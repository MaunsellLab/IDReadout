"""Coordinated increase-only search, October 2026.

Uses actual physical stimulus windows and the original IDS/IDQ gain analysis.
The three width classes are a quadrature of a direction-uniform population;
they share the direction-dependent noise field. Epochs and patches are independent.
No decreases are pooled, no offset-dependent parameters, no decision noise.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,time,argparse
import numpy as np
from scipy.special import erf
from scipy.stats import qmc
from scipy.fft import rfft,irfft
from scipy.optimize import brentq
from concurrent.futures import ProcessPoolExecutor
from research import ROOT,load_data
from search import choose_trials,probabilities,evaluate

OUT=ROOT/'positive-search'
OUT.mkdir(exist_ok=True)

def width_classes(c):
    cv=c.get('width_cv',0.)
    if cv==0:return np.array([c['sigma']]),np.array([1.])
    # Three equal-number classes, preserving the arithmetic mean sigma.
    a=np.sqrt(1.5*np.log1p(cv*cv))
    scale=np.exp(np.array([-a,0,a]));scale/=scale.mean()
    return c['sigma']*scale,np.full(3,1/3)

def response(data,indices,c,spacing):
    widths,weights=width_classes(c)
    theta=np.arange(0,360,spacing,dtype=np.float32)
    delta=(theta[None,None,None,:]-data['dirs'][indices,...,None]+180)%360-180
    out=[]
    for sigma in widths:
        b=np.sqrt(2*np.pi)*sigma/360*erf(180/(np.sqrt(2)*sigma))
        k=c.get('norm_strength',0.)
        input_sigma=sigma;input_b=b
        if k:
            s=c['norm_sigma'];bs=np.sqrt(2*np.pi)*s/360*erf(180/(np.sqrt(2)*s))
            q=1/bs-1
            input_b=1/(1+(1/b-1)*(1+k*.15*q)+k*q)
            input_sigma=brentq(lambda v:np.sqrt(2*np.pi)*v/360*erf(180/(np.sqrt(2)*v))-input_b,.01,180)
            suppress_gauss=np.exp(-.5*(delta/s)**2)*data['mask'][indices,...,None]
        gauss=np.exp(-.5*(delta/input_sigma)**2)*data['mask'][indices,...,None]
        epochs=[]
        for epoch in ['pre','post']:
            coh=data[epoch][indices]
            mean=(np.sum(coh[...,None]*gauss,axis=2)+(100-coh.sum(-1))[...,None]*input_b)*(b/input_b)
            if k:
                suppress=np.sum(coh[...,None]*suppress_gauss,axis=2)+(100-coh.sum(-1))[...,None]*bs
                mean=mean/(1+k*(suppress/(100*bs)-1))
            mean=mean*(c['sigma']/sigma)**c.get('width_gain',0.)
            epochs.append(mean)
        out.append(np.stack(epochs,axis=2))
    return np.stack(out,axis=-2).astype(np.float32),weights

def fields(n,draws,spacing,c,seed):
    rng=np.random.default_rng(seed);d=int(round(360/spacing));a=np.arange(d)*2*np.pi/d
    white=rng.standard_normal((n,draws//2,2,2,d),dtype=np.float32)
    cov=np.exp(-2*np.sin(a/2)**2/np.deg2rad(c['corr'])**2)
    spec=np.sqrt(np.maximum(rfft(cov).real,0)).astype(np.float32)
    z=irfft(rfft(white,axis=-1)*spec,n=d,axis=-1)
    if c.get('balanced',False):
        z-=z.mean(-1,keepdims=True);z/=np.sqrt(1-cov.mean())
    ab=rng.standard_normal((n,draws//2,2,2,3),dtype=np.float32)
    order=c.get('harmonic_order',1)
    harmonic=ab[...,0,None]*np.cos(order*a)+ab[...,1,None]*np.sin(order*a)
    h=c.get('harmonic',0.);g=c.get('global_fraction',0.)
    z=np.sqrt(max(0,1-h-g))*z+np.sqrt(h)*harmonic+np.sqrt(g)*ab[...,2,None]
    z=z.astype(np.float32)
    return np.concatenate([z,-z],axis=1)

def noisy_rates(r,z,c):
    mean=r[:,None]
    variance_mean=mean[:,:,:,0:1] if c.get('freeze_variance',False) else mean
    scale=c['sd']*np.maximum(variance_mean/30,1e-12)**(c['variance_exp']/2)
    if c.get('noise','gaussian')=='lognormal':
        s=np.sqrt(np.log1p((scale/np.maximum(mean,1e-12))**2))
        raw=mean*np.exp(s*z[...,None,:]-.5*s*s)
    else:raw=mean+scale*z[...,None,:]
    return raw

def evidence(r,z,weights,c,diagnostic=False):
    raw=noisy_rates(r,z,c)
    x=np.maximum(raw,0)
    delta=np.maximum(x[:,:,:,1]-x[:,:,:,0],0)
    p=c['p']
    if np.isinf(p):e=delta.max(axis=(-1,-2))
    else:e=100*np.sum((delta/100)**p*weights[None,None,None,:,None],axis=(-1,-2))**(1/p)
    if diagnostic:return e,dict(floored=float(np.mean(raw<0)),mean_rate=float(r.mean()))
    return e

def predict(data,indices,c,draws=128,spacing=2,batch=64,seed=81029):
    out=np.empty(len(indices))
    for start in range(0,len(indices),batch):
        ii=indices[start:start+batch]
        r,w=response(data,ii,c,spacing)
        z=fields(len(ii),draws,spacing,c,seed+start)
        out[start:start+len(ii)]=probabilities(evidence(r,z,w,c),data['change'][ii])
    return out

def score_rows(rows):
    score=0.;shape=0.;acc=0.
    scales={'IDS_10':.15,'IDS_25':.18,'IDS_45':.15,'IDQ_120':.15,'IDS_155':.22,'IDS_170':.25}
    for r in rows:
        if r['label'] not in scales:continue
        if not np.isfinite(r['ratio']) or r['drift']<.15:return 1e5,1e5,1e5
        shape+=((r['ratio']-r['observed_ratio'])/scales[r['label']])**2
        acc+=((r['acc']-r['observed_acc'])/.08)**2
    return shape+acc,shape,acc

_STATE=None
def initialize(n,draws,spacing,seed):
    global _STATE
    data,groups=load_data();indices,specs=choose_trials(groups,n=n,seed=seed)
    _STATE=(data,indices,specs,draws,spacing)

def worker(item):
    number,c=item;data,idx,specs,draws,spacing=_STATE
    # Identical random primitives for parameter comparisons.
    r,w=response(data,idx,c,spacing);z=fields(len(idx),draws,spacing,c,22473)
    pc=probabilities(evidence(r,z,w,c),data['change'][idx])
    _,_,rows=evaluate(pc,specs);score,shape,acc=score_rows(rows)
    return dict(number=number,config=c,score=score,shape_score=shape,accuracy_score=acc,groups=rows)

def configs(power=10):
    points=qmc.Sobol(9,scramble=True,seed=428261).random_base2(power)
    out=[]
    # Equal-sized searches: empirical width / heterogeneity / free mean width.
    for branch in ['fixed_width','heterogeneous','free_width']:
        for i,v in enumerate(points):
            sigma=37.5 if branch!='free_width' else 18+32*v[0]
            cv=0 if branch=='fixed_width' else .65*v[1]
            b=[0,1,1.2,1.6,2][i%5]
            h=.98*v[4];g=(1-h)*.3*v[5]
            c=dict(sigma=float(sigma),width_cv=float(cv),corr=float(10*9**v[2]),
                   sd=float(1.25*12**v[3]),harmonic=float(h),global_fraction=float(g),
                   balanced=bool(v[6]>.5),variance_exp=b,
                   p=float('inf') if v[7]>.85 else float(1.5*8**(v[7]/.85)),
                   noise='lognormal' if v[8]>.5 else 'gaussian',branch=branch)
            out.append(c)
    # Replicate earlier models and isolate broad covariance trends at p=5.
    for sigma in [25,37.5]:
        for b in [1.2,2]:
            for h in [0,.5,.9,1]:
                for p in [3,5,8,float('inf')]:
                    out.append(dict(sigma=sigma,width_cv=0,corr=30,sd=5,harmonic=h,
                                    global_fraction=0,balanced=False,variance_exp=b,p=p,
                                    noise='gaussian',branch='reference'))
    return out

def run(power=8,n=220,draws=32,workers=3):
    cc=configs(power);start=time.monotonic();out=[]
    with open(OUT/'coordinated.jsonl','w') as f,ProcessPoolExecutor(workers,initializer=initialize,
            initargs=(n,draws,6,20261028)) as pool:
        for rec in pool.map(worker,enumerate(cc),chunksize=1):
            out.append(rec);f.write(json.dumps(rec)+'\n')
            if len(out)%48==0:
                f.flush();best=min(out,key=lambda x:x['score'])
                print(len(out),'/',len(cc),'seconds',round(time.monotonic()-start,1),
                      'best',round(best['score'],2),best['config'],flush=True)
    (OUT/'coordinated_best.json').write_text(json.dumps(sorted(out,key=lambda r:r['score'])[:60],indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--power',type=int,default=8);p.add_argument('--n',type=int,default=220)
    p.add_argument('--draws',type=int,default=32);p.add_argument('--workers',type=int,default=3)
    a=p.parse_args();run(a.power,a.n,a.draws,a.workers)
