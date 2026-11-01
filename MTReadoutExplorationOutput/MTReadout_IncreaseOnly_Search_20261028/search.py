"""Exploratory population-model grid; all noise precedes directional pooling.

The coarse grid uses 5-degree direction spacing and 32 antithetic draws.
Finalists must be checked with finer grids, more draws and all trials.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
from pathlib import Path
import json,csv,time,argparse
import numpy as np
from scipy.special import erf, ndtri
from scipy.fft import rfft,irfft
from research import ROOT,load_data,gain_fit

def response(data,indices,sigma,spacing=5,background=1):
    theta=np.arange(0,360,spacing,dtype=np.float32)
    dirs=data['dirs'][indices].astype(np.float32)
    delta=(theta[None,None,None,:]-dirs[...,None]+180)%360-180
    gauss=np.exp(-.5*(delta/sigma)**2)*data['mask'][indices,...,None]
    gbar=np.sqrt(2*np.pi)*sigma/360*erf(180/(np.sqrt(2)*sigma))
    out=[]
    for epoch in ['pre','post']:
        c=data[epoch][indices].astype(np.float32)
        out.append(np.sum(c[...,None]*gauss,axis=2)+(100-c.sum(-1))[...,None]*gbar*background)
    # N x physical patch x epoch x direction
    return np.stack(out,axis=2).astype(np.float32)

def noise_fields(n,draws,spacing,corr,seed,style='gp'):
    rng=np.random.default_rng(seed);d=int(round(360/spacing))
    white=rng.standard_normal((n,draws//2,2,2,d),dtype=np.float32)
    a=np.arange(d)*2*np.pi/d;ell=np.deg2rad(corr)
    cov=np.exp(-2*np.sin(a/2)**2/ell**2)
    spec=np.sqrt(np.maximum(rfft(cov).real,0)).astype(np.float32)
    z=irfft(rfft(white,axis=-1)*spec,n=d,axis=-1)
    if style=='balanced':
        z-=z.mean(axis=-1,keepdims=True)
        z/=np.sqrt(1-cov.mean())
    elif style=='opponent':
        ab=rng.standard_normal((n,draws//2,2,2,2),dtype=np.float32)
        harmonic=ab[...,0,None]*np.cos(a)+ab[...,1,None]*np.sin(a)
        z=(z+harmonic.astype(np.float32))/np.sqrt(2.)
    elif style!='gp':raise ValueError(style)
    return np.concatenate([z,-z],axis=1)

def pool(x,p):
    if np.isinf(p):return np.max(x,axis=-1)
    # Scaled to avoid overflow; scaling is identical for both RDKs.
    return np.sum((x/100)**p,axis=-1)**(1/p)*100

def evidence(r,z,sd,p,family,variance='constant'):
    # Genuine rate rectification precedes all readout variants.
    factor=1
    if variance=='sqrt':factor=np.sqrt(np.maximum(r[:,None],0)/30)
    elif variance=='proportional':factor=r[:,None]/30
    elif variance!='constant':raise ValueError(variance)
    x=np.maximum(r[:,None]+sd*z*factor,0)
    if family.startswith('adapted_'):
        _,mode,strength=family.split('_');strength=float(strength)
        # Prestep activity supplies a direction-specific divisive template.
        # It is held fixed for the comparison, including its estimation noise.
        weight=1/(1+strength*x[:,:,:,0,:]/30)
        x=x*weight[:,:,:,None,:]
        family='pooled_delta' if mode=='pooled' else 'rectified_delta'
    if family in ['pooled_delta','pooled_square_delta']:
        e=pool(x,p)
        if family=='pooled_square_delta':e=(e/100)**2
        return e[:,:,:,1]-e[:,:,:,0]
    if family=='post_only':return pool(x[:,:,:,1],p)
    a=.5 if family=='half_reference' else 1.
    delta=x[:,:,:,1]-a*x[:,:,:,0]
    e=pool(np.maximum(delta,0),p)
    if family=='opponent_delta':e-=pool(np.maximum(-delta,0),p)
    return e

def probabilities(e,change):
    # All cross-patch draw pairs are legitimate because patch noise is independent.
    a=np.take_along_axis(e,change[:,None,None],axis=2)[:,:,0]
    b=np.take_along_axis(e,1-change[:,None,None],axis=2)[:,:,0]
    # Exact same all-pairs statistic, including half-credit ties, without
    # constructing a draws-by-draws array. This makes large-draw checks cheap.
    if not np.all(np.isfinite(e)):raise ValueError('Nonfinite simulated evidence')
    b=np.sort(b,axis=1)
    return np.array([(np.searchsorted(bb,aa,side='left').sum()+
                      np.searchsorted(bb,aa,side='right').sum())/(2*len(aa)*len(bb))
                     for aa,bb in zip(a,b)])

def choose_trials(groups,n=180,seed=20261025):
    rng=np.random.default_rng(seed);selection=[];specs=[]
    for g in groups:
        local=np.sort(rng.choice(len(g['y']),min(n,len(g['y'])),replace=False))
        start=len(selection);selection.extend(g['indices'][local]);stop=len(selection)
        specs.append((g,local,slice(start,stop)))
    return np.array(selection),specs

def evaluate(pc,specs):
    rows=[];score=0.;maxerr=0
    for g,local,sl in specs:
        pred=pc[sl];gg,nll,ok=gain_fit(g['x'][local],g['step'][local],g['alpha'][local],g['beta'],g['lapse'],pred)
        ratio=float(gg[1]/gg[0]) if gg[0]>.02 else float('nan')
        target=float(g['original'][1]/g['original'][0]);acc=float(pred.mean());obs=float(g['y'].mean())
        # Native 120 has only one session; retain it as a diagnostic, not a target.
        use=g['label']!='IDS_120'
        if use:
            score+=((ratio-target)/.15)**2 if np.isfinite(ratio) else 1000
            score+=((acc-obs)/.05)**2
            maxerr=max(maxerr,abs(acc-obs))
        rows.append(dict(label=g['label'],ratio=ratio,drift=float(gg[0]),probe=float(gg[1]),
                         other=gg[2:].tolist(),acc=acc,observed_acc=obs,observed_ratio=target,fit_ok=ok))
    return score,maxerr,rows

def run_grid(n=180,draws=32):
    data,groups=load_data();indices,specs=choose_trials(groups,n)
    np.save(ROOT/'coarse_indices.npy',indices)
    families=['pooled_delta','rectified_delta','opponent_delta','half_reference','post_only']
    out=[];start=time.monotonic()
    with open(ROOT/'coarse_search.jsonl','w') as f:
        for corr in [15,30,60]:
            z=noise_fields(len(indices),draws,5,corr,1729)
            for sigma in [25,37.5,50]:
                r=response(data,indices,sigma)
                for sd in [2.5,5,10,20]:
                    for p in [2,3,5,8,float('inf')]:
                        for family in families:
                            e=evidence(r,z,sd,p,family)
                            pc=probabilities(e,data['change'][indices])
                            score,maxerr,rows=evaluate(pc,specs)
                            rec=dict(family=family,sigma=sigma,corr=corr,sd=sd,p=p,score=score,
                                     max_accuracy_error=maxerr,groups=rows)
                            f.write(json.dumps(rec)+'\n');out.append(rec)
                f.flush()
                print(f'{len(out)} models; sigma={sigma}, corr={corr}; {time.monotonic()-start:.1f}s',flush=True)
    best=sorted(out,key=lambda r:r['score'])[:20]
    (ROOT/'coarse_best.json').write_text(json.dumps(best,indent=2))
    print(json.dumps(best[:2],indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--n',type=int,default=180);parser.add_argument('--draws',type=int,default=32)
    args=parser.parse_args();run_grid(args.n,args.draws)
