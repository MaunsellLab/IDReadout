"""Joint IDS/IDQ and IDR feasibility search; fixed MT/noise, task-specific p.

IDR gains are matched-window probit diagnostics, not the unavailable historical
kernel gain analysis. Observed choices are fitted with the identical diagnostic.
Variance exponents 3 and 4 are deliberately relaxed mechanistic diagnostics.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,time
import numpy as np
from scipy.stats import qmc
from concurrent.futures import ProcessPoolExecutor
from positive_search import OUT,response,fields,evidence,score_rows
from positive_normalization import calibration
from research import load_data
from search import choose_trials,probabilities,evaluate
from refine import probit_slopes

_STATE=None
def init():
    global _STATE
    d,g=load_data();ids,spec=choose_trials(g,n=300,seed=20261028)
    rng=np.random.default_rng(40319);idr=[];idr_spec=[]
    for o in [10,45,90,135,180]:
        i=np.flatnonzero((d['domain']=='IDR')&d['keep']&(d['offset']==o));i=np.sort(rng.choice(i,300,replace=False))
        noise=d['post'][i]-d['detpost'][i];x=noise[np.arange(len(i)),d['change'][i]];x=np.c_[x[:,0],x[:,1:].sum(-1)]
        sl=slice(len(idr),len(idr)+len(i));idr.extend(i);idr_spec.append((o,i,x,sl))
    _STATE=d,ids,spec,np.array(idr),idr_spec

def worker(item):
    j,c=item;d,ids,spec,idr,idr_spec=_STATE
    r,w=response(d,ids,c,6);z=fields(len(ids),32,6,c,22473)
    pc=probabilities(evidence(r,z,w,c),d['change'][ids]);_,_,rows=evaluate(pc,spec);ss,shape,ac=score_rows(rows)
    ci=c|dict(p=c['idr_p']);r,w=response(d,idr,ci,6);z=fields(len(idr),32,6,ci,22841)
    pc=probabilities(evidence(r,z,w,ci),d['change'][idr]);rr=[];idrscore=0
    targets={10:(.9,.2),45:(.25,.15),90:(.25,.15),135:(.2,.15),180:(.5,.2)}
    for o,i,x,sl in idr_spec:
        coeff=probit_slopes(x,d['step'][i],pc[sl]);ratio=float(coeff[1]/coeff[0])
        target,scale=targets[o]
        idrscore+=((ratio-target)/scale)**2 if coeff[0]>.001 and np.isfinite(ratio) else 1e4
        idrscore+=((float(pc[sl].mean())-float(d['y'][i].mean()))/.1)**2
        rr.append(dict(offset=o,n=len(i),ratio=ratio,drift_slope=float(coeff[0]),probe_slope=float(coeff[1]),acc=float(pc[sl].mean())))
    return dict(number=j,config=c,score=ss+idrscore,ids_score=ss,idr_score=idrscore,
                groups=rows,idr=rr,calibration=calibration(c) if c.get('norm_strength',0) else None)

def configs():
    out=[]
    points=qmc.Sobol(10,scramble=True,seed=69203).random_base2(10)
    for j,v in enumerate(points):
        k=[0,0,.15,.6,.85][j%5]
        c=dict(sigma=float(18+27*v[0]),width_cv=0 if v[1]<.5 else float(.65*(2*v[1]-1)),
               corr=float(10*9**v[2]),sd=float(1.5*10**v[3]),
               harmonic=1. if v[4]>.75 else float(v[4]/.75),global_fraction=0,
               balanced=bool(v[5]>.5),variance_exp=[0,1.2,2,3,4][(j//5)%5],
               p=5 if v[6]<.5 else float('inf'),idr_p=[1.5,2,3][j%3],
               noise='gaussian' if v[7]<.5 else 'lognormal',branch='joint',
               norm_strength=k,norm_sigma=float(40+100*v[8]))
        out.append(c)
    # Dense local sweep of the promising pure first-harmonic covariance.
    for sigma in [18,22,25,30,37.5]:
        for b in [1.2,2,3,4]:
            for sd in [3,5,7.5]:
                for p in [5,float('inf')]:
                    out.append(dict(sigma=sigma,width_cv=0,corr=30,sd=sd,harmonic=1,global_fraction=0,
                                    balanced=False,variance_exp=b,p=p,idr_p=2,noise='gaussian',branch='joint_local'))
    return out

if __name__=='__main__':
    cc=configs();start=time.monotonic();out=[]
    with open(OUT/'joint.jsonl','w') as f,ProcessPoolExecutor(3,initializer=init) as pool:
        for r in pool.map(worker,enumerate(cc),chunksize=1):
            out.append(r);f.write(json.dumps(r)+'\n')
            if len(out)%32==0:
                f.flush();best=min(out,key=lambda v:v['score'])
                print(len(out),'/',len(cc),round(time.monotonic()-start,1),'best',round(best['score'],2),best['config'],flush=True)
    (OUT/'joint_best.json').write_text(json.dumps(sorted(out,key=lambda v:v['score'])[:40],indent=2))
