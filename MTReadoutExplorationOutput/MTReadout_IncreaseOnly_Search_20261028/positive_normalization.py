"""Restricted, calibrated dot-distribution normalization diagnostic.

Denominator is a Gaussian-weighted sum of the physical dot distribution.
Baseline and preferred response at 15% coherence are calibrated exactly;
remaining single-direction deviations are measured, not assumed absent.
Internal noise is applied to the resulting MT output, before delta/readout.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from positive_search import OUT,initialize,worker,response

def calibration(c):
    cohs=np.linspace(0,100,101)
    d=dict(dirs=np.zeros((len(cohs),2,3)),mask=np.zeros((len(cohs),2,3),bool),
           pre=np.zeros((len(cohs),2,3)),post=np.zeros((len(cohs),2,3)))
    d['mask'][:,:,0]=True;d['post'][:,:,0]=cohs[:,None]
    r,_=response(d,np.arange(len(cohs)),c,1)
    ref,_=response(d,np.arange(len(cohs)),c|dict(norm_strength=0),1)
    pref=r[:,0,1,:,0];null=r[:,0,1,:,180]
    pref_ref=ref[:,0,1,:,0];null_ref=ref[:,0,1,:,180]
    delta=abs(r-ref)
    # Relative to the FULL 0--100% response excursion, plus low-coherence metric.
    return dict(preferred_max_fraction_dynamic_range=float(np.max(abs(pref-pref_ref)/(pref_ref[-1]-pref_ref[0]))),
                null_max_fraction_dynamic_range=float(np.max(abs(null-null_ref)/(null_ref[0]-null_ref[-1]))),
                low_coherence_0_30_max_response_error=float(delta[:31].max()),
                preferred_15_error=float(np.max(abs(pref[15]-pref_ref[15]))),
                output_tuning_15_rms_error=float(np.sqrt(np.mean((r[15,0,1]-ref[15,0,1])**2))))

def run():
    cc=[]
    for sigma in [25,37.5]:
        for k in [.05,.15,.35,.6,.85]:
            for s in [40,65,100,140]:
                for h in [0,1]:
                    for b in [1.2,2]:
                        for p in [3,5,float('inf')]:
                            cc.append(dict(sigma=sigma,width_cv=0,corr=30,sd=5,
                                           harmonic=h,global_fraction=0,balanced=False,
                                           variance_exp=b,p=p,noise='gaussian',branch='normalization',
                                           norm_strength=k,norm_sigma=s))
    start=time.monotonic();out=[]
    with open(OUT/'normalization.jsonl','w') as f,ProcessPoolExecutor(2,initializer=initialize,
                initargs=(220,32,6,20261028)) as pool:
        for r in pool.map(worker,enumerate(cc),chunksize=1):
            r['calibration']=calibration(r['config']);out.append(r);f.write(json.dumps(r)+'\n')
            if len(out)%48==0:
                f.flush();best=min(out,key=lambda v:v['score'])
                print('normalization',len(out),'/',len(cc),round(time.monotonic()-start,1),
                      'best',round(best['score'],2),best['config'],flush=True)
    (OUT/'normalization_best.json').write_text(json.dumps(sorted(out,key=lambda v:v['score'])[:30],indent=2))

if __name__=='__main__':run()
