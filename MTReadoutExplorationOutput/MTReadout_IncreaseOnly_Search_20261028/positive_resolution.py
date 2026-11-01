"""Separate common-noise direction-grid checks and independent Monte Carlo checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
import numpy as np
from positive_search import OUT,response,evidence,predict
from research import load_data
from search import choose_trials,probabilities,evaluate

def run():
    d,g=load_data();idx,spec=choose_trials(g,300,seed=420219)
    allc=json.loads((OUT/'pilot_configs.json').read_text())+json.loads((OUT/'control_configs.json').read_text())+json.loads((OUT/'norm_confirm_configs.json').read_text())
    names=['pilot_6','reference_p5','norm_moderate_p5','norm_broad_strong_p5'];out=[]
    for name in names:
        c=next(v for v in allc if v['name']==name)
        rng=np.random.default_rng(421205);p1=[];p2=[]
        for start in range(0,len(idx),32):
            ii=idx[start:start+32];ab=rng.standard_normal((len(ii),32,2,2,2)).astype(np.float32);pcs=[]
            for spacing in [1,2]:
                angle=np.deg2rad(np.arange(0,360,spacing))
                z=ab[...,0,None]*np.cos(angle)+ab[...,1,None]*np.sin(angle);z=np.concatenate([z,-z],axis=1).astype(np.float32)
                r,w=response(d,ii,c,spacing);e=evidence(r,z,w,c);pcs.append(probabilities(e,d['change'][ii]))
            p1.extend(pcs[0]);p2.extend(pcs[1])
        p1=np.array(p1);p2=np.array(p2);delta=p1-p2
        pcnew=predict(d,idx,c,draws=512,spacing=2,batch=32,seed=721742)
        _,_,a=evaluate(p1,spec);_,_,b=evaluate(p2,spec);_,_,mc=evaluate(pcnew,spec)
        out.append(dict(config=c,n=len(idx),common_noise_grid=dict(mean_pc_difference=float(delta.mean()),
                       max_pc_difference=float(abs(delta).max()),rmse=float(np.sqrt(np.mean(delta**2))),
                       gain_1deg=a,gain_2deg=b),independent_512_draws=dict(groups=mc,
                       mean_pc_difference=float((pcnew-p2).mean()),rmse=float(np.sqrt(np.mean((pcnew-p2)**2))))))
        (OUT/'separated_resolution_checks.json').write_text(json.dumps(out,indent=2))
        print(name,'grid',out[-1]['common_noise_grid']['max_pc_difference'],
              'MC mean',out[-1]['independent_512_draws']['mean_pc_difference'],flush=True)
if __name__=='__main__':run()
