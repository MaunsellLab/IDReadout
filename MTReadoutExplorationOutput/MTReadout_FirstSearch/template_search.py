"""Noisy prestep templates used as divisive adaptation, not just subtraction."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,time
import numpy as np
from research import ROOT,load_data
from search import choose_trials,noise_fields,response,evidence,probabilities,evaluate
from refine import predict,gain_own_scale

def run():
    data,groups=load_data();idx,specs=choose_trials(groups,180);out=[];t=time.monotonic()
    z=noise_fields(len(idx),32,5,30,81151)
    with open(ROOT/'template_search.jsonl','w') as f:
        for sigma in [25,37.5]:
            r=response(data,idx,sigma)
            for strength in [.3,1,3,10]:
                for mode in ['pooled','rectified']:
                    for sd in [1,2.5,5,10]:
                        for p in [3,5,float('inf')]:
                            fam=f'adapted_{mode}_{strength}'
                            c=dict(sigma=sigma,corr=30,sd=sd,p=p,family=fam)
                            e=evidence(r,z,sd,p,fam);pc=probabilities(e,data['change'][idx])
                            score,err,rows=evaluate(pc,specs)
                            rec=dict(config=c,score=score,groups=rows);out.append(rec);f.write(json.dumps(rec)+'\n')
                f.flush();print(len(out),'models',sigma,strength,'best',min(r['score'] for r in out),'seconds',time.monotonic()-t,flush=True)
    rec=min(out,key=lambda r:r['score']);idx=np.concatenate([g['indices'] for g in groups]);specs=[];start=0
    for g in groups:
        specs.append((g,np.arange(len(g['y'])),slice(start,start+len(g['y']))));start+=len(g['y'])
    pc=predict(data,idx,rec['config']);score,err,rows=evaluate(pc,specs)
    for row,(g,local,sl) in zip(rows,specs):row['own_scale']=gain_own_scale(g,pc[sl])
    np.savez_compressed(ROOT/'template_predictions.npz',indices=idx,p_correct=pc)
    (ROOT/'template_confirmed.json').write_text(json.dumps(dict(config=rec['config'],score=score,groups=rows),indent=2))
    print('confirmed',rec['config'],score,flush=True)

if __name__=='__main__':run()
