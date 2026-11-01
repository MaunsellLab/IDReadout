"""Additional covariance and response-dependent noise hypotheses."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,time
import numpy as np
from research import ROOT,load_data
from search import choose_trials,noise_fields,response,evidence,probabilities,evaluate
from refine import predict,gain_own_scale

def run():
    data,groups=load_data();idx,specs=choose_trials(groups,180);out=[];start=time.monotonic()
    with open(ROOT/'noise_search.jsonl','w') as f:
        for style in ['balanced','opponent','gp']:
            z=noise_fields(len(idx),32,5,30,31415,style)
            for variance in ['constant','sqrt','proportional']:
                for sigma in [25,37.5]:
                    r=response(data,idx,sigma)
                    for sd in [1,2.5,5,10]:
                        for p in [3,5,float('inf')]:
                            for family in ['pooled_delta','pooled_square_delta','rectified_delta','opponent_delta']:
                                c=dict(style=style,variance=variance,sigma=sigma,corr=30,sd=sd,p=p,family=family)
                                e=evidence(r,z,sd,p,family,variance);pc=probabilities(e,data['change'][idx])
                                score,err,rows=evaluate(pc,specs)
                                rec=dict(config=c,score=score,groups=rows);out.append(rec);f.write(json.dumps(rec)+'\n')
                    f.flush();print(len(out),'models',style,variance,sigma,'best',min(r['score'] for r in out),'seconds',time.monotonic()-start,flush=True)
    # Confirm the best overall and best with the requested p=5/sigma=37.5.
    candidates=[min(out,key=lambda r:r['score']),
                min((r for r in out if r['config']['p']==5 and r['config']['sigma']==37.5),key=lambda r:r['score'])]
    indices=np.concatenate([g['indices'] for g in groups]);specs=[];start=0
    for g in groups:
        specs.append((g,np.arange(len(g['y'])),slice(start,start+len(g['y']))));start+=len(g['y'])
    confirmed=[]
    for j,rec in enumerate(candidates):
        c=rec['config'];pc=predict(data,indices,c,draws=128,spacing=2)
        score,err,rows=evaluate(pc,specs)
        for row,(g,loc,sl) in zip(rows,specs):row['own_scale']=gain_own_scale(g,pc[sl])
        confirmed.append(dict(config=c,score=score,groups=rows))
        np.savez_compressed(ROOT/f'noise_finalist_{j}_predictions.npz',indices=indices,p_correct=pc)
        (ROOT/'noise_confirmed.json').write_text(json.dumps(confirmed,indent=2))
        print('confirmed',c,score,flush=True)
    checks=[]
    for j,rec in enumerate(confirmed):
        old=np.load(ROOT/f'noise_finalist_{j}_predictions.npz');take=np.arange(0,len(old['indices']),20)
        new=predict(data,old['indices'][take],rec['config'],draws=256,spacing=1,seed=91781)
        delta=new-old['p_correct'][take]
        checks.append(dict(config=rec['config'],n=len(delta),mean_difference=float(delta.mean()),
                           rmse=float(np.sqrt(np.mean(delta**2)))))
    (ROOT/'noise_resolution_check.json').write_text(json.dumps(checks,indent=2))

if __name__=='__main__':run()
