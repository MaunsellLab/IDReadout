"""Check whether directional pooling before temporal integration matters.

The 50-ms exported bins limit this to an approximation of online pooling.
Internal noise is a window-level rate-estimation field, independent between
the two windows and RDKs, held over bins within a window.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,time
import numpy as np
from scipy.special import erf
from research import ROOT,BUNDLE,load_data
from read_exploration import load_bundle
from search import choose_trials,noise_fields,pool,probabilities,evaluate,response

def temporal_inputs():
    path=ROOT/'temporal_cache.npz'
    if path.exists(): return np.load(path)['cohs']
    chunks=[]
    for s,row in load_bundle(BUNDLE):
        n=int(s['meta']['n_trials']);b=len(np.atleast_1d(s['bin_start_ms']))
        noise=np.asarray(s['noise_pc']).reshape(n,2,3,b,order='F')
        mask=np.asarray(s['stream_mask'],bool).reshape(n,2,3,order='F')
        starts=np.asarray(s['bin_start_ms']);epochs=[]
        for basefield,a,z in [('deterministic_pre_pc',500,750),('deterministic_step_pc',750,1000)]:
            which=(starts>=a)&(starts<z);assert which.sum()==5
            base=np.asarray(s[basefield]).reshape(n,2,3,order='F').copy();base[~mask]=0
            x=base[...,None]+noise[...,which];x=np.where(mask[...,None],x,0)
            epochs.append(x.transpose(0,1,3,2))
        chunks.append(np.stack(epochs,axis=2))
    cohs=np.concatenate(chunks).astype(np.float32);np.savez_compressed(path,cohs=cohs)
    return cohs

def temporal_response(data,cohs,indices,sigma,spacing):
    theta=np.arange(0,360,spacing,dtype=np.float32);dirs=data['dirs'][indices].astype(np.float32)
    delta=(theta[None,None,None,:]-dirs[...,None]+180)%360-180
    g=np.exp(-.5*(delta/sigma)**2)*data['mask'][indices,...,None]
    c=cohs[indices];gbar=np.sqrt(2*np.pi)*sigma/360*erf(180/(np.sqrt(2)*sigma))
    r=np.einsum('npetk,npkd->npetd',c,g,optimize=True)+(100-c.sum(-1))[...,None]*gbar
    return r.astype(np.float32)

def predict_temporal(data,cohs,indices,config,draws=32,spacing=5,batch=96):
    out=np.empty(len(indices));p=config['p'];sd=config['sd'];family=config['family']
    for start in range(0,len(indices),batch):
        ii=indices[start:start+batch];r=temporal_response(data,cohs,ii,config['sigma'],spacing)
        z=noise_fields(len(ii),draws,spacing,config['corr'],13411+start)
        e=0
        if family=='pooled_delta':
            for b in range(5):
                x=np.maximum(r[:,None,:,:,b,:]+sd*z,0);t=pool(x,p)
                e=e+(t[:,:,:,1]-t[:,:,:,0])/5
        elif family=='rectified_delta':
            pre=np.maximum(r[:,:,0].mean(axis=2)[:,None]+sd*z[:,:,:,0],0)
            for b in range(5):
                post=np.maximum(r[:,None,:,1,b,:]+sd*z[:,:,:,1],0)
                e=e+pool(np.maximum(post-pre,0),p)/5
        else:raise ValueError(family)
        out[start:start+len(ii)]=probabilities(e,data['change'][ii])
    return out

def run():
    data,groups=load_data();cohs=temporal_inputs();idx,specs=choose_trials(groups,180)
    t=time.monotonic();out=[]
    for family in ['pooled_delta','rectified_delta']:
        for p in [2,3,5,8,float('inf')]:
            for sd in [0,1,2.5,5,10]:
                c=dict(family=family,p=p,sd=sd,sigma=37.5,corr=30)
                pc=predict_temporal(data,cohs,idx,c)
                score,err,rows=evaluate(pc,specs)
                out.append(dict(config=c,score=score,groups=rows))
            (ROOT/'temporal_search.json').write_text(json.dumps(out,indent=2))
            print(family,p,'best',min(r['score'] for r in out),'elapsed',time.monotonic()-t,flush=True)
    best=min(out,key=lambda r:r['score']);c=best['config'];idx=np.concatenate([g['indices'] for g in groups])
    specs=[];start=0
    for g in groups:
        specs.append((g,np.arange(len(g['y'])),slice(start,start+len(g['y']))));start+=len(g['y'])
    pc=predict_temporal(data,cohs,idx,c,draws=128,spacing=2)
    score,err,rows=evaluate(pc,specs)
    (ROOT/'temporal_confirmed.json').write_text(json.dumps(dict(config=c,score=score,groups=rows),indent=2))
    np.savez_compressed(ROOT/'temporal_predictions.npz',indices=idx,p_correct=pc)

if __name__=='__main__':run()
