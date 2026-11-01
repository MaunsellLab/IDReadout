"""All-trial confirmation and tests of transfer to the other tasks."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import json,time
import numpy as np
from scipy.optimize import minimize
from scipy.special import ndtr
from research import ROOT,load_data,gain_fit
from search import response,noise_fields,evidence,probabilities,evaluate

def predict(data,indices,config,draws=128,spacing=2,batch=96,seed=40219):
    out=np.empty(len(indices))
    for start in range(0,len(indices),batch):
        ii=indices[start:start+batch]
        r=response(data,ii,config['sigma'],spacing)
        z=noise_fields(len(ii),draws,spacing,config['corr'],seed+start,config.get('style','gp'))
        e=evidence(r,z,config['sd'],config['p'],config['family'],config.get('variance','constant'))
        out[start:start+len(ii)]=probabilities(e,data['change'][ii])
    return out

def gain_own_scale(g,pc):
    x=g['x'];step=g['step'];alpha=g['alpha'];beta=g['beta'];lapse=g['lapse']
    def fg(par):
        raw=step+x@par[:-1];e=np.maximum(raw,0);a=alpha*np.exp(par[-1]);z=(e/a)**beta
        ez=np.exp(-z);p=np.clip(.5+(.5-lapse)*(1-ez),1e-12,1-1e-12)
        dp=(.5-lapse)*ez*beta/a*(np.maximum(e,1e-30)/a)**(beta-1);dp[raw<=0]=0
        weight=(p-pc)/(p*(1-p));grad=np.r_[x.T@(weight*dp),np.sum(weight*(-(.5-lapse)*ez*beta*z))]
        return -np.sum(pc*np.log(p)+(1-pc)*np.log1p(-p)),grad
    opt=minimize(fg,np.r_[1,np.zeros(x.shape[1]-1),0.],jac=True,method='SLSQP',
                 bounds=[(-5,5)]*x.shape[1]+[(-2,2)],options={'ftol':1e-9,'maxiter':1000})
    return dict(gains=opt.x[:-1].tolist(),scale=float(np.exp(opt.x[-1])),
                ratio=float(opt.x[1]/opt.x[0]),success=bool(opt.success))

def finalists():
    data,groups=load_data();indices=np.concatenate([g['indices'] for g in groups]);specs=[];start=0
    for g in groups:
        specs.append((g,np.arange(len(g['y'])),slice(start,start+len(g['y']))));start+=len(g['y'])
    coarse=[json.loads(l) for l in open(ROOT/'coarse_search.jsonl')]
    configs=[dict(name='original_p5',family='pooled_delta',sigma=37.5,corr=30,sd=5.29842,p=5)]
    for family in ['pooled_delta','rectified_delta','opponent_delta']:
        best=min((r for r in coarse if r['family']==family),key=lambda r:r['score'])
        configs.append({k:best[k] for k in ['family','sigma','corr','sd','p']}|{'name':'coarse_best_'+family})
    out=[]
    for c in configs:
        t=time.monotonic();pc=predict(data,indices,c)
        score,maxerr,rows=evaluate(pc,specs)
        for row,(g,local,sl) in zip(rows,specs):
            row['own_scale']=gain_own_scale(g,pc[sl])
            pp=np.clip(pc[sl],1e-6,1-1e-6)
            row['choice_nll']=float(-np.sum(g['y']*np.log(pp)+(1-g['y'])*np.log1p(-pp)))
        out.append(dict(config=c,score=score,groups=rows,draws=128,spacing=2))
        np.savez_compressed(ROOT/(c['name']+'_predictions.npz'),indices=indices,p_correct=pc)
        (ROOT/'confirmed.json').write_text(json.dumps(out,indent=2))
        print(c['name'],round(time.monotonic()-t,1),'sec',[(r['label'],round(r['ratio'],3),round(r['own_scale']['ratio'],3),round(r['acc'],3)) for r in rows],flush=True)
    # An independent-seed, 1-degree check on the exact same selected trials.
    take=np.arange(0,len(indices),20);check=[]
    for c in configs:
        a=predict(data,indices[take],c,draws=256,spacing=1,seed=92823)
        old=np.load(ROOT/(c['name']+'_predictions.npz'))['p_correct'][take]
        check.append(dict(name=c['name'],n=len(a),mean_difference=float(np.mean(a-old)),
                          rmse=float(np.sqrt(np.mean((a-old)**2)))))
    (ROOT/'resolution_check.json').write_text(json.dumps(check,indent=2))

def idj_transfer(alternative=False):
    data,_=load_data();indices=np.flatnonzero((data['domain']=='IDJ')&data['keep'])
    allocation=np.max(data['allocation'][indices],axis=1);ref=allocation>.999
    # Same sigma/noise covariance as the initial agreed model. Internal SD is
    # calibrated only on single-stream allocation trials, then transferred.
    grid=[]
    for family in (['opponent_delta'] if alternative else ['pooled_delta','rectified_delta']):
        for p in ([3,5,float('inf')] if alternative else [2,3,5,float('inf')]):
            records=[]
            base=dict(family=family,p=p,sigma=37.5,corr=30)
            if alternative:base.update(sigma=25,style='balanced',variance='proportional')
            for sd in [0,.5,1,1.5,2.5,4,6,9,14]:
                c=base|dict(sd=sd)
                pc=predict(data,indices[ref],c,draws=64,spacing=5,seed=14317)
                pp=np.clip(pc,1e-5,1-1e-5);y=data['y'][indices[ref]]
                nll=float(-np.sum(y*np.log(pp)+(1-y)*np.log1p(-pp)))
                records.append((nll,sd))
            nll,sd=min(records);c=base|dict(sd=sd)
            pc=predict(data,indices,c,draws=128,spacing=2,seed=14317)
            rows=[]
            for a in [1.,.7,.5,1/3]:
                m=np.abs(allocation-a)<1e-4;y=data['y'][indices[m]];pp=np.clip(pc[m],1e-6,1-1e-6)
                rows.append(dict(allocation=a,n=int(m.sum()),observed_accuracy=float(y.mean()),
                                 predicted_accuracy=float(pp.mean()),
                                 choice_nll=float(-np.sum(y*np.log(pp)+(1-y)*np.log1p(-pp)))))
            rec=dict(config=c,calibration='single-stream trials only; shared SD across sessions',noise_grid=records,groups=rows)
            grid.append(rec);(ROOT/('idj_alternative.json' if alternative else 'idj_transfer.json')).write_text(json.dumps(grid,indent=2))
            print('IDJ',family,p,sd,[(r['allocation'],round(r['observed_accuracy'],3),round(r['predicted_accuracy'],3)) for r in rows],flush=True)

def probit_slopes(x,offset,y):
    """Diagnostic only: regression of full 2AFC probabilities on window noise.

    Each offset has an intercept. This is not claimed to reproduce the
    historical IDR kernel-weighted gain analysis, which is absent in the export.
    """
    unique=np.unique(offset);columns=[offset==o for o in unique]
    design=np.column_stack(columns+[x[:,0],x[:,1]])
    def fg(b):
        z=design@b;p=np.clip(ndtr(z),1e-9,1-1e-9);dp=np.exp(-z*z/2)/np.sqrt(2*np.pi)
        return -np.sum(y*np.log(p)+(1-y)*np.log1p(-p)),design.T@((p-y)/(p*(1-p))*dp)
    opt=minimize(fg,np.zeros(design.shape[1]),jac=True,method='BFGS',options={'gtol':1e-6})
    return opt.x[-2:]

def idr_transfer(alternative=False):
    data,_=load_data();rng=np.random.default_rng(911)
    groups=[];allidx=[]
    for offset in [10,25,45,90,135,179,180]:
        k=np.flatnonzero((data['domain']=='IDR')&data['keep']&(data['offset']==offset))
        k=np.sort(rng.choice(k,min(800,len(k)),replace=False));groups.append((offset,k));allidx.extend(k)
    indices=np.array(allidx);noise=data['post'][indices]-data['detpost'][indices]
    x=noise[np.arange(len(indices)),data['change'][indices]]
    x=np.column_stack([x[:,0],x[:,1:].sum(-1)])
    rows=[]
    for family in (['opponent_delta'] if alternative else ['pooled_delta','rectified_delta']):
        for p in [2,3,5,float('inf')]:
            c=dict(family=family,p=p,sigma=37.5,corr=30,sd=5.)
            if alternative:c.update(sigma=25,style='balanced',variance='proportional')
            pc=predict(data,indices,c,draws=128,spacing=2,seed=2113);start=0;result=[]
            for offset,k in groups:
                sl=slice(start,start+len(k));start+=len(k)
                # Also condition on nominal step to avoid its sessions being a
                # confound for sensitivity to the external fluctuations.
                level=data['step'][k]
                coeff=probit_slopes(x[sl],level,pc[sl])
                result.append(dict(offset=offset,n=len(k),accuracy=float(pc[sl].mean()),
                                   drift_slope=float(coeff[0]),probe_slope=float(coeff[1]),
                                   ratio=float(coeff[1]/coeff[0])))
            rows.append(dict(config=c,groups=result));(ROOT/('idr_alternative.json' if alternative else 'idr_transfer.json')).write_text(json.dumps(rows,indent=2))
            print('IDR',family,p,[(r['offset'],round(r['ratio'],3)) for r in result],flush=True)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('task',choices=['finalists','idj','idr','idj-alt','idr-alt']);a=p.parse_args()
    {'finalists':finalists,'idj':idj_transfer,'idr':idr_transfer,
     'idj-alt':lambda:idj_transfer(True),'idr-alt':lambda:idr_transfer(True)}[a.task]()
