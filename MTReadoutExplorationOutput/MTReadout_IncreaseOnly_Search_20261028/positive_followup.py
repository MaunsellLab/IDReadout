"""Larger-sample screening, all-trial confirmation, and task transfer."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,time,argparse
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from scipy.special import ndtr
from scipy.optimize import minimize_scalar
from positive_search import OUT,predict,score_rows,width_classes
from research import load_data
from search import choose_trials,evaluate
from refine import gain_own_scale,probit_slopes

def read_records(path):
    with open(path) as f:return [json.loads(line) for line in f if line.strip()]

def unique_configs(records):
    seen=set();result=[]
    for r in records:
        key=json.dumps({k:v for k,v in r['config'].items() if k not in ['branch','name']},sort_keys=True)
        if key not in seen:seen.add(key);result.append(r['config'])
    return result

_STATE=None
def init(stage):
    global _STATE
    data,groups=load_data()
    if stage=='screen':idx,specs=choose_trials(groups,n=1000,seed=281042)
    else:
        idx=np.concatenate([g['indices'] for g in groups]);specs=[];start=0
        for g in groups:
            specs.append((g,np.arange(len(g['y'])),slice(start,start+len(g['y']))));start+=len(g['y'])
    _STATE=data,groups,idx,specs,stage

def work(item):
    j,c=item;data,groups,idx,specs,stage=_STATE;t=time.monotonic()
    pc=predict(data,idx,c,draws=64 if stage=='screen' else 128,spacing=3 if stage=='screen' else 2,
               seed=71123 if stage=='screen' else 81029)
    _,_,rows=evaluate(pc,specs);score,shape,acc=score_rows(rows)
    if stage!='screen':
        for r,(g,local,sl) in zip(rows,specs):r['own_scale']=gain_own_scale(g,pc[sl])
        np.savez_compressed(OUT/f"{c.get('name','final_'+str(j))}_predictions.npz",indices=idx,p_correct=pc)
    return dict(config=c,score=score,shape_score=shape,accuracy_score=acc,groups=rows,
                n=len(idx),draws=64 if stage=='screen' else 128,spacing=3 if stage=='screen' else 2,
                elapsed=time.monotonic()-t)

def run_stage(stage,configfile,outfile,workers):
    cc=json.loads((OUT/configfile).read_text());out=[]
    with ProcessPoolExecutor(workers,initializer=init,initargs=(stage,)) as pool:
        for r in pool.map(work,enumerate(cc),chunksize=1):
            out.append(r);(OUT/outfile).write_text(json.dumps(out,indent=2))
            print(stage,len(out),'/',len(cc),'score',round(r['score'],2),r['config'],
                  [(g['label'],round(g['ratio'],3),round(g['acc'],3)) for g in r['groups']],flush=True)

def transfer(configfile):
    data,_=load_data();cc=json.loads((OUT/configfile).read_text());rng=np.random.default_rng(40319)
    groups=[];allidx=[]
    for offset in [10,25,45,90,135,179,180]:
        k=np.flatnonzero((data['domain']=='IDR')&data['keep']&(data['offset']==offset))
        k=np.sort(rng.choice(k,min(1500,len(k)),replace=False));groups.append((offset,k));allidx.extend(k)
    idx=np.array(allidx);noise=data['post'][idx]-data['detpost'][idx]
    x=noise[np.arange(len(idx)),data['change'][idx]];x=np.column_stack([x[:,0],x[:,1:].sum(-1)])
    observed=[];start=0
    for offset,k in groups:
        sl=slice(start,start+len(k));start+=len(k)
        slope=probit_slopes(x[sl],data['step'][k],data['y'][k])
        observed.append(dict(offset=offset,n=len(k),accuracy=float(data['y'][k].mean()),
                             drift_slope=float(slope[0]),probe_slope=float(slope[1]),ratio=float(slope[1]/slope[0])))
    (OUT/'idr_observed_diagnostic.json').write_text(json.dumps(observed,indent=2))
    resultfile=OUT/'idr_positive_transfer.json'
    rows=json.loads(resultfile.read_text()) if resultfile.exists() else []
    current_names={v['name'] for v in cc}
    rows=[v for v in rows if v['config'].get('name') not in current_names]
    for c in cc:
        for p in c.get('transfer_p',[1.5,2,3,5,10,float('inf')]):
            cfg=c|dict(p=p);pc=predict(data,idx,cfg,draws=128,spacing=2,seed=92743);start=0;result=[]
            for offset,k in groups:
                sl=slice(start,start+len(k));start+=len(k)
                slope=probit_slopes(x[sl],data['step'][k],pc[sl])
                result.append(dict(offset=offset,n=len(k),accuracy=float(pc[sl].mean()),
                                   drift_slope=float(slope[0]),probe_slope=float(slope[1]),ratio=float(slope[1]/slope[0])))
            rows.append(dict(config=cfg,groups=result))
            resultfile.write_text(json.dumps(rows,indent=2))
            (OUT/f"idr_{c['name']}.json").write_text(json.dumps([v for v in rows if v['config'].get('name')==c['name']],indent=2))
            print('IDR',c.get('name'),p,[(g['offset'],round(g['ratio'],3),round(g['accuracy'],3)) for g in result],flush=True)
        # IDJ, same MT and internal noise: no task-specific noise calibration.
        ji=np.flatnonzero((data['domain']=='IDJ')&data['keep']);allocation=np.max(data['allocation'][ji],axis=1)
        out=[]
        for p in c.get('idj_p',sorted(set([3,5,c['p']]))):
            cfg=c|dict(p=p);pc=predict(data,ji,cfg,draws=128,spacing=2,seed=13403);result=[]
            for a in [1.,.7,.5,1/3]:
                mask=abs(allocation-a)<1e-4
                result.append(dict(allocation=a,n=int(mask.sum()),observed_accuracy=float(data['y'][ji[mask]].mean()),
                                   predicted_accuracy=float(pc[mask].mean())))
            out.append(dict(config=cfg,groups=result))
            print('IDJ',c.get('name'),p,[(g['allocation'],round(g['predicted_accuracy'],3)) for g in result],flush=True)
        (OUT/f"idj_{c['name']}.json").write_text(json.dumps(out,indent=2))

def checks(configfile):
    data,_=load_data();cc=json.loads((OUT/configfile).read_text());out=[]
    for c in cc:
        old=np.load(OUT/f"{c['name']}_predictions.npz");take=np.arange(0,len(old['indices']),12)
        pc=predict(data,old['indices'][take],c,draws=512,spacing=1,seed=932035)
        delta=pc-old['p_correct'][take];widths,_=width_classes(c)
        out.append(dict(config=c,n=len(take),new_draws=512,new_spacing=1,
                        mean_difference=float(delta.mean()),rmse=float(np.sqrt(np.mean(delta**2))),
                        widths=widths.tolist(),actual_width_cv=float(widths.std()/widths.mean())))
        (OUT/'resolution_checks.json').write_text(json.dumps(out,indent=2));print('check',out[-1],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['screen','confirm','transfer','checks']);
    p.add_argument('configfile');p.add_argument('--outfile',default='screened.json');p.add_argument('--workers',type=int,default=3)
    a=p.parse_args()
    if a.stage in ['transfer','checks']:globals()[a.stage](a.configfile)
    else:run_stage(a.stage,a.configfile,a.outfile,a.workers)
