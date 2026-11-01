"""Independent full-data Monte Carlo replication of the main p=5 results."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
import numpy as np
from concurrent.futures import ProcessPoolExecutor
from positive_followup import init,_STATE
from positive_search import OUT,predict,score_rows
from research import load_data
from search import evaluate
from refine import gain_own_scale

def worker(c):
    data,groups=load_data();idx=np.concatenate([g['indices'] for g in groups]);spec=[];start=0
    for g in groups:
        spec.append((g,np.arange(len(g['y'])),slice(start,start+len(g['y']))));start+=len(g['y'])
    t=time.monotonic();pc=predict(data,idx,c,draws=256,spacing=2,batch=48,seed=550729)
    _,_,rows=evaluate(pc,spec)
    for r,(g,local,sl) in zip(rows,spec):r['own_scale']=gain_own_scale(g,pc[sl])
    np.savez_compressed(OUT/f"{c['name']}_independent_predictions.npz",indices=idx,p_correct=pc)
    return dict(config=c,n=len(idx),draws=256,spacing=2,seed=550729,groups=rows,elapsed=time.monotonic()-t)

if __name__=='__main__':
    a=json.loads((OUT/'control_configs.json').read_text())[0]
    b=json.loads((OUT/'norm_confirm_configs.json').read_text())[-1]
    with ProcessPoolExecutor(2) as pool:
        out=list(pool.map(worker,[a,b]))
    (OUT/'independent_full_data.json').write_text(json.dumps(out,indent=2))
    for r in out:print(r['config']['name'],r['elapsed'],[(g['label'],round(g['ratio'],4),round(g['acc'],4)) for g in r['groups']],flush=True)
