import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
from concurrent.futures import ProcessPoolExecutor
from positive_search import OUT,initialize,worker
cc=json.loads((OUT/'harmonic_configs.json').read_text())
if __name__=='__main__':
    start=time.monotonic();out=[]
    with open(OUT/'harmonics.jsonl','w') as f,ProcessPoolExecutor(2,initializer=initialize,
            initargs=(220,32,6,20261028)) as pool:
        for r in pool.map(worker,enumerate(cc),chunksize=1):
            out.append(r);f.write(json.dumps(r)+'\n')
            if len(out)%48==0:
                f.flush();best=min(out,key=lambda x:x['score'])
                print(len(out),'/',len(cc),round(time.monotonic()-start,1),'best',round(best['score'],2),best['config'],flush=True)
