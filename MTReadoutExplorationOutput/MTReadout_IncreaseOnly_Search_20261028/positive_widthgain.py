"""Width heterogeneity with fixed peak versus fixed isotropic rate."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time
from concurrent.futures import ProcessPoolExecutor
from positive_joint import init,worker
from positive_search import OUT
if __name__=='__main__':
    cc=[]
    for cv in [.25,.5,.65]:
        for gain in [0,.5,1]:
            for b in [1.2,2]:
                for sd in [2.5,5,7.5]:
                    for h in [0,1]:
                        for p in [5,float('inf')]:
                            cc.append(dict(sigma=37.5,width_cv=cv,width_gain=gain,corr=30,sd=sd,harmonic=h,
                                           global_fraction=0,balanced=False,variance_exp=b,p=p,idr_p=2,
                                           noise='gaussian',branch='width_gain'))
    start=time.monotonic();out=[]
    with open(OUT/'width_gain_joint.jsonl','w') as f,ProcessPoolExecutor(2,initializer=init) as pool:
        for r in pool.map(worker,enumerate(cc),chunksize=1):
            out.append(r);f.write(json.dumps(r)+'\n')
            if len(out)%36==0:
                f.flush();best=min(out,key=lambda r:r['score'])
                print(len(out),'/',len(cc),round(time.monotonic()-start,1),'best',round(best['score'],2),best['config'],flush=True)
