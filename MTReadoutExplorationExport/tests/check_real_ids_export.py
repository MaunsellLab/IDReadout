"""Independent schema/reconstruction check using scipy-readable IDS inputs.

This checks numerical interchange using real recordings; it is not execution
of the MATLAB exporter. MATLAB integration is tested by its supplied fixture.
"""
import argparse,csv,io,json,sys,zipfile
from pathlib import Path
import numpy as np
from scipy.io import loadmat,savemat
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from read_exploration import validate_bundle


def value(x):
    while isinstance(x,dict) and 'data' in x: x=x['data']
    return x


def trace(tr,prefix,nf,frame_ms):
    ts=np.asarray(value(tr[prefix+'TimesMS']),float).ravel(); out=np.full((nf,3),np.nan)
    for s,role in enumerate(('Drift','Probe0','Probe1')):
        v=np.asarray(value(tr[prefix+role+'CohsPC']),float).ravel()
        assert len(v)==len(ts)
        for q,t in enumerate(ts):
            a=max(0,int(np.floor(t/frame_ms)))
            b=min(nf,int(np.floor(ts[q+1]/frame_ms)) if q+1<len(ts) else nf)
            if b>a: out[a:b,s]=v[q]
    assert np.isfinite(out).all()
    return out


def convert(path):
    A=loadmat(path,simplify_cells=True); H=A['sessionHeader']; trs=np.atleast_1d(A['trials'])
    used=[]
    for k,tr in enumerate(trs):
        if not all(f in tr for f in ['trial','trialCertify','extendedEOT','targetChosen']): continue
        D=tr['trial']['data']
        if value(tr['trialCertify'])!=0 or value(tr['extendedEOT']) not in (0,1) or D['instructTrial']: continue
        used.append((k,tr))
    n=len(used); frame_ms=1000/H['frameRateHz']; nf=round(1000/frame_ms)
    t=np.arange(nf)*frame_ms; noise=np.zeros((n,2,3,nf),np.float32)
    pre=np.zeros((n,2,3),np.float32); post=pre.copy(); dirs=pre.copy()
    parent=[]; change=[]; choice=[]; correct=[]; has=[]; offset=[]
    for i,(k,tr) in enumerate(used):
        D=tr['trial']['data']; c=int(D['stepSide'])
        noise[i,c]=trace(tr,'change',nf,frame_ms).T
        noise[i,1-c]=trace(tr,'noChange',nf,frame_ms).T
        pre[i]=[D['baseCohPC'],D['cohNoisePC'],D['cohNoisePC']]
        post[i]=pre[i]; post[i,c,0]=D['stepCohPC']
        dirs[i]=(D['driftDirDeg']+np.array([0,D['probeOffsetDeg'],-D['probeOffsetDeg']]))%360
        parent.append(k+1); change.append(c+1); choice.append(int(value(tr['targetChosen']))+1)
        correct.append(value(tr['extendedEOT'])==0); has.append(bool(D['stepHasNoise'])); offset.append(D['probeOffsetDeg'])
    starts=np.arange(0,1000,50.); stops=starts+50
    bins=np.stack([noise[..., (t>=a)&(t<z)].mean(-1) for a,z in zip(starts,stops)],-1)
    source=dict(meta=dict(schema_version=1,domain='IDS',animal='Meetz',session_id=H['fileName'],
        n_trials=n,n_fit_trials=int(sum(has)),pre_window_ms=[500,750],step_window_ms=[750,1000]),
        parent_trial_index=np.asarray(parent),change_side=np.asarray(change,np.uint8),choice_side=np.asarray(choice,np.uint8),
        correct=np.asarray(correct,np.uint8),has_step_noise=np.asarray(has,np.uint8),topup_trial=np.zeros(n,np.uint8),
        randomized_base_direction=np.ones(n),probe_offset_deg=np.array(offset),
        allocation_fraction=np.tile([1.,0,0],(n,1)),nominal_allocation_fraction=np.tile([1.,0,0],(n,1)),
        directions_deg=dirs,stream_mask=np.ones((n,2,3),np.uint8),yoking_group=np.tile([1,2,3],(n,2,1)).astype(np.uint8),
        deterministic_pre_pc=pre,deterministic_step_pc=post,noise_pc=bins,
        pre_window_total_pc=pre+noise[..., (t>=500)&(t<750)].mean(-1),
        step_window_total_pc=post+noise[..., (t>=750)&(t<1000)].mean(-1),
        bin_start_ms=starts,bin_stop_ms=stops,bin_frame_count=np.full(20,6))
    return source


def main():
    parser=argparse.ArgumentParser(); parser.add_argument('--files',nargs='+',required=True)
    parser.add_argument('--output',required=True); args=parser.parse_args(); rows=[]
    with zipfile.ZipFile(args.output,'w',zipfile.ZIP_DEFLATED) as z:
        for k,p in enumerate(args.files):
            source=convert(p); name=f'sources/source_{k:03d}.mat'
            stream=io.BytesIO(); savemat(stream,dict(source=source),do_compression=True)
            z.writestr(name,stream.getvalue()); m=source['meta']
            rows.append(dict(bundle_file=name,domain='IDS',n_trials=m['n_trials']))
        text=io.StringIO(); writer=csv.DictWriter(text,fieldnames=['bundle_file','domain','n_trials']);writer.writeheader();writer.writerows(rows)
        z.writestr('SourceManifest.csv',text.getvalue());z.writestr('Schema.json',json.dumps(dict(schema_version=1)))
    print(json.dumps(validate_bundle(args.output),indent=2))


if __name__=='__main__': main()
