"""Numerical access, audit, and gain recovery for the 2026-10-25 MT export.

No MATLAB engine is required. Saved summaries supply the original predictors
and psychometric calibration; source shards supply physical stimulus inputs.
"""
from pathlib import Path
import sys, json, csv,os,zipfile
import numpy as np
from scipy.io import loadmat
from scipy.optimize import minimize
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'MTReadoutExplorationExport'))
from read_exploration import load_bundle, validate_bundle

ROOT=Path(__file__).resolve().parent
_default_bundle=ROOT/'MTReadoutExplorationInputs.zip'
if not _default_bundle.exists():_default_bundle=ROOT.parent/'upload/MTReadoutExplorationInputs.zip'
BUNDLE=Path(os.environ.get('MT_INPUT_BUNDLE',str(_default_bundle)))

def gain_fit(x, step, alpha, beta, lapse, y, initial=None):
    """Same bounded Weibull gain objective as fitIDSNoiseGain.m.

    y may be observed correctness or a simulated probability of correctness.
    Analytic gradient avoids finite-difference optimization noise.
    """
    x=np.asarray(x,float); y=np.asarray(y,float)
    def fg(g):
        raw=step+x@g; e=np.maximum(raw,0); z=(e/alpha)**beta
        ez=np.exp(-z); p=np.clip(.5+(.5-lapse)*(1-ez),1e-12,1-1e-12)
        dp=(.5-lapse)*ez*beta/alpha*(np.maximum(e,1e-30)/alpha)**(beta-1)
        dp[raw<=0]=0
        nll=-np.sum(y*np.log(p)+(1-y)*np.log1p(-p))
        grad=x.T@((p-y)/(p*(1-p))*dp)
        return nll,grad
    if initial is None: initial=np.r_[1.,np.zeros(x.shape[1]-1)]
    opt=minimize(fg,initial,jac=True,bounds=[(-5,5)]*x.shape[1],method='SLSQP',
                 options={'ftol':1e-10,'maxiter':1500})
    return opt.x, float(opt.fun), bool(opt.success)

def summaries():
    if not (ROOT/'analyses').exists():
        (ROOT/'analyses').mkdir()
        with zipfile.ZipFile(BUNDLE) as archive:
            for name in archive.namelist():
                if name.startswith('analyses/') and name.endswith('.mat'):
                    (ROOT/'analyses'/Path(name).name).write_bytes(archive.read(name))
    return [loadmat(f,simplify_cells=True)['analysis'] for f in sorted((ROOT/'analyses').glob('*.mat'))]

def aligned_alpha(t,psych):
    sf=psych['sessionFits']; wb=psych['alignedWeibull']; beta=wb['betaWeibull']; lapse=wb['lapse']
    thresholds={int(i):v for i,v in zip(np.atleast_1d(sf['sessionIndex']),np.atleast_1d(sf['threshold']))}
    threshold=np.array([thresholds[int(i)] for i in np.atleast_1d(t['sessionIndex'])])
    scale=(-np.log(1-(psych['targetPerformance']-.5)/(.5-lapse)))**(1/beta)
    return threshold/scale, beta,lapse

def gain_groups():
    groups=[]
    for a in summaries():
        name=a['export_source_name']; s=a.get('acrossSummary',{})
        native=name.startswith('IDS_Probe') and 'noiseGain' in s
        idq=name=='IDS120CN_IDQRandom_All_AcrossSessionSummary.mat'
        if not (native or idq): continue
        t=s['trialTable']; alpha,beta,lapse=aligned_alpha(t,s['psychFit'])
        if native:
            x=np.column_stack([t['noisePredDrift'],np.asarray(t['noisePredPlus'])+t['noisePredMinus']])
            original=np.atleast_1d(s['noiseGain']['driftProbe']['gain']); label=f"IDS_{int(s['probeOffsetDeg'])}"
        else:
            x=np.column_stack([t['noisePredCDrift'],np.asarray(t['noisePredCPlus'])+t['noisePredCMinus'],
                               t['noisePredNDrift'],np.asarray(t['noisePredNPlus'])+t['noisePredNMinus']])
            original=np.atleast_1d(s['noiseGain']['sideRole']['gain']);label='IDQ_120'
        good=np.asarray(t['hasStepNoise'],bool)
        idx=np.flatnonzero(good)
        group=dict(label=label,x=x[idx],step=np.asarray(t['stepCoh'])[idx],alpha=alpha[idx],
                   beta=beta,lapse=lapse,y=np.asarray(t['correct'],float)[idx],original=original,
                   session=np.asarray(t['sourceSessionID'])[idx],parent=np.asarray(t['parentTrialIdx'])[idx])
        g,nll,ok=gain_fit(group['x'],group['step'],group['alpha'],beta,lapse,group['y'])
        group['recovered']=g;group['recovered_nll']=nll
        assert np.max(np.abs(g-original))<2e-4,(label,g,original,ok)
        groups.append(group)
    return sorted(groups,key=lambda g: (g['label'].split('_')[0],int(g['label'].split('_')[1])))

def build_data():
    pieces=[]; audit=[]; keys={};pos=0
    for source,row in load_bundle(BUNDLE):
        s=source;m=s['meta'];n=int(m['n_trials']);d=m['domain']
        pre=np.asarray(s['pre_window_total_pc'],float).reshape(n,2,3,order='F')
        post=np.asarray(s['step_window_total_pc'],float).reshape(n,2,3,order='F')
        mask=np.asarray(s['stream_mask'],bool).reshape(n,2,3,order='F')
        dirs=np.asarray(s['directions_deg'],float).reshape(n,2,3,order='F')
        ref=dirs[:,0,0];dirs=(dirs-ref[:,None,None]+180)%360-180
        dirs[~mask]=0;pre[~mask]=0;post[~mask]=0
        parent=np.asarray(s['parent_trial_index']).reshape(-1);side=np.asarray(s['change_side'],int).reshape(-1)-1
        detpre=np.asarray(s['deterministic_pre_pc'],float).reshape(n,2,3,order='F')
        detpost=np.asarray(s['deterministic_step_pc'],float).reshape(n,2,3,order='F')
        detpre[~mask]=0;detpost[~mask]=0
        step=(detpost-detpre)[np.arange(n),side].sum(-1)
        keep=np.asarray(s['has_step_noise'],bool).reshape(-1)&np.asarray(s['is_increment'],bool).reshape(-1)
        pieces.append(dict(pre=pre,post=post,detpre=detpre,detpost=detpost,dirs=dirs,mask=mask,step=step,change=side,
                           y=np.asarray(s['correct'],float).reshape(-1),keep=keep,
                           offset=np.asarray(s['probe_offset_deg'],float).reshape(-1),
                           allocation=np.asarray(s['nominal_allocation_fraction']).reshape(n,3,order='F'),
                           domain=np.repeat(d,n),session=np.repeat(m['session_id'],n),parent=parent))
        for i,p in enumerate(parent): keys[(d,m['session_id'],int(p))]=pos+i
        pos+=n
        audit.append(dict(domain=d,session=m['session_id'],source=m['source_name'],n_trials=n,n_fit=int(keep.sum()),
                          min_coherence=float(np.minimum(pre,post)[mask].min()),
                          max_budget=float(np.maximum(pre.sum(-1),post.sum(-1)).max()),
                          offsets=';'.join(map(str,np.unique(s['probe_offset_deg'])))))
    data={k:np.concatenate([p[k] for p in pieces]) for k in pieces[0]}
    groups=gain_groups()
    for g in groups:
        domain=g['label'].split('_')[0]
        g['indices']=np.array([keys[(domain,str(s),int(p))] for s,p in zip(g['session'],g['parent'])])
        assert np.allclose(data['step'][g['indices']],g['step'],atol=2e-5)
        assert np.array_equal(data['y'][g['indices']],g['y'])
    with open(ROOT/'source_audit.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(audit[0]));w.writeheader();w.writerows(audit)
    np.savez_compressed(ROOT/'trial_cache.npz',**data)
    return data,groups

def load_data():
    if (ROOT/'trial_cache.npz').exists():
        data=dict(np.load(ROOT/'trial_cache.npz'));groups=gain_groups()
        keys={(str(d),str(s),int(p)):i for i,(d,s,p) in enumerate(zip(data['domain'],data['session'],data['parent']))}
        for g in groups:
            domain=g['label'].split('_')[0]
            g['indices']=np.array([keys[(domain,str(s),int(p))] for s,p in zip(g['session'],g['parent'])])
        return data,groups
    return build_data()

if __name__=='__main__':
    validation=validate_bundle(BUNDLE)
    (ROOT/'validation.json').write_text(json.dumps(validation,indent=2))
    data,groups=build_data()
    out=[]
    for g in groups:
        out.append(dict(label=g['label'],n=len(g['y']),gains=g['recovered'].tolist(),
                        original=g['original'].tolist(),ratio=float(g['recovered'][1]/g['recovered'][0])))
    (ROOT/'observed_gains.json').write_text(json.dumps(out,indent=2))
    print(json.dumps(out,indent=2))
    print('negative coherence trials',int(((data['pre']<-.001)|(data['post']<-.001)).any((1,2)).sum()))
    print('overbudget trials',int(((data['pre'].sum(-1)>100.001)|(data['post'].sum(-1)>100.001)).any(1).sum()))
