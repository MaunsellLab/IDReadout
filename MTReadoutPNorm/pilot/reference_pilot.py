"""Independent Python numerical reference for the MATLAB model.

Run with converted v5 MATLAB IDS recordings (scipy.io.loadmat):
 python reference_pilot.py --files SESSION.mat ... --out pilot --draws 64
This pilot does not replace gain recovery through the existing IDS pipeline.
"""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
from scipy.io import loadmat
from scipy.linalg import cholesky
from scipy.optimize import minimize_scalar
from scipy.special import erf


def value(x):
    while isinstance(x, dict) and 'data' in x:
        x = x['data']
    return x


def trace(tr, prefix, nframe, frame_ms):
    times = np.asarray(value(tr[prefix+'TimesMS']), float).ravel()
    out = np.full((nframe, 3), np.nan)
    for s, role in enumerate(('Drift', 'Probe0', 'Probe1')):
        c = np.asarray(value(tr[prefix+role+'CohsPC']), float).ravel()
        assert len(c) == len(times)
        for q, t in enumerate(times):
            a = max(0, int(np.floor(t/frame_ms)))
            b = min(nframe, int(np.floor(times[q+1]/frame_ms)) if q+1<len(times) else nframe)
            if b>a:
                out[a:b, s] = c[q]
    assert np.isfinite(out).all()
    return out


def load_trials(files):
    rows, coh = [], []
    for file in files:
        s = loadmat(file, simplify_cells=True)
        h = s['sessionHeader']
        assert h['createdBy'] == 'convertIDSData' and h['animal'] == 'Meetz'
        for k, tr in enumerate(np.atleast_1d(s['trials'])):
            if not all(f in tr for f in ('trial', 'trialCertify', 'extendedEOT', 'targetChosen')):
                continue
            td = tr['trial']['data']
            if value(tr['trialCertify']) != 0 or value(tr['extendedEOT']) not in (0,1) or td['instructTrial']:
                continue
            if any(bool(td.get(f, False)) for f in ('topUpTrial','isTopUpTrial','topUp','isTopUp','topupTrial','isTopupTrial')):
                continue
            assert td['randomizeDriftDir']
            assert td['preStepMS'] == 750 and td['stepMS'] == 250
            frame_ms = 1000/h['frameRateHz']; nframe = round(1000/frame_ms)
            t = np.arange(nframe)*frame_ms
            ip = (t>=500)&(t<750); istep = (t>=750)&(t<1000)
            assert ip.sum() == istep.sum()
            change = int(td['stepSide']); choice = int(value(tr['targetChosen']))
            correct = value(tr['extendedEOT']) == 0
            assert correct == (choice == change)
            c = np.zeros((2,2,3))
            for j,prefix in ((change,'change'),(1-change,'noChange')):
                noise = trace(tr,prefix,nframe,frame_ms)
                c[j,0] = noise[ip].mean(0) + [td['baseCohPC'],td['cohNoisePC'],td['cohNoisePC']]
                endpoint = td['stepCohPC'] if j == change else td['baseCohPC']
                c[j,1] = noise[istep].mean(0) + [endpoint,td['cohNoisePC'],td['cohNoisePC']]
            rows.append(dict(sessionID=h['fileName'],parentTrialIdx=k+1,
                probeOffsetDeg=float(td['probeOffsetDeg']),driftDirDeg=float(td['driftDirDeg']),
                stepMagnitudePC=float(abs(td['stepCohPC']-td['baseCohPC'])),
                changeRDK=change,choiceRDK=choice,correct=bool(correct),
                hasStepNoise=bool(td['stepHasNoise']),isIncrement=bool(td['stepSign']==1)))
            coh.append(c)
    assert rows
    return rows, np.array(coh)


def responses(rows,coh,spacing=1):
    theta = np.arange(0,360,spacing)
    gbar = np.sqrt(2*np.pi)*37.5/360*erf(180/(np.sqrt(2)*37.5))
    assert coh.min()>-1e-4 and coh.sum(-1).max()<100.0001
    c = np.maximum(coh,0)
    r = np.zeros((len(rows),2,2,len(theta)))
    for k,tr in enumerate(rows):
        ds = tr['driftDirDeg']+np.array([0,tr['probeOffsetDeg'],-tr['probeOffsetDeg']])
        d = (theta[:,None]-ds+180)%360-180
        g = np.exp(-.5*(d/37.5)**2)
        r[k] = (100-c[k].sum(-1))[...,None]*gbar + c[k]@g.T
    return r


def fields(draws=64,seed=1729,spacing=1):
    theta = np.deg2rad(np.arange(0,360,spacing))
    d = theta[:,None]-theta
    k = np.exp(-2*np.sin(d/2)**2/np.deg2rad(30)**2)
    l = cholesky(k+1e-10*np.eye(len(theta)),lower=True)
    rng = np.random.default_rng(seed)
    z = np.zeros((2,2,len(theta),draws))
    for j in range(2):
        for e in range(2):
            a = l@rng.standard_normal((len(theta),draws//2))
            z[j,e] = np.concatenate([a,-a],axis=1)
    return z


def pool(r,p,axis=-2):
    m = r.max(axis=axis,keepdims=True)
    if np.isinf(p):
        return m.squeeze(axis)
    den = np.where(m==0,1,m)
    return (m*(np.power(r/den,p).sum(axis=axis,keepdims=True))**(1/p)).squeeze(axis)


def predict(r,p,sd,z,batch=32):
    prob = np.empty(len(r)); clipped=0; total=0
    for a in range(0,len(r),batch):
        rb = r[a:a+batch]; delta=np.zeros((len(rb),2,z.shape[-1]))
        for j in range(2):
            for e in range(2):
                x = rb[:,j,e,:,None]+sd*z[j,e]
                clipped += np.count_nonzero(x<0); total += x.size
                delta[:,j] += (1 if e==1 else -1)*pool(np.maximum(x,0),p)
        for q in range(len(rb)):
            x=delta[q,0]; y=np.sort(delta[q,1])
            less=np.searchsorted(y,x,side='left'); leq=np.searchsorted(y,x,side='right')
            prob[a+q]=np.mean((less+leq)/2)/len(y)
    return prob,clipped/max(total,1)


def checks():
    assert np.isclose(pool(np.array([[3.],[4.]]),2)[0],5)
    assert pool(np.zeros((3,2)),5).sum() == 0
    r=np.zeros((1,2,2,2)); z=np.zeros((2,2,2,4))
    r[0,0,0]=[20,10]; r[0,0,1]=[19,9]
    r[0,1,0]=[20,10]; r[0,1,1]=[15,5]
    assert predict(r,5,0,z)[0][0] == 1
    r[:,1]=r[:,0]
    assert predict(r,5,0,z)[0][0] == .5
    # The nonnegative raw-response p=1 sum conserves total dot drive.
    rows=[dict(driftDirDeg=0,probeOffsetDeg=120)]
    c=np.array([[[[7,7,7],[25,7,7]],[[7,7,7],[7,7,7]]]],float)
    a=responses(rows,c)
    assert np.ptp(a.sum(-1)) < 1e-4
    # Rotate stimuli AND a realization of the noise: decisions are invariant.
    z=fields(32); p,_=predict(a,5,2,z)
    p2,_=predict(np.roll(a,41,-1),5,2,np.roll(z,41,-2))
    assert np.array_equal(p,p2)
    # Pooled changes retain pedestal-dependent nonlinear interactions.
    zero=np.array([[[[7,0,0],[25,0,0]],[[7,0,0],[7,0,0]]]],float)
    b=responses(rows,zero)
    da=pool(a[0,0,1,:,None],5)[0]-pool(a[0,0,0,:,None],5)[0]
    db=pool(b[0,0,1,:,None],5)[0]-pool(b[0,0,0,:,None],5)[0]
    assert not np.isclose(da,db)
    print('Analytic reference checks passed.',flush=True)


def run(args):
    checks()
    rows,c=load_trials(args.files); r=responses(rows,c)
    mask=np.array([t['hasStepNoise'] and t['isIncrement'] for t in rows])
    y=np.array([t['choiceRDK']==0 for t in rows],float)
    change=np.array([t['changeRDK'] for t in rows])
    z=fields(args.draws); zf=fields(args.final_draws,seed=1730)
    out=Path(args.out); out.mkdir(parents=True,exist_ok=True)
    result=dict(specification='full response, signed pooled pre/post differences, four independent directional fields',
                nTrials=len(rows),nFitTrials=int(mask.sum()),fitDraws=args.draws,finalDraws=args.final_draws,
                absoluteRateFloor='zero; pooled differences remain signed',fits=[],sessions=[])
    for sid in sorted(set(t['sessionID'] for t in rows)):
        q=np.array([t['sessionID']==sid for t in rows]); qm=q&mask
        result['sessions'].append(dict(sessionID=sid,nValid=int(q.sum()),nFit=int(qm.sum()),
            probeOffsetDeg=sorted(set(t['probeOffsetDeg'] for t in rows if t['sessionID']==sid)),
            fitAccuracy=float(np.mean([rows[i]['correct'] for i in np.flatnonzero(qm)]))))
    for p in (5,np.inf):
        history=[]
        def fun(sd):
            pr,cf=predict(r[mask],p,sd,z)
            pr=np.clip(pr,1e-8,1-1e-8)
            nll=-np.sum(y[mask]*np.log(pr)+(1-y[mask])*np.log(1-pr))
            history.append(dict(noiseSD=float(sd),negLogLikelihood=float(nll),negativeRateFraction=float(cf)))
            print('p',p,'sd',round(sd,4),'nll',round(nll,3),flush=True)
            return nll
        grid=np.array([0,.5,1,2,4,8,16,32,40])
        vals=np.array([fun(sd) for sd in grid]); b=vals.argmin()
        opt=minimize_scalar(fun,bounds=(grid[max(b-1,0)],grid[min(b+1,len(grid)-1)]),
            method='bounded',options=dict(xatol=.05,maxiter=12))
        sd=float(grid[b]) if vals[b] < opt.fun else float(opt.x)
        prob,cf=predict(r,p,sd,zf); pc=np.where(change==0,prob,1-prob)
        clamped=np.clip(prob,1e-8,1-1e-8)
        ll=np.sum(y[mask]*np.log(clamped[mask])+(1-y[mask])*np.log(1-clamped[mask]))
        result['fits'].append(dict(p='Inf' if np.isinf(p) else p,noiseSD=sd,finalLogLikelihood=float(ll),
            predictedAccuracy=float(pc[mask].mean()),negativeRateFraction=float(cf),
            profile=history,perSession=[dict(sessionID=sid,predictedAccuracy=float(pc[np.array([t['sessionID']==sid for t in rows])&mask].mean()))
               for sid in sorted(set(t['sessionID'] for t in rows))]))
        np.savez_compressed(out/f"prediction_p{'Inf' if np.isinf(p) else int(p)}.npz",pRDK1=prob,pCorrect=pc)
    (out/'Meetz_IDS_pilot.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='fits'},indent=2),flush=True)
    for fit in result['fits']:
        print({k:v for k,v in fit.items() if k!='profile'},flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--files',nargs='+',required=True)
    parser.add_argument('--out',required=True)
    parser.add_argument('--draws',type=int,default=64)
    parser.add_argument('--final-draws',type=int,default=256)
    run(parser.parse_args())
