"""Physical invariances and covariance/moment checks for the expanded model."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from positive_search import OUT,response,fields,evidence,noisy_rates
from research import load_data
from search import response as old_response,evidence as old_evidence,probabilities
from positive_normalization import calibration

def run():
    tied=np.array([[[0,0],[1,1],[2,2]],[[0,0],[1,0],[2,0]],[[0,1],[0,1],[0,1]]],float)
    assert np.allclose(probabilities(tied,np.zeros(3,int)),[.5,5/6,0])
    d,_=load_data();i=np.arange(20)
    c=dict(sigma=37.5,width_cv=0,corr=30,sd=2.5,harmonic=.7,global_fraction=.1,
           balanced=True,variance_exp=1.2,p=5,noise='gaussian')
    r,w=response(d,i,c,3);z=fields(len(i),64,3,c,810)
    old=old_response(d,i,37.5,3)
    meanerr=float(np.max(abs(r[...,0,:]-old)));assert meanerr<2e-5
    errs=[]
    for b,v in [(0,'constant'),(1,'sqrt'),(2,'proportional')]:
        a=evidence(r,z,w,c|dict(variance_exp=b));bval=old_evidence(old,z,c['sd'],5,'rectified_delta',v)
        err=float(np.max(abs(a-bval)));assert err<5e-5;errs.append(err)
    # Same physical distribution represented as one or two coincident streams.
    dd={k:v[:20].copy() for k,v in d.items()};dd['dirs'][:]=0;dd['mask'][:]=True
    for ep in ['pre','post']:
        dd[ep][:]=0;dd[ep][:,:,0]=10;dd[ep][:,:,1]=5
    c=c|dict(norm_strength=.35,norm_sigma=140)
    r1,_=response(dd,i,c,3)
    for ep in ['pre','post']:dd[ep][:,:,0]+=dd[ep][:,:,1];dd[ep][:,:,1]=0
    dd['mask'][:,:,1]=False;r2,_=response(dd,i,c,3)
    doterr=float(np.max(abs(r1-r2)));assert doterr<2e-5
    # Rotations commute with all mean and readout operations.
    base={k:v[:20].copy() for k,v in d.items()};rr,ww=response(base,i,c,3)
    base['dirs']=(base['dirs']+90+180)%360-180;rot,_=response(base,i,c,3)
    rotationerr=float(np.max(abs(rot-np.roll(rr,30,axis=-1))));assert rotationerr<3e-5
    ee=evidence(rr,z,ww,c);ee2=evidence(rot,np.roll(z,30,axis=-1),ww,c)
    assert np.max(abs(ee-ee2))<1e-4
    # Stationary covariance remains PSD by construction; sample check plus
    # independence of the two patches and two epochs.
    z=fields(1200,16,6,c,301)[:,::2] # retain a subsample, both signs still valid
    flat=z.reshape(-1,60);co=np.mean(flat*flat[:,0,None],axis=0)
    a=np.arange(60)*2*np.pi/60;cov=np.exp(-2*np.sin(a/2)**2/np.deg2rad(c['corr'])**2)
    cov=(cov-cov.mean())/(1-cov.mean())
    expected=.2*cov+.7*np.cos(a)+.1
    coverr=float(np.max(abs(co-expected)));assert coverr<.04
    epochcorr=float(np.mean(z[:,:,:,0]*z[:,:,:,1]));patchcorr=float(np.mean(z[:,:,0]*z[:,:,1]));assert abs(epochcorr)<.03 and abs(patchcorr)<.03
    # Noise moments, away from the Gaussian floor.
    moments=[]
    for style in ['gaussian','lognormal']:
        rng=np.random.default_rng(193);q=rng.standard_normal((4000,16,2,2,1));mu=20.;sd=5.;b=1.2
        scale=sd*(mu/30)**(b/2)
        x=noisy_rates(np.full((4000,2,2,1,1),mu),q,c|dict(noise=style,sd=sd,variance_exp=b))
        moments.append(dict(style=style,mean=float(x.mean()),variance=float(x.var()),expected_mean=mu,expected_variance=scale**2))
        assert abs(x.mean()-mu)<.05 and abs(x.var()/scale**2-1)<.025
    out=dict(homogeneous_mean_error=meanerr,old_evidence_errors=errs,coincident_stream_error=doterr,
             rotation_error=rotationerr,covariance_max_error=coverr,epoch_covariance=epochcorr,patch_covariance=patchcorr,
             noise_moments=moments,normalization_calibration=calibration(c),rank_tie_cases_passed=True,passed=True)
    (OUT/'expanded_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))

if __name__=='__main__':run()
