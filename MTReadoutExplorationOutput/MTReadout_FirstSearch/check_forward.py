"""Numerical checks relevant to the scientific interpretation."""
import numpy as np
from search import noise_fields,response,evidence,probabilities,pool
from research import load_data

def main():
    # Independent patches/windows, with the intended circular covariance.
    z=noise_fields(1000,32,5,30,918)
    v=z[:,:,0,0,:].reshape(-1,72)
    assert abs(v.var()-1)<.02
    expected=np.exp(-2*np.sin(np.deg2rad(30)/2)**2/np.deg2rad(30)**2)
    empirical=np.mean(v[:,0]*v[:,6])
    assert abs(empirical-expected)<.03,(empirical,expected)
    for a in [(0,1),(1,0),(1,1)]:
        cross=np.mean(z[:,:,0,0,:]*z[:,:,a[0],a[1],:]);assert abs(cross)<.01
    # A full uniform MT bank conserves directional total response. At p=1
    # response to a coherence step therefore vanishes (apart from quadrature).
    data,_=load_data();idx=np.flatnonzero(data['keep']&(data['domain']=='IDS'))[:50]
    r=response(data,idx,37.5,1)
    assert np.max(np.abs(r[:,:,1,:].sum(-1)-r[:,:,0,:].sum(-1)))<.02
    # Two identical patches have chance performance, including exact ties.
    e=np.repeat(np.arange(32)[None,:,None],2,axis=2)
    assert probabilities(e,np.array([0]))[0]==.5
    # A noiseless single-stream increment improves its patch under signed
    # pooled differences. Uniform rotation must not alter the continuous bank.
    idx=np.flatnonzero((data['domain']=='IDJ')&(np.max(data['allocation'],axis=1)>.999))[:40]
    r1=response(data,idx,37.5,1);r5=response(data,idx,37.5,5)
    for p in [2,5,np.inf]:
        e1=pool(r1,p);e5=pool(r5,p)
        # Remove the trivial number-of-units factor for a quadrature check.
        if not np.isinf(p):e5*=5**(1/p)
        assert np.max(abs(e1-e5))/np.max(abs(e1))<.001
    from temporal import temporal_inputs,temporal_response
    cohs=temporal_inputs();idx=np.flatnonzero(data['keep'])[:100]
    rb=temporal_response(data,cohs,idx,37.5,5)
    rw=response(data,idx,37.5,5)
    assert np.allclose(rb.mean(axis=3),rw,atol=2e-5)
    zz=noise_fields(len(idx),32,5,30,111)
    e1=evidence(rw,zz,5,5,'pooled_delta')
    e2=evidence(rw,zz,5,5,'adapted_pooled_0')
    assert np.array_equal(e1,e2)
    print('Covariance, temporal/patch independence, p=1 conservation, 2AFC ties, and direction quadrature checks passed.')
    print('Temporal reconstruction and zero-strength adaptation checks passed.')

if __name__=='__main__':main()
