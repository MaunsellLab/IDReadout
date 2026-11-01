"""Read and validate plain-numeric MT input bundles (numpy + scipy only).

python read_exploration.py MTReadoutExplorationInputs.zip
"""
from __future__ import annotations
import argparse
import csv
import io
import json
from pathlib import PurePosixPath
import zipfile
import numpy as np
from scipy.io import loadmat


def array(source, name, shape, dtype=None):
    # MATLAB/scipy may squeeze the trial dimension for a one-trial source.
    return np.asarray(source[name], dtype=dtype).reshape(shape, order='F')


def validate_source(source):
    meta=source['meta']; n=int(meta['n_trials'])
    starts=np.asarray(source['bin_start_ms'],float).reshape(-1)
    stops=np.asarray(source['bin_stop_ms'],float).reshape(-1)
    counts=np.asarray(source['bin_frame_count'],float).reshape(-1)
    b=len(starts)
    assert n>0 and len(stops)==b and len(counts)==b
    assert np.all(stops>starts) and np.all(counts>0)
    assert np.allclose(stops[:-1],starts[1:])
    mask=array(source,'stream_mask',(n,2,3),bool)
    dirs=array(source,'directions_deg',(n,2,3),float)
    pre=array(source,'deterministic_pre_pc',(n,2,3),float)
    post=array(source,'deterministic_step_pc',(n,2,3),float)
    noise=array(source,'noise_pc',(n,2,3,b),float)
    assert np.isfinite(dirs[mask]).all()
    assert np.isfinite(pre[mask]).all() and np.isfinite(post[mask]).all()
    active=np.broadcast_to(mask[...,None],noise.shape)
    assert np.isfinite(noise[active]).all()
    assert not np.any(np.asarray(source['topup_trial'],bool))
    change=np.asarray(source['change_side'],int).reshape(-1)
    choice=np.asarray(source['choice_side'],int).reshape(-1)
    correct=np.asarray(source['correct'],bool).reshape(-1)
    assert len(change)==len(choice)==len(correct)==n
    assert np.isin(change,[1,2]).all() and np.isin(choice,[1,2]).all()
    assert np.array_equal(correct,change==choice)
    parent=np.asarray(source['parent_trial_index'],int).reshape(-1)
    assert len(np.unique(parent))==n
    for epoch,baseline,field in [
        (meta['pre_window_ms'],pre,'pre_window_total_pc'),
        (meta['step_window_ms'],post,'step_window_total_pc')]:
        a,z=np.asarray(epoch,float).reshape(-1)
        idx=(starts>=a-1e-8)&(stops<=z+1e-8)
        assert idx.any() and np.isclose(starts[idx][0],a) and np.isclose(stops[idx][-1],z)
        weights=counts[idx]/counts[idx].sum()
        reconstructed=baseline+np.sum(noise[...,idx]*weights,axis=-1)
        saved=array(source,field,(n,2,3),float)
        assert np.allclose(reconstructed[mask],saved[mask],atol=2e-5,rtol=1e-6),field
    allocations=array(source,'allocation_fraction',(n,3),float)
    nominal=array(source,'nominal_allocation_fraction',(n,3),float)
    assert np.allclose(allocations.sum(-1),1,atol=1e-4)
    if meta['domain']=='IDQ':
        assert np.all(np.asarray(source['randomized_base_direction'])==1)
    if meta['domain']=='IDJ':
        equal=(np.abs(allocations-1/3).max(-1)<1e-4)
        assert np.allclose(nominal[equal],1/3)
    yoking=array(source,'yoking_group',(n,2,3),int)
    for a in range(3):
        for c in range(a+1,3):
            same=mask[:,:,a]&mask[:,:,c]&(yoking[:,:,a]==yoking[:,:,c])
            if same.any():
                assert np.array_equal(noise[:,:,a,:][same],noise[:,:,c,:][same]),'Yoked traces differ'
    return dict(domain=meta['domain'],animal=meta['animal'],session_id=meta['session_id'],
                n_trials=n,n_fit_trials=int(meta['n_fit_trials']),n_bins=b)


def load_bundle(path):
    """Yield (source, manifest row); never extract arbitrary ZIP paths."""
    with zipfile.ZipFile(path) as archive:
        manifest=list(csv.DictReader(io.StringIO(archive.read('SourceManifest.csv').decode('utf-8-sig'))))
        for row in manifest:
            name=row['bundle_file'].replace('\\','/')
            p=PurePosixPath(name)
            assert not p.is_absolute() and '..' not in p.parts and p.parts[0]=='sources'
            values=loadmat(io.BytesIO(archive.read(name)),simplify_cells=True)
            yield values['source'],row


def validate_bundle(path):
    with zipfile.ZipFile(path) as archive:
        schema=json.loads(archive.read('Schema.json'))
    assert schema['schema_version']==1
    result=[]; keys=set()
    for source,row in load_bundle(path):
        info=validate_source(source)
        assert info['n_trials']==int(row['n_trials'])
        assert info['domain']==row['domain']
        for parent in np.asarray(source['parent_trial_index']).reshape(-1):
            key=(info['domain'],info['animal'],info['session_id'],int(parent))
            assert key not in keys, f'Duplicate source trial: {key}'
            keys.add(key)
        result.append(info)
    assert result
    return dict(n_sources=len(result),n_trials=sum(r['n_trials'] for r in result),
                n_fit_trials=sum(r['n_fit_trials'] for r in result),
                domains=sorted(set(r['domain'] for r in result)),sources=result)


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('bundle')
    args=parser.parse_args()
    print(json.dumps(validate_bundle(args.bundle),indent=2))
