"""Assemble the reproducible increase-only search report and scientific figure."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import csv,json,sys,platform,hashlib,zipfile
from pathlib import Path
import numpy as np
import scipy
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from positive_search import OUT,width_classes
from positive_normalization import calibration
from research import ROOT,BUNDLE

def load(name):return json.loads((OUT/name).read_text())
def rows(records,label):return next(r for r in records if r['label']==label)
def table(header,rr):
    return '| '+' | '.join(header)+' |\n| '+' | '.join(['---']*len(header))+' |\n'+''.join('| '+' | '.join(map(str,r))+' |\n' for r in rr)
def f(x):return f'{x:.2f}'

def run():
    logs=['pilot_gaussian.jsonl','coordinated.jsonl','normalization.jsonl','harmonics.jsonl',
          'joint.jsonl','strong_normalization_joint.jsonl','width_gain_joint.jsonl']
    counts={n:sum(1 for l in open(OUT/n) if l.strip()) for n in logs}
    allrecords=[r for name in ['pilot_confirmed.json','control_confirmed.json','norm_confirmed.json','joint_confirmed.json'] for r in load(name)]
    assert len(allrecords)==19
    byname={r['config']['name']:r for r in allrecords}
    base=byname['reference_p5'];moderate=byname['norm_moderate_p5'];strong=byname['norm_broad_strong_p5']
    # Merge completed transfer checkpoints by configuration and exponent.
    transfers={}
    for name in ['idr_transfer_checkpoint.json','idr_positive_transfer.json','idr_all_transfers.json']:
        for r in load(name):transfers[(r['config']['name'],r['config']['p'])]=r
    assert all((n,p) in transfers for n,p in [('harmonic_reference',2),('norm_broad_strong_p5',2),('joint_final_0',2),('joint_final_1',3)])
    (OUT/'idr_all_transfers.json').write_text(json.dumps(list(transfers.values()),indent=2))
    observed=load('idr_observed_all_diagnostic.json')
    # A single flat table preserves absolute gains as well as normalized gains.
    flat=[]
    for record in allrecords:
        for g in record['groups']:
            flat.append(dict(model=record['config']['name'],condition=g['label'],drift_gain=g['drift'],probe_gain=g['probe'],
                             normalized_probe_gain=g['ratio'],accuracy=g['acc'],observed_accuracy=g['observed_acc'],
                             observed_normalized_probe_gain=g['observed_ratio'],own_scale_ratio=g['own_scale']['ratio']))
    with open(OUT/'confirmed_gains.csv','w') as out:
        w=csv.DictWriter(out,fieldnames=list(flat[0]));w.writeheader();w.writerows(flat)
    for c in [moderate['config'],strong['config']]:
        (OUT/(c['name']+'_calibration.json')).write_text(json.dumps(calibration(c),indent=2))
    labels=['IDS_10','IDS_25','IDS_45','IDQ_120','IDS_155','IDS_170'];angles=[10,25,45,120,155,170]
    modeltable=table(['Probe offset','Observed','Increase-only, p=5','+ moderate normalization','+ stronger normalization'],
       [[str(a)+'°',f(rows(base['groups'],l)['observed_ratio']),f(rows(base['groups'],l)['ratio']),
         f(rows(moderate['groups'],l)['ratio']),f(rows(strong['groups'],l)['ratio'])] for a,l in zip(angles,labels)])
    accuracytable=table(['Offset','Observed','Increase-only','Moderate norm.','Stronger norm.'],
       [[str(a)+'°',f(rows(base['groups'],l)['observed_acc']),f(rows(base['groups'],l)['acc']),
         f(rows(moderate['groups'],l)['acc']),f(rows(strong['groups'],l)['acc'])] for a,l in zip(angles,labels)])
    ablations=['pilot_6','reference_p5','reference_b12','reference_gp','reference_frozen_variance','reference_width375']
    ablnames=['Reference: p=∞, σ=25°, pure first harmonic','Change p to 5','Change variance exponent to 1.2',
              'Replace harmonic with smooth GP','Hold variance amplitude at prestep value','Change σ to 37.5°']
    ablationtable=table(['Comparison','120° gain','170° gain'],[[n,f(rows(byname[k]['groups'],'IDQ_120')['ratio']),
                       f(rows(byname[k]['groups'],'IDS_170')['ratio'])] for k,n in zip(ablations,ablnames)])
    idrbase=transfers[('harmonic_reference',2)]['groups'];idrnorm=transfers[('norm_broad_strong_p5',2)]['groups']
    idrtable=table(['IDR offset','Observed window diagnostic','Same MT/noise, p=2','Stronger normalization, p=2'],
       [[str(a)+'°',f(next(v['ratio'] for v in observed if v['offset']==a)),
         f(next(v['ratio'] for v in idrbase if v['offset']==a)),f(next(v['ratio'] for v in idrnorm if v['offset']==a))] for a in [10,25,45,90,135,179,180]])
    idj=load('idj_harmonic_reference.json');idjp5=next(v for v in idj if v['config']['p']==5)['groups']
    idjb3=load('idj_joint_final_0.json')[0]['groups']
    idjtable=table(['Largest allocation fraction','Observed accuracy','Reference p=5','Stronger variance scaling'],
       [[f(a['allocation']),f(a['observed_accuracy']),f(a['predicted_accuracy']),f(b['predicted_accuracy'])] for a,b in zip(idjp5,idjb3)])
    numerical=load('separated_resolution_checks.json');maxgrid=max(v['common_noise_grid']['max_pc_difference'] for v in numerical)
    maxmc=max(abs(v['independent_512_draws']['mean_pc_difference']) for v in numerical)
    rms=[v['independent_512_draws']['rmse'] for v in numerical]
    independent=load('independent_full_data.json') if (OUT/'independent_full_data.json').exists() else []
    replication=''
    if independent:
        comparisons=[]
        for r in independent:
            old=byname[r['config']['name']]
            for l in labels:
                comparisons.append(abs(rows(old['groups'],l)['ratio']-rows(r['groups'],l)['ratio']))
        replication=f'Two main p=5 candidates were additionally repeated on all 23,038 trials with an independent seed and 256 draws. The largest change in normalized gain across the six target conditions was {max(comparisons):.3f}. Full replication results are retained in independent_full_data.json.'
    # Figure: actual gain recovery, within-model IDR exponent change, joint search.
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140})
    fig,axs=plt.subplots(1,3,figsize=(15,4.7));colors=['#176eaa','#21a179','#d77a21']
    ax=axs[0]
    ax.plot(angles,[rows(base['groups'],l)['observed_ratio'] for l in labels],'ko',label='Observed (120°: IDQ)',zorder=5)
    for rec,color,name in zip([base,moderate,strong],colors,['Increase-only, p=5','+ moderate normalization','+ stronger normalization']):
        ax.plot(angles,[rows(rec['groups'],l)['ratio'] for l in labels],'-o',color=color,ms=4,label=name)
    ax.axhline(0,color='.7',lw=.8);ax.set(xlabel='Probe offset (degrees)',ylabel='Probe gain / drift gain',title='A   IDS/IDQ: useful partial matches',ylim=(-.35,1.15),xticks=[0,45,90,135,180]);ax.legend(fontsize=8,loc='upper right')
    ax=axs[1];aa=[v['offset'] for v in observed]
    ax.plot(aa,[v['ratio'] for v in observed],'ko',label='Observed window diagnostic',zorder=5)
    for p,ls in [(5,'--'),(2,'-')]:
        v=transfers[('harmonic_reference',p)]['groups'];ax.plot([r['offset'] for r in v],[r['ratio'] for r in v],ls+'o',color=colors[0],ms=4,label=f'Same MT/noise, p={p}')
    ax.axhline(0,color='.7',lw=.8);ax.set(xlabel='Probe offset (degrees)',ylabel='Probe / drift slope',title='B   IDR: lowering p misses the baseline',ylim=(-.4,1.15),xticks=[0,45,90,135,180]);ax.legend(fontsize=8)
    ax=axs[2];joint=[json.loads(l) for name in ['joint.jsonl','strong_normalization_joint.jsonl','width_gain_joint.jsonl'] for l in open(OUT/name)]
    xx=[];yy=[];good=[]
    for r in joint:
        gg={g['label']:g for g in r['groups']};jj={g['offset']:g for g in r['idr']}
        x=gg['IDS_170']['ratio'];y=jj[90]['ratio']
        if not np.isfinite(x+y):continue
        xx.append(x);yy.append(y)
        good.append(all(.65<=g['acc']<=.9 for g in gg.values() if g['label']!='IDS_120') and all(.65<=g['acc']<=.9 for g in jj.values()) and
                    .75<gg['IDS_10']['ratio']<1.25 and -.25<gg['IDQ_120']['ratio']<.1 and -.05<gg['IDS_45']['ratio']<.3)
    xx=np.array(xx);yy=np.array(yy);good=np.array(good)
    ax.scatter(xx,yy,s=9,color='.75',alpha=.35,label='Joint-search settings');ax.scatter(xx[good],yy[good],s=15,color=colors[0],alpha=.7,label='Reasonable accuracy + near/120° gains')
    ax.plot([.718],[.289],'k*',ms=11,label='Observed diagnostic targets');ax.axhline(0,color='.7',lw=.8);ax.axvline(0,color='.7',lw=.8)
    ax.set(xlabel='IDS 170° normalized gain',ylabel='IDR 90° slope ratio',title='C   Remaining tradeoff (coarse search)',xlim=(-.65,1.05),ylim=(-.7,.8));ax.legend(fontsize=7.5,loc='upper left')
    fig.text(.01,.015,'All readouts pool only positive pre-to-step changes. IDR uses matched window-average diagnostics; historical IDR kernel fitting was not available.',fontsize=8,color='.3')
    fig.tight_layout(rect=(0,.05,1,1));fig.savefig(OUT/'MTReadout_IncreaseOnly_Search.png',dpi=180);plt.close(fig)
    report=f'''# MT readout: increase-only parameter exploration

28 October 2026 · Meetz · compact exported stimuli and behavioral summaries

The increase-only model can approach the IDS/IDQ gain profile without a negative-change pool. I did not find a shared MT/noise model that also produces the positive middle-angle IDR gains by lowering only the exponent. The clearest obstacle is the IDR response around 90°, rather than the existence of any near-null rebound.

![Simulation results](MTReadout_IncreaseOnly_Search.png)

## What was explored

{sum(counts.values()):,} screening evaluations, including overlapping grids and reference settings; 22 larger-sample finalist screens; {len(allrecords)} all-trial IDS/IDQ confirmations. All-trial confirmations used 23,038 trials, 128 internal-noise draws and a 2° bank. No new raw session files were needed. This is a feasibility search, not a significance test or a claim of exhaustive parameter coverage.

Every model used the full physical dot distribution, equal 500–750 and 750–1000 ms windows, independent internal noise across windows and patches, and the full 2AFC comparison. Its evidence was

`E = p-norm(max(Rstep − Rpre, 0))`.

The signed unit difference is computed before rectification. No negative-change pool, offset threshold, downstream decision noise, or angle-specific exponent was used. A single exponent applies uniformly over directions and patches; IDR can use a lower exponent. Recorded external stimulus fluctuations remain in the inputs.

The search varied mean Gaussian width (18–50°, with a separate 37.5° branch), direction-uniform three-class width mixtures, noise amplitude, angular correlation length (10–90°), global noise, smooth versus mean-balanced fields, and first-harmonic noise. Additional covariance checks used second through fourth harmonics. Width mixtures shared their directional field; their peak versus isotropic-rate scaling was also varied. Width-mixture `width_cv` is a latent spread control; `width_classes()` returns the actual widths and their actual CV.

Noise SD followed `sd × (mean response / 30)^(b/2)`. The main grid used variance exponents b = 0, 1, 1.2, 1.6 and 2. Joint diagnostics extended b to 3 and 4. Gaussian rates were floored at zero; a lognormal alternative preserved the intended mean and variance without a rate floor. The familiar b≈1.2 provides a physiological reference from the supplied Britten et al. (1993), p.1161; it is not treated as an exact bound for these shorter windows. Exponents ranged from 1.5 through 12 plus infinity; none required p=1.

## IDS/IDQ partial matches

All values below are probe gain divided by drift gain, recovered using the original predictors and bounded Weibull gain fit. The 120° target is randomized-drift IDQ. Sparse native IDS120 is retained in the machine results but excluded as a target.

{modeltable}

The three displayed models share a nominal reference σ=25°, SD=5 response units at rate 30, b=2, and a pure first-harmonic noise field. Normalization changes the output tuning as described below. The field correlation is cos(direction separation), including −1 at 180°: a strong assumption needing biological justification. The MT mean response still uses the actual isotropic remainder and coherent dot directions. These are partial fits, not a restoration of the usual homogeneous 37.5° model.

The moderate and stronger normalization models additionally use a 140° suppressive kernel with strength 0.35 or 0.85. Their largest gains require less suppression from remote directions. The moderate case improves the rebound while leaving 45° close to zero. The stronger case reaches about 0.42 at 170° but makes 45° negative. Increasing p toward a hard max brings relatively little improvement.

Performance remains uneven:

{accuracytable}

In particular, the near-drift conditions are underpredicted. Independent refits allowing an overall psychometric scale preserve the main gain effects; absolute gains and these refits are in `confirmed_gains.csv` and the JSON results.

A separate standardized endpoint check reassigned probe directions on the same 1,500 IDS10 trials. The baseline p=5 model gave approximately 1.00 at 0° and 0.15 at 180°; the stronger normalization p=5 model gave 1.01 and 0.36. Those are counterfactual stimulus checks, not observed IDS180 fits and not points to splice into the cross-session curve above.

## What drives the rebound

These comparisons each change one item from the hard-max, narrow-width reference:

{ablationtable}

Narrower tuning and opposing noise fluctuations matter most in this comparison. Response-dependent noise amplitude contributes but is not the entire effect: freezing it at prestep values leaves a smaller positive rebound. Noise changes which units contribute positive evidence; there is no explicit decrease penalty.

Width heterogeneity can sharpen the near-drift profile while preserving a 37.5° mean, but did not solve the full profile. Changing peak versus isotropic-rate scaling across width classes did not remove the joint IDR obstacle. The oscillatory higher-harmonic covariances also failed to supply a joint solution in the tested settings.

## IDR is the stronger constraint

The historical IDR kernel-gain source was absent from the export. I therefore fitted the same window-average probit diagnostic to observed choices and simulated probabilities, conditioning on nominal step. These diagnostic ratios must not be equated with the historical kernel-weighted estimates. Observed values below use all 52,163 noisy increment trials; simulations use up to 1,500 per offset (10,293 total). The 179° yoked pair and 180° single probe were preserved; those adjacent points have different physical stimulus constructions.

{idrtable}

For the unnormalized candidate, decreasing p from 5 to 2 changes 90° from about −0.19 to −0.22, and 135° from −0.12 to −0.15. This goes in the wrong direction for the positive IDR baseline. Stronger broad normalization creates a large near-null rise and some gain at 135°, but still misses 90° and can overshoot the 180° gain.

The joint search included IDR conditions from the outset and extended to stronger noise scaling and strong local normalization. None of 2,080 joint screening settings passed a loose descriptive combination of reasonable performance, near-drift/120°/45° IDS gains, positive 170° gain, and positive middle-angle IDR gains. This finite screen is evidence about these implemented families, not an impossibility result.

## IDJ catches an additional tradeoff

The MT/noise parameters were transferred without a new noise fit. Actual allocations, external noise and dot budgets were preserved:

{idjtable}

The reference has a modest allocation cost. A stronger-noise-scaling candidate (σ=18°, b=3, SD=7.5, p=5) raises the full-data IDS170 gain to about 0.49, but its IDQ120 gain becomes positive (~0.16), and IDJ performance improves with divided allocation. That reverses the observed allocation cost. Extra variance nonlinearity does not provide a clean rescue.

## Directional normalization: useful but constrained

This was one explicit mean-response normalization family, using a Gaussian-weighted sum of the actual dot distribution in the denominator. It operates on physical directions, so coincident streams are equivalent to one combined stream. Noise is added at the normalized MT output, not inside a simulated recurrent normalizer.

For each strength, the excitatory kernel and response gain were calibrated to preserve the isotropic baseline and preferred response at 15% coherence. Other responses were checked rather than assumed to remain linear. For the displayed moderate and stronger models, the maximum single-direction deviations from the unnormalized preferred coherence curve are approximately 7.6% and 16.3% of its full 0–100% response excursion; corresponding null deviations are 4.2% and 11.7%. Up to 30% coherence the largest absolute differences are about 0.90 and 2.14 response units. These are calibration costs to evaluate biologically, not numerical fit tolerances inferred from Britten et al.

## Numerical and scope checks

Rotation invariance, coincident-stream invariance, agreement with the original homogeneous implementation, noise moments, and patch/window independence passed. With exactly matched first-harmonic noise, switching 2° to 1° changes individual trial probability by at most {maxgrid:.6f} in the tested finalists. Independent 512-draw runs change mean probability by at most {maxmc:.6f} on the matched 2,100-trial check set, with per-trial RMS differences of {min(rms):.3f}–{max(rms):.3f} relative to 64 draws. Small-subset fitted gains can consequently shift by around 0.1. {replication} These checks support numerical stability; they do not validate a biological mechanism or guarantee that the parameter search found every solution.

Unexplored mechanisms include independent noise among additional cell subtypes, more realistic stimulus-dependent covariance, recurrent normalization with noise inside its denominator, and richer prestep adaptation. All present results use the specified window-average representation.

## Recommendation

Keep the positive-change readout as a serious candidate: it can supply a substantial near-null rebound without explicitly pooling decreases. Retain the moderate-normalization p=5 model as a concrete reference, with the stronger case illustrating the tradeoff. Further broad tuning of the same noise parameters now has lower value. The next mechanism should first explain how oblique probes create positive IDR evidence while preserving single-direction coherence calibration, IDS120, and the sign of the IDJ allocation cost. Prestep-dependent gain control or recurrent competition could be tested for that specific job. This round does not establish that either will work.

The code, all screening records, finalist predictions, absolute gains, endpoint checks, task transfers, and numerical checks accompany this report. The compact input ZIP is not duplicated.
'''
    (OUT/'MTReadout_IncreaseOnly_Report.md').write_text(report)
    summary=dict(date='2026-10-28',screening_counts=counts,screening_evaluations=sum(counts.values()),
                 larger_screens=22,all_trial_confirmations=len(allrecords),all_trial_n=23038,
                 input_sha256=hashlib.sha256(BUNDLE.read_bytes()).hexdigest(),python=sys.version,
                 numpy=np.__version__,scipy=scipy.__version__,matplotlib=matplotlib.__version__,
                 numerical_validation=load('expanded_validation.json'),main_model_names=['reference_p5','norm_moderate_p5','norm_broad_strong_p5'])
    (OUT/'run_summary.json').write_text(json.dumps(summary,indent=2))
    readme='''# MT readout increase-only exploration, 28 October 2026

Start with positive-search/MTReadout_IncreaseOnly_Report.md and its PNG figure.
confirmed_gains.csv preserves both absolute and normalized gains for all 19 all-trial confirmations.
This Python research package consumes the existing compact four-domain export; it does not open raw ProbeSessions or alter the MATLAB file-selection workflow.

Reproduce using Python with NumPy, SciPy and Matplotlib. Place MTReadoutExplorationInputs.zip next to research.py, or set MT_INPUT_BUNDLE to its absolute path. First run `python research.py` once to validate the bundle and build the derived trial cache before starting parallel searches.

Commands, in order:
```
python positive_check.py
python positive_search.py --power 10 --workers 3
python positive_normalization.py
python positive_harmonics.py
python positive_joint.py
python positive_strongnorm.py
python positive_widthgain.py
python positive_followup.py screen screen_configs.json --outfile screened.json --workers 3
python positive_followup.py confirm pilot_configs.json --outfile pilot_confirmed.json --workers 2
python positive_followup.py confirm control_configs.json --outfile control_confirmed.json --workers 2
python positive_followup.py confirm norm_confirm_configs.json --outfile norm_confirmed.json --workers 2
python positive_followup.py confirm joint_final_configs.json --outfile joint_confirmed.json --workers 2
python positive_followup.py transfer transfer_initial_config.json
python positive_followup.py transfer transfer_norm_config.json
python positive_followup.py transfer joint_final_configs.json
python positive_resolution.py
python positive_independent.py
python positive_report.py
```
Transfer commands should run sequentially because they maintain an aggregate file. Configuration JSON files are under positive-search; the CLI resolves those filenames there. Infinite p is represented by Python JSON's Infinity token. Original stochastic seeds, sample sizes, and direction spacings are in the scripts. The coarse pilot predates the lognormal branch and is retained as a search record; the full main grid supersedes it.

Large derived trial/temporal caches and the original input ZIP are excluded. Output NPZ files contain global exported-trial indices and predicted probabilities; run research.py to regenerate their trial identities. Source MATLAB gain-analysis snapshots are in the original compact input ZIP. Historical IDR kernel-gain code was not present, so IDR transfers are matched window-average diagnostics.

Endpoint records use IDS10 stimuli with reassigned directions and are counterfactual checks. All-trial confirmations use the original physical stimuli. The report figure does not splice synthetic endpoints into empirical offset curves.
'''
    (ROOT/'README_IncreaseOnly.md').write_text(readme)
    print('screening counts',counts,'total',sum(counts.values()),'confirmed',len(allrecords))
    print(OUT/'MTReadout_IncreaseOnly_Report.md')

if __name__=='__main__':run()
