"""Render the completed first exploration batch; no new fits are performed."""
from pathlib import Path
import hashlib,json,csv,zipfile
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import erf
from research import ROOT,BUNDLE
from search import pool

def main():
    original=json.loads((ROOT/'confirmed.json').read_text())[0]
    alt=json.loads((ROOT/'noise_confirmed.json').read_text())[0]
    constrained=json.loads((ROOT/'noise_confirmed.json').read_text())[1]
    observed={g['label']:g for g in json.loads((ROOT/'observed_gains.json').read_text())}
    labels=['IDS_10','IDS_25','IDS_45','IDQ_120','IDS_155','IDS_170'];angles=[10,25,45,120,155,170]
    def rows(model):return {r['label']:r for r in model['groups']}
    base=rows(original);alternative=rows(alt);cp5=rows(constrained)
    fig,axes=plt.subplots(1,3,figsize=(13.5,4.4),layout='constrained')
    colors=['#2274a5','#bd5319','#6650a4']
    for ax,key,ylabel in [(axes[0],'ratio','Probe / drift gain'),(axes[1],'acc','Fraction correct')]:
        oy=[observed[l]['ratio'] if key=='ratio' else base[l]['observed_acc'] for l in labels]
        ax.plot(angles,oy,'ko-',lw=1.6,ms=5,label='Observed')
        for data,color,label in [(base,colors[0],'Starting MT p=5'),(cp5,colors[2],'p=5, changed noise/readout'),
                                 (alternative,colors[1],'Exploratory alternative')]:
            ax.plot(angles,[data[l][key] for l in labels],'o--',color=color,lw=1.3,ms=4,label=label)
        ax.set(xlabel='Probe offset (degrees)',ylabel=ylabel,xticks=[0,45,90,120,180],xlim=(0,180))
        ax.spines[['top','right']].set_visible(False);ax.grid(alpha=.16)
    axes[0].axhline(0,color='gray',lw=.7);axes[0].set_ylim(-.55,1.15)
    axes[1].set_ylim(.48,.9);axes[1].axhline(.5,color='gray',lw=.7)
    axes[0].set_title('A. Gain recovery: original predictors',loc='left',fontsize=11)
    axes[1].set_title('B. Same noisy increment trials',loc='left',fontsize=11)
    axes[0].legend(fontsize=8,loc='lower left',frameon=False)
    # Illustrative noiseless calculation, not a behavioral exponent estimate.
    th=np.arange(360);directions=np.array([0,120,-120]);sigma=37.5
    gauss=np.exp(-.5*(((th[:,None]-directions+180)%360-180)/sigma)**2)
    gbar=np.sqrt(2*np.pi)*sigma/360*erf(180/(np.sqrt(2)*sigma))
    pedestal=np.ones(3)*7.071;pre=gauss@pedestal+(100-pedestal.sum())*gbar
    allocations=np.array([[1,0,0],[.7,.15,.15],[.5,.25,.25],[1/3]*3])
    d=[]
    for a in allocations:
        c=pedestal+15*a;post=gauss@c+(100-c.sum())*gbar
        d.append(float(pool(post[None],5)[0]-pool(pre[None],5)[0]))
    relative=np.array(d)/d[0];stream=(allocations**5).sum(1)**.2
    diagnostic=dict(allocations=allocations.tolist(),total_increment_pc=15,pedestal_pc=7.071,
                    population_delta=d,relative_population_delta=relative.tolist(),stream_p5_norm=stream.tolist())
    (ROOT/'allocation_diagnostic.json').write_text(json.dumps(diagnostic,indent=2))
    axes[2].plot(range(4),stream,'o-',color='black',label='p=5 on stimulus streams')
    axes[2].plot(range(4),relative,'o-',color=colors[0],label='p=5 on full MT responses')
    axes[2].set(xticks=range(4),xticklabels=['100/0/0','70/15/15','50/25/25','Equal thirds'],
                xlabel='Allocation of a fixed total increment',ylabel='Evidence relative to single stream',ylim=(0,1.08))
    axes[2].tick_params(axis='x',labelsize=8);axes[2].spines[['top','right']].set_visible(False)
    axes[2].grid(alpha=.16);axes[2].legend(fontsize=8,frameon=False)
    axes[2].set_title('C. Exponents do not transfer directly',loc='left',fontsize=11)
    fig.suptitle('MT readout: first population-model search',fontsize=14)
    fig.supxlabel('Panels A/B: 120° uses randomized IDQ; other offsets use IDS. Altered-noise curves use positive-minus-negative unit-change pools.',fontsize=8)
    fig.savefig(ROOT/'MTReadout_FirstSearch.png',dpi=180)
    fig.savefig(ROOT/'MTReadout_FirstSearch.pdf')
    table=[]
    for l,a in zip(labels,angles):
        table.append(dict(condition=l,offset=a,observed_ratio=observed[l]['ratio'],
                          original_p5_ratio=base[l]['ratio'],alternative_ratio=alternative[l]['ratio'],
                          observed_accuracy=base[l]['observed_acc'],original_p5_accuracy=base[l]['acc'],
                          alternative_accuracy=alternative[l]['acc']))
    with open(ROOT/'main_comparison.csv','w') as f:
        w=csv.DictWriter(f,fieldnames=list(table[0]));w.writeheader();w.writerows(table)
    coverage=(ROOT/'Coverage.csv').read_text().strip()
    comparison='\n'.join(f"| {r['condition']} | {r['observed_ratio']:.2f} | {r['original_p5_ratio']:.2f} | {r['alternative_ratio']:.2f} |" for r in table)
    idj=json.loads((ROOT/'idj_transfer.json').read_text())
    idjbase=next(r for r in idj if r['config']['family']=='pooled_delta' and r['config']['p']==5)
    jt='\n'.join(f"| {r['allocation']:.3f} | {r['n']} | {r['observed_accuracy']:.3f} | {r['predicted_accuracy']:.3f} |" for r in idjbase['groups'])
    digest=hashlib.sha256(BUNDLE.read_bytes()).hexdigest()
    text=f'''# MT readout: first exploration batch — 25 October 2026

The export is sufficient for local exploration. No joint solution was found in this first bounded search. The starting p=5 population model misses the far-offset rebound and predicts an excessive cost of divided signal allocation. Several modifications improve individual features, but no confirmed setting satisfies the whole pattern.

## Data and validation

All 190 exported sources passed numerical validation, including reconstruction of both 250-ms windows, trial identity, 2AFC side/outcome consistency, direction reordering, IDR yoking, and top-up exclusion. There are 147,265 valid trials and 82,143 noisy increment trials. All sources are Meetz.

```csv
{coverage}
```

Random-base-direction IDQ supplies 7,124 noisy increment trials at 120 degrees. The fixed-base IDQ sources were deliberately excluded. One IDR source, `IDReadout_Meetz_20260305_Probe45.mat`, failed export because its noise/time event counts differ; this does not block the current analysis. No exported window has materially negative physical coherence or exceeds the 100% coherence budget.

The experimental gain fits were reproduced from the saved predictors and Weibull calibration, to approximately 0.000002 in the gain coefficients. Native IDS uses its original change-RDK drift and pooled probe predictors. Random IDQ uses the saved four-predictor model (change/no-change RDK × drift/probe), retaining the distinction between RDK roles. No new behavioral confidence intervals or significance tests were computed.

The native IDS 120-degree estimate has only one source/session and a very broad interval. It is retained as a diagnostic but is not a search target. The randomized IDQ summary is the principal target at 120 degrees. The uploaded native IDS data have 10 and 170 degrees as the nearest sampled endpoints; these are not measurements at exactly 0 and 180 degrees.

## Starting model

For unit i, the mean response is the Gaussian-tuned response to every coherent stream plus the response to the remaining isotropic dots:

`R_i = (100 - sum(c)) * mean(G) + sum_k c_k * G_i(direction_k)`.

Gaussian sigma is initially 37.5 degrees, with a dense circular bank. Independent internal Gaussian fields are added in each RDK and each time window, with correlation over preferred direction. Absolute firing rates are floored at zero. The starting evidence for each RDK is `norm_p(R_step) - norm_p(R_pre)`, retaining its sign. The decision selects the RDK with larger evidence. Both windows last 250 ms: 500–750 and 750–1000 ms.

The confirmed starting p=5 run uses the earlier IDS-calibrated noise SD 5.29842, a 30-degree correlation length, and zero spontaneous offset. This is a reference replication, not a newly optimized fit to all four tasks.

## What was searched

The initial grids evaluated **2,006 settings** against the IDS/random-IDQ gains and accuracy. Every exponent was uniform over direction, RDK, and window within a model. No angle-specific exponents or temporally correlated noise between the two windows were introduced. Nothing below p=2 was needed for this batch.

- 900 settings: p=2,3,5,8,infinity; tuning sigma=25,37.5,50 degrees; directional correlation lengths=15,30,60 degrees; internal SD=2.5,5,10,20. Readouts were pooled signed changes, rectified unit changes, positive-minus-negative change pools, half-strength unit reference subtraction, and a post-only control.
- 864 settings: internal SD=1,2.5,5,10; p=3,5,infinity; sigma=25 or37.5. Compared ordinary smooth fields, fields with their uniform mode removed, and fields mixing smooth fluctuations with a first-harmonic opponent component. Noise SD was constant, proportional to sqrt(response), or proportional to response. Also tested squaring the pooled response before pre/post subtraction.
- 192 settings: a direction-specific noisy prestep response supplied divisive adaptation, at strengths0.3,1,3,10; pooled-change and rectified-unit-change variants.
- 50 settings: directional pooling in each exported50-ms bin followed by temporal averaging, versus averaging MT responses before directional pooling. This is an approximation to online processing; activity within a50-ms bin is not retained in this export. Internal noise was a window-estimation field held across bins within its window, independent between windows.

Search scores combined deviations in probe/drift gain ratios and accuracy, using fixed navigation scales0.15 and0.05 respectively. They are optimization aids, not statistical evidence. Native IDS120 did not contribute to the score.

Coarse runs used up to180 randomly selected trials per gain group, a5-degree bank,32 antithetic draws per patch/window, and all cross-patch draw pairs for 2AFC probabilities. Selected settings were then evaluated on **all23,038 noisy IDS/IDQ trials**, using a2-degree bank and128 draws, with independent Monte Carlo noise. Four initial finalists also received an independent-seed1-degree/256-draw check on1,152 trials. Mean probability differences were below0.002; single-trial RMS differences0.016–0.024 reflect grid and Monte Carlo changes together. The two later noise finalists were checked separately in `noise_resolution_check.json`.

## Gain results

All values below are pooled probe gain divided by drift gain, recovered with the original experimental predictors and calibration.

| Condition | Observed | Starting MT p=5 | Exploratory alternative |
|---|---:|---:|---:|
{comparison}

The alternative uses sigma25 degrees, hard maximum, positive-minus-negative unit-change pools, response-proportional noise SD, and a latent directional field with its uniform mode removed. Its partial far-offset recovery is informative, but its 25/45-degree gains are too large, its 120-degree gain is not negative, and its far-offset accuracy remains too low. It also explicitly uses negative unit changes, departing from the preferred increase-only physiological interpretation. It is a diagnostic candidate, not an accepted explanation.

Keeping p=5 and sigma37.5 in that noise/readout family yields approximately0.16 at170 degrees and−0.28 at120, still far from the combination required. The all-trial rectified-unit-change finalist also fails the far-offset rebound. Changing subtraction/rectification order alone was insufficient in this grid. The tested divisive templates and temporal-pooling order did not solve the joint problem.

Allowing a separate overall Weibull sensitivity scale when recovering each model's gains changes their magnitudes somewhat but does not remove these pattern failures. These checks are saved with the confirmed results under `own_scale`.

![First-search comparison](MTReadout_FirstSearch.png)

## Why population p=5 is not stimulus-stream p=5

With sigma37.5, three streams separated by120 degrees, a7.071% pedestal in each stream, and a15-percentage-point total increment, the noiseless signed population-norm difference for equal thirds is **{relative[-1]:.3f}** of that for a single-stream increment. Applying p=5 directly to the three increments gives **{stream[-1]:.3f}** instead. Rectifying unit changes first gives approximately0.115 in this example.

These are evidence ratios for an illustrative stimulus, not fitted behavioral exponents or predicted accuracy ratios. The full MT response includes coherent-dot replacement of isotropic dots. This replacement introduces a broad negative response component, while the equal pre-step pedestals affect the population norm. Consequently, a p inferred between a few stimulus streams cannot simply be inserted as the exponent over full MT firing rates and expected to preserve its behavioral meaning.

## IDJ and IDR checks

IDJ noise SD was calibrated by choice likelihood on the single-stream allocation trials only, then transferred unchanged to divided-allocation trials. The exponent and other parameters were fixed within each check. A shared SD across sessions is a deliberately restricted calibration; it does not guarantee matching reference accuracy.

| Largest allocation fraction | Trials | Observed accuracy | Starting family, p=5 |
|---|---:|---:|---:|
{jt}

This starting family has too much allocation cost, and even its reference accuracy is low. The alternative noise/readout family reduces the allocation penalty; `idj_alternative.json` records its transfer, including remaining reference mismatch. No joint IDS/IDQ/IDJ fit is claimed.

For IDR, eight starting-family and four alternative-family transfer checks used actual stimulus windows and yoking. Lowering p alone did not supply the desired broadly elevated off-axis influence. **These IDR gain ratios are exploratory probit regressions on uniform step-window change-RDK noise, with nominal-step intercepts.** The export does not include the historical IDR gain-analysis function/summary, so they must not be presented as exact replications of the published/current kernel-weighted IDR gains. Separate results at179 and180 preserve the different physical stream construction. These checks used at most800 trials per offset, not all104,050 IDR trials.

## Interpretation and next priorities

The completed batch provides a bounded negative result for the simple starting implementation, with two constructive findings: the mapping from population p to behavioral p depends strongly on reference/signal representation, and response-dependent directional noise can change the far-offset rebound substantially. Neither finding establishes a biologically acceptable joint account yet.

Further work should compare explicit prestep-template/reference mechanisms and response-dependent directional covariance while keeping the full physical dot construction. IDJ allocation behavior should constrain those mechanisms alongside the gain curves, rather than being checked only after fitting them. Tuning-width heterogeneity and finer searches near partial candidates remain untested. Session-specific nuisance calibration and the original IDR gain-analysis code would improve later comparisons. A failure of the finite grids here is not an impossibility result for MT p-norm readout.

## Reproduction

Python3 with numpy, scipy and matplotlib is sufficient. Place the original uploaded ZIP anywhere and set its path, for example:

```bash
export MT_INPUT_BUNDLE=/absolute/path/MTReadoutExplorationInputs.zip
export OPENBLAS_NUM_THREADS=1
python research.py
python check_forward.py
python search.py
python refine.py finalists
python noise_search.py
python template_search.py
python temporal.py
python refine.py idj
python refine.py idj-alt
python refine.py idr
python refine.py idr-alt
python make_report.py
```

Input summaries are unpacked from the ZIP when first needed. Derived trial/temporal caches are reproducible and omitted from the results package. Original input SHA-256: `{digest}`.

The package contains code, the exported analysis-function snapshots, the validation/source audit, all grid records, confirmed results and predictions, and this figure. It does not duplicate the raw behavioral input ZIP.
'''
    # Repair visual spacing in prose without altering scientific expressions.
    for a,b in [('or37','or 37'),('strengths0','strengths 0'),('exported50','exported 50'),('a50','a 50'),
                ('scales0','scales 0'),('and0','and 0'),('used up to180','used up to 180'),('a5','a 5'),
                ('bank,32','bank, 32'),('all23','all 23'),('a2-degree','a 2-degree'),('and128','and 128'),
                ('independent-seed1','independent-seed 1'),('on1,152','on 1,152'),('below0','below 0'),
                ('differences0','differences 0'),('sigma25','sigma 25'),('25/45','25/45'),('sigma37','sigma 37'),
                ('approximately0','approximately 0'),('at170','at 170'),('and−0','and −0'),('at120','at 120'),
                ('by120','by 120'),('a7.071','a 7.071'),('a15-','a 15-'),('gives approximately0','gives approximately 0'),
                ('at179','at 179'),('and180','and 180'),('most800','most 800'),('all104','all 104')]:text=text.replace(a,b)
    (ROOT/'MTReadout_FirstSearch.md').write_text(text)

if __name__=='__main__':main()
