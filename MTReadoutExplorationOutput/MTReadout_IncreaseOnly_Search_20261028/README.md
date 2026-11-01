# MT readout increase-only exploration, 28 October 2026

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
