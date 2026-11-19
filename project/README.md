# Hjort et al. 2026 meta-RPE model applied to the Dec-2025 opioid cohort (2026-11-17)

Paper: Hjort, ... Stuber (2026) Prefrontal to ventral tegmental area dynamics drive contingency degradation.
Nature 655:174. doi:10.1038/s41586-026-10443-5. Code: github.com/stuberlab/Hjort-et-al.-2026-PFC-and-reversal-learning
(commit e06ae22, MIT licence; their 40-mouse sample data copied to `data/hjort2026_sample_choice/`).

## Scripts (run in order; raw data on K: read only)
| script | does | main output |
|---|---|---|
| `m01_parse_trials.py` | raw controller logs -> one row per PR trial, every lick and drop | `data/trials_pr.csv`, `data/licks_all.csv` (rewards = raw table, max diff 0) |
| `m02_hjort_models.py` | the paper's 6 models, ported line by line (numba) | module |
| `m03_validate_on_paper_data.py` | refit the authors' 40 mice | `data/m03_paper_data_summary.json` (reproduces Fig. 1k,l) |
| `m04_fit_user_data.py` | same learning rules on our data, per mouse | `data/m04_fits.csv`, `m04_summary.json` |
| `m05_recovery_and_group_fits.py` | model recovery; group fits with shared learning rates; bootstrap | `data/m05_recovery.json`, `m05_group_fits.json` |
| `m06_modelfree_switch.py` | change in licking on day 14 / day 17 by 3-min block | `data/m06_modelfree.json` |
| `m07_fig1_paper_explained.py` | Fig. 1 of the paper redrawn from the authors' data | `figs/fig1_paper_fig1_explained.png` |
| `m08_bootstrap_asymmetry.py` | RPE 2-alpha group fit, 1000 bootstraps over mice | `data/m08_bootstrap_asymmetry.json` |
| `m09_figures_user_data.py` | figures 2-4 | `figs/fig2...png`, `fig3...png`, `fig4...png` |
| `m13_trial_by_trial.py` | (18 Nov) trial-by-trial checks as in the paper's Fig. 1d,e,m: within day 14 vs overnight, sigmoid midpoints, start vs models, posterior predictive check | `data/m13_summary.json`, `m13_per_mouse.csv`, `m13_bins.csv`, `m13_log.txt` |
| `m14_make_trial_slides.py` | (18 Nov) two new slides 7-8 (EN + KR) + v3 decks (v1 of slides 7-8, superseded) | `slides_20261118_trial_by_trial\` |
| `m15_trial_checks_v2.py` | (18 Nov, after the review) noise-matched model check (2 calibrations), placebo day pair, mouse-level t, equal windows, leave-one-out | `data/m15_summary.json`, `m15_log.txt` |
| `m16_make_trial_slides_v2.py` | (18 Nov) slides 7-8 v2 (EN + KR) + v4 decks | `slides_20261118_trial_by_trial\` (`*_v2_*`, `*_v4_*`) |

## Slides (2026-11-17, `m10_make_slides.py`, python-pptx only, lab template)
`slides_20261117\` - `mrpe_hjort2026_opioid_6slides_EN_20261117.pptx` (English) and `..._KR_20261117.pptx`
(Korean slide text + a full Korean explanation of each slide in the speaker notes; figures keep English labels).
1 paper Fig. 1 (authors' data, our port) | 2 mapping table + timeline | 3 behaviour at the two switches |
4 group model fits | 5 per-mouse recovery + learning asymmetry | 6 model variables per drop (neural regressors).
Every printed number is asserted against the JSON outputs and every text box is measured (fit_report.json);
previews in `preview_EN\` / `preview_KR\` (LibreOffice headless, own profile); one-slide files in `single_slides\`.

### Slides v2 (slides 1-6 still current; the CURRENT 8-slide decks are v4 = these 6 + trial-by-trial v2 slides 7-8; after one independent review; `m11_extra_checks.py` + `m12_make_slides_v2.py`)
`slides_20261117\mrpe_hjort2026_opioid_6slides_KR_v2_20261117.pptx` / `..._EN_v2_20261117.pptx`, previews in
`preview_KR_v2\` / `preview_EN_v2\` (v1 kept). Review fixes: whole-session switch tests with Active-vs-Passive Welch
test and CI (no post-hoc 9-min window); +0.32 -> +0.31 (double rounding); titles state model fits, not mechanisms;
leave-one-mouse-out CV and group-level recovery added; Passive beta < 0 stated; meta-RPE prediction framed as
conditional; KR wording (피팅, 마우스) and footers match EN. 291 checks pass.
New numbers (m11): day 14 vs 13 Active -0.06 [-0.27, 0.15] (p = 0.53), Passive +0.40 (p = 0.0002), Welch p = 0.001;
Active days 15/16 -0.27 / -0.35 (p = 0.018 / 0.015); day 17 Active +0.35 (p = 0.014); asymmetry 0.18 vs 0.35, p = 0.21.
LOMO-CV Active: RPE 2a best, meta-RPE +1.3, RPE +8.3, instant +15.5 negLL; Passive all within 2.6 of instant switch.
Group recovery (Active): RPE 2a never picked for RPE- or instant-made data; meta-RPE recovered 16%.

### Trial by trial v2 (CURRENT, 18 Nov 2026; after one independent review; `m15` + `m16`)
The review reproduced every m13 number from the raw lick table; the problems were inferential. v1 (below) is kept.
- m13's model check drew every 30-s bin independently from the fitted p: real licking varies ~2x more from day to day
  and is autocorrelated, so its p values were too small (not 'conservative'); it even rejected a no-change placebo
  (Passive day 12 -> 13, p < 0.001). Now: mouse-level paired t (mice minus model, n = 6) as the main test, and a
  noise-matched simulation (mouse-day offset + AR(1) logit noise fitted to day-to-day residual differences on days
  11-13; two calibrations as a range; placebo passes for Active).
- Day 14 (Active): change -0.06 vs every model predicting more (mice minus model +0.05 to +0.17). Mouse-level p:
  instant switch 0.073, RPE 0.15, RPE 2a 0.30, meta-RPE 0.58; simulations: instant 0.021-0.050, RPE 0.078-0.14,
  RPE 2a 0.23-0.31, meta-RPE 0.52-0.59. No model is robustly rejected; only instant switch is near rejection.
  Robust observation: no clear fall in the first 7.5 min (+0.02) although 88% of the day's drops were earned by then
  (drops are front-loaded on every PR day: day 13 58% by 2.5 min; first 2 min at ceiling, P(lick) = 1).
- Within day 14 (-0.16, 95% CI -0.60 to +0.29) vs overnight (-0.02, -0.53 to +0.49): cannot be separated
  (mirror images across mice; 5-min windows split it 50/50). v1's 'no extra change overnight' is withdrawn.
- Timing: group-curve midpoint 36.5 vs 4.5 (later in > 99% of bootstraps; equal 60-trial windows 49.0 vs 4.5;
  leave-one-out 31.5-44 vs 4.0-4.5), but per mouse p = 0.11 (60-trial windows 0.12), n = 6. Start 14.4 (95% CI
  0-47, leave-one-out 6-42) excludes no model (noise-matched p >= 0.068). v1's 'rules out instant switch' is withdrawn.
- v1 notes said slide 8a was normalised: wrong (raw change in P(lick)); fixed in v2. 'r2' in m13 is uncentred.
Slides v2: `slides_20261118_trial_by_trial\mrpe_trial_by_trial_2slides_{EN,KR}_v2_20261118.pptx`, `single_slides_v2\`,
previews `preview_v2_{EN,KR}\`; decks `mrpe_hjort2026_opioid_8slides_{EN,KR}_v4_20261118.pptx` (v2 slides 1-6
byte-identical + slides 7-8 v2). 134 checks pass (`fit_report_v2.json`).

### Trial by trial v1 (18 Nov 2026; `m13_trial_by_trial.py` + `m14_make_trial_slides.py`) - SUPERSEDED by v2 above
Kept for the record; the overnight, 'rules out' and 'conservative' statements below were overstated (see v2).
Asked: can we do the paper's trial-level analyses? Trial = the same 30-s bins the models were fitted on; each bin is
compared with the same bin on the day before the switch (day 13 for drug -> water, day 16 for water -> drug).
Models = the m05 group fits (shared learning rates, readout refitted per mouse). Model check = posterior predictive,
2000 simulated data sets per model, not refitted (conservative).
- (A) Day 14: Active drank 66% of the day's 28.3 water drops in the first 2.5 min (88% by 7.5 min), yet licking was
  unchanged in the first half (+0.02) and lower in the second half (-0.13; within-day change -0.16, p = 0.41);
  day 14 -> 15 overnight -0.02 (p = 0.92): no extra overnight step. Whole day 14 -0.06 vs model predictions
  -0.22 (instant switch, p < 0.001), -0.19 (RPE, p = 0.005), -0.15 (RPE 2a, p = 0.050), -0.11 (meta-RPE, p = 0.28).
  Passive changed at once (day 14 +0.40, first half +0.31). Lick counts give the same picture.
- (B) Fig. 1d,e: Active group-curve midpoint trial 36.5 after drug -> water (95% CI 13.5-68.0) vs 4.5 after
  water -> drug (3.5-15.0); later in 99.7% of bootstraps; per mouse 38.6 vs 12.4, paired t p = 0.11 (n = 6).
- (C) Fig. 1m: start of the drug -> water change trial 14.4; vs simulated data: instant switch p = 0.011, RPE 0.054,
  RPE 2a 0.11, meta-RPE 0.79. Paper-style noise-free curves all start by trial 1.8 vs mice 37.9 (Holm p = 0.072),
  but that comparison is biased by noise in yes/no licking.
- Meaning: the paper's signature (no change when the errors are largest, change later) is present; the delay rules
  out an instant switch and one learning rate, is borderline for two fixed rates and fits meta-RPE. Not proof of
  error averaging: n = 6, one switch each way, and a delayed oral drug effect gives the same lag.
Slides: `slides_20261118_trial_by_trial\mrpe_trial_by_trial_2slides_{EN,KR}_20261118.pptx` (the two new slides only),
`single_slides\`, previews `preview_{EN,KR}\`; v3 decks `mrpe_hjort2026_opioid_8slides_{EN,KR}_v3_20261118.pptx` =
v2 slides 1-6 (byte-identical parts, checked) + slides 7-8. 107 checks pass (fit_report.json).

## Mapping (paper -> our task)
trial (cue, ~20 s) -> 30-s bin of the 15-min PR session; choice (anticipatory lick) -> any lick in the bin;
outcome (sucrose 1 / omission 0) -> each drop: morphine 1 / water 0; value -> expected drug from the spout
(starts 0 on day 3, updated by every drop the mouse drank, Passive replay drops included, carried over days);
reversal -> day 14 morphine -> water, day 17 water -> morphine. Likelihood on days 11-18 (8 x 30 bins per mouse).
Readout for every model alike: logistic(b0 + beta*Q + time in session + first-3-min term), beta may be negative.
Learning-rule bounds = the paper's. mRPE step alpha*|M| capped at 1 (never reached in any fit).

## Results
- Port check: authors' 40 mice, BIC change from best: meta-RPE 2.6 vs RPE 5.7, RPE 2a 6.2, RPE t-1 6.2, PH 12.4,
  eMack 18.6; meta-RPE better than each (paired t, Holm p <= 0.010). negLL: meta-RPE best in 30/40.
- Model-free (our data): day 14 first 9 min, change vs day 13: Active 0.00 (p = 1.00), Passive +0.35 (p = 0.005).
  Day 17 first 9 min: Active +0.32 (p = 0.038), Passive -0.25 (p = 0.093).
- Per-mouse model comparison (the paper's method) has no power here: data simulated from each mouse's fitted
  meta-RPE are picked as meta-RPE by BIC in 2% (Active) / 1% (Passive) of sets.
- Group fits (shared learning rates): Active best = RPE 2a (alpha- 0.01 [0.01, 0.03], alpha+ 0.12 [0.05, 1.00];
  alpha- < alpha+ in 99.5% of bootstraps); meta-RPE dBIC +9 (d = 0.01 < e = 0.10; d < e in 98.3%).
  Passive best = today's outcome (instant switch); meta-RPE dBIC +20.

## Caveats
Different task (free-operant PR, outcome identity not probability); in-session drug effect and withdrawal state can
mimic fast enhancement / slow degradation; 6 Active + 7 Passive mice with one switch each way; Passive instant
switch may be taste (morphine is bitter), not learning; 6099_white (died day 13) left out of comparisons.
