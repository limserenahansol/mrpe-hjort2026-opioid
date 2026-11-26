# 5. Trial-by-trial analyses (analogues of the paper's Fig. 1d, e, m)

Asked on 18 Nov 2026: "based on the paper, can we also do trial by trial?"

- **v1** (`m13`, slides 7-8 v1) was reviewed once, independently. The review reproduced every number from the raw
  lick table, but found inferential errors.
- **v2** (`m15` + `m16`, slides 7-8 v2) is **CURRENT**. This file describes v2 and says exactly what v1 got wrong
  (section 5.8).

## 5.1 The three questions

- **(A) Day 14.** When does Active licking change after morphine → water: within day 14, or only after a night? The
  paper's point (the pattern behind its Fig. 1h) is that behaviour does *not* change when the errors are largest. The Fig. 1h correlation itself is not computed
  here.
- **(B) Fig. 1d,e.** Is the sigmoid midpoint of the change later after morphine → water than after water → morphine?
- **(C) Fig. 1m.** At which trial does the morphine → water change start, and do the models start at the same trial?

## 5.2 Definitions

- **Trial** = one 30-s bin. These are the same bins the models were fitted on.
- **x** = trials from the switch:
  - morphine → water: x = 0-89 (days 14, 15, 16 × 30 bins);
  - water → morphine: x = 0-59 (days 17, 18).
- **D(day, bin) = lick(day, bin) - lick(reference day, same bin)**:
  - the reference day is **day 13** (last morphine day) for morphine → water, and **day 16** (last water day) for
    water → morphine;
  - comparing the same bin removes the fall in licking within every session.
- **Model predictions:**
  - learning rates come from the group fits (`m05`);
  - the readout is refitted per mouse on days 11-18 (`logit_fit`: Newton iterations with a 1e-6 ridge, the same as
    in `m10`);
  - this gives a predicted P(lick) for every bin, and the same statistics are computed on the predictions.

### Statistic keys

These are the keys of `a_stats` in `m13`, used in every JSON. Each is a mean of D over the bins named.

| Key | Definition |
|---|---|
| `day14`, `day15`, `day16` | mean of D over the 30 bins of that day |
| `d14_h1` / `d14_h2` | day 14, bins 0-14 (0-7.5 min) / bins 15-29 (7.5-15 min) |
| `d15_h1` | day 15, bins 0-14 |
| `within14_halves` | `d14_h2 - d14_h1` (change within day 14) |
| `night_halves` | `d15_h1 - d14_h2` (change over the night) |
| `within14_5min` / `night_5min` | the same, using the last and first 5 min (bins 20-29, 0-9) |
| `d15_minus_d14` | `day15 - day14` |

### Sigmoid timing (Fig. 1d analogue)

- The fitted curve is `D(x) ≈ A / (1 + exp(-(x - x0)/s))`.
- It is fitted by exact least squares on a grid:
  - x0 runs from -10 to (window + 10) in steps of 0.5;
  - s runs from 0.3 to 40 in 40 log-spaced steps;
  - A is solved in closed form and clipped to |A| ≤ 1.
- **Midpoint** = x0, clipped to [0, window].
- **Start** = max(0, x0 - s·ln 9). This is the trial at which the fitted curve has made 10% of its change; the
  paper's version is "decays more than 0.1 normalized licks".
- **Unlike the paper's,** our curves are the raw change in P(lick). They are not normalised or smoothed.
- The curve is fitted per mouse (the paper's way) and to the group-mean curve.

### Tests, in order of trust

1. **Mouse-level paired t (PRIMARY).** For each mouse, take the observed statistic minus the model's predicted
   statistic. Test the n = 6 (or 7) differences against 0. No noise model is needed.
2. **Noise-matched simulation (`m15`).** See 5.4. It is reported as a range over two calibrations, "A (B)".
3. **Bootstrap over mice, leave-one-mouse-out, equal windows.** These check robustness.
4. **Naive check (`m13`, v1): withdrawn.** Every bin is drawn independently, Bernoulli(P). Its p values are far too
   small; see 5.4. **Never cite them as evidence.**

## 5.3 Results (A): day 14 (Active, n = 6)

### When the drops came

- Active mice earned **28.3 water drops** on day 14 (range 17 to 33).
- **66% of them came by 2.5 min and 88% by 7.5 min.**
- This front-loading is the PR schedule: the ratio restarts at 1 every session. By 2.5 min, 58% of the drops had come
  on day 13 (morphine), 77% on day 15 and 48% on day 17.
- Every model updates at each drop, so **every model expects most of the day-14 change early**.

### Ceiling

P(lick) in the first 2 min (bins 0-3) is **1.0 on days 13 and 14**, so a fall cannot show there. On days 16 and 17 it
is 0.88.

### Observed changes

| Statistic | Mean [95% CI] | p | Mice down / up |
|---|---|---|---|
| `day14` | **-0.056** [-0.27, +0.15] | 0.53 | 4 / 2 |
| `d14_h1` (0-7.5 min) | **+0.022** [-0.18, +0.22] | 0.79 | 1 / 3 |
| `d14_h2` (7.5-15 min) | -0.133 [-0.51, +0.25] | 0.41 | 4 / 2 |
| `d15_h1` | -0.156 [-0.40, +0.09] | 0.16 | 5 / 1 |
| `day15` | -0.267 [-0.46, -0.07] | 0.018 | 6 / 0 |
| `day16` | -0.350 [-0.60, -0.10] | 0.015 | 5 / 1 |
| `within14_halves` | **-0.156** [-0.60, +0.29] | 0.41 | 3 / 3 |
| `night_halves` | **-0.022** [-0.53, +0.49] | 0.92 | 3 / 3 |
| `within14_5min` / `night_5min` | -0.083 / -0.083 | 0.64 / 0.69 | 3 / 3 each |

### Within day 14 or overnight

**This cannot be separated.** The two changes are near mirror images across mice, and 5-min windows split the change
exactly in half. v1's "no extra change overnight" is withdrawn.

### Models versus mice, whole day 14 (slide 7 panel c)

"Mice − model" is the observed day-14 change minus the predicted change, averaged over mice, with its 95% t CI.
Holm correction is over the 4 main models.

| Model | Predicted | Mice − model [95% CI] | Mouse-level p | Holm | Noise-matched p, A (B) | Naive p (do not use) |
|---|---|---|---|---|---|---|
| instant switch | -0.224 | +0.169 [-0.02, +0.36] | **0.073** | 0.29 | 0.050 (0.021) | 0.0005 |
| RPE | -0.190 | +0.134 [-0.07, +0.34] | 0.15 | 0.46 | 0.14 (0.078) | 0.002 |
| RPE 2α | -0.147 | +0.092 [-0.11, +0.29] | 0.30 | 0.59 | 0.31 (0.23) | 0.053 |
| meta-RPE | -0.107 | +0.052 [-0.17, +0.27] | 0.58 | 0.59 | 0.59 (0.52) | 0.29 |
| RPE t-1 | -0.211 | +0.156 | 0.10 | - | 0.089 (0.041) | 0.001 |
| PH | -0.150 | +0.094 | 0.28 | - | 0.34 (0.23) | 0.054 |
| eMack | -0.129 | +0.074 | 0.39 | - | 0.47 (0.33) | 0.12 |
| no learning | 0.000 | -0.056 | 0.53 | - | 0.56 (0.48) | 0.25 |

- **Every model except "no learning" predicted a larger fall than the mice showed.** None is robustly rejected; only the instant switch
  comes near.
- After Holm correction, no mouse-level p is below 0.29.
- "No learning" passes day 14 because it predicts no change. It fails days 15 and 16 (mouse-level p = 0.018 and
  0.015).
- The learning models (RPE, RPE 2α, meta-RPE) match days 15 and 16 closely (mouse-level p 0.56 to 0.99). The
  instant switch gives 0.41 and 0.15.
- The naive p values above are from `m15`. The v1 slides printed the `m13` naive run (instant switch < 0.001,
  RPE 0.005, RPE 2α 0.050, meta-RPE 0.28). It is the same invalid check with a different random seed.

### Supplementary, not on the slides: the first half of day 14

This is the 0-7.5 min window of slide 7's KEY box. The instant switch predicts -0.198 against the observed +0.022:

- mouse-level p = 0.010; noise-matched p 0.040 (A) and 0.019 (B);
- RPE is borderline: 0.052, and 0.12 (0.083);
- RPE 2α (0.12) and meta-RPE (0.36) pass.

These p values are **uncorrected**, and the half-session window is one of several examined. Treat this as
descriptive support for "only instant switch is near rejection", not as a rejection.

### Passive (n = 7), with a caveat added while writing this guide

- **Observed:** licking rose at once on day 14 (+0.405; all 7 mice up; first half +0.314). The models predict +0.20
  (PH) to +0.32 (instant switch); mouse-level p = 0.004 to 0.008.
- **Do not interpret this as a model failure:**
  - Passive licking fell from day 12 to day 13 (0.47 → 0.28, a change of -0.19) with no change in outcome. Every model
    predicts 0 for that pair.
  - This **placebo pair fails the mouse-level test** (p = 0.046 to 0.053; section 5.4).
  - Measured from day 12, the day-14 rise is +0.22, **smaller** than the +0.30 to +0.32 the main models predict.
- The noise-matched checks do not reject any main model for Passive day 14 (A 0.37-0.44, B 0.23-0.31), and the
  first half of day 14 matches the models (observed +0.314 vs predicted +0.27 to +0.30; mouse-level p 0.62-0.90).
- So the Passive result depends on an unusually low reference day. Slide 7's notes ("more than every model
  predicted ... likely taste") do not carry this caveat.

## 5.4 Why the naive check is wrong, and the noise-matched check (`m15`)

**The problem.** The naive check treats every 30-s bin as an independent coin flip. On the stable days 11-13, real
licking is far more variable from day to day than that, and it is autocorrelated. The table compares the day-to-day
differences of residuals at the same bin (all three equal 1, 1, 0 for independent bins).

| Feature (days 11-13; range over the 4 main models) | Active | Passive |
|---|---|---|
| T1: variance of whole-session mean differences / independent-bin variance | 1.99-2.12 | 3.35-3.36 |
| T2: the same for half-session means | 2.83-3.00 | 3.23-3.24 |
| T3: lag-1 autocorrelation of the bin differences | 0.37 | 0.42 |

The naive check **rejected a no-change placebo**: Passive day 12 → 13, p = 0.0005 for every model.

**Noise model.** Licks are simulated as
`lick ~ Bernoulli( logistic( k·logit(P) + σu·z + a_t ) )`, where:

- z ~ N(0, 1) is one offset per mouse-day;
- `a_t` is AR(1) noise within the session, with stationary sd σa and coefficient φ;
- `k = sqrt(1 + π(σu² + σa²)/8)` keeps the mean close to the model's P (the logistic-normal approximation).

**Calibration.** The noise model is fitted per model and group on **days 11-13** (stable morphine days, before any
switch), to T1-T3 of the day-to-day residual differences at the same bin (11→12, 12→13).

- **Why differences:** every slide statistic compares the same bin across days. A within-session misfit that repeats
  every day cancels in these differences.
- **Grid:** σu from 0 to 2.0 in steps of 0.1; σa from 0 to 3.0 in steps of 0.1; φ ∈ {0, 0.3, 0.5, 0.7, 0.9}.
- 200 simulated sets per grid point for the calibration, then 2,000 sets for the tests. Seed 20261119.

No single setting matches T1, T2 and T3 together, so two calibrations are reported as a range. Reading "A (B)" as a
bracket of the true p is an assumption, not a proven property. For Passive, calibration A is also at the grid edge
(σa = 3.0) for 6 of 8 models.

| Calibration | Fits | Behaviour | Effect on p values |
|---|---|---|---|
| **A** | T1 + T2 + T3 (loss = squared log-ratios + ((T3 - T3obs)/0.05)²) | matches T2 and T3 but **overshoots T1** (3.5-3.9 simulated vs 1.99-2.12 observed, Active, 4 main models); σa = 3.0 is at the grid edge for all Active models | whole-session p too **large** (conservative) |
| **B** | T1 + T2 only | **undershoots T3** (0.24-0.28 simulated vs 0.37, Active, 4 main models); for Passive σa = 3.0 is at the grid edge and T2/T3 are under-matched | timing p may be too **small** |

### The placebo: day 13 minus day 12 (no outcome change)

A valid test must not reject it.

| | Active (observed -0.006) | Passive (observed -0.186) |
|---|---|---|
| naive | 0.76-0.96 (passes) | **0.0005 (fails)** |
| noise-matched A | 0.86-0.97 (passes) | 0.08-0.11 (passes) |
| noise-matched B | 0.81-0.97 (passes) | **0.02-0.03 (fails)** |
| mouse-level t | 0.72-0.90 (passes) | **0.046-0.053 (fails or borderline)** |

- **Active:** every robust test is validated by the placebo.
- **Passive:** no test here can be trusted at α = 0.05.

## 5.5 Results (B): timing of the change (Fig. 1d,e; slide 8)

### Active group-mean curve

| | Morphine → water | Water → morphine |
|---|---|---|
| Midpoint | **trial 36.5** (bootstrap 95% CI 13.5 to 68.0) | **trial 4.5** (3.5 to 15.0) |
| Fitted amplitude A, slope s | A = -0.38, s = 10.1 | A = +0.35, s = 0.3 (a step) |

- Trial 30 is the start of day 15.
- The morphine → water midpoint is later in **99.65% of 2,000 bootstraps** over mice. This is a share, not a p value.
- **Equal 60-trial windows** (days 14-15 only): 49.0 vs 4.5. The morphine → water fit is flat for midpoints of about
  28-50.
- **Leave one mouse out:** 31.5 to 44.0 vs 4.0 to 4.5. The direction is robust.

### Models, noise-free group curves

| | instant switch | RPE | RPE 2α | meta-RPE |
|---|---|---|---|---|
| Morphine → water midpoint | 1.0 | 10.0 | 30.0 | 31.5 |
| Water → morphine midpoint | 1.0 | 16.0 | 3.0 | 4.5 |
| Start (morphine → water) | 0.0 | 0.0 | 0.0 | 0.0 |

RPE 2α and meta-RPE predict a late midpoint after morphine → water and an early one after water → morphine.

### Noise-matched check of the group-curve midpoint

p shown as A (B): instant switch 0.069 (0.047), RPE 0.30 (0.22), RPE 2α 0.87 (0.84), meta-RPE 0.93 (0.95).

### Per mouse, the paper's test (Active)

| | Midpoints |
|---|---|
| Morphine → water | 42.5, 33.0, 15.0, 64.5, 73.5, 3.0 (mean 38.6) |
| Water → morphine | 1.5, 47.0, 4.5, 4.0, 11.5, 6.0 (mean 12.4) |

- Paired t p = **0.11**, Wilcoxon 0.22; later in 4 of 6 mice. With equal 60-trial windows p = 0.12.
- **Not significant.**
- All 6 morphine → water per-mouse fits and 5 of 6 water → morphine fits are steps (s at its 0.3-trial floor), so
  their start is the midpoint minus 0.66. Exception: 0910_black water → morphine (s = 1.53, start 0).

### Caveats on timing

- The water → morphine midpoint is an **upper bound**: P(lick) is near ceiling in the first 2 min.
- Every model predicts a fast water → morphine change.
- The two switches are 3 days apart, with different histories.

### Passive

- Per mouse: 12.9 vs 20.1 (p = 0.58).
- Morphine → water happens at once (group-curve midpoint 1.0).
- The water → morphine group curve (midpoint 49.0) hits its amplitude bound (A = -1.0). At face value the change is
  gradual (slide 8 notes), but with the fit at its bound and the Passive placebo failing, it is not interpretable.

## 5.6 Results (C): start of the change (Fig. 1m; slide 8 panel c)

- The Active group-curve start is **trial 14.4**: bootstrap 95% CI 0 to 46.8; leave one mouse out 6.4 to 41.8.
- Compared with noise-matched simulations of the same group curve, A (B):

| Model | Simulated median start | 95% range (A) | p, A (B) |
|---|---|---|---|
| instant switch | 0.3 | 0 to 36.9 | 0.094 (0.068) |
| RPE | 4.8 | 0 to 40.4 | 0.60 (0.52) |
| RPE 2α | 6.2 | 0 to 52.9 | 0.92 (0.81) |
| meta-RPE | 12.8 | 0 to 51.8 | 0.94 (0.99) |

- **The start excludes no model**; the smallest p is 0.068. v1's "the start rules out instant switch" is withdrawn.
- v1 also compared the mice's per-mouse starts (mean 37.9) with noise-free model curves, which all start by
  trial 1.8 (Holm p = 0.072). That comparison is biased: noise in yes/no licking makes observed curves start later
  than noise-free predictions. Do not use it.

## 5.7 What this adds, in one paragraph

Descriptively, the Active mice show the paper's pattern. Licking did not fall while most of the water drops, and
therefore most of the negative errors, arrived. It fell later, and the morphine → water change is later than the
water → morphine change on the group curve. Statistically, with 6 mice and one switch per direction:

- no model is rejected (only the instant switch comes near);
- the per-mouse timing test is not significant;
- within-day versus overnight change cannot be told apart.

A delayed oral-drug effect, a withdrawal state, or taste could produce the same lag.

## 5.8 v1 → v2: what changed and why

| v1 statement or method (18 Nov, superseded) | v2 (current) |
|---|---|
| Naive check: bins independent, "not refitted, so conservative" | **Wrong**: anti-conservative; rejects a placebo. Replaced by mouse-level paired t (primary) and the noise-matched check (A/B range) |
| "Instant switch rejected" (p < 0.001), RPE p = 0.005 | Mouse-level p 0.073 and 0.15; noise-matched 0.021-0.050 and 0.078-0.14 |
| "No extra change overnight" | Cannot be separated from change within day 14 (95% CIs -0.60 to +0.29 and -0.53 to +0.49) |
| "The start rules out instant switch" (p = 0.011) | Start excludes no model (p ≥ 0.068) |
| Midpoints compared over unequal windows (90 vs 60 trials) only | Equal 60-trial windows added (49.0 vs 4.5; per mouse p = 0.12) |
| Slide 8a notes said the curves were normalised | Raw change in P(lick); fixed |
| "Drank" water drops | **Earned**; drinking was not measured |
| `r2` in `m13` called R² | It is uncentred (1 - SSE/Σy²) |
