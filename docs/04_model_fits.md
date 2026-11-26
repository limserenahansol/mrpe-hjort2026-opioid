# 4. Model-free tests, model fits, recovery and asymmetry

All results in this file come from the first analysis round (17 Nov 2026), after one independent review. They are
shown on slides 2-6. Unless stated, Active n = 6 and Passive n = 7.

## 4.1 Model-free switch tests (`m06`, `m11`)

**Measure.** `p(lick)` is the share of a session's 30-s bins that contain at least one lick. The **change** is a
mouse's session mean minus its reference day: day 13 for days 14-16, and day 16 for days 17-18.

**Tests.** One-sample t against 0 (two-sided); Welch t for Active versus Passive. Source:
`project/data/m11_modelfree.json`.

| Comparison | Active: mean [95% CI], p | Passive: mean [95% CI], p | Active vs Passive (Welch p) |
|---|---|---|---|
| day 14 vs 13 (first water day) | **-0.056** [-0.265, +0.154], 0.53 | **+0.405** [+0.277, +0.532], 0.0002 | 0.0011 |
| day 15 vs 13 | -0.267 [-0.465, -0.069], 0.018 | +0.367 [+0.256, +0.478], 0.0002 | < 0.0001 |
| day 16 vs 13 | -0.350 [-0.597, -0.103], 0.015 | +0.429 [+0.182, +0.675], 0.0053 | 0.0002 |
| day 17 vs 16 (first morphine day) | **+0.350** [+0.107, +0.593], 0.014 | -0.281 [-0.477, -0.085], 0.013 | 0.0004 |
| day 18 vs 16 | +0.289 [+0.102, +0.476], 0.011 | -0.467 [-0.636, -0.297], 0.0005 | < 0.0001 |

In words:

- **Active:** no detectable change on the first water day. Licking fell on days 15-16 and rose at once when morphine
  returned.
- **Passive:** licking changed on the first day in both directions, with the opposite signs (β < 0).

**Asymmetry test (Active).** Compare the *size* of the first-day change at the two switches:

- |day 14 change| = 0.18 and |day 17 change| = 0.35 (mean of each mouse's absolute change, not the absolute
  group mean);
- 4 of 6 mice changed more on day 17;
- paired t p = 0.21 (Wilcoxon 0.31).

The direction matches the paper's "H→L slower than L→H", but this test alone is **not significant**.

**Superseded.** `m06` also reports changes over the first 9 min of each switch day (Active day 14: 0.00, p = 1.00;
day 17: +0.32, p = 0.038). That window was chosen post hoc. Use the whole-session values from `m11`.

## 4.2 Per-mouse model comparison (`m04`) is not informative here

Each of the 8 models was fitted to each mouse by maximum likelihood: L-BFGS-B from a grid of starting values, with
BIC computed on n = 240 bins. The mean BIC change from the best model:

| | Active | Passive |
|---|---|---|
| instant switch (today's outcome) | 2.19 | 0.50 |
| RPE | 3.62 | 3.94 |
| RPE t-1 | 4.28 | 4.83 |
| RPE 2α | 5.78 | 8.47 |
| meta-RPE | 9.91 | 13.88 |
| no learning | 12.30 | 22.60 |
| eMack | 16.27 | 19.27 |
| PH | 22.38 | 33.12 |

**Why this ranking means little: model recovery (`m05` part a).** We simulated licking from each mouse's own fitted
model, keeping its real drop sequence; this is exact, because value depends only on the drops. There were 20 sets per
mouse for each generating model (meta-RPE, RPE, instant switch). We then refitted those three models and asked which
one BIC picks.

| Generating model | Active: picked correctly | Passive: picked correctly |
|---|---|---|
| meta-RPE | **1.7%** (47.5% → RPE, 50.8% → instant switch) | **0.7%** |
| RPE | 42.5% | 23.6% |
| instant switch | 99.2% | 97.9% |

By per-mouse negLL (no complexity penalty), meta-RPE is best in 4 of 6 Active and 3 of 7 Passive mice. Per-mouse
BIC leans towards the simplest model. So "meta-RPE lost, per mouse" means nothing here. The recovered
rates also track the true ones poorly: the correlation is 0.18 for d and 0.41 for e.

## 4.3 Group fits with shared learning rates (`m05` part b)

**Method.** The learning-rule parameters are shared by all mice of a group. The readout (b0, β, g1, g2) is fitted per
mouse, solved exactly by Newton logistic regression for every candidate learning parameter (a profile likelihood).
Only the 1-4 learning parameters are searched.

BIC uses k = (learning parameters) + 4 × (number of mice) and n = all bins of the group (1,440 Active, 1,680
Passive). **BIC counts 30-s bins as independent evidence, so these ΔBIC values overstate the evidence.**

| Active model | ΔBIC | negLL | Learning parameters |
|---|---|---|---|
| **RPE 2α** | **0** | 809.3 | α- 0.010 (lower bound), α+ 0.119 |
| meta-RPE | +9.1 | 810.3 | α 0.181, d 0.010, e 0.100 (both at the paper's bounds) |
| RPE | +11.3 | 818.6 | α 0.039 |
| eMack | +12.7 | 808.4 | 4 parameters |
| RPE t-1 | +13.5 | 819.7 | α 0.065 |
| instant switch | +22.0 | 827.6 | none |
| PH | +70.5 | 840.9 | 3 parameters |
| no learning | +71.9 | 874.4 | none |

| Passive model | ΔBIC | negLL | Learning parameters |
|---|---|---|---|
| **instant switch** | **0** | 1001.9 | none |
| RPE | +7.3 | 1001.8 | α 0.607 |
| RPE t-1 | +9.1 | 1002.7 | α 0.500 |
| RPE 2α | +13.7 | 1001.3 | α- 0.170, α+ 0.625 |
| meta-RPE | +20.2 | 1000.8 | α 3.29, d 0.010, e 0.100 |
| eMack | +28.5 | 1001.3 | 4 parameters |
| PH | +139.3 | 1060.4 | 3 parameters |
| no learning | +141.1 | 1098.4 | none |

- Meta-RPE fits Active mice slightly *worse* than RPE 2α (negLL 810.3 vs 809.3) and has one parameter more.
- In Passive mice, all learning rules reach almost the same negLL as the instant switch, which has no learning
  parameters.

**Bootstrap over mice (1,000 resamples).**

- Active: d < e in 98.3%; meta-RPE beats RPE by BIC in 53%.
- Passive: d < e in 90.8%; meta-RPE beats RPE in 0.1%.

## 4.4 Held-out mice: leave-one-mouse-out cross-validation (`m11`)

**Method.** The learning rates are fitted on the other mice. The held-out mouse's readout is then fitted on that mouse,
and we score its negLL. This tests the **learning rule**, not the readout, and it does not treat bins as independent
evidence for the learning rates.

Total held-out negLL, relative to the best model:

| Active (best: RPE 2α, 812.16) | Δ | Passive (best: instant switch, 1001.87) | Δ |
|---|---|---|---|
| meta-RPE | +1.31 | RPE t-1 | +0.84 |
| RPE | +8.29 (RPE 2α better in 5/6 mice; Wilcoxon p = 0.16) | meta-RPE | +1.45 |
| RPE t-1 | +8.86 | RPE | +1.80 |
| instant switch | +15.47 | RPE 2α | +2.61 |
| no learning | +62.24 | no learning | +96.54 |

RPE 2α versus meta-RPE on held-out mice is a tie: RPE 2α is better in only 2 of 6 mice (Wilcoxon p = 0.69). The
clear held-out result is that both beat one-rate RPE and the instant switch.

## 4.5 Group-level model recovery (`m11`)

**Method.** We simulated 50 data sets per generating model from the group fits, using each mouse's own readout, then
refitted four models at the group level (RPE 2α, RPE, meta-RPE, instant switch). The table shows the share picked by
BIC.

| Generating model | Active | Passive |
|---|---|---|
| RPE | RPE 100% | RPE 38%, instant switch 62% |
| RPE 2α | **RPE 2α 56%**, RPE 40%, instant switch 4% | RPE 42%, instant switch 58% |
| meta-RPE | RPE 2α 54%, RPE 30%, **meta-RPE 16%** | RPE 48%, instant switch 52% |
| instant switch | instant switch 100% | instant switch 100% |

What this means:

- **Active:** RPE 2α is never picked for data made by RPE or by instant switch (0%). So its win is not a BIC false
  positive.
- **Meta-RPE is not recoverable** (16%). Its data are mostly called RPE 2α.
- **Passive:** data made with learning are called "instant switch" 52-62% of the time. So the Passive instant-switch
  result does **not** exclude fast learning.

## 4.6 Learning-rate asymmetry (`m08`)

The RPE 2α group fit was bootstrapped over mice 1,000 times.

Values are the full-data fit [bootstrap 95% CI].

| | α- (negative errors) | α+ (positive errors) | Share with α- < α+ |
|---|---|---|---|
| Active | 0.010 [0.010, 0.032], at the lower bound | 0.119 [0.053, 1.00] | **99.5%** |
| Passive | 0.170 [0.010, 1.00] | 0.625 [0.122, 1.00] | 78.1% (unclear) |

- In Active mice, log10(α+/α-) has a 95% interval of 0.29 to 2.0.
- Active minus Passive α-: [-0.99, 0.00]; Active is lower in 88.5% of bootstraps.
- **Robust part:** in Active mice, updating is slow after morphine is lost and fast after it returns.
- **The bound matters.** α- sits at its lower bound (0.01). The true value could be smaller, so the bound limits the
  estimate.

## 4.7 Model variables per drop (slide 6): a prediction, not a result

From the group-fitted meta-RPE we computed, at every drop, the RPE δ, CD, CE, the gain M and the value Q. These are
the regressors a neural GLM would use, as in the paper's Fig. 2 (with a circular-shift null).

- **Example pair 4** (Active 0910_black, Passive 0910_orange):
  - Active: CD is lowest after 54 of 68 water drops on days 14-16, so a CD cell would **ramp slowly**.
  - Passive: CD is lowest 16 drops after the switch, so a CD cell would show a **step**.
- With shared parameters the counts are nearly the same in every mouse: Active 45-54, Passive 16.
- These counts were computed when slide 6 was built (`project/m12_make_slides_v2.py`). They are not stored in a JSON.
- **Caveats:**
  - meta-RPE lost the group comparison (ΔBIC +9 Active, +20 Passive);
  - d and e sit at the paper's bounds, and d is poorly recovered (r = 0.18), so the CD time course is a model
    assumption;
  - the winning RPE 2α gives δ and Q but no averaged error;
  - the paper found CD cells in mPFC (prelimbic) and in mPFC→VTA cells; ORBm and BMAp were not studied.
