# 1. The reference paper and its six learning models

**Paper.** Hjort, ..., Stuber (2026). *Prefrontal to ventral tegmental area dynamics drive contingency degradation.*
Nature 655, 174 (2 July 2026; online 6 May 2026). doi:10.1038/s41586-026-10443-5.

**Code and data.** github.com/stuberlab/Hjort-et-al.-2026-PFC-and-reversal-learning, commit e06ae22, MIT licence.
The models are in the notebook `BehaviorModels.ipynb`. The authors' 40-mouse sample data are copied in
`project/data/hjort2026_sample_choice/` (A01 to A40, with their licence).

Statements in this file were checked against the paper text unless marked otherwise.

---

## 1.1 The paper's task

- Head-fixed mice, four odour cues. Two cues are followed by sucrose (2.5 µl of 10% sucrose) with 85% probability;
  the other two with 15% probability.
- After the mice master these contingencies, there is **one unsignalled reversal**:
  - one 85% cue becomes 15% (**H→L**, contingency degradation);
  - one 15% cue becomes 85% (**L→H**, contingency enhancement);
  - the other two cues do not change (H→H and L→L controls).
- Behaviour is **anticipatory licking** between cue and reward (0 to 2.5 s after cue onset).

### Behavioural facts from Fig. 1 that we use

| Panel | Finding | Statistics |
|---|---|---|
| Fig. 1d,e | Sigmoid fits to licking across trials; the inflection (midpoint) is later for H→L than for L→H | P = 0.0273, n = 20 mice |
| Fig. 1h | The single-trial RPE does not correlate with the next trial's behaviour change during the reversal. This motivates an *averaged* error | P = 0.3299, n = 20 |
| Fig. 1k,l | Meta-RPE fits best by negLL and by BIC | n = 40 |
| Fig. 1m | Trial at which value degradation begins. Methods: "decays more than 0.1 normalized licks or value units after reversal". Only meta-RPE does not differ from behaviour | mRPE P = 0.2086; RPE2α 0.0255; RPE 0.0002; n = 40 |

### Neural findings (context only; we do not use them)

- 13.1 ± 1.1% of mPFC neurons significantly encode the averaged error CD. This was tested with a GLM against a
  circular-shift null.
- 61.5 ± 9.45% of VTA-projecting mPFC cells encode CD.
- The authors *speculate* that CD-encoding cells may be compromised after extended volitional access to drugs of
  abuse. The paper does not test this. It is the motivation for applying the model to our opioid data.

---

## 1.2 How the authors fit the models

- The models are fitted to the **H→L cue only**: 160 trials per mouse, with the reversal at trial 60, as set in the
  notebook.
- Per trial t:
  - choice `c_t` = 1 if the mouse licked in anticipation, 0 if not;
  - reward `r_t` = 1 or 0.
- Two option values: `Q[0]` (no lick) and `Q[1]` (lick). For the H→L cue they start at `Q = [0.15, 0.85]`.
- **Choice rule (softmax):** `P(lick) = 1 / (1 + exp(-β (Q[1] - Q[0])))`.
- **Only the chosen option is updated**, in every model except Pearce–Hall and eMack.
- **Fitting:** `scipy.optimize.minimize` (L-BFGS-B), started from the notebook's grid of starting values, within
  the notebook's bounds. **BIC = k ln T + 2 negLL**, where T is the number of trials.

### The six models

In every model the prediction error is `δ_t = r_t - Q[c_t]`.

| Model (code name) | Update of the chosen value | Free parameters (bounds) |
|---|---|---|
| RPE (`RPE`), Rescorla–Wagner | `Q[c] += α δ` | α (0.01-1), β (0-50) |
| RPE with two rates (`RPE2a`, "RPE 2α") | `Q[c] += α- δ` if δ ≤ 0, else `α+ δ` | α- (0.01-1), α+ (0.01-1), β (0-10) |
| RPE gated by the previous error (`RPEt-1`) | `Q[c] += α δ |δ_{t-1}|` | α (0.01-0.5), β (0-50) |
| **Meta-RPE** (`mRPE`) | `CD ← (1-d) CD + d δ`; `CE ← (1-e) CE + e δ`; `M = max(0, -CD) + max(0, CE)`; `Q[c] += α δ M` | α (0.01-10), β (1-10), d (0.01-0.1), e (0.01-0.1) |
| Pearce–Hall (`PH`), notebook form | Two values `vy`, `vn`; associability `a` (starts 0.05); `R = r - (vy - vn)`. If R > 0: `vy += a·Be·r`. If R < 0: `vn += a·Bi·|R|`. Then `a ← γ|r - (vy - vn)| + (1-γ) a`. Clip `vy ≤ 1`, `vn ≥ 0`. Choice uses `vy - vn` | Be (0.01-1), Bi (0.05-1), γ (0.05-1), β (1-10) |
| Extended Mackintosh (`eMack`), notebook form | Like PH with `a` starting at 0.1. If R > 0: `vy += a·Be·(1 - vy + vn)·|R|`. If R < 0: `vn += a·Bi·(1 - vn + vy)·|R|`. The learning rate `a` can only **grow** (clipped to 0.05-1) | Be (0.01-0.3), Bi (0.009-0.3), θe (0.002-0.2), θi (0.002-0.19), β (1-5) |

Quirks kept on purpose from the notebook (see `project/m02_hjort_models.py`):

- In PH, the first parameter is unpacked as Be, although the notebook's start grid names it γ.
- In eMack, the start grid requires Be > Bi and θe > θi.

### Meta-RPE in words

- **CD** ("contingency degradation") and **CE** ("contingency enhancement") are two running averages of δ, with rates
  d and e.
- After an H→L reversal, δ is negative on every omitted reward. CD drifts negative over several trials, and the gain
  `M = max(0, -CD)` grows.
- Learning, and with it the drop in value, therefore **accelerates only after errors have accumulated**. This gives a
  delayed behavioural change.
- CE does the same for positive errors, which drive enhancement.
- Only a negative CD or a positive CE raises M. A negative CE does not.
- In the code, the rates are the third and fourth parameters, named `m` (the rate d of CD) and `n` (the rate e of CE).

**d and e have no fixed order.** The model does not say that CD is the slow average and CE the fast one.

- In the authors' 40 mice, the medians are d = 0.026 and e = 0.010, and d < e in only 8 of 40 mice.
- In our opioid group fits, d = 0.01 and e = 0.10, both at their bounds. This makes CD the slow average there.

"CD slow, CE fast" is a property of our fits, not of the model.

---

## 1.3 Our port and its check (scripts m02, m03, m07)

- **`m02_hjort_models.py`** ports the six models line by line. The equations, starting values, softmax, bounds,
  start grid and BIC are the same. The trial loops are compiled with numba.
- **`m03_validate_on_paper_data.py`** refits the authors' 40 mice. The paper's ranking reproduces, as the table
  below shows.

| Model | mean BIC change from best (n = 40) |
|---|---|
| meta-RPE | **2.57** |
| RPE | 5.69 |
| RPE 2α | 6.20 |
| RPE t-1 | 6.24 |
| PH | 12.35 |
| eMack | 18.63 |

- Meta-RPE is better than each other model (paired t on BIC, Holm-corrected; largest p = 0.0104).
- By negLL, meta-RPE is the best model in 30 of 40 mice.
- The paper's own Fig. 1l p values (RM ANOVA + Holm–Šidák) differ (RPE2α 0.0064, RPEt-1 0.0258, RPE 0.0339), but the
  ranking is the same.
- Median meta-RPE parameters on the authors' data: α 1.24, β 4.18, d 0.026, e 0.010.

**`m07_fig1_paper_explained.py`** redraws Fig. 1 from the authors' data. All values below are normalised licking
or model value.

| | Licking | RPE value | Meta-RPE value |
|---|---|---|---|
| Trials 0-2 after the reversal | 1.03 | 0.86 | 0.96 |
| Trials 3-5 after the reversal | 0.93 | 0.72 | 0.72 |

- So meta-RPE delays the fall in value by about 3 trials. That is less than the mice delay their licking.
- The gain M peaks 11 trials after the reversal.
- The single-trial RPE is most negative at trial 1.

The port is correct. Everything that follows is about our data.
