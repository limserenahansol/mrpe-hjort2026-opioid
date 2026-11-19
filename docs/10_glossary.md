# 10. Glossary

| Term | Meaning here |
|---|---|
| **Active** | Mouse that earns drops by licking on the PR schedule, including the Yoked phase (days 6-10) |
| **Passive** | Yoked partner. On days 6-10 it receives a replay of its Active partner's drop times, and its licks earn nothing. From day 11 on it earns drops on PR like the Active mice |
| **Yoke pair** | One Active mouse plus one or two Passive partners (6 pairs; see docs/02) |
| **Replay** | A Passive session on days 6-10 (`MicroLogReplayExperiment`) that delivers drops at the partner's logged times, run at a different clock time |
| **PR (progressive ratio)** | Trial k needs k + 1 licks for one drop. The ratio restarts at 1 every session |
| **Drop** | One fluid delivery at the spout (`reward_cmd`). Earned on PR days, delivered on replay days. Consumption was not measured |
| **Baseline / Yoked / Equal Access / Withdrawal / Re-exposure** | Days 1-5 / 6-10 / 11-13 / 14-16 / 17-18. Equal Access was formerly called "Post" |
| **Switch** | Day 14 (morphine → water) and day 17 (water → morphine). There is one of each per mouse |
| **H→L / L→H** | The paper's reversal directions: a high-probability cue becomes low (degradation), or the reverse (enhancement). Our analogues are day 14 and day 17, for Active mice only |
| **Trial (here)** | One 30-s bin of a session; 30 bins per 15-min session |
| **Choice (here)** | 1 if the bin contains at least one logged lick, else 0 |
| **p(lick), P(lick)** | The share of bins with a lick (observed), or the model's probability of a lick in a bin |
| **Q (value)** | Expected drug from the spout. Starts at 0 on day 3, is updated at every drop, and carries over between days |
| **δ (RPE)** | Reward prediction error at a drop: `r - Q` (r = 1 for morphine, 0 for water) |
| **CD, CE** | Running averages of δ at rates d and e ("contingency degradation" / "contingency enhancement" error) |
| **M (gain)** | `max(0, -CD) + max(0, CE)`. It multiplies the learning step in meta-RPE |
| **meta-RPE (mRPE)** | The paper's model: `Q += α·δ·M`, so learning speeds up once errors accumulate |
| **RPE, RPE 2α, RPE t-1, PH, eMack** | The paper's comparison models (docs/01) |
| **Instant switch (`today's outcome`)** | Reference model: Q = the current day's outcome from the first bin. No learning lag |
| **No learning** | Reference model with β = 0 (no value term in the readout) |
| **Readout** | `logistic(b0 + β·Q + g1·τ + g2·early)`, fitted per mouse. τ = bin/29; early = first 3 min |
| **β** | Readout weight on value. It is > 0 in all Active mice and < 0 in all Passive mice |
| **Group fit** | Learning rates shared by the mice of a group; readout per mouse (profile likelihood) |
| **BIC, ΔBIC** | k ln n + 2 negLL. ΔBIC is the difference from the best model. It overstates evidence here because bins are not independent |
| **LOMO-CV** | Leave-one-mouse-out cross-validation: the learning rates are fitted on the other mice and scored on the held-out mouse |
| **Model recovery** | Simulate data from a known model, refit, and count how often the right model is picked |
| **D** | Lick (0/1) in a bin minus the same bin on the reference day (day 13 or day 16) |
| **Reference day** | Day 13 (last morphine day) for morphine → water; day 16 (last water day) for water → morphine |
| **Within day 14 / overnight** | `d14_h2 - d14_h1` / `d15_h1 - d14_h2` (halves of the 15-min session) |
| **Sigmoid, midpoint, start** | `A/(1 + exp(-(x - x0)/s))`. Midpoint = x0. Start = x0 - s·ln 9, the point where 10% of the change is reached |
| **Group curve** | The sigmoid fitted to the mean D across mice (as opposed to per mouse) |
| **Mouse-level paired t** | Per mouse, the observed statistic minus the model's predicted statistic; t-test of those n differences against 0. The primary test |
| **Naive check** | Simulation with every bin drawn independently from the model's P. **Withdrawn**: too narrow |
| **Noise-matched check** | Simulation with a mouse-day offset and AR(1) noise calibrated on days 11-13. Calibrations A and B are reported as a range "A (B)" |
| **T1, T2, T3** | Calibration features of day-to-day residual differences: whole-session dispersion, half-session dispersion, lag-1 autocorrelation |
| **Placebo** | Day 13 minus day 12: two morphine days with no outcome change. A valid test must not reject it |
| **Holm** | Holm-Bonferroni correction for multiple comparisons |
| **v1 / v2 (trial by trial)** | v1 = `m13`/`m14` (superseded); v2 = `m15`/`m16` (current) |
| **v4 deck** | The CURRENT 8-slide deck: v2 slides 1-6 + v2 slides 7-8 |
| **ORBm, BMAp** | The two brain regions for the planned calcium imaging (not part of this behavioural analysis) |
