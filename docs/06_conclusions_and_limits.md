# 6. Conclusions, limits and next designs

## 6.1 Claims table

Each claim gets one of four labels:

- **Supported**: robust to the checks in docs/04-05.
- **Descriptive**: the pattern is there, but there is no significant test behind it.
- **Not supported**: tested, and not significant or not robust.
- **Untestable**: this design cannot answer it.

| # | Claim | Label | Numbers (Active n = 6 unless noted) | Where |
|---|---|---|---|---|
| 1 | Our port reproduces the paper's model ranking on the authors' data | Supported | meta-RPE BIC change from best 2.6 vs 5.7-18.6; largest Holm p 0.010 | docs/01, slide 1 |
| 2 | Active licking does not change on the first water day; Passive licking does | Supported | Active -0.06 [-0.27, +0.15], p 0.53; Passive +0.40, p 0.0002; Welch p 0.001 | docs/04 §4.1, slide 3 |
| 3 | Active licking falls one to two days after morphine is removed and rises on the first morphine day | Supported | days 15/16: -0.27 / -0.35 (p 0.018 / 0.015); day 17: +0.35 (p 0.014) | slide 3 |
| 4 | Active: the first-day change is smaller after morphine → water than after water → morphine | Not supported, as a direct test | 0.18 vs 0.35, paired p 0.21 | slide 3 footer |
| 5 | Active: two learning rates (slow after loss, fast after return) describe licking best | Supported at group level vs one rate and instant switch; RPE 2α vs meta-RPE is a tie on held-out mice (meta-RPE better in 4/6) | ΔBIC: meta-RPE +9, RPE +11, instant switch +22; held-out +1.3 / +8.3 / +15.5; α- < α+ in 99.5% of bootstraps | docs/04, slides 4-5 |
| 6 | RPE 2α's win is not a BIC artefact | Supported | never picked for RPE-made or instant-switch-made data (0%) | §4.5 |
| 7 | Meta-RPE (error averaging) explains our data | Untestable | recovered in 2% (per mouse) / 16% (group) of its own simulations | §4.2, §4.5, slide 5 |
| 8 | Passive mice switch without a learning lag | Not supported as stated | instant switch best, but learning-made data are also called "instant switch" 52-62% of the time | §4.5 |
| 9 | Day 14: no clear fall in the first 7.5 min although 88% of the drops had been earned | Supported (observation) | +0.02 [-0.18, +0.22]; first 2 min at ceiling | §5.3, slide 7 |
| 10 | Every model predicted a larger day-14 fall than observed | Descriptive | mice − model +0.05 to +0.17; no model robustly rejected | §5.3 |
| 11 | The instant switch is rejected | Not supported | mouse-level p 0.073 (Holm 0.29); noise-matched 0.021-0.050; first half alone p 0.010, uncorrected | §5.3 |
| 12 | The fall happened within day 14 rather than overnight, or the reverse | Untestable here | within -0.16 [-0.60, +0.29]; overnight -0.02 [-0.53, +0.49] | §5.3 |
| 13 | The morphine → water change is later than the water → morphine change | Descriptive (group curve only) | midpoint 36.5 vs 4.5, later in 99.65% of bootstraps; per mouse p 0.11 | §5.5, slide 8 |
| 14 | The start of the change separates the models | Not supported | start 14.4 (CI 0-47); p ≥ 0.068 | §5.6 |
| 15 | Passive mice rose more than every model predicted on day 14 | Not supported | depends on the unusually low reference day 13; the placebo fails (p 0.046) | §5.3 |
| 16 | CD cells exist in ORBm or BMAp and ramp in Active mice | Untestable (no neural data yet); a prediction only | Active CD lowest after 54 of 68 water drops; Passive after 16 | §4.7, slide 6 |

## 6.2 Other explanations for the same behaviour

These cannot be separated from learning with behaviour alone.

1. **Delayed pharmacology.** Oral morphine is absorbed over minutes. The drug's effect, not the drop, may be the real
   "outcome". It is delayed and smeared over bins, and that alone creates a lag. The fast water → morphine change
   happened within about 2.5 min, which is before an oral dose can act. Taste may act as a learned drug cue; this was
   not tested.
2. **Withdrawal state** on days 14-16. It can slow stopping, or change licking for reasons unrelated to value.
3. **Taste.** Morphine is bitter. Passive mice (β < 0) may simply prefer water. This would explain their
   instant changes in both directions.
4. **Different histories.** Before day 14 there were 8 morphine days; before day 17 there were 3 water days. The
   paper's two cues reverse in the same session.
5. **In-session drug action** on day 17 can produce a fast return of licking.

## 6.3 What a real test needs

- **Several switches per mouse.** For example, repeated morphine ↔ water blocks. This allows per-mouse model
  comparison and recovery.
- **Probabilistic outcomes.** For example, 85% → 15% morphine drops. This is the paper's actual manipulation.
- **A cue that signals morphine versus water.** This makes trials cue-locked, as in the paper.
- **Neural recordings** (ORBm and BMAp calcium imaging), analysed with a GLM. The regressors would be δ, CD, CE and Q
  at each drop (from `m04`/`m05`), tested against a circular-shift null as in the paper's Fig. 2.
  - A slowly ramping CD-like signal in Active mice only would support meta-RPE.
  - A step would not.

## 6.4 One-sentence summary

The Active mice show the paper's behavioural signature descriptively: they do not change while errors pile up, then
change, and drug loss is slower than drug return. With 6 mice and one switch each way, however, no learning model is
rejected and error averaging cannot be tested. A delayed oral-drug effect, a withdrawal state or taste could produce
the same pattern.
