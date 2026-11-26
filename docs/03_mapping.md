# 3. Mapping the paper's model onto PR licking

The **learning-rule equations and their bounds are the paper's**. What changed is the definition of a trial, a
choice and an outcome, and the readout. Source: `project/m04_fit_user_data.py` (docstring and `build_events`). Its docstring says "every drop the mouse
drank"; read it as "received" (consumption was not measured).

## 3.1 Mapping table

| | Hjort 2026 (Pavlovian) | Our task (free-operant PR licking) |
|---|---|---|
| Trial | one cue presentation (~20 s apart) | one **30-s bin** of the 15-min session (30 bins, 0-900 s) |
| Choice | anticipatory lick: yes / no | **any logged lick in the bin**: yes / no |
| Outcome r | sucrose 1 / none 0; the cue's *probability* changes 85% ↔ 15% | **each drop: morphine = 1, water = 0**; the drop's *identity* changes, and every completed ratio gives a drop |
| Value Q | expected reward of the cue | **expected drug from the spout**: one number, updated at every drop |
| Switch | one reversal inside a session | **day 14** morphine → water; **day 17** water → morphine |
| Data fitted | 160 trials (H→L cue) per mouse | **days 11-18**: 240 bins per mouse |
| Readout | softmax(β·(Q1 - Q0)), β > 0 | **logistic(b0 + β·Q + g1·τ + g2·early)**; β may take either sign |

## 3.2 Value bookkeeping, exactly

- Q starts at **0 on day 3**, the first PR day. Only water had been given before.
- Q is updated at **every drop on days 3-18**. This includes Active PR drops on days 6-10 and the Passive replay
  drops on days 6-10: both groups received morphine then.
- Q is **carried over between days unchanged**. There is no overnight decay.
- **Event order.** Bins and drops are interleaved in time. Bin b covers `[30b, 30b + 30)` s. Its choice is predicted
  from Q **after all drops that happened before the bin started**; a drop inside bin b updates Q only after bin b has
  been scored. Drops after 900 s still update Q before the next day.
- **Likelihood:**
  - only bins of PR sessions on **days 11-18** are fitted (`use = 1`);
  - PR bins on days 3-10 (Active days 3-10, Passive days 3-5) exist but are not fitted;
  - Passive replay sessions (days 6-10) contribute drops only, with no bins.

- **Sessions shorter than 900 s.** Every session is scored over 30 bins (0-900 s), so bins after an early session
  end count as "no lick". 21 fitted sessions end before 900 s. Example: 6100_red (Passive) day 12 ends at 795.5 s,
  so bins 26-29 are zeros by construction; day 12 is one side of the Passive placebo pair. This is a known
  simplification, not corrected in the results.

## 3.3 Readout, the same for every model

`P(lick in bin) = 1 / (1 + exp(-(b0 + β·q + g1·τ + g2·early)))`

- `τ = bin / 29`, which runs from 0 to 1 across the session. This captures the fall in licking within every session.
- `early = 1` for bins 0-5 (the first 3 min).
- `q = Q` for RPE-type models. For PH and eMack, `q = vy - vn`.
- For the instant-switch model, q = the current day's outcome (1 on morphine days, 0 on water days).
- Readout bounds: b0 from -10 to 10; β, g1 and g2 from -20 to 20.

Why the readout differs from the paper's:

1. In an operant task, a drop only ever follows licking. So "update the chosen option" means updating the lick value,
   and the no-lick value stays at its start. The paper's two-option softmax then reduces to a logistic in Q with an
   intercept.
2. Licking falls within every session, whatever the outcome. The time terms absorb this.
3. A mouse may lick *less* when it expects morphine. So β is free in sign.

**One guard.** The meta-RPE step α·M is capped at 1. Our outcomes are certain, so M can approach 1, compared with
about 0.04 in the paper; uncapped, α up to 10 would overshoot. The cap was never reached in any fit
(`n_steps_capped = 0` for every mouse).

## 3.4 Models fitted to our data

Model codes in `m04` are 0 mRPE, 1 RPE2a, 2 PH, 3 eMack, 4 RPEt-1, 5 RPE, 6 no learning, 7 today's outcome.

| Name in code | Name on slides | What it is |
|---|---|---|
| `mRPE` | meta-RPE | the paper's model (α, d, e) |
| `RPE2a` | RPE 2α | separate rates for negative and positive errors (α-, α+) |
| `RPE` | RPE | one learning rate |
| `RPEt-1` | RPE t-1 | rate gated by the previous error |
| `PH`, `eMack` | Pearce–Hall, eMack | the paper's dynamic-rate models |
| `today's outcome` | **instant switch** | Q = 1 on morphine days and 0 on water days, from the first bin. No learning lag and no learning parameters |
| `no learning` | no learning | β = 0: the readout has no value term |

## 3.5 Sign of β, and why the H→L analogy holds only for Active mice

In the per-mouse meta-RPE fits (`project/data/m04_summary.json`):

- **β > 0 in all 6 Active mice** (1.25 to 3.61);
- **β < 0 in all 7 Passive mice** (-1.25 to -2.53).

The RPE fits give the same signs.

- Active mice lick more when they expect morphine. For them, morphine → water is a **reward omission**, the
  analogue of the paper's H→L degradation.
- Passive mice lick less when they expect morphine. For them, morphine → water is **relief, not an omission**.
- Passive "learning" results should therefore not be read as contingency degradation. Their fast changes may be
  taste: morphine is bitter.

## 3.6 What the mapping cannot reproduce

- **Probability.** Our outcomes change in identity (every drop switches), not in probability.
- **Cues.** There is no cue-locked trial structure. A "trial" is a time bin.
- **One switch per direction.** Each mouse has one switch each way, 3 days apart, with different histories. The
  paper has one reversal per session with value-matched control cues.
- **Delayed outcome.** Oral morphine is absorbed over minutes, so the pharmacological "outcome" of a drop is delayed
  and smeared across bins. The model treats it as immediate. This alone can create a lag that looks like slow
  learning.
