# 12. Analysis plan for the next cohort (written 2026-11-26, before its data)

**Status: fixed on 2026-11-26**, before any data of the next cohort (expected at the end of November 2026) were seen.
The GitHub commit that added this file records the date. A change made after the data are opened goes in section 10
with its reason; sections 1-9 are not edited.

**Scope:** the model-free licking tests that decide whether the Dec-2025 pattern replicates. Model fits stay
exploratory (section 8). Code: `project/m18_planned_tests.py`, written and tested on the Dec-2025 cohort before the
new data arrived (section 9). The same file is on the lab PC as
`mrpe_hjort2026_20261117\ANALYSIS_PLAN_next_cohort_20261126.md`, where the paths below have no `project/` prefix.

## 1. Question

In the Dec-2025 cohort, Active mice did not change their licking on the first water day and lowered it on the next
two days (README items 4 and 13; docs/04). **Primary question:** does a new cohort show the same lag?

## 2. Design (to be confirmed when the data arrive)

- The same 18-day design as Dec-2025, with 15-min PR sessions: Baseline days 1-5 (water); Yoked days 6-10 and Equal
  Access days 11-13 (morphine); Withdrawal days 14-16 (water); Re-exposure days 17-18 (morphine).
- Days are named relative to the switches: **M0** = last Equal Access day (day 13); **W1-W3** = Withdrawal days 1-3
  (days 14-16); **R1** = first Re-exposure day (day 17); **M1** = the Equal Access day before M0 (day 12). If the new
  schedule differs, only these day numbers change. They are set from the schedule before any licking data are read,
  and the change is listed in section 10.
- Groups and yoke pairs come from the raw controller logs. A Passive mouse is one whose Yoked-day sessions are replays
  (`MicroLogReplayExperiment`; `replay_source` = the partner's log), as in
  `yoked_stats_20261116/s02_yoke_from_raw_logs.py`.
- Pipeline: copies of `s02`, `s03` and `m01` in which only the raw-data root and the data folders change, then `m18`.
  No other script changes for the tests in this plan. The Dec-2025 files stay as they are; new results go into new
  files next to them.

## 3. Outcome measure

- p(lick) of a session = the share of its first 30 bins of 30 s (0-15 min) that hold at least one lick (as in m04 and
  m13). A change = p(lick) on a day minus p(lick) on its reference day.
- A session counts only if its controller log ends at 14 min (840 s) or later; otherwise that mouse-day is missing.
  In Dec-2025 every session on days 13-17 lasted at least 889.3 s, so this rule removes nothing
  there.

## 4. Inclusion rules

- A mouse enters a test only if every day that the test uses is a valid session (section 3). A mouse that dies stays
  in the tests whose days it completed.
- There are no other exclusions: no outlier removal and no exclusion based on licking. Every excluded mouse is listed
  with its reason.
- In Dec-2025 these rules give the mice used so far: 6 Active and 7 Passive. 6099 white died on day
  13 and has no day-14 session. For the S2 sensitivity test, 6100 red (Passive) drops out because its day-12 session
  ended at 795.5 s.

## 5. Tests

All tests are two-sided at alpha = 0.05 and use one value per mouse.

| ID | Role | Value per mouse | Mice | Test | Dec-2025 (reference, not evidence) |
|---|---|---|---|---|---|
| P1 | primary | lag = change on W1 minus the mean change on W2-W3 (both vs M0) = p(W1) - mean p(W2, W3) | Active, new cohort | one-sample t vs 0 | +0.25 (95% CI +0.08 to +0.43), p = 0.014; 6 of 6 mice > 0 |
| S1 | secondary | lag, as P1 | Active, both cohorts | OLS: lag ~ 1 + cohort (sum-coded); test of the intercept | (needs the new cohort) |
| S2 | secondary | change on W1 (vs M0) | Active vs Passive, new cohort | Welch t | -0.06 vs +0.40, p = 0.0011 |
| S3 | secondary | change on R1 (vs W3) | Active, new cohort | one-sample t vs 0 | +0.35, p = 0.014 |
| S4 | secondary | size of the R1 change minus size of the W1 change | Active, new cohort | one-sample t vs 0 (= paired t) | 0.35 vs 0.18, p = 0.21 |

- P1 is primary because it is the quantity behind the sample-size calculation (README item 13) and because p(M0)
  cancels: an unusual reference day (such as the low Passive day 13 of Dec-2025, README item 11) cannot move it.
- S1-S4 are corrected together with Holm's method.
- **S2 sensitivity** (reported, no claim): S2 with the mean of M1 and M0 as the reference. Dec-2025:
  -0.06 vs +0.35, p = 0.003 (n = 6 / 6).
- Pooled versions of S2-S4 (cohort as a sum-coded factor; HC3 standard errors for S2) are reported as estimates with
  95% CIs, without claims.
- The Dec-2025 values were computed while writing this plan. The P1 contrast was chosen after seeing these data, so
  its Dec-2025 p value is not evidence. The new cohort is the test.

## 6. Decision rules

- **P1 replicates** if p < 0.05 and the mean lag is above 0. Otherwise the estimate and its 95% CI are reported as
  "not replicated at this n", never as "no lag".
- A secondary result counts only after the Holm correction.
- The analysis runs once, after every mouse of the new cohort has finished the design (day 18) or has stopped. There
  are no interim looks, and no mice are added after looking. Mice that arrive later form a further cohort, analysed
  with this same plan.
- The wording rules of docs/07 still apply: drops are earned (not drunk); the phase names; the pharmacology caveat; no
  independent-bin p values.

## 7. Sample size

The planning effect is d = 1.18 (README item 13, from Dec-2025). For the within-mouse P1 test,
dz = d / sqrt(2(1 - r)), where r is the across-mouse correlation between the W1 change and the mean W2-W3 change.

| Active mice in the new cohort | 4 | 5 | 6 | 7 | 8 | 10 | 12 |
|---|---|---|---|---|---|---|---|
| power, r = 0.3 | 0.29 | 0.40 | 0.51 | 0.60 | 0.68 | 0.80 | 0.88 |
| power, r = 0.5 | 0.37 | 0.52 | 0.64 | 0.74 | 0.82 | 0.91 | 0.96 |
| power, r = 0.7 | 0.55 | 0.72 | 0.84 | 0.92 | 0.96 | 0.99 | 1.00 |

- 80% power needs 10.0, 7.7 and 5.6 Active mice for
  r = 0.3, 0.5 and 0.7.
- In Dec-2025 the observed r was 0.61 and the observed dz 1.52 (80% power at
  5.6 mice). These come from 6 mice, and from the data that suggested the effect, so they are
  probably too optimistic. Plan with r = 0.5: 8 Active mice.
- At the Dec-2025 effects, S2 (Hedges' g = -2.54) has 98% power with 6 mice per
  group, and S3 (dz = 1.51) has 84% with 6 Active mice. S4 (dz =
  0.59) has 22% with 6 mice and needs about 25
  Active mice for 80%, so it will probably stay inconclusive.

## 8. Exploratory analyses (no claims without these checks)

- The model fits, cross-validation, learning-rate bootstrap and trial-by-trial analyses (m04-m16) are rerun on the
  new cohort and on both cohorts together. They are reported as exploratory.
- Before any model claim, the group-level recovery check (m11) is rerun at the new n. Error averaging (meta-RPE) is
  claimed only if data made by meta-RPE are recovered as meta-RPE in at least 50% of simulated sets (Dec-2025:
  16% at group level, 2% per mouse).
- No p values from the independent-bin simulation, and no per-mouse BIC model choice.
- Sex and yoke pair are reported descriptively only; there are too few mice to test them.

## 9. Code

- `project/m18_planned_tests.py` computes P1, S1-S4, the S2 sensitivity test, the Holm correction and the pooled
  estimates. On the Dec-2025 cohort it reproduces the m11 changes, 95% CIs, one-sample and Welch p values and the
  asymmetry test exactly (asserted). `--selftest` checks the pooled code path.
- Dec-2025 run: `project/data/m18_planned_tests_dec2025.json` holds the reference values of sections 5 and 7.
- New cohort: `python m18_planned_tests.py <new data folder> <cohort name>` writes `m18_planned_tests.json` in that
  folder. Its field `primary.replicated` is the P1 decision.

## 10. Changes after the data were opened

None yet.
