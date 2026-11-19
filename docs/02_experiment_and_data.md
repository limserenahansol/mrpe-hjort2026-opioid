# 2. Our experiment and data

## 2.1 Cohort and groups

The Dec-2025 opioid cohort has 14 mice in 6 **yoke pairs**: 6 Active and 8 Passive. Two Actives have two Passive
partners each.

| Pair | Active | Passive partner(s) | Sex |
|---|---|---|---|
| 1 | 6100_black | 6100_orange, 6100_red | F |
| 2 | 0911_red | 0911_orange | F |
| 3 | 0911_white | 0911_black | F |
| 4 | 0910_black | 0910_orange, 0910_red | M |
| 5 | 6099_orange | 6099_red | M |
| 6 | 6099_black | 6099_white (died on day 13) | M |

- **Exclusion:** 6099_white has no data after day 13, so it has no switch day. It is left out of every comparison.
  All comparisons therefore use **n = 6 Active and n = 7 Passive**.
- The pairs come from the raw controller logs. On days 6-10, every Passive log starts with
  `"paradigm": "MicroLogReplayExperiment", "replay_source": "<partner>.jsonl"`. The pair table is
  `project/data/pairs_from_raw.csv`.
- **Yoking means replay.** A Passive session replays the partner's logged drop times; it is run at a different clock
  time, not at the same moment. The Passive mouse's own log records the drops it actually received, and the model
  uses those.

## 2.2 Task

- Free-operant licking at a spout on a **progressive ratio (PR)** schedule. Trial k needs k + 1 licks: the first
  drop comes after 1 lick, the next after 2, and so on (`req = trial + 1` holds for every trial). Each completed
  ratio gives a reward sound and **one drop** at the spout.
- The ratio restarts at 1 every session. So drops come quickly early in a session and slowly later: they are
  **front-loaded**.
- Sessions end at a fixed time of about 15 min. Over 179 PR sessions the end time averages 918 s (range 796 to
  1026 s); one session has no end tag. Every PR session ends with exactly one unfinished trial.
- **Fluid:** water on days 1-5 and 14-16; morphine solution on days 6-13 and 17-18. The concentration is not
  recorded in this project; do not state one.
- **Consumption was not measured.** A drop is a delivery. Say "earned" on PR days and "delivered" on replay days,
  not "drunk". The model assumes that every drop is consumed.

## 2.3 The 18 days: phase names

**Use these names in all new text.** They were decided on 2026-11-17. The code still writes the old labels in its
`phase` column.

| Days | Phase (current name) | Old label in the code's `phase` column | Fluid | Active | Passive |
|---|---|---|---|---|---|
| 1-2 | Baseline (training) | `Train` | water | fixed-ratio training | fixed-ratio training |
| 3-5 | Baseline | `Pre` | water | PR | PR |
| 6-10 | Yoked | `During` | morphine | PR (earns drops) | replay of the partner's drops; its licks are logged as `observed_lick` and earn nothing |
| 11-13 | **Equal Access** (formerly "Post") | `Post` | morphine | PR | PR |
| 14-16 | Withdrawal | `Withdrawal` | water | PR | PR |
| 17-18 | Re-exposure | `Re-exposure` | morphine | PR | PR |

- From day 11 on, both groups run the same PR schedule with the same fluid sequence.
- **The two switches** used by every analysis here are:
  - day 13 → day 14: morphine → water. This is the analogue of the paper's H→L, for Active mice only (see
    docs/03).
  - day 16 → day 17: water → morphine. This is the analogue of L→H.
- The histories before the two switches differ. Before day 14 there were 8 morphine days (6-13); before day 17 there
  were 3 water days (14-16). The paper's reversal happens inside one session.

## 2.4 Data provenance

- **Raw data** are controller logs, one `*.jsonl` per mouse-day, in `K:\addiction_concate_Dec_2025\dayN\<cage>\<color>\`.
  They are read only and are **not** in this package. The session list is
  `yoked_stats_20261116/data/raw_sessions_jsonl.csv`, which is not included either.
- **Event tags** used by `m01_parse_trials.py`:
  - `lick` and `observed_lick`, with `elapsed`, `trial` and `valid`;
  - `reward_cmd` (a drop);
  - `trial_start` (with `req`);
  - `experiment_end`.
- Every logged lick counts, including the 1,990 of 101,912 licks (2.0%) whose `valid` flag is False.
- **Checks.** Rewards counted from the logs equal the raw reward table in all 179 PR sessions (maximum
  difference 0). There are 4,929 PR trials, of which 4,750 were completed (`project/data/m01_check.json`).

## 2.5 Tables in `project/data/`

| File | One row per | Columns |
|---|---|---|
| `licks_all.csv` | lick or drop, in PR and replay sessions (days 3-18; training days 1-2 are not included) | `mouse_key, group, pair, day, phase` (old label), `outcome` (water/morphine), `paradigm, t` (s from session start), `kind` (lick/drop), `trial, valid` |
| `trials_pr.csv` | PR trial | `... trial, req, t_start, t_end, complete, n_licks, dur, max_gap, first_lick_lat, session_end` |
| `sessions.csv` | session (PR and replay) | `... end, n_licks, n_drops, last_lick` |
| `m13_bins.csv` | mouse × day (11-18) × 30-s bin | `mouse` (not `mouse_key`), `group, day, bin, lick` (0/1), `n_licks, drops_in_bin, p_<model>` (fitted P(lick) under each group-fitted model) |
| `m13_per_mouse.csv` | mouse | per-mouse statistics behind slides 7-8 (observed, predicted, sigmoid fits) |
| `m06_switch_blocks.csv` | mouse × day × 3-min block | `p` = share of 30-s bins with ≥ 1 lick |
| `pairs_from_raw.csv` | mouse | `mouse_key, group, pair, sex, cage` |

There are 1,162 drops on replay days (Passive, days 6-10) and 4,750 drops on PR days.

Model fitting uses days 11-18 only:

- 30 bins × 8 days = **240 bins per mouse**;
- 1,440 bins for Active (6 mice) and 1,680 for Passive (7 mice).

## 2.6 Mouse order

Per-mouse lists in the result files use this order:

- **Active:** 0910_black, 0911_red, 0911_white, 6099_black, 6099_orange, 6100_black.
- **Passive:** 0910_orange, 0910_red, 0911_black, 0911_orange, 6099_red, 6100_orange, 6100_red.

The lists are in `project/data/m13_summary.json → meta.mice`.
