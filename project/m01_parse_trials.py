"""Step 1: one row per PR trial (= one reward and the licks that earned it) from the raw controller logs.

Raw source (read only): every *.jsonl listed in yoked_stats_20261116/data/raw_sessions_jsonl.csv
(K:/addiction_concate_Dec_2025/dayN/<cage>/<color>/). One log per mouse-day (already checked there).
PR sessions (days 3-18; Active also 6-10): trial k needs k+1 licks, a reward sound + drop at completion,
session ends at a fixed ~15 min (experiment_end). Passive days 6-10 are replays (MicroLogReplayExperiment):
drops at the partner's times, the mouse's own licks are logged as observed_lick and earn nothing.

Outcome per day (study design): water on days 1-5 and 14-16, morphine on days 6-13 and 17-18.
Output: data/trials_pr.csv (one row per PR trial), data/licks_all.csv (every lick and drop time),
        data/sessions.csv (one row per session), data/m01_check.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
D = HERE / "data"
Y = Path(r"C:\Users\hsollim\Research_Projects\01_Opioid_Behavior\yoked_stats_20261116\data")
R = pd.read_csv(Y / "raw_sessions_jsonl.csv")
P = pd.read_csv(Y / "pairs_from_raw.csv")
R = R.merge(P[["mouse_key", "group", "pair", "sex"]], on="mouse_key", how="left")


def phase_of(d):
    return ("Train" if d <= 2 else "Pre" if d <= 5 else "During" if d <= 10 else "Post" if d <= 13
            else "Withdrawal" if d <= 16 else "Re-exposure")


def outcome_of(d):
    return "morphine" if (6 <= d <= 13 or d >= 17) else "water"


trials, licks, sess = [], [], []
for r in R.itertuples():
    if r.paradigm not in ("ProgressiveRatioExperiment", "MicroLogReplayExperiment"):
        continue
    ev = []
    for ln in open(r.path, encoding="utf-8"):
        try:
            ev.append(json.loads(ln))
        except Exception:
            pass
    end = next((e["elapsed"] for e in ev if e.get("tag") == "experiment_end"), np.nan)
    lk = [(e["elapsed"], e.get("trial"), bool(e.get("valid", True))) for e in ev if e.get("tag") in ("lick", "observed_lick")]
    rw = [(e["elapsed"], e.get("trial")) for e in ev if e.get("tag") == "reward_cmd"]
    ts = {}
    for e in ev:
        if e.get("tag") == "trial_start" and e.get("trial") not in ts:
            ts[e["trial"]] = (e["elapsed"], e.get("req"))
    base = dict(mouse_key=r.mouse_key, group=r.group, pair=r.pair, day=r.day, phase=phase_of(r.day),
                outcome=outcome_of(r.day), paradigm=r.paradigm)
    for t, k, v in lk:
        licks.append(dict(base, t=t, kind="lick", trial=k, valid=v))
    for t, k in rw:
        licks.append(dict(base, t=t, kind="drop", trial=k, valid=True))
    sess.append(dict(base, end=end, n_licks=len(lk), n_drops=len(rw), last_lick=max([x[0] for x in lk], default=np.nan)))
    if r.paradigm != "ProgressiveRatioExperiment":
        continue
    lick_t = np.array(sorted(x[0] for x in lk))
    rw_by_trial = {k: t for t, k in rw}
    for k in sorted(ts):
        t0, req = ts[k]
        done = k in rw_by_trial
        t1 = rw_by_trial[k] if done else end
        inside = lick_t[(lick_t > t0) & (lick_t <= t1 + 1e-9)]
        pts = np.r_[t0, inside, t1]
        gaps = np.diff(pts)
        trials.append(dict(base, trial=k, req=req, t_start=t0, t_end=t1, complete=done, n_licks=len(inside),
                           dur=t1 - t0, max_gap=float(gaps.max()) if len(gaps) else np.nan,
                           first_lick_lat=float(inside[0] - t0) if len(inside) else np.nan, session_end=end))

T = pd.DataFrame(trials).sort_values(["mouse_key", "day", "trial"])
L = pd.DataFrame(licks).sort_values(["mouse_key", "day", "t"])
S = pd.DataFrame(sess).sort_values(["mouse_key", "day"])
T.to_csv(D / "trials_pr.csv", index=False)
L.to_csv(D / "licks_all.csv", index=False)
S.to_csv(D / "sessions.csv", index=False)
chk = dict(n_pr_sessions=int(T.groupby(["mouse_key", "day"]).ngroups), n_trials=int(len(T)),
           n_complete=int(T.complete.sum()), incomplete_per_session=T.groupby(["mouse_key", "day"]).complete.apply(lambda s: int((~s).sum())).value_counts().to_dict(),
           req_equals_trial_plus_1=bool((T.req == T.trial + 1).all()),
           session_end_s=S[S.paradigm == "ProgressiveRatioExperiment"].end.describe().round(1).to_dict(),
           rewards_vs_raw_table=None)
# cross-check against the raw table used for the breakpoint statistics
n_rw = T[T.complete].groupby(["mouse_key", "day"]).size().rename("n").reset_index()
m = n_rw.merge(R[["mouse_key", "day", "n_rewards"]], on=["mouse_key", "day"])
chk["rewards_vs_raw_table"] = dict(n=int(len(m)), max_absdiff=int((m.n - m.n_rewards).abs().max()))
json.dump(chk, open(D / "m01_check.json", "w"), indent=1, default=str)
print(json.dumps(chk, indent=1, default=str))
