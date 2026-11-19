"""Step 6: model-free view of the two switches (the analogue of the paper's Fig. 1d,e,m: when does behaviour change?).

Per mouse and 3-min block (5 blocks per 15-min session): p(lick in a 30-s bin) on the switch day minus the day
before. Day 14 = morphine -> water (H->L analogue), day 17 = water -> morphine (L->H analogue).
'Delayed onset' test (per group): is the change in the first 9 min different from the change in the last 6 min?
Paired t-test and exact Wilcoxon signed-rank across mice (n = 6 Active, 7 Passive; 6099_white died on day 13).
Output: data/m06_switch_blocks.csv, data/m06_modelfree.json
"""
import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
D = HERE / "data"
EV = pickle.load(open(D / "m04_events.pkl", "rb"))
F = pd.read_csv(D / "m04_fits.csv")
GRP = F.groupby("mouse_key").group.first().to_dict()

rows = []
for m, A in EV.items():
    B = A[(A.typ == 1)].copy()
    B["block"] = (B.bin // 6).astype(int)
    for d, g in B.groupby("day"):
        for b, h in g.groupby("block"):
            rows.append(dict(mouse=m, group=GRP[m], day=int(d), block=int(b), p=float(h.c.mean())))
X = pd.DataFrame(rows)
X.to_csv(D / "m06_switch_blocks.csv", index=False)

out = {}
for sw, (d0, d1) in {"day14_morphine_to_water": (13, 14), "day17_water_to_morphine": (16, 17)}.items():
    W = X[X.day.isin([d0, d1])].pivot_table(index=["group", "mouse", "block"], columns="day", values="p").dropna()
    W["chg"] = W[d1] - W[d0]
    W = W.reset_index()
    res = {}
    for g in ("Active", "Passive"):
        G = W[W.group == g].pivot(index="mouse", columns="block", values="chg")
        early = G[[0, 1, 2]].mean(axis=1)
        late = G[[3, 4]].mean(axis=1)
        diff = late - early
        res[g] = dict(n=int(len(G)), mean_change_by_block=G.mean().round(3).tolist(),
                      sem_by_block=(G.std(ddof=1) / np.sqrt(len(G))).round(3).tolist(),
                      change_first9min=round(float(early.mean()), 3), change_last6min=round(float(late.mean()), 3),
                      late_minus_early=round(float(diff.mean()), 3),
                      p_t_late_vs_early=float(stats.ttest_rel(late, early).pvalue),
                      p_wilcoxon_exact=float(stats.wilcoxon(late, early, method="exact").pvalue),
                      n_mice_late_change_more_negative=int((diff < 0).sum()),
                      p_t_first9min_vs0=float(stats.ttest_1samp(early, 0).pvalue),
                      p_t_last6min_vs0=float(stats.ttest_1samp(late, 0).pvalue),
                      per_mouse_first9=early.round(3).to_dict(), per_mouse_last6=late.round(3).to_dict())
    out[sw] = res
# asymmetry inside each group (paper Fig. 1e: degradation slower than enhancement):
# size of the change in the first 9 min, day 14 (morphine -> water) vs day 17 (water -> morphine), paired by mouse
out["asymmetry_first9min_abs_change_day17_vs_day14"] = {}
for g in ("Active", "Passive"):
    a = pd.Series(out["day14_morphine_to_water"][g]["per_mouse_first9"]).abs()
    b = pd.Series(out["day17_water_to_morphine"][g]["per_mouse_first9"]).abs()
    k = a.index.intersection(b.index)
    out["asymmetry_first9min_abs_change_day17_vs_day14"][g] = dict(
        n=int(len(k)), mean_abs_day14=round(float(a[k].mean()), 3), mean_abs_day17=round(float(b[k].mean()), 3),
        n_mice_day17_larger=int((b[k] > a[k]).sum()), p_t=float(stats.ttest_rel(b[k], a[k]).pvalue),
        p_wilcoxon_exact=float(stats.wilcoxon(b[k], a[k], method="exact").pvalue))
print(json.dumps(out["asymmetry_first9min_abs_change_day17_vs_day14"], indent=1))
json.dump(out, open(D / "m06_modelfree.json", "w"), indent=1)
for sw, res in out.items():
    if sw.startswith("asymmetry"):
        continue
    for g, r in res.items():
        print(sw, g, "blocks", r["mean_change_by_block"], "first9", r["change_first9min"], "last6", r["change_last6min"],
              "p_t", round(r["p_t_late_vs_early"], 4), "p_w", round(r["p_wilcoxon_exact"], 4),
              "first9 vs 0 p", round(r["p_t_first9min_vs0"], 4), "last6 vs 0 p", round(r["p_t_last6min_vs0"], 4))
