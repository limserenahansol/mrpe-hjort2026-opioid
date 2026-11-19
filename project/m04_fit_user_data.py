"""Step 4: the six Hjort et al. 2026 learning rules applied to the Dec-2025 opioid cohort.

How the paper's trial maps onto our task (free-operant PR licking, 15-min sessions):
  paper                                   here
  trial (cue, ~20 s apart)            ->  30-s bin of the session (30 bins = 15 min)
  choice: anticipatory lick yes/no    ->  any lick in the bin yes/no ("seeking")
  outcome r: sucrose 1 / omission 0   ->  every drop: morphine 1 / water 0 (outcome identity)
  value Q: expected reward            ->  expected drug from the spout, updated at every drop
  reversal H->L / L->H                ->  day 14 morphine->water / day 17 water->morphine (same for both groups)
Value starts at 0 on day 3 (only water before), is updated by every drop the mouse drank (Active PR days 6-10 and
Passive replay days 6-10 included: both drank morphine), and carries over between days. The likelihood uses days
11-18 (8 x 30 bins), where both groups work on the same PR schedule with the same outcome sequence.

Choice rule. In an operant task a drop only follows licking, so 'update only the chosen option' (paper) means the
lick value is updated at each drop and the no-lick value stays at its start; the paper's softmax then becomes
p(lick) = logistic(b0 + beta*Q). Added for every model alike: time in session (linear) and a first-3-min term,
because licking falls within every session whatever the outcome. beta may be negative (a mouse that licks less
when it expects morphine). Learning-rule bounds = the paper's bounds. One guard: the mRPE step alpha*|M| is capped at 1
(our outcomes are certain, so |M| reaches ~1, far above the paper's ~0.04; uncapped, alpha up to 10 overshoots).
Extra reference models: 'no learning' (beta = 0) and 'today's outcome' (Q = 1 on morphine days, 0 on water days,
known from the first bin: an instant switch with no learning lag).
Output: data/m04_events.pkl, data/m04_fits.csv, data/m04_traces.pkl, data/m04_summary.json
"""
import itertools
import json
import pickle
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from numba import njit
from scipy import stats
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
D = HERE / "data"
BW, NB = 30.0, 30
FIT_DAYS = list(range(11, 19))


# ------------------------------------------------------------------ events per mouse
def build_events():
    L = pd.read_csv(D / "licks_all.csv")
    ev = {}
    for m, g in L.groupby("mouse_key"):
        rows = []
        for d, s in g.groupby("day"):
            if d < 3:
                continue
            r_day = 1.0 if (6 <= d <= 13 or d >= 17) else 0.0
            drops = np.sort(s[s.kind == "drop"].t.values)
            licks = np.sort(s[s.kind == "lick"].t.values)
            pr = s.paradigm.iloc[0] == "ProgressiveRatioExperiment"
            use = pr and d in FIT_DAYS
            k = 0
            for b in range(NB):
                t0 = b * BW
                while k < len(drops) and drops[k] < t0:
                    rows.append((0, d, drops[k], r_day, 0, 0.0, 0.0, 0, r_day, -1))
                    k += 1
                if pr:
                    c = int(((licks >= t0) & (licks < t0 + BW)).any())
                    rows.append((1, d, t0, np.nan, c, b / (NB - 1), float(b < 6), int(use), r_day, b))
            while k < len(drops):
                rows.append((0, d, drops[k], r_day, 0, 0.0, 0.0, 0, r_day, -1))
                k += 1
        A = pd.DataFrame(rows, columns=["typ", "day", "t", "r", "c", "tau", "early", "use", "dayout", "bin"])
        ev[m] = A
    return ev


def arrays(A):
    return (A.typ.values.astype(np.int64), np.nan_to_num(A.r.values).astype(np.float64), A.c.values.astype(np.int64),
            A.tau.values.astype(np.float64), A.early.values.astype(np.float64), A.use.values.astype(np.int64),
            A.dayout.values.astype(np.float64))


# ------------------------------------------------------------------ learning rules (one value Q, updated at drops)
# model codes: 0 mRPE, 1 RPE2a, 2 PH, 3 eMack, 4 RPEt-1, 5 RPE, 6 no learning, 7 today's outcome
@njit(cache=True)
def run(code, lp, ro, typ, r, c, tau, early, use, dayout, want_trace):
    b0, beta, g1, g2 = ro[0], ro[1], ro[2], ro[3]
    Q = 0.0
    md = 0.0
    me = 0.0
    dprev = 0.0
    vy = 0.0
    vn = 0.0
    a_ph = 0.05
    a_mk = 0.1
    capped = 0
    nll = 0.0
    n = typ.shape[0]
    tr = np.zeros((n, 6)) if want_trace else np.zeros((1, 6))
    for i in range(n):
        if typ[i] == 1:                                  # a 30-s bin: choice likelihood
            q = Q
            if code == 2 or code == 3:
                q = vy - vn
            elif code == 7:
                q = dayout[i]
            x = b0 + g1 * tau[i] + g2 * early[i]
            if code != 6:
                x += beta * q
            p = 1.0 / (1.0 + np.exp(-x))
            if want_trace:
                tr[i, 0] = q
                tr[i, 5] = p
            if use[i] == 1:
                pp = p if c[i] == 1 else 1.0 - p
                nll -= np.log(max(pp, 1e-300))
            continue
        # a drop: outcome r[i]
        if code == 0:
            delta = r[i] - Q
            md = (1 - lp[1]) * md + lp[1] * delta        # CD, rate d
            me = (1 - lp[2]) * me + lp[2] * delta        # CE, rate e
            M = (-md + abs(md)) / 2 + (me + abs(me)) / 2
            step = lp[0] * abs(M)
            if step > 1.0:
                step = 1.0
                capped += 1
            Q = Q + step * delta
            if want_trace:
                tr[i, 0] = Q
                tr[i, 1] = delta
                tr[i, 2] = md
                tr[i, 3] = me
                tr[i, 4] = M
        elif code == 1:
            delta = r[i] - Q
            Q = Q + (lp[0] if delta <= 0 else lp[1]) * delta
            if want_trace:
                tr[i, 0] = Q
                tr[i, 1] = delta
        elif code == 4:
            delta = r[i] - Q
            Q = Q + lp[0] * delta * dprev
            dprev = abs(delta)
            if want_trace:
                tr[i, 0] = Q
                tr[i, 1] = delta
        elif code == 5:
            delta = r[i] - Q
            Q = Q + lp[0] * delta
            if want_trace:
                tr[i, 0] = Q
                tr[i, 1] = delta
        elif code == 2:                                   # Pearce-Hall (notebook form)
            Be, Bi, gam = lp[0], lp[1], lp[2]
            vnet = vy - vn
            R = r[i] - (vy - vn)
            dvy = 0.0
            dvn = 0.0
            if R > 0:
                dvy = a_ph * Be * r[i]
            elif R < 0:
                dvn = a_ph * Bi * abs(R)
            a_ph = gam * abs(r[i] - vnet) + (1 - gam) * a_ph
            vy = min(vy + dvy, 1.0)
            vn = max(vn + dvn, 0.0)
            if want_trace:
                tr[i, 0] = vy - vn
                tr[i, 1] = R
        elif code == 3:                                   # extended Mackintosh (notebook form)
            Be, Bi, te, ti = lp[0], lp[1], lp[2], lp[3]
            R = r[i] - (vy - vn)
            dvy = 0.0
            dvn = 0.0
            if R > 0:
                dvy = a_mk * Be * (1 - vy + vn) * abs(R)
                da = -te * abs(r[i] - vy + vn)
            elif R < 0:
                dvn = a_mk * Bi * (1 - vn + vy) * abs(R)
                da = -ti * abs(R * abs(-vn + vy))
            else:
                da = 0.0
            vy = vy + dvy
            vn = vn + dvn
            a_mk = min(max(a_mk - da, 0.05), 1.0)
            if want_trace:
                tr[i, 0] = vy - vn
                tr[i, 1] = R
    return nll, capped, tr


RO_BOUNDS = [(-10, 10), (-20, 20), (-20, 20), (-20, 20)]
SPECS = {  # name: (code, learning bounds = paper's, learning start grid)
    "mRPE": (0, [(0.01, 10), (0.01, 0.1), (0.01, 0.1)], [[0.1, 1, 5], [0.01, 0.05, 0.1], [0.01, 0.05, 0.1]]),
    "RPE2a": (1, [(0.01, 1), (0.01, 1)], [[0.02, 0.2, 0.8], [0.02, 0.2, 0.8]]),
    "PH": (2, [(0.01, 1), (0.05, 1), (0.05, 1)], [[0.05, 0.5, 0.9], [0.05, 0.5, 0.9], [0.05, 0.5, 0.9]]),
    "eMack": (3, [(0.01, 0.3), (0.009, 0.3), (0.002, 0.2), (0.002, 0.19)],
              [[0.01, 0.13, 0.25], [0.009, 0.12, 0.23], [0.01, 0.1], [0.009, 0.09]]),
    "RPEt-1": (4, [(0.01, 0.5)], [[0.02, 0.1, 0.3, 0.5]]),
    "RPE": (5, [(0.01, 1)], [[0.01, 0.03, 0.1, 0.3, 0.9]]),
    "no learning": (6, [], []),
    "today's outcome": (7, [], []),
}
ORDER = ["mRPE", "RPE2a", "PH", "eMack", "RPEt-1", "RPE", "today's outcome", "no learning"]


def fit_one(name, arr):
    code, lb, grid = SPECS[name]
    typ, r, c, tau, early, use, dayout = arr
    nb = len(lb)
    bounds = lb + RO_BOUNDS
    cu = c[(typ == 1) & (use == 1)]
    pm = min(max(cu.mean(), 0.02), 0.98)
    b0 = np.log(pm / (1 - pm))
    starts = []
    for g in itertools.product(*grid) if nb else [()]:
        for bsign in ((-3.0, 3.0) if code != 6 else (0.0,)):
            starts.append(np.r_[np.array(g, float), [b0 - bsign / 2, bsign, 0.0, 0.0]])

    def f(x):
        return run(code, x[:nb], x[nb:], typ, r, c, tau, early, use, dayout, False)[0]

    best = (np.inf, None)
    for x0 in starts:
        res = minimize(f, x0, bounds=bounds, method="L-BFGS-B")
        if res.fun < best[0]:
            best = (float(res.fun), res.x.copy())
    k = nb + 4 - (1 if code == 6 else 0)
    n = int(cu.size)
    capped = int(run(code, best[1][:nb], best[1][nb:], typ, r, c, tau, early, use, dayout, False)[1])
    return dict(model=name, negLL=best[0], BIC=k * np.log(n) + 2 * best[0], k=k, n_bins=n,
                params=best[1].tolist(), n_steps_capped=capped)


def fit_mouse(args):
    m, A = args
    arr = arrays(A)
    return [dict(mouse_key=m, **fit_one(nm, arr)) for nm in ORDER]


def holm(ps):
    ps = np.asarray(ps, float)
    o = np.argsort(ps)
    adj = np.empty_like(ps)
    run_ = 0.0
    for i, j in enumerate(o):
        run_ = max(run_, min(1.0, (len(ps) - i) * ps[j]))
        adj[j] = run_
    return adj


if __name__ == "__main__":
    EV = build_events()
    pickle.dump(EV, open(D / "m04_events.pkl", "wb"))
    P = pd.read_csv(Path(r"C:\Users\hsollim\Research_Projects\01_Opioid_Behavior\yoked_stats_20261116\data") / "pairs_from_raw.csv")
    grp = dict(zip(P.mouse_key, P.group))
    with ProcessPoolExecutor(max_workers=14) as ex:
        rows = [x for res in ex.map(fit_mouse, list(EV.items())) for x in res]
    F = pd.DataFrame(rows)
    F["group"] = F.mouse_key.map(grp)
    F.to_csv(D / "m04_fits.csv", index=False)

    # traces of the fitted models (value, RPE, CD, CE, M per drop; value and p(lick) per bin)
    TR = {}
    for m, A in EV.items():
        arr = arrays(A)
        for nm in ORDER:
            row = F[(F.mouse_key == m) & (F.model == nm)].iloc[0]
            code, lb, _ = SPECS[nm]
            x = np.array(row.params)
            nb = len(lb)
            TR[(m, nm)] = run(code, x[:nb], x[nb:], *arr, True)[2]
    pickle.dump(TR, open(D / "m04_traces.pkl", "wb"))

    # model comparison per group (6099_white left out: died on day 13, so no switch day in its data)
    out = {"bins_per_mouse": F.groupby("mouse_key").n_bins.first().to_dict()}
    K = F[F.mouse_key != "6099_white"]
    for g in ("Active", "Passive"):
        G = K[K.group == g]
        res = {}
        for metric in ("negLL", "BIC"):
            W = G.pivot(index="mouse_key", columns="model", values=metric)[ORDER]
            chg = W.sub(W.min(axis=1), axis=0)
            vs = {}
            ps = []
            for nm in ORDER[1:]:
                d = W[nm] - W["mRPE"]
                t = stats.ttest_rel(W[nm], W["mRPE"])
                w = stats.wilcoxon(W[nm], W["mRPE"]) if (d != 0).any() else None
                vs[nm] = dict(mean_diff=round(float(d.mean()), 2), n_mRPE_better=int((d > 0).sum()),
                              p_t=float(t.pvalue), p_wilcoxon_exact=float(w.pvalue) if w else None)
                ps.append(t.pvalue)
            for nm, a in zip(ORDER[1:], holm(ps)):
                vs[nm]["p_t_holm"] = float(a)
            res[metric] = dict(mean_change_from_best=chg.mean().round(2).to_dict(),
                               sem=(chg.std(ddof=1) / np.sqrt(len(chg))).round(2).to_dict(),
                               group_sum=W.sum().round(1).to_dict(), n_best=W.idxmin(axis=1).value_counts().to_dict(),
                               vs_mRPE=vs)
        Pm = G[G.model == "mRPE"].set_index("mouse_key").params.apply(pd.Series)
        Pm.columns = ["alpha", "d_CD", "e_CE", "b0", "beta", "g_time", "g_first3min"]
        res["mRPE_params"] = Pm.round(4).to_dict(orient="index")
        Pr = G[G.model == "RPE"].set_index("mouse_key").params.apply(pd.Series)
        Pr.columns = ["alpha", "b0", "beta", "g_time", "g_first3min"]
        res["RPE_params"] = Pr.round(4).to_dict(orient="index")
        res["mRPE_steps_capped"] = G[G.model == "mRPE"].set_index("mouse_key").n_steps_capped.to_dict()
        out[g] = res
    json.dump(out, open(D / "m04_summary.json", "w"), indent=1)
    for g in ("Active", "Passive"):
        print(g, "BIC change from best:", out[g]["BIC"]["mean_change_from_best"])
        print(g, "BIC n_best:", out[g]["BIC"]["n_best"])
        print(g, "negLL change from best:", out[g]["negLL"]["mean_change_from_best"])
        print(g, "vs mRPE (BIC):", {k: (v["mean_diff"], v["n_mRPE_better"], round(v["p_t_holm"], 4))
                                    for k, v in out[g]["BIC"]["vs_mRPE"].items()})
        print(g, "mRPE params:", out[g]["mRPE_params"])
        print(g, "capped steps:", out[g]["mRPE_steps_capped"])
