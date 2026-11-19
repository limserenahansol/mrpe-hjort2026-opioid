"""Step 5: (a) can this design tell mRPE from simpler models at all? (b) group-level fits with shared learning rates.

(a) Model recovery, per mouse as in the paper: simulate licking (yes/no per 30-s bin) from each mouse's fitted
    mRPE, RPE and 'today's outcome' models, keeping the mouse's real drop sequence (value depends only on drops,
    so this is exact for the model), 20 data sets per mouse x model; refit the three models; count which one BIC picks.
    Also: do the fitted CD/CE rates come back (d < e)?
(b) Group fits: learning-rule parameters shared by all mice of a group, readout (b0, beta, time, first 3 min) per
    mouse. The readout is a logistic regression, solved exactly for every candidate learning rule (profile likelihood),
    so only the 1-4 learning parameters are searched. BIC with n = all bins of the group. Bootstrap over mice (2000)
    for the CD and CE rates and their difference.
Output: data/m05_recovery.json, data/m05_group_fits.json
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
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m04_fit_user_data as m4  # noqa: E402

D = HERE / "data"
EV = pickle.load(open(D / "m04_events.pkl", "rb"))
F = pd.read_csv(D / "m04_fits.csv")
F["params"] = F.params.apply(lambda s: np.array(json.loads(s)))
GRP = F.groupby("mouse_key").group.first().to_dict()
MICE = {g: [m for m in EV if GRP[m] == g and m != "6099_white"] for g in ("Active", "Passive")}
ARR = {m: m4.arrays(A) for m, A in EV.items()}


# ------------------------------------------------------------------ (a) model recovery
def recov_one(args):
    m, seed = args
    rng = np.random.default_rng(seed)
    typ, r, c, tau, early, use, dayout = ARR[m]
    out = []
    for gen in ("mRPE", "RPE", "today's outcome"):
        row = F[(F.mouse_key == m) & (F.model == gen)].iloc[0]
        code, lb, _ = m4.SPECS[gen]
        nb = len(lb)
        p = m4.run(code, row.params[:nb], row.params[nb:], typ, r, c, tau, early, use, dayout, True)[2][:, 5]
        for k in range(20):
            cs = c.copy()
            sel = typ == 1
            cs[sel] = (rng.random(sel.sum()) < p[sel]).astype(np.int64)
            arr = (typ, r, cs, tau, early, use, dayout)
            fits = {nm: m4.fit_one(nm, arr) for nm in ("mRPE", "RPE", "today's outcome")}
            best = min(fits, key=lambda nm: fits[nm]["BIC"])
            out.append(dict(mouse=m, gen=gen, rep=k, best=best,
                            dBIC_mRPE_minus_best_other=fits["mRPE"]["BIC"] - min(fits["RPE"]["BIC"], fits["today's outcome"]["BIC"]),
                            d_hat=fits["mRPE"]["params"][1], e_hat=fits["mRPE"]["params"][2],
                            d_true=row.params[1] if gen == "mRPE" else np.nan,
                            e_true=row.params[2] if gen == "mRPE" else np.nan))
    return out


# ------------------------------------------------------------------ (b) group fits, profile likelihood
@njit(cache=True)
def q_at_bins(code, lp, typ, r, dayout):
    """value seen by each bin under a learning rule (same rules as m04.run)."""
    Q = 0.0
    md = 0.0
    me = 0.0
    dprev = 0.0
    vy = 0.0
    vn = 0.0
    a_ph = 0.05
    a_mk = 0.1
    out = np.zeros(typ.shape[0])
    for i in range(typ.shape[0]):
        if typ[i] == 1:
            if code == 2 or code == 3:
                out[i] = vy - vn
            elif code == 7:
                out[i] = dayout[i]
            else:
                out[i] = Q
            continue
        if code == 0:
            delta = r[i] - Q
            md = (1 - lp[1]) * md + lp[1] * delta
            me = (1 - lp[2]) * me + lp[2] * delta
            M = (-md + abs(md)) / 2 + (me + abs(me)) / 2
            step = min(lp[0] * abs(M), 1.0)
            Q = Q + step * delta
        elif code == 1:
            delta = r[i] - Q
            Q = Q + (lp[0] if delta <= 0 else lp[1]) * delta
        elif code == 4:
            delta = r[i] - Q
            Q = Q + lp[0] * delta * dprev
            dprev = abs(delta)
        elif code == 5:
            Q = Q + lp[0] * (r[i] - Q)
        elif code == 2:
            vnet = vy - vn
            R = r[i] - (vy - vn)
            dvy = 0.0
            dvn = 0.0
            if R > 0:
                dvy = a_ph * lp[0] * r[i]
            elif R < 0:
                dvn = a_ph * lp[1] * abs(R)
            a_ph = lp[2] * abs(r[i] - vnet) + (1 - lp[2]) * a_ph
            vy = min(vy + dvy, 1.0)
            vn = max(vn + dvn, 0.0)
        elif code == 3:
            R = r[i] - (vy - vn)
            dvy = 0.0
            dvn = 0.0
            if R > 0:
                dvy = a_mk * lp[0] * (1 - vy + vn) * abs(R)
                da = -lp[2] * abs(r[i] - vy + vn)
            elif R < 0:
                dvn = a_mk * lp[1] * (1 - vn + vy) * abs(R)
                da = -lp[3] * abs(R * abs(-vn + vy))
            else:
                da = 0.0
            vy = vy + dvy
            vn = vn + dvn
            a_mk = min(max(a_mk - da, 0.05), 1.0)
    return out


@njit(cache=True)
def logit_nll(X, y):
    """max-likelihood logistic regression by Newton steps (tiny ridge for stability); returns min negLL."""
    k = X.shape[1]
    w = np.zeros(k)
    for it in range(60):
        z = X @ w
        p = 1.0 / (1.0 + np.exp(-z))
        g = X.T @ (p - y) + 1e-6 * w
        Wd = p * (1 - p)
        H = (X.T * Wd) @ X + 1e-6 * np.eye(k)
        step = np.linalg.solve(H, g)
        w = w - step
        if np.max(np.abs(step)) < 1e-8:
            break
    z = X @ w
    nll = 0.0
    for i in range(y.shape[0]):
        nll += np.log1p(np.exp(-abs(z[i]))) + max(-z[i], 0.0) if y[i] == 1 else np.log1p(np.exp(-abs(z[i]))) + max(z[i], 0.0)
    return nll


def group_nll(code, lp, mice, no_q=False):
    tot = 0.0
    for m in mice:
        typ, r, c, tau, early, use, dayout = ARR[m]
        q = q_at_bins(code, np.asarray(lp, float), typ, r, dayout)
        sel = (typ == 1) & (use == 1)
        cols = [np.ones(sel.sum()), tau[sel], early[sel]] + ([] if no_q else [q[sel]])
        X = np.column_stack(cols).astype(np.float64)
        tot += logit_nll(X, c[sel].astype(np.float64))
    return tot


def fit_group(name, mice):
    code, lb, grid = m4.SPECS[name]
    n = sum(int(((ARR[m][0] == 1) & (ARR[m][5] == 1)).sum()) for m in mice)
    k_read = (3 if code == 6 else 4) * len(mice)
    if not lb:
        nll = group_nll(code, [], mice, no_q=(code == 6))
        return dict(model=name, negLL=nll, k=k_read, n=n, BIC=k_read * np.log(n) + 2 * nll, params=[])
    best = (np.inf, None)
    for g in itertools.product(*grid):
        res = minimize(lambda x: group_nll(code, x, mice), np.array(g, float), bounds=lb, method="L-BFGS-B")
        if res.fun < best[0]:
            best = (float(res.fun), res.x.copy())
    k = len(lb) + k_read
    return dict(model=name, negLL=best[0], k=k, n=n, BIC=k * np.log(n) + 2 * best[0], params=best[1].tolist())


def boot_one(args):
    g, seed = args
    rng = np.random.default_rng(seed)
    ms = list(rng.choice(MICE[g], size=len(MICE[g]), replace=True))
    f = fit_group("mRPE", ms)
    fr = fit_group("RPE", ms)
    return dict(group=g, seed=seed, alpha=f["params"][0], d=f["params"][1], e=f["params"][2],
                dBIC_RPE_minus_mRPE=fr["BIC"] - f["BIC"])


if __name__ == "__main__":
    # (b) group fits first (fast)
    G = {}
    for g in ("Active", "Passive"):
        with ProcessPoolExecutor(max_workers=8) as ex:
            fits = list(ex.map(fit_group, m4.ORDER, [MICE[g]] * len(m4.ORDER)))
        best = min(f["BIC"] for f in fits)
        G[g] = {f["model"]: dict(f, dBIC_from_best=round(f["BIC"] - best, 2)) for f in fits}
        print(g, {k: (round(v["negLL"], 1), v["dBIC_from_best"], [round(x, 4) for x in v["params"]]) for k, v in G[g].items()})
    with ProcessPoolExecutor(max_workers=30) as ex:
        B = pd.DataFrame(list(ex.map(boot_one, [(g, s) for g in ("Active", "Passive") for s in range(1000)])))
    G["bootstrap_over_mice"] = {}
    for g in ("Active", "Passive"):
        b = B[B.group == g]
        G["bootstrap_over_mice"][g] = dict(
            n=int(len(b)), d_CI=b.d.quantile([.025, .5, .975]).round(4).tolist(), e_CI=b.e.quantile([.025, .5, .975]).round(4).tolist(),
            frac_d_less_than_e=float((b.d < b.e).mean()), dBIC_RPE_minus_mRPE_CI=b.dBIC_RPE_minus_mRPE.quantile([.025, .5, .975]).round(2).tolist(),
            frac_mRPE_better_BIC=float((b.dBIC_RPE_minus_mRPE > 0).mean()))
    print(json.dumps(G["bootstrap_over_mice"], indent=1))
    json.dump(G, open(D / "m05_group_fits.json", "w"), indent=1)

    # (a) model recovery
    jobs = [(m, 1000 + i) for i, m in enumerate(MICE["Active"] + MICE["Passive"])]
    with ProcessPoolExecutor(max_workers=14) as ex:
        R = pd.DataFrame([x for res in ex.map(recov_one, jobs) for x in res])
    R["group"] = R.mouse.map(GRP)
    conf = {g: R[R.group == g].groupby("gen").best.value_counts(normalize=True).round(3).unstack(fill_value=0).to_dict(orient="index")
            for g in ("Active", "Passive")}
    rec = R[R.gen == "mRPE"]
    out = dict(confusion_BIC=conf, n_sets_per_mouse_and_model=20,
               mRPE_generated_d_less_than_e_true=float((rec.d_true < rec.e_true).mean()),
               mRPE_generated_d_less_than_e_recovered=float((rec.d_hat < rec.e_hat).mean()),
               corr_d=float(np.corrcoef(rec.d_true, rec.d_hat)[0, 1]) if rec.d_true.std() > 0 else None,
               corr_e=float(np.corrcoef(rec.e_true, rec.e_hat)[0, 1]) if rec.e_true.std() > 0 else None)
    json.dump(out, open(D / "m05_recovery.json", "w"), indent=1)
    R.to_csv(D / "m05_recovery_sets.csv", index=False)
    print(json.dumps(out, indent=1))
