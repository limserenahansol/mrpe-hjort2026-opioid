"""Step 11: checks asked for by the independent slide review (2026-11-17).

(1) Model-free switch tests with unrounded values: change in p(lick) vs the day before, whole session and first 9 min,
    one-sample t per group with 95% CI, Active vs Passive (Welch t and Mann-Whitney), and days 15/16 vs 13.
(2) Group-level model recovery: simulate licking from the group fits (shared learning rates, each mouse's own readout)
    of RPE 2alpha, RPE, meta-RPE and today's outcome; refit those four at the group level; count which one BIC picks.
(3) Leave-one-mouse-out cross-validation (LOMO): learning rates fitted on the other mice, readout fitted on the held-out
    mouse; held-out negLL per mouse and model. This does not treat 30-s bins as independent evidence for the learning
    rates (they are never fitted on the held-out mouse).
Output: data/m11_modelfree.json, data/m11_group_recovery.json, data/m11_lomo_cv.json
"""
import itertools
import json
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.optimize import minimize

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m04_fit_user_data as m4  # noqa: E402
import m05_recovery_and_group_fits as m5  # noqa: E402

D = HERE / "data"
G5 = json.load(open(D / "m05_group_fits.json"))
MODELS = ["RPE2a", "RPE", "mRPE", "today's outcome"]
CV_MODELS = ["RPE2a", "RPE", "mRPE", "RPEt-1", "today's outcome", "no learning"]


# ------------------------------------------------------------------ (1) model-free
def modelfree():
    X = pd.read_csv(D / "m06_switch_blocks.csv")
    X = X[X.mouse != "6099_white"]
    out = {}

    def tci(v):
        v = np.asarray(v, float)
        t = stats.ttest_1samp(v, 0)
        h = stats.t.ppf(0.975, len(v) - 1) * v.std(ddof=1) / np.sqrt(len(v))
        return dict(mean=float(v.mean()), lo=float(v.mean() - h), hi=float(v.mean() + h), p=float(t.pvalue), n=len(v))

    for d0, d1 in ((13, 14), (13, 15), (13, 16), (16, 17), (16, 18)):
        key = f"day{d1}_vs_{d0}"
        out[key] = {}
        per = {}
        for g in ("Active", "Passive"):
            W = X[(X.group == g) & X.day.isin([d0, d1])].pivot_table(index="mouse", columns=["day", "block"], values="p")
            whole = W[d1].mean(axis=1) - W[d0].mean(axis=1)
            first9 = W[d1][[0, 1, 2]].mean(axis=1) - W[d0][[0, 1, 2]].mean(axis=1)
            per[g] = dict(whole=whole, first9=first9)
            out[key][g] = dict(whole=tci(whole), first9=tci(first9))
        for win in ("whole", "first9"):
            a, p = per["Active"][win], per["Passive"][win]
            out[key][f"Active_vs_Passive_{win}"] = dict(
                welch_p=float(stats.ttest_ind(a, p, equal_var=False).pvalue),
                mannwhitney_p=float(stats.mannwhitneyu(a, p, alternative="two-sided").pvalue))
    json.dump(out, open(D / "m11_modelfree.json", "w"), indent=1)
    return out


# ------------------------------------------------------------------ helpers on arbitrary choice arrays
def nll_mice(code, lp, arrs):
    tot = 0.0
    for typ, r, c, tau, early, use, dayout in arrs:
        q = m5.q_at_bins(code, np.asarray(lp, float), typ, r, dayout)
        sel = (typ == 1) & (use == 1)
        cols = [np.ones(sel.sum()), tau[sel], early[sel]] + ([] if code == 6 else [q[sel]])
        tot += m5.logit_nll(np.column_stack(cols).astype(np.float64), c[sel].astype(np.float64))
    return tot


def fit_group_arrs(name, arrs):
    code, lb, grid = m4.SPECS[name]
    n = sum(int(((a[0] == 1) & (a[5] == 1)).sum()) for a in arrs)
    k_read = (3 if code == 6 else 4) * len(arrs)
    if not lb:
        nll = nll_mice(code, [], arrs)
        return dict(model=name, negLL=nll, k=k_read, BIC=k_read * np.log(n) + 2 * nll, params=[])
    best = (np.inf, None)
    for g in itertools.product(*grid):
        res = minimize(lambda x: nll_mice(code, x, arrs), np.array(g, float), bounds=lb, method="L-BFGS-B")
        if res.fun < best[0]:
            best = (float(res.fun), res.x.copy())
    k = len(lb) + k_read
    return dict(model=name, negLL=best[0], k=k, BIC=k * np.log(n) + 2 * best[0], params=best[1].tolist())


def readout_w(code, lp, a):
    typ, r, c, tau, early, use, dayout = a
    q = m5.q_at_bins(code, np.asarray(lp, float), typ, r, dayout)
    sel = (typ == 1) & (use == 1)
    Xa = np.column_stack([np.ones(len(typ)), tau, early, q])
    w = np.zeros(4)
    Xs, y = Xa[sel], c[sel].astype(float)
    for _ in range(60):
        p = 1 / (1 + np.exp(-(Xs @ w)))
        H = (Xs.T * (p * (1 - p))) @ Xs + 1e-6 * np.eye(4)
        st = np.linalg.solve(H, Xs.T @ (p - y) + 1e-6 * w)
        w -= st
        if np.abs(st).max() < 1e-8:
            break
    return Xa, sel, w


# ------------------------------------------------------------------ (2) group-level recovery
def recov_job(args):
    g, gen, seed = args
    rng = np.random.default_rng(seed)
    code = m4.SPECS[gen][0]
    lp = G5[g][gen]["params"]
    arrs = []
    for m in m5.MICE[g]:
        a = m5.ARR[m]
        Xa, sel, w = readout_w(code, lp, a)
        cs = a[2].copy()
        p = 1 / (1 + np.exp(-(Xa[sel] @ w)))
        cs[sel] = (rng.random(sel.sum()) < p).astype(np.int64)
        arrs.append((a[0], a[1], cs, a[3], a[4], a[5], a[6]))
    fits = {nm: fit_group_arrs(nm, arrs) for nm in MODELS}
    best = min(fits, key=lambda k: fits[k]["BIC"])
    return dict(group=g, gen=gen, seed=seed, best=best, **{f"BIC_{k}": v["BIC"] for k, v in fits.items()})


# ------------------------------------------------------------------ (3) leave-one-mouse-out CV
def cv_job(args):
    g, name, held = args
    code = m4.SPECS[name][0]
    train = [m5.ARR[m] for m in m5.MICE[g] if m != held]
    f = fit_group_arrs(name, train)
    test_nll = nll_mice(code, f["params"], [m5.ARR[held]])
    return dict(group=g, model=name, held=held, test_negLL=float(test_nll), params=f["params"])


if __name__ == "__main__":
    MF = modelfree()
    for k, v in MF.items():
        print(k, {g: (round(v[g]["whole"]["mean"], 3), round(v[g]["whole"]["p"], 4), round(v[g]["first9"]["mean"], 4),
                      round(v[g]["first9"]["lo"], 3), round(v[g]["first9"]["hi"], 3)) for g in ("Active", "Passive")},
              "AvsP whole", round(v["Active_vs_Passive_whole"]["welch_p"], 4), "first9", round(v["Active_vs_Passive_first9"]["welch_p"], 4))
    jobs = [(g, gen, 7000 + 100 * i + s) for i, (g, gen) in enumerate(itertools.product(("Active", "Passive"), MODELS))
            for s in range(50)]
    with ProcessPoolExecutor(max_workers=30) as ex:
        R = pd.DataFrame(list(ex.map(recov_job, jobs)))
    R.to_csv(D / "m11_group_recovery_sets.csv", index=False)
    conf = {g: R[R.group == g].groupby("gen").best.value_counts(normalize=True).round(3).unstack(fill_value=0)
            .to_dict(orient="index") for g in ("Active", "Passive")}
    json.dump(dict(confusion=conf, n_sets=50, refit_models=MODELS), open(D / "m11_group_recovery.json", "w"), indent=1)
    print(json.dumps(conf, indent=1))
    cjobs = [(g, nm, m) for g in ("Active", "Passive") for nm in CV_MODELS for m in m5.MICE[g]]
    with ProcessPoolExecutor(max_workers=30) as ex:
        C = pd.DataFrame(list(ex.map(cv_job, cjobs)))
    C.to_csv(D / "m11_lomo_cv_rows.csv", index=False)
    out = {}
    for g in ("Active", "Passive"):
        W = C[C.group == g].pivot(index="held", columns="model", values="test_negLL")[CV_MODELS]
        tot = W.sum()
        ref = tot.idxmin()
        res = dict(total_test_negLL=tot.round(2).to_dict(), best=ref,
                   delta_from_best=(tot - tot.min()).round(2).to_dict(), vs_best={})
        for nm in CV_MODELS:
            if nm == ref:
                continue
            dlt = W[nm] - W[ref]
            res["vs_best"][nm] = dict(sum=round(float(dlt.sum()), 2), n_mice_best_better=int((dlt > 0).sum()),
                                      n=int(len(dlt)), p_wilcoxon=float(stats.wilcoxon(W[nm], W[ref]).pvalue)
                                      if (dlt != 0).any() else None)
        out[g] = res
    json.dump(out, open(D / "m11_lomo_cv.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


def asymmetry():
    """Active (and Passive): size of the whole-session change on the first day after each switch, paired by mouse."""
    X = pd.read_csv(D / "m06_switch_blocks.csv")
    X = X[X.mouse != "6099_white"]
    out = {}
    for g in ("Active", "Passive"):
        W = X[X.group == g].groupby(["mouse", "day"]).p.mean().unstack()
        d14 = (W[14] - W[13]).abs()
        d17 = (W[17] - W[16]).abs()
        out[g] = dict(n=int(len(W)), abs_day14=float(d14.mean()), abs_day17=float(d17.mean()),
                      n_day17_larger=int((d17 > d14).sum()), p_t=float(stats.ttest_rel(d17, d14).pvalue),
                      p_wilcoxon=float(stats.wilcoxon(d17, d14).pvalue))
    return out


if __name__ == "__main__" and "--asym" in sys.argv:
    A = asymmetry()
    MFj = json.load(open(D / "m11_modelfree.json"))
    MFj["asymmetry_whole_session_abs_change_day17_vs_day14"] = A
    json.dump(MFj, open(D / "m11_modelfree.json", "w"), indent=1)
    print(json.dumps(A, indent=1))
