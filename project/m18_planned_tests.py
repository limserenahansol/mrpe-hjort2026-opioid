"""Step 18: the pre-specified tests for the next cohort (plan: ANALYSIS_PLAN_next_cohort_20261126.md).

Written on 2026-11-26, before any data of the next cohort were seen. The input is one data folder per cohort holding
the m01 outputs `licks_all.csv` and `sessions.csv`. On the Dec-2025 cohort this script must reproduce the m11
per-mouse changes, Welch p values and the asymmetry test exactly (asserted), which checks the code.

Usage (from this folder):
    python m18_planned_tests.py                  -> Dec-2025 cohort (data/) only: reference values + sample size
                                                    -> data/m18_planned_tests_dec2025.json
    python m18_planned_tests.py NEW_DIR [NAME]   -> the new cohort alone (primary + secondary tests, Holm) and pooled
                                                    with Dec-2025 (cohort as a factor) -> NEW_DIR/m18_planned_tests.json
    python m18_planned_tests.py --selftest       -> pooled code path on Dec-2025 + an identical copy (no output file)
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.stats.power import TTestIndPower, TTestPower

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
D0 = HERE / "data"
BW, NB = 30.0, 30          # p(lick) = share of the first 30 x 30-s bins (0-15 min) with >= 1 lick, as in m04 / m13
MIN_END = 840.0            # a session counts only if its controller log ends at >= 14 min (else: missing day)
M0, W1, W2, W3, R1 = 13, 14, 15, 16, 17   # last Equal Access day, Withdrawal days 1-3, first Re-exposure day
M1 = 12                    # the Equal Access day before M0 (only for the S2 sensitivity reference)
ALPHA = 0.05
NS = (4, 5, 6, 7, 8, 10, 12)
RS = (0.3, 0.5, 0.7)


# ------------------------------------------------------------------ outcome measure
def plick(data_dir):
    """one row per PR session: p(lick), session end, valid flag"""
    L = pd.read_csv(Path(data_dir) / "licks_all.csv")
    L = L[L.kind == "lick"]
    S = pd.read_csv(Path(data_dir) / "sessions.csv")
    S = S[S.paradigm == "ProgressiveRatioExperiment"]
    lk = {k: g.t.values for k, g in L.groupby(["mouse_key", "day"])}
    rows = []
    for r in S.itertuples():
        t = lk.get((r.mouse_key, r.day), np.array([]))
        b = np.unique(np.floor(t[(t >= 0) & (t < NB * BW)] / BW).astype(int))
        rows.append(dict(mouse=r.mouse_key, group=r.group, day=int(r.day), end=float(r.end), p=len(b) / NB,
                         valid=bool(r.end >= MIN_END)))
    return pd.DataFrame(rows)


def table(data_dir):
    X = plick(data_dir)
    P = X[X.valid].pivot(index="mouse", columns="day", values="p")
    grp = X.groupby("mouse").group.first()
    short = X[~X.valid & X.day.isin([M0, W1, W2, W3, R1])][["mouse", "day", "end"]].to_dict("records")
    return P, grp.reindex(P.index), X, short


# ------------------------------------------------------------------ tests
def onesample(x):
    x = np.asarray(x, float)
    t = stats.ttest_1samp(x, 0)
    ci = t.confidence_interval(0.95)
    sd = float(x.std(ddof=1))
    return dict(n=int(len(x)), mean=float(x.mean()), lo=float(ci.low), hi=float(ci.high), p=float(t.pvalue), sd=sd,
                dz=float(x.mean() / sd), n_pos=int((x > 0).sum()))


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    t = stats.ttest_ind(a, b, equal_var=False)
    ci = t.confidence_interval(0.95)
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    g = (a.mean() - b.mean()) / sp * (1 - 3 / (4 * (len(a) + len(b)) - 9))      # Hedges' g
    return dict(n_a=int(len(a)), n_b=int(len(b)), mean_a=float(a.mean()), mean_b=float(b.mean()),
                diff=float(a.mean() - b.mean()), lo=float(ci.low), hi=float(ci.high), p=float(t.pvalue), g=float(g))


def contrasts(P, grp):
    """per-mouse values each test uses (a mouse enters a test only if all of that test's days are valid)"""
    A, Pa = P[grp == "Active"], P[grp == "Passive"]
    lag = (A[W1] - (A[W2] + A[W3]) / 2).dropna()                      # P1: = change on W1 minus mean change on W2-W3
    d1a, d1p = (A[W1] - A[M0]).dropna(), (Pa[W1] - Pa[M0]).dropna()  # S2: change on the first water day
    r1 = (A[R1] - A[W3]).dropna()                                     # S3: change on the first morphine day
    asym = ((A[R1] - A[W3]).abs() - (A[W1] - A[M0]).abs()).dropna()   # S4: |R1 change| - |W1 change|, paired
    sa = (A[W1] - (A[M1] + A[M0]) / 2).dropna()                       # S2 sensitivity: reference = mean of M1, M0
    sp = (Pa[W1] - (Pa[M1] + Pa[M0]) / 2).dropna()
    return dict(lag=lag, d1a=d1a, d1p=d1p, r1=r1, asym=asym, s2sens_a=sa, s2sens_p=sp,
                abs_w1=(A[W1] - A[M0]).abs()[asym.index], abs_r1=(A[R1] - A[W3]).abs()[asym.index])


def run_tests(P, grp):
    C = contrasts(P, grp)
    S4 = onesample(C["asym"])
    S4.update(abs_w1=float(C["abs_w1"].mean()), abs_r1=float(C["abs_r1"].mean()),
              p_wilcoxon=float(stats.wilcoxon(C["abs_r1"], C["abs_w1"]).pvalue))
    return dict(P1_lag=onesample(C["lag"]), S2_W1_active_vs_passive=welch(C["d1a"], C["d1p"]),
                S2_W1_active=onesample(C["d1a"]), S2_W1_passive=onesample(C["d1p"]),
                S2_sensitivity_ref_M1_M0=welch(C["s2sens_a"], C["s2sens_p"]),
                S3_R1_active=onesample(C["r1"]), S4_asymmetry=S4), C


def pooled(Cs):
    """cohort as a sum-coded factor: the intercept is the unweighted mean of the cohort means"""
    # (no local variable may be called C here: patsy would take it for its C() function)
    def stack(key):
        rows = []
        for coh, cv in Cs.items():
            for m, v in cv[key].items():
                rows.append(dict(cohort=coh, mouse=m, y=float(v)))
        return pd.DataFrame(rows)

    def icpt(df):
        f = smf.ols("y ~ C(cohort, Sum)", df).fit()
        lo, hi = f.conf_int().loc["Intercept"]
        cf = [k for k in f.params.index if k.startswith("C(cohort")]
        return dict(n=int(len(df)), mean=float(f.params["Intercept"]), lo=float(lo), hi=float(hi),
                    p=float(f.pvalues["Intercept"]), cohort_coef=float(f.params[cf[0]]), cohort_p=float(f.pvalues[cf[0]]))

    out = dict(S1_lag=icpt(stack("lag")), S3_R1_active=icpt(stack("r1")), S4_asymmetry=icpt(stack("asym")))
    rows = []
    for coh, cv in Cs.items():
        rows += [dict(cohort=coh, group="Active", y=float(v)) for v in cv["d1a"]]
        rows += [dict(cohort=coh, group="Passive", y=float(v)) for v in cv["d1p"]]
    df = pd.DataFrame(rows)
    f = smf.ols("y ~ C(group, Treatment('Passive')) + C(cohort, Sum)", df).fit(cov_type="HC3")
    k = "C(group, Treatment('Passive'))[T.Active]"
    lo, hi = f.conf_int().loc[k]
    out["S2_W1_active_vs_passive"] = dict(n=int(len(df)), diff=float(f.params[k]), lo=float(lo), hi=float(hi),
                                          p=float(f.pvalues[k]), se="HC3")
    return out


# ------------------------------------------------------------------ sample size (planning d as in m17)
def planning_d():
    MF = json.load(open(D0 / "m11_modelfree.json"))
    A = {k: MF[k]["Active"]["whole"] for k in ("day14_vs_13", "day15_vs_13", "day16_vs_13", "day17_vs_16")}
    n = A["day14_vs_13"]["n"]
    sd = float(np.mean([(v["hi"] - v["lo"]) / 2 * np.sqrt(n) / stats.t.ppf(0.975, n - 1) for v in A.values()]))
    late = (A["day15_vs_13"]["mean"] + A["day16_vs_13"]["mean"]) / 2
    return (abs(late) - abs(A["day14_vs_13"]["mean"])) / sd, MF


def sample_size(ref, C):
    d, _ = planning_d()
    one, ind = TTestPower(), TTestIndPower()
    plan = {f"r={r}": dict(dz=float(d / np.sqrt(2 * (1 - r))),
                           power={f"n={n}": float(one.power(d / np.sqrt(2 * (1 - r)), nobs=n, alpha=ALPHA)) for n in NS},
                           n80=float(one.solve_power(d / np.sqrt(2 * (1 - r)), alpha=ALPHA, power=0.8)))
            for r in RS}
    A_ = C["lag"].index
    r_obs = float(np.corrcoef((C["d1a"][A_]).values, (C["d1a"][A_] - C["lag"]).values)[0, 1])  # W1 vs mean W2-W3 change
    p1, s2, s3, s4 = ref["P1_lag"], ref["S2_W1_active_vs_passive"], ref["S3_R1_active"], ref["S4_asymmetry"]
    obs = dict(P1=dict(dz=p1["dz"], r=r_obs, power_n6=float(one.power(p1["dz"], nobs=6, alpha=ALPHA)),
                       n80=float(one.solve_power(p1["dz"], alpha=ALPHA, power=0.8))),
               S2=dict(g=s2["g"], power_n6=float(ind.power(s2["g"], nobs1=6, alpha=ALPHA)),
                       n80=float(ind.solve_power(s2["g"], alpha=ALPHA, power=0.8))),
               S3=dict(dz=s3["dz"], power_n6=float(one.power(s3["dz"], nobs=6, alpha=ALPHA)),
                       n80=float(one.solve_power(s3["dz"], alpha=ALPHA, power=0.8))),
               S4=dict(dz=s4["dz"], power_n6=float(one.power(s4["dz"], nobs=6, alpha=ALPHA)),
                       n80=float(one.solve_power(s4["dz"], alpha=ALPHA, power=0.8))))
    return dict(planning_d=float(d), P1_within_mouse=plan, observed_dec2025=obs)


# ------------------------------------------------------------------ the Dec-2025 cohort must reproduce m11
def check_dec2025(P, grp, ref):
    _, MF = planning_d()
    A, Pa = P[grp == "Active"], P[grp == "Passive"]
    for key, d1, d0 in (("day14_vs_13", W1, M0), ("day15_vs_13", W2, M0), ("day16_vs_13", W3, M0), ("day17_vs_16", R1, W3)):
        for g, X in (("Active", A), ("Passive", Pa)):
            o, m = onesample((X[d1] - X[d0]).dropna()), MF[key][g]["whole"]
            assert o["n"] == m["n"] and all(abs(o[k] - m[k]) < 1e-9 for k in ("mean", "lo", "hi", "p")), (key, g)
        w = welch((A[d1] - A[d0]).dropna(), (Pa[d1] - Pa[d0]).dropna())["p"]
        assert abs(w - MF[key]["Active_vs_Passive_whole"]["welch_p"]) < 1e-9, key
    a = MF["asymmetry_whole_session_abs_change_day17_vs_day14"]["Active"]
    s4 = ref["S4_asymmetry"]
    assert abs(s4["p"] - a["p_t"]) < 1e-9 and abs(s4["abs_w1"] - a["abs_day14"]) < 1e-9 \
        and abs(s4["abs_r1"] - a["abs_day17"]) < 1e-9 and abs(s4["p_wilcoxon"] - a["p_wilcoxon"]) < 1e-9
    return "Dec-2025 reproduces m11 (changes, CIs, one-sample and Welch p, asymmetry test)"


def describe(X, P, grp, short):
    used = sorted(set(P.index))
    return dict(n_active=int((grp == "Active").sum()), n_passive=int((grp == "Passive").sum()), mice=used,
                min_end_days13_17=float(X[X.day.isin([M0, W1, W2, W3, R1])].end.min()), short_sessions=short)


def main():
    P0, g0, X0, sh0 = table(D0)
    ref, C0 = run_tests(P0, g0)
    msg = check_dec2025(P0, g0, ref)
    print(msg)
    out = dict(plan="ANALYSIS_PLAN_next_cohort_20261126.md", written="2026-11-26",
               rules=dict(bin_s=BW, n_bins=NB, min_end_s=MIN_END, days=dict(M0=M0, W1=W1, W2=W2, W3=W3, R1=R1)),
               check=msg, cohorts=dict(dec2025=describe(X0, P0, g0, sh0)), tests=dict(dec2025=ref))
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    if "--selftest" in sys.argv:
        pl = pooled({"a": C0, "b": C0})
        assert abs(pl["S1_lag"]["mean"] - ref["P1_lag"]["mean"]) < 1e-12 and abs(pl["S1_lag"]["cohort_coef"]) < 1e-12
        assert abs(pl["S3_R1_active"]["mean"] - ref["S3_R1_active"]["mean"]) < 1e-12
        print("selftest OK: pooled intercept = cohort mean, cohort coefficient = 0")
        return
    if not args:
        out["sample_size"] = sample_size(ref, C0)
        dest = D0 / "m18_planned_tests_dec2025.json"
    else:
        new = Path(args[0])
        name = args[1] if len(args) > 1 else new.name
        P1_, g1, X1, sh1 = table(new)
        res, C1 = run_tests(P1_, g1)
        pl = pooled({"dec2025": C0, name: C1})
        sec = [pl["S1_lag"]["p"], res["S2_W1_active_vs_passive"]["p"], res["S3_R1_active"]["p"], res["S4_asymmetry"]["p"]]
        rej, ph, _, _ = multipletests(sec, alpha=ALPHA, method="holm")
        p1 = res["P1_lag"]
        out["cohorts"][name] = describe(X1, P1_, g1, sh1)
        out["tests"][name] = res
        out["pooled"] = pl
        out["primary"] = dict(test="P1_lag, new cohort alone", p=p1["p"], mean=p1["mean"],
                              replicated=bool(p1["p"] < ALPHA and p1["mean"] > 0))
        out["secondary_holm"] = {k: dict(p=float(p), p_holm=float(q), significant=bool(s)) for k, p, q, s in
                                 zip(["S1_lag_pooled", "S2_W1_active_vs_passive", "S3_R1_active", "S4_asymmetry"],
                                     sec, ph, rej)}
        dest = new / "m18_planned_tests.json"
    json.dump(out, open(dest, "w", encoding="utf-8"), indent=1)
    p1 = ref["P1_lag"]
    print(f"Dec-2025 P1 lag: {p1['mean']:+.3f} (95% CI {p1['lo']:+.3f} to {p1['hi']:+.3f}), p = {p1['p']:.4f}, "
          f"dz = {p1['dz']:.2f}, n = {p1['n']}")
    print("wrote", dest)


if __name__ == "__main__":
    main()
