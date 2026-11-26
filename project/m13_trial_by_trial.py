"""Step 13 (asked 18 Nov 2026): trial-by-trial checks of the two switches, as in Hjort et al. 2026 Fig. 1d,e,m.

CORRECTION (18 Nov 2026, after the one independent review; outputs of this script are unchanged = v1 slides 7-8):
the posterior predictive check below draws every 30-s bin independently from the fitted p. Real licking varies about
twice as much from day to day and is autocorrelated, so its p values are TOO SMALL (anti-conservative), not
'conservative' as written below; a no-change placebo day pair (Passive day 12 -> 13) is rejected by every model.
The 'r2' columns are uncentred (1 - SSE / sum y^2), not the usual R^2. See m15_trial_checks_v2.py (calibrated check,
placebo, mouse-level tests, equal windows, leave-one-out) and m16 (slides v2).

Trial = one 30-s bin: the same 30 trials per 15-min session that every model was fitted on (m04, m05).
Measure = lick yes/no in the bin (the models' choice). Licking falls within every session, so each bin is compared
with the SAME bin on the last day before the switch:
    drug -> water: D(day, bin) = lick(day, bin) - lick(day 13, bin), days 14-16 = trials 0-89 from the switch
    water -> drug: D(day, bin) = lick(day, bin) - lick(day 16, bin), days 17-18 = trials 0-59 from the switch
Models: the group fits of m05 (learning rates shared in the group; readout b0, beta, time in session, first 3 min
refitted per mouse on days 11-18, exactly as m10). Every model updates value at each drop only; value is carried
over the night unchanged.

(A) Within day 14 or overnight?
    - whole sessions (no window choice): mean D on days 14, 15, 16, observed and predicted by each model;
      'lag' = observed minus predicted on day 14, 'night' = (D15 - D14) observed minus predicted.
    - day 14 first half (0-7.5 min) vs second half (7.5-15 min) and the night from day 14's second half to
      day 15's first half; 5-min windows as a sensitivity check; 2.5-min profile (descriptive).
    - water drops through day 14 (when the errors happened).
    - reference check: within-session fall on day 13 vs days 11-12 (morphine) and days 4-5 (Baseline, water).
(B) Fig. 1d,e: sigmoid D(x) = A / (1 + exp(-(x - x0) / s)); midpoint x0 = trials from the switch.
    Per mouse (as the paper) and for the group-mean curve.
(C) Fig. 1m: start = first trial at which the fitted curve has made 10% of its change = max(0, x0 - s ln 9)
    (the paper: 'decays more than 0.1 normalized licks'); mice vs each model's predicted licking through the
    identical pipeline (noise-free prediction: paired t-test, Holm) and by posterior predictive simulation.
Posterior predictive check: 2000 data sets per model (lick yes/no per bin drawn from the fitted model's p; the
mouse's real drops); every statistic computed the same way; p = (1 + #{|sim - centre| >= |obs - centre|}) / 2001
with centre = mean of the simulated values. Not refitted per data set: refitting can only shrink the observed
misfit, so these p values are conservative.
Sigmoid fit: exact least squares over a grid (x0 from -10 to window end + 10 in 0.5-trial steps, s from 0.3 to 40
trials in 40 log steps, A in closed form and bounded to |A| <= 1, the largest change yes/no licking can show),
identical for mice, models and simulations. A curve whose amplitude |A| < 0.10 has no clear change; its timing is
reported but flagged. Midpoints are clipped to the window (0 = at once, window end = not reached).
Output: data/m13_bins.csv, data/m13_per_mouse.csv, data/m13_summary.json, data/m13_log.txt
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m04_fit_user_data as m4  # noqa: E402
import m05_recovery_and_group_fits as m5  # noqa: E402

D = HERE / "data"
G5 = json.load(open(D / "m05_group_fits.json"))
MODELS = ["today's outcome", "RPE", "RPE2a", "mRPE", "PH", "eMack", "RPEt-1", "no learning"]
MAIN = ["today's outcome", "RPE", "RPE2a", "mRPE"]
DAYS = list(range(11, 19))                      # array index 0..7
DI = {d: i for i, d in enumerate(DAYS)}
NSIM = 2000
SEED = 20261118
LN9 = np.log(9.0)
AMIN = 0.10
AMAX = 1.0
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


def logit_fit(Xm, y):
    """same readout fit as m10.logit_fit (Newton, tiny ridge)"""
    w = np.zeros(Xm.shape[1])
    for _ in range(60):
        p = 1 / (1 + np.exp(-(Xm @ w)))
        H = (Xm.T * (p * (1 - p))) @ Xm + 1e-6 * np.eye(len(w))
        step = np.linalg.solve(H, Xm.T @ (p - y) + 1e-6 * w)
        w -= step
        if np.abs(step).max() < 1e-8:
            break
    return w


# ============================================================================ data: observed, predicted, drops
def build():
    OBS, PRED, DROPS, REF = {}, {k: {} for k in MODELS}, {}, {}
    group = {}
    for g in ("Active", "Passive"):
        for m in m5.MICE[g]:
            group[m] = g
            typ, r, c, tau, early, use, dayout = m5.ARR[m]
            A = m5.EV[m]
            sel = (typ == 1) & (use == 1)
            dd, bb = A.day.values[sel].astype(int), A.bin.values[sel].astype(int)
            assert sel.sum() == 240 and set(dd) == set(DAYS), m
            o = np.full((8, 30), np.nan)
            o[dd - 11, bb] = c[sel]
            assert not np.isnan(o).any()
            OBS[m] = o
            for k in MODELS:
                code = m4.SPECS[k][0]
                q = m5.q_at_bins(code, np.asarray(G5[g][k]["params"], float), typ, r, dayout)
                cols = [np.ones(len(typ)), tau, early] + ([] if k == "no learning" else [q])
                Xa = np.column_stack(cols)
                w = logit_fit(Xa[sel], c[sel].astype(float))
                pm = np.full((8, 30), np.nan)
                pm[dd - 11, bb] = expit(Xa[sel] @ w)
                PRED[k][m] = pm
            # drops per 30-s bin (bins 0-29; drops after 15 min counted as 'after')
            dr = A[(A.typ == 0) & A.day.isin(DAYS)]
            cnt = np.zeros((8, 31))
            for d, t in zip(dr.day.values.astype(int), dr.t.values):
                cnt[d - 11, min(int(t // 30), 30)] += 1
            DROPS[m] = cnt
            # reference days for the within-session check: all PR bins of days 3-5 and 11-13
            rb = A[(A.typ == 1) & A.day.isin([3, 4, 5, 11, 12, 13])]
            REF[m] = {int(d): s.sort_values("bin").c.values.astype(float) for d, s in rb.groupby("day")}
    return OBS, PRED, DROPS, REF, group


def lick_counts(mice):
    L = pd.read_csv(D / "licks_all.csv")
    L = L[(L.kind == "lick") & (L.paradigm == "ProgressiveRatioExperiment") & L.day.isin(DAYS) & L.mouse_key.isin(mice)]
    L = L.assign(b=(L.t // 30).astype(int))
    L = L[L.b < 30]
    out = {}
    for m in mice:
        a = np.zeros((8, 30))
        x = L[L.mouse_key == m]
        for (d, b), n in x.groupby(["day", "b"]).size().items():
            a[int(d) - 11, int(b)] = n
        out[m] = a
    return out


# ============================================================================ (A) statistics (vectorised)
def a_stats(C):
    """C: (..., 8, 30) lick yes/no or probability for days 11-18. Returns dict of (...) arrays."""
    ref = C[..., DI[13], :]
    D14, D15, D16 = (C[..., DI[d], :] - ref for d in (14, 15, 16))
    o = dict(day14=D14.mean(-1), day15=D15.mean(-1), day16=D16.mean(-1),
             d14_h1=D14[..., :15].mean(-1), d14_h2=D14[..., 15:].mean(-1),
             d15_h1=D15[..., :15].mean(-1), d15_h2=D15[..., 15:].mean(-1),
             d14_e5=D14[..., :10].mean(-1), d14_l5=D14[..., 20:].mean(-1), d15_e5=D15[..., :10].mean(-1))
    o["within14_halves"] = o["d14_h2"] - o["d14_h1"]
    o["night_halves"] = o["d15_h1"] - o["d14_h2"]
    o["within14_5min"] = o["d14_l5"] - o["d14_e5"]
    o["night_5min"] = o["d15_e5"] - o["d14_l5"]
    o["d15_minus_d14"] = o["day15"] - o["day14"]
    return o


A_KEYS = ["day14", "day15", "day16", "d14_h1", "d14_h2", "d15_h1", "within14_halves", "night_halves",
          "within14_5min", "night_5min", "d15_minus_d14"]


def profile25(C):
    """2.5-min windows (5 bins) of D for days 14-16 vs day 13 -> (..., 18)"""
    ref = C[..., DI[13], :]
    segs = []
    for d in (14, 15, 16):
        x = C[..., DI[d], :] - ref
        segs.append(x.reshape(x.shape[:-1] + (6, 5)).mean(-1))
    return np.concatenate(segs, -1)


# ============================================================================ (B, C) sigmoid by grid least squares
class Grid:
    def __init__(self, n):
        x = np.arange(n, dtype=float)
        x0 = np.arange(-10.0, n + 10.0 + 1e-9, 0.5)
        s = np.geomspace(0.3, 40.0, 40)
        X0, S = np.meshgrid(x0, s, indexing="ij")
        self.x0, self.s = X0.ravel(), S.ravel()
        self.B = expit((x[None, :] - self.x0[:, None]) / self.s[:, None])      # (G, n)
        self.BB = (self.B ** 2).sum(1)
        self.n = n

    def fit(self, Y, chunk=800):
        """Y: (N, n) trajectories -> dict of (N,) arrays (processed in chunks to bound memory)"""
        Y = np.atleast_2d(np.asarray(Y, float))
        out = {k: [] for k in ("A", "x0", "s", "start", "r2", "mid")}
        for i0 in range(0, Y.shape[0], chunk):
            Yc = Y[i0:i0 + chunk]
            DS = self.B @ Yc.T                                                # (G, N)
            yy = (Yc ** 2).sum(1)
            Aall = np.clip(DS / self.BB[:, None], -AMAX, AMAX)                 # |A| <= 1: D of yes/no data
            sse = yy[None, :] - 2 * Aall * DS + Aall ** 2 * self.BB[:, None]
            k = np.argmin(sse, 0)
            cols = np.arange(Yc.shape[0])
            A = Aall[k, cols]
            x0, s = self.x0[k], self.s[k]
            out["A"].append(A)
            out["x0"].append(x0)
            out["s"].append(s)
            out["start"].append(np.maximum(0.0, x0 - s * LN9))
            out["r2"].append(np.where(yy > 0, 1 - sse[k, cols] / np.where(yy > 0, yy, 1), np.nan))
            out["mid"].append(np.clip(x0, 0, self.n))
        return {k: np.concatenate(v) for k, v in out.items()}


GD, GE = Grid(90), Grid(60)


def traj_deg(C):
    """(..., 8, 30) -> (..., 90): days 14-16 minus day 13, bin by bin"""
    ref = C[..., DI[13], :]
    return np.concatenate([C[..., DI[d], :] - ref for d in (14, 15, 16)], -1)


def traj_enh(C):
    ref = C[..., DI[16], :]
    return np.concatenate([C[..., DI[d], :] - ref for d in (17, 18)], -1)


# ============================================================================ tests
def ttest1(x):
    x = np.asarray(x, float)
    t = stats.ttest_1samp(x, 0.0)
    se = x.std(ddof=1) / np.sqrt(len(x))
    h = stats.t.ppf(0.975, len(x) - 1) * se
    return dict(mean=float(x.mean()), lo=float(x.mean() - h), hi=float(x.mean() + h), p=float(t.pvalue),
                n=int(len(x)), n_neg=int((x < 0).sum()), n_pos=int((x > 0).sum()), signflip_p=signflip(x))


def signflip(x):
    """exact two-sided sign-flip permutation p for the mean (all 2^n sign patterns)"""
    x = np.asarray(x, float)
    n = len(x)
    pats = ((np.arange(2 ** n)[:, None] >> np.arange(n)) & 1) * 2 - 1
    sims = (pats * x).mean(1)
    return float((np.abs(sims) >= abs(x.mean()) - 1e-12).mean())


def paired(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    t = stats.ttest_rel(a, b)
    d = a - b
    w = stats.wilcoxon(a, b) if (d != 0).any() else None
    return dict(mean_a=float(a.mean()), mean_b=float(b.mean()), mean_diff=float(d.mean()), p_t=float(t.pvalue),
                p_wilcoxon=float(w.pvalue) if w is not None else None, n=int(len(a)), n_a_greater=int((d > 0).sum()))


def welch(a, b):
    t = stats.ttest_ind(a, b, equal_var=False)
    return dict(p=float(t.pvalue), mean_a=float(np.mean(a)), mean_b=float(np.mean(b)))


def holm(ps):
    return m4.holm(ps).tolist()


def ppc_p(sim, obs):
    sim = np.asarray(sim, float)
    sim = sim[np.isfinite(sim)]
    c = sim.mean()
    k = (np.abs(sim - c) >= abs(obs - c) - 1e-12).sum()
    return float((k + 1) / (len(sim) + 1))


# ============================================================================ main
def main():
    OBS, PRED, DROPS, REF, GROUP = build()
    MICE = {g: [m for m in OBS if GROUP[m] == g] for g in ("Active", "Passive")}
    NL = lick_counts(list(OBS))
    for m in OBS:          # lick yes/no must equal 'any lick counted'
        assert np.array_equal(OBS[m] > 0.5, NL[m] > 0), m
    rng = np.random.default_rng(SEED)
    S = {"meta": dict(trial="30-s bin", n_sim=NSIM, seed=SEED, amp_min=AMIN, models=MODELS,
                      mice={g: MICE[g] for g in MICE},
                      note="posterior predictive without refitting (conservative for misfit)")}
    rows_pm, rows_bins = [], []

    # ---------------------------------------------------------------- long table of bins (observed + predicted)
    for m in OBS:
        for d in DAYS:
            for b in range(30):
                row = dict(mouse=m, group=GROUP[m], day=d, bin=b, lick=int(OBS[m][DI[d], b]), n_licks=int(NL[m][DI[d], b]),
                           drops_in_bin=int(DROPS[m][DI[d], b]))
                for k in MODELS:
                    row["p_" + k] = float(PRED[k][m][DI[d], b])
                rows_bins.append(row)
    pd.DataFrame(rows_bins).to_csv(D / "m13_bins.csv", index=False)

    for g in ("Active", "Passive"):
        ms = MICE[g]
        Og = np.stack([OBS[m] for m in ms])                                    # (n, 8, 30)
        Pg = {k: np.stack([PRED[k][m] for m in ms]) for k in MODELS}
        res = {}

        # ---------------- water drops through day 14 (when the errors happened)
        dr = np.stack([DROPS[m][DI[14]] for m in ms])                         # (n, 31)
        cum = np.cumsum(dr[:, :30], 1)
        tot = dr.sum(1)
        res["day14_drops"] = dict(total_mean=float(tot.mean()), total_min=int(tot.min()), total_max=int(tot.max()),
                                  by_2p5min=[float(x) for x in cum[:, 4::5].mean(0)],
                                  frac_by_2p5min=[float(x) for x in (cum[:, 4::5] / tot[:, None]).mean(0)],
                                  after_15min_mean=float(dr[:, 30].mean()))
        dr15 = np.stack([DROPS[m][DI[15]] for m in ms])
        res["day15_drops_total_mean"] = float(dr15.sum(1).mean())

        # ---------------- (A) observed, predicted, posterior predictive
        ao = a_stats(Og)
        res["A_observed"] = {k: ttest1(ao[k]) for k in A_KEYS}
        res["A_observed_counts_log1p"] = {}
        Cg = np.log1p(np.stack([NL[m] for m in ms]))
        ac = a_stats(Cg)
        for k in ("day14", "day15", "d14_h1", "d14_h2", "d15_h1", "within14_halves", "night_halves"):
            res["A_observed_counts_log1p"][k] = ttest1(ac[k])
        res["profile25_observed"] = dict(mean=profile25(Og).mean(0).tolist(),
                                         sem=(profile25(Og).std(0, ddof=1) / np.sqrt(len(ms))).tolist())
        res["A_models"] = {}
        res["profile25_models"] = {}
        res["profile25_ppc_band"] = {}
        for k in MODELS:
            ap = a_stats(Pg[k])
            sims = (rng.random((NSIM,) + Pg[k].shape) < Pg[k][None]).astype(float)   # (NSIM, n, 8, 30)
            asim = a_stats(sims)
            mk = {}
            for key in A_KEYS:
                gm = asim[key].mean(1)                                          # group mean per simulated set
                mk[key] = dict(pred_mean=float(ap[key].mean()), sim_lo=float(np.percentile(gm, 2.5)),
                               sim_hi=float(np.percentile(gm, 97.5)), ppc_p=ppc_p(gm, ao[key].mean()),
                               obs_minus_pred=paired(ao[key], ap[key]))
            res["A_models"][k] = mk
            prof = profile25(Pg[k]).mean(0)
            psim = profile25(sims).mean(1)
            res["profile25_models"][k] = prof.tolist()
            res["profile25_ppc_band"][k] = dict(lo=np.percentile(psim, 2.5, 0).tolist(),
                                                hi=np.percentile(psim, 97.5, 0).tolist())
            # ---------------- (B, C) sigmoid timing under the model: noise-free and simulated
            fd_pred = GD.fit(traj_deg(Pg[k]))
            fe_pred = GE.fit(traj_enh(Pg[k]))
            res["A_models"][k]["sig_pred"] = dict(
                deg_mid=fd_pred["mid"].tolist(), deg_start=fd_pred["start"].tolist(), deg_A=fd_pred["A"].tolist(),
                enh_mid=fe_pred["mid"].tolist(), enh_start=fe_pred["start"].tolist(), enh_A=fe_pred["A"].tolist())
            tds = traj_deg(sims).reshape(-1, 90)
            tes = traj_enh(sims).reshape(-1, 60)
            fds, fes = GD.fit(tds), GE.fit(tes)
            n = len(ms)
            res["A_models"][k]["sig_sim"] = dict(
                deg_mid_gm=fds["mid"].reshape(NSIM, n).mean(1).tolist(),
                deg_start_gm=fds["start"].reshape(NSIM, n).mean(1).tolist(),
                enh_mid_gm=fes["mid"].reshape(NSIM, n).mean(1).tolist(),
                enh_start_gm=fes["start"].reshape(NSIM, n).mean(1).tolist())
            # group-mean curve under the model, fitted the same way as the mice's mean curve
            gcurve = GD.fit(traj_deg(sims).mean(1))
            res["A_models"][k]["sig_sim_groupcurve"] = dict(deg_mid=gcurve["mid"].tolist(),
                                                            deg_start=gcurve["start"].tolist())
            gcurve_e = GE.fit(traj_enh(sims).mean(1))
            res["A_models"][k]["sig_sim_groupcurve"].update(enh_mid=gcurve_e["mid"].tolist(),
                                                            enh_start=gcurve_e["start"].tolist())
            gp = GD.fit(traj_deg(Pg[k]).mean(0))
            gpe = GE.fit(traj_enh(Pg[k]).mean(0))
            res["A_models"][k]["sig_pred_groupcurve"] = dict(deg_mid=float(gp["mid"][0]), deg_start=float(gp["start"][0]),
                                                             enh_mid=float(gpe["mid"][0]), enh_start=float(gpe["start"][0]))
            log(g, k, "pred day14", round(float(ap["day14"].mean()), 3), "obs", round(float(ao["day14"].mean()), 3),
                "ppc p", round(mk["day14"]["ppc_p"], 4), "| night_halves pred", round(float(ap["night_halves"].mean()), 3),
                "obs", round(float(ao["night_halves"].mean()), 3), "ppc", round(mk["night_halves"]["ppc_p"], 4))

        # ---------------- (B) mice: per-mouse sigmoids + the group-mean curve
        fd, fe = GD.fit(traj_deg(Og)), GE.fit(traj_enh(Og))
        gd_, ge_ = GD.fit(traj_deg(Og).mean(0)), GE.fit(traj_enh(Og).mean(0))
        res["sig_obs"] = dict(
            deg=dict(A=fd["A"].tolist(), mid=fd["mid"].tolist(), x0=fd["x0"].tolist(), s=fd["s"].tolist(),
                     start=fd["start"].tolist(), r2=fd["r2"].tolist()),
            enh=dict(A=fe["A"].tolist(), mid=fe["mid"].tolist(), x0=fe["x0"].tolist(), s=fe["s"].tolist(),
                     start=fe["start"].tolist(), r2=fe["r2"].tolist()),
            groupcurve=dict(deg=dict(A=float(gd_["A"][0]), mid=float(gd_["mid"][0]), x0=float(gd_["x0"][0]), s=float(gd_["s"][0]),
                                     start=float(gd_["start"][0]), r2=float(gd_["r2"][0])),
                            enh=dict(A=float(ge_["A"][0]), mid=float(ge_["mid"][0]), x0=float(ge_["x0"][0]), s=float(ge_["s"][0]),
                                     start=float(ge_["start"][0]), r2=float(ge_["r2"][0]))))
        ok = (np.abs(fd["A"]) >= AMIN) & (np.abs(fe["A"]) >= AMIN)
        res["sig_tests"] = dict(
            deg_vs_enh_mid_all=paired(fd["mid"], fe["mid"]),
            deg_vs_enh_start_all=paired(fd["start"], fe["start"]),
            n_clear_both=int(ok.sum()),
            deg_vs_enh_mid_clear=paired(fd["mid"][ok], fe["mid"][ok]) if ok.sum() >= 3 else None)
        # bootstrap over mice for the group-curve midpoints (2000)
        bm, be, bs_, bse = [], [], [], []
        Td, Te = traj_deg(Og), traj_enh(Og)
        for _ in range(2000):
            ii = rng.integers(0, len(ms), len(ms))
            f1, f2 = GD.fit(Td[ii].mean(0)), GE.fit(Te[ii].mean(0))
            bm.append(f1["mid"][0]); be.append(f2["mid"][0]); bs_.append(f1["start"][0]); bse.append(f2["start"][0])
        bm, be, bs_, bse = map(np.array, (bm, be, bs_, bse))
        res["groupcurve_boot"] = dict(deg_mid_CI=np.percentile(bm, [2.5, 50, 97.5]).tolist(),
                                      enh_mid_CI=np.percentile(be, [2.5, 50, 97.5]).tolist(),
                                      deg_start_CI=np.percentile(bs_, [2.5, 50, 97.5]).tolist(),
                                      enh_start_CI=np.percentile(bse, [2.5, 50, 97.5]).tolist(),
                                      frac_deg_mid_later=float((bm > be).mean()))

        # ---------------- (C) start of change: mice vs models (paired, noise-free) + posterior predictive
        cmp_, ps = {}, []
        for k in MAIN:
            sp = np.array(res["A_models"][k]["sig_pred"]["deg_start"])
            mp = np.array(res["A_models"][k]["sig_pred"]["deg_mid"])
            pr = paired(fd["start"], sp)
            cmp_[k] = dict(start_vs_model=pr, mid_vs_model=paired(fd["mid"], mp),
                           ppc_start_gm=ppc_p(res["A_models"][k]["sig_sim"]["deg_start_gm"], fd["start"].mean()),
                           ppc_mid_gm=ppc_p(res["A_models"][k]["sig_sim"]["deg_mid_gm"], fd["mid"].mean()),
                           ppc_groupcurve_start=ppc_p(res["A_models"][k]["sig_sim_groupcurve"]["deg_start"],
                                                      gd_["start"][0]),
                           ppc_groupcurve_mid=ppc_p(res["A_models"][k]["sig_sim_groupcurve"]["deg_mid"], gd_["mid"][0]),
                           sim_groupcurve_start_CI=np.percentile(res["A_models"][k]["sig_sim_groupcurve"]["deg_start"],
                                                                 [2.5, 50, 97.5]).tolist(),
                           sim_groupcurve_mid_CI=np.percentile(res["A_models"][k]["sig_sim_groupcurve"]["deg_mid"],
                                                               [2.5, 50, 97.5]).tolist(),
                           ppc_enh_groupcurve_mid=ppc_p(res["A_models"][k]["sig_sim_groupcurve"]["enh_mid"], ge_["mid"][0]),
                           sim_enh_groupcurve_mid_CI=np.percentile(res["A_models"][k]["sig_sim_groupcurve"]["enh_mid"],
                                                                   [2.5, 50, 97.5]).tolist())
            ps.append(pr["p_t"])
        for k, a in zip(MAIN, holm(ps)):
            cmp_[k]["start_vs_model"]["p_t_holm"] = a
        res["C_start_vs_models"] = cmp_

        # ---------------- reference days: within-session fall (second half minus first half)
        refd = {}
        for d in (3, 4, 5, 11, 12, 13):
            v = [REF[m][d][15:30].mean() - REF[m][d][:15].mean() for m in ms if d in REF[m] and len(REF[m][d]) == 30]
            refd[str(d)] = dict(mean=float(np.mean(v)), n=len(v)) if v else None
        base = [np.mean([REF[m][d][15:30].mean() - REF[m][d][:15].mean() for d in (4, 5) if d in REF[m]]) for m in ms]
        d13 = [REF[m][13][15:30].mean() - REF[m][13][:15].mean() for m in ms]
        d14 = [OBS[m][DI[14], 15:].mean() - OBS[m][DI[14], :15].mean() for m in ms]
        res["within_session_reference"] = dict(by_day=refd, day14_vs_baseline45=paired(d14, base),
                                               day13_vs_baseline45=paired(d13, base))
        # per-mouse rows
        for i, m in enumerate(ms):
            row = dict(mouse=m, group=g)
            for k in A_KEYS:
                row["obs_" + k] = float(ao[k][i])
            for k in MAIN:
                ap = a_stats(Pg[k][i])
                row["pred_day14_" + k] = float(ap["day14"])
                row["pred_night_halves_" + k] = float(ap["night_halves"])
                row["pred_deg_start_" + k] = res["A_models"][k]["sig_pred"]["deg_start"][i]
                row["pred_deg_mid_" + k] = res["A_models"][k]["sig_pred"]["deg_mid"][i]
            row.update(deg_A=float(fd["A"][i]), deg_mid=float(fd["mid"][i]), deg_s=float(fd["s"][i]),
                       deg_start=float(fd["start"][i]), deg_r2=float(fd["r2"][i]),
                       enh_A=float(fe["A"][i]), enh_mid=float(fe["mid"][i]), enh_s=float(fe["s"][i]),
                       enh_start=float(fe["start"][i]), enh_r2=float(fe["r2"][i]),
                       day14_water_drops=int(DROPS[m][DI[14]].sum()))
            rows_pm.append(row)
        # Active vs Passive is added after both groups
        S[g] = res
        log(g, "obs day14 / day15 / day16:", [round(res["A_observed"][k]["mean"], 3) for k in ("day14", "day15", "day16")])
        log(g, "obs halves: d14 h1", round(res["A_observed"]["d14_h1"]["mean"], 3), "h2", round(res["A_observed"]["d14_h2"]["mean"], 3),
            "d15 h1", round(res["A_observed"]["d15_h1"]["mean"], 3), "| within14", round(res["A_observed"]["within14_halves"]["mean"], 3),
            "p", round(res["A_observed"]["within14_halves"]["p"], 4), "| night", round(res["A_observed"]["night_halves"]["mean"], 3),
            "p", round(res["A_observed"]["night_halves"]["p"], 4))
        log(g, "per-mouse deg mid", np.round(fd["mid"], 1).tolist(), "A", np.round(fd["A"], 2).tolist(),
            "start", np.round(fd["start"], 1).tolist())
        log(g, "per-mouse enh mid", np.round(fe["mid"], 1).tolist(), "A", np.round(fe["A"], 2).tolist())
        log(g, "group curve deg", {k: round(v, 2) for k, v in res["sig_obs"]["groupcurve"]["deg"].items()},
            "enh", {k: round(v, 2) for k, v in res["sig_obs"]["groupcurve"]["enh"].items()})
        log(g, "deg vs enh mid (all):", {k: (round(v, 4) if isinstance(v, float) else v)
                                          for k, v in res["sig_tests"]["deg_vs_enh_mid_all"].items()})
        for k in MAIN:
            c = res["C_start_vs_models"][k]
            log(g, "start vs", k, "obs", round(c["start_vs_model"]["mean_a"], 1), "model", round(c["start_vs_model"]["mean_b"], 1),
                "p_t", round(c["start_vs_model"]["p_t"], 4), "holm", round(c["start_vs_model"]["p_t_holm"], 4),
                "| ppc groupcurve start", round(c["ppc_groupcurve_start"], 4), "mid", round(c["ppc_groupcurve_mid"], 4))
        log(g, "within-session reference:", res["within_session_reference"]["by_day"])

    # ---------------- Active vs Passive
    PM = pd.DataFrame(rows_pm)
    a, p = PM[PM.group == "Active"], PM[PM.group == "Passive"]
    S["Active_vs_Passive"] = {k: welch(a["obs_" + k], p["obs_" + k]) for k in A_KEYS}
    S["Active_vs_Passive"]["deg_mid"] = welch(a.deg_mid, p.deg_mid)
    S["Active_vs_Passive"]["enh_mid"] = welch(a.enh_mid, p.enh_mid)
    S["Active_vs_Passive"]["deg_start"] = welch(a.deg_start, p.deg_start)
    PM.to_csv(D / "m13_per_mouse.csv", index=False)
    json.dump(S, open(D / "m13_summary.json", "w"), indent=1)
    (D / "m13_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()
