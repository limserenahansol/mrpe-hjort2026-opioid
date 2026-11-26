"""Step 15 (18 Nov 2026): extra checks after the one independent review of the trial-by-trial slides (m13 + m14 = v1).

The review found (all numbers of m13 reproduce; the problems are inferential):
  1. m13's model check is a plug-in posterior predictive check that draws every 30-s bin independently from the fitted
     p. Real licking varies about twice as much from day to day and is autocorrelated within a session, so those p
     values are TOO SMALL (anti-conservative), not 'conservative' (m13's docstring was wrong). A placebo day pair with
     no outcome change (Passive day 12 -> 13) is rejected by every model under the naive check.
  2. 'No extra change overnight' and 'the start rules out instant switch' are not robust (wide CIs, n = 6).
  3. The midpoint comparison uses unequal windows (90 vs 60 trials).
This script adds, next to m13 (m13 outputs unchanged):
  (a) a CALIBRATED check: lick yes/no simulated with a mouse-day logit offset (sd sig_u) and AR(1) logit noise within
      each session (stationary sd sig_a, coefficient phi); the mean is kept at the model p (logistic-normal
      approximation, logit scaled by sqrt(1 + pi*var/8)). sig_u, sig_a and phi are fitted, per model and group, on the
      stable days 11-13 to three features of the DAY-TO-DAY DIFFERENCES of residuals at the same bin (11->12, 12->13),
      because every slide statistic compares the same bin across days (a within-session misfit that repeats every day
      cancels in these differences; a first attempt that fitted raw residuals could not match both features):
      T1 = dispersion of whole-session mean differences, T2 = dispersion of half-session mean differences (both as a
      ratio to independent bins), T3 = lag-1 autocorrelation of the bin differences. 2000 data sets per model.
      No single setting matches all three (the mice's bin-to-bin correlation is strong but their session totals
      vary little), so two calibrations are reported as a range: A fits T1+T2+T3 (overshoots T1 about 2x, so its
      whole-session p values are too large) and B fits T1+T2 only (undershoots T3, so its timing p values may be
      too small). The mouse-level paired t (c) needs no noise model and is the primary test on the slides.
  (b) placebo: day 13 minus day 12 (no outcome change), naive and calibrated check, and mouse-level t.
  (c) mouse-level paired t of observed minus predicted (n = 6 / 7), the most robust test (no noise model needed).
  (d) equal windows: drug -> water midpoint on days 14-15 only (60 trials, like water -> drug).
  (e) leave-one-mouse-out group-curve start and midpoint; drop front-loading on day 13; ceiling in the first 2 min.
Output: data/m15_summary.json, data/m15_log.txt
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
from scipy import stats
from scipy.special import expit, logit

warnings.filterwarnings("ignore")
sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m13_trial_by_trial as m13  # noqa: E402  (functions only; m13 outputs are not rewritten)

D = HERE / "data"
NSIM, NCAL, SEED = 2000, 200, 20261119
SU = np.round(np.arange(0.0, 2.01, 0.1), 2)
SA = np.round(np.arange(0.0, 3.01, 0.1), 2)
PHIS = (0.0, 0.3, 0.5, 0.7, 0.9)
MODELS, MAIN = m13.MODELS, m13.MAIN
DI = m13.DI
LOG = []


def log(*a):
    s = " ".join(str(x) for x in a)
    LOG.append(s)
    print(s)


# ============================================================================ noise model
def latent_draws(rng, shape_days, nsim):
    """common random numbers: offsets, AR innovations, uniforms for (nsim, n, days, 30)"""
    n, nd = shape_days
    return (rng.normal(size=(nsim, n, nd, 1)), rng.normal(size=(nsim, n, nd, 30)), rng.random((nsim, n, nd, 30)))


def simulate(P, su, sa, draws, phi):
    """P: (n, nd, 30) model p -> (nsim, n, nd, 30) yes/no with mouse-day offset and AR(1) noise, mean kept at P"""
    zu, e, U = draws
    var = su ** 2 + sa ** 2
    k = np.sqrt(1 + np.pi * var / 8)
    L = k * logit(np.clip(P, 1e-6, 1 - 1e-6))[None]
    a = np.empty_like(e)
    a[..., 0] = sa * e[..., 0]
    sin = sa * np.sqrt(1 - phi ** 2)
    for t in range(1, 30):
        a[..., t] = phi * a[..., t - 1] + sin * e[..., t]
    return (U < expit(L + su * zu + a)).astype(float)


def disp_targets(C, P):
    """features of day-to-day residual differences at the same bin on the stable days (C: (..., n, 3, 30), P: (n, 3, 30)):
    T1 whole-session mean difference, T2 half-session mean difference (both: mean square / independent-bin variance),
    T3 pooled lag-1 autocorrelation of the bin differences"""
    r = C - P
    dr = r[..., 1:, :] - r[..., :-1, :]
    v = P * (1 - P)
    vd = v[:, 1:, :] + v[:, :-1, :]
    T1 = (dr.mean(-1) ** 2).mean((-1, -2)) / (vd.sum(-1) / 900.0).mean()
    hm = np.stack([dr[..., :15].mean(-1), dr[..., 15:].mean(-1)], -1)
    hv = np.stack([vd[..., :15].sum(-1) / 225.0, vd[..., 15:].sum(-1) / 225.0], -1)
    T2 = (hm ** 2).mean((-1, -2, -3)) / hv.mean()
    num = (dr[..., 1:] * dr[..., :-1]).sum((-1, -2, -3))
    den = np.sqrt((dr[..., 1:] ** 2).sum((-1, -2, -3)) * (dr[..., :-1] ** 2).sum((-1, -2, -3)))
    return T1, T2, num / den


def calibrate(O, P, rng):
    """fit (sig_u, sig_a, phi) on days 11-13; A = match T1, T2, T3; B = match T1, T2 only"""
    Os, Ps = O[:, 0:3], P[:, 0:3]
    t1o, t2o, t3o = disp_targets(Os, Ps)
    dr = latent_draws(rng, Ps.shape[:2], NCAL)
    rows = []
    for phi in PHIS:
        for su in SU:
            for sa in SA:
                if phi > 0 and sa == 0:
                    continue
                C = simulate(Ps, su, sa, dr, phi)
                t1, t2, t3 = (float(np.mean(x)) for x in disp_targets(C, Ps))
                rows.append((su, sa, phi, t1, t2, t3))
    R = np.array(rows)
    fA = np.log(R[:, 3] / t1o) ** 2 + np.log(R[:, 4] / t2o) ** 2 + ((R[:, 5] - t3o) / 0.05) ** 2
    fB = np.log(R[:, 3] / t1o) ** 2 + np.log(R[:, 4] / t2o) ** 2
    ind = R[(R[:, 0] == 0) & (R[:, 1] == 0)][0]
    out = {}
    for name, f in (("A", fA), ("B", fB)):
        b = R[int(np.argmin(f))]
        out[name] = dict(sig_u=float(b[0]), sig_a=float(b[1]), phi=float(b[2]), T1_sim=float(b[3]), T2_sim=float(b[4]),
                         T3_sim=float(b[5]), loss=float(f.min()),
                         at_grid_edge=bool(b[0] >= SU.max() or b[1] >= SA.max()))
    out.update(T1_obs=float(t1o), T2_obs=float(t2o), T3_obs=float(t3o), T1_indep=float(ind[3]), T2_indep=float(ind[4]),
               T3_indep=float(ind[5]))
    return out


# ============================================================================ statistics
def placebo(C):
    """day 13 minus day 12 (no outcome change), whole session and halves"""
    d = C[..., DI[13], :] - C[..., DI[12], :]
    return dict(whole=d.mean(-1), h1=d[..., :15].mean(-1), h2=d[..., 15:].mean(-1))


def ppc(sim_gm, obs):
    return dict(p=m13.ppc_p(sim_gm, obs), lo=float(np.percentile(sim_gm, 2.5)), hi=float(np.percentile(sim_gm, 97.5)),
                centre=float(np.mean(sim_gm)))


def timing_stats(C, grids):
    """per-mouse and group-curve sigmoid timing for (..., n, 8, 30)"""
    GD, GE, GD60 = grids
    td, te = m13.traj_deg(C), m13.traj_enh(C)
    n = C.shape[-3]
    lead = C.shape[:-3]
    fd = GD.fit(td.reshape(-1, 90))
    fe = GE.fit(te.reshape(-1, 60))
    gd = GD.fit(td.mean(-2).reshape(-1, 90))
    ge = GE.fit(te.mean(-2).reshape(-1, 60))
    g60 = GD60.fit(td[..., :60].mean(-2).reshape(-1, 60))
    rs = lambda v: v.reshape(lead + (n,)) if lead else v
    rg = lambda v: v.reshape(lead) if lead else v[0]
    return dict(deg_start_pm=rs(fd["start"]), deg_mid_pm=rs(fd["mid"]), enh_mid_pm=rs(fe["mid"]),
                g_deg_start=rg(gd["start"]), g_deg_mid=rg(gd["mid"]), g_enh_mid=rg(ge["mid"]), g_enh_start=rg(ge["start"]),
                g_deg60_mid=rg(g60["mid"]))


# ============================================================================ main
def main():
    rng = np.random.default_rng(SEED)
    OBS, PRED, DROPS, REF, GROUP = m13.build()
    grids = (m13.GD, m13.GE, m13.Grid(60))
    S = {"meta": dict(nsim=NSIM, ncal=NCAL, seed=SEED, phis=list(PHIS), grid_sig_u=SU.tolist(), grid_sig_a=SA.tolist(),
                      note="calibrated = mouse-day logit offset + AR(1) logit noise, fitted per model on day-to-day residual differences of days 11-13")}
    for g in ("Active", "Passive"):
        ms = [m for m in OBS if GROUP[m] == g]
        O = np.stack([OBS[m] for m in ms])
        res = {"mice": ms}
        # ---------------- observed timing incl. equal windows and leave-one-out
        to = timing_stats(O, grids)
        res["obs_timing"] = {k: (v.tolist() if isinstance(v, np.ndarray) else float(v)) for k, v in to.items()}
        f60 = grids[2].fit(m13.traj_deg(O)[:, :60])
        fe = grids[1].fit(m13.traj_enh(O))
        res["equal_window"] = dict(deg60_mid_pm=f60["mid"].tolist(), deg60_A_pm=f60["A"].tolist(),
                                   paired_deg60_vs_enh=m13.paired(f60["mid"], fe["mid"]),
                                   group_deg60_mid=float(to["g_deg60_mid"]), group_enh_mid=float(to["g_enh_mid"]))
        loo = []
        for i in range(len(ms)):
            keep = [j for j in range(len(ms)) if j != i]
            t_ = timing_stats(O[keep], grids)
            loo.append(dict(left_out=ms[i], deg_start=float(t_["g_deg_start"]), deg_mid=float(t_["g_deg_mid"]),
                            enh_mid=float(t_["g_enh_mid"])))
        res["leave_one_out"] = loo
        # ---------------- drops front-loaded on every PR day; ceiling in the first 2 min
        res["drops_frac_by_2p5min"] = {str(d): float(np.mean([DROPS[m][DI[d], :5].sum() / max(DROPS[m][DI[d]].sum(), 1)
                                                               for m in ms])) for d in (13, 14, 15, 17)}
        res["p_lick_bins0_3"] = {str(d): float(O[:, DI[d], :4].mean()) for d in (13, 14, 16, 17)}
        # ---------------- per model: calibration, calibrated + naive checks, mouse-level tests
        ao = m13.a_stats(O)
        po = placebo(O)
        res["models"] = {}
        for k in MODELS:
            P = np.stack([PRED[k][m] for m in ms])
            cal = calibrate(O, P, rng)
            dr = latent_draws(rng, P.shape[:2], NSIM)
            Cc = simulate(P, cal["A"]["sig_u"], cal["A"]["sig_a"], dr, cal["A"]["phi"])   # calibrated A
            Cb = simulate(P, cal["B"]["sig_u"], cal["B"]["sig_a"], dr, cal["B"]["phi"])   # calibrated B
            Cn = (rng.random((NSIM,) + P.shape) < P[None]).astype(float)                  # naive (as m13)
            ap, ac, ab, an = m13.a_stats(P), m13.a_stats(Cc), m13.a_stats(Cb), m13.a_stats(Cn)
            pp, pc_, pb, pn = placebo(P), placebo(Cc), placebo(Cb), placebo(Cn)
            mk = dict(calibration=cal, A={}, placebo={})
            for key in ("day14", "day15", "day16", "d14_h1", "d14_h2", "d15_h1", "within14_halves", "night_halves"):
                mk["A"][key] = dict(obs=float(ao[key].mean()), pred=float(ap[key].mean()),
                                    calibrated=ppc(ac[key].mean(1), ao[key].mean()),
                                    calibrated_B=ppc(ab[key].mean(1), ao[key].mean()),
                                    naive=ppc(an[key].mean(1), ao[key].mean()),
                                    mouse_t=m13.paired(ao[key], ap[key]))
            for key in ("whole", "h1", "h2"):
                mk["placebo"][key] = dict(obs=float(po[key].mean()), pred=float(pp[key].mean()),
                                          calibrated=ppc(pc_[key].mean(1), po[key].mean()),
                                          calibrated_B=ppc(pb[key].mean(1), po[key].mean()),
                                          naive=ppc(pn[key].mean(1), po[key].mean()),
                                          mouse_t=m13.paired(po[key], pp[key]))
            prof_c = m13.profile25(Cc).mean(1)
            mk["profile25_band_calibrated"] = dict(lo=np.percentile(prof_c, 2.5, 0).tolist(),
                                                   hi=np.percentile(prof_c, 97.5, 0).tolist())
            if k in MAIN:
                tc = timing_stats(Cc, grids)
                tb = timing_stats(Cb, grids)
                tn = timing_stats(Cn, grids)
                mk["timing"] = {}
                for key, obs in (("g_deg_start", to["g_deg_start"]), ("g_deg_mid", to["g_deg_mid"]),
                                 ("g_enh_mid", to["g_enh_mid"]), ("g_deg60_mid", to["g_deg60_mid"])):
                    mk["timing"][key] = dict(obs=float(obs), calibrated=ppc(tc[key], obs), calibrated_B=ppc(tb[key], obs),
                                             naive=ppc(tn[key], obs), cal_median=float(np.median(tc[key])),
                                             calB_median=float(np.median(tb[key])))
                for key, obs in (("deg_start_pm", to["deg_start_pm"]), ("deg_mid_pm", to["deg_mid_pm"])):
                    mk["timing"][key + "_mean"] = dict(obs=float(np.mean(obs)), calibrated=ppc(tc[key].mean(-1), np.mean(obs)),
                                                       calibrated_B=ppc(tb[key].mean(-1), np.mean(obs)),
                                                       naive=ppc(tn[key].mean(-1), np.mean(obs)))
            res["models"][k] = mk
            a = mk["A"]["day14"]
            ca, cb = cal["A"], cal["B"]
            log(g, f"{k:16s}", f"A u {ca['sig_u']:.1f} a {ca['sig_a']:.1f} phi {ca['phi']:.1f} T1 {cal['T1_obs']:.2f}/{ca['T1_sim']:.2f} "
                f"T3 {cal['T3_obs']:.2f}/{ca['T3_sim']:.2f} | B u {cb['sig_u']:.1f} a {cb['sig_a']:.1f} phi {cb['phi']:.1f} "
                f"T1 {cb['T1_sim']:.2f} T2 {cal['T2_obs']:.2f}/{cb['T2_sim']:.2f} T3 {cb['T3_sim']:.2f} | day14 obs {a['obs']:+.3f} "
                f"pred {a['pred']:+.3f} p naive {a['naive']['p']:.4f} A {a['calibrated']['p']:.4f} B {a['calibrated_B']['p']:.4f} "
                f"mouse-t {a['mouse_t']['p_t']:.3f} | placebo naive {mk['placebo']['whole']['naive']['p']:.4f} "
                f"A {mk['placebo']['whole']['calibrated']['p']:.4f} B {mk['placebo']['whole']['calibrated_B']['p']:.4f}"
                + (f" | start A {mk['timing']['g_deg_start']['calibrated']['p']:.3f} B {mk['timing']['g_deg_start']['calibrated_B']['p']:.3f}"
                   f" mid A {mk['timing']['g_deg_mid']['calibrated']['p']:.3f} B {mk['timing']['g_deg_mid']['calibrated_B']['p']:.3f}"
                   if k in MAIN else ""))
        # Holm across the 4 main models (calibrated day-14 check and mouse-level t)
        for src in ("calibrated", "calibrated_B", "mouse_t"):
            ps = [res["models"][k]["A"]["day14"][src]["p" if src != "mouse_t" else "p_t"] for k in MAIN]
            for k, a in zip(MAIN, m13.holm(ps)):
                res["models"][k]["A"]["day14"][src]["p_holm"] = a
        S[g] = res
        log(g, "equal windows: group deg60 mid", round(res["equal_window"]["group_deg60_mid"], 1), "enh mid",
            round(res["equal_window"]["group_enh_mid"], 1), "| per mouse deg60",
            np.round(res["equal_window"]["deg60_mid_pm"], 1).tolist(), "paired p",
            round(res["equal_window"]["paired_deg60_vs_enh"]["p_t"], 4))
        log(g, "LOO group start", [round(x["deg_start"], 1) for x in loo], "mid", [round(x["deg_mid"], 1) for x in loo],
            "enh", [round(x["enh_mid"], 1) for x in loo])
        log(g, "drops by 2.5 min", {d: round(v, 2) for d, v in res["drops_frac_by_2p5min"].items()},
            "P(lick) bins 0-3", {d: round(v, 2) for d, v in res["p_lick_bins0_3"].items()})
    json.dump(S, open(D / "m15_summary.json", "w"), indent=1)
    (D / "m15_log.txt").write_text("\n".join(LOG), encoding="utf-8")


if __name__ == "__main__":
    main()
