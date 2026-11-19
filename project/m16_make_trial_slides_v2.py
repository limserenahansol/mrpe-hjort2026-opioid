"""Step 16 (18 Nov 2026): slides 7-8 v2 after the one independent review (v1 = m14 outputs, kept unchanged).

What changed against v1 (review findings; numbers from m13 + m15):
- slide 7: title states only the robust observation (no clear fall in the first 7.5 min); 'no extra change overnight'
  removed (within day 14 vs overnight cannot be separated: CIs about +-0.5, mirror images across mice); panel c now
  shows mice minus each model per mouse with the mouse-level paired t (no noise model needed); the independent-bin
  simulation check is dropped (too narrow: it rejected a no-change placebo day pair); noise-matched simulation p values
  (two calibrations, m15) are given as a range in the footer and notes; drops are 'earned', front-loaded on every PR
  day (day 13 shown); the first 2 min are at ceiling.
- slide 8: title says 'group curve'; per-mouse p = 0.11 and equal 60-trial windows printed; the start 'rules out' claim
  removed (start 95% CI 0-47; leave-one-out 6-42; noise-matched checks); curves described as raw change (v1 notes
  wrongly said normalised); per-mouse fits are steps; water -> drug midpoint is an upper bound and predicted by every
  model; the paper's Fig. 1m used model value curves.
Outputs (slides_20261118_trial_by_trial\\): mrpe_trial_by_trial_2slides_{EN,KR}_v2_20261118.pptx, single_slides_v2\\,
  mrpe_hjort2026_opioid_8slides_{EN,KR}_v4_20261118.pptx (v2 slides 1-6 byte-identical + slides 7-8 v2),
  preview_v2_{EN,KR}\\, preview_v4_{EN,KR}\\, figs_v2\\, fit_report_v2.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m10_make_slides as b  # noqa: E402
import m14_make_trial_slides as v1  # noqa: E402  (data loaders, traj, restore_v2_parts)

plt, Line2D, Presentation = b.plt, b.Line2D, b.Presentation
P, Rn, tbox, measure, at, place, key_box, common = b.P, b.Rn, b.tbox, b.measure, b.at, b.place, b.key_box, b.common
frac, sf, fp, peq, reg, check, T = b.frac, b.sf, b.fp, b.peq, b.reg, b.check, b.T
D = b.D
OUT = v1.OUT
FIGD = OUT / "figs_v2"
for p_ in (FIGD, OUT / "single_slides_v2"):
    p_.mkdir(parents=True, exist_ok=True)
V2 = v1.V2
RED, BLUE, GREY, GREEN, ORANGE, PURPLE = b.RED, b.BLUE, b.GREY, b.GREEN, b.ORANGE, b.PURPLE
S = v1.S
S15 = json.load(open(D / "m15_summary.json"))
PM, BINS = v1.PM, v1.BINS
MAIN, LAB, CM = v1.MAIN, v1.LAB, v1.CM
LAB_KR = {"today's outcome": "즉시 전환", "RPE": "RPE", "RPE2a": "RPE 2α", "mRPE": "meta-RPE"}
NA, NP = v1.NA, v1.NP
NSIM = S15["meta"]["nsim"]


def save(fig, name):
    p_ = FIGD / name
    fig.savefig(p_, facecolor="white")
    plt.close(fig)
    return p_


def fp2(p):
    """p for a simulation check: smallest possible value is 1/(NSIM+1)"""
    return "< 0.001" if p <= 1.0 / (NSIM + 1) + 1e-12 else fp(p)


def diff_rows(g="Active"):
    """per-mouse mice-minus-model day-14 change, mean, 95% t CI, paired p (= one-sample t of the differences)"""
    x = PM[PM.group == g]
    out = {}
    for k in MAIN:
        d = (x.obs_day14 - x["pred_day14_" + k]).values
        se = d.std(ddof=1) / np.sqrt(len(d))
        h = stats.t.ppf(0.975, len(d) - 1) * se
        out[k] = dict(d=d, mean=float(d.mean()), lo=float(d.mean() - h), hi=float(d.mean() + h),
                      p=float(stats.ttest_1samp(d, 0).pvalue))
        check(f"mouse-level p for {k} equals m15 mouse_t",
              abs(out[k]["p"] - S15[g]["models"][k]["A"]["day14"]["mouse_t"]["p_t"]) < 1e-9)
    return out


DR = diff_rows()


def calp(k, key, src, g="Active", part="A"):
    return S15[g]["models"][k][part][key][src]["p"]


# ============================================================================================ figure 7 v2
def fig7():
    A = S["Active"]
    fig = plt.figure(figsize=(12.7, 4.85))
    # ---- a: drops earned through the session: day 14 (water) both groups, day 13 (morphine) Active
    a = fig.add_axes([0.062, 0.17, 0.175, 0.70])
    for g, d, col, ls, lab in (("Active", 14, RED, "-", "Active d14 (water)"),
                               ("Passive", 14, BLUE, "-", "Passive d14 (water)"),
                               ("Active", 13, RED, "--", "Active d13 (morphine)")):
        x = BINS[(BINS.group == g) & (BINS.day == d)].pivot_table(index="mouse", columns="bin", values="drops_in_bin")
        cum = np.c_[np.zeros(len(x)), np.cumsum(x.values, 1)]
        tt = np.arange(31) * 0.5
        mu, se = cum.mean(0), cum.std(0, ddof=1) / np.sqrt(len(cum))
        a.plot(tt, mu, color=col, lw=2.2 if ls == "-" else 1.6, ls=ls, label=lab)
        if ls == "-":
            a.fill_between(tt, mu - se, mu + se, color=col, alpha=0.16, lw=0)
    a.axvline(7.5, color="#999999", lw=1, ls=":")
    a.set_xlim(0, 15)
    a.set_ylim(0, 36)
    a.set_xticks([0, 2.5, 7.5, 15], ["0", "2.5", "7.5", "15"])
    a.set_xlabel("Minute of the session")
    a.set_ylabel("Drops earned so far\n(mean ± s.e.m.)")
    a.set_title(f"a  Drops earned (n = {NA} / {NP})", loc="left", fontweight="bold")
    a.legend(loc="lower right", frameon=False, fontsize=9.4, handlelength=1.6, borderaxespad=0.2)
    # ---- b: 2.5-min profile, Active, observed vs models (as v1) + ceiling of the first 2 min
    bx = fig.add_axes([0.315, 0.17, 0.405, 0.70])
    ob = np.array(A["profile25_observed"]["mean"])
    se = np.array(A["profile25_observed"]["sem"])
    off = {14: 0.0, 15: 17.0, 16: 34.0}
    xs = np.concatenate([off[d] + np.arange(6) * 2.5 + 1.25 for d in (14, 15, 16)])
    bx.axvspan(0, 2.0, color="#F3E6C8", lw=0, zorder=0)
    bx.text(1.0, -1.21, "ceiling", ha="center", va="bottom", fontsize=9, color="#8A6D1F", rotation=90)
    for k in MAIN:
        pr = np.array(A["profile25_models"][k])
        for j in range(3):
            sl = slice(6 * j, 6 * j + 6)
            bx.plot(xs[sl], pr[sl], color=CM[k], lw=2.0 if k in ("RPE2a", "mRPE") else 1.6,
                    ls="--" if k == "mRPE" else "-", zorder=2)
    for j, d in enumerate((14, 15, 16)):
        sl = slice(6 * j, 6 * j + 6)
        bx.errorbar(xs[sl], ob[sl], yerr=se[sl], color="black", marker="o", ms=5.5, lw=1.6, capsize=2.5, zorder=4)
        if d < 16:
            bx.axvspan(off[d] + 15.0, off[d] + 17.0, color="#EEEEEE", lw=0)
            bx.text(off[d] + 16.0, -1.21, "night", ha="center", va="bottom", fontsize=9.5, color="#666666", rotation=90)
    bx.axhline(0, color="#777777", lw=1)
    bx.axvline(7.5, color="#BBBBBB", lw=1, ls=":")
    bx.set_xlim(-0.5, 49.5)
    bx.set_ylim(-1.25, 0.58)
    bx.set_yticks([-0.75, -0.5, -0.25, 0, 0.25, 0.5])
    bx.set_xticks([0, 7.5, 15, 17, 24.5, 32, 34, 41.5, 49],
                  ["0", "7.5\nday 14", "15", "0", "7.5\nday 15", "15", "0", "7.5\nday 16", "15"])
    bx.set_xlabel("Minute of the session (2.5-min means)", labelpad=2)
    bx.set_ylabel("Licking minus the same minute\non day 13 (Δ P(lick), Active)")
    bx.set_title(f"b  Licking vs day 13, Active (n = {NA}): mice and models", loc="left", fontweight="bold")
    hs = [Line2D([], [], color="black", marker="o", lw=1.6, label="mice ± s.e.m.")]
    hs += [Line2D([], [], color=CM[k], lw=2, ls="--" if k == "mRPE" else "-", label=LAB[k]) for k in MAIN]
    bx.legend(handles=hs, loc="lower left", bbox_to_anchor=(0.045, 0.0), frameon=True, fontsize=9.6, ncol=1,
              borderpad=0.35, framealpha=0.95, edgecolor="#DDDDDD", handlelength=1.6, borderaxespad=0.3)
    # ---- c: whole day 14, mice minus each model, per mouse
    c = fig.add_axes([0.80, 0.17, 0.185, 0.70])
    for i, k in enumerate(MAIN):
        r = DR[k]
        jit = np.linspace(-0.17, 0.17, len(r["d"]))
        c.scatter(r["d"], np.full(len(r["d"]), i) + jit, s=20, color=CM[k], alpha=0.7, lw=0, zorder=3)
        c.errorbar(r["mean"], i, xerr=[[r["mean"] - r["lo"]], [r["hi"] - r["mean"]]], color="black", marker="D",
                   ms=6.5, capsize=4, lw=2, zorder=4)
        c.text(0.86, i, peq(r["p"]), ha="right", va="center", fontsize=10.5, color=CM[k])
    lo = min(min(r["d"].min(), r["lo"]) for r in DR.values())
    hi = max(max(r["d"].max(), r["hi"]) for r in DR.values())
    check("fig7c: all mice and CIs inside the axis", lo > -0.45 and hi < 0.62)
    c.axvline(0, color="#777777", lw=1.2)
    c.set_yticks(range(4), ["instant\nswitch", "RPE", "RPE 2α", "meta-RPE"], fontsize=11)
    c.set_ylim(3.6, -0.6)
    c.set_xlim(-0.45, 0.88)
    c.set_xticks([-0.25, 0, 0.25, 0.5])
    c.text(0.02, -0.45, "mice changed less →", fontsize=9.5, color="#555555", va="center")
    c.set_xlabel("Mice minus model, whole day 14\n(Δ P(lick); dots = mice)")
    c.set_title("c  Day 14: mice − model", loc="left", fontweight="bold")
    tg = {"key_tl": frac(fig, bx, 18.4, 0.55), "key_br": frac(fig, bx, 49.0, 0.06),
          "tip": frac(fig, bx, 4.6, ob[1] + se[1] + 0.03)}
    return save(fig, "s7_day14_v2.png"), fig.get_size_inches(), tg


# ============================================================================================ figure 8 v2
def fig8():
    A = S["Active"]
    fig = plt.figure(figsize=(12.7, 4.85))
    a = fig.add_axes([0.065, 0.17, 0.33, 0.70])
    gc = A["sig_obs"]["groupcurve"]
    YL = (-0.90, 0.78)
    for kind, col, mk, lab in (("deg", "#222222", "o", "drug → water"), ("enh", "#8A8A8A", "^", "water → drug")):
        Tm = v1.traj("Active", kind)
        n = Tm.shape[1]
        blk = Tm.reshape(len(Tm), n // 5, 5).mean(2)
        mu, se = blk.mean(0), blk.std(0, ddof=1) / np.sqrt(len(blk))
        xc = np.arange(n // 5) * 5 + 2.0
        check(f"fig8a {kind}: block means inside the axis", (mu - se).min() > YL[0] + 0.1 and (mu + se).max() < YL[1])
        a.errorbar(xc, mu, yerr=se, color=col, marker=mk, ms=5.5, lw=0, elinewidth=1.3, capsize=2.5, zorder=3, label=lab)
        Am, x0, s_ = gc[kind]["A"], gc[kind]["x0"], gc[kind]["s"]
        xf = np.linspace(0, n - 1, 400)
        a.plot(xf, Am * expit((xf - x0) / s_), color=col, lw=2.8, zorder=2)
        a.axvline(x0, color=col, lw=1.2, ls="--")
        a.text(x0 + 1.2, -0.82, f"midpoint {x0:.1f}", color=col, fontsize=10.5, fontweight="bold", va="center", ha="left")
    for xv in (30, 60):
        a.axvline(xv, color="#BBBBBB", lw=1, ls=":")
        a.text(xv - 0.6, 0.76, "next day", fontsize=9.5, color="#777777", va="top", ha="right", rotation=90)
    a.axhline(0, color="#777777", lw=1)
    a.set_xlim(-1, 90)
    a.set_ylim(*YL)
    a.set_xlabel("Trial (30-s bin) from the switch")
    a.set_ylabel("Licking minus the day before the switch\n(Δ P(lick), mean ± s.e.m. of 2.5 min)")
    a.set_title(f"a  Licking by trial, Active (n = {NA}); lines = fits", loc="left", fontweight="bold")
    a.legend(loc="upper right", frameon=True, fontsize=9.8, framealpha=0.95, edgecolor="#DDDDDD", borderaxespad=0.3)
    # ---- b: per-mouse midpoints (as v1)
    bx = fig.add_axes([0.47, 0.17, 0.23, 0.70])
    for j, (g, col) in enumerate((("Active", RED), ("Passive", BLUE))):
        x = PM[PM.group == g]
        x0s = j * 2.6
        for _, r in x.iterrows():
            bx.plot([x0s, x0s + 1], [r.deg_mid, r.enh_mid], color=col, alpha=0.55, lw=1.2, marker="o", ms=5)
        bx.plot([x0s - 0.22, x0s + 0.22], [x.deg_mid.mean()] * 2, color="black", lw=2.6)
        bx.plot([x0s + 0.78, x0s + 1.22], [x.enh_mid.mean()] * 2, color="black", lw=2.6)
        t = S[g]["sig_tests"]["deg_vs_enh_mid_all"]
        bx.text(x0s + 0.5, 96, f"{g} (n = {len(x)})\npaired t {peq(t['p_t'])}", ha="center", va="bottom", fontsize=10.5,
                color=col, fontweight="bold")
    bx.set_xticks([0, 1, 2.6, 3.6], ["drug→\nwater", "water→\ndrug", "drug→\nwater", "water→\ndrug"], fontsize=10.5)
    for yv in (30, 60):
        bx.axhline(yv, color="#BBBBBB", lw=1, ls=":")
    bx.set_xlim(-0.55, 4.15)
    bx.set_ylim(-3, 116)
    bx.set_yticks([0, 30, 60, 90])
    bx.set_ylabel("Midpoint (trial from the switch)")
    bx.set_title("b  Midpoint per mouse (bar = mean)", loc="left", fontweight="bold")
    # ---- c: start, mice vs noise-matched simulations (calibration A bars; p for A and B)
    c = fig.add_axes([0.80, 0.17, 0.185, 0.70])
    obs = A["sig_obs"]["groupcurve"]["deg"]["start"]
    lo, med, hi = A["groupcurve_boot"]["deg_start_CI"]
    c.axvspan(lo, hi, color="#DDDDDD", lw=0, zorder=0)
    c.axvline(obs, color="black", lw=2, zorder=3)
    for i, k in enumerate(MAIN):
        tm = S15["Active"]["models"][k]["timing"]["g_deg_start"]
        q0, q2, qm = tm["calibrated"]["lo"], tm["calibrated"]["hi"], tm["cal_median"]
        c.plot([q0, q2], [i, i], color=CM[k], lw=5, solid_capstyle="butt", alpha=0.85)
        c.plot(qm, i, marker="o", color="white", mec=CM[k], mew=2, ms=8)
        c.text(89, i - 0.22, f"p = {fp(tm['calibrated']['p'])}", ha="right", va="center", fontsize=10, color=CM[k])
        c.text(89, i + 0.2, f"({fp(tm['calibrated_B']['p'])})", ha="right", va="center", fontsize=9.5, color=CM[k])
    c.set_yticks(range(4), ["instant\nswitch", "RPE", "RPE 2α", "meta-RPE"], fontsize=11)
    c.set_ylim(3.6, -0.75)
    c.set_xlim(-2, 90)
    c.set_xticks([0, 30, 60, 90])
    c.text(obs + 1.5, -0.45, f"mice {obs:.1f}", fontsize=10.5, fontweight="bold")
    c.text(hi - 1.0, 3.45, "mice 95% CI", fontsize=9, color="#666666", ha="right", va="bottom")
    c.set_xlabel("Start of the change (trial)\ndrug → water, group curve")
    c.set_title("c  Start: mice vs models", loc="left", fontweight="bold")
    return save(fig, "s8_switch_timing_v2.png"), fig.get_size_inches(), None


# ============================================================================================ slide 7 v2
def slide7(prs, F):
    b.SLIDE = 7
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    A, Pv = S["Active"], S["Passive"]
    o = A["A_observed"]
    f75 = reg("w7_f75", f"{A['day14_drops']['frac_by_2p5min'][2] * 100:.0f}", A["day14_drops"]["frac_by_2p5min"][2] * 100)
    h1 = reg("w7_h1", sf(o["d14_h1"]["mean"]), o["d14_h1"]["mean"])
    d14 = reg("w7_d14", sf(o["day14"]["mean"]), o["day14"]["mean"])
    means = [DR[k]["mean"] for k in MAIN]
    ps = [DR[k]["p"] for k in MAIN]
    check("w7: every model predicted a larger fall than the mice (mice - model > 0)", min(means) > 0)
    check("w7: no model rejected at the mouse level (all p >= 0.05)", min(ps) >= 0.05)
    dmin, dmax = reg("w7_dmin", f"{min(means):.2f}", min(means)), reg("w7_dmax", f"{max(means):.2f}", max(means))
    pmin, pmax = reg("w7_pmin", fp(min(ps)), min(ps)), reg("w7_pmax", fp(max(ps)), max(ps))
    wi, wlo, whi = o["within14_halves"]["mean"], o["within14_halves"]["lo"], o["within14_halves"]["hi"]
    ng, nlo, nhi = o["night_halves"]["mean"], o["night_halves"]["lo"], o["night_halves"]["hi"]
    w_, wl_, wh_ = reg("w7_w", sf(wi), wi), reg("w7_wl", sf(wlo), wlo), reg("w7_wh", sf(whi), whi)
    n_, nl_, nh_ = reg("w7_n", sf(ng), ng), reg("w7_nl", sf(nlo), nlo), reg("w7_nh", sf(nhi), nhi)
    sim = {k: (calp(k, "day14", "calibrated"), calp(k, "day14", "calibrated_B")) for k in MAIN}
    sim_txt_en = ", ".join(f"{LAB[k]} {fp2(min(v))}–{fp2(max(v))}" for k, v in sim.items())
    sim_txt_kr = ", ".join(f"{LAB_KR[k]} {fp2(min(v))}–{fp2(max(v))}" for k, v in sim.items())
    for k, v in sim.items():
        for tag, val in (("lo", min(v)), ("hi", max(v))):
            if val > 1 / (NSIM + 1) + 1e-12:
                reg(f"w7_sim{tag}_{k}", fp2(val), val)
    pi_m = reg("w7_pim", fp(DR[MAIN[0]]["p"]), DR[MAIN[0]]["p"])
    pi_lo = reg("w7_pilo", fp2(min(sim[MAIN[0]])), min(sim[MAIN[0]]))
    pi_hi = reg("w7_pihi", fp2(max(sim[MAIN[0]])), max(sim[MAIN[0]]))
    check("w7: instant switch has the smallest mouse-level and simulation p",
          DR[MAIN[0]]["p"] == min(ps) and max(sim[MAIN[0]]) <= min(max(v) for v in sim.values()))
    check("w7: every other main model has simulation p >= 0.05", all(min(sim[k]) >= 0.05 for k in MAIN[1:]))
    plac = S15["Passive"]["models"]["RPE2a"]["placebo"]["whole"]["naive"]["p"]
    check("w7: naive check rejects the Passive placebo (motivates the change)", plac < 0.01)
    pv14 = reg("w7_pv14", sf(Pv["A_observed"]["day14"]["mean"]), Pv["A_observed"]["day14"]["mean"])
    ks = 12.0 if b.LANG == "EN" else 11.5
    key_box(s, box, tg, [P(Rn(T(f"By minute 7.5 the mice had earned {f75}% of the day's water drops, yet licking was "
                                f"unchanged ({h1}).",
                                f"7.5분까지 그날 물 drop의 {f75}%를 얻었는데도 licking은 그대로였다 ({h1})."),
                              ks, True, "1F4E79"))],
            tip=tg["tip"], name="KEY", anchor="left")
    common(s, T("Day 14: no clear fall in licking in the first 7.5 min",
                "Day 14: 처음 7.5분 동안 licking 감소가 뚜렷하지 않았다"),
           T("compared each 30-s trial with the same trial on day 13 in mice and fitted models; tested per mouse.",
             "30초 trial마다 day 13 같은 trial과 비교 (마우스, 피팅 모델); 마우스별 검정."),
           T(f"day-14 licking changed {d14}, less than all models predicted (by {dmin}–{dmax}); only instant switch nears "
             f"rejection (mouse p = {pi_m}, simulated {pi_lo}–{pi_hi}); within day 14 vs overnight cannot be separated.",
             f"day 14 변화 {d14}: 모든 모델 예측보다 작음 (차이 {dmin}–{dmax}); 기각에 가까운 것은 즉시 전환뿐 (마우스별 "
             f"p = {pi_m}, 시뮬레이션 {pi_lo}–{pi_hi}); 세션 안 vs 밤사이는 구분 불가."),
           T(f"Trial = 30-s bin; Δ = any lick minus the same bin on day 13. n = {NA} Active, {NP} Passive. Models = group "
             f"fits. c: dots = mice, ◆ = mean ± 95% CI, paired t. Noise-matched simulations, p: {sim_txt_en}. First 2 min "
             f"at ceiling. 95% CI: within {wl_} to {wh_}, overnight {nl_} to {nh_}. Passive day 14 {pv14}.",
             f"Trial = 30초 bin; Δ = lick 유무 − day 13 같은 bin. n = Active {NA}, Passive {NP}. 모델 = 그룹 피팅. c: 점 = 마우스, "
             f"◆ = 평균 ± 95% CI, 대응 t. 잡음 맞춘 시뮬레이션 p: {sim_txt_kr}. 처음 2분은 천장. 95% CI: 세션 안 {wl_} ~ {wh_}, "
             f"밤사이 {nl_} ~ {nh_}. Passive day 14 {pv14}."),
           T(notes7_en(), notes7_kr()))
    return s


def _sim(k, key="day14", g="Active"):
    a = S15[g]["models"][k]["A"][key]
    return a["calibrated"]["p"], a["calibrated_B"]["p"], a["naive"]["p"]


def notes7_en():
    A, Pv, A15 = S["Active"], S["Passive"], S15["Active"]
    o, dr = A["A_observed"], A["day14_drops"]
    ca = A15["models"]["RPE2a"]["calibration"]
    return "\n".join([
        "WHAT: the paper's behavioural point is that mice do not change right after the reversal, when the errors are "
        "largest (the pattern it uses to explain Fig. 1h; the Fig. 1h correlation itself is not computed here). Same "
        "question for the drug -> water switch, trial by trial (30-s bins = the trials the models were fitted on).",
        f"a: drops are reward deliveries earned by licking (drinking was not measured). Active earned {dr['total_mean']:.1f} "
        f"water drops on day 14, {dr['frac_by_2p5min'][0] * 100:.0f}% by 2.5 min and {dr['frac_by_2p5min'][2] * 100:.0f}% by "
        f"7.5 min. This front-loading is the PR schedule (the ratio starts at 1 lick): on day 13 (morphine) "
        f"{S15['Active']['drops_frac_by_2p5min']['13'] * 100:.0f}% came by 2.5 min (dashed). Every model updates at each drop, "
        "so the models expect most of the day-14 change early.",
        f"b: P(lick) is 1.0 in the first 2 min on days 13 and 14 (ceiling), so a fall cannot show there. From 2.5 to 7.5 min "
        f"licking was unchanged; one low block at 7.5-10 min ({A['profile25_observed']['mean'][3]:+.2f}); the last 5 min were "
        f"back to {A['profile25_observed']['mean'][4]:+.2f} / {A['profile25_observed']['mean'][5]:+.2f}; day 15 was lower "
        "throughout.",
        "c: mice minus model for the whole day-14 change, per mouse, paired t (n = 6): "
        + "; ".join(f"{LAB[k]} {DR[k]['mean']:+.3f} (95% CI {DR[k]['lo']:+.2f} to {DR[k]['hi']:+.2f}), p {fp(DR[k]['p'])}"
                    for k in MAIN)
        + ". Noise-matched simulations (two calibrations of day-to-day and within-session noise, fitted to days 11-13; "
          "range = calibration A, B): "
        + "; ".join(f"{LAB[k]} p {fp2(_sim(k)[0])}, {fp2(_sim(k)[1])}" for k in MAIN)
        + f". The v1 check with independent bins gave p {fp2(_sim(MAIN[0])[2])} for instant switch but was too narrow: it "
          f"rejected a no-change placebo (Passive day 12 -> 13, p {fp2(S15['Passive']['models']['RPE2a']['placebo']['whole']['naive']['p'])}). "
          "Not on the slide: eMack and no learning also pass day 14; no learning fails days 15-16.",
        f"Within vs overnight: within day 14 {o['within14_halves']['mean']:+.3f} (95% CI {o['within14_halves']['lo']:+.2f} to "
        f"{o['within14_halves']['hi']:+.2f}), overnight {o['night_halves']['mean']:+.3f} ({o['night_halves']['lo']:+.2f} to "
        f"{o['night_halves']['hi']:+.2f}); the two are near mirror images across mice, and 5-min windows split the change "
        f"equally ({o['within14_5min']['mean']:+.3f} / {o['night_5min']['mean']:+.3f}). This design cannot place the fall.",
        f"Passive: licking rose at once on day 14 ({Pv['A_observed']['day14']['mean']:+.2f}, first half "
        f"{Pv['A_observed']['d14_h1']['mean']:+.2f}), more than every model predicted (mouse-level p "
        f"{fp(S15['Passive']['models']['RPE2a']['A']['day14']['mouse_t']['p_t'])} for RPE 2alpha); likely taste (beta < 0).",
        "MEANING: descriptively, the fall came after most water drops (the paper's pattern), but with 6 mice and one switch "
        "no model is rejected, and a delayed oral drug effect would give the same lag. Noise model details: m15 "
        f"(e.g. RPE 2alpha calibration A: sd {ca['A']['sig_u']:.1f} per day, {ca['A']['sig_a']:.1f} within session, phi "
        f"{ca['A']['phi']:.1f})."])


def notes7_kr():
    A, Pv = S["Active"], S["Passive"]
    o, dr = A["A_observed"], A["day14_drops"]
    return "\n".join([
        "[이 슬라이드] 논문의 행동 포인트는 '전환 직후, 오차가 가장 클 때는 행동이 바뀌지 않는다'이다 (논문이 Fig. 1h를 설명할 때 "
        "쓰는 모양; Fig. 1h의 상관 자체는 여기서 계산하지 않음). 같은 질문을 약 → 물 전환에 trial 단위로 던졌다 (30초 bin = "
        "모델을 피팅한 trial).",
        f"[a 읽는 법] drop은 licking으로 얻은 보상 전달이다 (마셨는지는 측정 안 함). Active는 day 14에 물 drop을 "
        f"{dr['total_mean']:.1f}개 얻었고, 2.5분까지 {dr['frac_by_2p5min'][0] * 100:.0f}%, 7.5분까지 "
        f"{dr['frac_by_2p5min'][2] * 100:.0f}%였다. 이렇게 앞에 몰리는 것은 PR 때문이다 (lick 1번부터 시작): day 13(모르핀)도 "
        f"2.5분까지 {S15['Active']['drops_frac_by_2p5min']['13'] * 100:.0f}% (점선). 모든 모델은 drop마다 배우므로 day 14 변화를 "
        "대부분 앞에서 예상한다.",
        f"[b 읽는 법] day 13과 14의 처음 2분은 P(lick) = 1.0 (천장)이라 감소가 보일 수 없다. 2.5–7.5분에는 그대로였고, "
        f"7.5–10분에 한 구간만 낮았다 ({A['profile25_observed']['mean'][3]:+.2f}); 마지막 5분은 다시 "
        f"{A['profile25_observed']['mean'][4]:+.2f} / {A['profile25_observed']['mean'][5]:+.2f}; day 15는 내내 낮았다.",
        "[c 읽는 법] day 14 전체 변화에서 마우스 − 모델, 마우스별, 대응 t (n = 6): "
        + "; ".join(f"{LAB_KR[k]} {DR[k]['mean']:+.3f} (95% CI {DR[k]['lo']:+.2f} ~ {DR[k]['hi']:+.2f}), p {fp(DR[k]['p'])}"
                    for k in MAIN)
        + ". 잡음 맞춘 시뮬레이션 (day 11–13에 맞춘 날 사이·세션 안 잡음, 보정 A, B 두 가지): "
        + "; ".join(f"{LAB_KR[k]} p {fp2(_sim(k)[0])}, {fp2(_sim(k)[1])}" for k in MAIN)
        + f". v1에서 쓴 독립 bin 검정은 즉시 전환 p {fp2(_sim(MAIN[0])[2])}였지만 너무 좁았다: 변화가 없어야 할 placebo "
          f"(Passive day 12 → 13)도 기각했다 (p {fp2(S15['Passive']['models']['RPE2a']['placebo']['whole']['naive']['p'])}). "
          "슬라이드에 없는 모델: eMack과 no learning도 day 14는 통과; no learning은 day 15–16에서 탈락.",
        f"[세션 안 vs 밤사이] day 14 안 {o['within14_halves']['mean']:+.3f} (95% CI {o['within14_halves']['lo']:+.2f} ~ "
        f"{o['within14_halves']['hi']:+.2f}), 밤사이 {o['night_halves']['mean']:+.3f} ({o['night_halves']['lo']:+.2f} ~ "
        f"{o['night_halves']['hi']:+.2f}); 두 값은 마우스마다 거의 반대로 움직이고, 5분 구간으로 나누면 변화가 반반 "
        f"({o['within14_5min']['mean']:+.3f} / {o['night_5min']['mean']:+.3f}). 이 설계로는 감소 시점을 정할 수 없다.",
        f"[Passive] day 14에 licking이 바로 올랐다 ({Pv['A_observed']['day14']['mean']:+.2f}, 전반 "
        f"{Pv['A_observed']['d14_h1']['mean']:+.2f}); 모든 모델 예측보다 컸다 (RPE 2α 마우스별 p "
        f"{fp(S15['Passive']['models']['RPE2a']['A']['day14']['mouse_t']['p_t'])}). 맛 때문으로 보인다 (beta < 0).",
        "[의미] 기술적으로는 감소가 물 drop 대부분 뒤에 왔다 (논문과 같은 모양). 하지만 마우스 6마리, 전환 한 번으로는 어떤 모델도 "
        "기각되지 않고, 늦게 오는 약효(마신 모르핀은 몇 분 뒤 작용)도 같은 지연을 만든다. 잡음 모델 자세한 내용은 m15."])


# ============================================================================================ slide 8 v2
def slide8(prs, F):
    b.SLIDE = 8
    path, size, _ = F
    s = b.new_slide(prs)
    place(s, path, size, 0.30, 0.98, 12.70)
    A, Pv, A15 = S["Active"], S["Passive"], S15["Active"]
    gc = A["sig_obs"]["groupcurve"]
    md = reg("w8_md", f"{gc['deg']['mid']:.1f}", gc["deg"]["mid"])
    me = reg("w8_me", f"{gc['enh']['mid']:.1f}", gc["enh"]["mid"])
    check("w8: group-curve midpoint later in > 99% of bootstraps", A["groupcurve_boot"]["frac_deg_mid_later"] > 0.99)
    m60 = reg("w8_m60", f"{A15['equal_window']['group_deg60_mid']:.1f}", A15["equal_window"]["group_deg60_mid"])
    tA = A["sig_tests"]["deg_vs_enh_mid_all"]
    pmA = reg("w8_pmA", fp(tA["p_t"]), tA["p_t"])
    p60 = reg("w8_p60", fp(A15["equal_window"]["paired_deg60_vs_enh"]["p_t"]), A15["equal_window"]["paired_deg60_vs_enh"]["p_t"])
    st = reg("w8_st", f"{gc['deg']['start']:.1f}", gc["deg"]["start"])
    slo, shi = A["groupcurve_boot"]["deg_start_CI"][0], A["groupcurve_boot"]["deg_start_CI"][2]
    s_lo, s_hi = reg("w8_slo", f"{slo:.0f}", slo), reg("w8_shi", f"{shi:.0f}", shi)
    pst = {k: (A15["models"][k]["timing"]["g_deg_start"]["calibrated"]["p"],
               A15["models"][k]["timing"]["g_deg_start"]["calibrated_B"]["p"]) for k in MAIN}
    pmin_st = min(min(v) for v in pst.values())
    ex_en = (f"excludes no model (p ≥ {fp(pmin_st)})" if pmin_st >= 0.05 else
             f"is later than most instant-switch simulations (p = {fp(pst[MAIN[0]][0])} / {fp(pst[MAIN[0]][1])}), "
             f"not robust")
    ex_kr = (f": 배제되는 모델 없음 (p ≥ {fp(pmin_st)})" if pmin_st >= 0.05 else
             f": 즉시 전환 시뮬레이션 대부분보다 늦지만 (p = {fp(pst[MAIN[0]][0])} / {fp(pst[MAIN[0]][1])}) 견고하지 않음")
    if pmin_st >= 0.05:
        reg("w8_pminst", fp(pmin_st), pmin_st)
    loo = A15["leave_one_out"]
    l_lo, l_hi = min(x["deg_start"] for x in loo), max(x["deg_start"] for x in loo)
    ll, lh = reg("w8_llo", f"{l_lo:.0f}", l_lo), reg("w8_lhi", f"{l_hi:.0f}", l_hi)
    enh_pred = {k: A["A_models"][k]["sig_pred_groupcurve"]["enh_mid"] for k in ("RPE2a", "mRPE")}
    e2a, emr = reg("w8_e2a", f"{enh_pred['RPE2a']:.1f}", enh_pred["RPE2a"]), reg("w8_emr", f"{enh_pred['mRPE']:.1f}", enh_pred["mRPE"])
    pmP = reg("w8_pmP", fp(Pv["sig_tests"]["deg_vs_enh_mid_all"]["p_t"]), Pv["sig_tests"]["deg_vs_enh_mid_all"]["p_t"])
    common(s, T("Timing: drug → water later than water → drug (group curve)",
                "전환 시점: 그룹 곡선에서 약 → 물이 물 → 약보다 늦다"),
           T("sigmoid fits by trial from each switch (paper Fig. 1d): midpoint (1e), 10% start (1m); same on model "
             "simulations.",
             "trial별 sigmoid 피팅 (논문 Fig. 1d): midpoint (1e), 10% 시작점 (1m); 모델 시뮬레이션도 같게."),
           T(f"group-curve midpoint trial {md} vs {me} (> 99% of bootstraps; 60-trial windows {m60} vs {me}), per mouse "
             f"p = {pmA} (n = {NA}); the start ({st}, 95% CI {s_lo}–{s_hi}) {ex_en}.",
             f"그룹 곡선 midpoint {md} vs {me} trial (bootstrap 99% 이상; 60 trial 창 {m60}), 마우스별 p = {pmA}; "
             f"시작점 {st} (95% CI {s_lo}–{s_hi}){ex_kr}."),
           T(f"Sigmoid A/(1 + e^(−(x − x0)/s)), |A| ≤ 1; midpoint = x0, start = x0 − s·ln 9. a: raw Δ P(lick), not "
             f"normalised; water → drug midpoint is an upper bound (first 2 min near ceiling) and predicted by the models "
             f"(RPE 2α {e2a}, meta-RPE {emr}). b: paired t (60-trial windows p = {p60}; Passive p = {pmP}); per-mouse fits are "
             f"steps, start ≈ midpoint. c: bars = 95% range in noise-matched simulations (A), ○ = median; p = A (B). "
             f"Leave-one-mouse-out start {ll}–{lh}.",
             f"Sigmoid A/(1 + e^(−(x − x0)/s)), |A| ≤ 1; midpoint = x0, 시작 = x0 − s·ln 9. a: 정규화 안 한 Δ P(lick); 물 → 약 "
             f"midpoint는 상한값 (천장)이고 모델도 예측 (RPE 2α {e2a}, meta-RPE {emr}). b: 대응 t (60 trial 창 p = {p60}; "
             f"Passive p = {pmP}); 마우스별 피팅은 계단 (시작 ≈ midpoint). c: 막대 = 잡음 맞춘 시뮬레이션(A) 95% 범위, "
             f"○ = 중앙값; p = A (B). 한 마리씩 빼면 시작 {ll}–{lh}."),
           T(notes8_en(), notes8_kr()))
    return s


def notes8_en():
    A, Pv, A15 = S["Active"], S["Passive"], S15["Active"]
    gc, gb = A["sig_obs"]["groupcurve"], A["groupcurve_boot"]
    tA, tP = A["sig_tests"]["deg_vs_enh_mid_all"], Pv["sig_tests"]["deg_vs_enh_mid_all"]
    pmA = PM[PM.group == "Active"]
    loo = A15["leave_one_out"]
    tim = lambda k, key: A15["models"][k]["timing"][key]
    return "\n".join([
        "WHAT: the paper's Fig. 1d,e,m on our trials. d: a sigmoid fitted to licking across trials from the reversal; e: its "
        "midpoint per mouse (H->L later than L->H, P = 0.027, n = 20, both cues reversed in the same session); m: the trial "
        "where the change starts, compared with each model's value curve (only meta-RPE was not significantly different, "
        "P = 0.21, n = 40).",
        f"a: Active group-mean curves, raw change in P(lick) against the day before the switch (not normalised; the paper's "
        f"Fig. 1d is normalised and smoothed). Drug -> water midpoint trial {gc['deg']['mid']:.1f} (95% CI "
        f"{gb['deg_mid_CI'][0]:.1f}-{gb['deg_mid_CI'][2]:.1f}; trial 30 = start of day 15), water -> drug "
        f"{gc['enh']['mid']:.1f} ({gb['enh_mid_CI'][0]:.1f}-{gb['enh_mid_CI'][2]:.1f}). Drug -> water was later in "
        f"{gb['frac_deg_mid_later'] * 100:.1f}% of 2000 bootstraps over mice (a share, not a p value). Windows differ "
        f"(90 vs 60 trials); with equal 60-trial windows the drug -> water midpoint is {A15['equal_window']['group_deg60_mid']:.1f} "
        f"(the fit is flat for midpoints about 28-50). Leave-one-out midpoints "
        f"{min(x['deg_mid'] for x in loo):.1f}-{max(x['deg_mid'] for x in loo):.1f} vs "
        f"{min(x['enh_mid'] for x in loo):.1f}-{max(x['enh_mid'] for x in loo):.1f}: the direction is robust. The water -> "
        "drug midpoint is an upper bound (licking near ceiling in the first 2 min) and every model predicts it.",
        f"b: per mouse, the paper's test. Active drug -> water {', '.join(f'{v:.1f}' for v in pmA.deg_mid)}; water -> drug "
        f"{', '.join(f'{v:.1f}' for v in pmA.enh_mid)} (paired t p {fp(tA['p_t'])}, Wilcoxon p {fp(tA['p_wilcoxon'])}, "
        f"{tA['n_a_greater']}/{tA['n']} later; equal windows p {fp(A15['equal_window']['paired_deg60_vs_enh']['p_t'])}). "
        f"Passive {tP['mean_a']:.1f} vs {tP['mean_b']:.1f} (p {fp(tP['p_t'])}): at once at drug -> water; the water -> drug "
        "change is gradual (group fit at its amplitude bound). Per-mouse fits are steps (s at its 0.3-trial floor), so a "
        "mouse's start is its midpoint minus 0.66.",
        "c: start of the drug -> water change in the group curve vs simulated data from each fitted model with noise matched "
        "to days 11-13 (calibration A; B in brackets): " + "; ".join(
            f"{LAB[k]} median {tim(k, 'g_deg_start')['cal_median']:.1f}, 95% range {tim(k, 'g_deg_start')['calibrated']['lo']:.1f}-"
            f"{tim(k, 'g_deg_start')['calibrated']['hi']:.1f}, p {fp(tim(k, 'g_deg_start')['calibrated']['p'])} "
            f"({fp(tim(k, 'g_deg_start')['calibrated_B']['p'])})" for k in MAIN)
        + ". Midpoint, same check: " + "; ".join(
            f"{LAB[k]} p {fp(tim(k, 'g_deg_mid')['calibrated']['p'])} ({fp(tim(k, 'g_deg_mid')['calibrated_B']['p'])})" for k in MAIN)
        + f". The mice's own start is uncertain: 95% CI {gb['deg_start_CI'][0]:.0f}-{gb['deg_start_CI'][2]:.0f}, leave-one-out "
          f"{min(x['deg_start'] for x in loo):.1f}-{max(x['deg_start'] for x in loo):.1f}.",
        "MEANING: at the group level the drug -> water change is later than the water -> drug change (as in the paper), but "
        "the paper's per-mouse test is not significant here (n = 6), and the start does not separate the models. Our two "
        "switches are 3 days apart with different histories (8 drug days vs 3 water days), unlike the paper's single "
        "session. The fast water -> drug change happened within ~2.5 min, before an oral dose acts; taste may act as a "
        "learned drug cue (not tested)."])


def notes8_kr():
    A, Pv, A15 = S["Active"], S["Passive"], S15["Active"]
    gc, gb = A["sig_obs"]["groupcurve"], A["groupcurve_boot"]
    tA, tP = A["sig_tests"]["deg_vs_enh_mid_all"], Pv["sig_tests"]["deg_vs_enh_mid_all"]
    pmA = PM[PM.group == "Active"]
    loo = A15["leave_one_out"]
    tim = lambda k, key: A15["models"][k]["timing"][key]
    return "\n".join([
        "[이 슬라이드] 논문 Fig. 1d,e,m을 우리 trial에 적용했다. d: 전환부터 trial별 licking에 sigmoid 피팅; e: 마우스별 "
        "midpoint (논문: H→L이 L→H보다 늦음, P = 0.027, n = 20, 두 cue가 같은 세션에서 반전); m: 변화가 시작되는 trial을 각 "
        "모델의 가치 곡선과 비교 (논문: meta-RPE만 유의한 차이 없음, P = 0.21, n = 40).",
        f"[a 읽는 법] Active 그룹 평균 곡선. 전환 전날과의 P(lick) 차이 그대로이며 정규화하지 않았다 (논문 Fig. 1d는 정규화·스무딩). "
        f"약 → 물 midpoint trial {gc['deg']['mid']:.1f} (95% CI {gb['deg_mid_CI'][0]:.1f}–{gb['deg_mid_CI'][2]:.1f}; trial 30 = "
        f"day 15 시작), 물 → 약 {gc['enh']['mid']:.1f} ({gb['enh_mid_CI'][0]:.1f}–{gb['enh_mid_CI'][2]:.1f}). 마우스 bootstrap "
        f"2000번 중 {gb['frac_deg_mid_later'] * 100:.1f}%에서 약 → 물이 더 늦다 (비율이지 p값이 아님). 창 길이가 다르다 (90 vs 60 "
        f"trial); 같은 60 trial 창이면 약 → 물 midpoint는 {A15['equal_window']['group_deg60_mid']:.1f} (midpoint 28–50 사이에서 "
        f"피팅이 거의 같음). 한 마리씩 빼면 midpoint {min(x['deg_mid'] for x in loo):.1f}–{max(x['deg_mid'] for x in loo):.1f} vs "
        f"{min(x['enh_mid'] for x in loo):.1f}–{max(x['enh_mid'] for x in loo):.1f}: 방향은 견고하다. 물 → 약 midpoint는 상한값이고 "
        "(처음 2분 licking이 천장 근처) 모든 모델이 예측한다.",
        f"[b 읽는 법] 논문의 마우스별 검정. Active 약 → 물 {', '.join(f'{v:.1f}' for v in pmA.deg_mid)}; 물 → 약 "
        f"{', '.join(f'{v:.1f}' for v in pmA.enh_mid)} (대응 t p {fp(tA['p_t'])}, Wilcoxon p {fp(tA['p_wilcoxon'])}, "
        f"{tA['n']}마리 중 {tA['n_a_greater']}마리가 더 늦음; 같은 창 p {fp(A15['equal_window']['paired_deg60_vs_enh']['p_t'])}). "
        f"Passive {tP['mean_a']:.1f} vs {tP['mean_b']:.1f} (p {fp(tP['p_t'])}): 약 → 물은 바로, 물 → 약은 천천히 (그룹 피팅이 "
        "진폭 한계에 걸림). 마우스별 피팅은 계단 모양이라 (s가 최소값 0.3) 시작점 = midpoint − 0.66.",
        "[c 읽는 법] 그룹 곡선의 약 → 물 시작점을, day 11–13에 맞춘 잡음을 넣어 각 모델로 만든 데이터와 비교 (보정 A; 괄호 B): "
        + "; ".join(
            f"{LAB_KR[k]} 중앙값 {tim(k, 'g_deg_start')['cal_median']:.1f}, 95% 범위 {tim(k, 'g_deg_start')['calibrated']['lo']:.1f}–"
            f"{tim(k, 'g_deg_start')['calibrated']['hi']:.1f}, p {fp(tim(k, 'g_deg_start')['calibrated']['p'])} "
            f"({fp(tim(k, 'g_deg_start')['calibrated_B']['p'])})" for k in MAIN)
        + ". midpoint로 같은 검정: " + "; ".join(
            f"{LAB_KR[k]} p {fp(tim(k, 'g_deg_mid')['calibrated']['p'])} ({fp(tim(k, 'g_deg_mid')['calibrated_B']['p'])})" for k in MAIN)
        + f". 마우스 시작점 자체가 불확실하다: 95% CI {gb['deg_start_CI'][0]:.0f}–{gb['deg_start_CI'][2]:.0f}, 한 마리씩 빼면 "
          f"{min(x['deg_start'] for x in loo):.1f}–{max(x['deg_start'] for x in loo):.1f}.",
        "[의미] 그룹 수준에서는 약 → 물 변화가 물 → 약보다 늦다 (논문과 같은 방향). 하지만 논문의 마우스별 검정은 여기서 유의하지 "
        "않고 (n = 6), 시작점으로는 모델을 구분할 수 없다. 우리 두 전환은 3일 떨어져 있고 이력도 다르다 (약 8일 vs 물 3일); 논문은 "
        "한 세션 안. 물 → 약 변화는 약 2.5분 안에 일어났는데, 마신 모르핀이 작용하기 전이다. 맛이 학습된 약물 신호일 수 있다 "
        "(검정 안 함)."])


# ============================================================================================ main
def main():
    F7, F8 = fig7(), fig8()
    outs, v4 = {}, {}
    for lang in ("EN", "KR"):
        b.LANG = lang
        two = Presentation(str(b.TEMPLATE))
        slide7(two, F7)
        slide8(two, F8)
        outs[lang] = OUT / f"mrpe_trial_by_trial_2slides_{lang}_v2_20261118.pptx"
        two.save(outs[lang])
        for fn, F, nm in ((slide7, F7, "7_day14"), (slide8, F8, "8_switch_timing")):
            one = Presentation(str(b.TEMPLATE))
            fn(one, F)
            one.save(OUT / "single_slides_v2" / f"{lang}_slide{nm}_v2_20261118.pptx")
        v2p = V2 / f"mrpe_hjort2026_opioid_6slides_{lang}_v2_20261117.pptx"
        deck = Presentation(str(v2p))
        n0 = len(deck.slides)
        slide7(deck, F7)
        slide8(deck, F8)
        check(f"[{lang}] v4 deck = 6 v2 slides + 2 new", n0 == 6 and len(deck.slides) == 8)
        v4[lang] = OUT / f"mrpe_hjort2026_opioid_8slides_{lang}_v4_20261118.pptx"
        deck.save(v4[lang])
        check(f"[{lang}] v4 deck: every v2 part except the 3 index parts is byte-identical", v1.restore_v2_parts(v2p, v4[lang]))
    for k, (txt, val) in b.PRINTED.items():
        if not k.startswith("w"):
            continue
        t_ = txt.replace("−", "-").replace("+", "").strip()
        if t_.startswith("<"):
            ok = val < float(t_[1:])
        else:
            nd = len(t_.split(".")[1]) if "." in t_ else 0
            ok = abs(float(t_) - val) <= 0.5 * 10 ** -nd + 1e-12
        check(f"v2 printed {k} = {txt} matches {val:.6g}", ok)
    json.dump(dict(checks=b.CHECKS, fit=b.FIT), open(OUT / "fit_report_v2.json", "w", encoding="utf-8"), indent=1,
              ensure_ascii=False)
    bad = [c_ for c_ in b.CHECKS if not c_[1]]
    print(f"{len(b.CHECKS)} checks, {len(bad)} failed")
    for c_ in bad:
        print("  FAILED", c_[0])
    for lang in ("EN", "KR"):
        print(lang, [p_.name for p_ in b.render(outs[lang], OUT / f"preview_v2_{lang}")])
        print(lang, "v4", [p_.name for p_ in b.render(v4[lang], OUT / f"preview_v4_{lang}")])


if __name__ == "__main__":
    main()
