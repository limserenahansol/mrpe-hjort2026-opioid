"""Step 14 (18 Nov 2026): two NEW slides (7, 8) with the trial-by-trial results of m13, EN + KR, python-pptx only.

Slides 1-6 (v2) are not changed. Outputs (slides_20261118_trial_by_trial\\):
  mrpe_trial_by_trial_2slides_{EN,KR}_20261118.pptx   the two new slides only
  single_slides\\{EN,KR}_slide7_..., _slide8_...        one slide per file
  mrpe_hjort2026_opioid_8slides_{EN,KR}_v3_20261118.pptx  = the v2 deck (slides 1-6 untouched) + slides 7-8
  preview_{EN,KR}\\slide7.png, slide8.png; preview_v3_{EN,KR}\\; figs\\; fit_report.json
Helpers (text boxes measured with the real fonts, KEY boxes, footer, number registry, LibreOffice render) are those of
m10_make_slides.py, used exactly as in m12 (v2). Every printed number is registered and checked against m13 outputs.
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m10_make_slides as b  # noqa: E402

plt, Line2D, Presentation = b.plt, b.Line2D, b.Presentation
P, Rn, tbox, measure, at, place, key_box, common = b.P, b.Rn, b.tbox, b.measure, b.at, b.place, b.key_box, b.common
frac, sf, fp, peq, reg, check, T = b.frac, b.sf, b.fp, b.peq, b.reg, b.check, b.T
D = b.D
OUT = HERE / "slides_20261118_trial_by_trial"
FIGD = OUT / "figs"
for p_ in (OUT, FIGD, OUT / "single_slides"):
    p_.mkdir(parents=True, exist_ok=True)
V2 = HERE / "slides_20261117"
RED, BLUE, GREY, GREEN, ORANGE, PURPLE = b.RED, b.BLUE, b.GREY, b.GREEN, b.ORANGE, b.PURPLE
S = json.load(open(D / "m13_summary.json"))
PM = pd.read_csv(D / "m13_per_mouse.csv")
BINS = pd.read_csv(D / "m13_bins.csv")
NSIM = S["meta"]["n_sim"]
MAIN = ["today's outcome", "RPE", "RPE2a", "mRPE"]
LAB = {"today's outcome": "instant switch", "RPE": "RPE", "RPE2a": "RPE 2α", "mRPE": "meta-RPE"}
LAB_KR = LAB
CM = {"today's outcome": PURPLE, "RPE": GREY, "RPE2a": ORANGE, "mRPE": GREEN}
NA, NP = len(S["meta"]["mice"]["Active"]), len(S["meta"]["mice"]["Passive"])


def ppc_txt(p):
    """posterior predictive p: the smallest possible value is 1/(NSIM+1)"""
    return "p < 0.001" if p <= 1.0 / (NSIM + 1) + 1e-12 else peq(p)


def save(fig, name):
    p_ = FIGD / name
    fig.savefig(p_, facecolor="white")
    plt.close(fig)
    return p_


def traj(g, kind):
    """per-mouse D trajectories (n, 90 or 60) from the bins table, same definition as m13"""
    x = BINS[BINS.group == g].pivot_table(index=["mouse", "day"], columns="bin", values="lick")
    ms = S["meta"]["mice"][g]
    out = []
    for m in ms:
        X = x.loc[m]
        if kind == "deg":
            out.append(np.concatenate([X.loc[d].values - X.loc[13].values for d in (14, 15, 16)]))
        else:
            out.append(np.concatenate([X.loc[d].values - X.loc[16].values for d in (17, 18)]))
    return np.array(out)


# ============================================================================================ figure 7
def fig7():
    A = S["Active"]
    fig = plt.figure(figsize=(12.7, 4.85))
    # ---- a: water drops through day 14
    a = fig.add_axes([0.062, 0.17, 0.175, 0.70])
    for g, col in (("Active", RED), ("Passive", BLUE)):
        x = BINS[(BINS.group == g) & (BINS.day == 14)].pivot_table(index="mouse", columns="bin", values="drops_in_bin")
        cum = np.c_[np.zeros(len(x)), np.cumsum(x.values, 1)]
        tt = np.arange(31) * 0.5
        mu, se = cum.mean(0), cum.std(0, ddof=1) / np.sqrt(len(cum))
        a.plot(tt, mu, color=col, lw=2.2, label=f"{g} (n = {len(cum)})")
        a.fill_between(tt, mu - se, mu + se, color=col, alpha=0.18, lw=0)
    fa = A["day14_drops"]["frac_by_2p5min"][0]
    nA25 = A["day14_drops"]["by_2p5min"][0]
    a.axvline(2.5, color="#999999", lw=1, ls=":")
    a.set_xlim(0, 15)
    a.set_ylim(0, 34)
    a.set_xticks([0, 2.5, 5, 10, 15], ["0", "2.5", "5", "10", "15"])
    a.set_xlabel("Minute of the day-14 session")
    a.set_ylabel("Water drops drunk so far\n(day 14, mean ± s.e.m.)")
    a.set_title("a  When the water drops came", loc="left", fontweight="bold")
    a.legend(loc="lower right", frameon=False, fontsize=10.5)
    # ---- b: 2.5-min profile, Active, observed vs models
    bx = fig.add_axes([0.315, 0.17, 0.405, 0.70])
    ob = np.array(A["profile25_observed"]["mean"])
    se = np.array(A["profile25_observed"]["sem"])
    off = {14: 0.0, 15: 17.0, 16: 34.0}
    xs = np.concatenate([off[d] + np.arange(6) * 2.5 + 1.25 for d in (14, 15, 16)])
    for k in MAIN:
        pr = np.array(A["profile25_models"][k])
        for j, d in enumerate((14, 15, 16)):
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
    bx.legend(handles=hs, loc="lower left", frameon=True, fontsize=9.6, ncol=1, borderpad=0.35, framealpha=0.95,
              edgecolor="#DDDDDD", handlelength=1.6, borderaxespad=0.3)
    # ---- c: whole day 14, mice vs models
    c = fig.add_axes([0.80, 0.17, 0.185, 0.70])
    ob14 = A["A_observed"]["day14"]
    ys = {"mice": 0}
    ys.update({k: i + 1 for i, k in enumerate(MAIN)})
    pm = PM[PM.group == "Active"].obs_day14.values
    jit = np.linspace(-0.16, 0.16, len(pm))
    c.scatter(pm, np.zeros(len(pm)) + jit, s=22, color=RED, alpha=0.75, zorder=3, lw=0)
    c.errorbar(ob14["mean"], 0, xerr=[[ob14["mean"] - ob14["lo"]], [ob14["hi"] - ob14["mean"]]], color="black",
               marker="D", ms=7, capsize=4, lw=2, zorder=4)
    for k in MAIN:
        r = A["A_models"][k]["day14"]
        y = ys[k]
        c.errorbar(r["pred_mean"], y, xerr=[[r["pred_mean"] - r["sim_lo"]], [r["sim_hi"] - r["pred_mean"]]],
                   color=CM[k], marker="s", ms=7, capsize=4, lw=2.2)
        c.text(0.58, y, ppc_txt(r["ppc_p"]), ha="right", va="center", fontsize=10.5, color=CM[k])
    c.axvline(ob14["mean"], color="black", lw=1, ls=":")
    c.axvline(0, color="#777777", lw=1)
    c.set_yticks(list(ys.values()), ["mice", "instant\nswitch", "RPE", "RPE 2α", "meta-RPE"], fontsize=11)
    c.set_ylim(4.6, -0.6)
    c.set_xlim(-0.62, 0.6)
    c.set_xticks([-0.5, -0.25, 0, 0.25])
    c.set_xlabel("Day-14 licking minus day 13\n(whole session, Δ P(lick))")
    c.set_title("c  Whole day 14, Active", loc="left", fontweight="bold")
    tg = {"key_tl": frac(fig, bx, 18.4, 0.55), "key_br": frac(fig, bx, 49.0, 0.06),
          "tip": frac(fig, bx, 4.6, ob[1] + se[1] + 0.03)}
    return save(fig, "s7_day14_vs_overnight.png"), fig.get_size_inches(), tg


# ============================================================================================ figure 8
def fig8():
    A = S["Active"]
    fig = plt.figure(figsize=(12.7, 4.85))
    # ---- a: change by trial from the switch (group mean, Active): 2.5-min blocks +- s.e.m. and the sigmoid fits
    a = fig.add_axes([0.065, 0.17, 0.33, 0.70])
    gc = A["sig_obs"]["groupcurve"]
    YL = (-0.78, 0.78)
    for kind, col, mk, lab, ytxt in (("deg", "#222222", "o", "drug → water", -0.66),
                                     ("enh", "#8A8A8A", "^", "water → drug", -0.66)):
        Tm = traj("Active", kind)
        n = Tm.shape[1]
        blk = Tm.reshape(len(Tm), n // 5, 5).mean(2)
        mu, se = blk.mean(0), blk.std(0, ddof=1) / np.sqrt(len(blk))
        xc = np.arange(n // 5) * 5 + 2.0
        check(f"fig8a {kind}: block means inside the axis", (mu - se).min() > YL[0] and (mu + se).max() < YL[1])
        a.errorbar(xc, mu, yerr=se, color=col, marker=mk, ms=5.5, lw=0, elinewidth=1.3, capsize=2.5, zorder=3, label=lab)
        Am, x0, s_ = gc[kind]["A"], gc[kind]["x0"], gc[kind]["s"]
        xf = np.linspace(0, n - 1, 400)
        a.plot(xf, Am * expit((xf - x0) / s_), color=col, lw=2.8, zorder=2)
        a.axvline(x0, color=col, lw=1.2, ls="--")
        a.text(x0 + 1.2, ytxt, f"midpoint {x0:.1f}", color=col, fontsize=10.5, fontweight="bold", va="center")
    for xv in (30, 60):
        a.axvline(xv, color="#BBBBBB", lw=1, ls=":")
        a.text(xv - 0.6, 0.76, "next day", fontsize=9.5, color="#777777", va="top", ha="right", rotation=90)
    a.axhline(0, color="#777777", lw=1)
    a.set_xlim(-1, 90)
    a.set_ylim(*YL)
    a.set_xlabel("Trial (30-s bin) from the switch")
    a.set_ylabel("Licking minus the day before the switch\n(Δ P(lick), mean ± s.e.m. of 2.5 min)")
    a.set_title(f"a  Licking by trial, Active (n = {NA}); lines = fits", loc="left", fontweight="bold")
    a.legend(loc="upper right", frameon=True, fontsize=9.8, framealpha=0.95, edgecolor="#DDDDDD",
             borderaxespad=0.3)
    # ---- b: per-mouse midpoints (Fig. 1e)
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
    # ---- c: start of the change: mice vs models (Fig. 1m), group curve
    c = fig.add_axes([0.80, 0.17, 0.185, 0.70])
    obs = A["sig_obs"]["groupcurve"]["deg"]["start"]
    lo, med, hi = A["groupcurve_boot"]["deg_start_CI"]
    c.axvspan(lo, hi, color="#DDDDDD", lw=0, zorder=0)
    c.axvline(obs, color="black", lw=2, zorder=3)
    for i, k in enumerate(MAIN):
        sim = np.array(A["A_models"][k]["sig_sim_groupcurve"]["deg_start"])
        q = np.percentile(sim, [2.5, 50, 97.5])
        y = i
        c.plot([q[0], q[2]], [y, y], color=CM[k], lw=5, solid_capstyle="butt", alpha=0.85)
        c.plot(q[1], y, marker="o", color="white", mec=CM[k], mew=2, ms=8)
        c.text(89, y, ppc_txt(A["C_start_vs_models"][k]["ppc_groupcurve_start"]), ha="right", va="center",
               fontsize=10.5, color=CM[k])
    c.set_yticks(range(4), ["instant\nswitch", "RPE", "RPE 2α", "meta-RPE"], fontsize=11)
    c.set_ylim(3.6, -0.75)
    c.set_xlim(-2, 90)
    c.set_xticks([0, 30, 60, 90])
    c.text(obs + 1.5, -0.45, f"mice {obs:.1f}", fontsize=10.5, fontweight="bold")
    c.text(hi - 1.0, 3.45, "mice 95% CI", fontsize=9, color="#666666", ha="right", va="bottom")
    c.set_xlabel("Start of the change (trial)\ndrug → water, group curve")
    c.set_title("c  Start: mice vs models", loc="left", fontweight="bold")
    return save(fig, "s8_switch_timing.png"), fig.get_size_inches(), None


# ============================================================================================ slide 7
def slide7(prs, F):
    b.SLIDE = 7
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    A, Pv = S["Active"], S["Passive"]
    o = A["A_observed"]
    f25 = reg("t7_f25", f"{A['day14_drops']['frac_by_2p5min'][0] * 100:.0f}", A["day14_drops"]["frac_by_2p5min"][0] * 100)
    f75 = reg("t7_f75", f"{A['day14_drops']['frac_by_2p5min'][2] * 100:.0f}", A["day14_drops"]["frac_by_2p5min"][2] * 100)
    d14 = reg("t7_d14", sf(o["day14"]["mean"]), o["day14"]["mean"])
    mlo = reg("t7_mlo", sf(A["A_models"]["mRPE"]["day14"]["pred_mean"]), A["A_models"]["mRPE"]["day14"]["pred_mean"])
    mhi = reg("t7_mhi", sf(A["A_models"][MAIN[0]]["day14"]["pred_mean"]), A["A_models"][MAIN[0]]["day14"]["pred_mean"])
    check("t7 meta-RPE predicts the smallest and instant switch the largest day-14 fall",
          min(A["A_models"][k]["day14"]["pred_mean"] for k in MAIN) == A["A_models"][MAIN[0]]["day14"]["pred_mean"]
          and max(A["A_models"][k]["day14"]["pred_mean"] for k in MAIN) == A["A_models"]["mRPE"]["day14"]["pred_mean"])
    pk = {k: A["A_models"][k]["day14"]["ppc_p"] for k in MAIN}
    p_rpe = reg("t7_prpe", fp(pk["RPE"]), pk["RPE"])
    p_2a = reg("t7_p2a", fp(pk["RPE2a"]), pk["RPE2a"])
    p_m = reg("t7_pm", fp(pk["mRPE"]), pk["mRPE"])
    check("t7 instant-switch ppc p is the minimum (1/2001)", pk[MAIN[0]] <= 1 / (NSIM + 1) + 1e-12)
    h1 = reg("t7_h1", sf(o["d14_h1"]["mean"]), o["d14_h1"]["mean"])
    w = reg("t7_w", sf(o["within14_halves"]["mean"]), o["within14_halves"]["mean"])
    pw = reg("t7_pw", fp(o["within14_halves"]["p"]), o["within14_halves"]["p"])
    ng = reg("t7_ng", sf(o["night_halves"]["mean"]), o["night_halves"]["mean"])
    pn = reg("t7_pn", fp(o["night_halves"]["p"]), o["night_halves"]["p"])
    pv14 = reg("t7_pv14", sf(Pv["A_observed"]["day14"]["mean"]), Pv["A_observed"]["day14"]["mean"])
    ppv = reg("t7_ppv", fp(Pv["A_observed"]["day14"]["p"]), Pv["A_observed"]["day14"]["p"])
    ks = 12.0 if b.LANG == "EN" else 11.5
    key_box(s, box, tg, [P(Rn(T(f"By minute 7.5 the mice had drunk {f75}% of the day's water drops, yet licking was "
                                f"unchanged ({h1}).",
                                f"7.5분까지 그날 물 drop의 {f75}%를 마셨는데도 licking은 그대로였다 ({h1})."),
                              ks, True, "1F4E79"))],
            tip=tg["tip"], name="KEY", anchor="left")
    common(s, T("Day 14: licking lagged the water drops; no extra change overnight",
                "Day 14: licking 변화는 물 drop보다 늦었고, 밤사이 추가 변화는 없었다"),
           T("compared each 30-s trial with the same trial on day 13, in the mice and in each fitted model.",
             "30초 trial마다 day 13의 같은 trial과 비교 (마우스와 피팅된 각 모델)."),
           T(f"{f25}% of water drops came in the first 2.5 min, yet day-14 licking changed only {d14} (models {mlo} to "
             f"{mhi}; p: instant switch < 0.001, RPE {p_rpe}, RPE 2α {p_2a}, meta-RPE {p_m}); overnight {ng} (p = {pn}).",
             f"물 drop의 {f25}%가 처음 2.5분에 왔지만 day 14 변화는 {d14}뿐 (모델 {mlo} ~ {mhi}; p: 즉시 전환 < 0.001, "
             f"RPE {p_rpe}, RPE 2α {p_2a}, meta-RPE {p_m}); 밤사이 {ng} (p = {pn})."),
           T(f"Trial = 30-s bin, 30 per session; Δ = any-lick minus the same bin on day 13. n = {NA} Active, {NP} Passive "
             f"(6099_white left out). Models = group fits (shared learning rates, readout per mouse, days 11–18). Model "
             f"bars = 95% range of the group mean in 2000 sets simulated from each model; p = two-sided posterior predictive, "
             f"not refitted (conservative). Halves (one-sample t): within day 14 {w} (p = {pw}). Passive: day 14 {pv14} "
             f"(p = {ppv}).",
             f"Trial = 30초 bin, 세션당 30개; Δ = lick 유무 − day 13 같은 bin. n = Active {NA}, Passive {NP} (6099_white 제외). "
             f"모델 = 그룹 피팅 (학습률 공유, 마우스별 출력, day 11–18). 모델 막대 = 모델로 만든 2000세트의 그룹 평균 95% 범위; "
             f"p 양측 posterior predictive, 재피팅 없음 (보수적). 절반 비교 (단일표본 t): day 14 안 {w} (p = {pw}). "
             f"Passive: day 14 {pv14} (p = {ppv})."),
           T(notes7_en(), notes7_kr()))
    return s


def _m(k, key="day14"):
    r = S["Active"]["A_models"][k][key]
    return r["pred_mean"], r["sim_lo"], r["sim_hi"], r["ppc_p"]


def notes7_en():
    A, Pv = S["Active"], S["Passive"]
    o = A["A_observed"]
    dr = A["day14_drops"]
    lines = ["WHAT: the paper's key behavioural observation (Fig. 1h) is that mice do not change right after the reversal, "
             "when the errors are largest. Here we ask the same of the drug -> water switch, trial by trial (30-s bins, "
             "the trials the models were fitted on), and whether the change came within day 14 or only overnight.",
             f"a: Active mice drank {dr['total_mean']:.1f} water drops on day 14 ({dr['total_min']}-{dr['total_max']}); "
             f"{dr['frac_by_2p5min'][0] * 100:.0f}% of them in the first 2.5 min, {dr['frac_by_2p5min'][2] * 100:.0f}% by "
             "7.5 min (the PR ratio starts at 1 lick). Every model here updates at each drop, so the errors were "
             "concentrated at the start of the session.",
             "b: licking minus the same 2.5 min on day 13 (black), and each model's prediction (lines). The models expect "
             "the fall in the first minutes of day 14; the mice held until about minute 7.5. Day 15 is lower from its "
             "first minutes.",
             f"c and halves: whole day 14 {o['day14']['mean']:+.3f} (95% CI {o['day14']['lo']:+.2f} to "
             f"{o['day14']['hi']:+.2f}). Model predictions and posterior predictive p: "
             + "; ".join(f"{LAB[k]} {_m(k)[0]:+.3f} (range {_m(k)[1]:+.2f} to {_m(k)[2]:+.2f}, p {fp(_m(k)[3])})" for k in MAIN)
             + f". First half {o['d14_h1']['mean']:+.3f}, second half {o['d14_h2']['mean']:+.3f}, day 15 first half "
             f"{o['d15_h1']['mean']:+.3f}: within day 14 {o['within14_halves']['mean']:+.3f} (p {fp(o['within14_halves']['p'])}); "
             f"overnight {o['night_halves']['mean']:+.3f} (p {fp(o['night_halves']['p'])}); models predict "
             + ", ".join(f"{LAB[k]} {_m(k, 'night_halves')[0]:+.3f}" for k in MAIN) + " overnight. "
             f"5-min windows instead of halves: within {o['within14_5min']['mean']:+.3f} (p {fp(o['within14_5min']['p'])}), "
             f"overnight {o['night_5min']['mean']:+.3f} (p {fp(o['night_5min']['p'])}).",
             f"Checks: lick counts (log) give the same picture (day 14 {A['A_observed_counts_log1p']['day14']['mean']:+.2f}, "
             f"p {fp(A['A_observed_counts_log1p']['day14']['p'])}; overnight "
             f"{A['A_observed_counts_log1p']['night_halves']['mean']:+.2f}, p "
             f"{fp(A['A_observed_counts_log1p']['night_halves']['p'])}). Day 13 is a fair reference: its within-session "
             f"fall ({A['within_session_reference']['day13_vs_baseline45']['mean_a']:+.2f}) is like Baseline days 4-5 "
             f"({A['within_session_reference']['day13_vs_baseline45']['mean_b']:+.2f}, p "
             f"{fp(A['within_session_reference']['day13_vs_baseline45']['p_t'])}).",
             f"Passive: licking rose at once on day 14 ({Pv['A_observed']['day14']['mean']:+.2f}, p "
             f"{fp(Pv['A_observed']['day14']['p'])}; first half {Pv['A_observed']['d14_h1']['mean']:+.2f}), as expected "
             "from taste (beta < 0).",
             "MEANING: the change lagged the drops (the paper's Fig. 1h pattern) and then carried over the night. The lag is "
             "too long for an instant switch or one learning rate (RPE), borderline for two fixed rates (RPE 2alpha) and "
             "fits meta-RPE. CAVEATS: n = 6, one switch; posterior predictive without refitting; an oral morphine dose acts "
             "minutes after drinking, so the missing-drug error may only be felt later in the session (a delayed outcome "
             "gives the same lag as error averaging)."]
    return "\n".join(lines)


def notes7_kr():
    A, Pv = S["Active"], S["Passive"]
    o = A["A_observed"]
    dr = A["day14_drops"]
    lines = ["[이 슬라이드] 논문의 핵심 행동 관찰(Fig. 1h)은 '전환 직후, 오차가 가장 클 때는 행동이 바뀌지 않는다'이다. 같은 질문을 "
             "약 → 물 전환에 trial 단위로 던졌다 (30초 bin = 모델을 피팅한 trial). 그리고 변화가 day 14 세션 안에서 생겼는지, "
             "밤사이에만 생겼는지 봤다.",
             f"[a 읽는 법] Active는 day 14에 물 drop을 평균 {dr['total_mean']:.1f}개 마셨다 ({dr['total_min']}–{dr['total_max']}개). "
             f"그중 {dr['frac_by_2p5min'][0] * 100:.0f}%를 처음 2.5분에, {dr['frac_by_2p5min'][2] * 100:.0f}%를 7.5분까지 마셨다 "
             "(PR은 lick 1번부터 시작). 모든 모델은 drop마다 배우므로, 오차는 세션 처음에 몰려 있었다.",
             "[b 읽는 법] 검정 점 = day 13의 같은 2.5분 구간과의 licking 차이, 선 = 각 모델의 예측. 모델들은 day 14 처음 몇 분에 "
             "감소를 예상하지만, 마우스는 7.5분쯤까지 그대로였다. day 15는 처음 몇 분부터 낮다.",
             f"[c와 절반 비교] day 14 전체 {o['day14']['mean']:+.3f} (95% CI {o['day14']['lo']:+.2f} ~ {o['day14']['hi']:+.2f}). "
             "모델 예측과 posterior predictive p: "
             + "; ".join(f"{LAB[k]} {_m(k)[0]:+.3f} (범위 {_m(k)[1]:+.2f} ~ {_m(k)[2]:+.2f}, p {fp(_m(k)[3])})" for k in MAIN)
             + f". day 14 전반 {o['d14_h1']['mean']:+.3f}, 후반 {o['d14_h2']['mean']:+.3f}, day 15 전반 {o['d15_h1']['mean']:+.3f}: "
             f"day 14 안의 변화 {o['within14_halves']['mean']:+.3f} (p {fp(o['within14_halves']['p'])}), 밤사이 "
             f"{o['night_halves']['mean']:+.3f} (p {fp(o['night_halves']['p'])}). 모델의 밤사이 예측: "
             + ", ".join(f"{LAB[k]} {_m(k, 'night_halves')[0]:+.3f}" for k in MAIN) + ". "
             f"절반 대신 5분 구간: 세션 안 {o['within14_5min']['mean']:+.3f} (p {fp(o['within14_5min']['p'])}), 밤사이 "
             f"{o['night_5min']['mean']:+.3f} (p {fp(o['night_5min']['p'])}).",
             f"[확인] lick 수(log)로 봐도 같다 (day 14 {A['A_observed_counts_log1p']['day14']['mean']:+.2f}, p "
             f"{fp(A['A_observed_counts_log1p']['day14']['p'])}; 밤사이 {A['A_observed_counts_log1p']['night_halves']['mean']:+.2f}, "
             f"p {fp(A['A_observed_counts_log1p']['night_halves']['p'])}). day 13을 기준으로 써도 된다: day 13의 세션 안 감소 "
             f"({A['within_session_reference']['day13_vs_baseline45']['mean_a']:+.2f})는 Baseline day 4–5 "
             f"({A['within_session_reference']['day13_vs_baseline45']['mean_b']:+.2f})와 비슷하다 (p "
             f"{fp(A['within_session_reference']['day13_vs_baseline45']['p_t'])}).",
             f"[Passive] day 14에 licking이 바로 올랐다 ({Pv['A_observed']['day14']['mean']:+.2f}, p "
             f"{fp(Pv['A_observed']['day14']['p'])}; 전반 {Pv['A_observed']['d14_h1']['mean']:+.2f}). 맛 때문으로 보인다 (beta < 0).",
             "[의미] 변화는 drop보다 늦게 왔고 (논문 Fig. 1h와 같은 모양), 그 뒤 밤을 넘어 그대로 이어졌다. 이 지연은 즉시 전환이나 "
             "학습률 하나(RPE)로는 너무 길고, 고정 학습률 두 개(RPE 2α)로는 경계선이며, meta-RPE와 맞는다. [주의] n = 6, 전환 한 번; "
             "posterior predictive는 재피팅 없음; 마신 모르핀은 몇 분 뒤에 작용하므로 '약이 없다'는 오차는 세션 후반에야 느낄 수 "
             "있다 (늦게 오는 결과도 오차 평균과 같은 지연을 만든다)."]
    return "\n".join(lines)


# ============================================================================================ slide 8
def slide8(prs, F):
    b.SLIDE = 8
    path, size, _ = F
    s = b.new_slide(prs)
    place(s, path, size, 0.30, 0.98, 12.70)
    A, Pv = S["Active"], S["Passive"]
    gc = A["sig_obs"]["groupcurve"]
    md = reg("t8_md", f"{gc['deg']['mid']:.1f}", gc["deg"]["mid"])
    me = reg("t8_me", f"{gc['enh']['mid']:.1f}", gc["enh"]["mid"])
    fl = reg("t8_fl", f"{A['groupcurve_boot']['frac_deg_mid_later'] * 100:.1f}", A["groupcurve_boot"]["frac_deg_mid_later"] * 100)
    tA = A["sig_tests"]["deg_vs_enh_mid_all"]
    tP = Pv["sig_tests"]["deg_vs_enh_mid_all"]
    pmA = reg("t8_pmA", fp(tA["p_t"]), tA["p_t"])
    pmP = reg("t8_pmP", fp(tP["p_t"]), tP["p_t"])
    st = reg("t8_st", f"{gc['deg']['start']:.1f}", gc["deg"]["start"])
    pc = {k: A["C_start_vs_models"][k]["ppc_groupcurve_start"] for k in MAIN}
    p_i = reg("t8_pi", fp(pc[MAIN[0]]), pc[MAIN[0]])
    check("t8 only instant switch is below 0.05 for the start", pc[MAIN[0]] < 0.05 and all(pc[k] >= 0.05 for k in MAIN[1:]))
    blo, bhi = A["groupcurve_boot"]["deg_mid_CI"][0], A["groupcurve_boot"]["deg_mid_CI"][2]
    b_lo, b_hi = reg("t8_blo", f"{blo:.1f}", blo), reg("t8_bhi", f"{bhi:.1f}", bhi)
    elo, ehi = A["groupcurve_boot"]["enh_mid_CI"][0], A["groupcurve_boot"]["enh_mid_CI"][2]
    e_lo, e_hi = reg("t8_elo", f"{elo:.1f}", elo), reg("t8_ehi", f"{ehi:.1f}", ehi)
    nf = [A["C_start_vs_models"][k]["start_vs_model"] for k in MAIN]
    nf_hi = reg("t8_nfhi", f"{max(x['mean_b'] for x in nf):.1f}", max(x["mean_b"] for x in nf))
    nf_obs = reg("t8_nfobs", f"{nf[0]['mean_a']:.1f}", nf[0]["mean_a"])
    nf_h = reg("t8_nfh", fp(max(x["p_t_holm"] for x in nf)), max(x["p_t_holm"] for x in nf))
    check("t8 Holm p equal for the four models", len({round(x["p_t_holm"], 6) for x in nf}) == 1)
    common(s, T("Switch timing: drug → water late, water → drug fast",
                "전환 시점: 약 → 물은 늦고, 물 → 약은 빠르다"),
           T("fitted a sigmoid to licking by trial from each switch (paper Fig. 1d) for the midpoint (1e) and 10% start "
             "(1m); same on 2000 sets per model.",
             "trial별 licking에 sigmoid 피팅 (논문 Fig. 1d) → midpoint (1e), 10% 시작점 (1m); 모델별 2000세트도 같게."),
           T(f"Active midpoint trial {md} after drug → water vs {me} after water → drug (later in {fl}% of bootstraps; per "
             f"mouse p = {pmA}, n = {NA}); the start (trial {st}) rules out only instant switch (p = {p_i}).",
             f"Active midpoint: 약 → 물 trial {md} vs 물 → 약 {me} (bootstrap {fl}%; 마우스별 p = {pmA}); 시작점 trial {st}: "
             f"즉시 전환만 배제 (p = {p_i})."),
           T(f"Sigmoid A/(1 + e^(−(x − x0)/s)), grid least squares, |A| ≤ 1; midpoint = x0, start = x0 − s·ln 9 (10% of the "
             f"change). a: group-mean curve; midpoint 95% CI {b_lo}–{b_hi} vs {e_lo}–{e_hi} (2000 bootstraps over mice). b: "
             f"paired t (Passive p = {pmP}, n = {NP}). c: bars = 95% range of the start in 2000 sets per model, grey = mice 95% "
             f"CI; two-sided posterior predictive p. Noise-free model curves (paper style) start by trial {nf_hi} vs mice "
             f"{nf_obs} (Holm p = {nf_h}).",
             f"Sigmoid A/(1 + e^(−(x − x0)/s)), 격자 최소제곱, |A| ≤ 1; midpoint = x0, 시작 = x0 − s·ln 9 (변화의 10%). a: 그룹 "
             f"평균 곡선; midpoint 95% CI {b_lo}–{b_hi} vs {e_lo}–{e_hi} (bootstrap 2000번). b: 대응 t (Passive p = {pmP}, "
             f"n = {NP}). c: 막대 = 모델별 2000세트의 시작점 95% 범위, 회색 = 마우스 95% CI; p 양측. 잡음 없는 모델 곡선 "
             f"(논문 방식): trial {nf_hi} 이내 vs 마우스 {nf_obs} (Holm p = {nf_h})."),
           T(notes8_en(), notes8_kr()))
    return s


def notes8_en():
    A, Pv = S["Active"], S["Passive"]
    gc, gb = A["sig_obs"]["groupcurve"], A["groupcurve_boot"]
    tA, tP = A["sig_tests"]["deg_vs_enh_mid_all"], Pv["sig_tests"]["deg_vs_enh_mid_all"]
    pmA = PM[PM.group == "Active"]
    pmP = PM[PM.group == "Passive"]
    C = A["C_start_vs_models"]
    return "\n".join([
        "WHAT: the paper's Fig. 1d,e,m done on our trials. d: a sigmoid fitted to licking across trials from the reversal; "
        "e: its midpoint (H->L later than L->H, P = 0.027, n = 20); m: the trial where the change starts (10%), compared with "
        "each model (only meta-RPE matched the mice).",
        f"a: Active group-mean curve, normalised so 0 = before and 1 = after. Drug -> water midpoint trial {gc['deg']['mid']:.1f} "
        f"(95% CI {gb['deg_mid_CI'][0]:.1f}-{gb['deg_mid_CI'][2]:.1f}; trial 30 = start of day 15), start "
        f"{gc['deg']['start']:.1f}; water -> drug midpoint {gc['enh']['mid']:.1f} (CI {gb['enh_mid_CI'][0]:.1f}-"
        f"{gb['enh_mid_CI'][2]:.1f}). Drug -> water later in {gb['frac_deg_mid_later'] * 100:.1f}% of 2000 bootstraps over mice.",
        f"b: per mouse, as in the paper. Active drug -> water {', '.join(f'{v:.1f}' for v in pmA.deg_mid)}; water -> drug "
        f"{', '.join(f'{v:.1f}' for v in pmA.enh_mid)} (mean {tA['mean_a']:.1f} vs {tA['mean_b']:.1f}, paired t p "
        f"{fp(tA['p_t'])}, Wilcoxon p {fp(tA['p_wilcoxon'])}, {tA['n_a_greater']}/{tA['n']} later). Passive "
        f"{tP['mean_a']:.1f} vs {tP['mean_b']:.1f} (p {fp(tP['p_t'])}): Passive switch at once at drug -> water (group midpoint "
        f"{Pv['sig_obs']['groupcurve']['deg']['mid']:.1f}) and more slowly at water -> drug (licking fell further on day 18).",
        "c: start of the drug -> water change in the group curve vs 2000 data sets simulated from each fitted model and "
        "fitted the same way: " + "; ".join(
            f"{LAB[k]} median {np.percentile(A['A_models'][k]['sig_sim_groupcurve']['deg_start'], 50):.1f} (95% range "
            f"{np.percentile(A['A_models'][k]['sig_sim_groupcurve']['deg_start'], 2.5):.1f}-"
            f"{np.percentile(A['A_models'][k]['sig_sim_groupcurve']['deg_start'], 97.5):.1f}), p {fp(C[k]['ppc_groupcurve_start'])}"
            for k in MAIN) + ". Midpoint, same test: " + "; ".join(f"{LAB[k]} p {fp(C[k]['ppc_groupcurve_mid'])}" for k in MAIN) + ".",
        f"Paper-style comparison (noise-free model curve vs each mouse): models start by trial "
        f"{max(C[k]['start_vs_model']['mean_b'] for k in MAIN):.1f}, mice {C[MAIN[0]]['start_vs_model']['mean_a']:.1f} "
        f"(paired t, Holm p {fp(C[MAIN[0]]['start_vs_model']['p_t_holm'])} for every model). This is biased by noise: "
        "fitting noisy yes/no licking already moves the start later, which is why c uses simulated data.",
        "MEANING: the asymmetry of slide 5 is visible trial by trial (late after morphine loss, fast after its return), and "
        "the start is too late for an instant switch. These data cannot separate meta-RPE from two fixed rates or from a "
        "delayed drug effect (n = 6, one switch each way). Water -> drug rose within ~2.5 min, before an oral dose acts; "
        "the taste may act as a learned drug cue (not tested)."])


def notes8_kr():
    A, Pv = S["Active"], S["Passive"]
    gc, gb = A["sig_obs"]["groupcurve"], A["groupcurve_boot"]
    tA, tP = A["sig_tests"]["deg_vs_enh_mid_all"], Pv["sig_tests"]["deg_vs_enh_mid_all"]
    pmA = PM[PM.group == "Active"]
    C = A["C_start_vs_models"]
    return "\n".join([
        "[이 슬라이드] 논문 Fig. 1d,e,m을 우리 trial에 그대로 적용했다. d: 전환부터 trial별 licking에 sigmoid 피팅; e: 그 "
        "midpoint (논문: H→L이 L→H보다 늦음, P = 0.027, n = 20); m: 변화가 시작되는(10%) trial을 각 모델과 비교 (논문: "
        "meta-RPE만 마우스와 맞음).",
        f"[a 읽는 법] Active 그룹 평균 곡선, 0 = 전환 전, 1 = 전환 후가 되게 정규화. 약 → 물 midpoint trial {gc['deg']['mid']:.1f} "
        f"(95% CI {gb['deg_mid_CI'][0]:.1f}–{gb['deg_mid_CI'][2]:.1f}; trial 30 = day 15 시작), 시작 {gc['deg']['start']:.1f}; "
        f"물 → 약 midpoint {gc['enh']['mid']:.1f} (CI {gb['enh_mid_CI'][0]:.1f}–{gb['enh_mid_CI'][2]:.1f}). 마우스 bootstrap "
        f"2000번 중 {gb['frac_deg_mid_later'] * 100:.1f}%에서 약 → 물이 더 늦다.",
        f"[b 읽는 법] 논문처럼 마우스별. Active 약 → 물 {', '.join(f'{v:.1f}' for v in pmA.deg_mid)}; 물 → 약 "
        f"{', '.join(f'{v:.1f}' for v in pmA.enh_mid)} (평균 {tA['mean_a']:.1f} vs {tA['mean_b']:.1f}, 대응 t p {fp(tA['p_t'])}, "
        f"Wilcoxon p {fp(tA['p_wilcoxon'])}, {tA['n']}마리 중 {tA['n_a_greater']}마리가 더 늦음). Passive {tP['mean_a']:.1f} vs "
        f"{tP['mean_b']:.1f} (p {fp(tP['p_t'])}): Passive는 약 → 물에서 바로 바뀌고 (그룹 midpoint "
        f"{Pv['sig_obs']['groupcurve']['deg']['mid']:.1f}), 물 → 약에서는 더 천천히 바뀐다 (day 18에 licking이 더 줄었다).",
        "[c 읽는 법] 그룹 곡선의 약 → 물 시작점을, 피팅된 각 모델로 만든 2000세트에 같은 방법을 적용한 값과 비교: " + "; ".join(
            f"{LAB[k]} 중앙값 {np.percentile(A['A_models'][k]['sig_sim_groupcurve']['deg_start'], 50):.1f} (95% 범위 "
            f"{np.percentile(A['A_models'][k]['sig_sim_groupcurve']['deg_start'], 2.5):.1f}–"
            f"{np.percentile(A['A_models'][k]['sig_sim_groupcurve']['deg_start'], 97.5):.1f}), p {fp(C[k]['ppc_groupcurve_start'])}"
            for k in MAIN) + ". midpoint로 같은 검정: " + "; ".join(f"{LAB[k]} p {fp(C[k]['ppc_groupcurve_mid'])}" for k in MAIN) + ".",
        f"[논문 방식 비교] 잡음 없는 모델 곡선 vs 마우스별 값: 모델은 모두 trial "
        f"{max(C[k]['start_vs_model']['mean_b'] for k in MAIN):.1f} 이내에 시작, 마우스는 {C[MAIN[0]]['start_vs_model']['mean_a']:.1f} "
        f"(대응 t, 모든 모델 Holm p {fp(C[MAIN[0]]['start_vs_model']['p_t_holm'])}). 하지만 이 비교는 잡음 때문에 치우친다: 잡음 "
        "있는 yes/no licking에 피팅하면 시작점이 원래 늦게 나온다. 그래서 c는 시뮬레이션 데이터와 비교했다.",
        "[의미] 슬라이드 5의 비대칭이 trial 단위로도 보인다 (모르핀이 빠지면 늦게, 돌아오면 빠르게). 시작점은 즉시 전환으로는 "
        "너무 늦다. 하지만 meta-RPE를 고정 학습률 두 개나 늦게 오는 약효와 구분하지는 못한다 (n = 6, 전환 각 방향 한 번). "
        "물 → 약은 약 2.5분 안에 올랐는데, 마신 모르핀이 작용하기 전이다. 맛이 학습된 약물 큐로 작용했을 수 있다 (검정 안 함)."])


# ============================================================================================ main

INDEX_PARTS = {"[Content_Types].xml", "ppt/presentation.xml", "ppt/_rels/presentation.xml.rels"}


def restore_v2_parts(v2, v3):
    """python-pptx re-serialises every slide it loads and drops whitespace-only runs (v2 slide 2 has one), so put the
    original v2 bytes back for every v2 part except the three index parts that must list the new slides."""
    import zipfile
    with zipfile.ZipFile(v2) as za:
        orig = {n: za.read(n) for n in za.namelist()}
    with zipfile.ZipFile(v3) as zb:
        items = [(i, zb.read(i.filename)) for i in zb.infolist()]
    tmp = v3.with_suffix(".tmp")
    with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zo:
        for info, data in items:
            n = info.filename
            zo.writestr(n, orig[n] if (n in orig and n not in INDEX_PARTS) else data)
    tmp.replace(v3)
    with zipfile.ZipFile(v3) as zb:
        names = set(zb.namelist())
        same = all(zb.read(n) == orig[n] for n in orig if n not in INDEX_PARTS)
    return same and set(orig) <= names


def build_into(prs, F7, F8):
    slide7(prs, F7)
    slide8(prs, F8)


def main():
    F7, F8 = fig7(), fig8()
    outs, v3 = {}, {}
    for lang in ("EN", "KR"):
        b.LANG = lang
        two = Presentation(str(b.TEMPLATE))
        build_into(two, F7, F8)
        outs[lang] = OUT / f"mrpe_trial_by_trial_2slides_{lang}_20261118.pptx"
        two.save(outs[lang])
        for fn, F, nm in ((slide7, F7, "7_day14_vs_overnight"), (slide8, F8, "8_switch_timing")):
            one = Presentation(str(b.TEMPLATE))
            fn(one, F)
            one.save(OUT / "single_slides" / f"{lang}_slide{nm}_20261118.pptx")
        deck = Presentation(str(V2 / f"mrpe_hjort2026_opioid_6slides_{lang}_v2_20261117.pptx"))
        n_before = len(deck.slides)
        build_into(deck, F7, F8)
        check(f"[{lang}] v3 deck = 6 v2 slides + 2 new", n_before == 6 and len(deck.slides) == 8)
        v3[lang] = OUT / f"mrpe_hjort2026_opioid_8slides_{lang}_v3_20261118.pptx"
        deck.save(v3[lang])
        v2p = V2 / f"mrpe_hjort2026_opioid_6slides_{lang}_v2_20261117.pptx"
        check(f"[{lang}] v3 deck: every v2 part except the 3 index parts is byte-identical",
              restore_v2_parts(v2p, v3[lang]))
    for k, (txt, val) in b.PRINTED.items():
        t_ = txt.replace("−", "-").replace("+", "").strip()
        if t_.startswith("<"):
            ok = val < float(t_[1:])
        else:
            nd = len(t_.split(".")[1]) if "." in t_ else 0
            ok = abs(float(t_) - val) <= 0.5 * 10 ** -nd + 1e-12
        check(f"printed {k} = {txt} matches {val:.6g}", ok)
    json.dump(dict(checks=b.CHECKS, fit=b.FIT), open(OUT / "fit_report.json", "w", encoding="utf-8"), indent=1,
              ensure_ascii=False)
    bad = [c_ for c_ in b.CHECKS if not c_[1]]
    print(f"{len(b.CHECKS)} checks, {len(bad)} failed")
    for c_ in bad:
        print("  FAILED", c_[0])
    for lang in ("EN", "KR"):
        print(lang, [p_.name for p_ in b.render(outs[lang], OUT / f"preview_{lang}")])
        print(lang, "v3", [p_.name for p_ in b.render(v3[lang], OUT / f"preview_v3_{lang}")])


if __name__ == "__main__":
    main()
