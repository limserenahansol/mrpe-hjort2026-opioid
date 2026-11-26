"""Step 17 (26 Nov 2026): the whole meta-RPE analysis on 2 summary slides, English (EN) and Korean (KR) decks.

User request (26 Nov 2026): "make 2 slide summarizing all results". python-pptx only (never Office/COM).
Slide 1: behaviour at the two switches (m11) + group model fits and held-out CV (m05, m11) + the two learning rates (m08).
Slide 2: what holds / what does not / what this design cannot test (m03, m05, m11, m13, m15; the claims table of
         teaching_package_20261119/docs/06) + power for a follow-up cohort (statsmodels, from the between-mouse SD in m11;
         same numbers as deep-research-output/dopamine-rpe-drug-cessation/hypotheses/power_results.json, asserted).
Every number printed in slide text is registered (reg) and checked against the data at the end; text boxes are measured
with the real fonts (helpers from m10). Phase names: Equal Access (days 11-13), Withdrawal (14-16), Re-exposure (17-18).
Earlier decks are not touched.
Output (slides_20261126_summary\\): mrpe_summary_2slides_{EN,KR}_20261126.pptx, figs\\, preview_{EN,KR}\\,
        fit_report.json
"""
import json
import sys
from pathlib import Path

import matplotlib.transforms as mtrans
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.power import TTestIndPower, TTestPower

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m10_make_slides as b  # noqa: E402  (fonts, text boxes, fit checks, render)

plt, Line2D, Presentation, Patch = b.plt, b.Line2D, b.Presentation, b.Patch
P, Rn, tbox, measure, common = b.P, b.Rn, b.tbox, b.measure, b.common
sf, fp, peq, reg, check, T = b.sf, b.fp, b.peq, b.reg, b.check, b.T
D = b.D
OUT = HERE / "slides_20261126_summary"
FIGD = OUT / "figs"
FIGD.mkdir(parents=True, exist_ok=True)
RED, BLUE, GREEN, ORANGE, PURPLE, GREY = b.RED, b.BLUE, b.GREEN, b.ORANGE, b.PURPLE, b.GREY
CMOD = b.CMOD

MF = json.load(open(D / "m11_modelfree.json"))
G5 = json.load(open(D / "m05_group_fits.json"))
LO = json.load(open(D / "m11_lomo_cv.json"))
LOR = pd.read_csv(D / "m11_lomo_cv_rows.csv")
BS = json.load(open(D / "m08_bootstrap_asymmetry.json"))
REC = json.load(open(D / "m05_recovery.json"))
GR = json.load(open(D / "m11_group_recovery.json"))
S3 = json.load(open(D / "m03_paper_data_summary.json"))
S13 = json.load(open(D / "m13_summary.json"))
S15 = json.load(open(D / "m15_summary.json"))
BINS = pd.read_csv(D / "m13_bins.csv")
PM = BINS.groupby(["group", "mouse", "day"])["lick"].mean().unstack("day")
NA, NP = len(PM.loc["Active"]), len(PM.loc["Passive"])
POWER_JSON = Path(r"C:\Users\hsollim\deep-research-output\dopamine-rpe-drug-cessation\hypotheses\power_results.json")
SW = [("day14_vs_13", 14, 13), ("day15_vs_13", 15, 13), ("day16_vs_13", 16, 13), ("day17_vs_16", 17, 16),
      ("day18_vs_16", 18, 16)]
TO = "today's outcome"
MODELS = ["RPE2a", "mRPE", "RPE", TO]


def TL(en, kr):
    """T() for a list of strings"""
    return en if b.LANG == "EN" else [T(e, k) for e, k in zip(en, kr)]


def common2(s, title, did, found, foot, notes):
    """b.common with the same boxes; KR bottom line 13.5 pt and footer 9.5 pt so that 3 lines fit"""
    ts = 28 if b.LANG == "EN" else 26
    tbox(s, 0.35, 0.10, 12.63, 0.80, [P(Rn(title, ts, True, b.NAVY), align=b.PP_ALIGN.CENTER)], "Title",
         anchor=b.MSO_ANCHOR.MIDDLE)
    bs = 16 if b.LANG == "EN" else 13.5
    tbox(s, 0.45, 5.90, 12.45, 0.92, [P(Rn(T("Did: ", "한 것: "), bs, True), Rn(did, bs),
                                         Rn(T("  Found: ", "  찾은 것: "), bs, True), Rn(found, bs))], "Bottom line")
    fs = 11 if b.LANG == "EN" else 9.5
    tbox(s, 0.45, 6.83, 12.45, 0.665, [P(Rn(foot, fs, color="555555"))], "Footer (n and test)")
    s.notes_slide.notes_text_frame.text = notes

# ------------------------------------------------------------------ per-mouse changes (must equal m11 exactly)
CH = {}
for key, d_, ref in SW:
    for g in ("Active", "Passive"):
        x = (PM.loc[g][d_] - PM.loc[g][ref]).dropna().values
        CH[(key, g)] = x
        m = MF[key][g]["whole"]
        half = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
        check(f"{key} {g}: per-mouse mean, CI and one-sample p equal m11",
              abs(x.mean() - m["mean"]) < 1e-9 and abs(x.mean() - half - m["lo"]) < 1e-9
              and abs(stats.ttest_1samp(x, 0).pvalue - m["p"]) < 1e-9 and len(x) == m["n"])
    w = stats.ttest_ind(CH[(key, "Active")], CH[(key, "Passive")], equal_var=False).pvalue
    check(f"{key}: Welch p equals m11", abs(w - MF[key]["Active_vs_Passive_whole"]["welch_p"]) < 1e-9)
check("n = 6 Active, 7 Passive", NA == 6 and NP == 7)

# ------------------------------------------------------------------ held-out CV: mice where meta-RPE beats RPE 2a
_r = LOR[LOR.group == "Active"].pivot(index="held", columns="model", values="test_negLL")
N_MRPE_BETTER = int((_r["mRPE"] < _r["RPE2a"]).sum())
check("held-out: RPE 2a total negLL is the best, meta-RPE second", LO["Active"]["best"] == "RPE2a"
      and sorted(LO["Active"]["delta_from_best"].items(), key=lambda kv: kv[1])[1][0] == "mRPE")

# ------------------------------------------------------------------ power (same inputs and code path as power_analysis.py)
A_ = {k: MF[k]["Active"]["whole"] for k in ("day14_vs_13", "day15_vs_13", "day16_vs_13", "day17_vs_16")}
sdA = {k: (v["hi"] - v["lo"]) / 2 * np.sqrt(NA) / stats.t.ppf(0.975, NA - 1) for k, v in A_.items()}
SD = float(np.mean(list(sdA.values())))
LATE = (A_["day15_vs_13"]["mean"] + A_["day16_vs_13"]["mean"]) / 2
DEFF = (abs(LATE) - abs(A_["day14_vs_13"]["mean"])) / SD
IND, ONE = TTestIndPower(), TTestPower()
N80 = IND.solve_power(DEFF, alpha=0.05, power=0.8)
N90 = IND.solve_power(DEFF, alpha=0.05, power=0.9)
PW6 = IND.power(DEFF, nobs1=6, alpha=0.05)
NHALF = IND.solve_power(DEFF / 2, alpha=0.05, power=0.8)
RS = (0.3, 0.5, 0.7)
DZ = {r: DEFF / np.sqrt(2 * (1 - r)) for r in RS}
NW = {r: ONE.solve_power(DZ[r], alpha=0.05, power=0.8) for r in RS}
PW6W = {r: ONE.power(DZ[r], nobs=6, alpha=0.05) for r in RS}
_rng = np.random.default_rng(1)


def _sim(n, mu1, mu2, s1, s2, reps=20000):
    a = _rng.normal(mu1, s1, (reps, n))
    c = _rng.normal(mu2, s2, (reps, n))
    return float((stats.ttest_ind(a, c, axis=1, equal_var=False).pvalue < 0.05).mean())


SIMN = (4, 6, 8, 10, 12, 16)
SIM = {n: _sim(n, A_["day14_vs_13"]["mean"], LATE, sdA["day14_vs_13"], np.mean([sdA["day15_vs_13"], sdA["day16_vs_13"]]))
       for n in SIMN}
if POWER_JSON.exists():
    PJ = json.load(open(POWER_JSON))["plans"]
    check("power: d, n80, n90, power at 6 equal power_results.json",
          all(abs(u - v) < 1e-9 for u, v in ((DEFF, PJ["E1_FR1_vs_PR_two_sample"]["d"]),
                                             (N80, PJ["E1_FR1_vs_PR_two_sample"]["n_per_group_80"]),
                                             (N90, PJ["E1_FR1_vs_PR_two_sample"]["n_per_group_90"]),
                                             (PW6, PJ["E1_FR1_vs_PR_two_sample"]["power_at_n6"]),
                                             (NHALF, PJ["interaction_group_x_manipulation"]["n_per_cell_80"]))))
    check("power: Welch simulation equals power_results.json",
          all(abs(SIM[n] - PJ["E1_simulated_welch"][f"n={n}"]) < 1e-12 for n in SIMN))
print(f"SD {SD:.4f}  d {DEFF:.4f}  n80 {N80:.2f}  n90 {N90:.2f}  power@6 {PW6:.3f}  half {NHALF:.1f}")
print("within n80:", {r: round(v, 2) for r, v in NW.items()}, " power@6:", {r: round(v, 3) for r, v in PW6W.items()})
print("sim:", SIM, " meta-RPE better held-out:", N_MRPE_BETTER)


def save(fig, name):
    p_ = FIGD / name
    fig.savefig(p_, facecolor="white")
    plt.close(fig)
    return p_


# ============================================================================================ figure, slide 1
def fig_s1():
    W, H = 12.70, 4.85
    fig = plt.figure(figsize=(W, H))

    def ax_in(x0, y0, w, h):
        return fig.add_axes([x0 / W, y0 / H, w / W, h / H])

    a = ax_in(1.30, 1.40, 4.10, 2.85)
    # ---- a: per-mouse change at each switch
    a.axvspan(-0.5, 2.5, color="#E8F1FB", zorder=0, lw=0)
    a.axvspan(2.5, 4.5, color="#F3ECF8", zorder=0, lw=0)
    for i, (key, _, _) in enumerate(SW):
        for g, off, col in (("Active", -0.19, RED), ("Passive", 0.19, BLUE)):
            x = np.sort(CH[(key, g)])
            jit = np.linspace(-0.07, 0.07, len(x))
            a.scatter(i + off + jit, x, s=16, color=col, alpha=0.40, lw=0, zorder=3)
            m = MF[key][g]["whole"]
            a.errorbar(i + off, m["mean"], yerr=[[m["mean"] - m["lo"]], [m["hi"] - m["mean"]]], fmt="o", color=col,
                       ms=6.5, capsize=3.5, lw=1.9, zorder=4)
    a.axhline(0, color="#444444", lw=0.9, zorder=2)
    a.text(1.0, 1.04, "Withdrawal (water)\nvs day 13 (Equal Access)", ha="center", va="top", fontsize=10.5,
           color="#1F4E79")
    a.text(3.5, 1.04, "Re-exposure (morphine)\nvs day 16", ha="center", va="top", fontsize=10.5, color="#6A3D9A")
    a.set_xlim(-0.5, 4.5)
    a.set_ylim(-1.28, 1.08)
    a.set_xticks(range(5), ["day 14", "day 15", "day 16", "day 17", "day 18"])
    a.set_ylabel("Change in p(lick) vs reference day")
    hs = [Line2D([], [], marker="o", color=RED, lw=1.9, ms=6.5, label=f"Active (n = {NA})"),
          Line2D([], [], marker="o", color=BLUE, lw=1.9, ms=6.5, label=f"Passive (n = {NP})"),
          Line2D([], [], marker="o", color="#999999", lw=0, ms=4.5, alpha=0.6, label="dot = mouse; bar = mean ± 95% CI")]
    a.legend(handles=hs, loc="lower center", ncol=2, frameon=True, framealpha=0.95, edgecolor="#DDDDDD", fontsize=9.5,
             borderaxespad=0.25, handlelength=1.5, columnspacing=1.2)
    a.set_yticks([-1.0, -0.5, 0.0, 0.5, 1.0])
    tr = mtrans.blended_transform_factory(a.transData, a.transAxes)
    rows = [("Active vs 0", [MF[k]["Active"]["whole"]["p"] for k, _, _ in SW], RED),
            ("Passive vs 0", [MF[k]["Passive"]["whole"]["p"] for k, _, _ in SW], BLUE),
            ("Active vs Passive", [MF[k]["Active_vs_Passive_whole"]["welch_p"] for k, _, _ in SW], "#333333")]
    for j, (lab, ps, col) in enumerate(rows):
        yy = -0.175 - 0.088 * j
        a.text(-0.62, yy, lab, transform=tr, ha="right", va="center", fontsize=10, color=col)
        for i, p_ in enumerate(ps):
            a.text(i, yy, fp(p_), transform=tr, ha="center", va="center", fontsize=10, color=col,
                   fontweight="bold" if p_ < 0.05 else "normal")
    fig.text(0.20 / W, 4.62 / H, "a  Change in licking at each switch", fontsize=13.5, fontweight="bold", va="center")
    fig.text(0.20 / W, 0.17 / H, "p: one-sample t vs 0, Welch t (bold: p < 0.05)", fontsize=9.5, color="#555555",
             va="center")

    # ---- b: group fits and held-out CV, Active
    b1 = ax_in(6.80, 1.40, 1.15, 2.85)
    b2 = ax_in(8.25, 1.40, 1.15, 2.85)
    names = ["RPE 2α\n(two rates)", "meta-RPE", "RPE\n(one rate)", "instant\nswitch"]
    yv = np.arange(4)[::-1]
    dbic = [G5["Active"][m]["dBIC_from_best"] for m in MODELS]
    dho = [LO["Active"]["delta_from_best"][m] for m in MODELS]
    for ax, vals, lim, ttl, xl in ((b1, dbic, 30, "group fit", "ΔBIC"), (b2, dho, 21, "held-out mice", "Δ held-out negLL")):
        ax.barh(yv, vals, height=0.62, color=[CMOD[m] for m in MODELS])
        for y_, v in zip(yv, vals):
            ax.text(v + lim * 0.03, y_, "0 best" if v == 0 else f"+{v:.1f}", va="center", fontsize=10.5)
        ax.set_xlim(0, lim)
        ax.set_ylim(-0.6, 3.6)
        ax.set_title(ttl, fontsize=11.5)
        ax.set_xlabel(xl, fontsize=11)
        ax.tick_params(axis="x", labelsize=10.5)
    b1.set_yticks(yv, names, fontsize=11)
    b2.set_yticks(yv, [])
    fig.text(5.72 / W, 4.62 / H, f"b  Model fits, Active (n = {NA})", fontsize=13.5, fontweight="bold", va="center")
    fig.text(5.72 / W, 0.42 / H, f"Lower = better. Held-out: RPE 2α and meta-RPE tie\n(meta-RPE better in {N_MRPE_BETTER} "
             f"of {NA} mice)", fontsize=10, color="#555555", va="center")

    # ---- c: RPE 2a learning rates, 95% bootstrap CI
    c = ax_in(10.60, 1.40, 1.95, 2.85)
    rows_c = [("Active α−", "Active", "a_minus", RED, "v", 3.3), ("Active α+", "Active", "a_plus", RED, "^", 2.4),
              ("Passive α−", "Passive", "a_minus", BLUE, "v", 1.1), ("Passive α+", "Passive", "a_plus", BLUE, "^", 0.2)]
    for lab, g, k, col, mk, y_ in rows_c:
        v, ci = BS[g]["full_fit"][k], BS[g][f"{k}_CI"]
        c.plot([ci[0], ci[2]], [y_, y_], color=col, lw=2.0)
        c.plot([ci[0], ci[2]], [y_, y_], "|", color=col, ms=8, mew=1.8)
        c.plot(v, y_, marker=mk, color=col, ms=9, zorder=4)
        c.text(0.0105, y_ + 0.30, f"{v:.2f} [{ci[0]:.2f}, {ci[2]:.2f}]", fontsize=10, color=col, va="bottom")
    c.text(0.0105, 1.78, f"α− < α+: {BS['Active']['frac_a_minus_less_than_a_plus'] * 100:.1f}%", fontsize=10,
           color=RED, fontweight="bold")
    c.text(0.0105, -0.45, f"α− < α+: {BS['Passive']['frac_a_minus_less_than_a_plus'] * 100:.1f}%", fontsize=10,
           color=BLUE)
    c.set_xscale("log")
    c.set_xlim(0.0085, 1.4)
    c.set_ylim(-0.7, 4.1)
    c.set_yticks([r[5] for r in rows_c], [r[0] for r in rows_c], fontsize=11)
    c.set_xlabel("Learning rate (log), 95% CI", fontsize=11)
    c.tick_params(axis="x", labelsize=10.5)
    fig.text(9.55 / W, 4.62 / H, "c  RPE 2α learning rates", fontsize=13.5, fontweight="bold", va="center")
    fig.text(9.55 / W, 0.42 / H, "α−: morphine → water (negative error)\nα+: water → morphine (positive error)",
             fontsize=10, color="#555555", va="center")
    return save(fig, "s1_summary_results.png"), (W, H)


# ============================================================================================ figure, slide 2 (power)
def fig_s2():
    W, H = 5.30, 4.85
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0.82 / W, 1.00 / H, 4.33 / W, 3.25 / H])
    ns = np.linspace(2, 50, 400)
    lo_ = [ONE.power(DZ[0.3], nobs=n, alpha=0.05) for n in ns]
    hi_ = [ONE.power(DZ[0.7], nobs=n, alpha=0.05) for n in ns]
    ax.fill_between(ns, lo_, hi_, color=GREEN, alpha=0.13, lw=0)
    ax.plot(ns, [ONE.power(DZ[0.5], nobs=n, alpha=0.05) for n in ns], color=GREEN, lw=2.2,
            label=f"within mouse, d = {DEFF:.2f} (r = 0.5; band 0.3–0.7)")
    ax.plot(ns, [IND.power(DEFF, nobs1=n, alpha=0.05) for n in ns], color="#222222", lw=2.2,
            label=f"between groups, d = {DEFF:.2f}")
    ax.plot(ns, [IND.power(DEFF / 2, nobs1=n, alpha=0.05) for n in ns], color=ORANGE, lw=2.2,
            label=f"half effect or interaction, d = {DEFF / 2:.2f}")
    ax.plot(list(SIM), list(SIM.values()), "o", mfc="white", mec="#222222", mew=1.4, ms=6.5,
            label="Welch t simulation (20,000 draws)")
    ax.axhline(0.8, color="#888888", ls="--", lw=1)
    ax.axvline(6, color="#888888", ls=":", lw=1.2)
    ax.text(6, 1.03, "current n = 6", ha="center", va="bottom", fontsize=9.5, color="#555555", clip_on=False)
    ax.plot(6, PW6, "s", color="#222222", ms=5, zorder=5)
    ax.text(6.8, 0.42, f"power {PW6:.2f}\nat n = 6", fontsize=9.5, color="#222222", va="top")
    for n, col, txt, tx, ty, ha in ((NW[0.5], GREEN, f"{int(np.ceil(NW[0.5]))} mice", 2.6, 0.935, "left"),
                                    (N80, "#222222", f"{int(np.ceil(N80))} per group", N80 + 0.9, 0.725, "left"),
                                    (NHALF, ORANGE, f"{int(np.ceil(NHALF))} per group", NHALF - 1.0, 0.845, "right")):
        ax.plot(n, 0.8, "D", color=col, ms=6.5, zorder=5)
        ax.text(tx, ty, txt, color=col, fontsize=10.5, fontweight="bold", ha=ha, va="center")
    ax.set_xlim(2, 50)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("Mice per group (within mouse: mice in total)")
    ax.set_ylabel("Power (two-sided t, α = 0.05)")
    ax.legend(loc="lower right", frameon=True, framealpha=0.95, edgecolor="#DDDDDD", fontsize=9.3, borderaxespad=0.3)
    fig.text(0.12 / W, 4.62 / H, "b  Mice needed to confirm the day-14 lag", fontsize=13.5, fontweight="bold",
             va="center")
    return save(fig, "s2_power.png"), (W, H)


# ============================================================================================ slide 1
def slide1(prs, F):
    b.SLIDE = 1
    path, size = F
    s = b.new_slide(prs)
    b.place(s, path, size, 0.32, 0.98, 12.70)
    m14a, m14p = MF["day14_vs_13"]["Active"]["whole"], MF["day14_vs_13"]["Passive"]["whole"]
    d14a = reg("z1_d14a", sf(m14a["mean"]), m14a["mean"])
    p14a = reg("z1_p14a", fp(m14a["p"]), m14a["p"])
    d14p = reg("z1_d14p", sf(m14p["mean"]), m14p["mean"])
    w14 = reg("z1_w14", fp(MF["day14_vs_13"]["Active_vs_Passive_whole"]["welch_p"]),
              MF["day14_vs_13"]["Active_vs_Passive_whole"]["welch_p"])
    am = reg("z1_am", f"{BS['Active']['full_fit']['a_minus']:.2f}", BS["Active"]["full_fit"]["a_minus"])
    ap = reg("z1_ap", f"{BS['Active']['full_fit']['a_plus']:.2f}", BS["Active"]["full_fit"]["a_plus"])
    check("slide 1: RPE 2a group fit rates equal m05 (m08 stores them rounded to 4 decimals)",
          abs(G5["Active"]["RPE2a"]["params"][0] - BS["Active"]["full_fit"]["a_minus"]) < 5e-5
          and abs(G5["Active"]["RPE2a"]["params"][1] - BS["Active"]["full_fit"]["a_plus"]) < 5e-5)
    check("slide 1: RPE 2a best by group BIC and by held-out negLL", G5["Active"]["RPE2a"]["dBIC_from_best"] == 0
          and LO["Active"]["delta_from_best"]["RPE2a"] == 0)
    common2(s, T("Active licking lags after morphine removal; Passive does not",
                "모르핀 제거 후 Active는 하루 늦게 바뀌고, Passive는 바로 바뀐다"),
           T("fitted the Hjort 2026 learning models (our port reproduces the paper's ranking) to PR licking on days "
             "11–18; tested every effect at the mouse level.",
             "Hjort 2026 학습 모델을 day 11–18 PR licking에 맞추고 모든 효과를 마우스 단위로 검정했다 (포트는 논문의 "
             "순위를 재현)."),
           T(f"Active licking held on the first water day ({d14a}, p = {p14a}; Passive {d14p}, Welch p = {w14}) and "
             f"fell on days 15–16; two learning rates fit Active best (α− {am} vs α+ {ap}).",
             f"Active는 첫 물 날 licking이 그대로였고 ({d14a}, p = {p14a}; Passive {d14p}, Welch p = {w14}) day 15–16에 "
             f"줄었다. 학습률 두 개가 Active에 가장 잘 맞았다 (α− {am} vs α+ {ap})."),
           T(f"n = {NA} Active, {NP} Passive (6099 white died on day 13). p(lick) = share of 30-s bins with ≥ 1 lick; "
             "change = per-mouse session mean minus the reference day (13 or 16). a: one-sample t vs 0 and Welch t, "
             "two-sided. b: group-shared learning rates + each mouse's readout; held-out = leave-one-mouse-out CV. "
             "c: 1000 bootstraps over mice; bounds 0.01–1. A delayed oral-drug effect, Withdrawal or taste can also "
             "slow stopping.",
             f"n = Active {NA}, Passive {NP} (6099 white는 day 13에 사망). p(lick) = lick이 1번 이상 있는 30초 bin의 비율; "
             "변화 = 마우스별 세션 평균 − 기준일 (13 또는 16). a: 0 대비 one-sample t, Welch t, 양측. b: 그룹이 공유하는 "
             "학습률 + 마우스별 readout; held-out = 한 마리씩 빼고 교차검증. c: 마우스 bootstrap 1000번; 범위 0.01–1. "
             "늦게 작용하는 경구 약물 효과, Withdrawal, 맛도 멈춤을 늦출 수 있다."),
           T(notes1_en(), notes1_kr()))
    return s


def _n1():
    A, G, L = MF, G5["Active"], LO["Active"]["delta_from_best"]
    return dict(
        d14a=sf(A["day14_vs_13"]["Active"]["whole"]["mean"]), p14a=fp(A["day14_vs_13"]["Active"]["whole"]["p"]),
        d15a=sf(A["day15_vs_13"]["Active"]["whole"]["mean"]), d16a=sf(A["day16_vs_13"]["Active"]["whole"]["mean"]),
        d14p=sf(A["day14_vs_13"]["Passive"]["whole"]["mean"]),
        d17a=sf(A["day17_vs_16"]["Active"]["whole"]["mean"]), p17a=fp(A["day17_vs_16"]["Active"]["whole"]["p"]),
        bm=f"{G['mRPE']['dBIC_from_best']:.1f}", br=f"{G['RPE']['dBIC_from_best']:.1f}",
        bi=f"{G[TO]['dBIC_from_best']:.1f}",
        hm=f"{L['mRPE']:.1f}", hr=f"{L['RPE']:.1f}", hi=f"{L[TO]:.1f}",
        k=N_MRPE_BETTER,
        am=f"{BS['Active']['full_fit']['a_minus']:.2f}", al=f"{BS['Active']['a_minus_CI'][0]:.2f}",
        ah=f"{BS['Active']['a_minus_CI'][2]:.2f}", ap=f"{BS['Active']['full_fit']['a_plus']:.2f}",
        pl=f"{BS['Active']['a_plus_CI'][0]:.2f}", ph=f"{BS['Active']['a_plus_CI'][2]:.2f}",
        fa=f"{BS['Active']['frac_a_minus_less_than_a_plus'] * 100:.1f}",
        fpv=f"{BS['Passive']['frac_a_minus_less_than_a_plus'] * 100:.1f}",
        m3=f"{S3['BIC']['mean_change_from_best']['mRPE']:.1f}",
        m3lo=f"{min(v for k_, v in S3['BIC']['mean_change_from_best'].items() if k_ != 'mRPE'):.1f}",
        m3hi=f"{max(S3['BIC']['mean_change_from_best'].values()):.1f}")


def notes1_en():
    v = _n1()
    return (
        "WHAT: the behaviour and model results on one slide (details: slides 3-5 and 7-8 of the 8-slide v4 deck).\n"
        "a: per-mouse change in p(lick) (share of 30-s bins with at least one lick) from the last day before each "
        "switch. Days 14-16 (Withdrawal, water) are compared with day 13 (Equal Access, morphine); days 17-18 "
        f"(Re-exposure, morphine) with day 16. Active licking did not change on day 14 ({v['d14a']}, p = {v['p14a']}) "
        f"and fell on days 15-16 ({v['d15a']}, {v['d16a']}); it rose again on day 17 ({v['d17a']}, p = {v['p17a']}). "
        f"Passive licking rose on day 14 already ({v['d14p']}). Rows under the axis: one-sample t vs 0 per group and "
        "Welch t between groups.\n"
        "b: group fits to Active licking on days 11-18 (learning rates shared by the group, each mouse's own readout). "
        f"RPE 2alpha (separate rates after negative and positive errors) is best by BIC (meta-RPE +{v['bm']}, RPE "
        f"+{v['br']}, instant switch +{v['bi']}). On held-out mice (leave-one-mouse-out) RPE 2alpha and meta-RPE tie "
        f"(+{v['hm']}; meta-RPE better in {v['k']} of 6 mice); one rate (+{v['hr']}) and the instant switch "
        f"(+{v['hi']}) are worse.\n"
        f"c: RPE 2alpha rates with 95% CIs from 1000 bootstraps over mice. Active: alpha- {v['am']} [{v['al']}, "
        f"{v['ah']}] vs alpha+ {v['ap']} [{v['pl']}, {v['ph']}]; alpha- < alpha+ in {v['fa']}% of bootstraps. "
        f"Passive rates are wide ({v['fpv']}%).\n"
        "CAVEATS: one switch each way per mouse, n = 6 / 7. Oral morphine acts minutes after a drop, so a delayed drug "
        "effect, the Withdrawal state or taste could produce the same lag. Port check: on the authors' 40 mice our "
        f"Python port gives the paper's ranking (meta-RPE BIC change from best {v['m3']} vs {v['m3lo']}-{v['m3hi']}).")


def notes1_kr():
    v = _n1()
    return (
        "무엇: 행동 결과와 모델 결과를 한 장에 모은 요약입니다 (자세한 내용은 8장짜리 v4 덱의 슬라이드 3-5, 7-8).\n"
        "a: 각 전환 직전 날과 비교한 마우스별 p(lick) 변화입니다 (p(lick) = lick이 1번 이상 있는 30초 bin의 비율). "
        "Day 14-16 (Withdrawal, 물)은 day 13 (Equal Access, 모르핀)과, day 17-18 (Re-exposure, 모르핀)은 day 16과 "
        f"비교했습니다. Active는 day 14에 licking이 변하지 않았고 ({v['d14a']}, p = {v['p14a']}) day 15-16에 줄었으며 "
        f"({v['d15a']}, {v['d16a']}), day 17에 다시 늘었습니다 ({v['d17a']}, p = {v['p17a']}). Passive는 day 14에 이미 "
        f"늘었습니다 ({v['d14p']}). 축 아래 숫자는 그룹별 0 대비 one-sample t와 그룹 간 Welch t의 p입니다.\n"
        "b: day 11-18 Active licking에 대한 그룹 피팅입니다 (학습률은 그룹이 공유하고 readout은 마우스마다 따로). "
        f"음의 오류와 양의 오류에 학습률을 따로 두는 RPE 2alpha가 BIC로 가장 좋습니다 (meta-RPE +{v['bm']}, RPE "
        f"+{v['br']}, 즉시 전환 +{v['bi']}). 한 마리씩 빼고 맞춘 교차검증에서는 RPE 2alpha와 meta-RPE가 비깁니다 "
        f"(+{v['hm']}; 6마리 중 {v['k']}마리에서 meta-RPE가 나음). 학습률 하나짜리 RPE (+{v['hr']})와 즉시 전환 "
        f"(+{v['hi']})은 더 나쁩니다.\n"
        f"c: RPE 2alpha의 두 학습률과 마우스 bootstrap 1000번의 95% CI입니다. Active: alpha- {v['am']} [{v['al']}, "
        f"{v['ah']}] vs alpha+ {v['ap']} [{v['pl']}, {v['ph']}]; bootstrap의 {v['fa']}%에서 alpha- < alpha+. "
        f"Passive는 CI가 넓습니다 ({v['fpv']}%).\n"
        "주의: 마우스마다 방향별 전환이 한 번뿐이고 n = 6 / 7입니다. 경구 모르핀은 drop 뒤 몇 분 지나서 작용하므로, "
        "늦게 나타나는 약물 효과, Withdrawal 상태, 맛도 같은 지연을 만들 수 있습니다. 포트 검증: 저자들의 40마리 "
        f"데이터에서 우리 Python 포트가 논문의 순위를 재현합니다 (meta-RPE BIC 차이 {v['m3']} vs {v['m3lo']}-{v['m3hi']}).")


# ============================================================================================ slide 2
def bullet_box(s, x, y, w, head, items, col, fill, name):
    hs, bs = (12, 10.5) if b.LANG == "EN" else (11.5, 10)
    paras = [P(Rn(head, hs, True, col), after=2)] + [P(Rn("• " + it, bs), after=1.5) for it in items]
    h = measure(paras, w - 0.14) + 0.08
    tb = tbox(s, x, y, w, h, paras, name, fill=fill, line=col)
    for p_ in tb.text_frame.paragraphs[1:]:
        pPr = p_._p.get_or_add_pPr()
        pPr.set("marL", str(int(0.11 * 914400)))
        pPr.set("indent", str(int(-0.11 * 914400)))
    return y + h


def claims():
    A = MF
    m15i = S15["Active"]["models"]["today's outcome"]["A"]["day14"]
    plac = S15["Passive"]["models"]["today's outcome"]["placebo"]["whole"]["mouse_t"]
    asy = MF["asymmetry_whole_session_abs_change_day17_vs_day14"]["Active"]
    gc, gt = S13["Active"]["sig_obs"]["groupcurve"], S13["Active"]["sig_tests"]["deg_vs_enh_mid_all"]
    G = G5["Active"]
    r = dict(
        m3=reg("z2_m3", f"{S3['BIC']['mean_change_from_best']['mRPE']:.1f}", S3["BIC"]["mean_change_from_best"]["mRPE"]),
        m3lo=reg("z2_m3lo", f"{S3['BIC']['mean_change_from_best']['RPE']:.1f}", S3["BIC"]["mean_change_from_best"]["RPE"]),
        m3hi=reg("z2_m3hi", f"{S3['BIC']['mean_change_from_best']['eMack']:.1f}",
                 S3["BIC"]["mean_change_from_best"]["eMack"]),
        d14a=reg("z2_d14a", sf(A["day14_vs_13"]["Active"]["whole"]["mean"]), A["day14_vs_13"]["Active"]["whole"]["mean"]),
        p14a=reg("z2_p14a", fp(A["day14_vs_13"]["Active"]["whole"]["p"]), A["day14_vs_13"]["Active"]["whole"]["p"]),
        d14p=reg("z2_d14p", sf(A["day14_vs_13"]["Passive"]["whole"]["mean"]),
                 A["day14_vs_13"]["Passive"]["whole"]["mean"]),
        w14=reg("z2_w14", fp(A["day14_vs_13"]["Active_vs_Passive_whole"]["welch_p"]),
                A["day14_vs_13"]["Active_vs_Passive_whole"]["welch_p"]),
        h1=reg("z2_h1", sf(S13["Active"]["A_observed"]["d14_h1"]["mean"]), S13["Active"]["A_observed"]["d14_h1"]["mean"]),
        f75=reg("z2_f75", f"{S13['Active']['day14_drops']['frac_by_2p5min'][2] * 100:.0f}",
                S13["Active"]["day14_drops"]["frac_by_2p5min"][2] * 100),
        d15=reg("z2_d15", sf(A["day15_vs_13"]["Active"]["whole"]["mean"]), A["day15_vs_13"]["Active"]["whole"]["mean"]),
        d16=reg("z2_d16", sf(A["day16_vs_13"]["Active"]["whole"]["mean"]), A["day16_vs_13"]["Active"]["whole"]["mean"]),
        p15=reg("z2_p15", fp(A["day15_vs_13"]["Active"]["whole"]["p"]), A["day15_vs_13"]["Active"]["whole"]["p"]),
        p16=reg("z2_p16", fp(A["day16_vs_13"]["Active"]["whole"]["p"]), A["day16_vs_13"]["Active"]["whole"]["p"]),
        d17=reg("z2_d17", sf(A["day17_vs_16"]["Active"]["whole"]["mean"]), A["day17_vs_16"]["Active"]["whole"]["mean"]),
        p17=reg("z2_p17", fp(A["day17_vs_16"]["Active"]["whole"]["p"]), A["day17_vs_16"]["Active"]["whole"]["p"]),
        br=reg("z2_br", f"+{G['RPE']['dBIC_from_best']:.0f}", G["RPE"]["dBIC_from_best"]),
        bi=reg("z2_bi", f"+{G[TO]['dBIC_from_best']:.0f}", G[TO]["dBIC_from_best"]),
        a14=reg("z2_a14", f"{asy['abs_day14']:.2f}", asy["abs_day14"]),
        a17=reg("z2_a17", f"{asy['abs_day17']:.2f}", asy["abs_day17"]),
        pa=reg("z2_pa", fp(asy["p_t"]), asy["p_t"]),
        mid_d=reg("z2_midd", f"{gc['deg']['mid']:.1f}", gc["deg"]["mid"]),
        mid_e=reg("z2_mide", f"{gc['enh']['mid']:.1f}", gc["enh"]["mid"]),
        fl=reg("z2_fl", f"{S13['Active']['groupcurve_boot']['frac_deg_mid_later'] * 100:.2f}",
               S13["Active"]["groupcurve_boot"]["frac_deg_mid_later"] * 100),
        pm=reg("z2_pm", fp(gt["p_t"]), gt["p_t"]),
        pi=reg("z2_pi", fp(m15i["mouse_t"]["p_t"]), m15i["mouse_t"]["p_t"]),
        ph=reg("z2_ph", fp(m15i["mouse_t"]["p_holm"]), m15i["mouse_t"]["p_holm"]),
        nlo=reg("z2_nlo", fp(min(m15i["calibrated"]["p"], m15i["calibrated_B"]["p"])),
                min(m15i["calibrated"]["p"], m15i["calibrated_B"]["p"])),
        nhi=reg("z2_nhi", fp(max(m15i["calibrated"]["p"], m15i["calibrated_B"]["p"])),
                max(m15i["calibrated"]["p"], m15i["calibrated_B"]["p"])),
        k=N_MRPE_BETTER,
        pp=reg("z2_pp", fp(plac["p_t"]), plac["p_t"]),
        r1=reg("z2_r1", f"{REC['confusion_BIC']['Active']['mRPE']['mRPE'] * 100:.0f}",
               REC["confusion_BIC"]["Active"]["mRPE"]["mRPE"] * 100),
        r2=reg("z2_r2", f"{GR['confusion']['Active']['mRPE']['mRPE'] * 100:.0f}",
               GR["confusion"]["Active"]["mRPE"]["mRPE"] * 100),
    )
    check("slide 2: port range = smallest and largest non-meta-RPE BIC change",
          S3["BIC"]["mean_change_from_best"]["RPE"] == min(v for k_, v in S3["BIC"]["mean_change_from_best"].items()
                                                          if k_ != "mRPE")
          and S3["BIC"]["mean_change_from_best"]["eMack"] == max(S3["BIC"]["mean_change_from_best"].values()))
    check("slide 2: two rates never picked for one-rate or instant-switch data (group recovery)",
          GR["confusion"]["Active"]["RPE"]["RPE2a"] == 0 and GR["confusion"]["Active"]["today's outcome"]["RPE2a"] == 0)
    check("slide 2: day-14 first half is the first 7.5 min (15 of 30 bins); drops by 7.5 min = 3rd 2.5-min step",
          len(S13["Active"]["day14_drops"]["frac_by_2p5min"]) == 6)
    return r


def slide2(prs, F):
    b.SLIDE = 2
    path, size = F
    s = b.new_slide(prs)
    b.place(s, path, size, 7.90, 0.98, 5.30)
    r = claims()
    hd = 13.5 if b.LANG == "EN" else 13
    y = 0.98
    tbox(s, 0.30, y, 7.50, 0.37, [P(Rn(T(f"a  What the data show (mouse-level tests, n = {NA} / {NP})",
                                          f"a  데이터가 보여 주는 것 (마우스 단위 검정, n = {NA} / {NP})"),
                                        hd, True))], "Panel a header")
    y += 0.40
    sup = TL([f"Port check: on the paper's 40 mice meta-RPE fits best (ΔBIC {r['m3']} vs {r['m3lo']}–{r['m3hi']}), "
              "as in its Fig. 1k,l",
              f"Day 14: Active licking unchanged ({r['d14a']}, p = {r['p14a']}), Passive rose ({r['d14p']}); "
              f"Welch p = {r['w14']}",
              f"First 7.5 min of day 14: no fall ({r['h1']}) although {r['f75']}% of the day's water drops were earned",
              f"Active fell on days 15–16 ({r['d15']}, {r['d16']}; p = {r['p15']}, {r['p16']}), rose on day 17 "
              f"({r['d17']}, p = {r['p17']})",
              f"Group level: two rates beat one rate and instant switch (ΔBIC {r['br']}, {r['bi']}); 0% false picks"],
             [f"포트 검증: 논문 40마리에서 meta-RPE가 최적 (ΔBIC {r['m3']} vs {r['m3lo']}–{r['m3hi']}), 논문 Fig. 1k,l과 같음",
              f"Day 14: Active licking 변화 없음 ({r['d14a']}, p = {r['p14a']}), Passive는 증가 ({r['d14p']}); Welch "
              f"p = {r['w14']}",
              f"Day 14 처음 7.5분: 그날 물 drop의 {r['f75']}%를 얻었는데도 감소 없음 ({r['h1']})",
              f"Active는 day 15–16에 감소 ({r['d15']}, {r['d16']}; p = {r['p15']}, {r['p16']}), day 17에 증가 "
              f"({r['d17']}, p = {r['p17']})",
              f"그룹 단위: 학습률 두 개가 하나짜리와 즉시 전환보다 나음 (ΔBIC {r['br']}, {r['bi']}); 잘못 고른 비율 0%"])
    nos = TL([f"Slower loss than return, tested directly: first-day change {r['a14']} vs {r['a17']}, paired p = {r['pa']}",
              f"By timing: midpoint {r['mid_d']} vs {r['mid_e']} trials in the group curve, but per mouse p = {r['pm']}",
              f"Instant switch rejected: mouse-level p = {r['pi']} (Holm {r['ph']}); noise-matched {r['nlo']}–{r['nhi']}",
              f"Two rates better than meta-RPE: a tie on held-out mice (meta-RPE better in {r['k']} of {NA})",
              f"Passive above every model: depends on the reference day (placebo pair fails, p = {r['pp']})"],
             [f"제거가 복귀보다 느린지 직접 검정: 첫날 변화 {r['a14']} vs {r['a17']}, 대응 t p = {r['pa']}",
              f"타이밍으로: 그룹 곡선 중간점 {r['mid_d']} vs {r['mid_e']} trial, 그러나 마우스별 p = {r['pm']}",
              f"즉시 전환 기각: 마우스 단위 p = {r['pi']} (Holm {r['ph']}); 잡음 맞춘 시뮬레이션 {r['nlo']}–{r['nhi']}",
              f"학습률 두 개 vs meta-RPE: 남겨 둔 마우스에서 동점 ({NA}마리 중 {r['k']}마리는 meta-RPE가 나음)",
              f"Passive가 모든 모델보다 더 증가: 기준일에 따라 달라짐 (placebo 날 쌍 실패, p = {r['pp']})"])
    unt = TL([f"Error averaging (meta part): meta-RPE data recovered in {r['r1']}% (per mouse) / {r['r2']}% (group) of "
              "simulations",
              "Learning vs a delayed oral-drug effect, Withdrawal state or taste: not separable by behaviour"],
             [f"오류 평균 (meta 부분): meta-RPE로 만든 데이터도 {r['r1']}% (마우스별) / {r['r2']}% (그룹)만 회복",
              "학습 vs 늦게 작용하는 경구 약물 효과, Withdrawal 상태, 맛: 행동만으로 구분 불가"])
    y = bullet_box(s, 0.30, y, 7.50, T("Supported", "지지됨"), sup, "2E7D32", "F1F8F1", "Supported") + 0.07
    y = bullet_box(s, 0.30, y, 7.50, T("Not supported (tested: not significant or not robust)",
                                       "지지되지 않음 (검정했지만 유의하지 않거나 견고하지 않음)"), nos, "B35A00", "FDF3EA",
                   "Not supported") + 0.07
    y = bullet_box(s, 0.30, y, 7.50, T("Untestable here (one switch each way per mouse)",
                                       "이 설계로는 검정 불가 (마우스마다 방향별 전환 1번)"), unt, "5F5F5F", "F4F4F4",
                   "Untestable")
    check(f"[{b.LANG}] slide 2: the three claim boxes end above the bottom line ({y:.2f} in <= 5.86 in)", y <= 5.86)
    sd = reg("z2_sd", f"{SD:.2f}", SD)
    pw = reg("z2_pw6", f"{PW6 * 100:.0f}", PW6 * 100)
    n80 = reg("z2_n80c", f"{int(np.ceil(N80))}", np.ceil(N80))
    nw = reg("z2_nw5c", f"{int(np.ceil(NW[0.5]))}", np.ceil(NW[0.5]))
    d_ = reg("z2_d", f"{DEFF:.2f}", DEFF)
    d14 = reg("z2_pd14", sf(A_["day14_vs_13"]["mean"]), A_["day14_vs_13"]["mean"])
    late = reg("z2_late", sf(LATE), LATE)
    nh = reg("z2_nh", f"{NHALF:.0f}", NHALF)
    nwr = (reg("z2_nw7c", f"{int(np.ceil(NW[0.7]))}", np.ceil(NW[0.7])),
           reg("z2_nw3c", f"{int(np.ceil(NW[0.3]))}", np.ceil(NW[0.3])))
    common2(s, T("What 6 mice can show, and how many a follow-up needs",
                "6마리로 보일 수 있는 것과 다음 실험에 필요한 마우스 수"),
           T(f"sorted every claim by its mouse-level test; computed power from our between-mouse SD ({sd}).",
             f"모든 주장을 마우스 단위 검정 결과로 분류하고, 우리 데이터의 마우스 간 SD ({sd})로 검정력을 계산했다."),
           T(f"the day-14 hold and the two group-level rates hold, but error averaging is untestable and n = 6 gives "
             f"{pw}% power; confirming the lag needs {n80} mice per group, or {nw} if each mouse is tested in both "
             "conditions.",
             f"day 14에 licking이 유지된 것과 그룹 단위의 두 학습률은 지지되지만, 오류 평균은 검정할 수 없고 n = 6의 "
             f"검정력은 {pw}%다. 지연을 확인하려면 그룹당 {n80}마리, 같은 마우스에 두 조건을 모두 주면 {nw}마리가 필요하다."),
           T(f"Effect = Active change on day 14 ({d14}) vs days 15–16 (mean {late}) ÷ between-mouse SD {sd}: d = {d_}. "
             "Two-sample t (statsmodels), two-sided α = 0.05; dots = Welch t simulation. Within mouse: paired t, "
             f"dz = d / √(2(1 − r)); r = 0.7 / 0.3 gives {nwr[0]} / {nwr[1]} mice. Half the effect or a group × "
             f"condition interaction (~4× the n): {nh} per group. No independent-bin p values are used.",
             f"효과 = Active의 day 14 변화 ({d14}) vs day 15–16 (평균 {late}) ÷ 마우스 간 SD {sd}: d = {d_}. 두 표본 t "
             "(statsmodels), 양측 α = 0.05; 점 = Welch t 시뮬레이션. 마우스 내 설계: 대응 t, dz = d / √(2(1 − r)); "
             f"r = 0.7 / 0.3이면 {nwr[0]} / {nwr[1]}마리. 효과가 절반이거나 그룹 × 조건 상호작용 (n 약 4배): 그룹당 "
             f"{nh}마리. 독립 bin p값은 쓰지 않았다."),
           T(notes2_en(r), notes2_kr(r)))
    return s


def _n2():
    return dict(sd=f"{SD:.2f}", d14=sf(A_["day14_vs_13"]["mean"]).replace("\u2212", "-"),
                late=sf(LATE).replace("\u2212", "-"), d=f"{DEFF:.2f}", n80=f"{N80:.1f}", n90=f"{N90:.1f}",
                pw6=f"{PW6:.2f}", s12=f"{SIM[12]:.2f}", w3=f"{NW[0.3]:.1f}", w5=f"{NW[0.5]:.1f}", w7=f"{NW[0.7]:.1f}",
                nh=f"{NHALF:.0f}", pw6w=f"{PW6W[0.5]:.2f}")


def notes2_en(r):
    v = _n2()
    c = {k_: (x.replace("\u2212", "-") if isinstance(x, str) else x) for k_, x in r.items()}
    return (
        "WHAT: every claim with the label it earned in the checks (claims table: teaching package docs/06), and the "
        "sample size a follow-up cohort needs.\n"
        f"SUPPORTED: the port reproduces the paper's ranking on its 40 mice; day 14: Active licking unchanged "
        f"({c['d14a']}, p = {c['p14a']}) while Passive rose ({c['d14p']}); the first 7.5 min of day 14 did not fall "
        f"({c['h1']}) although {c['f75']}% of that day's drops had been earned (PR front-loads drops); Active fell on "
        f"days 15-16 and rose on day 17; at group level two rates beat one rate and the instant switch, and two rates "
        "are never picked for data made with one rate or the instant switch.\n"
        f"NOT SUPPORTED: the asymmetry as a direct test ({c['a14']} vs {c['a17']}, paired t p = {c['pa']}); the later "
        f"group-curve midpoint after drug loss ({c['mid_d']} vs {c['mid_e']} trials) is not significant per mouse "
        f"(p = {c['pm']}); the instant switch is not rejected at the mouse level (p = {c['pi']}, Holm {c['ph']}; "
        f"noise-matched simulations {c['nlo']}-{c['nhi']}); RPE 2alpha vs meta-RPE is a tie on held-out mice; "
        f"'Passive rose more than every model' depends on the unusually low day-13 reference (placebo day 13 - 12 "
        f"fails, p = {c['pp']}).\n"
        f"UNTESTABLE: error averaging: data made by meta-RPE are recovered in {c['r1']}% (per mouse) / {c['r2']}% "
        "(group) of simulations, because each mouse has one switch each way; a delayed oral-drug effect, the "
        "Withdrawal state and taste cannot be separated from learning with behaviour alone.\n"
        f"POWER (statsmodels; inputs from our data only): between-mouse SD of the change in p(lick) = {v['sd']} (mean "
        f"of 4 SDs recovered from the 95% CIs). Effect = the day-14 lag in Active mice: day-14 change {v['d14']} vs "
        f"mean of days 15-16 {v['late']} -> d = {v['d']}. Two-sample t, two-sided alpha 0.05: {v['n80']} per group for "
        f"80% ({v['n90']} for 90%); power at n = 6 = {v['pw6']}; the Welch simulation (20,000 draws) agrees: n = 12 "
        f"gives {v['s12']}. Within mouse (each mouse tested in both conditions; paired t, dz = d / sqrt(2(1 - r))): "
        f"80% at {v['w3']} / {v['w5']} / {v['w7']} mice for r = 0.3 / 0.5 / 0.7 (power at n = 6 with r = 0.5: "
        f"{v['pw6w']}). Half the effect, or a group x condition interaction (about 4x the n): {v['nh']} per group or "
        "cell.")


def notes2_kr(r):
    v = _n2()
    c = {k_: (x.replace("\u2212", "-") if isinstance(x, str) else x) for k_, x in r.items()}
    return (
        "무엇: 모든 주장을 검정 결과에 따라 붙인 라벨과 함께 정리했고 (주장 표: teaching package docs/06), 다음 코호트에 "
        "필요한 마우스 수를 계산했습니다.\n"
        f"지지됨: 우리 포트가 논문의 40마리에서 모델 순위를 재현합니다. Day 14에 Active licking은 변하지 않았고 "
        f"({c['d14a']}, p = {c['p14a']}) Passive는 늘었습니다 ({c['d14p']}). Day 14 처음 7.5분 동안 그날 drop의 "
        f"{c['f75']}%를 이미 얻었는데도 licking이 줄지 않았습니다 ({c['h1']}; PR은 drop이 세션 앞쪽에 몰림). Active는 "
        "day 15-16에 줄고 day 17에 늘었습니다. 그룹 단위에서 학습률 두 개가 하나짜리와 즉시 전환보다 낫고, 학습률 "
        "하나나 즉시 전환으로 만든 데이터에서는 학습률 두 개 모델이 한 번도 선택되지 않았습니다.\n"
        f"지지되지 않음: 비대칭을 직접 검정하면 유의하지 않습니다 ({c['a14']} vs {c['a17']}, 대응 t p = {c['pa']}). "
        f"약물 제거 후 그룹 곡선 중간점이 더 늦지만 ({c['mid_d']} vs {c['mid_e']} trial) 마우스별로는 유의하지 않습니다 "
        f"(p = {c['pm']}). 즉시 전환은 마우스 단위에서 기각되지 않습니다 (p = {c['pi']}, Holm {c['ph']}; 잡음 맞춘 "
        f"시뮬레이션 {c['nlo']}-{c['nhi']}). 남겨 둔 마우스에서 RPE 2alpha와 meta-RPE는 비깁니다. 'Passive가 모든 "
        f"모델보다 더 늘었다'는 day 13 기준일이 유난히 낮은 탓이 큽니다 (placebo day 13 - 12에서 실패, p = {c['pp']}).\n"
        f"검정 불가: 오류 평균은 meta-RPE로 만든 데이터조차 {c['r1']}% (마우스별) / {c['r2']}% (그룹)만 회복됩니다. "
        "마우스마다 방향별 전환이 한 번뿐이기 때문입니다. 늦게 작용하는 경구 약물 효과, Withdrawal 상태, 맛은 행동만으로 "
        "학습과 구분할 수 없습니다.\n"
        f"검정력 (statsmodels; 입력은 우리 데이터만): p(lick) 변화의 마우스 간 SD = {v['sd']} (95% CI에서 되돌린 SD "
        f"4개의 평균). 효과 = Active의 day 14 지연: day 14 변화 {v['d14']} vs day 15-16 평균 {v['late']} -> d = "
        f"{v['d']}. 두 표본 t, 양측 alpha 0.05: 80%에 그룹당 {v['n80']}마리 (90%는 {v['n90']}); n = 6의 검정력 = "
        f"{v['pw6']}; Welch 시뮬레이션 (20,000번)도 일치 (n = 12에서 {v['s12']}). 같은 마우스에 두 조건을 모두 주는 "
        f"설계 (대응 t, dz = d / sqrt(2(1 - r))): r = 0.3 / 0.5 / 0.7이면 80%에 {v['w3']} / {v['w5']} / {v['w7']}마리 "
        f"(r = 0.5, n = 6의 검정력 {v['pw6w']}). 효과가 절반이거나 그룹 x 조건 상호작용을 보려면 (n 약 4배) 그룹 또는 "
        f"셀당 {v['nh']}마리.")


# ============================================================================================ main
def main():
    F1, F2 = fig_s1(), fig_s2()
    outs = {}
    for lang in ("EN", "KR"):
        b.LANG = lang
        prs = Presentation(str(b.TEMPLATE))
        slide1(prs, F1)
        slide2(prs, F2)
        check(f"[{lang}] deck has 2 slides", len(prs.slides) == 2)
        outs[lang] = OUT / f"mrpe_summary_2slides_{lang}_20261126.pptx"
        prs.save(outs[lang])
    for k_, (txt, val) in b.PRINTED.items():
        if not k_.startswith("z"):
            continue
        t_ = txt.replace("\u2212", "-").replace("+", "").strip()
        if t_.startswith("<"):
            ok = val < float(t_[1:])
        else:
            nd = len(t_.split(".")[1]) if "." in t_ else 0
            ok = abs(float(t_) - val) <= 0.5 * 10 ** -nd + 1e-12
        check(f"printed {k_} = {txt} matches {val:.6g}", ok)
    json.dump(dict(checks=b.CHECKS, fit=b.FIT), open(OUT / "fit_report.json", "w", encoding="utf-8"), indent=1,
              ensure_ascii=False)
    bad = [c_ for c_ in b.CHECKS if not c_[1]]
    print(f"{len(b.CHECKS)} checks, {len(bad)} failed")
    for c_ in bad:
        print("  FAILED", c_[0])
    for lang in ("EN", "KR"):
        print(lang, [p_.name for p_ in b.render(outs[lang], OUT / f"preview_{lang}")])


if __name__ == "__main__":
    main()
