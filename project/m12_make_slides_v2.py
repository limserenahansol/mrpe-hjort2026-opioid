"""Step 12: slides v2 (EN + KR) after the independent review of v1 (v1 decks are kept unchanged).

What changed against v1 (m10_make_slides.py, whose helpers are reused here):
- slide 1: KEY says meta-RPE delays the value drop by only ~3 trials (trials 3-5: both models 0.72, licking 0.93);
  legend boxed; footer names the test difference with the paper (paired t + Holm vs RM ANOVA + Holm-Sidak); notes fixed.
- slide 2: Did line says the readout was replaced; readout row lists every term; n and "beta > 0 in all Active,
  < 0 in all Passive" on the slide; Korean KEY reworded; two arrows (both switches).
- slide 3: whole-session tests (no post-hoc 9-min window), Active vs Passive Welch test, CI for the Active null,
  days 15-16, and the non-significant asymmetry test printed; unrounded values (m11_modelfree.json).
- slide 4: title states a model fit, not a mechanism; leave-one-mouse-out CV added; caveats on the slide;
  Passive panel shows the best model on top; Did counts 6 learning rules + 2 controls.
- slide 5: group-level recovery panel (incl. RPE 2alpha); KEY gives both groups; percentages rounded half up.
- slide 6: prediction framed as conditional on meta-RPE (which lost), with its caveats; "or" for the gain M;
  ORBm/BMAp not studied in the paper; larger labels.
- KR: fitting = 피팅, mouse = 마우스, footers no longer drop content, notes match EN.
Output (slides_20261117\\): mrpe_hjort2026_opioid_6slides_EN_v2_20261117.pptx, ..._KR_v2_20261117.pptx,
        preview_EN_v2\\, preview_KR_v2\\, single_slides_v2\\, figs_v2\\, fit_report_v2.json
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m10_make_slides as b  # noqa: E402  (helpers, data loaders, LANG / SLIDE / checks live there)

plt, Rectangle, Line2D, Presentation = b.plt, b.Rectangle, b.Line2D, b.Presentation
P, Rn, tbox, measure, arrow, at, place, key_box, common = (b.P, b.Rn, b.tbox, b.measure, b.arrow, b.at, b.place,
                                                         b.key_box, b.common)
frac, sf, fp, peq, reg, check, T = b.frac, b.sf, b.fp, b.peq, b.reg, b.check, b.T
hm, m4, m5 = b.hm, b.m4, b.m5
D, OUT = b.D, b.OUT
FIG2 = OUT / "figs_v2"
for p in (FIG2, OUT / "single_slides_v2"):
    p.mkdir(parents=True, exist_ok=True)
RED, BLUE, GREY, GREEN, WBG, CMOD = b.RED, b.BLUE, b.GREY, b.GREEN, b.WBG, b.CMOD
S3, N7, G5, REC, BS, X6 = b.S3, b.N7, b.G5, b.REC, b.BS, b.X6
M11 = json.load(open(D / "m11_modelfree.json"))
GR = json.load(open(D / "m11_group_recovery.json"))
CV = json.load(open(D / "m11_lomo_cv.json"))
F4 = pd.read_csv(D / "m04_fits.csv")
F4["params"] = F4.params.apply(json.loads)
REPO = b.REPO


def save2(fig, name):
    p = FIG2 / name
    fig.savefig(p, facecolor="white")
    plt.close(fig)
    return p


def pct(v):
    """fraction -> whole percent, rounded half up (0.425 -> 43)"""
    return str(int(np.floor(v * 100 + 0.5)))


def beta_signs():
    """sign of the readout weight beta in the per-mouse meta-RPE fits, compared mice only (6099 white: no switch)."""
    out = {}
    for g in ("Active", "Passive"):
        bb = [F4[(F4.mouse_key == m) & (F4.model == "mRPE")].params.iloc[0][4] for m in m5.MICE[g]]
        out[g] = dict(n=len(bb), n_pos=int(sum(x > 0 for x in bb)), n_neg=int(sum(x < 0 for x in bb)))
    return out


BETA = beta_signs()
check("beta > 0 in all Active, < 0 in all Passive (compared mice)",
      BETA["Active"]["n_pos"] == BETA["Active"]["n"] == 6 and BETA["Passive"]["n_neg"] == BETA["Passive"]["n"] == 7)


# ============================================================================================ figures v2
def fig_paper():
    L, QR, QM, TM = b.paper_curves()
    Ln, Rn_, Mn = b.nz(L), b.nz(QR), b.nz(QM)
    x = np.arange(160) - 60
    win = (x >= -20) & (x <= 60)
    fig = plt.figure(figsize=(12.7, 4.85))
    a = fig.add_axes([0.062, 0.145, 0.375, 0.735])
    a.axvline(0, color="k", lw=1, ls="--")
    for v, col, lw, lab in ((Ln, "k", 2.2, "licking (running mean, 5 trials)"), (Rn_, GREY, 2.6, "RPE model value"),
                            (Mn, GREEN, 3.0, "meta-RPE model value")):
        m_ = v.mean(0)
        se = v.std(0, ddof=1) / np.sqrt(len(v))
        a.fill_between(x[win], (m_ - se)[win], (m_ + se)[win], color=col, alpha=0.15, lw=0)
        a.plot(x[win], m_[win], color=col, lw=lw, label=lab)
    a.axhline(1, color="#BBBBBB", lw=0.8, ls=":")
    a.axhline(0, color="#BBBBBB", lw=0.8, ls=":")
    a.set_xlim(-22, 62)
    a.set_ylim(-0.25, 1.32)
    a.set_yticks([0, 0.5, 1.0])
    a.set_xlabel("Trials from reversal (cue: 85% → 15% reward)")
    a.set_ylabel("Normalized (pre = 1, post = 0)")
    a.set_title("a  Licking and model value, 40 mice", loc="left", fontweight="bold")
    leg = a.legend(loc="lower left", frameon=True, framealpha=1.0, facecolor="white", edgecolor="none")
    leg.set_zorder(6)
    b1 = fig.add_axes([0.515, 0.66, 0.235, 0.22])
    d_mu = TM[:, :, 1].mean(0)
    b1.axvline(0, color="k", lw=1, ls="--")
    b1.axhline(0, color="k", lw=0.6)
    b1.fill_between(x[win], 0, d_mu[win], color="#C8C8C8", lw=0)
    b1.plot(x[win], d_mu[win], color="#555555", lw=1.6)
    b1.set_xlim(-22, 62)
    b1.set_ylabel("RPE δ")
    b1.tick_params(labelbottom=False)
    b1.set_title("b  Inside the meta-RPE model, 40 mice", loc="left", fontweight="bold")
    b2 = fig.add_axes([0.515, 0.145, 0.235, 0.47])
    cd, mm = TM[:, :, 2].mean(0), TM[:, :, 4].mean(0)
    b2.axvline(0, color="k", lw=1, ls="--")
    b2.axhline(0, color="k", lw=0.6)
    b2.plot(x[win], cd[win], color=RED, lw=2.4, label="CD = running mean of δ")
    b2.plot(x[win], mm[win], color=GREEN, lw=2.6, ls="--", label="gain M")
    pk = int(x[win][np.argmax(mm[win])])
    check("v2: M peak trial = m07 value", pk == N7["gain_M_peak_trial_from_reversal"])
    b2.annotate(f"M peaks at trial {pk}", xy=(pk, mm[60 + pk]), xytext=(pk + 9, mm[60 + pk] + 0.018), fontsize=11,
                color=GREEN, arrowprops=dict(arrowstyle="-", color=GREEN, lw=1))
    b2.set_xlim(-22, 62)
    b2.set_ylim(-0.215, 0.155)
    b2.set_xlabel("Trials from reversal")
    b2.set_ylabel("CD, gain M")
    b2.legend(loc="lower right", frameon=False, fontsize=10.5)
    c = fig.add_axes([0.872, 0.145, 0.123, 0.735])
    names = hm.ORDER
    lab = {"mRPE": "meta-RPE", "RPE2a": "RPE 2α", "PH": "PH", "eMack": "eMack", "RPEt-1": "RPE t−1", "RPE": "RPE"}
    vals = [S3["BIC"]["mean_change_from_best"][n] for n in names]
    sems = [S3["BIC"]["sem"][n] for n in names]
    yy = np.arange(len(names))
    c.barh(yy, vals, xerr=sems, color=[GREEN] + ["#A8A8A8"] * 5, capsize=3, height=0.65)
    for i, n in enumerate(names[1:], start=1):
        c.text(vals[i] + sems[i] + 0.6, i, "*" if S3["BIC"]["vs_mRPE"][n]["p_t_holm"] < 0.05 else "ns", va="center",
               fontsize=13)
    c.set_yticks(yy, [lab[n] for n in names])
    c.invert_yaxis()
    c.set_xlim(0, 25)
    c.set_xlabel("ΔBIC from best\n(lower = better)")
    c.set_title("c  Fit, 40 mice", loc="left", fontweight="bold")
    tg = {"key_tl": frac(fig, a, 13, 1.30), "key_br": frac(fig, a, 61, 0.68),
          "tip": frac(fig, a, 3.4, float(Mn.mean(0)[63]) + 0.02)}
    return save2(fig, "s1_paper_fig1.png"), fig.get_size_inches(), tg


def fig_map():
    fig = plt.figure(figsize=(5.55, 4.85))
    t = fig.add_axes([0.15, 0.625, 0.82, 0.25])
    days = list(range(3, 19))
    t.axvspan(10.5, 18.5, color="#FFF4CC", zorder=0)
    for row, g in ((1, "Active"), (0, "Passive")):
        for d in days:
            morph = (6 <= d <= 13) or d >= 17
            replay = g == "Passive" and 6 <= d <= 10
            t.add_patch(Rectangle((d - 0.43, row - 0.34), 0.86, 0.68, facecolor="#F2C4D6" if morph else "#CFE2F5",
                                  edgecolor="white", hatch="////" if replay else None, lw=0.5))
    t.text(8, -0.02, "replay", ha="center", va="center", fontsize=10, color="#333333",
           bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9))
    for xs, nm in ((13.5, "switch 1"), (16.5, "switch 2")):
        t.plot([xs, xs], [-0.5, 1.5], color="k", lw=1.6, ls="--")
        t.text(xs, 1.62, nm, ha="center", va="bottom", fontsize=10.5, fontweight="bold")
    t.text(12.0, -0.72, "fitted", ha="center", va="top", fontsize=10, color="#7A5C00")
    t.set_xlim(2.4, 18.6)
    t.set_ylim(-1.0, 1.55)
    t.set_xticks(days, [str(d) for d in days], fontsize=9.5)
    t.set_yticks([1, 0], ["Active", "Passive"], fontsize=11)
    for lbl, col in zip(t.get_yticklabels(), (RED, BLUE)):
        lbl.set_color(col)
        lbl.set_fontweight("bold")
    t.spines["left"].set_visible(False)
    t.tick_params(axis="y", length=0)
    t.set_xlabel("Day", fontsize=10.5, labelpad=1)
    fig.text(0.15, 0.965, "Outcome of every drop (pink = morphine, blue = water)", fontsize=11.5, fontweight="bold",
             va="top")
    s = fig.add_axes([0.15, 0.095, 0.82, 0.19])
    A14 = m5.EV["0910_black"]
    A14 = A14[A14.day == 14]
    Lk = pd.read_csv(D / "licks_all.csv")
    L14 = Lk[(Lk.mouse_key == "0910_black") & (Lk.day == 14)]
    bins = A14[A14.typ == 1]
    for _, rr in bins.iterrows():
        if rr.c == 1:
            s.add_patch(Rectangle((rr.t / 60, 0), 0.5, 1, facecolor="#D0D0D0", edgecolor="white", lw=0.6))
    s.vlines(L14[L14.kind == "lick"].t.values / 60, 0.50, 0.95, color="k", lw=0.5)
    dt = L14[L14.kind == "drop"].t.values / 60
    s.scatter(dt, np.full(len(dt), 0.25), marker="v", s=16, color="#1565C0", zorder=3)
    s.set_xlim(0, 15)
    s.set_ylim(0, 1)
    s.set_yticks([])
    s.spines["left"].set_visible(False)
    s.set_xticks([0, 3, 6, 9, 12, 15])
    s.set_xlabel("Minutes into the session", fontsize=10.5, labelpad=1)
    s.set_title("One session (Active 0910 black, day 14): grey = 30-s bin with a lick,\n| = lick, ▼ = water drop",
                fontsize=10, loc="left")
    tg = {"key_tl": (0.15, 1 - 0.535), "key_br": (0.97, 1 - 0.368), "sw1": frac(fig, t, 13.5, -0.5),
          "sw2": frac(fig, t, 16.5, -0.5)}
    return save2(fig, "s2_mapping.png"), fig.get_size_inches(), tg


def fig_behaviour():
    DAYS = list(range(11, 19))
    Y = X6[X6.day.isin(DAYS) & (X6.mouse != "6099_white")]
    fig = plt.figure(figsize=(12.7, 4.85))
    a = fig.add_axes([0.058, 0.145, 0.50, 0.735])
    for i, d in enumerate(DAYS):
        if 14 <= d <= 16:
            a.axvspan(i * 5 - 0.5, i * 5 + 4.5, color=WBG, zorder=0, lw=0)
        a.axvline(i * 5 - 0.5, color="#DDDDDD", lw=0.7, zorder=0)
        a.text(i * 5 + 2, 1.585, "water" if 14 <= d <= 16 else "morphine", ha="center", va="top", fontsize=10.5,
               color="#1565C0" if 14 <= d <= 16 else "#7B1FA2")
    blk = {}
    for g, col in (("Active", RED), ("Passive", BLUE)):
        W = Y[Y.group == g].pivot_table(index="mouse", columns=["day", "block"], values="p")
        m_, se = W.mean(), W.std(ddof=1) / np.sqrt(W.notna().sum())
        blk[g] = m_
        xx = np.arange(len(m_))
        for i in range(len(DAYS)):
            sl = slice(i * 5, i * 5 + 5)
            a.fill_between(xx[sl], (m_ - se).values[sl], (m_ + se).values[sl], color=col, alpha=0.18, lw=0)
            a.plot(xx[sl], m_.values[sl], color=col, lw=2.2, marker="o", ms=3.5,
                   label=f"{g} (n = {W.shape[0]})" if i == 0 else None)
    a.set_xticks([i * 5 + 2 for i in range(len(DAYS))], [f"day {d}" for d in DAYS])
    a.set_xlim(-0.7, len(DAYS) * 5 - 0.3)
    a.set_ylim(0, 1.62)
    a.set_yticks([0, 0.5, 1.0])
    a.set_ylabel("p(lick in a 30-s bin)")
    a.set_title("a  Licking in 3-min blocks of each 15-min session (mean ± SEM)", loc="left", fontweight="bold")
    a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.905), frameon=False)
    y14 = float(blk["Active"][(14, 2)])
    for k, (key, ref, d1, ttl, old) in enumerate((("day14_vs_13", 13, 14, "b  Day 14 vs 13\n(morphine → water)",
                                                   "day14_morphine_to_water"),
                                                  ("day17_vs_16", 16, 17, "c  Day 17 vs 16\n(water → morphine)",
                                                   "day17_water_to_morphine"))):
        bx = fig.add_axes([0.645 + k * 0.2, 0.145, 0.15, 0.735])
        bx.axhline(0, color="k", lw=0.8)
        for g, col, off in (("Active", RED, -0.08), ("Passive", BLUE, 0.08)):
            r_ = b.MF[old][g]
            bx.errorbar(np.arange(5) + off, r_["mean_change_by_block"], yerr=r_["sem_by_block"], color=col, lw=2,
                        marker="o", ms=5, capsize=3, label=f"{g} (n = {r_['n']})")
        wa, wp = M11[key]["Active"]["whole"], M11[key]["Passive"]["whole"]
        wel = M11[key]["Active_vs_Passive_whole"]["welch_p"]
        txt = (f"whole session vs 0:\nActive {sf(wa['mean'])}, {peq(wa['p'])}\nPassive {sf(wp['mean'])}, "
               f"{peq(wp['p'])}\nActive vs Passive: {peq(wel)}")
        bx.text(0.03, 0.03 if k == 0 else 0.97, txt, transform=bx.transAxes, fontsize=9.8,
                va="bottom" if k == 0 else "top")
        bx.set_xticks(range(5), ["3", "6", "9", "12", "15"])
        bx.set_xlim(-0.5, 4.5)
        bx.set_ylim(-0.82, 1.08)
        bx.set_xlabel("Minute (block end)")
        bx.set_ylabel(f"Change in p(lick) vs day {ref}")
        bx.set_title(ttl, loc="left", fontweight="bold", fontsize=12.5)
        if k == 0:
            bx.legend(loc="upper right", frameon=False, fontsize=10)
    tg = {"key_tl": frac(fig, a, 9.6, 1.47), "key_br": frac(fig, a, 34.4, 1.03), "tip": frac(fig, a, 17, y14 + 0.07)}
    return save2(fig, "s3_behaviour_switches.png"), fig.get_size_inches(), tg


def fig_models(PRED):
    DD = list(range(13, 19))
    fig = plt.figure(figsize=(12.7, 4.85))
    tg = {}
    for k, (g, col) in enumerate((("Active", RED), ("Passive", BLUE))):
        a = fig.add_axes([0.055 + k * 0.345, 0.20, 0.30, 0.68])
        Yg = PRED[(PRED.group == g) & PRED.day.isin(DD)]
        for i, d in enumerate(DD):
            if 14 <= d <= 16:
                a.axvspan(i * 5 - 0.5, i * 5 + 4.5, color=WBG, zorder=0, lw=0)
            a.axvline(i * 5 - 0.5, color="#DDDDDD", lw=0.7, zorder=0)
        obs = Yg[Yg.model == "RPE"].pivot_table(index="mouse", columns=["day", "block"], values="obs")
        m_, se = obs.mean(), obs.std(ddof=1) / np.sqrt(len(obs))
        xx = np.arange(len(m_))
        curves = {mm: Yg[Yg.model == mm].pivot_table(index="mouse", columns=["day", "block"], values="pred").mean()
                  for mm in ("RPE", "mRPE", "RPE2a", "today's outcome")}
        # the group's best model is drawn last (on top); in Passive it is also thicker, since all lines overlap
        sty = {"RPE": ("RPE (one learning rate)", 1.9, "-"), "mRPE": ("meta-RPE (paper)", 2.6, "-"),
               "RPE2a": ("RPE 2α (rates for − / + errors)", 2.6, "-"),
               "today's outcome": ("today's outcome (instant switch)", 2.0, (0, (4, 2)))}
        best = "RPE2a" if g == "Active" else "today's outcome"
        if g == "Passive":
            sty["today's outcome"] = ("today's outcome (instant switch)", 3.0, (0, (3, 1.6)))
        order = tuple((mm,) + sty[mm] for mm in ("RPE", "mRPE", "RPE2a", "today's outcome") if mm != best) + \
            ((best,) + sty[best],)
        for i in range(len(DD)):
            sl = slice(i * 5, i * 5 + 5)
            a.errorbar(xx[sl], m_.values[sl], yerr=se.values[sl], color="k", lw=0, elinewidth=1.2, marker="o", ms=4.5,
                       capsize=2, zorder=6, label="observed (mean ± SEM)" if i == 0 else None)
            for z, (model, lab, lw, ls) in enumerate(order):
                a.plot(xx[sl], curves[model].values[sl], color=CMOD[model], lw=lw, ls=ls, zorder=3 + z * 0.1,
                       label=lab if i == 0 else None)
        a.set_xticks([i * 5 + 2 for i in range(len(DD))],
                     [f"{d}\n{'water' if 14 <= d <= 16 else 'morphine'}" for d in DD], fontsize=10.5)
        a.set_xlim(-0.7, len(DD) * 5 - 0.3)
        a.set_ylim(0, 1.62)
        a.set_yticks([0, 0.5, 1.0])
        a.set_ylabel("p(lick in a 30-s bin)" if k == 0 else "")
        a.set_xlabel("Day (3-min blocks)", fontsize=12)
        a.set_title(f"{'ab'[k]}  {g} (n = {len(obs)})", loc="left", fontweight="bold")
        if k == 1:
            h, l_ = a.get_legend_handles_labels()
            keep = ["observed (mean ± SEM)", "today's outcome (instant switch)", "RPE (one learning rate)",
                    "RPE 2α (rates for − / + errors)", "meta-RPE (paper)"]
            hh = [h[l_.index(x)] for x in keep]
            a.legend(hh, keep, loc="upper center", bbox_to_anchor=(0.5, 1.005), ncol=2, frameon=False, fontsize=10,
                     columnspacing=1.0, handlelength=2.2)
            a.text(14.5, 1.10, "all 4 model lines overlap; best = instant switch (purple, on top)", ha="center",
                   va="center", fontsize=9.5, color="#555555")
        else:
            y = float(curves["RPE2a"][(14, 2)])
            tg = {"key_tl": frac(fig, a, -0.4, 1.60), "key_br": frac(fig, a, 29.4, 1.05), "tip": frac(fig, a, 7, y + 0.06)}
    c = fig.add_axes([0.862, 0.20, 0.133, 0.68])
    models = ["RPE2a", "mRPE", "eMack", "RPE", "RPEt-1", "today's outcome"]
    nm = {"RPE2a": "RPE 2α", "mRPE": "meta-RPE", "eMack": "eMack", "RPE": "RPE", "RPEt-1": "RPE t−1",
          "today's outcome": "today's outcome"}
    yy = np.arange(len(models))
    for g, col, off in (("Active", RED, -0.2), ("Passive", BLUE, 0.2)):
        vals = [G5[g][mm]["dBIC_from_best"] for mm in models]
        c.barh(yy + off, [min(v, 33) for v in vals], height=0.38, color=col, label=g)
        for yv, v in zip(yy + off, vals):
            c.text(min(v, 33) + 0.8, yv, f"{v:.0f}", va="center", fontsize=9.5, color=col)
    c.axvline(10, color="k", lw=0.9, ls=":")
    c.text(10.6, len(models) - 0.45, "10", fontsize=9.5, color="#333333", va="center")
    c.set_yticks(yy, [nm[mm] for mm in models], fontsize=11)
    c.invert_yaxis()
    c.set_xlim(0, 38)
    c.set_xticks([0, 10, 20, 30])
    c.set_xlabel("ΔBIC from best\n(lower = better)", fontsize=12)
    c.set_title("c  Group fits", loc="left", fontweight="bold")
    c.legend(loc="upper right", frameon=False, fontsize=10)
    return save2(fig, "s4_group_model_fits.png"), fig.get_size_inches(), tg


def fig_recovery():
    fig = plt.figure(figsize=(12.7, 4.85))
    a = fig.add_axes([0.075, 0.20, 0.255, 0.68])
    conf = REC["confusion_BIC"]
    gens = ["mRPE", "RPE", "today's outcome"]
    for j, (g, col) in enumerate((("Active", RED), ("Passive", BLUE))):
        vals = [conf[g][gen].get(gen, 0.0) for gen in gens]
        xs = np.arange(3) + (j - 0.5) * 0.38
        a.bar(xs, [v * 100 for v in vals], width=0.38, color=col, label=g)
        for xv, v in zip(xs, vals):
            a.text(xv, v * 100 + 1.5, f"{pct(v)}%", ha="center", va="bottom", fontsize=10.5)
    a.set_xticks(range(3), ["meta-\nRPE", "RPE", "instant\nswitch"])
    a.set_ylim(0, 132)
    a.set_yticks([0, 25, 50, 75, 100])
    a.set_ylabel("% of simulated sets where\nBIC picks the true model")
    a.set_xlabel("Model that made the data")
    a.set_title("a  Per mouse (20 sets each)", loc="left", fontweight="bold")
    a.legend(loc="upper right", ncol=1, frameon=False, fontsize=10.5)
    # group-level confusion, Active
    h = fig.add_axes([0.445, 0.20, 0.185, 0.68])
    rows = ["RPE2a", "mRPE", "RPE", "today's outcome"]
    lab = {"RPE2a": "RPE 2α", "mRPE": "meta-RPE", "RPE": "RPE", "today's outcome": "instant"}
    Cm = np.array([[GR["confusion"]["Active"][gen].get(pk, 0.0) for pk in rows] for gen in rows])
    h.imshow(Cm, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(4):
        for j in range(4):
            v = Cm[i, j]
            h.text(j, i, f"{pct(v)}%", ha="center", va="center", fontsize=10.5, color="white" if v > 0.55 else "#222222",
                   fontweight="bold" if i == j else "normal")
        h.add_patch(Rectangle((i - 0.5, i - 0.5), 1, 1, fill=False, edgecolor="k", lw=1.4))
    h.set_xticks(range(4), [lab[r] for r in rows], rotation=35, ha="right", fontsize=10.5)
    h.set_yticks(range(4), [lab[r] for r in rows], fontsize=10.5)
    h.set_ylabel("Model that made the data", fontsize=11.5)
    h.set_xlabel("Model BIC picked", fontsize=11.5)
    h.set_title("b  Group level, Active (50 sets)", loc="left", fontweight="bold", fontsize=12.5)
    for sp in h.spines.values():
        sp.set_visible(False)
    cx = fig.add_axes([0.80, 0.20, 0.19, 0.68])
    ylab = []
    for j, (g, col) in enumerate((("Active", RED), ("Passive", BLUE))):
        r_ = BS[g]
        for k2, (key, lab2, mk) in enumerate((("a_minus", "α− (drug → water)", "v"),
                                              ("a_plus", "α+ (water → drug)", "^"))):
            lo, med, hi = r_[f"{key}_CI"]
            full = r_["full_fit"][key]
            y = j * 2.7 + k2
            cx.errorbar(full, y, xerr=[[max(full - lo, 0)], [max(hi - full, 0)]], color=col, marker=mk, ms=8,
                        capsize=3, lw=2)
            cx.text(max(lo, 0.0095), y - 0.3, f"{full:.2f} [{lo:.2f}, {hi:.2f}]", ha="left", fontsize=9.5, color=col)
            ylab.append((y, f"{g} {lab2}"))
        cx.text(0.0092, j * 2.7 + 1.62, f"α− < α+: {r_['frac_a_minus_less_than_a_plus'] * 100:.1f}%",
                fontsize=10, color=col)
    cx.set_xscale("log")
    cx.set_xlim(0.008, 1.6)
    cx.set_ylim(4.55, -0.85)
    cx.set_yticks([y for y, _ in ylab], [t_ for _, t_ in ylab], fontsize=9.8)
    cx.set_xlabel("Learning rate (log scale)", fontsize=12)
    cx.set_title("c  RPE 2α group fit, 95% CI", loc="left", fontweight="bold", fontsize=12.5)
    tg = {"key_tl": frac(fig, a, -0.5, 128), "key_br": frac(fig, a, 1.5, 58),
          "tip": frac(fig, a, -0.19, conf["Active"]["mRPE"].get("mRPE", 0) * 100 + 7)}
    return save2(fig, "s5_recovery_and_asymmetry.png"), fig.get_size_inches(), tg


def fig_vars():
    fig = plt.figure(figsize=(9.25, 4.85))
    tg, info = {}, {}
    for k, (m, g) in enumerate((("0910_black", "Active"), ("0910_orange", "Passive"))):
        Tb = b.drop_variables(m, g)
        xx = np.arange(len(Tb))
        first14 = int(Tb.index[Tb.day == 14][0])
        last16 = int(Tb.index[Tb.day == 16][-1])
        imin = int(Tb.CD.values[first14:last16 + 1].argmin()) + first14
        info[g] = dict(mouse=m, n_water=last16 - first14 + 1, cd_min_after=imin - first14 + 1)
        for col_i, kind in enumerate(("val", "cd")):
            ax = fig.add_axes([0.085 + col_i * 0.50, 0.565 - k * 0.455, 0.395, 0.285])
            for d in sorted(Tb.day.unique()):
                idx = Tb.index[Tb.day == d]
                if 14 <= d <= 16:
                    ax.axvspan(idx[0] - 0.5, idx[-1] + 0.5, color=WBG, zorder=0, lw=0)
                ax.axvline(idx[0] - 0.5, color="#CCCCCC", lw=0.6, zorder=0)
                ax.text((idx[0] + idx[-1]) / 2, 1.03, str(d), transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=9.5, color="#444444")
            if kind == "val":
                ax.bar(xx, Tb.delta, color="#BDBDBD", width=1.0, label="RPE δ per drop")
                ax.plot(xx, Tb.Q, color="k", lw=2, label="value Q (expected drug)")
                ax.set_ylim(-1.1, 1.15)
                ax.set_ylabel("Q, δ", fontsize=11.5)
                ax.legend(loc="lower left", frameon=False, fontsize=9.8)
                ttl = f"{'ac'[k]}  {g} {m.replace('_', ' ')}: value and RPE (day above)"
            else:
                ax.axhline(0, color="k", lw=0.6)
                ax.plot(xx, Tb.CD, color=RED, lw=2.2, label="CD (slow mean of δ)")
                ax.plot(xx, Tb.CE, color=BLUE, lw=1.5, label="CE (fast mean of δ)")
                ax.plot(xx, Tb.M, color=GREEN, lw=2.2, ls="--", label="gain M")
                ax.scatter([imin], [Tb.CD[imin]], s=40, color=RED, zorder=5, edgecolors="k", linewidths=0.6)
                ax.set_ylabel("CD, CE, M", fontsize=11.5)
                ax.legend(loc="lower left" if g == "Active" else "lower right", frameon=False, fontsize=9.8)
                ttl = f"{'bd'[k]}  {g}: averaged errors and gain"
                lo_, hi_ = ax.get_ylim()
                ax.set_ylim(lo_ - 0.12 * (hi_ - lo_), hi_ + 0.10 * (hi_ - lo_))
                tg[g] = frac(fig, ax, imin, float(Tb.CD[imin]))
            ax.set_title(ttl, loc="left", fontsize=11.5, fontweight="bold", pad=21)
            ax.tick_params(labelsize=10.5)
            ax.set_xlabel("Drops in order, days 11–18" if k == 1 else "", fontsize=11)
    # the same drop counts in every mouse of a group (shared parameters)?
    rng = {}
    for g in ("Active", "Passive"):
        v = []
        for m in m5.MICE[g]:
            Tb = b.drop_variables(m, g)
            f14 = int(Tb.index[Tb.day == 14][0])
            l16 = int(Tb.index[Tb.day == 16][-1])
            v.append(int(Tb.CD.values[f14:l16 + 1].argmin()) + 1)
        rng[g] = (min(v), max(v))
    info["range"] = rng
    return save2(fig, "s6_model_variables.png"), fig.get_size_inches(), tg, info


# ============================================================================================ slides v2
def slide1(prs, F):
    b.SLIDE = 1
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    B = S3["BIC"]
    others = [B["mean_change_from_best"][n] for n in hm.ORDER[1:]]
    hp = max(B["vs_mRPE"][n]["p_t_holm"] for n in hm.ORDER[1:])
    mb = reg("v2s1_bic", f"{B['mean_change_from_best']['mRPE']:.1f}", B["mean_change_from_best"]["mRPE"])
    lo_ = reg("v2s1_lo", f"{min(others):.1f}", min(others))
    hi_ = reg("v2s1_hi", f"{max(others):.1f}", max(others))
    hps = reg("v2s1_hp", f"{hp:.3f}", hp)
    nb = S3["negLL"]["n_best"]["mRPE"]
    km, kr_, kl = "meta-RPE model value", "RPE model value", "licking (running mean, 5 trials)"
    v = {f"{w}_{t}": N7[w2][t] for w, w2 in (("m", km), ("r", kr_), ("l", kl)) for t in ("trials_0_2", "trials_3_5")}
    s_ = {k_: reg(f"v2s1_{k_}", f"{val:.2f}", val) for k_, val in v.items()}
    check("v2 slide 1: both models 0.72 in trials 3-5", s_["m_trials_3_5"] == s_["r_trials_3_5"] == "0.72")
    pk = N7["gain_M_peak_trial_from_reversal"]
    ks = 13.5 if b.LANG == "EN" else 12
    key_box(s, box, tg, [P(Rn(T("Meta-RPE delays the value drop by ~3 trials", "meta-RPE는 가치 하락을 ~3 trial 늦춘다"),
                              ks + 0.5, True, "1F4E79"), after=3),
                         P(Rn(T(f"Trials 0–2: meta-RPE {s_['m_trials_0_2']}, RPE {s_['r_trials_0_2']}, licking "
                                f"{s_['l_trials_0_2']}. Trials 3–5: both {s_['m_trials_3_5']}, licking {s_['l_trials_3_5']}.",
                                f"trial 0–2: meta-RPE {s_['m_trials_0_2']}, RPE {s_['r_trials_0_2']}, lick "
                                f"{s_['l_trials_0_2']}. trial 3–5: 둘 다 {s_['m_trials_3_5']}, lick {s_['l_trials_3_5']}."),
                              ks))], tip=tg["tip"], anchor="left")
    common(s, T("Hjort et al. 2026: behaviour changes only after errors pile up",
                "Hjort et al. 2026: 오류가 쌓인 뒤에야 행동이 바뀐다"),
           T("re-fitted the authors' 40 mice with their own code (ported to Python).",
             "저자 코드(Python으로 옮김)로 저자 데이터 40마리를 다시 피팅함."),
           T(f"the paper's ranking reproduces: meta-RPE fits best (ΔBIC {mb} vs {lo_}–{hi_}; better than each model, "
             f"largest Holm p = {hps}).",
             f"논문의 순위가 재현됨: meta-RPE가 가장 잘 맞음 (ΔBIC {mb} vs {lo_}–{hi_}; 모든 모델보다 좋음, 가장 큰 "
             f"Holm p = {hps})."),
           T(f"Hjort et al., Nature 655:174 (2026), Fig. 1; code + data: {REPO} (H→L cue, 160 trials, reversal at trial "
             f"60). a: per mouse, trials 20–59 = 1, 120–159 = 0; mean ± SEM. b: step = α·δ·M, "
             f"M = max(0, −CD) + max(0, CE). c: BIC = k ln T + 2 negLL; * paired t vs meta-RPE, Holm over 5 (paper: RM ANOVA "
             f"+ Holm–Šídák; same ranking, p values differ).",
             f"Hjort et al., Nature 655:174 (2026), Fig. 1; 코드·데이터: {REPO} (H→L 큐, 160 trial, trial 60에서 "
             f"reversal). a: 마우스마다 trial 20–59 = 1, 120–159 = 0; 평균 ± SEM. b: 학습 크기 = α·δ·M, "
             f"M = max(0, −CD) + max(0, CE). c: BIC = k ln T + 2 negLL; * meta-RPE와 paired t, Holm 보정 (논문은 RM ANOVA: 순위 같고 p값 다름)."),
           T(notes1_en(mb, lo_, hi_, hps, nb, s_, pk), notes1_kr(mb, lo_, hi_, hps, nb, s_, pk)))
    return s


def notes1_en(mb, lo_, hi_, hps, nb, s_, pk):
    return (
        "WHAT: Fig. 1 of Hjort et al. 2026 (Nature 655:174) redrawn from the authors' 40 sample mice with their own model code "
        "(BehaviorModels.ipynb, ported line by line to Python + numba; same equations, bounds, start grid and BIC).\n"
        "TASK: head-fixed mice, 4 odour cues (2 at 85% sucrose, 2 at 15%). After stable licking, one unsignalled reversal: one "
        "85% cue becomes 15% (H->L, contingency degradation) and one 15% cue becomes 85% (L->H, enhancement); the other two "
        "(H->H, L->L) are controls. Behaviour = anticipatory licks 0-2.5 s after cue onset. The model is fitted to the H->L cue: "
        "160 trials (the authors' notebook: 4 days of 40 trials, reversal at trial 60).\n"
        f"a: licking stays high for ~6 trials. In trials 0-2 the RPE value is already {s_['r_trials_0_2']} while the meta-RPE "
        f"value is {s_['m_trials_0_2']} (licking {s_['l_trials_0_2']}); by trials 3-5 both models are at {s_['m_trials_3_5']} "
        f"while licking is still {s_['l_trials_3_5']}. So meta-RPE delays the drop by ~3 trials, not as long as the mice.\n"
        f"b: the single-trial RPE is largest right after the reversal; the averaged error CD needs several trials to build, so "
        f"the gain M peaks at trial {pk}. Learning step = alpha x delta x M; M rises when CD < 0 or CE > 0.\n"
        f"c: our port gives the paper's ranking (Fig. 1k,l): BIC change from best {mb} (meta-RPE) vs {lo_}-{hi_}; paired t vs "
        f"meta-RPE with Holm correction, largest p = {hps} (the paper used RM ANOVA + Holm-Sidak, so its p values differ); "
        f"by negLL meta-RPE is best in {nb}/40 mice.\n"
        "Rest of the paper: 13% of mPFC neurons encode CD (GLM + circular-shift null); holographic activation of CD cells "
        "reduces licking; mPFC dopamine carries RPE/CE signals; 61% of mPFC->VTA cells encode CD; D1/D2 block in mPFC lowers CD "
        "encoding and delays reversal; mPFC->VTA stimulation speeds degradation. The authors speculate that CD-encoding cells "
        "might be impaired after long-term voluntary drug taking.")


def notes1_kr(mb, lo_, hi_, hps, nb, s_, pk):
    return (
        "[이 슬라이드] Hjort et al. 2026 (Nature 655:174) Figure 1의 핵심을 저자 샘플 데이터 40마리와 저자 코드로 다시 그린 것. "
        "저자 노트북(BehaviorModels.ipynb)을 Python으로 한 줄씩 옮겼고 식, 파라미터 범위, 시작값, BIC 계산은 모두 같다.\n"
        "[과제] 머리 고정 마우스에게 냄새 큐 4개: 2개는 85% 확률로 설탕물, 2개는 15%. lick이 안정되면 예고 없이 한 번 바꾼다: "
        "85% 큐 하나는 15%로(H→L, contingency degradation), 15% 큐 하나는 85%로(L→H, enhancement) 바뀌고 나머지 두 큐"
        "(H→H, L→L)는 대조로 그대로 둔다. 행동 지표는 큐 뒤 0–2.5초의 미리 핥기(anticipatory lick). 모델은 H→L 큐의 160 "
        "trial (저자 노트북: 하루 40 trial씩 4일, trial 60에서 reversal)에 피팅한다.\n"
        f"[a 읽는 법] 검은 선 = 실제 lick, 회색 = 일반 RPE 모델의 가치, 초록 = meta-RPE 모델의 가치 (마우스마다 reversal 전 = 1, "
        f"후 = 0으로 정규화). reversal 직후 trial 0–2: RPE 가치 {s_['r_trials_0_2']}, meta-RPE 가치 {s_['m_trials_0_2']}, lick "
        f"{s_['l_trials_0_2']}. trial 3–5에는 두 모델 모두 {s_['m_trials_3_5']}, lick은 아직 {s_['l_trials_3_5']}. 즉 meta-RPE는 "
        f"가치 하락을 약 3 trial 늦출 뿐이고, 실제 마우스(약 6 trial)만큼 오래 버티지는 않는다.\n"
        f"[b 읽는 법] 위: 한 trial의 오류 δ = 보상 − 가치. reversal 바로 다음 trial에서 가장 크다. 아래: CD = δ의 이동평균(빨강), "
        f"초록 점선 = gain M = max(0, −CD) + max(0, CE). 한 번에 배우는 양 = α × δ × M. CD < 0이거나 CE > 0일 때 M이 커지는데, "
        f"오류가 쌓여야 CD가 커지므로 M은 trial {pk}에서 최대가 된다.\n"
        f"[c 읽는 법] 6개 모델의 BIC 차이 (낮을수록 좋음). meta-RPE {mb} vs 나머지 {lo_}–{hi_}; 모든 모델보다 좋음 (paired t, "
        f"Holm 보정, 가장 큰 p = {hps}; 논문은 RM ANOVA + Holm–Šídák이라 p값은 다르지만 순위는 같다). negLL로는 40마리 중 "
        f"{nb}마리에서 meta-RPE가 최선. 우리 Python 코드가 논문 결과를 재현한다는 확인이다.\n"
        "[논문의 나머지] Fig. 2: mPFC 뉴런의 13%가 CD를 표현 (GLM + circular shift 검정). Fig. 3: CD 세포를 홀로그래픽 광유전학으로 "
        "켜면 lick 감소. Fig. 4: mPFC 도파민은 RPE/CE 신호를 담고, mPFC→VTA 뉴런의 61%가 CD를 표현; mPFC에 D1/D2 차단제를 넣으면 "
        "CD 표현이 줄고 reversal이 늦어짐. Fig. 5: mPFC→VTA 자극 → degradation이 빨라짐. 저자들은 약물을 오래 자발적으로 섭취하면 "
        "CD 세포의 기능이 떨어질 수 있다고 추측한다.")


ROWS = [
    (("Trial", "시행 (trial)"), ("cue every ~20 s", "~20초마다 큐"), ("30-s bin of the 15-min session", "15분 세션의 30초 bin")),
    (("Choice", "선택"), ("anticipatory lick: yes / no", "미리 핥기: 있음 / 없음"), ("any lick in the bin: yes / no",
                                                                            "bin 안의 lick: 있음 / 없음")),
    (("Outcome", "결과"), ("sucrose 1 / none 0 (85% ↔ 15%)", "설탕물 1 / 없음 0 (85% ↔ 15%)"),
     ("each drop: morphine 1 / water 0", "매 drop: 모르핀 1 / 물 0")),
    (("Value Q", "가치 Q"), ("expected reward of the cue", "큐의 보상 기대"),
     ("expected drug from the spout; updated at every drop, from 0 on day 3",
      "spout의 약물 기대; drop마다 업데이트, day 3에 0에서 시작")),
    (("Switch", "전환"), ("one reversal: H→L, L→H", "reversal 1번: H→L, L→H"),
     ("day 14 morphine → water; day 17 water → morphine", "day 14 모르핀 → 물; day 17 물 → 모르핀")),
    (("Fitted", "피팅 데이터"), ("160 trials (4 days) per mouse", "마우스당 160 trial (4일)"),
     ("days 11–18: 240 bins per mouse", "day 11–18: 마우스당 240 bin")),
    (("Readout", "행동 출력"), ("softmax(β·Q), β > 0", "softmax(β·Q), β > 0"),
     ("logistic(b0 + β·Q + time + first 3 min); β of either sign",
      "logistic(b0 + β·Q + 시간 + 처음 3분); β 부호 자유")),
]


def slide2(prs, F):
    b.SLIDE = 2
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 7.55, 0.98, 5.55)
    xs, ws = [0.35, 1.70, 4.20], [1.35, 2.50, 2.95]
    fs = 13 if b.LANG == "EN" else 12
    head = ("", "Hjort 2026 (Pavlovian)", T("Our task (PR licking)", "우리 과제 (PR licking)"))
    y = 1.05
    hh = max(measure([P(Rn(h_, fs, True))], w) for h_, w in zip(head, ws)) + 0.06
    for x, w, h_ in zip(xs, ws, head):
        tbox(s, x, y, w, hh, [P(Rn(h_, fs, True, "FFFFFF"))], f"table head {h_}", fill="1F4E79", anchor=b.MSO_ANCHOR.MIDDLE)
    y += hh + 0.04
    for i, row in enumerate(ROWS):
        cells = [T(*row[0]), T(*row[1]), T(*row[2])]
        rh = max(measure([P(Rn(cl, fs, j == 0))], w) for j, (cl, w) in enumerate(zip(cells, ws))) + 0.06
        for j, (x, w, cl) in enumerate(zip(xs, ws, cells)):
            tbox(s, x, y, w, rh, [P(Rn(cl, fs, j == 0, "1F4E79" if j == 0 else "222222"))], f"table r{i}c{j}",
                 fill="EEF3F9" if i % 2 == 0 else "FFFFFF", anchor=b.MSO_ANCHOR.MIDDLE)
        y += rh + 0.04
    check(f"[{b.LANG}] v2 mapping table ends above the bottom line (ends {y:.2f} in)", y <= 5.88)
    ks = 12 if b.LANG == "EN" else 11.5
    x0, y0 = at(box, tg["key_tl"])
    x1, y1 = at(box, tg["key_br"])
    paras = [P(Rn(T("Days 11–18: same drops for both groups on the same days, so both are compared at the same "
                    "two switches.",
                    "day 11–18: 두 그룹 모두 같은 날 같은 종류의 drop(모르핀/물)을 받으므로 같은 두 전환에서 비교한다."),
                  ks, True, "1F4E79"))]
    tbox(s, x0, y0, x1 - x0, None, paras, "KEY", fill="F2F6FB", line="1F4E79", region_bottom=y1)
    for k_ in ("sw1", "sw2"):
        tx, ty = at(box, tg[k_])
        arrow(s, tx, y0, tx, ty + 0.05)
    nb = int(G5["Active"]["RPE"]["n"] / len(m5.MICE["Active"]))
    reg("v2s2_bins", str(nb), nb)
    na, npv = BETA["Active"]["n"], BETA["Passive"]["n"]
    common(s, T("Same equations, different task: the model on our PR licking", "식은 같고 과제는 다르다: 우리 PR licking에 맞춘 방법"),
           T("kept the paper's learning-rule equations and bounds; replaced the readout (logistic with intercept, time terms, "
             "β of either sign); redefined trial, choice and outcome.",
             "논문의 학습 규칙 식과 범위는 유지; 행동 출력을 logistic(절편, 시간 항, β 부호 자유)으로 바꾸고 trial, 선택, "
             "결과를 다시 정의함."),
           T(f"each mouse has one morphine → water and one water → morphine switch ({nb} bins); β > 0 in all {na} "
             f"Active and β < 0 in all {npv} Passive mice.",
             f"마우스마다 모르핀 → 물, 물 → 모르핀 전환이 한 번씩 ({nb} bin); β는 Active {na}마리 모두 > 0, "
             f"Passive {npv}마리 모두 < 0."),
           T("n = 6 Active, 8 Passive (7 after day 13: 6099 white died, left out of model comparisons). PR = progressive ratio "
             "(+1 lick per drop), ~15-min sessions. Passive days 6–10 = replay of the partner's drops (morphine drunk, licks "
             "earn nothing): drops update value, licking is not fitted. β < 0 = licks less when expecting morphine, so for "
             "Passive mice morphine → water is not a reward omission. Raw controller logs (K:, read only).",
             "n = Active 6, Passive 8 (day 13 이후 7: 6099 white 사망, 모델 비교에서 제외). PR: drop마다 lick +1, 세션 "
             "~15분. Passive day 6–10 = 짝의 drop 재생(replay): 모르핀은 마셨지만 lick은 보상과 무관 → "
             "drop은 가치 업데이트에 넣고 lick은 피팅하지 않음. β < 0 = 모르핀을 기대할수록 덜 핥음: Passive에게 모르핀 "
             "→ 물은 보상 상실이 아님. 원자료: K: 컨트롤러 로그 (읽기 전용)."),
           T(notes2_en(nb, na, npv), notes2_kr(nb, na, npv)))
    return s


def notes2_en(nb, na, npv):
    return (
        "WHAT: what had to change to use the paper's model on our data.\n"
        "SAME: the learning-rule equations (RPE, meta-RPE, RPE 2alpha, RPE t-1, Pearce-Hall, eMack) and their parameter "
        "bounds, maximum-likelihood fitting and BIC.\n"
        "CHANGED: (1) the task: the paper changes the reward PROBABILITY of a cue (Pavlovian); in our free-operant PR task every "
        "completed ratio gives a drop and its IDENTITY changes (morphine <-> water). So trial = 30-s bin, choice = any lick in "
        "the bin, outcome = drop identity (morphine 1, water 0). (2) the readout: the paper uses softmax(beta*Q) with beta > 0, "
        "start values 0.15 / 0.85 and updates only the chosen option; we use logistic(b0 + beta*Q + time in session + first "
        "3 min) with beta of either sign, value 0 on day 3, updated at every drop (in an operant task a drop only follows "
        "licking).\n"
        f"VALUE Q = expected drug from the spout, carried over between days; Passive replay drops (days 6-10) update it too, "
        f"but their replay-day licking is not fitted. Fitted days 11-18 ({nb} bins per mouse), where both groups run the same "
        f"PR schedule with the same outcome sequence.\n"
        f"SIGN OF BETA: positive in all {na} Active mice, negative in all {npv} Passive mice (per-mouse fits): Passive mice lick "
        "less when they expect morphine, so for them day 14 (morphine -> water) is relief rather than a reward omission, and "
        "the H->L analogy holds only for Active mice.\n"
        "Right, bottom: one real session (0910 black, day 14): grey bins had at least one lick; triangles = water drops.")


def notes2_kr(nb, na, npv):
    return (
        "[이 슬라이드] 논문 모델을 우리 데이터에 쓰려면 무엇을 바꿔야 했는지.\n"
        "[같은 것] 학습 규칙 식 6개(RPE, meta-RPE, RPE 2α, RPE t−1, Pearce–Hall, eMack)와 파라미터 범위, 최대우도 피팅, BIC 비교.\n"
        "[바꾼 것] (1) 과제: 논문은 큐의 '보상 확률'이 바뀌는 Pavlovian 과제. 우리는 free-operant PR 과제라서 정해진 수만큼 "
        "핥으면 drop이 항상 나오고, 그 drop의 '종류'(모르핀 ↔ 물)가 바뀐다. 그래서 trial = 30초 bin, 선택 = bin 안에 lick이 "
        "있었는지, 결과 = drop의 종류(모르핀 1, 물 0). (2) 행동 출력: 논문은 softmax(β·Q), β > 0, 시작값 0.15/0.85, 선택한 쪽만 "
        "업데이트. 우리는 logistic(b0 + β·Q + 세션 내 시간 + 처음 3분), β 부호 자유, day 3에 가치 0, drop마다 업데이트 "
        "(operant 과제에서는 drop이 lick 뒤에만 나오기 때문).\n"
        f"[가치 Q] 'spout에서 기대하는 약물', 날이 바뀌어도 이어진다. Passive의 replay drop(day 6–10)도 마셨으므로 업데이트에 "
        f"넣었지만 그날 lick은 피팅하지 않는다. 피팅은 day 11–18 (마우스당 {nb} bin): 두 그룹이 같은 PR 일정과 같은 결과 순서를 "
        f"겪는 기간이다.\n"
        f"[β의 부호] Active {na}마리 모두 양수, Passive {npv}마리 모두 음수 (마우스별 피팅). Passive는 모르핀을 기대할수록 덜 "
        "핥으므로, 이들에게 day 14 (모르핀 → 물)는 보상 상실이 아니라 오히려 해방에 가깝다. 즉 논문의 H→L 비유는 Active에만 "
        "맞는다.\n"
        "[오른쪽 아래] 실제 한 세션 (Active 0910 black, day 14): 회색 칸 = lick이 있었던 30초 bin, ▼ = 물 drop. 모델은 bin마다 "
        "'lick이 있을 확률'을 예측한다.")


def slide3(prs, F):
    b.SLIDE = 3
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    a14, p14 = M11["day14_vs_13"]["Active"]["whole"], M11["day14_vs_13"]["Passive"]["whole"]
    a15, a16 = M11["day15_vs_13"]["Active"]["whole"], M11["day16_vs_13"]["Active"]["whole"]
    a17, p17 = M11["day17_vs_16"]["Active"]["whole"], M11["day17_vs_16"]["Passive"]["whole"]
    w14 = M11["day14_vs_13"]["Active_vs_Passive_whole"]["welch_p"]
    asy = M11["asymmetry_whole_session_abs_change_day17_vs_day14"]["Active"]
    v = dict(a14=reg("v2s3_a14", sf(a14["mean"]), a14["mean"]), a14lo=reg("v2s3_a14lo", sf(a14["lo"]), a14["lo"]),
             a14hi=reg("v2s3_a14hi", sf(a14["hi"]), a14["hi"]), p14=reg("v2s3_p14", sf(p14["mean"]), p14["mean"]),
             a15=reg("v2s3_a15", sf(a15["mean"]), a15["mean"]), a16=reg("v2s3_a16", sf(a16["mean"]), a16["mean"]),
             a17=reg("v2s3_a17", sf(a17["mean"]), a17["mean"]), p17=reg("v2s3_p17", sf(p17["mean"]), p17["mean"]),
             asy14=reg("v2s3_asy14", f"{asy['abs_day14']:.2f}", asy["abs_day14"]),
             asy17=reg("v2s3_asy17", f"{asy['abs_day17']:.2f}", asy["abs_day17"]))
    pv = dict(a14=fp(a14["p"]), p14=fp(p14["p"]), a15=fp(a15["p"]), a16=fp(a16["p"]), a17=fp(a17["p"]), p17=fp(p17["p"]),
              w14=fp(w14), asy=fp(asy["p_t"]))
    for k_, val in (("a14", a14["p"]), ("p14", p14["p"]), ("a15", a15["p"]), ("a16", a16["p"]), ("a17", a17["p"]),
                    ("p17", p17["p"]), ("w14", w14), ("asy", asy["p_t"])):
        reg(f"v2s3_p_{k_}", pv[k_], val)
    pe = {k_: (f"p {x}" if x.startswith("<") else f"p = {x}") for k_, x in pv.items()}
    ks = 13.5 if b.LANG == "EN" else 12.5
    key_box(s, box, tg, [P(Rn(T(f"First water day: Active {v['a14']} (95% CI {v['a14lo']} to {v['a14hi']}, n = 6) vs "
                                f"Passive {v['p14']} (n = 7); Welch t {pe['w14']}.",
                                f"물로 바뀐 첫날: Active {v['a14']} (95% CI {v['a14lo']} ~ {v['a14hi']}, n = 6) vs "
                                f"Passive {v['p14']} (n = 7); Welch t {pe['w14']}."), ks, True, "1F4E79"))],
            tip=tg["tip"])
    common(s, T("First water day: Active licking holds, Passive licking rises", "물로 바뀐 첫날: Active의 lick은 유지, Passive는 바로 증가"),
           T("compared each day with the last day before the switch (session means, within mice).",
             "전환 전 마지막 날과 각 날을 비교함 (세션 평균, 같은 마우스 안에서)."),
           T(f"Active fell only on days 15–16 ({v['a15']}, {v['a16']}; {pe['a15']}, {pe['a16']}) and rose on the first morphine "
             f"day ({v['a17']}, {pe['a17']}); Passive changed on the first day both ways ({v['p14']}, {v['p17']}).",
             f"Active는 day 15–16에야 감소 ({v['a15']}, {v['a16']}; {pe['a15']}, {pe['a16']})했고 모르핀 첫날 바로 증가 "
             f"({v['a17']}, {pe['a17']}); Passive는 양방향 모두 첫날 바뀜 ({v['p14']}, {v['p17']})."),
           T(f"p(lick) = share of 30-s bins with ≥ 1 lick; change = per-mouse session mean minus the reference day (13 or 16). "
             f"n = 6 Active, 7 Passive (6099 white died on day 13). One-sample t vs 0 and Welch t (Active vs Passive), two-sided. "
             f"Active, size of the first-day change, day 14 vs day 17: {v['asy14']} vs {v['asy17']}, paired t {pe['asy']}: the "
             f"asymmetry is not significant on its own. Drug action within the session and withdrawal can also give slow "
             f"stopping and fast resumption.",
             f"p(lick) = lick이 1번 이상 있었던 30초 bin의 비율; 변화 = 마우스별 세션 평균 − 기준일(13 또는 16). n = Active 6, "
             f"Passive 7 (6099 white는 day 13 사망). one-sample t (vs 0)와 Welch t (Active vs Passive), 양측. Active 첫날 변화 "
             f"크기, day 14 vs day 17: {v['asy14']} vs {v['asy17']}, paired t {pe['asy']}: 비대칭 자체는 유의하지 않음. 세션 중 "
             f"약효나 금단도 '늦게 멈추고 빨리 재개'를 만들 수 있음."),
           T(notes3_en(v, pe), notes3_kr(v, pe)))
    return s


def notes3_en(v, pe):
    return (
        "WHAT: behaviour without any model. a: p(lick) per 3-min block of each 15-min PR session, days 11-18; light blue = water "
        "days. b: day 14 minus day 13 (morphine -> water); c: day 17 minus day 16 (water -> morphine), per 3-min block; the "
        "printed tests use whole-session means.\n"
        f"Day 14: Active {v['a14']} (95% CI {v['a14lo']} to {v['a14hi']}, {pe['a14']}): no detectable change on the first water "
        f"day; Passive {v['p14']} ({pe['p14']}); Active vs Passive, Welch t {pe['w14']}.\n"
        f"Days 15-16 vs 13: Active {v['a15']} ({pe['a15']}) and {v['a16']} ({pe['a16']}): Active licking falls one day after the "
        f"switch. Day 17: Active {v['a17']} ({pe['a17']}), Passive {v['p17']} ({pe['p17']}).\n"
        "MEANING: Active mice changed one day late when morphine was removed but on the first day when it returned - the same "
        "direction as the paper's 'H->L slower than L->H'. Passive mice changed on the first day both ways (possibly taste: "
        "morphine is bitter; and beta < 0, so for them the valence is reversed).\n"
        f"CAVEATS: the direct test of the asymmetry in Active (size of the first-day change, {v['asy14']} vs {v['asy17']}, "
        f"{pe['asy']}, n = 6) is not significant. Drug action within the session (day 17) and a withdrawal state (day 14) could "
        "produce the same pattern; behaviour alone cannot separate these from learning.")


def notes3_kr(v, pe):
    return (
        "[이 슬라이드] 모델 없이 본 행동.\n"
        "[a 읽는 법] 15분 세션을 3분 블록 5개로 나눠 p(lick) (lick이 있었던 30초 bin의 비율)을 day 11–18에 걸쳐 그렸다. 빨강 = "
        "Active, 파랑 = Passive, 연파랑 배경 = 물 날.\n"
        "[b, c 읽는 법] b: day 14 − day 13 (모르핀 → 물), c: day 17 − day 16 (물 → 모르핀), 3분 블록별. 그림 안의 검정은 세션 "
        "전체 평균으로 했다.\n"
        f"[숫자] day 14: Active {v['a14']} (95% CI {v['a14lo']} ~ {v['a14hi']}, {pe['a14']}) → 물로 바뀐 첫날 변화가 검출되지 "
        f"않음; Passive {v['p14']} ({pe['p14']}); Active vs Passive Welch t {pe['w14']}. day 15–16 (vs day 13): Active {v['a15']} "
        f"({pe['a15']}), {v['a16']} ({pe['a16']}) → 하루 늦게 감소. day 17: Active {v['a17']} ({pe['a17']}), Passive {v['p17']} "
        f"({pe['p17']}).\n"
        "[의미] Active: 모르핀을 뺐을 때는 하루 늦게, 돌려줬을 때는 첫날 바뀐다 = 논문의 'H→L이 L→H보다 느리다'와 같은 방향. "
        "Passive: 양쪽 다 첫날 바뀐다 (학습보다 맛일 수 있음 — 모르핀은 쓰다; 또 β < 0이라 이들에게는 좋고 나쁨이 반대).\n"
        f"[주의] Active 안에서 두 첫날 변화의 크기를 직접 비교하면 {v['asy14']} vs {v['asy17']}, {pe['asy']} (n = 6)로 유의하지 "
        "않다. day 17의 빠른 회복은 세션 중 마신 모르핀의 약효, day 14에 그대로인 것은 금단 상태 때문일 수도 있다 — 행동만으로는 "
        "학습과 구분할 수 없다.")


def slide4(prs, F):
    b.SLIDE = 4
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    d = {f"{g}_{m}": reg(f"v2s4_{g}_{m}", f"{G5[g][m]['dBIC_from_best']:.0f}", G5[g][m]["dBIC_from_best"])
         for g in ("Active", "Passive") for m in ("RPE2a", "mRPE", "RPE", "today's outcome", "eMack", "RPEt-1")}
    a_to = d["Active_today's outcome"]
    phmin = min(G5[g][m]["dBIC_from_best"] for g in ("Active", "Passive") for m in ("PH", "no learning"))
    check(f"v2: PH / no learning dBIC >= 70 (min {phmin:.1f})", phmin >= 70)
    cva = CV["Active"]["delta_from_best"]
    check("v2: LOMO best = RPE2a (Active), today's outcome (Passive)",
          CV["Active"]["best"] == "RPE2a" and CV["Passive"]["best"] == "today's outcome")
    c_m = reg("v2s4_cv_m", f"{cva['mRPE']:.1f}", cva["mRPE"])
    c_r = reg("v2s4_cv_r", f"{cva['RPE']:.1f}", cva["RPE"])
    c_t = reg("v2s4_cv_t", f"{cva[chr(116) + 'oday' + chr(39) + 's outcome']:.1f}", cva["today's outcome"])
    cvp = CV["Passive"]["delta_from_best"]
    pmax = max(v_ for k_, v_ in cvp.items() if k_ not in ("today's outcome", "no learning"))
    c_p = reg("v2s4_cv_p", f"{pmax:.1f}", pmax)
    nA, nP = G5["Active"]["RPE"]["n"], G5["Passive"]["RPE"]["n"]
    am = reg("v2s4_am", f"{G5['Active']['RPE2a']['params'][0]:.2f}", G5["Active"]["RPE2a"]["params"][0])
    ks = 13.5 if b.LANG == "EN" else 12.5
    key_box(s, box, tg, [P(Rn(T("Best fit (BIC and held-out mice): two learning rates, slow after morphine → water, fast "
                                "after water → morphine.",
                                "최적 모델 (BIC와 held-out 마우스): 학습률 2개, 모르핀 → 물 뒤엔 느리게, 물 → 모르핀 "
                                "뒤엔 빠르게."), ks, True, "1F4E79"))], tip=tg["tip"])
    common(s, T("Active: slow updating after morphine loss fits best", "Active: 모르핀이 빠진 뒤 느리게 업데이트하는 모델이 최적"),
           T("fitted 6 learning rules + 2 no-learning controls (learning rates shared within a group); BIC and "
             "leave-one-mouse-out CV.",
             "학습 규칙 6개 + 대조 2개를 그룹 단위로 피팅; BIC, leave-one-mouse-out."),
           T(f"Active: two rates best (ΔBIC meta-RPE +{d['Active_mRPE']}, one rate +{d['Active_RPE']}, instant switch "
             f"+{a_to}; held-out negLL +{c_m}, +{c_r}, +{c_t}); Passive: learning rules ≈ instant switch (held-out ≤ {c_p}).",
             f"Active: 학습률 2개가 최선 (ΔBIC meta-RPE +{d['Active_mRPE']}, 학습률 1개 +{d['Active_RPE']}, 즉시 전환 "
             f"+{a_to}; held-out +{c_m}, +{c_r}, +{c_t}); Passive: 학습 규칙 ≈ 즉시 전환 (held-out ≤ {c_p})."),
           T(f"Lines = mean model predictions; dots = observed mean ± SEM (3-min blocks). Readout logistic(b0 + β·Q + "
             f"time + first 3 min) per mouse. BIC (n = {nA} / {nP} bins) treats bins as independent and overstates the evidence; "
             f"held-out = learning rates fitted on the other mice. n = 6 Active, one switch each way; α− at its lower "
             f"bound ({am}); withdrawal or drug action within the session can mimic the pattern. PH, no learning: ΔBIC ≥ 70.",
             f"선 = 모델 예측 평균; 점 = 실제 평균 ± SEM (3분 블록). 행동 출력 logistic(b0 + β·Q + 시간 + 처음 "
             f"3분)은 마우스별. BIC (n = {nA} / {nP} bin)는 bin을 독립으로 봐서 증거를 과대평가; held-out = 다른 마우스로 학습률을 "
             f"피팅. Active 6마리, 방향마다 전환 1번; α−는 하한 ({am}); 금단이나 세션 중 약효도 같은 패턴을 만들 수 있음. "
             f"PH, no learning: ΔBIC ≥ 70."),
           T(notes4_en(d, nA, nP, c_m, c_r, c_t, c_p), notes4_kr(d, nA, nP, c_m, c_r, c_t, c_p)))
    return s


def _gr(g, gen, pk):
    return pct(GR["confusion"][g][gen].get(pk, 0.0))


def notes4_en(d, nA, nP, c_m, c_r, c_t, c_p):
    pa = G5["Active"]["mRPE"]["params"]
    p2 = G5["Active"]["RPE2a"]["params"]
    a_to = d["Active_today's outcome"]
    nl = G5["Active"]["mRPE"]["negLL"], G5["Active"]["RPE2a"]["negLL"]
    return (
        "WHAT: the learning rules fitted per group: learning rates shared by the mice of a group, readout (b0, beta, time terms) "
        f"per mouse, solved exactly for every candidate learning rate. BIC with n = {nA} (Active) / {nP} (Passive) bins.\n"
        "a: dots = observed, lines = mean model predictions. Orange (RPE 2alpha) and green (meta-RPE) fall slowly on day 14 and "
        "rise fast on day 17; dashed purple (today's outcome = switch at once) predicts day 14 far too low. b: in Passive all "
        "four model lines overlap.\n"
        f"c: dBIC. Active: RPE 2alpha 0, meta-RPE {d['Active_mRPE']}, RPE {d['Active_RPE']}, RPE t-1 {d['Active_RPEt-1']}, eMack "
        f"{d['Active_eMack']}, today's outcome {a_to}. Passive: today's outcome 0, RPE {d['Passive_RPE']}, RPE t-1 "
        f"{d['Passive_RPEt-1']}, RPE 2alpha {d['Passive_RPE2a']}, meta-RPE {d['Passive_mRPE']}. Because 30-s bins are not "
        "independent, these dBIC values overstate the evidence.\n"
        f"HELD-OUT MICE (leave-one-mouse-out): Active total held-out negLL relative to RPE 2alpha: meta-RPE +{c_m}, RPE +{c_r}, "
        f"today's outcome +{c_t} (RPE 2alpha better than RPE in 5/6 held-out mice). Passive: every learning rule within {c_p} of "
        "today's outcome - indistinguishable.\n"
        f"GROUP-LEVEL RECOVERY (50 simulated sets per model): RPE 2alpha was picked in {_gr('Active', 'RPE', 'RPE2a')}% of "
        f"RPE-made and {_gr('Active', chr(116) + 'oday' + chr(39) + 's outcome', 'RPE2a')}% of instant-switch-made sets, so its win "
        "is not a false positive of BIC; meta-RPE-made sets were called RPE 2alpha in "
        f"{_gr('Active', 'mRPE', 'RPE2a')}% (meta-RPE itself {_gr('Active', 'mRPE', 'mRPE')}%).\n"
        f"Parameters (Active): RPE 2alpha alpha- = {p2[0]:.2f} (lower bound), alpha+ = {p2[1]:.2f}; meta-RPE alpha = {pa[0]:.2f}, "
        f"d = {pa[1]:.2f}, e = {pa[2]:.2f} (both at the paper's bounds). meta-RPE also fits slightly worse (negLL {nl[0]:.1f} vs "
        f"{nl[1]:.1f}) with one parameter more.\n"
        "MEANING: Active licking is best described by slow updating after morphine loss and fast updating after its return; "
        "whether errors are averaged (the meta part) cannot be told. Passive mice switch without a detectable learning lag. "
        "All Passive mice have beta < 0 (they lick less when expecting morphine).")


def notes4_kr(d, nA, nP, c_m, c_r, c_t, c_p):
    pa = G5["Active"]["mRPE"]["params"]
    p2 = G5["Active"]["RPE2a"]["params"]
    a_to = d["Active_today's outcome"]
    nl = G5["Active"]["mRPE"]["negLL"], G5["Active"]["RPE2a"]["negLL"]
    tdo = "today's outcome"
    return (
        "[이 슬라이드] 학습 규칙을 그룹 단위로 피팅한 결과.\n"
        f"[방법] 학습률은 그룹 안의 마우스들이 공유, 행동 출력(b0, β, 시간 항)은 마우스별. 후보 학습률마다 행동 출력은 로지스틱 "
        f"회귀로 정확히 풀었다. BIC의 n = Active {nA}, Passive {nP} bin.\n"
        "[a, b 읽는 법] 점 = 실제, 선 = 모델 예측 평균. 주황(RPE 2α)과 초록(meta-RPE)은 day 14에 천천히 내려가고 day 17에 빨리 "
        "올라간다. 보라 점선(today's outcome = 그날 결과로 바로 전환)은 day 14를 너무 낮게 예측한다. Passive(b)에서는 네 모델 선이 "
        "모두 겹친다.\n"
        f"[c 읽는 법] ΔBIC. Active: RPE 2α 0, meta-RPE {d['Active_mRPE']}, RPE {d['Active_RPE']}, RPE t−1 {d['Active_RPEt-1']}, "
        f"eMack {d['Active_eMack']}, today's outcome {a_to}. Passive: today's outcome 0, RPE {d['Passive_RPE']}, RPE t−1 "
        f"{d['Passive_RPEt-1']}, RPE 2α {d['Passive_RPE2a']}, meta-RPE {d['Passive_mRPE']}. 30초 bin은 서로 독립이 아니라서 이 "
        "ΔBIC는 증거를 과대평가한다.\n"
        f"[held-out 마우스] leave-one-mouse-out (학습률은 나머지 마우스로 피팅): Active held-out negLL 합, RPE 2α 대비 meta-RPE "
        f"+{c_m}, RPE +{c_r}, today's outcome +{c_t} (held-out 6마리 중 5마리에서 RPE 2α가 RPE보다 좋음). Passive: 모든 학습 규칙이 "
        f"today's outcome과 {c_p} 이내 → 구분 불가.\n"
        f"[그룹 단위 복원] 모델마다 50세트 시뮬레이션: RPE로 만든 데이터에서 RPE 2α가 뽑힌 비율 {_gr('Active', 'RPE', 'RPE2a')}%, "
        f"즉시 전환 데이터에서 {_gr('Active', tdo, 'RPE2a')}% → RPE 2α의 승리는 BIC의 거짓 양성이 아님. meta-RPE로 만든 데이터는 "
        f"{_gr('Active', 'mRPE', 'RPE2a')}%가 RPE 2α로 분류됨 (meta-RPE 자신은 {_gr('Active', 'mRPE', 'mRPE')}%).\n"
        f"[파라미터] Active RPE 2α: α− = {p2[0]:.2f} (하한), α+ = {p2[1]:.2f}. meta-RPE: α = {pa[0]:.2f}, d = {pa[1]:.2f}, "
        f"e = {pa[2]:.2f} (둘 다 논문 범위 끝). meta-RPE는 같은 비대칭을 보이지만 우도도 약간 나쁘고 (negLL {nl[0]:.1f} vs "
        f"{nl[1]:.1f}) 파라미터가 하나 더 많다.\n"
        "[의미] Active의 lick은 '모르핀이 빠지면 느리게, 돌아오면 빠르게 업데이트'로 가장 잘 설명된다. 오류를 평균하는지(meta "
        "부분)는 알 수 없다. Passive는 검출 가능한 학습 지연 없이 바뀐다. Passive는 모두 β < 0 (모르핀을 기대할수록 덜 핥음).")


def slide5(prs, F):
    b.SLIDE = 5
    path, size, tg = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    conf = REC["confusion_BIC"]
    va, vp = conf["Active"]["mRPE"].get("mRPE", 0), conf["Passive"]["mRPE"].get("mRPE", 0)
    ra, rp = reg("v2s5_ra", pct(va), va * 100), reg("v2s5_rp", pct(vp), vp * 100)
    gm = GR["confusion"]["Active"]["mRPE"].get("mRPE", 0.0)
    g_m = reg("v2s5_gm", pct(gm), gm * 100)
    fp_rpe2a = max(GR["confusion"]["Active"]["RPE"].get("RPE2a", 0.0), GR["confusion"]["Active"]["today's outcome"].get("RPE2a", 0.0))
    g_fp = reg("v2s5_gfp", pct(fp_rpe2a), fp_rpe2a * 100)
    A = BS["Active"]
    am = reg("v2s5_am", f"{A['full_fit']['a_minus']:.2f}", A["full_fit"]["a_minus"])
    ap = reg("v2s5_ap", f"{A['full_fit']['a_plus']:.2f}", A["full_fit"]["a_plus"])
    fa = reg("v2s5_fa", f"{A['frac_a_minus_less_than_a_plus'] * 100:.1f}", A["frac_a_minus_less_than_a_plus"] * 100)
    pl = [GR["confusion"]["Passive"][gen].get("today's outcome", 0.0) for gen in ("RPE", "RPE2a", "mRPE")]
    plo, phi = reg("v2s5_plo", pct(min(pl)), min(pl) * 100), reg("v2s5_phi", pct(max(pl)), max(pl) * 100)
    ks = 12.5 if b.LANG == "EN" else 11.5
    key_box(s, box, tg, [P(Rn(T(f"meta-RPE made the data, yet per-mouse BIC found it in only {ra}% (Active) / {rp}% "
                                f"(Passive) of sets.",
                                f"meta-RPE로 만든 데이터에서도 마우스별 BIC가 찾은 비율은 {ra}% (Active) / {rp}% (Passive)뿐."),
                              ks, True, "1F4E79"))],
            tip=tg["tip"])
    common(s, T("The data show the asymmetry but cannot test error averaging", "비대칭은 보이지만, 오류 '평균'이 필요한지는 검정할 수 없다"),
           T("refitted simulated licking per mouse (20 sets each) and per group (50 sets each); bootstrapped the group learning "
             "rates.",
             "시뮬레이션 데이터를 마우스별(20세트)·그룹별(50세트)로 다시 피팅; 학습률은 bootstrap."),
           T(f"meta-RPE is not recoverable ({ra}% per mouse, {g_m}% per group), but two rates are never picked for one-rate or "
             f"instant-switch data ({g_fp}%); Active α− {am} vs α+ {ap} (α− < α+ in {fa}% of "
             f"bootstraps).",
             f"meta-RPE는 복원 안 됨 (마우스별 {ra}%, 그룹 {g_m}%); 학습률 1개·즉시 전환 데이터가 학습률 2개로 분류된 적 없음 "
             f"({g_fp}%); Active α− {am} vs α+ {ap} ({fa}%)."),
           T(f"a: each mouse's own fitted model + its real drop sequence; 3 models refitted; 6 + 7 mice x 20 sets. b: shared "
             f"learning rates + each mouse's readout; 4 models refitted; 50 sets per model (Passive: sets made with learning are "
             f"called 'instant switch' in {plo}–{phi}%). c: point = all mice, bar = 95% CI of 1000 bootstraps over mice; bounds "
             f"0.01–1. A per-mouse test needs several switches per mouse (e.g. repeated morphine ↔ water, or 85% → 15% "
             f"morphine drops).",
             f"a: 마우스별 모델 + 실제 drop 순서, 후보 3개, 6 + 7마리 x 20세트. b: 공유 학습률 + "
             f"마우스별 행동 출력; 후보 4개; 모델마다 50세트 (Passive: 학습 데이터도 {plo}–{phi}%가 '즉시 전환'으로 "
             f"분류). c: 점 = 전체 마우스 피팅, 막대 = 마우스 1000번 bootstrap의 95% CI; 범위 0.01–1. 마우스별 검정에는 "
             f"마우스당 전환이 여러 번 필요 (예: 모르핀 ↔ 물 반복, 또는 모르핀 drop 85% → 15%)."),
           T(notes5_en(ra, rp, g_m, g_fp, am, ap, fa, plo, phi), notes5_kr(ra, rp, g_m, g_fp, am, ap, fa, plo, phi)))
    return s


def notes5_en(ra, rp, g_m, g_fp, am, ap, fa, plo, phi):
    conf = REC["confusion_BIC"]
    A, Pv = BS["Active"], BS["Passive"]
    tdo = "today's outcome"
    return (
        "WHAT: can these data tell the models apart, and how sure is the learning asymmetry?\n"
        f"a (per mouse, like the paper): from each mouse's fitted models we simulated 20 licking sets and refitted three models. "
        f"When meta-RPE made the data, BIC picked meta-RPE in {ra}% (Active) / {rp}% (Passive) of sets; RPE-made sets were "
        f"recognised in {pct(conf['Active']['RPE'].get('RPE', 0))}% / {pct(conf['Passive']['RPE'].get('RPE', 0))}%, instant-switch "
        f"sets in {pct(conf['Active'][tdo].get(tdo, 0))}% / {pct(conf['Passive'][tdo].get(tdo, 0))}%. Per-mouse BIC leans to the "
        "simplest model, so 'meta-RPE lost' per mouse means nothing here.\n"
        f"b (group level, Active): RPE 2alpha-made sets were recovered in {_gr('Active', 'RPE2a', 'RPE2a')}% (otherwise called RPE); "
        f"RPE- and instant-switch-made sets were never called RPE 2alpha ({g_fp}%); meta-RPE-made sets were recovered in only "
        f"{g_m}% (mostly called RPE 2alpha). Passive: sets made with learning were called instant switch in {plo}-{phi}%, so the "
        "Passive 'instant switch' result does not exclude fast learning.\n"
        f"c: RPE 2alpha group fit, 1000 bootstraps over mice. Active alpha- {am} [{A['a_minus_CI'][0]:.2f}, {A['a_minus_CI'][2]:.2f}] "
        f"vs alpha+ {ap} [{A['a_plus_CI'][0]:.2f}, {A['a_plus_CI'][2]:.2f}]; alpha- < alpha+ in {fa}%. Passive alpha- "
        f"{Pv['full_fit']['a_minus']:.2f} vs alpha+ {Pv['full_fit']['a_plus']:.2f}, "
        f"{Pv['frac_a_minus_less_than_a_plus'] * 100:.0f}% (unclear).\n"
        "MEANING: slow updating after morphine loss / fast after return is the robust part in Active mice; whether errors are "
        "averaged (the meta part) cannot be decided with this design. NEXT DESIGN: several switches per mouse, probabilistic "
        "outcomes (e.g. 85% -> 15% morphine drops) or a cue that signals morphine vs water.")


def notes5_kr(ra, rp, g_m, g_fp, am, ap, fa, plo, phi):
    conf = REC["confusion_BIC"]
    A, Pv = BS["Active"], BS["Passive"]
    tdo = "today's outcome"
    return (
        "[이 슬라이드] 이 데이터로 모델을 구분할 수 있는가, 그리고 학습률 비대칭은 얼마나 확실한가.\n"
        f"[a 읽는 법] (마우스별, 논문 방식) 각 마우스의 피팅 모델로 시뮬레이션 데이터를 20세트씩 만들고 후보 3개를 다시 피팅했다. "
        f"meta-RPE로 만든 데이터인데도 BIC가 meta-RPE를 고른 비율은 {ra}% (Active), {rp}% (Passive). RPE로 만든 데이터는 "
        f"{pct(conf['Active']['RPE'].get('RPE', 0))}% / {pct(conf['Passive']['RPE'].get('RPE', 0))}%, 즉시 전환으로 만든 데이터는 "
        f"{pct(conf['Active'][tdo].get(tdo, 0))}% / {pct(conf['Passive'][tdo].get(tdo, 0))}% 맞게 찾았다. 마우스별 BIC는 가장 단순한 "
        "모델 쪽으로 기울기 때문에, 마우스별 비교에서 'meta-RPE가 졌다'는 여기서는 의미가 없다.\n"
        f"[b 읽는 법] (그룹 단위, Active) RPE 2α로 만든 데이터는 {_gr('Active', 'RPE2a', 'RPE2a')}% 복원 (나머지는 RPE로 분류). RPE나 "
        f"즉시 전환으로 만든 데이터가 RPE 2α로 분류된 적은 없음 ({g_fp}%). meta-RPE로 만든 데이터는 {g_m}%만 복원 (대부분 RPE 2α로 "
        f"분류). Passive: 학습으로 만든 데이터도 {plo}–{phi}%는 '즉시 전환'으로 분류되므로, Passive의 '즉시 전환' 결과가 빠른 학습을 "
        "배제하지는 못한다.\n"
        f"[c 읽는 법] RPE 2α 그룹 피팅의 학습률과 마우스 1000번 bootstrap의 95% CI. Active α−(모르핀 → 물) {am} "
        f"[{A['a_minus_CI'][0]:.2f}, {A['a_minus_CI'][2]:.2f}] vs α+(물 → 모르핀) {ap} [{A['a_plus_CI'][0]:.2f}, "
        f"{A['a_plus_CI'][2]:.2f}]; {fa}%에서 α− < α+. Passive α− {Pv['full_fit']['a_minus']:.2f} vs α+ "
        f"{Pv['full_fit']['a_plus']:.2f}, {Pv['frac_a_minus_less_than_a_plus'] * 100:.0f}% (불확실).\n"
        "[의미] Active에서 확실한 부분은 '모르핀이 빠지면 느리게, 돌아오면 빠르게 업데이트'이다. 오류를 평균하는지(meta 부분)는 이 "
        "설계로는 판단할 수 없다. [다음 설계] 마우스마다 전환 여러 번, 확률적 결과 (예: 모르핀 drop 85% → 15%), 또는 모르핀/물을 "
        "미리 알려 주는 큐가 있으면 구분할 수 있다.")


def slide6(prs, F):
    b.SLIDE = 6
    path, size, tg, info = F
    s = b.new_slide(prs)
    box = place(s, path, size, 0.25, 0.98, 9.25)
    ia, ip = info["Active"], info["Passive"]
    na = reg("v2s6_na", str(ia["cd_min_after"]), ia["cd_min_after"])
    nw = reg("v2s6_nw", str(ia["n_water"]), ia["n_water"])
    npv = reg("v2s6_np", str(ip["cd_min_after"]), ip["cd_min_after"])
    ra, rp = info["range"]["Active"], info["range"]["Passive"]
    rtxt_a = f"{ra[0]}–{ra[1]}" if ra[0] != ra[1] else f"{ra[0]}"
    rtxt_p = f"{rp[0]}–{rp[1]}" if rp[0] != rp[1] else f"{rp[0]}"
    ks = 12.5 if b.LANG == "EN" else 11.5
    paras = [P(Rn(T("If CD cells exist (paper: mPFC), meta-RPE predicts", "CD 세포가 있다면 (논문: mPFC), meta-RPE의 예측"),
                  ks + 1.5, True, "1F4E79"), after=6),
             P(Rn("Active: ", ks, True, "D62728"),
               Rn(T(f"slow ramp; CD lowest after {na} of {nw} water drops (days 14–16).",
                    f"천천히 ramp; 물 drop {nw}개 중 {na}번째에 CD 최소 (day 14–16)."), ks), after=6),
             P(Rn("Passive: ", ks, True, "1F77B4"),
               Rn(T(f"fast step; CD lowest {npv} drops after the switch (β < 0: for them morphine → water is not an "
                    f"omission).",
                    f"빠른 계단; 전환 후 {npv}번째 drop에 CD 최소 (β < 0: 이들에게 모르핀 → 물은 보상 상실이 아님)."),
                  ks), after=6),
             P(Rn(T("Test as in the paper: GLM with δ, CD, CE and Q per drop, circular-shift null.",
                    "검정은 논문처럼: drop마다 δ, CD, CE, Q를 넣은 GLM, circular shift로 p값."), ks, False, "444444"))]
    kx, ky, kw = 9.62, 1.30, 3.45
    tbox(s, kx, ky, kw, None, paras, "KEY", fill="F2F6FB", line="1F4E79")
    h = measure(paras, kw) + 0.06
    check(f"[{b.LANG}] v2 slide 6 KEY ends above the bottom line", ky + h <= 5.85)
    tops = np.cumsum([0.05] + [measure([pp], kw) - 0.10 for pp in paras])
    for g, k in (("Active", 1), ("Passive", 2)):
        tx, ty = at(box, tg[g])
        ay = ky + (tops[k] + tops[k + 1]) / 2 if g == "Active" else ky + h - 0.15
        arrow(s, kx, ay, tx + 0.06, ty)
    pa, pp_ = G5["Active"]["mRPE"]["params"], G5["Passive"]["mRPE"]["params"]
    dA, dP = G5["Active"]["mRPE"]["dBIC_from_best"], G5["Passive"]["mRPE"]["dBIC_from_best"]
    rd = REC["corr_d"]
    common(s, T("Model variables per drop: regressors for ORBm / BMAp neurons", "drop마다 나오는 모델 변수: ORBm / BMAp 뉴런 회귀변수 후보"),
           T("computed RPE δ, CD, CE, gain M and value Q at every drop with the group-fitted meta-RPE (one yoked pair "
             "shown).",
             "그룹에 피팅한 meta-RPE로 drop마다 RPE δ, CD, CE, gain M, 가치 Q를 계산 (yoked pair 하나)."),
           T(f"Active CD is lowest {na} drops after the switch, Passive after {npv} (shared parameters: Active {rtxt_a}, Passive "
             f"{rtxt_p} across mice): a prediction, not a result.",
             f"Active CD는 전환 후 {na}번째, Passive는 {npv}번째 drop에서 최소 (공유 파라미터라 마우스 간 Active {rtxt_a}, "
             f"Passive {rtxt_p}): 결과가 아닌 예측."),
           T(f"Pair 4: Active 0910 black, Passive 0910 orange. meta-RPE lost the group comparison (ΔBIC +{dA:.0f} Active, "
             f"+{dP:.0f} Passive); d = {pa[1]:.2f} and e = {pa[2]:.2f} sit at the paper's bounds and d is poorly recovered "
             f"(r = {rd:.2f}), so the CD time course is a model assumption. Only a negative CD or a positive CE raises the gain M. "
             f"Paper: CD cells in mPFC (prelimbic) and mPFC→VTA; ORBm and BMAp were not studied.",
             f"Pair 4: Active 0910 black, Passive 0910 orange. meta-RPE는 그룹 비교에서 짐 (ΔBIC Active +{dA:.0f}, Passive "
             f"+{dP:.0f}); d = {pa[1]:.2f}, e = {pa[2]:.2f}는 논문 범위 끝이고 d는 복원이 안 됨 (r = {rd:.2f}) → CD 시간 "
             f"경과는 모델 가정. CD < 0 이거나 CE > 0 일 때만 gain M이 커짐. 논문: mPFC(prelimbic)와 mPFC→VTA의 CD 세포; "
             f"ORBm, BMAp는 연구되지 않음."),
           T(notes6_en(na, nw, npv, pa, pp_, rtxt_a, rtxt_p, dA, dP, rd), notes6_kr(na, nw, npv, pa, pp_, rtxt_a, rtxt_p, dA, dP, rd)))
    reg("v2s6_dA", f"{dA:.0f}", dA)
    reg("v2s6_dP", f"{dP:.0f}", dP)
    reg("v2s6_rd", f"{rd:.2f}", rd)
    return s


def notes6_en(na, nw, npv, pa, pp_, rtxt_a, rtxt_p, dA, dP, rd):
    return (
        "WHAT: the variables meta-RPE computes at every drop - what a neural GLM would use once ORBm / BMAp calcium imaging "
        "exists (as in the paper's Fig. 2, with a circular-shift null for p values).\n"
        "a, c: black = value Q (expected drug), grey bars = RPE delta per drop; from day 14 the drops are water, so delta < 0.\n"
        "b, d: red CD = slow running mean of errors, blue CE = fast running mean, green dashed = gain M (sets how much is "
        "learned). M rises when CD < 0 or CE > 0. CE also dips fast at the switch in both groups, but a negative CE does not raise "
        "M.\n"
        f"PREDICTION (only if meta-RPE is right): Active 0910 black: CD builds over days 14-16 and is lowest after {na} of {nw} "
        f"water drops -> a CD cell would ramp slowly. Passive 0910 orange: CD is lowest {npv} drops after the switch -> a step. "
        f"With shared group parameters the counts are nearly the same in every mouse (Active {rtxt_a}, Passive {rtxt_p}).\n"
        f"CAVEATS: meta-RPE lost the group comparison (dBIC +{dA:.0f} Active, +{dP:.0f} Passive); d and e sit at the paper's "
        f"bounds (Active alpha = {pa[0]:.2f}, d = {pa[1]:.2f}, e = {pa[2]:.2f}; Passive alpha = {pp_[0]:.2f}, d = {pp_[1]:.2f}, "
        f"e = {pp_[2]:.2f}) and d is poorly recovered (r = {rd:.2f}). The winning RPE 2alpha model gives value Q and delta but no "
        "averaged error. Passive mice have beta < 0, so 'drug gone' is not an omission for them. The paper found CD cells in "
        "mPFC (prelimbic; 13% of cells) and among mPFC->VTA cells (61%); ORBm and BMAp were not studied.")


def notes6_kr(na, nw, npv, pa, pp_, rtxt_a, rtxt_p, dA, dP, rd):
    return (
        "[이 슬라이드] meta-RPE가 drop마다 계산하는 변수들. 나중에 ORBm / BMAp 칼슘 이미징을 하면 GLM 회귀변수로 쓸 것이다 "
        "(논문 Fig. 2와 같은 방식, circular shift로 p값).\n"
        "[a, c 읽는 법] 검은 선 = 가치 Q (기대 약물), 회색 막대 = drop마다의 오류 δ. day 14부터 물 drop이라 δ가 음수.\n"
        "[b, d 읽는 법] 빨강 CD = 오류의 느린 이동평균, 파랑 CE = 빠른 이동평균, 초록 점선 = gain M (한 번에 얼마나 배울지 결정). "
        "CD < 0 이거나 CE > 0 일 때 M이 커진다. 전환 때 CE도 두 그룹 모두 빠르게 내려가지만, 음수 CE는 M을 키우지 않는다.\n"
        f"[예측, meta-RPE가 맞을 때만] Active (0910 black): CD가 day 14–16 동안 천천히 쌓여 물 drop {nw}개 중 {na}번째에 최소 → "
        f"'CD 세포'라면 천천히 ramp. Passive (0910 orange): 전환 후 {npv}번째 drop에서 최소 → 계단처럼 빠르게. 그룹 공유 파라미터라 "
        f"마우스마다 거의 같다 (Active {rtxt_a}, Passive {rtxt_p}).\n"
        f"[주의] meta-RPE는 그룹 비교에서 졌다 (ΔBIC Active +{dA:.0f}, Passive +{dP:.0f}). d, e는 논문 범위 끝에 붙어 있고 (Active "
        f"α = {pa[0]:.2f}, d = {pa[1]:.2f}, e = {pa[2]:.2f}; Passive α = {pp_[0]:.2f}, d = {pp_[1]:.2f}, e = {pp_[2]:.2f}) d는 복원이 "
        f"안 된다 (r = {rd:.2f}). 이긴 RPE 2α 모델은 가치 Q와 δ는 주지만 평균 오류는 없다. Passive는 β < 0이라 '약이 없어짐'이 보상 "
        "상실이 아니다. 논문은 mPFC(prelimbic, 뉴런의 13%)와 mPFC→VTA 뉴런(61%)에서 CD를 찾았고 ORBm, BMAp는 연구하지 않았다.")


# ============================================================================================ build
def main():
    PRED = b.group_predictions()
    FIGS = [fig_paper(), fig_map(), fig_behaviour(), fig_models(PRED), fig_recovery(), fig_vars()]
    BUILD = [slide1, slide2, slide3, slide4, slide5, slide6]
    NAMES = ["1_paper_fig1", "2_mapping", "3_behaviour_switches", "4_group_model_fits", "5_recovery_asymmetry",
             "6_model_variables"]
    decks = {}
    for lang in ("EN", "KR"):
        b.LANG = lang
        allp = Presentation(str(b.TEMPLATE))
        for fn, F, nm in zip(BUILD, FIGS, NAMES):
            one = Presentation(str(b.TEMPLATE))
            fn(one, F)
            one.save(OUT / "single_slides_v2" / f"{lang}_slide{nm}_v2_20261117.pptx")
            fn(allp, F)
        decks[lang] = OUT / f"mrpe_hjort2026_opioid_6slides_{lang}_v2_20261117.pptx"
        allp.save(decks[lang])
    for k, (txt, val) in b.PRINTED.items():
        if not k.startswith("v2"):
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
    bad = [c for c in b.CHECKS if not c[1]]
    print(f"{len(b.CHECKS)} checks, {len(bad)} failed")
    for c in bad:
        print("  FAILED", c[0])
    for lang, deck in decks.items():
        print(lang, deck, [p.name for p in b.render(deck, OUT / f"preview_{lang}_v2")])


if __name__ == "__main__":
    main()
