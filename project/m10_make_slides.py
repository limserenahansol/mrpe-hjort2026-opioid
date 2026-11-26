"""Step 10: the meta-RPE analysis as 6 slides, an English (EN) and a Korean (KR) deck, python-pptx only (never Office/COM).

Slides: 1 the paper's Fig. 1 (authors' data, our port) | 2 how the model maps onto our task | 3 behaviour at the two
switches (model-free) | 4 group model fits | 5 per-mouse recovery + learning asymmetry | 6 model variables per drop.
Every number printed on a slide comes from the m03-m08 outputs and is asserted (PRINTED). Text boxes are measured with
the real fonts (Arial; Malgun Gothic for KR) and must fit. KR deck: Korean slide text + a full Korean explanation of each
slide in the speaker notes (figures keep English labels). Previews: LibreOffice headless (own profile) -> PDF -> PNG.
Output (slides_20261117\\): mrpe_hjort2026_opioid_6slides_EN_20261117.pptx, ..._KR_20261117.pptx, preview_EN\\,
        preview_KR\\, single_slides\\, figs\\, fit_report.json
"""
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pymupdf  # noqa: E402
from lxml import etree  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402
from PIL import ImageFont  # noqa: E402
from pptx import Presentation  # noqa: E402
from pptx.dml.color import RGBColor  # noqa: E402
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE  # noqa: E402
from pptx.enum.text import MSO_ANCHOR, MSO_AUTO_SIZE, PP_ALIGN  # noqa: E402
from pptx.oxml.ns import qn  # noqa: E402
from pptx.util import Inches, Pt  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m02_hjort_models as hm  # noqa: E402
import m04_fit_user_data as m4  # noqa: E402
import m05_recovery_and_group_fits as m5  # noqa: E402

D = HERE / "data"
OUT = HERE / "slides_20261117"
FIG = OUT / "figs"
for p in (OUT, FIG, OUT / "single_slides"):
    p.mkdir(parents=True, exist_ok=True)
TEMPLATE = Path(r"C:\Users\hsollim\Research_Projects\CLAUDE_SCIENCE\labmeeting_NOV26_2026\round2"
                r"\_template_from_real_labmeeting_20261102.pptx")
SOFFICE = r"C:\Program Files\LibreOffice\program\soffice.exe"
LO_PROFILE = Path(tempfile.gettempdir()) / "lo_profile_mrpe_slides"

S3 = json.load(open(D / "m03_paper_data_summary.json"))
N7 = json.load(open(D / "m07_fig1_numbers.json"))
MF = json.load(open(D / "m06_modelfree.json"))
G5 = json.load(open(D / "m05_group_fits.json"))
REC = json.load(open(D / "m05_recovery.json"))
BS = json.load(open(D / "m08_bootstrap_asymmetry.json"))
X6 = pd.read_csv(D / "m06_switch_blocks.csv")

RED, BLUE, INK, GREY, MUT, NAVY = "#D62728", "#1F77B4", "#222222", "#8A8A8A", "#555555", "000090"
GREEN, ORANGE, PURPLE, WBG = "#2E7D32", "#FF7F0E", "#9467BD", "#E8F1FB"
CMOD = {"today's outcome": PURPLE, "RPE": GREY, "RPE2a": ORANGE, "mRPE": GREEN}
DPI = 200
mpl.rcParams.update({"font.family": "Arial", "font.size": 13, "axes.labelsize": 13.5, "xtick.labelsize": 12,
                     "ytick.labelsize": 12, "axes.titlesize": 13.5, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.edgecolor": "#444444", "axes.unicode_minus": False,
                     "legend.fontsize": 11, "figure.dpi": DPI, "savefig.dpi": DPI})
CHECKS, FIT = [], []
PRINTED = {}
LANG = "EN"
SLIDE = 0


def check(name, ok):
    CHECKS.append((name, bool(ok)))
    if not ok:
        print("CHECK FAILED:", name)


def T(en, kr):
    """slide text in the deck's language; Malgun Gothic has no U+2212 minus or U+2248."""
    return en if LANG == "EN" else kr.replace("\u2212", "-").replace("\u2248", "~")


def fp(p):
    if p < 0.0001:
        return "< 0.0001"
    if p < 0.001:
        return f"{p:.4f}"
    if p < 0.01:
        return f"{p:.3f}"
    if p < 0.1:
        return f"{p:.3f}"
    return f"{p:.2f}"


def peq(p):
    s = fp(p)
    return f"p {s}" if s.startswith("<") else f"p = {s}"


def sf(x, nd=2):
    if round(x, nd) == 0:
        return f"{0:.{nd}f}"
    return f"{x:+.{nd}f}".replace("-", "\u2212")


def reg(key, txt, val):
    """remember a printed number and the value it stands for (checked at the end)"""
    PRINTED[key] = (txt, float(val))
    return txt


def frac(fig, ax, x, y):
    fig.canvas.draw()
    X, Y = ax.transData.transform((x, y))
    W, H = fig.get_size_inches() * fig.dpi
    return X / W, 1 - Y / H


def save(fig, name):
    p = FIG / name
    fig.savefig(p, facecolor="white")
    plt.close(fig)
    return p


# ============================================================================================ data for the figures
def rw_q(p, c, r):
    Q = [0.15, 0.85]
    out = []
    for t in range(len(c)):
        out.append(Q[1])
        Q[c[t]] += p[0] * (r[t] - Q[c[t]])
    return np.array(out)


def paper_curves():
    F3 = pd.read_csv(D / "m03_paper_data_fits.csv")
    F3["params"] = F3.params.apply(json.loads)
    L, QR, QM, TM = [], [], [], []
    for a in sorted(F3.animal.unique()):
        c = np.load(D / "hjort2026_sample_choice" / a / "HL_choice.npy", allow_pickle=False).astype(np.int64)
        r = np.load(D / "hjort2026_sample_choice" / a / "HL_outcome.npy", allow_pickle=False).astype(np.float64)
        pm = np.array(F3[(F3.animal == a) & (F3.model == "mRPE")].params.iloc[0])
        pr = np.array(F3[(F3.animal == a) & (F3.model == "RPE")].params.iloc[0])
        t = hm.sim_mrpe(pm, c, r, 0.15, 0.85)
        L.append(pd.Series(c.astype(float)).rolling(5, center=True, min_periods=2).mean().values)
        QR.append(rw_q(pr, c, r))
        QM.append(t[:, 0])
        TM.append(t)
    return map(np.array, (L, QR, QM, TM))


def nz(v):
    pre = v[:, 20:60].mean(1, keepdims=True)
    post = v[:, 120:].mean(1, keepdims=True)
    return (v - post) / (pre - post)


def logit_fit(Xm, y):
    w = np.zeros(Xm.shape[1])
    for _ in range(60):
        p = 1 / (1 + np.exp(-(Xm @ w)))
        H = (Xm.T * (p * (1 - p))) @ Xm + 1e-6 * np.eye(len(w))
        step = np.linalg.solve(H, Xm.T @ (p - y) + 1e-6 * w)
        w -= step
        if np.abs(step).max() < 1e-8:
            break
    return w


def group_predictions():
    rows = []
    for g in ("Active", "Passive"):
        for model in ("today's outcome", "RPE", "RPE2a", "mRPE"):
            for m in m5.MICE[g]:
                typ, r, c, tau, early, use, dayout = m5.ARR[m]
                q = m5.q_at_bins(m4.SPECS[model][0], np.asarray(G5[g][model]["params"], float), typ, r, dayout)
                sel = (typ == 1) & (use == 1)
                Xa = np.column_stack([np.ones(len(typ)), tau, early, q])
                w = logit_fit(Xa[sel], c[sel].astype(float))
                A = m5.EV[m][sel].copy()
                A["p"] = 1 / (1 + np.exp(-(Xa[sel] @ w)))
                A["block"] = (A.bin // 6).astype(int)
                rows.append(A.groupby(["day", "block"]).agg(obs=("c", "mean"), pred=("p", "mean")).reset_index()
                            .assign(mouse=m, group=g, model=model))
    return pd.concat(rows)


def drop_variables(m, g):
    """value, RPE, CD, CE, gain M at every drop with the group's meta-RPE parameters (days 11-18)."""
    typ, r, c, tau, early, use, dayout = m5.ARR[m]
    A = m5.EV[m]
    lp = np.asarray(G5[g]["mRPE"]["params"], float)
    Q = md = me = 0.0
    rows = []
    for i in range(len(typ)):
        if typ[i] != 0:
            continue
        delta = r[i] - Q
        md = (1 - lp[1]) * md + lp[1] * delta
        me = (1 - lp[2]) * me + lp[2] * delta
        M = (-md + abs(md)) / 2 + (me + abs(me)) / 2
        Q = Q + min(lp[0] * abs(M), 1.0) * delta
        rows.append(dict(day=int(A.day.iloc[i]), delta=delta, CD=md, CE=me, M=M, Q=Q))
    Tb = pd.DataFrame(rows)
    return Tb[Tb.day >= 11].reset_index(drop=True)


# ============================================================================================ figures
def fig_paper():
    L, QR, QM, TM = paper_curves()
    keep = (L[:, 20:60].mean(1) - L[:, 120:].mean(1)) > 0.1
    check("paper: all 40 mice changed their licking (normalisation defined)", keep.sum() == 40)
    Ln, Rn, Mn = nz(L), nz(QR), nz(QM)
    x = np.arange(160) - 60
    win = (x >= -20) & (x <= 60)
    fig = plt.figure(figsize=(12.7, 4.85))
    a = fig.add_axes([0.062, 0.145, 0.375, 0.735])
    a.axvline(0, color="k", lw=1, ls="--")
    mu = {}
    for key, v, col, lw, lab in (("lick", Ln, "k", 2.2, "licking (running mean, 5 trials)"),
                                 ("rpe", Rn, GREY, 2.6, "RPE model value"),
                                 ("mrpe", Mn, GREEN, 3.0, "meta-RPE model value")):
        m_ = v.mean(0)
        se = v.std(0, ddof=1) / np.sqrt(len(v))
        a.fill_between(x[win], (m_ - se)[win], (m_ + se)[win], color=col, alpha=0.15, lw=0)
        a.plot(x[win], m_[win], color=col, lw=lw, label=lab)
        mu[key] = m_
    for k7, key in (("licking (running mean, 5 trials)", "lick"), ("RPE model value", "rpe"),
                    ("meta-RPE model value", "mrpe")):
        check(f"paper curve {key}: trials 0-2 mean = m07 value",
              abs(float(np.nanmean(nz({"lick": L, "rpe": QR, "mrpe": QM}[key])[:, 60:63])) - N7[k7]["trials_0_2"]) < 1e-3)
    a.axhline(1, color="#BBBBBB", lw=0.8, ls=":")
    a.axhline(0, color="#BBBBBB", lw=0.8, ls=":")
    a.set_xlim(-22, 62)
    a.set_ylim(-0.25, 1.32)
    a.set_yticks([0, 0.5, 1.0])
    a.set_xlabel("Trials from reversal (cue: 85% \u2192 15% reward)")
    a.set_ylabel("Normalized (pre = 1, post = 0)")
    a.set_title("a  Licking and model value, 40 mice", loc="left", fontweight="bold")
    a.legend(loc="lower left", frameon=False)
    b1 = fig.add_axes([0.515, 0.66, 0.235, 0.22])
    d_mu = TM[:, :, 1].mean(0)
    b1.axvline(0, color="k", lw=1, ls="--")
    b1.axhline(0, color="k", lw=0.6)
    b1.fill_between(x[win], 0, d_mu[win], color="#C8C8C8", lw=0)
    b1.plot(x[win], d_mu[win], color="#555555", lw=1.6)
    b1.set_xlim(-22, 62)
    b1.set_ylabel("RPE \u03b4")
    b1.tick_params(labelbottom=False)
    b1.set_title("b  Inside the meta-RPE model, 40 mice", loc="left", fontweight="bold")
    b2 = fig.add_axes([0.515, 0.145, 0.235, 0.47])
    cd, mm = TM[:, :, 2].mean(0), TM[:, :, 4].mean(0)
    b2.axvline(0, color="k", lw=1, ls="--")
    b2.axhline(0, color="k", lw=0.6)
    b2.plot(x[win], cd[win], color=RED, lw=2.4, label="CD = running mean of \u03b4")
    b2.plot(x[win], mm[win], color=GREEN, lw=2.6, ls="--", label="gain M")
    pk = int(x[win][np.argmax(mm[win])])
    check("M peak trial = m07 value", pk == N7["gain_M_peak_trial_from_reversal"])
    b2.annotate(f"M peaks at trial {pk}", xy=(pk, mm[60 + pk]), xytext=(pk + 9, mm[60 + pk] + 0.018),
                fontsize=11, color=GREEN, arrowprops=dict(arrowstyle="-", color=GREEN, lw=1))
    b2.set_xlim(-22, 62)
    b2.set_ylim(-0.215, 0.155)
    b2.set_xlabel("Trials from reversal")
    b2.set_ylabel("CD, gain M")
    b2.legend(loc="lower right", frameon=False, fontsize=10.5)
    c = fig.add_axes([0.872, 0.145, 0.123, 0.735])
    names = hm.ORDER
    lab = {"mRPE": "meta-RPE", "RPE2a": "RPE 2\u03b1", "PH": "PH", "eMack": "eMack", "RPEt-1": "RPE t\u22121", "RPE": "RPE"}
    vals = [S3["BIC"]["mean_change_from_best"][n] for n in names]
    sems = [S3["BIC"]["sem"][n] for n in names]
    yy = np.arange(len(names))
    c.barh(yy, vals, xerr=sems, color=[GREEN] + ["#A8A8A8"] * 5, capsize=3, height=0.65)
    for i, n in enumerate(names[1:], start=1):
        c.text(vals[i] + sems[i] + 0.6, i, "*" if S3["BIC"]["vs_mRPE"][n]["p_t_holm"] < 0.05 else "ns",
               va="center", fontsize=13)
    c.set_yticks(yy, [lab[n] for n in names])
    c.invert_yaxis()
    c.set_xlim(0, 25)
    c.set_xlabel("\u0394BIC from best\n(lower = better)")
    c.set_title("c  Fit, 40 mice", loc="left", fontweight="bold")
    tg = {"key_tl": frac(fig, a, 13, 1.30), "key_br": frac(fig, a, 61, 0.70),
          "tip": frac(fig, a, 3.4, float(Mn.mean(0)[63]) + 0.02)}
    return save(fig, "s1_paper_fig1.png"), fig.get_size_inches(), tg


def fig_map():
    fig = plt.figure(figsize=(5.55, 4.85))
    t = fig.add_axes([0.15, 0.625, 0.82, 0.25])
    days = list(range(3, 19))
    for row, g in ((1, "Active"), (0, "Passive")):
        for d in days:
            morph = (6 <= d <= 13) or d >= 17
            replay = g == "Passive" and 6 <= d <= 10
            t.add_patch(Rectangle((d - 0.43, row - 0.34), 0.86, 0.68, facecolor="#F2C4D6" if morph else "#CFE2F5",
                                  edgecolor="white", hatch="////" if replay else None, lw=0.5))
    t.text(8, -0.02, "replay", ha="center", va="center", fontsize=9.5, color="#333333",
           bbox=dict(facecolor="white", edgecolor="none", pad=0.6, alpha=0.9))
    t.axvspan(10.5, 18.5, ymin=0.0, ymax=1.0, color="#FFF4CC", zorder=0)
    for xs, nm in ((13.5, "switch 1"), (16.5, "switch 2")):
        t.plot([xs, xs], [-0.5, 1.5], color="k", lw=1.6, ls="--")
        t.text(xs, 1.62, nm, ha="center", va="bottom", fontsize=10, fontweight="bold")
    t.text(12.0, -0.72, "fitted", ha="center", va="top", fontsize=9.5, color="#7A5C00")
    t.set_xlim(2.4, 18.6)
    t.set_ylim(-1.0, 1.55)
    t.set_xticks(days, [str(d) for d in days], fontsize=8.5)
    t.set_yticks([1, 0], ["Active", "Passive"], fontsize=11)
    for lbl, col in zip(t.get_yticklabels(), (RED, BLUE)):
        lbl.set_color(col)
        lbl.set_fontweight("bold")
    t.spines["left"].set_visible(False)
    t.tick_params(axis="y", length=0)
    t.set_xlabel("Day", fontsize=10.5, labelpad=1)
    fig.text(0.15, 0.965, "Outcome of every drop (pink = morphine, blue = water)", fontsize=11.5, fontweight="bold",
             va="top")
    # one real session: 0910 black, day 14
    s = fig.add_axes([0.15, 0.095, 0.82, 0.19])
    A = m5.EV["0910_black"]
    A14 = A[(A.day == 14)]
    Lk = pd.read_csv(D / "licks_all.csv")
    L14 = Lk[(Lk.mouse_key == "0910_black") & (Lk.day == 14)]
    lick_t = L14[L14.kind == "lick"].t.values / 60
    drop_t = L14[L14.kind == "drop"].t.values / 60
    bins = A14[A14.typ == 1]
    for _, rr in bins.iterrows():
        if rr.c == 1:
            s.add_patch(Rectangle((rr.t / 60, 0), 0.5, 1, facecolor="#D0D0D0", edgecolor="white", lw=0.6))
    s.vlines(lick_t, 0.50, 0.95, color="k", lw=0.5)
    s.scatter(drop_t, np.full(len(drop_t), 0.25), marker="v", s=16, color="#1565C0", zorder=3)
    s.set_xlim(0, 15)
    s.set_ylim(0, 1)
    s.set_yticks([])
    s.spines["left"].set_visible(False)
    s.set_xticks([0, 3, 6, 9, 12, 15])
    s.set_xlabel("Minutes into the session", fontsize=10.5, labelpad=1)
    s.set_title("One session (Active 0910 black, day 14): grey = 30-s bin with a lick,\n| = lick, \u25bc = water drop",
                fontsize=10, loc="left")
    check("map: 0910 black day 14 has 30 bins", len(bins) == 30)
    tg = {"key_tl": (0.15, 1 - 0.535), "key_br": (0.97, 1 - 0.385), "switch": frac(fig, t, 13.5, -0.5)}
    return save(fig, "s2_mapping.png"), fig.get_size_inches(), tg


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
    nn = {}
    blk = {}
    for g, col in (("Active", RED), ("Passive", BLUE)):
        W = Y[Y.group == g].pivot_table(index="mouse", columns=["day", "block"], values="p")
        m_, se = W.mean(), W.std(ddof=1) / np.sqrt(W.notna().sum())
        nn[g] = W.shape[0]
        blk[g] = m_
        xx = np.arange(len(m_))
        for i in range(len(DAYS)):
            sl = slice(i * 5, i * 5 + 5)
            a.fill_between(xx[sl], (m_ - se).values[sl], (m_ + se).values[sl], color=col, alpha=0.18, lw=0)
            a.plot(xx[sl], m_.values[sl], color=col, lw=2.2, marker="o", ms=3.5,
                   label=f"{g} (n = {W.shape[0]})" if i == 0 else None)
    check("behaviour: n = 6 Active, 7 Passive", nn == {"Active": 6, "Passive": 7})
    a.set_xticks([i * 5 + 2 for i in range(len(DAYS))], [f"day {d}" for d in DAYS])
    a.set_xlim(-0.7, len(DAYS) * 5 - 0.3)
    a.set_ylim(0, 1.62)
    a.set_yticks([0, 0.5, 1.0])
    a.set_ylabel("p(lick in a 30-s bin)")
    a.set_title("a  Licking in 3-min blocks of each 15-min session (mean \u00b1 SEM)", loc="left", fontweight="bold")
    a.legend(loc="upper left", bbox_to_anchor=(0.0, 0.905), frameon=False)
    y14 = float(blk["Active"][(14, 1)])
    for k, (sw, d0, d1, ttl) in enumerate((("day14_morphine_to_water", 13, 14, "b  Day 14 vs 13\n(morphine \u2192 water)"),
                                           ("day17_water_to_morphine", 16, 17, "c  Day 17 vs 16\n(water \u2192 morphine)"))):
        b = fig.add_axes([0.645 + k * 0.2, 0.145, 0.15, 0.735])
        b.axhline(0, color="k", lw=0.8)
        b.axvspan(-0.45, 2.45, color="#F1F1F1", zorder=0)
        for g, col, off in (("Active", RED, -0.08), ("Passive", BLUE, 0.08)):
            r_ = MF[sw][g]
            b.errorbar(np.arange(5) + off, r_["mean_change_by_block"], yerr=r_["sem_by_block"], color=col, lw=2,
                       marker="o", ms=5, capsize=3, label=f"{g} (n = {r_['n']})")
        ra, rp = MF[sw]["Active"], MF[sw]["Passive"]
        txt = (f"first 9 min (grey) vs 0:\nActive {sf(ra['change_first9min'])}, {peq(ra['p_t_first9min_vs0'])}\n"
               f"Passive {sf(rp['change_first9min'])}, {peq(rp['p_t_first9min_vs0'])}")
        b.text(0.03, 0.03 if k == 0 else 0.97, txt, transform=b.transAxes, fontsize=10, va="bottom" if k == 0 else "top")
        b.set_xticks(range(5), ["3", "6", "9", "12", "15"])
        b.set_xlim(-0.5, 4.5)
        b.set_ylim(-0.78, 0.98)
        b.set_xlabel("Minute (block end)")
        b.set_ylabel(f"Change in p(lick) vs day {d0}" if k == 0 else "")
        b.set_title(ttl, loc="left", fontweight="bold", fontsize=12.5)
        if k == 0:
            b.legend(loc="upper right", frameon=False, fontsize=10)
    tg = {"key_tl": frac(fig, a, 9.6, 1.47), "key_br": frac(fig, a, 34.4, 1.04), "tip": frac(fig, a, 16, y14 + 0.07)}
    return save(fig, "s3_behaviour_switches.png"), fig.get_size_inches(), tg


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
        curves = {}
        for model in ("RPE", "mRPE", "RPE2a", "today's outcome"):
            curves[model] = Yg[Yg.model == model].pivot_table(index="mouse", columns=["day", "block"], values="pred").mean()
        for i in range(len(DD)):
            sl = slice(i * 5, i * 5 + 5)
            a.errorbar(xx[sl], m_.values[sl], yerr=se.values[sl], color="k", lw=0, elinewidth=1.2, marker="o", ms=4.5,
                       capsize=2, zorder=4, label="observed (mean \u00b1 SEM)" if i == 0 else None)
            for model, lab, lw, ls in (("today's outcome", "today's outcome (instant switch)", 1.9, (0, (4, 2))),
                                       ("RPE", "RPE (one learning rate)", 1.9, "-"),
                                       ("mRPE", "meta-RPE (paper)", 2.6, "-"),
                                       ("RPE2a", "RPE 2\u03b1 (rates for \u2212 / + errors)", 2.6, "-")):
                a.plot(xx[sl], curves[model].values[sl], color=CMOD[model], lw=lw, ls=ls, zorder=3,
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
            a.legend(loc="upper center", bbox_to_anchor=(0.5, 1.005), ncol=2, frameon=False, fontsize=10,
                     columnspacing=1.0, handlelength=2.2)
        else:
            y = float(curves["RPE2a"][(14, 2)])
            tg = {"key_tl": frac(fig, a, -0.4, 1.60), "key_br": frac(fig, a, 29.4, 1.05), "tip": frac(fig, a, 7, y + 0.06)}
    c = fig.add_axes([0.862, 0.20, 0.133, 0.68])
    models = ["RPE2a", "mRPE", "eMack", "RPE", "RPEt-1", "today's outcome"]
    nm = {"RPE2a": "RPE 2\u03b1", "mRPE": "meta-RPE", "eMack": "eMack", "RPE": "RPE", "RPEt-1": "RPE t\u22121",
          "today's outcome": "today's outcome"}
    yy = np.arange(len(models))
    for g, col, off in (("Active", RED, -0.2), ("Passive", BLUE, 0.2)):
        vals = [G5[g][mm]["dBIC_from_best"] for mm in models]
        c.barh(yy + off, [min(v, 35) for v in vals], height=0.38, color=col, label=g)
        for yv, v in zip(yy + off, vals):
            c.text(min(v, 35) + 0.8, yv, f"{v:.0f}", va="center", fontsize=9.5, color=col)
    c.axvline(10, color="k", lw=0.9, ls=":")
    c.set_yticks(yy, [nm[mm] for mm in models], fontsize=11)
    c.invert_yaxis()
    c.set_xlim(0, 40)
    c.set_xlabel("\u0394BIC from best\n(lower = better)", fontsize=12)
    c.set_title("c  Group fits", loc="left", fontweight="bold")
    c.legend(loc="upper right", frameon=False, fontsize=10)
    return save(fig, "s4_group_model_fits.png"), fig.get_size_inches(), tg


def fig_recovery():
    fig = plt.figure(figsize=(12.7, 4.85))
    a = fig.add_axes([0.085, 0.20, 0.35, 0.68])
    conf = REC["confusion_BIC"]
    gens = ["mRPE", "RPE", "today's outcome"]
    for j, (g, col) in enumerate((("Active", RED), ("Passive", BLUE))):
        vals = [conf[g][gen].get(gen, 0.0) * 100 for gen in gens]
        xs = np.arange(3) + (j - 0.5) * 0.38
        a.bar(xs, vals, width=0.38, color=col, label=g)
        for xv, v in zip(xs, vals):
            a.text(xv, v + 1.5, f"{v:.0f}%", ha="center", va="bottom", fontsize=11)
    a.set_xticks(range(3), ["meta-RPE", "RPE", "today's\noutcome"])
    a.set_ylim(0, 132)
    a.set_yticks([0, 25, 50, 75, 100])
    a.set_ylabel("% of simulated sets where\nper-mouse BIC picks the true model")
    a.set_xlabel("Model that generated the simulated licking")
    a.set_title("a  Can a per-mouse test find the true model?", loc="left", fontweight="bold")
    a.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), ncol=1, frameon=False)
    b = fig.add_axes([0.655, 0.20, 0.335, 0.68])
    ylab = []
    for j, (g, col) in enumerate((("Active", RED), ("Passive", BLUE))):
        r_ = BS[g]
        for k2, (key, lab, mk) in enumerate((("a_minus", "\u03b1\u2212 (morphine \u2192 water)", "v"),
                                             ("a_plus", "\u03b1+ (water \u2192 morphine)", "^"))):
            lo, med, hi = r_[f"{key}_CI"]
            full = r_["full_fit"][key]
            y = j * 2.7 + k2
            b.errorbar(full, y, xerr=[[max(full - lo, 0)], [max(hi - full, 0)]], color=col, marker=mk, ms=9, capsize=4,
                       lw=2.2)
            b.text(np.sqrt(lo * hi), y - 0.3, f"{full:.2f} [{lo:.2f}, {hi:.2f}]", ha="center", fontsize=10.5, color=col)
            ylab.append((y, f"{g} {lab}"))
        b.text(0.0092, j * 2.7 + 1.62, f"{g}: \u03b1\u2212 < \u03b1+ in {r_['frac_a_minus_less_than_a_plus'] * 100:.1f}% "
               f"of bootstraps", fontsize=10.5, color=col)
    b.set_xscale("log")
    b.set_xlim(0.008, 1.6)
    b.set_ylim(4.55, -0.85)
    b.set_yticks([y for y, _ in ylab], [t for _, t in ylab], fontsize=10.5)
    b.set_xlabel("Learning rate (log scale)")
    b.set_title("b  Learning rates, RPE 2\u03b1 group fit (95% CI)", loc="left", fontweight="bold")
    tg = {"key_tl": frac(fig, a, -0.5, 128), "key_br": frac(fig, a, 1.45, 58),
          "tip": frac(fig, a, -0.19, conf["Active"]["mRPE"].get("mRPE", 0) * 100 + 7)}
    return save(fig, "s5_recovery_and_asymmetry.png"), fig.get_size_inches(), tg


def fig_vars():
    fig = plt.figure(figsize=(9.25, 4.85))
    tg, info = {}, {}
    for k, (m, g) in enumerate((("0910_black", "Active"), ("0910_orange", "Passive"))):
        Tb = drop_variables(m, g)
        xx = np.arange(len(Tb))
        first14 = int(Tb.index[Tb.day == 14][0])
        last16 = int(Tb.index[Tb.day == 16][-1])
        imin = int(Tb.CD.values[first14:last16 + 1].argmin()) + first14
        info[g] = dict(mouse=m, n_water=last16 - first14 + 1, cd_min_after=imin - first14 + 1, cd_min=float(Tb.CD[imin]))
        for col_i, kind in enumerate(("val", "cd")):
            ax = fig.add_axes([0.085 + col_i * 0.50, 0.565 - k * 0.455, 0.395, 0.285])
            for d in sorted(Tb.day.unique()):
                idx = Tb.index[Tb.day == d]
                if 14 <= d <= 16:
                    ax.axvspan(idx[0] - 0.5, idx[-1] + 0.5, color=WBG, zorder=0, lw=0)
                ax.axvline(idx[0] - 0.5, color="#CCCCCC", lw=0.6, zorder=0)
                ax.text((idx[0] + idx[-1]) / 2, 1.04, f"d{d}", transform=ax.get_xaxis_transform(), ha="center",
                        va="bottom", fontsize=8.5, color="#444444")
            if kind == "val":
                ax.bar(xx, Tb.delta, color="#BDBDBD", width=1.0, label="RPE \u03b4 per drop")
                ax.plot(xx, Tb.Q, color="k", lw=2, label="value Q (expected drug)")
                ax.set_ylim(-1.1, 1.15)
                ax.set_ylabel("Q, \u03b4", fontsize=11.5)
                ax.legend(loc="lower left", frameon=False, fontsize=9)
                ttl = f"{'ac'[k]}  {g} {m.replace('_', ' ')}: value and RPE"
            else:
                ax.axhline(0, color="k", lw=0.6)
                ax.plot(xx, Tb.CD, color=RED, lw=2.2, label="CD (\u2212 errors, slow)")
                ax.plot(xx, Tb.CE, color=BLUE, lw=1.5, label="CE (fast)")
                ax.plot(xx, Tb.M, color=GREEN, lw=2.2, ls="--", label="gain M")
                ax.scatter([imin], [Tb.CD[imin]], s=40, color=RED, zorder=5, edgecolors="k", linewidths=0.6)
                ax.set_ylabel("CD, CE, M", fontsize=11.5)
                ax.legend(loc="lower left" if g == "Active" else "lower right", frameon=False, fontsize=9)
                ttl = f"{'bd'[k]}  {g}: averaged errors and gain"
                lo_, hi_ = ax.get_ylim()
                ax.set_ylim(lo_ - 0.12 * (hi_ - lo_), hi_)
                tg[g] = frac(fig, ax, imin, float(Tb.CD[imin]))
            ax.set_title(ttl, loc="left", fontsize=11.5, fontweight="bold", pad=21)
            ax.tick_params(labelsize=10)
            ax.set_xlabel("Drops in order, days 11\u201318" if k == 1 else "", fontsize=10.5)
    return save(fig, "s6_model_variables.png"), fig.get_size_inches(), tg, info


# ============================================================================================ pptx helpers
FONTS = {"EN": {False: r"C:\Windows\Fonts\arial.ttf", True: r"C:\Windows\Fonts\arialbd.ttf"},
         "KR": {False: r"C:\Windows\Fonts\malgun.ttf", True: r"C:\Windows\Fonts\malgunbd.ttf"}}
FACE = {"EN": "Arial", "KR": "Malgun Gothic"}
LH = {"EN": 1.22, "KR": 1.40}
_FC = {}


def _f(size, bold):
    k = (LANG, size, bold)
    if k not in _FC:
        _FC[k] = ImageFont.truetype(FONTS[LANG][bold], int(round(size * 300 / 72)))
    return _FC[k]


def lines_needed(runs, width_in):
    W, cur, lines, maxsz = width_in * 300, 0.0, 1, 0
    for text, st in runs:
        f = _f(st.get("size", 18), st.get("bold", False))
        maxsz = max(maxsz, st.get("size", 18))
        for tok in text.split(" "):
            w = f.getlength(tok + " ")
            if cur > 0 and cur + f.getlength(tok) > W:
                lines += 1
                cur = w
            else:
                cur += w
    return lines, maxsz


def rgb(h):
    h = h.lstrip("#")
    return RGBColor(int(h[0:2], 16), int(h[2:4], 16), int(h[4:6], 16))


def measure(paras, w):
    need = 0.10
    for runs, align, after in paras:
        nl, msz = lines_needed(runs, w - 0.30)
        need += nl * LH[LANG] * msz / 72 + after / 72
    return need


def tbox(s, x, y, w, h, paras, name, anchor=MSO_ANCHOR.TOP, fill=None, line=None, region_bottom=None):
    if h is None:
        h = measure(paras, w) + 0.06
    if fill or line:
        tb = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
        tb.adjustments[0] = 0.08
        if fill:
            tb.fill.solid()
            tb.fill.fore_color.rgb = rgb(fill)
        else:
            tb.fill.background()
        if line:
            tb.line.color.rgb = rgb(line)
            tb.line.width = Pt(1.25)
        else:
            tb.line.fill.background()
        tb.shadow.inherit = False
    else:
        tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    name = f"S{SLIDE} {name}"
    tb.name = name
    tf = tb.text_frame
    tf.word_wrap = True
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.vertical_anchor = anchor
    tf.margin_left = tf.margin_right = Inches(0.1)
    tf.margin_top = tf.margin_bottom = Inches(0.05)
    for i, (runs, align, after) in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        p.alignment = align
        p.space_after = Pt(after)
        p.line_spacing = 1.0
        for text, st in runs:
            r = p.add_run()
            r.text = text
            fo = r.font
            fo.size = Pt(st.get("size", 18))
            fo.bold = st.get("bold", False)
            fo.color.rgb = rgb(st.get("color", "222222"))
            fo.name = FACE[LANG]
            rPr = r._r.get_or_add_rPr()
            if LANG == "KR":
                rPr.set("lang", "ko-KR")
                ea = rPr.find(qn("a:ea"))
                if ea is None:
                    ea = etree.SubElement(rPr, qn("a:ea"))
                ea.set("typeface", FACE[LANG])
    need = measure(paras, w)
    ok = need <= h + 0.01 and x >= 0 and y >= 0 and x + w <= 13.334 and y + h <= 7.5
    if region_bottom is not None:
        ok = ok and y + h <= region_bottom + 0.02
    FIT.append(dict(lang=LANG, name=name, need_in=round(need, 3), box_in=round(h, 3), ok=bool(ok)))
    check(f"[{LANG}] text fits: {name} ({need:.2f} in needed, box {h:.2f} in"
          + (f", must end above {region_bottom:.2f} in, ends {y + h:.2f}" if region_bottom else "") + ")", ok)
    return tb


def P(*runs, align=PP_ALIGN.LEFT, after=0):
    return (list(runs), align, after)


def Rn(text, size=18, bold=False, color="222222"):
    return (text, dict(size=size, bold=bold, color=color))


def arrow(s, x0, y0, x1, y1, color="333333", width=2.25):
    c = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x0), Inches(y0), Inches(x1), Inches(y1))
    c.line.color.rgb = rgb(color)
    c.line.width = Pt(width)
    ln = c.line._get_or_add_ln()
    te = etree.SubElement(ln, qn("a:tailEnd"))
    te.set("type", "triangle")
    te.set("w", "med")
    te.set("len", "med")
    c.name = "Key arrow"
    return c


def new_slide(prs):
    lay = next(lay for lay in prs.slide_layouts if lay.name == "1_Title Slide")
    s = prs.slides.add_slide(lay)
    bg = s.background.fill
    bg.solid()
    bg.fore_color.rgb = rgb("FFFFFF")
    return s


def place(s, path, size_in, x, y, w):
    pic = s.shapes.add_picture(str(path), Inches(x), Inches(y), Inches(w))
    pic.name = path.stem
    return x, y, w, w * size_in[1] / size_in[0]


def at(box, f):
    x, y, w, h = box
    return x + w * f[0], y + h * f[1]


def key_box(s, box, tg, paras, tip=None, name="KEY", anchor="bottom"):
    """KEY box in an empty region of the figure; the arrow leaves the box bottom (straight down to the tip's x) or its
    left edge, so it never crosses the box."""
    x0, y0 = at(box, tg["key_tl"])
    x1, y1 = at(box, tg["key_br"])
    tb = tbox(s, x0, y0, x1 - x0, None, paras, name, fill="F2F6FB", line="1F4E79", region_bottom=y1)
    if tip is not None:
        h = measure(paras, x1 - x0) + 0.06
        tx, ty = at(box, tip)
        if anchor == "bottom":
            check(f"[{LANG}] {name}: arrow tip lies under the box", x0 + 0.05 <= tx <= x1 - 0.05 and ty > y0 + h)
            arrow(s, tx, y0 + h, tx, ty)
        else:
            check(f"[{LANG}] {name}: arrow tip lies left of the box", tx < x0)
            arrow(s, x0, y0 + 0.75 * h, tx, ty)
    return tb


def common(s, title, did, found, foot, notes):
    ts = 28 if LANG == "EN" else 26
    tbox(s, 0.35, 0.10, 12.63, 0.80, [P(Rn(title, ts, True, NAVY), align=PP_ALIGN.CENTER)], "Title",
         anchor=MSO_ANCHOR.MIDDLE)
    bs = 16 if LANG == "EN" else 14.5
    tbox(s, 0.45, 5.90, 12.45, 0.92, [P(Rn(T("Did: ", "한 것: "), bs, True), Rn(did, bs),
                                         Rn(T("  Found: ", "  찾은 것: "), bs, True), Rn(found, bs))], "Bottom line")
    fs = 11 if LANG == "EN" else 10
    tbox(s, 0.45, 6.83, 12.45, 0.665, [P(Rn(foot, fs, color="555555"))], "Footer (n and test)")
    s.notes_slide.notes_text_frame.text = notes


# ============================================================================================ slides
REPO = "github.com/stuberlab/Hjort-et-al.-2026-PFC-and-reversal-learning"


def slide1(prs, F):
    global SLIDE
    SLIDE = 1
    path, size, tg = F
    s = new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    B = S3["BIC"]
    others = [B["mean_change_from_best"][n] for n in hm.ORDER[1:]]
    hp = max(B["vs_mRPE"][n]["p_t_holm"] for n in hm.ORDER[1:])
    mb = reg("s1_bic_mrpe", f"{B['mean_change_from_best']['mRPE']:.1f}", B["mean_change_from_best"]["mRPE"])
    lo_ = reg("s1_bic_lo", f"{min(others):.1f}", min(others))
    hi_ = reg("s1_bic_hi", f"{max(others):.1f}", max(others))
    hps = reg("s1_holm_max", f"{hp:.3f}", hp)
    nb = S3["negLL"]["n_best"]["mRPE"]
    v_m = N7["meta-RPE model value"]["trials_0_2"]
    v_r = N7["RPE model value"]["trials_0_2"]
    v_l = N7["licking (running mean, 5 trials)"]["trials_0_2"]
    reg("s1_vm", f"{v_m:.2f}", v_m)
    reg("s1_vr", f"{v_r:.2f}", v_r)
    reg("s1_vl", f"{v_l:.2f}", v_l)
    pk = N7["gain_M_peak_trial_from_reversal"]
    ks = 13.5 if LANG == "EN" else 12
    key_box(s, box, tg, [P(Rn(T("Value waits for errors to pile up", "오류가 쌓일 때까지 가치는 기다린다"), ks + 0.5, True,
                                "1F4E79"), after=3),
                         P(Rn(T(f"Trials 0\u20132 after reversal: meta-RPE {v_m:.2f}, RPE {v_r:.2f}, licking {v_l:.2f}.",
                                f"reversal 뒤 0\u20132 trial: meta-RPE {v_m:.2f}, RPE {v_r:.2f}, 핥기 {v_l:.2f}"),
                              ks))], tip=tg["tip"], anchor="left")
    common(s, T("Hjort et al. 2026: behaviour changes only after errors pile up",
                "Hjort et al. 2026: 오류가 쌓인 뒤에야 행동이 바뀐다"),
           T("re-fitted the authors' 40 mice with their own code (ported to Python).",
             "저자 코드(Python으로 옮김)로 저자 데이터 40마리를 다시 맞춤."),
           T(f"the paper's result reproduces: meta-RPE fits best (\u0394BIC {mb} vs {lo_}\u2013{hi_}; better than each "
             f"model, largest Holm p = {hps}).",
             f"논문 결과 재현: meta-RPE가 가장 잘 맞음 (\u0394BIC {mb} vs {lo_}\u2013{hi_}; 모든 모델보다 좋음, 가장 큰 "
             f"Holm p = {hps})."),
           T(f"Hjort et al., Nature 655:174 (2026), Fig. 1; code + 40-mouse data: {REPO} (H\u2192L cue, 160 trials, "
             f"reversal at trial 60). a: per mouse, trials 20\u201359 = 1 and 120\u2013159 = 0; mean \u00b1 SEM. b: learning step = "
             f"\u03b1\u00b7\u03b4\u00b7M, M = max(0, \u2212CD) + max(0, CE); CD, CE = running means of \u03b4 at rates d, e. c: BIC = "
             f"k ln T + 2 negLL, mean \u00b1 SEM; * paired t vs meta-RPE, Holm over 5.",
             f"Hjort et al., Nature 655:174 (2026), Fig. 1; 코드, 데이터: {REPO} (H\u2192L 큐, 160 trial). "
             f"a: 마우스마다 trial 20\u201359 = 1, 120\u2013159 = 0; 평균 \u00b1 SEM. b: 학습 크기 = "
             f"\u03b1\u00b7\u03b4\u00b7M, M = max(0, \u2212CD) + max(0, CE); CD, CE = \u03b4의 이동평균(속도 d, e). c: BIC = "
             f"k ln T + 2 negLL, 평균 \u00b1 SEM; * meta-RPE와 paired t, 5개 Holm 보정."),
           T(notes1_en(mb, lo_, hi_, hps, nb, v_m, v_r, v_l, pk), notes1_kr(mb, lo_, hi_, hps, nb, v_m, v_r, v_l, pk)))
    return s


def notes1_en(mb, lo_, hi_, hps, nb, v_m, v_r, v_l, pk):
    return (
        "WHAT: Fig. 1 of Hjort et al. 2026 (Nature 655:174) redrawn from the authors' own 40 mice with their own model code "
        "(BehaviorModels.ipynb, ported line by line to Python + numba; same equations, bounds, start grid and BIC).\n"
        "TASK: head-fixed mice, 4 odour cues (2 at 85% sucrose, 2 at 15%); after stable licking one unsignalled reversal "
        "(H->L = contingency degradation, L->H = enhancement). Behaviour = anticipatory licks 0-2.5 s after cue onset; the "
        "model is fitted to the H->L cue (160 trials over 4 days, reversal at trial 60).\n"
        f"a: licking stays high for ~6 trials after the reversal. In trials 0-2 the RPE (Rescorla-Wagner) value has already "
        f"dropped to {v_r:.2f} while the meta-RPE value is {v_m:.2f} (licking {v_l:.2f}).\n"
        f"b: the single-trial RPE is largest right after the reversal, but the averaged error CD needs several trials to "
        f"build, so the gain M peaks at trial {pk}; learning step = alpha x delta x M.\n"
        f"c: our port gives the paper's Fig. 1k,l: BIC change from best {mb} (meta-RPE) vs {lo_}-{hi_}; paired t vs "
        f"meta-RPE, Holm-corrected, largest p = {hps}; by negLL meta-RPE is best in {nb}/40 mice.\n"
        "Rest of the paper: 13% of mPFC neurons encode CD (GLM + circular-shift null); holographic activation of CD cells "
        "reduces licking; mPFC dopamine encodes RPE/CE; 61% of mPFC->VTA cells encode CD; D1/D2 block in mPFC lowers CD and "
        "delays reversal; mPFC->VTA stimulation speeds degradation. The authors speculate that volitional drug use impairs CD.")


def notes1_kr(mb, lo_, hi_, hps, nb, v_m, v_r, v_l, pk):
    return (
        "[이 슬라이드] Hjort et al. 2026 (Nature 655:174) Figure 1의 핵심을 저자 데이터 40마리와 저자 코드로 다시 그린 것. "
        "저자 노트북(BehaviorModels.ipynb)을 Python으로 한 줄씩 옮겼고 식, 파라미터 범위, 시작값, BIC 계산은 모두 같다.\n"
        "[과제] 머리 고정 쥐에게 냄새 큐 4개: 2개는 85% 확률로 설탕물(H), 2개는 15%(L). 핥기가 안정되면 예고 없이 한 번 "
        "바꾼다 (H→L = contingency degradation, L→H = enhancement). 행동 지표는 큐 뒤 0–2.5초의 미리 핥기(anticipatory lick). "
        "모델은 H→L 큐의 160 trial(4일, trial 60에서 reversal)에 맞춘다.\n"
        f"[a 읽는 법] 검은 선 = 실제 핥기, 회색 = 일반 RPE(Rescorla–Wagner) 모델의 가치, 초록 = meta-RPE 모델의 가치 "
        f"(마우스마다 reversal 전 = 1, 후 = 0으로 정규화). reversal 직후 0–2 trial에서 핥기 {v_l:.2f}, RPE 가치 {v_r:.2f}, "
        f"meta-RPE 가치 {v_m:.2f}. 즉 일반 RPE 모델은 첫 오류부터 가치를 깎지만, 실제 쥐와 meta-RPE는 몇 trial 버틴다.\n"
        f"[b 읽는 법] 위: 한 trial의 오류 δ = 보상 − 가치. reversal 바로 다음 trial에서 가장 크다. 아래: CD = δ의 이동평균 "
        f"(빨강), 초록 점선 = gain M = max(0, −CD) + max(0, CE). 한 번에 배우는 양 = α × δ × M. 오류가 한 번 나면 M이 작아 "
        f"거의 안 배우고, 같은 방향 오류가 쌓여야 M이 커진다 (trial {pk}에서 최대). CD(음의 오류 누적)와 CE(양의 오류 누적)의 "
        f"속도 d, e를 따로 두어 degradation과 enhancement 속도가 다를 수 있다.\n"
        f"[c 읽는 법] 6개 모델의 BIC 차이 (낮을수록 좋음). meta-RPE {mb} vs 나머지 {lo_}–{hi_}; 모든 모델보다 좋음 "
        f"(paired t, Holm 보정, 가장 큰 p = {hps}). negLL로는 40마리 중 {nb}마리에서 meta-RPE가 최선. 우리 Python 코드가 "
        f"논문 Fig. 1k,l을 재현한다는 확인이다.\n"
        "[논문의 나머지] Fig. 2: mPFC 뉴런의 13%가 CD를 표현 (GLM + circular shift 검정). Fig. 3: CD 세포를 홀로그래픽 "
        "광유전학으로 켜면 핥기 감소. Fig. 4: mPFC 도파민은 RPE/CE를 표현, mPFC→VTA 뉴런의 61%가 CD를 표현, mPFC에 D1/D2 "
        "차단제를 넣으면 CD가 줄고 reversal이 늦어짐. Fig. 5: mPFC→VTA 자극 → degradation이 빨라짐. 저자들은 자발적 약물 "
        "사용 뒤 CD 세포가 손상될 수 있다고 추측한다.")


ROWS = [
    (("Trial", "Trial"), ("cue every ~20 s", "~20초마다 큐"), ("30-s bin of the 15-min session", "15분 세션의 30초 bin")),
    (("Choice", "선택"), ("anticipatory lick: yes / no", "미리 핥기: 있음 / 없음"), ("any lick in the bin: yes / no",
                                                                                "bin 안의 lick: 있음 / 없음")),
    (("Outcome", "결과"), ("sucrose 1 / none 0 (85% \u2194 15%)", "설탕물 1 / 없음 0 (85% \u2194 15%)"),
     ("each drop: morphine 1 / water 0", "매 drop: 모르핀 1 / 물 0")),
    (("Value Q", "가치 Q"), ("expected reward of the cue", "큐의 보상 기대"),
     ("expected drug from the spout; updated at every drop, from 0 on day 3",
      "spout의 약물 기대; drop마다 업데이트, day 3에 0에서 시작")),
    (("Switch", "전환"), ("one reversal: H\u2192L, L\u2192H", "reversal 1번: H\u2192L, L\u2192H"),
     ("day 14 morphine \u2192 water; day 17 water \u2192 morphine", "day 14 모르핀 \u2192 물; day 17 물 \u2192 모르핀")),
    (("Fitted", "맞춘 데이터"), ("160 trials (4 days) per mouse", "마우스당 160 trial (4일)"),
     ("days 11\u201318: 240 bins per mouse", "day 11\u201318: 마우스당 240 bin")),
    (("Readout", "행동 출력"), ("softmax(\u03b2\u00b7Q)", "softmax(\u03b2\u00b7Q)"),
     ("logistic(b0 + \u03b2\u00b7Q + time in session); \u03b2 may be negative",
      "logistic(b0 + \u03b2\u00b7Q + 세션 내 시간); \u03b2는 음수 가능")),
]


def slide2(prs, F):
    global SLIDE
    SLIDE = 2
    path, size, tg = F
    s = new_slide(prs)
    box = place(s, path, size, 7.55, 0.98, 5.55)
    xs, ws = [0.35, 1.70, 4.20], [1.35, 2.50, 2.95]
    fs = 13 if LANG == "EN" else 12
    head = (T("", ""), T("Hjort 2026 (Pavlovian)", "Hjort 2026 (Pavlovian)"), T("Our task (PR licking)", "우리 과제 (PR licking)"))
    y = 1.05
    hh = max(measure([P(Rn(h, fs, True))], w) for h, w in zip(head, ws)) + 0.06
    for x, w, h in zip(xs, ws, head):
        tbox(s, x, y, w, hh, [P(Rn(h, fs, True, "FFFFFF"))], f"table head {h}", fill="1F4E79", anchor=MSO_ANCHOR.MIDDLE)
    y += hh + 0.04
    for i, row in enumerate(ROWS):
        cells = [T(*row[0]), T(*row[1]), T(*row[2])]
        rh = max(measure([P(Rn(cl, fs, j == 0))], w) for j, (cl, w) in enumerate(zip(cells, ws))) + 0.06
        for j, (x, w, cl) in enumerate(zip(xs, ws, cells)):
            tbox(s, x, y, w, rh, [P(Rn(cl, fs, j == 0, "1F4E79" if j == 0 else "222222"))], f"table r{i}c{j}",
                 fill="EEF3F9" if i % 2 == 0 else "FFFFFF", anchor=MSO_ANCHOR.MIDDLE)
        y += rh + 0.04
    check(f"[{LANG}] mapping table ends above the bottom line (ends {y:.2f} in)", y <= 5.88)
    ks = 12.5 if LANG == "EN" else 12
    x0, y0 = at(box, tg["key_tl"])
    x1, y1 = at(box, tg["key_br"])
    paras = [P(Rn(T("Days 11\u201318: same outcomes for both groups, so they are compared at the same two switches.",
                    "day 11\u201318: 두 그룹의 결과가 같으므로 같은 두 전환에서 비교한다."), ks, True, "1F4E79"))]
    tbox(s, x0, y0, x1 - x0, None, paras, "KEY", fill="F2F6FB", line="1F4E79", region_bottom=y1)
    tx, ty = at(box, tg["switch"])
    arrow(s, tx, y0, tx, ty + 0.05)
    nb = int(G5["Active"]["RPE"]["n"] / len(m5.MICE["Active"]))
    reg("s2_bins", str(nb), nb)
    common(s, T("Same equations, different task: the model on our PR licking", "식은 같고 과제는 다르다: 우리 PR licking에 맞춘 방법"),
           T("kept the paper's learning rules, bounds and fitting; redefined trial, choice and outcome.",
             "논문의 학습 규칙, 파라미터 범위, 맞추는 방법은 그대로 두고 trial, 선택, 결과만 다시 정의함."),
           T(f"days 11\u201318 give each mouse one morphine \u2192 water and one water \u2192 morphine switch ({nb} bins per mouse).",
             f"day 11\u201318에서 마우스마다 모르핀 \u2192 물, 물 \u2192 모르핀 전환이 한 번씩 있음 (마우스당 {nb} bin)."),
           T("PR = progressive ratio (+1 lick per drop), sessions end at ~15 min. Passive days 6\u201310 = replay of the partner's "
             "drops (morphine drunk, licks earn nothing): their drops update value, their licking is not fitted. Passive 6099 white "
             "died on day 13 (no switch): left out of model comparisons. Raw controller logs, K:\\addiction_concate_Dec_2025 (read only).",
             "PR = progressive ratio (drop마다 lick +1), 세션은 ~15분. Passive day 6\u201310 = 짝의 drop을 재생(replay): 모르핀은 "
             "마셨지만 lick은 보상과 무관 \u2192 drop은 가치 업데이트에 넣고 lick은 맞추지 않음. Passive 6099 white는 day 13에 사망 "
             "(전환 없음) \u2192 모델 비교에서 제외. 원자료: K:\\addiction_concate_Dec_2025 컨트롤러 로그 (읽기만)."),
           T(notes2_en(nb), notes2_kr(nb)))
    return s


def notes2_en(nb):
    return (
        "WHAT: what has to change to use the paper's model on our data.\n"
        "SAME: the learning rules (RPE, meta-RPE, RPE 2alpha, RPE t-1, Pearce-Hall, eMack), their bounds, maximum-likelihood "
        "fitting and BIC.\n"
        "DIFFERENT: the paper changes the reward PROBABILITY of a cue (Pavlovian); in our free-operant PR task every completed "
        "ratio gives a drop and its IDENTITY changes (morphine <-> water). So trial = 30-s bin, choice = any lick in the bin, "
        "outcome = drop identity (morphine 1, water 0).\n"
        "VALUE Q = expected drug from the spout: 0 on day 3 (only water before), updated by every drop the mouse drank (Active "
        "PR days 6-10 and Passive replay days 6-10 included), carried over between days. In an operant task a drop only follows "
        "licking, so the paper's 'update only the chosen option' becomes 'update the lick value at each drop', and its softmax "
        "becomes logistic(b0 + beta*Q). Added for every model alike: time in session and a first-3-min term (licking falls "
        "within every session). beta may be negative (a mouse that licks less when it expects morphine).\n"
        f"SWITCHES: day 14 morphine -> water = the paper's H->L; day 17 water -> morphine = L->H. Likelihood on days 11-18 "
        f"({nb} bins per mouse), where both groups run the same PR schedule with the same outcome sequence.\n"
        "Right, bottom: one real session (0910 black, day 14): grey bins had at least one lick; triangles = water drops.")


def notes2_kr(nb):
    return (
        "[이 슬라이드] 논문 모델을 우리 데이터에 쓰려면 무엇을 바꿔야 하는지.\n"
        "[같은 것] 학습 규칙 6개(RPE, meta-RPE, RPE 2α, RPE t−1, Pearce–Hall, eMack), 파라미터 범위, 최대우도로 맞추기, "
        "BIC 비교는 논문 그대로.\n"
        "[다른 것] 논문은 큐의 '보상 확률'이 바뀌는 Pavlovian 과제. 우리는 free-operant PR 과제라서 정해진 수만큼 핥으면 "
        "drop이 항상 나오고, 그 drop의 '정체'(모르핀 ↔ 물)가 바뀐다. 그래서 trial = 30초 bin, 선택 = bin 안에 lick이 있었는지, "
        "결과 = drop의 정체(모르핀 1, 물 0)로 다시 정의했다.\n"
        "[가치 Q] 'spout에서 기대하는 약물'. day 3(그전엔 물만 경험)에 0에서 시작, 마신 drop마다 업데이트, 날이 바뀌어도 "
        "이어진다. Passive의 day 6–10 replay drop도 마셨으므로 업데이트에 넣었다 (그날 lick은 보상과 무관하니 모델 맞추기에선 "
        "제외). operant 과제에서는 drop이 lick 뒤에만 나오므로, 논문의 '선택한 쪽만 업데이트'는 'drop마다 lick 가치를 "
        "업데이트'가 되고 softmax는 logistic(b0 + β·Q)가 된다. 모든 모델에 똑같이 세션 내 시간 항과 처음 3분 항을 넣었다 "
        "(어느 날이든 세션 안에서 lick이 줄기 때문). β는 음수도 가능 (모르핀을 기대하면 덜 핥는 쥐 = Passive).\n"
        f"[전환] day 14 모르핀 → 물 = 논문의 H→L, day 17 물 → 모르핀 = L→H. 모델은 day 11–18에 맞춤 (마우스당 {nb} bin). "
        "이 기간에는 두 그룹이 같은 PR 일정과 같은 결과 순서를 겪으므로 같은 전환에서 비교할 수 있다.\n"
        "[오른쪽 아래] 실제 한 세션 (Active 0910 black, day 14): 회색 칸 = lick이 있었던 30초 bin, ▼ = 물 drop. 모델은 "
        "bin마다 'lick이 있을 확률'을 예측한다.")


def slide3(prs, F):
    global SLIDE
    SLIDE = 3
    path, size, tg = F
    s = new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    a14, p14 = MF["day14_morphine_to_water"]["Active"], MF["day14_morphine_to_water"]["Passive"]
    a17, p17 = MF["day17_water_to_morphine"]["Active"], MF["day17_water_to_morphine"]["Passive"]
    asy = MF["asymmetry_first9min_abs_change_day17_vs_day14"]["Active"]
    va = reg("s3_a14", sf(a14["change_first9min"]), a14["change_first9min"])
    vp = reg("s3_p14", sf(p14["change_first9min"]), p14["change_first9min"])
    va17 = reg("s3_a17", sf(a17["change_first9min"]), a17["change_first9min"])
    pa = reg("s3_pa14", fp(a14["p_t_first9min_vs0"]), a14["p_t_first9min_vs0"])
    pp = reg("s3_pp14", fp(p14["p_t_first9min_vs0"]), p14["p_t_first9min_vs0"])
    pa17 = reg("s3_pa17", fp(a17["p_t_first9min_vs0"]), a17["p_t_first9min_vs0"])
    ple_a = reg("s3_le_a", fp(a14["p_t_late_vs_early"]), a14["p_t_late_vs_early"])
    ple_p = reg("s3_le_p", fp(p14["p_t_late_vs_early"]), p14["p_t_late_vs_early"])
    as14 = reg("s3_asy14", f"{asy['mean_abs_day14']:.2f}", asy["mean_abs_day14"])
    as17 = reg("s3_asy17", f"{asy['mean_abs_day17']:.2f}", asy["mean_abs_day17"])
    pasy = reg("s3_asy_p", fp(asy["p_t"]), asy["p_t"])
    ks = 14.5 if LANG == "EN" else 13.5
    key_box(s, box, tg, [P(Rn(T(f"First water day, first 9 min: Active unchanged ({va}); Passive up at once ({vp}).",
                                f"첫 물 날(day 14) 처음 9분: Active 변화 없음 ({va}); Passive 바로 증가 ({vp})."), ks, True,
                              "1F4E79"))], tip=tg["tip"])
    common(s, T("After morphine is removed, Active mice keep licking for 9 minutes",
                "모르핀을 빼도 Active 쥐는 9분 동안 그대로 핥는다"),
           T("compared licking on each switch day with the day before, in 3-min blocks, within mice.",
             "전환 날의 lick을 전날과 3분 단위로 비교함 (같은 쥐 안에서)."),
           T(f"morphine \u2192 water: Active {va} ({peq(a14['p_t_first9min_vs0'])}), Passive {vp} "
             f"({peq(p14['p_t_first9min_vs0'])}); water \u2192 morphine: Active {va17} ({peq(a17['p_t_first9min_vs0'])}): "
             f"slow to stop, fast to resume.",
             f"모르핀 \u2192 물: Active {va} ({peq(a14['p_t_first9min_vs0'])}), Passive {vp} ({peq(p14['p_t_first9min_vs0'])}); "
             f"물 \u2192 모르핀: Active {va17} ({peq(a17['p_t_first9min_vs0'])}): 멈추는 건 느리고 다시 시작은 빠름."),
           T(f"p(lick) = share of 30-s bins with \u2265 1 lick. n = 6 Active, 7 Passive (6099 white died on day 13, left out). "
             f"b, c: change vs the day before per mouse, mean \u00b1 SEM; first 9 min vs 0: one-sample t-test, two-sided. Last 6 "
             f"vs first 9 min (paired t): Active p = {ple_a}, Passive p = {ple_p}. Active size of the first-9-min change, day 14 vs "
             f"day 17: {as14} vs {as17}, paired t p = {pasy}.",
             f"p(lick) = lick이 1번 이상 있었던 30초 bin의 비율. n = Active 6, Passive 7 (6099 white는 day 13 사망으로 제외). "
             f"b, c: 마우스마다 전날 대비 변화, 평균 \u00b1 SEM; 처음 9분 vs 0: one-sample t-test, 양측. 마지막 6분 vs 처음 9분 "
             f"(paired t): Active p = {ple_a}, Passive p = {ple_p}. Active 처음 9분 변화 크기, day 14 vs day 17: {as14} vs "
             f"{as17}, paired t p = {pasy}."),
           T(notes3_en(va, vp, va17, pa, pp, pa17, ple_a, ple_p, as14, as17, pasy, p17),
             notes3_kr(va, vp, va17, pa, pp, pa17, ple_a, ple_p, as14, as17, pasy, p17)))
    return s


def notes3_en(va, vp, va17, pa, pp, pa17, ple_a, ple_p, as14, as17, pasy, p17):
    return (
        "WHAT: behaviour without any model. a: p(lick) per 3-min block of each 15-min PR session, days 11-18; light blue = "
        "water days. b: day 14 minus day 13 (morphine -> water); c: day 17 minus day 16 (water -> morphine).\n"
        f"Day 14, first 9 min (grey): Active {va} (p = {pa}): on the first water day they lick as on a morphine day; Passive "
        f"{vp} (p = {pp}): they increase at once. Active licking falls only at the end of day 14 and on days 15-16.\n"
        f"Day 17, first 9 min: Active {va17} (p = {pa17}): back at once; Passive {sf(p17['change_first9min'])} "
        f"(p = {fp(p17['p_t_first9min_vs0'])}).\n"
        "MEANING: Active mice stop slowly when the drug is removed but resume fast when it returns - the same direction as the "
        "paper's 'H->L slower than L->H'. Passive mice switch at once both ways (possibly taste: morphine is bitter).\n"
        f"CAVEATS: Active size of change day 14 vs day 17 {as14} vs {as17}, p = {pasy} (6 mice); late vs early on day 14: "
        f"Active p = {ple_a}, Passive p = {ple_p}. The fast return on day 17 could be the drug acting within the session, and "
        "the persistence on day 14 could be a withdrawal state: behaviour alone cannot separate these from learning.")


def notes3_kr(va, vp, va17, pa, pp, pa17, ple_a, ple_p, as14, as17, pasy, p17):
    return (
        "[이 슬라이드] 모델 없이 본 행동.\n"
        "[a 읽는 법] 15분 세션을 3분 블록 5개로 나눠 p(lick) (lick이 있었던 30초 bin의 비율)을 day 11–18에 걸쳐 그렸다. "
        "빨강 = Active, 파랑 = Passive, 연파랑 배경 = 물 날.\n"
        f"[b 읽는 법] day 14 − day 13 (모르핀 → 물). 처음 9분(회색 영역): Active {va} (p = {pa}) → 물인데도 모르핀 날처럼 "
        f"핥는다. Passive {vp} (p = {pp}) → 바로 증가. Active는 세션 끝무렵과 다음 날(day 15–16)에야 줄어든다.\n"
        f"[c 읽는 법] day 17 − day 16 (물 → 모르핀). 처음 9분: Active {va17} (p = {pa17}) → 바로 회복. Passive "
        f"{sf(p17['change_first9min'])} (p = {fp(p17['p_t_first9min_vs0'])}).\n"
        "[의미] Active: 약을 뺐을 때는 늦게, 돌려줬을 때는 빨리 바뀐다 = 논문의 'H→L이 L→H보다 느리다'와 같은 방향. "
        "Passive: 양쪽 다 바로 바뀐다 (학습이라기보다 맛일 수 있음, 모르핀은 쓰다).\n"
        f"[주의] Active 안에서 두 변화의 크기를 직접 비교하면 day 14 {as14} vs day 17 {as17}, p = {pasy} (6마리로는 유의하지 "
        f"않음). day 14의 마지막 6분 vs 처음 9분: Active p = {ple_a}, Passive p = {ple_p}. day 17의 빠른 회복은 세션 중 마신 "
        "모르핀의 약효일 수 있고, day 14에 그대로인 것은 금단 상태 때문일 수도 있다 — 행동만으로는 학습과 구분할 수 없다.")


def slide4(prs, F):
    global SLIDE
    SLIDE = 4
    path, size, tg = F
    s = new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    A, Pv = G5["Active"], G5["Passive"]
    d = {f"{g}_{m}": reg(f"s4_{g}_{m}", f"{G5[g][m]['dBIC_from_best']:.0f}", G5[g][m]["dBIC_from_best"])
         for g in ("Active", "Passive") for m in ("RPE2a", "mRPE", "RPE", "today's outcome", "eMack", "RPEt-1")}
    phmin = min(G5[g][m]["dBIC_from_best"] for g in ("Active", "Passive") for m in ("PH", "no learning"))
    check(f"PH / no learning dBIC >= 70 (min {phmin:.1f})", phmin >= 70)
    a_to = d["Active_today's outcome"]
    nA, nP = A["RPE"]["n"], Pv["RPE"]["n"]
    check("best models: Active RPE2a, Passive today's outcome",
          A["RPE2a"]["dBIC_from_best"] == 0 and Pv["today's outcome"]["dBIC_from_best"] == 0)
    ks = 13.5 if LANG == "EN" else 12.5
    key_box(s, box, tg, [P(Rn(T("Active best fit: two learning rates, slow after morphine \u2192 water, fast after "
                                "water \u2192 morphine.", "Active 최적 모델: 학습률 2개 (모르핀 \u2192 물 뒤엔 느리게, 물 \u2192 모르핀 뒤엔 빠르게)"),
                              ks, True, "1F4E79"))],
            tip=tg["tip"])
    common(s, T("Active mice learn slowly from the loss of morphine", "Active 쥐는 모르핀이 없어진 것을 천천히 배운다"),
           T("fitted 8 learning rules, learning rates shared within each group; compared by BIC.",
             "학습 규칙 8개를 맞춤 (학습률은 그룹 안에서 공유); BIC로 비교."),
           T(f"Active: two rates best (\u0394BIC meta-RPE +{d['Active_mRPE']}, one rate +{d['Active_RPE']}, instant switch "
             f"+{a_to}); Passive: instant switch best (meta-RPE +{d['Passive_mRPE']}).",
             f"Active: 학습률 2개가 최선 (\u0394BIC meta-RPE +{d['Active_mRPE']}, 학습률 1개 +{d['Active_RPE']}, 즉시 전환 "
             f"+{a_to}); Passive: 즉시 전환이 최선 (meta-RPE +{d['Passive_mRPE']})."),
           T(f"Lines = model predictions averaged over mice; dots = observed mean \u00b1 SEM (3-min blocks). Learning rules and "
             f"bounds as in Hjort et al.; readout logistic(b0 + \u03b2\u00b7Q + time in session + first 3 min) per mouse, solved "
             f"exactly for each candidate learning rate. BIC with n = all bins of the group (Active {nA}, Passive {nP}). PH and "
             f"no learning: \u0394BIC \u2265 70 (not shown).",
             f"선 = 마우스 평균 모델 예측; 점 = 실제 평균 \u00b1 SEM (3분 블록). 학습 규칙과 범위는 Hjort et al.과 같음; 행동 출력 "
             f"logistic(b0 + \u03b2\u00b7Q + 세션 내 시간 + 처음 3분)은 마우스별로, 후보 학습률마다 정확히 풂. BIC의 n = 그룹 "
             f"전체 bin 수 (Active {nA}, Passive {nP}). PH와 no learning은 \u0394BIC \u2265 70 (그림에 없음)."),
           T(notes4_en(d, nA, nP), notes4_kr(d, nA, nP)))
    return s


def notes4_en(d, nA, nP):
    pa = G5["Active"]["mRPE"]["params"]
    p2 = G5["Active"]["RPE2a"]["params"]
    a_to = d["Active_today's outcome"]
    return (
        "WHAT: the 8 learning rules fitted per group: learning rates shared by the mice of a group, readout (b0, beta, time "
        "terms) per mouse. For each candidate learning rate the readout is an exact logistic regression (profile "
        f"likelihood), so only the 1-4 learning parameters are searched. BIC with n = {nA} (Active) / {nP} (Passive) bins.\n"
        "a: dots = observed, lines = mean model predictions. Orange (RPE 2alpha) and green (meta-RPE) fall slowly on day 14 "
        "and rise fast on day 17; the dashed purple line (today's outcome = switch at once to the day's outcome) predicts "
        "day 14 far too low.\n"
        f"c: dBIC (0 = best, >10 = strong). Active: RPE 2alpha 0, meta-RPE {d['Active_mRPE']}, RPE {d['Active_RPE']}, "
        f"RPE t-1 {d['Active_RPEt-1']}, eMack {d['Active_eMack']}, today's outcome {a_to}. "
        f"Passive: today's outcome 0, RPE {d['Passive_RPE']}, RPE t-1 {d['Passive_RPEt-1']}, RPE 2alpha {d['Passive_RPE2a']}, "
        f"meta-RPE {d['Passive_mRPE']}.\n"
        f"Active group parameters: RPE 2alpha alpha- = {p2[0]:.2f}, alpha+ = {p2[1]:.2f}; meta-RPE alpha = {pa[0]:.2f}, "
        f"d = {pa[1]:.2f}, e = {pa[2]:.2f} (d at its lower bound, e at its upper bound: the same asymmetry, one parameter more).\n"
        "MEANING: Active mice learn slowly from negative errors (drug gone) and fast from positive errors (drug back). The "
        "error-averaging part of meta-RPE is not needed to explain it. Passive mice switch without a learning lag.")


def notes4_kr(d, nA, nP):
    pa = G5["Active"]["mRPE"]["params"]
    p2 = G5["Active"]["RPE2a"]["params"]
    a_to = d["Active_today's outcome"]
    return (
        "[이 슬라이드] 학습 규칙 8개를 그룹 단위로 맞춘 결과.\n"
        "[방법] 학습률은 그룹 안의 쥐들이 공유하고, 행동 출력(b0, β, 시간 항)은 마우스별. 후보 학습률마다 행동 출력은 "
        f"로지스틱 회귀로 정확히 풀고(profile likelihood) 학습률 1–4개만 탐색했다. BIC의 n = Active {nA}, Passive {nP} bin.\n"
        "[a 읽는 법] 점 = 실제 (평균 ± SEM), 선 = 모델 예측의 마우스 평균. 주황(RPE 2α)과 초록(meta-RPE)은 day 14에 천천히 "
        "내려가고 day 17에 빨리 올라간다. 보라 점선(today's outcome = 그날 결과로 바로 전환)은 day 14를 너무 낮게 예측한다.\n"
        f"[c 읽는 법] ΔBIC (0 = 최선, 10 이상 = 강한 증거). Active: RPE 2α 0, meta-RPE {d['Active_mRPE']}, RPE {d['Active_RPE']}, "
        f"RPE t−1 {d['Active_RPEt-1']}, eMack {d['Active_eMack']}, today's outcome {a_to}. "
        f"Passive: today's outcome 0, RPE {d['Passive_RPE']}, RPE t−1 {d['Passive_RPEt-1']}, RPE 2α {d['Passive_RPE2a']}, "
        f"meta-RPE {d['Passive_mRPE']}.\n"
        f"[파라미터] Active RPE 2α: α− = {p2[0]:.2f}, α+ = {p2[1]:.2f}. meta-RPE: α = {pa[0]:.2f}, d = {pa[1]:.2f}, e = {pa[2]:.2f} "
        "(d는 하한, e는 상한에 붙음 = 같은 비대칭을 파라미터 하나 더 써서 표현).\n"
        "[의미] Active는 '음의 오류(약이 없어짐)에서는 천천히, 양의 오류(약이 돌아옴)에서는 빨리' 배운다. meta-RPE도 같은 "
        "비대칭을 보이지만 파라미터가 더 많아 BIC에서 진다. 즉 오류를 '평균'하는 부분까지는 필요 없다. Passive는 학습 지연 "
        "없이 바로 전환한다.")


def slide5(prs, F):
    global SLIDE
    SLIDE = 5
    path, size, tg = F
    s = new_slide(prs)
    box = place(s, path, size, 0.30, 0.98, 12.70)
    conf = REC["confusion_BIC"]
    ra = reg("s5_rec_a", f"{conf['Active']['mRPE'].get('mRPE', 0) * 100:.0f}", conf["Active"]["mRPE"].get("mRPE", 0) * 100)
    rp = reg("s5_rec_p", f"{conf['Passive']['mRPE'].get('mRPE', 0) * 100:.0f}", conf["Passive"]["mRPE"].get("mRPE", 0) * 100)
    A = BS["Active"]
    am = reg("s5_am", f"{A['full_fit']['a_minus']:.2f}", A["full_fit"]["a_minus"])
    ap = reg("s5_ap", f"{A['full_fit']['a_plus']:.2f}", A["full_fit"]["a_plus"])
    fa = reg("s5_fa", f"{A['frac_a_minus_less_than_a_plus'] * 100:.1f}", A["frac_a_minus_less_than_a_plus"] * 100)
    fpv = reg("s5_fp", f"{BS['Passive']['frac_a_minus_less_than_a_plus'] * 100:.0f}",
              BS["Passive"]["frac_a_minus_less_than_a_plus"] * 100)
    ks = 14 if LANG == "EN" else 13
    key_box(s, box, tg, [P(Rn(T(f"Even when meta-RPE made the data, per-mouse BIC picked it in only {ra}% of sets.",
                                f"meta-RPE로 만든 데이터에서도 마우스별 BIC가 meta-RPE를 고른 건 {ra}%뿐."), ks, True,
                              "1F4E79"), after=2),
                         P(Rn(T("The paper's per-mouse test has no power with one switch each way.",
                                "방향마다 전환이 한 번뿐이라 논문식 마우스별 비교는 검정력이 없음."), ks))], tip=tg["tip"])
    common(s, T("The data show the asymmetry, not the error averaging", "비대칭은 보이지만, 오류 '평균'이 필요한지는 못 보인다"),
           T("refitted licking simulated from each mouse's models (20 sets each); bootstrapped the group learning rates.",
             "마우스마다 맞춘 모델로 만든 가짜 데이터(20세트)를 다시 맞춤; 그룹 학습률은 bootstrap."),
           T(f"per-mouse BIC finds meta-RPE in {ra}% / {rp}% of sets; Active \u03b1\u2212 {am} vs \u03b1+ {ap} "
             f"(\u03b1\u2212 < \u03b1+ in {fa}% of bootstraps; Passive {fpv}%).",
             f"마우스별 BIC가 meta-RPE를 찾은 비율 {ra}% / {rp}%; Active \u03b1\u2212 {am} vs \u03b1+ {ap} "
             f"(bootstrap {fa}%에서 \u03b1\u2212 < \u03b1+; Passive {fpv}%)."),
           T("a: licking simulated from each mouse's own fitted model with its real drop sequence; 3 candidate models refitted per "
             "set; 6 Active + 7 Passive mice x 20 sets x 3 generating models. b: RPE 2\u03b1 group fit (\u03b1\u2212 after negative, "
             "\u03b1+ after positive errors), point = fit to all mice, bar = 95% CI of 1000 bootstraps over mice; bounds 0.01\u20131. "
             "A per-mouse test needs several switches per mouse (e.g. repeated morphine \u2194 water, or 85% \u2192 15% morphine drops).",
             "a: 각 마우스의 맞춘 모델 + 실제 drop 순서로 lick을 시뮬레이션; 세트마다 후보 모델 3개를 다시 맞춤; Active 6 + Passive 7"
             "마리 x 20세트 x 생성 모델 3개. b: RPE 2\u03b1 그룹 fit (\u03b1\u2212 = 음의 오류 뒤, \u03b1+ = 양의 오류 뒤), 점 = 전체 "
             "마우스 fit, 막대 = 마우스 1000번 bootstrap의 95% CI; 범위 0.01\u20131. 마우스별 검정에는 마우스당 전환이 여러 번 필요."),
           T(notes5_en(ra, rp, am, ap, fa, fpv), notes5_kr(ra, rp, am, ap, fa, fpv)))
    return s


def notes5_en(ra, rp, am, ap, fa, fpv):
    conf = REC["confusion_BIC"]
    A, Pv = BS["Active"], BS["Passive"]
    tdo = "today's outcome"
    return (
        "WHAT: two checks. (a) Can a per-mouse model comparison like the paper's work on this data? (b) How sure is the learning "
        "asymmetry?\n"
        f"a: from each mouse's fitted models we simulated 20 licking data sets and refitted. When meta-RPE generated the data, "
        f"BIC picked meta-RPE in {ra}% (Active) / {rp}% (Passive) of sets; RPE-generated data were recognised in "
        f"{conf['Active']['RPE'].get('RPE', 0) * 100:.0f}% / {conf['Passive']['RPE'].get('RPE', 0) * 100:.0f}%, today's-outcome data in "
        f"{conf['Active'][tdo].get(tdo, 0) * 100:.0f}% / {conf['Passive'][tdo].get(tdo, 0) * 100:.0f}%. "
        "With one switch each way per mouse, the data cannot pay for meta-RPE's extra parameters, so 'meta-RPE lost' per mouse "
        "means nothing.\n"
        f"b: RPE 2alpha group fit, 1000 bootstraps over mice. Active alpha- {am} [{A['a_minus_CI'][0]:.2f}, {A['a_minus_CI'][2]:.2f}] "
        f"vs alpha+ {ap} [{A['a_plus_CI'][0]:.2f}, {A['a_plus_CI'][2]:.2f}]; alpha- < alpha+ in {fa}%. Passive alpha- "
        f"{Pv['full_fit']['a_minus']:.2f} vs alpha+ {Pv['full_fit']['a_plus']:.2f}, {fpv}% (unclear).\n"
        "MEANING: 'learn slowly that the drug is gone, fast that it is back' is solid in Active mice; whether errors are AVERAGED "
        "(the meta part) cannot be decided with this design.\n"
        "NEXT DESIGN: several switches per mouse (repeated morphine <-> water), probabilistic outcomes (e.g. 85% -> 15% morphine "
        "drops), or a cue that signals morphine vs water.")


def notes5_kr(ra, rp, am, ap, fa, fpv):
    conf = REC["confusion_BIC"]
    A, Pv = BS["Active"], BS["Passive"]
    tdo = "today's outcome"
    return (
        "[이 슬라이드] 두 가지 확인: (a) 논문처럼 마우스마다 모델을 비교하는 것이 이 데이터에서 가능한가? (b) 학습률 비대칭은 "
        "얼마나 확실한가?\n"
        f"[a 읽는 법] 각 마우스의 맞춘 모델로 가짜 lick 데이터를 20세트씩 만들고 다시 맞췄다. meta-RPE로 만든 데이터인데도 "
        f"BIC가 meta-RPE를 고른 비율은 {ra}% (Active), {rp}% (Passive). RPE로 만든 데이터는 "
        f"{conf['Active']['RPE'].get('RPE', 0) * 100:.0f}% / {conf['Passive']['RPE'].get('RPE', 0) * 100:.0f}%, today's outcome으로 "
        f"만든 데이터는 {conf['Active'][tdo].get(tdo, 0) * 100:.0f}% / {conf['Passive'][tdo].get(tdo, 0) * 100:.0f}% 맞게 찾았다. "
        "마우스당 전환이 방향마다 한 번뿐이라 정보가 부족해서, meta-RPE의 추가 파라미터를 정당화할 수 없다. 그래서 마우스별 "
        "비교에서 'meta-RPE가 졌다'는 아무 의미가 없다.\n"
        f"[b 읽는 법] RPE 2α 그룹 fit의 학습률과 1000번 bootstrap(마우스 재표집)의 95% CI. Active α−(모르핀 → 물) {am} "
        f"[{A['a_minus_CI'][0]:.2f}, {A['a_minus_CI'][2]:.2f}] vs α+(물 → 모르핀) {ap} [{A['a_plus_CI'][0]:.2f}, "
        f"{A['a_plus_CI'][2]:.2f}]; {fa}%에서 α− < α+. Passive α− {Pv['full_fit']['a_minus']:.2f} vs α+ "
        f"{Pv['full_fit']['a_plus']:.2f}, {fpv}% (불확실).\n"
        "[의미] '약이 없어진 것은 천천히, 돌아온 것은 빨리 배운다'는 Active에서 꽤 확실하다. 하지만 오류를 '평균'하는 meta "
        "부분이 꼭 필요한지는 이 설계로는 판단할 수 없다.\n"
        "[다음 설계] 마우스마다 전환을 여러 번 (모르핀 ↔ 물 반복), 확률적인 결과 (예: 모르핀 drop 85% → 15%), 또는 모르핀/물을 "
        "미리 알려 주는 큐가 있으면 마우스별 비교에 검정력이 생긴다.")


def slide6(prs, F):
    global SLIDE
    SLIDE = 6
    path, size, tg, info = F
    s = new_slide(prs)
    box = place(s, path, size, 0.25, 0.98, 9.25)
    ia, ip = info["Active"], info["Passive"]
    na = reg("s6_na", str(ia["cd_min_after"]), ia["cd_min_after"])
    nw = reg("s6_nw", str(ia["n_water"]), ia["n_water"])
    npv = reg("s6_np", str(ip["cd_min_after"]), ip["cd_min_after"])
    ks = 13.5 if LANG == "EN" else 12.5
    paras = [P(Rn(T("Prediction for imaging", "이미징에 대한 예측"), ks + 1.5, True, "1F4E79"), after=6),
             P(Rn(T("Active: ", "Active: "), ks, True, "D62728"),
               Rn(T(f"CD builds slowly; minimum after {na} of {nw} water drops (days 14\u201316). A 'CD cell' should ramp.",
                    f"CD가 천천히 쌓여 물 drop {nw}개 중 {na}번째에 최소 (day 14\u201316). 'CD 세포'라면 천천히 ramp."), ks), after=6),
             P(Rn(T("Passive: ", "Passive: "), ks, True, "1F77B4"),
               Rn(T(f"CD minimum {npv} drops after the switch. A 'CD cell' should step.",
                    f"전환 후 {npv}번째 drop에서 CD 최소. 'CD 세포'라면 계단처럼 바뀜."), ks), after=6),
             P(Rn(T("Test as in the paper: GLM with \u03b4, CD, CE and Q per drop, circular-shift null.",
                    "검정은 논문처럼: drop마다 \u03b4, CD, CE, Q를 넣은 GLM, circular shift로 p값."), ks, False, "444444"))]
    kx, ky, kw = 9.62, 1.30, 3.45
    tbox(s, kx, ky, kw, None, paras, "KEY", fill="F2F6FB", line="1F4E79")
    h = measure(paras, kw) + 0.06
    check(f"[{LANG}] slide 6 KEY ends above the bottom line", ky + h <= 5.85)
    # arrows leave the box's left edge at the middle of the Active / Passive paragraph
    tops = np.cumsum([0.05] + [measure([pp], kw) - 0.10 for pp in paras])
    for g, k in (("Active", 1), ("Passive", 2)):
        tx, ty = at(box, tg[g])
        # Active arrow from its paragraph; Passive arrow from the box's lower edge so it passes under panel titles
        ay = ky + (tops[k] + tops[k + 1]) / 2 if g == "Active" else ky + h - 0.15
        check(f"[{LANG}] slide 6 arrow {g} starts on the box edge", ky < ay < ky + h)
        arrow(s, kx, ay, tx + 0.06, ty)
    pa, pp_ = G5["Active"]["mRPE"]["params"], G5["Passive"]["mRPE"]["params"]
    common(s, T("Model variables per drop: regressors for ORBm / BMAp neurons", "drop마다 나오는 모델 변수: ORBm / BMAp 뉴런 회귀변수 후보"),
           T("computed RPE \u03b4, CD, CE, gain M and value Q at every drop with the group-fitted meta-RPE (one yoked pair).",
             "그룹에 맞춘 meta-RPE로 drop마다 RPE \u03b4, CD, CE, gain M, 가치 Q를 계산함 (yoked pair 하나)."),
           T(f"Active CD reaches its minimum {na} drops after the switch, Passive after {npv}: two time courses to look for in "
             f"neurons.",
             f"Active CD는 전환 후 {na}번째, Passive는 {npv}번째 drop에서 최소: 뉴런에서 찾아볼 두 가지 시간 경과."),
           T(f"Pair 4: Active 0910 black, Passive 0910 orange. Group meta-RPE parameters (Active \u03b1 = {pa[0]:.2f}, d = "
             f"{pa[1]:.2f}, e = {pa[2]:.2f}; Passive \u03b1 = {pp_[0]:.2f}, d = {pp_[1]:.2f}, e = {pp_[2]:.2f}); value starts at 0 on "
             f"day 3. Only CD < 0 and CE > 0 raise the gain M. A prediction, not a test (Hjort et al. Fig. 2: 13% of mPFC cells "
             f"encode CD).",
             f"Pair 4: Active 0910 black, Passive 0910 orange. 그룹 meta-RPE 파라미터 (Active \u03b1 = {pa[0]:.2f}, d = {pa[1]:.2f}, "
             f"e = {pa[2]:.2f}; Passive \u03b1 = {pp_[0]:.2f}, d = {pp_[1]:.2f}, e = {pp_[2]:.2f}); 가치는 day 3에 0에서 시작. "
             f"CD < 0, CE > 0일 때만 gain M이 커짐. 검정이 아닌 예측 (Hjort et al. Fig. 2: mPFC 뉴런의 13%가 CD 표현)."),
           T(notes6_en(na, nw, npv, pa, pp_), notes6_kr(na, nw, npv, pa, pp_)))
    return s


def notes6_en(na, nw, npv, pa, pp_):
    return (
        "WHAT: the variables the model computes at every drop - what a neural GLM would use once ORBm / BMAp calcium imaging "
        "exists (as in the paper's Fig. 2, with a circular-shift null for p values).\n"
        "a, c: black = value Q (expected drug), grey bars = RPE delta per drop; from day 14 the drops are water, so delta < 0.\n"
        "b, d: red CD = running mean of errors (slow), blue CE (fast), green dashed = gain M (sets how much is learned). Only "
        "CD < 0 or CE > 0 raises M.\n"
        f"PREDICTION: Active 0910 black: CD builds over days 14-16 and is lowest after {na} of {nw} water drops -> a CD cell would "
        f"ramp slowly. Passive 0910 orange: CD is lowest {npv} drops after the switch -> a CD cell would step.\n"
        f"Parameters (group fits): Active alpha = {pa[0]:.2f}, d = {pa[1]:.2f}, e = {pa[2]:.2f}; Passive alpha = {pp_[0]:.2f}, "
        f"d = {pp_[1]:.2f}, e = {pp_[2]:.2f}. Paper: 13% of mPFC neurons and 61% of mPFC->VTA neurons encode CD.")


def notes6_kr(na, nw, npv, pa, pp_):
    return (
        "[이 슬라이드] 모델이 drop마다 계산하는 변수들. 나중에 ORBm / BMAp 칼슘 이미징을 하면 GLM 회귀변수로 쓸 것이다 "
        "(논문 Fig. 2와 같은 방식, circular shift로 p값).\n"
        "[a, c 읽는 법] 검은 선 = 가치 Q (기대 약물), 회색 막대 = drop마다의 오류 δ. day 14부터 물 drop이라 δ가 음수.\n"
        "[b, d 읽는 법] 빨강 CD = 오류의 이동평균 (느림), 파랑 CE (빠름), 초록 점선 = gain M (한 번에 얼마나 배울지 결정). "
        "CD < 0이거나 CE > 0일 때만 M이 커진다.\n"
        f"[예측] Active (0910 black): CD가 day 14–16 동안 천천히 쌓여 물 drop {nw}개 중 {na}번째에 최소 → 'CD 세포'가 있다면 "
        f"천천히 ramp. Passive (0910 orange): 전환 후 {npv}번째 drop에서 CD 최소 → 계단처럼 빠르게 바뀜.\n"
        f"[주의] 검정이 아니라 예측이다. 파라미터는 그룹 fit (Active α = {pa[0]:.2f}, d = {pa[1]:.2f}, e = {pa[2]:.2f}; Passive "
        f"α = {pp_[0]:.2f}, d = {pp_[1]:.2f}, e = {pp_[2]:.2f}). 논문에서는 mPFC 뉴런의 13%, mPFC→VTA 뉴런의 61%가 CD를 표현했다.")


# ============================================================================================ build + render
def render(pptx, outdir):
    outdir.mkdir(parents=True, exist_ok=True)
    LO_PROFILE.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        r = subprocess.run([SOFFICE, f"-env:UserInstallation={LO_PROFILE.as_uri()}", "--headless", "--convert-to", "pdf",
                            "--outdir", td, str(pptx)], capture_output=True, text=True, timeout=300)
        pdf = Path(td) / (pptx.stem + ".pdf")
        if not pdf.exists():
            raise RuntimeError(f"LibreOffice failed: {r.stdout} {r.stderr}")
        doc = pymupdf.open(pdf)
        paths = []
        for i, pg in enumerate(doc):
            out = outdir / f"slide{i + 1}.png"
            pg.get_pixmap(dpi=110).save(out)
            paths.append(out)
        doc.close()
        shutil.copy(pdf, outdir / (pptx.stem + ".pdf"))
    return paths


def main():
    global LANG
    PRED = group_predictions()
    FIGS = [fig_paper(), fig_map(), fig_behaviour(), fig_models(PRED), fig_recovery(), fig_vars()]
    BUILD = [slide1, slide2, slide3, slide4, slide5, slide6]
    NAMES = ["1_paper_fig1", "2_mapping", "3_behaviour_switches", "4_group_model_fits", "5_recovery_asymmetry",
             "6_model_variables"]
    decks = {}
    for LANG in ("EN", "KR"):
        allp = Presentation(str(TEMPLATE))
        for fn, F, nm in zip(BUILD, FIGS, NAMES):
            one = Presentation(str(TEMPLATE))
            fn(one, F)
            one.save(OUT / "single_slides" / f"{LANG}_slide{nm}_20261117.pptx")
            fn(allp, F)
        decks[LANG] = OUT / f"mrpe_hjort2026_opioid_6slides_{LANG}_20261117.pptx"
        allp.save(decks[LANG])
    for k, (txt, val) in PRINTED.items():
        t = txt.replace("\u2212", "-").replace("+", "").strip()
        if t.startswith("<"):
            ok = val < float(t[1:])
        else:
            nd = len(t.split(".")[1]) if "." in t else 0
            ok = abs(float(t) - val) <= 0.5 * 10 ** -nd + 1e-12
        check(f"printed {k} = {txt} matches {val:.6g}", ok)
    json.dump(dict(checks=CHECKS, fit=FIT), open(OUT / "fit_report.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    bad = [c for c in CHECKS if not c[1]]
    print(f"{len(CHECKS)} checks, {len(bad)} failed")
    for c in bad:
        print("  FAILED", c[0])
    for lang, deck in decks.items():
        print(lang, deck, [p.name for p in render(deck, OUT / f"preview_{lang}")])


if __name__ == "__main__":
    main()
