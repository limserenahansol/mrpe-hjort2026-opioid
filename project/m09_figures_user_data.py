"""Figures 2-4 for the Dec-2025 opioid cohort (numbers come from m04-m08 outputs; nothing typed in).
fig2: behaviour around the two switches (model-free); fig3: models with group-shared learning rates + model recovery;
fig4: the model's trial-by-trial variables for one yoked pair (what would be regressed onto neurons) + asymmetry bootstrap.
"""
import json
import pickle
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import m04_fit_user_data as m4  # noqa: E402
import m05_recovery_and_group_fits as m5  # noqa: E402

D = HERE / "data"
FIGS = HERE / "figs"
plt.rcParams.update({"font.family": "Arial", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
CA, CP = "#d62728", "#1f77b4"              # Active red, Passive blue
CMOD = {"today's outcome": "#9467bd", "RPE": "#8a8a8a", "RPE2a": "#ff7f0e", "mRPE": "#2e7d32"}
LAB = {"today's outcome": "today's outcome (instant switch)", "RPE": "RPE (one learning rate)",
       "RPE2a": "RPE 2α (separate rates for −/+ errors)", "mRPE": "meta-RPE (paper)"}
X = pd.read_csv(D / "m06_switch_blocks.csv")
MF = json.load(open(D / "m06_modelfree.json"))
G = json.load(open(D / "m05_group_fits.json"))
REC = json.load(open(D / "m05_recovery.json"))
BS = json.load(open(D / "m08_bootstrap_asymmetry.json"))
NA, NP = len(m5.MICE["Active"]), len(m5.MICE["Passive"])


def pstr(p):
    return "p < 0.001" if p < 0.001 else f"p = {p:.3f}" if p < 0.1 else f"p = {p:.2f}"


def sgn(v):
    return "0.00" if abs(v) < 0.005 else f"{v:+.2f}"


def shade_days(ax, days, y0, y1, blocks=5):
    for i, d in enumerate(days):
        if 14 <= d <= 16:
            ax.axvspan(i * blocks - 0.5, (i + 1) * blocks - 0.5, color="#e8f1fb", zorder=0)
        ax.axvline(i * blocks - 0.5, color="#ddd", lw=0.6, zorder=0)


# ------------------------------------------------------------------ Fig 2: behaviour (model-free)
DAYS = list(range(11, 19))
fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8), gridspec_kw=dict(width_ratios=[2.0, 1, 1]))
a = ax[0]
shade_days(a, DAYS, 0, 1)
for g, col in (("Active", CA), ("Passive", CP)):
    Y = X[(X.group == g) & X.day.isin(DAYS) & (X.mouse != "6099_white")]
    W = Y.pivot_table(index="mouse", columns=["day", "block"], values="p")
    mu, se = W.mean(), W.std(ddof=1) / np.sqrt(W.notna().sum())
    xx = np.arange(len(mu))
    for i in range(len(DAYS)):
        s = slice(i * 5, i * 5 + 5)
        a.fill_between(xx[s], (mu - se).values[s], (mu + se).values[s], color=col, alpha=0.18, lw=0)
        a.plot(xx[s], mu.values[s], color=col, lw=2.2, marker="o", ms=3.5,
               label=f"{g} (n = {W.shape[0]})" if i == 0 else None)
a.set_xticks([i * 5 + 2 for i in range(len(DAYS))], [f"day {d}" for d in DAYS])
for i, d in enumerate(DAYS):
    a.text(i * 5 + 2, 1.04, "water" if 14 <= d <= 16 else "morphine", ha="center", fontsize=9,
           color="#1565c0" if 14 <= d <= 16 else "#7b1fa2")
a.set(ylim=(0, 1.1), ylabel="p(lick in a 30-s bin)", xlim=(-0.7, len(DAYS) * 5 - 0.3),
      title="a  Licking in 3-min blocks of each 15-min PR session (mean ± SEM over mice)")
a.legend(loc="lower left", fontsize=9, frameon=False)
for k, (sw, d0, d1, title) in enumerate((("day14_morphine_to_water", 13, 14, "b  Day 14: morphine → water"),
                                         ("day17_water_to_morphine", 16, 17, "c  Day 17: water → morphine"))):
    b = ax[1 + k]
    b.axhline(0, color="k", lw=0.8)
    b.axvspan(-0.4, 2.4, color="#f3f3f3", zorder=0)
    b.text(1.0, 0.8, "first 9 min", ha="center", fontsize=8.5, color="#666")
    for g, col, off in (("Active", CA, -0.08), ("Passive", CP, 0.08)):
        r = MF[sw][g]
        b.errorbar(np.arange(5) + off, r["mean_change_by_block"], yerr=r["sem_by_block"], color=col, lw=2, marker="o",
                   capsize=3, label=f"{g} (n = {r['n']})")
    ra, rp = MF[sw]["Active"], MF[sw]["Passive"]
    b.text(0.02, 0.03, f"first 9 min vs 0 (one-sample t):\nActive {sgn(ra['change_first9min'])}, {pstr(ra['p_t_first9min_vs0'])}\n"
           f"Passive {sgn(rp['change_first9min'])}, {pstr(rp['p_t_first9min_vs0'])}", transform=b.transAxes, fontsize=8.6)
    b.set_xticks(range(5), ["0-3", "3-6", "6-9", "9-12", "12-15"])
    b.set(xlabel="Minutes into the session", ylabel=f"Change in p(lick): day {d1} − day {d0}", ylim=(-0.75, 0.85),
          title=title)
    b.legend(fontsize=8.6, loc="upper right", frameon=False)
fig.tight_layout()
fig.savefig(FIGS / "fig2_your_data_switches.png", dpi=170)
plt.close(fig)


# ------------------------------------------------------------------ group-fit predictions
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


def predict(model, params, m):
    typ, r, c, tau, early, use, dayout = m5.ARR[m]
    code = m4.SPECS[model][0]
    q = m5.q_at_bins(code, np.asarray(params, float), typ, r, dayout)
    sel = (typ == 1) & (use == 1)
    Xa = np.column_stack([np.ones(len(typ)), tau, early, q])
    w = logit_fit(Xa[sel], c[sel].astype(float))
    A = m5.EV[m]
    out = A[sel].copy()
    out["p"] = 1 / (1 + np.exp(-(Xa[sel] @ w)))
    out["q"] = q[sel]
    return out


PRED = []
for g in ("Active", "Passive"):
    for model in ("today's outcome", "RPE", "RPE2a", "mRPE"):
        for m in m5.MICE[g]:
            P = predict(model, G[g][model]["params"], m)
            P["block"] = (P.bin // 6).astype(int)
            PRED.append(P.groupby(["day", "block"]).agg(obs=("c", "mean"), pred=("p", "mean")).reset_index()
                        .assign(mouse=m, group=g, model=model))
PRED = pd.concat(PRED)

# ------------------------------------------------------------------ Fig 3: models
fig = plt.figure(figsize=(16.5, 8.6))
gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.35, 1.0], hspace=0.42, wspace=0.28)
DD = list(range(13, 19))
for k, (g, col) in enumerate((("Active", CA), ("Passive", CP))):
    a = fig.add_subplot(gs[k, 0:2])
    shade_days(a, DD, 0, 1)
    Y = PRED[(PRED.group == g) & PRED.day.isin(DD)]
    obs = Y[Y.model == "RPE"].pivot_table(index="mouse", columns=["day", "block"], values="obs")
    mu, se = obs.mean(), obs.std(ddof=1) / np.sqrt(len(obs))
    xx = np.arange(len(mu))
    for i in range(len(DD)):
        s = slice(i * 5, i * 5 + 5)
        a.errorbar(xx[s], mu.values[s], yerr=se.values[s], color="k", lw=0, elinewidth=1.2, marker="o", ms=5,
                   capsize=2, label="observed (mean ± SEM)" if i == 0 else None)
        for model in ("RPE", "mRPE", "RPE2a", "today's outcome"):
            pr = Y[Y.model == model].pivot_table(index="mouse", columns=["day", "block"], values="pred").mean()
            a.plot(xx[s], pr.values[s], color=CMOD[model], lw=2.6 if model in ("mRPE", "RPE2a") else 1.8,
                   ls="-" if model != "today's outcome" else (0, (4, 2)), zorder=3,
                   label=LAB[model] if (i == 0 and k == 0) else None)
    a.set_xticks([i * 5 + 2 for i in range(len(DD))], [f"day {d}\n{'water' if 14 <= d <= 16 else 'morphine'}" for d in DD])
    dl = ", ".join(f"{nm} {G[g][mm]['dBIC_from_best']:.0f}" for mm, nm in
                   (("RPE2a", "RPE 2α"), ("mRPE", "meta-RPE"), ("RPE", "RPE"), ("today's outcome", "today's outcome")))
    a.set(ylim=(0, 1.05), ylabel="p(lick in a 30-s bin)", xlim=(-0.7, len(DD) * 5 - 0.3),
          title=f"{'ab'[k]}  {g} (n = {len(obs)}), 3-min blocks; group-level ΔBIC: {dl}")
    if k == 0:
        a.legend(fontsize=9, loc="lower center", bbox_to_anchor=(0.5, 1.13), frameon=False, ncol=5)
c = fig.add_subplot(gs[0, 2])
models = ["RPE2a", "mRPE", "eMack", "RPE", "RPEt-1", "today's outcome", "PH", "no learning"]
yy = np.arange(len(models))
for g, col, off in (("Active", CA, -0.2), ("Passive", CP, 0.2)):
    vals = [G[g][mm]["dBIC_from_best"] for mm in models]
    c.barh(yy + off, [min(v, 80) for v in vals], height=0.38, color=col, label=g)
    for yv, v in zip(yy + off, vals):
        c.text(min(v, 80) + 1, yv, f"{v:.0f}", va="center", fontsize=7.8, color=col)
c.set_yticks(yy, [{"RPE2a": "RPE 2α", "mRPE": "meta-RPE", "RPEt-1": "RPE t−1", "today's outcome": "today's outcome",
                   "no learning": "no learning"}.get(mm, mm) for mm in models])
c.invert_yaxis()
c.axvline(10, color="k", lw=0.8, ls=":")
c.set_xlim(0, 92)
c.set(xlabel="ΔBIC from the best model (lower = better; bars capped at 80)",
      title="c  Group-level model comparison\n(dotted line: ΔBIC 10 = strong evidence)")
c.legend(fontsize=9, frameon=False, loc="upper right")
d = fig.add_subplot(gs[1, 2])
conf = REC["confusion_BIC"]
gens = ["mRPE", "RPE", "today's outcome"]
for j, (g, col) in enumerate((("Active", CA), ("Passive", CP))):
    vals = [conf[g][gen].get(gen, 0.0) * 100 for gen in gens]
    d.bar(np.arange(3) + (j - 0.5) * 0.38, vals, width=0.38, color=col, label=g)
    for i, v in enumerate(vals):
        d.text(i + (j - 0.5) * 0.38, v + 1.5, f"{v:.0f}%", ha="center", fontsize=8.5)
d.set_xticks(range(3), ["meta-RPE", "RPE", "today's\noutcome"])
d.set(ylim=(0, 112), xlabel="Model that generated the simulated licking",
      ylabel="% of sets where BIC picks\nthe true model (per mouse)",
      title=f"d  Can a per-mouse comparison find meta-RPE?\n20 simulated sets per mouse; n = {NA} + {NP} mice")
d.legend(fontsize=9, frameon=False, loc="upper left")
fig.savefig(FIGS / "fig3_models_on_your_data.png", dpi=170, bbox_inches="tight")
plt.close(fig)

# ------------------------------------------------------------------ Fig 4: model variables for one yoked pair + asymmetry
fig = plt.figure(figsize=(16.5, 7.4))
gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1.45, 0.9], hspace=0.38, wspace=0.45)
PAIR = ("0910_black", "0910_orange")
for k, (m, g, col) in enumerate(((PAIR[0], "Active", CA), (PAIR[1], "Passive", CP))):
    A = m5.EV[m]
    typ, r, c_, tau, early, use, dayout = m5.ARR[m]
    lp = np.asarray(G[g]["mRPE"]["params"], float)
    # walk the drops with the group's meta-RPE parameters to get delta, CD, CE, M, Q per drop
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
        rows.append(dict(day=A.day.iloc[i], delta=delta, CD=md, CE=me, M=M, Q=Q))
    T = pd.DataFrame(rows)
    T = T[T.day >= 11].reset_index(drop=True)
    a = fig.add_subplot(gs[k, 0])
    xx = np.arange(len(T))
    for d in sorted(T.day.unique()):
        idx = T.index[T.day == d]
        if 14 <= d <= 16:
            a.axvspan(idx[0] - 0.5, idx[-1] + 0.5, color="#e8f1fb", zorder=0)
        a.axvline(idx[0] - 0.5, color="#ccc", lw=0.6)
        a.text((idx[0] + idx[-1]) / 2, 1.12, f"d{d}", ha="center", fontsize=8.5)
    a.bar(xx, T.delta, color="#bbb", width=1.0, label="RPE δ per drop")
    a.plot(xx, T.Q, color="k", lw=2, label="value Q (expected drug)")
    a.set(ylim=(-1.1, 1.25), ylabel="Q, δ", xlabel="Drops, in order (days 11-18)",
          title=f"{'ac'[k]}  {g} {m}: value and single-drop RPE")
    a.legend(fontsize=8.4, loc="lower left", frameon=False)
    b = fig.add_subplot(gs[k, 1], sharex=a)
    for d in sorted(T.day.unique()):
        idx = T.index[T.day == d]
        if 14 <= d <= 16:
            b.axvspan(idx[0] - 0.5, idx[-1] + 0.5, color="#e8f1fb", zorder=0)
        b.axvline(idx[0] - 0.5, color="#ccc", lw=0.6)
    b.axhline(0, color="k", lw=0.6)
    b.plot(xx, T.CD, color="#c62828", lw=2.2, label=f"CD (rate d = {lp[1]:.2f})")
    b.plot(xx, T.CE, color="#1565c0", lw=1.6, label=f"CE (rate e = {lp[2]:.2f})")
    b.plot(xx, T.M, color="#2e7d32", lw=2.2, ls="--", label="gain M")
    b.set(xlabel="Drops, in order (days 11-18)", ylabel="CD, CE, M",
          title=f"{'bd'[k]}  {g}: averaged errors (candidate regressors for neurons)")
    b.legend(fontsize=8.4, loc="lower left", frameon=False)
    if k == 0:
        b.text(0.02, 0.36, "only CD < 0 and CE > 0\nraise the gain M", transform=b.transAxes, ha="left", fontsize=8.2,
               color="#555")
c = fig.add_subplot(gs[:, 2])
ylab = []
for j, (g, col) in enumerate((("Active", CA), ("Passive", CP))):
    r_ = BS[g]
    for k2, (key, lab, mk) in enumerate((("a_minus", "α−\n(drug → water)", "v"),
                                         ("a_plus", "α+\n(water → drug)", "^"))):
        lo, med, hi = r_[f"{key}_CI"]
        y = j * 3 + k2
        c.errorbar(med, y, xerr=[[med - lo], [hi - med]], color=col, marker=mk, ms=9, capsize=4, lw=2.2)
        c.text(med, y - 0.32, f"{med:.2f} [{lo:.2f}, {hi:.2f}]", ha="center", fontsize=8.6, color=col)
        ylab.append((y, f"{g} {lab}"))
    c.text(0.0105, j * 3 + 1.75, f"α− < α+ in {r_['frac_a_minus_less_than_a_plus'] * 100:.1f}% of bootstraps",
           fontsize=8.8, color=col)
c.set_xscale("log")
c.set_xlim(0.008, 1.5)
c.set_ylim(5.4, -0.9)
c.set_yticks([y for y, _ in ylab], [t for _, t in ylab], fontsize=9)
c.set(xlabel="Learning rate (RPE 2α, group fit; log scale)",
      title="e  Slow learning when morphine is removed,\nfast when it returns (median [95% CI];\n1000 bootstraps over mice)")
fig.savefig(FIGS / "fig4_model_variables_and_asymmetry.png", dpi=170, bbox_inches="tight")
plt.close(fig)
print("saved fig2, fig3, fig4")
