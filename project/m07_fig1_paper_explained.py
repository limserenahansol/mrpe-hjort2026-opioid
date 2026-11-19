"""Figure 1: the paper's Fig. 1 idea redrawn from the authors' own data and code (ported in m02, checked in m03).
All 40 mice, H->L cue, aligned to the reversal (trial 0 = first 15% trial), each mouse normalized pre = 1, post = 0:
(a) licking vs the value curve of the RPE model and of the meta-RPE model (mean +/- SEM);
(b) inside the meta-RPE model: single-trial RPE (top) and its rolling averages CD, CE and the gain M (bottom);
(c) BIC change from the best model (mean +/- SEM), paired t vs meta-RPE with Holm correction.
Output: figs/fig1_paper_fig1_explained.png, data/m07_fig1_numbers.json
"""
import json
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
import m02_hjort_models as hm  # noqa: E402

D = HERE / "data"
FIGS = HERE / "figs"
plt.rcParams.update({"font.family": "Arial", "font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
C_RW, C_M, C_CD, C_CE = "#8a8a8a", "#2e7d32", "#c62828", "#1565c0"

F = pd.read_csv(D / "m03_paper_data_fits.csv")
F["params"] = F.params.apply(json.loads)
S = json.load(open(D / "m03_paper_data_summary.json"))


def rw_q(p, c, r):
    Q = [0.15, 0.85]
    out = []
    for t in range(len(c)):
        out.append(Q[1])
        Q[c[t]] += p[0] * (r[t] - Q[c[t]])
    return np.array(out)


L, QR, QM, TM = [], [], [], []
for a in sorted(F.animal.unique()):
    c = np.load(D / "hjort2026_sample_choice" / a / "HL_choice.npy", allow_pickle=False).astype(np.int64)
    r = np.load(D / "hjort2026_sample_choice" / a / "HL_outcome.npy", allow_pickle=False).astype(np.float64)
    pm = np.array(F[(F.animal == a) & (F.model == "mRPE")].params.iloc[0])
    pr = np.array(F[(F.animal == a) & (F.model == "RPE")].params.iloc[0])
    t = hm.sim_mrpe(pm, c, r, 0.15, 0.85)
    L.append(pd.Series(c.astype(float)).rolling(5, center=True, min_periods=2).mean().values)
    QR.append(rw_q(pr, c, r))
    QM.append(t[:, 0])
    TM.append(t)
L, QR, QM, TM = map(np.array, (L, QR, QM, TM))


def nz(v):  # per mouse: pre-reversal stable part (trials 20-59) = 1, post stable day (120-159) = 0
    pre = v[:, 20:60].mean(1, keepdims=True)
    post = v[:, 120:].mean(1, keepdims=True)
    return (v - post) / (pre - post)


keep = (L[:, 20:60].mean(1) - L[:, 120:].mean(1)) > 0.1     # mice whose licking changed at all
x = np.arange(160) - 60
win = (x >= -20) & (x <= 60)
nums = {"n_mice_normalized": int(keep.sum())}

fig = plt.figure(figsize=(16, 5.0))
gs = fig.add_gridspec(2, 3, width_ratios=[1.3, 1.2, 0.95], height_ratios=[1, 1.3], hspace=0.14, wspace=0.33)
a = fig.add_subplot(gs[:, 0])
a.axvline(0, color="k", lw=1, ls="--")
for v, col, lw, lab in ((nz(L)[keep], "k", 2.0, "licking (running mean, 5 trials)"),
                        (nz(QR)[keep], C_RW, 2.6, "RPE model value"),
                        (nz(QM)[keep], C_M, 2.8, "meta-RPE model value")):
    mu = np.nanmean(v, 0)
    se = np.nanstd(v, 0, ddof=1) / np.sqrt(keep.sum())
    a.fill_between(x[win], (mu - se)[win], (mu + se)[win], color=col, alpha=0.15, lw=0)
    a.plot(x[win], mu[win], color=col, lw=lw, label=lab)
    nums[lab] = {k: round(float(np.nanmean(v[:, 60 + lo:60 + hi])), 3) for k, (lo, hi) in
                 {"trials_0_2": (0, 3), "trials_3_5": (3, 6), "trials_6_9": (6, 10), "trials_10_19": (10, 20)}.items()}
a.axhline(1, color="#bbb", lw=0.8, ls=":")
a.axhline(0, color="#bbb", lw=0.8, ls=":")
a.set(xlabel="Trials from reversal (cue H\u2192L: 85% \u2192 15% reward)", ylabel="Normalized (pre-reversal = 1, post = 0)",
      ylim=(-0.25, 1.3), title=f"a  {keep.sum()} mice, authors' data: licking holds ~6 trials; RPE value\n"
                               "falls from trial 1, meta-RPE value holds ~3 trials, then falls steeply")
a.legend(fontsize=9, loc="lower left", frameon=False)
a.text(0.98, 0.97, "mean \u00b1 SEM over mice", transform=a.transAxes, ha="right", va="top", fontsize=8.5)

b1 = fig.add_subplot(gs[0, 1])
b1.axvline(0, color="k", lw=1, ls="--")
b1.axhline(0, color="k", lw=0.6)
d_mu = TM[:, :, 1].mean(0)
b1.plot(x[win], d_mu[win], color="#555", lw=2)
b1.fill_between(x[win], 0, d_mu[win], color="#bbb", alpha=0.6)
b1.set(ylabel="RPE  \u03b4", title="b  Inside the meta-RPE model (40 mice)\nsingle-trial RPE is largest right after the reversal")
b1.tick_params(labelbottom=False)
b2 = fig.add_subplot(gs[1, 1], sharex=b1)
b2.axvline(0, color="k", lw=1, ls="--")
b2.axhline(0, color="k", lw=0.6)
b2.plot(x[win], TM[:, :, 2].mean(0)[win], color=C_CD, lw=2.4, label="CD = rolling mean of \u03b4 (rate d)")
b2.plot(x[win], TM[:, :, 3].mean(0)[win], color=C_CE, lw=1.8, label="CE = rolling mean of \u03b4 (rate e)")
b2.plot(x[win], TM[:, :, 4].mean(0)[win], color=C_M, lw=2.6, ls="--", label="gain M = max(0, \u2212CD) + max(0, CE)")
b2.set(xlabel="Trials from reversal", ylabel="CD, CE, gain M")
lo_, hi_ = b2.get_ylim()
b2.set_ylim(lo_, hi_ + 0.45 * (hi_ - lo_))
b2.legend(fontsize=8.2, loc="lower left", frameon=False)
b2.text(0.98, 0.97, "learning step = α × δ × M\nM peaks ~10 trials after the reversal,\nonce the errors have piled up",
        transform=b2.transAxes, fontsize=9, va="top", ha="right")
mM = TM[:, :, 4].mean(0)
nums["gain_M_peak_trial_from_reversal"] = int(x[win][np.argmax(mM[win])])
nums["rpe_most_negative_trial"] = int(x[win][np.argmin(d_mu[win])])

d = fig.add_subplot(gs[:, 2])
names = hm.ORDER
mu = [S["BIC"]["mean_change_from_best"][n] for n in names]
se = [S["BIC"]["sem"][n] for n in names]
d.bar(range(6), mu, yerr=se, color=[C_M] + ["#9e9e9e"] * 5, capsize=3)
lab = {"RPE2a": "RPE\n2\u03b1", "RPEt-1": "RPE\nt\u22121", "eMack": "eMack", "PH": "PH", "RPE": "RPE", "mRPE": "meta-\nRPE"}
d.set_xticks(range(6), [lab[n] for n in names])
for i, n in enumerate(names[1:], start=1):
    p = S["BIC"]["vs_mRPE"][n]["p_t_holm"]
    d.text(i, mu[i] + se[i] + 0.4, "*" if p < 0.05 else "ns", ha="center", fontsize=12)
d.set_ylim(0, max(mu) + max(se) + 4)
d.set(ylabel="BIC change from the best model (lower = better)",
      title=f"c  {S['n_animals']} mice: meta-RPE fits best\n(our port of the authors' code)")
d.text(0.02, 0.97, "mean \u00b1 SEM; * paired t vs meta-RPE,\nHolm p < 0.05 (all five)", transform=d.transAxes,
       ha="left", va="top", fontsize=8.5)
fig.savefig(FIGS / "fig1_paper_fig1_explained.png", dpi=170, bbox_inches="tight")
json.dump(nums, open(D / "m07_fig1_numbers.json", "w"), indent=1)
print(json.dumps(nums, indent=1))
