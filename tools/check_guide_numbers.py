"""Check that the key numbers quoted in README.md and docs/*.md equal the values in project/data/.

Each check formats a value read from the result files exactly as the guide prints it and asserts that the resulting
text occurs in the named document (whitespace collapsed, '**' bold markers removed). Run from the package root:
    python tools/check_guide_numbers.py        -> prints ALL OK or the failing checks
"""
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "project" / "data"
J = lambda n: json.load(open(D / n, encoding="utf-8"))
m03, m07, m04, m05g, m05r = (J("m03_paper_data_summary.json"), J("m07_fig1_numbers.json"), J("m04_summary.json"),
                             J("m05_group_fits.json"), J("m05_recovery.json"))
m08, m11, m11cv, m11r = (J("m08_bootstrap_asymmetry.json"), J("m11_modelfree.json"), J("m11_lomo_cv.json"),
                         J("m11_group_recovery.json"))
m13, m15 = J("m13_summary.json"), J("m15_summary.json")
PM = pd.read_csv(D / "m13_per_mouse.csv")
BINS = pd.read_csv(D / "m13_bins.csv")

DOCS = {}
for p in [ROOT / "README.md"] + sorted((ROOT / "docs").glob("*.md")):
    t = p.read_text(encoding="utf-8").replace("**", "")
    DOCS[p.name] = re.sub(r"\s+", " ", t)
FAIL, N = [], [0]


def s2(x): return f"{x:+.2f}"
def s3(x): return f"{x:+.3f}"
def f0(x): return f"{x:.0f}"
def f1(x): return f"{x:.1f}"
def f2(x): return f"{x:.2f}"
def f3(x): return f"{x:.3f}"
def f4(x): return f"{x:.4f}"
def pc0(x): return f"{100 * x:.0f}%"
def pc1(x): return f"{100 * x:.1f}%"
def p2(x):            # p value as printed in the guide: 2 decimals from 0.095, else 2 significant digits
    if x >= 0.095:
        return f"{x:.2f}"
    return f"{x:.2g}" if x >= 0.001 else f"{x:.1g}"


def ptab(p):         # p in tables: 2 significant digits; 4 decimals below 0.001; '< 0.0001' below that
    if p >= 0.001:
        return f"{p:.2g}"
    return f"{p:.4f}" if p >= 0.0001 else "< 0.0001"


def chk(doc, text):
    N[0] += 1
    t = re.sub(r"\s+", " ", text)
    if t not in DOCS[doc]:
        FAIL.append((doc, t))


A13 = {g: m13[g]["A_observed"] for g in ("Active", "Passive")}
MF = m11
M15 = {g: m15[g]["models"] for g in ("Active", "Passive")}
MAIN = ["today's outcome", "RPE", "RPE2a", "mRPE"]
mt = lambda g, m, k="day14": M15[g][m]["A"][k]["mouse_t"]

# ============================================================ README
others = [m03["BIC"]["mean_change_from_best"][k] for k in ("RPE", "RPE2a", "RPEt-1", "PH", "eMack")]
chk("README.md", f"BIC change from best of {f1(m03['BIC']['mean_change_from_best']['mRPE'])} versus {f1(min(others))} to {f1(max(others))}")
w = MF["day14_vs_13"]["Active"]["whole"]
chk("README.md", f"{s2(w['mean'])} (95% CI {s2(w['lo'])} to {s2(w['hi'])}, p = {f2(w['p'])})")
chk("README.md", f"Passive licking rose by {s2(MF['day14_vs_13']['Passive']['whole']['mean'])}")
chk("README.md", f"fell on days 15 and 16 ({s2(MF['day15_vs_13']['Active']['whole']['mean'])}, "
                 f"{s2(MF['day16_vs_13']['Active']['whole']['mean'])}) and rose on the first morphine day "
                 f"({s2(MF['day17_vs_16']['Active']['whole']['mean'])})")
chk("README.md", f"alpha- = {f2(m08['Active']['full_fit']['a_minus'])} (its lower bound) and alpha+ = "
                 f"{f2(m08['Active']['full_fit']['a_plus'])}, with alpha- < alpha+ in {pc1(m08['Active']['frac_a_minus_less_than_a_plus'])}")
cv = m11cv["Active"]["delta_from_best"]
chk("README.md", f"+{f1(cv['mRPE'])} for meta-RPE, +{f1(cv['RPE'])} for one-rate RPE and +{f1(cv['today' + chr(39) + 's outcome'])} for instant switch")
chk("README.md", f"only {pc0(m05r['confusion_BIC']['Active']['mRPE']['mRPE'])} of per-mouse simulations and "
                 f"{pc0(m11r['confusion']['Active']['mRPE']['mRPE'])} of group simulations")
chk("README.md", f"no clear fall in licking ({s2(A13['Active']['d14_h1']['mean'])}), although "
                 f"{pc0(m13['Active']['day14_drops']['frac_by_2p5min'][2])} of the day's water drops")
chk("README.md", f"The whole-day change was {s2(A13['Active']['day14']['mean'])}. The models predicted "
                 f"{s2(M15['Active']['mRPE']['A']['day14']['pred'])} (meta-RPE) to "
                 f"{s2(M15['Active'][MAIN[0]]['A']['day14']['pred'])} (instant switch)")
chk("README.md", f"instant switch p = {p2(mt('Active', MAIN[0])['p_t'])}, RPE {p2(mt('Active', 'RPE')['p_t'])}, "
                 f"RPE 2α {f2(mt('Active', 'RPE2a')['p_t'])}, meta-RPE {p2(mt('Active', 'mRPE')['p_t'])}")
a, n_ = A13["Active"]["within14_halves"], A13["Active"]["night_halves"]
chk("README.md", f"Within day 14 the change was {s2(a['mean'])} (95% CI {s2(a['lo'])} to {s2(a['hi'])}) and overnight "
                 f"{s2(n_['mean'])} (95% CI {s2(n_['lo'])} to {s2(n_['hi'])})")
gc = m13["Active"]["sig_obs"]["groupcurve"]
chk("README.md", f"the midpoint of the morphine-to-water change is trial {f1(gc['deg']['mid'])} and of the water-to-morphine change trial {f1(gc['enh']['mid'])}")
chk("README.md", f"per-mouse test gives p = {f2(m13['Active']['sig_tests']['deg_vs_enh_mid_all']['p_t'])} (n = 6)")
chk("README.md", f"(trial {f1(gc['deg']['start'])}, 95% CI 0 to {f0(m13['Active']['groupcurve_boot']['deg_start_CI'][2])})")
pl_p = [M15["Passive"][m]["placebo"]["whole"]["mouse_t"]["p_t"] for m in M15["Passive"]]
chk("README.md", f"(day 13 minus day 12: {s2(M15['Passive']['RPE2a']['placebo']['whole']['obs'])})")
chk("README.md", f"(p = {f3(min(pl_p))} to {f3(max(pl_p))})")
Pb = BINS[BINS.group == "Passive"]
obs_1412 = (Pb[Pb.day == 14].groupby("mouse").lick.mean() - Pb[Pb.day == 12].groupby("mouse").lick.mean()).mean()
pred_1412 = [(Pb[Pb.day == 14].groupby("mouse")["p_" + m].mean() - Pb[Pb.day == 12].groupby("mouse")["p_" + m].mean()).mean()
             for m in MAIN]
chk("README.md", f"the Passive rise is {s2(obs_1412)}, smaller than the {s2(min(pred_1412))} to {s2(max(pred_1412))} the main models")

# ============================================================ docs/01
bic = m03["BIC"]["mean_change_from_best"]
for nm, k in (("meta-RPE", "mRPE"), ("RPE", "RPE"), ("RPE 2α", "RPE2a"), ("RPE t-1", "RPEt-1"), ("PH", "PH"), ("eMack", "eMack")):
    chk("01_paper_and_model.md", f"| {nm} | {f2(bic[k])} |")
chk("01_paper_and_model.md", f"largest p = {f4(max(v['p_t_holm'] for v in m03['BIC']['vs_mRPE'].values()))}")
chk("01_paper_and_model.md", f"meta-RPE is the best model in {m03['negLL']['n_best']['mRPE']} of 40 mice")
pm_ = m03["mRPE_params_median"]
chk("01_paper_and_model.md", f"α {f2(pm_['alpha'])}, β {f2(pm_['beta'])}, d {f3(pm_['d_CD'])}, e {f3(pm_['e_CE'])}")
chk("01_paper_and_model.md", f"d < e in only {m03['mRPE_d_less_than_e_n']} of 40 mice")
L, R_, M_ = m07["licking (running mean, 5 trials)"], m07["RPE model value"], m07["meta-RPE model value"]
chk("01_paper_and_model.md", f"| Trials 0-2 after the reversal | {f2(L['trials_0_2'])} | {f2(R_['trials_0_2'])} | {f2(M_['trials_0_2'])} |")
chk("01_paper_and_model.md", f"| Trials 3-5 after the reversal | {f2(L['trials_3_5'])} | {f2(R_['trials_3_5'])} | {f2(M_['trials_3_5'])} |")
chk("01_paper_and_model.md", f"The gain M peaks {m07['gain_M_peak_trial_from_reversal']} trials after the reversal")

# ============================================================ docs/03 (beta signs)
bA = [v["beta"] for v in m04["Active"]["mRPE_params"].values()]
bP = [v["beta"] for v in m04["Passive"]["mRPE_params"].values()]
assert all(b > 0 for b in bA) and all(b < 0 for b in bP)
chk("03_mapping.md", f"β > 0 in all 6 Active mice ({f2(min(bA))} to {f2(max(bA))})")
chk("03_mapping.md", f"β < 0 in all 7 Passive mice ({f2(max(bP))} to {f2(min(bP))})")
assert all(v == 0 for g in ("Active", "Passive") for v in m04[g]["mRPE_steps_capped"].values())

# ============================================================ docs/04
rows = (("day 14 vs 13 (first water day)", "day14_vs_13"), ("day 15 vs 13", "day15_vs_13"), ("day 16 vs 13", "day16_vs_13"),
        ("day 17 vs 16 (first morphine day)", "day17_vs_16"), ("day 18 vs 16", "day18_vs_16"))
for lab, k in rows:
    def cell(g):
        w = MF[k][g]["whole"]
        return f"{s3(w['mean'])} [{s3(w['lo'])}, {s3(w['hi'])}], {ptab(w['p'])}"
    wps = ptab(MF[k]["Active_vs_Passive_whole"]["welch_p"])
    chk("04_model_fits.md", f"| {lab} | {cell('Active')} | {cell('Passive')} | {wps} |")
asy = MF["asymmetry_whole_session_abs_change_day17_vs_day14"]["Active"]
chk("04_model_fits.md", f"|day 14 change| = {f2(asy['abs_day14'])} and |day 17 change| = {f2(asy['abs_day17'])}")
chk("04_model_fits.md", f"{asy['n_day17_larger']} of 6 mice changed more on day 17")
chk("04_model_fits.md", f"paired t p = {f2(asy['p_t'])} (Wilcoxon {f2(asy['p_wilcoxon'])})")
pmb = {g: m04[g]["BIC"]["mean_change_from_best"] for g in ("Active", "Passive")}
for nm, k in (("instant switch (today's outcome)", "today's outcome"), ("RPE", "RPE"), ("RPE t-1", "RPEt-1"), ("RPE 2α", "RPE2a"),
              ("meta-RPE", "mRPE"), ("no learning", "no learning"), ("eMack", "eMack"), ("PH", "PH")):
    chk("04_model_fits.md", f"| {nm} | {f2(pmb['Active'][k])} | {f2(pmb['Passive'][k])} |")
c = m05r["confusion_BIC"]
chk("04_model_fits.md", f"| meta-RPE | {pc1(c['Active']['mRPE']['mRPE'])} ({pc1(c['Active']['mRPE']['RPE'])} → RPE, "
                        f"{pc1(c['Active']['mRPE'][chr(116) + 'oday' + chr(39) + 's outcome'])} → instant switch) | {pc1(c['Passive']['mRPE']['mRPE'])} |")
chk("04_model_fits.md", f"| RPE | {pc1(c['Active']['RPE']['RPE'])} | {pc1(c['Passive']['RPE']['RPE'])} |")
to = "today's outcome"
chk("04_model_fits.md", f"| instant switch | {pc1(c['Active'][to][to])} | {pc1(c['Passive'][to][to])} |")
chk("04_model_fits.md", f"the correlation is {f2(m05r['corr_d'])} for d and {f2(m05r['corr_e'])} for e")
for g, order in (("Active", [("RPE 2α", "RPE2a"), ("meta-RPE", "mRPE"), ("RPE", "RPE"), ("eMack", "eMack"), ("RPE t-1", "RPEt-1"),
                             ("instant switch", to), ("PH", "PH"), ("no learning", "no learning")]),
                 ("Passive", [("instant switch", to), ("RPE", "RPE"), ("RPE t-1", "RPEt-1"), ("RPE 2α", "RPE2a"),
                              ("meta-RPE", "mRPE"), ("eMack", "eMack"), ("PH", "PH"), ("no learning", "no learning")])):
    for nm, k in order:
        r = m05g[g][k]
        d = r["dBIC_from_best"]
        ds = "0" if d == 0 else f"+{d:.1f}"
        chk("04_model_fits.md", f"| {nm} | {ds} | {r['negLL']:.1f} |")
pa = m05g["Active"]
chk("04_model_fits.md", f"α- {f3(pa['RPE2a']['params'][0])} (lower bound), α+ {f3(pa['RPE2a']['params'][1])}")
chk("04_model_fits.md", f"α {f3(pa['mRPE']['params'][0])}, d {f3(pa['mRPE']['params'][1])}, e {f3(pa['mRPE']['params'][2])}")
chk("04_model_fits.md", f"| RPE | +{m05g['Active']['RPE']['dBIC_from_best']:.1f} | {pa['RPE']['negLL']:.1f} | α {f3(pa['RPE']['params'][0])} |")
pp = m05g["Passive"]
chk("04_model_fits.md", f"α- {f3(pp['RPE2a']['params'][0])}, α+ {f3(pp['RPE2a']['params'][1])}")
chk("04_model_fits.md", f"α {f2(pp['mRPE']['params'][0])}, d {f3(pp['mRPE']['params'][1])}, e {f3(pp['mRPE']['params'][2])}")
bo = m05g["bootstrap_over_mice"]
chk("04_model_fits.md", f"Active: d < e in {pc1(bo['Active']['frac_d_less_than_e'])}; meta-RPE beats RPE by BIC in {pc0(bo['Active']['frac_mRPE_better_BIC'])}")
chk("04_model_fits.md", f"Passive: d < e in {pc1(bo['Passive']['frac_d_less_than_e'])}; meta-RPE beats RPE in {pc1(bo['Passive']['frac_mRPE_better_BIC'])}")
cvA, cvP = m11cv["Active"], m11cv["Passive"]
chk("04_model_fits.md", f"| Active (best: RPE 2α, {cvA['total_test_negLL']['RPE2a']:.2f}) | Δ | Passive (best: instant switch, {cvP['total_test_negLL'][to]:.2f}) | Δ |")
for (na, ka), (np_, kp) in zip((("meta-RPE", "mRPE"), ("RPE", "RPE"), ("RPE t-1", "RPEt-1"), ("instant switch", to), ("no learning", "no learning")),
                               (("RPE t-1", "RPEt-1"), ("meta-RPE", "mRPE"), ("RPE", "RPE"), ("RPE 2α", "RPE2a"), ("no learning", "no learning"))):
    extra = (f" (RPE 2α better in {cvA['vs_best']['RPE']['n_mice_best_better']}/6 mice; Wilcoxon p = {f2(cvA['vs_best']['RPE']['p_wilcoxon'])})"
             if ka == "RPE" else "")
    chk("04_model_fits.md", f"| {na} | +{cvA['delta_from_best'][ka]:.2f}{extra} | {np_} | +{cvP['delta_from_best'][kp]:.2f} |")
gr = m11r["confusion"]
chk("04_model_fits.md", f"| RPE 2α | RPE 2α {pc0(gr['Active']['RPE2a']['RPE2a'])}, RPE {pc0(gr['Active']['RPE2a']['RPE'])}, instant switch "
                        f"{pc0(gr['Active']['RPE2a'][to])} | RPE {pc0(gr['Passive']['RPE2a']['RPE'])}, instant switch {pc0(gr['Passive']['RPE2a'][to])} |")
chk("04_model_fits.md", f"| meta-RPE | RPE 2α {pc0(gr['Active']['mRPE']['RPE2a'])}, RPE {pc0(gr['Active']['mRPE']['RPE'])}, meta-RPE "
                        f"{pc0(gr['Active']['mRPE']['mRPE'])} | RPE {pc0(gr['Passive']['mRPE']['RPE'])}, instant switch {pc0(gr['Passive']['mRPE'][to])} |")
chk("04_model_fits.md", f"| RPE | RPE {pc0(gr['Active']['RPE']['RPE'])} | RPE {pc0(gr['Passive']['RPE']['RPE'])}, instant switch {pc0(gr['Passive']['RPE'][to])} |")
for g in ("Active", "Passive"):
    b = m08[g]
    lo_m, _, hi_m = b["a_minus_CI"]
    lo_p, _, hi_p = b["a_plus_CI"]
    chk("04_model_fits.md", f"| {g} | {f3(b['full_fit']['a_minus'])} [{f3(lo_m)}, {f3(hi_m) if hi_m < 1 else '1.00'}]")
    chk("04_model_fits.md", f"| {f3(b['full_fit']['a_plus'])} [{f3(lo_p)}, {f2(hi_p)}] | {pc1(b['frac_a_minus_less_than_a_plus'])}")
lr = m08["Active"]["log10_ratio_plus_over_minus_CI"]
chk("04_model_fits.md", f"log10(α+/α-) has a 95% interval of {f2(lr[0])} to {f1(lr[2])}")
chk("04_model_fits.md", f"Active minus Passive α-: [{f2(m08['Active_minus_Passive_a_minus_CI'][0])}, {f2(abs(m08['Active_minus_Passive_a_minus_CI'][2]))}]; "
                        f"Active is lower in {pc1(m08['frac_Active_a_minus_lower'])}")

# ============================================================ docs/05
dr = m13["Active"]["day14_drops"]
chk("05_trial_by_trial.md", f"Active mice earned {f1(dr['total_mean'])} water drops on day 14 (range {dr['total_min']} to {dr['total_max']})")
chk("05_trial_by_trial.md", f"{pc0(dr['frac_by_2p5min'][0])} of them came by 2.5 min and {pc0(dr['frac_by_2p5min'][2])} by 7.5 min")
fr = m15["Active"]["drops_frac_by_2p5min"]
chk("05_trial_by_trial.md", f"By 2.5 min, {pc0(fr['13'])} of the drops had come on day 13 (morphine), {pc0(fr['15'])} on day 15 and {pc0(fr['17'])} on day 17")
pl = m15["Active"]["p_lick_bins0_3"]
chk("05_trial_by_trial.md", f"P(lick) in the first 2 min (bins 0-3) is {f1(pl['13'])} on days 13 and 14")
assert pl["14"] == pl["13"]
chk("05_trial_by_trial.md", f"On days 16 and 17 it is {f2(pl['16'])}")
assert abs(pl["16"] - pl["17"]) < 1e-9
for lab, k in (("`day14`", "day14"), ("`d14_h1` (0-7.5 min)", "d14_h1"), ("`d14_h2` (7.5-15 min)", "d14_h2"), ("`d15_h1`", "d15_h1"),
               ("`day15`", "day15"), ("`day16`", "day16"), ("`within14_halves`", "within14_halves"), ("`night_halves`", "night_halves")):
    a = A13["Active"][k]
    chk("05_trial_by_trial.md", f"| {lab} | {s3(a['mean'])} [{s2(a['lo'])}, {s2(a['hi'])}] | {p2(a['p'])} | {a['n_neg']} / {a['n_pos']} |")
a5, n5 = A13["Active"]["within14_5min"], A13["Active"]["night_5min"]
chk("05_trial_by_trial.md", f"| `within14_5min` / `night_5min` | {s3(a5['mean'])} / {s3(n5['mean'])} | {f2(a5['p'])} / {f2(n5['p'])} |")
# models table
PMA = PM[PM.group == "Active"]
from scipy import stats  # noqa: E402
for nm, k in (("instant switch", to), ("RPE", "RPE"), ("RPE 2α", "RPE2a"), ("meta-RPE", "mRPE")):
    d = PMA["obs_day14"] - PMA["pred_day14_" + k]
    h = stats.t.ppf(0.975, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
    A_ = M15["Active"][k]["A"]["day14"]
    mtt = A_["mouse_t"]
    assert abs(mtt["mean_diff"] - d.mean()) < 1e-9
    pA, pB, pn = A_["calibrated"]["p"], A_["calibrated_B"]["p"], A_["naive"]["p"]
    fmtp = lambda x: f"{x:.3f}" if x < 0.06 and x >= 0.01 else (f"{x:.4f}" if x < 0.001 else (f"{x:.3f}" if x < 0.01 else f"{x:.2g}"))
    chk("05_trial_by_trial.md", f"| {nm} | {s3(A_['pred'])} | {s3(d.mean())} [{s2(d.mean() - h)}, {s2(d.mean() + h)}] | "
                                f"{p2(mtt['p_t'])} | {f2(mtt['p_holm'])} | {fmtp(pA)} ({fmtp(pB)}) | {fmtp(pn)} |")
for nm, k in (("RPE t-1", "RPEt-1"), ("PH", "PH"), ("eMack", "eMack"), ("no learning", "no learning")):
    A_ = M15["Active"][k]["A"]["day14"]
    pr = A_["pred"]
    prs = "0.000" if abs(pr) < 5e-4 else s3(pr)
    fmtp = lambda x: f"{x:.3f}" if x < 0.06 and x >= 0.01 else (f"{x:.3f}" if x < 0.01 else f"{x:.2g}")
    chk("05_trial_by_trial.md", f"| {nm} | {prs} | {s3(A_['mouse_t']['mean_diff'])} | {p2(A_['mouse_t']['p_t'])} | - | "
                                f"{fmtp(A_['calibrated']['p'])} ({fmtp(A_['calibrated_B']['p'])}) | {fmtp(A_['naive']['p'])} |")
hol = [M15["Active"][k]["A"]["day14"]["mouse_t"]["p_holm"] for k in MAIN]
chk("05_trial_by_trial.md", f"no mouse-level p is below {f2(min(hol))}")
nl = M15["Active"]["no learning"]["A"]
chk("05_trial_by_trial.md", f"It fails days 15 and 16 (mouse-level p = {f3(nl['day15']['mouse_t']['p_t'])} and {f3(nl['day16']['mouse_t']['p_t'])})")
lp = [M15["Active"][k]["A"][d]["mouse_t"]["p_t"] for k in ("RPE", "RPE2a", "mRPE") for d in ("day15", "day16")]
chk("05_trial_by_trial.md", f"(mouse-level p {f2(min(lp))} to {f2(max(lp))})")
chk("05_trial_by_trial.md", f"The instant switch gives {f2(M15['Active'][to]['A']['day15']['mouse_t']['p_t'])} and {f2(M15['Active'][to]['A']['day16']['mouse_t']['p_t'])}")
m13n = [m13["Active"]["A_models"][k]["day14"]["ppc_p"] for k in MAIN]
chk("05_trial_by_trial.md", f"(instant switch < 0.001, RPE {f3(m13n[1])}, RPE 2α {f3(m13n[2])}, meta-RPE {f2(m13n[3])})")
assert m13n[0] < 0.001
h1 = {k: M15["Active"][k]["A"]["d14_h1"] for k in MAIN}
chk("05_trial_by_trial.md", f"The instant switch predicts {s3(h1[to]['pred'])} against the observed {s3(A13['Active']['d14_h1']['mean'])}")
chk("05_trial_by_trial.md", f"mouse-level p = {f3(h1[to]['mouse_t']['p_t'])}; noise-matched p {f3(h1[to]['calibrated']['p'])} (A) and {f3(h1[to]['calibrated_B']['p'])} (B)")
chk("05_trial_by_trial.md", f"RPE is borderline: {f3(h1['RPE']['mouse_t']['p_t'])}, and {f2(h1['RPE']['calibrated']['p'])} ({f3(h1['RPE']['calibrated_B']['p'])})")
chk("05_trial_by_trial.md", f"RPE 2α ({f2(h1['RPE2a']['mouse_t']['p_t'])}) and meta-RPE ({f2(h1['mRPE']['mouse_t']['p_t'])}) pass")
# Passive
AP = A13["Passive"]
chk("05_trial_by_trial.md", f"licking rose at once on day 14 ({s3(AP['day14']['mean'])}; all 7 mice up; first half {s3(AP['d14_h1']['mean'])})")
assert AP["day14"]["n_pos"] == 7
pp_ = [M15["Passive"][k]["A"]["day14"]["pred"] for k in M15["Passive"] if k != "no learning"]
chk("05_trial_by_trial.md", f"The models predict {s2(min(pp_))} (PH) to {s2(max(pp_))} (instant switch)")
mp_ = [M15["Passive"][k]["A"]["day14"]["mouse_t"]["p_t"] for k in M15["Passive"] if k != "no learning"]
chk("05_trial_by_trial.md", f"mouse-level p = {f3(min(mp_))} to {f3(max(mp_))}")
dm = BINS[BINS.group == "Passive"].groupby(["mouse", "day"]).lick.mean().unstack().mean()
chk("05_trial_by_trial.md", f"Passive licking fell from day 12 to day 13 ({f2(dm[12])} → {f2(dm[13])}, a change of {s2(dm[13] - dm[12])})")
chk("05_trial_by_trial.md", f"(p = {f3(min(pl_p))} to {f3(max(pl_p))}; section 5.4)")
chk("05_trial_by_trial.md", f"Measured from day 12, the day-14 rise is {s2(obs_1412)}, smaller than the {s2(min(pred_1412))} to {s2(max(pred_1412))}")
# calibration features (ranges over the 4 main models)
rg = lambda g, key, nd: "-".join(f"{x:.{nd}f}" for x in (min(M15[g][k]["calibration"][key] for k in MAIN),
                                                          max(M15[g][k]["calibration"][key] for k in MAIN)))
one = lambda g, key: f2(np.mean([M15[g][k]["calibration"][key] for k in MAIN]))
chk("05_trial_by_trial.md", f"| T1: variance of whole-session mean differences / independent-bin variance | {rg('Active', 'T1_obs', 2)} | {rg('Passive', 'T1_obs', 2)} |")
chk("05_trial_by_trial.md", f"| T2: the same for half-session means | {rg('Active', 'T2_obs', 2)} | {rg('Passive', 'T2_obs', 2)} |")
chk("05_trial_by_trial.md", f"| T3: lag-1 autocorrelation of the bin differences | {one('Active', 'T3_obs')} | {one('Passive', 'T3_obs')} |")
for g in ("Active", "Passive"):
    v = [M15[g][k]["calibration"]["T3_obs"] for k in MAIN]
    assert max(v) - min(v) < 0.005, v
t1A = [M15["Active"][k]["calibration"]["A"]["T1_sim"] for k in MAIN]
chk("05_trial_by_trial.md", f"overshoots T1 ({f1(min(t1A))}-{f1(max(t1A))} simulated vs {rg('Active', 'T1_obs', 2)} observed, Active, 4 main models)")
assert all(M15["Active"][k]["calibration"]["A"]["at_grid_edge"] for k in M15["Active"])
t3B = [M15["Active"][k]["calibration"]["B"]["T3_sim"] for k in MAIN]
chk("05_trial_by_trial.md", f"undershoots T3 ({f2(min(t3B))}-{f2(max(t3B))} simulated vs {one('Active', 'T3_obs')}, Active, 4 main models)")
assert all(M15["Passive"][k]["calibration"]["B"]["sig_a"] == 3.0 for k in M15["Passive"])


def rng_(g, src):
    v = [M15[g][k]["placebo"]["whole"][src]["p"] if src != "mouse_t" else M15[g][k]["placebo"]["whole"]["mouse_t"]["p_t"]
         for k in M15[g]]
    return min(v), max(v)


chk("05_trial_by_trial.md", f"| | Active (observed {s3(M15['Active']['RPE2a']['placebo']['whole']['obs'])}) | Passive (observed {s3(M15['Passive']['RPE2a']['placebo']['whole']['obs'])}) |")
lo, hi = rng_("Active", "naive"); chk("05_trial_by_trial.md", f"| naive | {f2(lo)}-{f2(hi)} (passes)")
lo, hi = rng_("Active", "calibrated"); chk("05_trial_by_trial.md", f"| noise-matched A | {f2(lo)}-{f2(hi)} (passes) | {f2(rng_('Passive', 'calibrated')[0])}-{f2(rng_('Passive', 'calibrated')[1])} (passes) |")
lo, hi = rng_("Active", "calibrated_B"); chk("05_trial_by_trial.md", f"| noise-matched B | {f2(lo)}-{f2(hi)} (passes) | {f2(rng_('Passive', 'calibrated_B')[0])}-{f2(rng_('Passive', 'calibrated_B')[1])} (fails) |")
lo, hi = rng_("Active", "mouse_t"); chk("05_trial_by_trial.md", f"| mouse-level t | {f2(lo)}-{f2(hi)} (passes) | {f3(rng_('Passive', 'mouse_t')[0])}-{f3(rng_('Passive', 'mouse_t')[1])} (fails or borderline) |")
assert rng_("Passive", "naive")[1] <= 0.0005 + 1e-12
# timing
bt = m13["Active"]["groupcurve_boot"]
chk("05_trial_by_trial.md", f"| Midpoint | trial {f1(gc['deg']['mid'])} (bootstrap 95% CI {f1(bt['deg_mid_CI'][0])} to {f1(bt['deg_mid_CI'][2])}) | "
                            f"trial {f1(gc['enh']['mid'])} ({f1(bt['enh_mid_CI'][0])} to {f1(bt['enh_mid_CI'][2])}) |")
chk("05_trial_by_trial.md", f"| Fitted amplitude A, slope s | A = {s2(gc['deg']['A'])}, s = {f1(gc['deg']['s'])} | A = {s2(gc['enh']['A'])}, s = {f1(gc['enh']['s'])} (a step) |")
chk("05_trial_by_trial.md", f"later in {bt['frac_deg_mid_later'] * 100:.2f}% of 2,000 bootstraps")
ew = m15["Active"]["equal_window"]
chk("05_trial_by_trial.md", f"(days 14-15 only): {f1(ew['group_deg60_mid'])} vs {f1(ew['group_enh_mid'])}")
loo = m15["Active"]["leave_one_out"]
chk("05_trial_by_trial.md", f"{f1(min(x['deg_mid'] for x in loo))} to {f1(max(x['deg_mid'] for x in loo))} vs {f1(min(x['enh_mid'] for x in loo))} to {f1(max(x['enh_mid'] for x in loo))}")
sp = {k: m13["Active"]["A_models"][k]["sig_pred_groupcurve"] for k in MAIN}
chk("05_trial_by_trial.md", "| Morphine → water midpoint | " + " | ".join(f1(sp[k]["deg_mid"]) for k in MAIN) + " |")
chk("05_trial_by_trial.md", "| Water → morphine midpoint | " + " | ".join(f1(sp[k]["enh_mid"]) for k in MAIN) + " |")
tm = {k: M15["Active"][k]["timing"] for k in MAIN}
pf = lambda x: f"{x:.3f}" if x < 0.1 else f"{x:.2f}"
chk("05_trial_by_trial.md", "p shown as A (B): " + ", ".join(
    f"{nm} {pf(tm[k]['g_deg_mid']['calibrated']['p'])} ({pf(tm[k]['g_deg_mid']['calibrated_B']['p'])})"
    for nm, k in (("instant switch", to), ("RPE", "RPE"), ("RPE 2α", "RPE2a"), ("meta-RPE", "mRPE"))))
so = m13["Active"]["sig_obs"]
chk("05_trial_by_trial.md", "| Morphine → water | " + ", ".join(f1(x) for x in so["deg"]["mid"]) + f" (mean {f1(np.mean(so['deg']['mid']))}) |")
chk("05_trial_by_trial.md", "| Water → morphine | " + ", ".join(f1(x) for x in so["enh"]["mid"]) + f" (mean {f1(np.mean(so['enh']['mid']))}) |")
st = m13["Active"]["sig_tests"]["deg_vs_enh_mid_all"]
chk("05_trial_by_trial.md", f"Paired t p = {f2(st['p_t'])}, Wilcoxon {f2(st['p_wilcoxon'])}; later in {st['n_a_greater']} of 6 mice. With equal 60-trial windows p = {f2(ew['paired_deg60_vs_enh']['p_t'])}")
assert all(abs(x - 0.3) < 1e-9 for x in so["deg"]["s"])
sP = m13["Passive"]["sig_tests"]["deg_vs_enh_mid_all"]
chk("05_trial_by_trial.md", f"Per mouse: {f1(sP['mean_a'])} vs {f1(sP['mean_b'])} (p = {f2(sP['p_t'])})")
gP = m13["Passive"]["sig_obs"]["groupcurve"]
chk("05_trial_by_trial.md", f"(group-curve midpoint {f1(gP['deg']['mid'])})")
chk("05_trial_by_trial.md", f"(midpoint {f1(gP['enh']['mid'])}) hits its amplitude bound (A = {f1(gP['enh']['A'])})")
chk("05_trial_by_trial.md", f"The Active group-curve start is trial {f1(gc['deg']['start'])}: bootstrap 95% CI 0 to {f1(bt['deg_start_CI'][2])}; "
                            f"leave one mouse out {f1(min(x['deg_start'] for x in loo))} to {f1(max(x['deg_start'] for x in loo))}")
for nm, k in (("instant switch", to), ("RPE", "RPE"), ("RPE 2α", "RPE2a"), ("meta-RPE", "mRPE")):
    t = tm[k]["g_deg_start"]
    pa_, pb_ = t["calibrated"]["p"], t["calibrated_B"]["p"]
    chk("05_trial_by_trial.md", f"| {nm} | {f1(t['cal_median'])} | 0 to {f1(t['calibrated']['hi'])} | {pf(pa_)} ({pf(pb_)}) |")
mins = min(min(tm[k]["g_deg_start"]["calibrated"]["p"], tm[k]["g_deg_start"]["calibrated_B"]["p"]) for k in MAIN)
chk("05_trial_by_trial.md", f"the smallest p is {f3(mins)}")
cs = m13["Active"]["C_start_vs_models"]
chk("05_trial_by_trial.md", f"(mean {f1(cs[to]['start_vs_model']['mean_a'])}) with noise-free model curves, which all start by trial "
                            f"{f1(max(cs[k]['start_vs_model']['mean_b'] for k in MAIN))} (Holm p = {f3(cs[to]['start_vs_model']['p_t_holm'])})")

# ============================================================ report
if FAIL:
    print(f"{len(FAIL)} of {N[0]} checks FAILED:")
    for d, t in FAIL:
        print(f"  [{d}] expected text: {t}")
    sys.exit(1)
print(f"ALL OK ({N[0]} checks)")
