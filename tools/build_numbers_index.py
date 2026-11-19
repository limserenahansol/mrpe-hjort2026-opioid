"""Build docs/11_numbers_index.md and docs/numbers_index.json from the result files in project/data/.

Every number is read from a JSON file; nothing is typed in. Run from the package root:
    python tools/build_numbers_index.py
"""
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "project" / "data"
J = {n: json.load(open(D / n, encoding="utf-8")) for n in (
    "m03_paper_data_summary.json", "m07_fig1_numbers.json", "m04_summary.json", "m05_group_fits.json",
    "m05_recovery.json", "m08_bootstrap_asymmetry.json", "m11_modelfree.json", "m11_lomo_cv.json",
    "m11_group_recovery.json", "m13_summary.json", "m15_summary.json")}
GROUPS = ("Active", "Passive")
MAIN = ["today's outcome", "RPE", "RPE2a", "mRPE"]
ALLM = ["today's outcome", "RPE", "RPE2a", "mRPE", "PH", "eMack", "RPEt-1", "no learning"]
NICE = {"today's outcome": "instant switch", "RPE": "RPE", "RPE2a": "RPE 2a", "mRPE": "meta-RPE", "PH": "PH",
        "eMack": "eMack", "RPEt-1": "RPE t-1", "no learning": "no learning"}
ROWS = []          # (section, label, value, file, path, note)


def get(f, path):
    x = J[f]
    for k in path:
        x = x[k]
    return x


def fmt(v):
    if isinstance(v, bool) or v is None:
        return str(v)
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if v != 0 and abs(v) < 0.001:
            return f"{v:.2e}"
        return f"{round(v, 4):g}" if abs(v) < 1e6 else f"{v:.4g}"
    if isinstance(v, list):
        return "[" + ", ".join(fmt(x) for x in v) + "]"
    if isinstance(v, dict):
        return "{" + ", ".join(f"{k}: {fmt(x)}" for k, x in v.items()) + "}"
    return str(v)


def add(section, label, f, path, note=""):
    ROWS.append((section, label, get(f, path), f, path, note))


# ---------------------------------------------------------------- port check
S = "1. Port check on the authors' 40 mice (m03, m07)"
f = "m03_paper_data_summary.json"
for m in ("mRPE", "RPE", "RPE2a", "RPEt-1", "PH", "eMack"):
    add(S, f"BIC change from best, {NICE[m]}", f, ["BIC", "mean_change_from_best", m])
for m in ("RPE", "RPE2a", "RPEt-1", "PH", "eMack"):
    add(S, f"paired t on BIC vs meta-RPE, Holm p, {NICE[m]}", f, ["BIC", "vs_mRPE", m, "p_t_holm"])
add(S, "negLL: mice where meta-RPE is best (of 40)", f, ["negLL", "n_best", "mRPE"])
add(S, "median meta-RPE parameters", f, ["mRPE_params_median"])
add(S, "mice with d < e (of 40)", f, ["mRPE_d_less_than_e_n"])
f = "m07_fig1_numbers.json"
for k in ("licking (running mean, 5 trials)", "RPE model value", "meta-RPE model value"):
    add(S, f"Fig. 1 redrawn: {k}", f, [k])
add(S, "gain M peak, trials after reversal", f, ["gain_M_peak_trial_from_reversal"])

# ---------------------------------------------------------------- model-free
S = "2. Model-free switch tests (m11; whole session; change vs reference day)"
f = "m11_modelfree.json"
for key in ("day14_vs_13", "day15_vs_13", "day16_vs_13", "day17_vs_16", "day18_vs_16"):
    for g in GROUPS:
        for q in ("mean", "lo", "hi", "p"):
            add(S, f"{key} {g} {q}", f, [key, g, "whole", q])
    add(S, f"{key} Active vs Passive Welch p", f, [key, "Active_vs_Passive_whole", "welch_p"])
for q in ("abs_day14", "abs_day17", "n_day17_larger", "p_t", "p_wilcoxon"):
    add(S, f"asymmetry (Active) {q}", f, ["asymmetry_whole_session_abs_change_day17_vs_day14", "Active", q])

# ---------------------------------------------------------------- per-mouse fits and recovery
S = "3. Per-mouse fits (m04) and per-mouse recovery (m05a)"
for g in GROUPS:
    add(S, f"{g} per-mouse BIC change from best (mean)", "m04_summary.json", [g, "BIC", "mean_change_from_best"])
    add(S, f"{g} meta-RPE readout beta per mouse", "m04_summary.json", [g, "mRPE_params"],
        "beta = 5th entry; see JSON")
    add(S, f"{g} meta-RPE steps capped", "m04_summary.json", [g, "mRPE_steps_capped"])
    for gen in ("mRPE", "RPE", "today's outcome"):
        add(S, f"{g} recovery: data made by {NICE[gen]} -> BIC picks", "m05_recovery.json", ["confusion_BIC", g, gen])
add(S, "recovery corr(d true, d fitted)", "m05_recovery.json", ["corr_d"])
add(S, "recovery corr(e true, e fitted)", "m05_recovery.json", ["corr_e"])

# ---------------------------------------------------------------- group fits
S = "4. Group fits (m05b), LOMO-CV and group recovery (m11), asymmetry bootstrap (m08)"
f = "m05_group_fits.json"
for g in GROUPS:
    for m in ALLM:
        add(S, f"{g} {NICE[m]}: dBIC from best", f, [g, m, "dBIC_from_best"])
        add(S, f"{g} {NICE[m]}: negLL", f, [g, m, "negLL"])
        add(S, f"{g} {NICE[m]}: learning parameters", f, [g, m, "params"])
    for q in ("d_CI", "e_CI", "frac_d_less_than_e", "frac_mRPE_better_BIC"):
        add(S, f"{g} bootstrap (meta-RPE) {q}", f, ["bootstrap_over_mice", g, q])
for g in GROUPS:
    add(S, f"{g} LOMO-CV held-out negLL, delta from best", "m11_lomo_cv.json", [g, "delta_from_best"])
    add(S, f"{g} LOMO-CV best model", "m11_lomo_cv.json", [g, "best"])
    add(S, f"{g} group recovery (50 sets per model)", "m11_group_recovery.json", ["confusion", g])
f = "m08_bootstrap_asymmetry.json"
for g in GROUPS:
    for q in ("full_fit", "a_minus_CI", "a_plus_CI", "frac_a_minus_less_than_a_plus", "log10_ratio_plus_over_minus_CI"):
        add(S, f"{g} RPE 2a {q}", f, [g, q])
add(S, "Active minus Passive alpha- CI", f, ["Active_minus_Passive_a_minus_CI"])

# ---------------------------------------------------------------- trial by trial, observed (m13)
S = "5. Trial by trial: observed statistics, drops, timing (m13; v1 file, observed values are valid)"
f = "m13_summary.json"
for g in GROUPS:
    add(S, f"{g} day-14 drops: mean, min, max", f, [g, "day14_drops"], "frac_by_2p5min = share earned by 2.5, 5, ... 15 min")
    for k in ("day14", "day15", "day16", "d14_h1", "d14_h2", "d15_h1", "within14_halves", "night_halves",
              "within14_5min", "night_5min", "d15_minus_d14"):
        r = get(f, [g, "A_observed", k])
        ROWS.append((S, f"{g} observed {k}: mean [lo, hi], p, n down/up",
                     [r["mean"], r["lo"], r["hi"], r["p"], r["n_neg"], r["n_pos"]], f, [g, "A_observed", k], ""))
    add(S, f"{g} group-curve sigmoid, morphine -> water", f, [g, "sig_obs", "groupcurve", "deg"])
    add(S, f"{g} group-curve sigmoid, water -> morphine", f, [g, "sig_obs", "groupcurve", "enh"])
    add(S, f"{g} bootstrap over mice (2000): CIs and share later", f, [g, "groupcurve_boot"])
    add(S, f"{g} per-mouse midpoints, morphine -> water", f, [g, "sig_obs", "deg", "mid"])
    add(S, f"{g} per-mouse midpoints, water -> morphine", f, [g, "sig_obs", "enh", "mid"])
    add(S, f"{g} per-mouse midpoint test (paired)", f, [g, "sig_tests", "deg_vs_enh_mid_all"])
    for m in MAIN:
        add(S, f"{g} {NICE[m]} noise-free group-curve midpoints/starts", f, [g, "A_models", m, "sig_pred_groupcurve"])
add(S, "Active vs Passive Welch p (observed statistics)", f, ["Active_vs_Passive"])

# ---------------------------------------------------------------- trial by trial v2 (m15)
S = "6. Trial by trial v2: model checks, placebo, calibration, timing (m15; CURRENT)"
f = "m15_summary.json"
for g in GROUPS:
    add(S, f"{g} share of drops earned by 2.5 min, by day", f, [g, "drops_frac_by_2p5min"])
    add(S, f"{g} P(lick) in the first 2 min (bins 0-3), by day", f, [g, "p_lick_bins0_3"])
    add(S, f"{g} equal 60-trial windows", f, [g, "equal_window"])
    add(S, f"{g} leave one mouse out (group curve)", f, [g, "leave_one_out"])
    for m in ALLM:
        base = [g, "models", m]
        a = get(f, base + ["A", "day14"])
        ROWS.append((S, f"{g} {NICE[m]} day 14: obs, pred, mice-model, mouse p, Holm p",
                     [a["obs"], a["pred"], a["mouse_t"]["mean_diff"], a["mouse_t"]["p_t"], a["mouse_t"].get("p_holm")],
                     f, base + ["A", "day14", "mouse_t"], "Holm over the 4 main models only"))
        ROWS.append((S, f"{g} {NICE[m]} day 14: noise-matched p A, B; naive p (WITHDRAWN)",
                     [a["calibrated"]["p"], a["calibrated_B"]["p"], a["naive"]["p"]],
                     f, base + ["A", "day14"], "naive = independent bins; never cite"))
        h = get(f, base + ["A", "d14_h1"])
        ROWS.append((S, f"{g} {NICE[m]} first half of day 14: pred, mouse p, A, B",
                     [h["pred"], h["mouse_t"]["p_t"], h["calibrated"]["p"], h["calibrated_B"]["p"]],
                     f, base + ["A", "d14_h1"], "uncorrected; descriptive"))
        pl = get(f, base + ["placebo", "whole"])
        ROWS.append((S, f"{g} {NICE[m]} placebo (day 13 - day 12): obs, pred, mouse p, A, B, naive",
                     [pl["obs"], pl["pred"], pl["mouse_t"]["p_t"], pl["calibrated"]["p"], pl["calibrated_B"]["p"],
                      pl["naive"]["p"]], f, base + ["placebo", "whole"], ""))
        cal = get(f, base + ["calibration"])
        ROWS.append((S, f"{g} {NICE[m]} calibration: T1/T2/T3 observed; A (sig_u, sig_a, phi, edge); B (...)",
                     [cal["T1_obs"], cal["T2_obs"], cal["T3_obs"],
                      [cal["A"]["sig_u"], cal["A"]["sig_a"], cal["A"]["phi"], cal["A"]["at_grid_edge"]],
                      [cal["B"]["sig_u"], cal["B"]["sig_a"], cal["B"]["phi"], cal["B"]["at_grid_edge"]]],
                     f, base + ["calibration"], ""))
        if m in MAIN:
            for key in ("g_deg_start", "g_deg_mid", "g_enh_mid", "g_deg60_mid"):
                t = get(f, base + ["timing", key])
                ROWS.append((S, f"{g} {NICE[m]} timing {key}: obs, p A, p B, sim median A, 95% range A",
                             [t["obs"], t["calibrated"]["p"], t["calibrated_B"]["p"], t["cal_median"],
                              [t["calibrated"]["lo"], t["calibrated"]["hi"]]], f, base + ["timing", key], ""))

# ---------------------------------------------------------------- write
md = ["# 11. Numbers index (auto-generated; do not edit by hand)", "",
      "Generated by `tools/build_numbers_index.py` from `project/data/*.json`. Values are rounded to 4 decimals "
      "(very small p values in scientific notation). 'Path' is the key path inside the JSON file. Group order of "
      "per-mouse lists: see docs/02 section 2.6. Model names: instant switch = `today's outcome`, RPE 2a = `RPE2a`, "
      "meta-RPE = `mRPE`.", ""]
sec = None
for s, label, v, fname, path, note in ROWS:
    if s != sec:
        md += ["", f"## {s}", "", "| Quantity | Value | File | Path | Note |", "|---|---|---|---|---|"]
        sec = s
    p = " → ".join(str(x) for x in path)
    md.append(f"| {label} | {fmt(v)} | `{fname}` | `{p}` | {note} |")
out_md = ROOT / "docs" / "11_numbers_index.md"
out_md.write_text("\n".join(md) + "\n", encoding="utf-8")
out_js = ROOT / "docs" / "numbers_index.json"
json.dump([dict(section=s, quantity=label, value=v, file=fname, path=path, note=note)
           for s, label, v, fname, path, note in ROWS], open(out_js, "w", encoding="utf-8"), indent=1)
print(f"wrote {out_md.name} ({len(ROWS)} rows) and {out_js.name}")
