# 8. Code map and how to reproduce

The scripts in `project/` are **byte-identical copies** of the originals; `CHECKSUMS.sha256` lets you verify this.
Each script reads from and writes to `data/` next to itself (`D = HERE / "data"`), so run them from inside
`project/`. Run them in numeric order.

## 8.1 Environment used

| | |
|---|---|
| Machine | Windows 11 |
| Python | 3.13.3 |
| Libraries | numpy 2.2.5, scipy 1.17.1, pandas 2.3.3, numba 0.66.0, matplotlib 3.10.8, python-pptx 1.0.2, statsmodels 0.14.6 (m17, m18 and the docs checker only) |
| Previews (slides only) | LibreOffice headless (`C:\Program Files\LibreOffice\program\soffice.exe`), run with its own profile |
| Fonts (slides only) | Arial; Malgun Gothic for Korean |

The fitting scripts use process pools with up to 30 workers.

## 8.2 Script by script

| Script | Does | Main inputs | Main outputs (in `project/data/` unless stated) | Seeds |
|---|---|---|---|---|
| `m01_parse_trials.py` | Turns the raw controller logs into one row per PR trial and every lick and drop | raw `*.jsonl` on `K:` (not included); session list and `pairs_from_raw.csv` | `trials_pr.csv`, `licks_all.csv`, `sessions.csv`, `m01_check.json` | - |
| `m02_hjort_models.py` | Module: the paper's six models ported line by line (numba); `fit()` with the notebook's bounds and start grid | - | - | - |
| `m03_validate_on_paper_data.py` | Refits the authors' 40 mice (H→L cue, 160 trials) | `hjort2026_sample_choice/` | `m03_paper_data_fits.csv`, `m03_paper_data_summary.json` | - |
| `m04_fit_user_data.py` | Builds the event stream per mouse (`build_events`: drops and 30-s bins in time order) and fits 8 models per mouse | `licks_all.csv`, `pairs_from_raw.csv` | `m04_events.pkl`, `m04_fits.csv`, `m04_traces.pkl`, `m04_summary.json` | deterministic |
| `m05_recovery_and_group_fits.py` | (b) group fits with shared learning rates (profile likelihood); bootstrap over mice; (a) per-mouse model recovery | m04 outputs | `m05_group_fits.json`, `m05_recovery.json`, `m05_recovery_sets.csv`, `m05_log.txt` | bootstrap 0-999; recovery 1000 + mouse index |
| `m06_modelfree_switch.py` | p(lick) per 3-min block, switch day minus the day before (the first-9-min window is superseded) | `m04_events.pkl` | `m06_switch_blocks.csv`, `m06_modelfree.json` | - |
| `m07_fig1_paper_explained.py` | Redraws the paper's Fig. 1 from the authors' data | authors' data | `figs/fig1_paper_fig1_explained.png`, `m07_fig1_numbers.json` | - |
| `m08_bootstrap_asymmetry.py` | RPE 2α group fit, 1,000 bootstraps over mice | m05 | `m08_bootstrap_asymmetry.json` | 5000 + s |
| `m09_figures_user_data.py` | Figures 2-4 of the first round | m04-m08 | `figs/fig2..fig4*.png` | - |
| `m10_make_slides.py` | Slides 1-6 v1 (EN, KR) and the slide helpers (measured text boxes, KEY boxes, number registry, checks, LibreOffice render) | m03-m08, lab template | `slides_20261117/` (v1) | - |
| `m11_extra_checks.py` | Whole-session model-free tests, Welch t, group-level recovery, LOMO-CV; `--asym` adds the asymmetry test | m04, m05, m06 | `m11_modelfree.json`, `m11_group_recovery.json`, `m11_lomo_cv.json` and their CSVs | 7000 + 100·i + s |
| `m12_make_slides_v2.py` | Slides 1-6 v2 (CURRENT) | as m10 + m11 | `slides_20261117/*_v2_*` | - |
| `m13_trial_by_trial.py` | Trial-by-trial v1: statistics, sigmoids, bootstrap, naive check (withdrawn) | m04, m05 | `m13_bins.csv`, `m13_per_mouse.csv`, `m13_summary.json`, `m13_log.txt` | 20261118 |
| `m14_make_trial_slides.py` | Slides 7-8 v1 and v3 decks (superseded) | m13 | `slides_20261118_trial_by_trial/` (v1, v3) | - |
| `m15_trial_checks_v2.py` | Noise-matched check (calibrations A and B), placebo, mouse-level t, equal windows, leave one out | m13 functions + m04/m05 | `m15_summary.json`, `m15_log.txt` | 20261119 |
| `m16_make_trial_slides_v2.py` | Slides 7-8 v2 and v4 decks (CURRENT); 134 checks | m13, m15 | `slides_20261118_trial_by_trial/*_v2_*`, `*_v4_*` | - |
| `m17_make_summary_2slides.py` | (26 Nov) 2-slide summary deck (EN, KR) and the sample-size calculation (statsmodels: two-sample and paired t; Welch-t simulation, 20,000 draws per n); 100 number checks + 20 text-fit checks | m03, m05, m08, m11, m13, m15 outputs; m10 helpers | `slides_20261126_summary/` (copied to `slides/summary_20261126/`) | simulation `default_rng(1)` |
| `m18_planned_tests.py` | (26 Nov) the pre-specified tests of docs/12 (P1, S1-S4, S2 sensitivity, Holm, pooled estimates); on the Dec-2025 cohort it reproduces the m11 changes, CIs, p values and the asymmetry test exactly (asserted); `--selftest` checks the pooled code | `licks_all.csv`, `sessions.csv`, `m11_modelfree.json` | `m18_planned_tests_dec2025.json`; for a new cohort, `m18_planned_tests.json` in its data folder | - |

### Key functions to know

**`m04.run(code, lp, ro, ...)`** — the event loop.

- A drop updates the value according to the learning rule.
- A bin adds `-log P(choice)` to the likelihood when `use == 1`.
- It optionally returns traces: Q, δ, CD, CE, M and P(lick).

**`m05.q_at_bins(code, lp, typ, r, dayout)`** — the value seen by each bin. `m05.logit_nll(X, y)` is an exact
logistic-regression readout (Newton steps with a 1e-6 ridge).

**`m05.fit_group(name, mice)`** — group fit by profile likelihood.

**`m13`:**

- `build()` returns per-mouse arrays (8 days × 30 bins) of observed licks (`OBS`), model predictions (`PRED`),
  drops per bin (`DROPS`) and reference days.
- `a_stats()` computes the day-14 statistics.
- `traj_deg()` / `traj_enh()` build the trajectories from each switch.
- `Grid.fit()` is the exact grid least-squares sigmoid.
- `paired()`, `ttest1()`, `signflip()`, `holm()` are the tests.

**`m15`:**

- `simulate(P, su, sa, draws, phi)` is the noise model.
- `disp_targets()` computes T1-T3.
- `calibrate()` returns calibrations A and B.
- `placebo()` computes day 13 minus day 12.
- `timing_stats()` returns per-mouse and group-curve timing.

### Results-file structure (most-used keys)

**`m13_summary.json` → `<group>`:**

- `day14_drops`: totals and the share of drops by each 2.5 min.
- `A_observed.<key>`: observed statistic with `mean, lo, hi, p, n_neg, n_pos, signflip_p`.
- `A_models.<model>.<key>`: model check; contains the **naive** `ppc_p` (withdrawn).
- `sig_obs`: per-mouse and group-curve sigmoid fits.
- `groupcurve_boot`: bootstrap CIs.
- `sig_tests`: per-mouse paired tests.
- `C_start_vs_models`: start of the change versus each model.

**`m15_summary.json` → `<group>`:**

- `obs_timing`, `equal_window`, `leave_one_out`, `drops_frac_by_2p5min`, `p_lick_bins0_3`.
- `models.<model>`:
  - `calibration.{A,B}`, including `sig_u, sig_a, phi, T*_sim` and `at_grid_edge`;
  - `A.<key>.{obs, pred, calibrated, calibrated_B, naive, mouse_t}`;
  - `placebo.<whole|h1|h2>`;
  - `timing.<key>`.

## 8.3 How to rerun from this package

1. **`m01` cannot be rerun.** It needs the raw logs on `K:`. Its outputs (`licks_all.csv`, `trials_pr.csv`,
   `sessions.csv`) are included.
2. **Hard-coded Windows paths** to edit on another machine:
   - `m04_fit_user_data.py` line 255 reads `pairs_from_raw.csv` from
     `C:\Users\hsollim\Research_Projects\01_Opioid_Behavior\yoked_stats_20261116\data`. Point it to
     `project/data/pairs_from_raw.csv` (included).
   - `m01_parse_trials.py` line 23 has the same folder (only needed with raw data).
   - `m10_make_slides.py` lines 48-50 give the lab slide template (not included) and LibreOffice; lines 554-555 give
     the font files. The slide scripts (m10, m12, m14, m16, m17) are not needed to understand or check the analysis.
   - `m17_make_summary_2slides.py` line 52 points to a power-results file outside the package (not included). It is
     used only for two extra cross-checks of the power numbers, which are skipped when the file is absent (line 124),
     so a rerun elsewhere reports 98 checks.
3. **Pickles are not shipped.** Run `m04` first; it rebuilds `m04_events.pkl` and `m04_traces.pkl` deterministically
   from `licks_all.csv`.
4. Then run m05 → m06 → m08 → m11 (and `m11 --asym`) → m13 → m15. m03 needs only the authors' data; m07 needs
   the authors' data and the m03 outputs. `m18` needs only `licks_all.csv`, `sessions.csv` and `m11_modelfree.json`,
   has no hard-coded paths, and runs directly from `project/` (`python m18_planned_tests.py`).
5. **Work on a copy.** Every script writes into `project/data/` and overwrites the shipped results.
6. **Expected result:** the JSON files match the shipped ones, up to floating-point differences between machines.
   The seeds are fixed (see the table).

## 8.4 Tools in this package

- **`tools/build_numbers_index.py`** regenerates `docs/11_numbers_index.md` and `docs/numbers_index.json` from the
  JSON files, with each number's source path.
- **`tools/check_guide_numbers.py`** checks that the key numbers quoted in `README.md` and `docs/*.md` equal the
  values in the JSON files. Run `python tools/check_guide_numbers.py` from the package root; it should print
  `ALL OK`. Since 2026-11-26 it also recomputes the sample-size numbers (README item 13 and the last section of
  docs/09) with statsmodels.
