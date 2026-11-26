# Meta-RPE model (Hjort et al. 2026) applied to the opioid Active/Passive licking cohort

**Teaching package for another analyst or AI agent.** Version 2026-11-19. Status: CURRENT.
The trial-by-trial results are version 2. Version 1 is kept in the code but its claims are withdrawn.
**Added on 2026-11-26:** a 2-slide summary deck (English, Korean) with a sample-size calculation, in
`slides/summary_20261126/` (docs/09, last section; item 13 below). Earlier results and decks are unchanged.

This package explains, from first principles, a re-analysis of our oral-morphine licking experiment with the learning
models of a 2026 Nature paper. It covers what was done, every result with its number, and what must **not** be
claimed. Everything here is a copy. The original project lives on the lab PC at
`C:\Users\hsollim\Research_Projects\01_Opioid_Behavior\mrpe_hjort2026_20261117\`.

> If you are an AI agent: read `docs/07_rules_for_agents.md` before you write anything about this work.
> Take numbers from `docs/11_numbers_index.md`, which is generated from the result files. Never take them from
> memory or from superseded (v1) text.

---

## Reading order

| # | File | What it teaches |
|---|---|---|
| 1 | [docs/01_paper_and_model.md](docs/01_paper_and_model.md) | The reference paper, its task, and the six learning models with their equations |
| 2 | [docs/02_experiment_and_data.md](docs/02_experiment_and_data.md) | Our cohort, the 18-day design (phase names), yoke pairs, and the data tables |
| 3 | [docs/03_mapping.md](docs/03_mapping.md) | How the paper's trial, choice, outcome and value were mapped onto free-operant PR licking |
| 4 | [docs/04_model_fits.md](docs/04_model_fits.md) | Model-free switch tests, per-mouse and group fits, cross-validation, recovery, learning-rate asymmetry |
| 5 | [docs/05_trial_by_trial.md](docs/05_trial_by_trial.md) | The paper-style trial-by-trial analyses (Fig. 1d,e,m analogues), the noise-matched check and the placebo |
| 6 | [docs/06_conclusions_and_limits.md](docs/06_conclusions_and_limits.md) | What holds, what does not, what cannot be tested, alternative explanations, next designs |
| 7 | [docs/07_rules_for_agents.md](docs/07_rules_for_agents.md) | Binding wording rules and known traps |
| 8 | [docs/08_code_and_reproduce.md](docs/08_code_and_reproduce.md) | Script-by-script map (m01 to m17), inputs and outputs, how to rerun |
| 9 | [docs/09_slides.md](docs/09_slides.md) | The current 8-slide deck, slide by slide; the 2-slide summary of 2026-11-26 at the end |
| 10 | [docs/10_glossary.md](docs/10_glossary.md) | Every term and symbol |
| 11 | [docs/11_numbers_index.md](docs/11_numbers_index.md) | Every key number with its source file and JSON key (auto-generated) |

---

## The answer in brief

Unless noted, numbers are for Active mice (n = 6) and p values are two-sided.

1. **The port is faithful.** We refitted the authors' 40 mice with their own code (ported to Python). This
   reproduces their Fig. 1k,l: meta-RPE fits best, with a BIC change from best of 2.6 versus 5.7 to 18.6 for the
   other models.
2. **The task had to be mapped.** In our task:
   - one trial is a 30-s bin;
   - the choice is whether the mouse licked at all in that bin;
   - the outcome is the identity of each drop (morphine = 1, water = 0);
   - the value is the expected drug from the spout.

   Each mouse has exactly one switch in each direction: morphine to water on day 14, and water to morphine on day 17.
3. **Passive mice are not a reward-omission case.** The readout weight on value is β > 0 in all 6 Active mice and
   β < 0 in all 7 Passive mice. Passive mice lick less when they expect morphine. So the paper's "H to L"
   (degradation) analogy applies to **Active mice only**.
4. **Model-free:**
   - On the first water day, Active licking did not change detectably: -0.06 (95% CI -0.27 to +0.15, p = 0.53).
     Passive licking rose by +0.40.
   - Active licking fell on days 15 and 16 (-0.27, -0.35) and rose on the first morphine day (+0.35).
5. **Group fits** use learning rates shared within a group and a readout fitted per mouse.
   - In Active mice the best model is RPE with two learning rates: alpha- = 0.01 (its lower bound) and
     alpha+ = 0.12, with alpha- < alpha+ in 99.5% of bootstraps.
   - Held-out mice agree in total. Relative to the best model, the held-out negative log-likelihood is +1.3 for
     meta-RPE, +8.3 for one-rate RPE and +15.5 for instant switch. RPE 2α and meta-RPE are effectively tied
     (meta-RPE is better in 4 of 6 held-out mice); the clear part is that both beat one-rate RPE and instant switch.
   - In Passive mice, every learning rule is about as good as an instant switch.
6. **Meta-RPE's defining part (error averaging) cannot be tested with this design.** When we simulated data from
   meta-RPE, it was recovered as meta-RPE in only 2% of per-mouse simulations and 16% of group simulations.
7. **Trial by trial, day 14.**
   - In the first 7.5 min there was no clear fall in licking (+0.02), although 88% of the day's water drops had
     been earned by then.
   - The whole-day change was -0.06. The models predicted -0.11 (meta-RPE) to -0.22 (instant switch).
   - With n = 6, no model is robustly rejected. Mouse-level paired t: instant switch p = 0.073, RPE 0.15,
     RPE 2α 0.30, meta-RPE 0.58.
8. **Within day 14 or overnight cannot be separated.** Within day 14 the change was -0.16 (95% CI -0.60 to +0.29)
   and overnight -0.02 (95% CI -0.53 to +0.49).
9. **Timing.**
   - On the group-mean curve, the midpoint of the morphine-to-water change is trial 36.5 and of the water-to-morphine
     change trial 4.5. The first is later in more than 99% of bootstraps.
   - The paper's per-mouse test gives p = 0.11 (n = 6), which is not significant.
   - The start of the change (trial 14.4, 95% CI 0 to 47) excludes no model.
10. **Statistical lesson.** A simulation check that draws 30-s bins independently is far too narrow: it rejected a
    no-change placebo day pair. Use mouse-level tests, and validate any simulation check on a placebo first.
11. **Caveat added on 2026-11-19 while writing this guide (not on the slides).** In Passive mice, even the
    placebo pair (day 13 minus day 12: -0.19) differs from every model at the mouse level (p = 0.046 to 0.053).
    Passive day 13 was unusually low. So "Passive rose more than every model predicted on day 14" depends on the
    reference day. Measured from day 12, the Passive rise is +0.22, smaller than the +0.30 to +0.32 the main models
    predict.
12. **Bottom line.** The paper's behavioural signature is present descriptively. Active mice did not change while
    the errors were largest, and changed later; the morphine-to-water change was slower than the water-to-morphine
    change on the group curve. But with 6 mice and one switch per direction, no learning model is rejected and error
    averaging is untestable. A delayed oral-drug effect, a withdrawal state, or taste could produce the same lag.
13. **Sample size for a follow-up (added 2026-11-26; `m17`, statsmodels).** The effect is the Active day-14 lag:
    the day-14 change (-0.06) versus the mean change on days 15-16 (-0.31), divided by the between-mouse SD of the
    change (0.21, the mean of four SDs recovered from the 95% CIs in `m11_modelfree.json`): d = 1.18.
    - Two groups, two-sided t test at alpha = 0.05: 80% power needs 12.3 mice per group, so 13 (90% power: 16.1).
      The current n = 6 gives 46% power.
    - Each mouse tested in both conditions (paired t): 80% power needs 7.7 mice, so 8, if the two measures correlate
      at r = 0.5 (10.0 at r = 0.3, 5.6 at r = 0.7).
    - Half the effect, or a group x condition interaction: 46 mice per group.
    - d is estimated from these 6 mice and is uncertain; a smaller true effect needs more mice.

## Ten rules (full list in docs/07)

1. Write that drops were **earned** (or delivered, on replay days), not "drunk". Drinking was not measured.
2. Phase names are Baseline (days 1-5), Yoked (6-10), **Equal Access** (11-13), Withdrawal (14-16) and Re-exposure
   (17-18). The code's `phase` column still uses old labels: Train/Pre = Baseline, During = Yoked,
   **Post = Equal Access**. Never write "Post" in new material.
3. Never cite p values from the independent-bin ("naive") simulation check in `m13` as evidence.
4. Always say whether a timing result is for the **group curve** or **per mouse**. The per-mouse test is p = 0.11.
5. BIC treats 30-s bins as independent and overstates the evidence. Prefer held-out mice and mouse-level tests.
6. Never claim that meta-RPE or error averaging was shown or rejected. It is untestable with one switch per
   direction.
7. Never write "learned" without the pharmacology caveat: oral morphine acts minutes later, which is a delayed
   outcome. Withdrawal state and taste are alternatives.
8. The H-to-L degradation analogy holds for Active mice only, because Passive β < 0.
9. Group sizes are n = 6 Active and 7 Passive in every comparison. Mouse 6099_white (Passive) died on day 13; there
   are 8 Passive mice in total.
10. Three version-1 statements are withdrawn: "no extra change overnight", "the start rules out instant switch", and
    "the simulation p values are conservative".

---

## Folder map

```
README.md               this file
docs/                   the guide (01-11)
project/                byte-identical copy of the analysis scripts m01-m17 + the project log (README.md)
  project/data/         derived tables and JSON results (no raw data, no pickles)
  project/figs/         figures 1-4 of the first analysis round
slides/                 CURRENT 8-slide decks (EN, KR), slide previews (PNG), figure PNGs, slide check reports
  slides/summary_20261126/  2-slide summary deck (EN, KR) of 2026-11-26: previews (PNG, PDF), figures, check report
release/                zip of the package as of 2026-11-19 (made before the summary deck was added)
tools/                  build_numbers_index.py (writes docs/11) and check_guide_numbers.py (checks the docs)
CHECKSUMS.sha256        sha256 of every file in the package
```

Raw controller logs (`.jsonl`, one per mouse-day, on the lab drive `K:`, read only) are not included.
`project/data/licks_all.csv` holds every lick and drop time and is enough to rerun everything from `m04` onward. See
docs/08.

## Reference

Hjort, ..., Stuber (2026). Prefrontal to ventral tegmental area dynamics drive contingency degradation.
*Nature* 655, 174. doi:10.1038/s41586-026-10443-5. Code and sample data: github.com/stuberlab/Hjort-et-al.-2026-PFC-and-reversal-learning
(commit e06ae22, MIT licence). The authors' 40-mouse sample data and licence are in `project/data/hjort2026_sample_choice/`.
