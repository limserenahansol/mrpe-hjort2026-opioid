# 7. Rules for agents (binding)

Read this file before writing any text, slide, figure label, e-mail or code comment about this work. Each rule has a
reason, so it can be applied to cases not listed.

## A. Wording

1. **"Earned" (PR days) or "delivered" (replay days), never "drank" or "drunk".**
   *Why:* consumption was not measured. A drop is a delivery at the spout.
2. **Phase names:** Baseline (days 1-5), Yoked (6-10), Equal Access (11-13), Withdrawal (14-16), Re-exposure (17-18).
   **Never write "Post", "POST" or "Post-Yoked"** in new material. Write "Equal Access (formerly Post)" once if the
   reader may know the old name.
   *Why:* decided on 2026-11-17 for all new material. The code's `phase` column still holds the old labels
   (`Train`/`Pre` = Baseline, `During` = Yoked, `Post` = Equal Access). Do not copy them into text.
3. **"Instant switch"** is the slide name of the `today's outcome` model. Use one name per document and define it
   once.
4. **Say "group curve" or "per mouse"** for every timing result.
   *Why:* the group-curve midpoint difference (36.5 vs 4.5) is robust, but the per-mouse test is not significant
   (p = 0.11).
5. **Do not write "learned", "learning rate shows" or "the mice expected"** without the caveat that oral morphine
   acts minutes later (a delayed outcome), and that withdrawal state and taste are alternatives.
6. **"Model fits best"**, not "the mice use this algorithm". The titles state model fits, not mechanisms.
7. **Passive is not reward omission.** Do not describe the Passive morphine → water switch as contingency degradation.
   *Why:* β < 0 in all 7 Passive mice.

## B. Statistics

8. **Never cite the naive (independent-bin) simulation p values** from `m13` (keys `ppc_p`, `A_models.*.ppc_p`,
   `C_start_vs_models.*.ppc_*`) or the `naive` entries in `m15` as evidence.
   *Why:* they are anti-conservative and reject a no-change placebo.
9. **Primary test = mouse-level paired t** (n = 6 / 7). Report the noise-matched check as a range "A (B)".
   *Why:* A and B bracket the true noise; neither matches all three features.
10. **Any new simulation-based test must first pass the placebo** (day 13 minus day 12, the same pipeline).
11. **BIC overstates the evidence** because 30-s bins are not independent. Pair every ΔBIC with the held-out-mouse
    result (LOMO-CV) and with recovery.
12. **Per-mouse model comparison has no power here.** Meta-RPE is recovered in 2% of its own simulations. Never
    report "meta-RPE lost in most mice" as evidence against it.
13. **Do not claim within-day versus overnight.** It is untestable: the CIs are about 1 unit wide. v1's "no extra
    change overnight" is withdrawn.
14. **Do not claim the start of the change rules out any model** (p ≥ 0.068). v1's statement is withdrawn.
15. **Passive model tests are unreliable**, because the Passive placebo fails at the mouse level (p = 0.046). Do
    not report Passive "the models under-predict" as a finding.
16. **n = 6 Active and 7 Passive** in every comparison. 6099_white (Passive, pair 6) died on day 13 and is left out;
    there are 8 Passive mice in total.
17. **Bounds.** α- = 0.01, d = 0.01 and e = 0.10 sit at their bounds. Say so whenever you quote them.

## C. Numbers

18. Take numbers from `docs/11_numbers_index.md` or the JSON paths given there. Changes are given to 2 or 3 decimals and p values to
    2 significant digits; keep the precision of the source you quote.
19. **"88% of drops by 7.5 min"** refers to Active mice on day 14 (28.3 drops on average). Drops are front-loaded on
    every PR day (day 13: 58% by 2.5 min). Do not present front-loading as special to day 14.
20. **Do not mix versions.** Trial by trial: v1 = `m13` + `m14` (`*_2slides_*_20261118.pptx`, `*_8slides_*_v3_*`);
    v2 = `m15` + `m16` (`*_2slides_*_v2_20261118.pptx`). Slides 1-6: v1 = `m10`, v2 = `m12`
    (`*_6slides_*_v2_20261117.pptx`). The CURRENT deck is the 8-slide v4.

## D. Files and code (lab conventions)

21. **Only add, never delete or overwrite** earlier versions. Rerunning a script overwrites `project/data/`, so
    work on a copy of the package. A new analysis, check or figure goes *next to* the old
    one, with a new version suffix.
22. **Raw data are read only.** Never write to `K:` or `H:`.
23. **Slides:** build with python-pptx only. Never drive PowerPoint through Office COM: it closes the user's open
    documents. Render previews with LibreOffice headless, using its own profile.
24. **Slide fixes:** send only the edited slide (a 1-slide pptx + PNG), not a rebuilt deck.
25. **Every graph** has a title, axis labels, a legend, n and the statistical test. Each slide has its title on top
    and a 1-2 sentence "Did / Found" at the bottom.
26. **Do not invent** a value that is not in the files, such as a morphine concentration, drinking volumes, or
    session-to-session yoking at the same moment.

## E. Known traps (each of these happened once)

| Trap | Correct statement |
|---|---|
| "The simulation check is conservative because it is not refitted" | It is anti-conservative: bins are not independent |
| "Instant switch rejected, p < 0.001" | Mouse-level p = 0.073; noise-matched 0.021-0.050 |
| "No extra change overnight" | Within day 14 vs overnight cannot be separated |
| "The start rules out instant switch" | The start excludes no model |
| "Slide 8a is normalised" | Raw change in P(lick) |
| "Passive rose more than every model predicted" | Depends on a low day-13 reference; placebo fails |
| "CD is the slow average, CE the fast one" (as a model rule) | That is only our fit (d = 0.01, e = 0.10 at bounds); in the authors' mice d < e in only 8 of 40 |
| "r2 = R²" in `m13` | It is uncentred (1 - SSE/Σy²) |
