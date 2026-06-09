# Verification checklist

After the first-pass model file is written, re-read the source independently and walk through every item below. Each item has a failure mode that has caused real translation errors in the past. Fix everything you can; **flag anything ambiguous to the user** using this format:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

Do not silently resolve ambiguity. Do not tune parameters to make a validation output match a target — if the validation disagrees with the paper, investigate the source, not the parameters.

## H. Source identity

- [ ] On-disk PDF / XML title, first author, journal, and year match the task's `Paper metadata` block. Filename PMID matches the actual PMID in the source (catches mislabelled drops like `PMID_23436260.pdf` actually being Frey 2010 / PMID 20097931).
- [ ] Drug named in the task metadata matches the molecule the paper actually models.
- [ ] **Species recorded.** `population$species` is set to the species the final model was fit to (e.g., `"human"`, `"rat (Sprague-Dawley)"`, `"beagle dog"`, `"in vitro (SKBR3 cell line)"`, `"human + rat"` for pooled). Non-human models additionally prepend the species to `description` (e.g., `"Preclinical (rat). ..."`). Filename carries a species suffix when the same drug has both human and animal extractions (e.g., `Geldof_2008_fluvoxamine_rat.R`). Preclinical and in-vitro models are first-class and are extracted without sidecar-asking — the species check is a labelling step, not a gating decision. See SKILL.md Phase 1 step 3.
- [ ] **Model class recorded** (PBPK / QSP / MBMA / mechanistic-systems papers). `description` prefixed with the framework label (e.g. `"PBPK (whole-body, SimCYP V12.1). ..."`, `"QSP. 33 ODEs + 113 parameters ..."`, `"MBMA. Between-study variance only ..."`). Filename uses a `_pbpk` / `_qsp` / `_mbma` suffix when a popPK extraction of the same paper / drug already exists. See `references/pbpk-qsp-mbma.md`.
- [ ] **No class-typical-default substitutions.** For PBPK / QSP / MBMA models, every parameter in `ini()` is sourced from the paper text, an on-disk supplement, or an on-disk upstream paper that the current paper explicitly cites. NO parameter was filled in from training-data knowledge of "what a similar PBPK paper for this drug typically uses", default SimCYP / GastroPlus / OSP library entries, or class-typical literature ranges. Missing parameters went through the Phase 4 missing-parameter sidecar, not silent substitution. See `references/pbpk-qsp-mbma.md` — this rule is stricter for mechanism models than for popPK.

## A. Parameter values

- [ ] **Errata / corrigenda checked.** Confirmed search for published corrections to the source; any erratum value supersedes the main publication for conflicting values, with the most recent erratum winning when multiple exist. Erratum citation recorded in the model file's `reference` field, and each affected `ini()` comment points to the erratum.
- [ ] Every parameter line in `ini()` has a clear provenance comment if it did not come from the paper's text/tables (author correspondence, figure-digitisation, upstream-task model file). The comment cites the followup-register entry (`tracking/operator_followups.md F<n>`) when applicable.
- [ ] Every parameter in `ini()` has an in-file trailing comment pointing to the source location (table, equation, page, or figure). Re-check each comment matches what the source actually says.
- [ ] Values are **final estimates**, not initial estimates. Supplement NONMEM control streams often list initial values in `$THETA` and `$OMEGA`; the final values come from the `$TABLE` output, an accompanying `.lst` `FINAL PARAMETER ESTIMATE` block, or the main paper. If the only source is a control stream, confirm the values match any published point estimates.
- [ ] **Log-vs-linear reporting.** NONMEM often reports THETAs on the estimation scale (already log), but tables in the paper usually show the back-transformed value. A `log()` wrapper in `ini()` must match what the paper reports: `lcl <- log(0.0388)` is correct when the paper says "CL = 0.0388 L/day."
- [ ] **CV% vs. variance.** `omega²` in NONMEM output is the variance on the internal scale. For log-normal parameters, CV% relates via `omega² = log(CV² + 1)`. Do not paste CV% directly into `ini()` as if it were a variance.
- [ ] **Correlated IIV.** If the paper reports a correlation `r` and individual CV%, the covariance is `cov = r × sqrt(var_1 × var_2)`. Verify the block matrix entries match this formula.
- [ ] **Fixed parameters** the source holds fixed are wrapped in `fixed(...)` in `ini()` — applies to ALL parameter types (THETAs, allometric exponents, IIVs, residual errors, covariate effects, bioavailability anchors), not just IIVs. Source signals: explicit "fixed at <value>" prose, NONMEM `FIX` flags on `$THETA`/`$OMEGA`/`$SIGMA`, allometric exponents reported without uncertainty, bioavailability `F1=1` set as structural anchor, parameters inherited from upstream papers without re-fitting. If a parameter is reported without uncertainty but the paper does not explicitly say "fixed", sidecar-ask before guessing — see `references/parameter-names.md` § "Fixed parameters" for the encoding examples.

## B. Structural model

- [ ] The number of compartments matches the source.
- [ ] ODEs (or `linCmt()` shortcut) reproduce the equations shown in the paper.
- [ ] **Reference weight / age** for allometric and maturation terms matches the source (70 kg adult? 75 kg? 5 kg infant? 40 weeks PMA? term birth?).
- [ ] **Allometric exponents** apply to the right parameters (usually 0.75 for CL/Q and 1.0 for volumes, but papers often fit custom exponents — use whatever the paper reports).
- [ ] **Maturation form** matches the paper (sigmoidal asymptotic vs. Hill-type vs. exponential). The functional form determines whether `beta_cl` is "fraction mature at birth" or something else.
- [ ] **Bioavailability** (`f(depot) <- ...`) applied to the correct compartment. `F1` in NONMEM sometimes targets a different compartment than you'd expect.
- [ ] **Lag time / Tlag**, transit compartments, and zero-order absorption phases match the source.
- [ ] **Units consistency.** Dose units × bioavailability ÷ volume units must yield the concentration units declared in `units`. Walk through the dimensional analysis at least once.
- [ ] **Full dimensional analysis (mandatory for endogenous / mechanistic models).** For every ODE term and every derived rate / flux, write down units of each symbol and multiply them out. The right-hand side of `d/dt(state)` must equal `[state]/[time]`. Past endogenous-model bugs caught only here: `igg_kim_2006` V1 mislabeled `(mg/kg)` instead of `(mL/kg)`; `phenylalanine_charbonneau_2021` `daily_phe_intake` augmentation line carried a stray `vd` factor, reporting `(L/kg)·mg/day` instead of `mg/day`. See `references/endogenous-validation.md`.

## C. Covariate effects

- [ ] Every covariate used in `model()` is registered in `inst/references/covariate-columns.md` with a canonical name, or the PR adds a new entry.
- [ ] No `## Change log` / `## Summary` section or per-extraction history line was added to `inst/references/covariate-columns.md`. Per-entry context (derivation rules, scope-promotion rationale, naming-decision sidecars) goes in the H3 entry's Description / Notes / Source aliases. Chronological history is read from `git log`.
- [ ] Source column names different from the canonical names are recorded in `covariateData[[name]]$source_name` and any value transformation (e.g., `SEXM → SEXF` inverts values and flips the effect sign) is documented in `notes`.
- [ ] **Reference categories** for categorical effects match the paper (especially after composite race groups like `RACE_BLACK_OTH` — the reference is everyone NOT in the composite).
- [ ] **Effect form** is correct: multiplicative (`1 + e × COV`), power (`COV^e`), or exponential (`exp(e × COV)`). The form determines what `e` means.
- [ ] **Continuous covariates** are centered / normalized the way the paper describes (e.g., `WT / 70`, `AGE / 40`, `PAGE - 40/4.35`). The skill uses the paper's convention even when it's not "round."
- [ ] **Time-varying vs. time-fixed.** Covariates declared time-varying in the source (e.g., `WT`, `PAGE`, `PNA`) must be supplied at every time point in the event dataset. Fixed covariates (e.g., `GA`) are one-per-subject.

## D. Error model

- [ ] Residual error form matches the paper: proportional, additive, combined, or log-additive. NONMEM "additive on log-scale" = proportional in linear space for nlmixr2.
- [ ] Error magnitude units match the concentration units (e.g., `addSd = 0.231 ug/mL` only if `units$concentration = "ug/mL"`).
- [ ] Conditional error models (different error by study or assay, like `Cirincione_2017_exenatide`, `Kyhl_2016_nalmefene`) use the `| condition` syntax and the condition variable is in `covariateData`.

## E. File plumbing

- [ ] File path is `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`.
- [ ] Function name inside the file **equals** the filename minus `.R`. `buildModelDb()` rejects mismatches.
- [ ] `description`, `reference`, `vignette`, `units`, `covariateData`, `population` all present before `ini()`.
- [ ] Validation vignette lives at `vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd`, and the model file's `vignette <- "..."` value matches that basename (no path, no extension). Confirm `readModelDb("<model>")$meta$vignette` returns the basename after `devtools::load_all()`.
- [ ] `population` uses the extensible schema documented in `references/parameter-names.md` § "File-level metadata"; any paper-specific keys are allowed.
- [ ] No stray `#!` instruction comments from the template remain.
- [ ] Vignette path is `vignettes/articles/<...>.Rmd` (not top-level `vignettes/`), matching the pkgdown "articles" convention used by every non-legacy vignette in the package.
- [ ] If this task `depends_on` an upstream-PK task, the model file's `reference` field cites the upstream model (e.g., `"... PK structure adapted from <Author> <Year>; see modellib('<Upstream>')"`). See `references/model-file-template.md` for the lineage snippet.

## F. Sanity simulations

- [ ] Running `readModelDb("<model>")` returns without error.
- [ ] `rxode2::rxSolve(mod, events)` produces non-NaN, non-negative concentrations across the relevant time window.
- [ ] Simulated Cmax, AUC, and half-life are within ~20% of published values for a typical dose in a typical subject. Larger discrepancies: investigate, don't tune. **Skip this check** when the source publication does not report NCA values (count / Markov / IRT / dropout / TTE modalities, or endogenous models — see § F.1).
- [ ] Rendered vignette contains **no** `Requesting an AUC range starting (0) before the first measurement` warning. If it does, the PKNCA input filter is dropping the time-zero row (commonly `time > 0` or `Cc > 0`); see `pknca-recipes.md` § "Time-zero records (mandatory)" for the fix.
- [ ] The published-vs-simulated comparison is a single combined kable built via `nlmixr2lib::ncaComparisonTable()`; the parameter column header is `NCA parameter`; cells contain friendly labels (`Cmax`, `AUC0-∞ (obs)`, `t½`, …) not raw PKNCA codes (`cmax`, `aucinf.obs`, `half.life`); no "see above" cross-references appear anywhere in the comparison section.
- [ ] A simulated VPC visually resembles the paper's VPC (dose-proportional scaling, right terminal slope, reasonable spread). Skip when the publication has no VPC figure to compare against.
- [ ] If the event table was built from multiple cohorts via `bind_rows()`, ID ranges are disjoint (`anyDuplicated(events[, c("id","time","evid")]) == 0`). See `vignette-template.md`'s `make_cohort(..., id_offset = )` snippet for the pattern.

### F.1 Endogenous / mechanistic models

For models with no dosing (endogenous, mechanistic, steady-state turnover), replace the PK sanity checks above with:

- [ ] Steady-state hold: simulate with `<state>(0) <- <baseline>` and no perturbation; the state stays at baseline within numerical tolerance. Example: `igg_kim_2006` must hold `igg = 12.1` across the full simulation horizon.
- [ ] Perturbation recovery: initialize at `0.5 × baseline` and `2 × baseline`; the trajectory monotonically returns to baseline.
- [ ] Mass-balance / flux check: at steady state, sum every production and elimination flux; the sum equals zero symbolically (not just numerically).
- [ ] All augmentation outputs (e.g., `daily_phe_intake`) have correct units — verified by dimensional analysis, not just plausible magnitude.

See `references/endogenous-validation.md` for full recipes.

### F.2 Count / Markov / IRT / dropout / TTE models

- [ ] Typical-value (no-IIV, no-residual-error) simulation reproduces the published expected-count / hazard / probability trajectory within ~5% at canonical time points.
- [ ] For TTE models: simulated event rate over the observation window matches the published Kaplan-Meier or expected-hazard curve within ~10% across decile time points.
- [ ] For count models: simulated mean count at each scheduled visit matches the published mean within ~10%.
- [ ] For IRT / dropout: predicted item-response or dropout probability at the published reference cohort time-points matches the published value within ~5%.

## G. Registration

- [ ] **Vignette renders end-to-end** (mandatory pre-push gate; HIGHEST priority). Run the procedure in `SKILL.md` Phase 6 Step 2 BEFORE anything else in this section. Record the verbatim `RENDER_GATE stem=... exit=0 html_bytes=...` line; paste it into the PR body. Past PRs that fabricated this evidence broke CI for downstream merges — do not.

- [ ] `nlmixr2lib::buildModelDb()` runs to completion.
- [ ] The new model appears in `modellib()`.
- [ ] `nlmixr2lib::checkModelConventions(model = "<...>")` returns clean (or all deviations are justified in the vignette's Assumptions and deviations section and noted in the PR body).
- [ ] **No non-ASCII characters in the new model file or vignette** (mandatory pre-push gate). Run the bash gate in `SKILL.md` Phase 6 Step 5. The model file's `description <-` string is reified into `data/modeldb.rda`; a single non-ASCII character triggers `R CMD check` `WARNING: found non-ASCII strings`, breaking the R-CMD-check matrix. The vignette doesn't break R CMD check directly but can cause downstream encoding surprises (cross-platform, locale, knitr, pkgdown), so apply the gate to both. When the gate matches, substitute using this table:

   | Non-ASCII | ASCII replacement |
   | --- | --- |
   | `—` (em dash, U+2014) | `--` |
   | `–` (en dash, U+2013) | `-` |
   | `×` (multiplication sign) | `x` (or `*` when between numeric quantities) |
   | `·` (middle dot) | `*` (multiplication) or `.` (decimal-style) |
   | `≈` | `~=` |
   | `→` | `->` |
   | `∈` | `in` |
   | `≤` | `<=` |
   | `≡` | `==` |
   | `…` | `...` |
   | `²` | `^2` |
   | `µ` (Greek mu) | `u` (e.g. `µg/mL` → `ug/mL`) |
   | Greek letters in prose (`η`, `λ`, etc.) | spell out (`eta`, `lambda`) |
   | `§` | `Section ` |

- [ ] `devtools::check()` passes (warnings OK to discuss; errors are blocking). The vignette-render gate above (the top item in § G) is the cheap, fast, vignette-only version of the same check; this `devtools::check()` is the full-package gate that also covers R CMD check warnings and tests. A C-level segfault during `check()` or vignette rendering is a red flag — stop, sidecar-ask the operator to investigate the environment, and do **not** work around with `--no-build-vignettes` or similar flags. See SKILL.md Phase 6 step 2.
- [ ] `NEWS.md` entry added.

### Verifying against the **worktree's** nlmixr2lib, not the system install

When this task runs under `claude_runner`, your `working_dir` is a fresh
git worktree — NOT the system-installed `nlmixr2lib`. Any verification
that uses `library(nlmixr2lib)` or `system.file("modeldb.qs2", package =
"nlmixr2lib")` will read from the system-installed package, which is
almost always **stale** (it predates your extraction and has no
awareness of the new model). You will get a spurious `FALSE` even when
the worktree is correct.

Always verify against the worktree:

```r
# Wrong — reads from the stale system install:
# library(nlmixr2lib)
# "My_2024_drug" %in% nlmixr2lib::modellib()

# Right — loads the worktree's in-development package:
devtools::load_all(".")                       # cwd is the worktree root
"My_2024_drug" %in% nlmixr2lib::modellib()    # now reads worktree modeldb
```

`devtools::check()` already does the right thing (it installs and loads
from the worktree temp build), so use that as the canonical
registration check. The `library(...)` / `qs2::qs_read(system.file(...))`
path is fine only for a quick *sanity smoke* against the old installed
state — don't rely on its `TRUE`/`FALSE` result when the worktree is
under development.

## Common escalations

Ask the user when:

- Parameters look like initial estimates, not final values, and no final-estimate table is available.
- The source paper uses a covariate encoding that conflicts with an existing registered name (e.g., inverted sex convention) and a value transformation is required.
- A new covariate is needed that isn't in `inst/references/covariate-columns.md`.
- The source defines multiple non-hierarchical models and it's unclear which to implement.
- The simulation disagrees with a published figure or NCA table by more than ~20%.
- Units in the source are ambiguous (e.g., "CL = 38.8" without units stated in the immediate context).
- The source uses a parameterization that doesn't map cleanly to nlmixr2 (e.g., custom residual-error shapes that require `add() + prop() + lnorm()` combinations).
