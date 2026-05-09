# Verification checklist

After the first-pass model file is written, re-read the source independently and walk through every item below. Each item has a failure mode that has caused real translation errors in the past. Fix everything you can; **flag anything ambiguous to the user** using this format:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

Do not silently resolve ambiguity. Do not tune parameters to make a validation output match a target — if the validation disagrees with the paper, investigate the source, not the parameters.

## H. Source identity

- [ ] On-disk PDF / XML title, first author, journal, and year match the task's `Paper metadata` block. Filename PMID matches the actual PMID in the source (catches mislabelled drops like `PMID_23436260.pdf` actually being Frey 2010 / PMID 20097931).
- [ ] Drug named in the task metadata matches the molecule the paper actually models.

## A. Parameter values

- [ ] **Errata / corrigenda checked.** Confirmed search for published corrections to the source; any erratum value supersedes the main publication for conflicting values, with the most recent erratum winning when multiple exist. Erratum citation recorded in the model file's `reference` field, and each affected `ini()` comment points to the erratum.
- [ ] Every parameter line in `ini()` has a clear provenance comment if it did not come from the paper's text/tables (author correspondence, figure-digitisation, upstream-task model file). The comment cites the followup-register entry (`tracking/operator_followups.md F<n>`) when applicable.
- [ ] Every parameter in `ini()` has an in-file trailing comment pointing to the source location (table, equation, page, or figure). Re-check each comment matches what the source actually says.
- [ ] Values are **final estimates**, not initial estimates. Supplement NONMEM control streams often list initial values in `$THETA` and `$OMEGA`; the final values come from the `$TABLE` output, an accompanying `.lst` `FINAL PARAMETER ESTIMATE` block, or the main paper. If the only source is a control stream, confirm the values match any published point estimates.
- [ ] **Log-vs-linear reporting.** NONMEM often reports THETAs on the estimation scale (already log), but tables in the paper usually show the back-transformed value. A `log()` wrapper in `ini()` must match what the paper reports: `lcl <- log(0.0388)` is correct when the paper says "CL = 0.0388 L/day."
- [ ] **CV% vs. variance.** `omega²` in NONMEM output is the variance on the internal scale. For log-normal parameters, CV% relates via `omega² = log(CV² + 1)`. Do not paste CV% directly into `ini()` as if it were a variance.
- [ ] **Correlated IIV.** If the paper reports a correlation `r` and individual CV%, the covariance is `cov = r × sqrt(var_1 × var_2)`. Verify the block matrix entries match this formula.
- [ ] **Fixed parameters** the source holds fixed are wrapped in `fixed(...)` in `ini()` — applies to ALL parameter types (THETAs, allometric exponents, IIVs, residual errors, covariate effects, bioavailability anchors), not just IIVs. Source signals: explicit "fixed at <value>" prose, NONMEM `FIX` flags on `$THETA`/`$OMEGA`/`$SIGMA`, allometric exponents reported without uncertainty, bioavailability `F1=1` set as structural anchor, parameters inherited from upstream papers without re-fitting. If a parameter is reported without uncertainty but the paper does not explicitly say "fixed", sidecar-ask before guessing — see `SKILL.md` § "Fixed parameters in `ini()`" for the encoding examples.

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
- [ ] `population` uses the extensible schema documented in `naming-conventions.md`; any paper-specific keys are allowed.
- [ ] No stray `#!` instruction comments from the template remain.
- [ ] Vignette path is `vignettes/articles/<...>.Rmd` (not top-level `vignettes/`), matching the pkgdown "articles" convention used by every non-legacy vignette in the package.
- [ ] If this task `depends_on` an upstream-PK task, the model file's `reference` field cites the upstream model (e.g., `"... PK structure adapted from <Author> <Year>; see modellib('<Upstream>')"`). See `references/model-file-template.md` for the lineage snippet.

## F. Sanity simulations

- [ ] Running `readModelDb("<model>")` returns without error.
- [ ] `rxode2::rxSolve(mod, events)` produces non-NaN, non-negative concentrations across the relevant time window.
- [ ] Simulated Cmax, AUC, and half-life are within ~20% of published values for a typical dose in a typical subject. Larger discrepancies: investigate, don't tune. **Skip this check** when the source publication does not report NCA values (count / Markov / IRT / dropout / TTE modalities, or endogenous models — see § F.1).
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

- [ ] `nlmixr2lib::buildModelDb()` runs to completion.
- [ ] The new model appears in `modellib()`.
- [ ] `nlmixr2lib::checkModelConventions(model = "<...>")` returns clean (or all deviations are justified in the vignette's Assumptions and deviations section and noted in the PR body).
- [ ] **Verify no non-ASCII characters in the new model file or vignette** (mandatory pre-push gate). The model file's `description <-` string is reified into `data/modeldb.rda`, and a single non-ASCII character there triggers `R CMD check` `WARNING: found non-ASCII strings`, breaking the R-CMD-check matrix. The vignette doesn't break R CMD check directly but can cause downstream encoding surprises (cross-platform, locale, knitr, pkgdown), so apply the gate to both.

   ```bash
   for f in inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R \
            vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd; do
     LC_ALL=C grep -nP "[^\x00-\x7F]" "$f" && { echo "FAIL: non-ASCII in $f"; exit 1; }
   done
   ```

   No matches → pass; replace any matches with ASCII (`—` → `--`, `–` → `-`, `×` → `x`, `·` → `*`, `≈` → `~=`, `→` → `->`, `≤` → `<=`, `…` → `...`, `²` → `^2`, `µ` → `u`, Greek letters in prose spelled out, `§` → `Section `). See SKILL.md Phase 6 step 5 for the full substitution table and rationale.

- [ ] **Render the new vignette end-to-end** (mandatory pre-push gate). Run the exact command below on the worktree, with a 5-minute wall-clock budget — do not skip it, do not interpret silence as success, and do not push if it errors. This is the single highest-value gate the worker performs: it catches the failure modes that pkgdown CI also surfaces (missing data columns, time-varying covariates assigned to rxEt objects that get silently dropped, PKNCA formulas referencing absent columns, simulation crashes, vignettes exceeding the time budget).

   ```bash
   timeout 300 Rscript --vanilla -e "
     pkgload::load_all('.', quiet = TRUE)
     rmarkdown::render(
       'vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd',
       quiet = FALSE
     )
   "
   ```

   Interpret the exit code:
   - `0` with output HTML present → gate passed; proceed.
   - non-zero R error → fix the vignette before committing. The traceback identifies the failing chunk; the most common causes are listed in SKILL.md Phase 6 step 4.
   - `124` (timeout) → reduce simulation size; do not skip the gate.
   - C-level segfault → broken environment; stop and sidecar-ask. Do **not** work around with `--no-build-vignettes`.

- [ ] `devtools::check()` passes (warnings OK to discuss; errors are blocking). A C-level segfault during `check()` or vignette rendering is a red flag — stop, sidecar-ask the operator to investigate the environment, and do **not** work around with `--no-build-vignettes` or similar flags. See SKILL.md Phase 6 step 4.
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
