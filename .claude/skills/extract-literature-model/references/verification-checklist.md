# Verification checklist

After the first-pass model file is written, re-read the source independently and walk through every item below. Each item has a failure mode that has caused real translation errors in the past. Fix everything you can; **flag anything ambiguous to the user** using this format:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

Do not silently resolve ambiguity. Do not tune parameters to make a validation output match a target — if the validation disagrees with the paper, investigate the source, not the parameters.

## A. Parameter values

- [ ] Every parameter in `ini()` has an in-file trailing comment pointing to the source location (table, equation, page, or figure). Re-check each comment matches what the source actually says.
- [ ] Values are **final estimates**, not initial estimates. Supplement NONMEM control streams often list initial values in `$THETA` and `$OMEGA`; the final values come from the `$TABLE` output or the main paper. If the only source is a control stream, confirm the values match any published point estimates.
- [ ] **Log-vs-linear reporting.** NONMEM often reports THETAs on the estimation scale (already log), but tables in the paper usually show the back-transformed value. A `log()` wrapper in `ini()` must match what the paper reports: `lcl <- log(0.0388)` is correct when the paper says "CL = 0.0388 L/day."
- [ ] **CV% vs. variance.** `omega²` in NONMEM output is the variance on the internal scale. For log-normal parameters, CV% relates via `omega² = log(CV² + 1)`. Do not paste CV% directly into `ini()` as if it were a variance.
- [ ] **Correlated IIV.** If the paper reports a correlation `r` and individual CV%, the covariance is `cov = r × sqrt(var_1 × var_2)`. Verify the block matrix entries match this formula.
- [ ] **Fixed parameters** the source holds fixed are wrapped in `fixed(...)` in `ini()`.

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

- [ ] Every covariate used in `model()` is registered in `references/covariate-columns.md` with a canonical name, or the PR adds a new entry.
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
- [ ] `description`, `reference`, `units`, `covariateData`, `population` all present before `ini()`.
- [ ] `population` uses the extensible schema documented in `naming-conventions.md`; any paper-specific keys are allowed.
- [ ] No stray `#!` instruction comments from the template remain.

## F. Sanity simulations

- [ ] Running `readModelDb("<model>")` returns without error.
- [ ] `rxode2::rxSolve(mod, events)` produces non-NaN, non-negative concentrations across the relevant time window.
- [ ] Simulated Cmax, AUC, and half-life are within ~20% of published values for a typical dose in a typical subject. Larger discrepancies: investigate, don't tune.
- [ ] A simulated VPC visually resembles the paper's VPC (dose-proportional scaling, right terminal slope, reasonable spread).

### F.1 Endogenous / mechanistic models

For models with no dosing (endogenous, mechanistic, steady-state turnover), replace the PK sanity checks above with:

- [ ] Steady-state hold: simulate with `<state>(0) <- <baseline>` and no perturbation; the state stays at baseline within numerical tolerance. Example: `igg_kim_2006` must hold `igg = 12.1` across the full simulation horizon.
- [ ] Perturbation recovery: initialize at `0.5 × baseline` and `2 × baseline`; the trajectory monotonically returns to baseline.
- [ ] Mass-balance / flux check: at steady state, sum every production and elimination flux; the sum equals zero symbolically (not just numerically).
- [ ] All augmentation outputs (e.g., `daily_phe_intake`) have correct units — verified by dimensional analysis, not just plausible magnitude.

See `references/endogenous-validation.md` for full recipes.

## G. Registration

- [ ] `nlmixr2lib::buildModelDb()` runs to completion.
- [ ] The new model appears in `modellib()`.
- [ ] `devtools::check()` passes (warnings OK to discuss; errors are blocking).
- [ ] Vignette builds cleanly.
- [ ] `NEWS.md` entry added.

## Common escalations

Ask the user when:

- Parameters look like initial estimates, not final values, and no final-estimate table is available.
- The source paper uses a covariate encoding that conflicts with an existing registered name (e.g., inverted sex convention) and a value transformation is required.
- A new covariate is needed that isn't in `references/covariate-columns.md`.
- The source defines multiple non-hierarchical models and it's unclear which to implement.
- The simulation disagrees with a published figure or NCA table by more than ~20%.
- Units in the source are ambiguous (e.g., "CL = 38.8" without units stated in the immediate context).
- The source uses a parameterization that doesn't map cleanly to nlmixr2 (e.g., custom residual-error shapes that require `add() + prop() + lnorm()` combinations).
