---
name: extract-literature-model
description: This skill should be used when the user wants to "add a model from a paper", "extract a pharmacometric model from the literature", "implement a published PK/PD model in nlmixr2lib", or provides a scientific article, conference poster, supplement, or regulatory review document and asks to add the model to nlmixr2lib. Guides source review, standardized model file creation under inst/modeldb/, in-file source-trace verification, and validation vignette with PKNCA NCA checks.
---

# Extract a pharmacometric model from the literature

Input: a scientific source describing a pharmacometric model (journal article, supplement, conference poster, or regulatory review).
Output: a packaged nlmixr2lib model file under `inst/modeldb/`, a validation vignette under `vignettes/`, and updated registry artifacts — opened as a pull request against `main`.

Work through the six phases below. Stop and ask the user at any of the decision points called out explicitly; ambiguity is the main failure mode for this workflow, and silent assumptions are what get shipped as bugs.

## References

Read these on demand; don't load them up front.

- `references/naming-conventions.md` — parameter, compartment, IIV, and error-model names.
- `references/covariate-columns.md` — authoritative register of covariate column names. Consult before introducing any new covariate.
- `references/model-file-template.md` — skeleton for the `.R` file.
- `references/vignette-template.md` — skeleton for the validation vignette.
- `references/pknca-recipes.md` — PKNCA setups for single-dose, steady-state, and multi-dose NCA.
- `references/endogenous-validation.md` — validation strategy for endogenous / mechanistic / turnover models where PKNCA is not the right check.
- `references/verification-checklist.md` — checklist to walk after the first-pass implementation.

## Phase 1 — Source acquisition and scoping

1. Confirm the source type (journal article, supplement, poster, regulatory document).
2. **Always search for supplementary information.** Supplements frequently contain the NONMEM control stream and parameter tables that disambiguate model structure. If the user provided only a main article, ask whether a supplement exists and request it.
3. **Always search for errata, corrigenda, or author corrections.** Check the journal's landing page for the article, the publisher's "corrections" / "notices" feed, and a search like `"<first author> <year> <drug>" erratum` on PubMed and Google Scholar. Ask the user whether they are aware of any corrections if the source is paywalled or the search is inconclusive. **When an erratum revises a value used in the model (parameter estimate, covariate effect, equation, units), the erratum value takes precedence over the main publication.** If multiple errata exist, the most recent supersedes earlier ones. Record the erratum citation in the model file's `reference` metadata alongside the main paper, and in every in-file source-trace comment whose value comes from the erratum, point to the erratum (not the original table).
4. **Verify parameters are final estimates, not initial estimates.** Supplement control streams usually list initial values in `$THETA` / `$OMEGA`; final values come from the paper's results table or `$TABLE` output. If only a control stream is available, confirm values against any published point estimates before treating them as final.
5. **Multiple-model handling.**
   - Base model + final model → extract only the final.
   - Any other "multiple model" case (per-subpopulation, per-endpoint, sensitivity analyses) → list the candidates to the user and ask which to extract. Offer "one," "all," or "a subset."
6. Confirm the target subdirectory under `inst/modeldb/` (usually `specificDrugs/`; endogenous, therapeuticArea, pharmacokinetics, and pharmacodynamics are also valid).

## Phase 2 — Sync with origin/main and branch

Do this before any file changes:

```bash
git fetch origin
git checkout -b <firstauthor>-<year>-<drug> origin/main
```

Local `main` may be stale. The regenerated `data/modeldb.rda` / `inst/modeldb.qs2` must reflect current origin/main or they will clobber models added upstream. Never push directly to `main`; always open a PR.

## Phase 3 — Model file

File path: `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`.

The function name **must equal the filename** minus `.R`; `buildModelDb()` enforces this.

Use `references/model-file-template.md` as the starting skeleton and the two best-formed existing models as anchors:

- `inst/modeldb/specificDrugs/Clegg_2024_nirsevimab.R` — covariates, maturation, correlated IIV, exported race-derivation helper.
- `inst/modeldb/specificDrugs/Hu_2026_clesrovimab.R` — simpler case with Hill-type maturation.

The file body has this shape:

1. `description`, `reference`, `units`, `covariateData`, `population` — metadata before `ini()`.
2. `ini()` — parameters with `label()` and a trailing **in-file comment pointing to the source location** for every value.
3. `model()` — derived terms → individual parameters → micro-constants → ODEs → bioavailability → observation and error.

Follow `references/naming-conventions.md` strictly:

- Structural PK parameters log-transformed: `lka`, `lcl`, `lvc`, `lvp`, `lvp2`, `lq`, `lq2`, `lfdepot`.
- IIV: `eta` + transformed name, e.g., `etalcl` (not `etacl`). Block correlations via `etalcl + etalvc ~ c(var, cov, var)`.
- Residual error: `propSd`, `addSd`. Multi-output: `CcpropSd`, `tumorSizeaddSd`, etc.
- Compartments: `depot`, `central`, `peripheral1`, `peripheral2`, `effect`. Observation: `Cc`.

Covariate columns come from `references/covariate-columns.md`. Before writing any covariate into the file:

- If the canonical name exists, use it and record the source column name in `covariateData[[name]]$source_name`.
- If the source name is an alias of an existing canonical name (e.g., source uses `SEXM`, canonical is `SEXF`), use the canonical name, note the required value transformation (`SEXF = 1 - SEXM`), and **ask the user to confirm the effect-coefficient sign and reference-category implications** before committing.
- If the concept isn't in the register at all, propose a new entry (canonical name, description, units, type, reference category, source aliases) and ask the user to confirm before adding it. The new entry is committed alongside the model.

`population` uses the extensible schema in `naming-conventions.md`. Common fields: `n_subjects`, `n_studies`, `age_range`, `weight_range`, `sex_female_pct`, `race_ethnicity`, `disease_state`, `dose_range`, `regions`, `notes`. Add any additional keys the paper describes (e.g., `ga_range`, `renal_function`, `co_medication`) — do not force facts into the common schema.

## Phase 4 — Verification (re-read the source)

After the first pass, re-read the source independently and walk through `references/verification-checklist.md`. Common pitfalls:

- NONMEM THETA log-vs-linear reporting and `omega² = log(CV² + 1)` for log-normal variance.
- NONMEM "additive on log-scale" ≡ proportional in nlmixr2's linear space.
- Reference weight / age for allometric and maturation terms (70 kg? 5 kg? 40 weeks PMA?).
- Reference category for categorical effects after composite-group renames.
- Bioavailability target compartment.
- Units consistency (dose × F ÷ V must give the declared concentration units).

Every parameter's in-file source-trace comment must be verified — the comment states where it came from, so checking is a line-by-line audit. 

When any uncertainty remains, ask the user using this format:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

Never silently resolve ambiguity. Never tune parameter values to match a validation target.

## Phase 5 — Validation vignette

File path: `vignettes/<FirstAuthor>_<Year>_<drug>.Rmd`, with matching `VignetteIndexEntry`.

Use `references/vignette-template.md`. Required sections, in order:

1. **Header and setup** — libraries include `nlmixr2lib`, `PKNCA`, `rxode2`, `dplyr`, `ggplot2`.
2. **Population** — narrative reproducing the `population` metadata; cite the source table listing baseline demographics.
3. **Source trace** — a dedicated table listing the source location (page / table / equation / figure) for every model equation and every `ini()` parameter. This is in addition to the in-file comments; the vignette gives reviewers a single place to audit provenance.
4. **Virtual cohort** — covariate distributions match the population metadata. Use WHO weight-for-age curves for pediatric models.
5. **Simulation** — `rxode2::rxSolve(mod, events)` for stochastic VPCs; `rxode2::zeroRe()` + `rxSolve` for typical-value replications.
6. **Replicate published figures** — one code chunk per figure, caption linking to the source figure number ("Replicates Figure 4 of <Author Year>").
7. **PKNCA validation** — required; no inline trapezoidal NCA. See `references/pknca-recipes.md`. The PKNCA formula **must include a treatment grouping variable** (`conc ~ time | id/treatment`) so per-group results can be compared against the paper.
8. **Comparison against published NCA** — if the source paper reports Cmax / Tmax / AUC / half-life, render a side-by-side comparison table. Flag differences > 20% in the narrative and investigate the source — do not tune.
9. **Assumptions and deviations** — explicit list of what you had to assume because the paper didn't say (race distribution, z-score stability, etc.).

For **endogenous / turnover models** where NCA isn't the right validation, replace the PKNCA section with the steady-state / perturbation-recovery / mass-balance checks described in `references/endogenous-validation.md`. For **multi-output models**, run one PKNCA block per output.

### Endogenous and mechanistic models

Papers that describe endogenous turnover, steady-state balance, or mechanistic enzyme kinetics (e.g., Kim 2006 IgG FcRn recycling, Charbonneau 2021 phenylalanine) have a different shape than drug PK models:

- Parameters are mechanistic constants (`Vmax`, `Km`, `kint`, `kcat`, `kpro`, baseline concentrations `bl_<species>`, fractional-activity scalars like `f_<enzyme>`) rather than log-transformed CL/V.
- `ini()` usually has **no IIV etas and no residual error** — the model describes population-typical mechanism, not variability.
- `model()` has **no dosing events**; the state starts at a biological baseline (`<state>(0) <- bl_<state>` or `<state>(0) <- css`).
- Validation is **not** PKNCA. Use steady-state / perturbation-recovery / mass-balance checks. See `references/endogenous-validation.md`.
- Dimensional analysis is load-bearing. These models often mix `mg/mL`, `mg/kg`, `L/kg`, `1/day`; a single unit slip silently corrupts the balance. Walk through every term in every ODE and the derived rates.

Naming conventions for mechanistic parameters are documented in `references/naming-conventions.md` under "Endogenous / mechanistic parameters."

## Phase 6 — Registration, tests, docs, PR

1. Re-confirm the branch is on top of fresh `origin/main` (`git fetch origin && git rebase origin/main` if needed).
2. Run `nlmixr2lib::buildModelDb()` to regenerate `data/modeldb.rda` and `inst/modeldb.qs2`. Confirm the new model appears in `modellib()`.
3. Run `devtools::check()`. Vignettes must build cleanly.
4. Add a `NEWS.md` entry under the current development version:
   ```
   - Added <Author> <Year> population PK model for <drug> (#PR).
   ```
5. Commit the model file, the vignette, the regenerated `modeldb.rda` / `modeldb.qs2`, the `NEWS.md` entry, and any updates to `references/covariate-columns.md` (if a new covariate was registered) together on the feature branch.
6. Push the branch and open a PR against `main`. Use `gh pr create` with a title like `Add <Author> <Year> <drug> model`.

## Stop-and-ask triggers (consolidated)

Don't guess — ask the user when:

- The source has multiple non-hierarchical models and it's not obvious which to extract.
- Parameter values look like initial estimates rather than final.
- Covariate encoding isn't fully specified (reference category, units, transformation).
- A source column name is not in `references/covariate-columns.md` (propose a new entry and confirm).
- A source column is an alias of an existing canonical name and the mapping involves value inversion or a reference-category flip.
- A parameter name deviates from the nlmixr2lib standard (propose the canonical name and confirm).
- PKNCA output disagrees with a published NCA table by more than ~20% after careful review.
- The source is paywalled and the user hasn't supplied the text.
- An erratum search is inconclusive (e.g., paywalled journal, ambiguous correction notice) — ask the user to confirm whether any corrections apply.

Use this fixed format for ambiguities:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

## Constraints

- Never invent parameter values; if it's not in the source, ask.
- Never tune parameters to make a validation output match a target.
- This skill **only adds** new models. Retrofitting existing vignettes to use PKNCA, or renaming covariates in existing files, is a future separate skill.
- Never push directly to `main`. Open a PR.
