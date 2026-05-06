---
name: extract-literature-model
description: This skill should be used when the user wants to "add a model from a paper", "extract a pharmacometric model from the literature", "implement a published PK/PD model in nlmixr2lib", or provides a scientific article, conference poster, supplement, regulatory review document, or DDMORE Foundation Model Repository bundle (NONMEM `.mod` / `.ctl` control stream + `.lst` listing + RDF metadata + Model_Accomodations file) and asks to add the model to nlmixr2lib. Guides source review, standardized model file creation under inst/modeldb/, in-file source-trace verification, and validation vignette with PKNCA NCA checks.
---

# Extract a pharmacometric model from the literature

Input: a scientific source describing a pharmacometric model. Two source shapes are supported:

- **Paper-source**: journal article, supplement, conference poster, or regulatory review.
- **DDMORE-source**: a directory from the DDMORE Foundation Model Repository (typically a NONMEM `.mod` / `.ctl` control stream plus `.lst` listing, with optional RDF metadata and `Model_Accomodations` text file). DDMORE-source extractions land under `inst/modeldb/ddmore/`. See `references/ddmore-source.md` for the full DDMORE workflow.

Output: a packaged nlmixr2lib model file under `inst/modeldb/`, a validation vignette under `vignettes/articles/`, and updated registry artifacts — opened as a pull request against `main`.

> **DDMORE-source support is a temporary one-time skill extension.**
> The current revision adds extraction of the
> `dpastoor/ddmore_scraping` bundle (~58 DDMORE Foundation Model
> Repository models) into `inst/modeldb/ddmore/`. Once that batch is
> done, the DDMORE-aware sections of this file and its references —
> together with `references/ddmore-source.md` — are scheduled for
> removal. See `references/ddmore-source.md` § "Cleanup recipe" for
> the rollback steps. Until removal, the skill works for both
> paper-source and DDMORE-source extractions.

Work through the six phases below. Stop and ask the user at any of the decision points called out explicitly; ambiguity is the main failure mode for this workflow, and silent assumptions are what get shipped as bugs.

## References

Read these on demand; don't load them up front.

- `references/naming-conventions.md` — parameter, compartment, IIV, and error-model names. Includes a NONMEM → nlmixr2 syntax-translation section consulted by DDMORE-source extractions.
- `inst/references/covariate-columns.md` — authoritative register of covariate column names. Consult before introducing any new covariate. (Installed with the package so R code like `checkModelConventions()` can parse it at runtime.)
- `references/model-file-template.md` — skeleton for the `.R` file.
- `references/vignette-template.md` — skeleton for the validation vignette.
- `references/pknca-recipes.md` — PKNCA setups for single-dose, steady-state, and multi-dose NCA.
- `references/endogenous-validation.md` — validation strategy for endogenous / mechanistic / turnover models where PKNCA is not the right check.
- `references/verification-checklist.md` — checklist to walk after the first-pass implementation.
- `references/ddmore-source.md` — extra guidance for DDMORE Foundation Model Repository sources (NONMEM `.mod` / `.ctl` / MDL XML / `.txt` NMTRAN executables plus `.lst` listings, `DDMODEL00000<id>.rdf`, `Model_Accomodations.text|.txt`, `<id>.json`). **Read only when the source is a DDMORE bundle directory.**

## Phase 1 — Source acquisition and scoping

1. Confirm the source type. Two shapes are supported:
   - **Paper-source** — journal article, supplement, poster, regulatory document. Continue with steps 2–11 below.
   - **DDMORE-source** — directory from the DDMORE Foundation Model Repository. Telltale files: `Executable_*.{mod,ctl,xml,txt}`, `Output_real_*.lst`, `Output_simulated_*.lst`, optionally `DDMODEL00000<id>.rdf`, `Model_Accomodations.text` (or `.txt`), `<id>.json`. **Switch to `references/ddmore-source.md`** for the DDMORE-specific Phase 1 steps; come back to Phase 2 onward in this file when you reach the model-file-creation phase. The DDMORE flow replaces step 2 (paper identity), step 3 (species check), step 4 (trimmed markdown), step 5 (full-text check), and step 7 (supplement search) with DDMORE-tailored equivalents documented in `ddmore-source.md`. Steps 6 (upstream-popPK), 8 (errata), 9 (final estimates), 10 (multiple-model handling), and 11 (target subdirectory — DDMORE-source models always go to `ddmore/`) still apply, with DDMORE-aware tweaks.
2. **(Paper-source)** **Verify the on-disk file is the paper the task names.** Open the source file (or its trimmed `.md` companion) and read the title + first-author line + journal + year. Compare against the task's `Paper metadata` block. If any of {first-author, year, journal, drug} disagrees, stop and sidecar-ask:

   > The on-disk file `<filename>` reports `<actual-author> <actual-year>, <actual-journal>: "<actual-title>"` but the task names `<expected-author> <expected-year>, ... <expected-drug>`. Options: (A) the task metadata is wrong and the file is the right paper — please correct the task metadata, (B) the file is mislabelled — please drop the correct PDF, (C) skip this task. Which applies?

   This catches both "wrong-PMID-in-filename" (e.g. Frey 2010 saved as `PMID_23436260.pdf`) and "wrong-drug-guessed" task-generator errors. Do not proceed with extraction until the mismatch is resolved.

3. **Check the study population species.** Skim the Methods / Subjects section. If every PK dataset contributing to the final model is non-human (rat, mouse, cyno, dog, etc.), stop and sidecar-ask the operator before drafting anything:

   > This paper reports a preclinical-only (<species>) population PK model. nlmixr2lib is a library of human population-PK models. Should I (A) extract it anyway with clear preclinical metadata, (B) extract only a human-scaled projection if the paper includes one, or (C) skip this paper?

   If the paper has both preclinical and human cohorts and the final model is fit to pooled or human-only data, proceed normally — this trigger is specifically for animal-only final models. Record the operator's decision in the PR body.

4. **Prefer trimmed markdown when available.** The preprocessor at `mab_human_consensus/tracking/preprocess_papers.py` writes a `<stem>_trimmed.md` next to each source file (PMC XML, PDF, DOCX, XLSX) containing only the sections the extraction actually needs: Title + Abstract + Methods + Results + Tables + Figure captions. The Introduction, Discussion, Conclusions, References, Acknowledgments, and publisher boilerplate are stripped. If `PMID_<pmid>_pmc_trimmed.md` (or `PMID_<pmid>_trimmed.md` for a PDF, or `<stem>_trimmed.md` for a supplement) exists, read it **instead of** the raw `.xml` / `.pdf` / `.docx` — it's typically 40-95% smaller with no loss of extractable content. Full-text sanity check on the trimmed file: ~15 KB+ (full-text trim) vs < 3 KB (abstract-only trim). Fall back to the raw source only if the `_trimmed.md` doesn't exist, the trim appears to have lost a specific piece of information you need (rare — only when the paper is structurally unusual), or you explicitly need the discarded sections (e.g., to quote a Discussion claim in the vignette narrative).

5. **Verify the source contains full text, not just the abstract.** Wiley / BJCP and some other publishers serve PMC XML containing only front matter + abstract. Before reading for model structure, run a quick sanity check (against the trimmed file if present, otherwise the raw source):

   - Trimmed `.md` file size ≥ ~15 KB, or raw PMC XML ≥ ~40 KB (full-text XML is typically 100 KB+; abstract-only is usually < 20 KB).
   - The file contains a materially-present Methods section (not just a "Methods" heading followed by one abstract paragraph).
   - If only a PDF is on disk, confirm it runs past the abstract (multi-page, Methods / Results / Tables present).

   If only the abstract is available, sidecar-ask:

   > The source on disk for <paper> contains only the abstract and front matter; full text appears to be blocked by the publisher. Options: (A) pause this task until full text is provided, (B) proceed only if a supplement / regulatory review on disk contains the model equations and parameter tables, (C) skip this paper. Which applies?

   Never attempt extraction from an abstract alone — population-PK parameter values, covariate effects, and equations are not in an abstract.

6. **Detect upstream-popPK dependencies.** Skim Methods for phrases like "PK was described using the popPK model previously developed from <study/phase>", "the structural PK model was fixed from <reference>", "covariate effects were carried over from <author> et al.", or "the PK model from <prior publication> was used as a backbone." If the current paper's PD model fixes its PK parameters from a separate publication that is not on disk:

   1. Try to identify the upstream paper from the references list (PMID, DOI, or full citation).
   2. If identifiable: sidecar-ask whether the operator wants this task to acquire `depends_on: [<upstream-task>]` and pause until the upstream task completes — versus extracting only the PD layer with the upstream PK parameters reproduced inline.
   3. If unidentifiable (e.g. "the popPK model from internal Phase 1/2 studies" with no specific citation): sidecar-ask whether to (A) skip the task, (B) proceed with parameters fixed inline as reported in the current paper, with a clear "upstream PK source not located" note in the vignette Errata, or (C) defer pending operator investigation.

   Never silently fabricate upstream PK parameters from training data.

7. **Always search for supplementary information.** Supplements frequently contain the NONMEM control stream and parameter tables that disambiguate model structure. If the user provided only a main article, ask whether a supplement exists and request it.
8. **Always search for errata, corrigenda, or author corrections.** Check the journal's landing page for the article, the publisher's "corrections" / "notices" feed, and a search like `"<first author> <year> <drug>" erratum` on PubMed and Google Scholar. Ask the user whether they are aware of any corrections if the source is paywalled or the search is inconclusive. **When an erratum revises a value used in the model (parameter estimate, covariate effect, equation, units), the erratum value takes precedence over the main publication.** If multiple errata exist, the most recent supersedes earlier ones. Record the erratum citation in the model file's `reference` metadata alongside the main paper, and in every in-file source-trace comment whose value comes from the erratum, point to the erratum (not the original table).
9. **Verify parameters are final estimates, not initial estimates.** Supplement control streams usually list initial values in `$THETA` / `$OMEGA`; final values come from the paper's results table or `$TABLE` output. If only a control stream is available, confirm values against any published point estimates before treating them as final. **For DDMORE-source extractions**, final estimates are in the `Output_*.lst` listing (the `FINAL PARAMETER ESTIMATE` block after `MINIMIZATION SUCCESSFUL`); the `.mod` `$THETA` / `$OMEGA` / `$SIGMA` blocks are initial values. See `ddmore-source.md` for the parsing recipe.
10. **Multiple-model handling.**
    - Base model + final model → extract only the final.
    - Any other "multiple model" case (per-subpopulation, per-endpoint, sensitivity analyses) → list the candidates to the user and ask which to extract. Offer "one," "all," or "a subset."
11. Confirm the target subdirectory under `inst/modeldb/`:
    - **Paper-source** — usually `specificDrugs/`; `endogenous/`, `therapeuticArea/`, `pharmacokinetics/`, `pharmacodynamics/` also valid.
    - **DDMORE-source** — always `ddmore/` (parallel to `specificDrugs/`).

## Phase 2 — Sync with origin/main and branch

Do this before any file changes:

```bash
git fetch origin
git checkout -b <firstauthor>-<year>-<drug> origin/main
```

Local `main` may be stale. The regenerated `data/modeldb.rda` / `inst/modeldb.qs2` must reflect current origin/main or they will clobber models added upstream. Never push directly to `main`; always open a PR.

(When this skill runs as a task under `claude_runner`, the runner's preamble
will note that a fresh worktree has already been set up — verify with
`git branch --show-current` and skip the `git fetch` / `git checkout -b`
in that case. The runner provides the authoritative instructions for its
own environment.)

**Worktree resumption — handle pre-existing WIP.** Before drafting anything, run `git status -s` in the worktree. If the working tree is **not clean** (uncommitted modifications or untracked files), this worktree is a resumption of a prior task instance that crashed or was interrupted:

1. Read every modified / untracked file under `inst/modeldb/specificDrugs/`, `inst/modeldb/<other-categories>/`, `vignettes/articles/`, `data/`, `inst/`, and `NEWS.md`.
2. Decide whether the prior WIP is salvageable (sound enough to continue from where it stopped) or unsalvageable (revert to clean state via `git restore .` and `git clean -fd`, then start over).
3. If salvageable: continue from where the prior run left off and note the resumption in the PR body so the reviewer is aware.
4. If unsalvageable: clean and restart. No operator approval is required for the cleanup itself; the reset is automatic for an unsalvageable WIP.

If the worktree's branch was **already committed and pushed in a prior run** (i.e. `git log origin/<branch>..HEAD` shows no new commits and the branch exists on origin with task content), that is a different case — the task was already completed previously. Sidecar-ask:

> Worktree branch `<branch>` is already pushed with prior task content. Options: (A) verify the pushed branch and exit cleanly without re-extracting, (B) tear down and re-extract from scratch (operator may have new source files / new skill version). Which applies?

## Phase 3 — Model file

File path:

- **Paper-source**: `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`.
- **DDMORE-source**: `inst/modeldb/ddmore/<FirstAuthor>_<Year>_<drug>.R`.

If the chosen `<FirstAuthor>_<Year>_<drug>` name collides with an existing file in `ddmore/` (rare — e.g., two Themans 2019 meropenem entries with different scenarios), append a lowercase letter to the year: `Themans_2019a_meropenem.R`, `Themans_2019b_meropenem.R`. Use the same year-letter for both files in the pair so the chronological ordering is preserved.

The function name **must equal the filename** minus `.R`; `buildModelDb()` enforces this.

Use `references/model-file-template.md` as the starting skeleton and the two best-formed existing models as anchors:

- `inst/modeldb/specificDrugs/Clegg_2024_nirsevimab.R` — covariates, maturation, correlated IIV, exported race-derivation helper.
- `inst/modeldb/specificDrugs/Hu_2026_clesrovimab.R` — simpler case with Hill-type maturation.

The file body has this shape:

1. `description`, `reference`, `vignette`, `units`, `covariateData`, `population` — metadata before `ini()`. `vignette` is the basename of the validation vignette in `vignettes/articles/` (e.g., `"Clegg_2024_nirsevimab"`, no path, no extension); `buildModelDb()` extracts it so the list-of-models table can link to the rendered vignette on the pkgdown site.

   **DDMORE-source models additionally set:**
   - `ddmore_id <- "DDMODEL00000<NNN>"` — the DDMORE Foundation ID, for repo traceability. Always present on `ddmore/` models.
   - `replicate_of <- "inst/modeldb/specificDrugs/<Other>_<Year>_<drug>.R"` — relative path to a paper-derived counterpart in `specificDrugs/` when one exists (i.e., the same publication is also extracted there). Reciprocal — the `specificDrugs/` counterpart adds `replicate_of <- "inst/modeldb/ddmore/<This>_<Year>_<drug>.R"` pointing back. Omit when there is no counterpart.

2. `ini()` — parameters with `label()` and a trailing **in-file comment pointing to the source location** for every value.
3. `model()` — derived terms → individual parameters → micro-constants → ODEs → bioavailability → observation and error.

Follow `references/naming-conventions.md` strictly:

- Structural PK parameters log-transformed: `lka`, `lcl`, `lvc`, `lvp`, `lvp2`, `lq`, `lq2`, `lfdepot`.
- IIV: `eta` + transformed name, e.g., `etalcl` (not `etacl`). Block correlations via `etalcl + etalvc ~ c(var, cov, var)`.
- Residual error: `propSd`, `addSd`. Multi-output: `CcpropSd`, `tumorSizeaddSd`, etc.
- Compartments: `depot`, `central`, `peripheral1`, `peripheral2`, `effect`. Observation: `Cc`.

Covariate columns come from `inst/references/covariate-columns.md`. Before writing any covariate into the file:

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

**Missing-parameter / author-correspondence pathway.** If a parameter you need to populate the model is **not present anywhere on disk** (not in the main paper, not in any supplement, not in any associated regulatory review), do not substitute a training-data value or a "typical" textbook value. Sidecar-ask:

> Parameter `<name>` (e.g. `kdeg`, `V_DXd`) is required for the <full-TMDD / payload-PK / structural> model but is not reported in any source on disk for <paper>. Options: (A) draft an author-correspondence email and pause this task pending reply (operator handles the email), (B) approximate with <QSS / steady-state / fixed-from-class> and document the approximation in vignette Errata, (C) skip this paper. Which applies?

Operator-followup tracking: unresolved-parameter cases are recorded in `tracking/operator_followups.md` (the `F1`, `F2`, ... numbered register) so the operator can batch-send author emails. When a reply arrives with a numeric value, the operator records it in the followups file; this skill does NOT email authors.

**Non-paper provenance — annotate inline.** When a parameter value did not come from the paper's text or tables — e.g. the operator read it off a graphical figure, an author supplied it via email correspondence, or it was carried from an upstream-task model file — record the provenance as an inline comment on the parameter line so the source-trace is unambiguous:

```r
ini({
  # paper-derived (Table 4)
  lcl <- log(0.225) ; label("Clearance (L/day)")

  # author correspondence (J. Almquist email 2026-04-29);
  # see tracking/operator_followups.md F12
  lkdeg <- log(0.0231) ; label("Receptor degradation rate (1/day)")

  # operator-extracted from Figure 2 (digitised); ±10% uncertainty
  lvdxd <- log(0.038) ; label("DXd payload volume (L/kg) — figure-derived")
})
```

This is mandatory for any parameter not from the paper text/tables. The vignette's Assumptions and deviations / Errata section must also list the non-paper provenance (see `references/vignette-template.md`).

## Phase 5 — Validation vignette

File path: `vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd`, with matching `VignetteIndexEntry`. Drug-specific vignettes live under `vignettes/articles/` so pkgdown builds them for the site but CRAN does not — `.Rbuildignore` excludes that directory. The basename (without `.Rmd`) must match the `vignette <- "..."` field in the model file. Only the legacy `PK_2cmt_mAb_Davda_2014.Rmd` remains at top-level `vignettes/`.

**Validation strategy is gated on model type, not source shape:**

- **PK / PD models with published NCA** (paper-source or DDMORE-source) — full PKNCA section with side-by-side comparison against the publication. Default path.
- **Endogenous / mechanistic / steady-state models** — replace PKNCA with the recipes in `references/endogenous-validation.md`.
- **Count / Markov / IRT / dropout / TTE models** (most common for DDMORE-source non-PK entries) — replace PKNCA with mechanistic sanity checks: simulation reproduces the typical-value trajectory, baseline / steady-state probability matches the source, and (for TTE) the simulated event rate matches the published Kaplan-Meier or hazard. See `references/ddmore-source.md` § "Validation for non-PK/PD model classes" for templates. Skipping PKNCA without one of these substitutes is not allowed; if you can't decide which substitute applies, sidecar-ask the operator.
- **DDMORE-source models with no linked publication** — comparison against published figures is not possible. Validation reduces to (a) the model parses, (b) `rxode2::rxSolve()` runs to completion on the simulated dataset shipped in the DDMORE bundle, (c) the simulated trajectory visually matches the bundle's `Output_simulated_*.lst`. Document the absence of a paper in the vignette's Assumptions and deviations section.

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
2. Run `nlmixr2lib::buildModelDb()` to regenerate `data/modeldb.rda` and `inst/modeldb.qs2`. Confirm the new model appears in `modellib()`. **When verifying in R, do `devtools::load_all(".")` first so `modellib()` reads the worktree's in-development package, not the stale system install** — see `references/verification-checklist.md` § "Verifying against the worktree's nlmixr2lib" for why a bare `library(nlmixr2lib)` can return a misleading `FALSE`.
3. Run `nlmixr2lib::checkModelConventions(model = "<FirstAuthor>_<Year>_<drug>")` and review the output. Any deviations from the canonical parameter / IIV / residual-error / covariate / compartment conventions (see `references/naming-conventions.md` and `inst/references/covariate-columns.md`) that the function flags should be either fixed in the model file before committing, or explicitly justified in the vignette's Assumptions and deviations section. `buildModelDb()` runs `checkModelConventions()` implicitly at package-build time, but running it explicitly on your new model makes drift visible before commit, not after CI. Paste the key lines of the output into the PR body so a reviewer can see what was checked.
4. **Render the vignette locally and verify it completes without error in under 5 minutes of wall-clock time.** This is a mandatory pre-push gate — it catches missing columns, wrong variable names, PKNCA formula errors, and simulation crashes before CI does.

   ```bash
   timeout 300 Rscript --vanilla -e "
     pkgload::load_all('.', quiet = TRUE)
     rmarkdown::render(
       'vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd',
       quiet = FALSE
     )
   "
   ```

   Interpret the result:

   - **Exit 0, output HTML present** → vignette is clean; proceed to commit.
   - **Exit non-zero with an R error** → fix the vignette before committing. Read the traceback carefully; the most common causes are (a) a `select()` / PKNCA formula referencing a column that was never created (add it in the preceding `mutate()`), (b) a variable name used before it is defined (reorder chunks), (c) a simulation that errors because covariate values are out of range for the model.
   - **Exit 124 (timeout)** → the vignette exceeds 5 minutes. Reduce simulation size: cut `nSub` (stochastic VPC subjects), shorten the observation grid, or move expensive chunks behind `eval = FALSE` with a note. Do **not** skip the time budget — pkgdown CI has strict wall-time limits and a slow vignette breaks the build for everyone.
   - **C-level segfault (`*** caught segfault ***`)** → broken R / rxode2 / nlmixr2 install, not a model-file problem. Stop, sidecar-ask the operator to fix the environment; do not work around with `--no-build-vignettes`.

5. Run `devtools::check()`. Vignettes must build cleanly.

   A C-level segfault (`*** caught segfault ***`) during `check()` or vignette rendering is a red flag — it indicates a broken R / rxode2 / nlmixr2 install in the environment, not a model-file problem. Stop, sidecar-ask the operator to investigate and fix the environment, and do not work around it with `--no-build-vignettes` or similar flags.
6. Add a short, single-line `NEWS.md` entry under the current development
   version. The goal is a scannable changelog — the model file and vignette
   already contain the full detail, so NEWS should only mention:
   - the drug,
   - a minimal reference (author + year + DOI link — not a full citation),
   - and the population studied.

   Format:
   ```
   - Add <Author> <Year> <drug> ([doi:<doi>](https://doi.org/<doi>)) — <population phrase>.
   ```

   Examples:
   ```
   - Add Xu 2019 sarilumab ([doi:10.1007/s40262-019-00765-1](https://doi.org/10.1007/s40262-019-00765-1)) — adults with rheumatoid arthritis.
   - Add Clegg 2024 nirsevimab ([doi:10.1007/s40262-024-01387-y](https://doi.org/10.1007/s40262-024-01387-y)) — preterm and term infants.
   ```

   For DDMORE-source models, append `[DDMODEL00000<NNN>]` after the
   citation:

   ```
   - Add Themans 2019 meropenem ([doi:10.1111/bcp.14025](https://doi.org/10.1111/bcp.14025)) [DDMODEL00000301] — adults with severe pneumonia.
   ```

   Do NOT include: covariate list, compartment count, IIV structure,
   residual-error form, data origin, study counts, PKNCA sentence, or
   anything else that lives in the model file's metadata or vignette. A
   reviewer who wants those details clicks through to the model file.
7. Commit the model file, the vignette under `vignettes/articles/`, the regenerated `modeldb.rda` / `modeldb.qs2` / `modeldb.Rd`, the `NEWS.md` entry, and any updates to `inst/references/covariate-columns.md` (if a new covariate was registered) together on the feature branch.
8. Push the branch and open a PR against `main`. Use `gh pr create` with a title like `Add <Author> <Year> <drug> model`.

## Stop-and-ask triggers (consolidated)

Don't guess — ask the user when:

- The source has multiple non-hierarchical models and it's not obvious which to extract.
- Parameter values look like initial estimates rather than final.
- Covariate encoding isn't fully specified (reference category, units, transformation).
- A source column name is not in `inst/references/covariate-columns.md` (propose a new entry and confirm).
- A source column is an alias of an existing canonical name and the mapping involves value inversion or a reference-category flip.
- A parameter name deviates from the nlmixr2lib standard (propose the canonical name and confirm).
- PKNCA output disagrees with a published NCA table by more than ~20% after careful review.
- The source is paywalled and the user hasn't supplied the text.
- An erratum search is inconclusive (e.g., paywalled journal, ambiguous correction notice) — ask the user to confirm whether any corrections apply.
- The paper's final model was fit to animal-only data (see Phase 1 species-check step).
- The PMC XML / PDF on disk contains only the paper's abstract (see Phase 1 full-text-check step).
- The on-disk file's title / first-author / journal / year / drug disagrees with the task's `Paper metadata` block (see Phase 1 paper-identity-check step).
- The paper depends on an upstream popPK model that is not on disk — either the upstream paper is identifiable (queue dependency) or unidentifiable (operator decision needed). See Phase 1 upstream-popPK detection step.
- A required structural-model parameter is absent from every on-disk source (paper + supplements + regulatory review). Use the Phase 4 missing-parameter sidecar template; do not fall back to training data.
- The worktree at dispatch has a pre-existing pushed branch for this task ID (prior completed run). See Phase 2 worktree-resumption step.
- `devtools::check()` or vignette rendering produces a C-level segfault — the environment is broken; do not paper over with `--no-build-vignettes` or similar.
- The `timeout 300 rmarkdown::render()` gate exits 124 — the vignette exceeds 5 minutes; reduce simulation size before committing.
- (DDMORE-source) The DDMORE bundle directory has no `Executable_*` file (only `Output_*.lst`) — the model is output-only and cannot be re-extracted faithfully. Sidecar-ask whether to skip or chase the source `.mod` separately.
- (DDMORE-source) The `Model_Accomodations.text` is missing or content-light, the `.mod` `$PROBLEM` line and `;;` author comments don't disambiguate authorship/drug, and a PubMed lookup by `$PROBLEM` keywords returns nothing. Sidecar-ask the operator to identify the publication or flag "DDMORE repo entry only" with no linked publication.
- (DDMORE-source) The `.mod` references a covariate not present in the simulated dataset shipped with the bundle (i.e., the simulated dataset is incomplete relative to the model). Sidecar-ask whether to (A) reconstruct the covariate from the bundle's metadata, (B) skip the covariate effect, or (C) defer the task pending a complete dataset.

Use this fixed format for ambiguities:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

When running interactively, use `AskUserQuestion`. When running under
`claude_runner`, follow the runner's injected preamble instructions for
the sidecar stop-and-ask protocol (the runner provides the file paths and
schema).

## Refusal handling

If at any phase you find yourself producing a content-policy refusal in response to legitimate clinical-trial content (drug + dose + indication + adverse-event language is normal in pharmacometric papers and is not a safety concern), do not silently degrade the extraction. Sidecar-ask:

> Phase <N> step <M> hit a content-policy refusal while reading <paper section / table>. The content is standard published clinical pharmacology (drug + dose + indication + AE language). Options: (A) rephrase the extraction prompt and continue, (B) skip the affected section and document the gap in vignette Errata, (C) defer the task pending operator review.

A refusal is operator-actionable signal, not a fatal error. Treat it the same as any other stop-and-ask trigger.

## Constraints

- Never invent parameter values; if it's not in the source, ask.
- Never tune parameters to make a validation output match a target.
- This skill **only adds** new models. Retrofitting existing vignettes to use PKNCA, or renaming covariates in existing files, is a future separate skill.
- Never push directly to `main`. Open a PR.
