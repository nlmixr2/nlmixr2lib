# DDMORE Foundation Model Repository extractions

> **Temporary skill extension.** The DDMORE-source flow documented here
> is a one-time batch addition for extracting the
> `dpastoor/ddmore_scraping` bundle (~58 models) into
> `inst/modeldb/ddmore/`. **Once the batch is complete and merged, this
> file and the DDMORE-aware sections of `SKILL.md`,
> `naming-conventions.md`, `verification-checklist.md`, and
> `model-file-template.md` should be removed**, leaving the skill back
> at its paper-source-only shape. The completed
> `inst/modeldb/ddmore/` model files stay; only the skill scaffolding
> goes.
>
> Cleanup recipe (when the batch is done):
>
> 1. Delete `references/ddmore-source.md` (this file).
> 2. In `SKILL.md`: revert the frontmatter `description` to its
>    paper-only form; remove the "Two source shapes…" paragraph;
>    remove the `references/ddmore-source.md` bullet from References;
>    revert Phase 1 step 1 to its paper-source-only form; remove the
>    DDMORE-source addendum on step 9; revert step 11; remove the
>    DDMORE-source variant of Phase 3 file path + ddmore_id /
>    replicate_of; remove the DDMORE-aware Phase 5 validation gating;
>    remove the DDMORE entry from Phase 6 NEWS examples; remove the
>    three DDMORE-specific stop-and-ask triggers at the bottom.
> 3. In `naming-conventions.md`: remove the year-letter collision
>    bullet (or keep it — it's general-purpose); remove the entire
>    "NONMEM → nlmixr2 syntax translation" section.
> 4. In `verification-checklist.md`: revert the DDMORE-source
>    sentence on the final-estimates item; revert the Phase E
>    DDMORE bullet on file path; remove the two DDMORE-only Phase E
>    bullets on `ddmore_id` / `replicate_of`; revert F first item to
>    drop the "Skip when…" qualifier; remove F.2 and F.3; revert the
>    Phase H source-identity DDMORE bullet.
> 5. In `model-file-template.md`: revert the lead paragraph; remove
>    the `ddmore_id` / `replicate_of` placeholders in the template;
>    remove the "DDMORE-source metadata" section.
>
> A `git revert` of the single commit that introduced this support
> (`add-ddmore-extraction-support` branch's first commit) achieves
> almost the whole rollback in one step — re-apply by hand only the
> general-purpose year-letter collision rule in `naming-conventions.md`
> if it's worth keeping.

This reference is loaded only when the source for `extract-literature-model` is a DDMORE Foundation Model Repository bundle (NONMEM `.mod` / `.ctl` control stream + `.lst` listing + RDF metadata + `Model_Accomodations` text). For paper-source extractions, follow `SKILL.md` directly without reading this file.

The DDMORE Foundation Model Repository ([repository.ddmore.eu](https://repository.ddmore.eu); formerly `repository.ddmore.foundation`) hosts a curated set of executable pharmacometric models. Each entry is identified by a 12-character ID of the form `DDMODEL00000<NNN>`, where `<NNN>` is a zero-padded 5-digit number. The website is sometimes flaky; the canonical local source for these models is the GitHub mirror [`dpastoor/ddmore_scraping`](https://github.com/dpastoor/ddmore_scraping), which has 58 numbered directories spanning IDs 155-301 (with gaps) and includes the per-model RDF metadata files the headline `repository.ddmore.eu` API does not expose.

DDMORE-source models always land under `inst/modeldb/ddmore/`. Filename convention is the same as for `specificDrugs/`: `<FirstAuthor>_<Year>_<drug>.R`. On collision (e.g., two scenarios from the same first author / year / drug), append a lowercase letter to the year: `Themans_2019a_meropenem.R`, `Themans_2019b_meropenem.R`.

## Per-directory file inventory

A typical dpastoor directory `<id>/` (where the numeric `id` corresponds to `DDMODEL00000<id>`) contains:

| File | Purpose |
|------|---------|
| `Executable_<name>.{mod,ctl,xml,txt}` | NONMEM control stream (`.mod` / `.ctl` typical), MDL XML (rare), or NMTRAN text (`.txt`). Re-runnable. |
| `Output_real_<name>.lst` | NONMEM listing from running the model on the original real dataset. **Final parameter estimates** live in the `FINAL PARAMETER ESTIMATE` block after `MINIMIZATION SUCCESSFUL`. |
| `Output_simulated_<name>.lst` | Listing from running the model on a shipped simulated dataset (a regression-style smoke test). |
| `Output_<name>.{ext,cor,cov,phi,coi,shk,shm}` | Optional NONMEM run artifacts (some dirs only): `.ext` is the iteration log, `.cor` / `.cov` are correlation / covariance matrices, `.phi` is per-subject empirical Bayes, `.shk` / `.coi` / `.shm` are shrinkage / condition-number / mixture diagnostics. Useful for cross-checks; not directly used in the model file. |
| `*tab.simulated_<name>` | NONMEM `$TABLE` output for the simulated run. |
| `Simulated_<name>.csv` (or `.CSV`) | Simulated event dataset. The `Output_simulated_*` listing is reproduced from this — useful for the F.2 self-consistency vignette check. |
| `DDMODEL00000<id>.rdf` | RDF metadata: `model-has-description`, `model-field-purpose` URI (PK / PD / count / etc.), `model-tasks-in-scope`, `model-research-stage`, `model-modelling-question`, `model-origin-of-code-in-literature-controlled` flag. |
| `Model_Accomodations.text` (or `.txt`) | Plain-text reference citation: paper title, journal, year, authors, scenario notes, model differences from the publication. **This is the primary publication-mapping source** when the `.mod` `$PROBLEM` line is terse. |
| `Command.txt` | The NONMEM run command used (typically just `execute Executable_<name>.mod`). |
| `<id>.json` | Scraper metadata: file list, repo version (e.g., `"version": 52` for DDMODEL00000301). The version field can be cross-checked against the live DDMORE repository. |

When extracting, prefer files in this order: `Output_real_*.lst` (final estimates) → `Executable_*.mod` (structural model + initial estimates + comments) → `Model_Accomodations.text` (publication identification) → `DDMODEL00000<id>.rdf` (model classification).

## Phase 1 (DDMORE-source variant)

Replaces steps 2-5 and 7 of `SKILL.md` Phase 1; steps 6, 8, 9, 10, 11 still apply (DDMORE-aware tweaks below).

1. **Verify the bundle directory matches the task ID.** The task names a specific DDMORE folder (e.g., `.../ddmore_scraping/301/`). Confirm the directory exists, contains an `Executable_*` file, and that any `<id>.json` field matches the task ID. If the directory is missing or has no executable, sidecar-ask:

   > The DDMORE bundle directory `<path>/<id>/` either doesn't exist or contains no `Executable_*` file (only `Output_*.lst` listings). The model is not re-extractable from this bundle. Options: (A) skip this task, (B) chase the source `.mod` from the live DDMORE repository (the operator handles), (C) defer pending operator investigation. Which applies?

2. **Identify the publication from `Model_Accomodations.text|.txt`.** Read the file (it is short — typically under 30 lines). Parse for:
   - **Title** of the linked publication (often the first non-`##` line).
   - **Authors** (a line starting with `## Authors:` or similar).
   - **Journal + year** (look for "Br J Clin Pharmacol", "submitted - July 2019", etc.).
   - **Model differences from publication** (lines after the `######` separator; flag any non-trivial deviation in the vignette's Errata section).

   If `Model_Accomodations.text|.txt` is missing or content-light:

   - Fall back to the `.mod` header: read `head -40 Executable_*.mod` and look for `;;` author comments (e.g., `;; Pauline Thémans, Joseph J. Winkin, Flora T. Musuamba`) and the `$PROBLEM` line.
   - Try a PubMed E-utilities search by `<first-author> <year> <drug>` keywords; if a match returns within the first three results, use it.
   - If still ambiguous, sidecar-ask:

     > The DDMORE bundle for `DDMODEL00000<id>` does not have a `Model_Accomodations` file or sufficient `.mod` header metadata to identify the publication. PubMed search by `<keywords>` returns: `<top 3 hits>`. Options: (A) use hit #1, (B) use a different hit (specify), (C) extract with `reference = "DDMORE repo entry only; no linked publication"` and skip the publication-comparison vignette section. Which applies?

3. **Classify the model type from RDF + structural inspection.** Read `DDMODEL00000<id>.rdf` and look for `j.2:model-field-purpose rdf:resource="..."`. The URI fragment encodes the modality (e.g., `pkpd_0001024` ≡ pharmacokinetics). Cross-check by inspecting the `.mod`:
   - `$SUBROUTINE ADVAN1/2/3/4/11/12` → linear PK; `linCmt()`-compatible.
   - `$SUBROUTINE ADVAN5/6/8/9/13` → ODE-based; will need explicit `d/dt(...)` translation.
   - `$PRED` only (no `$SUBROUTINE`) → typically count, Markov, IRT, dropout, or TTE; the validation strategy is mechanistic-sanity (verification-checklist.md § F.3), **not** PKNCA.
   - Mixed `$SUBROUTINE ADVAN... + $DES + biomarker covariate equations` → PD or PK/PD; verification strategy depends on what the publication reports.

   The model type drives the vignette validation strategy (PKNCA vs mechanistic sanity vs F.2 self-consistency).

4. **Detect `replicate_of` candidates** by grepping `inst/modeldb/specificDrugs/*.R` for the same first-author + year + drug combination. If a paper-derived counterpart already exists, plan to add reciprocal `replicate_of` pointers in both files. (See `model-file-template.md` § "DDMORE-source metadata" for the schema.) The same task should update the counterpart in `specificDrugs/` — either in the same PR or, if the counterpart's PR is already merged, in a follow-up edit (flag in PR body).

5. **Multiple-scenario handling.** Some DDMORE entries ship two `Executable_*` files representing different scenarios (e.g., `280/` has `Executable_real_*` and `Executable_simulated_*` with materially different structures). Treat each as a separate extraction (`<author>_<year>a_<drug>.R`, `<author>_<year>b_<drug>.R`) unless the scenarios are non-mechanistic differences (e.g., dataset-only changes); in that case extract once and document the scenario range in the vignette's Errata.

## Reading final estimates from `.lst`

NONMEM listing files have a fixed structure. The key block is `FINAL PARAMETER ESTIMATE` after `MINIMIZATION SUCCESSFUL` (or a similar success marker — `MINIMIZATION TERMINATED` indicates non-convergence and is a stop-and-ask trigger). The block looks like:

```
 #OBJV:********************************************    -1234.567       **************************************************
 ************************************************************************************************************************
 ************************************************                      ************************************************
 ************************************************      FINAL PARAMETER ESTIMATE           ************************************************
 ************************************************                      ************************************************
 ************************************************************************************************************************


 THETA - VECTOR OF FIXED EFFECTS PARAMETERS   ********


         TH 1      TH 2      TH 3      ...
         7.94E+00  7.22E-01  1.36E+01  ...

 OMEGA - COV MATRIX FOR RANDOM EFFECTS - ETAS  ********


            ETA1      ETA2      ETA3
 ETA1
 +        1.26E-01
 ETA2
 +        0.00E+00  1.40E-01
 ETA3
 +        0.00E+00  0.00E+00  1.76E+00
 ...
```

Translate by mapping each `TH <i>` to its labelled `THETA(<i>)` slot in the `.mod` (the `;` comments after each `$THETA` line in the executable identify what each slot is — `; CL`, `; GFRCL`, etc.). The `OMEGA` block diagonals are the IIV variances; off-diagonals are covariances (for `BLOCK` `$OMEGA` declarations).

`SIGMA - COV MATRIX FOR RANDOM EFFECTS - EPSILONS` follows the same structure for residual error.

Some listings also have a `STANDARD ERROR OF ESTIMATE` block after the final estimates — useful for confidence intervals but not used in the model file (nlmixr2 doesn't store them).

If the listing reports `MINIMIZATION TERMINATED` or `R MATRIX ALGORITHMICALLY SINGULAR`, sidecar-ask:

> DDMODEL00000<id>'s `.lst` reports `<status>` rather than `MINIMIZATION SUCCESSFUL`. Final estimates may be unreliable. Options: (A) extract the reported values with a clear caveat in the vignette Errata, (B) skip this task, (C) defer pending operator investigation. Which applies?

## Validation strategy by model type

Decision tree for the vignette validation section:

1. Is there a linked publication on disk (or accessible via DOI)?
   - **No** → use F.2 self-consistency check (re-simulate vs `Output_simulated_*.lst`). Skip publication comparison.
   - **Yes** → continue.
2. Does the publication report PKNCA-amenable PK/PD with NCA values (Cmax, AUC, t½)?
   - **Yes** → standard PKNCA vignette section (see `pknca-recipes.md`).
   - **No** → continue.
3. Is the model count / Markov / IRT / dropout / TTE?
   - **Yes** → use F.3 mechanistic-sanity recipes (verification-checklist.md § F.3): typical-value simulation reproduces published expected count / hazard / probability trajectory.
   - **No** → continue.
4. Is the model endogenous / mechanistic / steady-state turnover?
   - **Yes** → use endogenous-validation.md recipes.
   - **No** → sidecar-ask the operator for a custom validation strategy. Do not commit a vignette without one.

Vignette template stays the same (`vignette-template.md`); only the validation section's body changes.

## When the live DDMORE repository is needed

The dpastoor mirror lags the live repository by an unknown amount. If the operator suspects the local bundle is out of date — or wants a model not in dpastoor — the live URL pattern is:

```
https://repository.ddmore.eu/model/DDMODEL00000<id>#Files
```

Files are downloaded via:

```
http://repository.ddmore.eu/model/download/DDMODEL00000<id>.<version>?filename=<name>
```

Where `<version>` is the latest version listed under the model's "History" tab on the web UI. The `<id>.json` `version` field in the dpastoor bundle is the version that was current when dpastoor scraped; if the live UI shows a higher version, the bundle is stale.

This skill does NOT scrape the live repository directly — that path is operator-handled. If the bundle is stale, sidecar-ask:

> The dpastoor bundle for DDMODEL00000<id> reports version <X>; the live DDMORE repository web UI shows version <Y>. Options: (A) extract from the dpastoor bundle anyway and note the version skew in the vignette Errata, (B) defer this task pending an operator-pulled refresh of the bundle, (C) skip. Which applies?

## Common DDMORE pitfalls

1. **`.mod` $THETA / $OMEGA / $SIGMA are initial values.** Always read the `Output_real_*.lst` for final estimates (see "Reading final estimates from `.lst`" above). The verification-checklist § A item is binding.

2. **`$PROBLEM` lines are sometimes terse.** Many DDMORE entries have `$PROBLEM   model_run` or similar minimal text. Do not rely on `$PROBLEM` as the primary metadata source; always cross-check against `Model_Accomodations.text|.txt`.

3. **Authors may be in `;;` comments at the very top of the `.mod`.** Conventional NONMEM uses `;` for comments and `;;` for double-emphasis. DDMORE conventions sometimes use `;;` for author / metadata lines. Read the first 5-10 lines of `Executable_*.mod` even when `Model_Accomodations.text` is present.

4. **Multi-output models with conditional `$ERROR`.** When the `.mod` has `IF (CMT.EQ.1) Y = ...` / `IF (CMT.EQ.2) Y = ...`, the nlmixr2 translation produces one residual-error line per output (`Cc ~ prop(...)` and `<other_output> ~ ...`), with output-suffixed residual-error parameter names. See `naming-conventions.md` § "Observation variable" and § "Residual error."

5. **MDL XML format (a small minority of entries — e.g., DDMODEL00000213).** MDL is a different language from NONMEM; the auto-translation tooling embedded in DDMORE renders MDL → NMTRAN at submission, but the raw `.mod` is sometimes the rendered NMTRAN form, sometimes the original MDL XML. If `Executable_*.xml` is the only executable, sidecar-ask the operator before attempting a translation.

6. **Bundle-shipped simulated datasets are intentionally minimal.** A `Simulated_<name>.csv` may contain a single dose, one subject, one observation per visit — enough to verify the `.lst` reproduces but **not** representative of the publication's clinical study. Do **not** treat it as a representative population for the vignette virtual cohort; build a virtual cohort from the publication's reported demographics (Phase 5, vignette-template.md).

## Worked example

The reference worked example (committed alongside this skill update) is `inst/modeldb/ddmore/Themans_2019_meropenem.R` — a 3-compartment meropenem PK model with covariates on CL (GFR) and V1, V2 (WT). Its source bundle is `.../ddmore_scraping/301/`. Read both side-by-side when starting your first DDMORE-source extraction.
