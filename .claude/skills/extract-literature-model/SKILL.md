---
name: extract-literature-model
description: This skill should be used when the user wants to "add a model from a paper", "extract a pharmacometric model from the literature", "implement a published PK/PD model in nlmixr2lib", or provides a scientific article, conference poster, supplement, or regulatory review document and asks to add the model to nlmixr2lib. Guides source review, standardized model file creation under inst/modeldb/, in-file source-trace verification, and validation vignette with PKNCA NCA checks.
---

# Extract a pharmacometric model from the literature

Input: a scientific source describing a pharmacometric model (journal article, supplement, conference poster, or regulatory review).
Output: a packaged nlmixr2lib model file under `inst/modeldb/`, a validation vignette under `vignettes/articles/`, and updated registry artifacts — opened as a pull request against `main`.

Work through the six phases below. Stop and ask the user at any of the decision points called out explicitly; ambiguity is the main failure mode for this workflow, and silent assumptions are what get shipped as bugs.

## References

Read on demand based on the paper's model class. The most-used references are listed first.

**Always:**
- `references/pre-flight-checklist.md` — consolidated stop-and-ask triggers. Read this once at dispatch.
- `references/compartment-names.md` — compartment and observation-variable names.
- `references/parameter-names.md` — structural PK, IIV, residual-error, covariate-effect, fixed-parameter encoding, file-level metadata, file naming.
- `inst/references/covariate-columns.md` — authoritative register of covariate column names. Consult before introducing any new covariate. (Installed with the package so R code like `checkModelConventions()` can parse it at runtime.)
- `references/model-file-template.md` — skeleton for the `.R` file.
- `references/vignette-template.md` — skeleton for the validation vignette.
- `references/verification-checklist.md` — checklist to walk after the first-pass implementation.

**Conditional (load only when the model class applies):**
- `references/nonmem-translation.md` — NONMEM → nlmixr2 syntax conversion. Load only when the source is a NONMEM control stream (`.mod` / `.ctl`, supplement with `$PK` / `$DES` / `$THETA` blocks).
- `references/pknca-recipes.md` — PKNCA setups for single-dose, steady-state, and multi-dose NCA. Load for popPK validation; skip for endogenous / mechanistic models.
- `references/endogenous-validation.md` — validation strategy for endogenous / mechanistic / turnover models. Load only when the source describes turnover, enzyme-kinetic, or steady-state-balance models where PKNCA is not the right check.
- `references/oa-acquisition.md` — detailed OA-PDF ladder. The Phase 1 Step 0 script (`scripts/acquire-paper.R`) is the entry point; this reference is for understanding the script's failure modes or extending the source list.

**Scripts:**
- `scripts/acquire-paper.R` — open-access PDF acquisition (Phase 1 Step 0).
- `scripts/lint-conventions.R` — pretty-print `checkModelConventions()` results (Phase 6 step 3).

## Phase 1 — Source acquisition and scoping

### Step 0 — Source-file presence and self-acquisition

Confirm every source file the task references is on disk and is a valid PDF. The task block includes `Lead PDF`, `Supplements`, `Model files`, and `Source dir` paths; check each path that is set.

For each missing or invalid file, run the acquisition script:

```bash
Rscript .claude/skills/extract-literature-model/scripts/acquire-paper.R \
  --doi "<DOI>" \
  --out "<expected-on-disk-path.pdf>" \
  --title "<expected-first-author or short title>" \
  --log "/tmp/acquire-log.json"
```

Interpret the exit code:

- **0** → PDF at `--out` is valid; proceed.
- **2** → all 5 sources tried, none produced a valid PDF whose title matched. Read `/tmp/acquire-log.json` for the per-source outcome, then sidecar-ask the operator with the structured attempts list (Options: A operator drops the PDF and re-dispatches, B operator emails corresponding author, C skip). Do **not** attempt to fabricate the missing source.
- **3** → environment problem (missing `curl`, bad DOI). Fix the environment; do not work around.

For supplements / errata that the script doesn't cover, see `references/oa-acquisition.md` for the full ladder; if a supplement is unobtainable but the lead PDF is on disk, decide based on whether the missing supplement contains parameter values you need (sidecar-ask if it does; document the gap in vignette Errata if it doesn't).

### Numbered steps

1. Confirm the source type (journal article, supplement, poster, regulatory document).
2. **Verify the on-disk file is the paper the task names.** Open the source file (or its trimmed `.md` companion) and read the title + first-author line + journal + year. Compare against the task's `Paper metadata` block. If any of {first-author, year, journal, drug} disagrees, stop and sidecar-ask:

   > The on-disk file `<filename>` reports `<actual-author> <actual-year>, <actual-journal>: "<actual-title>"` but the task names `<expected-author> <expected-year>, ... <expected-drug>`. Options: (A) the task metadata is wrong and the file is the right paper — please correct the task metadata, (B) the file is mislabelled — please drop the correct PDF, (C) skip this task. Which applies?

   This catches both "wrong-PMID-in-filename" (e.g. Frey 2010 saved as `PMID_23436260.pdf`) and "wrong-drug-guessed" task-generator errors. Do not proceed with extraction until the mismatch is resolved.

3. **Record the study population species.** Skim the Methods / Subjects section to identify which species the final model was fit to (human, rat, mouse, cyno, dog, beagle, sheep, swine, in-vitro cell line, etc.). Preclinical and in-vitro models are first-class members of nlmixr2lib — extract them the same way you would a human popPK model, with the species recorded explicitly so downstream users can filter by it.

   - **Required**: set `population$species` to a human-readable label (e.g., `"rat (Sprague-Dawley)"`, `"mouse (HBCx-9 PDX)"`, `"beagle dog"`, `"in vitro (SKBR3)"`, `"human"`). For pooled human+animal cohorts, list each species (e.g., `"human + rat"`).
   - **Required for non-human**: prepend the species label to `description` (e.g., `description = "Preclinical (rat). Two-compartment popPK ..."` or `description = "In vitro (SKBR3 cell line). Mechanistic HER2 trafficking ..."`) so the species is visible in `modellib()` listings without inspecting `population`.
   - **Recommended**: if the paper provides a human-scaled projection alongside the preclinical fit (e.g., allometric forward projection for first-in-human dosing), extract the preclinical and the human-scaled projection as **two separate model files** (`Author_Year_drug_<species>.R` and `Author_Year_drug_human.R`) so each is self-contained.
   - **No sidecar needed for the species choice itself** — extract preclinical models without asking. Sidecar-ask only when the paper has *other* triggers (missing parameters, ambiguous covariate encoding, mixture model, multiple non-hierarchical models, etc.).

   The drug field name should still match the on-disk PDF (per step 2). Animal cohorts often use the same INN as the human drug; in that case, the file is `Author_Year_drug_rat.R` (or similar species suffix) so it doesn't collide with a sibling human extraction.

### Step 3a — Model class: PBPK, QSP, MBMA, mechanistic-systems

Whole-body PBPK, semi-mechanistic PBPK, quantitative systems pharmacology (QSP), and model-based meta-analysis (MBMA) papers are first-class members of nlmixr2lib when the source provides enough information to reproduce the model end-to-end. Extract them the same way you would a popPK model, with the following discipline points specific to non-popPK model classes:

- **`description` prefix.** Prepend the modeling framework to `description` so users filtering by mechanism class see it in `modellib()` listings. Examples:
  - `"PBPK (whole-body, SimCYP V12.1). APAP disposition in pregnancy ..."`
  - `"QSP. 5-lipoxygenase pathway with 33 ODEs, 64 rate laws, 113 parameters ..."`
  - `"MBMA. Trastuzumab + T-DXd from 103 published trials; between-study variability only ..."`

- **Filename suffix when a popPK extraction already exists for the same drug.** Use `Author_Year_<drug>_pbpk.R`, `_qsp.R`, or `_mbma.R` so the mechanism class is visible at file-tree scan and so they don't collide with a sibling compartmental popPK extraction of the same paper / drug.

- **Parameter sourcing — NEVER substitute from another paper unless that paper is on disk.** PBPK and QSP models often have 50-200+ parameters, many of which live in supplements rather than the main text. If a parameter is not reported in any on-disk source (paper text, on-disk supplement, on-disk upstream citation), DO NOT fill in the gap from training-data knowledge, class-typical values, default SimCYP / GastroPlus / Open Systems Pharmacology library entries, or "what a similar PBPK paper used for the same drug class." Instead, fire the Phase 4 missing-parameter sidecar and let the operator decide between (A) author correspondence, (B) operator drops the supplement, or (C) skip the paper.

  This rule is **stricter** for PBPK / QSP / MBMA than for popPK because the literature is full of mechanism papers that cite each other circularly — substituting from "a representative PBPK paper for the same drug" introduces silent provenance ambiguity (which paper's parameters? which body-composition assumptions? which permeability model? which between-study covariance prior?) that a downstream user cannot audit. The acceptable substitution path is exactly one: the current paper explicitly states "we used the <Author Year> body-composition / metabolism / permeability model" AND that paper is on disk in the source directory; the inherited parameters then get an inline `# inherited from <Author Year> Table N` comment pointing to the upstream paper.

  Class-typical training-data defaults are NOT an acceptable substitute. "kdeg for trastuzumab is typically 0.5/day in the literature" is not a citation; it's a guess. Sidecar instead.

- **Reproducibility check before drafting.** Walk the paper for these provenance markers:
  - Are all structural parameters (Vmax, Km, kcat, organ volumes, blood flows, permeabilities, partition coefficients, between-study covariance, etc.) tabulated explicitly with units?
  - Are the ODEs written out, or only described in prose ("we used the standard SimCYP whole-body model")?
  - Are the dosing-event handling and observation-mapping rules unambiguous?
  - For MBMA: are the per-study weights / variance estimates explicit, or only summary parameter point estimates?

  If any of these are insufficient, sidecar-ask the operator before drafting. Do not "fill in the gaps" from PBPK / QSP / MBMA class knowledge.

- **PBPK / QSP platforms vs nlmixr2.** SimCYP, GastroPlus, PK-Sim/MoBi, OSP, and similar commercial platforms come with built-in whole-body models that users parameterise. nlmixr2 does not have those platforms — when extracting a SimCYP-published model, you must encode the platform's whole-body ODEs explicitly in `model()`. If the paper only states "we used SimCYP's standard whole-body model" without writing out the ODEs and parameter values, the model is not reproducible from on-disk sources and the appropriate action is sidecar-ask or skip — not fill in from training-data knowledge of how SimCYP's standard model is typically parameterised.

- **MBMA `between-study variability` ≠ popPK BSV.** MBMA papers typically report between-study variance (a study-level random effect on summary metrics) rather than between-subject variance (individual-level $\eta$ on $\theta$). When extracting an MBMA model:
  - Encode between-study variance as a study-level eta (e.g., `eta_study_lcl ~ <var>`) clearly labelled to distinguish from the popPK pattern; do NOT silently relabel it as `etalcl` (between-subject).
  - The model is suitable for simulating study-level summary outcomes (per-arm mean response) but NOT individual-level concentrations.
  - Document the simulation scope prominently in `description` and in the vignette's Assumptions and deviations section.

4. **Prefer trimmed markdown when available.** The preprocessor at `mab_human_consensus/tracking/preprocess_papers.py` writes a `<stem>_trimmed.md` next to each source file (PMC XML, PDF, DOCX, XLSX) containing only the sections the extraction actually needs: Title + Abstract + Methods + Results + Tables + Figure captions. The Introduction, Discussion, Conclusions, References, Acknowledgments, and publisher boilerplate are stripped. If `PMID_<pmid>_pmc_trimmed.md` (or `PMID_<pmid>_trimmed.md` for a PDF, or `<stem>_trimmed.md` for a supplement) exists, read it **instead of** the raw `.xml` / `.pdf` / `.docx` — it's typically 40-95% smaller with no loss of extractable content. Full-text sanity check on the trimmed file: ~15 KB+ (full-text trim) vs < 3 KB (abstract-only trim). Fall back to the raw source only if the `_trimmed.md` doesn't exist, the trim appears to have lost a specific piece of information you need (rare — only when the paper is structurally unusual), or you explicitly need the discarded sections (e.g., to quote a Discussion claim in the vignette narrative).

5. **Verify the source contains full text, not just the abstract.** Wiley / BJCP and some other publishers serve PMC XML containing only front matter + abstract. Before reading for model structure, run a quick sanity check (against the trimmed file if present, otherwise the raw source):

   - Trimmed `.md` file size ≥ ~15 KB, or raw PMC XML ≥ ~40 KB (full-text XML is typically 100 KB+; abstract-only is usually < 20 KB).
   - The file contains a materially-present Methods section (not just a "Methods" heading followed by one abstract paragraph).
   - If only a PDF is on disk, confirm it runs past the abstract (multi-page, Methods / Results / Tables present).

   If only the abstract is available, sidecar-ask:

   > The source on disk for <paper> contains only the abstract and front matter; full text appears to be blocked by the publisher. Options: (A) pause this task until full text is provided, (B) proceed only if a supplement / regulatory review on disk contains the model equations and parameter tables, (C) skip this paper. Which applies?

   Never attempt extraction from an abstract alone — population-PK parameter values, covariate effects, and equations are not in an abstract.

6. **Detect upstream-popPK dependencies — DEFAULT is auto-handle, not sidecar.** Skim Methods for phrases like "PK was described using the popPK model previously developed from <study/phase>", "the structural PK model was fixed from <reference>", "covariate effects were carried over from <author> et al.", or "the PK model from <prior publication> was used as a backbone." If the current paper's PD model fixes its PK parameters from a separate publication that is not on disk:

   1. **Try to identify** the upstream paper from the references list (PMID, DOI, or full citation).

   2. **If identifiable, run the helper script** to handle the dependency automatically (no sidecar):

      ```sh
      python3 .claude/skills/extract-literature-model/scripts/handle-upstream-dependency.py \\
          --queue-dir <queue path>             `# the runner queue dir` \\
          --current-task-id "$TASK_ID"         `# the downstream task` \\
          --upstream-pmid <PMID>               `# or --upstream-doi <DOI>` \\
          --upstream-citation '<verbatim citation>' \\
          --upstream-drug <drug name>          `# e.g. clozapine, warfarin`
      ```

      The script:
        a. **Tries to acquire** the upstream PDF via `acquire-paper.R`'s OA-PDF ladder.
        b. **If acquired**: places the PDF at `papers/PMID_<id>/PMID_<id>.pdf` and drops a `trim_queue` marker so the trim daemon picks it up.
        c. **If NOT acquired** (paywall, no OA): writes `papers/PMID_<id>/PMID_<id>_needs_acquisition.flag` with clear "operator drops the PDF here" instructions.
        d. **Either way**: queues the upstream as a new task in `todo/` (next free `NNN-<author>_<year>_<drug>.yaml` slot), with a prompt that re-uses this same `/extract-literature-model` skill.
        e. **Edits the current task's YAML** to add `depends_on: [<upstream_task_id>]`.
        f. **Writes a deferral report** and exits cleanly.

      The current task does NOT inline-extract the upstream model and does NOT extract the current paper's PD layer in isolation. Wait for the upstream extraction to commit; on re-dispatch, the current task imports the upstream PK parameters with `# inherited from <Author Year> Table N` provenance comments.

   3. **If unidentifiable** (e.g. "the popPK model from internal Phase 1/2 studies" with no specific citation, an unpublished company internal model, or an in-house simulator output), sidecar-ask: (A) skip the task, (B) proceed with parameters fixed inline as reported in the current paper, with a clear "upstream PK source not located" note in the vignette Errata, or (C) defer pending operator investigation. The auto-handle helper above is only safe when the upstream is uniquely identifiable.

   Never silently fabricate upstream PK parameters from training data.

7. **Always search for supplementary information.** Supplements frequently contain the NONMEM control stream and parameter tables that disambiguate model structure. If the user provided only a main article, ask whether a supplement exists and request it.
8. **Always search for errata, corrigenda, or author corrections.** Check the journal's landing page for the article, the publisher's "corrections" / "notices" feed, and a search like `"<first author> <year> <drug>" erratum` on PubMed and Google Scholar. Ask the user whether they are aware of any corrections if the source is paywalled or the search is inconclusive. **When an erratum revises a value used in the model (parameter estimate, covariate effect, equation, units), the erratum value takes precedence over the main publication.** If multiple errata exist, the most recent supersedes earlier ones. Record the erratum citation in the model file's `reference` metadata alongside the main paper, and in every in-file source-trace comment whose value comes from the erratum, point to the erratum (not the original table).
9. **Verify parameters are final estimates, not initial estimates.** Supplement control streams usually list initial values in `$THETA` / `$OMEGA`; final values come from the paper's results table or `$TABLE` output. If only a control stream is available, confirm values against any published point estimates before treating them as final.
10. **Multiple-model handling.**
    - Base model + final model → extract only the final.
    - Any other "multiple model" case (per-subpopulation, per-endpoint, sensitivity analyses) → list the candidates to the user and ask which to extract. Offer "one," "all," or "a subset."
11. **Systematic review / meta-analysis handling.** If the source is a systematic review or meta-analysis that catalogs other authors' published popPK / PD models without developing an original model of its own, **the default action is to skip the task and queue the cited primary papers for future extraction.** Do not extract a cataloged model from the review's summary table — extracting from a secondary source loses model-selection rationale, covariate-encoding details (reference categories, units, allometric exponents, scaling normalisations), and the full residual-error / IIV structure that only the primary papers contain. Recognise systematic reviews by:

    - An explicit "Review" / "Systematic Review" tag on page 1 or in the article-type header.
    - A Methods section describing a literature-search protocol (PubMed / EMBASE / Scopus / Web of Science search query, PRISMA-style screening flowchart, inclusion / exclusion criteria).
    - A Results section that tabulates other authors' models in side-by-side tables (one row per cited study, columns for structural model / parameter estimates / covariates).
    - No `$THETA` / `$OMEGA` / `$SIGMA` block, no original VPC / GOF figures, no original NONMEM control stream attributable to the review's authors.

    Sidecar-ask the operator with the list of cited primary papers and three options:

    > `<paper-name>` is a systematic review of <N> previously published popPK / PD models. Per the standing policy, the recommended action is to skip this task and queue the cited primary papers individually for future extraction. Cited primary models: <numbered list with first-author + year + drug + journal>. Confirm: (A) skip this task and queue all <N> cited primary papers for future extraction (operator-followups register), (B) skip without queueing references (e.g. the references duplicate already-queued tasks), (C) extract one or more models from the review's tables (operator names which; the review becomes the transcription source with provenance noted prominently in vignette Errata).

    Option (A) is the recommended default. The cited-papers list is committed in the report so the operator can add them to the queue in batch.

12. Confirm the target subdirectory under `inst/modeldb/` (usually `specificDrugs/`; endogenous, therapeuticArea, pharmacokinetics, and pharmacodynamics are also valid).

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

File path: `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`.

If the chosen `<FirstAuthor>_<Year>_<drug>` name collides with an existing file (rare — e.g., two same-author/year/drug entries with different scenarios), append a lowercase letter to the year: `Author_2019a_drug.R`, `Author_2019b_drug.R`. Use the same year-letter for both files in the pair so the chronological ordering is preserved.

The function name **must equal the filename** minus `.R`; `buildModelDb()` enforces this.

Use `references/model-file-template.md` as the starting skeleton and the two best-formed existing models as anchors:

- `inst/modeldb/specificDrugs/Clegg_2024_nirsevimab.R` — covariates, maturation, correlated IIV, exported race-derivation helper.
- `inst/modeldb/specificDrugs/Hu_2026_clesrovimab.R` — simpler case with Hill-type maturation.

The file body has this shape:

1. `description`, `reference`, `vignette`, `units`, `covariateData`, `population` — metadata before `ini()`. `vignette` is the basename of the validation vignette in `vignettes/articles/` (e.g., `"Clegg_2024_nirsevimab"`, no path, no extension); `buildModelDb()` extracts it so the list-of-models table can link to the rendered vignette on the pkgdown site.
2. `ini()` — parameters with `label()` and a trailing **in-file comment pointing to the source location** for every value. **Wrap fixed parameters in `fixed()`** — see the "Fixed parameters" subsection below.
3. `model()` — derived terms → individual parameters → micro-constants → ODEs → bioavailability → observation and error.

### Fixed parameters in `ini()`

When the source paper holds a parameter at a known value rather than estimating it, encode that fact by wrapping the value in `fixed()` in `ini()`. This applies to **every** parameter type — structural THETAs, allometric exponents, IIV variances/covariances, residual-error magnitudes, covariate-effect coefficients — not just IIV. Failing to encode the fixed status loses load-bearing provenance: a downstream user cannot tell whether the value was estimated and reported as a point estimate, or whether the source authors held it constant; re-fitting the model under one assumption vs the other gives different results.

See `references/parameter-names.md` § "Fixed parameters" for the full list of source signals (prose, NONMEM `FIX` flag, etc.) and the encoding examples. The key syntactic point is that `log()` goes **inside** `fixed()` for log-transformed parameters (`lcl <- fixed(log(2))`, never `log(fixed(2))`).

When in doubt — a $THETA reported with an RSE of 0% but no FIX flag, or a parameter reported to three decimal places with no uncertainty — sidecar-ask the operator before guessing.

Follow `references/compartment-names.md` and `references/parameter-names.md` strictly. The quick-reference rules:

- Structural PK parameters log-transformed: `lka`, `lcl`, `lvc`, `lvp`, `lvp2`, `lq`, `lq2`, `lfdepot`.
- IIV: `eta` + transformed name, e.g., `etalcl` (not `etacl`). Block correlations via `etalcl + etalvc ~ c(var, cov, var)`.
- Residual error: `propSd`, `addSd`, `expSd`. Multi-output: `propSd_<output>`, `addSd_<output>`, etc.
- Compartments: `depot`, `central`, `peripheral1`, `peripheral2`, `effect`, `target`, `complex`, `csf`, `isf`, etc. Observation: `Cc`.

Covariate columns come from `inst/references/covariate-columns.md`. Before writing any covariate into the file:

- If the canonical name exists, use it and record the source column name in `covariateData[[name]]$source_name`.
- If the source name is an alias of an existing canonical name (e.g., source uses `SEXM`, canonical is `SEXF`), use the canonical name, note the required value transformation (`SEXF = 1 - SEXM`), and **ask the user to confirm the effect-coefficient sign and reference-category implications** before committing.
- If the concept isn't in the register at all, propose a new entry (canonical name, description, units, type, reference category, source aliases) and ask the user to confirm before adding it. The new entry is committed alongside the model.
- Do **not** add a change-log / history / "## Change log" / "## Summary" section or per-extraction history line to `inst/references/covariate-columns.md`. The H3 entry itself is the authoritative record; chronological history of when an entry was added or modified is read from `git log`. Any context that someone needs to use the covariate (derivation rules, scope-promotion rationale, name-collision history, naming-decision sidecars) belongs in the H3 entry's Description / Notes / Source aliases — not in a separate change log.

`population` uses the extensible schema documented in `references/parameter-names.md` § "File-level metadata". Common fields: `n_subjects`, `n_studies`, `age_range`, `weight_range`, `sex_female_pct`, `race_ethnicity`, `disease_state`, `dose_range`, `regions`, `notes`. Add any additional keys the paper describes (e.g., `ga_range`, `renal_function`, `co_medication`) — do not force facts into the common schema.

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

File path: `vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd`. Drug-specific vignettes live under `vignettes/articles/` so pkgdown builds them for the site but CRAN does not — `.Rbuildignore` excludes that directory. The basename (without `.Rmd`) must match the `vignette <- "..."` field in the model file. Only the legacy `PK_2cmt_mAb_Davda_2014.Rmd` remains at top-level `vignettes/`.

The YAML `title:` and `\VignetteIndexEntry{...}` use the **human form** `<Drug> (<FirstAuthor> <Year>)` (e.g. `Ustekinumab (Aguiar 2021)`), matching the navbar dropdown label generated by `R/modeldb.R::.parseModelLabel()`. The filename stays in the machine form. See `references/vignette-template.md` for the full header and rules around longer descriptive titles.

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

Naming conventions for mechanistic parameters are documented in `references/parameter-names.md` § "Endogenous / mechanistic parameters."

## Phase 6 — Registration, tests, docs, PR

1. Re-confirm the branch is on top of fresh `origin/main` (`git fetch origin && git rebase origin/main` if needed).
2. Run `nlmixr2lib::buildModelDb()` to regenerate `data/modeldb.rda` and `inst/modeldb.qs2`. Confirm the new model appears in `modellib()`. **When verifying in R, do `devtools::load_all(".")` first so `modellib()` reads the worktree's in-development package, not the stale system install** — see `references/verification-checklist.md` § "Verifying against the worktree's nlmixr2lib" for why a bare `library(nlmixr2lib)` can return a misleading `FALSE`.
3. Run the convention lint and review the output:

   ```bash
   Rscript .claude/skills/extract-literature-model/scripts/lint-conventions.R \
     "<FirstAuthor>_<Year>_<drug>"
   ```

   The script wraps `nlmixr2lib::checkModelConventions()` and emits PR-paste-ready output. Exit 0 = clean; exit 1 = warnings only; exit 2 = errors. Any deviations from the canonical parameter / IIV / residual-error / covariate / compartment conventions (see `references/compartment-names.md`, `references/parameter-names.md`, and `inst/references/covariate-columns.md`) should be either fixed in the model file before committing, or explicitly justified in the vignette's Assumptions and deviations section. `buildModelDb()` runs `checkModelConventions()` implicitly at package-build time, but running the lint explicitly on your new model surfaces drift before commit. Paste the lint output into the PR body so a reviewer can see what was checked.
4. **Render the vignette locally and verify it completes without error in under 5 minutes of wall-clock time. This is a mandatory pre-push gate — the worker MUST run this command and verify exit 0 / HTML present before any `git commit` / `git push` of the new vignette.** It catches the same failure modes pkgdown CI does — missing data columns, time-varying covariates assigned to rxEt objects that get silently dropped, PKNCA formulas referencing absent columns, simulation crashes — and surfaces them in seconds rather than after a CI cycle.

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

   - **Exit 0, output HTML present** → vignette is clean; the gate is passed; proceed to commit.
   - **Exit non-zero with an R error** → fix the vignette before committing. Read the traceback carefully; the most common causes are (a) a `select()` / PKNCA formula referencing a column that was never created (add it in the preceding `mutate()`), (b) a variable name used before it is defined (reorder chunks), (c) a simulation that errors because covariate values are out of range for the model, (d) `event_table$col <- vector` syntax against an `rxEt` object — rxode2 silently drops these assignments, so always materialize via `as.data.frame()` and add covariate columns there before passing to `rxSolve()`.
   - **Exit 124 (timeout)** → the vignette exceeds 5 minutes. Reduce simulation size: cut `nSub` (stochastic VPC subjects), shorten the observation grid, or move expensive chunks behind `eval = FALSE` with a note. Do **not** skip the time budget — pkgdown CI has strict wall-time limits and a slow vignette breaks the build for everyone.
   - **C-level segfault (`*** caught segfault ***`)** → broken R / rxode2 / nlmixr2 install, not a model-file problem. Stop, sidecar-ask the operator to fix the environment; do not work around with `--no-build-vignettes`.

   Do not interpret silence as success — re-run with `quiet = FALSE` and confirm an HTML output file landed next to the `.Rmd`. Do not assume "checkModelConventions was clean, vignette must be fine" — those are independent failure modes. The render gate is the only one that exercises the full data path the way pkgdown CI does.

5. **Verify no non-ASCII characters in the new model file or vignette (mandatory pre-push gate).** R CMD check warns on non-ASCII strings in package data — the `description <-` field of every model file is reified into `data/modeldb.rda`, so a single em-dash, en-dash, multiplication sign, or Greek letter in the description triggers a `WARNING: found non-ASCII strings` on every platform. Apply the gate to both the `.R` file and the `.Rmd` vignette.

   ```bash
   for f in inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R \
            vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd; do
     if LC_ALL=C grep -nP "[^\x00-\x7F]" "$f"; then
       echo "FAIL: non-ASCII in $f"; exit 1
     fi
   done
   echo "OK: all checked files are ASCII-only"
   ```

   See `references/verification-checklist.md` for the full ASCII substitution table when anything matches.

6. Run `devtools::check()`. Vignettes must build cleanly.

   A C-level segfault (`*** caught segfault ***`) during `check()` or vignette rendering is a red flag — it indicates a broken R / rxode2 / nlmixr2 install in the environment, not a model-file problem. Stop, sidecar-ask the operator to investigate and fix the environment, and do not work around it with `--no-build-vignettes` or similar flags.
7. Add a short, single-line `NEWS.md` entry under the current development
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

   Do NOT include: covariate list, compartment count, IIV structure,
   residual-error form, data origin, study counts, PKNCA sentence, or
   anything else that lives in the model file's metadata or vignette. A
   reviewer who wants those details clicks through to the model file.
8. Commit the model file, the vignette under `vignettes/articles/`, the regenerated `modeldb.rda` / `modeldb.qs2` / `modeldb.Rd`, the `NEWS.md` entry, and any updates to `inst/references/covariate-columns.md` (if a new covariate was registered) together on the feature branch.
9. Push the branch and open a PR against `main`. Use `gh pr create` with a title like `Add <Author> <Year> <drug> model`. **Before pushing, confirm steps 4 (vignette render) and 5 (non-ASCII check) were both run on the final state of the model and vignette files and both exited clean.** If either file was edited after the last gate run — even for a typo fix or a cosmetic comment change — re-run both gates; do not push on the strength of an earlier check against an older copy.

## Stop-and-ask triggers

The full consolidated list is in `references/pre-flight-checklist.md`. Review it once at dispatch (before starting Phase 1) so each trigger is visible up-front rather than discovered inline.

Use this fixed format for any ambiguity:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

When running interactively, use `AskUserQuestion` and wait for the answer. When running under a task runner, use the runner's documented stop-and-ask protocol. Either way: **stop, ask, wait — do not guess past the trigger.**

## Refusal handling

If at any phase you find yourself producing a content-policy refusal in response to legitimate clinical-trial content (drug + dose + indication + adverse-event language is normal in pharmacometric papers and is not a safety concern), do not silently degrade the extraction. Sidecar-ask:

> Phase <N> step <M> hit a content-policy refusal while reading <paper section / table>. The content is standard published clinical pharmacology (drug + dose + indication + AE language). Options: (A) rephrase the extraction prompt and continue, (B) skip the affected section and document the gap in vignette Errata, (C) defer the task pending operator review.

A refusal is operator-actionable signal, not a fatal error. Treat it the same as any other stop-and-ask trigger.

## Constraints

- Never invent parameter values; if it's not in the source, ask.
- Never tune parameters to make a validation output match a target.
- This skill **only adds** new models. Retrofitting existing vignettes to use PKNCA, or renaming covariates in existing files, is a future separate skill.
- Never push directly to `main`. Open a PR.
