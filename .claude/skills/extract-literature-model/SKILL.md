---
name: extract-literature-model
description: This skill should be used when the user wants to "add a model from a paper", "extract a pharmacometric model from the literature", "implement a published PK/PD model in nlmixr2lib", or provides a scientific article, conference poster, supplement, or regulatory review document and asks to add the model to nlmixr2lib. Guides source review, standardized model file creation under inst/modeldb/, in-file source-trace verification, and validation vignette with PKNCA NCA checks.
---

# Extract a pharmacometric model from the literature

Input: a scientific source describing a pharmacometric model (journal article, supplement, conference poster, or regulatory review).
Output: a packaged nlmixr2lib model file under `inst/modeldb/`, a validation vignette under `vignettes/articles/`, and updated registry artifacts — opened as a pull request against `main`.

Work through the six phases below. Stop and ask the user at any of the decision points called out explicitly; ambiguity is the main failure mode for this workflow, and silent assumptions are what get shipped as bugs.

## Tooling: `nlmixr2libingest` (assumed installed)

This skill uses the helper package **`nlmixr2libingest`** for token-efficient
canonical-name resolution (Phase 3) and a source-trace pre-check (Phase 4). It is
assumed to be installed; install from
<https://github.com/billdenney/nlmixr2libingest> if missing:

```r
remotes::install_github("billdenney/nlmixr2libingest")
```

Resolve its CLI directory **once** at the start of the task and reuse it (if the
result is empty the package is not installed — install it, then continue):

```bash
NLI=$(Rscript -e 'cat(system.file("scripts", package="nlmixr2libingest"))')
```

You will use `"$NLI/prebrief.R"` (Phase 3), `"$NLI/validate.R"` (Phase 4
source-trace pre-check, and Phase 6 as the combined per-iteration gate
`validate.R … --model`), and `"$NLI/lint_vignette.R"` (Phase 6, a static
pre-render lint). All are **priors / assists** — they never gate what you read of
the paper, and the quality firewall (source-trace every final value and any
nonstandard equation against the source) is unchanged.

## References

Each reference is loaded just in time. Most tasks load only the "Always" set plus one or two conditionals.

**Always — read at dispatch:**
- `references/pre-flight-checklist.md` — consolidated stop-and-ask triggers.
- `references/replicate-author-structure.md` — standing default policy on multi-model handling (N files vs 1 file, always 1 vignette per paper) with worked examples.

**Always — load by phase:**
- `references/compartment-names.md`, `references/parameter-names.md`, `inst/references/covariate-columns.md`, `references/model-file-template.md` — Phase 3 (drafting the model file). `parameter-names.md` is also the canonical home for author-surname normalization and fixed-parameter encoding rules. `inst/references/covariate-columns.md` is installed with the package so `checkModelConventions()` parses it at runtime.
- `references/verification-checklist.md` — Phase 4 (verifying the first-pass implementation).
- `references/vignette-template.md` — Phase 5 (drafting the validation vignette).
- `references/known-vignette-failure-patterns.md` — Phase 5 step 2 (before running the render gate). Catalogue of the recurring vignette-render failures that have shipped because the gate was skipped, fabricated, or run with events too simple to exercise the failure. Read this before your render, not after CI fails.

**Conditional — load only when the model class applies:**
- `references/nonmem-translation.md` — NONMEM → nlmixr2 syntax conversion. Load when the source is a NONMEM control stream (`.mod` / `.ctl`, supplement with `$PK` / `$DES` / `$THETA` blocks).
- `references/pknca-recipes.md` — PKNCA setups for single-dose, steady-state, and multi-dose NCA. Load for popPK validation; skip for endogenous / mechanistic models.
- `references/endogenous-validation.md` — validation strategy for endogenous / mechanistic / turnover models. Load when the source describes turnover, enzyme-kinetic, or steady-state-balance models where PKNCA is not the right check.
- `references/pbpk-qsp-mbma.md` — extra discipline for whole-body PBPK, QSP, and MBMA papers (description prefix, filename suffix, on-disk-only parameter sourcing, reproducibility check, platform-vs-nlmixr2, MBMA between-study variance ≠ popPK BSV).

**Troubleshooting — read on failure:**
- `references/oa-acquisition.md` — OA-PDF ladder source-by-source detail. Read when `scripts/acquire-paper.R` exits non-zero and you need to interpret the per-source failure log, or when extending the source list.

**Scripts:**
- `scripts/acquire-paper.R` — open-access PDF acquisition (Phase 1 Step 0).
- `scripts/lint-conventions.R` — pretty-print `checkModelConventions()` results (Phase 6 step 4).

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

If the paper describes a whole-body PBPK, semi-mechanistic PBPK, quantitative systems pharmacology (QSP), or model-based meta-analysis (MBMA) model, load `references/pbpk-qsp-mbma.md` now. That reference covers the `description` prefix and filename suffix conventions, the strict on-disk-only parameter-sourcing rule (stricter than popPK — never substitute from training-data defaults or "a typical SimCYP entry"), the reproducibility check to walk before drafting, the platform-vs-nlmixr2 ODE-encoding rule, and the MBMA between-study-variance ≠ popPK BSV distinction. Extract these models the same way you would a popPK model, with the extra discipline that reference imposes.

### Step 3b — Target file naming (author-surname normalization)

Apply normalization silently — no sidecar — when choosing the filename, function name, vignette basename, branch, and PR title. The full table and rules live in `references/parameter-names.md` § "Author-surname normalization" (canonical). Quick examples: `Lohy Das` → `LohyDas` (drop spaces); `Ait-Oudhia` → `AitOudhia` (drop hyphens); `O'Brien` → `OBrien` (drop apostrophes); `Câmara` → `Camara` (transliterate accents); `van Rongen` → `vanRongen` (lowercase particle preserved per published form); `Von Bonin` → `VonBonin` (capitalised particle preserved per published form).

The same normalised stem is used for the model filename (`<Stem>_<Year>_<drug>.R`), the function name (which `buildModelDb()` enforces must equal the filename minus `.R`), the `vignette <- "..."` field, the vignette basename under `vignettes/articles/`, the worktree branch name in Phase 2, and the PR title. Year-letter collision suffixes (`Author_2019a_drug.R`, `Author_2019b_drug.R`) for genuinely-same-author/year/drug pairs are handled separately — see Phase 3 and `references/parameter-names.md` § "Year-letter collision suffix".

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
9. **Verify parameters are final estimates, not initial estimates.** Supplement control streams usually list initial values in `$THETA` / `$OMEGA`; final values come from the paper's results table or `$TABLE` output. If only a control stream is available, confirm values against any published point estimates before treating them as final.
10. **Multiple-model handling — replicate the author's structure (default).** Extract the model exactly as the authors built it. Apply the N-files vs. 1-file mapping silently — no sidecar. See `references/replicate-author-structure.md` for the full decision table (independent vs. coupled, base-vs-final, sensitivity analyses), the two narrow sidecar conditions (infeasible to express in rxode2/nlmixr2; genuinely ambiguous paper text), and worked examples (de Vries Schultink 2018 → 2 `.R` files + 1 vignette; Yoshida 2021 ipatasertib → 1 `.R` file + 1 vignette). Regardless of how many `.R` files come out, produce one vignette per paper (Phase 5). Do NOT sidecar to ask "the paper has multiple models — which do you want?"
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

`<firstauthor>` here is the **normalised** surname stem from Phase 1 Step 3b (CamelCase across hyphens/spaces/apostrophes, accents transliterated, lowercase particles preserved per the published form). Use the same stem the model file and vignette use, just in the lowercase form typical for branch names (e.g. `lohydas-2018-mefloquine`, `vanrongen-2017-midazolam`).

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

File path: `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`. `<FirstAuthor>` is the **normalised** surname stem — see Phase 1 Step 3b for the quick-reference examples and `references/parameter-names.md` § "Author-surname normalization" for the full table and rules (hyphens/spaces/apostrophes dropped with CamelCase, accents transliterated, lowercase particles preserved per the published form). Do not stop to ask the operator to confirm the normalised form; the policy is a hard default.

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

Wrap fixed values in `fixed()` in `ini()` for ALL parameter types (structural THETAs, allometric exponents, IIVs, residual error, covariate effects). The `log()` goes **inside** `fixed()` for log-transformed parameters (`lcl <- fixed(log(2))`, never `log(fixed(2))`). Failing to encode the fixed status loses load-bearing provenance — a downstream user cannot tell whether the value was estimated or held constant.

See `references/parameter-names.md` § "Fixed parameters" for source signals (prose, NONMEM `FIX` flag, allometric exponents without uncertainty, F1 anchors, upstream-inherited parameters) and encoding examples. When in doubt — a $THETA reported with an RSE of 0% but no FIX flag, or a parameter reported to three decimal places with no uncertainty — sidecar-ask the operator before guessing.

Follow `references/compartment-names.md` and `references/parameter-names.md` strictly. The quick-reference rules:

- Structural PK parameters log-transformed: `lka`, `lcl`, `lvc`, `lvp`, `lvp2`, `lq`, `lq2`, `lfdepot`.
- IIV: `eta` + transformed name, e.g., `etalcl` (not `etacl`). Block correlations via `etalcl + etalvc ~ c(var, cov, var)`.
- Residual error: `propSd`, `addSd`, `expSd`. Multi-output: `propSd_<output>`, `addSd_<output>`, etc.
- Compartments: `depot`, `central`, `peripheral1`, `peripheral2`, `effect`, `target`, `complex`, `csf`, `isf`, etc. Observation: `Cc`.

**Stop-and-ask before introducing any new canonical parameter or compartment name.** The authoritative lists live in `references/parameter-names.md` (~263 lines covering structural PK params, transform prefixes, IIV/eta names, residual error, covariate-effect parameters, file-level metadata) and `references/compartment-names.md` (~71 lines covering canonical compartments and drug-suffix patterns). Treat these the same way as `inst/references/covariate-columns.md`:

- If the paper's local name is just notation for an existing canonical (e.g., paper uses `Kel`, canonical is `kel`; paper uses `V1/V2/V3`, canonical is `vc/vp/vp2`), translate to canonical and proceed — no sidecar needed for trivial casing or notation differences explicitly covered by the references.
- If the paper introduces a structural concept that ISN'T in the references at all (e.g., a new clearance component suffix, a new compartment role, a new transform prefix), **sidecar-ask before writing the model file**. Propose the new canonical name with: its role (one sentence), its source-paper local name(s) it would replace, why it isn't an alias of an existing canonical, and any cross-precedent in the existing registered models (cite filenames). Wait for operator approval before committing the model file and the new entry to `parameter-names.md` / `compartment-names.md`.
- Do NOT pick a "seems obvious" extension silently (e.g., choosing `kel_distinct_<suffix>` when nothing in the registry uses that pattern, or registering `compartment_NN` numbered names instead of role-based names). Even when the new name looks like a natural extension of an existing one, file the sidecar so the operator can decide whether it's actually a new canonical or a synonym of something that already exists.
- After the operator approves, the new entry is committed alongside the model file in the same PR (same convention as `covariate-columns.md`). Append the new name to the appropriate H2 section in the reference file with a one-paragraph description and a `Founding example: <Author_Year_drug.R>` citation. Do not add a change-log or history section — `git log` is the record.

Covariate columns must use the canonical names in `inst/references/covariate-columns.md` (~1.1 MB ≈ 284k tokens — **do NOT `Read` it whole; one whole-file read can exhaust the task's token budget**). While reading the paper you identified the covariates it uses; resolve them **all at once** in a single batched pre-brief rather than a lookup per name (resolve parameter and compartment names the same way, with `parameter` / `compartment`, instead of reading those reference files in full):

```bash
Rscript "$NLI/prebrief.R" covariate "body weight" "serum albumin" "creatinine clearance"
# first arg = kind (covariate | parameter | compartment); the rest are the paper's terms (quote multi-word phrases)
```

It prints a `term -> CANONICAL [units, scope]` line per term, or `UNMATCHED` when a term has no canonical entry. Use the printed canonical names. For an `UNMATCHED` term (a possible NEW canonical → stop-and-ask), or to double-check a single name, fall back to the per-name lookup:

```bash
Rscript --vanilla "$NLI/lookup.R" "body weight" covariate
```

Before writing any covariate into the file:

- If the canonical name exists, use it and record the source column name in `covariateData[[name]]$source_name`.
- If the source name is an alias of an existing canonical name (e.g., source uses `SEXM`, canonical is `SEXF`), use the canonical name, note the required value transformation (`SEXF = 1 - SEXM`), and **ask the user to confirm the effect-coefficient sign and reference-category implications** before committing.
- If the concept isn't in the register at all, propose a new entry (canonical name, description, units, type, reference category, source aliases) and ask the user to confirm before adding it. The new entry is committed alongside the model.
- Do **not** add a change-log / history / "## Change log" / "## Summary" section or per-extraction history line to `inst/references/covariate-columns.md`. The H3 entry itself is the authoritative record; chronological history of when an entry was added or modified is read from `git log`. Any context that someone needs to use the covariate (derivation rules, scope-promotion rationale, name-collision history, naming-decision sidecars) belongs in the H3 entry's Description / Notes / Source aliases — not in a separate change log.

**Documented-but-unused covariates → `covariatesDataExcluded`.** When the source paper *screens* a covariate (e.g. a Full Random Effects Model / FREM forest-plot sweep of demographics) but does **not** retain it in the final model — because the effect is not clinically meaningful and/or the paper reports only graphical estimates with no usable point estimate — do not put it in `covariateData`. A `covariateData` entry that is never referenced in `model()` triggers a "declared but not referenced" convention warning. Instead, document these in a separate metadata list named `covariatesDataExcluded`, with the identical entry shape (`covariatesDataExcluded <- list(WT = list(description=..., units=..., type=..., notes="screened in the FREM but not retained; see ..."), ...)`). `checkModelConventions()` treats this list as documentation only: it does not require the entries to appear in `model()` and does not flag them as unused. The one rule the checker enforces is that a name in `covariatesDataExcluded` must **not** also be referenced in `model()` (if it is actually used, it belongs in `covariateData`). Use this to preserve the provenance of a paper's covariate screen without carrying convention warnings.

`population` uses the extensible schema documented in `references/parameter-names.md` § "File-level metadata". Common fields: `n_subjects`, `n_studies`, `age_range`, `weight_range`, `sex_female_pct`, `race_ethnicity`, `disease_state`, `dose_range`, `regions`, `notes`. Add any additional keys the paper describes (e.g., `ga_range`, `renal_function`, `co_medication`) — do not force facts into the common schema.

## Phase 4 — Verification (re-read the source)

First run the **source-trace pre-check** — it parses the model, runs
`checkModelConventions()`, and flags every final `ini()` value with no supporting
number anywhere in the paper (back-transform- and rounding-tolerant), plus the
nonstandard `model()` equations and hardcoded constants to confirm:

```bash
Rscript "$NLI/validate.R" inst/modeldb/<category>/<Stem>.R <trimmed-paper.md>
```

Treat its output as a **worklist, not a verdict**: a "found" match can be
coincidental, and an unverified flag can be a legitimate derived value (e.g. a
unit conversion) — it focuses your manual source-trace, it does not replace it,
and the firewall (every value and nonstandard equation traced to the source) is
unchanged. Then re-read the source independently and walk through
`references/verification-checklist.md`. Common pitfalls:

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

**Non-paper provenance — annotate inline.** When a parameter value did not come from the paper's text or tables (e.g. the operator read it off a graphical figure, an author supplied it via email correspondence, or it was carried from an upstream-task model file), record the provenance as an inline comment on the parameter line. See `references/model-file-template.md` § "Documenting non-paper-derived parameter values" for the comment forms. This is mandatory for any parameter not from the paper text/tables. The vignette's Assumptions and deviations / Errata section must also list the non-paper provenance (see `references/vignette-template.md`).

## Phase 5 — Validation vignette

File path: `vignettes/articles/<FirstAuthor>_<Year>_<topic>.Rmd`. Drug-specific vignettes live under `vignettes/articles/` so pkgdown builds them for the site but CRAN does not — `.Rbuildignore` excludes that directory. The basename (without `.Rmd`) must match the `vignette <- "..."` field in every model file the paper contributed. Only the legacy `PK_2cmt_mAb_Davda_2014.Rmd` remains at top-level `vignettes/`.

**One vignette per paper, regardless of N model files.** When a paper produces N `.R` files (per the Phase 1 step 10 policy — independent models extracted as N separate files), produce **a single vignette** that walks the paper's narrative and uses each `modellib()` call at the appropriate point. Do NOT produce N vignettes (`Author_Year_drug1.Rmd` / `Author_Year_drug2.Rmd`) for a paper that contributed N models — that fragments the reviewer's view of the paper-as-a-unit.

- **Naming**: `<FirstAuthor>_<Year>_<topic>.Rmd`. For a single-model paper, `<topic>` is the drug (e.g. `Hu_2026_clesrovimab.Rmd`). For a multi-model paper, `<topic>` is a phrase that names what the paper is about as a whole, not one of its drugs / endpoints (e.g. `de_vries_schultink_2018_anthracycline_trastuzumab_cardiotoxicity.Rmd`, NOT `de_vries_schultink_2018_anthracyclines.Rmd` plus a second vignette). When in doubt, pick the most prominent shared subject of the paper (a disease, a mechanism, a combination therapy).
- **`vignette <- "..."` in every contributing model file** points to the same basename. `buildModelDb()` records the vignette per model so all N models in `modellib()` link to the same rendered page.
- **YAML `title:` and `\VignetteIndexEntry{...}`** use the **human form** of the paper's subject and citation, e.g. `Anthracycline + trastuzumab cardiotoxicity (de Vries Schultink 2018)` — drug list at the start, citation in parentheses. For single-model papers the existing form `<Drug> (<FirstAuthor> <Year>)` (e.g. `Ustekinumab (Aguiar 2021)`) continues unchanged. The filename stays in the machine form. See `references/vignette-template.md` for header rules.

Use `references/vignette-template.md`. Required sections, in order:

1. **Header and setup** — libraries include `nlmixr2lib`, `PKNCA`, `rxode2`, `dplyr`, `ggplot2`.
2. **Population** — narrative reproducing the `population` metadata; cite the source table listing baseline demographics.
3. **Source trace** — a dedicated table listing the source location (page / table / equation / figure) for every model equation and every `ini()` parameter. This is in addition to the in-file comments; the vignette gives reviewers a single place to audit provenance.
4. **Virtual cohort** — covariate distributions match the population metadata. Use WHO weight-for-age curves for pediatric models.
5. **Simulation** — `rxode2::rxSolve(mod, events)` for stochastic VPCs; `rxode2::zeroRe()` + `rxSolve` for typical-value replications. **Keep cohorts small: never simulate more than 200 participants per arm** (per treatment / dose group). 200/arm is ample for a VPC; larger cohorts add no validation value and are the top render-timeout and token-cost cause. `references/vignette-template.md` defaults to this cap, and `lint_vignette.R` (Phase 6) flags any cohort over it.
6. **Replicate published figures** — one code chunk per figure, caption linking to the source figure number ("Replicates Figure 4 of <Author Year>").
7. **PKNCA validation** — required; no inline trapezoidal NCA. See `references/pknca-recipes.md`. The PKNCA formula **must include a treatment grouping variable** (`conc ~ time | id/treatment`) so per-group results can be compared against the paper. The PKNCA input filter must be `dplyr::filter(!is.na(Cc))` only — adding `time > 0` or `Cc > 0` drops the time-zero row and triggers the "Requesting an AUC range starting (0) before the first measurement" warning across every subject; if the simulation grid doesn't naturally produce a time-zero row, add one defensively (see `pknca-recipes.md` § "Time-zero records (mandatory)").
8. **Comparison against published NCA** — if the source paper reports Cmax / Tmax / AUC / half-life, render **one** combined side-by-side table via `nlmixr2lib::ncaComparisonTable()` with header column `NCA parameter` and friendly parameter labels (`Cmax`, `AUC0-∞ (obs)`, `t½`, …). Do **not** split simulated and reference values across separate tables or use "see above" cross-references. Flag any starred (>20% difference) rows in the narrative and investigate the source — do not tune.
9. **Assumptions and deviations** — explicit list of what you had to assume because the paper didn't say (race distribution, z-score stability, etc.).

For **endogenous / turnover models** where NCA isn't the right validation, replace the PKNCA section with the steady-state / perturbation-recovery / mass-balance checks described in `references/endogenous-validation.md`. For **multi-output models**, run one PKNCA block per output.

### Endogenous and mechanistic models

For papers that describe endogenous turnover, steady-state-balance, or mechanistic enzyme kinetics (e.g., Kim 2006 IgG FcRn recycling, Charbonneau 2021 phenylalanine), load `references/endogenous-validation.md` for vignette structure and the four validation patterns (steady-state hold, perturbation recovery, mass-balance / flux check, dimensional analysis). Use the mechanistic parameter names from `references/parameter-names.md` § "Endogenous / mechanistic parameters" (`Vmax`, `Km`, `kint`, `kcat`, baseline `bl_<species>`, fractional-activity `f_<enzyme>`). These models typically have no IIV, no residual error, and no dosing events — the state starts at a biological baseline. Replace the PKNCA section of the vignette with the steady-state / perturbation-recovery / mass-balance / dimensional-analysis checks.

## Phase 6 — Registration, tests, docs, PR

**Hard rule for the entire phase: if the vignette does not render cleanly, the task is BLOCKED.** No `git add` of the model file or vignette, no `git commit`, no `git push`, no PR description, no "renders in ~Xs" claim in any commit message or PR body — until step 2 below returns exit 0 with an HTML file ≥ 1 KB on disk. A failed render is not a soft warning that gets resolved in CI; it is a stop-the-line condition. CI catching the failure for you is a process violation, not an acceptable workflow. Past extractions (e.g. HillMcManus 2017 in PR #433, which failed CI on `Assertion on 'amt' failed: May not have names.` despite a commit message claiming "vignette renders in ~11 s") have demonstrated that fabricating this evidence breaks the entire downstream merge / release pipeline. Do not do it.

1. Re-confirm the branch is on top of fresh `origin/main` (`git fetch origin && git rebase origin/main` if needed).

2. **Render the vignette locally and verify it completes without error in under 5 minutes of wall-clock time. This is the single highest-value pre-push gate.** Run this BEFORE `buildModelDb()`, BEFORE the convention lint, BEFORE any `git add` of the model file or vignette — if the render is broken, fixing it is the only thing that matters. The gate catches the failure modes pkgdown CI surfaces (missing data columns, time-varying covariates assigned to rxEt objects that get silently dropped, PKNCA formulas referencing absent columns, simulation crashes, named-scalar arguments to `rxode2::et()`, vignettes exceeding the time budget) and surfaces them in seconds.

   **Before you run the gate, read `references/known-vignette-failure-patterns.md`.** Those are the recurring shapes that have shipped broken because the gate was skipped, fabricated, or run with events too simple to exercise the failure: singular OMEGA (`chol(): decomposition failed`), missing explicit `cmt()` declarations when the model has algebraic observables and the vignette references them on observation rows (12-of-15 failures in the 2026-06-17 consolidation), dplyr `unique()` on a varying column, PKNCA zero-row filters, callr timeouts under parallel build. Scan that doc, apply the prophylactic fixes, then render. The consolidation merge (`runner-merge-claude-branches`) re-runs every vignette in parallel as a HARD gate — anything you ship broken WILL be caught there, but the cost of un-stacking a fix from a 130-branch merge is much higher than fixing it in your one-paper worktree right now.

   **Cover the failure modes in your event table, not just the dose path.** A render that only doses (no observations) or only observes with `cmt = NA_character_` will pass the gate and still ship a broken vignette. Make sure your event table includes observation rows **on the relevant ODE state** (`cmt = "central"`, etc.) — rxode2 returns each algebraic observable (e.g. `Cc`) as a column at those rows, so observing the state exercises the observable-computation path. Do **NOT** write `cmt = "<observable-name>"` (e.g. `cmt = "Cc"`): referencing an observable as a compartment auto-injects a `cmt()` slot for it AFTER the ODE states and renumbers every compartment — that IS the slot-renumbering bug, not a way to test it (see `references/known-vignette-failure-patterns.md` patterns; `lint_vignette.R`, next paragraph, flags it).

   **Pre-lint, then iterate with the combined gate — both before the evidence render below.** Two `nlmixr2libingest` helpers make the fix loop cheap:

   - **Static pre-lint (milliseconds, no render).** Catches the most common render-killers up front — `cmt =` on an algebraic observable, a named-vector character subscript to `amt =`, a PKNCA `time > 0` / `Cc > 0` filter, and a cohort over the 200-per-arm cap:
     ```bash
     Rscript "$NLI/lint_vignette.R" vignettes/articles/${STEM}.Rmd inst/modeldb/<category>/${STEM}.R
     ```
     Fix what it flags. It is a pre-check, not a substitute — a clean lint does not guarantee a clean render.

   - **Combined model-scoped gate (one turn, no whole-package check).** During the fix loop, run parse + `checkModelConventions()` + source-trace + `load_all` + the vignette render in a SINGLE call instead of separate render / convention-lint turns (each of which reloads the package in its own turn):
     ```bash
     Rscript "$NLI/validate.R" inst/modeldb/<category>/${STEM}.R <trimmed-paper.md> --model --pkg . --vignette vignettes/articles/${STEM}.Rmd
     ```
     It returns one terse Success / fix-list. When it finally reports Success, run the explicit `RENDER_GATE` command below **once** to capture the verbatim evidence line that step 9 requires (the combined gate is for iterating; the `RENDER_GATE` line is the PR evidence).

   **Run the render SYNCHRONOUSLY, in one tool call (as below). NEVER launch it with `run_in_background`, and NEVER poll a running render's log (`cat …log`, `tail -f`, repeated reads) across multiple tool calls.** Each poll re-reads the entire accumulated context, so waiting on a slow render costs far more in cache-read tokens than the extraction itself. (A single 3-model statin vignette did exactly this — one backgrounded render polled ~98 times across ~250 turns → 71M cache-read / ~$42 for one task, ~25× a normal extraction.) If the render is slow, **fix it** (timeout bullet below); do not wait on it.

   Run the command, capture stdout/stderr to a log file, then check BOTH the exit code AND the HTML byte size:

   ```bash
   STEM="<FirstAuthor>_<Year>_<drug>"
   LOG="/tmp/vignette-render-${STEM}.log"
   timeout 300 Rscript --vanilla -e "
     pkgload::load_all('.', quiet = TRUE)
     rmarkdown::render(
       'vignettes/articles/${STEM}.Rmd',
       quiet = FALSE
     )
   " 2>&1 | tee "$LOG"
   RC=${PIPESTATUS[0]}
   HTML="vignettes/articles/${STEM}.html"
   BYTES=$(stat -c%s "$HTML" 2>/dev/null || echo MISSING)
   echo "RENDER_GATE stem=$STEM exit=$RC html_bytes=$BYTES"
   Rscript --vanilla .claude/skills/extract-literature-model/scripts/filter-r-log.R "$LOG" > "${LOG%.log}.filt.log" 2>/dev/null
   echo "RENDER_GATE_FILTERED=${LOG%.log}.filt.log"
   ```

   **Read the filtered log (`${LOG%.log}.filt.log`), not the raw `$LOG`.** It compacts the render log to its errors/warnings/traceback (the raw log is often hundreds–thousands of lines that get re-read every turn of the fix loop). Open the raw `$LOG` only if the filtered view is insufficient.

   Interpret the printed `RENDER_GATE` line:

   - **`exit=0` AND `html_bytes` is a number ≥ 1024** → gate passed; proceed to step 3.
   - **`exit` non-zero (R error)** → STOP. Do not `git add` anything. Read the traceback in `$LOG`; the most common causes are (a) a `select()` / PKNCA formula referencing a column that was never created (add it in the preceding `mutate()`), (b) a variable name used before it is defined (reorder chunks), (c) a simulation that errors because covariate values are out of range for the model, (d) `event_table$col <- vector` syntax against an `rxEt` object — rxode2 silently drops these assignments, so always materialize via `as.data.frame()` and add covariate columns there before passing to `rxSolve()`, (e) a named scalar passed as `amt =` to `rxode2::et()` (e.g. `doses["fbx"]` instead of `doses[["fbx"]]`) — rxode2 rejects names with `Assertion on 'amt' failed: May not have names.`. Fix the vignette, re-run this step, and only proceed when both `exit=0` and `html_bytes ≥ 1024` come back clean.
   - **PKNCA "Requesting an AUC range starting (0) before the first measurement" warning** in `$LOG` (often repeated once per subject) → the concentration frame passed to PKNCA has no `time = 0` row. Two causes: (1) the PKNCA input filter is too aggressive (typically `dplyr::filter(time > 0, ...)` or `Cc > 0`) — use only `!is.na(Cc)`; (2) the simulation grid never produced a time-zero observation — add one defensively via the bind_rows + distinct pattern in `pknca-recipes.md` § "Time-zero records (mandatory)". Render again until the warning is gone.
   - **`exit=124` (timeout)** → the vignette exceeds 5 minutes. Reduce simulation size: cut `nSub` (stochastic VPC subjects) — **the cap is 200 participants per arm; if any cohort is above it, that is the first thing to cut** — shorten the observation grid, or move expensive chunks behind `eval = FALSE` with a note. Do **not** skip the time budget — pkgdown CI has strict wall-time limits and a slow vignette breaks the build for everyone. **A timeout means SHRINK the simulation, never background-and-poll the render. For a multi-model (N-drug) vignette the combined cohort simulation is the usual culprit — cut the per-model cohort size first.** Re-run the gate synchronously after shrinking.
   - **`html_bytes=MISSING`** → render exited but no HTML was written. Treat as a failure even if `exit=0`; pandoc may have aborted silently. Re-run with `quiet = FALSE` (already set above) and read the bottom of `$LOG`.
   - **C-level segfault (`*** caught segfault ***`)** → broken R / rxode2 / nlmixr2 install, not a model-file problem. Stop, sidecar-ask the operator to fix the environment; do not work around with `--no-build-vignettes`.

   Do not interpret silence as success — read `$LOG`, confirm the final lines show a successful pandoc invocation, and confirm the HTML file ≥ 1 KB landed next to the `.Rmd`. Do not assume "checkModelConventions was clean, the vignette must be fine" — those are independent failure modes. Do not claim "vignette renders in ~Xs" in the commit message or PR body unless you actually observed `exit=0` with `html_bytes ≥ 1024` in the current shell on the current file contents; fabricating that line is a process violation that masks real CI failures (and HAS done so — see the hard rule above). Always re-run the gate after any edit to the vignette or the model file, even cosmetic edits — pkgdown CI does, and you must too.

   **Record the result.** Copy the printed `RENDER_GATE stem=... exit=... html_bytes=...` line verbatim — it goes into the PR body in step 9 as evidence the gate actually ran on the final file contents.

3. Run `nlmixr2lib::buildModelDb()` to regenerate `data/modeldb.rda` and `inst/modeldb.qs2`. Confirm the new model appears in `modellib()`. **When verifying in R, do `devtools::load_all(".")` first so `modellib()` reads the worktree's in-development package, not the stale system install** — see `references/verification-checklist.md` § "Verifying against the worktree's nlmixr2lib" for why a bare `library(nlmixr2lib)` can return a misleading `FALSE`. (`load_all()` is required in every fresh `Rscript` session for this reason; that is correct setup, not overhead.)

   **Run `devtools::document()` at most once, at the very end** — after the model file and vignette are final, just before the final `devtools::check()`. Roxygen reprocesses the whole package on every call, so running it inside the edit/fix loop is wasted work; the model and vignette do not need re-documenting after each tweak.

4. Run the convention lint and review the output:

   ```bash
   Rscript .claude/skills/extract-literature-model/scripts/lint-conventions.R \
     "<FirstAuthor>_<Year>_<drug>"
   ```

   The script wraps `nlmixr2lib::checkModelConventions()` and emits PR-paste-ready output. Exit 0 = clean; exit 1 = warnings only; exit 2 = errors. Any deviations from the canonical parameter / IIV / residual-error / covariate / compartment conventions (see `references/compartment-names.md`, `references/parameter-names.md`, and `inst/references/covariate-columns.md`) should be either fixed in the model file before committing, or explicitly justified in the vignette's Assumptions and deviations section. `buildModelDb()` runs `checkModelConventions()` implicitly at package-build time, but running the lint explicitly on your new model surfaces drift before commit. Paste the lint output into the PR body so a reviewer can see what was checked.

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

6. Run `devtools::check()`. Vignettes must build cleanly. Capture the output and read the **filtered** view — a full check log is ~1,000–5,000 lines and gets re-read every turn:

   ```bash
   Rscript --vanilla -e 'devtools::check(args = "--no-build-vignettes", error_on = "never")' 2>&1 | tee /tmp/check.log
   Rscript --vanilla .claude/skills/extract-literature-model/scripts/filter-r-log.R /tmp/check.log
   ```

   Read the filtered errors/warnings/NOTEs; open `/tmp/check.log` only if the summary is insufficient.

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
9. Push the branch and open a PR against `main`. Use `gh pr create` with a title like `Add <Author> <Year> <drug> model`. **Before pushing, confirm steps 2 (vignette render) and 5 (non-ASCII check) were both run on the final state of the model and vignette files and both exited clean.** If either file was edited after the last gate run — even for a typo fix or a cosmetic comment change — re-run both gates; do not push on the strength of an earlier check against an older copy.

   The PR body MUST include the verbatim `RENDER_GATE stem=... exit=0 html_bytes=...` line printed by step 2 (and a matching line for any re-runs after edits). Reviewers use this line to confirm the gate actually ran on the file contents in the PR head — a missing or fabricated line is grounds for the reviewer to request a fresh render before merge.

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
