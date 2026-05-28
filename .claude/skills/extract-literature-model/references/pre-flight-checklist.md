# Pre-flight checklist — stop-and-ask triggers

Read this once at dispatch, before starting the 6-phase workflow. Every trigger below is a hard stop: sidecar-ask the operator and wait. A trigger is **not advisory** — silent best-guesses ship as bugs the operator cannot retroactively correct without re-running the whole extraction.

## Source acquisition (Phase 1)

- **Lead PDF missing and OA acquisition failed.** All 5 ladder sources tried, no valid PDF whose title matches the task expectation. Sidecar with the structured attempts log from `scripts/acquire-paper.R` (`acquire-log.json`).
- **Paper-identity mismatch.** The on-disk file's title / first-author / journal / year / drug disagrees with the task's `Paper metadata` block. Sidecar with both actual and expected metadata.
- **Abstract only on disk.** The trimmed `.md` is < ~3 KB or the raw PMC XML / PDF only contains front matter + abstract. Population-PK parameter values are not in an abstract.
- **Upstream popPK dependency missing.** The current paper's PD model fixes PK parameters from a separate publication that is not on disk. Sidecar with the identifiable upstream paper (queue as dependency) or with the unidentifiable case (operator decision).
- **Systematic review / meta-analysis source.** The source catalogs other authors' models without developing an original one. Default action: skip and queue cited primary papers. Sidecar with the cited-primary list.

## Multi-model and ambiguity (Phase 1)

- **Multiple non-hierarchical models in one source.** Base + final → extract final; per-subpopulation, per-endpoint, sensitivity-analysis → list candidates and ask which to extract.
- **Parameter values look like initial estimates, not final.** Confirm against published point-estimate tables before treating any `$THETA` / `$OMEGA` initial as final.
- **Covariate encoding ambiguous.** Reference category, units, transformation not fully specified.
- **Source column name not in `inst/references/covariate-columns.md`.** Propose a new entry and confirm before adding.
- **Source column is an alias of an existing canonical name with value inversion or reference-category flip.** Confirm sign and reference-category implications.
- **Parameter name not in `references/parameter-names.md`.** Trivial notation differences (case, NONMEM `V1`/`V2`/`V3` → `vc`/`vp`/`vp2`, paper's `Kel` → canonical `kel`, etc.) translate silently. Any *new structural concept* (new clearance-component suffix like `cl_renal_intermittent`, new transform prefix, new endogenous-system parameter family, etc.) requires a sidecar BEFORE drafting the model. Propose: canonical name, role (one sentence), source paper's local name(s) it replaces, and any cross-precedent in existing registered models. Do not invent a "seems obvious" name silently.
- **Compartment name not in `references/compartment-names.md`.** Trivial casing differences translate silently. Any *new compartment role* (new endogenous-state compartment, new drug-suffix pattern beyond the registered list, new effect-compartment variant) requires a sidecar BEFORE drafting. Never introduce numbered `cmt1` / `cmt2` / `compartment_3` placeholders silently — propose a canonical role-based name and confirm. Same applies to drug-specific suffix patterns: don't extend `_tfvdp` → `_newdrugX` without operator sign-off on the suffix.

## Worktree state (Phase 2)

- **Worktree branch already pushed in a prior run.** `git log origin/<branch>..HEAD` shows no new commits and the branch exists upstream with task content. Sidecar to decide: verify and exit, or tear down and re-extract.

## Parameter sourcing (Phase 4)

- **Required parameter absent from every on-disk source.** Don't substitute from training data. Sidecar with options: (A) author-correspondence email, (B) approximate (QSS / steady-state / fixed-from-class) and document in vignette Errata, (C) skip.
- **`$THETA` reported with 0% RSE but no FIX flag.** Could be fixed or estimated; the encoding choice (`fixed()` wrapper or not) materially affects re-fits.
- **NCA disagreement > ~20%.** PKNCA output disagrees with published NCA after careful review. Do not tune — confirm the source has been correctly transcribed first.

## PBPK / QSP / MBMA specifics (Phase 1 Step 3a)

- **PBPK / QSP parameter not in any on-disk source.** Stricter than popPK: never substitute from "a representative PBPK paper for the same drug class" unless that paper is explicitly cited AND on disk. Sidecar instead.
- **Source only states "we used SimCYP's standard whole-body model" without writing out ODEs.** The model is not reproducible from on-disk sources. Sidecar or skip — do not fill in from training-data knowledge of typical SimCYP parameterizations.
- **MBMA paper reports only summary parameter point estimates without per-study weights / variance estimates.** Sidecar before drafting.

## Environment failures (Phase 6)

- **C-level segfault during `devtools::check()` or vignette rendering.** Broken R / rxode2 / nlmixr2 install. Sidecar — do not paper over with `--no-build-vignettes`.
- **`timeout 300 rmarkdown::render()` exits 124.** Vignette > 5 minutes. Reduce simulation size before committing; do not skip the time budget.

## Erratum search inconclusive

- **Paywalled erratum journal or ambiguous correction notice.** Ask the user whether any corrections apply before treating the main paper's values as final.

## Refusal handling

- **Content-policy refusal on legitimate clinical-trial content.** Drug + dose + indication + AE language is normal in pharmacometric papers and is not a safety concern. Sidecar with: (A) rephrase and continue, (B) skip the affected section and document in Errata, (C) defer pending operator review.

## Ambiguity format

When a stop-and-ask trigger fires, use this fixed format:

> Ambiguity at [source location]. Two plausible interpretations: (A) …, (B) …. Which applies?

Stop, ask, wait. Do not guess past the trigger. If you find yourself thinking "I'll just pick the safest option and move on," that itself is a stop-and-ask signal.
