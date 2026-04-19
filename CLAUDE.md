# nlmixr2lib — Claude project notes

This file gives AI coding assistants a fast orientation to the
repository. Human contributors should read
`vignettes/create-model-library.Rmd` for the full conventions and the
README for a project overview.

## What this package is

`nlmixr2lib` is a model library for `nlmixr2`. It ships curated
pharmacometric models (PK, PD, endogenous, therapeutic-area-specific)
that users can load by name and combine or modify.

## Layout highlights

- `inst/modeldb/` — every distributed model, organized by category
  (`specificDrugs/`, `endogenous/`, `pharmacokinetics/`,
  `pharmacodynamics/`, `therapeuticArea/`). The function name inside
  each `.R` file must match the filename.
- `R/modeldb.R` — the registry-building machinery (`buildModelDb()`,
  [`readModelDb()`](https://nlmixr2.github.io/nlmixr2lib/reference/readModelDb.md),
  [`modellib()`](https://nlmixr2.github.io/nlmixr2lib/reference/modellib.md)).
- `vignettes/` — validation vignettes, one per published model.
- `vignettes/create-model-library.Rmd` — the canonical naming
  conventions document.

## Skills

Repo-scoped skills live under `.claude/skills/`:

- **`extract-literature-model`** — guided workflow for extracting a
  pharmacometric model from a scientific source (journal article,
  supplement, poster, regulatory document) into the package. Use when a
  user provides a paper and asks to add the model. See
  `.claude/skills/extract-literature-model/SKILL.md`.

The skill’s `references/` folder contains the templates and standards it
enforces. The authoritative covariate-column register is
`inst/references/covariate-columns.md` — consult it before introducing
any new covariate column. The register is installed with the package so
[`checkModelConventions()`](https://nlmixr2.github.io/nlmixr2lib/reference/checkModelConventions.md)
can parse it at runtime; update the file (not R code) when ratifying a
new canonical covariate.

## Conventions (quick reference)

For full details see `vignettes/create-model-library.Rmd` and
`.claude/skills/extract-literature-model/references/naming-conventions.md`.

- Compartments: `depot`, `central`, `peripheral1`, `peripheral2`,
  `effect`. Observation: `Cc`.
- PK parameters (log-scale): `lka`, `lcl`, `lvc`, `lvp`, `lq`,
  `lfdepot`; derived `ka`, `cl`, `vc`, `vp`, `q`, `kel`, `k12`, etc.
- IIV: `eta` + transformed name (`etalcl`, `etalvc`). Use `etalcl` even
  when the paper used `etacl`.
- Residual error: `propSd`, `addSd`; multi-output prefixes with the
  output (`CcpropSd`).
- Covariate columns: canonical names in
  `inst/references/covariate-columns.md`. Standardized choices include
  `SEXF` (1 = female), `ADA_POS` (1 = ADA-positive), and a
  `RACE_<GROUP>` prefix for race indicators.

## Git workflow

- Never push directly to `main`. Always branch and open a PR.
- Before regenerating `data/modeldb.rda` / `inst/modeldb.qs2`, sync with
  `origin/main` so upstream additions aren’t clobbered.
- `nlmixr2lib::buildModelDb()` regenerates the registry; commit the
  regenerated artifacts alongside the model file and vignette.

## Testing

- `devtools::check()` must pass. Model files are validated via
  `buildModelDb()` (filename ↔︎ function-name match, parseable `ini()` /
  `model()`).
- Vignettes must build cleanly; they use `PKNCA` (in `Suggests`) for NCA
  validation.
