# Validation vignette template

File path: `vignettes/articles/<FirstAuthor>_<Year>_<topic>.Rmd` — the pkgdown "articles" directory, not top-level `vignettes/`. Drug-specific vignettes live under `vignettes/articles/` so pkgdown builds them for the site but CRAN skips them (`.Rbuildignore` excludes the directory). The top-level `vignettes/` directory is reserved for the one legacy vignette (`PK_2cmt_mAb_Davda_2014.Rmd`). The basename must match the `vignette <- "..."` field in every model file the paper contributed.

**One vignette per paper.** Even when a paper contributes N model files (per the `replicate-author-structure.md` policy — independent models extracted as N separate `.R` files), produce **a single vignette** that walks the paper's narrative as a unit and uses each `modellib()` at the appropriate point. Do NOT produce N vignettes for a paper that contributed N models — that fragments the reviewer's view and duplicates the population / methods narrative.

**Filename vs. title.** The **filename** stays in the machine form `<FirstAuthor>_<Year>_<topic>.Rmd`. For a single-model paper, `<topic>` is the drug (e.g. `Hu_2026_clesrovimab.Rmd`). For a multi-model paper, `<topic>` names what the paper is about as a whole (e.g. `de_vries_schultink_2018_anthracycline_trastuzumab_cardiotoxicity.Rmd`), not any one of its drugs / endpoints. It is what `buildModelDb()` matches against the `vignette <- "..."` field in each contributing model file, what the modeldb `vignette` column stores, and what sorts the navbar dropdowns. The YAML `title:` field and `\VignetteIndexEntry{...}` use the **human form**:

- **Single-model paper**: `<Drug> (<FirstAuthor> <Year>)` (e.g. `Ustekinumab (Aguiar 2021)`).
- **Multi-model paper**: `<paper-subject> (<FirstAuthor> <Year>)` (e.g. `Anthracycline + trastuzumab cardiotoxicity (de Vries Schultink 2018)`) — the drug list or shared subject at the start, citation in parentheses.

The two surfaces — link text and page title — therefore agree, and the list-of-models vignette renders the same form. For papers that warrant a richer descriptive title (mechanistic ADC + payload-coupled models, intrathecal-delivery scenarios, multi-agent regimens), expand the human title — but always keep `(<FirstAuthor> <Year>)` as the trailing parenthetical and keep `title:` and `\VignetteIndexEntry{...}` byte-identical.

A validation vignette has two jobs:

1. **Verification trail** — document exactly where each model equation and parameter came from.
2. **Validation** — show that a simulation from the packaged model reproduces published results (usually one or more figures and an NCA table) using **PKNCA** for NCA parameters.

Use this template as a starting point. Every section is required unless marked optional.

## Template

````markdown
---
title: "<Drug> (<FirstAuthor> <Year>)"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{<Drug> (<FirstAuthor> <Year>)}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(nlmixr2lib)
library(PKNCA)
library(rxode2)
library(dplyr)
library(tidyr)
library(ggplot2)
```

# Model and source

* Citation: `r rxode2::rxode(readModelDb("<FirstAuthor>_<Year>_<drug>"))$reference`
* Description: `r rxode2::rxode(readModelDb("<FirstAuthor>_<Year>_<drug>"))$description`
* Article: [<Journal short citation>](<direct DOI or publisher URL>)

<!--
Always include a direct link to the original article (DOI preferred) so readers
can jump to the source without searching for the citation. If a free
supplementary document exists (open access supplement, regulatory review,
author-hosted preprint), add additional lines for each. Examples:

* Article: <https://doi.org/10.1002/cpt.3447>
* Supplement: <https://doi.org/10.1002/cpt.3447-sup-0001>
-->


# Population

<One or two paragraphs summarizing the study population that informed the model:
N subjects, age / weight range, sex balance, race / ethnicity distribution,
disease state, dose levels, regions, and any relevant baseline characteristics.
Cite the source Table that lists baseline demographics.>

The same information is available programmatically via the model's `population`
metadata (`readModelDb("<model>")$population` after the model is loaded).

# Source trace

The per-parameter origin is recorded as an in-file comment next to each `ini()`
entry in `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R`. The table
below collects them in one place for review.

| Equation / parameter | Value | Source location |
|----------------------|-------|-----------------|
| `<param1>`           | `<value>` | <Table / equation / page / figure of the source> |
| `<param2>`           | `<value>` | <...> |
| `<equation: d/dt(central)>` | n/a | <Equation 3, page X> |
| ...                  | ...     | ... |

# Virtual cohort

Original observed data are not publicly available. The figures below use virtual
populations whose covariate distributions approximate the published trial
demographics.

```{r cohort}
set.seed(<seed>)

# Helper: build one cohort as a self-contained event table. `id_offset`
# shifts subject IDs so multiple cohorts can be bind_rows()-ed without
# colliding. rxSolve treats ID as the subject key; duplicate IDs across
# cohorts silently collapse into single (wrong) subjects, so offsetting
# is mandatory, not decorative, for any multi-cohort simulation.
make_cohort <- function(n, ..., id_offset = 0L) {
  tibble(
    id = id_offset + seq_len(n),
    # ... covariate columns by canonical name (see inst/references/covariate-columns.md)
    # ... cohort/treatment/regimen label as a column so it can ride through rxSolve via `keep = `
  ) |>
    # expand into dosing + observation rows (evid, amt, cmt, time, ...)
    ...
}

# Single-cohort case:
events <- make_cohort(n = <integer>, ...)

# Multi-cohort case — pass distinct id_offset per call so IDs are disjoint:
# events <- dplyr::bind_rows(
#   make_cohort(200, ..., id_offset =   0L) |> mutate(cohort = "A"),
#   make_cohort(200, ..., id_offset = 200L) |> mutate(cohort = "B"),
#   make_cohort(200, ..., id_offset = 400L) |> mutate(cohort = "C")
# )
# stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

# Simulation

```{r simulate}
mod <- readModelDb("<FirstAuthor>_<Year>_<drug>")
# Prefer `keep = c("col1", "col2")` to carry source columns (cohort,
# treatment, dose group, regimen) through to the simulation output.
# This is cleaner than a post-hoc left_join back from `events`, and
# avoids row-alignment surprises when rxSolve drops or expands rows.
sim <- rxode2::rxSolve(mod, events = events, keep = c("cohort"))
```

For deterministic replication (reproducing Figure N without between-subject
variability), zero out the random effects:

```{r simulate-typical, eval = FALSE}
mod_typical <- mod |> rxode2::zeroRe()
sim_typical <- rxode2::rxSolve(mod_typical, events = events)
```

# Replicate published figures

Name each figure block with the source figure number so reviewers can pair them
up at a glance.

```{r figure-4}
# Replicates Figure 4 of <Author Year>: VPC of Cc vs. time by dose group.
sim |>
  group_by(time, treatment) |>
  summarise(
    Q05 = quantile(Cc, 0.05, na.rm = TRUE),
    Q50 = quantile(Cc, 0.50, na.rm = TRUE),
    Q95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  ggplot(aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line() +
  facet_wrap(~treatment) +
  scale_y_log10() +
  labs(x = "Time (<unit>)", y = "Cc (<unit>)",
       title = "Figure 4 — VPC by dose group",
       caption = "Replicates Figure 4 of <Author Year>.")
```

# PKNCA validation

Use PKNCA for Cmax, Tmax, AUC, and half-life. **Always include a treatment
grouping variable** (dose group, cohort, regimen) so the results can be
compared against any per-group values reported in the source paper. See
`references/pknca-recipes.md` for recipes covering single-dose,
steady-state, multiple-dose, and sparse-sampling cases.

```{r pknca}
# Concentrations — keep the column named Cc (nlmixr2lib convention).
# IMPORTANT: do NOT add `time > 0` or `Cc > 0` to this filter — both drop the
# time-zero row that PKNCA needs to anchor AUC0-*. Use only `!is.na(Cc)`.
sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, treatment)

# Guarantee a time=0 row per (id, treatment); for extravascular pre-dose Cc=0
# is the correct value. (See `pknca-recipes.md` § "Time-zero guarantee".)
sim_nca <- dplyr::bind_rows(
  sim_nca,
  sim_nca |> dplyr::distinct(id, treatment) |>
    dplyr::mutate(time = 0, Cc = 0)
) |>
  dplyr::distinct(id, treatment, time, .keep_all = TRUE) |>
  dplyr::arrange(id, treatment, time)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id)

# Doses — one row per dose event per subject; keep the column named amt
dose_df <- events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, treatment)

dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id)

# Intervals: what parameters to compute, over what time window
intervals <- data.frame(
  start       = 0,
  end         = Inf,
  cmax        = TRUE,
  tmax        = TRUE,
  aucinf.obs  = TRUE,
  half.life   = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
```

## Comparison against published NCA

When the source paper reports NCA values (Cmax, Tmax, AUC, half-life, …),
render a single combined side-by-side table using
`nlmixr2lib::ncaComparisonTable()`. Do **not** split simulated and reference
values across separate tables or refer to "see above"; the comparison only
works if the reader can see both numbers and the discrepancy in one row.

```{r nca_compare}
# Transcribe the paper's table into a wide tibble, one row per group, one
# column per PKNCA parameter code (cmax, tmax, aucinf.obs, auclast,
# half.life, …).
published <- tibble::tribble(
  ~treatment, ~cmax,  ~tmax, ~aucinf.obs, ~half.life,
  "50 mg",    14.8,   2.0,   125.0,       6.5,
  "100 mg",   28.5,   2.1,   250.0,       6.7
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated = nca_res,
  reference = published,
  by        = "treatment",
  units     = c(cmax = "ng/mL", aucinf.obs = "ng*h/mL",
                tmax = "h", half.life = "h"),
  tolerance_pct = 20
)

knitr::kable(
  cmp,
  caption = "Simulated vs. published NCA. * differs from reference by >20%.",
  align   = c("l", "l", "r", "r", "r")
)
```

Flag any starred rows in the narrative and investigate the source — do not
tune parameters to match.

# Assumptions and deviations

- <Any distributional assumption you had to make because the paper did not
  specify, e.g., "Race distribution set to X% White / Y% Black to match Table 1."
- <Any simplification, e.g., "Time-varying weight held constant at the baseline
  z-score over the follow-up window."
- <Any parameter the paper did not publish, and what was used instead.>
- **Non-paper-derived parameter values** — call out any value that came from a
  source other than the paper's text or tables. Cite the followup-register
  entry (`tracking/operator_followups.md F<n>`) when applicable. Examples:
  - "`kdeg = 0.0231 1/day` was supplied by the corresponding author
    (J. Almquist, email 2026-04-29; see `tracking/operator_followups.md F12`);
    the paper's Wiley supplement holds the value but is subscription-blocked."
  - "`V_DXd = 0.038 L/kg` was extracted by the operator from Figure 2 of
    <Author Year> via on-screen digitisation; ±10% uncertainty estimated from
    the figure's gridline resolution."
  - "PK structural parameters (`CL`, `Vc`, `Q`, `Vp`) were carried from the
    upstream model `<Upstream_Year_drug>` (this paper fixes its PK from that
    publication and reports only PD parameters)."
````

## Notes

- PKNCA formulas must include a treatment grouping (`| treatment + id`, or
  `| cohort + id`, `| regimen + id` — whatever stratification the source paper
  uses). The treatment grouping variable goes **before** `id` so summaries roll
  up per treatment. Omitting the grouping collapses all subjects into a single
  group and defeats the comparison against per-group published values.
- Keep the column named `Cc` (observation variable in nlmixr2 models), not
  `conc`. Keep dose named `amt`, not `dose`. This matches the rest of the
  nlmixr2lib / rxode2 pipeline.
- For endogenous / turnover models (where NCA isn't the right validation),
  replace the PKNCA section with the four validation patterns in
  `references/endogenous-validation.md` (steady-state hold,
  perturbation recovery, mass-balance / flux check, dimensional
  analysis), and use the endogenous vignette outline there.
- For multi-output models (e.g., parent + metabolite, or plasma + tissue),
  run one PKNCA block per output, each with its own `conc ~ time | id/treatment`.
- `rxode2::zeroRe()` is helpful when the published figure is a typical-value
  prediction rather than a VPC. Simulating with between-subject variability
  always, then plotting percentiles, matches VPC-style figures.
- **Multi-cohort simulations — disjoint IDs are mandatory.** If you build
  the event table by `bind_rows()`-ing several cohorts (dose groups,
  regimens, age strata, GA strata, trial panels, etc.), each cohort's `id`
  column must span a **disjoint integer range**. `rxSolve` treats `id` as
  the subject key; duplicate IDs across cohorts are silently merged into a
  single "Frankenstein" subject that receives the *summed* dose, producing
  predictions that can be 2-4× too high (the bug is silent — no warning).
  The `make_cohort(..., id_offset = )` helper in the cohort chunk above
  shows the pattern. Worked examples: `Robbie_2012_palivizumab.Rmd`,
  `Clegg_2024_nirsevimab.Rmd` (Figure 4), `Lin_2024_casirivimab.Rmd`,
  `Kuchimanchi_2018_evolocumab.Rmd`. Always add the assertion

  ```r
  stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
  ```

  immediately after the `bind_rows()` as a cheap regression guard.

- **Always carry grouping labels through `rxSolve(..., keep = ...)` —
  never with a post-hoc `left_join`.** `keep` attaches source columns
  directly in the rxode2 output aligned per row. The post-hoc pattern

  ```r
  # AVOID — fans out rows when an id has multiple labels (multi-cohort
  # bug surface), and even when ids are unique it's a footgun that
  # silently breaks if a future refactor introduces collisions.
  sim <- rxode2::rxSolve(mod, events = events) |>
    dplyr::left_join(events |> dplyr::select(id, treatment) |>
                       dplyr::distinct(), by = "id")
  ```

  collapses to wrong group totals the moment any `id` is repeated across
  cohort labels (e.g. the Clegg 2024 nirsevimab Figure 4 bug, where the
  `distinct()` step produced four trial-rows per subject and inflated
  every panel's median ~3-fold). Use `keep` instead:

  ```r
  # PREFER
  sim <- rxode2::rxSolve(mod, events = events,
                         keep = c("treatment", "cohort", "WT")) |>
    as.data.frame()
  ```

  Use `keep` for `cohort`, `treatment`, `regimen`, `trial`, dose-group
  labels, any per-subject covariate the plots or PKNCA formulas need
  downstream. The column-name is up to you (`"cohort"`, `"treatment"`,
  `"regimen"`, `"trial"`, etc. — whatever the paper's figure legends use);
  the skill does not mandate one. Character columns are preserved — they
  may come back as factor or character but the value is intact.

  The only legitimate post-rxSolve `left_join` use cases are:
  - joining NCA output (e.g. PKNCA results, one row per id) onto a
    per-subject summary table — that's not a per-time-point fan-out;
  - attaching reference values (published Cmax / AUC) by a label that
    was already carried through `keep`.

  When you see a `left_join(events |> select(id, X) |> distinct())` in
  any vignette, that's a bug surface to remove, not a pattern to copy.
