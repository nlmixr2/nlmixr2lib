# Validation vignette template

File path: `vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd` — the pkgdown "articles" directory, not top-level `vignettes/`. Drug-specific vignettes live under `vignettes/articles/` so pkgdown builds them for the site but CRAN skips them (`.Rbuildignore` excludes the directory). The top-level `vignettes/` directory is reserved for the one legacy vignette (`PK_2cmt_mAb_Davda_2014.Rmd`). The basename must match the `vignette <- "..."` field in the model file.

A validation vignette has two jobs:

1. **Verification trail** — document exactly where each model equation and parameter came from.
2. **Validation** — show that a simulation from the packaged model reproduces published results (usually one or more figures and an NCA table) using **PKNCA** for NCA parameters.

Use this template as a starting point. Every section is required unless marked optional.

## Template

````markdown
---
title: "<FirstAuthor>_<Year>_<drug>"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{<FirstAuthor>_<Year>_<drug>}
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
# Concentrations — keep the column named Cc (nlmixr2lib convention)
sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, treatment)

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

nca_summary <- summary(nca_res)
knitr::kable(nca_summary, caption = "Simulated NCA parameters by dose group.")
```

## Comparison against published NCA

<If the source paper reports NCA parameters (e.g., mean Cmax and AUC by dose
group in a Table), reproduce that table side-by-side with the simulated values
here. Note any differences > 20% and investigate — do not tune parameters to
match.>

# Assumptions and deviations

- <Any distributional assumption you had to make because the paper did not
  specify, e.g., "Race distribution set to X% White / Y% Black to match Table 1."
- <Any simplification, e.g., "Time-varying weight held constant at the baseline
  z-score over the follow-up window."
- <Any parameter the paper did not publish, and what was used instead.>
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
  replace the PKNCA section with:
  - baseline recovery (simulate with drug removed and confirm return to
    `css` within the reported time)
  - turnover / steady-state check (integrate long enough to verify the
    reported steady-state concentration)
- For multi-output models (e.g., parent + metabolite, or plasma + tissue),
  run one PKNCA block per output, each with its own `conc ~ time | id/treatment`.
- `rxode2::zeroRe()` is helpful when the published figure is a typical-value
  prediction rather than a VPC. Simulating with between-subject variability
  always, then plotting percentiles, matches VPC-style figures.
- **Multi-cohort simulations.** If you build the event table by
  `bind_rows()`-ing several cohorts (dose groups, regimens, age strata, etc.),
  each cohort's `id` column must span a disjoint integer range. The
  `make_cohort(..., id_offset = )` helper in the cohort chunk above shows the
  pattern. Duplicate IDs across cohorts are silently merged by `rxSolve` into
  a single subject and are a recurring bug — the `Robbie_2012_palivizumab`
  vignette has a worked example. Add the assertion
  `stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))`
  after the bind_rows as a cheap regression guard.
- **Prefer `rxSolve(..., keep = c("col1", "col2"))`** over a post-hoc
  `left_join` from `events` back into the simulation output. `keep` attaches
  source columns directly in the rxode2 output, aligned per row. Use it for
  `cohort`, `treatment`, dose-group labels, and any grouping column the plots
  or PKNCA formulas need. The column-name is up to you (`"cohort"`,
  `"treatment"`, `"regimen"` — whatever the paper's figure legends use); the
  skill does not mandate one.
