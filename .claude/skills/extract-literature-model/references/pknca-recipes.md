# PKNCA recipes for validation vignettes

Every validation vignette uses PKNCA for NCA parameters (Cmax, Tmax, AUC, half-life, …) rather than inline trapezoidal calculations. All recipes include a **treatment grouping variable** in the formula so results can be compared against per-group values reported in the source paper.

PKNCA reference: `?PKNCA::PKNCA` and `vignette("Introduction-and-Usage", package = "PKNCA")`.

## Data shape required

- **Concentration data:** one row per subject × time. Columns: `id`, `time`, `Cc` (the simulated concentration from the model — keep the column named `Cc` to match nlmixr2lib conventions), plus the grouping column (`treatment`, `cohort`, `regimen`, or similar).
- **Dose data:** one row per dose event. Columns: `id`, `time`, `amt` (dose amount — keep the rxode2/nlmixr2 column name), plus the same grouping column.

Both frames must agree on `id` and the grouping column.

The formula is always `Cc ~ time | treatment + id` (and `amt ~ time | treatment + id` for dose). The **treatment grouping variable goes first**, before `id`, so PKNCA summaries roll up per treatment as reported in most source papers.

## Recipe 1 — Single-dose, dense sampling (Cmax, Tmax, AUC0-inf, half-life)

Use when the paper reports NCA after a single dose with enough sampling to characterize the terminal phase.

```r
library(PKNCA)

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, treatment)

dose_df <- events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, treatment)

# Always declare units: `concu` + `timeu` on PKNCAconc, `doseu` on PKNCAdose.
# PKNCA uses these to do automatic unit-consistent AUC / clearance calculations;
# without them the summary tables carry only numbers and downstream consumers
# can't tell ng/mL from ug/mL. Use the same units that appear in the model
# file's `units` metadata and the paper's Table footnotes.
conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id,
                             concu = "<conc unit, e.g. ug/mL>",
                             timeu = "<time unit, e.g. day>")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id,
                             doseu = "<dose unit, e.g. mg>")

intervals <- data.frame(
  start      = 0,
  end        = Inf,
  cmax       = TRUE,
  tmax       = TRUE,
  aucinf.obs = TRUE,
  aucinf.pred = TRUE,
  half.life  = TRUE,
  clast.obs  = TRUE,
  lambda.z   = TRUE
)

res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
summary(res)
```

## Recipe 2 — Single-dose, AUClast only (sparse terminal data)

Use when terminal sampling is too sparse to estimate `lambda.z` reliably.

```r
intervals <- data.frame(
  start     = 0,
  end       = max(sim_nca$time),
  cmax      = TRUE,
  tmax      = TRUE,
  auclast   = TRUE,
  clast.obs = TRUE
)
res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
```

## Recipe 3 — Steady state (AUC0-tau, Cmax,ss, Cmin,ss, Cavg,ss)

Use when the paper reports steady-state NCA over a dosing interval `tau`.

```r
tau <- <dosing interval, e.g., 24>  # same units as time

# Pick the AUC0-tau interval at steady state, e.g., the last dosing interval
start_ss <- max(dose_df$time)          # time of the final dose
end_ss   <- start_ss + tau

intervals <- data.frame(
  start    = start_ss,
  end      = end_ss,
  cmax     = TRUE,
  tmax     = TRUE,
  cmin     = TRUE,
  auclast  = TRUE,
  cav      = TRUE,   # average concentration over the interval
  ctau     = TRUE    # concentration at end of interval
)

res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
```

## Recipe 4 — Multiple-dose with accumulation

Compute AUC0-tau on the first dosing interval, AUC0-tau at steady state, and the accumulation ratio.

```r
intervals <- data.frame(
  start = c(0, start_ss),
  end   = c(tau, end_ss),
  cmax     = TRUE,
  auclast  = TRUE
)

res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
# Accumulation ratio = AUC0-tau at SS / AUC0-tau after first dose
```

## Recipe 5 — Per-subgroup summaries beyond the grouping variable

PKNCA's `summary()` collapses within the group levels defined by the formula. If
the source paper reports NCA stratified by an additional covariate (e.g., by
baseline weight band), the cleanest path is to **carry the stratifier through
`rxSolve(..., keep = ...)`** so it lands directly in the simulation output and
in the PKNCA grouping formula — no post-summarise join is needed:

```r
sim <- rxode2::rxSolve(mod, events = events,
                       keep = c("treatment", "weight_band"))

conc_obj <- PKNCA::PKNCAconc(
  data = as.data.frame(sim),
  formula = Cc ~ time | treatment + weight_band + id,
  concu = "ug/mL", timeu = "day"
)
```

If you cannot route the stratifier through `keep` (e.g. it was derived inside
the simulation loop and only exists in `cohort`), join it onto the **PKNCA
result table** (one row per subject per parameter), not onto the per-time-point
simulation rows — the result table is a 1:1 join by `id` and never fans out:

```r
res_tbl <- as.data.frame(res$result)

res_joined <- res_tbl |>
  dplyr::left_join(cohort |> dplyr::select(id, weight_band), by = "id") |>
  dplyr::group_by(treatment, weight_band, PPTESTCD) |>
  dplyr::summarise(
    median_value = median(PPORRES),
    q05          = quantile(PPORRES, 0.05),
    q95          = quantile(PPORRES, 0.95),
    .groups      = "drop"
  )
```

## Comparing against the published table

When the source paper reports NCA values, render a **single combined**
side-by-side table using `nlmixr2lib::ncaComparisonTable()`. The helper
accepts a `PKNCAresults` object (or its `$result` data frame) for the
simulated side, and a wide tibble (one row per group, one column per
PKNCA code) for the reference side — convenient when transcribing values
from a paper table.

```r
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

Output: one row per parameter × group, with columns `NCA parameter`, the
grouping variable(s), `Reference`, `Simulated`, and `% diff`. Rows whose
discrepancy exceeds `tolerance_pct` get a trailing `*`, and the function
attaches a footnote string at `attr(cmp, "footnote")`. Friendly parameter
labels (`Cmax`, `Tmax`, `AUC0-∞ (obs)`, `t½`, …) come from
[`nlmixr2lib::ncaParamLabel()`] — never render `PPTESTCD` codes
(`cmax`, `auclast`, …) directly to the reader.

**Never split the comparison into separate simulated / reference tables**
or use "see above" cross-references — the reader can only verify the
comparison when both values and the discrepancy sit in the same row.
`ncaComparisonTable()` enforces this shape automatically.

Flag any starred rows in the narrative; do not tune parameters to match.

## Common pitfalls

- **Missing treatment grouping** — if the formula is `Cc ~ time | id` with no
  treatment, PKNCA aggregates across dose groups and the Cmax / AUC results
  are uninterpretable. Always include the treatment grouping, and put it
  **before** `id` (`Cc ~ time | treatment + id`).
- **Missing unit declarations** — pass `concu` + `timeu` to `PKNCAconc()` and
  `doseu` to `PKNCAdose()`. Without them, PKNCA's output is unit-blind: AUC
  values are just numbers, clearance is not derivable, and the summary
  tables cannot be cross-checked against the paper's own NCA values (which
  always report units). Match the units to the model file's `units`
  metadata (`time = "day"`, `dosing = "mg"`, `concentration = "ug/mL"` —
  whatever the paper uses). Note that the units are strings, not R
  objects — e.g. `concu = "ug/mL"`, not `concu = units::as_units("ug/mL")`.
- **Renaming `Cc` to `conc`** — keep the column named `Cc` (same as the
  observation variable in the nlmixr2 model) rather than renaming it to
  `conc`. Same for dose: keep `amt`, not `dose`.
- **Dose units ≠ concentration units** — PKNCA doesn't check. Confirm dose is
  in the same mass unit as `Cc × volume` implied by the model (e.g., mg vs.
  ng × L).
- **`lambda.z` warnings** — PKNCA will emit warnings when there aren't enough
  post-peak points. That's informational; if the paper used a specific
  regression window, set `lambda.z.time.first` and `lambda.z.n.points` in the
  `intervals` data frame.
- **BLQ handling** — the model emits continuous concentrations (no BLQ), so
  BLQ handling in PKNCA is usually a no-op. If the paper applied a specific
  BLQ rule (e.g., M3 method), that's outside the NCA step.
- **Time-zero records (mandatory)** — PKNCA's AUC integration starts at the
  interval's `start` (`= 0` by default), and emits the warning
  `Requesting an AUC range starting (0) before the first measurement
  (<t1>) is not allowed` — repeated once per subject — when the
  concentration frame passed to it has no row at `time = 0`. Two failure
  modes recur:

  1. **Over-eager filters.** The PKNCA input chunk uses
     `dplyr::filter(time > 0, !is.na(Cc), Cc > 0)` — borrowed from a
     log-scale plotting filter — which drops the time-zero row outright.
     The Archary 2019 abacavir vignette exhibited this; the fix is to use
     **only** `!is.na(Cc)` in the PKNCA input filter.
  2. **Sim grid never produced a time-zero observation.** Some
     `make_cohort()` patterns omit `time = 0` from the observation grid
     (especially when sampling is sparse), so even an unfiltered output
     lacks the row.

  The defensive pattern guarantees a row regardless. For extravascular
  models, pre-dose `Cc = 0` is the correct value; for IV bolus models,
  set `Cc = dose / Vc` (or accept the back-extrapolated value PKNCA
  computes when `lambda.z` is fit):

  ```r
  sim_nca <- sim |>
    dplyr::filter(!is.na(Cc)) |>
    dplyr::select(id, time, Cc, treatment)

  # Add time=0 with Cc=0 (extravascular); existing time=0 rows win via
  # .keep_all = TRUE on the first occurrence.
  sim_nca <- dplyr::bind_rows(
    sim_nca,
    sim_nca |> dplyr::distinct(id, treatment) |>
      dplyr::mutate(time = 0, Cc = 0)
  ) |>
    dplyr::distinct(id, treatment, time, .keep_all = TRUE) |>
    dplyr::arrange(id, treatment, time)
  ```

  Verify the rendered vignette is free of the
  `Requesting an AUC range starting (0)` warning before pushing — it is
  the single most common PKNCA defect in this repo.
- **Cohort-ID safety and grouping via `keep =`** — when `events` is built from
  multiple `make_cohort()` calls and `bind_rows`-ed, and when the PKNCA
  formula needs a `treatment` / `cohort` / `regimen` label, follow the
  patterns in `references/vignette-template.md` § "Notes" (`id_offset`
  for disjoint IDs; `rxSolve(..., keep = ...)` instead of post-hoc
  `left_join`). Both errors silently produce wrong NCA values — the
  Clegg 2024 nirsevimab Figure 4 bug (~3-fold inflation) was a
  `left_join` after an id collision.
