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
baseline weight band), join the group metadata after summarizing:

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

When the source paper reports NCA values (e.g., geometric mean with 95% CI), the
vignette should render a side-by-side table. A simple pattern:

```r
published <- tibble::tibble(
  treatment   = c("50 mg", "100 mg"),
  Cmax_pub    = c(<value>, <value>),
  AUCinf_pub  = c(<value>, <value>)
)

simulated <- res_tbl |>
  dplyr::filter(PPTESTCD %in% c("cmax", "aucinf.obs")) |>
  dplyr::group_by(treatment, PPTESTCD) |>
  dplyr::summarise(value = median(PPORRES), .groups = "drop") |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = value)

comparison <- published |> dplyr::left_join(simulated, by = "treatment")
knitr::kable(comparison)
```

Flag any differences > 20% in the narrative; do not tune parameters to match.

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
- **Time-zero records** — PKNCA expects a pre-dose record for IV or a `time = 0`
  observation for extravascular. Ensure the simulation grid includes `time = 0`.
- **Duplicate IDs across cohorts** — when `events` is assembled from multiple
  `make_cohort()` calls and `bind_rows`-ed together, confirm
  `anyDuplicated(sim[, c("id", "time")]) == 0` before handing `sim` to PKNCA.
  `rxSolve` silently collapses duplicated-ID rows into a single subject;
  PKNCA then aggregates a single (wrong) subject as if it were the whole
  group. Use the `id_offset` pattern in the `make_cohort` snippet in
  `references/vignette-template.md` to keep ID ranges disjoint.
- **Carry grouping via `rxSolve(..., keep = ...)`, not a post-hoc
  `left_join`.** `rxSolve` accepts a `keep = c("col1", "col2")` argument
  that attaches source columns (cohort, treatment, dose group, regimen)
  directly to the simulation output. This is aligned per row and far
  cleaner than joining back from `events` — a post-hoc `left_join` will
  produce multiplied rows if any IDs collide or if rxSolve expanded
  observation times. Round-trip `treatment` / `cohort` through `keep` and
  PKNCA's grouping formula picks it up automatically.
