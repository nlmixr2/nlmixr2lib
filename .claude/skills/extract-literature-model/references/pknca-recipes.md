# PKNCA recipes for validation vignettes

Every validation vignette uses PKNCA for NCA parameters (Cmax, Tmax, AUC, half-life, …) rather than inline trapezoidal calculations. All recipes include a **treatment grouping variable** in the formula so results can be compared against per-group values reported in the source paper.

PKNCA reference: `?PKNCA::PKNCA` and `vignette("Introduction-and-Usage", package = "PKNCA")`.

## Data shape required

- **Concentration data:** one row per subject × time. Columns: `id`, `time`, `conc` (the `Cc` output from simulation), plus the grouping column (`treatment`, `cohort`, `regimen`, or similar).
- **Dose data:** one row per dose event. Columns: `id`, `time`, `dose` (amount), plus the same grouping column.

Both frames must agree on `id` and the grouping column.

## Recipe 1 — Single-dose, dense sampling (Cmax, Tmax, AUC0-inf, half-life)

Use when the paper reports NCA after a single dose with enough sampling to characterize the terminal phase.

```r
library(PKNCA)

conc_df <- sim |>
  filter(!is.na(Cc)) |>
  transmute(id, time, conc = Cc, treatment)

dose_df <- events |>
  filter(evid == 1) |>
  transmute(id, time, dose = amt, treatment)

conc_obj <- PKNCAconc(conc_df, conc ~ time | id / treatment)
dose_obj <- PKNCAdose(dose_df, dose ~ time | id / treatment)

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

res <- pk.nca(PKNCAdata(conc_obj, dose_obj, intervals = intervals))
summary(res)
```

## Recipe 2 — Single-dose, AUClast only (sparse terminal data)

Use when terminal sampling is too sparse to estimate `lambda.z` reliably.

```r
intervals <- data.frame(
  start     = 0,
  end       = max(conc_df$time),
  cmax      = TRUE,
  tmax      = TRUE,
  auclast   = TRUE,
  clast.obs = TRUE
)
res <- pk.nca(PKNCAdata(conc_obj, dose_obj, intervals = intervals))
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

res <- pk.nca(PKNCAdata(conc_obj, dose_obj, intervals = intervals))
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

res <- pk.nca(PKNCAdata(conc_obj, dose_obj, intervals = intervals))
# Accumulation ratio = AUC0-tau at SS / AUC0-tau after first dose
```

## Recipe 5 — Per-subgroup summaries beyond the grouping variable

PKNCA's `summary()` collapses within the group levels defined by the formula. If
the source paper reports NCA stratified by an additional covariate (e.g., by
baseline weight band), join the group metadata after summarizing:

```r
res_tbl <- as.data.frame(res$result)

res_joined <- res_tbl |>
  left_join(cohort |> select(id, weight_band), by = "id") |>
  group_by(treatment, weight_band, PPTESTCD) |>
  summarise(
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
published <- tibble(
  treatment   = c("50 mg", "100 mg"),
  Cmax_pub    = c(<value>, <value>),
  AUCinf_pub  = c(<value>, <value>)
)

simulated <- res_tbl |>
  filter(PPTESTCD %in% c("cmax", "aucinf.obs")) |>
  group_by(treatment, PPTESTCD) |>
  summarise(value = median(PPORRES), .groups = "drop") |>
  pivot_wider(names_from = PPTESTCD, values_from = value)

comparison <- published |> left_join(simulated, by = "treatment")
knitr::kable(comparison)
```

Flag any differences > 20% in the narrative; do not tune parameters to match.

## Common pitfalls

- **Missing treatment grouping** — if the formula is `conc ~ time | id` with no
  `/ treatment`, PKNCA aggregates across dose groups and the Cmax / AUC results
  are uninterpretable. Always include the grouping.
- **Dose units ≠ concentration units** — PKNCA doesn't check. Confirm dose is
  in the same mass unit as `conc × volume` implied by the model (e.g., mg vs.
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
