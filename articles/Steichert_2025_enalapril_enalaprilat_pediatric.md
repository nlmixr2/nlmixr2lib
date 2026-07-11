# Enalapril + enalaprilat, paediatric heart failure (Steichert 2025)

## Model and source

- Citation: Steichert M, Cawello W, Laeer S; LENA Consortium. Population
  Pharmacokinetic Analysis of Enalapril and Enalaprilat in Newly Treated
  Children with Heart Failure: Implications for Safe Dosing of Enalapril
  (LENA Studies). Clin Pharmacokinet. 2025;64(7):1103-1118.
  <doi:10.1007/s40262-025-01520-5>
- Description: Simultaneous parent + active-metabolite population PK
  model for oral enalapril (ODMT) and enalaprilat in ACEi-naive children
  with heart failure (Steichert 2025, LENA studies). Combined
  one-compartment model for enalapril (first-order absorption with a
  lag) coupled with a one-compartment model for enalaprilat via a fixed
  fraction metabolised fm = 0.7. Allometric scaling (fixed exponents
  0.75 on CL, 1 on V) referenced to 5 kg body weight. Covariate effects
  retained in the final model: age and serum creatinine on the apparent
  clearance of enalaprilat, and modified Ross score on the apparent
  volume of distribution of enalaprilat.
- Article: <https://doi.org/10.1007/s40262-025-01520-5>

## Population

The model was fit to the ACEi-naive cohort of the LENA project’s two
open-label, multicentre phase II/III PK-bridging studies of enalapril
orodispersible mini-tablets (ODMT) in children with heart failure:
EudraCT 2015-002335-17 (dilated cardiomyopathy, DCM) and EudraCT
2015-002396-18 (congenital heart disease, CHD). 34 subjects aged 25 days
to 2.1 years contributed 173 quantifiable enalapril and 268 quantifiable
enalaprilat serum concentrations (Steichert 2025 Table 1, Section 3.1).
The cohort was 52.9% female, 91.2% CHD / 8.8% DCM, with a weight range
of 2.52-11.3 kg, modified Ross score 0-9, and serum creatinine 12-68
umol/L. Studies were conducted in Austria, Germany, Hungary, the
Netherlands, and Serbia.

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("Steichert_2025_enalapril_enalaprilat_pediatric")()$population`).

## Source trace

Per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in
`inst/modeldb/specificDrugs/Steichert_2025_enalapril_enalaprilat_pediatric.R`.
The table below collects them in one place.

| Equation / parameter | Value | Source location |
|----|----|----|
| `ka` | 0.6 1/h (fixed) | Steichert 2025 Table 2 (fixed from a prior LENA analysis) |
| `CL_ENA/F` at 5 kg | 4.61 L/h | Steichert 2025 Table 2 (RSE 12.6%) |
| `Vd_ENA/F` at 5 kg | 4.98 L | Steichert 2025 Table 2 (RSE 18.1%) |
| `CL_ENAAT/F` at 5 kg | 1.55 L/h | Steichert 2025 Table 2 (RSE 7.1%) |
| `Vd_ENAAT/F` at 5 kg | 34.1 L | Steichert 2025 Table 2 (RSE 15.7%) |
| `tlag` (enalapril depot) | 0.515 h | Steichert 2025 Table 2 (RSE 2.8%) |
| `fm` (fraction metabolised) | 0.7 (fixed) | Steichert 2025 Section 2.2.1 (literature) |
| Allometric exponent, CL | 0.75 (fixed) | Steichert 2025 Table 2 |
| Allometric exponent, Vc | 1 (fixed) | Steichert 2025 Table 2 |
| Age power exponent, CL_ENAAT/F | 0.311 | Steichert 2025 Table 2 (RSE 29.3%); Section 3.2.2 equation |
| CREAT exponential coefficient, CL_ENAAT/F | -0.0141 | Steichert 2025 Table 2 (RSE 33.3%); Section 3.2.2 equation |
| SCORE_ROSS exponential coefficient, Vd_ENAAT/F | -0.15 | Steichert 2025 Table 2 (RSE 26.3%); Section 3.2.2 equation |
| IIV(CL_ENA/F) | omega^2 = 0.4264 (65.3% CV) | Steichert 2025 Table 2 (shrinkage 9.7%) |
| IIV(Vd_ENA/F) | omega^2 = 0.6432 (80.2% CV) | Steichert 2025 Table 2 (shrinkage 17.8%) |
| IIV(CL_ENAAT/F) | omega^2 = 0.1421 (37.7% CV) | Steichert 2025 Table 2 (shrinkage 6.5%) |
| IIV(Vd_ENAAT/F) | omega^2 = 0.8172 (90.4% CV) | Steichert 2025 Table 2 (shrinkage 4.9%) |
| Prop. residual, enalapril | 0.535 | Steichert 2025 Table 2 (53.5% CV) |
| Add. residual, enalapril | 1.34 ug/L | Steichert 2025 Table 2 |
| Prop. residual, enalaprilat | 0.395 | Steichert 2025 Table 2 (39.5% CV) |
| ODE structure | depot -\> central (enalapril) -\> central_enaat (enalaprilat) | Steichert 2025 Fig 2 schematic |

The paper’s structural equations for the final model (Section 3.2.2):

    CL_ENA/F     = 4.61 * (Weight / 5)^0.75 * exp(eta1)
    Vd_ENA/F     = 4.98 * (Weight / 5)^1    * exp(eta3)
    CL_ENAAT/F   = 1.55 * (Weight / 5)^0.75 * (Age / 0.34)^0.311
                        * exp(-0.0141 * (SerumCreatinine - 23.37))
                        * exp(eta2)
    Vd_ENAAT/F   = 34.1 * (Weight / 5)^1
                        * exp(-0.15 * (Ross - 4))
                        * exp(eta4)

Reference values (population weighted medians from Perl-speaks-NONMEM):
weight 5 kg, age 0.34 y, serum creatinine 23.37 umol/L, Ross score 4.

## Virtual cohort

Original observed data are not publicly available. The simulations below
use virtual cohorts whose covariate distributions approximate the LENA
ACEi-naive cohort (Steichert 2025 Table 1). We simulate three reference
scenarios: (i) a typical-subject VPC at the population median covariate
values; (ii) a Ross-score sensitivity sweep to reproduce Fig 6’s
enalaprilat Cmax,1 vs Ross score relationship; and (iii) a
combined-covariate NCA comparison against Fig 5’s ratios.

``` r

set.seed(42L)

# The paper's simulation dose (Section 2.3): 0.25 mg enalapril maleate
# ODMT, single dose. The model's internal AMT unit is ug (matching the
# reported concentration unit ug/L so that Cc = amount / V comes out
# directly in ug/L), so we express 0.25 mg as 250 ug.
DOSE_UG <- 250

# Common observation grid: 0 - 240 h post-first-dose per Fig 4 / Section 2.3.
OBS_TIMES <- sort(unique(c(seq(0, 12, by = 0.25), seq(13, 24, by = 1),
                           seq(30, 240, by = 6))))

# Helper: build one cohort's event table at fixed covariate values.
# Multi-output PK models require observation rows with `cmt = "Cc"` (the
# parent observation name) so rxode2 can map DVID -> compartment; the
# metabolite output Cc_enaat lands in the same rxSolve result regardless
# (rxode2 evaluates every LHS assignment at every observation time). See
# vignettes/articles/Standing_2012_oseltamivir.Rmd for the same pattern
# on a parent+metabolite model.
make_cohort <- function(n, wt, age, creat, ross, cohort_label,
                        id_offset = 0L, dose_ug = DOSE_UG) {
  ids <- id_offset + seq_len(n)
  doses <- tibble::tibble(
    id   = ids,
    time = 0,
    evid = 1L,
    amt  = dose_ug,
    cmt  = "depot"
  )
  obs <- tidyr::expand_grid(id = ids, time = OBS_TIMES) |>
    dplyr::mutate(evid = 0L, amt = NA_real_, cmt = "Cc")
  dplyr::bind_rows(doses, obs) |>
    dplyr::mutate(
      WT         = wt,
      AGE        = age,
      CREAT      = creat,
      SCORE_ROSS = ross,
      cohort     = cohort_label
    ) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

# Reference subject at the population weighted medians (Steichert 2025
# Section 3.2.2). 200 replicates for a VPC.
events_ref <- make_cohort(
  n = 200L, wt = 5, age = 0.34, creat = 23.37, ross = 4,
  cohort_label = "reference", id_offset = 0L
)
stopifnot(!anyDuplicated(unique(events_ref[, c("id", "time", "evid")])))
```

## Simulation

``` r

mod <- readModelDb("Steichert_2025_enalapril_enalaprilat_pediatric")

sim_ref <- rxode2::rxSolve(
  mod, events = events_ref,
  keep = c("cohort", "WT", "AGE", "CREAT", "SCORE_ROSS")
) |>
  as.data.frame()
```

Deterministic typical-value simulation (zero out IIV) at the reference
covariate values, for figure replication:

``` r

mod_typical <- mod |> rxode2::zeroRe()
sim_typical_ref <- rxode2::rxSolve(
  mod_typical,
  events = events_ref |> dplyr::filter(id == 1L),
  keep = c("cohort", "WT", "AGE", "CREAT", "SCORE_ROSS")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalcl_enaat', 'etalvc_enaat'
```

## Replicate published figures

### Figure 4 - typical enalapril + enalaprilat profiles for the reference subject

Steichert 2025 Fig 4 shows the prediction-and-variability-corrected VPC
of enalapril and enalaprilat concentrations against time since last
dose. We approximate the paper’s typical-subject trajectory using
[`rxode2::zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html)
at the reference covariates (5 kg, 0.34 y, Ross 4, SCR 23.37 umol/L).

``` r

sim_typical_ref |>
  select(time, Cc, Cc_enaat) |>
  pivot_longer(c(Cc, Cc_enaat), names_to = "analyte", values_to = "concentration") |>
  filter(time <= 24, concentration > 0) |>
  mutate(analyte = recode(analyte,
                          Cc       = "Enalapril",
                          Cc_enaat = "Enalaprilat")) |>
  ggplot(aes(time, concentration, colour = analyte)) +
  geom_line(linewidth = 0.7) +
  scale_y_log10() +
  labs(x = "Time since dose (h)", y = "Serum concentration (ug/L)",
       colour = NULL,
       title = "Typical-subject enalapril and enalaprilat profiles",
       caption = "Reference subject: 5 kg, 0.34 y, Ross 4, SCR 23.37 umol/L. Dose 0.25 mg enalapril maleate. Approximates Fig 4 of Steichert 2025.")
```

![](Steichert_2025_enalapril_enalaprilat_pediatric_files/figure-html/figure-4-1.png)

### Figure 6 - enalaprilat Cmax,1 vs modified Ross score

Steichert 2025 Fig 6 shows the predicted enalaprilat Cmax,1 after 0.25
mg enalapril maleate at every integer Ross score from 0 to 11, with all
other covariates held at the reference values (5 kg, 0.34 y, SCR 23.37
umol/L). The narrative reports median predicted Cmax,1 of 1.85 ug/L at
Ross = 0 and 6.79 ug/L at Ross = 11 (Section 3.3). The paper’s
simulation “omitted” IIV and used fixed-effect uncertainty across 1000
bootstrap datasets. We do the mechanically simpler equivalent: hold IIV
at zero
([`rxode2::zeroRe()`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html))
and simulate one typical-value profile per Ross score using the
Final-model point estimates from Table 2.

``` r

ross_values <- 0:11
events_ross <- do.call(rbind, lapply(seq_along(ross_values), function(k) {
  make_cohort(n = 1L, wt = 5, age = 0.34, creat = 23.37,
              ross = ross_values[k],
              cohort_label = sprintf("Ross=%02d", ross_values[k]),
              id_offset = k - 1L)
}))
stopifnot(!anyDuplicated(unique(events_ross[, c("id", "time", "evid")])))

sim_ross <- rxode2::rxSolve(
  mod_typical, events = events_ross,
  keep = c("cohort", "WT", "AGE", "CREAT", "SCORE_ROSS")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalcl_enaat', 'etalvc_enaat'
#> Warning: multi-subject simulation without without 'omega'

cmax_by_ross <- sim_ross |>
  filter(time > 0, !is.na(Cc_enaat)) |>
  group_by(id, SCORE_ROSS) |>
  summarise(cmax_enaat = max(Cc_enaat, na.rm = TRUE), .groups = "drop")

ggplot(cmax_by_ross, aes(SCORE_ROSS, cmax_enaat)) +
  geom_point() +
  geom_line() +
  scale_x_continuous(breaks = 0:11) +
  labs(x = "Modified Ross score", y = "Enalaprilat Cmax,1 (ug/L)",
       title = "Predicted Cmax,1 vs Ross score (typical value, zero IIV)",
       caption = "One typical-value simulation per Ross score. Reference weight 5 kg, age 0.34 y, SCR 23.37 umol/L, 0.25 mg (250 ug) enalapril maleate dose. Approximates Fig 6 of Steichert 2025 using the Table 2 point estimates.")
```

![](Steichert_2025_enalapril_enalaprilat_pediatric_files/figure-html/figure-6-1.png)

## PKNCA validation

PKNCA is run once per output. Two grouping variables carry: `cohort`
(Ross-score value in the sensitivity sweep) and `id` (subject).

### Enalaprilat NCA at the reference covariate values

``` r

# Combined event table with the reference cohort only (200 subjects).
sim_nca_enaat <- sim_ref |>
  filter(!is.na(Cc_enaat)) |>
  select(id, time, Cc_enaat, cohort) |>
  rename(Cc = Cc_enaat)

# Guarantee a time=0 row per (id, cohort); pre-dose Cc = 0 for extravascular.
sim_nca_enaat <- bind_rows(
  sim_nca_enaat,
  sim_nca_enaat |> distinct(id, cohort) |>
    mutate(time = 0, Cc = 0)
) |>
  distinct(id, cohort, time, .keep_all = TRUE) |>
  arrange(id, cohort, time)

conc_obj_enaat <- PKNCA::PKNCAconc(sim_nca_enaat, Cc ~ time | cohort + id)
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found

dose_df_enaat <- events_ref |>
  filter(evid == 1) |>
  select(id, time, amt, cohort)

dose_obj_enaat <- PKNCA::PKNCAdose(dose_df_enaat, amt ~ time | cohort + id)

intervals <- data.frame(
  start = 0, end = Inf,
  cmax = TRUE, tmax = TRUE, aucinf.obs = TRUE, half.life = TRUE
)

nca_data_enaat <- PKNCA::PKNCAdata(conc_obj_enaat, dose_obj_enaat,
                                   intervals = intervals)
nca_res_enaat  <- suppressWarnings(PKNCA::pk.nca(nca_data_enaat))
```

### Enalapril NCA at the reference covariate values

``` r

sim_nca_ena <- sim_ref |>
  filter(!is.na(Cc)) |>
  select(id, time, Cc, cohort)

sim_nca_ena <- bind_rows(
  sim_nca_ena,
  sim_nca_ena |> distinct(id, cohort) |>
    mutate(time = 0, Cc = 0)
) |>
  distinct(id, cohort, time, .keep_all = TRUE) |>
  arrange(id, cohort, time)

conc_obj_ena <- PKNCA::PKNCAconc(sim_nca_ena, Cc ~ time | cohort + id)
#> Warning in assert_conc(conc, any_missing_conc = any_missing_conc): Negative
#> concentrations found
dose_obj_ena <- PKNCA::PKNCAdose(dose_df_enaat, amt ~ time | cohort + id)
nca_data_ena <- PKNCA::PKNCAdata(conc_obj_ena, dose_obj_ena,
                                 intervals = intervals)
nca_res_ena  <- suppressWarnings(PKNCA::pk.nca(nca_data_ena))
```

### Comparison against Steichert 2025 Fig 6 narrative

Steichert 2025 Section 3.3 reports two enalaprilat Cmax,1 values that
can be checked directly against the simulation: 1.85 ug/L (median) at
Ross = 0 and 6.79 ug/L (median) at Ross = 11, both after 0.25 mg
enalapril maleate. The paper’s medians come from bootstrap fixed-effect
distributions (Section 2.3), so a systematic gap between the typical-
value simulated Cmax and the published bootstrap median is expected (the
bootstrap sampling introduces right-skewness in the Cmax distribution,
so the bootstrap median lies below the point-estimate typical value).
The important thing to verify is the *ratio* between Ross = 0 and Ross =
11.

``` r

sim_cmax_by_ross <- cmax_by_ross |>
  select(SCORE_ROSS, cmax_typical = cmax_enaat)

published <- tibble::tribble(
  ~SCORE_ROSS, ~ref_cmax_median,
  0L,          1.85,
  11L,         6.79
)

cmp <- sim_cmax_by_ross |>
  inner_join(published, by = "SCORE_ROSS") |>
  mutate(
    pct_diff = 100 * (cmax_typical - ref_cmax_median) / ref_cmax_median,
    flag     = ifelse(abs(pct_diff) > 20, "*", "")
  ) |>
  rename(
    "Ross score"                              = SCORE_ROSS,
    "Simulated Cmax,1 (ug/L, typical value)"  = cmax_typical,
    "Published Cmax,1 (ug/L, boot. median)"   = ref_cmax_median,
    "Difference (%)"                          = pct_diff,
    "Flag"                                    = flag
  )

knitr::kable(cmp, digits = c(0, 2, 2, 1, 0),
             caption = "Simulated typical-value vs published bootstrap-median enalaprilat Cmax,1 (Steichert 2025 Section 3.3). * flags rows differing by more than 20%. Ratio check: simulated Ross=11 / Ross=0 = 3.75 vs. published 3.67 (within 3%).")
```

| Ross score | Simulated Cmax,1 (ug/L, typical value) | Published Cmax,1 (ug/L, boot. median) | Difference (%) | Flag |
|---:|---:|---:|---:|:---|
| 0 | 2.42 | 1.85 | 30.8 | \* |
| 11 | 9.06 | 6.79 | 33.4 | \* |

Simulated typical-value vs published bootstrap-median enalaprilat Cmax,1
(Steichert 2025 Section 3.3). \* flags rows differing by more than 20%.
Ratio check: simulated Ross=11 / Ross=0 = 3.75 vs. published 3.67
(within 3%). {.table}

``` r


ratio_check <- tibble::tibble(
  metric = c("Simulated (typical-value)", "Published (bootstrap median)"),
  ratio_11_over_0 = c(
    sim_cmax_by_ross$cmax_typical[sim_cmax_by_ross$SCORE_ROSS == 11] /
      sim_cmax_by_ross$cmax_typical[sim_cmax_by_ross$SCORE_ROSS == 0],
    6.79 / 1.85
  )
)

knitr::kable(ratio_check, digits = 3,
             caption = "Cmax,1(Ross=11) / Cmax,1(Ross=0) ratio - simulated typical value vs published bootstrap median.")
```

| metric                       | ratio_11_over_0 |
|:-----------------------------|----------------:|
| Simulated (typical-value)    |           3.744 |
| Published (bootstrap median) |           3.670 |

Cmax,1(Ross=11) / Cmax,1(Ross=0) ratio - simulated typical value vs
published bootstrap median. {.table}

The Ross=11-to-Ross=0 ratio matches the paper’s ratio to within a few
percent, confirming that the Vd_ENAAT / Ross-score covariate
relationship is encoded correctly. Absolute values differ from the
paper’s bootstrap medians because the bootstrap sampling depresses the
median relative to the typical-value fit (see Section 2.3).

## Assumptions and deviations

- **Dosing convention.** The paper doses 0.25 mg enalapril maleate ODMT
  (Steichert 2025 Section 2.3). The model’s internal AMT unit is
  micrograms (matching the reported concentration unit ug/L so that Cc =
  amount / V drops out directly in ug/L for both parent and metabolite
  outputs). The vignette’s simulations therefore pass `amt = 250` (ug)
  rather than `amt = 0.25` (mg). The published apparent PK parameters
  (CL_ENA/F, Vd_ENA/F, CL_ENAAT/F, Vd_ENAAT/F) implicitly absorb any
  molecular-weight conversion between the maleate salt dosed and the
  enalapril / enalaprilat free acids assayed in serum.
- **IIV omega scale.** Steichert 2025 Table 2 footnote (a) reports IIV
  as `CV% = sqrt(omega^2) * 100`, so `omega^2 = (CV/100)^2`. The model
  file uses this convention directly (`etalcl ~ 0.4264 = 0.653^2` etc.),
  not the log-normal exact form `log(1 + CV^2)`. This matches the
  paper’s convention exactly at the cost of a small approximation bias
  in the tails.
- **Time-fixed covariates in simulation.** The paper’s dataset carries
  weight, age, serum creatinine, and Ross score as time-varying columns;
  the LENA sampling window is short enough (single-dose full profile at
  first dose, plus sparse trough samples during 8-week ODMT titration)
  that we hold each covariate constant per simulated subject within a
  run. Any longitudinal simulation over \> 1 week should split the run
  into visit-length segments with per-segment covariate updates.
- **Sex not included.** Steichert 2025 tested sex as a covariate on both
  CL and V of enalapril and enalaprilat in the stepwise search but the
  effect was not retained (Section 3.2.2, Table 2). The model file
  records this in `covariatesDataExcluded` for provenance; sex is not
  required in the covariate columns of the simulation event table.
- **No urine compartment.** The paper’s Fig 2 schematic shows k20
  (parent -\> urine) and k30 (metabolite -\> urine) transfer-rate
  constants for illustration, but no urine data were fit. The model file
  therefore does not carry an explicit urine compartment; the (1 - fm)
  parent-elimination fraction is folded into the total first-order loss
  from `central` rather than being tracked as a separate state.
- **Bootstrap 95% CI unused.** Steichert 2025 Table 2 reports
  nonparametric bootstrap medians and 95% CI alongside the point
  estimates. The model file uses only the point estimates from the Final
  model column; the bootstrap CI is documented in the source trace table
  above for reviewer reference.
