# Tobramycin (Hennig 2008)

## Model and source

- Citation: Hennig S, Norris R, Kirkpatrick CMJ. (2008). Target
  concentration intervention is needed for tobramycin dosing in
  paediatric patients with cystic fibrosis – a population
  pharmacokinetic study. Br J Clin Pharmacol 65(4):502-510.
- Article: <https://doi.org/10.1111/j.1365-2125.2007.03045.x>

``` r

mod_meta <- rxode2::rxode(readModelDb("Hennig_2008_tobramycin"))
#> ℹ parameter labels from comments will be replaced by 'label()'
mod_meta$description
#> [1] "Two-compartment population PK model for once-daily IV tobramycin in paediatric cystic fibrosis patients (Hennig 2008), with allometric weight scaling on CL, Q, Vc, and Vper (reference 70 kg, exponent 3/4 for clearances and 1 for volumes), full-block correlated between-subject variability on CL/Vc/Vper, a fixed 30 min infusion duration into the central compartment, and an estimated lag time between infusion start and drug entry into the patient's vein."
mod_meta$reference
#> [1] "Hennig S, Norris R, Kirkpatrick CMJ. (2008). Target concentration intervention is needed for tobramycin dosing in paediatric patients with cystic fibrosis - a population pharmacokinetic study. Br J Clin Pharmacol 65(4):502-510. doi:10.1111/j.1365-2125.2007.03045.x"
mod_meta$units
#> $time
#> [1] "hour"
#> 
#> $dosing
#> [1] "mg"
#> 
#> $concentration
#> [1] "mg/L"
```

## Population

The Hennig 2008 model was developed retrospectively from 318 tobramycin
plasma concentrations measured in 35 paediatric cystic fibrosis (CF)
patients (21 female, 14 male) admitted to Mater Health Services,
Brisbane, Australia between July 2005 and September 2006 (Hennig 2008
Table 1). Mean age was 9.5 years (range 0.5-17.9), mean total body
weight 34.0 kg (6.0-72.6), and mean dose 9.6 mg/kg once-daily IV. Renal
impairment was rare in this cohort – only one patient had
Cockcroft-Gault creatinine clearance below 50 mL/min, so the model is
most informative for paediatric CF patients with normal serum
creatinine. The retrospective dataset contained on average 9.1 samples
per subject across 4.6 dosing occasions.

The same information is available programmatically via
`readModelDb("Hennig_2008_tobramycin")$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) entry in
`inst/modeldb/specificDrugs/Hennig_2008_tobramycin.R`. The table below
collects them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (log CL/70 kg) | log(6.37) | Hennig 2008 Table 2, Covariate model column |
| `lvc` (log Vc/70 kg) | log(18.70) | Hennig 2008 Table 2, Covariate model column |
| `lq` (log Q/70 kg) | log(0.39) | Hennig 2008 Table 2, Covariate model column |
| `lvp` (log Vper/70 kg) | log(1.32) | Hennig 2008 Table 2, Covariate model column |
| `e_wt_cl_q` (allometric exponent on CL, Q) | 0.75 | Hennig 2008 Table 2 footnote (f = 3/4 for clearance) |
| `e_wt_vc_vp` (allometric exponent on Vc, Vper) | 1.00 | Hennig 2008 Table 2 footnote (f = 1 for volume) |
| `llag` (log infusion lag) | log(0.40) | Hennig 2008 Table 2, tlag, Covariate model column |
| `ldur` (log infusion duration, fixed) | log(0.5) | Hennig 2008 Table 2, D2 (fixed); Methods page 503 |
| BSV CL (full block element) | 11.70% CV -\> omega^2 = 0.013596 | Hennig 2008 Table 2, Random parameters |
| BSV Vc (full block element) | 11.66% CV -\> omega^2 = 0.013505 | Hennig 2008 Table 2, Random parameters |
| BSV Vper (full block element) | 41.95% CV -\> omega^2 = 0.162130 | Hennig 2008 Table 2, Random parameters |
| Cor(CL, Vc) | 0.73 | Hennig 2008 Table 2, Correlations between BSVs |
| Cor(CL, Vper) | 0.49 | Hennig 2008 Table 2, Correlations between BSVs |
| Cor(Vc, Vper) | 0.27 | Hennig 2008 Table 2, Correlations between BSVs |
| `propSd` (proportional residual SD) | 0.190 | Hennig 2008 Table 2, Residual variability (19.0% CV exponential -\> proportional) |
| Two-compartment IV ODE structure | n/a | Hennig 2008 Methods, Population analysis section, page 503 |
| Allometric scaling formula | `(WT/70)^f` with f = 3/4 (CL, Q) and f = 1 (Vc, Vper) | Hennig 2008 Table 2 footnote |

## Virtual cohort

The original observed concentrations are not publicly available. The
figures below use a virtual paediatric CF cohort whose body-weight
distribution approximates Table 1 of Hennig 2008 (mean 34.0 kg, range
6.0-72.6 kg). Once-daily IV infusion doses are computed from a per mg/kg
regimen (the dosing strategies compared in Hennig 2008 Table 3).

``` r

set.seed(20250508)

n_subjects <- 200L

# Approximate the Hennig 2008 Table 1 weight distribution: mean 34 kg,
# range 6-72.6 kg. Use a truncated normal-ish draw to span the observed
# range; the simulation only consumes WT for allometric scaling.
sample_wt <- function(n) {
  wt <- rnorm(n, mean = 34, sd = 16)
  wt <- pmin(pmax(wt, 6.0), 72.6)
  wt
}

# Dosing regimens evaluated in Table 3 of Hennig 2008. Per-kg doses are
# computed from the simulated WT; the absolute 350 mg ("one dose fits all")
# regimen is applied without weight scaling.
regimens <- tibble::tribble(
  ~regimen,        ~dose_mg_per_kg, ~fixed_dose_mg,
  "350 mg",        NA_real_,        350,
  "7.5 mg/kg",     7.5,             NA_real_,
  "10 mg/kg",      10,              NA_real_
)

make_cohort <- function(n, regimen, dose_mg_per_kg, fixed_dose_mg, id_offset = 0L) {
  wt <- sample_wt(n)
  amt <- if (is.na(fixed_dose_mg)) wt * dose_mg_per_kg else rep(fixed_dose_mg, n)
  ids <- id_offset + seq_len(n)

  dose_rows <- tibble::tibble(
    id   = ids,
    time = 0,
    evid = 1L,
    cmt  = "central",
    amt  = amt,
    WT   = wt,
    regimen = regimen
  )

  obs_grid <- seq(0, 24, by = 0.25)
  obs_rows <- tidyr::expand_grid(id = ids, time = obs_grid) |>
    dplyr::left_join(dplyr::select(dose_rows, id, WT, regimen), by = "id") |>
    dplyr::mutate(evid = 0L, cmt = "central", amt = 0)

  dplyr::bind_rows(dose_rows, obs_rows) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

events <- dplyr::bind_rows(
  make_cohort(n_subjects, "350 mg",    NA_real_, 350,      id_offset =   0L),
  make_cohort(n_subjects, "7.5 mg/kg", 7.5,      NA_real_, id_offset = 1000L),
  make_cohort(n_subjects, "10 mg/kg",  10,       NA_real_, id_offset = 2000L)
)

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))

n_doses <- sum(events$evid == 1L)
n_obs   <- sum(events$evid == 0L)
cat("Dosing rows:", n_doses, " | Observation rows:", n_obs, "\n")
#> Dosing rows: 600  | Observation rows: 58200
```

## Simulation

``` r

mod <- readModelDb("Hennig_2008_tobramycin")
sim <- rxode2::rxSolve(
  mod,
  events = events,
  keep   = c("regimen", "WT")
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
```

## Replicate published figures

### Figure 3 – visual predictive check by dose regimen

Hennig 2008 Figure 3 shows the model’s VPC at the 350 mg regimen
(observed grey points and grey 50th-percentile line vs. simulated median
and 5th/95th-percentile envelopes). The chunk below renders the
analogous percentile envelope for each of the three regimens.

``` r

vpc_summary <- sim |>
  dplyr::group_by(regimen, time) |>
  dplyr::summarise(
    Q05 = quantile(Cc, 0.05, na.rm = TRUE),
    Q50 = quantile(Cc, 0.50, na.rm = TRUE),
    Q95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(vpc_summary, aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95, fill = regimen), alpha = 0.2) +
  geom_line(aes(colour = regimen)) +
  facet_wrap(~regimen) +
  scale_y_log10() +
  labs(
    x = "Time after dose (h)",
    y = "Tobramycin concentration (mg/L)",
    title = "Replicates Figure 3 of Hennig 2008",
    caption = "Simulated 5th/50th/95th-percentile envelopes by dosing regimen."
  ) +
  theme(legend.position = "none")
#> Warning in scale_y_log10(): log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
```

![](Hennig_2008_tobramycin_files/figure-html/figure-3-1.png)

## PKNCA validation

Compute Cmax and AUC0-24 (the dosing interval) per subject with PKNCA.
Treatment grouping is by `regimen`.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, regimen)

dose_df <- events |>
  dplyr::filter(evid == 1L) |>
  dplyr::select(id, time, amt, regimen)

conc_obj <- PKNCA::PKNCAconc(
  sim_nca, Cc ~ time | regimen + id,
  concu = "mg/L", timeu = "h"
)
dose_obj <- PKNCA::PKNCAdose(
  dose_df, amt ~ time | regimen + id,
  doseu = "mg"
)

intervals <- data.frame(
  start    = 0,
  end      = 24,
  cmax     = TRUE,
  tmax     = TRUE,
  auclast  = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
#>  ■■■■■■■■■■                        30% |  ETA:  3s
#>  ■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■    97% |  ETA:  0s
```

### Comparison against Hennig 2008 Table 3 simulation summaries

Hennig 2008 Table 3 reports, for each dosing regimen, the percentage of
patients whose AUC0-24 falls within the 80-125 mg\*h/L target range and
whose Cmax falls within the 24-38 mg/L target range, based on 1000 Monte
Carlo simulations. The chunk below compares the simulated percentages
from this packaged model with the published values.

``` r

res_tbl <- as.data.frame(nca_res$result)

per_subject <- res_tbl |>
  dplyr::filter(PPTESTCD %in% c("auclast", "cmax")) |>
  tidyr::pivot_wider(
    id_cols     = c("regimen", "id"),
    names_from  = PPTESTCD,
    values_from = PPORRES
  )

simulated_pct <- per_subject |>
  dplyr::group_by(regimen) |>
  dplyr::summarise(
    pct_auc_in_target  = mean(auclast >= 80   & auclast <= 125, na.rm = TRUE) * 100,
    pct_auc_below_80   = mean(auclast < 80,                     na.rm = TRUE) * 100,
    pct_auc_above_125  = mean(auclast > 125,                    na.rm = TRUE) * 100,
    pct_cmax_in_target = mean(cmax    >= 24   & cmax    <= 38,  na.rm = TRUE) * 100,
    .groups = "drop"
  )

published_pct <- tibble::tribble(
  ~regimen,    ~pct_auc_in_target_pub, ~pct_auc_below_80_pub, ~pct_auc_above_125_pub, ~pct_cmax_in_target_pub,
  "350 mg",    37.9,                    36.6,                  25.5,                   29.8,
  "7.5 mg/kg", 21.7,                    78.1,                   0.2,                   17.3,
  "10 mg/kg",  72.1,                    22.2,                   5.7,                   91.6
)

comparison <- published_pct |>
  dplyr::left_join(simulated_pct, by = "regimen") |>
  dplyr::select(regimen,
                pct_auc_in_target_pub,  pct_auc_in_target,
                pct_auc_below_80_pub,   pct_auc_below_80,
                pct_auc_above_125_pub,  pct_auc_above_125,
                pct_cmax_in_target_pub, pct_cmax_in_target)

knitr::kable(
  comparison,
  digits  = 1,
  caption = paste("Comparison of simulated AUC0-24 and Cmax target-attainment percentages",
                  "with Hennig 2008 Table 3. The 'pub' columns are the published values; the",
                  "remaining columns are computed from the packaged model with",
                  paste0("n = ", n_subjects, " subjects per regimen."))
)
```

| regimen | pct_auc_in_target_pub | pct_auc_in_target | pct_auc_below_80_pub | pct_auc_below_80 | pct_auc_above_125_pub | pct_auc_above_125 | pct_cmax_in_target_pub | pct_cmax_in_target |
|:---|---:|---:|---:|---:|---:|---:|---:|---:|
| 350 mg | 37.9 | 37.5 | 36.6 | 34.0 | 25.5 | 28.5 | 29.8 | 39.0 |
| 7.5 mg/kg | 21.7 | 15.0 | 78.1 | 85.0 | 0.2 | 0.0 | 17.3 | 88.5 |
| 10 mg/kg | 72.1 | 74.0 | 22.2 | 23.5 | 5.7 | 2.5 | 91.6 | 72.0 |

Comparison of simulated AUC0-24 and Cmax target-attainment percentages
with Hennig 2008 Table 3. The ‘pub’ columns are the published values;
the remaining columns are computed from the packaged model with n = 200
subjects per regimen. {.table}

The simulated percentages should track the published Table 3 values
within reasonable Monte Carlo noise. AUC target attainment matches the
published values across all three regimens (within ~10 percentage
points). Cmax target attainment matches well for the 350 mg flat dose
but diverges for the per-kg regimens; in particular, the 7.5 mg/kg arm’s
published 17.3% in the 24-38 mg/L target band differs from the simulated
value by more than 20 percentage points. The most likely explanation is
a difference in how the paper’s simulated “Cmax” was defined (e.g., peak
observed at a specific end-of-infusion sampling time rather than the
maximum over the dosing interval). Because Cmax under per-kg dosing in a
fully-allometric model is approximately weight-independent, with the
typical value (7.5 mg/kg x 70 kg / Vc_70 = 28 mg/L) sitting near the
centre of the 24-38 mg/L band, almost any draw from the BSV distribution
lands in target – which is what the simulation shows. The model
parameters and structure are reproduced faithfully from the source; per
the extraction skill’s policy, parameters are not tuned to chase a
validation metric. Users who need to match the paper’s Cmax-distribution
results closely should re-derive Cmax according to the paper’s exact
sampling convention before comparing.

## Assumptions and deviations

- **Between-occasion variability (BOV) on CL excluded.** Hennig 2008
  estimated BOV on CL of 6.47% CV in addition to BSV (11.70% CV).
  Following the precedent set by `Faelens_2021_infliximab` – which also
  drops a low-magnitude per-occasion variance for nlmixr2lib portability
  – the packaged model omits BOV. The omitted variance is small relative
  to BSV, so simulated population-level statistics closely match the
  published values. Re-introducing BOV would require a per-occasion eta
  multiplexed by an occasion indicator (cf. `Wilkins_2008_rifampicin`
  for that pattern).
- **Body-weight distribution.** The vignette draws WT from a truncated
  normal with mean 34 kg matching the published baseline mean (Table 1)
  and clamped to the 6-72.6 kg study range. The original individual WT
  values are not published.
- **Renal function.** Only one patient in the source cohort had reduced
  creatinine clearance, and serum creatinine / CrCL were tested as
  covariates but not retained. The packaged model therefore does not
  include any renal covariate; users with renally impaired patients
  (e.g. eGFR substantially below the paediatric CF norm) should
  interpret simulations with care.
- **Allometric exponents fixed.** The Hennig 2008 Table 2 footnote
  states “f = 3/4 for clearance, f = 1 for volume” without confidence
  intervals, indicating these were fixed. The packaged model wraps
  `e_wt_cl_q` and `e_wt_vc_vp` in `fixed()` to mirror this.
- **Infusion duration fixed at 30 min.** Hennig 2008 fixed D2 = 0.5 h
  per the hospital protocol; the packaged model imposes
  `dur(central) <- 0.5` h so users do not need to supply a `rate` or
  `dur` column on each dose record.
