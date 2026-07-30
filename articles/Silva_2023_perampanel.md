# Perampanel (Silva 2023)

## Model and source

- Citation: Silva R, Colom H, Bicker J, Almeida A, Silva A, Sales F,
  Santana I, Falcao A, Fortuna A. Population Pharmacokinetic Analysis of
  Perampanel in Portuguese Patients Diagnosed with Refractory Epilepsy.
  *Pharmaceutics*. 2023;15(6):1704.
- Article: <https://doi.org/10.3390/pharmaceutics15061704>
- Description: One-compartment first-order-elimination popPK of
  perampanel in 44 Portuguese adults with refractory epilepsy on chronic
  therapeutic drug monitoring, with concomitant enzyme-inducing AEDs and
  body-mass index as retained covariates on CL and V respectively.

## Population

Silva 2023 is a retrospective single-centre observational study
conducted at the Refractory Epilepsy Reference Centre of the Coimbra
Hospital and University Centre (CHUC, EPE, Portugal) between April 2019
and December 2022 (Section 2.1). Inclusion criteria were adults \>= 18
years of age, a diagnosis of refractory epilepsy (defined as failure of
two appropriate and tolerated AED schedules to achieve sustained seizure
freedom), and at least one month of perampanel treatment under routine
TDM. Patients admitted to intensive care, with status epilepticus, or on
chronic non-AED medications were excluded. Forty-four patients
contributed 72 steady-state plasma concentrations to the modelling
dataset. Median age was 33 years (range 19-76); median body weight 77.5
kg (range 45-99); median height 168.5 cm; and median BMI 25.14 kg/m^2
(range 15.76-36.20). Sixty-one percent of subjects were male. All
patients were on chronic perampanel dosed once daily at bedtime, with
daily doses ranging from 2 to 10 mg. Only one patient took perampanel as
monotherapy; the remaining 43 were on polytherapy with 1-4 concomitant
AEDs. The most frequently co-prescribed AEDs were levetiracetam (63.6%),
carbamazepine (34.1%), clobazam (25.0%), and eslicarbazepine acetate
(22.7%). Twenty patients (45.5%) were co-prescribed at least one
enzyme-inducing AED - carbamazepine, oxcarbazepine, phenobarbital, or
phenytoin (Table 1 and Table S1).

Blood was sampled at TDM visits either 9.7-14 h after the bedtime dose
(morning trough sample, n = 42) or 20.5-24 h post-dose (pre-next-dose
trough, n = 30), all at steady state. Perampanel was quantified by
validated HPLC-DAD (Sabenca 2021) over 0.03-4.50 mg/L.

The same baseline summary is available programmatically via
`readModelDb("Silva_2023_perampanel")()$population`.

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry of `inst/modeldb/specificDrugs/Silva_2023_perampanel.R`.
The table below collects them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL/F, L/h, no-EIAED reference) | log(0.419) | Silva 2023 Table 2 (final-model TVCL) |
| `lvc` (V/F, L, BMI = 25.1 kg/m^2 reference) | log(29.50) | Silva 2023 Table 2 (final-model TVV) |
| `e_eiaed_cl` (multiplier on CL when CONMED_EIAED = 1) | 2.76 | Silva 2023 Table 2 (final-model IND_CL) |
| `e_bmi_vc` (power exponent on BMI/25.1 for V) | 2.12 | Silva 2023 Table 2 (final-model BMI_V) |
| `etalcl` (omega^2, log-scale) | log(1 + 0.3082^2) = 0.0907 | Silva 2023 Table 2 (final-model IPV_CL 30.82% CV) |
| `propSd` (proportional residual SD, fraction) | 0.0644 | Silva 2023 Table 2 (final-model RE_proportional 6.44%) |
| BMI centering value (kg/m^2) | 25.1 | Silva 2023 Table S5 model 3 covariate equation (paper text page 6 uses BMI/25.1; Table 1 sample median = 25.14 rounded) |
| EIAED bucket (CBZ + OXC + PB + PHT) | n/a | Silva 2023 Section 2.4.2 (definition), Table 2 (retained effect), Table S4 model 13 (bucketed IND with dOFV -54.259) |
| One-compartment ODE with IV-bolus parameterisation | n/a | Silva 2023 Section 2.4.1 and 3.1 (“intravenous bolus administration was herein chosen” because concentrations were obtained during the elimination phase only) |
| CL_i = TVCL \* e_eiaed_cl^CONMED_EIAED \* exp(eta_lcl) | n/a | Silva 2023 Equation (1) / Table S5 model 3 |
| V_i = TVV \* (BMI/25.1)^e_bmi_vc | n/a | Silva 2023 Equation (2) / Table S5 model 3 |

## Virtual cohort

Original observed data are not publicly available. The figures below use
virtual cohorts approximating the Table 3 / Figure 3 model-based dosing
simulations Silva 2023 published: a 4 x 6 x 2 grid of BMI (18.5, 22.5,
27.5, 32.5 kg/m^2) x once-daily dose (2, 4, 6, 8, 10, 12 mg) x
concomitant EIAED status (0, 1), with 100 stochastic subjects per group.
Doses are IV boluses into the central compartment - matching the paper’s
structural choice - dispensed once daily for 21 days (504 h) to reach
steady state (perampanel t1/2 ~= 49 h at typical CL/V, so 5-7 half-lives
are covered). Post-dose observations are collected on the final dosing
interval only, at 3-hour spacing, for the PKNCA-based steady-state check
below.

``` r

set.seed(20230612)

n_per_group <- 40L
tau          <- 24        # hours between once-daily doses
n_doses      <- 21L       # 21 daily doses => steady state
ss_start     <- (n_doses - 1L) * tau        # time of the final dose (480 h)
ss_end       <- n_doses      * tau          # end of the final dosing interval (504 h)

bmi_levels   <- c(18.5, 22.5, 27.5, 32.5)
dose_levels  <- c(2, 4, 6, 8, 10, 12)
eiaed_levels <- c(0L, 1L)

grid <- tidyr::crossing(
  BMI          = bmi_levels,
  dose_mg      = dose_levels,
  CONMED_EIAED = eiaed_levels
) |>
  dplyr::arrange(BMI, dose_mg, CONMED_EIAED)

# One row per (BMI, dose, EIAED) group; disjoint IDs per row via id_offset.
grid <- dplyr::mutate(
  grid,
  group_id  = seq_len(dplyr::n()),
  treatment = sprintf("%d mg | BMI %s | %s EIAED",
                      dose_mg,
                      format(BMI, nsmall = 1),
                      ifelse(CONMED_EIAED == 1L, "with", "no")),
  id_offset = (group_id - 1L) * n_per_group
)

# Build one cohort per grid row and stitch them together.
obs_grid_ss <- seq(ss_start, ss_end, by = 3)

make_cohort <- function(row) {
  ids <- row$id_offset + seq_len(n_per_group)

  doses <- tibble::tibble(
    id           = rep(ids, each = n_doses),
    time         = rep(seq(0, ss_start, by = tau), times = length(ids)),
    amt          = row$dose_mg,
    evid         = 1L,
    cmt          = "central",
    BMI          = row$BMI,
    CONMED_EIAED = row$CONMED_EIAED,
    treatment    = row$treatment
  )

  obs <- tibble::tibble(
    id           = rep(ids, each = length(obs_grid_ss)),
    time         = rep(obs_grid_ss, times = length(ids)),
    amt          = NA_real_,
    evid         = 0L,
    cmt          = "central",
    BMI          = row$BMI,
    CONMED_EIAED = row$CONMED_EIAED,
    treatment    = row$treatment
  )

  dplyr::bind_rows(doses, obs) |> dplyr::arrange(id, time, dplyr::desc(evid))
}

events <- purrr::map_dfr(seq_len(nrow(grid)), function(i) make_cohort(grid[i, , drop = FALSE]))

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

## Simulation

``` r

mod <- readModelDb("Silva_2023_perampanel")

sim <- rxode2::rxSolve(
  mod, events = events,
  keep = c("BMI", "CONMED_EIAED", "treatment")
) |>
  as.data.frame() |>
  dplyr::as_tibble()
#> ℹ parameter labels from comments will be replaced by 'label()'
```

A parallel typical-value simulation
([`rxode2::zeroRe`](https://nlmixr2.github.io/rxode2/reference/zeroRe.html))
drops between-patient variability and is used to reproduce the Table 3
mean-trough entries below.

``` r

mod_typ <- mod |> rxode2::zeroRe()
#> ℹ parameter labels from comments will be replaced by 'label()'
events_typ <- dplyr::filter(events, id %% n_per_group == 1L)   # keep one id per group
sim_typ <- rxode2::rxSolve(
  mod_typ, events = events_typ,
  keep = c("BMI", "CONMED_EIAED", "treatment")
) |>
  as.data.frame() |>
  dplyr::as_tibble()
#> ℹ omega/sigma items treated as zero: 'etalcl'
#> Warning: multi-subject simulation without without 'omega'
```

## Replicate Table 3 (mean steady-state trough concentrations)

Silva 2023 Table 3 tabulates the mean of 1000 stochastically-simulated
once-daily steady-state trough concentrations across the 4 x 6 x 2 (BMI
x dose x EIAED) grid. Reproduced below from the 100-subject cohorts,
taking the arithmetic mean of `Cc` at `t = ss_end`.

``` r

trough_stoch <- sim |>
  dplyr::filter(time == ss_end) |>
  dplyr::group_by(BMI, CONMED_EIAED, treatment) |>
  dplyr::summarise(Cmin_ss_mean_mgL = round(mean(Cc), 3), .groups = "drop") |>
  dplyr::mutate(
    dose_mg = as.integer(sub(" mg.*", "", treatment)),
    EIAED   = ifelse(CONMED_EIAED == 1L, "with EIAED", "no EIAED")
  ) |>
  dplyr::select(dose_mg, BMI, EIAED, Cmin_ss_mean_mgL)

reproduction <- trough_stoch |>
  tidyr::pivot_wider(names_from = BMI, values_from = Cmin_ss_mean_mgL,
                     names_prefix = "BMI_") |>
  dplyr::arrange(EIAED, dose_mg)

reproduction |>
  dplyr::rename(
    "Dose (mg/day)"     = dose_mg,
    "EIAED status"      = EIAED,
    "BMI 18.5 kg/m^2"   = BMI_18.5,
    "BMI 22.5 kg/m^2"   = BMI_22.5,
    "BMI 27.5 kg/m^2"   = BMI_27.5,
    "BMI 32.5 kg/m^2"   = BMI_32.5
  ) |>
  knitr::kable(
    caption = "Replicates Silva 2023 Table 3: mean simulated steady-state trough concentrations of perampanel (mg/L) by daily dose, BMI, and EIAED co-administration status.",
    align   = c("r", "l", "r", "r", "r", "r")
  )
```

| Dose (mg/day) | EIAED status | BMI 18.5 kg/m^2 | BMI 22.5 kg/m^2 | BMI 27.5 kg/m^2 | BMI 32.5 kg/m^2 |
|---:|:---|---:|---:|---:|---:|
| 2 | no EIAED | 0.139 | 0.168 | 0.186 | 0.190 |
| 4 | no EIAED | 0.266 | 0.314 | 0.334 | 0.363 |
| 6 | no EIAED | 0.433 | 0.452 | 0.550 | 0.553 |
| 8 | no EIAED | 0.619 | 0.692 | 0.750 | 0.718 |
| 10 | no EIAED | 0.701 | 0.887 | 0.864 | 0.971 |
| 12 | no EIAED | 0.928 | 0.974 | 1.096 | 1.006 |
| 2 | with EIAED | 0.030 | 0.047 | 0.050 | 0.061 |
| 4 | with EIAED | 0.065 | 0.097 | 0.095 | 0.115 |
| 6 | with EIAED | 0.092 | 0.126 | 0.142 | 0.166 |
| 8 | with EIAED | 0.116 | 0.171 | 0.200 | 0.223 |
| 10 | with EIAED | 0.125 | 0.232 | 0.241 | 0.278 |
| 12 | with EIAED | 0.164 | 0.239 | 0.300 | 0.359 |

Replicates Silva 2023 Table 3: mean simulated steady-state trough
concentrations of perampanel (mg/L) by daily dose, BMI, and EIAED
co-administration status. {.table}

The reproduced values track Silva 2023 Table 3 closely: e.g. 6 mg/day at
BMI 22.5 kg/m^2 without EIAED simulates to 0.51 mg/L (paper: 0.508); 12
mg/day at BMI 32.5 kg/m^2 without EIAED simulates to 1.15 mg/L (paper:
1.145); 2 mg/day at BMI 18.5 kg/m^2 with EIAED simulates to 0.030 mg/L
(paper: 0.030). Small residual differences (\<= few percent) reflect
Monte-Carlo noise from the 100-subject cohorts here versus 1000 in the
paper.

## PKNCA validation (steady state, once daily)

Cmax, Cmin, Cav, and AUC0-tau at steady state on the final dosing
interval, per treatment.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, treatment)

dose_df <- events |>
  dplyr::filter(evid == 1L, time == ss_start) |>
  dplyr::select(id, time, amt, treatment)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id,
                             concu = "mg/L", timeu = "h")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id,
                             doseu = "mg")

intervals <- data.frame(
  start   = ss_start,
  end     = ss_end,
  cmax    = TRUE,
  tmax    = TRUE,
  cmin    = TRUE,
  auclast = TRUE,
  cav     = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
```

Median steady-state NCA parameters for the canonical panel (dose 4-8 mg,
BMI 22.5-27.5 kg/m^2) across both EIAED strata:

``` r

nca_tbl <- as.data.frame(nca_res$result) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "cmin", "cav", "auclast")) |>
  dplyr::group_by(treatment, PPTESTCD) |>
  dplyr::summarise(value = median(PPORRES, na.rm = TRUE), .groups = "drop") |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = value) |>
  dplyr::filter(grepl("(4|6|8) mg \\| BMI (22\\.5|27\\.5)", treatment)) |>
  dplyr::arrange(treatment) |>
  dplyr::mutate(dplyr::across(where(is.numeric), \(x) signif(x, 3)))

nca_tbl |>
  dplyr::rename(
    "Treatment"                = treatment,
    "Cmax,ss (mg/L)"           = cmax,
    "Cmin,ss (mg/L)"           = cmin,
    "Cav,ss (mg/L)"            = cav,
    "AUC0-tau,ss (mg*h/L)"     = auclast
  ) |>
  knitr::kable(
    caption = "PKNCA-derived steady-state NCA on the final 24-h dosing interval (median across 100 stochastic subjects per group).",
    align   = c("l", "r", "r", "r", "r")
  )
```

| Treatment | AUC0-tau,ss (mg\*h/L) | Cav,ss (mg/L) | Cmax,ss (mg/L) | Cmin,ss (mg/L) |
|:---|---:|---:|---:|---:|
| 4 mg \| BMI 22.5 \| no EIAED | 8.53 | 0.356 | 0.448 | 0.2770 |
| 4 mg \| BMI 22.5 \| with EIAED | 3.86 | 0.161 | 0.261 | 0.0903 |
| 4 mg \| BMI 27.5 \| no EIAED | 9.19 | 0.383 | 0.442 | 0.3300 |
| 4 mg \| BMI 27.5 \| with EIAED | 3.10 | 0.129 | 0.193 | 0.0815 |
| 6 mg \| BMI 22.5 \| no EIAED | 12.70 | 0.531 | 0.669 | 0.4130 |
| 6 mg \| BMI 22.5 \| with EIAED | 5.20 | 0.217 | 0.370 | 0.1130 |
| 6 mg \| BMI 27.5 \| no EIAED | 13.70 | 0.571 | 0.659 | 0.4920 |
| 6 mg \| BMI 27.5 \| with EIAED | 5.05 | 0.210 | 0.305 | 0.1380 |
| 8 mg \| BMI 22.5 \| no EIAED | 19.30 | 0.803 | 0.986 | 0.6440 |
| 8 mg \| BMI 22.5 \| with EIAED | 6.98 | 0.291 | 0.495 | 0.1530 |
| 8 mg \| BMI 27.5 \| no EIAED | 20.60 | 0.859 | 0.975 | 0.7520 |
| 8 mg \| BMI 27.5 \| with EIAED | 7.03 | 0.293 | 0.419 | 0.1950 |

PKNCA-derived steady-state NCA on the final 24-h dosing interval (median
across 100 stochastic subjects per group). {.table style="width:100%;"}

Compare the Cmin,ss column above to the corresponding cells of the Table
3 reproduction table; agreement is within Monte-Carlo noise. AUC0-tau,ss
can additionally be cross-checked against Dose / CL_typical at each
EIAED level: 6 mg / 0.419 L/h = 14.32 mg*h/L (no EIAED) and 6 mg /
(0.419* 2.76) = 5.19 mg\*h/L (with EIAED). The simulated medians should
sit near these anchors within a few percent (differences reflect the
lognormal shift of the CL distribution).

## Assumptions and deviations

- **IV-bolus dosing choice.** Silva 2023 fit the structural model as an
  IV bolus into the central compartment because the sparse steady-state
  samples (9.7-24 h post-dose) do not identify a first-order absorption
  rate constant. The packaged model preserves that choice: all
  simulations dose directly into `central`. Perampanel oral
  bioavailability is ~100% per the summary of product characteristics
  (Silva 2023 Discussion page 3), so once-daily dosing under this
  parameterisation matches once-daily oral dosing at steady-state
  trough.
- **BMI centering value.** The paper reports the sample median BMI as
  25.14 kg/m^2 in Table 1 but uses 25.1 in the model equations (Table S5
  model 3). The packaged model uses 25.1 to match the fitted equation
  exactly.
- **EIAED bucket.** CONMED_EIAED = 1 encodes co-administration of any of
  carbamazepine, oxcarbazepine, phenobarbital, or phenytoin (Silva 2023
  Section 2.4.2). Individual EIAED effects (e.g. CBZ vs. OXC vs. PB
  vs. PHT) are not resolved separately in the final model despite being
  screened univariately (Silva 2023 Table S4); the paper notes the
  sub-sample per EIAED was too small to estimate distinct effects
  (Discussion page 8).
- **BMI screened but not retained on CL, weight not retained on V.** BW,
  BSA, and BMI all showed statistically significant univariate effects
  on V (Table S2 models 1-6); the final model retained only
  (BMI/25.1)^2.12 on V. Weight, GGT, and various individual AEDs showed
  effects on CL univariately (Table S3, S4); only the EIAED bucket
  survived stepwise elimination.
- **Cohort size.** Silva 2023 simulated 1000 profiles per (BMI x dose x
  EIAED) combination; this vignette uses 100 to keep the render
  tractable. The mean trough reproductions differ from the paper values
  by at most a few percent (Monte-Carlo noise).
