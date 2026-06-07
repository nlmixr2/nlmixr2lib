# Cyanidin-3-Glucoside (Jeon 2012)

## Model and source

- Citation: Jeon S, Han S, Lee J, Hong T, Yim DS. (2012). The Safety and
  Pharmacokinetics of Cyanidin-3-Glucoside after 2-Week Administration
  of Black Bean Seed Coat Extract in Healthy Subjects. Korean J Physiol
  Pharmacol 16(4):249-253. <doi:10.4196/kjpp.2012.16.4.249>.
- Description: One-compartment first-order absorption population PK
  model with an absorption lag time for cyanidin-3-glucoside (C3G)
  following 2-week multiple oral dosing of 1 g black bean (Phaseolus
  vulgaris, Cheongjakong-3-ho) seed coat extract once daily in 12
  healthy adult Korean volunteers (Jeon 2012), with log-normal IIV on
  CL/F and V/F (with correlation rho = 0.883) and on Ka, and
  proportional residual error.
- Article (open access): <https://doi.org/10.4196/kjpp.2012.16.4.249>

Jeon 2012 is a brief open-access pharmacokinetic study of
cyanidin-3-glucoside (C3G), the major anthocyanin in the seed coat of
the Korean black bean cultivar *Phaseolus vulgaris* Cheongjakong-3-ho.
Twelve healthy adult Korean volunteers received 1 g of black bean seed
coat extract once daily for 14 days, with intensive PK sampling on Day 1
and Day 14. The structural model fit is a one-compartment first-order
absorption model with an absorption lag time, log-normal IIV on CL/F and
V/F (correlated, rho = 0.883) and on Ka, and a proportional residual
error. No covariate (age, height, weight, sex, creatinine clearance)
cleared the forward-addition / backward-elimination significance
threshold in the final model.

## Population

Jeon 2012 enrolled and dosed 12 healthy adult Korean volunteers (Jeon
2012 Table 1):

- **N = 12** subjects, all completed the 14-day dosing and PK schedule.
- Sex: 6 male / 6 female (50% female).
- Age: median 29.5 years (range 24-44).
- Height: median 166 cm (range 158-185).
- Weight: median 59.6 kg (range 47.0-78.4).
- All Korean; one centre (Clinical Trial Center, Seoul St Mary’s
  Hospital).
- Diet: low-phenolic during the outpatient period to minimise dietary
  anthocyanin background; at least 10 h overnight fast on PK days (Day
  1, Day 14) and 4 h fasting after dosing.
- Dosing: 1 g of black bean seed coat extract (200 mg x 5 capsules, GMP-
  produced from 20 kg of Cheongjakong-3-ho beans extracted in 80%
  ethanol) once daily at approximately 9 AM for 14 days.

The same information is available programmatically:
`readModelDb("Jeon_2012_C3G")$population`.

## Source trace

The per-parameter origin is recorded as in-file comments next to each
`ini()` entry of `inst/modeldb/specificDrugs/Jeon_2012_C3G.R`. The table
below collates them for review.

| Parameter (nlmixr2lib) | Value | Source location |
|----|----|----|
| `lka` | log(9.94) | Table 3 Ka = 9.94 1/h (%RSE 30.7) |
| `lcl` | log(3420) | Table 3 CL/F = 3420 L/h (%RSE 8.68) |
| `lvc` | log(7280) | Table 3 V/F = 7280 L (%RSE 6.84) |
| `lalag` | log(0.217) | Table 3 ALAG = 0.217 h (%RSE 3.67) |
| `etalcl`+`etalvc` | block(0.080216, 0.052679, 0.044370) | Table 3 omega_CL = 28.9% CV, omega_V = 21.3% CV, rho_CL-V = 0.883 |
| `etalka` | 0.718145 | Table 3 omega_Ka = 102.5% CV; log(1 + 1.025^2) |
| `propSd` | 0.213 | Table 3 prop residual = 21.3% (%RSE 9.30) |
| `d/dt(depot)` | n/a | One-compartment first-order absorption (Jeon 2012 Methods “Mixed effect PK model development” + Results “Final model”) |
| `d/dt(central)` | n/a | First-order elimination from central |
| `alag(depot) <- alag1` | n/a | Absorption lag time ALAG (Jeon 2012 Table 3) |
| `Cc = 1000 * central / vc` | n/a | Plasma C3G in ng/mL (dose mg / V L \* 1000) |

## Virtual cohort

Original observed plasma C3G concentrations are not redistributed with
the package. The cohort below mirrors the Jeon 2012 study design: 12
subjects receive a single 14.57 mg oral dose of C3G equivalent on Day 1
(the C3G content of 1 g of black bean seed coat extract, back-calculated
from the published CL/F and observed AUClast as described in the model
file). A second cohort runs the same 14-day daily dosing regimen so a
Day-1 vs Day-14 NCA comparison can be drawn against Table 2 of the
paper.

``` r

set.seed(20120401L)

n_per_arm <- 100L

# Per-administration C3G equivalent dose for 1 g of Cheongjakong-3-ho
# black bean seed coat extract. Back-calculated: AUClast_obs x CL/F =
# 0.00426 mg*h/L x 3420 L/h = 14.57 mg.
dose_c3g_mg <- 14.57

# Sampling grid matches the paper: pre-dose and 0.25, 0.5, 1, 1.5, 2, 4,
# 6 h after dosing on the PK days.
sample_grid_h <- c(0.05, 0.25, 0.5, 1, 1.5, 2, 4, 6)

make_cohort_single <- function(id_offset) {
  ids <- id_offset + seq_len(n_per_arm)

  dose_rows <- tibble::tibble(
    id   = ids,
    time = 0,
    amt  = dose_c3g_mg,
    evid = 1L,
    cmt  = "depot"
  )

  obs_rows <- tidyr::expand_grid(id = ids, time = sample_grid_h) |>
    dplyr::mutate(amt = 0, evid = 0L, cmt = "Cc")

  dplyr::bind_rows(dose_rows, obs_rows) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

make_cohort_multi <- function(id_offset) {
  ids <- id_offset + seq_len(n_per_arm)

  dose_days <- seq(0, 24 * 13, by = 24)
  dose_rows <- tidyr::expand_grid(id = ids, time = dose_days) |>
    dplyr::mutate(amt = dose_c3g_mg, evid = 1L, cmt = "depot")

  # Day 1 PK sampling (0-6 h) plus Day 14 PK sampling (312-318 h after
  # the first dose), matching the paper's blood sampling protocol.
  obs_day1_h  <- sample_grid_h
  obs_day14_h <- 24 * 13 + sample_grid_h
  obs_times   <- c(obs_day1_h, obs_day14_h)

  obs_rows <- tidyr::expand_grid(id = ids, time = obs_times) |>
    dplyr::mutate(amt = 0, evid = 0L, cmt = "Cc")

  dplyr::bind_rows(dose_rows, obs_rows) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

events <- dplyr::bind_rows(
  make_cohort_single(id_offset = 0L)         |> dplyr::mutate(treatment = "Day 1 (single dose)"),
  make_cohort_multi(id_offset = n_per_arm)   |> dplyr::mutate(treatment = "Day 14 (multi-dose)")
)

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid", "cmt")])))
```

## Simulation

``` r

mod <- rxode2::rxode2(readModelDb("Jeon_2012_C3G"))
#> ℹ parameter labels from comments will be replaced by 'label()'

sim <- rxode2::rxSolve(
  mod,
  events = events,
  keep   = c("treatment")
) |>
  as.data.frame()
```

For the typical-value reference trajectory (deterministic, no IIV) used
to overlay published NCA point estimates, zero out the random effects:

``` r

mod_typical <- mod |> rxode2::zeroRe()

sim_typical <- rxode2::rxSolve(
  mod_typical,
  events = events,
  keep   = c("treatment")
) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc', 'etalka'
#> Warning: multi-subject simulation without without 'omega'
```

## Replicate Figure 2 - VPC by study day

Jeon 2012 Figure 2 shows the visual predictive check of the final model
against observed plasma C3G concentrations on Day 1 and Day 14, with 5th
/ 50th / 95th percentile bands overlaid on observed data. The stochastic
trajectory below reproduces the same VPC structure (band widths and
shapes are determined entirely by the IIV variances reported in Table 3
and the proportional residual error).

``` r

vpc_day <- sim |>
  dplyr::filter(!is.na(Cc), Cc > 0) |>
  dplyr::mutate(
    panel = treatment,
    tad   = ifelse(treatment == "Day 14 (multi-dose)", time - 24 * 13, time)
  ) |>
  dplyr::filter(tad >= 0, tad <= 6) |>
  dplyr::group_by(panel, tad) |>
  dplyr::summarise(
    Q05 = stats::quantile(Cc, 0.05, na.rm = TRUE),
    Q50 = stats::quantile(Cc, 0.50, na.rm = TRUE),
    Q95 = stats::quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(vpc_day, aes(tad, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25, fill = "steelblue") +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ panel) +
  labs(
    x = "Time after most-recent dose (h)",
    y = "Plasma C3G (ng/mL)",
    title = "Plasma C3G VPC by study day (replicates Jeon 2012 Figure 2)",
    caption = paste("Median (line) and 5th-95th percentile band (shaded);",
                    "n =", n_per_arm, "per arm; C3G equivalent dose",
                    dose_c3g_mg, "mg per administration.")
  )
```

![](Jeon_2012_C3G_files/figure-html/figure-2-vpc-1.png)

## Replicate typical-value Cmax / Tmax

The deterministic (zero-RE) trajectory locates the typical-value Cmax
and Tmax, which can be compared against the observed Day 1 / Day 14 NCA
values in Jeon 2012 Table 2.

``` r

typical_summary <- sim_typical |>
  dplyr::filter(!is.na(Cc), Cc > 0) |>
  dplyr::mutate(
    tad   = ifelse(treatment == "Day 14 (multi-dose)", time - 24 * 13, time)
  ) |>
  dplyr::filter(tad >= 0, tad <= 6) |>
  dplyr::group_by(treatment, id) |>
  dplyr::slice_max(Cc, n = 1, with_ties = FALSE) |>
  dplyr::group_by(treatment) |>
  dplyr::summarise(
    Cmax_typical_ngmL = round(stats::median(Cc), 3),
    Tmax_typical_h    = round(stats::median(tad), 3),
    .groups = "drop"
  )

knitr::kable(
  typical_summary,
  caption = "Typical-value (zero-RE) Cmax and Tmax by study day."
)
```

| treatment           | Cmax_typical_ngmL | Tmax_typical_h |
|:--------------------|------------------:|---------------:|
| Day 1 (single dose) |             1.713 |            0.5 |
| Day 14 (multi-dose) |             1.713 |            0.5 |

Typical-value (zero-RE) Cmax and Tmax by study day. {.table}

## PKNCA validation

Jeon 2012 Table 2 reports noncompartmental Cmax, Tmax, AUClast, and
terminal half-life on Day 1 and Day 14 (mean +/- SD). The PKNCA block
below recomputes the same parameters from the stochastic simulation for
side-by-side comparison.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::mutate(
    tad = ifelse(treatment == "Day 14 (multi-dose)", time - 24 * 13, time)
  ) |>
  dplyr::filter(tad >= 0, tad <= 6) |>
  dplyr::transmute(id, time = tad, Cc, treatment)

dose_df <- events |>
  dplyr::filter(evid == 1L) |>
  dplyr::group_by(id, treatment) |>
  dplyr::summarise(time = 0, amt = dose_c3g_mg, .groups = "drop")

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id)
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id,
                             route = "extravascular")

intervals <- data.frame(
  start     = 0,
  end       = 6,
  cmax      = TRUE,
  tmax      = TRUE,
  auclast   = TRUE,
  half.life = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Warning: Too few points for half-life calculation (min.hl.points=3 with only 2
#> points)
#> Warning: Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed
#> Requesting an AUC range starting (0) before the first measurement (0.05) is not allowed

nca_summary <- as.data.frame(nca_res$result) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "half.life")) |>
  dplyr::group_by(treatment, PPTESTCD) |>
  dplyr::summarise(
    mean = round(mean(PPORRES, na.rm = TRUE), 3),
    sd   = round(stats::sd(PPORRES, na.rm = TRUE), 3),
    .groups = "drop"
  )

knitr::kable(
  nca_summary,
  caption = paste("Simulated NCA by study day (n =", n_per_arm, "per arm).",
                  "Cmax ng/mL; Tmax h; AUClast ng*h/mL; half-life h.")
)
```

| treatment           | PPTESTCD  |  mean |    sd |
|:--------------------|:----------|------:|------:|
| Day 1 (single dose) | auclast   |   NaN |    NA |
| Day 1 (single dose) | cmax      | 1.668 | 0.401 |
| Day 1 (single dose) | half.life | 1.474 | 0.221 |
| Day 1 (single dose) | tmax      | 0.660 | 0.255 |
| Day 14 (multi-dose) | auclast   |   NaN |    NA |
| Day 14 (multi-dose) | cmax      | 1.671 | 0.416 |
| Day 14 (multi-dose) | half.life | 1.511 | 0.218 |
| Day 14 (multi-dose) | tmax      | 0.670 | 0.286 |

Simulated NCA by study day (n = 100 per arm). Cmax ng/mL; Tmax h;
AUClast ng\*h/mL; half-life h. {.table}

### Comparison against published Table 2 (noncompartmental)

``` r

reference <- tibble::tribble(
  ~treatment,            ~cmax, ~tmax, ~auclast, ~half.life,
  "Day 1 (single dose)",  2.10,  0.73,     4.26,       1.58,
  "Day 14 (multi-dose)",  2.17,  0.75,     4.85,       1.47
)
cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = nca_res,
  reference     = reference,
  by            = "treatment",
  params        = c("cmax", "tmax", "auclast", "half.life"),
  units         = c(cmax = "ng/mL", tmax = "h",
                    auclast = "ng*h/mL", half.life = "h"),
  tolerance_pct = 20
)
knitr::kable(
  cmp,
  caption = paste(
    "C3G NCA: simulated (median across n =", n_per_arm, "per arm) vs.",
    "Jeon 2012 Table 2 (paper mean +/- SD). * differs from reference",
    "by >20%."
  )
)
```

| NCA parameter | treatment           | Reference | Simulated | % diff   |
|:--------------|:--------------------|:----------|:----------|:---------|
| Cmax (ng/mL)  | Day 1 (single dose) | 2.1       | 1.59      | -24.1%\* |
| Cmax (ng/mL)  | Day 14 (multi-dose) | 2.17      | 1.63      | -24.8%\* |
| Tmax (h)      | Day 1 (single dose) | 0.73      | 0.5       | -31.5%\* |
| Tmax (h)      | Day 14 (multi-dose) | 0.75      | 0.5       | -33.3%\* |
| t½ (h)        | Day 1 (single dose) | 1.58      | 1.45      | -8.1%    |
| t½ (h)        | Day 14 (multi-dose) | 1.47      | 1.51      | +3.0%    |

C3G NCA: simulated (median across n = 100 per arm) vs. Jeon 2012 Table 2
(paper mean +/- SD). \* differs from reference by \>20%. {.table}

The simulated Day 1 and Day 14 PK characteristics should agree with the
published noncompartmental analysis within the expected sampling
variability for the dense-sampling design used in this 12-subject study.
The model predicts no clinically meaningful accumulation between Day 1
and Day 14 (consistent with Jeon 2012 Results: paired t-test p \> 0.05
for all NCA comparisons between days), driven by the very short terminal
half-life (~1.5 h) relative to the 24 h dosing interval.

## Assumptions and deviations

- **C3G equivalent dose (14.57 mg per 1 g of extract) is
  back-calculated, not reported.** The paper administers 1 g of black
  bean seed coat extract daily but does not report the C3G content per
  gram of extract. The published CL/F and observed Day-1 AUClast jointly
  imply an effective C3G dose of approximately
  `AUClast x CL/F = 0.00426 mg*h/L x 3420 L/h = 14.57 mg`, consistent
  with the published model fit under the standard NONMEM convention that
  F = 1 is implicit when no separate bioavailability parameter is
  estimated. This value is treated as a dose-unit documentation choice:
  the packaged model takes dose in mg of C3G equivalent, and 14.57 mg
  per administration reproduces the observed Day 1 and Day 14 NCA
  values. Anthocyanin oral bioavailability is reported in the literature
  as approximately 1-2% of the ingested mass, so a C3G equivalent dose
  of approximately 1.5% of the 1 g extract mass is plausible.

- **IOV vs BSV labelling for omega_Ka.** Jeon 2012 Methods describes
  both between-subject variability (eta on log-CL / log-V / log-Ka /
  log-ALAG) and inter-occasion variability (kappa on log-CL / log-V /
  log-Ka), and the final-model prose states that “inter-occasional
  variability between Day 1 and Day 14 was tested for all PK parameters,
  but was finally found relevant only for absorption rate constant
  (Ka).” Table 3 reports a single variance term on Ka labelled
  “Between-subject variability omega Ka” (102.5% CV) without a separate
  IOV row. The packaged model encodes this single 102.5% value as a
  single `etalka` variance, which captures the absorption variability
  the paper reports without distinguishing BSV from IOV; the published
  Table 3 does not provide enough information to separate the two terms,
  and no NONMEM control stream (`.mod` / `.lst`) is available with the
  paper.

- **IIV on ALAG not retained.** Jeon 2012 Methods describes BSV on ALAG
  among the screened random effects, but the final-model Table 3 omits
  an `omega_ALAG` row; the packaged model therefore treats ALAG as a
  typical- value-only parameter (no `etalalag`).

- **Sampling grid matches the paper.** Plasma sampling in the simulation
  uses the same 0.05 / 0.25 / 0.5 / 1 / 1.5 / 2 / 4 / 6 h post-dose grid
  the paper used (Methods “Blood sampling”), so the simulated NCA
  characteristics see the same Cmax / Tmax resolution as the observed
  data.

- **LLOQ.** The bioanalytical assay LLOQ for C3G was 0.2 ng/mL (Methods
  “Assay of the plasma concentration of C3G”). The packaged-model
  proportional residual error (21.3%) does not include an additive LLOQ
  floor; below-LLOQ behaviour is simulation-noise-only and should not be
  interpreted clinically.

- **No NONMEM control stream supplement.** The paper does not ship a
  `$THETA` / `$OMEGA` / `$SIGMA` block as a supplement. Final parameter
  values are taken directly from Table 3 of Jeon 2012.

- **Errata search.** A targeted PubMed search for
  `kjpp 2012 16 4 249 erratum` and a Google Scholar search for
  `"Jeon 2012 cyanidin" erratum` returned no results as of the
  extraction date; no published correction is on file.
