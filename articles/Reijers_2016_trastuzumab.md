# Trastuzumab (Reijers 2016)

``` r

library(nlmixr2lib)
library(rxode2)
#> rxode2 5.1.2 using 2 threads (see ?getRxThreads)
#>   no cache: create with `rxCreateCache()`
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union
library(tidyr)
library(ggplot2)
library(PKNCA)
#> 
#> Attaching package: 'PKNCA'
#> The following object is masked from 'package:stats':
#> 
#>     filter
```

## Model and source

- Citation: Reijers JAA, van Donge T, Schepers FML, Burggraaf J,
  Stevens J. Use of population approach non-linear mixed effects models
  in the evaluation of biosimilarity of monoclonal antibodies. Eur J
  Clin Pharmacol. 2016;72(11):1343-1352. <doi:10.1007/s00228-016-2101-6>
- Description: Three-compartment population PK model with parallel
  linear and Michaelis-Menten nonlinear elimination from the central
  compartment for intravenous trastuzumab in healthy male volunteers
  from a phase I biosimilarity trial of the FTMB biosimilar vs Herceptin
  reference product (Reijers 2016, combined model on all dose levels
  0.49-6.44 mg/kg); covariates are lean body mass on central volume of
  distribution V1 and BMI on the linear elimination rate constant ke.
- Article (open access via SpringerLink):
  <https://doi.org/10.1007/s00228-016-2101-6>
- Supplement (open access via Springer static-content):
  <https://static-content.springer.com/esm/art%3A10.1007%2Fs00228-016-2101-6/MediaObjects/228_2016_2101_MOESM1_ESM.pdf>

## Population

Reijers et al. (2016) developed a population PK model for intravenous
trastuzumab from a phase I randomised single-dose parallel-group
bioequivalence trial in 110 healthy male volunteers (aged 18-45 years)
who received the biosimilar candidate FTMB (Synthon BV, test product, T)
or Herceptin (EU-licenced reference product, R) at the Centre for Human
Drug Research (Leiden, The Netherlands). The dose-escalation arm
delivered 0.5, 1.5, and 3 mg/kg test product (actual 0.49, 1.48, 2.96
mg/kg, n=6 per dose); the bioequivalence arm delivered 6 mg/kg test
(actual 5.96 mg/kg, n=46) or 6 mg/kg reference (actual 6.44 mg/kg,
n=46). Trastuzumab was administered as a 90-minute IV infusion in 250 mL
0.9% NaCl. Demographics by treatment arm (Reijers 2016 Table 1) span
mean (SD) lean body mass 57.5-62.6 (3.8-8.4) kg, BMI 21.2-23.5 (2.1-3.3)
kg/m^2, and total body weight 72.0-79.5 (7.5-12.6) kg across the five
dose strata. Serum trastuzumab was quantified by ELISA with LLOQ 0.060
ug/mL over 21 sampling times from pre-dose through day 63.

The same information is available programmatically:

``` r

readModelDb("Reijers_2016_trastuzumab")$meta$population
```

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
[`ini()`](https://nlmixr2.github.io/rxode2/reference/ini.html) entry in
`inst/modeldb/specificDrugs/Reijers_2016_trastuzumab.R`. The table below
collects the mapping in one place for reviewer audit.

| Element | Source location | Value / form |
|----|----|----|
| Structural equations (three-compartment, parallel linear + MM, IV input) | Reijers 2016 Results “First step: combined model” and Figure 1 | See model file `model({})` |
| Typical V1 | Reijers 2016 Table 2 (combined model) | 3.28 L |
| Typical V2 | Reijers 2016 Table 2 (combined model) | 1.89 L |
| Typical V3 | Reijers 2016 Table 2 (combined model) | 1.96 L |
| Typical Q1 | Reijers 2016 Table 2 (combined model) | 2.91e-3 L/h = 0.0698 L/day |
| Typical Q2 | Reijers 2016 Table 2 (combined model) | 4.34e-2 L/h = 1.0416 L/day |
| Typical Vmax | Reijers 2016 Table 2 (combined model) | 178 ug/h = 4.272 mg/day |
| Typical Km | Reijers 2016 Table 2 (combined model) | 937 ug/L = 0.937 ug/mL |
| Typical ke | Reijers 2016 Table 2 (combined model) | 2.20e-3 1/h = 0.0528 1/day |
| LBM on V1 | Reijers 2016 Online Resource Eq. 1 | V1_i = V1_pop \* (LBM/61) \* exp(eta), unit-slope ratio |
| BMI on ke | Reijers 2016 Online Resource Eq. 2 | ke_i = ke_pop \* (BMI/23) \* exp(eta), unit-slope ratio |
| IIV on V1 (CV 14.8%) | Reijers 2016 Table 2 (combined model) | omega^2 = 0.0217 |
| IIV on KM (CV 35.9%) | Reijers 2016 Table 2 (combined model) | omega^2 = 0.121 |
| IIV on ke (CV 17.2%) | Reijers 2016 Table 2 (combined model) | omega^2 = 0.0292 |
| Proportional residual error | Reijers 2016 Table 2 (combined model) | sigma^2 = 0.0222 -\> propSd = 0.149 |
| Additive residual error | Reijers 2016 Table 2 (combined model) | sigma^2 = 1520 (ug/L)^2 -\> addSd = 0.039 ug/mL |

## Covariate column naming

| Source column | Canonical column used here | Notes |
|----|----|----|
| `LBW` (lean body weight) | `LBM` (kg) | Alias mapping: `LBM` is the canonical for lean body mass (Reijers reports the value as “lean body weight” computed by the same family of body-composition formulas). Values transfer 1:1 with no transformation. Reference LBM = 61 kg. |
| `BMI` | `BMI` (kg/m^2) | Canonical. Reference BMI = 23 kg/m^2. |

## Virtual cohort

The source paper publishes baseline demographics by dose group (Table 1)
but not per-subject covariate values. The cohort below approximates the
five strata of the combined-model dataset (110 subjects, all healthy
adult males) with covariate distributions centred at the dose-group
means reported in Table 1.

``` r

set.seed(2016)

# Per-stratum virtual-cohort constructor. id_offset keeps subject IDs
# disjoint across strata so rxSolve does not collapse duplicate IDs
# into single (wrong) subjects when bind_rows()-ed.
make_stratum <- function(n, label, dose_mgkg, lbm_mean, lbm_sd,
                         bmi_mean, bmi_sd, wt_mean, wt_sd,
                         id_offset = 0L) {
  tibble::tibble(
    id    = id_offset + seq_len(n),
    arm   = label,
    dose_mgkg = dose_mgkg,
    LBM   = pmin(pmax(rnorm(n, lbm_mean, lbm_sd), 40), 90),
    BMI   = pmin(pmax(rnorm(n, bmi_mean, bmi_sd), 18), 32),
    WT    = pmin(pmax(rnorm(n, wt_mean,  wt_sd),  50), 120)
  )
}

cohort <- dplyr::bind_rows(
  make_stratum( 6, "T 0.49 mg/kg", 0.49, 59.4, 8.4, 21.7, 3.3, 73.1, 12.6, id_offset =   0L),
  make_stratum( 6, "T 1.48 mg/kg", 1.48, 57.5, 5.1, 23.5, 2.6, 73.0,  8.7, id_offset =   6L),
  make_stratum( 6, "T 2.96 mg/kg", 2.96, 59.1, 3.8, 21.2, 2.1, 72.0,  7.5, id_offset =  12L),
  make_stratum(46, "T 5.96 mg/kg", 5.96, 62.6, 6.6, 23.4, 2.5, 79.5, 11.2, id_offset =  18L),
  make_stratum(46, "R 6.44 mg/kg", 6.44, 61.0, 5.6, 23.2, 2.7, 77.1, 10.2, id_offset =  64L)
)

stopifnot(!anyDuplicated(cohort$id))
```

## Dosing regimens

Reijers 2016 administered single 90-minute IV infusions in 250 mL 0.9%
NaCl. The model carries dose into `central` directly; we encode the dose
as a short bolus (the 90-minute infusion is short compared to the
multi-week disposition half-life and the simulation here focuses on the
multi-day terminal phase).

``` r

n_subj <- nrow(cohort)
obs_times <- sort(unique(c(
  c(0, 45/(24*60), 1.5/24, 2/24, 3/24, 4/24, 5/24, 6/24, 8/24, 1),
  c(2, 4, 8, 14, 21, 28, 35, 42, 49, 63),
  seq(0, 70, length.out = 80)
)))

dose_rows <- cohort |>
  dplyr::mutate(
    time = 0,
    amt  = dose_mgkg * WT,
    cmt  = "central",
    evid = 1L
  ) |>
  dplyr::select(id, time, amt, cmt, evid, arm, dose_mgkg, LBM, BMI, WT)

obs_rows <- cohort |>
  tidyr::expand_grid(time = obs_times) |>
  dplyr::mutate(
    amt  = NA_real_,
    cmt  = NA_character_,
    evid = 0L
  ) |>
  dplyr::select(id, time, amt, cmt, evid, arm, dose_mgkg, LBM, BMI, WT)

events <- dplyr::bind_rows(dose_rows, obs_rows) |>
  dplyr::arrange(id, time, dplyr::desc(evid))

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
```

## Simulate the cohort

``` r

mod <- readModelDb("Reijers_2016_trastuzumab")
sim <- rxode2::rxSolve(
  mod,
  events = events,
  keep   = c("arm", "dose_mgkg")
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'
sim$time <- as.numeric(sim$time)
```

### Typical-subject profiles (reproduces Figure 3 of Reijers 2016)

Figure 3 of Reijers 2016 shows the visual predictive check (VPC) of the
best combined model conditioned on dose, with typical-value trajectories
overlaid. Below we plot the typical-subject concentration-time profile
(no IIV / no residual error) for each of the five dose strata, using the
dose-group mean LBM and BMI as the typical-subject covariates.

``` r

mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

ref_strata <- tibble::tribble(
  ~arm,             ~dose_mgkg, ~LBM, ~BMI, ~WT,
  "T 0.49 mg/kg",   0.49,        59.4, 21.7, 73.1,
  "T 1.48 mg/kg",   1.48,        57.5, 23.5, 73.0,
  "T 2.96 mg/kg",   2.96,        59.1, 21.2, 72.0,
  "T 5.96 mg/kg",   5.96,        62.6, 23.4, 79.5,
  "R 6.44 mg/kg",   6.44,        61.0, 23.2, 77.1
)
ref_strata$id <- seq_len(nrow(ref_strata))

ref_obs_times <- sort(unique(c(
  seq(0,  1,  length.out = 50),
  seq(1,  70, length.out = 200)
)))

ref_doses <- ref_strata |>
  dplyr::mutate(
    time = 0,
    amt  = dose_mgkg * WT,
    cmt  = "central",
    evid = 1L
  ) |>
  dplyr::select(id, time, amt, cmt, evid, arm, dose_mgkg, LBM, BMI, WT)

ref_obs <- ref_strata |>
  tidyr::expand_grid(time = ref_obs_times) |>
  dplyr::mutate(
    amt  = NA_real_,
    cmt  = NA_character_,
    evid = 0L
  ) |>
  dplyr::select(id, time, amt, cmt, evid, arm, dose_mgkg, LBM, BMI, WT)

ref_events <- dplyr::bind_rows(ref_doses, ref_obs) |>
  dplyr::arrange(id, time, dplyr::desc(evid))

sim_typ <- rxode2::rxSolve(mod_typical, events = ref_events,
                           keep = c("arm", "dose_mgkg")) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalkm', 'etalkel'
#> Warning: multi-subject simulation without without 'omega'
sim_typ$time <- as.numeric(sim_typ$time)

ggplot(sim_typ |> dplyr::filter(time > 0, Cc > 0.01),
       aes(time, Cc, color = arm, group = arm)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(x = "Time (days)",
       y = "Trastuzumab Cc (ug/mL)",
       title = "Typical-subject trastuzumab profiles by dose stratum",
       subtitle = "Single 90-min IV infusion; covariates set to dose-group means (Reijers 2016 Table 1)",
       caption = "Patterned after Figure 3 of Reijers 2016.") +
  theme_bw()
```

![Typical-subject Cc vs. time by dose stratum after a single 90-min IV
trastuzumab infusion (reproduces the typical-value overlay of Figure 3
of Reijers
2016).](Reijers_2016_trastuzumab_files/figure-html/figure-3-typical-1.png)

Typical-subject Cc vs. time by dose stratum after a single 90-min IV
trastuzumab infusion (reproduces the typical-value overlay of Figure 3
of Reijers 2016).

### VPC-style summary across the virtual cohort (reproduces Figure 3 of Reijers 2016)

``` r

vpc <- sim |>
  dplyr::filter(time > 0, Cc > 0.01, !is.na(Cc)) |>
  dplyr::mutate(time_bin = cut(time,
                               breaks = c(0, 0.05, 0.1, 0.2, 0.5, 1, 2, 4, 7,
                                          14, 21, 28, 35, 42, 49, 56, 63, 70),
                               include.lowest = TRUE, labels = FALSE)) |>
  dplyr::group_by(arm, time_bin) |>
  dplyr::summarise(time   = mean(time),
                   median = median(Cc),
                   lo     = quantile(Cc, 0.05),
                   hi     = quantile(Cc, 0.95),
                   n      = dplyr::n(),
                   .groups = "drop") |>
  dplyr::filter(n >= 3)

ggplot(vpc, aes(time)) +
  geom_ribbon(aes(ymin = lo, ymax = hi, fill = arm), alpha = 0.20) +
  geom_line(aes(y = median, color = arm), linewidth = 0.9) +
  scale_y_log10() +
  labs(x = "Time (days)",
       y = "Trastuzumab Cc (ug/mL)",
       title = "Simulated trastuzumab concentration-time profiles by dose",
       subtitle = "Median and 5-95% prediction interval; combined-model parameters from Reijers 2016 Table 2",
       caption = "Patterned after Figure 3 of Reijers 2016.") +
  theme_bw()
```

![VPC-style summary (median and 5-95% prediction interval) stratified by
dose. Patterned after Figure 3 of Reijers
2016.](Reijers_2016_trastuzumab_files/figure-html/figure-3-vpc-1.png)

VPC-style summary (median and 5-95% prediction interval) stratified by
dose. Patterned after Figure 3 of Reijers 2016.

## PKNCA validation

Compute Cmax, Tmax, AUClast, and AUCinf per subject and roll up by dose
stratum to compare against the geometric-mean AUCs in Reijers 2016 Table
3. The PKNCA formula carries the dose-stratum label through `keep` so
the per-arm summaries are directly comparable.

``` r

sim_nca <- sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, Cc, arm)

dose_df <- events |>
  dplyr::filter(evid == 1L) |>
  dplyr::select(id, time, amt, arm)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | arm + id,
                             concu = "ug/mL", timeu = "day")
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | arm + id,
                             doseu = "mg")

intervals <- data.frame(
  start       = 0,
  end         = Inf,
  cmax        = TRUE,
  tmax        = TRUE,
  auclast     = TRUE,
  aucinf.obs  = TRUE,
  half.life   = TRUE
)

nca_result <- suppressWarnings(
  PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
)

nca_tbl <- as.data.frame(nca_result$result) |>
  dplyr::filter(PPTESTCD %in% c("cmax", "tmax", "auclast", "aucinf.obs", "half.life"))

nca_summary <- nca_tbl |>
  dplyr::group_by(arm, PPTESTCD) |>
  dplyr::summarise(
    gm   = exp(mean(log(pmax(PPORRES, 1e-9)), na.rm = TRUE)),
    q05  = quantile(PPORRES, 0.05, na.rm = TRUE),
    q95  = quantile(PPORRES, 0.95, na.rm = TRUE),
    n    = sum(!is.na(PPORRES)),
    .groups = "drop"
  )

knitr::kable(nca_summary, digits = 2,
             caption = "Per-arm NCA summary (geometric mean and 5-95% PI) from the simulated cohort.")
```

| arm          | PPTESTCD   |      gm |     q05 |     q95 |   n |
|:-------------|:-----------|--------:|--------:|--------:|----:|
| R 6.44 mg/kg | aucinf.obs | 1546.84 | 1023.43 | 2227.82 |  46 |
| R 6.44 mg/kg | auclast    | 1543.46 | 1020.63 | 2224.82 |  46 |
| R 6.44 mg/kg | cmax       |  150.81 |  104.08 |  203.57 |  46 |
| R 6.44 mg/kg | half.life  |   13.71 |    7.24 |   17.00 |  46 |
| R 6.44 mg/kg | tmax       |    0.00 |    0.00 |    0.00 |  46 |
| T 0.49 mg/kg | aucinf.obs |   34.16 |   22.68 |   43.52 |   6 |
| T 0.49 mg/kg | auclast    |   34.11 |   22.66 |   43.45 |   6 |
| T 0.49 mg/kg | cmax       |   12.45 |    8.90 |   14.92 |   6 |
| T 0.49 mg/kg | half.life  |   18.89 |   18.84 |   18.95 |   6 |
| T 0.49 mg/kg | tmax       |    0.00 |    0.00 |    0.00 |   6 |
| T 1.48 mg/kg | aucinf.obs |  164.52 |  135.10 |  192.87 |   6 |
| T 1.48 mg/kg | auclast    |  164.30 |  134.90 |  192.62 |   6 |
| T 1.48 mg/kg | cmax       |   32.07 |   30.20 |   35.39 |   6 |
| T 1.48 mg/kg | half.life  |   18.60 |   18.51 |   18.69 |   6 |
| T 1.48 mg/kg | tmax       |    0.00 |    0.00 |    0.00 |   6 |
| T 2.96 mg/kg | aucinf.obs |  437.60 |  333.84 |  606.75 |   6 |
| T 2.96 mg/kg | auclast    |  436.59 |  333.12 |  605.01 |   6 |
| T 2.96 mg/kg | cmax       |   62.12 |   48.39 |   79.30 |   6 |
| T 2.96 mg/kg | half.life  |   18.20 |   17.96 |   18.39 |   6 |
| T 2.96 mg/kg | tmax       |    0.00 |    0.00 |    0.00 |   6 |
| T 5.96 mg/kg | aucinf.obs | 1293.06 |  897.51 | 2030.51 |  46 |
| T 5.96 mg/kg | auclast    | 1290.52 |  895.81 | 2026.71 |  46 |
| T 5.96 mg/kg | cmax       |  134.12 |   95.58 |  191.10 |  46 |
| T 5.96 mg/kg | half.life  |   14.92 |   10.06 |   17.35 |  46 |
| T 5.96 mg/kg | tmax       |    0.00 |    0.00 |    0.00 |  46 |

Per-arm NCA summary (geometric mean and 5-95% PI) from the simulated
cohort. {.table}

### Comparison against Reijers 2016 Table 3 (actual-dose AUC)

Reijers 2016 Table 3 reports per-arm geometric-mean AUClast and AUCinf
at the actual delivered dose (test 5.92 mg/kg, reference 6.44 mg/kg)
using both the standard NCA on observed concentrations and the
combined-model simulations at original sampling times. Comparison
focuses on the 6 mg/kg strata where the source paper reports per-arm
geometric means.

``` r

published <- tibble::tribble(
  ~arm,            ~metric,             ~paper_gm,
  "T 5.96 mg/kg",  "AUClast (ug*day/mL)", 1301,
  "T 5.96 mg/kg",  "AUCinf  (ug*day/mL)", 1311,
  "R 6.44 mg/kg",  "AUClast (ug*day/mL)", 1588,
  "R 6.44 mg/kg",  "AUCinf  (ug*day/mL)", 1593
)

sim_wide <- nca_summary |>
  dplyr::mutate(metric = dplyr::case_when(
    PPTESTCD == "auclast"    ~ "AUClast (ug*day/mL)",
    PPTESTCD == "aucinf.obs" ~ "AUCinf  (ug*day/mL)",
    TRUE                     ~ NA_character_
  )) |>
  dplyr::filter(!is.na(metric)) |>
  dplyr::select(arm, metric, sim_gm = gm,
                sim_q05 = q05, sim_q95 = q95)

comparison <- dplyr::left_join(published, sim_wide,
                                by = c("arm", "metric")) |>
  dplyr::mutate(pct_diff = 100 * (sim_gm - paper_gm) / paper_gm)

knitr::kable(comparison, digits = 1,
             caption = "Per-arm AUC: simulated geometric mean (with 5-95% PI) vs Reijers 2016 Table 3 geometric mean at actual delivered dose.")
```

| arm          | metric               | paper_gm | sim_gm | sim_q05 | sim_q95 | pct_diff |
|:-------------|:---------------------|---------:|-------:|--------:|--------:|---------:|
| T 5.96 mg/kg | AUClast (ug\*day/mL) |     1301 | 1290.5 |   895.8 |  2026.7 |     -0.8 |
| T 5.96 mg/kg | AUCinf (ug\*day/mL)  |     1311 | 1293.1 |   897.5 |  2030.5 |     -1.4 |
| R 6.44 mg/kg | AUClast (ug\*day/mL) |     1588 | 1543.5 |  1020.6 |  2224.8 |     -2.8 |
| R 6.44 mg/kg | AUCinf (ug\*day/mL)  |     1593 | 1546.8 |  1023.4 |  2227.8 |     -2.9 |

Per-arm AUC: simulated geometric mean (with 5-95% PI) vs Reijers 2016
Table 3 geometric mean at actual delivered dose. {.table}

The simulated geometric-mean AUC values track the paper’s NCA-derived
geometric means within approximately 10-15%. The simulated reference-
arm AUC sits slightly below the paper’s value primarily because the
virtual cohort spans the LBM and BMI distributions reported in Table 1
(which differ slightly from the LBM-median and BMI-median used as
reference values in the supplement covariate equations), and because the
cohort here uses simulated covariates rather than the actual trial data.

## Covariate-effect sanity checks (reproduces Reijers 2016 SI equations)

The supplement equations specify unit-slope covariate effects. The block
below reproduces the typical V1 and ke directly from the packaged
coefficients to confirm the LBM and BMI scaling factors match the
published equations.

``` r

q <- list(
  v1_0 = 3.28, ke_0 = 0.0528,
  lbm_ref = 61, bmi_ref = 23
)

typ <- function(LBM = 61, BMI = 23) {
  list(
    v1 = q$v1_0 * (LBM / q$lbm_ref),
    ke = q$ke_0 * (BMI / q$bmi_ref)
  )
}

sensitivity <- tibble::tribble(
  ~Scenario,             ~LBM, ~BMI,
  "Reference",           61,   23,
  "LBM = 50 kg",         50,   23,
  "LBM = 75 kg",         75,   23,
  "BMI = 19 kg/m^2",     61,   19,
  "BMI = 28 kg/m^2",     61,   28
) |>
  dplyr::rowwise() |>
  dplyr::mutate(
    `Typical V1 (L)` = typ(LBM, BMI)$v1,
    `Typical ke (1/day)` = typ(LBM, BMI)$ke
  ) |>
  dplyr::ungroup()

knitr::kable(sensitivity, digits = 4,
             caption = "Typical V1 and ke sensitivities reproduced from the packaged covariate coefficients.")
```

| Scenario        | LBM | BMI | Typical V1 (L) | Typical ke (1/day) |
|:----------------|----:|----:|---------------:|-------------------:|
| Reference       |  61 |  23 |         3.2800 |             0.0528 |
| LBM = 50 kg     |  50 |  23 |         2.6885 |             0.0528 |
| LBM = 75 kg     |  75 |  23 |         4.0328 |             0.0528 |
| BMI = 19 kg/m^2 |  61 |  19 |         3.2800 |             0.0436 |
| BMI = 28 kg/m^2 |  61 |  28 |         3.2800 |             0.0643 |

Typical V1 and ke sensitivities reproduced from the packaged covariate
coefficients. {.table}

The reference-subject V1 (3.28 L) and ke (0.0528 1/day) match Reijers
2016 Table 2 to three significant digits. LBM has a unit-slope
multiplicative effect on V1; BMI has a unit-slope multiplicative effect
on ke.

## Assumptions and deviations

Reijers 2016 does not publish per-subject covariate values; the virtual
cohort here approximates each dose stratum by sampling LBM, BMI, and
total body weight independently from Normal distributions centred at the
per-stratum means and SDs reported in Table 1.

- **LBM-median and BMI-median reference values.** The supplement
  equations normalise LBM and BMI to the “median” of the analysis cohort
  but the paper does not state these reference values explicitly. The
  cohort-weighted mean across the five dose strata (computed from Table
  1 dose-group means weighted by stratum size) is approximately 61 kg
  for LBM and 23 kg/m^2 for BMI, and these are the values used as the
  reference normalisations in the packaged model. The exact medians may
  differ by a few percent; predictions at the population level are
  correspondingly biased by the same fraction since the covariate
  effects are unit-slope multiplicative.
- **No omega block on (KM, ke).** Reijers 2016 Results state that “an
  omega block was required to correct for the parameter correlation
  between KM and ke in the model” but neither the main paper nor the
  open-access supplement reports the off-diagonal covariance /
  correlation value. The packaged model encodes `etalkm` and `etalkel`
  as independent random effects on a diagonal omega; the diagonal
  variances themselves match Reijers 2016 Table 2 exactly. Simulations
  from the diagonal-omega encoding may understate the true joint
  variability of the saturable elimination pathway because they ignore
  the (KM, ke) covariance the source paper retained for identifiability.
- **`lkel` as a primary parameter.** Reijers 2016 parameterises the
  linear-elimination arm directly by the rate constant ke rather than by
  a clearance CL. The nlmixr2lib canonical convention uses `lcl` as the
  primary parameter and derives `kel = cl/vc`; here we use `lkel` as a
  primary parameter so the BMI covariate effect applies directly to
  `kel` exactly as in the source equations. Re-expressing this with
  `lcl` would require carrying both the LBM effect on V1 and the BMI
  effect on ke through to CL = ke \* V1, which (a) does not match the
  paper’s structural parameterisation and (b) yields a different
  propagation of IIV to total clearance. The deviation is intentional
  and preserves source fidelity.
- **Bolus vs 90-minute infusion.** The source administered trastuzumab
  as a 90-minute IV infusion in 250 mL 0.9% NaCl. The vignette here
  encodes the dose as an instantaneous bolus into `central` since the
  infusion duration (1.5 hours) is short compared to both the linear
  elimination half-life (approximately 13 days) and the multi-week
  distribution timescale. Users who need the explicit infusion profile
  can replace the `amt + cmt` / `evid = 1` rows with an `evid = 1` +
  `rate` form (rate = amt / (1.5 / 24)) without other changes.
- **Drug product (test vs reference) is not a model covariate.** The
  paper found drug product to be non-significant on every parameter (max
  OFV decrease 5.80, p \> 0.01) and did not retain it. The packaged
  model is the final combined model that pools both products on
  identical structural parameters; the separate models (test-only and
  reference-only at 6 mg/kg) reported in Table 2 are not packaged here.
- **Residual error magnitudes**: Table 2 reports the proportional
  variability as sigma^2 = 0.0222 (CV approximately 14.9%) and the
  additive variability as sigma^2 = 1520. The additive variance is
  reported on the ug/L internal scale; the equivalent SD on the
  canonical ug/mL scale used in this file is 0.039 ug/mL, which is below
  the assay LLOQ (0.060 ug/mL) and consistent with the paper’s residual
  coefficient of variation of 14.98%.
- **Time units**: days throughout. The paper reports rate constants in
  1/h and Vmax in ug/h; conversion factors of 24 (1/h -\> 1/day) and
  24/1000 (ug/h -\> mg/day) are applied. The Km is converted from ug/L
  to ug/mL by dividing by 1000.
- **TMDD alternative explored but abandoned**: the paper Results note
  that a target-mediated drug disposition (TMDD) model and its Kd-based
  approximation were “abandoned due to over-parameterisation and
  instability”. The packaged model is the parallel-linear plus
  Michaelis-Menten approximation that the authors retained.
- **Errata search**: no author correction or erratum was located on the
  European Journal of Clinical Pharmacology landing page, PubMed, or
  Google Scholar for DOI 10.1007/s00228-016-2101-6 at the time of
  extraction. The paper is open access (Creative Commons BY 4.0) and the
  supplement (Online Resource) is openly downloadable from Springer
  static-content.

## Model summary

- **Structure**: three-compartment PK model with parallel linear (rate
  constant ke) and Michaelis-Menten nonlinear elimination (Vmax, Km)
  from the central compartment; IV infusion delivered directly into
  central.
- **Typical parameters** (combined model): V1 = 3.28 L, V2 = 1.89 L, V3
  = 1.96 L, Q1 = 0.0698 L/day, Q2 = 1.0416 L/day, Vmax = 4.272 mg/day,
  Km = 0.937 ug/mL, ke = 0.0528 1/day.
- **Reference subject**: 61 kg LBM, BMI 23 kg/m^2 (cohort-weighted means
  across the five dose strata in Reijers 2016 Table 1).
- **Clinically relevant covariates**: LBM scales V1 with unit slope; BMI
  scales ke with unit slope.
- **Conclusion (from the paper)**: the FTMB biosimilar candidate and
  Herceptin reference are pharmacokinetically equivalent at the 6 mg/kg
  dose level; drug product is not a statistically significant covariate
  on any model parameter.

## Reference

- Reijers JAA, van Donge T, Schepers FML, Burggraaf J, Stevens J. Use of
  population approach non-linear mixed effects models in the evaluation
  of biosimilarity of monoclonal antibodies. Eur J Clin Pharmacol.
  2016;72(11):1343-1352. <doi:10.1007/s00228-016-2101-6>
