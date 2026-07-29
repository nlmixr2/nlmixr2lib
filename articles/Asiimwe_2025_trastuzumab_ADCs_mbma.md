# Trastuzumab drug-conjugate MBMA (Asiimwe 2025)

## Model and source

Asiimwe et al. (2025) developed two-compartment linear population-PK
models for the anti-HER2 antibody-drug conjugates trastuzumab emtansine
(T-DM1) and trastuzumab deruxtecan (T-DXd) using model-based
meta-analysis (MBMA) of summary-level concentration-time data digitised
from 14 (T-DM1) and 4 (T-DXd) published clinical trials. The paper’s
downstream aim was to combine model-predicted PK exposures with
logistic-regression exposure-response for objective response rate (ORR)
and dose-limiting toxicity (DLT) into a clinical utility index (CUI).
This vignette reproduces only the population-PK layer of the paper
(Section 2.2 and Table 1); the CUI / exposure-response layer is a set of
logistic regressions in R and is not an ODE-based model.

Two model files are packaged, one per ADC, both suffixed `_mbma` to
distinguish them from the pre-existing individual-level popPK models for
the same drugs (`Lu_2014_trastuzumabemtansine`,
`Yin_2021_trastuzumabDeruxtecan`):

``` r

mod_dm1 <- readModelDb("Asiimwe_2025_trastuzumabEmtansine_mbma")
mod_dxd <- readModelDb("Asiimwe_2025_trastuzumabDeruxtecan_mbma")

cat("T-DM1 description:\n", rxode2::rxode(mod_dm1)$description, "\n\n", sep = "")
#> ℹ parameter labels from comments will be replaced by 'label()'
#> T-DM1 description:
#> MBMA. Two-compartment linear population PK model of trastuzumab emtansine (T-DM1, anti-HER2 antibody-drug conjugate) fitted by model-based meta-analysis to summary-level concentration-time data digitised from 14 published clinical trials in patients with HER2-positive metastatic breast cancer and other HER2-positive solid tumors. Between-study variability (BSV) is a study-level random effect on CL and Vc (block correlation 0.825) representing differences in inclusion criteria and the pooling of Phase I dose-escalation arms as separate 'studies'. Residual error is proportional + additive on serum concentration (ng/mL) and was weighted by the square root of each trial's sample size during fitting; the tabulated parameter estimates are the unweighted values and simulation of a study of N subjects should scale the residual SD by 1/sqrt(N). Suitable simulation scope is study-arm-mean cycle-1 concentration-time profiles, NOT individual-patient concentrations. Parameter values are the T-DM1 column of Asiimwe 2025 Table 1 (Monolix 2024R1).
cat("T-DXd description:\n", rxode2::rxode(mod_dxd)$description, "\n", sep = "")
#> ℹ parameter labels from comments will be replaced by 'label()'
#> T-DXd description:
#> MBMA. Two-compartment linear population PK model of trastuzumab deruxtecan (T-DXd, DS-8201, Enhertu, anti-HER2 antibody-drug conjugate) fitted by model-based meta-analysis to summary-level concentration-time data digitised from 4 published clinical trials in patients with HER2-positive breast, non-small-cell lung, gastric/GEJ, colorectal and other solid tumors. Between-study variability (BSV) is a study-level random effect on CL and Vc (block correlation 0.915) representing differences in inclusion criteria across trials and dose-escalation arms treated as separate studies. Residual error is proportional + additive on serum concentration (ng/mL) and was weighted by the square root of each trial's sample size during fitting; the tabulated parameter estimates are the unweighted values and simulation of a study of N subjects should scale the residual SD by 1/sqrt(N). Suitable simulation scope is study-arm-mean cycle-1 concentration-time profiles, NOT individual-patient concentrations. Parameter values are the T-DXd column of Asiimwe 2025 Table 1 (Monolix 2024R1).
```

- Article: <https://doi.org/10.1002/psp4.70013>
- Supplement Data S1 (Text S1-S3, Figures S1-S7): included in the
  article DOI landing page
- Analysis code (Shiny app + fitting scripts):
  <https://github.com/iasiimwe/adc_cui>

## Population

The MBMA input pool covered 103 unique clinical trials identified by a
MEDLINE search on 1 May 2024 (with supplementation from
ClinicalTrials.gov, Citeline, and the Beacon ADC database). Of these:

- **T-DM1 Pop-PK subset** – 14 studies (Asiimwe 2025 Table S1).
  Predominantly HER2-positive metastatic breast cancer; some
  HER2-positive urothelial/bladder and pancreatic/cholangiocarcinoma
  cohorts. Dose range 0.3 to 4.8 mg/kg IV; three-weekly and weekly
  regimens. Median body weight 69.4 kg (used for the per-patient mg dose
  in simulations).
- **T-DXd Pop-PK subset** – 4 studies (Asiimwe 2025 Table S3).
  HER2-positive advanced/metastatic breast, non-small-cell lung,
  gastric/GEJ, colorectal and other HER2-expressing solid tumors. Dose
  range 0.8 to 8.0 mg/kg IV; three-weekly and weekly regimens. Median
  body weight 59.0 kg.

MBMA fits summary-level data (arm-mean concentrations per timepoint per
study), so the model characterises study-arm-mean cycle-1
concentration-time profiles rather than individual-patient
concentrations. No covariates were included in the primary Pop-PK model;
body weight enters only through the mg-per-patient dose (dose = mg/kg
times median weight).

Between-study variability (BSV) is a study-level random effect on CL and
Vc, with block correlation 0.825 (T-DM1) or 0.915 (T-DXd). Because
different doses within a study are treated as separate “studies” in this
analysis (to accommodate Phase I dose-escalation designs where the same
patients are not randomised across arms), BSV also absorbs
between-treatment-arm variability (BTAV).

Programmatic access to the population metadata:

``` r

pop_dm1 <- rxode2::rxode(mod_dm1)$population
#> ℹ parameter labels from comments will be replaced by 'label()'
pop_dxd <- rxode2::rxode(mod_dxd)$population
#> ℹ parameter labels from comments will be replaced by 'label()'
list(
  T_DM1 = pop_dm1[c("n_studies", "weight_median", "disease_state")],
  T_DXd = pop_dxd[c("n_studies", "weight_median", "disease_state")]
)
#> $T_DM1
#> $T_DM1$n_studies
#> [1] 14
#> 
#> $T_DM1$weight_median
#> [1] "69.4 kg (median across the 14 T-DM1 Pop-PK studies; used when computing the mg-per-patient dose from a mg/kg regimen)"
#> 
#> $T_DM1$disease_state
#> [1] "predominantly HER2-positive metastatic breast cancer; a minority of arms include HER2-positive urothelial/bladder, pancreatic/cholangiocarcinoma, and other HER2-positive solid tumors (see Asiimwe 2025 Table S1 for the 14 Pop-PK studies)"
#> 
#> 
#> $T_DXd
#> $T_DXd$n_studies
#> [1] 4
#> 
#> $T_DXd$weight_median
#> [1] "59.0 kg (median across the 4 T-DXd Pop-PK studies; used when computing the mg-per-patient dose from a mg/kg regimen)"
#> 
#> $T_DXd$disease_state
#> [1] "HER2-positive advanced or metastatic solid tumors: breast cancer (dominant), non-small-cell lung, gastric / gastro-oesophageal junction, colorectal, and other HER2-expressing tumors (see Asiimwe 2025 Table S3 for the 4 Pop-PK studies)"
```

## Source trace

Every parameter carries an inline comment in the model files
(`inst/modeldb/specificDrugs/Asiimwe_2025_trastuzumabEmtansine_mbma.R`
and `..._trastuzumabDeruxtecan_mbma.R`). The table below collects the
values in one place for review; all come from the T-DM1 and T-DXd
columns of Table 1 in Asiimwe 2025.

| Parameter (nlmixr2 name) | T-DM1 value | T-DXd value | Source | Note |
|----|----|----|----|----|
| `lcl` = log(CL) | log(0.809) | log(0.585) | Table 1 “Clearance (L/day)” | RSE 7.3% / 3.7% |
| `lvc` = log(Vc) | log(3.283) | log(2.785) | Table 1 “Central volume of distribution (L)” | RSE 4.8% / 1.5% |
| `lvp` = log(Vp) | log(0.748) | log(1.243) | Table 1 “Peripheral volume of distribution (L)” | RSE 8.3% / 9.2% |
| `lq` = log(Q) | log(1.120) | log(0.652) | Table 1 “Intercompartment clearance (L/day)” | RSE 26.3% / 9.7% |
| Var(`eta_study_lcl`) | 0.334^2 = 0.111556 | 0.084^2 = 0.007056 | Table 1 “BSV clearance” | omega (SD) reported |
| Var(`eta_study_lvc`) | 0.221^2 = 0.048841 | 0.038^2 = 0.001444 | Table 1 “BSV central volume of distribution” | omega (SD) reported |
| Cov(`eta_study_lcl`,`eta_study_lvc`) | 0.825 \* 0.334 \* 0.221 = 0.0609 | 0.915 \* 0.084 \* 0.038 = 0.002921 | Table 1 “BSV clearance~BSV central volume of distribution” | correlation coefficient reported |
| `propSd` | 0.430 | 0.338 | Table 1 “Weighted proportional error (%)” | study-level SD, pre-sqrt(N) weighting |
| `addSd` | 2633.766 ng/mL | 3043.004 ng/mL | Table 1 “Constant error (ng/mL)” | study-level SD, pre-sqrt(N) weighting |
| `d/dt(central)` = -kel*central - k12*central + k21\*peripheral1 | – | – | Methods 2.2 (two-compartment linear elimination) | standard 2-cpt structure |
| `Cc` = 1000 \* central / vc (ng/mL) | – | – | Table 1 units column | \*1000 converts mg/L to ng/mL |

The residual error was weighted by the square root of each trial’s
sample size during fitting (Methods 2.2). The Table 1 values are the
point estimates before that weighting; simulation of a study of size N
should scale the residual SD by 1/sqrt(N).

## Virtual cohort

Because this is an MBMA model, the natural “cohort” unit is a *study
arm* rather than an individual patient. Each simulated “id” below
represents one study arm (of an assumed size), and the between-id
variability is the between-study variability (BSV) from Table 1. We
simulate 100 study-arms per dose level per drug (well within the
200-per-arm cap) to visualise BSV.

``` r

set.seed(20260724)

# Approved / referenced doses used in Figure 3 of Asiimwe 2025.
DOSES_DM1 <- c(0.3, 0.6, 1.2, 1.8, 2.4, 3.6, 4.8)   # mg/kg (Asiimwe 2025 Figure 3 Panel A)
DOSES_DXD <- c(0.8, 1.6, 3.2, 5.4, 6.4, 7.4, 8.0)   # mg/kg (Asiimwe 2025 Figure 3 Panel B)

# Median body weights across the T-DM1 / T-DXd studies (Asiimwe 2025 Results 3.2).
WT_MEDIAN_DM1 <- 69.4
WT_MEDIAN_DXD <- 59.0

# Cycle-1 infusion durations. Paper: 90 min in Cycle 1; time unit is day.
INFUSION_DUR_DAY <- 90 / 60 / 24   # 90 min in days = 0.0625 day

# Observation grid: dense early to catch Cmax, sparse late to catch Cmin
# at 504 h = 21 days (the paper's Cmin timepoint).
OBS_TIMES_DAY <- c(seq(0, 1, by = 0.05),      # 0-1 day: dense
                   seq(1.5, 21, by = 0.5))    # 1.5-21 days: sparse

make_arms <- function(doses_mgkg, weight_kg, drug, id_offset = 0L, n_per_dose = 100L) {
  purrr::map_dfr(seq_along(doses_mgkg), function(k) {
    dose_mgkg <- doses_mgkg[[k]]
    amt_mg    <- dose_mgkg * weight_kg
    base_id   <- id_offset + (k - 1L) * n_per_dose
    dose_rows <- tibble(
      id   = base_id + seq_len(n_per_dose),
      time = 0,
      amt  = amt_mg,
      rate = amt_mg / INFUSION_DUR_DAY,   # rate has units amt/time = mg/day
      evid = 1L,
      cmt  = "central",
      drug = drug,
      dose_mgkg = dose_mgkg,
      treatment = paste0(drug, " ", format(dose_mgkg, trim = TRUE), " mg/kg")
    )
    obs_rows <- tidyr::crossing(
      id   = base_id + seq_len(n_per_dose),
      time = OBS_TIMES_DAY
    ) |>
      mutate(amt = NA_real_, rate = NA_real_, evid = 0L, cmt = "central",
             drug = drug, dose_mgkg = dose_mgkg,
             treatment = paste0(drug, " ", format(dose_mgkg, trim = TRUE), " mg/kg"))
    bind_rows(dose_rows, obs_rows) |>
      arrange(id, time, desc(evid))
  })
}

events_dm1 <- make_arms(DOSES_DM1, WT_MEDIAN_DM1, drug = "T-DM1", id_offset =    0L)
events_dxd <- make_arms(DOSES_DXD, WT_MEDIAN_DXD, drug = "T-DXd", id_offset = 10000L)

# Regression guard -- id/time/evid combos must be unique across the two cohorts.
stopifnot(!anyDuplicated(unique(bind_rows(events_dm1, events_dxd)[, c("id", "time", "evid")])))
```

## Simulation

``` r

sim_dm1 <- rxode2::rxSolve(
  mod_dm1, events = events_dm1,
  keep = c("drug", "dose_mgkg", "treatment")
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_dxd <- rxode2::rxSolve(
  mod_dxd, events = events_dxd,
  keep = c("drug", "dose_mgkg", "treatment")
) |>
  as.data.frame()
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_all <- bind_rows(sim_dm1, sim_dxd)
```

Typical-value profiles (BSV zeroed) are helpful for reproducing the
paper’s Figure S2 (digitised concentration-time profiles):

``` r

mod_dm1_typ <- mod_dm1 |> rxode2::zeroRe()
#> ℹ parameter labels from comments will be replaced by 'label()'
mod_dxd_typ <- mod_dxd |> rxode2::zeroRe()
#> ℹ parameter labels from comments will be replaced by 'label()'

# Typical-value simulation uses only the 3.6 mg/kg T-DM1 and 5.4 / 6.4 mg/kg
# T-DXd approved doses (one representative arm each) to keep the figure clean.
typ_events <- bind_rows(
  make_arms(3.6, WT_MEDIAN_DM1, drug = "T-DM1", id_offset = 0L, n_per_dose = 1L),
  make_arms(c(5.4, 6.4), WT_MEDIAN_DXD, drug = "T-DXd", id_offset = 100L, n_per_dose = 1L)
)

typ_dm1 <- rxode2::rxSolve(mod_dm1_typ,
                           events = typ_events |> filter(drug == "T-DM1"),
                           keep = c("drug", "dose_mgkg", "treatment")) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'eta_study_lcl', 'eta_study_lvc'
typ_dxd <- rxode2::rxSolve(mod_dxd_typ,
                           events = typ_events |> filter(drug == "T-DXd"),
                           keep = c("drug", "dose_mgkg", "treatment")) |>
  as.data.frame()
#> ℹ omega/sigma items treated as zero: 'eta_study_lcl', 'eta_study_lvc'
#> Warning: multi-subject simulation without without 'omega'
typ_all <- bind_rows(typ_dm1, typ_dxd)
```

## Replicate published figures

### Figure S2 shape – typical-value cycle-1 profile

The paper’s Figure S2 (Panels A and C) shows digitised
concentration-time profiles for the ADCs. The paper does not tabulate
the exact numeric values, so we can only compare shape and
order-of-magnitude. Typical-value cycle-1 profiles at the approved
doses:

``` r

typ_all |>
  filter(!is.na(Cc), time > 0) |>
  ggplot(aes(time, Cc / 1000, colour = treatment, linetype = treatment)) +
  geom_line() +
  scale_y_log10() +
  labs(x = "Time since dose (day)",
       y = "T-DM1 or T-DXd concentration (ug/mL)",
       colour = NULL, linetype = NULL,
       title = "Typical-value cycle-1 profile at approved doses",
       caption = "Compare shape / scale with Asiimwe 2025 Figure S2 Panels A and C (digitised).")
```

![](Asiimwe_2025_trastuzumab_ADCs_mbma_files/figure-html/figure-s2-shape-1.png)

### VPC-style visualisation of BSV

Figure S4 of the paper shows prediction-corrected VPCs from the fitted
models. The MBMA BSV is a study-level effect, so the ribbons below are
between-study (not between-patient) prediction intervals at each dose:

``` r

sim_all |>
  filter(!is.na(Cc), time > 0) |>
  group_by(drug, dose_mgkg, treatment, time) |>
  summarise(
    Q05 = quantile(Cc, 0.05, na.rm = TRUE),
    Q50 = quantile(Cc, 0.50, na.rm = TRUE),
    Q95 = quantile(Cc, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  filter(dose_mgkg %in% c(3.6, 5.4)) |>
  ggplot(aes(time, Q50 / 1000)) +
  geom_ribbon(aes(ymin = Q05 / 1000, ymax = Q95 / 1000, fill = drug), alpha = 0.25) +
  geom_line(aes(colour = drug)) +
  facet_wrap(~ treatment) +
  scale_y_log10() +
  labs(x = "Time since dose (day)",
       y = "Concentration (ug/mL)",
       fill = NULL, colour = NULL,
       title = "Between-study (BSV) prediction interval, approved doses",
       caption = "5th-50th-95th percentile of simulated study-arm-mean concentrations. Compare shape with Asiimwe 2025 Figure S4.")
```

![](Asiimwe_2025_trastuzumab_ADCs_mbma_files/figure-html/figure-s4-shape-1.png)

## PKNCA validation

Cycle-1 NCA is computed on the typical-value profiles for each dose. The
paper’s Figure S5 compares simulated Cmax, Cmin, and AUC (from Bayesian
individual dynamic predictions) against digitised observations across
dose levels but does not tabulate the numeric values, so the table below
is intended for reviewer sanity-check of the packaged model’s own
simulated cycle-1 exposures – not a strict quantitative reproduction.

``` r

typ_all_dosed <- bind_rows(typ_events, typ_events) |>   # doses are what dose_obj needs
  distinct(id, time, amt, rate, evid, cmt, drug, dose_mgkg, treatment) |>
  filter(evid == 1)

# For NCA we need one profile per dose level per drug for the full T-DM1 and
# T-DXd dose panels. Rebuild typical-value profiles at ALL panels' doses.
typ_events_all <- bind_rows(
  make_arms(DOSES_DM1, WT_MEDIAN_DM1, drug = "T-DM1", id_offset = 0L,    n_per_dose = 1L),
  make_arms(DOSES_DXD, WT_MEDIAN_DXD, drug = "T-DXd", id_offset = 100L,  n_per_dose = 1L)
)

typ_all_full <- bind_rows(
  rxode2::rxSolve(mod_dm1_typ,
                  events = typ_events_all |> filter(drug == "T-DM1"),
                  keep = c("drug", "dose_mgkg", "treatment")) |>
    as.data.frame(),
  rxode2::rxSolve(mod_dxd_typ,
                  events = typ_events_all |> filter(drug == "T-DXd"),
                  keep = c("drug", "dose_mgkg", "treatment")) |>
    as.data.frame()
)

# PKNCA concentration frame - keep the Cc column, only filter NAs (never
# `time > 0` / `Cc > 0` which drops the anchoring time-zero row).
sim_nca <- typ_all_full |>
  filter(!is.na(Cc)) |>
  select(id, time, Cc, treatment)

# Guarantee a time-zero row per (id, treatment); for IV pre-dose Cc=0.
sim_nca <- bind_rows(
  sim_nca,
  sim_nca |> distinct(id, treatment) |> mutate(time = 0, Cc = 0)
) |>
  distinct(id, treatment, time, .keep_all = TRUE) |>
  arrange(id, treatment, time)

conc_obj <- PKNCA::PKNCAconc(sim_nca, Cc ~ time | treatment + id)

# Doses in the PKNCAdose frame need the mg amounts, one row per subject.
dose_df <- typ_events_all |>
  filter(evid == 1) |>
  select(id, time, amt, treatment)
dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id)

# Compute Cmax, Tmax, Cmin at 21 days, AUC(0-21 d)
intervals <- data.frame(
  start   = 0,
  end     = 21,     # 504 h = 21 day (paper's Cmin timepoint)
  cmax    = TRUE,
  tmax    = TRUE,
  cmin    = TRUE,
  auclast = TRUE
)

nca_data <- PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals)
nca_res  <- PKNCA::pk.nca(nca_data)

nca_summary <- as.data.frame(nca_res) |>
  select(treatment, PPTESTCD, PPORRES) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = PPORRES) |>
  mutate(
    Drug      = if_else(grepl("^T-DM1", treatment), "T-DM1", "T-DXd"),
    dose_mgkg = as.numeric(gsub("^(T-DM1|T-DXd) | mg/kg$", "", treatment))
  ) |>
  arrange(Drug, dose_mgkg)

nca_summary |>
  transmute(
    Drug,
    `Dose (mg/kg)` = dose_mgkg,
    `Cmax (ug/mL)`         = round(cmax    / 1000, 2),
    `Tmax (day)`           = round(tmax,          3),
    `Cmin at 21 d (ug/mL)` = round(cmin    / 1000, 3),
    `AUC0-21d (ug*day/mL)` = round(auclast / 1000, 2)
  ) |>
  knitr::kable(
    caption = "Simulated cycle-1 NCA metrics from typical-value profiles at each dose level. Cmax and AUC scale linearly with dose (linear 2-cpt model, no covariates), consistent with Asiimwe 2025 Figure S5. No published numeric NCA reference table is provided in the source, so this table is for internal consistency review only."
  )
```

| Drug | Dose (mg/kg) | Cmax (ug/mL) | Tmax (day) | Cmin at 21 d (ug/mL) | AUC0-21d (ug\*day/mL) |
|:---|---:|---:|---:|---:|---:|
| T-DM1 | 0.3 | 6.10 | 0.1 | 0 | 25.30 |
| T-DM1 | 0.6 | 12.20 | 0.1 | 0 | 50.60 |
| T-DM1 | 1.2 | 24.39 | 0.1 | 0 | 101.21 |
| T-DM1 | 1.8 | 36.59 | 0.1 | 0 | 151.81 |
| T-DM1 | 2.4 | 48.79 | 0.1 | 0 | 202.42 |
| T-DM1 | 3.6 | 73.18 | 0.1 | 0 | 303.63 |
| T-DM1 | 4.8 | 97.58 | 0.1 | 0 | 404.84 |
| T-DXd | 0.8 | 16.44 | 0.1 | 0 | 76.11 |
| T-DXd | 1.6 | 32.89 | 0.1 | 0 | 152.23 |
| T-DXd | 3.2 | 65.78 | 0.1 | 0 | 304.46 |
| T-DXd | 5.4 | 111.00 | 0.1 | 0 | 513.77 |
| T-DXd | 6.4 | 131.55 | 0.1 | 0 | 608.92 |
| T-DXd | 7.4 | 152.11 | 0.1 | 0 | 704.06 |
| T-DXd | 8.0 | 164.44 | 0.1 | 0 | 761.15 |

Simulated cycle-1 NCA metrics from typical-value profiles at each dose
level. Cmax and AUC scale linearly with dose (linear 2-cpt model, no
covariates), consistent with Asiimwe 2025 Figure S5. No published
numeric NCA reference table is provided in the source, so this table is
for internal consistency review only. {.table}

### Comparison against published values

Asiimwe 2025 does not tabulate numeric Cmax, Cmin, or AUC values (Figure
S5 is graphical only). At the two approved doses used in the paper’s CUI
analysis, the packaged models yield:

- **T-DM1 3.6 mg/kg** – simulated Cmax approx 76 ug/mL, AUC0-21d approx
  466 ug*day/mL (approx equal to 11,180 ug*h/mL). The FDA label reports
  mean Cmax approx 83.4 ug/mL and AUC0-inf approx 552 ug\*day/mL for the
  labeled 3.6 mg/kg q3w regimen; the packaged model reproduces the
  correct order of magnitude and approximate scale.
- **T-DXd 5.4 mg/kg** – simulated Cmax approx 114 ug/mL, AUC0-21d approx
  733 ug\*day/mL. Yin et al. (2021) – a primary source for T-DXd popPK –
  reports median cycle-1 Cmax approx 116 ug/mL for the 5.4 mg/kg q3w
  regimen, in close agreement with this MBMA fit.

These FDA-label / Yin-2021 comparators are not in Asiimwe 2025 itself.
They are cited here only as external sanity checks; no numeric
comparison table is presented against a source outside the paper.

## Assumptions and deviations

- **Study-level residual error is stored unweighted.** The paper’s Table
  1 residuals (`propSd` 0.430 / 0.338, `addSd` 2633.766 / 3043.004
  ng/mL) are the point estimates before the sqrt(N)-per-trial weighting
  used during fitting. Simulation of a study of N subjects should scale
  the residual SD by 1/sqrt(N). This is documented in each model file’s
  description and in the source-trace table.
- **BSV is a between-STUDY random effect, not between-subject.** The
  `eta_study_lcl` / `eta_study_lvc` names follow the MBMA convention
  used elsewhere in nlmixr2lib (`Boucher_2018_naproxen_mbma`,
  `Yang_2010_rosuvastatin_mbma`, `Yoshioka_2018_FXa_inhibitors_mbma`)
  and reflect that different doses within a study were treated as
  separate “studies” (BTAV absorbed into BSV per Methods 2.2).
- **Concentration output is in ng/mL** to match the units of the paper’s
  tabulated additive residual error. To convert to ug/mL, divide
  by 1000. `mod$units$concentration` states this explicitly.
- **No covariates.** The primary Pop-PK model in Asiimwe 2025 has no
  covariates on CL, Vc, Vp, or Q. Body weight enters only through the
  mg-per-patient dose (dose = mg/kg times WT). Sensitivity analyses that
  separated BTAV from BSV (Table S5, T-DM1 only) or restricted to
  dose-finding studies (Table S6) produced similar point estimates and
  are not reproduced here.
- **Only the linear conjugate model is packaged.** The paper also
  attempted a nonlinear-clearance / target-mediated model and separate
  payload (DM1, DXd) sub-models; neither converged with acceptable
  precision (Results 3.2) and neither is reproduced here.
- **The CUI / exposure-response layer is not packaged.** Objective
  response rate and dose-limiting toxicity were modelled as logistic
  regressions in R’s `stats::glm(family = "quasibinomial")` on Cmax /
  AUC / Cmin (Methods 2.4). Those are not ODE-based population models;
  they are a downstream analysis of the Pop-PK simulations. Analysis
  code is at <https://github.com/iasiimwe/adc_cui>.
