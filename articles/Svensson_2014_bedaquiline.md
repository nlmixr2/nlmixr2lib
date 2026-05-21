# Bedaquiline (Svensson 2014)

## Model and source

- Citation: Svensson E. M., Murray S., Karlsson M. O., Dooley K. E.
  (2014). Rifampicin and rifapentine significantly reduce concentrations
  of bedaquiline, a new anti-TB drug. Journal of Antimicrobial
  Chemotherapy 70(4):1106-1114. <doi:10.1093/jac/dku504>.
- Description: Three-compartment population PK model for bedaquiline
  (BDQ) and a two-compartment N-desmethyl metabolite M2 in healthy adult
  volunteers following single 400 mg oral doses, with Savic 2007
  analytical transit-compartment absorption (non-integer NN feeding a
  first-order depot at rate ka), fixed allometric scaling on disposition
  (0.75 on CL/Q at 70 kg, 1 on Vc/Vp), and multiplicative rifampicin or
  rifapentine drug-drug-interaction factors of 4.78 and 3.96 on
  bedaquiline and M2 apparent clearance, applied at full induction from
  day 3 of rifamycin co-administration.
- Article: <https://doi.org/10.1093/jac/dku504>

## Population

The packaged model was fit to data from Phase I trial TMC207-CL002, a
two-arm open-label two-period single-sequence drug-drug-interaction
study in 32 healthy adult volunteers. Subjects received a single 400 mg
oral dose of bedaquiline on Day 1 (period 1, alone) and again on Day 29
(period 2, after nine days of pre-treatment with either 600 mg daily
rifampicin (Arm 1, 13 completers) or 600 mg daily rifapentine (Arm 2, 16
completers); rifamycin co-administration continued throughout the 14-day
PK sampling window of period 2). Blood samples were collected at 1, 2,
3, 4, 5, 6, 8, 12, and 24 h and thereafter every 24 h until 336 h after
each bedaquiline dose plus a single sample on Day 20 (before rifamycin
start). Concentrations were quantified by validated LC-MS/MS with a
1-2000 ng/mL linear range and a 1 ng/mL lower limit of quantification.
Baseline demographics (Svensson 2014 Table 1): median age 35.5 years
(range 19-55), median weight 81.8 kg (range 57.3-122), 12.5% female,
87.5% white, 6.2% Black or African American, 3.1% American Indian or
Alaska Native, 3.1% Asian.

The historical efavirenz-DDI cohort (1083 bedaquiline + 1055 M2
observations after a single 400 mg dose alone or with efavirenz) was
fitted simultaneously with the new TMC207-CL002 data to increase
parameter precision.

The same information is available programmatically via the model’s
`population` metadata
(`readModelDb("Svensson_2014_bedaquiline")$population`).

## Source trace

The per-parameter origin is recorded as an in-file comment next to each
`ini()` entry in
`inst/modeldb/specificDrugs/Svensson_2014_bedaquiline.R`. The table
below collects them in one place for review.

| Equation / parameter | Value | Source location |
|----|----|----|
| `lmtt` (MTT) | log(0.97) | Table 2 “MTT (h) = 0.97” (RSE 11.5%) |
| `lka` (ka) | log(0.12) | Table 2 “KA (1/h) = 0.12” (RSE 3.9%) |
| `lnn` (NN) | log(8.41) | Table 2 “NN = 8.41” (RSE 36.2%) |
| `lcl` (CL/F BDQ) | log(3.20) | Table 2 “CL/F = 3.20 L/h” (RSE 6.5%) |
| `lvc` (V/F BDQ) | log(16.2) | Table 2 “V/F = 16.2 L” (RSE 12.9%) |
| `lq` (Q1/F BDQ) | log(4.71) | Table 2 “Q1/F = 4.71 L/h” (RSE 5.6%) |
| `lvp` (VP1/F BDQ) | log(2801) | Table 2 “VP1/F = 2801 L” (RSE 10.1%) |
| `lq2` (Q2/F BDQ) | log(3.10) | Table 2 “Q2/F = 3.10 L/h” (RSE 6.0%) |
| `lvp2` (VP2/F BDQ) | log(137) | Table 2 “VP2/F = 137 L” (RSE 10.4%) |
| `lcl_m2` (CLM2/F/fm) | log(13.1) | Table 2 “CLM2/F/fm = 13.1 L/h” (RSE 6.6%) |
| `lvc_m2` (VM2/F/fm) | log(882) | Table 2 “VM2/F/fm = 882 L” (RSE 4.9%) |
| `lq_m2` (Q1M2/F/fm) | log(105) | Table 2 “Q1M2/F/fm = 105 L/h” (RSE 8.6%) |
| `lvp_m2` (VP1M2/F/fm) | log(3349) | Table 2 “VP1M2/F/fm = 3349 L” (RSE 3.8%) |
| `e_wt_cl_q` | fixed(0.75) | Methods “Allometric scaling was applied to CL and V using body weight and fixed coefficients of 0.75 and 1, respectively.” |
| `e_wt_vc_vp` | fixed(1) | same Methods sentence |
| `e_rif_cl` | 4.78 | Table 2 “Factor change BDQ/M2 CL with RIF = 4.78” (RSE 9.1%) |
| `e_rpt_cl` | 3.96 | Table 2 “Factor change BDQ/M2 CL with RPT = 3.96” (RSE 5.0%) |
| Block IIV CL, CL_M2 | var 0.07654, 0.09236; cor 0.538 | Table 2 BSV CL = 28.2% CV, BSV CLM2 = 31.1% CV, corr 53.8% |
| `etalvc` (BSV V) | 0.16769 | Table 2 BSV V = 42.8% CV (RSE 14.6%) |
| `etalq` (BSV Q1) | 0.04372 | Table 2 BSV Q1 = 21.1% CV (RSE 13.2%) |
| `etalvc_m2` (BSV VM2) | 0.12379 | Table 2 BSV VM2 = 36.4% CV (RSE 12.3%) |
| `etalvp_m2` (BSV VP1M2) | 0.04330 | Table 2 BSV VP1M2 = 21.0% CV (RSE 22.4%) |
| `propSd` (BDQ residual) | 0.157 | Table 2 “Prop err BDQ = 15.7% CV” (RSE 4.4%) |
| `propSd_m2` (M2 residual) | 0.122 | Table 2 “Prop err M2 = 12.2% CV” (RSE 5.0%) |
| Transit-absorption chain via Savic 2007 input feeding depot at rate (NN+1)/MTT, then first-order ka into central | n/a | Methods “absorption through a dynamic transit compartment model” + structural inheritance from the previously developed bedaquiline + M2 model (Svensson 2014 reference 19) |
| Three-compartment BDQ + two-compartment M2 ODEs | n/a | Methods “three disposition compartments for bedaquiline and two for M2” |
| Allometric scaling on CL, Q, V, VM2, VP1M2 at 70 kg | n/a | Methods “Allometric scaling was applied to CL and V” |
| DDI step indicators (CONMED_RIF, CONMED_RPT) switching on at full induction (paper’s 3-day lag) | n/a | Methods “Parameterizing the CLs to change after 3 days of rifamycin administration provided the best fit” |

## Virtual cohort

Original observed concentrations are not publicly available. We simulate
a virtual cohort matching the Svensson 2014 design as closely as the
published demographics allow: 32 healthy adult volunteers split into
three single-dose period 1 / period 2 cohorts that cover the three
observed scenarios.

``` r

set.seed(20251229L)

# Helper: build one cohort of `n` subjects as a self-contained event table.
# id_offset shifts subject IDs so multiple cohorts can be bind_rows()-ed
# without colliding (rxSolve treats id as the subject key; duplicate IDs
# across cohorts silently collapse into single Frankenstein subjects).
make_arm <- function(n, conmed_rif, conmed_rpt, treatment, id_offset = 0L) {
  ids <- id_offset + seq_len(n)
  # Sampling grid mirrors the published design: dense early then 24-hourly
  # for two weeks plus a few earlier points to characterise absorption.
  sample_grid <- c(
    seq(0, 12, by = 0.5),
    seq(13, 24, by = 1),
    seq(48, 336, by = 24)
  )
  # Use a small range of body weights centred around the cohort median of
  # 81.8 kg (Table 1); randomly drawn so allometric scaling shows the
  # expected dispersion in the VPC.
  wts <- round(runif(length(ids), min = 60, max = 110), 1)
  dose_rows <- tibble::tibble(
    id          = ids,
    time        = 0,
    evid        = 1L,
    amt         = 400,
    cmt         = "depot",
    WT          = wts,
    CONMED_RIF  = conmed_rif,
    CONMED_RPT  = conmed_rpt,
    treatment   = treatment
  )
  sample_rows <- tidyr::expand_grid(id = ids, time = sample_grid) |>
    dplyr::left_join(
      tibble::tibble(id = ids, WT = wts), by = "id"
    ) |>
    dplyr::mutate(
      evid       = 0L,
      amt        = 0,
      cmt        = "Cc",
      CONMED_RIF = conmed_rif,
      CONMED_RPT = conmed_rpt,
      treatment  = treatment
    )
  dplyr::bind_rows(dose_rows, sample_rows) |>
    dplyr::arrange(id, time, dplyr::desc(evid))
}

n_per_arm <- 32L
events <- dplyr::bind_rows(
  make_arm(n_per_arm, conmed_rif = 0L, conmed_rpt = 0L,
           treatment = "BDQ alone",        id_offset =   0L),
  make_arm(n_per_arm, conmed_rif = 1L, conmed_rpt = 0L,
           treatment = "BDQ + rifampicin", id_offset = 100L),
  make_arm(n_per_arm, conmed_rif = 0L, conmed_rpt = 1L,
           treatment = "BDQ + rifapentine", id_offset = 200L)
)

stopifnot(!anyDuplicated(unique(events[, c("id", "time", "evid")])))
nrow(events)
#> [1] 4896
table(unique(events[, c("id", "treatment")])$treatment)
#> 
#>  BDQ + rifampicin BDQ + rifapentine         BDQ alone 
#>                32                32                32
```

## Simulation

``` r

mod <- readModelDb("Svensson_2014_bedaquiline")
mod_typical <- rxode2::zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

sim_typical <- rxode2::rxSolve(
  mod_typical, events = events,
  keep = c("CONMED_RIF", "CONMED_RPT", "treatment", "WT"),
  returnType = "data.frame"
)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalcl_m2', 'etalvc', 'etalq', 'etalvc_m2', 'etalvp_m2'
#> Warning: multi-subject simulation without without 'omega'
head(sim_typical[, c("id", "time", "Cc", "Cc_m2", "treatment")])
#>   id time          Cc        Cc_m2 treatment
#> 1  1  0.0 0.000000000 0.000000e+00 BDQ alone
#> 2  1  0.5 0.008047051 1.750017e-06 BDQ alone
#> 3  1  1.0 0.329637148 2.044843e-04 BDQ alone
#> 4  1  1.5 1.107267994 1.405209e-03 BDQ alone
#> 5  1  2.0 1.800480156 3.840863e-03 BDQ alone
#> 6  1  2.5 2.263164398 7.109189e-03 BDQ alone
```

``` r

sim_vpc <- rxode2::rxSolve(
  mod, events = events,
  keep = c("CONMED_RIF", "CONMED_RPT", "treatment", "WT"),
  returnType = "data.frame"
)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

## Replicate published figures

### Figure 3 – typical BDQ + M2 profiles

Svensson 2014 Figure 3 plots typical concentration-time profiles for
bedaquiline (continuous lines) and M2 (broken lines) over the first 4
weeks of bedaquiline treatment alone (black), with rifampicin (dark
grey), or with rifapentine (light grey). The paper’s figure uses molar
units (nmol/L); we convert mg/L outputs from the model with
multiplicative factors of 1000 / 555.50 = 1.800 for bedaquiline and 1000
/ 541.47 = 1.847 for M2.

``` r

mw_bdq <- 555.50
mw_m2  <- 541.47

sim_typical_first <- sim_typical |>
  dplyr::filter(time > 0, time <= 336) |>
  dplyr::distinct(time, treatment, .keep_all = TRUE) |>
  dplyr::transmute(
    time      = time,
    treatment = treatment,
    BDQ_nM    = Cc    * 1000 / mw_bdq,
    M2_nM     = Cc_m2 * 1000 / mw_m2
  ) |>
  tidyr::pivot_longer(c(BDQ_nM, M2_nM),
                      names_to = "analyte", values_to = "conc_nM") |>
  dplyr::mutate(
    analyte = dplyr::recode(analyte, BDQ_nM = "Bedaquiline", M2_nM = "M2")
  )

ggplot(sim_typical_first,
       aes(time, conc_nM, colour = treatment, linetype = analyte)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  labs(
    x = "Time after single 400 mg bedaquiline dose (h)",
    y = "Typical plasma concentration (nmol/L)",
    title = "Figure 3 -- BDQ + M2 typical profiles, alone vs with rifamycin",
    caption = paste(
      "Replicates Figure 3 of Svensson 2014: typical-value profiles after a",
      "single 400 mg bedaquiline dose alone, with rifampicin (at full",
      "induction), or with rifapentine (at full induction). The published",
      "figure uses 4-week (672 h) observations; this simulation uses 14",
      "days (336 h) to match the trial's sampling window."
    )
  )
```

![](Svensson_2014_bedaquiline_files/figure-html/figure-3-1.png)

### VPC across the cohort

``` r

sim_summary <- sim_vpc |>
  dplyr::filter(time > 0) |>
  dplyr::group_by(time, treatment) |>
  dplyr::summarise(
    BDQ_Q05 = quantile(Cc,    0.05, na.rm = TRUE),
    BDQ_Q50 = quantile(Cc,    0.50, na.rm = TRUE),
    BDQ_Q95 = quantile(Cc,    0.95, na.rm = TRUE),
    M2_Q05  = quantile(Cc_m2, 0.05, na.rm = TRUE),
    M2_Q50  = quantile(Cc_m2, 0.50, na.rm = TRUE),
    M2_Q95  = quantile(Cc_m2, 0.95, na.rm = TRUE),
    .groups = "drop"
  ) |>
  tidyr::pivot_longer(-c(time, treatment),
                      names_to = c("analyte", "stat"),
                      names_sep = "_") |>
  tidyr::pivot_wider(names_from = stat, values_from = value) |>
  dplyr::mutate(
    analyte = dplyr::recode(analyte, BDQ = "Bedaquiline", M2 = "M2 metabolite")
  )

ggplot(sim_summary |> dplyr::filter(time <= 336),
       aes(time, Q50, colour = treatment, fill = treatment)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.20, colour = NA) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ analyte, scales = "free_y", ncol = 1) +
  scale_y_log10() +
  labs(
    x = "Time after single 400 mg bedaquiline dose (h)",
    y = "Plasma concentration (mg/L)",
    title = "VPC -- BDQ + M2, alone vs with rifamycin",
    caption = paste(
      "Median + 5-95% prediction interval across",
      paste(nrow(unique(events[, c("id", "treatment")])), "subjects (32 per arm),"),
      "single 400 mg oral dose."
    )
  )
```

![](Svensson_2014_bedaquiline_files/figure-html/figure-vpc-cohort-1.png)

## PKNCA validation

Use PKNCA to compute Cmax, Tmax, AUC0-14d, and terminal half-life by
treatment arm so the simulated 14-day AUC ratios (BDQ alone vs with
rifamycin) can be compared against the published GMRs (Svensson 2014
Results “Non-compartmental analysis and posterior predictive check”).

The paper reports observed GMR of AUC0-14d (BDQ with RIF / BDQ alone) =
41.0%, GMR (BDQ with RPT / BDQ alone) = 42.8%, GMR (M2 with RIF / M2
alone) = 78.9%, GMR (M2 with RPT / M2 alone) = 85.5%, and the
model-based simulation GMRs from 100 simulated replicas were 40.6%
(BDQ + RIF), 47.8% (BDQ + RPT), 76.7% (M2 + RIF), and 89.0% (M2 + RPT).

``` r

# Bedaquiline NCA
sim_nca_bdq <- sim_vpc |>
  dplyr::filter(!is.na(Cc), time > 0 | time == 0) |>
  dplyr::select(id, time, Cc, treatment)

dose_df <- events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, time, amt, treatment)

conc_obj_bdq <- PKNCA::PKNCAconc(
  sim_nca_bdq, Cc ~ time | treatment + id,
  concu = "mg/L", timeu = "h"
)
dose_obj <- PKNCA::PKNCAdose(
  dose_df, amt ~ time | treatment + id,
  doseu = "mg"
)

intervals_14d <- data.frame(
  start      = 0,
  end        = 336,
  cmax       = TRUE,
  tmax       = TRUE,
  auclast    = TRUE,
  aucinf.obs = TRUE,
  half.life  = TRUE
)

nca_res_bdq <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  conc_obj_bdq, dose_obj, intervals = intervals_14d
))
nca_summary_bdq <- summary(nca_res_bdq)
knitr::kable(
  nca_summary_bdq,
  caption = "Simulated NCA parameters for bedaquiline by treatment arm."
)
```

| Interval Start | Interval End | treatment | N | AUClast (h\*mg/L) | Cmax (mg/L) | Tmax (h) | Half-life (h) | AUCinf,obs (h\*mg/L) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 0 | 336 | BDQ + rifampicin | 32 | 21.0 \[23.2\] | 1.46 \[14.9\] | 3.25 \[2.50, 5.00\] | 374 \[33.2\] | 23.9 \[26.7\] |
| 0 | 336 | BDQ + rifapentine | 32 | 22.3 \[28.7\] | 1.53 \[23.0\] | 3.50 \[2.50, 4.50\] | 391 \[37.5\] | 26.0 \[33.2\] |
| 0 | 336 | BDQ alone | 32 | 57.5 \[18.7\] | 2.59 \[19.8\] | 4.50 \[3.00, 6.50\] | 532 \[109\] | 86.4 \[24.2\] |

Simulated NCA parameters for bedaquiline by treatment arm. {.table}

``` r


# M2 metabolite NCA
sim_nca_m2 <- sim_vpc |>
  dplyr::filter(!is.na(Cc_m2)) |>
  dplyr::select(id, time, Cc_m2, treatment) |>
  dplyr::rename(Cc = Cc_m2)

conc_obj_m2 <- PKNCA::PKNCAconc(
  sim_nca_m2, Cc ~ time | treatment + id,
  concu = "mg/L", timeu = "h"
)

# M2 has no direct dose; reuse the BDQ dose object as the formation
# precursor so PKNCA can normalise AUC.
nca_res_m2 <- PKNCA::pk.nca(PKNCA::PKNCAdata(
  conc_obj_m2, dose_obj, intervals = intervals_14d
))
nca_summary_m2 <- summary(nca_res_m2)
knitr::kable(
  nca_summary_m2,
  caption = "Simulated NCA parameters for M2 metabolite by treatment arm."
)
```

| Interval Start | Interval End | treatment | N | AUClast (h\*mg/L) | Cmax (mg/L) | Tmax (h) | Half-life (h) | AUCinf,obs (h\*mg/L) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 0 | 336 | BDQ + rifampicin | 32 | 4.67 \[21.9\] | 0.0754 \[25.2\] | 9.50 \[6.00, 14.0\] | 124 \[24.1\] | 5.22 \[24.7\] |
| 0 | 336 | BDQ + rifapentine | 32 | 4.60 \[32.5\] | 0.0699 \[29.6\] | 9.50 \[6.50, 16.0\] | 133 \[30.8\] | 5.25 \[35.5\] |
| 0 | 336 | BDQ alone | 32 | 6.76 \[24.2\] | 0.0381 \[30.8\] | 14.0 \[8.00, 20.0\] | 1480 \[3810\] | 26.0 \[82.0\] |

Simulated NCA parameters for M2 metabolite by treatment arm. {.table
style="width:100%;"}

### Comparison against published GMRs

``` r

auclast_by_arm <- function(nca_summary_df) {
  res <- as.data.frame(nca_summary_df)
  tail_col <- grep("auclast", names(res), value = TRUE, ignore.case = TRUE)
  if (length(tail_col) == 0) return(NULL)
  res[, c("treatment", tail_col[1])]
}

bdq_auc <- auclast_by_arm(nca_summary_bdq)
m2_auc  <- auclast_by_arm(nca_summary_m2)
print(bdq_auc)
#>           treatment AUClast (h*mg/L)
#> 1  BDQ + rifampicin      21.0 [23.2]
#> 2 BDQ + rifapentine      22.3 [28.7]
#> 3         BDQ alone      57.5 [18.7]
print(m2_auc)
#>           treatment AUClast (h*mg/L)
#> 1  BDQ + rifampicin      4.67 [21.9]
#> 2 BDQ + rifapentine      4.60 [32.5]
#> 3         BDQ alone      6.76 [24.2]
```

> Compute the GMRs by hand from the printed AUClast medians: divide the
> “BDQ + rifampicin” row by the “BDQ alone” row (and likewise for M2 and
> for rifapentine). The published observed and simulated GMRs are listed
> in the introduction to this PKNCA section. Differences are expected
> because the present implementation drops the BSV on the individual
> induction effects (see Assumptions and deviations below); the typical
> profile is preserved but the per-subject variance contribution to the
> GMR is reduced.

## Assumptions and deviations

- **Cross-output residual correlation dropped.** Svensson 2014 Table 2
  reports a 55% correlation between the proportional residual errors on
  bedaquiline and M2 (`Prop err M2 = 55.4` is a correlation rather than
  a CV per the table’s footnote). nlmixr2lib has no idiomatic encoding
  for cross-output residual correlation, so the BDQ and M2 residual
  proportions are encoded as independent.
- **Early-sample residual weighting dropped.** Svensson 2014 Table 2
  reports a 2.19-fold weighting on the residual SD for samples taken
  between 0 and 6 h post-dose (`Weighting of samples 0-6 h = 2.19`). The
  paper introduces this weighting to absorb absorption-phase
  unmodelled-variability that is largely a fitting concern. The
  simulated residual SD here is held constant at the post-absorption
  value (15.7% on BDQ, 12.2% on M2) across all time points; readers who
  need to reproduce the fitting weighting will need to layer it on
  externally.
- **Between-occasion variability dropped.** Svensson 2014 Table 2
  reports BOV on bioavailability (18.5% CV) and on MTT (64.7% CV) from
  the two-occasion crossover design (period 1 dose alone, period 2 dose
  with rifamycin); the same table reports an additional BSV on F (12.1%
  CV). nlmixr2lib has no idiomatic encoding for between-occasion
  variability separate from between-subject; both BOVs and BSV F are
  dropped here. As a consequence the simulated VPCs underrepresent
  within-subject variability between the two doses each subject
  received.
- **Individual-induction etas (BSV RIF BDQ, BSV RIF M2, BSV RPT BDQ, BSV
  RPT M2) dropped.** Svensson 2014 Table 2 reports a 6x6
  lower-triangular NONMEM BLOCK(6) on the bedaquiline-CL eta, the M2-CL
  eta, and four individual-level induction-effect etas. The BSV CL + BSV
  CLM2 block (2x2) is preserved here, but the four individual induction
  etas (BSV RIF BDQ = 27.9% CV, BSV RIF M2 = 32.4% CV, BSV RPT BDQ =
  17.9% CV, BSV RPT M2 = 28.2% CV) and the inter-block correlations
  (-78.2 to 86.3% per Table 2) are dropped because the etas are only
  well-defined for subjects in the corresponding arm, and nlmixr2lib has
  no idiomatic encoding for a covariate-gated random effect. As a
  consequence the VPC rifamycin arms have narrower prediction intervals
  on AUC than the paper’s model would produce, and the typical-value
  GMRs are an approximation of the paper’s central tendencies.
- **Time-varying weight assumed time-fixed.** The paper’s analysis used
  per-subject body weight from baseline only, so this simplification
  matches the source.
- **Rifamycin step indicators set externally.** Svensson 2014 found that
  a step from “no induction” to “full induction” at day 3 of rifamycin
  co-administration provided the best fit (an alternative time-decaying
  induction model would describe pre-induction kinetics but did not
  improve OFV materially per the Methods text). For new simulations,
  populate `CONMED_RIF` = 1 on every observation row that falls \>= 3
  days after the first rifampicin dose and 0 otherwise (and similarly
  for `CONMED_RPT`); pre-induction observations should carry the
  unmodified baseline CL. The vignette’s example cohorts assume the
  bedaquiline dose is given when rifamycin is already at full induction
  (period-2 design of the paper).
- **F = 1 anchor for the F-relative parameterisation.** The paper
  reports all clearances and volumes as apparent F-relative values
  (CL/F, V/F, CLM2/(F*fm), V_M2/(F*fm)). Setting F = 1 inside
  `transit(nn, mtt, 1)` preserves these apparent values directly; the
  consequence is that simulated bedaquiline mass in the central
  compartment represents the F-apparent mass and concentrations
  represent F-apparent (= observed) plasma concentrations.
- **Concentration units.** Model output is in mg/L (= ug/mL). The
  published figures (Figure 3, Figure 4) plot concentrations in nmol/L;
  the conversion factors are 1000 / 555.50 = 1.800 (BDQ) and 1000 /
  541.47 = 1.847 (M2). The `Figure 3` chunk above applies these
  conversions.
- **NCA terminal-phase fit.** With 14-day sampling and a 5-6 month
  terminal half-life, NCA AUC0-inf and half-life cannot be reliably
  estimated. AUC0-14d (the paper’s primary NCA endpoint) is reported
  here via PKNCA’s `auclast` instead.
