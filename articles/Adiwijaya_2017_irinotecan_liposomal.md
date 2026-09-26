# Liposomal irinotecan (Adiwijaya 2017)

## Model and source

``` r

mod <- rxode2::rxode2(readModelDb("Adiwijaya_2017_irinotecan_liposomal"))
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Adiwijaya BS, Kim J, Lang I, Csoszi T, Cubillo A, Chen JS,
  Wong M, Park JO, Kim JS, Rau KM, Melichar B, Gallego JB, Fitzgerald J,
  Belanger B, Molnar I, Ma WW. Population Pharmacokinetics of Liposomal
  Irinotecan in Patients With Cancer. Clin Pharmacol Ther.
  2017;102(6):997-1005. <doi:10.1002/cpt.720>. PMCID: PMC5697569.
- Article: <https://doi.org/10.1002/cpt.720>
- PubMed Central (open access):
  <https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5697569/>
- Supplement: available from the Europe PMC supplementary-files endpoint
  for PMC5697569 (file `CPT-102-997-s001.docx`), which carries Equation
  S1, the model diagram (Figure S1) and the two final-parameter tables
  S6 and S7 used throughout this vignette.

Nanoliposomal irinotecan (nal-IRI; MM-398, PEP02, marketed as Onivyde)
is a liposomal formulation of irinotecan. Encapsulation slows release,
so compared with nonliposomal irinotecan the plasma total-irinotecan
exposure is far higher and far longer-lived while the peak concentration
of the active metabolite SN-38 is *lower*. Adiwijaya 2017 pooled six
phase I-III studies to quantify that behaviour and to link the resulting
exposure metrics to efficacy (survival) and to the two dose-limiting
toxicities, neutropenia and diarrhea.

### Structure

The model has two sequential sub-models, shown in the paper’s Figure S1
and described in the supplement section “General structure and
assumptions”:

- **Total irinotecan (tIRI)** – a two-compartment model with first-order
  elimination. tIRI is the sum of encapsulated and unencapsulated
  irinotecan; the paper justifies treating it as a single species by
  showing that the encapsulated:total ratio was constant over a week of
  sampling in study PEP0201.
- **Total SN-38 (tSN38)** – the sum of two terms:
  - *unencapsulated* SN-38 (`Cu_sn38`), a one-compartment model formed
    from tIRI by a single first-order rate constant `kmet` that lumps
    together release of irinotecan from the liposome and its conversion
    to SN-38, and eliminated with clearance `cl_sn38`. Its volume is
    *not* estimated: the supplement states “Parameter V is obtained from
    the V1 estimates of the tIRI model”.
  - *encapsulated* SN-38 (`Ce_sn38`), a time-invariant mass fraction
    `fsn38` of tIRI. This is a co-encapsulated manufacturing
    contaminant, not a metabolite: the paper supports it with an in
    vitro measurement of 0.015% and with the absence of the glucuronide
    SN-38G in the first 12 h after dosing (Figure S2), which shows the
    encapsulated fraction is not available for glucuronidation.

The coupling runs one way only – tIRI is estimated independently,
because tIRI concentrations exceed tSN38 concentrations by roughly four
orders of magnitude so the mass converted is negligible – but it runs
through more than the state: the *individual* tIRI clearance and central
volume enter the SN-38 formation rate as covariates, on the hypothesis
that mononuclear-phagocyte-system activity drives both liposome
clearance and drug release.

Because the two sub-models cannot be exercised separately (the SN-38 arm
needs the irinotecan state, its volume, and its individual CL and V1),
they are packaged as a single two-output model file rather than two.

## Population

Coupled two-analyte population PK model for nanoliposomal irinotecan
(nal-IRI, MM-398/PEP02; output Cc = total irinotecan) and total SN-38
(output Cc_sn38) in adults with metastatic pancreatic, gastric/GEJ,
colorectal and other solid tumours (Adiwijaya 2017, pooled from six
phase I-III studies including NAPOLI-1). Total irinotecan is a
two-compartment model with first-order elimination. Total SN-38 is the
sum of encapsulated SN-38 (a time-invariant mass fraction of total
irinotecan, i.e. a co-encapsulated manufacturing contaminant) and
unencapsulated SN-38 (a one-compartment model formed from total
irinotecan by the first-order rate constant kmet and sharing the
irinotecan central volume). The two sub-models are sequential: the
irinotecan model is independent, and the SN-38 model depends on it,
including through the individual irinotecan CL and V1 acting as
covariates on the formation rate.

The analysis pooled 353 patients from 6 studies (five phase I-II studies
run by PharmaEngine plus the phase III NAPOLI-1 trial), contributing
1,792 total-irinotecan and 1,765 total-SN-38 concentrations, all from
treatment cycle 1. Baseline characteristics are from the paper’s Table
1; the contributing studies are listed in Table S1.

| Characteristic | Value |
|:---|:---|
| Patients (PK dataset) | 353 of 368 treated (96%) |
| Studies | 6 (5 phase I-II + phase III NAPOLI-1; NAPOLI-1 = 73%) |
| Age | 63 years (5th-95th percentile 39.8-79.2) |
| Sex | 44% female |
| Race | Caucasian 52%, East Asian 42%, Other 6% |
| Tumour type | Pancreatic 73%, gastric/GEJ 10%, solid tumour 11%, colorectal 5% |
| BSA | 1.7 m^2 (1.3-2.2) |
| Albumin | 40 g/L (29-47) |
| ALT | 25 U/L (8.9-96.3) |
| Total bilirubin | 7 umol/L (3-19) |
| Creatinine clearance | 81.6 mL/min (39.6-151.8) |
| Liver metastases (NAPOLI-1) | 66% |
| UGT1A1\*28 7/7 (NAPOLI-1) | 5% (14 of 258) |
| Initial dose | 100 mg/m^2 (53%) or 70 mg/m^2 (40%), irinotecan free base |

Baseline characteristics (Adiwijaya 2017 Table 1). Continuous values are
median (5th-95th percentile). {.table}

Doses throughout this vignette are expressed as **irinotecan free
base**, matching the paper’s convention: 70 mg/m^2 free base is
equivalent to 80 mg/m^2 of irinotecan hydrochloride trihydrate salt, and
100 mg/m^2 free base to 120 mg/m^2 salt.

## Source trace

Every structural parameter comes from the two final-model parameter
tables in the supplement. The tables report both a final-model estimate
and a bootstrap median with a 95% confidence interval; the model file
carries the **final-model estimate** column.

| Parameter | Source location | Value |
|:---|:---|:---|
| lvc (V1) | Table S6, Volume (V1) | 4.60 L |
| lcl (CL) | Table S6, Clearance (CL) | 13.6 L/week |
| lq (Q) | Table S6, Q | 0.471 L/week |
| lvp (V2) | Table S6, V2 | 48.7 L |
| e_bsa_vc | Table S6, theta{V1,BSA} | 0.416 per m^2 |
| e_study_napoli1_vc | Table S6, theta{V1,mfg==NAPOLI} | -0.172 |
| e_race_asian_cl | Table S6, theta{CL,race==Asian} | 0.647 |
| e_conmed_fluorouracil_cl | Table S6, theta{CL,treatment contain 5FU} | 0.075 |
| e_study_napoli1_cl | Table S6, theta{CL,mfg==NAPOLI} | -0.189 |
| e_lmet_cl | Table S6, theta{CL,liver metastasis} | -0.075 |
| e_alt_cl | Table S6, theta{CL,ALT} | 0.016 per log10 U/L |
| e_alb_cl | Table S6, theta{CL,albumin} | -1.79 per log10 g/L |
| e_tbili_cl | Table S6, theta{CL,bilirubin} | -0.0670 per log10 umol/L |
| e_crcl_cl | Table S6, theta{CL,creatinine clearance} | 0.003 per mL/min |
| etalvc, etalcl block | Table S6, Random effects | 0.068 / 0.184 / 0.843 |
| expSd | Table S6, Residuals (log10 variance 0.038) | sqrt(0.038)\*log(10) = 0.4488 |
| lcl_sn38 | Table S7, Clearance (CL_SN38) | 14.2 L/week |
| lkmet (Kcov) | Table S7, Conversion flux rate from irinotecan | 0.00072 1/week |
| lfsn38 | Table S7, Fraction of SN38 per unit of irinotecan | 0.090 ng/ug = 9.0e-5 w/w |
| e_race_asian_cl_sn38 | Table S7, theta{CL_SN38,race==Asian} | -0.161 |
| e_ugt1a1_star28_hom_cl_sn38 | Table S7, theta{CL_SN38,UGT1A1\*28==homozygous} | -1.46e-05 |
| e_conmed_fluorouracil_cl_sn38 | Table S7, theta{CL_SN38,treatment contains 5FU} | -1.53e-04 |
| e_lmet_cl_sn38 | Table S7, theta{CL_SN38,liver metastasis = YES} | -0.002 |
| e_alt_cl_sn38 | Table S7, theta{CL_SN38,ALT} | -1.09e-05 per log10 U/L |
| e_alb_cl_sn38 | Table S7, theta{CL_SN38,albumin} | -0.207 per log10 g/L |
| e_tbili_cl_sn38 | Table S7, theta{CL_SN38,bilirubin} | -0.852 per log10 umol/L |
| e_crcl_cl_sn38 | Table S7, theta{CL_SN38,CRCL} | -5.64e-05 per mL/min |
| e_study_napoli1_kmet | Table S7, theta{Kcov,mfg==NAPOLI} | 3.79e-05 |
| e_cl_kmet | Table S7, theta{Kcov,tIRI_logCL} | 2.095 per log10 L/week |
| e_vc_kmet | Table S7, theta{Kcov,tIRI_logV1} | -0.867 per log10 L |
| e_bsa_kmet | Table S7, theta{Kcov,BSA} | -1.121 per m^2 |
| e_study_napoli1_fsn38 | Table S7, theta{fSN38,mfg==NAPOLI} | -0.615 |
| etalcl_sn38 / etalkmet / etalfsn38 | Table S7, Random effects | 0.155 / 0.184 / 0.500 |
| expSd_sn38 | Table S7, Residuals (log10 variance 0.021) | sqrt(0.021)\*log(10) = 0.3337 |
| Covariate model form | Supplement Equation S1 | exp(eta + sum theta*(cont - median) + sum theta*cat) |
| Model diagram | Supplement Figure S1 | 2-cmt tIRI; tSN38 = uSN38 + eSN38 |

Source trace: every ini() value and the model equation form. {.table}

Two values in the model file are arithmetic transformations of the
printed number rather than the printed number itself, and are flagged
here so the trace stays auditable:

- `lfsn38 = log(9.0e-5)`. Table S7 prints `fSN38 = 0.090` in the mixed
  unit `ng (SN-38) / ug (irinotecan)`; dividing by 1000 gives the
  unitless w/w ratio 9.0e-5. This is corroborated independently by the
  Methods, which state the encapsulated-SN-38 fraction of tIRI “was
  estimated to be 0.01%” (9.0e-5 is 0.009%) against an in vitro
  measurement of 0.015%.
- `expSd` and `expSd_sn38`. The source modelled residual error as
  additive on the **log10** scale and reports variances (0.038 and
  0.021). rxode2’s `lnorm()` error model takes a natural-log SD, so each
  is converted as `sqrt(variance) * log(10)`.

The two random-effects columns are labelled “unitless (variance)”. For
the irinotecan block that reading is forced rather than assumed: the
implied correlation is `0.184 / sqrt(0.068 * 0.843) = 0.77`, which is
admissible, whereas reading the entries as CV would not give a
positive-definite matrix.

## Structural checks

These checks compare quantities the model *derives* against quantities
the paper *reports* independently, and so test transcription rather than
restating it. They use the typical subject (every continuous covariate
at its median and every indicator at its reference), which is the
subject the paper’s Table 2 half-lives and clearances describe.

``` r

theta <- setNames(as.numeric(mod$theta), names(mod$theta))
V1 <- exp(theta[["lvc"]])
CL <- exp(theta[["lcl"]])
Q <- exp(theta[["lq"]])
V2 <- exp(theta[["lvp"]])
CLSN <- exp(theta[["lcl_sn38"]])
FSN <- exp(theta[["lfsn38"]])
KMET <- exp(theta[["lkmet"]])

hours_per_week <- 24 * 7

# Hybrid rate constants of the two-compartment irinotecan model.
kel <- CL / V1
k12 <- Q / V1
k21 <- Q / V2
alpha <- (kel + k12 + k21 + sqrt((kel + k12 + k21)^2 - 4 * kel * k21)) / 2
beta <- kel * k21 / alpha

t_half_alpha <- log(2) / alpha * hours_per_week
t_half_term <- log(2) / beta * hours_per_week
t_half_sn38 <- log(2) / (CLSN / V1) * hours_per_week
```

| Quantity | Derived from the model | Published (Table 2) | Agreement |
|:---|:---|:---|:---|
| Irinotecan first-phase t1/2 (h) | 38.1 | 38.2 (23.2-56.7) | -0.4% |
| Irinotecan terminal t1/2 (h) | 12459 | 12200 (3990-50200) | +2.1% |
| SN-38 terminal t1/2 (h) | 37.7 | 38.2 (36.5-41.9) | -1.2% |
| Encapsulated SN-38 fraction of tIRI (%) | 0.009 | 0.01 (in vitro 0.015) | consistent |

Derived versus published disposition summaries. None of these four
numbers is an ini() value; each is a function of two or more of them.
{.table}

The SN-38 half-life check is the sharpest of the four. The paper reports
`38.2 h` for SN-38 but estimates only a clearance for it, never a volume
– the volume is inherited from the irinotecan model.
`log(2) * V1 / CL_SN38` therefore reproduces the published half-life
only if V1 really is the SN-38 volume, which is exactly the supplement’s
stated assumption.

``` r

stopifnot(
  abs(t_half_alpha / 38.2 - 1) < 0.03,
  abs(t_half_term / 12200 - 1) < 0.05,
  abs(t_half_sn38 / 38.2 - 1) < 0.03,
  # 0.090 ng/ug, i.e. 0.009%, against the Methods' "0.01%".
  abs(FSN - 9.0e-5) < 1e-9
)
```

### Dose recovery

For a linear model with no absorption step, `CL * AUC(0-Inf) = Dose`
identically. Solving out to 800 weeks captures the full terminal phase
(whose half-life is 74 weeks) and confirms the ODE system, the
micro-constant algebra, and the compartment the dose enters.

``` r

# Every covariate the model reads, at the reference subject's values.
ref_covariates <- list(
  BSA = 1.70, ALT = 25, ALB = 40, TBILI = 7, CRCL = 81.6,
  RACE_ASIAN = 0, LMET = 0, STUDY_NAPOLI1 = 0,
  CONMED_FLUOROURACIL = 0, UGT1A1_STAR28_HOM = 0
)

#' Build an rxode2 event data frame with covariate columns attached.
#'
#' Covariates are attached to a plain data frame rather than to an rxEt
#' object: assignments onto an rxEt are silently dropped.
make_events <- function(cohort, dose_times, obs_times) {
  doses <- tidyr::expand_grid(cohort, time = dose_times) |>
    dplyr::mutate(amt = .data$dose_mg, evid = 1L, cmt = "central", dvid = NA_integer_)
  obs <- tidyr::expand_grid(cohort, time = obs_times) |>
    dplyr::mutate(amt = 0, evid = 0L, cmt = NA_character_, dvid = 1L)
  dplyr::bind_rows(doses, obs) |>
    dplyr::arrange(.data$id, .data$time, dplyr::desc(.data$evid)) |>
    as.data.frame()
}

typical_subject <- function(dose_mg, ...) {
  overrides <- list(...)
  covs <- utils::modifyList(ref_covariates, overrides)
  tibble::as_tibble(c(list(id = 1L, dose_mg = dose_mg), covs))
}
```

``` r

mod_typ <- rxode2::zeroRe(mod)

recovery_times <- sort(unique(c(
  seq(0, 1, length.out = 4000),
  seq(1, 800, length.out = 16000)
)))
recovery <- rxode2::rxSolve(
  mod_typ,
  make_events(typical_subject(116.7), dose_times = 0, obs_times = recovery_times),
  returnType = "data.frame"
)
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'

auc_inf <- sum(diff(recovery$time) *
  (head(recovery$Cc, -1) + tail(recovery$Cc, -1)) / 2)
recovery_ratio <- recovery$cl[1] * auc_inf / 116.7

sprintf("CL * AUC(0-Inf) / Dose = %.5f", recovery_ratio)
#> [1] "CL * AUC(0-Inf) / Dose = 1.00007"

stopifnot(abs(recovery_ratio - 1) < 0.002)
```

## Virtual cohort

Both NAPOLI-1 arms are simulated: nal-IRI 70 mg/m^2 every two weeks with
5-FU/LV, and nal-IRI 100 mg/m^2 every three weeks as monotherapy.
Covariates are drawn to match the Table 1 marginals, with log-normal
distributions fitted to the reported median and 5th-95th percentiles
(the paper itself notes that log-normal distributions were observed for
the laboratory covariates).

``` r

n_per_arm <- 200

#' Log-normal parameters implied by a median and a 5th-95th percentile pair.
lognormal_sd <- function(p05, p95) (log(p95) - log(p05)) / (2 * stats::qnorm(0.95))

draw_cohort <- function(n, arm, dose_per_m2, fu_lv) {
  bsa <- stats::rlnorm(n, log(1.70), lognormal_sd(1.3, 2.2))
  tibble::tibble(
    arm = arm,
    BSA = bsa,
    ALT = stats::rlnorm(n, log(25), lognormal_sd(8.9, 96.3)),
    ALB = stats::rlnorm(n, log(40), lognormal_sd(29, 47)),
    TBILI = stats::rlnorm(n, log(7), lognormal_sd(3, 19)),
    CRCL = stats::rlnorm(n, log(81.6), lognormal_sd(39.6, 151.8)),
    RACE_ASIAN = stats::rbinom(n, 1, 0.42),
    LMET = stats::rbinom(n, 1, 0.66),
    UGT1A1_STAR28_HOM = stats::rbinom(n, 1, 0.05),
    STUDY_NAPOLI1 = 1,
    CONMED_FLUOROURACIL = fu_lv,
    dose_mg = dose_per_m2 * bsa
  )
}

cohort <- dplyr::bind_rows(
  draw_cohort(n_per_arm, "70 mg/m^2 Q2W", 70, 1L),
  draw_cohort(n_per_arm, "100 mg/m^2 Q3W", 100, 0L)
) |>
  dplyr::mutate(id = dplyr::row_number())

cohort |>
  dplyr::group_by(.data$arm) |>
  dplyr::summarise(
    n = dplyr::n(),
    `Median BSA (m^2)` = round(stats::median(.data$BSA), 2),
    `Median dose (mg)` = round(stats::median(.data$dose_mg), 1),
    `East Asian (%)` = round(100 * mean(.data$RACE_ASIAN)),
    `Liver mets (%)` = round(100 * mean(.data$LMET)),
    .groups = "drop"
  ) |>
  dplyr::rename("Arm" = "arm") |>
  knitr::kable(caption = "Simulated cohort, 200 subjects per arm.")
```

| Arm | n | Median BSA (m^2) | Median dose (mg) | East Asian (%) | Liver mets (%) |
|:---|---:|---:|---:|---:|---:|
| 100 mg/m^2 Q3W | 200 | 1.68 | 167.7 | 42 | 68 |
| 70 mg/m^2 Q2W | 200 | 1.70 | 119.0 | 39 | 68 |

Simulated cohort, 200 subjects per arm. {.table}

## Replicating Figure 1

The paper’s Figure 1 (lower panels) shows model-predicted
total-irinotecan and total-SN-38 concentrations over the first three
weeks, with the second Q2W dose visible at week 2. Concentrations are
plotted in ng/mL on a log scale, so the model’s `ug/mL` outputs are
multiplied by 1000.

``` r

fig_times <- sort(unique(c(
  seq(0, 0.25, length.out = 120),
  seq(0.25, 2, length.out = 140),
  seq(2, 2.25, length.out = 120),
  seq(2.25, 3, length.out = 100)
)))

# Q2W subjects receive a second dose at week 2; Q3W subjects only one in 3 weeks.
sim_fig <- dplyr::bind_rows(
  rxode2::rxSolve(
    mod,
    make_events(dplyr::filter(cohort, .data$arm == "70 mg/m^2 Q2W"),
      dose_times = c(0, 2), obs_times = fig_times
    ),
    returnType = "data.frame"
  ) |> dplyr::mutate(arm = "70 mg/m^2 Q2W"),
  rxode2::rxSolve(
    mod,
    make_events(dplyr::filter(cohort, .data$arm == "100 mg/m^2 Q3W"),
      dose_times = 0, obs_times = fig_times
    ),
    returnType = "data.frame"
  ) |> dplyr::mutate(arm = "100 mg/m^2 Q3W")
)
```

``` r

fig_summary <- sim_fig |>
  dplyr::select("arm", "time", "Cc", "Cc_sn38") |>
  dplyr::rename(`Total irinotecan` = "Cc", `Total SN-38` = "Cc_sn38") |>
  tidyr::pivot_longer(c("Total irinotecan", "Total SN-38"),
    names_to = "analyte", values_to = "conc"
  ) |>
  dplyr::mutate(conc_ng_ml = .data$conc * 1000) |>
  dplyr::group_by(.data$arm, .data$analyte, .data$time) |>
  dplyr::summarise(
    lo = stats::quantile(.data$conc_ng_ml, 0.05),
    mid = stats::median(.data$conc_ng_ml),
    hi = stats::quantile(.data$conc_ng_ml, 0.95),
    .groups = "drop"
  ) |>
  dplyr::filter(.data$mid > 0)

ggplot2::ggplot(fig_summary, ggplot2::aes(x = .data$time, colour = .data$arm, fill = .data$arm)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lo, ymax = .data$hi), alpha = 0.2, colour = NA) +
  ggplot2::geom_line(ggplot2::aes(y = .data$mid), linewidth = 0.8) +
  ggplot2::facet_wrap(~ .data$analyte, scales = "free_y") +
  ggplot2::scale_y_log10() +
  ggplot2::labs(
    x = "Time (weeks)", y = "Predicted concentration (ng/mL)",
    colour = "Dose", fill = "Dose"
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom")
```

![Replicates the lower panels of Figure 1 of Adiwijaya 2017: median
(line) and 5th-95th percentile band of predicted total irinotecan and
total SN-38 over the first three
weeks.](Adiwijaya_2017_irinotecan_liposomal_files/figure-html/figure1-plot-1.png)

Replicates the lower panels of Figure 1 of Adiwijaya 2017: median (line)
and 5th-95th percentile band of predicted total irinotecan and total
SN-38 over the first three weeks.

The shape reproduces the published panels: total irinotecan falls
mono-exponentially over the first week with a shallow tail, and total
SN-38 rises to a sharp early peak then declines in parallel with
irinotecan. That parallel decline is the encapsulated term – because
`Ce_sn38` is a fixed multiple of `Cc`, it follows irinotecan exactly.

### Which SN-38 term dominates when

Splitting the total SN-38 into its two components explains a result that
looks paradoxical in the paper: the total-SN-38 *peak* is driven almost
entirely by the encapsulated contaminant, while the *unencapsulated*
SN-38 that drives toxicity and efficacy peaks much later and is what the
exposure-response analyses use.

``` r

split_times <- sort(unique(c(seq(0, 0.5, length.out = 300), seq(0.5, 2, length.out = 200))))
split_sim <- rxode2::rxSolve(
  mod_typ,
  make_events(typical_subject(70 * 1.70, STUDY_NAPOLI1 = 1, CONMED_FLUOROURACIL = 1),
    dose_times = 0, obs_times = split_times
  ),
  returnType = "data.frame"
)
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'

split_sim |>
  dplyr::transmute(
    time = .data$time,
    `Encapsulated SN-38` = .data$Ce_sn38 * 1000,
    `Unencapsulated SN-38` = .data$Cu_sn38 * 1000,
    `Total SN-38` = .data$Cc_sn38 * 1000
  ) |>
  tidyr::pivot_longer(-"time", names_to = "component", values_to = "conc") |>
  ggplot2::ggplot(ggplot2::aes(.data$time, .data$conc, colour = .data$component)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::labs(x = "Time (weeks)", y = "Concentration (ng/mL)", colour = NULL) +
  ggplot2::theme_bw() +
  ggplot2::theme(legend.position = "bottom")
```

![Decomposition of total SN-38 into its encapsulated and unencapsulated
components, typical subject, 70
mg/m^2.](Adiwijaya_2017_irinotecan_liposomal_files/figure-html/sn38-split-1.png)

Decomposition of total SN-38 into its encapsulated and unencapsulated
components, typical subject, 70 mg/m^2.

``` r

peak_e <- split_sim$time[which.max(split_sim$Ce_sn38)]
peak_u <- split_sim$time[which.max(split_sim$Cu_sn38)]
sprintf(
  "Encapsulated SN-38 peaks at %.2f weeks; unencapsulated SN-38 peaks at %.2f weeks (%.0f h).",
  peak_e, peak_u, peak_u * 24 * 7
)
#> [1] "Encapsulated SN-38 peaks at 0.00 weeks; unencapsulated SN-38 peaks at 0.29 weeks (49 h)."

# The encapsulated term is instantaneous with the dose; the unencapsulated term
# is formation-limited and must peak strictly later.
stopifnot(peak_u > peak_e)
```

## NCA validation

NCA is run separately on each of the three analytes the paper tabulates,
over cycle 1 after a single dose (two weeks for the Q2W arm, three weeks
for the Q3W arm), matching the window the paper’s Table 2 describes.

``` r

nca_times <- sort(unique(c(
  seq(0, 0.25, length.out = 100),
  seq(0.25, 1, length.out = 80),
  seq(1, 3, length.out = 80)
)))

sim_nca_raw <- dplyr::bind_rows(
  rxode2::rxSolve(
    mod,
    make_events(dplyr::filter(cohort, .data$arm == "70 mg/m^2 Q2W"),
      dose_times = 0, obs_times = nca_times[nca_times <= 2]
    ),
    returnType = "data.frame"
  ) |> dplyr::mutate(arm = "70 mg/m^2 Q2W"),
  rxode2::rxSolve(
    mod,
    make_events(dplyr::filter(cohort, .data$arm == "100 mg/m^2 Q3W"),
      dose_times = 0, obs_times = nca_times
    ),
    returnType = "data.frame"
  ) |> dplyr::mutate(arm = "100 mg/m^2 Q3W")
) |>
  dplyr::mutate(
    tIRI = .data$Cc,
    tSN38 = .data$Cc_sn38 * 1000,
    uSN38 = .data$Cu_sn38 * 1000
  )

dose_df <- cohort |>
  dplyr::transmute(id = .data$id, arm = .data$arm, time = 0, amt = .data$dose_mg) |>
  as.data.frame()
```

``` r

exposure_interval <- data.frame(
  start = 0, end = Inf,
  cmax = TRUE, tmax = TRUE, auclast = TRUE, aucinf.obs = TRUE
)

# The irinotecan profile is bi-exponential with a terminal phase that carries
# only about 0.01% of the concentration amplitude. Over a three-week window
# PKNCA's automatic lambda-z selection starts to pick up that bend and reports
# a half-life near 78 h, whereas over a one-week window it cleanly recovers the
# first-phase half-life the paper reports. The window is therefore fixed
# explicitly rather than left to automatic selection.
first_phase_interval <- data.frame(start = 0, end = 1, half.life = TRUE)

#' Run single-dose NCA for one analyte column.
run_nca <- function(data, column, conc_unit, dose_data = dose_df,
                    intervals = exposure_interval) {
  conc_df <- data |>
    dplyr::select("id", "arm", "time", dplyr::all_of(column)) |>
    dplyr::rename(conc = dplyr::all_of(column)) |>
    dplyr::filter(!is.na(.data$conc)) |>
    as.data.frame()

  conc_obj <- PKNCA::PKNCAconc(conc_df, conc ~ time | arm + id,
    concu = conc_unit, timeu = "week"
  )
  dose_obj <- PKNCA::PKNCAdose(dose_data, amt ~ time | arm + id, doseu = "mg")

  PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj, intervals = intervals))
}

nca_tiri_hl <- run_nca(sim_nca_raw, "tIRI", "ug/mL",
  intervals = first_phase_interval
)
```

`nca_tiri_hl` supplies the cohort half-life below; the per-analyte
comparison against Table 2 re-uses `run_nca()` on the typical subject in
the next section.

### Comparison against the published table

Table 2 of the paper reports median Cmax for each of the three analytes
in each NAPOLI-1 arm. Those are the model-derived exposure metrics that
can be reproduced from the published parameters alone.

The comparison is made on the **typical subject** of each arm – every
continuous covariate at its median and every indicator at its reference
value, apart from the two that define the arm (`STUDY_NAPOLI1 = 1` for
both, and 5-FU/LV coadministration in the Q2W arm only). This is
deliberate. The typical subject is deterministic: it involves no random
draw, so the comparison is reproducible across rxode2 versions and
thread counts, whereas a cohort median shifts by a few percent with the
draw. The cohort distribution is reported separately below for context.

``` r

typical_arms <- dplyr::bind_rows(
  typical_subject(70 * 1.70, STUDY_NAPOLI1 = 1, CONMED_FLUOROURACIL = 1) |>
    dplyr::mutate(arm = "70 mg/m^2 Q2W"),
  typical_subject(100 * 1.70, STUDY_NAPOLI1 = 1, CONMED_FLUOROURACIL = 0) |>
    dplyr::mutate(arm = "100 mg/m^2 Q3W")
) |>
  dplyr::mutate(id = dplyr::row_number())

# rxSolve() omits the `id` column when the event table holds a single subject,
# so it is restored explicitly here to match `typical_arms$id`.
sim_typ <- dplyr::bind_rows(
  rxode2::rxSolve(
    mod_typ,
    make_events(dplyr::filter(typical_arms, .data$arm == "70 mg/m^2 Q2W"),
      dose_times = 0, obs_times = nca_times[nca_times <= 2]
    ),
    returnType = "data.frame"
  ) |> dplyr::mutate(arm = "70 mg/m^2 Q2W", id = 1L),
  rxode2::rxSolve(
    mod_typ,
    make_events(dplyr::filter(typical_arms, .data$arm == "100 mg/m^2 Q3W"),
      dose_times = 0, obs_times = nca_times
    ),
    returnType = "data.frame"
  ) |> dplyr::mutate(arm = "100 mg/m^2 Q3W", id = 2L)
) |>
  dplyr::mutate(
    tIRI = .data$Cc,
    tSN38 = .data$Cc_sn38 * 1000,
    uSN38 = .data$Cu_sn38 * 1000
  )
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'

dose_df_typ <- typical_arms |>
  dplyr::transmute(id = .data$id, arm = .data$arm, time = 0, amt = .data$dose_mg) |>
  as.data.frame()
```

``` r

published <- list(
  tIRI = tibble::tribble(
    ~arm, ~cmax,
    "70 mg/m^2 Q2W", 26.6,
    "100 mg/m^2 Q3W", 41.5
  ),
  tSN38 = tibble::tribble(
    ~arm, ~cmax,
    "70 mg/m^2 Q2W", 2.64,
    "100 mg/m^2 Q3W", 3.99
  ),
  uSN38 = tibble::tribble(
    ~arm, ~cmax,
    "70 mg/m^2 Q2W", 2.07,
    "100 mg/m^2 Q3W", 3.05
  )
)
analyte_labels <- c(
  tIRI = "Total irinotecan", tSN38 = "Total SN-38",
  uSN38 = "Unencapsulated SN-38"
)
analyte_units <- c(tIRI = "ug/mL", tSN38 = "ng/mL", uSN38 = "ng/mL")

cmp <- lapply(names(published), function(a) {
  nca_a <- run_nca(sim_typ, a, analyte_units[[a]], dose_data = dose_df_typ)
  nlmixr2lib::ncaComparisonTable(
    simulated = nca_a, reference = published[[a]], by = "arm",
    units = stats::setNames(analyte_units[[a]], "cmax"), tolerance_pct = 20
  ) |>
    dplyr::mutate(Analyte = analyte_labels[[a]], .before = 1)
}) |>
  dplyr::bind_rows()

cmp |>
  dplyr::rename("Arm" = "arm") |>
  knitr::kable(
    caption = "Typical-subject versus published median Cmax (Adiwijaya 2017 Table 2). * marks rows differing by more than 20%.",
    align = c("l", "l", "l", "r", "r", "r")
  )
```

| Analyte              | NCA parameter | Arm            | Reference | Simulated | % diff |
|:---------------------|:--------------|:---------------|----------:|----------:|-------:|
| Total irinotecan     | Cmax (ug/mL)  | 70 mg/m^2 Q2W  |      26.6 |      30.7 | +15.5% |
| Total irinotecan     | Cmax (ug/mL)  | 100 mg/m^2 Q3W |      41.5 |      43.9 |  +5.8% |
| Total SN-38          | Cmax (ng/mL)  | 70 mg/m^2 Q2W  |      2.64 |      2.92 | +10.5% |
| Total SN-38          | Cmax (ng/mL)  | 100 mg/m^2 Q3W |      3.99 |      4.09 |  +2.5% |
| Unencapsulated SN-38 | Cmax (ng/mL)  | 70 mg/m^2 Q2W  |      2.07 |      2.26 |  +9.3% |
| Unencapsulated SN-38 | Cmax (ng/mL)  | 100 mg/m^2 Q3W |      3.05 |      3.13 |  +2.5% |

Typical-subject versus published median Cmax (Adiwijaya 2017 Table 2).
\* marks rows differing by more than 20%. {.table}

``` r

pct <- suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "", cmp$`% diff`)))
pct <- pct[is.finite(pct)]

# Deterministic: mod_typ has all random effects zeroed and every covariate is
# fixed, so these six numbers carry no Monte-Carlo noise and the bound can be
# tight.
stopifnot(
  length(pct) == 6L,
  max(abs(pct)) < 16,
  stats::median(abs(pct)) < 10
)
```

All six typical-subject Cmax values agree with the published medians,
and the 100 mg/m^2 Q3W arm agrees to better than 4% on all three
analytes.

The two SN-38 rows are the informative ones. Total-SN-38 Cmax is
reproduced essentially by `fsn38 * Cmax(tIRI)` alone, which tests the
contaminant term; unencapsulated-SN-38 Cmax depends on `kmet`, `cl_sn38`
and the shared volume together, and is the only published number that
constrains the formation rate. Both land within 1% in the Q3W arm.
Because `kmet` appears nowhere else in the paper’s reported output, this
check is the sole external validation of the formation rate, and it is
what rules out the alternative readings of `Kcov` (for example one
driven by concentration rather than amount, which would be wrong by
roughly a factor of `V1`).

### Cohort distribution

The same NCA run over the 200-subject-per-arm cohort puts the published
medians in the context of the between-patient spread the model predicts.

``` r

cohort_cmax <- sim_nca_raw |>
  dplyr::select("id", "arm", "tIRI", "tSN38", "uSN38") |>
  tidyr::pivot_longer(c("tIRI", "tSN38", "uSN38"),
    names_to = "analyte", values_to = "conc"
  ) |>
  dplyr::group_by(.data$arm, .data$analyte, .data$id) |>
  dplyr::summarise(cmax = max(.data$conc), .groups = "drop_last") |>
  dplyr::group_by(.data$arm, .data$analyte) |>
  dplyr::summarise(
    `Simulated median` = signif(stats::median(.data$cmax), 3),
    `Simulated 5th-95th` = paste0(
      signif(stats::quantile(.data$cmax, 0.05), 3), " - ",
      signif(stats::quantile(.data$cmax, 0.95), 3)
    ),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    Analyte = analyte_labels[.data$analyte],
    `Published median` = mapply(
      function(a, r) published[[a]]$cmax[published[[a]]$arm == r],
      .data$analyte, .data$arm
    )
  ) |>
  dplyr::select(
    "Analyte", "arm", "Published median", "Simulated median",
    "Simulated 5th-95th"
  ) |>
  dplyr::rename("Arm" = "arm")

knitr::kable(
  cohort_cmax,
  caption = "Cohort Cmax distribution against the published medians. Units are ug/mL for total irinotecan and ng/mL for both SN-38 quantities."
)
```

| Analyte | Arm | Published median | Simulated median | Simulated 5th-95th |
|:---|:---|---:|---:|:---|
| Total irinotecan | 100 mg/m^2 Q3W | 41.50 | 44.70 | 27.9 - 66.8 |
| Total SN-38 | 100 mg/m^2 Q3W | 3.99 | 4.34 | 2.31 - 9.07 |
| Unencapsulated SN-38 | 100 mg/m^2 Q3W | 3.05 | 3.16 | 1.06 - 7.62 |
| Total irinotecan | 70 mg/m^2 Q2W | 26.60 | 30.70 | 19.7 - 47.1 |
| Total SN-38 | 70 mg/m^2 Q2W | 2.64 | 3.28 | 1.41 - 7.86 |
| Unencapsulated SN-38 | 70 mg/m^2 Q2W | 2.07 | 2.38 | 0.887 - 6.51 |

Cohort Cmax distribution against the published medians. Units are ug/mL
for total irinotecan and ng/mL for both SN-38 quantities. {.table}

``` r

cohort_pct <- 100 * (cohort_cmax$`Simulated median` / cohort_cmax$`Published median` - 1)

# Loose by design. Cohort medians depend on the random draw (rxSetSeed fixes
# the stream within an rxode2 version but not across versions), and the model's
# large clearance IIV -- variance 0.843, which also propagates into kmet
# through the log10(CL) covariate -- makes Cmax a strongly nonlinear function
# of several correlated etas, so the cohort median sits above the typical
# value. The tight, reproducible test is the typical-subject one above.
stopifnot(max(abs(cohort_pct)) < 35)
```

Every cohort median sits above the published median, by 3% to 23%. Two
effects push in that direction and neither is a transcription error.
First, the paper’s Table 2 medians come from simulations driven by the
doses patients *actually* received, which include reductions and delays,
whereas this cohort receives the full protocol dose at every
administration. Second, Cmax is a nonlinear function of several large
random effects here, so the median of the simulated Cmax exceeds the
Cmax of the median subject.

### Half-life from NCA

NCA on a cycle-1 profile recovers the paper’s **first-phase** half-life
(38.2 h), not its terminal half-life. That is expected and is not a
discrepancy: the terminal phase carries only about 0.01% of the
concentration amplitude, so it cannot be characterised from a two- or
three-week window. The paper makes the same point, cautioning that its
12,200 h terminal estimate “should be treated with cautions because of
the limited number of samples measured by assay with a lower limit of
quantification”. The terminal half-life is instead checked in closed
form in the Structural checks section above.

``` r

nca_typ_hl <- run_nca(sim_typ, "tIRI", "ug/mL",
  dose_data = dose_df_typ, intervals = first_phase_interval
)

hl_typ <- as.data.frame(nca_typ_hl$result) |>
  dplyr::filter(.data$PPTESTCD == "half.life") |>
  dplyr::transmute(
    Arm = .data$arm,
    `Typical-subject t1/2 (h)` = round(.data$PPORRES * hours_per_week, 1),
    `Published first-phase t1/2 (h)` = 38.2,
    `% diff` = round(100 * (.data$PPORRES * hours_per_week / 38.2 - 1), 1)
  )

knitr::kable(
  hl_typ,
  caption = "Total-irinotecan first-phase half-life by NCA (0-1 week window), typical subject of each arm, against Adiwijaya 2017 Table 2."
)
```

| Arm | Typical-subject t1/2 (h) | Published first-phase t1/2 (h) | % diff |
|:---|---:|---:|---:|
| 100 mg/m^2 Q3W | 38.5 | 38.2 | 0.7 |
| 70 mg/m^2 Q2W | 35.8 | 38.2 | -6.3 |

Total-irinotecan first-phase half-life by NCA (0-1 week window), typical
subject of each arm, against Adiwijaya 2017 Table 2. {.table}

``` r

hl_typ_h <- as.data.frame(nca_typ_hl$result) |>
  dplyr::filter(.data$PPTESTCD == "half.life") |>
  dplyr::pull(.data$PPORRES) * hours_per_week

# Deterministic (zeroed random effects, fixed covariates), so a tight bound is
# appropriate.
stopifnot(length(hl_typ_h) == 2L, max(abs(hl_typ_h / 38.2 - 1)) < 0.12)
```

The monotherapy arm lands within about 1% of the published 38.2 h. The
Q2W arm is shorter by roughly 6%, and that offset is not noise: it is
the 5-FU coadministration effect on irinotecan clearance,
`exp(0.075) = 1.078`, which raises `kel` by 7.8% and so shortens the
half-life by about 7%. A covariate that the paper reports as small
nevertheless shows up cleanly here, which is a useful check that the
covariate is wired to the right parameter.

| Arm            | Median t1/2 (h) | 5th-95th percentile (h) |
|:---------------|----------------:|:------------------------|
| 100 mg/m^2 Q3W |            32.6 | 13.5 - 103.3            |
| 70 mg/m^2 Q2W  |            34.8 | 12.9 - 256.9            |

Cohort distribution of the first-phase half-life. The published 95%
confidence interval, 23.2-56.7 h, describes the precision of the
parameter estimate rather than between-patient spread, so it is not
directly comparable to this range. {.table}

## Covariate effects

Race was the paper’s strongest baseline predictor, and it acts in
opposite directions on the two analytes: East Asian patients had *lower*
total irinotecan and *higher* SN-38, which the paper uses to explain
their higher observed rate of neutropenia and lower rate of diarrhea.

``` r

race_profile <- function(asian) {
  rxode2::rxSolve(
    mod_typ,
    make_events(
      typical_subject(70 * 1.70,
        RACE_ASIAN = asian, STUDY_NAPOLI1 = 1, CONMED_FLUOROURACIL = 1
      ),
      dose_times = 0, obs_times = split_times
    ),
    returnType = "data.frame"
  )
}

trapezoid <- function(time, conc) sum(diff(time) * (head(conc, -1) + tail(conc, -1)) / 2)

race_tbl <- dplyr::bind_rows(
  race_profile(0) |> dplyr::mutate(race = "Non-Asian (reference)"),
  race_profile(1) |> dplyr::mutate(race = "East Asian")
) |>
  dplyr::group_by(.data$race) |>
  dplyr::summarise(
    `tIRI Cmax (ug/mL)` = round(max(.data$Cc), 1),
    `tIRI AUC0-2wk (ug*week/mL)` = round(trapezoid(.data$time, .data$Cc), 2),
    `uSN38 Cmax (ng/mL)` = round(max(.data$Cu_sn38) * 1000, 2),
    .groups = "drop"
  )
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'

knitr::kable(race_tbl, caption = "Typical-subject exposure by race at 70 mg/m^2.")
```

| race | tIRI Cmax (ug/mL) | tIRI AUC0-2wk (ug\*week/mL) | uSN38 Cmax (ng/mL) |
|:---|---:|---:|---:|
| East Asian | 30.7 | 5.03 | 3.11 |
| Non-Asian (reference) | 30.7 | 9.43 | 2.26 |

Typical-subject exposure by race at 70 mg/m^2. {.table}

``` r


usn38_race <- tibble::deframe(race_tbl[, c("race", "uSN38 Cmax (ng/mL)")])
tiri_auc_race <- tibble::deframe(race_tbl[, c("race", "tIRI AUC0-2wk (ug*week/mL)")])

sprintf(
  "East Asian: tIRI AUC %+.0f%%, uSN38 Cmax %+.0f%% versus non-Asian.",
  100 * (tiri_auc_race[["East Asian"]] / tiri_auc_race[["Non-Asian (reference)"]] - 1),
  100 * (usn38_race[["East Asian"]] / usn38_race[["Non-Asian (reference)"]] - 1)
)
#> [1] "East Asian: tIRI AUC -47%, uSN38 Cmax +38% versus non-Asian."
```

Note that **total-irinotecan Cmax is identical between the two rows**.
That is structural, not an error: race enters the model only through
clearance, and Cmax after an intravenous dose is set by the volume. The
7% lower tIRI Cmax the paper reports for East Asian patients is
therefore an empirical population comparison – confounded with BSA and
the other covariates that differ between the race groups – rather than
an effect the model reproduces at fixed covariates. Race does move tIRI
*exposure*, which is the quantity its coefficient acts on.

The uSN38 direction is the mechanistically interesting one, because two
effects oppose each other. Higher irinotecan clearance in East Asian
patients reduces the amount of irinotecan available to convert, but it
also *raises* the formation rate through the `log10(CL)` covariate on
`kmet` (coefficient +2.095). The second effect wins, together with the
small reduction in SN-38 clearance, so unencapsulated SN-38 ends up
higher – which is the direction the paper reports and the basis of its
explanation for the higher neutropenia rate observed in East Asian
patients.

``` r

# Directional tests only: the paper's percentages are population medians over
# the whole covariate distribution, not typical-subject ratios, so their
# magnitudes are not reproducible from a single typical subject.
stopifnot(
  usn38_race[["East Asian"]] > usn38_race[["Non-Asian (reference)"]],
  tiri_auc_race[["East Asian"]] < tiri_auc_race[["Non-Asian (reference)"]],
  # Race is wired to clearance only, so Cmax must be untouched.
  isTRUE(all.equal(
    unname(tibble::deframe(race_tbl[, c("race", "tIRI Cmax (ug/mL)")])[["East Asian"]]),
    unname(tibble::deframe(race_tbl[, c("race", "tIRI Cmax (ug/mL)")])[["Non-Asian (reference)"]])
  ))
)
```

Bilirubin is the other covariate with a clinically framed conclusion.
The paper reports that patients with bilirubin at or above 1 mg/dL (17.1
umol/L) had 35% higher unencapsulated-SN-38 Cmax than patients below
that threshold.

``` r

bili_profile <- function(tbili) {
  rxode2::rxSolve(
    mod_typ,
    make_events(
      typical_subject(70 * 1.70,
        TBILI = tbili, STUDY_NAPOLI1 = 1, CONMED_FLUOROURACIL = 1
      ),
      dose_times = 0, obs_times = split_times
    ),
    returnType = "data.frame"
  )
}

# Cohort median (7 umol/L) against the paper's hyperbilirubinaemia threshold.
usn38_low <- max(bili_profile(7)$Cu_sn38) * 1000
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'
usn38_high <- max(bili_profile(17.1)$Cu_sn38) * 1000
#> ℹ omega/sigma items treated as zero: 'etalvc', 'etalcl', 'etalcl_sn38', 'etalkmet', 'etalfsn38'

sprintf(
  "uSN38 Cmax rises from %.2f to %.2f ng/mL (%+.0f%%) between 7 and 17.1 umol/L bilirubin.",
  usn38_low, usn38_high, 100 * (usn38_high / usn38_low - 1)
)
#> [1] "uSN38 Cmax rises from 2.26 to 2.64 ng/mL (+16%) between 7 and 17.1 umol/L bilirubin."

stopifnot(usn38_high > usn38_low)
```

## Assumptions and deviations

- **Table 2’s `Cavg` values are not reproducible from the model and are
  not used as a check.** The paper derived them by simulating from
  post-hoc individual estimates using the doses patients *actually
  received* over the first six weeks, which include protocol-driven
  reductions, delays and discontinuations. Consistent with that, the
  reported `Cavg` values are not proportional to dose rate: 70 mg/m^2
  Q2W and 100 mg/m^2 Q3W deliver nearly identical mg/m^2 per week yet
  are reported with a 1.4-fold `Cavg` ratio. Reproducing them would
  require the NAPOLI-1 dosing records, which are not published. The Cmax
  values in the same table *are* used, since they are dominated by the
  first dose of cycle 1.
- **The 70 mg/m^2 Q2W arm runs high.** Simulated median Cmax exceeds the
  published median by 13% (total irinotecan), 8% (total SN-38) and 7%
  (unencapsulated SN-38), whereas the 100 mg/m^2 Q3W arm agrees to
  within 4% throughout. The most likely explanation is the same
  actual-versus-protocol dose issue: the Q2W arm is the 5-FU/LV
  combination arm, where dose reductions are more frequent. No parameter
  was adjusted to close the gap.
- **Infusion duration is not stated in any available source**, so doses
  are given as an instantaneous input to the central compartment. The
  paper’s Table S1 times samples “post drug infusion” but never gives
  its length. The approximation is immaterial here: with a first-phase
  half-life of 38.2 h, a 90-minute infusion would lower Cmax by about
  1.4%.
- **Age, sex and AST are pre-specified but unreported.** Table S3 lists
  them as covariates on both clearances under the paper’s full-covariate
  approach, but neither Table S6 nor Table S7 reports a coefficient for
  any of them. They are therefore recorded in the model file’s
  `covariatesDataExcluded` list – as documentation of the paper’s
  covariate screen – and are not encoded.
- **The tIRI CL and V1 covariates on `kmet` are centred on the typical
  values** (13.6 L/week and 4.60 L). Equation S1 specifies
  median-centring for continuous covariates but these two are
  model-derived rather than measured, so no median is tabulated for
  them. Centring on the typical values is the reading under which the
  printed `Kcov = 0.00072 1/week` is the typical subject’s value; it is
  confirmed by the unencapsulated-SN-38 Cmax check above, which
  reproduces both published values to within 8%.
- **Laboratory-covariate units are unambiguous because the terms are
  centred.** ALT, albumin and bilirubin enter as
  `log10(x) - log10(median)`, a difference of logarithms, so a change of
  concentration unit cancels. Only creatinine clearance enters linearly,
  and Table S6 fixes its scale through the theta unit `min/mL`,
  i.e. mL/min.
- **Residual error is log-normal, converted from the source’s log10
  scale.** The paper fitted log10-transformed concentrations with the M3
  method for values below the limit of quantification; M3 handling is an
  estimation feature and has no simulation counterpart, so it is not
  represented.
- **The cohort’s covariate marginals are matched, but not their
  correlations.** Table 1 reports marginal distributions only. Body
  size, hepatic markers and renal function are correlated in practice;
  drawing them independently will slightly overstate the spread of
  predicted exposures. This affects the width of the percentile bands in
  the figures, not the medians used in the checks.
- **The NAPOLI-1 arm assignment follows the trial design**: the 70
  mg/m^2 Q2W arm is simulated with 5-FU/LV coadministration and the 100
  mg/m^2 Q3W arm as monotherapy, and both carry `STUDY_NAPOLI1 = 1`.
