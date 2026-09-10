# Sitafloxacin (Paiboonvong 2025)

## Model and source

``` r

mod <- readModelDb("Paiboonvong_2025_sitafloxacin")
ui  <- rxode2::rxode(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Paiboonvong T, Montakantikul P, Panjasawatwong N, Singkham
  N, Punyawudho B. Population pharmacokinetics and pharmacodynamics of
  sitafloxacin in plasma and alveolar epithelial lining fluid of
  critically ill Thai patients with pneumonia. Pharmacol Res Perspect.
  2025;13(2):e70081. <doi:10.1002/prp2.70081>. PMCID: PMC11930543. Final
  parameter estimates: Table 2. Model structure, the transit-absorption
  chain with KA = KTR, the ELF compartment and the age-on-F
  relationship: Results section 3.3 and Figure 1. Structural encoding
  including the placement of the unbound fraction inside the
  central-to-ELF micro-rate constant and the omega / sigma values on
  their estimation scales: Data S1 (Supporting Information), NONMEM
  control streams Run56 (estimation) and Run59 (simulation). Underlying
  clinical study, sampling design and bioanalysis: reference 7 of the
  paper. See also modellib(‘Wu_2025_sitafloxacin’) and
  modellib(‘Rodjun_2023_sitafloxacin’) for independent sitafloxacin
  popPK models in non-critically-ill populations.
- Description: Population PK model for oral sitafloxacin in plasma and
  pulmonary epithelial lining fluid (ELF) of 12 critically ill Thai
  patients with pneumonia, developed by Paiboonvong 2025 to drive Monte
  Carlo probability-of-target-attainment simulation against
  Streptococcus pneumoniae. Plasma disposition is one-compartment;
  absorption is a Savic transit chain of one transit compartment in
  which the absorption rate constant is constrained to equal the transit
  rate constant (ka = ktr = (ntr + 1) / mtt with ntr fixed to 1). An ELF
  compartment hangs off the central compartment as a peripheral
  compartment with a fixed physiological volume of 0.025 L and an
  inter-compartmental clearance q_elf. IMPORTANT: in the authors’ own
  control stream the central-to-ELF micro-rate constant carries the
  plasma unbound fraction, k_central_elf = q_elf \* fu \* pcelf / vc
  while k_elf_central = q_elf / v_elf, so the estimated partition
  coefficient pcelf = 0.772 is the ELF-to-UNBOUND-plasma ratio and the
  ELF-to-total-plasma concentration (and AUC) ratio the model actually
  produces is fu \* pcelf = 0.63 \* 0.772 = 0.486. Body weight is an
  a-priori allometric covariate on CL/F (exponent 0.75 fixed), V/F
  (exponent 1 fixed) and Q/F_ELF (exponent 0.75 fixed), all normalized
  to the cohort median 52 kg. Age is a LINEAR covariate on relative
  bioavailability, F = 1 + 0.0258 \* (AGE - 57), which is only usable
  over the adult age range the authors simulated (30-70 years): the
  expression reaches zero at 18.2 years and is negative below it. All
  disposition parameters are apparent (/F): only oral data were analysed
  and F was fixed to 1 with its inter-individual variability estimated.
  Residual variability was additive on the natural-log concentration
  scale for both plasma and ELF, encoded here as lnorm().
- Article: <https://doi.org/10.1002/prp2.70081>
- Supporting Information (Data S1, NONMEM control streams Run56 and
  Run59): <https://europepmc.org/article/MED/40128137>

## Population

Twelve critically ill Thai patients with pneumonia, admitted to an
intensive care unit, each given a single 200 mg oral dose of
sitafloxacin under fasting conditions (Methods 2.1). Six of the twelve
were male. Median (IQR) age was 57 (40-65) years, body weight 52 (44-68)
kg and Cockcroft-Gault creatinine clearance 68 (30-96) mL/min; median
APACHE II score was 21 (18-33) and median serum albumin 2.0 (1.8-2.2)
g/dL (Results 3.1). Two thirds of the cohort had mild-to-moderate renal
impairment and one patient showed augmented renal clearance at 235
mL/min.

Plasma was sampled pre-dose and at 0.5, 1, 2, 3, 8 and 12 h; each
patient additionally contributed one bronchoalveolar-lavage sample,
randomly assigned to the 0.5-2, 3-4, 5-6 or 7-9 h window. The model was
built on 83 plasma and 12 epithelial lining fluid (ELF) concentrations,
one plasma sample below the 0.025 mg/L limit of quantification having
been excluded (Results 3.3).

The same information is available programmatically from the model’s
`population` metadata:

``` r

str(ui$population)
#> List of 10
#>  $ species       : chr "human"
#>  $ n_subjects    : int 12
#>  $ n_studies     : int 1
#>  $ age_median    : chr "57 years (IQR 40-65)"
#>  $ weight_median : chr "52 kg (IQR 44-68)"
#>  $ sex_female_pct: num 50
#>  $ disease_state : chr "Critically ill patients with pneumonia admitted to an intensive care unit. Median APACHE II score 21 (IQR 18-33"| __truncated__
#>  $ dose_range    : chr "Sitafloxacin 200 mg orally as a single dose under fasting conditions"
#>  $ regions       : chr "Thailand"
#>  $ notes         : chr "Baseline demographics: Results 3.1. The model was built on 83 plasma and 12 ELF concentrations; one plasma conc"| __truncated__
```

## Source trace

Every `ini()` entry in
`inst/modeldb/specificDrugs/Paiboonvong_2025_sitafloxacin.R` carries an
in-file comment naming its origin. They are collected here for review.
“Table 2” and “Results” refer to the main article; “Run59” refers to the
simulation control stream printed in Data S1 of the Supporting
Information, which carries the final parameter estimates (Run56, the
estimation stream, carries earlier initial values for two of them).

| Equation / parameter | Value | Source location |
|----|----|----|
| `lcl` (CL/F at 52 kg) | 7.03 L/h | Table 2, RSE 24.9%, SIR 95% CI 4.36-11.0; Run59 `THETA(1)` |
| `lvc` (V/F at 52 kg) | 116 L | Table 2, RSE 13.0%, SIR 95% CI 93.9-153; Run59 `THETA(2)` |
| `lfdepot` (F) | 1, fixed | Table 2 “1 fix”; Methods 2.2; Run59 `THETA(3) (1) FIX` |
| `lmtt` (MTT) | 1.48 h | Table 2, RSE 27.1%, SIR 95% CI 0.858-2.44; Run59 `THETA(4)` |
| `lq_elf` (Q/F ELF at 52 kg) | 0.0441 L/h | Table 2, RSE 36.5%, SIR 95% CI 0.0194-0.0812; Run59 `THETA(5)` |
| `lv_elf` (V/F ELF) | 0.025 L, fixed | Table 2 “0.025 fix”; Results 3.3 (lung physiological value, paper reference 11); Run59 `THETA(6) (0.025) FIX` |
| `lpcelf` (partition coefficient) | 0.772 | Table 2, RSE 18.2%, SIR 95% CI 0.556-1.08; Run59 `THETA(7)` |
| `fu` (fraction unbound) | 0.63, fixed | Table 2 “0.63 fix”; Methods 2.4; used inside Run59 `K24 = QELF*0.63*PC/VBD` |
| `e_age_fdepot` (age on F) | 0.0258 /year | Table 2 “AGE on F 2.58”, RSE 7.64%; Results 3.3 “F = 1 + (Age - 57) \* 0.0258”; Run59 `THETA(8)` |
| `e_wt_cl` | 0.75, fixed | Methods 2.2; Run59 `TVCLBD = THETA(1)*((WT/52)**0.75)` |
| `e_wt_vc` | 1, fixed | Methods 2.2; Run59 `TVVBD = THETA(2)*(WT/52)` |
| `e_wt_q_elf` | 0.75, fixed | Run59 `TVQELF = THETA(5)*((WT/52)**0.75)` and its header comment |
| `etalcl` | 0.566 (87.2 %CV) | Run59 `$OMEGA` 1; Table 2 IIV CL/F 87.2 %CV |
| `etalfdepot` | 0.0285 (17.0 %CV) | Run59 `$OMEGA` 3; Table 2 IIV F 17.0 %CV |
| `etalmtt` | 0.626 (93.3 %CV) | Run59 `$OMEGA` 4; Table 2 IIV MTT 93.3 %CV |
| `etalpcelf` | 0.0528 (23.3 %CV) | Run59 `$OMEGA` 7; Table 2 IIV PC 23.3 %CV |
| `expSd` (plasma) | sqrt(0.122) | Run59 `$SIGMA` 1; Table 2 sigma plasma 36.0 |
| `expSd_Celf` (ELF) | sqrt(0.249) | Run59 `$SIGMA` 2; Table 2 sigma ELF 53.2 |
| `ktr = (ntr + 1) / mtt`, `ntr = 1` | n/a | Results 3.3 “one fixed transit absorption compartment, which set Ka to be identical to Ktr”; Run59 `NN = 1`, `KTR = (NN+1)/MTT`, `K13 = KTR`, `K32 = KTR` |
| `d/dt(depot)`, `d/dt(transit1)`, `d/dt(central)`, `d/dt(elf)` | n/a | Run59 `$DES` `DADT(1)`, `DADT(3)`, `DADT(2)`, `DADT(4)`; Figure 1 schematic |
| `k_central_elf = q_elf * fu * pcelf / vc` | n/a | Figure 1 arrow labelled `(Q_ELF x F_u x PC)/V`; Run59 `K24 = QELF*0.63*PC/VBD` |
| `k_elf_central = q_elf / v_elf` | n/a | Figure 1 arrow labelled `Q_ELF/V_ELF`; Run59 `K42 = QELF/VELF` |
| `Cc ~ lnorm()`, `Celf ~ lnorm()` | n/a | Results 3.3 “additive residual error model on the logarithmic scale”; Run59 `$ERROR` `Y = IPRED + EPS(n)` with `IPRED = LOG(IPRED)` |

The IIV column of Table 2 is printed as %CV while the control stream
carries log-scale variances; the two agree exactly under
`CV = sqrt(exp(v) - 1)`, which is checked below. The same identity
converts the two sigma rows.

``` r

th <- setNames(ui$iniDf$est, ui$iniDf$name)
etas <- ui$iniDf |> dplyr::filter(!is.na(neta1), neta1 == neta2)
omega_chk <- tibble::tibble(
  parameter = etas$name,
  variance  = etas$est,
  cv_pct    = 100 * sqrt(exp(etas$est) - 1),
  published = c(87.2, 17.0, 93.3, 23.3)
)
sig_chk <- tibble::tibble(
  parameter = c("expSd", "expSd_Celf"),
  variance  = c(th[["expSd"]], th[["expSd_Celf"]])^2,
  cv_pct    = 100 * sqrt(exp(c(th[["expSd"]], th[["expSd_Celf"]])^2) - 1),
  published = c(36.0, 53.2)
)
scale_chk <- dplyr::bind_rows(omega_chk, sig_chk) |>
  dplyr::mutate(abs_diff = abs(cv_pct - published))
knitr::kable(scale_chk, digits = 3,
             caption = "Log-scale variances re-expressed as %CV against Table 2.")
```

| parameter  | variance | cv_pct | published | abs_diff |
|:-----------|---------:|-------:|----------:|---------:|
| etalcl     |    0.566 | 87.247 |      87.2 |    0.047 |
| etalfdepot |    0.028 | 17.003 |      17.0 |    0.003 |
| etalmtt    |    0.626 | 93.280 |      93.3 |    0.020 |
| etalpcelf  |    0.053 | 23.285 |      23.3 |    0.015 |
| expSd      |    0.122 | 36.021 |      36.0 |    0.021 |
| expSd_Celf |    0.249 | 53.173 |      53.2 |    0.027 |

Log-scale variances re-expressed as %CV against Table 2. {.table}

``` r


# Deterministic identity: no simulation, no RNG. A variance/SD mix-up on any
# row moves the %CV by tens of points, so 0.1 is a tight but safe bound.
stopifnot(max(scale_chk$abs_diff) < 0.1)
```

## Structural verification

Two exact identities follow from the ODE system and hold for every
subject, independently of the parameter values. They are the cheapest
available test that the compartments are wired the way the control
stream wires them.

1.  Because `elf` returns all of its drug to `central` (it is a true
    peripheral compartment, not a sink), total plasma exposure over one
    steady-state dosing interval is fixed by mass balance:
    `AUC(0-24, plasma) = 24 h / tau * dose * F / CL`.
2.  Because the central-to-ELF micro-constant carries `fu` while the
    ELF-to-central one does not, ELF equilibrates to `fu * pcelf` times
    the *total* plasma concentration, so
    `AUC(ELF) / AUC(plasma) = fu * pcelf` exactly.

``` r

make_ss_events <- function(subj, dose, tau = 12, n_doses = 16L,
                           obs = seq(168, 192, by = 0.5)) {
  dosing <- subj |>
    dplyr::mutate(time = 0, amt = dose, evid = 1L, cmt = "depot",
                  dvid = NA_integer_, ii = tau, addl = n_doses - 1L)
  # Observation rows address the ODE state `central` and select the endpoint
  # with `dvid = 1L`. rxode2 returns BOTH observables (Cc and Celf) as columns
  # on those rows; naming an algebraic observable in the compartment column
  # instead would inject a new compartment slot and renumber the ODE states.
  observing <- subj |>
    tidyr::crossing(time = obs) |>
    dplyr::mutate(amt = 0, evid = 0L, cmt = "central",
                  dvid = 1L, ii = 0, addl = 0L)
  dplyr::bind_rows(dosing, observing) |>
    dplyr::arrange(id, time, dplyr::desc(evid)) |>
    as.data.frame()
}

trapz <- function(t, y) sum(diff(t) * (utils::head(y, -1) + utils::tail(y, -1)) / 2)

# Typical-value subjects spanning the age range the paper simulated.
typ_subj <- tibble::tibble(id = 1:5, WT = 52, AGE = c(30, 40, 50, 60, 70),
                           treatment = "100 mg q12h")
typ_sim <- rxode2::rxSolve(
  rxode2::zeroRe(mod), events = make_ss_events(typ_subj, dose = 100),
  keep = c("WT", "AGE"), useLinCmt = FALSE, returnType = "data.frame"
) |>
  dplyr::filter(!is.na(Cc))
#> ℹ parameter labels from comments will be replaced by 'label()'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalfdepot', 'etalmtt', 'etalpcelf'
#> Warning: multi-subject simulation without without 'omega'

ident <- typ_sim |>
  dplyr::group_by(id, WT, AGE) |>
  dplyr::summarise(
    auc_plasma = trapz(time, Cc),
    auc_elf    = trapz(time, Celf),
    cl = dplyr::first(cl), fdepot = dplyr::first(fdepot),
    pcelf = dplyr::first(pcelf), .groups = "drop"
  ) |>
  dplyr::mutate(
    auc_closed_form = 2 * 100 * fdepot / cl,
    pct_diff_auc    = 100 * (auc_plasma / auc_closed_form - 1),
    ratio_elf       = auc_elf / auc_plasma,
    ratio_expected  = th[["fu"]] * pcelf
  )

ident |>
  dplyr::select(AGE, auc_plasma, auc_closed_form, pct_diff_auc,
                ratio_elf, ratio_expected) |>
  dplyr::rename(
    "Age (years)"                = AGE,
    "AUC0-24 plasma, solved"     = auc_plasma,
    "AUC0-24 plasma, closed form" = auc_closed_form,
    "% difference"               = pct_diff_auc,
    "AUC(ELF)/AUC(plasma)"       = ratio_elf,
    "fu * pcelf"                 = ratio_expected
  ) |>
  knitr::kable(digits = 5,
               caption = "Steady-state mass-balance identities, typical-value solve.")
```

| Age (years) | AUC0-24 plasma, solved | AUC0-24 plasma, closed form | % difference | AUC(ELF)/AUC(plasma) | fu \* pcelf |
|---:|---:|---:|---:|---:|---:|
| 30 | 8.63123 | 8.63158 | -0.00409 | 0.48638 | 0.48636 |
| 40 | 15.97090 | 15.97155 | -0.00409 | 0.48638 | 0.48636 |
| 50 | 23.31057 | 23.31152 | -0.00409 | 0.48638 | 0.48636 |
| 60 | 30.65024 | 30.65149 | -0.00409 | 0.48638 | 0.48636 |
| 70 | 37.98991 | 37.99147 | -0.00409 | 0.48638 | 0.48636 |

Steady-state mass-balance identities, typical-value solve. {.table
style="width:100%;"}

``` r


# Both sides of each comparison use the SAME drawn parameters, so the residual
# is pure trapezoidal / solver error and a tight bound is correct here. The ELF
# ratio carries slightly more of it than the plasma AUC because the ELF profile
# lags plasma, so the trapezoid over a finite interval is not exactly
# proportional; the realised relative error is about 4e-5. Dropping `fu` from
# the micro-constant would move this ratio from 0.486 to 0.772, a 59% error.
stopifnot(
  max(abs(ident$pct_diff_auc)) < 0.05,
  max(abs(ident$ratio_elf / ident$ratio_expected - 1)) < 1e-3
)
```

The second identity is the reason this vignette treats the published
partition coefficient carefully: the model’s ELF-to-*total*-plasma
exposure ratio is `fu * pcelf` = 0.63 \* 0.772 = 0.486, not the 0.772
the abstract quotes. See Assumptions and deviations.

## Virtual cohort

Original observed data are not publicly available. Three virtual cohorts
are used, each of 200 subjects (the per-arm cap for validation
vignettes).

- **Study cohort** - a 200 mg single dose, with age and weight drawn to
  match the medians and IQRs of Results 3.1. Used for the Figure 3
  visual predictive check.
- **100 mg q12h** and **50 mg q12h** - the two approved regimens, dosed
  to steady state. Ages cycle through the five values the paper
  simulated (30, 40, 50, 60 and 70 years) and weight is held at the 52
  kg reference, which is what reproduces the published median exposures
  (see Assumptions).

``` r

# set.seed() seeds R's RNG, which is what draws the covariates below. It does
# NOT seed rxode2's simulation RNG, and rxode2's streams are partitioned per
# solver thread, so the etas drawn downstream differ between a 2-core CI runner
# and a 16-thread workstation. Every assertion below is therefore either
# deterministic (no cohort involved) or written with headroom for that spread.
set.seed(20250926)
n_arm <- 200L

# Study cohort: lognormal weight with median 52 kg and IQR 44-68 kg; normal age
# with median 57 y and IQR 40-65 y, truncated to 25-85 y so that the LINEAR
# age-on-F relationship stays positive (it reaches zero at 18.2 years).
study_subj <- tibble::tibble(
  id  = seq_len(n_arm),
  WT  = 52 * exp(stats::rnorm(n_arm, 0, log(68 / 44) / (2 * stats::qnorm(0.75)))),
  AGE = pmin(85, pmax(25, stats::rnorm(n_arm, 57, 25 / (2 * stats::qnorm(0.75))))),
  treatment = "200 mg single dose"
)

make_ss_subj <- function(label, id_offset) {
  tibble::tibble(
    id  = id_offset + seq_len(n_arm),
    WT  = 52,
    AGE = rep(c(30, 40, 50, 60, 70), length.out = n_arm),
    treatment = label
  )
}
ss_subj <- dplyr::bind_rows(
  make_ss_subj("100 mg q12h", id_offset = 1000L),
  make_ss_subj("50 mg q12h",  id_offset = 2000L)
)

ss_events <- dplyr::bind_rows(
  make_ss_events(dplyr::filter(ss_subj, treatment == "100 mg q12h"), dose = 100),
  make_ss_events(dplyr::filter(ss_subj, treatment == "50 mg q12h"),  dose = 50)
)
stopifnot(!anyDuplicated(ss_events[, c("id", "time", "evid")]))

# Study cohort event table: single dose, sampled over the observed 0-12 h window.
study_events <- dplyr::bind_rows(
  study_subj |>
    dplyr::mutate(time = 0, amt = 200, evid = 1L, cmt = "depot",
                  dvid = NA_integer_, ii = 0, addl = 0L),
  study_subj |>
    tidyr::crossing(time = seq(0, 12, by = 0.25)) |>
    dplyr::mutate(amt = 0, evid = 0L, cmt = "central",
                  dvid = 1L, ii = 0, addl = 0L)
) |>
  dplyr::arrange(id, time, dplyr::desc(evid)) |>
  as.data.frame()
stopifnot(!anyDuplicated(study_events[, c("id", "time", "evid")]))
```

## Simulation

``` r

rxode2::rxSetSeed(20250926)

study_sim <- rxode2::rxSolve(
  mod, events = study_events, keep = c("WT", "AGE", "treatment"),
  useLinCmt = FALSE, returnType = "data.frame"
)
#> ℹ parameter labels from comments will be replaced by 'label()'

ss_sim <- rxode2::rxSolve(
  mod, events = ss_events, keep = c("WT", "AGE", "treatment"),
  useLinCmt = FALSE, returnType = "data.frame"
)

stopifnot(nrow(study_sim) > 0, nrow(ss_sim) > 0,
          all(study_sim$Cc[!is.na(study_sim$Cc)] >= 0),
          all(ss_sim$Celf[!is.na(ss_sim$Celf)] >= 0))
```

## Replicate published figures

``` r

vpc_dat <- study_sim |>
  dplyr::filter(!is.na(Cc)) |>
  tidyr::pivot_longer(c(Cc, Celf), names_to = "matrix", values_to = "conc") |>
  dplyr::mutate(matrix = dplyr::recode(matrix,
                                       Cc = "A. Plasma", Celf = "B. ELF")) |>
  dplyr::group_by(matrix, time) |>
  dplyr::summarise(
    Q05 = stats::quantile(conc, 0.05),
    Q50 = stats::quantile(conc, 0.50),
    Q95 = stats::quantile(conc, 0.95),
    .groups = "drop"
  )

ggplot(vpc_dat, aes(time, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~matrix) +
  scale_y_log10() +
  labs(x = "Time after dose (h)", y = "Sitafloxacin concentration (mg/L)",
       caption = paste("Replicates Figure 3 of Paiboonvong 2025.",
                       "Line = median, band = 5th-95th percentile,",
                       "200 virtual subjects given 200 mg orally."))
#> Warning in scale_y_log10(): log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
#> log-10 transformation introduced infinite values.
```

![Replicates Figure 3 of Paiboonvong 2025: visual predictive check of
plasma (A) and ELF (B) sitafloxacin concentrations after a single 200 mg
oral
dose.](Paiboonvong_2025_sitafloxacin_files/figure-html/figure-3-1.png)

Replicates Figure 3 of Paiboonvong 2025: visual predictive check of
plasma (A) and ELF (B) sitafloxacin concentrations after a single 200 mg
oral dose.

Figure 3 of the paper plots concentrations on a ug/L axis. Its plasma
panel shows a broad, flat peak: the median climbs to roughly 1000 ug/L
(1 mg/L) between about 2 and 5 h and is still near 700-900 ug/L at 12 h,
which is the delayed, prolonged absorption the Discussion attributes to
gastrointestinal dysmotility in critical illness. Its ELF panel sits
near 500-600 ug/L across the 1-9 h bronchoalveolar-lavage window,
roughly half of plasma. Both features are reproduced above.

``` r

peak_chk <- study_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::group_by(time) |>
  dplyr::summarise(med_plasma = stats::median(Cc),
                   med_elf = stats::median(Celf), .groups = "drop")
tmax_med <- peak_chk$time[which.max(peak_chk$med_plasma)]
cmax_med <- max(peak_chk$med_plasma)
# ELF-to-plasma concentration ratio once the ELF compartment has equilibrated
# (its half-life of equilibration is about 0.3 h, so 4 h is well past it).
ratio_med <- stats::median(
  study_sim$Celf[!is.na(study_sim$Cc) & study_sim$time >= 4] /
    study_sim$Cc[!is.na(study_sim$Cc) & study_sim$time >= 4]
)
tibble::tibble(
  quantity = c("Median Cmax (mg/L)", "Time of median peak (h)",
               "Median Celf/Cc after 4 h"),
  value    = c(cmax_med, tmax_med, ratio_med)
) |>
  dplyr::rename("Quantity" = quantity, "Value" = value) |>
  knitr::kable(digits = 3,
               caption = "Summary of the simulated single-dose median profile.")
```

| Quantity                 | Value |
|:-------------------------|------:|
| Median Cmax (mg/L)       | 1.102 |
| Time of median peak (h)  | 4.000 |
| Median Celf/Cc after 4 h | 0.501 |

Summary of the simulated single-dose median profile. {.table}

``` r


# Cohort-derived, so the bounds carry headroom for the per-thread RNG spread.
# The typical-value profile peaks at 1.44 mg/L at 3.6 h and the paper's own
# Figure 3 median peaks near 1 mg/L between 2 and 5 h; a mis-transcribed
# volume, dose or MTT moves these by a factor, not by a few percent.
stopifnot(cmax_med > 0.4, cmax_med < 3.0, tmax_med > 1, tmax_med < 9)
# Structural, and therefore tight: the equilibrated ratio is fu * pcelf with
# only the 23.3 %CV on pcelf around it. Dropping fu would put it at 0.772.
stopifnot(ratio_med > 0.42, ratio_med < 0.56)
```

## PKNCA validation

Non-compartmental analysis of the last steady-state dosing interval
(168-192 h), run once over a stacked frame in which each matrix is
carried as its own treatment level so that plasma and ELF appear in a
single comparison table.

``` r

# Only `!is.na(Cc)` -- a `time > 0` or `Cc > 0` filter would drop rows PKNCA
# needs to anchor the interval.
ss_long <- ss_sim |>
  dplyr::filter(!is.na(Cc)) |>
  dplyr::select(id, time, treatment, Cc, Celf) |>
  tidyr::pivot_longer(c(Cc, Celf), names_to = "matrix", values_to = "conc") |>
  dplyr::mutate(
    matrix    = dplyr::recode(matrix, Cc = "plasma", Celf = "ELF"),
    treatment = paste(treatment, matrix, sep = " - ")
  ) |>
  dplyr::rename(Cc = conc) |>
  dplyr::select(id, time, Cc, treatment) |>
  dplyr::arrange(treatment, id, time)

conc_obj <- PKNCA::PKNCAconc(ss_long, Cc ~ time | treatment + id,
                             concu = "mg/L", timeu = "h")

dose_df <- ss_events |>
  dplyr::filter(evid == 1) |>
  dplyr::select(id, amt, treatment, ii, addl) |>
  tidyr::crossing(dose_index = 0:15) |>
  dplyr::mutate(time = dose_index * ii) |>
  dplyr::filter(dose_index <= addl) |>
  dplyr::select(id, time, amt, treatment) |>
  tidyr::crossing(matrix = c("plasma", "ELF")) |>
  dplyr::mutate(treatment = paste(treatment, matrix, sep = " - ")) |>
  dplyr::select(id, time, amt, treatment)

dose_obj <- PKNCA::PKNCAdose(dose_df, amt ~ time | treatment + id, doseu = "mg")

intervals <- data.frame(
  start = 168, end = 192,
  cmax = TRUE, tmax = TRUE, cmin = TRUE, auclast = TRUE, cav = TRUE
)

nca_res <- PKNCA::pk.nca(PKNCA::PKNCAdata(conc_obj, dose_obj,
                                          intervals = intervals))
knitr::kable(summary(nca_res),
             caption = "Steady-state NCA over the 168-192 h dosing interval.")
```

| Interval Start | Interval End | treatment | N | AUClast (h\*mg/L) | Cmax (mg/L) | Cmin (mg/L) | Tmax (h) | Cav (mg/L) |
|---:|---:|:---|:---|:---|:---|:---|:---|:---|
| 168 | 192 | 100 mg q12h - ELF | 200 | 9.44 \[110\] | 0.496 \[92.8\] | 0.269 \[173\] | 15.5 \[1.50, 20.0\] | 0.393 \[110\] |
| 168 | 192 | 100 mg q12h - plasma | 200 | 18.5 \[105\] | 0.996 \[86.8\] | 0.506 \[175\] | 15.0 \[0.500, 19.5\] | 0.772 \[105\] |
| 168 | 192 | 50 mg q12h - ELF | 200 | 5.12 \[126\] | 0.265 \[108\] | 0.150 \[184\] | 15.5 \[2.00, 20.0\] | 0.213 \[126\] |
| 168 | 192 | 50 mg q12h - plasma | 200 | 10.4 \[116\] | 0.552 \[97.3\] | 0.294 \[178\] | 15.0 \[1.00, 19.5\] | 0.435 \[116\] |

Steady-state NCA over the 168-192 h dosing interval. {.table
style="width:100%;"}

### Comparison against published NCA

Results 3.4 reports the median steady-state 24 h AUC in plasma and ELF
for both regimens. Note that the paper labels these values “fAUC” but,
as shown under Assumptions and deviations, they are the model’s *total*
plasma AUC and the ELF AUC with no further unbound-fraction
multiplication.

``` r

published <- tibble::tribble(
  ~treatment,             ~auclast,
  "100 mg q12h - plasma", 21.02,
  "100 mg q12h - ELF",    10.17,
  "50 mg q12h - plasma",  10.42,
  "50 mg q12h - ELF",      5.12
)

cmp <- nlmixr2lib::ncaComparisonTable(
  simulated     = nca_res,
  reference     = published,
  by            = "treatment",
  units         = c(auclast = "mg.h/L"),
  tolerance_pct = 20
)

knitr::kable(cmp, caption = paste(
  "Simulated (200 subjects per arm) vs. published median steady-state AUC0-24.",
  "* differs from the reference by more than 20%."
), align = c("l", "l", "r", "r", "r"))
```

| NCA parameter    | treatment            | Reference | Simulated | % diff |
|:-----------------|:---------------------|----------:|----------:|-------:|
| AUClast (mg.h/L) | 100 mg q12h - plasma |        21 |      18.1 | -13.7% |
| AUClast (mg.h/L) | 100 mg q12h - ELF    |      10.2 |      9.18 |  -9.7% |
| AUClast (mg.h/L) | 50 mg q12h - plasma  |      10.4 |      11.3 |  +8.8% |
| AUClast (mg.h/L) | 50 mg q12h - ELF     |      5.12 |      5.37 |  +4.9% |

Simulated (200 subjects per arm) vs. published median steady-state
AUC0-24. \* differs from the reference by more than 20%. {.table}

``` r


# The simulated side is a 200-subject median of a 5-component age mixture whose
# log-scale spread is about 0.93, giving a median standard error near 8%. The
# bound below is roughly 3.5 of those and still goes red on any mis-transcribed
# clearance, dose or unit, which move the median by tens of percent. Realised
# values on one 16-thread run were +0.3 / -4.3 / -3.2 / -5.5 percent; do NOT
# tighten to those, they are one draw. The deterministic reproduction in the
# next section is the tight gate.
cmp_pct <- suppressWarnings(as.numeric(gsub("[^0-9.eE+-]", "", cmp$`% diff`)))
stopifnot(max(abs(cmp_pct), na.rm = TRUE) < 30)
```

## Deterministic reproduction of the published exposures

The Monte Carlo medians of Results 3.4 can be reproduced exactly,
without any simulation, because the steady-state 24 h plasma AUC of
subject *i* is `24/tau * dose * F_i / CL_i`, which is log-normal within
each of the five simulated age groups. The population median is then the
median of an equally weighted five-component log-normal mixture,
obtained below by root-finding on the mixture CDF. Every input is read
out of the packaged model, so the check goes red if any shipped value
drifts.

``` r

cl_typ <- exp(th[["lcl"]])
fu_typ <- th[["fu"]]
pc_typ <- exp(th[["lpcelf"]])
age_slope <- th[["e_age_fdepot"]]
var_of <- function(nm) etas$est[match(nm, etas$name)]

ages <- c(30, 40, 50, 60, 70)
f_age <- 1 + age_slope * (ages - 57)
# AUC = 24/tau * dose * F / CL varies with eta_F and eta_CL; the ELF AUC picks
# up eta_PC as well.
sd_plasma <- sqrt(var_of("etalfdepot") + var_of("etalcl"))
sd_elf    <- sqrt(var_of("etalfdepot") + var_of("etalcl") + var_of("etalpcelf"))

mixture_median <- function(mu, sdlog) {
  stats::uniroot(function(x) mean(stats::pnorm((x - mu) / sdlog)) - 0.5,
                 c(-20, 20))$root |> exp()
}

predicted <- function(dose_per_24h, matrix) {
  scale <- if (matrix == "plasma") 1 else fu_typ * pc_typ
  sdlog <- if (matrix == "plasma") sd_plasma else sd_elf
  mixture_median(log(dose_per_24h * f_age / cl_typ) + log(scale), sdlog)
}

exposure_chk <- tibble::tribble(
  ~regimen,       ~matrix,  ~dose24, ~published,
  "100 mg q12h",  "plasma", 200,     21.02,
  "100 mg q12h",  "ELF",    200,     10.17,
  "50 mg q12h",   "plasma", 100,     10.42,
  "50 mg q12h",   "ELF",    100,      5.12
) |>
  dplyr::rowwise() |>
  dplyr::mutate(predicted = predicted(dose24, matrix)) |>
  dplyr::ungroup() |>
  dplyr::mutate(pct_diff = 100 * (predicted / published - 1))

exposure_chk |>
  dplyr::select(-dose24) |>
  dplyr::rename(
    "Regimen"                        = regimen,
    "Matrix"                         = matrix,
    "Published median (mg.h/L)"      = published,
    "Model median (mg.h/L)"          = predicted,
    "% difference"                   = pct_diff
  ) |>
  knitr::kable(digits = c(0, 0, 2, 3, 2), caption = paste(
    "Published median steady-state AUC0-24 (Results 3.4) against the exact",
    "population median implied by the packaged parameters. No simulation is",
    "involved, so this comparison is identical on every machine."
  ))
```

| Regimen | Matrix | Published median (mg.h/L) | Model median (mg.h/L) | % difference |
|:---|:---|---:|---:|---:|
| 100 mg q12h | plasma | 21.02 | 21.067 | 0.22 |
| 100 mg q12h | ELF | 10.17 | 10.230 | 0.59 |
| 50 mg q12h | plasma | 10.42 | 10.533 | 1.09 |
| 50 mg q12h | ELF | 5.12 | 5.115 | -0.10 |

Published median steady-state AUC0-24 (Results 3.4) against the exact
population median implied by the packaged parameters. No simulation is
involved, so this comparison is identical on every machine. {.table}

``` r


# Deterministic; the realised maximum is about 1.1%. 3% still admits the
# paper's own rounding while going red on any parameter drift.
stopifnot(max(abs(exposure_chk$pct_diff)) < 3)
```

All four published medians are reproduced to within about 1%. That is
the strongest available evidence that the parameter values, the age
covariate, the allometric reference weight and the placement of `fu`
inside the ELF micro-constant have all been transcribed correctly:
getting any one of them wrong moves at least one of these four numbers
by tens of percent.

## Reproduction of the target-attainment table

Table 3 reports the probability of attaining `fAUC/MIC > 30`. The same
log-normal mixture gives each cell in closed form.

``` r

pta <- function(dose24, mic, matrix) {
  scale <- if (matrix == "plasma") 1 else fu_typ * pc_typ
  sdlog <- if (matrix == "plasma") sd_plasma else sd_elf
  mu <- log(dose24 * f_age / cl_typ) + log(scale)
  100 * mean(1 - stats::pnorm((log(30 * mic) - mu) / sdlog))
}

pta_chk <- tidyr::expand_grid(
  regimen = c("100 mg q12h", "50 mg q12h"),
  matrix  = c("plasma", "ELF"),
  mic     = c(0.032, 0.0625, 0.125, 0.25)
) |>
  dplyr::mutate(
    dose24    = ifelse(regimen == "100 mg q12h", 200, 100),
    published = c(100.00, 100.00, 94.60, 85.72,
                   99.18,  95.36, 79.44, 62.50,
                  100.00,  96.34, 81.42, 63.48,
                   95.34,  84.84, 55.70, 33.70)
  ) |>
  dplyr::rowwise() |>
  dplyr::mutate(model = pta(dose24, mic, matrix)) |>
  dplyr::ungroup() |>
  dplyr::mutate(
    abs_diff  = abs(model - published),
    # Both MIC = 0.125 cells of each matrix are recorded deviations; see below.
    deviation = mic == 0.125
  )

pta_chk |>
  dplyr::select(regimen, matrix, mic, published, model, abs_diff, deviation) |>
  dplyr::rename(
    "Regimen"                 = regimen,
    "Matrix"                  = matrix,
    "MIC (mg/L)"              = mic,
    "Published %PTA"          = published,
    "Model %PTA"              = model,
    "Absolute difference (pp)" = abs_diff,
    "Recorded deviation"      = deviation
  ) |>
  knitr::kable(digits = c(0, 0, 4, 2, 2, 2, 0),
               caption = "Table 3 of Paiboonvong 2025 against the model.")
```

| Regimen | Matrix | MIC (mg/L) | Published %PTA | Model %PTA | Absolute difference (pp) | Recorded deviation |
|:---|:---|---:|---:|---:|---:|:---|
| 100 mg q12h | plasma | 0.0320 | 100.00 | 99.95 | 0.05 | FALSE |
| 100 mg q12h | plasma | 0.0625 | 100.00 | 99.45 | 0.55 | FALSE |
| 100 mg q12h | plasma | 0.1250 | 94.60 | 96.33 | 1.73 | TRUE |
| 100 mg q12h | plasma | 0.2500 | 85.72 | 85.73 | 0.01 | FALSE |
| 100 mg q12h | ELF | 0.0320 | 99.18 | 99.21 | 0.03 | FALSE |
| 100 mg q12h | ELF | 0.0625 | 95.36 | 95.67 | 0.31 | FALSE |
| 100 mg q12h | ELF | 0.1250 | 79.44 | 84.43 | 4.99 | TRUE |
| 100 mg q12h | ELF | 0.2500 | 62.50 | 62.42 | 0.08 | FALSE |
| 50 mg q12h | plasma | 0.0320 | 100.00 | 99.41 | 0.59 | FALSE |
| 50 mg q12h | plasma | 0.0625 | 96.34 | 96.33 | 0.01 | FALSE |
| 50 mg q12h | plasma | 0.1250 | 81.42 | 85.73 | 4.31 | TRUE |
| 50 mg q12h | plasma | 0.2500 | 63.48 | 63.89 | 0.41 | FALSE |
| 50 mg q12h | ELF | 0.0320 | 95.34 | 95.44 | 0.10 | FALSE |
| 50 mg q12h | ELF | 0.0625 | 84.84 | 84.43 | 0.41 | FALSE |
| 50 mg q12h | ELF | 0.1250 | 55.70 | 62.42 | 6.72 | TRUE |
| 50 mg q12h | ELF | 0.2500 | 33.70 | 34.53 | 0.83 | FALSE |

Table 3 of Paiboonvong 2025 against the model. {.table}

``` r


# Deterministic. Realised maximum over the gated cells is 0.83 pp.
stopifnot(max(pta_chk$abs_diff[!pta_chk$deviation]) < 2)
```

Twelve of the sixteen cells reproduce to within 0.9 percentage points.
The four that do not are exactly the four `MIC = 0.125 mg/L` cells, and
Table 3 is internally inconsistent on those cells regardless of any
model: because the attainment probability for `fAUC/MIC > 30` depends
only on the ratio of dose to MIC, halving the dose and halving the MIC
must give the identical value.

``` r

pairs_chk <- tibble::tribble(
  ~matrix,  ~cell_a,                   ~a,     ~cell_b,                  ~b,
  "plasma", "100 mg q12h, MIC 0.125",  94.60,  "50 mg q12h, MIC 0.0625", 96.34,
  "plasma", "100 mg q12h, MIC 0.25",   85.72,  "50 mg q12h, MIC 0.125",  81.42,
  "ELF",    "100 mg q12h, MIC 0.125",  79.44,  "50 mg q12h, MIC 0.0625", 84.84,
  "ELF",    "100 mg q12h, MIC 0.25",   62.50,  "50 mg q12h, MIC 0.125",  55.70
) |>
  dplyr::mutate(published_gap = abs(a - b))

pairs_chk |>
  dplyr::rename(
    "Matrix"                     = matrix,
    "Cell A"                     = cell_a, "%PTA A" = a,
    "Cell B (identical dose/MIC)" = cell_b, "%PTA B" = b,
    "Published gap (pp)"         = published_gap
  ) |>
  knitr::kable(digits = 2, caption = paste(
    "Pairs of Table 3 cells that carry the same dose-to-MIC ratio and must",
    "therefore be equal. Every discrepant pair involves a MIC = 0.125 cell."
  ))
```

| Matrix | Cell A | %PTA A | Cell B (identical dose/MIC) | %PTA B | Published gap (pp) |
|:---|:---|---:|:---|---:|---:|
| plasma | 100 mg q12h, MIC 0.125 | 94.60 | 50 mg q12h, MIC 0.0625 | 96.34 | 1.74 |
| plasma | 100 mg q12h, MIC 0.25 | 85.72 | 50 mg q12h, MIC 0.125 | 81.42 | 4.30 |
| ELF | 100 mg q12h, MIC 0.125 | 79.44 | 50 mg q12h, MIC 0.0625 | 84.84 | 5.40 |
| ELF | 100 mg q12h, MIC 0.25 | 62.50 | 50 mg q12h, MIC 0.125 | 55.70 | 6.80 |

Pairs of Table 3 cells that carry the same dose-to-MIC ratio and must
therefore be equal. Every discrepant pair involves a MIC = 0.125 cell.
{.table}

``` r


# Every published pair whose gap is large contains a MIC = 0.125 cell; the one
# pair that contains none of them (not shown, MIC 0.0625 vs 0.032 has no
# matching partner) cannot be formed from the four tabulated MICs. This
# assertion pins the observation that the discrepancy is confined to MIC 0.125.
stopifnot(all(grepl("0.125", paste(pairs_chk$cell_a, pairs_chk$cell_b))))
```

## Assumptions and deviations

- **The unbound fraction sits inside the ELF distribution rate, so the
  published partition coefficient is not the ELF:total-plasma ratio.**
  The abstract states that “the partition coefficient which relates drug
  exposure in ELF to drug exposure in plasma was estimated to be 0.77”,
  and the Discussion adds that “approximately 80% of the plasma
  concentration of sitafloxacin can be diffused into ELF”. The
  structural diagram in Figure 1 of the paper labels the central-to-ELF
  arrow `(Q_ELF x F_u x PC)/V` against a plain `Q_ELF/V_ELF` in the
  return direction, and the control stream codes the same thing as
  `K24 = QELF*0.63*PC/VBD` against `K42 = QELF/VELF`, so the
  concentration ratio the model actually produces is `fu * PC` = 0.63 \*
  0.772 = 0.486. The model file reproduces the control stream verbatim;
  the two prose statements should be read as the ELF-to-*unbound*-plasma
  ratio. The published median exposures corroborate the code rather than
  the prose: the reported ELF-to-plasma AUC ratios are 10.17/21.02 =
  0.484 and 5.12/10.42 = 0.491, both at 0.486 and neither at 0.772.
- **The values Results 3.4 labels “fAUC” are not multiplied by the
  unbound fraction on the plasma side.** Reproducing them requires
  taking the total plasma AUC (Deterministic reproduction section,
  agreement within 1%); multiplying by `fu = 0.63` would put the plasma
  value at 13.3 rather than the published 21.02. The ELF values do carry
  `fu`, but only through the `K24` micro-constant described above, not
  as a second multiplication. The same reading reproduces twelve of the
  sixteen `%PTA` cells of Table 3 to within 0.9 percentage points.
- **Recorded deviation: the four MIC = 0.125 mg/L cells of Table 3.**
  The model gives 96.33 / 85.73 (plasma, 100 and 50 mg q12h) and 84.43 /
  62.42 (ELF) where the paper reports 94.60 / 81.42 and 79.44 / 55.70.
  As shown above, Table 3’s own MIC = 0.125 entries contradict its MIC =
  0.0625 and MIC = 0.25 entries under the exact dose-to-MIC scaling that
  the target `fAUC/MIC > 30` implies, so the discrepancy is in the
  published row rather than in the transcription. These cells are
  excluded from the assertion and left visible in the table.
- **Simulation covariate distributions.** Methods 2.4 states the ages
  simulated (30, 40, 50, 60 and 70 years) but not the body weights.
  Holding weight at the 52 kg allometric reference reproduces all four
  published median exposures to within about 1%, so that is what the
  steady-state cohorts use. The study cohort’s weight and age
  distributions were constructed to match the medians and IQRs of
  Results 3.1; the paper does not publish the underlying distributions.
- **The age-on-F relationship must not be extrapolated below
  adulthood.** It is linear rather than exponential, so
  `F = 1 + 0.0258 * (AGE - 57)` reaches zero at 18.2 years and is
  negative below it. The virtual cohorts here are truncated at 25 years.
  The model file records the same caution in the `AGE` covariate
  metadata.
- **No inter-individual variability on V/F, Q/F ELF or V/F ELF.** The
  paper states that the data did not support estimating them (Results
  3.3 and the Discussion’s limitations paragraph) and the control stream
  fixes the corresponding `$OMEGA` elements to 0. They are therefore
  omitted from the model rather than carried as zero-variance etas.
- **Final estimates come from the simulation control stream.** Data S1
  contains two streams. Run56 (estimation) carries earlier initial
  values for the ELF inter-compartmental clearance (0.0714 vs 0.0441),
  the partition coefficient (0.733 vs 0.772), clearance (7.01 vs 7.03),
  central volume (115 vs 116), mean transit time (1.49 vs 1.48) and the
  ELF residual variance (0.274 vs 0.249). Run59 (simulation) carries
  values that match Table 2 exactly, and it is Run59 that the model file
  follows.
- **Residual variability is encoded as `lnorm()`.** The paper’s
  “additive residual error model on the logarithmic scale” is a
  log-normal residual on the linear concentration scale;
  `Cc ~ lnorm(expSd)` is nlmixr2’s spelling of the control stream’s
  `Y = LOG(IPRED) + EPS(1)`.
- **Related sitafloxacin models.** `modellib("Wu_2025_sitafloxacin")`
  and `modellib("Rodjun_2023_sitafloxacin")` package independent
  sitafloxacin popPK models fitted in non-critically-ill populations;
  neither carries an ELF compartment.
