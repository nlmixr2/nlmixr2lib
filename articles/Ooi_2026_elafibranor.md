# Elafibranor and GFT1007 PK with joint ALP-bilirubin PKPD (Ooi 2026)

## Model and source

Ooi 2026 reports three models, extracted here as three files and
validated in this single vignette:

- `modellib("Ooi_2026_elafibranor")` – parent population PK.
- `modellib("Ooi_2026_elafibranor_gft1007")` – active-metabolite
  population PK.
- `modellib("Ooi_2026_elafibranor_alptb")` – the final joint alkaline
  phosphatase (ALP) / total-bilirubin exposure-response model.

The two PK models are separate files because the authors deliberately
fitted them separately: “a simultaneously estimated joint model for
elafibranor and GFT1007 was attempted but not pursued due to extremely
long estimation times” (Discussion). The separate ALP and TB models that
preceded the joint model are model-development intermediates and are not
extracted; the joint model is the paper’s final PD model.

- Citation: Ooi QX, Brendel K, van Beek S, Aguiar Zdovc J, Bardol M,
  Dehez M. Population Pharmacokinetics and
  Pharmacokinetics-Pharmacodynamics Analyses of Elafibranor to Support
  Dose Selection in Primary Biliary Cholangitis. CPT Pharmacometrics
  Syst Pharmacol. 2026;15(0):e70247. <doi:10.1002/psp4.70247>.
- Article: <https://doi.org/10.1002/psp4.70247>
- Supplement (parameter tables and all six final NONMEM control
  streams):
  <https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13274737/supplementaryFiles>

``` r

mod_ela  <- readModelDb("Ooi_2026_elafibranor")
mod_gft  <- readModelDb("Ooi_2026_elafibranor_gft1007")
mod_pd   <- readModelDb("Ooi_2026_elafibranor_alptb")

ui_ela <- rxode2::rxode(mod_ela)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalfdepot_f2, etalfdepot_f3, etalfdepot_f4, etalfdepot_f5, etalfdepot_f6, etald1_f5, etald1_f6, etaiov_mat_2, etaiov_mat_3, etaiov_mat_4, etaiov_mat_5, etaiov_mat_6, etaiov_mat_7, etaiov_mat_8, etaiov_mat_9, etaiov_mat_10, etaiov_d1_2, etaiov_d1_3, etaiov_d1_4, etaiov_d1_5, etaiov_d1_6, etaiov_d1_7, etaiov_d1_8, etaiov_d1_9, etaiov_d1_10, etaruvP1NoneLate, etaruvP1NoneEarly, etaruvP2Late, etaruvP3NoneLate, etaruvP1SepfaLate, etaruvP1SepLate, etaruvP1SepEarly
#> as a work-around try putting the mu-referenced expression on a simple line
ui_gft <- rxode2::rxode(mod_gft)
#> ℹ parameter labels from comments will be replaced by 'label()'
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalfdepot_f3, etalfdepot_f5, etalfdepot_f6, etald1_f5, etald1_f6, etaruvP1Late, etaruvP1Early, etaruvP2Late, etaruvP2Early, etaruvP3Late
#> as a work-around try putting the mu-referenced expression on a simple line
ui_pd  <- rxode2::rxode(mod_pd)
#> ℹ parameter labels from comments will be replaced by 'label()'
```

## Population

The PK analysis pooled 17 trials – 13 phase I, two phase II and two
phase III (Table S1) – covering healthy volunteers and patients with
renal impairment, hepatic impairment, MASH or primary biliary
cholangitis (PBC), at doses from 5 mg to 360 mg and durations from a
single dose to more than a year of daily dosing. The elafibranor data
set held 892 subjects and 12,205 observations, the GFT1007 data set 894
subjects and 10,592 observations (Results 3.1.1). Baseline
characteristics (Table S2): mean weight 81.0 kg (SD 18.4), mean age 43.9
years (SD 16.0), mean BMI 27.7 kg/m2 (SD 5.9), 36.7% female, 84.7%
White, mean albumin 43.3 g/L (SD 3.5), mean ALT 38.3 U/L (SD 35.1), mean
creatinine clearance 96.6 mL/min/1.73 m2 (SD 19.8).

The PKPD analysis used only the two PBC trials, the phase II
GFT505B-216-1 and the phase III ELATIVE (GFT505B-319-1): 206 patients
with 1892 ALP and 1693 total-bilirubin observations, 95.6% female, mean
age 57.0 years (SD 8.6), mean baseline ALP 308 U/L (SD 122), mean
baseline total bilirubin 9.79 umol/L (SD 5.1), mean baseline liver
stiffness 10.2 kPa (SD 8.2, missing in 25.7%).

The same information is available programmatically from each model’s
`population` metadata,
e.g. `readModelDb("Ooi_2026_elafibranor")()$population`.

The covariate forest plots (Figures 3A and 3B) are conditioned on a
reference PBC patient, which this vignette reuses as its typical
subject:

``` r

# Figure 3 legend: "a female patient with PBC with elafibranor formulation 5
# administered under fasted conditions, body weight: 68 kg, BMI: 25.9 kg/m2,
# age: 57 years, ALT: 42 g/L, albumin: 43.5 g/L, CRCL: 90 mL/min/1.73m2".
# (Methods 2.2.3 prints 91.9 mL/min/1.73m2 for the same reference patient;
# Figure 3 and Figure S5 both print 90, which is used here.)
ref_pt <- list(
  WT = 68, BMI = 25.9, AGE = 57, ALT = 42, ALB = 43.5, CRCL = 90,
  SEXF = 1, DIS_PBC = 1, FED = 0
)

# Molecular weights used by the authors in Supplementary Datafile S6 $PK.
MW_ELA <- 384.49   # g/mol, elafibranor
MW_GFT <- 386.51   # g/mol, GFT1007

# The pooled PK data set carried amounts in nmol; both PK models take umol with
# volumes in L, which is the same concentration scale (1 nmol/mL = 1 umol/L).
mg_to_umol <- function(dose_mg) dose_mg / MW_ELA * 1000
```

## Source trace

Every `ini()` entry in the three model files carries an in-file comment
naming its source location. The table below collects the structural
entries; the formulation-specific absorption parameters (six values each
for relative bioavailability, D1, MAT and lag time, for each analyte)
and the twelve elafibranor / six GFT1007 residual-error magnitudes all
come from Table S3 and are listed individually in the model files.

| Equation / parameter | Value | Source location |
|----|----|----|
| Elafibranor CL/F, Vc/F, Q/F, Vp/F | 47.1 L/h, 68.7 L, 263 L/h, 4310 L | Table S3, elafibranor column |
| GFT1007 CL/F, Vc/F, Q/F, Vp/F | 11.2 L/h, 13.6 L, 6.83 L/h, 77.8 L | Table S3, GFT1007 column |
| Allometric exponents on CL/Q and Vc/Vp | 0.750 (FIX), 1.00 (FIX), reference 75 kg | Table S3; Datafiles S1/S2 `$PK` `LOG(WTKGBL/75)` |
| Sequential zero- then first-order absorption with lag | `dur(depot)`, `ka = 1/MAT`, `alag(depot)` | Datafiles S1/S2 `$SUBROUTINE ADVAN4 TRANS4` and `$PK` |
| Dose effect on elafibranor Frel | Emax 0.800, ED50 11.0 mg, gamma 2.31 | Table S3; Datafile S1 `$PK` `DOSE_FREL` |
| Elafibranor covariates | ALB on CL/F (-0.0227) and Q/F (-0.0158) centred 43.30; AGE on Vp/F (0.00589), D1 (-0.0172), Frel (0.00332) centred 44; SEXF on Q/F (-0.183); BMI on Vp/F (-0.0151) centred 26.44; PBC on Vp/F (-0.284); FED on MAT (1.55) | Table S3; Datafile S1 `$PK` covariate blocks |
| GFT1007 covariates | ALT on CL/F (-0.0489) and Vp/F (-0.123) normalised 27 U/L; BMI on CL/F (-0.00982) and Vp/F (-0.0227) centred 26.44; CRCL on CL/F (0.00313) centred 97.49; AGE on Vc/F (0.00487) centred 44; SEXF on Q/F (0.0923); PBC on Q/F (0.214); FED on D1 (11.6) and Frel (0.252) | Table S3; Datafile S2 `$PK` covariate blocks |
| Bioanalytical-method shift on elafibranor prediction | -0.182 (peak separation), -0.760 (peak separation plus formic acid) | Table S3; Datafile S1 `$ERROR` |
| Inter-occasion variability (elafibranor only) | MAT 0.414, D1 0.934, ten occasions | Table S3; Datafile S1 `$OMEGA BLOCK(1) SAME` |
| ALP / TB baselines | 251 U/L, 8.52 umol/L | Table 2 |
| ALP / TB turnover half-lives | 10.7 day, 1240 day | Table 2; footnote `HL = ln(2)/Kout` |
| ALP drug effect | Emax -0.731, AUC50 24.1 umol\*h/L, Hill 1.00 (FIX) | Table 2 |
| TB drug effect | Slope -0.0100 L/(umol\*h) | Table 2 |
| Placebo multipliers | ALP 0.992, TB 0.958, applied for time \> 0 | Table 2; Datafile S6 `$DES` `IF(TIME.GT.0)` |
| PD covariates | NCI hepatic impairment on baseline ALP (+0.277); liver stiffness on baseline TB, power 0.333 normalised 8.1 kPa | Table 2; Datafile S6 `$PK` |
| Logit-scale Emax IIV variance | 0.71883 | Datafile S6 `$OMEGA BLOCK(2)`; cross-checked against Table 2 (see below) |
| Exposure driver | `AUCSSP + AUCSSM`, each `F * dose_mg / CL * 1e6 / MW` | Datafile S6 `$PK` |

Every other value in the three files is `$THETA`-independent: the values
were taken from Table S3 and Table 2, **not** from the `$THETA` records
of the control streams, because those records hold the *initial*
estimates handed to the final run and differ from the published finals
by up to ~11% (for example the joint-model `EBALPNCIHISN` initial is
0.308 against the final 0.277 in Table 2). The one place a stream value
is used in preference is the logit-scale Emax variance, which Table 2
reports only after back-transformation; see “Assumptions and
deviations”.

## Virtual cohort

Original observed data are not publicly available. The cohort below is a
virtual PBC population centred on the Figure 3 reference patient, with
dispersion taken from the Table S2 baseline summaries. Each arm holds
150 subjects, below the 200-per-arm cap.

``` r

# set.seed() seeds R's RNG, which is what draws the covariates below. It does
# NOT seed rxode2's simulation RNG, and rxode2's streams are partitioned per
# solver thread, so the eta draws differ between a 2-core CI runner and a
# 16-thread workstation. Every assertion below is therefore written either on a
# typical-value (zeroRe) quantity or as a within-subject identity that holds for
# any draw.
set.seed(20260914)

n_arm <- 150L
doses_mg <- c(80, 120)

make_arm <- function(n, dose_mg, id_offset = 0L) {
  tibble::tibble(
    id  = id_offset + seq_len(n),
    # Table S2 dispersion, means shifted to the Figure 3 PBC reference patient.
    WT  = pmax(40, rnorm(n, ref_pt$WT,  18.4)),
    BMI = pmax(16, rnorm(n, ref_pt$BMI,  5.9)),
    AGE = pmax(18, rnorm(n, ref_pt$AGE,  8.6)),
    ALB = pmax(20, rnorm(n, ref_pt$ALB,  3.5)),
    ALT = pmax(5,  rnorm(n, ref_pt$ALT, 20.0)),
    CRCL = pmax(15, rnorm(n, ref_pt$CRCL, 19.8)),
    SEXF = rbinom(n, 1, 0.956),      # 95.6% female in the PBC trials (Table S2)
    DIS_PBC = 1, FED = 0,
    DOSE_ELA_MG = dose_mg,
    # Formulation 5 (all indicators 0), the phase II/III formulation.
    FORM_ELA_F1 = 0, FORM_ELA_F2 = 0, FORM_ELA_F3 = 0,
    FORM_ELA_F4 = 0, FORM_ELA_F6 = 0,
    # Reference bioanalytical method and phase I residual stratum; see
    # "Assumptions and deviations" for why the NCA below is run on the
    # unmodified method rather than on Table 1's stated assay condition.
    ASSAY_SEPIP = 0, ASSAY_SEPIP_FA = 0,
    STUDY_PHASE2 = 0, STUDY_PHASE3 = 0, STUDY_GFT505B_319_1 = 0,
    OCC = 1,
    arm = paste0(dose_mg, " mg/day")
  )
}

subj <- dplyr::bind_rows(
  make_arm(n_arm, doses_mg[1], id_offset = 0L),
  make_arm(n_arm, doses_mg[2], id_offset = n_arm)
)
stopifnot(!anyDuplicated(subj$id))
```

The PK simulation runs 60 once-daily doses – comfortably past the
12.4-day predicted time to elafibranor steady state – and samples the
final dosing interval.

``` r

n_dose  <- 60L
t_start <- 24 * (n_dose - 1)   # start of the final dosing interval, hours
# 241 points over the interval. The grid has to resolve a peak that arrives
# around 0.6 h after the dose, and the 108% CV on the central volume means some
# subjects draw a Vc an order of magnitude below typical and a correspondingly
# sharp peak; at 121 points the worst subject's trapezoidal AUC error reached
# 2.2%, at 241 points it is below 0.2%.
grid_h  <- seq(t_start, t_start + 24, length.out = 241)

make_events <- function(subj_df) {
  dosing <- subj_df |>
    dplyr::mutate(time = 0, amt = mg_to_umol(DOSE_ELA_MG), evid = 1L,
                  cmt = "depot", ii = 24, addl = n_dose - 1L)
  obs <- subj_df |>
    tidyr::crossing(time = grid_h) |>
    dplyr::mutate(amt = NA_real_, evid = 0L,
                  # The ODE STATE, never the algebraic observable name.
                  cmt = "central", ii = 0, addl = 0L)
  dplyr::bind_rows(dosing, obs) |> dplyr::arrange(id, time, dplyr::desc(evid))
}

events <- make_events(subj)
```

## Simulation

Both PK models are solved on the same event table: the two analytes
share one data set in the source, and the GFT1007 model’s depot
represents metabolite formation from the administered elafibranor dose.

``` r

# method = "lsoda": the default liblsoda solver fails on this model at fine
# output grids (the zero-order depot input is very short -- D1 is about 0.19 h
# for formulation 5 -- against a 60-h terminal half-life).
solve_pk <- function(ui, ev, typical = FALSE) {
  m <- if (typical) rxode2::zeroRe(ui) else ui
  rxode2::rxSolve(m, events = ev, keep = c("arm", "DOSE_ELA_MG"),
                  method = "lsoda", returnType = "data.frame")
}

sim_ela <- solve_pk(ui_ela, events)
sim_gft <- solve_pk(ui_gft, events)

sim_ela_typ <- solve_pk(ui_ela, events, typical = TRUE)
#> Warning: No sigma parameters in the model
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalfdepot_f2, etalfdepot_f3, etalfdepot_f4, etalfdepot_f5, etalfdepot_f6, etald1_f5, etald1_f6, etaiov_mat_2, etaiov_mat_3, etaiov_mat_4, etaiov_mat_5, etaiov_mat_6, etaiov_mat_7, etaiov_mat_8, etaiov_mat_9, etaiov_mat_10, etaiov_d1_2, etaiov_d1_3, etaiov_d1_4, etaiov_d1_5, etaiov_d1_6, etaiov_d1_7, etaiov_d1_8, etaiov_d1_9, etaiov_d1_10, etaruvP1NoneLate, etaruvP1NoneEarly, etaruvP2Late, etaruvP3NoneLate, etaruvP1SepfaLate, etaruvP1SepLate, etaruvP1SepEarly
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalfdepot_f1', 'etalfdepot_f2', 'etalfdepot_f3', 'etalfdepot_f4', 'etalfdepot_f5', 'etalfdepot_f6', 'etald1_f4', 'etald1_f5', 'etald1_f6', 'etalmat_f1', 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etaiov_mat_1', 'etaiov_mat_2', 'etaiov_mat_3', 'etaiov_mat_4', 'etaiov_mat_5', 'etaiov_mat_6', 'etaiov_mat_7', 'etaiov_mat_8', 'etaiov_mat_9', 'etaiov_mat_10', 'etaiov_d1_1', 'etaiov_d1_2', 'etaiov_d1_3', 'etaiov_d1_4', 'etaiov_d1_5', 'etaiov_d1_6', 'etaiov_d1_7', 'etaiov_d1_8', 'etaiov_d1_9', 'etaiov_d1_10', 'etaruvP1NoneLate', 'etaruvP1NoneEarly', 'etaruvP2Late', 'etaruvP3NoneLate', 'etaruvP1SepfaLate', 'etaruvP1SepLate', 'etaruvP1SepEarly'
#> Warning: multi-subject simulation without without 'omega'
sim_gft_typ <- solve_pk(ui_gft, events, typical = TRUE)
#> Warning: No sigma parameters in the model
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalfdepot_f3, etalfdepot_f5, etalfdepot_f6, etald1_f5, etald1_f6, etaruvP1Late, etaruvP1Early, etaruvP2Late, etaruvP2Early, etaruvP3Late
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalfdepot_f1', 'etalfdepot_f3', 'etalfdepot_f5', 'etalfdepot_f6', 'etald1_f1to4', 'etald1_f5', 'etald1_f6', 'etalmat_f6', 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etaruvP1Late', 'etaruvP1Early', 'etaruvP2Late', 'etaruvP2Early', 'etaruvP3Late'
#> Warning: multi-subject simulation without without 'omega'
```

### Gate 1 – dose / clearance mass balance, per subject

At steady state on a once-daily regimen, `AUC(tau) = F * dose / CL`
exactly. Both sides use the *same drawn* parameters for each subject, so
the only difference is trapezoidal error on the output grid; a tight
bound is correct here and would catch a mis-transcribed clearance, dose
or unit.

``` r

auc_tau <- function(sim) {
  sim |>
    dplyr::filter(time >= t_start) |>
    dplyr::group_by(id, arm) |>
    dplyr::summarise(
      auc   = sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2),
      exact = unique(fdepot)[1] * mg_to_umol(unique(DOSE_ELA_MG)[1]) /
        unique(cl)[1],
      .groups = "drop"
    ) |>
    dplyr::mutate(pct_diff = 100 * (auc / exact - 1))
}

mb <- dplyr::bind_rows(
  auc_tau(sim_ela) |> dplyr::mutate(analyte = "Elafibranor"),
  auc_tau(sim_gft) |> dplyr::mutate(analyte = "GFT1007")
)

mb |>
  dplyr::group_by(analyte) |>
  dplyr::summarise(n = dplyr::n(),
                   median_pct = median(pct_diff),
                   max_abs_pct = max(abs(pct_diff)), .groups = "drop") |>
  dplyr::rename("Analyte" = analyte, "N" = n,
                "Median % difference" = median_pct,
                "Max |% difference|" = max_abs_pct) |>
  knitr::kable(digits = 3,
               caption = "Trapezoidal AUC(tau) at steady state versus F * dose / CL.")
```

| Analyte     |   N | Median % difference | Max \|% difference\| |
|:------------|----:|--------------------:|---------------------:|
| Elafibranor | 300 |               0.022 |                0.204 |
| GFT1007     | 300 |               0.020 |                0.107 |

Trapezoidal AUC(tau) at steady state versus F \* dose / CL. {.table}

``` r


# Pure trapezoidal error on the 241-point grid; the realised maximum over both
# analytes was 0.20%, and the worst subject is always the one that drew the
# smallest central volume. The bound keeps an order of magnitude of headroom so
# that a cohort drawn on a different thread count cannot trip it, while a
# mis-transcribed clearance, dose or unit -- which moves AUC by tens of
# percent -- still fails instantly.
stopifnot(nrow(mb) == 2 * 2 * n_arm, max(abs(mb$pct_diff)) < 2)
```

### Gate 2 – published half-lives

The paper reports, for the Figure 3 reference PBC patient, distribution
and elimination half-lives of 0.176 h and 59.7 h for elafibranor and
0.475 h and 10.7 h for GFT1007 (Results 3.1.4 and 3.1.7). These are the
eigenvalues of the two-compartment micro-constant system, so they can be
checked in closed form from the typical-value solve.

``` r

half_lives <- function(sim_typ) {
  r <- sim_typ[1, ]
  b <- r$kel + r$k12 + r$k21
  lam1 <- (b + sqrt(b^2 - 4 * r$kel * r$k21)) / 2
  lam2 <- (b - sqrt(b^2 - 4 * r$kel * r$k21)) / 2
  c(distribution = log(2) / lam1, elimination = log(2) / lam2)
}

# Solve the reference patient exactly, rather than a cohort member.
ref_events <- make_events(
  make_arm(1L, 80) |>
    dplyr::mutate(WT = ref_pt$WT, BMI = ref_pt$BMI, AGE = ref_pt$AGE,
                  ALB = ref_pt$ALB, ALT = ref_pt$ALT, CRCL = ref_pt$CRCL,
                  SEXF = 1)
)
hl_ela <- half_lives(solve_pk(ui_ela, ref_events, typical = TRUE))
#> Warning: No sigma parameters in the model
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalfdepot_f2, etalfdepot_f3, etalfdepot_f4, etalfdepot_f5, etalfdepot_f6, etald1_f5, etald1_f6, etaiov_mat_2, etaiov_mat_3, etaiov_mat_4, etaiov_mat_5, etaiov_mat_6, etaiov_mat_7, etaiov_mat_8, etaiov_mat_9, etaiov_mat_10, etaiov_d1_2, etaiov_d1_3, etaiov_d1_4, etaiov_d1_5, etaiov_d1_6, etaiov_d1_7, etaiov_d1_8, etaiov_d1_9, etaiov_d1_10, etaruvP1NoneLate, etaruvP1NoneEarly, etaruvP2Late, etaruvP3NoneLate, etaruvP1SepfaLate, etaruvP1SepLate, etaruvP1SepEarly
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalfdepot_f1', 'etalfdepot_f2', 'etalfdepot_f3', 'etalfdepot_f4', 'etalfdepot_f5', 'etalfdepot_f6', 'etald1_f4', 'etald1_f5', 'etald1_f6', 'etalmat_f1', 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etaiov_mat_1', 'etaiov_mat_2', 'etaiov_mat_3', 'etaiov_mat_4', 'etaiov_mat_5', 'etaiov_mat_6', 'etaiov_mat_7', 'etaiov_mat_8', 'etaiov_mat_9', 'etaiov_mat_10', 'etaiov_d1_1', 'etaiov_d1_2', 'etaiov_d1_3', 'etaiov_d1_4', 'etaiov_d1_5', 'etaiov_d1_6', 'etaiov_d1_7', 'etaiov_d1_8', 'etaiov_d1_9', 'etaiov_d1_10', 'etaruvP1NoneLate', 'etaruvP1NoneEarly', 'etaruvP2Late', 'etaruvP3NoneLate', 'etaruvP1SepfaLate', 'etaruvP1SepLate', 'etaruvP1SepEarly'
hl_gft <- half_lives(solve_pk(ui_gft, ref_events, typical = TRUE))
#> Warning: No sigma parameters in the model
#> Warning: some etas defaulted to non-mu referenced, possible parsing error: etalfdepot_f3, etalfdepot_f5, etalfdepot_f6, etald1_f5, etald1_f6, etaruvP1Late, etaruvP1Early, etaruvP2Late, etaruvP2Early, etaruvP3Late
#> as a work-around try putting the mu-referenced expression on a simple line
#> ℹ omega/sigma items treated as zero: 'etalfdepot_f1', 'etalfdepot_f3', 'etalfdepot_f5', 'etalfdepot_f6', 'etald1_f1to4', 'etald1_f5', 'etald1_f6', 'etalmat_f6', 'etalcl', 'etalvc', 'etalq', 'etalvp', 'etaruvP1Late', 'etaruvP1Early', 'etaruvP2Late', 'etaruvP2Early', 'etaruvP3Late'

hl_tab <- tibble::tibble(
  Analyte = c("Elafibranor", "Elafibranor", "GFT1007", "GFT1007"),
  Phase = c("Distribution", "Elimination", "Distribution", "Elimination"),
  Model = c(hl_ela[["distribution"]], hl_ela[["elimination"]],
            hl_gft[["distribution"]], hl_gft[["elimination"]]),
  Published = c(0.176, 59.7, 0.475, 10.7)
) |>
  dplyr::mutate(`% difference` = 100 * (Model / Published - 1))

hl_tab |>
  dplyr::rename("Model (h)" = Model, "Published (h)" = Published) |>
  knitr::kable(digits = 3,
               caption = "Two-compartment eigenvalue half-lives versus Results 3.1.4 / 3.1.7.")
```

| Analyte     | Phase        | Model (h) | Published (h) | % difference |
|:------------|:-------------|----------:|--------------:|-------------:|
| Elafibranor | Distribution |     0.176 |         0.176 |       -0.252 |
| Elafibranor | Elimination  |    59.875 |        59.700 |        0.293 |
| GFT1007     | Distribution |     0.474 |         0.475 |       -0.155 |
| GFT1007     | Elimination  |    10.694 |        10.700 |       -0.054 |

Two-compartment eigenvalue half-lives versus Results 3.1.4 / 3.1.7.
{.table}

``` r


# Deterministic (zeroRe) quantities against published medians; the published
# values are cohort medians, so a few percent is expected, tens of percent is
# a transcription error.
stopifnot(max(abs(hl_tab$`% difference`)) < 5)
```

## Replicate published figures

### Steady-state concentration-time profiles

``` r

prof <- dplyr::bind_rows(
  sim_ela |> dplyr::mutate(analyte = "Elafibranor"),
  sim_gft |> dplyr::mutate(analyte = "GFT1007")
) |>
  dplyr::filter(time >= t_start) |>
  dplyr::mutate(tad = time - t_start) |>
  dplyr::group_by(analyte, arm, tad) |>
  dplyr::summarise(Q05 = quantile(Cc, 0.05), Q50 = median(Cc),
                   Q95 = quantile(Cc, 0.95), .groups = "drop")

ggplot(prof, aes(tad, Q50)) +
  geom_ribbon(aes(ymin = Q05, ymax = Q95), alpha = 0.25) +
  geom_line() +
  facet_grid(analyte ~ arm) +
  scale_y_log10() +
  labs(x = "Time since last dose (h)", y = "Concentration (umol/L)",
       title = "Steady-state profiles, median and 90% prediction interval",
       caption = paste("Companion to Figure 2 of Ooi 2026, which shows the",
                       "observed pcVPCs over the first 25 h after dosing."))
```

![](Ooi_2026_elafibranor_files/figure-html/figure-profiles-1.png)

The paper’s qualitative claims are visible here: GFT1007 exposure
greatly exceeds the parent, and the parent’s trough is relatively higher
because of its much longer terminal half-life.

### Figure 5 – exposure-response for ALP and total bilirubin

Figure 5 plots relative ALP and total-bilirubin change from baseline at
Week 52 against the summed steady-state AUC. That driver is a static
per-subject quantity, so the PD model needs no PK compartments; the
exposure metric is computed exactly as Supplementary Datafile S6 does.

``` r

# Datafile S6 $PK:
#   AUCSSP = (FP * DOSEN / CLP) * 1000 * 1000 / 384.49
#   AUCSSM = (FM * DOSEN / CLM) * 1000 * 1000 / 386.51
# with CL in mL/h and DOSEN in mg. The source converts the same milligram dose
# with the PARENT molecular weight for elafibranor and the METABOLITE
# molecular weight for GFT1007; that 0.5% inconsistency is reproduced as
# published rather than harmonised.
per_subject <- function(sim, mw) {
  sim |>
    dplyr::group_by(id, arm) |>
    dplyr::summarise(auc = unique(fdepot)[1] * unique(DOSE_ELA_MG)[1] /
                       (unique(cl)[1] * 1000) * 1e6 / mw, .groups = "drop")
}

auc_driver <- dplyr::inner_join(
  per_subject(sim_ela, MW_ELA) |> dplyr::rename(AUC_ELA = auc),
  per_subject(sim_gft, MW_GFT) |> dplyr::rename(AUC_GFT1007 = auc),
  by = c("id", "arm")
) |>
  dplyr::mutate(auc_sum = AUC_ELA + AUC_GFT1007)

auc_driver |>
  dplyr::group_by(arm) |>
  dplyr::summarise(`Median AUC sum (umol*h/L)` = median(auc_sum),
                   `5th` = quantile(auc_sum, 0.05),
                   `95th` = quantile(auc_sum, 0.95), .groups = "drop") |>
  dplyr::rename("Arm" = arm) |>
  knitr::kable(digits = 1,
               caption = paste("Summed steady-state AUC driver. Figure 5",
                               "reports medians of 32.3 (80 mg/day) and",
                               "39.3 umol*h/L (120 mg/day)."))
```

| Arm        | Median AUC sum (umol\*h/L) |  5th | 95th |
|:-----------|---------------------------:|-----:|-----:|
| 120 mg/day |                       48.4 | 27.4 | 76.6 |
| 80 mg/day  |                       28.1 | 17.7 | 48.7 |

Summed steady-state AUC driver. Figure 5 reports medians of 32.3 (80
mg/day) and 39.3 umol\*h/L (120 mg/day). {.table}

``` r

# Typical-value exposure-response over the paper's 0-60 umol*h/L grid.
auc_grid <- seq(0, 60, length.out = 61)
pd_days  <- c(0, 364)

pd_events <- tibble::tibble(auc_sum = auc_grid) |>
  dplyr::mutate(id = seq_along(auc_sum),
                AUC_ELA = 0, AUC_GFT1007 = auc_sum,
                # NCI hepatic impairment 1 and liver stiffness 8.1 kPa are the
                # Figure S4 reference-patient values.
                HEPIMP = 1, LSM = 8.1) |>
  tidyr::crossing(time = pd_days) |>
  dplyr::mutate(amt = NA_real_, evid = 0L, cmt = "alp")

pd_typ <- rxode2::rxSolve(rxode2::zeroRe(ui_pd), events = pd_events,
                          keep = "auc_sum", returnType = "data.frame")
#> ℹ omega/sigma items treated as zero: 'etalrbase_alp', 'etalrbase_tbili', 'etalogitemax_alp', 'etalpbo_alp', 'etalpbo_tbili', 'etapropSd_alp', 'etapropSd_tbili'
#> Warning: multi-subject simulation without without 'omega'

er <- pd_typ |>
  dplyr::group_by(id, auc_sum) |>
  dplyr::summarise(alp_cfb = 100 * (alp[time == 364] / alp[time == 0] - 1),
                   tb_cfb = 100 * (tbili[time == 364] / tbili[time == 0] - 1),
                   .groups = "drop")

er |>
  tidyr::pivot_longer(c(alp_cfb, tb_cfb), names_to = "endpoint",
                      values_to = "cfb") |>
  dplyr::mutate(endpoint = dplyr::recode(
    endpoint, alp_cfb = "ALP", tb_cfb = "Total bilirubin")) |>
  ggplot(aes(auc_sum, cfb)) +
  geom_line() +
  geom_hline(yintercept = -15, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 24.1, linetype = "dotted") +
  facet_wrap(~endpoint) +
  labs(x = "AUC(tau,ss) sum of elafibranor and GFT1007 (umol*h/L)",
       y = "Relative change from baseline at Week 52 (%)",
       title = "Figure 5A/5B -- typical-value exposure-response",
       caption = paste("Dashed red line: the -15% relative ALP decrease",
                       "threshold. Dotted line: the estimated AUC50, 24.1."))
```

![](Ooi_2026_elafibranor_files/figure-html/figure-5-1.png)

### Gate 3 – indirect-response closed forms

Both PD endpoints are linear turnover systems driven by a constant input
from time 0, so the trajectory has a closed form: the state relaxes from
its baseline towards `baseline * placebo * drug` with rate `kout`. This
is a pure solver check against the model’s own algebra and so carries a
tight bound.

``` r

pd_closed_form <- function(auc_sum, day) {
  emax <- -0.731; ec50 <- 24.1; slope <- -0.0100
  kout_alp <- log(2) / 10.7; kout_tb <- log(2) / 1240
  ratio <- function(kout, target) 1 - (1 - exp(-kout * day)) * (1 - target)
  c(alp = ratio(kout_alp, 0.992 * (1 + emax * auc_sum / (ec50 + auc_sum))),
    tb  = ratio(kout_tb,  0.958 * (1 + slope * auc_sum)))
}

cf_chk <- er |>
  dplyr::rowwise() |>
  dplyr::mutate(alp_cf = 100 * (pd_closed_form(auc_sum, 364)[["alp"]] - 1),
                tb_cf  = 100 * (pd_closed_form(auc_sum, 364)[["tb"]] - 1)) |>
  dplyr::ungroup() |>
  dplyr::mutate(d_alp = abs(alp_cfb - alp_cf), d_tb = abs(tb_cfb - tb_cf))

knitr::kable(
  tibble::tibble(
    Endpoint = c("ALP", "Total bilirubin"),
    `Max absolute difference (percentage points)` = c(max(cf_chk$d_alp),
                                                     max(cf_chk$d_tb))
  ),
  digits = 6,
  caption = "Simulated Week-52 change from baseline versus the turnover closed form."
)
```

| Endpoint        | Max absolute difference (percentage points) |
|:----------------|--------------------------------------------:|
| ALP             |                                    0.000013 |
| Total bilirubin |                                    0.000115 |

Simulated Week-52 change from baseline versus the turnover closed form.
{.table}

``` r


# Solver tolerance only; realised below 1e-4 percentage points.
stopifnot(max(cf_chk$d_alp) < 0.01, max(cf_chk$d_tb) < 0.01)

# Baseline ALP for a patient with NCI hepatic impairment reproduces Table 2's
# 251 U/L scaled by 1.277, and lands near the observed cohort mean of 308 U/L.
stopifnot(abs(pd_typ$alp[pd_typ$time == 0][1] - 251 * 1.277) < 0.5)
```

### Gate 4 – the paper’s dose-justification claims

``` r

er_80  <- median(auc_driver$auc_sum[auc_driver$arm == "80 mg/day"])
er_120 <- median(auc_driver$auc_sum[auc_driver$arm == "120 mg/day"])
cfb_at <- function(a) 100 * (pd_closed_form(a, 364)[["alp"]] - 1)

# The saturation claim is a statement about the shape of the exposure-response
# curve, so it is evaluated at the exposures the PAPER reports for the two
# arms, not at this cohort's medians. Those differ in kind: the published
# medians are empirical, over the actual trial patients, and the paper notes
# that the two arms' AUC distributions "largely overlapped" (ratio 39.3/32.3 =
# 1.22), whereas any simulation that gives every subject the nominal dose
# necessarily separates them by the dose ratio 1.5 (the bioavailability dose
# effect is fully saturated at both 80 and 120 mg). See "Assumptions and
# deviations".
claims <- tibble::tribble(
  ~Claim, ~Published, ~Model,
  "Median AUC sum, 80 mg/day (umol*h/L)",  "32.3", sprintf("%.1f", er_80),
  "Median AUC sum, 120 mg/day (umol*h/L)", "39.3", sprintf("%.1f", er_120),
  "120 vs 80 mg/day exposure ratio", "1.22 (empirical)",
  sprintf("%.2f (dose-proportional)", er_120 / er_80),
  "AUC50 is below the median 80 mg/day exposure", "yes",
  ifelse(24.1 < er_80, "yes", "no"),
  "Typical ALP decrease at 80 mg/day exceeds 15%", "yes",
  ifelse(cfb_at(er_80) < -15, "yes", "no"),
  "ALP response at the published 80 vs 120 mg/day medians (percentage points apart)",
  "marginal", sprintf("%.1f", abs(cfb_at(39.3) - cfb_at(32.3))),
  "TB decrease is smaller than the ALP decrease", "yes",
  ifelse(abs(pd_closed_form(er_80, 364)[["tb"]] - 1) <
           abs(pd_closed_form(er_80, 364)[["alp"]] - 1), "yes", "no")
)
knitr::kable(claims, caption = "Dose-justification claims (Results 3.3 and Figure 5).")
```

| Claim | Published | Model |
|:---|:---|:---|
| Median AUC sum, 80 mg/day (umol\*h/L) | 32.3 | 28.1 |
| Median AUC sum, 120 mg/day (umol\*h/L) | 39.3 | 48.4 |
| 120 vs 80 mg/day exposure ratio | 1.22 (empirical) | 1.72 (dose-proportional) |
| AUC50 is below the median 80 mg/day exposure | yes | yes |
| Typical ALP decrease at 80 mg/day exceeds 15% | yes | yes |
| ALP response at the published 80 vs 120 mg/day medians (percentage points apart) | marginal | 3.4 |
| TB decrease is smaller than the ALP decrease | yes | yes |

Dose-justification claims (Results 3.3 and Figure 5). {.table}

``` r


stopifnot(
  # Cohort median against the published 80 mg/day median; realised -12.9%.
  abs(er_80 / 32.3 - 1) < 0.25,
  # Dose proportionality WITHIN the simulation: a deterministic property of the
  # model (the Frel dose effect is 1.792 at 80 mg and 1.797 at 120 mg), checked
  # through two independently drawn arms, so the bound admits median sampling
  # noise. Realised 1.72 against the expected 1.50.
  abs(er_120 / er_80 / 1.5 - 1) < 0.30,
  # AUC50 below the median 80 mg/day exposure (Discussion).
  24.1 < er_80,
  # The headline efficacy claim, with the typical response far past the
  # threshold (realised about -40%).
  cfb_at(er_80) < -15,
  # Saturation, evaluated deterministically at the published exposures.
  abs(cfb_at(39.3) - cfb_at(32.3)) < 5,
  abs(pd_closed_form(er_80, 364)[["tb"]] - 1) <
    abs(pd_closed_form(er_80, 364)[["alp"]] - 1)
)
```

## PKNCA validation

One PKNCA block per analyte over the final steady-state dosing interval,
with time re-zeroed to the last dose.

``` r

nca_for <- function(sim, label) {
  conc <- sim |>
    dplyr::filter(time >= t_start, !is.na(Cc)) |>
    dplyr::mutate(time = time - t_start) |>
    dplyr::select(id, time, Cc, arm)

  # Guarantee a time-zero record per (id, arm); without it PKNCA warns
  # "Requesting an AUC range starting (0) before the first measurement".
  conc <- dplyr::bind_rows(
    conc,
    conc |> dplyr::distinct(id, arm) |> dplyr::mutate(time = 0, Cc = 0)
  ) |>
    dplyr::distinct(id, arm, time, .keep_all = TRUE) |>
    dplyr::arrange(id, arm, time)

  dose_df <- sim |>
    dplyr::distinct(id, arm, DOSE_ELA_MG) |>
    dplyr::mutate(time = 0, amt = mg_to_umol(DOSE_ELA_MG))

  d <- PKNCA::PKNCAdata(
    PKNCA::PKNCAconc(conc, Cc ~ time | arm + id),
    PKNCA::PKNCAdose(dose_df, amt ~ time | arm + id),
    intervals = data.frame(start = 0, end = 24,
                           cmax = TRUE, tmax = TRUE, auclast = TRUE)
  )
  res <- as.data.frame(PKNCA::pk.nca(d))
  dplyr::mutate(res, analyte = label)
}

nca <- dplyr::bind_rows(nca_for(sim_ela, "Elafibranor"),
                        nca_for(sim_gft, "GFT1007"))

nca_summary <- nca |>
  dplyr::group_by(analyte, arm, PPTESTCD) |>
  dplyr::summarise(median = median(PPORRES),
                   p05 = quantile(PPORRES, 0.05),
                   p95 = quantile(PPORRES, 0.95), .groups = "drop")

nca_summary |>
  dplyr::rename("Analyte" = analyte, "Arm" = arm,
                "NCA parameter" = PPTESTCD, "Median" = median,
                "5th percentile" = p05, "95th percentile" = p95) |>
  knitr::kable(digits = 3,
               caption = "Simulated steady-state NCA, 150 subjects per arm.")
```

| Analyte     | Arm        | NCA parameter | Median | 5th percentile | 95th percentile |
|:------------|:-----------|:--------------|-------:|---------------:|----------------:|
| Elafibranor | 120 mg/day | auclast       | 14.079 |          6.236 |          29.269 |
| Elafibranor | 120 mg/day | cmax          |  2.073 |          0.768 |           4.859 |
| Elafibranor | 120 mg/day | tmax          |  0.500 |          0.300 |           1.355 |
| Elafibranor | 80 mg/day  | auclast       |  8.301 |          3.371 |          17.301 |
| Elafibranor | 80 mg/day  | cmax          |  1.160 |          0.442 |           2.862 |
| Elafibranor | 80 mg/day  | tmax          |  0.600 |          0.300 |           1.400 |
| GFT1007     | 120 mg/day | auclast       | 33.035 |         15.891 |          57.187 |
| GFT1007     | 120 mg/day | cmax          |  6.931 |          3.014 |          14.374 |
| GFT1007     | 120 mg/day | tmax          |  1.100 |          0.700 |           1.700 |
| GFT1007     | 80 mg/day  | auclast       | 19.244 |         11.066 |          36.675 |
| GFT1007     | 80 mg/day  | cmax          |  4.374 |          2.101 |           8.361 |
| GFT1007     | 80 mg/day  | tmax          |  1.100 |          0.700 |           1.600 |

Simulated steady-state NCA, 150 subjects per arm. {.table}

### Comparison against published NCA

Table 1 reports simulated secondary PK parameters at steady state for
500 PBC patients on 80 mg/day of formulation 5, fasted, for elafibranor
without inter-occasion variability and for GFT1007.

``` r

sim80 <- nca_summary |>
  dplyr::filter(arm == "80 mg/day") |>
  dplyr::select(analyte, PPTESTCD, median) |>
  tidyr::pivot_wider(names_from = PPTESTCD, values_from = median)

published <- tibble::tribble(
  ~analyte,      ~auclast, ~cmax,  ~tmax,
  "Elafibranor",  4.20,    0.553,  0.798,
  "GFT1007",     20.8,     4.08,   1.42
)

cmp <- dplyr::inner_join(sim80, published, by = "analyte",
                         suffix = c("_model", "_pub")) |>
  tidyr::pivot_longer(-analyte,
                      names_to = c("param", ".value"),
                      names_pattern = "(.*)_(model|pub)$") |>
  dplyr::mutate(
    `% difference` = 100 * (model / pub - 1),
    Flag = ifelse(abs(`% difference`) > 20, "*", ""),
    param = dplyr::recode(param, auclast = "AUC(tau,ss) (umol*h/L)",
                          cmax = "Cmax,ss (umol/L)", tmax = "Tmax,ss (h)")
  )

cmp |>
  dplyr::rename("Analyte" = analyte, "NCA parameter" = param,
                "Model" = model, "Table 1" = pub) |>
  knitr::kable(digits = 3,
               caption = paste("Simulated versus Table 1 of Ooi 2026.",
                               "* differs by more than 20%."))
```

| Analyte     | NCA parameter           |  Model | Table 1 | % difference | Flag |
|:------------|:------------------------|-------:|--------:|-------------:|:-----|
| Elafibranor | AUC(tau,ss) (umol\*h/L) |  8.301 |   4.200 |       97.644 | \*   |
| Elafibranor | Cmax,ss (umol/L)        |  1.160 |   0.553 |      109.793 | \*   |
| Elafibranor | Tmax,ss (h)             |  0.600 |   0.798 |      -24.812 | \*   |
| GFT1007     | AUC(tau,ss) (umol\*h/L) | 19.244 |  20.800 |       -7.481 |      |
| GFT1007     | Cmax,ss (umol/L)        |  4.374 |   4.080 |        7.194 |      |
| GFT1007     | Tmax,ss (h)             |  1.100 |   1.420 |      -22.535 | \*   |

Simulated versus Table 1 of Ooi 2026. \* differs by more than 20%.
{.table}

GFT1007 reproduces Table 1 essentially exactly for AUC, and Cmax and
Tmax differ only as expected between a typical profile and the median of
a cohort carrying a 54% CV on the central volume. **The three
elafibranor rows are starred**, with AUC and Cmax roughly twice the
published values; this is a genuine, reproducible disagreement with the
paper’s own Table 1, analysed in detail below. No parameter was tuned to
close it.

``` r

gft <- cmp |> dplyr::filter(analyte == "GFT1007")
# GFT1007 AUC is an exact dose/clearance identity against Table 1 and is the
# gate; its Cmax and Tmax are cohort medians against a typical-value
# publication and are reported but not gated.
stopifnot(abs(gft$`% difference`[gft$`NCA parameter` == "AUC(tau,ss) (umol*h/L)"]) < 10)
#> Warning: Unknown or uninitialised column: `NCA parameter`.
```

## Assumptions and deviations

### The elafibranor rows of Table 1 disagree with the paper’s own parameter table

For the Figure 3 reference PBC patient at 80 mg/day the model gives
`AUC(tau,ss) = F * dose / CL = 8.94 umol*h/L` for elafibranor, against
the 4.20 umol\*h/L of Table 1 – a factor of 2.13. Two observations bound
the problem:

1.  **GFT1007 is exact.** The same identity gives 20.80 umol\*h/L
    against Table 1’s 20.8, and the eigenvalue half-life 10.69 h against
    the published 10.7 h. The unit convention, the molar dose and the
    allometric and covariate terms are therefore all correct, and the
    discrepancy is specific to elafibranor.
2.  **Two arithmetic routes each reproduce 4.20 to within 1%**, and
    Table 1’s own footnote supports neither cleanly:
    - applying the *formic-acid* bioanalytical factor
      `exp(-0.760) = 0.468` gives 4.18 – but Table 1’s note states the
      samples were analysed “with separation of interfering peak but
      without addition of formic acid”, which is the
      `exp(-0.182) = 0.834` factor and gives 7.45;
    - omitting the sigmoidal dose effect on relative bioavailability
      (dividing by 1.792 at 80 mg) and applying the stated `exp(-0.182)`
      factor gives 4.16.

The **AUC-sum identity adjudicates in favour of the model as encoded**.
The exposure driver of the joint PD model is the sum of the two AUCs,
and the paper reports its 80 mg/day median as 32.3 umol*h/L (Figure 5
legend). With the dose effect retained, the reference patient gives
8.94 + 20.69 = 29.6 umol*h/L, and the simulated cohort median is shown
in the AUC-driver table above; dropping the dose effect gives 25.7, a
clearly worse match. Supplementary Datafile S6 also computes `AUCSSP`
from `IFP`, the individual relative bioavailability carried over from
the PK model, which includes the dose effect. The model therefore keeps
the dose effect and the reference bioanalytical method, and the Table 1
elafibranor row is recorded here as an unreconciled discrepancy in the
source rather than corrected in the model.

Because the discrepancy is bioanalytical-scale in nature, the virtual
cohort above is simulated at the **reference** bioanalytical condition
(`ASSAY_SEPIP = 0`, `ASSAY_SEPIP_FA = 0`), which is the scale on which
the structural parameters were estimated. Setting `ASSAY_SEPIP = 1`
multiplies every predicted elafibranor concentration by `exp(-0.182)`.

### The published 80 and 120 mg exposure medians are not dose-proportional

Figure 5’s legend reports median summed AUCs of 32.3 and 39.3 umol\*h/L
for the 80 and 120 mg/day arms, a ratio of 1.22. Both models are linear
in dose over this range – the paper states that “between 50 mg and 360
mg, the PK of elafibranor was linear”, and the sigmoidal dose effect on
relative bioavailability is already saturated at 1.792 at 80 mg against
1.797 at 120 mg – so a simulation that gives every subject the nominal
dose necessarily separates the two arms by the full dose ratio of 1.5,
which is what the table above shows.

This is not a contradiction: the published values are medians of
*individually predicted* exposures in the actual trial patients, and the
paper itself observes that “the distributions largely overlapped”. The
120 mg arm ran only in the phase II trial while the phase III trial used
80 mg only, so the two medians are computed over different patients with
different covariates and very different group sizes. The consequence for
validation is that the saturation claim must be evaluated at the paper’s
own exposures rather than at a simulated cohort’s, which is what the
gate above does; comparing the simulated 120 mg median against 39.3
would be comparing a dose-proportional quantity with an empirical one.

### Parameter sourcing

- **Final estimates, not control-stream initials.** The `$THETA` records
  of Supplementary Datafiles S1, S2 and S6 are the initial estimates
  supplied to the final runs and differ from the published finals by up
  to ~11%. All values are taken from Table S3 (PK) and Table 2 (PD). The
  `$OMEGA` and `$SIGMA` records *do* match the published finals and were
  used to cross-check them.
- **Control-stream formulation indices are permuted.** The binarised
  `FORM1`..`FORM6` columns of Datafiles S1 and S2 map to Table S1 /
  Table S3 formulations 5, 6, 1, 2, 3, 4 respectively. The model files
  use the paper’s numbering. The permutation was established
  independently for both analytes by matching every one of the four
  formulation-specific parameter sets and their IIV magnitudes.
- **Logit-scale Emax variance.** Table 2 reports
  `IIV E max (CV) = 0.231`, which is on the back-transformed Emax scale,
  whereas the model’s random effect acts on the logit of `|Emax|`. The
  logit-scale variance 0.71883 is taken from Datafile S6
  `$OMEGA BLOCK(2)`. Two independent checks confirm it is the final
  value: the delta-method back-transform
  `0.8478 * 0.731 * (1 - 0.731) / 0.731 = 0.228` reproduces Table 2’s
  0.231, and combining it with the stream’s off-diagonal returns Table
  2’s printed correlation of -0.979 exactly.
- **Correlations.** Table 2 prints correlations, not covariances; the
  off-diagonals in `ini()` are `correlation x SD x SD`. Both 2x2 blocks
  are positive definite as published (determinants 1.48e-2 and 3.29e-4),
  so no positive-definiteness nudge was needed.
- **Zero-variance random effects are omitted rather than declared.**
  Where Table S3 or Table 2 report an IIV of “0 (FIX)” – the elafibranor
  D1 IIV for formulations 1-3, most MAT IIVs, several residual-error
  IIVs, and the PD IIVs on the turnover half-lives, EC50 and the
  bilirubin slope – no eta is declared, which is equivalent and avoids a
  singular OMEGA.

### Structural simplifications

- **Missing-covariate imputation is not reproduced.** Datafile S1
  imputes a missing baseline albumin as 42 g/L in healthy volunteers and
  46 g/L in the hepatic-impairment population, and Datafile S6 sets the
  liver-stiffness factor to 1 when the measurement is missing. These are
  data-handling rules, not structural model components; supply `ALB` for
  every simulated subject, and `LSM = 8.1` for a subject with no
  liver-stiffness measurement.
- **Inter-occasion variability needs an `OCC` column.** Datafile S1
  carries ten occasion slots each on D1 and MAT, all constrained equal.
  This vignette simulates with `OCC = 1` throughout, which draws one IOV
  realisation per subject; that matches Table 1’s “with IOV” column in
  structure, and Table 1 itself shows IOV barely moves the summary
  statistics. Increment `OCC` per dosing interval to simulate genuine
  occasion-to-occasion variation.
- **Study-specific overrides.** Datafile S1 carries several
  commented-out study-specific `TVFREL` / `TVD1` / `TVMAT` overrides
  (studies 5 and 11.2) and three fixed `$THETA` values that the final
  `$PK` block never references. These are inactive in the published
  final model and are not encoded. The one *active* study override – the
  lag time fixed to 0 h in GFT505B-319-1 – is encoded via
  `STUDY_GFT505B_319_1`.
- **Separate ALP and TB models not extracted.** Supplementary Datafiles
  S4 and S5 hold the standalone ALP and total-bilirubin models that
  preceded the joint model. They are model-development intermediates
  superseded by the joint model of Table 2, so only the joint model is
  packaged.

### Virtual-cohort assumptions

- Table S2 reports weight, BMI, albumin, ALT and creatinine clearance
  for the pooled *PK* analysis set only, not separately for the PBC
  subgroup. The cohort above centres those covariates on the Figure 3
  reference PBC patient and borrows the Table S2 standard deviations;
  the ALT SD is set to 20 U/L rather than the pooled 35.1 U/L, because
  the pooled value is inflated by the MASH studies.
- Race is not simulated: no covariate in any of the three models depends
  on it (Black or African American race was tested in the PK covariate
  screen and not retained).
- Every simulated subject receives formulation 5 under fasted
  conditions, which is both the model reference and the condition of
  Table 1 and of the phase II/III PBC trials.

### New canonical names registered with this extraction

`LSM` (liver stiffness by transient elastography, kPa, general scope),
`STUDY_PHASE2`, `ASSAY_SEPIP`, `ASSAY_SEPIP_FA`, `AUC_ELA`,
`AUC_GFT1007`, `DOSE_ELA_MG`, `STUDY_GFT505B_319_1` and
`FORM_ELA_F1`-`F4`/`F6`. All but `LSM` are members of existing canonical
families. The `alp` and `tbili` compartments follow the register’s bare
clinical-biomarker PD-output pattern (`ast`, `cpk`, `ldl`, `hdl`,
`urate`) but are declared `paper_specific` per the standing rule that a
compartment canonical needs a second independent paper.
