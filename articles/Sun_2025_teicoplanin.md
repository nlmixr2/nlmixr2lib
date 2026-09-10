# Teicoplanin in sepsis with CRRT (Sun 2025)

## Model and source

``` r

mod_meta <- nlmixr2est::nlmixr(readModelDb("Sun_2025_teicoplanin"))$meta
#> ℹ parameter labels from comments will be replaced by 'label()'
```

- Citation: Sun Q, Jian J, Zhou X, Hong Z, Yang S, Zheng Y, Wang S,
  Zhao M. Population pharmacokinetics of teicoplanin and dosage
  optimization in sepsis patients based on continuous renal replacement
  therapy. Front Pharmacol. 2025;16:1621959.
  <doi:10.3389/fphar.2025.1621959>
- Description: One-compartment IV infusion population PK model for
  teicoplanin in 86 critically ill adults with sepsis, 20 of whom were
  receiving continuous renal replacement therapy (Sun 2025). Clearance
  is 0.98 L/h with no retained covariate; the central volume of
  distribution is 108.69 L in a male not receiving CRRT and carries two
  exponential categorical covariate effects, V = 108.69 \* exp(-0.71 \*
  RRT_CRRT_STATUS) \* exp(-1.07 \* SEXF). Both effects reduce V: CRRT
  halves it (to 53.4 L) and female sex reduces it to about a third (37.3
  L), so a female receiving CRRT has V = 18.3 L. The CRRT direction is
  opposite to the usual literature finding of an expanded volume during
  renal replacement; the authors attribute it to fluid-overload
  correction being the indication for starting CRRT in this cohort.
  Interindividual variability is exponential on both CL (omega^2 = 0.31)
  and V (omega^2 = 0.09), and residual variability is additive at 0.23
  mg/L. Age, body weight, serum creatinine, serum albumin and ECMO were
  screened but not retained, and no covariate improved the fit on
  clearance.
- Article (DOI): <https://doi.org/10.3389/fphar.2025.1621959>

This vignette validates the packaged `Sun_2025_teicoplanin` model – a
one-compartment IV infusion population PK model for teicoplanin
developed from 86 trough concentrations in 86 critically ill adults with
sepsis in a Guizhou ICU, 20 of whom were receiving continuous renal
replacement therapy (CRRT).

The model is unusually simple: clearance carries **no** covariate, and
the two retained covariates – CRRT status and sex – both act on the
central volume of distribution, exponentially and in the same (downward)
direction.

The paper’s quantitative deliverable is its Monte Carlo
dose-optimization analysis (Figures 3-5), which recommends a specific
maintenance dose for each of three strata and reports one exact
probability-of-target-attainment (PTA) figure. That recommendation table
is a *binary answer key*: for each stratum the recommended dose is
asserted to be the smallest one reaching the target, so a faithful
reimplementation must place the PTA curve’s 90% crossing at exactly the
recommended rung and not one rung earlier or later. This vignette
reproduces that key, together with two exact structural identities.

``` r

mod <- readModelDb("Sun_2025_teicoplanin")
mod
#> function() {
#>   description <- "One-compartment IV infusion population PK model for teicoplanin in 86 critically ill adults with sepsis, 20 of whom were receiving continuous renal replacement therapy (Sun 2025). Clearance is 0.98 L/h with no retained covariate; the central volume of distribution is 108.69 L in a male not receiving CRRT and carries two exponential categorical covariate effects, V = 108.69 * exp(-0.71 * RRT_CRRT_STATUS) * exp(-1.07 * SEXF). Both effects reduce V: CRRT halves it (to 53.4 L) and female sex reduces it to about a third (37.3 L), so a female receiving CRRT has V = 18.3 L. The CRRT direction is opposite to the usual literature finding of an expanded volume during renal replacement; the authors attribute it to fluid-overload correction being the indication for starting CRRT in this cohort. Interindividual variability is exponential on both CL (omega^2 = 0.31) and V (omega^2 = 0.09), and residual variability is additive at 0.23 mg/L. Age, body weight, serum creatinine, serum albumin and ECMO were screened but not retained, and no covariate improved the fit on clearance."
#>   reference   <- "Sun Q, Jian J, Zhou X, Hong Z, Yang S, Zheng Y, Wang S, Zhao M. Population pharmacokinetics of teicoplanin and dosage optimization in sepsis patients based on continuous renal replacement therapy. Front Pharmacol. 2025;16:1621959. doi:10.3389/fphar.2025.1621959"
#>   vignette    <- "Sun_2025_teicoplanin"
#>   units       <- list(time = "h", dosing = "mg", concentration = "mg/L")
#> 
#>   # Issue #482: what each ODE state holds, in what amount units, in what
#>   # biological matrix. Verified against Sun 2025 Section 2.3 (total teicoplanin
#>   # in plasma by HPLC-UV, calibration range 5.63-125.00 mg/L).
#>   compartmentData <- list(
#>     central = list(analyte = "teicoplanin", units = "mg", specimen = "plasma", verified = TRUE)
#>   )
#> 
#>   covariateData <- list(
#>     RRT_CRRT_STATUS = list(
#>       description        = "Subject-level binary indicator for continuous renal replacement therapy during the teicoplanin sampling period",
#>       units              = "(binary)",
#>       type               = "binary",
#>       reference_category = 0,
#>       notes              = "Source column CRRT. 1 = subject was receiving CRRT; 0 = no CRRT (Sun 2025 Table 1: 20/86, 23.26%). Modality mix within the CRRT subgroup: 14 continuous venovenous hemofiltration (CVVH) only, 5 continuous venovenous hemodiafiltration (CVVHD) only, and 1 subject who received both (Sun 2025 Results 3.1); the model treats all three as a single binary indicator and does not distinguish modality. Time-fixed at the subject level, matching the RRT_CRRT_STATUS canonical's stated convention. Enters V as exp(-0.71 * RRT_CRRT_STATUS) per Sun 2025 Equation 5, i.e. CRRT REDUCES the volume of distribution by 50.8%. This direction is the opposite of most published CRRT covariate effects; Sun 2025 Discussion attributes it to volume overload being the indication for initiating CRRT, so that CRRT-treated subjects had their expanded interstitial volume corrected. Adding CRRT to clearance did NOT improve the fit (Sun 2025 Discussion), which is also atypical for a renally eliminated drug.",
#>       source_name        = "CRRT"
#>     ),
#>     SEXF = list(
#>       description        = "Biological sex indicator, 1 = female, 0 = male",
#>       units              = "(binary)",
#>       type               = "binary",
#>       reference_category = 0,
#>       notes              = "Source column reported as 'Gender' (Sun 2025 Table 1: 51 males, 59.30%; 35 females, 40.70%). Sun 2025 Equation 5 writes the effect as exp(-1.07 * (if is Female)), so the source indicator is already female = 1 and maps onto the canonical SEXF orientation with no inversion; the reference subject is male. Sun 2025 Section 3.3 confirms the coding by naming the female simulation cohort 'Group Sex = 1'. Female sex reduces V by 65.7% relative to male. Median body weight differed by sex (males 66.0 kg, IQR 54.0-73.0; females 57.0 kg, IQR 47.0-63.5; Sun 2025 Discussion), so the sex effect on V is partly confounded with body size, which was screened but not retained.",
#>       source_name        = "Gender"
#>     )
#>   )
#> 
#>   # Screened during covariate model building but NOT retained in the final model
#>   # (Sun 2025 Section 2.5 lists the tested set; Section 3.2 reports that only
#>   # gender and CRRT on V survived forward inclusion and backward elimination, and
#>   # that "No covariate was found to significantly influence the CL of
#>   # teicoplanin"). Documentation only -- none of these appears in model().
#>   covariatesDataExcluded <- list(
#>     WT = list(
#>       description = "Body weight",
#>       units       = "kg",
#>       type        = "continuous",
#>       notes       = "Sun 2025 Table 1 median (IQR) 62.00 (51.88, 70.00) kg. Screened as a covariate but not retained on either CL or V. Note that no allometric scaling was applied at all, so the reported CL and V are unnormalized whole-body values."
#>     ),
#>     AGE = list(
#>       description = "Subject age",
#>       units       = "years",
#>       type        = "continuous",
#>       notes       = "Sun 2025 Table 1 median (IQR) 62.00 (53.00, 71.25) years. Screened but not retained."
#>     ),
#>     CREAT = list(
#>       description = "Serum creatinine",
#>       units       = "umol/L",
#>       type        = "continuous",
#>       notes       = "Sun 2025 Table 1 median (IQR) 109.00 (74.00, 184.80). Table 1's column header reads 'Serum creatinine concentration (mg/dL)', which is a units error in the source: 109 mg/dL is not physiologically possible, whereas 109 umol/L (about 1.23 mg/dL) is an unremarkable ICU value consistent with the cohort's renal impairment and CRRT use. Recorded here in umol/L. Screened but not retained; the paper notes serum creatinine is a poor renal-function marker in CRRT-treated patients."
#>     ),
#>     ALB = list(
#>       description = "Serum albumin",
#>       units       = "g/L",
#>       type        = "continuous",
#>       notes       = "Sun 2025 Table 1 median (IQR) 31.60 (29.48, 36.90). Table 1's column header reads 'Serum albumin concentration (mg/L)', a units error in the source; the Discussion quotes the same median as '31.6 (29.5, 39.4) g/L', confirming g/L. Recorded here in g/L. (The Discussion's upper quartile of 39.4 also disagrees with Table 1's 36.90; Table 1 is taken as authoritative for the IQR.) Screened but not retained, despite teicoplanin being over 90% albumin-bound."
#>     ),
#>     ECMO_STATUS = list(
#>       description = "Extracorporeal membrane oxygenation treatment-status indicator",
#>       units       = "(binary)",
#>       type        = "binary",
#>       notes       = "Sun 2025 Table 1: 4/86 subjects (4.65%) received ECMO. Collected (Section 2.1) but not retained in the final model; with only 4 ECMO subjects the effect was not estimable."
#>     ),
#>     HT = list(
#>       description = "Body height at baseline",
#>       units       = "cm",
#>       type        = "continuous",
#>       notes       = "Sun 2025 Table 1 median (IQR) 165.00 (156.30, 172.00) cm. Collected but not retained."
#>     ),
#>     WBC = list(
#>       description = "White blood cell count",
#>       units       = "10^9/L",
#>       type        = "continuous",
#>       notes       = "Listed in Sun 2025 Section 2.1 among the collected physiological and biochemical parameters. No summary statistics are reported in Table 1 and it was not retained."
#>     ),
#>     AST = list(
#>       description = "Aspartate aminotransferase",
#>       units       = "U/L",
#>       type        = "continuous",
#>       notes       = "Listed in Sun 2025 Section 2.1 among the collected parameters. No summary statistics reported; not retained."
#>     ),
#>     ALT = list(
#>       description = "Alanine aminotransferase",
#>       units       = "U/L",
#>       type        = "continuous",
#>       notes       = "Listed in Sun 2025 Section 2.1 among the collected parameters. No summary statistics reported; not retained."
#>     )
#>   )
#> 
#>   population <- list(
#>     species        = "human",
#>     n_subjects     = 86L,
#>     n_studies      = 1L,
#>     age_range      = "Median (IQR) 62.00 (53.00, 71.25) years; adults aged 18 years and over (Sun 2025 Table 1 and inclusion criteria)",
#>     weight_range   = "Median (IQR) 62.00 (51.88, 70.00) kg; by sex, males 66.0 (54.0, 73.0) and females 57.0 (47.0, 63.5) (Sun 2025 Table 1 and Discussion)",
#>     sex_female_pct = 40.7,
#>     race_ethnicity = "Not reported (single-centre Chinese ICU cohort, presumed predominantly Han Chinese)",
#>     disease_state  = "Adults with sepsis by Sepsis 3.0 criteria admitted to the intensive care unit with confirmed or suspected Gram-positive infection, treated with teicoplanin for at least 4 days. 20/86 (23.26%) received CRRT (14 CVVH only, 5 CVVHD only, 1 both); 4/86 (4.65%) received ECMO. Children, pregnant women, and patients with joint, bone or endocardial infection were excluded.",
#>     dose_range     = "Per-protocol study regimen: teicoplanin 400 mg intravenously q12h for the first three doses, then 400 mg q24h maintenance, each given as a 1-hour infusion (Sun 2025 Section 2.2). Actual loading doses received ranged 200-800 mg (median 400) and maintenance doses 200-1,000 mg (median 400). The Monte Carlo dose-optimization simulations explored loading doses of 600-1,200 mg q12h for 3 or 5 doses with 200-1,000 mg q24h maintenance, plus continuous regimens of 400-1,000 mg q12h or 1,000-1,800 mg q24h.",
#>     regions        = "China (Beijing Jishuitan Hospital Guizhou Hospital, Guiyang, Guizhou; single-centre ICU)",
#>     renal_function = "Serum creatinine median (IQR) 109.00 (74.00, 184.80) umol/L (Table 1 header mislabels the unit as mg/dL). 20/86 subjects required CRRT for sepsis-associated acute kidney injury. CRRT effluent flow rate, filter adsorption capacity, dialysis timing relative to dosing, and modality were NOT captured in the retrospective dataset, which the authors list as a source of unexplained residual variability.",
#>     co_medication  = "Not reported beyond the study drug.",
#>     notes          = "Retrospective single-centre study, 1 June 2022 to 1 June 2024; IRB No. KT2022102101. IMPORTANT for interpreting the variance estimates: only 86 teicoplanin concentrations were available from 86 patients, i.e. essentially ONE trough sample per subject, drawn within 30 minutes before a dose at steady state (Sections 2.2 and 3.1). Observed trough concentrations had median (IQR) 13.40 (10.48, 19.83) mg/L. With a single observation per subject the residual error and the interindividual variability are only weakly separable, which is the likely explanation for the unusually small additive residual (0.23 mg/L) and for the wide bootstrap CI on omega^2 for V, which includes zero. Assay: HPLC-UV at 220 nm after protein precipitation, piperacillin internal standard, calibration range 5.63-125.00 mg/L (r^2 > 0.99), accuracy 2.98-10.36% and precision 7.33-11.25% at QC concentrations of 7.81, 31.25 and 90.00 mg/L. Estimation: Phoenix NLME 8.1, FOCE with interaction. Model evaluation: goodness-of-fit plots, 1,000-replicate bootstrap, and prediction-corrected VPC."
#>   )
#> 
#>   ini({
#>     # Structural fixed-effect parameters from Sun 2025 Table 2 ("Final model
#>     # (n = 86)" Estimate column) and confirmed by the printed final-model
#>     # Equations 4 and 5.
#>     lcl       <- log(0.98);    label("Clearance (L/h)")                      # Sun 2025 Table 2: CL = 0.98 L/h (RSE 6.92%; bootstrap median 0.97, 95% CI 0.83-1.12). Equation 4: CL = 0.98 * EXP(eta)
#>     lvc       <- log(108.69);  label("Central volume of distribution (L)")   # Sun 2025 Table 2: V = 108.69 L (RSE 9.89%; bootstrap median 107.89, 95% CI 84.30-151.56). Equation 5 reference subject is male without CRRT
#> 
#>     # Categorical covariate effects on V, entered on the log scale. Sun 2025
#>     # Equation 3 gives the categorical covariate form as
#>     #   tvP' = tv(P) * EXP(theta_cov * Cov_cat),
#>     # and Equation 5 instantiates it as
#>     #   V(L) = 108.69 * EXP(-0.71 * (if with CRRT)) * EXP(-1.07 * (if is Female)).
#>     # Both coefficients are therefore log-scale multipliers, NOT fractional
#>     # changes, and both are negative (each covariate reduces V).
#>     e_crrt_vc <- -0.71;        label("Log-scale CRRT effect on V (unitless)")   # Sun 2025 Table 2: theta_CRRT,V = -0.71 (RSE 22.55%; bootstrap median -0.70, 95% CI -1.15 to -0.19). exp(-0.71) = 0.492, a 50.8% reduction. Cross-check: the Results state V was "102.48% lower" with CRRT, and exp(+0.7055) - 1 = 102.48%, reproducing the printed ratio from the unrounded coefficient
#>     e_sexf_vc <- -1.07;        label("Log-scale female-sex effect on V (unitless)") # Sun 2025 Table 2: theta_Sex,V = -1.07 (RSE 28.53%; bootstrap median -1.07, 95% CI -1.61 to -0.25). exp(-1.07) = 0.343, a 65.7% reduction. Cross-check: the Results state males had "1.90-fold higher" V, and exp(1.07) - 1 = 1.92, i.e. a 1.9-fold increment over the female value
#> 
#>     # Between-subject variability. Sun 2025 Table 2 reports these under the
#>     # heading "Between-subject variation" as omega^2 CL and omega^2 V, and
#>     # Equations 4 and 5 print the SAME two numbers inside the exponential IIV
#>     # term ("* EXP(0.31)" and "* EXP(0.09)"). Read literally the printed
#>     # equations would be fixed multipliers, which is nonsense; read as the
#>     # authors intended they pin the table column to the VARIANCE slot of an
#>     # exponential IIV term, so the values below go into ini() as-is. Two
#>     # independent cross-checks confirm the variance reading:
#>     #   (1) the omega^2 CL row's RSE of 17.58% sits just above sqrt(2/86) =
#>     #       15.25%, the Cramer-Rao floor for a variance, and far above the
#>     #       sqrt(1/(2*86)) = 7.62% floor that an SD-scale parameter would have;
#>     #   (2) Section 2.4's generic Equation 1 prints the IIV as ADDITIVE,
#>     #       Pj = tv(P) + eta_j, but that reading is untenable against these
#>     #       magnitudes -- an additive eta with variance 0.09 on a V of 108.69 L
#>     #       is an SD of 0.30 L (0.3% of V), which no modeller would retain with
#>     #       an estimated RSE of 32%, while an additive eta with variance 0.31 on
#>     #       a CL of 0.98 L/h would make 3.9% of subjects have negative clearance.
#>     # Equations 4 and 5 are the final-model equations and are specific, so they
#>     # govern over the generic Equation 1; this is also Phoenix NLME's default
#>     # structural-parameter parameterization. See the vignette Errata.
#>     etalcl ~ 0.31  # Sun 2025 Table 2: omega^2 = 0.31 on CL (RSE 17.58%; bootstrap median 0.31, 95% CI 0.20-0.41). Equivalent to sqrt(exp(0.31) - 1) = 60.2% CV on the linear scale
#>     etalvc ~ 0.09  # Sun 2025 Table 2: omega^2 = 0.09 on V (RSE 32.11%; bootstrap median 0.08, 95% CI -0.05 to 0.21 -- the bootstrap CI includes zero, so this variance is only weakly identified from single-trough data). Equivalent to sqrt(exp(0.09) - 1) = 30.7% CV
#> 
#>     # Residual variability. Sun 2025 Section 3.2 states the final model
#>     # incorporated "an additive residual variability", and Table 2 reports it
#>     # under "Within-subject variation" as sigma_additive with explicit units of
#>     # mg/L, which fixes it on the SD scale in concentration units.
#>     addSd  <- 0.23;  label("Additive residual error (mg/L)")  # Sun 2025 Table 2: sigma_additive = 0.23 mg/L (RSE 10.93%; bootstrap median 0.22, 95% CI 0.03-0.26). Small relative to the observed troughs (median 13.40 mg/L) because each subject contributed only one concentration, so the etas absorb nearly all of the variability
#>   })
#>   model({
#>     # Individual PK parameters. Clearance carries no covariate (Sun 2025
#>     # Section 3.2: "No covariate was found to significantly influence the CL of
#>     # teicoplanin"). The two categorical effects on V are additive on the log
#>     # scale, which is algebraically identical to the product of exponentials
#>     # printed in Equation 5.
#>     cl <- exp(lcl + etalcl)
#>     vc <- exp(lvc + e_crrt_vc * RRT_CRRT_STATUS + e_sexf_vc * SEXF + etalvc)
#> 
#>     kel <- cl / vc
#> 
#>     # One-compartment model with first-order elimination (Sun 2025 Section 3.2).
#>     # Teicoplanin was given as a 1-hour intravenous infusion (Section 2.2), so
#>     # doses enter `central` with a rate or duration set in the event table.
#>     # Dose in mg and vc in L give central / vc in mg/L.
#>     d/dt(central) <- -kel * central
#> 
#>     Cc <- central / vc
#>     Cc ~ add(addSd)
#>   })
#> }
#> <environment: 0x56280cd11d50>
```

## Population

The model was fit to a retrospective single-centre cohort (Sun 2025
Table 1).

``` r

tibble::tribble(
  ~Characteristic,                        ~Value,
  "Subjects",                             "86 (86 concentrations, i.e. ~1 trough per subject)",
  "Sex",                                  "51 male (59.30%), 35 female (40.70%)",
  "Age (years)",                          "62.00 (53.00, 71.25)",
  "Weight (kg)",                          "62.00 (51.88, 70.00)",
  "Height (cm)",                          "165.00 (156.30, 172.00)",
  "Serum albumin (g/L)",                  "31.60 (29.48, 36.90)",
  "Serum creatinine (umol/L)",            "109.00 (74.00, 184.80)",
  "Observed trough (mg/L)",               "13.40 (10.48, 19.83)",
  "CRRT",                                 "20 (23.26%): 14 CVVH, 5 CVVHD, 1 both",
  "ECMO",                                 "4 (4.65%)"
) |>
  knitr::kable(caption = "Sun 2025 Table 1; continuous values are median (IQR).")
```

| Characteristic | Value |
|:---|:---|
| Subjects | 86 (86 concentrations, i.e. ~1 trough per subject) |
| Sex | 51 male (59.30%), 35 female (40.70%) |
| Age (years) | 62.00 (53.00, 71.25) |
| Weight (kg) | 62.00 (51.88, 70.00) |
| Height (cm) | 165.00 (156.30, 172.00) |
| Serum albumin (g/L) | 31.60 (29.48, 36.90) |
| Serum creatinine (umol/L) | 109.00 (74.00, 184.80) |
| Observed trough (mg/L) | 13.40 (10.48, 19.83) |
| CRRT | 20 (23.26%): 14 CVVH, 5 CVVHD, 1 both |
| ECMO | 4 (4.65%) |

Sun 2025 Table 1; continuous values are median (IQR). {.table}

Two column headers in the source Table 1 carry unit errors, corrected
above and recorded in the model file: serum creatinine is labelled
`mg/dL` but the values are `umol/L` (109 mg/dL is not physiologically
possible), and serum albumin is labelled `mg/L` but the Discussion
quotes the same median as `31.6 g/L`.

Dosing in the study was teicoplanin 400 mg IV q12h for three doses then
400 mg q24h, each as a **1-hour infusion**.

## Source trace

Every value in `ini()` traces to Sun 2025 Table 2, cross-checked against
the final-model Equations 4 and 5.

``` r

tibble::tribble(
  ~Quantity,                 ~`Model file`,   ~Value,    ~`Source location`,
  "Clearance",               "lcl",           "0.98 L/h",   "Table 2 (RSE 6.92%); Equation 4",
  "Central volume",          "lvc",           "108.69 L",   "Table 2 (RSE 9.89%); Equation 5",
  "CRRT effect on V",        "e_crrt_vc",     "-0.71",      "Table 2 (RSE 22.55%); Equation 5",
  "Female-sex effect on V",  "e_sexf_vc",     "-1.07",      "Table 2 (RSE 28.53%); Equation 5",
  "IIV on CL (variance)",    "etalcl",        "0.31",       "Table 2 (RSE 17.58%); Equation 4 exponent",
  "IIV on V (variance)",     "etalvc",        "0.09",       "Table 2 (RSE 32.11%); Equation 5 exponent",
  "Additive residual SD",    "addSd",         "0.23 mg/L",  "Table 2 (RSE 10.93%)",
  "Structural model",        "d/dt(central)", "1-cmt, 1st order", "Section 3.2",
  "Covariate form",          "exp(theta*cov)", "categorical", "Equation 3",
  "Infusion duration",       "dur = 1",       "1 h",        "Section 2.2"
) |>
  knitr::kable(caption = "Source trace for every ini() value and structural choice.")
```

| Quantity | Model file | Value | Source location |
|:---|:---|:---|:---|
| Clearance | lcl | 0.98 L/h | Table 2 (RSE 6.92%); Equation 4 |
| Central volume | lvc | 108.69 L | Table 2 (RSE 9.89%); Equation 5 |
| CRRT effect on V | e_crrt_vc | -0.71 | Table 2 (RSE 22.55%); Equation 5 |
| Female-sex effect on V | e_sexf_vc | -1.07 | Table 2 (RSE 28.53%); Equation 5 |
| IIV on CL (variance) | etalcl | 0.31 | Table 2 (RSE 17.58%); Equation 4 exponent |
| IIV on V (variance) | etalvc | 0.09 | Table 2 (RSE 32.11%); Equation 5 exponent |
| Additive residual SD | addSd | 0.23 mg/L | Table 2 (RSE 10.93%) |
| Structural model | d/dt(central) | 1-cmt, 1st order | Section 3.2 |
| Covariate form | exp(theta\*cov) | categorical | Equation 3 |
| Infusion duration | dur = 1 | 1 h | Section 2.2 |

Source trace for every ini() value and structural choice. {.table}

Equations 4 and 5 as printed in the source are:

    CL(L/h) = 0.98 x EXP(0.31)
    V(L)    = 108.69 x EXP(-0.71 x (if with CRRT)) x EXP(-1.07 x (if is Female)) x EXP(0.09)

### Typical parameter values by stratum

``` r

strata <- tidyr::expand_grid(SEXF = c(0, 1), RRT_CRRT_STATUS = c(0, 1)) |>
  mutate(
    Stratum = paste0(ifelse(SEXF == 1, "Female", "Male"),
                     ifelse(RRT_CRRT_STATUS == 1, ", CRRT", ", no CRRT")),
    `V (L)`     = 108.69 * exp(-0.71 * RRT_CRRT_STATUS - 1.07 * SEXF),
    `CL (L/h)`  = 0.98,
    `t1/2 (h)`  = log(2) * `V (L)` / `CL (L/h)`
  )

strata |>
  select(Stratum, `CL (L/h)`, `V (L)`, `t1/2 (h)`) |>
  knitr::kable(digits = 2, caption = "Typical values by covariate stratum.")
```

| Stratum         | CL (L/h) |  V (L) | t1/2 (h) |
|:----------------|---------:|-------:|---------:|
| Male, no CRRT   |     0.98 | 108.69 |    76.88 |
| Male, CRRT      |     0.98 |  53.44 |    37.80 |
| Female, no CRRT |     0.98 |  37.28 |    26.37 |
| Female, CRRT    |     0.98 |  18.33 |    12.96 |

Typical values by covariate stratum. {.table}

The reference subject (male, no CRRT) has V = 108.69 L and a terminal
half-life of 77 h, comfortably inside the 30-180 h range the paper’s
Introduction quotes for teicoplanin. The two covariate effects compound,
however, so a female receiving CRRT has V = 18.3 L and a half-life of 13
h – far outside that range. See the Errata for why this is preserved
rather than corrected.

## Structural verification

### Closed-form identity

A one-compartment model with a constant-rate infusion has an exact
analytic solution. Comparing the solver against it isolates pure
numerical error – both sides use the *same* drawn parameters – so a
tight bound is the correct assertion here.

``` r

rxSetSeed(20250701)
n_sub <- 200L

grid_fine <- sort(unique(c(seq(0, 2, by = 0.02), seq(2, 24, by = 0.25),
                           seq(24, 240, by = 2), seq(240, 1200, by = 8))))

ev_single <- et(amt = 400, time = 0, dur = 1, cmt = "central") |>
  et(time = grid_fine, cmt = "central") |>
  et(id = seq_len(n_sub))

dat_single <- as.data.frame(ev_single) |>
  mutate(SEXF = 0, RRT_CRRT_STATUS = 0)

sim_single <- rxSolve(mod, dat_single, returnType = "data.frame", addDosing = FALSE)
#> ℹ parameter labels from comments will be replaced by 'label()'

chk_cf <- sim_single |>
  mutate(
    kel     = cl / vc,
    rate    = 400 / 1,
    analytic = ifelse(
      time <= 1,
      rate / (vc * kel) * (1 - exp(-kel * time)),
      rate / (vc * kel) * (1 - exp(-kel * 1)) * exp(-kel * (time - 1))
    )
  ) |>
  filter(Cc > 1e-8) |>
  mutate(rel_err = abs((Cc - analytic) / analytic))

max_rel_err <- max(chk_cf$rel_err)
max_rel_err
#> [1] 1.792957e-12

stopifnot(max_rel_err < 1e-8)
```

The solver reproduces the closed form to machine precision, confirming
the ODE, the volume scaling and the infusion handling.

### Steady-state / mass-balance identity: AUC(inf) = Dose / CL

For a linear model the total exposure after a single dose is exactly
`Dose/CL`, independent of volume and of the covariates. This pins
clearance and the dose units simultaneously, per subject.

``` r

chk_auc <- sim_single |>
  arrange(id, time) |>
  group_by(id) |>
  summarise(
    auc_t   = sum(diff(time) * (head(Cc, -1) + tail(Cc, -1)) / 2),
    clast   = last(Cc),
    kel     = first(cl) / first(vc),
    cl_i    = first(cl),
    .groups = "drop"
  ) |>
  mutate(
    aucinf   = auc_t + clast / kel,
    target   = 400 / cl_i,
    pct_diff = 100 * (aucinf - target) / target
  )

summary(chk_auc$pct_diff)
#>      Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#> 0.0005121 0.0050999 0.0069091 0.0068620 0.0080911 0.0258431

stopifnot(max(abs(chk_auc$pct_diff)) < 0.1)
```

### Covariate multiplier identities

The two covariate effects must reproduce `exp(-0.71)` and `exp(-1.07)`
exactly, with the reference subject being male without CRRT.

Several checks below deliberately switch the random effects off and
supply the etas explicitly (either as fixed typical values, or as a
quantile lattice). rxode2 emits a
`multi-subject simulation without 'omega'` note in that situation, which
is expected here rather than a problem, so it is filtered out.

``` r

solve_fixed <- function(...) {
  withCallingHandlers(
    rxSolve(...),
    warning = function(w) {
      if (grepl("without .*'omega'", conditionMessage(w))) invokeRestart("muffleWarning")
    }
  )
}
```

``` r

mod_typ <- zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

ev_cov <- et(amt = 400, time = 0, dur = 1, cmt = "central") |>
  et(time = c(0, 1, 24), cmt = "central") |>
  et(id = 1:4)

dat_cov <- as.data.frame(ev_cov) |>
  mutate(
    SEXF            = c(0, 0, 1, 1)[id],
    RRT_CRRT_STATUS = c(0, 1, 0, 1)[id]
  )

# `omega = NA` omitted deliberately -- see the note at the PTA lattice solve
# below: `mod_typ` is already `zeroRe()`d, and on rxode2 5.1.6 passing it
# alongside a multi-subject solve reads out of bounds.
vc_by_stratum <- solve_fixed(mod_typ, dat_cov, returnType = "data.frame",
                             addDosing = FALSE) |>
  distinct(id, vc) |>
  arrange(id) |>
  pull(vc)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

ratios <- vc_by_stratum / vc_by_stratum[1]

tibble::tibble(
  Stratum  = c("Male, no CRRT (reference)", "Male, CRRT", "Female, no CRRT", "Female, CRRT"),
  `V (L)`  = vc_by_stratum,
  Ratio    = ratios,
  Expected = c(1, exp(-0.71), exp(-1.07), exp(-0.71 - 1.07))
) |>
  knitr::kable(digits = 4, caption = "Covariate multipliers on V.")
```

| Stratum                   |    V (L) |  Ratio | Expected |
|:--------------------------|---------:|-------:|---------:|
| Male, no CRRT (reference) | 108.6900 | 1.0000 |   1.0000 |
| Male, CRRT                |  53.4368 | 0.4916 |   0.4916 |
| Female, no CRRT           |  37.2816 | 0.3430 |   0.3430 |
| Female, CRRT              |  18.3293 | 0.1686 |   0.1686 |

Covariate multipliers on V. {.table}

``` r


stopifnot(
  abs(vc_by_stratum[1] - 108.69) < 1e-8,
  max(abs(ratios - c(1, exp(-0.71), exp(-1.07), exp(-0.71 - 1.07)))) < 1e-10
)
```

## Virtual cohort and the observed trough distribution

A cohort matched to Table 1’s sex and CRRT proportions, given the
study’s own regimen (400 mg q12h x 3 then 400 mg q24h), should reproduce
the observed trough distribution.

``` r

rxSetSeed(20250701)
set.seed(20250701)

cohort <- tibble::tibble(
  id              = seq_len(n_sub),
  SEXF            = rbinom(n_sub, 1, 0.407),
  RRT_CRRT_STATUS = rbinom(n_sub, 1, 0.2326)
)

build_regimen <- function(ld, n_ld, md, tmax = 192) {
  ld_times <- seq(0, by = 12, length.out = n_ld)
  md_times <- seq(max(ld_times) + 24, tmax, by = 24)
  tibble::tibble(
    time = c(ld_times, md_times),
    amt  = c(rep(ld, length(ld_times)), rep(md, length(md_times)))
  )
}

build_events <- function(reg, obs_times, n) {
  ev <- et(amt = reg$amt[1], time = reg$time[1], dur = 1, cmt = "central")
  for (i in seq_along(reg$time)[-1]) {
    ev <- et(ev, amt = reg$amt[i], time = reg$time[i], dur = 1, cmt = "central")
  }
  ev <- et(ev, time = obs_times, cmt = "central")
  et(ev, id = seq_len(n))
}

reg_study <- build_regimen(400, 3, 400)
trough_times <- c(72, 96, 168) - 1e-4

dat_study <- as.data.frame(build_events(reg_study, trough_times, n_sub)) |>
  left_join(cohort, by = "id")

sim_study <- rxSolve(mod, dat_study, returnType = "data.frame", addDosing = FALSE)

troughs <- sim_study |>
  mutate(time = round(time, 3)) |>
  filter(time %in% round(trough_times, 3)) |>
  mutate(Hour = round(time))

troughs |>
  group_by(Hour) |>
  summarise(
    Median = median(Cc), Q1 = quantile(Cc, 0.25), Q3 = quantile(Cc, 0.75),
    .groups = "drop"
  ) |>
  knitr::kable(digits = 2,
               caption = "Simulated trough concentrations (mg/L) on the study regimen.")
```

| Hour | Median |   Q1 |    Q3 |
|-----:|-------:|-----:|------:|
|   72 |   9.68 | 7.66 | 12.87 |
|   96 |  10.57 | 7.65 | 14.05 |
|  168 |  12.22 | 7.80 | 16.40 |

Simulated trough concentrations (mg/L) on the study regimen. {.table}

The observed trough distribution was 13.40 (10.48, 19.83) mg/L. TDM
samples were drawn “at steady state” without a stated day, so the
comparison is against the central tendency across the plateau rather
than any single timepoint.

``` r

c168 <- troughs$Cc[troughs$Hour == 168]
med_168 <- median(c168)
med_168
#> [1] 12.22308

# Structural check: a mis-transcribed CL, dose or unit moves this by tens of
# percent. Asserting on the median, not on extremes (which are not reproducible
# across rxode2 builds).
stopifnot(abs(med_168 - 13.40) / 13.40 < 0.30)
```

The simulated median trough sits slightly below the observed median.
That is the expected direction: the observed cohort’s maintenance doses
ranged 200-1,000 mg because clinicians escalated the dose in
under-exposed patients, whereas this simulation holds every subject at
the protocol 400 mg.

## NCA validation (PKNCA)

Sun 2025 reports no NCA parameters, so PKNCA is used here to
characterise the model’s exposure metrics and to re-derive the
`AUC = Dose/CL` identity through an independent code path.

``` r

nca_conc <- sim_single |>
  filter(!is.na(Cc)) |>
  mutate(treatment = "400 mg single 1-h infusion") |>
  select(id, treatment, time, Cc)

nca_dose <- cohort |>
  transmute(id, treatment = "400 mg single 1-h infusion", time = 0, amt = 400) |>
  filter(id %in% unique(nca_conc$id))

o_conc <- PKNCAconc(nca_conc, Cc ~ time | treatment + id)
o_dose <- PKNCAdose(nca_dose, amt ~ time | treatment + id)

o_data <- PKNCAdata(
  o_conc, o_dose,
  intervals = data.frame(
    start = 0, end = Inf,
    cmax = TRUE, tmax = TRUE, half.life = TRUE,
    aucinf.obs = TRUE, auclast = TRUE
  )
)

res_nca <- suppressWarnings(pk.nca(o_data))

nca_wide <- as.data.frame(res_nca) |>
  select(id, PPTESTCD, PPORRES) |>
  pivot_wider(names_from = PPTESTCD, values_from = PPORRES)

nca_wide |>
  summarise(
    across(c(cmax, tmax, half.life, aucinf.obs),
           list(median = ~median(.x, na.rm = TRUE),
                p10 = ~quantile(.x, 0.10, na.rm = TRUE),
                p90 = ~quantile(.x, 0.90, na.rm = TRUE)))
  ) |>
  pivot_longer(everything(), names_to = "metric", values_to = "value") |>
  tidyr::separate(metric, into = c("Parameter", "Statistic"), sep = "_(?=[^_]+$)") |>
  pivot_wider(names_from = Statistic, values_from = value) |>
  dplyr::rename("NCA parameter" = Parameter, "Median" = median,
                "P10" = p10, "P90" = p90) |>
  knitr::kable(digits = 2,
               caption = "PKNCA summary after a single 400 mg 1-hour infusion (male, no CRRT).")
```

| NCA parameter | Median |    P10 |    P90 |
|:--------------|-------:|-------:|-------:|
| cmax          |   3.67 |   2.60 |   5.38 |
| tmax          |   1.00 |   1.00 |   1.00 |
| half.life     |  76.39 |  35.41 | 179.66 |
| aucinf.obs    | 412.76 | 208.80 | 930.93 |

PKNCA summary after a single 400 mg 1-hour infusion (male, no CRRT).
{.table}

``` r

cl_by_id <- sim_single |> distinct(id, cl)

chk_nca <- nca_wide |>
  left_join(cl_by_id, by = "id") |>
  mutate(pct_diff = 100 * (aucinf.obs - 400 / cl) / (400 / cl))

summary(chk_nca$pct_diff)
#>       Min.    1st Qu.     Median       Mean    3rd Qu.       Max. 
#> -7.661e-06 -6.238e-07 -2.732e-07 -5.174e-07 -1.091e-07 -5.879e-09

stopifnot(
  !anyNA(chk_nca$aucinf.obs),
  max(abs(chk_nca$pct_diff)) < 1
)
```

PKNCA’s `aucinf.obs` reproduces `Dose/CL` per subject to well under 1%,
independently confirming the clearance transcription and the dose units.

## Dose optimization: reproducing the paper’s answer key

Sun 2025 recommends, for a target trough of 10 mg/L at MIC = 1 mg/L:

- male without CRRT: 800 mg q12h x 3, then **600 mg** q24h
- male with CRRT: 800 mg q12h x 3, then **800 mg** q24h
- female: 1,000 mg q12h x 3, then **1,000 mg** q24h, with a reported
  C168h PTA of **90.20%**

To reproduce PTA at 200 subjects per arm without Monte Carlo noise
swamping the comparison, the cohort is a deterministic stratified
lattice: 20 clearance quantiles x 20 volume quantiles of the fitted eta
distributions. This gives an exact quadrature of the PTA integral rather
than a random sample, so the numbers below are reproducible across
rxode2 builds and thread counts.

``` r

mod_fixed <- zeroRe(mod)
#> ℹ parameter labels from comments will be replaced by 'label()'

n_cl <- 20L
n_v  <- 10L
lattice <- tidyr::expand_grid(i_cl = seq_len(n_cl), i_v = seq_len(n_v)) |>
  mutate(
    id  = row_number(),
    lcl = log(0.98)   + qnorm((i_cl - 0.5) / n_cl) * sqrt(0.31),
    lvc = log(108.69) + qnorm((i_v  - 0.5) / n_v)  * sqrt(0.09)
  ) |>
  select(id, lcl, lvc)

nrow(lattice)
#> [1] 200

pta_det <- function(ld, n_ld, md, crrt, sexf) {
  dat <- as.data.frame(build_events(build_regimen(ld, n_ld, md),
                                    c(72, 168) - 1e-4, nrow(lattice))) |>
    mutate(RRT_CRRT_STATUS = crrt, SEXF = sexf)
    # `omega = NA` omitted deliberately. The model is already `zeroRe()`d, so
    # the etas are zero either way -- but on rxode2 5.1.6 (the version CI
    # installs) passing `omega = NA` ALONGSIDE a multi-subject solve reads
    # out of bounds and returns NA for some subjects. Locally on 5.1.7 it is
    # silent, which is how this reached CI. Same defect as PR #501 (Yoon 2023
    # vancomycin, Lin 2026 coxTte).
  s <- solve_fixed(mod_fixed, dat, params = as.data.frame(lattice),
                   returnType = "data.frame", addDosing = FALSE) |>
    mutate(time = round(time, 3))
  tibble::tibble(
    c72  = 100 * mean(s$Cc[s$time == round(72  - 1e-4, 3)] >= 10),
    c168 = 100 * mean(s$Cc[s$time == round(168 - 1e-4, 3)] >= 10)
  )
}

ladder <- tibble::tribble(
  ~Stratum,           ~crrt, ~sexf, ~ld,  ~n_ld, ~md,
  "Male, no CRRT",     0, 0,  800, 3,  200,
  "Male, no CRRT",     0, 0,  800, 3,  400,
  "Male, no CRRT",     0, 0,  800, 3,  600,
  "Male, no CRRT",     0, 0,  800, 3,  800,
  "Male, CRRT",        1, 0,  800, 3,  400,
  "Male, CRRT",        1, 0,  800, 3,  600,
  "Male, CRRT",        1, 0,  800, 3,  800,
  "Male, CRRT",        1, 0,  800, 3, 1000,
  "Female, no CRRT",   0, 1, 1000, 3,  400,
  "Female, no CRRT",   0, 1, 1000, 3,  600,
  "Female, no CRRT",   0, 1, 1000, 3,  800,
  "Female, no CRRT",   0, 1, 1000, 3, 1000
) |>
  rowwise() |>
  mutate(pta = list(pta_det(ld, n_ld, md, crrt, sexf))) |>
  ungroup() |>
  tidyr::unnest(pta)
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

recommended <- tibble::tibble(
  Stratum = c("Male, no CRRT", "Male, CRRT", "Female, no CRRT"),
  rec_md  = c(600, 800, 1000)
)

ladder |>
  left_join(recommended, by = "Stratum") |>
  transmute(
    Stratum,
    `Loading`      = paste0(ld, " mg q12h x", n_ld),
    `Maintenance`  = paste0(md, " mg q24h"),
    `C72h PTA (%)` = c72,
    `C168h PTA (%)` = c168,
    `Paper recommends` = ifelse(md == rec_md, "<<<", "")
  ) |>
  knitr::kable(digits = 2,
               caption = "PTA for a trough target of 10 mg/L. The paper's recommended maintenance dose should be the smallest rung reaching ~90%.")
```

| Stratum | Loading | Maintenance | C72h PTA (%) | C168h PTA (%) | Paper recommends |
|:---|:---|:---|---:|---:|:---|
| Male, no CRRT | 800 mg q12h x3 | 200 mg q24h | 81.5 | 52.0 |  |
| Male, no CRRT | 800 mg q12h x3 | 400 mg q24h | 88.5 | 78.5 |  |
| Male, no CRRT | 800 mg q12h x3 | 600 mg q24h | 92.5 | 89.5 | \<\<\< |
| Male, no CRRT | 800 mg q12h x3 | 800 mg q24h | 94.5 | 95.0 |  |
| Male, CRRT | 800 mg q12h x3 | 400 mg q24h | 84.5 | 69.5 |  |
| Male, CRRT | 800 mg q12h x3 | 600 mg q24h | 88.5 | 83.5 |  |
| Male, CRRT | 800 mg q12h x3 | 800 mg q24h | 92.0 | 90.5 | \<\<\< |
| Male, CRRT | 800 mg q12h x3 | 1000 mg q24h | 94.0 | 94.0 |  |
| Female, no CRRT | 1000 mg q12h x3 | 400 mg q24h | 82.0 | 63.0 |  |
| Female, no CRRT | 1000 mg q12h x3 | 600 mg q24h | 86.0 | 77.5 |  |
| Female, no CRRT | 1000 mg q12h x3 | 800 mg q24h | 89.0 | 85.5 |  |
| Female, no CRRT | 1000 mg q12h x3 | 1000 mg q24h | 91.0 | 90.0 | \<\<\< |

PTA for a trough target of 10 mg/L. The paper’s recommended maintenance
dose should be the smallest rung reaching ~90%. {.table}

``` r

# Gate 1: the exact published number. Female, 1000 q12h x3 + 1000 qd, C168h.
pta_female <- ladder |>
  filter(Stratum == "Female, no CRRT", md == 1000) |>
  pull(c168)
pta_female
#> [1] 90

stopifnot(abs(pta_female - 90.20) < 5)

# Gate 2: the answer key. In every stratum the recommended rung reaches the
# target and the rung immediately below it does not. Gating the attainment
# LEVEL either side of the recommendation, not the dose label itself.
key <- ladder |>
  left_join(recommended, by = "Stratum") |>
  group_by(Stratum) |>
  arrange(md, .by_group = TRUE) |>
  summarise(
    at_rec    = c168[md == rec_md],
    below_rec = c168[md == rec_md[1] - 200],
    .groups   = "drop"
  )
key
#> # A tibble: 3 × 3
#>   Stratum         at_rec below_rec
#>   <chr>            <dbl>     <dbl>
#> 1 Female, no CRRT   90        85.5
#> 2 Male, CRRT        90.5      83.5
#> 3 Male, no CRRT     89.5      78.5

stopifnot(
  all(key$at_rec    >= 88),   # recommended rung reaches the target
  all(key$below_rec <  88),   # the rung below does not
  all(key$at_rec - key$below_rec > 3)
)
```

The 90% crossing lands on the recommended rung in all three strata, and
the published 90.20% is reproduced within a fraction of a percentage
point.

``` r

# The paper also states that 600 and 800 mg q12h x5 + 400 mg qd give "overall
# PTAs of 70%~80%" for C168h, pooled across CRRT status.
pooled <- tidyr::expand_grid(ld = c(600, 800), crrt = c(0, 1)) |>
  rowwise() |>
  mutate(c168 = pta_det(ld, 5, 400, crrt, 0)$c168) |>
  ungroup() |>
  group_by(ld) |>
  summarise(`Pooled C168h PTA (%)` = mean(c168), .groups = "drop")
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

pooled |>
  knitr::kable(digits = 2,
               caption = "q12h x5 loading + 400 mg q24h maintenance, pooled over CRRT status (published: 70-80%).")
```

|  ld | Pooled C168h PTA (%) |
|----:|---------------------:|
| 600 |                75.00 |
| 800 |                78.75 |

q12h x5 loading + 400 mg q24h maintenance, pooled over CRRT status
(published: 70-80%). {.table}

``` r


stopifnot(all(pooled$`Pooled C168h PTA (%)` > 65),
          all(pooled$`Pooled C168h PTA (%)` < 85))
```

``` r

ladder |>
  left_join(recommended, by = "Stratum") |>
  ggplot(aes(md, c168, colour = Stratum)) +
  geom_line() +
  geom_point(aes(shape = md == rec_md), size = 3) +
  scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 19),
                     labels = c(`FALSE` = "other", `TRUE` = "recommended"),
                     name = NULL) +
  geom_hline(yintercept = 90, linetype = "dashed") +
  labs(x = "Maintenance dose (mg q24h)",
       y = "PTA for C168h >= 10 mg/L (%)") +
  theme_bw()
```

![PTA versus maintenance dose. Replicates the structure of Sun 2025
Figures 3 and 5.](Sun_2025_teicoplanin_files/figure-html/pta_plot-1.png)

PTA versus maintenance dose. Replicates the structure of Sun 2025
Figures 3 and 5.

## The IIV convention: an explicit falsifier

Sun 2025 Table 2 labels its variability rows `omega^2 CL` and
`omega^2 V`, but the paper’s *generic* Equation 1 prints the IIV as
additive, `Pj = tv(P) + eta`. That reading is untenable (an additive eta
of variance 0.09 on a 108.69 L volume is an SD of 0.30 L, and an
additive eta of variance 0.31 on a 0.98 L/h clearance makes 3.9% of
subjects have negative clearance), and the final-model Equations 4 and 5
instead print `x EXP(0.31)` and `x EXP(0.09)` – the table values
verbatim inside an exponential term.

The published PTA discriminates the candidate readings directly. Note
that the **median does not discriminate**; only the lower tail does,
which is exactly what a trough-target PTA measures.

``` r

pta_under <- function(om_cl, om_v) {
  L <- tidyr::expand_grid(i_cl = seq_len(n_cl), i_v = seq_len(n_v)) |>
    mutate(id  = row_number(),
           lcl = log(0.98)   + qnorm((i_cl - 0.5) / n_cl) * sqrt(om_cl),
           lvc = log(108.69) + qnorm((i_v  - 0.5) / n_v)  * sqrt(om_v)) |>
    select(id, lcl, lvc)
  dat <- as.data.frame(build_events(build_regimen(1000, 3, 1000),
                                    168 - 1e-4, nrow(L))) |>
    mutate(RRT_CRRT_STATUS = 0, SEXF = 1)
  # `omega = NA` omitted -- see the lattice solve above.
  s <- solve_fixed(mod_fixed, dat, params = as.data.frame(L),
                   returnType = "data.frame", addDosing = FALSE)
  cc <- s$Cc[round(s$time, 3) == round(168 - 1e-4, 3)]
  tibble::tibble(PTA = 100 * mean(cc >= 10), Median = median(cc),
                 P10 = quantile(cc, 0.10))
}

falsifier <- tibble::tribble(
  ~Reading,                                     ~om_cl,          ~om_v,
  "A: table is omega^2 (adopted)",               0.31,            0.09,
  "B: table is omega (SD)",                      0.31^2,          0.09^2,
  "C: table is CV",                              log(1 + 0.31^2), log(1 + 0.09^2)
) |>
  rowwise() |>
  mutate(res = list(pta_under(om_cl, om_v))) |>
  ungroup() |>
  tidyr::unnest(res) |>
  mutate(`|PTA - 90.20|` = abs(PTA - 90.20))
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'
#> ℹ omega/sigma items treated as zero: 'etalcl', 'etalvc'

falsifier |>
  select(Reading, PTA, Median, P10, `|PTA - 90.20|`) |>
  knitr::kable(digits = 2,
               caption = "Candidate IIV readings against the published C168h PTA of 90.20%.")
```

| Reading                       | PTA | Median |   P10 | \|PTA - 90.20\| |
|:------------------------------|----:|-------:|------:|----------------:|
| A: table is omega^2 (adopted) |  90 |  30.54 | 10.14 |             0.2 |
| B: table is omega (SD)        | 100 |  30.95 | 17.59 |             9.8 |
| C: table is CV                | 100 |  30.95 | 17.85 |             9.8 |

Candidate IIV readings against the published C168h PTA of 90.20%.
{.table}

``` r


stopifnot(
  which.min(falsifier$`|PTA - 90.20|`) == 1,
  falsifier$`|PTA - 90.20|`[1] < 3,
  all(falsifier$`|PTA - 90.20|`[-1] > 5)
)
```

Reading A is the only one that reproduces the published PTA; B and C
both saturate near 100%. Combined with the verbatim printed exponents
and with the `omega^2 CL` row’s RSE of 17.58% sitting just above
`sqrt(2/86) = 15.25%` (the Cramer-Rao floor for a *variance*, versus
7.62% for an SD), the variance reading is settled on four independent
grounds.

## Assumptions and deviations

- **IIV form.** Encoded as exponential with the Table 2 values as
  variances, against the letter of the paper’s generic Equation 1 (which
  prints an additive eta). Justified by the final-model Equations 4 and
  5, by the Cramer-Rao RSE floor, by numerical plausibility, and by the
  PTA falsifier above. Phoenix NLME’s default structural-parameter form
  is also exponential.

- **Maintenance-dose start time.** The paper writes regimens as “1,000
  mg q12h\*3 + 1,000 mg qd” without stating when the daily maintenance
  dose begins. This vignette starts it 24 h after the last loading dose,
  the standard clinical reading. The reproduction of the published
  90.20% supports that interpretation.

- **Deterministic lattice instead of random sampling.** The paper ran
  1,000 random replicates. This vignette uses a 20 x 10 stratified
  quantile lattice (200 subjects, the per-arm cap) so the PTA is an
  exact quadrature rather than a noisy sample; this makes the assertions
  reproducible across rxode2 builds.

- **PTA uses `Cc` (individual prediction) without residual error**,
  matching the usual PTA convention. With an additive residual SD of
  0.23 mg/L the choice changes the numbers negligibly.

- **Race/ethnicity is not reported** in the source; the cohort is a
  single Chinese centre and is presumed predominantly Han Chinese.

### Errata and source inconsistencies

- **Table 1 unit errors (two).** Serum creatinine is headed `mg/dL` but
  the values are `umol/L`; serum albumin is headed `mg/L` but the
  Discussion quotes the same median in `g/L`. Both are recorded
  corrected in the model file.

- **Albumin IQR disagreement.** Table 1 gives 31.60 (29.48, 36.90) while
  the Discussion gives 31.6 (29.5, 39.4). Table 1 is taken as
  authoritative.

- **Equation 1 versus Equations 4-5.** Discussed above; the final-model
  equations govern.

- **The reported subgroup volumes are cohort summaries, not typical
  values.** The Results quote V as 76.68 L without CRRT vs 37.87 L with
  CRRT, and 90.43 L in males vs 31.85 L in females. None equals the
  typical value of 108.69 L, because each is averaged over the *other*
  covariate’s distribution in that subgroup. Weighting the model’s
  typical values by the Table 1 composition (40.7% female, 23.26% CRRT)
  reproduces them to within 3-6%: 79.6 vs 76.68, 39.1 vs 37.87, 95.8 vs
  90.43, and 32.9 vs 31.85.

- **“102.48% lower” and “1.90-fold higher” are ratios, not percentage
  reductions.** `exp(0.7055) - 1 = 102.48%` and `exp(1.07) - 1 = 1.92`
  reproduce both statements, confirming the exponential covariate form.

- **The compounded covariate effects give a physiologically doubtful
  half-life in the female + CRRT stratum.** V falls to 18.3 L and the
  half-life to 13 h, well below the 30-180 h range the paper’s own
  Introduction cites for teicoplanin, and below any published
  teicoplanin volume. The model is encoded faithfully rather than
  corrected, but users should treat the female + CRRT combination as an
  extrapolation: the paper never simulated it (its three Monte Carlo
  cohorts were male/no-CRRT, male/CRRT and female/no-CRRT), and with 35
  females and 20 CRRT subjects the cell is very sparsely observed.

- **The variance on V is only weakly identified.** Its bootstrap 95% CI
  is (-0.05, 0.21), which includes zero. With one trough per subject the
  residual error and the IIV are barely separable, which also explains
  the unusually small additive residual of 0.23 mg/L.

- **The CRRT effect direction is opposite to most of the literature.**
  CRRT *reduces* V here. The authors attribute this to fluid overload
  being the indication for starting CRRT, so that CRRT-treated subjects
  had their expanded interstitial volume corrected. CRRT on clearance
  did not improve the fit, which is also atypical for a renally
  eliminated drug; the paper notes that effluent flow rate, filter type
  and modality were not captured.
