Liu_2026_riluzole <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption, an absorption lag time and linear elimination for oral riluzole in healthy adult volunteers (Liu 2026), with a power body-weight effect on the apparent peripheral volume Vp/F normalised to 75 kg (exponent 0.993)."
  reference <- paste(
    "Liu C, Tang Q, Yu L, Yu Z. Handling missed dose of riluzole by population",
    "pharmacokinetic model informed remedial regimen. Front Pharmacol.",
    "2026;17:1890355. doi:10.3389/fphar.2026.1890355.",
    sep = " "
  )
  vignette <- "Liu_2026_riluzole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # What each ODE state holds, in what amount units, in what biological matrix.
  # Verified against Liu 2026 Methods 'Data collection and PPK model
  # development for riluzole' (single oral tablet doses; plasma
  # concentration-time observations).
  compartmentData <- list(
    depot = list(analyte = "riluzole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "riluzole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "riluzole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the apparent peripheral volume only, normalised to 75 kg: Vp/F = 376 * (BW/75)^0.993 (Liu 2026 Results Equation 1; Methods: 'Body weight was normalized to 75 kg in the final model'). The 75 kg reference is above the cohort median of 63.3 kg (Supplementary Table S1). Cohort weight range 45.9-80.1 kg.",
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on Vp/F (Supplementary Table S2, models 5, 10, 13) and dropped in backward elimination (model 18, delta OFV 6.098, below the 10.83 retention threshold). Only 7 of 72 subjects (9.7%) were female."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL/F and Vc/F (Supplementary Table S2, models 6, 7, 11, 12, 14, 15); never significant in forward inclusion. Cohort age range 18-50 years."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 72L,
    n_studies = 1L,
    age_range = "18-50 years (median 27)",
    age_median = "27 years",
    weight_range = "45.9-80.1 kg (median 63.3)",
    weight_median = "63.3 kg",
    sex_female_pct = 9.7,
    disease_state = "Healthy adult volunteers. The model is intended by the authors for exposure prediction in patients with amyotrophic lateral sclerosis (ALS) without major organ dysfunction, but was not fitted to ALS-patient data.",
    dose_range = "Single 50 mg oral riluzole tablet of each of two formulations (test and reference), fasting, in a two-period crossover with washout.",
    regions = "China",
    notes = "Data from a previously published riluzole tablet bioequivalence study (Yu F et al. 2025, BMC Pharmacol Toxicol 26:105). 42 plasma samples per subject, 3,024 observations in total, no missing data. Concentrations below the 0.700 ng/mL LLOQ were set to LLOQ/2 (0.35 ng/mL) and kept in the analysis. Fitted in NONMEM 7.6.0 with FOCE-I. Demographics from Supplementary Table S1; parameter estimates from Table 1. The paper reports no formulation effect, so both tablet formulations share one set of parameters."
  )

  ini({
    # Structural PK parameters -- Liu 2026 Table 1 'Final model results'.
    # Vp/F is the typical value at the 75 kg reference weight.
    ltlag <- log(0.157); label("Absorption lag time (h)") # Table 1 'ALAG' = 0.157 h (RSE 1.3%)
    lka <- log(5.25); label("First-order absorption rate constant (1/h)") # Table 1 'Ka' = 5.25 h-1 (RSE 19%)
    lcl <- log(67); label("Apparent clearance CL/F (L/h)") # Table 1 'CL/F' = 67 L/h (RSE 4.9%)
    lvc <- log(242); label("Apparent central volume Vc/F (L)") # Table 1 'Vc/F' = 242 L (RSE 6.3%)
    lq <- log(33.1); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 'Q/F' = 33.1 L/h (RSE 5.4%)
    lvp <- log(376); label("Apparent peripheral volume Vp/F at WT = 75 kg (L)") # Table 1 'Vp/F' = 376 L (RSE 8.6%); Equation 1

    # Power exponent of body weight on Vp/F, estimated (RSE reported).
    e_wt_vp <- 0.993; label("Power exponent of (WT/75) on Vp/F (unitless)") # Table 1 'theta BW on Vp' = 0.993 (RSE 33.9%); Equation 1

    # Inter-individual variability, exponential (log-normal). Table 1 labels
    # the rows 'omega CL/Vc/Ka'; the values are NONMEM $OMEGA variances: the
    # base model in Supplementary Table S3 prints the same quantities as
    # 'omega2' (0.142, 0.214, 0.989), and the variance reading reproduces the
    # 5th-95th percentile band of Figure 4 (see the vignette). IIV on Q and Vp
    # was tested and removed (Results 'PPK model development for riluzole').
    etalcl ~ 0.141 # Table 1 'omega CL' = 0.141 (RSE 15%, shrinkage 1%)
    etalvc ~ 0.218 # Table 1 'omega Vc' = 0.218 (RSE 16.8%, shrinkage 5%)
    etalka ~ 0.994 # Table 1 'omega Ka' = 0.994 (RSE 14%, shrinkage 5%)

    # Combined proportional + additive residual error. Table 1 prints the
    # NONMEM $SIGMA variances; nlmixr2 takes SDs, so each value is the square
    # root of the printed variance. sqrt(0.127) = 0.356 matches the base
    # model's proportional error of 35.5 %CV (Supplementary Table S3).
    propSd <- 0.3564; label("Proportional residual SD (fraction)") # Table 1 'Proportional' = 0.127 (variance; RSE 7%); sqrt(0.127) = 0.3564
    addSd <- 0.2512; label("Additive residual SD (ng/mL)") # Table 1 'Additive' = 0.0631 (variance; RSE 30.1%); sqrt(0.0631) = 0.2512
  })

  model({
    tlag <- exp(ltlag)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp) * (WT / 75)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
