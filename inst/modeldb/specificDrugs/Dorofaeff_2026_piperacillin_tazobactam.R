Dorofaeff_2026_piperacillin_tazobactam <- function() {
  description <- paste(
    "Joint one-compartment population PK model for piperacillin and",
    "tazobactam in 31 critically ill children aged 1 month to 12 years in a",
    "paediatric intensive care unit (Dorofaeff 2026). Each drug has its own",
    "first-order-elimination central compartment with clearance and volume",
    "estimated separately. Both clearances carry allometric weight scaling",
    "(fixed exponent 0.75, 70 kg reference) and power terms on",
    "Schwartz-estimated GFR (reference 117.9 mL/min/1.73 m^2) and age",
    "(reference 1.8 years); both volumes scale linearly with weight.",
    "Between-subject variability is correlated across the two drugs on",
    "clearance and on volume. Residual error is combined additive plus",
    "proportional for piperacillin and proportional for tazobactam.",
    "Piperacillin uses the unsuffixed canonical compartment / parameter set;",
    "tazobactam carries the sibling-drug suffix _taz throughout."
  )
  reference <- paste(
    "Dorofaeff T, Hernandez-Mitre MP, Liu X, Guerra Valero YC, Coulthard MG,",
    "Lipman J, Wallis SC, Maresco-Pennisi DL, Parker SL, Roberts JA, Tamma PD.",
    "Population pharmacokinetics of piperacillin/tazobactam in critically",
    "ill pediatric patients. Antimicrob Agents Chemother. 2026;70(10):e00272-26.",
    "doi:10.1128/aac.00272-26."
  )
  vignette <- "Dorofaeff_2026_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    central_taz = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  dosing <- c("central", "central_taz")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of clearance (fixed exponent 0.75) and volume",
        "(fixed exponent 1) for both drugs, standardised to 70 kg and",
        "included a priori (Methods, 'Population PK modeling and",
        "validation'). Cohort median 10.80 kg, range 3-68 (Table 1)."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (AGE / 1.8) on the clearance of both drugs, exponent",
        "0.2 for piperacillin and 0.17 for tazobactam (Table 2 and the",
        "Results display equations). The reference 1.8 years is the cohort",
        "median (Table 1: median 1.8, range 0.1-11.8 years). Age entered in",
        "addition to allometric weight scaling, so it acts as a maturation",
        "term over the 1-month to 12-year range studied; it is not",
        "supported outside that range."
      ),
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (modified bedside Schwartz equation), BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power term (eGFR / 117.9) on the clearance of both drugs (Results",
        "display equations; exponents from Table 2). eGFR was computed by",
        "the modified Schwartz equation eGFR = 36.5 x height (cm) / serum",
        "creatinine (umol/L) (Methods, 'Data collection'). The reference",
        "117.9 mL/min/1.73 m^2 is the value printed in the equations; Table 1",
        "gives the cohort median as 120.5 (range 54.3-289) among the 29",
        "patients with eGFR available. 44.8% of those patients had",
        "augmented renal clearance (eGFR > 130)."
      ),
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    PELOD2 = list(
      description = "Pediatric Logistic Organ Dysfunction 2 score",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous power covariate on piperacillin and tazobactam PK parameters (Methods) but not retained. Cohort median 5, range 0-12 (Table 1). Documented for provenance only; not referenced in model()."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous power covariate (Methods) but not retained. Cohort median 28 g/L, range 19-48 (Table 1). Documented for provenance only; not referenced in model()."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a continuous power covariate (Methods) but not retained in its own right; it enters the model only through the Schwartz eGFR. Cohort median 28 umol/L, range 10-88 (Table 1). Documented for provenance only; not referenced in model()."
    ),
    MECH_VENT = list(
      description = "Mechanical ventilation (1 = ventilated)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Tested as a categorical covariate (Methods) but not retained. 19/31 (61.3%) ventilated (Table 1). Documented for provenance only; not referenced in model()."
    ),
    CONMED_INOTROPE = list(
      description = "Vasopressor / inotrope use (1 = receiving)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Tested as a categorical covariate (Methods) but not retained. 12/31 (38.7%) receiving (Table 1). Documented for provenance only; not referenced in model()."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31L,
    n_studies = 1L,
    age_range = "0.1-11.8 years (median 1.8); 11 infants (1 month to < 1 year) and 20 children",
    age_median = "1.8 years",
    weight_range = "3-68 kg (median 10.80)",
    weight_median = "10.80 kg",
    height_range = "47-165 cm (median 84)",
    sex_female_pct = 54.8,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Critically ill children receiving piperacillin-tazobactam in a",
      "paediatric intensive care unit. 48.4% admitted post-operatively",
      "(general/trauma, cardiac, neurosurgical, liver transplant and other",
      "surgery); others for respiratory distress, congenital heart disease,",
      "altered consciousness, trauma, apnoea and other causes. Median",
      "PELOD-2 5 (range 0-12); 61.3% mechanically ventilated; 38.7%",
      "receiving vasopressors or inotropes (Table 1)."
    ),
    renal_function = paste(
      "Schwartz eGFR median 120.5 mL/min/1.73 m^2 (range 54.3-289); eGFR",
      "< 90 in 34.5%, 90-130 in 20.7% and > 130 (augmented renal clearance)",
      "in 44.8% of the 29 patients with eGFR available. Serum creatinine",
      "median 28 umol/L (range 10-88). Albumin median 28 g/L (range 19-48)."
    ),
    dose_range = paste(
      "Piperacillin-tazobactam 100/12.5 mg/kg every 6 h (maximum 4 g",
      "piperacillin per dose) infused over 30 min, at the clinician's",
      "discretion. Each dose is split 8:1 into the two central compartments."
    ),
    regions = "Australia (Queensland Children's Hospital, Brisbane)",
    notes = paste(
      "Prospective open-label observational study, June 2019 to February",
      "2023. Samples before the dose and 60, 120, 240 and 360 min after the",
      "start of one infusion; 102 piperacillin and 97 tazobactam total",
      "plasma concentrations by UHPLC-MS/MS (five tazobactam samples below",
      "the 0.625 mg/L limit of quantification excluded). Estimation used",
      "Monolix 2024R1 (SAEM) and the target-attainment simulations Simulx",
      "2024R1. The model predicts TOTAL plasma concentrations; the source",
      "paper's target-attainment analysis assumes 30% protein binding for",
      "both drugs (unbound fraction 0.70) and scores 100% fT > MIC for",
      "piperacillin together with 85% fT > 2 mg/L for tazobactam."
    )
  )

  ini({
    # Typical values at 70 kg, eGFR 117.9 mL/min/1.73 m^2 and age 1.8 years.
    lcl <- log(10.31); label("Piperacillin clearance at 70 kg, eGFR 117.9, age 1.8 y (L/h)")  # Table 2 CL_pip = 10.31 L/h (RSE 5.35%; bootstrap median 10.36, 95% CI 9.28-11.4)
    lvc <- log(24.9);  label("Piperacillin volume of distribution at 70 kg (L)")              # Table 2 V_pip = 24.9 L (RSE 8.67%; bootstrap median 25.01, 95% CI 21.19-29.39)

    lcl_taz <- log(9.97);  label("Tazobactam clearance at 70 kg, eGFR 117.9, age 1.8 y (L/h)") # Table 2 CL_taz = 9.97 L/h (RSE 5.15%; bootstrap median 10, 95% CI 8.94-10.99)
    lvc_taz <- log(28.61); label("Tazobactam volume of distribution at 70 kg (L)")             # Table 2 V_taz = 28.61 L (RSE 8.03%; bootstrap median 28.53, 95% CI 24.17-33.24)

    # Allometric exponents: fixed a priori (Methods, 'Population PK modeling and validation').
    e_wt_cl     <- fixed(0.75); label("Allometric exponent of (WT/70) on piperacillin clearance (unitless)") # Methods: fixed allometric exponent 0.75 for clearance
    e_wt_vc     <- fixed(1);    label("Allometric exponent of (WT/70) on piperacillin volume (unitless)")    # Methods: fixed allometric exponent 1 for volume
    e_wt_cl_taz <- fixed(0.75); label("Allometric exponent of (WT/70) on tazobactam clearance (unitless)")  # Methods: fixed allometric exponent 0.75 for clearance
    e_wt_vc_taz <- fixed(1);    label("Allometric exponent of (WT/70) on tazobactam volume (unitless)")     # Methods: fixed allometric exponent 1 for volume

    # Covariate exponents (power models; Results display equations).
    e_crcl_cl     <- 0.65; label("Power exponent of (eGFR/117.9) on piperacillin clearance (unitless)") # Table 2 eGFR_CL_pip = 0.65 (RSE 17.4%; bootstrap median 0.65, 95% CI 0.46-0.84); the display equation prints 0.76, see vignette Errata
    e_age_cl      <- 0.2;  label("Power exponent of (AGE/1.8) on piperacillin clearance (unitless)")    # Table 2 AGE_CL_pip = 0.2 (RSE 18.3%; bootstrap median 0.2, 95% CI 0.12-0.27)
    e_crcl_cl_taz <- 0.69; label("Power exponent of (eGFR/117.9) on tazobactam clearance (unitless)")  # Table 2 eGFR_CL_taz = 0.69 (RSE 16.5%; bootstrap median 0.69, 95% CI 0.52-0.88)
    e_age_cl_taz  <- 0.17; label("Power exponent of (AGE/1.8) on tazobactam clearance (unitless)")     # Table 2 AGE_CL_taz = 0.17 (RSE 19.6%; bootstrap median 0.17, 95% CI 0.095-0.25)

    # Between-subject variability. Table 2 reports omega as SDs; variance
    # = omega^2 and covariance = r * omega1 * omega2.
    # CL block: omega 0.27 (pip), 0.27 (taz), r = 0.92 -> 0.0729, 0.067068, 0.0729
    etalcl + etalcl_taz ~ c(0.0729, 0.067068, 0.0729) # Table 2 IIV_CL_pip SD 0.27, IIV_CL_taz SD 0.27, correlation CL_pip and CL_taz 0.92
    # V block: omega 0.44 (pip), 0.41 (taz), r = 0.96 -> 0.1936, 0.173184, 0.1681
    etalvc + etalvc_taz ~ c(0.1936, 0.173184, 0.1681) # Table 2 IIV_V_pip SD 0.44, IIV_V_taz SD 0.41, correlation V_pip and V_taz 0.96

    # Residual error ('combined constant-proportional' for piperacillin,
    # proportional for tazobactam; Results).
    addSd      <- 1.33; label("Piperacillin additive residual SD (mg/L)")             # Table 2 Constant_pip = 1.33 (RSE 31.1%; bootstrap median 1.29, 95% CI 0.78-1.85)
    propSd     <- 0.09; label("Piperacillin proportional residual SD (fraction)")     # Table 2 Proportional_pip = 0.09 (RSE 15.9%; bootstrap median 0.09, 95% CI 0.06-0.12)
    propSd_taz <- 0.11; label("Tazobactam proportional residual SD (fraction)")      # Table 2 Proportional_taz = 0.11 (RSE 10.0%; bootstrap median 0.11, 95% CI 0.08-0.14)
  })

  model({
    # Individual PK parameters (Results display equations):
    # CL = CLpop * (WT/70)^0.75 * (eGFR/117.9)^theta_eGFR * (AGE/1.8)^theta_AGE
    # V  = Vpop  * (WT/70)^1
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CRCL / 117.9)^e_crcl_cl * (AGE / 1.8)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    cl_taz <- exp(lcl_taz + etalcl_taz) * (WT / 70)^e_wt_cl_taz * (CRCL / 117.9)^e_crcl_cl_taz * (AGE / 1.8)^e_age_cl_taz
    vc_taz <- exp(lvc_taz + etalvc_taz) * (WT / 70)^e_wt_vc_taz

    kel     <- cl / vc
    kel_taz <- cl_taz / vc_taz

    # One-compartment IV disposition per drug; no interaction between the
    # two drugs. Each dose is split 8:1 (piperacillin:tazobactam) into the
    # two central compartments.
    d/dt(central)     <- -kel * central
    d/dt(central_taz) <- -kel_taz * central_taz

    # Total plasma concentrations; mg / L = mg/L.
    Cc     <- central / vc
    Cc_taz <- central_taz / vc_taz

    # Monolix combined1 form, SD = a + b * f (see vignette Errata).
    Cc     ~ add(addSd) + prop(propSd) + combined1()
    Cc_taz ~ prop(propSd_taz)
  })
}
