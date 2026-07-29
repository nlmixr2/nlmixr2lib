Li_2006_meropenem <- function() {
  description <- "Two-compartment population PK model for meropenem in adult patients (Li 2006), as packaged in DDMORE Foundation Model Repository entry DDMODEL00000213."
  reference <- paste(
    "Li C, Kuti JL, Nightingale CH, Nicolau DP (2006).",
    "Population pharmacokinetic analysis and dosing regimen optimization of meropenem in adult patients.",
    "J Clin Pharmacol 46(10):1171-1178.",
    "doi:10.1177/0091270006291035.",
    "DDMORE Foundation Model Repository: DDMODEL00000213.",
    sep = " "
  )
  vignette <- "Li_2006_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  ddmore_id    <- "DDMODEL00000213"
  replicate_of <- NULL

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Used in power-form effect on CL with reference 35 years (sample median per DDMORE Model_Accommodations and Li 2006 demographics).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Used in power-form effect on V1 (central volume) with reference 70 kg (sample median per DDMORE Model_Accommodations and Li 2006 demographics).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw, measured; not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Power-form effect on CL with reference 83 mL/min",
        "(sample median per DDMORE Model_Accommodations and Li 2006 demographics).",
        "Deviation from canonical CRCL: the canonical register entry is BSA-normalized",
        "(mL/min/1.73 m^2); Li 2006 reports raw measured CrCl in mL/min, and the CRCL",
        "exponent (0.62, FIXED) was estimated under that raw-mL/min parameterization.",
        "Source data column 'CLCR' is mapped to canonical 'CRCL' on input."
      ),
      source_name        = "CLCR"
    )
  )

  population <- list(
    n_subjects     = 79,
    n_studies      = 1,
    age_range      = "Not extractable from DDMORE bundle (Li 2006 PDF not on disk).",
    age_median     = "35 years (DDMORE Model_Accommodations)",
    age_sd         = "18.2 years (DDMORE Model_Accommodations; treated as lognormal sd in the bundle's simulated cohort)",
    weight_range   = "Not extractable from DDMORE bundle (Li 2006 PDF not on disk).",
    weight_median  = "70 kg (DDMORE Model_Accommodations)",
    weight_sd      = "16.1 kg (DDMORE Model_Accommodations; treated as lognormal sd in the bundle's simulated cohort)",
    sex_female_pct = "Not extractable from DDMORE bundle (Li 2006 PDF not on disk).",
    race_ethnicity = "Not extractable from DDMORE bundle (Li 2006 PDF not on disk).",
    disease_state  = "Adult patients receiving meropenem for clinical infection. Specific indication not extractable from DDMORE bundle alone; Li 2006 reports a population PK and dosing-regimen-optimization analysis in adult patients.",
    dose_range     = "500-2000 mg meropenem IV. Six dosing groups in the DDMORE simulated cohort: 500/1000/2000 mg given as 0.5-h infusions (rates 1000/2000/4000 mg/h, all corresponding to 0.5-h infusion duration) and as 3-h infusions (rates 166.7/333.3/666.7 mg/h). Note: Model_Accommodations.txt states the rate for the 2000 mg short-infusion group as 3000 mg/h, but the simulated dataset (Simulated_DatasetMeropenem.csv) uses 4000 mg/h, consistent with a 0.5-h infusion. The dataset value (4000 mg/h) is used here.",
    crcl_median    = "83 mL/min (DDMORE Model_Accommodations; raw measured CrCl, not BSA-normalized)",
    regions        = "Not extractable from DDMORE bundle (Li 2006 PDF not on disk).",
    notes          = "Demographic medians come from DDMORE Model_Accommodations.txt, which states the model 'was used as described in publication (Li et al)'. The Li 2006 PDF is not on disk under the literature tree; full demographics, study design, indication, and inclusion criteria could not be cross-checked. The bundle's simulated dataset uses lognormal AGE (median 35, sd 18.2) and lognormal WT (median 70, sd 16.1), with all CLCR fixed at 83 mL/min because no CLCR distribution was reported in the publication."
  )

  ini({
    # Structural typical values; from DDMODEL00000213 .mdl parObj, which Model_Accommodations.txt
    # certifies as Li 2006's published estimates ("Original model was used as described in publication").
    lcl <- log(14.6); label("Clearance for AGE = 35 yr, CLCR = 83 mL/min, typical individual (CL, L/h)")  # DDMODEL00000213 parObj POP_CL
    lvc <- log(10.8); label("Central volume of distribution at WT = 70 kg, typical individual (V1, L)")  # DDMODEL00000213 parObj POP_V1
    lq  <- log(18.6); label("Inter-compartmental clearance, typical individual (Q, L/h)")  # DDMODEL00000213 parObj POP_Q
    lvp <- log(12.6); label("Peripheral volume of distribution, typical individual (V2, L)")  # DDMODEL00000213 parObj POP_V2

    # Covariate effects (power form via the .mdl INDIVIDUAL_VARIABLES log-linear block)
    e_age_cl  <- -0.34; label("AGE effect on CL: exponent of (AGE / 35) (unitless)")  # DDMODEL00000213 parObj COV_CL_AGE
    # COV_CL_CLCR was FIXED at the publication value because no CLCR distribution was reported in Li 2006.
    e_crcl_cl <- fixed(0.62); label("CRCL effect on CL: exponent of (CRCL / 83), FIXED (unitless)")  # DDMODEL00000213 parObj COV_CL_CLCR (fix=true)
    e_wt_vc   <- 0.99; label("WT effect on V1: exponent of (WT / 70) (unitless)")  # DDMODEL00000213 parObj COV_V1_WT

    # Inter-individual variability; DDMORE parObj reports OMEGA as variance directly (type is var).
    etalcl ~ 0.118  # DDMODEL00000213 parObj PPV_CL (variance)
    etalvc ~ 0.143  # DDMODEL00000213 parObj PPV_V1 (variance)
    etalq  ~ 0.290  # DDMODEL00000213 parObj PPV_Q  (variance)
    etalvp ~ 0.102  # DDMODEL00000213 parObj PPV_V2 (variance)

    # Residual error: combined additive + proportional, per the NMTRAN-rendered $ERROR block
    #   W = sqrt(RUV_ADD^2 + RUV_PROP^2 * IPRED^2).
    propSd <- 0.19; label("Proportional residual error (fraction)")  # DDMODEL00000213 parObj RUV_PROP
    addSd  <- 0.47; label("Additive residual error (mg/L)")  # DDMODEL00000213 parObj RUV_ADD
  })
  model({
    # Individual PK parameters with power-form covariate effects (per DDMODEL00000213 .mdl
    # INDIVIDUAL_VARIABLES block: ln(CL) = ln(POP_CL) + COV_CL_AGE*ln(AGE/35) +
    # COV_CL_CLCR*ln(CRCL/83) + ETA_CL, etc.).
    cl <- exp(lcl + etalcl) * (AGE / 35)^e_age_cl * (CRCL / 83)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT  / 70)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants for the 2-compartment IV ODE
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system (per DDMODEL00000213 .mdl DEQ block)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation (mg / L) and combined error model
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
