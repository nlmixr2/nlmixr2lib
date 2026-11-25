Dumas_2026_immunoglobulin_allometryFixed <- function() {
  description <- "Two-compartment population PK model for total serum IgG after intravenous immune globulin 10% (BIVIGAM) in children (2-16 years) and adults with primary immunodeficiency (Dumas 2026), the fixed-allometry sensitivity variant with the standard allometric exponents FIXED to 0.75 on CL and 1 on Vc (reference weight 70 kg; Q and Vp are not weight-scaled), which the authors qualified and used for simulation alongside the preferred estimated-exponent model. Linear elimination with no endogenous IgG production term: both compartments are initialised at the subject's pre-treatment total IgG (assumed steady state on prior replacement therapy). Log-normal IIV on CL, Vc and Vp; proportional residual error whose magnitude depends on whether the IgG result came from a dedicated PK sample or from routine central-laboratory monitoring. The companion Dumas_2026_immunoglobulin_allometryEstimated model estimates the exponents and has the lower OFV (delta OFV -50.398 for 2 extra parameters)."
  reference <- paste(
    "Dumas T, Mehr S, Avila R, Doman A, Kestenberg K, Melamed I.",
    "Pharmacokinetic Model Based Sensitivity Analysis to Lower Recruitment Burden for Young Children",
    "Requiring Intravenous Immunoglobulin G Replacement.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(6):e70275.",
    "doi:10.1002/psp4.70275",
    sep = " "
  )
  vignette <- "Dumas_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed) body weight; fixed power scaling (WT/70)^0.75 on CL and (WT/70)^1 on Vc only (Dumas 2026 Methods 'Fixed-allometry standard exponents were applied to CL and Vc'; Results 'standard exponents were applied to the central compartment'). Q and Vp are not weight-scaled in this variant, unlike the estimated-exponent model, whose Listing S1 code shares the exponents with Q and Vp. Cohort range 16.4-135 kg (Table 1).",
      source_name = "BWT"
    ),
    IGG = list(
      description = "Pre-treatment (baseline) serum total IgG concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a covariate effect: sets the initial condition of both compartments (central(0) = IGG * vc, peripheral1(0) = IGG * vp), assuming pre-dose steady-state distribution because every subject entered on stable replacement therapy (Dumas 2026 Methods; Listing S1 A_0(1) = BIGG*V1, A_0(2) = BIGG*V2). The source dataset carried mg/dL (Listing S1 'BIGG = baseline total IGG LEVEL MEAN=989 mg/dL'); divide by 100 to get g/L. Set IGG = 0 to simulate exogenous IgG alone from an empty system.",
      source_name = "BIGG"
    ),
    SAMPLE_CENTRAL_LAB = list(
      description = "Per-observation indicator: 1 = total IgG result from a routine-monitoring blood draw assayed by the clinical central laboratory (SDTM LB domain); 0 = dedicated PK sample assayed by the bioanalytical laboratory (SDTM PC domain)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dedicated PK sample, SDTM PC)",
      notes = "Switches the proportional residual-error magnitude only (Dumas 2026 Table 2 footnote: 'Residual error was estimated separately for IgG results from SDTM.PC and SDTM.LB'; Listing S1 DVID = 1 'PK draws/ba lab - pk profiles and sparse', DVID = 2 'Routine monitoring draws/central lab - sparse'). SAMPLE_CENTRAL_LAB = DVID - 1. Does not change the prediction.",
      source_name = "DVID"
    )
  )

  compartmentData <- list(
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 79L,
    n_studies = 2L,
    age_range = "3-75 years",
    age_groups = "2 to <6 years: 3; 6 to <12 years: 9; 12 to 16 years: 13; >16 years: 54",
    weight_range = "16.4-135 kg",
    sex_female_pct = round(100 * 32 / 79, 1),
    race_ethnicity = "Not reported",
    disease_state = "Primary immunodeficiency (PID) on stable intravenous IgG replacement therapy",
    dose_range = "BIVIGAM (IVIG 10%) 300-800 mg/kg IV every 3 or 4 weeks (administered 265-1100 mg/kg after trough-guided adjustment)",
    regions = "Not reported",
    notes = "Pooled Phase 3 study Nabi-7101 (NCT00538915, n = 63, about 1 year) and Phase 4 pediatric study 994 (NCT03164967, n = 16 evaluable, 5 months); 1243 total IgG observations (Dumas 2026 Table 1, Table S1). Every subject had at least 3 months of stable dosing with a trough of at least 5 g/L before entry. Baseline total IgG mean 7.95-10.40 g/L by age group (Table 1)."
  )

  ini({
    # Structural parameters for a 70 kg subject. Dumas 2026 Table 2,
    # 'Fixed allometry model'. Structure as Listing S1 (printed for the
    # estimated-exponent model), but the fixed standard exponents scale only
    # CL and Vc (Methods: 'applied to CL and Vc'); Q and Vp are unscaled.
    lcl <- log(0.0979); label("Clearance for a 70 kg subject (L/day)")                      # Table 2: CL = 0.0979 L/day (RSE 3.16%)
    lvc <- log(3.36);   label("Central volume of distribution for a 70 kg subject (L)")     # Table 2: Vc = 3.36 L (RSE 3.29%)
    lvp <- log(3.62);   label("Peripheral volume of distribution for a 70 kg subject (L)")  # Table 2: Vp = 3.62 L (RSE 19.8%)
    lq  <- log(0.4872); label("Intercompartmental clearance for a 70 kg subject (L/day)")   # Table 2: Q = 0.4872 L/day (RSE 14.9%)

    # Allometric exponents, fixed to the standard values.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")   # Table 2: 'Allometric clearance (fixed)' = 0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on Vc (unitless)")   # Table 2: 'Allometric volume (fixed)' = 1.0

    # IIV. Table 2 reports %CV; omega^2 = log(1 + CV^2), the conversion that
    # reproduces the Listing S1 $OMEGA values for the estimated-exponent
    # model. No IIV on Q (Listing S1 $OMEGA 0 FIX on ETA(4)).
    etalcl ~ 0.06157  # Table 2: IIV on CL 25.2 %CV; log(1 + 0.252^2)
    etalvc ~ 0.05876  # Table 2: IIV on Vc 24.6 %CV; log(1 + 0.246^2)
    etalvp ~ 1.73441  # Table 2: IIV on Vp 216 %CV; log(1 + 2.16^2)

    # Residual error: proportional, magnitude by data source (Listing S1
    # $ERROR: W = THETA * IPRED with $SIGMA 1 FIX, so THETA is the SD).
    propSd    <- 0.119;  label("Proportional residual error, dedicated PK samples (fraction)")             # Table 2: Proportional error (SDTM.PC) = 0.119
    propSd_lb <- 0.0820; label("Proportional residual error, routine central-laboratory samples (fraction)") # Table 2: Proportional error (SDTM.LB) = 0.0820
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # No endogenous production: the states hold total IgG, initialised at
    # the pre-treatment level in both compartments (Listing S1 A_0).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    central(0)     <- IGG * vc
    peripheral1(0) <- IGG * vp

    Cc <- central / vc
    propSd_eff <- propSd * (1 - SAMPLE_CENTRAL_LAB) + propSd_lb * SAMPLE_CENTRAL_LAB
    Cc ~ prop(propSd_eff)
  })
}
