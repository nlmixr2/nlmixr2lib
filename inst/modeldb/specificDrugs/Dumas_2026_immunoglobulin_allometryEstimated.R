Dumas_2026_immunoglobulin_allometryEstimated <- function() {
  description <- "Two-compartment population PK model for total serum IgG after intravenous immune globulin 10% (BIVIGAM) in children (2-16 years) and adults with primary immunodeficiency (Dumas 2026), the variant with the allometric exponents ESTIMATED (0.939 shared by CL and Q, 0.528 shared by Vc and Vp; reference weight 70 kg). Linear elimination with no endogenous IgG production term: both compartments are initialised at the subject's pre-treatment total IgG (assumed steady state on prior replacement therapy). Log-normal IIV on CL, Vc and Vp; proportional residual error whose magnitude depends on whether the IgG result came from a dedicated PK sample or from routine central-laboratory monitoring. The companion Dumas_2026_immunoglobulin_allometryFixed model fixes the exponents to 0.75 and 1; this estimated-exponent model has the lower OFV (delta OFV -50.398 for 2 extra parameters) and is the preferred model."
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
      notes = "Baseline (time-fixed) body weight; power scaling (WT/70)^0.939 on CL and Q and (WT/70)^0.528 on Vc and Vp (Dumas 2026 Table 2; Listing S1 $PK). Cohort range 16.4-135 kg (Table 1).",
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
    # 'Estimated allometry model'. Listing S1 carries the same values in
    # dL/h and dL (0.0404 dL/h = 0.0969 L/day; 32.4 dL; 45.5 dL; 0.18 dL/h).
    lcl <- log(0.0969); label("Clearance for a 70 kg subject (L/day)")                      # Table 2: CL = 0.0969 L/day (RSE 3.13%)
    lvc <- log(3.25);   label("Central volume of distribution for a 70 kg subject (L)")     # Table 2: Vc = 3.25 L (RSE 2.59%)
    lvp <- log(4.55);   label("Peripheral volume of distribution for a 70 kg subject (L)")  # Table 2: Vp = 4.55 L (RSE 21.3%)
    lq  <- log(0.432);  label("Intercompartmental clearance for a 70 kg subject (L/day)")   # Table 2: Q = 0.432 L/day (RSE 9.63%)

    # Allometric exponents. Listing S1: THETA(7) is shared by CL and Q,
    # THETA(8) is shared by V1 and V2.
    e_wt_cl_q  <- 0.939; label("Allometric exponent on CL and Q (unitless)")    # Table 2: 'Allometric clearance' = 0.939 (RSE 6.69%)
    e_wt_vc_vp <- 0.528; label("Allometric exponent on Vc and Vp (unitless)")   # Table 2: 'Allometric volume' = 0.528 (RSE 11.3%)

    # IIV. Table 2 reports %CV; omega^2 = log(1 + CV^2) reproduces the
    # Listing S1 $OMEGA values (0.0549, 0.0308, 2.37). No IIV on Q
    # (Listing S1 $OMEGA 0 FIX on ETA(4)).
    etalcl ~ 0.05469  # Table 2: IIV on CL 23.7 %CV; log(1 + 0.237^2)
    etalvc ~ 0.03085  # Table 2: IIV on Vc 17.7 %CV; log(1 + 0.177^2)
    etalvp ~ 2.37339  # Table 2: IIV on Vp 312 %CV; log(1 + 3.12^2)

    # Residual error: proportional, magnitude by data source (Listing S1
    # $ERROR: W = THETA * IPRED with $SIGMA 1 FIX, so THETA is the SD).
    propSd    <- 0.119;  label("Proportional residual error, dedicated PK samples (fraction)")           # Table 2: Proportional error (SDTM.PC) = 0.119
    propSd_lb <- 0.0808; label("Proportional residual error, routine central-laboratory samples (fraction)") # Table 2: Proportional error (SDTM.LB) = 0.0808
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq) * (WT / 70)^e_wt_cl_q

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
