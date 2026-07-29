Wahlby_2004_paclitaxel_myelosuppression <- function() {
  description <- "Semi-mechanistic Friberg-Karlsson myelosuppression PD model for paclitaxel-induced neutropenia in 45 cancer patients, demonstrating Wahlby 2004's extended covariate-model formulation. Final-model adds a time-varying bilirubin (TBILI) effect on mean transit time and a per-occasion delta-from-baseline-bilirubin effect on the linear drug-effect Slope, with inter-individual variability in the delta-bilirubin-Slope coefficient (Wahlby 2004 Eq 3 demonstrated). Paclitaxel PK is supplied via per-subject empirical-Bayes columns (CL_INDIV, VC_INDIV, VP_INDIV) following the Friberg 2002 paclitaxel convention; users can also pair this PD model with the Friberg_2002_paclitaxel PK structure directly via the modellib registry."
  reference <- paste(
    "Wahlby U, Thomson AH, Milligan PA, Karlsson MO.",
    "Models for time-varying covariates in population pharmacokinetic-pharmacodynamic analysis.",
    "Br J Clin Pharmacol 2004;58(4):367-377.",
    "doi:10.1111/j.1365-2125.2004.02170.x.",
    "Underlying Friberg-Karlsson semi-mechanistic myelosuppression model from Friberg LE, Henningsson A, Maas H, Nguyen L, Karlsson MO.",
    "Model of chemotherapy-induced myelosuppression with parameter consistency across drugs.",
    "J Clin Oncol 2002;20(24):4713-4719; see modellib('Friberg_2002_paclitaxel').",
    "Encodes Wahlby 2004 Table 8 Final-Model column.",
    sep = " "
  )
  vignette <- "Wahlby_2004_time_varying_covariates"
  units    <- list(time = "hour", dosing = "umol", concentration = "umol/L", neutrophils = "10^9/L")

  covariateData <- list(
    CL_INDIV = list(
      description        = "Per-subject empirical-Bayes paclitaxel clearance",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Individual paclitaxel CL EBE supplied per-subject as a data column following the Friberg 2002 paclitaxel convention. Reference values: median ~285 L/h, range ~160-540 L/h.",
      source_name        = "CLI"
    ),
    VC_INDIV = list(
      description        = "Per-subject empirical-Bayes paclitaxel central volume of distribution",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Individual paclitaxel V1 EBE supplied per-subject as a data column. Reference values: median ~290 L.",
      source_name        = "V1I"
    ),
    VP_INDIV = list(
      description        = "Per-subject empirical-Bayes paclitaxel peripheral volume of distribution",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Individual paclitaxel V2 EBE supplied per-subject as a data column. Reference values: median ~995 L.",
      source_name        = "V2I"
    ),
    TBILI = list(
      description        = "Total serum bilirubin, time-varying.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Table 4: mean 8.8, median 6, range 2-41 umol/L. The BIL effect on MTT is centered at the BIL median = 6 umol/L per Wahlby 2004 Eq 1 convention.",
      source_name        = "BIL"
    ),
    TBILI_BASE = list(
      description        = "Per-subject baseline total serum bilirubin (BBIL in the source paper), time-fixed.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-subject constant. Median across 45 subjects = 6 umol/L (Table 4, BBIL row). Used together with TBILI to compute the within-subject delta (paper's DBIL term) inside model() as (TBILI - TBILI_BASE).",
      source_name        = "BBIL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adult cancer patients receiving paclitaxel chemotherapy as a 3-hour intravenous infusion every third week. Wahlby 2004 re-analyses a previously-described cohort.",
    dose_range     = NA_character_,
    n_observations = 530L,
    n_courses      = 196L,
    follow_up      = "One to 18 courses per patient (median three).",
    regions        = NA_character_,
    notes          = "One individual with a large and constant DBIL trend over time was identified as influential (Sadray, Jonsson, Karlsson - Pharm Res 1999;16(8):1260-9) and was excluded from the BIL-relationship analyses; the final-model parameter estimates in Table 8 include all 45 individuals."
  )

  ini({
    # Structural typical-value parameters (Wahlby 2004 Table 8 Final-Model column)
    lcirc0 <- log(5.47); label("Baseline circulating neutrophil count Circ0 (10^9 cells/L)")  # Table 8: theta_Circ0 = 5.47 (RSE 4.4%)
    lmtt   <- log(131);  label("Mean transit time MTT through proliferation->circulation (h)") # Table 8: theta_MTT = 131 h (RSE 3.6%)
    lslope <- log(48.6); label("Linear drug-effect Slope on paclitaxel concentration (L/umol)") # Table 8: theta_Slope = 48.6 (RSE 8%)
    gamma  <- 0.22;      label("Feedback exponent gamma on (Circ0 / circ) (unitless)")          # Table 8: theta_gamma = 0.22 (RSE 8.9%)

    # Covariate effects (Wahlby 2004 Eq 1 form for BIL-MTT; Eq 3 form for DBIL-Slope with IIV)
    e_bil_mtt    <- -0.012; label("TBILI effect on MTT (per umol/L)")          # Table 8: theta_BIL-MTT = -0.012 (RSE 13%); enters as MTT = theta_MTT * (1 + e_bil_mtt * (TBILI - 6))
    e_dbil_slope <- -0.034; label("Delta-TBILI effect on Slope (per umol/L)")  # Table 8: theta_DBIL-Slope = -0.034 (RSE 25%); enters as Slope = theta_Slope * (1 + e_dbil_slope_i * (TBILI - TBILI_BASE))

    # Inter-individual variability on structural parameters (Table 8 Final-Model column)
    etalcirc0 ~ 0.36^2  # Table 8: omega_Circ0 = 0.36 (RSE 33% rel. to variance)
    etalmtt   ~ 0.17^2  # Table 8: omega_MTT   = 0.17 (RSE 32% rel. to variance)
    etalslope ~ 0.38^2  # Table 8: omega_Slope = 0.38 (RSE 27% rel. to variance)

    # IIV in the DBIL-Slope covariate-effect coefficient (Wahlby 2004 Eq 3 demonstrated)
    etae_dbil_slope ~ 0.89^2  # Table 8: omega_DBIL-Slope = 0.89 (RSE 41% rel. to variance)

    # Residual error (Wahlby 2004 Table 8 Final-Model column)
    propSd <- 0.38; label("Proportional residual error on neutrophil count (fraction)")  # Table 8: sigma = 0.38 (RSE 9.9%)
  })

  model({
    # Individual covariate-effect coefficient (Wahlby 2004 Eq 3) for DBIL on Slope
    e_dbil_slope_i <- e_dbil_slope * exp(etae_dbil_slope)

    # Individual myelosuppression parameters with the Wahlby 2004 covariate effects.
    # MTT has a time-varying BIL effect (Eq 1, centered at BIL median = 6 umol/L);
    # Slope has a within-subject delta-BIL effect (Eq 3 with IIV).
    circ0 <- exp(lcirc0 + etalcirc0)
    mtt   <- exp(lmtt   + etalmtt) * (1 + e_bil_mtt * (TBILI - 6))
    slope <- exp(lslope + etalslope) * (1 + e_dbil_slope_i * (TBILI - TBILI_BASE))

    # PK driver: paclitaxel pharmacokinetics is supplied via per-subject EBE PK estimates
    # (CL_INDIV, VC_INDIV, VP_INDIV) carried as data columns. Intercompartmental clearance Q
    # is fixed at 204 L/h per the Friberg 2002 paclitaxel convention.
    q  <- 204
    Cc <- central / VC_INDIV  # paclitaxel central concentration (umol/L)

    # Drug effect on proliferation (linear in paclitaxel concentration) and feedback term
    edrug <- 1 - slope * Cc
    feed  <- (circ0 / circ)^gamma

    # Transit-rate constant; the Friberg-Karlsson chain has 3 transit compartments plus a
    # single proliferation compartment, so KTR = (NN + 1) / MTT = 4 / MTT.
    ktr <- 4 / mtt

    # Paclitaxel 2-compartment IV PK driven by the per-subject EBE columns
    d/dt(central)     <- -q / VC_INDIV * central - CL_INDIV / VC_INDIV * central + q / VP_INDIV * peripheral1
    d/dt(peripheral1) <-  q / VC_INDIV * central - q / VP_INDIV * peripheral1
    # Friberg-Karlsson myelosuppression chain
    d/dt(circ)        <-  ktr * precursor4 - ktr * circ
    d/dt(precursor1)  <-  ktr * precursor1 * edrug * feed - ktr * precursor1
    d/dt(precursor2)  <-  ktr * precursor1 - ktr * precursor2
    d/dt(precursor3)  <-  ktr * precursor2 - ktr * precursor3
    d/dt(precursor4)  <-  ktr * precursor3 - ktr * precursor4

    # Initial conditions: all five myelosuppression compartments start at the individual
    # baseline neutrophil count.
    circ(0)       <- circ0
    precursor1(0) <- circ0
    precursor2(0) <- circ0
    precursor3(0) <- circ0
    precursor4(0) <- circ0

    # Observation: absolute neutrophil count (10^9/L)
    ANC <- circ
    ANC ~ prop(propSd)
  })
}
