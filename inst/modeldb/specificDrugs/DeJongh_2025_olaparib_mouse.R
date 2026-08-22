DeJongh_2025_olaparib_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female SCID and athymic nude).",
    "Population PK model for the PARP inhibitor olaparib after oral",
    "dosing, with first-order absorption (ka fixed at a high value for",
    "lack of absorption-phase samples), two-compartment disposition and",
    "a single linear clearance. Co-treatment with the DNA-PK inhibitor",
    "AZD7648 reduces olaparib clearance through an Emax-shaped",
    "drug-drug-interaction term driven by the AZD7648 total daily dose",
    "(CmD50 = the daily dose halving olaparib CL). Relative oral",
    "bioavailability varied substantially between source studies and is",
    "carried as a per-study scaling factor, anchored to the absolute",
    "bioavailability estimated from the intravenous / oral study S11448.",
    "Amounts are molar (umol/kg) and concentrations are uM. Parameter",
    "values from DeJongh 2025 Table 2."
  )
  reference <- paste(
    "DeJongh J, Cadogan E, Davies M, Ramos-Montoya A, Smith A,",
    "van Steeg T, Richards R. (2025).",
    "Defining preclinical efficacy with the DNAPK inhibitor AZD7648",
    "in combination with olaparib: a minimal systems",
    "pharmacokinetic-pharmacodynamic model.",
    "J Pharmacokinet Pharmacodyn 52:17.",
    "doi:10.1007/s10928-025-09962-x.",
    sep = " "
  )
  vignette <- "DeJongh_2025_azd7648_olaparib_xenograft"

  units <- list(
    time          = "h",
    dosing        = "umol/kg",
    concentration = "uM"
  )

  compartmentData <- list(
    depot       = list(analyte = "olaparib", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "olaparib", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "olaparib", units = "umol/kg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    STRAIN_NUDE = list(
      source_name        = "STRN",
      description        = paste(
        "Mouse strain indicator: 1 = Hsd:Athymic Nude-Foxn1nu, 0 =",
        "C.B-17/IcrHan(R)Hsd-Prkdcscid (SCID). Used only to split the",
        "study S1734 relative-bioavailability factor, which the source",
        "estimated separately in each strain."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Control stream (supplementary file 6) THETA(10) for STRN = 1 (SCID) and THETA(11) for STRN = 2 (nude) within STUDY 1734."
    ),
    STUDY_S1143 = list(
      source_name        = "STUDY",
      description        = "1 = animal enrolled in study S1143 (2, 5, 10, 50 or 100 mg/kg single oral dose, nude mice); 0 = reference study S11448.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Selects the study-specific relative bioavailability lfdepot_s1143."
    ),
    STUDY_S1721 = list(
      source_name        = "STUDY",
      description        = "1 = animal enrolled in study S1721 (50, 75, 100 mg/kg single and multiple oral dose, nude mice); 0 = reference study S11448.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Selects the study-specific relative bioavailability lfdepot_s1721."
    ),
    STUDY_S1734 = list(
      source_name        = "STUDY",
      description        = paste(
        "1 = animal enrolled in study S1734 (100 mg/kg single and",
        "multiple oral dose in SCID and nude mice, with and without",
        "AZD7648 co-treatment); 0 = reference study S11448. This is the",
        "only source study that informs the AZD7648-olaparib",
        "drug-drug-interaction term."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Combined with STRAIN_NUDE to select lfdepot_s1734scid or lfdepot_s1734nude."
    ),
    STUDY_S1770 = list(
      source_name        = "STUDY",
      description        = "1 = animal enrolled in study S1770 (multiple oral dose with and without AZD7648 co-treatment, SCID mice); 0 = reference study S11448.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Selects the study-specific relative bioavailability lfdepot_s1770."
    ),
    STUDY_S1816 = list(
      source_name        = "STUDY",
      description        = "1 = animal enrolled in study S1816 (multiple oral dose with AZD7648 co-treatment, SCID mice); 0 = reference study S11448.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Selects the study-specific relative bioavailability lfdepot_s1816; this is also the factor carried into the DeJongh 2025 PK-PD tumour model."
    ),
    DOSE_AZD7648_MGKGD = list(
      source_name        = "CMDDOS",
      description        = paste(
        "Total daily dose of co-administered AZD7648 in mg/kg/day.",
        "Drives the Emax-shaped reduction of olaparib clearance;",
        "0 when olaparib is given as monotherapy."
      ),
      units              = "mg/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "DeJongh 2025 Equations 7-8: CmS = 1 when DoseAZ = 0 and",
        "CmS = 1 - DoseAZ/(DoseAZ + CmD50) when DoseAZ > 0. The two",
        "branches agree at DoseAZ = 0, so the model uses the single",
        "continuous form."
      )
    )
  )

  population <- list(
    species        = "mouse (female SCID C.B-17/IcrHan(R)Hsd-Prkdcscid and female Hsd:Athymic Nude-Foxn1nu)",
    n_studies      = 6L,
    age_range      = "at least 6 weeks of age",
    weight_range   = "(not reported in the source publication)",
    sex_female_pct = 100,
    disease_state  = "non-tumour-bearing and FaDu ATM-knockout xenograft-bearing satellite PK groups",
    dose_range     = paste(
      "Olaparib 2-100 mg/kg orally as single and multiple doses, plus a",
      "20 mg/kg intravenous arm in study S11448 that anchors absolute",
      "bioavailability; given alone and with 0-150 mg/kg/day AZD7648.",
      "Doses enter the model in umol/kg (100 mg/kg = 229.8428 umol/kg,",
      "i.e. a molar mass of 435.1 g/mol, read from the AMT and DOSE_Ola",
      "columns of the supplied NONMEM dataset, supplementary file 7)."
    ),
    regions        = "preclinical (AstraZeneca, Cambridge UK)",
    notes          = paste(
      "Olaparib plasma PK pooled across six mouse studies (S11448,",
      "S1143, S1721, S1734, S1770, S1816). Below-quantification samples",
      "were excluded after confirming they did not influence the",
      "parameter estimates. Estimation was FOCEI in NONMEM 7.4.3 with",
      "ADVAN6. Condition number 271.1. Only study S1734 contained",
      "olaparib arms with and without AZD7648 co-treatment, so the",
      "drug-drug-interaction parameter CmD50 is informed by that study",
      "alone. See DeJongh 2025 Materials and methods and Table 2."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK -- DeJongh 2025 Table 2, cross-checked against the
    # theta list carried into the PK-PD control stream (supplementary
    # file 4, 10928_2025_9962_MOESM4_ESM.txt lines 102-114), which
    # reports the same estimates to five significant figures.
    #
    # Note the unit convention the source uses: CL and Q are absolute
    # (L/h) while Vc and Vp are weight-normalised (L/kg). The resulting
    # rate constants are what the model actually uses, and the scale
    # mismatch is absorbed into the per-study bioavailability factors
    # (which is why some exceed 1). This is transcribed as published.
    # ------------------------------------------------------------------
    lka  <- fixed(log(10.0)); label("Absorption rate constant (ka, 1/h)")                  # Table 2: 10.0, fixed by assumption for lack of post-dose absorption data; checked by likelihood profiling
    lcl  <- log(1.9460);      label("Clearance (CL, L/h)")                                 # Table 2 CL = 1.95 (RSE 11.6%); control stream THETA(3) = 1.9460
    lvc  <- log(1.1656);      label("Central volume of distribution (Vc, L/kg)")           # Table 2 Vc = 1.17 (RSE 16.5%); control stream THETA(4) = 1.1656
    lq   <- log(0.6676);      label("Intercompartmental clearance (Q, L/h)")               # Table 2 Q = 0.668 (RSE 12.8%); control stream THETA(5) = 0.6676
    lvp  <- log(2.6608);      label("Peripheral volume of distribution (Vp, L/kg)")        # Table 2 Vp = 2.66 (RSE 13.3%); control stream THETA(6) = 2.6608

    # ------------------------------------------------------------------
    # AZD7648 -> olaparib drug-drug interaction (DeJongh 2025 Eq. 7-8).
    # ------------------------------------------------------------------
    lcmd50 <- log(82.538); label("AZD7648 daily dose halving olaparib clearance (CmD50, mg/kg/day)")  # Table 2: 82.5 (RSE 19.7%); control stream THETA(14) = 82.538

    # ------------------------------------------------------------------
    # Per-study relative oral bioavailability. The source estimated one
    # factor per study (and, in S1734, one per strain) as a descriptive
    # inter-study scaling on F1; every stratum therefore carries an
    # explicit suffix. S11448 is the reference study: it contained both
    # intravenous and oral arms and so anchors the absolute
    # bioavailability of the whole set.
    # ------------------------------------------------------------------
    lfdepot_s11448    <- log(0.42536); label("Relative oral bioavailability, study S11448, nude mice (F1, fraction)")  # Table 2 "F1 11,448 (Nude)" = 0.425 (RSE 15.1%); control stream THETA(7) = 0.42536
    lfdepot_s1143     <- log(0.17443); label("Relative oral bioavailability, study S1143, nude mice (F1, fraction)")   # Table 2 "F1 1143 (Nude)" = 0.174 (RSE 25.2%); control stream THETA(8) = 0.17443
    lfdepot_s1721     <- log(0.48873); label("Relative oral bioavailability, study S1721, nude mice (F1, fraction)")   # Table 2 "F1 1721 (Nude)" = 0.489 (RSE 16.1%); control stream THETA(9) = 0.48873
    lfdepot_s1734scid <- log(0.59754); label("Relative oral bioavailability, study S1734, SCID mice (F1, fraction)")   # Table 2 "Absolute bio-avail F1" = 0.598 (RSE 15.5%); control stream THETA(10) = 0.59754 (labelled "F1 1734 SCID")
    lfdepot_s1734nude <- log(0.23833); label("Relative oral bioavailability, study S1734, nude mice (F1, fraction)")   # Table 2 "F1 1734 (Nude)" = 0.238 (RSE 17.1%); control stream THETA(11) = 0.23833
    lfdepot_s1770     <- log(1.4785);  label("Relative oral bioavailability, study S1770, SCID mice (F1, fraction)")   # Table 2 "F1 1770 (SCID)" = 1.48 (RSE 21.9%); control stream THETA(12) = 1.4785
    lfdepot_s1816     <- log(4.1817);  label("Relative oral bioavailability, study S1816, SCID mice (F1, fraction)")   # Table 2 "F1 1816 (SCID)" = 4.18 (RSE 36.8%); control stream THETA(13) = 4.1817

    # ------------------------------------------------------------------
    # Inter-individual variability -- one level only, on Vc, identifiable
    # because most studies had longitudinal sampling.
    # ------------------------------------------------------------------
    etalvc ~ 0.591                                                                        # Table 2: Om_1 (Vc) = 0.591 (RSE 23.9%)

    # ------------------------------------------------------------------
    # Residual error. Table 2 reports the proportional residual error as a
    # VARIANCE (0.399); the control stream confirms the scale with
    # `W = IPRED*SQRT(THETA(1))` and `$SIGMA 1 FIX`, so the residual SD is
    # sqrt(0.399) = 0.6317 (63.2% CV).
    # ------------------------------------------------------------------
    propSd <- 0.63166; label("Proportional residual error (fraction)")                     # sqrt(Table 2 proportional residual error variance 0.399)
  })

  model({
    # --- Per-study relative bioavailability -----------------------------
    # Exactly one of the study indicators is 1; when all are 0 the
    # reference study S11448 applies, so a downstream simulation that
    # supplies no study covariate gets the absolute-bioavailability
    # anchor rather than zero drug.
    isS11448 <- 1 - STUDY_S1143 - STUDY_S1721 - STUDY_S1734 -
                    STUDY_S1770 - STUDY_S1816
    fdepot <- exp(lfdepot_s11448)    * isS11448 +
              exp(lfdepot_s1143)     * STUDY_S1143 +
              exp(lfdepot_s1721)     * STUDY_S1721 +
              exp(lfdepot_s1734scid) * STUDY_S1734 * (1 - STRAIN_NUDE) +
              exp(lfdepot_s1734nude) * STUDY_S1734 * STRAIN_NUDE +
              exp(lfdepot_s1770)     * STUDY_S1770 +
              exp(lfdepot_s1816)     * STUDY_S1816

    # --- Individual structural parameters --------------------------------
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # --- AZD7648 interaction on olaparib clearance (Eq. 7-8) -------------
    cms <- 1 - DOSE_AZD7648_MGKGD / (DOSE_AZD7648_MGKGD + exp(lcmd50))

    Cc <- central / vc

    # --- Disposition (DeJongh 2025 Equations 4-6) ------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot -
                         (cl / vc) * cms * central +
                         (q / vp) * peripheral1 - (q / vc) * central
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    f(depot) <- fdepot

    Cc ~ prop(propSd)
  })
}
