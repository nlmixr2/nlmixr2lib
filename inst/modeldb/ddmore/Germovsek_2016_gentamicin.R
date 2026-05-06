Germovsek_2016_gentamicin <- function() {
  description <- "Three-compartment population PK model for gentamicin in neonates and infants (Germovsek 2016), as packaged in DDMORE Foundation Model Repository entry DDMODEL00000238."
  reference <- paste(
    "Germovsek E, Kent A, Metsvaht T, Lutsar I, Klein N, Turner MA,",
    "Sharland M, Nielsen EI, Heath PT, Standing JF (2016).",
    "Development and Evaluation of a Gentamicin Pharmacokinetic Model that Facilitates",
    "Opportunistic Gentamicin Therapeutic Drug Monitoring in Neonates and Infants.",
    "Antimicrobial Agents and Chemotherapy 60(8):4869-4877.",
    "doi:10.1128/AAC.00577-16.",
    "DDMORE Foundation Model Repository: DDMODEL00000238.",
    sep = " "
  )
  vignette <- "Germovsek_2016_gentamicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  ddmore_id    <- "DDMODEL00000238"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Allometric scaling with reference 70 kg, exponent 0.632",
        "on CL (Germovsek 2016) and 1 on V1/V2/V3, 0.75 on Q/Q2.",
        "Source data column WT was reported in grams (g); the canonical convention is kg,",
        "so the model expects WT in kg (i.e., source_g / 1000)."
      ),
      source_name        = "WT (g, divide by 1000 to obtain kg)"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the Hill-type maturation function on CL.",
        "Source paper / DDMORE bundle reported PMA in weeks; the canonical PAGE is in months,",
        "so the model converts internally as PMA_weeks = PAGE_months * 4.345.",
        "T50 = 55.4 weeks (fixed) and Hill = 3.33 (fixed) are kept on the source-paper weeks scale",
        "for traceability against the .lst final estimates."
      ),
      source_name        = "PMA (weeks; multiply by 1/4.345 to obtain canonical PAGE in months)"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives a saturable postnatal-age effect on CL via PNAF = PNA / (P50 + PNA).",
        "Source paper / DDMORE bundle reported PNA in days; the canonical PNA is in months,",
        "so the model converts internally as PNA_days = PNA_months * 30.4375.",
        "P50 = 1.70 days is kept on the source-paper days scale for traceability against the .lst.",
        "When postnatal age is < 1 day (immediately post-partum), supply a small positive value",
        "(e.g., the time fraction in days / 30.4375) rather than 0 to avoid PNAF = 0."
      ),
      source_name        = "PNA (days; multiply by 1/30.4375 to obtain canonical PNA in months)"
    ),
    CREAT = list(
      description        = "Serum creatinine (measured)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Power-form effect on CL relative to a typical PMA-dependent serum",
        "creatinine TCREA = -2.8488 * PMA_weeks + 166.48 (Cuzzolin 2006, Rudd 1983);",
        "the model recomputes TCREA internally from PAGE so users only need to supply CREAT.",
        "Source rule: when CREAT was missing in the original dataset (coded -99), the typical",
        "TCREA was substituted; users should impute missing serum creatinine with the same",
        "TCREA formula before passing CREAT to the model."
      ),
      source_name        = "CREAT"
    )
  )

  population <- list(
    n_subjects     = 205,
    n_studies      = 3,
    age_range      = "Neonates and infants; postnatal age (PNA) and postmenstrual age (PMA) ranges not extractable from the DDMORE bundle (Germovsek 2016 PDF not on disk).",
    weight_range   = "Not extractable from DDMORE bundle (Germovsek 2016 PDF not on disk). The bundle's simulated dataset uses 2.12 kg as a representative weight (preterm neonate).",
    sex_female_pct = "Not extractable from DDMORE bundle.",
    race_ethnicity = "Not extractable from DDMORE bundle.",
    disease_state  = "Neonates and infants receiving gentamicin (typical clinical indication: suspected or confirmed neonatal sepsis). Pooled across three studies: Glasgow (Thomson 1988), Uppsala (Nielsen 2009), and Estonia (unpublished).",
    dose_range     = paste(
      "Gentamicin given as short IV infusions; doses in the bundle's simulated dataset",
      "range from approximately 4-7 mg per dose (about 2-3.5 mg/kg) for a 2.12 kg neonate,",
      "delivered as 5-min infusions (rate ~ 50-90 mg/h). Real clinical regimens are once-daily",
      "to every-36-hours dosing depending on PMA and renal function."
    ),
    regions        = "United Kingdom (Glasgow), Sweden (Uppsala), Estonia.",
    notes          = paste(
      "Population description is reconstructed from the .mod $PROBLEM/$DATA comments",
      "(`data from Nielsen2009, Thomson1988, Estonia_unpub`) plus the PubMed abstract of",
      "Germovsek 2016 AAC (PMID 27270281), which states 1,325 concentrations from 205 patients",
      "in the meta-analysis phase. The full PDF is not on disk under the literature tree;",
      "detailed demographics (age range, weight range, sex distribution, race) could not be",
      "cross-checked. The bundle's simulated dataset uses median covariates (WT = 2.12 kg,",
      "PMA = 33 weeks, PNA = 5.4 days, CREAT = 78 umol/L) from the original analysis."
    )
  )

  ini({
    # Structural typical values: from FINAL PARAMETER ESTIMATE block of
    # Output_real_run35b.lst (DDMODEL00000238), THETA vector. Reference subject:
    # WT = 70 kg, PMA -> infinity (fully mature), CREAT = TCREA at the subject's PMA.
    lcl  <- log(6.21);    label("Typical clearance at WT = 70 kg, fully mature, typical SCr (CL, L/h)")  # .lst TH1 = 6.21
    lvc  <- log(26.5);    label("Typical central volume at WT = 70 kg (V1, L)")  # .lst TH2 = 26.5
    lq   <- log(2.15);    label("Typical inter-compartmental CL to peripheral1 at WT = 70 kg (Q, L/h)")  # .lst TH3 = 2.15
    lvp  <- log(21.2);    label("Typical first peripheral volume at WT = 70 kg (V2, L)")  # .lst TH4 = 21.2
    lq2  <- log(0.271);   label("Typical inter-compartmental CL to peripheral2 at WT = 70 kg (Q2, L/h)")  # .lst TH5 = 0.271
    lvp2 <- log(148);     label("Typical second peripheral volume at WT = 70 kg (V3, L)")  # .lst TH6 = 148

    # Maturation parameters for CL (Hill-type function of postmenstrual age, weeks).
    # T50 and Hill are FIXED in the source ($THETA 55.4 FIX, 3.33 FIX); values are
    # carried over with fixed() so estimation will not perturb them.
    t50_pma  <- fixed(55.4);  label("Postmenstrual-age maturation T50 (weeks PMA)")  # .lst TH7 = 55.4 (FIXED)
    hill_pma <- fixed(3.33);  label("Postmenstrual-age maturation Hill exponent (unitless)")  # .lst TH8 = 3.33 (FIXED)

    # Postnatal-age (PNA) saturable effect on CL: PNAF = PNA_days / (P50 + PNA_days)
    p50_pna  <- 1.70;  label("Postnatal-age half-saturation for residual postnatal CL maturation (days)")  # .lst TH10 = 1.70

    # Power-form serum-creatinine effect on CL: OF = (CREAT / TCREA)^e_creat_cl
    e_creat_cl <- -0.130; label("Serum-creatinine effect on CL: exponent of (CREAT / TCREA), unitless")  # .lst TH9 = -0.130

    # Allometric exponents (hard-coded in the source $PK; not estimated)
    e_wt_cl  <- fixed(0.632); label("Body-weight effect on CL: exponent of (WT/70) (unitless)")  # .mod $PK TVCL line
    e_wt_q   <- fixed(0.75);  label("Body-weight effect on Q (and Q2): exponent of (WT/70) (unitless)")  # .mod $PK TVQ / TVQ2 lines
    e_wt_vc  <- fixed(1);     label("Body-weight effect on V1 (and V2/V3): exponent of (WT/70) (unitless)")  # .mod $PK TVV1 / TVV2 / TVV3 lines

    # Inter-individual variability. NONMEM $OMEGA reports variance on the
    # internal log-normal scale; final OMEGA values are read from the OMEGA
    # block of FINAL PARAMETER ESTIMATE in Output_real_run35b.lst.
    # BLOCK(2) on (CL, V1) -> correlated etalcl, etalvc
    etalcl + etalvc ~ c(0.175, 0.116, 0.112)  # .lst OMEGA(1,1)=0.175; OMEGA(2,1)=0.116; OMEGA(2,2)=0.112
    # No IIV on Q (ETA(3) FIXED at 0 in source)
    etalvp  ~ 0.132   # .lst OMEGA(4,4) = 0.132 (V2)
    # No IIV on Q2 (ETA(5) FIXED at 0 in source)
    etalvp2 ~ 0.177   # .lst OMEGA(6,6) = 0.177 (V3)

    # Residual error: combined proportional + additive on the linear scale.
    # NONMEM $SIGMA reports variances; nlmixr2 expects standard deviations.
    propSd <- sqrt(0.0360); label("Proportional residual error (fraction; SD scale)")  # .lst SIGMA(1,1) = 0.0360
    addSd  <- sqrt(0.0164); label("Additive residual error (mg/L; SD scale)")          # .lst SIGMA(2,2) = 0.0164
  })
  model({
    # 1. Covariate transforms back to source-paper units (weeks, days)
    pma_weeks <- PAGE * 4.345
    pna_days  <- PNA  * 30.4375

    # 2. Hill-type postmenstrual-age maturation factor on CL (approaches 1 at full maturity)
    mf <- pma_weeks^hill_pma / (pma_weeks^hill_pma + t50_pma^hill_pma)

    # 3. Saturable postnatal-age factor on CL (approaches 1 as PNA grows beyond P50)
    pnaf <- pna_days / (p50_pna + pna_days)

    # 4. Typical PMA-dependent serum creatinine (umol/L) per Cuzzolin 2006 and
    #    Rudd 1983 (TCREA = -2.8488 * PMA_weeks + 166.48); used as the reference
    #    in the power-form CREAT effect on CL.
    tcrea <- -2.8488 * pma_weeks + 166.48
    of    <- (CREAT / tcrea)^e_creat_cl

    # 5. Individual PK parameters with allometric and physiological covariate effects
    cl <- exp(lcl + etalcl) * mf * (WT / 70)^e_wt_cl * pnaf * of
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)           * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    q2 <- exp(lq2)          * (WT / 70)^e_wt_q
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc

    # 6. Micro-constants for the 3-compartment IV ODE
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # 7. ODE system (3-compartment IV; gentamicin given as short infusion to central)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # 8. Observation (mg/L; dose in mg, volume in L) and combined residual error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
