Janssen_2023_doxorubicin <- function() {
  description <- "Semi-physiological enriched two-compartment population PK model for intravenous doxorubicin in pregnant cancer patients, applying gestational changes in albumin binding, glomerular filtration, hepatic plasma flow, CYP3A4 activity and body fluid volumes to the non-pregnant Joerger 2007 base model"
  reference <- paste(
    "Janssen JM, Damoiseaux D, van Hasselt JGC, Amant FCH, van Calsteren K,",
    "Beijnen JH, Huitema ADR, Dorlo TPC. Semi-physiological Enriched Population",
    "Pharmacokinetic Modelling to Predict the Effects of Pregnancy on the",
    "Pharmacokinetics of Cytotoxic Drugs. Clin Pharmacokinet. 2023;62(8):1157-1167.",
    "doi:10.1007/s40262-023-01263-1.",
    "Non-pregnant base model from Joerger M et al. Clin Pharmacokinet.",
    "2007;46:1051-1068 (Janssen 2023 reference [15], reprinted in Janssen Table 1).",
    "Gestational physiology relations from Abduljalil K et al.",
    "Clin Pharmacokinet. 2012;51:365-396 (Janssen 2023 reference [3]),",
    "as reprinted in Janssen 2023 Eqs 1-17.",
    sep = " "
  )
  vignette <- "Janssen_2023_pregnancy_cytotoxics"
  units <- list(time = "hr", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    EGA = list(
      description = "Maternal estimated gestational age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives every gestational physiology relation (Janssen 2023 Eqs 1-17).",
        "EGA = 0 is the non-pregnant anchor at which the model collapses exactly",
        "to the Joerger 2007 non-pregnant base model. Time-varying across a",
        "pregnancy but effectively fixed within a single chemotherapy cycle.",
        "Study range 15.0-36.3 weeks (Janssen 2023 Table 2).",
        sep = " "
      ),
      source_name = "EGA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22,
    n_cycles = 27,
    n_studies = 1,
    disease_state = "pregnant patients with cancer receiving doxorubicin chemotherapy",
    ega_range = "15.0-36.3 weeks (median 28.7)",
    bsa_range = "1.56-2.49 m^2 (median 1.78)",
    sex_female_pct = 100,
    dose_range = "25, 50 and 60 mg/m^2 intravenous",
    regions = "International (INCIP registry; 137 participating institutions)",
    notes = paste(
      "Evaluation cohort from Janssen 2023 Table 2 (INCIP registry, de Haan 2018",
      "Lancet Oncol / Janssen 2021). The non-pregnant structural parameters come",
      "from the Joerger 2007 base model in breast cancer patients; the pregnancy",
      "layer was evaluated, not estimated, against these 22 patients. Janssen 2023",
      "Sect. 4 notes the sparse-sampling origin of the Joerger 2007 model and that",
      "a three-compartment model may better describe the richly sampled pregnant data.",
      sep = " "
    )
  )

  ini({
    # --- Non-pregnant (EGA = 0) structural parameters -----------------------
    # All fixed: Janssen 2023 re-uses the published Joerger 2007 base model
    # without re-estimation (Janssen 2023 Sect. 2.2 "Prediction").
    lcl <- fixed(log(47.6)); label("Clearance at EGA = 0 (L/hr)")                                  # Table 1, doxorubicin column, CL
    lvc <- fixed(log(12.3)); label("Central volume at EGA = 0 (L)")                                # Table 1, doxorubicin column, V1
    lq <- fixed(log(60.3)); label("Intercompartmental clearance to peripheral1 at EGA = 0 (L/hr)") # Table 1, doxorubicin column, Q1
    lvp <- fixed(log(421)); label("Peripheral volume at EGA = 0 (L)")                              # Table 1, doxorubicin column, V2

    # --- Drug-specific disposition constants --------------------------------
    fu_ref <- fixed(0.25); label("Unbound fraction in plasma at EGA = 0 (unitless)")               # Table 1: 75% protein binding (albumin; see vignette Errata)
    f_renal <- fixed(0.05); label("Fraction of total clearance that is renal (unitless)")          # Table 1, doxorubicin column, CL_R = 5%

    # --- IIV (Table 1 reports CV% for all four structural parameters) -------
    # omega^2 = log(CV^2 + 1)
    etalcl ~ fixed(0.058753)   # Table 1: CL IIV 24.6% CV -> log(0.246^2 + 1)
    etalvc ~ fixed(0.013828)   # Table 1: V1 IIV 11.8% CV -> log(0.118^2 + 1)
    etalq ~ fixed(0.041955)    # Table 1: Q1 IIV 20.7% CV -> log(0.207^2 + 1)
    etalvp ~ fixed(0.060625)   # Table 1: V2 IIV 25.0% CV -> log(0.250^2 + 1)

    # --- Residual error -----------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; FIXED AT ZERO - not reported)")  # Janssen 2023 reports predictions, never a residual-error model
  })

  model({
    # ---- System constants --------------------------------------------------
    qhblood <- 109 # hepatic blood flow (L/hr), fixed to the non-pregnant value  # Sect. 2.1.2, from Nakai 2002 (reference [11])

    # ---- Gestational physiology (EGA in weeks; EGA = 0 = non-pregnant) -----
    calb <- 45.8 + -0.177 * EGA + -0.0033 * EGA^2                          # Eq 1, serum albumin (g/L)
    gfr <- 114 + 3.236 * EGA + -0.0572 * EGA^2                             # Eq 7, GFR (mL/min)
    hct <- 39.1 + -0.054 * EGA + -0.0021 * EGA^2                           # Eq 8; see vignette Errata - printed coefficient is -0.0098, which is incompatible with Tables 3-4
    cyp3a4 <- 100 + 2.9826 * EGA + -0.0741 * EGA^2                         # Eq 13, CYP3A4 activity (% of non-pregnant)
    vplasma <- 2.5 + -0.0223 * EGA + 0.0042 * EGA^2 + -0.00007 * EGA^3     # Eq 15, plasma volume (L)
    ecw <- 11.86 + 0.0187 * EGA + 0.0016 * EGA^2                           # Eq 17, extracellular water (L)
    tbw <- 31.67 + 0.275 * EGA + 0.0024 * EGA^2                            # Eq 16, total body water (L)
    qhp <- (1 - hct / 100) * qhblood                                       # Eq 9, hepatic plasma flow (L/hr)

    # Non-pregnant anchors: the same polynomials evaluated at EGA = 0
    calb0 <- 45.8
    gfr0 <- 114
    qhp0 <- (1 - 39.1 / 100) * qhblood
    vplasma0 <- 2.5
    ecw0 <- 11.86
    tbw0 <- 31.67

    # ---- Protein binding (albumin) -----------------------------------------
    kd <- calb0 * fu_ref / (1 - fu_ref)                                    # Eq 3, dissociation constant, constant over pregnancy
    fu <- 1 / (1 + calb / kd)                                              # Eq 4, unbound fraction at EGA

    # ---- Non-pregnant parameter values -------------------------------------
    cl0 <- exp(lcl + etalcl)
    vc0 <- exp(lvc + etalvc)
    q0 <- exp(lq + etalq)
    vp0 <- exp(lvp + etalvp)

    # ---- Clearance during pregnancy (Eqs 5-13) -----------------------------
    clr0 <- cl0 * f_renal
    clh0 <- cl0 - clr0
    clint0 <- -(clh0 * qhp0) / (clh0 * fu_ref - qhp0 * fu_ref)             # Eq 11, well-stirred rearrangement
    clint <- clint0 * cyp3a4 / 100                                         # Eq 12 with Eq 13 (CYP3A4 only; doxorubicin Table 1)
    clh <- qhp * clint * fu / (qhp + clint * fu)                           # Eq 10, well-stirred liver model
    clr <- clr0 * (gfr / gfr0) * (fu / fu_ref)                             # Eq 6
    cl <- clr + clh                                                        # Eq 5

    # ---- Volumes during pregnancy (Eq 14, nested fluid shells) ------------
    # fu/ft is back-calculated once at EGA = 0 and held constant (Sect. 2.1.3).
    # The central compartment scales between plasma volume and ECW. The sole
    # peripheral of this two-compartment model spans the whole extravascular
    # space and scales between plasma volume and TBW (see vignette Errata).
    vc <- vplasma + (ecw - vplasma) * ((vc0 - vplasma0) / (ecw0 - vplasma0))
    vp <- vplasma + (tbw - vplasma) * ((vp0 - vplasma0) / (tbw0 - vplasma0))
    q <- q0                                                                # Q is not scaled: no gestational relation is given

    # ---- Micro-constants and ODEs ------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
