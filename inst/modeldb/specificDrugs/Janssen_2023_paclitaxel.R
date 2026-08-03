Janssen_2023_paclitaxel <- function() {
  description <- "Semi-physiological enriched three-compartment population PK model for intravenous paclitaxel with saturable elimination and saturable distribution in pregnant cancer patients, applying gestational changes in albumin binding, glomerular filtration, hepatic plasma flow, CYP3A4 activity and body fluid volumes to the non-pregnant Crombag 2019 base model"
  reference <- paste(
    "Janssen JM, Damoiseaux D, van Hasselt JGC, Amant FCH, van Calsteren K,",
    "Beijnen JH, Huitema ADR, Dorlo TPC. Semi-physiological Enriched Population",
    "Pharmacokinetic Modelling to Predict the Effects of Pregnancy on the",
    "Pharmacokinetics of Cytotoxic Drugs. Clin Pharmacokinet. 2023;62(8):1157-1167.",
    "doi:10.1007/s40262-023-01263-1.",
    "Non-pregnant base model from Crombag MRBS et al. Pharm Res.",
    "2019;36:33 (Janssen 2023 reference [14], reprinted in Janssen Table 1).",
    "Gestational physiology relations from Abduljalil K et al.",
    "Clin Pharmacokinet. 2012;51:365-396 (Janssen 2023 reference [3]),",
    "as reprinted in Janssen 2023 Eqs 1-17.",
    sep = " "
  )
  vignette <- "Janssen_2023_pregnancy_cytotoxics"
  units <- list(time = "hr", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    EGA = list(
      description = "Maternal estimated gestational age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives every gestational physiology relation (Janssen 2023 Eqs 1-17).",
        "EGA = 0 is the non-pregnant anchor at which the model collapses exactly",
        "to the Crombag 2019 non-pregnant base model. Time-varying across a",
        "pregnancy but effectively fixed within a single chemotherapy cycle.",
        "Study range 16.7-35.7 weeks (Janssen 2023 Table 2).",
        sep = " "
      ),
      source_name = "EGA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20,
    n_cycles = 25,
    n_studies = 1,
    disease_state = "pregnant patients with cancer receiving paclitaxel chemotherapy",
    ega_range = "16.7-35.7 weeks (median 31.0)",
    bsa_range = "1.74-2.27 m^2 (median 1.92)",
    sex_female_pct = 100,
    dose_range = "60, 80 and 175 mg/m^2 intravenous",
    regions = "International (INCIP registry; 137 participating institutions)",
    notes = paste(
      "Evaluation cohort from Janssen 2023 Table 2 (INCIP registry, de Haan 2018",
      "Lancet Oncol / Janssen 2021). The non-pregnant structural parameters come",
      "from the Crombag 2019 base model; the pregnancy layer was evaluated, not",
      "estimated, against these 20 patients. Doses are expressed in umol: the",
      "Crombag 2019 saturable-kinetics parameters are in umol/hr and umol/L, so a",
      "mg dose must be converted with the paclitaxel molar mass before use.",
      sep = " "
    )
  )

  ini({
    # --- Non-pregnant (EGA = 0) structural parameters -----------------------
    # All fixed: Janssen 2023 re-uses the published Crombag 2019 base model
    # without re-estimation (Janssen 2023 Sect. 2.2 "Prediction").
    lvc <- fixed(log(12)); label("Central volume at EGA = 0 (L)")                                   # Table 1, paclitaxel column, V1
    lvmax_el <- fixed(log(33.8)); label("Maximum elimination rate VM_EL at EGA = 0 (umol/hr)")      # Table 1, paclitaxel column, VM_EL
    lkm_el <- fixed(log(0.44)); label("Concentration at half VM_EL, KM_EL (umol/L)")                # Table 1, paclitaxel column, KM_EL
    lvmax_tr <- fixed(log(177)); label("Maximum transport rate to peripheral1, VM_TR (umol/hr)")    # Table 1, paclitaxel column, VM_TR
    lkm_tr <- fixed(log(1.61)); label("Concentration at half VM_TR, KM_TR (umol/L)")                # Table 1, paclitaxel column, KM_TR
    lk21 <- fixed(log(1.21)); label("First-order rate constant peripheral1 to central (1/hr)")      # Table 1, paclitaxel column, K21
    lq2 <- fixed(log(16.8)); label("Intercompartmental clearance to peripheral2 at EGA = 0 (L/hr)") # Table 1, paclitaxel column, Q2
    lvp2 <- fixed(log(268)); label("Second peripheral volume at EGA = 0 (L)")                       # Table 1, paclitaxel column, V3

    # --- Drug-specific disposition constants --------------------------------
    fu_ref <- fixed(0.05); label("Unbound fraction in plasma at EGA = 0 (unitless)")                # Table 1: 95% protein binding (paclitaxel binds albumin; Sect. 3.2)
    f_renal <- fixed(0.06); label("Fraction of total elimination that is renal (unitless)")         # Table 1, paclitaxel column, CL_R = 6%

    # --- IIV (Table 1 reports CV% for V3, VM_EL, VM_TR, KM_TR and Q2) ------
    # omega^2 = log(CV^2 + 1)
    etalvp2 ~ fixed(0.123141)      # Table 1: V3    IIV 36.2% CV -> log(0.362^2 + 1)
    etalvmax_el ~ fixed(0.070363)  # Table 1: VM_EL IIV 27.0% CV -> log(0.270^2 + 1)
    etalvmax_tr ~ fixed(0.068860)  # Table 1: VM_TR IIV 26.7% CV -> log(0.267^2 + 1)
    etalkm_tr ~ fixed(0.361545)    # Table 1: KM_TR IIV 66.0% CV -> log(0.660^2 + 1)
    etalq2 ~ fixed(0.219101)       # Table 1: Q2    IIV 49.5% CV -> log(0.495^2 + 1)

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

    # ---- Protein binding (paclitaxel binds albumin; Sect. 3.2) -------------
    kd <- calb0 * fu_ref / (1 - fu_ref)                                    # Eq 3, dissociation constant, constant over pregnancy
    fu <- 1 / (1 + calb / kd)                                              # Eq 4, unbound fraction at EGA

    # ---- Non-pregnant parameter values -------------------------------------
    vmax_el0 <- exp(lvmax_el + etalvmax_el)
    km_el <- exp(lkm_el)
    vmax_tr <- exp(lvmax_tr + etalvmax_tr)                                 # VM_TR is not scaled: Table 3 reports a gestational change for VM_EL only
    km_tr <- exp(lkm_tr + etalkm_tr)
    k21 <- exp(lk21)
    vc0 <- exp(lvc)
    q20 <- exp(lq2 + etalq2)
    vp20 <- exp(lvp2 + etalvp2)

    # ---- VM_EL during pregnancy (Sect. 3.2: "VM_EL was scaled according to
    # Eqs. 10-13"). The numeric value of VM_EL (umol/hr) is substituted for
    # CL (L/hr) in the Eq 5-13 clearance cascade. This is dimensionally
    # inconsistent but is what reproduces Table 3 (see vignette Errata).
    clr0 <- vmax_el0 * f_renal
    clh0 <- vmax_el0 - clr0
    clint0 <- -(clh0 * qhp0) / (clh0 * fu_ref - qhp0 * fu_ref)             # Eq 11, well-stirred rearrangement
    clint <- clint0 * cyp3a4 / 100                                         # Eq 12 with Eq 13 (CYP3A4 only; CYP2C8 held at the non-pregnant value, Sect. 2.1.2)
    clh <- qhp * clint * fu / (qhp + clint * fu)                           # Eq 10, well-stirred liver model
    clr <- clr0 * (gfr / gfr0) * (fu / fu_ref)                             # Eq 6
    vmax_el <- clr + clh                                                   # Eq 5, scaled VM_EL (umol/hr)

    # ---- Volumes during pregnancy (Eq 14, nested fluid shells) ------------
    # fu/ft is back-calculated once at EGA = 0 and held constant (Sect. 2.1.3).
    # The central compartment scales between plasma volume and ECW; the deep
    # peripheral scales between ECW and TBW (see vignette Errata). The first
    # peripheral compartment is an amount compartment with no volume, so it
    # carries no gestational scaling (Table 3 reports no V2 for paclitaxel).
    vc <- vplasma + (ecw - vplasma) * ((vc0 - vplasma0) / (ecw0 - vplasma0))
    vp2 <- ecw + (tbw - ecw) * ((vp20 - ecw0) / (tbw0 - ecw0))
    q2 <- q20                                                              # Q is not scaled: no gestational relation is given

    # ---- ODEs --------------------------------------------------------------
    # Saturable elimination from central and saturable distribution into
    # peripheral1 with first-order return (Crombag 2019 / Henningsson
    # parameterisation); linear distribution to peripheral2.
    Cc <- central / vc
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -vmax_el * Cc / (km_el + Cc) -
      vmax_tr * Cc / (km_tr + Cc) + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- vmax_tr * Cc / (km_tr + Cc) - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc ~ prop(propSd)
  })
}
