Janssen_2023_docetaxel <- function() {
  description <- "Semi-physiological enriched three-compartment population PK model for intravenous docetaxel in pregnant cancer patients, applying gestational changes in AAG binding, glomerular filtration, hepatic plasma flow, CYP3A4 activity and body fluid volumes to the non-pregnant Koolen 2010 base model"
  reference <- paste(
    "Janssen JM, Damoiseaux D, van Hasselt JGC, Amant FCH, van Calsteren K,",
    "Beijnen JH, Huitema ADR, Dorlo TPC. Semi-physiological Enriched Population",
    "Pharmacokinetic Modelling to Predict the Effects of Pregnancy on the",
    "Pharmacokinetics of Cytotoxic Drugs. Clin Pharmacokinet. 2023;62(8):1157-1167.",
    "doi:10.1007/s40262-023-01263-1.",
    "Non-pregnant base model from Koolen SLW et al. Br J Clin Pharmacol.",
    "2010;69:465-474 (Janssen 2023 reference [13], reprinted in Janssen Table 1).",
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
        "to the Koolen 2010 non-pregnant base model. Time-varying across a",
        "pregnancy but effectively fixed within a single chemotherapy cycle.",
        "Study range 26.1-35.0 weeks (Janssen 2023 Table 2).",
        sep = " "
      ),
      source_name = "EGA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 9,
    n_cycles = 10,
    n_studies = 1,
    disease_state = "pregnant patients with cancer receiving docetaxel chemotherapy",
    ega_range = "26.1-35.0 weeks (median 31.8)",
    bsa_range = "1.66-2.06 m^2 (median 1.91)",
    sex_female_pct = 100,
    dose_range = "100 mg/m^2 intravenous",
    regions = "International (INCIP registry; 137 participating institutions)",
    notes = paste(
      "Evaluation cohort from Janssen 2023 Table 2 (INCIP registry, de Haan 2018",
      "Lancet Oncol / Janssen 2021). The non-pregnant structural parameters come",
      "from the Koolen 2010 base model in patients with advanced cancer; the",
      "pregnancy layer was evaluated, not estimated, against these 9 patients.",
      sep = " "
    )
  )

  ini({
    # --- Non-pregnant (EGA = 0) structural parameters -----------------------
    # All fixed: Janssen 2023 re-uses the published Koolen 2010 base model
    # without re-estimation (Janssen 2023 Sect. 2.2 "Prediction").
    lcl <- fixed(log(44.1)); label("Clearance at EGA = 0 (L/hr)")                                  # Table 1, docetaxel column, CL
    lvc <- fixed(log(8.9)); label("Central volume at EGA = 0 (L)")                                 # Table 1, docetaxel column, V1
    lq <- fixed(log(6.1)); label("Intercompartmental clearance to peripheral1 at EGA = 0 (L/hr)")  # Table 1, docetaxel column, Q1
    lvp <- fixed(log(7.3)); label("First peripheral volume at EGA = 0 (L)")                        # Table 1, docetaxel column, V2
    lq2 <- fixed(log(14.4)); label("Intercompartmental clearance to peripheral2 at EGA = 0 (L/hr)") # Table 1, docetaxel column, Q2
    lvp2 <- fixed(log(388)); label("Second peripheral volume at EGA = 0 (L)")                      # Table 1, docetaxel column, V3

    # --- Drug-specific disposition constants --------------------------------
    fu_ref <- fixed(0.06); label("Unbound fraction in plasma at EGA = 0 (unitless)")               # Table 1: 94% protein binding (docetaxel binds AAG)
    f_renal <- fixed(0.06); label("Fraction of total clearance that is renal (unitless)")          # Table 1, docetaxel column, CL_R = 6%

    # --- IIV (Table 1 reports CV% for CL, V1 and Q2 only) -------------------
    # omega^2 = log(CV^2 + 1)
    etalcl ~ fixed(0.082906)   # Table 1: CL  IIV 29.4% CV -> log(0.294^2 + 1)
    etalvc ~ fixed(0.133560)   # Table 1: V1  IIV 37.8% CV -> log(0.378^2 + 1)
    etalq2 ~ fixed(0.040381)   # Table 1: Q2  IIV 20.3% CV -> log(0.203^2 + 1)

    # --- Residual error -----------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; FIXED AT ZERO - not reported)")  # Janssen 2023 reports predictions, never a residual-error model
  })

  model({
    # ---- System constants --------------------------------------------------
    qhblood <- 109 # hepatic blood flow (L/hr), fixed to the non-pregnant value  # Sect. 2.1.2, from Nakai 2002 (reference [11])

    # ---- Gestational physiology (EGA in weeks; EGA = 0 = non-pregnant) -----
    caag <- 0.74 + -0.0088 * EGA + 0.0001 * EGA^2                          # Eq 2, serum AAG (g/L)
    gfr <- 114 + 3.236 * EGA + -0.0572 * EGA^2                             # Eq 7, GFR (mL/min)
    hct <- 39.1 + -0.054 * EGA + -0.0021 * EGA^2                           # Eq 8; see vignette Errata - printed coefficient is -0.0098, which is incompatible with Tables 3-4
    cyp3a4 <- 100 + 2.9826 * EGA + -0.0741 * EGA^2                         # Eq 13, CYP3A4 activity (% of non-pregnant)
    vplasma <- 2.5 + -0.0223 * EGA + 0.0042 * EGA^2 + -0.00007 * EGA^3     # Eq 15, plasma volume (L)
    ecw <- 11.86 + 0.0187 * EGA + 0.0016 * EGA^2                           # Eq 17, extracellular water (L)
    tbw <- 31.67 + 0.275 * EGA + 0.0024 * EGA^2                            # Eq 16, total body water (L)
    qhp <- (1 - hct / 100) * qhblood                                       # Eq 9, hepatic plasma flow (L/hr)

    # Non-pregnant anchors: the same polynomials evaluated at EGA = 0
    caag0 <- 0.74
    gfr0 <- 114
    qhp0 <- (1 - 39.1 / 100) * qhblood
    vplasma0 <- 2.5
    ecw0 <- 11.86
    tbw0 <- 31.67

    # ---- Protein binding (docetaxel is mainly AAG bound; Sect. 3.1) --------
    kd <- caag0 * fu_ref / (1 - fu_ref)                                    # Eq 3, dissociation constant, constant over pregnancy
    fu <- 1 / (1 + caag / kd)                                              # Eq 4, unbound fraction at EGA

    # ---- Non-pregnant parameter values -------------------------------------
    cl0 <- exp(lcl + etalcl)
    vc0 <- exp(lvc + etalvc)
    q0 <- exp(lq)
    vp0 <- exp(lvp)
    q20 <- exp(lq2 + etalq2)
    vp20 <- exp(lvp2)

    # ---- Clearance during pregnancy (Eqs 5-13) -----------------------------
    clr0 <- cl0 * f_renal
    clh0 <- cl0 - clr0
    clint0 <- -(clh0 * qhp0) / (clh0 * fu_ref - qhp0 * fu_ref)             # Eq 11, well-stirred rearrangement
    clint <- clint0 * cyp3a4 / 100                                         # Eq 12 with Eq 13 (CYP3A4 only; docetaxel Table 1)
    clh <- qhp * clint * fu / (qhp + clint * fu)                           # Eq 10, well-stirred liver model
    clr <- clr0 * (gfr / gfr0) * (fu / fu_ref)                             # Eq 6
    cl <- clr + clh                                                        # Eq 5

    # ---- Volumes during pregnancy (Eq 14, nested fluid shells) ------------
    # fu/ft is back-calculated once at EGA = 0 and held constant (Sect. 2.1.3).
    # Central and shallow peripheral scale between plasma volume and ECW; the
    # deep peripheral scales between ECW and TBW (see vignette Errata).
    vc <- vplasma + (ecw - vplasma) * ((vc0 - vplasma0) / (ecw0 - vplasma0))
    vp <- vplasma + (ecw - vplasma) * ((vp0 - vplasma0) / (ecw0 - vplasma0))
    vp2 <- ecw + (tbw - ecw) * ((vp20 - ecw0) / (tbw0 - ecw0))
    q <- q0                                                                # Q is not scaled: no gestational relation is given
    q2 <- q20

    # ---- Micro-constants and ODEs ------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
