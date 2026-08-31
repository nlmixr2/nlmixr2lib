Mukker_2026_tuvusertib_hematology <- function() {
  description <- paste(
    "Semi-mechanistic population PK/PD model of tuvusertib-induced anemia in 55",
    "patients with advanced solid tumors (DDRiver Solid Tumors 301 Part A1).",
    "Three serial four-compartment transit chains carry erythroid progenitor",
    "cells, reticulocytes and red blood cells; hemoglobin is proportional to the",
    "summed reticulocyte and red-cell counts. Tuvusertib inhibits progenitor",
    "cell production through an Emax function of plasma concentration,",
    "and two power-law negative-feedback loops on the reticulocyte and red-cell",
    "counts stabilise the system. The tuvusertib PK driving the effect is the",
    "companion Mukker 2026 population PK model, carried here with its parameters",
    "fixed because the PK/PD analysis was sequential."
  )
  reference <- paste(
    "Mukker JK, Diderichsen PM, Hellmann F, Yap TA, Plummer R, Tolcher AW,",
    "de Bono JS, Gounaris I, Szucs Z, Zimmermann A, Kareva I, Bolleddula J,",
    "Seithel-Keuth A, Locatelli G, Enderlin M, Hicking C, Zutshi A, Gao W,",
    "Strotmann R, Benincosa L, Venkatakrishnan K.",
    "Integrated Population Pharmacokinetic, Pharmacodynamic, and Safety Analyses",
    "to Inform Dosage Selection in the Clinical Development of the ATR Inhibitor",
    "Tuvusertib. Clin Pharmacol Ther. 2026;119(3):618-628. doi:10.1002/cpt.70029.",
    "Structural PK/PD framework (transit chains, feedback, SHB) inherited from",
    "Zhang X, Chua L, Ernest C II, Macias W, Rooney T, Tham LS.",
    "Dose/exposure-response modeling to support dosing recommendation for phase",
    "III development of baricitinib in patients with rheumatoid arthritis.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(12):804-813. doi:10.1002/psp4.12251",
    "(Equations 8-20 and 21-33).",
    sep = " "
  )
  vignette <- "Mukker_2026_tuvusertib"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot              = list(analyte = "tuvusertib", units = "mg", specimen = "administration site", verified = TRUE),
    central            = list(analyte = "tuvusertib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1        = list(analyte = "tuvusertib", units = "mg", specimen = "tissue", verified = TRUE),
    clearance_capacity = list(analyte = "none", units = "unitless (relative to drug-free baseline of 1)", specimen = "not applicable", verified = TRUE),
    prol               = list(analyte = "erythroid progenitor cells", units = "10^9/mL", specimen = "tissue", verified = TRUE),
    precursor1         = list(analyte = "erythroid progenitor cells", units = "10^9/mL", specimen = "tissue", verified = TRUE),
    precursor2         = list(analyte = "erythroid progenitor cells", units = "10^9/mL", specimen = "tissue", verified = TRUE),
    precursor3         = list(analyte = "erythroid progenitor cells", units = "10^9/mL", specimen = "tissue", verified = TRUE),
    reticulocytes1     = list(analyte = "reticulocytes", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    reticulocytes2     = list(analyte = "reticulocytes", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    reticulocytes3     = list(analyte = "reticulocytes", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    reticulocytes4     = list(analyte = "reticulocytes", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    erythrocytes1      = list(analyte = "red blood cells", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    erythrocytes2      = list(analyte = "red blood cells", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    erythrocytes3      = list(analyte = "red blood cells", units = "10^9/mL", specimen = "blood cell", verified = TRUE),
    erythrocytes4      = list(analyte = "red blood cells", units = "10^9/mL", specimen = "blood cell", verified = TRUE)
  )

  # Mukker 2026 Results: "No discernible effect of the evaluated covariates
  # (age, body weight, sex, race, ECOG status) on POPPK/PD parameters'
  # variability for tuvusertib could be identified."
  covariateData <- list()
  covariatesDataExcluded <- list(
    AGE = list(description = "Age", units = "y", type = "continuous",
               notes = "Screened in the stepwise PK/PD covariate search; not retained."),
    WT = list(description = "Baseline body weight", units = "kg", type = "continuous",
              notes = "Screened; not retained."),
    SEXF = list(description = "Female sex indicator", units = "(binary)", type = "binary",
                reference_category = "male", notes = "Screened; not retained."),
    ECOG = list(description = "Baseline ECOG performance status", units = "(score)", type = "categorical",
                notes = "Screened; not retained."),
    RACE_ASIAN = list(description = "Asian race indicator", units = "(binary)", type = "binary",
                      reference_category = "non-Asian",
                      notes = "Screened; not retained. Assessed post hoc by overlaying observed HGB-time profiles for Asian and Black or African American patients on the model 90% PI (Figure 4f).")
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 1L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = c(Asian = 5, `Black or African American` = 2, `non-Asian` = 48),
    disease_state  = "advanced / metastatic solid tumors",
    dose_range     = "5-270 mg once daily, plus intermittent regimens (180 mg QD 2w on/1w off, 220 mg QD 2w on/1w off, 150 mg BID 4d on/3d off)",
    regions        = "multicenter first-in-human trial (DDRiver Solid Tumors 301, NCT04170153, Part A1)",
    notes          = paste(
      "Same 55-patient cohort as the companion population PK model. Reticulocyte,",
      "red-blood-cell and hemoglobin time courses were fitted jointly, driven by",
      "the tuvusertib concentration predicted by the population PK model (a",
      "sequential PK-then-PD fit), so the PK parameters carried in this file are",
      "fixed rather than estimated. The system is started in equilibrium from the",
      "estimated baseline reticulocyte count (Zhang 2017 Equations 21-32).",
      "Anemia was the dose-limiting toxicity: 36% Grade >=3 across Part A1.",
      "Estimation in NONMEM; 250-replicate nonparametric bootstrap."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # PK parameters -- FIXED. Mukker 2026 Table 1 (final POPPK model). The
    # PK/PD analysis was sequential ("in response to the time course of
    # tuvusertib plasma concentration predicted by the POPPK model"), so these
    # are not re-estimated here.
    # ---------------------------------------------------------------------
    lka   <- fixed(log(0.441));  label("Absorption rate constant (KA, 1/h)")                        # Table 1: KA = 0.441 1/h
    lcl   <- fixed(log(55.7));   label("Apparent clearance at the drug-free baseline (CL/F, L/h)")  # Table 1: CL/F = 55.7 L/h
    lvc   <- fixed(log(30.0));   label("Apparent central volume of distribution (VC/F, L)")         # Table 1: VC/F = 30.0 L
    lq    <- fixed(log(3.59));   label("Apparent intercompartmental clearance (Q/F, L/h)")          # Table 1: Q/F = 3.59 L/h
    lvp   <- fixed(log(136));    label("Apparent peripheral volume of distribution (VP/F, L)")      # Table 1: VP/F = 136 L
    ltlag <- fixed(log(0.369));  label("Absorption lag time (ALAG1, h)")                            # Table 1: ALAG1 = 0.369 h
    lkcl  <- fixed(log(0.0878)); label("Clearance-compartment turnover rate constant (KCL, 1/h)")   # Table 1: KCL = 0.0878 1/h
    slp   <- fixed(0.303);       label("Tuvusertib effect on clearance-compartment loss (SLP, %/(ng/mL))") # Table 1: SLP = 0.303 %/(ng/mL)

    # ---------------------------------------------------------------------
    # PK/PD system parameters -- Mukker 2026 Table 2 (final POPPK/PD model).
    # Rate constants are tabulated per DAY; the model time unit is hours, so
    # model() divides each by 24. ini() keeps the published per-day value so
    # the source trace is a direct read of Table 2.
    # ---------------------------------------------------------------------
    lretbl <- log(0.0645); label("Baseline TOTAL reticulocyte count (RETBL, 10^9/mL)")           # Table 2: RETBL = 0.0645 10^9/mL (0.0583, 0.0718)
    lktr1  <- log(3.85);   label("Progenitor transit / proliferation rate constant (KTR1, 1/day)") # Table 2: KTR1 = 3.85 1/day (1.51, 67.6)
    lktr2  <- log(3.33);   label("Reticulocyte transit rate constant (KTR2, 1/day)")             # Table 2: KTR2 = 3.33 1/day (2.58, 4.17)
    lkcir  <- log(0.0557); label("Red-blood-cell circulation rate constant (KCIR, 1/day)")       # Table 2: KCIR = 0.0557 1/day (0.0463, 0.0670)
    lshb   <- log(29.4);   label("Hemoglobin proportionality factor (SHB, (g/L)/(10^9/mL))")     # Table 2: SHB = 29.4 (g/L)/(10^9/mL) (28.9, 30.0)

    # Feedback exponents. Table 2 reports GAM1 / GAM2 not as the exponent but
    # as the percent change in progenitor production for a 10% increase in the
    # respective cell count (footnotes a and b). The feedback enters as
    # FB = (RETBL/RET)^GAM1 * (RBCBL/RBC)^GAM2 (Zhang 2017 Eq. 33), so
    #   GAM = ln(1 + pct/100) / ln(1/1.1).
    #   GAM1: ln(1 - 0.00883) / ln(1/1.1) = 0.09306
    #   GAM2: ln(1 - 0.01930) / ln(1/1.1) = 0.20447
    # The published bootstrap CIs are printed in reversed order
    # (-0.0533, -2.67) because the transform is monotone decreasing in GAM,
    # which independently confirms the exponent reading.
    lgam1 <- log(0.09306); label("Reticulocyte negative-feedback exponent (GAM1, unitless)")  # Table 2: GAM1 = -0.883% per 10% rise in RET (-0.0533, -2.67) -> exponent 0.09306
    lgam2 <- log(0.20447); label("Red-blood-cell negative-feedback exponent (GAM2, unitless)") # Table 2: GAM2 = -1.93% per 10% rise in RBC (-0.129, -3.46) -> exponent 0.20447

    lemaxpd <- log(59.8); label("Maximum fractional inhibition of progenitor production (EMAXPD, %)") # Table 2: EMAXPD = 59.8% (10.9, 1040)
    lec50pd <- log(736);  label("Tuvusertib concentration at half-maximal effect (EC50PD, ng/mL)")    # Table 2: EC50PD = 736 ng/mL (69.7, 68,700)

    # ---------------------------------------------------------------------
    # Interindividual variability. Table 2 reports %CV; omega^2 = log(CV^2 + 1).
    # ---------------------------------------------------------------------
    etalretbl  ~ 0.010965; label("IIV on baseline reticulocyte count")   # Table 2: IIV RETBL = 10.5 %CV (4.00, 24.4), shrinkage 32%;   log(0.105^2 + 1) = 0.010965
    etalkcir   ~ 0.024045; label("IIV on KCIR")                          # Table 2: IIV KCIR = 15.6 %CV (10.9, 26.4), shrinkage 12.8%;  log(0.156^2 + 1) = 0.024045
    etalshb    ~ 0.003775; label("IIV on SHB")                           # Table 2: IIV SHB = 6.15 %CV (2.68, 8.79), shrinkage 19.8%;   log(0.0615^2 + 1) = 0.003775
    etalemaxpd ~ 0.148420; label("IIV on EMAXPD")                        # Table 2: IIV EMAXPD = 40.0 %CV (17.2, 88.4), shrinkage 33.5%; log(0.400^2 + 1) = 0.148420

    # PK IIV carried at the Table 1 values, fixed (sequential PK/PD fit).
    etalcl ~ fixed(0.05922); label("IIV on CL/F")   # Table 1: IIV CL/F = 24.7 %CV; log(0.247^2 + 1)
    etalvc ~ fixed(0.86243); label("IIV on VC/F")   # Table 1: IIV VC/F = 117 %CV;  log(1.17^2 + 1)

    # ---------------------------------------------------------------------
    # Residual unexplained variability -- all three PD endpoints additive.
    # ---------------------------------------------------------------------
    expSd     <- fixed(0.714); label("Residual SD on the natural-log tuvusertib concentration scale (log ng/mL)") # Table 1: SD of RUV = 0.714 log ng/mL; carried fixed with the rest of the PK model
    addSd_RET <- 0.0190; label("Additive residual SD on reticulocyte count (10^9/mL)")   # Table 2: Add. RUV RET = 0.0190 10^9/mL (0.0111, 0.0241)
    addSd_RBC <- 0.221;  label("Additive residual SD on red-blood-cell count (10^9/mL)") # Table 2: Add. RUV RBC = 0.221 10^9/mL (0.183, 0.256)
    addSd_hb  <- 6.86;   label("Additive residual SD on hemoglobin (g/L)")               # Table 2: Add. RUV HGB = 6.86 g/L (5.22, 8.13)
  })

  model({
    # --- PK parameters --------------------------------------------------------
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    tlag <- exp(ltlag)
    kcl  <- exp(lkcl)

    # --- PK/PD system parameters ---------------------------------------------
    # Table 2 tabulates KTR1 / KTR2 / KCIR per DAY; the model time unit is the
    # hour, so each is divided by 24 here.
    retbl  <- exp(lretbl + etalretbl)
    ktr1   <- exp(lktr1) / 24
    ktr2   <- exp(lktr2) / 24
    kcir   <- exp(lkcir + etalkcir) / 24
    shb    <- exp(lshb + etalshb)
    gam1   <- exp(lgam1)
    gam2   <- exp(lgam2)
    emaxpd <- exp(lemaxpd + etalemaxpd)
    ec50pd <- exp(lec50pd)

    # --- Tuvusertib PK (companion Mukker 2026 population PK model) ------------
    Cc <- 1000 * central / vc

    d/dt(clearance_capacity) <- kcl - kcl * (1 + slp / 100 * Cc) * clearance_capacity
    clearance_capacity(0)    <- 1

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          (cl / vc) * clearance_capacity * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    alag(depot) <- tlag

    # --- Baselines (Zhang 2017 Equations 21-32) -------------------------------
    # RETBL is the observed TOTAL baseline reticulocyte count, partitioned
    # equally across the four reticulocyte transit bins. The progenitor and
    # red-cell bins follow from the steady-state flux balance:
    #   progenitor bin = ret_bin * KTR2 / KTR1
    #   red-cell  bin  = ret_bin * KTR2 / KCIR
    # Cross-check performed while extracting: the implied baseline hemoglobin
    # SHB * (RETBL + RETBL*KTR2/KCIR) = 29.4 * 0.0645 * (1 + 3.33/0.0557)
    # = 115.3 g/L, which matches the ~115 g/L observed median at day 0 in the
    # Figure 5 prediction-corrected VPC. Reading RETBL as a PER-BIN value
    # instead would give ~461 g/L and is excluded.
    ret_bin  <- retbl / 4
    prol_bin <- ret_bin * ktr2 / ktr1
    rbc_bin  <- ret_bin * ktr2 / kcir
    rbcbl    <- 4 * rbc_bin

    # --- Observed totals ------------------------------------------------------
    RET <- reticulocytes1 + reticulocytes2 + reticulocytes3 + reticulocytes4
    RBC <- erythrocytes1 + erythrocytes2 + erythrocytes3 + erythrocytes4

    # --- Drug effect and negative feedback (Zhang 2017 Eq. 33) ----------------
    # EMAXPD is tabulated as a percent, so it enters as a fraction.
    edrug <- (emaxpd / 100) * Cc / (ec50pd + Cc)
    fb    <- (retbl / RET)^gam1 * (rbcbl / RBC)^gam2

    # --- Erythropoiesis cascade (Zhang 2017 Equations 8-19) -------------------
    # Three serial four-compartment chains: progenitors (prol + precursor1-3),
    # reticulocytes, red blood cells. The proliferation constant of the head
    # compartment equals KTR1 (Mukker 2026 Methods).
    #
    # The drug effect multiplies the PROLIFERATION term of the head compartment
    # only, leaving every downstream maturation flux untouched. Two main-text
    # statements say exactly this -- Methods, "inhibition of the production of
    # progenitor cells by an Emax model", and Results, "an Emax drug effect on
    # the production of the progenitor cells" -- and Figure 2b draws the
    # "Negative effect (EMAXPD; EC50PD)" arrow into B1 alone. Supplement S1
    # instead says "production and maturation of progenitor cells" once; that
    # reading is excluded because applying (1 - edrug) to the KTR1 maturation
    # fluxes as well makes cells pile up in the progenitor chain and flattens
    # the dose response, giving ~3% Grade >=3 anemia at every dose against the
    # 17.2-34.5% of Table S2 (RMSE 28.5 percentage points over the nine
    # published cells, versus 3.5 for the form encoded here). See the
    # vignette's Assumptions and deviations section.
    d/dt(prol)       <- ktr1 * prol * (1 - edrug) * fb - ktr1 * prol
    d/dt(precursor1) <- ktr1 * (prol - precursor1)
    d/dt(precursor2) <- ktr1 * (precursor1 - precursor2)
    d/dt(precursor3) <- ktr1 * (precursor2 - precursor3)

    d/dt(reticulocytes1) <- ktr1 * precursor3 - ktr2 * reticulocytes1
    d/dt(reticulocytes2) <- ktr2 * (reticulocytes1 - reticulocytes2)
    d/dt(reticulocytes3) <- ktr2 * (reticulocytes2 - reticulocytes3)
    d/dt(reticulocytes4) <- ktr2 * (reticulocytes3 - reticulocytes4)

    d/dt(erythrocytes1) <- ktr2 * reticulocytes4 - kcir * erythrocytes1
    d/dt(erythrocytes2) <- kcir * (erythrocytes1 - erythrocytes2)
    d/dt(erythrocytes3) <- kcir * (erythrocytes2 - erythrocytes3)
    d/dt(erythrocytes4) <- kcir * (erythrocytes3 - erythrocytes4)

    prol(0)           <- prol_bin
    precursor1(0)     <- prol_bin
    precursor2(0)     <- prol_bin
    precursor3(0)     <- prol_bin
    reticulocytes1(0) <- ret_bin
    reticulocytes2(0) <- ret_bin
    reticulocytes3(0) <- ret_bin
    reticulocytes4(0) <- ret_bin
    erythrocytes1(0)  <- rbc_bin
    erythrocytes2(0)  <- rbc_bin
    erythrocytes3(0)  <- rbc_bin
    erythrocytes4(0)  <- rbc_bin

    # --- Hemoglobin (Zhang 2017 Eq. 20) ---------------------------------------
    # "The HGB content was calculated as proportional to the sum of the amounts
    # in the RET count and RBC count compartments, based on SHB."
    hb <- shb * (RET + RBC)

    Cc  ~ lnorm(expSd)
    RET ~ add(addSd_RET)
    RBC ~ add(addSd_RBC)
    hb  ~ add(addSd_hb)
  })
}
