DeJongh_2025_azd7648_olaparib_xenograft_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, female SCID with FaDu ATM-knockout xenograft).",
    "Minimal systems PK-PD model for tumour growth inhibition by the",
    "DNA-PK inhibitor AZD7648 combined with the PARP inhibitor olaparib.",
    "The two oral PK models (AZD7648 with parallel linear and",
    "Michaelis-Menten elimination; olaparib with a dose-driven",
    "drug-drug-interaction term on clearance) are carried over as fixed",
    "typical values and feed separate biophase compartments. The tumour",
    "is a four-state cell chain: proliferating cells grow logistically",
    "(Verhulst) and are pushed into a reversible DNA-damage-induced",
    "quiescent state by olaparib, from which AZD7648 blocks non-homologous",
    "end-joining repair back to proliferation; cells that do not recover",
    "pass through two dying-cell transit states before removal. An",
    "Allee-effect term multiplies proliferation, so tumours whose",
    "proliferating mass falls far enough regress permanently. Amounts are",
    "molar (umol/kg), plasma concentrations uM, and the observed tumour",
    "volume is in cm3. Parameter values from DeJongh 2025 Tables 1-3."
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
    concentration = "uM (plasma); the tumour output tumor_size is a volume in cm3"
  )

  compartmentData <- list(
    depot                 = list(analyte = "AZD7648",  units = "umol/kg", specimen = "administration site",  verified = TRUE),
    central               = list(analyte = "AZD7648",  units = "umol/kg", specimen = "plasma",                verified = TRUE),
    peripheral1           = list(analyte = "AZD7648",  units = "umol/kg", specimen = "tissue",     verified = TRUE),
    effect                = list(analyte = "AZD7648",  units = "umol/kg", specimen = "tumor",       verified = TRUE),
    depot_olaparib        = list(analyte = "olaparib", units = "umol/kg", specimen = "administration site",  verified = TRUE),
    central_olaparib      = list(analyte = "olaparib", units = "umol/kg", specimen = "plasma",                verified = TRUE),
    peripheral1_olaparib  = list(analyte = "olaparib", units = "umol/kg", specimen = "tissue",     verified = TRUE),
    effect_olaparib       = list(analyte = "olaparib", units = "umol/kg", specimen = "tumor",       verified = TRUE),
    cycling_cells         = list(analyte = "cells",    units = "cm3",     specimen = "tumor",             verified = TRUE),
    damaged_cells1        = list(analyte = "cells",    units = "cm3",     specimen = "tumor",             verified = TRUE),
    damaged_cells2        = list(analyte = "cells",    units = "cm3",     specimen = "tumor",             verified = TRUE),
    damaged_cells3        = list(analyte = "cells",    units = "cm3",     specimen = "tumor",             verified = TRUE)
  )

  covariateData <- list(
    DOSE_AZD7648_MGKGD = list(
      source_name        = "AZD7648",
      description        = paste(
        "Total daily dose of AZD7648 in mg/kg/day. Drives the Emax-shaped",
        "reduction of olaparib clearance; 0 for olaparib monotherapy and",
        "for vehicle controls. This is the total across the day, so a",
        "75 mg/kg twice-daily regimen carries the value 150."
      ),
      units              = "mg/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "DeJongh 2025 Equations 7-8. Carried as a covariate column rather",
        "than derived from the AZD7648 event records because the source",
        "NONMEM dataset supplies it that way (DOSE_AZD column of",
        "supplementary file 7), and because the interaction is a function",
        "of the daily dose rather than of the instantaneous amount."
      )
    )
  )

  population <- list(
    species        = "mouse (female SCID C.B-17/IcrHan(R)Hsd-Prkdcscid bearing a subcutaneous FaDu ATM-knockout xenograft)",
    n_subjects     = 120L,
    n_studies      = 1L,
    age_range      = "at least 6 weeks of age at implantation",
    weight_range   = "(not reported in the source publication)",
    sex_female_pct = 100,
    disease_state  = "subcutaneous FaDu ATM-knockout head-and-neck squamous-carcinoma xenograft in the dorsal left flank; randomised to treatment at a mean tumour volume of about 0.2 cm3",
    dose_range     = paste(
      "Study S1822 (model calibration): eight groups of 15 mice dosed",
      "orally once daily for 28 days with vehicle, AZD7648 50-100 mg/kg,",
      "olaparib 50-100 mg/kg, or five AZD7648 + olaparib combinations,",
      "followed by washout of up to eight weeks. Olaparib was given 1 h",
      "after the morning AZD7648 dose. Doses enter the model in umol/kg",
      "(AZD7648 100 mg/kg = 262.8812 umol/kg; olaparib 100 mg/kg =",
      "229.8428 umol/kg)."
    ),
    regions        = "preclinical (AstraZeneca, Cambridge UK; modelling by LAP&P Consultants, Leiden)",
    notes          = paste(
      "Parameters in Table 3 were fitted simultaneously to all eight",
      "treatment groups of study S1822 (120 mice). Tumour volume was",
      "measured by bilateral Vernier calliper (length x width) and",
      "computed with Mousetrap software. Only the fixed effects of the",
      "two population PK models were carried into the PK-PD analysis,",
      "because most PK data came from satellite groups; the PK-PD model",
      "therefore has no PK random effects. Estimation was FOCEI in",
      "NONMEM 7.4.3 with ADVAN13 (TOL=6). Condition number 40.12.",
      "Eta shrinkage was 30-40% for the growth rate and the Allee",
      "coefficient. Studies S1842 (six-week intermittent schedules) and",
      "S1826 were held out for a priori forecast evaluation and did not",
      "contribute to these estimates. See DeJongh 2025 Table 3 and",
      "supplementary Table S1."
    )
  )

  ini({
    # ==================================================================
    # PHARMACOKINETICS -- carried over as FIXED typical values.
    # DeJongh 2025 Results: "For both compounds, the estimated values for
    # PK model parameters were carried over to the PK-PD analysis. [...]
    # As most PK data was obtained from satellite groups, only fixed
    # effects from the population PK models for olaparib and AZD7648 were
    # taken forward." The PK-PD control stream (supplementary file 4)
    # accordingly hard-codes every PK parameter with its random effects
    # commented out.
    # ==================================================================

    # --- AZD7648 (DeJongh 2025 Table 1) -------------------------------
    # The TGI studies used SCID mice only, so the SCID absorption rate
    # applies and relative bioavailability is 1 (the control stream keeps
    # the nude-mouse branch "for later ref" but it is never active here).
    lka   <- fixed(log(2.7726));  label("AZD7648 absorption rate constant (ka, 1/h)")                     # Table 1 Ka* = 2.77; PK-PD control stream THETA(1) = 2.7726
    lcl   <- fixed(log(0.25119)); label("AZD7648 linear clearance (CL, L/kg/h)")                          # Table 1 CL = 0.251
    lvmax <- fixed(log(4.2741));  label("AZD7648 maximum saturable elimination rate (Vmax, umol/kg/h)")   # Table 1 Vmax = 4.27
    lkm   <- fixed(log(3.7268));  label("AZD7648 concentration at half-maximal elimination (Km, uM)")     # Table 1 Km = 3.73
    lvc   <- fixed(log(3.4475));  label("AZD7648 central volume of distribution (Vc, L/kg)")              # Table 1 Vc = 3.45; Vp is set equal to Vc per Table 1
    lq    <- fixed(log(0.93210)); label("AZD7648 intercompartmental clearance (Q, L/kg/h)")               # Table 1 Q = 0.932

    # --- Olaparib (DeJongh 2025 Table 2) ------------------------------
    lka_olaparib     <- fixed(log(10.0));   label("Olaparib absorption rate constant (ka, 1/h)")                  # Table 2: 10.0, fixed by assumption for lack of absorption-phase data
    lcl_olaparib     <- fixed(log(1.9460)); label("Olaparib clearance (CL, L/h)")                                 # Table 2 CL = 1.95
    lvc_olaparib     <- fixed(log(1.1656)); label("Olaparib central volume of distribution (Vc, L/kg)")           # Table 2 Vc = 1.17
    lq_olaparib      <- fixed(log(0.6676)); label("Olaparib intercompartmental clearance (Q, L/h)")               # Table 2 Q = 0.668
    lvp_olaparib     <- fixed(log(2.6608)); label("Olaparib peripheral volume of distribution (Vp, L/kg)")        # Table 2 Vp = 2.66
    lfdepot_olaparib <- fixed(log(4.1817)); label("Olaparib relative oral bioavailability in the xenograft studies (F1, fraction)")  # Table 2 "F1 1816 (SCID)" = 4.18; PK-PD control stream F11 = 4.1817 ("rel. F1, based on S_11448/S_1816")
    lcmd50           <- fixed(log(82.538)); label("AZD7648 daily dose halving olaparib clearance (CmD50, mg/kg/day)")                # Table 2: 82.5

    # ==================================================================
    # TUMOUR PK-PD -- DeJongh 2025 Table 3 (estimated).
    # ==================================================================
    lrbase      <- log(0.115);  label("Xenograft volume at baseline (BL, cm3)")                                   # Table 3 BL = 0.115 (RSE 3.9%)
    lkge        <- log(0.329);  label("Tumour growth rate coefficient (GC, 1/day)")                               # Table 3 GC = 0.329 (RSE 5.6%)
    ltsmax      <- log(2.25);   label("Carrying capacity (KC, cm3)")                                              # Table 3 KC = 2.25 (RSE 5.8%)
    lec50_allee <- log(1.93);   label("Proliferating-cell volume giving 50% Allee effect (AEC, mm3)")             # Table 3 AEC = 1.93 mm3 (RSE 11%); converted to cm3 in model()
    lhill       <- log(0.637);  label("Power coefficient on the Allee function (G, unitless)")                    # Table 3 G = 0.637 (RSE 9.4%)
    lke0        <- log(0.961);  label("Biophase equilibration rate constant (Ktr, 1/day)")                        # Table 3 Ktr = 0.961 (RSE 17%); control stream Ke0 = THETA(6)/24
    lec50       <- log(8.99);   label("AZD7648 concentration giving 50% inhibition of quiescent-to-proliferating transition (EC50az, uM)")  # Table 3 EC50az = 8.99 (RSE 7.4%)

    lslope_olaparib <- log(0.0154); label("Olaparib linear effect coefficient on proliferating-to-quiescent transition (EColap, 1/uM)")     # Table 3 EColap = 0.0154 (RSE 4.5%)

    # Kto, the turnover rate constant shared by every cell-state
    # transition, is NOT a separate parameter: DeJongh 2025 Table 3
    # footnote c states "Kturn-over was assumed to be equal to the Growth
    # rate Coefficient GC, as these parameters could not be identified
    # independently from the model fit to the data". The control stream
    # implements this as `Kto = GC`. It is set in model() below.

    # ------------------------------------------------------------------
    # Inter-individual variability -- DeJongh 2025 Table 3.
    # A 2x2 block on baseline volume and the Allee coefficient, plus a
    # diagonal element on the growth coefficient.
    # ------------------------------------------------------------------
    etalrbase + etalec50_allee ~ c(0.0951,
                                   0.135, 1.83)                                                                   # Table 3: Om_1 (BL) = 0.0951, Om_1 x Om_2 = 0.135, Om_2 (AEC) = 1.83
    etalkge ~ 0.0362                                                                                              # Table 3: Om_3 (GC) = 0.0362

    # ------------------------------------------------------------------
    # Residual error. Unlike Tables 1 and 2, which report the residual
    # error as a variance, Table 3's 0.243 is a standard deviation: the
    # PK-PD control stream uses `W = IPRED*THETA(9)` (no square root) with
    # `$SIGMA 1 FIX`. It is therefore used directly, as 24.3% CV.
    # ------------------------------------------------------------------
    propSd_tumor_size <- 0.243; label("Proportional residual error on tumour volume (fraction)")                   # Table 3 proportional residual error = 0.243 (RSE 2.0%)
  })

  model({
    # --- Unit conversions (source reports per-day rates and mm3) --------
    hoursPerDay  <- 24      # the model runs on the hour timescale of the PK
    cm3PerMm3    <- 0.001   # control stream: AEC = 0.001*THETA(4)

    # --- AZD7648 pharmacokinetics (DeJongh 2025 Eq. 1-3) ---------------
    # Michaelis-Menten term transcribed from the NONMEM $DES block,
    # `Vmax*A(2)/(A(2)/Vc + Km)`; printed Equation 1 omits the amount in
    # the numerator (see vignette Errata).
    ka   <- exp(lka)
    cl   <- exp(lcl)
    vmax <- exp(lvmax)
    km   <- exp(lkm)
    vc   <- exp(lvc)
    vp   <- vc
    q    <- exp(lq)

    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <- ka * depot -
                         (q / vc) * central + (q / vp) * peripheral1 -
                         (cl / vc) * central - vmax * central / (Cc + km)
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1

    # AZD7648 bioavailability is 1: the xenograft studies used SCID mice,
    # the reference strain of the Table 1 relative-F covariate.

    # --- Olaparib pharmacokinetics (DeJongh 2025 Eq. 4-8) --------------
    ka_olaparib     <- exp(lka_olaparib)
    cl_olaparib     <- exp(lcl_olaparib)
    vc_olaparib     <- exp(lvc_olaparib)
    q_olaparib      <- exp(lq_olaparib)
    vp_olaparib     <- exp(lvp_olaparib)
    fdepot_olaparib <- exp(lfdepot_olaparib)

    # Emax-shaped reduction of olaparib clearance by the AZD7648 daily
    # dose. Equations 7 and 8 agree at DoseAZ = 0, so one branch suffices.
    cms <- 1 - DOSE_AZD7648_MGKGD / (DOSE_AZD7648_MGKGD + exp(lcmd50))

    Cc_olaparib <- central_olaparib / vc_olaparib

    d/dt(depot_olaparib)       <- -ka_olaparib * depot_olaparib
    d/dt(central_olaparib)     <- ka_olaparib * depot_olaparib -
                                  (cl_olaparib / vc_olaparib) * cms * central_olaparib +
                                  (q_olaparib / vp_olaparib) * peripheral1_olaparib -
                                  (q_olaparib / vc_olaparib) * central_olaparib
    d/dt(peripheral1_olaparib) <- (q_olaparib / vc_olaparib) * central_olaparib -
                                  (q_olaparib / vp_olaparib) * peripheral1_olaparib

    f(depot_olaparib) <- fdepot_olaparib

    # --- Biophase (effect-site) compartments (Eq. 14-15) ---------------
    # One shared equilibration rate constant for both drugs. The transfer
    # is written on amounts, exactly as in the control stream, and the
    # biophase concentration then divides by the same central volume.
    ke0 <- exp(lke0) / hoursPerDay

    d/dt(effect)          <- ke0 * (central - effect)
    d/dt(effect_olaparib) <- ke0 * (central_olaparib - effect_olaparib)

    Ce          <- effect / vc
    Ce_olaparib <- effect_olaparib / vc_olaparib

    # --- Drug effects (DeJongh 2025 Eq. 16-17) -------------------------
    # azEff multiplies the quiescent-to-proliferating (repair) flux and
    # falls from 1 to 0 as AZD7648 rises, i.e. Emax is fixed at -1
    # (complete inhibition of DNA repair).
    #
    # olEff multiplies the proliferating-to-quiescent flux. The control
    # stream reads `OLEFF = 1 + ECol*COL`: the leading 1 is the untreated
    # baseline turnover between the two states, without which the system
    # would have no cell cycling at all in vehicle animals. Printed
    # Equation 16 gives only the drug term; the control stream governs
    # (see vignette Errata).
    azEff <- 1 - Ce / (Ce + exp(lec50))
    olEff <- 1 + exp(lslope_olaparib) * Ce_olaparib

    # --- Tumour system parameters --------------------------------------
    rbase <- exp(lrbase + etalrbase)
    kge   <- exp(lkge + etalkge) / hoursPerDay
    tsmax <- exp(ltsmax)
    aec   <- cm3PerMm3 * exp(lec50_allee + etalec50_allee)
    hill  <- exp(lhill)

    # Table 3 footnote c: turnover equals the growth rate coefficient.
    kto <- kge

    # Total xenograft volume (Eq. 9), also the model observation.
    tumor_size <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    # Allee term (Eq. 18): proliferation is progressively switched off as
    # the proliferating mass falls below aec.
    alleeEff <- cycling_cells^hill / (cycling_cells^hill + aec^hill)

    # --- Cell-state system (DeJongh 2025 Eq. 10-13) --------------------
    #   cycling_cells  = Prol    (proliferating; paper A9)
    #   damaged_cells1 = Quiesc  (reversible DNA-damage-induced quiescence; A8)
    #   damaged_cells2 = ApopI   (first dying-cell transit; A7)
    #   damaged_cells3 = ApopII  (second dying-cell transit; A6)
    #
    # The logistic limitation term is transcribed from the control stream
    # as `1 - Atot/(2*KC)`; with Kto = GC this makes the untreated steady
    # state exactly KC (supplementary file 2, Remark 1).
    d/dt(cycling_cells)  <- cycling_cells * kge * (1 - tumor_size / (2 * tsmax)) * alleeEff -
                            olEff * kto * cycling_cells +
                            azEff * kto * damaged_cells1
    d/dt(damaged_cells1) <- olEff * kto * cycling_cells -
                            azEff * kto * damaged_cells1 -
                            kto * damaged_cells1
    d/dt(damaged_cells2) <- kto * damaged_cells1 - kto * damaged_cells2
    d/dt(damaged_cells3) <- kto * damaged_cells2 - kto * damaged_cells3

    # Initial conditions: all cells proliferating at study start
    # (control stream BL8 = BL7 = BL6 = 0).
    cycling_cells(0)  <- rbase
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0

    tumor_size ~ prop(propSd_tumor_size)
  })
}
