Bian_2024_lefamulin_original_ppb <- function() {
  description <- paste(
    "Three-compartment population PK model for lefamulin with parallel",
    "immediate and delayed oral absorption, saturable plasma protein binding,",
    "and an epithelial lining fluid (ELF) effect compartment, using the",
    "ORIGINAL (diluted) plasma protein binding parameterisation (Bian 2024",
    "Supplementary Table S1). The system is parameterised on UNBOUND",
    "concentrations: all clearances and volumes are free-drug referenced, so",
    "the central state divided by vc is the unbound plasma concentration and",
    "the measured TOTAL plasma concentration is recovered by dividing by the",
    "concentration-dependent unbound fraction fu. An oral dose is split by a",
    "logit-scale total bioavailability into an immediate depot (rate ka) and a",
    "delayed depot2 feeding a three-compartment transit chain (rate ka2),",
    "sharing one lag time; food reduces ka, ka2 and total bioavailability. An",
    "intravenous dose enters central directly with bioavailability 1. ELF is a",
    "concentration-valued effect compartment driven by the unbound plasma",
    "concentration through kin_elf with first-order loss kout_elf, so the",
    "steady-state AUC(ELF) / fAUC(plasma) ratio is kin_elf / kout_elf = 5.3.",
    "Covariates are serum albumin and study phase on CL, study phase on CLd1,",
    "and study phase plus body weight on Vp1. The companion model",
    "Bian_2024_lefamulin_higher_ppb refits the same structure to the",
    "non-diluted (higher) protein binding assumption and replaces the ELF",
    "compartment with the FDA review's algebraic power penetration function."
  )
  reference <- paste(
    "Bian X, Li N, Li Y, Zhu X, Yu J, Hu Y, Yang H, Wei Q, Wu X, Wang J,",
    "Cao G, Wu J, Wang Y, Zhang J. Lefamulin dosing optimization using",
    "population pharmacokinetic and pharmacokinetic/pharmacodynamic assessment",
    "in Chinese patients with community-acquired bacterial pneumonia.",
    "Front Pharmacol. 2024;15:1456741. doi:10.3389/fphar.2024.1456741.",
    "Structural equations from Equations 1-20; parameter values from",
    "Supplementary Table S1. The underlying popPK plus ELF framework was",
    "developed in Zhang L, Wicha SG, Bhavnani SM, et al.",
    "J Antimicrob Chemother. 2019;74(Suppl 3):iii27-iii34.",
    "doi:10.1093/jac/dkz088, and in the FDA XENLETA Multi-Discipline Review,",
    "NDA 211672/211673 (2019).",
    sep = " "
  )
  vignette <- "Bian_2024_lefamulin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Oral doses target depot (immediate) and depot2 (delayed); intravenous
  # doses target central. Declared explicitly because buildModelDb()'s
  # heuristic only recognises literal "depot" / "central".
  dosing <- c("depot", "depot2", "central")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: analyte and specimen checked against
  # Bian 2024 Figure 2 (model diagram) and the Methods text.
  compartmentData <- list(
    depot       = list(analyte = "lefamulin", units = "mg",   specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "lefamulin", units = "mg",   specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "lefamulin", units = "mg",   specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "lefamulin", units = "mg",   specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "lefamulin", units = "mg",   specimen = "administration site", verified = TRUE),
    central     = list(analyte = "lefamulin", units = "mg",   specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lefamulin", units = "mg",   specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "lefamulin", units = "mg",   specimen = "plasma", verified = TRUE),
    # The ELF state is a Sheiner-style effect compartment and therefore holds a
    # CONCENTRATION (mg/L), not an amount: Bian 2024 Methods, "ELF, as an
    # effect compartment, was linked to the central compartment of the PopPK
    # model", with no ELF volume reported anywhere in the paper, Zhang 2019, or
    # the FDA review.
    elf         = list(analyte = "lefamulin", units = "mg/L", specimen = "epithelial lining fluid", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear (non-power) proportional shift on CL:",
        "CL = [1 + e_alb_cl * (ALB_gdL - 4.1)] * clphase * theta1 (Bian 2024",
        "Equation 1). The structural coefficient was calibrated in g/dL, so",
        "model() applies the inline SI-to-US conversion alb_gdL <- ALB * 0.1",
        "required by the covariate register's ALB Units note (same construction",
        "as Fasanmade_2009_infliximab.R and Roepcke_2023_rezafungin.R).",
        "Reference 4.1 g/dL is the pooled model-building cohort median albumin",
        "(Supplementary Table S5, Total column). Supplementary Table S1 prints",
        "the row 'The effect of albumin on CL' as 1.214, which is (1 + theta20)",
        "rather than theta20 itself; see the ini() source-trace comment for the",
        "three independent confirmations. Baseline (time-fixed) per subject."
      ),
      source_name        = "ALB"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear proportional shift on the first peripheral volume:",
        "Vp1 = theta4 * vp1phase * (1 + e_wt_vp * (WT - 78)) (Bian 2024",
        "Equation 4). Reference 78 kg is the pooled model-building cohort",
        "median weight (Supplementary Table S5, Total column; range",
        "31-174.6 kg). Supplementary Table S1 prints the row 'The effect of",
        "WTKG on Vp1' as 1.0129, i.e. (1 + theta27) with theta27 = 0.0129 per",
        "kg. Baseline (time-fixed) per subject."
      ),
      source_name        = "WTKG"
    ),
    FED = list(
      description        = "Fed state at the time of the oral dose",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = fasted (the FAST = 1 stratum of Bian 2024 Equations 7-9)",
      notes              = paste(
        "1 = dose administered fed, 0 = fasted. Bian 2024 Equations 7-9 write",
        "each affected parameter as theta * FAST + theta * theta_fed * FEDD",
        "with FAST = 1 - FED and FEDD = FED, i.e. the tabulated fed values are",
        "MULTIPLIERS applied to the fasted typical value, not absolute fed",
        "estimates -- despite Supplementary Table S1 labelling theta17 and",
        "theta19 with 1/hr units. Food multiplies ka by 0.0541, ka2 by 0.445",
        "and total oral bioavailability by 0.763. Applies to oral doses only;",
        "intravenous doses bypass both depots."
      ),
      source_name        = "FEDD"
    ),
    STUDY_LEFAMULIN_PHASE1 = list(
      description        = "Phase 1 study cohort indicator in the pooled lefamulin popPK analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with STUDY_LEFAMULIN_PHASE2 also 0 (Phase 3 studies; the reference stratum)",
      notes              = paste(
        "1 = subject enrolled in one of the Phase 1 studies (n = 98) of the",
        "pooled model-building dataset, 0 = otherwise. Paired with",
        "STUDY_LEFAMULIN_PHASE2; both 0 selects Phase 3 (n = 622), the",
        "reference stratum with all three phase multipliers equal to 1 (Bian",
        "2024 Equations 1, 3, 4: 'When PHASE = 3, CLPHASE = 1'). Shifts CL,",
        "CLd1 and Vp1. Subject-level (time-fixed)."
      ),
      source_name        = "PHASE"
    ),
    STUDY_LEFAMULIN_PHASE2 = list(
      description        = "Phase 2 study cohort indicator in the pooled lefamulin popPK analysis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 with STUDY_LEFAMULIN_PHASE1 also 0 (Phase 3 studies; the reference stratum)",
      notes              = paste(
        "1 = subject enrolled in the Phase 2 study (n = 129) of the pooled",
        "model-building dataset, 0 = otherwise. Paired with",
        "STUDY_LEFAMULIN_PHASE1. Shifts CL, CLd1 and Vp1.",
        "Subject-level (time-fixed)."
      ),
      source_name        = "PHASE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 849L,
    n_studies      = 3L,
    age_range      = "18-97 years",
    age_median     = "57 years",
    weight_range   = "31-174.6 kg",
    weight_median  = "78 kg",
    height_range   = "133-200 cm",
    height_median  = "170 cm",
    sex_female_pct = 38.8,
    race_ethnicity = c(White = 78.2, Black = 8.13, Asian = 8.48,
                       `American Indian/Alaskan Native` = 3.65,
                       `Native Hawaiian/Other Pacific Islander` = 0.471,
                       Other = 1.06),
    disease_state  = paste(
      "Pooled Phase 1 healthy adults (n = 98), Phase 2 acute bacterial skin",
      "and skin structure infection patients (n = 129), and Phase 3",
      "community-acquired bacterial pneumonia (CABP) patients (n = 622)."
    ),
    dose_range     = paste(
      "150 mg q12h intravenous infusion over 1, 1.5 or 2 h, and 600 mg q12h",
      "oral under fasted or fed conditions (the regimens simulated in",
      "Supplementary Table S6)."
    ),
    regions        = "Foreign (predominantly North American and European) model-building studies; externally validated in China.",
    albumin_range  = "2-5.6 g/dL",
    albumin_median = "4.1 g/dL",
    renal_function = "Creatinine clearance 6.2-192.4 mL/min (median 80.2 mL/min); CLcr was not retained as a covariate.",
    notes          = paste(
      "Baseline demographics from Bian 2024 Supplementary Table S5. The final",
      "model was taken from the previously published foreign analysis (Zhang",
      "2019; FDA XENLETA review) and externally validated against 33 Chinese",
      "CABP patients from a Phase III study (Supplementary Table S5, External",
      "validation column: median age 56.12 years, weight 69.38 kg, albumin",
      "3.88 g/dL, 100% Asian, 18.18% female). Monte Carlo simulations were run",
      "on 5000 virtual CABP patients resampled from 125 Chinese CABP patients."
    )
  )

  ini({
    # ================================================================
    # Structural fixed effects -- Bian 2024 Supplementary Table S1
    # ("Parameter list for the lefamulin final PopPK model with
    # original plasma protein binding"). Every clearance and volume is
    # FREE-DRUG referenced: Zhang 2019 Table 2 footnote states
    # "Volumes and clearances are scaled to free-drug concentrations"
    # and its Results that "all PK parameters (e.g. clearances and
    # volumes) were conditioned on lefamulin unbound concentrations".
    # Self-check: Dose / CL = 300 mg per day / 79.4 L/h = 3.78 mg*h/L,
    # which reproduces the Table S6 steady-state fAUC(0-24) of
    # 3.87 mg*h/L for 150 mg q12h IV.
    # ================================================================
    lcl  <- log(79.4); label("Clearance, free-drug referenced (L/h)")                      # theta1, Table S1
    lvc  <- log(46.3); label("Central volume of distribution, free-drug referenced (L)")   # theta2, Table S1
    lq   <- log(40.6); label("Distributional clearance to peripheral1, free-drug referenced (L/h)") # theta3, Table S1
    lvp  <- log(249);  label("Peripheral1 volume of distribution, free-drug referenced (L)")        # theta4, Table S1
    lq2  <- log(199);  label("Distributional clearance to peripheral2, free-drug referenced (L/h)") # theta5, Table S1
    lvp2 <- log(259);  label("Peripheral2 volume of distribution, free-drug referenced (L)")        # theta6, Table S1

    # Absorption. Two parallel oral routes share one lag time
    # (Equation 15). ka drives the immediate depot, ka2 drives the
    # delayed depot2 -> transit1 -> transit2 -> transit3 chain
    # (Figure 2: Depot2 -Ka2-> Abs1 -Ka2-> ... -Ka2-> Abs3 -Ka2->
    # Plasma; the diagram elides Abs2 behind a dotted line).
    lka   <- log(1.2);  label("First-order absorption rate constant, immediate route (1/h)") # theta7, Table S1
    lka2  <- log(2.12); label("First-order absorption rate constant, delayed transit route (1/h)") # theta8, Table S1
    ltlag <- log(0.15); label("Absorption lag time, both oral routes (h)")                   # theta11 ALAG, Table S1

    # Oral bioavailability. Equation 10 puts the IIV on the LOGIT of
    # the total bioavailability and Equations 11-12 do the same for
    # the slow-route fraction, so both typical values are carried on
    # the logit scale here. frac is the canonical generic bounded
    # fraction (parameter register); it is the fraction of the
    # ABSORBED dose routed through the delayed pathway, so that
    # F1 = fdepot * (1 - frac) and F2 = fdepot * frac
    # (Equations 13, 14).
    logitfdepot <- log(0.244 / (1 - 0.244)); label("Logit of total oral bioavailability Ftot (unitless)") # theta9 = 0.244, Table S1
    logitfrac   <- log(0.802 / (1 - 0.802)); label("Logit of the fraction of absorbed dose taking the delayed route FS (unitless)") # theta10 = 0.802, Table S1

    # ================================================================
    # Saturable plasma protein binding -- Equations 16-19, an Emax
    # (Michaelis-Menten) rise of the unbound fraction with the UNBOUND
    # plasma concentration:
    #   fu = fumin + fumax * Cu / (cup50 + Cu)
    # All three were held FIXED (Table S1 marks each "Fixed"). The FDA
    # XENLETA Multi-Discipline Review (NDA 211672/211673, p. 245)
    # states the submodel and the driver explicitly: "fu, min =
    # population minimum unbound fraction fixed at 0.0997; fu, max =
    # population maximum unbound fraction fixed at 0.259; Cu, plasma =
    # unbound plasma concentration; Cu, plasma 50 = population
    # Cu,plasma where fu is increased by half".
    # ================================================================
    fumin <- fixed(0.0997); label("Minimum (zero-concentration) unbound fraction in plasma (unitless)") # theta12, Table S1, Fixed
    fumax <- fixed(0.259);  label("Maximum increment in unbound fraction at saturating concentration (unitless)") # theta13, Table S1, Fixed
    cup50 <- fixed(1.35);   label("Unbound plasma concentration giving half the maximal increment in fu (mg/L)")  # theta14 fu50, Table S1, Fixed

    # ================================================================
    # Food effect. Equations 7-9 have the shape
    #   parameter = theta * FAST + theta * theta_fed * FEDD,
    # so each tabulated "fed" theta is a MULTIPLIER on the fasted
    # typical value, not an absolute fed estimate. Table S1's unit
    # annotations "(1/hr)" on theta17 and theta19 are inconsistent with
    # the printed equations; per the standing rule the equation wins.
    # ================================================================
    e_fed_ka     <- 0.0541; label("Fed-state multiplier on ka (unitless)")                     # theta19 Kafed, Table S1
    e_fed_ka2    <- 0.445;  label("Fed-state multiplier on ka2 (unitless)")                    # theta17 Ka2fed, Table S1
    e_fed_fdepot <- 0.763;  label("Fed-state multiplier on total oral bioavailability (unitless)") # theta18 Ftot,fed, Table S1

    # ================================================================
    # Covariate effects. Table S1 prints these rows as the MULTIPLIER
    # (1 + theta), not as theta, so each value below is the tabulated
    # number minus 1. Three independent confirmations:
    #  (a) theta27 is printed as 1.0129 and theta27 of the higher-ppb
    #      companion table as 1.00637 -- four- and five-significant-
    #      figure values that are only interpretable as 1 + 0.0129 and
    #      1 + 0.00637;
    #  (b) reading them as theta directly makes CL negative at the
    #      observed albumin minimum of 2.0 g/dL
    #      (1 + 1.214 * (2 - 4.1) = -1.55) and Vp1 negative at 31 kg;
    #  (c) the upstream Zhang 2019 Phase-1-only estimates are
    #      reproduced by the multiplier reading: 40.6 * 2.12 = 86.1 vs
    #      Zhang's CLd1 = 86.6 L/h, and 79.4 * 1.766 * [1 + 0.214 *
    #      (4.5 - 4.1)] = 152 vs Zhang's CL = 159 L/h. The paper's own
    #      prose ("AUC0-24 and Cmax in CABP patients was approximately
    #      1.73- and 1.3-fold greater compared to adults without
    #      pneumonia") matches the 1.766 Phase 1 CL multiplier.
    # Phase 3 is the reference stratum (all multipliers 1).
    # ================================================================
    e_alb_cl    <- 0.214; label("Additive effect of albumin on CL, per g/dL above 4.1 g/dL (1/(g/dL))") # theta20; Table S1 prints 1 + theta20 = 1.214
    e_phase1_cl <- 0.766; label("Proportional shift in CL for Phase 1 subjects (unitless)")   # theta22; Table S1 prints 1 + theta22 = 1.766
    e_phase2_cl <- 0.827; label("Proportional shift in CL for Phase 2 subjects (unitless)")   # theta21; Table S1 prints 1 + theta21 = 1.827
    e_phase1_q  <- 1.12;  label("Proportional shift in CLd1 for Phase 1 subjects (unitless)") # theta24; Table S1 prints 1 + theta24 = 2.12
    e_phase2_q  <- 0.44;  label("Proportional shift in CLd1 for Phase 2 subjects (unitless)") # theta23; Table S1 prints 1 + theta23 = 1.44
    e_phase1_vp <- 1.75;  label("Proportional shift in Vp1 for Phase 1 subjects (unitless)")  # theta26; Table S1 prints 1 + theta26 = 2.75
    e_phase2_vp <- 0.985; label("Proportional shift in Vp1 for Phase 2 subjects (unitless)")  # theta25; Table S1 prints 1 + theta25 = 1.985
    e_wt_vp     <- 0.0129; label("Additive effect of body weight on Vp1, per kg above 78 kg (1/kg)") # theta27; Table S1 prints 1 + theta27 = 1.0129

    # ================================================================
    # Epithelial lining fluid effect compartment. Bian 2024 Methods:
    # "the distribution rate constant from the central compartment to
    # the ELF compartment (KELF,in) and the elimination rate constant
    # from the ELF compartment (KEFL,out) ... are derived from the
    # previously established model" (Figure 2; Supplementary Table S1).
    # The driver is the UNBOUND plasma concentration, which is what
    # makes the steady-state exposure ratio exactly
    # kin_elf / kout_elf = 2.71 / 0.51 = 5.31 -- corroborated by Zhang
    # 2019 Results ("the median lefamulin total-drug ELF AUC(0-24):
    # free-drug plasma AUC(0-24) ratio was ~5:1") and by Table S6,
    # where 20.51 / 3.87 = 5.30 for the 150 mg q12h IV Day 3 arm.
    # ================================================================
    lkin_elf  <- log(2.71); label("Plasma-to-ELF distribution rate constant (1/h)") # theta28 KELF,in, Table S1
    lkout_elf <- log(0.51); label("ELF elimination rate constant (1/h)")            # theta29 KELF,out, Table S1

    # ================================================================
    # Inter-individual variability -- Table S1 "Interindividual
    # variability" block, reported as omega^2 with the sqrt shown as a
    # %CV (e.g. 0.171 -> 41.4%CV), so the tabulated numbers are
    # VARIANCES and are used directly. All etas are exponential except
    # eta9 and eta10, which are additive on the logit scale
    # (Equations 10, 12).
    #
    # Equations 5 and 6 additionally carry exp(eta5) on CLd2 and
    # exp(eta5) * exp(eta6) on Vp2, but Table S1 reports no eta5 or
    # eta6 variance, so both were not estimated (variance 0) and are
    # omitted here rather than encoded as a singular fixed(0) omega
    # diagonal.
    # ================================================================
    etalcl         ~ 0.171; label("IIV on CL (variance)")                   # eta1, Table S1 (41.4%CV)
    etalvc         ~ 0.39;  label("IIV on Vc (variance)")                   # eta2, Table S1 (62.4%CV)
    etalq          ~ 0.119; label("IIV on CLd1 (variance)")                 # eta3, Table S1 (34.5%CV)
    etalvp         ~ 0.623; label("IIV on Vp1 (variance)")                  # eta4, Table S1 (78.9%CV)
    etalka         ~ 0.800; label("IIV on ka (variance)")                   # eta7, Table S1 (89.4%CV)
    etalka2        ~ 0.400; label("IIV on ka2 (variance)")                  # eta8, Table S1 (63.2%CV)
    etalogitfdepot ~ 0.100; label("IIV on logit total oral bioavailability (variance)") # eta9, Table S1 (31.6%CV)
    etalogitfrac   ~ 0.170; label("IIV on logit delayed-route fraction (variance)")     # eta10, Table S1 (41.2%CV)

    # ================================================================
    # Residual unexplained variability -- Equation 20:
    #   Y = IPRED + IPRED*eps1*PLASMA + eps2*PLASMA + IPRED*eps3*ELF
    # i.e. combined proportional + additive on plasma and proportional
    # only on ELF. Table S1 reports sigma^2; the SDs below are the
    # square roots, which reproduce the table's own parenthetical
    # conversions (sqrt(0.103) = 0.321 -> 32.0%CV;
    # sqrt(0.0000343) = 0.00586 mg/L; sqrt(0.372) = 0.610 -> 61.0%CV).
    # ================================================================
    propSd      <- 0.3209;  label("Proportional residual error, plasma (fraction)") # sqrt(sigma^2 eps1 = 0.103), Table S1
    addSd       <- 0.005857; label("Additive residual error, plasma (mg/L)")        # sqrt(sigma^2 eps2 = 0.0000343), Table S1
    propSd_Celf <- 0.6099;  label("Proportional residual error, ELF (fraction)")    # sqrt(sigma^2 eps3 = 0.372), Table S1
  })

  model({
    # ---- Covariate preparation ------------------------------------
    # The register's canonical ALB unit is g/L (SI); Equation 1 was
    # calibrated against g/dL, so convert inline.
    alb_gdL <- ALB * 0.1

    # Study-phase multipliers. Phase 3 is the reference stratum, so
    # both indicators 0 gives a multiplier of 1 (Equations 1, 3, 4).
    clphase <- 1 + e_phase1_cl * STUDY_LEFAMULIN_PHASE1 +
                   e_phase2_cl * STUDY_LEFAMULIN_PHASE2
    qphase  <- 1 + e_phase1_q  * STUDY_LEFAMULIN_PHASE1 +
                   e_phase2_q  * STUDY_LEFAMULIN_PHASE2
    vpphase <- 1 + e_phase1_vp * STUDY_LEFAMULIN_PHASE1 +
                   e_phase2_vp * STUDY_LEFAMULIN_PHASE2

    # ---- Individual PK parameters ---------------------------------
    # Equations 1-6. Equations 3 and 4 as printed omit theta3 and
    # theta4 from the right-hand side; Supplementary Table S1 supplies
    # them as CLd1 = 40.6 L/hr and Vp1 = 249 L, and the units of the
    # printed left-hand sides confirm the omission is a typesetting
    # slip (a bare multiplier cannot have units of L/hr or L).
    cl  <- exp(lcl + etalcl) * (1 + e_alb_cl * (alb_gdL - 4.1)) * clphase
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + etalq) * qphase
    vp  <- exp(lvp + etalvp) * vpphase * (1 + e_wt_vp * (WT - 78))
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    # Absorption rate constants with the multiplicative food effect
    # (Equations 7-8; FAST = 1 - FED, FEDD = FED).
    ka   <- exp(lka + etalka)   * (1 - FED + e_fed_ka  * FED)
    ka2  <- exp(lka2 + etalka2) * (1 - FED + e_fed_ka2 * FED)
    tlag <- exp(ltlag)

    # Oral bioavailability. Equation 9 applies the food effect on the
    # LINEAR scale to theta9; Equation 10 then logit-transforms the
    # result and adds eta9 on the logit scale. Equations 11-12 do the
    # same for the delayed-route fraction, whose typical value is
    # already carried on the logit scale in ini().
    fdepot_typ <- exp(logitfdepot) / (1 + exp(logitfdepot))
    bio        <- fdepot_typ * (1 - FED + e_fed_fdepot * FED)
    fpo        <- log(bio / (1 - bio))
    fdepot     <- exp(fpo + etalogitfdepot) / (1 + exp(fpo + etalogitfdepot))
    frac       <- exp(logitfrac + etalogitfrac) /
                  (1 + exp(logitfrac + etalogitfrac))

    # ELF micro-constants.
    kin_elf  <- exp(lkin_elf)
    kout_elf <- exp(lkout_elf)

    # ---- Concentrations -------------------------------------------
    # vc is free-drug referenced, so central / vc is the UNBOUND
    # plasma concentration. Equation 16 then gives the unbound
    # fraction and the measured TOTAL plasma concentration is Cu / fu.
    Cu <- central / vc
    fu <- fumin + fumax * Cu / (cup50 + Cu)

    # ---- ODE system (Figure 2) ------------------------------------
    # Immediate oral route: depot -ka-> central.
    # Delayed oral route:   depot2 -ka2-> transit1 -ka2-> transit2
    #                       -ka2-> transit3 -ka2-> central.
    # Intravenous doses are administered directly to central.
    # All disposition fluxes are driven by the unbound concentration
    # because the clearances and volumes are free-drug referenced.
    d/dt(depot)       <- -ka * depot
    d/dt(depot2)      <- -ka2 * depot2
    d/dt(transit1)    <-  ka2 * depot2   - ka2 * transit1
    d/dt(transit2)    <-  ka2 * transit1 - ka2 * transit2
    d/dt(transit3)    <-  ka2 * transit2 - ka2 * transit3
    d/dt(central)     <-  ka * depot + ka2 * transit3 -
                          cl * Cu -
                          q  * (Cu - peripheral1 / vp) -
                          q2 * (Cu - peripheral2 / vp2)
    d/dt(peripheral1) <-  q  * (Cu - peripheral1 / vp)
    d/dt(peripheral2) <-  q2 * (Cu - peripheral2 / vp2)
    # ELF is an effect compartment: it is driven by the unbound plasma
    # concentration and takes no mass out of central, which is what
    # makes AUC(ELF) / fAUC(plasma) equal kin_elf / kout_elf exactly.
    d/dt(elf)         <-  kin_elf * Cu - kout_elf * elf

    # ---- Bioavailability and lag time (Equations 13-15) -----------
    f(depot)     <- fdepot * (1 - frac)
    f(depot2)    <- fdepot * frac
    alag(depot)  <- tlag
    alag(depot2) <- tlag

    # ---- Observations and residual error (Equation 20) ------------
    Cc   <- Cu / fu
    Celf <- elf
    Cc   ~ add(addSd) + prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
