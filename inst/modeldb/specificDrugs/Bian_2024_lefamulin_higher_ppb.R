Bian_2024_lefamulin_higher_ppb <- function() {
  description <- paste(
    "Three-compartment population PK model for lefamulin with parallel",
    "immediate and delayed oral absorption, saturable plasma protein binding,",
    "and an algebraic epithelial lining fluid (ELF) power penetration",
    "function, using the HIGHER (non-diluted) plasma protein binding",
    "parameterisation (Bian 2024 Supplementary Table S2). The structure is",
    "identical to the companion Bian_2024_lefamulin_original_ppb apart from the",
    "ELF layer: the whole system is refitted under the non-diluted binding",
    "assumption, and ELF concentrations are obtained not from an effect",
    "compartment but from the FDA XENLETA review's lung-penetration-ratio power",
    "model, C_ELF = penratio_elf * (Cc * 0.0379)^pwr_penratio_elf, in which the",
    "total plasma concentration is first scaled to free drug by the review's",
    "fixed linear unbound fraction of 0.0379. All clearances and volumes remain",
    "free-drug referenced, so the central state divided by vc is the unbound",
    "plasma concentration and the measured TOTAL plasma concentration is",
    "recovered by dividing by the concentration-dependent unbound fraction fu."
  )
  reference <- paste(
    "Bian X, Li N, Li Y, Zhu X, Yu J, Hu Y, Yang H, Wei Q, Wu X, Wang J,",
    "Cao G, Wu J, Wang Y, Zhang J. Lefamulin dosing optimization using",
    "population pharmacokinetic and pharmacokinetic/pharmacodynamic assessment",
    "in Chinese patients with community-acquired bacterial pneumonia.",
    "Front Pharmacol. 2024;15:1456741. doi:10.3389/fphar.2024.1456741.",
    "Structural equations from Equations 1-21; parameter values from",
    "Supplementary Table S2. The ELF power penetration function and its",
    "0.0379 fixed unbound fraction are from the FDA XENLETA Multi-Discipline",
    "Review, NDA 211672/211673 (2019), section 16.3.2.4.1 'ELF PK Model'.",
    "The underlying popPK framework was developed in Zhang L, Wicha SG,",
    "Bhavnani SM, et al. J Antimicrob Chemother. 2019;74(Suppl 3):iii27-iii34.",
    "doi:10.1093/jac/dkz088.",
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
  # Bian 2024 Figure 2 (model diagram) and the Methods text. There is no ELF
  # state in this variant -- ELF is algebraic (Equation 21).
  compartmentData <- list(
    depot       = list(analyte = "lefamulin", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "lefamulin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "lefamulin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "lefamulin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "lefamulin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "lefamulin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lefamulin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "lefamulin", units = "mg", specimen = "plasma", verified = TRUE)
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
        "required by the covariate register's ALB Units note. Reference",
        "4.1 g/dL is the pooled model-building cohort median albumin",
        "(Supplementary Table S5, Total column). Supplementary Table S2 prints",
        "the row 'The effect of albumin on CL' as 1.2, which is (1 + theta20)",
        "rather than theta20 itself. Baseline (time-fixed) per subject."
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
        "median weight (Supplementary Table S5, Total column). Supplementary",
        "Table S2 prints the row 'The effect of WTKG on Vp1' as 1.00637, i.e.",
        "(1 + theta27) with theta27 = 0.00637 per kg -- a five-significant-",
        "figure value that is only interpretable as a multiplier. Baseline",
        "(time-fixed) per subject."
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
        "estimates. Under the higher-ppb fit food multiplies ka by 0.0525, ka2",
        "by 0.513 and total oral bioavailability by 0.802. Applies to oral",
        "doses only; intravenous doses bypass both depots."
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
        "2024 Equations 1, 3, 4). Shifts CL, CLd1 and Vp1. Subject-level",
        "(time-fixed)."
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
      "Same model-building population as the companion original-ppb model",
      "(Bian 2024 Supplementary Table S5); the difference is the plasma",
      "protein binding assumption used when fitting. 'Higher ppb' refers to",
      "the non-diluted plasma protein binding rate, which the FDA XENLETA",
      "review adopted after finding that the diluted assay understated",
      "binding; it reduces free-drug exposure roughly 3.5-fold relative to the",
      "original-ppb fit (Supplementary Table S6)."
    )
  )

  ini({
    # ================================================================
    # Structural fixed effects -- Bian 2024 Supplementary Table S2
    # ("Parameter list for the lefamulin final PopPK model with higher
    # plasma protein binding"). As in the original-ppb companion, every
    # clearance and volume is FREE-DRUG referenced (Zhang 2019 Table 2
    # footnote, "Volumes and clearances are scaled to free-drug
    # concentrations").
    # Self-check: Dose / CL = 300 mg per day / 282 L/h = 1.06 mg*h/L,
    # which reproduces the Table S6 higher-ppb steady-state fAUC(0-24)
    # of 1.08 mg*h/L for 150 mg q12h IV.
    # ================================================================
    lcl  <- log(282);  label("Clearance, free-drug referenced (L/h)")                      # theta1, Table S2
    lvc  <- log(138);  label("Central volume of distribution, free-drug referenced (L)")   # theta2, Table S2
    lq   <- log(187);  label("Distributional clearance to peripheral1, free-drug referenced (L/h)") # theta3, Table S2
    lvp  <- log(1300); label("Peripheral1 volume of distribution, free-drug referenced (L)")        # theta4, Table S2
    lq2  <- log(421);  label("Distributional clearance to peripheral2, free-drug referenced (L/h)") # theta5, Table S2
    lvp2 <- log(449);  label("Peripheral2 volume of distribution, free-drug referenced (L)")        # theta6, Table S2

    # Absorption -- see the original-ppb companion for the structure.
    lka   <- log(0.686); label("First-order absorption rate constant, immediate route (1/h)")        # theta7, Table S2
    lka2  <- log(1.34);  label("First-order absorption rate constant, delayed transit route (1/h)")  # theta8, Table S2
    ltlag <- log(0.124); label("Absorption lag time, both oral routes (h)")                          # theta11 ALAG, Table S2

    # Oral bioavailability, carried on the logit scale per Equations
    # 10-12. frac is the fraction of the ABSORBED dose routed through
    # the delayed pathway: F1 = fdepot * (1 - frac),
    # F2 = fdepot * frac (Equations 13, 14).
    logitfdepot <- log(0.293 / (1 - 0.293)); label("Logit of total oral bioavailability Ftot (unitless)") # theta9 = 0.293, Table S2
    logitfrac   <- log(0.622 / (1 - 0.622)); label("Logit of the fraction of absorbed dose taking the delayed route FS (unitless)") # theta10 = 0.622, Table S2

    # ================================================================
    # Saturable plasma protein binding -- Equations 16-19, refitted
    # under the non-diluted binding assumption. All three fixed
    # (Table S2 marks each "Fixed"). The maximum attainable unbound
    # fraction is fumin + fumax = 0.220, against 0.359 for the
    # original-ppb fit.
    # ================================================================
    fumin <- fixed(0.0260); label("Minimum (zero-concentration) unbound fraction in plasma (unitless)") # theta12, Table S2, Fixed
    fumax <- fixed(0.194);  label("Maximum increment in unbound fraction at saturating concentration (unitless)") # theta13, Table S2, Fixed
    cup50 <- fixed(0.814);  label("Unbound plasma concentration giving half the maximal increment in fu (mg/L)")  # theta14 fu50, Table S2, Fixed

    # ---- Food effect (multipliers; Equations 7-9) -----------------
    e_fed_ka     <- 0.0525; label("Fed-state multiplier on ka (unitless)")                         # theta19 Kafed, Table S2
    e_fed_ka2    <- 0.513;  label("Fed-state multiplier on ka2 (unitless)")                        # theta17 Ka2fed, Table S2
    e_fed_fdepot <- 0.802;  label("Fed-state multiplier on total oral bioavailability (unitless)") # theta18 Ftot,fed, Table S2

    # ================================================================
    # Covariate effects. As in Table S1, Table S2 prints these rows as
    # the MULTIPLIER (1 + theta); each value below is the tabulated
    # number minus 1. The five-significant-figure theta27 row (1.00637)
    # is the clearest internal evidence for the multiplier reading.
    # Phase 3 is the reference stratum.
    # ================================================================
    e_alb_cl    <- 0.2;     label("Additive effect of albumin on CL, per g/dL above 4.1 g/dL (1/(g/dL))") # theta20; Table S2 prints 1 + theta20 = 1.2
    e_phase1_cl <- 0.710;   label("Proportional shift in CL for Phase 1 subjects (unitless)")   # theta22; Table S2 prints 1 + theta22 = 1.710
    e_phase2_cl <- 0.707;   label("Proportional shift in CL for Phase 2 subjects (unitless)")   # theta21; Table S2 prints 1 + theta21 = 1.707
    e_phase1_q  <- 0.788;   label("Proportional shift in CLd1 for Phase 1 subjects (unitless)") # theta24; Table S2 prints 1 + theta24 = 1.788
    e_phase2_q  <- 0.192;   label("Proportional shift in CLd1 for Phase 2 subjects (unitless)") # theta23; Table S2 prints 1 + theta23 = 1.192
    e_phase1_vp <- 0.889;   label("Proportional shift in Vp1 for Phase 1 subjects (unitless)")  # theta26; Table S2 prints 1 + theta26 = 1.889
    e_phase2_vp <- 0.28;    label("Proportional shift in Vp1 for Phase 2 subjects (unitless)")  # theta25; Table S2 prints 1 + theta25 = 1.28
    e_wt_vp     <- 0.00637; label("Additive effect of body weight on Vp1, per kg above 78 kg (1/kg)") # theta27; Table S2 prints 1 + theta27 = 1.00637

    # ================================================================
    # Algebraic ELF penetration -- Bian 2024 Equation 21, adopted
    # verbatim from the FDA XENLETA Multi-Discipline Review
    # (NDA 211672/211673, section 16.3.2.4.1 "ELF PK Model", p. 242):
    #
    #   C_ELF = LPR(1 mg/L) * [C_P(t) * 0.0379]^power
    #
    # "where LPR (1 mg/L) is the LPR at a plasma concentration of
    # 1 mg/L, and the power parameter allows the penetration ratio to
    # change with plasma concentration. The plasma concentration was
    # adjusted by 0.0379, fraction unbound of lefamulin in plasma."
    # C_P(t) is the TOTAL plasma concentration; the review's 0.0379 is
    # a separate FIXED LINEAR unbound fraction used only inside this
    # link function, and is deliberately not the saturable fu above.
    #
    # NON-PAPER PROVENANCE -- penratio_elf = 3.45 is BACK-SOLVED, not
    # tabulated. Supplementary Table S2 prints LPR = 2.71 and Power =
    # 0.51, which are byte-identical to Table S1's KELF,in = 2.71 and
    # KELF,out = 0.51 and appear to have been copied across. Celf is
    # directly proportional to LPR, so the two readings differ by the
    # factor 2.71 / 3.45 = 0.786 exactly. With the power held at the
    # tabulated 0.51, LPR = 2.71 undershoots the paper's own higher-ppb
    # ELF exposures for the IV 1 h, oral fasted and oral fed regimens
    # (Table S6: 16.12, 18.10 and 16.49 mg*h/L) by 23.1%, 21.9% and
    # 20.0%, whereas 3.45 lands within -2.1% to +1.9% across
    # all five simulated regimens. Encoding the back-solved value is an explicit operator
    # ruling (sidecar oare_PMC11544005 request-001 q3, option B) that
    # overrides the usual "values from the paper" rule; see the
    # vignette's Assumptions and deviations section.
    # ================================================================
    lpenratio_elf     <- log(3.45); label("Lung penetration ratio at 1 mg/L unbound plasma, ELF (unitless)") # BACK-SOLVED from Table S6; Table S2 prints 2.71 (see comment above)
    pwr_penratio_elf  <- 0.51;      label("Power on the unbound plasma concentration in the ELF penetration function (unitless)") # theta29 Power, Table S2

    # ================================================================
    # Inter-individual variability -- Table S2, reported as omega^2
    # with sqrt shown as %CV (e.g. 2.65 -> 162.8%CV), so the tabulated
    # numbers are VARIANCES and are used directly. eta9 and eta10 are
    # additive on the logit scale (Equations 10, 12); the rest are
    # exponential. As in the original-ppb fit, Table S2 reports no
    # eta5 or eta6 variance, so the exp(eta5) on CLd2 and the
    # exp(eta5) * exp(eta6) on Vp2 of Equations 5-6 are omitted.
    # ================================================================
    etalcl         ~ 0.151;  label("IIV on CL (variance)")                   # eta1, Table S2 (38.9%CV)
    etalvc         ~ 0.195;  label("IIV on Vc (variance)")                   # eta2, Table S2 (44.2%CV)
    etalq          ~ 0.064;  label("IIV on CLd1 (variance)")                 # eta3, Table S2 (25.3%CV)
    etalvp         ~ 0.0875; label("IIV on Vp1 (variance)")                  # eta4, Table S2 (29.6%CV)
    etalka         ~ 2.65;   label("IIV on ka (variance)")                   # eta7, Table S2 (162.8%CV)
    etalka2        ~ 0.264;  label("IIV on ka2 (variance)")                  # eta8, Table S2 (51.4%CV)
    etalogitfdepot ~ 0.791;  label("IIV on logit total oral bioavailability (variance)") # eta9, Table S2 (88.9%CV)
    etalogitfrac   ~ 0.716;  label("IIV on logit delayed-route fraction (variance)")     # eta10, Table S2 (87.2%CV)

    # ================================================================
    # Residual unexplained variability -- Equation 20. Table S2 reports
    # sigma^2; the SDs below are the square roots. sqrt(0.0958) = 0.310
    # and sqrt(0.192) = 0.438 reproduce the table's own 31.0%CV and
    # 43.8%CV. The additive row's parenthetical "0. 0000167 mg/L"
    # simply repeats the variance instead of converting it, unlike
    # Table S1's matching row (0.0000343 -> 0.00586 mg/L); the square
    # root 0.00409 mg/L is used here.
    # ================================================================
    propSd      <- 0.30952;  label("Proportional residual error, plasma (fraction)") # sqrt(sigma^2 eps1 = 0.0958), Table S2
    addSd       <- 0.004087; label("Additive residual error, plasma (mg/L)")         # sqrt(sigma^2 eps2 = 0.0000167), Table S2
    propSd_Celf <- 0.43818;  label("Proportional residual error, ELF (fraction)")    # sqrt(sigma^2 eps3 = 0.192), Table S2
  })

  model({
    # ---- Covariate preparation ------------------------------------
    alb_gdL <- ALB * 0.1

    clphase <- 1 + e_phase1_cl * STUDY_LEFAMULIN_PHASE1 +
                   e_phase2_cl * STUDY_LEFAMULIN_PHASE2
    qphase  <- 1 + e_phase1_q  * STUDY_LEFAMULIN_PHASE1 +
                   e_phase2_q  * STUDY_LEFAMULIN_PHASE2
    vpphase <- 1 + e_phase1_vp * STUDY_LEFAMULIN_PHASE1 +
                   e_phase2_vp * STUDY_LEFAMULIN_PHASE2

    # ---- Individual PK parameters (Equations 1-6) -----------------
    # Equations 3 and 4 as printed omit theta3 and theta4; Table S2
    # supplies them as CLd1 = 187 L/hr and Vp1 = 1300 L.
    cl  <- exp(lcl + etalcl) * (1 + e_alb_cl * (alb_gdL - 4.1)) * clphase
    vc  <- exp(lvc + etalvc)
    q   <- exp(lq + etalq) * qphase
    vp  <- exp(lvp + etalvp) * vpphase * (1 + e_wt_vp * (WT - 78))
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    ka   <- exp(lka + etalka)   * (1 - FED + e_fed_ka  * FED)
    ka2  <- exp(lka2 + etalka2) * (1 - FED + e_fed_ka2 * FED)
    tlag <- exp(ltlag)

    fdepot_typ <- exp(logitfdepot) / (1 + exp(logitfdepot))
    bio        <- fdepot_typ * (1 - FED + e_fed_fdepot * FED)
    fpo        <- log(bio / (1 - bio))
    fdepot     <- exp(fpo + etalogitfdepot) / (1 + exp(fpo + etalogitfdepot))
    frac       <- exp(logitfrac + etalogitfrac) /
                  (1 + exp(logitfrac + etalogitfrac))

    penratio_elf <- exp(lpenratio_elf)

    # ---- Concentrations -------------------------------------------
    Cu <- central / vc
    fu <- fumin + fumax * Cu / (cup50 + Cu)

    # ---- ODE system (Figure 2) ------------------------------------
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

    # ---- Bioavailability and lag time (Equations 13-15) -----------
    f(depot)     <- fdepot * (1 - frac)
    f(depot2)    <- fdepot * frac
    alag(depot)  <- tlag
    alag(depot2) <- tlag

    # ---- Observations and residual error --------------------------
    # Equation 21: the total plasma concentration is scaled to free
    # drug by the FDA review's fixed linear unbound fraction 0.0379
    # BEFORE the power is applied. The constant is hardcoded because it
    # is part of the published link function rather than an estimate of
    # this analysis.
    Cc   <- Cu / fu
    Celf <- penratio_elf * (Cc * 0.0379)^pwr_penratio_elf
    Cc   ~ add(addSd) + prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
