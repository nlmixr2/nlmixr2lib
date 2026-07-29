vandenBerg_2021_uprifosbuvir_pbpk <- function() {
  description <- paste(
    "PBPK-PD (minimal PBPK; well-stirred liver + Qgut gut wall + intrahepatic",
    "UXP/CXP nucleotide cycle) for uprifosbuvir (MK-3682), an HCV NS5B",
    "nucleoside polymerase inhibitor, and its plasma metabolites M5 and M6, with",
    "a sigmoid Emax link from intrahepatic UXP to viral inhibition (epsilon).",
    "22 ODE compartments (gut absorption + transit chains, hepatic + portal-vein +",
    "central + two peripherals for parent, intrahepatic UXP/CXP cycle, M5 and M6",
    "central + M5 gut transit chain, pseudo-M4 gut pool). All parameters fixed at",
    "the final NONMEM point estimates (Tables 1-6 of van den Berg 2021).",
    "Supports i.v. (25 mg), tablet oral (150-750 mg) and capsule oral (50-400 mg)",
    "dosing in HV and HCV patients, with optional itraconazole DDI encoded via",
    "covariate multipliers. All amounts are in nmol; concentrations in nmol/L.",
    "Framework adapted from the Brill et al. 2016 midazolam mPBPK with metabolite."
  )
  reference <- paste(
    "van den Berg P, Gao W, Ahsman MJ, Arrington L, Kesisoglou F, Miller R,",
    "Post TM, Rizk ML (2021). Understanding effect site pharmacology of",
    "uprifosbuvir, a hepatitis C virus nucleoside inhibitor: Case study of a",
    "multidisciplinary modeling approach in drug development.",
    "CPT Pharmacometrics Syst Pharmacol 10(7):658-670.",
    "doi:10.1002/psp4.12644.",
    sep = " "
  )
  vignette <- "vandenBerg_2021_uprifosbuvir"
  units    <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  paper_specific_compartments <- c(
    "gut", "gut_m5", "gut_m6",
    "liver", "portal_vein",
    "uxp", "cxp",
    "central_m5", "central_m6",
    "transit_cxp_m5",
    "transit_gut_m5_1", "transit_gut_m5_2",
    "transit_gut_m5_3", "transit_gut_m5_4",
    "transit_gut_par",
    "pseudo_m4_gut"
  )

  covariateData <- list(
    DOSE = list(
      description        = paste(
        "Nominal dose level in mg per administration. Enters the model as a per-",
        "dose covariate: it centres the dose-dependent slope on absorption rate",
        "KA2 (all formulations) and, for the capsule formulation only, on",
        "absorption rate KA1 and on the gut-M6 elimination rate KelM6g. Reference",
        "dose is 150 mg; the effect is `exp(coeff * (DOSE - 150))`. Set DOSE per",
        "dose event in the event table."
      ),
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference dose 150 mg per van den Berg 2021 supplement Table 1 note (a)",
        "and NONMEM code KA2 = TVKA2 * EXP(DOSF * (DOSE - 150))."
      ),
      source_name        = "DOSE"
    ),
    FORM_CAPSULE = list(
      description        = "Capsule formulation indicator (1 = capsule, 0 = tablet or IV).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tablet or IV = reference)",
      notes              = paste(
        "The capsule formulation used in Phase 1 studies P001/P003 had different",
        "absorption kinetics: 91% of the capsule dose is absorbed via the fast",
        "KA1 route (vs 67% for the tablet), KA1 is 2.02-fold higher, apparent",
        "bioavailability is 1.17-fold higher, and gut-M6 elimination follows a",
        "separate parameter set. Set FORM_CAPSULE = 0 for tablet and IV records."
      ),
      source_name        = "FORM (2 = capsule; recode)"
    ),
    ROUTE_IV = list(
      description        = "Intravenous administration indicator (1 = IV, 0 = oral).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral = reference)",
      notes              = paste(
        "When 1, the absorbed fractions F1/F2 (fast/slow oral) are set to 0 so",
        "the depot compartments contribute no mass; the dose must be placed in",
        "the central compartment via cmt = 'central' on the dose event. When 0,",
        "the oral absorption fractions apply per the FORM_CAPSULE / itraconazole",
        "logic in the model. Corresponds to FORM = 10 in the NONMEM control",
        "stream (i.v.)."
      ),
      source_name        = "FORM (10 = IV; recode to binary)"
    ),
    CONMED_ITRACONAZOLE = list(
      description        = "Co-administered itraconazole indicator (1 = with itraconazole, 0 = alone).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (uprifosbuvir alone = reference)",
      notes              = paste(
        "Itraconazole 200 mg QD reduces uprifosbuvir gut CYP3A4 / P-gp metabolism,",
        "delaying and reducing early gut-mediated M6 formation, and it delays the",
        "second uprifosbuvir absorption peak. Encoded as multipliers on CLiG,",
        "CLintH, ALAG2, KA2, Ktr1, KM6g, FrM6g, and KelM6g per Table 3."
      ),
      source_name        = "INTR"
    ),
    HCV_POS = list(
      description        = "HCV infection status (1 = HCV-seropositive patient, 0 = healthy volunteer).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer = reference)",
      notes              = paste(
        "In HCV-positive patients M5 clearance is 28.4% lower than in HV",
        "(multiplicative factor 0.716 on CLM5). Attributed to a higher prevalence",
        "of impaired hepatic function among chronic HCV patients."
      ),
      source_name        = "HCVSTAT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 220L,
    n_studies      = 6L,
    age_range      = "adults",
    disease_state  = paste(
      "Pooled healthy volunteers and chronic HCV-infected patients across six",
      "Phase 1/1-2a studies (PN001, PN002, PN003, PN010, PN013, PN026).",
      "HCV cohort included multiple HCV genotypes (GT1a/1b/2/3), no cirrhosis,",
      "no HIV co-infection."
    ),
    dose_range     = paste(
      "Uprifosbuvir 25 mg IV (single dose, PN026); tablet 150-750 mg single or",
      "multiple oral doses (PN001/PN002/PN003/PN010/PN013/PN026); capsule",
      "50-400 mg single or multiple oral doses (PN001/PN003). Itraconazole",
      "200 mg QD co-administered in the DDI arms (PN001 Group F, PN010)."
    ),
    regions        = "PN001-026 multinational Phase 1 programme",
    notes          = paste(
      "See van den Berg 2021 Table 1 and supplementary Table S1 for the full",
      "per-study breakdown (25 mg IV: N=14; tablet SD/MD: N=88; capsule SD/MD:",
      "N=111; itraconazole DDI: N=19). The model was fit simultaneously to",
      "uprifosbuvir + M5 + M6 plasma concentrations; a K-PD viral inhibition",
      "estimate (epsilon) from a prior K-PD analysis of Phase 1b PN001 data was",
      "used as the PD endpoint in the final PBPK-PD step (N=88 HCV patients)."
    )
  )

  ini({
    # ==================== STRUCTURAL PARAMETERS ====================
    # All THETAs fixed at the point estimates of the final minimal PBPK-PD
    # model (van den Berg 2021 supplement Tables 1, 3, 5, 6). NONMEM index
    # is retained in the trailing comment to make provenance auditable.

    # ---------- Gut absorption / bioavailability ----------
    # FDOS1 is estimated as the pre-logit for the fast-fraction of dose;
    # the tablet reference fraction on natural scale is 0.669 (Table 1
    # note a: TVF1 = exp(FDOS1) / (1 + exp(FDOS1)) so for TH1 = 0.704
    # this gives 0.669). For the capsule the logit-scale value is
    # 0.704 * 3.23 = 2.274, corresponding to 0.907 on natural scale.
    fdos1_logit    <- fixed(0.704);    label("Logit-scale fast-fraction of dose (tablet baseline); F1 = expit(fdos1_logit)")   # TH(1); Table 1
    lka1_150       <- fixed(log(0.670));   label("Uprifosbuvir absorption rate KA1 at 150 mg (log 1/h)")                      # TH(2); Table 1
    lka2_150       <- fixed(log(0.497));   label("Uprifosbuvir absorption rate KA2 at 150 mg (log 1/h)")                      # TH(3); Table 1
    e_dose_ka2     <- fixed(-0.00178); label("Dose slope on KA2 (linear on log-KA2, 1/mg)")                                   # TH(4); Table 1
    ltlag_fast     <- fixed(log(0.125));   label("Lag time for the fast absorption peak (log h)")                             # TH(38); Table 1
    ltlag_slow     <- fixed(log(1.87));    label("Lag time for the slow absorption peak (log h)")                             # TH(5); Table 1
    lktr1          <- fixed(log(0.0598));  label("Rate constant into transit gut parent (log 1/h)")                           # TH(6); Table 1

    # ---------- Gut wall metabolism ----------
    lcl_int_g      <- fixed(log(15.5));   label("Uprifosbuvir intrinsic gut wall clearance (log L/h, unbound-blood basis)")   # TH(8); Table 1
    frm5_logit     <- fixed(-1.56);   label("Pre-logit for gut fraction to M5 (tablet baseline); F_M5 = exp(x)/SUM")          # TH(9); Table 1 (natural scale 0.110)
    frm4_logit     <- fixed(-0.345);  label("Pre-logit for gut fraction to pseudo-M4 (tablet baseline); F_M4 = exp(x)/SUM")   # TH(39); Table 1 (natural scale 0.369)
    # FrM6g = 1 - FrM5g - FrM4g = 0.521 (Table 1 note f).
    lkm4g          <- fixed(log(0.119));  label("Pseudo-M4 gut absorption rate (log 1/h)")                                    # TH(40); Table 1
    lkm6g          <- fixed(log(0.106));  label("Gut M6 absorption rate constant (log 1/h)")                                  # TH(10); Table 1
    lk2m5g         <- fixed(log(0.402));  label("Rate gut M5 -> M5 transit chain (log 1/h; K3M5g = K2M5g fixed)")             # TH(12); Table 1
    lkelm6g_50000  <- fixed(log(6.10));   label("Conc-dep gut M6 elimination rate at 50,000 nmol gut M6 (log 1/h; tablet)")   # TH(11); Table 1
    e_pow_kelm6g   <- fixed(1.39);        label("Power exponent on the conc-dep gut M6 elimination (tablet)")                 # TH(41); Table 1

    # ---------- Systemic (parent) ----------
    lcl_int_h      <- fixed(log(82.6));   label("Uprifosbuvir intrinsic hepatic clearance (log L/h, unbound-blood basis)")    # TH(14); Table 1
    lvc            <- fixed(log(3.92));   label("Uprifosbuvir central volume of distribution (log L)")                        # TH(15); Table 1
    lvp            <- fixed(log(7.19));   label("Uprifosbuvir first peripheral volume of distribution (log L; V_pp2 = V_pp)") # TH(16); Table 1
    lq             <- fixed(log(17.3));   label("Uprifosbuvir first inter-compartmental clearance (log L/h)")                 # TH(17); Table 1
    lq2            <- fixed(log(1.11));   label("Uprifosbuvir second inter-compartmental clearance (log L/h)")                # TH(19); Table 1

    # ---------- UXP - CXP intrahepatic cycle ----------
    lkunuc1        <- fixed(log(13.3));   label("Rate constant UXP -> M6 (log; NONMEM parameterises as THETA/1000)")          # TH(20); Table 1 note i (natural 0.0133 1/h)
    lkucxp         <- fixed(log(0.785));  label("Rate constant UXP -> CXP (log 1/h; kCUXP = kUCXP fixed)")                    # TH(21); Table 1
    lkeluxp_50000  <- fixed(log(0.481));  label("Conc-dep UXP -> CXP additional rate at 50,000 nmol UXP (log 1/h)")           # TH(43); Table 1
    e_pow_keluxp   <- fixed(1.34);        label("Power exponent on the conc-dep UXP -> CXP rate")                             # TH(44); Table 1
    lkcnuc         <- fixed(log(33.8));   label("Rate constant CXP -> M5 transit (log; NONMEM parameterises as THETA/1000)")  # TH(23); Table 1 note k (natural 0.0338 1/h)
    lkcnuc2        <- fixed(log(0.122));  label("Rate constant transit CXP-M5 -> central M5 (log 1/h)")                       # TH(37); Table 1

    # ---------- M5 and M6 systemic clearance ----------
    lcl_m6         <- fixed(log(2.60));   label("M6 systemic clearance (log L/h; V_M6 = V_c fixed)")                          # TH(24); Table 1
    lcl_m5         <- fixed(log(79.5));   label("M5 systemic clearance in HV (log L/h; V_M5 = V_c fixed)")                    # TH(26); Table 1

    # ---------- HCV covariate on CL_M5 ----------
    e_hcv_cl_m5    <- fixed(0.716);       label("HCV factor on CL_M5 (multiplicative; 28.4% lower in HCV patients)")          # TH(53); Table 3

    # ---------- Capsule formulation covariates ----------
    e_capsule_ka1     <- fixed(2.02);     label("Capsule factor on KA1 (multiplicative; tablet reference)")                   # TH(45); Table 3
    e_capsule_fdos1   <- fixed(3.23);     label("Capsule factor on FDOS1 (logit scale; corresponds to F1 = 0.91 for capsule)")# TH(46); Table 3
    e_capsule_bioav   <- fixed(1.17);     label("Capsule factor on bioavailability (multiplicative on F1 and F2)")            # TH(54); Table 3
    lkelm6g_50000_cap <- fixed(log(5.76));    label("Capsule conc-dep gut M6 elim at 50,000 nmol gut M6 (log 1/h)")           # TH(63); Table 3
    e_dose_kelm6g_cap <- fixed(-0.00420); label("Capsule dose slope on KelM6g (linear on log; 1/mg)")                         # TH(64); Table 3
    e_dose_ka1_cap    <- fixed(-0.00189); label("Capsule dose slope on KA1 (linear on log; 1/mg)")                            # TH(65); Table 3
    e_capsule_tlag_fast <- fixed(1.65);   label("Capsule factor on ALAG1 (multiplicative)")                                   # TH(66); Table 3

    # ---------- Itraconazole DDI covariates ----------
    e_itra_cl_int_g <- fixed(0.444);      label("Itraconazole factor on CL_int,G (multiplicative)")                           # TH(55); Table 3
    e_itra_cl_int_h <- fixed(0.797);      label("Itraconazole factor on CL_int,H (multiplicative)")                           # TH(56); Table 3
    e_itra_tlag_slow <- fixed(0.145);     label("Itraconazole factor on ALAG2 (multiplicative)")                              # TH(57); Table 3
    e_itra_ka2      <- fixed(0.499);      label("Itraconazole factor on KA2 (multiplicative)")                                # TH(58); Table 3
    e_itra_ktr1     <- fixed(0.166);      label("Itraconazole factor on Ktr1 (multiplicative)")                               # TH(59); Table 3
    e_itra_km6g     <- fixed(0.294);      label("Itraconazole factor on KM6g (multiplicative)")                               # TH(60); Table 3
    e_itra_frm6g    <- fixed(0.174);      label("Itraconazole factor on FrM6g (multiplicative)")                              # TH(61); Table 3
    e_itra_kelm6g   <- fixed(0.0290);     label("Itraconazole factor on KelM6g (multiplicative)")                             # TH(62); Table 3

    # ---------- PBPK-PD viral inhibition (Table 6) ----------
    # UXP -> viral inhibition: EFFW = 1 + Emax * UXP^Hill / (EC50^Hill + UXP^Hill)
    # NONMEM: EC50 = exp(THETA(68) + ETA(1)) so typical EC50 in UXP units is
    # ~exp(9.77) = 17500 nmol/L (Table 6 note a).
    lec50          <- fixed(9.77);        label("Log EC50 for viral inhibition (log nmol/L UXP; typical EC50 = 17,500)")      # TH(68); Table 6 note a
    lhill          <- fixed(log(2.80));   label("Log Hill coefficient for viral inhibition (unitless)")                       # Table 6
    lemax          <- fixed(log(1.0));    label("Log Emax for viral inhibition (fixed to 1; unitless multiplier)")            # Table 6

    # ==================== IIV (Table 4) ====================
    # Full 15-eta variance-covariance matrix (with the kCnuc/CLM5 block) is
    # reported in Table 4; the four etas with the greatest information content
    # (CL_int,H, Vc, CL_M6, CL_M5) are encoded here as fixed() variances so the
    # model supports population VPC-type simulations. The remaining 11 etas
    # (KA2, FDOS1, CL_int,G, KelM6g, FrM5g, FrM4g, KA1_capsule, kUnuc1,
    # F_capsule, kCnuc, ALAG1 factor) are not encoded; downstream users needing
    # the full block should extend this file per Table 4.
    etalcl_int_h ~ fixed(0.00693)                                                                                             # Table 4 ETA(1); omega^2
    etalvc       ~ fixed(0.137)                                                                                               # Table 4 ETA(2); omega^2
    etalcl_m6    ~ fixed(0.0222)                                                                                              # Table 4 ETA(3); omega^2
    etalcl_m5    ~ fixed(0.0722)                                                                                              # Table 4 ETA(5); omega^2

    # ==================== Residual error ====================
    # Table 5: sigma^2 values converted to SD = sqrt(variance). Reported for
    # each formulation x route x analyte combination; the tablet oral values
    # for uprifosbuvir, M5 and M6 are encoded as the default outputs (the
    # dominant use case in the paper's Figures 3 and 4). IV, capsule and the
    # additive terms are reproduced in the vignette Errata.
    propSd        <- fixed(sqrt(0.155));  label("Uprifosbuvir tablet oral proportional residual SD (fraction)")               # TH(29); Table 5
    addSd         <- fixed(sqrt(3.36));   label("Uprifosbuvir tablet oral additive residual SD (nmol/L)")                     # TH(30); Table 5
    propSd_Cc_m5  <- fixed(sqrt(0.0294)); label("M5 tablet oral proportional residual SD (fraction)")                         # TH(32); Table 5
    propSd_Cc_m6  <- fixed(sqrt(0.0153)); label("M6 tablet oral proportional residual SD (fraction)")                         # TH(35); Table 5
    addSd_Cc_m6   <- fixed(sqrt(342));    label("M6 tablet oral additive residual SD (nmol/L)")                               # TH(36); Table 5
  })

  model({
    # =============== PHYSIOLOGICAL PARAMETERS (fixed) ===============
    # Values from van den Berg 2021 supplement Table 2 / NONMEM code.
    # All flows in L/h; all reference volumes fixed to 1 L (the compartments
    # are book-keeping volumes for the well-stirred / Qgut extraction models).
    qhep   <- 88.0                    # Hepatic blood flow (Davies & Morris 1993)
    qpv    <- 0.75 * qhep             # Portal vein blood flow (Williams & Leggett 1989)
    qha    <- 0.25 * qhep             # Hepatic artery blood flow (Williams & Leggett 1989)
    vpv    <- 1.0                     # Portal vein volume (fixed to 1 L)
    vh     <- 1.0                     # Liver volume (fixed to 1 L)
    vgw    <- 1.0                     # Gut wall volume (fixed to 1 L)
    fu_b   <- 0.577                   # Uprifosbuvir plasma unbound fraction (Sponsor)
    bpr    <- 0.730                   # Uprifosbuvir blood:plasma ratio (Sponsor)
    fug    <- 1.0                     # Uprifosbuvir gut unbound fraction (Yang 2007)
    qin    <- 0.40 * qhep             # Small-intestine blood flow (Williams 1989)
    qmu    <- 0.80 * qin              # Mucosa blood flow (Yang 2007)
    qvi    <- 0.60 * qmu              # Villous blood flow (Yang 2007)
    cl_per <- 2.9                     # Enterocyte permeability CL (L/h; Sponsor)

    # =============== INDIVIDUAL PK PARAMETERS ===============
    # Formulation and DDI multipliers switched by binary covariates.
    # (Reference: TABLET oral, HV, alone, DOSE = 150 mg.)

    # Absorption fast fraction (F1) via logit; F1 baseline = expit(fdos1_logit).
    # Capsule multiplies the pre-logit scalar by e_capsule_fdos1.
    fdos_scaled <- fdos1_logit * (1 + FORM_CAPSULE * (e_capsule_fdos1 - 1))
    f1_oral     <- exp(fdos_scaled) / (1 + exp(fdos_scaled))
    f2_oral     <- 1 - f1_oral

    # Itraconazole shifts absorption entirely to the slow route (per NONMEM:
    # F1 = 0 and F2 = 1 when INTR == 1). IV administration bypasses both
    # oral routes (F1 = F2 = 0; dose goes to `central` via the event table).
    f1 <- (1 - ROUTE_IV) * ((1 - CONMED_ITRACONAZOLE) * f1_oral) *
          (1 + FORM_CAPSULE * (e_capsule_bioav - 1))
    f2 <- (1 - ROUTE_IV) * ((1 - CONMED_ITRACONAZOLE) * f2_oral + CONMED_ITRACONAZOLE * 1.0) *
          (1 + FORM_CAPSULE * (e_capsule_bioav - 1))

    # KA1 baseline is estimated at 150 mg. Capsule scales KA1 by the factor
    # e_capsule_ka1 and adds a dose slope. Tablet has no dose slope on KA1.
    ka1_tab <- exp(lka1_150)
    ka1_cap <- exp(lka1_150) * e_capsule_ka1 * exp(e_dose_ka1_cap * (DOSE - 150))
    ka1     <- (1 - FORM_CAPSULE) * ka1_tab + FORM_CAPSULE * ka1_cap

    # KA2 dose slope applies for both formulations; ITRA scales KA2 by factor.
    ka2_base <- exp(lka2_150) * exp(e_dose_ka2 * (DOSE - 150))
    ka2      <- ka2_base * (1 + CONMED_ITRACONAZOLE * (e_itra_ka2 - 1))

    # Lag times
    tlag_fast_tab <- exp(ltlag_fast)
    tlag_fast_cap <- exp(ltlag_fast) * e_capsule_tlag_fast
    tlag_fast     <- (1 - FORM_CAPSULE) * tlag_fast_tab + FORM_CAPSULE * tlag_fast_cap
    tlag_slow     <- exp(ltlag_slow) * (1 + CONMED_ITRACONAZOLE * (e_itra_tlag_slow - 1))

    # Transit rates gut -> gut wall
    ktr1 <- exp(lktr1) * (1 + CONMED_ITRACONAZOLE * (e_itra_ktr1 - 1))
    ktr2 <- ktr1 * 1.0    # K3M5g = K2M5g in NONMEM (TH(13) = 1 FIX); labelled ktr2 here for symmetry

    # Gut wall extraction (Qgut hybrid model; Yang 2007)
    cl_int_g <- exp(lcl_int_g + etalcl_int_h * 0) * (1 + CONMED_ITRACONAZOLE * (e_itra_cl_int_g - 1))
    qgut     <- (qvi * cl_per) / (qvi + cl_per)
    egut     <- (fug * cl_int_g) / (qgut + fug * cl_int_g)
    fracg    <- 1 - egut

    # Fractions to M4, M5, M6 in the gut (softmax over pre-logit; FrM6g = 1/SUM)
    tsum   <- 1 + exp(frm5_logit) + exp(frm4_logit)
    frm5g  <- exp(frm5_logit) / tsum
    frm4g  <- exp(frm4_logit) / tsum
    frm6g_base <- 1 / tsum
    frm6g  <- frm6g_base * (1 + CONMED_ITRACONAZOLE * (e_itra_frm6g - 1))

    # Gut M4, M5, M6 absorption/conversion rates
    km4g   <- exp(lkm4g)
    km6g   <- exp(lkm6g) * (1 + CONMED_ITRACONAZOLE * (e_itra_km6g - 1))
    k2m5g  <- exp(lk2m5g)
    k3m5g  <- k2m5g * 1.0    # K3M5g = K2M5g fixed

    # Gut-M6 conc-dep elimination rate (formulation-specific)
    tkelm6g_tab <- exp(lkelm6g_50000) *
      (1 + CONMED_ITRACONAZOLE * (e_itra_kelm6g - 1))
    tkelm6g_cap <- exp(lkelm6g_50000_cap) *
      exp(e_dose_kelm6g_cap * (DOSE - 150))
    tkelm6g     <- (1 - FORM_CAPSULE) * tkelm6g_tab + FORM_CAPSULE * tkelm6g_cap
    pow_kelm6g  <- e_pow_kelm6g

    # Hepatic extraction (parent, well-stirred model)
    cl_int_h_tv <- exp(lcl_int_h + etalcl_int_h) *
      (1 + CONMED_ITRACONAZOLE * (e_itra_cl_int_h - 1))
    eh   <- (cl_int_h_tv * fu_b) / (qhep + cl_int_h_tv * fu_b)
    fhp  <- 1 - eh

    # Volumes and inter-compartmental rates for parent
    vc  <- exp(lvc + etalvc)
    vp  <- exp(lvp)
    vp2 <- vp * 1.0                # V_pp2 = V_pp fixed (TH(18) = 1 FIX)
    q   <- exp(lq)
    q2  <- exp(lq2)
    kcp  <- q  / vc
    kpc  <- q  / vp
    kcp2 <- q2 / vc
    kp2c <- q2 / vp2

    # UXP <-> CXP intrahepatic cycle
    kunuc1 <- exp(lkunuc1) / 1000    # NONMEM parameterises as THETA/1000
    kucxp  <- exp(lkucxp)
    kcuxp  <- kucxp * 1.0            # kCUXP = kUCXP fixed
    kcnuc  <- exp(lkcnuc) / 1000     # NONMEM parameterises as THETA/1000
    kcnuc2 <- exp(lkcnuc2)
    tkeluxp <- exp(lkeluxp_50000)
    pow_keluxp <- e_pow_keluxp

    # Metabolite clearances (V = Vc for both M5 and M6; TH(25) = TH(27) = 1 FIX)
    hcv_cl_m5 <- 1.0 + HCV_POS * (e_hcv_cl_m5 - 1)
    cl_m6 <- exp(lcl_m6 + etalcl_m6)
    cl_m5 <- exp(lcl_m5 + etalcl_m5) * hcv_cl_m5
    v_m6  <- vc
    v_m5  <- vc
    kem6  <- cl_m6 / v_m6
    kem5  <- cl_m5 / v_m5

    # =============== ODE SYSTEM (22 compartments) ===============
    # Compartment map (van den Berg 2021 supplement, $MODEL block):
    #   1  depot            (fast oral depot;      DEPOT1)
    #   2  depot2           (slow oral depot;      DEPOT2)
    #   3  gut              (gut wall parent;      GUT)
    #   4  gut_m6           (gut M6 transit;       GUTM6)
    #   5  liver            (liver parent;         LIVER)
    #   6  central          (parent central;       CENTPAR)
    #   7  peripheral1      (parent peripheral 1;  PERIPAR)
    #   8  uxp              (intrahepatic UXP;     UXP)
    #   9  central_m6       (M6 central;           M6CENTR)
    #  10  transit_cxp_m5   (CXP -> M5 transit;    CTRM5)  (NONMEM slot 11)
    #  11  cxp              (intrahepatic CXP;     CXP)   (NONMEM slot 12)
    #  12  central_m5       (M5 central;           M5CENTR)(NONMEM slot 13)
    #  13  portal_vein      (portal vein parent;   PVPAR) (NONMEM slot 14)
    #  14  gut_m5           (gut M5;               GUTM5) (NONMEM slot 15)
    #  15  transit_gut_m5_1 (M5 gut transit 1;     TRGTM51)
    #  16  transit_gut_m5_2 (M5 gut transit 2;     TRGTM52)
    #  17  transit_gut_m5_3 (M5 gut transit 3;     TRGTM53)
    #  18  transit_gut_m5_4 (M5 gut transit 4;     TRGTM54)
    #  19  transit_gut_par  (parent gut transit;   TRPAR1)
    #  20  peripheral2      (parent peripheral 2;  PERIPAR2)
    #  21  pseudo_m4_gut    (pseudo-M4 gut;        PSM4GUT)
    #
    # NONMEM slot 10 (DUMMY) is unused (d/dt = 0) and omitted here.
    # Amounts in nmol; concentrations in nmol/L.

    # Conc-dependent kinetics (guard against t=0 zero-amount to avoid 0^POW=NaN)
    cm6_gut  <- gut_m6
    kelm6g   <- 0
    if (cm6_gut > 0) kelm6g <- tkelm6g * (cm6_gut / 50000)^pow_kelm6g

    c_uxp    <- uxp
    keluxp   <- 0
    if (c_uxp > 0) keluxp <- tkeluxp * (c_uxp / 50000)^pow_keluxp

    d/dt(depot)  <- -ka1 * depot
    d/dt(depot2) <- -ka2 * depot2 - ktr1 * depot2
    d/dt(gut)    <- ka1 * depot + ka2 * depot2 + ktr2 * transit_gut_par -
                    fracg * (qvi / vgw) * gut - egut * (cl_int_g / vgw) * gut
    d/dt(gut_m6) <- frm6g * egut * (cl_int_g / vgw) * gut -
                    km6g * gut_m6 - kelm6g * gut_m6
    d/dt(liver)  <- (qpv / vpv) * portal_vein + (qha / vc) * central -
                    eh * (cl_int_h_tv / vh) * liver -
                    fhp * (qhep / vh) * liver
    d/dt(central) <- -(qpv / vc) * central - (qha / vc) * central +
                    fhp * (qhep / vh) * liver -
                    kcp * central + kpc * peripheral1 -
                    kcp2 * central + kp2c * peripheral2
    d/dt(peripheral1) <- kcp * central - kpc * peripheral1
    d/dt(uxp)    <- eh * (cl_int_h_tv / vh) * liver -
                    kunuc1 * uxp - kucxp * uxp + kcuxp * cxp +
                    km4g * pseudo_m4_gut - keluxp * uxp
    d/dt(central_m6) <- kunuc1 * uxp + km6g * gut_m6 - kem6 * central_m6
    d/dt(transit_cxp_m5) <- kcnuc * cxp - kcnuc2 * transit_cxp_m5
    d/dt(cxp)    <- kucxp * uxp - kcuxp * cxp - kcnuc * cxp + keluxp * uxp
    d/dt(central_m5) <- kcnuc2 * transit_cxp_m5 + k3m5g * transit_gut_m5_4 -
                    kem5 * central_m5
    d/dt(portal_vein) <- (qpv / vc) * central - (qpv / vpv) * portal_vein +
                    fracg * (qvi / vgw) * gut
    d/dt(gut_m5) <- frm5g * egut * (cl_int_g / vgw) * gut - k2m5g * gut_m5
    d/dt(transit_gut_m5_1) <- k2m5g * gut_m5 - k3m5g * transit_gut_m5_1
    d/dt(transit_gut_m5_2) <- k3m5g * transit_gut_m5_1 - k3m5g * transit_gut_m5_2
    d/dt(transit_gut_m5_3) <- k3m5g * transit_gut_m5_2 - k3m5g * transit_gut_m5_3
    d/dt(transit_gut_m5_4) <- k3m5g * transit_gut_m5_3 - k3m5g * transit_gut_m5_4
    d/dt(transit_gut_par) <- ktr1 * depot2 - ktr2 * transit_gut_par
    d/dt(peripheral2) <- kcp2 * central - kp2c * peripheral2
    d/dt(pseudo_m4_gut) <- frm4g * egut * (cl_int_g / vgw) * gut -
                    km4g * pseudo_m4_gut

    # Bioavailabilities apply on the fast (depot) and slow (depot2) routes
    f(depot)  <- f1
    f(depot2) <- f2

    # Absorption lag times
    alag(depot)  <- tlag_fast
    alag(depot2) <- tlag_slow

    # =============== OBSERVATIONS ===============
    Cc     <- central     / vc           # uprifosbuvir plasma concentration (nmol/L)
    Cc_m5  <- central_m5  / v_m5         # M5 plasma concentration (nmol/L)
    Cc_m6  <- central_m6  / v_m6         # M6 plasma concentration (nmol/L)

    # Viral inhibition scaling factor from UXP (van den Berg 2021 Eq. 1;
    # Table 6). eff = 1 + Emax * UXP^Hill / (EC50^Hill + UXP^Hill).
    ec50   <- exp(lec50)
    hill   <- exp(lhill)
    emax   <- exp(lemax)
    eff    <- 1 + emax * (uxp^hill) / ((ec50^hill) + (uxp^hill))

    # Multi-output residual error. rxode2 auto-injects DVID slots for the
    # three algebraic observables (Cc, Cc_m5, Cc_m6) AFTER the 21 ODE-state
    # slots -- so DVID 1 -> compartment 22 (Cc), 2 -> 23 (Cc_m5), 3 -> 24
    # (Cc_m6). When building an event table for this model use numeric
    # `cmt = 1L` on the dose row (targets the `depot` ODE state) and either
    # numeric `cmt = 22L/23L/24L` on observation rows to select the specific
    # DV, or set `dvid` per row. See the vignette for a worked example.
    Cc    ~ add(addSd)       + prop(propSd)
    Cc_m5 ~ prop(propSd_Cc_m5)
    Cc_m6 ~ add(addSd_Cc_m6) + prop(propSd_Cc_m6)
  })
}
