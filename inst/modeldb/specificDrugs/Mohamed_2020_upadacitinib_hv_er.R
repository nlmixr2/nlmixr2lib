Mohamed_2020_upadacitinib_hv_er <- function() {
  description <- paste(
    "Two-compartment population PK model for oral upadacitinib immediate-release",
    "(IR) and extended-release (ER) formulations in healthy adults, fitted to",
    "the phase 1 formulation study M14-680 (Mohamed 2020, Table S4).",
    "Disposition, its covariate effects (sex and creatinine clearance on CL/F,",
    "sex and body weight on Vc/F) and their interindividual variability are",
    "fixed to the final IR model (Mohamed_2020_upadacitinib); only absorption",
    "was estimated. IR: first-order absorption after a lag time. ER: a",
    "first-order input for 44.6% of the absorbed dose and a parallel zero-order",
    "input into the central compartment for the rest, both after a common lag",
    "time, with 74.2% bioavailability relative to IR."
  )
  reference <- paste(
    "Mohamed M-EF, Klunder B, Lacerda AP, Othman AA. Exposure-Response",
    "Analyses for Upadacitinib Efficacy and Safety in the Crohn's Disease",
    "CELEST Study and Bridging to the Extended-Release Formulation. Clin",
    "Pharmacol Ther. 2020;107(3):639-649. doi:10.1002/cpt.1668"
  )
  vignette <- "Mohamed_2020_upadacitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc/F, exponent 0.490 fixed, reference 74 kg (Table S4 footnote b).",
      source_name = "WTKG"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male). Table S4 footnotes a and b set the sex factors to 1 for males.",
      notes = "Female CL/F ratio 0.850 and Vc/F ratio 0.748, both fixed (Table S4). Source SEX.EQ.1 = female.",
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance by the Cockcroft-Gault equation (raw mL/min).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F, exponent 0.313 fixed; the control stream divides by 106.9 mL/min (Table S4 footnote a prints 107).",
      source_name = "CRCL"
    ),
    FORM_UPA_ER = list(
      description = "Upadacitinib extended-release formulation indicator (1 = extended-release tablet, 0 = immediate-release reference).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (immediate release).",
      notes = "Selects the absorption branch. Each ER administration needs two dose records at the same time and amount: one to 'depot' and one to 'central' with rate = -2. Source FORM is 1 = IR, 2 = ER, so FORM_UPA_ER = FORM - 1.",
      source_name = "FORM (1 = IR, 2 = ER)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 71L,
    n_studies = 1L,
    age_range = "19-53 years (mean 36.0)",
    weight_range = "55-101 kg (mean 78.2)",
    sex_female_pct = 14,
    race_ethnicity = c(White = 51, Black = 42, Other = 7),
    disease_state = "Healthy adults",
    dose_range = "Single doses of 12 mg IR and 15 mg ER or 24 mg IR and 30 mg ER (crossover); 15 or 30 mg ER once daily for 7 days; 6 mg IR twice daily vs 15 mg ER once daily and 12 mg IR twice daily vs 30 mg ER once daily (crossover).",
    regions = "United States",
    notes = "Study M14-680 (Table S1, parts 1-5; Table S2 'Study 4' column). Fitted in NONMEM with SAEM; ETA(4)-ETA(8) of the ER parameters were estimated as a 5 x 5 block whose off-diagonal elements are not reported."
  )

  ini({
    # ---- Disposition, fixed to the final IR model (Table S4) ----------------
    lcl <- fixed(log(39.6))
    label("Apparent clearance CL/F, male at CrCL 106.9 mL/min (L/h)") # Table S4 'CL/F (L/h)' = 39.6 fix
    lvc <- fixed(log(141))
    label("Apparent central volume Vc/F, male at 74 kg (L)") # Table S4 'Vc/F (L)' = 141 fix
    lvp <- fixed(log(66.3))
    label("Apparent peripheral volume Vp/F (L)") # Table S4 'Vp/F (L)' = 66.3 fix
    lq <- fixed(log(3.39))
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table S4 'Q/F (L/h)' = 3.39 fix
    e_sexf_cl <- fixed(log(0.850))
    label("Log ratio of CL/F, female vs male (unitless)") # Table S4 'Ratio of CL/F in Females Compared to Males' = 0.850 fix
    e_crcl_cl <- fixed(0.313)
    label("Power exponent of creatinine clearance on CL/F (unitless)") # Table S4 'Exponent for the Effect of CrCL on CL/F' = 0.313 fix
    e_sexf_vc <- fixed(log(0.748))
    label("Log ratio of Vc/F, female vs male (unitless)") # Table S4 'Ratio of Vc/F in Females Compared to Males' = 0.748 fix
    e_wt_vc <- fixed(0.490)
    label("Power exponent of body weight on Vc/F (unitless)") # Table S4 'Exponent for the Effect of Bodyweight on Vc/F' = 0.490 fix

    # ---- IR absorption (estimated in this study) ------------------------------
    lka <- log(9.68)
    label("Absorption rate constant Ka, IR formulation (1/h)") # Table S4 'KA IR (1/h)' = 9.68 (95% CI 4.53-20.7)
    ltlag <- log(0.247)
    label("Absorption lag time, IR formulation (h)") # Table S4 'ALAG1 IR (h)' = 0.247 (95% CI 0.242-0.251)

    # ---- ER absorption ---------------------------------------------------------
    lka_er <- log(0.100)
    label("First-order absorption rate constant Ka, ER formulation (1/h)") # Table S4 'KA ER (1/h)' = 0.100 (95% CI 0.0494-0.203)
    ltlag_er <- log(0.225)
    label("Absorption lag time of both ER input routes (h)") # Table S4 'ALAG1 ER (h)' = 0.225 (95% CI 0.183-0.278)
    logitffo <- log(0.446 / (1 - 0.446))
    label("Logit of the fraction of the ER dose absorbed by the first-order route (unitless)") # Table S4 'Fraction of Extended-release Dose Absorbed through First-Order Process (%)' = 44.6 (95% CI 26.7-64.1)
    ld2 <- log(6.11)
    label("Duration of the zero-order ER input into the central compartment (h)") # Table S4 'Zero-Order Infusion duration D2 (h)' = 6.11 (95% CI 2.37-15.8)
    lfrel_er <- log(0.742)
    label("Bioavailability of the ER formulation relative to IR (fraction)") # Table S4 relative bioavailability (%) = 74.2 (95% CI 69.0-79.8)

    # ---- Interindividual variability ---------------------------------------------
    # Table S4 '%ISV' = SQRT(omega^2) x 100, so variance = (%ISV / 100)^2. The
    # disposition ISVs are the healthy-subject values of the final model
    # (control stream THETA(14) = 0.162 FIX and THETA(12) = 0.134 FIX, the
    # standard deviations of EXP(THETA * ETA) with OMEGA 1 FIX).
    etalcl ~ fixed(0.026244) # Table S4 CL/F '%ISV' = 16 fix; control stream THETA(14) = 0.162 FIX; 0.162^2
    etalvc ~ fixed(0.017956) # Table S4 Vc/F '%ISV' = 13 fix; control stream THETA(12) = 0.134 FIX; 0.134^2
    etalka ~ 1.4884 # Table S4 'KA IR' %ISV = 122; 1.22^2
    etaltlag ~ 0.00011236 # Table S4 'ALAG1 IR' %ISV = 1.06; 0.0106^2
    etalka_er ~ 1.2321 # Table S4 'KA ER' %ISV = 111; 1.11^2
    etaltlag_er ~ 0.218089 # Table S4 'ALAG1 ER' %ISV = 46.7; 0.467^2
    etalogitffo ~ 1.49 # Table S4 'Variance of the ISV for the Fraction of Extended-release Dose Absorbed through First-Order Process' = 1.49
    etald2 ~ 1.4641 # Table S4 'Zero-Order Infusion duration D2' %ISV = 121; 1.21^2
    etalfrel_er ~ 0.0289 # Table S4 relative bioavailability %ISV = 17.0; 0.170^2

    # ---- Residual error ------------------------------------------------------------
    propSd <- 0.288
    label("Proportional residual error SD (fraction)") # Table S4 'Proportional Residual Error SD (%)' = 28.8
    addSd <- 0.122
    label("Additive residual error SD (ng/mL)") # Table S4 'Additive Residual Error SD (ng/mL)' = 0.122
  })

  model({
    ref_crcl <- 106.9
    ref_wt <- 74

    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF + e_crcl_cl * log(CRCL / ref_crcl))
    vc <- exp(lvc + etalvc + e_sexf_vc * SEXF + e_wt_vc * log(WT / ref_wt))
    vp <- exp(lvp)
    q <- exp(lq)

    # Control stream: KA = EXP(MU_3 + ETA(3)); IF(FORM.EQ.2) KA = EXP(MU_5 +
    # ETA(5)); ALAG1 = EXP(MU_9 + ETA(9)); IF(FORM.EQ.2) ALAG1 = ALAG2 = LAGT2;
    # F1 = Frel * Frac and F2 = Frel * (1 - Frac) for ER; D2 = EXP(MU_4 + ETA(4)).
    ka_ir <- exp(lka + etalka)
    tlag_ir <- exp(ltlag + etaltlag)
    ka_er <- exp(lka_er + etalka_er)
    tlag_er <- exp(ltlag_er + etaltlag_er)
    ffo_er <- expit(logitffo + etalogitffo)
    frel_er <- exp(lfrel_er + etalfrel_er)

    ka <- ka_er * FORM_UPA_ER + ka_ir * (1 - FORM_UPA_ER)
    tlag <- tlag_er * FORM_UPA_ER + tlag_ir * (1 - FORM_UPA_ER)
    frel <- frel_er * FORM_UPA_ER + (1 - FORM_UPA_ER)
    ffo <- ffo_er * FORM_UPA_ER + (1 - FORM_UPA_ER)
    d2 <- exp(ld2 + etald2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # IR doses: one record to depot. ER doses: one record to depot plus one to
    # central with rate = -2 at the same time and amount.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- frel * ffo
    f(central) <- frel * (1 - ffo)
    dur(central) <- d2
    alag(depot) <- tlag
    alag(central) <- tlag

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
