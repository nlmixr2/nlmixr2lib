Kim_2025_infliximab_ternant <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in adults with inflammatory bowel disease, with anti-drug-antibody-stratified clearance and an additively decomposed central volume (Ternant 2008 model, as re-coded and externally validated by Kim 2025 in Korean IBD patients)"
  reference <- paste(
    "Kim Y, Baek SH, Jang IJ, Chung JY. Model-Informed Precision Dosing of",
    "Infliximab in Korean Inflammatory Bowel Disease Patients: External",
    "Validation of Population Pharmacokinetic Models.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14:1682-1694.",
    "doi:10.1002/psp4.70089.",
    "The structural model and parameter values were originally developed by",
    "Ternant D, Aubourg A, Magdelaine-Beuzelin C, et al. Infliximab",
    "pharmacokinetics in inflammatory bowel disease patients.",
    "Ther Drug Monit. 2008;30(4):523-529.",
    "doi:10.1097/FTD.0b013e318180e300",
    "(reference 19 of the Kim 2025 main-text reference list; the same model is",
    "labelled [16] in the separately numbered Kim 2025 Tables S1 and S2).",
    "Kim 2025 reproduces the model verbatim as a NONMEM control stream in its",
    "Data S1 supplement (Supplementary Material 1, model #3) with MAXEVAL=0,",
    "i.e. every THETA / OMEGA / SIGMA is held at the originally published value.",
    "Values here are transcribed from that control stream and cross-checked",
    "against Kim 2025 Table S2.",
    sep = " "
  )
  vignette <- "Kim_2025_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central     = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    ADA_POS = list(
      description        = "Anti-drug antibody positivity (antibodies toward infliximab)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "no reference level: a separate typical clearance is estimated in each ADA stratum",
      notes              = "This is the only model in the Kim 2025 panel that estimates a wholly separate typical clearance -- and a separate clearance eta -- in each ADA stratum, rather than applying a multiplicative offset to a single typical value: CL = 0.288 L/day when ADA-negative and 0.768 L/day when ADA-positive, a 2.67-fold increase. Source paper labels this covariate 'ATI' (antibodies toward infliximab); renamed to canonical ADA_POS per covariate-columns.md. Time-varying in the Kim 2025 validation dataset. Kim 2025 showed (Figure 3) that assuming every patient to be ATI-positive improved predictive performance for the five ATI-carrying models, which they attributed to a positive bias inherent in the population PK models.",
      source_name        = "ATI"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "no reference level: a separate sex-specific intercept of the central volume is estimated for each sex",
      notes              = "Contributes an additive, weight-independent term to the central volume: 1.1 L in females and 2.3 L in males, each carrying its own eta. The source control stream codes SEX = 1 for female and SEX = 0 for male, which matches the canonical SEXF polarity directly (no value inversion needed). Time-invariant.",
      source_name        = "SEX"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the weight-proportional component of the central volume, normalized to a 67 kg reference: Vc gains 1.7 * (WT/67) L. Note this is a LINEAR proportionality, not the power function used by every other model in the Kim 2025 panel, and it is added to -- not multiplied by -- the sex-specific intercept. Clearance carries no weight effect in this model. Kim 2025 treated body weight as time-invariant (recorded at the last infliximab concentration measurement).",
      source_name        = "WGT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    age_range      = "Adults; median age 33 years (Kim 2025 Table S1).",
    weight_range   = "Not reported; median body weight 67 kg (Kim 2025 Table S1), which is the reference weight used by the model.",
    sex_female_pct = 55,
    race_ethnicity = "Not specified; developed in a French IBD cohort.",
    disease_state  = "Inflammatory bowel disease: Crohn's disease (n = 30) and ulcerative colitis (n = 3).",
    dose_range     = "Intravenous infliximab during both induction and maintenance phases.",
    regions        = "France.",
    notes          = paste(
      "Development-population characteristics are as summarised by Kim 2025",
      "Table S1 for the Ternant model: Crohn's disease (n = 30) and ulcerative",
      "colitis (n = 3), adults, induction and maintenance phases, peak,",
      "intermediate and trough sampling, 15% ATI positive, median age 33 years,",
      "55% female, median body weight 67 kg. Albumin, ESR and immunomodulator",
      "use were not specified.",
      "This is the smallest development cohort in the Kim 2025 panel (n = 33),",
      "which is worth keeping in view alongside its seven random effects.",
      "The external-validation cohort in Kim 2025 itself (142 Korean adult IBD",
      "patients with 294 concentrations) was used to assess predictive",
      "performance, not to develop or update these parameters; in that cohort",
      "this model gave MPE 37.0%, MAE 0.8 mg/L and rRMSE 294.5%",
      "(Kim 2025 Table 2) -- the best mean absolute error but the worst relative",
      "root-mean-square error of the three adult models.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Clearance is estimated separately in each ADA
    # stratum rather than as a reference value plus an offset, so both strata
    # carry an explicit suffix (no bare lcl).
    lcl_ada_neg <- log(0.288); label("Typical clearance in ADA-negative patients (CL, L/day)")  # Data S1 model #3 $THETA1 = 0.288 L/day; Table S2 "ATI negative: 0.288"
    lcl_ada_pos <- log(0.768); label("Typical clearance in ADA-positive patients (CL, L/day)")  # Data S1 model #3 $THETA2 = 0.768 L/day; Table S2 "ATI positive: 0.768"

    # The central volume is an ADDITIVE decomposition, not the usual
    # multiplicative covariate model:
    #   Vc = lvc_wt * (WT/67) + [lvc_female (female) or lvc_male (male)]
    # so each of the three terms is a volume in litres in its own right and each
    # carries its own eta. The sex terms are weight-independent intercepts; the
    # weight term is shared by both sexes.
    lvc_female <- log(1.1); label("Sex-specific intercept of the central volume, females (L)")  # Data S1 model #3 $THETA3 = 1.1 L; Table S2 "F: 1.1"
    lvc_male   <- log(2.3); label("Sex-specific intercept of the central volume, males (L)")    # Data S1 model #3 $THETA4 = 2.3 L; Table S2 "M: 2.3"
    lvc_wt     <- log(1.7); label("Weight-proportional component of the central volume at 67 kg (L)")  # Data S1 model #3 $THETA5 = 1.7 L
    lvp        <- log(1.9);   label("Peripheral volume of distribution (Vp, L)")        # Data S1 model #3 $THETA6 = 1.9 L; Table S2 "1.9"
    lq         <- log(0.130); label("Inter-compartmental clearance (Q, L/day)")         # Data S1 model #3 $THETA7 = 0.130 L/day; Table S2 "0.130"

    # Inter-individual variability. Following the panel-wide convention stated
    # in the Aubourg stream ("OMEGA = (%CV/100)**2"), each variance is the
    # squared coefficient of variation; Table S2 reports the matching %CV for
    # six of the seven etas. Every structural term in this model carries its
    # own eta, including each ADA stratum of CL and each additive component of
    # Vc -- seven random effects on a development cohort of 33 patients.
    etalcl_ada_neg ~ 0.051  # Data S1 model #3 $OMEGA ETA1 = 0.051; sqrt = 22.6%, matching Table S2 "0.288 (22.6%)"
    etalcl_ada_pos ~ 0.075  # Data S1 model #3 $OMEGA ETA2 = 0.075; sqrt = 27.4%, matching Table S2 "0.768 (27.4%)"
    etalvc_female  ~ 0.013  # Data S1 model #3 $OMEGA ETA3 = 0.013; sqrt = 11.4%, matching Table S2 "F: 1.1 (11.4%)"
    etalvc_male    ~ 0.020  # Data S1 model #3 $OMEGA ETA4 = 0.020; sqrt = 14.1%, matching Table S2 "M: 2.3 (14.1%)"
    etalvc_wt      ~ 0.037  # Data S1 model #3 $OMEGA ETA5 = 0.037 (sqrt = 19.2%); the only eta of the seven that Table S2 does not tabulate, because Table S2 reports a single Vc %CV per sex
    etalvp         ~ 0.023  # Data S1 model #3 $OMEGA ETA6 = 0.023; sqrt = 15.2%, matching Table S2 "1.9 (15.2%)"
    etalq          ~ 0.010  # Data S1 model #3 $OMEGA ETA7 = 0.010; sqrt = 10.0%, matching Table S2 "0.130 (10%)"

    # Residual error -- proportional only. This is the only model in the panel
    # whose $ERROR block carries no additive term (Y = IPRED*(1+ERR(1))), and
    # Table S2 lists its error model as "Proportional".
    propSd <- 0.04; label("Proportional residual error (fraction)")  # Data S1 model #3 $SIGMA = 0.0016 (variance); sqrt(0.0016) = 0.04
  })
  model({
    # Individual PK parameters. ADA_POS selects both the typical clearance and
    # the clearance eta; SEXF selects the volume intercept and its eta.
    #   CL (L/day) = 0.288 (ADA-negative) or 0.768 (ADA-positive)
    #   Vc (L)     = 1.7 * (WT/67) + [1.1 (female) or 2.3 (male)]
    #   Vp (L)     = 1.9
    #   Q  (L/day) = 0.130
    cl <- exp(lcl_ada_neg * (1 - ADA_POS) + lcl_ada_pos * ADA_POS +
                etalcl_ada_neg * (1 - ADA_POS) + etalcl_ada_pos * ADA_POS)

    # Additive volume decomposition: a weight-proportional term plus a
    # sex-specific intercept. Unlike every other model in this panel the weight
    # term is linear in WT/67 rather than a power function.
    vc <- exp(lvc_wt + etalvc_wt) * (WT / 67) +
      exp(lvc_female + etalvc_female) * SEXF +
      exp(lvc_male + etalvc_male) * (1 - SEXF)

    vp <- exp(lvp + etalvp)
    q  <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment model; infliximab is given as an IV infusion into central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L, numerically identical to ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
