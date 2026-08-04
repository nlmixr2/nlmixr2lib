Zhang_2015_dolutegravir <- function() {
  description <- "One-compartment oral population PK model with first-order absorption, absorption lag time, and first-order elimination for dolutegravir (integrase strand transfer inhibitor) in HIV-1-infected treatment-naive adults (Zhang 2015)"
  reference <- "Zhang J, Hayes S, Sadler BM, Minto I, Brandt J, Piscitelli S, Min S, Song IH. Population pharmacokinetics of dolutegravir in HIV-infected treatment-naive patients. Br J Clin Pharmacol. 2015;80(3):502-514. doi:10.1111/bcp.12639"
  vignette <- "Zhang_2015_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot   = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F and V/F, reference 70 kg (Zhang 2015 Table 3 footnote). Range in the analysis population 39-135 kg (Zhang 2015 Table 2).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F, reference 40 years (Zhang 2015 Table 3 footnote). Range 18-68 years.",
      source_name        = "AGE"
    ),
    TBILI = list(
      description        = "Baseline total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F, reference 9 umol/L (Zhang 2015 Table 3 footnote). Range 3-38 umol/L. Higher bilirubin decreases CL, plausibly through competition with dolutegravir for UGT1A1 metabolism (Zhang 2015 Discussion).",
      source_name        = "BILI"
    ),
    SEXF = list(
      description        = "Sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on bioavailability F; female subjects have 21% higher oral bioavailability than males (Zhang 2015 Table 3 footnote: F = 1.21^GEND, GEND = 1 for female). Source column is GEND (1 = female, 0 = male) in the Zhang 2015 dataset, renamed to canonical SEXF.",
      source_name        = "GEND"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-current smoker; pools never and former smokers)",
      notes              = "Exponential effect on CL/F; current smokers have 16% higher CL than non-current smokers (Zhang 2015 Table 3 footnote: CL/F factor = 1.16^SMOK). Attributed to tobacco-induced UGT1A1 / CYP1A1 / CYP1A2 / CYP2B6 / CYP3A4 activity (Zhang 2015 Discussion).",
      source_name        = "SMOK"
    ),
    STUDY_ING111521 = list(
      description        = "Study ING111521 (proof-of-concept, phase 2a monotherapy) cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (SPRING-1 phase 2b or SPRING-2 phase 3 combination cohort)",
      notes              = "Exponential effect on CL/F; typical CL/F is 35% higher in study ING111521 than in SPRING-1 / SPRING-2 (Zhang 2015 Table 3 footnote: CL/F factor = 1.35^POC). The paper attributes this unexplained study difference to the smaller sample size and less diverse patient population in the phase 2a study (Zhang 2015 Discussion). Source column POC (= 1 for ING111521, 0 otherwise) renamed to canonical STUDY_ING111521.",
      source_name        = "POC"
    ),
    DOSE_10MG = list(
      description        = "10 mg dose-level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (25 mg or 50 mg dose)",
      notes              = "Exponential effect on bioavailability F; the 10 mg dose has 24% higher relative oral bioavailability than the pooled 25 / 50 mg reference (Zhang 2015 Table 3 footnote: F = 1.24^DOSE, DOSE = 1 for 10 mg). Continuous dose was screened and was not significant; 25 mg vs 50 mg were indistinguishable. Attributed to better tablet dispersion at the lower strength (Zhang 2015 Discussion). Source column DOSE renamed to canonical DOSE_10MG.",
      source_name        = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    RACE = list(
      description = "Race",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened in the full-model covariate step but not significant; race did not influence dolutegravir PK (Zhang 2015 Results)."
    ),
    ETHNIC = list(
      description = "Hispanic / Latino ethnicity",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    ),
    HBV = list(
      description = "Hepatitis B virus co-infection at baseline",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    ),
    HCV = list(
      description = "Hepatitis C virus co-infection at baseline",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    ),
    CDC = list(
      description = "CDC HIV disease-severity classification (A / B / C)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    ),
    ALB = list(
      description = "Serum albumin at baseline",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not significant; the model preferred age over albumin as the CL/F predictor (Zhang 2015 Discussion)."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    ),
    AST = list(
      description = "Aspartate aminotransferase at baseline",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not significant (Zhang 2015 Results)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 563L,
    n_studies      = 3L,
    n_observations = 3357L,
    age_range      = "18-68 years",
    age_median     = "37 years",
    weight_range   = "39.0-135 kg",
    weight_median  = "74.5 kg",
    sex_female_pct = 15,
    race_ethnicity = c(Caucasian = 83, Black = 12, Asian = 1, Other = 4),
    disease_state  = "HIV-1 infected antiretroviral treatment-naive adults",
    dose_range     = "10, 25, or 50 mg dolutegravir once daily orally, alone (ING111521) or in combination with ABC/3TC or TDF/FTC (SPRING-1, SPRING-2)",
    regions        = "Multinational",
    notes          = "Pooled data from three phase 2/3 studies (Zhang 2015 Table 1): ING111521 (phase 2a monotherapy, n = 19), SPRING-1 (phase 2b, n = 141), SPRING-2 (phase 3, n = 403). Baseline demographics per Zhang 2015 Table 2. Race percentages summed across 'Other' cases per the paper's category assignments; 'Other' includes multiracial / unreported.",
    smoking_pct    = "42% current, 13% former, 42% never (3% unknown; Zhang 2015 Table 2)",
    tbili_median   = "9 umol/L (range 3-38)"
  )

  ini({
    # Structural parameters -- reference covariates:
    #   WT = 70 kg, AGE = 40 y, TBILI = 9 umol/L, SEXF = 0 (male),
    #   SMOKE = 0 (non-current smoker), STUDY_ING111521 = 0 (SPRING),
    #   DOSE_10MG = 0 (25 or 50 mg tablet strength).
    # From Zhang 2015 Table 3 final PK model.
    lka     <- log(2.24);       label("Absorption rate constant ka (1/h)")           # Zhang 2015 Table 3: ka = 2.24 h^-1 (95% CI 1.56, 2.92)
    lcl     <- log(0.901);      label("Apparent oral clearance CL/F (L/h)")           # Zhang 2015 Table 3: CL/F = 0.901 L/h (95% CI 0.864, 0.938)
    lvc     <- log(17.4);       label("Apparent central volume V/F (L)")              # Zhang 2015 Table 3: V/F = 17.4 L (95% CI 16.5, 18.3)
    ltlag   <- log(0.263);      label("Absorption lag time tlag (h)")                 # Zhang 2015 Table 3: ALAG = 0.263 h (95% CI 0.0942, 0.432)
    lfdepot <- fixed(log(1));   label("Reference relative bioavailability F (unitless; anchor at the 25/50 mg dose in males)")

    # Continuous covariate power exponents (Zhang 2015 Table 3 footnote).
    #   CL/F: (WT/70)^e_wt_cl * (AGE/40)^e_age_cl * (TBILI/9)^e_tbili_cl
    #   V/F:  (WT/70)^e_wt_vc
    e_wt_cl    <-  0.438; label("Power exponent of body weight on CL/F (unitless)")            # Zhang 2015 Table 3: CL/F ~ WT = 0.438 (95% CI 0.293, 0.583)
    e_wt_vc    <-  0.768; label("Power exponent of body weight on V/F (unitless)")             # Zhang 2015 Table 3: V/F ~ WT = 0.768 (95% CI 0.605, 0.931)
    e_age_cl   <-  0.193; label("Power exponent of age on CL/F (unitless)")                    # Zhang 2015 Table 3: CL ~ AGE = 0.193 (95% CI 0.103, 0.283)
    e_tbili_cl <- -0.211; label("Power exponent of total bilirubin on CL/F (unitless)")        # Zhang 2015 Table 3: CL ~ BILIRUBIN = -0.211 (95% CI -0.269, -0.153)

    # Categorical / binary covariate effects encoded as exp(theta * indicator);
    # theta = log(fold-change) so exp(theta) reproduces the Zhang 2015 fold-change.
    e_smoke_cl         <- log(1.16); label("log-scale multiplier of current smoking on CL/F (unitless; exp(e_smoke_cl) = 1.16)")           # Zhang 2015 Table 3: CL ~ SMOKING = 1.16 (95% CI 1.10, 1.22)
    e_ing111521_cl     <- log(1.35); label("log-scale multiplier of study ING111521 on CL/F (unitless; exp(e_ing111521_cl) = 1.35)")       # Zhang 2015 Table 3: CL/F ~ proof-of-concept = 1.35 (95% CI 1.22, 1.48)
    e_sexf_fdepot      <- log(1.21); label("log-scale multiplier of female sex on F (unitless; exp(e_sexf_fdepot) = 1.21)")                # Zhang 2015 Table 3: F ~ GENDER = 1.21 (95% CI 1.13, 1.29)
    e_dose_10mg_fdepot <- log(1.24); label("log-scale multiplier of the 10 mg dose level on F (unitless; exp(e_dose_10mg_fdepot) = 1.24)") # Zhang 2015 Table 3: F ~ 10 mg = 1.24 (95% CI 1.17, 1.31)

    # IIV -- omega^2 on the log scale per Zhang 2015 Table 3 (paper's own
    # notation uses omega^2_P). CL / V / ka were reported with IIV; the
    # correlation between CL and V (0.375) was removed from the final
    # model per Zhang 2015 Results (11-point OFV increase considered non-
    # material). IOV on CL/F (omega^2 = 0.0296, 17.2% CV) is dropped per
    # nlmixr2lib convention (BOV is dropped whenever a BSV term is
    # reported for the same parameter); a note in the vignette Assumptions
    # and deviations section documents the drop.
    etalcl ~ 0.0551   # Zhang 2015 Table 3: omega^2_CL = 0.0551 (95% CI 0.0451, 0.0651; CV 23.5%)
    etalvc ~ 0.0188   # Zhang 2015 Table 3: omega^2_V  = 0.0188 (95% CI 0.00794, 0.0297; CV 13.7%)
    etalka ~ 0.224    # Zhang 2015 Table 3: omega^2_ka = 0.224 (95% CI 0.0535, 0.395; CV 50.1%)

    # Residual error -- proportional only (additive was not significant
    # per Zhang 2015 Results). Table 3 reports sigma^2_prop = 0.0704 as
    # a variance; propSd = sqrt(sigma^2_prop) = 0.2653, consistent with
    # the CV = 26.5% figure in the Table 3 footnote.
    propSd <- sqrt(0.0704); label("Proportional residual error (fraction; CV 26.5%)")  # Zhang 2015 Table 3: sigma^2_prop = 0.0704 (95% CI 0.0602, 0.0806)
  })

  model({
    # Individual PK parameters. Reference covariates: WT = 70 kg, AGE = 40 y,
    # TBILI = 9 umol/L, SEXF = 0 (male), SMOKE = 0 (non-current), STUDY_ING111521
    # = 0 (SPRING-1 / SPRING-2), DOSE_10MG = 0 (25 or 50 mg dose).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) *
      (WT    / 70)^e_wt_cl *
      (AGE   / 40)^e_age_cl *
      (TBILI /  9)^e_tbili_cl *
      exp(e_smoke_cl     * SMOKE) *
      exp(e_ing111521_cl * STUDY_ING111521)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # Bioavailability -- reference F fixed to 1 (log-scale anchor); female sex
    # and the 10 mg dose are exponential multipliers on that anchor.
    fdepot <- exp(lfdepot) *
      exp(e_sexf_fdepot      * SEXF) *
      exp(e_dose_10mg_fdepot * DOSE_10MG)

    tlag <- exp(ltlag)
    kel  <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Dose is in mg, V in L => central / vc is in mg/L = ug/mL (matches the
    # paper's plasma concentration units in Table 4).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
