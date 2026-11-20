Gaurav_2026_favezelimab <- function() {
  description <- paste(
    "Two-compartment population PK model of favezelimab (anti-LAG-3 IgG4",
    "monoclonal antibody) with parallel linear and Michaelis-Menten",
    "(target-mediated) elimination from the central compartment, in adults",
    "with advanced solid tumors or PD-L1-positive metastatic colorectal",
    "cancer receiving 7-800 mg IV every 3 weeks alone or with pembrolizumab",
    "(Gaurav 2026). Full-model covariate structure: body weight, albumin,",
    "sex, Asian race, ECOG performance status, hepatic and renal impairment",
    "on CL and Vc, age on CL, and time-varying anti-drug antibody status on",
    "the linear component of CL."
  )
  reference <- paste(
    "Gaurav M, Barcomb H, Maxwell KF, Kandala B, Chatterjee MS.",
    "Population Pharmacokinetic, Exposure-Response Efficacy and Safety",
    "Analyses of Favezelimab in Patients With Solid Tumors.",
    "Clin Transl Sci. 2026;19(7):e70658. doi:10.1111/cts.70658.",
    "PMCID PMC13330136."
  )
  vignette <- "Gaurav_2026_favezelimab"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL (exponent 0.599) and Vc (exponent 0.353) normalised to 71 kg, the reference printed in the Table 1 footnote equations (analysis median 70.9 kg, Table S2). Baseline (time-fixed).",
      source_name = "WTKG"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL (exponent -1.11) and Vc (exponent -0.116) normalised to 40 g/L, the analysis median. Table S2 and the Table 1 footnote label albumin 'g/dL', but the values (median 40, range 21-70.6) and the Results text ('lowest quartile of albumin (<= 37 g/L)') and Figure 1 axis ('Albumin (g/L)') are in g/L, so no conversion is applied. Missing albumin was imputed with the median (Table S2 footnote).",
      source_name = "ALB"
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL only (exponent -0.115) normalised to the analysis median of 59 years (Table S2). Retained in the full model although its 95% CI includes zero; it is omitted from the Table 1 footnote equation, which lists only effects whose CI excludes zero.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male) in this library's encoding; the source estimates the effect for MALE relative to female",
      notes = "The Table 1 footnote equations carry exp(0.2 * [SEX = Male]) on CL and exp(0.26 * [SEX = Male]) on Vc, and the Results typical-value sentence (male, 71 kg: CL 0.023 L/h = 0.0195 * exp(0.200)) confirms the effect applies to males. Encoded as exp(e_male_cl * (1 - SEXF)) so that SEXF = 1 (female) returns the Table 1 typical value. Consistent with the Figure 1 forest plot, where females have higher exposure (GMR 1.39).",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = Caucasian, Black or African American, American Indian or Alaska Native, Native Hawaiian or Other Pacific Islander, multiple or unknown)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (pooled Caucasian / other / missing group)",
      notes = "Exponential effects on CL (-0.245) and Vc (-0.0974). Table S2 footnote: races other than Asian and Caucasian were under 5% of the population and were merged with Caucasian as the reference group.",
      source_name = "RACE"
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 (1 = ECOG 1, 0 = ECOG 0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG performance status 0)",
      notes = "Exponential effects on CL (0.0742) and Vc (-0.0031). Only ECOG 0 (48.1%) and 1 (51.9%) occur in the analysis set (Table S2), so the source ECOGBLN indicator equals ECOG_GE1.",
      source_name = "ECOGBLN"
    ),
    HEPIMP = list(
      description = "Hepatic impairment indicator (1 = mild, moderate or severe; 0 = normal or missing)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal or missing hepatic function)",
      notes = "Exponential effects on CL (0.0461) and Vc (0.0123). Table S2 footnote: mild, moderate and severe were combined because the moderate and severe groups were too small, and missing was pooled with normal. The grading system is not stated in the paper (NCI ODWG is the usual oncology convention).",
      source_name = "hepatic impairment category"
    ),
    RENALIMP = list(
      description = "Renal impairment indicator (1 = mild, moderate or severe; 0 = normal or missing)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal or missing renal function)",
      notes = "Exponential effects on CL (-0.0926) and Vc (-0.0184). The Table 1 footnote equation names the source column RFCAT = mild/moderate/severe. Table S2 footnote: the impairment strata were pooled and missing was grouped with normal. The renal-function thresholds are not stated in the paper.",
      source_name = "RFCAT"
    ),
    ADA_POS = list(
      description = "Anti-drug antibody positive status (time-varying; 1 = positive, 0 = negative)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA-negative)",
      notes = "Time-varying regressor on the LINEAR component of CL only; the source carried the last observed ADA status forward. Multiplicative factor 1.11 on linear CL for ADA-positive samples (Table 1 'ADA Effect on CL' 1.11). The paper does not print the functional form; the factor reading (CL * 1.11^ADA_POS) is the maintainers' interpretation, chosen because the exponential reading exp(1.11) = 3.03 contradicts the steady-state AUC GMR of 0.807 for ADA-positive patients in Figure 1 (see the vignette).",
      source_name = "ADA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "favezelimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "favezelimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 724L,
    n_studies = 2L,
    n_observations = "8980 favezelimab serum concentrations",
    age_range = "25-85 years (median 59)",
    weight_range = "33.0-145.9 kg (median 70.9)",
    sex_female_pct = 38.1,
    race_ethnicity = "Caucasian 68.2%, Asian 27.1%, Black or African American 2.5%, other or unknown 2.1%",
    disease_state = "Advanced solid tumors (P001: colorectal 60.4%, gastric 17.5%, head and neck, other solid tumors) and PD-L1-positive previously treated metastatic colorectal adenocarcinoma (KEYFORM-007)",
    dose_range = "Favezelimab 7, 21, 70, 200, 210, 700 and 800 mg IV every 3 weeks as monotherapy, with pembrolizumab 200 mg Q3W (separate or co-formulated), with or without chemotherapy or lenvatinib",
    regions = "Multinational (MK-4280-001 / NCT02720068 and KEYFORM-007 / NCT05064059)",
    albumin = "median 40 g/L (range 21-70.6)",
    ecog = "ECOG 0 48.1%, ECOG 1 51.9%",
    hepatic_function = "normal 68.5%, mild 29.3%, moderate 0.4%, severe 0.1%, missing 1.7%",
    renal_function = "normal 60.2%, mild 34.0%, moderate 5.5%, severe 0.1%, missing 0.1%",
    ada_incidence = "28.9% ADA-positive overall; 26.2% at 800 mg Q3W; median onset about 3 weeks (Table S3)",
    notes = "Demographics from Gaurav 2026 Results 'Patient Characteristics' and Supplementary Table S2. P001 contributed 475 patients and KEYFORM-007 249."
  )

  ini({
    # ---- Structural PK (Gaurav 2026 Table 1, Monolix SAEM) ------------------
    # Typical values for the reference patient: female, non-Asian, ECOG 0,
    # normal/missing hepatic and renal function, 71 kg, albumin 40 g/L,
    # age 59 years, ADA-negative.
    lcl <- log(0.0195); label("Linear clearance (L/h)") # Table 1 CL (L/h) 0.0195 (%RSE 3.93)
    lvc <- log(3.08); label("Central volume of distribution (L)") # Table 1 Vc (L) 3.08 (%RSE 2.16)
    lq <- log(0.0322); label("Intercompartmental clearance (L/h)") # Table 1 Q (L/h) 0.0322 (%RSE 7.23)
    lvp <- log(1.22); label("Peripheral volume of distribution (L)") # Table 1 Vp (L) 1.22 (%RSE 3.15)
    lvmax <- log(13.6); label("Maximum Michaelis-Menten elimination rate (ug/h)") # Table 1 Vmax (ug/h) 13.6 (%RSE 12.1)
    lkm <- log(49.5); label("Michaelis-Menten constant (ug/L, equal to ng/mL)") # Table 1 Km (ug/L) 49.5 (%RSE 65.4)

    # ---- Covariate effects on CL -------------------------------------------
    # Continuous covariates: power model (cov / median)^theta.
    # Categorical covariates: exponential model exp(theta * indicator).
    e_wt_cl <- 0.599; label("Power exponent of body weight on CL (unitless)") # Table 1 'Body weight effect on CL' 0.599
    e_alb_cl <- -1.11; label("Power exponent of albumin on CL (unitless)") # Table 1 'Albumin effect on CL' -1.11
    e_age_cl <- -0.115; label("Power exponent of age on CL (unitless)") # Table 1 'Age effect on CL' -0.115
    e_male_cl <- 0.200; label("Exponential effect of male sex on CL (unitless)") # Table 1 'Sex effect on CL' 0.200
    e_race_asian_cl <- -0.245; label("Exponential effect of Asian race on CL (unitless)") # Table 1 'Asian race effect on CL' -0.245
    e_ecog_ge1_cl <- 0.0742; label("Exponential effect of ECOG >= 1 on CL (unitless)") # Table 1 'ECOG Effect on CL' 0.0742
    e_hepimp_cl <- 0.0461; label("Exponential effect of hepatic impairment on CL (unitless)") # Table 1 'Hepatic impairment effect on CL' 0.0461
    e_renalimp_cl <- -0.0926; label("Exponential effect of renal impairment on CL (unitless)") # Table 1 'Renal impairment effect on CL' -0.0926
    e_ada_pos_cl <- 1.11; label("Multiplicative factor on linear CL when ADA-positive (unitless)") # Table 1 'ADA Effect on CL' 1.11

    # ---- Covariate effects on Vc -------------------------------------------
    e_wt_vc <- 0.353; label("Power exponent of body weight on Vc (unitless)") # Table 1 'Body weight effect on Vc' 0.353
    e_alb_vc <- -0.116; label("Power exponent of albumin on Vc (unitless)") # Table 1 'Albumin effect on Vc' -0.116
    e_male_vc <- 0.255; label("Exponential effect of male sex on Vc (unitless)") # Table 1 'Sex effect on Vc' 0.255
    e_race_asian_vc <- -0.0974; label("Exponential effect of Asian race on Vc (unitless)") # Table 1 'Asian race effect on Vc' -0.0974
    e_ecog_ge1_vc <- -0.0031; label("Exponential effect of ECOG >= 1 on Vc (unitless)") # Table 1 'ECOG effect on Vc' -0.0031
    e_hepimp_vc <- 0.0123; label("Exponential effect of hepatic impairment on Vc (unitless)") # Table 1 'Hepatic impairment effect on Vc' 0.0123
    e_renalimp_vc <- -0.0184; label("Exponential effect of renal impairment on Vc (unitless)") # Table 1 'Renal impairment effect on Vc' -0.0184

    # ---- Between-subject variability ---------------------------------------
    # Table 1 reports SDs omega_CL = 0.409 and omega_Vc = 0.184 and the
    # correlation corr_Vc_CL = 0.242. Variances: 0.409^2 = 0.167281,
    # 0.184^2 = 0.033856; covariance 0.242 * 0.409 * 0.184 = 0.01821195.
    etalcl + etalvc ~ c(0.167281, 0.01821195, 0.033856) # Table 1 omega_CL 0.409, omega_Vc 0.184, corr_Vc_CL 0.242

    # ---- Residual error ------------------------------------------------------
    addSd <- 30.3; label("Additive residual error (ng/mL)") # Table 1 'Additive Error' 30.3 (ng/mL per Results text)
    propSd <- 0.322; label("Proportional residual error (fraction)") # Table 1 'Proportional Error' 0.322
  })

  model({
    # ---- Individual parameters (Table 1 footnote equations) ----------------
    # The source reference category for sex is FEMALE (the effect is
    # estimated for males), so male = 1 - SEXF.
    cl <- exp(lcl + etalcl) *
      (WT / 71)^e_wt_cl *
      (ALB / 40)^e_alb_cl *
      (AGE / 59)^e_age_cl *
      exp(e_male_cl * (1 - SEXF)) *
      exp(e_race_asian_cl * RACE_ASIAN) *
      exp(e_ecog_ge1_cl * ECOG_GE1) *
      exp(e_hepimp_cl * HEPIMP) *
      exp(e_renalimp_cl * RENALIMP) *
      e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) *
      (WT / 71)^e_wt_vc *
      (ALB / 40)^e_alb_vc *
      exp(e_male_vc * (1 - SEXF)) *
      exp(e_race_asian_vc * RACE_ASIAN) *
      exp(e_ecog_ge1_vc * ECOG_GE1) *
      exp(e_hepimp_vc * HEPIMP) *
      exp(e_renalimp_vc * RENALIMP)
    q <- exp(lq)
    vp <- exp(lvp)
    vmax <- exp(lvmax)
    km <- exp(lkm)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Amounts in mg and volumes in L give mg/L = ug/mL; x 1000 gives ng/mL,
    # the scale of the additive residual error and of Km (ug/L = ng/mL).
    Cc <- 1000 * central / vc

    # Michaelis-Menten arm: vmax (ug/h) * Cc / (km + Cc) is in ug/h; divide
    # by 1000 for mg/h to match the mg amount in central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      (vmax / 1000) * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
