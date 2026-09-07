Yang_2025_iruplinalkib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "first-order elimination for oral iruplinalkib (WX-0593, a selective",
    "ALK/ROS1 tyrosine kinase inhibitor approved in China for ALK-positive",
    "non-small-cell lung cancer), pooled over four Chinese trials in 392",
    "subjects: 16 healthy volunteers and 376 patients with solid tumors.",
    "Apparent oral clearance carries power effects of baseline body weight,",
    "time-varying serum albumin, time-varying creatinine clearance and",
    "time-varying lactate dehydrogenase; the apparent central volume carries a",
    "power effect of baseline body weight. Food slows absorption without",
    "changing exposure: dosing in the fed state multiplies the absorption rate",
    "constant by 0.588 and adds a 0.472 h absorption lag. Residual error is",
    "proportional, with a separate magnitude for the healthy-volunteer",
    "food-effect study WX-0593-002 (55.7%) and for the three solid-tumor",
    "studies WX-0593-001, -003 and -004 (35.7%)."
  )

  reference <- paste(
    "Yang G, Wang Y, Zhao H, Jiang Z, Zheng S, Ge M, Si M, Kang X.",
    "Population pharmacokinetics of iruplinalkib in healthy volunteers and",
    "patients with solid tumors.",
    "Clin Transl Sci. 2025;18(1):e70099.",
    "doi:10.1111/cts.70099"
  )

  vignette <- "Yang_2025_iruplinalkib"

  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "iruplinalkib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "iruplinalkib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "iruplinalkib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT_BASE = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed BASELINE body weight, which is what Yang 2025 fitted",
        "(Results and Eqs. 2-3 both name the covariate 'BBWT', baseline body",
        "weight, in contrast to the three laboratory covariates on CL/F which",
        "the same paper describes as time-varying). Power effects on both",
        "CL/F (exponent 0.441) and V1/F (exponent 1.36), each normalised to",
        "63 kg. The 63 kg reference is the pooled cohort median (Table 1) and",
        "is also named as the typical-subject value in Results 'Covariate",
        "effects on iruplinalkib steady-state exposure'. Pooled range",
        "35.0-98.9 kg; the paper's Figure 2 forest plot spans the 5th-95th",
        "percentile window 48-85 kg. Body weight was strongly correlated with",
        "BMI, so the two were screened separately and only weight retained",
        "(Results 'PopPK analysis')."
      ),
      source_name        = "BBWT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING serum albumin (Yang 2025 abstract, Results and",
        "Discussion all say 'time-varying albumin'; the Methods screened each",
        "laboratory covariate 'at baseline and at the time of measurements').",
        "Power effect on CL/F with exponent 1.05, normalised to 43.89 g/L -",
        "the typical-subject value named in Results 'Covariate effects on",
        "iruplinalkib steady-state exposure'. Note that 43.89 g/L is NOT the",
        "Table 1 pooled BASELINE median of 40.9 g/L; the reference belongs to",
        "the longitudinal record, which sits higher. Table 1 baseline range",
        "24.3-53.8 g/L; the Figure 2 forest plot spans 34.7-49.6 g/L.",
        "The exponent is POSITIVE, i.e. lower albumin lowers CL/F and so",
        "RAISES exposure - the opposite of what a free-fraction mechanism",
        "would predict. The Discussion flags this explicitly and proposes that",
        "low albumin marks a poorer disease state and impaired hepatic",
        "metabolism rather than acting through protein binding. Iruplinalkib",
        "is 76% protein bound."
      ),
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Creatinine clearance, NOT BSA-normalised",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw (absolute) creatinine clearance in mL/min, not normalised to",
        "1.73 m^2 - the same convention as Wada_2023_sparsentan.R,",
        "Delattre_2010_amikacin.R and Chen_2023_nemonoxacin.R. Table 1 of",
        "Yang 2025 lists creatinine clearance (mL/min) and eGFR",
        "(mL/min/1.73 m^2) as two separate rows, which settles the units:",
        "this column is the un-normalised one. The two were strongly",
        "correlated and therefore screened separately, with creatinine",
        "clearance retained (Results 'PopPK analysis'). TIME-VARYING.",
        "Power effect on CL/F with exponent 0.22, normalised to",
        "84.79 mL/min - the typical-subject value named in Results",
        "'Covariate effects on iruplinalkib steady-state exposure', which is",
        "again NOT the Table 1 pooled baseline median of 97.9 mL/min.",
        "Table 1 baseline range 30.8-243 mL/min; the Figure 2 forest plot",
        "spans 48.87-141.7 mL/min. The paper attributes the effect to the",
        "20.23% of a dose recovered in urine. The paper does not name the",
        "estimating equation."
      ),
      source_name        = "CRCL"
    ),
    LDH = list(
      description        = "Serum lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING serum lactate dehydrogenase. Power effect on CL/F with",
        "exponent -0.225, normalised to 242.63 U/L - the typical-subject",
        "value named in Results 'Covariate effects on iruplinalkib",
        "steady-state exposure', not the Table 1 pooled baseline median of",
        "219 U/L. Table 1 baseline range 86.6-1160 U/L; the Figure 2 forest",
        "plot spans 161-388.45 U/L. The exponent is NEGATIVE, so a rising LDH",
        "lowers CL/F and raises exposure. LDH is a general marker of tissue",
        "turnover and tumour burden in solid-tumour cohorts; the Discussion",
        "groups it with creatinine clearance under the partial renal",
        "elimination route rather than giving it an independent mechanism."
      ),
      source_name        = "LDH"
    ),
    FED = list(
      description        = "Fed-versus-fasted dose-record indicator, 1 = fed",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "Two effects on absorption only, both from Yang 2025 Table 2 and",
        "Eqs. 1 and 6: the absorption rate constant is multiplied by 0.588",
        "when fed (Ka 1.06 -> 0.623 /h), and a 0.472 h absorption lag is",
        "added when fed (there is NO lag in the fasted state - Table 2 names",
        "the row 'ALAG for fed subjects' and Eq. 6 writes",
        "'ALAG = 0.472 (if fed)'). Food therefore slows and delays absorption",
        "without touching CL/F or bioavailability, which is exactly the",
        "paper's finding that 'there was no difference in exposure of",
        "iruplinalkib between the fasted and fed states'. The only fed data",
        "come from the two-period crossover study WX-0593-002, in which 16",
        "healthy volunteers received a single 120 mg dose fasted then fed, or",
        "in the reverse order, with a 7-day washout (Table S1); the indicator",
        "is therefore time-varying per dose record within a subject. The",
        "paper does not state the meal composition, so the general FED",
        "canonical is used rather than FED_HIGHFAT."
      ),
      source_name        = "Food"
    ),
    STUDY_WX0593_002 = list(
      description        = "WX-0593-002 study indicator: 1 = the food-effect crossover study in healthy volunteers, 0 = studies WX-0593-001, -003 and -004 in patients with solid tumors",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the pooled solid-tumor studies WX-0593-001, -003 and -004)",
      notes              = paste(
        "Switches the proportional residual-error magnitude only; it carries",
        "no structural or covariate effect. Yang 2025 Table 2 reports",
        "'Residual proportional errors for WX-0593-001, 003 and 004' = 35.7%",
        "against 'Residual proportional errors for WX-0593-002' = 55.7%, so",
        "the healthy-volunteer crossover is the noisier of the two despite",
        "its far denser sampling schedule (21 timepoints per period against",
        "sparse trough sampling in the phase 2/3 studies - Table S1).",
        "WX-0593-002 contributed 16 of the 392 subjects. Same construction as",
        "STUDY_NIPOCALIMAB_PHASE1 in Valenzuela_2025_nipocalimab.R."
      ),
      source_name        = "Study"
    )
  )

  # Covariates that Yang 2025 screened (Methods 'Covariate model development')
  # but did not retain in the final model. Documented here so the provenance of
  # the covariate screen is preserved without carrying unused covariateData
  # entries.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened; strongly correlated with baseline body weight, so the two",
        "were entered separately and body weight - the greater effect - was",
        "retained (Results 'PopPK analysis'). Pooled baseline median",
        "23.6 kg/m^2, range 16.5-35.2 (Table 1)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened at baseline and time-varying; strongly correlated with AST",
        "and not retained. Pooled baseline median 19.4 U/L, range 2.00-172",
        "(Table 1)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened at baseline and time-varying; strongly correlated with ALT",
        "and not retained. Pooled baseline median 21.0 U/L, range 8.00-178",
        "(Table 1)."
      )
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened at baseline and time-varying; not retained. Pooled baseline median 10.2 umol/L, range 3.20-32.6 (Table 1)."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and not retained; also compared",
        "post hoc as >=65 vs <65 years with no apparent effect on exposure",
        "(Results 'Effect of other demographic factors'; Figure S3). Pooled",
        "median 52.0 years, range 25.0-76.0; 14.5% were 65 or older",
        "(Table 1)."
      )
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained. Cohort 52.0% female (Table 1)."
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator, 1 = healthy subject, 0 = patient with a solid tumor",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as the 'subject type' covariate and not retained on any",
        "structural parameter - only the residual-error magnitude differs",
        "between the healthy-volunteer study and the patient studies, and",
        "that is carried by STUDY_WX0593_002. 16 of 392 subjects (4.1%) were",
        "healthy volunteers (Table 1)."
      )
    ),
    WHO_PS = list(
      description = "Eastern Cooperative Oncology Group performance status",
      units       = "(score 0-5)",
      type        = "continuous",
      notes       = paste(
        "Screened and not retained. Pooled distribution 0: 28.1%, 1: 69.4%,",
        "2: 1.8%, missing 0.8% (Table 1); no subject scored above 2."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 392L,
    n_studies      = 4L,
    n_observations = 3788L,
    age_range      = "25.0-76.0 years",
    age_median     = "52.0 years",
    weight_range   = "35.0-98.9 kg",
    weight_median  = "63.0 kg",
    sex_female_pct = 52.0,
    race_ethnicity = c(Han = 95.7, Other = 4.3),
    disease_state  = paste(
      "Pooled: 376 Chinese patients with ALK/ROS1-positive advanced solid",
      "tumors (predominantly non-small-cell lung cancer) across a phase 1",
      "dose-escalation/expansion trial, a single-arm phase 2 trial and a",
      "randomised phase 3 trial, plus 16 Chinese healthy volunteers in a",
      "two-period food-effect crossover"
    ),
    renal_function = paste(
      "Baseline creatinine clearance median 97.9 mL/min, range",
      "30.8-243 mL/min; eGFR median 98.7 mL/min/1.73 m^2, range 46.9-234.",
      "Only normal and mild-to-moderate impairment were represented; the",
      "paper states the effect of moderate or severe impairment remains to",
      "be defined"
    ),
    hepatic_function = paste(
      "Baseline albumin median 40.9 g/L (range 24.3-53.8), ALT median",
      "19.4 U/L, AST median 21.0 U/L, total bilirubin median 10.2 umol/L.",
      "Only normal and mild impairment (NCI-ODWG) were meaningfully",
      "represented"
    ),
    co_medication  = paste(
      "Prior ALK inhibitor: never 54.1%, crizotinib only 44.6%, other 1.3%.",
      "Crizotinib resistance: no 42.9%, yes 39.3%, missing 17.9%"
    ),
    dose_range     = paste(
      "WX-0593-001: single and once-daily oral doses of 30, 60, 90, 120,",
      "180, 240 and 300 mg. WX-0593-002: a single 120 mg dose fasted and",
      "fed. WX-0593-003 and -004: 180 mg once daily after a 7-day 60 mg",
      "lead-in (Table S1)"
    ),
    regions        = "China",
    notes          = paste(
      "Baseline demographics are Yang 2025 Table 1; the study list, dosing",
      "regimens and PK sampling schedules are Table S1. Bioanalysis was",
      "validated LC-MS/MS with a 2-800 ng/mL calibration range and",
      "WX-0593-d6 as internal standard. Estimation was FOCEI in NONMEM 7.5",
      "with PsN 4.8.1. The final model had a condition number of 27.75 and",
      "was confirmed by 1000 bootstrap replicates (Table 2).",
      "Four further covariates were screened (Methods 'Covariate model",
      "development') and not retained but are not listed in",
      "covariatesDataExcluded because the register has no canonical column",
      "for them, or because the only fitting canonical is already in use:",
      "total bile acid; eGFR (which maps onto the same CRCL canonical as the",
      "retained raw creatinine-clearance column - the two were strongly",
      "correlated, screened separately, and creatinine clearance won, and",
      "eGFR was additionally used post hoc to band renal function as normal /",
      "mild / moderate for Figure S3, with no apparent effect on exposure);",
      "ALK mutational status (72.7% positive, 19.9% negative, 7.4% missing);",
      "and ROS1 mutational status (17.1% positive, 25.5% negative, 57.4%",
      "missing - above the paper's 15% threshold, so under the stated Methods",
      "rule those records formed their own category rather than falling to",
      "the reference). ALK and ROS1 status and crizotinib resistance were",
      "also compared post hoc with no apparent effect on exposure",
      "(Figure S3)."
    )
  )

  ini({
    # Structural parameters. Every value is the typical value for the
    # reference subject named in Yang 2025 Results 'Covariate effects on
    # iruplinalkib steady-state exposure': baseline body weight 63 kg,
    # albumin 43.89 g/L, LDH 242.63 U/L and creatinine clearance
    # 84.79 mL/min, dosed FASTED (FED = 0, so no absorption lag and no
    # slowing of Ka). All disposition parameters are apparent (per unit
    # bioavailability); iruplinalkib was given only orally, so F is not
    # separately identifiable and no bioavailability parameter is estimated.

    lka <- log(1.06); label("First-order absorption rate constant Ka, fasted (1/h)")            # Yang 2025 Table 2 final model: Ka = 1.06 /h (RSE 8.4%; bootstrap median 1.05, 95% CI 0.873-1.27); Eq. 1
    lcl <- log(18.9); label("Apparent oral clearance CL/F (L/h)")                               # Yang 2025 Table 2 final model: CL/F = 18.9 L/h (RSE 2.2%; bootstrap median 18.9, 95% CI 18.1-19.7); Eq. 2
    lvc <- log(348);  label("Apparent central volume of distribution V1/F (L)")                 # Yang 2025 Table 2 final model: V1/F = 348 L (RSE 6.7%; bootstrap median 346, 95% CI 304-399); Eq. 3
    lq  <- log(15.5); label("Apparent intercompartmental clearance Q/F (L/h)")                  # Yang 2025 Table 2 final model: Q/F = 15.5 L/h (RSE 11.4%; bootstrap median 15.5, 95% CI 12.0-20.6); Eq. 5
    lvp <- log(295);  label("Apparent peripheral volume of distribution V2/F (L)")              # Yang 2025 Table 2 final model: V2/F = 295 L (RSE 8.4%; bootstrap median 296, 95% CI 252-359); Eq. 4

    # Absorption lag. Yang 2025 Eq. 6 is "ALAG = 0.472 (if fed)" and Table 2
    # names the row "ALAG for fed subjects"; there is no lag in the fasted
    # state. The value below is therefore the FED-state lag, gated by the FED
    # indicator inside model().
    ltlag <- log(0.472); label("Absorption lag time in the fed state (h)")                      # Yang 2025 Table 2 final model: ALAG for fed subjects = 0.472 h (RSE 2.5%; bootstrap median 0.473, 95% CI 0.436-0.880); Eq. 6

    # Food effect on absorption rate. Yang 2025 Eq. 1 writes
    #   Ka = 1.06 * 0.588 (if fed) * exp(eta1)
    # i.e. a multiplicative FACTOR with the fasted state as reference, not a
    # log-scale coefficient. Applied as e_fed_ka^FED so the verbatim 0.588 is
    # preserved and Ka is unchanged when fasted.
    e_fed_ka <- 0.588; label("Multiplicative factor of the fed state on Ka, fasted reference (unitless)")  # Yang 2025 Table 2 final model: Food on Ka = 0.588 (RSE 21.3%; bootstrap median 0.600, 95% CI 0.389-1.02); Eq. 1

    # Covariate effects on CL/F (Yang 2025 Eq. 2). All four are power terms on
    # a covariate normalised to the typical-subject value.
    e_wt_base_cl <-  0.441; label("Power exponent of baseline body weight (/63 kg) on CL/F (unitless)")            # Yang 2025 Table 2 final model: Baseline body weight on CL/F = 0.441 (RSE 31.1%; bootstrap median 0.431, 95% CI 0.160-0.709); Eq. 2 (BBWT/63)^0.441
    e_ldh_cl     <- -0.225; label("Power exponent of lactate dehydrogenase (/242.63 U/L) on CL/F (unitless)")      # Yang 2025 Table 2 final model: LDH on CL/F = -0.225 (RSE 19.2%; bootstrap median -0.228, 95% CI -0.306 to -0.144); Eq. 2 (LDH/242.63)^-0.225
    e_crcl_cl    <-  0.22;  label("Power exponent of creatinine clearance (/84.79 mL/min) on CL/F (unitless)")     # Yang 2025 Table 2 final model: Creatinine clearance on CL/F = 0.22 (RSE 30.1%; bootstrap median 0.22, 95% CI 0.097-0.348); Eq. 2 (CRCL/84.79)^0.22
    e_alb_cl     <-  1.05;  label("Power exponent of serum albumin (/43.89 g/L) on CL/F (unitless)")               # Yang 2025 Table 2 final model: Albumin on CL/F = 1.05 (RSE 12.6%; bootstrap median 1.04, 95% CI 0.778-1.32); Eq. 2 (ALB/43.89)^1.05

    # Covariate effect on V1/F (Yang 2025 Eq. 3).
    e_wt_base_vc <- 1.36; label("Power exponent of baseline body weight (/63 kg) on V1/F (unitless)")              # Yang 2025 Table 2 final model: Baseline body weight on V1/F = 1.36 (RSE 25%; bootstrap median 1.33, 95% CI 0.562-2.00); Eq. 3

    # Inter-individual variability. Yang 2025 Methods: "Individual variability
    # in the PK parameters was estimated using an exponential relationship for
    # all PK parameters", i.e. theta_i = theta_typical * exp(eta_i) with
    # eta ~ N(0, omega^2) - the form written out in Eqs. 1-3. Table 2 reports
    # the variability as "IIV (%)" and its abbreviation list defines
    # "CV, coefficient of variation", so the tabulated percentages are
    # log-normal CVs and are converted here with omega^2 = log(1 + CV^2).
    # No IIV is reported on Q/F, V2/F or the absorption lag (Table 2 shows "-"
    # for those rows) and no eta correlations are reported, so OMEGA is
    # diagonal.
    etalka ~ 0.444368  # Yang 2025 Table 2 final model: IIV Ka   = 74.8% (IIV_RSE 10.9%, shrinkage 58.6%); log(1 + 0.748^2)
    etalcl ~ 0.059687  # Yang 2025 Table 2 final model: IIV CL/F = 24.8% (IIV_RSE 6.1%,  shrinkage 14.8%); log(1 + 0.248^2)
    etalvc ~ 0.250880  # Yang 2025 Table 2 final model: IIV V1/F = 53.4% (IIV_RSE 14.1%, shrinkage 39.4%); log(1 + 0.534^2)

    # Residual error: proportional only, with a separate magnitude per study.
    # Table 2 gives the two rows as percentages; Table S2 heads the same two
    # rows "CV of residual errors ... %", which confirms they are proportional
    # SDs on the linear scale.
    propSdPatient <- 0.357; label("Proportional residual error SD for studies WX-0593-001, -003 and -004 (fraction)")  # Yang 2025 Table 2 final model: 35.7% (RSE 4.6%; bootstrap median 35.4, 95% CI 32.3-39.0)
    propSdHv      <- 0.557; label("Proportional residual error SD for study WX-0593-002 (fraction)")                   # Yang 2025 Table 2 final model: 55.7% (RSE 6.4%; bootstrap median 55.0, 95% CI 47.5-62.3)
  })

  model({
    # ---- 1. Individual PK parameters --------------------------------------
    # Yang 2025 Eq. 1: Ka = 1.06 * 0.588 (if fed) * exp(eta1)
    ka <- exp(lka + etalka) * e_fed_ka^FED

    # Yang 2025 Eq. 2:
    #   CL/F = 18.9 * (BBWT/63)^0.441 * (LDH/242.63)^-0.225 *
    #                 (CRCL/84.79)^0.22 * (ALB/43.89)^1.05 * exp(eta2)
    cl <-
      exp(lcl + etalcl) *
      (WT_BASE / 63)^e_wt_base_cl *
      (LDH / 242.63)^e_ldh_cl *
      (CRCL / 84.79)^e_crcl_cl *
      (ALB / 43.89)^e_alb_cl

    # Yang 2025 Eq. 3: V1/F = 348 * (BBWT/63)^1.36 * exp(eta3)
    vc <- exp(lvc + etalvc) * (WT_BASE / 63)^e_wt_base_vc

    # Yang 2025 Eqs. 4-5: V2/F and Q/F carry no covariates and no IIV.
    q  <- exp(lq)
    vp <- exp(lvp)

    # Yang 2025 Eq. 6: ALAG = 0.472 (if fed). Gating by FED gives a zero lag
    # in the fasted state, which is the model's reference condition.
    tlag <- exp(ltlag) * FED

    # ---- 2. Micro-constants ------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. ODE system -----------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <-
      ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 4. Absorption lag -------------------------------------------------
    alag(depot) <- tlag

    # ---- 5. Observation and error model ------------------------------------
    # Study-specific proportional residual magnitude (Yang 2025 Table 2), the
    # same indicator-weighted switch used in Valenzuela_2025_nipocalimab.R.
    propSd <-
      propSdHv * STUDY_WX0593_002 +
      propSdPatient * (1 - STUDY_WX0593_002)

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
