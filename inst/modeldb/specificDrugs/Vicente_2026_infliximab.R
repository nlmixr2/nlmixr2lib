Vicente_2026_infliximab <- function() {
  description <- "One-compartment population PK model with first-order absorption for subcutaneous infliximab (CT-P13 SC) in adults with inflammatory bowel disease in a Spanish real-world cohort of anti-TNF-naive and IV-to-SC switch patients (Vicente 2026). The SC absorption rate constant is fixed to the pivotal Hanzel 2021 value; SC bioavailability is estimated on the logit scale and applies to the depot only, so intravenous doses of the switch patients enter the central compartment unscaled. Clearance increases linearly with total body weight centred on 73 kg and by 70 percent for observations coded as suspected immunogenicity (a time-varying, observation-level composite of confirmed anti-drug antibodies or an otherwise unexplained trough drop of more than 33 percent). Inter-individual variability on clearance only (IIV on V and F fixed to zero); proportional residual error."
  reference <- paste(
    "Vicente B, Zazo H, Sanchez-Hernandez JG, Revilla N, Teixeira-da-Silva P.",
    "Population Pharmacokinetics of Subcutaneous Infliximab: A Real-World Study",
    "of Naive and Switch Patients. Pharmaceutics. 2026;18(9):1162.",
    "doi:10.3390/pharmaceutics18091162.",
    sep = " "
  )
  vignette <- "Vicente_2026_infliximab"
  dosing <- c("depot", "central")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds. Vicente 2026 Methods 2.3: "a one-compartment
  # model with first-order absorption and elimination was implemented (ADVAN2,
  # TRANS2)"; Methods 2.2: serum IFX measured by the Promonitor ELISA.
  compartmentData <- list(
    depot = list(analyte = "infliximab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL centred on the population median of 73 kg: CL = 0.017 * (1 + 0.00888 * (TBW - 73)) (Vicente 2026 Eq. 1 and Results 3.2.1, 'TBW is centered on the population median of 73 kg'). Table 1 reports the cohort median as 74 kg (IQR 65-82); the equation's 73 kg centring value is used as printed. Table 1 footnote: non-longitudinal variables are a single measurement per patient, so weight was effectively time-fixed in the analysis.",
      source_name = "TBW"
    ),
    ADA_SUSPECTED = list(
      description = "Suspected immunogenicity indicator, observation-level and time-varying: 1 = the concentration record meets either (i) analytically confirmed anti-drug antibodies or (ii) a trough decrease of more than 33 percent relative to the previous stable level not explained by non-adherence, recent dose change, delayed administration or an acute flare; 0 otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no suspected immunogenicity)",
      notes = "Multiplicative proportional effect on CL: (1 + 0.700 * ADA_SUSPECTED), i.e. +70.0% CL (Vicente 2026 Eq. 1 and Table 3 'CLSUSP' = 0.700). Definition from Vicente 2026 Methods 2.2. The paper stresses that SUSP is an exploratory PK-informed surrogate rather than confirmed immunogenicity, and that because a concentration drop is part of its definition the effect is partly circular with the PK outcome. It is a superset of assay-confirmed ADA positivity (ADA_POS): only 4 of 138 assayed samples were ADA-positive while 32 of 250 samples were SUSP-positive (Results 3.1).",
      source_name = "SUSP"
    )
  )

  # Covariates screened by stepwise covariate modelling (Methods 2.3) but not
  # retained in the final model (Discussion: none of these 'exert a
  # statistically significant or clinically relevant effect on the PK
  # parameters in this cohort'). The paper reports no estimate for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the SCM; not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the SCM; not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened in the SCM ('HGT', reported in m in Table 1); not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the SCM; not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    SCORE_CALPRO = list(
      description = "Fecal calprotectin",
      units = "ug/g",
      type = "continuous",
      notes = "Screened in the SCM ('FCP'); not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened in the SCM ('CREA'); not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened in the SCM (reported in g/dL in Table 1; canonical unit g/L); not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    FERRITIN_BL = list(
      description = "Serum ferritin",
      units = "ug/L",
      type = "continuous",
      notes = "Screened in the SCM ('FER', measured longitudinally, ng/mL); not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened in the SCM ('HB'); not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened in the SCM; not retained (Vicente 2026 Methods 2.3, Discussion)."
    ),
    ADA_POS = list(
      description = "Assay-confirmed anti-drug antibody positivity",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the SCM ('ADA'); not retained. Its information is carried by the composite ADA_SUSPECTED covariate (Vicente 2026 Methods 2.3, Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    n_observations = 250L,
    age_range = "Adults; median 56 years (IQR 38-61)",
    weight_range = "Median 74 kg (IQR 65-82)",
    sex_female_pct = 31.1,
    race_ethnicity = "Not reported (Spanish multicentre cohort).",
    disease_state = "Inflammatory bowel disease: ulcerative colitis 64.4%, Crohn's disease 35.6%; clinically stable maintenance phase (median fecal calprotectin 69 ug/g, CRP 1 mg/L).",
    dose_range = "SC infliximab 120 mg every 2 weeks (some patients every 1 week, 10 days or 3 weeks under TDM); prior IV infliximab 5-10 mg/kg every 4-8 weeks in switch patients.",
    regions = "Spain (University Hospitals of Burgos and Salamanca), January 2020 - June 2024.",
    treatment_status = "16 (35.5%) anti-TNF-naive patients starting SC after IV induction; 29 (64.5%) switched from IV maintenance.",
    co_medication = "Concomitant azathioprine or methotrexate in 30 (66.7%); prior biologic in 10 (22.2%).",
    notes = "Vicente 2026 Table 1. Retrospective real-world data, almost exclusively trough (pre-dose) samples; median 6 (IQR 4-7) samples per patient. ADA-positive in 3 patients (4/138 assayed samples); suspected immunogenicity in 12 patients (32/250 samples)."
  )

  ini({
    # Final-model column of Vicente 2026 Table 3 (NONMEM 7.5.1, FOCEI).
    lka <- fixed(log(0.0113)); label("SC absorption rate constant (1/h)") # Table 3 Ka = 0.0113 h-1, fixed (Hanzel 2021 value; Table S1 sensitivity analysis)
    lcl <- log(0.0170); label("Clearance at 73 kg without suspected immunogenicity (L/h)") # Table 3 CL = 0.0170 L/h
    lvc <- log(9.14); label("Volume of distribution (L)") # Table 3 V = 9.14 L
    logitfdepot <- log(0.750 / (1 - 0.750)); label("Logit of SC bioavailability (unitless)") # Table 3 F_SC = 0.750 (estimated on the logit scale, Methods 2.3)

    e_wt_cl <- 0.00888; label("Linear effect of total body weight on CL per kg from 73 kg (1/kg)") # Table 3 CLTBW = 0.00888; Eq. 1
    e_ada_suspected_cl <- 0.700; label("Fractional increase in CL with suspected immunogenicity (unitless)") # Table 3 CLSUSP = 0.700; Eq. 1

    # IIV: exponential (log-normal). The 25.7% is converted as a CV,
    # omega^2 = log(1 + 0.257^2) = 0.06393.
    etalcl ~ 0.06393 # Table 3 IIV_CL = 25.7%
    etalvc ~ fixed(0) # Table 3 IIV_V = 0 (held at zero for stability)
    etalogitfdepot ~ fixed(0) # Table 3 IIV_F_SC = 0 (held at zero; 28.4% in the basic model only)

    propSd <- 0.391; label("Proportional residual error (fraction)") # Table 3 RUV = 39.1%, proportional (Results 3.2.1)
  })
  model({
    # Vicente 2026 Eq. 1: CL = 0.017 * (1 + 0.70 * SUSP) * (1 + 0.00888 * (TBW - 73))
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (1 + e_ada_suspected_cl * ADA_SUSPECTED) * (1 + e_wt_cl * (WT - 73))
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Bioavailability applies to SC (depot) doses only; IV doses of the switch
    # patients go to central with F = 1.
    f(depot) <- expit(logitfdepot + etalogitfdepot)

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
