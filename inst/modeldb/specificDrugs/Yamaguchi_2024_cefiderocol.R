Yamaguchi_2024_cefiderocol <- function() {
  description <- "Two-compartment pediatric population PK model for intravenous cefiderocol from neonates to adults, with body-weight allometric scaling (beta-lactam class exponents 0.578 on CL and Q, 0.677 on V1 and V2) and a sigmoidal postmenstrual-age renal-maturation function on CL (TM50 47.7 weeks, Hill 3.4); the simulation model of a stochastic simulation and estimation study of pediatric sample size (Yamaguchi 2024)"
  reference <- paste(
    "Yamaguchi D, Katsube T, Wajima T.",
    "Evaluation and optimization of sample size of neonates and infants for",
    "pediatric clinical studies on cefiderocol using a model-based approach.",
    "Pharmacol Res Perspect. 2024;12(5):e70001. doi:10.1002/prp2.70001.",
    "The pediatric model is that of Katsube T, Echols R, Wajima T.",
    "Prediction of cefiderocol pharmacokinetics and probability of target",
    "attainment in pediatric subjects for proposing dose regimens.",
    "IDWeek 2019 abstract, Washington, DC; allometric exponents from",
    "Shimamura K, Wajima T, Yano Y. J Pharm Sci. 2007;96(11):3125-3139.",
    "doi:10.1002/jps.20977; maturation constants from Rhodin MM et al.",
    "Pediatr Nephrol. 2009;24(1):67-76. doi:10.1007/s00467-008-0997-5.",
    sep = " "
  )
  vignette <- "Yamaguchi_2024_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Enters UNNORMALISED as WT^0.578 on CL and Q and",
        "WT^0.677 on V1 and V2 (Yamaguchi 2024 Equations 1-4), so the",
        "structural coefficients are the values at 1 kg body weight."
      ),
      source_name = "BW"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "WEEKS, not the register default of months, because the Rhodin 2009",
        "maturation constants (TM50 = 47.7 weeks, Hill = 3.4) are on that",
        "scale (same convention as Germovsek_2018_meropenem.R and",
        "Riccobene_2017_ceftaroline.R). Required for every subject: for",
        "children past the neonatal period and for adults supply GA (40 weeks",
        "for term births) plus postnatal age in weeks; the maturation factor",
        "is then effectively 1."
      ),
      source_name = "PMA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 516,
    n_studies = 5,
    age_range = "Neonates (gestational age 26 weeks or more, postnatal age from 0 days) to adults (18 to <93 years)",
    weight_range = "Not reported; body weight was simulated from age (Sumpter and Holford 2011; Sy 2014)",
    sex_female_pct = 50,
    race_ethnicity = "Not reported",
    disease_state = "Pediatric patients with Gram-negative infections (planned pediatric trials NCT04335539, NCT04215991 and a neonate/infant trial); adult data from healthy subjects and patients",
    dose_range = "Single 3-h IV infusion: 30-60 mg/kg by gestational and postnatal age below 3 months; 60 mg/kg (body weight <34 kg) or 2 g (34 kg or more) from 3 months to <18 years (Yamaguchi 2024 Table S1)",
    regions = "Not reported",
    pma_range = "26 to <57 weeks across the neonate and infant groups (Yamaguchi 2024 Tables 1 and 2)",
    notes = paste(
      "No pediatric subject was observed. The model is the pediatric",
      "projection of the adult cefiderocol popPK (Kawaguchi 2021) built",
      "by allometric scaling with beta-lactam class exponents and a",
      "renal-maturation factor, and serves as the true model of the",
      "stochastic simulation and estimation study. n_subjects and",
      "n_studies count the 516 adults from five clinical studies whose",
      "observed concentrations formed Group 11 of every design. The",
      "virtual pediatric cohorts were 1:1 female:male with uniformly",
      "distributed gestational and postnatal age per group (Methods",
      "'Simulation'; Tables 1 and 2)."
    )
  )

  ini({
    # Structural coefficients: Yamaguchi 2024 Equations 1-4, and Table 4
    # 'True Value'. Body weight is not normalised, so each is the typical
    # value at 1 kg (and, for CL, at full renal maturation). Not estimated in
    # this paper; these are the simulation model's own typical values.
    lcl <- log(0.508); label("Clearance coefficient at 1 kg body weight, fully mature (L/h)") # Eq. 1 and Table 4: CL = 0.508 L/h
    lvc <- log(0.460); label("Central volume V1 coefficient at 1 kg body weight (L)") # Eq. 2 and Table 4: V1 = 0.460 L
    lq <- log(0.495); label("Intercompartmental clearance coefficient at 1 kg body weight (L/h)") # Eq. 3 and Table 4: Q = 0.495 L/h
    lvp <- log(0.318); label("Peripheral volume V2 coefficient at 1 kg body weight (L)") # Eq. 4 and Table 4: V2 = 0.318 L

    # Allometric exponents taken from Shimamura 2007 (paper ref 21), shared
    # by CL and Q and by V1 and V2.
    e_wt_cl_q <- fixed(0.578); label("Allometric exponent on body weight for CL and Q (unitless)") # Eqs. 1 and 3; Table 4 'Power on CL and Q' = 0.578
    e_wt_vc_vp <- fixed(0.677); label("Allometric exponent on body weight for V1 and V2 (unitless)") # Eqs. 2 and 4; Table 4 'Power on V1 and V2' = 0.677

    # Renal-maturation constants set to the glomerular-filtration-rate values
    # of Rhodin 2009 (paper refs 22-23).
    tmat50 <- fixed(47.7); label("Postmenstrual age at 50 percent renal maturation of CL (weeks)") # Methods 'Pharmacokinetic model': TM50 = 47.7 weeks; Table 4
    hill_mat <- fixed(3.4); label("Hill coefficient of renal maturation of CL (unitless)") # Methods 'Pharmacokinetic model': HILL = 3.4; Table 4

    # Between-subject variability: log-normal, 20% CV on CL, V1 and V2
    # (Methods), printed as variance 0.040 (= 0.20^2) in Table 4.
    etalcl ~ 0.040 # Table 4 'BSV_CL' = 0.040
    etalvc ~ 0.040 # Table 4 'BSV_V1' = 0.040
    etalvp ~ 0.040 # Table 4 'BSV_V2' = 0.040

    # Table 4 prints RUVprop on the same variance scale as the BSV rows
    # (0.040); the proportional SD is sqrt(0.040) = 0.20.
    propSd <- 0.20; label("Proportional residual error (fraction)") # Table 4 'RUVprop' = 0.040 (variance)
  })

  model({
    # Sigmoidal renal maturation on CL, PAGE in weeks (Eq. 1)
    fmat <- PAGE^hill_mat / (PAGE^hill_mat + tmat50^hill_mat)

    cl <- exp(lcl + etalcl) * WT^e_wt_cl_q * fmat
    vc <- exp(lvc + etalvc) * WT^e_wt_vc_vp
    q <- exp(lq) * WT^e_wt_cl_q
    vp <- exp(lvp + etalvp) * WT^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L give mg/L (= ug/mL, the paper's unit)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
