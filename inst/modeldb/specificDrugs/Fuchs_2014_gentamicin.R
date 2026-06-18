Fuchs_2014_gentamicin <- function() {
  description <- "Two-compartment IV population PK model for gentamicin in 1449 preterm and term neonates (Fuchs 2014) with fixed allometric body-weight scaling (0.75 on CL/Q, 1 on Vc/Vp), linear centred-on-median effects of gestational age on CL and Vc, postnatal age on CL, and dopamine co-administration on CL."
  reference <- paste(
    "Fuchs A, Guidi M, Giannoni E, Werner D, Buclin T, Widmer N, Csajka C.",
    "Population pharmacokinetic study of gentamicin in a large cohort of premature",
    "and term neonates.",
    "Br J Clin Pharmacol. 2014;78(5):1090-1101.",
    "doi:10.1111/bcp.12444.",
    sep = " "
  )
  vignette <- "Fuchs_2014_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the time of blood sampling.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL and Q (fixed exponent 0.75) and on Vc and Vp",
        "(fixed exponent 1.0). Reference weight is the cohort median 2.170 kg",
        "(Methods 'Model-based pharmacokinetic analysis' and Table 2 footnote).",
        "Source column BW (Table 1)."
      ),
      source_name        = "BW"
    ),
    GA = list(
      description        = "Gestational age at birth.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters CL and Vc as a linear effect centred on",
        "the cohort median 34 weeks: f_ga = 1 + theta * (GA/34 - 1). Range in",
        "the model-building cohort 24-42 weeks (Table 1)."
      ),
      source_name        = "GA"
    ),
    PNA = list(
      description        = "Postnatal age (chronological time since birth).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Fuchs 2014 reports postnatal age in DAYS (cohort median 1 day, range",
        "0-94 days; Table 1). The canonical PNA column is in MONTHS, so the",
        "linear-effect equation 1 + theta * (PNA_days/1 - 1) is reparameterised",
        "inside model() as 1 + theta * (PNA / pna_ref_months - 1) with",
        "pna_ref_months = 1 / 30.4375 (the cohort median 1 day expressed in",
        "months). Both numerator and denominator carry the same units factor",
        "and cancel, leaving the paper's theta_CL_PNA = 0.054 unchanged.",
        "Users should supply PNA in months in the dataset."
      ),
      source_name        = "PNA"
    ),
    CONMED_DOPA = list(
      description        = "Concomitant dopamine administration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant dopamine).",
      notes              = paste(
        "1 = subject is receiving dopamine concomitantly at the PK observation;",
        "0 = no concomitant dopamine. 9.4% of subjects (136/1449) in the",
        "model-building cohort (Table 1). Enters CL as a multiplicative linear",
        "effect: 1 + theta_CL_DOPA * CONMED_DOPA, with theta_CL_DOPA = -0.120",
        "(Table 2). Treated as time-varying per observation in the source paper;",
        "users should supply the appropriate flag per observation record."
      ),
      source_name        = "DOPA"
    )
  )

  covariatesDataExcluded <- list(
    PMA = list(
      description = "Postmenstrual age (= GA + PNA). Screened during covariate model building.",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained in the final model. The univariate effect on CL",
        "(DeltaOF = -658.1) and on Vc (DeltaOF = -214.6) was strong, but the",
        "paper showed that GA + PNA as two distinct covariates on CL outperformed",
        "PMA alone (DeltaOF = -750.8 in favour of GA + PNA); for Vc the",
        "AIC-preferred model was GA alone (Results 'Covariate model'). PMA is",
        "therefore tracked for population description only."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male). Screened.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained. Model-building cohort 57.5% male, 42.5% female",
        "(Table 1). Results 'Covariate model': 'No other covariates showed any",
        "significant effect on gentamicin disposition (DeltaOF > -6.1, P > 0.01)'."
      )
    ),
    CONMED_INDOMETH = list(
      description = "Concomitant indomethacin co-administration indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained. Univariate CL reduction of 18% (DeltaOF = -7.8,",
        "P = 0.005) but the bootstrap 95% CI included 0 and the indomethacin",
        "coefficient was omitted from the final model (Results 'Model validation",
        "and simulation'). 1.9% of subjects (27/1449) in the model-building",
        "cohort (Table 1)."
      )
    ),
    CONMED_FUROSEMIDE = list(
      description = "Concomitant furosemide co-administration indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained. Univariate CL reduction of 34% but DeltaOF",
        "= -6.3 (P = 0.012), not significant at the P < 0.01 threshold used",
        "during covariate inclusion. Only 5 subjects received furosemide,",
        "limiting power (Results 'Covariate model' and Discussion). 0.3% of",
        "subjects in the model-building cohort (Table 1)."
      )
    ),
    PDA = list(
      description = "Patent ductus arteriosus indicator. Screened.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained. 10.6% of subjects (153/1449) in the",
        "model-building cohort (Table 1)."
      )
    ),
    IV_VENT = list(
      description = "Invasive ventilation indicator. Screened.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained. 20.8% of subjects in the model-building cohort.",
        "Results 'Covariate model' and Discussion paragraph on respiratory",
        "support."
      )
    ),
    NIV_VENT = list(
      description = "Non-invasive ventilation indicator. Screened.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened; not retained. 59.4% of subjects in the model-building cohort."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1449L,
    n_studies      = 1L,
    n_observations = 3039L,
    age_range      = "GA 24-42 weeks at birth; PNA 0-94 days (median 1 day) at first sample",
    age_median     = "GA 34 weeks; PMA 34.4 weeks; PNA 1 day",
    weight_range   = "0.440-5.510 kg",
    weight_median  = "2.170 kg",
    sex_female_pct = 42.5,
    race_ethnicity = "Not reported (single-centre Swiss cohort).",
    disease_state  = paste(
      "Preterm and term neonates (994 preterm, 455 term) admitted to a Neonatal",
      "Intensive Care Unit and receiving gentamicin for suspected or proven",
      "infection. Within a routine therapeutic drug monitoring (TDM) programme."
    ),
    dose_range     = paste(
      "3 mg/kg per dose IV (December 2006-April 2011) and 4 mg/kg per dose IV",
      "(May 2011-October 2011), each given as a 30-minute IV infusion, mostly",
      "in combination with amoxicillin."
    ),
    regions        = "Switzerland (Centre Hospitalier Universitaire Vaudois, Lausanne; December 2006-October 2011).",
    co_medications = paste(
      "9.4% dopamine; 1.9% indomethacin; 0.3% furosemide; 20.8% invasive",
      "ventilation; 59.4% non-invasive ventilation; 10.6% patent ductus",
      "arteriosus (Table 1)."
    ),
    validation_cohort = paste(
      "Independent external validation set of 69 preterm and term newborns",
      "recruited through TDM January-April 2013 (Table 1), used for accuracy /",
      "precision assessment and dosage-adjustment-method comparison."
    ),
    notes          = paste(
      "Retrospective TDM cohort, December 2006-October 2011. Routine sampling",
      "comprised two concentrations after the first dose: peak (between 0.5",
      "and 1.5 h after the start of infusion) and a 12 h sample (between 11.5",
      "and 12.5 h). 86% of measurements were after the first gentamicin dose",
      "and only 3% beyond 72 h; 98% of measurements were within the first",
      "week of life. Concentrations ranged 0.5-22.1 mg/L (one subject had 29",
      "mg/L after an accidental 10x overdose). Baseline demographics from",
      "Table 1; PK estimates from Table 2 (final model)."
    )
  )

  ini({
    # Structural PK parameters (Fuchs 2014 Table 2 'Final model' column),
    # reported at the cohort median body weight 2.170 kg. All log-transformed.
    lcl <- log(0.089);  label("Clearance at reference BW = 2.170 kg, GA = 34 weeks, PNA = 1 day, no dopamine (L/h)")  # Fuchs 2014 Table 2 final-model CL = 0.089 L/h (SE 1%)
    lvc <- log(0.908);  label("Central volume of distribution at reference BW = 2.170 kg, GA = 34 weeks (L)")          # Fuchs 2014 Table 2 final-model Vc = 0.908 L (SE 2%)
    lq  <- log(0.157);  label("Intercompartmental clearance at reference BW = 2.170 kg (L/h)")                          # Fuchs 2014 Table 2 final-model Q = 0.157 L/h (SE 7%)
    lvp <- log(0.560);  label("Peripheral volume of distribution at reference BW = 2.170 kg (L)")                       # Fuchs 2014 Table 2 final-model Vp = 0.560 L (SE 4%)

    # Allometric body-weight exponents. Fuchs 2014 reports these as fixed
    # (no SE reported in Table 2; canonical allometric exponents 0.75 on
    # clearances and 1 on volumes, per Methods 'Covariate model' citation
    # of the allometric scaling approach).
    e_wt_cl_q  <- fixed(0.75); label("Allometric body-weight exponent shared across CL and Q (unitless)")  # Fuchs 2014 Table 2 theta_CL_BW = theta_Q_BW = 0.75 (fixed)
    e_wt_vc_vp <- fixed(1.0);  label("Allometric body-weight exponent shared across Vc and Vp (unitless)") # Fuchs 2014 Table 2 theta_Vc_BW = theta_Vp_BW = 1 (fixed)

    # Linear (centred-on-median) covariate effects.
    # Form: f = 1 + theta * (cov / cov_ref - 1).
    e_ga_cl         <-  1.870; label("Linear GA effect on CL centred on GA = 34 weeks (unitless)")    # Fuchs 2014 Table 2 theta_CL_GA = 1.870 (SE 3%)
    e_pna_cl        <-  0.054; label("Linear PNA effect on CL centred on PNA = 1 day (unitless)")     # Fuchs 2014 Table 2 theta_CL_PNA = 0.054 (SE 6%)
    e_conmed_dopa_cl <- -0.120; label("Multiplicative dopamine effect on CL (unitless)")              # Fuchs 2014 Table 2 theta_CL_DOPA = -0.120 (SE 22%)
    e_ga_vc         <- -0.922; label("Linear GA effect on Vc centred on GA = 34 weeks (unitless)")    # Fuchs 2014 Table 2 theta_Vc_GA = -0.922 (SE 8%)

    # IIV (log-normal). Fuchs 2014 Table 2 reports BSV as CV%; converted to
    # log-scale variance via omega^2 = log(1 + CV^2). The CL-Vc correlation
    # of 87% is preserved as an off-diagonal element of the OMEGA block.
    #   IIV CL  : 28% CV  -> log(1 + 0.28^2) = 0.0754785
    #   IIV Vc  : 18% CV  -> log(1 + 0.18^2) = 0.0318862
    #   cov(CL, Vc) = 0.87 * sqrt(0.0754785 * 0.0318862) = 0.0426808
    etalcl + etalvc ~ c(0.0754785,
                        0.0426808, 0.0318862)  # Fuchs 2014 Table 2 BSV CL 28% CV, BSV Vc 18% CV, correlation CL-Vc 87%

    # Residual error. Fuchs 2014 Table 2 final-model residual is a combined
    # additive (0.10 mg/L) and proportional (18%) model.
    addSd  <- 0.10;  label("Additive residual error (mg/L)")             # Fuchs 2014 Table 2 final-model additive residual error 0.10 mg/L (SE 24%)
    propSd <- 0.18;  label("Proportional residual error (fraction)")     # Fuchs 2014 Table 2 final-model proportional residual error 18% (SE 1%)
  })

  model({
    # ----- Derived covariate terms (centred on cohort medians) -----
    # Reference values from Methods 'Model-based pharmacokinetic analysis'
    # and Table 1: BW median 2.170 kg, GA median 34 weeks, PNA median 1 day.
    # PNA is supplied in canonical months; the reference 1 day expressed in
    # months is 1 / 30.4375, so the ratio PNA / pna_ref_months reproduces
    # the paper's days-scaled ratio PNA_days / 1.
    pna_ref_months <- 1 / 30.4375

    f_ga_cl  <- 1 + e_ga_cl  * (GA  / 34            - 1)
    f_ga_vc  <- 1 + e_ga_vc  * (GA  / 34            - 1)
    f_pna_cl <- 1 + e_pna_cl * (PNA / pna_ref_months - 1)
    f_dopa_cl <- 1 + e_conmed_dopa_cl * CONMED_DOPA

    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (WT / 2.170) ^ e_wt_cl_q  * f_ga_cl * f_pna_cl * f_dopa_cl
    vc <- exp(lvc + etalvc) * (WT / 2.170) ^ e_wt_vc_vp * f_ga_vc
    q  <- exp(lq)           * (WT / 2.170) ^ e_wt_cl_q
    vp <- exp(lvp)          * (WT / 2.170) ^ e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system: two-compartment IV PK -----
    # Gentamicin was delivered as a 30-min IV infusion in the source paper;
    # the model does not hard-code the infusion duration, so users specify
    # rate / dur per dose in their event table.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central                 - k21 * peripheral1

    # ----- Observation and error -----
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
