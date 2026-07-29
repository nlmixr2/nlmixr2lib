`Boer-Perez_2026_piperacillin` <- function() {
  description <- "One-compartment population PK model for piperacillin in preterm and term neonates with severe infections (Boer-Perez 2026); body-weight allometric scaling, sigmoidal postmenstrual-age maturation on CL fixed from Rhodin 2009, and a power effect of serum creatinine on CL."
  reference <- paste(
    "Boer-Perez FS, Lima-Rogel V, Romano-Moreno S, Mejia-Elizondo AR,",
    "Medellin-Garibay SE, Schaiquevich P, Noyola-Cherpitel DE,",
    "Rodriguez-Baez AS, Rodriguez-Pinal CJ, Milan-Segovia RC.",
    "Population pharmacokinetics and dose optimization of",
    "piperacillin-tazobactam in premature and term neonates with severe",
    "infections.",
    "Antimicrob Agents Chemother. 2026;70(1):e00998-25.",
    "doi:10.1128/aac.00998-25.",
    "Maturation parameters (TM50 = 47.7 weeks, Hill = 3.4) fixed from",
    "Rhodin MM, Anderson BJ, Peters AM, Coulthard MG, Wilkins B,",
    "Cole M, Chatelut E, Grubb A, Veal GJ, Keir MJ, Holford NHG.",
    "Human renal function maturation: a quantitative description using",
    "weight and postmenstrual age.",
    "Pediatr Nephrol. 2009;24(1):67-76. doi:10.1007/s00467-008-0997-5.",
    sep = " "
  )
  vignette <- "Boer-Perez_2026_piperacillin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Allometric scaling on CL (exponent 0.75 fixed) and on V",
        "(exponent 1 fixed) with reference body weight 1.76 kg, the median of",
        "the model-development cohort (Boer-Perez 2026 Table 1)."
      ),
      source_name        = "BW"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the Rhodin 2009 sigmoidal renal-maturation",
        "function on CL (TM50 = 47.7 weeks, Hill = 3.4, both fixed). The",
        "source paper reports PMA in weeks; this model converts the canonical",
        "PAGE (months) back to weeks via pma_wk = PAGE * 4.35 so the published",
        "TM50 and Hill apply unchanged."
      ),
      source_name        = "PMA (weeks)"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on CL: (CREAT / 0.4)^theta_SCr with theta_SCr =",
        "-0.635 (Boer-Perez 2026 Table 2). Reference 0.4 mg/dL is the median",
        "of the model-development cohort (Table 1). Source paper column SCr,",
        "quantified by an enzymatic assay on an automated analyzer."
      ),
      source_name        = "SCr"
    )
  )

  population <- list(
    n_subjects     = 25,
    n_studies      = 1,
    age_range      = "Postnatal age 6-28 days; postmenstrual age 28.1-43.3 weeks; gestational age 26-41.1 weeks",
    age_median     = "Postnatal age 14 days; postmenstrual age 36.0 weeks; gestational age 34.2 weeks",
    weight_range   = "0.89-3.57 kg",
    weight_median  = "1.76 kg",
    sex_female_pct = 56,
    race_ethnicity = "Not reported",
    disease_state  = "Preterm and term neonates with severe infections (most commonly late-onset sepsis, necrotising enterocolitis Bell stage II or higher, and healthcare-associated pneumonia)",
    dose_range     = "100 mg/kg piperacillin (with 12.5 mg/kg tazobactam in 8:1 ratio) every 8 or 12 hours by IV infusion (0.5 or 1 h), per Neofax recommendations",
    regions        = "Mexico (single neonatal intensive care unit at Hospital Central, San Luis Potosi)",
    ga_range       = "26-41.1 weeks (4% extremely preterm, 52% moderate-late preterm, 44% term)",
    pma_range      = "28.1-43.3 weeks postmenstrual age",
    pna_range      = "6-28 days postnatal age",
    creat_range    = "0.2-0.9 mg/dL serum creatinine (cohort median 0.4 mg/dL)",
    crcl_range     = "15-103 mL/min/1.73m^2 (Schwartz formula; cohort median 32 mL/min/1.73m^2)",
    notes          = paste(
      "25 neonates from a single neonatal intensive care unit, all admitted",
      "between September 2020 and May 2024 with PNA <= 29 days, BW >= 850 g",
      "and hematocrit >= 30%. 65 piperacillin plasma concentrations (median",
      "3 samples per patient, opportunistic sampling, range 5.7-302.6 mg/L)",
      "with 60% obtained during the first two weeks of life. An additional",
      "8 neonates (17 plasma samples) formed the external-validation cohort.",
      "Boer-Perez 2026 Tables 1 and demographics paragraph in Methods."
    )
  )

  ini({
    # ===== Structural PK (Boer-Perez 2026 Table 2 final model) =====
    # Reference subject: BW = 1.76 kg, SCr = 0.4 mg/dL, full PMA-driven CL maturation.
    lcl <- log(0.748);  label("Typical CL at BW = 1.76 kg, SCr = 0.4 mg/dL, full PMA maturation (L/h)") # Boer-Perez 2026 Table 2: theta_CL = 0.748 L/h (RSE 8%)
    lvc <- log(0.866);  label("Typical V at BW = 1.76 kg (L)")                                          # Boer-Perez 2026 Table 2: theta_V  = 0.866 L (RSE 8%)

    # Allometric / covariate exponents (Boer-Perez 2026 Table 2 covariate equation; eq. 3 in Methods)
    e_wt_cl    <- 0.75;   label("Allometric exponent on CL (unitless, fixed)")  # Methods: k = 0.75 for CL (eq. 3)
    e_wt_vc    <- 1.0;    label("Allometric exponent on V (unitless, fixed)")   # Methods: k = 1 for V (eq. 3)
    e_creat_cl <- -0.635; label("Power exponent on (CREAT / 0.4) for CL (unitless)")  # Boer-Perez 2026 Table 2: theta_SCr = -0.635 (RSE 36%)

    # Maturation Hill function constants (Boer-Perez 2026 Methods, fixed from Rhodin 2009, ref 35)
    pma50_cl <- 47.7; label("Postmenstrual age at half-maximal CL maturation (weeks)")  # Boer-Perez 2026 Methods + Discussion: TM50 = 47.7 weeks (Rhodin 2009)
    h_pma_cl <- 3.4;  label("Hill coefficient on PMA for CL maturation (unitless)")     # Boer-Perez 2026 Methods + Discussion: Hill = 3.4 (Rhodin 2009)

    # Reference covariate values (Boer-Perez 2026 Table 1 cohort medians, used in Table 2 covariate equation)
    bw_ref    <- 1.76; label("Reference body weight (kg, cohort median)")
    creat_ref <- 0.4;  label("Reference serum creatinine (mg/dL, cohort median)")

    # ===== IIV (Boer-Perez 2026 Table 2 final model) =====
    # Exponential (log-normal) IIV; CV%-to-variance: omega^2 = log(1 + CV^2)
    # 38.3% -> log(1 + 0.383^2) = 0.13688
    # 37.7% -> log(1 + 0.377^2) = 0.13288
    etalcl ~ 0.13688  # Table 2: omega_CL 38.3% CV (shrinkage 2%)
    etalvc ~ 0.13288  # Table 2: omega_V  37.7% CV (shrinkage 6%)

    # ===== Residual error (Boer-Perez 2026 Table 2: 11.4% CV proportional) =====
    propSd <- 0.114; label("Proportional residual error (fraction)")  # Table 2: sigma_proportional 11.4% CV
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # Rhodin 2009 TM50 = 47.7 weeks and Hill = 3.4 apply unchanged.
    pma_wk <- PAGE * 4.35

    # Sigmoidal renal-maturation factor on CL (Boer-Perez 2026 eq. 4):
    #   F_PMA = PMA^Hill / (PMA^Hill + TM50^Hill)
    fpma <- pma_wk^h_pma_cl / (pma50_cl^h_pma_cl + pma_wk^h_pma_cl)

    # Power effect of SCr on CL (Boer-Perez 2026 Table 2 covariate equation):
    #   (SCr / SCr_median)^theta_SCr  with SCr_median = 0.4 mg/dL.
    # theta_SCr = -0.635, so higher SCr (worse renal function) reduces CL.
    creat_factor <- (CREAT / creat_ref)^e_creat_cl

    # ----- Individual PK parameters (Boer-Perez 2026 Table 2 final model) -----
    cl <- exp(lcl + etalcl) * (WT / bw_ref)^e_wt_cl * fpma * creat_factor
    vc <- exp(lvc + etalvc) * (WT / bw_ref)^e_wt_vc

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # IV piperacillin-tazobactam dosed into the central compartment (no depot).
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Plasma piperacillin concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
