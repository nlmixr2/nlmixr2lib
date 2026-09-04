Li_2025_amikacin <- function() {
  description <- "One-compartment population PK model for amikacin in Chinese premature infants receiving therapeutic drug monitoring (Li 2025); body-weight allometric scaling on CL and V, a linear postmenstrual-age maturation factor on CL, and a creatinine-production-rate renal-function factor on CL."
  reference <- paste(
    "Li J-h, Xu J-h, Huang Y, Zhou M, Yang Z-m, Li J-j, Feng Z-t, Zhang Q,",
    "Yu Y-x, Duan L-f, Tang L.",
    "Therapeutic drug monitoring of amikacin in Chinese premature infant:",
    "a population pharmacokinetic analysis and dosage optimization.",
    "BMC Infect Dis. 2026;26:43. doi:10.1186/s12879-025-11747-z.",
    "(Accepted 18 September 2025; copyright 2025; assigned to the 2026",
    "volume. Indexed by EuropePMC as PMC12797814.)",
    "The postmenstrual-age maturation form (Eq. 3) is cited by Li 2025 to",
    "Tod M, Jullien V, Pons G. Facilitation of drug evaluation in children",
    "by population methods and modelling.",
    "Clin Pharmacokinet. 2008;47:231-243.",
    "doi:10.2165/00003088-200847040-00002.",
    "The creatinine-production-rate renal-function form (Eqs. 4 and 5) is",
    "cited by Li 2025 to",
    "Allegaert K, Scheers I, Cossey V, Anderson BJ.",
    "Covariates of amikacin clearance in neonates: the impact of postnatal",
    "age on predictability.",
    "Drug Metab Lett. 2008;2:286-289. doi:10.2174/187231208786734120.",
    sep = " "
  )
  vignette <- "Li_2025_amikacin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE -- amikacin was quantified in SERUM by
  # LC-MS/MS (Li 2025 Methods, "Blood sampling and concentration
  # determination"), and doses are given in mg by 0.5-h IV infusion.
  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Allometric scaling on CL (exponent 0.75) and on V",
        "(exponent 1), both standardised to an adult body mass of 70 kg",
        "(Li 2025 Eq. 2 and its legend: 'Weight_standard is 70 kg for adults;",
        "PWR is the empirical coefficient, for CL is 0.75 and V is 1').",
        "Cohort weight at amikacin administration: median 1.36 kg",
        "(range 0.80-4.00 kg), Li 2025 Table 1."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age in weeks / 4.35 + postnatal age in months)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives two separate terms on CL: (a) the linear",
        "maturation factor F_PMA = 1 + SLPCL * (PMA - 40) (Li 2025 Eq. 3),",
        "and (b) the age dependence of the creatinine production rate",
        "CPR = 516 * exp(Kage * ((PMA - 40) / 52 - 40)) (Li 2025 Eq. 5).",
        "The source paper reports PMA in WEEKS; the canonical PAGE column",
        "carries months, so model() converts back with pma_wk = PAGE * 4.35",
        "before evaluating either term, exactly as in",
        "BoerPerez_2026_piperacillin.R and CohenWolkowiez_2014_piperacillin.R.",
        "Cohort PMA: median 32.1 weeks (range 29.1-39.1), Li 2025 Table 1",
        "(equivalently PAGE median 7.38 months, range 6.69-8.99)."
      ),
      source_name        = "PMA (weeks)"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline serum creatinine at the start of amikacin therapy",
        "(Li 2025 Methods: 'to assess renal function by using the baseline",
        "Scr at the beginning of amikacin administration'), measured by the",
        "Jaffe method. Enters CL as the denominator of the creatinine",
        "clearance CLcr = CPR / Scr (Li 2025 Eq. 4), which is then calibrated",
        "to the adult normal value of 6 L/h per 70 kg. umol/L is required:",
        "the calibration only balances when CPR is in umol/h per 70 kg and",
        "Scr is in umol/L, since 516 / 86 = 6.0 (86 umol/L = 1.0 mg/dL is the",
        "adult normal serum creatinine). Cohort Scr: mean 31.72 +/- 11.06",
        "umol/L (Li 2025 Table 1); Monte Carlo simulation strata spanned",
        "15-22, 23-36 and 37-60 umol/L (Li 2025 Table 5)."
      ),
      source_name        = "Scr"
    )
  )

  # Covariates screened by the stepwise covariate model (Li 2025 Methods,
  # "Pharmacokinetic analysis and model evaluation") but NOT retained in the
  # final model: "Other covariates were systematically analyzed but were not
  # found to be statistically significant predictors of PK parameters"
  # (Li 2025 Results, "Population PK analysis"). No point estimates are
  # reported for any of them, so none can be encoded in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'gender' in the SCM; not retained. Cohort 13/23 (56.52%) male (Li 2025 Table 1)."
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort mean 28.90 +/- 2.53 weeks (Li 2025 Table 1). PMA (PAGE) was retained instead."
    ),
    PNA = list(
      description = "Postnatal age",
      units       = "months",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort mean 29.56 +/- 13.53 days (Li 2025 Table 1)."
    ),
    WT_BIRTH = list(
      description = "Birth weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort median 1.10 kg (range 0.70-4.30) (Li 2025 Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort mean 29.00 +/- 3.79 g/L (Li 2025 Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort median 7 U/L (range 4-410) (Li 2025 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort median 21 U/L (range 10-920) (Li 2025 Table 1)."
    ),
    TBIL = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort mean 80.59 +/- 54.79 umol/L (Li 2025 Table 1)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened as 'urea nitrogen' in the SCM; not retained. Cohort mean 5.31 +/- 3.47 mmol/L (Li 2025 Table 1)."
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort mean 121.04 +/- 19.58 g/L (Li 2025 Table 1)."
    ),
    PLT = list(
      description = "Platelet count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened in the SCM; not retained. Cohort mean 127.00 +/- 91.92 10^9/L (Li 2025 Table 1)."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate (Schwartz equation)",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = paste(
        "Screened as 'GFR' in the SCM; not retained. Li 2025 Eq. 1 computes",
        "eGFR = K * L / SCR with K = 0.33 (premature infants under 1 year),",
        "L = body length (cm) and SCR in mg/dL. eGFR was used to define acute",
        "kidney injury, not to scale CL; the retained renal-function term is",
        "the CPR / Scr creatinine clearance of Eqs. 4-5. Cohort mean",
        "43.49 +/- 3.58 mL/min/1.73m^2 (Li 2025 Table 1)."
      )
    ),
    APGAR1 = list(
      description = "Apgar score at 1 minute",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened as 'apgar score' in the SCM; not retained. Cohort median 7 (range 0-10) (Li 2025 Table 1)."
    ),
    APGAR5 = list(
      description = "Apgar score at 5 minutes",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened as 'apgar score' in the SCM; not retained. Cohort median 8 (range 5-10) (Li 2025 Table 1)."
    ),
    CONMED_IBUPROFEN = list(
      description = "Concomitant ibuprofen administration indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'co-administration of ibuprofen' in the SCM; not retained. No counts reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 23,
    n_studies      = 1,
    n_observations = 54,
    age_range      = "Postmenstrual age 29.1-39.1 weeks (median 32.1); postnatal age 29.56 +/- 13.53 days; gestational age 28.90 +/- 2.53 weeks",
    weight_range   = "0.80-4.00 kg at amikacin administration (median 1.36); 0.92-4.01 kg at the time of sampling (median 1.48)",
    sex_female_pct = 43.5,
    race_ethnicity = "Chinese (not otherwise stratified)",
    disease_state  = "Premature infants (gestational age < 37 weeks) in neonatal intensive care with carbapenem-resistant-organism nosocomial pneumonia; 82.6% also had bloodstream infection and 34.8% suppurative meningitis",
    dose_range     = "10.34-19.70 mg/kg/day (median 14.32), most commonly 15 mg/kg/day split as 7.5 mg/kg every 12 h by 0.5-h IV infusion; median 12 days of therapy",
    regions        = "China (two neonatal intensive care units in Suzhou, Jiangsu: the Affiliated Suzhou Hospital of Nanjing Medical University and the Children's Hospital of Soochow University)",
    ga_range       = "28.90 +/- 2.53 weeks (all < 37 weeks by inclusion criteria)",
    pma_range      = "29.1-39.1 weeks postmenstrual age (median 32.1)",
    pna_range      = "29.56 +/- 13.53 days postnatal age",
    birth_weight   = "0.70-4.30 kg (median 1.10)",
    renal_function = "Serum creatinine 31.72 +/- 11.06 umol/L; Schwartz eGFR 43.49 +/- 3.58 mL/min/1.73m^2; 3 of 23 developed acute kidney injury on days 5-7 of therapy",
    co_medication  = "Vancomycin 43.5%, meropenem 34.8%, piperacillin-tazobactam 21.7%, other 21.7% (all intravenous, given separately from amikacin)",
    notes          = paste(
      "Two-centre retrospective therapeutic-drug-monitoring study; clinical",
      "and demographic records collected January 2021 - December 2022.",
      "54 amikacin serum concentrations from 23 premature infants. Samples",
      "were drawn 1 h after the end of an infusion and 30 min before the next",
      "dose, all after the fifth dose, with additional random inter-dose",
      "times. Assay LLOQ 0.8 ug/mL by LC-MS/MS with isepamicin internal",
      "standard. Observed peak concentration median 17.45 ug/mL (range",
      "4.56-164.72) and trough median 4.07 ug/mL (range 1.01-30.99)",
      "(Li 2025 Tables 1 and 2)."
    )
  )

  ini({
    # ===== Structural PK (Li 2025 Table 4, final model; Eqs. 6 and 7) =====
    # Both typical values are standardised to a 70 kg adult, so they are the
    # value a hypothetical 70 kg subject would have at full PMA maturation
    # and normal adult renal function -- NOT a value for a neonate.
    lcl <- log(1.43);   label("Typical CL standardised to 70 kg (L/h)")  # Li 2025 Table 4: tvCL = 1.43 L/h/70kg (SE 0.14, CV 9.58%, 95% CI 1.15-1.70; bootstrap median 1.42)
    lvc <- log(30.97);  label("Typical V standardised to 70 kg (L)")     # Li 2025 Table 4: tvV  = 30.97 L/70kg (SE 2.89, CV 9.32%, 95% CI 25.17-36.77; bootstrap median 31.09)

    # Allometric exponents (Li 2025 Eq. 2 legend: "PWR is the empirical
    # coefficient, for CL is 0.75 and V is 1"). Declared empirical constants,
    # reported without SE / CI / bootstrap.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT / 70) for CL (unitless)")
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on (WT / 70) for V (unitless)")

    # Postmenstrual-age and creatinine-production constants (Li 2025 Table 4).
    # Both are reported as bare point estimates with "-" in every uncertainty
    # column (SE, CV, 95% CI) and are absent from the bootstrap columns, in
    # contrast to tvCL, tvV, both omegas and stdev0, which all carry a full
    # SE / CV / 95% CI and a bootstrap median. That reporting split is the
    # evidence they were held constant rather than estimated; Kage is in
    # addition a literature constant of the creatinine-production model that
    # Li 2025 cites to Allegaert 2008. See the vignette Errata section.
    e_page_cl  <- fixed(0.032);   label("Linear slope of (PMA - 40 weeks) on CL, SLPCL (per week)")
    e_page_cpr <- fixed(0.00823); label("Age scaling constant of the creatinine production rate, Kage (per year)")

    # ===== IIV (Li 2025 Table 4) =====
    # Exponential (log-normal) random effects: V = tvV * exp(etaV),
    # CL = tvCL * exp(etaCL) (Li 2025 Eqs. 6-7 and Table 3 basic model).
    # The table reports these as omega^2, i.e. VARIANCES on the log scale.
    # 0.16 -> CV = sqrt(exp(0.16) - 1) = 41.7%; 0.15 -> 40.4%.
    etalcl ~ 0.16  # Li 2025 Table 4: omega^2 CL = 0.16 (SE 0.046, 95% CI 0.085-0.23, shrinkage 4.61%; bootstrap median 0.15)
    etalvc ~ 0.15  # Li 2025 Table 4: omega^2 V  = 0.15 (SE 0.047, 95% CI 0.073-0.23, shrinkage 10.35%; bootstrap median 0.15)

    # ===== Residual error (Li 2025 Table 4) =====
    # "residual variability was fitted with an additive residual error model"
    # (Li 2025 Results, "Population PK analysis"); the Phoenix NLME additive
    # residual SD is reported as stdev0, in the concentration unit ug/mL.
    addSd <- 0.92; label("Additive residual SD (mg/L)")  # Li 2025 Table 4: stdev0 = 0.92 (SE 0.10, CV 11.14%, 95% CI 0.72-1.13; bootstrap median 0.91)

    # ===== Reference / centering constants =====
    wt_ref    <- 70;  label("Adult standard body mass for allometric scaling (kg)")                  # Li 2025 Eq. 2 legend: "Weight_standard is 70 kg for adults"
    pma_ref   <- 40;  label("Reference postmenstrual age at term birth (weeks)")                     # Li 2025 Eqs. 3 and 5: (PMA - 40)
    age_ref   <- 40;  label("Adult reference age of the creatinine production rate (years)")         # Li 2025 Eq. 5: the "- 40" outside the (PMA - 40)/52 years conversion
    cpr_adult <- 516; label("Adult male creatinine production rate (umol/h per 70 kg)")              # Li 2025 Methods: "the CPR was 516 umol.h-1" for "an individual mass of 70 kg in adult males"
    clcr_ref  <- 6;   label("Adult normal creatinine clearance used to calibrate renal function (L/h per 70 kg)")  # Li 2025 Eq. 2-5 legend: "calibrated by CLcr of 6 L/h.per (70 kg)-1"
  })

  model({
    # ----- Derived covariate terms -----
    # Canonical PAGE is postmenstrual age in MONTHS; Li 2025 writes every
    # equation in postmenstrual WEEKS, so convert back before use.
    pma_wk <- PAGE * 4.35

    # Li 2025 Eq. 3 -- linear postmenstrual-age maturation factor on CL:
    #   F_PMA = [1 + SLPCL * (PMA - 40)]
    # F_PMA = 1 at term birth (PMA = 40 weeks) and falls below 1 in preterm
    # infants (e.g. 0.747 at the cohort median PMA of 32.1 weeks).
    fpma <- 1 + e_page_cl * (pma_wk - pma_ref)

    # Li 2025 Eq. 5 -- creatinine production rate (umol/h per 70 kg):
    #   CPR = 516 * exp(Kage * ((PMA - 40) / 52 - 40))
    # The inner term is age in YEARS relative to term birth, (PMA - 40) / 52,
    # minus the 40-year adult reference age at which CPR = 516.
    age_yr <- (pma_wk - pma_ref) / 52
    cpr    <- cpr_adult * exp(e_page_cpr * (age_yr - age_ref))

    # Li 2025 Eq. 4 -- renal function. CPR (umol/h per 70 kg) divided by serum
    # creatinine (umol/L) is creatinine clearance in L/h per 70 kg; Li 2025
    # states this is "calibrated by CLcr of 6 L/h per (70 kg)^-1", so the
    # dimensionless factor that multiplies CL is CLcr / 6, equal to 1 at
    # normal adult renal function (516 / 86 umol/L = 6.0 L/h per 70 kg).
    # Printed Eq. 6 abbreviates this to "CPR / Creatinine" and omits the
    # calibration; see the vignette Errata for the numeric evidence that the
    # calibration is required.
    rf <- (cpr / CREAT) / clcr_ref

    # ----- Individual PK parameters (Li 2025 Eqs. 6 and 7) -----
    #   CL = 1.43 * (WT / 70)^0.75 * F_PMA * RF * exp(etaCL)
    #   V  = 30.97 * (WT / 70)^1                * exp(etaV)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl * fpma * rf
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # One compartment with first-order elimination; amikacin is given as a
    # 0.5-h IV infusion, so dosing goes directly into central.
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Serum amikacin concentration: dose in mg, vc in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
