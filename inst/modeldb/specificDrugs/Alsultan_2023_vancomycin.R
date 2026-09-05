Alsultan_2023_vancomycin <- function() {
  description <- "One-compartment population PK model for intravenous vancomycin in very low birth weight neonates, with allometric body weight on CL and V (exponents fixed at 0.75 and 1), sigmoidal Hill maturation on CL by postmenstrual age, and a power effect of serum creatinine on CL (Alsultan 2023)."
  reference <- "Alsultan A, Al Munjem MF, Atiq KM, Aljehani ZK, Al Muqati H, Almohaizeie A, Ballal DA, Refaei TM, Al Jeraisy M, Assiri A, Abouelkheir M. Population pharmacokinetics of vancomycin in very low birth weight neonates. Front Pediatr. 2023;11:1093171. doi:10.3389/fped.2023.1093171"
  vignette <- "Alsultan_2023_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin is dosed as a short intravenous infusion
  # directly into the sampled plasma compartment (Methods "Patients and data
  # collection": trough concentrations are serum vancomycin).
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Current total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Scaled to the cohort MEDIAN bodyweight of 0.93 kg (Results 'Population pharmacokinetics': 'Bodyweight was scaled to the median bodyweight of 0.93 kg'), not to the 70 kg adult standard, so the reported CL and V thetas are already neonate-sized. Allometric exponents were FIXED, not estimated: CL ~ (WT/0.93)^0.75 and V ~ (WT/0.93)^1 (Methods 'Covariates', citing Anderson & Holford). Cohort current weight 1.0 kg (SD 0.29, range 0.46-1.7) in training and 1.1 kg (SD 0.3, range 0.5-2.2) in validation (Table 2); the simulation dataset spanned 0.46-2.2 kg (Results 'Simulation').",
      source_name        = "Weight"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months. Alsultan 2023 writes the sigmoidal Hill maturation function on CL directly in weeks (TMA50 = 26.3 weeks, Hill = 4.42; Results equation and Table 3), and the register explicitly permits a weeks-scaled PAGE for models whose source equations are written that way. The maturation term is NOT normalized to a reference PMA -- it is the bare Hill fraction PMA^4.42/(PMA^4.42 + 26.3^4.42), which equals 0.487 at PMA 26 weeks and approaches 1 only well beyond the studied range. Cohort PMA 29.8 weeks (SD 3.15, range 22-39) in training and 30.7 (SD 3.4, range 24-42) in validation (Table 2); the simulation dataset spanned 22-42 weeks (Results 'Simulation'). Table 1 footnote a defines PMA as gestational age plus postnatal age.",
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "mg/dL, NOT the SI umol/L (divide umol/L by 88.4 to convert). Time-varying; enters CL as the RECIPROCAL power ratio (0.6 / CREAT)^0.48 with reference 0.6 mg/dL, so clearance FALLS as creatinine rises. This is equivalent to (CREAT / 0.6)^(-0.48); the model keeps the paper's printed orientation so that the Table 3 estimate 0.48 appears verbatim. Values from the two centers using the enzymatic assay were converted to the Jaffe scale by the authors before modelling, as Jaffe = 0.122 + enzymatic / 1.05 (Methods 'Analytical assay', citing Srivastava 2009); a user supplying enzymatic creatinine should apply the same conversion. Cohort 0.65 mg/dL (SD 0.22, range 0.2-1.5) in training and 0.62 (SD 0.23, range 0.15-1.4) in validation (Table 2). The model does not apply above 1.2 mg/dL: only five patients exceeded it, the Monte Carlo simulations were capped there, and the Discussion states 'our model does not apply to this population'.",
      source_name        = "Scr"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units       = "days",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate search (Methods 'Covariates') but not retained; postmenstrual age carried the maturation signal. Cohort 10.7 days (SD 7.5, range 1-30) in training (Table 2)."
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "weeks",
      type        = "continuous",
      notes       = "Screened but not retained (Methods 'Covariates'); it is a component of the retained PAGE. Cohort 28.0 weeks (SD 2.9, range 22-35) in training (Table 2)."
    ),
    WT_BIRTH = list(
      description = "Birth weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened via the VLBW-vs-ELBW (birth weight < 1.0 kg) contrast (Methods 'Covariates') but not retained; current weight was the retained size descriptor. Cohort 0.95 kg (SD 0.27, range 0.46-1.5) in training (Table 2).",
      source_name = "Birth weight"
    ),
    HT = list(
      description = "Body length / height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Methods 'Covariates'); no cohort summary is reported."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "",
      type        = "binary",
      notes       = "Screened but not retained (Methods 'Covariates'). Sex was missing for 22% of the training and 18% of the validation cohort (Table 2), which limited its power as a covariate."
    ),
    DIS_CHD = list(
      description = "Congenital heart disease indicator",
      units       = "",
      type        = "binary",
      notes       = "Screened but not retained (Methods 'Covariates'). Present in 26% of training and 37% of validation patients (Table 2)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 236L,
    n_studies      = 1L,
    age_range      = "1-30 days postnatal; postmenstrual age 22-42 weeks",
    age_median     = "10.7 days postnatal (mean), 29.8 weeks postmenstrual (mean), training set",
    weight_range   = "0.46-2.2 kg current body weight; 0.46-1.5 kg birth weight",
    weight_median  = "1.0 kg current body weight (mean), training set",
    sex_female_pct = 32,
    race_ethnicity = "Not reported; single-country cohort recruited in Saudi Arabia",
    disease_state  = "Very low birth weight neonates (birth weight < 1.5 kg, 58% of them extremely low birth weight < 1.0 kg) admitted to a NICU and treated with intravenous vancomycin for proven or suspected MRSA or methicillin-resistant coagulase-negative staphylococcal infection; 36% had a culture-confirmed infection. Patients with renal failure or on haemodialysis were excluded.",
    dose_range     = "Total daily dose 22 mg/kg (SD 8, range 7.5-55) in the training set and 24.1 mg/kg (SD 11.5, range 9-68) in the validation set; initial doses followed NeoFax or Lexicomp neonatal guidance (Table 1) with subsequent trough-guided adjustment",
    regions        = "Saudi Arabia (six centers: King Saud University Medical City, King Faisal Specialist Hospital Riyadh and Jeddah, National Guard Hospital Riyadh, Prince Sultan Military Medical City, Armed Forces Hospital Southern Region)",
    renal_function = "Renal failure and haemodialysis were exclusion criteria. Serum creatinine 0.65 mg/dL (SD 0.22, range 0.2-1.5) in training. Only five patients had creatinine > 1.2 mg/dL, so the simulations -- and the model's stated domain of applicability -- were capped at 1.2 mg/dL (Results 'Simulation'; Discussion).",
    notes          = "Retrospective multicenter observational study. The 236 neonates were split randomly 70/30 into a training set (n = 162, 214 concentrations) and an external validation set (n = 74, 97 concentrations); the parameters carried here are from the training-set fit (Table 3), which the authors then evaluated on the validation set by VPC and NPDE (mean NPDE 0.043, SD 1.1). Baseline demographics per Table 2. Sampling was sparse -- 1-2 samples per patient, predominantly steady-state troughs drawn 30 min before the next dose, with peaks (1 h after end of infusion) collected only at one of the six centers -- so the volume of distribution is far less well informed than clearance; the authors list this among the study limitations."
  )

  ini({
    # Structural parameters -- Alsultan 2023 Table 3, in the parameterisation of
    # the two Results equations. The reference subject is a neonate at the cohort
    # MEDIAN weight of 0.93 kg with serum creatinine 0.6 mg/dL.
    #
    # NOTE on lcl: 0.09 L/h is the theta multiplying the Hill maturation
    # FRACTION, i.e. the fully-mature (PMA -> infinity) clearance at 0.93 kg and
    # Scr 0.6, not the clearance at any PMA in the studied range. The Results
    # text calls 0.09 L/h "the typical Cl value for a VLBW neonate weighing
    # 0.93 kg, PMA equal to 26 weeks, and Scr of 0.6 mg/dl", which conflicts with
    # the printed equation directly below it: that equation gives
    # 0.09 * 26^4.42/(26^4.42 + 26.3^4.42) = 0.0439 L/h at PMA 26 weeks. The
    # equation is the correct reading -- it reproduces the Table 4 Monte Carlo
    # AUC0-24 means to a median 11% (the text reading is uniformly 46% low, a
    # factor of exactly 1/0.487). See the vignette Errata.
    lcl <- log(0.09); label("Fully mature clearance at 0.93 kg body weight and 0.6 mg/dL serum creatinine (L/h)")  # Alsultan 2023 Table 3 (Cl = 0.09 L/h, RSE 11%) and Results Cl equation
    lvc <- log(0.81); label("Volume of distribution at 0.93 kg body weight (L)")                                   # Alsultan 2023 Table 3 (V = 0.81 L, RSE 6.6%) and Results V equation

    # Allometric exponents on current body weight, referenced to the cohort
    # median 0.93 kg. Both were FIXED, not estimated -- Results: "the exponents
    # for the effect of weight on V and Cl were fixed at 1 and 0.75"; Methods
    # "Covariates" cites Anderson & Holford. Neither appears in Table 3, which
    # lists only estimated parameters with an RSE.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT / 0.93 kg) for CL (unitless)")  # Alsultan 2023 Methods "Covariates" and Results (fixed)
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT / 0.93 kg) for V (unitless)")   # Alsultan 2023 Methods "Covariates" and Results (fixed)

    # Sigmoidal Hill maturation of clearance with postmenstrual age, in WEEKS.
    # Methods "Covariates" gives the form F_PMA = PMA^hill / (PMA^hill + TMA50^hill).
    pma_tm50 <- 26.3; label("Postmenstrual age at 50% of mature CL (TMA50, weeks)")     # Alsultan 2023 Table 3 (TMA50 = 26.3 weeks, RSE 7%)
    pma_hill <- 4.42; label("Hill coefficient for sigmoidal CL maturation (unitless)")  # Alsultan 2023 Table 3 (Hill coefficient for clearance = 4.42, RSE 19%)

    # Serum-creatinine effect on clearance. The printed Results equation uses the
    # RECIPROCAL ratio (0.6 / Scr)^0.48, so a POSITIVE exponent means clearance
    # falls as creatinine rises; this is identical to (CREAT / 0.6)^(-0.48).
    e_creat_cl <- 0.48; label("Power exponent on (0.6 mg/dL / CREAT) for CL (unitless)")  # Alsultan 2023 Table 3 (Exponent of Scr on Cl = 0.48, RSE 19%) and Results Cl equation

    # IIV. Table 3 reports "IIV V 24%" and "IIV Cl 28%" as percentages, taken as
    # the coefficient of variation of the log-normal random effect declared in
    # Methods "Base model" (p_i = p * exp(eta_i), eta ~ N(0, omega^2)); converted
    # with omega^2 = log(1 + CV^2). The alternative reading -- that the printed
    # percentages ARE omega x 100 -- gives omega 0.2400 / 0.2800 instead of
    # 0.2367 / 0.2747, a difference of under 1.5% relative, so nothing downstream
    # turns on the choice. No CL-V correlation is reported, so the etas are
    # independent.
    etalcl ~ 0.075474  # Alsultan 2023 Table 3 (IIV Cl = 28%, RSE 11%); log(1 + 0.28^2) = 0.075474
    etalvc ~ 0.056011  # Alsultan 2023 Table 3 (IIV V  = 24%, RSE 26%); log(1 + 0.24^2) = 0.056011

    # Residual error. Results: "best described using ... a proportional error
    # model"; Table 3 footnote b defines the reported 0.3 as the proportional
    # error, on Monolix's y = f * (1 + b * eps) scale.
    propSd <- 0.3; label("Proportional residual error (fraction)")  # Alsultan 2023 Table 3 (Residual variability b = 0.3, RSE 9.9%)
  })
  model({
    # 1. Derived covariate terms -- Alsultan 2023 Results, Cl and V equations.
    #    Reference constants: median bodyweight 0.93 kg, reference serum
    #    creatinine 0.6 mg/dL. PAGE is in weeks for this model (see
    #    covariateData$PAGE$notes).
    maturation_cl <- PAGE^pma_hill / (PAGE^pma_hill + pma_tm50^pma_hill)
    creat_cl      <- (0.6 / CREAT)^e_creat_cl

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (WT / 0.93)^e_wt_cl * creat_cl * maturation_cl
    vc <- exp(lvc + etalvc) * (WT / 0.93)^e_wt_vc

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system -- one compartment, linear elimination, intravenous only
    #    (Results: "best described using a one-compartment model with linear
    #    elimination"; vancomycin was given by intermittent IV infusion).
    d/dt(central) <- -kel * central

    # 5. Observation and error. Dose in mg, V in L -> mg/L = ug/mL.
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
