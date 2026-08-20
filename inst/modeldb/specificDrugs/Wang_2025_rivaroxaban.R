Wang_2025_rivaroxaban <- function() {
  description <- "One-compartment population PK model for rivaroxaban in Chinese patients with non-valvular atrial fibrillation, with an AST/ALT-ratio power effect on CL/F and V/F (Wang 2025)"
  reference <- "Wang F, Li Z, Huang Y, Liu Q, Zhao L, Wang H, Gao H, Chen M, Lin Y, Li X, Chen M. Effect of ABCB1 SNP polymorphisms on the plasma concentrations and clinical outcomes of rivaroxaban in Chinese NVAF patients: a population pharmacokinetic-based study. Front Pharmacol. 2025;16:1574949. doi:10.3389/fphar.2025.1574949"
  vignette <- "Wang_2025_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Rivaroxaban is given orally (5-20 mg once daily,
  # Methods "Subjects and therapeutic interventions") and measured in plasma
  # by HPLC-MS/MS (Methods "Measurement of the collected rivaroxaban plasma
  # concentrations").
  compartmentData <- list(
    depot   = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description        = "Serum aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters the model only through the AST/ALT ratio (the De Ritis ratio), which was the single covariate retained on both CL/F and V/F by forward inclusion / backward elimination. The ratio is normalized by its population median of 1.188 (Results, text following Equations 7 and 8). Baseline AST 21 U/L (range 6.2-197) per Table 1; the observed AST/ALT ratio ranged 0.37-6.5 (Discussion).",
      source_name        = "AST"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Denominator of the AST/ALT ratio. Baseline ALT 18 U/L (range 1.5-82) per Table 1. AST and ALT are carried as two separate canonical columns rather than a single pre-computed ratio column so that the ratio is formed inside model() exactly as Equations 7 and 8 write it.",
      source_name        = "ALT"
    )
  )

  # Screened in the covariate analysis but NOT retained in the final model
  # (Results, "Population pharmacokinetic model": "The covariates analyzed
  # included age, sex, body weight (BW), body mass index (BMI), albumin (ALB),
  # bilirubin (BIL), alanine aminotransferase (ALT), aspartate aminotransferase
  # (AST), the AST/ALT ratio, serum creatinine (SCR), and estimated glomerular
  # filtration rate (eGFR)." ... "This iterative process identified AST/ALT
  # ratios as the sole covariate significantly influencing CL/F and V/F").
  # No point estimate is published for any of these, so none can be encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; showed an effect on CL/F or V/F in the diagnostic covariate plots and entered the initial model, but was eliminated during backward elimination. Median 73 years (range 36-94), Table 1."
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened but not retained; 88 of 228 participants (38.6%) were female, Table 1."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained. Median 65 kg (range 33.5-99), Table 1. Unlike prior rivaroxaban popPK analyses, this study did not exclude patients below 45 kg (Discussion)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained. Median 23.7 kg/m^2 (range 13.6-36), Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained. Median 41 g/L (range 27-54), Table 1."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Reported as 'BIL' in Table 1. Screened but not retained. Median 11.6 umol/L (range 2.1-53.6)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Reported as 'SCR' in Table 1 in mg/dL (not the SI umol/L). Screened but not retained. Median 0.89 mg/dL (range 0.31-4.21)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Reported as 'eGFR' in Table 1 and computed with the Chronic Kidney Disease Epidemiology Collaboration equation (Kong 2013). Screened; entered the initial model but was eliminated during backward elimination. Median 79.2 mL/min (range 13.3-130.4). Reported in mL/min rather than the CRCL canonical mL/min/1.73 m^2."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 228,
    n_studies      = 1,
    n_observations = 287,
    age_range      = "36-94 years",
    age_median     = "73 years",
    weight_range   = "33.5-99 kg",
    weight_median  = "65 kg",
    sex_female_pct = 38.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "non-valvular atrial fibrillation (NVAF)",
    dose_range     = "5, 7.5, 10, 15 or 20 mg orally once daily",
    regions        = "China (single center: Fujian Provincial Hospital, Fuzhou)",
    renal_function = "eGFR (CKD-EPI) median 79.2 mL/min, range 13.3-130.4 mL/min",
    hepatic_function = "AST median 21 U/L (6.2-197); ALT median 18 U/L (1.5-82); AST/ALT ratio median 1.188, range 0.37-6.5",
    notes          = "Prospective single-center study; Table 1 baseline demographics. Sampling was sparse and opportunistic: 287 concentrations from 228 patients, median 1 sample per patient (range 1-3; mean 1.26 +/- 0.54), taken from residual blood after routine biochemistry. Risk scores: CHA2DS2-VASc median 4 (2-10), HAS-BLED median 2 (1-5). Four ABCB1 SNPs (3435C>T, 1236C>T, 2677G>T/A, c.2482-2236C>T) were genotyped but were NOT covariates in the PK model; their effects were assessed post hoc on model-simulated Cmax/D and Ctrough/D (Tables 6-8)."
  )

  ini({
    # Structural parameters. Time in h, dose in mg, CL/F in L/h, V/F in L.
    # ka was not estimated: "the absorption rate constant (k_a) was fixed at
    # 0.617 h-1 based on a prior PPK study involving Japanese patients
    # (Kaneko et al., 2013)" (Results); Table 5 marks it "(Freeze)".
    lka <- fixed(log(0.617)); label("Absorption rate constant (ka, 1/h), taken from Kaneko 2013")  # Table 5 row "0.617 ka (h-1) (Freeze)"; Results, Population pharmacokinetic model
    lcl <- log(5.64);         label("Apparent clearance (CL/F, L/h)")                              # Table 5 CL/F estimate 5.64 (RSE 5.49%); Equation 7
    lvc <- log(41.7);         label("Apparent volume of distribution (V/F, L)")                    # Table 5 V/F estimate 41.7 (RSE 7.58%); Equation 8

    # Covariate effects. Equation 6 gives the continuous-covariate form
    #   P_i = P * (Cov / Cov_median)^theta * exp(eta_i)
    # and Equations 7-8 instantiate it with Cov = AST/ALT and
    # Cov_median = 1.188:
    #   CL/F (L/h) = 5.64 * (AST/ALT / 1.188)^-0.074
    #   V/F  (L)   = 41.7 * (AST/ALT / 1.188)^ 0.213
    e_astalt_cl <- -0.074; label("AST/ALT ratio power exponent on CL/F (unitless)")  # Table 5 fCL/F-AST/ALT = -0.074 (RSE -14.74%); Equation 7
    e_astalt_vc <-  0.213; label("AST/ALT ratio power exponent on V/F (unitless)")   # Table 5 fV/F-AST/ALT  =  0.213 (RSE  15.87%); Equation 8

    # IIV, log-normal (Methods: "a log-normal distribution was assumed for the
    # interindividual variability (IIV) of the PK parameters", Equation 1).
    # Table 5 reports IIV as CV%; mapped via omega^2 = log(CV^2 + 1):
    #   CL/F CV 34.64% -> omega^2 = log(0.3464^2 + 1) = 0.113322
    #   V/F  CV 19.81% -> omega^2 = log(0.1981^2 + 1) = 0.038490
    # A full OMEGA block was evaluated and rejected: "the correlation between
    # CL/F and V/F was found to be negligible ... Thus, a diagonal OMEGA matrix
    # was retained in the final model to maintain parsimony" (Results).
    etalcl ~ 0.113322  # Table 5 IIV(CV%) on CL/F = 34.64 (eta shrinkage 15.8%)
    etalvc ~ 0.038490  # Table 5 IIV(CV%) on V/F  = 19.81 (eta shrinkage 22.4%)

    # Residual error. Table 4 selected the one-compartment + proportional-error
    # combination (lowest OFV, 2661 before the covariate was added).
    propSd <- 0.71; label("Proportional residual error (fraction)")  # Table 5 sigma = 0.71 (RSE 9.61%, epsilon shrinkage 18.7%)
  })
  model({
    # Covariate ratio, exactly as written in Equations 7 and 8: the individual
    # AST/ALT ratio normalized by the population median AST/ALT of 1.188
    # (Results, text following Equations 7 and 8).
    astalt_ratio <- (AST / ALT) / 1.188

    # Individual parameters (Equation 6: power covariate x log-normal IIV).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * astalt_ratio^e_astalt_cl
    vc <- exp(lvc + etalvc) * astalt_ratio^e_astalt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x1000 converts to the reported ng/mL
    # (equivalently ug/L, the unit used in the goodness-of-fit narrative).
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
