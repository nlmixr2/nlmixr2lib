Carral_2023_levofloxacin <- function() {
  description <- "Two-compartment population PK simulation model for oral levofloxacin in virtual adult male outpatients with community-acquired pneumonia (Carral 2023). First-order absorption; micro-constant parameterisation in which apparent clearance is a linear regression on Cockcroft-Gault creatinine clearance (computed from ideal body weight) and age, central volume and k21 are linear in age. Reparameterised by the authors from the Preston 1998 demographic model for their Monte Carlo adherence / relative-forgiveness analysis."
  reference <- paste(
    "Carral N, Lukas JC, Estrade O, Jauregizar N, Morillas H, Suarez E. Non-Adherence in Adult Male",
    "Patients with Community-Acquired Pneumonia: Relative Forgiveness of",
    "Amoxicillin versus Respiratory Fluoroquinolones. Antibiotics (Basel).",
    "2023;12(5):838. doi:10.3390/antibiotics12050838. Structural model and",
    "covariate regressions taken from Preston SL et al. Antimicrob Agents",
    "Chemother. 1998;42(5):1098-1104. doi:10.1128/AAC.42.5.1098.",
    sep = " "
  )
  vignette <- "Carral_2023_cap_nonadherence"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the Cockcroft-Gault creatinine clearance and the linear regressions for CL/F, V2/F and K32 (Appendix B). Generated in the control stream as SIMAGE = 45 * exp(ETA(2)) with omega^2 = 0.006; here it is a data column and the vignette reproduces the draw.",
      source_name = "SIMAGE"
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Used as the weight term of the Cockcroft-Gault creatinine clearance (Appendix B). Generated in the control stream as SIMIBW = 70.11 * exp(ETA(1)) with omega^2 = 0.0003.",
      source_name = "SIMIBW"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters only through the Cockcroft-Gault creatinine clearance CLcr = (140 - AGE) * IBW / (72 * CREAT) (male form). The printed Appendix B stream hard-codes 0.9 mg/dL (the middle renal subgroup); the paper's subgroups use 0.7, 0.9 and 1.3 mg/dL (Methods 4.1; Supplementary Table S2), so creatinine is a data column here.",
      source_name = "SCr"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    age_range = "31-65 years (virtual; typical 45 years)",
    weight_range = "typical ideal body weight 70.11 kg (virtual)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "Virtual adult male outpatients with community-acquired pneumonia (CAP) caused by Streptococcus pneumoniae",
    dose_range = "Oral levofloxacin 750 mg once daily for 10 days",
    renal_function = "Three virtual subgroups by serum creatinine 0.7 / 0.9 / 1.3 mg/dL: Cockcroft-Gault CLcr 131 +/- 10, 100 +/- 8 and 70 +/- 6 mL/min (Supplementary Table S2)",
    notes = "Simulation-only study: no patients were enrolled. The PK model is the authors' reparameterisation of the Preston 1998 levofloxacin population model (Methods 4.3), run in NONMEM ONLYSIM mode with 10,000 virtual patients per subgroup."
  )

  ini({
    # Appendix B control stream ($THETA, ADVAN4 TRANS1). THETA(1) = 0.7
    # 'Blood creatine mean value' is not referenced in the printed $PK.
    lcl <- log(5.94); label("Constant term of the linear CL/F regression, THETA(2) (L/h)") # Appendix B THETA(2) = 5.94 'Constant for prediction of CL'
    e_crcl_cl <- fixed(0.07); label("Linear slope of CL/F on Cockcroft-Gault CLcr ((L/h)/(mL/min))") # Appendix B TVCL = THETA(2) - 1.486 + 0.07*CLcr - 0.032*SIMAGE (hard-coded)
    e_age_cl <- fixed(-0.032); label("Linear slope of CL/F on age ((L/h)/year)") # Appendix B TVCL regression (hard-coded)
    lvc <- log(72.096); label("Constant term of the linear V2/F regression, THETA(3) (L)") # Appendix B THETA(3) = 72.096 'Regression coefficient for V'
    e_age_vc <- fixed(-0.332); label("Linear slope of V2/F on age (L/year)") # Appendix B TVV2 = THETA(3) - 0.332*SIMAGE + 16.51 (hard-coded)
    lk12 <- log(0.308); label("Constant term of the K23 (central to peripheral) regression, THETA(4) (1/h)") # Appendix B THETA(4) = 0.308; TVK23 = THETA(4) + 0.083
    lk21 <- log(0.430); label("Constant term of the K32 (peripheral to central) regression, THETA(5) (1/h)") # Appendix B THETA(5) = 0.430; TVK32 = THETA(5) + SIMAGE*0.005 - 0.066
    e_age_k21 <- fixed(0.005); label("Linear slope of K32 on age ((1/h)/year)") # Appendix B TVK32 regression (hard-coded)
    lka <- log(2.38); label("First-order absorption rate constant KA (1/h)") # Appendix B THETA(6) = 2.38 'Ka'

    # Appendix B $OMEGA, in ETA order. ETA(1) (ideal body weight, 0.0003)
    # and ETA(2) (age, 0.006) generate the virtual covariates and are
    # reproduced in the vignette cohort, not here.
    etalcl ~ 0.001 # Appendix B $OMEGA ETA(3)
    etalvc ~ 0.04 # Appendix B $OMEGA ETA(4)
    etalk12 ~ 0.04 # Appendix B $OMEGA ETA(5)
    etalk21 ~ 0.04 # Appendix B $OMEGA ETA(6)

    propSd <- sqrt(0.029); label("Proportional residual error (fraction)") # Appendix B $SIGMA 0.029 (variance); Y = F + F*EPS(1)
  })

  model({
    # Cockcroft-Gault creatinine clearance on ideal body weight, male form
    # (Appendix B, with the hard-coded 0.9 mg/dL replaced by CREAT).
    crcl_cg <- (140 - AGE) * IBW / (CREAT * 72)

    # The additive constants -1.486, +16.51, +0.083 and -0.066 are written
    # as literals in the Appendix B $PK block.
    cl <- (exp(lcl) - 1.486 + e_crcl_cl * crcl_cg + e_age_cl * AGE) * exp(etalcl)
    vc <- (exp(lvc) + e_age_vc * AGE + 16.51) * exp(etalvc)
    k12 <- (exp(lk12) + 0.083) * exp(etalk12)
    k21 <- (exp(lk21) + e_age_k21 * AGE - 0.066) * exp(etalk21)
    ka <- exp(lka)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
