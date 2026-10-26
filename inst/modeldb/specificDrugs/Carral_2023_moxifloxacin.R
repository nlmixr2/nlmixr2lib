Carral_2023_moxifloxacin <- function() {
  description <- "Two-compartment population PK simulation model for oral moxifloxacin in virtual adult male outpatients with community-acquired pneumonia (Carral 2023). First-order absorption; CL/F and Q/F scale with (LBM/60)^0.75 and V2/F and V3/F linearly with lean body mass. Taken by the authors from the Grosjean and Urien 2012 population model for their Monte Carlo adherence / relative-forgiveness analysis."
  reference <- paste(
    "Carral N, Lukas JC, Estrade O, Jauregizar N, Morillas H, Suarez E. Non-Adherence in Adult Male",
    "Patients with Community-Acquired Pneumonia: Relative Forgiveness of",
    "Amoxicillin versus Respiratory Fluoroquinolones. Antibiotics (Basel).",
    "2023;12(5):838. doi:10.3390/antibiotics12050838. Structural model taken",
    "from Grosjean P, Urien S. J Clin Pharmacol. 2012;52(3):329-338.",
    "doi:10.1177/0091270011398361.",
    sep = " "
  )
  vignette <- "Carral_2023_cap_nonadherence"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 60 kg. Scales CL/F and Q/F by (LBM/60)^0.75 and V2/F and V3/F by LBM/60 (Appendix C). Generated in the control stream as SIMLBM = 57.5 * exp(ETA(1)) with omega^2 = 0.00035; here it is a data column and the vignette reproduces the draw.",
      source_name = "SIMLBM"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "moxifloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "moxifloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    age_range = "31-65 years (virtual)",
    weight_range = "typical 70 kg body weight, 57.5 kg lean body mass (virtual)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "Virtual adult male outpatients with community-acquired pneumonia (CAP) caused by Streptococcus pneumoniae",
    dose_range = "Oral moxifloxacin 400 mg once daily for 10 days",
    renal_function = "Single virtual group spanning CLcr 70-131 mL/min; moxifloxacin clearance does not depend on CLcr (Methods 4.3; Supplementary Table S2)",
    notes = "Simulation-only study: no patients were enrolled. The PK model is the Grosjean and Urien 2012 moxifloxacin population model as used by the authors (Methods 4.3), run in NONMEM ONLYSIM mode with 10,000 virtual patients."
  )

  ini({
    # Appendix C control stream ($THETA, ADVAN4 TRANS4).
    lcl <- log(10); label("Apparent clearance CL/F at LBM = 60 kg (L/h)") # Appendix C THETA(1) = 10 '1.CL L/h'
    lvc <- log(131); label("Apparent central volume V2/F at LBM = 60 kg (L)") # Appendix C THETA(2) = 131 '2.V2 L'
    lvp <- log(44.2); label("Apparent peripheral volume V3/F at LBM = 60 kg (L)") # Appendix C THETA(3) = 44.2 '3.V3 L'
    lq <- log(4.91); label("Apparent intercompartmental clearance Q/F at LBM = 60 kg (L/h)") # Appendix C THETA(4) = 4.91 '4. Q'
    lka <- log(5.97); label("First-order absorption rate constant KA (1/h)") # Appendix C THETA(5) = 5.97 '6.KA 1/h'
    e_lbm_cl_q <- fixed(0.75); label("Exponent of LBM/60 on CL/F and Q/F (unitless)") # Appendix C TVCL/F = THETA(1)*((SIMLBM/60)**0.75); TVQ likewise (hard-coded)
    e_lbm_vc_vp <- fixed(1); label("Exponent of LBM/60 on V2/F and V3/F (unitless)") # Appendix C TV2/F = THETA(2)*(SIMLBM/60); TVV3/F = THETA(3)*SIMLBM/60

    # Appendix C $OMEGA, in ETA order. ETA(1) (lean body mass, 0.00035)
    # generates the virtual covariate and is reproduced in the vignette
    # cohort; V3/F carries no ETA.
    etalcl ~ 0.045 # Appendix C $OMEGA ETA(2)
    etalvc ~ 0.021 # Appendix C $OMEGA ETA(3)
    etalq ~ 0.25 # Appendix C $OMEGA ETA(4), printed '0 .25'

    propSd <- sqrt(0.029); label("Proportional residual error (fraction)") # Appendix C $SIGMA 0.029 (variance); Y = F + F*EPS(1)
  })

  model({
    cl <- exp(lcl + etalcl) * (LBM / 60)^e_lbm_cl_q
    vc <- exp(lvc + etalvc) * (LBM / 60)^e_lbm_vc_vp
    vp <- exp(lvp) * (LBM / 60)^e_lbm_vc_vp
    q <- exp(lq + etalq) * (LBM / 60)^e_lbm_cl_q
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
