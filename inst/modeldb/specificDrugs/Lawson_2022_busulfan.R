Lawson_2022_busulfan <- function() {
  description <- "Two-compartment IV PK model for once-daily busulfan in pediatric hematopoietic stem cell transplant recipients with allometric normal-fat-mass (NFM) scaling, postmenstrual-age maturation on CL, and a time-associated within-treatment-course CL decline (Lawson 2022)."
  reference <- "Lawson R, Staatz CE, Fraser CJ, et al. Population pharmacokinetic model for once-daily intravenous busulfan in pediatric subjects describing time-associated clearance. CPT Pharmacometrics Syst Pharmacol. 2022;11(8):1002-1017. doi:10.1002/psp4.12809"
  vignette <- "Lawson_2022_busulfan"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight (TBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline in the source analysis. Combined with FFM to compute parameter-specific normal fat mass (NFM = FFM + Ffat * (WT - FFM)) for allometric scaling. Reference adult TBW corresponds to 70 kg (Lawson 2022 Methods, Equation 2 and Table 3 footnote).",
      source_name        = "TBW"
    ),
    FFM = list(
      description        = "Fat-free mass (Al-Sallami 2015 pediatric extension of the Janmahasatian 2005 semi-mechanistic model, derived from TBW, height, age, and sex)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Combined with WT to compute parameter-specific NFM. NFM reference values for a 70 kg TBW adult are 62 kg for CL (Ffat = 0.509), 59 kg for V1 and V2 (Ffat = 0.203), and 56.1 kg for Q (Ffat = 0); Ffat values are fixed a priori per McCune 2014 (Lawson 2022 Methods + Table 3).",
      source_name        = "FFM"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Subjects in the source dataset were assumed to have been born at term, so PAGE = postnatal age + 40 weeks. Drives the Hill maturation function on CL with TM50 = 45.6 weeks and Hill = 2.3 (both fixed; Lawson 2022 Methods Equation 3 and Table 3). PAGE is canonical in months and is converted to weeks inside the model.",
      source_name        = "PMA"
    )
  )

  population <- list(
    n_subjects     = 95,
    n_studies      = 1,
    n_centers      = 4,
    age_range      = "0.735-17.2 years",
    age_median     = "4.20 years",
    weight_range   = "7.77-83.3 kg",
    weight_median  = "17.0 kg",
    bmi_range      = "13.35-32.39 kg/m^2",
    bmi_median     = "18.2 kg/m^2",
    sex_female_pct = 48.4,
    disease_state  = "Pediatric hematopoietic stem cell transplant (HSCT) recipients (80 with malignant disease, 15 non-malignant) receiving once-daily IV busulfan as part of conditioning",
    dose_range     = "Initial 3.2-4.8 mg/kg/dose IV over ~3 h per Australian Busulfex product information weight bands (<9 kg: 4 mg/kg; 9 to <16 kg: 4.8 mg/kg; 16-23 kg: 4.4 mg/kg; >23-34 kg: 3.8 mg/kg; >34 kg: 3.2 mg/kg). Subsequent doses across 4 days were individualized to a cumulative AUC target near 90 mg.h/L.",
    regions        = "Australia (Brisbane, Perth, Sydney) and New Zealand (Auckland)",
    notes          = "Lawson 2022 Tables 1-2 and Subjects/data section; data collected 2016-2021 across 4 children's hospitals (80/95 from Queensland Children's Hospital). Conditioning regimens included Bu/Flu (25.3%), Bu/Flu/TT (25.3%), Bu/Mel (22.1%), Bu/Flu/Mel (20%), Bu/Cy (5.3%), and other (2.1%). Final dataset: 379 dosing days, 2491 plasma busulfan concentrations."
  )

  ini({
    # Structural parameters at adult reference NFM (70 kg TBW corresponds to NFM 62 kg for CL,
    # 59 kg for V1 and V2, 56.1 kg for Q; allometric exponents 0.75 on CL and Q, 1 on V1 and V2).
    lcl <- log(14.5);  label("Typical clearance at adult reference (CL, L/h)")                              # Lawson 2022 Table 3
    lvc <- log(40.6);  label("Typical central volume of distribution at adult reference (V1, L)")          # Lawson 2022 Table 3
    lq  <- log(1.92);  label("Typical intercompartmental clearance at adult reference (Q, L/h)")           # Lawson 2022 Table 3
    lvp <- log(3.57);  label("Typical peripheral volume of distribution at adult reference (V2, L)")       # Lawson 2022 Table 3

    # Time-associated CL decline within a treatment course (Lawson 2022 Equation 1, simplified
    # to gamma = 1; CL_t = CL * exp(dCLmax * t / (tm50_time + t)) where t is time since the
    # start of the first infusion).
    dCLmax    <- -0.198;  label("Maximum fractional change in CL over a treatment course (unitless)")      # Lawson 2022 Table 3
    tm50_time <-  50.6;   label("Time at which 50% of dCLmax is attained (h)")                             # Lawson 2022 Table 3

    # Fat-mass fractions for parameter-specific NFM (fixed a priori per McCune 2014).
    ffat_cl <- fixed(0.509);  label("Fraction of fat mass contributing to NFM on CL (unitless, fixed)")    # Lawson 2022 Table 3
    ffat_v1 <- fixed(0.203);  label("Fraction of fat mass contributing to NFM on V1 (unitless, fixed)")    # Lawson 2022 Table 3
    ffat_q  <- fixed(0.0);    label("Fraction of fat mass contributing to NFM on Q (unitless, fixed)")     # Lawson 2022 Table 3
    ffat_v2 <- fixed(0.203);  label("Fraction of fat mass contributing to NFM on V2 (unitless, fixed)")    # Lawson 2022 Table 3

    # Maturation parameters fixed a priori (Lawson 2022 Methods + Table 3).
    hill_mat <- fixed(2.3);   label("Hill coefficient for the CL maturation function (unitless, fixed)")   # Lawson 2022 Table 3 + Methods
    tm50_mat <- fixed(45.6);  label("Postmenstrual age at 50% of adult CL via maturation (weeks, fixed)")  # Lawson 2022 Table 3 + Methods

    # Inter-individual variability on CL and V1 (block).
    # Per Lawson 2022 Table 3 footnote: 'CV% are calculated as the Square root of variance
    # (OMEGA from NONMEM) x 100', so omega^2 = (CV/100)^2. Reported correlation coefficient
    # between CL and V1 is 0.0295, giving cov = 0.0295 * sqrt(omega2_lcl * omega2_lvc).
    etalcl + etalvc ~ c(0.021609,
                        0.001513, 0.121801)                                                                # Lawson 2022 Table 3 (IIV CL 14.7%, IIV V1 34.9%, corr 0.0295)

    # Combined residual error (Lawson 2022 Methods + Table 3; proportional RUV reported as
    # standard deviation per the table footnote).
    propSd <- 0.243;  label("Proportional residual error (fraction)")                                      # Lawson 2022 Table 3
    addSd  <- 0.030;  label("Additive residual error (mg/L)")                                              # Lawson 2022 Table 3
  })

  model({
    # Parameter-specific normal fat mass (Lawson 2022 Methods + Table 3 footnote): NFM is a
    # weighted sum of FFM and the fat compartment (TBW - FFM). With Ffat = 0 the size term
    # collapses to FFM; with Ffat = 1 it collapses to TBW.
    nfm_cl <- FFM + ffat_cl * (WT - FFM)
    nfm_v1 <- FFM + ffat_v1 * (WT - FFM)
    nfm_q  <- FFM + ffat_q  * (WT - FFM)
    nfm_v2 <- FFM + ffat_v2 * (WT - FFM)

    # Allometric size factors (Lawson 2022 Methods Equation 2). Reference NFMs correspond to
    # a 70 kg TBW adult; allometric exponents are 0.75 on CL and Q, 1 on V1 and V2.
    size_cl <- (nfm_cl / 62)^0.75
    size_v1 <- (nfm_v1 / 59)
    size_q  <- (nfm_q  / 56.1)^0.75
    size_v2 <- (nfm_v2 / 59)

    # CL maturation as a sigmoid Hill function of postmenstrual age in weeks (Lawson 2022
    # Methods Equation 3). PAGE is canonical in months; convert to weeks via the average
    # 365.25/12/7 ~= 4.348 weeks/month.
    pma_weeks <- PAGE * (365.25 / 12 / 7)
    maturation_cl <- 1 / (1 + (pma_weeks / tm50_mat)^(-hill_mat))

    # Time-associated CL multiplier over the treatment course (Lawson 2022 Equation 1 with
    # gamma = 1, Table 3). The simulation clock t must equal time since the first dose (i.e.,
    # dosing in the event table starts at t = 0); at t = 0 the multiplier is 1.
    time_factor_cl <- exp(dCLmax * t / (tm50_time + t))

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * size_cl * maturation_cl * time_factor_cl
    vc <- exp(lvc + etalvc) * size_v1
    q  <- exp(lq)           * size_q
    vp <- exp(lvp)          * size_v2

    # Two-compartment IV ODE in mass-balance form. The corresponding rate-of-change of central
    # concentration follows from CONC = central / vc and PERI = peripheral1 / vp; this matches
    # Lawson 2022 Table 3 final-model equations after correcting the printed sign of the
    # central-to-peripheral term.
    conc <- central     / vc
    peri <- peripheral1 / vp
    d/dt(central)     <- q * peri - q * conc - cl * conc
    d/dt(peripheral1) <- q * conc - q * peri

    Cc <- conc
    Cc ~ add(addSd) + prop(propSd)
  })
}
