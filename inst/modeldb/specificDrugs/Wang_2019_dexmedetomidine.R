Wang_2019_dexmedetomidine <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "intranasally administered dexmedetomidine in Chinese children aged 4-10",
    "years undergoing surgery (Wang 2019; n = 13, 145 arterial plasma",
    "concentrations after a single 1 ug/kg intranasal dose given as drops).",
    "No intravenous arm was studied, so all disposition parameters are",
    "apparent (CL/F, V1/F, Q2/F, V2/F). Body weight enters as fixed",
    "allometric scaling standardised to 70 kg (exponent 0.75 on clearances,",
    "1 on volumes). Inter-individual variability on Ka, CL/F, V1/F and Q2/F",
    "(V2/F variance fixed to zero); additive residual error only."
  )
  reference <- paste(
    "Wang CY, Ihmsen H, Hu ZY, Chen J, Ye XF, Chen F, Lu Y, Schuttler J,",
    "Lian QQ, Liu HC. Pharmacokinetics of Intranasally Administered",
    "Dexmedetomidine in Chinese Children. Front Pharmacol. 2019;10:756.",
    "doi:10.3389/fphar.2019.00756"
  )
  vignette <- "Wang_2019_dexmedetomidine"
  units <- list(time = "min", dosing = "ug", concentration = "pg/mL")

  compartmentData <- list(
    depot = list(
      analyte = "dexmedetomidine",
      units = "ug",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "dexmedetomidine",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "dexmedetomidine",
      units = "ug",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of all disposition parameters, Methods",
        "'Pharmacokinetic Modeling': Pi = PSTD * (BWi/70)^PBW, with PBW",
        "fixed at 0.75 for clearances (CL/F, Q2/F) and 1 for volumes",
        "(V1/F, V2/F). Estimating the exponents did not improve the fit",
        "(Results, OFV 1286.6 vs 1287.1). Ka is not weight-scaled. Study",
        "range 14-36 kg (Table 1, mean 24.8 kg)."
      ),
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened after weight scaling and not retained (Results: 'there",
        "were no further significant effects of age, sex, and laboratory",
        "parameters detected'). Table 1: 6.5 +/- 1.7 years (4-10)."
      ),
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. Table 1: 8 female / 5 male.",
      source_name = "sex"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened and not retained. Table 1: 37.3 +/- 7.3 umol/L (23.8-47.0).",
      source_name = "Cr"
    ),
    BUN_CREAT_RATIO = list(
      description = "Ratio of blood urea nitrogen to serum creatinine",
      units = "(ratio)",
      type = "continuous",
      notes = "Screened and not retained. Table 1: 0.11 +/- 0.03 (0.04-0.15).",
      source_name = "BUN/Cr"
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault, not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened and not retained. Table 1: 99.2 +/- 19.0 mL/min (78.2-150).",
      source_name = "CrCL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 13,
    n_studies = 1,
    n_observations = 145,
    age_range = "4-10 years",
    age_mean = "6.5 years (SD 1.7)",
    weight_range = "14-36 kg",
    weight_mean = "24.8 kg (SD 5.8)",
    sex_female_pct = 61.5,
    race_ethnicity = "Chinese (100%)",
    disease_state = paste(
      "Children under general anaesthesia (propofol/fentanyl induction,",
      "sevoflurane maintenance) for orthopaedic, lower abdominal, urologic",
      "or plastic surgery"
    ),
    dose_range = paste(
      "1 ug/kg dexmedetomidine (100 ug/mL, undiluted) intranasally as",
      "drops from a 1-mL syringe, single dose 30 min before the end of",
      "surgery"
    ),
    regions = "China (Wenzhou; single centre; NCT03337581)",
    notes = paste(
      "Table 1 demographics. 16 enrolled, 1 withdrew, 2 excluded for blocked",
      "arterial catheters. Arterial samples to 180 min (shorter in 5",
      "subjects); no concentrations below the 0.01 ng/mL LLOQ. NONMEM 7.4.3",
      "FOCE-I."
    )
  )

  ini({
    # Table 2 'Final population pharmacokinetic parameters', typical values for a 70 kg subject.
    lka <- log(0.0378); label("Absorption rate constant Ka (1/min)")                         # Table 2 'Ka (1/min)' = 0.0378 (RSE 10.2%)
    lcl <- log(0.319);  label("Apparent clearance CL/F (L/min) at WT = 70 kg")                # Table 2 'CL/F (L/min/70 kg)' = 0.319 (RSE 6.6%)
    lvc <- log(34.2);   label("Apparent central volume V1/F (L) at WT = 70 kg")               # Table 2 'V1/F (L/70 kg)' = 34.2 (RSE 14.4%)
    lq  <- log(10.0);   label("Apparent intercompartmental clearance Q2/F (L/min) at WT = 70 kg") # Table 2 'Q2/F (L/min/70 kg)' = 10.0 (RSE 22.2%)
    lvp <- log(34.9);   label("Apparent peripheral volume V2/F (L) at WT = 70 kg")            # Table 2 'V2/F (L/70 kg)' = 34.9 (RSE 6.5%)

    # Methods 'Pharmacokinetic Modeling' and Results: allometric power model with
    # PBW = 0.75 for clearances and 1 for volumes, not estimated.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of body weight on CL/F and Q2/F (unitless)") # Methods/Results: 'exponent of 1 for all volumes and of 0.75 for all clearances'
    e_wt_vc_vp  <- fixed(1);    label("Allometric exponent of body weight on V1/F and V2/F (unitless)") # Methods/Results: 'exponent of 1 for all volumes and of 0.75 for all clearances'

    # IIV: Methods Pi = PTV * exp(eta_i); Table 2 reports omega^2 (variances).
    # The V2/F variance was < 0.0001 and fixed to zero (Table 2 '0 (fixed)'),
    # so no eta is carried on lvp.
    etalka ~ 0.0112  # Table 2 'Interindividual variability omega2, Ka' = 0.0112 (RSE 50.0%)
    etalcl ~ 0.0207  # Table 2 'Interindividual variability omega2, CL/F' = 0.0207 (RSE 38.7%)
    etalvc ~ 0.00569 # Table 2 'Interindividual variability omega2, V1/F' = 0.00569 (RSE 82.7%)
    etalq  ~ 1.29    # Table 2 'Interindividual variability omega2, Q2/F' = 1.29 (RSE 50.6%)

    # Residual: the proportional term was removed (sigma2_prop < 0.0001), leaving
    # additive error only. Table 2 sigma2_add = 1752 (pg/mL)^2 -> SD = sqrt(1752).
    addSd <- 41.857; label("Additive residual error SD (pg/mL)") # Table 2 'sigma2 add' = 1752 (RSE 14.6%); sqrt(1752) = 41.857
  })

  model({
    # Allometric scaling standardised to 70 kg (Methods 'Pharmacokinetic Modeling').
    allom_cl <- (WT / 70)^e_wt_cl_q
    allom_v <- (WT / 70)^e_wt_vc_vp

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * allom_cl
    vc <- exp(lvc + etalvc) * allom_v
    q <- exp(lq + etalq) * allom_cl
    vp <- exp(lvp) * allom_v

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intranasal dose into depot; parameters are apparent (/F), so no
    # bioavailability term is applied.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in ug and volume in L give ug/L = ng/mL; the paper reports pg/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
