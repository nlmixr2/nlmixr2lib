Hornik_2025_furosemide <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for furosemide in adults with chronic heart failure and volume overload, allometrically scaled to adolescents (Hornik 2025)"
  reference <- "Hornik CP, Foote HP, Kendig E, Mohr J. An Adult Population Pharmacokinetic Model to Simulate Subcutaneous Administration of a Fixed Dose of Furosemide in Adolescents with Heart Failure and Volume Overload. Clin Pharmacokinet. 2025;64:899-908. doi:10.1007/s40262-025-01515-2"
  vignette <- "Hornik_2025_furosemide"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "furosemide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "furosemide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "furosemide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Renal function. Hornik 2025 Sect. 2.4 calls this 'eGFR' and states it was computed with the Cockcroft-Gault equation; the ETA-versus-covariate panels in the electronic supplementary material (Figures S1-S5) label the same axis 'Creatinine Clearance, mL/min/1.73 m^2'. The median value in the adult analysis population, and the reference used by Eq. 15, is 86.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL normalized to the adult population median of 86 (Eq. 15, Table 1). The paper labels this covariate inconsistently between the main text (Cockcroft-Gault eGFR, mL/min) and the supplement (BSA-normalized creatinine clearance, mL/min/1.73 m^2); the register's CRCL units are used here. During covariate screening the authors also evaluated a saturable form in which eGFR was capped at 120 mL/min (Sect. 2.4); that form was not retained.",
      source_name        = "eGFR"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling to the 70 kg standardized adult (Eqs. 9-13). Weight was NOT retained as a covariate in the adult fit itself; the allometric terms were applied post hoc to extrapolate the adult model to adolescents, so all size terms evaluate to 1 at WT = 70 kg and the model then reduces exactly to the published adult model (Eqs. 14-19).",
      source_name        = "WT"
    )
  )

  # Screened in the adult covariate search (Sect. 2.4, Table S1, Figures S1-S5)
  # but not retained in the final model, so they are documented rather than
  # implemented. AGE on Vp was selected in the second forward step (dOFV
  # -6.229; it had scored -6.313 in the first step, which eGFR on CL won at
  # -15.546), then removed in backward elimination, where dropping it cost
  # only +6.383 (Sect. 3.2; Table S1).
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Selected on peripheral volume (Vp) in the second forward step but removed in backward elimination (Sect. 3.2; Table S1). Normalized to the population median of 69 years when tested."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened graphically against all five ETAs (Figures S1-S5); not carried into forward selection. Superseded by CRCL, which is derived from it."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened graphically (Figures S1-S5); not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a dichotomous covariate via Eq. 8 (Figures S1-S5); not retained. The analysis population was 93.3% male (14 of 15), so the female stratum was a single subject."
    ),
    NYHA3 = list(
      description = "New York Heart Association heart failure class III indicator (0 = class II)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a dichotomous covariate via Eq. 8 (Figures S1-S5); not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 15,
    n_studies      = 1,
    age_range      = "52-83 years",
    age_median     = "69 years",
    weight_range   = "approximately 60-120 kg (not tabulated in the paper; read from the weight axis of ETA-versus-covariate Figures S1-S5)",
    weight_median  = "not reported",
    sex_female_pct = 6.7,
    race_ethnicity = "not reported",
    disease_state  = "chronic heart failure, NYHA class II/III, with chronic volume overload requiring oral furosemide >= 40 mg/day for >= 30 days before entry",
    renal_function = "median eGFR 86 mL/min (Cockcroft-Gault); observed range approximately 60-140 mL/min/1.73 m^2 from Figures S1-S5",
    dose_range     = "80 mg subcutaneous (30 mg over the first 60 min then 12.5 mg/h for 4 h) and 80 mg intravenous (40 mg over 2 min, repeated 2 h later), in a randomized crossover with a 7-day washout",
    regions        = "USA (single site)",
    notes          = paste(
      "Adult analysis population of NCT02329834 (Sect. 2.1, Sect. 3.1).",
      "Sixteen adults received furosemide; one was excluded for high pre-dose concentrations.",
      "The paper reports no baseline demographics table, so only age and sex are quoted directly.",
      "The model was subsequently applied to virtual adolescents aged 12-17 years with body weight >= 42.5 kg",
      "in three weight bands (42.5-50.0, >50-60, >60-70 kg) via the allometric terms in Eqs. 9-13;",
      "no maturation function was used because renal maturation is complete in this age range (Sect. 2.6)."
    )
  )

  ini({
    # Structural parameters. Typical values refer to a 70 kg adult with
    # CRCL 86 mL/min/1.73 m^2 (Sect. 2.6 definition of Kastd/CLstd/Vstd).
    lka <- log(1.30); label("Subcutaneous first-order absorption rate constant (Ka, 1/h)")             # Table 1, Ka; Eqs. 9, 14
    lcl <- log(6.51); label("Clearance in a 70 kg adult at the median CRCL of 86 (CL, L/h)")           # Table 1, CL; Eqs. 10, 15
    lvc <- log(5.37); label("Central volume of distribution in a 70 kg adult (Vc, L)")                 # Table 1, Vc; Eqs. 11, 16
    lq  <- log(3.74); label("Intercompartmental clearance in a 70 kg adult (Q, L/h)")                  # Table 1, Q; Eqs. 12, 17
    lvp <- log(5.95); label("Peripheral volume of distribution in a 70 kg adult (Vp, L)")              # Table 1, Vp; Eqs. 13, 18

    # Sect. 3.2: "bioavailability following subcutaneous furosemide
    # administration was estimated at 0.96". The paper estimated it on the
    # logit scale, F = exp(theta(6)) / (1 + exp(theta(6))) (Eq. 19); only the
    # back-transformed 0.96 is published, and Table 1 omits theta(6) entirely.
    # See the vignette Errata: the paper's own Table 2 exposures are only
    # reproducible with F = 1.00, so this value is internally inconsistent with
    # the published simulations. The published estimate is used here.
    lfdepot <- log(0.96); label("Subcutaneous bioavailability (fraction)")                             # Sect. 3.2 text; Eq. 19

    # Covariate effect
    e_crcl_cl <- 0.74; label("Power exponent on (CRCL / 86) for CL (unitless)")                        # Table 1, "eGFR on CL"; Eq. 15

    # Allometric size terms. Sect. 2.6 and Sect. 4 state the CL exponent was
    # held at the conventional 0.75 and the volume exponent at 1 rather than
    # estimated; neither is reported with an RSE or bootstrap interval.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT / 70) for CL and Q (unitless)")       # Eqs. 10, 12
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT / 70) for Vc and Vp (unitless)")      # Eqs. 11, 13

    # Inter-individual variability, exponential (Eq. 1). Table 1 reports CV%;
    # converted to the variance scale as omega^2 = log(CV^2 + 1). Table 1
    # reports diagonal elements only -- no off-diagonal was reported for the
    # final model, so the etas are independent here even though Eq. 2 defines
    # how a CL-V correlation would have been computed.
    etalka ~ 0.046014  # 21.7% CV (shrinkage 10%), Table 1
    etalcl ~ 0.033652  # 18.5% CV (shrinkage 0%),  Table 1
    etalvc ~ 0.095755  # 31.7% CV (shrinkage 0%),  Table 1
    etalq  ~ 0.103964  # 33.1% CV (shrinkage 6%),  Table 1
    etalvp ~ 0.035098  # 18.9% CV (shrinkage 3%),  Table 1

    # Combined residual error (Eq. 5)
    propSd <- 0.126; label("Proportional residual error (fraction)")                                   # Table 1, 12.6%
    addSd  <- 0.100; label("Additive residual error (ug/mL)")                                          # Table 1, 100 ng/mL = 0.100 ug/mL
  })
  model({
    # Allometric scaling to the 70 kg standardized adult (Eqs. 9-13). Both
    # terms equal 1 at WT = 70 kg, where the model reduces to the published
    # adult model of Eqs. 14-19.
    allom_cl <- (WT / 70)^e_wt_cl_q
    allom_v  <- (WT / 70)^e_wt_vc_vp

    ka <- exp(lka + etalka)                                          # Eq. 14; Ka is not weight-scaled (Eq. 9)
    cl <- exp(lcl + etalcl) * (CRCL / 86)^e_crcl_cl * allom_cl       # Eqs. 15, 10
    vc <- exp(lvc + etalvc) * allom_v                                # Eqs. 16, 11
    q  <- exp(lq  + etalq)  * allom_cl                               # Eqs. 17, 12
    vp <- exp(lvp + etalvp) * allom_v                                # Eqs. 18, 13

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability applies to the subcutaneous depot only; intravenous
    # doses are administered directly into central and are unaffected.
    f(depot) <- exp(lfdepot)                                         # Eq. 19

    # Dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
