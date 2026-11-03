vanderVeer_2025_amikacin <- function() {
  description <- "Two-compartment IV-infusion population PK model for amikacin in neonates (van der Veer 2025 'optimized' model: the De Cock 2012 structure re-estimated on pooled Amsterdam UMC and University Hospitals Leuven NICU data). Clearance scales as a power of birth weight (reference 1685 g) times a linear function of postnatal age (1 + n * PNA / 5 days); central volume scales as a power of current weight (reference 1700 g); peripheral volume equals central volume and intercompartmental clearance is proportional to clearance (Q = 0.468 x CL). IIV on CL, and a 30% IIV on V1 set by the authors; combined additive and proportional residual error."
  reference <- "van der Veer MAA, Smits A, de Haan TR, Franken LGW, van Kaam AH, Hodiamont CJ, Bijleveld YA, Allegaert K, Mathot RAA. Amikacin dosing in neonates: evaluation of target attainment using a simplified and complex pharmacokinetic model-derived dosing regimen in clinical practice. Antimicrob Agents Chemother. 2025;69(4):e01118-24. doi:10.1128/aac.01118-24"
  vignette <- "vanderVeer_2025_amikacin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT_BIRTH = list(
      description = "Body weight at birth",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL centred on the pooled-dataset median birth weight of 1685 g (Supplement Table S5 footnote b). The paper reports birth weight (bBW) in grams; encoded here as WT_BIRTH / 1.685 with WT_BIRTH in kg.",
      source_name = "bBW"
    ),
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the central volume V1 (and, through V2 = V1, on the peripheral volume) centred on the pooled-dataset median current weight of 1700 g (Supplement Table S5 footnote b). The paper reports current weight (cBW) in grams; encoded here as WT / 1.700 with WT in kg.",
      source_name = "cBW"
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL, (1 + n * PNA / 5) with PNA in days centred on the pooled-dataset median of 5 days (Supplement Table S5 footnote b). model() converts the canonical months to days with PNA * 30.4375 before forming the ratio.",
      source_name = "PNA"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_IBUPROFEN = list(
      description = "Concomitant ibuprofen at the time of amikacin administration (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Carried in the De Cock 2012 model as a multiplicative factor o = 0.838 on CL. In the van der Veer 2025 optimized model the factor equals 1 'given low ibuprofen concurrent use' (Supplement Table S5 footnote c; 2% and 5% of the two cohorts, Table 3), so ibuprofen has no effect on CL and is not referenced in model().",
      source_name = "ibuprofen"
    )
  )

  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 710L,
    n_studies = 2L,
    n_concentrations = 1412L,
    age_range = "PNA 0-31 days; PMA 24-54 weeks (pooled Table 3 ranges)",
    age_median = "PNA 5 days (pooled model-building dataset, Supplement Table S5 footnote b)",
    ga_range = "24-42.1 weeks (Amsterdam median 28.7, Leuven median 34)",
    weight_range = "Current weight 420-5500 g; birth weight 420-5200 g (pooled Table 3 ranges)",
    weight_median = "Current weight 1700 g; birth weight 1685 g (pooled model-building dataset, Supplement Table S5 footnote b)",
    sex_female_pct = 42.9,
    disease_state = "Neonates in the NICU treated with amikacin for suspected or proven early- or late-onset sepsis",
    dose_range = "IV infusion over 20-30 min; simplified regimen 15 mg/kg every 36 h (GA < 30 weeks) or every 24 h (GA >= 30 weeks); complex regimen 15-20 mg/kg every 20-48 h by postnatal age and current weight (Tables 1 and 2)",
    regions = "Netherlands (Amsterdam UMC) and Belgium (University Hospitals Leuven)",
    renal_function = "Serum creatinine median 42 umol/L (Amsterdam) and 88 umol/L (Leuven)",
    co_medication = "Ibuprofen in 2% (Amsterdam) and 5% (Leuven) of neonates",
    notes = "Model-building dataset = a random 75% of the 945 neonates pooled from Amsterdam UMC (n = 366, simplified regimen, 2018-2024) and University Hospitals Leuven (n = 579, complex regimen, 2011-2012); the remaining 25% formed the validation dataset. Demographics are reported per centre in Table 3 (sex from 216/366 and 324/579 males); the model-building subset is described only by its covariate medians. NONMEM; routine TDM samples."
  )

  ini({
    # Structural parameters: Supplement Table S5, column 'Model-building
    # dataset (optimized PK model)'. Covariate centring values: birth weight
    # 1685 g, current weight 1700 g, PNA 5 days (Table S5 footnote b).
    lcl <- log(0.044)
    label("Clearance at birth (PNA 0) for birth weight 1685 g (L/h)") # Table S5: CLp = 0.044 (RSE 3.2%); the PNA term (1 + n x PNA/5) equals 1 at PNA 0
    lvc <- log(0.794)
    label("Central volume V1 at current weight 1700 g (L)") # Table S5: Vp = 0.794 (RSE 1.5%)
    qclratio <- 0.468
    label("Intercompartmental clearance as a fraction of CL, Q / CL (unitless)") # Table S5: Q = r x CL, r = 0.468 (RSE 7.8%)

    # Covariate effects (Table S5 equations:
    # CL = CLp x (bBW/median)^m x (1 + n x [PNA/median]);
    # V1 = Vp x (cBW/median)^p; V2 = V1).
    e_wt_birth_cl <- 1.34
    label("Power exponent of birth weight on CL (unitless)") # Table S5: m = 1.34 (RSE 2.4%)
    e_pna_cl <- 0.597
    label("Linear coefficient of PNA / 5 days on CL (unitless)") # Table S5: n = 0.597 (RSE 7.8%)
    e_wt_vc <- 0.903
    label("Power exponent of current weight on V1 (unitless)") # Table S5: p = 0.903 (RSE 2.6%)

    # IIV: Table S5 reports omega^2 (the 'IIV on V1' entry 0.09 is the
    # text's 30% fixed IIV, 0.3^2).
    etalcl ~ 0.0856 # Table S5: IIV on CL = 0.0856 (RSE 11.3%)
    etalvc ~ fixed(0.09) # Table S5: IIV on V1 = 0.09, set to a 30% CV by the authors (Results; Supplementary Data 4)

    # Residual error: Table S5 reports NONMEM sigma^2; SD = sqrt(variance).
    propSd <- 0.16703
    label("Proportional residual error (fraction)") # Table S5: proportional error variance 0.0279 (RSE 13.7%) -> sqrt = 0.167
    addSd <- 0.61482
    label("Additive residual error (mg/L)") # Table S5: additive error variance 0.378 (RSE 17.5%) -> sqrt = 0.615
  })
  model({
    # Postnatal age in days (canonical PNA is in months) and weights in kg
    # against the gram medians of Table S5 footnote b.
    pna_days <- PNA * 30.4375

    cl <- exp(lcl + etalcl) *
      (WT_BIRTH / 1.685)^e_wt_birth_cl *
      (1 + e_pna_cl * (pna_days / 5))
    vc <- exp(lvc + etalvc) * (WT / 1.700)^e_wt_vc
    vp <- vc
    q <- qclratio * cl

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
