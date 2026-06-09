SanchezPena_2005_enoxaparin <- function() {
  description <- paste(
    "One-compartment population pharmacokinetic model of anti-factor Xa",
    "activity after intravenous enoxaparin in 546 adults undergoing elective",
    "percutaneous coronary intervention (Sanchez-Pena 2005). The IV bolus is",
    "modelled as a brief zero-order input phase of duration T0 with linear",
    "elimination. Body weight is the only retained covariate, applied as",
    "estimated allometric exponents on clearance (0.9) and volume (0.7) with",
    "reference 75 kg. A fixed basal anti-Xa activity (0.0725 IU/mL) is added",
    "to the dose-driven concentration to account for the endogenous pre-dose",
    "background measured by the chromogenic anti-Xa assay. Doses must be",
    "entered in IU (1 mg enoxaparin = 100 IU anti-Xa); the typical 0.5 mg/kg",
    "clinical dose corresponds to 3830 IU for a 76 kg patient."
  )
  reference <- paste(
    "Sanchez-Pena P, Hulot JS, Urien S, Ankri A, Collet JP, Choussat R,",
    "Lechat P, Montalescot G. Anti-factor Xa kinetics after intravenous",
    "enoxaparin in patients undergoing percutaneous coronary intervention:",
    "a population model analysis. Br J Clin Pharmacol. 2005;60(4):364-373.",
    "doi:10.1111/j.1365-2125.2005.02452.x"
  )
  vignette <- "SanchezPena_2005_enoxaparin"
  units    <- list(
    time          = "hour",
    dosing        = "IU",
    concentration = "IU/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Estimated allometric scaling exponents on CL (0.9) and V (0.7);",
        "reference weight 75 kg per the covariate equation in",
        "Sanchez-Pena 2005 Results (immediately above Table 2):",
        "CL = TV_CL * (BW/75)^0.9 and V = TV_V * (BW/75)^0.7.",
        "Time-fixed (baseline only) in the source data set."
      ),
      source_name        = "BW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age (years)",
      units       = "year",
      type        = "continuous",
      notes       = paste(
        "Screened in the covariate sub-model (Sanchez-Pena 2005 Methods,",
        "Population pharmacokinetic modelling section). Not retained in the",
        "final model: no significant improvement in the population PK model",
        "when included."
      )
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      reference_category = "0 (male)",
      notes       = paste(
        "Screened (paper text: 'gender'). Not retained in the final model.",
        "Cohort was 79% male / 21% female (Table 1)."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened. The Results note CL was 'influenced to a lesser degree by",
        "serum creatinine', but the effect did not survive the OFV-7-unit",
        "backward elimination step. Not retained in the final model."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (raw Cockcroft-Gault, NOT BSA-normalized)",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Screened. The Results note CL was 'influenced to a lesser degree by",
        "creatinine clearance', but the effect did not survive the OFV-7-unit",
        "backward elimination step. Discussion attributes the lack of CrCl",
        "retention to the short anticoagulation window after a single IV bolus,",
        "for which V (not CL) is the main determinant of plasma concentration.",
        "Note this is raw Cockcroft-Gault in mL/min, not the canonical",
        "BSA-normalised CRCL (mL/min/1.73 m^2)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 546L,
    n_studies      = 1L,
    n_observations = 1978L,
    age_range      = "21-93 years",
    age_median     = "63 years (SD 11)",
    weight_range   = "35-153 kg",
    weight_median  = "76 kg (SD 15)",
    sex_female_pct = 21,
    race_ethnicity = "Not reported (single-centre Paris cohort)",
    disease_state  = paste(
      "Adults referred for elective PCI of coronary or vein-graft stenosis",
      "> 70%. Excluded: primary PCI for ST-elevation MI, LMWH or UFH within",
      "the prior 48 h, or GPIIb/IIIa antagonist within the prior 2 weeks.",
      "Renal function distribution: CrCl < 30 mL/min in 4%, 31-59 mL/min in",
      "33%, >= 60 mL/min in 62%. 15% aged > 75 years."
    ),
    dose_range     = paste(
      "0.5 mg/kg single IV bolus (mean dose 38 +/- 7 mg = 3830 +/- 730 IU)",
      "immediately before PCI. Simulations in the source paper also explored",
      "0.75 mg/kg and 1 mg/kg single doses, and second-bolus regimens."
    ),
    sampling       = paste(
      "Five samples per patient: pre-bolus, 10 min post-bolus (start of PCI),",
      "end of PCI (mean 45 min), 3 h post-PCI, and the morning after PCI;",
      "1978 quantifiable anti-Xa concentrations after exclusions."
    ),
    regions        = "Paris, France (Pitie-Salpetriere University Hospital).",
    notes          = paste(
      "Demographics from Sanchez-Pena 2005 Table 1. 556 patients enrolled;",
      "10 (1.8%) excluded for probable misadministration leaving 546 in the",
      "final analysis. Concomitant GPIIb/IIIa inhibitors: eptifibatide (146",
      "patients), abciximab (29 patients). All patients received aspirin",
      "(500 mg IV loading then 75 mg/day PO)."
    )
  )

  ini({
    # Structural PK parameters (Sanchez-Pena 2005 Table 2, "Final model"
    # column). Typical values are reported for a patient with BW = 75 kg
    # (the reference weight in the covariate equation immediately preceding
    # Table 2: CL = TV_CL * (BW/75)^0.9, V = TV_V * (BW/75)^0.7).
    lcl <- log(1.20)
    label("Clearance at BW = 75 kg (L/h)")  # Sanchez-Pena 2005 Table 2: TV.CL = 1.20 L/h (SE 0.03)
    lvc <- log(2.9)
    label("Central volume at BW = 75 kg (L)")  # Sanchez-Pena 2005 Table 2: TV.V = 2.9 L (SE 0.1)

    # Zero-order IV input duration T0 (paper's nomenclature). Untransformed
    # because the paper modelled inter-individual variability on T0 as
    # additive rather than exponential (Sanchez-Pena 2005 Results, page on
    # population PK: "ISV for zero-order input (T0) and residual variability
    # were best described by an additive error model").
    tdur <- 0.25
    label("Zero-order IV input duration T0 (h)")  # Sanchez-Pena 2005 Table 2: T0 = 0.25 h (SE 0.01)

    # Allometric exponents on body weight. Estimated (SE reported for both),
    # NOT fixed. Reference weight 75 kg.
    e_wt_cl <- 0.9
    label("Allometric exponent of BW on CL (unitless)")  # Sanchez-Pena 2005 Table 2: BW effect on CL = 0.9 (SE 0.1)
    e_wt_vc <- 0.7
    label("Allometric exponent of BW on V (unitless)")  # Sanchez-Pena 2005 Table 2: BW effect on V = 0.7 (SE 0.1)

    # Fixed basal anti-Xa activity. The chromogenic assay returns a non-zero
    # pre-dose reading; the paper added this constant to the dose-driven
    # prediction and held it fixed during estimation.
    bl_antiXa <- fixed(0.0725)
    label("Basal anti-Xa activity (IU/mL); fixed")  # Sanchez-Pena 2005 Table 2: Basal anti-Xa = 0.0725 IU/mL (FIXED)

    # Inter-individual variability. CL and V use an exponential (log-normal)
    # model so the etas attach to the log-transformed structural parameters;
    # the CV percentages from Table 2 convert to the internal variance scale
    # via omega^2 = log(CV^2 + 1). T0 uses an additive model so its eta is
    # an additive shift with the reported SD (0.06 h); variance = 0.06^2.
    etalcl  ~ log(0.33^2 + 1)  # Sanchez-Pena 2005 Table 2: ISV(CL) = 33% CV; var = log(0.33^2 + 1) = 0.10336
    etalvc  ~ log(0.30^2 + 1)  # Sanchez-Pena 2005 Table 2: ISV(V)  = 30% CV; var = log(0.30^2 + 1) = 0.08618
    etatdur ~ 0.06^2           # Sanchez-Pena 2005 Table 2: ISV(T0) = 0.06 h (additive SD); var = 0.0036

    # Residual error: additive on linear (anti-Xa activity) scale.
    addSd <- 0.09
    label("Additive residual SD (IU/mL)")  # Sanchez-Pena 2005 Table 2: residual variability = 0.09 IU/mL (SE 0.02)
  })

  model({
    # Individual PK parameters with the BW power-model covariate. Reference
    # weight 75 kg, with separately estimated exponents on CL and V.
    cl     <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    vc     <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    tdur_i <- tdur + etatdur

    # First-order elimination rate constant
    kel <- cl / vc

    # ODE: single central compartment receiving zero-order input
    d/dt(central) <- -kel * central

    # Engage the rxode2 model-defined zero-order input duration. Dose
    # records must specify rate = -2 to trigger duration-based delivery.
    dur(central) <- tdur_i

    # Observation. Dose is entered in IU and V is in L, so central / vc has
    # units IU/L; divide by 1000 to obtain IU/mL (the assay unit reported
    # throughout the paper) and add the fixed endogenous basal level.
    Cc <- central / vc / 1000 + bl_antiXa

    Cc ~ add(addSd)
  })
}
