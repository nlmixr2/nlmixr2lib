Ngara_2020_atazanavir <- function() {
  description <- paste(
    "One-compartment plasma PK model with first-order absorption plus a",
    "non-eliminating hair compartment for atazanavir given as ritonavir-boosted",
    "atazanavir (ATV/r 300/100 mg once daily) to HIV-infected Zimbabwean",
    "adolescents (Ngara 2020). Only hair concentrations were measured: the",
    "plasma parameters (ka 0.44 1/h, CL/F 10 L/h, Vc 63.4 L, allometric",
    "exponents 0.75 / 1 on body weight) are fixed from published adult and",
    "paediatric atazanavir models, and a fraction Frac (0.16) of the drug",
    "eliminated from plasma is deposited in hair (dA3/dt = Frac * k20 * A2),",
    "read out as Chair = A3 / Vh with Vh fixed to 1. Frac carries all",
    "covariate effects (PsN linear form): 30% lower at the enrolment visit",
    "than at the 90-day follow-up, +2% per unit of self-reported adherence",
    "(visual analogue scale), +54% / -21% for thin / overweight",
    "BMI-for-age versus normal, and +53% / +12% / -54% for parent /",
    "aunt-or-uncle / sibling primary caregiver versus grandparent. Combined",
    "additive plus proportional residual error. The printed parameters",
    "predict hair levels about 18-fold above the paper's own observations",
    "(see the vignette). The companion ritonavir model from the same paper",
    "is Ngara_2020_ritonavir."
  )
  reference <- paste(
    "Ngara B, Zvada S, Chawana TD, Stray-Pedersen B, Nhachi CFB,",
    "Rusakaniko S (2020). A population pharmacokinetic model is beneficial",
    "in quantifying hair concentrations of ritonavir-boosted atazanavir: a",
    "study of HIV-infected Zimbabwean adolescents.",
    "BMC Pharmacol Toxicol 21:58. doi:10.1186/s40360-020-00437-y",
    sep = " "
  )
  vignette <- "Ngara_2020_atazanavir_ritonavir_hair"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = paste(
      "ng/mg hair (Chair, the fitted output; numerically the hair amount in",
      "mg divided by Vh in L, as the paper's Equation 4 writes it); mg/L",
      "(plasma Cc, not observed)"
    )
  )

  # Issue #482. `hair` accumulates drug and never loses it: Ngara 2020
  # Equation 3 has no outflow term.
  compartmentData <- list(
    depot = list(
      analyte = "atazanavir",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "atazanavir",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    hair = list(
      analyte = "atazanavir",
      units = "mg",
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
        "Allometric scaling of CL/F (exponent 0.75) and Vc (exponent 1), both",
        "fixed (Ngara 2020 Results, 'Atazanavir model'). The paper does not",
        "print the reference weight or the cohort weights; 70 kg is assumed",
        "here. The Discussion says some plasma values were 'adjusted by the",
        "median body weight observed in our study', which is not printed.",
        "See the vignette Errata."
      ),
      source_name = "WT"
    ),
    OCC = list(
      description = "Hair-sampling occasion (1 = enrolment, 2 = 90-day follow-up)",
      units = "(count)",
      type = "categorical",
      reference_category = "2 (follow-up visit)",
      notes = paste(
        "Two occasions per participant. Used as a fixed effect on Frac, not",
        "for inter-occasion variability: Frac is 30% lower at enrolment",
        "(OCC = 1) than at follow-up (Ngara 2020 Table 3,",
        "'Occasion (Enrolment)_Frac' = -0.30, follow-up the reference).",
        "The authors attribute the rise at follow-up to better adherence while",
        "under study rather than to a pharmacokinetic change."
      ),
      source_name = "Occasion"
    ),
    ADHERENCE_VAS = list(
      description = "Self-reported adherence on a 0-100% visual analogue scale",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on Frac, +0.02 per unit (Ngara 2020 Table 3,",
        "'Adherence score_Frac'). The centring value is not printed; the",
        "effect is centred on 84, the rounded cohort mean (84.2, SD 18.1,",
        "range 30-100; Table 1). With this centring the factor is zero at",
        "a VAS of 34 and NEGATIVE below it, so the model as printed predicts a",
        "negative hair concentration for the least adherent participants; it",
        "is not clamped. See the vignette Errata."
      ),
      source_name = "Adherence VAS"
    ),
    BMI_FOR_AGE_THIN = list(
      description = "WHO BMI-for-age category: thinness (1 = thin, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal BMI-for-age, when BMI_FOR_AGE_OVERWEIGHT is also 0)",
      notes = "Ngara 2020 Table 3 'Body Mass Index-for-age (Thin)_Frac' = +0.54; 14 of 46 categorised participants (30%, Table 1).",
      source_name = "Body Mass Index-for-age"
    ),
    BMI_FOR_AGE_OVERWEIGHT = list(
      description = "WHO BMI-for-age category: overweight (1 = overweight, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal BMI-for-age, when BMI_FOR_AGE_THIN is also 0)",
      notes = "Ngara 2020 Table 3 'Body Mass Index-for-age (Over-weight)_Frac' = -0.21; 7 of 46 categorised participants (15%, Table 1).",
      source_name = "Body Mass Index-for-age"
    ),
    CAREGIVER_PARENT = list(
      description = "Primary caregiver is a parent (1 = yes, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (grandparent caregiver, when the other CAREGIVER_ indicators are also 0)",
      notes = "Ngara 2020 Table 3 'Guardian (Parent)_Frac' = +0.53; 10 of 50 participants (20%, Table 1).",
      source_name = "Guardian"
    ),
    CAREGIVER_AUNT_UNCLE = list(
      description = "Primary caregiver is an aunt or uncle (1 = yes, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (grandparent caregiver, when the other CAREGIVER_ indicators are also 0)",
      notes = "Ngara 2020 Table 3 'Guardian (Uncle/Aunt)_Frac' = +0.12; 15 of 50 participants (30%, Table 1).",
      source_name = "Guardian"
    ),
    CAREGIVER_SIBLING = list(
      description = "Primary caregiver is a sibling (1 = yes, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (grandparent caregiver, when the other CAREGIVER_ indicators are also 0)",
      notes = "Ngara 2020 Table 3 'Guardian (Sibling)_Frac' = -0.54; 5 of 50 participants (10%, Table 1), siblings younger than 19 years.",
      source_name = "Guardian"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on Frac by PsN SCM (Ngara 2020 Methods: all Table 1 covariates tested) and not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Frac by PsN SCM and not retained. Education level, WHO disease stage, ART regimen and hair length were also screened and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = "10-18 years (Table 1 range 11-18)",
    age_mean = "15.8 years (SD 1.8)",
    weight_range = "not reported",
    sex_female_pct = 54,
    race_ethnicity = "Zimbabwean (Black African); not tabulated",
    disease_state = paste(
      "HIV-infected adolescents on ritonavir-boosted-atazanavir-based",
      "second-line antiretroviral therapy for at least 6 months, with viral",
      "load >= 1000 copies/mL within the previous 2 months; participants of a",
      "randomised trial of a home-based modified directly administered",
      "adherence intervention. Backbones: tenofovir/lamivudine (82% of",
      "visits), abacavir/didanosine (7%), zidovudine/lamivudine (7%),",
      "abacavir/lamivudine (3%), tenofovir/emtricitabine (2%)."
    ),
    dose_range = "Atazanavir 300 mg once daily, boosted with ritonavir 100 mg",
    regions = "Zimbabwe (Harare)",
    n_observations = "82 quantifiable hair samples of 100 (18 below the 0.05 ng/mg limit), at enrolment and at 90-day follow-up",
    notes = paste(
      "Ngara 2020 Methods 'Source of data', Results 'Study participant",
      "details' and Table 1. Secondary analysis of data collected January",
      "2015 to May 2016. Hair cut closest to the scalp, median length 1 cm",
      "(0.5-1.5 cm); atazanavir in hair mean 2.5 ng/mg (SD 2.0, range",
      "0.07-8.65). No plasma concentrations were measured."
    )
  )

  ini({
    # Plasma PK -- fixed from Ngara 2020 refs 36/37 (Results 'Atazanavir model'; Table 3)
    lka <- fixed(log(0.44))
    label("Absorption rate constant (k12, 1/h)") # Table 3 'k12' 0.44 fixed; printed unit 'litres hour-1' is a typo for 1/h
    lcl <- fixed(log(10))
    label("Apparent clearance at the reference weight (CL/F, L/h)") # Table 3 'CL/F' 10 fixed
    lvc <- fixed(log(63.4))
    label("Apparent central volume at the reference weight (Vc, L)") # Table 3 'Vc' 63.4 fixed
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)") # Results 'Atazanavir model': exponent 0.75
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc (unitless)") # Results 'Atazanavir model': exponent 1

    # Hair layer (Equations 3-4; Table 3)
    lfhair <- log(0.16)
    label("Fraction of plasma elimination deposited in hair, reference covariates (Frac, unitless)") # Table 3 'Frac' 0.16 (SE 16%)
    lvhair <- fixed(log(1))
    label("Apparent hair volume (Vh, L)") # Table 3 'Vh' 1 fixed; Results: estimate was close to 1

    # Covariate effects on Frac, PsN linear form 1 + theta (Table 3)
    e_occ_enrol_fhair <- -0.30
    label("Fractional change in Frac at enrolment vs follow-up (unitless)") # Table 3 'Occasion (Enrolment)_Frac' -0.30
    e_adherence_vas_fhair <- 0.02
    label("Fractional change in Frac per unit of adherence score (1/%)") # Table 3 'Adherence score_Frac' 0.02
    e_bmi_for_age_thin_fhair <- 0.54
    label("Fractional change in Frac for thin vs normal BMI-for-age (unitless)") # Table 3 'Body Mass Index-for-age (Thin)_Frac' 0.54
    e_bmi_for_age_overweight_fhair <- -0.21
    label("Fractional change in Frac for overweight vs normal BMI-for-age (unitless)") # Table 3 'Body Mass Index-for-age (Over-weight)_Frac' -0.21
    e_caregiver_parent_fhair <- 0.53
    label("Fractional change in Frac for parent vs grandparent caregiver (unitless)") # Table 3 'Guardian (Parent)_Frac' 0.53
    e_caregiver_aunt_uncle_fhair <- 0.12
    label("Fractional change in Frac for aunt/uncle vs grandparent caregiver (unitless)") # Table 3 'Guardian (Uncle/Aunt)_Frac' 0.12
    e_caregiver_sibling_fhair <- -0.54
    label("Fractional change in Frac for sibling vs grandparent caregiver (unitless)") # Table 3 'Guardian (Sibling)_Frac' -0.54

    # IIV -- Table 3 'Variability' column, read as NONMEM omega^2
    etalka ~ fixed(0.45) # Table 3 'k12' Variability 0.45, held constant
    etalcl ~ 1.04 # Table 3 'CL/F' Variability 1.04 (SE 99%)
    etalvc ~ fixed(0.50) # Table 3 'Vc' Variability 0.50, held constant

    # Residual error -- Table 3 'eps ADD' / 'eps PROP' thetas with sigma fixed to 1
    addSd_Chair <- 0.30
    label("Additive residual SD of hair concentration (ng/mg)") # Table 3 'eps ADD' 0.30
    propSd_Chair <- 0.50
    label("Proportional residual SD of hair concentration (fraction)") # Table 3 'eps PROP' 0.50
  })

  model({
    # Occasion indicator: OCC 1 = enrolment, 2 = follow-up (reference)
    enrol <- 0
    if (OCC == 1) {
      enrol <- 1
    }

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    kel <- cl / vc

    # Frac with PsN-linear covariate factors; adherence centred on 84%
    fhair <- exp(lfhair) *
      (1 + e_occ_enrol_fhair * enrol) *
      (1 + e_adherence_vas_fhair * (ADHERENCE_VAS - 84)) *
      (1 +
        e_bmi_for_age_thin_fhair * BMI_FOR_AGE_THIN +
        e_bmi_for_age_overweight_fhair * BMI_FOR_AGE_OVERWEIGHT) *
      (1 +
        e_caregiver_parent_fhair * CAREGIVER_PARENT +
        e_caregiver_aunt_uncle_fhair * CAREGIVER_AUNT_UNCLE +
        e_caregiver_sibling_fhair * CAREGIVER_SIBLING)
    vhair <- exp(lvhair)

    # Equations 1-3: depot -> central -> hair (no loss from hair)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(hair) <- fhair * kel * central

    Cc <- central / vc
    # Equation 4
    Chair <- hair / vhair
    Chair ~ add(addSd_Chair) + prop(propSd_Chair)
  })
}
