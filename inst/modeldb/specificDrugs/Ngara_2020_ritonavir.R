Ngara_2020_ritonavir <- function() {
  description <- paste(
    "One-compartment plasma PK model with first-order absorption plus a",
    "non-eliminating hair compartment for ritonavir given as the booster of",
    "ritonavir-boosted atazanavir (ATV/r 300/100 mg once daily) to",
    "HIV-infected Zimbabwean adolescents (Ngara 2020). Only hair",
    "concentrations were measured: the plasma parameters (ka 2.31 1/h, CL/F",
    "12.8 L/h, Vc 105 L) are fixed from published ritonavir models, and a",
    "fraction Frac (0.18) of the drug eliminated from plasma is deposited in",
    "hair (dA3/dt = Frac * k20 * A2), read out as Chair = A3 / Vh with Vh",
    "fixed to 1. Frac carries the covariate effects (PsN linear form): 42%",
    "lower at the enrolment visit than at the 90-day follow-up, and +2% per",
    "unit of self-reported adherence (visual analogue scale). Combined",
    "additive plus proportional residual error. The printed parameters",
    "predict hair levels about 34-fold above the paper's own observations",
    "(see the vignette). The companion atazanavir model from the same paper",
    "is Ngara_2020_atazanavir."
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
      analyte = "ritonavir",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "ritonavir",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    hair = list(
      analyte = "ritonavir",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    OCC = list(
      description = "Hair-sampling occasion (1 = enrolment, 2 = 90-day follow-up)",
      units = "(count)",
      type = "categorical",
      reference_category = "2 (follow-up visit)",
      notes = paste(
        "Two occasions per participant. Used as a fixed effect on Frac, not",
        "for inter-occasion variability: Frac is 42% lower at enrolment",
        "(OCC = 1) than at follow-up (Ngara 2020 Table 5,",
        "'Occasion (Enrolment)_Frac' = -0.42, follow-up the reference)."
      ),
      source_name = "Occasion"
    ),
    ADHERENCE_VAS = list(
      description = "Self-reported adherence on a 0-100% visual analogue scale",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear effect on Frac, +0.02 per unit (Ngara 2020 Table 5,",
        "'Adherence score_Frac'). The centring value is not printed; the",
        "effect is centred on 84, the rounded cohort mean (84.2, SD 18.1,",
        "range 30-100; Table 1). With this centring the factor is zero at",
        "a VAS of 34 and NEGATIVE below it, so the model as printed predicts a",
        "negative hair concentration for the least adherent participants; it",
        "is not clamped. See the vignette Errata."
      ),
      source_name = "Adherence VAS"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Ngara 2020 states allometric weight scaling of the fixed plasma",
        "parameters only for the atazanavir model; the ritonavir paragraph",
        "lists the fixed values without a weight term, so none is applied",
        "here. See the vignette Errata."
      )
    ),
    BMI_FOR_AGE_THIN = list(
      description = "WHO BMI-for-age category: thinness",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Frac by PsN SCM and not retained for ritonavir (Table 4); retained for atazanavir."
    ),
    BMI_FOR_AGE_OVERWEIGHT = list(
      description = "WHO BMI-for-age category: overweight",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Frac by PsN SCM and not retained for ritonavir (Table 4); retained for atazanavir."
    ),
    CAREGIVER_PARENT = list(
      description = "Primary caregiver is a parent",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Frac by PsN SCM and not retained for ritonavir (Table 4); retained for atazanavir."
    ),
    CAREGIVER_AUNT_UNCLE = list(
      description = "Primary caregiver is an aunt or uncle",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Frac by PsN SCM and not retained for ritonavir (Table 4); retained for atazanavir."
    ),
    CAREGIVER_SIBLING = list(
      description = "Primary caregiver is a sibling",
      units = "(binary)",
      type = "binary",
      notes = "Screened on Frac by PsN SCM and not retained for ritonavir (Table 4); retained for atazanavir."
    ),
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
      "adherence intervention."
    ),
    dose_range = "Ritonavir 100 mg once daily as the booster of atazanavir 300 mg",
    regions = "Zimbabwe (Harare)",
    n_observations = "88 quantifiable hair samples of 100 (12 below the 0.01 ng/mg limit), at enrolment and at 90-day follow-up",
    notes = paste(
      "Ngara 2020 Methods 'Source of data', Results 'Study participant",
      "details' and Table 1. Ritonavir in hair mean 0.5 ng/mg (SD 0.4, range",
      "0.01-1.39). No plasma concentrations were measured."
    )
  )

  ini({
    # Plasma PK -- fixed from Ngara 2020 refs 36/37 (Results 'Ritonavir model'; Table 5)
    lka <- fixed(log(2.31))
    label("Absorption rate constant (k12, 1/h)") # Table 5 'k12' 2.31 fixed; printed unit 'litres hour-1' is a typo for 1/h
    lcl <- fixed(log(12.8))
    label("Apparent clearance (CL/F, L/h)") # Table 5 'CL/F' 12.8 fixed
    lvc <- fixed(log(105))
    label("Apparent central volume (Vc, L)") # Table 5 'Vc' 105 fixed

    # Hair layer (Equations 3-4; Table 5)
    lfhair <- log(0.18)
    label("Fraction of plasma elimination deposited in hair, reference covariates (Frac, unitless)") # Table 5 'Frac' 0.18 (SE 16%)
    lvhair <- fixed(log(1))
    label("Apparent hair volume (Vh, L)") # Table 5 'Vh' 1 fixed; Results: estimate was close to 1

    # Covariate effects on Frac, PsN linear form 1 + theta (Table 5)
    e_occ_enrol_fhair <- -0.42
    label("Fractional change in Frac at enrolment vs follow-up (unitless)") # Table 5 'Occasion (Enrolment)_Frac' -0.42
    e_adherence_vas_fhair <- 0.02
    label("Fractional change in Frac per unit of adherence score (1/%)") # Table 5 'Adherence score_Frac' 0.02

    # IIV -- Table 5 'Variability' column, read as NONMEM omega^2
    etalka ~ fixed(0.45) # Table 5 'k12' Variability 0.45, held constant
    etalcl ~ 0.28 # Table 5 'CL/F' Variability 0.28 (SE 31%)
    etalvc ~ fixed(0.50) # Table 5 'Vc' Variability 0.50, held constant

    # Residual error -- Table 5 'eps ADD' / 'eps PROP' thetas (sigma fixed to 1, as Table 3)
    addSd_Chair <- 0.34
    label("Additive residual SD of hair concentration (ng/mg)") # Table 5 'eps ADD' 0.34
    propSd_Chair <- 0.26
    label("Proportional residual SD of hair concentration (fraction)") # Table 5 'eps PROP' 0.26
  })

  model({
    # Occasion indicator: OCC 1 = enrolment, 2 = follow-up (reference)
    enrol <- 0
    if (OCC == 1) {
      enrol <- 1
    }

    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # Frac with PsN-linear covariate factors; adherence centred on 84%
    fhair <- exp(lfhair) *
      (1 + e_occ_enrol_fhair * enrol) *
      (1 + e_adherence_vas_fhair * (ADHERENCE_VAS - 84))
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
