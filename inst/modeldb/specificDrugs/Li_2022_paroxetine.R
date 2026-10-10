Li_2022_paroxetine <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral paroxetine immediate-release and sustained-release tablets in Chinese psychiatric inpatients receiving therapeutic drug monitoring, with daily dose on CL/F, formulation on V/F and sex on relative bioavailability (Li 2022)."
  reference <- "Li X-l, Huang S-q, Xiao T, Wang X-p, Kong W, Liu S-j, Zhang Z, Yang Y, Huang S-s, Ni X-j, Lu H-y, Zhang M, Wen Y-g, Shang D-w. Pharmacokinetics of immediate and sustained-release formulations of paroxetine: Population pharmacokinetic approach to guide paroxetine personalized therapy in Chinese psychotic patients. Front Pharmacol. 2022;13:966622. doi:10.3389/fphar.2022.966622"
  vignette <- "Li_2022_paroxetine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description = "Total daily paroxetine dose at the time of the record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (a) of the DOSE register entry, time-varying within a patient because the TDM cohort was titrated (20-75 mg/day). Power effect on CL/F centred at 40 mg/day: CL/F = 21.2 * (DOSE/40)^-1.03 (Li 2022 final-model equation, Section 3.2.1). Supply the total daily dose, not the per-administration amount; for once-daily regimens these are equal.",
      source_name = "dose"
    ),
    FORM_PAROXETINE_IR = list(
      description = "Paroxetine immediate-release tablet indicator (1 = immediate-release tablet, 0 = sustained-release tablet)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sustained-release tablet)",
      notes = "Li 2022 Methods 2.3.2: 'COV = 1 for immediate-release tablet; COV = 0 for sustained release tablet'. Enters V/F as 8850 * (1 - 0.666 * formulation), so the immediate-release tablet has the smaller V/F (2956 L) and the sustained-release tablet the reference V/F (8850 L). Per record (formulation was recorded per concentration in Table 1).",
      source_name = "formulation"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Li 2022 Methods 2.3.2: 'COV for the sex covariate is 1 (representing female) and 0 (representing male)'. Enters relative bioavailability as F1 = 1 + 0.475 * SEX, so females have 47.5% higher F than males (the reference).",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "paroxetine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "paroxetine", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 184L,
    n_studies = 1L,
    n_observations = 372L,
    age_range = "15-90 years",
    age_median = "37.5 years",
    weight_range = "58-96 kg (as printed in Table 1)",
    weight_median = "60 kg",
    height_median = "164.5 cm (range 142-180)",
    sex_female_pct = 56,
    race_ethnicity = c(Chinese = 100),
    disease_state = "Psychiatric inpatients (described as patients with psychosis) treated with paroxetine and monitored by therapeutic drug monitoring.",
    dose_range = "20-75 mg/day orally (immediate-release tablets and/or 25 mg sustained-release tablets); 85.5% of the 372 concentrations were on immediate-release and 14.5% on sustained-release tablets.",
    regions = "China (Affiliated Brain Hospital of Guangzhou Medical University, Guangzhou; retrospective TDM data 1 January 2019 to 31 May 2021).",
    co_medication = "Olanzapine 16.8%, tandospirone 14.1%, risperidone 13.0%, metoprolol 5.4% (none retained as a covariate).",
    notes = "Baseline demographics per Li 2022 Table 1. Serum samples were collected 10-22.5 h after the most recent dose (sparse trough-type TDM data, about 2 samples per subject), so Ka was fixed to a literature value and no absorption lag could be estimated. CYP2D6 genotype was examined but not retained."
  )

  ini({
    # Structural parameters - Li 2022 Table 2 and the final-model equations
    # (Section 3.2.1). Ka was fixed to the literature value used by
    # Venkatakrishnan & Obach 2005 and Nishimura 2016.
    lka <- fixed(log(0.908))
    label("Absorption rate constant Ka (1/h)") # Li 2022 Table 2 'Ka (1/h) 0.908 fixed'
    lcl <- log(21.2)
    label("Apparent clearance CL/F at DOSE = 40 mg/day (L/h)") # Li 2022 Table 2 'CL/F (L/h) 21.2', RSE 7.2%
    lvc <- log(8850)
    label("Apparent volume of distribution V/F for the sustained-release tablet (L)") # Li 2022 Table 2 'V/F (L) 8,850', RSE 17.2%

    # Covariate effects (Li 2022 Section 3.2.1 equations; Table 2 estimates)
    #   CL/F = 21.2 * (dose/40)^-1.03
    #   V/F  = 8850 * (1 - 0.666 * formulation)
    #   F1   = 1 + 0.475 * SEX
    e_dose_cl <- -1.03
    label("Power exponent of daily dose on CL/F (reference 40 mg/day; unitless)") # Li 2022 Table 2 'theta CL-dosage -1.03', RSE 5.7%
    e_form_ir_vc <- 0.666
    label("Fractional decrease in V/F for the immediate-release tablet (unitless)") # Li 2022 Table 2 'theta V-formulation 0.666', RSE 9.7%
    e_sexf_fdepot <- 0.475
    label("Fractional increase in relative bioavailability for females (unitless)") # Li 2022 Table 2 'theta F-sex 0.475', RSE 25.3%

    # IIV - exponential random effects (Li 2022 Eq. 1). Table 2 reports IIV
    # as CV%; converted with omega^2 = log(CV^2 + 1).
    etalcl ~ 0.2385136 # Li 2022 Table 2 'IIV (CV%)' CL/F 51.9 -> log(1 + 0.519^2)
    etalvc ~ 0.5289946 # Li 2022 Table 2 'IIV (CV%)' V/F 83.5 -> log(1 + 0.835^2)

    # Proportional residual error. Table 2 prints 'PRO (CV%) 0.0929' as a
    # fraction while the IIV rows of the same table are percentages; read as
    # the raw NONMEM $SIGMA variance, so the SD is sqrt(0.0929).
    propSd <- 0.3047950
    label("Proportional residual error (fraction)") # Li 2022 Table 2 'PRO (CV%) 0.0929' -> sqrt(0.0929)
  })

  model({
    # Individual parameters (Li 2022 Section 3.2.1 final-model equations)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (DOSE / 40)^e_dose_cl
    vc <- exp(lvc + etalvc) * (1 - e_form_ir_vc * FORM_PAROXETINE_IR)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Relative bioavailability: females 47.5% higher than males (F1 = 1 + 0.475 * SEX)
    f(depot) <- 1 + e_sexf_fdepot * SEXF

    # Dose in mg and vc in L give mg/L; x 1000 gives ng/mL (the paper's unit)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
