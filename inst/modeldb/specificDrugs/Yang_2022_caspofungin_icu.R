Yang_2022_caspofungin_icu <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous caspofungin in critically ill adults",
    "in the intensive care unit (the 'ICU-subgroup PPK model' of Yang 2022), pooling 51 ICU",
    "patients from Xi'an, China with 85 ICU patients digitised from seven published caspofungin",
    "studies. Clearance carries one dichotomised covariate, body weight >= 70 kg, entering as a",
    "(1 + theta x indicator) factor; volume has no covariate. Log-normal IIV on CL and V and a",
    "proportional residual error. The general-population model from the same paper is",
    "Yang_2022_caspofungin.",
    sep = " "
  )
  reference <- paste(
    "Yang Q, Zhang T, Zhang Y, Sun D, Zheng X, Du Q, Wang X, Cheng X, Xing J, Dong Y. The",
    "recommended dosage regimen for caspofungin in patients with higher body weight or",
    "hypoalbuminaemia will result in low exposure: five years of data based on a population",
    "pharmacokinetic model and Monte-Carlo simulations. Front Pharmacol. 2022;13:993330.",
    "doi:10.3389/fphar.2022.993330. PMCID: PMC9669616.",
    sep = " "
  )
  vignette <- "Yang_2022_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only DICHOTOMISED: model() forms wt_high <- (WT >= 70), the paper's WT* (Table 3",
        "footnote d uses WT* as defined in footnote c: WT* = 1 if weight >= 70 kg, otherwise 0).",
        "The reference (WT* = 0) is weight < 70 kg. Cohort 68.1 +/- 12.7 kg (30-139; Table 1).",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Significant on CL in the forward step (Table 2, dOFV -4.8) but not retained under the backward criterion (dOFV > 10.83); the Discussion attributes this to > 90% of ICU patients having albumin < 35 g/L. Cohort 28.0 +/- 6.9 g/L (Table 1)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy",
      units = "(binary)",
      type = "categorical",
      notes = "Collected for the ICU-subgroup model (Methods 2.2; 50 of 136 patients, Table 1) but not retained."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Collected (Methods 2.2) but not retained. Cohort 60.5 +/- 16.2 years (23-99; Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "categorical",
      notes = "Collected (Methods 2.2) but not retained. 53 of 136 patients female (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 136L,
    n_studies = 8L,
    age_range = "23-99 years",
    age_mean = "60.5 years (SD 16.2)",
    weight_range = "30-139 kg",
    weight_mean = "68.1 kg (SD 12.7)",
    sex_female_pct = 39.0,
    race_ethnicity = "Chinese (51 patients) and the populations of seven published ICU studies (85 patients, mostly European)",
    disease_state = paste(
      "Critically ill adults in the ICU receiving caspofungin; 50 of 136 on continuous renal",
      "replacement therapy. Mean serum albumin 28.0 g/L, > 90% below 35 g/L.",
      sep = " "
    ),
    dose_range = paste(
      "1 h intravenous infusion; 70 mg loading then 50 mg daily, with the same 70/35 and 70/70",
      "mg adjustments as the general cohort.",
      sep = " "
    ),
    regions = "China (Xi'an) plus the published cohorts of Spriet 2009, Spriet 2011, Weiler 2013, Muilwijk 2014, Sinnollareddy 2015, Roger 2017 and van der Elst 2017",
    notes = paste(
      "51 ICU patients from the First Affiliated Hospital of Xi'an Jiaotong University (sparse TDM",
      "samples) and 85 ICU patients whose concentrations were taken from seven published papers",
      "(dense and sparse sampling; Supplementary 3.1.2). Only sex, age, weight, dosing, albumin",
      "and CRRT status were available for the literature patients.",
      sep = " "
    )
  )

  ini({
    # Final model (Table 3, footnote d):
    #   CL = theta1 * (1 + theta3 * WT*) * exp(eta1)
    #   V  = theta2 * exp(eta2)
    lcl <- log(0.42); label("Clearance for weight < 70 kg (L/h)") # Table 3, theta1 = 0.42 L/h; the Abstract's 0.40 L/h is the bootstrap mean
    lvc <- log(10.20); label("Volume of distribution (L)") # Table 3, theta2 = 10.20 L (bootstrap 95% CI 7.70-19.00)
    e_wt_cl <- 0.34; label("Fractional change in CL for weight >= 70 kg (unitless)") # Table 3, theta3 = 0.34 (bootstrap 95% CI 0.11-0.70)

    # IIV: Table 3 'Inter-individual variability (%)' read as CV%, omega^2 = log(CV^2 + 1)
    etalcl ~ 0.018061 # Table 3, IIV CL 13.5 -> log(0.135^2 + 1)
    etalvc ~ 0.037697 # Table 3, IIV V 19.6 -> log(0.196^2 + 1)

    propSd <- 0.243; label("Proportional residual error (fraction)") # Table 3, residual variability 'Proportional (%)' = 24.3
  })
  model({
    # Dichotomised weight (Table 3 footnotes c and d)
    wt_high <- (WT >= 70)

    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * wt_high)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # Caspofungin is given as a 1 h intravenous infusion into central.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
