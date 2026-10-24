Gu_2023_naloxone <- function() {
  description <- paste(
    "One-compartment population PK model for naloxone after the sublingual",
    "buprenorphine/naloxone (4:1) tablet in Chinese male healthy volunteers",
    "and patients with opioid use disorder (Gu 2023). First-order absorption",
    "(no IIV on ka), bioavailability fixed to 0.01 at the 8 mg dose level",
    "with a power effect of dose (exponent -0.36), and combined additive and",
    "proportional residual error. The dose is entered as the labelled",
    "buprenorphine strength of the combination dose (8 mg for four",
    "2 mg/0.5 mg tablets carrying 2 mg naloxone), as in the published",
    "model. The",
    "buprenorphine and norbuprenorphine component is a separately fitted",
    "model (Gu_2023_buprenorphine)."
  )
  reference <- paste(
    "Gu M, Li A, Mak W, Dong F, Xu N, Zhang J, Shi Y, Zheng N, Tang Z, He Q,",
    "Ruan C, Guo W, Xiang X, Wang C, Han B, Zhu X (2023). Population",
    "pharmacokinetics of buprenorphine and naloxone sublingual combination",
    "in Chinese healthy volunteers and patients with opioid use disorder:",
    "Model-based dose optimization. Front. Pharmacol. 14:1089862.",
    "doi:10.3389/fphar.2023.1089862."
  )
  vignette <- "Gu_2023_buprenorphine_naloxone"
  units <- list(
    time = "h",
    dosing = "mg (labelled buprenorphine strength of the buprenorphine/naloxone dose; the naloxone content is one quarter of it)",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(
      analyte = "naloxone",
      units = "mg (labelled buprenorphine strength)",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "naloxone",
      units = "mg (labelled buprenorphine strength)",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    DOSE_BPN_SL_MG = list(
      description = "Sublingual buprenorphine dose of the buprenorphine/naloxone tablets administered (mg buprenorphine)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on naloxone bioavailability",
        "F = 0.01 * (DOSE_BPN_SL_MG / 8)^-0.36 (supplementary control stream",
        "'TVF1 = THETA(4)*(DOSE/8)**THETA(7)'; Table 3 'Exponent for dose on",
        "F' = -0.36). The same column carries the dose-group level (4, 8 or",
        "16 mg buprenorphine) as in Gu_2023_buprenorphine; it is also the",
        "value to supply as amt (see the model description)."
      ),
      source_name = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the stepwise covariate search (Section 4.4) and not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in the stepwise covariate search (Section 4.4) and not retained."
    ),
    NSWALLOWS = list(
      description = "Number of swallows performed after sublingual administration of the dose",
      units = "(count)",
      type = "continuous",
      notes = "Tested in the stepwise covariate search (Section 4.4) and not retained; dose rather than swallowing explained the fall in naloxone bioavailability."
    ),
    DIS_OUD = list(
      description = "Opioid use disorder patient (1) vs healthy volunteer (0)",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the stepwise covariate search (Section 4.4) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 46L,
    n_studies = 1L,
    age_range = "20-44 years",
    age_median = "27.5 years",
    weight_range = "53-88 kg",
    weight_median = "65.0 kg",
    sex_female_pct = 0,
    race_ethnicity = "Han Chinese (100%)",
    disease_state = "34 healthy volunteers (single dose) and 12 patients with opioid use disorder (multiple dose)",
    dose_range = paste(
      "Sublingual buprenorphine/naloxone tablets (2 mg buprenorphine +",
      "0.5 mg naloxone per tablet): single 4, 8 or 16 mg buprenorphine doses",
      "(1, 2 or 4 mg naloxone) in healthy volunteers and 8 mg buprenorphine",
      "(2 mg naloxone) once daily for 8 days in OUD subjects"
    ),
    regions = "China (Beijing Anding Hospital)",
    notes = paste(
      "Table 1. 484 naloxone plasma concentrations; 20.3% of single-dose and",
      "56.2% of multiple-dose samples were below the LLOQ of 20 ng/L and were",
      "handled with the M3 likelihood method."
    )
  )

  ini({
    # Table 3 (final naloxone model).
    lcl <- log(205)
    label("Naloxone clearance CL (L/h)") # Table 3 'CL NLX , L/h' = 205
    lvc <- log(104)
    label("Naloxone central volume Vc (L)") # Table 3 'V C NLX , L' = 104
    lka <- log(0.585)
    label("First-order absorption rate constant (1/h)") # Table 3 'ka,/h' = 0.585
    lfdepot <- fixed(log(0.01))
    label("Sublingual naloxone bioavailability at the 8 mg dose level (unitless)") # Table 3 'F' = 0.01 FIX (McDonald et al., 2018)
    e_dose_fdepot <- -0.36
    label("Power exponent of DOSE/8 mg on bioavailability (unitless)") # Table 3 'Exponent for dose on F' = -0.36

    # IIV printed as sqrt(omega^2) * 100, the convention verified for the
    # buprenorphine model of the same paper, so omega^2 = (CV/100)^2. No IIV
    # on ka (control stream $OMEGA 0 FIX).
    etalcl ~ 0.097344 # Table 3 IIV 'CL NLX' = 31.2% -> 0.312^2
    etalvc ~ 1.050625 # Table 3 IIV 'V C NLX' = 102.5% -> 1.025^2

    # Combined error W = sqrt(THETA(5)^2 * IPRED^2 + THETA(6)^2) with the
    # concentration in ng/L in the control stream; the additive SD is
    # converted to ng/mL here (9.24 ng/L = 0.00924 ng/mL).
    propSd <- 0.208
    label("Proportional residual error (fraction)") # Table 3 'sigma prop (%CV)' = 20.8%
    addSd <- 0.00924
    label("Additive residual error (ng/mL)") # Table 3 'sigma additive (sd)' = 9.24 ng/L
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)
    kel <- cl / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot) * (DOSE_BPN_SL_MG / 8)^e_dose_fdepot

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
