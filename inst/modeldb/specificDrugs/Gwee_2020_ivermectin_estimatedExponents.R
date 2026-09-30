Gwee_2020_ivermectin_estimatedExponents <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for a",
    "single oral dose of ivermectin in Indigenous Australian children aged 5",
    "to 15 years and weighing more than 15 kg who were treated for scabies",
    "(Gwee 2020, alternative model with ESTIMATED weight exponents). Body",
    "weight scales CL/F (exponent 0.944) and Vc/F (exponent 2.16), referenced",
    "to 37.55 kg; Q/F and Vp/F are not weight-scaled; ka is fixed to 0.5 /h.",
    "The authors preferred the fixed-exponent model (Gwee_2020_ivermectin) as",
    "final because both gave the same dose inference; this model fitted the",
    "data better and is provided for completeness."
  )
  reference <- paste(
    "Gwee A, Duffull S, Zhu X, Tong SYC, Cranswick N, McWhinney B, Ungerer J,",
    "Francis J, Steer AC (2020). Population pharmacokinetics of ivermectin for",
    "the treatment of scabies in Indigenous Australian children.",
    "PLoS Negl Trop Dis 14(12):e0008886. doi:10.1371/journal.pntd.0008886.",
    "Parameter values are from the supporting information S2 Text and S2 Table."
  )
  vignette <- "Gwee_2020_ivermectin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed (single-dose study). Power model referenced to 37.55 kg on",
        "CL/F and Vc/F only (S2 Text: 'CL/F(L/h) = 6.73 *",
        "(Weight/37.55)^0.944', 'Vc/F(L) = 160 * (Weight/37.55)^2.16',",
        "'Q/F(L/h) = 2.95', 'Vp/F(L) = 447'). The study enrolled children",
        "weighing 18.5 to 74.5 kg."
      ),
      source_name = "Weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested (Methods) and not retained (Results: 'No other covariates were found to improve the fit')."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Tested with females as the reference category (Methods); not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "ivermectin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ivermectin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ivermectin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 26,
    n_studies = 1,
    n_observations = 48,
    age_range = "5.5-14.9 years",
    age_median = "10.9 years",
    weight_range = "18.5-74.5 kg",
    weight_median = "37.6 kg",
    sex_female_pct = 57.7,
    race_ethnicity = c(`Indigenous Australian` = 100),
    disease_state = paste(
      "Scabies that failed topical therapy (24 children), suspected",
      "permethrin allergy (1) or severe scabies with secondary bacterial",
      "infection (1)."
    ),
    dose_range = paste(
      "Single oral dose of 200 ug/kg rounded to the nearest whole or half 3 mg",
      "tablet, given with food (median 190 ug/kg, range 120-230 ug/kg)."
    ),
    regions = "Remote community in Arnhem Land, Northern Territory, Australia",
    notes = paste(
      "Same cohort and data as Gwee_2020_ivermectin (ITCH study, April 2016 to",
      "March 2018; two samples per child spanning 2 hours to 14 days after",
      "dosing)."
    )
  )

  ini({
    lka <- fixed(log(0.5)); label("Absorption rate constant ka (1/h)") # S2 Table row 'Absorption rate constant (TVka, 1/h)' = '0.5 fixed'
    lcl <- log(6.73); label("Apparent clearance CL/F at 37.55 kg (L/h)") # S2 Table row 'TVCL, L/h' = 6.73 (RSE 23%); S2 Text CL/F equation
    lvc <- log(160); label("Apparent central volume Vc/F at 37.55 kg (L)") # S2 Table row 'TVVc, L' = 160 (RSE 28%); S2 Text Vc/F equation
    lq <- log(2.95); label("Apparent intercompartmental clearance Q/F (L/h; not weight-scaled)") # S2 Table row 'TVQ, L/h' = 2.95 (RSE 37%); S2 Text 'Q/F(L/h) = 2.95'
    lvp <- log(447); label("Apparent peripheral volume Vp/F (L; not weight-scaled)") # S2 Table row 'TVVp, L' = 447 (RSE 10%); S2 Text 'Vp/F(L) = 447'

    e_wt_cl <- 0.944; label("Allometric exponent of body weight on CL/F (unitless)") # S2 Table row 'Exponent for body weight on clearance' = 0.944 (RSE 3%)
    e_wt_vc <- 2.16; label("Allometric exponent of body weight on Vc/F (unitless)") # S2 Table row 'Exponent for body weight on central volume' = 2.16 (RSE 4%)

    # S2 Table 'Between subject variability (%CV)' block, printed one row
    # above its labels exactly as in S1 Table (see Gwee_2020_ivermectin).
    # Read in order: CV CL 68.7%, Vc 86.9%, Q 108.7%, Vp 29.8%; correlations
    # CL-Vc -0.042, CL-Q 0.399, CL-Vp -0.958, Vc-Q -0.934, Vc-Vp -0.239,
    # Q-Vp -0.125. Variances are (CV/100)^2 and covariances r * SD_i * SD_j.
    # The printed (rounded) correlation matrix is non-positive-definite
    # (smallest eigenvalue -0.001), so every off-diagonal is multiplied by
    # 0.99; the variances are unchanged.
    etalcl + etalvc + etalq + etalvp ~ c(
      0.471969,
      0.99 * -0.042 * 0.687 * 0.869, 0.755161,
      0.99 * 0.399 * 0.687 * 1.087, 0.99 * -0.934 * 0.869 * 1.087, 1.181569,
      0.99 * -0.958 * 0.687 * 0.298, 0.99 * -0.239 * 0.869 * 0.298, 0.99 * -0.125 * 1.087 * 0.298, 0.088804
    )

    expSd <- 0.062; label("Residual error, additive on the log scale (SD)") # S2 Table row 'Residual error, sigma additive in log domain' = 0.062
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 37.55)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 37.55)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; x 1000 gives the paper's ug/L.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
