Gwee_2020_ivermectin <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for a",
    "single oral dose of ivermectin in Indigenous Australian children aged 5",
    "to 15 years and weighing more than 15 kg who were treated for scabies",
    "(Gwee 2020, final 'biological prior' model). Body weight scales CL/F and",
    "Q/F with an exponent fixed to 0.75 and Vc/F and Vp/F with an exponent",
    "fixed to 1, all referenced to 37.55 kg; ka is fixed to 0.5 /h.",
    "Between-subject variability is a full 4 x 4 block on CL/F, Vc/F, Q/F and",
    "Vp/F; the residual error is additive on the log scale."
  )
  reference <- paste(
    "Gwee A, Duffull S, Zhu X, Tong SYC, Cranswick N, McWhinney B, Ungerer J,",
    "Francis J, Steer AC (2020). Population pharmacokinetics of ivermectin for",
    "the treatment of scabies in Indigenous Australian children.",
    "PLoS Negl Trop Dis 14(12):e0008886. doi:10.1371/journal.pntd.0008886.",
    "Parameter values are from the supporting information S1 Text and S1 Table."
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
        "Time-fixed (single-dose study). Power model referenced to 37.55 kg,",
        "the cohort median (S1 Text: 'CL/F(L/h) = 6.94 * (Weight/37.55)^0.75',",
        "'Vc/F(L) = 88 * (Weight/37.55)', 'Q/F(L/h) = 13.7 *",
        "(Weight/37.55)^0.75', 'Vp/F(L) = 344 * (Weight/37.55)'). The study",
        "enrolled children weighing 18.5 to 74.5 kg; the authors used the",
        "model to extrapolate to children weighing 10 to 15 kg."
      ),
      source_name = "Weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested as a power-model covariate (Methods: 'Covariates were tested",
        "that were deemed biologically relevant (weight, age, sex)') and not",
        "retained (Results: 'No other covariates were found to improve the fit')."
      )
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested with females as the reference category and the effect of male",
        "sex estimated (Methods); not retained in the final model."
      )
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
      "Ivermectin Therapy in Children (ITCH) study, April 2016 to March 2018.",
      "Two plasma samples per child (median 2, range 1-2) drawn in one of five",
      "sampling-design groups spanning 2 hours to 14 days after dosing (Table 1).",
      "Demographics from the first paragraph of the Results; 11 of 26 (42%)",
      "children were male."
    )
  )

  ini({
    lka <- fixed(log(0.5)); label("Absorption rate constant ka (1/h)") # S1 Table row 'Absorption rate constant (TVka, 1/h)' = '0.5 fixed'
    lcl <- log(6.94); label("Apparent clearance CL/F at 37.55 kg (L/h)") # S1 Table row 'TVCL, L/h' = 6.94 (RSE 18%); S1 Text CL/F equation
    lvc <- log(88); label("Apparent central volume Vc/F at 37.55 kg (L)") # S1 Table row 'TVVc, L' = 88 (RSE 53%); S1 Text Vc/F equation
    lq <- log(13.7); label("Apparent intercompartmental clearance Q/F at 37.55 kg (L/h)") # S1 Table row 'TVQ, L/h' = 13.7 (RSE 14%); S1 Text Q/F equation
    lvp <- log(344); label("Apparent peripheral volume Vp/F at 37.55 kg (L)") # S1 Table row 'TVVp, L' = 344 (RSE 14%); S1 Text Vp/F equation

    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # S1 Table row 'Exponent for body weight on clearance' = '0.75 fixed'
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F (unitless)") # S1 Table row 'Exponent for body weight on central volume (theta WT,Vc)' = '1 fixed'
    e_wt_q <- fixed(0.75); label("Allometric exponent of body weight on Q/F (unitless)") # S1 Table row labelled 'clearance' but with symbol 'theta WT,Q' = '0.75 fixed'
    e_wt_vp <- fixed(1); label("Allometric exponent of body weight on Vp/F (unitless)") # S1 Table row labelled 'central volume' but with symbol 'theta WT,Vp' = '1 fixed'

    # S1 Table 'Between subject variability (%CV)' block. The value column is
    # printed one row above its labels: 10 values (4 CVs then 6 correlations
    # in the order CL-Vc, CL-Q, CL-Vp, Vc-Q, Vc-Vp, Q-Vp) sit against 11 label
    # rows, the first value on the block-heading row and the 'Correlation
    # (Q, Vp)' row empty. Read in order: CV CL 54%, Vc 237.1%, Q 40.8%,
    # Vp 38.9%; correlations -0.30, 0.70, 0.863, 0.467, 0.215, 0.964.
    # The %CV is taken as 100 * omega (SD scale), so each variance is
    # (CV/100)^2 and each covariance is r * SD_i * SD_j; the vignette gives
    # the evidence for this reading. The printed (rounded) correlation matrix
    # is slightly non-positive-definite (smallest eigenvalue -2e-5), so every
    # off-diagonal is multiplied by 0.99; the variances are unchanged.
    etalcl + etalvc + etalq + etalvp ~ c(
      0.2916,
      0.99 * -0.300 * 0.540 * 2.371, 5.621641,
      0.99 * 0.700 * 0.540 * 0.408, 0.99 * 0.467 * 2.371 * 0.408, 0.166464,
      0.99 * 0.863 * 0.540 * 0.389, 0.99 * 0.215 * 2.371 * 0.389, 0.99 * 0.964 * 0.408 * 0.389, 0.151321
    )

    expSd <- 0.198; label("Residual error, additive on the log scale (SD)") # S1 Table row 'Residual error, sigma additive in log domain' = 0.198
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 37.55)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 37.55)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 37.55)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 37.55)^e_wt_vp

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
