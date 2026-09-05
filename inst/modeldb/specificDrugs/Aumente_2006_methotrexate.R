Aumente_2006_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (1.23-5.23 g/m^2) in Spanish children with acute",
    "lymphoblastic leukaemia (Aumente 2006; n = 37). Linear first-order",
    "elimination from the central compartment, with the peripheral",
    "distribution published as the micro-constants K12 and K21 rather than as",
    "Q and Vp. Clearance and central volume are STRATIFIED on a 10-year age",
    "cut-off: below 10 years clearance scales with body weight to the power",
    "0.876, at or above 10 years it is linear in weight, and the central",
    "volume is linear in weight in both strata with different coefficients.",
    "Exponential between-subject variability on CL, V1, K12 and K21, and a",
    "combined proportional plus additive residual error. Transcribed from the",
    "Wang 2023 external evaluation, which found this model unbiased for",
    "individual predictions (median PE 6.52%) but markedly under-predicting",
    "for population predictions (median PE -31.70%).",
    sep = " "
  )
  reference <- paste(
    "Aumente D, Buelga DS, Lukas JC, Gomez P, Torres A, Garcia MJ (2006).",
    "Population pharmacokinetics of high-dose methotrexate in children with",
    "acute lymphoblastic leukaemia. Clin Pharmacokinet 45(12):1227-1238.",
    "doi:10.2165/00003088-200645120-00007.",
    "TRANSCRIPTION SOURCE: the Aumente 2006 primary was not available; every",
    "value here is transcribed from Table 2 of the external evaluation",
    "Wang S, Yin Q, Yang M, Cheng Z, Xie F (2023). External Evaluation of",
    "Population Pharmacokinetic Models of Methotrexate for Model-Informed",
    "Precision Dosing in Pediatric Patients with Acute Lymphoid Leukemia.",
    "Pharmaceutics 15(2):569. doi:10.3390/pharmaceutics15020569.",
    "Re-extract from the primary when it is obtained; see",
    "vignette('Wang_2023_methotrexate') Errata.",
    sep = " "
  )
  vignette <- "Wang_2023_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482. Wang 2023 Table 2 records this model as '2 CMT' and names the
  # states V1 (central) and, implicitly through K12/K21, one peripheral
  # compartment. Amount units are umol because Wang reports the additive
  # residual error for this model as 0.0035 umol/L and all methotrexate
  # concentrations throughout the paper in umol/L. specimen is 'plasma':
  # Wang 2023 Methods 2.2 describes plasma sampling for the evaluation cohort
  # and Table 2 lists 'V1: central volume of distribution'. verified = FALSE
  # for the specimen because the Aumente 2006 primary, which would state the
  # assayed matrix for THIS cohort, is not on disk.
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column Weight. Enters BOTH clearance and central volume, with",
        "a different form in each age stratum, per Wang 2023 Table 2:",
        "'CL (L/h, age > 10 years) = 0.149 x Weight',",
        "'CL (L/h, age < 10 years) = 0.287 x Weight^0.876',",
        "'V1 (L, age > 10 years) = 0.437 x Weight',",
        "'V1 (L, age < 10 years) = 0.465 x Weight'. Note the weight terms are",
        "NOT normalised to a reference weight -- the coefficients are absolute",
        "(L/h per kg, or per kg^0.876; and L per kg). Wang 2023 Table 1",
        "records median 24.2 kg (range 7.5-80.0) for this cohort."
      ),
      source_name        = "Weight"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column age. Used ONLY as a binary stratifier at 10 years, not",
        "as a continuous effect: Wang 2023 Table 2 gives one pair of CL and V1",
        "equations for 'age > 10 years' and another for 'age < 10 years'.",
        "The table leaves age EXACTLY 10 undefined (it prints strict '>' and",
        "strict '<'). This implementation assigns age exactly 10 to the older",
        "stratum, i.e. the split is AGE < 10 vs AGE >= 10; the alternative",
        "assignment changes the prediction only for a subject whose recorded",
        "age is exactly 10.0 years. Wang 2023 Table 1 records median 5.0 years",
        "(range 0.5-17.0) for this cohort. See the vignette Errata."
      ),
      source_name        = "age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 37L,
    n_studies      = 1L,
    age_range      = "0.5 to 17.0 years",
    age_median     = "5.0 years",
    weight_range   = "7.5 to 80.0 kg",
    weight_median  = "24.2 kg",
    height_range   = "69 to 174 cm (median 115)",
    bsa_range      = "0.30 to 1.90 m^2 (median 1.07)",
    sex_female_pct = 100 * 20 / 37,
    disease_state  = "Childhood acute lymphoblastic leukaemia (ALL) receiving high-dose methotrexate.",
    renal_function = "Serum creatinine median 0.5 mg/dL (range 0.3-0.8).",
    dose_range     = "1.23 to 5.23 g/m^2 intravenous high-dose methotrexate.",
    regions        = "Spain (Reina Sofia University Hospital, Cordoba).",
    notes          = paste(
      "All demographics are from Wang 2023 Table 1, which tabulates the six",
      "evaluated cohorts side by side; the Aumente 2006 primary is not on",
      "disk. Wang 2023 Results reports that this model, together with the",
      "Medellin-Garibay and Zhang models, 'showed the best predictive",
      "performance across all the different tests', while Figure 2 panel A",
      "shows it is the one model whose individual predictions retain a",
      "systematic bias."
    )
  )

  ini({
    # Structural parameters -- Wang 2023 Table 2, Aumente row. The paper
    # publishes a stratified pair for CL and for V1, so each stratum takes an
    # explicit suffix per parameter-names.md 'Stratum-suffixed parameters'.
    # Neither stratum keeps the bare canonical name.
    #
    # No uncertainty (RSE, CI or SE) is reproduced in Wang 2023 for any of
    # these values, so none is recorded here.
    lcl_agelt10 <- log(0.287); label("Clearance coefficient on weight^0.876, age < 10 years (L/h/kg^0.876)")  # Wang 2023 Table 2: 'CL (L/h, age < 10 years) = 0.287 x Weight^0.876'
    lcl_agege10 <- log(0.149); label("Clearance coefficient on weight, age >= 10 years (L/h/kg)")             # Wang 2023 Table 2: 'CL (L/h, age > 10 years) = 0.149 x Weight'
    lvc_agelt10 <- log(0.465); label("Central volume coefficient on weight, age < 10 years (L/kg)")           # Wang 2023 Table 2: 'V1 (L, age < 10 years) = 0.465 x Weight'
    lvc_agege10 <- log(0.437); label("Central volume coefficient on weight, age >= 10 years (L/kg)")          # Wang 2023 Table 2: 'V1 (L, age > 10 years) = 0.437 x Weight'

    # Only the younger stratum's clearance carries a non-unit weight exponent;
    # the older stratum's CL and both strata's V1 are linear in weight by the
    # published model form, so there is nothing to estimate for those and only
    # this one exponent takes a stratum suffix.
    e_wt_cl_agelt10 <- 0.876; label("Power exponent of body weight on CL, age < 10 years (unitless)")  # Wang 2023 Table 2: exponent in '0.287 x Weight^0.876'

    # Peripheral distribution. Aumente publishes MICRO-CONSTANTS rather than a
    # Q / Vp pair, and puts an eta on each of them, so lk12 / lk21 are kept as
    # the stored fixed effects (preserving the published numbers and giving
    # etalk12 / etalk21 a matching fixed effect). model() converts them to
    # q / vp before writing the ODEs -- see the comment there.
    lk12 <- log(0.0155); label("Central to peripheral transfer rate constant (1/h)")  # Wang 2023 Table 2: 'K12 (1/h) = 0.0155'
    lk21 <- log(0.0724); label("Peripheral to central transfer rate constant (1/h)")  # Wang 2023 Table 2: 'K21 (1/h) = 0.0724'

    # Between-subject variability, exponential. Wang 2023 Table 2 reports these
    # as CV percentages in its 'IIV (%)' column; converted with the log-normal
    # identity omega^2 = log(1 + CV^2). A single eta per parameter is shared
    # across both age strata -- Wang reports one IIV for CL and one for V1, not
    # one per stratum.
    etalcl  ~ 0.160322  # Wang 2023 Table 2 IIV CL  = 41.7% -> log(1 + 0.417^2)
    etalvc  ~ 0.159612  # Wang 2023 Table 2 IIV V1  = 41.6% -> log(1 + 0.416^2)
    etalk12 ~ 0.042354  # Wang 2023 Table 2 IIV K12 = 20.8% -> log(1 + 0.208^2)
    etalk21 ~ 0.116808  # Wang 2023 Table 2 IIV K21 = 35.2% -> log(1 + 0.352^2)

    # Residual error -- combined proportional and additive.
    propSd <- 0.162;  label("Proportional residual error (fraction)")  # Wang 2023 Table 2 RUV: 'Prop = 16.2%'
    addSd  <- 0.0035; label("Additive residual error (umol/L)")        # Wang 2023 Table 2 RUV: 'Add = 0.0035 umol/L'
  })

  model({
    # Age stratifier. Wang 2023 Table 2 splits at 10 years; age exactly 10 is
    # assigned to the older stratum (see covariateData$AGE notes).
    ageLt10 <- (AGE < 10)

    # Individual PK parameters. The stratum switch is applied on the LOG scale
    # so that the etas stay mu-referenced: within each branch the linear
    # predictor is log(coefficient) + exponent * log(WT).
    lclInd <- ageLt10 * (lcl_agelt10 + e_wt_cl_agelt10 * log(WT)) +
      (1 - ageLt10) * (lcl_agege10 + log(WT))
    lvcInd <- ageLt10 * lvc_agelt10 + (1 - ageLt10) * lvc_agege10 + log(WT)

    cl <- exp(lclInd + etalcl)
    vc <- exp(lvcInd + etalvc)

    # The published micro-constants are converted to flow / volume form before
    # the ODEs are written. rxSolve() defaults to useLinCmt = TRUE, and a
    # constant-coefficient two-compartment system written directly in k12/k21
    # is silently rewritten to a ONE-compartment closed form, dropping
    # peripheral1 entirely. Going through q and vp is exact
    # (Q = K12 x V1, Vp = Q / K21) and converts correctly with and without the
    # flag. The intermediates are deliberately named kcp / kpc rather than
    # k12 / k21 so the canonical names below are derived last.
    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)
    q   <- kcp * vc
    vp  <- q / kpc

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition; methotrexate is infused intravenously
    # straight into the central compartment, so there is no absorption phase
    # and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma methotrexate concentration in umol/L (amounts in umol, volumes L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
