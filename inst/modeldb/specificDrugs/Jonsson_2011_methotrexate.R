Jonsson_2011_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (5-8 g/m^2) in Nordic children with acute lymphoblastic",
    "leukaemia (Jonsson 2011; n = 304 children, 1,284 high-dose courses).",
    "Linear first-order elimination from the central compartment. Body weight",
    "is the ONLY covariate and enters all four disposition parameters",
    "strictly linearly (per kg), which is the paper's central claim -- that",
    "weight-based rather than body-surface-area-based dose calculation is the",
    "better normalisation. Exponential between-subject variability on all",
    "four parameters, with a very large 109% CV on clearance. Transcribed",
    "from the Wang 2023 external evaluation, where it was by a wide margin",
    "the worst performer of the six models tested (population median PE",
    "442.04%, MPE 780.87%, RMSE 1182.24%).",
    sep = " "
  )
  reference <- paste(
    "Jonsson P, Skarby T, Heldrup J, Schroder H, Hoglund P (2011).",
    "High dose methotrexate treatment in children with acute lymphoblastic",
    "leukaemia may be optimised by a weight-based dose calculation.",
    "Pediatr Blood Cancer 57(1):41-46. doi:10.1002/pbc.22999.",
    "TRANSCRIPTION SOURCE: the full primary was not available; every value",
    "here is transcribed from Table 2 of the external evaluation",
    "Wang S, Yin Q, Yang M, Cheng Z, Xie F (2023). Pharmaceutics 15(2):569.",
    "doi:10.3390/pharmaceutics15020569. The residual error in particular is",
    "NOT from the primary -- see the propSd comment in ini().",
    "Re-extract from the primary when it is obtained; see",
    "vignette('Wang_2023_methotrexate') Errata.",
    sep = " "
  )
  vignette <- "Wang_2023_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482. Wang 2023 Table 2 records this model as '2 CMT' with states V1
  # (central) and V2 (peripheral). Amount units are umol because every
  # methotrexate concentration in Wang 2023 is reported in umol/L.
  # verified = FALSE for the specimen: the Jonsson 2011 primary, which would
  # state the assayed matrix, is not on disk, and neither Wang 2023 Table 2 nor
  # the primary's abstract names it.
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
        "Source column Weight. The ONLY covariate in the model, and it enters",
        "all four disposition parameters strictly linearly (per kg), per Wang",
        "2023 Table 2: 'CL (L/h) = Weight x 0.185 x e^eta_CL',",
        "'V1 (L) = Weight x 1.27 x e^eta_V1',",
        "'Q (L/h) = Weight x 0.017 x e^eta_Q',",
        "'V2 = Weight x 1.02 x e^eta_V2'. No exponent is estimated and no",
        "reference weight is used, so every coefficient is an absolute",
        "per-kilogram value. That linear-in-weight form is the paper's",
        "conclusion, per its abstract: 'Body weight improved the population",
        "pharmacokinetic model significantly more than any of the other",
        "patient characteristics, indicating that body weight may be the",
        "better way of dose normalisation.' Wang 2023 Table 1 records median",
        "19.0 kg (range 5.8-93.3) for this cohort, at which the typical",
        "clearance is 0.185 x 19 = 3.52 L/h and the typical central volume is",
        "1.27 x 19 = 24.1 L -- exactly the lowest CL and highest V1 that Wang",
        "2023 Results quotes across the six evaluated models."
      ),
      source_name        = "Weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened and not retained. The primary's abstract lists the candidate",
        "set: 'A population pharmacokinetic model was developed with data from",
        "1284 HDMTX courses in 304 children evaluating age, height, weight,",
        "body surface area (BSA), sex, serum creatinine and serum alanine",
        "aminotransferase as potential covariates.' Wang 2023 Table 1 records",
        "median 5.0 years (range 0.4-17.8)."
      )
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and not retained (same abstract sentence as AGE). Wang 2023 Table 1 median 110 cm (range 63-192)."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened and not retained -- deliberately, and this is the paper's",
        "headline finding: body weight 'improved the population",
        "pharmacokinetic model significantly more'. BSA nonetheless remains",
        "the DOSING descriptor for the protocol (5-8 g/m^2), so a simulation",
        "still needs it to compute the dose amount even though it enters no",
        "model equation. Wang 2023 Table 1 median 0.76 m^2 (range 0.31-2.19)."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and not retained (same abstract sentence as AGE). Wang 2023 Table 1 records M/F as 'NA' for this study."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened and not retained (same abstract sentence as AGE). Wang 2023 Table 1 records Scr as 'NA' for this study."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened and not retained (same abstract sentence as AGE). Wang 2023 Table 1 median 31.2 U/L (range 6.0-228.6)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 304L,
    n_studies      = 1L,
    age_range      = "0.4 to 17.8 years",
    age_median     = "5.0 years",
    weight_range   = "5.8 to 93.3 kg",
    weight_median  = "19.0 kg",
    height_range   = "63 to 192 cm (median 110)",
    bsa_range      = "0.31 to 2.19 m^2 (median 0.76)",
    disease_state  = "Childhood acute lymphoblastic leukaemia (ALL) receiving high-dose methotrexate.",
    hepatic_function = "ALT median 31.2 U/L (range 6.0-228.6).",
    dose_range     = "5 to 8 g/m^2 intravenous high-dose methotrexate -- the highest dose range of the six cohorts Wang 2023 evaluated.",
    regions        = "Sweden and Denmark (NOPHO ALL protocol centres).",
    notes          = paste(
      "Demographics from Wang 2023 Table 1; n = 304 children contributing",
      "1,284 high-dose methotrexate courses, confirmed from the primary's own",
      "abstract, which also notes the model was developed alongside an outcome",
      "analysis in 340 patients. Sex distribution and serum creatinine are not",
      "available (Wang 2023 Table 1 records both as 'NA').",
      "Wang 2023 attributes this model's poor external performance to its",
      "sparse covariate model and very large unexplained clearance",
      "variability: 'a large variability (109%) was found for MTX clearance in",
      "this model, and probably the most influential covariates were not",
      "identified to reduce the variability', and 'the Jonsson et al. model",
      "only considered the effect of body weight on a drug's PK, which seems",
      "to be significantly insufficient for our data'. Note too that this",
      "cohort's 5-8 g/m^2 dose range does not overlap the 1-5 g/m^2 of Wang's",
      "evaluation cohort, so the external evaluation extrapolated below the",
      "doses the model was built on."
    )
  )

  ini({
    # Structural parameters -- Wang 2023 Table 2, Jonsson row. Every one is a
    # per-kilogram coefficient; there is no reference weight and no exponent.
    # No uncertainty (RSE, CI or SE) is reproduced in Wang 2023, so none is
    # recorded here.
    lcl <- log(0.185); label("Clearance per kg body weight (L/h/kg)")                     # Wang 2023 Table 2: 'CL (L/h) = Weight x 0.185 x e^eta_CL'
    lvc <- log(1.27);  label("Central volume of distribution per kg body weight (L/kg)")  # Wang 2023 Table 2: 'V1 (L) = Weight x 1.27 x e^eta_V1'
    lq  <- log(0.017); label("Intercompartmental clearance per kg body weight (L/h/kg)")  # Wang 2023 Table 2: 'Q (L/h) = Weight x 0.017 x e^eta_Q'
    lvp <- log(1.02);  label("Peripheral volume of distribution per kg body weight (L/kg)")  # Wang 2023 Table 2: 'V2 = Weight x 1.02 x e^eta_V2'

    # Between-subject variability, exponential. Wang 2023 Table 2 reports these
    # as CV percentages in its 'IIV (%)' column; converted with the log-normal
    # identity omega^2 = log(1 + CV^2). The 109% CV on clearance is the largest
    # of the six models Wang evaluated and is what Wang blames for this model's
    # extreme population-prediction error.
    etalcl ~ 0.783034  # Wang 2023 Table 2 IIV CL = 109.0% -> log(1 + 1.090^2)
    etalvc ~ 0.065413  # Wang 2023 Table 2 IIV V1 = 26.0%  -> log(1 + 0.260^2)
    etalq  ~ 0.047265  # Wang 2023 Table 2 IIV Q  = 22.0%  -> log(1 + 0.220^2)
    etalvp ~ 0.176974  # Wang 2023 Table 2 IIV V2 = 44.0%  -> log(1 + 0.440^2)

    # Residual error. IMPORTANT PROVENANCE NOTE: this value is NOT a Jonsson
    # 2011 estimate. Wang 2023 Table 2 records this model's RUV as 'NR' (not
    # reported), and Wang supplied one for the evaluation per its Methods 2.3:
    # 'In addition, a proportional error of 30% was assumed for the residual
    # variability if the models did not report this information.' The 30%
    # assumption is therefore what Wang actually evaluated, and it is carried
    # here so this file reproduces the evaluated model, but it must not be read
    # as a published estimate. Re-extraction from the primary should replace
    # it. See the vignette Errata. Wrapped in fixed() because it was held at
    # 30% by assumption rather than estimated from data -- nothing in either
    # source estimated this quantity.
    propSd <- fixed(0.30); label("Proportional residual error (fraction) -- ASSUMED, not published")  # Wang 2023 Methods 2.3 assumption; Table 2 RUV = 'NR'
  })

  model({
    # Individual PK parameters. All four are strictly linear in body weight.
    cl <- exp(lcl + etalcl) * WT
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq  + etalq)  * WT
    vp <- exp(lvp + etalvp) * WT

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
    Cc ~ prop(propSd)
  })
}
