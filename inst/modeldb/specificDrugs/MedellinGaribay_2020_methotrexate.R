MedellinGaribay_2020_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (2.89 +/- 0.9 g/m^2 as a 24 h continuous infusion) in",
    "Mexican children with acute lymphoblastic leukaemia",
    "(Medellin-Garibay 2020; n = 41 in the model-building cohort). Linear",
    "first-order elimination from the central compartment. Body surface area",
    "raises clearance through an un-normalised power term (exponent 0.62) and",
    "body weight scales the central volume linearly; intercompartmental",
    "clearance and peripheral volume carry no covariate. Exponential",
    "between-subject variability on CL, Vc and Vp, and a proportional",
    "residual error. Transcribed from the Wang 2023 external evaluation,",
    "where it was the only model whose population predictions met the",
    "pre-specified +/-20% median-bias criterion (median PE -7.56%).",
    sep = " "
  )
  reference <- paste(
    "Medellin-Garibay SE, Hernandez-Villa N, Correa-Gonzalez LC,",
    "Morales-Barragan MN, Valero-Rivera KP, Resendiz-Galvan JE,",
    "Ortiz-Zamudio JJ, Milan-Segovia RD, Romano-Moreno S (2020).",
    "Population pharmacokinetics of methotrexate in Mexican pediatric",
    "patients with acute lymphoblastic leukemia.",
    "Cancer Chemother Pharmacol 85(1):21-31. doi:10.1007/s00280-019-03977-1.",
    "TRANSCRIPTION SOURCE: the full primary was not available; the values here",
    "are transcribed from Table 2 of the external evaluation",
    "Wang S, Yin Q, Yang M, Cheng Z, Xie F (2023). Pharmaceutics 15(2):569.",
    "doi:10.3390/pharmaceutics15020569, and the four structural values were",
    "independently confirmed against the primary's own published abstract.",
    "Re-extract from the primary when it is obtained; see",
    "vignette('Wang_2023_methotrexate') Errata.",
    sep = " "
  )
  vignette <- "Wang_2023_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482. Wang 2023 Table 2 records this model as '2 CMT' with states V1
  # (central) and V2 (peripheral). Amount units are umol because every
  # methotrexate concentration in Wang 2023 is reported in umol/L. The
  # primary's abstract states 'Plasma concentrations of MTX were determined in
  # blood samples collected at 24, 36, 42 or 48 h post-infusion, by means of
  # the CMIA immunoassay', which fixes the specimen as plasma.
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BSA. Power-form effect on clearance, per Wang 2023",
        "Table 2: 'CL (L/h) = 6.5 x BSA^0.62 x e^eta_CL'. Confirmed verbatim",
        "in the primary's own abstract: 'The population pharmacokinetic model",
        "obtained was: CL (L/h) = 6.5 x BSA^0.62, Vc (L) = 0.36 x Weight,",
        "Q (L/h) = 0.41 and Vp (L) = 3.2', which also states that 'body",
        "surface area (BSA) was the covariate that influences on MTX total",
        "CL'. Note the term is NOT normalised to a reference BSA -- the",
        "coefficient 6.5 is clearance at BSA = 1 m^2, and at the cohort median",
        "BSA of 0.79 m^2 the typical clearance is 6.5 x 0.79^0.62 = 5.6 L/h.",
        "Wang 2023 Table 1 records BSA median 0.79 m^2 (range 0.41-1.6)."
      ),
      source_name        = "BSA"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column Weight. Strictly linear (per-kg) effect on the central",
        "volume, per Wang 2023 Table 2: 'V1 (L) = 0.36 x Weight x e^eta_V1',",
        "confirmed verbatim in the primary's abstract ('Vc (L) = 0.36 x",
        "Weight'). Not normalised to a reference weight; 0.36 is L per kg.",
        "Wang 2023 Table 1 records median 21.2 kg (range 8.0-57.3), giving a",
        "typical central volume of 7.6 L -- the smallest of the six models",
        "Wang evaluated, which Wang's Results quotes as '7.5 [L]'."
      ),
      source_name        = "Weight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "1.0 to 15.0 years",
    age_median     = "5.0 years",
    weight_range   = "8.0 to 57.3 kg",
    weight_median  = "21.2 kg",
    height_range   = "115 +/- 24 cm (mean +/- SD)",
    bsa_range      = "0.41 to 1.6 m^2 (median 0.79)",
    disease_state  = "Childhood acute lymphoblastic leukaemia (ALL) receiving high-dose methotrexate as a 24 h continuous intravenous infusion.",
    renal_function = "Serum creatinine 0.37 +/- 0.11 mg/dL (mean +/- SD).",
    hepatic_function = "ALT median 19.1 U/L (range 7.2-98.2); AST median 25.8 U/L (range 13.9-112.7).",
    dose_range     = "2.89 +/- 0.9 g/m^2 (mean +/- SD) intravenous high-dose methotrexate over 24 h.",
    regions        = "Mexico (Hospital Central 'Dr. Ignacio Morones Prieto', San Luis Potosi).",
    notes          = paste(
      "Demographics from Wang 2023 Table 1, which reports n = 41 for the",
      "model-building cohort; the primary's abstract describes a prospective",
      "study in 50 children aged 1-15 years, with predictive performance",
      "'evaluated by external validation in a different group of patients', so",
      "the 41 counts the subset used to build the model reported here.",
      "Sex distribution is not available (Wang 2023 Table 1 records M/F as",
      "'NA' for this study). Wang 2023 Results groups this model with the",
      "Aumente and Zhang models as those that 'showed the best predictive",
      "performance across all the different tests', and it is the only one of",
      "the six whose POPULATION median prediction error fell inside the",
      "pre-specified +/-20% band."
    )
  )

  ini({
    # Structural parameters -- Wang 2023 Table 2, Medellin-Garibay row. All
    # four reproduce verbatim in the primary's own abstract. No uncertainty
    # (RSE, CI or SE) is reproduced in either source, so none is recorded here.
    lcl <- log(6.5);  label("Clearance at 1 m^2 body surface area (L/h)")      # Wang 2023 Table 2: 'CL (L/h) = 6.5 x BSA^0.62 x e^eta_CL'; primary abstract 'CL (L/h) = 6.5 x BSA^0.62'
    lvc <- log(0.36); label("Central volume of distribution per kg (L/kg)")    # Wang 2023 Table 2: 'V1 (L) = 0.36 x Weight x e^eta_V1'; primary abstract 'Vc (L) = 0.36 x Weight'
    lq  <- log(0.41); label("Intercompartmental clearance (L/h)")              # Wang 2023 Table 2: 'Q (L/h) = 0.41'; primary abstract 'Q (L/h) = 0.41'
    lvp <- log(3.2);  label("Peripheral volume of distribution (L)")           # Wang 2023 Table 2: 'V2 (L) = 3.2 x e^eta_V2'; primary abstract 'Vp (L) = 3.2'

    # Covariate effect. Only clearance carries a power exponent; the central
    # volume is linear in weight by the published model form and so needs no
    # exponent parameter.
    e_bsa_cl <- 0.62; label("Power exponent of body surface area on CL (unitless)")  # Wang 2023 Table 2 and primary abstract: exponent in '6.5 x BSA^0.62'

    # Between-subject variability, exponential. Wang 2023 Table 2 reports these
    # as CV percentages in its 'IIV (%)' column; converted with the log-normal
    # identity omega^2 = log(1 + CV^2). Q carries no eta in the published
    # equation and none is added here.
    etalcl ~ 0.006701  # Wang 2023 Table 2 IIV CL = 8.2%  -> log(1 + 0.082^2)
    etalvc ~ 0.064927  # Wang 2023 Table 2 IIV V1 = 25.9% -> log(1 + 0.259^2)
    etalvp ~ 0.068863  # Wang 2023 Table 2 IIV V2 = 26.7% -> log(1 + 0.267^2)

    # Residual error.
    propSd <- 0.201; label("Proportional residual error (fraction)")  # Wang 2023 Table 2 RUV: 'Prop = 20.1%'
  })

  model({
    # Individual PK parameters. Both covariate terms are un-normalised: cl is
    # clearance per (m^2)^0.62 and vc is volume per kg.
    cl <- exp(lcl + etalcl) * BSA^e_bsa_cl
    vc <- exp(lvc + etalvc) * WT
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

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
