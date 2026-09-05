Hui_2019_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (2-5 g/m^2) in Chinese children with acute lymphoblastic",
    "leukaemia (Hui 2019; n = 36 ALL patients, 354 observations). Linear",
    "first-order elimination from the central compartment. Body surface area",
    "raises clearance and the central volume through power terms normalised",
    "to a 0.735 m^2 typical child, BSA-normalised eGFR raises clearance",
    "further, and age raises intercompartmental clearance. Exponential",
    "between-subject variability on CL and the peripheral volume, plus",
    "between-OCCASION variability on clearance, and a proportional residual",
    "error. Transcribed from the Wang 2023 external evaluation, where its",
    "numerical bias was acceptable (MPE -8.24% population / 8.78%",
    "individual) but the Bland-Altman plots showed an obvious trend in",
    "prediction error across the concentration range.",
    sep = " "
  )
  reference <- paste(
    "Hui KH, Chu HM, Fong PS, Cheng WTF, Lam TN (2019).",
    "Population Pharmacokinetic Study and Individual Dose Adjustments of",
    "High-Dose Methotrexate in Chinese Pediatric Patients with Acute",
    "Lymphoblastic Leukemia or Osteosarcoma.",
    "J Clin Pharmacol 59(4):566-577. doi:10.1002/jcph.1349.",
    "TRANSCRIPTION SOURCE: the full primary was not available; every value",
    "here is transcribed from Table 2 of the external evaluation",
    "Wang S, Yin Q, Yang M, Cheng Z, Xie F (2023). Pharmaceutics 15(2):569.",
    "doi:10.3390/pharmaceutics15020569. This file encodes only the ALL model;",
    "the primary also reports a separate osteosarcoma model, which Wang 2023",
    "did not evaluate and which is therefore not extracted here.",
    "Re-extract from the primary when it is obtained; see",
    "vignette('Wang_2023_methotrexate') Errata.",
    sep = " "
  )
  vignette <- "Wang_2023_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482. Wang 2023 Table 2 records this model as '2 CMT' with states V1
  # (central) and V2 (peripheral). Amount units are umol because every
  # methotrexate concentration in Wang 2023 is reported in umol/L.
  # verified = FALSE for the specimen: the Hui 2019 primary, which would state
  # the assayed matrix, is not on disk, and neither Wang 2023 Table 2 nor the
  # primary's abstract names it.
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BSA. Power-form effects on BOTH clearance and the",
        "central volume, each normalised to 0.735 m^2, per Wang 2023 Table 2:",
        "'CL (L/h) = 7.73 x (BSA/0.735)^0.721 x ...' and",
        "'V1 (L) = 19.0 x (BSA/0.735)^0.985'. BSA ALSO appears inside the",
        "renal-function term on clearance (see the CRCL entry). The 0.735 m^2",
        "reference is the cohort median -- Wang 2023 Table 1 records BSA",
        "median 0.74 m^2 (range 0.47-1.64) for this study."
      ),
      source_name        = "BSA"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. Power-form effect on clearance, per Wang 2023",
        "Table 2: 'CL (L/h) = 7.73 x (BSA/0.735)^0.721 x",
        "(eGFR x 1.73/192 x BSA)^0.256 x e^(eta_CL + IOV)'. Only two of the",
        "six models Wang evaluated include a renal-function descriptor on",
        "clearance (this one and Gao 2021), which Wang's Discussion attributes",
        "to most cohorts having normal renal function.",
        "PARENTHESISATION -- the bracket is encoded exactly as typeset,",
        "i.e. (eGFR x 1.73 x BSA / 192), reading the inline solidus with the",
        "usual left-to-right precedence. A 400 dpi render of Wang 2023 Table 2",
        "confirms an inline '/' rather than a fraction bar, so no denominator",
        "grouping is implied by the typesetting. The reading is also the one",
        "that behaves like a normalised covariate: at the evaluation cohort's",
        "median eGFR of 149.9 mL/min/1.73 m^2 (Wang 2023 Results 3.2) and this",
        "cohort's median BSA of 0.735 m^2, the bracket evaluates to 0.993, so",
        "the term is centred at approximately 1 as a normalised term should",
        "be. The alternative grouping (eGFR x 1.73)/(192 x BSA) -- which would",
        "read eGFR as an absolute mL/min value and 192 as a reference",
        "normalised eGFR -- cannot be excluded without the primary, and would",
        "raise typical clearance by about 17%. Flagged in the vignette Errata;",
        "the primary is queued for re-extraction."
      ),
      source_name        = "eGFR"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column Age. Power-form effect on intercompartmental clearance",
        "only, normalised to 5.29 years, per Wang 2023 Table 2:",
        "'Q (L/h) = 0.283 x (Age/5.29)^0.278'. The 5.29-year reference is the",
        "cohort median -- Wang 2023 Table 1 records age median 5.3 years",
        "(range 1.3-15.8) for this study. Age does NOT enter clearance or",
        "either volume in this model."
      ),
      source_name        = "Age"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator (one high-dose methotrexate course per occasion)",
      units              = "(integer)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Time-varying between courses, constant within a course. Wang 2023",
        "Table 2 records this model's variability as 'CL = 14.3 (14.9)' under",
        "the header 'IIV (%) (IOV (%))', i.e. a 14.3% between-subject CV and a",
        "14.9% between-occasion CV on clearance, and writes the exponential",
        "term as 'e^(eta_CL + IOV)'. Decomposed inside model() into binary",
        "indicators that multiplex per-occasion IOV etas onto log-clearance,",
        "following the Barnett_2018_coproporphyrin_I.R idiom. Wang 2023 does",
        "NOT state how many occasions the primary modelled; four slots are",
        "provided here as an implementation choice. This invents no parameter",
        "value -- Wang reports a SINGLE shared IOV variance, so every slot",
        "carries the same variance and only the first is estimable; the slot",
        "count merely bounds how many distinct courses a user can encode.",
        "Set OCC = 1 for a single-course simulation."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    age_range      = "1.3 to 15.8 years",
    age_median     = "5.3 years",
    weight_range   = "10.4 to 57.8 kg",
    weight_median  = "18.4 kg",
    height_range   = "77 to 176 cm (median 107)",
    bsa_range      = "0.47 to 1.64 m^2 (median 0.74)",
    sex_female_pct = 100 * 13 / 36,
    disease_state  = "Childhood acute lymphoblastic leukaemia (ALL) receiving high-dose methotrexate. The primary also studied 16 osteosarcoma patients under a separate model, not extracted here.",
    renal_function = "Serum creatinine median 0.362 mg/dL (range 0.11-1.07); eGFR is a model covariate but its cohort distribution is not reproduced in Wang 2023.",
    hepatic_function = "ALT median 16.5 U/L (range 5.0-247.0).",
    dose_range     = "2 to 5 g/m^2 intravenous high-dose methotrexate.",
    regions        = "Hong Kong / China.",
    notes          = paste(
      "Demographics from Wang 2023 Table 1; the primary's abstract records",
      "'36 ALL (354 observations) and 16 osteosarcoma (585 observations)",
      "patients' with covariate model building and parameter estimation in",
      "NONMEM and Perl-speaks-NONMEM, and the paper's purpose was to",
      "demonstrate popPK-based individual dose optimisation in R and shiny.",
      "This is the smallest of the six cohorts Wang 2023 evaluated. Wang's",
      "Results notes that this model, like Gao 2021, showed 'obvious trends'",
      "in the Bland-Altman prediction-error plots despite acceptable summary",
      "bias, and 'significant underprediction of peak MTX concentrations'.",
      "Wang's Results also mis-reads this row once, describing the 14.9% as an",
      "inter-individual effect ('implemented an inter-individual variability",
      "on MTX clearance with an estimated effect of 14.9%'); Table 2's own",
      "column header makes 14.9% the inter-OCCASION value and 14.3% the",
      "inter-individual one, and Table 2 is followed here."
    )
  )

  ini({
    # Structural parameters -- Wang 2023 Table 2, Hui row, for a typical child
    # of BSA 0.735 m^2 and age 5.29 years at the reference renal function.
    # No uncertainty (RSE, CI or SE) is reproduced in Wang 2023, so none is
    # recorded here.
    lcl <- log(7.73); label("Clearance at BSA 0.735 m^2 and reference renal function (L/h)")  # Wang 2023 Table 2: 'CL (L/h) = 7.73 x (BSA/0.735)^0.721 x (eGFR x 1.73/192 x BSA)^0.256'
    lvc <- log(19.0); label("Central volume of distribution at BSA 0.735 m^2 (L)")            # Wang 2023 Table 2: 'V1 (L) = 19.0 x (BSA/0.735)^0.985'
    lq  <- log(0.283); label("Intercompartmental clearance at age 5.29 years (L/h)")          # Wang 2023 Table 2: 'Q (L/h) = 0.283 x (Age/5.29)^0.278'
    lvp <- log(6.63); label("Peripheral volume of distribution (L)")                          # Wang 2023 Table 2: 'V2 (L) = 6.63 x e^eta_V2'

    # Covariate effects, all power exponents.
    e_bsa_cl  <- 0.721; label("Power exponent of body surface area on CL (unitless)")           # Wang 2023 Table 2: exponent in '(BSA/0.735)^0.721'
    e_crcl_cl <- 0.256; label("Power exponent of BSA-normalized eGFR on CL (unitless)")         # Wang 2023 Table 2: exponent in '(eGFR x 1.73/192 x BSA)^0.256'
    e_bsa_vc  <- 0.985; label("Power exponent of body surface area on Vc (unitless)")           # Wang 2023 Table 2: exponent in '(BSA/0.735)^0.985'
    e_age_q   <- 0.278; label("Power exponent of age on Q (unitless)")                          # Wang 2023 Table 2: exponent in '(Age/5.29)^0.278'

    # Between-subject variability, exponential. Wang 2023 Table 2 'IIV (%)'
    # column; converted with the log-normal identity omega^2 = log(1 + CV^2).
    # Neither Vc nor Q carries an eta in the published equations and none is
    # added here.
    etalcl ~ 0.020243  # Wang 2023 Table 2 IIV CL = 14.3% -> log(1 + 0.143^2)
    etalvp ~ 0.113075  # Wang 2023 Table 2 IIV V2 = 34.6% -> log(1 + 0.346^2)

    # Between-occasion variability on clearance. Wang 2023 Table 2 reports a
    # SINGLE shared IOV CV of 14.9% (the parenthesised value in
    # 'CL = 14.3 (14.9)'), so the same variance is carried on every occasion
    # slot and only the first is estimable -- the fix()ed siblings mirror the
    # NONMEM 'OMEGA BLOCK(1) SAME' idiom, as in
    # Barnett_2018_coproporphyrin_I.R. The number of slots (four) is an
    # implementation choice, not a published quantity; see covariateData$OCC.
    etaiov_cl_1 ~ 0.021958        # Wang 2023 Table 2 IOV CL = 14.9% -> log(1 + 0.149^2)
    etaiov_cl_2 ~ fix(0.021958)   # same shared IOV variance
    etaiov_cl_3 ~ fix(0.021958)   # same shared IOV variance
    etaiov_cl_4 ~ fix(0.021958)   # same shared IOV variance

    # Residual error.
    propSd <- 0.302; label("Proportional residual error (fraction)")  # Wang 2023 Table 2 RUV: 'Prop = 30.2%'
  })

  model({
    # Decompose the integer-valued OCC column into binary indicators and
    # multiplex the per-occasion IOV etas onto log-clearance.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iovCl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 +
      oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # Individual PK parameters. The renal term is written exactly as typeset in
    # Wang 2023 Table 2 -- (eGFR x 1.73 / 192 x BSA), i.e. eGFR scaled by 1.73
    # and by BSA and divided by 192. See covariateData$CRCL for why this
    # grouping is used and what the alternative would be.
    cl <- exp(lcl + etalcl + iovCl) *
      (BSA / 0.735)^e_bsa_cl *
      (CRCL * 1.73 * BSA / 192)^e_crcl_cl
    vc <- exp(lvc) * (BSA / 0.735)^e_bsa_vc
    q  <- exp(lq)  * (AGE / 5.29)^e_age_q
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
