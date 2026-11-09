Sato_2025_midazolam_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment population PK model for intravenous",
    "bolus midazolam in anaesthetised male Sprague-Dawley rats with or",
    "without extracorporeal membrane oxygenation (ECMO)",
    "(Sato 2025). The only retained covariate is a binary ECMO indicator",
    "that multiplies the central volume 3.02-fold; clearance, the",
    "peripheral volume and the intercompartmental clearance are shared by",
    "both groups. Parameters are absolute (L, L/h) for 400-450 g rats, not",
    "per kilogram."
  )
  reference <- paste(
    "Sato Y, Motoishi E, Kakuba R, Nagatsuka Y, Abe T, Fujii Y.",
    "Pharmacokinetics of midazolam and hepatic cytochrome P450 3A activity",
    "in an in vivo extracorporeal membrane oxygenation rat model.",
    "Sci Rep. 2025;15:19363. doi:10.1038/s41598-025-03949-x.",
    sep = " "
  )
  vignette <- "Sato_2025_midazolam_rat"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    ECMO_STATUS = list(
      description = "Binary indicator: 1 = rat on veno-arterial ECMO, 0 = sham-operated control",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source covariate 'dic (ECMO)' (Methods, PK analysis): a",
        "dichotomous 0/1 subject-level indicator. The circuit drained the",
        "right internal jugular vein and used the left common carotid",
        "artery cannula as its outflow (a veno-arterial configuration),",
        "with a membrane oxygenator and roller pump primed with 7 mL saline",
        "plus 1 mL heparin (Supplementary Information 2); auxiliary flow",
        "60-70 mL/kg/min (Methods, Experimental design).",
        "Midazolam was given 2 min BEFORE ECMO started, and the first",
        "sample (2 min) was drawn before ECMO initiation; the paper treats",
        "the indicator as constant per rat, and so does this model.",
        "Effect: V1 = theta_V1 * 3.02^ECMO_STATUS."
      ),
      source_name = "dic (ECMO)"
    )
  )

  compartmentData <- list(
    central = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "midazolam", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "14-16 weeks",
    weight_range = "400-450 g",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy male rats under isoflurane anaesthesia and mechanical",
      "ventilation; 3 on veno-arterial ECMO and 3 sham-operated controls"
    ),
    dose_range = "Midazolam 400 ug/kg single IV bolus via the jugular vein",
    regions = "Japan (Fukuyama University; Niigata University of Health and Welfare)",
    notes = paste(
      "Methods, Experimental design: 12 male Sprague-Dawley rats were",
      "studied, of which 6 contributed the midazolam PK data (3 ECMO, 3",
      "sham); 48 plasma samples at 2 min (pre-ECMO) and 5-180 min after",
      "ECMO initiation (Results, PK analysis). LC-MS/MS assay, linear range",
      "10-1000 ng/mL. NONMEM 7.5.1 with PsN; bootstrap (500) and VPC (500)."
    )
  )

  ini({
    # Typical values, Table 3 (final model, 'Mean' column). Absolute
    # volumes and clearances for the 400-450 g rats; no weight scaling.
    lcl <- log(0.69); label("Clearance CL (L/h)") # Table 3: theta CL = 0.69 L/hr (RSE 18.0%)
    lvc <- log(0.22); label("Central volume V1 in sham rats (L)") # Table 3: theta V1 = 0.22 L (RSE 41.2%)
    lvp <- log(0.49); label("Peripheral volume V2 (L)") # Table 3: theta V2 = 0.49 L (RSE 34.4%)
    lq <- log(0.29); label("Intercompartmental clearance Q (L/h)") # Table 3: theta Q = 0.29 L/hr (RSE 20.5%)

    # Power-of-indicator ECMO multiplier on V1: the printed equation
    # theta_i = theta_pop * dic * exp(eta) read as theta_pop * dic^ECMO.
    # Discussion: 'ECMO induced an increase in V1 (0.66 vs. 0.22 L)'.
    e_ecmo_status_vc <- 3.02; label("Multiplicative factor on V1 for ECMO_STATUS = 1 (unitless)") # Table 3: dic (ECMO) = 3.02 (RSE 41.7%)

    # IIV reported as CV%; omega^2 = log(CV^2 + 1) for the exponential model.
    etalcl ~ 0.617225 # Table 3: omega CL = 92.4 CV% -> log(0.924^2 + 1)
    etalvc ~ 0.215983 # Table 3: omega V1 = 49.1 CV% -> log(0.491^2 + 1)

    propSd <- 0.395; label("Proportional residual error (fraction)") # Table 3: Proportional = 39.5% (RSE 25.8%)
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * e_ecmo_status_vc^ECMO_STATUS
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in ug, volumes in L -> ug/L = ng/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
