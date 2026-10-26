Edwina_2023_fosfomycin <- function() {
  description <- paste(
    "Two-compartment IV population PK model for fosfomycin in healthy adult",
    "men, fit jointly to plasma and urine concentrations after intermittent",
    "(8 g q8h) and continuous (8 g loading + 1 g/h) infusion (Edwina 2023).",
    "All elimination is renal: the clearance from the central compartment",
    "feeds a urine compartment and is a power function of absolute",
    "(not BSA-indexed) CKD-EPI estimated GFR (reference 120 mL/min). The",
    "urine compartment is read as a urinary concentration over the urine",
    "volume collected in each interval."
  )
  reference <- paste(
    "Edwina AE, Koch BCP, Muller AE, al Jalali V, Matzneller P,",
    "Zeitlinger M, Sassen SDT.",
    "Population plasma and urine pharmacokinetics and the probability of",
    "target attainment of fosfomycin in healthy male volunteers.",
    "Eur J Clin Pharmacol. 2023;79(6):775.",
    "doi:10.1007/s00228-023-03477-5."
  )
  vignette <- "Edwina_2023_fosfomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "fosfomycin", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate by the CKD-EPI creatinine",
        "equation, in ABSOLUTE mL/min (NOT indexed to 1.73 m^2)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column GFR (supplementary NONMEM control stream). Methods",
        "list 'eGFR ... calculated with the CKD-EPI in millilitres per",
        "minute' and, separately, 'eGFR calculated with the CKD-EPI formula",
        "normalised to 1.73 m2'; Supplementary Table 4 retains the",
        "un-normalised form (model 11, dOFV -12.909) over the BSA-normalised",
        "one (model 12, dOFV -4.186). Enters clearance as (CRCL/120)^0.774",
        "(control stream COVCLRGFR = (GFR/120)**THETA(7)); cohort median",
        "119.63 mL/min, range 89-138 (Table 1)."
      ),
      source_name = "GFR"
    ),
    URINE_VOL_INTERVAL = list(
      description = "Urine volume collected in the urine collection interval containing the current urinary observation",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Control stream: S3 = UVL with UVL = UVOL, so the urinary prediction",
        "is the amount in the urine compartment divided by the collected",
        "volume. The source does not state the UVOL units; stored here in mL",
        "and divided by 1000 in model() so that the mg amount yields mg/L,",
        "the assay's reporting unit. Reset the urine state at each",
        "collection-interval boundary (evid = 5, amt = 0 placed just after",
        "the boundary). Only the urine observation uses it; it has no effect",
        "on plasma concentrations."
      ),
      source_name = "UVOL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Supplementary Table 4 model 7, dOFV -0.985); not retained.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Supplementary Table 4 model 8, dOFV -2.369); not retained.",
      source_name = "WEIGHT"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Supplementary Table 4 model 9, dOFV -2.989); not retained.",
      source_name = "BMI"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on CL (Supplementary Table 4 model 6, dOFV -7.456; forward",
        "step model 15, dOFV -3.764); not retained. Table 1 labels the unit",
        "g/dL but the values (median 44) are g/L."
      ),
      source_name = "ALB"
    ),
    TPRO = list(
      description = "Serum total protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on CL (Supplementary Table 4 model 10, dOFV -7.471); added",
        "in forward inclusion (model 16) but removed at backward elimination",
        "(model 17, dOFV +6.124 < 6.63 at p < 0.01). Table 1 labels the unit",
        "g/dL but the values (median 69.45) are g/L."
      ),
      source_name = "PRO"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL (Supplementary Table 4 model 13, dOFV -1.551); not retained.",
      source_name = "SCR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8,
    n_studies = 1,
    age_range = "23-42 years (median 36.5)",
    weight_range = "53-93 kg (median 83.5)",
    sex_female_pct = 0,
    disease_state = "Healthy volunteers",
    dose_range = paste(
      "Crossover: 8 g IV over 30 min every 8 h for three doses, and an 8 g",
      "30-min loading dose followed by 1 g/h continuous infusion for 18 h;",
      "washout of at least 48 h between regimens."
    ),
    regions = "Austria (Medical University of Vienna)",
    renal_function = "CKD-EPI eGFR 89-138 mL/min (median 119.63); Cockcroft-Gault CrCL 97-151 mL/min (Table 1)",
    notes = paste(
      "Table 1 and Methods: eight healthy men, 176 plasma and 86 urine",
      "samples. Urine collected over 0-3, 3-8, 8-16, 16-18, 18-20, 20-22",
      "and 22-24 h. Urine amounts exceeding the administered dose were",
      "truncated to the dose before fitting. NONMEM 7.4, FOCE-I, ADVAN6."
    )
  )

  ini({
    # Structural parameters -- Edwina 2023 Table 2 (final model).
    lcl <- log(6.34)
    label("Renal (= total) clearance at CRCL = 120 mL/min (L/h)") # Table 2: CLpop = 6.34 L/h (RSE 3.3%)
    lvc <- log(10.2)
    label("Central volume of distribution (L)") # Table 2: Vc = 10.2 L (RSE 7.7%)
    lq <- log(11.8)
    label("Intercompartmental clearance (L/h)") # Table 2: Q = 11.8 L/h (RSE 10.9%)
    lvp <- log(9.74)
    label("Peripheral volume of distribution (L)") # Table 2: Vp = 9.74 L (RSE 6.2%)

    # Covariate effect -- Table 2, CL = CLpop * (eGFR/120)^beta.
    e_crcl_cl <- 0.774
    label("Exponent of CRCL (absolute CKD-EPI eGFR) on clearance (unitless)") # Table 2: beta = 0.774 (RSE 14.2%)

    # IIV -- Table 2 'Interindividual variability of CL (%CV)' = 4.8;
    # omega^2 = log(CV^2 + 1) = log(0.048^2 + 1).
    etalcl ~ 0.0023014 # Table 2: IIV CL 4.8 %CV (shrinkage 4%)

    # Residual error -- control stream W = IPRED * THETA with SIGMA 1 FIX,
    # so the THETA is the proportional SD.
    propSd <- 0.3817
    label("Proportional residual SD, plasma (fraction)") # Table 2: proportional error in plasma = 38.17 %CV
    propSd_Curine <- 0.8715
    label("Proportional residual SD, urine (fraction)") # Table 2: proportional error in urine = 87.15 %CV
  })

  model({
    # Supplementary control stream $PK: CLR = THETA(3) * (GFR/120)**THETA(7)
    # * EXP(ETA(1)); K13 = CLR/V1 is the only elimination route.
    cl <- exp(lcl + etalcl) * (CRCL / 120)^e_crcl_cl
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(urine) <- kel * central

    Cc <- central / vc
    urine_volume <- max(URINE_VOL_INTERVAL / 1000, 0.001)
    Curine <- urine / urine_volume

    Cc ~ prop(propSd)
    Curine ~ prop(propSd_Curine)
  })
}
