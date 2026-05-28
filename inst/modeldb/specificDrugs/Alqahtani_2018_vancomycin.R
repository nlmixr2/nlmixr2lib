Alqahtani_2018_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin used as prophylactic antibiotic in 28 adult patients undergoing open heart surgery with cardiopulmonary bypass (Alqahtani 2018). Clearance scales by power exponent with Cockcroft-Gault creatinine clearance (raw mL/min, reference 83.5) and serum albumin (g/L, reference 35.5); central volume scales by power exponent with body weight (kg, reference 79.6)."
  reference <- "Alqahtani SA, Alsultan AS, Alqattan HM, Eldemerdash A, Albacker TB. Population pharmacokinetic model for vancomycin used in open heart surgery: model-based evaluation of standard dosing regimens. Antimicrob Agents Chemother. 2018;62(7):e00088-18. doi:10.1128/AAC.00088-18"
  vignette <- "Alqahtani_2018_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Alqahtani 2018 Table 1: mean 79.6 kg (SD 17; range 52-111.8). Reference 79.6 kg used in the V1 power-scaling term in Table 3 footnote.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Alqahtani 2018 Methods: CLCR estimated using the Cockcroft-Gault equation. Table 1 cohort mean 83.5 mL/min (SD 29.3; range 33.4-125); reference value 83.5 mL/min used in the CL power-scaling term in Table 3 footnote (CL = 6.13 * (CLCR/83.5)^0.514 * (albumin/35.5)^0.854). Stored under canonical CRCL with raw mL/min (per inst/references/covariate-columns.md, CRCL accepts raw Cockcroft-Gault when the source paper does not apply BSA normalization).",
      source_name        = "CLCR"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Alqahtani 2018 Table 1: mean 35.5 g/L (SD 4.5; range 25-44). Reference 35.5 g/L used in the CL power-scaling term in Table 3 footnote.",
      source_name        = "albumin"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 28L,
    n_studies        = 1L,
    age_range        = "18-78 years",
    age_median       = "51.7 years (mean, SD 15.9)",
    weight_range     = "52-111.8 kg",
    weight_median    = "79.6 kg (mean, SD 17)",
    sex_female_pct   = 39,
    race_ethnicity   = "Not reported (single-center study in Riyadh, Saudi Arabia)",
    disease_state    = "Adults >=18 years undergoing open heart surgery with cardiopulmonary bypass; receiving vancomycin as prophylactic antibiotic. Patients excluded if vancomycin was administered in the 72 h before surgery.",
    dose_range       = "1 g vancomycin IV infusion over 30 min, first dose 2 h before skin incision then q12h for 48 h; extra (redosing) dose if surgery lasted >4 h (about 12 of 28 patients received an extra dose).",
    regions          = "Saudi Arabia (single center, King Fahad Cardiac Center, King Saud University Medical City, Riyadh)",
    renal_function   = "Cockcroft-Gault CLCR mean 83.5 mL/min (range 33.4-125)",
    n_concentrations = 168L,
    notes            = "Prospective open-label observational PK study (Alqahtani 2018 Table 1). 168 plasma vancomycin concentrations from 28 patients; 6 nominal sampling times per patient (before skin incision, at start of CPB, 1 h after starting CPB, before skin closure, 24 h and 48 h after first dose). Vancomycin assay: Architect iVancomycin chemiluminescent microparticle immunoassay (analytical range 0.24-100 ug/mL). Model fit using Monolix 4.4 SAEM. Body mass index mean 29.8 (range 20.2-41.9); serum creatinine mean 77.2 umol/L (range 41-134). The paper notes that CPB itself was tested as a covariate on CL and V and was not retained in the final model."
  )

  ini({
    # Structural parameters (Alqahtani 2018 Table 3). Reference subject is 79.6
    # kg, with CLCR = 83.5 mL/min and albumin = 35.5 g/L.
    lcl <- log(6.13); label("Clearance at CLCR=83.5 mL/min, albumin=35.5 g/L (L/h)")  # Alqahtani 2018 Table 3: CL = 6.13 L/h (RSE 19%)
    lvc <- log(40);   label("Central volume at WT=79.6 kg (L)")                       # Alqahtani 2018 Table 3: V1 = 40 L (RSE 15%)
    lq  <- log(0.22); label("Intercompartmental clearance Q (L/h)")                   # Alqahtani 2018 Table 3: Q = 0.22 L/h (RSE 10%)
    lvp <- log(3.88); label("Peripheral volume V2 (L)")                               # Alqahtani 2018 Table 3: V2 = 3.88 L (RSE 16%)

    # Covariate effects (Alqahtani 2018 Table 3 footnote):
    #   CL = 6.13 * (CLCR/83.5)^0.514 * (albumin/35.5)^0.854
    #   V1 = 40   * (WT/79.6)^0.466
    e_crcl_cl <- 0.514; label("Power exponent on (CRCL/83.5) for CL")    # Alqahtani 2018 Table 3 footnote
    e_alb_cl  <- 0.854; label("Power exponent on (ALB/35.5) for CL")     # Alqahtani 2018 Table 3 footnote
    e_wt_vc   <- 0.466; label("Power exponent on (WT/79.6) for V1")      # Alqahtani 2018 Table 3 footnote

    # Inter-individual variability (Alqahtani 2018 Table 3 IIV rows, %CV).
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.047686  # log(0.221^2 + 1); 22.1% CV on CL
    etalvc ~ 0.004012  # log(0.0634^2 + 1); 6.34% CV on V1
    etalq  ~ 0.288245  # log(0.578^2 + 1); 57.8% CV on Q
    etalvp ~ 0.318122  # log(0.612^2 + 1); 61.2% CV on V2

    # Combined additive + proportional residual error (Alqahtani 2018 Table 3).
    addSd  <- 0.055; label("Additive residual error (mg/L)")        # Alqahtani 2018 Table 3: a = 0.055 mg/L (RSE 11%)
    propSd <- 0.152; label("Proportional residual error (fraction)") # Alqahtani 2018 Table 3: b = 15.2% (RSE 7%)
  })
  model({
    # Individual PK parameters. CL scales by (CRCL/83.5)^0.514 and by
    # (ALB/35.5)^0.854; V1 scales by (WT/79.6)^0.466; Q and V2 have no
    # retained covariates.
    cl <- exp(lcl + etalcl) * (CRCL / 83.5)^e_crcl_cl * (ALB / 35.5)^e_alb_cl
    vc <- exp(lvc + etalvc) * (WT / 79.6)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
