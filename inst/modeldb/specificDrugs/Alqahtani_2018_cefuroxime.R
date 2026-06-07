Alqahtani_2018_cefuroxime <- function() {
  description <- "Two-compartment IV population PK model for cefuroxime in adults undergoing coronary artery bypass graft (CABG) surgery with cardiopulmonary bypass (Alqahtani 2018), with a power-form creatinine-clearance (Cockcroft-Gault) effect on clearance."
  reference <- paste(
    "Alqahtani SA, Alsultan AS, Alqattan HM, Eldemerdash A, Albacker TB.",
    "Population pharmacokinetic model-based evaluation of standard dosing",
    "regimens for cefuroxime used in coronary artery bypass graft surgery",
    "with cardiopulmonary bypass.",
    "Antimicrob Agents Chemother. 2018;62(6):e02241-17.",
    "doi:10.1128/AAC.02241-17.",
    sep = " "
  )
  vignette <- "Alqahtani_2018_cefuroxime"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw mL/min, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CL_CR estimated by the Cockcroft-Gault equation (Alqahtani 2018 Methods, Study design and settings). Raw mL/min, NOT BSA-normalized to mL/min/1.73 m^2. Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form -- same convention as Delattre_2010_amikacin and Shi_2018_ceftazidime). Reference value 78.5 mL/min (cohort mean, Alqahtani 2018 Table 1 and the Table 2 footnote b CL covariate formula). The effect on CL is a power form: CL = 3.43 * (CRCL / 78.5)^0.56.",
      source_name        = "CL_CR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 78L,
    n_studies      = 1L,
    age_range      = "18-80 years",
    age_mean       = "54.2 years (SD 13.2)",
    weight_range   = "41-111.8 kg",
    weight_mean    = "76.7 kg (SD 14.7)",
    bmi_range      = "16.6-43.2 kg/m^2 (mean 28.6, SD 5.2)",
    sex_female_pct = 24,
    race_ethnicity = "Not reported (single-centre Saudi Arabian cohort; King Fahad Cardiac Center, King Saud University Medical City, Riyadh)",
    disease_state  = "Adults scheduled to undergo cardiac surgical procedures (coronary artery bypass graft surgery with cardiopulmonary bypass); patients allergic to beta-lactams, with prior systemic infections, or who received antibiotic therapy in the 72 h before surgery were excluded",
    dose_range     = "Cefuroxime 1.5 g IV infusion over 30 min administered 30-60 min before skin incision; an extra 1.5 g dose was mixed into the CPB solution if surgery lasted >4 h; subsequent prophylactic doses for 48 h were either 1.5 g every 12 h or 1 g every 8 h",
    regions        = "Saudi Arabia (single-centre prospective open-label study at King Fahad Cardiac Center, King Saud University Medical City)",
    renal_function = "Cockcroft-Gault creatinine clearance mean 78.5 mL/min (SD 24, range 28.3-125); serum creatinine mean 85 umol/L (SD 29.7, range 41-245); raw mL/min, not BSA-normalized",
    notes          = "Baseline demographics per Alqahtani 2018 Table 1. 78 adults; 468 plasma samples analyzed by validated HPLC (linearity 0.5-200 ug/mL, equivalent to 0.5-200 mg/L). Cefuroxime PK fit using Monolix v4.4 (SAEM). Eight covariates were screened (age, weight, serum creatinine, CL_CR, gender, height, albumin, BMI); only CL_CR was retained on CL after stepwise log-likelihood ratio testing. Indication was antibiotic prophylaxis against postoperative surgical-site infection."
  )

  ini({
    # Structural parameters at the reference subject (CRCL = 78.5 mL/min);
    # Alqahtani 2018 Table 2 final-model column.
    lcl <- log(3.43); label("Clearance at 78.5 mL/min CRCL (L/h)")              # Alqahtani 2018 Table 2: CL = 3.43 L/h
    lvc <- log(3.88); label("Central volume of distribution V1 (L)")             # Alqahtani 2018 Table 2: V1 = 3.88 L
    lq  <- log(22.2); label("Intercompartmental clearance Q (L/h)")              # Alqahtani 2018 Table 2: Q  = 22.2 L/h
    lvp <- log(5.7);  label("Peripheral volume of distribution V2 (L)")          # Alqahtani 2018 Table 2: V2 = 5.7 L

    # Covariate effect: power-form CRCL on CL,
    #   CL = 3.43 * (CRCL / 78.5)^0.56
    # (Alqahtani 2018 Table 2 footnote b).
    e_crcl_cl <- 0.56; label("CRCL exponent on CL (unitless)")                   # Alqahtani 2018 Table 2 footnote b: CL = 3.43*(CL_CR/78.5)^0.56

    # Inter-individual variability (Alqahtani 2018 Table 2 final-model CV%);
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.18590 # log(0.452^2 + 1); 45.2% CV on CL  (Alqahtani 2018 Table 2)
    etalvc ~ 0.31901 # log(0.613^2 + 1); 61.3% CV on V1  (Alqahtani 2018 Table 2)
    etalq  ~ 0.27448 # log(0.562^2 + 1); 56.2% CV on Q   (Alqahtani 2018 Table 2)
    etalvp ~ 0.04643 # log(0.218^2 + 1); 21.8% CV on V2  (Alqahtani 2018 Table 2)

    # Combined residual error (Alqahtani 2018 Table 2 final-model rows a and b).
    # Monolix combined-error convention: SD(eps) = a + b*pred, with
    #   a = additive SD in concentration units (mg/L; cefuroxime HPLC range
    #       0.5-200 ug/mL == mg/L per Methods, Analytical method), and
    #   b = proportional fraction (unitless). This matches the encoding used in
    #   Delattre_2010_amikacin and Shi_2018_ceftazidime for the same Monolix /
    #   NONMEM "combined additive + proportional" form.
    addSd  <- 6.45;  label("Additive residual error (mg/L)")          # Alqahtani 2018 Table 2: a = 6.45
    propSd <- 0.092; label("Proportional residual error (fraction)")  # Alqahtani 2018 Table 2: b = 0.092
  })
  model({
    # Individual PK parameters. CL has a power-form CRCL effect with reference
    # 78.5 mL/min; V1, Q and V2 have no retained covariates.
    cl <- exp(lcl + etalcl) * (CRCL / 78.5)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L (equivalent to ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
