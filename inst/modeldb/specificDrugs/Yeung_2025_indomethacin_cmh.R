Yeung_2025_indomethacin_cmh <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous indomethacin in",
    "preterm neonates (CYP2C9-reduced IND popPK model, re-estimated on the",
    "Children's Mercy Hospital development dataset). Total clearance is",
    "split into renal and non-renal (presumably hepatic) clearance, each a",
    "power function of time-varying postnatal age centred on 1.25 days;",
    "the renal arm feeds a cumulative urine-amount state. Plasma and dried",
    "blood spot (DBS, converted to plasma-equivalent) concentrations share",
    "the same prediction with matrix-specific log-scale residual errors.",
    "Structure of Krzyzanski 2023 with the CYP2C9 covariate on renal",
    "clearance removed.",
    sep = " "
  )
  reference <- paste(
    "Yeung CHT, Sekulich DC, Scott A, Nolte WM, Gibson K, Su R,",
    "Alrifai MW, Lopata SM, Lewis T, Reese J, Maharaj A (2025).",
    "The Relationship of Indomethacin Exposure With Efficacy and Renal",
    "Toxicity Outcomes for Preterm Infants in the Neonatal Intensive Care",
    "Unit. Clin Transl Sci 18(5):e70251. doi:10.1111/cts.70251.",
    "Structural model from Krzyzanski W, Stockard B, Gaedigk A, et al.",
    "(2023) Developmental pharmacokinetics of indomethacin in preterm",
    "neonates: severely decreased drug clearance in the first week of",
    "life. CPT Pharmacometrics Syst Pharmacol 12:110-121.",
    "doi:10.1002/psp4.12881.",
    sep = " "
  )
  vignette <- "Yeung_2025_indomethacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    PNA = list(
      description = "Postnatal age, time-varying (updated at every record).",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Krzyzanski 2023 (the structural model Yeung 2025 retained) enters",
        "PNA in DAYS as a power term (PNA / PNAmean)^theta with PNAmean =",
        "1.25 days, the mean postnatal age at the first dose of the",
        "development cohort, and updates PNA at every observation. The",
        "canonical PNA column is in months, so model() recovers days with",
        "PNA * 30.4375. Supply PNA > 0: the renal exponent is negative,",
        "so the renal clearance diverges at PNA = 0."
      ),
      source_name = "PNA (days)"
    ),
    SAMPLE_DBS = list(
      description = "Per-observation sampling-matrix indicator (1 = concentration from a dried blood spot converted to plasma; 0 = plasma).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (plasma)",
      notes = paste(
        "Selects the matrix-specific log-scale residual error only; the",
        "DBS value was converted to a plasma-equivalent concentration",
        "before modelling (Yeung 2025 Equation 1), so the structural",
        "prediction is the same plasma concentration for both matrices.",
        "Has no effect on urine records."
      ),
      source_name = "DBS"
    )
  )

  compartmentData <- list(
    central = list(analyte = "indomethacin", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "indomethacin", units = "mg", specimen = "urine", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 38,
    n_studies = 1,
    age_range = "postnatal age at first dose 0.23 h to 53.7 days across the CMH cohorts (Yeung 2025 Table 1)",
    weight_range = "360-2630 g at first dose across the CMH IVH and PDA cohorts (Yeung 2025 Table 1)",
    ga_range = "22.6-34 weeks gestational age at birth (Yeung 2025 Table 1)",
    disease_state = "Preterm neonates receiving IV indomethacin for intraventricular hemorrhage (IVH) prophylaxis or patent ductus arteriosus (PDA) closure",
    dose_range = "0.1 mg/kg IV q24h x3 (IVH prophylaxis) or 0.1 mg/kg IV q12h x3 per course (PDA, up to 3 courses, later courses could be increased); CMH median 0.107 (0.085-0.40) mg/kg",
    regions = "USA (Children's Mercy Hospital, Kansas City)",
    notes = paste(
      "Re-estimation of the Krzyzanski 2023 model on the CMH development",
      "dataset with the CYP2C9 genotype covariate removed (Yeung 2025",
      "Methods 'PopPK Modeling', Table S2). Yeung 2025 does not restate",
      "the subject count of that dataset; 38 is the development cohort of",
      "Krzyzanski 2023 Table 1 (median weight 779 g, GA 26.1 weeks, PNA",
      "at first dose 0.26 days). Plasma, DBS and 12-h pooled urine samples."
    )
  )

  ini({
    lvc <- log(0.435); label("Volume of distribution (L)") # Table S2, CMH 'theta V [L]' = 0.435 (RSE 7.2%)
    lcl_renal <- log(0.0000172); label("Renal clearance at PNA = 1.25 days (L/h)") # Table S2, CMH 'theta CLR [L/h]' = 0.0000172 (RSE 11.6%)
    lcl_nonren <- log(0.000647); label("Non-renal (hepatic) clearance at PNA = 1.25 days (L/h)") # Table S2, CMH 'theta CLH [L/h]' = 0.000647 (RSE 63.5%)
    e_pna_cl_renal <- -0.192; label("Power exponent of PNA on renal clearance (unitless)") # Table S2, CMH 'theta PNA-CLR' = -0.192
    e_pna_cl_nonren <- 0.901; label("Power exponent of PNA on non-renal clearance (unitless)") # Table S2, CMH 'theta PNA-CLH' = 0.901 (RSE 25.4%)

    etalvc ~ log(1 + 0.385^2) # Table S2, CMH 'IIV V (%)' = 38.5
    etalcl_renal ~ log(1 + 0.690^2) # Table S2, CMH 'IIV CLR (%)' = 69.0
    etalcl_nonren ~ log(1 + 2.154^2) # Table S2, CMH 'IIV CLH (%)' = 215.4

    expSd <- 0.312730; label("Log-scale residual SD, plasma concentrations (unitless)") # Table S2, CMH 'epsilon2 log-normal plasma' = 0.0978; sqrt(0.0978)
    expSdDbs <- 0.521536; label("Log-scale residual SD, dried blood spot concentrations (unitless)") # Table S2, CMH 'epsilon2 log-normal DBS' = 0.272; sqrt(0.272)
    expSd_Aurine <- 0.310483; label("Log-scale residual SD, cumulative urine amount (unitless)") # Table S2, CMH 'epsilon2 log-normal urine' = 0.0964; sqrt(0.0964)
  })

  model({
    # PNA power terms centred on the development-cohort mean PNA at first
    # dose, 1.25 days (Krzyzanski 2023 Equation 5 and Results).
    pna_days <- PNA * 30.4375

    vc <- exp(lvc + etalvc)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (pna_days / 1.25)^e_pna_cl_renal
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * (pna_days / 1.25)^e_pna_cl_nonren

    # Krzyzanski 2023 Equations 1-3: IV bolus into the plasma compartment;
    # the renal outflow accumulates in urine.
    d/dt(central) <- -(cl_renal + cl_nonren) / vc * central
    d/dt(urine) <- cl_renal / vc * central

    Cc <- central / vc
    Aurine <- urine

    # Krzyzanski 2023 Equations 7-8: constant residual error on the log
    # scale, with a separate variance for DBS-derived concentrations.
    expSdCc <- expSd * (1 - SAMPLE_DBS) + expSdDbs * SAMPLE_DBS
    Cc ~ lnorm(expSdCc)
    Aurine ~ lnorm(expSd_Aurine)
  })
}
