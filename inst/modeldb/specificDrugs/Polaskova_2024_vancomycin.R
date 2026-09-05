Polaskova_2024_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin during the INITIAL phase of therapy (first 3 days) in 138 adult obese (BMI >= 30 kg/m2) patients (Polaskova 2024). Clearance is 1.32 L/h multiplied by two UNCENTERED exponential covariate terms, exp(0.61 x eGFR) and exp(0.011 x LBM), where eGFR is the creatinine-based CKD-EPI estimate in mL/s/1.73 m2 and LBM is Boer-formula lean body mass in kg; at the cohort median covariates (eGFR 1.51, LBM 68 kg) this gives 7 L/h. The central volume is 75.0 L with NO retained covariate: the authors screened body weight, LBM, BSA and BMI against Vd and found none of them a reliable predictor in this obese cohort. Residual error is additive (constant) at 2.9 mg/L."
  reference <- "Polaskova L, Murinova I, Gregorova J, Slanar O, Sima M. Vancomycin population pharmacokinetics and dosing proposal for the initial treatment in obese adult patients. Front Pharmacol. 2024;15:1364681. doi:10.3389/fphar.2024.1364681"
  vignette <- "Polaskova_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine-based CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Polaskova 2024 Methods 2.2: 'Both creatinine- and cystatin C-based (if available) glomerular",
        "filtration rates (eGFR) were estimated according to the Chronic Kidney Disease Epidemiology",
        "Collaboration (CKD-EPI) formula for each patient (Inker et al., 2012).' The creatinine-based",
        "variant is the one retained -- Results 3.2 and Figure 2B report that creatinine-based eGFR",
        "predicted vancomycin CL better than the cystatin C-based estimate in the n = 46 subgroup where",
        "both were measured, which the Discussion attributes to cystatin C being falsely elevated in",
        "obesity.",
        "UNITS: the paper reports eGFR throughout in mL/s/1.73 m^2 (SI, the Czech clinical convention),",
        "NOT in the canonical mL/min/1.73 m^2 of this register. Table 1 gives median 1.51, IQR 1.12-1.72,",
        "range 0.17-2.47 mL/s/1.73 m^2, i.e. median 90.6, IQR 67.2-103.2, range 10.2-148.2",
        "mL/min/1.73 m^2. This column carries the CANONICAL mL/min/1.73 m^2 value and model() divides it",
        "by 60 before applying the published coefficient 0.61, which is per (mL/s/1.73 m^2). A user who",
        "supplies mL/s/1.73 m^2 here will understate clearance by a large factor.",
        "The effect is EXPONENTIAL and UNCENTERED: exp(0.61 x eGFR_SI), so it is 1 at eGFR = 0 and 3.10",
        "at the cohort median 1.51 mL/s/1.73 m^2. The cohort skews preserved-to-augmented -- 52% of",
        "patients had eGFR >= 1.5 and 5% >= 2.13 mL/s/1.73 m^2 (>= 90 and >= 128 mL/min/1.73 m^2) --",
        "which is why the paper is able to propose dosing for augmented renal clearance. Patients on",
        "renal replacement therapy or extracorporeal life support were excluded, so the model carries no",
        "information about dialysis. The Discussion explicitly rejects Cockcroft-Gault for this cohort",
        "because it takes body weight as an input and therefore overestimates filtration in obesity.",
        collapse = " "
      ),
      source_name        = "eGFR"
    ),
    LBM = list(
      description        = "Lean body mass by the Boer formula",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Polaskova 2024 Methods 2.2: 'The body surface area (BSA) and lean body mass (LBM) were",
        "calculated using Du Bois and Boer formulas, respectively (Boer, 1984; Du Bois and Du Bois,",
        "1989).' The Boer formulae are LBM (male) = 0.407 x WT + 0.267 x HT - 19.2 and LBM (female) =",
        "0.252 x WT + 0.473 x HT - 48.3, with WT in kg and HT in cm; the register's LBM note warns that",
        "the Boer, James and Hume formulae differ by several kg at a given height and weight, so a",
        "downstream user must use Boer to reproduce this model's covariate distribution.",
        "Table 1: median 68 kg, IQR 55-76, range 41-104.",
        "The effect is EXPONENTIAL and UNCENTERED: exp(0.011 x LBM), so it is 1 at LBM = 0 and 2.11 at",
        "the cohort median 68 kg. LBM enters CL only; the paper screened LBM (with BW, BSA and BMI)",
        "against Vd and retained none of them.",
        collapse = " "
      ),
      source_name        = "LBM"
    )
  )

  # Covariates the paper screened but did not retain in the final model. No
  # point estimate was published for any of them, so none can be encoded. The
  # Vd screen is the notable one: it is a NEGATIVE result the paper argues for
  # at length rather than an omission (Results 3.2, Figure 2A, Discussion).
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight. Screened against both Vd and CL; not retained.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table 1: median 104 kg, IQR 95-120, range 73-190. Results 3.2: against Vd, 'the other body",
        "size descriptors (BW, LBM, and BSA) were found to be without statistical significance",
        "(Figure 2A)'; against CL, BW was one of several variables positively related in the",
        "preliminary graphical assessment but 'vancomycin CL was best predicted using eGFR and LBM'",
        "after stepwise covariate modelling. The Discussion defends the absent weight effect on Vd:",
        "'it is physiologically plausible that the Vd of hydrophilic compounds does not increase",
        "proportionally with body weight in obese patients, in whom weight gain is mainly due to the",
        "deposition of adipose tissue.'",
        collapse = " "
      )
    ),
    HT = list(
      description        = "Height. Screened as a continuous covariate on CL; not retained.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1: median 1.74 m, IQR 1.65-1.80, range 1.50-1.95 (reported in m; this column is cm). Results 3.2 lists height among the variables positively related to CL in the preliminary graphical assessment, but it was not retained."
    ),
    BSA = list(
      description        = "Body surface area by the Du Bois formula. Screened against both Vd and CL; not retained.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1: median 2.17 m^2, IQR 2.02-2.35, range 1.68-2.77. Methods 2.2 names the Du Bois formula. Screened against Vd (Figure 2A, not significant) and against CL (positively related graphically, not retained after stepwise modelling)."
    ),
    BMI = list(
      description        = "Body mass index. Screened against Vd; not retained.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1: median 34.3 kg/m^2, IQR 32.5-38.3, range 30.1-65.7 (BMI >= 30 was the inclusion criterion). Results 3.2: 'The preliminary graphical assessment showed only a very weak relationship between BMI and vancomycin Vd', and covariate diagnostics on the final model 'found that none of the covariates tested reliably predicted vancomycin Vd.' No point estimate published."
    ),
    AGE = list(
      description        = "Subject age. Screened as a continuous covariate; not retained.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1: median 65 years, IQR 54-72, range 26-86. Listed in Methods 2.3 among the continuous covariates tested. No point estimate published."
    ),
    CREAT = list(
      description        = "Serum creatinine. Screened as a continuous covariate on CL; not retained.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1: median 72.3 umol/L, IQR 55.8-98.1, range 23.2-374.1. Results 3.2 reports CL was negatively related to serum creatinine in the preliminary graphical assessment, but the CKD-EPI eGFR derived from it was retained instead. No point estimate published."
    ),
    BUN = list(
      description        = "Serum urea. Screened as a continuous covariate on CL; not retained.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Table 1: median 5.6 mmol/L, IQR 3.8-7.9, range 1.2-24.4. Reported as serum UREA (mmol/L), not blood urea nitrogen; the register admits both under BUN with the unit documented per model. Results 3.2 reports CL was negatively related to urea graphically; not retained. No point estimate published."
    ),
    CYSC = list(
      description        = "Serum cystatin C, measured in a subgroup only. Screened as the basis of an alternative eGFR; not retained.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Measured in n = 46 of 138 patients (Results 3.2). The cystatin C-based CKD-EPI eGFR was compared head-to-head against the creatinine-based estimate in that subgroup (Figure 2B) and lost. No value distribution and no point estimate are published."
    ),
    SEXF = list(
      description        = "Female sex indicator. Screened as a categorical covariate; not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Results 3.1: 82 males and 56 females (40.6% female). Methods 2.3 lists gender among the categorical covariates tested. No point estimate published. Note that sex is still an INPUT to this model indirectly, through the sex-specific Boer LBM formula."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 138L,
    n_studies      = 1L,
    n_observations = 147L,
    age_range      = "26-86 years",
    age_median     = "65 years",
    weight_range   = "73-190 kg",
    weight_median  = "104 kg",
    sex_female_pct = 40.6,
    race_ethnicity = NULL,
    disease_state  = "Adult obese (BMI >= 30 kg/m2) inpatients treated for suspected or proven invasive Gram-positive infection: CNS infection 33%, sepsis 18%, orthopaedic 15%, ocular 9%, skin 6%, other (pneumonia, bacteriuria, bacteraemia, endocarditis, intra-abdominal) 19%",
    dose_range     = "Loading dose 1-4 g (median 2.5 g) as a 0.5-7 h (median 5 h) IV infusion in 122 of 138 patients; maintenance dose 0.5-1.5 g (median 1 g) q6h, q8h, q12h or q24h as a 1-3 h (median 2 h) IV infusion, or 1-4 g/day (median 2 g/day) by continuous infusion in 9 patients",
    regions        = "Czechia (single centre: Military University Hospital Prague)",
    renal_function = "eGFR (CKD-EPI, creatinine) median 1.51 mL/s/1.73 m2 (IQR 1.12-1.72, range 0.17-2.47), i.e. median 90.6 mL/min/1.73 m2 (IQR 67.2-103.2, range 10.2-148.2). 52% of patients had eGFR >= 1.5 mL/s/1.73 m2 and 5% >= 2.13. Patients on renal replacement therapy or extracorporeal life support were EXCLUDED.",
    notes          = paste(
      "Retrospective open-label observational study of routine therapeutic drug monitoring data,",
      "January 2013 to December 2022 (Methods 2.1). Baseline demographics are Polaskova 2024 Table 1.",
      "Only concentrations from the INITIAL phase of therapy -- the first 3 days of treatment -- were",
      "included, which is the paper's whole point: most published vancomycin popPK models describe",
      "maintenance dosing at steady state. 147 serum concentrations from 138 patients (1-2 per",
      "patient): 11 (7.5%) peaks taken up to 2 h after the end of infusion, 124 (84.4%) troughs taken",
      "0-1 h before the next dose, and 12 (8.2%) mid-interval samples with recorded sampling times.",
      "Assay: immunoturbidimetric KIMS on a Roche Cobas 8000, LLOQ 4.0 mg/L, measuring range 4.0-80.0",
      "mg/L. Estimation was by SAEM in Monolix 2021R2. Model stability was checked by a 250-replicate",
      "bootstrap (Table 2); the bootstrap median for Vd_pop, 79.3 L, sits about 6% above the 75.0 L",
      "final estimate, and every other bootstrap median is within a few percent of its point estimate.",
      collapse = " "
    )
  )

  ini({
    # Structural parameters -- Polaskova 2024 Table 2 "Fixed effects".
    #
    # NOTE ON PARAMETERIZATION: both covariate effects on CL are exponential in
    # the RAW, UNCENTERED covariate (Results 3.2 equation block, and the
    # Discussion's explicit form "CL = 1.32 x e^(0.61 x eGFR) x e^(0.011 x
    # LBM)"). lcl is therefore the clearance extrapolated to eGFR = 0 AND
    # LBM = 0, not a clearance at any covariate reference value, and it is not
    # physically interpretable on its own. The Discussion's worked example
    # reproduces exactly: at LBM 68 kg and eGFR 1.51 mL/s/1.73 m^2 (the cohort
    # medians), 1.32 * exp(0.61 * 1.51) * exp(0.011 * 68) = 7.01 L/h, which the
    # paper rounds to "7 L/h", and 0.693 * 75 / 7.01 = 7.42 h against the
    # paper's stated t1/2 of 7.4 h.
    lvc <- log(75.0); label("Volume of distribution (L)")                            # Table 2, Vd_pop = 75.0 L (R.S.E. 8.66%; bootstrap median 79.3, 95% CI 77.8-80.9)
    lcl <- log(1.32); label("Clearance at eGFR = 0 and LBM = 0 (L/h)")                # Table 2, CL_pop = 1.32 L/h (R.S.E. 19.3%; bootstrap median 1.27, 95% CI 1.24-1.31)

    # Covariate effects on CL -- exponential, uncentered.
    e_crcl_cl <- 0.61; label("Exponential eGFR effect on CL (per mL/s/1.73 m^2)")     # Table 2, beta_CL_eGFR = 0.61 (R.S.E. 11.8%; bootstrap median 0.61, 95% CI 0.60-0.63)
    e_lbm_cl <- 0.011; label("Exponential lean-body-mass effect on CL (per kg)")      # Table 2, beta_CL_LBM = 0.011 (R.S.E. 21.6%; bootstrap median 0.011, 95% CI 0.0108-0.0113)

    # IIV. Table 2 reports these under "Standard deviation of the random
    # effects" and the table footnote defines Omega as "standard deviation of
    # the random effects", so the published 0.31 / 0.28 are SDs on the log
    # scale and are squared here to give the variances nlmixr2 expects.
    # Monolix's log-normal parameter distribution matches the exp(l* + eta*)
    # form used in model().
    etalvc ~ 0.0961; label("IIV on volume of distribution (variance)")                # Table 2, Omega_Vd = 0.31 SD -> 0.31^2 = 0.0961 (R.S.E. 16.6%; CV 31.7%)
    etalcl ~ 0.0784; label("IIV on clearance (variance)")                             # Table 2, Omega_CL = 0.28 SD -> 0.28^2 = 0.0784 (R.S.E. 9.92%; CV 28.5%)

    # Residual error. Results 3.2: "A constant error model was the most
    # accurate for the description of residual [error]" -- i.e. purely
    # additive, in the assay's mg/L.
    addSd <- 2.9; label("Additive residual error (mg/L)")                             # Table 2, "Constant" = 2.9 mg/L (R.S.E. 20.7%; bootstrap median 2.2, 95% CI 2.1-2.3)
  })

  model({
    # CRCL is supplied in the canonical mL/min/1.73 m^2; the published
    # coefficient 0.61 is per mL/s/1.73 m^2, so convert first.
    eGFRsi <- CRCL / 60

    # Polaskova 2024 Results 3.2:
    #   Vd = Vd_pop
    #   CL = CL_pop * exp(beta_CL_eGFR * eGFR) * exp(beta_CL_LBM * LBM)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * exp(e_crcl_cl * eGFRsi) * exp(e_lbm_cl * LBM)

    kel <- cl / vc

    # One compartment, linear elimination, IV infusion only (Results 3.2:
    # "A one-compartment model with linear elimination kinetics best-fitted
    # vancomycin concentration-time data"). Vancomycin is not absorbed orally,
    # so there is no depot; dose all events into `central` with a rate or
    # duration.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
