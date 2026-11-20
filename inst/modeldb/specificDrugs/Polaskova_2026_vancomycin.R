Polaskova_2026_vancomycin <- function() {
  description <- "Two-compartment (serum + cerebrospinal fluid) population PK model for intravenous vancomycin in 64 adult neurosurgical patients with secondary post-neurosurgical meningitis (Polaskova 2026). Drug distributes from the central (serum) compartment, Vd 115.1 L, into a CSF compartment whose volume is fixed at 0.15 L, through a bidirectional inter-compartmental clearance Q of 0.0079 L/h; it leaves the body by systemic clearance from serum (6.50 L/h at an absolute CKD-EPI eGFR of 120 mL/min, power exponent 1.20) and by a separate clearance from the CSF compartment (0.019 L/h). Serum and CSF each carry a proportional residual error. The typical steady-state CSF-to-serum AUC ratio is Q / (Q + CLCSF) = 0.29."
  reference <- "Polaskova L, Moravec T, Krejci V, Netuka D, Slanar O, Sima M. Population pharmacokinetics of vancomycin in serum and cerebrospinal fluid: variability in cerebrospinal fluid concentrations despite targeted serum exposure. Antimicrob Agents Chemother. 2026;70(8):e00569-26. doi:10.1128/aac.00569-26"
  vignette <- "Polaskova_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE),
    csf = list(analyte = "vancomycin", units = "mg", specimen = "CSF", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate by the creatinine CKD-EPI",
        "equation, in ABSOLUTE mL/min (NOT indexed to 1.73 m^2)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source name eGFR. Methods 'Data retrieval': the creatinine-based GFR",
        "was estimated with CKD-EPI 'both as a normalized value (eGFRnorm:",
        "mL/min/1.73 m2) and as an absolute value (eGFR: mL/min)'; both were",
        "screened and Results 'Population pharmacokinetic analysis' names the",
        "ABSOLUTE eGFR as 'the most influential covariate on vancomycin CL,",
        "described using a power function' (dOFV -23.5). Table 2 equation:",
        "CL = CL_pop x (eGFR/120)^beta_CL_eGFR x exp(eta_CL). Table 1:",
        "median 119.8 mL/min, IQR 101.0-136.0, range 56.2-185.5. The paper",
        "does not say which BSA formula de-indexed the CKD-EPI value; Methods",
        "name Du Bois as the BSA formula used in the study. Supplying the",
        "BSA-normalized value (Table 1 median 106.5 mL/min/1.73 m^2) in this",
        "column under-predicts clearance in a typical adult with BSA > 1.73",
        "m^2 (cohort median BSA 2.03 m^2)."
      ),
      source_name = "eGFR"
    )
  )

  # Covariates the paper screened (Methods 'Population pharmacokinetic
  # analysis') but did not retain: Results state 'No other tested demographic
  # or laboratory variables showed a statistically significant effect on the
  # pharmacokinetic parameters.' No point estimates were published.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age. Screened; not retained.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 55 years, IQR 45-67.5, range 18-82."
    ),
    HT = list(
      description = "Height. Screened; not retained.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 175 cm, IQR 167-180, range 154-197."
    ),
    WT = list(
      description = "Total body weight. Screened; not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 85 kg, IQR 70-95, range 48-140."
    ),
    IBW = list(
      description = "Ideal body weight. Screened; not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 66 kg, IQR 58-75, range 47-87."
    ),
    ABW = list(
      description = "Adjusted body weight. Screened; not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 75 kg, IQR 65-82, range 48-100."
    ),
    LBM = list(
      description = "Lean body mass (Boer formula). Screened; not retained.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 57 kg, IQR 49-65, range 35-80."
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula). Screened; not retained.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1: median 2.03 m^2, IQR 1.85-2.18, range 1.51-2.62."
    ),
    SEXF = list(
      description = "Female sex indicator. Screened as a categorical covariate; not retained.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 1: 36 male / 28 female (43.8% female; the table prints 44.8%, but 28/64 = 43.75%)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 1L,
    n_observations = 233L,
    age_range = "18-82 years",
    age_median = "55 years",
    weight_range = "48-140 kg",
    weight_median = "85 kg",
    sex_female_pct = 43.8,
    race_ethnicity = NULL,
    disease_state = "Adult neurosurgical inpatients treated with vancomycin for secondary (post-neurosurgical) meningitis; no encephalitis or isolated ventriculitis. CSF sampled via external ventricular or lumbar drains or clinically indicated lumbar puncture.",
    dose_range = "Loading dose 1,500-4,000 mg in 53 patients; intermittent doses 500-3,000 mg every 6, 8, 12 or 24 h (53 patients intermittent only), 4 patients continuous infusion only, 7 both.",
    regions = "Czechia (single centre: Military University Hospital Prague)",
    renal_function = "Absolute CKD-EPI eGFR median 119.8 mL/min (IQR 101.0-136.0, range 56.2-185.5); BSA-normalized median 106.5 mL/min/1.73 m^2 (range 53.4-138). No patient was on haemodialysis or ECMO.",
    notes = paste(
      "Retrospective observational study of routine therapeutic drug monitoring,",
      "January 2017 to January 2026. Only concentrations from the first 4 days of",
      "therapy were analysed: 133 serum concentrations (1-4 per patient) and 100",
      "CSF observations (62 quantifiable, 38 below the 4.0 mg/L LOQ, handled as",
      "censored data in Monolix, an M3-like method). Assay: KIMS",
      "immunoturbidimetry on a Roche Cobas 8000 for both matrices. Estimation by",
      "SAEM in Monolix 2024R1; 500-replicate bootstrap (Table 2). Baseline",
      "demographics: Table 1.",
      collapse = " "
    )
  )

  ini({
    # Structural parameters -- Polaskova 2026 Table 2 'Fixed effects'. The
    # structure follows the Monolix model file printed in the Supplementary
    # data: cmt 1 = serum (volume Vd), cmt 2 = CSF (volume VdCSF = 0.15),
    # k12 = Q/Vd, k21 = Q/VdCSF, k10 = CL/Vd, k20 = CLCSF/VdCSF.
    lvc <- log(115.1); label("Volume of distribution of the central (serum) compartment (L)") # Table 2: Vd_pop = 115.1 L (R.S.E. 10.7%)
    lcl <- log(6.50); label("Systemic clearance from serum at eGFR = 120 mL/min (L/h)") # Table 2: CL_pop = 6.50 L/h (R.S.E. 4.61%)
    e_crcl_cl <- 1.20; label("Power exponent of absolute eGFR (mL/min, reference 120) on CL (unitless)") # Table 2: beta_CL_eGFR = 1.20 (R.S.E. 17.5%)
    lqcsf <- log(0.0079); label("Inter-compartmental clearance between serum and CSF (L/h)") # Table 2: Q_pop = 0.0079 L/h (R.S.E. 13.1%)
    lcl_csf <- log(0.019); label("Clearance out of the body from the CSF compartment (L/h)") # Table 2: CLCSF_pop = 0.019 L/h (R.S.E. 13.8%)
    lvcsf <- fixed(log(0.15)); label("Volume of the CSF compartment (L)") # Table 2 and Methods: VdCSF fixed to 0.15 L

    # IIV. Monolix parameterises log-normal random effects by their standard
    # deviation omega; Methods define shrinkage as 1 - Var(eta_i)/omega^2, so
    # the Table 2 'Between-subject variability (%)' values are 100 x omega.
    # They are squared here to the variances nlmixr2 expects.
    etalvc ~ 0.1681 # Table 2: omega_Vd = 41.0% -> 0.41^2 (R.S.E. 19.5%)
    etalcl ~ 0.09 # Table 2: omega_CL = 30.0% -> 0.30^2 (R.S.E. 11.4%)
    etalqcsf ~ 0.0961 # Table 2: omega_Q = 31.0% -> 0.31^2 (R.S.E. 43.9%)
    etalcl_csf ~ 0.1764 # Table 2: omega_CLCSF = 42.0% -> 0.42^2 (R.S.E. 22.4%)

    # Residual error: proportional in both matrices (Results 'Population
    # pharmacokinetic analysis').
    propSd <- 0.15; label("Proportional residual error, serum (fraction)") # Table 2: Proportional (serum) = 0.15 (R.S.E. 10.7%)
    propSd_Ccsf <- 0.16; label("Proportional residual error, CSF (fraction)") # Table 2: Proportional (CSF) = 0.16 (R.S.E. 18.4%)
  })

  model({
    # Table 2 covariate and parameter equations
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (CRCL / 120)^e_crcl_cl
    qcsf <- exp(lqcsf + etalqcsf)
    cl_csf <- exp(lcl_csf + etalcl_csf)
    vcsf <- exp(lvcsf)

    # Supplementary Monolix model file rate constants
    kel <- cl / vc
    k_central_csf <- qcsf / vc
    k_csf_central <- qcsf / vcsf
    kel_csf <- cl_csf / vcsf

    # Intravenous doses go to central (Monolix iv(cmt=1)); elimination from
    # both compartments to the outside (elimination(cmt=1), elimination(cmt=2)).
    d/dt(central) <- -kel * central - k_central_csf * central + k_csf_central * csf
    d/dt(csf) <- k_central_csf * central - k_csf_central * csf - kel_csf * csf

    Cc <- central / vc
    Ccsf <- csf / vcsf

    Cc ~ prop(propSd)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
