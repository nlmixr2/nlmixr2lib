Chen_2015_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order elimination for intravenous voriconazole in Chinese adult critically ill patients with pulmonary disease (Chen 2015); direct bilirubin enters as a power-form covariate on clearance."
  reference <- "Chen W, Xie H, Liang F, Meng D, Rui J, Yin X, Zhang T, Xiao X, Cai S, Liu X, Li Y. Population pharmacokinetics in China: the dynamics of intravenous voriconazole in critically ill patients with pulmonary disease. Biol Pharm Bull. 2015;38(7):996-1004. doi:10.1248/bpb.b14-00768"
  vignette <- "Chen_2015_voriconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DBIL = list(
      description        = "Direct (conjugated) serum bilirubin concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL: CL_typical = TVCL * (DBIL / 2.6)^e_dbil_cl with reference 2.6 umol/L (Chen 2015 Table 2 typical-value equation and Table 1 cohort median direct bilirubin 3.16 +/- 2.21 umol/L; the reference 2.6 is the value used inside the model equation as printed in Table 2 footnote rather than the Table 1 mean). The estimated exponent is -0.40 (95% CI -0.69 to -0.11), so a doubling of DBIL above the reference scales CL by 2^-0.40 = 0.76 (~24% decrease).",
      source_name        = "DBIL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 62L,
    n_studies      = 1L,
    n_observations = 240L,
    age_range      = "19-90 years",
    age_mean       = "59.71 +/- 16.67 years",
    weight_range   = "41-84 kg",
    weight_mean    = "60.13 +/- 10.03 kg",
    sex_female_pct = 32.3,
    race_ethnicity = c(Chinese = 100),
    disease_state  = "Adult intensive-care-unit patients with pulmonary diseases (clinical indication for invasive fungal infection; 84% with positive microbiology, predominantly aspergillosis 60% and candidiasis 19%). Hepatic and renal function ranged across normal and mildly impaired; severe renal impairment (CLCR < 50 mL/min) was an exclusion criterion under the hospital's voriconazole-stopping policy.",
    dose_range     = "Loading dose 300 mg intravenous voriconazole followed by 200 mg intravenous infusion every 12 h; therapeutic drug monitoring started at trough before the sixth maintenance dose (steady state after 72 h).",
    regions        = "Single center: The First Affiliated Hospital of Guangzhou Medical University, Guangzhou, P.R. China.",
    notes          = "Prospective observational study, March 2012 - May 2013. 240 plasma concentration samples collected at 0.5, 1, 1.5, 2, 4, 6, 9, and 12 h after the start of an infusion (or any two of those points per occasion). Bioanalytical: HPLC-UV with LLOQ 70 ng/mL, calibration range 208-20800 ng/mL. Baseline demographics per Chen 2015 Table 1; final population PK parameter estimates per Chen 2015 Table 2. APACHE II 21.6 +/- 13.8, SOFA median 4, observed mortality 14.5%. Hepatotoxicity in 24.2% of patients."
  )

  ini({
    # Structural parameters -- typical-value reference is a patient with
    # direct bilirubin 2.6 umol/L (Chen 2015 Table 2 typical-value
    # equation CL = theta_CL * (DBIL/2.6)^theta_DBIL). Voriconazole is
    # administered intravenously only, so absorption and bioavailability
    # are not part of the model.
    lcl <- log(4.28); label("Clearance at the DBIL = 2.6 umol/L reference (L/h)")  # Chen 2015 Table 2 final model: theta_CL = 4.28 L/h (95% CI 3.48-5.08)
    lvc <- log(93.4); label("Central volume of distribution (L)")                  # Chen 2015 Table 2 final model: theta_Vd = 93.4 L (95% CI 79.39-107.41)

    # Covariate effect on CL. Power-form (Chen 2015 Table 2 typical-value
    # equation): CL = theta_CL * (DBIL / 2.6)^theta_DBIL with theta_DBIL
    # = -0.40 (95% CI -0.69 to -0.11). The minus sign means higher DBIL
    # decreases CL, consistent with the paper's interpretation that
    # patients with poor liver function have impaired voriconazole
    # clearance (Discussion section).
    e_dbil_cl <- -0.40; label("Power exponent for direct bilirubin on CL (unitless)")  # Chen 2015 Table 2 final model: theta_DBIL = -0.40

    # IIV. Chen 2015 reports exponential IIV in NONMEM convention
    # (CL_i = TVCL * exp(eta_CL_i)) with reported approximate-CV%
    # values 72.94% on CL and 26.50% on Vc (final model, Table 2). Using
    # the standard NONMEM convention CV% ~= sqrt(omega^2) * 100, the
    # internal variance is omega^2 = (CV/100)^2: 0.7294^2 = 0.5320 for
    # CL and 0.2650^2 = 0.0702 for Vc. The exact log-normal variance
    # log(1 + CV^2) is 0.4234 for CL and 0.0680 for Vc; the approximate
    # form is used here to align with the bootstrap symmetry shown in
    # Table 2 (CV-as-omega: 72.67% median, 58.30-85.92% 95% CI for CL).
    etalcl ~ 0.5320                                                             # Chen 2015 Table 2 final model: omega_CL = 72.94% CV; see comment above
    etalvc ~ 0.0702                                                             # Chen 2015 Table 2 final model: omega_Vd = 26.50% CV

    # Residual error. Chen 2015 describes the residual model as a
    # "constant coefficient" model with delta = 13.0%, i.e., a
    # proportional residual error on the linear-concentration scale.
    propSd <- 0.130; label("Proportional residual error (fraction)")             # Chen 2015 Table 2 final model: delta = 13.0%
  })

  model({
    # Reference direct-bilirubin value taken from Chen 2015 Table 2
    # typical-value equation (the value 2.6 appearing in the printed
    # formula CL = theta_CL * (DBIL/2.6)^theta_DBIL).
    ref_dbil <- 2.6

    # Individual parameters. Power-form covariate on CL.
    cl <- exp(lcl + etalcl) * (DBIL / ref_dbil)^e_dbil_cl
    vc <- exp(lvc + etalvc)

    # One-compartment IV PK with first-order elimination -- linCmt()
    # resolves to the analytical solution given cl and vc (NONMEM
    # ADVAN1 TRANS2 equivalent).
    Cc <- linCmt()

    # Concentration units: dose mg / volume L = mg/L = ug/mL, matching
    # the source paper's plasma-concentration units (ug/mL throughout
    # Results and Tables).
    Cc ~ prop(propSd)
  })
}
