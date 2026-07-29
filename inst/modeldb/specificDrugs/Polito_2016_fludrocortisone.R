Polito_2016_fludrocortisone <- function() {
  description <- paste0(
    "One-compartment population PK model for oral fludrocortisone with ",
    "first-order absorption, an absorption lag time, and first-order ",
    "elimination, estimated in 14 adults with septic shock (out of 21 ",
    "enrolled; 7 had undetectable plasma concentrations) receiving a ",
    "single 50 ug oral dose of fludrocortisone acetate via naso-gastric ",
    "tube (Polito 2016). The Simplified Acute Physiology Score II ",
    "(SAPS II) is retained as a positive power covariate on both ",
    "apparent oral clearance CL/F (exponent 0.019) and absorption lag ",
    "time Tlag (exponent 0.036), normalised to the cohort median ",
    "SAPS II = 53. Inter-individual variability is exponential on every ",
    "PK parameter (ka, V/F, CL/F, Tlag) with a diagonal OMEGA matrix; ",
    "residual error is proportional."
  )
  reference <- paste0(
    "Polito A, Hamitouche N, Ribot M, Polito A, Laviolle B, Bellissant E, ",
    "Annane D, Alvarez JC. Pharmacokinetics of oral fludrocortisone in ",
    "septic shock. Br J Clin Pharmacol. 2016;82(6):1509-1516. ",
    "doi:10.1111/bcp.13065."
  )
  vignette <- "Polito_2016_fludrocortisone"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ug/L"
  )

  covariateData <- list(
    SAPS_II = list(
      description        = "New Simplified Acute Physiology Score II at ICU admission",
      units              = "points",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed per subject (computed during the first 24 hours of ",
        "ICU admission per Le Gall, Lemeshow & Saulnier 1993). Power ",
        "effect on CL/F and Tlag: theta_i = theta_pop * exp(eta_i) * ",
        "(SAPS_II / 53)^beta with reference 53 = the cohort median ",
        "SAPS II among the 14 patients with detectable plasma ",
        "fludrocortisone concentrations (Table 1, 'Yes (n = 14)' ",
        "column). Exponents are positive (Table 2): higher SAPS II is ",
        "associated with a longer absorption lag and a faster apparent ",
        "oral clearance. RSE 39 % on the Tlag exponent and 35 % on the ",
        "CL/F exponent."
      ),
      source_name        = "SAPSII"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 14L,
    n_enrolled      = 21L,
    n_studies       = 1L,
    age_range       = "Adults (detectable subset: median 65 years, IQR 57-75)",
    age_median      = "65 years (detectable subset)",
    weight_range    = "Detectable subset: median 71 kg, IQR 60-84",
    weight_median   = "71 kg (detectable subset)",
    sex_female_pct  = 43,
    saps_ii_median  = 53,
    disease_state   = paste0(
      "Adults with septic shock (clinically or microbiologically ",
      "documented source of infection; at least two signs of tissue ",
      "hypoperfusion / organ dysfunction such as urinary output ",
      "< 20 mL/h, Glasgow coma score < 10, mechanical ventilation, or ",
      "arterial lactate > 2 mmol/L; need for vasopressor therapy to ",
      "maintain SBP > 90 mmHg or MAP > 65 mmHg). Patients with known ",
      "endocrine disorders or any condition / treatment that may affect ",
      "cortisol synthesis or metabolism were excluded."
    ),
    dose_range      = paste0(
      "Single 50 ug oral dose of fludrocortisone acetate administered ",
      "via naso-gastric tube within the first 3 hours of septic shock ",
      "onset and prior to any other corticotherapy (Methods 'Patients ",
      "and settings')."
    ),
    regions         = "France (single-centre, Raymond Poincare Hospital, Garches).",
    sampling_design = paste0(
      "Arterial blood sampled pre-dose and every 30 min for 6 h, then ",
      "hourly to 18 h. Plasma fludrocortisone quantified by LC-MS/MS ",
      "(Oasis-HLB solid-phase extraction); LLOQ 0.10 ug/L. BLQ ",
      "observations handled in the likelihood (censored-data SAEM)."
    ),
    notes           = paste0(
      "Ancillary study to the CRISTAL trial (NCT00318942) at Raymond ",
      "Poincare Hospital (Garches, France), recruiting December 2010 to ",
      "May 2012. Of 21 enrolled patients, 7 (33 %) had undetectable ",
      "plasma fludrocortisone (< 0.10 ug/L) and were excluded from the ",
      "structural PK model; the remaining 14 patients informed all ",
      "parameter estimates in Table 2. Hospital mortality 50 % in the ",
      "detectable subset. Modelling software: MONOLIX 4.3.0 (LIXOFT, ",
      "France) with the Stochastic Approximation Expectation Maximization ",
      "(SAEM) algorithm. Diagonal OMEGA matrix (no IIV covariances ",
      "retained). Covariate selection: forward inclusion at delta-LL > ",
      "3.84 (p < 0.05); backward elimination at delta-LL > 6.64 ",
      "(p < 0.01). Of nine screened covariates (age, sex, weight, ",
      "SAPS II, total protein, albumin, creatinine, diuresis, proton ",
      "pump inhibitor use), only SAPS II was retained, on CL/F and on ",
      "Tlag."
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Adult age (years)",
      units       = "years",
      type        = "continuous",
      notes       = paste0(
        "Screened in the covariate analysis but not retained in the ",
        "Polito 2016 final model (Methods 'Population pharmacokinetic ",
        "analysis'; Results 'Population pharmacokinetic analysis')."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened as 'gender' in the covariate analysis but not retained."
      )
    ),
    WT = list(
      description = "Total body weight (kg)",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Screened as 'total body weight' in the covariate analysis but ",
        "not retained."
      )
    ),
    PROTEIN = list(
      description = "Serum total protein (g/L)",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    ALBUMIN = list(
      description = "Serum albumin (g/L)",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    CREATININE = list(
      description = "Serum creatinine (umol/L)",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    DIURESIS = list(
      description = "24-hour urine output (mL/24 h)",
      units       = "mL/24h",
      type        = "continuous",
      notes       = "Screened but not retained."
    ),
    PPI = list(
      description = "Proton pump inhibitor co-medication indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened but not retained as a fitted covariate in the final ",
        "model. Discussion suggests PPI co-medication may reduce ",
        "fludrocortisone absorption (9/14 in the detectable group vs ",
        "7/7 in the undetectable group, p = 0.06)."
      )
    )
  )

  ini({
    # ===== Structural parameters (Polito 2016 Table 2, 'Population parameter') =====
    # Reference subject: SAPS II = 53 (cohort median for the 14 detectable
    # patients; Table 1 'Yes (n = 14)' column). The paper reports CL/F in
    # L/h, V/F in L, ka in 1/h and Tlag in h. RSEs reported in Table 2.
    lka   <- log(0.67); label("First-order absorption rate constant ka (1/h)")              # Table 2 (RSE 23 %)
    lvc   <- log(78);   label("Apparent volume of distribution V/F (L)")                    # Table 2 (RSE 28 %)
    lcl   <- log(40);   label("Apparent oral clearance CL/F at SAPS II = 53 (L/h)")         # Table 2 (RSE 15 %)
    ltlag <- log(0.65); label("Absorption lag time Tlag at SAPS II = 53 (h)")               # Table 2 (RSE 34 %)

    # ===== Covariate effects on CL/F and Tlag (Polito 2016 Table 2 footnote b) =====
    # Multiplicative power-of-(SAPS_II / median) form per Methods
    # 'Population pharmacokinetic analysis' (continuous-covariate
    # equation) and Table 2 footnote b:
    #   theta_i = theta_pop * exp(eta_i) * (SAPS_II_i / 53)^beta
    # The exponents are positive: higher SAPS II is associated with
    # longer absorption lag and faster apparent oral clearance.
    e_saps_ii_cl   <- 0.019; label("Power exponent of (SAPS_II / 53) on CL/F (unitless)")   # Table 2 (RSE 35 %, p = 0.004)
    e_saps_ii_tlag <- 0.036; label("Power exponent of (SAPS_II / 53) on Tlag (unitless)")   # Table 2 (RSE 39 %, p = 0.037)

    # ===== Inter-individual variability (Polito 2016 Table 2 'Inter-individual variability') =====
    # Exponential IIV per Methods 'Population pharmacokinetic analysis':
    #   theta_i = theta_pop * exp(eta_i)  with  eta_i ~ N(0, omega_theta)
    # Table 2 reports omega in the column 'omega (%)' with values 42, 75,
    # 49 and 98; following the MONOLIX convention (and the Abboud 2009
    # BJCP septic-shock precedent in nlmixr2lib) these are interpreted as
    # the standard deviation of eta on the log scale expressed as a
    # percentage, i.e., omega_ka = 0.42, omega_V/F = 0.75, omega_CL/F =
    # 0.49 and omega_Tlag = 0.98. Variance = omega^2. Diagonal OMEGA
    # matrix (no covariances retained).
    etalka   ~ 0.1764  # Table 2: omega_ka   = 0.42; variance = 0.42^2 = 0.1764 (RSE 36 %)
    etalvc   ~ 0.5625  # Table 2: omega_V/F  = 0.75; variance = 0.75^2 = 0.5625 (RSE 23 %)
    etalcl   ~ 0.2401  # Table 2: omega_CL/F = 0.49; variance = 0.49^2 = 0.2401 (RSE 23 %)
    etaltlag ~ 0.9604  # Table 2: omega_Tlag = 0.98; variance = 0.98^2 = 0.9604 (RSE 24 %)

    # ===== Residual error (Polito 2016 Table 2 'Residual variability') =====
    # Proportional residual model per Methods and Table 2 footnote:
    #   C_obs_ij = C_pred_ij + C_pred_ij * eps_prop_ij  with
    #   eps_prop_ij ~ N(0, sigma_prop)
    # sigma_prop = 0.20 is the SD on the proportional scale (RSE 25 %).
    propSd <- 0.20; label("Proportional residual error (fraction)")                          # Table 2 (RSE 25 %)
  })

  model({
    # Individual PK parameters. Reference subject: SAPS_II = 53.
    # Final CL/F and Tlag relationships per Methods Eq. for continuous
    # covariates and Table 2 footnote b:
    #   CL/F_i = 40   * exp(eta_lcl_i)   * (SAPS_II_i / 53)^0.019  (L/h)
    #   Tlag_i = 0.65 * exp(eta_ltlag_i) * (SAPS_II_i / 53)^0.036  (h)
    # ka and V/F have no covariate effect retained in the final model.
    ka   <- exp(lka   + etalka)
    vc   <- exp(lvc   + etalvc)
    cl   <- exp(lcl   + etalcl)   * (SAPS_II / 53)^e_saps_ii_cl
    tlag <- exp(ltlag + etaltlag) * (SAPS_II / 53)^e_saps_ii_tlag

    kel <- cl / vc

    # One-compartment open model with first-order absorption + lag.
    # Bioavailability F is not separately identifiable from CL and V in
    # this oral-only design; both CL and V are reported as apparent
    # (CL/F and V/F). No explicit f(depot) term is needed (F implicitly
    # = 1 in the structural model).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag time on the depot.
    alag(depot) <- tlag

    # Observation and proportional residual error. Cc has units ug/L
    # (depot amt in ug, vc in L) matching the LC-MS/MS assay scale.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
