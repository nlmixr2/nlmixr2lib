Renaud_2020_glenzocimab <- function() {
  description <- "Two-compartment population PK model with a direct (immediate) Imax model of ex vivo collagen-induced platelet aggregation for glenzocimab (ACT017), an anti-GPVI Fab, in healthy volunteers (Renaud 2020)"
  reference <- "Renaud L, Lebozec K, Voors-Pette C, Dogterom P, Billiald P, Jandrot Perrus M, Pletan Y, Machacek M. Population Pharmacokinetic/Pharmacodynamic Modeling of Glenzocimab (ACT017) a Glycoprotein VI Inhibitor of Collagen-Induced Platelet Aggregation. J Clin Pharmacol. 2020;60(9):1198-1208. doi:10.1002/jcph.1616"
  vignette <- "Renaud_2020_glenzocimab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "glenzocimab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "glenzocimab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight; power effect on CL, V1 and Q normalized to 70 kg (Renaud 2020 Table 2 footnotes d, f, h).",
      source_name = "BW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline age; power effect on CL and Q normalized to 50 years (Renaud 2020 Table 2 footnotes c, g).",
      source_name = "AGE"
    ),
    CREAT = list(
      description = "Serum creatinine concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline plasma creatinine, used by the authors as a proxy for eGFR; power effect on CL normalized to 0.79 mg/dL (Renaud 2020 Table 2 footnote e).",
      source_name = "CRE"
    ),
    PLT = list(
      description = "Baseline platelet count",
      units = "10^9/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline platelet count; power effect on Imax normalized to 220 x 10^9/L (Renaud 2020 Table 2 footnote j).",
      source_name = "PLT"
    ),
    DOSE = list(
      description = "Total glenzocimab dose administered over the infusion",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (a): per-subject assigned total dose (62.5-2000 mg in the phase I study). Power effect on IC50 normalized to 500 mg (Renaud 2020 Table 2 footnote i); the authors could not identify a mechanism for the dose-dependent potency.",
      source_name = "DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 36,
    n_studies = 1,
    age_range = "22-63 years",
    age_median = "56 years",
    weight_range = "52-107 kg",
    weight_median = "74 kg",
    sex_female_pct = 36.1,
    race_ethnicity = c(White = 91.7, Other = 8.3),
    disease_state = "Healthy volunteers (81% normal renal function, 19% mild renal impairment by eGFR)",
    dose_range = "62.5, 125, 250, 500, 1000 or 2000 mg single IV infusion over 6 h (25% of the dose in the first 15 min, 75% over the remaining 5 h 45 min)",
    regions = "Netherlands",
    creatinine_range = "0.46-1.2 mg/dL (median 0.77)",
    platelet_range = "159-323 x 10^9/L (median 207)",
    notes = "Phase I first-in-human single-ascending-dose study (6 active + 2 placebo per dose group); 390 PK and 404 PD observations from the 36 actively treated subjects (Renaud 2020 Table 1 and Methods)."
  )

  ini({
    lcl <- log(2.67); label("Clearance for a 70 kg, 50-year-old subject with creatinine 0.79 mg/dL (L/h)")  # Table 2 'Cl (L/h)' = 2.67
    lvc <- log(4.1); label("Central volume of distribution V1 for a 70 kg subject (L)")  # Table 2 'V1 (L)' = 4.1
    lq <- log(0.626); label("Intercompartmental clearance for a 70 kg, 50-year-old subject (L/h)")  # Table 2 'Q (L/h)' = 0.626
    lvp <- log(6.89); label("Peripheral volume of distribution V2 (L)")  # Table 2 'V2 (L)' = 6.89

    lrbase <- log(79.8); label("Baseline ex vivo platelet aggregation Base_PPA (%)")  # Table 2 'Base_PPA (%)' = 79.8
    lic50 <- log(0.924); label("Glenzocimab concentration for half-maximal inhibition at a 500 mg dose (ug/mL)")  # Table 2 'IC50 (ug/mL)' = 0.924
    limax <- log(72.9); label("Maximum reduction in platelet aggregation at a platelet count of 220 x 10^9/L (percentage points)")  # Table 2 'Imax (%)' = 72.9

    e_age_cl <- -0.304; label("Power exponent of age (/50 years) on CL (unitless)")  # Table 2 'beta Cl_tAGE' = -0.304
    e_wt_cl <- 1.09; label("Power exponent of body weight (/70 kg) on CL (unitless)")  # Table 2 'beta Cl_tBW' = 1.09
    e_creat_cl <- -0.566; label("Power exponent of creatinine (/0.79 mg/dL) on CL (unitless)")  # Table 2 'beta Cl_tCRE' = -0.566
    e_wt_vc <- 0.694; label("Power exponent of body weight (/70 kg) on V1 (unitless)")  # Table 2 'beta V1_tBW' = 0.694
    e_age_q <- -0.318; label("Power exponent of age (/50 years) on Q (unitless)")  # Table 2 'beta Q_tAGE' = -0.318
    e_wt_q <- 0.812; label("Power exponent of body weight (/70 kg) on Q (unitless)")  # Table 2 'beta Q_tBW' = 0.812
    e_dose_ic50 <- -0.989; label("Power exponent of total dose (/500 mg) on IC50 (unitless)")  # Table 2 'beta IC50_tDOSE' = -0.989
    e_plt_imax <- 0.17; label("Power exponent of platelet count (/220 x 10^9/L) on Imax (unitless)")  # Table 2 'beta Imax_tPLT' = 0.17

    # Monolix log-normal SDs (Table 2 'Standard deviations') squared to variances;
    # covariances = r * omega_i * omega_j from Table 2 'Correlations'.
    # omega Cl = 0.182, omega V1 = 0.148, omega Q = 0.144;
    # r V1_Cl = 0.626, r Q_Cl = 0.796, r V1_Q = 0.84
    etalcl + etalvc + etalq ~ c(
      0.033124,
      0.016862, 0.021904,
      0.020862, 0.017902, 0.020736
    )
    etalvp ~ 0.070225 # Table 2 'omega V2' = 0.265 (SD), squared
    etalic50 ~ 1.8225 # Table 2 'omega IC50' = 1.35 (SD), squared

    addSd <- 0.0869; label("Additive residual error for glenzocimab concentration (ug/mL)")  # Table 2 'a1 (ug/mL)' = 0.0869
    propSd <- 0.0514; label("Proportional residual error for glenzocimab concentration (fraction)")  # Table 2 'b1 (-)' = 0.0514
    addSd_PPA <- 0.778; label("Additive residual error of platelet aggregation on the logit scale (unitless)")  # Table 2 'a2 (-)' = 0.778
  })
  model({
    # Individual PK parameters (Renaud 2020 Eq. 3 power covariate model with
    # log-normal IIV; references 70 kg, 50 years, 0.79 mg/dL)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (AGE / 50)^e_age_cl * (CREAT / 0.79)^e_creat_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_q * (AGE / 50)^e_age_q
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc

    # Immediate-response Imax model (Renaud 2020 Eq. 1), dose covariate on
    # IC50 (reference 500 mg) and platelet-count covariate on Imax (reference
    # 220 x 10^9/L); Base_PPA has no IIV.
    rbase <- exp(lrbase)
    ic50 <- exp(lic50 + etalic50) * (DOSE / 500)^e_dose_ic50
    imax <- exp(limax) * (PLT / 220)^e_plt_imax
    PPA <- rbase - imax * Cc / (ic50 + Cc)

    # Combined error (Renaud 2020 Eq. 5, Monolix combined1: sd = a1 + b1 * f)
    Cc ~ add(addSd) + prop(propSd) + combined1()
    # Logit-normal constant error on the 0-100 % scale (Renaud 2020 Eq. 6)
    PPA ~ logitNorm(addSd_PPA, 0, 100)
  })
}
