Moes_2022_tocilizumab <- function() {
  description <- "One-compartment population PK model for intravenous tocilizumab in ICU-admitted adults with COVID-19 co-treated with dexamethasone (Moes 2022), with parallel first-order linear and Michaelis-Menten elimination from the central compartment and correlated IIV on CL and V; no covariates retained."
  reference <- "Moes DJAR, van Westerloo DJ, Arend SM, Swen JJ, de Vries A, Guchelaar HJ, Joosten SA, de Boer MGJ, van Gelder T, van Paassen J. Towards Fixed Dosing of Tocilizumab in ICU-Admitted COVID-19 Patients: Results of an Observational Population Pharmacokinetic and Descriptive Pharmacodynamic Study. Clin Pharmacokinet. 2022;61(2):231-247. doi:10.1007/s40262-021-01074-2"
  vignette <- "Moes_2022_tocilizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "tocilizumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight (baseline)",
      units = "kg",
      type = "continuous",
      notes = "Screened, not retained. Tested as a power (allometric) effect on CL (Section 3.5, Eq. 1); the exponent was estimated at 0.002, did not significantly improve the fit and did not reduce the IIV on CL. Cohort mean 96.4 kg (range 58-130; Table 1)."
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Screened, not retained (Section 3.3, all p > 0.05). Cohort mean 64 years (range 45-80; Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, not retained (Section 3.3). 72.4% of the cohort was male (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (baseline and time-varying), not retained (Section 3.3). Cohort mean 32 g/L (range 23-37; Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Screened (baseline and time-varying), not retained (Section 3.3). Cohort mean 146 mg/L (range 13.4-309; Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (baseline and time-varying), not retained (Section 3.3). Cohort mean 87 umol/L (range 48-167; Table 1). eGFR (CKD-EPI), urea, LDH, ASAT, ALAT, GGT, total bilirubin, BMI, BSA and height were also screened and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 29L,
    n_observations = 139L,
    n_studies = 1L,
    age_range = "45-80 years",
    age_mean = "64 years",
    weight_range = "58-130 kg",
    weight_mean = "96.4 kg",
    sex_female_pct = 27.6,
    disease_state = "PCR-confirmed COVID-19 admitted to the ICU with respiratory organ support (82.8% mechanically ventilated; mean SOFA 7.5, APACHE IV 71.6). All co-treated with dexamethasone 6 mg once daily for up to 10 days; 24.1% also received methylprednisolone rescue.",
    dose_range = "Single intravenous dose of 8 mg/kg (maximum 800 mg) within 24 h of starting ICU organ support; administered doses 472-1552 mg (mean 781 mg; one patient accidentally received a double dose, one received the dose in two steps within 12 h).",
    regions = "Single centre, Leiden University Medical Center, The Netherlands (December 2020 - March 2021).",
    notes = paste(
      "Baseline characteristics from Moes 2022 Table 1.",
      "Free tocilizumab measured in EDTA plasma by ELISA (LLOQ 0.2 ug/mL) on 1-11 samples per patient (median 5) from day 1 to day 20 after dosing.",
      "Mean BMI 30.9 kg/m^2 (20.1-46.1); mean baseline CRP 146 mg/L; mean albumin 32 g/L."
    )
  )

  ini({
    # Structural parameters - Moes 2022 Table 3 final-model estimates; the same
    # values appear as the $THETA block of the final NONMEM control stream in
    # ESM 1. Time unit is days (CL in L/day, Vmax in mg/day); the $DES block
    # 'DADT(1) = -K*A(1) - (VM*C1)/(KM+C1)' with C1 = A(1)/V1 makes Vmax an
    # amount rate (mg/day) and Km a concentration (mg/L = ug/mL).
    lcl <- log(0.725); label("Linear clearance CL (L/day)") # Moes 2022 Table 3, CL (L/day) = 0.725
    lvc <- log(4.34); label("Volume of distribution Vd (L)") # Moes 2022 Table 3, Vd (L) = 4.34
    lvmax <- log(4.19); label("Maximum Michaelis-Menten elimination rate Vmax (mg/day)") # Moes 2022 Table 3, Vmax (mg/day) = 4.19
    lkm <- log(0.22); label("Michaelis-Menten constant Km (ug/mL)") # Moes 2022 Table 3, Km (mg/L) = 0.22

    # IIV - ESM 1 control stream $OMEGA BLOCK(2): 0.0351 (CL), 0.0355 (CL-V1
    # covariance), 0.043 (V1). The variances reproduce Table 3 exactly as
    # sqrt(exp(omega^2) - 1): 18.9% CV for CL and 21.0% CV for Vd. Table 3
    # does not print the covariance (correlation 0.914). No IIV on Vmax or Km
    # ('0 FIX' in the stream; Table 2 'No IIV identifiable for Km and Vmax').
    etalcl + etalvc ~ c(0.0351, 0.0355, 0.043)

    # Residual error - Table 3 / ESM 1 $ERROR: W = SQRT(THETA(5)^2*IPRED^2 +
    # THETA(6)^2) with $SIGMA 1 FIX, i.e. the combined2 form (nlmixr2 default).
    propSd <- 0.171; label("Proportional residual error (fraction)") # Moes 2022 Table 3, Proportional (CV%) = 17.1
    addSd <- 0.139; label("Additive residual error (ug/mL)") # Moes 2022 Table 3, Additive (ug/mL) = 0.139
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vmax <- exp(lvmax)
    km <- exp(lkm)

    # One-compartment IV disposition with parallel linear and Michaelis-Menten
    # elimination (Moes 2022 Fig. 1; ESM 1 $DES). central holds mg, so
    # Cc = central / vc is mg/L = ug/mL.
    Cc <- central / vc
    d/dt(central) <- -(cl / vc) * central - vmax * Cc / (km + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
