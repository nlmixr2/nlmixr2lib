Zavrelova_2025_linezolid <- function() {
  description <- "One-compartment population PK model for intravenous linezolid in hematooncological patients with suspected or proven Gram-positive sepsis (Zavrelova 2025); clearance and volume of distribution both decline exponentially with age, and clearance is 33% lower on day 4 of therapy than on day 1"
  reference <- "Zavrelova A, Merdita S, Zak P, Radocha J, Visek B, Lanska M, Malakova J, Michalek P, Slanar O, Sima M. Population Pharmacokinetic Model-Based Optimization of Linezolid Dosing in Hematooncological Patients With Suspected or Proven Gram-Positive Sepsis. Clin Transl Sci. 2025;18:e70346. doi:10.1111/cts.70346"
  vignette <- "Zavrelova_2025_linezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Linezolid was assayed in serum by LC-MS/MS (Methods 2.1).
  compartmentData <- list(
    central = list(analyte = "linezolid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Entered UNCENTRED on the log scale, exactly as printed in the source: Vd = 99.44 * exp(-0.013 * AGE) and CL = 46.93 * exp(-0.023 * AGE). Vd_pop and CL_pop are therefore the extrapolated age-0 intercepts, not values at a reference age. Cohort median 59 years (IQR 48-70, range 26-82; Table 1), at which Vd = 46.2 L and CL = 12.1 L/h as stated in Results 3.2.",
      source_name        = "Age"
    ),
    DAY4 = list(
      description        = "Day-4-of-therapy landmark indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (day 1 of therapy)",
      notes              = "Within-subject step indicator entered as a regressor (time-varying covariate) in Monolix: 0 for day-1 observations, 1 for day-4 observations. Gates the exp(-0.40) = 0.67 fold change in CL reported in Table 2, i.e. the 33% CL reduction between days 1 and 4 (Results 3.2, Discussion). The source samples only days 1 and 4 and states it could not localise when the change occurs (Discussion: 'the described reduction in CL probably occurs' by about day 2), so the cutoff is a sampling landmark, not a mechanistic threshold.",
      source_name        = "4th day"
    )
  )

  # Screened in the covariate analysis (Methods 2.3) but not retained in the
  # final model -- "The other tested variables exerted no influence on the
  # linezolid pharmacokinetics" (Results 3.2). Documented here for provenance;
  # these names are deliberately absent from model() and from covariateData.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      notes = "Median 84 kg (IQR 74-91, range 60-122; Table 1). Not correlated with Vd (p = 0.425) or CL (p = 0.513) in the graphical screen, and adding it did not significantly decrease OFV (Discussion). The nomogram-validation sub-set (Figure S6) likewise showed no trend of serum level on body weight."
    ),
    HT = list(
      description = "Body height", units = "cm", type = "continuous",
      notes = "Median 176 cm (IQR 167-185, range 158-190; Table 1). Screened, not retained."
    ),
    BSA = list(
      description = "Body surface area", units = "m^2", type = "continuous",
      notes = "Median 2.02 m^2 (IQR 1.85-2.14, range 1.62-2.44; Table 1). Screened, not retained."
    ),
    CREAT = list(
      description = "Serum creatinine", units = "umol/L", type = "continuous",
      notes = "Median 79 umol/L (IQR 56-123, range 36-459; Table 1). Correlation with linezolid CL was not significant (p = 0.162; Discussion)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI)", units = "mL/s", type = "continuous",
      notes = "Median 1.53 mL/s (IQR 0.88-1.76, range 0.16-2.27; Table 1). Correlated with CL in the graphical screen but more weakly than age (p = 0.0099 vs p = 0.0000125); once age entered the model the eGFR-CL relationship disappeared, which the authors attribute to age being a component of the CKD-EPI formula (Discussion). Reported in mL/s, not the mL/min/1.73 m^2 more usual for this canonical."
    ),
    TBILI = list(
      description = "Serum total bilirubin", units = "umol/L", type = "continuous",
      notes = "Median 12 umol/L (IQR 9-16, range 2-112; Table 1). Screened, not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "ukat/L", type = "continuous",
      notes = "Median 0.40 ukat/L (IQR 0.31-0.55, range 0.11-2.02; Table 1). Screened, not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "ukat/L", type = "continuous",
      notes = "Median 0.50 ukat/L (IQR 0.34-1.02, range 0.13-1.78; Table 1). Screened, not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase", units = "ukat/L", type = "continuous",
      notes = "Median 2.02 ukat/L (IQR 1.36-2.77, range 0.73-7.50; Table 1). Screened, not retained."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase", units = "ukat/L", type = "continuous",
      notes = "Median 2.81 ukat/L (IQR 1.05-4.09, range 0.35-13.91; Table 1). Screened, not retained."
    ),
    INR_BASE = list(
      description = "International normalized ratio", units = "(unitless)", type = "continuous",
      notes = "Median 1.15 (IQR 1.04-1.20, range 0.92-1.66; Table 1). Screened, not retained."
    ),
    SEXF = list(
      description = "Female sex indicator", units = "(binary)", type = "binary",
      notes = "5 of 22 patients female (22.7%; Table 1). Tested as a categorical covariate (Methods 2.3), not retained."
    ),
    NEUT = list(
      description = "Absolute neutrophil count", units = "cells/mm^3", type = "continuous",
      notes = "Tested as the binary covariate 'neutropenia', defined as an absolute neutrophil count < 500 cells/mm^3 (Methods 2.1, 2.3). Not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22,
    n_studies      = 1,
    n_observations = 197,
    age_range      = "26-82 years",
    age_median     = "59 years",
    weight_range   = "60-122 kg",
    weight_median  = "84 kg",
    sex_female_pct = 22.7,
    race_ethnicity = "Single-centre Czech cohort; race/ethnicity not reported in the source.",
    disease_state  = "Adult hematooncological patients (most frequently acute leukaemia, n = 14, and lymphoma, n = 3) treated in a haematology intensive care unit for suspected or proven Gram-positive sepsis; bloodstream infection documented in 12 of 22 (55%). Patients receiving renal replacement therapy were excluded.",
    renal_function = "eGFR (CKD-EPI) median 1.53 mL/s (IQR 0.88-1.76, range 0.16-2.27); serum creatinine median 79 umol/L (range 36-459).",
    dose_range     = "Intravenous linezolid; 19 of 22 patients received the off-label high dose of 600 mg every 8 h as a 4-h infusion, and 3 patients received 600 mg every 6 h as a 4-h infusion.",
    regions        = "Czech Republic (University Hospital Hradec Kralove).",
    notes          = "Retrospective cross-sectional therapeutic-drug-monitoring analysis, July 2023 to September 2024. Samples were drawn at 1, 2, 4, 6 and 8 h after the start of the first infusion (day 1) and again on day 4 (20 of 22 patients); 4-10 concentrations per patient, 9 on average. Seven trough levels were excluded because the next infusion had already started. Demographics in Table 1. A separate prospective validation sub-set of 40 patients dosed by the age-based nomogram (Table 4) is described in Results 3.4 but was not used to fit this model."
  )

  ini({
    # Structural typical values. Age enters UNCENTRED (see covariateData$AGE),
    # so these are the extrapolated age-0 intercepts of the source equations
    # Vd = 99.44 * exp(-0.013 * AGE) and CL = 46.93 * exp(-0.023 * AGE),
    # not values at a reference age.
    lvc <- log(99.44); label("Volume of distribution, age-0 intercept of the uncentred age model (Vd, L)")  # Table 2 (Vd_pop = 99.44, RSE 23.1%)
    lcl <- log(46.93); label("Clearance, age-0 intercept of the uncentred age model (CL, L/h)")             # Table 2 (CL_pop = 46.93, RSE 25.2%)

    # Covariate effects, all on the natural-log scale of the parameter.
    e_age_vc  <- -0.013; label("Effect of age on log Vd (1/year)")                    # Table 2 (beta_Vd_Age = -0.013, RSE 28.8%)
    e_age_cl  <- -0.023; label("Effect of age on log CL (1/year)")                    # Table 2 (beta_CL_Age = -0.023, RSE 18.4%)
    e_day4_cl <- -0.40;  label("Effect of day 4 of therapy on log CL (unitless)")     # Table 2 (beta_CL_4th day = -0.40, RSE 13.3%); exp(-0.40) = 0.67, the reported 33% CL reduction

    # IIV. Table 2 reports Omega as the standard deviation of the normally
    # distributed random effect on the log scale, so the variance is Omega^2.
    etalvc       ~ 0.23^2  # Table 2 (Omega_Vd = 0.23, RSE 21.2%)
    etalcl       ~ 0.26^2  # Table 2 (Omega_CL = 0.26, RSE 18.4%)
    # Between-subject variability on the day-4 effect itself: Table 2 reports a
    # fixed effect (beta_CL_4th day) and a matching random-effect SD
    # (Omega_CL_4th day) for the same term, which is how Monolix reports an
    # individual parameter that multiplies a regressor. See the vignette
    # "Assumptions and deviations" section for the reading and the alternative.
    etae_day4_cl ~ 0.21^2  # Table 2 (Omega_CL_4th day = 0.21, RSE 20.7%)

    # Combined residual error (Results 3.2: "a combined error model was found to
    # best fit ... and explain the residual variability").
    addSd  <- 0.63;  label("Additive residual error (mg/L)")            # Table 2 (Constant = 0.63, RSE 20.2%)
    propSd <- 0.052; label("Proportional residual error (fraction)")    # Table 2 (Proportional = 0.052, RSE 33.3%)
  })
  model({
    # Individual parameters. Source equations (Results 3.2):
    #   Log(Vd) = log(Vd_pop) + beta_Vd_Age * Age + eta_Vd
    #   Log(CL) = log(CL_pop) + beta_CL_Age * Age + beta_CL_4th day [if 4th day] + eta_CL
    vc <- exp(lvc + e_age_vc * AGE + etalvc)
    cl <- exp(lcl + e_age_cl * AGE + (e_day4_cl + etae_day4_cl) * DAY4 + etalcl)

    kel <- cl / vc

    # One compartment with linear elimination; linezolid is given intravenously,
    # so there is no absorption compartment and F is structurally 1.
    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L (= ug/mL)
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
