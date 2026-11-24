Peters_2016_ilofotaseAlfa <- function() {
  description <- paste(
    "Four-compartment (central + 3 peripheral) population PK model for",
    "ilofotase alfa (human recombinant alkaline phosphatase, recAP) given",
    "as 1-h IV infusions to healthy adult volunteers (Peters 2016), in a",
    "first-in-human single- (200-2000 U/kg) and multiple-ascending-dose",
    "(500 or 1000 U/kg once daily for 3 days) trial. Parameterised as the",
    "elimination and inter-compartmental rate constants plus the central",
    "volume; correlated IIV on the elimination rate constant and the",
    "central volume; proportional residual error; no covariates. Dosing",
    "is in mass units (ug); doses prescribed in U/kg must be converted",
    "(about 625 U per mg; see the vignette)."
  )
  reference <- paste(
    "Peters E, Heuberger JAAC, Tiessen R, van Elsas A, Masereeuw R,",
    "Arend J, Stevens J, Pickkers P.",
    "Pharmacokinetic Modeling and Dose Selection in a Randomized,",
    "Double-Blind, Placebo-Controlled Trial of a Human Recombinant",
    "Alkaline Phosphatase in Healthy Volunteers.",
    "Clin Pharmacokinet. 2016;55(10):1227-1237.",
    "doi:10.1007/s40262-016-0399-y.",
    sep = " "
  )
  vignette <- "Peters_2016_ilofotaseAlfa"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "ilofotase alfa", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "ilofotase alfa", units = "ug", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "ilofotase alfa", units = "ug", specimen = "tissue", verified = TRUE),
    peripheral3 = list(analyte = "ilofotase alfa", units = "ug", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age in years",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on V1 and K10 (Peters 2016 Methods, covariate analysis; Results, Pharmacokinetic Modeling). No candidate improved the fit, so none was retained."
    ),
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on V1 and K10 and not retained (Peters 2016 Results, Pharmacokinetic Modeling). Weight enters only through the U/kg dose."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened on V1 and K10 and not retained (Peters 2016 Results, Pharmacokinetic Modeling)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on V1 and K10 and not retained (Peters 2016 Results, Pharmacokinetic Modeling)."
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on V1 and K10 and not retained (Peters 2016 Results, Pharmacokinetic Modeling)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on V1 and K10 and not retained (Peters 2016 Results, Pharmacokinetic Modeling)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 36L,
    n_studies = 1L,
    age_range = "18-55 years (eligibility); active-arm group means 21-33 years",
    weight_range = "active-arm group means 69.4-73.8 kg (SD 8.2-16.3 kg)",
    sex_female_pct = 48.6,
    race_ethnicity = c(White = 94.6, Other = 5.4),
    disease_state = "healthy adult volunteers",
    dose_range = paste(
      "Single 1-h IV infusion of 200, 500, 1000 or 2000 U/kg (part A,",
      "6 PK-evaluable subjects per dose), or 1-h IV infusions of 500 or",
      "1000 U/kg once daily on 3 consecutive days (part B, 6 per dose)."
    ),
    regions = "Netherlands (single centre, Zuidlaren)",
    notes = paste(
      "Peters 2016 Tables 1 and 2. 37 subjects received recAP; one subject",
      "in the 200 U/kg group was dosed short because the infusion set",
      "leaked and was excluded from the PK analysis, leaving 36. Female",
      "fraction and race are for the 37 actively dosed subjects (Table 2:",
      "18/37 female, 35/37 white). EudraCT 2013-002694-21."
    )
  )

  ini({
    # Table 3 reports the rate constants and V1 as the estimated parameters;
    # CL, V2, V3 and V4 are derived from them (CL = K10 * V1, Vx = K1x * V1 / Kx1).
    # Paper compartment 2 -> peripheral1, 3 -> peripheral2, 4 -> peripheral3.
    lkel <- log(0.666); label("Elimination rate constant from the central compartment K10 (1/h)") # Table 3 K10 = 0.666 /h (RSE 4.97%)
    lvc <- log(3.71); label("Central volume of distribution V1 (L)") # Table 3 V1 = 3.71 L (RSE 5.23%)
    lk12 <- log(0.311); label("Central to peripheral1 rate constant K12 (1/h)") # Table 3 K12 = 0.311 /h (RSE 11.6%)
    lk21 <- log(0.0143); label("Peripheral1 to central rate constant K21 (1/h)") # Table 3 K21 = 0.0143 /h (RSE 7.64%)
    lk13 <- log(0.545); label("Central to peripheral2 rate constant K13 (1/h)") # Table 3 K13 = 0.545 /h (RSE 14.1%)
    lk31 <- log(0.571); label("Peripheral2 to central rate constant K31 (1/h)") # Table 3 K31 = 0.571 /h (RSE 20.7%)
    lk14 <- log(0.564); label("Central to peripheral3 rate constant K14 (1/h)") # Table 3 K14 = 0.564 /h (RSE 9%)
    lk41 <- log(0.0795); label("Peripheral3 to central rate constant K41 (1/h)") # Table 3 K41 = 0.0795 /h (RSE 11.3%)

    # IIV reported as %CV = sqrt(exp(omega^2) - 1) * 100 (ESM formula image).
    # K10 16.51% -> 0.026893; V1 28.65% -> 0.078887 (Table 3).
    # Covariance is not printed; the text says it was estimated. It is
    # back-solved from Table 3's derived IIV on CL = K10 * V1 (18.48% ->
    # 0.033581): cov = (0.033581 - 0.026893 - 0.078887) / 2 = -0.036100
    # (correlation -0.78).
    etalkel + etalvc ~ c(0.026893, -0.036100, 0.078887) # Table 3 IIV K10 16.51 pct CV, V1 28.65 pct CV, CL 18.48 pct CV

    propSd <- sqrt(0.02); label("Proportional residual error (fraction)") # Table 3 residual error variance 0.02 (printed as omega^2)
  })

  model({
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    k14 <- exp(lk14)
    k41 <- exp(lk41)

    d/dt(central) <- -(kel + k12 + k13 + k14) * central +
      k21 * peripheral1 + k31 * peripheral2 + k41 * peripheral3
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(peripheral3) <- k14 * central - k41 * peripheral3

    # central in ug, vc in L -> ug/L = ng/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
