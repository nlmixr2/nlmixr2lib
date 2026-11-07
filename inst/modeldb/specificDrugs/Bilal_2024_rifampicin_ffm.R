Bilal_2024_rifampicin_ffm <- function() {
  description <- paste(
    "One-compartment population PK model for a single 600 mg oral dose of",
    "rifampicin in healthy Caucasian adult volunteers (Bilal 2024), the",
    "fat-free-mass (FFM) covariate model that the authors ranked best.",
    "Absorption is a zero-order input of duration Tk0 (0.470 h) into the",
    "central compartment after a lag time (0.340 h); elimination is",
    "Michaelis-Menten (Vmax/F 190 mg/h, Km 20.1 mg/L) from a central volume",
    "V/F of 36.2 L. FFM scales V/F linearly and Vmax/F with a fixed 0.75",
    "exponent, both normalised to the cohort mean FFM of about 50.1 kg.",
    "Between-subject variability is on Vmax/F only (18.5%); inter-occasion",
    "variability over the two crossover periods is on lag time, Tk0, V/F and",
    "Vmax/F, with correlated lag-time and Tk0 occasion effects. Residual",
    "error is the Monolix combined1 form (additive 0.060 mg/L plus 11.0%",
    "proportional, summed as standard deviations). Dose records into",
    "central must carry rate = -2 so the modelled zero-order duration",
    "applies."
  )
  reference <- paste(
    "Bilal M, Ullah S, Jaehde U, Trueck C, Zaremba D, Wachall B, Wargenau M,",
    "Scheidel B, Wiesen MHJ, Gazzaz M, Chen C, Buesker S, Fuhr U, Taubert M,",
    "Dokos C. Assessment of body mass-related covariates for rifampicin",
    "pharmacokinetics in healthy Caucasian volunteers. Eur J Clin Pharmacol.",
    "2024;80(9):1271-1283. doi:10.1007/s00228-024-03697-3.",
    "Parameter estimates from Table 3 (FFM covariate model column); the",
    "structural model and random-effect structure from the Monolix code of",
    "the base model in the supplementary material. Fat-free mass follows",
    "Janmahasatian et al. (2005) Clin Pharmacokinet 44(10):1051-1065",
    "doi:10.2165/00003088-200544100-00004.",
    sep = " "
  )
  vignette <- "Bilal_2024_rifampicin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "rifampicin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass computed from body weight, BMI and sex by the Janmahasatian (2005) formula",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power model normalised to the observation-weighted mean FFM",
        "(Methods, Covariate analysis: 'Continuous covariates were modeled",
        "using power models normalized by weighted means'; Eq. 1). The",
        "exponents are fixed: 1 on V/F and 0.75 on Vmax/F (Table 3).",
        "FFM is computed as in Eqs. 2-3: male",
        "FFM = 9.27e3 * WT / (6.68e3 + 216 * BMI), female",
        "FFM = 9.27e3 * WT / (8.78e3 + 244 * BMI), with WT in kg and BMI in",
        "kg/m^2. The weighted-mean reference is not printed. Every subject",
        "contributed 38 observations (912 / 24; 19 post-dose samples in each",
        "of two periods), so it equals the arithmetic mean FFM, which is also",
        "not printed. 50.1 kg is the sex-composition-weighted average of the",
        "per-sex Table 1 medians, (13 * 39.2 + 11 * 63) / 24. Cross-check:",
        "with this reference the FFM model gives a typical female (FFM",
        "39.2 kg) Vmax/F of 190 * (39.2 / 50.1)^0.75 = 158 mg/h, matching the",
        "female-reference Vmax/F of 157 mg/h in the sex + body weight model",
        "(Bilal_2024_rifampicin_sexwt); matching it exactly would need a",
        "reference of 50.6 kg. Cohort FFM median (range)",
        "39.2 (33.7-44.8) kg in females and 63 (47.09-68.9) kg in males; the",
        "two ranges do not overlap (Table 1)."
      ),
      source_name = "FFM"
    ),
    OCC = list(
      description = "Crossover period (1 or 2), the occasion for the inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "The two periods of the bioequivalence crossover were treated as two",
        "occasions (Methods, Population PK analysis); each period is one",
        "single 600 mg dose and its 24-h sampling. Set OCC = 1 for the first",
        "period and OCC = 2 for the second. The Monolix code declares the",
        "variability levels id and id*occ."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Not a covariate of this model, but required upstream to compute FFM (Eqs. 2-3). BW on V/F and Vmax/F lowered the OFV by 37.1 points (Table 4, model 1), less than FFM. BW on V/F is retained in the alternative sex + body weight model, Bilal_2024_rifampicin_sexwt. Median (range) 63.8 (48.5-73.1) kg in females and 82.9 (57.8-91.0) kg in males (Table 1)."
    ),
    SEXF = list(
      description = "1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "Not a covariate of this model, but selects the sex-specific Janmahasatian coefficients when computing FFM. Adding sex in addition to FFM did not improve the model (Results, Covariate modeling). Sex is retained in the alternative Bilal_2024_rifampicin_sexwt model. 13 of 24 subjects female."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Input to the Janmahasatian FFM formula. Significant on Tlag (dOFV -4.02) and V/F (dOFV -5.11) in forward inclusion but failed backward elimination (Results). Median (range) 23.09 (18.7-26.7) in females and 24.7 (20.1-29.7) in males (Table 1)."
    ),
    HT = list(
      description = "Body height.",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "BH on V/F and Vmax/F lowered the OFV by 42.9 points (Table 4, model 3), less than FFM. Reported in metres: median (range) 1.65 (1.55-1.76) m in females and 1.81 (1.63-1.89) m in males (Table 1)."
    ),
    BSA = list(
      description = "Body surface area by the Mosteller formula.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "BSA on V/F and Vmax/F lowered the OFV by 48.2 points (Table 4, model 5), less than FFM. Median (range) 1.66 (1.47-1.85) m^2 in females and 2.05 (1.62-2.18) m^2 in males (Table 1)."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Methods, Covariate analysis) and not retained. Median (range) 37.0 (21.0-58.0) years in females and 43.0 (22.0-64.0) years in males (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "21-64 years (median 39.5)",
    weight_range = "48.5-91.0 kg (median 68)",
    ffm_range = "33.7-68.9 kg (female median 39.2, male median 63)",
    sex_female_pct = 54.2,
    race_ethnicity = "Caucasian (100%)",
    disease_state = "Healthy volunteers, BMI 18.5-30 kg/m^2, non-smokers",
    dose_range = "Single oral 600 mg rifampicin tablet in each of two periods (test tablet and the EREMFAT reference tablet), separated by a washout of at least 9 days",
    regions = "Germany (Cologne)",
    notes = paste(
      "Phase I/IV randomised, two-period crossover, open-label",
      "bioequivalence study (EudraCT 2017-004418-24). 25 enrolled, one",
      "dropped out before the first dose; 24 completed (11 men, 13 women).",
      "Plasma sampled predose and at 0.16, 0.33, 0.5, 0.75, 1, 1.25, 1.5,",
      "1.75, 2, 2.25, 2.5, 3, 3.5, 4, 6, 9, 12, 16 and 24 h. 912",
      "concentrations, 99 (10.8%) below the 0.1 mg/L LLOQ and treated as",
      "interval-censored. The two tablets were bioequivalent by NCA and the",
      "preparation had no significant effect on any parameter, so the data",
      "were pooled. Fit in Monolix 2023R1 (SAEM). Demographics are Table 1."
    )
  )

  ini({
    # Structural parameters: Table 3, 'FFM covariate model' estimates column.
    ltlag <- log(0.340)
    label("Absorption lag time Tlag (h)") # Table 3 FFM model Tlag = 0.340 h (RSE 5.53%)
    ld1 <- log(0.470)
    label("Duration of the zero-order absorption Tk0 (h)") # Table 3 FFM model Tk0 = 0.470 h (RSE 8.99%)
    lvc <- log(36.2)
    label("Apparent volume of distribution V/F at the reference FFM of 50.1 kg (L)") # Table 3 FFM model V/F = 36.2 L (RSE 2.21%)
    lvmax <- log(190)
    label("Apparent maximum elimination rate Vmax/F at the reference FFM of 50.1 kg (mg/h)") # Table 3 FFM model Vmax/F = 190 mg/h (RSE 4.39%)
    lkm <- log(20.1)
    label("Michaelis-Menten constant Km (mg/L)") # Table 3 FFM model Km = 20.1 mg/L (RSE 2.44%)

    # Covariate exponents, both fixed in the final model.
    e_ffm_vc <- fixed(1)
    label("Power exponent of FFM on V/F (unitless)") # Table 3 'beta FFM on V/F' = 1.00 (fixed)
    e_ffm_vmax <- fixed(0.75)
    label("Power exponent of FFM on Vmax/F (unitless)") # Table 3 'FFM on Vmax/F' = 0.750 (fixed)

    # Monolix reports random effects as omega, the SD of the log-normal eta;
    # Table 3 prints them x 100 under a 'CV %' heading. Variances are omega^2.
    etalvmax ~ 0.034225 # Table 3 FFM model IIV Vmax/F = 18.5; 0.185^2

    # Inter-occasion variability, two occasions (crossover periods). Monolix
    # gives one gamma per parameter, shared by both occasions; occasion 2
    # repeats the occasion-1 value as fixed. Lag-time and Tk0 occasion
    # effects are correlated: cov = 0.410 * 0.390 * 0.655 = 0.1047345.
    etaiov_tlag_1 + etaiov_d1_1 ~ c(
      0.1521,
      0.1047345,
      0.429025
    ) # Table 3 FFM model IOV Tlag = 39.0, IOV Tk0 = 65.5, Corr. IOV Tlag and IOV Tk0 = 0.410
    etaiov_tlag_2 + etaiov_d1_2 ~ fixed(c(
      0.1521,
      0.1047345,
      0.429025
    ))
    etaiov_vc_1 ~ 0.022201 # Table 3 FFM model IOV V = 14.9; 0.149^2
    etaiov_vc_2 ~ fixed(0.022201)
    etaiov_vmax_1 ~ 0.00881721 # Table 3 FFM model IOV Vmax = 9.39; 0.0939^2
    etaiov_vmax_2 ~ fixed(0.00881721)

    # Residual error: Monolix combined1, SD = a + b * f.
    addSd <- 0.060
    label("Additive residual error (mg/L)") # Table 3 FFM model additive residual error = 0.060 mg/L (RSE 7.69%)
    propSd <- 0.110
    label("Proportional residual error (fraction)") # Table 3 FFM model proportional residual error = 11.0% (RSE 3.31%)
  })
  model({
    # Occasion indicators: crossover period 1 or 2.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2
    iov_d1 <- oc1 * etaiov_d1_1 + oc2 * etaiov_d1_2
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2
    iov_vmax <- oc1 * etaiov_vmax_1 + oc2 * etaiov_vmax_2

    # Eq. 1 power model, centred on the cohort mean FFM (see covariateData).
    tlag <- exp(ltlag + iov_tlag)
    d1 <- exp(ld1 + iov_d1)
    vc <- exp(lvc + iov_vc) * (FFM / 50.1)^e_ffm_vc
    vmax <- exp(lvmax + etalvmax + iov_vmax) * (FFM / 50.1)^e_ffm_vmax
    km <- exp(lkm)

    # Figure 1: one compartment with Michaelis-Menten elimination
    # Vmax * Cp / (Km + Cp); the Monolix pkmodel(Tlag, Tk0, V, Vm, Km).
    Cc <- central / vc
    d/dt(central) <- -vmax * Cc / (km + Cc)

    # Zero-order input of duration Tk0 starting Tlag after the dose; dose
    # records need rate = -2.
    alag(central) <- tlag
    dur(central) <- d1

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
