Li_2023_granisetron <- function() {
  description <- "One-compartment population PK model for the granisetron transdermal delivery system (GTDS, Sancuso 34.3 mg/52 cm2 patch) in healthy Caucasian volunteers, with first-order absorption from the patch (no lag time) and first-order elimination; no covariates retained (Li 2023)"
  reference <- paste(
    "Li J, Hu P, Zhou L, Nagahama F, Chen R.",
    "Population pharmacokinetic analysis of transdermal granisetron in",
    "healthy Chinese and Caucasian volunteers.",
    "Front Pharmacol. 2023;14:1154026.",
    "doi:10.3389/fphar.2023.1154026"
  )
  vignette <- "Li_2023_granisetron"
  # Li 2023 Table 3 reports V in mL and CL in mL/h with concentrations in ng/mL
  # (Figures 2-4 axis labels). Encoding dose in ug with vc in L reproduces the
  # authors' concentration scale exactly, because 1 ug/L == 1 ng/mL. A 34.3 mg
  # patch is therefore dosed as amt = 34300.
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    # Li 2023 Section 2.2 + 3.2: age, weight, height, BMI and sex were each
    # screened by stepwise forward inclusion (dOFV > 3.84, p < 0.05) and
    # backward elimination (dOFV > 6.64, p < 0.01). None was retained, so the
    # final model carries no covariates. Recorded here to preserve the
    # provenance of the covariate screen.
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a power function normalised to the population median (Li 2023 Eq 2) but not retained in the final model (Li 2023 Section 3.2).",
      source_name        = "age"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a power function normalised to the population median (Li 2023 Eq 2) but not retained in the final model (Li 2023 Section 3.2). No allometric scaling is applied a priori either.",
      source_name        = "weight"
    ),
    HT = list(
      description        = "Height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a power function normalised to the population median (Li 2023 Eq 2) but not retained in the final model (Li 2023 Section 3.2).",
      source_name        = "height"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a power function normalised to the population median (Li 2023 Eq 2) but not retained in the final model (Li 2023 Section 3.2).",
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = "Screened as an exponential categorical effect (Li 2023 Eq 3) but not retained in the final model (Li 2023 Section 3.2). Li 2023 Discussion notes the same negative finding in Howell 2009 (n = 48 healthy subjects + 793 cancer patients).",
      source_name        = "gender"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 112L,                               # Li 2023 Section 3.1: 1372 concentrations from 112 Caucasian healthy subjects used for model building
    n_studies      = 4L,                                 # Li 2023 Table 1: 392MD/11/C, 392MD/26/C, 392MD/40/C, 392MD/43/C
    n_observations = 1372L,                              # Li 2023 Section 3.1
    age_range      = NULL,                               # Li 2023 Table 2 reports mean +/- SD per study, not ranges
    age_mean       = "43.01 +/- 17.80 years",            # Li 2023 Section 3.1 (pooled Caucasian)
    weight_mean    = "70.54 +/- 16.13 kg",               # Li 2023 Section 3.1 (pooled Caucasian)
    height_mean    = "169.51 +/- 9.69 cm",               # Li 2023 Section 3.1 (pooled Caucasian)
    bmi_mean       = "24.39 +/- 4.43 kg/m^2",            # Li 2023 Section 3.1 (pooled Caucasian)
    sex_female_pct = 42.86,                              # Li 2023 Section 3.1 (pooled Caucasian)
    race_ethnicity = c(White = 100, Black = 0, Asian = 0, Other = 0), # Li 2023 Table 2: all four model-building studies enrolled Caucasian healthy volunteers
    disease_state  = "Healthy volunteers (no chemotherapy-induced nausea and vomiting; GTDS is indicated for CINV prophylaxis).",
    dose_range     = "Single 34.3 mg/52 cm2 granisetron transdermal patch worn 6-9 days depending on study (Li 2023 Table 1).",
    regions        = "Four Caucasian healthy-volunteer studies (sponsor studies 392MD/11/C, 392MD/26/C, 392MD/40/C, 392MD/43/C).",
    sampling       = "Study-specific schedules from pre-dose out to 120-216 h post-application; 96-660 samples per study, 1372 in total (Li 2023 Table 1).",
    software       = "Phoenix NLME 1.30, FOCE ELS (Li 2023 Section 2.2).",
    notes          = paste(
      "The model was fit to the four Caucasian studies only. A fifth study",
      "(SP-0102, 24 Chinese healthy male volunteers, mean age 27.13 +/- 4.07",
      "years, weight 65.07 +/- 5.67 kg, sampling to 240 h) was held out and used",
      "as the external comparison cohort: the Caucasian model was simulated",
      "under the Chinese trial design and the observed Chinese concentrations",
      "were compared against the 5th-95th percentile prediction interval",
      "(Li 2023 Sections 2.3.2 and 3.4, Table 4). No dose adjustment for the",
      "Chinese population was indicated."
    )
  )

  ini({
    # Structural parameters - Li 2023 Table 3, "Fixed effect" block.
    # CL and V are apparent (CL/F, V/F): the delivered fraction of the 34.3 mg
    # patch content is not identifiable from plasma data alone, so no
    # bioavailability term is estimated and the nominal patch content is dosed.
    lka <- log(0.0179879); label("Absorption rate constant from the patch (1/h)")  # Li 2023 Table 3: Ka = 0.0179879 1/h (RSE 3.63%); bootstrap mean 0.017966084, 95% CI 0.016469625-0.019216464. Printed as "0.0179,879" in Table 3 because Frontiers' typesetting inserts a thousands separator; the bootstrap mean and CI confirm 0.0179879.
    lcl <- log(31.3163);   label("Apparent clearance CL/F (L/h)")                  # Li 2023 Table 3: CL = 31316.3 mL/h = 31.3163 L/h (RSE 5.52%); bootstrap mean 31042.673 mL/h, 95% CI 26407.14-36914.471 mL/h. Abstract: "apparent systemic clearance was determined to be 31316.3 mL/h".
    lvc <- log(6299.03);   label("Apparent central volume of distribution V/F (L)") # Li 2023 Table 3: V = 6299030 mL = 6299.03 L (RSE 13.12%); bootstrap mean 6351196.5 mL, 95% CI 4962557.7-8105176.7 mL. Abstract: "the central compartment volume of distribution was 6299.03 L".

    # IIV - Li 2023 Table 3, "Random effect" block; exponential (log-normal) IIV
    # per Li 2023 Eqs 8-10 (ka = tvka * exp(eta_ka), V = tvV * exp(eta_V),
    # CL = tvCL * exp(eta_CL)). Li 2023 Section 2.2 states eta is drawn from a
    # normal distribution "with mean zero and a variance omega^2", so the
    # tabulated random-effect values are variances and are used unchanged.
    etalka ~ 4.74078e-11 # Li 2023 Table 3: Random effect Ka = 4.74078E-11 (RSE 149.11%); bootstrap 95% CI 0.00000-0.00000. Estimated but numerically indistinguishable from zero, i.e. no between-subject variability in patch absorption rate was supported.
    etalcl ~ 0.25423392  # Li 2023 Table 3: Random effect CL = 0.25423392 (RSE 20.19%); bootstrap 95% CI 0.1536138-0.3548540. Equivalent to 53.8% CV.
    etalvc ~ 1.5462161   # Li 2023 Table 3: Random effect V = 1.5462161 (RSE 12.11%); bootstrap 95% CI 1.179257-1.9131749. Equivalent to 192% CV.

    # Residual error - Li 2023 Section 3.2 and Eq 7 (Cobs = C + epsilon): the
    # additive model was selected over proportional and mixed models.
    addSd <- 1.18094; label("Additive residual error (ng/mL)") # Li 2023 Table 3: sigma = 1.18094 (RSE 10.17%); bootstrap mean 1.1706386, 95% CI 0.94536484-1.4335495.
  })
  model({
    # Individual parameters - Li 2023 Eqs 8-10
    ka  <- exp(lka + etalka)
    cl  <- exp(lcl + etalcl)
    vc  <- exp(lvc + etalvc)

    kel <- cl / vc

    # ODE system - Li 2023 Eqs 4-5 and Figure 1
    d/dt(depot)   <- -ka * depot            # Eq 4: dAa/dt = -ka * Aa
    d/dt(central) <- ka * depot - kel * central # Eq 5: dA1/dt = -CL * C + ka * Aa

    # Observation - Li 2023 Eq 6 (C = A1/V) and Eq 7 (Cobs = C + epsilon).
    # With amt in ug and vc in L, central/vc is in ug/L == ng/mL.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
