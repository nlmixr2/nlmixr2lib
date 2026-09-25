Larizza_2018_cpeptideKinetics_M0 <- function() {
  description <- paste(
    "Population regression model (Van Cauter 1992 / Magni 2000 'M0') predicting the four",
    "macro constants of C-peptide two-compartment kinetics -- short half-life, amplitude",
    "fraction, long half-life and central distribution volume -- from health status",
    "(normal / obese / diabetic), sex, body surface area and age in 207 adults. This is the",
    "maximum-likelihood (NONMEM) variant with independent additive residual errors on each",
    "macro constant and no estimated interindividual variability; the Bayesian variant with a",
    "full 4x4 interindividual covariance matrix is Larizza_2018_cpeptideKinetics_M4. The macro",
    "constants convert algebraically to the micro rate constants k01, k12 and k21 that",
    "Larizza_2018_insulinMinimalModel requires. Re-estimated by Larizza 2018 via the DDMoRe",
    "Interoperability Framework; the estimates are identical to the original publication.",
    sep = " "
  )
  reference <- paste(
    "Larizza C, Borella E, Pasotti L, Tartaglione P, Smith M, Moodie S, Magni P. (2018).",
    "Complex Bayesian Modeling Workflows Encoding and Execution Made Easy With a Novel",
    "WinBUGS Plugin of the Drug Disease Model Resources Interoperability Framework.",
    "CPT Pharmacometrics Syst Pharmacol 7(5):298-308. doi:10.1002/psp4.12285.",
    "Underlying models: Van Cauter E, Mestrez F, Sturis J, Polonsky KS. Diabetes 1992;41:368-377;",
    "Magni P, Bellazzi R, Sparacino G, Cobelli C. Ann Biomed Eng 2000;28:812-823.",
    sep = " "
  )
  vignette <- "Larizza_2018_cpeptide_insulin_secretion"

  units <- list(
    time = "min",
    dosing = "n/a (covariate regression model; no drug input)",
    concentration = "n/a (outputs are C-peptide kinetic macro constants: min, unitless and L)"
  )

  covariateData <- list(
    DIS_OBESE = list(
      description = "Obese health-status stratum indicator (1 = obese, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal health status, when DIS_DIAB is also 0)",
      notes = paste(
        "One of the three mutually exclusive levels of the source paper's HSTATUS column",
        "(normal / obese / diabetic). Normal is the reference level and is encoded as",
        "DIS_OBESE = 0 and DIS_DIAB = 0; the normal-stratum indicator is recovered inside",
        "model() as (1 - DIS_OBESE - DIS_DIAB). Exactly one of the three levels must hold per",
        "subject. The source cohort is the 207-subject Van Cauter 1992 dataset, whose obese",
        "group is defined by body habitus rather than by a stated BMI threshold, so this is",
        "deliberately NOT the morbidly-obese canonical DIS_OBESE_MORBID.",
        sep = " "
      ),
      source_name = "HSTATUS (level 'obese')"
    ),
    DIS_DIAB = list(
      description = "Diabetic health-status stratum indicator (1 = diabetic, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal health status, when DIS_OBESE is also 0)",
      notes = paste(
        "Second non-reference level of the source paper's HSTATUS column. Van Cauter 1992",
        "enrolled non-insulin-dependent (type 2) diabetic subjects; per the DIS_DIAB register",
        "entry the type-2-specific detail is recorded here rather than on a parallel canonical.",
        sep = " "
      ),
      source_name = "HSTATUS (level 'diabetic')"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Selects between the male and female intercept / BSA-slope pair of the volume",
        "regression. Larizza 2018 Mathematical Models section: theta_V = aVm + bVm * BSA if SEX",
        "is male, aVf + bVf * BSA if SEX is female. Both strata are estimated, so neither is an",
        "offset from the other.",
        sep = " "
      ),
      source_name = "SEX"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Uncentred linear predictor of the C-peptide central volume. Larizza 2018 computes it",
        "with the Du Bois-family formula stated in the Mathematical Models section:",
        "BSA = 0.20247 * Height(m)^0.725 * Weight(kg)^0.425. Supply BSA directly; the model",
        "does not recompute it from height and weight.",
        sep = " "
      ),
      source_name = "BSA"
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Uncentred linear predictor of the C-peptide long half-life:",
        "theta_tl = atl + btl * AGE. Not centred in the source, so the intercept atl is the",
        "long half-life extrapolated to age 0 and is not itself interpretable.",
        sep = " "
      ),
      source_name = "AGE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 207,
    n_studies = 1,
    age_range = "not reported (adults)",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    disease_state = "normal, obese, or diabetic (type 2) health status",
    dose_range = "n/a (no drug administered; biosynthetic C-peptide kinetic study)",
    regions = "not reported",
    notes = paste(
      "Larizza 2018 Datasets section: 'a large dataset, including information about health",
      "status, sex, age, BSA, and corresponding CP kinetics macro constants of 207 subjects',",
      "citing Van Cauter 1992 and Magni 2000. The individual macro constants were themselves",
      "obtained in the source studies by fitting a two-compartment model to biosynthetic",
      "human C-peptide decay curves. Per-stratum counts, age and weight ranges are not",
      "reported in Larizza 2018.",
      sep = " "
    )
  )

  ini({
    # ---- Short half-life ts, per health-status stratum (Larizza 2018 Table 1, M0 column) ----
    # Natural scale, NOT log-transformed: the source regression is additive-normal on the
    # original scale (theta = mu_stratum, phi_i = theta_i + gamma_i), so a log
    # parameterisation would misstate the structure the random effects sit on.
    thalfShort_normal <- 5.000 ; label("Typical C-peptide short half-life, normal health status (min)") # Table 1 M0: mtsn = 5.000 min, 2.088 %RSE
    thalfShort_obese <- 4.554 ; label("Typical C-peptide short half-life, obese health status (min)") # Table 1 M0: mtso = 4.554 min, 3.758 %RSE
    thalfShort_diab <- 4.594 ; label("Typical C-peptide short half-life, diabetic health status (min)") # Table 1 M0: mtsd = 4.594 min, 3.727 %RSE

    # ---- Amplitude fraction F, per health-status stratum (Larizza 2018 Table 1, M0 column) ----
    famp_normal <- 0.764 ; label("Typical C-peptide amplitude fraction, normal health status (unitless)") # Table 1 M0: mFn = 0.764, 0.546 %RSE
    famp_obese <- 0.782 ; label("Typical C-peptide amplitude fraction, obese health status (unitless)") # Table 1 M0: mFo = 0.782, 0.629 %RSE
    famp_diab <- 0.780 ; label("Typical C-peptide amplitude fraction, diabetic health status (unitless)") # Table 1 M0: mFd = 0.780, 0.771 %RSE

    # ---- Long half-life tl: linear in AGE (Larizza 2018 Table 1, M0 column) ----
    thalfLong_intercept <- 27.797 ; label("Intercept of the C-peptide long half-life vs age regression (min)") # Table 1 M0: atl = 27.797 min, 4.802 %RSE
    e_age_thalfLong <- 0.177 ; label("Effect of AGE on the C-peptide long half-life (min/year)") # Table 1 M0: btl = 0.177 min/years, 22.728 %RSE

    # ---- Central volume V: linear in BSA, separate male and female strata ----
    # Both sexes carry their own intercept AND slope, so both are stratum-suffixed; neither
    # sex is an offset from the other (Larizza 2018 Mathematical Models, theta_V equation).
    vc_intercept_male <- 0.495 ; label("Intercept of the C-peptide central volume vs BSA regression, male (L)") # Table 1 M0: aVm = 0.495 L, 181.584 %RSE
    e_bsa_vc_male <- 1.982 ; label("Effect of BSA on the C-peptide central volume, male (L/m^2)") # Table 1 M0: bVm = 1.982 L/m^2, 22.630 %RSE
    vc_intercept_female <- 1.520 ; label("Intercept of the C-peptide central volume vs BSA regression, female (L)") # Table 1 M0: aVf = 1.520 L, 48.365 %RSE
    e_bsa_vc_female <- 1.432 ; label("Effect of BSA on the C-peptide central volume, female (L/m^2)") # Table 1 M0: bVf = 1.432 L/m^2, 28.790 %RSE

    # ---- Additive residual errors, one per macro constant ----
    # Larizza 2018 Table 1 footnote c: 'SDs of the additive residual errors' -- these are
    # standard deviations, not variances, and enter directly as add() magnitudes.
    # M0 estimates residual error and NO interindividual variability; M4 does the reverse.
    # With one observation per macro constant per subject the two are not simultaneously
    # identifiable, which is why the paper reports them in alternation.
    addSd_thalfShort <- 1.143 ; label("Additive residual SD on the short half-life (min)") # Table 1 M0: sigma_ADD_ts = 1.143, 5.394 %RSE
    addSd_famp <- 0.041 ; label("Additive residual SD on the amplitude fraction (unitless)") # Table 1 M0: sigma_ADD_F = 0.041, 5.625 %RSE
    addSd_thalfLong <- 5.778 ; label("Additive residual SD on the long half-life (min)") # Table 1 M0: sigma_ADD_tl = 5.778, 6.167 %RSE
    addSd_vc <- 0.846 ; label("Additive residual SD on the central volume (L)") # Table 1 M0: sigma_ADD_V = 0.846, 5.800 %RSE
  })

  model({
    # Larizza 2018 Mathematical Models, 'A population regression model to estimate CP kinetic
    # parameters'. Four linear regressions predict the macro constants of C-peptide kinetics:
    #
    #   theta_ts = mtsn | mtso | mtsd        by HSTATUS (normal | obese | diabetic)
    #   theta_F  = mFn  | mFo  | mFd         by HSTATUS
    #   theta_V  = aVm + bVm * BSA   (male)  |  aVf + bVf * BSA   (female)
    #   theta_tl = atl + btl * AGE
    #
    # The three health-status levels are mutually exclusive, so the normal-stratum indicator
    # is (1 - DIS_OBESE - DIS_DIAB). Writing the selector this way keeps the paper's three
    # typical values verbatim rather than re-parameterising two of them as offsets.
    isNormal <- 1 - DIS_OBESE - DIS_DIAB

    thalfShort <- thalfShort_normal * isNormal + thalfShort_obese * DIS_OBESE + thalfShort_diab * DIS_DIAB
    famp <- famp_normal * isNormal + famp_obese * DIS_OBESE + famp_diab * DIS_DIAB
    thalfLong <- thalfLong_intercept + e_age_thalfLong * AGE
    vc <-
      (vc_intercept_male + e_bsa_vc_male * BSA) * (1 - SEXF) +
      (vc_intercept_female + e_bsa_vc_female * BSA) * SEXF

    # The micro rate constants of the linear two-compartment C-peptide model follow from the
    # macro constants by the algebraic relations in Larizza 2018 Mathematical Models:
    #
    #   k12 = ln(2) * (F / tl + (1 - F) / ts)
    #   k01 = (ln(2) / ts) * (ln(2) / tl) * (1 / k12)
    #   k21 =  ln(2) / ts  +  ln(2) / tl  - k12 - k01
    #
    # They are deliberately NOT computed inside model(): this model is fitted to the macro
    # constants, the conversion is a pure post-processing step, and defining a k12 / k21 / vc
    # set inside an ODE-free model() risks rxode2 interpreting it as a solved linear system.
    # The vignette applies the conversion to the simulated macro constants and feeds the
    # result to Larizza_2018_insulinMinimalModel; see that model's ini() for the values.

    # Four observed macro constants, each with its own additive residual error.
    thalfShort ~ add(addSd_thalfShort)
    famp ~ add(addSd_famp)
    thalfLong ~ add(addSd_thalfLong)
    vc ~ add(addSd_vc)
  })
}
