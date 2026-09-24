Larizza_2018_cpeptideKinetics_M4 <- function() {
  description <- paste(
    "Population regression model (Magni 2000 'M4') predicting the four macro constants of",
    "C-peptide two-compartment kinetics -- short half-life, amplitude fraction, long half-life",
    "and central distribution volume -- from health status (normal / obese / diabetic), sex,",
    "body surface area and age in 207 adults. This is the Bayesian variant, estimated in",
    "WinBUGS with a full 4x4 interindividual covariance matrix and no separate residual error;",
    "the maximum-likelihood variant with independent additive residual errors is",
    "Larizza_2018_cpeptideKinetics_M0. The macro constants convert algebraically to the micro",
    "rate constants k01, k12 and k21 that Larizza_2018_insulinMinimalModel requires.",
    sep = " "
  )
  reference <- paste(
    "Larizza C, Borella E, Pasotti L, Tartaglione P, Smith M, Moodie S, Magni P. (2018).",
    "Complex Bayesian Modeling Workflows Encoding and Execution Made Easy With a Novel",
    "WinBUGS Plugin of the Drug Disease Model Resources Interoperability Framework.",
    "CPT Pharmacometrics Syst Pharmacol 7(5):298-308. doi:10.1002/psp4.12285.",
    "Underlying model: Magni P, Bellazzi R, Sparacino G, Cobelli C.",
    "Bayesian identification of a population compartmental model of C-peptide kinetics.",
    "Ann Biomed Eng 2000;28:812-823.",
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
      "Same 207-subject dataset as Larizza_2018_cpeptideKinetics_M0 (Larizza 2018 Datasets",
      "section, citing Van Cauter 1992 and Magni 2000). M4 was estimated in WinBUGS with one",
      "chain, 1,000 burn-in iterations, 100,000 updates and a thin of 10, repeated three times",
      "using the last chain values as the next run's initial values, giving 300,000 chain",
      "samples (Larizza 2018 Figure 2 caption).",
      sep = " "
    )
  )

  ini({
    # ---- Short half-life ts, per health-status stratum (Larizza 2018 Table 1, M4 column) ----
    # Natural scale, NOT log-transformed: the source regression is additive-normal on the
    # original scale (phi_i = theta_i + gamma_i, gamma_i ~ N(0, Xi)), so a log
    # parameterisation would misstate the structure the random effects sit on.
    thalfShort_normal <- 4.991 ; label("Typical C-peptide short half-life, normal health status (min)") # Table 1 M4: mtsn = 4.991 min, 1.942 %RSE
    thalfShort_obese <- 4.496 ; label("Typical C-peptide short half-life, obese health status (min)") # Table 1 M4: mtso = 4.496 min, 2.921 %RSE
    thalfShort_diab <- 4.693 ; label("Typical C-peptide short half-life, diabetic health status (min)") # Table 1 M4: mtsd = 4.693 min, 2.957 %RSE

    # ---- Amplitude fraction F, per health-status stratum (Larizza 2018 Table 1, M4 column) ----
    famp_normal <- 0.766 ; label("Typical C-peptide amplitude fraction, normal health status (unitless)") # Table 1 M4: mFn = 0.766, 0.562 %RSE
    famp_obese <- 0.781 ; label("Typical C-peptide amplitude fraction, obese health status (unitless)") # Table 1 M4: mFo = 0.781, 0.784 %RSE
    famp_diab <- 0.778 ; label("Typical C-peptide amplitude fraction, diabetic health status (unitless)") # Table 1 M4: mFd = 0.778, 0.858 %RSE

    # ---- Long half-life tl: linear in AGE (Larizza 2018 Table 1, M4 column) ----
    thalfLong_intercept <- 26.705 ; label("Intercept of the C-peptide long half-life vs age regression (min)") # Table 1 M4: atl = 26.705 min, 3.854 %RSE
    e_age_thalfLong <- 0.209 ; label("Effect of AGE on the C-peptide long half-life (min/year)") # Table 1 M4: btl = 0.209 min/years, 13.197 %RSE

    # ---- Central volume V: linear in BSA, separate male and female strata ----
    vc_intercept_male <- 0.344 ; label("Intercept of the C-peptide central volume vs BSA regression, male (L)") # Table 1 M4: aVm = 0.344 L, 131.067 %RSE
    e_bsa_vc_male <- 2.061 ; label("Effect of BSA on the C-peptide central volume, male (L/m^2)") # Table 1 M4: bVm = 2.061 L/m^2, 10.730 %RSE
    vc_intercept_female <- 0.795 ; label("Intercept of the C-peptide central volume vs BSA regression, female (L)") # Table 1 M4: aVf = 0.795 L, 59.352 %RSE
    e_bsa_vc_female <- 1.819 ; label("Effect of BSA on the C-peptide central volume, female (L/m^2)") # Table 1 M4: bVf = 1.819 L/m^2, 13.898 %RSE

    # ---- Full 4x4 interindividual covariance matrix Xi (Larizza 2018 Table 1, M4 column) ----
    # Table 1 footnote b: 'Elements of the full matrix Xi'. These are variances and covariances
    # on the ORIGINAL (natural) scale, applied additively: phi_i = theta_i + gamma_i with
    # gamma_i ~ N(0, Xi). Ordering below is the lower triangle of Xi read row-major over
    # (ts, F, tl, V):
    #   var(ts)        = 1.295   (9.781 %RSE)
    #   cov(ts, F)     = 0.006   (60.110 %RSE)   var(F)  = 0.002   (9.822 %RSE)
    #   cov(ts, tl)    = 3.250   (15.488 %RSE)   cov(F, tl) = 0.071 (26.966 %RSE)   var(tl) = 33.044 (9.847 %RSE)
    #   cov(ts, V)     = 0.596   (13.214 %RSE)   cov(F, V) = -0.006 (49.064 %RSE)   cov(tl, V) = 1.915 (19.095 %RSE)   var(V) = 0.713 (9.796 %RSE)
    # Implied correlations are 0.118 (ts,F), 0.497 (ts,tl), 0.620 (ts,V), 0.276 (F,tl),
    # -0.159 (F,V) and 0.394 (tl,V); the matrix is positive definite (smallest eigenvalue
    # 0.0016). Table 1 prints the unit of xi_tl as '-' where min^2 is meant, and the unit of
    # xi_F as 'min^2' where the amplitude fraction is unitless; both are label slips and do
    # not affect the numbers.
    etathalfShort + etafamp + etathalfLong + etavc ~ c(
      1.295,
      0.006, 0.002,
      3.250, 0.071, 33.044,
      0.596, -0.006, 1.915, 0.713
    )

    # ---- Residual error ----
    # M4 reports NO residual-error terms (Larizza 2018 Table 1: the sigma_ADD rows are blank in
    # the M4 column). With a single observation of each macro constant per subject, Xi and a
    # residual error are not simultaneously identifiable, and M4 places all of the stochastic
    # structure on Xi. nlmixr2 requires an error model per endpoint, so each is pinned at zero
    # rather than invented; see the vignette Errata.
    addSd_thalfShort <- fixed(0) ; label("Additive residual SD on the short half-life (min; not estimated in M4)") # Table 1: M4 column blank for sigma_ADD_ts
    addSd_famp <- fixed(0) ; label("Additive residual SD on the amplitude fraction (unitless; not estimated in M4)") # Table 1: M4 column blank for sigma_ADD_F
    addSd_thalfLong <- fixed(0) ; label("Additive residual SD on the long half-life (min; not estimated in M4)") # Table 1: M4 column blank for sigma_ADD_tl
    addSd_vc <- fixed(0) ; label("Additive residual SD on the central volume (L; not estimated in M4)") # Table 1: M4 column blank for sigma_ADD_V
  })

  model({
    # Larizza 2018 Mathematical Models, 'A population regression model to estimate CP kinetic
    # parameters'. Identical structural regressions to M0; M4 differs only in carrying the
    # full interindividual covariance matrix Xi and a Bayesian prior specification:
    #
    #   theta ~ N(theta_0, Sigma_0^-1) with theta_0 = [5 5 5 1 1 1 30 1 1 1 1 1] and Sigma_0
    #     diagonal with the squared elements of theta_0
    #   Xi^-1 ~ Wishart(q, R) with q = 10 and R = q^-1 * (0.01 * diag([5 1 30 4]))^-1
    #
    # The priors govern estimation only and are not part of the simulation model, so they are
    # recorded here rather than in ini().
    isNormal <- 1 - DIS_OBESE - DIS_DIAB

    thalfShort <-
      thalfShort_normal * isNormal + thalfShort_obese * DIS_OBESE + thalfShort_diab * DIS_DIAB +
      etathalfShort
    famp <- famp_normal * isNormal + famp_obese * DIS_OBESE + famp_diab * DIS_DIAB + etafamp
    thalfLong <- thalfLong_intercept + e_age_thalfLong * AGE + etathalfLong
    vc <-
      (vc_intercept_male + e_bsa_vc_male * BSA) * (1 - SEXF) +
      (vc_intercept_female + e_bsa_vc_female * BSA) * SEXF +
      etavc

    # The micro rate constants k01 / k12 / k21 follow from the macro constants by the algebraic
    # relations in Larizza 2018 Mathematical Models (see Larizza_2018_cpeptideKinetics_M0 for
    # the equations). They are applied as a post-processing step in the vignette rather than
    # inside model(), so that the sampled joint distribution of (ts, F, tl, V) propagates
    # through to a joint distribution of (k01, k12, k21) -- which is exactly what Larizza 2018
    # approach 2 uses as the empirical prior for the insulin minimal model.

    # Four observed macro constants. The additive residuals are pinned at zero (see ini()), so
    # the sampled observations equal phi_i = theta_i + gamma_i, the quantity M4 models.
    thalfShort ~ add(addSd_thalfShort)
    famp ~ add(addSd_famp)
    thalfLong ~ add(addSd_thalfLong)
    vc ~ add(addSd_vc)
  })
}
