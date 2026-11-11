Rohleff_2025_neonatalWeight <- function() {
  description <- "Endogenous turnover model of body-weight change over the first 7 days of life in 2425 neonates (Rohleff 2025, a smoothed version of the Wilbaux 2016/2019 neonatal weight model), fitted without covariates by a variational-autoencoder NLME estimator. Weight starts at the birth weight W0 and follows dW/dt = kprod(t) - kel(t) * W, where weight-gain production rises sigmoidally around a delay Tlag (kprod = kin / (1 + exp(-2 (t - Tlag)))) and the fractional loss rate declines hyperbolically (kel = kout * (1 - t / (T50 + t))). Log-normal IIV on all five parameters; additive residual error in grams. Time is postnatal age in days."
  reference <- paste(
    "Rohleff J, Bachmann F, Nahum U, Bram D, Steffens B, Pfister M, Koch G,",
    "Schropp J. Redefining Parameter Estimation and Covariate Selection via",
    "Variational Autoencoders: One Run Is All You Need.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):2232-2243.",
    "doi:10.1002/psp4.70129.",
    "Parameter values: Supporting Information S1, Table S3 (VAE column).",
    sep = " "
  )
  vignette <- "Rohleff_2025_theophylline_neonatalWeight"
  units <- list(
    time = "day",
    dosing = "n/a (no dosing; endogenous body-weight model)",
    concentration = "g (body weight)"
  )

  covariateData <- list()

  # The paper's covariate-selection run (Section 3.2, Figures 4 and 5) retains
  # covariate effects, but only the covariate-FREE population fit is reported
  # with its fixed effects, omegas and residual error (Table S3). The fixed
  # effects and omegas of the covariate model are not printed anywhere, so the
  # covariate model cannot be reconstructed; the selected covariates are
  # documented here only.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Neonate sex",
      units = "(binary)",
      type = "binary",
      reference_category = "Not stated",
      notes = "Source column 'sex' (Section 3.2). The paper does not state which sex is coded 1. Selected by the VAE on W0 and kin (Figure 4); final coefficients appear only as convergence traces in Figure 5.",
      source_name = "sex"
    ),
    DelM = list(
      description = "Mode of delivery (categorical; paper column name, not a registered canonical)",
      units = "(binary)",
      type = "binary",
      reference_category = "Not stated",
      notes = "Source column 'DelM' (Section 3.2). Coding not stated in the paper. Selected by the VAE on W0 and kin (Figure 4).",
      source_name = "DelM"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column 'GA' (Section 3.2), entered as log(GA / mean GA). Units are not stated in the paper; weeks is assumed from the clinical convention. Selected by the VAE on W0, kin and T50 (Figure 4).",
      source_name = "GA"
    ),
    Mage = list(
      description = "Maternal age (paper column name, not a registered canonical)",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column 'Mage' (Section 3.2), entered as log(Mage / mean Mage). Selected by the VAE on kin, Tlag and kout (Figure 4).",
      source_name = "Mage"
    ),
    Para2 = list(
      description = "Parity indicator (paper column name 'Para2', not a registered canonical)",
      units = "(binary)",
      type = "binary",
      reference_category = "Not stated",
      notes = "Source column 'Para2' (Section 3.2). The paper does not define the categories. Selected by the VAE on W0, kin and Tlag (Figure 4).",
      source_name = "Para2"
    )
  )

  compartmentData <- list(
    bw = list(analyte = "body weight", units = "g", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2425L,
    n_studies = 1L,
    age_range = "Birth to 7 days postnatal age",
    disease_state = "Neonates followed over the first week of life (the cohort of Wilbaux et al., J Pediatr 2016 and Clin Nutr 2019, the paper's references 19 and 20, both of which describe term neonates).",
    dose_range = "n/a (no drug)",
    notes = paste(
      "Section 3.2 'Data': repeated body weights from 2425 neonates over the",
      "first 7 days of life, with sex, mode of delivery, gestational age,",
      "maternal age and parity as candidate covariates. The paper reports no",
      "demographic table; the original data are not shared (the code archive",
      "ships a simulated stand-in dataset instead)."
    )
  )

  ini({
    lrbase <- log(3383.2); label("Birth weight W0, the weight at t = 0 (g)") # Table S3, VAE column, row 'W0,pop' = 3383.2
    lkin <- log(106.8); label("Maximal zero-order weight-gain rate kin (g/day)") # Table S3, VAE column, row 'kin,pop' = 106.8
    lkin_t50 <- log(1.56); label("Delay Tlag: postnatal age at which weight gain reaches half of kin (day)") # Table S3, VAE column, row 'Tlag,pop' = 1.56
    lkout <- log(0.072); label("Initial fractional weight-loss rate kout (1/day)") # Table S3, VAE column, row 'kout,pop' = 0.072
    lkout_t50 <- log(1.18); label("T50: postnatal age at which the fractional loss rate has halved (day)") # Table S3, VAE column, row 'T50,pop' = 1.18

    # Table S3 reports omega as the standard deviation of the log-normal random
    # effect, so each variance below is omega^2.
    etalrbase ~ 0.0196 # Table S3, VAE column, row 'omega_W0' = 0.14 (SD); 0.14^2
    etalkin ~ 0.0529 # Table S3, VAE column, row 'omega_kin' = 0.23 (SD); 0.23^2
    etalkin_t50 ~ 0.0144 # Table S3, VAE column, row 'omega_Tlag' = 0.12 (SD); 0.12^2
    etalkout ~ 0.01 # Table S3, VAE column, row 'omega_kout' = 0.10 (SD); 0.10^2
    etalkout_t50 ~ 0.0576 # Table S3, VAE column, row 'omega_T50' = 0.24 (SD); 0.24^2

    addSd <- 33.0; label("Additive residual error on body weight (g)") # Table S3, VAE column, row 'Error Model a' = 33.0
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    kin <- exp(lkin + etalkin)
    kin_t50 <- exp(lkin_t50 + etalkin_t50)
    kout <- exp(lkout + etalkout)
    kout_t50 <- exp(lkout_t50 + etalkout_t50)

    # Equation 15. The factor 2 (1/day) is the fixed steepness of the smoothed
    # lag; Rohleff 2025 replaced the step lag of the original model with this
    # sigmoid so that the model is differentiable.
    kprod <- kin / (1 + exp(-2 * (t - kin_t50)))
    kelim <- kout * (1 - t / (kout_t50 + t))

    bw(0) <- rbase
    d/dt(bw) <- kprod - kelim * bw

    bw ~ add(addSd)
  })
}
