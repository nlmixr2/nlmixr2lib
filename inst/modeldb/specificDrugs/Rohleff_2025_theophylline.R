Rohleff_2025_theophylline <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for single-dose oral theophylline in 12 subjects (the classic Boeckmann / Pinheiro-Bates theophylline dataset), fitted by Rohleff 2025 with a variational-autoencoder NLME estimator that also selected the covariates. Log-normal IIV on ka, ke and V; body weight (normalised to the cohort mean) on ka and V; additive residual error."
  reference <- paste(
    "Rohleff J, Bachmann F, Nahum U, Bram D, Steffens B, Pfister M, Koch G,",
    "Schropp J. Redefining Parameter Estimation and Covariate Selection via",
    "Variational Autoencoders: One Run Is All You Need.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(8):2232-2243.",
    "doi:10.1002/psp4.70129.",
    sep = " "
  )
  vignette <- "Rohleff_2025_theophylline_neonatalWeight"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as log(WT / wpop) on log(ka) and log(V) (Section 3.1,",
        "'Model Description'), i.e. a power model (WT / wpop)^beta. wpop is",
        "'the mean weight over the N individuals'; the paper does not print",
        "it, so 69.583 kg is computed from the 12 subject weights in the",
        "authors' deposited dataset (Supporting Information code archive,",
        "Data/theophylline_data.tab), which are those of the R datasets::Theoph",
        "data (79.6, 72.4, 70.5, 72.7, 54.6, 80.0, 64.6, 70.5, 86.4, 58.2,",
        "65.0, 60.5 kg; sum 835.0 kg)."
      ),
      source_name = "w"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "The paper codes sex as c_sex = 1 for male and 0 for female (Section",
        "3.1), i.e. SEXF = 1 - c_sex. Sex was offered on ka, ke and V, and the",
        "VAE (and COSSAC, SAMBA and SCM) set every sex effect to zero",
        "(Table 1, rows beta_ka,sex / beta_ke,sex / beta_V,sex = '-'; Figure 2)."
      ),
      source_name = "sex"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "theophylline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "theophylline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    weight_range = "54.6-86.4 kg (mean 69.6 kg; authors' deposited dataset)",
    sex_female_pct = 41.7,
    disease_state = "Not stated in the paper ('N = 12 patients'); the dataset is the classic theophylline study distributed with NONMEM and the R nlme package (Pinheiro and Bates 2000, the paper's reference 18).",
    dose_range = "Single oral theophylline loading dose at t = 0; 267.8-320.7 mg total (3.10-5.86 mg/kg) in the deposited dataset.",
    n_observations = 120L,
    notes = paste(
      "Section 3.1 'Data': plasma theophylline concentrations from 12",
      "subjects after a single oral dose at t = 0, with baseline weight and",
      "sex recorded. 10 post-dose concentrations per subject (120 in all) in",
      "the deposited dataset; the time-zero pre-dose samples of",
      "datasets::Theoph are not in it. Sex count (5 female, 7 male) is from",
      "the deposited dataset using the paper's c_sex = 1 for male coding."
    )
  )

  ini({
    lka <- log(1.63); label("Absorption rate constant ka at the mean weight (1/h)") # Table 1, VAE column, row 'ka,pop' = 1.63
    lkel <- log(0.086); label("Elimination rate constant ke at the mean weight (1/h)") # Table 1, VAE column, row 'ke,pop' = 0.086
    lvc <- log(31.97); label("Apparent volume of distribution V/F at the mean weight (L)") # Table 1, VAE column, row 'Vpop' = 31.97

    e_wt_ka <- 2.55; label("Exponent of (WT / mean WT) on ka (unitless)") # Table 1, VAE column, row 'beta_ka,w' = 2.55
    e_wt_vc <- 0.57; label("Exponent of (WT / mean WT) on V/F (unitless)") # Table 1, VAE column, row 'beta_V,w' = 0.57

    # Table 1 reports omega as the standard deviation of the log-normal random
    # effect, so each variance below is omega^2.
    etalka ~ 0.2809 # Table 1, VAE column, row 'omega_ka' = 0.53 (SD); 0.53^2
    etalkel ~ 0.0225 # Table 1, VAE column, row 'omega_ke' = 0.15 (SD); 0.15^2
    etalvc ~ 0.0169 # Table 1, VAE column, row 'omega_V' = 0.13 (SD); 0.13^2

    addSd <- 0.71; label("Additive residual error (mg/L)") # Table 1, VAE column, row 'Error model a' = 0.71 (constant error model, Equation 2)
  })

  model({
    # Covariate: log(w_i / w_pop) with w_pop the mean weight of the 12 subjects
    # (Section 3.1); 69.583 kg is computed from the deposited dataset.
    wt_ref <- 69.583
    ka <- exp(lka + etalka) * (WT / wt_ref)^e_wt_ka
    kel <- exp(lkel + etalkel)
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc

    # Equation 14: d/dt Abs = -ka Abs, d/dt A = ka Abs - ke A, Abs(0) = D.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
