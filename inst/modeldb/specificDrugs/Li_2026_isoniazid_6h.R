Li_2026_isoniazid_6h <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for oral isoniazid given as daily monotherapy (6H regimen, 300 mg once daily) for latent tuberculosis infection preventive therapy in Chinese adults (Li 2026). NAT2 acetylator phenotype (intermediate reference, with log-linear shifts for slow and fast acetylators) and a power function of body weight normalised to the regimen median act on apparent oral clearance; body weight also acts on the apparent volume. Log-normal between-subject variability on CL/F and V/F and a log-normal (exponential) residual error. One of three independently fitted regimen-specific models from the same paper (see Li_2026_isoniazid_3hr and Li_2026_isoniazid_3h2p2)."
  reference <- paste(
    "Li Z, Xiao X, Xu C, Liu Y, Gu L, Li X, Shen X, Hu Y. (2026).",
    "Regimen-Specific Population Pharmacokinetics of Isoniazid with and",
    "Without Rifamycin: A Bayesian Modeling Analysis of 6H, 3HR, and 3H2P2",
    "Regimens. Pharmaceutics 18(8):937. doi:10.3390/pharmaceutics18080937.",
    sep = " "
  )
  vignette <- "Li_2026_isoniazid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "isoniazid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isoniazid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power function on CL/F and V/F normalised to 66.4 kg, the median body weight of the 6H cohort (Table 1). The authors' Stan code (Supplementary Code, prepare_regimen_data) centres weight as log(Weight / median(Weight)) within each regimen's own dataset, so each regimen-specific model has its own reference weight.",
      source_name = "Weight"
    ),
    NAT2_SLOW = list(
      description = "NAT2 slow-acetylator phenotype indicator (1 = slow acetylator, two reduced-function alleles; 0 = intermediate or fast). Paired with NAT2_RAPID to encode the three-level NAT2 phenotype.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the joint state NAT2_SLOW = 0 and NAT2_RAPID = 0 denotes an intermediate acetylator, the reference phenotype of the source model.",
      notes = "Methods 2.3: NAT2 genotype by PCR; phenotype assigned by the number of reduced-function alleles (two = slow, one = intermediate, none = rapid). 6H cohort: 20 slow (32.3%), 18 intermediate (29.0%), 24 fast (38.7%) (Table 1).",
      source_name = "NAT2"
    ),
    NAT2_RAPID = list(
      description = "NAT2 rapid (fast) acetylator phenotype indicator (1 = rapid acetylator, no reduced-function allele; 0 = intermediate or slow). Paired with NAT2_SLOW to encode the three-level NAT2 phenotype.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the joint state NAT2_SLOW = 0 and NAT2_RAPID = 0 denotes an intermediate acetylator, the reference phenotype of the source model.",
      notes = "Same phenotyping source as NAT2_SLOW. The paper labels this level 'Fast' (Table 1) and 'Fast/Rapid' (Table 3).",
      source_name = "NAT2"
    )
  )

  # Screened in the stepwise covariate analysis (Methods 2.3) but not retained
  # in the final model; no point estimates are reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods 2.3) and not retained. 6H cohort median 38 years (IQR 32.3-45.3) (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 2.3) and not retained. 6H cohort 31 of 62 (50.0%) male (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 62L,
    n_studies = 1L,
    n_observations = 372L,
    age_median = "38 years (IQR 32.3-45.3)",
    weight_median = "66.4 kg (IQR 59.0-72.0)",
    sex_female_pct = 50.0,
    race_ethnicity = "Chinese (Sichuan Province); ancestry not otherwise reported",
    disease_state = "Latent tuberculosis infection receiving tuberculosis preventive therapy",
    dose_range = "Isoniazid 300 mg orally once daily for 6 months (6H), without a rifamycin",
    regions = "Sichuan Province, China",
    nat2_phenotype = "20 slow (32.3%), 18 intermediate (29.0%), 24 fast (38.7%)",
    notes = "Prospective observational study, July 2024 to June 2025 (Methods 2.1). Plasma sampled after two weeks of treatment at predose and 1, 2, 4, 6 and 8 h post-dose; 94 of 372 observations (25.3%) below the LLOQ, handled by the M3 likelihood. Baseline characteristics are Table 1. Bayesian estimation in Stan (CmdStanR 0.9.0); each regimen was fitted independently with no parameter shared across regimens (Supplementary Code)."
  )

  ini({
    # Structural parameters. Posterior medians for the 6H regimen from Table 2
    # (also Results 3.2.1). They are the typical values for an intermediate
    # acetylator at the regimen median weight of 66.4 kg.
    lka <- log(2.88); label("First-order absorption rate constant ka (1/h)") # Table 2: Ka 6H = 2.88 1/h (94% HDI 1.17-7.36; RSE 62.09%)
    lcl <- log(21.83); label("Apparent oral clearance CL/F, NAT2 intermediate acetylator at 66.4 kg (L/h)") # Table 2: CL/F 6H = 21.83 L/h (94% HDI 18.53-25.23; RSE 8.18%)
    lvc <- log(51.03); label("Apparent volume of distribution V/F at 66.4 kg (L)") # Table 2: V/F 6H = 51.03 L (94% HDI 42.50-60.62; RSE 9.32%)

    # Body-weight power exponents. Table 2 reports one value per parameter
    # 'Across regimens' (Results 3.2); the paper does not print the
    # regimen-specific exponents.
    e_wt_cl <- 0.728; label("Power exponent of body weight on CL/F (unitless)") # Table 2: Weight on CL/F = 0.728 (94% HDI 0.412-1.070; RSE 23.90%)
    e_wt_vc <- 0.999; label("Power exponent of body weight on V/F (unitless)") # Table 2: Weight on V/F = 0.999 (94% HDI 0.606-1.436; RSE 21.95%)

    # NAT2 phenotype effects on log CL/F (Stan beta_nat_CL; effect ratio =
    # exp(beta)). Table 3 and Results 3.2.4 report the effect ratios versus
    # the intermediate reference: slow 0.646, fast 1.484.
    e_nat2_slow_cl <- log(0.646); label("Log effect of NAT2 slow vs intermediate acetylator on CL/F (unitless)") # Table 3 / Results 3.2.4: slow effect ratio 0.646 (94% HDI 0.573-0.719), -35.4%
    e_nat2_rapid_cl <- log(1.484); label("Log effect of NAT2 fast vs intermediate acetylator on CL/F (unitless)") # Table 3 / Results 3.2.4: fast effect ratio 1.484 (94% HDI 1.327-1.642), +48.4%

    # IIV. Table 2 reports %CV computed in the Supplementary Code as
    # sqrt(exp(omega^2) - 1) * 100, so omega^2 = log(1 + CV^2).
    etalcl ~ 0.05157 # Table 2: IIV CL/F = 23.00 %CV (94% HDI 15.84-29.51); log(1 + 0.23^2)
    etalvc ~ 0.006348 # Table 2: IIV V/F = 7.98 %CV (94% HDI 0.001-22.79); log(1 + 0.0798^2)

    # Residual error: log(conc) ~ normal(log(Cpred), sigma) in the Stan model.
    expSd <- 0.895; label("Residual error SD on the log-concentration scale (unitless)") # Table 2: Proportional error (log-scale SD) = 0.895 (94% HDI 0.850-0.941)
  })

  model({
    # Covariate model (Supplementary Code Stan block): log CL/F = log_CL +
    # beta_nat_CL + beta_wt_CL * log(WT / median WT) + eta; log V/F = log_V +
    # beta_wt_V * log(WT / median WT) + eta. Ka has no IIV.
    ka <- exp(lka)
    cl <- exp(lcl + e_nat2_slow_cl * NAT2_SLOW + e_nat2_rapid_cl * NAT2_RAPID + etalcl) * (WT / 66.4)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 66.4)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose mg / volume L = mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
