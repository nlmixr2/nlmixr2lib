Glazar_2025_glioma_tgi <- function() {
  description <- paste(
    "Claret tumor growth inhibition (TGI) model of T2/FLAIR non-enhancing",
    "tumor volume coupled with a deterministic threshold model of",
    "progression-free survival in recurrent high-grade glioma treated with",
    "hypofractionated stereotactic radiotherapy (HFSRT) plus bevacizumab plus",
    "an immune checkpoint inhibitor (Glazar 2025). Tumor volume grows",
    "exponentially at rate p and, from treatment initiation (time 0), is",
    "killed at rate kkill * exp(-res * t), where res is the rate of evolution",
    "of resistance. The individual time to progression ttp is the time the",
    "tumor reaches kappa = 1.4 times its analytic nadir (or its baseline when",
    "it never shrinks), computed in closed form with the Lambert W function.",
    "No pharmacokinetics and no dose input: the three treatment modalities are",
    "not separated. Log-normal IIV on all four parameters with a diagonal",
    "omega matrix; no covariates.",
    sep = " "
  )
  reference <- paste(
    "Glazar DJ, Sahebjam S, Yu HHM, Chen DT, Bhandari M, Enderling H.",
    "A sample size analysis of a mathematical model of longitudinal tumor",
    "volume and progression-free survival for Bayesian individual dynamic",
    "predictions in recurrent high-grade glioma.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(3):495.",
    "doi:10.1002/psp4.13290.",
    "Parameter values are from Table 1 and the Stan code of Data S3; the",
    "time-to-progression solution and the Lambert W approximation are from",
    "Data S1 sections 3-5.",
    sep = " "
  )
  vignette <- "Glazar_2025_glioma_tgi"

  units <- list(
    time = "day",
    dosing = "none (no drug input; the treatment effect starts at time 0)",
    concentration = "cm^3 (the observed variable is T2/FLAIR tumor volume, not a drug concentration)"
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "T2/FLAIR non-enhancing tumor volume",
      units = "cm^3",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 52L,
    n_studies = 2L,
    disease_state = "recurrent high-grade glioma",
    dose_range = "hypofractionated stereotactic radiotherapy combined with bevacizumab and an immune checkpoint inhibitor (pembrolizumab in NCT02313272; nivolumab plus ipilimumab in NCT02829931); the model does not resolve individual modalities or doses",
    regions = "United States (two single-centre phase I trials)",
    notes = "Population parameters were estimated from 52 patients with recurrent high-grade glioma pooled from two phase I trials (NCT02313272 and NCT02829931; Glazar 2025 Methods, Data generation). Glazar 2025 does not tabulate demographics for this cohort. The paper itself applies the model to 2000 simulated patients imaged every 6 weeks with administrative censoring at day 714."
  )

  ini({
    # Glazar 2025 Table 1 reports the typical value mu and the IIV omega of
    # each parameter, and Data S3 (the Stan estimation code) hard-codes the
    # same values as TS0_pop, kge_pop, kkill_pop and lambda_pop. Equation 2:
    # ln(psi_i) = ln(mu) + eta_i, eta_i ~ N(0, omega), with omega the SD of
    # the log-transformed parameter.
    lrbase <- log(101.12)
    label("Tumor volume at treatment initiation TS0 (cm^3)") # Table 1 TS0 typical value 101.12 cm3; Data S3 TS0_pop = 101.12

    lp <- log(0.0020)
    label("Exponential tumor growth rate kge (1/day)") # Table 1 kge typical value 0.0020 day-1; Data S3 kge_pop = 0.0020

    lkkill <- log(0.0245)
    label("Maximum tumor kill rate kkill at treatment initiation (1/day)") # Table 1 kkill typical value 0.0245 day-1; Data S3 kkill_pop = 0.0245

    lres <- log(0.0474)
    label("Rate of evolution of resistance lambda (1/day)") # Table 1 lambda typical value 0.0474 day-1; Data S3 lambda_pop = 0.0474

    # Relative progression threshold, derived in Data S1 section 2 by
    # extrapolating the RANO 25% increase in the product of perpendicular
    # diameters to a spherical volume: (sqrt(5)/2)^3 = 5*sqrt(5)/8 = 1.398,
    # rounded to 1.4 in the text and the Stan code.
    kappa <- fixed(1.4)
    label("Progression threshold: tumor volume relative to nadir at progression (ratio)") # Methods, Progression-free survival: kappa = 1.4; Data S2 transformed data kappa = 1.4

    # IIV: Table 1 omega is the SD of the log-parameter (Equation 2); the
    # variances below are omega^2. Diagonal omega (Methods: 'the covariance
    # matrix is assumed to be diagonal').
    etalrbase ~ 0.4225 # Table 1 TS0 omega 0.65 -> 0.65^2
    etalp ~ 0.3136 # Table 1 kge omega 0.56 -> 0.56^2
    etalkkill ~ 0.04 # Table 1 kkill omega 0.20 -> 0.20^2
    etalres ~ 0.5329 # Table 1 lambda omega 0.73 -> 0.73^2

    # Residual error, Equation 3: ln(y) = ln(TS) + eps, eps ~ N(0, a). The
    # paper sets a = 0.01 (printed with unit cm3, but it acts on the log
    # scale); Data S3 uses it as the SD in y ~ lognormal(log(TS), sigma).
    expSd <- fixed(0.01)
    label("Log-scale residual SD on tumor volume (unitless)") # Equation 3 a = 0.01; Data S3 sigma = 0.01
  })
  model({
    rbase <- exp(lrbase + etalrbase)
    p <- exp(lp + etalp)
    kkill <- exp(lkkill + etalkkill)
    res <- exp(lres + etalres)

    # Equation 1 (Claret TGI). Time 0 is treatment initiation, so the
    # treatment indicator [t >= 0] of Equation 1 is always 1 here. The
    # pre-treatment course is TS0 * exp(p * t) for t < 0 (Data S2 TGI()).
    tumor_size(0) <- rbase
    d/dt(tumor_size) <- p * tumor_size - kkill * exp(-res * t) * tumor_size

    # Time to progression, Equation 4 / Data S1 Equation S8: the time at
    # which tumor_size reaches kappa times its nadir (kge < kkill) or kappa
    # times TS0 (kge >= kkill, no initial shrinkage). The principal-branch
    # Lambert W is evaluated with the authors' closed-form approximation
    # (Data S1 section 4, epsilon = 1e-12), transcribed from the TTE()
    # function of Data S2. Data S2 additionally rounds TTP up to whole days
    # with ceil(); ttp here is the unrounded time.
    if (p < kkill) {
      w_arg <- -exp(-(res / p) * log(kappa) - 1)
    } else {
      w_arg <- -kkill / p * exp((-res * log(kappa) - kkill) / p)
    }
    w_r <- exp(-1) - 1e-12
    w_s <- w_arg + w_r + 1e-12
    w_z0 <- sqrt(2 * w_r * w_s)
    w_z <- w_z0 * sqrt(1 + w_z0 / (3 * w_r))
    lambert_w <- w_arg * log(1 + exp(1) * w_z) / (w_s + w_z)
    if (p < kkill) {
      ttp <- (lambert_w + log(kkill / p) + 1) / res + log(kappa) / p
    } else {
      ttp <- lambert_w / res + log(kappa) / p + kkill / (res * p)
    }

    tumor_size ~ lnorm(expSd)
  })
}
