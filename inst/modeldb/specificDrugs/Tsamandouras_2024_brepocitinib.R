Tsamandouras_2024_brepocitinib <- function() {
  description <- paste(
    "Exposure-response model for the Psoriasis Area and Severity Index (PASI)",
    "under oral brepocitinib (PF-06700841, a dual TYK2/JAK1 inhibitor) in adults",
    "with moderate-to-severe plaque psoriasis (Tsamandouras 2024). The PASI",
    "score divided by 72 follows a zero-inflated beta distribution: a point",
    "mass at PASI = 0 whose logit falls linearly with the beta mean, and a",
    "beta density with mean mu and precision phi otherwise. The expected",
    "score is the baseline (logit-normal between 12 and 72) minus an",
    "exponential placebo effect that removes at most a fraction Pmax of",
    "baseline, minus a drug effect on the remaining fraction driven by a type",
    "I indirect-response latent variable whose production the individual",
    "average steady-state brepocitinib concentration (CAV) inhibits."
  )
  reference <- paste(
    "Tsamandouras N, Qiu R, Hughes JH, Sweeney K, Prybylski JP, Banfield C,",
    "Nicholas T. Employing zero-inflated beta distribution in an",
    "exposure-response analysis of TYK2/JAK1 inhibitor brepocitinib in",
    "patients with plaque psoriasis. J Pharmacokinet Pharmacodyn.",
    "2024;51(3):265-277. doi:10.1007/s10928-024-09901-2.",
    "The CAV input comes from the brepocitinib population PK model of Hughes",
    "JH et al. Clin Pharmacol Drug Dev. 2022;11(12):1447-1456.",
    "doi:10.1002/cpdd.1163; see modellib('Hughes_2022_brepocitinib')."
  )
  vignette <- "Tsamandouras_2024_brepocitinib"
  units <- list(
    time = "day",
    dosing = "n/a (no dose events; exposure enters as the covariate CAV)",
    concentration = "pasi (PASI score on its 0-72 scale; CAV and IC50 are in ng/mL)"
  )

  covariateData <- list(
    CAV = list(
      description = paste(
        "Individual average steady-state brepocitinib plasma concentration",
        "over the dosing interval of the current treatment period."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tsamandouras 2024 Eq. 1: CAV = Dose / (CL/F * tau), with tau = 24 h",
        "for once-daily and 168 h for once-weekly dosing and CL/F the",
        "empirical Bayes estimate from the Hughes 2022 population PK model",
        "(modellib('Hughes_2022_brepocitinib'); CL/F = 18.7 L/h at 70 kg,",
        "weight exponent 0.75, 78% CV). Time-varying: one value for the",
        "4-week induction period and one for the 8-week maintenance period,",
        "because brepocitinib (half-life 3.8-7.5 h) is assumed to be at",
        "steady state at every PASI assessment. Set to 0 during placebo",
        "periods. With dose in mg and CL/F in L/h, multiply by 1000 to get",
        "ng/mL."
      ),
      source_name = "CAVE"
    )
  )

  compartmentData <- list(
    latent_pasi = list(
      analyte = "latent drug-effect variable R(t)",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 212L,
    n_studies = 1L,
    age_range = "18-75 years (median 48)",
    weight_range = "45.1-204.3 kg (median 91.6)",
    sex_female_pct = 30.2,
    race_ethnicity = "White 89.2%, Black or African American 5.2%, Asian 1.9%, other 3.8%",
    baseline_pasi = "median 18.2 (range 12-54) on the 0-72 scale; arm means 19.1-23.8 (Table 1)",
    disease_state = "Adults with moderate-to-severe plaque psoriasis (PASI >= 12, PGA 3 or 4, >= 10% body surface area involved) for at least 6 months.",
    dose_range = "Phase 2a NCT02969018: 4-week induction with brepocitinib 30 or 60 mg QD or placebo, then 8-week maintenance with 10 or 30 mg QD, 100 mg once weekly or placebo (eight arms, Table 1).",
    regions = "Multicenter (countries not stated in the paper)",
    notes = "Tsamandouras 2024 Results 'Data description': 212 randomized patients, 148 men (69.8%), 170 completed the 12-week treatment period. Only PASI scores from weeks 0-12 entered the analysis; no covariates were tested."
  )

  ini({
    # Tsamandouras 2024 Table 2 (final model, NONMEM 7.3 Laplacian) and
    # Online Resource 1 control stream run6.mod.

    # Baseline PASI: THETA(1), logit-normal between the bounds 12 (the
    # inclusion criterion) and 72 (the scale maximum), Eq. 14.
    logitrbase <- logit((18.00 - 12) / (72 - 12)); label("Logit of the typical baseline PASI rescaled to its (12, 72) bounds (logit scale; BSL = 18.00 on the PASI scale)") # Table 2 BSL = 18.00 (RSE 2.69%); Online Resource 1 LGBASE = LOG((TVBASE-12)/(72-TVBASE))

    # Drug effect: THETA(2)-THETA(4), Eq. 11-12.
    imax <- fixed(1); label("Maximum fractional inhibition of the production of the latent variable R (unitless)") # Online Resource 1 $THETA '1 FIX ; IMAX' (fixed to 1, included 'just for completeness')
    lic50 <- log(53.24); label("Average brepocitinib concentration giving half-maximal inhibition of the production of R, IC50 (ng/mL)") # Table 2 IC50 = 53.24 ng/mL (RSE 21.43%)
    lkout <- log(0.067); label("First-order turnover rate constant of the latent variable R, kout = kin (1/day)") # Table 2 kout = 0.067 Day-1 (RSE 13.42%); footnote 'kin = kout'

    # Placebo effect: THETA(5)-THETA(6), Eq. 10.
    logitpmax <- logit(0.489); label("Logit of the maximum placebo effect Pmax, a fraction of baseline PASI (logit scale; Pmax = 0.489)") # Table 2 Pmax = 0.489 (RSE 8.41%); Online Resource 1 LGPMAX = LOG(TVPMAX/(1-TVPMAX))
    lkpbo <- log(0.047); label("First-order rate constant of onset of the placebo effect, kp (1/day)") # Table 2 kp = 0.047 Day-1 (RSE 6.48%)

    # Zero-inflated beta likelihood: THETA(7)-THETA(9), Eq. 3-8 and 13.
    lphi <- log(31.25); label("Precision phi of the beta distribution of PASI/72 (unitless)") # Table 2 phi = 31.25 (RSE 8.74%)
    logite0_pasi0 <- 2.228; label("Intercept zeta1 of the logit of the probability of PASI = 0 (logit scale)") # Table 2 zeta1 = 2.228 (RSE 39.95%)
    e_mu_pasi0 <- 62.67; label("Slope zeta2 by which the logit of the probability of PASI = 0 falls per unit of the beta mean mu (logit units)") # Table 2 zeta2 = 62.67 (RSE 20.24%)

    # Online Resource 1 $OMEGA: only ETA(1) on BASE is estimated; the etas on
    # IMAX, IC50, KOUT, PMAX and KP are '0 FIX' and are omitted here.
    etalogitrbase ~ 1.655 # Table 2 'etaBSL omega2' = 1.655 (RSE 12.07%); additive on the logit scale, Eq. 14
  })
  model({
    # Individual parameters (Online Resource 1 $PK). The baseline bounds 12
    # and 72 are fixed design constants of Eq. 14, not estimates.
    logitrbase_ind <- logitrbase + etalogitrbase
    rbase <- 12 + (72 - 12) * expit(logitrbase_ind)
    ic50 <- exp(lic50)
    kout <- exp(lkout)
    kin <- kout
    pmax <- expit(logitpmax)
    kpbo <- exp(lkpbo)
    phi <- exp(lphi)

    # Type I indirect response on the latent variable R, Eq. 12, with
    # R(0) = 1 and kin = kout. CAV is in ng/mL and t is days since the
    # first dose.
    latent_pasi(0) <- 1
    d/dt(latent_pasi) <- kin * (1 - imax * CAV / (ic50 + CAV)) - kout * latent_pasi

    # Placebo (Eq. 10) and drug (Eq. 11) effects in PASI units; Eq. 9 gives
    # the beta mean on the 0-1 scale. pasi_pred is the paper's IPRED * 72.
    placebo_eff <- rbase * pmax * (1 - exp(-kpbo * t))
    drug_eff <- rbase * (1 - pmax) * (1 - latent_pasi)
    pasi_pred <- rbase - placebo_eff - drug_eff
    mu_pasi <- pasi_pred / 72

    # Beta shapes (Eq. 6-7) and the probability of PASI = 0 (Eq. 13).
    alpha_pasi <- mu_pasi * phi
    beta_pasi <- (1 - mu_pasi) * phi
    p0_pasi <- expit(logite0_pasi0 - e_mu_pasi0 * mu_pasi)

    # Log of the gamma function by the Nemes approximation (Eq. 5), written
    # as in Online Resource 1, including its 3.1415 for pi. alpha + beta
    # equals phi.
    lgam_phi <- phi * (log(phi) - 1) + 0.5 * (log(2 * 3.1415) - log(phi)) +
      1.25 * phi * log(1 + 1 / (15 * phi^2))
    lgam_alpha <- alpha_pasi * (log(alpha_pasi) - 1) +
      0.5 * (log(2 * 3.1415) - log(alpha_pasi)) +
      1.25 * alpha_pasi * log(1 + 1 / (15 * alpha_pasi^2))
    lgam_beta <- beta_pasi * (log(beta_pasi) - 1) +
      0.5 * (log(2 * 3.1415) - log(beta_pasi)) +
      1.25 * beta_pasi * log(1 + 1 / (15 * beta_pasi^2))

    # Zero-inflated beta log-likelihood of the observed PASI score (Eq. 3-4;
    # Online Resource 1 $ERROR). DV is the PASI score on its 0-72 scale and
    # y = DV / 72 is the transformed value the source fitted (Eq. 2), so the
    # objective function equals the source's -2LL.
    y_pasi <- DV / 72
    if (DV == 0) {
      ll_pasi <- log(p0_pasi)
    } else {
      ll_pasi <- lgam_phi - lgam_alpha - lgam_beta +
        (alpha_pasi - 1) * log(y_pasi) + (beta_pasi - 1) * log(1 - y_pasi) +
        log(1 - p0_pasi)
    }
    ll(pasi) ~ ll_pasi
  })
}
