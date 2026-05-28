Zhang_2018_flurbiprofen <- function() {
  description <- "One-compartment IV population PK plus Holford-Sheiner effect-compartment for cerebrospinal fluid (CSF) disposition of flurbiprofen, the active metabolite of flurbiprofen axetil, in Chinese adults with postoperative pain receiving 1 mg/kg IV flurbiprofen axetil (Zhang 2018, Tables 1-2, Eq. 3 covariate form). Final-model typical values CL = 1.55 L/h, Vd = 7.91 L, plasma-CSF equilibration rate Ke = 0.0015/h; linear-multiplicative covariate effects of weight and height on Ke centered on the population medians (68.5 kg, 165 cm)."
  reference <- "Zhang J, Zhang H, Zhao L, Gu J, Feng Y, An H. Population pharmacokinetic modeling of flurbiprofen, the active metabolite of flurbiprofen axetil, in Chinese patients with postoperative pain. J Pain Res. 2018;11:3061-3070. doi:10.2147/JPR.S176475"
  vignette <- "Zhang_2018_flurbiprofen"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Centered on population median 68.5 kg (Zhang 2018 Table 1) inside Eq. 3 linear-multiplicative form on Ke.",
      source_name        = "WT"
    ),
    HT = list(
      description        = "Body height at baseline",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Centered on population median 165 cm (Zhang 2018 Table 1) inside Eq. 3 linear-multiplicative form on Ke.",
      source_name        = "HT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72L,
    n_studies      = 1L,
    age_range      = "18-72 years",
    age_median     = "50 years (mean 50.50, SD 13.63)",
    weight_range   = "45-96 kg",
    weight_median  = "68.50 kg (mean 67.47, SD 11.46)",
    height_range   = "149-185 cm",
    height_median  = "165 cm (mean 164.50, SD 9.08)",
    bmi_range      = "16.94-34.08 kg/m^2",
    bmi_median     = "24.33 kg/m^2 (mean 24.89, SD 3.49)",
    sex_female_pct = 62.5,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Adults undergoing surgery under subarachnoid anesthesia with postoperative pain. Approved by Peking University People's Hospital Medical Ethics Committee (ChiCTR-TRC-11001791).",
    dose_range     = "Single IV injection of 1 mg/kg flurbiprofen axetil (5050E; Tide Pharmaceutical) -- observed dose range 45-96 mg.",
    regions        = "China (single centre, Peking University People's Hospital, Beijing).",
    n_observations = "144 therapeutic drug-monitoring samples (72 plasma + 72 CSF) from the same 72 patients; each subject contributed one plasma sample and one CSF sample drawn simultaneously at a single time point. Patients were randomised into nine groups of 8 patients, and each group sampled at a fixed nominal time (5, 10, 15, 20, 25, 30, 35, 40, or 45 minutes after dose). Plasma concentrations 3.48-14.56 ug/mL; CSF concentrations 0-20.80 ng/mL.",
    notes          = "Sparse single-time-point joint plasma+CSF sampling design with stratified allocation across nine post-dose times. Sex 45 female / 27 male. CYP2C9 genotype was not collected; the authors note *1/*1 is the dominant Chinese genotype with *1/*3 and *1/*13 each <10%."
  )

  ini({
    # Structural PK parameters - final model, Zhang 2018 Table 2
    lcl  <- log(1.55);   label("Clearance CL (L/h)")                                   # Zhang 2018 Table 2 final model: theta_CL = 1.55 L/h (95% CI 1.543-1.550)
    lvc  <- log(7.91);   label("Volume of distribution Vd (L)")                        # Zhang 2018 Table 2 final model: theta_Vd = 7.91 L (95% CI 7.83-7.99)
    lkeo <- log(0.0015); label("Plasma-CSF equilibration rate constant Ke (1/h)")      # Zhang 2018 Table 2 final model: theta_Ke = 0.0015 /h (95% CI 0.0015-0.0015); effect-compartment equilibration rate (paper labels this Ke but it is the ke0 of the effect-compartment link model, not plasma elimination -- plasma elimination = CL/Vd = 0.196 /h)

    # Covariate-effect coefficients on Ke - linear-multiplicative form (Zhang 2018 Eq. 3 and final model equations)
    # Ke = theta_Ke * [1 + theta_weight * (WT - 68.5)] * [1 + theta_height * (HT - 165)] * exp(eta_Ke)
    e_wt_keo <- -0.0080; label("Linear coefficient on (WT - 68.5 kg) for Ke (1/kg)")   # Zhang 2018 Table 2 final model: theta_weight = -0.0080 (95% CI -0.0081 to -0.0080)
    e_ht_keo <- -0.0162; label("Linear coefficient on (HT - 165 cm) for Ke (1/cm)")    # Zhang 2018 Table 2 final model: theta_height = -0.0162 (95% CI -0.0162 to -0.0162)

    # IIV - exponential (log-normal) per Zhang 2018 Methods Eq. 1: P_i = theta * exp(eta_i).
    # Paper reports omega values as CV (fraction) under the Estimate column in Table 2; the row labels include "(%)".
    # Convert to internal log-scale variance via omega^2 = log(1 + CV^2).
    etalcl  ~ 0.5217   # CV(CL) = 82.75% (Zhang 2018 Table 2 final omega_CL = 0.8275)        -> omega^2 = log(1 + 0.8275^2) = 0.5217
    etalvc  ~ 0.00789  # CV(Vd) =  8.90% (Zhang 2018 Table 2 final omega_Vd = 0.0890)        -> omega^2 = log(1 + 0.0890^2) = 0.00789
    etalkeo ~ 0.000484 # CV(Ke) =  2.20% (Zhang 2018 Table 2 final omega_Ke = 0.0220)        -> omega^2 = log(1 + 0.0220^2) = 0.000484

    # Residual error - additive (Zhang 2018 Methods Eq. 2): C_obs = C_pred + eps, eps ~ N(0, sigma^2).
    # A single sigma = 0.0023 mg/L was reported for the joint plasma+CSF fit; the paper labels
    # eps on the "serum flurbiprofen concentration" but the same additive sigma was applied to
    # the CSF observations in the effect-compartment link fit. nlmixr2 requires distinct
    # residual-error parameters per endpoint, so the single source sigma is replicated as both
    # addSd (plasma) and addSd_Ccsf (CSF). Magnitude (0.0023 mg/L = 2.3 ng/mL) is commensurate
    # with the CSF observation range (0-20.8 ng/mL) and tiny relative to plasma (3.48-14.56 mg/L),
    # which is why plasma IPRED vs DV (Zhang 2018 Figure 1A left) lies almost exactly on the
    # identity line.
    addSd      <- 0.0023; label("Additive residual error on plasma Cc (mg/L)")   # Zhang 2018 Table 2 final model: sigma = 0.0023 mg/L (95% CI 0.0023-0.0023); shared with addSd_Ccsf
    addSd_Ccsf <- 0.0023; label("Additive residual error on CSF Ccsf (mg/L)")    # Zhang 2018 Table 2 final model: sigma = 0.0023 mg/L; same source sigma applied to CSF endpoint
  })
  model({
    # Individual PK parameters
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)

    # Effect-site equilibration rate with linear-multiplicative covariate effects (Zhang 2018 Eq. 3 + final equations).
    # Median weight = 68.5 kg, median height = 165 cm (Zhang 2018 Table 1).
    keo <- exp(lkeo + etalkeo) *
           (1 + e_wt_keo * (WT - 68.5)) *
           (1 + e_ht_keo * (HT - 165))

    # Micro-rate constant
    kel <- cl / vc

    # One-compartment IV disposition (dosing into central as IV bolus)
    d/dt(central) <- -kel * central
    Cc            <- central / vc

    # Effect compartment representing CSF (Holford-Sheiner link form):
    # dCe/dt = keo * (Cc - Ce); Ce is the apparent CSF concentration.
    d/dt(effect) <- keo * (Cc - effect)
    Ccsf         <- effect

    # Observation models - single additive sigma applied to both outputs per Zhang 2018 Eq. 2;
    # nlmixr2 requires per-endpoint residual-error parameters, so addSd and addSd_Ccsf carry
    # the same paper-reported initial value.
    Cc   ~ add(addSd)
    Ccsf ~ add(addSd_Ccsf)
  })
}
