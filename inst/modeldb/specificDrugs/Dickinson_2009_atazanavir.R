Dickinson_2009_atazanavir <- function() {
  description <- "One-compartment first-order-absorption population PK model with absorption lag-time for oral ritonavir-boosted atazanavir in HIV-infected adults and healthy volunteers; ritonavir AUC0-24 (median 7.52 mg*h/L) enters CL/F via a power function (Dickinson 2009)."
  reference <- "Dickinson L, Boffito M, Back D, Waters L, Else L, Davies G, Khoo S, Pozniak A, Aarons L. Population pharmacokinetics of ritonavir-boosted atazanavir in HIV-infected patients and healthy volunteers. J Antimicrob Chemother. 2009;63(6):1233-1243. doi:10.1093/jac/dkp102"
  vignette <- "Dickinson_2009_atazanavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AUC_RTV = list(
      description        = "Ritonavir AUC over the 0-24 h dosing interval (per-subject, time-fixed within an evaluated regimen)",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters atazanavir CL/F via the power form CL = exp(lcl) * (AUC_RTV / 7.52)^e_aucrtv_cl, centred at the cohort median 7.52 mg*h/L (Dickinson 2009 Table 1 / Results page 1236). Computed by the source authors using non-compartmental methods (WinNonlin 5.2) on the ritonavir concentration-time data. Set per-subject; for simulation users without observed ritonavir AUC, the cohort median (7.52) reproduces typical-value behaviour. Healthy volunteers 7.36 mg*h/L (range 4.31-13.42); HIV-infected 7.59 mg*h/L (range 2.41-22.05); pooled 7.52 (2.41-22.05) per Table 1.",
      source_name        = "RTVAUC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 46L,
    n_studies      = 3L,
    n_observations = 538L,
    age_range      = "22-62 years",
    age_median     = "43 years",
    weight_range   = "46-115 kg",
    weight_median  = "76 kg",
    sex_female_pct = 19.6,
    race_ethnicity = c(Caucasian = 72, BlackAfrican = 15, Hispanic = 13),
    disease_state  = "HIV-infected adults (n=30) and healthy volunteers (n=16); stable on atazanavir/ritonavir for >=2 weeks prior to PK sampling.",
    dose_range     = "Oral atazanavir/ritonavir 300/100 mg once daily, fed (16-20 g fat). 18 of 46 also received saquinavir 1600 mg once daily; 6 of 46 received tenofovir 300 mg once daily. Lower-dose regimens (200/100 and 150/100 mg once daily) used only for external validation.",
    regions        = "United Kingdom (St Stephen's Centre, Chelsea and Westminster Foundation Trust, London)",
    notes          = "Three pooled clinical studies in adults. Sampling pre-dose and 0.5, 1, 2, 3, 4, 6, 8, 10, 12, 24 h post-dose; healthy volunteers had additional samples at 16 and 20 h. Plasma atazanavir and ritonavir quantified by HPLC-MS/MS. Median (range) ritonavir AUC0-24 was 7.52 mg*h/L (2.41-22.05) (Table 1). Model fit with NONMEM VI 2.0 (FOCE-I) per Methods / Data analysis."
  )

  ini({
    # Structural parameters from Dickinson 2009 Table 2 (final-model column).
    lcl   <- log(7.7);  label("Apparent oral clearance at AUC_RTV = 7.52 mg*h/L (CL/F, L/h)")  # Table 2 final: CL/F = 7.7 L/h (RSE 5%)
    lvc   <- log(103);  label("Apparent volume of distribution (V/F, L)")                      # Table 2 final: V/F = 103 L (RSE 13%)
    lka   <- log(3.4);  label("First-order absorption rate constant (ka, 1/h)")                # Table 2 final: ka = 3.4 1/h (RSE 34%)
    ltlag <- log(0.96); label("Absorption lag-time (Tlag, h)")                                 # Table 2 final: Lag-time = 0.96 h (RSE 1%)

    # Covariate effect: ritonavir AUC0-24 on atazanavir CL/F via power function.
    # Paper Results page 1236 / Table 3 final equation: CL/F = theta1 * (RTV/7.52)^theta2
    # with the cohort-median ritonavir AUC0-24 of 7.52 mg*h/L as the centring reference.
    e_aucrtv_cl <- -0.8; label("Power exponent of ritonavir AUC0-24 on atazanavir CL/F (unitless)")  # Table 2 final: -0.8 (RSE 13%)

    # Inter-individual variability (exponential model on CL/F, V/F, ka).
    # Paper reports IIV (%) as CV; convert to log-normal variance via omega^2 = log(1 + CV^2).
    # Lag-time IIV was not significant (delta-OFV -0.8) and was not included in the final model.
    etalcl ~ 0.0807  # log(1 + 0.29^2) = 0.0807; CV = 29% per Table 2 final (RSE 59%)
    etalvc ~ 0.2074  # log(1 + 0.48^2) = 0.2074; CV = 48% per Table 2 final (RSE 37%)
    etalka ~ 1.2155  # log(1 + 1.54^2) = 1.2155; CV = 154% per Table 2 final (RSE 51%)

    # Residual error: combined proportional + additive (paper Results page 1236).
    propSd <- 0.23; label("Proportional residual error (CV, fraction)")  # Table 2 final: 23% (RSE 27%)
    addSd  <- 0.08; label("Additive residual error (mg/L)")              # Table 2 final: 0.08 mg/L (RSE 38%)
  })

  model({
    # Individual PK parameters
    cl   <- exp(lcl + etalcl) * (AUC_RTV / 7.52)^e_aucrtv_cl
    vc   <- exp(lvc + etalvc)
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    # Micro-constant
    kel <- cl / vc

    # ODE system (1-compartment with first-order absorption from depot)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Absorption lag-time on the depot compartment
    alag(depot) <- tlag

    # Observation and combined proportional + additive residual error
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
