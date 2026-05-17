Koloskoff_2025_ganciclovir <- function() {
  description <- "Indirect-response viral turnover PD model for cytomegalovirus (CMV) viral load decline in pediatric solid-organ and hematopoietic-stem-cell transplant recipients receiving (val)ganciclovir (Koloskoff 2025). The model treats the q12h-interval ganciclovir AUC (AUC_0-12) as a time-varying covariate input AUC_GCV that stimulates first-order viral degradation through an Emax-EC50 relationship. The upstream popPK that produces AUC_0-12 (Franck 2021 Bayesian estimator) is NOT included here; AUC_GCV must be supplied per record by the user, either from the Franck 2021 model or any other AUC source."
  reference   <- paste(
    "Koloskoff K, Franck B, Benito S, Welzel J, Autmizguine J, Theoret Y, Briand A,",
    "Ovetchkine P, Woillard J-B.",
    "Pharmacokinetic/Pharmacodynamic Modelling and Monte Carlo Simulations to Predict",
    "Cytomegalovirus Viral Load in Pediatric Transplant Recipients Treated with",
    "(val)Ganciclovir.",
    "Clin Pharmacokinet. 2025.",
    "doi:10.1007/s40262-025-01526-z.",
    "Upstream popPK used by the source authors to compute AUC_0-12 inputs:",
    "Franck B, Autmizguine J, Asberg A, Theoret Y, Marquet P, Ovetchkine P, et al.",
    "Thoroughly validated Bayesian estimator and limited sampling strategy for dose",
    "individualization of ganciclovir and valganciclovir in pediatric transplant",
    "recipients. Clin Pharmacokinet. 2021;60:1449-1462.",
    "doi:10.1007/s40262-021-01034-w.",
    sep = " "
  )
  vignette    <- "Koloskoff_2025_ganciclovir"
  units       <- list(
    time          = "hour",
    dosing        = "n/a (no drug dosing events; ganciclovir exposure enters as the time-varying AUC_GCV covariate)",
    concentration = "log10 copies/mL (CMV viral load output; not a drug concentration)",
    AUC_GCV       = "mg*h/L (per q12h dosing interval; q24h regimens are entered as AUC_0-24 / 2 so all data live in a q12h framework)"
  )

  covariateData <- list(
    AUC_GCV = list(
      description        = "Time-varying ganciclovir AUC over a q12h dosing interval (AUC_0-12) used as the drug-exposure input to the viral-load PD model.",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying drug-exposure covariate. The source authors compute AUC_0-12 from",
        "an upstream popPK model (Franck 2021, doi:10.1007/s40262-021-01034-w) and feed",
        "it as a varying input ('amount') into the Monolix PD model. This nlmixr2lib",
        "model expects pre-computed AUC values supplied as a data column; AUC_GCV is not",
        "computed inside the model. Q24h dosing intervals are entered as AUC_0-24 / 2 so",
        "all data live in a q12h framework (Koloskoff 2025 Methods, Section 2.1). Set to",
        "0 in pre-treatment or off-treatment records so the drug-stimulation term",
        "vanishes and the viral load returns to the kin / kout steady-state baseline."
      ),
      source_name        = "AUC_0-12"
    )
  )

  population <- list(
    species                    = "human",
    n_subjects                 = 29L,
    n_occurrences              = 36L,
    n_observations             = 184L,
    n_studies                  = 1L,
    age_range                  = "0.5-15 years",
    age_median                 = "8.2 years",
    weight_range               = "6.3-95.2 kg",
    weight_median              = "29.4 kg",
    height_range               = "41-172 cm",
    height_median              = "128 cm",
    sex_female_pct             = 41.7,
    serum_creatinine_median    = "46.0 umol/L (range 9-414)",
    crcl_median                = "111 mL/min/1.73 m^2 (Schwartz-modified; range 24.8-243)",
    disease_state              = "Pediatric solid-organ transplant (SOT, n = 18 occurrences) or hematopoietic stem cell transplant (HSCT, n = 18 occurrences) recipients monitored for CMV reactivation via weekly DNAemia screening; 6 occurrences had graft-versus-host disease (GVHD).",
    dose_range                 = "Pre-emptive treatment with IV ganciclovir 5 mg/kg q12h or oral valganciclovir 10 mg/kg q12h, adjusted by TDM. Initial AUC_0-12 median 23.5 mg*h/L (range 4.30-47.3); post-TDM AUC_0-12 median 23.7 mg*h/L (range 6.04-83.7).",
    regions                    = "Canada (CHU Sainte-Justine, Montreal, QC).",
    baseline_viral_load        = "Median 3.61 log10 copies/mL (range 2.57-4.85).",
    treatment_duration_median  = "22 days (range 6-76)",
    transplant_types           = "SOT (liver, kidney, heart) and HSCT (allogeneic). Six occurrences had GVHD; binary and time-dependent GVHD were tested as covariates and not retained.",
    notes                      = paste(
      "Retrospective single-center cohort, January 2007 - December 2015. Inclusion:",
      "SOT or HSCT receiving valganciclovir and/or IV ganciclovir for CMV disease",
      "prevention with at least one full PK profile and at least two CMV viral loads.",
      "Patients on foscarnet were excluded. 29 children met criteria; 5 had multiple",
      "distinct CMV episodes (>=3 month gap) counted as separate occurrences, giving",
      "36 occurrences and 184 viral-load observations of which 42 were below the",
      "LLOQ (200 copies/mL) and handled with the NONMEM M4 method (Bergstrand 2009).",
      "AUC_0-12 was computed externally from the Franck 2021 popPK model and fed into",
      "the Monolix PD model as a varying input. Type of transplant, time since",
      "transplantation, and GVHD were tested as covariates and not retained. Sources:",
      "Koloskoff 2025 Table 1 (demographics) and Methods Section 2.1 (cohort)."
    )
  )

  ini({
    # PD parameters from Koloskoff 2025 Table 2 'Final model estimate' column
    # (Monolix indirect viral turnover model with drug stimulation of viral
    # degradation; n = 36 occurrences / 29 children / 184 observations). Time
    # in hours, viral load R in log10 copies/mL, AUC_GCV in mg*h/L (per q12h
    # interval). The Eq. 1 ODE is
    #   dR/dt = kin - kout * (1 + Emax * AUC / (EC50 + AUC)) * R
    # with R(0) = kin / kout. The typical-value steady-state baseline is
    # 0.00087 / 0.00023 = 3.78 log10 copies/mL, consistent with the
    # observed median baseline of 3.61 log10 copies/mL (Table 1).
    lkin   <- log(0.00087); label("Zero-order viral production rate kin (log10 copies/mL per hour)")          # Koloskoff 2025 Table 2 final kin = 0.00087 (RSE 4.80%)
    lkout  <- log(0.00023); label("First-order viral elimination rate constant kout (1/hour)")                  # Koloskoff 2025 Table 2 final kout = 0.00023 (RSE 5.19%)
    lemax  <- log(16.3);    label("Maximum drug-induced fold-increase in viral elimination Emax (unitless)")    # Koloskoff 2025 Table 2 final Emax = 16.3 (RSE 18.1%)
    lec50  <- log(23.5);    label("Ganciclovir AUC at half-maximal stimulation EC50 (mg*h/L)")                  # Koloskoff 2025 Table 2 final EC50 = 23.5 mg*h/L (RSE 46.7%)

    # IIV from Koloskoff 2025 Table 2 'Final' column. Monolix reports omega as
    # the standard deviation of the log-scale random effect; nlmixr2 stores
    # the variance, so var = omega^2. Koloskoff 2025 Results: 'The
    # interindividual variabilities maximum effect and k_in were removed as
    # they had no impact on individual fit', so only kout and EC50 carry
    # between-subject variability.
    etalkout ~ 0.0196    # Koloskoff 2025 Table 2 omega_kout = 0.14 -> var = 0.14^2 = 0.0196 (RSE 19.9%)
    etalec50 ~ 1.8496    # Koloskoff 2025 Table 2 omega_EC50 = 1.36 -> var = 1.36^2 = 1.8496 (RSE 24.5%)

    # Residual error: Monolix proportional ('b') parameterisation
    # Y = ypred * (1 + b * eps) with eps ~ N(0, 1). Koloskoff 2025 Table 2
    # reports b = 0.15 on the log10 viral-load observations (RSE 8.66%).
    propSd <- 0.15; label("Proportional residual error on log10 viral load (fraction)") # Koloskoff 2025 Table 2 final b = 0.15 (RSE 8.66%)
  })

  model({
    # Individual structural parameters. IIV is present on kout and EC50
    # only; kin and Emax are typical-value-only per Koloskoff 2025 Results.
    kin  <- exp(lkin)
    kout <- exp(lkout + etalkout)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)

    # Indirect-response viral turnover with drug stimulation of viral
    # degradation (Koloskoff 2025 Eq. 1; reproduced from Cojutti 2018 model
    # structure with AUC_0-12 replacing instantaneous concentration). The
    # AUC_GCV covariate must be present in the input dataset as a
    # time-varying column in mg*h/L. The state viralLoad is on the
    # log10 copies/mL scale; the typical-value steady-state baseline at
    # AUC_GCV = 0 is kin / kout = 3.78 log10 copies/mL.
    d/dt(viralLoad) <- kin - kout * (1 + emax * AUC_GCV / (ec50 + AUC_GCV)) * viralLoad
    viralLoad(0)    <- kin / kout

    viralLoad ~ prop(propSd)
  })
}
