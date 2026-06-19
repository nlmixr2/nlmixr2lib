Ahmed_2016_lorenzosOil <- function() {
  description <- "Population pharmacodynamic model of Lorenzo's oil effect on plasma C26:0 in asymptomatic boys with X-linked adrenoleukodystrophy: inhibitory fractional Emax model relating observed plasma erucic acid concentration to plasma C26:0. The paper does not develop a PK model for erucic acid; observed erucic acid plasma concentration is supplied as a time-varying covariate."
  reference <- "Ahmed MA, Kartha RV, Brundage RC, Cloyd J, Basu C, Carlin BP, Jones RO, Moser AB, Fatemi A, Raymond GV. A model-based approach to assess the exposure-response relationship of Lorenzo's oil in adrenoleukodystrophy. Br J Clin Pharmacol. 2016 Jun;81(6):1058-1065. doi:10.1111/bcp.12897"
  vignette <- "Ahmed_2016_lorenzosOil"
  units <- list(
    time = "year",
    dosing = "mg/kg/day",
    concentration = "mg/L"
  )

  covariateData <- list(
    CP_ER_MGL = list(
      description        = "Instantaneous plasma erucic acid concentration as a time-varying PD driver",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying input. Observed plasma erucic acid concentration in X-ALD patients receiving Lorenzo's oil. The paper does not develop a PK model for erucic acid; observed values are supplied directly per observation record. Median pretreatment 0.5 mg/L (range 0.22-1.94, n=97); median post-treatment 18.63 mg/L (range 0.21-336.1) (Ahmed 2016 Table 1). Values > 30 mg/L were censored in the paper's time-to-event analysis as transient peak post-dose concentrations not reflective of steady-state brain exposure.",
      source_name        = "ER"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 104L,
    n_studies       = 1L,
    age_range       = "0.068-8.92 years",
    age_median      = "2.79 years",
    weight_range    = "9.40-40.60 kg",
    weight_median   = "14.90 kg",
    sex_female_pct  = 0,
    disease_state   = "Asymptomatic X-linked adrenoleukodystrophy (X-ALD); diagnosis confirmed by SVLCFAs assay. Subjects with baseline neurological or radiological abnormality, brain tumour, another peroxisomal disorder, or non-adherence to LO were excluded.",
    dose_range      = "Approximately 2-3 mg/kg/day of Lorenzo's oil providing 20% of caloric intake; halted and replaced by glyceryl trioleate at the same dosage on clinically significant platelet drop, then resumed at a lower dose with gradual return to target",
    regions         = "John Hopkins Research Hospital (United States); enrollment 2000-2014",
    notes           = "Open-label single-arm trial (ClinicalTrials.gov NCT02233257). 2384 paired C26:0 and erucic acid plasma measurements over a mean follow-up of 4.88 +/- 2.76 years (range 0-10.26). Sampling: baseline plus monthly for the first 6 months, then every 3-6 months. Pretreatment and post-treatment plasma C26:0 medians were 1.06 mg/L and 0.402 mg/L respectively (Table 1). Estimation in NONMEM 7.3 with FOCE-I; final model evaluated by prediction- and variability-corrected VPC (Figure 1) and 1000-replicate non-parametric bootstrap (Table 2). The accompanying Weibull time-to-event analyses (hazard of brain MRI abnormality vs. LAUCER and LAUCC26:0) were conducted in SAS PROC LIFEREG using per-subject summary LAUC values; those analyses are an auxiliary regression on a summary statistic and are not extracted as a separate nlmixr2 model file. See the vignette for the AFT-Weibull point estimates (Tables 3, 4) and a discussion of why they sit outside the popPD scope."
  )

  ini({
    # Structural fixed effects (Ahmed 2016 Table 2). Final estimates; bootstrap medians match the
    # point estimates within rounding (E0 1.44 vs 1.44; Emax 0.76 vs 0.761; EC50 0.734 vs 0.733).
    le0   <- log(1.44);  label("Baseline C26:0 plasma concentration at zero erucic acid (E0, mg/L)")    # Ahmed 2016 Table 2
    lemax <- log(0.76);  label("Fractional maximum reduction of C26:0 by erucic acid (Emax, unitless in [0,1])")  # Ahmed 2016 Table 2
    lec50 <- log(0.734); label("Erucic acid concentration at half-maximum effect (EC50, mg/L)")          # Ahmed 2016 Table 2

    # IIV: log-normal (exponential) model. Ahmed 2016 Table 2 reports CV% computed as
    # sqrt(exp(omega^2) - 1) (footnote *), so omega^2 = log(CV^2 + 1).
    # Correlation(E0, EC50) = -0.877 estimated; correlations between Emax and the other
    # fixed effects were "very small and, therefore, were fixed to zero in the final model"
    # (Results, Population PD analyses). The block on (E0, EC50) carries the off-diagonal
    # covariance; Emax sits on its own diagonal eta.
    etale0 + etalec50 ~ c(
      0.09464,                  # omega^2_E0   = log(0.315^2 + 1)                                        # Ahmed 2016 Table 2 BSV_E0 31.5%
     -0.31579,                  # cov(E0, EC50) = -0.877 * sqrt(0.09464 * 1.36998)                       # Ahmed 2016 Table 2 corr(E0,EC50) -0.877
      1.36998                   # omega^2_EC50 = log(1.713^2 + 1)                                        # Ahmed 2016 Table 2 BSV_EC50 171.3%
    )
    etalemax ~ 0.003836         # omega^2_Emax = log(0.062^2 + 1)                                        # Ahmed 2016 Table 2 BSV_Emax 6.2%

    # Residual error: proportional only (Results: "A proportional error model best described
    # the residual errors for these data. The combined error model did not provide any
    # significant improvement of fit (Delta OFV < 3.84); the additive error model provided a
    # worse fit (Delta AIC ~ 671)."). Table 2 footnote: CV% = sqrt(sigma^2), so propSd = sigma.
    propSd <- 0.276;     label("Proportional residual error on C26:0 (fraction)")                        # Ahmed 2016 Table 2 (footnote dagger)
  })

  model({
    # Inhibitory fractional Emax (Ahmed 2016, Methods/Population PD model building):
    #   C26:0 = E0 * (1 - Emax * ER / (EC50 + ER))
    # where ER is the observed erucic acid plasma concentration (mg/L), supplied per
    # observation record via the time-varying covariate CP_ER_MGL. The model has no
    # ODE state, no compartment, and no dosing event in nlmixr2 terms; LO administration
    # enters the model only through its measured effect on erucic acid plasma exposure.
    e0   <- exp(le0   + etale0)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)

    Cc26 <- e0 * (1 - emax * CP_ER_MGL / (ec50 + CP_ER_MGL))

    Cc26 ~ prop(propSd)
  })
}
