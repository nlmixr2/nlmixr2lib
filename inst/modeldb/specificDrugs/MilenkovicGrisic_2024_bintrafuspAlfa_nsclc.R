MilenkovicGrisic_2024_bintrafuspAlfa_nsclc <- function() {
  description <- "Empirical tumor-size model (sum of longest diameters of target lesions, mm) for adults with previously untreated advanced non-small cell lung cancer with high PD-L1 expression (at least 80 percent PD-L1-positive tumor cells) receiving first-line bintrafusp alfa 1200 mg IV every 2 weeks (bintrafusp alfa arm of the randomised phase III INTR@PID Lung 037, NCT03631706). Closed-form profile TS(t) = BTS * exp(-ks * t) + kg * t with t in weeks from the first dose: exponential shrinkage of the baseline tumor plus linear (zero-order) growth. Log-normal IIV on BTS, ks and kg with a full 3 x 3 correlation block; residual error is additive on the log scale (log-transform-both-sides). Covariates: CRP on BTS and lymphocytes as a percentage of total leukocytes on kg; bintrafusp alfa exposure was not significant. The paper prints the coefficients but not the covariate equations or reference values; the power form and the reference values are maintainers' assumptions documented in the vignette."
  reference <- paste(
    "Milenkovic-Grisic A-M, Terranova N, Mould DR, Vugmeyster Y, Mrowiec T,",
    "Machl A, Girard P, Venkatakrishnan K, Khandelwal A. Tumor growth inhibition",
    "modeling in patients with second line biliary tract cancer and first line",
    "non-small cell lung cancer based on bintrafusp alfa trials.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(1):143-153.",
    "doi:10.1002/psp4.13068",
    sep = " "
  )
  vignette <- "MilenkovicGrisic_2024_tumor_dynamics"
  units <- list(time = "week", dosing = "mg", concentration = "mm (tumor size)")

  covariateData <- list(
    CRP = list(
      description = "Baseline C-reactive protein, standard assay (time-fixed).",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters BTS as a power effect (CRP / 10)^0.125. The paper prints neither the form, the units nor the reference value; the power form, mg/L and the 10 mg/L reference are maintainers' assumptions. Must be > 0.",
      source_name = "CRP"
    ),
    LYMPH_PCT = list(
      description = "Baseline lymphocytes as a percentage of total leukocytes (time-fixed).",
      units = "percent",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters kg as a power effect (LYMPH_PCT / 20)^-0.623. The paper prints neither the form nor the reference value; the power form and the 20 percent reference are maintainers' assumptions. Must be > 0.",
      source_name = "Lymphocytes as percentage of total leukocytes"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 139L,
    n_studies = 1L,
    age_range = "median 67 years (range 42-90)",
    weight_range = "median 66.5 kg (range 35.0-104.3)",
    sex_female_pct = 20.1,
    disease_state = "Previously untreated advanced non-small cell lung cancer with high PD-L1 tumor expression (at least 80 percent PD-L1-positive tumor cells); first line.",
    dose_range = "Bintrafusp alfa 1200 mg IV every 2 weeks.",
    regions = "Multicenter international adaptive phase III randomised controlled trial INTR@PID Lung 037 (NCT03631706), bintrafusp alfa arm.",
    sampling_window = "Tumor assessments every 8 weeks up to 12 months, then every 12 weeks (Table S1).",
    assay = "Tumor size = RECIST 1.1 sum of longest diameters of target lesions (mm), mainly CT or X-ray.",
    notes = "139 of 152 patients randomised to bintrafusp alfa (91.4%) had evaluable PK and tumor-size data; 111 men (79.9%). Hyperprogression (Matos criteria) in 11 patients. Missing covariates were imputed with the median (continuous) or mode (categorical)."
  )

  ini({
    # Typical values: Milenkovic-Grisic 2024 Table 2, bintrafusp alfa columns,
    # cross-checked against the Results text for ks and kg.
    lrbase <- log(81.1)   ; label("Typical baseline tumor size BTS (mm)")                         # Table 2 bintrafusp alfa: Baseline TS (mm) = 81.1, RSE 5.75
    lkse   <- log(0.0368) ; label("Typical exponential tumor-shrinkage rate constant ks (1/week)") # Table 2 bintrafusp alfa: ks (week-1) = 0.0368, RSE 12.8; Results text 0.0368
    lp     <- log(0.728)  ; label("Typical linear tumor-growth rate constant kg (mm/week)")        # Table 2 bintrafusp alfa: kg (mm/week) = 0.728, RSE 17.9; Results text 0.728

    e_crp_rbase    <-  0.125 ; label("Power exponent of CRP on BTS (unitless)")                                 # Table 2 bintrafusp alfa: Effect of CRP on baseline TS = 0.125, RSE 23.8
    e_lymphpct_p   <- -0.623 ; label("Power exponent of lymphocytes as a percentage of leukocytes on kg (unitless)") # Table 2 bintrafusp alfa: Effect of lymphocyte as percentage of total leukocytes on kg = -0.623, RSE 43.7

    # IIV: Table 2 reports %CV; omega^2 = log(1 + CV^2). Off-diagonals are
    # corr * sqrt(omega_i^2 * omega_j^2) from the Table 2 correlations.
    #   BTS 55.4 %CV -> 0.2677; ks 62.8 %CV -> 0.3325; kg 122.9 %CV -> 0.9205
    #   cov(BTS, ks) = -0.0411 * sqrt(0.2677 * 0.3325) = -0.01226
    #   cov(BTS, kg) =  0.423  * sqrt(0.2677 * 0.9205) =  0.2100
    #   cov(ks, kg)  = -0.0639 * sqrt(0.3325 * 0.9205) = -0.03535
    etalrbase + etalkse + etalp ~ c(0.2677,
                                    -0.01226, 0.3325,
                                    0.2100, -0.03535, 0.9205) # Table 2 bintrafusp alfa: IIV BTS 55.4, ks 62.8, kg 122.9 %CV; Corr(BTS,ks) -0.0411, Corr(BTS,kg) 0.423, Corr(ks,kg) -0.0639

    expSd <- 0.128 ; label("Additive residual error SD on the log scale (unitless)") # Table 2 bintrafusp alfa: Proportional residual unexplained variability = 12.8 %CV, RSE 2.29
  })

  model({
    # Power covariate effects with maintainers' reference values (CRP 10 mg/L,
    # lymphocytes 20 percent of leukocytes).
    rbase <- exp(lrbase + etalrbase) * (CRP / 10)^e_crp_rbase
    kse <- exp(lkse + etalkse)
    p <- exp(lp + etalp) * (LYMPH_PCT / 20)^e_lymphpct_p

    # Milenkovic-Grisic 2024 Eq 1: TS(t) = BTS * exp(-ks * t) + kg * t, with t
    # in weeks from the first bintrafusp alfa dose.
    tumor_size <- rbase * exp(-kse * t) + p * t
    tumor_size ~ lnorm(expSd)
  })
}
