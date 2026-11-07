Maskawa_2025_vancomycin_mouse_firstorder <- function() {
  description <- paste(
    "Preclinical (mouse).",
    "Gastrointestinal transit model for fecal concentrations of orally",
    "administered, non-absorbed vancomycin in male C57BL/6N mice, with",
    "first-order gastric emptying (the paper's supplementary Model 2,",
    "the alternative to the zero-order main model",
    "Maskawa_2025_vancomycin_mouse). The dose empties from the stomach",
    "into the cecum with rate constant k_a, is eliminated from the cecum",
    "by first-order transit (k_e), and the fecal concentration is the",
    "cecal amount divided by the mass of cecal contents A. The",
    "stomach-to-cecum and cecum-to-feces transit delays are only",
    "identifiable as their sum t_tot, applied as a lag on the oral dose",
    "into stomach. Fitted per mouse by penalised least squares; typical",
    "values are the arithmetic means over the six treated mice and the",
    "between-mouse variances are derived from the reported coefficients",
    "of variation.",
    sep = " "
  )
  reference <- paste(
    "Maskawa R, Takayasu L, Takayasu H, Watanabe K, Takemine S, Kakimoto T,",
    "Takeshita K, Narushima S, Suda W, Takayasu M.",
    "High-resolution fecal pharmacokinetic modeling in mice with orally",
    "administered antibiotics.",
    "Sci Rep. 2025;15:24441.",
    "doi:10.1038/s41598-025-10293-7.",
    "Model 2 from the Supplementary Information, section 2.",
    sep = " "
  )
  vignette <- "Maskawa_2025_vancomycin_mouse"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/g"
  )

  compartmentData <- list(
    # Supplement section 2: M_a(t), antibiotic mass in the stomach.
    stomach = list(analyte = "vancomycin", units = "mg", specimen = "administration site", verified = TRUE),
    # Supplement section 2: M_b(t), antibiotic mass in the cecum. As in
    # Model 1, C(t) = M_b(t - dt_CF) / A, so with the lag set to t_tot the
    # state is the cecal mass time-shifted to the moment it is expelled,
    # and cecum / A is the measured fecal concentration.
    cecum = list(analyte = "vancomycin", units = "mg", specimen = "faeces", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (C57BL/6N, male, conventional)",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "8 weeks",
    sex_female_pct = 0,
    disease_state = "Healthy conventional mice (no infection model)",
    dose_range = paste(
      "Single oral dose of 0.5 mL vancomycin hydrochloride in water at",
      "1 mg/mL (q = 0.5 mg; low dose, mice L1-L3) or 20 mg/mL (q = 10 mg;",
      "high dose, mice H1-H3), given at 12:00 on day 0."
    ),
    regions = "Japan (Keio University School of Medicine)",
    notes = paste(
      "Same six mice and fecal series as Maskawa_2025_vancomycin_mouse.",
      "Three treated mice per dose group (plus three water-only controls",
      "not used in the model fit). Body weight not reported. Methods",
      "'Mice and fecal sample collection'; Fig. 1."
    )
  )

  ini({
    # Typical values: Table S1 'Mean (CV%)' row, arithmetic mean over the
    # six treated mice (L1-L3, H1-H3), each fitted individually.
    ltlag <- log(6.7); label("Total gastrointestinal transit time of the antibiotic t_tot (h)") # Table S1, t_tot mean 6.7 h (CV 46%)
    lkst <- log(0.32); label("First-order gastric emptying rate constant k_a (1/h)") # Table S1, k_a mean 0.32 1/h (CV 53%)
    lkfec <- log(0.36); label("First-order elimination rate constant from the cecum to feces k_e (1/h)") # Table S1, k_e mean 0.36 1/h (CV 57%)
    lv_cecum <- log(0.40); label("Mass of cecal contents A (g)") # Table S1, A mean 0.40 g (CV 91%)

    # Between-mouse variability: omega^2 = log(1 + CV^2) from the Table S1
    # CV% of the six individual estimates.
    etaltlag ~ 0.1919 # Table S1, t_tot CV 46%: log(1 + 0.46^2)
    etalkst ~ 0.2476 # Table S1, k_a CV 53%: log(1 + 0.53^2)
    etalkfec ~ 0.2813 # Table S1, k_e CV 57%: log(1 + 0.57^2)
    etalv_cecum ~ 0.6033 # Table S1, A CV 91%: log(1 + 0.91^2)

    # Unweighted least-squares fit (Methods Eq. 8, lambda = 0.9), so the
    # residual error is additive; its SD is the pooled RMSE of Model 2.
    addSd_Ccecum <- 0.60; label("Additive residual error on fecal concentration (mg/g)") # Table S2, Model2 RMSE 'All' 6.0e-1
  })
  model({
    tlag <- exp(ltlag + etaltlag)
    kst <- exp(lkst + etalkst)
    kfec <- exp(lkfec + etalkfec)
    v_cecum <- exp(lv_cecum + etalv_cecum)

    # Supplement section 2: dM_a/dt = -k_a M_a; dM_b/dt = k_a M_a - k_e M_b
    d/dt(stomach) <- -kst * stomach
    d/dt(cecum) <- kst * stomach - kfec * cecum
    alag(stomach) <- tlag

    # Supplement section 2: C(t) = M_b(t - dt_CF) / A
    Ccecum <- cecum / v_cecum
    Ccecum ~ add(addSd_Ccecum)
  })
}
