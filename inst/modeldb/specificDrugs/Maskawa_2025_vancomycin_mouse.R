Maskawa_2025_vancomycin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse).",
    "Gastrointestinal transit model for fecal concentrations of orally",
    "administered, non-absorbed vancomycin in male C57BL/6N mice (the",
    "paper's zero-order Model 1). The dose empties from the stomach into",
    "the cecum at a constant rate over a duration t_in, is eliminated from",
    "the cecum by first-order transit (k_e), and the fecal concentration is",
    "the cecal amount divided by the mass of cecal contents A. The",
    "stomach-to-cecum and cecum-to-feces transit delays are only",
    "identifiable as their sum t_tot, which is applied as a lag on the",
    "cecal input, so the cecum state holds the cecal amount shifted",
    "forward to the time it appears in feces. Doses go into cecum with",
    "rate = -2 so that dur(cecum) supplies the emptying duration. Fitted",
    "per mouse by penalised least squares; typical values are the",
    "arithmetic means over the six treated mice and the between-mouse",
    "variances are derived from the reported coefficients of variation.",
    sep = " "
  )
  reference <- paste(
    "Maskawa R, Takayasu L, Takayasu H, Watanabe K, Takemine S, Kakimoto T,",
    "Takeshita K, Narushima S, Suda W, Takayasu M.",
    "High-resolution fecal pharmacokinetic modeling in mice with orally",
    "administered antibiotics.",
    "Sci Rep. 2025;15:24441.",
    "doi:10.1038/s41598-025-10293-7.",
    sep = " "
  )
  vignette <- "Maskawa_2025_vancomycin_mouse"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/g"
  )

  compartmentData <- list(
    # Methods 'Mathematical model': M_b(t), antibiotic mass in the cecum.
    # The fecal concentration is C(t) = M_b(t - dt_CF) / A, so with the
    # lag set to t_tot this state is the cecal mass time-shifted to the
    # moment it is expelled, and cecum / A is the measured fecal
    # concentration (wet-weight basis, LC-MS/MS).
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
      "Three treated mice per dose group (plus three water-only controls",
      "not used in the model fit). Food (CL-2) and water ad libitum.",
      "Feces collected every 4 h from day -3 to day 3, at 12:00 and 16:00",
      "on day 4, then every 1-7 days to day 36; vancomycin quantified by",
      "solid-phase extraction and LC-MS/MS and normalised to fecal wet",
      "weight (lower limit of quantification 0.01 ug/mL in extract).",
      "Body weight not reported. Methods 'Mice and fecal sample",
      "collection'; Fig. 1."
    )
  )

  ini({
    # Typical values: Table 2 'Mean (CV%)' row, arithmetic mean over the
    # six treated mice (L1-L3, H1-H3), each fitted individually.
    ltlag <- log(5.4); label("Total gastrointestinal transit time of the antibiotic t_tot (h)") # Table 2, t_tot mean 5.4 h (CV 54%)
    ld1 <- log(12); label("Duration of zero-order stomach-to-cecum emptying t_in (h)") # Table 2, t_in mean 12 h (CV 67%)
    lkfec <- log(0.41); label("First-order elimination rate constant from the cecum to feces k_e (1/h)") # Table 2, k_e mean 0.41 1/h (CV 48%)
    lv_cecum <- log(0.31); label("Mass of cecal contents A (g)") # Table 2, A mean 0.31 g (CV 65%)

    # Between-mouse variability: omega^2 = log(1 + CV^2) from the Table 2
    # CV% of the six individual estimates (CV computed with the n
    # denominator, as numpy does).
    etaltlag ~ 0.2559 # Table 2, t_tot CV 54%: log(1 + 0.54^2)
    etald1 ~ 0.3708 # Table 2, t_in CV 67%: log(1 + 0.67^2)
    etalkfec ~ 0.2073 # Table 2, k_e CV 48%: log(1 + 0.48^2)
    etalv_cecum ~ 0.3524 # Table 2, A CV 65%: log(1 + 0.65^2)

    # Unweighted least-squares fit (Methods Eq. 8, lambda = 0.9), so the
    # residual error is additive; its SD is the pooled RMSE of Model 1.
    addSd_Ccecum <- 0.21; label("Additive residual error on fecal concentration (mg/g)") # Table S2, Model1 RMSE 'All' 2.1e-1
  })
  model({
    tlag <- exp(ltlag + etaltlag)
    d1 <- exp(ld1 + etald1)
    kfec <- exp(lkfec + etalkfec)
    v_cecum <- exp(lv_cecum + etalv_cecum)

    # Methods Eq. 3: dM_b/dt = q / t_in - k_e M_b while the stomach empties,
    # then -k_e M_b; the q / t_in input is the zero-order dose infusion.
    d/dt(cecum) <- -kfec * cecum
    alag(cecum) <- tlag
    dur(cecum) <- d1

    # Methods Eq. 5-7: C(t) = M_b(t - dt_CF) / A
    Ccecum <- cecum / v_cecum
    Ccecum ~ add(addSd_Ccecum)
  })
}
