Li_2026_PF06804103 <- function() {
  description <- paste(
    "Integrated longitudinal PK / repeated time-to-event peripheral neuropathy",
    "/ tumor-dynamics model for PF-06804103, a site-specific anti-HER2",
    "antibody-drug conjugate, in adults with HER2-expressing breast or gastric",
    "cancer (Li 2026). Two-compartment linear PK with IV infusion, body weight",
    "and sex on the central volume. Grade >= 2 peripheral neuropathy (PN)",
    "follows a Weibull hazard whose scale is a linear function of the",
    "concentration in an effect compartment reached through one transit",
    "compartment (shared rate KTR); the hazard is driven by the concentration",
    "in ng/mL. Tumor size (sum of longest diameters, mm) grows exponentially",
    "(rate KG, scaled by body weight and log baseline tumor size) and is",
    "killed at a rate proportional to the ADC concentration in ug/mL with",
    "an exponentially decaying resistance term (Claret form). The three",
    "sub-models were fitted sequentially (the PD fits used the PK empirical",
    "Bayes estimates) and integrated by the authors into one adaptive-dosing",
    "simulation framework, reproduced here as one model.",
    sep = " "
  )
  reference <- paste(
    "Li J, Zhou L, Leung ACF, Yin D, Chen X. Adaptive Dosing Trial Simulation",
    "Framework to Inform Dose Optimization: Case Study Using Longitudinal PKPD",
    "Models With PF-06804103 Phase 1 Data in Participants With Human Epidermal",
    "Growth Factor Receptor 2 (HER2) Positive and Low Solid Tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70261.",
    "doi:10.1002/psp4.70261. PMID 42118676.",
    "Final estimates from Tables 1-3; model structure, covariate reference",
    "values and ODEs from the NONMEM simulation code and equations in",
    "Supporting Information S1/S2.",
    sep = " "
  )
  vignette <- "Li_2026_PF06804103"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    central = list(analyte = "PF-06804103", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "PF-06804103", units = "mg", specimen = "serum", verified = TRUE),
    transit1 = list(
      analyte = "PF-06804103 (neuropathy delay chain, concentration-scaled, ng/mL)",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect = list(
      analyte = "PF-06804103 (neuropathy effect-compartment concentration, ng/mL)",
      units = "ng/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative hazard of grade >= 2 peripheral neuropathy",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_size = list(
      analyte = "tumor size (sum of longest diameters of target lesions)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed baseline. Power covariate on the central volume",
        "(Table 1 'ThetaSex-BWT' = -0.6198; see vignette for the reading of",
        "the label) and on the tumor growth rate KG (Table 3 'ThetaKG-BWT' =",
        "5.097). Reference weight 63.45 kg is the value the authors' NONMEM",
        "code (Supporting Information S1) uses for the KG covariate; Li 2026",
        "does not print a reference weight for the PK model, so the same",
        "63.45 kg is used there."
      ),
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = paste(
        "Time-fixed. Li 2026 Table 1 reports a male effect on V1",
        "('ThetaMale-V1' = -0.3186) under the paper's linear model for",
        "categorical covariates: vc = theta_vc * (1 + e_male_vc * (1 - SEXF)),",
        "so females are the reference."
      ),
      source_name = "SEX"
    ),
    TUM_SLD = list(
      description = "Baseline sum of longest diameters of target lesions (RECIST 1.1)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Initial condition of the tumor-size state",
        "(NONMEM A_0(7) = BSLD) and, after a natural-log transform, a power",
        "covariate on KG: (log(TUM_SLD) / 3.98)^e_tum_sld_p (NONMEM",
        "LBSLD / 3.98; the dummy simulation dataset has BSLD = 30 and",
        "LBSLD = 3.4 = log(30)). The reference 3.98 corresponds to",
        "exp(3.98) = 53.5 mm."
      ),
      source_name = "BSLD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 95L,
    n_studies = 1L,
    n_pk_records = 1508L,
    n_pn_events = 35L,
    n_tumor_records = 335L,
    weight_median = "63.45 kg (reference weight in the authors' code; inferred to be the analysis median)",
    disease_state = "HER2-expressing (HER2-positive and HER2-low) advanced breast or gastric cancer",
    dose_range = "PF-06804103 0.15-5 mg/kg IV on Day 1 of every 3-week cycle",
    regions = "Multi-regional phase 1 (NCT03284723)",
    notes = paste(
      "Phase 1 dose-escalation / expansion study NCT03284723. 95 participants",
      "contributed 1508 PK records (popPK), 35 grade >= 2 PN events (RTTE)",
      "and 335 tumor-size measurements (tumor dynamics). Li 2026 does not",
      "tabulate baseline demographics; see Meric-Bernstam 2023 (Mol Cancer",
      "Ther 22:1191) for the study population."
    )
  )

  ini({
    # ---- Population PK (Li 2026 Table 1) ----
    lcl <- log(0.03089); label("Clearance (L/h)") # Table 1 'ThetaCL (L/h)' = 0.03089 (RSE 9.211%)
    lvc <- log(1.681); label("Central volume of distribution for a 63.45-kg female (L)") # Table 1 'ThetaVc (L)' = 1.681 (RSE 9.406%)
    lvp <- log(3.005); label("Peripheral volume of distribution (L)") # Table 1 'ThetaVp (L)' = 3.005 (RSE 9.418%)
    lq <- log(24.86); label("Intercompartmental clearance (L/h)") # Table 1 'ThetaQ (L/h)' = 24.86 (RSE 21.5%)
    e_male_vc <- -0.3186; label("Fractional change in Vc for males, linear categorical effect (unitless)") # Table 1 'ThetaMale-V1' = -0.3186 (RSE 16.36%)
    e_wt_vc <- -0.6198; label("Power exponent of body weight (WT/63.45) on Vc (unitless)") # Table 1 'ThetaSex-BWT' = -0.6198 (RSE 18.24%)

    etalcl ~ 0.06063 # Table 1 'omega2 CL' = 0.06063 (CV 24.62%)
    etalvc ~ 0.01875 # Table 1 'omega2 Vc' = 0.01875 (CV 13.69%)
    etalvp ~ 0.009123 # Table 1 'omega2 Vp' = 0.009123 (CV 9.551%)
    etalq ~ 0.07759 # Table 1 'omega2 Q' = 0.07759 (CV 27.85%)

    expSd <- 0.3589; label("Residual error SD of log-transformed PF-06804103 concentration (log scale)") # Table 1 'ThetaRes. Error' = 0.3589 (RSE 6.126%); footnote: SD on the log scale

    # ---- Repeated time-to-event grade >= 2 peripheral neuropathy (Li 2026 Table 2) ----
    lktr <- -7.47; label("Log transit rate constant into and out of the neuropathy transit and effect compartments (log 1/h)") # Table 2 'Thetalog(KTR) (hr-1)' = -7.47 (RSE 12%)
    lalpha_haz <- -15.2; label("Log slope of the Weibull scale on effect-compartment concentration (log mL/ng/h)") # Table 2 'Thetalog(alpha)' = -15.2 (RSE 8%)
    lgam_haz <- -0.959; label("Log Weibull shape parameter (log unitless)") # Table 2 'Thetalog(gamma)' = -0.959 (RSE 10%)

    etalalpha_haz ~ 0.0205 # Table 2 'omega2 log(alpha)' = 0.0205 (CV 14%)

    # ---- Tumor dynamics (Li 2026 Table 3) ----
    lp <- log(1.958e-5); label("Tumor growth rate KG for a 63.45-kg patient with baseline SLD exp(3.98) mm (1/h)") # Table 3 'ThetaKG (hr-1)' = 0.00001958 (RSE 54%)
    lkkill <- log(4.079e-5); label("Tumor kill rate constant KD per unit concentration (1/h per ug/mL)") # Table 3 'ThetaKD (1/[hr * ug/mL])' = 0.00004079 (RSE 37%)
    lres <- log(0.0015); label("Resistance rate constant lambda, exponential decay of the kill rate (1/h)") # Table 3 'Thetalambda (hr-1)' = 0.0015 (RSE 77%)
    e_wt_p <- 5.097; label("Power exponent of body weight (WT/63.45) on KG (unitless)") # Table 3 'ThetaKG-BWT' = 5.097 (RSE 22%)
    e_tum_sld_p <- 4.536; label("Power exponent of log baseline SLD (log(TUM_SLD)/3.98) on KG (unitless)") # Table 3 'ThetaKG-BSLD' = 4.536 (RSE 38%)

    etalp ~ 2.04 # Table 3 'omega2 KG' = 2.04 (CV 143%)
    etalkkill ~ 0.894 # Table 3 'omega2 KD' = 0.894 (CV 95%)
    etalres ~ 1.31 # Table 3 'omega2 lambda' = 1.31 (CV 114%)

    expSd_tumor_size <- 0.132; label("Residual error SD of log-transformed tumor size (log scale)") # Table 3 'ThetaRes. Err' = 0.132 (RSE 20%); footnote: SD on the log scale
  })

  model({
    # ---- Population PK ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (1 + e_male_vc * (1 - SEXF)) * (WT / 63.45)^e_wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL
    Cc <- central / vc

    # ---- Grade >= 2 peripheral neuropathy (RTTE) ----
    # The hazard is driven by the concentration in ng/mL (the units of the
    # authors' simulation code, Supporting Information S1), i.e. 1000 * Cc.
    ktr <- exp(lktr)
    alpha_haz <- exp(lalpha_haz + etalalpha_haz)
    gam_haz <- exp(lgam_haz)
    cp_ngml <- 1000 * Cc

    d/dt(transit1) <- ktr * cp_ngml - ktr * transit1
    d/dt(effect) <- ktr * transit1 - ktr * effect

    # Weibull hazard (Supporting Information S2 eq 3-4); the 1e-6 offsets
    # and the zero hazard at t = 0 follow the authors' $DES block.
    lam_haz <- alpha_haz * effect
    hazard <- 0
    if (t > 0) {
      hazard <- gam_haz * (lam_haz + 1e-6)^gam_haz * (t + 1e-6)^(gam_haz - 1)
    }
    d/dt(cumhaz) <- hazard
    surv <- exp(-cumhaz)

    # ---- Tumor dynamics ----
    # KG covariate model from the authors' code: (BWT/63.45)^e * (LBSLD/3.98)^e,
    # with LBSLD = log(BSLD in mm). Time in the resistance term is time since
    # the first dose.
    p <- exp(lp + etalp) * (WT / 63.45)^e_wt_p * (log(TUM_SLD) / 3.98)^e_tum_sld_p
    kkill <- exp(lkkill + etalkkill)
    res <- exp(lres + etalres)

    tumor_size(0) <- TUM_SLD
    d/dt(tumor_size) <- p * tumor_size - kkill * exp(-res * t) * Cc * tumor_size

    Cc ~ lnorm(expSd)
    tumor_size ~ lnorm(expSd_tumor_size)
  })
}
