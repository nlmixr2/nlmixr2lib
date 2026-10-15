Wong_2025_efanesoctocog_bleed_pediatric <- function() {
  description <- "Repeated-time-to-event (RTTE) model for treated bleeding episodes in children (<12 years) with severe hemophilia A receiving once-weekly efanesoctocog alfa prophylaxis, from Wong 2025 (XTEND-Kids). The adult/adolescent model structure re-estimated on pediatric data only: Weibull base hazard with a power-model FVIII-activity effect and a log-normal individual random effect on the hazard, h(t) = lam * gam * (lam * t)^(gam - 1) * exp(beta * C^alpha) * exp(eta). Compared with adults, children have a lower base-hazard constant, a shape gamma < 1 (base hazard decreasing over time) and a shallower FVIII effect; no covariates were retained. FVIII activity C (IU/dL) is driven by the embedded one-compartment PopPK model of the same paper. The cumulative hazard is an ODE state and the survivor function is exposed. Sister model files from the same paper: modellib('Wong_2025_efanesoctocog') (PopPK) and modellib('Wong_2025_efanesoctocog_bleed_adult') (adults and adolescents)."
  reference <- paste(
    "Wong N, Bhagunde P, Nyberg J, Katragadda S, Demissie M, Willemze A,",
    "Benson C, Macha S. Efanesoctocog Alfa Population Pharmacokinetics and",
    "Repeated Time-To-Event Analysis of Bleeds in Adults, Adolescents, and",
    "Children with Severe Hemophilia A. J Clin Pharmacol. 2025;65(7):860-872.",
    "doi:10.1002/jcph.70008. Open Access under CC BY-NC.",
    "The hazard equation is in the main text (Results, Adult/Adolescent Model",
    "Development); the re-estimated pediatric RTTE estimates are in Table S6 and the embedded PopPK",
    "estimates in Table S4 of the supplement (Wiley file JCPH-65-860-s001.docx).",
    "Sister model files from the same paper:",
    "modellib('Wong_2025_efanesoctocog') and",
    "modellib('Wong_2025_efanesoctocog_bleed_adult').",
    sep = " "
  )
  vignette <- "Wong_2025_efanesoctocog"
  # The individual random effect multiplies the whole hazard, exp(eta), in
  # the printed hazard equation, not the base-hazard constant lambda inside
  # (lambda * t)^(gamma - 1); it has no 1-to-1 ini() typical value.
  paper_specific_etas <- "etalhaz"
  units <- list(time = "h", dosing = "IU", concentration = "IU/dL")

  covariateData <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Acts on the EMBEDDED PopPK model only (power effects on CL and V centred at 78.3 kg), and thereby on the hazard through FVIII activity. Not a covariate of the RTTE parameters. Time-varying in the source PopPK analysis.",
      source_name = "Body weight"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator: 1 = Asian, 0 = non-Asian.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Acts on the EMBEDDED PopPK model only (CL 10.4% lower in Asian patients). Not a covariate of the RTTE parameters.",
      source_name = "Asian"
    )
  )

  covariatesDataExcluded <- list(
    AGE_LT6 = list(
      description = "Age-group indicator: 1 = younger than 6 years, 0 = 6 to <12 years.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (6 to <12 years)",
      notes = "Screened on the re-estimated pediatric RTTE parameters but not retained ('No additional significant covariate-parameter relationships were identified using the pediatric dataset', Results). Observed time to bleed did not differ visibly between the two age groups (Figure S3B)."
    ),
    BLEED_HISTORY_HIGH = list(
      description = "Bleed-history indicator: 1 = above the median, 0 = at or below the median.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (at or below the median)",
      notes = "Screened but not retained (Methods, Pediatric Model Development; Results)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "efanesoctocog alfa (FVIII activity)", units = "IU", specimen = "plasma", verified = TRUE),
    cumhaz = list(
      analyte = "Cumulative hazard of a treated bleeding episode",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 74L,
    n_studies = 1L,
    age_range = "1.4-11 years (median 5)",
    weight_range = "11.4-66.5 kg (median 22.0)",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = "Previously treated children (<12 years) with severe hemophilia A (<1 IU/dL endogenous FVIII activity).",
    dose_range = "Intravenous efanesoctocog alfa 50 IU/kg once weekly for 52 weeks.",
    regions = "multinational",
    notes = "XTEND-Kids (NCT04759131): 63 treated bleeding episodes among 74 children, 38 aged <6 years and 36 aged 6 to <12 years (Results; Table S7). One patient contributed 21 inferred bleeding episodes because they received an intensive consolidation regimen after traumatic bleeds rather than protocol-specified weekly prophylaxis; the model was estimated with this patient included. The adult/adolescent model was first assessed against the pediatric data by VPC and its structural parameters were then re-estimated on the pediatric data alone, which improved the fit (Methods, Pediatric Model Development; Discussion)."
  )

  ini({
    # ------------------------------------------------------------------
    # Embedded PK -- identical to modellib('Wong_2025_efanesoctocog'),
    # Table S4.
    # ------------------------------------------------------------------
    lcl <- log(0.433); label("Clearance at 78.3 kg in a non-Asian patient (dL/h)") # Table S4 'Typical value of CL' = 0.433 dL/h
    lvc <- log(30.2); label("Central volume of distribution at 78.3 kg (dL)") # Table S4 'Typical value of V' = 30.2 dL
    e_wt_cl <- 0.677; label("Power exponent of body weight on CL, centred at 78.3 kg (unitless)") # Table S4 = 0.677
    e_wt_vc <- 0.791; label("Power exponent of body weight on V, centred at 78.3 kg (unitless)") # Table S4 = 0.791
    e_race_asian_cl <- log(0.896); label("Log fold change in CL for Asian vs non-Asian patients (unitless)") # Table S4 = 0.896

    # Table S4 log-scale variances and covariance; see the PopPK file.
    etalcl + etalvc ~ c(
      0.0354,
      0.0148, 0.0209
    ) # Table S4: IIV CL 0.0354; Corr-CL-V (covariance) 0.0148; IIV V 0.0209

    # ------------------------------------------------------------------
    # RTTE -- Table S6, pediatric (<12 years) re-estimated model column.
    # ------------------------------------------------------------------
    llam_haz <- log(0.0000429); label("Constant (scale) of the Weibull base hazard, lambda (1/h)") # Table S6 pediatric 'Constant of the base hazard (lambda)' = 0.0000429 1/h (RSE 58.2%)
    lgam_haz <- log(0.694); label("Shape of the Weibull base hazard, gamma (unitless)") # Table S6 pediatric 'Shape parameter of the base hazard (gamma)' = 0.694 (RSE 14.5%)
    beta_fviii <- -0.00535; label("Constant of the FVIII-activity effect on the hazard, beta (dL/IU)") # Table S6 pediatric 'Constant of the FVIII effect (beta)' = -0.00535 dL/IU (RSE 1.55%)
    lalpha_fviii <- log(1.27); label("Power exponent of FVIII activity in the hazard effect, alpha (unitless)") # Table S6 pediatric 'Shape parameter of FVIII effect (alpha)' = 1.27 (RSE 4.37%)

    # Table S6 pediatric 'IIV base hazard' = 186 CV% (RSE 260%), converted
    # with the log-normal convention of the paper's Table S4:
    # omega^2 = log(1 + 1.86^2) = 1.4951. The random effect multiplies the
    # whole hazard, exp(eta), as in the printed equation.
    etalhaz ~ 1.4951 # Table S6 pediatric 'IIV base hazard' = 186 CV%; log(1 + 1.86^2) = 1.4951

    # The source fits the RTTE likelihood directly, so there is no
    # observation-error model to translate. This placeholder additive
    # residual on the survivor function lets nlmixr2 accept the model for
    # forward simulation; it is not from the source.
    addSd <- 0.001; label("Placeholder additive residual error on the survivor-function output `sur` (unitless); not from the source")
  })
  model({
    # --- embedded PK (identical to modellib('Wong_2025_efanesoctocog')) ---
    cl <- exp(lcl + etalcl + e_race_asian_cl * RACE_ASIAN) * (WT / 78.3)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 78.3)^e_wt_vc
    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc

    # --- bleed hazard -----------------------------------------------------
    lam <- exp(llam_haz)
    gam <- exp(lgam_haz)
    alpha_fviii <- exp(lalpha_fviii)

    # Time since the start of efanesoctocog alfa treatment. With gamma < 1
    # the base hazard is infinite at t = 0 (an integrable singularity); the
    # floor keeps it finite, and the hazard integral it drops over
    # [0, 1e-4] h is below 1e-6.
    tw <- max(t, 1e-4)

    cfviii <- max(Cc, 0)
    fviii_effect <- exp(beta_fviii * cfviii^alpha_fviii)

    hazard <- lam * gam * (lam * tw)^(gam - 1) * fviii_effect * exp(etalhaz)

    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0

    # Probability of remaining free of treated bleeds since the start of
    # the integration window.
    sur <- exp(-cumhaz)

    sur ~ add(addSd)
  })
}
