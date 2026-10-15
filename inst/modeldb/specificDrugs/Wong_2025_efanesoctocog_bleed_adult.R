Wong_2025_efanesoctocog_bleed_adult <- function() {
  description <- "Repeated-time-to-event (RTTE) model for treated bleeding episodes in adults and adolescents (>=12 years) with severe hemophilia A receiving efanesoctocog alfa, from Wong 2025 (XTEND-1). Weibull base hazard with a power-model FVIII-activity effect and a log-normal individual random effect on the hazard: h(t) = lam * gam * (lam * t)^(gam - 1) * exp(beta * C^alpha) * exp(eta). Treatment strategy enters the base-hazard constant lam: on-demand treatment carries an exp(2.21) = 9.1-fold higher lam than once-weekly prophylaxis. FVIII activity C (IU/dL) is driven by the embedded one-compartment PopPK model of the same paper. The cumulative hazard is an ODE state and the survivor function is exposed. Sister model files from the same paper: modellib('Wong_2025_efanesoctocog') (PopPK) and modellib('Wong_2025_efanesoctocog_bleed_pediatric') (children <12 years)."
  reference <- paste(
    "Wong N, Bhagunde P, Nyberg J, Katragadda S, Demissie M, Willemze A,",
    "Benson C, Macha S. Efanesoctocog Alfa Population Pharmacokinetics and",
    "Repeated Time-To-Event Analysis of Bleeds in Adults, Adolescents, and",
    "Children with Severe Hemophilia A. J Clin Pharmacol. 2025;65(7):860-872.",
    "doi:10.1002/jcph.70008. Open Access under CC BY-NC.",
    "The hazard equation is in the main text (Results, Adult/Adolescent Model",
    "Development); the RTTE estimates are in Table S6 and the embedded PopPK",
    "estimates in Table S4 of the supplement (Wiley file JCPH-65-860-s001.docx).",
    "Sister model files from the same paper:",
    "modellib('Wong_2025_efanesoctocog') and",
    "modellib('Wong_2025_efanesoctocog_bleed_pediatric').",
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
    ),
    REGIMEN_ONDEMAND = list(
      description = "Treatment-strategy indicator: 1 = on-demand (episodic, bleed-triggered) dosing, 0 = scheduled once-weekly prophylaxis.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (once-weekly prophylaxis)",
      notes = "Exponential effect on the base-hazard constant lam: lam_i = lam * exp(2.21 * REGIMEN_ONDEMAND) (Table S6 'Arm effect: on-demand' = 2.21, RSE 14.3%; exp(2.21) = 9.12 reproduces the Results statement that 'the hazard (lambda) for a bleeding episode was ~9-fold higher for on-demand treatment than for prophylaxis'). TIME-VARYING for XTEND-1 Arm B, who received on-demand treatment for 26 weeks and then switched to once-weekly prophylaxis for 26 weeks: of their 276 treated bleeds, 268 occurred on demand and 8 on prophylaxis (Results, Adult/Adolescent Dataset). Arm A (133 patients) was on prophylaxis throughout. Both strategies receive active efanesoctocog alfa, so the FVIII-activity effect still acts in on-demand periods through the doses given to treat bleeds.",
      source_name = "Arm / treatment strategy"
    )
  )

  covariatesDataExcluded <- list(
    AGE_FIRST_PROPH = list(
      description = "Age at first prophylaxis.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on the RTTE parameters but not retained (Methods; Results 'No other covariates were identified as significant'). Median 1 year (range 0-69) overall (Table S5)."
    ),
    N_TARGET_JOINTS = list(
      description = "Number of target joints at baseline.",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = "Screened but not retained. Median 0 (range 0-8) overall; mean 0.60 in Arm A vs 2.88 in Arm B (Table S5)."
    ),
    HJHS = list(
      description = "Baseline Hemophilia Joint Health Score (0-124; lower is better).",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened but not retained. Median 14 (range 0-98) overall (Table S5)."
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
    n_subjects = 159L,
    n_studies = 1L,
    age_range = "12-72 years (median 35)",
    weight_range = "33.9-133.0 kg (median 78.0)",
    sex_female_pct = NULL,
    race_ethnicity = NULL,
    disease_state = "Previously treated adults and adolescents (>=12 years) with severe hemophilia A (<1 IU/dL endogenous FVIII activity).",
    dose_range = "Intravenous efanesoctocog alfa 50 IU/kg once weekly for 52 weeks (Arm A, N = 133), or 50 IU/kg on demand for 26 weeks followed by 50 IU/kg once weekly for 26 weeks (Arm B, N = 26).",
    regions = "multinational",
    notes = "XTEND-1 (NCT04161495): 361 treated bleeding episodes after the start of efanesoctocog alfa among 159 patients, 85 in the prophylaxis arm and 276 in Arm B (268 during on-demand treatment, 8 on prophylaxis) (Results; Table S5). Arm B patients had more treated bleeds in the 12 months before the study (median 30 vs 1). Individual FVIII activity in the RTTE dataset was predicted with the PopPK model (Methods, Adult/Adolescent Model Development). Estimated with the Laplacian method and the likelihood option in NONMEM."
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
    # RTTE -- Table S6, adult/adolescent model column.
    # ------------------------------------------------------------------
    llam_haz <- log(0.00247); label("Constant (scale) of the Weibull base hazard, lambda, under prophylaxis (1/h)") # Table S6 'Constant of the base hazard (lambda)' = 0.00247 1/h (RSE 36.2%)
    lgam_haz <- log(1.02); label("Shape of the Weibull base hazard, gamma (unitless)") # Table S6 'Shape parameter of the base hazard (gamma)' = 1.02 (RSE 6.31%)
    beta_fviii <- -1.73; label("Constant of the FVIII-activity effect on the hazard, beta (dL/IU)") # Table S6 'Constant of the FVIII effect (beta)' = -1.73 dL/IU (RSE 14.5%)
    lalpha_fviii <- log(0.247); label("Power exponent of FVIII activity in the hazard effect, alpha (unitless)") # Table S6 'Shape parameter of FVIII effect (alpha)' = 0.247 (RSE 11.2%)
    e_regimen_ondemand_lam <- 2.21; label("Log fold change in lambda for on-demand treatment vs once-weekly prophylaxis (unitless)") # Table S6 'Arm effect: on-demand' = 2.21 (RSE 14.3%); exp(2.21) = 9.12, the Results '~9-fold'

    # Table S6 'IIV base hazard' = 114 CV% (RSE 126%), converted with the
    # same log-normal convention the paper's own Table S4 uses for its %CV
    # column: omega^2 = log(1 + 1.14^2) = 0.8327. The random effect
    # multiplies the whole hazard, exp(eta), as in the printed equation.
    etalhaz ~ 0.8327 # Table S6 'IIV base hazard' = 114 CV%; log(1 + 1.14^2) = 0.8327

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
    lam <- exp(llam_haz + e_regimen_ondemand_lam * REGIMEN_ONDEMAND)
    gam <- exp(lgam_haz)
    alpha_fviii <- exp(lalpha_fviii)

    # Time since the start of efanesoctocog alfa treatment. The floor keeps
    # (lam * t)^(gam - 1) finite at t = 0 if gamma is re-estimated below 1;
    # the hazard integral it drops over [0, 1e-4] h is below 1e-6.
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
