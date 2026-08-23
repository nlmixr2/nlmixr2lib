Okada_2024_triazolam <- function() {
  description <- paste(
    "One-compartment first-order-absorption population PK model for oral",
    "triazolam in young and elderly adults (Okada 2024, Table 3), coupled to an",
    "effect-compartment linear direct-effect PK/PD layer for two parallel",
    "endpoints: sedation (visual analogue scale, mm) and cognitive function",
    "(percent change). Age enters apparent volume and apparent clearance as a",
    "power law centred at 30 years (Eq 1), and a continuous CYP3A",
    "drug-interaction AUC ratio (AUCR) divides apparent clearance, so the model",
    "reproduces the paper's benefit-risk simulations of triazolam dose in the",
    "elderly with and without CYP3A inhibitors. The model was fit to digitised",
    "mean young-adult and elderly concentration/effect profiles from Greenblatt",
    "1991 (N Engl J Med 324:1691-8); no inter-individual variability was",
    "estimated."
  )
  reference <- paste(
    "Okada A, Sera S, Nagai N. Appropriate use of triazolam in elderly patients",
    "considering a quantitative benefit-risk assessment based on the",
    "pharmacokinetic-pharmacodynamic modeling and simulation approach supported",
    "by real-world data. BMC Pharmacol Toxicol. 2024;25:60.",
    "doi:10.1186/s40360-024-00777-z"
  )
  vignette <- "Okada_2024_triazolam"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters both Vd/F and CL/F as a power law centred at 30 years",
        "(Okada 2024 Eq 1: parameter = theta * (age/30)^theta_cov).",
        "Only two ages informed the fit (30 years young, 69 years elderly),",
        "so the paper's own Limitations section states that predictions at",
        "other ages are extrapolations."
      ),
      source_name        = "age"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Okada 2024 Table 3 reports theta3 in L/hr/kg, i.e. CL/F is the",
        "per-kilogram value multiplied by body weight (a linear per-kg scaling,",
        "not an estimated allometric exponent). Vd/F (theta2 = 119 L) is",
        "absolute and is NOT weight-scaled. Weights used in the paper's",
        "simulations were 72 kg (young, 30 years) and 69 kg (elderly, 69 years)."
      ),
      source_name        = "weight"
    ),
    AUCR = list(
      description        = paste(
        "Fold increase in triazolam AUC caused by a co-administered",
        "CYP3A-inhibiting drug"
      ),
      units              = "(unitless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Okada 2024 Table 1 footnote: AUCR = triazolam AUC with precipitant /",
        "triazolam AUC without precipitant. Reference value 1.0 = no",
        "interacting drug. Enters CL/F as a divisor (Table 3:",
        "CL = theta3 * (age/30)^theta5 * 1/AUCR), so a larger AUCR lowers",
        "clearance and prolongs exposure proportionally. Table 1 tabulates the",
        "15 literature values used to drive the paper's simulations: ritonavir",
        "40.7, itraconazole 27.1, ketoconazole 22.4, mibefradil 8.36,",
        "clarithromycin 5.26, fluconazole 4.42, nefazodone 3.90, erythromycin",
        "3.80, troleandomycin 3.76, diltiazem 3.38, grapefruit juice 2.43,",
        "cimetidine 2.20, isoniazid 1.46, hormonal contraceptives 1.44,",
        "ranitidine 1.31."
      ),
      source_name        = "AUCR"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "triazolam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "triazolam", units = "mg", specimen = "plasma", verified = TRUE),
    effect  = list(analyte = "triazolam", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    age_range      = "30 and 69 years (two mean profiles only)",
    age_median     = NA_character_,
    weight_range   = "69 and 72 kg (two mean profiles only)",
    weight_median  = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Healthy young and elderly volunteers.",
    dose_range     = "0.25 mg oral triazolam, single dose (fitted); 0.0625-0.25 mg simulated",
    regions        = "United States (source PK/PD data); model built in Japan",
    notes          = paste(
      "Okada 2024 did not collect new subjects. Mean plasma triazolam",
      "concentration, sedation (visual analogue scale) and cognitive-function",
      "time courses were digitised with WebPlotDigitizer 4.3 from Greenblatt",
      "DJ, Harmatz JS, Shapiro L, Engelhardt N, Gouthro TA, Shader RI.",
      "'Sensitivity to triazolam in the elderly.' N Engl J Med.",
      "1991;324(24):1691-8 (reference 18). That source supplied only group",
      "means for a young cohort (30 years, 72 kg) and an elderly cohort",
      "(69 years, 69 kg); Okada 2024 does not restate the per-cohort subject",
      "counts, so n_subjects is left unset rather than inferred. Because only",
      "group means were fitted, inter-individual variability was deliberately",
      "not estimated (Methods, 'Pharmacokinetic-pharmacodynamic modeling').",
      "Estimation used Certara NLME 8.1. The paper's pharmacovigilance half",
      "(JADER, FAERS and NDB claims analyses) informed the clinical framing",
      "but not the model, and is not represented here."
    )
  )

  ini({
    # =================================================================
    # PK: one compartment, first-order oral absorption.
    # Time unit: hour. Concentration unit: ng/mL. Dose unit: mg.
    # No inter-individual variability: Okada 2024 Methods states
    # "inter-individual variability was not estimated" because only mean
    # profiles were available for fitting.
    # =================================================================
    lka <- log(3.16)
    label("Absorption rate constant (Ka, 1/h)")                                      # Table 3: Ka = theta1 = 3.16 1/hr (CV 5.27%)

    lvc <- log(119)
    label("Apparent volume of distribution at age 30 years (Vd/F, L)")               # Table 3: Vd = theta2 * (age/30)^theta4; theta2 = 119 L (CV 1.86%)

    e_age_vc <- -0.49
    label("Power-law exponent for age on apparent volume of distribution (unitless)") # Table 3: theta4 = -0.49 (CV 6.53%); form from Eq 1

    lcl <- log(0.47)
    label("Weight-normalised apparent clearance at age 30 years (CL/F, L/h/kg)")     # Table 3: CL = theta3 * (age/30)^theta5 * 1/AUCR; theta3 = 0.47 L/hr/kg (CV 1.92%)

    e_age_cl <- -0.53
    label("Power-law exponent for age on apparent clearance (unitless)")             # Table 3: theta5 = -0.53 (CV 4.30%); form from Eq 1

    propSd <- 0.0620
    label("Proportional residual error on plasma triazolam concentration (fraction)") # Table 3: sigma = 6.20% (CV 10.8%); proportional form from Eq 2

    # =================================================================
    # PK/PD: one effect compartment shared by two linear direct-effect
    # endpoints. An anticlockwise hysteresis loop between plasma
    # concentration and effect motivated the effect compartment (Results:
    # OFV 193.7 for the effect-compartment model vs 223.5 direct-response
    # linear and 225.2 direct-response Emax). Both endpoints are
    # change-from-baseline measures and the paper reports no intercept.
    # =================================================================
    lke0 <- log(3.42)
    label("Plasma-to-effect-compartment equilibration rate constant (Keo, 1/h)")     # Table 3: Keo = 3.42 1/hr (CV 12.6%)

    slope_sedation <- 12.1
    label("Slope of sedation VAS on effect-site triazolam concentration (mm per ng/mL)") # Table 3 footnote: "Sedation = 12.1 x concentration" (CV 1.05%)

    slope_cognitive <- -11.0
    label("Slope of cognitive function on effect-site triazolam concentration (% per ng/mL)") # Table 3 footnote: "Cognitive function = -11.0 x concentration" (CV 3.14%)

    addSd_sedation <- 6.81
    label("Additive residual error on sedation VAS (mm)")                            # Table 3: sigma Sedation = 6.81 mm (CV 21.8%); additive form from Eq 3

    addSd_cognitive <- 3.27
    label("Additive residual error on cognitive function (%)")                       # Table 3: sigma Cognitive function = 3.27 % (CV 14.5%); additive form from Eq 3
  })

  model({
    ka <- exp(lka)

    # Okada 2024 Eq 1 power law, centred at 30 years, instantiated by Table 3
    # as Vd = theta2 * (age/30)^theta4. Vd/F is absolute (L), not per kg.
    vc <- exp(lvc) * (AGE / 30)^e_age_vc

    # Table 3: CL = theta3 * (age/30)^theta5 * 1/AUCR. theta3 carries units of
    # L/hr/kg, so the per-kilogram value is multiplied by body weight. AUCR is
    # the fold increase in triazolam AUC produced by a co-administered CYP3A
    # inhibitor (Table 1 footnote); AUCR = 1 means no interacting drug.
    cl <- exp(lcl) * WT * (AGE / 30)^e_age_cl / AUCR

    kel <- cl / vc
    ke0 <- exp(lke0)

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # central (mg) / vc (L) = mg/L; 1 mg/L = 1000 ng/mL.
    # The state expression is written inline inside d/dt(effect) rather than
    # routed through the named Cc, because a named intermediate referencing an
    # ODE state can silently evaluate to zero inside d/dt() in the UI form.
    d/dt(effect) <- ke0 * (central / vc * 1000 - effect)

    Cc <- central / vc * 1000
    Cc ~ prop(propSd)

    # Both PD endpoints are linear in effect-site concentration with no
    # intercept: they are change-from-baseline measures (Table 3 footnote).
    sedation <- slope_sedation * effect
    sedation ~ add(addSd_sedation)

    cognitive <- slope_cognitive * effect
    cognitive ~ add(addSd_cognitive)
  })
}
