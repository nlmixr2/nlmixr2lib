Altcheh_2023_benznidazole <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination for oral benznidazole in neonates, infants",
    "and children (0-12 years) with Chagas disease treated with a new",
    "12.5 mg dispersible pediatric formulation (Altcheh 2023; PEDCHAGAS",
    "network, NCT01549236, n = 81). Apparent clearance CL/F = 2.1 L/h at",
    "the 70 kg allometric reference, scaled by (WT/70)^0.75 with the",
    "exponent fixed at the theory-based 3/4. Apparent volume of",
    "distribution V/F is a LINEAR function of body weight, V/F = 1.75 +",
    "0.804 * WT (L), not a power function. Absorption rate constant Ka =",
    "1.22 1/h. Inter-individual variability is on Ka (61.9% CV), V/F",
    "(47.4% CV) and CL/F (38.3% CV); residual error is proportional",
    "(37.1% CV). Weight was the only covariate retained; age, sex,",
    "postmenstrual age, dose, tablet formulation and study center were",
    "screened and not retained."
  )
  reference <- paste(
    "Altcheh J, Moscatelli G, Caruso M, Moroni S, Bisio M, Miranda MR,",
    "Monla C, Vaina M, Valdez M, Moran L, Ramirez T, Ledesma Patino O,",
    "Riarte A, Gonzalez N, Fernandes J, Alves F, Ribeiro I,",
    "Garcia-Bournissen F. Population pharmacokinetics of benznidazole in",
    "neonates, infants and children using a new pediatric formulation.",
    "PLoS Neglected Tropical Diseases. 2023;17(5):e0010850.",
    "doi:10.1371/journal.pntd.0010850.",
    sep = " "
  )
  vignette <- "Altcheh_2023_benznidazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Benznidazole was measured in whole dried blood
  # spots (Whatman 903 paper, LC-ESI-MS/MS, LLOQ 50 ng/mL) rather than
  # plasma (Altcheh 2023 Methods, "Measurement of benznidazole in blood
  # samples"), so the central compartment specimen is blood.
  compartmentData <- list(
    depot   = list(analyte = "benznidazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "benznidazole", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model (Altcheh 2023",
        "Results, 'Covariate modeling'). Enters CL/F allometrically as",
        "(WT/70)^0.75 with the exponent FIXED at the theory-based 3/4",
        "value, and enters V/F LINEARLY as 1.75 + 0.804 * WT -- note the",
        "V/F relationship is additive on the natural scale, NOT a power",
        "function, so V/F is not proportional to weight. Time-varying:",
        "because weight was recorded only at the start of the 60-day",
        "dosing period, intermediate weights for children under 18 months",
        "were interpolated from the WHO weight-for-age percentile curves",
        "following each child's own percentile (Altcheh 2023 Methods,",
        "'Covariate model')."
      ),
      source_name        = "Wt"
    )
  )

  # Covariates screened by Altcheh 2023 and NOT retained in the final
  # model (Methods 'Covariate model'; Results 'Covariate modeling':
  # "None of the remaining covariates ... had a significant effect on
  # the model fit"). Listed for provenance only; not referenced in
  # model(). Dose per kg, total daily benznidazole dose, study center
  # and tablet formulation (12.5 mg dispersible vs 100 mg) were also
  # screened and not retained -- see population$notes; they are not
  # given register entries here because no equation uses them.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on V/F and CL/F. Altcheh 2023 Results: 'Addition of age",
        "as a covariate for V/F did not improve the fit beyond the effect",
        "of weight.' CL/F correlates strongly with age (Fig 4) but only",
        "because age and weight are collinear in children; weight is the",
        "retained covariate. Register canonical carries years; the paper",
        "tabulates age in months (Table 1)."
      ),
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened as 'gender' and not retained. Neither Cmax nor trough",
        "concentrations differed between girls and boys (Mann-Whitney",
        "p > 0.05; Altcheh 2023 Results). Cohort was 47/81 (58.0%) female",
        "(Table 1)."
      ),
      source_name        = "gender"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on CL/F via an empirical sigmoidal PMA-CL/F relationship",
        "and not retained; the allometric weight model fit marginally",
        "better and was chosen for parsimony (Altcheh 2023 Results,",
        "'Covariate modeling'). Altcheh 2023 estimated PMA by adding 9",
        "months to the actual postnatal age in months at each blood",
        "sampling timepoint (Methods, 'Covariate model'), so this column",
        "is time-varying and carried in MONTHS, matching the register",
        "default."
      ),
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 81L,
    n_studies      = 1L,
    n_observations = 383L,
    age_range      = "0-12 years (14 newborns < 1 month of age; 41/81 under 2 years)",
    age_median     = "16 months (IQR 2.76-75.6); mean 42.3 months (SD 47.3)",
    weight_range   = "not reported as a range; median 12 kg (IQR 8.1-23)",
    weight_mean    = "16.69 kg (SD 12.8)",
    sex_female_pct = round(100 * 47 / 81, 1),
    disease_state  = paste(
      "Children 0-12 years with Chagas disease (Trypanosoma cruzi",
      "infection), treatment-naive, all asymptomatic without cardiac",
      "symptoms or other organ involvement. Diagnosis by at least two",
      "positive serologic tests (age > 8 months) or by direct",
      "microhematocrit parasite detection (age < 8 months). 69/80 (85%)",
      "were qPCR-positive for T. cruzi DNA at enrolment. Excluded:",
      "pregnancy, investigational drug within the prior month,",
      "cardiovascular / hepatic / neurologic / endocrine or other major",
      "systemic disease, immunocompromise."
    ),
    dose_range     = paste(
      "Oral benznidazole 7.5 mg/kg/day divided into two daily doses for",
      "60 days (actual administered dose mean 7.35 mg/kg, SD 0.65;",
      "Table 1). Formulation assigned by a 14 kg enrolment weight",
      "cut-off: 12.5 mg dispersible pediatric tablet below, 100 mg",
      "non-dispersible tablet at or above (LAFEPE, Brazil)."
    ),
    sampling       = paste(
      "Sparse design, 5 blood samples per child at random times within",
      "pre-specified windows: Day 0 at 2-5 h after the first dose; Day 7",
      "and Day 30 at steady state, each between pre-dose and 8 h",
      "post-dose; and two samples 12-24 h AFTER THE FINAL dose on Day 60",
      "(i.e. during washout, not steady-state troughs). 387 samples",
      "collected, 383 quantifiable; the 4 BLQ samples (1%) were excluded",
      "from the PK analysis after confirming they did not change the fit.",
      "Whole dried blood spots on Whatman 903 paper, LC-ESI-MS/MS,",
      "linear 50-20,000 ng/mL, LLOQ 50 ng/mL."
    ),
    regions        = paste(
      "Five pediatric centers in Argentina (PEDCHAGAS Network), all in",
      "areas with certified vector control to avoid confounding by",
      "re-infection. Enrolment May 2011 - August 2012."
    ),
    notes          = paste(
      "Open-label single-group interventional trial, ClinicalTrials.gov",
      "NCT01549236. Baseline demographics from Altcheh 2023 Table 1;",
      "final parameter estimates from Table 2 and the three displayed",
      "equations that immediately precede it. Covariates screened and NOT",
      "retained: age, sex, postmenstrual age, dose per kg, total daily",
      "benznidazole dose, study center, and tablet formulation (12.5 mg",
      "dispersible vs 100 mg). Table 1 and the Results text disagree on",
      "two demographic counts (median age 16 vs 22 months; 44 vs 46",
      "patients on the 12.5 mg formulation); Table 1 is used here because",
      "its formulation split, 44/81, is the one consistent with the 54%",
      "the text itself reports. Neither value is used by model()."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural PK parameters. Altcheh 2023 reports the final model as
    # three displayed equations between Fig 3 and Table 2:
    #   Ka (/h) = 1.22
    #   V/F (L) = 1.75 (L) + 0.804 (L/Kg) x Wt (Kg)
    #   CL/F (L/h) = 2.1 (L/h) x (Wt (Kg) / 70 (Kg))^0.75
    # with the same point estimates in Table 2. All parameters are
    # apparent (CL/F, V/F) because no intravenous formulation of
    # benznidazole exists, so F cannot be estimated (Discussion,
    # limitations paragraph).
    # ---------------------------------------------------------------
    lka <- log(1.22);  label("Absorption rate constant Ka (1/h)")                                        # Altcheh 2023 Table 2 and displayed equation: Ka = 1.22 [95% CI 0.91; 1.65]
    lcl <- log(2.1);   label("Apparent oral clearance CL/F (L/h) at the 70 kg allometric reference")     # Altcheh 2023 Table 2 and displayed equation: CL/F = 2.1 L/h at Wt = 70 kg [95% CI 1.89; 2.39]

    # V/F is LINEAR in weight: V/F = b0 + b1 * WT. lvc carries the
    # intercept b0 (the natural-scale offset, in L) and e_wt_vc carries
    # the slope b1 (L per kg). This is the additive-covariate encoding
    # established by Goggin_2004_emfilermin.R (tvp_vc <- exp(lvc) +
    # e_wt_vc * (WT - 62)); do NOT read lvc as a typical volume at a
    # reference weight, and do NOT exponentiate the weight term.
    lvc     <- log(1.75); label("Intercept b0 of the linear V/F-weight relationship (L, apparent volume at extrapolated WT = 0)")  # Altcheh 2023 Table 2: b0 = 1.75 [95% CI 0; 4.48]
    e_wt_vc <- 0.804;     label("Linear additive effect of WT on V/F (L per kg body weight)")                                      # Altcheh 2023 Table 2: b1 = 0.804 [95% CI 0.62; 0.99]

    # Allometric exponent on CL/F. Altcheh 2023 Methods ('Covariate
    # model') states the allometric scaling used "the theory-based
    # allometric exponent of 3/4", and Results confirms the final model
    # used an "empirical allometric influence of weight on CL/F
    # consistent with theory based allometry". The exponent is a
    # structural assumption, not an estimate: Table 2 reports a CI for
    # the 2.1 L/h coefficient but none for the 0.75 exponent, which is
    # printed as a literal in the displayed CL/F equation.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT on CL/F (unitless, the theory-based 3/4)")  # Altcheh 2023 Methods 'Covariate model' + displayed CL/F equation exponent 0.75

    # ---------------------------------------------------------------
    # Inter-individual variability. Altcheh 2023 Methods ('Population PK
    # analysis', structural model paragraph): "Inter-individual
    # variability in PK parameters was described using exponential error
    # models", i.e. P_i = TVP * exp(eta_i), and Table 2's
    # "Interindividual variability" column reports each as a %CV. The
    # internal log-scale variance is omega^2 = log(CV^2 + 1).
    # ---------------------------------------------------------------
    etalka ~ 0.324371  # Altcheh 2023 Table 2: IIV Ka   = 61.9% CV -> omega^2 = log(0.619^2 + 1) = 0.324371
    etalvc ~ 0.202676  # Altcheh 2023 Table 2: IIV V/F  = 47.4% CV -> omega^2 = log(0.474^2 + 1) = 0.202676
    etalcl ~ 0.136879  # Altcheh 2023 Table 2: IIV CL/F = 38.3% CV -> omega^2 = log(0.383^2 + 1) = 0.136879

    # ---------------------------------------------------------------
    # Residual error. Altcheh 2023 Methods: "The residual variability
    # was described using a proportional error model."
    # ---------------------------------------------------------------
    propSd <- 0.371; label("Proportional residual error on Cc (fraction CV)")  # Altcheh 2023 Table 2: Residual error (Proportional) = 37.1% [95% CI 34.6; 41.7%]
  })

  model({
    # Individual PK parameters. Ka and CL/F carry log-normal random
    # effects on the log scale; V/F carries its log-normal random effect
    # multiplicatively on the natural-scale linear weight model.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl

    # V/F = (b0 + b1 * WT) * exp(eta). The weight term is additive on
    # the natural scale (Altcheh 2023 displayed V/F equation), so the
    # bracket is formed before the exponential IIV multiplier is applied.
    vc <- (exp(lvc) + e_wt_vc * WT) * exp(etalvc)

    kel <- cl / vc

    # One-compartment system with first-order absorption from the depot
    # and first-order elimination from the central compartment. Altcheh
    # 2023 Results ('BNZ PK modeling and simulation'): two- and
    # three-compartment models failed to decrease the objective
    # function, so a one-compartment model with first-order oral
    # absorption was adopted.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration units: dose mg / volume L = mg/L, matching the
    # LC-ESI-MS/MS dried-blood-spot assay (LLOQ 0.05 mg/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
