Sager_2023_sotrovimab_progression <- function() {
  description <- paste0(
    "Logistic-regression exposure-response model for progression of ",
    "COVID-19 through day 29 -- hospitalization for more than 24 h for ",
    "acute management of illness due to any cause, or death -- in ",
    "non-hospitalized adults and adolescents with mild-to-moderate ",
    "COVID-19 at high risk of progression who received a single dose of ",
    "sotrovimab (Sager 2023, the COMET-TAIL study, NCT04913675; 902 ",
    "patients of whom 20 progressed). The probability of progression is ",
    "modelled on the logit scale as a linear function of the individual ",
    "serum sotrovimab concentration 168 h (7 days) after the dose, plus ",
    "an additive shift on the intercept for patients carrying more than ",
    "one protocol-defined risk factor for progression. The number of risk ",
    "factors moves only the intercept -- the model-estimated placebo ",
    "progression rate -- and has no effect on the exposure slope, so the ",
    "authors conclude that dose optimization on the basis of risk-factor ",
    "count is not indicated. The exposure metric is not an observation: ",
    "it is generated per patient from the companion population PK model, ",
    "modellib('Sager_2023_sotrovimab'), and supplied here as the ",
    "covariate column CONC_SOTRO_168H. Sotrovimab AUC through day 28 and the ",
    "96 h concentration were also significant predictors; the 168 h ",
    "concentration was carried forward because the median time from ",
    "randomization to a progression event was 5.5 days. The authors ",
    "caution that the small number of progressors, the absence of a ",
    "placebo arm in COMET-TAIL, and enrolment during a single ",
    "(delta-predominant) period limit generalization of this ",
    "relationship across SARS-CoV-2 variants."
  )
  reference <- paste(
    "Sager JE, El-Zailik A, Passarell J, Roepcke S, Li X, Aldinger M,",
    "Nader A, Skingsley A, Alexander EL, Yeh WW, Mogalian E, Garner C,",
    "Peppercorn A, Shapiro AE, Reyes M.",
    "Population pharmacokinetics and exposure-response analysis of a single",
    "dose of sotrovimab in the early treatment of patients with mild to",
    "moderate COVID-19.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(6):853-864.",
    "doi:10.1002/psp4.12958.",
    "Final exposure-response estimates are Table 3; the base-model",
    "estimates and the model-predicted progression rates by treatment arm",
    "are Tables S9 and S10 of the Supporting Information, and the",
    "NONMEM control stream is MODEL CODE S2 of the same document.",
    "Companion population PK model:",
    "modellib('Sager_2023_sotrovimab').",
    sep = " "
  )
  vignette <- "Sager_2023_sotrovimab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CONC_SOTRO_168H = list(
      description        = "Individual serum sotrovimab concentration 168 h (7 days) after a single intravenous or intramuscular dose.",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "NOT an observation. Sager 2023 (Methods, 'ER modeling') generated ",
        "it by integrating each patient's predicted concentration-time ",
        "profile from the final population PK model using that patient's ",
        "own empirical Bayes parameter estimates, then reading the profile ",
        "at 168 h. Generate it for simulation with ",
        "modellib('Sager_2023_sotrovimab'). Enters UNCENTRED, so the ",
        "intercept is the log-odds of progression at zero sotrovimab ",
        "exposure -- which is why the authors describe the intercept as ",
        "the model-estimated placebo response even though COMET-TAIL had ",
        "no placebo arm. The source column is CP168 (MODEL CODE S2 $INPUT) ",
        "and the slope is printed in units of 1/[ug/mL], fixing the ",
        "concentration scale. Distribution implied by Sager 2023 Table ",
        "S10 and reproduced in the vignette: roughly 60 ug/mL after 500 mg ",
        "i.v., 21 ug/mL after 500 mg i.m. and 11 ug/mL after 250 mg i.m."
      ),
      source_name        = "CP168"
    ),
    NRISK_GT1 = list(
      description        = "Indicator that the patient has more than one protocol-defined risk factor for progression of COVID-19; 1 = more than one, 0 = one or fewer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (one or fewer risk factors; 630 of 902 patients, 69.8%)",
      notes              = paste0(
        "The source column is RISKCATN (MODEL CODE S2 $INPUT, listed ",
        "under ';--catcovs-'), the dichotomised form of the risk-condition ",
        "count RFNUMCND. Sager 2023 Table S8 gives the split as 630 ",
        "(69.8%) with one or fewer risk factors and 272 (30.2%) with more ",
        "than one. The count is INCLUSIVE of age and body mass index ",
        "(Sager 2023 Methods, 'ER modeling'), which is what distinguishes ",
        "it from the separately screened 'number of other risk factors' ",
        "column NOTHRISK that excludes them. ",
        "DIRECTION OF CODING. The paper prints the covariate as 'number ",
        "of risk factors (<=1 vs. >1)' with a POSITIVE additive shift of ",
        "1.887 on the intercept for RISKCATN = 1, so RISKCATN = 1 is the ",
        "more-than-one group and the shift raises the progression risk, ",
        "which is the clinically coherent direction. The paper's own ",
        "numbers corroborate it: mixing the two final intercepts over the ",
        "published 69.8% / 30.2% split at zero exposure gives a mean ",
        "probability of 0.0386, against 0.0461 for the single-intercept ",
        "base model of Table S9; the reversed coding would give 0.0692, ",
        "which sits on the wrong side of the base model. The vignette ",
        "reproduces this arithmetic as an assertion."
      ),
      source_name        = "RISKCATN"
    )
  )

  # Covariates screened for this endpoint and NOT retained. Documentation
  # only; checkModelConventions() does not require these to be referenced
  # in model(). Sager 2023 Methods ('ER modeling') lists the full screen,
  # and MODEL CODE S2 tags them ';--contcovs-' and ';--catcovs-'.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate on the progression probability (source column AGE) and not retained. Median 50.0 years, range 15-92 (Sager 2023 Table S8). Age is already embedded in the retained NRISK_GT1 count, which is inclusive of age.",
      source_name = "AGE"
    ),
    BMI = list(
      description = "Body mass index at baseline.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a continuous covariate on the progression probability (source column BMIBL) and not retained. Median 31.0 kg/m^2, range 18-63 (Sager 2023 Table S8). Like age, BMI is already embedded in the retained NRISK_GT1 count. BMI IS retained in the companion population PK model, where it acts on the intramuscular absorption rate.",
      source_name = "BMIBL"
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical covariate on the progression probability (source column SEXF) and not retained. The COMET-TAIL exposure-response population was 409 of 902 male (45.3%; Sager 2023 Table S8). Sex IS retained in the companion population PK model, on both intramuscular bioavailability and the absorption rate.",
      source_name = "SEXF"
    ),
    SARS_VLOAD = list(
      description = "Baseline SARS-CoV-2 viral load.",
      units       = "log10 copies/mL",
      type        = "continuous",
      notes       = "Screened as a continuous covariate on the progression probability (source column VLCOVBL) and not retained. Median 6.09 log10 copies/mL, range 3.2-10.2 (Sager 2023 Table S8). Also screened, and likewise not retained, in the companion population PK model.",
      source_name = "VLCOVBL"
    ),
    T_SYMPTOM_ONSET = list(
      description = "Duration of COVID-19 symptoms before dosing, i.e. time since symptom onset.",
      units       = "days",
      type        = "continuous",
      notes       = "Screened as a continuous covariate on the progression probability (source column SYMDUR) and not retained. Median 4.0 days, range 0-7 (Sager 2023 Table S8). The paper also screened, and likewise did not retain, the categorical form of the same quantity (source column SYMDCTN; <=3 days 49.0%, 4-5 days 38.1%, >5 days 12.9%).",
      source_name = "SYMDUR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 902L,
    n_studies      = 1L,
    n_observations = "one record per patient: 2706 source observations from 902 patients, reduced from 2877 records on 959 patients by removing those outside the intent-to-treat population and those lacking PK data (Sager 2023 Results, 'ER modeling', and Table S7)",
    age_median     = "50.0 years (range 15-92)",
    bmi_median     = "31.0 kg/m^2 (range 18-63)",
    sex_female_pct = 54.7,
    disease_state  = "non-hospitalized mild-to-moderate COVID-19 at high risk of progression to hospitalization or death; 69.8% had one or fewer protocol-defined risk factors and 30.2% more than one",
    dose_range     = "single dose: 500 mg intravenous (n = 367), 500 mg intramuscular (n = 361) or 250 mg intramuscular (n = 174). COMET-TAIL had no placebo arm.",
    regions        = "COMET-TAIL was multinational but 85% of participants were recruited in Florida, USA",
    endpoint       = "Progression of COVID-19 through day 29: hospitalization for more than 24 h for acute management of illness due to any cause, or death. Observed rates by arm (Sager 2023 Table 2): 7 of 174 (4.0%) at 250 mg i.m., 8 of 361 (2.2%) at 500 mg i.m., 5 of 367 (1.4%) at 500 mg i.v., 20 of 902 (2.2%) overall.",
    notes          = paste0(
      "COMET-TAIL enrolled from June to September 2021, when the delta ",
      "(B.1.617.2) variant predominated; it was the variant detected in ",
      "88.2% of the 764 participants with sequencing data. Median time ",
      "from randomization to a progression event was 5.5 days, which is ",
      "why the 168 h concentration was carried forward over the 96 h ",
      "alternative. Additional covariates screened and not retained, ",
      "beyond those in covariatesDataExcluded, were the number of OTHER ",
      "risk factors excluding age and BMI (source column NOTHRISK: 67.0% ",
      "had none, 29.4% one, 3.2% two, 0.4% three), the categorical form ",
      "of symptom duration (SYMDCTN), the SARS-CoV-2 variant (VARIANTN), ",
      "the route of administration (ROUTEN; 40.7% i.v., 59.3% i.m.) and ",
      "the time to event (TTEVENT). None reached the alpha = 0.01 ",
      "forward-selection threshold."
    )
  )

  ini({
    # ==================================================================
    # All three values are the FINAL exposure-response estimates in
    # Sager 2023 Table 3 ("Parameter estimates and standard errors from
    # the final ER model ... sotrovimab concentrations at 168 h in the
    # COMET-TAIL study"). Minimum objective function value 170.913.
    #
    # The structural form is MODEL CODE S2 of the Supporting Information
    # (a NONMEM $PRED logistic regression, one record per patient):
    #     LOGIT = INT + SLP*CP168 + ETA(1)
    #     PROB  = EXP(LOGIT) / (1 + EXP(LOGIT))
    # with the likelihood Y = PROB if DV == 1 and Y = 1 - PROB if DV == 0.
    # MODEL CODE S2 is the BASE model, so it carries no risk-factor term;
    # the final model of Table 3 adds the additive intercept shift below,
    # described in Results ("Covariate analysis led to the addition of
    # the number of risk factors (<=1 vs. >1) as an additive shift on the
    # model intercept").
    #
    # NO RANDOM EFFECT. MODEL CODE S2 declares $OMEGA 0 FIX on ETA(1),
    # i.e. the between-patient random effect on the intercept is fixed at
    # zero, which is the only possibility with one record per patient.
    # It is therefore omitted here rather than carried as a degenerate
    # zero-variance eta. There is likewise no residual-error parameter:
    # the Bernoulli likelihood supplies the whole error model.
    #
    # CONCENTRATION UNIT. Table 3 prints the slope as "1/[ug/mL]", so
    # CONC_SOTRO_168H must be supplied in ug/mL -- the same units the
    # companion population PK model reports Cc in.
    # ==================================================================

    logitprog <- -4.169
    label("Logit of the probability of COVID-19 progression through day 29 at zero sotrovimab exposure, in patients with one or fewer risk factors (unitless)") # Sager 2023 Table 3, "Intercept / Overall response (logit)" (RSE 12.03%); expit(-4.169) = 0.01524

    e_nrisk_gt1_logitprog <- 1.887
    label("Additive shift on the progression logit for patients with more than one risk factor (unitless logit)") # Sager 2023 Table 3, "Additive shift in INT for RISKCATN = 1" (RSE 27.34%); expit(-4.169 + 1.887) = 0.09254

    e_conc_sotro_168h_logitprog <- -0.02037
    label("Effect of the 168 h serum sotrovimab concentration on the progression logit, per ug/mL") # Sager 2023 Table 3, "Slope for concentration at 168h (1/[ug/mL])" (RSE 54.43%)
  })

  model({
    # ----------------------------------------------------------------
    # Linear predictor on the logit scale. CONC_SOTRO_168H is UNCENTRED and
    # the risk-factor indicator shifts only the intercept -- Sager 2023
    # Results state explicitly that the number of risk factors "had no
    # impact on overall drug response", i.e. there is deliberately no
    # interaction term between NRISK_GT1 and CONC_SOTRO_168H.
    # ----------------------------------------------------------------
    lpprog <-
      logitprog +
      e_nrisk_gt1_logitprog * NRISK_GT1 +
      e_conc_sotro_168h_logitprog * CONC_SOTRO_168H

    # pprog is exposed as an output column so the fitted probabilities
    # can be compared against Tables 2, S10 and Figure 4 directly,
    # without having to draw Bernoulli samples.
    pprog <- expit(lpprog)

    prog ~ dbinom(1, pprog)
  })
}
