Shah_2026_rrmm_neutropenia_mbma <- function() {
  description <- paste0(
    "MBMA. Landmark model-based meta-analysis of the incidence of Grade >= 3 ",
    "neutropenia across treatment regimens for relapsed/refractory multiple ",
    "myeloma (RRMM), fitted in Monolix 2023R1 to study-arm-level event counts ",
    "from 130 published clinical trials (154 arms, 17,246 patients) in the ",
    "Certara CODEx RRMM database. The logit of the arm incidence is a fixed ",
    "reference intercept plus one regimen-specific effect (11 regimens, each ",
    "selected by its own mutually exclusive TRT_ indicator; lenalidomide ",
    "monotherapy enters as a linear slope in its daily dose) plus an additive ",
    "shift for a background alkylating agent, a between-study random effect ",
    "and a between-treatment-arm random effect scaled by 1/sqrt(arm size). ",
    "The source likelihood is binomial on the arm event count; this model ",
    "approximates it with an additive residual of SD sqrt(p * (1 - p) / ",
    "N_ARM). Simulation scope is the STUDY-ARM incidence, not an individual ",
    "patient's risk. Companion model of the same paper: ",
    "Shah_2026_rrmm_orr_mbma (overall response rate)."
  )

  reference <- paste(
    "Shah Z, Anderson C, McCormick K, Vallejo C, Duncan W, Hu C, Ratushny A,",
    "Zhou J, Kondic AG. Integrating MBMA and QSP to Identify Key Covariates",
    "and Predict Treatment Outcomes in Relapsed/Refractory Multiple Myeloma.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70145.",
    "doi:10.1002/psp4.70145. PMC12823780.",
    "Structural model: Methods 'MBMA Approach' Equations 1-3 and the",
    "Supporting Information Monolix code 'Model code for Grade >=3",
    "neutropenia'. Parameter estimates: Supporting Information Table S3.",
    "Reference intercept: digitised from Figure 3A (see the logite0 comment).",
    sep = " "
  )
  vignette <- "Shah_2026_rrmm_mbma"

  units <- list(
    time = paste0(
      "not applicable (LANDMARK meta-analysis: each observation is the ",
      "Grade >= 3 neutropenia incidence accumulated over a whole trial, so ",
      "the model has no time term)"
    ),
    dosing = paste0(
      "mg (lenalidomide daily dose of the study arm, supplied through the ",
      "DOSE_LENALIDOMIDE_MG covariate column and NOT as rxode2 dose events)"
    ),
    concentration = paste0(
      "fraction/arm (prob_neutropenia_grade3 is the study-arm proportion of ",
      "patients with Grade >= 3 neutropenia, on a 0-1 scale; it is NOT a ",
      "drug concentration)"
    )
  )

  trtNotes <- paste(
    "MBMA study-arm-level regimen indicator (a property of the published trial arm, not of an individual patient).",
    "The 11 regimen indicators of this model are MUTUALLY EXCLUSIVE: exactly one is 1 on every arm, because the source estimates one effect per regimen (Monolix 'trt_reg' switch, Supporting Information) and a combination regimen carries its own effect rather than the sum of its components.",
    "Background dexamethasone is not part of the regimen definition.",
    "With every indicator 0 the model returns the reference intercept only, which is not a clinically meaningful regimen."
  )

  covariateData <- list(
    TRT_BORTEZOMIB = list(
      description = "1 = study arm received bortezomib (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 1 (Monolix code); Table S3 'bortezomib'"
    ),
    TRT_CARFILZOMIB = list(
      description = "1 = study arm received carfilzomib (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 2 (Monolix code); Table S3 'carfilzomib'"
    ),
    TRT_CARFILZOMIB_LENALIDOMIDE = list(
      description = "1 = study arm received carfilzomib + lenalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 3 (Monolix code); Table S3 'carfilzomib +lenalidomide'"
    ),
    TRT_DARATUMUMAB = list(
      description = "1 = study arm received intravenous daratumumab (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 4 (Monolix code); Table S3 'daratumumab'"
    ),
    TRT_DARATUMUMAB_POMALIDOMIDE = list(
      description = "1 = study arm received daratumumab + pomalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 5 (Monolix code); Table S3 'daratumumab+pomalidomide'"
    ),
    TRT_DARATUMUMAB_HYALURONIDASE = list(
      description = "1 = study arm received subcutaneous daratumumab co-formulated with recombinant human hyaluronidase PH20 (rHuPH20), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste(
        trtNotes,
        "The source labels the regimen 'daratumumab+rHuPH20' and keeps it separate from intravenous daratumumab."
      ),
      source_name = "trt_reg == 6 (Monolix code); Table S3 'daratumumab+rHuPH20'"
    ),
    TRT_ISATUXIMAB_POMALIDOMIDE = list(
      description = "1 = study arm received isatuximab + pomalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 7 (Monolix code); Table S3 'isatuximab+pomalidomide'"
    ),
    TRT_IXAZOMIB = list(
      description = "1 = study arm received ixazomib (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 8 (Monolix code); Table S3 'ixazomib'"
    ),
    TRT_LENALIDOMIDE = list(
      description = "1 = study arm received lenalidomide (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste(
        trtNotes,
        "The lenalidomide regimen has no constant effect: its contribution is e_dose_lenalidomide_logite0 * DOSE_LENALIDOMIDE_MG, so a lenalidomide arm must also carry its dose."
      ),
      source_name = "trt_reg == 9 (Monolix code)"
    ),
    TRT_POMALIDOMIDE = list(
      description = "1 = study arm received pomalidomide (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 10 (Monolix code); Table S3 'pomalidomide'"
    ),
    TRT_TECLISTAMAB = list(
      description = "1 = study arm received teclistamab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 11 (Monolix code); Table S3 'teclistamab'"
    ),
    DOSE_LENALIDOMIDE_MG = list(
      description = "Lenalidomide daily dose of the study arm.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "MBMA study-arm-level dose covariate, read only on a lenalidomide-monotherapy arm (TRT_LENALIDOMIDE = 1); its value on other arms is ignored. The source normalised dose intensities to the FDA-approved dosing frequency (Methods 'MBMA Approach'), i.e. 25 mg once daily on days 1-21 of a 28-day cycle; Figures 3A and 4A evaluate lenalidomide at 25 mg (Figure 3A's 26.9% equals expit(-4 + 0.12 * 25)).",
      source_name = "dose (Monolix regressor, used for trt_reg == 9)"
    ),
    CONMED_ALKYLATOR = list(
      description = "1 = the study-arm regimen includes a background alkylating agent (e.g. melphalan, cyclophosphamide), 0 = no alkylating agent.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no background alkylating agent)",
      notes = "MBMA study-arm-level indicator. Adds 0.44 to the logit of the arm incidence (Table S3 'beta_backalkylating agent'; Figure S1A forest plot). The source does not list which alkylating agents were counted; the Discussion names melphalan and cyclophosphamide. Figure 4A evaluates every regimen with and without it.",
      source_name = "backalkylating agent (Table S3 'beta_backalkylating agent')"
    ),
    N_ARM = list(
      description = "Number of patients in the study arm (the denominator of the arm incidence).",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Used twice. (1) Equation 1: the between-treatment-arm random effect is divided by sqrt(N_ARM) (Monolix 'etaBAV/sqrt(N)'). (2) The binomial sampling SD of the observed arm proportion, sqrt(p * (1 - p) / N_ARM), is the residual of this model. The between-study random effect is not scaled.",
      source_name = "N (Monolix regressor); n_ij in Equation 1"
    )
  )

  covariatesDataExcluded <- list(
    LINE_2L_PCT = list(
      description = "Percentage of the arm's patients who had received exactly one prior line of therapy.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened for Grade >= 3 neutropenia (Table S5 'Prior.therapy.1', 119 arms, mean 29.1%) and not retained; it IS retained in the companion ORR model."
    ),
    AGE = list(
      description = "Study-arm mean age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S5 'Age', 149 arms, mean 64.97, range 54-76) and not retained. No estimate published."
    ),
    SEXM_PCT = list(
      description = "Percentage of male patients in the arm.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S5 'Males', 144 arms, mean 55.84%) and not retained."
    ),
    RACE_WHITE_PCT = list(
      description = "Percentage of white patients in the arm.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S5 'Race.white', 98 arms, mean 67.6%) and not retained."
    ),
    OTHER_SCREENED = list(
      description = "Other arm-level covariates screened in the stepwise BIC search.",
      units = "(various)",
      type = "continuous",
      reference_category = NULL,
      notes = "Table S5 also lists study year, treatment line, study phase, prior bortezomib %, prior lenalidomide %, prior anti-CD38 %, prior IMiD %, prior stem-cell transplant %, ISS stage 1/2/3 %, ECOG 0/1/2 %, high-risk cytogenetics % and disease duration. None was retained. Covariates with more than 40% missing data were not tested and missing values were median-imputed (Methods 'MBMA Approach')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 17246L,
    n_studies = 130L,
    n_arms = 154L,
    age_range = "study-arm mean age 54-76 years, mean 65 (Results; Table S5)",
    weight_range = "not reported",
    sex_female_pct = 44.2,
    race_ethnicity = "mean 67.6% white across arms (Results; Table S5)",
    disease_state = "relapsed/refractory multiple myeloma (second line or later); mean 29.1% of patients per arm had exactly one prior line of therapy",
    dose_range = "regimen doses as reported in each trial; dose intensities normalised to the FDA-approved dosing frequency",
    regimens = "11 regimens with an estimated effect: bortezomib, carfilzomib, carfilzomib + lenalidomide, daratumumab, daratumumab + pomalidomide, daratumumab + rHuPH20, isatuximab + pomalidomide, ixazomib, lenalidomide (linear dose-response), pomalidomide, teclistamab",
    regions = "not reported",
    notes = "MBMA at the study-arm level from the Certara CODEx RRMM clinical-trial outcomes database (352 RCTs, 512 arms overall), restricted to arms with more than 10 patients and to regimens supported by at least three independent studies. Table S5 reports 152 arms with covariate data; the Results text gives 154 arms. Prophylactic G-CSF use was not consistently reported and is not a covariate (Discussion), so regimens that used it may appear safer. The model simulates arm-level incidence for a trial of a given size; it is not an individual-patient model."
  )

  ini({
    # ================================================================
    # Shah 2026 Equations 1-2 and the Supporting Information Monolix
    # code ('Model code for Grade >=3 neutropenia'):
    #
    #   EFF_RE = EFF0 + EFF + etaBSV + etaBAV / sqrt(N)
    #   P      = invlogit(EFF_RE)
    #
    # EFF is the effect of the arm's regimen (one value per regimen;
    # lenalidomide monotherapy is slope_lena * dose). Every term is
    # additive on the logit scale, so no value is log-transformed.
    # ================================================================

    # REFERENCE INTERCEPT. Table S3's caption states 'EFF0 was fixed to
    # zero', but with EFF0 = 0 every regimen would have an incidence of
    # 82-98% (expit(1.55) to expit(4.05)), which Figures 2A, 3A and 4A
    # all contradict. Digitising the 11 point estimates of Figure 3A
    # (model-derived typical probabilities) and subtracting each
    # regimen's Table S3 value gives logit(p) - EFF = -3.99 to -4.02 for
    # every regimen, i.e. a common intercept of exactly -4. Figure 4A
    # (typical arm without alkylator: teclistamab ~50%, carfilzomib ~8%)
    # agrees. Encoded as fixed(-4); see the vignette Errata.
    logite0 <- fixed(-4)
    label("Reference logit intercept EFF0 of the arm Grade >= 3 neutropenia incidence (logit units)") # Back-solved from Figure 3A (digitised; all 11 regimens give -4.00 +/- 0.02); Table S3 caption says 'EFF0 was fixed to zero'

    e_trt_bortezomib_logite0 <- 1.72
    label("Effect of the bortezomib regimen on the logit of the arm incidence (logit units)") # Table S3 'bortezomib' = 1.72 (SE 0.098, RSE 5.67%)

    e_trt_carfilzomib_logite0 <- 1.55
    label("Effect of the carfilzomib regimen on the logit of the arm incidence (logit units)") # Table S3 'carfilzomib' = 1.55 (SE 0.15, RSE 9.97%)

    e_trt_carfilzomib_lenalidomide_logite0 <- 3.06
    label("Effect of the carfilzomib + lenalidomide regimen on the logit of the arm incidence (logit units)") # Table S3 'carfilzomib +lenalidomide' = 3.06 (SE 0.63, RSE 20.5%)

    e_trt_daratumumab_logite0 <- 2.14
    label("Effect of the daratumumab regimen on the logit of the arm incidence (logit units)") # Table S3 'daratumumab' = 2.14 (SE 0.24, RSE 11.3%)

    e_trt_daratumumab_pomalidomide_logite0 <- 3.97
    label("Effect of the daratumumab + pomalidomide regimen on the logit of the arm incidence (logit units)") # Table S3 'daratumumab+pomalidomide' = 3.97 (SE 0.45, RSE 11.3%)

    e_trt_daratumumab_hyaluronidase_logite0 <- 1.56
    label("Effect of the daratumumab + rHuPH20 (subcutaneous) regimen on the logit of the arm incidence (logit units)") # Table S3 'daratumumab+rHuPH20' = 1.56 (SE 0.69, RSE 44.3%)

    e_trt_isatuximab_pomalidomide_logite0 <- 3.99
    label("Effect of the isatuximab + pomalidomide regimen on the logit of the arm incidence (logit units)") # Table S3 'isatuximab+pomalidomide' = 3.99 (SE 0.28, RSE 6.96%)

    e_trt_ixazomib_logite0 <- 1.55
    label("Effect of the ixazomib regimen on the logit of the arm incidence (logit units)") # Table S3 'ixazomib' = 1.55 (SE 0.52, RSE 33.6%)

    e_dose_lenalidomide_logite0 <- 0.12
    label("Linear slope of the lenalidomide-monotherapy effect on the logit of the arm incidence per mg of daily lenalidomide dose (logit units/mg)") # Table S3 'slope_lenalidomide' = 0.12 (SE 0.0091, RSE 7.32%)

    e_trt_pomalidomide_logite0 <- 3.54
    label("Effect of the pomalidomide regimen on the logit of the arm incidence (logit units)") # Table S3 'pomalidomide' = 3.54 (SE 0.15, RSE 4.27%)

    e_trt_teclistamab_logite0 <- 4.05
    label("Effect of the teclistamab regimen on the logit of the arm incidence (logit units)") # Table S3 'teclistamab' = 4.05 (SE 0.35, RSE 8.71%)

    e_conmed_alkylator_logite0 <- 0.44
    label("Additive shift in the logit of the arm incidence for a background alkylating agent (logit units)") # Table S3 'beta_backalkylating agent' = 0.44 (SE 0.18, RSE 40.9%)

    # ================================================================
    # Two-level random effects on the logit (Equation 1). Monolix
    # reports random effects as standard deviations, so each ini()
    # value is the SD squared. The Figure 4A 90% PIs (e.g. lenalidomide
    # 25 mg about 13-52%) match an SD of about 0.6 on the logit scale,
    # not the 0.77 a variance reading would give.
    #   eta_study_logite0 : between-study variability (BSV), unscaled
    #   eta_arm_logite0   : between-treatment-arm variability (BTAV)
    #                       for an arm of ONE patient; model() divides
    #                       it by sqrt(N_ARM).
    # ================================================================
    eta_study_logite0 ~ 0.3481 # Table S3 'BSV' = 0.59 (SE 0.08, RSE 13.4%); variance 0.59^2
    eta_arm_logite0 ~ 10.5625 # Table S3 'BTAV' = 3.25 (SE 0.75, RSE 23.2%); variance 3.25^2, divided by sqrt(N_ARM) in model()
  })

  model({
    # Regimen effect (Monolix 'trt_reg' switch). Exactly one TRT_
    # indicator is 1 per arm. Lenalidomide monotherapy has no constant
    # term, only the linear dose slope.
    eff <- e_trt_bortezomib_logite0 * TRT_BORTEZOMIB +
      e_trt_carfilzomib_logite0 * TRT_CARFILZOMIB +
      e_trt_carfilzomib_lenalidomide_logite0 * TRT_CARFILZOMIB_LENALIDOMIDE +
      e_trt_daratumumab_logite0 * TRT_DARATUMUMAB +
      e_trt_daratumumab_pomalidomide_logite0 * TRT_DARATUMUMAB_POMALIDOMIDE +
      e_trt_daratumumab_hyaluronidase_logite0 * TRT_DARATUMUMAB_HYALURONIDASE +
      e_trt_isatuximab_pomalidomide_logite0 * TRT_ISATUXIMAB_POMALIDOMIDE +
      e_trt_ixazomib_logite0 * TRT_IXAZOMIB +
      e_dose_lenalidomide_logite0 * DOSE_LENALIDOMIDE_MG * TRT_LENALIDOMIDE +
      e_trt_pomalidomide_logite0 * TRT_POMALIDOMIDE +
      e_trt_teclistamab_logite0 * TRT_TECLISTAMAB

    # Background alkylating agent, additive on the logit (Table S3).
    effCov <- e_conmed_alkylator_logite0 * CONMED_ALKYLATOR

    # Between-arm term scaled by arm size; between-study term unscaled
    # (Equation 1).
    btav <- eta_arm_logite0 / sqrt(N_ARM)
    logitP <- logite0 + eff + effCov + eta_study_logite0 + btav

    # Study-arm incidence of Grade >= 3 neutropenia (Equation 2).
    prob_neutropenia_grade3 <- expit(logitP)

    # The source likelihood is binomial on the arm event count
    # (Equation 3). rxode2 cannot take the arm size as the binomial
    # size, so the observed arm proportion is given its binomial
    # sampling SD, sqrt(p * (1 - p) / N), as an additive residual.
    binomSd <- sqrt(prob_neutropenia_grade3 * (1 - prob_neutropenia_grade3) / N_ARM)
    prob_neutropenia_grade3 ~ add(binomSd)
  })
}
