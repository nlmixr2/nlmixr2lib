Shah_2026_rrmm_orr_mbma <- function() {
  description <- paste0(
    "MBMA. Landmark model-based meta-analysis of the overall response rate ",
    "(ORR, IMWG partial response or better) across treatment regimens for ",
    "relapsed/refractory multiple myeloma (RRMM), fitted in Monolix 2023R1 ",
    "to study-arm-level responder counts from 192 published clinical trials ",
    "(239 arms, 22,227 patients) in the Certara CODEx RRMM database. The ",
    "logit of the arm ORR is a fixed reference intercept plus one ",
    "regimen-specific effect (21 regimens, each selected by its own mutually ",
    "exclusive TRT_ indicator; carfilzomib monotherapy is linear in its ",
    "dose), an additive shift for background corticosteroid use, an Emax ",
    "function of the percentage of the arm's patients with exactly one ",
    "prior line of therapy, a between-study random effect and a ",
    "between-treatment-arm random effect scaled by 1/sqrt(arm size). The ",
    "source likelihood is binomial on the arm responder count; this model ",
    "approximates it with an additive residual of SD sqrt(p * (1 - p) / ",
    "N_ARM). Simulation scope is the STUDY-ARM response rate, not an ",
    "individual patient's response. Companion model of the same paper: ",
    "Shah_2026_rrmm_neutropenia_mbma (Grade >= 3 neutropenia)."
  )

  reference <- paste(
    "Shah Z, Anderson C, McCormick K, Vallejo C, Duncan W, Hu C, Ratushny A,",
    "Zhou J, Kondic AG. Integrating MBMA and QSP to Identify Key Covariates",
    "and Predict Treatment Outcomes in Relapsed/Refractory Multiple Myeloma.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70145.",
    "doi:10.1002/psp4.70145. PMC12823780.",
    "Structural model: Methods 'MBMA Approach' Equations 1-3 and the",
    "Supporting Information Monolix code 'Model code for ORR'. Parameter",
    "estimates: Supporting Information Table S4. Reference intercept:",
    "digitised from Figure 3B (see the logite0 comment).",
    sep = " "
  )
  vignette <- "Shah_2026_rrmm_mbma"

  units <- list(
    time = paste0(
      "not applicable (LANDMARK meta-analysis: each observation is the best ",
      "overall response rate of a whole trial arm, so the model has no time ",
      "term)"
    ),
    dosing = paste0(
      "mg/m^2 (carfilzomib dose per administration of the study arm, ",
      "supplied through the DOSE_CARFILZOMIB_MGM2 covariate column and NOT as ",
      "rxode2 dose events)"
    ),
    concentration = paste0(
      "fraction/arm (prob_orr is the study-arm proportion of patients with ",
      "an overall response, on a 0-1 scale; it is NOT a drug concentration)"
    )
  )

  trtNotes <- paste(
    "MBMA study-arm-level regimen indicator (a property of the published trial arm, not of an individual patient).",
    "The 21 regimen indicators of this model are MUTUALLY EXCLUSIVE: exactly one is 1 on every arm, because the source estimates one effect per regimen (Monolix 'trt_reg' switch, Supporting Information) and a combination regimen carries its own effect rather than the sum of its components.",
    "Background corticosteroid is not part of the regimen definition; it is the separate CONMED_STEROID covariate.",
    "With every indicator 0 the model returns the reference intercept only, which is not a clinically meaningful regimen."
  )

  covariateData <- list(
    TRT_BORTEZOMIB = list(
      description = "1 = study arm received bortezomib (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 1 (Monolix code); Table S4 'bortezomib'"
    ),
    TRT_BORTEZOMIB_LENALIDOMIDE = list(
      description = "1 = study arm received bortezomib + lenalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 2 (Monolix code); Table S4 'bortezomib+lenalidomide'"
    ),
    TRT_BORTEZOMIB_THALIDOMIDE = list(
      description = "1 = study arm received bortezomib + thalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 3 (Monolix code); Table S4 'bortezomib+thalidomide'"
    ),
    TRT_CARFILZOMIB = list(
      description = "1 = study arm received carfilzomib (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = paste(
        trtNotes,
        "The carfilzomib effect is inter_carfilzomib + slope_carfilzomib * DOSE_CARFILZOMIB_MGM2, so a carfilzomib arm must also carry its dose."
      ),
      source_name = "trt_reg == 4 (Monolix code); Table S4 'inter_carfilzomib' and 'slope_carfilzomib'"
    ),
    TRT_CARFILZOMIB_LENALIDOMIDE = list(
      description = "1 = study arm received carfilzomib + lenalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 5 (Monolix code); Table S4 'carfilzomib+lenalidomide'"
    ),
    TRT_DARATUMUMAB = list(
      description = "1 = study arm received intravenous daratumumab (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 6 (Monolix code); Table S4 'daratumumab'"
    ),
    TRT_BORTEZOMIB_DARATUMUMAB = list(
      description = "1 = study arm received daratumumab + bortezomib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 7 (Monolix code); Table S4 'daratumumab+bortezomib'"
    ),
    TRT_CARFILZOMIB_DARATUMUMAB = list(
      description = "1 = study arm received daratumumab + carfilzomib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 8 (Monolix code); Table S4 'daratumumab+carfilzomib'"
    ),
    TRT_DARATUMUMAB_LENALIDOMIDE = list(
      description = "1 = study arm received daratumumab + lenalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 9 (Monolix code); Table S4 'daratumumab+lenalidomide'"
    ),
    TRT_DARATUMUMAB_POMALIDOMIDE = list(
      description = "1 = study arm received daratumumab + pomalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 10 (Monolix code); Table S4 'daratumumab+pomalidomide'"
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
      source_name = "trt_reg == 11 (Monolix code); Table S4 'daratumumab+rHuPH20'"
    ),
    TRT_ELRANATAMAB = list(
      description = "1 = study arm received elranatamab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 12 (Monolix code); Table S4 'elranatamab'"
    ),
    TRT_ISATUXIMAB = list(
      description = "1 = study arm received isatuximab (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 13 (Monolix code); Table S4 'isatuximab'"
    ),
    TRT_CARFILZOMIB_ISATUXIMAB = list(
      description = "1 = study arm received isatuximab + carfilzomib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 14 (Monolix code); Table S4 'isatuximab+carfilzomib'"
    ),
    TRT_ISATUXIMAB_POMALIDOMIDE = list(
      description = "1 = study arm received isatuximab + pomalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 15 (Monolix code); Table S4 'isatuximab+pomalidomide'"
    ),
    TRT_IXAZOMIB = list(
      description = "1 = study arm received ixazomib (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 16 (Monolix code); Table S4 'ixazomib'"
    ),
    TRT_IXAZOMIB_POMALIDOMIDE = list(
      description = "1 = study arm received ixazomib + pomalidomide, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 17 (Monolix code); Table S4 'ixazomib+pomalidomide'"
    ),
    TRT_LENALIDOMIDE = list(
      description = "1 = study arm received lenalidomide (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 18 (Monolix code); Table S4 'lenalidomide'"
    ),
    TRT_POMALIDOMIDE = list(
      description = "1 = study arm received pomalidomide (without another regimen-defining agent), 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 19 (Monolix code); Table S4 'pomalidomide'"
    ),
    TRT_CARFILZOMIB_SELINEXOR = list(
      description = "1 = study arm received selinexor + carfilzomib, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 20 (Monolix code); Table S4 'selinexor+carfilzomib'"
    ),
    TRT_TECLISTAMAB = list(
      description = "1 = study arm received teclistamab, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = trtNotes,
      source_name = "trt_reg == 21 (Monolix code); Table S4 'teclistamab'"
    ),
    DOSE_CARFILZOMIB_MGM2 = list(
      description = "Carfilzomib dose per administration of the study arm.",
      units = "mg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "MBMA study-arm-level dose covariate, read only on a carfilzomib-monotherapy arm (TRT_CARFILZOMIB = 1); its value on other arms is ignored. The source normalised dose intensities to the FDA-approved dosing frequency (Methods 'MBMA Approach'). Figure 4 evaluates carfilzomib at 27 mg/m^2, and the digitised Figure 4C medians (about 41% without and 58% with background steroid at 12% one-prior-line patients) back-solve to 27-28 mg/m^2. Figure 3B's carfilzomib point (about 28.7%) corresponds to about 22 mg/m^2; the dose used for that figure is not stated.",
      source_name = "dose (Monolix regressor, used for trt_reg == 4)"
    ),
    CONMED_STEROID = list(
      description = "1 = the study-arm regimen includes a background corticosteroid (typically dexamethasone), 0 = no background corticosteroid.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no background corticosteroid)",
      notes = "MBMA STUDY-ARM-level grain of the general CONMED_STEROID indicator: a property of the published regimen, not of an individual patient. Adds 0.7 to the logit of the arm ORR (Table S4 'beta_backsteroid'; Figure S1B). The source's simulations hold it at 1 (steroid use, the most common value across trials) unless varied (Figure 4 caption).",
      source_name = "backsteroid (Table S4 'beta_backsteroid'; Figure 4C 'background steroid')"
    ),
    LINE_2L_PCT = list(
      description = "Percentage (0-100) of the arm's patients who had received exactly one prior line of therapy, i.e. who were treated in the second line.",
      units = "%",
      type = "continuous",
      reference_category = "0 (no patient in the second line; all had two or more prior lines)",
      notes = "MBMA study-arm-level covariate. Enters as Emax_PT * LINE_2L_PCT / (EC50_PT + LINE_2L_PCT) on the logit of the arm ORR (Monolix code). Supplied in PERCENT, not as a fraction: EC50_PT = 40.54 is on the 0-100 scale. Table S6 'Prior.therapy.1': 167 arms, mean 29.17%, median 12%, range 0-100; missing values were median-imputed (Methods). The Figure 4 simulations hold it at 12% (the median) unless varied, and Figure 4B contrasts 0% ('>= 2 prior lines') with 100% ('1 prior line'); both readings reproduce the digitised Figure 4B medians.",
      source_name = "ptherapy1 (Monolix regressor); Table S6 'Prior.therapy.1'; 'Ptherapy1%' in the Figure 4 caption"
    ),
    N_ARM = list(
      description = "Number of patients in the study arm (the denominator of the arm ORR).",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Used twice. (1) Equation 1: the between-treatment-arm random effect is divided by sqrt(N_ARM) (Monolix 'etaBAV/sqrt(N)'). (2) The binomial sampling SD of the observed arm proportion, sqrt(p * (1 - p) / N_ARM), is the residual of this model. The between-study random effect is not scaled.",
      source_name = "N (Monolix regressor); n_ij in Equation 1"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Study-arm mean age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S6 'Age', 225 arms, mean 64.72, range 48-76) and not retained. No estimate published."
    ),
    SEXM_PCT = list(
      description = "Percentage of male patients in the arm.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S6 'Males', 216 arms, mean 55.73%) and not retained."
    ),
    RACE_WHITE_PCT = list(
      description = "Percentage of white patients in the arm.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Table S6 'Race.white', 145 arms, mean 65.98%) and not retained."
    ),
    OTHER_SCREENED = list(
      description = "Other arm-level covariates screened in the stepwise BIC search.",
      units = "(various)",
      type = "continuous",
      reference_category = NULL,
      notes = "Table S6 also lists treatment line, study year, study phase, prior anti-CD38 %, prior bortezomib %, prior lenalidomide %, prior IMiD %, ISS stage 1/2/3 %, ECOG 0/1/2 % and high-risk cytogenetics %. None was retained. Covariates with more than 40% missing data were not tested and missing values were median-imputed (Methods 'MBMA Approach')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22227L,
    n_studies = 192L,
    n_arms = 239L,
    age_range = "study-arm mean age 48-76 years, mean 64.7 (Results; Table S6)",
    weight_range = "not reported",
    sex_female_pct = 44.3,
    race_ethnicity = "mean 66.0% white across arms (Results; Table S6)",
    disease_state = "relapsed/refractory multiple myeloma (second line or later); mean 29.2% of patients per arm had exactly one prior line of therapy",
    dose_range = "regimen doses as reported in each trial; dose intensities normalised to the FDA-approved dosing frequency",
    regimens = "21 regimens with an estimated effect: bortezomib, bortezomib + lenalidomide, bortezomib + thalidomide, carfilzomib (linear dose-response), carfilzomib + lenalidomide, daratumumab, daratumumab + bortezomib, daratumumab + carfilzomib, daratumumab + lenalidomide, daratumumab + pomalidomide, daratumumab + rHuPH20, elranatamab, isatuximab, isatuximab + carfilzomib, isatuximab + pomalidomide, ixazomib, ixazomib + pomalidomide, lenalidomide, pomalidomide, selinexor + carfilzomib, teclistamab",
    regions = "not reported",
    notes = "MBMA at the study-arm level from the Certara CODEx RRMM clinical-trial outcomes database (352 RCTs, 512 arms overall), restricted to arms with more than 10 patients and to regimens supported by at least three independent studies. Table S6 reports 237 arms with treatment-line data; the Results text gives 239 arms. ORR follows the IMWG response criteria. The model simulates arm-level response rates for a trial of a given size; it is not an individual-patient model."
  )

  ini({
    # ================================================================
    # Shah 2026 Equations 1-2 and the Supporting Information Monolix
    # code ('Model code for ORR'):
    #
    #   EFF_RE = EFF0 + EFF + Emax_PT * ptherapy1 / (EC50_PT + ptherapy1)
    #            + etaBSV + etaBAV / sqrt(N)
    #   P      = invlogit(EFF_RE)
    #
    # plus the background-steroid shift of Table S4, which the printed
    # code omits but Figures 4C and S1B show acting additively on the
    # logit. Every term is additive on the logit scale.
    # ================================================================

    # REFERENCE INTERCEPT. Table S4's caption states 'All effects are
    # modeled as deviations from the reference (EFF0 = 0)', but with
    # EFF0 = 0 every regimen would have an ORR of 97-100%, which
    # Figures 2B, 3B and 4 all contradict. Digitising the 20 non-dose
    # point estimates of Figure 3B (model-derived typical ORR, no
    # steroid and no one-prior-line patients) and subtracting each
    # regimen's Table S4 value gives logit(p) - EFF = -4.99 to -5.00
    # for every regimen, i.e. a common intercept of exactly -5. The
    # digitised Figure 4B/4C medians agree. Encoded as fixed(-5); see
    # the vignette Errata.
    logite0 <- fixed(-5)
    label("Reference logit intercept EFF0 of the arm overall response rate (logit units)") # Back-solved from Figure 3B (digitised; all 20 non-dose regimens give -5.00 +/- 0.01); Table S4 caption says 'EFF0 = 0'

    e_trt_bortezomib_logite0 <- 3.94
    label("Effect of the bortezomib regimen on the logit of the arm ORR (logit units)") # Table S4 'bortezomib' = 3.94 (SE 0.14, RSE 3.55%)

    e_trt_bortezomib_lenalidomide_logite0 <- 4.58
    label("Effect of the bortezomib + lenalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'bortezomib+lenalidomide' = 4.58 (SE 0.37, RSE 8.06%)

    e_trt_bortezomib_thalidomide_logite0 <- 5.5
    label("Effect of the bortezomib + thalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'bortezomib+thalidomide' = 5.5 (SE 0.38, RSE 6.90%)

    e_trt_carfilzomib_logite0 <- 3.55
    label("Intercept of the carfilzomib-monotherapy effect on the logit of the arm ORR (logit units)") # Table S4 'inter_carfilzomib' = 3.55 (SE 0.29, RSE 8.30%)

    e_dose_carfilzomib_logite0 <- 0.024
    label("Linear slope of the carfilzomib-monotherapy effect on the logit of the arm ORR per mg/m^2 of carfilzomib dose (logit units per mg/m^2)") # Table S4 'slope_carfilzomib' = 0.024 (SE 0.007, RSE 29.2%)

    e_trt_carfilzomib_lenalidomide_logite0 <- 5.1
    label("Effect of the carfilzomib + lenalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'carfilzomib+lenalidomide' = 5.1 (SE 0.23, RSE 4.52%)

    e_trt_daratumumab_logite0 <- 4.11
    label("Effect of the daratumumab regimen on the logit of the arm ORR (logit units)") # Table S4 'daratumumab' = 4.11 (SE 0.16, RSE 3.95%)

    e_trt_bortezomib_daratumumab_logite0 <- 4.9
    label("Effect of the daratumumab + bortezomib regimen on the logit of the arm ORR (logit units)") # Table S4 'daratumumab+bortezomib' = 4.9 (SE 0.24, RSE 4.93%)

    e_trt_carfilzomib_daratumumab_logite0 <- 5.3
    label("Effect of the daratumumab + carfilzomib regimen on the logit of the arm ORR (logit units)") # Table S4 'daratumumab+carfilzomib' = 5.3 (SE 0.26, RSE 4.87%)

    e_trt_daratumumab_lenalidomide_logite0 <- 5.52
    label("Effect of the daratumumab + lenalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'daratumumab+lenalidomide' = 5.52 (SE 0.29, RSE 5.21%)

    e_trt_daratumumab_pomalidomide_logite0 <- 4.54
    label("Effect of the daratumumab + pomalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'daratumumab+pomalidomide' = 4.54 (SE 0.22, RSE 4.90%)

    e_trt_daratumumab_hyaluronidase_logite0 <- 4.79
    label("Effect of the daratumumab + rHuPH20 (subcutaneous) regimen on the logit of the arm ORR (logit units)") # Table S4 'daratumumab+rHuPH20' = 4.79 (SE 0.3, RSE 6.33%)

    e_trt_elranatamab_logite0 <- 5.11
    label("Effect of the elranatamab regimen on the logit of the arm ORR (logit units)") # Table S4 'elranatamab' = 5.11 (SE 0.23, RSE 4.43%)

    e_trt_isatuximab_logite0 <- 3.33
    label("Effect of the isatuximab regimen on the logit of the arm ORR (logit units)") # Table S4 'isatuximab' = 3.33 (SE 0.34, RSE 10.1%)

    e_trt_carfilzomib_isatuximab_logite0 <- 5.48
    label("Effect of the isatuximab + carfilzomib regimen on the logit of the arm ORR (logit units)") # Table S4 'isatuximab+carfilzomib' = 5.48 (SE 0.29, RSE 5.36%)

    e_trt_isatuximab_pomalidomide_logite0 <- 4.79
    label("Effect of the isatuximab + pomalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'isatuximab+pomalidomide' = 4.79 (SE 0.2, RSE 4.23%)

    e_trt_ixazomib_logite0 <- 3.41
    label("Effect of the ixazomib regimen on the logit of the arm ORR (logit units)") # Table S4 'ixazomib' = 3.41 (SE 0.25, RSE 7.20%)

    e_trt_ixazomib_pomalidomide_logite0 <- 3.81
    label("Effect of the ixazomib + pomalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'ixazomib+pomalidomide' = 3.81 (SE 0.35, RSE 9.08%)

    e_trt_lenalidomide_logite0 <- 3.94
    label("Effect of the lenalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'lenalidomide' = 3.94 (SE 0.16, RSE 4.16%)

    e_trt_pomalidomide_logite0 <- 3.62
    label("Effect of the pomalidomide regimen on the logit of the arm ORR (logit units)") # Table S4 'pomalidomide' = 3.62 (SE 0.13, RSE 3.46%)

    e_trt_carfilzomib_selinexor_logite0 <- 4.84
    label("Effect of the selinexor + carfilzomib regimen on the logit of the arm ORR (logit units)") # Table S4 'selinexor+carfilzomib' = 4.84 (SE 0.37, RSE 7.72%)

    e_trt_teclistamab_logite0 <- 5.38
    label("Effect of the teclistamab regimen on the logit of the arm ORR (logit units)") # Table S4 'teclistamab' = 5.38 (SE 0.31, RSE 5.74%)

    e_conmed_steroid_logite0 <- 0.7
    label("Additive shift in the logit of the arm ORR for background corticosteroid use (logit units)") # Table S4 'beta_backsteroid' = 0.7 (SE 0.09, RSE 12.9%)

    emax_line_2l_pct <- 1.85
    label("Maximum additive effect on the logit of the arm ORR of the percentage of patients with exactly one prior line (logit units)") # Table S4 'Emax_Priortherapy1' = 1.85 (SE 0.41, RSE 22.2%)

    ec50_line_2l_pct <- 40.54
    label("Percentage of patients with exactly one prior line giving half the maximum effect (%)") # Table S4 'EC50_Priortherapy1' = 40.54 (SE 26.32, RSE 64.9%)

    # ================================================================
    # Two-level random effects on the logit (Equation 1). Monolix
    # reports random effects as standard deviations, so each ini()
    # value is the SD squared. The Figure 4C 90% PIs (e.g.
    # daratumumab + lenalidomide without steroid about 57-84%) match an
    # SD of about 0.42 on the logit scale, not the 0.64 a variance
    # reading would give.
    # ================================================================
    eta_study_logite0 ~ 0.1681 # Table S4 'BSV' = 0.41 (SE 0.052, RSE 12.7%); variance 0.41^2
    eta_arm_logite0 ~ 3.61 # Table S4 'BTAV' = 1.9 (SE 0.36, RSE 19%); variance 1.9^2, divided by sqrt(N_ARM) in model()
  })

  model({
    # Regimen effect (Monolix 'trt_reg' switch). Exactly one TRT_
    # indicator is 1 per arm. Carfilzomib monotherapy is linear in dose.
    eff <- e_trt_bortezomib_logite0 * TRT_BORTEZOMIB +
      e_trt_bortezomib_lenalidomide_logite0 * TRT_BORTEZOMIB_LENALIDOMIDE +
      e_trt_bortezomib_thalidomide_logite0 * TRT_BORTEZOMIB_THALIDOMIDE +
      (e_trt_carfilzomib_logite0 + e_dose_carfilzomib_logite0 * DOSE_CARFILZOMIB_MGM2) * TRT_CARFILZOMIB +
      e_trt_carfilzomib_lenalidomide_logite0 * TRT_CARFILZOMIB_LENALIDOMIDE +
      e_trt_daratumumab_logite0 * TRT_DARATUMUMAB +
      e_trt_bortezomib_daratumumab_logite0 * TRT_BORTEZOMIB_DARATUMUMAB +
      e_trt_carfilzomib_daratumumab_logite0 * TRT_CARFILZOMIB_DARATUMUMAB +
      e_trt_daratumumab_lenalidomide_logite0 * TRT_DARATUMUMAB_LENALIDOMIDE +
      e_trt_daratumumab_pomalidomide_logite0 * TRT_DARATUMUMAB_POMALIDOMIDE +
      e_trt_daratumumab_hyaluronidase_logite0 * TRT_DARATUMUMAB_HYALURONIDASE +
      e_trt_elranatamab_logite0 * TRT_ELRANATAMAB +
      e_trt_isatuximab_logite0 * TRT_ISATUXIMAB +
      e_trt_carfilzomib_isatuximab_logite0 * TRT_CARFILZOMIB_ISATUXIMAB +
      e_trt_isatuximab_pomalidomide_logite0 * TRT_ISATUXIMAB_POMALIDOMIDE +
      e_trt_ixazomib_logite0 * TRT_IXAZOMIB +
      e_trt_ixazomib_pomalidomide_logite0 * TRT_IXAZOMIB_POMALIDOMIDE +
      e_trt_lenalidomide_logite0 * TRT_LENALIDOMIDE +
      e_trt_pomalidomide_logite0 * TRT_POMALIDOMIDE +
      e_trt_carfilzomib_selinexor_logite0 * TRT_CARFILZOMIB_SELINEXOR +
      e_trt_teclistamab_logite0 * TRT_TECLISTAMAB

    # Background corticosteroid (Table S4) and the Emax effect of the
    # percentage of second-line patients (Monolix code), both additive
    # on the logit. LINE_2L_PCT is in percent (0-100).
    effCov <- e_conmed_steroid_logite0 * CONMED_STEROID +
      emax_line_2l_pct * LINE_2L_PCT / (ec50_line_2l_pct + LINE_2L_PCT)

    # Between-arm term scaled by arm size; between-study term unscaled
    # (Equation 1).
    btav <- eta_arm_logite0 / sqrt(N_ARM)
    logitP <- logite0 + eff + effCov + eta_study_logite0 + btav

    # Study-arm overall response rate (Equation 2).
    prob_orr <- expit(logitP)

    # The source likelihood is binomial on the arm responder count
    # (Equation 3). rxode2 cannot take the arm size as the binomial
    # size, so the observed arm proportion is given its binomial
    # sampling SD, sqrt(p * (1 - p) / N), as an additive residual.
    binomSd <- sqrt(prob_orr * (1 - prob_orr) / N_ARM)
    prob_orr ~ add(binomSd)
  })
}
