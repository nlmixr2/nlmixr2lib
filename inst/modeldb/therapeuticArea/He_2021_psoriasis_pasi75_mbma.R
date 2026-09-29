He_2021_psoriasis_pasi75_mbma <- function() {
  description <- paste0(
    "MBMA. Longitudinal model-based meta-analysis of the PASI75 responder ",
    "rate (proportion of patients achieving a 75% reduction from baseline in ",
    "the Psoriasis Area and Severity Index) over time in moderate-to-severe ",
    "plaque psoriasis, fitted to study-arm summary data from 80 randomised ",
    "trials (233 arms reporting PASI75, 40,323 patients) of 13 biologics and ",
    "4 small targeted molecules published up to July 2019. The logit of the ",
    "PASI75 rate is the sum of a placebo component that rises exponentially ",
    "to a plateau, E0 = BSL + A*(1 - exp(-kpbo*exp(eta)*t)), and a drug ",
    "component Edrug = Emax_d*(1 - exp(-k_d*t))*Dose/(Dose + ED50_d) with ",
    "Emax, ED50 and k estimated separately for each of 17 drugs (adalimumab, ",
    "infliximab, etanercept, certolizumab pegol, ustekinumab, briakinumab, ",
    "guselkumab, tildrakizumab, risankizumab, secukinumab, ixekizumab, ",
    "brodalumab, apremilast, tofacitinib, baricitinib, alefacept, ",
    "methotrexate). Arm-mean body weight acts as a power function centred at ",
    "90 kg on the placebo asymptote A. Dose is the per-administration ",
    "maintenance dose supplied in one CONMED_<drug>_DOSE covariate column per ",
    "drug; there is no PK layer and no rxode2 dose event. Between-study ",
    "variability is on the placebo onset rate. Simulation scope is ",
    "STUDY-ARM-MEAN responder trajectories, NOT individual patients. The ",
    "companion PASI90 model of the same paper is He_2021_psoriasis_pasi90_mbma."
  )
  reference <- paste(
    "He H, Wu W, Zhang Y, Zhang M, Sun N, Zhao L, Wang X.",
    "Model-Based Meta-Analysis in Psoriasis: A Quantitative Comparison of",
    "Biologics and Small Targeted Molecules.",
    "Front Pharmacol. 2021;12:586827. doi:10.3389/fphar.2021.586827.",
    "PMC8281289.",
    "Structural model: Methods 'Model Development' Equations 1-5.",
    "Covariate model: Equation 9. Residual model: Equations 7-8.",
    "Drug-specific estimates: Table 2. Placebo, covariate, random-effect and",
    "residual estimates: Supplementary Table S1 (Table1.docx of the",
    "Supplementary Material, EuropePMC supplementaryFiles for PMC8281289).",
    sep = " "
  )

  vignette <- "He_2021_psoriasis_mbma"

  # This model has no PK layer, no concentration and no rxode2 dose events:
  # dose is a per-arm covariate column. The placeholder `units` entries follow
  # the arm-level responder-rate MBMA convention of
  # Checchio_2017_psoriasis_pasi75_longitudinal_mbma, whose methodology He
  # 2021 states it follows (Methods 'Model Development').
  units <- list(
    time = "week (weeks since first dose; kpbo and k are in week^-1, and Tables 3 and 5 report Weeks 4, 8, 12, 16 and 24)",
    dosing = "mg/administration (maintenance dose per administration of the named agent, supplied in the CONMED_<drug>_DOSE covariate columns; infliximab is mg/kg per administration. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_pasi75 is the STUDY-ARM probability that a patient achieves a 75% reduction from baseline in the Psoriasis Area and Severity Index, on a 0-1 scale; it is NOT a drug concentration)"
  )

  covariateData <- list(
    WT = list(
      description = "Study-arm mean (or median) body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as the power ratio (WT/90)^e_wt_asym_pbo on the placebo asymptote A (He 2021 Equation 9 and Supplementary Table S1 row 'Body weight effect on A'). Equation 9 normalises by 'mean(covariate)'; 90 kg is the 'average body weight' named in the Figure 2-4 captions and the value at which every published prediction (Tables 3 and 5) is made. The dataset median is 89.6 kg (Table 1, 'Total' row), which changes the factor by 0.1%. The source set arms with unreported weight to the dataset median.",
      source_name = "body weight (He 2021 Equation 9; Supplementary Table S1)"
    ),
    N_ARM = list(
      description = "Number of patients contributing to the study-arm PASI75 proportion at a given timepoint.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity supplied per observation row. It is the meta-analytic weight of He 2021 Equation 8, Weight = sqrt(P*(1 - P)/N), which scales the residual of Equation 7. It does not affect the typical-value prediction.",
      source_name = "Nij (He 2021 Equations 1 and 8, 'the sample size in each arm of each trial')"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab maintenance dose per administration in the study arm; 0 if the arm did not receive adalimumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose PER ADMINISTRATION of the maintenance regimen (40 mg q2w), not the 80 mg loading dose: Dose = 40 reproduces Table 3 (Week 12: 68.1% vs published 67.8%), whereas 80 gives 83.0%. The Results text anchors the same reading ('the ED50 value was estimated to be 23.1 mg and the clinical dosage is 40 mg every 2 weeks').",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description = "Infliximab dose per administration in the study arm, per kilogram of body weight; 0 if the arm did not receive infliximab.",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = "UNIQUE UNITS within this model: infliximab is dosed in mg/kg (Table 1: '3 mg/kg 0, 2, 6' and '5 mg/kg 0, 2, 6, q8w') and its ED50 of 0.462 is therefore in mg/kg. Dose = 5 reproduces Table 3 (Week 12: 76.0% vs 75.65%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose per administration in the study arm; 0 if the arm did not receive etanercept.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (Table 1: 25 mg and 50 mg biw). Dose = 50 reproduces the Table 3 '50 mg biw' row (Week 12: 50.0% vs 49.6%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol maintenance dose per administration in the study arm; 0 if the arm did not receive certolizumab pegol.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Maintenance dose per administration. Table 3's regimen '400 mg 0, 2, 4 200 mg q2w' is reproduced by Dose = 200 (Week 12: 68.0% vs 67.7%); Dose = 400, the other Table 1 arm, gives 70.0%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description = "Ustekinumab dose per administration in the study arm; 0 if the arm did not receive ustekinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (Table 1: 45 and 90 mg at weeks 0, 4 then q12w). Dose = 45 reproduces Table 3 (Week 12: 65.7% vs 65.3%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_BRIAKINUMAB_DOSE = list(
      description = "Briakinumab maintenance dose per administration in the study arm; 0 if the arm did not receive briakinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Maintenance dose per administration (Table 1: 200 mg at weeks 0 and 4, then 100 mg q4w). Dose = 100 reproduces Table 3 (Week 12: 79.9% vs 79.5%). Briakinumab was discontinued in development (Table 1 footnote c).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_GUSELKUMAB_DOSE = list(
      description = "Guselkumab dose per administration in the study arm; 0 if the arm did not receive guselkumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (Table 1: 100 mg at weeks 0, 4 then q8w). Dose = 100 reproduces Table 3 (Week 12: 80.9% vs 80.4%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_TILDRAKIZUMAB_DOSE = list(
      description = "Tildrakizumab dose per administration in the study arm; 0 if the arm did not receive tildrakizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (Table 1: 100 and 200 mg at weeks 0, 4 then q12w). Dose = 100 reproduces Table 3 (Week 12: 61.7% vs 61.0%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_RISANKIZUMAB_DOSE = list(
      description = "Risankizumab dose per administration in the study arm; 0 if the arm did not receive risankizumab. Any positive value selects the full risankizumab effect.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the risankizumab ED50 is '0 FIX' in Table 2 ('the dose-response relationships were not obvious. Therefore, the ED50 for these drugs was fixed to 0'), so Dose/(Dose + 0) = 1 for every positive dose. The value records the arm's clinical dose (150 mg at weeks 0, 4 then q12w) for provenance and must not be read as a prediction for other doses. The companion PASI90 model DOES estimate a risankizumab ED50 (13.3 mg).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column; Table 2 'ED50 0 FIX'"
    ),
    CONMED_SECUKINUMAB_DOSE = list(
      description = "Secukinumab dose per administration in the study arm; 0 if the arm did not receive secukinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (Table 1: 150 and 300 mg at weeks 0, 1, 2, 3, 4 then q4w). Dose = 300 reproduces Table 3 (Week 12: 84.9% vs 84.5%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_IXEKIZUMAB_DOSE = list(
      description = "Ixekizumab maintenance dose per administration in the study arm; 0 if the arm did not receive ixekizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Maintenance dose per administration. Table 3 carries two ixekizumab rows that share one Emax / ED50 / k triplet and are separated only by this column: '160 mg 0, q4w' is Dose = 160 (Week 12: 86.2% vs 85.9%) and '160 mg 0 80 mg q4w' is Dose = 80 (82.8% vs 82.4%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_BRODALUMAB_DOSE = list(
      description = "Brodalumab dose per administration in the study arm; 0 if the arm did not receive brodalumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (Table 1: 140 and 210 mg at weeks 0, 1, 2 then q2w). Dose = 210 reproduces Table 3 (Week 12: 82.2% vs 81.7%). Unlike Checchio_2017_psoriasis_pasi_landmark_mbma, this model has no dosing-interval (REGI_Q4W) term.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_APREMILAST_DOSE = list(
      description = "Apremilast dose per administration in the study arm; 0 if the arm did not receive apremilast.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: Dose = 30 reproduces the Table 3 '30 mg b.i.d.' row (Week 12: 29.9% vs 29.5%) whereas Dose = 60 gives 60.9%. Consistent with Results: 'For all the drugs except apremilast, the dose regimen is higher than ED50' (ED50 98.5 mg).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose per administration in the study arm; 0 if the arm did not receive tofacitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: the Table 3 '5 mg b.i.d.' and '10 mg b.i.d.' rows share one parameter triplet and are reproduced by Dose = 5 (39.7% vs 39.2% at Week 12) and Dose = 10 (58.2% vs 57.7%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_BARICITINIB_DOSE = list(
      description = "Baricitinib daily dose in the study arm; 0 if the arm did not receive baricitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Once-daily dose (Table 1: 2, 4, 8 and 10 mg qd), so per administration equals per day. Dose = 10 reproduces Table 3 (Week 12: 48.9% vs 48.7%). Baricitinib data come from a single phase 2 trial (Table 1).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_ALEFACEPT_DOSE = list(
      description = "Alefacept dose per administration in the study arm; 0 if the arm did not receive alefacept. Any positive value selects the full alefacept effect.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED: alefacept ED50 is '0 FIX' in Table 2, so only the comparison > 0 is read. The value records the arm's clinical dose (Table 1: 10 or 15 mg i.m. qw) for provenance.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column; Table 2 'ED50 0 FIX'"
    ),
    CONMED_MTX_DOSE = list(
      description = "Methotrexate weekly dose in the study arm; 0 if the arm did not receive randomised methotrexate. Any positive value selects the full methotrexate effect.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED: methotrexate ED50 is '0 FIX' in Table 2, so only the comparison > 0 is read. The value records the arm's clinical dose (Table 1: 20 mg p.o. qw) for provenance.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column; Table 2 'ED50 0 FIX'"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 40323L,
    n_studies = 80L,
    age_range = "arm median ages 38.6-55.3 years (Table 1, 'Total' row: median 45)",
    weight_range = "arm median weights 66.6-99 kg (Table 1, 'Total' row: median 89.6 kg); predictions are made at 90 kg",
    sex_female_pct = 30.95,
    race_ethnicity = "not reported",
    disease_state = "adults with moderate to severe plaque psoriasis in randomised placebo- or active-controlled trials; arm median baseline PASI 11-33.1 and body surface area involved 19-50.2% (Table 1)",
    dose_range = "per Table 1: adalimumab 40 and 80 mg q2w; infliximab 3 and 5 mg/kg; etanercept 25 and 50 mg biw; certolizumab pegol 200 and 400 mg q2w; ustekinumab 45 and 90 mg q12w; briakinumab 100 mg q4w; guselkumab 100 mg q8w; tildrakizumab 100 and 200 mg q12w; risankizumab 150 mg q12w; secukinumab 150 and 300 mg q4w; ixekizumab 80 and 160 mg q4w; brodalumab 140 and 210 mg q2w; apremilast 10-30 mg b.i.d.; tofacitinib 5 and 10 mg b.i.d.; baricitinib 2-10 mg qd; alefacept 10 and 15 mg qw; methotrexate 20 mg qw",
    regions = "international; literature (PubMed, Cochrane, Embase, ClinicalTrials.gov) to 18 July 2019",
    notes = "MBMA at the STUDY-ARM level: 235 arms in total, 233 reporting PASI75 (Table 1). Female percentage is 100 minus the Table 1 'Total' median percentage male (69.05%). The random effect is BETWEEN-STUDY and this model must not be used to simulate individual patients."
  )

  ini({
    # ========================================================================
    # Structural model, He 2021 Methods 'Model Development'. The display
    # equations are dropped by a plain markdown conversion of the PDF; the
    # forms below were read with `pdftotext -layout` and the superscript
    # placement in Equation 4 was checked from the PDF glyph bounding boxes.
    #
    #   (1) N_response,ijt ~ binomial(N_ij, P_response,ijt)
    #   (2) P_response,ijt = g(E0 + Edrug)
    #   (3) g = 1 / (1 + exp(-(E0 + Edrug)))                 [inverse logit]
    #   (4) E0 = BSL + A * (1 - exp(-kpbo * time * exp(eta)))
    #   (5) Edrug = Emax * (1 - exp(-k * time)) * dose^c / (dose^c + ED50^c)
    #   (7) Obs_ijt = P_ijt + Weight * eps
    #   (8) Weight = sqrt(P_ijt * (1 - P_ijt) / N_ij)
    #   (9) Covariate effect = (Covariate / mean(covariate))^theta
    #
    # The Hill coefficient c is fixed to 1 (Methods: 'The Hill coefficient (c)
    # was fixed to 1 in the model because there was not sufficient
    # dose-response information available for each drug'), so it is not an
    # ini() parameter.
    #
    # Nothing below is tuned. With these values at 90 kg the typical-value
    # trajectory reproduces all 95 cells of Table 3 (19 regimens x 5 weeks),
    # which are NOT parameter estimates; the vignette tabulates every one.
    # ========================================================================

    # ---- Placebo component (Equation 4; Supplementary Table S1) ------------
    bsl_pbo <- -7.52
    label("Intercept of the placebo effect on the PASI75 logit scale (paper: BSL; unitless log-odds)")
    # Supplementary Table S1, 'Intercept of placebo effect (BSL)' = -7.52
    # (RSE 3.2%).

    asym_pbo <- 4.83
    label("Asymptote of the placebo effect on the PASI75 logit scale for a 90 kg arm (paper: A; unitless log-odds)")
    # Supplementary Table S1, 'Asymptote of placebo effect (A)' = 4.83
    # (RSE 4.5%). Placebo plateau logit -7.52 + 4.83 = -2.69, a 6.4% PASI75 rate.

    lkpbo <- log(0.25)
    label("Log placebo-effect onset rate constant (paper: kpbo; back-transform 0.25 /week)")
    # Supplementary Table S1, 'Rate of onset of placebo effect (kpbo)' = 0.25
    # (RSE 6.4%).

    e_wt_asym_pbo <- -0.245
    label("Power exponent on (WT/90) for the placebo asymptote A (unitless)")
    # Supplementary Table S1, 'Body weight effect on A' = -0.245 (RSE 44.1%).
    # Equation 9 power form; heavier arms reach a lower placebo plateau.

    # ========================================================================
    # DRUG-SPECIFIC PARAMETERS -- Table 2, 'Final parameter estimates of
    # PASI75 longitudinal model'. ED50 is in mg per administration except
    # infliximab (mg/kg); ED50 and k are held in log space so they stay
    # positive, the source reports both on the linear scale.
    #
    # ED50 '0 FIX' for risankizumab, alefacept and methotrexate (Results: 'the
    # dose-response relationships were not obvious. Therefore, the ED50 for
    # these drugs was fixed to 0'). An ED50 of exactly 0 makes
    # Dose/(Dose + ED50) = 1 at every positive dose, which model() writes as
    # the step (DOSE > 0); a log-scale ED50 cannot hold 0, so no ED50
    # parameter is declared for these three drugs.
    # ========================================================================

    # ---- TNF-alpha inhibitors ----
    emax_adalimumab <- 5.84
    label("Maximum adalimumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, adalimumab Emax = 5.84 (RSE 3.7%)
    led50_adalimumab <- log(23.1)
    label("Log adalimumab ED50 (log mg per administration); back-transform 23.1 mg") # Table 2, adalimumab ED50 = 23.1 (RSE 10%)
    lkdrug_adalimumab <- log(0.463)
    label("Log adalimumab drug-effect onset rate (log 1/week); back-transform 0.463 /week") # Table 2, adalimumab k = 0.463 (RSE 14.1%)

    emax_infliximab <- 4.46
    label("Maximum infliximab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, infliximab Emax = 4.46 (RSE 7%)
    led50_infliximab <- log(0.462)
    label("Log infliximab ED50 (log mg/kg per administration, NOT mg); back-transform 0.462 mg/kg") # Table 2, infliximab ED50 = 0.462 (RSE 48.9%)
    lkdrug_infliximab <- log(0.663)
    label("Log infliximab drug-effect onset rate (log 1/week); back-transform 0.663 /week") # Table 2, infliximab k = 0.663 (RSE 13%)

    emax_etanercept <- 4.34
    label("Maximum etanercept effect on the PASI75 logit scale (unitless log-odds)") # Table 2, etanercept Emax = 4.34 (RSE 10%)
    led50_etanercept <- log(21.5)
    label("Log etanercept ED50 (log mg per administration); back-transform 21.5 mg") # Table 2, etanercept ED50 = 21.5 (RSE 30.2%)
    lkdrug_etanercept <- log(0.282)
    label("Log etanercept drug-effect onset rate (log 1/week); back-transform 0.282 /week") # Table 2, etanercept k = 0.282 (RSE 10.7%)

    emax_certolizumab <- 3.95
    label("Maximum certolizumab pegol effect on the PASI75 logit scale (unitless log-odds)") # Table 2, certolizumab pegol Emax = 3.95 (RSE 3.4%)
    led50_certolizumab <- log(10.1)
    label("Log certolizumab pegol ED50 (log mg per administration); back-transform 10.1 mg") # Table 2, certolizumab pegol ED50 = 10.1 (RSE 70.2%)
    lkdrug_certolizumab <- log(0.327)
    label("Log certolizumab pegol drug-effect onset rate (log 1/week); back-transform 0.327 /week") # Table 2, certolizumab pegol k = 0.327 (RSE 21.7%)

    # ---- IL-12/23 inhibitors ----
    emax_ustekinumab <- 4.33
    label("Maximum ustekinumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, ustekinumab Emax = 4.33 (RSE 4.4%)
    led50_ustekinumab <- log(6.38)
    label("Log ustekinumab ED50 (log mg per administration); back-transform 6.38 mg") # Table 2, ustekinumab ED50 = 6.38 (RSE 41.5%)
    lkdrug_ustekinumab <- log(0.241)
    label("Log ustekinumab drug-effect onset rate (log 1/week); back-transform 0.241 /week") # Table 2, ustekinumab k = 0.241 (RSE 10.5%)

    emax_briakinumab <- 4.94
    label("Maximum briakinumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, briakinumab Emax = 4.94 (RSE 5%)
    led50_briakinumab <- log(12.1)
    label("Log briakinumab ED50 (log mg per administration); back-transform 12.1 mg") # Table 2, briakinumab ED50 = 12.1 (RSE 27%)
    lkdrug_briakinumab <- log(0.317)
    label("Log briakinumab drug-effect onset rate (log 1/week); back-transform 0.317 /week") # Table 2, briakinumab k = 0.317 (RSE 8.4%)

    # ---- IL-23 inhibitors ----
    emax_guselkumab <- 4.63
    label("Maximum guselkumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, guselkumab Emax = 4.63 (RSE 2.4%)
    led50_guselkumab <- log(2.75)
    label("Log guselkumab ED50 (log mg per administration); back-transform 2.75 mg") # Table 2, guselkumab ED50 = 2.75 (RSE 4%)
    lkdrug_guselkumab <- log(0.294)
    label("Log guselkumab drug-effect onset rate (log 1/week); back-transform 0.294 /week") # Table 2, guselkumab k = 0.294 (RSE 11.1%)

    emax_tildrakizumab <- 3.81
    label("Maximum tildrakizumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, tildrakizumab Emax = 3.81 (RSE 4%)
    led50_tildrakizumab <- log(4.63)
    label("Log tildrakizumab ED50 (log mg per administration); back-transform 4.63 mg") # Table 2, tildrakizumab ED50 = 4.63 (RSE 9.5%)
    lkdrug_tildrakizumab <- log(0.229)
    label("Log tildrakizumab drug-effect onset rate (log 1/week); back-transform 0.229 /week") # Table 2, tildrakizumab k = 0.229 (RSE 10.5%)

    emax_risankizumab <- 5.05
    label("Maximum risankizumab effect on the PASI75 logit scale, reached at any positive dose (unitless log-odds)") # Table 2, risankizumab Emax = 5.05 (RSE 1.8%); ED50 '0 FIX'
    lkdrug_risankizumab <- log(0.242)
    label("Log risankizumab drug-effect onset rate (log 1/week); back-transform 0.242 /week") # Table 2, risankizumab k = 0.242 (RSE 9.5%)

    # ---- IL-17 inhibitors ----
    emax_secukinumab <- 5.74
    label("Maximum secukinumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, secukinumab Emax = 5.74 (RSE 5.4%)
    led50_secukinumab <- log(69.2)
    label("Log secukinumab ED50 (log mg per administration); back-transform 69.2 mg") # Table 2, secukinumab ED50 = 69.2 (RSE 10.4%)
    lkdrug_secukinumab <- log(0.509)
    label("Log secukinumab drug-effect onset rate (log 1/week); back-transform 0.509 /week") # Table 2, secukinumab k = 0.509 (RSE 8%)

    emax_ixekizumab <- 5.05
    label("Maximum ixekizumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, ixekizumab Emax = 5.05 (RSE 2.3%)
    led50_ixekizumab <- log(9.74)
    label("Log ixekizumab ED50 (log mg per administration); back-transform 9.74 mg") # Table 2, ixekizumab ED50 = 9.74 (RSE 12%)
    lkdrug_ixekizumab <- log(1.16)
    label("Log ixekizumab drug-effect onset rate (log 1/week); back-transform 1.16 /week") # Table 2, ixekizumab k = 1.16 (RSE 9.2%)

    emax_brodalumab <- 7.69
    label("Maximum brodalumab effect on the PASI75 logit scale (unitless log-odds)") # Table 2, brodalumab Emax = 7.69 (RSE 8.4%)
    led50_brodalumab <- log(152)
    label("Log brodalumab ED50 (log mg per administration); back-transform 152 mg") # Table 2, brodalumab ED50 = 152 (RSE 16.6%)
    lkdrug_brodalumab <- log(1.14)
    label("Log brodalumab drug-effect onset rate (log 1/week); back-transform 1.14 /week") # Table 2, brodalumab k = 1.14 (RSE 7.5%)

    # ---- PDE4 inhibitor ----
    emax_apremilast <- fixed(9)
    label("Maximum apremilast effect on the PASI75 logit scale, held at the published value (unitless log-odds)") # Table 2, apremilast Emax = '9 FIX' (Results: 'The Emax for apremilast was fixed, otherwise the estimation for apremilast showed larger RSE%')
    led50_apremilast <- log(98.5)
    label("Log apremilast ED50 (log mg per administration); back-transform 98.5 mg") # Table 2, apremilast ED50 = 98.5 (RSE 6.1%)
    lkdrug_apremilast <- log(0.381)
    label("Log apremilast drug-effect onset rate (log 1/week); back-transform 0.381 /week") # Table 2, apremilast k = 0.381 (RSE 20.8%)

    # ---- JAK inhibitors ----
    emax_tofacitinib <- 4.65
    label("Maximum tofacitinib effect on the PASI75 logit scale (unitless log-odds)") # Table 2, tofacitinib Emax = 4.65 (RSE 6.2%)
    led50_tofacitinib <- log(4.23)
    label("Log tofacitinib ED50 (log mg per administration); back-transform 4.23 mg") # Table 2, tofacitinib ED50 = 4.23 (RSE 13.4%)
    lkdrug_tofacitinib <- log(0.498)
    label("Log tofacitinib drug-effect onset rate (log 1/week); back-transform 0.498 /week") # Table 2, tofacitinib k = 0.498 (RSE 13.1%)

    emax_baricitinib <- 5.63
    label("Maximum baricitinib effect on the PASI75 logit scale (unitless log-odds)") # Table 2, baricitinib Emax = 5.63 (RSE 19.2%)
    led50_baricitinib <- log(8.4)
    label("Log baricitinib ED50 (log mg per day); back-transform 8.4 mg") # Table 2, baricitinib ED50 = 8.4 (RSE 20.4%)
    lkdrug_baricitinib <- log(0.24)
    label("Log baricitinib drug-effect onset rate (log 1/week); back-transform 0.24 /week") # Table 2, baricitinib k = 0.24 (RSE 42.1%)

    # ---- CD2 antagonist (ED50 '0 FIX') ----
    emax_alefacept <- 1.38
    label("Maximum alefacept effect on the PASI75 logit scale, reached at any positive dose (unitless log-odds)") # Table 2, alefacept Emax = 1.38 (RSE 5.5%); ED50 '0 FIX'
    lkdrug_alefacept <- log(0.115)
    label("Log alefacept drug-effect onset rate (log 1/week); back-transform 0.115 /week") # Table 2, alefacept k = 0.115 (RSE 6.8%)

    # ---- Dihydrofolate reductase inhibitor (ED50 '0 FIX') ----
    emax_methotrexate <- 2.32
    label("Maximum methotrexate effect on the PASI75 logit scale, reached at any positive dose (unitless log-odds)") # Table 2, methotrexate Emax = 2.32 (RSE 4.1%); ED50 '0 FIX'
    lkdrug_methotrexate <- log(0.212)
    label("Log methotrexate drug-effect onset rate (log 1/week); back-transform 0.212 /week") # Table 2, methotrexate k = 0.212 (RSE 28.9%)

    # ========================================================================
    # BETWEEN-STUDY RANDOM EFFECT. An aggregate-data BETWEEN-STUDY variance,
    # NOT popPK between-subject variability: one draw describes one trial.
    #
    # PLACEMENT: on the placebo ONSET RATE kpbo, as Equation 4 prints it
    # (exp(eta) multiplies kpbo * time INSIDE the exponent; in the PDF the
    # glyphs '.time.exp(eta)' sit at the same superscript baseline as
    # '-kpbo'). Supplementary Table S1 labels the row 'omega(A), %', which
    # would put it on the asymptote. The paper's own VPC (Supplementary
    # Figure S3, digitised exactly from the vector PDF) settles it for the
    # equation: at Week 22 the risankizumab panel's 2.5th percentile is 84.0%
    # around a median of 91.3%. A 24.9% CV log-normal effect on A alone moves
    # the Week-22 typical value's 2.5th percentile to about 61% before any
    # residual error is added, which the VPC rules out; the same effect on
    # kpbo has almost vanished by Week 22 (the placebo term is at its plateau)
    # and leaves a band that residual error alone explains.
    #
    # SCALE: Supplementary Table S1 reports 'omega(A), %' = 24.9 (RSE 21.7%),
    # a CV%, converted to the log-normal variance by omega^2 = log(CV^2 + 1)
    # = log(1 + 0.249^2) = 0.06015. (Reading 24.9% as sqrt(omega^2) gives
    # 0.0620 instead; the difference is immaterial.)
    # ========================================================================
    eta_study_lkpbo ~ 0.06015

    # ========================================================================
    # RESIDUAL ERROR (Equations 7-8). Obs = P + Weight * eps, Weight =
    # sqrt(P*(1 - P)/N): the residual SD on the probability scale is the
    # binomial standard error of the arm proportion times sigma. Carried with
    # N_ARM as a covariate column, following
    # Checchio_2017_psoriasis_pasi75_longitudinal_mbma.
    # ========================================================================
    addSd_prob_pasi75 <- 1.33
    label("Multiplier on the binomial standard error sqrt(P*(1-P)/N_ARM) giving the residual SD of the arm PASI75 proportion (unitless)")
    # Supplementary Table S1, 'sigma' = 1.33 (RSE 7.4%). Read as the SD
    # (the table symbol is sigma, and Methods defines eps as having variance
    # sigma^2); reading it as the variance would give SD 1.153.
  })

  model({
    # Equation 9: power model on the placebo asymptote, centred at 90 kg.
    wtAsym <- (WT / 90)^e_wt_asym_pbo

    # Equation 4: placebo component; the between-study effect scales the
    # onset rate inside the exponent.
    kpbo <- exp(lkpbo + eta_study_lkpbo)
    e0 <- bsl_pbo + asym_pbo * wtAsym * (1 - exp(-kpbo * time))

    # Equation 5 with c = 1, one term per drug. Each arm supplies a positive
    # dose in exactly one CONMED_<drug>_DOSE column and zero in the rest; a
    # zero dose makes that drug's term exactly zero, so a placebo arm (all
    # columns zero) reduces to e0 alone.
    edAdalimumab <- emax_adalimumab * (1 - exp(-exp(lkdrug_adalimumab) * time)) *
      CONMED_ADALIMUMAB_DOSE / (CONMED_ADALIMUMAB_DOSE + exp(led50_adalimumab))
    edInfliximab <- emax_infliximab * (1 - exp(-exp(lkdrug_infliximab) * time)) *
      CONMED_INFLIXIMAB_DOSE / (CONMED_INFLIXIMAB_DOSE + exp(led50_infliximab))
    edEtanercept <- emax_etanercept * (1 - exp(-exp(lkdrug_etanercept) * time)) *
      CONMED_ETANERCEPT_DOSE / (CONMED_ETANERCEPT_DOSE + exp(led50_etanercept))
    edCertolizumab <- emax_certolizumab * (1 - exp(-exp(lkdrug_certolizumab) * time)) *
      CONMED_CERTOLIZUMAB_DOSE / (CONMED_CERTOLIZUMAB_DOSE + exp(led50_certolizumab))
    edUstekinumab <- emax_ustekinumab * (1 - exp(-exp(lkdrug_ustekinumab) * time)) *
      CONMED_USTEKINUMAB_DOSE / (CONMED_USTEKINUMAB_DOSE + exp(led50_ustekinumab))
    edBriakinumab <- emax_briakinumab * (1 - exp(-exp(lkdrug_briakinumab) * time)) *
      CONMED_BRIAKINUMAB_DOSE / (CONMED_BRIAKINUMAB_DOSE + exp(led50_briakinumab))
    edGuselkumab <- emax_guselkumab * (1 - exp(-exp(lkdrug_guselkumab) * time)) *
      CONMED_GUSELKUMAB_DOSE / (CONMED_GUSELKUMAB_DOSE + exp(led50_guselkumab))
    edTildrakizumab <- emax_tildrakizumab * (1 - exp(-exp(lkdrug_tildrakizumab) * time)) *
      CONMED_TILDRAKIZUMAB_DOSE / (CONMED_TILDRAKIZUMAB_DOSE + exp(led50_tildrakizumab))
    # ED50 '0 FIX': Dose/(Dose + 0) = 1 for any positive dose.
    edRisankizumab <- emax_risankizumab * (1 - exp(-exp(lkdrug_risankizumab) * time)) *
      (CONMED_RISANKIZUMAB_DOSE > 0)
    edSecukinumab <- emax_secukinumab * (1 - exp(-exp(lkdrug_secukinumab) * time)) *
      CONMED_SECUKINUMAB_DOSE / (CONMED_SECUKINUMAB_DOSE + exp(led50_secukinumab))
    edIxekizumab <- emax_ixekizumab * (1 - exp(-exp(lkdrug_ixekizumab) * time)) *
      CONMED_IXEKIZUMAB_DOSE / (CONMED_IXEKIZUMAB_DOSE + exp(led50_ixekizumab))
    edBrodalumab <- emax_brodalumab * (1 - exp(-exp(lkdrug_brodalumab) * time)) *
      CONMED_BRODALUMAB_DOSE / (CONMED_BRODALUMAB_DOSE + exp(led50_brodalumab))
    edApremilast <- emax_apremilast * (1 - exp(-exp(lkdrug_apremilast) * time)) *
      CONMED_APREMILAST_DOSE / (CONMED_APREMILAST_DOSE + exp(led50_apremilast))
    edTofacitinib <- emax_tofacitinib * (1 - exp(-exp(lkdrug_tofacitinib) * time)) *
      CONMED_TOFACITINIB_DOSE / (CONMED_TOFACITINIB_DOSE + exp(led50_tofacitinib))
    edBaricitinib <- emax_baricitinib * (1 - exp(-exp(lkdrug_baricitinib) * time)) *
      CONMED_BARICITINIB_DOSE / (CONMED_BARICITINIB_DOSE + exp(led50_baricitinib))
    # ED50 '0 FIX' for alefacept and methotrexate as well.
    edAlefacept <- emax_alefacept * (1 - exp(-exp(lkdrug_alefacept) * time)) *
      (CONMED_ALEFACEPT_DOSE > 0)
    edMethotrexate <- emax_methotrexate * (1 - exp(-exp(lkdrug_methotrexate) * time)) *
      (CONMED_MTX_DOSE > 0)

    edrug <- edAdalimumab + edInfliximab + edEtanercept + edCertolizumab +
      edUstekinumab + edBriakinumab + edGuselkumab + edTildrakizumab +
      edRisankizumab + edSecukinumab + edIxekizumab + edBrodalumab +
      edApremilast + edTofacitinib + edBaricitinib + edAlefacept +
      edMethotrexate

    # Equations 2-3: inverse logit of the summed placebo and drug components.
    lp_pasi75 <- e0 + edrug
    prob_pasi75 <- expit(lp_pasi75)

    # Equations 7-8: residual SD is sigma times the binomial standard error.
    sdArm <- addSd_prob_pasi75 * sqrt(prob_pasi75 * (1 - prob_pasi75) / N_ARM)
    prob_pasi75 ~ add(sdArm)
  })
}
