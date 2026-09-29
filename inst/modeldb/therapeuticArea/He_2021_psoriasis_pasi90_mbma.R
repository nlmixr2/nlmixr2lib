He_2021_psoriasis_pasi90_mbma <- function() {
  description <- paste0(
    "MBMA. Longitudinal model-based meta-analysis of the PASI90 responder ",
    "rate (proportion of patients achieving a 90% reduction from baseline in ",
    "the Psoriasis Area and Severity Index) over time in moderate-to-severe ",
    "plaque psoriasis, fitted to study-arm summary data from randomised ",
    "trials (224 arms reporting PASI90 out of 235 arms from 80 trials and ",
    "40,323 patients) of 12 biologics and 4 small targeted molecules ",
    "published up to July 2019. The logit of the PASI90 rate is the sum of a ",
    "placebo component that rises exponentially to a plateau, E0 = BSL + ",
    "A*(1 - exp(-kpbo*exp(eta)*t)), and a drug component Edrug = ",
    "Emax_d*(1 - exp(-k_d*t))*Dose/(Dose + ED50_d) with Emax, ED50 and k ",
    "estimated separately for each of 16 drugs (adalimumab, infliximab, ",
    "etanercept, certolizumab pegol, ustekinumab, briakinumab, guselkumab, ",
    "tildrakizumab, risankizumab, secukinumab, ixekizumab, brodalumab, ",
    "apremilast, tofacitinib, baricitinib, methotrexate; alefacept had too ",
    "few PASI90 data). Arm-mean body weight acts as a power function centred ",
    "at 90 kg on the placebo asymptote A. Dose is the per-administration ",
    "maintenance dose supplied in one CONMED_<drug>_DOSE covariate column per ",
    "drug; there is no PK layer and no rxode2 dose event. Between-study ",
    "variability is on the placebo onset rate. Simulation scope is ",
    "STUDY-ARM-MEAN responder trajectories, NOT individual patients. The ",
    "companion PASI75 model of the same paper is He_2021_psoriasis_pasi75_mbma."
  )
  reference <- paste(
    "He H, Wu W, Zhang Y, Zhang M, Sun N, Zhao L, Wang X.",
    "Model-Based Meta-Analysis in Psoriasis: A Quantitative Comparison of",
    "Biologics and Small Targeted Molecules.",
    "Front Pharmacol. 2021;12:586827. doi:10.3389/fphar.2021.586827.",
    "PMC8281289.",
    "Structural model: Methods 'Model Development' Equations 1-5.",
    "Covariate model: Equation 9. Residual model: Equations 7-8.",
    "Drug-specific estimates: Table 4. Placebo, covariate, random-effect and",
    "residual estimates: Supplementary Table S2 (Table2.docx of the",
    "Supplementary Material, PMC8281289).",
    sep = " "
  )

  vignette <- "He_2021_psoriasis_mbma"

  # Same arm-level responder-rate MBMA layout as the companion
  # He_2021_psoriasis_pasi75_mbma: no PK layer, no concentration, no rxode2
  # dose events; dose is a per-arm covariate column.
  units <- list(
    time = "week (weeks since first dose; kpbo and k are in week^-1, and Table 5 reports Weeks 4, 8, 12, 16 and 24)",
    dosing = "mg/administration (maintenance dose per administration of the named agent, supplied in the CONMED_<drug>_DOSE covariate columns; infliximab is mg/kg per administration. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_pasi90 is the STUDY-ARM probability that a patient achieves a 90% reduction from baseline in the Psoriasis Area and Severity Index, on a 0-1 scale; it is NOT a drug concentration)"
  )

  # The dose metric of every CONMED_<drug>_DOSE column is the one settled in
  # He_2021_psoriasis_pasi75_mbma against Table 3; each note below gives the
  # independent Week-12 check against this model's Table 5.
  covariateData <- list(
    WT = list(
      description = "Study-arm mean (or median) body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as the power ratio (WT/90)^e_wt_asym_pbo on the placebo asymptote A (He 2021 Equation 9 and Supplementary Table S2 row 'Body weight effect on A'). 90 kg is the 'average body weight' named in the Figure 4 caption and the value at which every Table 5 prediction is made ('assuming a typical body weight with 90 kg'). The dataset median is 89.6 kg (Table 1). The source set arms with unreported weight to the dataset median.",
      source_name = "body weight (He 2021 Equation 9; Supplementary Table S2)"
    ),
    N_ARM = list(
      description = "Number of patients contributing to the study-arm PASI90 proportion at a given timepoint.",
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
      notes = "Maintenance dose per administration (40 mg q2w). Dose = 40 gives a Week-12 PASI90 of 44.1% against the Table 5 '80 mg 0 40 mg 1, q2w' value of 43.5%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description = "Infliximab dose per administration in the study arm, per kilogram of body weight; 0 if the arm did not receive infliximab.",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = "UNIQUE UNITS within this model: mg/kg, so the ED50 of 0.689 is in mg/kg. Dose = 5 gives 49.9% at Week 12 against Table 5's 49.15%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose per administration in the study arm; 0 if the arm did not receive etanercept.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (25 or 50 mg biw). Dose = 50 gives 22.5% at Week 12 against Table 5's 22.0%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol maintenance dose per administration in the study arm; 0 if the arm did not receive certolizumab pegol.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Maintenance dose per administration. Dose = 200 gives 36.9% at Week 12 against Table 5's 36.5% for '400 mg 0, 2, 4 200 mg q2w'; Dose = 400 gives 38.6%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description = "Ustekinumab dose per administration in the study arm; 0 if the arm did not receive ustekinumab. Any positive value selects the full ustekinumab effect.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: Table 4 reports the ustekinumab ED50 as '0 FIX' (Results: 'The dose-response difference of ustekinumab was not significant and estimating ED50 resulted in a poor model estimation accuracy. Therefore, ED50 for ustekinumab was fixed to 0'). Only the comparison > 0 is read; the value records the arm's clinical dose (45 or 90 mg q12w). The companion PASI75 model DOES estimate an ustekinumab ED50 (6.38 mg). Week-12 check: 42.0% vs Table 5's 41.6%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column; Table 4 'ED50 0 FIX'"
    ),
    CONMED_BRIAKINUMAB_DOSE = list(
      description = "Briakinumab maintenance dose per administration in the study arm; 0 if the arm did not receive briakinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Maintenance dose per administration (100 mg q4w after 200 mg at weeks 0 and 4). Dose = 100 gives 55.9% at Week 12 against Table 5's 54.85%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_GUSELKUMAB_DOSE = list(
      description = "Guselkumab dose per administration in the study arm; 0 if the arm did not receive guselkumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (100 mg at weeks 0, 4 then q8w). Dose = 100 gives 61.6% at Week 12 against Table 5's 61.0%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_TILDRAKIZUMAB_DOSE = list(
      description = "Tildrakizumab dose per administration in the study arm; 0 if the arm did not receive tildrakizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (100 or 200 mg at weeks 0, 4 then q12w). Dose = 100 gives 35.0% at Week 12 against Table 5's 34.8%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_RISANKIZUMAB_DOSE = list(
      description = "Risankizumab dose per administration in the study arm; 0 if the arm did not receive risankizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED in this model (ED50 13.3 mg, Table 4), unlike the companion PASI75 model where the risankizumab ED50 is '0 FIX'. Dose = 150 gives 66.3% at Week 12 against Table 5's 65.5%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_SECUKINUMAB_DOSE = list(
      description = "Secukinumab dose per administration in the study arm; 0 if the arm did not receive secukinumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (150 or 300 mg at weeks 0-4 then q4w). Dose = 300 gives 63.4% at Week 12 against Table 5's 62.5%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_IXEKIZUMAB_DOSE = list(
      description = "Ixekizumab maintenance dose per administration in the study arm; 0 if the arm did not receive ixekizumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Maintenance dose per administration. Table 5's '160 mg 0, q4w' row is Dose = 160 (68.2% vs 67.2% at Week 12) and its '160 mg 0 80 mg q4w' row is Dose = 80 (63.7% vs 62.6%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_BRODALUMAB_DOSE = list(
      description = "Brodalumab dose per administration in the study arm; 0 if the arm did not receive brodalumab.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration (140 or 210 mg at weeks 0, 1, 2 then q2w). Dose = 210 gives 65.5% at Week 12 against Table 5's 64.5%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_APREMILAST_DOSE = list(
      description = "Apremilast dose per administration in the study arm; 0 if the arm did not receive apremilast.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose (30 mg b.i.d. is Dose = 30). Dose = 30 gives 6.5% at Week 12 against Table 5's 6.09%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose per administration in the study arm; 0 if the arm did not receive tofacitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Dose per administration, NOT total daily dose: Table 5's '5 mg b.i.d.' and '10 mg b.i.d.' rows are Dose = 5 (20.4% vs 19.5% at Week 12) and Dose = 10 (34.4% vs 33.6%).",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_BARICITINIB_DOSE = list(
      description = "Baricitinib daily dose in the study arm; 0 if the arm did not receive baricitinib.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Once-daily dose (2, 4, 8 or 10 mg qd). Dose = 10 gives 27.0% at Week 12 against Table 5's 26.2%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column"
    ),
    CONMED_MTX_DOSE = list(
      description = "Methotrexate weekly dose in the study arm; 0 if the arm did not receive randomised methotrexate. Any positive value selects the full methotrexate effect.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED: methotrexate ED50 is '0 FIX' in Table 4, so only the comparison > 0 is read. The value records the arm's clinical dose (20 mg p.o. qw). Week-12 check: 13.4% vs Table 5's 13.05%.",
      source_name = "dose (He 2021 Equation 5); Table 1 'Route (regimen)' column; Table 4 'ED50 0 FIX'"
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
    dose_range = "per Table 1 (alefacept arms excluded from this model): adalimumab 40 and 80 mg q2w; infliximab 3 and 5 mg/kg; etanercept 25 and 50 mg biw; certolizumab pegol 200 and 400 mg q2w; ustekinumab 45 and 90 mg q12w; briakinumab 100 mg q4w; guselkumab 100 mg q8w; tildrakizumab 100 and 200 mg q12w; risankizumab 150 mg q12w; secukinumab 150 and 300 mg q4w; ixekizumab 80 and 160 mg q4w; brodalumab 140 and 210 mg q2w; apremilast 10-30 mg b.i.d.; tofacitinib 5 and 10 mg b.i.d.; baricitinib 2-10 mg qd; methotrexate 20 mg qw",
    regions = "international; literature (PubMed, Cochrane, Embase, ClinicalTrials.gov) to 18 July 2019",
    notes = "MBMA at the STUDY-ARM level: 224 of the 235 arms report PASI90 (Table 1); n_subjects and n_studies are the whole-dataset counts, the paper does not give a PASI90-only patient total. Female percentage is 100 minus the Table 1 'Total' median percentage male (69.05%). The random effect is BETWEEN-STUDY and this model must not be used to simulate individual patients."
  )

  ini({
    # ========================================================================
    # Structural model: identical equations to He_2021_psoriasis_pasi75_mbma
    # (He 2021 Methods Equations 1-9, Hill coefficient c fixed to 1), fitted
    # separately to the PASI90 end point ('Two independent longitudinal
    # model-based meta-analyses', Discussion). See that file for the equation
    # listing.
    #
    # Nothing below is tuned. With these values at 90 kg the typical-value
    # trajectory reproduces all 90 cells of Table 5 (18 regimens x 5 weeks);
    # the vignette tabulates every one.
    # ========================================================================

    # ---- Placebo component (Equation 4; Supplementary Table S2) ------------
    bsl_pbo <- -10.2
    label("Intercept of the placebo effect on the PASI90 logit scale (paper: BSL; unitless log-odds)")
    # Supplementary Table S2, 'Intercept of placebo effect (BSL)' = -10.2
    # (RSE 6.3%).

    asym_pbo <- 6.32
    label("Asymptote of the placebo effect on the PASI90 logit scale for a 90 kg arm (paper: A; unitless log-odds)")
    # Supplementary Table S2, 'Asymptote of placebo effect (A)' = 6.32
    # (RSE 9.5%). Placebo plateau logit -10.2 + 6.32 = -3.88, a 2.0% PASI90 rate.

    lkpbo <- log(0.259)
    label("Log placebo-effect onset rate constant (paper: kpbo; back-transform 0.259 /week)")
    # Supplementary Table S2, 'Rate of onset of placebo effect (kpbo)' = 0.259
    # (RSE 7.4%).

    e_wt_asym_pbo <- -0.214
    label("Power exponent on (WT/90) for the placebo asymptote A (unitless)")
    # Supplementary Table S2, 'Body weight effect on A' = -0.214 (RSE 33.6%).

    # ========================================================================
    # DRUG-SPECIFIC PARAMETERS -- Table 4, 'Final parameter estimates of
    # PASI90 longitudinal model'. ED50 in mg per administration except
    # infliximab (mg/kg); ED50 and k held in log space.
    #
    # ED50 '0 FIX' for ustekinumab and methotrexate: written as the step
    # (DOSE > 0) in model(), with no ED50 parameter declared (a log-scale
    # ED50 cannot hold 0).
    # ========================================================================

    # ---- TNF-alpha inhibitors ----
    emax_adalimumab <- 5.74
    label("Maximum adalimumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, adalimumab Emax = 5.74 (RSE 7.6%)
    led50_adalimumab <- log(18.3)
    label("Log adalimumab ED50 (log mg per administration); back-transform 18.3 mg") # Table 4, adalimumab ED50 = 18.3 (RSE 25.7%)
    lkdrug_adalimumab <- log(0.487)
    label("Log adalimumab drug-effect onset rate (log 1/week); back-transform 0.487 /week") # Table 4, adalimumab k = 0.487 (RSE 26.5%)

    emax_infliximab <- 4.73
    label("Maximum infliximab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, infliximab Emax = 4.73 (RSE 4.6%)
    led50_infliximab <- log(0.689)
    label("Log infliximab ED50 (log mg/kg per administration, NOT mg); back-transform 0.689 mg/kg") # Table 4, infliximab ED50 = 0.689 (RSE 23.8%)
    lkdrug_infliximab <- log(0.802)
    label("Log infliximab drug-effect onset rate (log 1/week); back-transform 0.802 /week") # Table 4, infliximab k = 0.802 (RSE 35.2%)

    emax_etanercept <- 5.33
    label("Maximum etanercept effect on the PASI90 logit scale (unitless log-odds)") # Table 4, etanercept Emax = 5.33 (RSE 11.2%)
    led50_etanercept <- log(30.7)
    label("Log etanercept ED50 (log mg per administration); back-transform 30.7 mg") # Table 4, etanercept ED50 = 30.7 (RSE 25.3%)
    lkdrug_etanercept <- log(0.181)
    label("Log etanercept drug-effect onset rate (log 1/week); back-transform 0.181 /week") # Table 4, etanercept k = 0.181 (RSE 12.7%)

    emax_certolizumab <- 3.99
    label("Maximum certolizumab pegol effect on the PASI90 logit scale (unitless log-odds)") # Table 4, certolizumab pegol Emax = 3.99 (RSE 4.9%)
    led50_certolizumab <- log(8.1)
    label("Log certolizumab pegol ED50 (log mg per administration); back-transform 8.1 mg") # Table 4, certolizumab pegol ED50 = 8.1 (RSE 63.2%)
    lkdrug_certolizumab <- log(0.242)
    label("Log certolizumab pegol drug-effect onset rate (log 1/week); back-transform 0.242 /week") # Table 4, certolizumab pegol k = 0.242 (RSE 21.2%)

    # ---- IL-12/23 inhibitors ----
    emax_ustekinumab <- 4.06
    label("Maximum ustekinumab effect on the PASI90 logit scale, reached at any positive dose (unitless log-odds)") # Table 4, ustekinumab Emax = 4.06 (RSE 3.3%); ED50 '0 FIX'
    lkdrug_ustekinumab <- log(0.243)
    label("Log ustekinumab drug-effect onset rate (log 1/week); back-transform 0.243 /week") # Table 4, ustekinumab k = 0.243 (RSE 11.6%)

    emax_briakinumab <- 4.9
    label("Maximum briakinumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, briakinumab Emax = 4.9 (RSE 4%)
    led50_briakinumab <- log(9.47)
    label("Log briakinumab ED50 (log mg per administration); back-transform 9.47 mg") # Table 4, briakinumab ED50 = 9.47 (RSE 20.1%)
    lkdrug_briakinumab <- log(0.341)
    label("Log briakinumab drug-effect onset rate (log 1/week); back-transform 0.341 /week") # Table 4, briakinumab k = 0.341 (RSE 14.3%)

    # ---- IL-23 inhibitors ----
    emax_guselkumab <- 4.78
    label("Maximum guselkumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, guselkumab Emax = 4.78 (RSE 3.8%)
    led50_guselkumab <- log(2.95)
    label("Log guselkumab ED50 (log mg per administration); back-transform 2.95 mg") # Table 4, guselkumab ED50 = 2.95 (RSE 6.4%)
    lkdrug_guselkumab <- log(0.537)
    label("Log guselkumab drug-effect onset rate (log 1/week); back-transform 0.537 /week") # Table 4, guselkumab k = 0.537 (RSE 43%)

    emax_tildrakizumab <- 3.99
    label("Maximum tildrakizumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, tildrakizumab Emax = 3.99 (RSE 5.4%)
    led50_tildrakizumab <- log(7.16)
    label("Log tildrakizumab ED50 (log mg per administration); back-transform 7.16 mg") # Table 4, tildrakizumab ED50 = 7.16 (RSE 25.8%)
    lkdrug_tildrakizumab <- log(0.252)
    label("Log tildrakizumab drug-effect onset rate (log 1/week); back-transform 0.252 /week") # Table 4, tildrakizumab k = 0.252 (RSE 45.6%)

    emax_risankizumab <- 5.56
    label("Maximum risankizumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, risankizumab Emax = 5.56 (RSE 2.8%)
    led50_risankizumab <- log(13.3)
    label("Log risankizumab ED50 (log mg per administration); back-transform 13.3 mg") # Table 4, risankizumab ED50 = 13.3 (RSE 23.9%)
    lkdrug_risankizumab <- log(0.246)
    label("Log risankizumab drug-effect onset rate (log 1/week); back-transform 0.246 /week") # Table 4, risankizumab k = 0.246 (RSE 9.4%)

    # ---- IL-17 inhibitors ----
    emax_secukinumab <- 5.98
    label("Maximum secukinumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, secukinumab Emax = 5.98 (RSE 5.1%)
    led50_secukinumab <- log(79.2)
    label("Log secukinumab ED50 (log mg per administration); back-transform 79.2 mg") # Table 4, secukinumab ED50 = 79.2 (RSE 15.4%)
    lkdrug_secukinumab <- log(0.467)
    label("Log secukinumab drug-effect onset rate (log 1/week); back-transform 0.467 /week") # Table 4, secukinumab k = 0.467 (RSE 16.1%)

    emax_ixekizumab <- 5.14
    label("Maximum ixekizumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, ixekizumab Emax = 5.14 (RSE 3.3%)
    led50_ixekizumab <- log(7.02)
    label("Log ixekizumab ED50 (log mg per administration); back-transform 7.02 mg") # Table 4, ixekizumab ED50 = 7.02 (RSE 15.8%)
    lkdrug_ixekizumab <- log(1.12)
    label("Log ixekizumab drug-effect onset rate (log 1/week); back-transform 1.12 /week") # Table 4, ixekizumab k = 1.12 (RSE 22.5%)

    emax_brodalumab <- 6.92
    label("Maximum brodalumab effect on the PASI90 logit scale (unitless log-odds)") # Table 4, brodalumab Emax = 6.92 (RSE 7.7%)
    led50_brodalumab <- log(92.5)
    label("Log brodalumab ED50 (log mg per administration); back-transform 92.5 mg") # Table 4, brodalumab ED50 = 92.5 (RSE 20.4%)
    lkdrug_brodalumab <- log(1.19)
    label("Log brodalumab drug-effect onset rate (log 1/week); back-transform 1.19 /week") # Table 4, brodalumab k = 1.19 (RSE 16%)

    # ---- PDE4 inhibitor ----
    emax_apremilast <- fixed(8.8)
    label("Maximum apremilast effect on the PASI90 logit scale, held at the published value (unitless log-odds)") # Table 4, apremilast Emax = '8.8 FIX'
    led50_apremilast <- log(54.1)
    label("Log apremilast ED50 (log mg per administration); back-transform 54.1 mg") # Table 4, apremilast ED50 = 54.1 (RSE 52.3%)
    lkdrug_apremilast <- log(0.0541)
    label("Log apremilast drug-effect onset rate (log 1/week); back-transform 0.0541 /week") # Table 4, apremilast k = 0.0541 (RSE 50.1%)

    # ---- JAK inhibitors ----
    emax_tofacitinib <- 4.76
    label("Maximum tofacitinib effect on the PASI90 logit scale (unitless log-odds)") # Table 4, tofacitinib Emax = 4.76 (RSE 10.7%)
    led50_tofacitinib <- log(3.42)
    label("Log tofacitinib ED50 (log mg per administration); back-transform 3.42 mg") # Table 4, tofacitinib ED50 = 3.42 (RSE 31.9%)
    lkdrug_tofacitinib <- log(0.395)
    label("Log tofacitinib drug-effect onset rate (log 1/week); back-transform 0.395 /week") # Table 4, tofacitinib k = 0.395 (RSE 62.8%)

    emax_baricitinib <- 4.1
    label("Maximum baricitinib effect on the PASI90 logit scale (unitless log-odds)") # Table 4, baricitinib Emax = 4.1 (RSE 4.4%)
    led50_baricitinib <- log(2.9)
    label("Log baricitinib ED50 (log mg per day); back-transform 2.9 mg") # Table 4, baricitinib ED50 = 2.9 (RSE 14.1%)
    lkdrug_baricitinib <- log(0.481)
    label("Log baricitinib drug-effect onset rate (log 1/week); back-transform 0.481 /week") # Table 4, baricitinib k = 0.481 (RSE 24.7%)

    # ---- Dihydrofolate reductase inhibitor (ED50 '0 FIX') ----
    emax_methotrexate <- 2.56
    label("Maximum methotrexate effect on the PASI90 logit scale, reached at any positive dose (unitless log-odds)") # Table 4, methotrexate Emax = 2.56 (RSE 5.6%); ED50 '0 FIX'
    lkdrug_methotrexate <- log(0.19)
    label("Log methotrexate drug-effect onset rate (log 1/week); back-transform 0.19 /week") # Table 4, methotrexate k = 0.19 (RSE 32.6%)

    # ========================================================================
    # BETWEEN-STUDY RANDOM EFFECT on the placebo onset rate kpbo, as printed
    # in Equation 4 (exp(eta) inside the exponent). Supplementary Table S2
    # labels the row 'omega(A), %'; the equation placement is used for the
    # reason given in He_2021_psoriasis_pasi75_mbma (the paper's VPC rules
    # out a 26% CV effect on the asymptote).
    #
    # Supplementary Table S2, 'omega(A), %' = 26 (RSE 24.9%), a CV%:
    # omega^2 = log(1 + 0.26^2) = 0.06541.
    # ========================================================================
    eta_study_lkpbo ~ 0.06541

    # ---- Residual error (Equations 7-8) ----
    addSd_prob_pasi90 <- 1.33
    label("Multiplier on the binomial standard error sqrt(P*(1-P)/N_ARM) giving the residual SD of the arm PASI90 proportion (unitless)")
    # Supplementary Table S2, 'sigma' = 1.33 (RSE 9.7%), read as the SD.
  })

  model({
    # Equation 9: power model on the placebo asymptote, centred at 90 kg.
    wtAsym <- (WT / 90)^e_wt_asym_pbo

    # Equation 4: placebo component; the between-study effect scales the
    # onset rate inside the exponent.
    kpbo <- exp(lkpbo + eta_study_lkpbo)
    e0 <- bsl_pbo + asym_pbo * wtAsym * (1 - exp(-kpbo * time))

    # Equation 5 with c = 1, one term per drug; a zero dose makes the drug's
    # term exactly zero.
    edAdalimumab <- emax_adalimumab * (1 - exp(-exp(lkdrug_adalimumab) * time)) *
      CONMED_ADALIMUMAB_DOSE / (CONMED_ADALIMUMAB_DOSE + exp(led50_adalimumab))
    edInfliximab <- emax_infliximab * (1 - exp(-exp(lkdrug_infliximab) * time)) *
      CONMED_INFLIXIMAB_DOSE / (CONMED_INFLIXIMAB_DOSE + exp(led50_infliximab))
    edEtanercept <- emax_etanercept * (1 - exp(-exp(lkdrug_etanercept) * time)) *
      CONMED_ETANERCEPT_DOSE / (CONMED_ETANERCEPT_DOSE + exp(led50_etanercept))
    edCertolizumab <- emax_certolizumab * (1 - exp(-exp(lkdrug_certolizumab) * time)) *
      CONMED_CERTOLIZUMAB_DOSE / (CONMED_CERTOLIZUMAB_DOSE + exp(led50_certolizumab))
    # ED50 '0 FIX': Dose/(Dose + 0) = 1 for any positive dose.
    edUstekinumab <- emax_ustekinumab * (1 - exp(-exp(lkdrug_ustekinumab) * time)) *
      (CONMED_USTEKINUMAB_DOSE > 0)
    edBriakinumab <- emax_briakinumab * (1 - exp(-exp(lkdrug_briakinumab) * time)) *
      CONMED_BRIAKINUMAB_DOSE / (CONMED_BRIAKINUMAB_DOSE + exp(led50_briakinumab))
    edGuselkumab <- emax_guselkumab * (1 - exp(-exp(lkdrug_guselkumab) * time)) *
      CONMED_GUSELKUMAB_DOSE / (CONMED_GUSELKUMAB_DOSE + exp(led50_guselkumab))
    edTildrakizumab <- emax_tildrakizumab * (1 - exp(-exp(lkdrug_tildrakizumab) * time)) *
      CONMED_TILDRAKIZUMAB_DOSE / (CONMED_TILDRAKIZUMAB_DOSE + exp(led50_tildrakizumab))
    edRisankizumab <- emax_risankizumab * (1 - exp(-exp(lkdrug_risankizumab) * time)) *
      CONMED_RISANKIZUMAB_DOSE / (CONMED_RISANKIZUMAB_DOSE + exp(led50_risankizumab))
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
    # ED50 '0 FIX' for methotrexate as well.
    edMethotrexate <- emax_methotrexate * (1 - exp(-exp(lkdrug_methotrexate) * time)) *
      (CONMED_MTX_DOSE > 0)

    edrug <- edAdalimumab + edInfliximab + edEtanercept + edCertolizumab +
      edUstekinumab + edBriakinumab + edGuselkumab + edTildrakizumab +
      edRisankizumab + edSecukinumab + edIxekizumab + edBrodalumab +
      edApremilast + edTofacitinib + edBaricitinib + edMethotrexate

    # Equations 2-3: inverse logit of the summed placebo and drug components.
    lp_pasi90 <- e0 + edrug
    prob_pasi90 <- expit(lp_pasi90)

    # Equations 7-8: residual SD is sigma times the binomial standard error.
    sdArm <- addSd_prob_pasi90 * sqrt(prob_pasi90 * (1 - prob_pasi90) / N_ARM)
    prob_pasi90 ~ add(sdArm)
  })
}
