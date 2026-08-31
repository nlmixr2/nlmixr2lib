Mandema_2011_anticoagulants_mbma <- function() {
  description <- "MBMA. Joint dose-response meta-analysis of efficacy and bleeding for 23 anticoagulants (7 drug classes) used for venous thromboembolism (VTE) prophylaxis after total hip or total knee replacement surgery. Fit to per-arm event counts from 89 randomized controlled trials (92,543 patients) by binomial nonlinear regression in S-PLUS 6.2. Efficacy follows a sigmoid Emax dose-response on the logit scale (a dose-dependent REDUCTION in event log-odds) for three endpoints (clinical pulmonary embolism, major VTE, total VTE); bleeding follows a linear dose-response on the logit scale (a dose-dependent INCREASE in event log-odds) for three endpoints (major bleeding, major + clinically-relevant-non-major bleeding, total bleeding). Each drug carries its own ED50 for major VTE; the total-VTE ED50 ratio and the bleeding slope ratio (EDbld/ED50) are drug-CLASS properties, which is what makes the relative therapeutic index a class property. Per-arm daily dose enters through one CONMED_<drug>_DOSE column per drug (0 outside that drug's arm; all zero = placebo). Outputs are the six per-arm event probabilities, the corresponding six log odds ratios versus the trial control response, and the enoxaparin-equivalent dose. Suitable simulation scope is per-arm (study-arm-mean) event probability, NOT individual-subject events or concentrations; there is no PK layer and no between-subject variability. The six placebo-response intercepts are digitized from the published Figures 1 and 2 because the source estimated a separate nuisance intercept for each of the 89 trials and reports none of them."

  reference <- paste(
    "Mandema JW, Boyd RA, DiCarlo LA.",
    "Therapeutic index of anticoagulants for prevention of venous thromboembolism",
    "following orthopedic surgery: a dose-response meta-analysis.",
    "Clin Pharmacol Ther. 2011;90(6):820-827.",
    "doi:10.1038/clpt.2011.232.",
    "Dose-response equations and all parameter estimates are in the",
    "Supplementary Data section 'Parameter estimates [95% CI] of the joint",
    "analysis of all efficacy and bleeding endpoints'.",
    sep = " "
  )

  vignette <- "Mandema_2011_anticoagulants_mbma"

  # Algebraic MBMA dose-response model: no rxode2 dose events are consumed (the
  # per-arm daily dose enters through the CONMED_<drug>_DOSE covariate columns)
  # and the outputs are dimensionless event probabilities in [0, 1] rather than
  # drug concentrations. The `units` entries below follow the placeholder
  # convention already used by Vargo_2014_statins_ezetimibe_mbma so that
  # checkModelConventions() sees a parseable dosing / concentration pair.
  units <- list(
    time          = "day (placeholder; the model is a time-independent per-arm dose-response over the trial's prophylaxis period, not a time course)",
    dosing        = "mg/day (per-arm TOTAL DAILY dose supplied through the CONMED_<drug>_DOSE covariate columns, NOT as rxode2 dose events. Units are per-drug: mg/day for most agents, IU/kg/day for ardeparin / nadroparin / tinzaparin, kIU/day for bemiparin / dalteparin / heparin / reviparin, and achieved INR for warfarin -- see each covariateData entry)",
    concentration = "%/arm (per-arm event probability expressed as a percentage of patients in the study arm; the outputs p_* are on the 0-1 fraction scale. Output is NOT a drug concentration; the slash satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    # ---- Enoxaparin (its own drug class; the reference class) ----------------
    CONMED_ENOXAPARIN_DOSE = list(
      description        = "Per-arm total daily enoxaparin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an enoxaparin arm. Enoxaparin was the most frequent active control (39% of trials) and is treated by the source as its own drug class, separate from the other LMWHs, so that all relative effects are estimated against it. Reported daily-dose range across the 54 enoxaparin trials is 10-60 mg (Mandema 2011 Table 1); the approved regimens pooled here are 40 mg q.d. (Europe), 30 mg Q12H (North America, i.e. 60 mg/day) and 20 mg b.i.d. (Japan, i.e. 40 mg/day).",
      source_name        = "Dose (Mandema 2011 Supplementary Data dose-response equations; enoxaparin arms of Table 1)"
    ),
    # ---- LMWHs other than enoxaparin ----------------------------------------
    CONMED_ARDEPARIN_DOSE = list(
      description        = "Per-arm total daily ardeparin dose.",
      units              = "IU/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an ardeparin arm. Weight-normalized anti-Xa activity, NOT a mass dose; the matching ED50 (101 IU/kg/day) carries the same units. Daily-dose range 50-100 IU/kg across 3 trials (Mandema 2011 Table 1). Drug class: LMWH other than enoxaparin.",
      source_name        = "Dose (Mandema 2011 Table 1 ardeparin row; Supplementary Data ED50 Ardeparin)"
    ),
    CONMED_BEMIPARIN_DOSE = list(
      description        = "Per-arm total daily bemiparin dose.",
      units              = "kIU/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a bemiparin arm. Thousands of anti-Xa IU per day; the matching ED50 (3.62 kIU/day) carries the same units. Daily-dose range 3.5-3.5 kIU across 2 trials (Mandema 2011 Table 1). Drug class: LMWH other than enoxaparin.",
      source_name        = "Dose (Mandema 2011 Table 1 bemiparin row; Supplementary Data ED50 Bemiparin)"
    ),
    CONMED_DALTEPARIN_DOSE = list(
      description        = "Per-arm total daily dalteparin dose.",
      units              = "kIU/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a dalteparin arm. Thousands of anti-Xa IU per day; the matching ED50 (5.57 kIU/day) carries the same units. Daily-dose range 5-5 kIU across 8 trials (Mandema 2011 Table 1). Drug class: LMWH other than enoxaparin.",
      source_name        = "Dose (Mandema 2011 Table 1 dalteparin row; Supplementary Data ED50 Dalteparin)"
    ),
    CONMED_NADROPARIN_DOSE = list(
      description        = "Per-arm total daily nadroparin dose.",
      units              = "IU/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a nadroparin arm. Weight-normalized anti-Xa activity; the matching ED50 (78.4 IU/kg/day) carries the same units. Daily-dose range 43-60 IU/kg across 5 trials (Mandema 2011 Table 1). Drug class: LMWH other than enoxaparin.",
      source_name        = "Dose (Mandema 2011 Table 1 nadroparin row; Supplementary Data ED50 Nadroparin)"
    ),
    CONMED_REVIPARIN_DOSE = list(
      description        = "Per-arm total daily reviparin dose.",
      units              = "kIU/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a reviparin arm. Thousands of anti-Xa IU per day; the matching ED50 (6.29 kIU/day) carries the same units. Daily-dose range 2.8-4.2 kIU across 3 trials (Mandema 2011 Table 1). Drug class: LMWH other than enoxaparin.",
      source_name        = "Dose (Mandema 2011 Table 1 reviparin row; Supplementary Data ED50 Reviparin)"
    ),
    CONMED_TINZAPARIN_DOSE = list(
      description        = "Per-arm total daily tinzaparin dose.",
      units              = "IU/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a tinzaparin arm. Weight-normalized anti-Xa activity; the matching ED50 (105 IU/kg/day) carries the same units. Daily-dose range 50-75 IU/kg across 4 trials (Mandema 2011 Table 1). Drug class: LMWH other than enoxaparin.",
      source_name        = "Dose (Mandema 2011 Table 1 tinzaparin row; Supplementary Data ED50 Tinzaparin)"
    ),
    # ---- Direct FXa inhibitors ----------------------------------------------
    CONMED_APIXABAN_DOSE = list(
      description        = "Per-arm total daily apixaban dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an apixaban arm. Daily-dose range 5-20 mg across 4 trials, 6,686 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 apixaban row; Supplementary Data ED50 Apixaban)"
    ),
    CONMED_BETRIXABAN_DOSE = list(
      description        = "Per-arm total daily betrixaban dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a betrixaban arm. Daily-dose range 30-80 mg in 1 trial, 171 patients (Mandema 2011 Table 1); the ED50 is correspondingly imprecise (214 mg/day, 95% CI 41.6-1100). Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 betrixaban row; Supplementary Data ED50 Betrixaban)"
    ),
    CONMED_EDOXABAN_DOSE = list(
      description        = "Per-arm total daily edoxaban dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an edoxaban arm. Daily-dose range 5-90 mg across 2 trials (Mandema 2011 Table 1); appears in the source's Methods search list under its development code DU-176b. Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 edoxaban row; Supplementary Data ED50 Edoxaban)"
    ),
    CONMED_RAZAXABAN_DOSE = list(
      description        = "Per-arm total daily razaxaban dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a razaxaban arm. Daily-dose range 50-200 mg in 1 trial, 506 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 razaxaban row; Supplementary Data ED50 Razaxaban)"
    ),
    CONMED_RIVAROXABAN_DOSE = list(
      description        = "Per-arm total daily rivaroxaban dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a rivaroxaban arm. Daily-dose range 5-60 mg across 7 trials, 7,187 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor. Distinct from the canonical DOSE_RIV_MGKG (per-administration weight-normalized rivaroxaban dose in a paediatric popPK model); this column is a per-arm TOTAL DAILY mass dose in an MBMA.",
      source_name        = "Dose (Mandema 2011 Table 1 rivaroxaban row; Supplementary Data ED50 Rivaroxaban)"
    ),
    CONMED_LY517717_DOSE = list(
      description        = "Per-arm total daily LY517717 dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an LY517717 arm. Investigational oral direct FXa inhibitor with no INN; the development code is used as the column suffix. Daily-dose range 25-150 mg in 1 trial, 417 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 LY517717 row; Supplementary Data ED50 LY517717)"
    ),
    CONMED_PD0348292_DOSE = list(
      description        = "Per-arm total daily PD0348292 dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a PD0348292 arm. Investigational oral direct FXa inhibitor with no INN; the development code is used as the column suffix. Daily-dose range 0.1-10 mg in 1 trial, 992 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 PD0348292 row; Supplementary Data ED50 PD0348292)"
    ),
    CONMED_YM150_DOSE = list(
      description        = "Per-arm total daily YM150 dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a YM150 arm. Investigational oral direct FXa inhibitor (darexaban) identified in the source only by its development code. Daily-dose range 3-60 mg in 1 trial, 138 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 YM150 row; Supplementary Data ED50 YM150)"
    ),
    CONMED_AVE5026_DOSE = list(
      description        = "Per-arm total daily AVE5026 dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an AVE5026 arm. Investigational ultra-low-molecular-weight heparin (semuloparin) that the source classifies as a direct FXa inhibitor and identifies only by its development code. Daily-dose range 5-60 mg in 1 trial, 559 patients (Mandema 2011 Table 1). Drug class: direct FXa inhibitor per the source's Methods inclusion-criteria list.",
      source_name        = "Dose (Mandema 2011 Table 1 AVE5026 row; Supplementary Data ED50 AVE5026)"
    ),
    # ---- Indirect FXa inhibitor ---------------------------------------------
    CONMED_FONDAPARINUX_DOSE = list(
      description        = "Per-arm total daily fondaparinux dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a fondaparinux arm. Daily-dose range 0.8-8 mg across 6 trials, 6,121 patients (Mandema 2011 Table 1). Fondaparinux is the only member of the indirect-FXa-inhibitor class, so that class's parameters are identified by fondaparinux alone.",
      source_name        = "Dose (Mandema 2011 Table 1 fondaparinux row; Supplementary Data ED50 Fondaparinux)"
    ),
    # ---- Univalent thrombin inhibitors --------------------------------------
    CONMED_DABIGATRAN_DOSE = list(
      description        = "Per-arm total daily dabigatran dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a dabigatran arm. Dose is of dabigatran etexilate as reported by the trials. Daily-dose range 25-600 mg across 6 trials, 7,653 patients (Mandema 2011 Table 1). Drug class: univalent (direct) thrombin inhibitor.",
      source_name        = "Dose (Mandema 2011 Table 1 dabigatran row; Supplementary Data ED50 Dabigatran)"
    ),
    CONMED_XIMELAGATRAN_DOSE = list(
      description        = "Per-arm total daily ximelagatran dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a ximelagatran arm. Daily-dose range 12-72 mg across 10 trials, 9,249 patients (Mandema 2011 Table 1). Drug class: univalent (direct) thrombin inhibitor. Arms in which ximelagatran was combined with subcutaneous melagatran STARTED BEFORE SURGERY are flagged by FORM_XIMELAGATRAN_PRESURG and use a proportionally lower ED50; the source treated that regimen as a separate drug because its potency differed significantly (P < 0.001).",
      source_name        = "Dose (Mandema 2011 Table 1 ximelagatran row; Supplementary Data ED50 Ximelagatran)"
    ),
    FORM_XIMELAGATRAN_PRESURG = list(
      description        = "Indicator that a ximelagatran arm used the pre-surgery regimen: ximelagatran given in combination with subcutaneous melagatran started BEFORE surgery (1) versus ximelagatran administered alone (0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = ximelagatran administered alone (the reference regimen whose ED50 is the ximelagatran ED50 of 78.8 mg/day).",
      notes              = "Study-arm-level regimen indicator, not a subject-level covariate. Multiplies the ximelagatran ED50 by 0.365 when set to 1, i.e. the pre-surgery melagatran combination is about 2.7-fold more potent. The source states: 'Ximelagatran given in combination with subcutaneous melagatran started before surgery was treated as a separate drug from ximelagatran when administered alone because of a statistically significant difference in potency between these two treatments (P < 0.001).' Follows the FORM_<drug>_<regimen> family established by FORM_FLV_BID_XR / FORM_LOV_BID_XR, which likewise encode an ED50-ratio regimen switch in an MBMA dose-response model. Has no effect when CONMED_XIMELAGATRAN_DOSE is 0.",
      source_name        = "'Ximelagatran started before surgery / ED50 Ximelagatran' (Mandema 2011 Supplementary Data parameter table)"
    ),
    # ---- Bivalent thrombin inhibitor ----------------------------------------
    CONMED_DESIRUDIN_DOSE = list(
      description        = "Per-arm total daily desirudin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a desirudin arm. Daily-dose range 20-40 mg across 3 trials, 2,095 patients (Mandema 2011 Table 1). Desirudin is the only member of the bivalent-thrombin-inhibitor class. The source excluded this class from the Table 2 pairwise therapeutic-index comparison 'because only limited data are available'; its EDbld ratio (3.72) is estimated but very imprecise (95% CI 0.652-21.3).",
      source_name        = "Dose (Mandema 2011 Table 1 desirudin row; Supplementary Data ED50 Desirudin)"
    ),
    # ---- Synthetic mixed FXa and thrombin inhibitor --------------------------
    CONMED_SR123781A_DOSE = list(
      description        = "Per-arm total daily SR123781A dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an SR123781A arm. Investigational synthetic mixed FXa and thrombin inhibitor with no INN; the development code is used as the column suffix. Daily-dose range 0.2-4 mg in 1 trial, 843 patients (Mandema 2011 Table 1). Only member of its class, and excluded from the Table 2 pairwise therapeutic-index comparison for the same limited-data reason as desirudin.",
      source_name        = "Dose (Mandema 2011 Table 1 SR123781A row; Supplementary Data ED50 SR123781A)"
    ),
    # ---- Unfractionated heparin ---------------------------------------------
    CONMED_HEPARIN_DOSE = list(
      description        = "Per-arm total daily unfractionated heparin dose.",
      units              = "kIU/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an unfractionated-heparin arm. Thousands of IU per day; the matching ED50 (57.9 kIU/day) carries the same units. Daily-dose range 10-15 kIU across 18 trials, 3,493 patients (Mandema 2011 Table 1). Heparin is its own drug class and has by far the narrowest therapeutic index in the analysis. Distinct from the canonical DOSE_UFH_UH (a concomitant continuous-infusion heparin dose RATE in U/h); this column is a per-arm total DAILY prophylactic dose in an MBMA.",
      source_name        = "Dose (Mandema 2011 Table 1 heparin row; Supplementary Data ED50 Heparin)"
    ),
    # ---- Warfarin ------------------------------------------------------------
    CONMED_WARFARIN_DOSE = list(
      description        = "Per-arm achieved international normalized ratio (INR) for a warfarin arm, used as the warfarin dose regressor.",
      units              = "INR",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a warfarin arm. Warfarin is dosed to an INR target rather than to a fixed mass dose, so the source uses the arm's achieved INR as the dose variable: Mandema 2011 Table 1 reports the warfarin 'daily dose range' as 2.2-2.5 INR and the Supplementary Data reports 'ED50 Warfarin (INR)' = 4.05. Feeding a milligram warfarin dose into this column would be a units error. Distinct from the canonical INR_BASE, which is a subject-level PRE-medication baseline INR used as a PK/PD model covariate; this column is a per-arm ON-TREATMENT achieved INR used as an MBMA dose surrogate. NOTE: the source could not estimate a bleeding slope (EDbld) for warfarin, so a warfarin arm has NO bleeding dose-response -- the model reports this through the warfarin_bleeding_undefined output.",
      source_name        = "Dose in INR (Mandema 2011 Table 1 warfarin row 'Daily dose range 2.2-2.5 INR'; Supplementary Data 'ED50 Warfarin (INR)')"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Per-arm mean patient age. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Mandema 2011 Methods lists age among the extracted patient-population characteristics tested on the dose-response parameters. Results: 'None of the covariates significantly influenced the ED50 or Emax of the dose-response relationship' and 'None of the covariates significantly impacted the slope of the bleeding dose-response relationship.' No point estimate is reported for any screened covariate, so no effect parameter is included in ini(). Per-arm mean ages are tabulated in the Supplementary Data trial-overview table.",
      source_name        = "mean age (Mandema 2011 Methods covariate list; Supplementary Data trial-overview table)"
    ),
    WT = list(
      description        = "Per-arm mean patient body weight. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods ('patient-population characteristics (age, weight, gender, type of anesthesia, type of surgery ...)'); not retained and no point estimate reported. Note that three LMWHs in this model are dosed per kilogram (ardeparin, nadroparin, tinzaparin), so body weight enters those arms through the IU/kg dose units rather than as a covariate effect.",
      source_name        = "weight (Mandema 2011 Methods covariate list)"
    ),
    SEXF = list(
      description        = "Per-arm proportion female. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods ('gender'); not retained and no point estimate reported.",
      source_name        = "gender (Mandema 2011 Methods covariate list)"
    ),
    SURG_HIP = list(
      description        = "Per-arm proportion of patients undergoing total hip replacement rather than total knee replacement. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods ('type of surgery (hip vs. knee)') and tabulated per arm in the Supplementary Data trial-overview table as '% hip replacement'. NOT retained on the dose-response parameters -- the RELATIVE treatment effect was the same for hip and knee surgery. Surgery type does strongly affect the ABSOLUTE event rate, which the source absorbed into the per-trial placebo intercept rather than into the dose-response: 'knee surgery showing a higher absolute total VTE event rate (28.6% for enoxaparin) than hip surgery (13.5% for enoxaparin)'. A user who wants a knee-specific or hip-specific absolute prediction should therefore shift the e0_* intercepts, not add a covariate effect.",
      source_name        = "% hip replacement (Mandema 2011 Methods covariate list; Supplementary Data trial-overview table)"
    ),
    REGION = list(
      description        = "Primary geographic location of the trial (Asia, Australia, Europe, or North America). Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods trial-characteristics list; not retained and no point estimate reported. Regional counts per drug are in Mandema 2011 Table 1 and per trial in the Supplementary Data trial-overview table. Region matters clinically because the approved enoxaparin regimen differs by region (40 mg q.d. in Europe, 30 mg Q12H in North America, 20 mg b.i.d. in Japan), but the source handled that through the dose axis rather than through a covariate effect.",
      source_name        = "primary geographic location (Mandema 2011 Methods covariate list)"
    ),
    TRT_START_REL_SURGERY = list(
      description        = "Timing of treatment start relative to surgery. Screened as a dose-response covariate but NOT retained in the final model, EXCEPT for the ximelagatran pre-surgery regimen which was carried as a separate potency.",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods treatment-covariate list ('regimen, time of treatment start relative to surgery, treatment duration'). Not retained as a general covariate effect; no point estimate is reported. The one place the source did retain a start-timing effect is the ximelagatran + subcutaneous melagatran started-before-surgery regimen, which is modelled as an ED50 ratio through FORM_XIMELAGATRAN_PRESURG rather than as a continuous covariate.",
      source_name        = "time of treatment start relative to surgery (Mandema 2011 Methods covariate list)"
    ),
    TRT_DURATION = list(
      description        = "Treatment (prophylaxis) duration. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods treatment-covariate list; not retained and no point estimate reported. Trials of extended prophylaxis and trials whose treatment duration differed across arms were excluded at the screening stage, which limits the range over which duration could have acted.",
      source_name        = "treatment duration (Mandema 2011 Methods covariate list)"
    ),
    VENOGRAPHY_BILATERAL = list(
      description        = "Method of venography used for VTE ascertainment (bilateral versus unilateral). Screened as a trial-level dose-response covariate but NOT retained in the final model.",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods trial-characteristics list; not retained and no point estimate reported. Only trials using mandatory venography at the end of the evaluation period were eligible, so this covariate distinguishes ascertainment thoroughness rather than presence.",
      source_name        = "method of venography (unilateral vs. bilateral) (Mandema 2011 Methods covariate list)"
    ),
    TRIAL_START_YEAR = list(
      description        = "Calendar year in which the trial started. Screened as a trial-level dose-response covariate but NOT retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods trial-characteristics list; not retained and no point estimate reported.",
      source_name        = "year of trial start (Mandema 2011 Methods covariate list)"
    ),
    ANESTHESIA_TYPE = list(
      description        = "Type of anesthesia used for the surgery. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods patient-population covariate list ('type of anesthesia'); not retained and no point estimate reported.",
      source_name        = "type of anesthesia (Mandema 2011 Methods covariate list)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 92543L,
    n_studies        = 89L,
    age_range        = "Adults undergoing elective total hip or total knee replacement; per-arm mean ages (roughly 60-70 years for most arms) are tabulated in the Supplementary Data trial-overview table. Age was screened as a covariate and not retained.",
    weight_range     = "Not tabulated in the main paper. Body weight was screened as a covariate and not retained; three LMWHs (ardeparin, nadroparin, tinzaparin) are dosed per kilogram so weight enters those arms through the IU/kg dose units.",
    sex_female_pct   = NA_real_,
    disease_state    = "Prophylaxis of venous thromboembolism after elective total hip replacement (THR) or total knee replacement (TKR) surgery. 61 trials enrolled THR patients and 41 TKR (some trials both); trials required at least 75% of patients to have had hip or knee replacement. Efficacy endpoints: clinical pulmonary embolism (PE), major VTE (proximal DVT + clinical PE +/- all-cause or VTE-related death) and total VTE (distal DVT + major VTE), all ascertained by mandatory end-of-period venography. Safety endpoints: major bleeding (ISTH criteria adapted for surgery), major + clinically-relevant-non-major (CRNM) bleeding, and total bleeding.",
    dose_range       = "23 anticoagulants across 7 drug classes; per-drug daily-dose ranges are in Mandema 2011 Table 1 (e.g. enoxaparin 10-60 mg, dabigatran 25-600 mg, apixaban 5-20 mg, rivaroxaban 5-60 mg, fondaparinux 0.8-8 mg, heparin 10-15 kIU, warfarin INR 2.2-2.5). 15 trials (1,286 patients) were placebo-controlled.",
    regions          = "6 trials in Asia, 47 in Europe (including Australia), and 35 in North America (Mandema 2011 Table 1).",
    n_trials_efficacy = 87L,
    n_trials_bleeding = 74L,
    drug_classes     = "Enoxaparin (treated as its own class, and the reference for all relative effects); other LMWHs (ardeparin, bemiparin, dalteparin, nadroparin, reviparin, tinzaparin); direct FXa inhibitors (apixaban, betrixaban, edoxaban, razaxaban, rivaroxaban, LY517717, PD0348292, YM150, AVE5026); indirect FXa inhibitors (fondaparinux); univalent thrombin inhibitors (dabigatran, ximelagatran); bivalent thrombin inhibitors (desirudin); synthetic mixed FXa and thrombin inhibitors (SR123781A); unfractionated heparin; and warfarin. Certoparin appears in the Methods search list but no certoparin trial met the inclusion criteria, so it has no ED50 and is not a column in this model.",
    notes            = "Summary-level MBMA: the modelled observations were per-arm event COUNTS (number of patients with the event out of the arm sample size), assumed binomial, with a separate placebo-response intercept estimated for EACH of the 89 trials. Because of those per-trial intercepts the identified response variable is the log odds ratio between the active and control arms, not the absolute event rate; the source reports none of the 89 intercepts. Trial-specific random effects on Emax and ED50 were tested and rejected (estimated variances close to zero), so the model carries NO random effects and NO between-subject variability. Within-arm correlation across the multiple endpoints measured in the same arm was handled by a compound-symmetry structure that is a property of the estimation, not of the predictive model, and is therefore not encoded. Absolute event rates varied enormously across trials (1.7-46% for total VTE and 0-5.3% for major bleeding across all enoxaparin doses), so the e0_* intercepts in this model reproduce the single reference trial drawn in Figures 1 and 2, not any particular study."
  )

  ini({
    # ========================================================================
    # Placebo-response intercepts, one per endpoint, on the logit scale.
    #
    # PROVENANCE -- these six values are the ONLY parameters in this file that
    # are not printed in the source. Mandema 2011 estimated a separate
    # intercept for each of the 89 trials ("The trial-specific placebo response
    # (i.e., mean event rate in each trial) accounted for the trial-to-trial
    # variability in the overall event incidence. Because a different intercept
    # is estimated for each trial, the primary response variable is the log of
    # the odds ratio between the active and the control arms") and reports none
    # of them. Figures 1 and 2 each plot the PREDICTED dose-response curve for
    # enoxaparin against enoxaparin-equivalent dose for their three endpoints,
    # so the curve's dose -> 0 limit is the reference placebo response.
    #
    # Each value below was recovered by digitizing the published solid curve
    # (200 dpi render, axis calibration from the printed tick marks) and
    # solving E0 = logit(P(dose)) +/- g(dose) at every extracted pixel column
    # using the source's own published parameters. The recovered E0 was
    # essentially constant along each curve -- interquartile range 0.009-0.023
    # on the logit scale, i.e. under +/- 0.15 percentage points of incidence --
    # which is simultaneously a strong confirmation that the equations and
    # parameter values encoded below reproduce the published figures.
    #
    # These intercepts describe the single reference trial drawn in the
    # figures. Real trials varied from 1.7% to 46% for total VTE. Shift these
    # values, do NOT add a covariate effect, to move to a different absolute
    # event-rate setting (e.g. knee vs hip surgery).
    # ========================================================================
    e0_pe <- -4.768
    label("Control-arm log-odds of clinical pulmonary embolism (logit scale; 0.84% incidence)")  # Mandema 2011 Figure 1 PE panel, digitized dose -> 0 asymptote of the predicted enoxaparin curve

    e0_majorvte <- -1.964
    label("Control-arm log-odds of major VTE (logit scale; 12.3% incidence)")  # Mandema 2011 Figure 1 major VTE panel, digitized dose -> 0 asymptote of the predicted enoxaparin curve

    e0_totalvte <- -0.644
    label("Control-arm log-odds of total VTE (logit scale; 34.4% incidence)")  # Mandema 2011 Figure 1 VTE panel, digitized dose -> 0 asymptote of the predicted enoxaparin curve

    e0_majorbld <- -4.976
    label("Control-arm log-odds of major bleeding (logit scale; 0.69% incidence)")  # Mandema 2011 Figure 2 major bleeding panel, digitized dose -> 0 asymptote of the predicted enoxaparin curve

    e0_crnmbld <- -3.579
    label("Control-arm log-odds of major + CRNM bleeding (logit scale; 2.7% incidence)")  # Mandema 2011 Figure 2 major + CRNM bleeding panel, digitized dose -> 0 asymptote of the predicted enoxaparin curve

    e0_totalbld <- -2.953
    label("Control-arm log-odds of total bleeding (logit scale; 5.0% incidence)")  # Mandema 2011 Figure 2 total bleeding panel, digitized dose -> 0 asymptote of the predicted enoxaparin curve

    # ========================================================================
    # Efficacy dose-response shape (sigmoid Emax on the logit scale).
    #
    #   g_efficacy = (Emax + Emax_endpoint) * Dose^n
    #                / (Dose^n + (ED50_drug * ED50_endpoint * ED50_class_endpoint)^n)
    #
    # and logit(P_event) = E0 - g_efficacy, so a larger g is a LARGER reduction
    # in VTE risk. Emax and ED50 are shared between major VTE and PE ("A model
    # that assumed that the Emax and ED50 values for PE were the same as those
    # for major VTE for each drug best described the PE data"); total VTE gets
    # its own Emax shift and ED50 ratios.
    # ========================================================================
    emax <- 3.18
    label("Maximal reduction in event log-odds for major VTE and PE (logit units)")  # Mandema 2011 Supplementary Data parameter table, Emax = 3.18 [2.59 to 3.76]

    emax_totalvte <- -0.744
    label("Additive shift in Emax for total VTE relative to major VTE (logit units)")  # Mandema 2011 Supplementary Data parameter table, EmaxtotalVTE = -0.744 [-1.18 to -0.303]

    lhill <- log(1.33)
    label("Hill exponent of the efficacy dose-response (unitless)")  # Mandema 2011 Supplementary Data parameter table, n = 1.33 [1.09 to 1.61]

    # ========================================================================
    # ED50 for MAJOR VTE, one per drug. This is the drug's potency anchor: the
    # total-VTE ED50 and the bleeding slope are both expressed relative to it.
    # Units follow the drug (mg/day unless the label says otherwise) and match
    # the corresponding CONMED_<drug>_DOSE column.
    # All values: Mandema 2011 Supplementary Data parameter table.
    # ========================================================================
    led50_enoxaparin <- log(61.6)
    label("Enoxaparin ED50 for major VTE (mg/day)")  # ED50 Enoxaparin 61.6 [44.9 to 84.5]

    led50_ardeparin <- log(101)
    label("Ardeparin ED50 for major VTE (IU/kg/day)")  # ED50 Ardeparin 101 [66.6 to 154]

    led50_bemiparin <- log(3.62)
    label("Bemiparin ED50 for major VTE (kIU/day)")  # ED50 Bemiparin 3.62 [2.12 to 6.18]

    led50_dalteparin <- log(5.57)
    label("Dalteparin ED50 for major VTE (kIU/day)")  # ED50 Dalteparin 5.57 [3.81 to 8.14]

    led50_nadroparin <- log(78.4)
    label("Nadroparin ED50 for major VTE (IU/kg/day)")  # ED50 Nadroparin 78.4 [50.1 to 123]

    led50_reviparin <- log(6.29)
    label("Reviparin ED50 for major VTE (kIU/day)")  # ED50 Reviparin 6.29 [4.07 to 9.72]

    led50_tinzaparin <- log(105)
    label("Tinzaparin ED50 for major VTE (IU/kg/day)")  # ED50 Tinzaparin 105 [69.9 to 158]

    led50_apixaban <- log(4.31)
    label("Apixaban ED50 for major VTE (mg/day)")  # ED50 Apixaban 4.31 [2.85 to 6.51]

    led50_betrixaban <- log(214)
    label("Betrixaban ED50 for major VTE (mg/day)")  # ED50 Betrixaban 214 [41.6 to 1100]

    led50_edoxaban <- log(15.2)
    label("Edoxaban ED50 for major VTE (mg/day)")  # ED50 Edoxaban 15.2 [8.78 to 26.3]

    led50_razaxaban <- log(18.1)
    label("Razaxaban ED50 for major VTE (mg/day)")  # ED50 Razaxaban 18.1 [8.94 to 36.5]

    led50_rivaroxaban <- log(5.51)
    label("Rivaroxaban ED50 for major VTE (mg/day)")  # ED50 Rivaroxaban 5.51 [3.6 to 8.44]

    led50_ly517717 <- log(225)
    label("LY517717 ED50 for major VTE (mg/day)")  # ED50 LY517717 225 [94.9 to 534]

    led50_pd0348292 <- log(1.98)
    label("PD0348292 ED50 for major VTE (mg/day)")  # ED50 Pd0348292 1.98 [1.07 to 3.66]

    led50_ym150 <- log(26.1)
    label("YM150 ED50 for major VTE (mg/day)")  # ED50 YM150 26.1 [8.95 to 76.3]

    led50_ave5026 <- log(13.1)
    label("AVE5026 ED50 for major VTE (mg/day)")  # ED50 AVE5026 13.1 [7.49 to 23]

    led50_fondaparinux <- log(1.97)
    label("Fondaparinux ED50 for major VTE (mg/day)")  # ED50 Fondaparinux 1.97 [1.22 to 3.19]

    led50_dabigatran <- log(229)
    label("Dabigatran ED50 for major VTE (mg/day)")  # ED50 dabigatran 229 [159 to 330]

    led50_ximelagatran <- log(78.8)
    label("Ximelagatran ED50 for major VTE (mg/day; administered alone)")  # ED50 Ximelagatran 78.8 [54.3 to 114]

    led50_desirudin <- log(21.4)
    label("Desirudin ED50 for major VTE (mg/day)")  # ED50 Desirudin 21.4 [13.5 to 33.8]

    led50_sr123781a <- log(1.41)
    label("SR123781A ED50 for major VTE (mg/day)")  # ED50 SR123781A 1.41 [0.397 to 4.99]

    led50_heparin <- log(57.9)
    label("Unfractionated heparin ED50 for major VTE (kIU/day)")  # ED50 Heparin 57.9 [28.5 to 118]

    led50_warfarin <- log(4.05)
    label("Warfarin ED50 for major VTE (INR)")  # ED50 Warfarin 4.05 [2.71 to 6.04]

    # ------------------------------------------------------------------------
    # Regimen ratio: ximelagatran + subcutaneous melagatran started BEFORE
    # surgery is about 2.7-fold more potent than ximelagatran alone and was
    # treated by the source as a separate drug (P < 0.001).
    # ------------------------------------------------------------------------
    ratio_ed50_ximelagatran_presurg <- 0.365
    label("ED50 ratio for ximelagatran started before surgery relative to ximelagatran alone (unitless)")  # Mandema 2011 Supplementary Data parameter table, ED50 Ximelagatran started before surgery / ED50 Ximelagatran = 0.365 [0.286 to 0.466]

    # ========================================================================
    # Total-VTE ED50 ratios. The total-VTE ED50 of a drug is
    #   ED50_totalVTE = ED50_majorVTE * ratio_ed50_totalvte * ratio_ed50_totalvte_<class>
    # The first factor is the overall (enoxaparin-reference) endpoint ratio;
    # the second is the drug-class interaction, which is the finding that
    # "the difference in drug potency for major VTE vs. total VTE is dependent
    # on the mechanism of action only". Enoxaparin is the reference class and
    # its class ratio is 1 by definition (no parameter; see model()).
    # All values: Mandema 2011 Supplementary Data parameter table.
    # ========================================================================
    ratio_ed50_totalvte <- 0.721
    label("Total-VTE / major-VTE ED50 ratio for enoxaparin (unitless)")  # ED50,totalVTE 0.721 [0.549 to 0.947]

    ratio_ed50_totalvte_lmwh <- 1.26
    label("Additional total-VTE ED50 ratio for LMWHs other than enoxaparin (unitless)")  # ED50,totalVTE LMWH other than enoxaparin 1.26 [1.04 to 1.52]

    ratio_ed50_totalvte_dfxa <- 0.769
    label("Additional total-VTE ED50 ratio for direct FXa inhibitors (unitless)")  # ED50,totalVTE direct FXa inhibitor 0.769 [0.621 to 0.953]

    ratio_ed50_totalvte_ifxa <- 0.499
    label("Additional total-VTE ED50 ratio for indirect FXa inhibitors (unitless)")  # ED50,totalVTE indirect FXa inhibitor 0.499 [0.331 to 0.752]

    ratio_ed50_totalvte_uthr <- 1.28
    label("Additional total-VTE ED50 ratio for univalent thrombin inhibitors (unitless)")  # ED50,totalVTE univalent thrombin inhibitor 1.28 [1.1 to 1.48]

    ratio_ed50_totalvte_bthr <- 1.03
    label("Additional total-VTE ED50 ratio for bivalent thrombin inhibitors (unitless)")  # ED50,totalVTE bivalent thrombin inhibitor 1.03 [0.744 to 1.44]

    ratio_ed50_totalvte_mixed <- 1.97
    label("Additional total-VTE ED50 ratio for mixed FXa and thrombin inhibitors (unitless)")  # ED50,totalVTE mixed FXa and thrombin inhibitors 1.97 [0.626 to 6.2]

    ratio_ed50_totalvte_heparin <- 0.535
    label("Additional total-VTE ED50 ratio for unfractionated heparin (unitless)")  # ED50,totalVTE Heparin 0.535 [0.332 to 0.864]

    ratio_ed50_totalvte_warfarin <- 1.7
    label("Additional total-VTE ED50 ratio for warfarin (unitless)")  # ED50,totalVTE Warfarin 1.7 [1.34 to 2.15]

    # ========================================================================
    # Bleeding dose-response (linear on the logit scale).
    #
    #   g_bleeding = Dose / (ED50_majorVTE * EDbld_class) * (1 + Sc_endpoint)
    #
    # and logit(P_bleed) = E0 + g_bleeding, so a larger g is a LARGER bleeding
    # risk. EDbld is parameterised as a RATIO to the drug's own major-VTE ED50,
    # which is exactly why the relative therapeutic index reduces to a ratio of
    # these class parameters: relative TI = EDbld_class / EDbld_enoxaparin.
    # No EDbld is reported for warfarin ("The TI for warfarin could not be
    # estimated"), so warfarin has no bleeding dose-response in this model.
    # All values: Mandema 2011 Supplementary Data parameter table.
    # ========================================================================
    ratio_edbld_enox <- 1.7
    label("EDbld / ED50 ratio for enoxaparin (unitless)")  # EDbldEnoxaparin 1.7 [1.14 to 2.54]

    ratio_edbld_lmwh <- 1.58
    label("EDbld / ED50 ratio for LMWHs other than enoxaparin (unitless)")  # EDbld LMWH other than enoxaparin 1.58 [0.923 to 2.71]

    ratio_edbld_dfxa <- 3.98
    label("EDbld / ED50 ratio for direct FXa inhibitors (unitless)")  # EDbld direct FXa inhibitor 3.98 [2.51 to 6.29]

    ratio_edbld_ifxa <- 1.28
    label("EDbld / ED50 ratio for indirect FXa inhibitors (unitless)")  # EDbld indirect FXa inhibitor 1.28 [0.762 to 2.16]

    ratio_edbld_uthr <- 1.43
    label("EDbld / ED50 ratio for univalent thrombin inhibitors (unitless)")  # EDbld univalent thrombin inhibitor 1.43 [0.943 to 2.16]

    ratio_edbld_bthr <- 3.72
    label("EDbld / ED50 ratio for bivalent thrombin inhibitors (unitless)")  # EDbld bivalent thrombin inhibitor 3.72 [0.652 to 21.3]

    ratio_edbld_mixed <- 0.841
    label("EDbld / ED50 ratio for mixed FXa and thrombin inhibitors (unitless)")  # EDbld mixed FXa and thrombin inhibitors 0.841 [0.231 to 3.07]

    ratio_edbld_heparin <- 0.312
    label("EDbld / ED50 ratio for unfractionated heparin (unitless)")  # EDbld Heparin 0.312 [0.142 to 0.687]

    # ------------------------------------------------------------------------
    # Bleeding-endpoint slope scalars. Major bleeding is the reference
    # (Sc = 0 by definition; no parameter, see model()). A NEGATIVE value is a
    # SHALLOWER dose-response slope, matching the source's statement that "the
    # slope of the dose-response relationship for total bleeding was smaller
    # than that for major bleeding". The shift is common to all drugs.
    # ------------------------------------------------------------------------
    sc_bld_crnm <- -0.291
    label("Fractional change in the bleeding dose-response slope for major + CRNM bleeding (unitless)")  # Mandema 2011 Supplementary Data parameter table, Scbld,CRNM+major -0.291 [-0.501 to -0.0804]

    sc_bld_total <- -0.369
    label("Fractional change in the bleeding dose-response slope for total bleeding (unitless)")  # Mandema 2011 Supplementary Data parameter table, Scbld,total bleed -0.369 [-0.494 to -0.243]

    # ========================================================================
    # No random effects and no residual error model.
    #
    # Random effects: "Inclusion of a trial-specific random effect on Emax or
    # ED50 did not result in a significant improvement of fit, and the
    # estimated variances were close to zero", so the published final model
    # carries none. The trial-to-trial variability in the ABSOLUTE event rate
    # was absorbed by the 89 per-trial intercepts, not by a random effect.
    #
    # Residual error: the observations were per-arm event COUNTS with a
    # binomial likelihood (N_event,kij ~ binomial(P(event)kij, N_kij)); there
    # is no additive or proportional residual SD to encode. The outputs below
    # are therefore deterministic per-arm probabilities. A user who wants
    # simulated counts should draw rbinom(n = N_arm, prob = p_<endpoint>)
    # downstream, which is what the source's likelihood assumes.
    # ========================================================================
  })

  model({
    hill <- exp(lhill)

    # ---- Per-drug ED50 for major VTE (back-transformed) --------------------
    ed50_enoxaparin   <- exp(led50_enoxaparin)
    ed50_ardeparin    <- exp(led50_ardeparin)
    ed50_bemiparin    <- exp(led50_bemiparin)
    ed50_dalteparin   <- exp(led50_dalteparin)
    ed50_nadroparin   <- exp(led50_nadroparin)
    ed50_reviparin    <- exp(led50_reviparin)
    ed50_tinzaparin   <- exp(led50_tinzaparin)
    ed50_apixaban     <- exp(led50_apixaban)
    ed50_betrixaban   <- exp(led50_betrixaban)
    ed50_edoxaban     <- exp(led50_edoxaban)
    ed50_razaxaban    <- exp(led50_razaxaban)
    ed50_rivaroxaban  <- exp(led50_rivaroxaban)
    ed50_ly517717     <- exp(led50_ly517717)
    ed50_pd0348292    <- exp(led50_pd0348292)
    ed50_ym150        <- exp(led50_ym150)
    ed50_ave5026      <- exp(led50_ave5026)
    ed50_fondaparinux <- exp(led50_fondaparinux)
    ed50_dabigatran   <- exp(led50_dabigatran)
    ed50_desirudin    <- exp(led50_desirudin)
    ed50_sr123781a    <- exp(led50_sr123781a)
    ed50_heparin      <- exp(led50_heparin)
    ed50_warfarin     <- exp(led50_warfarin)

    # Ximelagatran + subcutaneous melagatran started before surgery is a
    # separate potency. FORM_XIMELAGATRAN_PRESURG is binary, so the algebraic
    # form (1 + (ratio - 1) * flag) gives the ratio when the flag is 1 and
    # leaves the ED50 unchanged when it is 0.
    ed50_ximelagatran <- exp(led50_ximelagatran) *
      (1 + (ratio_ed50_ximelagatran_presurg - 1) * FORM_XIMELAGATRAN_PRESURG)

    # ---- Potency-normalised dose, accumulated per DRUG CLASS ----------------
    # r_<class> = Dose / ED50_majorVTE for whichever drug in that class was
    # given. Every trial arm in the source received exactly ONE anticoagulant,
    # so at most one of these nine terms is non-zero and the sums below are
    # exact. A placebo arm sets every CONMED_<drug>_DOSE to 0, which makes
    # every r_<class> zero and collapses the dose-response to the intercept.
    # Coding two drugs simultaneously is outside the source's calibration and
    # would make the model add their normalised doses.
    r_enox <- CONMED_ENOXAPARIN_DOSE / ed50_enoxaparin

    r_lmwh <- CONMED_ARDEPARIN_DOSE  / ed50_ardeparin +
      CONMED_BEMIPARIN_DOSE  / ed50_bemiparin +
      CONMED_DALTEPARIN_DOSE / ed50_dalteparin +
      CONMED_NADROPARIN_DOSE / ed50_nadroparin +
      CONMED_REVIPARIN_DOSE  / ed50_reviparin +
      CONMED_TINZAPARIN_DOSE / ed50_tinzaparin

    r_dfxa <- CONMED_APIXABAN_DOSE    / ed50_apixaban +
      CONMED_BETRIXABAN_DOSE  / ed50_betrixaban +
      CONMED_EDOXABAN_DOSE    / ed50_edoxaban +
      CONMED_RAZAXABAN_DOSE   / ed50_razaxaban +
      CONMED_RIVAROXABAN_DOSE / ed50_rivaroxaban +
      CONMED_LY517717_DOSE    / ed50_ly517717 +
      CONMED_PD0348292_DOSE   / ed50_pd0348292 +
      CONMED_YM150_DOSE       / ed50_ym150 +
      CONMED_AVE5026_DOSE     / ed50_ave5026

    r_ifxa  <- CONMED_FONDAPARINUX_DOSE / ed50_fondaparinux

    r_uthr  <- CONMED_DABIGATRAN_DOSE   / ed50_dabigatran +
      CONMED_XIMELAGATRAN_DOSE / ed50_ximelagatran

    r_bthr  <- CONMED_DESIRUDIN_DOSE  / ed50_desirudin
    r_mixed <- CONMED_SR123781A_DOSE  / ed50_sr123781a
    r_hep   <- CONMED_HEPARIN_DOSE    / ed50_heparin
    r_warf  <- CONMED_WARFARIN_DOSE   / ed50_warfarin

    r_major <- r_enox + r_lmwh + r_dfxa + r_ifxa + r_uthr +
      r_bthr + r_mixed + r_hep + r_warf

    # ---- Efficacy: major VTE and clinical PE -------------------------------
    # Sigmoid Emax in the drug's own major-VTE ED50 units. PE shares Emax and
    # ED50 with major VTE, so the two differ only in their intercept.
    u_major    <- r_major^hill
    g_majorvte <- emax * u_major / (u_major + 1)
    g_pe       <- g_majorvte

    # ---- Efficacy: total VTE -----------------------------------------------
    # ED50_totalVTE = ED50_majorVTE * ratio_ed50_totalvte * class ratio, so the
    # normalised dose for total VTE is r_<class> divided by those two ratios.
    # Enoxaparin is the reference class: its class ratio is 1 by definition.
    q_enox  <- r_enox  / (ratio_ed50_totalvte * 1)
    q_lmwh  <- r_lmwh  / (ratio_ed50_totalvte * ratio_ed50_totalvte_lmwh)
    q_dfxa  <- r_dfxa  / (ratio_ed50_totalvte * ratio_ed50_totalvte_dfxa)
    q_ifxa  <- r_ifxa  / (ratio_ed50_totalvte * ratio_ed50_totalvte_ifxa)
    q_uthr  <- r_uthr  / (ratio_ed50_totalvte * ratio_ed50_totalvte_uthr)
    q_bthr  <- r_bthr  / (ratio_ed50_totalvte * ratio_ed50_totalvte_bthr)
    q_mixed <- r_mixed / (ratio_ed50_totalvte * ratio_ed50_totalvte_mixed)
    q_hep   <- r_hep   / (ratio_ed50_totalvte * ratio_ed50_totalvte_heparin)
    q_warf  <- r_warf  / (ratio_ed50_totalvte * ratio_ed50_totalvte_warfarin)

    u_total <- q_enox^hill + q_lmwh^hill + q_dfxa^hill + q_ifxa^hill +
      q_uthr^hill + q_bthr^hill + q_mixed^hill + q_hep^hill + q_warf^hill

    g_totalvte <- (emax + emax_totalvte) * u_total / (u_total + 1)

    # ---- Bleeding ----------------------------------------------------------
    # Linear in Dose / (ED50_majorVTE * EDbld_class); major bleeding is the
    # reference endpoint (Sc = 0 by definition) and the other two endpoints
    # scale the common slope by (1 + Sc). Warfarin is deliberately absent:
    # the source could not estimate its EDbld.
    s_bld <- r_enox  / ratio_edbld_enox +
      r_lmwh  / ratio_edbld_lmwh +
      r_dfxa  / ratio_edbld_dfxa +
      r_ifxa  / ratio_edbld_ifxa +
      r_uthr  / ratio_edbld_uthr +
      r_bthr  / ratio_edbld_bthr +
      r_mixed / ratio_edbld_mixed +
      r_hep   / ratio_edbld_heparin

    g_majorbld <- s_bld
    g_crnmbld  <- s_bld * (1 + sc_bld_crnm)
    g_totalbld <- s_bld * (1 + sc_bld_total)

    # ---- Log odds ratios versus the trial's control arm --------------------
    # These are the source's stated primary response variable and are the
    # quantities identified without reference to any trial intercept.
    # Efficacy enters the logit with a minus sign (dose REDUCES VTE risk),
    # bleeding with a plus sign (dose INCREASES bleeding risk; Figure 2).
    lor_pe       <- -g_pe
    lor_majorvte <- -g_majorvte
    lor_totalvte <- -g_totalvte
    lor_majorbld <- g_majorbld
    lor_crnmbld  <- g_crnmbld
    lor_totalbld <- g_totalbld

    # ---- Per-arm event probabilities ---------------------------------------
    p_pe       <- 1 / (1 + exp(-(e0_pe       + lor_pe)))
    p_majorvte <- 1 / (1 + exp(-(e0_majorvte + lor_majorvte)))
    p_totalvte <- 1 / (1 + exp(-(e0_totalvte + lor_totalvte)))
    p_majorbld <- 1 / (1 + exp(-(e0_majorbld + lor_majorbld)))
    p_crnmbld  <- 1 / (1 + exp(-(e0_crnmbld  + lor_crnmbld)))
    p_totalbld <- 1 / (1 + exp(-(e0_totalbld + lor_totalbld)))

    # ---- Derived reporting quantities --------------------------------------
    # The x-axis of the source's Figures 1 and 2: the dose of enoxaparin that
    # would give the same major-VTE effect as this arm's dose of its own drug.
    dose_enox_equiv <- r_major * ed50_enoxaparin

    # 1 when the arm is a warfarin arm, whose bleeding predictions are NOT
    # supported: the source reports no EDbld for warfarin, so p_majorbld /
    # p_crnmbld / p_totalbld collapse to the control-arm rate for these arms.
    warfarin_bleeding_undefined <- (r_warf > 0)
  })
}
