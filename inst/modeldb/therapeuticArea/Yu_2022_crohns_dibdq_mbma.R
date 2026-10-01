Yu_2022_crohns_dibdq_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of the change from baseline in the ",
    "Inflammatory Bowel Disease Questionnaire score (Delta IBDQ; positive ",
    "values are improvement) in adults with moderate-to-severe Crohn's ",
    "disease, fitted to study-arm summary data from 20 double-blind ",
    "randomised controlled trials (induction period) of biologics and small ",
    "targeted molecules published up to March 2020. The arm-mean change ",
    "from baseline in IBDQ is the sum of a placebo term and one drug term ",
    "per treatment, Y = E0 + Edrug, on the natural scale of the endpoint. ",
    "Each of 16 drugs carries its own maximum effect; the dose-response ",
    "drugs are adalimumab (linear slope per mg) and risankizumab (linear ",
    "slope per mg). Drug effects are constant in time. Arm-level disease ",
    "duration (T_DIAG_CD) acts as power functions (covariate / dataset ",
    "mean)^theta multiplying every drug effect. The authors estimated the ",
    "placebo response non-parametrically (one value per trial per visit); ",
    "this model holds it at the published Week-12 typical-trial placebo ",
    "value, which a user should replace with the placebo response of the ",
    "trial being simulated. Dose enters through one CONMED_<drug>_DOSE ",
    "covariate column per drug (0 = not received); there is no PK layer and ",
    "no rxode2 dose event. Simulation scope is STUDY-ARM-MEAN outcomes, NOT ",
    "individual patients. Companion models of the same paper: ",
    "Yu_2022_crohns_cdai150_mbma, Yu_2022_crohns_cdai100_mbma, ",
    "Yu_2022_crohns_cdai70_mbma, Yu_2022_crohns_dcdai_mbma, ",
    "Yu_2022_crohns_dcrp_mbma."
  )
  reference <- paste0(
    "Yu B, Zhao L, Jin S, He H, Zhang J, Wang X. Model-Based Meta-Analysis ",
    "on the Efficacy of Biologics and Small Targeted Molecules for Crohn's ",
    "Disease. Front Immunol. 2022;13:828219. doi:10.3389/fimmu.2022.828219. ",
    "PMC8967940. Structural model: Methods 'Model Development' Equations ",
    "1-3 and Supplementary Materials (DataSheet 2) Equations 1-12. ",
    "Parameter estimates: Supplementary Table 7 (DataSheet 2, 'Final ",
    "Models'); key parameters also in main-text Table 2. Model code: ",
    "DataSheet 2 'Model Code' (R nlme::gnls). Analysis dataset: ",
    "Supplementary DataSheet 1 (sheet 'IBDQ')."
  )
  vignette <- "Yu_2022_crohns_mbma"

  # No PK layer, no concentration and no rxode2 dose events: dose is a
  # per-arm covariate column, as in He_2021_psoriasis_pasi75_mbma (same
  # research group and the same MBMA methodology).
  units <- list(
    time = "week (weeks since the first dose of the induction period; the onset rate constants are in 1/week and the paper simulates Week 12)",
    dosing = "mg (the dose metric of each CONMED_<drug>_DOSE column is described in its covariateData notes; weight-based regimens are expressed at 70 kg per Methods 'Data Development'. This model consumes NO rxode2 dose events.)",
    concentration = "points/arm (d_ibdq is the STUDY-ARM mean change from baseline in IBDQ; it is NOT a drug concentration)"
  )

  covariateData <- list(
    T_DIAG_CD = list(
      description = "Study-arm mean duration of Crohn's disease at baseline (time since diagnosis).",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as (T_DIAG_CD / 9.679)^e_t_diag_cd_emax on every drug effect; 9.679 is the mean arm disease duration over the 'IBDQ' sheet of the deposited dataset with missing arm values replaced by the median arm value. The exponent is steep, so this reference value matters: the alternative imputation readings give a mean within about 1% of it, which moves the drug effects by up to about 20%. The paper's typical simulated trial has a disease duration of 9.52 years.",
      source_name = "dur (DataSheet 1 column; DataSheet 2 model code); 'disease duration' in Table 1 and Table 2"
    ),
    N_ARM = list(
      description = "Number of patients contributing to the study-arm mean at a given visit.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity supplied per observation row; the meta-analytic weight of the residual. The source weights each arm by its standard error SD/sqrt(N) (DataSheet 2 Equation 8; model code varFixed(~SD^2/no)). It does not affect the typical-value prediction.",
      source_name = "no (DataSheet 1 column); N in DataSheet 2 Equations 8 and 11"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description = "Infliximab dose in the study arm; 0 if the arm did not receive infliximab. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of infliximab ('(FLAG == k)' in the model code), so only CONMED_INFLIXIMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 350 mg (5 mg/kg at the 70 kg per-patient normalisation of Methods 'Data Development'); the IBDQ sheet of the deposited dataset records the same arms as 5 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = infliximab (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_NATALIZUMAB_DOSE = list(
      description = "Natalizumab dose in the study arm; 0 if the arm did not receive natalizumab. Drug class: integrin-alpha4 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of natalizumab ('(FLAG == k)' in the model code), so only CONMED_NATALIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210 or 420 mg (3 or 6 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = natalizumab (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CDP571_DOSE = list(
      description = "CDP571 dose in the study arm; 0 if the arm did not receive CDP571. Drug class: TNF-alpha inhibitor (humanised anti-TNF antibody; development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of CDP571 ('(FLAG == k)' in the model code), so only CONMED_CDP571_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 700 mg (10 mg/kg at 70 kg); the IBDQ sheet records 10 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = CDP571 (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose in the study arm; 0 if the arm did not receive etanercept. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of etanercept ('(FLAG == k)' in the model code), so only CONMED_ETANERCEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 25 mg per subcutaneous administration (25 mg twice weekly). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = etanercept (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol dose in the study arm; 0 if the arm did not receive certolizumab pegol. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of certolizumab pegol ('(FLAG == k)' in the model code), so only CONMED_CERTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 100, 200 or 400 mg per subcutaneous administration (q2w or q4w regimens). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = certolizumab pegol (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab dose in the study arm; 0 if the arm did not receive adalimumab. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_adalimumab * CONMED_ADALIMUMAB_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: the week-0 induction dose, 40, 80 or 160 mg (regimens 40/20, 80/40 and 160/80 mg at weeks 0/2). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = adalimumab (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VEDOLIZUMAB_DOSE = list(
      description = "Vedolizumab dose in the study arm; 0 if the arm did not receive vedolizumab. Drug class: integrin-alpha4beta7 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vedolizumab ('(FLAG == k)' in the model code), so only CONMED_VEDOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 35 or 140 mg (0.5 or 2 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vedolizumab (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_SEMAPIMOD_DOSE = list(
      description = "Semapimod dose in the study arm; 0 if the arm did not receive semapimod. Drug class: TNF-alpha inhibitor group of Table 1.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of semapimod ('(FLAG == k)' in the model code), so only CONMED_SEMAPIMOD_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 180 mg (60 mg i.v. daily for 3 days). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = semapimod (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_APILIMOD_DOSE = list(
      description = "Apilimod dose in the study arm; 0 if the arm did not receive apilimod. Drug class: IL-12/23 inhibitor (oral small molecule).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of apilimod ('(FLAG == k)' in the model code), so only CONMED_APILIMOD_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 50 or 100 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = apilimod (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FILGOTINIB_DOSE = list(
      description = "Filgotinib dose in the study arm; 0 if the arm did not receive filgotinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of filgotinib ('(FLAG == k)' in the model code), so only CONMED_FILGOTINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 200 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = filgotinib (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose in the study arm; 0 if the arm did not receive tofacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of tofacitinib ('(FLAG == k)' in the model code), so only CONMED_TOFACITINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 1, 5, 10 or 15 mg per administration of a twice-daily regimen. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = tofacitinib (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ABATACEPT_DOSE = list(
      description = "Abatacept dose in the study arm; 0 if the arm did not receive abatacept. Drug class: T-cell activation inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of abatacept ('(FLAG == k)' in the model code), so only CONMED_ABATACEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210, 700 or 2100 mg (3, 10 or 30 mg/kg at 70 kg); the IBDQ sheet records 3, 10 or 30 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = abatacept (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_RISANKIZUMAB_DOSE = list(
      description = "Risankizumab dose in the study arm; 0 if the arm did not receive risankizumab. Drug class: IL-23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_risankizumab * CONMED_RISANKIZUMAB_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: 200 or 600 mg per intravenous infusion (q4w). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = risankizumab (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_UPADACITINIB_DOSE = list(
      description = "Upadacitinib dose in the study arm; 0 if the arm did not receive upadacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of upadacitinib ('(FLAG == k)' in the model code), so only CONMED_UPADACITINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 3, 6, 12 or 24 mg per administration of a twice-daily regimen; the 24 mg once-daily arm is recorded as 12, its twice-daily equivalent. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = upadacitinib (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FONTOLIZUMAB_DOSE = list(
      description = "Fontolizumab dose in the study arm; 0 if the arm did not receive fontolizumab. Drug class: IFN-gamma inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of fontolizumab ('(FLAG == k)' in the model code), so only CONMED_FONTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 70, 280 or 700 mg (1, 4 or 10 mg/kg loading dose at 70 kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = fontolizumab (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VERCIRNON_DOSE = list(
      description = "Vercirnon dose in the study arm; 0 if the arm did not receive vercirnon. Drug class: CCR9 antagonist.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vercirnon ('(FLAG == k)' in the model code), so only CONMED_VERCIRNON_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: the total daily dose, 250, 500 or 1000 mg (250 mg bid is recorded as 500). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vercirnon (DataSheet 1 'IBDQ' sheet); Table 1 'Route (regimen)' column"
    )
  )

  # Screened in the source (Methods 'Covariate': age, percentage of male,
  # disease duration, smoking status, CDAI, CRP, IBDQ) but not retained in
  # this endpoint's final model. Baseline IBDQ was dropped from every model
  # because more than 40% of trials lacked it (Methods 'Data Development').
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female fraction of the study arm (the source screens the percentage of men).",
      units = "fraction",
      type = "continuous",
      notes = "Screened as 'percentage of male' (Methods 'Covariate'; Results 'Covariates') but not retained in this model; see the Yu_2022_crohns_mbma vignette."
    ),
    SMOKE = list(
      description = "Fraction of current smokers in the study arm.",
      units = "fraction",
      type = "continuous",
      notes = "Screened as 'smoking status' but not retained in any of the six final models; see the Yu_2022_crohns_mbma vignette."
    ),
    AGE = list(
      description = "Study-arm mean age.",
      units = "years",
      type = "continuous",
      notes = "Screened but retained only in the Delta CRP model (Yu_2022_crohns_dcrp_mbma); see the Yu_2022_crohns_mbma vignette."
    ),
    SCORE_CDAI = list(
      description = "Study-arm mean baseline CDAI.",
      units = "(score, 0-600)",
      type = "continuous",
      notes = "Screened but not retained in the Delta CRP and Delta IBDQ models; see the Yu_2022_crohns_mbma vignette."
    ),
    CRP = list(
      description = "Study-arm mean baseline CRP.",
      units = "mg/dL (see Yu_2022_crohns_dcrp_mbma)",
      type = "continuous",
      notes = "Screened but not retained in the Delta CDAI and Delta IBDQ models; see the Yu_2022_crohns_mbma vignette."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12846L,
    n_studies = 46L,
    age_range = "adults (>= 18 years); per-drug arm-mean ages 34.8-41.0 years (Table 1; overall mean 37.2)",
    weight_range = "not reported; weight-based doses are normalised to 70 kg",
    sex_female_pct = 53.91,
    race_ethnicity = "not reported",
    disease_state = "moderate-to-severe active Crohn's disease confirmed by radiologic, endoscopic or histologic criteria; induction-period data only; prior TNF-inhibitor exposure and concomitant 5-ASA, oral steroids and immunomodulators allowed. Table 1 overall means: disease duration 9.52 years, baseline CDAI 306.75, baseline CRP 1.67, baseline IBDQ 127.37",
    dose_range = "per Table 1 'Route (regimen)'; see the covariateData notes for the dose metric of each drug column",
    regions = "international; MEDLINE, CENTRAL, EMBASE and ClinicalTrials.gov to 14 March 2020 plus UEGW, ACG, DDW and ECCO abstracts to 2019",
    notes = "MBMA at the STUDY-ARM level: 46 trials, 146 treatment arms and 12,846 patients in total (Results 'Available Data'), of which 20 trials report change from baseline in IBDQ (IBDQ sheet of the deposited dataset). Female percentage is 100 minus the Table 1 'Total' percentage of men (46.09). The model has no between-study random effect: the source absorbs between-trial differences into its non-parametric placebo, one estimate per trial and visit."
  )

  ini({
    # ========================================================================
    # Structure (Methods 'Model Development'; DataSheet 2 Equations 1-12 and
    # 'Code for ... model'). The published gnls() code is the authority:
    #
    #   Y_ijt   = E0_it + Edrug_ijt                        (Equation 7)
    #   Edrug   = Emax_drug * f(time) * prod((X / mean(X))^theta)   (Eq. 2, 3, 12)
    #   Emax_drug = constant, slope * DOSE, or Emax * DOSE / (DOSE + ED50)
    #   (Supplementary table footnotes a and b; the footnote typesets the
    #   Emax denominator as 'ED50 * dose', the code as DOSE + exp(led50)).
    #
    # With the values below every study-arm prediction of the paper's own
    # Week-12 ranking figure is reproduced; see the vignette.
    # ========================================================================

    # ---- Placebo (the source estimates it non-parametrically) --------------
    e0_d_ibdq <- fixed(17.36)
    label("Placebo change from baseline at the typical simulated trial, Week 12 (points)")
    # Results 'Model Simulation': 'The placebo effect was simulated as 17.36
    # with a longitudinal model'.
    # The fitted model estimates one placebo value per trial and visit (gnls
    # params 'eo ~ -1 + I(group.week)') and the paper reports none of them.
    # The Week-12 simulation value is held constant in time here; replace it
    # with the observed placebo response of the trial being simulated.

    # ---- Covariate effects (power form, Equation 12) -----------------------
    e_t_diag_cd_emax <- -8.98
    label("Power exponent on disease duration over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 7 and main-text Table 2 row 'Covariate: Disease
    # duration' = -8.98 (95% CI -10.60, -7.36).

    # ---- Drug effects (Supplementary Table 7, 'Emax' block) ----
    emax_infliximab <- 15.72
    label("Effect of infliximab at any dose (points)")
    # Supplementary Table 7, 'Infliximab' = 15.72 (95% CI 7.59, 23.85).
    emax_natalizumab <- 2.72
    label("Effect of natalizumab at any dose (points)")
    # Supplementary Table 7, 'Natalizumab' = 2.72 (95% CI 0.66, 4.79).
    emax_cdp571 <- 3.93
    label("Effect of CDP571 at any dose (points)")
    # Supplementary Table 7, 'CDP571' = 3.93 (95% CI 0.55, 7.31).
    emax_etanercept <- -1.51
    label("Effect of etanercept at any dose (points)")
    # Supplementary Table 7, 'Etanercept' = -1.51 (95% CI -26.40, 23.38).
    emax_certolizumab <- 0.31
    label("Effect of certolizumab pegol at any dose (points)")
    # Supplementary Table 7, 'Certolizumab pegol' = 0.31 (95% CI -0.19, 0.82).
    slope_adalimumab <- 0.1
    label("Linear dose-response slope of adalimumab (points per mg)")
    # Supplementary Table 7 and main-text Table 2, 'Adalimumab (slope)' = 0.10
    # (95% CI 0.05, 0.15).
    emax_vedolizumab <- 0.97
    label("Effect of vedolizumab at any dose (points)")
    # Supplementary Table 7, 'Vedolizumab' = 0.97 (95% CI -1.80, 3.75).
    emax_semapimod <- -1.3
    label("Effect of semapimod at any dose (points)")
    # Supplementary Table 7, 'Semapimod' = -1.30 (95% CI -27.19, 24.60).
    emax_apilimod <- -8.26
    label("Effect of apilimod at any dose (points)")
    # Supplementary Table 7, 'Apilimod' = -8.26 (95% CI -41.19, 24.68).
    emax_filgotinib <- 6.74
    label("Effect of filgotinib at any dose (points)")
    # Supplementary Table 7, 'Filgotinib' = 6.74 (95% CI -0.05, 13.54).
    emax_tofacitinib <- 52.74
    label("Effect of tofacitinib at any dose (points)")
    # Supplementary Table 7, 'Tofacitinib' = 52.74 (95% CI 39.23, 66.25).
    emax_abatacept <- -3.54
    label("Effect of abatacept at any dose (points)")
    # Supplementary Table 7, 'Abatacept' = -3.54 (95% CI -4.35, -2.72).
    slope_risankizumab <- 0.05
    label("Linear dose-response slope of risankizumab (points per mg)")
    # Supplementary Table 7 and main-text Table 2, 'Risankizumab (slope)' =
    # 0.05 (95% CI 0.01, 0.10).
    emax_upadacitinib <- 4.7
    label("Effect of upadacitinib at any dose (points)")
    # Supplementary Table 7, 'Upadacitinib' = 4.70 (95% CI -1.57, 10.97).
    emax_fontolizumab <- -0.34
    label("Effect of fontolizumab at any dose (points)")
    # Supplementary Table 7, 'Fontolizumab' = -0.34 (95% CI -4.02, 3.35).
    emax_vercirnon <- 3.54
    label("Effect of vercirnon at any dose (points)")
    # Supplementary Table 7, 'Vercirnon' = 3.54 (95% CI -5.84, 12.93).

    # ---- Residual error ------------------------------------------------------
    addSd <- fixed(31.1)
    label("Between-patient SD of the change; arm-mean residual SD is addSd/sqrt(N_ARM) (points)")
    # NOT REPORTED in the paper. The source weights each arm by its own
    # reported SD/sqrt(N) (varFixed(~SD^2/no)) with an estimated scale
    # factor. Encoded at unit scale with the median arm SD of the 'IBDQ'
    # sheet of the deposited dataset (31.1); supply an arm's own SD by
    # overriding addSd. The source also fits a within-arm residual
    # autocorrelation (AR1), which a per-record rxode2 residual cannot
    # carry.
  })

  model({
    # Equation 12: every drug effect is scaled by the arm's covariates.
    covEff <- (T_DIAG_CD / 9.679)^e_t_diag_cd_emax

    # One term per drug. An arm supplies a positive dose in exactly one
    # CONMED_<drug>_DOSE column and zero in the rest, so a placebo arm (all
    # columns zero) reduces to the placebo term alone.
    edInfliximab <- emax_infliximab * (CONMED_INFLIXIMAB_DOSE > 0)
    edNatalizumab <- emax_natalizumab * (CONMED_NATALIZUMAB_DOSE > 0)
    edCdp571 <- emax_cdp571 * (CONMED_CDP571_DOSE > 0)
    edEtanercept <- emax_etanercept * (CONMED_ETANERCEPT_DOSE > 0)
    edCertolizumab <- emax_certolizumab * (CONMED_CERTOLIZUMAB_DOSE > 0)
    edAdalimumab <- slope_adalimumab * CONMED_ADALIMUMAB_DOSE
    edVedolizumab <- emax_vedolizumab * (CONMED_VEDOLIZUMAB_DOSE > 0)
    edSemapimod <- emax_semapimod * (CONMED_SEMAPIMOD_DOSE > 0)
    edApilimod <- emax_apilimod * (CONMED_APILIMOD_DOSE > 0)
    edFilgotinib <- emax_filgotinib * (CONMED_FILGOTINIB_DOSE > 0)
    edTofacitinib <- emax_tofacitinib * (CONMED_TOFACITINIB_DOSE > 0)
    edAbatacept <- emax_abatacept * (CONMED_ABATACEPT_DOSE > 0)
    edRisankizumab <- slope_risankizumab * CONMED_RISANKIZUMAB_DOSE
    edUpadacitinib <- emax_upadacitinib * (CONMED_UPADACITINIB_DOSE > 0)
    edFontolizumab <- emax_fontolizumab * (CONMED_FONTOLIZUMAB_DOSE > 0)
    edVercirnon <- emax_vercirnon * (CONMED_VERCIRNON_DOSE > 0)

    edrug <- (edInfliximab + edNatalizumab + edCdp571 + edEtanercept +
      edCertolizumab + edAdalimumab + edVedolizumab + edSemapimod +
      edApilimod + edFilgotinib + edTofacitinib + edAbatacept +
      edRisankizumab + edUpadacitinib + edFontolizumab + edVercirnon) * covEff

    # Equation 7: placebo plus drug term on the natural scale.
    d_ibdq <- e0_d_ibdq + edrug

    # Equation 8: the residual SD of an arm mean scales with 1/sqrt(N).
    sdArm <- addSd / sqrt(N_ARM)
    d_ibdq ~ add(sdArm)
  })
}
