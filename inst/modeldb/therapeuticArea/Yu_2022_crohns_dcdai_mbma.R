Yu_2022_crohns_dcdai_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of the change from baseline in the ",
    "Crohn's Disease Activity Index (Delta CDAI; negative values are ",
    "improvement) in adults with moderate-to-severe Crohn's disease, fitted ",
    "to study-arm summary data from 21 double-blind randomised controlled ",
    "trials (induction period) of biologics and small targeted molecules ",
    "published up to March 2020. The arm-mean change from baseline in CDAI ",
    "is the sum of a placebo term and one drug term per treatment, Y = E0 + ",
    "Edrug, on the natural scale of the endpoint. Each of 15 drugs carries ",
    "its own maximum effect; the dose-response drug is adalimumab (Emax ",
    "model in dose). Every drug effect rises as 1 - exp(-k * t) with one ",
    "shared onset rate. Arm-level baseline CDAI (SCORE_CDAI) acts as power ",
    "functions (covariate / dataset mean)^theta multiplying every drug ",
    "effect. The authors estimated the placebo response non-parametrically ",
    "(one value per trial per visit); this model holds it at the published ",
    "Week-12 typical-trial placebo value, which a user should replace with ",
    "the placebo response of the trial being simulated. Dose enters through ",
    "one CONMED_<drug>_DOSE covariate column per drug (0 = not received); ",
    "there is no PK layer and no rxode2 dose event. Simulation scope is ",
    "STUDY-ARM-MEAN outcomes, NOT individual patients. Companion models of ",
    "the same paper: Yu_2022_crohns_cdai150_mbma, ",
    "Yu_2022_crohns_cdai100_mbma, Yu_2022_crohns_cdai70_mbma, ",
    "Yu_2022_crohns_dcrp_mbma, Yu_2022_crohns_dibdq_mbma."
  )
  reference <- paste0(
    "Yu B, Zhao L, Jin S, He H, Zhang J, Wang X. Model-Based Meta-Analysis ",
    "on the Efficacy of Biologics and Small Targeted Molecules for Crohn's ",
    "Disease. Front Immunol. 2022;13:828219. doi:10.3389/fimmu.2022.828219. ",
    "PMC8967940. Structural model: Methods 'Model Development' Equations ",
    "1-3 and Supplementary Materials (DataSheet 2) Equations 1-12. ",
    "Parameter estimates: Supplementary Table 5 (DataSheet 2, 'Final ",
    "Models'); key parameters also in main-text Table 2. Model code: ",
    "DataSheet 2 'Model Code' (R nlme::gnls). Analysis dataset: ",
    "Supplementary DataSheet 1 (sheet 'CDAI')."
  )
  vignette <- "Yu_2022_crohns_mbma"

  # No PK layer, no concentration and no rxode2 dose events: dose is a
  # per-arm covariate column, as in He_2021_psoriasis_pasi75_mbma (same
  # research group and the same MBMA methodology).
  units <- list(
    time = "week (weeks since the first dose of the induction period; the onset rate constants are in 1/week and the paper simulates Week 12)",
    dosing = "mg (the dose metric of each CONMED_<drug>_DOSE column is described in its covariateData notes; weight-based regimens are expressed at 70 kg per Methods 'Data Development'. This model consumes NO rxode2 dose events.)",
    concentration = "points/arm (d_cdai is the STUDY-ARM mean change from baseline in CDAI; it is NOT a drug concentration)"
  )

  covariateData <- list(
    SCORE_CDAI = list(
      description = "Study-arm mean baseline Crohn's Disease Activity Index score.",
      units = "(score, 0-600)",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as (SCORE_CDAI / 306.285)^e_score_cdai_emax on every drug effect (DataSheet 2 Equation 12). 306.285 is the mean of the arm baseline CDAI over the 'CDAI' sheet of the deposited analysis dataset as hard-coded in the published model code. The paper's typical simulated trial has a baseline CDAI of 306.75 (Results 'Model Simulation'; Table 1 'Total' row).",
      source_name = "CDAI (DataSheet 1 column; DataSheet 2 model code)"
    ),
    N_ARM = list(
      description = "Number of patients contributing to the study-arm mean at a given visit.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity supplied per observation row; the meta-analytic weight of the residual. The source weights each arm by its standard error SD/sqrt(N) (DataSheet 2 Equation 8; model code varFixed(~SD^2/no)). It does not affect the typical-value prediction.",
      source_name = "no (DataSheet 1 column); N in DataSheet 2 Equations 8 and 11"
    ),
    CONMED_NATALIZUMAB_DOSE = list(
      description = "Natalizumab dose in the study arm; 0 if the arm did not receive natalizumab. Drug class: integrin-alpha4 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of natalizumab ('(FLAG == k)' in the model code), so only CONMED_NATALIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210 or 420 mg (3 or 6 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = natalizumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CDP571_DOSE = list(
      description = "CDP571 dose in the study arm; 0 if the arm did not receive CDP571. Drug class: TNF-alpha inhibitor (humanised anti-TNF antibody; development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of CDP571 ('(FLAG == k)' in the model code), so only CONMED_CDP571_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 700 mg (10 mg/kg at 70 kg); the IBDQ sheet records 10 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = CDP571 (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose in the study arm; 0 if the arm did not receive etanercept. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of etanercept ('(FLAG == k)' in the model code), so only CONMED_ETANERCEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 25 mg per subcutaneous administration (25 mg twice weekly). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = etanercept (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol dose in the study arm; 0 if the arm did not receive certolizumab pegol. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of certolizumab pegol ('(FLAG == k)' in the model code), so only CONMED_CERTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 100, 200 or 400 mg per subcutaneous administration (q2w or q4w regimens). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = certolizumab pegol (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab dose in the study arm; 0 if the arm did not receive adalimumab. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: Emax dose-response, effect = emax_adalimumab * DOSE / (DOSE + ED50) (model code), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: the week-0 induction dose, 40, 80 or 160 mg (regimens 40/20, 80/40 and 160/80 mg at weeks 0/2). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = adalimumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VEDOLIZUMAB_DOSE = list(
      description = "Vedolizumab dose in the study arm; 0 if the arm did not receive vedolizumab. Drug class: integrin-alpha4beta7 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vedolizumab ('(FLAG == k)' in the model code), so only CONMED_VEDOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 35 or 140 mg (0.5 or 2 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vedolizumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose in the study arm; 0 if the arm did not receive tofacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of tofacitinib ('(FLAG == k)' in the model code), so only CONMED_TOFACITINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 1, 5, 10 or 15 mg per administration of a twice-daily regimen. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = tofacitinib (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FILGOTINIB_DOSE = list(
      description = "Filgotinib dose in the study arm; 0 if the arm did not receive filgotinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of filgotinib ('(FLAG == k)' in the model code), so only CONMED_FILGOTINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 200 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = filgotinib (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ABATACEPT_DOSE = list(
      description = "Abatacept dose in the study arm; 0 if the arm did not receive abatacept. Drug class: T-cell activation inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of abatacept ('(FLAG == k)' in the model code), so only CONMED_ABATACEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210, 700 or 2100 mg (3, 10 or 30 mg/kg at 70 kg); the IBDQ sheet records 3, 10 or 30 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = abatacept (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_RISANKIZUMAB_DOSE = list(
      description = "Risankizumab dose in the study arm; 0 if the arm did not receive risankizumab. Drug class: IL-23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of risankizumab ('(FLAG == k)' in the model code), so only CONMED_RISANKIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 200 or 600 mg per intravenous infusion (q4w). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = risankizumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_PF04236921_DOSE = list(
      description = "PF-04236921 dose in the study arm; 0 if the arm did not receive PF-04236921. Drug class: IL-6 inhibitor (development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of PF-04236921 ('(FLAG == k)' in the model code), so only CONMED_PF04236921_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 10, 50 or 200 mg per subcutaneous administration (weeks 0 and 4). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = PF-04236921 (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_BRAZIKUMAB_DOSE = list(
      description = "Brazikumab dose in the study arm; 0 if the arm did not receive brazikumab. Drug class: IL-23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of brazikumab ('(FLAG == k)' in the model code), so only CONMED_BRAZIKUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 700 mg per intravenous infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = brazikumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FONTOLIZUMAB_DOSE = list(
      description = "Fontolizumab dose in the study arm; 0 if the arm did not receive fontolizumab. Drug class: IFN-gamma inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of fontolizumab ('(FLAG == k)' in the model code), so only CONMED_FONTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 70, 280 or 700 mg (1, 4 or 10 mg/kg loading dose at 70 kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = fontolizumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ABRILUMAB_DOSE = list(
      description = "Abrilumab dose in the study arm; 0 if the arm did not receive abrilumab. Drug class: integrin-alpha4beta7 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of abrilumab ('(FLAG == k)' in the model code), so only CONMED_ABRILUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 21, 70 or 210 mg per subcutaneous administration. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = abrilumab (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VERCIRNON_DOSE = list(
      description = "Vercirnon dose in the study arm; 0 if the arm did not receive vercirnon. Drug class: CCR9 antagonist.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vercirnon ('(FLAG == k)' in the model code), so only CONMED_VERCIRNON_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: the total daily dose, 250, 500 or 1000 mg (250 mg bid is recorded as 500). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vercirnon (DataSheet 1 'CDAI' sheet); Table 1 'Route (regimen)' column"
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
    T_DIAG_CD = list(
      description = "Study-arm mean duration of Crohn's disease.",
      units = "years",
      type = "continuous",
      notes = "Screened but retained only in the Delta CRP and Delta IBDQ models; see the Yu_2022_crohns_mbma vignette."
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
    notes = "MBMA at the STUDY-ARM level: 46 trials, 146 treatment arms and 12,846 patients in total (Results 'Available Data'), of which 21 trials report change from baseline in CDAI (CDAI sheet of the deposited dataset). Female percentage is 100 minus the Table 1 'Total' percentage of men (46.09). The model has no between-study random effect: the source absorbs between-trial differences into its non-parametric placebo, one estimate per trial and visit."
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
    e0_d_cdai <- fixed(-58.67)
    label("Placebo change from baseline at the typical simulated trial, Week 12 (points)")
    # Results 'Model Simulation': 'The placebo effect was estimated as
    # -58.67'.
    # The fitted model estimates one placebo value per trial and visit (gnls
    # params 'eo ~ -1 + I(group.week)') and the paper reports none of them.
    # The Week-12 simulation value is held constant in time here; replace it
    # with the observed placebo response of the trial being simulated.

    # ---- Onset of the drug effect -------------------------------------------
    lkdrug <- log(0.22)
    label("Log onset rate constant of the shared drug effect (log 1/week)")
    # Supplementary Table 5 row 'kgeneral, Rate constant for the onset of all
    # drugs' = 0.22 (95% CI 0.14, 0.34); also main-text Table 2. Results:
    # 'ET50 and ET90 were assumed to be 3.2 and 10.5 weeks in the Delta CDAI
    # model' = log(2)/0.22 and log(10)/0.22. The code writes the rate as
    # exp(k).

    # ---- Covariate effects (power form, Equation 12) -----------------------
    e_score_cdai_emax <- -2.05
    label("Power exponent on baseline CDAI over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 5 and main-text Table 2 row 'Covariate: Baseline
    # CDAI' = -2.05 (95% CI -4.93, 0.83).

    # ---- Drug effects (Supplementary Table 5, 'Emax' block) ----
    emax_natalizumab <- -51.36
    label("Effect of natalizumab at any dose (points)")
    # Supplementary Table 5, 'Natalizumab' = -51.36 (95% CI -70.89, -31.83).
    emax_cdp571 <- -18.26
    label("Effect of CDP571 at any dose (points)")
    # Supplementary Table 5, 'CDP571' = -18.26 (95% CI -33.68, -2.85).
    emax_etanercept <- -15.93
    label("Effect of etanercept at any dose (points)")
    # Supplementary Table 5, 'Etanercept' = -15.93 (95% CI -84.16, 52.30).
    emax_certolizumab <- -31.53
    label("Effect of certolizumab pegol at any dose (points)")
    # Supplementary Table 5, 'Certolizumab pegol' = -31.53 (95% CI -45.48,
    # -17.59).
    emax_adalimumab <- -151.24
    label("Maximum effect of adalimumab in the Emax dose-response (points)")
    # Supplementary Table 5 and main-text Table 2, 'Adalimumab (Emax)' =
    # -151.24 (95% CI -322.11, 19.64).
    led50_adalimumab <- log(112.23)
    label("Log dose giving half the maximum effect of adalimumab (log mg)")
    # Supplementary Table 5 and main-text Table 2, 'Adalimumab (ED50)' =
    # 112.23 (95% CI 9.02, 1.40x10^3); the interval is symmetric on the log
    # scale, matching the code's exp(led50) parameterisation.
    emax_vedolizumab <- -22.38
    label("Effect of vedolizumab at any dose (points)")
    # Supplementary Table 5, 'Vedolizumab' = -22.38 (95% CI -49.72, 4.97).
    emax_tofacitinib <- -26.23
    label("Effect of tofacitinib at any dose (points)")
    # Supplementary Table 5, 'Tofacitinib' = -26.23 (95% CI -49.46, -2.99).
    emax_filgotinib <- -32.16
    label("Effect of filgotinib at any dose (points)")
    # Supplementary Table 5, 'Filgotinib' = -32.16 (95% CI -69.00, 4.68).
    emax_abatacept <- -20.19
    label("Effect of abatacept at any dose (points)")
    # Supplementary Table 5, 'Abatacept' = -20.19 (95% CI -43.15, 2.76).
    emax_risankizumab <- -81.85
    label("Effect of risankizumab at any dose (points)")
    # Supplementary Table 5, 'Risankizumab' = -81.85 (95% CI -122.85, -40.86).
    emax_pf04236921 <- -31.11
    label("Effect of PF-04236921 at any dose (points)")
    # Supplementary Table 5, 'PF-04236921' = -31.11 (95% CI -60.23, -1.99).
    emax_brazikumab <- -39.99
    label("Effect of brazikumab at any dose (points)")
    # Supplementary Table 5, 'Brazikumab' = -39.99 (95% CI -85.67, 5.70).
    emax_fontolizumab <- -41
    label("Effect of fontolizumab at any dose (points)")
    # Supplementary Table 5, 'Fontolizumab' = -41.00 (95% CI -66.47, -15.52).
    emax_abrilumab <- -35.63
    label("Effect of abrilumab at any dose (points)")
    # Supplementary Table 5, 'Abrilumab' = -35.63 (95% CI -41.90, -29.36).
    emax_vercirnon <- -17.62
    label("Effect of vercirnon at any dose (points)")
    # Supplementary Table 5, 'Vercirnon' = -17.62 (95% CI -50.29, 15.06).

    # ---- Residual error ------------------------------------------------------
    addSd <- fixed(79.4)
    label("Between-patient SD of the change; arm-mean residual SD is addSd/sqrt(N_ARM) (points)")
    # NOT REPORTED in the paper. The source weights each arm by its own
    # reported SD/sqrt(N) (varFixed(~SD^2/no)) with an estimated scale
    # factor. Encoded at unit scale with the median arm SD of the 'CDAI'
    # sheet of the deposited dataset (79.4); supply an arm's own SD by
    # overriding addSd. The source also fits a within-arm residual
    # autocorrelation (AR1), which a per-record rxode2 residual cannot
    # carry.
  })

  model({
    # Equation 12: every drug effect is scaled by the arm's covariates.
    covEff <- (SCORE_CDAI / 306.285)^e_score_cdai_emax

    # Equation 3, one onset rate shared by every drug.
    onsetDrug <- 1 - exp(-exp(lkdrug) * time)

    # One term per drug. An arm supplies a positive dose in exactly one
    # CONMED_<drug>_DOSE column and zero in the rest, so a placebo arm (all
    # columns zero) reduces to the placebo term alone.
    edNatalizumab <- emax_natalizumab * (CONMED_NATALIZUMAB_DOSE > 0) * onsetDrug
    edCdp571 <- emax_cdp571 * (CONMED_CDP571_DOSE > 0) * onsetDrug
    edEtanercept <- emax_etanercept * (CONMED_ETANERCEPT_DOSE > 0) * onsetDrug
    edCertolizumab <- emax_certolizumab * (CONMED_CERTOLIZUMAB_DOSE > 0) * onsetDrug
    edAdalimumab <- emax_adalimumab *
      CONMED_ADALIMUMAB_DOSE / (CONMED_ADALIMUMAB_DOSE + exp(led50_adalimumab)) * onsetDrug
    edVedolizumab <- emax_vedolizumab * (CONMED_VEDOLIZUMAB_DOSE > 0) * onsetDrug
    edTofacitinib <- emax_tofacitinib * (CONMED_TOFACITINIB_DOSE > 0) * onsetDrug
    edFilgotinib <- emax_filgotinib * (CONMED_FILGOTINIB_DOSE > 0) * onsetDrug
    edAbatacept <- emax_abatacept * (CONMED_ABATACEPT_DOSE > 0) * onsetDrug
    edRisankizumab <- emax_risankizumab * (CONMED_RISANKIZUMAB_DOSE > 0) * onsetDrug
    edPf04236921 <- emax_pf04236921 * (CONMED_PF04236921_DOSE > 0) * onsetDrug
    edBrazikumab <- emax_brazikumab * (CONMED_BRAZIKUMAB_DOSE > 0) * onsetDrug
    edFontolizumab <- emax_fontolizumab * (CONMED_FONTOLIZUMAB_DOSE > 0) * onsetDrug
    edAbrilumab <- emax_abrilumab * (CONMED_ABRILUMAB_DOSE > 0) * onsetDrug
    edVercirnon <- emax_vercirnon * (CONMED_VERCIRNON_DOSE > 0) * onsetDrug

    edrug <- (edNatalizumab + edCdp571 + edEtanercept + edCertolizumab +
      edAdalimumab + edVedolizumab + edTofacitinib + edFilgotinib +
      edAbatacept + edRisankizumab + edPf04236921 + edBrazikumab +
      edFontolizumab + edAbrilumab + edVercirnon) * covEff

    # Equation 7: placebo plus drug term on the natural scale.
    d_cdai <- e0_d_cdai + edrug

    # Equation 8: the residual SD of an arm mean scales with 1/sqrt(N).
    sdArm <- addSd / sqrt(N_ARM)
    d_cdai ~ add(sdArm)
  })
}
