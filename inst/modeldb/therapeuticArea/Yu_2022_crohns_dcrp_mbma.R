Yu_2022_crohns_dcrp_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of the change from baseline in ",
    "C-reactive protein (Delta CRP; negative values are improvement) in ",
    "adults with moderate-to-severe Crohn's disease, fitted to study-arm ",
    "summary data from 26 double-blind randomised controlled trials ",
    "(induction period) of biologics and small targeted molecules published ",
    "up to March 2020. The arm-mean change from baseline in CRP is the sum ",
    "of a placebo term and one drug term per treatment, Y = E0 + Edrug, on ",
    "the natural scale of the endpoint. Each of 17 drugs carries its own ",
    "maximum effect; the dose-response drugs are certolizumab pegol (linear ",
    "slope per mg), PF-04236921 (Emax model in dose) and upadacitinib ",
    "(linear slope per mg). Drug effects are constant in time. Arm-level ",
    "age (AGE), disease duration (T_DIAG_CD) and baseline CRP (CRP) act as ",
    "power functions (covariate / dataset mean)^theta multiplying every ",
    "drug effect. The authors estimated the placebo response ",
    "non-parametrically (one value per trial per visit); this model holds ",
    "it at the published Week-12 typical-trial placebo value, which a user ",
    "should replace with the placebo response of the trial being simulated. ",
    "Dose enters through one CONMED_<drug>_DOSE covariate column per drug ",
    "(0 = not received); there is no PK layer and no rxode2 dose event. ",
    "Simulation scope is STUDY-ARM-MEAN outcomes, NOT individual patients. ",
    "Companion models of the same paper: Yu_2022_crohns_cdai150_mbma, ",
    "Yu_2022_crohns_cdai100_mbma, Yu_2022_crohns_cdai70_mbma, ",
    "Yu_2022_crohns_dcdai_mbma, Yu_2022_crohns_dibdq_mbma."
  )
  reference <- paste0(
    "Yu B, Zhao L, Jin S, He H, Zhang J, Wang X. Model-Based Meta-Analysis ",
    "on the Efficacy of Biologics and Small Targeted Molecules for Crohn's ",
    "Disease. Front Immunol. 2022;13:828219. doi:10.3389/fimmu.2022.828219. ",
    "PMC8967940. Structural model: Methods 'Model Development' Equations ",
    "1-3 and Supplementary Materials (DataSheet 2) Equations 1-12. ",
    "Parameter estimates: Supplementary Table 6 (DataSheet 2, 'Final ",
    "Models'); key parameters also in main-text Table 2. Model code: ",
    "DataSheet 2 'Model Code' (R nlme::gnls). Analysis dataset: ",
    "Supplementary DataSheet 1 (sheet 'CRP')."
  )
  vignette <- "Yu_2022_crohns_mbma"

  # No PK layer, no concentration and no rxode2 dose events: dose is a
  # per-arm covariate column, as in He_2021_psoriasis_pasi75_mbma (same
  # research group and the same MBMA methodology).
  units <- list(
    time = "week (weeks since the first dose of the induction period; the onset rate constants are in 1/week and the paper simulates Week 12)",
    dosing = "mg (the dose metric of each CONMED_<drug>_DOSE column is described in its covariateData notes; weight-based regimens are expressed at 70 kg per Methods 'Data Development'. This model consumes NO rxode2 dose events.)",
    concentration = "mg/dL (see notes)/arm (d_crp is the STUDY-ARM mean change from baseline in CRP; it is NOT a drug concentration)"
  )

  covariateData <- list(
    AGE = list(
      description = "Study-arm mean age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as (AGE / 36.945)^e_age_emax on every drug effect; 36.945 is the mean arm age over the 'CRP' sheet of the deposited dataset (no missing values). The paper does not state the age of its typical simulated trial; the Table 1 'Total' mean is 37.20 years.",
      source_name = "age (DataSheet 1 column; DataSheet 2 model code)"
    ),
    T_DIAG_CD = list(
      description = "Study-arm mean duration of Crohn's disease at baseline (time since diagnosis).",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as (T_DIAG_CD / 9.317)^e_t_diag_cd_emax on every drug effect; 9.317 is the mean arm disease duration over the 'CRP' sheet of the deposited dataset with missing arm values replaced by the median arm value. The exponent is steep, so this reference value matters: the alternative imputation readings give a mean within about 1% of it, which moves the drug effects by up to about 5%. The paper's typical simulated trial has a disease duration of 9.52 years.",
      source_name = "dur (DataSheet 1 column; DataSheet 2 model code); 'disease duration' in Table 1 and Table 2"
    ),
    CRP = list(
      description = "Study-arm mean baseline C-reactive protein, standard assay.",
      units = "mg/dL (see notes)",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, baseline only. Enters as (CRP / 1.786)^e_crp_emax on every drug effect (DataSheet 2 Equation 12); 1.786 is the mean arm baseline CRP over the 'CRP' sheet of the deposited analysis dataset with missing arm values replaced by the median arm value. UNITS: Methods 'Data Development' states CRP was standardised to mg/L, but every value in Table 1 (arm means 0.67-2.98) and in the dataset is on the mg/dL scale typical of moderate-to-severe Crohn's disease trials; supply CRP in the same unit as the reference value. The power form is scale-free only when column and reference share a unit. The paper's typical simulated trial has a baseline CRP of 1.67.",
      source_name = "CRP (DataSheet 1 column; DataSheet 2 model code)"
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
      source_name = "DOSE with FLAG = infliximab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_NATALIZUMAB_DOSE = list(
      description = "Natalizumab dose in the study arm; 0 if the arm did not receive natalizumab. Drug class: integrin-alpha4 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of natalizumab ('(FLAG == k)' in the model code), so only CONMED_NATALIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210 or 420 mg (3 or 6 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = natalizumab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CDP571_DOSE = list(
      description = "CDP571 dose in the study arm; 0 if the arm did not receive CDP571. Drug class: TNF-alpha inhibitor (humanised anti-TNF antibody; development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of CDP571 ('(FLAG == k)' in the model code), so only CONMED_CDP571_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 700 mg (10 mg/kg at 70 kg); the IBDQ sheet records 10 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = CDP571 (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol dose in the study arm; 0 if the arm did not receive certolizumab pegol. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_certolizumab * CONMED_CERTOLIZUMAB_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: 100, 200 or 400 mg per subcutaneous administration (q2w or q4w regimens). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = certolizumab pegol (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab dose in the study arm; 0 if the arm did not receive adalimumab. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of adalimumab ('(FLAG == k)' in the model code), so only CONMED_ADALIMUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: the week-0 induction dose, 40, 80 or 160 mg (regimens 40/20, 80/40 and 160/80 mg at weeks 0/2). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = adalimumab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VEDOLIZUMAB_DOSE = list(
      description = "Vedolizumab dose in the study arm; 0 if the arm did not receive vedolizumab. Drug class: integrin-alpha4beta7 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vedolizumab ('(FLAG == k)' in the model code), so only CONMED_VEDOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 35 or 140 mg (0.5 or 2 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vedolizumab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description = "Ustekinumab dose in the study arm; 0 if the arm did not receive ustekinumab. Drug class: IL-12/23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of ustekinumab ('(FLAG == k)' in the model code), so only CONMED_USTEKINUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 90 mg per subcutaneous administration (the CDAI-100 and CRP sheets record the same arm as 1.5). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = ustekinumab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_SEMAPIMOD_DOSE = list(
      description = "Semapimod dose in the study arm; 0 if the arm did not receive semapimod. Drug class: TNF-alpha inhibitor group of Table 1.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of semapimod ('(FLAG == k)' in the model code), so only CONMED_SEMAPIMOD_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 180 mg (60 mg i.v. daily for 3 days). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = semapimod (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose in the study arm; 0 if the arm did not receive tofacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of tofacitinib ('(FLAG == k)' in the model code), so only CONMED_TOFACITINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 1, 5, 10 or 15 mg per administration of a twice-daily regimen. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = tofacitinib (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ABATACEPT_DOSE = list(
      description = "Abatacept dose in the study arm; 0 if the arm did not receive abatacept. Drug class: T-cell activation inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of abatacept ('(FLAG == k)' in the model code), so only CONMED_ABATACEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210, 700 or 2100 mg (3, 10 or 30 mg/kg at 70 kg); the IBDQ sheet records 3, 10 or 30 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = abatacept (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_PF04236921_DOSE = list(
      description = "PF-04236921 dose in the study arm; 0 if the arm did not receive PF-04236921. Drug class: IL-6 inhibitor (development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: Emax dose-response, effect = emax_pf04236921 * DOSE / (DOSE + ED50) (model code), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: 10, 50 or 200 mg per subcutaneous administration (weeks 0 and 4). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = PF-04236921 (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_UPADACITINIB_DOSE = list(
      description = "Upadacitinib dose in the study arm; 0 if the arm did not receive upadacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_upadacitinib * CONMED_UPADACITINIB_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: 3, 6, 12 or 24 mg per administration of a twice-daily regimen; the 24 mg once-daily arm is recorded as 12, its twice-daily equivalent. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = upadacitinib (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FONTOLIZUMAB_DOSE = list(
      description = "Fontolizumab dose in the study arm; 0 if the arm did not receive fontolizumab. Drug class: IFN-gamma inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of fontolizumab ('(FLAG == k)' in the model code), so only CONMED_FONTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 70, 280 or 700 mg (1, 4 or 10 mg/kg loading dose at 70 kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = fontolizumab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VERCIRNON_DOSE = list(
      description = "Vercirnon dose in the study arm; 0 if the arm did not receive vercirnon. Drug class: CCR9 antagonist.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vercirnon ('(FLAG == k)' in the model code), so only CONMED_VERCIRNON_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: the total daily dose, 250, 500 or 1000 mg (250 mg bid is recorded as 500). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vercirnon (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_LAQUINIMOD_DOSE = list(
      description = "Laquinimod dose in the study arm; 0 if the arm did not receive laquinimod. Drug class: T-cell activation inhibitor group of Table 1.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of laquinimod ('(FLAG == k)' in the model code), so only CONMED_LAQUINIMOD_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 0.5, 1, 1.5 or 2 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = laquinimod (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_RISANKIZUMAB_DOSE = list(
      description = "Risankizumab dose in the study arm; 0 if the arm did not receive risankizumab. Drug class: IL-23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of risankizumab ('(FLAG == k)' in the model code), so only CONMED_RISANKIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 200 or 600 mg per intravenous infusion (q4w). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = risankizumab (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ONTAMALIMAB_DOSE = list(
      description = "Ontamalimab (PF-00547659) dose in the study arm; 0 if the arm did not receive ontamalimab (PF-00547659). Drug class: MAdCAM inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of ontamalimab (PF-00547659) ('(FLAG == k)' in the model code), so only CONMED_ONTAMALIMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 22.5, 75 or 225 mg per subcutaneous administration (q4w). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = ontamalimab (PF-00547659) (DataSheet 1 'CRP' sheet); Table 1 'Route (regimen)' column"
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
    SCORE_CDAI = list(
      description = "Study-arm mean baseline CDAI.",
      units = "(score, 0-600)",
      type = "continuous",
      notes = "Screened but not retained in the Delta CRP and Delta IBDQ models; see the Yu_2022_crohns_mbma vignette."
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
    notes = "MBMA at the STUDY-ARM level: 46 trials, 146 treatment arms and 12,846 patients in total (Results 'Available Data'), of which 26 trials report change from baseline in CRP (CRP sheet of the deposited dataset). Female percentage is 100 minus the Table 1 'Total' percentage of men (46.09). The model has no between-study random effect: the source absorbs between-trial differences into its non-parametric placebo, one estimate per trial and visit."
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
    e0_d_crp <- fixed(0.016)
    label("Placebo change from baseline at the typical simulated trial, Week 12 (mg/dL)")
    # Results 'Model Simulation': 'For CRP, the placebo effect was simulated
    # as 0.016 with a longitudinal placebo model'.
    # The fitted model estimates one placebo value per trial and visit (gnls
    # params 'eo ~ -1 + I(group.week)') and the paper reports none of them.
    # The Week-12 simulation value is held constant in time here; replace it
    # with the observed placebo response of the trial being simulated.

    # ---- Covariate effects (power form, Equation 12) -----------------------
    e_age_emax <- -7.69
    label("Power exponent on age over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 6 and main-text Table 2 row 'Covariate: Age' = -7.69
    # (95% CI -11.11, -4.26).
    e_t_diag_cd_emax <- 4.95
    label("Power exponent on disease duration over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 6 and main-text Table 2 row 'Covariate: Disease
    # duration' = 4.95 (95% CI 3.60, 6.29).
    e_crp_emax <- -0.87
    label("Power exponent on baseline CRP over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 6 and main-text Table 2 row 'Covariate: Baseline
    # CRP' = -0.87 (95% CI -1.24, -0.50).

    # ---- Drug effects (Supplementary Table 6, 'Emax' block) ----
    emax_infliximab <- -0.61
    label("Effect of infliximab at any dose (mg/dL)")
    # Supplementary Table 6, 'Infliximab' = -0.61 (95% CI -0.98, -0.24).
    emax_natalizumab <- -0.89
    label("Effect of natalizumab at any dose (mg/dL)")
    # Supplementary Table 6, 'Natalizumab' = -0.89 (95% CI -1.12, -0.65).
    emax_cdp571 <- -0.15
    label("Effect of CDP571 at any dose (mg/dL)")
    # Supplementary Table 6, 'CDP571' = -0.15 (95% CI -0.28, -0.02).
    slope_certolizumab <- -1.12e-3
    label("Linear dose-response slope of certolizumab pegol (mg/dL per mg)")
    # PRINTED '1.12x10-2' with 95% CI (-1.70x10-2, -5.47x10-4) in both
    # Supplementary Table 6 and main-text Table 2. The point estimate lies
    # outside its own interval, so the printed row is internally inconsistent.
    # The encoded -1.12x10-3 is the only reading that is (a) the exact
    # midpoint of a symmetric interval (-1.70x10-3, -5.47x10-4), i.e. the same
    # two exponent slips in the estimate and the lower bound plus a dropped
    # minus sign, and (b) consistent with Figure 5B, whose certolizumab
    # 100/200/400 mg points sit at about -0.10/-0.21/-0.46. The maintainers'
    # refit of the authors' deposited code and data gives -1.07x10-3.
    emax_adalimumab <- -0.47
    label("Effect of adalimumab at any dose (mg/dL)")
    # Supplementary Table 6, 'Adalimumab' = -0.47 (95% CI -0.67, -0.27).
    emax_vedolizumab <- 0.03
    label("Effect of vedolizumab at any dose (mg/dL)")
    # Supplementary Table 6, 'Vedolizumab' = 0.03 (95% CI -0.45, 0.52).
    emax_ustekinumab <- -0.22
    label("Effect of ustekinumab at any dose (mg/dL)")
    # Supplementary Table 6, 'Ustekinumab' = -0.22 (95% CI -0.54, 0.10).
    emax_semapimod <- -0.07
    label("Effect of semapimod at any dose (mg/dL)")
    # Supplementary Table 6, 'Semapimod' = -0.07 (95% CI -1.80, 1.66).
    emax_tofacitinib <- -0.07
    label("Effect of tofacitinib at any dose (mg/dL)")
    # Supplementary Table 6, 'Tofacitinib' = -0.07 (95% CI -0.10, -0.04).
    emax_abatacept <- 0.59
    label("Effect of abatacept at any dose (mg/dL)")
    # Supplementary Table 6, 'Abatacept' = 0.59 (95% CI 0.32, 0.87).
    emax_pf04236921 <- -7.47
    label("Maximum effect of PF-04236921 in the Emax dose-response (mg/dL)")
    # Supplementary Table 6 and main-text Table 2, 'PF-04236921 (Emax)' =
    # -7.47 (95% CI -12.76, -2.18).
    led50_pf04236921 <- log(93.69)
    label("Log dose giving half the maximum effect of PF-04236921 (log mg)")
    # Supplementary Table 6 and main-text Table 2, 'PF-04236921 (ED50)' =
    # 93.69 (95% CI 41.26, 241.86); the interval is symmetric on the log
    # scale, matching the code's exp(led50) parameterisation.
    slope_upadacitinib <- -7.3e-3
    label("Linear dose-response slope of upadacitinib (mg/dL per mg)")
    # PRINTED -0.22 (95% CI -0.43, -0.01) in Supplementary Table 6 and main-
    # text Table 2, which would give -5.3 at 24 mg and make upadacitinib the
    # most effective CRP-lowering regimen; the Results name PF-04236921 200 mg
    # (-5.52) as the most effective and Figure 5B plots upadacitinib 24/12/6/3
    # mg at about -0.17/-0.07/-0.03/0.00. The printed value is therefore wrong
    # by about 30-fold. The encoded -7.3x10-3 per mg is back-solved by the
    # maintainers from Figure 5B (the placebo-corrected 24 mg point relative
    # to the natalizumab, infliximab and PF-04236921 points of the same panel,
    # which cancels the covariate factor); the maintainers' refit of the
    # authors' deposited code and data gives -8.0x10-3.
    emax_fontolizumab <- -0.13
    label("Effect of fontolizumab at any dose (mg/dL)")
    # Supplementary Table 6, 'Fontolizumab' = -0.13 (95% CI -0.28, 0.03).
    emax_vercirnon <- 0.1
    label("Effect of vercirnon at any dose (mg/dL)")
    # Supplementary Table 6, 'Vercirnon' = 0.10 (95% CI -0.05, 0.25).
    emax_laquinimod <- 0.25
    label("Effect of laquinimod at any dose (mg/dL)")
    # Supplementary Table 6, 'Laquinimod' = 0.25 (95% CI -0.01, 0.50).
    emax_risankizumab <- -4.6e-3
    label("Effect of risankizumab at any dose (mg/dL)")
    # PRINTED '4.60x10-3' with 95% CI (5.70x10-2, 4.78x10-2), an interval that
    # excludes its estimate. The encoded -4.60x10-3 is the midpoint of the
    # symmetric interval (-5.70x10-2, 4.78x10-2), i.e. two dropped minus
    # signs; Figure 5B places risankizumab at the placebo line and the
    # maintainers' refit gives -7.0x10-3.
    emax_ontamalimab <- -0.05
    label("Effect of ontamalimab (PF-00547659) at any dose (mg/dL)")
    # Supplementary Table 6, 'Ontamalimab (PF-00547659)' = -0.05 (95% CI
    # -0.12, 0.02).

    # ---- Residual error ------------------------------------------------------
    addSd <- fixed(1.81)
    label("Between-patient SD of the change; arm-mean residual SD is addSd/sqrt(N_ARM) (mg/dL)")
    # NOT REPORTED in the paper. The source weights each arm by its own
    # reported SD/sqrt(N) (varFixed(~SD^2/no)) with an estimated scale
    # factor. Encoded at unit scale with the median arm SD of the 'CRP'
    # sheet of the deposited dataset (1.81); supply an arm's own SD by
    # overriding addSd. The source also fits a within-arm residual
    # autocorrelation (ARMA(2)), which a per-record rxode2 residual cannot
    # carry.
  })

  model({
    # Equation 12: every drug effect is scaled by the arm's covariates.
    covEff <- (AGE / 36.945)^e_age_emax *
      (T_DIAG_CD / 9.317)^e_t_diag_cd_emax *
      (CRP / 1.786)^e_crp_emax

    # One term per drug. An arm supplies a positive dose in exactly one
    # CONMED_<drug>_DOSE column and zero in the rest, so a placebo arm (all
    # columns zero) reduces to the placebo term alone.
    edInfliximab <- emax_infliximab * (CONMED_INFLIXIMAB_DOSE > 0)
    edNatalizumab <- emax_natalizumab * (CONMED_NATALIZUMAB_DOSE > 0)
    edCdp571 <- emax_cdp571 * (CONMED_CDP571_DOSE > 0)
    edCertolizumab <- slope_certolizumab * CONMED_CERTOLIZUMAB_DOSE
    edAdalimumab <- emax_adalimumab * (CONMED_ADALIMUMAB_DOSE > 0)
    edVedolizumab <- emax_vedolizumab * (CONMED_VEDOLIZUMAB_DOSE > 0)
    edUstekinumab <- emax_ustekinumab * (CONMED_USTEKINUMAB_DOSE > 0)
    edSemapimod <- emax_semapimod * (CONMED_SEMAPIMOD_DOSE > 0)
    edTofacitinib <- emax_tofacitinib * (CONMED_TOFACITINIB_DOSE > 0)
    edAbatacept <- emax_abatacept * (CONMED_ABATACEPT_DOSE > 0)
    edPf04236921 <- emax_pf04236921 *
      CONMED_PF04236921_DOSE / (CONMED_PF04236921_DOSE + exp(led50_pf04236921))
    edUpadacitinib <- slope_upadacitinib * CONMED_UPADACITINIB_DOSE
    edFontolizumab <- emax_fontolizumab * (CONMED_FONTOLIZUMAB_DOSE > 0)
    edVercirnon <- emax_vercirnon * (CONMED_VERCIRNON_DOSE > 0)
    edLaquinimod <- emax_laquinimod * (CONMED_LAQUINIMOD_DOSE > 0)
    edRisankizumab <- emax_risankizumab * (CONMED_RISANKIZUMAB_DOSE > 0)
    edOntamalimab <- emax_ontamalimab * (CONMED_ONTAMALIMAB_DOSE > 0)

    edrug <- (edInfliximab + edNatalizumab + edCdp571 + edCertolizumab +
      edAdalimumab + edVedolizumab + edUstekinumab + edSemapimod +
      edTofacitinib + edAbatacept + edPf04236921 + edUpadacitinib +
      edFontolizumab + edVercirnon + edLaquinimod + edRisankizumab +
      edOntamalimab) * covEff

    # Equation 7: placebo plus drug term on the natural scale.
    d_crp <- e0_d_crp + edrug

    # Equation 8: the residual SD of an arm mean scales with 1/sqrt(N).
    sdArm <- addSd / sqrt(N_ARM)
    d_crp ~ add(sdArm)
  })
}
