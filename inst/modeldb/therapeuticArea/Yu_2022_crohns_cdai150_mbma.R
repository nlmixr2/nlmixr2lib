Yu_2022_crohns_cdai150_mbma <- function() {
  description <- paste0(
    "MBMA. Model-based meta-analysis of the proportion of patients in ",
    "clinical remission (CDAI150: an absolute Crohn's Disease Activity ",
    "Index score below 150) in adults with moderate-to-severe Crohn's ",
    "disease, fitted to study-arm summary data from 38 double-blind ",
    "randomised controlled trials (induction period) of biologics and small ",
    "targeted molecules published up to March 2020. The logit of the arm ",
    "CDAI150 clinical-remission rate is the sum of a placebo term and one ",
    "drug term per treatment, P = expit(E0 + Edrug), so every drug effect ",
    "is a shift on the log-odds scale. Each of 22 drugs carries its own ",
    "maximum effect; the dose-response drugs are adalimumab (linear slope ",
    "per mg), risankizumab (linear slope per mg) and PF-04236921 (linear ",
    "slope per mg). The JAK-inhibitor effects (tofacitinib, filgotinib, ",
    "upadacitinib) rise as 1 - exp(-k_JAK * t); all other drug effects are ",
    "constant in time. Arm-level baseline CDAI (SCORE_CDAI) and baseline ",
    "CRP (CRP) act as power functions (covariate / dataset mean)^theta ",
    "multiplying every drug effect. The authors estimated the placebo ",
    "response non-parametrically (one value per trial per visit); this ",
    "model holds it at the published Week-12 typical-trial placebo value, ",
    "which a user should replace with the placebo response of the trial ",
    "being simulated. Dose enters through one CONMED_<drug>_DOSE covariate ",
    "column per drug (0 = not received); there is no PK layer and no rxode2 ",
    "dose event. Simulation scope is STUDY-ARM-MEAN outcomes, NOT ",
    "individual patients. Companion models of the same paper: ",
    "Yu_2022_crohns_cdai100_mbma, Yu_2022_crohns_cdai70_mbma, ",
    "Yu_2022_crohns_dcdai_mbma, Yu_2022_crohns_dcrp_mbma, ",
    "Yu_2022_crohns_dibdq_mbma."
  )
  reference <- paste0(
    "Yu B, Zhao L, Jin S, He H, Zhang J, Wang X. Model-Based Meta-Analysis ",
    "on the Efficacy of Biologics and Small Targeted Molecules for Crohn's ",
    "Disease. Front Immunol. 2022;13:828219. doi:10.3389/fimmu.2022.828219. ",
    "PMC8967940. Structural model: Methods 'Model Development' Equations ",
    "1-3 and Supplementary Materials (DataSheet 2) Equations 1-12. ",
    "Parameter estimates: Supplementary Table 2 (DataSheet 2, 'Final ",
    "Models'); key parameters also in main-text Table 2. Model code: ",
    "DataSheet 2 'Model Code' (R nlme::gnls). Analysis dataset: ",
    "Supplementary DataSheet 1 (sheet 'CDAI150')."
  )
  vignette <- "Yu_2022_crohns_mbma"

  # No PK layer, no concentration and no rxode2 dose events: dose is a
  # per-arm covariate column, as in He_2021_psoriasis_pasi75_mbma (same
  # research group and the same MBMA methodology).
  units <- list(
    time = "week (weeks since the first dose of the induction period; the onset rate constants are in 1/week and the paper simulates Week 12)",
    dosing = "mg (the dose metric of each CONMED_<drug>_DOSE column is described in its covariateData notes; weight-based regimens are expressed at 70 kg per Methods 'Data Development'. This model consumes NO rxode2 dose events.)",
    concentration = "probability/arm (prob_cdai150 is the STUDY-ARM proportion of patients reaching the endpoint, on a 0-1 scale; it is NOT a drug concentration)"
  )

  covariateData <- list(
    SCORE_CDAI = list(
      description = "Study-arm mean baseline Crohn's Disease Activity Index score.",
      units = "(score, 0-600)",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, not subject-level. Enters as (SCORE_CDAI / 303.19)^e_score_cdai_emax on every drug effect (DataSheet 2 Equation 12). 303.19 is the mean of the arm baseline CDAI over the 'CDAI150' sheet of the deposited analysis dataset, with missing arm values replaced by the median arm value as Methods 'Available Data' describes (the model code normalises by mean(<sheet>$CDAI)). The paper's typical simulated trial has a baseline CDAI of 306.75 (Results 'Model Simulation'; Table 1 'Total' row).",
      source_name = "CDAI (DataSheet 1 column; DataSheet 2 model code)"
    ),
    CRP = list(
      description = "Study-arm mean baseline C-reactive protein, standard assay.",
      units = "mg/dL (see notes)",
      type = "continuous",
      reference_category = NULL,
      notes = "TRIAL-ARM-LEVEL, baseline only. Enters as (CRP / 1.688)^e_crp_emax on every drug effect (DataSheet 2 Equation 12); 1.688 is the mean arm baseline CRP over the 'CDAI150' sheet of the deposited analysis dataset with missing arm values replaced by the median arm value. UNITS: Methods 'Data Development' states CRP was standardised to mg/L, but every value in Table 1 (arm means 0.67-2.98) and in the dataset is on the mg/dL scale typical of moderate-to-severe Crohn's disease trials; supply CRP in the same unit as the reference value. The power form is scale-free only when column and reference share a unit. The paper's typical simulated trial has a baseline CRP of 1.67.",
      source_name = "CRP (DataSheet 1 column; DataSheet 2 model code)"
    ),
    N_ARM = list(
      description = "Number of patients contributing to the study-arm proportion at a given visit.",
      units = "participants",
      type = "count",
      reference_category = NULL,
      notes = "Study-design quantity supplied per observation row; the meta-analytic weight of the residual. The source weights each arm by the binomial standard error sqrt(P*(1 - P)/N) (DataSheet 2 Equation 11; model code varPower(fixed = 0.5) on fitted*(1 - fitted)/no). It does not affect the typical-value prediction.",
      source_name = "no (DataSheet 1 column); N in DataSheet 2 Equations 8 and 11"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description = "Infliximab dose in the study arm; 0 if the arm did not receive infliximab. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of infliximab ('(FLAG == k)' in the model code), so only CONMED_INFLIXIMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 350 mg (5 mg/kg at the 70 kg per-patient normalisation of Methods 'Data Development'); the IBDQ sheet of the deposited dataset records the same arms as 5 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = infliximab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_NATALIZUMAB_DOSE = list(
      description = "Natalizumab dose in the study arm; 0 if the arm did not receive natalizumab. Drug class: integrin-alpha4 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of natalizumab ('(FLAG == k)' in the model code), so only CONMED_NATALIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 210 or 420 mg (3 or 6 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = natalizumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CDP571_DOSE = list(
      description = "CDP571 dose in the study arm; 0 if the arm did not receive CDP571. Drug class: TNF-alpha inhibitor (humanised anti-TNF antibody; development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of CDP571 ('(FLAG == k)' in the model code), so only CONMED_CDP571_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 700 mg (10 mg/kg at 70 kg); the IBDQ sheet records 10 (mg/kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = CDP571 (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description = "Etanercept dose in the study arm; 0 if the arm did not receive etanercept. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of etanercept ('(FLAG == k)' in the model code), so only CONMED_ETANERCEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 25 mg per subcutaneous administration (25 mg twice weekly). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = etanercept (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description = "Certolizumab pegol dose in the study arm; 0 if the arm did not receive certolizumab pegol. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of certolizumab pegol ('(FLAG == k)' in the model code), so only CONMED_CERTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 100, 200 or 400 mg per subcutaneous administration (q2w or q4w regimens). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = certolizumab pegol (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ADALIMUMAB_DOSE = list(
      description = "Adalimumab dose in the study arm; 0 if the arm did not receive adalimumab. Drug class: TNF-alpha inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_adalimumab * CONMED_ADALIMUMAB_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: the week-0 induction dose, 40, 80 or 160 mg (regimens 40/20, 80/40 and 160/80 mg at weeks 0/2). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = adalimumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ONERCEPT_DOSE = list(
      description = "Onercept dose in the study arm; 0 if the arm did not receive onercept. Drug class: TNF-alpha inhibitor (soluble p55 TNF receptor).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of onercept ('(FLAG == k)' in the model code), so only CONMED_ONERCEPT_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 10, 25, 35 or 50 mg per subcutaneous administration (three times weekly). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = onercept (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VEDOLIZUMAB_DOSE = list(
      description = "Vedolizumab dose in the study arm; 0 if the arm did not receive vedolizumab. Drug class: integrin-alpha4beta7 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vedolizumab ('(FLAG == k)' in the model code), so only CONMED_VEDOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 35 or 140 mg (0.5 or 2 mg/kg at 70 kg) or a flat 300 mg per infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vedolizumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_USTEKINUMAB_DOSE = list(
      description = "Ustekinumab dose in the study arm; 0 if the arm did not receive ustekinumab. Drug class: IL-12/23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of ustekinumab ('(FLAG == k)' in the model code), so only CONMED_USTEKINUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 90 mg per subcutaneous administration (the CDAI-100 and CRP sheets record the same arm as 1.5). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = ustekinumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_APILIMOD_DOSE = list(
      description = "Apilimod dose in the study arm; 0 if the arm did not receive apilimod. Drug class: IL-12/23 inhibitor (oral small molecule).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of apilimod ('(FLAG == k)' in the model code), so only CONMED_APILIMOD_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 50 or 100 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = apilimod (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ANDECALIXIMAB_DOSE = list(
      description = "Andecaliximab dose in the study arm; 0 if the arm did not receive andecaliximab. Drug class: MMP-9 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of andecaliximab ('(FLAG == k)' in the model code), so only CONMED_ANDECALIXIMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 75 (150 mg q2w), 150 or 300 mg (150 or 300 mg qw), i.e. the weekly-equivalent dose. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = andecaliximab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_TOFACITINIB_DOSE = list(
      description = "Tofacitinib dose in the study arm; 0 if the arm did not receive tofacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of tofacitinib ('(FLAG == k)' in the model code), so only CONMED_TOFACITINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 1, 5, 10 or 15 mg per administration of a twice-daily regimen. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = tofacitinib (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FILGOTINIB_DOSE = list(
      description = "Filgotinib dose in the study arm; 0 if the arm did not receive filgotinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of filgotinib ('(FLAG == k)' in the model code), so only CONMED_FILGOTINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 200 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = filgotinib (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_RISANKIZUMAB_DOSE = list(
      description = "Risankizumab dose in the study arm; 0 if the arm did not receive risankizumab. Drug class: IL-23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_risankizumab * CONMED_RISANKIZUMAB_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: 200 or 600 mg per intravenous infusion (q4w). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = risankizumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_PF04236921_DOSE = list(
      description = "PF-04236921 dose in the study arm; 0 if the arm did not receive PF-04236921. Drug class: IL-6 inhibitor (development code, no INN).",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "GRADED: linear dose-response, effect = slope_pf04236921 * CONMED_PF04236921_DOSE (model code 'k*DOSE'), so the dose metric matters. Dose metric of the deposited dataset's DOSE column for this drug: 10, 50 or 200 mg per subcutaneous administration (weeks 0 and 4). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = PF-04236921 (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_BRAZIKUMAB_DOSE = list(
      description = "Brazikumab dose in the study arm; 0 if the arm did not receive brazikumab. Drug class: IL-23 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of brazikumab ('(FLAG == k)' in the model code), so only CONMED_BRAZIKUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 700 mg per intravenous infusion. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = brazikumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_UPADACITINIB_DOSE = list(
      description = "Upadacitinib dose in the study arm; 0 if the arm did not receive upadacitinib. Drug class: JAK inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of upadacitinib ('(FLAG == k)' in the model code), so only CONMED_UPADACITINIB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 3, 6, 12 or 24 mg per administration of a twice-daily regimen; the 24 mg once-daily arm is recorded as 12, its twice-daily equivalent. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = upadacitinib (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_FONTOLIZUMAB_DOSE = list(
      description = "Fontolizumab dose in the study arm; 0 if the arm did not receive fontolizumab. Drug class: IFN-gamma inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of fontolizumab ('(FLAG == k)' in the model code), so only CONMED_FONTOLIZUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 70, 280 or 700 mg (1, 4 or 10 mg/kg loading dose at 70 kg). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = fontolizumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ABRILUMAB_DOSE = list(
      description = "Abrilumab dose in the study arm; 0 if the arm did not receive abrilumab. Drug class: integrin-alpha4beta7 inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of abrilumab ('(FLAG == k)' in the model code), so only CONMED_ABRILUMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 21, 70 or 210 mg per subcutaneous administration. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = abrilumab (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_VERCIRNON_DOSE = list(
      description = "Vercirnon dose in the study arm; 0 if the arm did not receive vercirnon. Drug class: CCR9 antagonist.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of vercirnon ('(FLAG == k)' in the model code), so only CONMED_VERCIRNON_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: the total daily dose, 250, 500 or 1000 mg (250 mg bid is recorded as 500). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = vercirnon (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_LAQUINIMOD_DOSE = list(
      description = "Laquinimod dose in the study arm; 0 if the arm did not receive laquinimod. Drug class: T-cell activation inhibitor group of Table 1.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of laquinimod ('(FLAG == k)' in the model code), so only CONMED_LAQUINIMOD_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 0.5, 1, 1.5 or 2 mg once daily. Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = laquinimod (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
    ),
    CONMED_ONTAMALIMAB_DOSE = list(
      description = "Ontamalimab (PF-00547659) dose in the study arm; 0 if the arm did not receive ontamalimab (PF-00547659). Drug class: MAdCAM inhibitor.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "STEP, NOT GRADED, in this model: the source estimates one effect for every regimen of ontamalimab (PF-00547659) ('(FLAG == k)' in the model code), so only CONMED_ONTAMALIMAB_DOSE > 0 is read and any positive value gives the same effect. Dose metric of the deposited dataset's DOSE column for this drug: 22.5, 75 or 225 mg per subcutaneous administration (q4w). Methods 'Data Development' describes normalisation 'by daily dose', which the dataset follows only for vercirnon and andecaliximab; the column carries the dataset's own values.",
      source_name = "DOSE with FLAG = ontamalimab (PF-00547659) (DataSheet 1 'CDAI150' sheet); Table 1 'Route (regimen)' column"
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
    notes = "MBMA at the STUDY-ARM level: 46 trials, 146 treatment arms and 12,846 patients in total (Results 'Available Data'), of which 38 trials report CDAI150 clinical-remission rate (CDAI150 sheet of the deposited dataset). Female percentage is 100 minus the Table 1 'Total' percentage of men (46.09). The model has no between-study random effect: the source absorbs between-trial differences into its non-parametric placebo, one estimate per trial and visit."
  )

  ini({
    # ========================================================================
    # Structure (Methods 'Model Development'; DataSheet 2 Equations 1-12 and
    # 'Code for ... model'). The published gnls() code is the authority:
    #
    #   P_ijt   = expit(E0_it + Edrug_ijt)                 (Equations 9-10)
    #   Edrug   = Emax_drug * f(time) * prod((X / mean(X))^theta)   (Eq. 2, 3, 12)
    #   Emax_drug = constant, slope * DOSE, or Emax * DOSE / (DOSE + ED50)
    #   (Supplementary table footnotes a and b; the footnote typesets the
    #   Emax denominator as 'ED50 * dose', the code as DOSE + exp(led50)).
    #
    # With the values below every study-arm prediction of the paper's own
    # Week-12 ranking figure is reproduced; see the vignette.
    # ========================================================================

    # ---- Placebo (the source estimates it non-parametrically) --------------
    e0_cdai150 <- fixed(-1.3093)
    label("Placebo response on the logit scale at the typical simulated trial, Week 12 (log-odds)")
    # Results 'Model Simulation': 'The model simulation of the CDAI150, with a
    # placebo effect estimated as 21.26%'. logit(0.2126) = -1.3093.
    # The fitted model estimates one placebo value per trial and visit (gnls
    # params 'eo ~ -1 + I(group.week)') and the paper reports none of them.
    # The Week-12 simulation value is held constant in time here; replace it
    # with the observed placebo response of the trial being simulated.

    # ---- Onset of the drug effect -------------------------------------------
    lkdrug_jak <- log(0.11)
    label("Log onset rate constant of the JAK-inhibitor drug effect (log 1/week)")
    # Supplementary Table 2 row 'kJAK, Rate constant for the onset of JAK
    # inhibitor' = 0.11 (95% CI 0.01, 1.40); also main-text Table 2. Results:
    # 'ET50 of JAK inhibitors was estimated to be about 6.3 weeks ... ET90 ...
    # 20.9 weeks' = log(2)/0.11 and log(10)/0.11. The code writes the rate as
    # exp(k), so the printed 0.11 is the rate itself.

    # ---- Covariate effects (power form, Equation 12) -----------------------
    e_score_cdai_emax <- -5.22
    label("Power exponent on baseline CDAI over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 2 and main-text Table 2 row 'Covariate: Baseline
    # CDAI' = -5.22 (95% CI -7.90, -2.54).
    e_crp_emax <- 0.51
    label("Power exponent on baseline CRP over its dataset mean, on every drug effect (unitless)")
    # Supplementary Table 2 and main-text Table 2 row 'Covariate: Baseline
    # CRP' = 0.51 (95% CI 0.20, 0.83).

    # ---- Drug effects (Supplementary Table 2, 'Emax' block) ----
    emax_infliximab <- 1.16
    label("Effect of infliximab at any dose (log-odds)")
    # Supplementary Table 2, 'Infliximab' = 1.16 (95% CI 0.87, 1.46).
    emax_natalizumab <- 0.44
    label("Effect of natalizumab at any dose (log-odds)")
    # Supplementary Table 2, 'Natalizumab' = 0.44 (95% CI 0.31, 0.57).
    emax_cdp571 <- 0.18
    label("Effect of CDP571 at any dose (log-odds)")
    # Supplementary Table 2, 'CDP571' = 0.18 (95% CI -0.14, 0.50).
    emax_etanercept <- -0.7
    label("Effect of etanercept at any dose (log-odds)")
    # Supplementary Table 2, 'Etanercept' = -0.70 (95% CI -1.68, 0.27).
    emax_certolizumab <- 0.47
    label("Effect of certolizumab pegol at any dose (log-odds)")
    # Supplementary Table 2, 'Certolizumab pegol' = 0.47 (95% CI 0.31, 0.63).
    slope_adalimumab <- 6e-3
    label("Linear dose-response slope of adalimumab (log-odds per mg)")
    # Supplementary Table 2 and main-text Table 2, 'Adalimumab (slope)' =
    # 6.00x10^-3 (95% CI 4.22x10^-3, 7.78x10^-3).
    emax_onercept <- 0.18
    label("Effect of onercept at any dose (log-odds)")
    # Supplementary Table 2, 'Onercept' = 0.18 (95% CI -0.56, 0.92).
    emax_vedolizumab <- 0.7
    label("Effect of vedolizumab at any dose (log-odds)")
    # Supplementary Table 2, 'Vedolizumab' = 0.70 (95% CI 0.45, 0.95).
    emax_ustekinumab <- 0.44
    label("Effect of ustekinumab at any dose (log-odds)")
    # Supplementary Table 2, 'Ustekinumab' = 0.44 (95% CI -0.45, 1.33).
    emax_apilimod <- -0.12
    label("Effect of apilimod at any dose (log-odds)")
    # Supplementary Table 2, 'Apilimod' = -0.12 (95% CI -0.55, 0.31).
    emax_andecaliximab <- -0.26
    label("Effect of andecaliximab at any dose (log-odds)")
    # Supplementary Table 2, 'Andecaliximab' = -0.26 (95% CI -1.40, 0.87).
    emax_tofacitinib <- 0.76
    label("Effect of tofacitinib at any dose (log-odds)")
    # Supplementary Table 2, 'Tofacitinib' = 0.76 (95% CI -1.14, 2.66).
    emax_filgotinib <- 1.5
    label("Effect of filgotinib at any dose (log-odds)")
    # Supplementary Table 2, 'Filgotinib' = 1.50 (95% CI -1.02, 4.01).
    slope_risankizumab <- 2.55e-3
    label("Linear dose-response slope of risankizumab (log-odds per mg)")
    # Supplementary Table 2 and main-text Table 2, 'Risankizumab (slope)' =
    # 2.55x10^-3 (95% CI 8.88x10^-4, 4.22x10^-3).
    slope_pf04236921 <- 8.33e-3
    label("Linear dose-response slope of PF-04236921 (log-odds per mg)")
    # Supplementary Table 2 and main-text Table 2, 'PF-04236921 (slope)' =
    # 8.33x10^-3 (95% CI 2.92x10^-3, 0.01).
    emax_brazikumab <- 0.69
    label("Effect of brazikumab at any dose (log-odds)")
    # Supplementary Table 2, 'Brazikumab' = 0.69 (95% CI -0.04, 1.42).
    emax_upadacitinib <- 1.09
    label("Effect of upadacitinib at any dose (log-odds)")
    # Supplementary Table 2, 'Upadacitinib' = 1.09 (95% CI -0.39, 2.58).
    emax_fontolizumab <- 0.78
    label("Effect of fontolizumab at any dose (log-odds)")
    # Supplementary Table 2, 'Fontolizumab' = 0.78 (95% CI 0.39, 1.16).
    emax_abrilumab <- 0.51
    label("Effect of abrilumab at any dose (log-odds)")
    # Supplementary Table 2, 'Abrilumab' = 0.51 (95% CI -0.03, 1.06).
    emax_vercirnon <- -0.19
    label("Effect of vercirnon at any dose (log-odds)")
    # Supplementary Table 2, 'Vercirnon' = -0.19 (95% CI -0.48, 0.11).
    emax_laquinimod <- 0.88
    label("Effect of laquinimod at any dose (log-odds)")
    # Supplementary Table 2, 'Laquinimod' = 0.88 (95% CI 0.45, 1.32).
    emax_ontamalimab <- 0.45
    label("Effect of ontamalimab (PF-00547659) at any dose (log-odds)")
    # Supplementary Table 2, 'Ontamalimab (PF-00547659)' = 0.45 (95% CI -0.13,
    # 1.02).

    # ---- Residual error ------------------------------------------------------
    addSd <- 0.837
    label("Multiplier on the binomial SE sqrt(P*(1-P)/N_ARM) giving the arm residual SD (unitless)")
    # NOT REPORTED in the paper. Residual standard error of the maintainers'
    # refit of the published gnls() code to the deposited dataset (DataSheet
    # 1, missing arm covariates replaced by the median arm value); the same
    # refit returns every published fixed effect of Supplementary Table 2 to within
    # about 2%. The source also fits a within-arm residual autocorrelation
    # (ARMA(2)), which a per-record rxode2 residual cannot carry.
  })

  model({
    # Equation 12: every drug effect is scaled by the arm's covariates.
    covEff <- (SCORE_CDAI / 303.19)^e_score_cdai_emax *
      (CRP / 1.688)^e_crp_emax

    # Equation 3, JAK inhibitors only (CLASS == 5 in the model code).
    onsetJak <- 1 - exp(-exp(lkdrug_jak) * time)

    # One term per drug. An arm supplies a positive dose in exactly one
    # CONMED_<drug>_DOSE column and zero in the rest, so a placebo arm (all
    # columns zero) reduces to the placebo term alone.
    edInfliximab <- emax_infliximab * (CONMED_INFLIXIMAB_DOSE > 0)
    edNatalizumab <- emax_natalizumab * (CONMED_NATALIZUMAB_DOSE > 0)
    edCdp571 <- emax_cdp571 * (CONMED_CDP571_DOSE > 0)
    edEtanercept <- emax_etanercept * (CONMED_ETANERCEPT_DOSE > 0)
    edCertolizumab <- emax_certolizumab * (CONMED_CERTOLIZUMAB_DOSE > 0)
    edAdalimumab <- slope_adalimumab * CONMED_ADALIMUMAB_DOSE
    edOnercept <- emax_onercept * (CONMED_ONERCEPT_DOSE > 0)
    edVedolizumab <- emax_vedolizumab * (CONMED_VEDOLIZUMAB_DOSE > 0)
    edUstekinumab <- emax_ustekinumab * (CONMED_USTEKINUMAB_DOSE > 0)
    edApilimod <- emax_apilimod * (CONMED_APILIMOD_DOSE > 0)
    edAndecaliximab <- emax_andecaliximab * (CONMED_ANDECALIXIMAB_DOSE > 0)
    edTofacitinib <- emax_tofacitinib * (CONMED_TOFACITINIB_DOSE > 0) * onsetJak
    edFilgotinib <- emax_filgotinib * (CONMED_FILGOTINIB_DOSE > 0) * onsetJak
    edRisankizumab <- slope_risankizumab * CONMED_RISANKIZUMAB_DOSE
    edPf04236921 <- slope_pf04236921 * CONMED_PF04236921_DOSE
    edBrazikumab <- emax_brazikumab * (CONMED_BRAZIKUMAB_DOSE > 0)
    edUpadacitinib <- emax_upadacitinib * (CONMED_UPADACITINIB_DOSE > 0) * onsetJak
    edFontolizumab <- emax_fontolizumab * (CONMED_FONTOLIZUMAB_DOSE > 0)
    edAbrilumab <- emax_abrilumab * (CONMED_ABRILUMAB_DOSE > 0)
    edVercirnon <- emax_vercirnon * (CONMED_VERCIRNON_DOSE > 0)
    edLaquinimod <- emax_laquinimod * (CONMED_LAQUINIMOD_DOSE > 0)
    edOntamalimab <- emax_ontamalimab * (CONMED_ONTAMALIMAB_DOSE > 0)

    edrug <- (edInfliximab + edNatalizumab + edCdp571 + edEtanercept +
      edCertolizumab + edAdalimumab + edOnercept + edVedolizumab +
      edUstekinumab + edApilimod + edAndecaliximab + edTofacitinib +
      edFilgotinib + edRisankizumab + edPf04236921 + edBrazikumab +
      edUpadacitinib + edFontolizumab + edAbrilumab + edVercirnon +
      edLaquinimod + edOntamalimab) * covEff

    # Equations 9-10: inverse logit of the summed placebo and drug terms.
    lp_cdai150 <- e0_cdai150 + edrug
    prob_cdai150 <- expit(lp_cdai150)

    # Equation 11: the residual SD is addSd times the binomial standard error.
    sdArm <- addSd * sqrt(prob_cdai150 * (1 - prob_cdai150) / N_ARM)
    prob_cdai150 ~ add(sdArm)
  })
}
