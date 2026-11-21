Baranda_2024_sar439459_mouse <- function() {
  description <- "Preclinical (mouse, MC38 syngeneic colon tumour). Direct Emax inhibition model of intratumoral active TGF-beta driven by the SAR439459 (anti-TGF-beta monoclonal antibody) plasma concentration, with a study-specific baseline and a shared Emax and EC50 (Baranda 2024)"
  reference <- paste(
    "Baranda JC, Robbrecht D, Sullivan R, Doger B, Santoro A, Barve M,",
    "Grob J-J, Bechter O, Vieito M, de Miguel MJ, Schadendorf D, Johnson M,",
    "Pouzin C, Cantalloube C, Wang R, Lee J, Chen X, Demers B, Amrate A,",
    "Abbadessa G, Hodi FS. Safety, pharmacokinetics, pharmacodynamics, and",
    "antitumor activity of SAR439459, a TGFbeta inhibitor, as monotherapy and",
    "in combination with cemiplimab in patients with advanced solid tumors:",
    "Findings from a phase 1/1b study. Clin Transl Sci. 2024;17(6):e13854.",
    "doi:10.1111/cts.13854",
    sep = " "
  )
  vignette <- "Baranda_2024_sar439459"
  units <- list(
    time = "day",
    dosing = "mg/kg",
    concentration = "ug/mL",
    tgfb_tumor = "pg/mg total protein",
    dosing_notes = "This model has no dosing compartment and no internal PK. The source reports no mouse PK model; the plasma SAR439459 concentration enters only through the time-varying CP_SAR439459_MGL covariate column (ug/mL, the EC50 units of Table S3). The `dosing` entry is nominal; the doses given in the four source studies are recorded in population$dose_range."
  )

  covariateData <- list(
    CP_SAR439459_MGL = list(
      description = "Instantaneous SAR439459 plasma concentration supplied per record as the PD driver",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters the direct Emax inhibition term emax * CP / (ec50 + CP) of supplement Eq. 1. The source analysis used the measured plasma concentration paired with each tumour sample; set to 0 for untreated control mice.",
      source_name = "C_SAR439459"
    ),
    STUDY_MC38_2 = list(
      description = "MC38 study 2 indicator (single dose, 10, 25 and 50 mg/kg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 1)",
      notes = "Selects the study-2 baseline E0 = 44.3 pg/mg protein (Table S3).",
      source_name = NA_character_
    ),
    STUDY_MC38_3 = list(
      description = "MC38 study 3 indicator (three doses of 25 mg/kg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 1)",
      notes = "Selects the study-3 baseline E0 = 58.5 pg/mg protein (Table S3).",
      source_name = NA_character_
    ),
    STUDY_MC38_4 = list(
      description = "MC38 study 4 indicator (single dose, 10, 25 and 50 mg/kg)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 1)",
      notes = "Selects the study-4 baseline E0 = 61.4 pg/mg protein (Table S3).",
      source_name = NA_character_
    )
  )

  population <- list(
    species = "mouse (MC38 syngeneic tumour model)",
    n_subjects = "Not reported",
    n_studies = 4L,
    disease_state = "Mice bearing MC38 colon-carcinoma tumours, each study with an untreated control group",
    dose_range = "SAR439459 1-50 mg/kg, single dose (studies 1, 2 and 4) or three doses (study 3): study 1 1 and 25 mg/kg, study 2 10, 25 and 50 mg/kg, study 3 3 x 25 mg/kg, study 4 10, 25 and 50 mg/kg (Figure S3A legend)",
    notes = "Plasma SAR439459 and intratumoral active TGF-beta were sampled at several times after dosing. The fitted model set the PK concentration thresholds (EC50 about 30 ug/mL, about 90% inhibition near 200 ug/mL) used to choose the 7.5 and 22.5 mg/kg Q3W clinical doses."
  )

  ini({
    # Table S3 (supplement), final PK/PD model estimates.
    le0_tgfb_tumor <- log(52.5); label("Baseline intratumoral active TGF-beta, study 1 (pg/mg protein)") # Table S3 'E0 study 1' = 52.5 (RSE 4.2%)
    le0_tgfb_tumor_study2 <- log(44.3); label("Baseline intratumoral active TGF-beta, study 2 (pg/mg protein)") # Table S3 'E0 study 2' = 44.3 (RSE 2.7%)
    le0_tgfb_tumor_study3 <- log(58.5); label("Baseline intratumoral active TGF-beta, study 3 (pg/mg protein)") # Table S3 'E0 study 3' = 58.5 (RSE 9.4%)
    le0_tgfb_tumor_study4 <- log(61.4); label("Baseline intratumoral active TGF-beta, study 4 (pg/mg protein)") # Table S3 'E0 study 4' = 61.4 (RSE 3.2%)
    lemax <- log(31.7); label("Maximum decrease in intratumoral active TGF-beta (pg/mg protein)") # Table S3 'Emax' = 31.7 (RSE 6.0%)
    lec50 <- log(29.4); label("Plasma SAR439459 concentration at half-maximal decrease (ug/mL)") # Table S3 'EC50' = 29.4 (RSE 25.5%)

    # No random effects are reported in Table S3; typical-value model only.
    propSd <- 0.311; label("Proportional residual error (fraction)") # Table S3 'Proportional error' = 31.1% (RSE 9.8%)
  })
  model({
    e0 <- exp(le0_tgfb_tumor) * (1 - STUDY_MC38_2 - STUDY_MC38_3 - STUDY_MC38_4) +
      exp(le0_tgfb_tumor_study2) * STUDY_MC38_2 +
      exp(le0_tgfb_tumor_study3) * STUDY_MC38_3 +
      exp(le0_tgfb_tumor_study4) * STUDY_MC38_4
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # Supplement Eq. 1: TGFbeta = E0 - Emax * C / (EC50 + C).
    tgfb_tumor <- e0 - emax * CP_SAR439459_MGL / (ec50 + CP_SAR439459_MGL)
    tgfb_tumor ~ prop(propSd)
  })
}
