Selig_2022_piperacillin_mbma <- function() {
  description <- paste(
    "MBMA. One-compartment IV population PK model for piperacillin in",
    "critically ill adults on continuous kidney replacement therapy (CKRT),",
    "fitted by FOCEI in Pumas to concentration-time curves digitised from",
    "20 arms of 10 published CKRT studies together with 14 individual",
    "concentrations from 4 Military Health System patients on CVVH (Selig",
    "2022). Total clearance is the estimated body (non-CKRT) clearance times",
    "exp(eta) plus the CKRT clearance, which is a per-subject data column",
    "(QEFF) rather than a parameter. The random effects are BETWEEN-STUDY",
    "(between-arm) variability of the arm-mean parameters; the source uses",
    "them as between-patient variability in its probability-of-target-",
    "attainment simulations. The proportional residual SD is for a single",
    "patient; the source scaled it by 1/sqrt(N) for an arm mean of N",
    "patients. No covariate was retained. The companion tazobactam model",
    "from the same paper is modellib('Selig_2022_tazobactam_mbma')."
  )
  reference <- paste(
    "Selig DJ, DeLuca JP, Chung KK, Pruskowski KA, Livezey JR, Nadeau RJ,",
    "Por ED, Akers KS. Pharmacokinetics of piperacillin and tazobactam in",
    "critically ill patients treated with continuous kidney replacement",
    "therapy: a mini-review and population pharmacokinetic analysis.",
    "J Clin Pharm Ther. 2022;47(8):1091-1102. doi:10.1111/jcpt.13657.",
    "Parameter values: Table 4A. Structure: Methods 2.3-2.5 (Equations 1-4)",
    "and the supplementary aggregate dataset (JCPT-47-1091-s002.xlsx), whose",
    "CLCKRT column is the per-arm CKRT clearance."
  )
  vignette <- "Selig_2022_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    QEFF = list(
      description = paste(
        "Piperacillin clearance by the CKRT circuit (L/h). For the Military",
        "Health System patients CL_CKRT = Qf * Sc * CF (Selig 2022 Equations",
        "1-3): ultrafiltrate flow times the measured sieving coefficient",
        "Sc = C_effluent / ((C_pre + C_post) / 2) times the prefilter",
        "replacement-fluid correction CF = Qb / (Qb + Qrep). For literature",
        "arms it was taken from the study, or assumed to be 70% of the",
        "reported total CKRT dose when the study did not report it."
      ),
      units = "L/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "A data column, not an estimated parameter: it is added to the",
        "eta-bearing body clearance after the random effect. Set it to 0",
        "for a patient off CKRT. Posted-dataset arms range 0.288-3 L/h",
        "(median 0.755 L/h); the literature median piperacillin CKRT",
        "clearance is 1.43 L/h, range 0.52-2.8 L/h (Table 2). The",
        "source's PTA simulations use 0, 1, 2 and 3 L/h (Figure 3). The",
        "additive placement is not printed as an equation; it is the form",
        "under which the posted dataset reproduces Table 4A (see the",
        "vignette)."
      ),
      source_name = "CLCKRT"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight (study mean or median)",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and V by forward addition (Methods 2.6), not retained. Literature median 78 kg, range 60.3-95.1 (Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance (study mean or median; Cockcroft-Gault where not reported)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened on CL (Methods 2.6), not retained; the source attributes this to CKRT clearing creatinine and to inconsistent timing of the CrCl measurement (Discussion). Literature median 40.91 mL/min (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin (study mean or median)",
      units = "g/L",
      type = "continuous",
      notes = "Screened on V (Methods 2.6), not retained. Reported in g/dL by the source (labelled mg/dL in Table 1); literature median 2.5 g/dL, range 2.11-2.93."
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 11L,
    age_range = "arm means/medians 44-77.7 years in the posted dataset; study medians 54-74 years (Table 1)",
    weight_range = "arm means/medians 52.5-90 kg in the posted dataset; study medians 60.3-95.1 kg (Table 1)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state = paste(
      "Critically ill adults (sepsis, burns, trauma) treated with continuous",
      "kidney replacement therapy (CVVH, CVVHD or CVVHDF). Study-level",
      "APACHE II median 23 (range 21-33.25) and mortality median 44% (range",
      "35-60%) (Table 1)."
    ),
    dose_range = paste(
      "Piperacillin 2000-4000 mg every 6-12 h infused over 0.33-4 h, or",
      "8000-10800 mg/day by continuous infusion (supplementary dataset)."
    ),
    regions = "Multinational literature cohorts plus the United States Military Health System",
    n_concentrations = 152L,
    notes = paste(
      "Aggregate-data model. 10 literature studies contributed 20 arms",
      "(each arm treated as its own trial; 145 patient-arm entries, some",
      "from crossover designs) and the Military Health System contributed",
      "4 individual CVVH patients (2 with burns) dosed 2250 or 3375 mg",
      "piperacillin-tazobactam every 6 h over 30 min with 14 steady-state",
      "concentrations (Methods 2.2; Table S1). 152 observations in total",
      "(Table 4A). Military Health System data are not in the posted dataset."
    )
  )

  ini({
    lcl <- log(2.7); label("Body (non-CKRT) clearance (L/h)") # Table 4A 'CL (L/hr)' 2.7 (RSE 13.48%, 95% CI 1.99-3.41)
    lvc <- log(25.83); label("Central volume of distribution (L)") # Table 4A 'Vc (L)' 25.83 (RSE 7.43%, 95% CI 22.07-29.59)

    # Between-study (between-arm) variability, exponential (Equation 4).
    # Table 4A prints the variances; sqrt(0.38) = 62% CV as quoted in the
    # Discussion. The 0.6 Pearson correlation in Table 4A is between the
    # post hoc etas; no covariance was reported.
    eta_study_lcl ~ 0.38 # Table 4A 'omega2 CL' 0.38 (RSE 32.58%); eta-shrinkage 2.4%
    eta_study_lvc ~ 0.067 # Table 4A 'omega2 Vc' 0.067 (RSE 31.01%); eta-shrinkage 26.8%

    propSd <- 0.42; label("Proportional residual SD for one patient (fraction)") # Table 4A 'Proportional Error' 0.42 (RSE 8.47%); scaled by 1/sqrt(N) for an arm of N patients (Methods 2.5)
  })

  model({
    cl_body <- exp(lcl + eta_study_lcl)
    # CKRT clearance is a per-subject data column, added after the random
    # effect.
    cl <- cl_body + QEFF
    vc <- exp(lvc + eta_study_lvc)
    kel <- cl / vc

    # IV infusion into central; the infusion duration comes from the event
    # table.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
