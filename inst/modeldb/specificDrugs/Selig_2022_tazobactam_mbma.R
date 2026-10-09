Selig_2022_tazobactam_mbma <- function() {
  description <- paste(
    "MBMA. One-compartment IV population PK model for tazobactam in",
    "critically ill adults on continuous kidney replacement therapy (CKRT),",
    "fitted by FOCEI in Pumas to concentration-time curves digitised from",
    "8 published piperacillin-tazobactam and ceftolozane-tazobactam CKRT",
    "studies together with 10 individual concentrations from 3 Military",
    "Health System patients on CVVH (Selig 2022). Total clearance is the",
    "estimated body (non-CKRT) clearance times exp(eta) plus the CKRT",
    "clearance, which is a per-subject data column (QEFF) rather than a",
    "parameter. The random effects are BETWEEN-STUDY (between-arm)",
    "variability of the arm-mean parameters. The proportional residual SD",
    "is for a single patient; the source scaled it by 1/sqrt(N) for an arm",
    "mean of N patients. No covariate was retained. The companion",
    "piperacillin model from the same paper is",
    "modellib('Selig_2022_piperacillin_mbma')."
  )
  reference <- paste(
    "Selig DJ, DeLuca JP, Chung KK, Pruskowski KA, Livezey JR, Nadeau RJ,",
    "Por ED, Akers KS. Pharmacokinetics of piperacillin and tazobactam in",
    "critically ill patients treated with continuous kidney replacement",
    "therapy: a mini-review and population pharmacokinetic analysis.",
    "J Clin Pharm Ther. 2022;47(8):1091-1102. doi:10.1111/jcpt.13657.",
    "Parameter values: Table 4B. Structure: Methods 2.3-2.5 (Equations 1-4)",
    "and the supplementary aggregate dataset (JCPT-47-1091-s002.xlsx), whose",
    "CLCKRT column is the per-arm CKRT clearance."
  )
  vignette <- "Selig_2022_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    QEFF = list(
      description = paste(
        "Tazobactam clearance by the CKRT circuit (L/h). For the Military",
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
        "for a patient off CKRT. Posted-dataset arms range 0.56-1.4 L/h",
        "(median 1.045 L/h); the literature median tazobactam CKRT",
        "clearance from piperacillin-tazobactam studies is 1.09 L/h, range",
        "0.56-2.73 L/h (Table 2). The additive placement is not printed as",
        "an equation; it is the form under which the posted dataset",
        "reproduces the piperacillin estimates of Table 4A (see the",
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
      notes = "Screened on CL and V by forward addition (Methods 2.6), not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance (study mean or median; Cockcroft-Gault where not reported)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened on CL (Methods 2.6), not retained; the Military Health System patients had the highest mean CrCl (77.34 mL/min) but the lowest tazobactam CL (0.72 L/h) (Results 3.2)."
    ),
    ALB = list(
      description = "Serum albumin (study mean or median)",
      units = "g/L",
      type = "continuous",
      notes = "Screened on V (Methods 2.6), not retained. Reported in g/dL by the source (labelled mg/dL in Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 9L,
    age_range = "arm means/medians 49.8-77.7 years in the posted dataset; Bremmer 47 years (Table 3)",
    weight_range = "arm means/medians 52.5-90 kg (supplementary dataset)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state = paste(
      "Critically ill adults (sepsis, burns, trauma) treated with continuous",
      "kidney replacement therapy (CVVH, CVVHD or CVVHDF), receiving",
      "piperacillin-tazobactam or ceftolozane-tazobactam."
    ),
    dose_range = paste(
      "Tazobactam 250-1000 mg every 6-12 h infused over 0.33-4 h",
      "(supplementary dataset; 1000 mg doses are from the 3 g",
      "ceftolozane-tazobactam case reports)."
    ),
    regions = "Multinational literature cohorts plus the United States Military Health System",
    n_concentrations = 112L,
    notes = paste(
      "Aggregate-data model. 8 literature studies contributed arm-level",
      "curves (Figure 1 legend: Aguilar, Arzuaga, Bremmer, Kohama, Kuti,",
      "Oliver, Valtonen, van der Werf; the posted dataset holds 14 arms from",
      "7 of them, 63 patient-arm entries) and the Military Health System",
      "contributed 3 individual CVVH patients with 10 steady-state",
      "concentrations (Methods 2.2). 112 observations in total (Table 4B)."
    )
  )

  ini({
    lcl <- log(2.49); label("Body (non-CKRT) clearance (L/h)") # Table 4B 'CL (L/hr)' 2.49 (RSE 19.43%, 95% CI 1.55-3.44)
    lvc <- log(30.62); label("Central volume of distribution (L)") # Table 4B 'Vc (L)' 30.62 (RSE 11.53%, 95% CI 23.7-37.54)

    # Between-study (between-arm) variability, exponential (Equation 4).
    # Table 4B prints the variances; sqrt(0.61) = 78% CV as quoted in the
    # Discussion. The 0.58 Pearson correlation in Table 4B is between the
    # post hoc etas; no covariance was reported.
    eta_study_lcl ~ 0.61 # Table 4B 'omega2 CL' 0.61 (RSE 55.91%); eta-shrinkage 1.09%
    eta_study_lvc ~ 0.12 # Table 4B 'omega2 Vc' 0.12 (RSE 30.47%); eta-shrinkage 19.76%

    propSd <- 0.3; label("Proportional residual SD for one patient (fraction)") # Table 4B 'Proportional Error' 0.3 (RSE 12.4%); scaled by 1/sqrt(N) for an arm of N patients (Methods 2.5)
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
