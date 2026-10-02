He_2022_paclitaxel <- function() {
  description <- "Two-compartment population PK model for oral paclitaxel co-administered with the P-glycoprotein inhibitor encequidar (Oraxol) in adults with advanced or metastatic solid tumors (He 2022). First-order absorption with a fixed lag time, linear elimination, bioavailability fixed at 0.119 for the capsule with a proportional increase for the oral solution, a proportional increase in the central volume for non-Asian patients, and a log-additive residual error. Pooled data from 197 patients in seven studies."
  reference <- "He J, Jackson CGCA, Deva S, Hung T, Clarke K, Segelov E, Chao TY, Dai MS, Yeh HT, Ma WW, Kramer D, Chan WK, Kwan R, Cutler D, Zhi J. Population pharmacokinetics for oral paclitaxel in patients with advanced/metastatic solid tumors. CPT Pharmacometrics Syst Pharmacol. 2022;11(7):867-879. doi:10.1002/psp4.12799"
  vignette <- "He_2022_paclitaxel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot = list(analyte = "paclitaxel", units = "mg", specimen = "administration site", verified = FALSE),
    central = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "paclitaxel", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = any other race).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (Asian; the typical-value central-volume reference in He 2022 Table 2).",
      notes = "The control stream (He 2022 Table S6) codes the indicator the other way round: RACB = 1 when RACE is 1-4 (Caucasian, African American, American Indian, Native Hawaiian or Other Pacific Islander) and 0 when RACE = 0 (Asian), and applies V2 = THETA(2) * (1 + THETA(8) * RACB). The model therefore uses (1 - RACE_ASIAN) in place of RACB, which leaves the coefficient 0.696 and its sign unchanged: non-Asian V2 = 50.7 * 1.696 = 86.0 L versus 50.7 L for Asian patients (He 2022 Results). The RACE = 0 coding of Asian was confirmed against the posted analysis dataset (Appendix S2): RACE = 0 gives 110 patients, matching the 110 Asian patients of Table 1. Race was also tested on F1 and on CL (Table S3); only the effect on V2 was retained.",
      source_name = "RACE"
    ),
    FORM_SOLUTION = list(
      description = "Oral-solution formulation indicator (1 = paclitaxel intravenous solution given orally, 0 = oral capsule).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (oral capsule; F1 = 0.119 fixed).",
      notes = "The source column FORM is 1 = solution and 2 = capsule (He 2022 Figure 3 legend); the control stream applies F1 = THETA(6) * (1 + THETA(9) * (2 - FORM)), so FORM_SOLUTION = 2 - FORM. The oral solution was used only in the first study (15 patients); all later studies and the commercial product use capsules or tablets (He 2022 Results). Per-dose-record covariate.",
      source_name = "FORM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = "Tested on CL, V2 and F1 (linear, exponential and power forms; He 2022 Table S3) and not retained in the final model."
    ),
    BSA = list(
      description = "Body surface area at baseline.",
      units = "m^2",
      type = "continuous",
      notes = "Tested on V2 and F1 (He 2022 Table S3) and not retained. Dosing in the source studies is BSA-based (mg/m^2), so BSA still sets the dose amount in a simulation; it has no effect on the PK parameters."
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units = "g/L",
      type = "continuous",
      notes = "Tested on V2 (He 2022 Table S3) and not retained."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL and V2 (He 2022 Table S3) and not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 197L,
    n_studies = 7L,
    n_observations = 4322L,
    age_range = "32-81 years (mean 59.6, SD 10.7)",
    weight_range = "38-139 kg (mean 67.2, SD 15.9)",
    bsa_range = "1.29-2.46 m^2 (mean 1.73, SD 0.23)",
    sex_female_pct = 52.3,
    race_ethnicity = c(
      Asian = 55.8,
      Caucasian = 41.1,
      `African American` = 2.0,
      `American Indian` = 0.5,
      `Native Hawaiian or Other Pacific Islander` = 0.5
    ),
    disease_state = "Adults with advanced or metastatic solid tumors (gastric cancer, metastatic breast cancer, lung cancer and others).",
    dose_range = "Oral paclitaxel 60-420 mg/m^2 (as the intravenous solution given orally or as liquid-filled capsules) in two early studies; 270 mg flat to 313 mg/m^2 once daily for 2-5 days in two dose-finding studies; 205 mg/m^2 once daily on days 1-3 weekly in two bioavailability studies. Every dose was given with the P-glycoprotein inhibitor encequidar (15 mg tablet in the later studies).",
    regions = "Republic of Korea, United States, New Zealand, Australia, Taiwan",
    renal_function = "Normal 40.1%, mild impairment 43.1%, moderate impairment 16.8% (Cockcroft-Gault creatinine clearance 33-230 mL/min)",
    hepatic_function = "Normal 89.8%, impaired 10.2% (19 mild, 1 moderate by NCI-ODWG criteria)",
    formulation = c(`Oral solution` = 7.6, `Oral capsule` = 92.4),
    notes = "Demographics from He 2022 Table 1. Studies HM-OXL-101, HM-OXL-201, ORAX-01-13-US, ORAX-01-14-NZ, KX-ORAX-002, KX-ORAX-003 and KX-ORAX-007 (He 2022 Table S1). 97% of patients were dosed fasted."
  )

  ini({
    # Structural parameters -- He 2022 Table 2 / Table S4 final estimates.
    # F1 was fixed at 0.119 for the capsule, so CL and the volumes are
    # systemic (not apparent) values given that F1 (Table 2 footnote a).
    lcl <- log(33.7); label("Clearance CL (L/h)") # Table 2 CL = 33.7 L/h (RSE 4.2%)
    lvc <- log(50.7); label("Central volume V2 in Asian patients (L)") # Table 2 V2 = 50.7 L (RSE 15.9%)
    lq <- log(40.6); label("Intercompartmental clearance Q (L/h)") # Table 2 Q = 40.6 L/h (RSE 5.1%)
    lvp <- log(855); label("Peripheral volume V3 (L)") # Table 2 V3 = 855 L (RSE 4.9%)
    lka <- log(0.724); label("First-order absorption rate constant KA (1/h)") # Table 2 KA = 0.724 1/h (RSE 5.2%)
    ltlag <- fixed(log(0.215)); label("Absorption lag time ALAG1 (h)") # Table 2 ALAG1 = 0.215 h, fixed; Table S6 THETA(7) 0.215 FIX
    lfdepot <- fixed(log(0.119)); label("Oral bioavailability F1 of the capsule (fraction)") # Table 2 F1 = 0.119, fixed; Table S6 THETA(6) 0.119 FIX

    # Covariate effects -- linear-proportional, Table S6 $PK.
    e_race_nonasian_vc <- 0.696; label("Proportional increase in V2 for non-Asian vs Asian patients (unitless)") # Table 2 'Race on V2 (proportional)' = 0.696 (RSE 36.6%)
    e_form_solution_fdepot <- 0.895; label("Proportional increase in F1 for the oral solution vs the capsule (unitless)") # Table 2 'Formulation on F1 (proportional)' = 0.895 (RSE 28.5%)

    # IIV -- exponential, diagonal (Table S6 $OMEGA; etas on Q, V3, KA and
    # ALAG1 are 0 FIX). Table 2 and Table S4 report CV%, converted here as
    # omega^2 = log(CV^2 + 1). The third IIV row is on F1: Table 2 labels it
    # 'ETA Q' but Table S4 ('BSV on F1') and the Table S6 control stream
    # (eta6 on F1, eta3 on Q fixed to 0) both place it on F1.
    etalcl ~ 0.1118 # Table 2 'ETA CL' = 34.4 CV% -> log(0.344^2 + 1)
    etalvc ~ 1.4104 # Table 2 'ETA V2' = 176 CV% -> log(1.76^2 + 1)
    etalfdepot ~ 0.2073 # Table S4 'BSV on F1' = 48.0 CV% -> log(0.480^2 + 1)

    # Residual error -- NONMEM Y = LOG(F) + ERR(1) (Table S6 $ERROR), i.e. a
    # log-normal residual. Table 2 'Log additive' = 0.208 is the SIGMA
    # variance, so the log-scale SD is sqrt(0.208).
    expSd <- 0.4561; label("Log-additive residual error SD (log scale)") # Table 2 'Log additive' = 0.208 (variance) -> sqrt(0.208)
  })

  model({
    # Individual parameters (Table S6 $PK). RACB in the control stream is the
    # non-Asian indicator, i.e. 1 - RACE_ASIAN.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (1 + e_race_nonasian_vc * (1 - RACE_ASIAN))
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_form_solution_fdepot * FORM_SOLUTION)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with first-order absorption (NONMEM ADVAN4 TRANS4).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL (Table S6 S2 = V2/1000).
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
