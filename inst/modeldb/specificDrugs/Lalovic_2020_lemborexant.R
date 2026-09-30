Lalovic_2020_lemborexant <- function() {
  description <- paste0(
    "Three-compartment population PK model for oral lemborexant (a dual orexin ",
    "receptor antagonist) in 1892 healthy adults, healthy elderly subjects and ",
    "adults and elderly subjects with insomnia disorder, pooled from 12 phase 1-3 ",
    "studies (Lalovic 2020). Absorption is sequential: a zero-order release of ",
    "duration D1 into the depot after an absorption lag, followed by first-order ",
    "absorption (ka) into the central compartment; elimination is linear from ",
    "central. Apparent clearance (CL/F) decreases with body mass index (power, ",
    "referenced to 25 kg/m^2), with alkaline phosphatase (power, referenced to 71 ",
    "U/L) and is 26% lower in the elderly (> 65 years). The absorption effects -- ",
    "tablet versus capsule and bedtime dosing on D1, tablet and a high-fat meal on ",
    "ka, and a high-fat meal on relative bioavailability -- were estimated on the ",
    "extensively sampled phase 1 data and FIXED in the final model, as were ka, D1, ",
    "the lag time and the IIV on D1, ka and F1. Diagonal IIV on CL/F, Vc/F, Q/F, ",
    "Vp/F, Q2/F, Vp2/F, D1, ka and F1. Combined additive + proportional residual ",
    "error with separate parameter pairs for samples taken up to 3 h after the ",
    "dose and later than 3 h after the dose. Companion landmark exposure-response ",
    "models for six treatment-emergent adverse events are packaged as ",
    "Lalovic_2020_lemborexant_<endpoint>."
  )
  reference <- paste(
    "Lalovic B, Majid O, Aluri J, Landry I, Moline M, Hussein Z. (2020).",
    "Population Pharmacokinetics and Exposure-Response Analyses for the Most",
    "Frequent Adverse Events Following Treatment With Lemborexant, an Orexin",
    "Receptor Antagonist, in Subjects With Insomnia Disorder.",
    "Journal of Clinical Pharmacology 60(12):1642-1654.",
    "doi:10.1002/jcph.1683.",
    sep = " "
  )
  vignette <- "Lalovic_2020_lemborexant"

  # Doses are in mg and volumes in L, so central / vc is mg/L; the 1000 factor in
  # model() reports Cc in ng/mL, the unit of the paper's concentrations and of
  # the additive residual error terms in Table 2.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BMI = list(
      description = "Body mass index at baseline",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (BMI / 25)^-0.428. Lalovic 2020 does not print the",
        "reference BMI. The value 25 kg/m^2 is back-solved from the three",
        "BMI effect sizes quoted in the Results ('Clinically Relevant",
        "Covariates'): relative to the reference, BMI 32 gives 11% and BMI 40",
        "gives 22% higher exposure and BMI 15 a 25% change in CL/F. A 25 kg/m^2",
        "reference reproduces all three (11.1%, 22.3%, 24.5%); the cohort median",
        "26.5 kg/m^2 gives 8.4%, 19.3% and 27.6% and does not. Cohort median 26.5,",
        "range 14.4-62.1 kg/m^2 (Table 1)."
      ),
      source_name = "BMI"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase at baseline",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (ALP / 71)^-0.118. The reference is not printed in",
        "Table 2; 71 U/L is the cohort median (Table 1) and is the 'population",
        "median alkaline phosphatase of 71 IU/L' Lalovic 2020 assumed for the",
        "Table 3 simulations. It reproduces the Results statement that ALP 150",
        "U/L raises and ALP 35 U/L lowers exposure by about 9%. Range 13-256 U/L",
        "(Table 1). Baseline value, time-fixed."
      ),
      source_name = "ALP"
    ),
    AGE_GT65 = list(
      description = "Elderly indicator, age > 65 years",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adults, 18-65 years)",
      notes = paste(
        "Multiplicative effect on CL/F, 0.739^AGE_GT65 (26% lower CL/F in the",
        "elderly). The Methods covariate definition, the Results text and the",
        "Table 3 footnote categories ('Elderly (> 65 years old)' versus 'Adult",
        "(18-65 years old)') put the cut point at > 65 years, which is the",
        "encoding used here. Table 1 and the Discussion instead write 'Elderly",
        "(>= 65 years)'; the two definitions differ only for subjects aged exactly",
        "65. Table 1 counts 547 elderly and 1345 adults."
      ),
      source_name = "Elderly"
    ),
    FORM_TABLET = list(
      description = "Tablet formulation indicator (per dose record)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the capsule used in the early phase 1 studies)",
      notes = paste(
        "NOTE the comparator is a CAPSULE, not the non-tablet oral liquid named as",
        "this canonical's default reference category (as for Zhu 2018 asunaprevir",
        "and Wada 2023 sparsentan). Multiplies D1 by 0.254 (0.467 h capsule to",
        "0.119 h tablet) and ka by 1.12; both effects FIXED in the final model.",
        "The tablet was used in phases 2 and 3. Table 1: tablet 1755, capsule 137",
        "subjects."
      ),
      source_name = "Tablet"
    ),
    FED_HIGHFAT = list(
      description = "Dose administered with a standard high-fat meal (per dose record)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Lalovic 2020 Methods describe the food covariate as 'coadministration",
        "with food [standard high-fat meal]'. Multiplies ka by 0.695 (30%",
        "slower) and relative bioavailability F1 by 1.21 (21% higher); both",
        "effects FIXED in the final model."
      ),
      source_name = "Food"
    ),
    DOSETIME_EVENING = list(
      description = "Bedtime (nighttime) dose indicator (per dose record)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (daytime dosing, as in the extensively sampled phase 1 studies)",
      notes = paste(
        "Lalovic 2020 calls this 'bedtime dosing' in the text and 'D1, nighttime",
        "dosing' in Table 2. It is the dose-record property of being taken at",
        "bedtime, which is how lemborexant is dosed in the phase 2 and 3 insomnia",
        "trials; the paper gives no clock-time window. Multiplies D1 by 2.33;",
        "FIXED in the final model."
      ),
      source_name = "Bedtime dosing"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_PPI = list(
      description = "Concomitant proton pump inhibitor use",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Significant on CL/F in forward addition but removed at backward",
        "elimination (P < .001); not in the final model (Lalovic 2020 Results,",
        "Final Model). 112 of 1892 subjects."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened; did not enter the final model (Results; Supplemental Figure S4)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; did not enter the final model -- BMI was retained instead (Results, Discussion)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened; did not enter the final model (Discussion)."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "lemborexant",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "lemborexant",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "lemborexant",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "lemborexant",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1892,
    n_studies = 12,
    n_observations = 12230,
    age_range = "18-88 years",
    age_median = "57 years",
    weight_range = "37-168 kg",
    weight_median = "74.1 kg",
    bmi_range = "14.4-62.1 kg/m^2 (median 26.5)",
    sex_female_pct = 66,
    race_ethnicity = c(
      White = 70.5,
      `Black or African American` = 17.7,
      Japanese = 8.2,
      `Other Asian` = 1.7,
      Chinese = 0.3,
      `American Indian / Alaskan / other / missing` = 1.6
    ),
    disease_state = paste(
      "Healthy adults and healthy elderly subjects (phase 1) and adults and",
      "elderly subjects with insomnia disorder (phases 2 and 3; SUNRISE 1 and",
      "SUNRISE 2), plus a phase 2 study in irregular sleep-wake rhythm disorder"
    ),
    dose_range = "1-100 mg once daily orally (phase 3: 5 and 10 mg at bedtime)",
    hepatic_function = "ALP median 71 U/L (range 13-256); ALT median 17 U/L; AST median 19 U/L",
    renal_function = "Creatinine clearance median 97.3 mL/min (range 26.8-319)",
    co_medication = "Concomitant PPI 112 subjects; weak CYP3A inhibitors 22 subjects",
    notes = paste(
      "Demographics from Lalovic 2020 Table 1 (n = 1892). 6 extensively and 3",
      "sparsely sampled phase 1 studies (407 subjects, 6543 observations), 1",
      "phase 2 study and 2 phase 3 studies (1485 subjects, 5687 observations);",
      "Supplemental Table S1. Race percentages are computed from the Table 1",
      "counts. Elderly (> 65 years) 547 subjects."
    )
  )

  ini({
    # --- Disposition (Lalovic 2020 Table 2, final model) ------------------------
    # Reference subject: BMI 25 kg/m^2, ALP 71 U/L, adult (<= 65 years), capsule,
    # fasted, daytime dosing. All disposition parameters are apparent (/F).
    lcl <- log(22.7); label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/h)' = 22.7
    lvc <- log(9.09); label("Apparent central volume V2/F (L)") # Table 2 'V2/F (L)' = 9.09
    lq <- log(32.1); label("Apparent intercompartmental clearance to the first peripheral compartment Q3/F (L/h)") # Table 2 'Q3/F (L/h)' = 32.1
    lvp <- log(278); label("Apparent first peripheral volume V3/F (L)") # Table 2 'V3/F (L)' = 278
    lq2 <- log(31.0); label("Apparent intercompartmental clearance to the second peripheral compartment Q4/F (L/h)") # Table 2 'Q4/F (L/h)' = 31.0
    lvp2 <- log(783); label("Apparent second peripheral volume V4/F (L)") # Table 2 'V4/F (L)' = 783

    # --- Absorption (Table 2; all FIXED in the final model) ---------------------
    lka <- fixed(log(0.532)); label("First-order absorption rate constant ka, capsule (1/h)") # Table 2 'Ka, capsule (h-1)' = 0.532, Fixed
    ld1 <- fixed(log(0.467)); label("Duration of the zero-order release into the depot D1, capsule (h)") # Table 2 'D1, capsule (h)' = 0.467, Fixed
    ltlag <- fixed(log(0.403)); label("Absorption lag time ALAG1 (h)") # Table 2 'ALAG1 (h)' = 0.403, Fixed
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1, fasted reference (unitless)") # Table 2 lists only the food effect on F1, so F1 = 1 in the fasted state (structural anchor)

    # --- Covariate effects on CL/F (Table 2; power and categorical forms per the Table 2 footnote) ---
    e_bmi_cl <- -0.428; label("Power exponent of BMI on CL/F, referenced to 25 kg/m^2 (unitless)") # Table 2 'CL/F, BMI' = -0.428
    e_alp_cl <- -0.118; label("Power exponent of ALP on CL/F, referenced to 71 U/L (unitless)") # Table 2 'CL/F,ALP' = -0.118
    e_age_gt65_cl <- 0.739; label("Multiplicative factor on CL/F for the elderly, age > 65 years (unitless)") # Table 2 'CL/F, elderly' = 0.739

    # --- Covariate effects on absorption (Table 2; all FIXED) ------------------
    e_form_tablet_d1 <- fixed(0.254); label("Multiplicative factor on D1 for the tablet versus the capsule (unitless)") # Table 2 'D1, tablet' = 0.254, Fixed; 0.467 * 0.254 = 0.119 h, the '0.118 hours' of the Results
    e_dosetime_evening_d1 <- fixed(2.33); label("Multiplicative factor on D1 for bedtime (nighttime) dosing (unitless)") # Table 2 'D1, nighttime dosing' = 2.33, Fixed
    e_form_tablet_ka <- fixed(1.12); label("Multiplicative factor on ka for the tablet versus the capsule (unitless)") # Table 2 'Ka, tablet' = 1.12, Fixed
    e_fed_highfat_ka <- fixed(0.695); label("Multiplicative factor on ka for dosing with a high-fat meal (unitless)") # Table 2 'Ka, food' = 0.695, Fixed
    e_fed_highfat_fdepot <- fixed(1.21); label("Multiplicative factor on relative bioavailability F1 for dosing with a high-fat meal (unitless)") # Table 2 'F1, food' = 1.21, Fixed

    # --- Inter-individual variability (Table 2, diagonal) ----------------------
    # Table 2 reports IIV as %CV and its abbreviation list defines 'CV, square
    # root of variance x 100', so each variance is (CV/100)^2.
    etalcl ~ 0.231361 # Table 2 IIV 'CL/F' 48.1 %CV -> 0.481^2
    etalvc ~ 2.0164 # Table 2 IIV 'V2/F' 142 %CV -> 1.42^2
    etalq ~ 0.322624 # Table 2 IIV 'Q3/F' 56.8 %CV -> 0.568^2
    etalvp ~ 0.6724 # Table 2 IIV 'V3/F' 82.0 %CV -> 0.820^2
    etalq2 ~ 0.216225 # Table 2 IIV 'Q4/F' 46.5 %CV -> 0.465^2
    etalvp2 ~ 0.171396 # Table 2 IIV 'V4/F' 41.4 %CV -> 0.414^2
    etald1 ~ fixed(2.7889) # Table 2 IIV 'D1' 167 %CV -> 1.67^2
    etalka ~ fixed(0.191844) # Table 2 IIV 'Ka' 43.8 %CV -> 0.438^2
    etalfdepot ~ fixed(0.463761) # Table 2 IIV 'F1' 68.1 %CV -> 0.681^2

    # --- Residual error (Table 2), switched at 3 h after the dose --------------
    propSd_early <- 0.329; label("Proportional residual SD for samples <= 3 h after the dose (fraction)") # Table 2 'Proportional (TAD <= 3 h),% CV' = 32.9
    addSd_early <- 2.62; label("Additive residual SD for samples <= 3 h after the dose (ng/mL)") # Table 2 'Additive (TAD <= 3 h), ng/mL' = 2.62
    propSd_late <- 0.143; label("Proportional residual SD for samples > 3 h after the dose (fraction)") # Table 2 'Proportional (TAD > 3 h),% CV' = 14.3
    addSd_late <- 0.0189; label("Additive residual SD for samples > 3 h after the dose (ng/mL)") # Table 2 'Additive (TAD > 3 h), ng/mL' = 0.0189
  })

  model({
    # 1. Individual parameters ----------------------------------------------------
    # Continuous covariates as power terms and dichotomous covariates as
    # theta^covariate (Table 2 footnote).
    cl <- exp(lcl + etalcl) *
      (BMI / 25)^e_bmi_cl *
      (ALP / 71)^e_alp_cl *
      e_age_gt65_cl^AGE_GT65
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    ka <- exp(lka + etalka) *
      e_form_tablet_ka^FORM_TABLET *
      e_fed_highfat_ka^FED_HIGHFAT
    d1 <- exp(ld1 + etald1) *
      e_form_tablet_d1^FORM_TABLET *
      e_dosetime_evening_d1^DOSETIME_EVENING
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot) * e_fed_highfat_fdepot^FED_HIGHFAT

    # 2. Micro-constants ------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system (Supplemental Figure S2: depot CMT 1, central CMT 2,
    #    peripheral CMT 3 and CMT 4) -------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Absorption input -------------------------------------------------------------
    # Mixed zero- then first-order absorption (NONMEM ADVAN12 with D1 on the
    # depot): after the lag the dose enters the depot at a constant rate over D1
    # hours and leaves it first-order at ka. Dose records MUST carry rate = -2 for
    # rxode2 to honour dur(depot); without it the dose is a bolus and D1 is
    # ignored.
    alag(depot) <- tlag
    dur(depot) <- d1
    f(depot) <- fdepot

    # 5. Observation and error -----------------------------------------------------
    Cc <- 1000 * central / vc

    # Residual error switches at 3 h after the most recent dose (Results, Base
    # Model). tad() is hoisted to its own line.
    tad_obs <- tad()
    early <- (tad_obs <= 3)
    propSd <- propSd_early * early + propSd_late * (1 - early)
    addSd <- addSd_early * early + addSd_late * (1 - early)
    Cc ~ add(addSd) + prop(propSd)
  })
}
