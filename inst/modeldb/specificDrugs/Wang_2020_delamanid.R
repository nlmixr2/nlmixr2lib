Wang_2020_delamanid <- function() {
  description <- "Two-compartment population PK model for oral delamanid in adults with pulmonary multidrug-resistant tuberculosis (Wang 2020): first-order absorption with lag time, morning doses (dosed into depot) and evening doses (dosed into depot2) with separate absorption rate constants, lag times and relative bioavailability, and relative bioavailability that also falls with dose and depends on inpatient versus outpatient setting and enrollment region."
  reference <- "Wang X, Mallikaarjun S, Gibiansky E. Population Pharmacokinetic Analysis of Delamanid in Patients with Pulmonary Multidrug-Resistant Tuberculosis. Antimicrob Agents Chemother. 2020;65(1):e01202-20. doi:10.1128/AAC.01202-20"
  vignette <- "Wang_2020_delamanid"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect with a single shared exponent on both the apparent central (V2/F) and apparent peripheral (V3/F) volumes (Table 4 theta13, row 'V2, V3,WT'). Reference weight 55 kg: the population median (Table 2) and the weight of the typical patient in Figure 1 panels A, C and D; the Figure 1 panel C/D ratios at 40, 75 and 90 kg reproduce (WT/55)^0.316. Apparent clearance is independent of weight.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative effect on V3/F: female V3/F is 1.65-fold the male value (Table 4 theta14, row 'V3,SEX'; Results: '65% ... higher than in male patients').",
      source_name = "SEX"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline (time-fixed) value. Canonical units are g/L; the model converts inline to g/dL (alb_gdl = ALB / 10) because the published effect is calibrated in g/dL. Hockey-stick power effect on CL/F (Table 4 footnote b): CL/F ~ IALB^theta17 with IALB = 1 when albumin >= 3.4 g/dL and IALB = ALB/3.4 when albumin < 3.4 g/dL, theta17 = -0.892, so clearance rises only in hypoalbuminemic patients.",
      source_name = "ALB"
    ),
    CONMED_EFV = list(
      description = "Concomitant efavirenz indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no efavirenz; HIV-negative patients and HIV-positive patients on other antiretroviral therapy)",
      notes = "Multiplicative effect on CL/F: 1.35-fold with efavirenz (Table 4 theta18, row 'CL efavirenz'). Tested in stage 2 among the 48 HIV-positive patients of the trial 213 HIV subtrial; 22 patients (3.0%) received efavirenz (Table 2). HIV status itself and lamivudine or tenofovir had no effect.",
      source_name = "efavirenz"
    ),
    DOSE_DELAMANID_MG = list(
      description = "Delamanid dose amount on this dose record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Used as a level selector for the dose-dependent relative bioavailability, not as a continuous dose-response driver. Reference dose 100 mg (F1 = 1). Table 4 theta7 (F1 for the 200-mg dose, 0.760) applies to doses above 100 mg up to 200 mg and theta8 (F1 for doses above 200 mg, 0.580) applies to the 250- and 300-mg doses. Only 100, 200, 250 and 300 mg doses were studied; the model does not interpolate between them.",
      source_name = "AMT"
    ),
    OUTPATIENT = list(
      description = "Indicator that the dose was taken in an outpatient setting",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (inpatient / hospitalized)",
      notes = "Time-varying at the dose record; patients were hospitalized all, part or none of the time depending on the trial (Discussion). Multiplies relative bioavailability by 1.09 (Table 4 theta12, row 'F1,out').",
      source_name = "OUT"
    ),
    REGION_EASTASIA = list(
      description = "Northeast Asian region indicator (patients enrolled in China, Japan or Korea)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian region, including the 2 Asian-race patients enrolled in Peru)",
      notes = "Multiplies relative bioavailability by 1.53 (Table 4 theta15, row 'F1,NEAsian'). Table 2 defines the group as Asian patients from Northeast Asia (China, Japan, or Korea), 99 patients (13.3%). Mutually exclusive with REGION_SOUTHEASTASIA.",
      source_name = "NE Asian"
    ),
    REGION_SOUTHEASTASIA = list(
      description = "Southeast Asian region indicator (patients enrolled in the Philippines)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian region)",
      notes = "Multiplies relative bioavailability by 1.40 (Table 4 theta16, row 'F1,SEAsian'). Table 2 defines the group as Asian patients from Southeast Asia (Philippines), 200 patients (26.9%). Mutually exclusive with REGION_EASTASIA.",
      source_name = "SE Asian"
    ),
    STUDY_242_07_208 = list(
      description = "Indicator that the observation record comes from trial 242-07-208 (NCT02573350)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (trials 242-07-204, 242-08-210 and 242-09-213)",
      notes = "Selects the larger proportional residual error of trial 208 (Table 4 SIGMA(2,2) = 0.174, CV 41.8%, vs SIGMA(1,1) = 0.0715, CV 26.7%, in the other trials). Trial 208 was the open-label extension that enrolled patients completing trial 204, with sparse sampling every 4 weeks (Table 1).",
      source_name = "STUDY"
    ),
    STUDY_242_09_213 = list(
      description = "Indicator that the observation record comes from phase III trial 242-09-213 (NCT01424670)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (trials 242-07-204, 242-07-208 and 242-08-210)",
      notes = "Selects the larger additive residual error of trial 213 (Table 4 SIGMA(3,3) = 1950 (ng/mL)^2, SD 44.2 ng/mL, vs SIGMA(4,4) = 2.39 (ng/mL)^2, SD 1.55 ng/mL, in the other trials). Trial 213 was the phase III trial (100 mg BID for 8 weeks then 200 mg QD for 18 weeks), added to the analysis in stage 2.",
      source_name = "STUDY"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "delamanid", units = "mg", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "delamanid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "delamanid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "delamanid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 744,
    n_studies = 4,
    n_observations = 20483,
    age_range = "18-64 years",
    age_median = "33 years",
    weight_range = "27-99.6 kg",
    weight_median = "55 kg",
    sex_female_pct = 30.5,
    race_ethnicity = c(Asian = 40.5, White = 23.8, Black = 6.9, Other = 28.9),
    disease_state = "Pulmonary multidrug-resistant tuberculosis (70.4% MDR-TB, 16.8% pre-XDR-TB, 12.8% XDR-TB); 4.2% HIV-coinfected",
    dose_range = "100 mg BID, 200 mg BID, 250 mg BID, 300 mg BID and 200 mg QD (morning) orally with food, for 8 to 28 weeks",
    regions = "Asia (Philippines, China, Japan, Korea), Peru and other global sites",
    renal_function = "15.3% CKD stage II and 1.6% CKD stage III by MDRD",
    hypoalbuminemia = "27.2% baseline albumin < 3.4 g/dL; 8.7% < 2.8 g/dL",
    co_medication = "Optimized background regimen in all but 10 patients (Table 3); efavirenz 3.0%, lamivudine 3.1%, tenofovir 2.2%",
    notes = "Pooled phase II trials 242-07-204, 242-07-208 and 242-08-210 and phase III trial 242-09-213 (Table 1). Baseline demographics from Table 2. The model was built in two stages: stage 1 on trials 204, 208 and 210; stage 2 added trial 213."
  )

  ini({
    # Structural parameters: Table 4 (final model). Reference: 100-mg
    # morning dose, inpatient, non-Asian region, male, 55 kg, albumin
    # >= 3.4 g/dL, no efavirenz.
    lcl <- log(37.1); label("Apparent clearance CL/F (L/h)") # Table 4 theta1 = 37.1 L/h
    lvc <- log(655); label("Apparent central volume V2/F at 55 kg (L)") # Table 4 theta2 = 655 L
    lq <- log(104); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 4 theta3 = 104 L/h
    lvp <- log(870); label("Apparent peripheral volume V3/F for a 55-kg male (L)") # Table 4 theta4 = 870 L
    lka_am <- log(0.397); label("Absorption rate constant after a morning dose (1/h)") # Table 4 theta5 = 0.397 1/h
    ltlag_am <- log(0.825); label("Absorption lag time after a morning dose (h)") # Table 4 theta6 = 0.825 h
    lka_pm <- log(0.248); label("Absorption rate constant after an evening dose (1/h)") # Table 4 theta9 = 0.248 1/h
    ltlag_pm <- log(1.38); label("Absorption lag time after an evening dose (h)") # Table 4 theta10 = 1.38 h
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 of the reference 100-mg morning inpatient dose in a non-Asian patient (unitless)") # Results: 100-mg morning inpatient non-Asian dose is the F1 reference; Figure 1B 'Typical F1 = 1'

    # Relative-bioavailability covariate effects (multiplicative ratios on F1)
    e_dose200_f <- 0.760; label("Relative bioavailability of a 200-mg dose vs a 100-mg dose (ratio)") # Table 4 theta7 = 0.760
    e_dosegt200_f <- 0.580; label("Relative bioavailability of a >200-mg (250 or 300 mg) dose vs a 100-mg dose (ratio)") # Table 4 theta8 = 0.580
    e_dosetime_evening_f <- 1.26; label("Relative bioavailability of an evening dose (depot2) vs a morning dose (depot) (ratio)") # Table 4 theta11 = 1.26
    e_outpatient_f <- 1.09; label("Relative bioavailability in an outpatient vs an inpatient setting (ratio)") # Table 4 theta12 = 1.09
    e_region_eastasia_f <- 1.53; label("Relative bioavailability in Northeast Asian vs non-Asian patients (ratio)") # Table 4 theta15 = 1.53
    e_region_southeastasia_f <- 1.40; label("Relative bioavailability in Southeast Asian vs non-Asian patients (ratio)") # Table 4 theta16 = 1.40

    # Volume and clearance covariate effects
    e_wt_vc_vp <- 0.316; label("Power exponent of (WT/55) shared by V2/F and V3/F (unitless)") # Table 4 theta13 = 0.316
    e_sexf_vp <- 1.65; label("V3/F in female vs male patients (ratio)") # Table 4 theta14 = 1.65
    e_alb_cl <- -0.892; label("Power exponent of IALB on CL/F, IALB = min(ALB/3.4 g/dL, 1) (unitless)") # Table 4 theta17 = -0.892 and footnote b
    e_conmed_efv_cl <- 1.35; label("CL/F with vs without concomitant efavirenz (ratio)") # Table 4 theta18 = 1.35

    # IIV: Table 4, diagonal OMEGA, log-normal (Methods: Pi = P exp(eta))
    etalcl ~ 0.056 # Table 4 OMEGA(1,1) = 0.056 (CV 23.7%)
    etalvp ~ 0.152 # Table 4 OMEGA(2,2) = 0.152 (CV 39.0%)
    etalka_am ~ 0.517 # Table 4 OMEGA(3,3) = 0.517 (CV 71.9%)
    etalfdepot ~ 0.0344 # Table 4 OMEGA(4,4) = 0.0344 (CV 18.5%)
    etalq ~ 0.456 # Table 4 OMEGA(5,5) = 0.456 (CV 67.5%)
    etalka_pm ~ 0.343 # Table 4 OMEGA(6,6) = 0.343 (CV 58.6%)

    # Residual error: combined proportional + additive with independent
    # epsilons (Methods: C = Chat(1 + w_pr eps1) + w_add eps2). Table 4
    # reports variances; the SDs below are their square roots.
    propSd <- 0.2674; label("Proportional residual SD, trials 204, 210 and 213 (fraction)") # Table 4 SIGMA(1,1) = 0.0715 -> sqrt = 0.2674 (CV 26.7%)
    propSd_t208 <- 0.4171; label("Proportional residual SD, trial 208 (fraction)") # Table 4 SIGMA(2,2) = 0.174 -> sqrt = 0.4171 (CV 41.8%)
    addSd <- 1.546; label("Additive residual SD, trials 204, 208 and 210 (ng/mL)") # Table 4 SIGMA(4,4) = 2.39 (ng/mL)^2 -> sqrt = 1.546 (SD 1.55)
    addSd_t213 <- 44.16; label("Additive residual SD, trial 213 (ng/mL)") # Table 4 SIGMA(3,3) = 1950 (ng/mL)^2 -> sqrt = 44.16 (SD 44.2)
  })

  model({
    # Covariate terms. Albumin is carried in g/L (canonical) and converted
    # to the g/dL scale on which theta17 was estimated. IALB is capped at 1,
    # so only hypoalbuminemia (< 3.4 g/dL) changes CL/F (Table 4 footnote b).
    alb_gdl <- ALB / 10
    ialb <- 1
    if (alb_gdl < 3.4) {
      ialb <- alb_gdl / 3.4
    }

    # Dose-level relative bioavailability (Table 4 theta7, theta8): the
    # 100-mg dose is the reference; 200 mg and doses above 200 mg each have
    # their own estimate.
    fdose <- 1
    if (DOSE_DELAMANID_MG > 200) {
      fdose <- e_dosegt200_f
    } else if (DOSE_DELAMANID_MG > 100) {
      fdose <- e_dose200_f
    }

    # Individual PK parameters
    cl <- exp(lcl + etalcl) * ialb^e_alb_cl * e_conmed_efv_cl^CONMED_EFV
    vc <- exp(lvc) * (WT / 55)^e_wt_vc_vp
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (WT / 55)^e_wt_vc_vp * e_sexf_vp^SEXF

    # Morning and evening doses are absorbed from separate depots, each with
    # its own absorption rate constant (with its own IIV), lag time and
    # relative bioavailability (Table 4 theta5/theta6 vs theta9/theta10/
    # theta11). Morning doses go into `depot` and evening doses into
    # `depot2`.
    ka <- exp(lka_am + etalka_am)
    ka2 <- exp(lka_pm + etalka_pm)
    tlag <- exp(ltlag_am)
    tlag2 <- exp(ltlag_pm)

    fdepot <- exp(lfdepot + etalfdepot) * fdose *
      e_outpatient_f^OUTPATIENT *
      e_region_eastasia_f^REGION_EASTASIA *
      e_region_southeastasia_f^REGION_SOUTHEASTASIA
    fdepot2 <- fdepot * e_dosetime_evening_f

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag
    f(depot2) <- fdepot2
    alag(depot2) <- tlag2

    # Dose in mg and volume in L give mg/L; x 1000 reports ng/mL.
    Cc <- 1000 * central / vc

    # Trial-specific residual error (Table 4): trial 208 has its own
    # proportional SD and trial 213 its own additive SD.
    propSdEff <- propSd * (1 - STUDY_242_07_208) + propSd_t208 * STUDY_242_07_208
    addSdEff <- addSd * (1 - STUDY_242_09_213) + addSd_t213 * STUDY_242_09_213
    Cc ~ add(addSdEff) + prop(propSdEff)
  })
}
