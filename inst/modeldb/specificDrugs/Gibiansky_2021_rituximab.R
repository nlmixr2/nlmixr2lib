Gibiansky_2021_rituximab <- function() {
  description <- "Two-compartment population PK model of intravenous and subcutaneous rituximab in adults with chronic lymphocytic leukemia (CLL); first-order subcutaneous absorption with bioavailability, and total clearance equal to a time-independent component cl_exp_inf plus a mono-exponentially decaying target-mediated component cl_exp_component * exp(-cl_exp_kdes * time). Covariates are body surface area on all clearances and volumes, baseline white blood cell count and baseline tumor size on the time-dependent clearance, body mass index on the absorption rate and bioavailability, and sex on the central volume; log-scale residual error whose SD falls with concentration, with study multipliers and IIV on the residual SD (Gibiansky 2021)."
  reference <- "Gibiansky E, Gibiansky L, Chavanne C, Frey N, Jamois C. Population pharmacokinetic and exposure-response analyses of intravenous and subcutaneous rituximab in patients with chronic lymphocytic leukemia. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):914-927. doi:10.1002/psp4.12665"
  vignette <- "Gibiansky_2021_rituximab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # The IIV on the residual-error SD (NONMEM ETA(7)) has no structural
  # typical-value partner; the residual SD is the concentration-dependent
  # expression built from sdL / sdH / sd50 in model().
  paper_specific_etas <- c("etaRUV")

  compartmentData <- list(
    depot = list(analyte = "rituximab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rituximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "rituximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 1.9 m^2 (Supplementary Table S5: IBSA = BSA/1.9). Power exponent 1.37 shared by CL_T, CL_inf and Q, and 0.8 shared by VC and VP (Table 3 theta14 and theta15). The BSA formula is not stated in the paper. The NONMEM code sets IBSA = 1 for the <1% of subjects with missing BSA (FLGBSA = 1); the library model has no missing-value flag, so supply an observed or imputed BSA for every subject.",
      source_name = "BSA"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 27 kg/m^2 (Supplementary Table S5: IBMI = BMI/27). Power exponent -1.01 on the subcutaneous absorption rate ka and -0.465 on the subcutaneous bioavailability F_SC (Table 3 theta18 and theta19). The NONMEM code sets IBMI = 1 when BMI is missing (FLGBMI = 1).",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Source column SEXF uses the canonical encoding directly (Supplementary Table S5 $INPUT). Multiplicative effect theta17^SEXF on VC with theta17 = 0.909 (VC 9.1% lower in women; Table 3 and Table 4).",
      source_name = "SEXF"
    ),
    WBC = list(
      description = "Baseline white blood cell count",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value (Table 2 footnote b). Reference 100 x 10^9/L (Supplementary Table S5: IWBC = WBC/100). Power exponent 0.223 on the time-dependent clearance CL_T (Table 3 theta16). In CLL the count is dominated by circulating leukaemic B cells, so it acts as a marker of CD20 target burden. Missing values (15%) were imputed to the mean and flagged (FLGWBC = 1 sets IWBC = 1).",
      source_name = "WBC"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (tumor load)",
      units = "mm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline value, source column BSIZ in mm^2 (Table 2; log-normally distributed, range 100-55,500 mm^2). Reference 7000 mm^2 (Supplementary Table S5: SIZ = BSIZ/7000). Power exponent 0.261 on the time-dependent clearance CL_T (Table 3 theta20). The paper coded missing values (20%) as -99 or 0, and the NONMEM code sets SIZ = 1 (reference value) when BSIZ <= 0; the library model keeps that rule, so TUMSZ = 0 gives the reference-subject CL_T.",
      source_name = "BSIZ"
    ),
    STUDY_REACH = list(
      description = "REACH study indicator (1 = record from the phase III REACH study BO17072, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SAWYER part 2)",
      notes = "Source NONMEM code: ST17072 = 1 if STU = 17072 (Supplementary Table S5). Multiplies the residual-error SD by theta12 = 1.42 (Table 3 sigma_17072). Affects the residual error only; set 0 for new simulations.",
      source_name = "STU"
    ),
    STUDY_SAWYER_PART1 = list(
      description = "SAWYER part 1 (dose-finding cohort) indicator (1 = record from SAWYER BO25341 part 1, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SAWYER part 2)",
      notes = "Source NONMEM code: COHA = 1 if STU = 25341 and ARM = 1 (Supplementary Table S5), labelled 'effect of SAWYER Part 1 on residual error' in Table 3. Multiplies the residual-error SD by theta13 = 0.568. Affects the residual error only; set 0 for new simulations.",
      source_name = "COHA"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in place of BSA as the body-size covariate (Supplementary Table S1 run 032) and rejected on objective function."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL_inf, VC, ka and F_SC in the full covariate model (run 031); removed during backward elimination (runs 041, 044, 056)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Source column BALB. Tested on CL_inf in the full covariate model (run 031); removed in run 043."
    ),
    ADA_POS = list(
      description = "Anti-drug antibody status",
      units = "(binary)",
      type = "binary",
      notes = "Evaluated by diagnostic plots only (Supplementary Table S2); ADA in 13 SAWYER patients had no influence on the concentration-time courses (Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 255L,
    n_studies = 2L,
    n_observations = 4739L,
    age_range = "25-78 years (mean 58.4)",
    weight_range = "47.0-124 kg (mean 79.3)",
    bsa_range = "1.41-2.42 m^2 (mean 1.91)",
    bmi_range = "16.7-41.1 kg/m^2 (mean 27.3)",
    sex_female_pct = 32.5,
    race_ethnicity = "Caucasian 94.5%, American Indian or Alaska native 1.2%, other 2.4%, missing 2.0% (Table 2)",
    disease_state = "Chronic lymphocytic leukemia: previously untreated (SAWYER, 234 patients) or previously treated (REACH, 21 patients); rituximab combined with fludarabine and cyclophosphamide.",
    dose_range = "Intravenous 375 mg/m^2 in cycle 1 then 500 mg/m^2 in cycles 2-6; subcutaneous 1400, 1600 or 1870 mg (up to 2200 mg) in cycle 6 of SAWYER part 1, and 1600 mg in cycles 2-6 of SAWYER part 2 after a 375 mg/m^2 intravenous cycle 1. Six 28-day cycles.",
    regions = "International (two multicentre trials)",
    studies = "SAWYER (BO25341, NCT01292603, phase Ib two-part; 4158 samples from 234 patients, 2381 i.v. and 1777 s.c.) and REACH (BO17072, NCT00090051, phase III; 581 samples from 21 patients).",
    wbc_summary = "Baseline WBC mean 95.2 (SD 75.8) x 10^9/L, range 4.03-436 (Table 2)",
    tumor_size_summary = "Baseline tumor size (BSIZ) mean 6810 (SD 7540) mm^2, range 100-55,500 mm^2 (Table 2)",
    excluded_obs = "Samples below the 0.5 ug/mL LLOQ (n = 348, 7.3%) were excluded.",
    notes = "45.1% of patients received intravenous rituximab only and 54.9% received both routes (Table 2)."
  )

  ini({
    # Structural PK parameters (Table 3, final covariate model). NONMEM units
    # are mL and mL/day; converted to L and L/day here (divide by 1000).
    # Reference subject: BSA 1.9 m^2, BMI 27 kg/m^2, male, WBC 100 x 10^9/L,
    # BSIZ 7000 mm^2 (Supplementary Table S5 $PK).
    lcl_exp_kdes <- log(0.0399); label("Decay rate constant of the time-dependent clearance k_des (1/day)") # Table 3 theta1 k_des = 0.0399 1/day
    lcl_exp_component <- log(1.55); label("Initial time-dependent clearance CL_T at time 0 (L/day)") # Table 3 theta2 CL_T = 1550 mL/day
    lcl_exp_inf <- log(0.207); label("Time-independent clearance CL_inf (L/day)") # Table 3 theta3 CL_inf = 207 mL/day
    lvc <- log(4.99); label("Central volume of distribution VC (L)") # Table 3 theta4 VC = 4990 mL
    lvp <- log(3.70); label("Peripheral volume of distribution VP (L)") # Table 3 theta5 VP = 3700 mL
    lq <- log(0.420); label("Intercompartmental clearance Q (L/day)") # Table 3 theta6 Q = 420 mL/day
    lka <- log(0.372); label("Subcutaneous absorption rate constant ka (1/day)") # Table 3 theta7 ka = 0.372 1/day
    lfdepot <- log(0.633); label("Absolute subcutaneous bioavailability F_SC (fraction)") # Table 3 theta8 F_SC = 0.633

    # Covariate effects (Table 3; functional forms from Supplementary Table S5)
    e_bsa_cl_q <- 1.37; label("Power exponent of (BSA/1.9) on CL_T, CL_inf and Q (unitless)") # Table 3 theta14 CL_BSA = Q_BSA = 1.37
    e_bsa_vc_vp <- 0.8; label("Power exponent of (BSA/1.9) on VC and VP (unitless)") # Table 3 theta15 VC_BSA = VP_BSA = 0.8
    e_wbc_cl_exp_component <- 0.223; label("Power exponent of (WBC/100) on CL_T (unitless)") # Table 3 theta16 CL_T,WBC = 0.223
    e_tumsz_cl_exp_component <- 0.261; label("Power exponent of (BSIZ/7000) on CL_T (unitless)") # Table 3 theta20 CL_T,BSIZ = 0.261
    e_sexf_vc <- 0.909; label("Multiplicative factor on VC for female sex, theta^SEXF (unitless)") # Table 3 theta17 VC,SEXF = 0.909
    e_bmi_ka <- -1.01; label("Power exponent of (BMI/27) on ka (unitless)") # Table 3 theta18 ka,BMI = -1.01
    e_bmi_fdepot <- -0.465; label("Power exponent of (BMI/27) on F_SC (unitless)") # Table 3 theta19 F_SC,BMI = -0.465

    # Residual error (Table 3; Supplementary Information residual error model
    # and Supplementary Table S5 $ERROR). Log-scale SD
    # W = sdL - (sdL - sdH) * C / (sd50 + C), multiplied by the study factors
    # and by exp(etaRUV); sigma^2 is fixed to 1 in NONMEM, so W is the SD.
    sdL <- 0.81; label("Log-scale residual SD at low concentrations (unitless)") # Table 3 theta9 sigma_L = 0.81
    sdH <- 0.134; label("Log-scale residual SD at high concentrations (unitless)") # Table 3 theta10 sigma_H = 0.134
    sd50 <- 6.35; label("Concentration at which the residual SD is (sdL + sdH)/2 (ug/mL)") # Table 3 theta11 sigma_50 = 6.35 ug/mL
    e_study_reach_ruv <- 1.42; label("Multiplicative factor on the residual SD for the REACH study (unitless)") # Table 3 theta12 sigma_17072 = 1.42
    e_study_sawyer1_ruv <- 0.568; label("Multiplicative factor on the residual SD for SAWYER part 1 (unitless)") # Table 3 theta13 sigma_CohA = 0.568

    # Inter-individual variability (Table 3; omega^2 on the log scale)
    etalcl_exp_kdes ~ 0.357 # Table 3 Omega(1,1) k_des = 0.357 (CV 59.7%)
    etalcl_exp_component ~ 0.691 # Table 3 Omega(2,2) CL_T = 0.691 (CV 83.1%)
    etalcl_exp_inf + etalvc ~ c(0.106, 0.0277, 0.0323) # Table 3 Omega(3,3) = 0.106, Omega(3,4) = 0.0277, Omega(4,4) = 0.0323
    etalka + etalfdepot ~ c(0.115, 0.0265, 0.0453) # Table 3 Omega(5,5) = 0.115, Omega(5,6) = 0.0265, Omega(6,6) = 0.0453
    etaRUV ~ 0.0929 # Table 3 Omega(7,7) on the residual SD = 0.0929 (CV 30.5%)
  })

  model({
    # Covariate ratios (Supplementary Table S5 $PK). A non-positive tumor size
    # is the paper's missing-value code and maps to the reference value.
    ibsa <- BSA / 1.9
    ibmi <- BMI / 27
    iwbc <- WBC / 100
    isiz <- TUMSZ / 7000
    if (TUMSZ <= 0) {
      isiz <- 1
    }

    # Individual parameters
    cl_exp_kdes <- exp(lcl_exp_kdes + etalcl_exp_kdes)
    cl_exp_component <- exp(lcl_exp_component + etalcl_exp_component) *
      iwbc^e_wbc_cl_exp_component * isiz^e_tumsz_cl_exp_component * ibsa^e_bsa_cl_q
    cl_exp_inf <- exp(lcl_exp_inf + etalcl_exp_inf) * ibsa^e_bsa_cl_q
    vc <- exp(lvc + etalvc) * e_sexf_vc^SEXF * ibsa^e_bsa_vc_vp
    vp <- exp(lvp) * ibsa^e_bsa_vc_vp
    q <- exp(lq) * ibsa^e_bsa_cl_q
    ka <- exp(lka + etalka) * ibmi^e_bmi_ka
    fdepot <- exp(lfdepot + etalfdepot) * ibmi^e_bmi_fdepot

    # Total clearance decays from CL_T + CL_inf towards CL_inf (Methods;
    # Supplementary Table S5 $DES: CL = CLT*EXP(-KDES*T)+CLINF). time is the
    # time since the first rituximab dose, as TIME in the NONMEM dataset.
    cl <- cl_exp_component * exp(-cl_exp_kdes * time) + cl_exp_inf

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Subcutaneous doses enter depot; intravenous doses go to central.
    f(depot) <- fdepot

    # Dose in mg and volume in L give ug/mL (NONMEM S2 = VC/1000 with VC in mL).
    Cc <- central / vc

    # Log-scale residual error: Y = log(IPRED) + W * EPS with EPS ~ N(0, 1).
    w <- (sdL - (sdL - sdH) * Cc / (sd50 + Cc)) *
      e_study_reach_ruv^STUDY_REACH * e_study_sawyer1_ruv^STUDY_SAWYER_PART1 * exp(etaRUV)
    Cc ~ lnorm(w)
  })
}
