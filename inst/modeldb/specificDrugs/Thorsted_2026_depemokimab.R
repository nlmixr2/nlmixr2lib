Thorsted_2026_depemokimab <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption for depemokimab (ultra-long-acting anti-IL-5 IgG1) in healthy volunteers and adults and adolescents with asthma or chronic rhinosinusitis with nasal polyps (CRSwNP), coupled to an indirect-response model of blood eosinophil count with sigmoid inhibition of production and a linear placebo drift of the baseline (Thorsted 2026). Allometric weight effects on CL/F and V/F; albumin, eGFR, CRSwNP and study effects on CL/F; Asian race and healthy-volunteer effects on V/F; age, abdominal injection site and study effects on ka; CRSwNP and study effects on relative bioavailability; baseline count, CRSwNP, age, weight and Asian race effects on the logit of Emax."
  reference <- "Thorsted A, Tanneau L, Lavalley-Morelle A, Follows R, Jacques L, Bird N, Gevaert P, Pavord I, Ribbing J, Gupta A, Howarth P, Schalkwijk S. Population Pharmacokinetics and Pharmacodynamics of Depemokimab in People with Asthma and Chronic Rhinosinusitis with Nasal Polyps. J Clin Pharmacol. 2026;66(8):e70187. doi:10.1002/jcph.70187"
  vignette <- "Thorsted_2026_depemokimab"
  paper_specific_etas <- c("etaruv_cc", "etaruv_eos")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "depemokimab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "depemokimab", units = "mg", specimen = "plasma", verified = TRUE),
    blood_eosinophils = list(analyte = "eosinophils", units = "cells/uL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK: power effects (WT/76)^0.841 on CL/F and (WT/76)^0.887 on V/F; 76 kg is the reference subject of Figure 2 (the median of the PK population). PD: linear effect -0.00471 per kg on the logit of Emax, centred on 77.4 kg (Table S4 footnote b).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "PK: exponential effect exp(-0.00735 * (AGE - 50)) on ka (50 years = Figure 2 reference subject). PD: linear effect +0.00660 per year on the logit of Emax, centred on 54 years (Table S4 footnote b).",
      source_name = "AGE"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect exp(-0.00684 * (ALB - 45)) on CL/F; 45 g/L is the Figure 2 reference subject.",
      source_name = "albumin"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "eGFR (estimating equation not stated in the paper). Exponential effect exp(0.00134 * (CRCL - 94)) on CL/F; 94 mL/min/1.73 m^2 is the Figure 2 reference subject.",
      source_name = "eGFR"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator (1 = Asian, 0 = other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "PK: fractional change +0.0788 on V/F. PD: additive effect on the logit of Emax that differs by population, +0.450 in asthma and +0.187 in CRSwNP (Table S4).",
      source_name = "Asian race"
    ),
    INJSITE_ARM = list(
      description = "Injection site = upper arm (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; the paper's reference site is the upper arm, so set INJSITE_ARM = 1 for the reference subject",
      notes = "The paper estimated only an abdomen effect on ka (+0.520, fractional change relative to upper arm or thigh). The model forms the abdomen indicator as 1 - INJSITE_ARM - INJSITE_THIGH, following the register convention that all site indicators equal to 0 means abdomen.",
      source_name = "injection site"
    ),
    INJSITE_THIGH = list(
      description = "Injection site = thigh (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Thigh injections share the upper-arm absorption rate (no thigh effect in Table S3); the column only enters through the abdomen indicator 1 - INJSITE_ARM - INJSITE_THIGH.",
      source_name = "injection site"
    ),
    DIS_CHRONIC_RHINOSINUSITIS_NASAL_POLYPS = list(
      description = "CRSwNP participant population (1 = CRSwNP, 0 = asthma or healthy volunteer)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (asthma participant population)",
      notes = "PK: fractional changes +0.0625 on CL/F and -0.0509 on Frel. PD: +0.320 on the logit of Emax, and selects the CRSwNP Asian-race effect on Emax (Table S3, Table S4). The ANCHOR-1/-2 participants.",
      source_name = "CRSwNP participant population"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer participant population (1 = healthy volunteer, 0 = patient)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (asthma participant population)",
      notes = "Fractional changes -0.0513 on V/F and +0.588 on ka (Table S3 'HV participant population'). Set to 1 for both the US HV and China HV studies; the China HV study adds its own STUDY_208021 effects on top (Table S1 footnote e). Healthy volunteers were not part of the PK/PD analysis.",
      source_name = "HV participant population"
    ),
    STUDY_205722 = list(
      description = "First-time-in-human asthma study 205722 (NCT03287310) participant (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (SWIFT-1/-2 within the asthma population)",
      notes = "PK: fractional changes +0.155 on CL/F, -0.173 on Frel, +0.627 on ka. PD: fractional change -0.561 on EC50, additive +0.357 on the placebo slope and fractional change -0.242 on the residual SD of the eosinophil count (Tables S3, S4).",
      source_name = "FTIH Asthma study population"
    ),
    STUDY_208021 = list(
      description = "China healthy-volunteer study 208021 (NCT05140200) participant (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (US HV study within the healthy-volunteer population)",
      notes = "Fractional changes -0.219 on CL/F and +0.405 on ka, applied on top of the DIS_HEALTHY effects (Table S3 'China HV study population').",
      source_name = "China HV study population"
    ),
    DOSE_LOW = list(
      description = "2 mg dose-level cohort of the first-time-in-human study (1 = yes, 0 = other dose)",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Fractional change -0.349 on CL/F (Table S3 '2 mg dose on CL'), tested because the 2 mg cohort had higher dose-normalised concentrations (Table S1 footnote f).",
      source_name = "2 mg dose level"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 961L,
    n_studies = 7L,
    age_range = "12-93 years (PK population mean 49.8 years)",
    weight_range = "34.6-161 kg (PK population mean 78.0 kg)",
    sex_female_pct = NA_real_,
    race_ethnicity = "19% Asian in the PK population (182 of 961)",
    disease_state = "Healthy volunteers (160), asthma (529; mild-to-moderate in the first-time-in-human study, type 2 asthma with exacerbations in SWIFT-1/-2) and inadequately controlled CRSwNP (272; ANCHOR-1/-2)",
    dose_range = "Single SC 2-300 mg (Phase I) or two SC 100 mg doses 26 weeks apart (Phase III)",
    regions = "UK, Germany, China, US and global Phase III sites",
    n_observations = "9383 plasma concentrations (PK); 14,273 blood eosinophil counts (PK/PD)",
    pd_population = "1324 participants (depemokimab or placebo) from the first-time-in-human asthma study, SWIFT-1/-2 and ANCHOR-1/-2; mean baseline blood eosinophil count 440 cells/uL (Table 3)",
    notes = "Pooled analysis of three Phase I studies (205722 first-time-in-human asthma, 208021 China HV, 214099 US HV) and four Phase III studies (SWIFT-1/-2, ANCHOR-1/-2); Table 1 (PK) and Table 3 (PK/PD). The PK/PD model was fitted sequentially on individual predicted concentrations."
  )

  ini({
    # ---- PK: Table S3 (final PK model, run 11) ----
    # Reference subject (Figure 2 caption): 76 kg, 50 years, non-Asian asthma
    # participant (SWIFT-1/-2), albumin 45 g/L, eGFR 94 mL/min/1.73 m^2,
    # upper-arm injection. Categorical effects are fractional differences,
    # TV * (1 + theta); continuous non-weight effects are exponential,
    # TV * exp(theta * (cov - ref)) (Supplementary Methods).
    lcl <- log(0.0920); label("Apparent clearance CL/F for the reference subject (L/day)") # Table S3 CL/F 0.0920 L/day
    lvc <- log(6.29); label("Apparent volume of distribution V/F for the reference subject (L)") # Table S3 V/F 6.29 L
    lka <- log(0.212); label("First-order absorption rate constant ka for the reference subject (1/day)") # Table S3 ka 0.212 /day
    lfdepot <- fixed(log(1)); label("Relative bioavailability Frel (fraction)") # Table S3 Frel 1.00 [FIX]

    e_wt_cl <- 0.841; label("Allometric exponent of body weight on CL/F (unitless)") # Table S3 'WT on CL' 0.841
    e_wt_vc <- 0.887; label("Allometric exponent of body weight on V/F (unitless)") # Table S3 'WT on V' 0.887
    e_alb_cl <- -0.00684; label("Exponential effect of albumin on CL/F (1/(g/L))") # Table S3 'Albumin on CL' -0.00684
    e_crcl_cl <- 0.00134; label("Exponential effect of eGFR on CL/F (1/(mL/min/1.73 m^2))") # Table S3 'eGFR on CL' 0.00134
    e_dis_chronic_rhinosinusitis_nasal_polyps_cl <- 0.0625; label("Fractional change in CL/F for the CRSwNP population (unitless)") # Table S3 'CRSwNP participant population on CL' 0.0625
    e_study_205722_cl <- 0.155; label("Fractional change in CL/F for the first-time-in-human asthma study (unitless)") # Table S3 'FTIH Asthma study population on CL' 0.155
    e_study_208021_cl <- -0.219; label("Fractional change in CL/F for the China HV study (unitless)") # Table S3 'China HV study population on CL' -0.219
    e_dose_low_cl <- -0.349; label("Fractional change in CL/F for the 2 mg dose cohort (unitless)") # Table S3 '2 mg dose on CL' -0.349
    e_race_asian_vc <- 0.0788; label("Fractional change in V/F for Asian race (unitless)") # Table S3 'Asian race on V' 0.0788
    e_dis_healthy_vc <- -0.0513; label("Fractional change in V/F for the healthy-volunteer population (unitless)") # Table S3 'HV participant population on V' -0.0513
    e_age_ka <- -0.00735; label("Exponential effect of age on ka (1/year)") # Table S3 'Age on ka' -0.00735
    e_injsite_abdomen_ka <- 0.520; label("Fractional change in ka for an abdominal injection (unitless)") # Table S3 'Abdomen injection site on ka' 0.520
    e_dis_healthy_ka <- 0.588; label("Fractional change in ka for the healthy-volunteer population (unitless)") # Table S3 'HV participant population on ka' 0.588
    e_study_205722_ka <- 0.627; label("Fractional change in ka for the first-time-in-human asthma study (unitless)") # Table S3 'FTIH Asthma study population on ka' 0.627
    e_study_208021_ka <- 0.405; label("Fractional change in ka for the China HV study (unitless)") # Table S3 'China HV study population on ka' 0.405
    e_dis_chronic_rhinosinusitis_nasal_polyps_fdepot <- -0.0509; label("Fractional change in Frel for the CRSwNP population (unitless)") # Table S3 'CRSwNP participant population on Frel' -0.0509
    e_study_205722_fdepot <- -0.173; label("Fractional change in Frel for the first-time-in-human asthma study (unitless)") # Table S3 'FTIH Asthma study population on Frel' -0.173

    # Table S3 reports IIV as SDs of exponential etas, so omega^2 = SD^2.
    etalcl ~ 0.101^2 # Table S3 'IIV CL' 0.101
    etalfdepot + etalka ~ c(
      0.167^2,
      0.638 * 0.167 * 0.351, 0.351^2
    ) # Table S3 'IIV Frel' 0.167, 'IIV ka' 0.351, 'Correlation IIV Frel - IIV ka' 0.638
    etaruv_cc ~ 0.488^2 # Table S3 'IIV RUV' 0.488 (exponential IIV on the residual SD)

    expSd <- 0.178; label("Residual SD of log concentration (log ug/mL)") # Table S3 'RUV' 0.178 (additive on the log scale)
    e_tadlt3h_expSd <- 5.18; label("Fractional change in the concentration residual SD within 3 h after a dose (unitless)") # Table S3 'TAD <3 h on RUV' 5.18

    # ---- PD: Table S4 (final blood eosinophil count model, run 18) ----
    # Reference subject (Table S4 footnote b): non-Asian SWIFT-1/-2 asthma
    # participant, predicted baseline 317 cells/uL, 54 years, 77.4 kg.
    lrbase <- log(317); label("Baseline blood eosinophil count (cells/uL)") # Table S4 'Baseline' 317 cells/uL
    slope_placebo <- -0.0466; label("Placebo drift of the baseline count (fractional change per year)") # Table S4 'Pslp' -0.0466 per year
    lkout <- log(log(2) / 3.08); label("Eosinophil turnover rate constant kout = ln2 / HLonset (1/day)") # Table S4 'HLonset' 3.08 day
    logitemax <- log(0.848 / (1 - 0.848)); label("Logit of the maximum fractional inhibition of eosinophil production Emax (logit, unitless)") # Table S4 'Emax' 0.848
    lec50 <- log(0.194); label("Depemokimab concentration at half-maximal effect EC50 (ug/mL)") # Table S4 'EC50' 0.194 ug/mL
    lhill <- log(1.64); label("Hill coefficient (unitless)") # Table S4 'Hill coefficient' 1.64

    e_lrbase_emax <- 0.834; label("Effect of log predicted baseline count on logit Emax (per log unit)") # Table S4 'Predicted baseline blood eosinophils on Emax' 0.834
    e_dis_chronic_rhinosinusitis_nasal_polyps_emax <- 0.320; label("Additive effect of the CRSwNP population on logit Emax") # Table S4 'CRSwNP participant population on Emax' 0.320
    e_age_emax <- 0.00660; label("Additive effect of age on logit Emax (per year)") # Table S4 'Age on Emax' 0.00660
    e_wt_emax <- -0.00471; label("Additive effect of body weight on logit Emax (per kg)") # Table S4 'Weight on Emax' -0.00471
    e_race_asian_asthma_emax <- 0.450; label("Additive effect of Asian race in the asthma population on logit Emax") # Table S4 'Asian race in Asthma population on Emax' 0.450
    e_race_asian_crswnp_emax <- 0.187; label("Additive effect of Asian race in the CRSwNP population on logit Emax") # Table S4 'Asian race in CRSwNP population on Emax' 0.187
    e_study_205722_ec50 <- -0.561; label("Fractional change in EC50 for the first-time-in-human asthma study (unitless)") # Table S4 'FTIH Asthma study population on EC50' -0.561
    e_study_205722_slope_placebo <- 0.357; label("Additive change in the placebo slope for the first-time-in-human asthma study (per year)") # Table S4 'FTIH Asthma study population on Pslp' 0.357

    etalrbase ~ 0.646^2 # Table S4 'IIV baseline' 0.646
    etalogitemax ~ 0.497^2 # Table S4 'IIV Emax' 0.497 (SD on the logit scale)
    etaruv_eos ~ 0.396^2 # Table S4 'IIV RUV' 0.396 (exponential IIV on the residual SD)

    expSd_blood_eosinophils <- 0.426; label("Residual SD of log blood eosinophil count (log cells/uL)") # Table S4 'RUV' 0.426
    e_study_205722_expSd_blood_eosinophils <- -0.242; label("Fractional change in the eosinophil residual SD for the first-time-in-human asthma study (unitless)") # Table S4 'Change in RUV for FTIH Asthma study population' -0.242
  })

  model({
    crswnp <- DIS_CHRONIC_RHINOSINUSITIS_NASAL_POLYPS
    asthma_pop <- 1 - crswnp
    injsite_abdomen <- 1 - INJSITE_ARM - INJSITE_THIGH

    # ---- PK ----
    cl <- exp(lcl + etalcl) * (WT / 76)^e_wt_cl *
      exp(e_alb_cl * (ALB - 45)) * exp(e_crcl_cl * (CRCL - 94)) *
      (1 + e_dis_chronic_rhinosinusitis_nasal_polyps_cl * crswnp) *
      (1 + e_study_205722_cl * STUDY_205722) *
      (1 + e_study_208021_cl * STUDY_208021) *
      (1 + e_dose_low_cl * DOSE_LOW)
    vc <- exp(lvc) * (WT / 76)^e_wt_vc *
      (1 + e_race_asian_vc * RACE_ASIAN) *
      (1 + e_dis_healthy_vc * DIS_HEALTHY)
    ka <- exp(lka + etalka) * exp(e_age_ka * (AGE - 50)) *
      (1 + e_injsite_abdomen_ka * injsite_abdomen) *
      (1 + e_dis_healthy_ka * DIS_HEALTHY) *
      (1 + e_study_205722_ka * STUDY_205722) *
      (1 + e_study_208021_ka * STUDY_208021)
    fdepot <- exp(lfdepot + etalfdepot) *
      (1 + e_dis_chronic_rhinosinusitis_nasal_polyps_fdepot * crswnp) *
      (1 + e_study_205722_fdepot * STUDY_205722)
    kel <- cl / vc

    # ---- PD (Figure S4): indirect response, inhibition of production ----
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    hill <- exp(lhill)
    ec50 <- exp(lec50) * (1 + e_study_205722_ec50 * STUDY_205722)
    # Logit-additive covariate model (Supplementary Methods Equation 2); the
    # baseline effect acts on the individual predicted baseline, centred on
    # log(317) (Table S4 footnote c).
    emax <- expit(logitemax + etalogitemax +
      e_lrbase_emax * (log(rbase) - log(317)) +
      e_dis_chronic_rhinosinusitis_nasal_polyps_emax * crswnp +
      e_age_emax * (AGE - 54) + e_wt_emax * (WT - 77.4) +
      e_race_asian_asthma_emax * RACE_ASIAN * asthma_pop +
      e_race_asian_crswnp_emax * RACE_ASIAN * crswnp)
    slope_pbo <- slope_placebo + e_study_205722_slope_placebo * STUDY_205722
    # Placebo count: linear proportional drift of the baseline, per year.
    eos_placebo <- rbase * (1 + slope_pbo * t / 365.25)
    d_eos_placebo <- rbase * slope_pbo / 365.25

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- fdepot

    Cc <- central / vc

    # The paper's state is the count as a ratio of placebo, R, with
    # dR/dt = kout * (1 - inh) - kout * R and R(0) = 1. The library state is
    # the absolute count B = eos_placebo * R, whose exact derivative is
    # dB/dt = (d_eos_placebo / eos_placebo) * B + kout * (eos_placebo * (1 - inh) - B).
    inh <- emax * Cc^hill / (Cc^hill + ec50^hill)
    d/dt(blood_eosinophils) <- d_eos_placebo / eos_placebo * blood_eosinophils +
      kout * (eos_placebo * (1 - inh) - blood_eosinophils)
    blood_eosinophils(0) <- rbase

    # Residual error: additive on the log scale, with exponential IIV on its
    # magnitude; the PK SD is 6.18-fold within 3 h of a dose.
    tad_depot <- tad(depot)
    early_obs <- (tad_depot < 3 / 24)
    sd_cc <- expSd * (1 + e_tadlt3h_expSd * early_obs) * exp(etaruv_cc)
    sd_eos <- expSd_blood_eosinophils *
      (1 + e_study_205722_expSd_blood_eosinophils * STUDY_205722) * exp(etaruv_eos)

    Cc ~ lnorm(sd_cc)
    blood_eosinophils ~ lnorm(sd_eos)
  })
}
