Wilkins_2019_avelumab <- function() {
  description <- "Two-compartment population PK model for avelumab (anti-PD-L1 IgG1) in adults with advanced solid tumors (14 tumor types, including metastatic Merkel cell carcinoma and urothelial carcinoma), with sigmoid time-varying clearance whose maximal change, T50 and shape differ for Merkel cell carcinoma and squamous cell carcinoma of the head and neck, and 13 covariate effects on clearance and 4 on central volume (Wilkins 2019)"
  reference <- "Wilkins JJ, Brockhaus B, Dai H, Vugmeyster Y, White JT, Brar S, Bello CL, Neuteboom B, Wade JR, Girard P, Khandelwal A. Time-Varying Clearance and Impact of Disease State on the Pharmacokinetics of Avelumab in Merkel Cell Carcinoma and Urothelial Carcinoma. CPT Pharmacometrics Syst Pharmacol. 2019;8(6):415-427. doi:10.1002/psp4.12406"
  vignette <- "Wilkins_2019_avelumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "avelumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "avelumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL and V1 normalised to the 71 kg population median (Data S1 '(WTT/71)'). A missing weight (<= 0 in the source data) is imputed to 71 kg, as in Data S1.",
      source_name = "WT_BL"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Linear categorical effects (1 + theta) on CL and V1 for SEX = 1, labelled 'Female sex' in Table 2.",
      source_name = "SEX"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL normalised to the 63-year population median (Data S1 '(AGE/63)'); no cap.",
      source_name = "AGE"
    ),
    RACE_BLACK = list(
      description = "Black or African American race, 1 = Black / African American, 0 = any other or missing race",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = all other races (White, Asian, American Indian / Alaska Native, Native Hawaiian / Pacific Islander, Other, missing)",
      notes = "Linear categorical effect on CL for RACE = 2 in Data S1 (n = 87, 4.8%).",
      source_name = "RACE"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL and V1 normalised to the 39 g/L median; capped at 61 g/L and a missing value (<= 0) imputed to 39 g/L (Data S1).",
      source_name = "ALB"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "eGFR (estimating equation not stated in the paper). Power effect on CL normalised to the 85.9 mL/min/1.73 m^2 median; capped at 226 and a missing value (<= 0) imputed to 85.9 (Data S1).",
      source_name = "EGFR"
    ),
    TUMSZ = list(
      description = "Baseline tumor size",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline tumor size in mm (Table 1; median 60 mm, range 10-750). Power effect on CL normalised to 60 mm; capped at 388 mm and a missing value (<= 0) imputed to 60 mm (Data S1).",
      source_name = "TS_BL"
    ),
    ADA_POS = list(
      description = "Treatment-emergent anti-drug antibody status ever positive, 1 = ADA-positive at any time point, 0 = never positive or missing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = never ADA-positive or ADA status missing",
      notes = "Time-fixed 'ever positive' flag (HAHAEVER in Data S1; n = 76, 4.2%). The Discussion notes the higher CL was evident from the first dose, before seroconversion.",
      source_name = "HAHAEVER"
    ),
    CRP = list(
      description = "Baseline C-reactive protein (standard assay)",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline-only value. Power effect on CL normalised to the 10.7 mg/L median; capped at 250 mg/L and a value <= 0 (missing, or a recorded 0) imputed to 10.7 mg/L (Data S1).",
      source_name = "CRP"
    ),
    PLT = list(
      description = "Baseline platelet count",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL normalised to the 247 x 10^9/L median; capped at 600 and a missing value (<= 0) imputed to 247 (Data S1).",
      source_name = "PLAT"
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL normalised to the 22 U/L median; capped at 105 U/L and a missing value (<= 0) imputed to 22 U/L (Data S1).",
      source_name = "AST"
    ),
    CONMED_OPIOID = list(
      description = "Concomitant opioid use, 1 = taking a concomitant opioid, 0 = not",
      units = "(binary)",
      type = "binary",
      reference_category = "1 = taking concomitant opioids (the most frequent category, 66.9%)",
      notes = "The source effect applies to patients NOT taking opioids (Data S1 'IF(CMOPI.EQ.0) CLCMOPI = 1 + THETA(21)'; Table 2 'Not taking concomitant opioids on CL'), so the model uses (1 - CONMED_OPIOID) to keep the published coefficient.",
      source_name = "CMOPI"
    ),
    PRIOR_BIO = list(
      description = "Previous use of biologics, 1 = previously treated with a biologic, 0 = biologic-naive",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no previous biologic",
      notes = "Linear categorical effects on CL and V1 (PMBIO = 1 in Data S1; n = 478, 26.2%).",
      source_name = "PMBIO"
    ),
    TUMTP_SCCHN = list(
      description = "Squamous cell carcinoma of the head and neck, 1 = SCCHN, 0 = other tumor type",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = any of the 12 tumor types without a time effect on CL (all TUMTP_* = 0)",
      notes = "Selects the SCCHN typical Imax (-0.284) and gamma (0.73) of the time-varying CL (Data S1 'TTYPE.EQ.6'); T50 is shared with the other tumor types. n = 153 (8.4%).",
      source_name = "TTYPE"
    ),
    TUMTP_MCC = list(
      description = "Metastatic Merkel cell carcinoma, 1 = mMCC, 0 = other tumor type",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = any of the 12 tumor types without a time effect on CL (all TUMTP_* = 0)",
      notes = "Selects the mMCC typical Imax (-0.387), T50 (131 days) and gamma (1.68) of the time-varying CL (Data S1 'TTYPE.EQ.14'). n = 88 (4.8%), all from JAVELIN Merkel 200.",
      source_name = "TTYPE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1827L,
    n_studies = 3L,
    n_observations = 10637L,
    age_range = "19-91 years (median 63)",
    weight_range = "30.4-204 kg (median 71)",
    sex_female_pct = 47.8,
    race_ethnicity = "White 76.2%, Asian 11.3%, Black or African American 4.8%, Other 7.1%, American Indian or Alaska Native 0.3%, Native Hawaiian or other Pacific Islander 0.2%, missing 0.2%",
    disease_state = "Metastatic or locally advanced solid tumors, 14 tumor types: NSCLC 18.6%, gastric / GEJ cancer 15.7%, urothelial carcinoma 13.6%, ovarian cancer 12.5%, metastatic breast cancer 9.2%, SCCHN 8.4%, Merkel cell carcinoma 4.8%, solid tumors (dose escalation) 3.8%, mesothelioma 2.9%, renal cell carcinoma 2.9%, melanoma 2.8%, adrenocortical carcinoma 2.7%, colorectal cancer 1.2%, castration-resistant prostate cancer 1.0%",
    dose_range = "1, 3, 10 or 20 mg/kg IV every 2 weeks (10 mg/kg in 97.3% of patients)",
    regions = "Global (JAVELIN Solid Tumor, JAVELIN Solid Tumor JPN in Japan, JAVELIN Merkel 200)",
    renal_function = "eGFR median 85.9 mL/min/1.73 m^2 (21.8-398); none 40.6%, mild 38.3%, moderate 20.7%, severe 0.2% renal impairment",
    notes = "Pooled from JAVELIN Solid Tumor (NCT01772004, n = 1688), JAVELIN Solid Tumor JPN (NCT01943461, n = 51) and JAVELIN Merkel 200 (NCT02155647, n = 88); Table 1. ADA ever positive 4.2%; concomitant opioids 66.9%; previous biologics 26.2%."
  )

  ini({
    # Structural PK (Table 2). Time in hours, as in the NONMEM data set.
    lcl <- log(0.0308); label("Baseline clearance CL for the reference patient (L/h)") # Table 2: CL = 0.0308 L/h
    lvc <- log(3.42); label("Central volume of distribution V1 for the reference patient (L)") # Table 2: Central volume = 3.42 L
    lvp <- log(0.918); label("Peripheral volume of distribution V2 (L)") # Table 2: Peripheral volume = 0.918 L
    lq <- log(0.0313); label("Intercompartmental clearance Q (L/h)") # Table 2: Intercompartmental CL = 0.0313 L/h

    # Sigmoid time-varying clearance (Eq. 1; Data S1):
    # CL(t) = CL * exp(Imax_i * t^gamma / (T50^gamma + t^gamma)), t = time after
    # first dose in days, with Imax_i = Imax + eta (additive IIV).
    cl_time_max <- fixed(0); label("Maximal change in log-CL, Imax, for the 12 tumor types without a time effect (unitless)") # Table 2: Imax for all other tumor types = 0, fixed
    cl_time_max_scchn <- -0.284; label("Maximal change in log-CL, Imax, for SCCHN (unitless)") # Table 2: Imax SCCHN = -0.284
    cl_time_max_mcc <- -0.387; label("Maximal change in log-CL, Imax, for mMCC (unitless)") # Table 2: Imax mMCC = -0.387
    lcl_t50 <- log(68.4); label("Time to half of the maximal CL change, T50, for all tumor types except mMCC (day)") # Table 2: T50 for all other tumor types = 68.4 days (also used for SCCHN, Data S1)
    lcl_t50_mcc <- log(131); label("Time to half of the maximal CL change, T50, for mMCC (day)") # Table 2: T50 mMCC = 131 days
    lcl_time_hill <- log(2.91); label("Shape of the time effect on CL, gamma, for the 12 other tumor types (unitless)") # Table 2: gamma for all other tumor types = 2.91
    lcl_time_hill_scchn <- log(0.73); label("Shape of the time effect on CL, gamma, for SCCHN (unitless)") # Table 2: gamma for SCCHN = 0.73
    lcl_time_hill_mcc <- log(1.68); label("Shape of the time effect on CL, gamma, for mMCC (unitless)") # Table 2: gamma for mMCC = 1.68

    # Covariate effects on CL (Table 2). Power exponents on (COV/median);
    # categorical effects are fractional changes, PAR * (1 + theta) (Eqs. 2-3).
    e_wt_cl <- 0.324; label("Power exponent of body weight on CL (unitless)") # Table 2: Body weight on CL = 0.324
    e_sexf_cl <- -0.15; label("Fractional change in CL for female sex (unitless)") # Table 2: Female sex on CL = -0.15
    e_age_cl <- -0.12; label("Power exponent of age on CL (unitless)") # Table 2: Age on CL = -0.12
    e_race_black_cl <- -0.12; label("Fractional change in CL for Black or African American race (unitless)") # Table 2: Black race on CL = -0.12
    e_alb_cl <- -0.438; label("Power exponent of albumin on CL (unitless)") # Table 2: Albumin on CL = -0.438
    e_crcl_cl <- 0.0579; label("Power exponent of eGFR on CL (unitless)") # Table 2: eGFR on CL = 0.0579
    e_tumsz_cl <- 0.0592; label("Power exponent of baseline tumor size on CL (unitless)") # Table 2: Tumor burden on CL = 0.0592
    e_ada_pos_cl <- 0.123; label("Fractional change in CL for ever ADA-positive (unitless)") # Table 2: Ever positive for ADA on CL = 0.123
    e_crp_cl <- 0.0677; label("Power exponent of C-reactive protein on CL (unitless)") # Table 2: CRP on CL = 0.0677
    e_plt_cl <- 0.0603; label("Power exponent of platelet count on CL (unitless)") # Table 2: Platelets on CL = 0.0603
    e_ast_cl <- -0.0514; label("Power exponent of AST on CL (unitless)") # Table 2: AST on CL = -0.0514
    e_no_opioid_cl <- -0.0445; label("Fractional change in CL when NOT taking concomitant opioids (unitless)") # Table 2: Not taking concomitant opioids on CL = -0.0445
    e_prior_bio_cl <- -0.0508; label("Fractional change in CL for previous use of biologics (unitless)") # Table 2: Previous use of biologics on CL = -0.0508

    # Covariate effects on V1 (Table 2).
    e_wt_vc <- 0.362; label("Power exponent of body weight on V1 (unitless)") # Table 2: Body weight on V1 = 0.362
    e_sexf_vc <- -0.160; label("Fractional change in V1 for female sex (unitless)") # Table 2: Female sex on V1 = -0.160
    e_alb_vc <- -0.278; label("Power exponent of albumin on V1 (unitless)") # Table 2: Albumin on V1 = -0.278
    e_prior_bio_vc <- -0.0526; label("Fractional change in V1 for previous use of biologics (unitless)") # Table 2: Previous use of biologics on V1 = -0.0526

    # IIV (Table 2). Log-normal BLOCK(3) on CL, V1, V2. The covariance values
    # are placed by the Data S1 $OMEGA BLOCK(3) positions: row 3 holds
    # cov(CL,V2) = 0.0633 and cov(V1,V2) = -0.014 (Table 2 prints these two
    # rows in NONMEM output order but labels them the other way round).
    etalcl + etalvc + etalvp ~ c(
      0.0535,
      0.0171, 0.0332,
      0.0633, -0.014, 0.858
    ) # Table 2: omega2 CL 0.0535, cov(CL,V1) 0.0171, omega2 V1 0.0332, cov 0.0633, cov -0.014, omega2 V2 0.858
    etacl_time_max ~ 0.0596 # Table 2: IIV on Imax (variance) = 0.0596; additive (normal) on Imax

    # Residual error (Table 2): combined, W = sqrt(prop^2 * F^2 + add^2), SIGMA 1 FIX
    propSd <- 0.162; label("Proportional residual error (fraction)") # Table 2: Proportional residual error = 0.162
    addSd <- 2.43; label("Additive residual error (ug/mL)") # Table 2: Additive residual error = 2.43 ug/mL
  })
  model({
    # Covariate handling as in Data S1: missing (<= 0) values imputed to the
    # population median and extreme values capped.
    wt_c <- WT
    if (WT <= 0) wt_c <- 71
    alb_c <- ALB
    if (ALB <= 0) alb_c <- 39
    if (ALB > 61) alb_c <- 61
    crcl_c <- CRCL
    if (CRCL <= 0) crcl_c <- 85.9
    if (CRCL > 226) crcl_c <- 226
    tumsz_c <- TUMSZ
    if (TUMSZ <= 0) tumsz_c <- 60
    if (TUMSZ > 388) tumsz_c <- 388
    crp_c <- CRP
    if (CRP <= 0) crp_c <- 10.7
    if (CRP > 250) crp_c <- 250
    plt_c <- PLT
    if (PLT <= 0) plt_c <- 247
    if (PLT > 600) plt_c <- 600
    ast_c <- AST
    if (AST <= 0) ast_c <- 22
    if (AST > 105) ast_c <- 105

    clcov <- (wt_c / 71)^e_wt_cl *
      (1 + e_sexf_cl * SEXF) *
      (AGE / 63)^e_age_cl *
      (1 + e_race_black_cl * RACE_BLACK) *
      (alb_c / 39)^e_alb_cl *
      (crcl_c / 85.9)^e_crcl_cl *
      (tumsz_c / 60)^e_tumsz_cl *
      (1 + e_ada_pos_cl * ADA_POS) *
      (crp_c / 10.7)^e_crp_cl *
      (plt_c / 247)^e_plt_cl *
      (ast_c / 22)^e_ast_cl *
      (1 + e_no_opioid_cl * (1 - CONMED_OPIOID)) *
      (1 + e_prior_bio_cl * PRIOR_BIO)
    vccov <- (wt_c / 71)^e_wt_vc *
      (1 + e_sexf_vc * SEXF) *
      (alb_c / 39)^e_alb_vc *
      (1 + e_prior_bio_vc * PRIOR_BIO)

    # Tumor-type-specific time-varying CL parameters (Data S1 TTYPE 6 = SCCHN,
    # TTYPE 14 = mMCC); T50 for SCCHN is the shared value.
    tv_cl_time_max <- cl_time_max
    tv_cl_t50 <- exp(lcl_t50)
    tv_cl_time_hill <- exp(lcl_time_hill)
    if (TUMTP_SCCHN == 1) {
      tv_cl_time_max <- cl_time_max_scchn
      tv_cl_time_hill <- exp(lcl_time_hill_scchn)
    }
    if (TUMTP_MCC == 1) {
      tv_cl_time_max <- cl_time_max_mcc
      tv_cl_t50 <- exp(lcl_t50_mcc)
      tv_cl_time_hill <- exp(lcl_time_hill_mcc)
    }
    cl_time_max_i <- tv_cl_time_max + etacl_time_max

    # Time after first dose in days (model time is hours from the first dose)
    tday <- t / 24
    td_cl <- exp(cl_time_max_i * tday^tv_cl_time_hill / (tv_cl_t50^tv_cl_time_hill + tday^tv_cl_time_hill))

    cl_base <- exp(lcl + etalcl) * clcov
    cl <- cl_base * td_cl
    vc <- exp(lvc + etalvc) * vccov
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and V1 in L, so central / vc is mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
