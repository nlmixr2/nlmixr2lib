Xiong_2022_tepotinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral tepotinib with sequential",
    "zero- then first-order absorption and first-order elimination, linked to a",
    "two-compartment model for its major circulating metabolite MSC2571109A",
    "(formed from parent clearance, fraction metabolised fixed to 1), in",
    "patients with cancer (including MET exon 14 skipping NSCLC from VISION)",
    "and healthy participants pooled from 12 studies. Dose-dependent relative",
    "bioavailability; food, formulation, hepatic dysfunction, eGFR, tumour",
    "type, opioid co-medication, INR and serum albumin covariates."
  )
  reference <- paste(
    "Xiong W, Papasouliotis O, Jonsson EN, Strotmann R, Girard P.",
    "Population pharmacokinetic analysis of tepotinib, an oral MET kinase",
    "inhibitor, including data from the VISION study.",
    "Cancer Chemother Pharmacol. 2022;89(5):655-669.",
    "doi:10.1007/s00280-022-04423-5. PMCID PMC9054876.",
    "Parameter values from Table 3; covariate-model functional forms from",
    "Electronic Supplementary Material ESM 13."
  )
  vignette <- "Xiong_2022_tepotinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description = "Labelled oral tepotinib dose level at the dose record (tepotinib hydrochloride hydrate)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying (use case a/b of the register: the current administered",
        "dose level). This is the LABELLED dose, e.g. 500 for the approved",
        "500 mg tablet, and it drives only the linear dose effect on relative",
        "bioavailability, (1 + e_dose_fdepot / 100 * (DOSE - 500)) per ESM 13.",
        "The dosing amount (amt) must instead be the free-base amount, 0.9 x",
        "DOSE (450 mg for a 500 mg tablet), because Table 3 footnote a states",
        "that CL/F, Vc/F, Q/F and Vp/F were multiplied by 0.9 to correct for",
        "the salt-to-base molar weight ratio. Studied range 30-1400 mg/day."
      ),
      source_name = "DOSE"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on relative bioavailability, reference 72 kg (ESM 13; Table 2 median 72.0 kg).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on parent central volume, reference 59 years (ESM 13).",
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Baseline estimated glomerular filtration rate (MDRD equation), BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "MDRD eGFR (Methods, Covariate model development). Power effects on",
        "parent CL/F and metabolite CL, reference 97.28 mL/min/1.73 m^2",
        "(ESM 13; Table 3 footnote rounds it to 97.3)."
      ),
      source_name = "eGFR"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on parent Q/F, reference 40 g/L (ESM 13). The Table 3",
        "footnote's '4 g/L' is a typographical slip for 40 g/L; ESM 13, the",
        "Figure 2 / Figure 3 captions and the ESM 2 note all give 40 g/L."
      ),
      source_name = "SALB"
    ),
    INR_BASE = list(
      description = "Baseline international normalized ratio of prothrombin time",
      units = "(unitless)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on parent Q/F, reference 1.06 (ESM 13).",
      source_name = "INR"
    ),
    HEPIMP = list(
      description = "Hepatic dysfunction indicator, NCI-ODWG class > 0",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NCI-ODWG class 0, normal hepatic function)",
      notes = paste(
        "Fractional effects on parent D1 and relative bioavailability and on",
        "metabolite Vc. Handled as time-varying in the covariate search; the",
        "baseline value was used for the forest plots (ESM 1)."
      ),
      source_name = "NCI ODG class > 0"
    ),
    DIS_HEALTHY = list(
      description = "Healthy participant (1) versus patient with cancer (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with cancer)",
      notes = paste(
        "Fractional effects on parent Vp/F and metabolite Vp; also selects the",
        "healthy-participant IIV variances for parent CL/F and F and for",
        "metabolite CL (Table 3)."
      ),
      source_name = "Patient/participant"
    ),
    TUMTP_HCC = list(
      description = "Hepatocellular carcinoma patient",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Fractional effects on parent CL/F and on the fraction metabolised to MSC2571109A (Table 3, ESM 13).",
      source_name = "HCC"
    ),
    TUMTP_CRC = list(
      description = "Colorectal cancer patient",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Fractional effect on parent CL/F (Table 3, ESM 13).",
      source_name = "Colorectal cancer"
    ),
    TUMTP_NSCLC = list(
      description = "Non-small cell lung cancer patient",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Fractional effects on parent Vc/F and metabolite CL (Table 3, ESM 13).",
      source_name = "NSCLC"
    ),
    CONMED_OPIOID = list(
      description = "Concomitant mu-opioid analgesic",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant opioid)",
      notes = "Time-varying (Methods). Fractional effect on parent CL/F (Table 3, ESM 13).",
      source_name = "mu-opioids"
    ),
    STUDY_MS200095_0028 = list(
      description = "Record from study MS200095-0028 (hepatic impairment study)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any other study)",
      notes = "Fractional effect on parent CL/F (Table 3, ESM 13). Set to 0 for simulation of new patients.",
      source_name = "Study MS200095-0028"
    ),
    FED = list(
      description = "Fed (1) versus fasted (0) at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (fed, non-high-fat standard breakfast) is the model reference",
      notes = paste(
        "Per dose record. The paper parameterises the FASTING effect, so",
        "fasting multiplies ka, D1 and F by (1 + effect) when FED = 0. The",
        "model reference is fed with a non-high-fat meal (Table 3 footnote)."
      ),
      source_name = "Fasting state"
    ),
    FED_HIGHFAT = list(
      description = "High-fat / high-calorie meal at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-high-fat meal or fasted)",
      notes = "Per dose record. Fractional effect on relative bioavailability (Table 3, ESM 13).",
      source_name = "High-fat meal"
    ),
    FORM_TEPOTINIB_CF1 = list(
      description = "Tepotinib capsule formulation 1 (non-micronised drug substance)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference formulation TF2; CF2 and the other unlisted formulations carry no effect)",
      notes = paste(
        "Time-varying per dose record. Fractional effects on F and ka, and",
        "selects the CF1-specific IIV variance on F (Table 3). Used only in",
        "the first-in-human study 001 at doses up to 230 mg."
      ),
      source_name = "CF1"
    ),
    FORM_TEPOTINIB_TF1 = list(
      description = "Tepotinib tablet formulation 1 (micronised drug substance)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference formulation TF2)",
      notes = "Time-varying per dose record. Fractional effect on ka (Table 3, ESM 13).",
      source_name = "TF1"
    ),
    FORM_TEPOTINIB_TF1FINE = list(
      description = "Tepotinib tablet formulation 1 made with finely micronised drug substance (TF1*)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference formulation TF2)",
      notes = "Time-varying per dose record. Fractional effect on ka (Table 3, ESM 13). Mutually exclusive with FORM_TEPOTINIB_TF1.",
      source_name = "TF1*"
    ),
    FORM_TEPOTINIB_TF3 = list(
      description = "Tepotinib tablet formulation 3 (the marketed tablet)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference formulation TF2)",
      notes = "Time-varying per dose record. Fractional effect on F (Table 3, ESM 13). Figure 3 simulates TF3 500 mg QD with food.",
      source_name = "TF3"
    ),
    RACE_ASIAN_NORTHEAST = list(
      description = "East Asian (Japanese or other East Asian) participant",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-East Asian)",
      notes = paste(
        "Fractional effect on metabolite Q (Table 3 'East Asian on Qmet').",
        "Table 2 race categories Japanese plus Other East Asian; the",
        "reference participant is 'non-East Asian' (Table 3 footnote)."
      ),
      source_name = "East Asian"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search, not retained (Methods; Discussion 'Other intrinsic factors')."
    ),
    CONMED_GEFITINIB = list(
      description = "Concomitant gefitinib (study 006)",
      units = "(binary)",
      type = "binary",
      notes = "Screened, no statistically significant effect on tepotinib PK (Results)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tepotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tepotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tepotinib", units = "mg", specimen = "tissue", verified = TRUE),
    central_msc2571109a = list(analyte = "MSC2571109A", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_msc2571109a = list(analyte = "MSC2571109A", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 613L,
    n_studies = 12L,
    age_range = "18-89 years",
    age_median = "58 years",
    weight_range = "35.5-136 kg",
    weight_median = "72.0 kg",
    sex_female_pct = 28.4,
    race_ethnicity = c(
      Caucasian = 59.1,
      Japanese = 4.6,
      `Other East Asian` = 23.7,
      `African origin` = 2.8,
      Hispanic = 4.1,
      `Other/missing` = 5.9
    ),
    disease_state = paste(
      "438 patients with cancer (NSCLC including MET exon 14 skipping NSCLC",
      "from the pivotal VISION study, hepatocellular carcinoma, colorectal,",
      "renal cell, gastroesophageal and other solid tumours) and 175 healthy",
      "participants, some with mild or moderate hepatic impairment."
    ),
    dose_range = paste(
      "30-1400 mg/day oral tepotinib (QD, three times weekly, or single dose)",
      "as capsules CF1/CF2 or tablets TF1/TF1*/TF2/TF3; therapeutic regimen",
      "500 mg QD (450 mg free base) with food."
    ),
    regions = "Multi-regional: Europe, United States, Japan, China and other Asian countries.",
    notes = paste(
      "Table 2: 613 participants with 10,788 tepotinib concentrations and 464",
      "participants with 7197 MSC2571109A concentrations; 439 male / 174",
      "female. Baseline eGFR mean 99.8 (SD 26.2, range 39.4-236)",
      "mL/min/1.73 m^2; 146 mild and 16 moderate NCI-ODWG hepatic impairment.",
      "Reference participant (Table 3 footnote): 59-year-old non-East Asian",
      "patient, 72 kg, eGFR 97.3, INR 1.06, NCI-ODG class 0, serum albumin",
      "40 g/L, no opioids, 500 mg TF2 with a non-high-fat meal."
    )
  )

  ini({
    # ---- Tepotinib structural parameters (Table 3; apparent values already
    #      multiplied by 0.9 for the salt-to-base ratio, Table 3 footnote a) ----
    lcl <- log(20.4); label("Apparent clearance CL/F (L/h)") # Table 3 'CLpar/F' = 20.4 L/h (RSE 2.07%); ESM 13 20.4
    lvc <- log(1020); label("Apparent central volume Vc/F (L)") # Table 3 'Vc,par/F' = 1020 L (RSE 2.00%); ESM 13 prints 1030
    lka <- log(0.278); label("First-order absorption rate constant ka (1/h)") # Table 3 'ka' = 0.278 1/h (RSE 6.16%); ESM 13 prints 1.47, falsified by Figure 3 (see vignette)
    lq <- log(1.32); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 3 'Qpar/F' = 1.32 L/h (RSE 4.22%)
    lvp <- log(1180); label("Apparent peripheral volume Vp/F (L)") # Table 3 'Vp,par/F' = 1180 L (RSE 16.6%)
    ld1 <- log(4.09); label("Zero-order absorption duration D1 (h)") # Table 3 'D1' = 4.09 h (RSE 5.34%)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F (fraction)") # Table 3 'Relative Fpar' = 1.00 (FIX)

    # ---- Tepotinib covariate effects (fractional unless noted; ESM 13) ----
    e_fasted_d1 <- -0.370; label("Fasting effect on D1 (fractional change)") # Table 3 'Fasting state covariate on D1' = -0.370
    e_hepimp_d1 <- -0.332; label("NCI-ODG class > 0 effect on D1 (fractional change)") # Table 3 'NCI-ODG class > 0 covariate (liver dysfunction) on D1' = -0.332
    e_dose_fdepot <- -0.0412; label("Linear dose effect on F (fractional change per 100 mg above 500 mg)") # Table 3 'DOSE covariate on Fpar (/100 mg)' = -0.0412; ESM 13 1 + (-0.0412/100)*(DOSE - 500)
    e_fasted_fdepot <- -0.209; label("Fasting effect on F (fractional change)") # Table 3 'Fasting state covariate on Fpar' = -0.209
    e_fed_highfat_fdepot <- 0.320; label("High-fat meal effect on F (fractional change)") # Table 3 'High-fat meal covariate on Fpar' = 0.320
    e_form_cf1_fdepot <- -0.656; label("CF1 capsule effect on F (fractional change)") # Table 3 'CF1 covariate on Fpar' = -0.656
    e_form_tf3_fdepot <- 0.154; label("TF3 tablet effect on F (fractional change)") # Table 3 'TF3 covariate on Fpar' = 0.154
    e_wt_fdepot <- -0.475; label("Power exponent of body weight on F (unitless)") # Table 3 'Body weight at baseline covariate on Fpar' = -0.475
    e_hepimp_fdepot <- -0.0729; label("NCI-ODG class > 0 effect on F (fractional change)") # Table 3 'NCI-ODG class > 0 covariate on Fpar' = -0.0729
    e_fasted_ka <- -0.561; label("Fasting effect on ka (fractional change)") # Table 3 'Fasting state covariate on ka' = -0.561
    e_form_cf1_ka <- -0.442; label("CF1 capsule effect on ka (fractional change)") # Table 3 'CF1 covariate on ka' = -0.442
    e_form_tf1_ka <- 0.305; label("TF1 tablet effect on ka (fractional change)") # Table 3 'TF1 covariate on ka' = 0.305
    e_form_tf1fine_ka <- 0.674; label("TF1* (finely micronised) tablet effect on ka (fractional change)") # Table 3 'TF1* covariate on ka' = 0.674
    e_crcl_cl <- 0.199; label("Power exponent of eGFR on CL/F (unitless)") # Table 3 'eGFR at baseline covariate on CLpar/F' = 0.199
    e_tumtp_hcc_cl <- 0.130; label("Hepatocellular carcinoma effect on CL/F (fractional change)") # Table 3 'Hepatocellular carcinoma covariate on CLpar/F' = 0.130
    e_tumtp_crc_cl <- -0.281; label("Colorectal cancer effect on CL/F (fractional change)") # Table 3 'Colorectal cancer covariate on CLpar/F' = -0.281
    e_conmed_opioid_cl <- -0.167; label("Concomitant mu-opioid effect on CL/F (fractional change)") # Table 3 'mu-Opioids covariate on CLpar/F' = -0.167
    e_study_0028_cl <- -0.115; label("Study MS200095-0028 effect on CL/F (fractional change)") # Table 3 'Study MS200095-0028 covariate on CLpar/F' = -0.115
    e_inr_q <- 3.81; label("Power exponent of baseline INR on Q/F (unitless)") # Table 3 'INR at baseline covariate on Qpar/F' = 3.81
    e_alb_q <- 4.14; label("Power exponent of baseline serum albumin on Q/F (unitless)") # Table 3 'Serum albumin at baseline covariate on Qpar/F' = 4.14
    e_age_vc <- 0.219; label("Power exponent of age on Vc/F (unitless)") # Table 3 'Age covariate on Vc,par/F' = 0.219
    e_tumtp_nsclc_vc <- -0.232; label("NSCLC effect on Vc/F (fractional change)") # Table 3 'Non-small cell lung cancer covariate on Vc,par/F' = -0.232
    e_dis_healthy_vp <- -0.810; label("Healthy-participant effect on Vp/F (fractional change)") # Table 3 'Patient/participant covariate on Vp,par/F' = -0.810; ESM 13 applies it to healthy volunteers

    # ---- MSC2571109A structural parameters (Table 3; footnote a salt factor) ----
    lcl_msc2571109a <- log(40.2); label("MSC2571109A apparent clearance (L/h)") # Table 3 'CLmet' = 40.2 L/h (RSE 2.40%)
    lvc_msc2571109a <- log(131); label("MSC2571109A apparent central volume (L)") # Table 3 'Vc,met' = 131 L (RSE 5.04%)
    lq_msc2571109a <- log(106); label("MSC2571109A apparent inter-compartmental clearance (L/h)") # Table 3 'Qmet' = 106 L/h (RSE 5.98%)
    lvp_msc2571109a <- log(152); label("MSC2571109A apparent peripheral volume (L)") # Table 3 'Vp,met' = 152 L (RSE 2.90%)
    lfm <- fixed(log(1)); label("Fraction of tepotinib clearance forming MSC2571109A (fraction)") # Results 'this parameter was fixed to 1'

    # ---- MSC2571109A covariate effects ----
    e_crcl_cl_msc2571109a <- 0.311; label("Power exponent of eGFR on MSC2571109A CL (unitless)") # Table 3 'eGFR at baseline covariate on CLmet' = 0.311
    e_wt_cl_msc2571109a <- -0.696; label("Power exponent of body weight on MSC2571109A CL (unitless)") # Table 3 'Body weight at baseline covariate on CLmet' = -0.696
    e_tumtp_nsclc_cl_msc2571109a <- 0.498; label("NSCLC effect on MSC2571109A CL (fractional change)") # Table 3 'Non-small cell lung cancer covariate on CLmet' = 0.498
    e_tumtp_hcc_fm <- -0.398; label("Hepatocellular carcinoma effect on fraction metabolised (fractional change)") # Table 3 'Hepatocellular carcinoma covariate on FM' = -0.398
    e_race_asian_northeast_q_msc2571109a <- 1.40; label("East Asian effect on MSC2571109A Q (fractional change)") # Table 3 'East Asian on Qmet' = 1.40
    e_dis_healthy_vp_msc2571109a <- 2.31; label("Healthy-participant effect on MSC2571109A Vp (fractional change)") # Table 3 'Patient/participant covariate on Vp,met' = 2.31
    e_hepimp_vc_msc2571109a <- 0.520; label("NCI-ODG class > 0 effect on MSC2571109A Vc (fractional change)") # Table 3 'NCI-ODG class > 0 covariate on Vc.met' = 0.520

    # ---- IIV. Methods: lognormal IIV 'with standard deviation omega'; the
    #      Table 3 '(CV)' column is read as omega, so variance = value^2.
    #      Subgroup-specific variances select one eta per subject (model()). ----
    etalcl_pt ~ 0.112225 # Table 3 'IIV CLpar (CV)' = 0.335 -> 0.335^2 (patients)
    etalcl_hv ~ 0.016384 # Table 3 'IIV CLpar for healthy participant (CV)' = 0.128 -> 0.128^2
    etalka ~ 0.426409 # Table 3 'IIV ka (CV)' = 0.653 -> 0.653^2
    etald1 ~ 0.425104 # Table 3 'IIV D1 (CV)' = 0.652 -> 0.652^2
    etalfdepot_pt ~ 0.080089 # Table 3 'IIV Fpar (CV)' = 0.283 -> 0.283^2 (patients, non-CF1 doses)
    etalfdepot_cf1 ~ 0.508369 # Table 3 'IIV Fpar for CF1 (CV)' = 0.713 -> 0.713^2
    etalfdepot_hv ~ 0.035344 # Table 3 'IIV Fpar for healthy participant (CV)' = 0.188 -> 0.188^2
    etalcl_msc2571109a_pt ~ 0.287296 # Table 3 'IIV CLmet (CV)' = 0.536 -> 0.536^2 (patients)
    etalcl_msc2571109a_hv ~ 0.065025 # Table 3 'IIV CLmet for healthy participants (CV)' = 0.255 -> 0.255^2
    etalvc_msc2571109a ~ 0.737881 # Table 3 'IIV Vc,met (CV)' = 0.859 -> 0.859^2
    etalq_msc2571109a ~ 0.625681 # Table 3 'IIV Qmet (CV)' = 0.791 -> 0.791^2
    etalvp_msc2571109a ~ 0.061504 # Table 3 'IIV Vp,met (CV)' = 0.248 -> 0.248^2
    etalfm ~ fixed(0) # Results: IIV on FM was estimated but its value is not reported in Table 3

    # ---- Residual error (additive on the log scale = proportional) ----
    propSd <- 0.337; label("Proportional residual error, tepotinib (fraction)") # Table 3 'Prop. RUV (CV)' = 0.337 (tepotinib); Discussion 'residual variability was 33.7%'
    propSd_msc2571109a <- 0.298; label("Proportional residual error, MSC2571109A (fraction)") # Table 3 (continued) 'Pro. RUV (CV)' = 0.298 (MSC2571109A)
  })

  model({
    # 1. Covariate terms (ESM 13: power models normalised to the reference
    #    participant; categorical effects as 1 + theta when the category applies)
    fasted <- 1 - FED
    cov_cl <- (CRCL / 97.28)^e_crcl_cl *
      (1 + e_tumtp_hcc_cl * TUMTP_HCC) *
      (1 + e_tumtp_crc_cl * TUMTP_CRC) *
      (1 + e_conmed_opioid_cl * CONMED_OPIOID) *
      (1 + e_study_0028_cl * STUDY_MS200095_0028)
    cov_vc <- (AGE / 59)^e_age_vc * (1 + e_tumtp_nsclc_vc * TUMTP_NSCLC)
    cov_vp <- 1 + e_dis_healthy_vp * DIS_HEALTHY
    cov_q <- (INR_BASE / 1.06)^e_inr_q * (ALB / 40)^e_alb_q
    cov_f <- (1 + e_dose_fdepot / 100 * (DOSE - 500)) *
      (WT / 72)^e_wt_fdepot *
      (1 + e_fasted_fdepot * fasted) *
      (1 + e_fed_highfat_fdepot * FED_HIGHFAT) *
      (1 + e_form_cf1_fdepot * FORM_TEPOTINIB_CF1) *
      (1 + e_form_tf3_fdepot * FORM_TEPOTINIB_TF3) *
      (1 + e_hepimp_fdepot * HEPIMP)
    cov_d1 <- (1 + e_fasted_d1 * fasted) * (1 + e_hepimp_d1 * HEPIMP)
    cov_ka <- (1 + e_fasted_ka * fasted) *
      (1 + e_form_cf1_ka * FORM_TEPOTINIB_CF1) *
      (1 + e_form_tf1_ka * FORM_TEPOTINIB_TF1) *
      (1 + e_form_tf1fine_ka * FORM_TEPOTINIB_TF1FINE)

    # 2. Individual tepotinib parameters. Clearance and bioavailability carry
    #    subgroup-specific IIV (Table 3): healthy participants draw the _hv eta;
    #    patients draw the CF1 eta on CF1 doses and the _pt eta otherwise.
    cl_pt <- exp(lcl + etalcl_pt)
    cl_hv <- exp(lcl + etalcl_hv)
    cl <- (cl_hv * DIS_HEALTHY + cl_pt * (1 - DIS_HEALTHY)) * cov_cl
    vc <- exp(lvc) * cov_vc
    q <- exp(lq) * cov_q
    vp <- exp(lvp) * cov_vp
    ka <- exp(lka + etalka) * cov_ka
    d1 <- exp(ld1 + etald1) * cov_d1
    f_pt <- exp(lfdepot + etalfdepot_pt)
    f_cf1 <- exp(lfdepot + etalfdepot_cf1)
    f_hv <- exp(lfdepot + etalfdepot_hv)
    f_pat <- f_cf1 * FORM_TEPOTINIB_CF1 + f_pt * (1 - FORM_TEPOTINIB_CF1)
    fdepot <- (f_hv * DIS_HEALTHY + f_pat * (1 - DIS_HEALTHY)) * cov_f

    # 3. Individual MSC2571109A parameters
    cl_msc_pt <- exp(lcl_msc2571109a + etalcl_msc2571109a_pt)
    cl_msc_hv <- exp(lcl_msc2571109a + etalcl_msc2571109a_hv)
    cl_msc2571109a <- (cl_msc_hv * DIS_HEALTHY + cl_msc_pt * (1 - DIS_HEALTHY)) *
      (CRCL / 97.28)^e_crcl_cl_msc2571109a *
      (WT / 72)^e_wt_cl_msc2571109a *
      (1 + e_tumtp_nsclc_cl_msc2571109a * TUMTP_NSCLC)
    vc_msc2571109a <- exp(lvc_msc2571109a + etalvc_msc2571109a) *
      (1 + e_hepimp_vc_msc2571109a * HEPIMP)
    q_msc2571109a <- exp(lq_msc2571109a + etalq_msc2571109a) *
      (1 + e_race_asian_northeast_q_msc2571109a * RACE_ASIAN_NORTHEAST)
    vp_msc2571109a <- exp(lvp_msc2571109a + etalvp_msc2571109a) *
      (1 + e_dis_healthy_vp_msc2571109a * DIS_HEALTHY)
    fm <- exp(lfm + etalfm) * (1 + e_tumtp_hcc_fm * TUMTP_HCC)

    # 4. ODEs (Figure 1). Dose records go to depot as a zero-order input of
    #    duration d1 (rate = -2), which then empties first-order into central.
    #    The metabolite is formed from parent clearance scaled by fm.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    d/dt(central_msc2571109a) <- fm * (cl / vc) * central -
      (cl_msc2571109a / vc_msc2571109a) * central_msc2571109a -
      (q_msc2571109a / vc_msc2571109a) * central_msc2571109a +
      (q_msc2571109a / vp_msc2571109a) * peripheral1_msc2571109a
    d/dt(peripheral1_msc2571109a) <- (q_msc2571109a / vc_msc2571109a) * central_msc2571109a -
      (q_msc2571109a / vp_msc2571109a) * peripheral1_msc2571109a

    # 5. Absorption: zero-order duration into depot and relative bioavailability
    dur(depot) <- d1
    f(depot) <- fdepot

    # 6. Observations. Amounts in mg (free base; parent-equivalent for the
    #    metabolite, since fm is fixed to 1), volumes in L -> x1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc_msc2571109a <- 1000 * central_msc2571109a / vc_msc2571109a
    Cc ~ prop(propSd)
    Cc_msc2571109a ~ prop(propSd_msc2571109a)
  })
}
