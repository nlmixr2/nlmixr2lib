Terranova_2021_berzosertib <- function() {
  description <- "Two-compartment population PK model (full covariate model) for intravenous berzosertib (M6620, VX-970), an ATR inhibitor, in adults with advanced solid tumors receiving 1-h infusions alone or combined with gemcitabine, cisplatin, carboplatin, or carboplatin plus paclitaxel. Linear elimination; a full 4x4 IIV block on CL, V1, Q and V2; combined additive plus proportional residual error. Body weight, age, albumin, platelet count and baseline tumor size (power, centred on the cohort medians) and sex, hepatic impairment, renal impairment, ECOG 0 and tumor type (NSCLC, TNBC, CRC; linear) act on CL, V1 and V2."
  reference <- paste(
    "Terranova N, Jansen M, Falk M, Hendriks BS.",
    "Population pharmacokinetics of ATR inhibitor berzosertib in phase I studies",
    "for different cancer types.",
    "Cancer Chemother Pharmacol. 2021;87(2):185-196.",
    "doi:10.1007/s00280-020-04184-z.",
    "Structural and random-effect estimates from Table 3; covariate coefficients",
    "recovered from the numeric ratios printed in the Figure 2 forest plots."
  )
  vignette <- "Terranova_2021_berzosertib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on CL, V1 and V2 centred on the analysis-population median 72.8 kg (Terranova 2021 Table 2; Methods 'Covariate analysis model development' power equation with CON_median = observed median). Baseline only. The paper chose body weight over BSA for the covariate model even though berzosertib is dosed per m^2.",
      source_name = "Body weight"
    ),
    AGE = list(
      description = "Baseline age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on CL, V1 and V2 centred on the median 60 years (Terranova 2021 Table 2).",
      source_name = "Age"
    ),
    ALB = list(
      description = "Baseline serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on CL, V1 and V2 centred on the median 38 g/L (Terranova 2021 Table 2).",
      source_name = "Albumin"
    ),
    PLT = list(
      description = "Baseline platelet count.",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power model on CL, V1 and V2 centred on the median 274 x 10^9/L (Terranova 2021 Table 2).",
      source_name = "Platelets"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (tumor burden).",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported as 'Tumor burden [mm]' in Terranova 2021 Table 2 and 'Tumor size' in Figure 2; the paper does not state the RECIST construct, so the pooled TUMSZ register is used. Power model on CL, V1 and V2 centred on the median 75 mm.",
      source_name = "Tumor burden / Tumor size"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Linear categorical effect (1 + theta * SEXF) on CL, V1 and V2. Figure 2 labels the plotted category 'Sex: Female', and the text states V1 in females is 15.8% lower than in males, so female is the indicator category and the Table 3 typical values refer to a male. The Methods say CAT = 0 for the most common category (female, 60.4%); the vignette shows that the male-reference reading is the one consistent with the base-model typical values (Table S2). Cohort 145 female / 95 male.",
      source_name = "Sex"
    ),
    HEPIMP = list(
      description = "Hepatic impairment indicator, NCI ODWG (1 = mild or worse, 0 = normal).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = normal hepatic function",
      notes = "The paper pooled mild (n = 43) and severe (n = 1) NCI ODWG impairment into one 'Mild/Severe' category (no moderate subjects), which equals the register's any-impairment HEPIMP. Missing values (n = 30) were set to the most common category, normal (Terranova 2021 Results 'Full covariate model').",
      source_name = "Hepatic impairment"
    ),
    RENALIMP = list(
      description = "Renal impairment indicator, FDA guidance categories (1 = mild or moderate, 0 = normal).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = normal renal function",
      notes = "Renal impairment was defined by the U.S. FDA guidance (Terranova 2021 Methods, ref 22). The paper pooled mild (n = 88) and moderate (n = 15) impairment; no subject had severe impairment, so the pool equals the register's any-impairment RENALIMP.",
      source_name = "Renal impairment"
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status >= 1 indicator (1 = ECOG >= 1, 0 = ECOG 0).",
      units = "(binary)",
      type = "binary",
      reference_category = "1 = ECOG >= 1 (the model reference, ECOG 1 with the rare ECOG 2-3 folded in)",
      notes = "The paper tested ECOG 0 against the most frequent ECOG 1 and folded ECOG >= 2 (n = 11) into ECOG 1, so the effect enters as (1 + theta * (1 - ECOG_GE1)).",
      source_name = "ECOG PS"
    ),
    TUMTP_NSCLC = list(
      description = "Non-small cell lung cancer tumor-type indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; all TUMTP_* = 0 is the 'other tumor types' reference",
      notes = "n = 48 (20.0%). Reference tumor type is the most common category 'Other' (n = 59) together with every tumor type without its own indicator (SCLC, prostate, non-TNBC breast, head and neck, ovarian, mesothelioma). At most one TUMTP_* indicator is 1 per subject.",
      source_name = "Tumor type: NSCLC"
    ),
    TUMTP_TNBC = list(
      description = "Triple-negative breast cancer tumor-type indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; all TUMTP_* = 0 is the 'other tumor types' reference",
      notes = "n = 33 (13.8%), mostly Study 001 Part C2 (basaloid TNBC). Distinct from the 11 non-TNBC breast-cancer patients, who are in the reference pool.",
      source_name = "Tumor type: TNBC"
    ),
    TUMTP_CRC = list(
      description = "Colorectal cancer tumor-type indicator.",
      units = "(binary)",
      type = "binary",
      reference_category = "0; all TUMTP_* = 0 is the 'other tumor types' reference",
      notes = "n = 44 (18.3%).",
      source_name = "Tumor type: CRC"
    )
  )

  covariatesDataExcluded <- list(
    RACE = list(
      description = "Race and ethnicity.",
      units = "(categorical)",
      type = "categorical",
      notes = "Explored graphically against the IIV estimates; no trend, so not included in the full covariate model (Terranova 2021 Results 'Full covariate model'). Cohort 93.3% White."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Explored graphically only; not included."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Explored graphically only; not included."
    ),
    TBILI = list(
      description = "Baseline serum bilirubin.",
      units = "umol/L",
      type = "continuous",
      notes = "Explored graphically only; not included."
    ),
    BSA = list(
      description = "Body surface area.",
      units = "m^2",
      type = "continuous",
      notes = "Used to compute the mg dose from the mg/m^2 dose level. A sensitivity run with BSA in place of body weight gave similar results (data not shown in Terranova 2021)."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "berzosertib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "berzosertib",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 240L,
    n_studies = 2L,
    n_observations = 2546L,
    age_range = "26-79 years",
    age_median = "60 years",
    weight_range = "46-150 kg",
    weight_median = "72.8 kg",
    bsa_median = "1.82 m^2 (range 1.4-2.59)",
    sex_female_pct = 60.4,
    race_ethnicity = c(White = 93.3, Black = 1.25, Asian = 2.08, Other = 2.5, Missing = 0.83),
    disease_state = "Advanced solid tumors (and DDR-defective lymphoma in Study 002 Part C): NSCLC 20.0%, CRC 18.3%, TNBC 13.8%, SCLC 6.7%, mesothelioma 5.0%, non-TNBC breast 4.6%, ovarian 3.8%, prostate 2.9%, head and neck 0.4%, other 24.6%.",
    dose_range = "Berzosertib 18-480 mg/m^2 as 1-h IV infusions (11 dose levels), alone (once or twice weekly) or on days 2 and 9 of 21-day cycles after gemcitabine, cisplatin, gemcitabine plus cisplatin, carboplatin, or carboplatin plus paclitaxel.",
    regions = "Two phase I studies: Study 001 (MS201923-0001, NCT02157792; n = 170) and Study 002 (VX13-970-002, EudraCT 2013-005100-34; n = 70); Terranova 2021 Introduction and Online Resource Table S1.",
    notes = "Baseline characteristics from Terranova 2021 Table 2. Renal impairment none 57.1% / mild 36.7% / moderate 6.25%; hepatic impairment none 69.2% / mild 17.9% / severe 0.4% / missing 12.5%; ECOG 0 27.1%, 1 68.3%, 2 2.1%, 3 2.5%. BLQ (< 10 ng/mL) samples were excluded."
  )

  ini({
    # Structural parameters: Terranova 2021 Table 3 (full covariate model).
    # Typical values are for a male with median continuous covariates, 'other'
    # tumor type, ECOG 1, and no hepatic or renal impairment.
    lcl <- log(65)
    label("Clearance CL (L/h)") # Table 3 'Clearance CL [L/h]' = 65 (RSE 5.2%)
    lvc <- log(118)
    label("Central volume V1 (L)") # Table 3 'Central volume V1 [L]' = 118 (RSE 12%)
    lq <- log(295)
    label("Intercompartmental clearance Q (L/h)") # Table 3 'Intercompartmental clearance Q [L/h]' = 295 (RSE 3.5%)
    lvp <- log(1030)
    label("Peripheral volume V2 (L)") # Table 3 'Peripheral volume V2 [L]' = 1030 (RSE 3.9%)

    # Continuous covariates, power form P = TV * (CON / median)^theta.
    # Terranova 2021 prints no covariate coefficients; each exponent is
    # back-solved from the two Figure 2 forest-plot ratios at the 2.5th and
    # 97.5th covariate percentiles as theta = log(r_high / r_low) /
    # log(x_high / x_low), which does not depend on the centring median.
    e_wt_cl <- 0.257
    label("Power exponent of WT/72.8 on CL (unitless)") # Fig 2A: Weight 51 kg 0.911, 114 kg 1.12 -> log(1.12/0.911)/log(114/51)
    e_age_cl <- -0.0929
    label("Power exponent of AGE/60 on CL (unitless)") # Fig 2A: Age 35 y 1.05, 76 y 0.977
    e_alb_cl <- 0.213
    label("Power exponent of ALB/38 on CL (unitless)") # Fig 2A: Albumin 26 g/L 0.921, 46 g/L 1.04
    e_plt_cl <- -0.132
    label("Power exponent of PLT/274 on CL (unitless)") # Fig 2A: Platelets 132 x10^9/L 1.1, 551 x10^9/L 0.911
    e_tumsz_cl <- 0.00122
    label("Power exponent of TUMSZ/75 on CL (unitless)") # Fig 2A: Tumor size 17.0 mm 0.997, 201 mm 1
    e_wt_vc <- -0.128
    label("Power exponent of WT/72.8 on V1 (unitless)") # Fig 2B: Weight 51 kg 1.05, 114 kg 0.947
    e_age_vc <- 0.0402
    label("Power exponent of AGE/60 on V1 (unitless)") # Fig 2B: Age 35 y 0.979, 76 y 1.01
    e_alb_vc <- 0.881
    label("Power exponent of ALB/38 on V1 (unitless)") # Fig 2B: Albumin 26 g/L 0.714, 46 g/L 1.18
    e_plt_vc <- -0.149
    label("Power exponent of PLT/274 on V1 (unitless)") # Fig 2B: Platelets 115 x10^9/L 1.13, 565 x10^9/L 0.891
    e_tumsz_vc <- 0.216
    label("Power exponent of TUMSZ/75 on V1 (unitless)") # Fig 2B: Tumor size 17.0 mm 0.727, 201 mm 1.24
    e_wt_vp <- 0.111
    label("Power exponent of WT/72.8 on V2 (unitless)") # Fig 2C: Weight 51 kg 0.96, 114 kg 1.05
    e_age_vp <- 0.335
    label("Power exponent of AGE/60 on V2 (unitless)") # Fig 2C: Age 35 y 0.833, 76 y 1.08
    e_alb_vp <- 0.259
    label("Power exponent of ALB/38 on V2 (unitless)") # Fig 2C: Albumin 26 g/L 0.906, 46 g/L 1.05
    e_plt_vp <- -0.0477
    label("Power exponent of PLT/274 on V2 (unitless)") # Fig 2C: Platelets 115 x10^9/L 1.04, 565 x10^9/L 0.964
    e_tumsz_vp <- 0.00810
    label("Power exponent of TUMSZ/75 on V2 (unitless)") # Fig 2C: Tumor size 17.0 mm 0.99, 201 mm 1.01

    # Categorical covariates, linear form P = TV * (1 + theta * CAT). The
    # Figure 2 ratio for the indicator category is 1 + theta.
    e_sexf_cl <- -0.121
    label("Fractional change in CL, female vs male (unitless)") # Fig 2A: Sex: Female 0.879
    e_hepimp_cl <- -0.034
    label("Fractional change in CL, mild/severe hepatic impairment (unitless)") # Fig 2A: Hepatic impairment Mild/Severe 0.966
    e_renalimp_cl <- -0.026
    label("Fractional change in CL, mild/moderate renal impairment (unitless)") # Fig 2A: Renal impairment Mild/Moderate 0.974
    e_ecog0_cl <- 0.02
    label("Fractional change in CL, ECOG 0 vs ECOG >= 1 (unitless)") # Fig 2A: ECOG: 0 1.02
    e_tumtp_nsclc_cl <- -0.12
    label("Fractional change in CL, NSCLC vs other tumor types (unitless)") # Fig 2A: Tumor type: NSCLC 0.88
    e_tumtp_tnbc_cl <- 0.42
    label("Fractional change in CL, TNBC vs other tumor types (unitless)") # Fig 2A: Tumor type: TNBC 1.42; text: TNBC had 42% higher CL
    e_tumtp_crc_cl <- -0.057
    label("Fractional change in CL, CRC vs other tumor types (unitless)") # Fig 2A: Tumor type: CRC 0.943
    e_sexf_vc <- -0.158
    label("Fractional change in V1, female vs male (unitless)") # Fig 2B: Sex: Female 0.842; text: V1 15.8% lower in females
    e_hepimp_vc <- 0.17
    label("Fractional change in V1, mild/severe hepatic impairment (unitless)") # Fig 2B: Hepatic impairment Mild/Severe 1.17
    e_renalimp_vc <- 0.12
    label("Fractional change in V1, mild/moderate renal impairment (unitless)") # Fig 2B: Renal impairment Mild/Moderate 1.12
    e_ecog0_vc <- -0.044
    label("Fractional change in V1, ECOG 0 vs ECOG >= 1 (unitless)") # Fig 2B: ECOG: 0 0.956
    e_tumtp_nsclc_vc <- -0.11
    label("Fractional change in V1, NSCLC vs other tumor types (unitless)") # Fig 2B: Tumor type: NSCLC 0.89
    e_tumtp_tnbc_vc <- 0.29
    label("Fractional change in V1, TNBC vs other tumor types (unitless)") # Fig 2B: Tumor type: TNBC 1.29
    e_tumtp_crc_vc <- 0.15
    label("Fractional change in V1, CRC vs other tumor types (unitless)") # Fig 2B: Tumor type: CRC 1.15
    e_sexf_vp <- -0.139
    label("Fractional change in V2, female vs male (unitless)") # Fig 2C: Sex: Female 0.861
    e_hepimp_vp <- -0.039
    label("Fractional change in V2, mild/severe hepatic impairment (unitless)") # Fig 2C: Hepatic impairment Mild/Severe 0.961
    e_renalimp_vp <- -0.042
    label("Fractional change in V2, mild/moderate renal impairment (unitless)") # Fig 2C: Renal impairment Mild/Moderate 0.958
    e_ecog0_vp <- 0
    label("Fractional change in V2, ECOG 0 vs ECOG >= 1 (unitless)") # Fig 2C: ECOG: 0 1 (95% CI 0.943-1.06)
    e_tumtp_nsclc_vp <- -0.109
    label("Fractional change in V2, NSCLC vs other tumor types (unitless)") # Fig 2C: Tumor type: NSCLC 0.891
    e_tumtp_tnbc_vp <- -0.307
    label("Fractional change in V2, TNBC vs other tumor types (unitless)") # Fig 2C: Tumor type: TNBC 0.693
    e_tumtp_crc_vp <- -0.052
    label("Fractional change in V2, CRC vs other tumor types (unitless)") # Fig 2C: Tumor type: CRC 0.948

    # IIV: full 4x4 block, Table 3 (variances and covariances on the log scale).
    # Order CL, V1, Q, V2 as in Table 3.
    etalcl + etalvc + etalq + etalvp ~ c(
      0.066,
      0.060, 0.32,
      0.054, 0.25, 0.24,
      0.041, 0.081, 0.090, 0.047
    ) # Table 3: IIV CL 0.066; cov(CL,V1) 0.060; IIV V1 0.32; cov(CL,Q) 0.054; cov(V1,Q) 0.25; IIV Q 0.24; cov(CL,V2) 0.041; cov(V1,V2) 0.081; cov(Q,V2) 0.090; IIV V2 0.047

    # Residual error: combined additive and proportional (Table 3).
    propSd <- 0.22
    label("Proportional residual error (SD, fraction)") # Table 3 'Proportional residual error [sd]' = 0.22 (RSE 4.6%)
    addSd <- 1.73
    label("Additive residual error (SD, ng/mL)") # Table 3 'Additive residual error [ng/mL]' = 1.73 (RSE 19%)
  })
  model({
    # Individual parameters (Terranova 2021 Methods, covariate equations:
    # power for continuous covariates centred on the Table 2 medians, linear
    # for categorical covariates; effects combine multiplicatively).
    cl <- exp(lcl + etalcl) *
      (WT / 72.8)^e_wt_cl *
      (AGE / 60)^e_age_cl *
      (ALB / 38)^e_alb_cl *
      (PLT / 274)^e_plt_cl *
      (TUMSZ / 75)^e_tumsz_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_hepimp_cl * HEPIMP) *
      (1 + e_renalimp_cl * RENALIMP) *
      (1 + e_ecog0_cl * (1 - ECOG_GE1)) *
      (1 + e_tumtp_nsclc_cl * TUMTP_NSCLC) *
      (1 + e_tumtp_tnbc_cl * TUMTP_TNBC) *
      (1 + e_tumtp_crc_cl * TUMTP_CRC)
    vc <- exp(lvc + etalvc) *
      (WT / 72.8)^e_wt_vc *
      (AGE / 60)^e_age_vc *
      (ALB / 38)^e_alb_vc *
      (PLT / 274)^e_plt_vc *
      (TUMSZ / 75)^e_tumsz_vc *
      (1 + e_sexf_vc * SEXF) *
      (1 + e_hepimp_vc * HEPIMP) *
      (1 + e_renalimp_vc * RENALIMP) *
      (1 + e_ecog0_vc * (1 - ECOG_GE1)) *
      (1 + e_tumtp_nsclc_vc * TUMTP_NSCLC) *
      (1 + e_tumtp_tnbc_vc * TUMTP_TNBC) *
      (1 + e_tumtp_crc_vc * TUMTP_CRC)
    # No covariates were tested on Q (Terranova 2021 Results).
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) *
      (WT / 72.8)^e_wt_vp *
      (AGE / 60)^e_age_vp *
      (ALB / 38)^e_alb_vp *
      (PLT / 274)^e_plt_vp *
      (TUMSZ / 75)^e_tumsz_vp *
      (1 + e_sexf_vp * SEXF) *
      (1 + e_hepimp_vp * HEPIMP) *
      (1 + e_renalimp_vp * RENALIMP) *
      (1 + e_ecog0_vp * (1 - ECOG_GE1)) *
      (1 + e_tumtp_nsclc_vp * TUMTP_NSCLC) *
      (1 + e_tumtp_tnbc_vp * TUMTP_TNBC) *
      (1 + e_tumtp_crc_vp * TUMTP_CRC)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L; x1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
