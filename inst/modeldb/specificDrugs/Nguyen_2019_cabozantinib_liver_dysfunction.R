Nguyen_2019_cabozantinib_liver_dysfunction <- function() {
  description <- "Two-compartment population PK model for oral cabozantinib carrying NCI-ODWG liver-dysfunction covariates on CL/F and Vc/F, in healthy volunteers and patients with hepatocellular carcinoma, renal cell carcinoma, castration-resistant prostate cancer, medullary thyroid cancer, glioblastoma multiforme, or other advanced malignancies (Nguyen 2019, n=2023 across 10 clinical studies). This is Table 3, column 2, the liver-dysfunction covariate model that motivated the analysis: it adds four parameters, for mild and for pooled moderate-or-severe hepatic impairment per the National Cancer Institute Organ Dysfunction Working Group criteria, on top of the final updated integrated model. The authors retained the model WITHOUT these covariates as their final model because the four extra parameters dropped the objective function by only about 10 units and changed the other parameter estimates by less than 15 percent; this file exists so the published liver-dysfunction effects are usable, and the companion Nguyen_2019_cabozantinib carries the final model. Mild hepatic impairment raises CL/F by 12 percent and pooled moderate-or-severe hepatic impairment lowers it by 2 percent, both clinically negligible, supporting the conclusion that no initial dosage adjustment is needed for cancer patients with mild liver dysfunction. Structure is otherwise identical to the final model: parallel first-order (fraction F1 into depot1, rate Ka, lag ALAG1) plus zero-order (duration D2) absorption, a dose power effect and a capsule formulation effect on Ka, a capsule effect on relative oral bioavailability, and two-compartment disposition with first-order elimination."
  reference <- paste(
    "Nguyen L, Chapel S, Tran BD, Lacy S.",
    "Updated population pharmacokinetic model of cabozantinib integrating",
    "various cancer types including hepatocellular carcinoma.",
    "J Clin Pharmacol. 2019;59(11):1551-1561.",
    "doi:10.1002/jcph.1467.",
    "Parameter values from Table 3, column 2",
    "('Including HCC Population and Liver Dysfunction Covariates').",
    "The companion final model is modellib('Nguyen_2019_cabozantinib').",
    sep = " "
  )
  vignette <- "Nguyen_2019_cabozantinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot1 = list(analyte = "cabozantinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cabozantinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Subject age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Power effect on CL/F (exponent -0.16) and Vc/F (exponent 0.077), centered at the updated-cohort median. Nguyen 2019 Methods 'Prior Integrated PPK Model': 'The relationship between continuous covariates and typical value of PK parameters was modeled using the power function with centering by median values.' Median age 64 years is the All Studies column of Nguyen 2019 Table 2.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Power effect on CL/F (exponent -0.0209, near-null) and Vc/F (exponent 1.2, near-linear scaling), centered at the updated-cohort median. Median weight 78 kg is the All Studies column of Nguyen 2019 Table 2.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male, the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.762) and Vc/F (1.08). Nguyen 2019 Table 3 column 2.",
      source_name = "SEXF"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black; White is the typical-value reference category when paired with RACE_ASIAN = 0 and RACE_OTHER = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.18) and Vc/F (1.07). Nguyen 2019 Table 3 column 2. Reference = White (77 percent of the pooled cohort).",
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian; White is the typical-value reference when paired with RACE_BLACK = 0 and RACE_OTHER = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.934) and Vc/F (0.739). Nguyen 2019 Table 3 column 2.",
      source_name = "RACE"
    ),
    RACE_OTHER = list(
      description = "Race category 'Other' indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Other; White is the typical-value reference when paired with RACE_BLACK = 0 and RACE_ASIAN = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.02, near-null) and Vc/F (0.965). Nguyen 2019 Table 3 column 2.",
      source_name = "RACE"
    ),
    TUMTP_HCC = list(
      description = "Hepatocellular carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-HCC; healthy volunteer is the typical-value reference when paired with all other TUMTP_* indicators = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.82) and Vc/F (0.81). Nguyen 2019 Table 3 column 2, 'HCC covariates' block. The HCC effect on CL/F is somewhat larger in this model (0.82) than in the final model without liver-dysfunction covariates (0.878), because the mild-impairment covariate (1.12) absorbs part of the HCC stratum: 65-68 percent of the HCC patients had mild liver dysfunction per Table 2. This covariate is distinct from the HEPIMP_* family, which classifies hepatic function by bilirubin and AST regardless of tumor type.",
      source_name = "POP"
    ),
    TUMTP_RCC = list(
      description = "Renal cell carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-RCC; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.862) and Vc/F (0.711). Nguyen 2019 Table 3 column 2.",
      source_name = "POP"
    ),
    TUMTP_HRPC = list(
      description = "Castration-resistant prostate cancer indicator (paper writes CRPC; canonical column TUMTP_HRPC covers both HRPC and CRPC wordings)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-CRPC; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.968) and Vc/F (0.721). Nguyen 2019 Table 3 column 2.",
      source_name = "POP"
    ),
    TUMTP_MTC = list(
      description = "Medullary thyroid carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-MTC; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.88, i.e. 88 percent higher than healthy volunteers) and Vc/F (0.912). Nguyen 2019 Table 3 column 2. MTC remains the only tumor type with an appreciable PK difference from healthy volunteers.",
      source_name = "POP"
    ),
    TUMTP_GLIO = list(
      description = "Glioblastoma multiforme indicator (canonical TUMTP_GLIO covers glioma of any grade including GB)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-glioma; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.2) and Vc/F (0.448). Nguyen 2019 Table 3 column 2.",
      source_name = "POP"
    ),
    TUMTP_OTHER = list(
      description = "Heterogeneous 'other malignancy' pool indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-other; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.14) and Vc/F (0.749). Nguyen 2019 Table 3 column 2.",
      source_name = "POP"
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator per NCI-ODWG criteria",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function, when paired with HEPIMP_MODSEV = 0)",
      notes = "Time-fixed baseline classification. Multiplicative effect on CL/F (1.12, i.e. 12 percent higher clearance) and Vc/F (1.04). Nguyen 2019 Table 3 column 2, 'Liver dysfunction (NCI-ODWG) covariates' block. Classification scheme is NCI-ODWG, defined in Methods 'Analysis of Data Files' from total bilirubin (TB) and aspartate aminotransferase (AST): mild = TB <= ULN and AST > ULN, or TB > 1-3 x ULN with any AST. Mutually exclusive with HEPIMP_MODSEV; both zero gives the normal-hepatic-function reference. Mild impairment covered 558 / 2023 = 28 percent of the cohort (Table 2). Results: 'Patients with mild or moderate/severe liver dysfunction were predicted to have minimal differences (12% or less) in CL/F and Vc/F relative to subjects with normal liver function.'",
      source_name = "NCI-ODWG liver dysfunction group"
    ),
    HEPIMP_MODSEV = list(
      description = "Composite moderate-or-severe hepatic impairment indicator per NCI-ODWG criteria",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function or mild impairment; paired with HEPIMP_MILD to identify the normal-function reference)",
      notes = "Time-fixed baseline classification. Multiplicative effect on CL/F (0.978) and Vc/F (1.06), both indistinguishable from no effect. Nguyen 2019 Table 3 column 2. Classification scheme is NCI-ODWG: moderate = TB > 1.5-3 x ULN with any AST; severe = TB > 3 x ULN with any AST. The paper pools the two because of sparse data: only 15 moderate and 1 severe subject out of 2023 (Table 2). Results: 'the effect of liver dysfunction per NCI-ODWG criteria on CL/F and Vc/F was assessed using the categorical covariates of mild and combined moderate/severe liver dysfunction per NCI-ODWG criteria due to limited data for these latter 2 groups.' The very wide 90 percent CI on Vc/F (0.658, 1.71) reflects that sample size; Conclusions accordingly state that limited data 'preclude providing any dosing recommendations for these subpopulations'.",
      source_name = "NCI-ODWG liver dysfunction group"
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet; the typical-value reference)",
      notes = "Per-dose-occasion indicator. Multiplicative effect on Ka (0.528) and on overall relative oral bioavailability (0.841). Nguyen 2019 Table 3 column 2. Reference = tablet (Cabometyx). Capsule data were 648 / 2023 = 32 percent of subjects.",
      source_name = "FORM"
    ),
    DOSE = list(
      description = "Administered cabozantinib dose level (free base equivalent)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-dose-occasion. Power covariate on the first-order absorption rate constant: Ka(DOSE) = Ka_ref * (DOSE / 60 mg)^0.564. Nguyen 2019 Table 3 column 2 row 'Dose-dependent Ka'. The reference dose is not printed alongside the exponent, but Methods 'Covariate Effects' fixes the model's reference condition as a subject 'receiving a 60-mg free base equivalent cabozantinib tablet dose once daily', so 60 mg is the paper's own reference dose.",
      source_name = "DOSE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2023L,
    n_observations = 9510L,
    n_studies = 10L,
    age_range = "18-87 years",
    age_median = "64 years (All Studies column, Nguyen 2019 Table 2)",
    weight_range = "30.4-190.7 kg",
    weight_median = "78 kg (All Studies column, Nguyen 2019 Table 2)",
    sex_female_pct = 15.7,
    race_ethnicity = "White 77%, Asian 10%, Black 3%, Other 2%, unknown 8% (Nguyen 2019 Table 2, All Studies)",
    disease_state = "Pooled population: healthy volunteers (7%); castration-resistant prostate cancer (41%); hepatocellular carcinoma (24%); renal cell carcinoma (14%); medullary thyroid cancer (10%); glioblastoma multiforme (2%); other malignancies (2%).",
    hepatic_function = "Per NCI-ODWG criteria: normal 1425 (70%); mild 558 (28%); moderate 15 (1%); severe 1 (<1%); missing 24 (1%). This stratification is the covariate that distinguishes this model from Nguyen_2019_cabozantinib. Of the hepatocellular carcinoma patients, 99% were Child-Pugh A.",
    dose_range = "Oral cabozantinib free base equivalent; 100 mg once daily capsule in the phase 2 RDT and 60 mg once daily tablet in CELESTIAL, pooled with the wider 20-200 mg/day range of the earlier integrated analysis.",
    regions = "Multinational; the two hepatocellular carcinoma studies enrolled 33-35% Asian subjects (Nguyen 2019 Table 2)",
    formulations = "Capsule 648 subjects (32%) and tablet 1375 subjects (68%) (Nguyen 2019 Table 2, All Studies).",
    notes = "Baseline demographics from Nguyen 2019 Table 2; liver-function strata from the same table's 'Liver Dysfunction' block. Bioanalysis by validated LC-MS/MS with a 0.5 ng/mL lower limit of quantification. This model was NOT the authors' final model: Results state that the four liver-dysfunction parameters 'dropped the objective function by ~10 units, and the difference in parameter estimates was <15% with and without liver dysfunction covariates', and that 'The initial model including the hepatocellular carcinoma population was considered the final updated integrated PPK model.'"
  )

  ini({
    # ---- Structural population parameters (Nguyen 2019 Table 3, column 2
    # 'Including HCC Population and Liver Dysfunction Covariates') ----
    #
    # Table 3 reports TRANSFORMED estimates throughout: the unnumbered footnote
    # states 'Transformed estimate is a PK parameter obtained by exponentiating
    # the original estimate.' The structural parameters below are therefore the
    # already-exponentiated values, and the log() here re-creates the estimated
    # THETA scale.
    #
    # Reference covariate set: tablet formulation, 60 mg dose, healthy
    # volunteer, male, White, normal hepatic function, age 64 y, weight 78 kg.
    lka <- log(1.23)
    label("First-order absorption rate constant from depot1 at the 60 mg tablet reference (1/h)") # Nguyen 2019 Table 3 col 2 Ka = 1.23 (90% CI 0.833, 1.82)
    ld2 <- log(2.53)
    label("Duration of the parallel zero-order absorption process (h)") # Nguyen 2019 Table 3 col 2 'Duration for zero-order absorption process' = 2.53 (90% CI 2.25, 2.84)
    lcl <- log(2.47)
    label("Apparent oral clearance at the reference covariate set (L/h)") # Nguyen 2019 Table 3 col 2 CL/F = 2.47 (90% CI 2.26, 2.7)
    lvc <- log(214)
    label("Apparent central volume of distribution at the reference covariate set (L)") # Nguyen 2019 Table 3 col 2 Vc/F = 214 (90% CI 181, 251)
    lq <- log(30.2)
    label("Apparent inter-compartmental clearance (L/h)") # Nguyen 2019 Table 3 col 2 Q/F = 30.2 (90% CI 27.6, 33.1)
    lvp <- log(179)
    label("Apparent peripheral volume of distribution (L)") # Nguyen 2019 Table 3 col 2 Vp/F = 179 (90% CI 167, 191)
    lalag1 <- log(0.82)
    label("Absorption lag time on the depot1 first-order process (h)") # Nguyen 2019 Table 3 col 2 ALAG1 = 0.82 (90% CI 0.795, 0.846)

    # F1 is the fraction of the dose routed to the first-order depot1; the
    # remaining (1 - F1) is delivered by the zero-order process. Table 3
    # footnote b: 'Anti-logit transformation was used to obtain F1', so the
    # estimated THETA is on the logit scale and 0.83 is the back-transformed
    # proportion. logit(0.83) = log(0.83 / 0.17) = 1.5856.
    logitffo <- 1.5856
    label("Logit of the fraction of dose absorbed via the first-order depot1 process (F1 = expit(logitffo))") # Nguyen 2019 Table 3 col 2 'Fraction of dose in first absorption depot F1' = 0.83 (90% CI 0.80, 0.87); logit(0.83) = 1.5856

    # Overall bioavailability at the tablet reference is fixed at 1: the paper
    # estimates only the capsule-relative availability covariate.
    lfdepot <- fixed(log(1.0))
    label("Overall relative oral bioavailability at the tablet reference (unitless)") # Nguyen 2019 Table 3 col 2 reports only 'Capsule on overall relative oral availability'; tablet is the reference formulation

    # ---- Dose-dependent Ka power exponent ----
    e_dose_ka <- 0.564
    label("Power exponent for dose on Ka, applied as (DOSE / 60 mg)^e_dose_ka (unitless)") # Nguyen 2019 Table 3 col 2 'Dose-dependent Ka' = 0.564 (90% CI 0.138, 0.989)

    # ---- Categorical covariate effects ----
    # Table 3 footnote c: 'For categorical covariates (eg, capsule),
    # transformed estimates correspond to multiplicative change from the
    # typical PK parameter.' Each value below is a direct MULTIPLIER applied as
    # factor^indicator, so an indicator of 0 contributes 1.

    # Formulation effects (capsule vs tablet reference)
    e_form_capsule_ka <- 0.528
    label("Capsule (vs tablet) multiplicative factor on Ka (unitless)") # Nguyen 2019 Table 3 col 2 'Capsule on Ka' = 0.528 (90% CI 0.28, 0.994)
    e_form_capsule_f <- 0.841
    label("Capsule (vs tablet) multiplicative factor on overall relative oral bioavailability (unitless)") # Nguyen 2019 Table 3 col 2 'Capsule on overall relative oral availability' = 0.841 (90% CI 0.824, 0.859)

    # Sex effects (female vs male reference)
    e_sexf_cl <- 0.762
    label("Female (vs male) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Female on CL/F' = 0.762 (90% CI 0.715, 0.811)
    e_sexf_vc <- 1.08
    label("Female (vs male) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Female on Vc/F' = 1.08 (90% CI 0.952, 1.22)

    # Race effects (Black / Asian / Other vs White reference)
    e_race_black_cl <- 1.18
    label("Race Black (vs White) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Black on CL/F' = 1.18 (90% CI 1.04, 1.33)
    e_race_black_vc <- 1.07
    label("Race Black (vs White) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Black on Vc/F' = 1.07 (90% CI 0.789, 1.44)
    e_race_asian_cl <- 0.934
    label("Race Asian (vs White) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Asian on CL/F' = 0.934 (90% CI 0.868, 1)
    e_race_asian_vc <- 0.739
    label("Race Asian (vs White) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Asian on Vc/F' = 0.739 (90% CI 0.595, 0.918)
    e_race_other_cl <- 1.02
    label("Race Other (vs White) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Other race on CL/F' = 1.02 (90% CI 0.9, 1.16)
    e_race_other_vc <- 0.965
    label("Race Other (vs White) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Other Race on Vc/F' = 0.965 (90% CI 0.681, 1.37)

    # Tumor-type effects (HCC / RCC / CRPC / MTC / GB / Other vs healthy-volunteer reference)
    e_hcc_cl <- 0.82
    label("Hepatocellular carcinoma (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'HCC on CL/F' = 0.82 (90% CI 0.738, 0.912)
    e_hcc_vc <- 0.81
    label("Hepatocellular carcinoma (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'HCC on Vc/F' = 0.81 (90% CI 0.652, 1.01)
    e_rcc_cl <- 0.862
    label("RCC (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'RCC on CL/F' = 0.862 (90% CI 0.778, 0.956)
    e_rcc_vc <- 0.711
    label("RCC (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'RCC on Vc/F' = 0.711 (90% CI 0.454, 1.11)
    e_hrpc_cl <- 0.968
    label("CRPC (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'CRPC on CL/F' = 0.968 (90% CI 0.874, 1.07)
    e_hrpc_vc <- 0.721
    label("CRPC (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'CRPC on Vc/F' = 0.721 (90% CI 0.583, 0.891)
    e_mtc_cl <- 1.88
    label("MTC (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'MTC on CL/F' = 1.88 (90% CI 1.69, 2.08)
    e_mtc_vc <- 0.912
    label("MTC (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'MTC on Vc/F' = 0.912 (90% CI 0.769, 1.08)
    e_glio_cl <- 1.2
    label("Glioblastoma multiforme (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'GB on CL/F' = 1.2 (90% CI 1, 1.44)
    e_glio_vc <- 0.448
    label("Glioblastoma multiforme (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'GB on Vc/F' = 0.448 (90% CI 0.304, 0.659)
    e_oth_cl <- 1.14
    label("Other malignancies (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Other malignancies on CL/F' = 1.14 (90% CI 0.971, 1.35)
    e_oth_vc <- 0.749
    label("Other malignancies (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Other malignancies on Vc/F' = 0.749 (90% CI 0.583, 0.962)

    # ---- NCI-ODWG liver-dysfunction effects ----
    # The four parameters that define this model. HEPIMP_MILD and
    # HEPIMP_MODSEV are mutually exclusive; both zero is the
    # normal-hepatic-function reference.
    e_hepimp_mild_cl <- 1.12
    label("Mild hepatic impairment (NCI-ODWG, vs normal) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Mild liver dysfunction on CL/F' = 1.12 (90% CI 1.06, 1.18)
    e_hepimp_mild_vc <- 1.04
    label("Mild hepatic impairment (NCI-ODWG, vs normal) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Mild liver dysfunction on Vc/F' = 1.04 (90% CI 0.904, 1.2)
    e_hepimp_modsev_cl <- 0.978
    label("Moderate or severe hepatic impairment (NCI-ODWG, vs normal) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 2 'Moderate and severe liver dysfunction on CL/F' = 0.978 (90% CI 0.781, 1.22)
    e_hepimp_modsev_vc <- 1.06
    label("Moderate or severe hepatic impairment (NCI-ODWG, vs normal) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 2 'Moderate and severe liver dysfunction on Vc/F' = 1.06 (90% CI 0.658, 1.71)

    # ---- Continuous covariate effects (power exponents, untransformed) ----
    e_age_cl <- -0.16
    label("Power exponent for age on CL/F, applied as (AGE / 64 y)^e_age_cl (unitless)") # Nguyen 2019 Table 3 col 2 'Age on CL/F' = -0.16 (90% CI -0.266, -0.0539)
    e_age_vc <- 0.077
    label("Power exponent for age on Vc/F, applied as (AGE / 64 y)^e_age_vc (unitless)") # Nguyen 2019 Table 3 col 2 'Age on Vc/F' = 0.077 (90% CI -0.136, 0.29)
    e_wt_cl <- -0.0209
    label("Power exponent for body weight on CL/F, applied as (WT / 78 kg)^e_wt_cl (unitless)") # Nguyen 2019 Table 3 col 2 'Weight on CL/F' = -0.0209 (90% CI -0.128, 0.0863)
    e_wt_vc <- 1.2
    label("Power exponent for body weight on Vc/F, applied as (WT / 78 kg)^e_wt_vc (unitless)") # Nguyen 2019 Table 3 col 2 'Weight on Vc/F' = 1.2 (90% CI 0.934, 1.46)

    # ---- Inter-individual variability ----
    # Nguyen 2019 Table 3 col 2 'Variance' block. The table footnote defines
    # omega^2 as the 'variance of population parameter modeled using
    # exponential model', so these are log-scale variances and the
    # colon-separated CL/F:Vc/F row is the OMEGA BLOCK off-diagonal covariance.
    # Cauchy-Schwarz holds: sqrt(0.210 * 0.430) = 0.3005 exceeds the reported
    # covariance 0.199, giving a valid correlation of 0.66.
    etalcl + etalvc ~ c(0.210, 0.199, 0.430) # Nguyen 2019 Table 3 col 2: omega^2 CL/F = 0.210; omega^2 CL/F:Vc/F = 0.199; omega^2 Vc/F = 0.430
    etalka ~ 2.21 # Nguyen 2019 Table 3 col 2 omega^2 Ka = 2.21 (90% CI 1.67, 2.75)
    etalogitffo ~ 2.73 # Nguyen 2019 Table 3 col 2 omega^2 F1 = 2.73 (90% CI 2.11, 3.36), on the logit scale per footnote b

    # ---- Residual error ----
    # Methods: 'The residual variability was modeled using the log-transformed
    # additive-error model', which is equivalent to a proportional residual in
    # the linear concentration space nlmixr2 works in. sqrt(0.127) = 0.35637.
    propSd <- 0.35637
    label("Proportional residual error (fraction; sqrt of the log-scale additive variance)") # Nguyen 2019 Table 3 col 2 sigma^2 = 0.127 (90% CI 0.123, 0.131)
  })

  model({
    # ---- Reference (centering) covariate values ----
    ref_age <- 64 # years; Nguyen 2019 Table 2 All Studies median age
    ref_wt <- 78 # kg; Nguyen 2019 Table 2 All Studies median body weight
    ref_dose <- 60 # mg; Nguyen 2019 Methods 'Covariate Effects' reference condition

    # ---- Categorical covariate multipliers ----
    # Applied as factor^indicator so that an indicator of 0 contributes a
    # multiplier of exactly 1 (Table 3 footnote c).
    cl_sex <- e_sexf_cl^SEXF
    vc_sex <- e_sexf_vc^SEXF

    cl_race <- e_race_black_cl^RACE_BLACK *
      e_race_asian_cl^RACE_ASIAN *
      e_race_other_cl^RACE_OTHER
    vc_race <- e_race_black_vc^RACE_BLACK *
      e_race_asian_vc^RACE_ASIAN *
      e_race_other_vc^RACE_OTHER

    cl_tumtp <- e_hcc_cl^TUMTP_HCC *
      e_rcc_cl^TUMTP_RCC *
      e_hrpc_cl^TUMTP_HRPC *
      e_mtc_cl^TUMTP_MTC *
      e_glio_cl^TUMTP_GLIO *
      e_oth_cl^TUMTP_OTHER
    vc_tumtp <- e_hcc_vc^TUMTP_HCC *
      e_rcc_vc^TUMTP_RCC *
      e_hrpc_vc^TUMTP_HRPC *
      e_mtc_vc^TUMTP_MTC *
      e_glio_vc^TUMTP_GLIO *
      e_oth_vc^TUMTP_OTHER

    # Liver dysfunction per NCI-ODWG; both indicators zero gives the
    # normal-hepatic-function reference multiplier of 1.
    cl_hepimp <- e_hepimp_mild_cl^HEPIMP_MILD *
      e_hepimp_modsev_cl^HEPIMP_MODSEV
    vc_hepimp <- e_hepimp_mild_vc^HEPIMP_MILD *
      e_hepimp_modsev_vc^HEPIMP_MODSEV

    ka_form <- e_form_capsule_ka^FORM_CAPSULE
    f_form <- e_form_capsule_f^FORM_CAPSULE

    # ---- Individual PK parameters ----
    ka <- exp(lka + etalka) *
      (DOSE / ref_dose)^e_dose_ka *
      ka_form
    cl <- exp(lcl + etalcl) *
      (AGE / ref_age)^e_age_cl *
      (WT / ref_wt)^e_wt_cl *
      cl_sex * cl_race * cl_tumtp * cl_hepimp
    vc <- exp(lvc + etalvc) *
      (AGE / ref_age)^e_age_vc *
      (WT / ref_wt)^e_wt_vc *
      vc_sex * vc_race * vc_tumtp * vc_hepimp
    q <- exp(lq)
    vp <- exp(lvp)

    d2 <- exp(ld2)
    alag1 <- exp(lalag1)
    f1 <- expit(logitffo + etalogitffo)
    fdepot <- exp(lfdepot) * f_form

    # ---- Micro-constants ----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system: parallel dual absorption + 2-compartment disposition ----
    # A fraction F1 of the dose enters depot1 and is absorbed first-order at
    # rate ka after the lag alag1; the remaining (1 - F1) is delivered
    # directly into central by a zero-order process of duration d2. Dosing
    # events must target BOTH compartments with the same nominal amount:
    # cmt = 'depot1' as an ordinary oral dose, and cmt = 'central' with
    # rate = -2 so that dur(central) is honoured rather than the row being
    # treated as a bolus.
    d/dt(depot1) <- -ka * depot1
    d/dt(central) <- ka * depot1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot1) <- fdepot * f1
    f(central) <- fdepot * (1 - f1)
    alag(depot1) <- alag1
    dur(central) <- d2

    # ---- Plasma concentration ----
    # central is in mg and vc in L, giving mg/L; the factor 1000 converts to
    # ng/mL, the unit of the 0.5 ng/mL assay lower limit of quantification.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
