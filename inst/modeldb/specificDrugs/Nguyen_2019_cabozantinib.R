Nguyen_2019_cabozantinib <- function() {
  description <- "Updated two-compartment population PK model for oral cabozantinib (tyrosine kinase inhibitor) in healthy volunteers and patients with hepatocellular carcinoma, renal cell carcinoma, castration-resistant prostate cancer, medullary thyroid cancer, glioblastoma multiforme, or other advanced malignancies (Nguyen 2019, n=2023 across 10 clinical studies, of whom 489 had hepatocellular carcinoma). This is the final updated integrated model of Table 3, column 1, which refits the earlier Lacy 2018 integrated model after adding hepatocellular carcinoma concentration data from the phase 2 randomized discontinuation trial XL184-203 and the phase 3 CELESTIAL trial XL184-309. Absorption is described by parallel dual processes: a fraction F1 of the dose enters depot1 and is absorbed first-order at rate Ka after a lag ALAG1, while the remaining (1-F1) is delivered to central by a zero-order process of duration D2. Capsule (vs tablet reference) formulation lowers both Ka and overall relative oral bioavailability, and Ka scales with dose through a power function. Two-compartment disposition (central + peripheral1) with first-order elimination from central. Covariates on CL/F and Vc/F are baseline age, body weight, female sex, race (Black/Asian/Other vs White reference) and tumor type (HCC/RCC/CRPC/MTC/GB/Other vs healthy-volunteer reference); medullary thyroid cancer remains the only tumor type with an appreciable effect, at 1.9-fold higher CL/F, whereas hepatocellular carcinoma lowers CL/F by only 12 percent. A companion model that additionally carries NCI-ODWG liver-dysfunction covariates is available as Nguyen_2019_cabozantinib_liver_dysfunction."
  reference <- paste(
    "Nguyen L, Chapel S, Tran BD, Lacy S.",
    "Updated population pharmacokinetic model of cabozantinib integrating",
    "various cancer types including hepatocellular carcinoma.",
    "J Clin Pharmacol. 2019;59(11):1551-1561.",
    "doi:10.1002/jcph.1467.",
    "Updates the prior integrated model of Lacy S, Yang B, Nielsen J, Miles D,",
    "Nguyen L, Hutmacher M. Cancer Chemother Pharmacol. 2018;81(6):1071-1082;",
    "doi:10.1007/s00280-018-3581-0; see modellib('Lacy_2018_cabozantinib').",
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
      notes = "Time-fixed baseline. Power effect on CL/F (exponent -0.157) and Vc/F (exponent 0.0644), centered at the updated-cohort median. Nguyen 2019 Methods 'Prior Integrated PPK Model': 'The relationship between continuous covariates and typical value of PK parameters was modeled using the power function with centering by median values.' Median age 64 years is the All Studies column of Nguyen 2019 Table 2.",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight at baseline",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Power effect on CL/F (exponent -0.0393, near-null) and Vc/F (exponent 1.19, near-linear scaling), centered at the updated-cohort median. Median weight 78 kg is the All Studies column of Nguyen 2019 Table 2. Weight was missing for 7 of the 2023 subjects (Table 2 footnote a).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male, the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.76, i.e. 24 percent lower in females) and Vc/F (1.1). Nguyen 2019 Table 3 column 1. The 24 percent figure is stated verbatim in Results: 'The CL/F estimate was 24% lower in women'. Females were 317 of 2023 subjects (16 percent).",
      source_name = "SEXF"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black; White is the typical-value reference category when paired with RACE_ASIAN = 0 and RACE_OTHER = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.18) and Vc/F (1.05). Nguyen 2019 Table 3 column 1. Reference = White (1556 / 2023 = 77 percent of the pooled cohort, Table 2). Black subjects were 53 / 2023 = 3 percent.",
      source_name = "RACE"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian; White is the typical-value reference when paired with RACE_BLACK = 0 and RACE_OTHER = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.935) and Vc/F (0.696). Nguyen 2019 Table 3 column 1. Asian subjects were 211 / 2023 = 10 percent, up sharply from the prior Lacy 2018 cohort because both hepatocellular carcinoma studies enrolled heavily in Asia (35 and 33 percent Asian, Table 2).",
      source_name = "RACE"
    ),
    RACE_OTHER = list(
      description = "Race category 'Other' indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Other; White is the typical-value reference when paired with RACE_BLACK = 0 and RACE_ASIAN = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.03, near-null) and Vc/F (0.882). Nguyen 2019 Table 3 column 1. 'Other' was 47 / 2023 = 2 percent of the cohort; a further 156 / 2023 = 8 percent were recorded as race unknown and are not separately estimated by the model.",
      source_name = "RACE"
    ),
    TUMTP_HCC = list(
      description = "Hepatocellular carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-HCC; healthy volunteer is the typical-value reference when paired with all other TUMTP_* indicators = 0)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.878, i.e. 12 percent lower than healthy volunteers) and Vc/F (0.847, not statistically significant). Nguyen 2019 Table 3 column 1, 'HCC covariates' block. This covariate is the addition that distinguishes the updated model from Lacy 2018. Results: 'The magnitude of the hepatocellular carcinoma population effect on CL/F was small (12% lower CL/F comparable to 13% lower in RCC) and not likely to be clinically meaningful.' HCC patients were 489 / 2023 = 24 percent of the pooled cohort, drawn from study XL184-203 (n=37) and CELESTIAL XL184-309 (n=452).",
      source_name = "POP"
    ),
    TUMTP_RCC = list(
      description = "Renal cell carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-RCC; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.87, 13 percent lower than healthy volunteers) and Vc/F (0.656). Nguyen 2019 Table 3 column 1. RCC patients were 282 / 2023 = 14 percent of the cohort.",
      source_name = "POP"
    ),
    TUMTP_HRPC = list(
      description = "Castration-resistant prostate cancer indicator (paper writes CRPC; canonical column TUMTP_HRPC covers both HRPC and CRPC wordings)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-CRPC; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (0.989, near-null) and Vc/F (0.743). Nguyen 2019 Table 3 column 1. CRPC patients were 823 / 2023 = 41 percent of the pooled cohort, the single largest stratum.",
      source_name = "POP"
    ),
    TUMTP_MTC = list(
      description = "Medullary thyroid carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-MTC; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.9, i.e. 90 percent higher than healthy volunteers) and Vc/F (0.936). Nguyen 2019 Table 3 column 1. The 90 percent figure is stated verbatim in Results: 'MTC patients are predicted to have a 90% larger CL/F'. This remains the load-bearing covariate finding: Abstract, 'Only medullary thyroid cancer had appreciable PK differences from healthy volunteers.' MTC patients were 210 / 2023 = 10 percent of the cohort, all dosed with the 140 mg capsule.",
      source_name = "POP"
    ),
    TUMTP_GLIO = list(
      description = "Glioblastoma multiforme indicator (canonical TUMTP_GLIO covers glioma of any grade including GB)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-glioma; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.2) and Vc/F (0.479). Nguyen 2019 Table 3 column 1. GB patients were 39 / 2023 = 2 percent of the cohort.",
      source_name = "POP"
    ),
    TUMTP_OTHER = list(
      description = "Heterogeneous 'other malignancy' pool indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-other; healthy volunteer is the typical-value reference)",
      notes = "Time-fixed. Multiplicative effect on CL/F (1.19) and Vc/F (0.762). Nguyen 2019 Table 3 column 1. 'Other malignancies' were 40 / 2023 = 2 percent of the cohort, from the first-in-human study in advanced malignancies; per-subject tumor composition is not enumerated.",
      source_name = "POP"
    ),
    FORM_CAPSULE = list(
      description = "Capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (tablet; the typical-value reference)",
      notes = "Per-dose-occasion indicator. Multiplicative effect on Ka (0.402, i.e. 60 percent lower absorption rate for capsule) and on overall relative oral bioavailability (0.847, i.e. 15 percent lower). Nguyen 2019 Table 3 column 1. Methods cites the earlier capsule-vs-tablet bioequivalence finding as 'Ka and relative oral bioavailability for the capsule formulation were 58% and 14% lower than the reference tablet formulation'. Reference = tablet (Cabometyx). Capsule data were 648 / 2023 = 32 percent of subjects; Methods Results: 'Approximately one third of the data were obtained with the capsule formulation and two thirds with the tablet formulation.'",
      source_name = "FORM"
    ),
    DOSE = list(
      description = "Administered cabozantinib dose level (free base equivalent)",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-dose-occasion. Power covariate on the first-order absorption rate constant: Ka(DOSE) = Ka_ref * (DOSE / 60 mg)^0.734. Nguyen 2019 Table 3 column 1 row 'Dose-dependent Ka'. The reference dose is not printed alongside the exponent, but Methods 'Covariate Effects' fixes the model's reference condition as a subject 'receiving a 60-mg free base equivalent cabozantinib tablet dose once daily', so 60 mg is the paper's own reference dose. The same 60 mg reference is used in the sibling model Lacy_2018_cabozantinib.",
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
    hepatic_function = "Per NCI-ODWG criteria: normal 1425 (70%); mild 558 (28%); moderate 15 (1%); severe 1 (<1%); missing 24 (1%). Of the hepatocellular carcinoma patients, 99% were Child-Pugh A.",
    dose_range = "Oral cabozantinib free base equivalent; 100 mg once daily capsule in the phase 2 RDT and 60 mg once daily tablet in CELESTIAL, pooled with the wider 20-200 mg/day range of the earlier integrated analysis.",
    regions = "Multinational; the two hepatocellular carcinoma studies enrolled 33-35% Asian subjects (Nguyen 2019 Table 2)",
    formulations = "Capsule 648 subjects (32%) and tablet 1375 subjects (68%) (Nguyen 2019 Table 2, All Studies).",
    studies = c(
      "XL184-203 (phase 2 randomized discontinuation trial, HCC cohort, 100 mg QD capsule, n=37)",
      "XL184-309 CELESTIAL (phase 3, HCC after prior sorafenib, 60 mg QD tablet, n=452)",
      "Eight further studies carried over from the prior integrated analysis (phase 1 in advanced malignancies, two phase 1 studies in healthy volunteers, phase 2 in glioblastoma and in CRPC, phase 3 in MTC, CRPC and RCC); enumerated in Nguyen 2019 Supplementary Table S1"
    ),
    notes = "Baseline demographics from Nguyen 2019 Table 2. Bioanalysis by validated LC-MS/MS with a 0.5 ng/mL lower limit of quantification. The two hepatocellular carcinoma studies contributed sparse sampling only: predose troughs at the end of even weeks in XL184-203, and samples 8 or more hours after the previous dose at the week 3, 5 and 9 visits in CELESTIAL (Table 1). Percent female is computed as 317 / 2023; Table 2 rounds this to 16%."
  )

  ini({
    # ---- Structural population parameters (Nguyen 2019 Table 3, column 1
    # 'Including HCC Population', which footnote a identifies as the final
    # updated integrated model) ----
    #
    # Table 3 reports TRANSFORMED estimates throughout: the unnumbered footnote
    # states 'Transformed estimate is a PK parameter obtained by exponentiating
    # the original estimate.' The structural parameters below are therefore the
    # already-exponentiated values, and the log() here re-creates the estimated
    # THETA scale.
    #
    # Reference covariate set: tablet formulation, 60 mg dose, healthy
    # volunteer, male, White, age 64 y, weight 78 kg. Apparent (oral)
    # parameters CL/F, Vc/F, Q/F, Vp/F are in L/h and L; Ka in 1/h; ALAG1 and
    # D2 in h.
    lka <- log(1.24)
    label("First-order absorption rate constant from depot1 at the 60 mg tablet reference (1/h)") # Nguyen 2019 Table 3 col 1 Ka = 1.24 (90% CI 0.849, 1.8)
    ld2 <- log(2.48)
    label("Duration of the parallel zero-order absorption process (h)") # Nguyen 2019 Table 3 col 1 'Duration for zero-order absorption process' = 2.48 (90% CI 2.2, 2.8)
    lcl <- log(2.48)
    label("Apparent oral clearance at the reference covariate set (L/h)") # Nguyen 2019 Table 3 col 1 CL/F = 2.48 (90% CI 2.27, 2.71); Results: 'For a white male subject, CL/F at steady state was estimated as 2.48 L/h'
    lvc <- log(212)
    label("Apparent central volume of distribution at the reference covariate set (L)") # Nguyen 2019 Table 3 col 1 Vc/F = 212 (90% CI 180, 250); Results: 'and Vc/F as 212 L'
    lq <- log(30.0)
    label("Apparent inter-compartmental clearance (L/h)") # Nguyen 2019 Table 3 col 1 Q/F = 30.0 (90% CI 27.3, 33)
    lvp <- log(177)
    label("Apparent peripheral volume of distribution (L)") # Nguyen 2019 Table 3 col 1 Vp/F = 177 (90% CI 165, 189)
    lalag1 <- log(0.821)
    label("Absorption lag time on the depot1 first-order process (h)") # Nguyen 2019 Table 3 col 1 ALAG1 = 0.821 (90% CI 0.795, 0.848)

    # F1 is the fraction of the dose routed to the first-order depot1; the
    # remaining (1 - F1) is delivered by the zero-order process. Table 3
    # footnote b: 'Anti-logit transformation was used to obtain F1', so the
    # estimated THETA is on the logit scale and 0.83 is the back-transformed
    # proportion. logit(0.83) = log(0.83 / 0.17) = 1.5856.
    logitffo <- 1.5856
    label("Logit of the fraction of dose absorbed via the first-order depot1 process (F1 = expit(logitffo))") # Nguyen 2019 Table 3 col 1 'Fraction of dose in first absorption depot F1' = 0.83 (90% CI 0.80, 0.87); logit(0.83) = 1.5856

    # Overall bioavailability at the tablet reference is fixed at 1: the paper
    # estimates only the capsule-relative availability covariate, never an
    # absolute F for the reference formulation.
    lfdepot <- fixed(log(1.0))
    label("Overall relative oral bioavailability at the tablet reference (unitless)") # Nguyen 2019 Table 3 col 1 reports only 'Capsule on overall relative oral availability'; tablet is the reference formulation

    # ---- Dose-dependent Ka power exponent ----
    # Methods 'Prior Integrated PPK Model': 'The first-order absorption process
    # included a lag time and a dose-dependent effect on the absorption rate
    # (Ka) that was characterized using a power model.' The exponent is NOT
    # marked with Table 3 footnote c (the multiplicative-categorical marker)
    # and is therefore the untransformed power exponent.
    e_dose_ka <- 0.734
    label("Power exponent for dose on Ka, applied as (DOSE / 60 mg)^e_dose_ka (unitless)") # Nguyen 2019 Table 3 col 1 'Dose-dependent Ka' = 0.734 (90% CI 0.331, 1.14)

    # ---- Categorical covariate effects ----
    # Table 3 footnote c: 'For categorical covariates (eg, capsule),
    # transformed estimates correspond to multiplicative change from the
    # typical PK parameter.' Each value below is therefore a direct MULTIPLIER
    # applied as factor^indicator, so an indicator of 0 contributes 1. This
    # differs from the sibling Lacy_2018_cabozantinib model, whose Table 3
    # reported the same effects as untransformed FRACTIONAL changes.
    #
    # Four independent cross-checks against the paper's own prose confirm the
    # multiplicative reading: Female on CL/F 0.76 vs 'CL/F estimate was 24%
    # lower in women'; MTC on CL/F 1.9 vs 'MTC patients are predicted to have a
    # 90% larger CL/F'; HCC on CL/F 0.878 and RCC on CL/F 0.87 vs '12% lower
    # CL/F comparable to 13% lower in RCC'.

    # Formulation effects (capsule vs tablet reference)
    e_form_capsule_ka <- 0.402
    label("Capsule (vs tablet) multiplicative factor on Ka (unitless)") # Nguyen 2019 Table 3 col 1 'Capsule on Ka' = 0.402 (90% CI 0.223, 0.725)
    e_form_capsule_f <- 0.847
    label("Capsule (vs tablet) multiplicative factor on overall relative oral bioavailability (unitless)") # Nguyen 2019 Table 3 col 1 'Capsule on overall relative oral availability' = 0.847 (90% CI 0.83, 0.865)

    # Sex effects (female vs male reference)
    e_sexf_cl <- 0.76
    label("Female (vs male) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'Female on CL/F' = 0.76 (90% CI 0.714, 0.81)
    e_sexf_vc <- 1.1
    label("Female (vs male) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'Female on Vc/F' = 1.1 (90% CI 0.973, 1.24)

    # Race effects (Black / Asian / Other vs White reference)
    e_race_black_cl <- 1.18
    label("Race Black (vs White) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'Black on CL/F' = 1.18 (90% CI 1.04, 1.33)
    e_race_black_vc <- 1.05
    label("Race Black (vs White) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'Black on Vc/F' = 1.05 (90% CI 0.773, 1.41)
    e_race_asian_cl <- 0.935
    label("Race Asian (vs White) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'Asian on CL/F' = 0.935 (90% CI 0.869, 1.01)
    e_race_asian_vc <- 0.696
    label("Race Asian (vs White) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'Asian on Vc/F' = 0.696 (90% CI 0.558, 0.867)
    e_race_other_cl <- 1.03
    label("Race Other (vs White) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'Other race on CL/F' = 1.03 (90% CI 0.903, 1.17)
    e_race_other_vc <- 0.882
    label("Race Other (vs White) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'Other Race on Vc/F' = 0.882 (90% CI 0.615, 1.26)

    # Tumor-type effects (HCC / RCC / CRPC / MTC / GB / Other vs healthy-volunteer reference)
    e_hcc_cl <- 0.878
    label("Hepatocellular carcinoma (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'HCC on CL/F' = 0.878 (90% CI 0.794, 0.971)
    e_hcc_vc <- 0.847
    label("Hepatocellular carcinoma (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'HCC on Vc/F' = 0.847 (90% CI 0.694, 1.03)
    e_rcc_cl <- 0.87
    label("RCC (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'RCC on CL/F' = 0.87 (90% CI 0.785, 0.965)
    e_rcc_vc <- 0.656
    label("RCC (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'RCC on Vc/F' = 0.656 (90% CI 0.41, 1.05)
    e_hrpc_cl <- 0.989
    label("CRPC (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'CRPC on CL/F' = 0.989 (90% CI 0.893, 1.09)
    e_hrpc_vc <- 0.743
    label("CRPC (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'CRPC on Vc/F' = 0.743 (90% CI 0.602, 0.917)
    e_mtc_cl <- 1.9
    label("MTC (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'MTC on CL/F' = 1.9 (90% CI 1.72, 2.11)
    e_mtc_vc <- 0.936
    label("MTC (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'MTC on Vc/F' = 0.936 (90% CI 0.79, 1.11)
    e_glio_cl <- 1.2
    label("Glioblastoma multiforme (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'GB on CL/F' = 1.2 (90% CI 0.997, 1.43)
    e_glio_vc <- 0.479
    label("Glioblastoma multiforme (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'GB on Vc/F' = 0.479 (90% CI 0.333, 0.689)
    e_oth_cl <- 1.19
    label("Other malignancies (vs healthy volunteer) multiplicative factor on CL/F (unitless)") # Nguyen 2019 Table 3 col 1 'Other malignancies on CL/F' = 1.19 (90% CI 1.01, 1.4)
    e_oth_vc <- 0.762
    label("Other malignancies (vs healthy volunteer) multiplicative factor on Vc/F (unitless)") # Nguyen 2019 Table 3 col 1 'Other malignancies on Vc/F' = 0.762 (90% CI 0.593, 0.979)

    # ---- Continuous covariate effects (power exponents, untransformed) ----
    # These rows carry no footnote c marker, and the negative age and weight
    # values could not be exponentiated estimates, so they are the raw power
    # exponents of the centered power model described in Methods.
    e_age_cl <- -0.157
    label("Power exponent for age on CL/F, applied as (AGE / 64 y)^e_age_cl (unitless)") # Nguyen 2019 Table 3 col 1 'Age on CL/F' = -0.157 (90% CI -0.264, -0.0509)
    e_age_vc <- 0.0644
    label("Power exponent for age on Vc/F, applied as (AGE / 64 y)^e_age_vc (unitless)") # Nguyen 2019 Table 3 col 1 'Age on Vc/F' = 0.0644 (90% CI -0.148, 0.277)
    e_wt_cl <- -0.0393
    label("Power exponent for body weight on CL/F, applied as (WT / 78 kg)^e_wt_cl (unitless)") # Nguyen 2019 Table 3 col 1 'Weight on CL/F' = -0.0393 (90% CI -0.147, 0.0679)
    e_wt_vc <- 1.19
    label("Power exponent for body weight on Vc/F, applied as (WT / 78 kg)^e_wt_vc (unitless)") # Nguyen 2019 Table 3 col 1 'Weight on Vc/F' = 1.19 (90% CI 0.934, 1.46)

    # ---- Inter-individual variability ----
    # Nguyen 2019 Table 3 col 1 'Variance' block. The table footnote defines
    # omega^2 as the 'variance of population parameter modeled using
    # exponential model', so these are log-scale variances and the
    # colon-separated CL/F:Vc/F row is the OMEGA BLOCK off-diagonal covariance.
    #
    # Consistency check against the paper's own prose: Results reports
    # 'approximately 46% for CL/F and 67% for Vc/F'; sqrt(0.213) = 0.462 and
    # sqrt(0.443) = 0.666, matching both figures.
    #
    # Cauchy-Schwarz holds for the CL/Vc block: sqrt(0.213 * 0.443) = 0.307 and
    # the reported covariance 0.211 is below it, giving a valid correlation of
    # 0.69. (The corresponding off-diagonal in the prior Lacy 2018 model was
    # not admissible and had to be dropped; here it can be carried faithfully.)
    etalcl + etalvc ~ c(0.213, 0.211, 0.443) # Nguyen 2019 Table 3 col 1: omega^2 CL/F = 0.213; omega^2 CL/F:Vc/F = 0.211; omega^2 Vc/F = 0.443
    etalka ~ 2.02 # Nguyen 2019 Table 3 col 1 omega^2 Ka = 2.02 (90% CI 1.59, 2.45)
    etalogitffo ~ 2.55 # Nguyen 2019 Table 3 col 1 omega^2 F1 = 2.55 (90% CI 1.99, 3.1), on the logit scale per footnote b

    # ---- Residual error ----
    # Methods: 'The residual variability was modeled using the log-transformed
    # additive-error model.' An additive residual on log-transformed
    # concentration is equivalent to a proportional residual in the linear
    # concentration space nlmixr2 works in. sqrt(0.127) = 0.35637.
    propSd <- 0.35637
    label("Proportional residual error (fraction; sqrt of the log-scale additive variance)") # Nguyen 2019 Table 3 col 1 sigma^2 = 0.127 (90% CI 0.123, 0.131)
  })

  model({
    # ---- Reference (centering) covariate values ----
    # Methods 'Prior Integrated PPK Model': continuous covariates enter as a
    # 'power function with centering by median values'. The medians are the
    # All Studies column of Table 2. The reference dose is the 60 mg free base
    # equivalent tablet named in Methods 'Covariate Effects'.
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

    ka_form <- e_form_capsule_ka^FORM_CAPSULE
    f_form <- e_form_capsule_f^FORM_CAPSULE

    # ---- Individual PK parameters ----
    ka <- exp(lka + etalka) *
      (DOSE / ref_dose)^e_dose_ka *
      ka_form
    cl <- exp(lcl + etalcl) *
      (AGE / ref_age)^e_age_cl *
      (WT / ref_wt)^e_wt_cl *
      cl_sex * cl_race * cl_tumtp
    vc <- exp(lvc + etalvc) *
      (AGE / ref_age)^e_age_vc *
      (WT / ref_wt)^e_wt_vc *
      vc_sex * vc_race * vc_tumtp
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
    # directly into central by a zero-order process of duration d2 (no holding
    # depot). The overall relative bioavailability fdepot scales both routes.
    # Dosing events must therefore target BOTH compartments with the same
    # nominal amount: cmt = 'depot1' as an ordinary oral dose, and
    # cmt = 'central' with rate = -2 so that dur(central) is honoured rather
    # than the row being treated as a bolus.
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
