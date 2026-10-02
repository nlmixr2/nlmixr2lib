Cheng_2021_avadomide <- function() {
  description <- "Two-compartment population PK model for oral avadomide (CC-122, a cereblon-modulating agent) in healthy adults, adults with renal impairment, and adults with advanced solid tumors, non-Hodgkin lymphoma or multiple myeloma (Cheng 2021, n=298 across 3 studies). First-order absorption with an absorption lag time and first-order elimination; the peripheral volume is fixed at 10 L. Linear creatinine-clearance and tumor-type (DLBCL, PCNSL, other solid tumor, MM) effects on CL/F; linear body-weight, female-sex and tumor-type (DLBCL, PCNSL, NHL, HCC, GBM, MM, brain cancer) effects on V2/F. Healthy subjects are the tumor-type reference."
  reference <- paste(
    "Cheng Y, Chen J, Pourdehnad M, Zhou S, Li Y.",
    "Population Pharmacokinetics of CC-122.",
    "Clin Pharmacol. 2021;13:61-71.",
    "doi:10.2147/CPAA.S310604.",
    sep = " "
  )
  vignette <- "Cheng_2021_avadomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "avadomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "avadomide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avadomide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Baseline creatinine clearance (estimating equation not stated; raw mL/min, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Linear centred effect on CL/F: (1 + 0.007 * (CRCL - 94.42)) per Table 2 footnote b. The reference 94.42 mL/min is the pooled-cohort median (Table 1: 94.4, range 9.0-321.2 mL/min). The paper does not state the estimating equation or a BSA normalization; values are carried in mL/min as printed in Table 1. The Results text writes the renal-function bands as 'mL/hr', a typo for mL/min (Table 1 and Figure 4A use mL/min).",
      source_name = "CLcr"
    ),
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Linear centred effect on V2/F: (1 + 0.009 * (WT - 74.5)) per Table 2 footnote c. The reference 74.5 kg is the pooled-cohort median (Table 1, range 39.8-159.0 kg).",
      source_name = "BW"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male, the most prevalent category and the typical-value reference)",
      notes = "Time-fixed. Multiplicative fractional effect on V2/F: (1 - 0.179) if female (Table 2 footnote c). Males were 62.1% of the cohort (Table 1).",
      source_name = "Sex"
    ),
    TUMTP_DLBCL = list(
      description = "Diffuse large B-cell lymphoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference when every TUMTP_* indicator = 0)",
      notes = "Time-fixed. Fractional effects: CL/F x (1 - 0.647), V2/F x (1 + 0.476) (Table 2). n = 60 (20.1%).",
      source_name = "Tumor"
    ),
    TUMTP_PCNSL = list(
      description = "Primary CNS lymphoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effects: CL/F x (1 - 0.692), V2/F x (1 + 0.846) (Table 2). n = 5 (1.7%). Coded separately from TUMTP_DLBCL and TUMTP_NHL; the paper's tumor-type levels are mutually exclusive.",
      source_name = "Tumor"
    ),
    TUMTP_OTHER = list(
      description = "Other solid tumor indicator (solid tumors outside the HCC, GBM and brain-cancer categories)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effect on CL/F x (1 - 0.364) (Table 2). The effect on V2/F was fixed to 0 in the final model (Results: 'effect of this tumor type (other solid tumor) was not included in the final model (fix to 0)'). n = 19 (6.4%). Named strata in this paper's decomposition: healthy (reference), DLBCL, PCNSL, NHL, HCC, GBM, MM, brain cancer; the histologies pooled into 'other solid tumor' are not enumerated.",
      source_name = "Tumor"
    ),
    TUMTP_MYELO = list(
      description = "Multiple myeloma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effects: CL/F x (1 + 0.309), V2/F x (1 + 0.682) (Table 2). n = 29 (9.7%).",
      source_name = "Tumor"
    ),
    TUMTP_NHL = list(
      description = "Non-Hodgkin lymphoma indicator (NHL other than DLBCL and PCNSL)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effect on V2/F x (1 + 0.344) (Table 2). The effect on CL/F was fixed to 0 (Results: 'effects of these tumor types were not incorporated in the model to ensure stability (fix to 0)'). n = 30 (10.1%). In this paper DLBCL and PCNSL carry their own indicators, so TUMTP_NHL = 1 only for the remaining NHL patients; their histologies are not enumerated.",
      source_name = "Tumor"
    ),
    TUMTP_HCC = list(
      description = "Hepatocellular carcinoma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effect on V2/F x (1 + 0.521) (Table 2). The effect on CL/F was fixed to 0 (Results). n = 27 (9.1%).",
      source_name = "Tumor"
    ),
    TUMTP_GLIO = list(
      description = "Glioblastoma multiforme indicator (canonical TUMTP_GLIO covers glioma including GBM)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effect on V2/F x (1 + 0.480) (Table 2). The effect on CL/F was fixed to 0 (Results). n = 44 (14.8%). The paper's separate 'brain cancer' category is TUMTP_BRAIN, not this indicator.",
      source_name = "Tumor"
    ),
    TUMTP_BRAIN = list(
      description = "Brain cancer indicator (brain tumors other than GBM and PCNSL)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy subject is the typical-value reference)",
      notes = "Time-fixed. Fractional effect on V2/F x (1 + 0.415) (Table 2). The effect on CL/F was fixed to 0 (Results: 'relatively limited number of patients but large variabilities in brain cancer (n=6) cohort'). n = 6 (2.0%). Histologies not enumerated; GBM (TUMTP_GLIO) and PCNSL (TUMTP_PCNSL) are separate categories.",
      source_name = "Tumor"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate model (SCM) but not retained; Results: 'no clinically meaningful effect ... age (20 to 91 years)'."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened in the SCM but not retained (Results)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened in the SCM but not retained (Results)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened in the SCM but not retained (Results)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened in the SCM but not retained (Results)."
    ),
    BILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened in the SCM but not retained (Results)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 298L,
    n_studies = 3L,
    age_range = "20-91 years",
    age_median = "59.5 years",
    weight_range = "39.8-159.0 kg",
    weight_median = "74.5 kg",
    sex_female_pct = 37.9,
    disease_state = "Healthy adults (26.2%; CC-122-CP-002 and the matched-control arm of CC-122-CP-005), adults with mild, moderate or severe renal impairment (CC-122-CP-005), and patients with advanced solid tumors, NHL or MM (CC-122-ST-001): DLBCL 20.1%, GBM 14.8%, NHL 10.1%, MM 9.7%, HCC 9.1%, other solid tumor 6.4%, brain cancer 2.0%, PCNSL 1.7%.",
    dose_range = "0.5-15 mg oral CC-122 (single doses in healthy and renally impaired subjects; daily dosing in patients).",
    renal_function = "Creatinine clearance median 94.4 mL/min (range 9.0-321.2).",
    studies = c(
      "CC-122-CP-002 Part 1 (phase 1 single ascending dose, healthy adults, n=30)",
      "CC-122-CP-005 (phase 1 single-dose renal impairment and matched healthy subjects, n=48)",
      "CC-122-ST-001 / NCT01421524 (phase 1a/b dose finding, advanced solid tumors, NHL, MM, n=220)"
    ),
    notes = "Demographics from Cheng 2021 Table 1. CC-122 is racemic; the achiral assay (CP-002, CP-005; LLOQ 1.0 ng/mL) measured total CC-122 and the chiral assay (ST-001; LLOQ 0.5 ng/mL per enantiomer) measured the R and S enantiomers. Race and region are not reported."
  )

  ini({
    # ---- Structural parameters (Cheng 2021 Table 2, final model) ----
    # Reference subject: healthy, male, CLcr 94.42 mL/min, body weight 74.5 kg.
    lka <- log(4.14); label("First-order absorption rate constant (1/h)") # Table 2 'TVKa' = 4.14 (RSE 20.5%)
    ltlag <- log(0.246); label("Absorption lag time (h)") # Table 2 'TVALAG' = 0.246 (RSE 0.6%)
    lcl <- log(3.63); label("Apparent clearance CL/F at the reference covariates (L/h)") # Table 2 'TVCL/F' = 3.63 (RSE 3.5%)
    lvc <- log(36.2); label("Apparent central volume V2/F at the reference covariates (L)") # Table 2 'TVV2/F' = 36.2 (RSE 5.7%)
    lq <- log(1.38); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 'TVQ/F' = 1.38 (RSE 49.9%)
    lvp <- fixed(log(10)); label("Apparent peripheral volume V3/F (L)") # Table 2 'TVV3/F' = 10 Fix; Results: V3 fixed to 10 L as estimated from the healthy-subject studies

    # ---- Covariate effects on CL/F (Table 2 footnote b) ----
    # CL/F = 3.63 * (1 + 0.007*(CLcr - 94.42)) * (1 - 0.647)^DLBCL *
    #        (1 - 0.692)^PCNSL * (1 - 0.364)^OtherSolidTumor * (1 + 0.309)^MM
    e_crcl_cl <- 0.007; label("Linear CLcr effect on CL/F (per mL/min above 94.42)") # Table 2 'CLcr on CL/F' = 0.007 (RSE 65%)
    e_tumtp_dlbcl_cl <- -0.647; label("DLBCL fractional change on CL/F (unitless)") # Table 2 'Tumor type (DLBCL) on CL/F' = -0.647
    e_tumtp_pcnsl_cl <- -0.692; label("PCNSL fractional change on CL/F (unitless)") # Table 2 'Tumor type (PCNSL) on CL/F' = -0.692
    e_tumtp_other_cl <- -0.364; label("Other solid tumor fractional change on CL/F (unitless)") # Table 2 'Tumor type (other solid tumor) on CL/F' = -0.364
    e_tumtp_myelo_cl <- 0.309; label("MM fractional change on CL/F (unitless)") # Table 2 'Tumor type (MM) on CL/F' = 0.309

    # ---- Covariate effects on V2/F (Table 2 footnote c) ----
    e_wt_vc <- 0.009; label("Linear body-weight effect on V2/F (per kg above 74.5)") # Table 2 'Body weight on V2/F' = 0.009 (RSE 47%)
    e_sexf_vc <- -0.179; label("Female fractional change on V2/F (unitless)") # Table 2 'Sex on V2/F' = -0.179; footnote c '(1 - 0.179) if female'
    e_tumtp_dlbcl_vc <- 0.476; label("DLBCL fractional change on V2/F (unitless)") # Table 2 'Tumor type (DLBCL) on V2/F' = 0.476
    e_tumtp_pcnsl_vc <- 0.846; label("PCNSL fractional change on V2/F (unitless)") # Table 2 'Tumor type (PCNSL) on V2/F' = 0.846
    e_tumtp_nhl_vc <- 0.344; label("NHL fractional change on V2/F (unitless)") # Table 2 'Tumor type (NHL) on V2/F' = 0.344
    e_tumtp_hcc_vc <- 0.521; label("HCC fractional change on V2/F (unitless)") # Table 2 'Tumor type (HCC) on V2/F' = 0.521
    e_tumtp_glio_vc <- 0.480; label("GBM fractional change on V2/F (unitless)") # Table 2 'Tumor type (GBM) on V2/F' = 0.480
    e_tumtp_myelo_vc <- 0.682; label("MM fractional change on V2/F (unitless)") # Table 2 'Tumor type (MM) on V2/F' = 0.682 (estimate column and footnote c; bootstrap median 0.662)
    e_tumtp_brain_vc <- 0.415; label("Brain cancer fractional change on V2/F (unitless)") # Table 2 'Tumor type (brain cancer) on V2/F' = 0.415

    # ---- Inter-individual variability ----
    # Table 2 reports CV% = sqrt(omega^2) * 100: the bootstrap-median column
    # prints the variances, and sqrt(0.068) = 26.1%, sqrt(2.05) = 143.2%,
    # sqrt(0.962) = 98.1% reproduce the printed CV% exactly. omega^2 is
    # therefore (CV%/100)^2. No covariances were reported (diagonal OMEGA).
    etalcl ~ 0.272484 # Table 2 'CV%, CL/F' = 52.2% -> 0.522^2 (bootstrap median 0.261)
    etalvc ~ 0.068121 # Table 2 'CV%, V2/F' = 26.1% -> 0.261^2 (bootstrap median 0.068)
    etalka ~ 2.050624 # Table 2 'CV%, Ka' = 143.2% -> 1.432^2 (bootstrap median 2.05)
    etalq ~ 0.962361 # Table 2 'CV%, Q' = 98.1% -> 0.981^2 (bootstrap median 0.962)

    # ---- Residual error ----
    # Methods: additive error on log-transformed concentrations; Table 2
    # 'sigma^2 (Log additive)' = 0.092 -> SD = sqrt(0.092) = 0.30332.
    expSd <- 0.30332; label("Additive residual SD on the log-concentration scale") # Table 2 'sigma^2 (Log additive)' = 0.092
  })

  model({
    # Reference covariate values (Table 2 footnotes b and c)
    ref_crcl <- 94.42
    ref_wt <- 74.5

    # Categorical covariate multipliers: P = theta * (1 + theta_cov * Z)
    # (Methods Eq. 3); the reference (healthy, male) gives multiplier 1.
    cl_tumtp <- (1 + e_tumtp_dlbcl_cl * TUMTP_DLBCL) *
      (1 + e_tumtp_pcnsl_cl * TUMTP_PCNSL) *
      (1 + e_tumtp_other_cl * TUMTP_OTHER) *
      (1 + e_tumtp_myelo_cl * TUMTP_MYELO)
    vc_tumtp <- (1 + e_tumtp_dlbcl_vc * TUMTP_DLBCL) *
      (1 + e_tumtp_pcnsl_vc * TUMTP_PCNSL) *
      (1 + e_tumtp_nhl_vc * TUMTP_NHL) *
      (1 + e_tumtp_hcc_vc * TUMTP_HCC) *
      (1 + e_tumtp_glio_vc * TUMTP_GLIO) *
      (1 + e_tumtp_myelo_vc * TUMTP_MYELO) *
      (1 + e_tumtp_brain_vc * TUMTP_BRAIN)

    # Individual parameters; continuous covariates enter linearly
    # (Methods Eq. 1: P = theta * (1 + theta_cov * (COV - COVm)))
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - ref_crcl)) * cl_tumtp
    vc <- exp(lvc + etalvc) * (1 + e_wt_vc * (WT - ref_wt)) *
      (1 + e_sexf_vc * SEXF) * vc_tumtp
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # central in mg, vc in L -> mg/L; x 1000 -> ng/mL (LLOQ 1.0 ng/mL)
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
