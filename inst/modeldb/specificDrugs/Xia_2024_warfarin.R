Xia_2024_warfarin <- function() {
  description <- "K-PD warfarin PK/PD model for adult Han Chinese (Alfalfa-Warfarin-PPK/PD; Xia 2024). PK parameters fixed from the Hamberg model; PD EC50 re-estimated, with VKORC1 -1639 G/A and CYP2C9 *1/*2/*3 allele-specific contributions, body-weight power scaling, and amiodarone effect on EC50. Two parallel coagulation-factor transit chains drive INR."
  reference <- paste(
    "Xia X, Cai X, Chen J, Jiang S, Zhang J.",
    "Construction of warfarin population pharmacokinetics and pharmacodynamics model",
    "in Han population based on Bayesian method.",
    "Sci Rep. 2024;14:14894. doi:10.1038/s41598-024-65048-7.",
    "PK structural parameters (CL per CYP2C9*1/*2/*3 allele, age effect on CL, V/F,",
    "Emax, gamma, MTT1, MTT2, INRmax) are fixed from the Hamberg model",
    "as reported in Table 3 of Xia 2024."
  )
  vignette <- "Xia_2024_warfarin"
  paper_specific_compartments <- c("coag_s1", "coag_s2", "coag_s3", "coag_l1", "coag_l2", "coag_l3")

  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline; power scaling on EC50 with reference 60 kg (cohort median; Xia 2024 supplement Section 1.3 Eq 8).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline; linear scale model on CL with reference age 71 years (Xia 2024 supplement Eq 7; Hamberg literature value).",
      source_name        = "AGE"
    ),
    CYP2C9_S1_COUNT = list(
      description        = "Count of CYP2C9*1 (wild-type) alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "CYP2C9_S1_COUNT + CYP2C9_S2_COUNT + CYP2C9_S3_COUNT = 2. Each *1 allele contributes 0.174 L/h to CL (Hamberg fixed value).",
      source_name        = "CYP2C9 (genotype string mapped to per-allele counts)"
    ),
    CYP2C9_S2_COUNT = list(
      description        = "Count of CYP2C9*2 alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Each *2 allele contributes 0.0879 L/h to CL (Hamberg fixed value). No *2 carriers were reported in the Xia 2024 Han cohort (Table 1) but the model retains the term for general use.",
      source_name        = "CYP2C9 (genotype string mapped to per-allele counts)"
    ),
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Each *3 allele contributes 0.0422 L/h to CL (Hamberg fixed value). Xia 2024 reported 5.7% of subjects as *1/*3 heterozygous (Table 1).",
      source_name        = "CYP2C9 (genotype string mapped to per-allele counts)"
    ),
    VKORC1_1639G_COUNT = list(
      description        = "Count of VKORC1 -1639G alleles per subject (0, 1, or 2). The complementary -1639A count is 2 - VKORC1_1639G_COUNT.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Each -1639G allele typical EC50 contribution is estimated at 4.3 mg/L; each -1639A allele contribution is 1.14 mg/L. Xia 2024 Han cohort distribution (Table 1): AA 80.3%, GA 18.7%, GG 0.9%.",
      source_name        = "VKORC1 (genotype string mapped to G allele count)"
    ),
    CONMED_AMIO = list(
      description        = "Concomitant amiodarone use indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant amiodarone)",
      notes              = "Multiplicative piecewise effect on EC50: ec50 *= (1 + e_amio_ec50 * CONMED_AMIO) per Xia 2024 supplement Section 1.3 Eq 9 with e_amio_ec50 = -0.602. Xia 2024 Han cohort prevalence (Table 1): 25.7% on amiodarone.",
      source_name        = "CM1"
    ),
    INR_BASE = list(
      description        = "Subject-specific baseline INR measured before warfarin administration",
      units              = "(unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (pre-medication INR). Added directly to the predicted INR equation; simulation default uses the Xia 2024 total-cohort mean of 1.13 (supplement Section 1.6).",
      source_name        = "INR_BASE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 646,
    n_studies      = 9,
    age_range      = "1-82 years",
    age_median     = "55.2 years (mean +/- SD 54.7 +/- 12.5 in the total cohort)",
    weight_range   = "7-106 kg",
    weight_median  = "59.7 kg (mean +/- SD 59.7 +/- 11.6 in the total cohort)",
    sex_female_pct = 48.9,
    race_ethnicity = "Han Chinese (100%)",
    disease_state  = "Adults requiring warfarin anticoagulation for valvular disease, post-valve replacement, non-valvular atrial fibrillation, or deep-vein thrombosis",
    dose_range     = "0.625-6.25 mg orally once daily",
    regions        = "China (Fujian Medical University Union Hospital plus 8 sub-centres, 2016-2020)",
    cyp2c9_freq    = "CYP2C9*1/*1 94.3%, *1/*3 5.7% (Xia 2024 Table 1, total cohort)",
    vkorc1_freq    = "VKORC1 -1639AA 80.3%, AG 18.7%, GG 0.9% (Xia 2024 Table 1, total cohort)",
    inr_baseline   = "1.13 +/- 0.59 (total cohort mean +/- SD; Xia 2024 Table 1)",
    n_inr_records  = 5467,
    notes          = "Pooled dataset across 9 Chinese centres; n=537 used for modelling, n=109 for external validation. Baseline demographics from Xia 2024 Table 1."
  )

  ini({
    # ---- PK parameters (all fixed from Hamberg model per Xia 2024 Table 3) ----
    # Per-CYP2C9-allele CL contributions; total subject CL is the sum across two alleles
    lcl_cyp2c9_s1 <- fixed(log(0.174))  ; label("CL per CYP2C9*1 allele (L/h)")              # Xia 2024 Table 3 (fixed, Hamberg literature value)
    lcl_cyp2c9_s2 <- fixed(log(0.0879)) ; label("CL per CYP2C9*2 allele (L/h)")              # Xia 2024 Table 3 (fixed, Hamberg literature value)
    lcl_cyp2c9_s3 <- fixed(log(0.0422)) ; label("CL per CYP2C9*3 allele (L/h)")              # Xia 2024 Table 3 (fixed, Hamberg literature value)
    e_age_cl      <- fixed(-0.00571)    ; label("Age effect on CL (fractional change per year, reference age 71 y)")  # Xia 2024 Table 3 / supplement Eq 7 (fixed)
    lvc           <- fixed(log(14.3))   ; label("Apparent volume of distribution V/F (L) at effect site")             # Xia 2024 Table 3 (fixed, Hamberg literature value)
    lemax         <- fixed(log(1))      ; label("Maximum anticoagulant effect Emax (fixed at 1 = 100%)")              # Xia 2024 supplement Section 1.1 (fixed)
    lhill        <- fixed(log(1.39))   ; label("Hill coefficient gamma for the sigmoid Emax effect (unitless)")      # Xia 2024 Table 3 (fixed)
    lmtt1         <- fixed(log(27.2))   ; label("Mean transit time of fast coagulation-factor chain MTT1 (h)")        # Xia 2024 Table 3 (fixed)
    lmtt2         <- fixed(log(110.9))  ; label("Mean transit time of slow coagulation-factor chain MTT2 (h)")        # Xia 2024 Table 3 (fixed)
    linrmax       <- fixed(log(20))     ; label("Maximum INR increment above baseline (unitless)")                    # Xia 2024 Table 3 / supplement Section 1.1 (fixed)

    # ---- PD parameters (re-estimated for Han Chinese cohort) ----
    lec50_vkorc1_g <- log(4.3)          ; label("Typical EC50 per VKORC1 -1639G allele (mg/L)")    # Xia 2024 Table 3 (estimated, RSE 16.3%)
    lec50_vkorc1_a <- log(1.14)         ; label("Typical EC50 per VKORC1 -1639A allele (mg/L)")    # Xia 2024 Table 3 (estimated, RSE 7.7%)
    e_wt_ec50      <- 1.34              ; label("Power exponent of body weight on EC50 (reference 60 kg)")  # Xia 2024 Table 3 (estimated, RSE 34.3%)
    e_amio_ec50    <- -0.602            ; label("Fractional multiplicative effect of amiodarone on EC50")   # Xia 2024 Table 3 (estimated, RSE 46.7%); piecewise per supplement Eq 9

    # ---- IIV: only EC50 has IIV in the final Xia 2024 model (Table 3) ----
    # omega^2 = log(CV^2 + 1) with CV(EC50) = 55.8% from Table 3
    etalec50 ~ 0.2712                                                                            # log(1 + 0.558^2) ~= 0.2712; Xia 2024 Table 3 reports CV = 55.8% (RSE 8.2%)

    # ---- Residual error (proportional on INR) ----
    propSd_INR <- 0.365                  ; label("Proportional residual error on INR (fraction)")  # Xia 2024 Table 3 (estimated, RSE 4.1%); supplement Eq 4
  })

  model({
    # Per-allele CL contributions (linear scale)
    cl_s1 <- exp(lcl_cyp2c9_s1)
    cl_s2 <- exp(lcl_cyp2c9_s2)
    cl_s3 <- exp(lcl_cyp2c9_s3)

    # Individual CL: per-allele sum scaled by age effect (Hamberg form; ref age 71 y)
    cl <- (CYP2C9_S1_COUNT * cl_s1 + CYP2C9_S2_COUNT * cl_s2 + CYP2C9_S3_COUNT * cl_s3) *
          (1 + e_age_cl * (AGE - 71))

    # Volume and derived constants
    vc      <- exp(lvc)
    emax    <- exp(lemax)
    hill <- exp(lhill)
    mtt1    <- exp(lmtt1)
    mtt2    <- exp(lmtt2)
    inrmax  <- exp(linrmax)

    # Individual typical EC50 from VKORC1 genotype (G allele dosage + complementary A count)
    vkorc1_a_count <- 2 - VKORC1_1639G_COUNT
    ec50_typ_geno  <- VKORC1_1639G_COUNT * exp(lec50_vkorc1_g) +
                      vkorc1_a_count     * exp(lec50_vkorc1_a)

    # Apply weight (power) and amiodarone (piecewise multiplicative) covariate effects, then IIV
    ec50 <- ec50_typ_geno *
            (WT / 60)^e_wt_ec50 *
            (1 + e_amio_ec50 * CONMED_AMIO) *
            exp(etalec50)

    # K-PD elimination of administered dose from a single virtual depot
    kde <- cl / vc

    # Effective drug-delivery-rate "concentration": cdr = depot/vc gives cdr/ec50 = DR/EDK50
    cdr <- depot / vc
    e_anticoag <- emax * cdr^hill / (ec50^hill + cdr^hill)

    # Hamberg transit-chain rate constants: ktr = N / MTT_total with N = 3 transit compartments
    ktr_s <- 3 / mtt1
    ktr_l <- 3 / mtt2

    # ODE system
    d/dt(depot)   <- -kde * depot
    d/dt(coag_s1) <- ktr_s * (1 - e_anticoag - coag_s1)
    d/dt(coag_s2) <- ktr_s * (coag_s1 - coag_s2)
    d/dt(coag_s3) <- ktr_s * (coag_s2 - coag_s3)
    d/dt(coag_l1) <- ktr_l * (1 - e_anticoag - coag_l1)
    d/dt(coag_l2) <- ktr_l * (coag_l1 - coag_l2)
    d/dt(coag_l3) <- ktr_l * (coag_l2 - coag_l3)

    # Initial coagulation-factor levels normalized to 1 in the no-drug state (Hamberg convention)
    coag_s1(0) <- 1
    coag_s2(0) <- 1
    coag_s3(0) <- 1
    coag_l1(0) <- 1
    coag_l2(0) <- 1
    coag_l3(0) <- 1

    # INR observation: subject-specific baseline plus inrmax x anticoagulation fraction
    INR <- INR_BASE + inrmax * (1 - (coag_s3 + coag_l3) / 2)
    INR ~ prop(propSd_INR)
  })
}
