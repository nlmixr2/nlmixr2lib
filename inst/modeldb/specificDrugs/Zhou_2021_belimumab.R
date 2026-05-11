Zhou_2021_belimumab <- function() {
  description <- "Linear two-compartment IV population PK model for belimumab in Chinese and non-Chinese adult and pediatric patients with systemic lupus erythematosus (Zhou 2021)"
  reference <- "Zhou L, Lee S, Zhu L, Roy A, Zhou H, Yang H. Prediction of Belimumab Pharmacokinetics in Chinese Pediatric Patients with Systemic Lupus Erythematosus. Drugs R D. 2021;21(4):407-417. doi:10.1007/s40268-021-00363-2"
  vignette <- "Zhou_2021_belimumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL and V1 with reference 40.69 kg (Zhou 2021 Table 2). Computed per Janmahasatian et al. Clin Pharmacokinet 2005;44:1051-1065 from height, weight, and sex.",
      source_name        = "FFM"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column BALB (baseline albumin). Power-form effect on CL with reference 40 g/L (Zhou 2021 Table 2).",
      source_name        = "BALB"
    ),
    IGG = list(
      description        = "Baseline serum immunoglobulin G",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column BIGG (baseline IgG). Power-form effect on CL with reference 14.8 g/L (Zhou 2021 Table 2).",
      source_name        = "BIGG"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Saturable maturation effect on V1: factor AGE/(AGE+1.58). Approaches 1 in adults, ~0.79 in 6-year-olds (Zhou 2021 Table 2).",
      source_name        = "AGE"
    ),
    RACE_NEAS = list(
      description        = "North East Asian race indicator (Chinese, Japanese, or Korean heritage)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-North East Asian)",
      notes              = "Multiplicative effect on V1 (factor 1.07) when RACE_NEAS = 1. Source column RAC4 (Zhou 2021 Table 2 footnote d).",
      source_name        = "RAC4"
    ),
    STUDY_LBSL = list(
      description        = "Indicator for early-phase belimumab studies LBSL01 or LBSL02",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any other belimumab study in the pooled analysis)",
      notes              = "Multiplicative effects on CL (factor 1.63) and V1 (factor 1.26) when STUDY_LBSL = 1. Source column INDR (Zhou 2021 Table 2 footnote: INDR=1 for LBSL01/LBSL02, 0 otherwise). The two early studies used a different ELISA assay than the later electrochemiluminescence-based assay.",
      source_name        = "INDR"
    )
  )

  population <- list(
    n_subjects     = 1783,
    n_studies      = 9,
    n_observations = 9650,
    age_range      = "6-80 years",
    age_median     = "38 years (adults); 14 years (pediatric)",
    weight_range   = "17-165.4 kg",
    weight_median  = "65.7 kg (adults); 52.5 kg (pediatric)",
    sex_female_pct = 93.1,
    race_ethnicity = "North East Asian 13.2% (235/1783); other 86.8% (1548/1783).",
    disease_state  = "Systemic lupus erythematosus (SLE), plus a small contribution from non-Chinese healthy volunteers in the early-phase studies.",
    dose_range     = "IV 1-10 mg/kg every 4 weeks (final-model adult/pediatric regimen: 10 mg/kg Q4W).",
    regions        = "Multinational pooled analysis combining Chinese and non-Chinese cohorts.",
    pediatric_n    = 53,
    adult_n        = 1730,
    studies        = "BEL114055 (PLUTO; NCT01649765), LBSL01 (NCT00657007), LBSL02 (NCT00071487), BEL113750 (BLISS-NEA; NCT01345253), BEL116119 (NCT01516450), BEL110751 (BLISS-76; NCT00410384), BEL110752 (BLISS-52; NCT00424476), 200909 (NCT02880852), BEL114448 (NCT01583530).",
    notes          = "Zhou 2021 Table 1 baseline demographics. Median (min, max) FFM 41.08 (24.32; 87.93) kg in adults and 34.45 (12.56; 57.16) kg in pediatric; median serum IgG 14.7 g/L (adults) / 14.5 g/L (pediatric); median albumin 40 g/L (adults) / 43 g/L (pediatric)."
  )

  ini({
    # Structural parameters — typical values from Zhou 2021 Table 2 (final model run613)
    # Published in mL/day and mL; converted to L/day and L by dividing by 1000.
    lcl <- log(0.238); label("Clearance for the reference adult (CL, L/day)")        # Table 2: theta1 = 238 mL/day
    lvc <- log(2.597); label("Central volume of distribution for the reference adult (V1, L)")  # Table 2: theta2 = 2597 mL
    lq  <- log(0.591); label("Intercompartmental clearance (Q, L/day)")              # Table 2: theta3 = 591 mL/day
    lvp <- log(2.318); label("Peripheral volume of distribution (V2, L)")            # Table 2: theta4 = 2318 mL

    # Covariate effects on CL (Zhou 2021 Table 2)
    e_ffm_cl  <-  0.673; label("FFM power exponent on CL (reference FFM 40.69 kg)")           # Table 2: theta7
    e_alb_cl  <- -1.12;  label("Baseline albumin power exponent on CL (reference 40 g/L)")    # Table 2: theta9
    e_igg_cl  <-  0.293; label("Baseline IgG power exponent on CL (reference 14.8 g/L)")      # Table 2: theta10
    e_stdy_cl <-  1.63;  label("Multiplicative factor on CL when STUDY_LBSL = 1 (LBSL01/02)")  # Table 2: theta11

    # Covariate effects on V1 (Zhou 2021 Table 2)
    e_ffm_vc  <-  0.891; label("FFM power exponent on V1 (reference FFM 40.69 kg)")               # Table 2: theta8
    e_stdy_vc <-  1.26;  label("Multiplicative factor on V1 when STUDY_LBSL = 1 (LBSL01/02)")      # Table 2: theta12
    e_neas_vc <-  1.07;  label("Multiplicative factor on V1 when RACE_NEAS = 1 (North East Asian)") # Table 2: theta13 (printed as theta12 due to a labelling typo)
    age50_vc  <-  1.58;  label("Age (years) at half-maximal V1 in the saturable AGE/(AGE+age50_vc) maturation term") # Table 2: theta14

    # IIV (omega^2 from Table 2; CL and V2 correlated, no IIV on Q)
    etalcl + etalvp ~ c(0.0682,
                        0.0125, 0.0781)
    etalvc ~ 0.0079

    # Residual error (Zhou 2021 Table 2: epsilon1 proportional, epsilon2 additive)
    propSd <- 0.247; label("Proportional residual error (fraction)") # Table 2: theta5
    addSd  <- 0.221; label("Additive residual error (mg/L)")         # Table 2: theta6
  })
  model({
    # Covariate effects (multiplicative on CL and V1)
    ffm_cl  <- (FFM / 40.69)^e_ffm_cl
    alb_cl  <- (ALB / 40)^e_alb_cl
    igg_cl  <- (IGG / 14.8)^e_igg_cl
    stdy_cl <- e_stdy_cl^STUDY_LBSL

    ffm_vc  <- (FFM / 40.69)^e_ffm_vc
    stdy_vc <- e_stdy_vc^STUDY_LBSL
    neas_vc <- e_neas_vc^RACE_NEAS
    age_vc  <- AGE / (AGE + age50_vc)

    # PK parameters
    cl <- exp(lcl + etalcl) * ffm_cl * alb_cl * igg_cl * stdy_cl
    vc <- exp(lvc + etalvc) * ffm_vc * stdy_vc * neas_vc * age_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Linear two-compartment IV model (NONMEM ADVAN3 TRANS4)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, vc in L => mg/L
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
