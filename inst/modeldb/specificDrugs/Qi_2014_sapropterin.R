Qi_2014_sapropterin <- function() {
  description <- "One-compartment population PK model with first-order oral absorption, an absorption lag, linear elimination, and an additive endogenous BH4 baseline for sapropterin dihydrochloride in pediatric and adult patients with phenylketonuria (Qi 2014)."
  reference <- "Qi Y, Mould DR, Zhou H, Merilainen M, Musson DG. A prospective population pharmacokinetic analysis of sapropterin dihydrochloride in infants and young children with phenylketonuria. Clinical Pharmacokinetics. 2015;54(2):195-207. doi:10.1007/s40262-014-0196-4"
  vignette <- "Qi_2014_sapropterin"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; constant within an individual in the source dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent clearance and apparent central volume, normalized to a reference weight of 70 kg (a typical adult). Source: Qi 2014 Equation 6.",
      source_name        = "WT"
    ),
    STUDY_PKU015 = list(
      description        = "Indicator for study PKU-015 (pediatric infants and young children 0-6 years) of the Qi 2014 pooled analysis.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (study PKU-004, adolescent / adult >= 9 years)",
      notes              = "Selects the PKU-015 log-scale residual-error magnitude (30.2% CV) under the LTBS approach; STUDY_PKU015 = 0 selects the PKU-004 residual (21.1% CV). Source: Qi 2014 Table 3 reports separate residual error per study.",
      source_name        = "STUDY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 156L,
    n_studies      = 2L,
    age_range      = "0.107-50 years",
    age_mean       = "12 years (SD 11.3)",
    weight_range   = "4.5-144 kg",
    weight_mean    = "40.9 kg (SD 30.3)",
    sex_female_pct = 51.3,
    race_ethnicity = c(Non_Hispanic = 97.4, Hispanic = 2.6),
    disease_state  = "Phenylketonuria (PKU); BH4-responsive subset enrolled in two BioMarin clinical studies.",
    dose_range     = "5 or 20 mg/kg oral sapropterin once daily.",
    regions        = "North America (USA, Canada) and Europe (France, Germany, Ireland, Italy, Poland, UK).",
    baseline_covariates = list(
      bsa_range_m2          = "0.259-2.65",
      bsa_mean_m2           = 1.17,
      alt_range_IU_L        = "11-127",
      alt_mean_IU_L         = 25.4,
      ast_range_IU_L        = "14-63",
      ast_mean_IU_L         = 32.2,
      albumin_range_g_dL    = "3.6-5.0",
      albumin_mean_g_dL     = 4.33,
      crcl_range_mL_min     = "9.39-276",
      crcl_mean_mL_min      = 87.8,
      baseline_phe_umol_L   = "53-2190 (mean 562, SD 378)"
    ),
    notes          = "Baseline demographics from Qi 2014 Table 2. 156 subjects (80 from PKU-015, 78 from PKU-004), 475 plasma concentrations. Age strata in the pooled dataset (Section 3.1): <1 year (n=10), 1-<2 y (n=14), 2-<4 y (n=28), 4-<7 y (n=28), 7-<12 y (n=10), >=12 y (n=66). BH4 plasma concentrations measured indirectly by oxidation to L-biopterin with a validated LC-MS/MS assay; LLOQ = 10.7 ng/mL for BH4 (5.00 ng/mL for L-biopterin) per Section 2.3. Only body weight was retained as a significant covariate; race, sex, age, and laboratory values were screened but not retained (Electronic Supplementary Material Table 1S)."
  )

  ini({
    # Structural parameters -- Qi 2014 Table 3, final population PK model
    # (reference weight 70 kg, typical adult; oral dose, apparent parameters).
    lka   <- log(0.235); label("First-order absorption rate (ka, 1/h)")            # Qi 2014 Table 3 (theta3): ka = 0.235 h^-1
    lcl   <- log(2710);  label("Apparent clearance at WT = 70 kg (CL/F, L/h)")     # Qi 2014 Table 3 (theta1): CL/F = 2710 L/h
    lvc   <- log(3010);  label("Apparent central volume at WT = 70 kg (Vc/F, L)")  # Qi 2014 Table 3 (theta2): Vc/F = 3010 L
    ltlag <- log(0.321); label("Absorption lag time (tlag, h)")                    # Qi 2014 Table 3 (theta4): tlag = 0.321 h
    lc0   <- log(16.6);  label("Endogenous BH4 baseline plasma concentration (C0, ug/L)")  # Qi 2014 Table 3 (theta5): C0 = 16.6 ug/L

    # Body-weight covariate effects: power form normalized to 70 kg (Eq. 6).
    e_wt_cl <- 0.864; label("Power exponent of body weight on CL/F (unitless, reference 70 kg)")  # Qi 2014 Table 3 (theta6)
    e_wt_vc <- 0.644; label("Power exponent of body weight on Vc/F (unitless, reference 70 kg)")  # Qi 2014 Table 3 (theta7)

    # Inter-individual variability. Qi 2014 reports IIV as %CV; the internal
    # variance on the log scale is omega^2 = log(1 + CV^2). CL/F and Vc/F are
    # correlated (Pearson correlation = 0.469 per Table 3).
    #   omega^2(CL) = log(1 + 0.4561^2) = 0.18899
    #   omega^2(Vc) = log(1 + 0.5657^2) = 0.27764
    #   cov(CL,Vc) = 0.469 * sqrt(0.18899 * 0.27764) = 0.10743
    #   omega^2(C0) = log(1 + 0.3647^2) = 0.12487
    etalcl + etalvc ~ c(0.18899, 0.10743, 0.27764)  # Qi 2014 Table 3: IIV CL/F CV = 45.61%, IIV Vc/F CV = 56.57%, corr(CL/F,Vc/F) = 0.469
    etalc0 ~ 0.12487                                # Qi 2014 Table 3: IIV C0 CV = 36.47%

    # Residual error -- log-transform-both-sides (LTBS) constant-CV model
    # (Eq. 5); separate magnitudes per study (Table 3 theta8 and theta9).
    expSdPKU004 <- 0.211; label("Log-scale residual SD, study PKU-004 (LTBS, fraction)")  # Qi 2014 Table 3 (theta8): 21.1% CV
    expSdPKU015 <- 0.302; label("Log-scale residual SD, study PKU-015 (LTBS, fraction)")  # Qi 2014 Table 3 (theta9): 30.2% CV
  })

  model({
    # 1. Derived covariate terms -- power-form WT effects on CL/F and Vc/F
    #    (Qi 2014 Eq. 6, reference weight 70 kg).
    wt_cl <- (WT / 70)^e_wt_cl
    wt_vc <- (WT / 70)^e_wt_vc

    # 2. Individual PK parameters.
    ka   <- exp(lka)
    cl   <- exp(lcl + etalcl) * wt_cl
    vc   <- exp(lvc + etalvc) * wt_vc
    tlag <- exp(ltlag)
    c0   <- exp(lc0 + etalc0)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. ODE system (oral dose into depot; apparent parameters).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Absorption lag on the depot compartment.
    alag(depot) <- tlag

    # 6. Observation: predicted concentration = drug-derived + endogenous BH4.
    #    Dose is in mg and Vc in L, so central/vc is in mg/L; multiply by
    #    1000 to express Cc in ug/L (= ng/mL), the paper's reporting unit.
    Cc <- (central / vc) * 1000 + c0

    # 7. Per-study residual under the LTBS approach: STUDY_PKU015 = 1
    #    selects the pediatric-cohort residual; 0 selects the PKU-004
    #    adolescent / adult residual.
    expSd <- expSdPKU015 * STUDY_PKU015 + expSdPKU004 * (1 - STUDY_PKU015)
    Cc ~ lnorm(expSd)
  })
}
