CirrincioneDall_2011_tacrolimus <- function() {
  description <- "One-compartment population PK model with first-order absorption for oral tacrolimus in pediatric liver transplant recipients (Cirrincione-Dall 2011 ACOP poster, Metrum Research Group). Apparent oral clearance CL/F (25.8 L/h at a 70 kg reference) and apparent volume V/F (2490 L at a 70 kg reference) are estimated; allometric body-weight scaling is fixed at exponent 0.75 on CL/F and 1.0 on V/F. The first-order absorption rate constant ka is fixed at 4.48 1/h from literature because the sparse therapeutic-drug-monitoring sampling could not identify it. CL/F additionally varies (full covariate model) with post-operative day as (POD/7)^0.409, with CYP3A5 expresser status as 1.24^CYP3A5_EXPR (missing genotype data imputed as non-expressers), with AST as (AST/510.5)^-0.0364, with albumin as (ALB/28)^-0.357, with hematocrit as 0.993^HCT (HCT entered as a fraction L/L, not as percent), and with age as (AGE/2)^-0.0310. Inter-individual random variation on CL/F and V/F was modeled exponentially with an estimated covariance of the two random effects per the poster text; the off-diagonal covariance value itself is not reported in the poster Table 2 so this implementation encodes uncorrelated diagonal IIVs and documents the gap in the vignette Errata. Residual error is a combined additive (SD 2.508 ng/mL) + proportional (SD 0.3674 fraction) model on whole-blood tacrolimus concentrations."
  reference <- "Cirrincione-Dall G, Gastonguay MR, Knebel W, Bergsma T, Zhang AY, Patel D, Barrett JS, van Schaik R, Soldin OP, Soldin SJ, Nulman I, Koren G, de Wildt SN. A Population Pharmacokinetic Model of Tacrolimus in Pediatric Liver Transplant Recipients. American Conference on Pharmacometrics (ACOP) 2011 poster, Metrum Research Group, Tariffville CT. https://metrumrg.com/wp-content/uploads/2018/07/acop_2011_tacrolimus.pdf"
  vignette <- "CirrincioneDall_2011_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject across the post-transplant follow-up. Allometric power scaling with reference 70 kg and theory-based fixed exponents: 0.75 on CL/F and 1.0 on V/F (Cirrincione-Dall 2011 Results, paragraph on the PK model). Cohort weight range 2.6-63.6 kg, median 10.6 kg (Table 1).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Centred at 2 years (the cohort median per Table 1) and enters CL/F as (AGE/2)^-0.0310. Cohort range 0.1-15 years.",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying clinical-laboratory covariate. Centred at 28 g/L (the cohort median per Table 1, range 19.4-42.5 g/L) and enters CL/F as (ALB/28)^-0.357. Cirrincione-Dall 2011 Table 1 reports albumin in g/L (cohort mean 29.1 g/L, median 28 g/L, range 19.4-42.5 g/L); the canonical-register ALB unit is g/L, matching the source.",
      source_name        = "ALB"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying clinical-laboratory covariate. Centred at 510.5 U/L (the cohort median per Table 1, range 42-3625 U/L) and enters CL/F as (AST/510.5)^-0.0364. Cohort median markedly elevated reflecting the immediate post-transplant period for liver-graft recipients.",
      source_name        = "AST"
    ),
    HCT = list(
      description        = "Hematocrit, expressed as a fraction of total blood volume (0-1, L/L)",
      units              = "fraction",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Cirrincione-Dall 2011 Table 1 reports HCT as a fraction (cohort mean 0.319, median 0.32, range 0.250-0.440), not as percent (0-100); the canonical-register HCT unit (%) is explicitly overridden here so the parameter value 0.993 reproduces the poster's equation 0.993^HCT directly. To use a dataset that records HCT in percent, multiply the column by 0.01 before passing it to this model.",
      source_name        = "HCT"
    ),
    POD = list(
      description        = "Days post-transplantation",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within subject; rises from 0 at the day of liver transplant surgery. Centred at 7 days (the cohort median per Table 1, range 0-15 days) and enters CL/F as (POD/7)^0.409. Cirrincione-Dall 2011 source name 'POD'. The poster reports CL/F was 55% lower at POD = 1 day and 36% higher at POD = 15 days relative to the POD = 7 day reference (consistent with (1/7)^0.409 = 0.45 and (15/7)^0.409 = 1.34).",
      source_name        = "POD"
    ),
    CYP3A5_EXPR = list(
      description        = "CYP3A5 expresser indicator: 1 if the patient carries at least one functional CYP3A5*1 allele (genotype *1/*1 or *1/*3), 0 if homozygous *3/*3 OR if the CYP3A5 genotype was not determined.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CYP3A5 *3/*3 nonexpresser, or genotype unknown -- Cirrincione-Dall 2011 imputed missing-genotype patients as non-expressers, the predominant category)",
      notes              = "Time-fixed germline genotype determined from rs776746 (CYP3A5 6986A>G). Cohort genotype distribution (n = 34 genotyped of 41 total): CYP3A5 *1/*1 = 1 (2.9%), *1/*3 = 12 (35.3%), *3/*3 = 21 (61.8%); 7 patients (14%) with no available genotype were imputed as non-expressers (the predominant category, per the Table 2 footnote) for the canonical model presented here. Enters CL/F as 1.24^CYP3A5_EXPR, so expressers have a 24% higher apparent oral clearance than *3/*3-or-unknown subjects with the same WT, POD, AST, ALB, HCT, and age. The poster also reports a sensitivity analysis with all missing data imputed as expressers (multiplier 1.14 instead of 1.24, 16.4% RSE); only the predominant-category-imputation case is encoded in this model file. The 95% bootstrap CI on the multiplier (1.05-1.52) does not contain 1.0, so the effect is statistically distinguishable from null in this small cohort.",
      source_name        = "CYP3A"
    )
  )

  population <- list(
    species                = "human",
    n_subjects             = 41L,
    n_studies              = 1L,
    n_observations         = 643L,
    age_range              = "0.1-15 years",
    age_median             = "2 years",
    weight_range           = "2.6-63.6 kg",
    weight_median          = "10.6 kg",
    sex_female_pct         = 53.7,
    race_ethnicity         = "Not reported in poster.",
    disease_state          = "Pediatric liver transplant recipients receiving routine clinical care; the analysis pooled observational therapeutic-drug-monitoring data over the post-transplant follow-up (POD 0-15 days).",
    dose_range             = "Not reported numerically in the poster; routine oral tacrolimus immunosuppression dosed per institutional protocol with concentrations measured by therapeutic drug monitoring (overall observed range 0.42-61.8 ng/mL tacrolimus whole blood).",
    regions                = "Multi-site collaboration: Metrum Research Group (Tariffville CT, USA); The Children's Hospital of Philadelphia (Philadelphia PA, USA); Erasmus MC Sophia Children's Hospital (Rotterdam, Netherlands); Georgetown University (Washington DC, USA); The Hospital for Sick Children (Toronto ON, Canada). Source data: pediatric liver-transplant cohort.",
    samples_per_patient    = "Mean 16, range 7-33 sparse therapeutic-drug-monitoring samples per patient (Table 1).",
    cyp3a5_distribution    = "Of 34 genotyped patients: *1/*1 = 1 (2.9%), *1/*3 = 12 (35.3%), *3/*3 = 21 (61.8%). 7 patients (14% of 41 total) were ungenotyped; the canonical model presented in poster Table 2 imputes ungenotyped patients as non-expressers (the predominant category).",
    notes                  = "Conference poster (ACOP 2011), no DOI. Patient population (Table 1): 19 males and 22 females; weights 2.6-63.6 kg (median 10.6); ages 0.1-15 years (median 2); albumin 19.4-42.5 g/L (median 28); AST 42-3625 U/L (median 510.5); hematocrit 0.250-0.440 (median 0.32); POD 0-15 days (median 7); tacrolimus concentrations 0.42-61.8 ng/mL (mean 15.7). Estimation: NONMEM VII with first-order conditional estimation; non-parametric bootstrap (1000 replicates stratified by age, 266 successful convergences) provided empirical 95% CIs. Covariates included in the full model (no stepwise testing): weight, age, ALB, AST, CYP3A genotype, HCT, and POD."
  )

  ini({
    # Fixed-effect parameter estimates from Cirrincione-Dall 2011 ACOP poster
    # Table 2 (case where missing CYP3A5 genotype data was imputed as
    # non-expressers, per the Table 2 footnote). Reference subject for the
    # typical-value structural parameters: WT = 70 kg adult, CYP3A5 = 0
    # (*3/*3-or-unknown), POD = 7 days, AST = 510.5 U/L, ALB = 28 g/L, HCT = 0
    # (so the 0.993^HCT factor evaluates to 1.0), AGE = 2 years. Note that the
    # 70 kg reference is allometric-theory anchored, not the cohort median (10.6
    # kg) -- typical values reported here therefore correspond to an adult-sized
    # reference subject and shrink for pediatric individuals via the
    # (WT/70)^0.75 factor inside model(). All apparent clearances in L/h;
    # apparent volumes in L; ka in 1/h.

    lka <- fixed(log(4.48)); label("Absorption rate constant ka (1/h) -- fixed from literature")  # Cirrincione-Dall 2011 Table 2 KA = 4.48 1/h, fixed (Results paragraph: "absorption rate constant could not be determined ... and was fixed to 4.48 hr^-1 [5]")
    lcl <- log(25.8)       ; label("Apparent oral clearance CL/F at the reference subject (L/h)") # Cirrincione-Dall 2011 Table 2 theta_1 = 25.8 L/h
    lvc <- log(2490)       ; label("Apparent volume V/F at WT = 70 kg (L)")                       # Cirrincione-Dall 2011 Table 2 theta_2 = 2490 L

    # Allometric exponents fixed by the authors at the theory-based values
    # (Results paragraph: "CL/F and V/F were allometrically scaled by weight
    # (CL/F: exponent fixed at 0.75, V/F: exponent fixed at 1)").
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/70) on CL/F (unitless; fixed at theory value)")  # Cirrincione-Dall 2011 Results -- CL/F exponent fixed at 0.75
    e_wt_vc <- fixed(1.0) ; label("Allometric exponent of (WT/70) on V/F  (unitless; fixed at theory value)")  # Cirrincione-Dall 2011 Results -- V/F  exponent fixed at 1.0

    # Covariate effects on CL/F -- Cirrincione-Dall 2011 Table 2 full-model
    # estimates (case: missing CYP3A5 imputed as non-expressers). The poster's
    # equation for CL/F (linear scale) is:
    #   CL/F = theta_1 * (WT/70)^0.75 * (POD/7)^theta_4 * theta_5^CYP3A *
    #          (AST/510.5)^theta_6 * (ALB/28)^theta_7 * theta_8^HCT *
    #          (AGE/2)^theta_9
    # V/F has no covariate effects other than the fixed allometric weight
    # scaling.
    e_pod_cl     <- 0.409  ; label("Power exponent of (POD/7) on CL/F (unitless)")                              # Cirrincione-Dall 2011 Table 2 theta_4 = 0.409
    e_cyp3a5_cl  <- 1.24   ; label("CYP3A5 expresser multiplicative base on CL/F (theta_5^CYP3A5_EXPR)")        # Cirrincione-Dall 2011 Table 2 theta_5 = 1.24
    e_ast_cl     <- -0.0364; label("Power exponent of (AST/510.5) on CL/F (unitless)")                          # Cirrincione-Dall 2011 Table 2 theta_6 = -0.0364
    e_alb_cl     <- -0.357 ; label("Power exponent of (ALB/28) on CL/F (unitless)")                             # Cirrincione-Dall 2011 Table 2 theta_7 = -0.357
    e_hct_cl     <- 0.993  ; label("Hematocrit multiplicative base on CL/F (theta_8^HCT; HCT as a fraction)")   # Cirrincione-Dall 2011 Table 2 theta_8 = 0.993
    e_age_cl     <- -0.0310; label("Power exponent of (AGE/2) on CL/F (unitless)")                              # Cirrincione-Dall 2011 Table 2 theta_9 = -0.0310

    # Inter-individual variability on CL/F and V/F, modeled exponentially per
    # the poster text. Cirrincione-Dall 2011 Table 2 reports the log-scale
    # variances directly (NONMEM Omega diagonals): Omega_1,1 = 0.554 for CL/F,
    # Omega_2,2 = 0.803 for V/F. The poster text states an inter-eta covariance
    # was estimated ("Inter-individual random variation on CL/F and V/F was
    # modeled exponentially with an estimated covariance of these random
    # effects") but the Omega_2,1 value is NOT reported in Table 2 -- this
    # implementation encodes uncorrelated diagonal IIVs (see vignette Errata).
    etalcl ~ 0.554  # Cirrincione-Dall 2011 Table 2 IIV CL/F (Omega_1,1) = 0.554
    etalvc ~ 0.803  # Cirrincione-Dall 2011 Table 2 IIV V/F  (Omega_2,2) = 0.803

    # Combined additive + proportional residual error on whole-blood tacrolimus
    # concentrations. Cirrincione-Dall 2011 Table 2 reports the NONMEM Sigma
    # diagonals as variances:
    #   Sigma_1,1 = 0.135 (proportional)  -> propSd = sqrt(0.135) = 0.3674
    #   Sigma_2,2 = 6.29  (additive ng/mL^2) -> addSd = sqrt(6.29) = 2.508 ng/mL
    propSd <- 0.3674 ; label("Proportional residual SD (fraction)")        # Cirrincione-Dall 2011 Table 2 prop err (Sigma_1,1) = 0.135 variance; sqrt(0.135) = 0.3674
    addSd  <- 2.508  ; label("Additive residual SD (ng/mL)")                # Cirrincione-Dall 2011 Table 2 add err  (Sigma_2,2) = 6.29 variance ng/mL^2; sqrt(6.29) = 2.508
  })

  model({
    # Body-weight scaling reference: 70 kg adult (Cirrincione-Dall 2011 Results).
    wt70 <- WT / 70
    f_wt_cl <- wt70 ^ e_wt_cl
    f_wt_vc <- wt70 ^ e_wt_vc

    # Covariate factors on CL/F per the poster equation
    # CL/F = theta_1 * (WT/70)^0.75 * (POD/7)^theta_4 * theta_5^CYP3A *
    #        (AST/510.5)^theta_6 * (ALB/28)^theta_7 * theta_8^HCT *
    #        (AGE/2)^theta_9
    f_pod    <- (POD   / 7)     ^ e_pod_cl
    f_cyp3a5 <- e_cyp3a5_cl     ^ CYP3A5_EXPR
    f_ast    <- (AST   / 510.5) ^ e_ast_cl
    f_alb    <- (ALB   / 28)    ^ e_alb_cl
    f_hct    <- e_hct_cl        ^ HCT
    f_age    <- (AGE   / 2)     ^ e_age_cl

    # Individual PK parameters. ka has no IIV in the source (the poster fixes
    # ka from literature because sparse TDM sampling could not identify it).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * f_wt_cl * f_pod * f_cyp3a5 * f_ast * f_alb * f_hct * f_age
    vc <- exp(lvc + etalvc) * f_wt_vc

    # One-compartment oral disposition; bioavailability is implicit in the
    # apparent CL/F and V/F parameterisation (the poster lacks reference IV
    # data, so F cannot be separated from CL or V).
    kel <- cl / vc
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Whole-blood tacrolimus concentration in ng/mL. Dose in mg, vc in L, so
    # central/vc is in mg/L = ug/mL; multiply by 1000 to convert to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
