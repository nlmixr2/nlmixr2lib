Kelman_1984_gentamicin_manchester2 <- function() {
  description <- "One-compartment population PK model for gentamicin in older children and adolescents (Kelman 1984, Manchester II cohort; n=68, age 6 months-15 years, body weight 5.7-62 kg). Clearance is a linear-additive function of body weight, age (in years), and serum creatinine (paper Equation 4, Model 1: CL = theta1*WT + theta2*AGE + theta3*CREAT); volume of distribution is proportional to body weight (Equation 5: V = theta4*WT). Encoded from Kelman 1984 Table 4 Model 1 (FULL model, paper's best by NONMEM objective function for this cohort). For Manchester II the paper's theta2 coefficient is in L/h per year (sign negative), implying age units are years rather than the days used in the Glasgow I and Manchester I cohorts; encoded here with the canonical AGE-in-years covariate."
  reference <- "Kelman AW, Thomson AH, Whiting B, Bryson SM, Steedman DA, Mawer GE, Samba-Donga LA. Estimation of gentamicin clearance and volume of distribution in neonates and young children. Br J Clin Pharmacol. 1984;18:685-692. doi:10.1111/j.1365-2125.1984.tb02530.x"
  vignette <- "Kelman_1984_gentamicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Manchester II cohort range 5.7-62 kg. Linear-additive slope on CL (e_wt_cl = 0.19 L/h per kg) AND sole determinant of V (e_wt_vc = 0.28 L per kg). Paper Table 1.",
      source_name        = "weight"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Manchester II cohort range 6 months to 15 years. Paper Table 4 Model 1 reports theta2 = -0.13 L/h per unit-of-AGE; the magnitude (a few L/h spread across the cohort) is only consistent with AGE in years, not days, since the cohort spans 0.5-15 years. Encoded here on the canonical AGE-in-years scale (see inst/references/covariate-columns.md). Note: the Glasgow I and Manchester I cohorts use PNA in months because their subjects are < 6 months old; the unit shift across cohorts is implicit in the paper.",
      source_name        = "age"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Manchester II cohort range 26-121 umol/L. Linear-additive slope on CL (e_creat_cl = -0.009 L/h per umol/L, paper Table 4 Model 1 theta3). Magnitude is ~10x larger (in absolute terms) than the neonatal Manchester I value, reflecting more developed renal function in the older cohort.",
      source_name        = "serum creatinine"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 68,
    n_studies      = 1,
    age_range      = "6 months to 15 years",
    weight_range   = "5.7-62 kg",
    disease_state  = "Older children and adolescents requiring gentamicin therapy (clinical therapeutic-drug-monitoring data; collected previously by Ashurst et al. 1977 and Samba-Donga 1977)",
    dose_range     = "Gentamicin IV bolus over 1-2 min (n=49 in pooled Manchester I+II; cohort II subset not separately reported) or IM (n=51 in pooled Manchester I+II); dose individualised per clinical practice; steady-state samples after at least 2 days of therapy",
    regions        = "Manchester, England (single centre)",
    notes          = "Manchester II subset (Kelman 1984 Table 1). The Manchester dataset was 'arbitrarily divided into two groups according to age, with the limit set at 6 months' (paper Methods, page 686); Manchester II is the subset with age >= 6 months. 272 total data points (136 peak, 136 trough)."
  )

  ini({
    # Structural anchors -- lcl / lvc fixed to log(1) so the linear-additive slopes below
    # are the primary CL / V drivers in model(); the etalcl / etalvc IIV terms then
    # apply on the derived CL / V via the canonical exp(lcl + etalcl) pattern.
    lcl <- fixed(log(1));    label("Structural log-CL anchor (L/h); CL is constructed in model() from linear slopes")  # see model() block
    lvc <- fixed(log(1));    label("Structural log-V  anchor (L);   V  is constructed in model() from e_wt_vc * WT")  # see model() block

    # Linear-additive covariate slopes (Kelman 1984 Equations 4-5; Table 4 Model 1, Manchester II).
    # The e_<cov>_<param> names are used here for LINEAR-ADDITIVE slope coefficients rather
    # than the canonical multiplicative power-law form (cov^e_<cov>_<param>); this faithfully
    # encodes the paper's Equation 4 (CL = theta1*WT + theta2*AGE + theta3*CREAT) and
    # Equation 5 (V = theta4*WT). See vignette Assumptions and deviations for details.
    e_wt_cl    <- 0.19;      label("Linear slope of CL on WT (L/h per kg)")            # Kelman 1984 Table 4 Model 1 (Manchester II): theta1
    e_age_cl   <- -0.13;     label("Linear slope of CL on AGE (L/h per year)")         # Kelman 1984 Table 4 Model 1 (Manchester II): theta2 (age in years)
    e_creat_cl <- -0.009;    label("Linear slope of CL on CREAT (L/h per umol/L)")     # Kelman 1984 Table 4 Model 1 (Manchester II): theta3
    e_wt_vc    <- 0.28;      label("Linear slope of V on WT (L per kg)")               # Kelman 1984 Table 4 Model 1 (Manchester II): theta4

    # Log-normal IIV on derived CL and V approximating the paper's per-slope additive
    # variances (Equation 6). At typical Manchester II covariates (WT = 20 kg, AGE = 5 years,
    # CREAT = 60 umol/L):
    #   TVCL = 0.19*20 + (-0.13)*5 + (-0.009)*60 = 2.61 L/h
    #   Var(CL|typical) = 400*5.5e-3 + 0 [VAR(eta2) indeterminate] + 3600*3.5e-5
    #                   = 2.2 + 0 + 0.126 = 2.326
    #   SD(CL) approx 1.525 L/h, CV(CL) approx 58%, omega^2(lcl) = log(1+0.584^2) = 0.296
    # V IIV: VAR(theta4) = 6.9e-3 (paper Table 4), so CV(theta4) = sqrt(6.9e-3)/0.28 = 29.6%,
    #        omega^2(lvc) = log(1+0.296^2) = 0.0843.
    # Note: the implied CV on CL is substantial (~58%); this reflects the large
    # paper-reported VAR(theta1) for the older Manchester II cohort, which is converted
    # to log-normal at typical covariates here. Per-individual IIV magnitudes will scale
    # with WT under the paper's true additive structure (heteroscedastic); the log-normal
    # approximation collapses this to homoscedastic IIV. See vignette Errata.
    etalcl ~ 0.296   # log-normal CV approx 58% on CL at typical Manchester II covariates (paper Table 4 VAR(eta_i) collapsed; VAR(eta2) indeterminate)
    etalvc ~ 0.0843  # log-normal CV approx 30% on V (paper Table 4 VAR(theta4) = 6.9e-3 converted)

    # Additive residual error (paper Equation 7; Table 4 reports VAR(eps) = 0.31 (mg/L)^2)
    addSd <- sqrt(0.31);     label("Additive residual SD (mg/L)")                       # Kelman 1984 Table 4 Model 1 (Manchester II): VAR(eps)
  })

  model({
    cl_tv <- e_wt_cl * WT + e_age_cl * AGE + e_creat_cl * CREAT
    v_tv  <- e_wt_vc * WT

    cl <- cl_tv * exp(lcl + etalcl)
    v  <- v_tv  * exp(lvc + etalvc)

    kel <- cl / v

    d/dt(central) <- -kel * central

    Cc <- central / v
    Cc ~ add(addSd)
  })
}
