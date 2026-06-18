Kelman_1984_gentamicin_glasgow <- function() {
  description <- "One-compartment population PK model for gentamicin in neonates and very young infants (Kelman 1984, Glasgow I cohort; n=43, postnatal age 2-120 days, body weight 0.8-3.7 kg). Clearance is a linear-additive function of body weight, postnatal age, and serum creatinine (paper Equation 4, Model 1: CL = theta1*WT + theta2*AGE + theta3*CREAT); volume of distribution is proportional to body weight (Equation 5: V = theta4*WT). Encoded from Kelman 1984 Table 2 Model 1 (FULL model, paper's best by NONMEM objective function)."
  reference <- "Kelman AW, Thomson AH, Whiting B, Bryson SM, Steedman DA, Mawer GE, Samba-Donga LA. Estimation of gentamicin clearance and volume of distribution in neonates and young children. Br J Clin Pharmacol. 1984;18:685-692. doi:10.1111/j.1365-2125.1984.tb02530.x"
  vignette <- "Kelman_1984_gentamicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Glasgow I cohort range 0.8-3.7 kg. Linear-additive slope on CL (e_wt_cl = 0.057 L/h per kg) AND sole determinant of V (e_wt_vc = 0.46 L per kg). Paper Table 1 records weight at the time of gentamicin sampling.",
      source_name        = "weight"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Glasgow I cohort range 2-120 days (approx 0.066-3.94 months). Paper Table 2 Model 1 reports theta2 = 0.00074 L/h per day-of-age; reparameterised here as e_pna_cl = 0.00074 * 30.4375 = 0.02252 L/h per month-of-PNA so the covariate column matches the canonical PNA-in-months scale (see inst/references/covariate-columns.md).",
      source_name        = "age (paper reports in days; converted to canonical PNA in months for covariate column)"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Glasgow I cohort range 39-120 umol/L. Linear-additive slope on CL (e_creat_cl = -0.00019 L/h per umol/L, paper Table 2 Model 1 theta3). The paper notes Glasgow I had no subjects with extreme creatinine values, so the magnitude of the CREAT effect is modest after adjusting for weight.",
      source_name        = "serum creatinine"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 43,
    n_studies      = 1,
    age_range      = "2-120 days postnatal",
    weight_range   = "0.8-3.7 kg",
    disease_state  = "Neonates and very young infants requiring gentamicin therapy (clinical therapeutic-drug-monitoring data)",
    dose_range     = "Gentamicin IV bolus over 1-2 min (n=39) or IM (n=4); dose individualised per clinical practice; steady-state samples after at least 2 days of therapy",
    regions        = "Glasgow, Scotland (single centre)",
    notes          = "Glasgow I cohort (Kelman 1984 Table 1). Estimated gestational age 26-39 weeks. 82 total data points (36 peak, 46 trough). Used as training data for the NONMEM Model 1 fit reported in Table 2; Glasgow II (n=26, 50 datapoints) was held out for prospective validation."
  )

  ini({
    # Structural anchors -- lcl / lvc fixed to log(1) so the linear-additive slopes below
    # are the primary CL / V drivers in model(); the etalcl / etalvc IIV terms then
    # apply on the derived CL / V via the canonical exp(lcl + etalcl) pattern.
    lcl <- fixed(log(1));    label("Structural log-CL anchor (L/h); CL is constructed in model() from linear slopes")  # see model() block
    lvc <- fixed(log(1));    label("Structural log-V  anchor (L);   V  is constructed in model() from e_wt_vc * WT")  # see model() block

    # Linear-additive covariate slopes (Kelman 1984 Equations 4-5; Table 2 Model 1, Glasgow I).
    # The e_<cov>_<param> names are used here for LINEAR-ADDITIVE slope coefficients rather
    # than the canonical multiplicative power-law form (cov^e_<cov>_<param>); this faithfully
    # encodes the paper's Equation 4 (CL = theta1*WT + theta2*AGE + theta3*CREAT) and
    # Equation 5 (V = theta4*WT). See vignette Assumptions and deviations for details.
    e_wt_cl    <- 0.057;     label("Linear slope of CL on WT (L/h per kg)")            # Kelman 1984 Table 2 Model 1 (Glasgow I): theta1
    e_pna_cl   <- 0.02252;   label("Linear slope of CL on PNA (L/h per month-of-PNA)") # Kelman 1984 Table 2 Model 1 theta2 = 0.00074 L/h/day * 30.4375 day/month
    e_creat_cl <- -0.00019;  label("Linear slope of CL on CREAT (L/h per umol/L)")     # Kelman 1984 Table 2 Model 1 (Glasgow I): theta3
    e_wt_vc    <- 0.46;      label("Linear slope of V on WT (L per kg)")               # Kelman 1984 Table 2 Model 1 (Glasgow I): theta4

    # Log-normal IIV on derived CL and V approximating the paper's per-slope additive
    # variances (Equation 6: theta_i,r = theta_i + eta_i,r). At typical Glasgow I covariates
    # (WT = 2 kg, PNA = 1 month, CREAT = 60 umol/L):
    #   TVCL = 0.057*2 + 0.02252*1 + (-0.00019)*60 = 0.125 L/h
    #   Var(CL|typical) = WT^2*VAR(eta1) + PNA^2*VAR(eta_pna) + CREAT^2*VAR(eta3)
    #                   = 4 * 2.5e-5 + 1 * (6.6e-7 * 30.4375^2) + 3600 * 1.1e-8
    #                   = 1.0e-4 + 6.12e-4 + 3.96e-5 = 7.51e-4
    #   SD(CL) approx 0.0274 L/h, CV(CL) approx 21.9%
    #   omega^2(lcl) = log(1 + 0.219^2) = 0.0469
    # Paper VAR(theta4) is reported as indeterminate ("-") in Glasgow I; use a small
    # placeholder for V IIV so the model is simulable.
    etalcl ~ 0.0469  # log-normal CV approx 22% on CL at typical Glasgow I covariates (paper Table 2 VAR(eta_i) collapsed)
    etalvc ~ 0.01    # placeholder; paper Table 2 VAR(theta4) indeterminate in Glasgow I

    # Additive residual error (paper Equation 7; Table 2 reports VAR(eps) = 0.88 (mg/L)^2)
    addSd <- sqrt(0.88);     label("Additive residual SD (mg/L)")                       # Kelman 1984 Table 2 Model 1 (Glasgow I): VAR(eps)
  })

  model({
    # Linear-additive typical CL (paper Equation 4) and V (Equation 5)
    cl_tv <- e_wt_cl * WT + e_pna_cl * PNA + e_creat_cl * CREAT
    v_tv  <- e_wt_vc * WT

    # Log-normal IIV on derived CL and V (approximation of paper Equation 6 per-slope
    # additive IIV; see vignette Assumptions and deviations). exp(lcl) and exp(lvc)
    # both evaluate to 1 (lcl = lvc = log(1) = 0, FIXED), so the linear-additive slopes
    # above are the structural drivers and the etas apply log-normally on top.
    cl <- cl_tv * exp(lcl + etalcl)
    v  <- v_tv  * exp(lvc + etalvc)

    kel <- cl / v

    # One-compartment model. IV bolus (n=39, 1-2 min) and rapid IM absorption (n=4)
    # both modelled as instantaneous input to central per paper Methods (page 686:
    # "assuming rapid intramuscular absorption").
    d/dt(central) <- -kel * central

    Cc <- central / v
    Cc ~ add(addSd)
  })
}
