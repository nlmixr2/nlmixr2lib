Kelman_1984_gentamicin_manchester1 <- function() {
  description <- "One-compartment population PK model for gentamicin in neonates and young infants (Kelman 1984, Manchester I cohort; n=32, postnatal age 1-153 days, body weight 1.6-9.1 kg). Clearance is a linear-additive function of body weight, postnatal age, and serum creatinine (paper Equation 4, Model 1: CL = theta1*WT + theta2*AGE + theta3*CREAT); volume of distribution is proportional to body weight (Equation 5: V = theta4*WT). Encoded from Kelman 1984 Table 3 Model 1 (FULL model, paper's best by NONMEM objective function for this cohort)."
  reference <- "Kelman AW, Thomson AH, Whiting B, Bryson SM, Steedman DA, Mawer GE, Samba-Donga LA. Estimation of gentamicin clearance and volume of distribution in neonates and young children. Br J Clin Pharmacol. 1984;18:685-692. doi:10.1111/j.1365-2125.1984.tb02530.x"
  vignette <- "Kelman_1984_gentamicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Manchester I cohort range 1.6-9.1 kg. Linear-additive slope on CL (e_wt_cl = 0.089 L/h per kg) AND sole determinant of V (e_wt_vc = 0.41 L per kg). Paper Table 1.",
      source_name        = "weight"
    ),
    PNA = list(
      description        = "Postnatal age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Manchester I cohort range 1-153 days (approx 0.033-5.03 months). Paper Table 3 Model 1 reports theta2 = 0.0013 L/h per day-of-age; reparameterised here as e_pna_cl = 0.0013 * 30.4375 = 0.03957 L/h per month-of-PNA so the covariate column matches the canonical PNA-in-months scale (see inst/references/covariate-columns.md).",
      source_name        = "age (paper reports in days; converted to canonical PNA in months for covariate column)"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Manchester I cohort range 26-165 umol/L. Linear-additive slope on CL (e_creat_cl = -0.0010 L/h per umol/L, paper Table 3 Model 1 theta3).",
      source_name        = "serum creatinine"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32,
    n_studies      = 1,
    age_range      = "1-153 days postnatal",
    weight_range   = "1.6-9.1 kg",
    disease_state  = "Neonates and young infants requiring gentamicin therapy (clinical therapeutic-drug-monitoring data; collected previously by Ashurst et al. 1977 and Samba-Donga 1977)",
    dose_range     = "Gentamicin IV bolus over 1-2 min (n=49 in pooled Manchester I+II; cohort I subset not separately reported) or IM (n=51 in pooled Manchester I+II); dose individualised per clinical practice; steady-state samples after at least 2 days of therapy",
    regions        = "Manchester, England (single centre)",
    notes          = "Manchester I subset (Kelman 1984 Table 1). The Manchester dataset was 'arbitrarily divided into two groups according to age, with the limit set at 6 months' (paper Methods, page 686); Manchester I is the subset with PNA < 6 months. 128 total data points (peak/trough split not separately reported by subset in Table 1)."
  )

  ini({
    # Structural anchors -- lcl / lvc fixed to log(1) so the linear-additive slopes below
    # are the primary CL / V drivers in model(); the etalcl / etalvc IIV terms then
    # apply on the derived CL / V via the canonical exp(lcl + etalcl) pattern.
    lcl <- fixed(log(1));    label("Structural log-CL anchor (L/h); CL is constructed in model() from linear slopes")  # see model() block
    lvc <- fixed(log(1));    label("Structural log-V  anchor (L);   V  is constructed in model() from e_wt_vc * WT")  # see model() block

    # Linear-additive covariate slopes (Kelman 1984 Equations 4-5; Table 3 Model 1, Manchester I).
    # The e_<cov>_<param> names are used here for LINEAR-ADDITIVE slope coefficients rather
    # than the canonical multiplicative power-law form (cov^e_<cov>_<param>); this faithfully
    # encodes the paper's Equation 4 (CL = theta1*WT + theta2*AGE + theta3*CREAT) and
    # Equation 5 (V = theta4*WT). See vignette Assumptions and deviations for details.
    e_wt_cl    <- 0.089;     label("Linear slope of CL on WT (L/h per kg)")            # Kelman 1984 Table 3 Model 1 (Manchester I): theta1
    e_pna_cl   <- 0.03957;   label("Linear slope of CL on PNA (L/h per month-of-PNA)") # Kelman 1984 Table 3 Model 1 theta2 = 0.0013 L/h/day * 30.4375 day/month
    e_creat_cl <- -0.001;    label("Linear slope of CL on CREAT (L/h per umol/L)")     # Kelman 1984 Table 3 Model 1 (Manchester I): theta3
    e_wt_vc    <- 0.41;      label("Linear slope of V on WT (L per kg)")               # Kelman 1984 Table 3 Model 1 (Manchester I): theta4

    # Log-normal IIV on derived CL and V approximating the paper's per-slope additive
    # variances (Equation 6). At typical Manchester I covariates (WT = 4 kg, PNA = 2 months,
    # CREAT = 60 umol/L):
    #   TVCL = 0.089*4 + 0.03957*2 + (-0.001)*60 = 0.375 L/h
    #   Var(CL|typical) = 16*2.4e-5 + 4*(6.0e-6*30.4375^2) + 0 [VAR(eta3) indeterminate]
    #                   = 3.84e-4 + 0.02224 = 0.02263
    #   SD(CL) approx 0.150 L/h, CV(CL) approx 40%, omega^2(lcl) = log(1+0.40^2) = 0.148
    # V IIV: VAR(theta4) = 1.6e-2 (paper Table 3), so CV(theta4) = sqrt(1.6e-2)/0.41 = 30.7%,
    #        omega^2(lvc) = log(1+0.307^2) = 0.0902.
    etalcl ~ 0.148   # log-normal CV approx 40% on CL at typical Manchester I covariates (paper Table 3 VAR(eta_i) collapsed; VAR(eta3) indeterminate)
    etalvc ~ 0.0902  # log-normal CV approx 31% on V (paper Table 3 VAR(theta4) = 1.6e-2 converted)

    # Additive residual error (paper Equation 7; Table 3 reports VAR(eps) = 0.53 (mg/L)^2)
    addSd <- sqrt(0.53);     label("Additive residual SD (mg/L)")                       # Kelman 1984 Table 3 Model 1 (Manchester I): VAR(eps)
  })

  model({
    cl_tv <- e_wt_cl * WT + e_pna_cl * PNA + e_creat_cl * CREAT
    v_tv  <- e_wt_vc * WT

    cl <- cl_tv * exp(lcl + etalcl)
    v  <- v_tv  * exp(lvc + etalvc)

    kel <- cl / v

    d/dt(central) <- -kel * central

    Cc <- central / v
    Cc ~ add(addSd)
  })
}
