Wickramasinghe_2025_pamiparib <- function() {
  description <- "One-compartment population PK model with first-order absorption (no lag time) and first-order elimination describing total and unbound plasma pamiparib, a PARP1/2 inhibitor, in 41 adults with newly diagnosed or recurrent glioblastoma receiving 60 mg orally twice daily (Wickramasinghe 2025). Total and unbound plasma concentrations were fitted simultaneously in Monolix 2024R1 with the unbound concentration derived algebraically as Cu = Cc * fu. This file encodes the total-drug parameterisation (Table 2 'Total Pamiparib' final covariate model), in which the apparent volume and apparent clearance are referenced to the total plasma concentration: V/F = 15 * exp(0.0094 * CRCL) L and CL/F = 6.76 * exp(-0.016 * AGE) L/h, giving typical values of 44 L and 2.59 L/h at the population median post-operative creatinine clearance (111.5 mL/min) and median age (60 years). Post-operative creatinine clearance on V/F and age on CL/F are the only retained covariates, explaining about 23% and 6% of the respective inter-individual variability. Log-normal IIV is carried on Ka (410% CV), V/F (41%), CL/F (50%) and fu (12%). A companion file Wickramasinghe_2025_pamiparib_unbound encodes the same joint fit reparameterised on the unbound concentration. The combined proportional-plus-additive residual error model the authors selected is declared but its magnitudes are not reported anywhere in the paper, so both residual SDs are carried at zero; see the vignette Errata."
  reference <- paste(
    "Wickramasinghe C, Kim S, Jiang Y, Bao X, Yue Y, Jiang J, Hong A, Sanai N, Li J.",
    "Population Pharmacokinetic Modeling of Total and Unbound Pamiparib in Glioblastoma Patients:",
    "Insights into Drug Disposition and Dosing Optimization.",
    "Pharmaceutics. 2025;17(4):524.",
    "doi:10.3390/pharmaceutics17040524.",
    sep = " "
  )
  vignette <- "Wickramasinghe_2025_pamiparib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "pamiparib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pamiparib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Post-operative creatinine clearance, measured after surgical tumour resection. Raw creatinine clearance in mL/min, NOT normalised to 1.73 m^2 body surface area. The paper reports post-operative creatinine clearance (its column 'PCC') separately from pre-dose creatinine clearance and from CKD-EPI-estimated GFR; the final covariate model uses the post-operative creatinine clearance value.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters V/F through an UNCENTERED exponential effect V/F = theta2 * exp(e_crcl_vc * CRCL), so exp(lvc) is the extrapolated V/F at CRCL = 0 rather than a typical value. Table 2 footnote gives the population median post-operative creatinine clearance as 111.5 mL/min, at which the typical V/F is 44 L; Table 1 reports the post-operative creatinine clearance median as 111 mL/min (range 59-169). Section 2.2.2 states that continuous covariates were centered on their medians, but the Table 2 footnote equation and the reported theta2 = 15 L are only mutually consistent in the uncentered form, so the uncentered form is encoded here. Eligibility required eGFR >= 30 mL/min/1.73 m^2 (CKD-EPI) and either serum creatinine <= 1.5 x ULN or creatinine clearance >= 60 mL/min, so the cohort excludes moderate-to-severe renal impairment.",
      source_name        = "PCC"
    ),
    AGE = list(
      description        = "Age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters CL/F through an UNCENTERED exponential effect CL/F = theta3 * exp(e_age_cl * AGE), so exp(lcl) is the extrapolated CL/F at age 0 rather than a typical value. Table 2 footnote gives the population median age as 60 years, at which the typical CL/F is 2.59 L/h. Table 1 reports age median 60 years (range 31-80); eligibility required age > 18 years.",
      source_name        = "Age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "31-80 years",
    age_median     = "60 years",
    weight_range   = "45-129 kg",
    weight_median  = "80 kg",
    height_range   = "155-193 cm",
    height_median  = "173 cm",
    bsa_range      = "1.41-2.53 m^2",
    bsa_median     = "1.99 m^2",
    sex_female_pct = 51.2,
    race_ethnicity = c(White = 92.7, `Non-white` = 7.3),
    disease_state  = "Newly diagnosed or recurrent glioblastoma; ECOG performance status <= 2; adequate bone marrow, hepatic and renal function.",
    dose_range     = "60 mg oral pamiparib twice daily for 3 days (9 doses) before surgical tumour resection. Plasma sampled pre-dose and at 0.5, 1, 2, 4, 7 and 24 h after the 9th dose.",
    regions        = "United States (Barrow Neurological Institute, St. Joseph's Hospital and Medical Center, Phoenix AZ; single centre).",
    renal_function = "Post-operative creatinine clearance median 111 mL/min (range 59-169); pre-dose median 98 mL/min (range 39-154). Post-operative GFR median 101 mL/min (range 48-117).",
    hepatic_function = "Pre-dose total bilirubin median 0.5 mg/dL (range 0.3-1.8); pre-dose AST median 19 IU/L, ALT median 24 IU/L; pre-dose plasma albumin median 4.1 mg/dL.",
    co_medication  = "36 of 41 patients received co-administered drugs during the trial; 34 of 41 received dexamethasone (median total dose 19 mg, range 0-66). No concomitant medication was retained as a significant covariate.",
    notes          = "Phase 0/II study NCT04614909. Baseline demographics from Table 1. Estimation by SAEM in Monolix 2024R1 on log-transformed concentrations; covariate screening by generalized additive modelling followed by stepwise covariate modelling (forward dOFV > 3.875, backward dOFV > 10.828). Parameter precision assessed by 500-replicate nonparametric bootstrap (Table 3)."
  )

  ini({
    # ============================================================
    # Structural parameters - Table 2, column "Total Pamiparib /
    # Final Model". The theta rows (theta1..theta4) are the model
    # coefficients; the TV_ rows above them are those coefficients
    # evaluated at the population median covariate values, per the
    # Table 2 footnotes. Because both covariate effects are
    # UNCENTERED exponentials, exp(lvc) and exp(lcl) are the
    # covariate-free intercepts theta2 and theta3, not typical
    # values.
    # ============================================================
    lka <- log(1.58);  label("Absorption rate constant Ka (1/h)")                                             # Table 2, Total Pamiparib final model, theta1 (Ka) = 1.58 (RSE 42%); Table 3 bootstrap mean 1.71 (95% CI 0.64-4.15)
    lvc <- log(15);    label("Apparent volume of distribution intercept theta2, V/F at CRCL = 0 (L)")         # Table 2, Total Pamiparib final model, theta2 (V/F) = 15 (RSE 28%); Table 3 bootstrap mean 17 (95% CI 8-30)
    lcl <- log(6.76);  label("Apparent oral clearance intercept theta3, CL/F at age 0 (L/h)")                 # Table 2, Total Pamiparib final model, theta3 (CL/F) = 6.76 (RSE 38%); Table 3 bootstrap mean 7.76 (95% CI 3.6-15.9)
    lfu <- log(0.041); label("Fraction of pamiparib unbound in plasma fu (unitless)")                         # Table 2, Total Pamiparib final model, theta4 (Fu) = 0.041 (RSE 2.8%); Table 3 bootstrap mean 0.041 (95% CI 0.039-0.044)

    # ============================================================
    # Covariate effects - Table 2, beta rows, and the Table 2
    # footnotes that give the functional form:
    #   TV_V  = theta2 * exp(beta1 * PCC),  median PCC  = 111.5 mL/min
    #   TV_CL = theta3 * exp(beta2 * Age),  median Age  = 60 years
    # Check: 15 * exp(0.0094 * 111.5)  = 42.8 L   (Table 2 TV_V/F  = 44)
    #        6.76 * exp(-0.016 * 60)   = 2.588 L/h (Table 2 TV_CL/F = 2.59)
    # The 42.8-versus-44 gap is a rounding artefact of the 2-significant-
    # figure theta2 = 15; it is NOT tuned away here. See vignette Errata.
    # ============================================================
    e_crcl_vc <- 0.0094; label("Exponential effect of post-operative creatinine clearance on V/F (per mL/min)")  # Table 2, Total Pamiparib final model, beta1 (PCC on V/F) = 0.0094 (RSE 25%); Table 3 bootstrap mean 0.0089 (95% CI 0.0035-0.015)
    e_age_cl  <- -0.016; label("Exponential effect of age on CL/F (per year)")                                   # Table 2, Total Pamiparib final model, beta2 (Age on CL/F) = -0.016 (RSE 37%); Table 3 bootstrap mean -0.017 (95% CI -0.03 to -0.0054)

    # ============================================================
    # Inter-individual variability. Section 2.2.1 states that IIV
    # was modelled with an exponential (log-normal) function on all
    # PK parameters and that the random effects are uncorrelated
    # (no covariance block is reported). Table 3 reports the
    # log-scale SDs of eta directly (Ka_SD, V/F_SD, CL/F_SD,
    # Fu_SD); nlmixr2 needs the variance, so each entry below is
    # that SD squared. Cross-check against the CV% figures in
    # Table 2 via CV = sqrt(exp(SD^2) - 1):
    #   SD 1.7  -> 412% (Table 2 IIV of Ka   = 410%)
    #   SD 0.39 ->  40.5% (Table 2 IIV of V/F  = 41%)
    #   SD 0.47 ->  49.7% (Table 2 IIV of CL/F = 50%)
    #   SD 0.12 ->  12.0% (Table 2 IIV of Fu   = 12%)
    # ============================================================
    etalka ~ 2.8900   # Table 3, Total Pamiparib, Ka_SD   = 1.7  (bootstrap 1.6, 95% CI 1.0-2.3)  -> variance 1.7^2
    etalvc ~ 0.1521   # Table 3, Total Pamiparib, V/F_SD  = 0.39 (bootstrap 0.35, 95% CI 0.16-0.5) -> variance 0.39^2
    etalcl ~ 0.2209   # Table 3, Total Pamiparib, CL/F_SD = 0.47 (bootstrap 0.46, 95% CI 0.35-0.56) -> variance 0.47^2
    etalfu ~ 0.0144   # Table 3, Total Pamiparib, Fu_SD   = 0.12 (bootstrap 0.11, 95% CI 0.06-0.16) -> variance 0.12^2

    # ============================================================
    # Residual error. Section 2.2.1 selected a combined error model
    # g = a + b * f with c fixed at 1, applied to both the total
    # and the unbound observations. The magnitudes a and b are
    # reported NOWHERE in the paper: Table 2 lists only the thetas,
    # betas and IIV percentages, and Table 3 lists only the eta
    # SDs. Carrying zero here records that the authors declared a
    # combined error structure while making clear that no sigma was
    # published; zero does NOT mean the authors estimated zero
    # residual variability. Simulations from this model reproduce
    # the structural and between-subject spread only. See the
    # vignette Errata.
    # ============================================================
    addSd     <- fixed(0); label("Additive residual SD on total plasma Cc (mg/L; not reported by the source)")           # Section 2.2.1 Equation (1): combined error g = a + b*f; a not reported
    propSd    <- fixed(0); label("Proportional residual SD on total plasma Cc (fraction; not reported by the source)")   # Section 2.2.1 Equation (1): combined error g = a + b*f; b not reported
    addSd_Cu  <- fixed(0); label("Additive residual SD on unbound plasma Cu (mg/L; not reported by the source)")         # Section 2.2.1 Equation (1) applied to the unbound observations; a not reported
    propSd_Cu <- fixed(0); label("Proportional residual SD on unbound plasma Cu (fraction; not reported by the source)") # Section 2.2.1 Equation (1) applied to the unbound observations; b not reported
  })

  model({
    # ---- Individual parameters -----------------------------------
    # Both covariate effects are uncentered exponentials, exactly as
    # printed in the Table 2 footnotes.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc + e_crcl_vc * CRCL)
    cl <- exp(lcl + etalcl + e_age_cl * AGE)
    fu <- exp(lfu + etalfu)

    kel <- cl / vc

    # ---- ODE system (Section 2.2.1) -------------------------------
    # One compartment, first-order absorption without lag time,
    # first-order elimination. Doses are oral and go to depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observations ---------------------------------------------
    # Cc is the total plasma concentration; Cu is the unbound plasma
    # concentration derived through the fixed binding relationship
    # Cu = Cc * fu (Section 2.2.1). Both were fitted simultaneously.
    Cc <- central / vc
    Cu <- fu * Cc

    Cc ~ add(addSd) + prop(propSd)
    Cu ~ add(addSd_Cu) + prop(propSd_Cu)
  })
}
