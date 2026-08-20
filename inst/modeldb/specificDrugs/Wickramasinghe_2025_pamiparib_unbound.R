Wickramasinghe_2025_pamiparib_unbound <- function() {
  description <- "One-compartment population PK model with first-order absorption (no lag time) and first-order elimination describing total and unbound plasma pamiparib, a PARP1/2 inhibitor, in 41 adults with newly diagnosed or recurrent glioblastoma receiving 60 mg orally twice daily (Wickramasinghe 2025). This file encodes the unbound-drug parameterisation (Table 2 'Unbound Pamiparib' final covariate model), in which the apparent volume and apparent clearance are referenced to the UNBOUND plasma concentration: Vu/F = 402 * exp(0.0087 * CRCL) L and CLu/F = 163 * exp(-0.016 * AGE) L/h, giving typical values of 1060 L and 62.5 L/h at the population median post-operative creatinine clearance (111.5 mL/min) and median age (60 years). The unbound concentration is the primary prediction (Cu = central / Vu) and the total plasma concentration is recovered as Cc = Cu / fu. This is the authors' second, separately estimated fit of the same simultaneously modelled total-plus-unbound dataset (OFV 9212 versus 9210); it is close to but not exactly the algebraic rescaling of the companion total-drug file Wickramasinghe_2025_pamiparib by fu, since each run was estimated independently by SAEM. Log-normal IIV is carried on Ka (380% CV), Vu/F (41%), CLu/F (52%) and fu (11%). The combined proportional-plus-additive residual error model the authors selected is declared but its magnitudes are not reported anywhere in the paper, so both residual SDs are carried at zero; see the vignette Errata."
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
      notes              = "Time-fixed per subject. Enters the unbound-referenced apparent volume through an UNCENTERED exponential effect Vu/F = theta2 * exp(e_crcl_vc * CRCL), so exp(lvc) is the extrapolated Vu/F at CRCL = 0 rather than a typical value. Table 2 footnote gives the population median post-operative creatinine clearance as 111.5 mL/min, at which the typical Vu/F is 1060 L; Table 1 reports the post-operative creatinine clearance median as 111 mL/min (range 59-169). Section 2.2.2 states that continuous covariates were centered on their medians, but the Table 2 footnote equation and the reported theta2 = 402 L are only mutually consistent in the uncentered form, so the uncentered form is encoded here. Eligibility required eGFR >= 30 mL/min/1.73 m^2 (CKD-EPI) and either serum creatinine <= 1.5 x ULN or creatinine clearance >= 60 mL/min, so the cohort excludes moderate-to-severe renal impairment.",
      source_name        = "PCC"
    ),
    AGE = list(
      description        = "Age at study entry",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters the unbound-referenced apparent clearance through an UNCENTERED exponential effect CLu/F = theta3 * exp(e_age_cl * AGE), so exp(lcl) is the extrapolated CLu/F at age 0 rather than a typical value. Table 2 footnote gives the population median age as 60 years, at which the typical CLu/F is 62.5 L/h. Table 1 reports age median 60 years (range 31-80); eligibility required age > 18 years.",
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
    # Structural parameters - Table 2, column "Unbound Pamiparib /
    # Final Model". Volume and clearance here are referenced to the
    # unbound plasma concentration, so they are larger than their
    # total-drug counterparts by roughly 1/fu. Both covariate
    # effects are UNCENTERED exponentials, so exp(lvc) and exp(lcl)
    # are the covariate-free intercepts theta2 and theta3, not
    # typical values.
    # ============================================================
    lka <- log(1.7);   label("Absorption rate constant Ka (1/h)")                                                     # Table 2, Unbound Pamiparib final model, theta1 (Ka) = 1.7 (RSE 44%); Table 3 bootstrap mean 1.86 (95% CI 0.67-4.07)
    lvc <- log(402);   label("Unbound-referenced apparent volume intercept theta2, Vu/F at CRCL = 0 (L)")             # Table 2, Unbound Pamiparib final model, theta2 (V/F) = 402 (RSE 28%); Table 3 bootstrap mean 441 (95% CI 219-787)
    lcl <- log(163);   label("Unbound-referenced apparent clearance intercept theta3, CLu/F at age 0 (L/h)")          # Table 2, Unbound Pamiparib final model, theta3 (CL/F) = 163 (RSE 39%); Table 3 bootstrap mean 185 (95% CI 83-454)
    lfu <- log(0.042); label("Fraction of pamiparib unbound in plasma fu (unitless)")                                 # Table 2, Unbound Pamiparib final model, theta4 (Fu) = 0.042 (RSE 2.6%); Table 3 bootstrap mean 0.041 (95% CI 0.039-0.044)

    # ============================================================
    # Covariate effects - Table 2, beta rows, with the functional
    # form from the Table 2 footnotes:
    #   TV_V  = theta2 * exp(beta1 * PCC),  median PCC = 111.5 mL/min
    #   TV_CL = theta3 * exp(beta2 * Age),  median Age = 60 years
    # Check: 402 * exp(0.0087 * 111.5) = 1060.4 L   (Table 2 TV_V/F  = 1060)
    #        163 * exp(-0.016 * 60)    =   62.4 L/h (Table 2 TV_CL/F =   62.5)
    # ============================================================
    e_crcl_vc <- 0.0087; label("Exponential effect of post-operative creatinine clearance on Vu/F (per mL/min)")  # Table 2, Unbound Pamiparib final model, beta1 (PCC on V/F) = 0.0087 (RSE 26%); Table 3 bootstrap mean 0.0083 (95% CI 0.0027-0.014)
    e_age_cl  <- -0.016; label("Exponential effect of age on CLu/F (per year)")                                   # Table 2, Unbound Pamiparib final model, beta2 (Age on CL/F) = -0.016 (RSE 38%); Table 3 bootstrap mean -0.016 (95% CI -0.032 to -0.0053)

    # ============================================================
    # Inter-individual variability. Table 3 reports the log-scale
    # SDs of eta directly; nlmixr2 needs the variance, so each
    # entry is that SD squared. Cross-check against the Table 2
    # CV% figures via CV = sqrt(exp(SD^2) - 1):
    #   SD 1.66 -> 384% (Table 2 IIV of Ka   = 380%)
    #   SD 0.39 ->  40.5% (Table 2 IIV of V/F  = 41%)
    #   SD 0.49 ->  52.1% (Table 2 IIV of CL/F = 52%)
    #   SD 0.11 ->  11.0% (Table 2 IIV of Fu   = 11%)
    # ============================================================
    etalka ~ 2.7556   # Table 3, Unbound Pamiparib, Ka_SD   = 1.66 (bootstrap 1.66, 95% CI 1.1-2.34) -> variance 1.66^2
    etalvc ~ 0.1521   # Table 3, Unbound Pamiparib, V/F_SD  = 0.39 (bootstrap 0.35, 95% CI 0.16-0.5) -> variance 0.39^2
    etalcl ~ 0.2401   # Table 3, Unbound Pamiparib, CL/F_SD = 0.49 (bootstrap 0.46) -> variance 0.49^2
    etalfu ~ 0.0121   # Table 3, Unbound Pamiparib, Fu_SD   = 0.11 (bootstrap 0.11, 95% CI 0.059-0.15) -> variance 0.11^2

    # ============================================================
    # Residual error. Section 2.2.1 selected a combined error model
    # g = a + b * f with c fixed at 1, applied to both the total
    # and the unbound observations. The magnitudes a and b are
    # reported NOWHERE in the paper. Carrying zero here records
    # that the authors declared a combined error structure while
    # making clear that no sigma was published; zero does NOT mean
    # the authors estimated zero residual variability. See the
    # vignette Errata.
    # ============================================================
    addSd     <- fixed(0); label("Additive residual SD on total plasma Cc (mg/L; not reported by the source)")           # Section 2.2.1 Equation (1): combined error g = a + b*f; a not reported
    propSd    <- fixed(0); label("Proportional residual SD on total plasma Cc (fraction; not reported by the source)")   # Section 2.2.1 Equation (1): combined error g = a + b*f; b not reported
    addSd_Cu  <- fixed(0); label("Additive residual SD on unbound plasma Cu (mg/L; not reported by the source)")         # Section 2.2.1 Equation (1) applied to the unbound observations; a not reported
    propSd_Cu <- fixed(0); label("Proportional residual SD on unbound plasma Cu (fraction; not reported by the source)") # Section 2.2.1 Equation (1) applied to the unbound observations; b not reported
  })

  model({
    # ---- Individual parameters -----------------------------------
    # vc and cl are referenced to the UNBOUND concentration in this
    # parameterisation.
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc + e_crcl_vc * CRCL)
    cl <- exp(lcl + etalcl + e_age_cl * AGE)
    fu <- exp(lfu + etalfu)

    kel <- cl / vc

    # ---- ODE system (Section 2.2.1) -------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observations ---------------------------------------------
    # The unbound concentration is the primary prediction here; the
    # total plasma concentration is recovered through the same fixed
    # binding relationship Cu = Cc * fu (Section 2.2.1), inverted.
    Cu <- central / vc
    Cc <- Cu / fu

    Cc ~ add(addSd) + prop(propSd)
    Cu ~ add(addSd_Cu) + prop(propSd_Cu)
  })
}
