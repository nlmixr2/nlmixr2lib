Hayashi_1998_epoetinBeta <- function() {
  description <- "One-compartment population PK model for subcutaneous recombinant human erythropoietin (epoetin beta) in healthy adult male Japanese volunteers with a constant endogenous EPO production rate carrying a fixed circadian sinusoid (acrophase near midnight) feeding the central compartment, and body weight as a power covariate on apparent absorption rate ka and apparent central volume V/F, plus serum creatinine and age as power covariates on the elimination rate constant k_e (reparameterised here onto canonical CL/F so the k_e covariates ride on CL/F together with the V/F weight exponent); apparent V/F and E/F throughout because bioavailability was not separately estimable from this SC-only study (Hayashi 1998)."
  reference <- paste(
    "Hayashi N, Kinoshita H, Yukawa E, Higuchi S.",
    "Pharmacokinetic analysis of subcutaneous erythropoietin administration",
    "with nonlinear mixed effect model including endogenous production.",
    "Br J Clin Pharmacol. 1998;46(1):11-19.",
    "doi:10.1046/j.1365-2125.1998.00043.x",
    sep = " "
  )
  vignette <- "Hayashi_1998_epoetinBeta"
  units <- list(
    time          = "h",
    dosing        = "IU",
    concentration = "IU/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at study entry",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on apparent ka (exponent -1.92) and on apparent CL/F and V/F (shared exponent 0.776 inherited from the paper's V/F WT effect). Reference WT = 62.0 kg (Hayashi 1998 Table 1 cohort mean). Healthy-adult-male range 51.0-79.0 kg.",
      source_name        = "Body weight"
    ),
    AGE = list(
      description        = "Age at study entry",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on apparent k_e (paper's parameterisation) with exponent -1.13; under the canonical CL/F reparameterisation here the same exponent rides on CL/F. Reference AGE = 22.7 y (Hayashi 1998 Table 1 cohort mean). Healthy-adult-male range 20-29 y; the strong negative exponent within this narrow window led the authors to discuss bone-marrow contribution to EPO clearance.",
      source_name        = "Age"
    ),
    CREAT = list(
      description        = "Serum creatinine at study entry",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on apparent k_e (paper's parameterisation) with exponent -0.542; under the canonical CL/F reparameterisation here the same exponent rides on CL/F. Reference CREAT = 0.98 mg/dL (Hayashi 1998 Table 1 cohort mean). Healthy-volunteer range 0.8-1.2 mg/dL.",
      source_name        = "Creatinine"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 48L,
    n_studies       = 1L,
    age_range       = "20-29 years",
    age_median      = "mean 22.7 (s.d. 2.1) years",
    weight_range    = "51.0-79.0 kg",
    weight_median   = "mean 62.0 (s.d. 5.7) kg",
    height_range    = "157.0-185.7 cm",
    height_median   = "mean 170.5 (s.d. 5.6) cm",
    sex_female_pct  = 0,
    race_ethnicity  = c(Asian = 100),
    creat_range     = "0.8-1.2 mg/dL (mean 0.98)",
    disease_state   = "Healthy adult male Japanese volunteers; normal complete blood count, platelet count, prothrombin time, and activated partial thromboplastin generation time; no clinically significant cardiovascular, renal, hepatic, metabolic, gastrointestinal, neurologic or endocrine disorders.",
    dose_range      = "1500 IU or 3000 IU subcutaneous (forearm) at 09:00 h, two administrations 2 weeks apart per subject. 16 subjects @ 1500 IU and 32 subjects @ 3000 IU.",
    regions         = "Japan",
    notes           = "Single-centre Phase I bioequivalence-style study at Kannondai Clinic (Ibaraki, Japan) comparing two Epoetin beta formulations; pooled for this popPK analysis. 1056 plasma EPO concentrations measured by radioimmunoassay (anti-recombinant-EPO antibody; LOD 2.1 IU/L; within-run CV 4.5-5.7%; between-run CV 1.9-6.7%). Sampling at -1, 3, 6, 9, 12, 15, 24, 36, 48, 72, 96 h after each administration. Pre-dose baseline EPO (Table 2) was 23.65 +/- 0.65 IU/L (3000 IU group) and 25.14 +/- 0.58 IU/L (1500 IU group). Demographics in Table 1; final pop-PK estimates in Table 4."
  )

  ini({
    # Structural parameters. The paper estimates the apparent absorption rate
    # k_a, the apparent elimination rate k_e, and the apparent central volume
    # V/F (Equation 2 and Table 4). To match the nlmixr2lib canonical CL/V
    # parameterisation, we reparameterise as (k_a, CL/F, V/F) using the exact
    # algebraic identity CL/F = k_e * V/F:
    #   CL/F_typ = 0.207 h^-1 * 14.4 L = 2.978 L/h at reference covariates.
    # Covariate effects on CL/F are the product of paper's effects on k_e and
    # V/F: CREAT and AGE on k_e carry to CL/F unchanged; the V/F WT exponent
    # 0.776 propagates from V/F to CL/F (the paper assigned WT only to V/F).
    # Reference covariates are the population means (Table 1):
    #   WT_ref = 62 kg, AGE_ref = 22.7 y, CREAT_ref = 0.98 mg/dL.
    lka <- log(0.0430); label("Apparent absorption rate constant ka (1/h) at WT=62 kg")  # Table 4 Results text (mean 0.0430 +/- 0.002); also reproduces from k_a = 118 * 62^-1.92.
    lcl <- log(2.978);  label("Apparent clearance CL/F (L/h) at reference covariates")    # Derived: k_e * V/F = 0.207 * 14.4 (Table 4)
    lvc <- log(14.4);   label("Apparent central volume V/F (L) at WT=62 kg")              # Table 4 (V/F = 0.585 * 62^0.776)

    # Covariate effects (power form: param * (cov/cov_ref)^exp).
    # All estimated by the paper with reported SE (Table 4); not held fixed.
    e_wt_ka    <- -1.92;  label("Power exponent of (WT/62) on apparent ka (unitless)")                                # Table 4 (theta_ka^WT; SE 0.95)
    e_wt_cl_vc <-  0.776; label("Shared power exponent of (WT/62) on apparent CL/F and V/F (unitless)")                # Table 4 (theta_V/F^WT; SE 0.235; propagated to CL/F via CL = ke * V/F)
    e_creat_cl <- -0.542; label("Power exponent of (CREAT/0.98) on apparent CL/F (unitless; from k_e covariate model)") # Table 4 (theta_ke^Cr; SE 0.288)
    e_age_cl   <- -1.13;  label("Power exponent of (AGE/22.7) on apparent CL/F (unitless; from k_e covariate model)")   # Table 4 (theta_ke^Age; SE 0.69)

    # Apparent endogenous EPO production rate (mesor) with fixed circadian
    # sinusoid. The paper estimates three period-specific E/F values:
    #   E1/F = 76.1 IU/h (around 1st administration; all subjects)
    #   E2/F = 75.5 IU/h (around 2nd administration; 1500 IU dose group)
    #   E3/F = 91.6 IU/h (around 2nd administration; 3000 IU dose group)
    # We use E1/F as the standing mesor here; the 3000 IU 2nd-period
    # elevation is a discussion point in the paper (interpreted as
    # increased endogenous production induced by the prior dose), is not
    # a feature of the structural single-dose model, and is documented
    # in the validation vignette.
    lrbase <- log(76.1);   label("Apparent endogenous EPO production rate E1/F (IU/h), 24-h mean")  # Table 4 (E1/F = 76.1; SE 5.1)
    lra    <- log(0.0986); label("Fractional amplitude of circadian endogenous production (unitless)")  # Table 4 (Amp = 9.86%; SE 0.0200)
    ltacro <- log(0.256);  label("Acrophase: clock-hour at which endogenous production peaks (h after midnight)")  # Table 4 (Phase t_0 = 3.86 rad; SE 0.18); acrophase = 15 - 24/(2*pi)*t_0 = 0.256 h (~ 00:15; "around 24.00 h" per Discussion).

    # Inter-individual variability. The paper's Methods (Equation 3) uses a
    # linear-scale "(1 + gamma)" eta on each parameter: param_ij = param_typ
    # * (1 + gamma_ij) with independent gamma_ij ~ N(0, omega^2). The
    # log-normal eta variance in nlmixr2 maps via omega^2 = log(CV^2 + 1).
    # The (k_a, k_e, V/F, E/F) etas are independent in the paper. Under our
    # canonical reparameterisation log(CL) = log(k_e) + log(V/F):
    #   Var(eta_lcl) = sigma^2_lke + sigma^2_lvc;
    #   Cov(eta_lcl, eta_lvc) = sigma^2_lvc; eta_lka and eta_lrbase remain independent.
    etalka          ~ 0.0889                              # log(0.305^2 + 1) = 0.0889; Table 4 IIV(k_a) = 30.5% (SE 3.7)
    etalcl + etalvc ~ c(0.01604, 0.00335, 0.00335)        # block: var_lcl = log(0.113^2 + 1) + log(0.0579^2 + 1) = 0.01269 + 0.00335 = 0.01604; cov = log(0.0579^2 + 1) = 0.00335; var_lvc = log(0.0579^2 + 1) = 0.00335. Table 4 IIV(k_e) = 11.3% (SE 4.5) and IIV(V/F) = 5.79% (SE 6.60).
    etalrbase       ~ 0.00710                             # log(0.0844^2 + 1) = 0.00710; Table 4 IIV(END/F) = 8.44% (SE 4.45)

    # Residual error (proportional, on linear-space concentration). Paper's
    # Methods Equation 3 form: C_ij = C_pred,ij * (1 + epsilon_ij) with
    # epsilon_ij ~ N(0, sigma^2). Reported as 13.9% in Table 4.
    propSd <- 0.139; label("Proportional residual error (fraction)")  # Table 4 (sigma = 13.9%; SE 0.7)
  })

  model({
    # Reference covariate values (Hayashi 1998 Table 1 cohort means).
    wt_ref    <- 62.0
    age_ref   <- 22.7
    creat_ref <- 0.98

    # Individual PK parameters. Reparameterised CL/F, V/F per the
    # algebraic identity CL/F = k_e * V/F (see ini-block comment).
    ka  <- exp(lka + etalka) * (WT / wt_ref)^e_wt_ka
    cl  <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl_vc *
           (CREAT / creat_ref)^e_creat_cl * (AGE / age_ref)^e_age_cl
    vc  <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_cl_vc
    kel <- cl / vc

    # Endogenous EPO production (IU/h) with circadian sinusoidal modulation.
    # END(clock_t) = rbase * (1 + ra * sin(2*pi/24 * (tacro + 6 - clock_t)))
    # The form is equivalent to Hayashi 1998 Equation 1:
    #   END_ij = E_j * (1 + Amp_j * sin(2*pi/24 * (6 - T_i + T_peak_j)))
    # The paper dosed every subject at 09:00 h (Methods), so the model's
    # rxode2 time t (h since dose) maps to clock time as clock_t = t + 9.
    # At clock_t = tacro the sinusoid argument is pi/2 and END = rbase * (1 + ra).
    rbase <- exp(lrbase + etalrbase)
    ra    <- exp(lra)
    tacro <- exp(ltacro)

    clock_t <- t + 9
    end_t   <- rbase * (1 + ra * sin(2 * pi / 24 * (tacro + 6 - clock_t)))

    # One-compartment open model with first-order absorption plus endogenous
    # input. The paper reports flip-flop kinetics: typical k_a (~0.043 h^-1)
    # is smaller than k_e (~0.207 h^-1), so the apparent terminal half-life
    # is governed by absorption (Methods, "flip-flop phenomenon"; Results).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central + end_t

    # Initial condition at t = 0 (clock 09:00): pre-dose endogenous steady
    # state. central(0) = end_t(0) / kel reproduces the observed pre-dose
    # baseline (Table 2: 23.65 +/- 0.65 IU/L for 3000 IU group). With typical
    # parameters the sinusoid at clock 09:00 evaluates to (tacro + 6 - 9) =
    # tacro - 3 = -2.744 h, giving end_t(0) = 76.1 * (1 + 0.0986 * sin(2*pi
    # * -2.744 / 24)) = 76.1 * (1 - 0.065) = 71.2 IU/h, central(0) = 71.2 /
    # 0.207 = 343.8 IU, and Cc(0) = 343.8 / 14.4 = 23.9 IU/L.
    central(0) <- rbase * (1 + ra * sin(2 * pi / 24 * (tacro - 3))) / kel

    # Observation: apparent concentration = central amount / apparent V/F.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
