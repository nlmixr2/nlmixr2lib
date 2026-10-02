Frederiksen_2021_tedatioxetine <- function() {
  description <- paste(
    "Joint parent + metabolite population PK model for oral tedatioxetine",
    "(Lu AA24530) and its major CYP2D6-dependent metabolite Lu AA37208,",
    "pooled across six phase I and one phase II study in 578 healthy",
    "subjects and patients with major depressive disorder (Frederiksen",
    "2021). Tedatioxetine is dosed into an absorption-site compartment",
    "that drains by two competing first-order routes: intact drug into a",
    "two-compartment tedatioxetine disposition model (ka), and",
    "pre-systemic first-pass metabolism into an amount-only precursor pool",
    "(the intermediate Lu AA37209; k_precursor_luaa37208_form). The",
    "precursor converts to Lu AA37208 (k_luaa37208_form), which also",
    "receives systemic formation from the full tedatioxetine clearance",
    "(CL_CYP2D6) and then follows its own two-compartment disposition. All",
    "of the parent's systemic clearance is assumed CYP2D6-mediated",
    "formation of Lu AA37208. Covariates retained: food state on the",
    "pre-systemic formation rate constant, and age (linear, centred at 37",
    "years) on the metabolite clearance. Residual error is proportional",
    "and shared between the two analytes. Fit in NONMEM 7.4 (SAEM +",
    "importance sampling). The paper's purpose was to quantify in vivo",
    "CYP2D6 activity per genotype from the individual formation-clearance",
    "estimates; the CYP2D6 genotype itself is NOT a covariate in the",
    "population PK model."
  )
  reference <- paste(
    "Frederiksen T, Areberg J, Schmidt E, Stage TB, Brosen K. Cytochrome",
    "P450 2D6 genotype-phenotype characterization through population",
    "pharmacokinetic modeling of tedatioxetine. CPT Pharmacometrics Syst",
    "Pharmacol. 2021;10(9):983-993. doi:10.1002/psp4.12635.",
    "PMCID: PMC8452298.",
    sep = " "
  )
  vignette <- "Frederiksen_2021_tedatioxetine"
  units <- list(
    time = "h",
    dosing = "mg",
    # A(3)/V3 and A(5)/V5 with dose in mg and volumes in L give mg/L; the
    # supplementary control stream applies no unit conversion. The metabolite
    # is carried in tedatioxetine-mass equivalents (no molar-mass correction on
    # the parent-to-metabolite flux; see compartmentData note).
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against the Frederiksen 2021 supplementary
  # NONMEM control stream ($MODEL / $DES, PSP4-10-983-s004.docx) and Figure 2.
  compartmentData <- list(
    depot = list(analyte = "tedatioxetine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tedatioxetine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tedatioxetine", units = "mg", specimen = "tissue", verified = TRUE),
    precursor_luaa37208 = list(
      analyte = "Lu AA37209 (Lu AA37208 precursor)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_luaa37208 = list(analyte = "Lu AA37208", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_luaa37208 = list(analyte = "Lu AA37208", units = "mg", specimen = "tissue", verified = TRUE)
  )

  # The second $MODEL COMP=(DEPOT) is an amount-only kinetic intermediate with
  # no volume and no measured concentration (the control stream sets only
  # S3=V3 and S5=V5), which is why it is `precursor_luaa37208` rather than a
  # `central_luaa37208`-style state with an invented volume. It holds the
  # intermediate Lu AA37209 (Frederiksen 2021, Figure 1: tedatioxetine ->
  # Lu AA37209 -> Lu AA37208). `precursor_<metab>` and the `luaa37208` suffix
  # are registered in inst/references/compartment-names.md.

  covariateData <- list(
    AGE = list(
      description = "Age at study entry",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear covariate on the Lu AA37208 clearance, centred at the",
        "population median of 37 years (Table 1): the control stream forms",
        "CLMET = (THETA(9) + (AGE - 37) * THETA(14)) * exp(eta), i.e. an",
        "ADDITIVE deviation on the linear clearance scale rather than the",
        "usual multiplicative power model. Age on CLMET was the second",
        "covariate retained in forward inclusion (OFV -14; Results,",
        "'PopPK analysis'). Older subjects have a lower metabolite clearance",
        "(THETA(14) = -0.0830 L/h per year)."
      ),
      source_name = "AGE"
    ),
    FED = list(
      description = "Fed state at the time of dosing (1 = fed, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Food state selects between two separately estimated values of the",
        "pre-systemic formation rate constant (control-stream KGMET):",
        "11.1 /h fasted and 0.0286 /h fed. Food was the single most",
        "significant covariate in forward inclusion (OFV -596; Results,",
        "'PopPK analysis'). NOTE: Table 2 of the paper labels this",
        "food-affected rate constant 'ka,met' and the 0.0972 /h",
        "precursor-to-metabolite conversion rate 'kg,met', but the",
        "executable supplementary control stream applies the food effect to",
        "its KGMET (the depot-to-precursor branch) and leaves KAMET (the",
        "precursor-to-central conversion) food-independent -- the two labels",
        "are swapped between the printed table and the code. This model",
        "follows the executable control stream."
      ),
      source_name = "FED"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 578L,
    n_studies = 7L,
    age_range = "18-80 years",
    age_median = "37 years (IQR 28-50)",
    weight_range = "42-140 kg",
    weight_median = "72 kg (IQR 61-82)",
    sex_female_pct = 46.2,
    race_ethnicity = c(
      Caucasian = 442,
      `Black or African American` = 31,
      Asian = 88,
      Other = 17
    ),
    disease_state = "220 healthy subjects and 358 patients with major depressive disorder",
    dose_range = "2-60 mg oral tedatioxetine (single and multiple dose across the seven pooled studies)",
    regions = "Multiregional phase I and phase II studies (H. Lundbeck A/S)",
    notes = paste(
      "Six phase I studies (dense PK sampling) plus one phase II",
      "dose-finding study in MDD patients (sparse sampling, <=3 samples over",
      "6 weeks). 5373 quantifiable tedatioxetine and 5449 quantifiable",
      "Lu AA37208 plasma concentrations. CYP2D6 phenotype distribution",
      "(UM/NM/IM/PM/missing) was 11/291/200/32/44; genotype was used only in",
      "the downstream activity-score analysis, not as a population-PK",
      "covariate. Creatinine clearance (Cockcroft-Gault) was screened but",
      "not retained. Tedatioxetine (Lu AA24530) development for MDD was",
      "terminated; the data are reused here to study CYP2D6",
      "genotype-phenotype relationships."
    )
  )

  ini({
    # ---- Tedatioxetine absorption and competing pre-systemic first-pass ----
    # Control stream $DES eq 1: dA(1)/dt = -(KA + KGMET) * A(1). The
    # absorption site drains by two parallel first-order routes: ka to
    # tedatioxetine central, and k_precursor_luaa37208_form (control-stream
    # KGMET, food-affected) to the Lu AA37209 precursor pool.
    lka <- log(0.195)
    label("Tedatioxetine absorption rate constant (ka, 1/h)") # Table 2, row 'Absorption rate constant, tedatioxetine (ka)': 0.195 1/h (6.4 %RSE)

    # Pre-systemic formation rate constant, control-stream KGMET (THETA(3)
    # fasted / THETA(13) fed). Table 2 mislabels this as 'ka,met'; the
    # executable control stream applies it to the depot->precursor branch.
    # The fasted value is the log-scale base; the food effect is encoded in
    # the canonical linear-deviation form (as in the sibling
    # Frederiksen_2023_brexpiprazole model): the fed value is
    # base * (1 + e_fed) with e_fed = 0.0286/11.1 - 1. This is algebraically
    # identical to the control stream's KGMET = (THETA3*(1-FED)+THETA13*FED)
    # * exp(eta).
    lk_precursor_luaa37208_form <- log(11.1)
    label("Pre-systemic Lu AA37208-precursor formation rate constant, fasted (1/h)") # Table 2, row 'Absorption rate constant, Lu AA37208 (ka,met) fasted' [control-stream KGMET fasting, THETA(3)]: 11.1 1/h (35.6 %RSE)
    e_fed_k_precursor_luaa37208_form <- -0.99742
    label("Food effect on the pre-systemic formation rate constant (fractional deviation)") # Table 2: fed value [control-stream KGMET fed, THETA(13)] 0.0286 1/h; e_fed = 0.0286/11.1 - 1 = -0.99742

    # ---- Precursor to Lu AA37208 conversion ----
    # Control stream $DES eq 2: A(2) drains by KAMET into central_luaa37208.
    # Table 2 mislabels this rate 'kg,met'.
    lk_luaa37208_form <- log(0.0972)
    label("Lu AA37209-precursor to Lu AA37208 conversion rate constant (1/h)") # Table 2, row 'Rate constant formation of Lu AA37208 (kg,met)' [control-stream KAMET, THETA(2)]: 0.0972 1/h (8.2 %RSE)

    ltlag <- log(0.652)
    label("Tedatioxetine absorption lag time (ALAG1, h)") # Table 2, row 'Lag-time (ALAG)': 0.652 h (0.7 %RSE)

    # ---- Tedatioxetine disposition (apparent, /F) ----
    lcl <- log(30.5)
    label("Tedatioxetine clearance = CYP2D6-mediated formation clearance (CL_CYP2D6, L/h)") # Table 2, row 'Clearance, tedatioxetine (CL_CYP2D6)': 30.5 L/h (6.6 %RSE)
    lvc <- log(1380)
    label("Tedatioxetine central volume (V3, L)") # Table 2, row 'Volume of distribution, central compartment, tedatioxetine (V3)': 1380 L (4.3 %RSE)
    lq <- log(39.1)
    label("Tedatioxetine inter-compartmental clearance (Q, L/h)") # Table 2, row 'Intercompartmental clearance, tedatioxetine (Q)': 39.1 L/h (0.8 %RSE)
    lvp <- log(507)
    label("Tedatioxetine peripheral volume (V4, L)") # Table 2, row 'Volume of distribution, peripheral compartment, tedatioxetine (V4)': 507 L (0.8 %RSE)

    # ---- Lu AA37208 disposition (apparent) and age effect on clearance ----
    # Control stream: CLMET = (THETA(9) + (AGE - 37) * THETA(14)) * exp(eta).
    # THETA(9) is the linear typical clearance; the age effect is an additive
    # deviation on that linear scale, so both are kept off the log scale.
    tvcl_luaa37208 <- 11.9
    label("Lu AA37208 clearance at age 37 (CLmet, L/h)") # Table 2, row 'Clearance, Lu AA37208 (CLmet)' [control-stream THETA(9)]: 11.9 L/h (3.4 %RSE)
    e_age_cl_luaa37208 <- -0.0830
    label("Age effect on Lu AA37208 clearance (linear, L/h per year)") # Table 2, row 'Age on CLmet' [control-stream THETA(14)]: -0.0830 (20.0 %RSE)
    lvc_luaa37208 <- log(33.1)
    label("Lu AA37208 central volume (V5, L)") # Table 2, row 'Volume of distribution, central compartment, Lu AA37208 (V5)': 33.1 L (5.9 %RSE)
    lq_luaa37208 <- log(0.940)
    label("Lu AA37208 inter-compartmental clearance (Qmet, L/h)") # Table 2, row 'Intercompartmental clearance, Lu AA37208 (Qmet)': 0.940 L/h (0.7 %RSE)
    lvp_luaa37208 <- log(12.2)
    label("Lu AA37208 peripheral volume (V6, L)") # Table 2, row 'Volume of distribution, peripheral compartment, Lu AA37208 (V6)': 12.2 L (8.1 %RSE)

    # ---- Between-subject variability ----
    # Table 2 reports IIV as %CV. The log-scale variance is omega^2 =
    # (CV/100)^2: this reading is forced by the CLmet-V5 covariance (0.372),
    # which yields a correlation of 0.986 (valid) under (CV/100)^2 but a
    # correlation > 1 (impossible) under the exact log-normal
    # omega^2 = log((CV/100)^2 + 1). See vignette Source trace.
    etalka ~ 0.4871 # Table 2, ka IIV 69.79% CV -> (0.6979)^2 = 0.4871
    etalk_precursor_luaa37208_form ~ 5.130 # Table 2, 'ka,met' IIV 226.50% CV -> (2.2650)^2 = 5.130 (control-stream KGMET eta)
    etalk_luaa37208_form ~ 0.8540 # Table 2, 'kg,met' IIV 92.41% CV -> (0.9241)^2 = 0.8540 (control-stream KAMET eta)
    # $OMEGA BLOCK(2) over (V3, CL): cov(CL,V3) = 0.079 (Table 2).
    etalvc + etalcl ~ c(
      0.18302, # V3 IIV 42.78% CV -> (0.4278)^2
      0.079, 0.69706 # cov(CL,V3) = 0.079; CL IIV 83.49% CV -> (0.8349)^2
    )
    # $OMEGA BLOCK(2) over (V5, CLmet): cov(CLmet,V5) = 0.372 (Table 2).
    etalvc_luaa37208 + etalcl_luaa37208 ~ c(
      0.47005, # V5 IIV 68.56% CV -> (0.6856)^2
      0.372, 0.30305 # cov(CLmet,V5) = 0.372; CLmet IIV 55.05% CV -> (0.5505)^2
    )

    # ---- Residual unexplained variability ----
    # $ERROR uses log-transform-both-sides additive error (Y = LOG(F+0.001) +
    # ERR(1)); Table 2 reports it as a proportional error of 23.6%. A SINGLE
    # sigma was estimated and applied to both analytes (Results, 'PopPK
    # analysis'); nlmixr2 requires one residual parameter per endpoint, so the
    # shared estimate is written into both propSd and propSd_luaa37208.
    propSd <- 0.236
    label("Tedatioxetine proportional residual SD (unitless)") # Table 2, row 'Residual error (proportional)': 23.6% (0.3 %RSE); single shared sigma
    propSd_luaa37208 <- 0.236
    label("Lu AA37208 proportional residual SD (unitless; same estimate as tedatioxetine)") # Table 2, row 'Residual error (proportional)': 23.6% (0.3 %RSE); single shared sigma
  })

  model({
    # 1. Individual parameters (mu-referenced: each log-parameter on its own
    #    simple line, then combined).
    ka <- exp(lka + etalka)
    # Food effect on the pre-systemic formation rate (control-stream KGMET),
    # in the canonical linear-deviation form: fasted base * (1 + e_fed * FED).
    k_precursor_luaa37208_form <- exp(lk_precursor_luaa37208_form + etalk_precursor_luaa37208_form) *
      (1 + e_fed_k_precursor_luaa37208_form * FED)
    k_luaa37208_form <- exp(lk_luaa37208_form + etalk_luaa37208_form)
    # The control stream declares etas on Q, V4, Qmet, V6 and ALAG1 but fixes
    # their variance to 0 ('0 FIX'), so they carry no IIV and are written here
    # without an eta.
    tlag <- exp(ltlag)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    # Linear age effect (additive) on the linear metabolite clearance scale,
    # then log-normal IIV (control stream: CLMET = (THETA9 + (AGE-37)*THETA14)
    # * exp(eta)).
    lcl_luaa37208 <- log(tvcl_luaa37208 + (AGE - 37) * e_age_cl_luaa37208)
    cl_luaa37208 <- exp(lcl_luaa37208 + etalcl_luaa37208)
    vc_luaa37208 <- exp(lvc_luaa37208 + etalvc_luaa37208)
    q_luaa37208 <- exp(lq_luaa37208)
    vp_luaa37208 <- exp(lvp_luaa37208)

    # 2. Concentrations entering the clearance-parameterised ODEs.
    Cc <- central / vc
    Cp <- peripheral1 / vp
    Cc_luaa37208 <- central_luaa37208 / vc_luaa37208
    Cp_luaa37208 <- peripheral1_luaa37208 / vp_luaa37208

    # 3. ODE system, supplementary control stream $DES verbatim.
    # A(1): tedatioxetine absorption site drains by both routes.
    d/dt(depot) <- -ka * depot - k_precursor_luaa37208_form * depot
    # A(2): Lu AA37209 precursor, fed pre-systemically from the absorption
    # site and drained by conversion to Lu AA37208.
    d/dt(precursor_luaa37208) <- k_precursor_luaa37208_form * depot -
      k_luaa37208_form * precursor_luaa37208
    # A(3): tedatioxetine central. All systemic clearance forms Lu AA37208.
    d/dt(central) <- ka * depot - q * Cc + q * Cp - cl * Cc
    # A(4): tedatioxetine peripheral.
    d/dt(peripheral1) <- q * Cc - q * Cp
    # A(5): Lu AA37208 central, fed by precursor conversion and by systemic
    # CYP2D6-mediated formation (cl * Cc).
    d/dt(central_luaa37208) <- k_luaa37208_form * precursor_luaa37208 +
      cl * Cc - q_luaa37208 * Cc_luaa37208 + q_luaa37208 * Cp_luaa37208 -
      cl_luaa37208 * Cc_luaa37208
    # A(6): Lu AA37208 peripheral.
    d/dt(peripheral1_luaa37208) <- q_luaa37208 * Cc_luaa37208 -
      q_luaa37208 * Cp_luaa37208

    alag(depot) <- tlag

    # 4. Observations. Proportional residual shared by both analytes.
    Cc ~ prop(propSd)
    Cc_luaa37208 ~ prop(propSd_luaa37208)
  })
}
