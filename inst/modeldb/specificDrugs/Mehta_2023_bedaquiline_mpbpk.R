Mehta_2023_bedaquiline_mpbpk <- function() {
  description <- paste(
    "mPBPK (minimal physiologically based, translational mouse-to-human).",
    "Bedaquiline and its N-desmethyl metabolite M2 in pulmonary tuberculosis patients, with",
    "explicit cavitary lung-lesion and uninvolved-lung site-of-action compartments.",
    "Parallel first-order absorption from a fast and a slow depot, a well-stirred liver",
    "compartment that converts bedaquiline to M2 via a hepatic extraction ratio, a lumped",
    "peripheral tissue pool for parent and for M2, and perfusion-limited lesion / uninvolved-lung",
    "compartments driven by penetration ratios. Physiological volumes and flows are computed from",
    "body weight (Brown 1997); clearance terms are allometrically scaled from the mouse fit with an",
    "exponent of 0.75. The lesion and uninvolved-lung states hold concentrations (mg/L), not",
    "amounts. Intended for typical-value and population simulation of site-of-action target",
    "attainment; the paper did not fit this model to individual human data.",
    sep = " "
  )
  reference <- paste(
    "Mehta K, Guo T, van der Graaf PH, van Hasselt JGC.",
    "Predictions of Bedaquiline and Pretomanid Target Attainment in Lung Lesions of",
    "Tuberculosis Patients using Translational Minimal Physiologically Based",
    "Pharmacokinetic Modeling. Clin Pharmacokinet. 2023;62(3):519-532.",
    "doi:10.1007/s40262-023-01217-7.",
    sep = " "
  )
  vignette <- "Mehta_2023_tb_lesion_mpbpk"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    weight = "kg"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot1         = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = FALSE),
    depot2         = list(analyte = "bedaquiline", units = "mg", specimen = "administration site", verified = FALSE),
    liver          = list(analyte = "bedaquiline", units = "mg", specimen = "tissue", verified = FALSE),
    blood          = list(analyte = "bedaquiline", units = "mg", specimen = "blood cell", verified = FALSE),
    peripheral1    = list(analyte = "bedaquiline", units = "mg", specimen = "plasma", verified = FALSE),
    blood_m2       = list(analyte = "M2", units = "mg", specimen = "blood cell", verified = FALSE),
    peripheral1_m2 = list(analyte = "M2", units = "mg", specimen = "plasma", verified = FALSE),
    lesion         = list(analyte = "bedaquiline", units = "mg", specimen = "tissue", verified = FALSE),
    lung           = list(analyte = "bedaquiline", units = "mg", specimen = "tissue", verified = FALSE),
    lesion_m2      = list(analyte = "M2", units = "mg", specimen = "tissue", verified = FALSE),
    lung_m2        = list(analyte = "M2", units = "mg", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives every physiological quantity in model(): cardiac output Qc = 312 * (WT/70)^0.75",
        "L/h, hepatic blood flow Qh = 0.227 * Qc, liver volume 0.0257 * WT L, lung volume",
        "0.0076 * WT L, blood-reservoir volume 0.0771 * WT L (0.0514 venous + 0.0257 arterial),",
        "and the lumped tissue volume as the residual. Also carries the 0.75 allometric exponent",
        "on CLint and CLM2 relative to a 70 kg reference. Body weights for the 500 virtual",
        "patients were sampled from the TB-PACTS clinical-trial weight distribution",
        "(Methods 2.2, ESM S1); the distribution itself is not tabulated in the paper.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 500L,
    disease_state  = paste(
      "Pulmonary (drug-susceptible and multidrug-resistant) tuberculosis with cavitary disease.",
      "Virtual population; cavity presence/absence and cavity size were sampled from the observed",
      "TB-PACTS distributions (Methods 2.2).",
      sep = " "
    ),
    dose_range     = paste(
      "Standard bedaquiline 400 mg once daily for 14 days followed by 200 mg three times weekly;",
      "alternative 200 mg once daily for 8 weeks followed by 100 mg once daily (Methods 2.5).",
      "Clinical validation data spanned 100-700 mg (ESM S1).",
      sep = " "
    ),
    notes          = paste(
      "The structural model and every drug-specific parameter were estimated on mouse data",
      "(plasma, liver, lesion and uninvolved lung after a single oral 25 mg/kg dose; ESM S1),",
      "then translated to humans by swapping in human physiology and allometrically scaling",
      "clearance. CLint and CLM2 were additionally calibrated to median rich concentration-time",
      "profiles from TB patients in one dose group (400 mg day 1, 300 mg day 2, 200 mg days 3-14)",
      "because direct allometric scaling under-predicted bedaquiline and over-predicted M2",
      "(Results 3.2). 50 trials x 500 virtual patients were simulated to propagate parameter",
      "uncertainty (RSE). No human lesion or uninvolved-lung PK data exist for bedaquiline.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Absorption -- two parallel first-order depots. The full dose is placed
    # in BOTH depot1 and depot2; the fdepot / (1 - fdepot) split is applied
    # at the transfer into the liver, exactly as in ESM S2.
    # ---------------------------------------------------------------------
    lka1    <- log(1.3)      ; label("Absorption rate from the fast depot, depot1 (1/h)")            # Table 1 bedaquiline ka1 = 1.3 (RSE 37.8%), same value mice and humans
    lka2    <- log(0.00501)  ; label("Absorption rate from the slow depot, depot2 (1/h)")            # Table 1 bedaquiline ka2 = 0.00501 (RSE 9.13%), same value mice and humans
    lfdepot <- log(0.609)    ; label("Fraction of dose absorbed through the fast depot, depot1 (unitless)") # Table 1 bedaquiline frc = 0.609 (RSE 11.9%)

    # ---------------------------------------------------------------------
    # Disposition. CLint is the intrinsic clearance entering the well-stirred
    # liver extraction ratio; CLM2 is the linear clearance of M2 from blood.
    # Both are the HUMAN column of Table 1 (mouse values 1.21 and 0.0119 L/h
    # were allometrically scaled and then calibrated to patient data).
    # ---------------------------------------------------------------------
    lclint  <- log(60.3)     ; label("Intrinsic hepatic clearance of bedaquiline, CLint (L/h, 70 kg)") # Table 1 bedaquiline CLint human = 60.3 (RSE 16.7%)
    lcl_m2  <- log(45.9)     ; label("Clearance of the M2 metabolite from blood (L/h, 70 kg)")       # Table 1 bedaquiline CLM2 human = 45.9 (RSE 13.6%)

    # ---------------------------------------------------------------------
    # Tissue:blood equilibrium concentration ratios. KpT / KpTM2 are the
    # partition coefficients of the lumped peripheral pool; R_le / R_ul are
    # the lesion and uninvolved-lung penetration ratios. All are unitless and
    # were assumed identical in mice and humans (Methods 2.2).
    # ---------------------------------------------------------------------
    lk_peripheral1    <- log(4.45)  ; label("Bedaquiline lumped-tissue partition coefficient, KpT (unitless)")    # Table 1 bedaquiline KpT = 4.45 (RSE 15.3%)
    lk_peripheral1_m2 <- log(9.54)  ; label("M2 lumped-tissue partition coefficient, KpTM2 (unitless)")           # Table 1 bedaquiline KpTM2 = 9.54 (RSE 18.4%). ESM S2 prints 'KpTM2=18.4', which is that row's %RSE mis-copied; see vignette Errata.
    lk_lesion         <- log(11)    ; label("Bedaquiline lung-lesion penetration ratio, R_le (unitless)")         # Table 1 bedaquiline R_le = 11 (RSE 10.7%); ESM S2 Rles = 11
    lk_lung           <- log(10.2)  ; label("Bedaquiline uninvolved-lung penetration ratio, R_ul (unitless)")     # Table 1 bedaquiline R_ul = 10.2 (RSE 10.9%); ESM S2 RUL = 10.2
    lk_lesion_m2      <- log(88.4)  ; label("M2 lung-lesion penetration ratio, R_le,M2 (unitless)")               # Table 1 bedaquiline R_le M2 = 88.4 (RSE 5.72%); ESM S2 RlesM2 = 88.4
    lk_lung_m2        <- log(88.8)  ; label("M2 uninvolved-lung penetration ratio, R_ul,M2 (unitless)")           # Table 1 bedaquiline R_ul M2 = 88.8 (RSE 5.53%); ESM S2 RULM2 = 88.8

    # ---------------------------------------------------------------------
    # Fixed physicochemical / scaling constants.
    # ---------------------------------------------------------------------
    fu       <- fixed(0.1)      ; label("Fraction unbound in plasma, fup (unitless)")                # Table 1 bedaquiline fup = 0.1, literature value from van Heeswijk 2014 (reference 27); ESM S2 fup = 0.1
    lbp      <- fixed(log(1))   ; label("Blood-to-plasma concentration ratio, B:P (unitless)")       # Table 1 bedaquiline BP = 1 (reference 34)
    e_wt_cl  <- fixed(0.75)     ; label("Allometric exponent on CLint and CLM2 (unitless)")          # Methods 2.2: 'a previously known allometric exponent of 0.75 for CL' (reference 20); ESM S2 (BW/70)^0.75
    lvlef    <- fixed(log(0.0216)); label("Cavitary lesion volume as a fraction of total lung volume, VF_le (unitless)") # Methods 2.1.1: 'Mean volume of lesions (VF_le) was assumed to be 0.0216, calculated based on the mean total lesion volume, approximately 14 mL' (reference 19)

    # ---------------------------------------------------------------------
    # Between-subject variability. Table 1 states 40% log-normal IIV for the
    # human simulations. ESM S2 applies exp(eta) to ka1, ka2, CLint and CLM2
    # only -- KpT and KpTM2 are plain constants in the final code -- so the
    # etas below follow the final model code (operator decision, 2026-07-27).
    # Encoded as variance = 0.40^2 on the log scale.
    # ---------------------------------------------------------------------
    etalka1   ~ 0.16  # Table 1 bedaquiline 'IIV for ka1, ka2, CLint, CLM2, KpT, and KpTM2 (%) = 40'; ESM S2 ka1 = TVka1*exp(eta.ka1)
    etalka2   ~ 0.16  # Table 1 bedaquiline IIV row = 40%; ESM S2 ka2 = TVka2*exp(eta.ka2)
    etalclint ~ 0.16  # Table 1 bedaquiline IIV row = 40%; ESM S2 CLint = TVCLint*((BW/70)^0.75)*exp(eta.CLint)
    etalcl_m2 ~ 0.16  # Table 1 bedaquiline IIV row = 40%; ESM S2 CLM2 = TVCLM2*((BW/70)^0.75)*exp(eta.CLM2)

    # ---------------------------------------------------------------------
    # Residual error. Estimated on the MOUSE fit; the paper states 'Residual
    # errors were not included in the human simulations' (Table 1 footnote).
    # Carried here so the model is fittable; zero it for typical-value runs.
    # ---------------------------------------------------------------------
    propSd <- 0.43 ; label("Proportional residual error (fraction)")  # Table 1 footnote: 'combined bedaquiline plasma, liver, and M2 plasma = proportional 43%'
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual drug-specific parameters
    # -------------------------------------------------------------------
    ka1    <- exp(lka1 + etalka1)
    ka2    <- exp(lka2 + etalka2)
    fdepot <- exp(lfdepot)
    clint  <- exp(lclint + etalclint) * (WT / 70)^e_wt_cl
    cl_m2  <- exp(lcl_m2 + etalcl_m2) * (WT / 70)^e_wt_cl
    kpt    <- exp(lk_peripheral1)
    kpt_m2 <- exp(lk_peripheral1_m2)
    r_le   <- exp(lk_lesion)
    r_ul   <- exp(lk_lung)
    r_le_m2 <- exp(lk_lesion_m2)
    r_ul_m2 <- exp(lk_lung_m2)
    bp     <- exp(lbp)
    vlef   <- exp(lvlef)

    # -------------------------------------------------------------------
    # 2. Human physiology from body weight (ESM S2; Brown 1997 = ref 15).
    #    Flows in L/h, volumes in L.
    # -------------------------------------------------------------------
    qc   <- 312 * (WT / 70)^0.75   # cardiac output; Brown 1997 Table 22 (5200 mL/min at 70 kg)
    qh   <- qc * 0.227             # hepatic blood flow as a fraction of Qc; Brown 1997 Table 23
    qt   <- qc - qh                # flow to the lumped peripheral tissue
    vliv <- 0.0257 * WT            # liver volume; Brown 1997 Table 7 (human)
    vlu  <- 0.0076 * WT            # total lung volume; Brown 1997 Table 7
    vbl  <- (0.0514 + 0.0257) * WT # blood reservoir = venous + arterial
    vt   <- WT - vlu - vbl - vliv  # residual (lumped tissue) volume

    # Lesion and uninvolved-lung volumes, and their first-order exchange
    # rate constants k_i = Qc / V_i (main-text equation 2).
    vle <- vlef * vlu
    vul <- vlu - vle
    kle <- qc / vle
    kul <- qc / vul

    # Well-stirred liver extraction ratio.
    eh <- (fu * clint) / (qh + fu * clint)

    # Blood concentrations driving the site-of-action compartments.
    cbld    <- blood / vbl
    cbld_m2 <- blood_m2 / vbl

    # -------------------------------------------------------------------
    # 3. ODE system (ESM S2 'Final Model Codes'). Amounts in mg for
    #    depot1/depot2/liver/blood/peripheral1 and their M2 counterparts;
    #    the lesion and lung states hold CONCENTRATIONS in mg/L.
    # -------------------------------------------------------------------
    d/dt(depot1) <- -ka1 * depot1
    d/dt(depot2) <- -ka2 * depot2

    d/dt(liver) <- fdepot * ka1 * depot1 + (1 - fdepot) * ka2 * depot2 -
      qh * eh * liver / vliv -        # conversion to M2
      qh * (1 - eh) * liver / vliv +  # outflow to blood
      qh * blood / vbl                # inflow from blood

    d/dt(blood) <- qh * (1 - eh) * liver / vliv -
      qh * blood / vbl +
      peripheral1 * qt / (vt * kpt) - blood * qt / vbl

    d/dt(peripheral1) <- blood * qt / vbl - peripheral1 * qt / (vt * kpt)

    d/dt(blood_m2) <- qh * eh * liver / vliv - cl_m2 * blood_m2 / vbl +
      peripheral1_m2 * qt / (vt * kpt_m2) - blood_m2 * qt / vbl

    d/dt(peripheral1_m2) <- blood_m2 * qt / vbl -
      peripheral1_m2 * qt / (vt * kpt_m2)

    # Site of action -- main-text equation 1: d/dt(Ci) = ki * (Cbld * Ri - Ci)
    d/dt(lesion)    <- kle * (cbld * r_le - lesion)
    d/dt(lung)      <- kul * (cbld * r_ul - lung)
    d/dt(lesion_m2) <- kle * (cbld_m2 * r_le_m2 - lesion_m2)
    d/dt(lung_m2)   <- kul * (cbld_m2 * r_ul_m2 - lung_m2)

    # -------------------------------------------------------------------
    # 4. Outputs, converted from mg/L to ng/mL (x1000).
    # -------------------------------------------------------------------
    Cc           <- (blood / vbl / bp) * 1000   # bedaquiline plasma
    Cc_m2        <- (blood_m2 / vbl) * 1000     # M2 plasma
    c_liver      <- (liver / vliv) * 1000       # bedaquiline liver
    c_lesion     <- lesion * 1000               # bedaquiline lung lesion
    c_lung       <- lung * 1000                 # bedaquiline uninvolved lung
    c_lesion_m2  <- lesion_m2 * 1000            # M2 lung lesion
    c_lung_m2    <- lung_m2 * 1000              # M2 uninvolved lung

    Cc ~ prop(propSd)
  })
}
