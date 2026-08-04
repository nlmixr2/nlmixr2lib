Mehta_2023_pyrazinamide_mpbpk <- function() {
  description <- paste(
    "mPBPK (minimal physiologically based, translational mouse-to-human).",
    "Pyrazinamide in pulmonary tuberculosis patients, with explicit cavitary lung-lesion and",
    "uninvolved-lung site-of-action compartments. This is the framework-qualification model of",
    "Mehta 2023: first-order oral absorption, first-order elimination from a blood compartment",
    "and perfusion-limited lesion / uninvolved-lung compartments driven by penetration ratios.",
    "No lumped peripheral tissue compartment -- the paper reports that adding one did not improve",
    "the fit. It is the only one of the paper's three models whose lesion and uninvolved-lung",
    "predictions could be checked against observed human site-of-action data. Physiological",
    "volumes and flows are computed from body weight (Brown 1997). The lesion and uninvolved-lung",
    "states hold concentrations (mg/L), not amounts.",
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
    depot  = list(analyte = "pyrazinamide", units = "mg", specimen = "administration site", verified = FALSE),
    blood  = list(analyte = "pyrazinamide", units = "mg", specimen = "blood cell", verified = FALSE),
    lesion = list(analyte = "pyrazinamide", units = "mg", specimen = "tumor", verified = FALSE),
    lung   = list(analyte = "pyrazinamide", units = "mg", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives cardiac output Qc = 312 * (WT/70)^0.75 L/h, total lung volume 0.0076 * WT L and",
        "blood-reservoir volume 0.0771 * WT L (0.0514 venous + 0.0257 arterial), and carries the",
        "0.75 allometric exponent on CL relative to a 70 kg reference. The absorption rate is NOT",
        "scaled by body weight within humans: the -0.25 exponent reported in the Table 1 footnote",
        "describes the mouse-to-human translation of ka (0.30 -> 0.05 1/h), and the human value is",
        "used directly. See the vignette Assumptions and deviations section.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 500L,
    disease_state  = paste(
      "Pulmonary tuberculosis with cavitary disease. Virtual population; cavity presence/absence",
      "and cavity size were sampled from the observed TB-PACTS distributions (Methods 2.2).",
      sep = " "
    ),
    dose_range     = "Simulated pyrazinamide 1500 mg oral daily dosing at steady state (Figure 1).",
    notes          = paste(
      "Developed on mouse plasma, lesion and uninvolved-lung data after a single oral 150 mg/kg",
      "dose digitised from the literature (ESM S1), then translated to humans using allometric",
      "exponents of -0.25 for ka and 0.75 for CL. Unlike bedaquiline and pretomanid, human lesion",
      "and uninvolved-lung concentrations ARE available for pyrazinamide (compiled from the",
      "literature and TB-PACTS), and the translated model reproduced them; this is the external",
      "qualification that the paper uses to justify applying the same framework to bedaquiline",
      "and pretomanid (Results 3.1). R_le and R_ul were estimated imprecisely (RSE 52% and 102%).",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Absorption and disposition -- HUMAN column of Table 1. The mouse
    # estimates were ka = 0.30 1/h and CL = 0.014 L/h.
    # ---------------------------------------------------------------------
    lka <- log(0.05) ; label("First-order oral absorption rate (1/h)")        # Table 1 pyrazinamide ka human = 0.05 (RSE 7.25%), footnote a: 'Allometrically scaled from mice to humans using an exponent of -0.25'
    lcl <- log(3.5)  ; label("Apparent clearance from blood (L/h, 70 kg)")    # Table 1 pyrazinamide CL human = 3.5 (RSE 3.04%)

    # ---------------------------------------------------------------------
    # Penetration ratios -- assumed identical in mice and humans
    # (Methods 2.2). Both were estimated with poor precision.
    # ---------------------------------------------------------------------
    lk_lesion <- log(1.37) ; label("Pyrazinamide lung-lesion penetration ratio, R_le (unitless)")     # Table 1 pyrazinamide R_le = 1.37 (RSE 52.2%)
    lk_lung   <- log(0.85) ; label("Pyrazinamide uninvolved-lung penetration ratio, R_ul (unitless)") # Table 1 pyrazinamide R_ul = 0.85 (RSE 102%)

    # ---------------------------------------------------------------------
    # Fixed physicochemical / scaling constants.
    # ---------------------------------------------------------------------
    lbp     <- fixed(log(0.79))   ; label("Blood-to-plasma concentration ratio, B:P (unitless)")  # Table 1 pyrazinamide BP = 0.79 (reference 42)
    e_wt_cl <- fixed(0.75)        ; label("Allometric exponent on CL (unitless)")                 # Methods 2.2 (reference 20); Results 3.1 'using allometric exponents of -0.25 for ka and 0.75 for CL'
    lvlef   <- fixed(log(0.0216)) ; label("Cavitary lesion volume as a fraction of total lung volume, VF_le (unitless)") # Methods 2.1.1: VF_le assumed 0.0216 from a mean total lesion volume of ~14 mL (reference 19)

    # ---------------------------------------------------------------------
    # Between-subject variability -- 40% log-normal on ka and CL for the
    # human simulations. Variance = 0.40^2 on the log scale.
    # ---------------------------------------------------------------------
    etalka ~ 0.16  # Table 1 pyrazinamide 'IIV for ka and CL (%) = 40'
    etalcl ~ 0.16  # Table 1 pyrazinamide 'IIV for ka and CL (%) = 40'

    # ---------------------------------------------------------------------
    # Residual error. Estimated on the MOUSE fit; the paper states 'Residual
    # errors were not included in the human simulations' (Table 1 footnote).
    # ---------------------------------------------------------------------
    propSd <- 0.35 ; label("Proportional residual error (fraction)")  # Table 1 footnote: 'pyrazinamide combined plasma, lungs, and lesions = proportional 35%'
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual drug-specific parameters
    # -------------------------------------------------------------------
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    bp   <- exp(lbp)
    vlef <- exp(lvlef)

    # -------------------------------------------------------------------
    # 2. Human physiology from body weight (Brown 1997 = ref 15), following
    #    the same construction as the bedaquiline and pretomanid ESM S2 code.
    # -------------------------------------------------------------------
    qc  <- 312 * (WT / 70)^0.75   # cardiac output; Brown 1997 Table 22
    vlu <- 0.0076 * WT            # total lung volume; Brown 1997 Table 7
    vbl <- (0.0514 + 0.0257) * WT # blood reservoir = venous + arterial

    # Lesion and uninvolved-lung volumes and exchange rate constants
    # k_i = Qc / V_i (main-text equation 2).
    vle <- vlef * vlu
    vul <- vlu - vle
    kle <- qc / vle
    kul <- qc / vul

    cbld <- blood / vbl

    # -------------------------------------------------------------------
    # 3. ODE system (ESM S3 panel A). Amounts in mg for depot and blood;
    #    the lesion and lung states hold CONCENTRATIONS in mg/L.
    # -------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(blood) <- ka * depot - cl * blood / vbl

    # Site of action -- main-text equation 1: d/dt(Ci) = ki * (Cbld * Ri - Ci)
    d/dt(lesion) <- kle * (cbld * exp(lk_lesion) - lesion)
    d/dt(lung)   <- kul * (cbld * exp(lk_lung) - lung)

    # -------------------------------------------------------------------
    # 4. Outputs, converted from mg/L to ng/mL (x1000).
    # -------------------------------------------------------------------
    Cc       <- (cbld / bp) * 1000   # pyrazinamide plasma
    c_lesion <- lesion * 1000        # pyrazinamide lung lesion
    c_lung   <- lung * 1000          # pyrazinamide uninvolved lung

    Cc ~ prop(propSd)
  })
}
