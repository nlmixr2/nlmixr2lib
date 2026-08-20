Mehta_2023_pretomanid_mpbpk <- function() {
  description <- paste(
    "mPBPK (minimal physiologically based, translational mouse-to-human).",
    "Pretomanid in pulmonary tuberculosis patients, with explicit cavitary lung-lesion and",
    "uninvolved-lung site-of-action compartments. First-order oral absorption with saturable,",
    "dose-dependent bioavailability F = 1 / (1 + DOSE/ED50), linear clearance from blood, two",
    "lumped peripheral tissue pools that split cardiac output and residual body volume by a",
    "fixed fraction, and perfusion-limited lesion / uninvolved-lung compartments driven by",
    "penetration ratios. Physiological volumes and flows are computed from body weight",
    "(Brown 1997); clearance is allometrically scaled from the mouse fit with an exponent of",
    "0.75. The lesion and uninvolved-lung states hold concentrations (mg/L), not amounts.",
    "Intended for typical-value and population simulation of site-of-action target attainment.",
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
    depot       = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = FALSE),
    blood       = list(analyte = "pretomanid", units = "mg", specimen = "blood cell", verified = FALSE),
    peripheral1 = list(analyte = "pretomanid", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "pretomanid", units = "mg", specimen = "plasma", verified = FALSE),
    lesion      = list(analyte = "pretomanid", units = "mg", specimen = "tissue", verified = FALSE),
    lung        = list(analyte = "pretomanid", units = "mg", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives cardiac output Qc = 312 * (WT/70)^0.75 L/h, total lung volume 0.0076 * WT L,",
        "blood-reservoir volume 0.0771 * WT L (0.0514 venous + 0.0257 arterial) and the residual",
        "body volume that is split between the two lumped peripheral pools. Also carries the 0.75",
        "allometric exponent on CL relative to a 70 kg reference. Body weights for the 500 virtual",
        "patients were sampled from the TB-PACTS clinical-trial weight distribution",
        "(Methods 2.2, ESM S1).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    DOSE_PTM_MG = list(
      description        = "Pretomanid dose amount for the current administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose-record covariate driving the saturable bioavailability",
        "F = Fmax / (1 + DOSE/ED50) with Fmax assumed to be 1 (Table 1 ED50 description).",
        "Must equal the amt of the corresponding dosing record; a 200 mg dose therefore gives",
        "F = 1 / (1 + 200/554) = 0.735. ESM S2 writes this as",
        "doseIn = Fmax*dose/(1 + dose/ED50) followed by f(depot) = doseIn/dose, which is",
        "algebraically identical to the form used here and avoids a division by the dose amount.",
        sep = " "
      ),
      source_name        = "dose"
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
      "Simulated standard pretomanid dosing of 200 mg once daily (Methods 2.5).",
      "Clinical validation data spanned 50-1200 mg (ESM S1).",
      sep = " "
    ),
    notes          = paste(
      "The structural model and every drug-specific parameter were estimated on mouse data:",
      "plasma after single oral doses of 6-486 mg/kg, plus plasma, lesion and uninvolved-lung",
      "PET-imaging data after intravenous F18-pretomanid (ESM S1). The model was then translated",
      "to humans by swapping in human physiology and allometrically scaling clearance; the human",
      "ED50 of 554 mg is a literature value (reference 41) rather than the mouse estimate of",
      "7.59 mg. Human plasma predictions agreed reasonably with TB-patient data, although the",
      "paper notes under-prediction at 1000 and 1200 mg (Results 3.2). No human lesion or",
      "uninvolved-lung PK data exist for pretomanid.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Absorption and saturable bioavailability.
    # ---------------------------------------------------------------------
    lka   <- log(0.3)   ; label("First-order oral absorption rate (1/h)")                 # Table 1 pretomanid ka = 0.3 (RSE 23.3%), same value mice and humans; ESM S2 TVka = 0.3
    led50 <- log(554)   ; label("Dose at which bioavailability is half of maximal, ED50 (mg)") # Table 1 pretomanid ED50 human = 554 (RSE 14.8%), literature value (reference 41); ESM S2 TVED50 = 554

    # ---------------------------------------------------------------------
    # Disposition. Human clearance from Table 1 (mouse estimate 0.016 L/h
    # allometrically scaled).
    # ---------------------------------------------------------------------
    lcl   <- log(4.42)  ; label("Apparent clearance from blood (L/h, 70 kg)")             # Table 1 pretomanid CL human = 4.42 (RSE 3.78%); ESM S2 TVCL = 4.42

    # ---------------------------------------------------------------------
    # Tissue:blood equilibrium concentration ratios and the peripheral split.
    # The paper's 'tissue compartment 1' (KpT) is peripheral1 and its 'tissue
    # compartment 2' (KpT1) is peripheral2; FT1 is the fraction of cardiac
    # output and of residual body volume assigned to peripheral2.
    # ---------------------------------------------------------------------
    lk_peripheral1   <- log(36.24) ; label("Pretomanid peripheral1 partition coefficient, KpT (unitless)")  # ESM S2 TVKpT = 36.24 (Table 1 prints the rounded 36.3, RSE 3.77%)
    lk_peripheral2   <- log(0.48)  ; label("Pretomanid peripheral2 partition coefficient, KpT1 (unitless)") # ESM S2 TVKpT1 = 0.48 (Table 1 prints 0.483, RSE 14%)
    frac_peripheral2 <- 0.975      ; label("Fraction of cardiac output and residual volume assigned to peripheral2, FT1 (unitless)") # ESM S2 FT1 = 0.975 (Table 1 prints the rounded 0.97, RSE 9.72%)
    lk_lesion        <- log(1.6)   ; label("Pretomanid lung-lesion penetration ratio, R_le (unitless)")     # ESM S2 Rles = 1.6. Table 1 prints 1.05 (RSE 145%), which cannot reproduce the paper's own reported lesion-to-plasma ratio of 2.6; see vignette Errata.
    lk_lung          <- log(1.76)  ; label("Pretomanid uninvolved-lung penetration ratio, R_ul (unitless)") # ESM S2 RUL = 1.76 (Table 1 prints the rounded 1.75, RSE 13.8%)

    # ---------------------------------------------------------------------
    # Fixed physicochemical / scaling constants.
    # ---------------------------------------------------------------------
    lbp     <- fixed(log(1.65))   ; label("Blood-to-plasma concentration ratio, B:P (unitless)")  # Table 1 pretomanid BP = 1.65 (reference 13); ESM S2 BP = 1.65
    e_wt_cl <- fixed(0.75)        ; label("Allometric exponent on CL (unitless)")                 # Methods 2.2 (reference 20); ESM S2 (BW/70)^0.75
    lvlef   <- fixed(log(0.0216)) ; label("Cavitary lesion volume as a fraction of total lung volume, VF_le (unitless)") # Methods 2.1.1: VF_le assumed 0.0216 from a mean total lesion volume of ~14 mL (reference 19)

    # ---------------------------------------------------------------------
    # Between-subject variability -- 40% log-normal (Table 1). ESM S2 applies
    # exp(eta) to ka, CL, KpT and KpT1; the ED50 eta is deliberately commented
    # out ('ED50 = TVED50 ; #*exp(eta.ED50) ;'), so ED50 carries no IIV here
    # (operator decision, 2026-07-27). Variance = 0.40^2 on the log scale.
    # ---------------------------------------------------------------------
    etalka             ~ 0.16  # Table 1 pretomanid 'IIV for ka, ED50, CL, KpT, and KpT2 (%) = 40'; ESM S2 ka = TVka*exp(eta.ka)
    etalcl             ~ 0.16  # Table 1 pretomanid IIV row = 40%; ESM S2 CL = TVCL*exp(eta.CL)*(BW/70)^0.75
    etalk_peripheral1  ~ 0.16  # Table 1 pretomanid IIV row = 40%; ESM S2 KpT = TVKpT*exp(eta.KpT)
    etalk_peripheral2  ~ 0.16  # Table 1 pretomanid IIV row = 40%; ESM S2 KpT1 = TVKpT1*exp(eta.KpT1)

    # ---------------------------------------------------------------------
    # Residual error. Estimated on the MOUSE fit; the paper states 'Residual
    # errors were not included in the human simulations' (Table 1 footnote).
    # ---------------------------------------------------------------------
    propSd <- 0.12 ; label("Proportional residual error (fraction)")  # Table 1 footnote: 'pretomanid plasma = proportional 12%'
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual drug-specific parameters
    # -------------------------------------------------------------------
    ka     <- exp(lka + etalka)
    ed50   <- exp(led50)
    cl     <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    kpt    <- exp(lk_peripheral1 + etalk_peripheral1)
    kpt1   <- exp(lk_peripheral2 + etalk_peripheral2)
    ft1    <- frac_peripheral2
    bp     <- exp(lbp)
    vlef   <- exp(lvlef)

    # -------------------------------------------------------------------
    # 2. Human physiology from body weight (ESM S2; Brown 1997 = ref 15).
    #    Pretomanid has no explicit liver compartment, so the residual body
    #    volume is body weight minus lung and blood only.
    # -------------------------------------------------------------------
    qc     <- 312 * (WT / 70)^0.75   # cardiac output; Brown 1997 Table 22
    qt     <- qc * (1 - ft1)         # flow to peripheral1
    qt1    <- qc * ft1               # flow to peripheral2
    vlu    <- 0.0076 * WT            # total lung volume; Brown 1997 Table 7
    vbl    <- (0.0514 + 0.0257) * WT # blood reservoir = venous + arterial
    vt_tot <- WT - vlu - vbl         # residual body volume
    vt     <- vt_tot * (1 - ft1)     # peripheral1 volume
    vt1    <- vt_tot * ft1           # peripheral2 volume

    # Lesion and uninvolved-lung volumes and exchange rate constants
    # k_i = Qc / V_i (main-text equation 2).
    vle <- vlef * vlu
    vul <- vlu - vle
    kle <- qc / vle
    kul <- qc / vul

    cbld <- blood / vbl

    # -------------------------------------------------------------------
    # 3. ODE system (ESM S2 'Final Model Codes'). Amounts in mg for
    #    depot/blood/peripheral1/peripheral2; the lesion and lung states
    #    hold CONCENTRATIONS in mg/L.
    # -------------------------------------------------------------------
    d/dt(depot) <- -ka * depot

    d/dt(blood) <- ka * depot - cl * blood / vbl +
      peripheral1 * qt / (vt * kpt) - blood * qt / vbl +
      peripheral2 * qt1 / (vt1 * kpt1) - blood * qt1 / vbl

    d/dt(peripheral1) <- blood * qt / vbl - peripheral1 * qt / (vt * kpt)
    d/dt(peripheral2) <- blood * qt1 / vbl - peripheral2 * qt1 / (vt1 * kpt1)

    # Site of action -- main-text equation 1: d/dt(Ci) = ki * (Cbld * Ri - Ci)
    d/dt(lesion) <- kle * (cbld * exp(lk_lesion) - lesion)
    d/dt(lung)   <- kul * (cbld * exp(lk_lung) - lung)

    # -------------------------------------------------------------------
    # 4. Saturable, dose-dependent bioavailability with Fmax fixed at 1.
    #    ESM S2: doseIn = Fmax*dose/(1 + dose/ED50); f(depot) = doseIn/dose.
    # -------------------------------------------------------------------
    f(depot) <- 1 / (1 + DOSE_PTM_MG / ed50)

    # -------------------------------------------------------------------
    # 5. Outputs, converted from mg/L to ng/mL (x1000).
    # -------------------------------------------------------------------
    Cc       <- (cbld / bp) * 1000   # pretomanid plasma
    c_lesion <- lesion * 1000        # pretomanid lung lesion
    c_lung   <- lung * 1000          # pretomanid uninvolved lung

    Cc ~ prop(propSd)
  })
}
