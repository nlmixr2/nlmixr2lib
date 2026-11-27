Jiang_2019_iopamidol <- function() {
  description <- paste(
    "Renal tracer-kinetic model of the iodinated contrast agent iopamidol measured by dynamic",
    "contrast-enhanced computed tomography (DCE-CT) in adults with essential hypertension or",
    "atherosclerotic renal artery stenosis. A modified two-compartment kidney model driven by",
    "two measured forcing functions: the aortic arterial input function (converted to plasma by",
    "an assumed haematocrit of 45%) enters the renal vascular space after a blood transit delay",
    "and with a first-order perfusion rate constant; the vascular contrast is filtered into the",
    "renal tubules with a normalized-GFR rate constant and leaves the tubules in proportion to",
    "the measured tubular (filtration) component of the renal papilla curve. The observed",
    "whole-kidney CT enhancement is the vascular-volume-fraction-weighted vascular",
    "concentration plus the tubular concentration. Single-kidney GFR, perfusion and renal blood",
    "flow are derived outputs. Parameters are the per-kidney-group summary values of",
    "individual least-squares fits (essential-hypertension, stenotic and contralateral kidneys);",
    "no between-kidney random effects or residual-error model were reported.",
    sep = " "
  )
  reference <- paste(
    "Jiang K, Ferguson CM, Abumoawad A, Saad A, Textor SC, Lerman LO. (2019).",
    "A modified two-compartment model for measurement of renal function using dynamic",
    "contrast-enhanced computed tomography. PLoS One 14(7):e0219605.",
    "doi:10.1371/journal.pone.0219605.",
    sep = " "
  )
  vignette <- "Jiang_2019_iopamidol"

  units <- list(
    time = "s",
    dosing = paste(
      "n/a (no dosing events; the kidney is driven by the measured aortic enhancement",
      "CP_IOPAMIDOL_HU and papilla tubular enhancement CPAPILLA_TUBULAR_IOPAMIDOL_HU)"
    ),
    concentration = paste(
      "HU (increase in CT attenuation above the pre-contrast baseline, used as a linear",
      "surrogate of iodine concentration)"
    )
  )

  covariateData <- list(
    CP_IOPAMIDOL_HU = list(
      description = paste(
        "Arterial (abdominal-aorta) whole-blood iodine enhancement after an iopamidol bolus:",
        "the aortic region-of-interest CT number minus its pre-contrast baseline"
      ),
      units = "HU",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying error-free regressor (the measured arterial input function); declared",
        "linear() so rxode2 interpolates linearly between scans, as the authors resampled the",
        "curves to 0.5 s by linear interpolation. The model converts it to plasma as",
        "CP_IOPAMIDOL_HU / (1 - hct) (Jiang 2019 Methods: Cp(t) = Ca(t)/(1 - Hct)). 0 before",
        "contrast arrival."
      ),
      source_name = "Ca(t)"
    ),
    CPAPILLA_TUBULAR_IOPAMIDOL_HU = list(
      description = paste(
        "Tubular-transit (filtration) component of the renal-papilla iodine enhancement,",
        "Cout''(t): the papilla CT enhancement minus its perfusion component"
      ),
      units = "HU",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying error-free regressor declared linear(). Drives the tubular outflow term",
        "kout * Cout''(t) of Jiang 2019 Eq. 2. The authors obtained it by fitting the perfusion",
        "component Cout'(t) of the papilla curve with Eq. 1 and subtracting it (Methods,",
        "Compartment model; Fig 2c left). 0 before filtered contrast reaches the papilla."
      ),
      source_name = "Cout''(t)"
    ),
    DIS_RAS_STENOTIC = list(
      description = paste(
        "1 = the kidney is supplied by a renal artery with atherosclerotic stenosis; 0 = otherwise.",
        "Kidney-level indicator"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (kidney of a patient with essential hypertension when DIS_RAS_CONTRALATERAL is also 0)",
      notes = paste(
        "Selects the 'Stenotic Kidney' column of Jiang 2019 Table 2. Both kidneys of a patient with",
        "bilateral stenosis are stenotic (5 of the 25 ARAS patients). Mutually exclusive with",
        "DIS_RAS_CONTRALATERAL."
      ),
      source_name = "Stenotic kidney"
    ),
    DIS_RAS_CONTRALATERAL = list(
      description = paste(
        "1 = the non-stenotic contralateral kidney of a patient with unilateral atherosclerotic",
        "renal artery stenosis; 0 = otherwise. Kidney-level indicator"
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (kidney of a patient with essential hypertension when DIS_RAS_STENOTIC is also 0)",
      notes = paste(
        "Selects the 'Contralateral Kidney' column of Jiang 2019 Table 2. Mutually exclusive with",
        "DIS_RAS_STENOTIC."
      ),
      source_name = "Contralateral kidney"
    ),
    ORGVOL_KIDNEY = list(
      description = "Volume of the modelled (single) kidney measured on helical CT",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "SINGLE-kidney volume, not the combined volume of both kidneys. Used only by the derived",
        "outputs: GFR = kGFR * V * 60 (Jiang 2019 Eq. 3) and renal blood flow = perfusion * V / 100.",
        "Table 2 group means: 148 mL (essential hypertension), 86 mL (stenotic), 170 mL",
        "(contralateral)."
      ),
      source_name = "V"
    )
  )

  compartmentData <- list(
    vp_kidney = list(
      analyte = "iopamidol",
      units = "HU",
      specimen = "plasma",
      verified = TRUE
    ),
    kidney = list(
      analyte = "iopamidol",
      units = "HU",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 38,
    n_studies = 1,
    n_kidneys = 76,
    age_range = "median 67.0 (IQR 62.0-71.0) years essential hypertension; 70.0 (58.0-74.0) years renal artery stenosis",
    weight_range = "BMI 28.5 +/- 3.7 (essential hypertension) and 27.6 +/- 4.0 (renal artery stenosis) kg/m2",
    sex_female_pct = 42.1,
    disease_state = paste(
      "Hypertension: essential hypertension (n = 13 patients, 26 kidneys) or renovascular",
      "hypertension from atherosclerotic renal artery stenosis (n = 25 patients; 30 stenotic and",
      "20 contralateral kidneys, 5 patients bilateral). Serum creatinine below 2.5 mg/dL."
    ),
    dose_range = "Single central venous bolus of iopamidol-370 0.5 mL/kg before dynamic CT",
    regions = "United States (Mayo Clinic, Rochester, MN)",
    renal_function = paste(
      "Iothalamate GFR median 89.2 (IQR 81.0-99.5) mL/min essential hypertension and 60.1",
      "(47.0-90.3) mL/min renal artery stenosis (Table 1)"
    ),
    co_medication = "ACE inhibitor or angiotensin receptor blocker; fixed 150 mEq/day sodium diet",
    notes = paste(
      "Jiang 2019 Table 1 (baseline characteristics). Sex: 7/13 and 15/25 male, so 16 of 38",
      "(42.1%) female. Every kidney was fitted separately (least squares, Levenberg-Marquardt,",
      "0-120 s of the dynamic scan); the packaged values are the per-group summaries of those fits."
    )
  )

  ini({
    # Jiang 2019 Table 2 'Kidney volume and model-fitted parameters'. Each kidney was fitted
    # individually; Table 2 summarises the fits per kidney group as mean +/- SD (f, kGFR) or
    # median (IQR) (kperf, Td, kout). The summary value is held fixed here; the spread across
    # kidneys is quoted in the comment and carries no random effect.

    # Essential-hypertension (EH) kidneys -- reference group
    kperf_kidney_eh <- fixed(0.15)
    label("Renal perfusion rate constant kperf, essential-hypertension kidney (1/s)") # Table 2 'kperf (s-1)' EH Kidney = 0.15 (0.11-0.19)
    tdelay_kidney_eh <- fixed(0.68)
    label("Aorta-to-kidney blood transit delay Td, essential-hypertension kidney (s)") # Table 2 'Td (s)' EH Kidney = 0.68 (0.34-1.16)
    fvasc_kidney_eh <- fixed(0.22)
    label("Renal vascular volume fraction f, essential-hypertension kidney (fraction)") # Table 2 'f' EH Kidney = 0.22 +/- 0.04
    kgfr_kidney_eh <- fixed(0.0048)
    label("Normalized GFR rate constant kGFR, essential-hypertension kidney (1/s)") # Table 2 'kGFR (s-1)' EH Kidney = 0.0048 +/- 0.0009
    kout_tubule_eh <- fixed(0.004)
    label("Tubular iodine outflow rate constant kout, essential-hypertension kidney (1/s)") # Table 2 'kout (s-1)' EH Kidney = 0.004 (0.003-0.011)

    # Stenotic kidneys (atherosclerotic renal artery stenosis)
    kperf_kidney_sten <- fixed(0.12)
    label("Renal perfusion rate constant kperf, stenotic kidney (1/s)") # Table 2 'kperf (s-1)' Stenotic Kidney = 0.12 (0.09-0.17)
    tdelay_kidney_sten <- fixed(0.58)
    label("Aorta-to-kidney blood transit delay Td, stenotic kidney (s)") # Table 2 'Td (s)' Stenotic Kidney = 0.58 (0.05-1.35)
    fvasc_kidney_sten <- fixed(0.19)
    label("Renal vascular volume fraction f, stenotic kidney (fraction)") # Table 2 'f' Stenotic Kidney = 0.19 +/- 0.06
    kgfr_kidney_sten <- fixed(0.0035)
    label("Normalized GFR rate constant kGFR, stenotic kidney (1/s)") # Table 2 'kGFR (s-1)' Stenotic Kidney = 0.0035 +/- 0.0013
    kout_tubule_sten <- fixed(0.005)
    label("Tubular iodine outflow rate constant kout, stenotic kidney (1/s)") # Table 2 'kout (s-1)' Stenotic Kidney = 0.005 (0.003-0.012)

    # Contralateral kidneys of unilateral renal artery stenosis
    kperf_kidney_contra <- fixed(0.17)
    label("Renal perfusion rate constant kperf, contralateral kidney (1/s)") # Table 2 'kperf (s-1)' Contralateral Kidney = 0.17 (0.11-0.18)
    tdelay_kidney_contra <- fixed(1.02)
    label("Aorta-to-kidney blood transit delay Td, contralateral kidney (s)") # Table 2 'Td (s)' Contralateral Kidney = 1.02 (0.47-1.57)
    fvasc_kidney_contra <- fixed(0.18)
    label("Renal vascular volume fraction f, contralateral kidney (fraction)") # Table 2 'f' Contralateral Kidney = 0.18 +/- 0.04
    kgfr_kidney_contra <- fixed(0.0041)
    label("Normalized GFR rate constant kGFR, contralateral kidney (1/s)") # Table 2 'kGFR (s-1)' Contralateral Kidney = 0.0041 +/- 0.0011
    kout_tubule_contra <- fixed(0.007)
    label("Tubular iodine outflow rate constant kout, contralateral kidney (1/s)") # Table 2 'kout (s-1)' Contralateral Kidney = 0.007 (0.002-0.009)

    # Fixed physiological constants (not estimated)
    hct <- fixed(0.45)
    label("Assumed haematocrit converting arterial blood to plasma enhancement (fraction)") # Methods, Compartment model: 'Hct the assumed hematocrit (45%)'
    pc_water_kidney <- fixed(0.9)
    label("Blood/tissue partition coefficient for water sigma in the perfusion formula (mL/g)") # Methods, Eq. 4: 'sigma ... blood/tissue partition coefficient for water, which is 0.9'

    # No residual-error model was reported (unweighted least squares on 0-120 s); a zero
    # additive SD keeps the model a typical-value predictor of the kidney enhancement.
    addSd <- fixed(0)
    label("Additive residual SD of the kidney CT enhancement (HU); not reported") # not reported (least-squares fit, Methods 'Image analysis')
  })

  model({
    # Measured forcing functions, interpolated linearly between scan times.
    linear(CP_IOPAMIDOL_HU)
    linear(CPAPILLA_TUBULAR_IOPAMIDOL_HU)

    # Kidney-group selection; exactly one weight is 1 for a valid covariate combination.
    sel_eh <- 1 - DIS_RAS_STENOTIC - DIS_RAS_CONTRALATERAL
    kperf_kidney <- sel_eh * kperf_kidney_eh + DIS_RAS_STENOTIC * kperf_kidney_sten +
      DIS_RAS_CONTRALATERAL * kperf_kidney_contra
    tdelay_kidney <- sel_eh * tdelay_kidney_eh + DIS_RAS_STENOTIC * tdelay_kidney_sten +
      DIS_RAS_CONTRALATERAL * tdelay_kidney_contra
    fvasc_kidney <- sel_eh * fvasc_kidney_eh + DIS_RAS_STENOTIC * fvasc_kidney_sten +
      DIS_RAS_CONTRALATERAL * fvasc_kidney_contra
    kgfr_kidney <- sel_eh * kgfr_kidney_eh + DIS_RAS_STENOTIC * kgfr_kidney_sten +
      DIS_RAS_CONTRALATERAL * kgfr_kidney_contra
    kout_tubule <- sel_eh * kout_tubule_eh + DIS_RAS_STENOTIC * kout_tubule_sten +
      DIS_RAS_CONTRALATERAL * kout_tubule_contra

    # Arterial input function in plasma, Cp(t) = Ca(t) / (1 - Hct).
    cp_aif <- CP_IOPAMIDOL_HU / (1 - hct)

    # Eq. 1, dCv/dt = kperf * (Cp(t - Td) - Cv). The equation is linear and time-invariant
    # with Cv(0) = 0, so Cv(t) equals the undelayed solution evaluated at t - Td. The state
    # vp_kidney is that undelayed solution (the renal vascular plasma enhancement on the aortic
    # clock), and delay() shifts it by the transit delay; rxode2 cannot delay a covariate
    # directly. delay() returns the initial value 0 for t < Td.
    d/dt(vp_kidney) <- kperf_kidney * (cp_aif - vp_kidney)
    vp_kidney(0) <- 0
    cv_kidney <- delay(vp_kidney, tdelay_kidney)

    # Eq. 2, dCt/dt = kGFR * Cv(t) - kout * Cout''(t): filtration into the tubules and outflow
    # proportional to the measured tubular component of the papilla curve.
    d/dt(kidney) <- kgfr_kidney * cv_kidney - kout_tubule * CPAPILLA_TUBULAR_IOPAMIDOL_HU
    kidney(0) <- 0

    # Derived single-kidney function. Eq. 3: GFR (mL/min) = kGFR (1/s) * V (mL) * 60.
    # Eq. 4: perfusion (mL/100 g/min) = kperf * f * 60 * 100 / (1 - Hct) * sigma.
    # Renal blood flow (mL/min) = perfusion * V / 100 (Methods, after Eq. 4).
    gfr <- kgfr_kidney * ORGVOL_KIDNEY * 60
    perfusion_kidney <- kperf_kidney * fvasc_kidney * 60 * 100 / (1 - hct) * pc_water_kidney
    rbf_kidney <- perfusion_kidney * ORGVOL_KIDNEY / 100

    # Whole-kidney enhancement, Ctot(t) = Cv(t) * f + Ct(t) (Methods, after Fig 2).
    Cc <- fvasc_kidney * cv_kidney + kidney
    Cc ~ add(addSd)
  })
}
