Themans_2019_meropenem <- function() {
  description <- "Three-compartment population PK model for meropenem in adults with severe pneumonia, with parallel ELF (epithelial lining fluid) sampling (Themans 2019), as packaged in DDMORE Foundation Model Repository entry DDMODEL00000301."
  reference <- paste(
    "Themans P, Winkin J J, Musuamba F T (2019).",
    "Towards a generic tool for prediction of meropenem systemic and infection-site exposure:",
    "a physiologically based pharmacokinetic model for adult patients with pneumonia.",
    "Drugs R D 19(4):339-355.",
    "doi:10.1007/s40268-019-0268-x.",
    "DDMORE Foundation Model Repository: DDMODEL00000301.",
    sep = " "
  )
  vignette <- "Themans_2019_meropenem"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  ddmore_id    <- "DDMODEL00000301"
  replicate_of <- NULL

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Power-form effect on V1 (central) with reference 75 kg",
        "(per .mod equation TVV1 = THETA(3)*((WT/75)**THETA(4))) and on V2 (the apparent ELF",
        "compartment) with reference 75 kg (TVV2 = THETA(6)*((WT/75)**THETA(9)))."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Renal function expressed as raw, measured glomerular filtration rate (not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Power-form effect on CL with reference 65 mL/min per the .mod equation",
        "TVCL = THETA(1)*((GFR/65)**THETA(2)).",
        "Deviation from canonical CRCL: the canonical register entry is BSA-normalized (mL/min/1.73 m^2);",
        "the Themans 2019 model uses raw measured GFR in mL/min (per the .mod $INPUT comment 'GFR [mL/min]')",
        "and the GFR exponent (0.722) was estimated under that raw-mL/min parameterization.",
        "Source data column 'GFR' (or 'GFR_valAbs' in the bundle's Simulated_dataset.csv) is mapped to canonical 'CRCL' on input."
      ),
      source_name        = "GFR"
    )
  )

  population <- list(
    n_subjects     = "Not extractable from DDMORE bundle (Themans 2019 PDF not on disk).",
    n_studies      = "Not extractable from DDMORE bundle.",
    age_range      = "Adult (specific range not extractable from bundle; the publication is reported in the bundle as covering adults with severe pneumonia).",
    weight_range   = "45-128 kg (range observed in DDMODEL00000301 Simulated_dataset.csv; n = 60 simulated subjects).",
    weight_median  = "approx. 78 kg (mean of DDMODEL00000301 Simulated_dataset.csv).",
    sex_female_pct = "Not extractable from DDMORE bundle.",
    race_ethnicity = "Not extractable from DDMORE bundle.",
    disease_state  = "Adult patients with severe pneumonia (per DDMORE Model_Accomodations.text).",
    dose_range     = paste(
      "1 g meropenem IV every 8 h at steady state in the bundle's Simulated_dataset.csv",
      "(AMT = 1 g, RATE = 2 g/h => 0.5 h infusion duration, II = 8 h, SS = 1).",
      "Clinical dosing in the indication is typically 1-2 g IV q8h."
    ),
    crcl_range     = "19-401 mL/min (range observed in DDMODEL00000301 Simulated_dataset.csv; raw GFR, not BSA-normalized).",
    regions        = "Not extractable from DDMORE bundle.",
    notes          = paste(
      "Demographics summarised from the DDMORE bundle's Simulated_dataset.csv (60 simulated subjects;",
      "WT 45-128 kg, GFR 19-401 mL/min). The Themans 2019 PDF is not on disk under the literature tree;",
      "full demographics, study design, sex / race / ethnicity, and inclusion criteria could not be cross-checked.",
      "The DDMORE Model_Accomodations.text references a British Journal of Pharmacology submission",
      "(July 2019, scenario = 4); the task metadata names a Drugs in R&D 2019 publication",
      "(DOI 10.1007/s40268-019-0268-x). The two are likely the same body of work at different stages of",
      "journal review; no external check was possible."
    )
  )

  ini({
    # Structural typical values; from DDMODEL00000301 Output_real_merop_PK_run3.lst
    # FINAL PARAMETER ESTIMATE block (post MINIMIZATION SUCCESSFUL, OBJV = 1488.719).
    # The .mod $THETA values in Executable_merop_PK_run3.mod are initial estimates and
    # must not be used as final values; the .lst's own embedded $THETA differs from the .mod's.
    lcl  <- log(7.94); label("Clearance for GFR = 65 mL/min, typical individual (CL, L/h)")  # DDMODEL00000301 Output_real .lst FINAL THETA(1)
    lvc  <- log(13.6); label("Central volume of distribution for WT = 75 kg, typical individual (V1, L)")  # DDMODEL00000301 .lst FINAL THETA(3)
    lq   <- log(6.73); label("Inter-compartmental clearance to ELF compartment, typical individual (Q2, L/h)")  # DDMODEL00000301 .lst FINAL THETA(5)
    lvp  <- log(4.08); label("Apparent volume of ELF compartment for WT = 75 kg, typical individual (V2, L)")  # DDMODEL00000301 .lst FINAL THETA(6)
    lq2  <- log(8.22); label("Inter-compartmental clearance to deep peripheral compartment, typical individual (Q3, L/h)")  # DDMODEL00000301 .lst FINAL THETA(7)
    lvp2 <- log(10.1); label("Volume of deep peripheral compartment, typical individual (V3, L)")  # DDMODEL00000301 .lst FINAL THETA(8)

    # Covariate effects (power form, all unitless exponents)
    e_crcl_cl <- 0.722; label("CRCL (GFR) effect on CL: exponent of (CRCL / 65) (unitless)")  # DDMODEL00000301 .lst FINAL THETA(2)
    e_wt_vc   <- 0.949; label("WT effect on V1 (central): exponent of (WT / 75) (unitless)")  # DDMODEL00000301 .lst FINAL THETA(4)
    e_wt_vp   <- 1.04;  label("WT effect on V2 (ELF compartment): exponent of (WT / 75) (unitless)")  # DDMODEL00000301 .lst FINAL THETA(9)

    # ELF / V2 partition factor (dimensionless). The source NONMEM model parameterised
    # this as THETA(10) = V2 / S2 = 249 under a mixed-unit convention where dose was
    # in g and concentration in mg/L (S1 = V1/1000 absorbs the g -> mg conversion). The
    # equivalent dimensionless penetration factor in an mg-dose, mg/L-concentration
    # convention is f_elf = THETA(10) / 1000 ~= 0.249, so ELF concentration is computed
    # as Celf = (peripheral1 / vp) * f_elf, equivalent to NONMEM's Celf = A2 * THETA(10) / V2
    # under its mixed-unit convention. See vignette Errata for the unit-trace derivation.
    lf_elf <- log(0.249); label("log meropenem ELF/V2 partition factor (unitless)")  # DDMODEL00000301 .lst FINAL THETA(10) = 249, rescaled by 1/1000 for mg-consistent units

    # Inter-individual variability. The .mod declared $OMEGA BLOCK(4) but every off-diagonal
    # was held at 0 (see both the .mod and the FINAL OMEGA matrix in the .lst), so this
    # maps to four independent etas in nlmixr2.
    etalcl ~ 0.126   # DDMODEL00000301 .lst FINAL OMEGA ETA1 (variance, ETA on CL)
    etalvc ~ 0.140   # DDMODEL00000301 .lst FINAL OMEGA ETA2 (variance, ETA on V1)
    etalvp ~ 1.76    # DDMODEL00000301 .lst FINAL OMEGA ETA3 (variance, ETA on V2)
    etalq2 ~ 0.187   # DDMODEL00000301 .lst FINAL OMEGA ETA4 (variance, ETA on Q3)

    # Multi-output residual error.
    # Plasma  (CMT 1, central):     Y = F * (1 + EPS1) + EPS2  -> combined proportional + additive on Cc.
    # ELF     (CMT 2, peripheral1): Y = F * (1 + EPS3)         -> proportional only on Celf.
    # SIGMA in the .lst is reported as variance; convert to SD for nlmixr2 prop()/add() arguments.
    propSd      <- sqrt(0.0240); label("Proportional residual error on plasma concentration Cc (fraction)")  # DDMODEL00000301 .lst FINAL SIGMA EPS1 (variance)
    addSd       <- sqrt(0.208);  label("Additive residual error on plasma concentration Cc (mg/L)")           # DDMODEL00000301 .lst FINAL SIGMA EPS2 (variance)
    propSd_Celf <- sqrt(0.404);  label("Proportional residual error on ELF concentration Celf (fraction)")   # DDMODEL00000301 .lst FINAL SIGMA EPS3 (variance)
  })
  model({
    # Individual PK parameters with power-form covariate effects (per DDMODEL00000301 $PK block:
    #   TVCL = THETA(1) * (GFR/65)^THETA(2)
    #   TVV1 = THETA(3) * (WT/75)^THETA(4)
    #   TVQ2 = THETA(5)
    #   TVV2 = THETA(6) * (WT/75)^THETA(9)
    #   TVQ3 = THETA(7)
    #   TVV3 = THETA(8)
    # IIV is on CL (ETA1), V1 (ETA2), V2 (ETA3), and Q3 (ETA4) only; Q2 and V3 have no IIV.
    cl  <- exp(lcl  + etalcl) * (CRCL / 65)^e_crcl_cl
    vc  <- exp(lvc  + etalvc) * (WT   / 75)^e_wt_vc
    q   <- exp(lq)
    vp  <- exp(lvp  + etalvp) * (WT   / 75)^e_wt_vp
    q2  <- exp(lq2  + etalq2)
    vp2 <- exp(lvp2)

    # ELF / V2 dimensionless penetration factor; inherits V2's WT scaling and IIV through vp below.
    f_elf <- exp(lf_elf)

    # Micro-constants for the 3-compartment IV ODE system (NONMEM ADVAN11 / TRANS4).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ODE system. central = plasma compartment (V1) where dosing occurs and plasma samples
    # are drawn (CMT = 1 in the source). peripheral1 = ELF compartment (V2) where ELF
    # samples are drawn (CMT = 2 in the source). peripheral2 = deeper peripheral (V3),
    # not directly sampled.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                       k13 * central - k31 * peripheral2

    # Plasma concentration (mg/L). Dose AMT in mg gives A_central in mg, divided by V1 in L.
    Cc   <- central / vc
    # ELF concentration (mg/L) = V2-equivalent concentration scaled by the meropenem ELF/V2
    # penetration factor. This reproduces NONMEM's Celf = A2 * THETA(10) / V2 in unit-consistent
    # mg/L when dose is supplied in mg (see notes in ini() block on f_elf).
    Celf <- (peripheral1 / vp) * f_elf

    Cc   ~ add(addSd) + prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
