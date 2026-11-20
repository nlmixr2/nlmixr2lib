Wong_2026_ginisortamab <- function() {
  description <- "Two-compartment population PK/PD model of ginisortamab (anti-gremlin-1 IgG monoclonal antibody) and total serum gremlin-1 in adults with advanced solid tumours (Wong 2026, ONC001 dose-escalation module). Linear clearance of free drug from the central compartment and quasi-steady-state target-mediated drug disposition (TMDD-QSS) binding to serum gremlin-1 with zero-order gremlin-1 synthesis, first-order free gremlin-1 degradation and first-order internalisation of the drug-gremlin-1 complex. Baseline body surface area (power, centred at 2 m^2) on the central and peripheral volumes. Outputs free ginisortamab, total gremlin-1, free gremlin-1 and percent gremlin-1 suppression."
  reference <- "Wong YC, Chavarro LR, Jones P, Egging D, Ungphakorn W, Valle M, Matthews I. A dose-finding population pharmacokinetic/pharmacodynamic model of ginisortamab, an anti-gremlin-1 monoclonal antibody, in patients with solid tumours. Br J Clin Pharmacol. 2026;92(8):2885-2898. doi:10.1002/bcp.70567"
  vignette <- "Wong_2026_ginisortamab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL",
    Ctotal_target = "ng/mL"
  )

  compartmentData <- list(
    central = list(
      analyte = "ginisortamab (total: free + gremlin-1-bound)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "ginisortamab (free)", units = "mg", specimen = "not applicable", verified = TRUE),
    total_target = list(
      analyte = "gremlin-1 (total: free + ginisortamab-bound)",
      units = "nmol/L",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    BSA = list(
      description = "Baseline body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline BSA (control-stream column BSABL), time-fixed per subject. Power scaling (BSA/2)^exponent on Vc (exponent 1.31) and Vp (exponent 2.18), centred at 2 m^2 as in the Supporting Information control stream and the Figure 2 equations. The paper calls 2 m^2 the participant median, although Table S1 reports a median of 1.8 m^2 (range 1.4-2.4); the centring constant of the fitted model is used. BSA formula not stated.",
      source_name = "BSABL"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 1L,
    n_observations = "371 ginisortamab and 270 total gremlin-1 serum concentrations",
    age_range = "median 64 years (range 37-74)",
    weight_range = "median 69 kg (range 46-116)",
    bmi_range = "median 24 kg/m^2 (range 19-36)",
    bsa_range = "median 1.8 m^2 (range 1.4-2.4)",
    sex_female_pct = 32,
    race_ethnicity = c(White = 92, Asian = 4, `Black or African American` = 4),
    disease_state = "Recurrent or metastatic advanced solid tumours associated with GREM1 mRNA expression, resistant or refractory to standard treatment, ECOG PS <= 1: colorectal adenocarcinoma 11, pancreatic adenocarcinoma 6, prostate adenocarcinoma 3, and one each of oesophageal carcinoma, gastric adenocarcinoma, breast invasive carcinoma and lung adenocarcinoma.",
    dose_range = "100, 250, 500, 1000 or 2000 mg as a 60-min IV infusion every 2 weeks (Q2W).",
    regions = "United Kingdom",
    ada = "All participants were ADA-negative at all visits.",
    notes = "ONC001 (NCT04393298) phase 1/2 first-in-human monotherapy dose-escalation module. Demographics from Wong 2026 Table S1; observation counts from Table S2."
  )

  ini({
    # ---- Structural PK (Wong 2026 Table 1, 'Estimated value' column = NONMEM final estimates) ----
    # Reference subject: baseline BSA = 2 m^2. Time unit: day.
    lcl <- log(0.433); label("Clearance of free ginisortamab CL (L/day)") # Table 1, CL 0.433 L/day (RSE 4.7%)
    lvc <- log(4.38); label("Central volume of distribution Vc at BSA 2 m^2 (L)") # Table 1, Vc 4.38 L (RSE 3.7%)
    lq <- log(0.849); label("Inter-compartmental clearance Q (L/day)") # Table 1, Q 0.849 L/day (RSE 12.8%)
    lvp <- log(3.23); label("Peripheral volume of distribution Vp at BSA 2 m^2 (L)") # Table 1, Vp 3.23 L (RSE 7.8%)

    # ---- Covariate effects (power on BSA/2; control stream THETA(9), THETA(10)) ----
    e_bsa_vc <- 1.31; label("Power exponent of BSA/2 on Vc (unitless)") # Table 1, covariate effect of BSA on Vc 1.31 (RSE 14.2%)
    e_bsa_vp <- 2.18; label("Power exponent of BSA/2 on Vp (unitless)") # Table 1, covariate effect of BSA on Vp 2.18 (RSE 21.3%)

    # ---- Gremlin-1 turnover and QSS binding (Table 1; Supporting Information $PK / $DES) ----
    lksyn <- log(0.378); label("Zero-order gremlin-1 synthesis rate Ksyn (nmol/L/day)") # Table 1, Ksyn 0.378 nmol/L/day (RSE 22.2%)
    lkdeg <- log(61.8); label("First-order free gremlin-1 degradation rate constant Kdeg (1/day)") # Table 1, Kdeg 61.8 1/day (RSE 28.5%)
    lkss <- log(7.59); label("Quasi-steady-state constant Kss = (Koff + Kint)/Kon (nmol/L)") # Table 1, Kss 7.59 nmol/L (RSE 35.7%)
    lkint <- log(3.23); label("First-order ginisortamab-gremlin-1 complex internalisation rate constant Kint (1/day)") # Table 1, Kint 3.23 1/day (RSE 31.3%)

    # ---- Inter-individual variability ----
    # Table 1 reports CV% = 100 * sqrt(exp(OMEGA) - 1), so OMEGA = log(1 + CV^2).
    # Vc and Vp share a BLOCK(2) (control stream $OMEGA BLOCK(2)) with correlation 0.97:
    # cov = 0.97 * sqrt(0.016504 * 0.047265) = 0.027092.
    etalcl ~ 0.051112 # Table 1, IIV CL 22.9% CV
    etalvc + etalvp ~ c(
      0.016504,
      0.027092, 0.047265
    ) # Table 1, IIV Vc 12.9% CV, IIV Vp 22.0% CV, correlation Vc-Vp 0.97
    etalkdeg ~ 0.13555 # Table 1, IIV Kdeg 38.1% CV
    etalkint ~ 0.26767 # Table 1, IIV Kint 55.4% CV

    # ---- Residual error (proportional; Table 1 note: CV% = 100 * THETA) ----
    propSd <- 0.108; label("Proportional residual error for free ginisortamab (fraction)") # Table 1, RUV ginisortamab 10.8%
    propSd_Ctotal_target <- 0.356; label("Proportional residual error for total gremlin-1 (fraction)") # Table 1, RUV total gremlin-1 35.6%
  })

  model({
    # Molecular weight used by the authors' control stream ($PK: MW = 150000) both
    # to convert free-drug nmol/L to ng/mL and to convert total gremlin-1 nmol/L to
    # ng/mL. The gremlin-1 parameters (Ksyn, Kss, baseline) were therefore fitted
    # on this 150 kDa scale and must be read back with the same factor.
    mw <- 150000 # g/mol (Supporting Information $PK)
    nmol_per_mg <- 1e6 / mw # 6.667 nmol per mg

    # ---- Individual parameters (Figure 2 equations; Supporting Information $PK) ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (BSA / 2)^e_bsa_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (BSA / 2)^e_bsa_vp
    ksyn <- exp(lksyn)
    kdeg <- exp(lkdeg + etalkdeg)
    kss <- exp(lkss)
    kint <- exp(lkint + etalkint)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Baseline free gremlin-1 (control stream BLR = KSYN / KDEG), nmol/L
    rbase <- ksyn / kdeg

    # ---- QSS algebra (Supporting Information $DES) ----
    # central holds the TOTAL ginisortamab amount (free + bound) in mg; ctot is
    # its concentration in nmol/L. total_target is the total gremlin-1
    # concentration in nmol/L. cfree is the control stream's CONC.
    ctot <- central * nmol_per_mg / vc
    cfree <- 0.5 * (ctot - total_target - kss) +
      0.5 * sqrt((ctot - total_target - kss)^2 + 4 * kss * ctot)

    # Control stream $DES, with DADT(1) and DADT(2) divided by nmol_per_mg so the
    # drug states carry mg:
    #   DADT(1) = -(K10+K12)*CONC*VC + K21*A(2) - KINT*A(3)*CONC*VC/(KSS+CONC)
    #   DADT(2) = K12*CONC*VC - K21*A(2)
    #   DADT(3) = KSYN - KDEG*A(3) - (KINT-KDEG)*CONC*A(3)/(KSS+CONC)
    d/dt(central) <- -(kel + k12) * cfree * vc / nmol_per_mg + k21 * peripheral1 -
      kint * total_target * cfree * vc / (kss + cfree) / nmol_per_mg
    d/dt(peripheral1) <- k12 * cfree * vc / nmol_per_mg - k21 * peripheral1
    d/dt(total_target) <- ksyn - kdeg * total_target -
      (kint - kdeg) * cfree * total_target / (kss + cfree)
    total_target(0) <- rbase

    # ---- Outputs ----
    # Free ginisortamab (the bridging assay measures free drug): nmol/L -> ug/mL
    Cc <- cfree / nmol_per_mg
    # Total gremlin-1, nmol/L -> ng/mL with the control stream's factor MW/1000
    Ctotal_target <- total_target * mw / 1000
    # Free gremlin-1 (control stream TFREE) and percent suppression relative to
    # the model baseline (Methods: change from the baseline free gremlin-1, BLR)
    tfree <- total_target * kss / (kss + cfree)
    target_suppression <- 100 * (1 - tfree / rbase)

    Cc ~ prop(propSd)
    Ctotal_target ~ prop(propSd_Ctotal_target)
  })
}
