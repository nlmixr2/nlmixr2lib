Kim_2019_sildenafil_dog_ivivc <- function() {
  description <- paste(
    "Preclinical (beagle dog). Level-A in vitro-in vivo correlation (IVIVC)",
    "model for oral sildenafil tablets with site-dependent dissolution (Kim",
    "2019). Identical to Kim_2019_sildenafil_dog (dose-normalised",
    "Michaelis-Menten tablet dissolution scaled by a time-varying Hill-type",
    "dissolved fraction, lumen -> gut -> two-compartment disposition, Table 4",
    "parameters) except that the in vivo maximum dissolution rate is not a",
    "fitted per-formulation value: it is computed from the in vitro maximum",
    "dissolution rate at pH 1.2 (Table 3) through the paper's power",
    "correlation, Vmax,in vivo/Dose = 2.1411 * (Vmax,in vitro/Dose)^0.4388 -",
    "0.6346 (Equation 13). This is the model the authors used to predict",
    "plasma profiles of the IR (20 mg) and three SR (60 mg) tablets from their",
    "in vitro dissolution (Figure 6, Table 5). To predict a new tablet,",
    "override the in vitro Vmax/Dose of the reference (IR) arm with",
    "ini(lvmax_invitro = log(<value>)); the correlation was built on in vitro",
    "Vmax/Dose of 0.124-6.58 1/h and returns a non-positive in vivo rate below",
    "about 0.062 1/h. As in the companion model, the printed Table 4 CL and",
    "CLD rows are transposed (CL = 9.97 L/h, CLD = 49.2 L/h are used) and the",
    "residual error is not reported (zero placeholder)."
  )
  reference <- paste(
    "Kim TH, Shin S, Jeong SW, Lee JB, Shin BS.",
    "Physiologically Relevant In Vitro-In Vivo Correlation (IVIVC) Approach",
    "for Sildenafil with Site-Dependent Dissolution.",
    "Pharmaceutics. 2019;11(6):251. doi:10.3390/pharmaceutics11060251"
  )
  vignette <- "Kim_2019_sildenafil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    FORM_SILDENAFIL_SRFAST = list(
      description = "1 = SRfast sustained-release tablet (10% HPMC 2208), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IR tablet when all three FORM_SILDENAFIL_SR* indicators are 0)",
      notes = paste(
        "Selects the SRfast in vitro maximum dissolution rate (Table 3), and",
        "with it the correlated in vivo rate, and its in vivo between-subject",
        "variability (Table 4). Mutually exclusive with the other two",
        "FORM_SILDENAFIL_SR* indicators."
      ),
      source_name = "SRfast"
    ),
    FORM_SILDENAFIL_SRMEDIUM = list(
      description = "1 = SRmedium sustained-release tablet (30% HPMC 2208), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IR tablet when all three FORM_SILDENAFIL_SR* indicators are 0)",
      notes = "Selects the SRmedium in vitro maximum dissolution rate (Table 3), hence the correlated in vivo rate, and its in vivo between-subject variability (Table 4).",
      source_name = "SRmedium"
    ),
    FORM_SILDENAFIL_SRSLOW = list(
      description = "1 = SRslow sustained-release tablet (60% HPMC 2208), 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IR tablet when all three FORM_SILDENAFIL_SR* indicators are 0)",
      notes = "Selects the SRslow in vitro maximum dissolution rate (Table 3), hence the correlated in vivo rate, and its in vivo between-subject variability (Table 4).",
      source_name = "SRslow"
    )
  )

  dosing <- "depot"

  # The in vivo maximum dissolution rate is not an ini() parameter here: it
  # is the Equation (13) transform of the in vitro rate. These etas carry
  # the Table 4 BSV of that in vivo rate, one per formulation.
  paper_specific_etas <- c("etalvmax", "etalvmax_srfast", "etalvmax_srmedium", "etalvmax_srslow")

  compartmentData <- list(
    depot = list(
      analyte = "sildenafil (undissolved, in the tablet)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    gut_lumen = list(
      analyte = "sildenafil (dissolved, gastrointestinal lumen)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    gut = list(
      analyte = "sildenafil (dissolved, gut compartment)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "sildenafil", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sildenafil", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "beagle dog",
    n_subjects = 12L,
    n_studies = 1L,
    age_range = "14-18 months",
    weight_range = "not reported",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range = paste(
      "Single oral dose after overnight fasting, with 30 mL water: IR tablet 20 mg",
      "or one of three SR tablets 60 mg (SRfast, SRmedium, SRslow); n = 3 dogs per",
      "formulation (parallel groups)."
    ),
    regions = "Korea (KNOTUS Co., Ltd.)",
    notes = paste(
      "Methods 2.4: 12 male beagle dogs (Orient Bio). Plasma sampled at 0, 0.25,",
      "0.5, 0.75, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 12, 24, 36 and 48 h.",
      "LC-MS/MS assay, LLOQ 5 ng/mL (Methods 2.5.2). Body weight is not reported",
      "and is not used by the model (doses are flat mg tablets)."
    )
  )

  ini({
    # ---- Systemic disposition (Kim 2019 Table 4) -------------------------
    # Table 4 prints 'Systemic clearance CL 49.2 (0.379)' and 'Distribution
    # clearance CLD 9.97 (0.2)'. These two rows are transposed: the IR tablet
    # dissolves completely (Figure 4B), so AUC(0-inf) = Dose/CL exactly; with
    # CL = 49.2 L/h a 20 mg dose gives 406 ng*h/mL against the observed
    # 2007 (Table 2 AUCinf) and the paper's own IVIVC prediction 1942 (Table
    # 5), whereas CL = 9.97 L/h gives 2006. Only CLD = 49.2 L/h then also
    # reproduces the Table 5 predicted Cmax (IR 351 vs 338.7; CLD = 9.97
    # gives 502). The whole rows (value and BSV) are swapped here.
    lcl <- log(9.97); label("Systemic clearance CL (L/h)") # Table 4 row printed as 'Distribution clearance CLD' = 9.97 L/h (rows transposed; see block comment)
    lq <- log(49.2); label("Distribution clearance CLD (L/h)") # Table 4 row printed as 'Systemic clearance CL' = 49.2 L/h (rows transposed; see block comment)
    lvc <- log(22); label("Central volume of distribution V1 (L)") # Table 4 V1 = 22 L
    lvp <- log(31); label("Peripheral volume of distribution V2 (L)") # Table 4 V2 = 31 L

    # ---- Post-dissolution transfer and absorption (Table 4) ---------------
    lklag <- log(3.13); label("Transfer rate constant of dissolved drug from lumen to gut klag (1/h)") # Table 4 klag = 3.13 1/h
    lka <- log(4.37); label("Absorption rate constant from gut to central ka (1/h)") # Table 4 ka = 4.37 1/h

    # ---- Time-varying dissolved fraction FDiss,total (Eqs 3-7, Table 4) --
    ltget <- log(0.73); label("Time of half-maximal decrease of FDiss,stomach (gastric emptying) TGET (h)") # Table 4 TGET = 0.73 h
    limax <- log(0.992); label("Maximum fractional decrease of FDiss,stomach Imax (fraction)") # Table 4 Imax = 0.992
    lhill_stomach <- log(16.7); label("Hill coefficient of the FDiss,stomach decrease (unitless)") # Table 4 Hillstomach = 16.7
    ltitt <- log(2.18); label("Intestinal transit time TITT (h)") # Table 4 TITT = 2.18 h
    ltctt <- log(4.1); label("Colon transit time TCTT (h)") # Table 4 TCTT = 4.1 h
    ldissmax <- log(0.115); label("Maximum FDiss,intestine Dissmax (fraction)") # Table 4 Dissmax = 0.115
    lhill_intestine <- log(13.1); label("Hill coefficient of the FDiss,intestine rise and fall (unitless)") # Table 4 Hillintestine = 13.1

    # ---- In vivo Michaelis-Menten dissolution (Eq 2, Table 4) -------------
    # Dose-normalised (Methods 2.6.3): AM50/dose is a fraction of the dose.
    lam50 <- log(0.315); label("Dose-normalised amount in the tablet at half-maximal in vivo dissolution rate AM50,in vivo/dose (fraction)") # Table 4 AM50,in vivo/dose = 0.315

    # ---- IVIVC inputs: in vitro Vmax/Dose at pH 1.2 (Table 3) -------------
    # Estimated in the separate in vitro dissolution fit (Kim_2019_sildenafil_invitro)
    # and carried into the IVIVC as fixed inputs.
    lvmax_invitro <- fixed(log(6.58)); label("Dose-normalised in vitro maximum dissolution rate at pH 1.2, IR tablet (1/h)") # Table 3 Vmax,in vitro,IR/dose = 6.58 1/h
    lvmax_invitro_srfast <- fixed(log(1.40)); label("Dose-normalised in vitro maximum dissolution rate at pH 1.2, SRfast tablet (1/h)") # Table 3 Vmax,in vitro,SRfast/dose = 1.40 1/h
    lvmax_invitro_srmedium <- fixed(log(0.358)); label("Dose-normalised in vitro maximum dissolution rate at pH 1.2, SRmedium tablet (1/h)") # Table 3 Vmax,in vitro,SRmedium/dose = 0.358 1/h
    lvmax_invitro_srslow <- fixed(log(0.124)); label("Dose-normalised in vitro maximum dissolution rate at pH 1.2, SRslow tablet (1/h)") # Table 3 Vmax,in vitro,SRslow/dose = 0.124 1/h

    # ---- IVIVC power correlation (Equation 13) ----------------------------
    # Fitted by regression (SigmaPlot) across the four tablets and carried
    # into the IVIVC model as fixed constants.
    ivivc_coef <- fixed(2.1411); label("IVIVC power-correlation coefficient (1/h per (1/h)^exponent)") # Equation (13) coefficient 2.1411
    ivivc_exp <- fixed(0.4388); label("IVIVC power-correlation exponent (unitless)") # Equation (13) exponent 0.4388
    ivivc_int <- fixed(-0.6346); label("IVIVC power-correlation intercept (1/h)") # Equation (13) intercept -0.6346

    # ---- Between-subject variability (Table 4 BSV column) -----------------
    # S-ADAPT exponential model (Methods 2.6.3). The BSV column is the SD of
    # the log-scale random effect (S-ADAPT 'apparent CV'), so omega^2 =
    # BSV^2; the 0.0316 entries are sqrt(0.001), the usual small-variance
    # floor of S-ADAPT fits.
    etalcl ~ 0.04 # Table 4 BSV 0.2 (printed on the CLD row; rows transposed) -> 0.2^2
    etalq ~ 0.143641 # Table 4 BSV 0.379 (printed on the CL row; rows transposed) -> 0.379^2
    etalvc ~ 0.024649 # Table 4 V1 BSV 0.157 -> 0.157^2
    etalvp ~ 0.564001 # Table 4 V2 BSV 0.751 -> 0.751^2
    etalklag ~ 0.247009 # Table 4 klag BSV 0.497 -> 0.497^2
    etalka ~ 0.906304 # Table 4 ka BSV 0.952 -> 0.952^2
    etaltget ~ 0.148225 # Table 4 TGET BSV 0.385 -> 0.385^2
    etalimax ~ 0.00000729 # Table 4 Imax BSV 0.0027 -> 0.0027^2
    etalhill_stomach ~ 0.00099856 # Table 4 Hillstomach BSV 0.0316 -> 0.0316^2
    etaltitt ~ 0.188356 # Table 4 TITT BSV 0.434 -> 0.434^2
    etaltctt ~ 0.0121 # Table 4 TCTT BSV 0.11 -> 0.11^2
    etaldissmax ~ 0.128881 # Table 4 Dissmax BSV 0.359 -> 0.359^2
    etalhill_intestine ~ 0.00099856 # Table 4 Hillintestine BSV 0.0316 -> 0.0316^2
    etalam50 ~ 0.011025 # Table 4 AM50,in vivo BSV 0.105 -> 0.105^2
    # The correlated in vivo Vmax carries the per-formulation Table 4 BSV of
    # the fitted in vivo Vmax; the paper does not state which variability its
    # IVIVC Monte Carlo simulations (Section 2.7) applied.
    etalvmax ~ 0.173056 # Table 4 Vmax,in vivo,IR BSV 0.416 -> 0.416^2
    etalvmax_srfast ~ 0.021025 # Table 4 Vmax,in vivo,SRfast BSV 0.145 -> 0.145^2
    etalvmax_srmedium ~ 0.087025 # Table 4 Vmax,in vivo,SRmedium BSV 0.295 -> 0.295^2
    etalvmax_srslow ~ 0.065025 # Table 4 Vmax,in vivo,SRslow BSV 0.255 -> 0.255^2

    # ---- Residual error -------------------------------------------------
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)") # not reported: Table 4 and Methods 2.6.3 give no residual-error model or estimate
  })

  model({
    # Individual disposition and absorption parameters
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    klag <- exp(lklag + etalklag)
    ka <- exp(lka + etalka)

    # Dissolved-fraction parameters
    tget <- exp(ltget + etaltget)
    imax <- exp(limax + etalimax)
    hill_stomach <- exp(lhill_stomach + etalhill_stomach)
    titt <- exp(ltitt + etaltitt)
    tctt <- exp(ltctt + etaltctt)
    dissmax <- exp(ldissmax + etaldissmax)
    hill_intestine <- exp(lhill_intestine + etalhill_intestine)

    # In vivo dissolution parameters; IR is the reference formulation
    am50 <- exp(lam50 + etalam50)
    # Equation (13): in vivo Vmax/Dose from in vitro Vmax/Dose at pH 1.2
    vmax_ir <- (ivivc_coef * exp(lvmax_invitro)^ivivc_exp + ivivc_int) * exp(etalvmax)
    vmax_srfast <- (ivivc_coef * exp(lvmax_invitro_srfast)^ivivc_exp + ivivc_int) * exp(etalvmax_srfast)
    vmax_srmedium <- (ivivc_coef * exp(lvmax_invitro_srmedium)^ivivc_exp + ivivc_int) * exp(etalvmax_srmedium)
    vmax_srslow <- (ivivc_coef * exp(lvmax_invitro_srslow)^ivivc_exp + ivivc_int) * exp(etalvmax_srslow)
    form_ir <- 1 - FORM_SILDENAFIL_SRFAST - FORM_SILDENAFIL_SRMEDIUM - FORM_SILDENAFIL_SRSLOW
    vmax <- form_ir * vmax_ir +
      FORM_SILDENAFIL_SRFAST * vmax_srfast +
      FORM_SILDENAFIL_SRMEDIUM * vmax_srmedium +
      FORM_SILDENAFIL_SRSLOW * vmax_srslow

    # Time-varying dissolved fraction, Equations (3)-(7). 'Time' is the time
    # since the tablet was given.
    tdiss <- tad(depot)
    if (is.na(tdiss)) tdiss <- 0
    fdiss_stomach <- 1 - imax * tdiss^hill_stomach / (tget^hill_stomach + tdiss^hill_stomach)
    ttrans1 <- tget + titt
    ttrans2 <- tget + titt + tctt
    fdiss_intestine <- dissmax *
      (tdiss^hill_intestine / (ttrans1^hill_intestine + tdiss^hill_intestine)) *
      (1 - tdiss^hill_intestine / (ttrans2^hill_intestine + tdiss^hill_intestine))
    fdiss_total <- fdiss_stomach + fdiss_intestine

    # Dose-normalised Michaelis-Menten dissolution, Equations (2) and (8):
    # d(X/Dose)/dt = -Vmax/Dose * (X/Dose) / (AM50/Dose + X/Dose) * FDiss,total
    # written in amounts. podo(depot) is the tablet dose.
    dose_tab <- podo(depot)
    rdiss <- 0
    if (dose_tab > 0) rdiss <- vmax * depot / (am50 + depot / dose_tab) * fdiss_total

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -rdiss
    d/dt(gut_lumen) <- rdiss - klag * gut_lumen
    d/dt(gut) <- klag * gut_lumen - ka * gut
    d/dt(central) <- ka * gut - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL; x 1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
