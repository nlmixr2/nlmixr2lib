Plachouras_2009_colistin <- function() {
  description <- "Two-compartment population PK model for colistin methanesulfonate (CMS, prodrug) and one-compartment model for the formed colistin (active metabolite) in critically ill adults receiving 3 MU q8h IV CMS for multidrug-resistant Gram-negative infections (Plachouras 2009). Colistin metabolite parameters are apparent values scaled by the unknown fraction (fm) of CMS that forms colistin (CL_col is CL/fm; Vc_col is V/fm)."
  reference <- "Plachouras D, Karvanen M, Friberg LE, Papadomichelakis E, Antoniadou A, Tsangaris I, Karaiskos I, Poulakou G, Kontopidou F, Armaganidis A, Cars O, Giamarellou H. Population pharmacokinetic analysis of colistin methanesulfonate and colistin after intravenous administration in critically ill patients with infections caused by gram-negative bacteria. Antimicrob Agents Chemother. 2009;53(8):3430-3436. doi:10.1128/AAC.01361-08."
  vignette <- "Plachouras_2009_colistin"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 1L,
    age_range      = "40-83 years",
    age_mean       = "63.6 years",
    weight_range   = "65-110 kg (actual); 60-85 kg (ideal)",
    sex_female_pct = 33.3,
    race_ethnicity = "Not reported; single-center Greek cohort.",
    disease_state  = "Critically ill adult ICU patients with probable or documented infections caused by multidrug-resistant Gram-negative bacteria (Acinetobacter baumannii, Pseudomonas aeruginosa, Enterobacteriaceae). Patients on continuous venovenous hemodiafiltration were excluded.",
    dose_range     = "CMS 3 MU (~240 mg) q8h as 15-min IV infusion. Two patients with CrCl < 50 mL/min received the empirically reduced regimen of 2 MU q8h (~160 mg q8h).",
    crcl_range     = "41-126 mL/min (Cockcroft-Gault); mean +- SD = 82.3 +- 24.35 mL/min on day 1.",
    apache_ii      = "median 13 (range 5-20).",
    regions        = "Single center: Attikon University General Hospital, Athens, Greece.",
    notes          = "Demographics from Table 1 of Plachouras 2009 (n = 18 patients; 12 males, 6 females). The fraction of CMS metabolised to colistin (fm) is not identifiable from CMS / colistin plasma data alone, so colistin CL and V are reported as the apparent values CL/fm and V/fm; the model file carries them as cl_col and vc_col without the explicit fm rescaling."
  )

  ini({
    # =====================================================================
    # Final population estimates from Plachouras 2009 Table 2 (page 3432).
    # Data fit in NONMEM VI by FOCE-INTER on log-transformed concentrations
    # in molar units (paper p. 3431 - Methods, "Population pharmacokinetic
    # modeling"). Molar masses used by the paper: 1743 g/mol CMS,
    # 1163 g/mol colistin. Additive residual SDs are converted from the
    # paper's nmol/L units to mg/L for the model below; the conversions
    # are spelled out next to each addSd line.
    # =====================================================================

    # ---- CMS structural PK (2-compartment, linear elimination) ----
    lcl  <- log(13.7);  label("CMS clearance CL (L/h)")                 # Table 2, "CMS" "CL (liters/h)" typical value = 13.7
    lvc  <- log(13.5);  label("CMS central volume of distribution V1 (L)") # Table 2, "CMS" "V1 (liters)" typical value = 13.5
    lq   <- log(133);   label("CMS intercompartmental clearance Q (L/h)")  # Table 2, "CMS" "Q (liters/h)" typical value = 133
    lvp  <- log(28.9);  label("CMS peripheral volume of distribution V2 (L)") # Table 2, "CMS" "V2 (liters)" typical value = 28.9

    # ---- Colistin structural PK (1-compartment; apparent values scaled by fm) ----
    lcl_col <- log(9.09); label("Apparent colistin clearance CL/fm (L/h); fm = unknown fraction of CMS metabolised to colistin") # Table 2, "Colistin" "CL/fm (liters/h)" typical value = 9.09
    lvc_col <- log(189);  label("Apparent colistin volume of distribution V/fm (L)")                                              # Table 2, "Colistin" "V/fm (liters)" typical value = 189

    # ---- IIV (omega^2 on the log scale; omega^2 = log(1 + CV^2)) ----
    # The paper reports IIV as CV%; the in-line `# CV%` comment converts
    # each value to the log-scale variance via omega^2 = log(1 + CV^2).
    # IIV on V1, Q, V2 (CMS) and on the colistin V/fm was not reported in
    # Table 2 (only CL CMS, CL/fm colistin, and colistin residual error
    # had IIV).
    etalcl     ~ 0.1289; label("IIV variance on log CMS CL (37% CV per Table 2)")        # log(1 + 0.37^2) = 0.1289
    etalcl_col ~ 0.2964; label("IIV variance on log apparent colistin CL/fm (59% CV per Table 2)") # log(1 + 0.59^2) = 0.2964

    # ---- Residual error ----
    # Combined additive + proportional error model on log-transformed
    # molar concentrations (paper p. 3431). Proportional values come
    # directly from Table 2 (% CV / 100); additive values are converted
    # from nmol/L to mg/L via the paper's molar masses.
    # CMS:        9.11 nmol/L * 1743 g/mol = 1.588e-5 g/L = 0.01588 mg/L
    # Colistin:   4.98 nmol/L * 1163 g/mol = 5.79e-6  g/L = 0.00579 mg/L
    propSd     <- 0.220;   label("CMS proportional residual SD (fraction)")        # Table 2, "CMS"      "Residual error" "Proportional (%)" = 22.0
    addSd      <- 0.01588; label("CMS additive residual SD (mg/L)")                # Table 2, "CMS"      "Residual error" "Additive (nmol/liter)" = 9.11; 9.11e-9 mol/L * 1743 g/mol
    propSd_col <- 0.0719;  label("Colistin proportional residual SD (fraction)")   # Table 2, "Colistin" "Residual error" "Proportional (%)" = 7.19
    addSd_col  <- 0.00579; label("Colistin additive residual SD (mg/L)")           # Table 2, "Colistin" "Residual error" "Additive (nmol/liters)" = 4.98; 4.98e-9 mol/L * 1163 g/mol
  })

  model({
    # ---- Molar-mass conversion CMS -> colistin ----
    # The paper fit on molar concentrations (paper p. 3431) so the CMS
    # clearance flux maps mole-for-mole to the colistin formation flux.
    # When the model runs in mass units (mg dose, mg/L plasma), 1 mg of
    # cleared CMS forms (mw_col / mw_cms) mg of colistin (= fm * that, but
    # fm is folded into the apparent colistin V/CL). Molar masses come
    # from paper p. 3431 ("molar masses are, on average, 1,743 g/mol for
    # CMS and 1,163 g/mol for colistin").
    mw_cms <- 1743                     # CMS g/mol (paper p. 3431)
    mw_col <- 1163                     # Colistin g/mol (paper p. 3431)
    mass_col_per_cms <- mw_col / mw_cms # = 0.6672

    # ---- Individual structural parameters ----
    cl     <- exp(lcl + etalcl)           # CMS CL (L/h); rate-determining for CMS elimination AND for the colistin formation flux
    vc     <- exp(lvc)                    # CMS V1 (L); no IIV reported
    q      <- exp(lq)                     # CMS Q (L/h); no IIV reported
    vp     <- exp(lvp)                    # CMS V2 (L); no IIV reported
    cl_col <- exp(lcl_col + etalcl_col)   # Colistin apparent CL/fm (L/h)
    vc_col <- exp(lvc_col)                # Colistin apparent V/fm (L); no IIV reported

    # ---- Micro-constants ----
    kel     <- cl  / vc        # CMS first-order elimination (1/h); the fraction fm of this clearance forms colistin
    k12     <- q   / vc        # CMS central -> CMS peripheral (1/h)
    k21     <- q   / vp        # CMS peripheral -> CMS central (1/h)
    kel_col <- cl_col / vc_col # Colistin first-order elimination (1/h); equivalent to CL/V since both apparent values share the same fm factor

    # ---- ODE system ----
    # CMS disposition is a standard 2-compartment IV system. Colistin is
    # formed by a first-order process from CMS at the rate fm * CL * Ccms
    # (mole basis). Because fm is non-identifiable from these data, the
    # colistin compartment carries the scaled mass
    #   central_col_modelled = A_col(t) / fm    (mg colistin / fm),
    # so the observed concentration is
    #   central_col / vc_col = (A_col/fm) / (V_col/fm) = A_col / V_col,
    # which is the true colistin concentration in mg/L.
    # The factor (mw_col / mw_cms) converts the CMS clearance flux
    # (mg CMS per hour) to colistin mass (mg colistin per hour), since the
    # paper's molar-unit fit is mass-blind to the species exchange:
    #   d A_col/dt = fm * (mw_col / mw_cms) * cl * C_cms - cl_col_true * C_col
    # Dividing by fm gives the equation below in central_col_modelled
    # (= A_col/fm) terms. Without this factor, the colistin compartment
    # would carry CMS-mass equivalents and Cc_col would over-predict by
    # mw_cms / mw_col = 1.499 (verified against the paper's typical Cmax
    # values: first-dose ~ 0.60 mg/L and steady-state ~ 2.3 mg/L per
    # paper p. 3433 / Fig. 4).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(central_col) <-  mass_col_per_cms * cl * central / vc - kel_col * central_col

    # ---- Observations and residual error ----
    # CMS concentration in mg CMS / L; colistin concentration in mg colistin / L.
    # The paper fits on log-transformed concentrations in molar units;
    # the typical-value prediction is unchanged on the linear (mg/L)
    # scale used here, and the combined additive + proportional residual
    # is the standard nlmixr2 mapping.
    Cc     <- central     / vc
    Cc_col <- central_col / vc_col
    Cc     ~ add(addSd)     + prop(propSd)
    Cc_col ~ add(addSd_col) + prop(propSd_col)
  })
}
