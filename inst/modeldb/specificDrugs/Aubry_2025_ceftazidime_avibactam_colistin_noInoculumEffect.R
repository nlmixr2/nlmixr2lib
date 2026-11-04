Aubry_2025_ceftazidime_avibactam_colistin_noInoculumEffect <- function() {
  description <- "In vitro (Klebsiella pneumoniae NARA864, KPC-3-producing clinical isolate). Semi-mechanistic PK/PD model of static time-kill experiments with ceftazidime (in the presence of a fixed 4 mg/L avibactam) and colistin, alone and in combination, fitted to the standard 5e5 CFU/mL inoculum only (the 'no-inoculum-effect model'). A single homogeneous bacterial population with logistic growth is killed by each drug through a sigmoidal Emax term whose Emax is progressively reduced by an empirical two-state adaptive-resistance model, combined under Bliss independence, with a general pharmacodynamic interaction (GPDI) term through which colistin lowers the ceftazidime/avibactam Emax. Unlike the companion inoculum-effect model it has no dependence of the kill parameters on the starting inoculum, so it serves as the comparator that quantifies the inoculum effect. Drug concentrations are static covariates (no dosing events). Companion inoculum-effect model fitted to all four inocula: Aubry_2025_ceftazidime_avibactam_colistin_inoculumEffect."
  reference <- "Aubry R, Buyck JM, Chauzy A, Prouvensier L, Decousser JW, Nordmann P, Wicha SG, Marchand S, Gregoire N. PKPD modeling of the inoculum effect of combined ceftazidime/avibactam and colistin against KPC-3 Klebsiella pneumoniae isolate. Antimicrob Agents Chemother. 2025;69(5):e01797-24. doi:10.1128/aac.01797-24. Structural equations: the no-inoculum-effect model reuses the structure of Aubry et al. 2023 (doi:10.1128/aac.00301-23), re-estimated on the present study's 5e5 CFU/mL data (Supplemental Text S1). Parameter estimates: Supplemental Table S1."
  vignette <- "Aubry_2025_ceftazidime_avibactam_colistin"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Adaptive-resistance states for ceftazidime/avibactam. `aron` / `aroff` are
  # canonical and accept a registered drug suffix (as `aron_col` does), but
  # `caz` is not a registered suffix token.
  paper_specific_compartments <- c("aron_caz", "aroff_caz")

  compartmentData <- list(
    cfu = list(
      analyte = "Klebsiella pneumoniae NARA864 (total viable bacteria)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron_caz = list(
      analyte = "fraction of bacteria adapted to ceftazidime/avibactam",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff_caz = list(
      analyte = "fraction of bacteria not adapted to ceftazidime/avibactam",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron_col = list(
      analyte = "fraction of bacteria adapted to colistin",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff_col = list(
      analyte = "fraction of bacteria not adapted to colistin",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONC_CAZ_MGL = list(
      description = "Static (time-invariant) ceftazidime concentration in the time-kill tube, always in the presence of a fixed 4 mg/L avibactam. Drives the sigmoidal Emax kill by ceftazidime/avibactam and switches on adaptive resistance to it.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the time-kill grid spanned 0.5-8 mg/L (Methods; MIC 2 mg/L)",
      notes = "Held constant for the 30-h experiment. The paper writes 'CZA' for the ceftazidime concentration and states that it always corresponds to ceftazidime with 4 mg/L avibactam; avibactam is assumed to inhibit the beta-lactamases fully without a bactericidal effect of its own, so it carries no covariate. Adaptive resistance develops only while the concentration is above zero (Code S1 $DES: IF(A(1).GT.0) ADAPT_CZA=KON_CZA). This no-inoculum-effect model was fitted to the 5e5 CFU/mL experiments only.",
      source_name = "CAZ_CONC"
    ),
    CONC_CST_MGL = list(
      description = "Static (time-invariant) colistin concentration in the time-kill tube. Drives the sigmoidal Emax kill by colistin, switches on adaptive resistance to it, and as the GPDI perpetrator lowers the ceftazidime/avibactam Emax.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the time-kill grid spanned 0.06-8 mg/L (Methods; MIC 0.25 mg/L)",
      notes = "Colistin (Sigma Aldrich), concentrations expressed as colistin. Held constant for the 30-h experiment. Adaptive resistance develops only while the concentration is above zero (Code S1 $DES: IF(A(2).GT.0) ADAPT_CST=KON_CST).",
      source_name = "CST_CONC"
    )
  )

  population <- list(
    species = "in vitro (Klebsiella pneumoniae NARA864 clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Clinical KPC-3-producing K. pneumoniae isolate NARA864, obtained from the National Reference Center for Emerging Antibiotic Resistance (Fribourg, Switzerland). Susceptible to ceftazidime/avibactam (MIC 2 mg/L with avibactam fixed at 4 mg/L) and to colistin (MIC 0.25 mg/L). The same isolate was used in the earlier interaction study (Aubry et al. 2023).",
    model_system = "Static time-kill experiments over 30 h at 35 degrees C with shaking (130 rpm), 10 mL cation-adjusted Mueller-Hinton broth II, at the standard 5e5 CFU/mL inoculum. Ceftazidime 0.5-8 mg/L (with avibactam 4 mg/L) or colistin 0.06-8 mg/L alone, and combinations from 0.5+0.06 to 4+0.5 mg/L. Bacterial counts at 0, 2, 4, 8, 24 and 30 h by spiral plating; lower limit of quantification 200 CFU/mL (2.3 log10 CFU/mL).",
    initial_inoculum = "Estimated 5.54 log10 CFU/mL (Table S1), nominal 5e5 CFU/mL",
    dose_range = "Static concentrations: ceftazidime 0.5-8 mg/L (with avibactam 4 mg/L), colistin 0.06-8 mg/L, alone and in combination (Methods)",
    notes = paste(
      "Estimated in NONMEM 7.4 with the Laplacian method (ADVAN13, TOL=9); data below the limit of quantification were handled with Beal's M3 method.",
      "The random effects are an estimated inter-experimental variability on the inoculum (etainoc, 9.40% CV) and an additive residual error on the log10 bacterial count.",
      "The no-inoculum-effect model reuses the structure of Aubry et al. 2023, re-estimated on this study's 5e5 CFU/mL data; it is the comparator against which the companion inoculum-effect model quantifies the inoculum effect.",
      "The colistin Emax and the interaction EC50 and Hill coefficient were not estimated (Table S1 'fixed')."
    )
  )

  ini({
    # --- Inoculum and bacterial growth (Table S1) ---
    log10inoc <- fixed(5.54); label("Estimated inoculum size at time 0 (log10 CFU/mL)")             # Table S1: INOC 5.54 fixed
    kg        <- 1.76; label("Apparent bacterial growth rate constant (1/h)")                       # Table S1: K_G 1.76 (RSE 8%)
    log10bmax <- 9.03; label("Maximal bacterial population supported by the system (log10 CFU/mL)") # Table S1: Bmax 9.03 (RSE 1%)

    # --- Ceftazidime/avibactam effect and adaptive resistance (Table S1) ---
    emax_caz  <- 7.53;  label("Maximum kill rate constant of ceftazidime/avibactam (1/h)")                     # Table S1: Emax_CZA 7.53 (RSE 3%)
    ec50_caz  <- 0.46;  label("Ceftazidime concentration giving half of the maximum kill rate (mg/L)")        # Table S1: EC50_CZA 0.46 (RSE 7%)
    hill_caz  <- 1.66;  label("Sigmoidicity of the ceftazidime/avibactam effect (unitless)")                   # Table S1: gamma_CZA 1.66 (RSE 10%)
    kon_caz   <- 0.49;  label("Rate constant for development of adaptive resistance to ceftazidime/avibactam (1/h)") # Table S1: K_ON_CZA 0.49 (RSE 4%)
    armax_caz <- 0.800; label("Maximal reduction of the ceftazidime/avibactam Emax when all bacteria are adapted (fraction)") # Table S1: ARmax_CZA 80.0 percent (RSE 3%), as a fraction

    # --- Reversal of adaptive resistance (Table S1), shared by both drugs ---
    koff <- fixed(0); label("Rate constant for reversal of adaptive resistance (1/h)") # Table S1: K_OFF 0 fixed

    # --- Colistin effect and adaptive resistance (Table S1) ---
    emax_col  <- fixed(8.94); label("Maximum kill rate constant of colistin (1/h)")                 # Table S1: Emax_CST 8.94 fixed
    ec50_col  <- 0.23;  label("Colistin concentration giving half of the maximum kill rate (mg/L)") # Table S1: EC50_CST 0.23 (RSE 13%)
    hill_col  <- 1.11;  label("Sigmoidicity of the colistin effect (unitless)")                     # Table S1: gamma_CST 1.11 (RSE 5%)
    kon_col   <- 0.77;  label("Rate constant for development of adaptive resistance to colistin (1/h)") # Table S1: K_ON_CST 0.77 (RSE 5%)
    armax_col <- 0.810; label("Maximal reduction of the colistin Emax when all bacteria are adapted (fraction)") # Table S1: ARmax_CST 81.0 percent (RSE 2%), as a fraction

    # --- GPDI interaction: colistin (perpetrator) lowers the ceftazidime/avibactam Emax (Table S1) ---
    int_col_caz     <- -0.60;       label("Maximal fractional change of the ceftazidime/avibactam Emax caused by colistin (fraction)") # Table S1: INT_Emax -60.0 percent (RSE 23%), as a fraction
    ec50int_col_caz <- fixed(0.38); label("Colistin concentration giving half of the maximal Emax change (mg/L)")                      # Table S1: EC50_INT_Emax 0.38 fixed
    hillint_col_caz <- fixed(1);    label("Sigmoidicity of the colistin-on-ceftazidime/avibactam interaction (unitless)")              # Table S1: gamma_INT_Emax 1 fixed

    # --- Inter-experimental variability on the inoculum (Table S1) ---
    etainoc ~ 0.008797 # IEV_INOC 9.40% CV estimated (Table S1, RSE 33%); log(0.094^2 + 1) = 0.008797; multiplies the log10-scale inoculum

    # --- Residual error (Table S1; Code S1 $ERROR) ---
    # Table S1 prints RES_ADD = 2.06 'log10 CFU/mL', but the value is the NONMEM
    # $SIGMA variance (W = SQRT(SIGMA(1))), as for the companion model. SD = sqrt(2.06).
    addSd <- sqrt(2.06); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # Table S1: RES_ADD 2.06 (RSE 2%) read as a variance; SD = sqrt(2.06) = 1.44
  })

  model({
    bmax <- 10^log10bmax

    # Inoculum with inter-experimental variability (Code S1 $PK:
    # INOC = INOC_TV * EXP(ETA(1))).
    inoc_i <- log10inoc * exp(etainoc)

    # GPDI (Text S1): the victim parameter is the ceftazidime/avibactam Emax,
    # theta_INT = theta_mono * (1 + INT * C_perp^gINT / (EC50_INT^gINT + C_perp^gINT)).
    fint_col_caz <- int_col_caz * CONC_CST_MGL^hillint_col_caz /
      (ec50int_col_caz^hillint_col_caz + CONC_CST_MGL^hillint_col_caz)
    emax_caz_i <- emax_caz * (1 + fint_col_caz)

    # Adaptation develops only while the drug is present (Code S1 $DES:
    # ADAPT_CZA = 0; IF(A(1).GT.0) ADAPT_CZA = KON_CZA, and likewise for colistin).
    kon_caz_i <- kon_caz * (CONC_CAZ_MGL > 0)
    kon_col_i <- kon_col * (CONC_CST_MGL > 0)

    # Sigmoidal Emax kill, Emax reduced in proportion to the adapted fraction.
    kill_caz <- emax_caz_i * (1 - armax_caz * aron_caz) * CONC_CAZ_MGL^hill_caz /
      (ec50_caz^hill_caz + CONC_CAZ_MGL^hill_caz)
    kill_col <- emax_col * (1 - armax_col * aron_col) * CONC_CST_MGL^hill_col /
      (ec50_col^hill_col + CONC_CST_MGL^hill_col)

    # Bliss independence (Code S1 $DES): each kill rate is normalised by the
    # larger of the two typical (base) Emax values, combined as Ea + Eb - Ea*Eb,
    # and rescaled; algebraically kill_caz + kill_col - kill_caz * kill_col / emax_norm.
    emax_norm <- max(emax_caz, emax_col)
    kill_comb <- kill_caz + kill_col - kill_caz * kill_col / emax_norm

    # Logistic growth on the linear CFU scale (Code S1 DADT(3)).
    d/dt(cfu) <- (kg * (1 - cfu / bmax) - kill_comb) * cfu

    # Two-state adaptive resistance per drug (Code S1 DADT(4)-DADT(7)).
    d/dt(aron_caz)  <- kon_caz_i * aroff_caz - koff * aron_caz
    d/dt(aroff_caz) <- koff * aron_caz - kon_caz_i * aroff_caz
    d/dt(aron_col)  <- kon_col_i * aroff_col - koff * aron_col
    d/dt(aroff_col) <- koff * aron_col - kon_col_i * aroff_col

    # Initial conditions (Code S1 A_0): inoculum 10^INOC, all bacteria non-adapted.
    cfu(0) <- 10^inoc_i
    aron_caz(0) <- 0
    aroff_caz(0) <- 1
    aron_col(0) <- 0
    aroff_col(0) <- 1

    # Observation: log10 CFU/mL with the 1e-6 CFU/mL floor of Code S1 $ERROR
    # (BTOT_ERR). The fitted data were censored below 200 CFU/mL (2.3 log10
    # CFU/mL) by the M3 method.
    log_cfu <- log10(max(cfu, 1e-6))
    log_cfu ~ add(addSd)
  })
}
