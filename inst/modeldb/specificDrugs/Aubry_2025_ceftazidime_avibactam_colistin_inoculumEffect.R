Aubry_2025_ceftazidime_avibactam_colistin_inoculumEffect <- function() {
  description <- "In vitro (Klebsiella pneumoniae NARA864, KPC-3-producing clinical isolate). Semi-mechanistic PK/PD model of static time-kill experiments with ceftazidime (in the presence of a fixed 4 mg/L avibactam) and colistin, alone and in combination, fitted across four starting inocula (10^4 to 10^8 CFU/mL) to describe the inoculum effect (the 'IE model'). A single homogeneous bacterial population with logistic growth is killed by each drug through a sigmoidal Emax term whose Emax is progressively reduced by an empirical two-state adaptive-resistance model, combined under Bliss independence, with a general pharmacodynamic interaction (GPDI) term through which colistin lowers the ceftazidime/avibactam Emax. The inoculum effect enters as a power function of the starting inoculum: increasing inoculum reduces the ceftazidime/avibactam Emax and increases the colistin EC50. Drug concentrations and the starting inoculum are static covariates (no dosing events). Companion no-inoculum-effect model fitted to the 5e5 CFU/mL experiments only: Aubry_2025_ceftazidime_avibactam_colistin_noInoculumEffect."
  reference <- "Aubry R, Buyck JM, Chauzy A, Prouvensier L, Decousser JW, Nordmann P, Wicha SG, Marchand S, Gregoire N. PKPD modeling of the inoculum effect of combined ceftazidime/avibactam and colistin against KPC-3 Klebsiella pneumoniae isolate. Antimicrob Agents Chemother. 2025;69(5):e01797-24. doi:10.1128/aac.01797-24. Structural equations: main-text Eq 1 and the NONMEM control stream of the observed interaction model with inoculum effect (Supplemental Code S1); PD interaction equation in Supplemental Text S1. Parameter estimates: Supplemental Table S2 (final estimates; the Code S1 $THETA block lists initial values)."
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
      notes = "Held constant for the 30-h experiment. The paper writes 'CZA' for the ceftazidime concentration and states that it always corresponds to ceftazidime with 4 mg/L avibactam; avibactam is assumed to inhibit the beta-lactamases fully without a bactericidal effect of its own, so it carries no covariate. Adaptive resistance develops only while the concentration is above zero (Code S1 $DES: IF(A(1).GT.0) ADAPT_CZA=KON_CZA). In the clinical-concentration simulation (Fig 7) this covariate is driven time-varying by the free ceftazidime PK profile.",
      source_name = "CAZ_CONC"
    ),
    CONC_CST_MGL = list(
      description = "Static (time-invariant) colistin concentration in the time-kill tube. Drives the sigmoidal Emax kill by colistin, switches on adaptive resistance to it, and as the GPDI perpetrator lowers the ceftazidime/avibactam Emax.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the time-kill grid spanned 0.06-8 mg/L (Methods; MIC 0.25 mg/L)",
      notes = "Colistin (Sigma Aldrich), concentrations expressed as colistin. Held constant for the 30-h experiment. Adaptive resistance develops only while the concentration is above zero (Code S1 $DES: IF(A(2).GT.0) ADAPT_CST=KON_CST). In the clinical-concentration simulation (Fig 7) this covariate is driven time-varying by the free colistin PK profile.",
      source_name = "CST_CONC"
    ),
    INOC_LOG10CFU = list(
      description = "Starting (time-zero) inoculum size of the in-vitro time-kill experiment, on the log10 CFU/mL scale. Sets the initial bacterial state (10^INOC_LOG10CFU CFU/mL) and drives the inoculum effect: it reduces the ceftazidime/avibactam Emax and increases the colistin EC50 through power functions of INOC_LOG10CFU / log10(5e5).",
      units = "log10 CFU/mL",
      type = "continuous",
      reference_category = "n/a -- the standard inoculum is log10(5e5) = 5.699 log10 CFU/mL, where the inoculum-effect power functions equal 1 (no effect)",
      notes = "Estimated per experiment from the measured time-zero CFU counts. The fit used four nominal inocula with estimated sizes 3.86, 5.54, 6.65 and 7.76 log10 CFU/mL (Table S2, INOC*; nominal 10^4, 5e5, 10^7 and 10^8 CFU/mL). Inter-experimental variability around the inoculum is carried by the fixed random effect etainoc (IEV_INOC = 4.30% CV, Code S1 $OMEGA 0.00188 FIX).",
      source_name = "INOC"
    )
  )

  population <- list(
    species = "in vitro (Klebsiella pneumoniae NARA864 clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Clinical KPC-3-producing K. pneumoniae isolate NARA864, obtained from the National Reference Center for Emerging Antibiotic Resistance (Fribourg, Switzerland). Susceptible to ceftazidime/avibactam (MIC 2 mg/L with avibactam fixed at 4 mg/L) and to colistin (MIC 0.25 mg/L). The same isolate was used in the earlier interaction study (Aubry et al. 2023).",
    model_system = "Static time-kill experiments over 30 h at 35 degrees C with shaking (130 rpm), 10 mL cation-adjusted Mueller-Hinton broth II, at four starting inocula (10^4, 5e5, 10^7 and 10^8 CFU/mL). Ceftazidime 0.5-8 mg/L (with avibactam 4 mg/L) or colistin 0.06-8 mg/L alone, and combinations from 0.5+0.06 to 4+0.5 mg/L. Bacterial counts at 0, 2, 4, 8, 24 and 30 h by spiral plating; lower limit of quantification 200 CFU/mL (2.3 log10 CFU/mL). Each condition tested one to four times.",
    inoculum_range = "Estimated 3.86, 5.54, 6.65 and 7.76 log10 CFU/mL (Table S2) for nominal 10^4, 5e5, 10^7 and 10^8 CFU/mL",
    dose_range = "Static concentrations: ceftazidime 0.5-8 mg/L (with avibactam 4 mg/L), colistin 0.06-8 mg/L, alone and in combination (Methods)",
    notes = paste(
      "Estimated in NONMEM 7.4 with the Laplacian method (ADVAN13, TOL=9); data below the limit of quantification were handled with Beal's M3 method.",
      "The only random effect is the fixed inter-experimental variability on the inoculum (etainoc); the residual is an additive error on the log10 bacterial count.",
      "Built in two steps: a no-inoculum-effect model fitted to the standard 5e5 CFU/mL data (reusing Aubry et al. 2023), then the inoculum effect added and fitted to all four inocula.",
      "The inoculum effect was best described as a reduction of the ceftazidime/avibactam Emax and an increase of the colistin EC50 with increasing inoculum; no interaction between the two inoculum effects was significant (E_CZA_IE_CST = E_CST_IE_CZA = 0, not included here).",
      "Table 1 and Table S2 give the final estimates; the Code S1 $THETA block lists initial values."
    )
  )

  ini({
    # --- Bacterial growth (Table S2) ---
    kg        <- 1.77; label("Apparent bacterial growth rate constant (1/h)")                       # Table S2: K_G 1.77 (RSE 6%)
    log10bmax <- 9.14; label("Maximal bacterial population supported by the system (log10 CFU/mL)") # Table S2: Bmax 9.14 (RSE 1%)

    # --- Ceftazidime/avibactam effect and adaptive resistance (Table S2; Eq 1) ---
    emax_caz  <- 5.59;  label("Maximum kill rate constant of ceftazidime/avibactam at the standard inoculum (1/h)") # Table S2: Emax_CZA 5.59 (RSE 3%)
    ec50_caz  <- 0.40;  label("Ceftazidime concentration giving half of the maximum kill rate (mg/L)")              # Table S2: EC50_CZA 0.40 (RSE 6%)
    hill_caz  <- 1.26;  label("Sigmoidicity of the ceftazidime/avibactam effect (unitless)")                        # Table S2: gamma_CZA 1.26 (RSE 8%)
    kon_caz   <- 0.31;  label("Rate constant for development of adaptive resistance to ceftazidime/avibactam (1/h)") # Table S2: K_ON_CZA 0.31 (RSE 7%)
    armax_caz <- 0.729; label("Maximal reduction of the ceftazidime/avibactam Emax when all bacteria are adapted (fraction)") # Table S2: ARmax_CZA 72.9 percent (RSE 2%), as a fraction

    # --- Reversal of adaptive resistance (Table S2), shared by both drugs ---
    koff <- fixed(0); label("Rate constant for reversal of adaptive resistance (1/h)") # Table S2: K_OFF 0 fixed; Code S1 THETA(11) (0) FIX

    # --- Colistin effect and adaptive resistance (Table S2; Eq 1) ---
    emax_col  <- 13.79; label("Maximum kill rate constant of colistin (1/h)")                       # Table S2: Emax_CST 13.79 (RSE 11%)
    ec50_col  <- 0.26;  label("Colistin concentration giving half of the maximum kill rate at the standard inoculum (mg/L)") # Table S2: EC50_CST 0.26 (RSE 9%)
    hill_col  <- 1.21;  label("Sigmoidicity of the colistin effect (unitless)")                     # Table S2: gamma_CST 1.21 (RSE 4%)
    kon_col   <- 1.28;  label("Rate constant for development of adaptive resistance to colistin (1/h)") # Table S2: K_ON_CST 1.28 (RSE 14%)
    armax_col <- 0.891; label("Maximal reduction of the colistin Emax when all bacteria are adapted (fraction)") # Table S2: ARmax_CST 89.1 percent (RSE 1%), as a fraction

    # --- Inoculum effect, power model (Eq 1; Code S1 $DES) ---
    # theta_INOC = theta_STD * (INOC / log10(5e5))^gamma_IE. A negative exponent
    # reduces the ceftazidime/avibactam Emax with increasing inoculum; a positive
    # exponent increases the colistin EC50 with increasing inoculum.
    ie_emax_caz <- -0.21; label("Inoculum-effect exponent on the ceftazidime/avibactam Emax (unitless)") # Table S2 / Table 1: gamma_IE_CZA -0.21 (RSE 13%)
    ie_ec50_col <- 1.53;  label("Inoculum-effect exponent on the colistin EC50 (unitless)")              # Table S2 / Table 1: gamma_IE_CST 1.53 (RSE 10%)

    # --- GPDI interaction: colistin (perpetrator) lowers the ceftazidime/avibactam Emax (Table S2; Text S1) ---
    int_col_caz     <- -0.63;       label("Maximal fractional change of the ceftazidime/avibactam Emax caused by colistin (fraction)") # Table S2: INT_Emax -63.0 percent (RSE 11%), as a fraction
    ec50int_col_caz <- fixed(0.30); label("Colistin concentration giving half of the maximal Emax change (mg/L)")                      # Table S2: EC50_INT_Emax 0.30 fixed
    hillint_col_caz <- fixed(1);    label("Sigmoidicity of the colistin-on-ceftazidime/avibactam interaction (unitless)")              # Table S2: gamma_INT_Emax 1 fixed

    # --- Inter-experimental variability on the inoculum (Table S2; Code S1 $OMEGA) ---
    etainoc ~ fixed(0.00188) # IEV_INOC 4.30% CV (Table S2); Code S1 $OMEGA 0.00188; sqrt(exp(0.00188)-1) = 4.34% CV

    # --- Residual error (Table S2; Code S1 $ERROR) ---
    # Table S2 prints RES_ADD = 1.73 'log10 CFU/mL', but the value is the NONMEM
    # $SIGMA variance: Code S1 uses W = SQRT(SIGMA(1)), and the 10th-90th
    # percentile bands of the supplementary VPCs (Figs S1, S3-S5) are about
    # 2.9 log10 wide, matching SD = sqrt(RES_ADD) = 1.32 and not SD = 1.73.
    addSd <- sqrt(1.73); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # Table S2: RES_ADD 1.73 (RSE 4%) read as a variance; SD = sqrt(1.73) = 1.32
  })

  model({
    bmax <- 10^log10bmax

    # Per-experiment inoculum with inter-experimental variability (Code S1 $PK:
    # INOC = INOC_TV * EXP(ETA(1)); the eta multiplies the log10-scale inoculum
    # and propagates to both the initial state and the inoculum-effect terms).
    inoc_i <- INOC_LOG10CFU * exp(etainoc)
    inoc_ratio <- inoc_i / log10(5e5)

    # Inoculum effect, power model (Eq 1; Code S1 $DES).
    emax_caz_ie <- emax_caz * inoc_ratio^ie_emax_caz
    ec50_col_ie <- ec50_col * inoc_ratio^ie_ec50_col

    # GPDI (Text S1): the victim parameter is the ceftazidime/avibactam Emax,
    # theta_INT = theta_mono * (1 + INT * C_perp^gINT / (EC50_INT^gINT + C_perp^gINT)).
    fint_col_caz <- int_col_caz * CONC_CST_MGL^hillint_col_caz /
      (ec50int_col_caz^hillint_col_caz + CONC_CST_MGL^hillint_col_caz)
    emax_caz_i <- emax_caz_ie * (1 + fint_col_caz)

    # Adaptation develops only while the drug is present (Code S1 $DES:
    # ADAPT_CZA = 0; IF(A(1).GT.0) ADAPT_CZA = KON_CZA, and likewise for colistin).
    kon_caz_i <- kon_caz * (CONC_CAZ_MGL > 0)
    kon_col_i <- kon_col * (CONC_CST_MGL > 0)

    # Sigmoidal Emax kill, Emax reduced in proportion to the adapted fraction;
    # the colistin EC50 carries the inoculum effect.
    kill_caz <- emax_caz_i * (1 - armax_caz * aron_caz) * CONC_CAZ_MGL^hill_caz /
      (ec50_caz^hill_caz + CONC_CAZ_MGL^hill_caz)
    kill_col <- emax_col * (1 - armax_col * aron_col) * CONC_CST_MGL^hill_col /
      (ec50_col_ie^hill_col + CONC_CST_MGL^hill_col)

    # Bliss independence (Code S1 $DES): each kill rate is normalised by the
    # larger of the two typical (base) Emax values, combined as Ea + Eb - Ea*Eb,
    # and rescaled; algebraically kill_caz + kill_col - kill_caz * kill_col / emax_norm.
    # The normaliser uses the mono-drug Emax values, as in Code S1, where EMAX is
    # chosen from the THETAs before any inoculum or interaction term is applied.
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
