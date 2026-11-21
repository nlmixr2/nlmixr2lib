Aubry_2023_ceftazidime_avibactam_colistin_NARA1295 <- function() {
  description <- "In vitro (Klebsiella pneumoniae NARA1295, KPC-3-producing ST512 clinical isolate). Semi-mechanistic PK/PD model of static time-kill experiments with ceftazidime (in the presence of a fixed 4 mg/L avibactam) and colistin, alone and in combination: a single homogeneous bacterial population with logistic growth, sigmoidal Emax kill by each drug whose Emax is progressively reduced by an empirical two-state adaptive-resistance model per drug, Bliss independence as the additivity criterion, and a general pharmacodynamic interaction (GPDI) term through which colistin lowers the ceftazidime EC50 by up to 51.4%. Drug concentrations are static covariates (no dosing events). Sibling models for the other three isolates: Aubry_2023_ceftazidime_avibactam_colistin_NARA1584, _NARA1182 and _NARA864."
  reference <- "Aubry R, Buyck J, Prouvensier L, Decousser JW, Nordmann P, Wicha SG, Marchand S, Gregoire N. An improved PKPD modeling approach to characterize the pharmacodynamic interaction over time between ceftazidime/avibactam and colistin from in vitro time-kill experiments against multidrug-resistant Klebsiella pneumoniae isolates. Antimicrob Agents Chemother. 2023;67(10):e00301-23. doi:10.1128/aac.00301-23. Structural equations: main-text Eqs 1-6 and the NONMEM control stream of this isolate's observed interaction model (Supplementary Code S2). Parameter estimates: main-text Table 2."
  vignette <- "Aubry_2023_ceftazidime_avibactam_colistin"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Adaptive-resistance states for ceftazidime/avibactam. `aron` / `aroff` are
  # canonical and accept a registered drug suffix (as `aron_col` does), but
  # `caz` is not a registered suffix token.
  paper_specific_compartments <- c("aron_caz", "aroff_caz")

  compartmentData <- list(
    cfu = list(
      analyte = "Klebsiella pneumoniae NARA1295 (total viable bacteria)",
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
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the NARA1295 time-kill grid spanned 2-64 mg/L (Fig. 3; MIC 8 mg/L, Table 1)",
      notes = "Held constant for the 30-h experiment. The paper writes 'CZA' for the ceftazidime concentration and states that it always corresponds to ceftazidime with 4 mg/L avibactam (Discussion); avibactam is assumed to inhibit the beta-lactamases fully without a bactericidal effect of its own, so it carries no covariate. Adaptive resistance develops only while the concentration is above zero (Code S2: IF(A(1).GT.0) ADAPT_CZA=KON_CZA).",
      source_name = "CAZ_CONC"
    ),
    CONC_CST_MGL = list(
      description = "Static (time-invariant) colistin concentration in the time-kill tube. Drives the sigmoidal Emax kill by colistin, switches on adaptive resistance to it, and as the GPDI perpetrator lowers the ceftazidime EC50.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the NARA1295 time-kill grid spanned 0.03-4 mg/L (Fig. 3; MIC 0.125 mg/L, Table 1)",
      notes = "Colistin sulfate (Sigma Aldrich), concentrations expressed as colistin. Held constant for the 30-h experiment. Adaptive resistance develops only while the concentration is above zero (Code S2: IF(A(2).GT.0) ADAPT_CST=KON_CST).",
      source_name = "CST_CONC"
    )
  )

  population <- list(
    species = "in vitro (Klebsiella pneumoniae NARA1295 clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Clinical K. pneumoniae isolate from the human strain collection (2019-2021) of the National Reference Center for Emerging Antibiotic Resistance (NARA, Fribourg, Switzerland). Sequence type ST512; resistance genes KPC-3, SHV-11 (35Q), TEM-1D, truncated OmpK35 (25%) and OmpK36 with a glycine-aspartate insertion (Table 1). MICs (broth microdilution, EUCAST): ceftazidime/avibactam 8 mg/L (avibactam fixed at 4 mg/L), colistin 0.125 mg/L; susceptible to both.",
    model_system = "Static time-kill experiments over 30 h at 35 degrees C with shaking (130 rpm), 10 mL cation-adjusted Mueller-Hinton broth II, inoculum about 5 x 10^5 CFU/mL from an exponential-phase culture. Bacterial counts at 0, 2, 4, 8, 24 and 30 h by spiral plating; lower limit of quantification 200 CFU/mL (2.3 log10 CFU/mL). Growth control, colistin alone at 0.125-32 x MIC, ceftazidime/avibactam alone at 0.125-8 x MIC, and at least 16 combinations from 0.25 + 0.25 to 2 + 2 x MIC; each condition tested one to four times.",
    initial_inoculum = "Estimated 10^5.63 CFU/mL (Table 2)",
    dose_range = "Static concentrations: ceftazidime 2-64 mg/L (with avibactam 4 mg/L), colistin 0.03-4 mg/L, alone and in combination (Fig. 3)",
    notes = paste(
      "Estimated in NONMEM 7.4 with the Laplacian method (ADVAN13, TOL=9); data below the limit of quantification were handled with Beal's M3 method.",
      "The only random effect is the additive residual error on the log10 bacterial count; there is no between-experiment or between-tube variability.",
      "Built in two steps: mono-drug effects and the Bliss expected-additivity model first, then the GPDI interaction terms added and selected by likelihood ratio test (alpha 5%) to give the observed interaction model encoded here.",
      "Table 2 reports the median of a 1000-run parametric bootstrap with relative standard errors.",
      "Sibling models for the other isolates: Aubry_2023_ceftazidime_avibactam_colistin_NARA1584, _NARA1182, _NARA864."
    )
  )

  ini({
    # --- Bacterial growth (Table 2) ---
    log10inoc <- 5.63; label("Inoculum, estimated (log10 CFU/mL)")                                   # Table 2 N1295: INOC 5.63 (RSE 2%)
    kg        <- 1.40; label("Apparent bacterial growth rate constant (1/h)")                       # Table 2 N1295: K_G 1.40 (RSE 13%)
    log10bmax <- 8.80; label("Maximal bacterial population supported by the system (log10 CFU/mL)") # Table 2 N1295: Bmax 8.80 (RSE 1%)

    # --- Ceftazidime/avibactam effect and adaptive resistance (Table 2; Eqs 1-4) ---
    emax_caz  <- 4.59;  label("Maximum kill rate constant of ceftazidime/avibactam (1/h)")                     # Table 2 N1295: Emax_CZA 4.59 (RSE 19%)
    ec50_caz  <- 13.39; label("Ceftazidime concentration giving half of the maximum kill rate (mg/L)")        # Table 2 N1295: EC50_CZA 13.39 (RSE 12%)
    hill_caz  <- 2.13;  label("Sigmoidicity of the ceftazidime/avibactam effect (unitless)")                   # Table 2 N1295: gamma_CZA 2.13 (RSE 16%)
    kon_caz   <- 0.24;  label("Rate constant for development of adaptive resistance to ceftazidime/avibactam (1/h)") # Table 2 N1295: K_ON_CZA 0.24 (RSE 51%)
    armax_caz <- 0.719; label("Maximal reduction of the ceftazidime/avibactam Emax when all bacteria are adapted (fraction)") # Table 2 N1295: ARmax_CZA 71.9 percent (RSE 8%), as a fraction

    # --- Reversal of adaptive resistance (Table 2), shared by both drugs ---
    koff <- fixed(0); label("Rate constant for reversal of adaptive resistance (1/h)") # Table 2: K_OFF 0 for every isolate; Code S2 THETA(8) (0) FIX

    # --- Colistin effect and adaptive resistance (Table 2; Eqs 1-4) ---
    emax_col  <- 13.12; label("Maximum kill rate constant of colistin (1/h)")                       # Table 2 N1295: Emax_CST 13.12 (RSE 41%)
    ec50_col  <- 0.60;  label("Colistin concentration giving half of the maximum kill rate (mg/L)") # Table 2 N1295: EC50_CST 0.60 (RSE 46%)
    hill_col  <- 0.96;  label("Sigmoidicity of the colistin effect (unitless)")                     # Table 2 N1295: gamma_CST 0.96 (RSE 17%)
    kon_col   <- 0.65;  label("Rate constant for development of adaptive resistance to colistin (1/h)") # Table 2 N1295: K_ON_CST 0.65 (RSE 56%)
    armax_col <- 0.918; label("Maximal reduction of the colistin Emax when all bacteria are adapted (fraction)") # Table 2 N1295: ARmax_CST 91.8 percent (RSE 3%), as a fraction

    # --- GPDI interaction: colistin (perpetrator) lowers the ceftazidime EC50 (Table 2; Eq 6) ---
    int_col_caz     <- -0.514;    label("Maximal fractional change of the ceftazidime EC50 caused by colistin (fraction)")       # Table 2 N1295: INT_EC50 -51.4 percent (RSE 10%), as a fraction
    ec50int_col_caz <- 0.06;      label("Colistin concentration giving half of the maximal EC50 change (mg/L)")                 # Table 2 N1295: EC50_INT_EC50 0.06 (RSE 5%)
    hillint_col_caz <- fixed(10); label("Sigmoidicity of the colistin-on-ceftazidime interaction (unitless)")                   # Table 2 N1295: gamma_INT_EC50 10, not estimated; Code S2 THETA(17) (10) FIX

    # --- Residual error (Table 2; Code S2 $ERROR) ---
    # Table 2 prints RES_ADD = 1.74 'log10CFU/mL', but the value is the NONMEM
    # $SIGMA variance: Code S2 uses W = SQRT(SIGMA(1)), and the 10th-90th
    # percentile bands of the supplementary VPCs (Figs. S2, S3) match
    # SD = sqrt(RES_ADD) and not SD = RES_ADD. See the vignette.
    addSd <- sqrt(1.74); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # Table 2 N1295: RES_ADD 1.74 (RSE 7%) read as a variance; SD = sqrt(1.74) = 1.319
  })

  model({
    bmax <- 10^log10bmax

    # GPDI (Eq 6; Code S2): the victim parameter is the ceftazidime EC50,
    # theta_INT = theta_mono * (1 + INT * C_perp^gINT / (EC50_INT^gINT + C_perp^gINT)).
    fint_col_caz <- int_col_caz * CONC_CST_MGL^hillint_col_caz /
      (ec50int_col_caz^hillint_col_caz + CONC_CST_MGL^hillint_col_caz)
    ec50_caz_i <- ec50_caz * (1 + fint_col_caz)

    # Adaptation develops only while the drug is present (Code S2 $DES:
    # ADAPT_CZA = 0; IF(A(1).GT.0) ADAPT_CZA = KON_CZA, and likewise for colistin).
    kon_caz_i <- kon_caz * (CONC_CAZ_MGL > 0)
    kon_col_i <- kon_col * (CONC_CST_MGL > 0)

    # Sigmoidal Emax kill, Emax reduced in proportion to the adapted fraction (Eq 4).
    kill_caz <- emax_caz * (1 - armax_caz * aron_caz) * CONC_CAZ_MGL^hill_caz /
      (ec50_caz_i^hill_caz + CONC_CAZ_MGL^hill_caz)
    kill_col <- emax_col * (1 - armax_col * aron_col) * CONC_CST_MGL^hill_col /
      (ec50_col^hill_col + CONC_CST_MGL^hill_col)

    # Bliss independence (Code S2 $DES): each kill rate is normalised by the
    # larger of the two typical Emax values, combined as Ea + Eb - Ea * Eb, and
    # rescaled; algebraically kill_caz + kill_col - kill_caz * kill_col / emax_norm.
    emax_norm <- max(emax_caz, emax_col)
    kill_comb <- kill_caz + kill_col - kill_caz * kill_col / emax_norm

    # Logistic growth on the linear CFU scale (Eq 5; Code S2 DADT(3)).
    d/dt(cfu) <- (kg * (1 - cfu / bmax) - kill_comb) * cfu

    # Two-state adaptive resistance per drug (Eqs 2-3; Code S2 DADT(4)-DADT(7)).
    d/dt(aron_caz)  <- kon_caz_i * aroff_caz - koff * aron_caz
    d/dt(aroff_caz) <- koff * aron_caz - kon_caz_i * aroff_caz
    d/dt(aron_col)  <- kon_col_i * aroff_col - koff * aron_col
    d/dt(aroff_col) <- koff * aron_col - kon_col_i * aroff_col

    # Initial conditions (Code S2 A_0): inoculum, all bacteria non-adapted.
    cfu(0) <- 10^log10inoc
    aron_caz(0) <- 0
    aroff_caz(0) <- 1
    aron_col(0) <- 0
    aroff_col(0) <- 1

    # Observation: log10 CFU/mL with the 1e-6 CFU/mL floor of Code S2 $ERROR
    # (BTOT_ERR). The fitted data were censored below 200 CFU/mL (2.3 log10
    # CFU/mL) by the M3 method.
    log_cfu <- log10(max(cfu, 1e-6))
    log_cfu ~ add(addSd)
  })
}
