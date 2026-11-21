Aubry_2023_ceftazidime_avibactam_colistin_NARA864 <- function() {
  description <- "In vitro (Klebsiella pneumoniae NARA864, KPC-3-producing ST348 clinical isolate). Semi-mechanistic PK/PD model of static time-kill experiments with ceftazidime (in the presence of a fixed 4 mg/L avibactam) and colistin, alone and in combination: a single homogeneous bacterial population with logistic growth, sigmoidal Emax kill by each drug whose Emax is progressively reduced by an empirical two-state adaptive-resistance model per drug, Bliss independence as the additivity criterion, and a general pharmacodynamic interaction (GPDI) term through which colistin lowers the ceftazidime/avibactam Emax by up to 65.8%. Drug concentrations are static covariates (no dosing events). Sibling models for the other three isolates: Aubry_2023_ceftazidime_avibactam_colistin_NARA1295, _NARA1584 and _NARA1182."
  reference <- "Aubry R, Buyck J, Prouvensier L, Decousser JW, Nordmann P, Wicha SG, Marchand S, Gregoire N. An improved PKPD modeling approach to characterize the pharmacodynamic interaction over time between ceftazidime/avibactam and colistin from in vitro time-kill experiments against multidrug-resistant Klebsiella pneumoniae isolates. Antimicrob Agents Chemother. 2023;67(10):e00301-23. doi:10.1128/aac.00301-23. Structural equations: main-text Eqs 1-6 and the NONMEM control stream of this isolate's observed interaction model (Supplementary Code S2). Parameter estimates: main-text Table 2."
  vignette <- "Aubry_2023_ceftazidime_avibactam_colistin"
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
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the NARA864 time-kill grid spanned 0.25-16 mg/L (Fig. S3; MIC 2 mg/L, Table 1)",
      notes = "Held constant for the 30-h experiment. The paper writes 'CZA' for the ceftazidime concentration and states that it always corresponds to ceftazidime with 4 mg/L avibactam (Discussion); avibactam is assumed to inhibit the beta-lactamases fully without a bactericidal effect of its own, so it carries no covariate. Adaptive resistance develops only while the concentration is above zero (Code S2, the NARA1295 control stream; the paper uses the same structure for all four isolates: IF(A(1).GT.0) ADAPT_CZA=KON_CZA).",
      source_name = "CAZ_CONC"
    ),
    CONC_CST_MGL = list(
      description = "Static (time-invariant) colistin concentration in the time-kill tube. Drives the sigmoidal Emax kill by colistin, switches on adaptive resistance to it, and as the GPDI perpetrator lowers the ceftazidime/avibactam Emax.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L is the drug-free growth control; the NARA864 time-kill grid spanned 0.03-8 mg/L (Fig. S3; MIC 0.25 mg/L, Table 1)",
      notes = "Colistin sulfate (Sigma Aldrich), concentrations expressed as colistin. Held constant for the 30-h experiment. Adaptive resistance develops only while the concentration is above zero (Code S2: IF(A(2).GT.0) ADAPT_CST=KON_CST).",
      source_name = "CST_CONC"
    )
  )

  population <- list(
    species = "in vitro (Klebsiella pneumoniae NARA864 clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    disease_state = "Clinical K. pneumoniae isolate from the human strain collection (2019-2021) of the National Reference Center for Emerging Antibiotic Resistance (NARA, Fribourg, Switzerland). Sequence type ST348; resistance genes KPC-3, SHV-11 (35Q), TEM-1D and OXA-9, with no porin alteration reported (Table 1). MICs (broth microdilution, EUCAST): ceftazidime/avibactam 2 mg/L (avibactam fixed at 4 mg/L), colistin 0.25 mg/L; susceptible to both.",
    model_system = "Static time-kill experiments over 30 h at 35 degrees C with shaking (130 rpm), 10 mL cation-adjusted Mueller-Hinton broth II, inoculum about 5 x 10^5 CFU/mL from an exponential-phase culture. Bacterial counts at 0, 2, 4, 8, 24 and 30 h by spiral plating; lower limit of quantification 200 CFU/mL (2.3 log10 CFU/mL). Growth control, colistin alone at 0.125-32 x MIC, ceftazidime/avibactam alone at 0.125-8 x MIC, and at least 16 combinations from 0.25 + 0.25 to 2 + 2 x MIC; each condition tested one to four times.",
    initial_inoculum = "Estimated 10^5.68 CFU/mL (Table 2)",
    dose_range = "Static concentrations: ceftazidime 0.25-16 mg/L (with avibactam 4 mg/L), colistin 0.03-8 mg/L, alone and in combination (Fig. S3)",
    notes = paste(
      "Estimated in NONMEM 7.4 with the Laplacian method (ADVAN13, TOL=9); data below the limit of quantification were handled with Beal's M3 method.",
      "The only random effect is the additive residual error on the log10 bacterial count; there is no between-experiment or between-tube variability.",
      "Built in two steps: mono-drug effects and the Bliss expected-additivity model first, then the GPDI interaction terms added and selected by likelihood ratio test (alpha 5%) to give the observed interaction model encoded here.",
      "Table 2 reports the median of a 1000-run parametric bootstrap with relative standard errors.",
      "Sibling models for the other isolates: Aubry_2023_ceftazidime_avibactam_colistin_NARA1295, _NARA1584, _NARA1182."
    )
  )

  ini({
    # --- Bacterial growth (Table 2) ---
    log10inoc <- 5.68; label("Inoculum, estimated (log10 CFU/mL)")                                   # Table 2 N864: INOC 5.68 (RSE 2%)
    kg        <- fixed(1.84); label("Apparent bacterial growth rate constant (1/h)")                # Table 2 N864: K_G 1.84, not estimated
    log10bmax <- 9.01; label("Maximal bacterial population supported by the system (log10 CFU/mL)") # Table 2 N864: Bmax 9.01 (RSE 2%)

    # --- Ceftazidime/avibactam effect and adaptive resistance (Table 2; Eqs 1-4) ---
    emax_caz  <- 9.17;  label("Maximum kill rate constant of ceftazidime/avibactam (1/h)")                     # Table 2 N864: Emax_CZA 9.17 (RSE 16%)
    ec50_caz  <- 0.46;  label("Ceftazidime concentration giving half of the maximum kill rate (mg/L)")        # Table 2 N864: EC50_CZA 0.46 (RSE 6%)
    hill_caz  <- 1.79;  label("Sigmoidicity of the ceftazidime/avibactam effect (unitless)")                   # Table 2 N864: gamma_CZA 1.79 (RSE 9%)
    kon_caz   <- 0.66;  label("Rate constant for development of adaptive resistance to ceftazidime/avibactam (1/h)") # Table 2 N864: K_ON_CZA 0.66 (RSE 25%)
    armax_caz <- 0.825; label("Maximal reduction of the ceftazidime/avibactam Emax when all bacteria are adapted (fraction)") # Table 2 N864: ARmax_CZA 82.5 percent (RSE 3%), as a fraction

    # --- Reversal of adaptive resistance (Table 2), shared by both drugs ---
    koff <- fixed(0); label("Rate constant for reversal of adaptive resistance (1/h)") # Table 2: K_OFF 0 for every isolate; Code S2 THETA(8) (0) FIX

    # --- Colistin effect and adaptive resistance (Table 2; Eqs 1-4) ---
    emax_col  <- 8.81;  label("Maximum kill rate constant of colistin (1/h)")                       # Table 2 N864: Emax_CST 8.81 (RSE 31%)
    ec50_col  <- 0.21;  label("Colistin concentration giving half of the maximum kill rate (mg/L)") # Table 2 N864: EC50_CST 0.21 (RSE 11%)
    hill_col  <- 1.16;  label("Sigmoidicity of the colistin effect (unitless)")                     # Table 2 N864: gamma_CST 1.16 (RSE 6%)
    kon_col   <- 0.75;  label("Rate constant for development of adaptive resistance to colistin (1/h)") # Table 2 N864: K_ON_CST 0.75 (RSE 47%)
    armax_col <- 0.803; label("Maximal reduction of the colistin Emax when all bacteria are adapted (fraction)") # Table 2 N864: ARmax_CST 80.3 percent (RSE 6%), as a fraction

    # --- GPDI interaction: colistin (perpetrator) lowers the ceftazidime/avibactam Emax (Table 2; Eq 6) ---
    # For this isolate the interaction acts on Emax rather than on EC50
    # (Results: 'In the K. pneumoniae NARA864, CST reduced the Emax of CZA by
    # 65% at maximum'). Table 2 gives no RSE for EC50_INT_Emax or gamma_INT_Emax.
    int_col_caz     <- -0.658;      label("Maximal fractional change of the ceftazidime/avibactam Emax caused by colistin (fraction)") # Table 2 N864: INT_Emax -65.8 percent (RSE 16%), as a fraction
    ec50int_col_caz <- fixed(0.39); label("Colistin concentration giving half of the maximal Emax change (mg/L)")                      # Table 2 N864: EC50_INT_Emax 0.39, not estimated
    hillint_col_caz <- fixed(1);    label("Sigmoidicity of the colistin-on-ceftazidime/avibactam interaction (unitless)")              # Table 2 N864: gamma_INT_Emax 1, not estimated

    # --- Residual error (Table 2; Code S2 $ERROR) ---
    # Table 2 prints RES_ADD = 2.25 'log10CFU/mL', but the value is the NONMEM
    # $SIGMA variance: Code S2 uses W = SQRT(SIGMA(1)), and the 10th-90th
    # percentile bands of the supplementary VPCs (Figs. S2, S3) match
    # SD = sqrt(RES_ADD) and not SD = RES_ADD. See the vignette.
    addSd <- sqrt(2.25); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # Table 2 N864: RES_ADD 2.25 (RSE 6%) read as a variance; SD = sqrt(2.25) = 1.5
  })

  model({
    bmax <- 10^log10bmax

    # GPDI (Eq 6): the victim parameter is the ceftazidime/avibactam Emax,
    # theta_INT = theta_mono * (1 + INT * C_perp^gINT / (EC50_INT^gINT + C_perp^gINT)).
    fint_col_caz <- int_col_caz * CONC_CST_MGL^hillint_col_caz /
      (ec50int_col_caz^hillint_col_caz + CONC_CST_MGL^hillint_col_caz)
    emax_caz_i <- emax_caz * (1 + fint_col_caz)

    # Adaptation develops only while the drug is present (Code S2 $DES:
    # ADAPT_CZA = 0; IF(A(1).GT.0) ADAPT_CZA = KON_CZA, and likewise for colistin).
    kon_caz_i <- kon_caz * (CONC_CAZ_MGL > 0)
    kon_col_i <- kon_col * (CONC_CST_MGL > 0)

    # Sigmoidal Emax kill, Emax reduced in proportion to the adapted fraction (Eq 4).
    kill_caz <- emax_caz_i * (1 - armax_caz * aron_caz) * CONC_CAZ_MGL^hill_caz /
      (ec50_caz^hill_caz + CONC_CAZ_MGL^hill_caz)
    kill_col <- emax_col * (1 - armax_col * aron_col) * CONC_CST_MGL^hill_col /
      (ec50_col^hill_col + CONC_CST_MGL^hill_col)

    # Bliss independence (Code S2 $DES): each kill rate is normalised by the
    # larger of the two typical Emax values, combined as Ea + Eb - Ea * Eb, and
    # rescaled; algebraically kill_caz + kill_col - kill_caz * kill_col / emax_norm.
    # The normaliser uses the mono-drug Emax values, as in Code S2, where
    # EMAX is chosen from the THETAs before any interaction term is applied.
    # The paper prints no control stream for an Emax-type interaction; see the
    # vignette for this reading.
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
