Lowe_2009_omalizumab <- function() {
  description <- "Mechanism-based binding population PK/PD model for omalizumab and free / total IgE in 1928 patients (1781 with severe persistent allergic asthma across four Phase III studies plus 152 healthy atopic volunteers in a single-dose bioequivalence study; Lowe 2009). Three serum entities (free omalizumab, free IgE, and the omalizumab-IgE complex) each carry their own clearance and apparent volume of distribution and are coupled through instantaneous-equilibrium binding (law of mass action) with a baseline-IgE-dependent dissociation constant that further varies with the instantaneous total-omalizumab-to-total-IgE molar ratio. Body weight modifies all clearances, all volumes, and the IgE production rate via allometric power covariates centred at 70 kg; baseline IgE modifies CL of free IgE, IgE production rate, and Kd via power covariates centred at 365 ng/mL. Subcutaneous absorption is first-order. Disposition parameters are reported as apparent (divided by SC bioavailability f). Extends Hayashi 2007 (modellib('Hayashi_2007_omalizumab')) with (i) IIV on Kd, (ii) baseline IgE as a covariate on Kd, and (iii) bodyweight covariates on the IgE production and clearance parameters. Three observed quantities: total omalizumab, total IgE, and free IgE (all in ng/mL)."
  reference <- "Lowe PJ, Tannenbaum S, Gautier A, Jimenez P. Relationship between omalizumab pharmacokinetics, IgE pharmacodynamics and symptoms in patients with severe persistent allergic (IgE-mediated) asthma. Br J Clin Pharmacol. 2009;68(1):61-76. doi:10.1111/j.1365-2125.2009.03401.x (PMID 19660004). Extends Hayashi N et al., Br J Clin Pharmacol. 2007;63(5):548-561; see modellib('Hayashi_2007_omalizumab')."
  vignette <- "Lowe_2009_omalizumab"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ng/mL (total omalizumab, total IgE, and free IgE)"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Pretreatment (baseline) body weight; used as a per-subject fixed covariate for dose selection and as a power covariate on all clearances and volumes (CL_X/F, CL_E/F, CL_C/F, V_X/f = V_E/f, V_C/f) and on IgE production rate (R/f). All effects centred at 70 kg (Lowe 2009 Table 3 footnote *).",
      source_name        = "Body weight"
    ),
    IGE = list(
      description        = "Baseline serum total IgE concentration (pretreatment)",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL of free IgE (CL_E/F), on IgE production rate (R/f), and on the equilibrium dissociation constant (Kd). All effects centred at 365 ng/mL (Lowe 2009 Table 3 footnote dagger). Also used as the initial value for the total-IgE state total_target at t = 0: total_target(0) = (IGE / MW_IgE) * V_E. Lowe 2009 Methods state 1 IU/mL = 2.42 ng/mL.",
      source_name        = "IgE0 (baseline IgE concentration)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 1928L,
    n_studies        = 5L,
    n_observations   = 23488L,
    n_omalizumab_obs = 5938L,
    n_total_ige_obs  = 11034L,
    n_free_ige_obs   = 6156L,
    age_range        = "12-79 years (Lowe 2009 Methods 'Study design and conduct'; Table 2)",
    weight_range     = "39-150 kg",
    weight_median    = "Reference 70 kg (Lowe 2009 Table 3 footnote *)",
    sex_female_pct   = NA_real_,
    disease_state    = "Severe persistent allergic (IgE-mediated) asthma (four Phase III studies) plus healthy atopic volunteers (single-dose bioequivalence study)",
    dose_range       = "Subcutaneous 75-375 mg per dose every 2 or 4 weeks (Phase III studies, indexed on bodyweight x baseline IgE per the dosing table) plus single-dose SC 150 or 300 mg (bioequivalence study)",
    regions          = "Multinational (INNOVATE + three other Phase III) plus USA (bioequivalence study)",
    baseline_ige_range = "19-1055 IU/mL (~46-2553 ng/mL using paper's 1 IU/mL = 2.42 ng/mL conversion); reference 365 ng/mL (~150.8 IU/mL)",
    notes            = "Five clinical studies: INNOVATE [13] (28-week treatment + 16-week follow-up, severe persistent allergic asthma), two 7-month Phase III parallel-group studies in adolescents and adults with moderate-to-severe allergic asthma on ICS [23, 24], a 32-week corticosteroid-reduction pilot [25], and a single-dose bioequivalence study in healthy atopic volunteers (Novartis, unpublished). Combined dataset 23 488 observations (5938 omalizumab, 11 034 total IgE, 6156 free IgE) from 1928 patients/volunteers. Lowe 2009 Table 2."
  )

  ini({
    # Structural parameters - reference 70 kg body weight, 365 ng/mL baseline IgE.
    # All CL and V are reported as apparent (divided by SC bioavailability f), so f = 1
    # is treated implicitly. Time unit is days (the unit used in the paper's Table 3).
    lka            <- log(0.458); label("Apparent SC absorption rate constant for omalizumab (ka, 1/day)")                              # Lowe 2009 Table 3 (ka 0.458 +/- 0.0626)
    lcl            <- log(0.208); label("Apparent CL of free omalizumab at 70 kg (CL_X/F, L/day)")                                       # Lowe 2009 Table 3 (CL_X/F 0.208 +/- 0.00338)
    lcl_ige        <- log(3.85);  label("Apparent CL of free IgE at 70 kg, 365 ng/mL baseline IgE (CL_E/F, L/day)")                      # Lowe 2009 Table 3 (CL_E/F 3.85 +/- 0.155)
    lcl_complex    <- log(0.832); label("Apparent CL of omalizumab-IgE complex at 70 kg (CL_C/F, L/day)")                                # Lowe 2009 Table 3 (CL_C/F 0.832 +/- 0.0344)
    lvc            <- log(9.33);  label("Apparent V_d of free omalizumab and free IgE at 70 kg (V_X/f = V_E/f, L)")                      # Lowe 2009 Table 3 (V_X/f & V_E/f 9.33 +/- 0.147)
    lvc_complex    <- log(6.31);  label("Apparent V_d of omalizumab-IgE complex at 70 kg (V_C/f, L)")                                    # Lowe 2009 Table 3 (V_C/f 6.31 +/- 0.196)
    lp_ige         <- log(1220 / 190.07); label("Apparent IgE production rate at 70 kg, 365 ng/mL baseline IgE (R/f, nmol/day; paper: 1220 ug/day / 190 kDa IgE)") # Lowe 2009 Table 3 footnote dagger (R_E/f 1220 +/- 49.9 ug/day)
    lkd0           <- log(1.81);  label("Equilibrium dissociation constant at central = total_target, 365 ng/mL baseline IgE (Kd, nmol/L)") # Lowe 2009 Table 3 (Kd 1.81 +/- 0.0808 nmol/L)

    # Covariate exponents (power form; Lowe 2009 page 64 Eq. 3 and Table 3 covariate block).
    e_wt_cl            <- 1.00;   label("Power exponent of body weight on CL_X/F (unitless)")                                            # Lowe 2009 Table 3 (1.00 +/- 0.0662)
    e_wt_cl_ige        <- 0.499;  label("Power exponent of body weight on CL_E/F (unitless)")                                            # Lowe 2009 Table 3 (0.499 +/- 0.114)
    e_wt_cl_complex    <- 0.671;  label("Power exponent of body weight on CL_C/F (unitless)")                                            # Lowe 2009 Table 3 (0.671 +/- 0.108)
    e_wt_vc            <- 0.828;  label("Power exponent of body weight on V_X/f (= V_E/f; unitless)")                                    # Lowe 2009 Table 3 (0.828 +/- 0.0635)
    e_wt_vc_complex    <- 0.549;  label("Power exponent of body weight on V_C/f (unitless)")                                             # Lowe 2009 Table 3 (0.549 +/- 0.0936)
    e_wt_p_ige         <- 0.491;  label("Power exponent of body weight on R/f (unitless)")                                               # Lowe 2009 Table 3 (0.491 +/- 0.116)
    e_ige_cl_ige       <- 0.372;  label("Power exponent of baseline IgE on CL_E/F (unitless)")                                           # Lowe 2009 Table 3 (0.372 +/- 0.0158)
    e_ige_p_ige        <- 0.594;  label("Power exponent of baseline IgE on R/f (unitless)")                                              # Lowe 2009 Table 3 (0.594 +/- 0.0156)
    e_ige_kd0          <- 0.115;  label("Power exponent of baseline IgE on Kd (unitless)")                                               # Lowe 2009 Table 3 (0.115 +/- 0.0142)
    alpha              <- 0.0902; label("Exponent on total omalizumab to total IgE molar ratio in Kd (unitless)")                        # Lowe 2009 Table 3 (alpha 0.0902 +/- 0.0108)

    # IIV - log-normal; variance values are the NONMEM omega^2 estimates reported in
    # Table 3 (Lowe 2009 page 67 footnote 'For convenience, inter- and intra-individual
    # variances are shown as %CV as well as the original values determined by NONMEM').
    # Block covariances: Cov(eta_CLX, eta_VX) = 0.103; Cov(eta_CLC, eta_R) = -0.0101.
    etalka                       ~ 2.01                              # Lowe 2009 Table 3 (ka variance 2.01 +/- 0.940; CV 141%)
    etalcl + etalvc              ~ c(0.162,
                                     0.103, 0.0901)                  # Lowe 2009 Table 3 (CL_X 0.162 +/- 0.0196; V_X = V_E 0.0901 +/- 0.00762; cov 0.103 +/- 0.0183)
    etalcl_ige                   ~ 0.0479                            # Lowe 2009 Table 3 (CL_E variance 0.0479 +/- 0.0208; CV 23%)
    etalcl_complex + etalp_ige   ~ c(0.0649,
                                     -0.0101, 0.0701)                # Lowe 2009 Table 3 (CL_C 0.0649 +/- 0.00921; R 0.0701 +/- 0.0186; cov -0.0101 +/- 0.00687)
    etalvc_complex               ~ 0.0519                            # Lowe 2009 Table 3 (V_C variance 0.0519 +/- 0.00829; CV 26%)
    etalkd0                      ~ 0.0991                            # Lowe 2009 Table 3 (Kd variance 0.0991 +/- 0.00722; CV 31%)

    # Residual error - log-transformed both sides (LTBS): observations and outputs were
    # natural-log transformed and additive residuals fit on the log scale (Lowe 2009 page
    # 65). For small sigma, Y = F * exp(eps) is well approximated by a proportional
    # model with propSd = sigma; encoded with prop() per nlmixr2 convention.
    propSd            <- sqrt(0.0568); label("Proportional residual error on total omalizumab (fraction)")   # Lowe 2009 Table 3 (omalizumab variance 0.0568 +/- 0.00429; CV 24%)
    propSd_totalIgE   <- sqrt(0.0671); label("Proportional residual error on total IgE (fraction)")          # Lowe 2009 Table 3 (total IgE variance 0.0671 +/- 0.00324; CV 26%)
    propSd_freeIgE    <- sqrt(0.0600); label("Proportional residual error on free IgE (fraction)")           # Lowe 2009 Table 3 (free IgE variance 0.0600 +/- 0.00352; CV 24%)
  })

  model({
    # ------------------------------------------------------------------
    # Constants. Molecular weights are in kDa, which equals ng/nmol; this
    # makes 1 nmol/L * MW_kDa = MW_kDa ng/mL and 1 mg = 1000 / MW_kDa nmol
    # for dose-mass-to-mole conversion. Values match Hayashi 2007
    # (modellib('Hayashi_2007_omalizumab')) which Lowe 2009 extends:
    # MW_omalizumab approximately 150 kDa; MW_IgE approximately 190 kDa.
    # ------------------------------------------------------------------
    MWX <- 150     # omalizumab molecular weight (kDa = ng/nmol)
    MWE <- 190.07  # IgE molecular weight (kDa); 190.07 chosen so 30.3 ug/h IgE production maps to 0.1595 nmol/h (Hayashi 2007 footnote dagger)

    # ------------------------------------------------------------------
    # 1. Individual parameters - power covariates on WT (centred 70 kg)
    # and IGE (centred 365 ng/mL); Lowe 2009 Eq. 3 page 64.
    # ------------------------------------------------------------------
    ka          <- exp(lka          + etalka)
    cl          <- exp(lcl          + etalcl)         * (WT  / 70 )^e_wt_cl
    cl_ige      <- exp(lcl_ige      + etalcl_ige)     * (WT  / 70 )^e_wt_cl_ige     * (IGE / 365)^e_ige_cl_ige
    cl_complex  <- exp(lcl_complex  + etalcl_complex) * (WT  / 70 )^e_wt_cl_complex
    vc          <- exp(lvc          + etalvc)         * (WT  / 70 )^e_wt_vc
    v_ige       <- vc                                                                # V_E/f = V_X/f (Lowe 2009 Table 3 row 'Volume, omalizumab and IgE')
    vc_complex  <- exp(lvc_complex  + etalvc_complex) * (WT  / 70 )^e_wt_vc_complex
    p_ige       <- exp(lp_ige       + etalp_ige)      * (WT  / 70 )^e_wt_p_ige      * (IGE / 365)^e_ige_p_ige
    kd0         <- exp(lkd0         + etalkd0)                                       * (IGE / 365)^e_ige_kd0

    # ------------------------------------------------------------------
    # 2. Concentration-dependent equilibrium dissociation constant Kd
    # (Lowe 2009 Eq. 2 page 64): Kd = Kd0 * (XT / ET)^alpha. At t = 0 with
    # central = 0 the ratio is 0 and Kd = 0 (0^alpha = 0 for alpha > 0),
    # so no complex forms initially.
    # ------------------------------------------------------------------
    kd <- kd0 * (central / total_target)^alpha

    # ------------------------------------------------------------------
    # 3. Equilibrium-binding solution for the complex amount X_C from the
    # mass-balance quadratic
    #   X_C^2 - X_C*(XT + ET + Kd * V_X * V_E / V_C) + XT*ET = 0
    # negative-discriminant root (Lowe 2009 Eq. 2; Hayashi 2007 page 552).
    # The discriminant equals (XT - ET)^2 + Kd_eff * (2*XT + 2*ET + Kd_eff)
    # and is non-negative for all physical XT, ET, Kd, V's.
    # ------------------------------------------------------------------
    S   <- central + total_target + kd * vc * v_ige / vc_complex
    X_C <- 0.5 * (S - sqrt(S * S - 4 * central * total_target))

    # Free / complex concentrations in nmol/L (= nM).
    C_fX <- (central       - X_C) / vc
    C_fE <- (total_target  - X_C) / v_ige
    C_C  <-  X_C                  / vc_complex

    # ------------------------------------------------------------------
    # 4. ODE system (Lowe 2009 Eq. 1 page 64).
    # State units: depot in mg, central / total_target / X_C in nmol.
    # ka * depot is in mg/day; the * (1000 / MWX) factor converts to
    # nmol/day. p_ige is in nmol/day. CL * concentration (L/day * nM)
    # = nmol/day.
    # ------------------------------------------------------------------
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot * (1000 / MWX) - cl     * C_fX - cl_complex * C_C
    d/dt(total_target) <-  p_ige                     - cl_ige * C_fE - cl_complex * C_C

    # Initial total IgE amount at t = 0 (no drug, so total IgE = free IgE):
    #   total_target(0) = (IGE [ng/mL] / MWE [kDa]) * V_E [L] = nmol.
    total_target(0) <- (IGE / MWE) * v_ige

    # ------------------------------------------------------------------
    # 5. Observation outputs in ng/mL (Lowe 2009 Methods page 64: 'the
    # free and complex concentrations were ... summed to attain total
    # omalizumab and IgE concentrations'). 1 nmol/L * MW_kDa = MW_kDa
    # ng/mL.
    # ------------------------------------------------------------------
    Cc        <- (C_fX + C_C) * MWX
    totalIgE  <- (C_fE + C_C) * MWE
    freeIgE   <-  C_fE        * MWE

    Cc        ~ prop(propSd)
    totalIgE  ~ prop(propSd_totalIgE)
    freeIgE   ~ prop(propSd_freeIgE)
  })
}
