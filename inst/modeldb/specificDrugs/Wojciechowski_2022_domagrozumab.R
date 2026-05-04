Wojciechowski_2022_domagrozumab <- function() {
  description <- "Quasi-steady-state TMDD population PK/PD model for domagrozumab (anti-myostatin IgG1) in healthy adult volunteers and pediatric patients with Duchenne muscular dystrophy (Wojciechowski 2022): two-compartment IV/SC drug disposition with parallel linear and Michaelis-Menten elimination, a synthesis-degradation total-myostatin compartment with drug-mediated internalization, and a study-population covariate (DIS_DMD) shifting myostatin baseline and turnover."
  reference <- "Wojciechowski J, Purohit VS, Nawarskas R, Marshall S, Charnas L, Bhattacharya I. Population pharmacokinetic-pharmacodynamic analysis of domagrozumab in pediatric patients with Duchenne muscular dystrophy. Clin Transl Sci. 2022;15(12):2939-2952. doi:10.1111/cts.13418"
  vignette <- "Wojciechowski_2022_domagrozumab"
  units <- list(time = "hour", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "All structural PK parameters in the source paper are reported per kilogram (CL, Q in L/hour/kg; V1, V2 in L/kg; Vmax in nM/hour/kg). They are scaled by WT in model() to convert to a per-subject basis so that doses can be supplied in mg. Reference body weight is implicit (linear weight scaling, exponent = 1); the population median was 76.9 kg in healthy adults and 28.9 kg in pediatric DMD patients (Table 1).",
      source_name        = "BWT"
    ),
    DIS_DMD = list(
      description        = "Study-population indicator: 1 = pediatric patient with Duchenne muscular dystrophy, 0 = healthy adult volunteer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy adult volunteer)",
      notes              = "Source paper covariate is SPOP (study population). Eq. 7 of Wojciechowski 2022 defines COVSPOP = 1 for healthy adult volunteers and COVSPOP = 1 + theta_SPOP for DMD pediatric patients, applied multiplicatively to the affected typical value (Eq. 8). DIS_DMD encodes the indicator in the orientation that matches the source equation (1 -> add the theta effect). Effects in the final model: theta_SPOP_BASE = -0.641 on baseline myostatin (BASE) and theta_SPOP_kdegkint = -0.900 jointly on kdeg and kint.",
      source_name        = "SPOP"
    )
  )

  population <- list(
    n_subjects     = 193L,
    n_studies      = 2L,
    age_range      = "6-61 years (healthy adults 19-61; DMD pediatric 6-15)",
    age_median     = "37.0 years (HV); 9.0 years (DMD)",
    weight_range   = "14.8-99.8 kg (HV 56.8-99.8; DMD 14.8-86.4)",
    weight_median  = "76.9 kg (HV); 28.9 kg (DMD)",
    sex_female_pct = 3.6,
    race_ethnicity = c(White = 68.4, Black = 12.4, Asian = 9.3, Other = 9.8),
    disease_state  = "Pooled healthy adult volunteers (Phase 1, NCT01616277, n = 73) and ambulatory boys with Duchenne muscular dystrophy (Phase 2, NCT02310763, n = 120). All DMD patients were male (100%); HV cohort 90.4% male.",
    dose_range     = "Healthy adults: single IV doses 1, 3, 10, 20, 40 mg/kg over 2-h infusion; single SC 3 mg/kg; or 10 mg/kg IV every 2 weeks x 3 doses. DMD patients: dose-escalation 5, 20, 40 mg/kg IV every 4 weeks x 4 doses each level (period 1, 16 weeks per level); period 2 (48 weeks) continued the highest tolerated dose.",
    regions        = "Multi-regional Phase 1 and Phase 2 trials (Pfizer-sponsored).",
    n_observations_drug      = 5181L,
    n_observations_myostatin = 8001L,
    notes          = "Demographics from Wojciechowski 2022 Table 1 (combined HV + DMD pediatric column). 5.79% of free domagrozumab and 0.075% of total myostatin observations were below LLOQ. Anti-drug antibodies were negative or below quantification in the entire cohort, so ADA was not tested as a covariate (paper Methods)."
  )

  ini({
    # ---- Structural PK (Wojciechowski 2022 Table 2, "Final model estimate" column) ----
    # Per-kg parameterization (the per-kg form is the paper's final model; the
    # allometric column in Table 2 is an alternative comparison and is not used
    # here). All clearances are L/hour/kg and volumes are L/kg; the model() block
    # multiplies by the WT covariate to obtain per-subject quantities so that
    # doses can be supplied in mg.
    lcl     <- log(0.0000982); label("Linear clearance per kilogram (L/hour/kg)")                  # Wojciechowski 2022 Table 2: CL = 0.0000982 L/hour/kg
    lvc     <- log(0.0415);    label("Central volume of distribution per kilogram (L/kg)")         # Wojciechowski 2022 Table 2: V1 = 0.0415 L/kg
    lq      <- log(0.000306);  label("Intercompartmental clearance per kilogram (L/hour/kg)")      # Wojciechowski 2022 Table 2: Q  = 0.000306 L/hour/kg
    lvp     <- log(0.0416);    label("Peripheral volume of distribution per kilogram (L/kg)")      # Wojciechowski 2022 Table 2: V2 = 0.0416 L/kg
    lka     <- log(0.00769);   label("First-order subcutaneous absorption rate (1/hour)")          # Wojciechowski 2022 Table 2: ka = 0.00769 1/hour
    lfdepot <- log(0.858);     label("Subcutaneous bioavailability (fraction)")                    # Wojciechowski 2022 Table 2: F  = 0.858

    # ---- Michaelis-Menten (target-mediated) elimination from central ----
    # Vmax in nM/hour/kg (paper Methods); km in nM. Inside model() the MM rate
    # is computed in nM/hour after WT scaling and then converted to mg/hour for
    # the amount-based ODE using the assumed domagrozumab molecular weight.
    lvmax   <- log(0.00251);   label("Maximum nonlinear elimination rate per kilogram (nM/hour/kg)") # Wojciechowski 2022 Table 2: Vmax = 0.00251 nM/hour/kg
    lkm     <- log(12.2);      label("Michaelis-Menten constant (nM)")                              # Wojciechowski 2022 Table 2: km   = 12.2 nM

    # ---- Total-myostatin compartment (Eqs 4-6, Table 2) ----
    # Total myostatin (free + drug-myostatin complex) tracked as a single state
    # with quasi-steady-state binding (kss). ksyn is derived as base * kdeg
    # (Eq. 5) so that the compartment starts at steady state.
    lbase   <- log(0.156);     label("Baseline total myostatin in healthy adult volunteers (nM)") # Wojciechowski 2022 Table 2: BASE = 0.156 nM
    lkdeg   <- log(0.0381);    label("First-order myostatin degradation rate (1/hour)")            # Wojciechowski 2022 Table 2: kdeg = 0.0381 1/hour
    lkint   <- log(0.00716);   label("First-order drug-myostatin complex internalization rate (1/hour)") # Wojciechowski 2022 Table 2: kint = 0.00716 1/hour
    lkss    <- log(7.76);      label("Quasi-steady-state binding constant for drug-myostatin (nM)") # Wojciechowski 2022 Table 2: kSS  = 7.76 nM

    # ---- Structural covariate effects ----
    # Eq. 7: COVSPOP = 1 (HV) or 1 + theta (DMD pediatric); Eq. 8 applies
    # COVSPOP multiplicatively to the typical value, so theta_SPOP_X is added
    # to 1 (not exponentiated) when DIS_DMD == 1.
    e_dmd_base     <- -0.641; label("Effect of DMD pediatric population on BASE: COVSPOP = 1 + theta when DIS_DMD = 1") # Wojciechowski 2022 Table 2: theta_SPOP_BASE = -0.641
    e_dmd_kdegkint <- -0.900; label("Effect of DMD pediatric population on kdeg and kint (shared): COVSPOP = 1 + theta when DIS_DMD = 1") # Wojciechowski 2022 Table 2: theta_SPOP_(kdeg,kint) = -0.900

    # Ratio of SD for eta_kint relative to eta_kdeg (Wojciechowski 2022 Table 2,
    # "Ratio of SD for eta_kint relative to eta_kdeg" row). The paper places a
    # single random effect on the kdeg-kint axis (column "omega kdeg-kint");
    # eta_kdeg enters kdeg directly and enters kint scaled by this ratio.
    e_ratio_kdegkint <- -0.295; label("Ratio of SD for eta_kint relative to eta_kdeg (unitless)") # Wojciechowski 2022 Table 2: ratio = -0.295

    # ---- IIV: 5x5 OMEGA block on (CL, V1, Vmax, BASE, kdeg-kint) ----
    # Diagonal variances are derived from the reported CV% via omega^2 = log(CV^2 + 1):
    #   CV_CL        = 24.3% -> 0.05727
    #   CV_V1        = 23.4% -> 0.05328
    #   CV_Vmax      =  104% -> 0.73304
    #   CV_BASE      = 31.8% -> 0.09634
    #   CV_kdegkint  = 23.3% -> 0.05292
    # Off-diagonals are taken directly from Table 2 ("Covariance" section); the
    # paper reports off-diagonal elements of the OMEGA variance-covariance
    # matrix (covariances in log-space) using the symbol rho. The implied
    # correlation between eta_CL and eta_V1 is 0.0363 / sqrt(0.05727 * 0.05328)
    # = 0.658, consistent with typical mAb CL-V1 correlations.
    etalcl + etalvc + etalvmax + etalbase + etalkdegkint ~ c(
       0.05727,
       0.0363,    0.05328,
      -0.0122,    0.00352,  0.73304,
       0.0048,   -0.00284, -0.0265,   0.09634,
      -0.024,    -0.018,   -0.00303, -0.0216,   0.05292
    )

    # ---- Residual unexplained variability (Wojciechowski 2022 Table 2) ----
    # Free domagrozumab observations were modelled with additive error in the
    # natural-log domain; in nlmixr2 linear-domain notation this is
    # proportional with SD = sigma_add (= 0.142). Total myostatin used a
    # proportional error model in the linear domain; sigma_pro = 20.6% CV.
    CcpropSd   <- 0.142; label("Proportional residual error on free domagrozumab concentration (fraction)")  # Wojciechowski 2022 Table 2: sigma_add = 0.142 SD log-additive
    propSd_Myo <- 0.206; label("Proportional residual error on total myostatin concentration (fraction)")     # Wojciechowski 2022 Table 2: sigma_pro = 20.6% CV
  })

  model({
    # ---- Physical constants ----
    # Domagrozumab is a humanized IgG1 monoclonal antibody. The paper does not
    # report a molecular weight, so a representative IgG1 mass of 145,000 g/mol
    # (145 kDa) is used for the in-model unit bridge between mg-dosing and the
    # paper's nM concentration scale. This assumption is documented as a
    # deviation in the validation vignette.
    MW_DOMA <- 145000

    # ---- Individual PK parameters ----
    # Per-kg fixed effects scaled by WT to give per-subject CL/V/Q/Vmax. The
    # source paper does not list a population reference weight; weight enters
    # multiplicatively as if it were an allometric exponent of 1 (linear).
    cl    <- exp(lcl   + etalcl)   * WT          # L/hour
    vc    <- exp(lvc   + etalvc)   * WT          # L
    q     <- exp(lq)               * WT          # L/hour
    vp    <- exp(lvp)              * WT          # L
    ka    <- exp(lka)
    fdepot <- exp(lfdepot)
    vmax  <- exp(lvmax + etalvmax) * WT          # nM/hour
    km    <- exp(lkm)                            # nM
    kss   <- exp(lkss)                           # nM

    # ---- Myostatin parameters (with study-population covariate effect) ----
    # Eqs. 7-8: COVSPOP = 1 + theta_SPOP * DIS_DMD; multiplicative on the
    # typical value. Single eta_kdegkint drives both kdeg and kint, with kint
    # scaled by the structural ratio e_ratio_kdegkint.
    cov_spop_base     <- 1 + e_dmd_base     * DIS_DMD
    cov_spop_kdegkint <- 1 + e_dmd_kdegkint * DIS_DMD

    base <- exp(lbase + etalbase) * cov_spop_base                                # nM
    kdeg <- exp(lkdeg + etalkdegkint) * cov_spop_kdegkint                        # 1/hour
    kint <- exp(lkint + e_ratio_kdegkint * etalkdegkint) * cov_spop_kdegkint     # 1/hour
    ksyn <- base * kdeg                                                          # nM/hour (Eq. 5)

    # ---- Drug concentration in nM (paper scale) ----
    # central holds the drug amount in mg (rxode2 dosing convention); vc is in
    # L; so central/vc is in mg/L = ug/mL. Multiplying by 1e6/MW gives nM.
    Cc_nM <- (central / vc) * 1e6 / MW_DOMA

    # ---- Saturable elimination rate ----
    # Vmax in nM/hour, Cc_nM in nM, km in nM -> mm_rate_nM in nM/hour.
    # Convert to mg/hour for the amount-based ODE: nM * L = nmol; nmol * MW
    # (g/mol) * 1e-6 = mg.
    mm_rate_nM   <- vmax * Cc_nM / (km + Cc_nM)
    mm_rate_mghr <- mm_rate_nM * vc * MW_DOMA * 1e-6

    # ---- ODE system (Eqs. 1-6 of Wojciechowski 2022) ----
    # Drug PK on amounts (mg). Linear elimination via cl, MM via mm_rate_mghr,
    # and standard 2-compartment distribution. SC bioavailability is applied
    # via f(depot).
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot -
                           (cl / vc) * central -
                           mm_rate_mghr -
                           (q  / vc) * central +
                           (q  / vp) * peripheral1
    d/dt(peripheral1)  <-  (q  / vc) * central -
                           (q  / vp) * peripheral1
    f(depot) <- fdepot

    # Total myostatin (free + drug-myostatin complex) compartment (Eq. 4),
    # initial condition Eq. 6: total_target(0) = BASE.
    d/dt(total_target) <-  ksyn - kdeg * total_target -
                           (kint - kdeg) * Cc_nM * total_target / (kss + Cc_nM)
    total_target(0)    <-  base

    # ---- Observations ----
    Cc  <- Cc_nM           # Free domagrozumab concentration (nM)
    Myo <- total_target    # Total myostatin concentration (nM)

    Cc  ~ prop(CcpropSd)
    Myo ~ prop(propSd_Myo)
  })
}
