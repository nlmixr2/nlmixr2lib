Zuo_2016_UDCA <- function() {
  description <- "Systems model. Enterohepatic recirculation of ursodeoxycholic acid (UDCA) and its glycine (GUDCA) and taurine (TUDCA) conjugates in healthy adults, with adaptation to primary biliary cirrhosis (PBC). 19 ODEs across stomach, intestine, portal vein, blood, liver, biliary system, and feces compartments per analyte; oral square-wave absorption (0.5 h) and meal/snack-modulated biliary-to-intestinal flux. No IIV or residual error - typical-value mechanistic simulation only."
  reference <- "Zuo P, Dobbins RL, O'Connor-Semmes RL, Young MA. A Systems Model for Ursodeoxycholic Acid Metabolism in Healthy and Patients With Primary Biliary Cirrhosis. CPT Pharmacometrics Syst Pharmacol. 2016 Aug;5(8):418-426. doi:10.1002/psp4.12100"
  vignette <- "Zuo_2016_UDCA"
  paper_specific_compartments <- c("portal_udca", "biliary_udca", "feces_udca", "portal_gudca", "biliary_gudca", "feces_gudca", "portal_tudca", "biliary_tudca", "feces_tudca")

  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  covariateData <- list(
    FRACABS = list(
      description       = "Fractional absorption of UDCA from stomach into circulation. Dose-dependent per the paper's regression on log(dose_mg): F=0.66 for 150 mg and F=0.31 for 1000 mg (R^2=0.99 over 200-2000 mg from Crosignani 1991 / Walker 1992).",
      units             = "fraction",
      type              = "continuous",
      reference_category = NULL,
      notes             = "Applied as bioavailability on the stomach compartment, combined with mg->mmol unit conversion via the UDCA molecular weight. Only modulates first-pass entry; recirculated bile acids re-enter the intestine from bile without an additional F scaling (per paper Methods).",
      source_name       = "F"
    ),
    MEAL_FLAG = list(
      description       = "Indicator that the current time falls within a meal window (lunch / dinner, 1 hour duration each). Scales biliary-to-intestine rate constants (K_BI for UDCA, GUDCA, and TUDCA) by E_meal during this interval to mimic gallbladder contraction.",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (no meal effect)",
      notes             = "Time-varying covariate. Set to 1 over [t_meal, t_meal + 1 h] for each meal, 0 otherwise. Snacks are encoded separately via SNACK_FLAG. Two meals per dosing day were modelled (lunch at +4 h, dinner at +10 h after the morning dose, per Xiang 2011 and Dilger 2012).",
      source_name       = "meal indicator"
    ),
    SNACK_FLAG = list(
      description       = "Indicator that the current time falls within a snack window (0.5 hour duration). Scales biliary-to-intestine rate constants by E_snack during this interval.",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (no snack effect)",
      notes             = "Time-varying covariate. Set to 1 over [t_snack, t_snack + 0.5 h] for each snack, 0 otherwise. Two snacks per dosing day were modelled (at +7 h after the morning dose for the Xiang 2011 single-dose study).",
      source_name       = "snack indicator"
    ),
    DIS_PBC = list(
      description       = "Primary biliary cirrhosis disease state. Multiplies liver-to-biliary rate constants by 0.10 (UDCA), 0.30 (GUDCA), and 0.10 (TUDCA) when DIS_PBC=1 to reproduce the Zuo 2016 Figure 3 disease-state simulation.",
      units             = "(binary)",
      type              = "binary",
      reference_category = "0 (healthy)",
      notes             = "Time-fixed per-subject covariate. PBC scaling factors come from the paper's PBC adaptation (Methods, 'Model for PBC' and Figure 3 legend): K_LB,0 reduced 90%, K_LB,1 reduced 70%, K_LB,2 reduced 90%.",
      source_name       = "PBC indicator"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 49L,
    n_studies       = 3L,
    age_range       = "adults (cohort mean age not reported per subject; PBC cohort 30-65 years typical for the disease)",
    weight_range    = "not reported per subject; cohort means imputed for the Dilger 2012 chronic-dosing study (1000 mg/day used as the average for 15 mg/kg body weight)",
    sex_female_pct  = NA_real_,
    disease_state   = "Healthy adults (calibration: Xiang 2011 n=27 single-dose + Dilger 2012 n=11 chronic; validation: Hess 2004 n=21 chronic). PBC adults (Dilger 2012 n=11) for the disease-adaptation simulations.",
    dose_range      = "150 mg single oral tablet (Xiang 2011), 15 mg/kg body weight daily for 21 days (Dilger 2012; modelled as 1000 mg/day for the mean weight), 900 mg twice daily for 21 days (Hess 2004 validation).",
    regions         = "Multi-source: Finland (Xiang 2011), Germany (Dilger 2012), United States (Hess 2004).",
    notes           = "Mechanistic systems model fit to mean published profiles, not individual-level data; no IIV, no residual error. The PBC adaptation reduces K_LB rate constants per the paper's hypothesis-testing simulations. Sensitivity analysis (Table 3) drove the choice of which parameters to adjust between healthy and PBC states."
  )

  ini({
    # Rate constants for UDCA (analyte index 0) - Table 1 columns 'UDCA'.
    # All values are 1/hour. Log-transformed because the paper reports them
    # as positive point estimates from least-squares parameter optimisation
    # (MATLAB Optimization Toolbox).
    lk_si    <- log(16.61); label("Stomach -> intestine rate constant K_SI (1/h, UDCA)")            # Table 1
    lk_ip0   <- log(2.70);  label("Intestine -> portal rate constant K_IP,0 for UDCA (1/h)")        # Table 1
    lk_pb0   <- log(0.61);  label("Portal -> blood rate constant K_PB,0 for UDCA (1/h)")            # Table 1
    lk_bp0   <- log(6.94);  label("Blood -> portal rate constant K_BP,0 for UDCA (1/h)")            # Table 1
    lk_pl0   <- log(0.82);  label("Portal -> liver rate constant K_PL,0 for UDCA (1/h)")            # Table 1
    lk_lb0   <- log(0.33);  label("Liver -> biliary system rate constant K_LB,0 for UDCA (1/h)")    # Table 1
    lk_bi0   <- log(0.64);  label("Biliary system -> intestine rate constant K_BI,0 for UDCA (1/h)") # Table 1
    lk_if0   <- log(0.79);  label("Intestine -> feces rate constant K_IF,0 for UDCA (1/h)")         # Table 1

    # Rate constants for GUDCA (analyte index 1). K_IP,2 / K_PB,2 / K_BP,2 /
    # K_LB,2 / K_BI,2 / K_IF,2 are reported as 'same as' the GUDCA values in
    # Table 1 and therefore share these parameters in the structural model.
    lk_ip1   <- log(0.32);  label("Intestine -> portal rate constant K_IP,1 for GUDCA (shared with TUDCA, 1/h)")     # Table 1
    lk_pb1   <- log(0.36);  label("Portal -> blood rate constant K_PB,1 for GUDCA (shared with TUDCA, 1/h)")          # Table 1
    lk_bp1   <- log(0.66);  label("Blood -> portal rate constant K_BP,1 for GUDCA (shared with TUDCA, 1/h)")          # Table 1
    lk_pl1   <- log(1.68);  label("Portal -> liver rate constant K_PL,1 for GUDCA (1/h)")           # Table 1
    lk_lp1   <- log(0.10);  label("Liver -> portal rate constant K_LP,1 for GUDCA (1/h)")           # Table 1
    lk_lb1   <- log(0.32);  label("Liver -> biliary system rate constant K_LB,1 for GUDCA (shared with TUDCA, 1/h)")  # Table 1
    lk_bi1   <- log(0.13);  label("Biliary system -> intestine rate constant K_BI,1 for GUDCA (shared with TUDCA, 1/h)") # Table 1
    lk_if1   <- log(0.21);  label("Intestine -> feces rate constant K_IF,1 for GUDCA (shared with TUDCA, 1/h)")        # Table 1

    # TUDCA-specific rate constants (analyte index 2). K_PL,2 and K_LP,2 are
    # distinct from the GUDCA values because the paper introduces them to
    # describe TUDCA's different hepatocyte uptake affinity.
    lk_pl2   <- log(3.98);  label("Portal -> liver rate constant K_PL,2 for TUDCA (1/h)")           # Table 1
    lk_lp2   <- log(0.03);  label("Liver -> portal rate constant K_LP,2 for TUDCA (1/h; CI lower bound 0)") # Table 1

    # Liver conjugation / deconjugation
    lk_conj1 <- log(54.61); label("UDCA -> GUDCA conjugation rate constant K_CONJ,1 (1/h)")         # Table 1
    lk_decj1 <- log(0.76);  label("GUDCA -> UDCA deconjugation rate constant K_DECONJ,1 (1/h)")     # Table 1
    lk_conj2 <- log(3.83);  label("UDCA -> TUDCA conjugation rate constant K_CONJ,2 (1/h)")         # Table 1
    lk_decj2 <- log(0.17);  label("TUDCA -> UDCA deconjugation rate constant K_DECONJ,2 (1/h)")     # Table 1

    # Food effects on biliary-system -> intestine rate constants (gallbladder
    # contraction surrogate). Apply to all three analytes simultaneously.
    le_meal  <- log(35.33); label("Meal scaling factor E_meal on K_BI (1 h duration; unitless)")    # Table 1
    le_snack <- log(9.53);  label("Snack scaling factor E_snack on K_BI (0.5 h duration; unitless)") # Table 1
  })

  model({
    # === Physiological constants - fixed per paper Methods (citing Hofmann 1983). ===
    v_plasma <- 2.5    # Plasma volume of distribution (L) -- Methods, paragraph on parameter optimisation
    v_bile   <- 0.075  # Bile (biliary system) volume of distribution (L) -- same paragraph

    # === Molecular weights (g/mol; equivalent to mg/mmol) - fixed physical constants. ===
    mw_udca  <- 392.572  # 3-alpha,7-beta-dihydroxy-5-beta-cholanic acid
    mw_gudca <- 449.626  # UDCA + glycine - H2O
    mw_tudca <- 499.703  # UDCA + taurine - H2O

    # === PBC adaptation scaling factors on liver -> biliary rate constants. ===
    # Methods 'Model for PBC' and Figure 3 legend: K_LB,0 reduced 90%, K_LB,1
    # reduced 70%, K_LB,2 reduced 90% when adapting to the PBC population.
    pbc_klb0_fr <- 0.10
    pbc_klb1_fr <- 0.30
    pbc_klb2_fr <- 0.10

    # === Rate constants (back-transformed from log scale). ===
    k_si    <- exp(lk_si)
    k_ip0   <- exp(lk_ip0);  k_ip1 <- exp(lk_ip1); k_ip2 <- k_ip1
    k_pb0   <- exp(lk_pb0);  k_pb1 <- exp(lk_pb1); k_pb2 <- k_pb1
    k_bp0   <- exp(lk_bp0);  k_bp1 <- exp(lk_bp1); k_bp2 <- k_bp1
    k_pl0   <- exp(lk_pl0);  k_pl1 <- exp(lk_pl1); k_pl2 <- exp(lk_pl2)
    k_lp1   <- exp(lk_lp1);  k_lp2 <- exp(lk_lp2)
    # K_LP,0 is not reported in Table 1 (UDCA does not return liver -> portal).
    k_lb0   <- exp(lk_lb0);  k_lb1 <- exp(lk_lb1); k_lb2 <- k_lb1
    k_bi0   <- exp(lk_bi0);  k_bi1 <- exp(lk_bi1); k_bi2 <- k_bi1
    k_if0   <- exp(lk_if0);  k_if1 <- exp(lk_if1); k_if2 <- k_if1
    k_conj1 <- exp(lk_conj1); k_conj2 <- exp(lk_conj2)
    k_decj1 <- exp(lk_decj1); k_decj2 <- exp(lk_decj2)
    e_meal  <- exp(le_meal)
    e_snack <- exp(le_snack)

    # === Disease-state scaling on K_LB rate constants. ===
    # When PBC=1 each K_LB,i is multiplied by pbc_klb<i>_fr; when PBC=0 the
    # scaling resolves to 1 (healthy). Encoded with an indicator so a virtual
    # cohort with mixed PBC status can be simulated in a single rxSolve call.
    pbc_scale0 <- 1 - DIS_PBC * (1 - pbc_klb0_fr)
    pbc_scale1 <- 1 - DIS_PBC * (1 - pbc_klb1_fr)
    pbc_scale2 <- 1 - DIS_PBC * (1 - pbc_klb2_fr)
    k_lb0_eff  <- k_lb0 * pbc_scale0
    k_lb1_eff  <- k_lb1 * pbc_scale1
    k_lb2_eff  <- k_lb2 * pbc_scale2

    # === Meal / snack scaling on biliary-system -> intestine rate constants. ===
    # food_mult = 1 outside meals and snacks, E_meal during a meal window,
    # E_snack during a snack window. The user supplies MEAL_FLAG and SNACK_FLAG
    # as time-varying covariates with the meal / snack on-times set to 1.
    food_mult <- 1 + (e_meal - 1) * MEAL_FLAG + (e_snack - 1) * SNACK_FLAG
    k_bi0_eff <- k_bi0 * food_mult
    k_bi1_eff <- k_bi1 * food_mult
    k_bi2_eff <- k_bi2 * food_mult

    # === ODE system - amounts in mmol. ===
    # Stomach receives the dose; zero-order release into intestine over a time
    # scale of 1 / K_SI (~3.6 minutes here), approximating the 0.5 h dissolution
    # square wave of the source paper (Methods).
    d/dt(stomach_udca)   <- -k_si * stomach_udca
    d/dt(intestine_udca) <-  k_si * stomach_udca +
                             k_bi0_eff * biliary_udca -
                             k_ip0 * intestine_udca -
                             k_if0 * intestine_udca
    d/dt(portal_udca)    <-  k_ip0 * intestine_udca +
                             k_bp0 * blood_udca -
                             k_pb0 * portal_udca -
                             k_pl0 * portal_udca
    d/dt(blood_udca)     <-  k_pb0 * portal_udca -
                             k_bp0 * blood_udca
    d/dt(liver_udca)     <-  k_pl0 * portal_udca +
                             k_decj1 * liver_gudca +
                             k_decj2 * liver_tudca -
                             k_conj1 * liver_udca -
                             k_conj2 * liver_udca -
                             k_lb0_eff * liver_udca
    d/dt(biliary_udca)   <-  k_lb0_eff * liver_udca -
                             k_bi0_eff * biliary_udca
    d/dt(feces_udca)     <-  k_if0 * intestine_udca

    # GUDCA - 6 compartments, no stomach (enters the system only via hepatic
    # conjugation of UDCA in liver_udca).
    d/dt(intestine_gudca) <- k_bi1_eff * biliary_gudca -
                             k_ip1 * intestine_gudca -
                             k_if1 * intestine_gudca
    d/dt(portal_gudca)    <- k_ip1 * intestine_gudca +
                             k_bp1 * blood_gudca +
                             k_lp1 * liver_gudca -
                             k_pb1 * portal_gudca -
                             k_pl1 * portal_gudca
    d/dt(blood_gudca)     <- k_pb1 * portal_gudca -
                             k_bp1 * blood_gudca
    d/dt(liver_gudca)     <- k_pl1 * portal_gudca +
                             k_conj1 * liver_udca -
                             k_lp1 * liver_gudca -
                             k_decj1 * liver_gudca -
                             k_lb1_eff * liver_gudca
    d/dt(biliary_gudca)   <- k_lb1_eff * liver_gudca -
                             k_bi1_eff * biliary_gudca
    d/dt(feces_gudca)     <- k_if1 * intestine_gudca

    # TUDCA - 6 compartments, same topology as GUDCA.
    d/dt(intestine_tudca) <- k_bi2_eff * biliary_tudca -
                             k_ip2 * intestine_tudca -
                             k_if2 * intestine_tudca
    d/dt(portal_tudca)    <- k_ip2 * intestine_tudca +
                             k_bp2 * blood_tudca +
                             k_lp2 * liver_tudca -
                             k_pb2 * portal_tudca -
                             k_pl2 * portal_tudca
    d/dt(blood_tudca)     <- k_pb2 * portal_tudca -
                             k_bp2 * blood_tudca
    d/dt(liver_tudca)     <- k_pl2 * portal_tudca +
                             k_conj2 * liver_udca -
                             k_lp2 * liver_tudca -
                             k_decj2 * liver_tudca -
                             k_lb2_eff * liver_tudca
    d/dt(biliary_tudca)   <- k_lb2_eff * liver_tudca -
                             k_bi2_eff * biliary_tudca
    d/dt(feces_tudca)     <- k_if2 * intestine_tudca

    # === Dose handling on stomach compartment. ===
    # Users dose UDCA tablets in mg. The bioavailability hook converts
    # mg -> mmol (divide by mw_udca, in g/mol = mg/mmol) and applies the
    # dose-dependent fractional absorption FRACABS (0.66 at 150 mg, 0.31 at
    # 1000 mg; see covariateData). Combined with dur = 0.5 h on the dose row,
    # this reproduces the 0.5 h square-wave dissolution of the source paper.
    f(stomach_udca) <- FRACABS / mw_udca

    # === Observation variables. ===
    # Plasma concentrations in umol/L. State amounts are in mmol and volumes
    # in L, so the amount/volume ratio is in mmol/L; multiply by 1000 to land
    # in umol/L (= uM), matching the reporting scale of Xiang 2011, Dilger
    # 2012, and Hess 2004.
    Cc_udca   <- blood_udca   * 1000 / v_plasma
    Cc_gudca  <- blood_gudca  * 1000 / v_plasma
    Cc_tudca  <- blood_tudca  * 1000 / v_plasma

    # Biliary (duodenal-bile aspirate) concentrations in umol/L. Used to
    # compare against Dilger 2012's baseline / peak / trough biliary
    # measurements.
    Cb_udca   <- biliary_udca  * 1000 / v_bile
    Cb_gudca  <- biliary_gudca * 1000 / v_bile
    Cb_tudca  <- biliary_tudca * 1000 / v_bile
  })
}
