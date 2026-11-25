Hasegawa_2018_denosumab_qsp <- function() {
  description <- paste(
    "QSP (reduced). Eight-state nonlinear bone-biology model obtained by inductive linearization and",
    "proper lumping of the 28-state Peterson and Riggs calcium-homeostasis / bone-remodelling model,",
    "driven by a two-compartment quasi-steady-state TMDD model of subcutaneous denosumab, with lumbar",
    "spine BMD (percent of baseline) as an indirect response to osteoblasts and osteoclasts. Lumped",
    "states: lumped1 = RANK-related state (RANK plus the original states pooled with it); lumped2 =",
    "RANKL-related state (RANKL, parathyroid max capacity and the OPG-RANKL complex). Unlumped states:",
    "responding, fast and slow osteoblasts, osteoclasts, active TGF-beta and the RANK-RANKL complex.",
    "Four parameters (osteoclast degradation, lumped RANKL production, the bone-side KSS and the BMD",
    "kout) were estimated from the first 12 months of mean BMD data in postmenopausal women; the",
    "rest are fixed. Typical-value (mean-data) model with an additive residual on BMD."
  )
  reference <- paste(
    "Hasegawa C, Duffull SB. Automated Scale Reduction of Nonlinear QSP Models With an Illustrative",
    "Application to a Bone Biology System. CPT Pharmacometrics Syst Pharmacol. 2018;7(9):562-572.",
    "doi:10.1002/psp4.12324. The authors' NONMEM control stream for the final nonlinear lumped model",
    "(Appendix S3, nonlinear_lumped.mod) and their MATLAB lumping code (Appendix S2) were used for",
    "the structure and fixed constants."
  )
  vignette <- "Hasegawa_2018_denosumab_bone_scale_reduction"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "denosumab", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "denosumab (total, free + bound)", units = "ug", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "denosumab", units = "ug", specimen = "tissue", verified = TRUE),
    total_target = list(
      analyte = "total RANKL (free + denosumab-bound)",
      units = "ng/mL",
      specimen = "serum",
      verified = TRUE
    ),
    ob_fast = list(
      analyte = "fast-turnover osteoblasts (FOB)",
      units = "model amount",
      specimen = "tissue",
      verified = TRUE
    ),
    ob_slow = list(
      analyte = "slow-turnover osteoblasts (SOB)",
      units = "model amount",
      specimen = "tissue",
      verified = TRUE
    ),
    oc = list(analyte = "osteoclasts", units = "model amount", specimen = "tissue", verified = TRUE),
    rob = list(analyte = "responding osteoblasts", units = "model amount", specimen = "tissue", verified = TRUE),
    lumped2 = list(
      analyte = "RANKL-related lumped state (L2)",
      units = "model amount",
      specimen = "tissue",
      verified = TRUE
    ),
    tgfb_act = list(analyte = "active TGF-beta", units = "model amount", specimen = "tissue", verified = TRUE),
    rank_rankl = list(
      analyte = "RANK-RANKL complex (CMX)",
      units = "model amount",
      specimen = "tissue",
      verified = TRUE
    ),
    lumped1 = list(
      analyte = "RANK-related lumped state (L1)",
      units = "model amount",
      specimen = "tissue",
      verified = TRUE
    ),
    BMD_LS = list(
      analyte = "lumbar spine bone mineral density",
      units = "% of baseline (100 at time 0)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # Unlumped bone-model states of the reduced model (Hasegawa 2018 Figure 4);
  # paper-mechanistic, amounts on the arbitrary scale of the Peterson and
  # Riggs model.
  paper_specific_compartments <- c("ob_fast", "ob_slow", "oc", "rob", "tgfb_act", "rank_rankl")

  population <- list(
    species = "human",
    n_subjects = 190L,
    n_studies = 1L,
    sex_female_pct = 100,
    weight_median = "66 kg (fixed in the denosumab absorption rate)",
    disease_state = "Postmenopausal women with low bone mineral density",
    dose_range = "Denosumab 14, 60, 100 or 210 mg subcutaneously every 6 months",
    regions = "United States (phase 2 denosumab dose-ranging study)",
    notes = paste(
      "Fitted to digitized arm-mean lumbar-spine BMD (percent change from baseline) from the phase 2",
      "denosumab dose-ranging study (McClung 2006; Miller 2008 for the 48-month extension), every-6-month",
      "arms of 54, 47, 42 and 47 women at 14, 60, 100 and 210 mg (Hasegawa 2018 Results; data in",
      "Table S2). Only months 0-12 were used for estimation; months 12-48 were predicted. The denosumab",
      "PK and the bone-model constants were not re-estimated: PK is the Sutjandra 2011 population model",
      "as coded by the authors, and the bone constants come from the lumping formulae (Table 1)."
    )
  )

  ini({
    # ---- Denosumab PK (fixed; Sutjandra 2011 as coded by the authors) ----
    # NONMEM stream $THETA 22-27 (per-hour values; the stream multiplies by
    # 24*30 to work in 30-day months) and MATLAB Model_parameter_values0.m.
    lcl <- fixed(log(0.00306))
    label("Linear clearance of denosumab (L/h)") # nonlinear_lumped.mod THETA(22) = 3.06 FIX (mL/h)
    lvc <- fixed(log(2.49))
    label("Central volume of distribution (L)") # nonlinear_lumped.mod THETA(23) = 2490 FIX (mL)
    lq <- fixed(log(0.0379))
    label("Intercompartmental clearance (L/h)") # nonlinear_lumped.mod THETA(24) = 37.9 FIX (mL/h)
    lvp <- fixed(log(1.36))
    label("Peripheral volume of distribution (L)") # nonlinear_lumped.mod THETA(25) = 1360 FIX (mL)
    lka <- fixed(log(0.212 / 24))
    label("Subcutaneous absorption rate constant at 71.5 kg (1/h)") # nonlinear_lumped.mod THETA(26) = 0.212 FIX (1/day)
    e_wt_ka <- fixed(-0.577)
    label("Power exponent of (WT/71.5) on ka, evaluated at a body weight of 66 kg (unitless)") # nonlinear_lumped.mod: KA = TVKA*(66/71.5)**(-0.577)
    lfdepot <- fixed(log(0.638))
    label("Subcutaneous bioavailability (fraction)") # nonlinear_lumped.mod THETA(27) = 0.638 FIX
    lkint <- fixed(log(0.00795))
    label("Internalization rate constant of the denosumab-RANKL complex (1/h)") # Table 1 kint = 0.00795; nonlinear_lumped.mod $DES KINT = 0.00795
    lkss <- fixed(log(1380))
    label("Quasi-steady-state constant of the denosumab PK model (ng/mL)") # MATLAB Model_parameter_values0.m Kss = 138*10; NONMEM read it from data column KSSINI (see vignette)
    lrbase_target <- fixed(log(614))
    label("Baseline total RANKL in the PK model (ng/mL)") # MATLAB Model_parameter_values0.m R0 = 614; NONMEM data column R00
    lkdeg <- fixed(log(0.00148))
    label("Degradation rate constant of free RANKL in the PK model (1/h)") # MATLAB Model_parameter_values0.m kdeg = 0.00148; NONMEM data column KDEG

    # ---- Initial conditions of the reduced bone model (Eqs 13-20) ----
    ob_fast0 <- fixed(0.0040)
    label("Initial fast-osteoblast amount FOB(t0)") # Eq 19 FOB(t0) = 0.0040; nonlinear_lumped.mod IC1
    ob_slow0 <- fixed(0.0010)
    label("Initial slow-osteoblast amount SOB(t0)") # Eq 20 SOB(t0) = 0.0010; nonlinear_lumped.mod IC2
    oc0 <- fixed(0.0012)
    label("Initial osteoclast amount OC(t0)") # Eq 16 OC(t0) = 0.0012; nonlinear_lumped.mod IC3
    rob0 <- fixed(0.0010)
    label("Initial responding-osteoblast amount ROB(t0)") # Eq 18 ROB(t0) = 0.0010; nonlinear_lumped.mod IC4
    lumped20 <- fixed(1.4001)
    label("Initial RANKL-related lumped state L2(t0)") # Eq 14 L2(t0) = 1.4001; nonlinear_lumped.mod IC5
    tgfb_act0 <- fixed(0.2281)
    label("Initial active TGF-beta TGF(t0)") # Eq 17 TGF(t0) = 0.2281; nonlinear_lumped.mod IC6
    rank_rankl0 <- fixed(0.00022283)
    label("Initial RANK-RANKL complex CMX(t0)") # Eq 15 CMX(t0) = 0.00022283; nonlinear_lumped.mod IC7
    lumped10 <- fixed(30082)
    label("Initial RANK-related lumped state L1(t0)") # Eq 13 L1(t0) = 30,082; nonlinear_lumped.mod IC8

    # ---- RANK-related lumped state L1 (Eq 21) ----
    lr_lumped1 <- fixed(log(75.0))
    label("Production rate of L1, R_L1 (amount/h)") # Table 1 R_L1 = 75.0; nonlinear_lumped.mod THETA(13) = 75 FIX
    lk_ob_lumped1 <- fixed(log(55.3))
    label("Rate constant of OB-driven production of L1, k_OB->L1 (1/h)") # Table 1 k_OB->L1 = 55.3; THETA(14) = 55.3 FIX
    lk_lumped2_lumped1 <- fixed(log(160))
    label("Rate constant from L2 to L1, k_L2->L1 (1/h)") # Table 1 k_L2->L1 = 160; THETA(15) = 160 FIX
    lk_cmx_l <- fixed(log(0.112))
    label("Rate constant from CMX to L1 and to L2, and CMX loss, k_CMX->L (1/h)") # Table 1 k_CMX->L = 0.112; THETA(10) = .1120 FIX
    ld_lumped1 <- fixed(log(0.970))
    label("Degradation rate constant of L1, d_L1 (1/h)") # Table 1 d_L1 = 0.970; THETA(16) = .97 FIX

    # ---- RANKL-related lumped state L2 (Eq 22) ----
    lr_lumped2 <- log(0.00337)
    label("Production rate of L2, R_L2 (amount/h)") # Table 1 R_L2 after estimation = 0.00337 (9.1% RSE)
    lk_ob_rankl <- fixed(log(0.234))
    label("Rate constant of OB-driven production of RANKL, k_OB->RANKL (1/h)") # Table 1 k_OB->RANKL = 0.234; THETA(17) = .234 FIX
    ld_lumped2 <- fixed(log(0.00110))
    label("Degradation rate constant of L2, d_L2 (1/h)") # Table 1 d_L2 = 0.00110; THETA(9) = .0011 FIX
    ld_rankl <- fixed(log(0.00290))
    label("Degradation rate constant of RANKL, d_RANKL (1/h)") # Table 1 d_RANKL = 0.00290; nonlinear_lumped.mod $DES KOUTL = 0.00290
    lkss_rankl <- log(63.4)
    label("Steady-state constant of denosumab-RANKL binding in the bone model, K_SS (ng/mL)") # Table 1 K_SS after estimation = 63.4 (63.7% RSE)

    # ---- RANK-RANKL complex CMX (Eq 23) ----
    lk_lumped2_cmx <- fixed(log(0.0000190))
    label("Rate constant from L2 to CMX, k_L2->CMX (1/h)") # Table 1 k_L2->CMX = 0.0000190; THETA(12) = .000019 FIX

    # ---- Osteoclasts (Eq 24) ----
    lr_oc <- fixed(log(0.00000298))
    label("Production rate of OC, R_OC (amount/h)") # Table 1 R_OC = 0.00000298; THETA(6) = .00298 FIX (x 1/1000)
    ld_oc <- log(0.0898)
    label("Degradation rate constant of OC, d_OC (1/h)") # Table 1 d_OC after estimation = 0.0898 (5.4% RSE)
    frho_tgf_oc <- fixed(0.878215)
    label("Minimum TGF effect on OC degradation as a fraction of TGF(t0), rho1/TGF(t0) (unitless)") # nonlinear_lumped.mod E0PICOC = 0.878215*PIC0 (Table 1 rho1 = 0.200)
    emax_tgf_oc <- fixed(2)
    label("Span of the TGF effect on OC degradation, a1 - rho1 (unitless)") # nonlinear_lumped.mod EMAXPICOC = 2 (Table 1 a1 = 2.18; see vignette)
    gam_tgf_oc <- fixed(1)
    label("Sigmoidicity of the TGF effect on OC degradation, gamma1 (unitless)") # Table 1 gamma1 = 1; stream PICOC has no exponent
    a_cmx_oc <- fixed(3.8034)
    label("Maximum CMX effect on OC degradation, a2 (unitless)") # Table 1 a2 = 3.80; nonlinear_lumped.mod E0RANKL = 3.8034
    rho_cmx_oc <- fixed(0.4698)
    label("Minimum CMX effect on OC degradation, rho2 (unitless)") # Table 1 rho2 = 0.470; nonlinear_lumped.mod EMAXL = 0.4698
    gam_cmx_oc <- fixed(3.0902)
    label("Sigmoidicity of the CMX effect on OC degradation, gamma2 (unitless)") # Table 1 gamma2 = 3.09; nonlinear_lumped.mod LSURVOCGAM = 3.0902

    # ---- Active TGF-beta (Eq 25) ----
    ld_tgf <- fixed(log(0.0298))
    label("Degradation rate constant of active TGF-beta, d_TGF (1/h)") # Table 1 d_TGF = 0.0298; THETA(11) = .0298 FIX

    # ---- Osteoblasts (Eqs 26-28) ----
    lr_rob <- fixed(log(0.000003))
    label("Production rate of ROB, R_ROB (amount/h)") # Table 1 R_ROB = 0.000003; THETA(1) = .003 FIX (x 1/1000)
    frho_tgf_rob <- fixed(0.8838)
    label("Minimum TGF effect on ROB production as a fraction of TGF(t0), rho3/TGF(t0) (unitless)") # nonlinear_lumped.mod E0PICROB = 0.8838*PIC0 (Table 1 rho3 = 0.202)
    emax_tgf_rob <- fixed(3.9745)
    label("Span of the TGF effect on ROB production, a3 - rho3 (unitless)") # nonlinear_lumped.mod EMAXPICROB = 3.9745 (Table 1 a3 = 4.18)
    gam_tgf_rob <- fixed(1.810)
    label("Sigmoidicity of the TGF effect on ROB production, gamma3 (unitless)") # Table 1 gamma3 = 1.81; nonlinear_lumped.mod GAMMA = 1.810
    lk_rob_ob <- fixed(log(0.003))
    label("Rate constant from ROB to OB, k_ROB->OB (1/h)") # Table 1 k_ROB->OB = 0.003; THETA(2) = .003 FIX
    f_sob <- fixed(0.06)
    label("Fraction of ROB converting to slow osteoblasts, f (fraction)") # Table 1 f = 0.06; THETA(3) = .06 FIX
    ld_fob <- fixed(log(0.01))
    label("Degradation rate constant of FOB, d_FOB (1/h)") # Table 1 d_FOB = 0.01; THETA(4) = .01 FIX
    ld_sob <- fixed(log(0.000001))
    label("Degradation rate constant of SOB, d_SOB (1/h)") # Table 1 d_SOB = 0.000001; THETA(5) = .001 FIX (x 1/1000)

    # ---- Lumbar spine BMD (Eq 4) ----
    lkout <- log(0.000533)
    label("First-order loss rate constant of BMD, k_out (1/h)") # Results: 'The estimate of kout in the reduced model was 0.000533/h (24.6%RSE)'
    gam_ob_bmd <- fixed(0.0739)
    label("Power of (OB/OB0) on BMD formation, gamma_OB (unitless)") # nonlinear_lumped.mod THETA(19) POWOB = .0739 FIX
    gam_oc_bmd <- fixed(0.0779)
    label("Power of (OC/OC0) on BMD loss, gamma_OC (unitless)") # nonlinear_lumped.mod THETA(20) POWOC = .0779 FIX

    # ---- Residual error ----
    addSd <- 0.6180615
    label("Additive residual error on BMD (% of baseline)") # Table 1 sigma^2 = 0.382 (24.0% RSE); sqrt(0.382)
  })

  model({
    # 1. Denosumab PK: two-compartment QSS TMDD (Appendix S1; stream $DES).
    #    The stream fixes body weight at 66 kg in the ka allometry.
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka) * (66 / 71.5)^e_wt_ka
    kint <- exp(lkint)
    kss <- exp(lkss)
    rbase_target <- exp(lrbase_target)
    kdeg <- exp(lkdeg)
    ksyn <- rbase_target * kdeg

    # Free serum denosumab (ng/mL): positive root of the QSS quadratic.
    ctot <- central / vc
    qssdisc <- ctot - total_target - kss
    Cc <- 0.5 * (qssdisc + sqrt(qssdisc * qssdisc + 4 * kss * ctot))
    cltot <- cl + kint * vc * total_target / (kss + Cc)

    # 2. Bone-model rate constants (Table 1; per hour).
    r_lumped1 <- exp(lr_lumped1)
    k_ob_lumped1 <- exp(lk_ob_lumped1)
    k_lumped2_lumped1 <- exp(lk_lumped2_lumped1)
    k_cmx_l <- exp(lk_cmx_l)
    d_lumped1 <- exp(ld_lumped1)
    r_lumped2 <- exp(lr_lumped2)
    k_ob_rankl <- exp(lk_ob_rankl)
    d_lumped2 <- exp(ld_lumped2)
    d_rankl <- exp(ld_rankl)
    kss_rankl <- exp(lkss_rankl)
    k_lumped2_cmx <- exp(lk_lumped2_cmx)
    r_oc <- exp(lr_oc)
    d_oc <- exp(ld_oc)
    d_tgf <- exp(ld_tgf)
    r_rob <- exp(lr_rob)
    k_rob_ob <- exp(lk_rob_ob)
    d_fob <- exp(ld_fob)
    d_sob <- exp(ld_sob)
    kout <- exp(lkout)

    # TGF production is tied to the initial TGF/OC ratio (stream K36 = K60*IC6/IC3;
    # Table 1 k_OC->TGF = 5.66).
    k_oc_tgf <- d_tgf * tgfb_act0 / oc0

    # 3. Hill-function constants, derived exactly as the stream derives them.
    #    TGF effect on ROB production (Eq 26): rho3, a3 = rho3 + emax, delta3.
    rho_tgf_rob <- frho_tgf_rob * tgfb_act0
    delta_tgf_rob <- emax_tgf_rob * tgfb_act0 / (tgfb_act0 - rho_tgf_rob) - tgfb_act0
    #    TGF effect on OC degradation (Eq 24): rho1, a1 = rho1 + emax, delta1.
    rho_tgf_oc <- frho_tgf_oc * tgfb_act0
    delta_tgf_oc <- emax_tgf_oc * tgfb_act0 / (tgfb_act0 - rho_tgf_oc) - tgfb_act0
    #    CMX effect on OC degradation (Eq 24): driver CMX/10, delta2 chosen so
    #    the effect equals 1 at CMX(t0).
    pil0 <- rank_rankl0 / 10
    delta_cmx_oc <- ((a_cmx_oc - rho_cmx_oc) * (pil0^gam_cmx_oc / (a_cmx_oc - 1)) -
      pil0^gam_cmx_oc)^(1 / gam_cmx_oc)

    eff_tgf_rob <- rho_tgf_rob + emax_tgf_rob * tgfb_act^gam_tgf_rob /
      (tgfb_act^gam_tgf_rob + delta_tgf_rob^gam_tgf_rob)
    eff_tgf_oc <- rho_tgf_oc + emax_tgf_oc * tgfb_act^gam_tgf_oc /
      (tgfb_act^gam_tgf_oc + delta_tgf_oc^gam_tgf_oc)
    pil <- rank_rankl / 10
    eff_cmx_oc <- a_cmx_oc - (a_cmx_oc - rho_cmx_oc) * pil^gam_cmx_oc /
      (pil^gam_cmx_oc + delta_cmx_oc^gam_cmx_oc)

    # Denosumab removal of the RANKL-related state (Eq 22; stream DRUG/3).
    drug_rankl <- (kint - d_rankl) * Cc / (kss_rankl + Cc)

    # 4. ODEs.
    d/dt(depot) <- -ka * depot
    f(depot) <- exp(lfdepot)
    d/dt(central) <- ka * depot - cltot * Cc - q * (Cc - peripheral1 / vp)
    d/dt(peripheral1) <- q * (Cc - peripheral1 / vp)
    d/dt(total_target) <- ksyn - (kdeg + (kint - kdeg) * Cc / (kss + Cc)) * total_target

    d/dt(ob_fast) <- k_rob_ob * (1 - f_sob) * rob - d_fob * ob_fast # Eq 27
    d/dt(ob_slow) <- k_rob_ob * f_sob * rob - d_sob * ob_slow # Eq 28
    d/dt(oc) <- r_oc - d_oc * eff_tgf_oc * eff_cmx_oc * oc # Eq 24
    d/dt(rob) <- r_rob * eff_tgf_rob - k_rob_ob * rob # Eq 26
    d/dt(lumped2) <- r_lumped2 + k_ob_rankl * (ob_fast + ob_slow) + k_cmx_l * rank_rankl -
      (d_lumped2 + drug_rankl / 3) * lumped2 # Eq 22
    d/dt(tgfb_act) <- k_oc_tgf * oc - d_tgf * tgfb_act # Eq 25
    d/dt(rank_rankl) <- k_lumped2_cmx * lumped2 - k_cmx_l * rank_rankl # Eq 23
    d/dt(lumped1) <- r_lumped1 + k_ob_lumped1 * (ob_fast + ob_slow) +
      k_lumped2_lumped1 * lumped2 + k_cmx_l * rank_rankl - d_lumped1 * lumped1 # Eq 21

    # Eq 4 with R_in = k_out * 100 (stream RIN_BMD and K_BMD).
    d/dt(BMD_LS) <- kout * 100 * ((ob_fast + ob_slow) / (ob_fast0 + ob_slow0))^gam_ob_bmd -
      kout * (oc / oc0)^gam_oc_bmd * BMD_LS

    total_target(0) <- rbase_target
    ob_fast(0) <- ob_fast0
    ob_slow(0) <- ob_slow0
    oc(0) <- oc0
    rob(0) <- rob0
    lumped2(0) <- lumped20
    tgfb_act(0) <- tgfb_act0
    rank_rankl(0) <- rank_rankl0
    lumped1(0) <- lumped10
    BMD_LS(0) <- 100

    # 5. Observation: BMD as percent of baseline; the fitted quantity was the
    #    percent change from baseline, BMD_LS - 100, with the same additive error.
    bmd_pcfb <- BMD_LS - 100
    BMD_LS ~ add(addSd)
  })
}
