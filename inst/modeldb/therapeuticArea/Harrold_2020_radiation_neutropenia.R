Harrold_2020_radiation_neutropenia <- function() {
  description <- "Preclinical (rhesus macaque). Semi-mechanistic ANC response model linking acute radiation injury to neutropenia, with a coupled time-to-event survival sub-model linking the ANC time course to overall survival in nonhuman primates exposed to 750 cGy whole-body irradiation on day 0 and treated with placebo or G-CSF (filgrastim or pegfilgrastim). Three components: (a) a K-PD radiation injury compartment (depot_kpd) that decays exponentially at k_PD,e and drives mitotic-cell killing via a power-law sensitivity coefficient k_kill = k_PD,kill * RAD^gamma (Eq. 2); (b) a Friberg-style granulopoiesis chain (precursor1 = N_SM stem cell -> precursor2 = N_MT mitotic -> precursor3 = N_PM1 -> precursor4 = N_PM2 -> circ = ANC) where the radiation kill effect operates on the mitotic compartment via the additive loss term -k_kill * precursor2 (Eqs. 3-7); and (c) a Box-Cox-transformed-ANCe-driven time-varying hazard model for overall survival (Eqs. 9-11), with study-specific parameter sets (Table III) selected by the binary indicator STUDY_HARROLD_PEG (0 = filgrastim pivotal study, 1 = pegfilgrastim pivotal study). The ANC response parameters (Table II) were jointly fit on the combined placebo cohorts (n = 45 NHPs); the OS sub-model parameters (Table III) were fit separately to each study because an unexplained study effect remained after covariate screening (Discussion p. 7)."
  reference <- paste(
    "Harrold JM, Olsson Gisleskog P, Delor I, Jacqmin P, Perez-Ruixo JJ, Narayanan A, Doshi S, Chow A, Yang B-B, Melhem M.",
    "Quantification of Radiation Injury on Neutropenia and the Link between Absolute Neutrophil Count Time Course",
    "and Overall Survival in Nonhuman Primates Treated with G-CSF.",
    "Pharm Res. 2020;37(7):102.",
    "doi:10.1007/s11095-020-02839-3 (PMID: 32440783; PMC: PMC7242243).",
    sep = " "
  )
  vignette <- "Harrold_2020_radiation_neutropenia"
  units <- list(
    time          = "day",
    dosing        = "K-PD radiation amount, numerical magnitude in cGy (the published 750 cGy pivotal dose enters the depot_kpd compartment as a 750-unit bolus at t = 0; Table II reports k_PD,kill = 2.14 d^-1 'KPD^-1' where KPD = cGy^gamma. Simulating with `amt = 7.5` (Gy units) produces an ANC nadir that is too shallow and too early to match the paper's Fig 2 / Fig 3 placebo VPC; `amt = 750` reproduces the deep nadir at ~14-15 days reported in the Discussion. See vignette Assumptions and deviations.)",
    concentration = "10^9 cells/L (ANC)"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric-like power effect on the maturation rate constant k_tr (Eq. 14 form), with reference body weight 5.9 kg (combined-studies median baseline weight from Harrold 2020 Table I). Fixed exponent e_wt_ktr = 0.629 (Table II: 'gamma_wt' Fixed).",
      source_name        = "BWT"
    ),
    STUDY_HARROLD_PEG = list(
      description        = "Binary indicator selecting the pegfilgrastim pivotal study OS parameter set",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (filgrastim pivotal study)",
      notes              = "1 = subject from the pegfilgrastim pivotal study (reference 13 in Harrold 2020; pegfilgrastim 300 ug/kg on days 1 and 8, SC); 0 = subject from the filgrastim pivotal study (reference 14 in Harrold 2020; filgrastim 10 ug/kg QD SC starting on day 1). Used only by the OS time-to-event sub-model (Eqs. 9-11) to select between the two parameter sets in Table III (lambda_ANC, lambda_BC, k_e0). The ANC response sub-model (Table II) does not depend on this indicator -- those parameters were fit jointly on the combined placebo cohorts of both studies.",
      source_name        = "STUDY"
    )
  )

  population <- list(
    species        = "nonhuman primate (rhesus macaque, Macaca mulatta)",
    n_subjects     = 92L,
    n_studies      = 2L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    weight_median  = "5.9 kg (combined studies; SD 0.86)",
    sex_female_pct = 8.7,
    disease_state  = "Hematopoietic syndrome of acute radiation syndrome (HS-ARS) following 750 cGy whole-body irradiation at 80 cGy/min on day 0.",
    dose_range     = "Whole-body irradiation 750 cGy (7.5 Gy) at 80 cGy/min on day 0 (administered to the radiation K-PD compartment as a 1 KPD-unit bolus at t = 0); placebo or filgrastim 10 ug/kg QD SC starting day 1 (until ANC recovered >= 1 x 10^9 cells/L for 3 consecutive days; median dosing days 1-19) or pegfilgrastim 300 ug/kg SC on days 1 and 8. Drug PK is not modelled here; the drug-treated cohorts are modelled implicitly through the OS sub-model's per-study parameter set.",
    n_observations = 1346L,
    follow_up      = "60 days; ANC measured every 1-2 days",
    notes          = "Two pivotal NHP studies: filgrastim pivotal (n = 22 placebo + n = 24 drug = 46; reference 14) and pegfilgrastim pivotal (n = 23 placebo + n = 23 drug = 46; reference 13). Median baseline ANC differs between studies (filgrastim 4.78, pegfilgrastim 1.77 x 10^9 cells/L; Table I and Results). Only the filgrastim pivotal study had female subjects (n = 4 per cohort). The ANC response model was fit to the combined 1346 placebo measurements; the OS model was fit separately per study (5000 bootstrap draws for the OS model)."
  )

  ini({
    # ----------------------------------------------------------------------
    # ANC response model (Table II of Harrold 2020).
    # Parameters were jointly fit on the combined placebo cohorts (n = 45,
    # 1346 ANC measurements). All structural parameters are log-transformed
    # in line with the paper's log-normal IIV (Eq. 12) and the standard
    # nlmixr2lib convention for strictly-positive rate constants and
    # baselines.
    # ----------------------------------------------------------------------
    lrbase    <- log(1.71)    ; label("Baseline circulating ANC at homeostasis (10^9 cells/L)")              # Table II: ANC_IC = 1.71 (95% CI 1.32, 2.10; RSE 11.6%)
    lktr      <- log(0.628)   ; label("Precursor maturation / transit rate k_tr (1/day)")                    # Table II: k_tr = 0.628 (95% CI 0.608, 0.648; RSE 1.64%)
    lkout     <- log(1.34)    ; label("Circulating ANC turnover rate k_c (1/day)")                           # Table II: k_c  = 1.34  (95% CI 0.892, 1.79; RSE 16.9%)
    lkel      <- log(0.312)   ; label("K-PD radiation elimination rate k_PD,e (1/day)")                      # Table II: k_PD,e = 0.312 (95% CI 0.296, 0.329; RSE 2.69%)
    lkdecay   <- log(2.14)    ; label("K-PD radiation kill coefficient k_PD,kill (1/(day * KPD^gamma))")     # Table II: k_PD,kill = 2.14 (95% CI 1.43, 2.86; RSE 19.2%)
    lgamma    <- log(1.79)    ; label("Radiation sensitivity exponent gamma (unitless; on RAD^gamma)")       # Table II: gamma = 1.79 (95% CI 1.61, 1.98; RSE 5.14%)
    e_wt_ktr  <- fixed(0.629) ; label("Allometric-like exponent of WT on k_tr (unitless; ref WT 5.9 kg)")    # Table II: gamma_wt = 0.629 (Fixed)

    # ----------------------------------------------------------------------
    # OS time-to-event sub-model (Table III of Harrold 2020).
    # Fit separately per pivotal study (filgrastim vs pegfilgrastim) because
    # an unexplained study effect remained after covariate screening
    # (Discussion p. 7). The binary STUDY_HARROLD_PEG indicator (0 = filg,
    # 1 = peg) selects between the two parameter sets inside model().
    # lambda_ANC and lambda_BC are not log-transformed because they can be
    # negative (lambda_ANC is a hazard slope; lambda_BC is the Box-Cox
    # power parameter). k_e0 is strictly positive and is reported on the
    # log scale (lke0_*).
    # ----------------------------------------------------------------------
    e_anc_haz_filg <- -2.15     ; label("Filgrastim study: hazard slope lambda_ANC on Box-Cox(ANCe) (Eq. 10)")  # Table III filgrastim: lambda_ANC = -2.15 (95% CI -3.32, -0.966; RSE 28%)
    e_anc_haz_peg  <- -0.229    ; label("Pegfilgrastim study: hazard slope lambda_ANC on Box-Cox(ANCe) (Eq. 10)") # Table III pegfilgrastim: lambda_ANC = -0.229 (95% CI -1.12, -0.00638; RSE 58%)
    lambda_bc_filg <- -0.347    ; label("Filgrastim study: Box-Cox power parameter lambda_BC (unitless)")        # Table III filgrastim: lambda_BC = -0.347 (95% CI -0.616, -0.0772; RSE 40%)
    lambda_bc_peg  <-  0.300    ; label("Pegfilgrastim study: Box-Cox power parameter lambda_BC (unitless)")     # Table III pegfilgrastim: lambda_BC = 0.300 (95% CI -0.107, 0.915; RSE 39%)
    lke0_filg      <- log(0.668); label("Filgrastim study: effect-compartment equilibration rate k_e0 (1/day)")  # Table III filgrastim: k_e0 = 0.668 (95% CI 0.592, 0.745; RSE 6%)
    lke0_peg       <- log(0.156); label("Pegfilgrastim study: effect-compartment equilibration rate k_e0 (1/day)") # Table III pegfilgrastim: k_e0 = 0.156 (95% CI 0.0221, 0.374; RSE 46%)

    # ----------------------------------------------------------------------
    # IIV on the ANC response parameters (Eq. 12; log-normal).
    # Table II reports %CV; the internal variance is omega^2 = log(1 + CV^2).
    # Paper notes that correlations between random effects were generally
    # low (r^2 < 0.16) "except for eta_kPD,e and eta_kPD,kill, which were
    # already incorporated as covariance components in the model" (ANC
    # Response, p. 6). The numerical covariance between eta_lkel and
    # eta_lkdecay is not published; this file declares the variances on
    # their own diagonal and notes the omission in the vignette Errata.
    # ----------------------------------------------------------------------
    etalrbase   ~ 0.0520   ; label("IIV variance on lrbase (log-normal; CV 23.1%)")    # Table II: omega_ANCIC %CV = 23.1; log(1 + 0.231^2) = 0.05201
    etalkel     ~ 0.01435  ; label("IIV variance on lkel (log-normal; CV 12.0%)")      # Table II: omega_kPD,e %CV = 12; log(1 + 0.12^2) = 0.01435
    etalkdecay  ~ 0.02380  ; label("IIV variance on lkdecay (log-normal; CV 15.5%)")   # Table II: omega_kPD,kill %CV = 15.5; log(1 + 0.155^2) = 0.02380
    etalkout    ~ 0.12715  ; label("IIV variance on lkout (log-normal; CV 36.8%)")     # Table II: omega_kc %CV = 36.8; log(1 + 0.368^2) = 0.12715

    # ----------------------------------------------------------------------
    # Residual error (Eq. 13). The paper reports an exponential error model
    # additive on the log scale; in nlmixr2 / rxode2 this is equivalent to
    # a proportional residual on the linear scale with SD equal to the
    # paper's reported CV fraction.
    # ----------------------------------------------------------------------
    propSd <- 0.620        ; label("Proportional residual error on ANC (fraction; 62% CV)")  # Table II: sigma_PE %CV = 62 (95% CI 59.1, 64.8; RSE 4.64%)
  })

  model({
    # Individual log-normal parameters
    rbase   <- exp(lrbase  + etalrbase)
    ktr     <- exp(lktr)                * (WT / 5.9) ^ e_wt_ktr   # WT effect on k_tr (Eq. 14), ref body weight 5.9 kg (Table I combined median)
    kout    <- exp(lkout   + etalkout)
    kel     <- exp(lkel    + etalkel)
    kdecay  <- exp(lkdecay + etalkdecay)
    gamma   <- exp(lgamma)

    # Zero-order stem-cell production rate kp such that the system is at
    # homeostasis when RAD = 0: at steady state Eq. 8 gives ANC(0) = kp/k_c
    # and N_SM(0) = N_MT(0) = N_PM1(0) = N_PM2(0) = kp/k_tr, so
    # kp = rbase * kout = ANC_IC * k_c.
    kp <- rbase * kout

    # Study-specific OS hazard parameters selected by the STUDY_HARROLD_PEG
    # indicator (0 = filgrastim, 1 = pegfilgrastim; Table III).
    lambda_anc <- e_anc_haz_filg  * (1 - STUDY_HARROLD_PEG) + e_anc_haz_peg  * STUDY_HARROLD_PEG
    lambda_bc  <- lambda_bc_filg  * (1 - STUDY_HARROLD_PEG) + lambda_bc_peg  * STUDY_HARROLD_PEG
    ke0        <- exp(lke0_filg) * (1 - STUDY_HARROLD_PEG) + exp(lke0_peg) * STUDY_HARROLD_PEG

    # ------------------------------------------------------------------
    # K-PD radiation compartment (Eq. 1). RAD enters as a 1 KPD-unit
    # bolus at t = 0 (cmt = "depot_kpd"); the depot decays exponentially
    # at k_PD,e and drives mitotic cell killing.
    # ------------------------------------------------------------------
    d/dt(depot_kpd) <- -kel * depot_kpd

    # K-PD-driven mitotic cell-loss rate (Eq. 2). The published Eq. 2
    # prints "k_kill = -k_PD,kill * RAD^gamma" with a leading minus
    # sign that is inconsistent with the positive Table II estimate
    # k_PD,kill = 2.14 (1/day * 1/KPD^gamma) and with the biological
    # direction (radiation increases mitotic cell loss). This file
    # follows the biologically- and Table-II-consistent interpretation
    # k_kill = +k_PD,kill * RAD^gamma; see vignette Errata for the
    # sign-convention discussion.
    rad_kill <- kdecay * depot_kpd ^ gamma

    # ------------------------------------------------------------------
    # Granulopoiesis chain (Eqs. 3-7).
    #   precursor1 = N_SM  : bone-marrow stem cells (Eq. 3; zero-order
    #                        input kp, first-order outflow k_tr)
    #   precursor2 = N_MT  : mitotic phase  (Eq. 4; radiation kill term)
    #   precursor3 = N_PM1 : maturation transit 1 (Eq. 5)
    #   precursor4 = N_PM2 : maturation transit 2 (Eq. 6)
    #   circ       = ANC   : circulating absolute neutrophils (Eq. 7)
    # ------------------------------------------------------------------
    d/dt(precursor1) <- kp - ktr * precursor1
    d/dt(precursor2) <- ktr * precursor1 - (ktr + rad_kill) * precursor2
    d/dt(precursor3) <- ktr * precursor2 - ktr * precursor3
    d/dt(precursor4) <- ktr * precursor3 - ktr * precursor4
    d/dt(circ)       <- ktr * precursor4 - kout * circ

    # Initial conditions at homeostasis (Eq. 8): N_SM(0) = N_MT(0) =
    # N_PM1(0) = N_PM2(0) = kp / k_tr; ANC(0) = kp / k_c = rbase.
    precursor1(0) <- kp / ktr
    precursor2(0) <- kp / ktr
    precursor3(0) <- kp / ktr
    precursor4(0) <- kp / ktr
    circ(0)       <- rbase

    # ------------------------------------------------------------------
    # Overall survival sub-model (Eqs. 9-11).
    #   - effect compartment ANCe equilibrates with circ at rate k_e0
    #     (Eq. 11); ANCe(0) = rbase so that the system starts at
    #     steady state with respect to OS as well.
    #   - hazard lambda(t) (Eq. 10, log-domain): log(lambda) is a
    #     linear function of the Box-Cox-transformed ANCe with slope
    #     lambda_ANC and power parameter lambda_BC. Both lambda_ANC
    #     and lambda_BC are study-specific (Table III).
    #   - the cumulative-hazard ODE integrates lambda(t) (Eq. 9) and
    #     the survival probability is sur = exp(-cumhaz).
    # ------------------------------------------------------------------
    d/dt(effect) <- ke0 * (circ - effect)
    effect(0)    <- rbase

    # Box-Cox transformation of ANCe (Eq. 10). The two reported
    # lambda_BC values (-0.347, +0.300) are well away from 0 so the
    # standard form (x^p - 1) / p is finite; no L'Hopital fallback is
    # needed at the reported parameter values.
    ancE_bc <- (effect ^ lambda_bc - 1) / lambda_bc

    # Hazard (log domain), cumulative hazard, and survival probability.
    log_hazard <- lambda_anc * ancE_bc
    hazard     <- exp(log_hazard)
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur <- exp(-cumhaz)

    # ------------------------------------------------------------------
    # Observation variable and residual error (Eq. 13).
    # ANC is the observed circulating neutrophil count (10^9 cells/L);
    # the paper's "exponential error additive on the log scale" maps
    # one-to-one onto rxode2's `prop(propSd)` linear-space proportional
    # error.
    # ------------------------------------------------------------------
    ANC <- circ
    ANC ~ prop(propSd)
  })
}
