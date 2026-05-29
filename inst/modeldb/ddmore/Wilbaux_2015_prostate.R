Wilbaux_2015_prostate <- function() {
  description <- "Joint semi-mechanistic kinetic-pharmacodynamic (K-PD) model of circulating tumour cell (CTC) count and prostate-specific antigen (PSA) longitudinal kinetics during chemotherapy and/or hormonotherapy in adults with metastatic castration-resistant prostate cancer (mCRPC) (Wilbaux 2015 / DDMODEL00000261). The structural model couples (i) two parallel first-order K-PD compartments for chemotherapy and hormonotherapy (no PK data -- virtual unit doses per cycle), (ii) a latent, dimensionless tumour-burden variable LV(t) governed by an indirect-response ODE with saturable Emax inhibition by both treatment compartments and a steady-state-anchored production rate, (iii) a CTC count that is the difference between two integrals of K0 * LV separated by the cell-lifespan delay LS (cell lifespan model implemented via a parallel delayed copy of the LV dynamics), and (iv) a PSA concentration following a non-steady-state indirect-response ODE driven by the same delayed LV. The published likelihood combines a negative-binomial count distribution for CTC observations (mean = alpha * CTC_total, overdispersion OVDP) with an exponential residual error on log-transformed PSA observations (W1 = 0.30); the full BLOCK(9) correlated $OMEGA across the nine PD parameters is preserved verbatim. nlmixr2's parser cannot natively express the .mod's F_FLAG=2 (-2*ln-likelihood) negative-binomial branch, so this implementation provides typical-value mechanistic outputs (NCTC = CTC * alpha, PSA, log_PSA) with placeholder additive residual errors on each -- see vignette Assumptions and deviations."
  reference <- paste(
    "Wilbaux M, Tod M, De Bono J, Lorente D, Mateo J, Freyer G, You B,",
    "Henin E. (2015).",
    "A joint model for the kinetics of CTC count and PSA concentration",
    "during treatment in metastatic castration-resistant prostate cancer.",
    "CPT Pharmacometrics Syst Pharmacol 4(5):277-285.",
    "doi:10.1002/psp4.34.",
    "DDMORE Foundation Model Repository: DDMODEL00000261.",
    sep = " "
  )
  vignette <- "Wilbaux_2015_prostate"
  paper_specific_compartments <- c("chemo", "hormo", "latent_tumor", "ctc", "chemo_d", "hormo_d", "latent_tumor_d")

  units <- list(
    time          = "day",
    dosing        = "AU (arbitrary unit per treatment cycle; the K-PD chemo/hormo compartments accept AMT = 1 per administered cycle)",
    concentration = "ng/mL (PSA endpoint); CTC count reported as cells per 7.5 mL aliquot (alpha = 0.0015 scales total-body CTC to aliquot count)"
  )

  ddmore_id    <- "DDMODEL00000261"
  replicate_of <- NULL

  covariateData <- list(
    # No subject-level covariates were retained in the final published model
    # (Wilbaux 2015 Results: "No covariates were identified in the final
    # model specification"). Treatment is encoded via dosing events in the
    # K-PD compartments (CMT 1 / 2 for current chemo / hormonotherapy and
    # CMT 5 / 6 for their lagged copies), not as covariate columns.
  )

  population <- list(
    n_subjects     = 223L,
    n_studies      = 1L,
    age_range      = "adults with metastatic castration-resistant prostate cancer (mCRPC); detailed age distribution not extracted in this run (PMC text quoted only the n = 223 cohort summary)",
    weight_range   = "not retained as a model covariate; not extracted",
    sex_female_pct = 0,
    disease_state  = "Metastatic castration-resistant prostate cancer (mCRPC) with frequent non-measurable bone metastases. Subjects received chemotherapy (most commonly docetaxel-based regimens) and/or hormonotherapy (most commonly abiraterone or enzalutamide) over a median follow-up of approximately six months (paper Methods; the .mod's K-PD compartments treat each cycle as an arbitrary AMT = 1 event with no chemical-distinction between agents).",
    dose_range     = "Chemotherapy and/or hormonotherapy administered per cycle. The K-PD parameterization is dose-/concentration-agnostic: each treatment cycle enters the model as AMT = 1 in the corresponding K-PD compartment (CMT 1 = chemo, CMT 2 = hormo) at the cycle start time, plus an identical dose at AMT = 1 to the lagged copies (CMT 5, CMT 6) that share the same lag time LS = 58 days. The bundle's Simulated_KPD_CTC.count_PSA.csv encodes 1-cycle / 2-cycle / 3-cycle administration patterns at days 0, 21-35, 42-77, 85, 108-119 across two simulated subjects.",
    biomarkers     = "Two longitudinal biomarkers: (i) CTC count per 7.5 mL aliquot (CellSearch system; non-negative integer; bundle CMT 4) and (ii) PSA concentration in plasma (ng/mL; bundle CMT 8 with DV recorded as log(PSA) in the simulated dataset). Latent tumour burden LV(t) is dimensionless and not directly observed.",
    notes          = "n_subjects = 223 from Wilbaux 2015 paper text (Methods, via PMC4452933 abstract / Methods extract). Detailed baseline demographics (age distribution, ECOG, Gleason, prior therapy) are in the linked publication (CPT Pharmacometrics Syst Pharmacol 2015;4(5):277-285, doi:10.1002/psp4.34); the paper PDF was not on disk in /home/bill/github/mab_human_consensus/literature at extraction time, so this population block reports only fields that are directly recoverable from the PMC full-text Methods extract performed during extraction. The DDMORE bundle's Simulated_KPD_CTC.count_PSA.csv contains only two simulated subjects with very sparse observation schedules (intentionally minimal regression dataset, not representative of the published cohort)."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Output_real_SAEM_KPD_CTC.count_PSA.lst,
    # FINAL PARAMETER ESTIMATE block (post-SAEM convergence).
    #
    # The bundle .mod's $THETA / $OMEGA blocks (initial values) match the
    # .lst's final estimates verbatim for this entry -- i.e. the .mod was
    # supplied with the published values. Cross-check against
    # Wilbaux 2015 Table 1 was performed via PMC4452933 full-text Methods
    # / Tables; all 14 thetas and the BLOCK(9) / BLOCK(2) OMEGAs match
    # within rounding (paper rounds A50c = 0.0003 vs .lst 2.61e-4;
    # KoutLV = 0.00513 vs .lst 5.13e-3; A50h = 0.004 vs .lst 3.97e-3;
    # other values agree to three significant figures). The paper's IIV
    # column is reported as `sqrt(omega^2) * 100%` (i.e. SD on the
    # estimation scale, in percent) rather than the back-transformed CV%.
    # See vignette Source trace for the full mapping.
    #
    # The .mod uses NONMEM SAEM with MU-referencing. THETA(7) (-> SFLV)
    # and THETA(11), THETA(12) (-> K0, LS) are NOT log-transformed via
    # MU: MU_7 = THETA(7), MU_11 = THETA(11), MU_12 = THETA(12). For
    # SFLV the EXP() wrapper still applies (TH = exp(MU_7 + ETA)), so
    # SFLV is on log scale internally and ETA(7) is on log scale. For
    # K0 and LS the EXP() wrapper is absent (K0 = MU_11 + ETA(11);
    # ALAG5 = MU_12 + ETA(12)) -- i.e. K0 and LS are estimated on the
    # **linear** scale and their ETAs are **additive** with linear-scale
    # variances (1430 (CTC*day/AU)^2 for K0; 61.9 day^2 for LS). The
    # nlmixr2 implementation preserves this: k0 / lag5 are plain
    # (un-logged) ini parameters and the etas are added linearly in
    # model() (`k0_ind = k0 + etak0`, `lag5_ind = lag5 + etalag5`).
    # ----------------------------------------------------------------------

    # --- Latent tumour-variable initial / steady-state (LV0; AU) ---
    # The paper fixes LV0 = 1 ("LV corresponded to the fractional change
    # in the latent variable from baseline"). The .mod fixes THETA(1) at
    # 1 and OMEGA(1,1) at 1e-7 FIX (effectively no IIV).
    lrbase <- fix(log(1)); label("Log of latent tumour-variable initial / baseline value LV0 (AU; fixed at 1 per Wilbaux 2015 Methods; .mod TH(1) = 1 FIX)")  # .lst FINAL TH 1 = 1.00E+00 (FIX)

    # --- K-PD elimination rates (1/day) ---
    lk1 <- log(0.248);  label("Log of chemotherapy K-PD elimination rate constant Kc (1/day)") # .lst FINAL TH 2 = 2.48E-01; paper Table 1 Kc = 0.248
    lk2 <- log(0.449);  label("Log of hormonotherapy K-PD elimination rate constant Kh (1/day)") # .lst FINAL TH 3 = 4.49E-01; paper Table 1 Kh = 0.449

    # --- Half-maximal-effect amounts in the K-PD compartments (AU) ---
    # Saturable Emax inhibition of LV production; Q501 = A50c, Q502 = A50h
    # (paper notation). With the AMT = 1 K-PD dosing convention, A50 < 1
    # means the inhibition saturates within a single cycle.
    lq501 <- log(2.61e-4); label("Log of half-maximal-effect amount of chemo K-PD compartment A50c (AU)")  # .lst FINAL TH 4 = 2.61E-04; paper Table 1 A50c = 0.0003 (rounded)
    lq502 <- log(3.97e-3); label("Log of half-maximal-effect amount of hormo K-PD compartment A50h (AU)")  # .lst FINAL TH 5 = 3.97E-03; paper Table 1 A50h = 0.004 (rounded)

    # --- Latent-variable elimination (1/day) ---
    lkoutts <- log(5.13e-3); label("Log of latent tumour-variable elimination rate KoutLV (1/day)") # .lst FINAL TH 6 = 5.13E-03; paper Table 1 KoutLV = 0.00513

    # --- Latent-variable steady-state scaling factor SFLV (AU*day, log scale) ---
    # Used inside model() to derive KinLV via
    #   KinLV = LV0 * KoutLV * (1 + SFLV) / SFLV
    # (NONMEM .mod $PK lines: KINTS = (TS0 * KOUTTS) / (TH / (1 + TH))).
    # The .mod parameterization is on the **log** scale because of the
    # `MU_7 = THETA(7); TH = EXP(MU_7 + ETA(7))` MU-reference pattern;
    # paper Table 1 reports SFLV = 6.33 (which is the THETA value, i.e.
    # log-scale; the back-transformed asymptote SFLV_lin = exp(6.33) is
    # large enough that KinLV ~= LV0 * KoutLV at typical values -- i.e.
    # the latent variable is approximately at steady state at LV0).
    lsflv <- 6.33; label("Log of latent-variable steady-state scaling factor SFLV (AU*day, log scale per .mod MU-reference pattern)") # .lst FINAL TH 7 = 6.33E+00; paper Table 1 SFLV = 6.33

    # --- PSA indirect-response turnover (ng/mL/day/AU and 1/day) ---
    lkinpsa  <- log(1.40);    label("Log of PSA production rate constant KinPSA (ng/mL/day per LV unit)")  # .lst FINAL TH 8 = 1.40E+00; paper Table 1 KinPSA = 1.40
    lkoutpsa <- log(8.13e-3); label("Log of PSA elimination rate constant KoutPSA (1/day)")              # .lst FINAL TH 9 = 8.13E-03; paper Table 1 KoutPSA = 0.00813

    # --- PSA initial concentration (ng/mL) ---
    lpsa0 <- log(153); label("Log of typical baseline PSA concentration PSA0 (ng/mL)")  # .lst FINAL TH10 = 1.53E+02; paper Table 1 PSA0 = 153

    # --- CTC production rate K0 and cell lifespan LS (LINEAR-scale ETAs) ---
    # The .mod estimates K0 and LS on the linear scale (no LOG() wrapper
    # in the MU-reference: MU_11 = THETA(11), K0 = MU_11 + ETA(11);
    # MU_12 = THETA(12), ALAG5 = MU_12 + ETA(12)). The ETAs are therefore
    # additive on the linear scale. Inside model() these are reconstructed
    # as `k0_ind = k0 + etak0` and `lag5_ind = lag5 + etalag5` (plain
    # addition, no exp()). The published CV% column uses CV = SD/typical,
    # e.g. K0: sqrt(1430)/308 ~ 0.123 = 12% CV.
    k0   <- 308;  label("Typical CTC zero-order production rate constant K0 (CTC*day^-1*AU^-1) -- linear-scale ini parameter; etak0 is additive in model()")  # .lst FINAL TH11 = 3.08E+02; paper Table 1 K0 = 308
    lag5 <- 57.7; label("Typical CTC cell lifespan LS (day); shared as the lag time on the chemo / hormo / latent-variable lagged-copy compartments -- linear-scale ini parameter; etalag5 is additive in model()") # .lst FINAL TH12 = 5.77E+01; paper Table 1 LS = 58

    # --- Negative-binomial overdispersion (AU, dimensionless) ---
    # OVDP > 0 increases variance of the count distribution above its mean
    # (Var = mu + OVDP * mu^2). OVDP -> 0 reduces to Poisson.
    lovdp <- log(4.89); label("Log of negative-binomial overdispersion parameter OVDP (unitless)")  # .lst FINAL TH13 = 4.89E+00; paper Table 1 OVDP = 4.9

    # --- Log-PSA exponential residual error SD (W1; unitless on log scale) ---
    # NONMEM .mod $ERROR: Y = LOG(PSA) + W1 * ERR(1) with $SIGMA 1 FIX,
    # i.e. additive residual on log(PSA) with SD = W1. Paper reports
    # PSA res error = 0.3.
    lw1 <- log(0.300); label("Log of exponential residual-error SD W1 on log-transformed PSA (unitless on log scale)") # .lst FINAL TH14 = 3.00E-01; paper Table 1 PSA res error = 0.3

    # ------------------------------------------------------------------
    # Inter-individual variability -- full BLOCK(9) correlated $OMEGA on
    # ETA(2) .. ETA(10) (mapped to {Kc, Kh, A50c, A50h, KoutLV, SFLV,
    # KinPSA, KoutPSA, PSA0}). NONMEM stores the lower-triangle row by
    # row; the c(...) vector below mirrors that order so the resulting
    # nlmixr2 covariance matrix matches the .lst's OMEGA(i,j) printout
    # element-for-element. Reading in (variance, off-diagonal) form:
    #   Row 2: var(K1)
    #   Row 3: cov(K2, K1), var(K2)
    #   Row 4: cov(Q501, K1), cov(Q501, K2), var(Q501)
    #   ... etc.
    # ------------------------------------------------------------------
    etalk1 + etalk2 + etalq501 + etalq502 + etalkoutts + etalsflv + etalkinpsa + etalkoutpsa + etalpsa0 ~ c(
      0.723,                                                                   # OMEGA(2,2)  -- var Kc
      -0.698, 1.83,                                                            # OMEGA(3,2..3)
      -0.638, -0.634, 4.76,                                                    # OMEGA(4,2..4)
      -0.799, -0.228, 2.93, 2.84,                                              # OMEGA(5,2..5)
      0.351, 2.25, -2.65, -4.24, 20.6,                                         # OMEGA(6,2..6)
      0.237, 0.04, -0.751, -0.954, 3.43, 0.786,                                # OMEGA(7,2..7)
      0.131, 0.974, -1.81, -1.45, -0.604, -0.0975, 2.60,                       # OMEGA(8,2..8)
      0.312, 0.0451, -2.29, -1.16, 1.84, 0.578, 0.0967, 1.53,                  # OMEGA(9,2..9)
      0.115, 0.693, -1.31, -1.25, -0.302, 0.153, 2.30, 0.0348, 2.40            # OMEGA(10,2..10)
    )

    # BLOCK(2) on ETA(11), ETA(12) -- additive (linear-scale) IIV on K0
    # and LS. Variances are in (CTC*day/AU)^2 and day^2.
    etak0 + etalag5 ~ c(1430, 294, 61.9)  # OMEGA(11,11), OMEGA(12,11), OMEGA(12,12)

    # Diagonal $OMEGA on ETA(13) -- log-scale IIV on OVDP.
    etalovdp ~ 2.25  # OMEGA(13,13)

    # The .mod also declares ETA(1) on LV0 (FIX 1e-7) and ETA(14) on W1
    # (FIX 1e-7). Both are treated as fixed-typical-only here: LV0 is
    # supplied as `lrbase <- fix(log(1))` above; W1 has no IIV and only
    # carries the population value `lw1 <- log(0.300)`.

    # ------------------------------------------------------------------
    # Placeholder additive residual-error SDs for the typical-value
    # mechanistic outputs (NCTC, log_PSA). The published likelihood is
    #   CTC ~ NegBin(mean = alpha * CTC_total, overdisp = OVDP)
    # for the count branch (.mod F_FLAG = 2, Y = -2 * log L), and
    #   log(PSA) ~ Normal(log(IPRED), W1)
    # for the continuous branch. nlmixr2's parser cannot natively
    # express the F_FLAG = 2 negative-binomial branch (custom -2*ln-L
    # forms outside dnbinom* / dpois are not parser-supported as
    # observation likelihoods), so the implementation here exposes
    # NCTC and log_PSA as deterministic typical-value outputs with
    # placeholder additive residual SDs. The log_PSA placeholder is
    # initialised to W1 = 0.30 (the published exponential residual
    # SD on log(PSA)); the NCTC placeholder is a small fixed value
    # the user can override at fit time. See vignette Assumptions
    # and deviations for the full deviation rationale.
    # ------------------------------------------------------------------
    addSd_NCTC    <- 1   ; label("Placeholder additive residual-error SD on NCTC typical-value output (cells per 7.5 mL aliquot); not the source likelihood -- see vignette Assumptions and deviations")
    addSd_log_PSA <- 0.300; label("Additive residual-error SD on log_PSA typical-value output (matches the published W1 exponential residual SD on log(PSA)); .lst FINAL TH14 = 3.00E-01")
  })

  model({
    # ------------------------------------------------------------------
    # Individual structural parameters. Log-scale parameters use the
    # standard exp(typical + eta) reconstruction; the two LINEAR-scale
    # parameters K0 and LS use additive (typical + eta) reconstruction
    # to match the .mod's MU-reference (no LOG() on THETA(11), THETA(12)).
    # ------------------------------------------------------------------
    rbase       <- exp(lrbase)
    k1        <- exp(lk1     + etalk1)
    k2        <- exp(lk2     + etalk2)
    q501      <- exp(lq501   + etalq501)
    q502      <- exp(lq502   + etalq502)
    koutts    <- exp(lkoutts + etalkoutts)
    sflv      <- exp(lsflv   + etalsflv)        # SFLV on linear scale (= TH in .mod)
    kinpsa    <- exp(lkinpsa + etalkinpsa)
    koutpsa   <- exp(lkoutpsa + etalkoutpsa)
    psa0      <- exp(lpsa0   + etalpsa0)
    k0_ind    <- k0   + etak0                   # additive ETA on linear scale
    lag5_ind  <- lag5 + etalag5                 # additive ETA on linear scale
    ovdp      <- exp(lovdp + etalovdp)
    w1        <- exp(lw1)                       # no IIV (OMEGA fix at 1e-7)

    # Latent-variable production rate. .mod $PK: KINTS = TS0 * KOUTTS /
    # (TH / (1 + TH)) = TS0 * KOUTTS * (1 + SFLV) / SFLV.
    kints <- rbase * koutts * (1 + sflv) / sflv

    # Saturable Emax inhibition by both K-PD compartments. Defined here
    # so the d/dt(latent_tumor) line below can reference them in CMT
    # order; rxode2 parses derivations top-to-bottom.
    inh_chemo   <- chemo   / (q501 + chemo)
    inh_hormo   <- hormo   / (q502 + hormo)
    inh_chemo_d <- chemo_d / (q501 + chemo_d)
    inh_hormo_d <- hormo_d / (q502 + hormo_d)

    # The .mod enforces that the "delayed" latent-variable copy returns
    # the initial LV0 (rather than the integrated A(7)) for t <= ALAG5,
    # via
    #   A7 = TS0
    #   IF (T .GT. ALAG5) A7 = A(7)
    # This ensures the cell-lifespan integral for CTC starts cleanly
    # from baseline. Reproduce by guarding the CTC ODE input.
    latent_tumor_d_eff <- ifelse(t > lag5_ind, latent_tumor_d, rbase)

    # ------------------------------------------------------------------
    # ODE system. Compartments are declared in CMT order (1..8) because
    # the simulated dataset references compartments by NUMBER:
    #   CMT 1 = chemo            CMT 2 = hormo
    #   CMT 3 = latent_tumor     CMT 4 = ctc
    #   CMT 5 = chemo_d          CMT 6 = hormo_d
    #   CMT 7 = latent_tumor_d   CMT 8 = psa
    # ------------------------------------------------------------------

    # CMT 1: chemotherapy K-PD compartment. AMT = 1 dose per cycle.
    d/dt(chemo) <- -k1 * chemo

    # CMT 2: hormonotherapy K-PD compartment. AMT = 1 dose per cycle.
    d/dt(hormo) <- -k2 * hormo

    # CMT 3: latent tumour variable LV(t). Indirect-response ODE with
    # saturable Emax inhibition by both K-PD compartments.
    d/dt(latent_tumor) <- kints * (1 - inh_chemo) * (1 - inh_hormo) - koutts * latent_tumor
    latent_tumor(0) <- rbase

    # CMT 4: CTC count. The .mod implements the cell-lifespan model via
    #   DADT(4) = K0 * A(3) - K0 * A7
    # where A7 is the latent variable delayed by LS (via the parallel
    # CMT 5 / 6 / 7 dynamics with lag time ALAG5). At t = 0 the
    # bioavailability F(ctc) = K0 * LS sets the steady-state initial
    # CTC count given LV(0) = LV0:
    #   CTC(0) = K0 * LV0 * LS    (= K0 * 1 * LS = K0 * LS).
    # The simulated dataset issues an AMT = 1 dose to CMT 4 at t = 0;
    # the bioavailability multiplier scales it to the steady-state
    # initial population.
    d/dt(ctc) <- k0_ind * latent_tumor - k0_ind * latent_tumor_d_eff
    f(ctc)    <- k0_ind * lag5_ind

    # CMT 5: lagged chemotherapy K-PD compartment. Same K-PD dynamics
    # as chemo, but doses arrive lag-shifted by lag5_ind. The simulated
    # dataset issues parallel AMT = 1 doses to CMT 5 alongside the CMT 1
    # doses (same time, same amount); the alag(chemo_d) <- lag5_ind
    # declaration delays each dose's arrival in the chemo_d compartment.
    d/dt(chemo_d) <- -k1 * chemo_d
    alag(chemo_d) <- lag5_ind

    # CMT 6: lagged hormonotherapy K-PD compartment. Same as CMT 5 for
    # hormo. Doses to CMT 6 arrive in the compartment after a delay of
    # lag5_ind (= LS).
    d/dt(hormo_d) <- -k2 * hormo_d
    alag(hormo_d) <- lag5_ind

    # CMT 7: lagged latent variable. Same dynamics as latent_tumor but
    # using the lagged K-PD compartments as inputs. Since CMT 7 has no
    # dose in the .mod / dataset, ALAG7 in the .mod is harmless. The
    # initial condition equals LV0; during 0 < t < lag5_ind the
    # inhibition terms are 0 (chemo_d = hormo_d = 0 because their doses
    # haven't arrived yet) and dlatent_tumor_d/dt = kints - koutts * LV
    # equals 0 at LV = LV0 (since kints / koutts = (1 + sflv) / sflv * LV0
    # ~= LV0 at typical sflv), so the lagged copy stays approximately at
    # LV0 until t > lag5_ind.
    d/dt(latent_tumor_d) <- kints * (1 - inh_chemo_d) * (1 - inh_hormo_d) - koutts * latent_tumor_d
    latent_tumor_d(0) <- rbase

    # CMT 8: PSA concentration. Indirect-response ODE driven by the
    # CURRENT (un-lagged) latent variable per .mod $DES line:
    #   DADT(8) = KINPSA * A(3) - KOUTPSA * A(8)
    # Initial condition at the published baseline.
    d/dt(psa) <- kinpsa * latent_tumor - koutpsa * psa
    psa(0) <- psa0

    # ------------------------------------------------------------------
    # Observation outputs.
    # ------------------------------------------------------------------
    # NCTC: CTC count per 7.5 mL aliquot (alpha = 0.0015 scales the
    # total-body CTC compartment to the per-aliquot count). The .mod
    # uses Y = -2 * log L from a negative-binomial likelihood with
    # mean = NCTC and overdispersion OVDP; nlmixr2's parser does not
    # natively express that custom F_FLAG = 2 form, so this output is
    # exposed as a typical-value mechanistic prediction with a
    # placeholder additive residual error.
    NCTC <- ctc * 0.0015
    NCTC ~ add(addSd_NCTC)

    # log_PSA: natural log of PSA. The .mod's published likelihood is
    # additive on log(PSA) with SD = W1 (= 0.30); the simulated dataset
    # records DV in log-PSA units. The addSd_log_PSA ini parameter is
    # initialised to the published W1 = 0.30 so this output's residual
    # error matches the publication for typical-value validation.
    log_PSA <- log(psa)
    log_PSA ~ add(addSd_log_PSA)
  })
}
