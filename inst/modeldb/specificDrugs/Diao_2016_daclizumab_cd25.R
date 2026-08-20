Diao_2016_daclizumab_cd25 <- function() {
  description <- "Kinetic-binding PK/PD model of CD25 receptor occupancy on peripheral CD4+ T cells following subcutaneous daclizumab high-yield process (HYP) in adults with relapsing-remitting multiple sclerosis (Diao 2016). The receptor-occupancy state occ_cd25 obeys dOcc/dt = kon*Cc*(1 - Occ) - koff*Occ; the PD output is the percentage of CD4+ T cells staining positive for unoccupied CD25, computed as E0 * (1 - occ_cd25). This is a REPARAMETERISATION of the paper's published sigmoidal Emax (Equation 1) with its two-parameter-set (saturation-phase / desaturation-phase) Hill function -- kon and koff are calibrated to reproduce the paper's phenomenology (rapid saturation within 7 h after first 150 mg SC dose per Figure 1A; return of unoccupied CD25 to baseline in ~24 weeks after last steady-state dose per Figure 1B) rather than transcribed from a paper table. See the companion vignette Assumptions and deviations for the calibration rationale. The PK backbone is inherited verbatim from Othman 2014."
  reference <- "Diao L, Hang Y, Othman AA, Nestorov I, Tran JQ, Mehta D, Amaravadi L. Population PK/PD analyses of CD25 occupancy, CD56 bright NK cell expansion and regulatory T cell reduction by daclizumab HYP in subjects with multiple sclerosis. Br J Clin Pharmacol. 2016;82(5):1333-1342. doi:10.1111/bcp.13051 (PMID 27333593). PK backbone: Othman AA, Tran JQ, Tang MT, Dutta S. Population Pharmacokinetics of Daclizumab High-Yield Process in Healthy Volunteers. Clin Pharmacokinet. 2014;53(10):907-918. doi:10.1007/s40262-014-0159-9."
  vignette <- "Diao_2016_daclizumab_cd25"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL", response = "% of CD4+ T cells unoccupied CD25")

  paper_specific_compartments <- c("occ_cd25")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "daclizumab", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "daclizumab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "daclizumab", units = "mg", specimen = "plasma", verified = FALSE),
    occ_cd25    = list(analyte = "CD25 receptor occupancy", units = "mg", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used for allometric scaling of the inherited Othman 2014 PK parameters (CL, Q, Vc, Vp) with reference 70 kg; exponents 0.54 on CL/Q and 0.64 on Vc/Vp. CD25 PD parameters do not carry weight covariates in Diao 2016.",
      source_name        = "WT"
    ),
    DOSE_50MG = list(
      description        = "Record-level indicator for the 50 mg SC dose (1 = 50 mg SC, 0 = any other SC dose or any IV dose)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (100, 150, 200, or 300 mg SC dose, or any IV dose)",
      notes              = "Inherited from the Othman 2014 PK backbone for the dose-dependent bioavailability term. The Diao 2016 RRMS regimens are 150 or 300 mg SC every 4 weeks, so leave DOSE_50MG = 0 in clinical simulations. See Othman_2014_daclizumab.R for the full rationale.",
      source_name        = "(derived from AMT)"
    )
  )

  population <- list(
    n_subjects     = 1459L,
    n_records      = 7622L,
    n_studies      = 4L,
    study_names    = c("205MS201 / SELECT (Phase 2, RRMS)",
                       "205MS202 / SELECTION (Phase 2 extension with washout cohort)",
                       "205MS302 / OBSERVE (immunogenicity / PK / PD with intensive substudy)",
                       "205MS301 / DECIDE (Phase 3 vs IFN beta-1a)"),
    disease_state  = "Relapsing-remitting multiple sclerosis (RRMS)",
    dose_range     = "Daclizumab HYP 150 or 300 mg SC every 4 weeks",
    notes          = paste0(
      "Pooled PK/PD dataset of 1459 RRMS subjects with 7622 CD25 occupancy ",
      "records from four daclizumab HYP clinical studies (Diao 2016 Table 2). ",
      "Subject-level demographics are reported in the companion population PK ",
      "analysis (paper reference [13]); they are not enumerated in Diao 2016 ",
      "itself."
    ),
    pd_subgroups = list(
      `205MS201/202 (SELECT/SELECTION)` = list(subjects = 580L, records = 5123L),
      `205MS302 (OBSERVE)`              = list(subjects = 113L, records =  974L,
                                                intensive_PK_PD_subgroup = 25L),
      `205MS301 (DECIDE)`                = list(subjects = 766L, records = 1525L)
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # PK backbone (copied verbatim from inst/modeldb/specificDrugs/Othman_2014_daclizumab.R).
    # The Diao 2016 PD analysis used a sequential approach: PK was fixed
    # to a previously published population PK model (Diao 2016 Methods,
    # ref [13]). For library coherence the Othman 2014 healthy-volunteer
    # PK is used here as the canonical daclizumab HYP PK backbone; see
    # vignette Assumptions and deviations for the small numerical
    # differences between Othman 2014 and the in-paper PK summary.
    # ----------------------------------------------------------------------
    lka      <- log(0.009 * 24); label("Absorption rate constant (Ka, 1/day; 0.009 /h SC)")           # Othman 2014 Table 2
    lcl      <- log(0.010 * 24); label("Clearance for a 70 kg adult (CL, L/day; 10 mL/h)")            # Othman 2014 Table 2
    lvc      <- log(3.89);       label("Central volume of distribution for a 70 kg adult (Vc, L)")    # Othman 2014 Table 2
    lvp      <- log(2.52);       label("Peripheral volume of distribution for a 70 kg adult (Vp, L)") # Othman 2014 Table 2
    lq       <- log(0.044 * 24); label("Inter-compartmental clearance for a 70 kg adult (Q, L/day)")  # Othman 2014 Table 2
    lfdepot  <- log(0.84);       label("SC bioavailability for 100-300 mg doses (F, fraction)")       # Othman 2014 Table 2
    ltlag    <- log(2 / 24);     label("Absorption lag time for SC doses (Tlag, day; 2 h)")           # Othman 2014 Table 2

    e_wt_cl_q <- 0.54; label("Allometric exponent on CL and Q (unitless)")  # Othman 2014 Table 2
    e_wt_vc_vp <- 0.64; label("Allometric exponent on Vc and Vp (unitless)") # Othman 2014 Table 2

    e_dose_50mg_f <- -0.32143; label("Relative change in F for 50 mg SC vs 100-300 mg SC (fraction)") # Othman 2014 Table 2

    # PK IIV (Othman 2014 Table 2 SC-cohort values).
    # ka  CV 58% -> omega^2 = log(1 + 0.58^2) = 0.29003
    # cl  CV 27% -> omega^2 = log(1 + 0.27^2) = 0.07038
    # corr(ka, cl) SC = -0.72 -> cov = -0.72 * sqrt(0.29003 * 0.07038) = -0.10290
    etalka + etalcl ~ c(0.29003,
                        -0.10290, 0.07038)                                                       # Othman 2014 Table 2
    etalvc ~ 0.09175                                                                              # Othman 2014 Table 2 (Vc CV 31%)

    # PK residual error (Othman 2014 Table 2): proportional 22% + additive 0.33 ug/mL.
    propSd <- 0.22; label("Proportional residual error on daclizumab HYP serum concentration (fraction)") # Othman 2014 Table 2
    addSd  <- 0.33; label("Additive residual error on daclizumab HYP serum concentration (ug/mL)")        # Othman 2014 Table 2

    # ----------------------------------------------------------------------
    # CD25 kinetic-binding parameters (REPARAMETERISATION of Diao 2016 Eq. 1).
    #
    # The paper reports two sigmoidal Emax parameter sets that cannot be
    # reconciled with a single instantaneous Hill function (Diao 2016 Table 3):
    #
    #   saturation phase   IC50 = 0.0135 mg/L, gamma = 1     (both FIXED
    #                      from the OBSERVE intensive substudy; describe
    #                      the rapid initial binding)
    #   desaturation phase IC50 = 2.07 mg/L,  gamma = 4.44   (estimated
    #                      from the SELECTION washout cohort; describe the
    #                      slow return to baseline)
    #
    # A single equilibrium Hill can only match one of the two limits at a
    # time (the earlier merged-to-main version of this model used the
    # desaturation set only, under-predicting the rapid initial saturation
    # reported in Diao 2016 Figure 1A). To reproduce BOTH the ~7 h
    # saturation timescale and the ~24 week baseline-return timescale in a
    # single simulation-ready model, the receptor occupancy is instead
    # tracked as an ODE state occ_cd25 in [0, 1] with mass-action kinetic
    # binding:
    #
    #     d(occ_cd25)/dt = kon * Cc * (1 - occ_cd25) - koff * occ_cd25
    #
    # kon is set high so kon*Cc dominates koff at any measurable clinical
    # Cc (fast association -> Figure 1A saturation within hours of the
    # first dose). koff is set low so once Cc drops well below the effective
    # equilibrium IC50 (koff/kon) during washout, occ_cd25 decays with
    # t1/2 = ln(2)/koff = 21 days (slow dissociation -> Figure 1B
    # ~24-week baseline return).
    #
    # kon and koff are CALIBRATION values, not fitted paper parameters.
    # They are wrapped in fixed() to make the provenance unambiguous.
    # The observation-scale IIV (baseline additive; koff CV) is retained
    # from the paper's Table 3 estimates on the E0 and IC50 (desaturation)
    # rows so that the marginal spread of the simulated occupancy tracks
    # the paper's reported variability.
    # ----------------------------------------------------------------------
    lkon      <- fixed(log(24));    label("Kinetic binding association rate (kon, L/(mg*day); calibrated -- see Assumptions)")   # Calibrated re-parameterisation; targets rapid Figure 1A saturation
    lkoff     <- fixed(log(0.033)); label("Kinetic binding dissociation rate (koff, 1/day; t1/2 off = 21 days; calibrated)")     # Calibrated re-parameterisation; targets Figure 1B ~24-week baseline return
    cd25E0    <- 56;                label("Typical baseline unoccupied CD25 (% of CD4+ T cells)")                                # Diao 2016 Table 3 (Baseline E0 = 56%)

    # IIV.
    # cd25E0 IIV is Diao 2016 Table 3 as "(additive) 11" (i.e., additive on
    # the linear percentage scale), not CV%. Encoded as an additive eta on
    # cd25E0 in linear (percentage-point) space; variance = 11^2 = 121.
    # koff IIV = 47% CV -> omega^2 = log(1 + 0.47^2) = 0.19770. This is
    # inherited from the paper's desaturation-IC50 IIV row: with kon fixed,
    # IIV on koff is equivalent to IIV on the effective equilibrium IC50
    # (which equals koff/kon), preserving the between-subject variability
    # the paper reported.
    etacd25E0 ~ 121                                                                              # Diao 2016 Table 3 (Baseline E0 IIV additive SD = 11 percentage points; variance = 11^2)
    etalkoff  ~ 0.19770                                                                          # Diao 2016 Table 3 (Desaturation IC50 IIV 47% CV, remapped to kinetic binding koff)

    # Residual error on the CD25 occupancy observation.
    # Diao 2016 Table 3 reports "Residual error (additive) = 4.02" (units:
    # percentage points of CD4+ T cells).
    addSd_cd25 <- 4.02; label("Additive residual error on unoccupied CD25 (% of CD4+ T cells)")  # Diao 2016 Table 3
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters (Othman 2014 PK backbone).
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 2. Two-compartment SC PK ODE system.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- exp(ltlag)
    f(depot)    <- exp(lfdepot) * (1 + e_dose_50mg_f * DOSE_50MG)

    Cc <- central / vc

    # ------------------------------------------------------------------
    # 3. Individual CD25 kinetic-binding parameters.
    # ------------------------------------------------------------------
    kon_i    <- exp(lkon)
    koff_i   <- exp(lkoff + etalkoff)
    cd25E0_i <- cd25E0 + etacd25E0

    # ------------------------------------------------------------------
    # 4. Kinetic binding ODE for the fraction of CD25 bound by daclizumab HYP.
    #    State occ_cd25 in [0, 1]; initial condition 0 (no drug present).
    #    See Assumptions and deviations for the calibration rationale.
    # ------------------------------------------------------------------
    d/dt(occ_cd25) <- kon_i * Cc * (1 - occ_cd25) - koff_i * occ_cd25
    occ_cd25(0)    <- 0

    # ------------------------------------------------------------------
    # 5. Unoccupied CD25 percentage (observable).
    # ------------------------------------------------------------------
    cd25 <- cd25E0_i * (1 - occ_cd25)

    # ------------------------------------------------------------------
    # 6. Observation and error model (PK + PD outputs).
    # ------------------------------------------------------------------
    Cc   ~ add(addSd) + prop(propSd)
    cd25 ~ add(addSd_cd25)
  })
}
