Chen_2017_TB_MTP_GPDI_mouse <- function() {
  description <- "Preclinical (BALB/c mouse). Multistate Tuberculosis Pharmacometric (MTP) model linked to General Pharmacodynamic Interaction (GPDI) model describing CFU/lungs dynamics in M. tuberculosis Beijing VN 2002-1585 infected BALB/c mice receiving oral monotherapy or combination therapy with rifampicin, isoniazid, ethambutol, and pyrazinamide. Four parallel population PK models (1-cmt RIF/INH/PZA, 2-cmt EMB) drive concentration-dependent drug effects on three bacterial states (fast-multiplying F, slow-multiplying S, non-multiplying N). Rifampicin and isoniazid interact antagonistically on the killing of S and N (INT_S_RH = 4.49; INT_N_RH = 0.32); rifampicin and ethambutol interact synergistically on the killing of N (INT_N_RE = -0.15). The natural-growth transfer rates kSF, kFN, kSN, kNS are fixed from the in vitro MTP fit of Clewe 2016; the absorption rate constants for isoniazid, ethambutol, and pyrazinamide are fixed from the upstream mouse popPK of Chen 2016."
  reference <- paste(
    "Chen C, Wicha SG, de Knegt GJ, Ortega F, Alameda L, Sousa V,",
    "de Steenwinkel JEM, Simonsson USH. (2017).",
    "Assessing Pharmacodynamic Interactions in Mice Using the Multistate",
    "Tuberculosis Pharmacometric and General Pharmacodynamic Interaction Models.",
    "CPT Pharmacometrics Syst. Pharmacol. 6(11):787-797.",
    "doi:10.1002/psp4.12226.",
    "Upstream parameter sources retained inline:",
    "MTP natural-growth transfer rates fixed from Clewe O, Aulin L, Hu Y,",
    "Coates AR, Simonsson US. J Antimicrob Chemother 71(4):964-974 (2016)",
    "doi:10.1093/jac/dkv416; INH/EMB/PZA absorption rate constants fixed",
    "from Chen C, Ortega F, Alameda L, Ferrer S, Simonsson US.",
    "Eur J Pharm Sci 93:319-333 (2016) doi:10.1016/j.ejps.2016.07.014."
  )
  vignette <- "Chen_2017_TB_MTP_GPDI_mouse"
  units <- list(
    time          = "h",
    dosing        = "mg/kg",
    concentration = "log(CFU/lungs) for the model observation; mg/L for the internal drug plasma trajectories Cc_rif / Cc_inh / Cc_emb / Cc_pza"
  )

  covariateData <- list(
    DOSE_INH_MGKG = list(
      description        = "Per-subject assigned isoniazid dose level (mg/kg/day)",
      units              = "mg/kg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used to compute the dose-dependent isoniazid clearance:",
        "CL_inh = CL_inh_lowest * (1 - slope_inh * (DOSE_INH_MGKG - 12.5))",
        "(Chen 2017 Table 1 footnote b; lowest dose = 12.5 mg/kg).",
        "Paper-specific covariate not in inst/references/covariate-columns.md",
        "because the canonical DOSE entry is single-drug and this model carries",
        "four parallel drugs in combination; only INH requires a per-subject",
        "dose-level covariate (the other three drugs have linear / fixed-dose",
        "PK and their assigned dose appears only through the event AMT).",
        "Set to 0 in regimens without isoniazid."
      ),
      source_name        = "INH dose level (Chen 2017 Methods / Table 1)"
    )
  )

  population <- list(
    species        = "mouse (BALB/c, female)",
    n_subjects     = 49L,
    n_studies      = 2L,
    age_range      = "13-15 weeks at infection",
    weight_range   = "20-25 g",
    sex_female_pct = 100,
    disease_state  = paste(
      "M. tuberculosis Beijing VN 2002-1585 genotype intratracheal infection",
      "(approx 9.5e4 CFU inoculum) plus a healthy-mouse rifampicin PK sub-study (n = 18)."
    ),
    dose_range     = paste(
      "Daily oral gavage, 5 days/week. Monotherapy: rifampicin 5/10/20 mg/kg",
      "(plus healthy-mouse 10/160 mg/kg supporting cohort), isoniazid",
      "12.5/25/50 mg/kg, ethambutol 50/100/200 mg/kg, pyrazinamide",
      "75/150/300 mg/kg. Combination therapies (R10H25, R10H25Z150,",
      "R10H25Z150E100) at fixed doses for up to 24 weeks."
    ),
    regions        = "Charles River Les Oncins, France; experiments at Erasmus MC, Rotterdam, NL",
    notes          = paste(
      "TB-infected mice contributed sparse PK (one plasma sample per mouse",
      "at 1, 4, or 8 h post-dose after 4 weeks of treatment for RIF/INH or",
      "after 1 week of treatment for EMB/PZA -- the EMB and PZA monotherapy",
      "groups did not survive beyond 1 week). The rifampicin PK was further",
      "supported by a healthy-mouse PK sub-study (n = 18) at 10 and 160 mg/kg.",
      "CFU sampling per Chen 2017 Methods: 1/2/4 weeks of RIF or INH",
      "monotherapy (9 mice/time point, 3 per dose); 1 week of EMB or PZA",
      "monotherapy (6 mice); 1/2/4/8/12/24 weeks of combination therapy",
      "(3 mice/time point). Natural-growth CFU collected at 1, 3, 7, 14, 21",
      "days post-infection. Interindividual variability in PK was not",
      "quantified in the source (one sample per mouse); the model is a",
      "typical-value mechanistic PK + MTP + GPDI without etas."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All structural values from Chen 2017 Table 1 (popPK) and Table 2
    # (MTP-GPDI PD). Upstream-fixed values noted inline. The model is a
    # typical-value mechanistic mouse PK + MTP-GPDI -- the source paper
    # explicitly states "interindividual variability in PK was not
    # quantified" (Methods, Population pharmacokinetic study and
    # modeling section), and the PD layer uses NONMEM FOCE on the
    # pooled CFU dataset without subject-level etas. There are
    # therefore no eta parameters in ini().
    # ------------------------------------------------------------------

    # ===== Rifampicin PK (1-compartment, TB-infected mice, R5/R10/R20) =====
    lka_rif <- log(6.23)
    label("Rifampicin absorption rate constant ka (1/h)")
    # Chen 2017 Table 1: ka = 6.23, RSE 21.0%, TB-infected.

    lcl_rif <- log(234)
    label("Rifampicin apparent clearance CL/F in TB-infected mice (mL/h/kg)")
    # Chen 2017 Table 1: CL/F = 234 mL/h/kg, RSE 22.4%, TB-infected R5/R10/R20.

    lvc_rif <- log(706)
    label("Rifampicin apparent central volume V/F (mL/kg)")
    # Chen 2017 Table 1: V/F = 706 mL/kg, RSE 32.4%, TB-infected R5/R10/R20.

    # ===== Isoniazid PK (1-compartment, dose-dependent CL) =====
    lka_inh <- fixed(log(12.6))
    label("Isoniazid absorption rate constant ka (1/h, fixed from Chen 2016)")
    # Chen 2017 Table 1: ka = 12.6 FIX, fixed per Chen 2016 EJPS upstream popPK.

    lcl_inh_lowdose <- log(612)
    label("Isoniazid apparent CL/F at the lowest dose 12.5 mg/kg (mL/h/kg)")
    # Chen 2017 Table 1: CL/F = 612 mL/h/kg at H12.5, RSE 5.4%.

    slope_inh <- 7.50e-3
    label("Isoniazid CL/F dose-dependence slope (1/(mg/kg); Slope sigma)")
    # Chen 2017 Table 1: Slope sigma = 7.50e-3, RSE 16.8%.
    # CL_inh = lcl_inh_lowdose * (1 - slope_inh * (DOSE_INH_MGKG - 12.5)).

    lvc_inh <- log(811)
    label("Isoniazid apparent central volume V/F (mL/kg)")
    # Chen 2017 Table 1: V/F = 811 mL/kg, RSE 8.0%.

    # ===== Ethambutol PK (2-compartment) =====
    lka_emb <- fixed(log(0.87))
    label("Ethambutol absorption rate constant ka (1/h, fixed from Chen 2016)")
    # Chen 2017 Table 1: ka = 0.87 FIX, fixed per Chen 2016 EJPS upstream popPK.

    lcl_emb <- log(3400)
    label("Ethambutol apparent clearance CL/F (mL/h/kg)")
    # Chen 2017 Table 1: CL/F = 3400 mL/h/kg, RSE 11.0%.

    lvc_emb <- log(1500)
    label("Ethambutol apparent central volume V/F (mL/kg)")
    # Chen 2017 Table 1: V/F = 1500 mL/kg, RSE 22.3%.

    lq_emb <- log(2530)
    label("Ethambutol apparent intercompartmental clearance Q/F (mL/h/kg)")
    # Chen 2017 Table 1: Q/F = 2530 mL/h/kg, RSE 42.7%.

    lvp_emb <- log(4690)
    label("Ethambutol apparent peripheral volume V2/F (mL/kg)")
    # Chen 2017 Table 1: V2/F = 4690 mL/kg, RSE 27.3%.

    # ===== Pyrazinamide PK (1-compartment, linear) =====
    lka_pza <- fixed(log(2.84))
    label("Pyrazinamide absorption rate constant ka (1/h, fixed from Chen 2016)")
    # Chen 2017 Table 1: ka = 2.84 FIX, fixed per Chen 2016 EJPS upstream popPK.

    lcl_pza <- log(273.71)
    label("Pyrazinamide apparent clearance CL/F (mL/h/kg)")
    # Chen 2017 Table 1: CL/F = 273.71 mL/h/kg, RSE 15.9%.

    lvc_pza <- log(525.1)
    label("Pyrazinamide apparent central volume V/F (mL/kg)")
    # Chen 2017 Table 1: V/F = 525.1 mL/kg, RSE 22.1%.

    # ===== MTP natural-growth structure =====
    f0_init <- log(20100)
    label("log(F0) -- initial fast-multiplying bacterial number (CFU/lungs)")
    # Chen 2017 Table 2: F0 = 20100, RSE 59.2%.

    s0_init <- log(119000)
    label("log(S0) -- initial slow-multiplying bacterial number (CFU/lungs)")
    # Chen 2017 Table 2: S0 = 119000, RSE 31.0%.

    lkg <- log(0.034)
    label("kG -- fast-multiplying bacterial growth rate (1/h)")
    # Chen 2017 Table 2: kG = 0.034 1/h, RSE 10.0% (exponential growth function).

    lkfs_lin <- log(6.65e-5)
    label("kFS_lin -- time-dependent F->S transfer slope (1/h^2)")
    # Chen 2017 Table 2: kFS_lin = 6.65e-5, RSE 20.3%. kFS(t) = kFS_lin * t.

    lksf <- fixed(log(6.03e-4))
    label("kSF -- first-order S->F transfer (1/h, fixed from Clewe 2016)")
    # Chen 2017 Table 2 footnote a: kSF = 6.03e-4, FIXED to in vitro Clewe 2016.

    lkfn <- fixed(log(3.74e-8))
    label("kFN -- first-order F->N transfer (1/h, fixed from Clewe 2016)")
    # Chen 2017 Table 2 footnote a: kFN = 3.74e-8, FIXED to in vitro Clewe 2016.

    lksn <- fixed(log(7.73e-3))
    label("kSN -- first-order S->N transfer (1/h, fixed from Clewe 2016)")
    # Chen 2017 Table 2 footnote a: kSN = 7.73e-3, FIXED to in vitro Clewe 2016.

    lkns <- fixed(log(5.11e-5))
    label("kNS -- first-order N->S transfer (1/h, fixed from Clewe 2016)")
    # Chen 2017 Table 2 footnote a: kNS = 5.11e-5, FIXED to in vitro Clewe 2016.

    # ===== Rifampicin monotherapy drug effects =====
    lkg_r <- log(0.0022)
    label("kG_R -- linear RIF inhibition of F growth (mL/(h*ug))")
    # Chen 2017 Table 2: kG_R = 0.0022, RSE 29.2%. FG_R = kG_R * Cc_rif.

    loo_f_r <- log(0.019)
    label("OO_F_R -- on/off RIF stimulation of F death (1/h)")
    # Chen 2017 Table 2: OO_F_R = 0.019, RSE 4.2%. FD_R = OO_F_R when Cc_rif > 0.

    lks_r <- log(0.0137)
    label("kS_R -- linear RIF stimulation of S death (mL/(h*ug))")
    # Chen 2017 Table 2: kS_R = 0.0137, RSE 15.8%. SD_R = kS_R * Cc_rif.

    lkn_r <- log(0.0033)
    label("kN_R -- linear RIF stimulation of N death (mL/(h*ug))")
    # Chen 2017 Table 2: kN_R = 0.0033, RSE 8.0%. ND_R = kN_R * Cc_rif.

    # ===== Isoniazid monotherapy drug effects =====
    lkf_h <- log(0.055)
    label("kF_H -- linear INH stimulation of F death (mL/(h*ug))")
    # Chen 2017 Table 2: kF_H = 0.055, RSE 15.6%. FD_H = kF_H * Cc_inh.

    lks_h <- log(0.00047)
    label("kS_H -- linear INH stimulation of S death (mL/(h*ug))")
    # Chen 2017 Table 2: kS_H = 0.00047, RSE 53.6%. SD_H = kS_H * Cc_inh.

    lkn_h <- log(0.0012)
    label("kN_H -- linear INH stimulation of N death (mL/(h*ug))")
    # Chen 2017 Table 2: kN_H = 0.0012, RSE 5.0% (identified only in
    # R10H25Z150 combination; carried as a monotherapy slope for simulation
    # so the GPDI on/off combination effect can be applied uniformly).

    # ===== GPDI interactions (joint INT_AB = INT_BA, EC50_INT -> 0, on/off) =====
    int_s_rh <- 4.49
    label("INT_S_RH -- joint RIF-INH GPDI interaction on S death (fractional change of slope)")
    # Chen 2017 Table 2: INT_S_RH = 4.49, RSE 28.1%. Positive => antagonism
    # (decreased slope of both RIF and INH on the killing of S when both
    # drugs are present). Applied as slope_eff = slope / (1 + INT) per the
    # canonical Wicha 2017 GPDI formulation (see vignette Errata).

    int_n_rh <- 0.32
    label("INT_N_RH -- joint RIF-INH GPDI interaction on N death (fractional change of slope)")
    # Chen 2017 Table 2: INT_N_RH = 0.32, RSE 38.0%. Positive => antagonism
    # (slight decrease in slope of both RIF and INH on the killing of N).

    int_n_re <- -0.15
    label("INT_N_RE -- joint RIF-EMB GPDI interaction on N death (fractional change of slope)")
    # Chen 2017 Table 2: INT_N_RE = -0.15, RSE 57.7%. Negative => synergism
    # (increase of RIF efficacy slope on N in presence of EMB). EMB has no
    # quantified monotherapy effect (Discussion), so this interaction enters
    # via the GPDI multiplier on RIF only.

    # ===== Residual error (additive on natural-log CFU/lungs) =====
    # Chen 2017 Table 2: "Additive residual variability on log scale
    # (variance)" sigma^2 = 1.23, RSE 16.7%. NONMEM's "additive on log
    # scale" with SIGMA expressed on log10(CFU) variance translates to a
    # natural-log SD of sqrt(1.23) * log(10). Per parameter-names.md
    # "Residual error" section, encoded as addSd on the log-transformed
    # CFU observation. (If the paper variance is on natural log instead of
    # log10, the natural-log SD becomes sqrt(1.23); see vignette Errata.)
    addSd <- fixed(sqrt(1.23) * log(10))
    label("Additive residual SD on natural-log CFU/lungs (sqrt(SIGMA on log10) * log(10))")
  })

  model({
    # ================================================================
    # 1. Individual PK parameters (typical values; no etas per Methods)
    # ================================================================
    ka_rif <- exp(lka_rif)
    cl_rif <- exp(lcl_rif)
    vc_rif <- exp(lvc_rif)
    kel_rif <- cl_rif / vc_rif

    ka_inh <- exp(lka_inh)
    cl_inh_lowdose <- exp(lcl_inh_lowdose)
    vc_inh <- exp(lvc_inh)
    cl_inh <- cl_inh_lowdose * (1 - slope_inh * (DOSE_INH_MGKG - 12.5))
    # Chen 2017 Table 1 footnote b. cl_inh evaluates to ~612 at H12.5,
    # ~554 at H25, and ~440 at H50; the simulation reference R10H25 uses
    # DOSE_INH_MGKG = 25.
    kel_inh <- cl_inh / vc_inh

    ka_emb <- exp(lka_emb)
    cl_emb <- exp(lcl_emb)
    vc_emb <- exp(lvc_emb)
    q_emb  <- exp(lq_emb)
    vp_emb <- exp(lvp_emb)
    k12_emb <- q_emb / vc_emb
    k21_emb <- q_emb / vp_emb
    kel_emb <- cl_emb / vc_emb

    ka_pza <- exp(lka_pza)
    cl_pza <- exp(lcl_pza)
    vc_pza <- exp(lvc_pza)
    kel_pza <- cl_pza / vc_pza

    # ================================================================
    # 2. Drug PK ODEs -- four parallel sub-trees, each fed by its own
    #    cmt= dose event. The covariates DOSE_*_MGKG are informational;
    #    the actual administered amount comes through the event table.
    # ================================================================
    d/dt(depot_rif)   <- -ka_rif * depot_rif
    d/dt(central_rif) <-  ka_rif * depot_rif - kel_rif * central_rif

    d/dt(depot_inh)   <- -ka_inh * depot_inh
    d/dt(central_inh) <-  ka_inh * depot_inh - kel_inh * central_inh

    d/dt(depot_emb)     <- -ka_emb * depot_emb
    d/dt(central_emb)   <-  ka_emb * depot_emb -
                            (kel_emb + k12_emb) * central_emb +
                            k21_emb * peripheral_emb
    d/dt(peripheral_emb) <- k12_emb * central_emb - k21_emb * peripheral_emb

    d/dt(depot_pza)   <- -ka_pza * depot_pza
    d/dt(central_pza) <-  ka_pza * depot_pza - kel_pza * central_pza

    Cc_rif <- central_rif / vc_rif
    Cc_inh <- central_inh / vc_inh
    Cc_emb <- central_emb / vc_emb
    Cc_pza <- central_pza / vc_pza

    # ================================================================
    # 3. MTP natural-growth rates (back-out log scale)
    # ================================================================
    kg      <- exp(lkg)
    kfs_lin <- exp(lkfs_lin)
    ksf     <- exp(lksf)
    kfn     <- exp(lkfn)
    ksn     <- exp(lksn)
    kns     <- exp(lkns)
    kfs     <- kfs_lin * t   # time-dependent F->S transfer (Chen 2017 Methods)

    # ================================================================
    # 4. Drug effect slopes / on-off rates
    # ================================================================
    kg_r   <- exp(lkg_r)
    oo_f_r <- exp(loo_f_r)
    ks_r   <- exp(lks_r)
    kn_r   <- exp(lkn_r)

    kf_h <- exp(lkf_h)
    ks_h <- exp(lks_h)
    kn_h <- exp(lkn_h)

    # Drug-presence indicators -- 1 when the drug is currently delivering
    # measurable plasma concentration, 0 otherwise. Used to gate the
    # on/off GPDI interactions per Chen 2017 (EC50_INT set to 1e-8 in the
    # source, collapsing the Emax interaction term to an on/off switch
    # whenever the partner drug concentration is non-zero).
    rif_on <- 0
    if (Cc_rif > 1e-8) rif_on <- 1
    inh_on <- 0
    if (Cc_inh > 1e-8) inh_on <- 1
    emb_on <- 0
    if (Cc_emb > 1e-8) emb_on <- 1

    # GPDI multipliers on the RIF-INH slope (S and N states) and the
    # RIF-EMB slope (N state). slope_eff = slope_mono / (1 + INT) per the
    # canonical Wicha 2017 GPDI formulation; INT > 0 -> decreased potency
    # (antagonism), INT < 0 -> increased potency (synergism). The
    # multiplier is applied only when the partner drug is on; the
    # monotherapy slope is recovered when the partner is absent.
    mult_rh_s <- 1 + int_s_rh * inh_on
    mult_rh_s_h <- 1 + int_s_rh * rif_on
    mult_rh_n <- 1 + int_n_rh * inh_on
    mult_rh_n_h <- 1 + int_n_rh * rif_on
    mult_re_n <- 1 + int_n_re * emb_on

    # ================================================================
    # 5. Cumulative drug effects on each bacterial state
    # ================================================================
    # F-state: RIF inhibits growth (linear in Cc_rif) and adds an on/off
    # death term; INH adds a linear death term in Cc_inh.
    fg_drug <- kg_r * Cc_rif
    fd_drug <- oo_f_r * rif_on + kf_h * Cc_inh

    # S-state: RIF and INH share the GPDI interaction (joint INT_S_RH).
    sd_drug <- (ks_r * Cc_rif) / mult_rh_s + (ks_h * Cc_inh) / mult_rh_s_h

    # N-state: RIF effect is modulated by both INH (GPDI N-RH) and EMB
    # (GPDI N-RE); INH effect is modulated by RIF (GPDI N-RH).
    nd_drug <- (kn_r * Cc_rif) / (mult_rh_n * mult_re_n) +
               (kn_h * Cc_inh) / mult_rh_n_h

    # ================================================================
    # 6. MTP bacterial-state ODEs (Chen 2017 Eq. 10-12) with initial
    #    conditions from Chen 2017 Table 2.
    # ================================================================
    fbugs(0) <- exp(f0_init)
    sbugs(0) <- exp(s0_init)
    nbugs(0) <- 0
    # Chen 2017 Results: no inoculum of non-multiplying bacteria was
    # supported by the data, fixed to 0.

    d/dt(fbugs) <- (kg - fg_drug) * fbugs - kfs * fbugs + ksf * sbugs -
                   fd_drug * fbugs - kfn * fbugs
    d/dt(sbugs) <- kfs * fbugs - ksf * sbugs + kns * nbugs - ksn * sbugs -
                   sd_drug * sbugs
    d/dt(nbugs) <- ksn * sbugs + kfn * fbugs - kns * nbugs - nd_drug * nbugs

    # ================================================================
    # 7. Observation: log(F + S + N) in CFU/lungs, additive residual
    #    error on the natural-log scale (Chen 2017 Methods: "additive
    #    error model on log scale" with variance 1.23).
    # ================================================================
    total_bugs <- fbugs + sbugs + nbugs
    # Floor at LLOQ (10 CFU/lungs per Chen 2017 Methods) to avoid log of
    # zero when drug effect drives the integrated total to ~0. The single
    # model observation Cc is the natural-log of the total CFU/lungs
    # (canonical single-output name; the "Cc" symbol is reused here for the
    # bacterial output to satisfy the nlmixr2lib convention even though the
    # underlying quantity is a count, not a concentration).
    if (total_bugs < 10) total_bugs <- 10
    Cc <- log(total_bugs)
    Cc ~ add(addSd)
  })
}
