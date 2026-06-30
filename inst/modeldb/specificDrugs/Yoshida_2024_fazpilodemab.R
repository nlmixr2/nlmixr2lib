Yoshida_2024_fazpilodemab <- function() {
  description <- "Population PK + longitudinal gastrointestinal-adverse-event (GIAE) discrete-time Markov + treatment-discontinuation logistic model for the anti-FGFR/KLB bispecific antibody fazpilodemab (BFKB8488A) in adults with type 2 diabetes mellitus or non-alcoholic fatty liver disease (Yoshida 2024). The PK structure is a 2-compartment disposition with parallel direct (ka) plus transit-mediated (ka2 via the depot2 absorption compartment, transit rate ktr) subcutaneous absorption; linear plus Michaelis-Menten elimination acting on free drug; target-mediated quasi-steady-state binding to a constant target concentration rmax; and a sigmoidal time-onset ADA-mediated clearance arm cl_ada that activates when ADA_POS = 1. The DTMM transition probabilities (grade 0 / grade 1 / grade 2-3) and the treatment-discontinuation probability are emitted as algebraic outputs conditional on the time-varying PREV_AE_SCORE covariate (the previous-day GIAE grade); the actual Markov-chain stochastic simulation is performed downstream in R using these probabilities (matching the mrgsolve plus Rcpp framework of the source paper supplement S2.3)."
  reference <- "Yoshida K, Poon V, Dash A, Kunder R, Chinn L, Kagedal M. Simulation-based evaluation of personalized dosing approaches for anti-FGFR/KLB bispecific antibody fazpilodemab. CPT Pharmacometrics Syst Pharmacol. 2024;13(4):544-550. doi:10.1002/psp4.13111"
  vignette <- "Yoshida_2024_fazpilodemab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    ADA_POS = list(
      description        = "Anti-drug antibody (ADA) positivity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Time-fixed per subject in the simulation framework. When ADA_POS = 1, an additional sigmoidal time-onset ADA-mediated clearance arm cl_ada activates (cl_ada(t) = cl_ada_max / (1 + exp(-k_ada * (t - t50_ada)))); when ADA_POS = 0 the arm is identically zero. Source column ADA in mrgsolve supplement S2.3.1 (zero default in [PARAM]).",
      source_name        = "ADA"
    ),
    PREV_AE_SCORE = list(
      description        = "Previous-day gastrointestinal AE (GIAE) grade",
      units              = "(ordinal score 0..2; 0 = no GIAE, 1 = grade 1, 2 = combined grade 2/3)",
      type               = "count",
      reference_category = "0 (no GIAE)",
      notes              = paste(
        "Time-varying (updated each day from the previous day's sampled grade). Used to condition the DTMM transition logits (B1, B2) and the discontinuation logit (B_dc) on the previous Markov state.",
        "Fazpilodemab pools grades 2 and 3 into a single category (PREV_AE_SCORE = 2) because of low event counts at grade 3 (Methods page 545).",
        "When constructing a simulation event table, set PREV_AE_SCORE = 0 at the first observation of every subject and update each subsequent observation to the previous observation's sampled grade -- matching the mrgsolve [TABLE] block carry-forward in supplement S2.3.1 (GR variable).",
        "Source column PDV in NONMEM control stream S2.2 (carried forward via IF (EVID.EQ.0.AND.TYPE.EQ.1) PDV=DV).",
        sep = " "
      ),
      source_name        = "PDV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 153L,
    n_studies      = 1L,
    age_range      = "Adults (specific age range not reported in main text or available supplement)",
    weight_range   = "Adults (specific weight range not reported in main text or available supplement)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in main text or available supplement",
    disease_state  = "Type 2 diabetes mellitus (T2DM) or non-alcoholic fatty liver disease (NAFLD)",
    dose_range     = "10-250 mg subcutaneous q1w, q2w, or q4w",
    regions        = "Not reported in main text or available supplement",
    notes          = paste(
      "Multiple ascending dose (MAD) study GC39547 (NCT03060538). 121 patients with T2DM or NAFLD received fazpilodemab and 32 patients received placebo (n = 153 total).",
      "Fazpilodemab was administered subcutaneously in the abdomen or thigh at 10-250 mg with intervals of q1w, q2w, or q4w.",
      "Detailed demographic breakdown (age / weight / sex / race) is not reported in the main text or in the available supplements; the available trimmed-markdown copies of Supplements 1 (figures) and 2 (model code) contain neither a Table 1 baseline-demographics summary nor an extended-text demographics paragraph.",
      "Phase I/II population characterisation; full popPK described as 'unpublished data' in the main text, with the typical-value parameter set and IIV/IOV/residual error provided in full as the mrgsolve [PARAM] / [OMEGA] / [SIGMA] block of Supplement S2.3.1.",
      sep = " "
    )
  )

  ini({
    # ---------- PK structural parameters ----------
    # All typical values from mrgsolve [PARAM] block of Supplement S2.3.1.
    # Source units there are ug / L / day; this model is in mg / L / day,
    # so VMAX (mass / time) and KM / RMAX / KSS (concentration) are
    # rescaled by 1e-3 (1 mg = 1000 ug). All others are unit-invariant.
    lka     <- log(0.243);   label("Direct absorption rate ka, depot -> central (1/day)")               # S2.3.1 TVKA    = 0.243
    lktr    <- log(0.115);   label("Transit rate ktr, depot -> depot2 (1/day)")                          # S2.3.1 TVKTR   = 0.115
    lka2    <- log(0.0687);  label("Secondary absorption rate ka2, depot2 -> central (1/day)")           # S2.3.1 TVKA2   = 0.0687
    lcl     <- log(0.349);   label("Linear clearance on free drug, cl (L/day)")                          # S2.3.1 TVCL    = 0.349
    lvc     <- log(2.98);    label("Central volume of distribution, vc (L)")                             # S2.3.1 TVV2    = 2.98
    lvp     <- log(1.83);    label("Peripheral volume of distribution, vp (L)")                          # S2.3.1 TVV3    = 1.83
    lq      <- log(3.07);    label("Intercompartmental clearance on free drug, q (L/day)")               # S2.3.1 TVQ     = 3.07
    lvmax   <- log(3.34);    label("Michaelis-Menten Vmax for saturable elimination on free drug (mg/day)")  # S2.3.1 TVVMAX = 3340 ug/day = 3.34 mg/day
    lkm     <- log(0.482);   label("Michaelis-Menten Km for saturable elimination on free drug (ug/mL)")  # S2.3.1 TVKM   = 482 ug/L = 0.482 ug/mL
    lrmax   <- fixed(log(2));     label("Steady-state target receptor concentration, rmax (ug/mL)")       # S2.3.1 TVRMAX = 2000 ug/L = 2 ug/mL; fixed structural constant for QSS-TMDD with no target turnover
    lkss    <- log(0.0885);  label("Quasi-steady-state dissociation constant, kss (ug/mL)")              # S2.3.1 TVKSS  = 88.5 ug/L = 0.0885 ug/mL
    lfdepot <- log(0.719);   label("Subcutaneous bioavailability F1 (fraction)")                          # S2.3.1 TVF1    = 0.719

    # ---------- ADA-mediated clearance arm (cl_ada family) ----------
    # When ADA_POS = 1, cl_ada(t) = cl_ada_max / (1 + exp(-k_ada * (t - t50_ada))).
    # The arm is identically zero when ADA_POS = 0.
    lcl_ada   <- log(3.39);  label("Maximum ADA-mediated clearance on free drug, cl_ada_max (L/day)")     # S2.3.1 TVCLADAMAX = 3.39
    lt50_ada  <- log(46.6);  label("Time at half-maximum ADA-mediated clearance, t50_ada (day)")          # S2.3.1 TVT50ADA   = 46.6
    lk_ada    <- log(0.195); label("Rate constant for sigmoidal ADA-effect onset, k_ada (1/day)")         # S2.3.1 TVKADA     = 0.195

    # ---------- Longitudinal GIAE discrete-time Markov model (DTMM) ----------
    # All typical values from mrgsolve [PARAM] block of Supplement S2.3.1
    # (these are the final estimates; the NONMEM control stream in S2.2 lists
    # the corresponding initial values). LGT1 = B1 + slp_ae * Cc, LGT2 = B2 +
    # slp_ae * Cc (with Cc in ug/mL); B1, B2 depend on PREV_AE_SCORE.
    # Time-effect term applies only when PREV_AE_SCORE = 0 and is capped at
    # time = 84 days (mrgsolve uses min(TIME + TIMEinit, 84.0)).
    lslp_ae     <- log(0.0927); label("Log of drug-effect slope on GIAE transition logit (per ug/mL)")    # S2.3.1 TVSLP    = 0.0927; matches main-text "9.72% increase in odds per ug/mL" = log(1.0972); IIV is exponential per mrgsolve SLP = TVSLP * exp(ETA_SLP)
    b1_g0_ae    <- -4.69;    label("DTMM intercept B1 conditional on PREV_AE_SCORE = 0 (logit grade >= 1)") # S2.3.1 TV0B1    = -4.69
    b2b1_g0_ae  <- -0.799;   label("DTMM offset B2-B1 conditional on PREV_AE_SCORE = 0 (logit grade >= 2 minus logit grade >= 1)") # S2.3.1 TV0B2B1 = -0.799
    b1_g1_ae    <-  2.4;     label("DTMM intercept B1 conditional on PREV_AE_SCORE = 1 (logit grade >= 1)") # S2.3.1 TV1B1    = 2.4
    b2b1_g1_ae  <- -14.8;    label("DTMM offset B2-B1 conditional on PREV_AE_SCORE = 1 (logit grade >= 2 minus logit grade >= 1)") # S2.3.1 TV1B2B1 = -14.8
    b1_g2_ae    <-  1.83;    label("DTMM intercept B1 conditional on PREV_AE_SCORE = 2 (logit grade >= 1)") # S2.3.1 TV2B1    = 1.83
    b2b1_g2_ae  <- -0.037;   label("DTMM offset B2-B1 conditional on PREV_AE_SCORE = 2 (logit grade >= 2 minus logit grade >= 1)") # S2.3.1 TV2B2B1 = -0.037
    tef_g0_ae   <- -0.0336;  label("Linear time effect on DTMM intercept when PREV_AE_SCORE = 0, capped at 84 days (per day)") # S2.3.1 TIME_EFF = -0.0336

    # ---------- Treatment-discontinuation logistic model ----------
    # logit(P_discon) = b_dc with b_dc depending on PREV_AE_SCORE.
    # In mrgsolve the discontinuation is evaluated once per dosing cycle on
    # day 14, gated to PREV_AE_SCORE > 0 (no discontinuation from grade 0).
    b_dc_g0     <- -4.31;    label("Discontinuation logit intercept conditional on PREV_AE_SCORE = 0")    # S2.3.1 B_discon_G0 = -4.31
    b_dc_g1     <- -1.69;    label("Discontinuation logit intercept conditional on PREV_AE_SCORE = 1")    # S2.3.1 B_discon_G1 = -1.69 (= -4.31 + 2.62)
    b_dc_g2     <- -1.33;    label("Discontinuation logit intercept conditional on PREV_AE_SCORE = 2")    # S2.3.1 B_discon_G2 = -1.33 (= -4.31 + 2.98)

    # ---------- Inter-individual variability ----------
    # All variances on the exponential / log scale (omega^2) from the mrgsolve
    # [OMEGA] @name IIV block of Supplement S2.3.1.
    etalcl       ~ 0.0174    # S2.3.1 E_CL       = 0.0174 (omega^2 on log CL)
    etalfdepot   ~ 0.0952    # S2.3.1 E_F1       = 0.0952 (omega^2 on log F1; note IOV on F1 of 0.0341 not implemented -- see vignette Errata)
    etalvc       ~ 0.286     # S2.3.1 E_V2       = 0.286
    etalka2      ~ 0.297     # S2.3.1 E_KA2      = 0.297
    etalcl_ada   ~ 5.21      # S2.3.1 E_CLADAMAX = 5.21 (very high; reflects bimodal ADA-response distribution)
    etalt50_ada  ~ 0.155     # S2.3.1 E_T50ADA   = 0.155
    etalslp_ae   ~ 1.85      # S2.3.1 E_SLP      = 1.85 (omega^2 on log slope of drug effect on AE)

    # ---------- Residual error ----------
    # mrgsolve [SIGMA] EPSP = 0.0437 PROPORTIONAL ERROR (variance). The
    # additive component (EPSA = 0) is omitted because it is identically zero.
    propSd       <- sqrt(0.0437); label("Proportional residual standard deviation (fraction)") # S2.3.1 EPSP = 0.0437 (variance); SD = sqrt(0.0437) ~ 0.209
  })

  model({
    # ---------- Individual PK parameters ----------
    ka      <- exp(lka)
    ktr     <- exp(lktr)
    ka2     <- exp(lka2 + etalka2)
    cl      <- exp(lcl + etalcl)
    vc      <- exp(lvc + etalvc)
    vp      <- exp(lvp)
    q       <- exp(lq)
    vmax    <- exp(lvmax)
    km      <- exp(lkm)
    rmax    <- exp(lrmax)
    kss     <- exp(lkss)
    fdepot  <- exp(lfdepot + etalfdepot)

    # ---------- ADA-mediated clearance (sigmoidal time-onset) ----------
    cl_ada_max <- exp(lcl_ada  + etalcl_ada)
    t50_ada    <- exp(lt50_ada + etalt50_ada)
    k_ada      <- exp(lk_ada)
    cl_ada     <- ADA_POS * cl_ada_max / (1 + exp(-k_ada * (time - t50_ada)))

    # ---------- QSS-TMDD with constant target concentration rmax ----------
    # ctot = total drug concentration in central compartment.
    # The QSS quadratic (Gibiansky 2008 Eq. 7 with constant target) yields
    # free drug concentration cfree; reported observation Cc = cfree because
    # the bioanalytical assay is documented as measuring free fazpilodemab
    # in the source paper's mrgsolve definition (IPRED = CFREE).
    ctot  <- central / vc
    disc  <- ctot - rmax - kss
    cfree <- 0.5 * (disc + sqrt(disc * disc + 4 * kss * ctot))
    Cc    <- cfree

    # ---------- ODE system (mass in mg; volumes in L; time in day) ----------
    # depot (A1): SC absorption depot; loses to central via ka and to
    #             depot2 via ktr (parallel direct + transit absorption).
    # depot2 (A4): secondary absorption compartment; gains from depot via
    #              ktr and loses to central via ka2.
    # central (A2): receives ka * depot + ka2 * depot2 + q * peripheral1 / vp;
    #               loses (cl + cl_ada + q) * cfree (linear + ADA + distribution
    #               all on free drug) and the Michaelis-Menten elimination term
    #               vmax * cfree / (km + cfree).
    # peripheral1 (A3): symmetric q * (cfree - peripheral1 / vp); free drug
    #                   exchanges with the peripheral compartment.
    d/dt(depot)       <- -ka * depot - ktr * depot
    d/dt(depot2)      <-  ktr * depot - ka2 * depot2
    d/dt(central)     <-  ka * depot + ka2 * depot2 + q * peripheral1 / vp -
                          (cl + cl_ada + q) * cfree -
                          vmax * cfree / (km + cfree)
    d/dt(peripheral1) <-  q * (cfree - peripheral1 / vp)

    f(depot) <- fdepot

    # ---------- DTMM transition probabilities (algebraic outputs) ----------
    # The previous-day AE grade conditions B1 / B2; the time effect applies
    # only when PREV_AE_SCORE = 0 and is capped at day 84 of treatment per
    # the mrgsolve [TABLE] block min(TIME + TIMEinit, 84.0).
    tef_now <- tef_g0_ae * min(time, 84.0)
    b1_dtmm <- (PREV_AE_SCORE == 0) * (b1_g0_ae + tef_now) +
               (PREV_AE_SCORE == 1) *  b1_g1_ae +
               (PREV_AE_SCORE == 2) *  b1_g2_ae
    b2_dtmm <- b1_dtmm +
               (PREV_AE_SCORE == 0) * b2b1_g0_ae +
               (PREV_AE_SCORE == 1) * b2b1_g1_ae +
               (PREV_AE_SCORE == 2) * b2b1_g2_ae

    slp_ae   <- exp(lslp_ae + etalslp_ae)
    lgt1     <- b1_dtmm + slp_ae * Cc
    lgt2     <- b2_dtmm + slp_ae * Cc

    pge1 <- exp(lgt1) / (1 + exp(lgt1))  # Pr(next grade >= 1 | PREV_AE_SCORE, Cc)
    pge2 <- exp(lgt2) / (1 + exp(lgt2))  # Pr(next grade >= 2 | PREV_AE_SCORE, Cc)
    p_ae_g0   <- 1    - pge1             # Pr(next grade = 0)
    p_ae_g1   <- pge1 - pge2             # Pr(next grade = 1)
    p_ae_g23  <- pge2                    # Pr(next grade in {2, 3})

    # ---------- Discontinuation probability (algebraic output) ----------
    # logit(P_discon) intercept conditional on PREV_AE_SCORE. In the source
    # mrgsolve workflow the discontinuation is evaluated once per dosing
    # cycle on day 14 and is gated to PREV_AE_SCORE > 0 (no discontinuation
    # from grade 0); the gate is enforced downstream in the simulator rather
    # than in the algebraic output here.
    b_dc <- (PREV_AE_SCORE == 0) * b_dc_g0 +
            (PREV_AE_SCORE == 1) * b_dc_g1 +
            (PREV_AE_SCORE == 2) * b_dc_g2
    p_dc <- exp(b_dc) / (1 + exp(b_dc))

    # ---------- Observation and residual error ----------
    # Observation is free fazpilodemab concentration in central (ug/mL).
    Cc ~ prop(propSd)
  })
}
