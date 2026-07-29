Iida_2008_nicorandil <- function() {
  description <- paste(
    "Two-compartment IV PK plus inhibitory Emax PD plus asymptotic-",
    "exponential disease-progression population PKPD model for",
    "nicorandil in acute heart failure (AHF) patients with pulmonary",
    "artery wedge pressure (PAWP, mmHg) as the haemodynamic biomarker",
    "(Iida 2008 model 10). Five clinical studies pooled: 11 healthy",
    "volunteer subjects (concentration only) and 94 AHF patients",
    "(concentration plus PAWP), 618 nicorandil and 559 PAWP",
    "observations. Allometric size scaling with a 70 kg reference",
    "(exponent 0.75 on CL and Q per paper Equation 5; canonical",
    "exponent 1.0 on V1 and V2 per Holford 1996). AHF-vs-healthy",
    "disease cohort modifies all four PK parameters multiplicatively",
    "(FCL = 1.94, FV1 = 1.39, FQ = 0.519, FV2 = 4.06; paper Table 2).",
    "PD layer: inhibitory Emax model on plasma concentration (no",
    "effect compartment, paper Results paragraph 'Pharmacodynamic",
    "analysis'), summed with an asymptotic-exponential disease-",
    "progression term that decreases PAWP from a baseline S0 to a",
    "steady-state Sss with half-life Tprog. Inhibitory maximum",
    "Emax = -11.7 mmHg, EC50 = 423 ng/mL, S0 = 25.6 mmHg, Sss = 19.5",
    "mmHg, Tprog = 5.83 h. Bootstrap median final estimates."
  )
  reference <- paste(
    "Iida S, Kinoshita H, Holford NHG. Population pharmacokinetic and",
    "pharmacodynamic modelling of the effects of nicorandil in the",
    "treatment of acute heart failure. Br J Clin Pharmacol.",
    "2008;66(3):352-365. doi:10.1111/j.1365-2125.2008.03257.x.",
    sep = " "
  )
  vignette <- "Iida_2008_nicorandil"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed body weight used for allometric scaling on a 70 kg",
        "reference: exponent 0.75 on CL and Q per paper Equation 5",
        "(CL_GRP = CL_POPSTD * (Wt / 70)^(3/4)); exponent 1.0 on V1 and",
        "V2 per the Holford 1996 size standard cited as reference [11]",
        "in the paper. Cohort weight (mean +/- SD) per Table 1: 62.6",
        "+/- 8.6 kg (healthy bolus), 63.6 +/- 6.4 kg (healthy bolus+inf),",
        "58.6 +/- 11.3 kg (AHF bolus), 60.3 +/- 10.7 kg (AHF bolus+inf),",
        "61.8 +/- 10.2 kg (AHF long-term infusion)."
      ),
      source_name        = "Weight (kg)"
    ),
    DIS_HEALTHY = list(
      description        = paste(
        "Healthy-volunteer cohort indicator (1 = healthy adult, 0 =",
        "acute heart failure (AHF) patient). Time-fixed per subject."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (AHF patient)",
      notes              = paste(
        "Iida 2008 reports POP_CL, POP_V1, POP_Q, POP_V2 in Table 2 as",
        "healthy-volunteer typicals with multiplicative fractional-",
        "change factors FCL = 1.94, FV1 = 1.39, FQ = 0.519, FV2 = 4.06",
        "encoding the AHF-vs-healthy contrast (paper Equation 6:",
        "CL_POP_patients = CL_POP_HV * FCL). The model file re-",
        "anchors the typical-value parameters on the AHF reference",
        "state (DIS_HEALTHY = 0) so the canonical DIS_HEALTHY",
        "orientation (0 = patient, 1 = healthy) applies: lcl, lvc, lq,",
        "lvp encode the AHF-typical value (POP_X * F_X), and the four",
        "e_dis_healthy_<pk> log-shifts equal log(1 / F_X), recovering",
        "the published healthy-volunteer typical when DIS_HEALTHY = 1.",
        "This re-expression follows the Bulitta_2010_ceftazidime.R",
        "precedent for the DIS_HEALTHY canonical when the paper's",
        "structural reference is the healthy cohort but the canonical",
        "reference category is the patient cohort. Only the PK",
        "parameters carry a DIS_HEALTHY effect; PD parameters",
        "(S0, Sss, Tprog, Emax, EC50) were estimated from AHF data",
        "only and have no healthy-vs-patient contrast in the source.",
        "The 11 healthy volunteers in the pooled dataset contributed",
        "only nicorandil concentrations; no PAWP was measured in",
        "healthy subjects."
      ),
      source_name        = "(paper text: 'healthy subjects' vs 'AHF patients'; Tables 1 and 2)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 105L,
    n_studies      = 5L,
    age_range      = paste(
      "Healthy bolus 33.0 +/- 7.7 y; healthy bolus+inf 41.2 +/- 5.9 y;",
      "AHF bolus 62.8 +/- 12.8 y; AHF bolus+inf 65.1 +/- 11.5 y;",
      "AHF long-term inf 70.8 +/- 8.4 y (Table 1, mean +/- SD)"
    ),
    age_median     = "(not reported as median; mean +/- SD per Table 1)",
    weight_range   = "58.6-63.6 kg cohort means (Table 1; individual range not tabulated)",
    weight_median  = "(not reported as median; mean +/- SD per Table 1)",
    sex_female_pct = 32.4,
    race_ethnicity = "(not reported in Table 1; conducted in Japan)",
    disease_state  = paste(
      "Acute heart failure (AHF) including acute exacerbation of",
      "chronic heart failure, NYHA functional class II to IV. AHF",
      "aetiologies pooled across studies (Table 1): ischaemic heart",
      "disease (33), valvular heart disease (22), hypertensive heart",
      "disease (13), dilated cardiomyopathy (22). Inclusion criteria:",
      "baseline PAWP (or diastolic pulmonary arterial pressure if",
      "PAWP was technically not measurable) > 18 mmHg in the bolus",
      "and bolus+infusion AHF studies; > 15 mmHg in the long-term",
      "infusion study. Healthy volunteers were enrolled in two of the",
      "five studies and contributed nicorandil concentrations only",
      "(no PAWP)."
    ),
    dose_range     = paste(
      "Healthy bolus: 6, 12, 18, 24 mg IV over 5 min (every subject",
      "received four ascending doses). Healthy bolus+infusion: 12 mg",
      "IV bolus over 5 min followed by 6, 9, or 12 mg/h infusion for",
      "235 min. AHF bolus: 4, 8, 12, or 18 mg IV bolus over 1-5 min.",
      "AHF bolus+infusion: 200 ug/kg IV bolus over 5 min followed by",
      "50, 100, 150, 200, or 250 ug/kg/h infusion for 6 or 24 h.",
      "AHF long-term infusion: 200 ug/kg IV bolus over 5 min followed",
      "by 200 ug/kg/h infusion for 48 h."
    ),
    regions        = "Japan (Chugai Clinical Research Center sponsor)",
    notes          = paste(
      "Population-PK estimation: NONMEM Version V with FOCE-INTER",
      "(Methods, 'Model estimation'). The PK model was developed",
      "first; the PD layer was developed with PK parameters held",
      "fixed (PPPD method of Zhang, Beal, Sheiner 2003), then the",
      "final PKPD model was refit simultaneously for the bootstrap",
      "runs (model 10 in Table 3, 'Simultaneous' estimation, OFV",
      "7600.48). Parameter uncertainty by non-parametric bootstrap",
      "(Table 2 / Table 4 'Bootstrap CV' and '5% / 95% CI' columns).",
      "Effect-compartment model was tested but rejected (paper",
      "Results, 'Pharmacodynamic analysis': estimated effect-",
      "compartment half-life 2.2 min with very small variability; the",
      "NONMEM run finished with rounding errors and less than three",
      "significant digits and an objective function fall of 8.9; the",
      "authors concluded there was little evidence to support the",
      "additional complexity)."
    )
  )

  ini({
    # ====================================================================
    # PK structural parameters -- Iida 2008 Table 2 final-model (model 10)
    # bootstrap-median estimates.
    #
    # Iida 2008 reports POP_CL, POP_V1, POP_Q, POP_V2 in Table 2 as
    # healthy-volunteer typical values at WT = 70 kg with multiplicative
    # AHF/healthy factors FCL = 1.94, FV1 = 1.39, FQ = 0.519, FV2 = 4.06
    # (paper Equation 6: CL_POP_patients = CL_POP_HV * FCL). To honour the
    # canonical DIS_HEALTHY orientation (0 = patient reference, 1 =
    # healthy), the typical-value parameters below are anchored on the
    # AHF cohort (DIS_HEALTHY = 0) by multiplying the published healthy
    # typicals by the corresponding F factor. The four e_dis_healthy_<pk>
    # log-shifts then equal log(1 / F_X) and recover the published
    # healthy-volunteer typicals when DIS_HEALTHY = 1. This is the
    # standing re-expression pattern documented on the DIS_HEALTHY entry
    # of inst/references/covariate-columns.md and used in
    # Bulitta_2010_ceftazidime.R, Yoneyama_2017_emicizumab.R, and
    # Goggin_2004_emfilermin.R (among others).
    # ====================================================================
    lcl <- log(26.3 * 1.94)
    label("Apparent clearance CL at WT=70 kg, AHF reference (L/h)")           # Table 2: POP_CL = 26.3 L/h/70kg (healthy); FCL = 1.94
    lvc <- log(18.1 * 1.39)
    label("Central volume of distribution V1 at WT=70 kg, AHF reference (L)") # Table 2: POP_V1 = 18.1 L/70kg (healthy); FV1 = 1.39
    lq  <- log(71.6 * 0.519)
    label("Intercompartmental clearance Q at WT=70 kg, AHF reference (L/h)")  # Table 2: POP_Q = 71.6 L/h/70kg (healthy); FQ = 0.519
    lvp <- log(24.1 * 4.06)
    label("Peripheral volume of distribution V2 at WT=70 kg, AHF reference (L)") # Table 2: POP_V2 = 24.1 L/70kg (healthy); FV2 = 4.06

    # AHF-vs-healthy disease cohort shifts on PK parameters (multiplicative
    # in linear space; log(1 / F_X) so DIS_HEALTHY = 1 multiplies the AHF
    # typical by 1 / F_X to recover the published healthy typical).
    e_dis_healthy_cl <- log(1 / 1.94)
    label("Healthy-vs-AHF effect on CL (log scale)")                          # Table 2: FCL = 1.94; log(1/1.94) ~ -0.6627
    e_dis_healthy_vc <- log(1 / 1.39)
    label("Healthy-vs-AHF effect on V1 (log scale)")                          # Table 2: FV1 = 1.39; log(1/1.39) ~ -0.3293
    e_dis_healthy_q  <- log(1 / 0.519)
    label("Healthy-vs-AHF effect on Q (log scale)")                           # Table 2: FQ = 0.519; log(1/0.519) ~ +0.6553
    e_dis_healthy_vp <- log(1 / 4.06)
    label("Healthy-vs-AHF effect on V2 (log scale)")                          # Table 2: FV2 = 4.06; log(1/4.06) ~ -1.4015

    # Allometric exponents on the 70 kg reference (paper Methods, 'Covariate
    # effects': Equation 5 with exponent 3/4 for CL; Holford 1996 size
    # standard cited as reference [11] supplies the canonical 1.0 for V).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL (unitless, FIXED)")                      # Paper Eq 5: (Wt / 70)^(3/4)
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on V1 (unitless, FIXED)")                      # Holford 1996 [11] canonical exponent on V
    e_wt_q  <- fixed(0.75)
    label("Allometric exponent on Q (unitless, FIXED)")                       # Paper Eq 5: (Wt / 70)^(3/4) applied to all clearance terms
    e_wt_vp <- fixed(1.0)
    label("Allometric exponent on V2 (unitless, FIXED)")                      # Holford 1996 [11] canonical exponent on V

    # ====================================================================
    # PD structural parameters -- Iida 2008 Table 2 final-model (model 10)
    # bootstrap-median estimates. PD layer was estimated from AHF patients
    # only; healthy subjects contributed no PAWP observations.
    # ====================================================================
    lrbase <- log(25.6)
    label("Baseline PAWP S0 at start of observation (mmHg)")                  # Table 2: POP_S0 = 25.6 mmHg
    lsss   <- log(19.5)
    label("Asymptotic (long-time) PAWP Sss after disease progression (mmHg)") # Table 2: POP_Sss = 19.5 mmHg
    ltprog <- log(5.83)
    label("Disease-progression half-life Tprog (h)")                          # Table 2: POP_TPROG = 5.83 h
    limax  <- log(11.7)
    label("Magnitude of inhibitory Emax of nicorandil on PAWP (mmHg)")        # Table 2: POP_Emax = -11.7 mmHg (magnitude 11.7; negative sign explicit in model())
    lec50  <- log(423)
    label("EC50 of nicorandil on PAWP (ng/mL)")                               # Table 2: POP_EC50 = 423 ug/L = 423 ng/mL

    # ====================================================================
    # PK between-subject variability block -- Iida 2008 Table 4 random-
    # effects final-model (model 10) bootstrap-median estimates.
    #
    # The paper text (Methods, 'Random effects') describes the random
    # effects as N(0, BSV^2) where BSV is the standard deviation, then
    # states 'The variance of BSV and WSV was estimated.' Table 4 values
    # are therefore taken as variances (omega^2) per the NONMEM $OMEGA
    # convention. Reported off-diagonals in Table 4: R78 (CL-V1) = 0.257,
    # R89 (V1-Q) = 0.957. R79 (CL-Q) is not reported and is fixed to 0 in
    # the block below. No BSV reported for V2 (no etalvp).
    #
    # Covariance derivations from Table 4 variances and correlations:
    # cov(CL, V1) = R78 * sqrt(var_CL * var_V1)
    #            = 0.257 * sqrt(0.686 * 0.301) = 0.1168
    # cov(V1, Q) = R89 * sqrt(var_V1 * var_Q)
    #            = 0.957 * sqrt(0.301 * 0.414) = 0.3378
    # ====================================================================
    # Lower-triangular row-major: var_CL; cov_CL_V1, var_V1; cov_CL_Q (=0), cov_V1_Q, var_Q.
    etalcl + etalvc + etalq ~ c(
      0.686,
      0.1168,    0.301,
      fixed(0),  0.3378,    0.414
    )

    # ====================================================================
    # PD between-subject variability block -- Iida 2008 Table 4 final-
    # model (model 10) random-effects variances on the log-transformed
    # S0, Sss, EC50. Reported off-diagonals in Table 4: R12 (S0-Sss) =
    # 0.864, R23 (Sss-EC50) = 0.054. R13 (S0-EC50) is not reported and is
    # fixed to 0 in the block below. No PPV reported for Tprog or Imax.
    #
    # Covariance derivations from Table 4 variances and correlations:
    # cov(S0, Sss)   = R12 * sqrt(var_S0 * var_Sss)
    #                = 0.864 * sqrt(0.210 * 0.320) = 0.2240
    # cov(Sss, EC50) = R23 * sqrt(var_Sss * var_EC50)
    #                = 0.054 * sqrt(0.320 * 0.879) = 0.02864
    # ====================================================================
    # Lower-triangular row-major: var_S0; cov_S0_Sss, var_Sss; cov_S0_EC50 (=0), cov_Sss_EC50, var_EC50.
    etalrbase + etalsss + etalec50 ~ c(
      0.210,
      0.2240,    0.320,
      fixed(0),  0.02864,   0.879
    )

    # ====================================================================
    # Residual error -- Iida 2008 Table 4 final-model (model 10).
    # Nicorandil concentration (Cc): combined proportional + additive
    # (paper Equation 8: SD = sqrt((PRED * RUV_CV)^2 + RUV_SD^2)).
    # PAWP: additive only.
    #
    # Iida 2008 also estimated a between-subject random scaling (eta on
    # the residual-error magnitude; paper Equation 8: SD * exp(eta_PPV_RUV),
    # citing Karlsson et al. 1998 [14]) with magnitudes PPV_RUVCP = 0.716
    # and PPV_RUVFX = 0.308 (Table 4). nlmixr2 has no native eta-on-sigma
    # parameterisation, so the per-subject random scaling is omitted in
    # this implementation; the typical residual-error magnitudes
    # (RUV_CVCP, RUV_SDCP, RUV_SDFX) are used. See vignette Errata.
    # ====================================================================
    propSd     <- 0.212
    label("Proportional residual error on nicorandil concentration (fraction)") # Table 4: RUV_CVCP = 0.212 (21.2% CV)
    addSd      <- 6.82
    label("Additive residual error on nicorandil concentration (ng/mL)")      # Table 4: RUV_SDCP = 6.82 ng/mL
    addSd_pawp <- 2.5
    label("Additive residual error on PAWP (mmHg)")                           # Table 4: RUV_SDFX = 2.5 mmHg
  })

  model({
    # ------------------------------------------------------------------
    # Individual PK parameters at the AHF reference state (DIS_HEALTHY =
    # 0). DIS_HEALTHY = 1 multiplies the AHF typical by exp(e_dis_healthy
    # _<pk>) = 1 / F_<pk>, recovering the published healthy-volunteer
    # typical. Allometric size scaling on a 70 kg reference.
    # ------------------------------------------------------------------
    cl <- exp(lcl + e_dis_healthy_cl * DIS_HEALTHY + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + e_dis_healthy_vc * DIS_HEALTHY + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + e_dis_healthy_q  * DIS_HEALTHY + etalq)  * (WT / 70)^e_wt_q
    vp <- exp(lvp + e_dis_healthy_vp * DIS_HEALTHY)          * (WT / 70)^e_wt_vp

    # ------------------------------------------------------------------
    # Individual PD parameters. PD layer was estimated from AHF patients
    # only; no DIS_HEALTHY effect on PD parameters in the source. Imax is
    # the magnitude of the inhibitory effect (positive); the negative
    # direction is encoded explicitly in the drug-effect equation below.
    # ------------------------------------------------------------------
    s0    <- exp(lrbase + etalrbase)
    sss   <- exp(lsss   + etalsss)
    ec50  <- exp(lec50  + etalec50)
    tprog <- exp(ltprog)
    imax  <- exp(limax)

    # ------------------------------------------------------------------
    # PK micro-rate constants and two-compartment IV ODE system. Input
    # arrives at the central compartment (cmt = 'central' on dose rows;
    # IV bolus or constant-rate IV infusion).
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma concentration of nicorandil (ng/mL). Dose units mg, vc units
    # L -> central / vc is mg/L; convert mg/L to ng/mL by multiplying by
    # 1000.
    Cc <- 1000 * central / vc

    # ------------------------------------------------------------------
    # Asymptotic-exponential disease progression on PAWP (paper Equation
    # 4): PAWP_base(t) = S0 + (Sss - S0) * (1 - exp(-log(2) / Tprog * t)),
    # i.e. an exponential approach from S0 at time 0 to Sss as t goes to
    # infinity, with half-life Tprog. (Sss - S0) is negative for typical
    # AHF parameters (19.5 - 25.6 = -6.1 mmHg), so the baseline decreases
    # over time.
    # ------------------------------------------------------------------
    progress <- (sss - s0) * (1 - exp(-log(2) / tprog * t))

    # Inhibitory Emax drug effect on PAWP (paper Equation 2 with E_max
    # negative). Imax is stored as a positive magnitude; the negative
    # sign is explicit in the drug-effect expression.
    pdeff <- -imax * Cc / (ec50 + Cc)

    # PAWP observation (mmHg) = baseline + disease progression + drug.
    pawp <- s0 + progress + pdeff

    # Multi-output residual error: Cc combined proportional + additive
    # (paper Equation 8); pawp additive only (paper Methods, 'Random
    # effects', 'Residual unidentified variability').
    Cc   ~ add(addSd)      + prop(propSd)
    pawp ~ add(addSd_pawp)
  })
}
