Kuchimanchi_2018_evolocumab_ldlc <- function() {
  description <- "Joint population PK + static Emax-on-AUC exposure-response model for evolocumab LDL-C lowering in adults with hypercholesterolemia (Kuchimanchi 2018). The PK layer (Table 3) is the one-compartment model with parallel linear and Michaelis-Menten elimination and SC bioavailability from the companion Kuchimanchi_2018_evolocumab.R file. The PD layer (Table 4) is an algebraic Emax model linking AUC over weeks 8-12 of dosing to the mean week-10-and-12 LDL-C reduction, with statin / ezetimibe / HeFH covariate effects on baseline LDL-C, a statin covariate effect on Emax, and a regimen-effect multiplier on EC50 distinguishing once-monthly (QM) from once-every-2-weeks (Q2W) dosing. AUC of evolocumab is integrated inside an extra rxode2 state over the 56-84-day window; the LDLC observable is meaningful only at t >= 84 (vignette documents the time-window discipline)."
  reference <- "Kuchimanchi M, Monine M, Kandadi Muralidharan K, Woodhead JL, Horner TJ. Population pharmacokinetics and exposure-response modeling and simulation for evolocumab in healthy volunteers and patients with hypercholesterolemia. J Pharmacokinet Pharmacodyn. 2019;46(2):133-148. doi:10.1007/s10928-018-9592-y"
  vignette <- "Kuchimanchi_2018_evolocumab"
  paper_specific_compartments <- c("auc_wk8_12")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL", ldlc = "mg/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate (WT/84)^exponent on CL (0.276), V (1.04), and Vmax (0.145) in the PK layer (Table 3). Reference 84 kg = mean body weight of the pooled phase 1-3 analysis population. Not retained as a covariate at the PD / exposure-response layer.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Female-sex exponent on V (1.11) in the PK layer (Table 3). Reference patient is male.",
      source_name        = "SEXF"
    ),
    CONMED_STATIN_MONO = list(
      description        = "Concomitant statin monotherapy indicator, 1 = patient on a statin only (no other lipid-lowering comedication), 0 = otherwise",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on statin monotherapy)",
      notes              = "Multiplicative exponent 1.13 on Vmax (PK layer, Table 3); multiplicative exponent 0.797 on baseline LDL-C and multiplicative exponent 0.937 on Emax in the exposure-response layer (Table 4). Reference patient is not on any lipid-lowering medication.",
      source_name        = "CONMED_STATIN_MONO"
    ),
    CONMED_EZE = list(
      description        = "Concomitant ezetimibe indicator, 1 = patient taking ezetimibe (with or without other lipid-lowering comedication), 0 = not on ezetimibe",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on ezetimibe)",
      notes              = "Kuchimanchi 2018 Methods defines this as 'all patients on ezetimibe, regardless of comedications'; in the popPK dataset ~79% of ezetimibe users were also on a statin, so the effect represents statin + ezetimibe combination therapy. Multiplicative exponent 1.20 on Vmax (PK layer, Table 3); multiplicative exponent 0.768 on baseline LDL-C in the exposure-response layer (Table 4).",
      source_name        = "CONMED_EZE"
    ),
    PCSK9 = list(
      description        = "Baseline unbound PCSK9 (proprotein convertase subtilisin/kexin type 9) serum concentration",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate (PCSK9/425)^0.194 on Vmax in the PK layer (Table 3). Reference 425 ng/mL (= 5.9 nM) is the population median. Baseline (time-fixed) covariate; patients with missing baseline PCSK9 were excluded from analyses that included PCSK9 as a covariate. Not retained at the PD / exposure-response layer.",
      source_name        = "PCSK9"
    ),
    DIS_HEFH = list(
      description        = "Heterozygous familial hypercholesterolemia indicator, 1 = patient with HeFH, 0 = non-HeFH",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-HeFH)",
      notes              = "Multiplicative exponent 1.28 on baseline LDL-C in the exposure-response layer (Table 4). HeFH patients have higher baseline LDL-C than the non-HeFH reference. Kuchimanchi 2018 Results notes that all HeFH patients in the dataset were taking a statin and/or ezetimibe, so any HeFH simulation should typically be combined with CONMED_STATIN_MONO = 1 or CONMED_EZE = 1.",
      source_name        = "HEFH"
    ),
    REGI_QM = list(
      description        = "Once-monthly dosing-regimen indicator, 1 = QM (e.g., 420 mg SC every 4 weeks), 0 = Q2W (e.g., 140 mg SC every 2 weeks)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Q2W; the structural reference for the Emax-on-AUC EC50 estimate)",
      notes              = "Per-subject (regimen-fixed) categorical indicator used by the static Emax-on-AUC exposure-response model to apply the QM-specific EC50 multiplier reg_qm (Table 4 REG = 2.30): ec50_eff = ec50 * reg_qm^REGI_QM. Required because the Emax-on-AUC formulation cannot represent the difference in target-saturation time courses between Q2W and QM regimens that produce similar AUC values. Set REGI_QM = 1 for QM simulations and REGI_QM = 0 for Q2W simulations; both values are valid inputs to the model regardless of the actual dosing interval used in the event table, but the canonical interpretation is to keep REGI_QM consistent with the simulated regimen.",
      source_name        = "REGI_QM"
    )
  )

  population <- list(
    species          = "human",
    n_subjects_pk    = 3414L,
    n_subjects_pkpd  = 1314L,
    n_observations   = 16179L,
    n_studies        = 11L,
    n_studies_pkpd   = 4L,
    age_range        = "18-80 years",
    age_median       = "57 years (mean; Table 2 reports SD 58 which appears to be a typographical error)",
    weight_range     = "41-175 kg",
    weight_median    = "84.2 kg (mean)",
    sex_female_pct   = 50,
    race_ethnicity   = c(White = 87, Black = 7, Asian = 4, Hispanic = 0, Other = 1, AmericanIndianAlaska = 0, NativeHawaiianPacific = 0, Multiple = 0),
    disease_state    = "Pooled adults: healthy volunteers (phase 1a) and patients with hypercholesterolemia (phase 1b, 2, and 3), including patients with heterozygous familial hypercholesterolemia (9%), diabetes (11%), and statin-intolerance cohorts. Most subjects received concomitant lipid-lowering therapy (statins 72%, ezetimibe 12%).",
    dose_range       = "Evolocumab 7-420 mg IV or SC, single- and multiple-dose across Q2W and QM regimens. Phase 3 studies used the commercial regimens 140 mg SC Q2W and 420 mg SC QM.",
    regions          = "Multi-regional (11 pooled clinical studies spanning phase 1, 2, and 3).",
    pcsk9_baseline   = "Mean 402 ng/mL (SD 375), range 15.5-1233 ng/mL; median used for reference patient = 425 ng/mL (= 5.9 nM).",
    ldlc_baseline    = "Typical baseline LDL-C used by the exposure-response model is 150 mg/dL for the reference patient (Table 4 typical value); the paper-predicted maximal LDL-C reduction is 99.7 mg/dL.",
    notes            = "PK parameters from Kuchimanchi 2018 Table 3 (updated phase 3 popPK model; N = 3414 evolocumab-treated patients pooled from 11 studies). Exposure-response parameters from Table 4 (n = 1314 patients pooled from 4 phase 2 studies). The PK layer fits the full pooled phase 1-3 dataset; the exposure-response layer is fit only on the phase 2 LDL-C measurements at the mean of weeks 10 and 12 of treatment, using individual-predicted AUCwk8-12 from the popPK model as the predictor."
  )

  ini({
    # ===========================================================================
    # PK structural parameters (Table 3, updated phase 3 popPK model)
    # ===========================================================================
    # Reference body weight = 84 kg; reference patient is male, not on lipid-
    # lowering medication, with baseline PCSK9 = 425 ng/mL (paper Methods).
    lka     <- fixed(log(0.319)); label("Absorption rate constant ka (1/day)")                       # Table 3: ka = 0.319 day^-1 (FIXED from phase 1a dense-data modeling)
    lcl     <- log(0.105);        label("Linear clearance CL for an 84 kg reference patient (L/day)") # Table 3: CL = 0.105 L/day (updated phase 3 popPK model)
    lvc     <- log(5.18);         label("Central volume of distribution V for an 84 kg male reference patient (L)") # Table 3: V = 5.18 L
    lfdepot <- fixed(log(0.72));  label("Subcutaneous bioavailability F (fraction; IV reference)")    # Table 3: F = 0.72 (FIXED; SC absolute bioavailability relative to IV)

    # ---- Michaelis-Menten (target-mediated) elimination from central ----
    # Paper parameterization: dC/dt-like rate = Vmax * C / (km + C), with both
    # Vmax and km reported in target-concentration units (nM). In the rxode2
    # mass-based ODE we multiply by V to get mass/day (see model() block).
    lvmax   <- fixed(log(9.85));  label("Nonlinear clearance capacity Vmax (nM/day; concentration-per-time form)") # Table 3: Vmax = 9.85 nM/day (FIXED in updated phase 3 popPK model)
    lkm     <- fixed(log(27.3));  label("Michaelis-Menten constant km (nM)")                          # Table 3: km = 27.3 nM (FIXED)

    # ---- PK covariate exponents (Table 3) ----
    e_wt_cl              <- 0.276;  label("Power exponent of WT/84 on CL (unitless)")                          # Table 3: Body weight exponent on CL = 0.276 (RSE 30.4%)
    e_wt_vc              <- 1.04;   label("Power exponent of WT/84 on Vc (unitless)")                          # Table 3: Body weight exponent on V = 1.04 (RSE 4.05%)
    e_sexf_vc            <- 1.11;   label("Female-vs-male multiplicative exponent on Vc (unitless)")           # Table 3: Female exponent on V = 1.11 (RSE 1.42%)
    e_wt_vmax            <- 0.145;  label("Power exponent of WT/84 on Vmax (unitless)")                        # Table 3: Body weight exponent on Vmax = 0.145 (RSE 33.0%)
    e_smono_vmax         <- 1.13;   label("Statin-monotherapy multiplicative exponent on Vmax (unitless)")     # Table 3: Statin exponent on Vmax = 1.13 (RSE 1.02%)
    e_conmed_eze_vmax    <- 1.20;   label("Ezetimibe (combination-therapy) multiplicative exponent on Vmax (unitless)") # Table 3: Statin + ezetimibe exponent on Vmax = 1.20 (RSE 1.59%)
    e_pcsk9_vmax         <- 0.194;  label("Power exponent of PCSK9/425 (ng/mL) on Vmax (unitless)")            # Table 3: PCSK9 baseline exponent on Vmax = 0.194 (RSE 7.47%)

    # ---- PK IIV (Table 3, omega^2 = log(CV^2 + 1)) ----
    #   CV 54.3% -> log(0.543^2 + 1) = 0.25839 (CL)
    #   CV 28.3% -> log(0.283^2 + 1) = 0.07704 (V)
    #   CV 31.1% -> log(0.311^2 + 1) = 0.09232 (Vmax)
    #   CV 74.6% -> log(0.746^2 + 1) = 0.44245 (ka; FIXED in the paper)
    # The paper states that CL, V, and Vmax share a full-block variance matrix
    # but only reports the diagonal CV values (the full omega correlation
    # structure is not published); the diagonal-only encoding matches what is
    # reproducible from the paper and is documented as a deviation in the
    # vignette (same approximation as the companion popPK-only model file).
    etalcl   ~ 0.25839                         # Table 3 IIV CV 54.3% on CL
    etalvc   ~ 0.07704                         # Table 3 IIV CV 28.3% on V
    etalvmax ~ 0.09232                         # Table 3 IIV CV 31.1% on Vmax
    etalka   ~ fix(0.44245)                    # Table 3 IIV CV 74.6% FIXED on ka

    # ---- PK residual error (Table 3) ----
    # Paper reports proportional error as a fraction (0.282 = 28.2% CV) and
    # additive error in nM (5.41 nM). The additive residual is converted to
    # ug/mL inside model() using the evolocumab molecular weight of 141,800
    # g/mol (MW_EVO) so the residual applies directly to Cc (ug/mL).
    propSd    <- 0.282;  label("Proportional residual error on Cc (fraction)")                       # Table 3: Residual proportional error = 0.282 (RSE 1.12%)
    addSd_nM  <- 5.41;   label("Additive residual error on Cc (nM; converted to ug/mL inside model)") # Table 3: Residual additive error = 5.41 nM (RSE 2.50%)

    # ===========================================================================
    # PD / exposure-response parameters (Table 4, Emax-on-AUC LDL-C model)
    # ===========================================================================
    # Reference patient (Methods, exposure-response section): 84 kg male,
    # baseline PCSK9 = 425 ng/mL, NOT on lipid-lowering medications, non-HeFH.
    # Reported parameter values are the population typical values from
    # Table 4 (n = 1314 patients pooled across 4 phase 2 studies).
    lrbase   <- log(150);   label("Typical baseline LDL-C for the reference patient (mg/dL)")        # Table 4: Baseline LDL-C = 150 mg/dL (RSE 0.92%)
    lemax    <- log(99.7);  label("Magnitude of the maximal LDL-C reduction (mg/dL; sign applied in model)") # Table 4: Emax = -99.7 mg/dL (RSE 2.17%); paper sign is negative (reduction from baseline)
    lec50    <- log(51.5);  label("AUC wk8-12 producing half-maximal LDL-C reduction at the Q2W reference (ug/mL*day)") # Table 4: EC50 = 51.5 (ug/mL)*day at Q2W (RSE 9.79%)
    lreg_qm  <- log(2.30);  label("Log of EC50 multiplier for the QM regimen (REGI_QM) relative to Q2W (unitless)") # Table 4: REG (QM/Q2W EC50 multiplier) = 2.30 (RSE 10.3%); paper formula EC50 * REG^i with i = REGI_QM

    # ---- Exposure-response covariate effects (Table 4) ----
    # Encoded as linear-scale multipliers raised to the binary indicator: paper
    # form `theta^cov` with cov in {0,1} so that the indicator selects the
    # multiplier when it is 1 and leaves the reference value unchanged when 0.
    e_conmed_statin_mono_rbase <- 0.797;  label("Statin-monotherapy multiplicative exponent on baseline LDL-C (unitless)") # Table 4: Statin exponent on baseline LDL-C = 0.797
    e_conmed_eze_rbase         <- 0.768;  label("Ezetimibe (combination-therapy) multiplicative exponent on baseline LDL-C (unitless)") # Table 4: Ezetimibe exponent on baseline LDL-C = 0.768
    e_dis_hefh_rbase           <- 1.28;   label("HeFH multiplicative exponent on baseline LDL-C (unitless)") # Table 4: HeFH exponent on baseline LDL-C = 1.28
    e_conmed_statin_mono_emax  <- 0.937;  label("Statin-monotherapy multiplicative exponent on |Emax| (unitless)") # Table 4: Statin exponent on Emax = 0.937

    # ---- LDL-C IIV (Table 4) ----
    # 20.0% CV on baseline LDL-C -> log(0.20^2 + 1) = 0.03922. No IIV on
    # Emax / EC50 / REG (the paper reports only baseline-LDL IIV).
    etalrbase ~ 0.03922                       # Table 4 IIV CV 20.0% on baseline LDL-C

    # ---- LDL-C residual error (Table 4) ----
    # Additive 19.3 mg/dL; proportional FIXED at 0 (Table 4).
    addSd_LDLC <- 19.3;  label("Additive residual error on LDL-C (mg/dL)")                            # Table 4: Residual additive error on LDL-C = 19.3 mg/dL (RSE 1.35%)
  })

  model({
    # ===========================================================================
    # Physical constants (shared by PK + PD layers)
    # ===========================================================================
    # Evolocumab is a fully human IgG2 monoclonal antibody; molecular weight
    # 141,800 g/mol per the FDA-approved Repatha prescribing information.
    # Used to convert nM (paper units for km, Vmax, additive PK residual) to
    # ug/mL (the rxode2 concentration scale with dose in mg and V in L).
    MW_EVO  <- 141800            # Evolocumab molecular weight, g/mol
    nM_per_ugmL <- 1e6 / MW_EVO  # 1 ug/mL = 1e6 ng/L / MW (g/mol) nM; = 7.052 nM/(ug/mL)

    # ===========================================================================
    # Individual PK parameters (Table 3 reference: WT = 84 kg, male, no
    # lipid-lowering medication, baseline PCSK9 = 425 ng/mL)
    # ===========================================================================
    ka    <- exp(lka + etalka)
    cl    <- exp(lcl + etalcl)   * (WT / 84)^e_wt_cl
    vc    <- exp(lvc + etalvc)   * (WT / 84)^e_wt_vc * e_sexf_vc^SEXF
    vmax  <- exp(lvmax + etalvmax) * (WT / 84)^e_wt_vmax *
             e_smono_vmax^CONMED_STATIN_MONO * e_conmed_eze_vmax^CONMED_EZE *
             (PCSK9 / 425)^e_pcsk9_vmax   # Vmax in nM/day
    km    <- exp(lkm)            # km in nM
    fdepot <- exp(lfdepot)

    # ---- Unit bridging between nM (paper) and ug/mL (model state) ----
    # central amount is in mg; vc in L; so Cc = central/vc is in mg/L = ug/mL.
    # The MM arithmetic uses the paper's nM scale: Cc_nM = Cc (ug/mL) * 1e6 / MW.
    Cc      <- central / vc                                 # ug/mL
    Cc_nM   <- Cc * nM_per_ugmL                              # nM
    mm_rate_nM <- vmax * Cc_nM / (km + Cc_nM)                # nM/day (concentration-per-time)
    # Convert the nonlinear elimination rate back to mass/time for the ODE:
    # V * (nM/day) * (MW / 1e6 (ug/nmol)) = V [L] * (nM/day) * (ug/mL / nM) * (1 mL/0.001 L)  ->  mg/day.
    mm_rate_mgday <- vc * mm_rate_nM / nM_per_ugmL           # ug/mL/day * L = mg/L/day * L = mg/day

    # ===========================================================================
    # ODE system (Figure 1a of Kuchimanchi 2018; PK layer)
    # ===========================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - cl * Cc - mm_rate_mgday

    # ---- Bioavailability on SC depot dosing (IV doses go directly to central) ----
    f(depot) <- fdepot

    # ===========================================================================
    # AUC over weeks 8-12 (days 56-84) -- exposure-response predictor
    # ===========================================================================
    # The exposure-response model uses the AUC of unbound evolocumab serum
    # concentration over the 4-week window from day 56 (start of week 8) to
    # day 84 (end of week 12) as the predictor of mean week-10-and-12 LDL-C.
    # This extra state integrates Cc only inside the window; auc_wk8_12(t) for
    # t >= 84 is the AUCwk8-12 used in the Emax expression below.
    d/dt(auc_wk8_12) <- Cc * (t >= 56) * (t <= 84)

    # ===========================================================================
    # Individual exposure-response parameters (Table 4 reference: no
    # lipid-lowering medication, non-HeFH)
    # ===========================================================================
    rbase    <- exp(lrbase + etalrbase)
    emax_mag <- exp(lemax)
    ec50     <- exp(lec50)
    reg_qm   <- exp(lreg_qm)

    # Covariate effects on baseline LDL-C and on Emax (Table 4). Each
    # multiplicative exponent acts only when its associated indicator is 1.
    rbase_cov <- rbase *
                 e_conmed_statin_mono_rbase^CONMED_STATIN_MONO *
                 e_conmed_eze_rbase^CONMED_EZE *
                 e_dis_hefh_rbase^DIS_HEFH
    emax_cov  <- -emax_mag * e_conmed_statin_mono_emax^CONMED_STATIN_MONO
    ec50_cov  <- ec50 * reg_qm^REGI_QM

    # Emax-on-AUC exposure-response formula (Kuchimanchi 2018 equations p. 509):
    #   Eff  = Emax * AUC / (EC50 * REG^i + AUC)
    #   LDLC = baseline_LDLC + Eff   (additive model selected in Results)
    # Note: LDLC is the modelled mean-of-weeks-10-and-12 LDL-C concentration;
    # because auc_wk8_12 only finishes accumulating at t = 84, the LDLC
    # observable is only meaningful at t >= 84. Observation events in user
    # data should target an observation time at or after day 84.
    Eff   <- emax_cov * auc_wk8_12 / (ec50_cov + auc_wk8_12)
    LDLC  <- rbase_cov + Eff                                  # mg/dL

    # ===========================================================================
    # Multi-output observation and residual error
    # ===========================================================================
    # Two algebraic observables: Cc (ug/mL, PK) and LDLC (mg/dL, PD). The
    # rxode2 DVID auto-mapping assigns dvid = 1 to the first residual line
    # (Cc) and dvid = 2 to the second (LDLC) in the order they appear below.
    # Observation rows in user event tables should carry dvid = 1L for PK
    # observations and dvid = 2L for LDLC observations (see the companion
    # vignette for the canonical event-table layout).

    # PK residual error: additive (converted from nM to ug/mL) + proportional.
    addSd <- addSd_nM / nM_per_ugmL                          # ug/mL
    Cc   ~ add(addSd) + prop(propSd)
    # LDL-C residual error: additive only (proportional FIXED at 0 in Table 4).
    LDLC ~ add(addSd_LDLC)
  })
}
