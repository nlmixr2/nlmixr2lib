Fu_2022_atenolol_qsp <- function() {
  description <- paste(
    "QSP (preclinical, beagle dog). Novel cardiovascular systems",
    "(CVS-CTR) model integrating heart rate (HR), left-ventricular end-",
    "diastolic volume (EDV), total peripheral resistance (TPR), and",
    "myocardial contractility (CTR) via pressure-volume-loop theory (Fu",
    "2022). Four indirect-response (turnover) ODEs -- one per",
    "hemodynamic state -- are coupled by mean arterial pressure (MAP)",
    "feedback and by algebraic PV-loop relations that derive stroke",
    "volume (SV), end-systolic volume (ESV), cardiac output (CO), and",
    "MAP from the four states plus a fixed literature-derived isovolumic",
    "contraction time (TIC = 0.256 s, Templeton 1979). The paper uses",
    "atenolol (0.3-30 mg/kg PO) as proof-of-concept beta1-blocker;",
    "Emax inhibition on HR-production and CTR-production is driven by",
    "the central-compartment atenolol concentration with EC50 fixed at",
    "the beta1 KD (58.3 ng/mL, Baker 2005). A positive TPR effect (beta2",
    "mediated) was rejected during model development (Emax_TPR -> 0) and",
    "is not included. Atenolol PK is fixed from a prior three-compartment",
    "beagle-dog popPK model (Snelder 2013 / Venkatasubramanian 2018 as",
    "cited in Fu 2022 references 16 and 18). Circadian rhythm modulates",
    "the production rates of HR / TPR / CTR (24 h period on HR and CTR;",
    "8 h effective period on TPR); amplitudes are shared within a study",
    "but differ between the two contributing studies (Servier -- Study 1,",
    "AstraZeneca -- Study 2). Four outputs are observed with log-normal",
    "residual error (HR, dP/dtmax); CO and MAP are derived algebraically.",
    "EDV is included structurally (Fu 2022 Discussion: 'baseline EDV was",
    "fixed to a reported literature value' 31.13 mL) to enable future",
    "extension to compounds with a primary effect on EDV.",
    sep = " "
  )
  reference <- paste(
    "Fu Y, Taghvafard H, Said MM, Rossman EI, Collins TA, Billiald-",
    "Desquand S, Leishman D, van der Graaf PH, van Hasselt JGC, Snelder",
    "N (2022). A novel cardiovascular systems model to quantify drugs",
    "effects on the inter-relationship between contractility and other",
    "hemodynamic variables. CPT Pharmacometrics Syst Pharmacol.",
    "11(5):640-652. doi:10.1002/psp4.12774.",
    "Atenolol PK structure fixed from Fu 2022 Methods (Pharmacokinetic",
    "model) as cited to Venkatasubramanian 2018 and Snelder 2013.",
    "TIC = 0.256 s from Templeton 1979 (Fu 2022 reference 21).",
    sep = " "
  )
  vignette <- "Fu_2022_atenolol_qsp"

  paper_specific_compartments <- c("hr", "edv", "tpr", "ctr")

  units <- list(
    time          = "h",
    dosing        = "mg/kg (per-kg atenolol oral dose; the PK typicals below are per-kg-normalized so amt in mg/kg with vc in L/kg yields amt/vc in mg/L; the model multiplies by 1000 to get Cc in ng/mL for comparison with the paper's EC50)",
    concentration = "ng/mL (central-compartment atenolol; the paper's Emax and EC50 values are on the ng/mL scale)"
  )

  covariateData <- list(
    STUDY_FU2022_AZ = list(
      description        = "Fu 2022 pooled-analysis study indicator: 1 = subject enrolled in Study 2 (AstraZeneca; Alderley Park, UK; 4 male beagle dogs, 14.2-14.6 kg, 17-22 months old; oral atenolol 0, 1, 3, 10 mg/kg; HR, dP/dtmax, and MAP measured; NO cardiac output measurement); 0 = Study 1 (Servier; France; 4 male beagle dogs, 10-15 kg; oral atenolol 0, 3, 10, 30 mg/kg; HR, dP/dtmax, CO, and MAP measured). Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Selects the study-specific baseline values (BSL_HR, V0, BSL_CTRM) and the study-specific circadian amplitude and horizontal displacements for HR and CTR. TPR-related circadian parameters (Amp_TPR = per-study Amp; Hor_TPR = 6.33 h) and all remaining system-specific parameters (BSL_TPR, BSL_EDV, Kout, FB) are shared across the two studies (Fu 2022 Table 2 footnotes). Study 3 (GlaxoSmithKline, 0.3/1/3 mg/kg) is external validation only and is not encoded as a covariate level; for external-validation simulations, use STUDY_FU2022_AZ = 0 with the Study 1 typical values.",
      source_name        = "SSID (in the NONMEM control stream; SSID = 1 -> Servier, SSID = 2 -> AstraZeneca, SSID = 3 -> GSK external validation)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight (documented in the supplement for each study cohort but NOT used as a covariate in the final Fu 2022 model). Screened only as a description of the animal cohort; the beagle-dog weight range is narrow (9-15 kg across the three studies) and no allometric scaling is applied.",
      units              = "kg",
      type               = "continuous",
      notes              = "Documented-but-unused: recorded in Fu 2022 supplement Section A per-study animal descriptions; not referenced in model()."
    )
  )

  population <- list(
    species        = "beagle dog",
    n_subjects     = 12L,
    n_studies      = 3L,
    age_range      = "17-72 months (Study 1: healthy naive adult; Study 2: 17-22 months; Study 3: 5-6 years)",
    weight_range   = "9-15 kg (Study 1: 10-15 kg; Study 2: 14.2-14.6 kg; Study 3: 9-13 kg)",
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state  = "Conscious chronically-instrumented healthy beagle dogs (male). Multi-site consortium: Servier (France; Study 1), AstraZeneca (Alderley Park, UK; Study 2), GlaxoSmithKline (Marshall Farms, NY, USA; Study 3, external validation). Hemodynamic markers monitored by telemetry: aortic and left-atrial pressures, LV pressure via solid-state micromanometer, aortic blood flow via transit-time flowmeter (Study 1) or DSI PhysioTel implants (Studies 2 and 3).",
    dose_range     = "Oral gavage; increasing atenolol doses with washout in between: 0/3/10/30 mg/kg (Study 1, vehicle = 0.5 percent methylcellulose); 0/1/3/10 mg/kg (Study 2, vehicle = water); 0/0.3/1/3 mg/kg (Study 3, vehicle = water). Studies 1 and 2 used for model development; Study 3 for external validation.",
    regions        = "France (Servier), United Kingdom (AstraZeneca Alderley Park), USA (GlaxoSmithKline).",
    notes          = "Data from three in vivo telemetry studies (Fu 2022 Table 1 and supplement Section A). Only HR, LV dP/dtmax, CO (Studies 1 and 3 only), and MAP time courses were used. Data from Studies 1 and 2 were simultaneously fit with NONMEM 7.4.3 (FOCE-INTER) via PsN 4.8.1. Model was initialised at 0 h and dosing began at 168 h in the fitting run so that the circadian rhythms were in oscillating steady state; for typical-user simulations the packaged model starts the CVS states at their baseline values with the circadian phase set so that t = 0 corresponds to the model's dosing time. Study 1 data used for external comparison of the developed model to CO-informative telemetry; the paper additionally fit a Model without CO data (with BSL_TPR fixed to 0.0743 mmHg*min/mL, IIV on BSL_TPR and CS_TPR removed, and Emax_TPR retained at zero) -- the packaged model corresponds to the primary Fu 2022 Table 2 final model with CO data."
  )

  ini({
    # =====================================================================
    # ATENOLOL PK -- three-compartment first-order absorption, all typical
    # values FIXED at the Fu 2022 Methods (Pharmacokinetic model) values
    # (cited to Venkatasubramanian 2018 / Snelder 2013). Volumes and
    # clearances are per-kg (L/kg, L/kg/h), so central-compartment amount
    # divided by vc yields Cc in mg/(L/kg) = ng/mL when amt is entered as
    # ug/kg or mg/kg * 1000; see units$dosing above.
    # =====================================================================
    lka      <- fixed(log(1.13))        ; label("Atenolol absorption rate ka (1/h; per Fu 2022 Methods)")             # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lcl      <- fixed(log(3.35))        ; label("Atenolol systemic clearance CL (L/kg/h; per Fu 2022 Methods)")       # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lvc      <- fixed(log(4.05))        ; label("Atenolol central volume Vc = V2 (L/kg; per Fu 2022 Methods)")        # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lq       <- fixed(log(8.85))        ; label("Atenolol inter-compartmental clearance Q = Q2 (L/kg/h)")            # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lvp      <- fixed(log(3.74))        ; label("Atenolol first-peripheral volume Vp = V3 (L/kg)")                   # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lq2      <- fixed(log(5.73))        ; label("Atenolol inter-compartmental clearance Q2 = Q3 (L/kg/h)")           # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lvp2     <- fixed(log(11.9))        ; label("Atenolol second-peripheral volume Vp2 = V4 (L/kg)")                 # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph
    lfdepot  <- fixed(log(0.783))       ; label("Atenolol oral bioavailability F1 (unitless)")                       # Fu 2022 Methods p.641 'Pharmacokinetic model' paragraph

    # =====================================================================
    # CVS-CTR SYSTEM -- SHARED parameters (Fu 2022 Table 2)
    # =====================================================================
    lbsl_tpr <- log(0.0743)             ; label("Baseline total peripheral resistance BSL_TPR (mmHg*min/mL)")                # Fu 2022 Table 2 (final model column) 0.0743 (RSE 4.55%)
    lbsl_edv <- fixed(log(31.13))       ; label("Baseline end-diastolic volume BSL_EDV (mL; literature value)")     # Fu 2022 Table 2 (fixed; the paper attributes the value to reference 29 and notes it must be fixed to avoid over-parameterization)
    lkout    <- log(0.830)              ; label("Shared turnover elimination rate Kout (1/h) for HR / EDV / TPR / CTR")      # Fu 2022 Table 2 (final model column) 0.830 (RSE 21.8%); Kout_HR = Kout_EDV = Kout_TPR = Kout_CTR = Kout (Table 2 footnote a)
    lfb      <- log(0.00558)            ; label("MAP-mediated negative-feedback strength FB (1/mmHg)")                       # Fu 2022 Table 2 (final model column) 0.00558 (RSE 12.5%)
    ltic     <- fixed(log(0.256))       ; label("Isovolumic contraction time TIC (s;, from Templeton 1979)")           # Fu 2022 Methods Eq 7 (TIC = 0.256 s reported for a typical young Beagle by Templeton 1979; kept in seconds so the CTRM = CTR * EDV / TIC identity retains the source's mmHg/s units for dP/dtmax)

    # Circadian: horizontal displacement of TPR is shared across studies (24-h period; but TPR uses a 3x-frequency cosine ==> effective period = 8 h per Fu 2022 Methods 'Circadian rhythms in hemodynamic variables' -- the paper tested 8/12/24 h and identified the optimal period per variable).
    hor_tpr  <- 6.33                    ; label("Circadian horizontal displacement Hor_TPR (h)")                             # Fu 2022 Table 2 (final model column) 6.33 (RSE 1.96%)

    # =====================================================================
    # CVS-CTR SYSTEM -- STUDY-SPECIFIC baselines (Fu 2022 Table 2)
    # STUDY_FU2022_AZ = 0 -> Study 1 (Servier) typical values
    # STUDY_FU2022_AZ = 1 -> Study 2 (AstraZeneca) typical values
    # =====================================================================
    lbsl_hr_s1   <- log(79.4)           ; label("Baseline heart rate BSL_HR Study 1 (beats/min)")                            # Fu 2022 Table 2 (final model column) 79.4 (RSE 5.78%)
    lbsl_hr_s2   <- log(77.0)           ; label("Baseline heart rate BSL_HR Study 2 (beats/min)")                            # Fu 2022 Table 2 (final model column) 77.0 (RSE 2.63%)
    lv0_s1       <- log(9.92)           ; label("PV-loop x-intercept V0 Study 1 (mL)")                                       # Fu 2022 Table 2 (final model column) 9.92 (RSE 9.81%)
    lv0_s2       <- log(9.15)           ; label("PV-loop x-intercept V0 Study 2 (mL)")                                       # Fu 2022 Table 2 (final model column) 9.15 (RSE 9.62%)
    lbsl_ctrm_s1 <- log(3777)           ; label("Baseline dP/dtmax BSL_CTRM Study 1 (mmHg/s)")                               # Fu 2022 Table 2 (final model column) 3777 (RSE 7.59%)
    lbsl_ctrm_s2 <- log(2422)           ; label("Baseline dP/dtmax BSL_CTRM Study 2 (mmHg/s)")                               # Fu 2022 Table 2 (final model column) 2422 (RSE 7.50%)

    # Study-specific circadian amplitude (shared across HR / CTR / TPR
    # within a study per the CTL: AMP2 = AMP1 x THETA(27) with THETA(27) =
    # 1 FIX; AMP3 = AMP1 x THETA(30) with THETA(30) = 1 FIX; Fu 2022
    # Results paragraph 6 confirms 'The amplitude parameters amp_HR, amp_CTR,
    # and amp_TPR could not be distinguished and were estimated to be
    # 0.0931 for study 1 and 0.168 for study 2').
    lamp_s1     <- log(0.0931)          ; label("Circadian amplitude Amp Study 1 (fraction; shared HR/CTR/TPR within study)") # Fu 2022 Table 2 (final model column) 0.0931 (RSE 41.3%)
    lamp_s2     <- log(0.168)           ; label("Circadian amplitude Amp Study 2 (fraction; shared HR/CTR/TPR within study)") # Fu 2022 Table 2 (final model column) 0.168 (RSE 13.6%)

    # Study-specific circadian horizontal displacements
    hor_hr_s1   <- 7.86                 ; label("Circadian horizontal displacement Hor_HR Study 1 (h)")                       # Fu 2022 Table 2 (final model column) 7.86 (RSE 28.7%)
    hor_hr_s2   <- 19.4                 ; label("Circadian horizontal displacement Hor_HR Study 2 (h)")                       # Fu 2022 Table 2 (final model column) 19.4 (RSE 1.60%)
    hor_ctr_s1  <- 9.82                 ; label("Circadian horizontal displacement Hor_CTR Study 1 (h)")                      # Fu 2022 Table 2 (final model column) 9.82 (RSE 18.9%)
    hor_ctr_s2  <- 21.8                 ; label("Circadian horizontal displacement Hor_CTR Study 2 (h)")                      # Fu 2022 Table 2 (final model column) 21.8 (RSE 1.84%)

    # =====================================================================
    # DRUG-SPECIFIC PARAMETERS (Fu 2022 Table 2). Emax on HR and CTR
    # inhibits their zero-order production rates (beta1-mediated). Emax
    # on TPR was fixed to 0 (rejected -- 'the final estimate of Emax for
    # TPR approached zero, and the OFV did not decrease significantly');
    # not included in the final model. EC50s are FIXED at the beta1 KD
    # (58.3 ng/mL) as reported by Baker 2005 (Fu 2022 reference 23).
    # =====================================================================
    lemax_hr <- log(0.415)              ; label("Emax inhibition of HR production by atenolol (fraction)")                   # Fu 2022 Table 2 (final model column) 0.415 (RSE 11.6%)
    lec50_hr <- fixed(log(58.3))        ; label("EC50 on HR (ng/mL; beta1 KD)")                                     # Fu 2022 Table 2 (fixed at Baker 2005 KD_beta1 = 58.3 ng/mL)
    lemax_ctr <- log(0.422)             ; label("Emax inhibition of CTR production by atenolol (fraction)")                  # Fu 2022 Table 2 (final model column) 0.422 (RSE 9.56%)
    lec50_ctr <- fixed(log(58.3))       ; label("EC50 on CTR (ng/mL; beta1 KD; = EC50_HR in the final model)")     # Fu 2022 Table 2 (fixed; CTL EC50_CTR = EC50_HR x THETA(11) with THETA(11) = 1 FIX)

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY (Fu 2022 Table 2)
    # CV% reported; internal variance = log(1 + (CV/100)^2) for a log-normal
    # random effect on the underlying baseline.
    # =====================================================================
    etalbsl_hr ~ 0.003690                                                                                                     # Fu 2022 Table 2 BSL_HR CV% = 6.08 -> var = log(1+0.0608^2) = 0.003690
    # IIV on V0 was fixed to 0 in NONMEM $OMEGA (OM2, `0 FIX`); not encoded here.
    etalbsl_tpr + etalbsl_ctrm ~ c(0.004190, 0.007854, 0.027119)                                                             # Fu 2022 Table 2 BSL_TPR CV% = 6.48, BSL_CTRM CV% = 16.58, off-diagonal CV% = 8.88 (interpretation: sqrt(cov) x 100; correlation approximately 0.74)

    # =====================================================================
    # RESIDUAL ERROR (Fu 2022 Table 2)
    # 'Exponential' NONMEM residual (Y = IPRED * exp(EPS)) maps to nlmixr2
    # `~ lnorm(expSd)`. TPR residual was reported as 'very small' and fixed
    # to zero in the CTL ($SIGMA `0 FIX`); no residual is encoded on TPR.
    # =====================================================================
    expSd_HR   <- 0.115                 ; label("Log-normal residual SD on HR (log-scale; CV% approximately 11.53)")           # Fu 2022 Table 2 Res.Error_HR CV% = 11.53 -> expSd = sqrt(log(1+0.1153^2)) = 0.115
    expSd_CTRM <- 0.104                 ; label("Log-normal residual SD on dP/dtmax (log-scale; CV% approximately 10.43)")      # Fu 2022 Table 2 Res.Error_dPdtmax CV% = 10.43 -> expSd = sqrt(log(1+0.1043^2)) = 0.104
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual PK parameters (typicals FIXED; no PK IIV encoded --
    #    PK was not measured in the CVS-CTR studies and is fixed from the
    #    Fu 2022 Methods 'Pharmacokinetic model' upstream reference).
    # ---------------------------------------------------------------------
    ka  <- exp(lka)
    cl  <- exp(lcl)
    vc  <- exp(lvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    q2  <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---------------------------------------------------------------------
    # 2. Study-selected typical baseline values (log-scale switch, then
    #    add IIV, then exp) -- STUDY_FU2022_AZ is a binary covariate:
    #    0 -> Study 1 (Servier); 1 -> Study 2 (AstraZeneca).
    # ---------------------------------------------------------------------
    lbsl_hr_typ   <- lbsl_hr_s1   * (1 - STUDY_FU2022_AZ) + lbsl_hr_s2   * STUDY_FU2022_AZ
    lv0_typ       <- lv0_s1       * (1 - STUDY_FU2022_AZ) + lv0_s2       * STUDY_FU2022_AZ
    lbsl_ctrm_typ <- lbsl_ctrm_s1 * (1 - STUDY_FU2022_AZ) + lbsl_ctrm_s2 * STUDY_FU2022_AZ
    lamp_typ      <- lamp_s1      * (1 - STUDY_FU2022_AZ) + lamp_s2      * STUDY_FU2022_AZ
    hor_hr_typ    <- hor_hr_s1    * (1 - STUDY_FU2022_AZ) + hor_hr_s2    * STUDY_FU2022_AZ
    hor_ctr_typ   <- hor_ctr_s1   * (1 - STUDY_FU2022_AZ) + hor_ctr_s2   * STUDY_FU2022_AZ

    bsl_hr   <- exp(lbsl_hr_typ + etalbsl_hr)
    v0       <- exp(lv0_typ)                       # IIV on V0 fixed at 0 (Fu 2022 Table 2 IIV row for BSL_HR only; V0 IIV set to 0 FIX in the CTL)
    bsl_tpr  <- exp(lbsl_tpr    + etalbsl_tpr)
    bsl_ctrm <- exp(lbsl_ctrm_typ + etalbsl_ctrm)
    bsl_edv  <- exp(lbsl_edv)                      # FIXED at 31.13 mL per Fu 2022 Table 2
    kout     <- exp(lkout)
    fb       <- exp(lfb)
    tic      <- exp(ltic)                          # 0.256 s -- Fu 2022 CTL keeps TIC in s; CTRM is mmHg/s and CTR is mmHg/mL so tic_s in seconds keeps the identity CTRM = CTR * EDV / TIC in mmHg/s
    amp      <- exp(lamp_typ)

    # ---------------------------------------------------------------------
    # 3. Derived baseline PV-loop quantities (Fu 2022 Methods Eqs 3-7)
    #    BSL_CTR = BSL_CTRM * TIC / BSL_EDV
    #    BSL_EA  = BSL_HR   * BSL_TPR
    #    BSL_ESV = (BSL_EDV * BSL_EA + V0 * BSL_CTR) / (BSL_EA + BSL_CTR)
    #    BSL_SV  = BSL_EDV - BSL_ESV
    #    BSL_CO  = BSL_HR  * BSL_SV
    #    BSL_MAP = BSL_TPR * BSL_CO
    # TIC is in seconds throughout (see ini()), so BSL_CTRM (mmHg/s) *
    # TIC (s) / BSL_EDV (mL) yields BSL_CTR (mmHg/mL) directly, matching
    # the source CTL. The CVS ODE time base is hours (Kout in 1/h) --
    # TIC appears only in the algebraic CTRM <-> CTR / EDV identity and
    # does not need to be in the same time unit as the ODEs.
    # ---------------------------------------------------------------------
    bsl_ctr <- bsl_ctrm * tic / bsl_edv
    bsl_ea  <- bsl_hr   * bsl_tpr
    bsl_esv <- (bsl_edv * bsl_ea + v0 * bsl_ctr) / (bsl_ea + bsl_ctr)
    bsl_sv  <- bsl_edv - bsl_esv
    bsl_co  <- bsl_hr   * bsl_sv
    bsl_map <- bsl_tpr  * bsl_co

    # ---------------------------------------------------------------------
    # 4. Steady-state production-rate constants (Fu 2022 Methods Eq 8 /
    #    'The system was forced to be in steady-state ...'). MAP feedback
    #    is baked into k_in for HR / TPR / CTR using the STEADY-STATE
    #    baseline MAP (bsl_map); EDV has no MAP feedback in the source CTL
    #    (KINEDV = KOUTEDV * BSLEDV, no `/ (1 - FB*MAP)` term).
    # ---------------------------------------------------------------------
    kin_hr  <- kout * bsl_hr  / (1 - fb * bsl_map)
    kin_tpr <- kout * bsl_tpr / (1 - fb * bsl_map)
    kin_ctr <- kout * bsl_ctr / (1 - fb * bsl_map)
    kin_edv <- kout * bsl_edv

    # ---------------------------------------------------------------------
    # 5. Emax drug-effects on HR and CTR production (Fu 2022 Methods
    #    'Modelling of the concentration-effect relationship' with Hill = 1).
    #    Cc is the atenolol central-compartment concentration in ng/mL --
    #    per-kg amt (mg/kg) / per-kg vc (L/kg) yields mg/L = ug/mL; the
    #    x1000 scales to ng/mL so the paper's EC50 = 58.3 ng/mL is used
    #    directly without a units conversion inline.
    # ---------------------------------------------------------------------
    emax_hr  <- exp(lemax_hr)
    ec50_hr  <- exp(lec50_hr)
    emax_ctr <- exp(lemax_ctr)
    ec50_ctr <- exp(lec50_ctr)

    Cc <- 1000 * central / vc

    effhr  <- emax_hr  * Cc / (ec50_hr  + Cc)
    effctr <- emax_ctr * Cc / (ec50_ctr + Cc)

    # ---------------------------------------------------------------------
    # 6. Circadian rhythms (Fu 2022 Methods 'Circadian rhythms in
    #    hemodynamic variables'; 24 h period on HR / CTR, 8 h effective
    #    period on TPR encoded as a 3x-frequency cosine of the 24 h base).
    # ---------------------------------------------------------------------
    cs_hr  <- amp * cos(2 * pi * (t + hor_hr_typ)  / 24)
    cs_ctr <- amp * cos(2 * pi * (t + hor_ctr_typ) / 24)
    cs_tpr <- amp * cos(3 * 2 * pi * (t + hor_tpr) / 24)

    # ---------------------------------------------------------------------
    # 7. Algebraic PV-loop derivations from the CURRENT CVS states
    #    (Fu 2022 Methods Eqs 1-6). MAP is computed here so it can be
    #    used in the feedback term below; CO / MAP / SV are also
    #    exposed as named outputs.
    # ---------------------------------------------------------------------
    EA   <- hr * tpr
    ESV  <- (edv * EA + v0 * ctr) / (EA + ctr)
    SV   <- edv - ESV
    CO   <- hr * SV                                # mL/min (cardiac output)
    MAP  <- tpr * CO                               # mmHg (mean arterial pressure)
    HR   <- hr                                     # named observation variable (beats/min)
    CTRM <- ctr * edv / tic                        # dP/dtmax (mmHg/s), Fu 2022 Methods Eq 7

    # ---------------------------------------------------------------------
    # 8. MAP feedback -- capped at (FB*MAP) < 0.99 to prevent division-by-
    #    zero or negative production if MAP spikes. The source CTL uses
    #    `IF (FB*MAP.GE.1) FDB = 0.99`; for realistic hemodynamic MAP the
    #    cap never triggers (fb approximately 0.006, MAP approximately 100
    #    mmHg -> fb*MAP approximately 0.6).
    # ---------------------------------------------------------------------
    fdb <- min(fb * MAP, 0.99)

    # ---------------------------------------------------------------------
    # 9. ODE system (Fu 2022 Methods Eqs 8 and 'Incorporation of drug
    #    effects in the systems model' box). Order of state declarations
    #    (depot / central / peripheral1 / peripheral2 / hr / edv / tpr /
    #    ctr) matches the NONMEM $MODEL COMP block.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    d/dt(hr)  <- kin_hr  * (1 + cs_hr)  * (1 - fdb) * (1 - effhr)  - kout * hr
    d/dt(edv) <- kin_edv                                             - kout * edv
    d/dt(tpr) <- kin_tpr * (1 + cs_tpr) * (1 - fdb)                  - kout * tpr
    d/dt(ctr) <- kin_ctr * (1 + cs_ctr) * (1 - fdb) * (1 - effctr) - kout * ctr

    # ---------------------------------------------------------------------
    # 10. Steady-state initial conditions (Fu 2022 CTL `A_0(...)` block).
    #     The circadian phase is taken relative to t = 0; users that want
    #     to replicate the paper's simulation (dosing at t = 168 h; system
    #     initialised at 0 h so the circadian rhythm has reached oscillating
    #     steady state) should set their dosing event time accordingly.
    # ---------------------------------------------------------------------
    hr(0)  <- bsl_hr
    edv(0) <- bsl_edv
    tpr(0) <- bsl_tpr
    ctr(0) <- bsl_ctr

    # ---------------------------------------------------------------------
    # 11. Bioavailability of the depot compartment (F1 = 0.783 FIXED).
    # ---------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # ---------------------------------------------------------------------
    # 12. Residual error on the two primary observations (HR and dP/dtmax).
    #     CO and MAP are simulated deterministically -- their published
    #     residual (via propagation of ERR(1)/ERR(2) through the algebraic
    #     PV-loop relations, Fu 2022 CTL $ERROR block) is not encoded in
    #     the packaged model. TPR residual was fixed to zero in the source.
    # ---------------------------------------------------------------------
    HR   ~ lnorm(expSd_HR)
    CTRM ~ lnorm(expSd_CTRM)
  })
}
