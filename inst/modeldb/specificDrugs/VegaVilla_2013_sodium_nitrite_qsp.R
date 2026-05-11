VegaVilla_2013_sodium_nitrite_qsp <- function() {
  description <- "QSP. Mechanistic systems pharmacology model of the NO metabolome (nitrite, nitrate) and methemoglobin (MetHb) in healthy adults receiving a 48-hour intravenous infusion of sodium nitrite. Nine ODEs covering plasma/RBC/tissue nitrite and nitrate, MetHb, NO and methemoglobin reductase activity; nonlinear nitrite/nitrate renal clearance (linear slope), entero-salivary nitrate-to-nitrite recycling, and indirect-response stimulation of MetHb reductase. Time in minutes; amounts in umol; concentrations in umol/L."
  reference <- "Vega-Villa K, Pluta R, Lonser R, Woo S. Quantitative Systems Pharmacology Model of NO Metabolome and Methemoglobin Following Long-Term Infusion of Sodium Nitrite in Humans. CPT Pharmacometrics Syst Pharmacol. 2013;2(8):e60. doi:10.1038/psp.2013.35"
  vignette <- "VegaVilla_2013_sodium_nitrite_qsp"
  units <- list(time = "min", dosing = "umol", concentration = "umol/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "21-56 years (mean 39, SD 9)",
    weight_range   = "49-115 kg (mean 77.8, SD 19)",
    sex_female_pct = NA_real_,
    disease_state  = "Healthy adult volunteers",
    dose_range     = "Sodium nitrite IV infusion 4.2-533.8 ug/kg/h for 48 h (escalating, single-dose-level per subject)",
    regions        = "United States",
    notes          = "Phase I dose-escalation cohort (Pluta et al. 2011 PLoS ONE 6:e14504; clinical protocol). Twelve subjects each receiving one of 4.2, 8.3, 16.7, 33.4, 66.8, 133.4, 266.9 (n=3, maximal-tolerated dose), 445.7 (n=2), or 533.8 (n=1) ug/kg/h IV for 48 h. Methods section cites Pluta 2011 for the demographics and study design. 333 plasma + 147 RBC nitrite/nitrate observations and 333 methemoglobin observations used for model development."
  )

  ini({
    # Structural rate constants (final estimates, Vega-Villa 2013 Table 1).
    # All times in minutes; concentrations in umol/L; second-order rate constants
    # in L/(umol*min) (reported as min^-1*umol^-1 in Table 1 with implicit
    # 1-L volume normalisation).

    # Inter-compartmental nitrite distribution. Supplement reparameterises
    # as Q/V1, Q/V4; Table 1 reports the derived rate constants directly,
    # which is what we use here.
    lkpt_no2 <- log(0.108);    label("Nitrite plasma-to-tissue rate constant (1/min)")  # Table 1
    lktp_no2 <- log(1.745);    label("Nitrite tissue-to-plasma rate constant (1/min)")  # Table 1

    # Inter-compartmental nitrate distribution (estimated directly, not via Q).
    lkpt_no3 <- log(0.160);    label("Nitrate plasma-to-tissue rate constant (1/min)")  # Table 1
    lktp_no3 <- log(0.515);    label("Nitrate tissue-to-plasma rate constant (1/min)")  # Table 1

    # Nitrate-to-nitrite recycling.
    lkno2_rdt <- log(1.65e-3); label("Nitrate-to-nitrite formation rate constant in plasma (entero-salivary; 1/min)")  # Table 1
    lkt_rdt   <- log(6.83e-5); label("Nitrate-to-nitrite formation rate constant in tissue (1/min)")                   # Table 1

    # Nitrite oxidation in tissue by myoglobin / heme carriers. Reported in
    # Table 1 as a stand-alone rate constant (0.234 1/min); supplement
    # derives it as kno3_r * baseline oxyHb. The Table 1 value is the
    # baseline typical value used here.
    lkmyo <- log(0.234);       label("Nitrite-to-nitrate formation rate constant in tissue (1/min, baseline typical)")  # Table 1

    # Renal/basal clearance and concentration-dependent slope factors
    # (Vega-Villa 2013 Eq. 5 and Eq. 6).
    lcl0_no2 <- log(0.382);    label("Basal clearance for endogenous nitrite (L/min)")  # Table 1
    ls_no2   <- log(4.524);    label("Slope factor on nitrite clearance (L/umol)")      # Table 1
    lvno2_p  <- log(12.418);   label("Volume of distribution of nitrite in central compartment (L)")  # Table 1

    lcl0_no3 <- log(9.41e-3);  label("Basal clearance for endogenous nitrate (L/min)")  # Table 1
    ls_no3   <- log(0.017);    label("Slope factor on nitrate clearance (L/umol)")      # Table 1
    lvno3_p  <- log(12.418);   label("Volume of distribution of nitrate in central compartment (L)")  # Table 1

    # Plasma <-> RBC exchange.
    lkpr_no2 <- log(0.019);    label("Nitrite uptake rate constant from plasma into RBC (1/min)")  # Table 1
    lkrp_no3 <- log(5.17e-3);  label("Nitrate diffusion rate constant from RBC into plasma (1/min)")  # Table 1

    # Nitrite-hemoglobin interactions in RBC (second-order).
    lkno3_r <- log(5.91e-6);   label("Nitrate+MetHb formation rate constant from nitrite (L/(umol*min))")  # Table 1
    lkno_r  <- log(1.42e-4);   label("NO+MetHb formation rate constant from nitrite (L/(umol*min))")      # Table 1

    # NO + Hb interactions in RBC. Table 1 reports kHbNO1 (NO+deoxyHb -> HbNO)
    # and kHbNO2 (NO+oxyHb -> MetHb+nitrate) as separate second-order rate
    # constants; the supplement parameterises them as a single kHbNO with a
    # fractional split FRHBNO, so kHbNO1 = kHbNO*(1-FRHBNO) and
    # kHbNO2 = kHbNO*FRHBNO. Table 1 numerical values imply FRHBNO ~ 1.4 %.
    # The paper's Discussion narrative quotes 4.61 % for this fraction; the
    # numerical Table 1 values are used here as the authoritative source.
    # See vignette Assumptions and deviations.
    lkhbno1 <- log(2.86e-5);   label("HbNO formation rate constant from NO+deoxyHb in RBC (L/(umol*min))")     # Table 1
    lkhbno2 <- log(4.07e-7);   label("MetHb+nitrate formation rate constant from NO+oxyHb in RBC (L/(umol*min))")  # Table 1

    # Methemoglobin reductase indirect-response model (Vega-Villa 2013 Eq. 12).
    lkdeg   <- log(0.016);     label("Degradation rate constant for methemoglobin reductase activity (1/min)")  # Table 1
    lsmethb <- log(0.045);     label("Slope factor on methemoglobin-induced reductase production (L/umol)")     # Table 1

    # RBC volumes of distribution.
    lvr_no2   <- log(5.816);   label("Volume of distribution of nitrite in RBC (L)")        # Table 1
    lvr_no3   <- log(4.391);   label("Volume of distribution of nitrate in RBC (L)")        # Table 1
    lvr_methb <- log(3.284);   label("Volume of distribution of methemoglobin in RBC (L)")  # Table 1

    # Direct plasma nitrite-to-nitrate conversion (supplement THETA(7) TVKNO3P,
    # fixed to zero in the final model; included for structural completeness).
    kno3_p <- fixed(0);        label("Direct plasma nitrite-to-nitrate conversion (1/min, fixed at 0)")  # supplement $THETA(7) FIX

    # ---- IIV (Table 1 %CV converted via omega^2 = log(1 + CV^2)) ----
    # Final %CV from Table 1, mapped to the corresponding nlmixr2 parameter.
    # Supplement parameterises kpt_no2/ktp_no2 via Q and V4 (correlated through
    # shared Q ETA); here we attach independent etas directly to the derived
    # kpt_no2 and ktp_no2, matching Table 1's marginal IIV reporting and
    # losing the implicit Q-mediated correlation. See vignette Assumptions
    # and deviations.
    etalkpt_no2  ~ 0.00807   # 9% CV (Table 1, on kPT_NO2 / supplement TVQ)
    etalktp_no2  ~ 1.0851    # 140% CV (Table 1, on kTP_NO2 / supplement TVV4)
    etalkno2_rdt ~ 1.4513    # 180% CV (Table 1, on kNO2_RDT / supplement KNO2RDT)
    etalvno2_p   ~ 0.5139    # 82% CV (Table 1, on VNO2_P / supplement V1)
    etalcl0_no3  ~ 1.1787    # 150% CV (Table 1, on CL(0)_NO3 / supplement CLR3)
    etalkpr_no2  ~ 0.7932    # 110% CV (Table 1, on kPR_NO2 / supplement KPRNO2)
    etalkrp_no3  ~ 0.6310    # 96% CV (Table 1, on kRP_NO3 / supplement KRPNO3)
    etalkno3_r   ~ 0.1907    # 46% CV (Table 1, on kNO3_R / supplement KNO3R)
    etalsmethb   ~ 1.0851    # 140% CV (Table 1, on SMetHb / supplement STIM)
    etalvr_no2   ~ 0.1284    # 37% CV (Table 1, on VR_NO2 / supplement VRBC_NO2)
    etalvr_methb ~ 0.2231    # 50% CV (Table 1, on VR_MethHb / supplement V5)

    # ---- Residual error (Vega-Villa 2013 Table 1, final estimates) ----
    # Combined additive + proportional on each endpoint; nitrate plasma and
    # RBC share a single pair per supplement $ERROR (same THETA(20)/THETA(21)
    # used for CMT 3 and CMT 4).
    addSd_Cc_nitrite_p  <- 0.024;  label("Additive residual SD on plasma nitrite (umol/L)")   # Table 1
    propSd_Cc_nitrite_p <- 0.5653; label("Proportional residual SD on plasma nitrite")        # Table 1, 56.53%
    addSd_Cc_nitrite_r  <- 0.217;  label("Additive residual SD on RBC nitrite (umol/L)")      # Table 1
    propSd_Cc_nitrite_r <- 0.4876; label("Proportional residual SD on RBC nitrite")           # Table 1, 48.76%
    # Supplement $ERROR shares THETA(20)/THETA(21) across plasma and RBC nitrate
    # (CMT 3 and CMT 4); nlmixr2 requires a unique residual-error parameter per
    # endpoint, so the same numerical values are declared twice.
    addSd_Cc_nitrate_p  <- 3.848;  label("Additive residual SD on plasma nitrate (umol/L)")   # Table 1 (shared)
    propSd_Cc_nitrate_p <- 0.3633; label("Proportional residual SD on plasma nitrate")        # Table 1, 36.33% (shared)
    addSd_Cc_nitrate_r  <- 3.848;  label("Additive residual SD on RBC nitrate (umol/L)")      # Table 1 (shared)
    propSd_Cc_nitrate_r <- 0.3633; label("Proportional residual SD on RBC nitrate")           # Table 1, 36.33% (shared)
    addSd_Cc_methb      <- 0.484;  label("Additive residual SD on methemoglobin (umol/L)")    # Table 1
    propSd_Cc_methb     <- 0.2614; label("Proportional residual SD on methemoglobin")         # Table 1, 26.14%
  })

  model({
    # ---- Constants ----
    # Typical adult total hemoglobin in the methemoglobin-RBC distribution
    # compartment (V5). Supplement $PK computes HBT = HBUM * V5 from the
    # per-subject covariate HBUM and the model's V5; the data column is not
    # available with the publication, so a typical-adult value is used here.
    # Chosen so that the supplement's baseline-typical kMYO = kno3_r * 0.77 *
    # HBT recovers Table 1's reported kMYO = 0.234 1/min when kno3_r equals
    # the Table 1 final (5.91e-6 L/(umol*min)). At this typical HBT and the
    # paper's 77/23 oxy/deoxy split, baseline oxyHb amount is ~39,580 umol
    # within the V5 compartment.
    hbt <- 51400  # umol total Hb in V5 compartment, typical adult (derived)

    # Typical adult plasma nitrite and methemoglobin baselines used as the
    # initial conditions of the ODE system. In the published fit these were
    # individual predose measurements read from the data file; for the
    # library model they are set to typical values consistent with the
    # paper's Results narrative (plasma nitrite baseline ~0.5 umol/L is
    # representative of Figure 4b predose, and ~1% MetHb on the typical
    # hbt above corresponds to ~22 umol/L in the V5 compartment).
    bl_nitrite_p_conc <- 0.5  # umol/L plasma nitrite baseline (Results, Figure 4b)
    bl_methb_conc     <- 22   # umol/L MetHb baseline (Results, ~1% of Hb total)

    # ---- Individual parameters ----
    kpt_no2  <- exp(lkpt_no2 + etalkpt_no2)
    ktp_no2  <- exp(lktp_no2 + etalktp_no2)
    kpt_no3  <- exp(lkpt_no3)
    ktp_no3  <- exp(lktp_no3)
    kno2_rdt <- exp(lkno2_rdt + etalkno2_rdt)
    kt_rdt   <- exp(lkt_rdt)
    kmyo     <- exp(lkmyo)
    cl0_no2  <- exp(lcl0_no2)
    s_no2    <- exp(ls_no2)
    vno2_p   <- exp(lvno2_p + etalvno2_p)
    cl0_no3  <- exp(lcl0_no3 + etalcl0_no3)
    s_no3    <- exp(ls_no3)
    vno3_p   <- exp(lvno3_p)
    kpr_no2  <- exp(lkpr_no2 + etalkpr_no2)
    krp_no3  <- exp(lkrp_no3 + etalkrp_no3)
    kno3_r   <- exp(lkno3_r + etalkno3_r)
    kno_r    <- exp(lkno_r)
    khbno1   <- exp(lkhbno1)
    khbno2   <- exp(lkhbno2)
    kdeg     <- exp(lkdeg)
    smethb   <- exp(lsmethb + etalsmethb)
    vr_no2   <- exp(lvr_no2 + etalvr_no2)
    vr_no3   <- exp(lvr_no3)
    vr_methb <- exp(lvr_methb + etalvr_methb)

    # ---- Baseline hemoglobin pools (computed once at t = 0) ----
    # Supplement $PK Eq. 13 and 14: oxyHb and deoxyHb are fractions of
    # (HBT - MetHb_amount). 77/23 split per Roberson 2012 (ref 45).
    bl_methb_amt   <- bl_methb_conc * vr_methb     # umol
    oxy_hb_0       <- (hbt - bl_methb_amt) * 0.77  # umol
    deoxy_hb_0     <- (hbt - bl_methb_amt) * 0.23  # umol

    # ---- Steady-state initial conditions (supplement $PK closed forms) ----
    bl_nitrite_p_amt <- bl_nitrite_p_conc * vno2_p
    bl_nitrite_r_amt <- kpr_no2 * bl_nitrite_p_amt /
      (kno_r * deoxy_hb_0 + kno3_r * oxy_hb_0)
    bl_no_r_amt <- kno_r * deoxy_hb_0 * bl_nitrite_r_amt /
      (khbno1 * deoxy_hb_0 + khbno2 * oxy_hb_0)
    bl_nitrate_r_amt <- (kno3_r * bl_nitrite_r_amt * oxy_hb_0 +
                         khbno2 * bl_no_r_amt * oxy_hb_0) / krp_no3

    kel_no2 <- cl0_no2 / vno2_p
    kel_no3 <- cl0_no3 / vno3_p
    ss_b <- kel_no3 + kno2_rdt + kpt_no3
    ss_d <- ktp_no2 + kmyo
    # Note: supplement E_ss has a KNO3P term that drops because KNO3P is
    # fixed at 0 in the final model.
    ss_e <- (kpt_no2 / ktp_no3 - ktp_no2 * kpt_no2 / (ktp_no3 * ss_d)) * bl_nitrite_p_amt
    ss_h <- kpt_no3 * krp_no3 / (ktp_no3 * ss_b) * bl_nitrate_r_amt
    bl_nitrate_t_amt <- (ss_e + ss_h) *
      (ss_b * ss_d * ktp_no3 /
        (ss_b * ss_d * ktp_no3 - (ss_d * ktp_no3 * kpt_no3 - ss_b * kt_rdt * ktp_no2)))
    bl_nitrite_t_amt <- (kpt_no2 * bl_nitrite_p_amt + kt_rdt * bl_nitrate_t_amt) /
      (ktp_no2 + kmyo)
    bl_nitrate_p_amt <- (krp_no3 * bl_nitrate_r_amt + ktp_no3 * bl_nitrate_t_amt) /
      (kno2_rdt + kpt_no3 + kel_no3)
    bl_kmr <- (kno_r * bl_nitrite_r_amt * deoxy_hb_0 +
               kno3_r * bl_nitrite_r_amt * oxy_hb_0 +
               khbno2 * bl_no_r_amt * oxy_hb_0) / bl_methb_amt

    # Endogenous zero-order nitrite production (Vega-Villa 2013 Eq. 7).
    kin_no2 <- kel_no2 * bl_nitrite_p_amt + kel_no3 * bl_nitrate_p_amt +
               khbno1 * bl_no_r_amt * deoxy_hb_0
    ksyn <- kdeg * bl_kmr

    # ---- Time-varying Hb pools ----
    oxy_hb   <- (hbt - methb) * 0.77
    deoxy_hb <- (hbt - methb) * 0.23

    # ---- Time-varying renal clearances (linear in concentration above baseline) ----
    nitrite_p_conc <- nitrite_p / vno2_p
    nitrate_p_conc <- nitrate_p / vno3_p
    cl_no2_t <- cl0_no2 * (1 + s_no2 * (nitrite_p_conc - bl_nitrite_p_conc))
    cl_no3_t <- cl0_no3 * (1 + s_no3 * (nitrate_p_conc - bl_nitrate_p_amt / vno3_p))

    # ---- Methemoglobin reductase stimulation ----
    methb_conc <- methb / vr_methb

    # ---- ODE system (Vega-Villa 2013 Eq. 1-12, supplement $DES; amounts in umol) ----
    d/dt(nitrite_p) <- kin_no2 + kno2_rdt * nitrate_p + ktp_no2 * nitrite_t -
                       kpt_no2 * nitrite_p - kno3_p * nitrite_p - kpr_no2 * nitrite_p -
                       cl_no2_t / vno2_p * nitrite_p
    d/dt(nitrite_r) <- kpr_no2 * nitrite_p - kno_r * nitrite_r * deoxy_hb -
                       kno3_r * nitrite_r * oxy_hb
    d/dt(nitrate_p) <- kno3_p * nitrite_p + krp_no3 * nitrate_r + ktp_no3 * nitrate_t -
                       kno2_rdt * nitrate_p - kpt_no3 * nitrate_p -
                       cl_no3_t / vno3_p * nitrate_p
    d/dt(nitrate_r) <- kno3_r * nitrite_r * oxy_hb - krp_no3 * nitrate_r +
                       khbno2 * no_r * oxy_hb
    d/dt(methb)     <- kno_r * nitrite_r * deoxy_hb + kno3_r * nitrite_r * oxy_hb -
                       kmr * methb + khbno2 * no_r * oxy_hb
    d/dt(nitrite_t) <- kpt_no2 * nitrite_p - ktp_no2 * nitrite_t +
                       kt_rdt * nitrate_t - kmyo * nitrite_t
    d/dt(nitrate_t) <- kpt_no3 * nitrate_p - ktp_no3 * nitrate_t -
                       kt_rdt * nitrate_t + kmyo * nitrite_t
    d/dt(kmr)       <- ksyn * (1 + smethb * (methb_conc - bl_methb_conc)) -
                       kdeg * kmr
    d/dt(no_r)      <- kno_r * nitrite_r * deoxy_hb - khbno2 * no_r * oxy_hb -
                       khbno1 * no_r * deoxy_hb

    # ---- Initial conditions ----
    nitrite_p(0) <- bl_nitrite_p_amt
    nitrite_r(0) <- bl_nitrite_r_amt
    nitrate_p(0) <- bl_nitrate_p_amt
    nitrate_r(0) <- bl_nitrate_r_amt
    methb(0)     <- bl_methb_amt
    nitrite_t(0) <- bl_nitrite_t_amt
    nitrate_t(0) <- bl_nitrate_t_amt
    kmr(0)       <- bl_kmr
    no_r(0)      <- bl_no_r_amt

    # ---- Observation variables (concentrations, umol/L) ----
    Cc_nitrite_p <- nitrite_p / vno2_p
    Cc_nitrite_r <- nitrite_r / vr_no2
    Cc_nitrate_p <- nitrate_p / vno3_p
    Cc_nitrate_r <- nitrate_r / vr_no3
    Cc_methb     <- methb / vr_methb

    # ---- Residual error ----
    Cc_nitrite_p ~ add(addSd_Cc_nitrite_p) + prop(propSd_Cc_nitrite_p)
    Cc_nitrite_r ~ add(addSd_Cc_nitrite_r) + prop(propSd_Cc_nitrite_r)
    Cc_nitrate_p ~ add(addSd_Cc_nitrate_p) + prop(propSd_Cc_nitrate_p)
    Cc_nitrate_r ~ add(addSd_Cc_nitrate_r) + prop(propSd_Cc_nitrate_r)
    Cc_methb     ~ add(addSd_Cc_methb)     + prop(propSd_Cc_methb)
  })
}
