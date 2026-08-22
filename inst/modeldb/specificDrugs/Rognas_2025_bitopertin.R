Rognas_2025_bitopertin <- function() {
  description <- paste(
    "QSP (semi-mechanistic erythropoiesis). Fifteen-state population PKPD model of human",
    "erythropoiesis and hemoglobin synthesis in healthy adults, driven by steady-state",
    "bitopertin (GlyT1 inhibitor) exposure. A bone-marrow erythroid precursor pool feeds a",
    "two-pathway four-compartment reticulocyte maturation structure (immature / mature x",
    "bone marrow / blood), which feeds a four-compartment erythrocyte age-transit chain with",
    "a parallel four-compartment corpuscular-hemoglobin chain. Total hemoglobin is the sum of",
    "the products across each erythrocyte / hemoglobin transit-compartment pair, and drives an",
    "exponential homeostatic feedback on both precursor recruitment and reticulocyte release.",
    "An empirical two-compartment moderator (tolerance) chain represents stem-cell reservoir",
    "depletion. Bitopertin inhibits the hemoglobin production rate via an Emax function of",
    "individual steady-state AUC. Simultaneous outputs: reticulocyte count, erythrocyte count,",
    "mean corpuscular hemoglobin, immature reticulocyte fraction, and total blood hemoglobin."
  )
  reference <- paste(
    "Rognas SV, Schaedeli Stark F, Marchesi M, Silber Baumann HE, Abrantes JA.",
    "A semi-mechanistic population pharmacokinetic-pharmacodynamic model to assess",
    "downstream drug-target effects on erythropoiesis.",
    "J Pharmacokinet Pharmacodyn. 2025;52(4):42. doi:10.1007/s10928-025-09990-7.",
    "Structural framework expanded from Schaedeli Stark F, Martin-Facklam M, Hofmann C,",
    "et al. (2012), semi-physiologic population PKPD model of bitopertin (RG1678) effects on",
    "hemoglobin turnover (reference 9 of the 2025 paper).",
    "Parameter values, the differential equation system and the steady-state parameterisation",
    "are taken from the 2025 paper and its Supplementary Appendix (Springer ESM MOESM1),",
    "which contains the full NONMEM control stream."
  )
  vignette <- "Rognas_2025_bitopertin"

  # The system time unit is days (control stream $INPUT: "TIME ; Days after first dose").
  # There is no drug compartment and no dosing event in this model: bitopertin enters only
  # through the AUC_BTP exposure covariate, so `dosing` is nominal. Cell counts are carried
  # in 10^9/L and corpuscular hemoglobin in pg, exactly as in the published control stream.
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SEXF = list(
      description        = "Sex indicator, 1 = female, 0 = male. Shifts the baseline erythrocyte count only.",
      units              = NA_character_,
      type               = "binary",
      reference_category = "male (SEXF = 0)",
      notes              = paste(
        "Control stream ($PK): TVRBC0 = (THETA(3) - THETA(10)**SEXF) * 1000, with the $INPUT",
        "column SEX coded 0 = male / 1 = female and re-derived as SEXF. Because THETA(10)**0 = 1,",
        "the male baseline is THETA(3) - 1 = 4.91 x 10^12/L (Table 2) and the female baseline is",
        "THETA(3) - THETA(10) = 4.32 x 10^12/L; the implied male-minus-female difference is",
        "THETA(10) - 1 = 0.59 x 10^12/L. Table 2 prints 1.59 for the row labelled",
        "RBCdiff,female, i.e. the raw THETA(10) rather than the difference its own footnote a",
        "defines. See the vignette Errata section.",
        "The only covariate effect retained in the final model."
      ),
      source_name        = "SEX"
    ),
    AUC_BTP = list(
      description        = paste(
        "Individual bitopertin steady-state area under the plasma concentration-time curve over",
        "the 24 h once-daily dosing interval (mg*h/L). Zero for placebo subjects and outside the",
        "treatment window."
      ),
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exposure metric driving the bitopertin drug effect (Methods, 'Individual AUCss values",
        "were used as the exposure metric driving the bitopertin drug effect'). In the control",
        "stream this is computed as DRUG = ADOS_BTP / ICL_BTP * TREAT_BTP -- the actual last dose",
        "divided by the empirical-Bayes apparent clearance from a bitopertin population PK model",
        "that the paper does not report ('data not shown'), multiplied by a treatment indicator.",
        "TREAT_BTP is 1 from the time of a dose until two days after the last dose and 0 otherwise,",
        "so the drug effect switches off two days after treatment stops; downstream users encode",
        "that carry-over in the AUC_BTP time course. Step-wise time-varying. The paper's own",
        "simulations (Fig. 5) specify exposure as multiples of AUC50 rather than as doses:",
        "0.5 x AUC50 = 8.25 mg*h/L gives 20% and 2 x AUC50 = 33 mg*h/L gives 40% inhibition."
      ),
      source_name        = "DRUG"
    )
  )

  compartmentData <- list(
    precursor1 = list(
      analyte = "erythroid precursor cell", units = "10^9/L",
      specimen = "tissue", verified = TRUE
    ),
    ret_imm_marrow = list(
      analyte = "immature reticulocyte", units = "10^9/L",
      specimen = "tissue", verified = TRUE
    ),
    ret_mat_marrow = list(
      analyte = "mature reticulocyte", units = "10^9/L",
      specimen = "tissue", verified = TRUE
    ),
    ret_imm_blood = list(
      analyte = "immature reticulocyte", units = "10^9/L",
      specimen = "blood cell", verified = TRUE
    ),
    ret_mat_blood = list(
      analyte = "mature reticulocyte", units = "10^9/L",
      specimen = "blood cell", verified = TRUE
    ),
    erythrocytes1 = list(
      analyte = "erythrocyte", units = "10^9/L",
      specimen = "blood cell", verified = TRUE
    ),
    erythrocytes2 = list(
      analyte = "erythrocyte", units = "10^9/L",
      specimen = "blood cell", verified = TRUE
    ),
    erythrocytes3 = list(
      analyte = "erythrocyte", units = "10^9/L",
      specimen = "blood cell", verified = TRUE
    ),
    erythrocytes4 = list(
      analyte = "erythrocyte", units = "10^9/L",
      specimen = "blood cell", verified = TRUE
    ),
    mch1 = list(
      analyte = "hemoglobin per erythrocyte", units = "pg",
      specimen = "blood cell", verified = TRUE
    ),
    mch2 = list(
      analyte = "hemoglobin per erythrocyte", units = "pg",
      specimen = "blood cell", verified = TRUE
    ),
    mch3 = list(
      analyte = "hemoglobin per erythrocyte", units = "pg",
      specimen = "blood cell", verified = TRUE
    ),
    mch4 = list(
      analyte = "hemoglobin per erythrocyte", units = "pg",
      specimen = "blood cell", verified = TRUE
    ),
    moderator1 = list(
      analyte = "erythroid precursor cell", units = "10^9/L",
      specimen = "not applicable", verified = TRUE
    ),
    moderator2 = list(
      analyte = "erythroid precursor cell", units = "10^9/L",
      specimen = "not applicable", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 62,
    n_studies      = 1,
    age_range      = "19-45 years (all subjects aged < 50 years by protocol)",
    age_median     = "31.0 years",
    weight_range   = "49.4-108 kg",
    weight_median  = "69.0 kg",
    sex_female_pct = 58,
    disease_state  = "Healthy male and female volunteers",
    dose_range     = "Bitopertin 10, 30 or 60 mg orally once daily for 120 days, or placebo, followed by a 120-day follow-up period",
    notes          = paste(
      "Phase 1 multicentre, randomised, double-blind, placebo-controlled, parallel-group trial.",
      "67 subjects enrolled; 62 included in the model (placebo n = 15, 10 mg n = 17, 30 mg n = 16,",
      "60 mg n = 14) -- Results 'Data used for model development' and Table 1.",
      "Baseline biomarkers (Table 1, all subjects, mean): reticulocyte count 40.2 x 10^9/L,",
      "erythrocyte count 4.70 x 10^12/L, total hemoglobin 140.4 g/L, mean corpuscular hemoglobin",
      "29.9 pg, immature reticulocyte fraction 0.052.",
      "Hematology sampled at baseline, weeks 1 and 2, then every 2 weeks to week 16, week 17",
      "(treatment end), week 18, then every 2 weeks to week 34.",
      "Fitted simultaneously to reticulocyte count, erythrocyte count, mean corpuscular hemoglobin",
      "and immature reticulocyte fraction; total hemoglobin observations were deliberately EXCLUDED",
      "from the fit ($INFN sets MDV = 1 for TYPE == 7; Table 1 footnote a) and are a pure",
      "out-of-sample prediction (Fig. 4).",
      "NONMEM 7.4.3, SAEM with a preceding ITS step and an IMP evaluation step; condition number",
      "272 (Table 2 footnote b). All IIV shrinkage below 25%.",
      "Age and body weight were carried in the dataset but no age or weight effect was retained."
    )
  )

  ini({
    # ---- Cell lifespans (mean transit times) -------------------------------
    # LS_PRE fixed to a literature value; LS_RBC estimated.
    lmtt_pre <- fixed(log(5))     ; label("Erythroid precursor transit time in bone marrow LS_PRE (day)")   # Table 2: 5 days, Fixed (literature refs 10, 11); control stream $PK LS_PRE = 5
    lmtt_rbc <- log(125)          ; label("Erythrocyte lifespan / mean transit time LS_RBC (day)")           # Table 2: 125 days (RSE 4.44%); $THETA(1) init 124

    # ---- Baselines (steady-state anchors) ----------------------------------
    lrbase_ret <- log(39.80)      ; label("Baseline total blood reticulocyte count RET0 (10^9/L)")           # Table 2: 39.80 x 10^9/L (RSE 3.37%); $THETA(11) init 39.9
    lrbase_rbc <- log(4.91)       ; label("Baseline erythrocyte count in males RBC0,male (10^12/L)")         # Table 2: 4.91 x 10^12/L (RSE 0.78%); = THETA(3) - 1 per $PK parameterisation
    lrbase_mch <- log(29.80)      ; label("Baseline mean corpuscular hemoglobin MCH0 (pg)")                  # Table 2: 29.80 pg (RSE 0.61%); $THETA(2) init 29.8
    lrbase_irf <- log(0.0471)     ; label("Baseline immature reticulocyte fraction IRF0 (fraction)")         # Table 2: 4.71% (RSE 4.31%); $THETA(14) init 0.045, bounded (0, 0.5)

    # ---- Covariate effect --------------------------------------------------
    e_sexf_rbase_rbc <- 0.59      ; label("Female minus male shift in baseline erythrocyte count (10^12/L, subtracted)")  # Table 2 RBCdiff,female row prints THETA(10) = 1.59 (RSE 3.18%); $PK TVRBC0 = (THETA(3) - THETA(10)**SEXF) * 1000 so the actual difference is THETA(10) - 1 = 0.59 (Table 2 footnote a definition)

    # ---- Homeostatic feedback and empirical tolerance ----------------------
    lfb   <- log(2.42)            ; label("Hemoglobin feedback amplification factor (unitless exponent multiplier)")  # Table 2 Feedback: 2.42 (RSE 8.84%); functional form in Table 2 footnote b / $DES HB_STIM
    lktol <- log(0.022)           ; label("Moderator (tolerance) chain turnover rate constant kTOL (1/day)")          # Table 2: 0.022 /day (RSE 12.86%); $THETA(17) init 0.0234

    # ---- Bitopertin exposure-response -------------------------------------
    limax  <- fixed(log(0.6))     ; label("Maximum fractional inhibition of hemoglobin synthesis by bitopertin Imax (fraction)")  # Table 2: 0.6, Fixed to a previously estimated value (Results: 'fixed to a previously estimated value to increase model stability'); $THETA(5) (0,0.6,1) FIX
    lauc50 <- log(16.50)          ; label("Bitopertin steady-state AUC giving half-maximal inhibition AUC50 (mg*h/L)")             # Table 2: 16.50 mg/L x h (RSE 8.91%); $THETA(6) init 16.3
    lhill  <- fixed(log(1))       ; label("Hill coefficient on the bitopertin exposure-response (unitless)")                       # $THETA(7) 1 FIX (GAM_BTP); not reported in Table 2 because it was held at 1

    # ---- Drug-target mechanism switches -----------------------------------
    # Not estimated parameters. The paper uses this model as a framework and
    # simulates four hypothetical drug-target interaction mechanisms (Methods
    # "Simulations", Fig. 2 boxes A / B / C / D1 / D2, Fig. 5). The defaults
    # below reproduce bitopertin itself, i.e. Mechanism A. Override via ini()
    # to place the Emax inhibition on a different pathway:
    #   Mechanism A (bitopertin): mech_mch 1, others 0                     -- inhibition of hemoglobin synthesis rate Rin,MCH
    #   Mechanism B:              mech_kinpre 1, others 0                  -- inhibition of reticulocyte precursor recruitment Rin,PRE
    #   Mechanism C:              mech_ktrpre 1, others 0                  -- inhibition of precursor-to-reticulocyte differentiation kPRE
    #   Mechanism D:              mech_mch 1 and mech_fb 1                 -- D2 inhibition of Rin,MCH plus D1 full inhibition of the Hbtot-driven feedback
    mech_mch    <- fixed(1)       ; label("Switch: place the inhibition on the hemoglobin input rate Rin,MCH (Mechanism A / D2; 1 = on)")   # Methods 'Simulations' Mechanism A, 'resembling the effect observed for bitopertin'; $DES applies INH_BTP to DADT(10) only
    mech_kinpre <- fixed(0)       ; label("Switch: place the inhibition on the precursor recruitment rate Rin,PRE (Mechanism B; 1 = on)")   # Methods 'Simulations' Mechanism B
    mech_ktrpre <- fixed(0)       ; label("Switch: place the inhibition on the precursor differentiation rate kPRE (Mechanism C; 1 = on)")  # Methods 'Simulations' Mechanism C
    mech_fb     <- fixed(0)       ; label("Switch: fully inhibit the hemoglobin-driven feedback while on treatment (Mechanism D1; 1 = on)") # Methods 'Simulations' Mechanism D, 'full inhibition of Hbtot-driven feedback'

    # ---- Between-subject variability ---------------------------------------
    # Table 2 reports IIV as %CV; converted to log-normal variance as
    # omega^2 = log(1 + CV^2). The paper estimated a full $OMEGA BLOCK(6) over
    # these six parameters but reports only the diagonal, so the packaged model
    # is diagonal -- see the vignette Errata section.
    etalmtt_rbc   ~ 0.0755824  # Table 2 IIV 28.02 %CV (shrinkage 16.92%)
    etalrbase_ret ~ 0.0655106  # Table 2 IIV 26.02 %CV (shrinkage 1.88%)
    etalrbase_rbc ~ 0.0028475  # Table 2 IIV 5.34 %CV (shrinkage 1.63%)
    etalrbase_mch ~ 0.0023205  # Table 2 IIV 4.82 %CV (shrinkage 0.47%)
    etalrbase_irf ~ 0.0980127  # Table 2 IIV 32.09 %CV (shrinkage 3.88%)
    etalauc50     ~ 0.2191556  # Table 2 IIV 49.50 %CV (shrinkage 24.97%)

    # ---- Residual unexplained variability ----------------------------------
    # NONMEM $ERROR uses W = sqrt((prop * IPRED)^2 + add^2) with $SIGMA 1 FIX,
    # i.e. the combined2 additive-plus-proportional form nlmixr2 uses by default.
    propSd_RET <- 0.22          ; label("Proportional residual error for reticulocyte count (fraction)")            # Table 2 eps(RET, proportional): 0.22 (RSE 2.35%); $THETA(12). $THETA(13) additive RET term is 0 FIX
    addSd_RBC  <- 0.18          ; label("Additive residual error for erythrocyte count (10^12/L)")                  # Table 2 eps(RBC, additive): 0.18 (RSE 2.22%); $THETA(8)
    addSd_MCH  <- 0.35          ; label("Additive residual error for mean corpuscular hemoglobin (pg)")             # Table 2 eps(MCH, additive): 0.35 (RSE 2.30%); $THETA(9)
    addSd_IRF  <- 0.0096        ; label("Additive residual error for immature reticulocyte fraction (fraction)")    # Table 2 eps(IRF, additive): 0.96% (RSE 16.48%); $THETA(15) init 0.00996
    propSd_IRF <- 0.38          ; label("Proportional residual error for immature reticulocyte fraction (fraction)")# Table 2 eps(IRF, proportional): 0.38 (RSE 5.41%); $THETA(16)
    addSd_thb  <- 0.53          ; label("Additive residual error for total blood hemoglobin (g/L)")                 # $ERROR for dv_hb: W = THETA(8) + THETA(9) = 0.18 + 0.35; hemoglobin was NOT fitted (MDV = 1) so this drives the Fig. 4 prediction interval only
  })

  model({
    # ------------------------------------------------------------------
    # Source trace. Equation numbers "Eq. N" are the main paper; "Suppl.
    # Eq. N" are the differential equation system in Supplementary
    # Appendix section 2; "$PK" / "$DES" / "$ERROR" are the NONMEM control
    # stream in Supplementary Appendix section 4.
    #
    # Cell states are carried in 10^9/L and corpuscular hemoglobin in pg,
    # matching the control stream. Reported erythrocyte counts (10^12/L)
    # and hemoglobin concentrations (g/L) are obtained by dividing by 1000
    # in the observation block, exactly as $ERROR does.
    # ------------------------------------------------------------------

    # 1. Individual parameters
    mttPre <- exp(lmtt_pre)
    mttRbc <- exp(lmtt_rbc + etalmtt_rbc)
    mch0   <- exp(lrbase_mch + etalrbase_mch)
    ret0   <- exp(lrbase_ret + etalrbase_ret)
    irf0   <- exp(lrbase_irf + etalrbase_irf)
    ktol   <- exp(lktol)
    fb     <- exp(lfb)
    imax   <- exp(limax)
    auc50  <- exp(lauc50 + etalauc50)
    hill   <- exp(lhill)

    # Baseline erythrocyte count: sex shift applied to the typical value
    # before the log-normal IIV, per $PK (RBC0 = EXP(LOG(TVRBC0) + ETA(3))).
    # Converted from 10^12/L to the 10^9/L working unit (the "* 1000" of $PK).
    rbc0 <- (exp(lrbase_rbc) - e_sexf_rbase_rbc * SEXF) * 1000 * exp(etalrbase_rbc)

    # 2. Transit-rate constants and steady-state-derived quantities
    nctr <- 4                                   # $PK NCTR = 4 (number of erythrocyte transit compartments)
    ktrRbc <- nctr / mttRbc                     # Eq. 3: kRBC = nCTR / LS_RBC
    ktrPre <- 1 / mttPre                        # Eq. 1: LS_PRE = 1 / kPRE

    # Probability that a reticulocyte is released from the bone marrow while
    # still immature. Suppl. Appendix section 3 Eq. 3: p_release = IRF0 / (1 - IRF0),
    # derived from the steady-state immature / mature blood balance; $PK writes the
    # algebraically identical p_release = -IRF0 / (IRF0 - 1).
    prelease <- irf0 / (1 - irf0)

    # Eq. 2: reticulocyte maturation / release rate constant, fixed by the
    # steady-state requirement that mature blood reticulocytes supply the
    # erythrocyte chain.
    ktrRet <- ktrRbc * (rbc0 / nctr) / (ret0 * (1 - irf0))

    # Zero-order input rates and the precursor baseline follow from homeostasis
    # at t = 0 ($PK RIN_PRE / PRE0 / RIN_MCH).
    kinPre <- ktrRbc * (rbc0 / nctr)
    pre0   <- kinPre / ktrPre
    kinMch <- ktrRbc * mch0
    thb0   <- mch0 * rbc0                       # $PK HB0 = MCH0 * RBC0

    # 3. Total hemoglobin, homeostatic feedback and drug effect
    # Eq. 4: Hbtot = sum over the four erythrocyte / MCH transit-compartment pairs.
    thbNow <- erythrocytes1 * mch1 + erythrocytes2 * mch2 +
      erythrocytes3 * mch3 + erythrocytes4 * mch4

    # Exponential stimulation driven by the fractional DECREASE of total
    # hemoglobin from baseline, acting on precursor recruitment and on
    # reticulocyte release from the bone marrow ($DES HB_CFB / HB_STIM_PRE /
    # HB_STIM_IRF, both computed with the same Feedback parameter -- Results:
    # "Estimating the feedback components as two separate parameters did not
    # significantly improve the fit"). NOTE: Table 2 footnote b prints the
    # exponent numerator as (Hbtot - Hb0); the control stream and the paper's
    # own prose ("a feedback mechanism increases precursor recruitment when
    # Hbtot decreased") both require (Hb0 - Hbtot). See vignette Errata.
    # Mechanism D1 fully switches the feedback off while on treatment; the
    # default mech_fb = 0 leaves it intact for bitopertin.
    trt <- (AUC_BTP > 0)
    stim <- exp((thb0 - thbNow) / thb0 * fb * (1 - mech_fb * trt))

    # Bitopertin inhibition of the hemoglobin production rate ($DES INH_BTP,
    # with GAM_BTP fixed to 1). AUC_BTP already carries the treatment-window
    # gating that $PK applies through TREAT_BTP. The mechanism switches route
    # the same Emax inhibition to the pathway being interrogated.
    inhFrac <- imax * AUC_BTP^hill / (auc50^hill + AUC_BTP^hill)
    inh     <- 1 - mech_mch * inhFrac      # on Rin,MCH  (Mechanism A / D2, bitopertin)
    inhPre  <- 1 - mech_kinpre * inhFrac   # on Rin,PRE  (Mechanism B)
    inhKpre <- 1 - mech_ktrpre * inhFrac   # on kPRE     (Mechanism C)

    # 4. ODE system (Suppl. Eq. 1-15 / $DES DADT(1)-DADT(15))
    # Mechanism C scales the precursor-to-reticulocyte differentiation rate;
    # the same scaled rate must appear on both sides of the PRE / RET_imm,bm
    # transfer so cell mass is conserved.
    ktrPreEff <- ktrPre * inhKpre

    # Precursors, modulated by the moderator chain (Suppl. Eq. 1)
    d/dt(precursor1) <- kinPre * inhPre * stim * (pre0 / moderator2) -
      ktrPreEff * precursor1

    # Empirical tolerance / moderator chain (Suppl. Eq. 14-15). Driven by the
    # precursor amount with NO cell transfer (Fig. 2 blue dotted line).
    d/dt(moderator1) <- ktol * (precursor1 - moderator1)
    d/dt(moderator2) <- ktol * (moderator1 - moderator2)

    # Reticulocytes in bone marrow (Suppl. Eq. 2-3). Suppl. Eq. 2 as typeset
    # carries k_LS,PRE on the immature-release term; the control stream
    # ($DES DADT(2)) uses KLS_RET there, which is what mass balance with
    # DADT(4) requires. The control stream form is used. See vignette Errata.
    d/dt(ret_imm_marrow) <- ktrPreEff * precursor1 -
      ktrRet * stim * prelease * ret_imm_marrow -
      ktrRet * (1 - prelease) * ret_imm_marrow
    d/dt(ret_mat_marrow) <- ktrRet * (1 - prelease) * ret_imm_marrow -
      ktrRet * stim * ret_mat_marrow

    # Reticulocytes in blood (Suppl. Eq. 4-5)
    d/dt(ret_imm_blood) <- ktrRet * stim * prelease * ret_imm_marrow -
      ktrRet * ret_imm_blood
    d/dt(ret_mat_blood) <- ktrRet * stim * ret_mat_marrow +
      ktrRet * ret_imm_blood - ktrRet * ret_mat_blood

    # Erythrocyte age-transit chain (Suppl. Eq. 6-9)
    d/dt(erythrocytes1) <- ktrRet * ret_mat_blood - ktrRbc * erythrocytes1
    d/dt(erythrocytes2) <- ktrRbc * (erythrocytes1 - erythrocytes2)
    d/dt(erythrocytes3) <- ktrRbc * (erythrocytes2 - erythrocytes3)
    d/dt(erythrocytes4) <- ktrRbc * (erythrocytes3 - erythrocytes4)

    # Parallel corpuscular-hemoglobin chain; bitopertin acts on the input rate
    # only, i.e. on hemoglobin loaded into newly produced cells (Suppl. Eq. 10-13).
    # The model assumes no loss of hemoglobin from an erythrocyte during aging.
    d/dt(mch1) <- kinMch * inh - ktrRbc * mch1
    d/dt(mch2) <- ktrRbc * (mch1 - mch2)
    d/dt(mch3) <- ktrRbc * (mch2 - mch3)
    d/dt(mch4) <- ktrRbc * (mch3 - mch4)

    # 5. Steady-state (homeostatic) initial conditions ($PK A_0 block)
    precursor1(0)     <- pre0
    moderator1(0)     <- pre0
    moderator2(0)     <- pre0
    ret_imm_marrow(0) <- ret0 * (irf0 / prelease)
    ret_mat_marrow(0) <- ret0 * (1 - 2 * irf0)
    ret_imm_blood(0)  <- ret0 * irf0
    ret_mat_blood(0)  <- ret0 * (1 - irf0)
    erythrocytes1(0)  <- rbc0 / nctr
    erythrocytes2(0)  <- rbc0 / nctr
    erythrocytes3(0)  <- rbc0 / nctr
    erythrocytes4(0)  <- rbc0 / nctr
    mch1(0)           <- mch0
    mch2(0)           <- mch0
    mch3(0)           <- mch0
    mch4(0)           <- mch0

    # 6. Observations ($ERROR). Unit conversions match $ERROR exactly.
    RET <- ret_imm_blood + ret_mat_blood                          # 10^9/L
    rbcTot <- erythrocytes1 + erythrocytes2 + erythrocytes3 + erythrocytes4
    RBC <- rbcTot / 1000                                          # 10^12/L
    MCH <- thbNow / rbcTot                                        # pg
    IRF <- ret_imm_blood / RET                                    # fraction
    thb <- thbNow / 1000                                          # g/L

    RET ~ prop(propSd_RET)
    RBC ~ add(addSd_RBC)
    MCH ~ add(addSd_MCH)
    IRF ~ add(addSd_IRF) + prop(propSd_IRF)
    thb ~ add(addSd_thb)
  })
}
