Suzuki_2024_mycophenolic_acid <- function() {
  description <- paste0(
    "Population pharmacokinetic model for mycophenolic acid (MPA), the active ",
    "metabolite of mycophenolate mofetil (MMF), in 42 adult renal transplant ",
    "recipients at a single Japanese centre (312 plasma MPA observations from ",
    "14 days pre- to 84 days post-transplantation; MMF 500-1000 mg twice ",
    "daily). The paper ",
    "compared six absorption structures plus an enterohepatic-circulation ",
    "(EHC) model; a two-compartment disposition model with SEQUENTIAL ",
    "zero- then first-order absorption and a lag time gave the lowest OFV, ",
    "and the EHC process did not improve the fit. Dose enters the depot as a ",
    "zero-order input over D1 = 1.95 h beginning after a lag of 0.283 h, then ",
    "drains first-order with ka = ln(2)/TABS (absorption half-life TABS = ",
    "0.0453 h, so ka = 15.3 /h and absorption is rate-limited by D1). ",
    "Apparent total clearance is split into an additive non-renal arm ",
    "(7.56 L/h) plus a renal arm (11.6 L/h) that scales linearly with renal ",
    "function RF = CLcr / 100, where CLcr is a Cockcroft-Gault estimate ",
    "standardised to a 70 kg body weight. Allometric size scaling uses the ",
    "Holford normal-fat-mass (NFM) framework, but Ffat was estimated as not ",
    "different from 1 and fixed there for both clearance and volume, so NFM ",
    "collapses exactly to total body weight and the size model reduces to ",
    "(WT/70)^0.75 on CL/F and Q/F and (WT/70)^1 on VC/F and VP/F (the ",
    "paper's own Equation 8). Post-transplantation day acts on relative ",
    "bioavailability as F1 = exp(0.956 * POD / 84) for POD >= 0 (1 otherwise), ",
    "a 2.8-fold rise in exposure by 90 days. The MMF-to-MPA molecular-weight ",
    "ratio 320.3/433.5 is carried inside F1, so doses are given as mg of MMF. ",
    "Between-subject variability is exponential on CL/F, VC/F, TABS, lag time ",
    "and D1, and additive-proportional on F1; BSV on Q/F and VP/F was fixed ",
    "to zero. Residual error is combined proportional (28.3%) plus additive ",
    "(0.0634 mg/L)."
  )
  reference <- paste(
    "Suzuki Y, Matsunaga N, Aoyama T, Ogami C, Hasegawa C, Iida S, To H,",
    "Kitahara T, Tsuji Y. (2024). Population pharmacokinetic analysis",
    "identifies an absorption process model for mycophenolic acid in patients",
    "with renal transplant. Clin Transl Sci 17(12):e70097.",
    "doi:10.1111/cts.70097",
    sep = " "
  )
  vignette <- "Suzuki_2024_mycophenolic_acid"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Suzuki 2024 Appendix S1 (NM-TRAN
  # control stream): ADVAN4 TRANS4, S2 = V2 with "MPA dosing is mg, DV is
  # mg/L", and F1 carries the MMF -> MPA molecular-weight conversion, so
  # every state holds an MPA-equivalent amount in mg.
  compartmentData <- list(
    depot       = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (TBW).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Allometric size scaling with reference TBW = 70 kg; exponent fixed ",
        "at 0.75 on CL/F and Q/F and at 1 on VC/F and VP/F (Suzuki 2024 ",
        "Methods 'Covariate model' Equation 1, Results Equation 8). The ",
        "paper's size descriptor is formally normal fat mass (NFM = FFM + ",
        "Ffat * FAT, Equation 2, with FFM from the Janmahasatian equation, ",
        "Equation 3), but Ffat for both clearance and volume was estimated ",
        "as not significantly different from 1 and fixed at 1 (Table 2, ",
        "Results 'Covariate model'; Appendix S1 'FFATCL' / 'FFATV' = 1 FIX). ",
        "With Ffat = 1, NFM = FFM + (TBW - FFM) = TBW and NFMstd = 56.1 + ",
        "(70 - 56.1) = 70 kg, so the NFM machinery collapses exactly to ",
        "TBW / 70 and height and sex drop out of the final model. Equation 8 ",
        "is written by the authors in the collapsed TBW form and that is what ",
        "is encoded here. Cohort TBW 2.5th / 50th / 97.5th percentiles ",
        "39.6 / 55.9 / 77.5 kg (Table 1)."
      ),
      source_name        = "TBW"
    ),
    CRCL = list(
      description        = paste0(
        "Creatinine clearance estimated by the Cockcroft-Gault equation with ",
        "body weight FIXED at 70 kg, i.e. standardised per 70 kg rather than ",
        "per 1.73 m^2 BSA. Appendix S1: CLCR = (140 - AGE) / SCR * 70 / 72 * ",
        "FCPR with FCPR = 1 for men and 0.85 for women, SCR in mg/dL and AGE ",
        "in years."
      ),
      units              = "mL/min/70 kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters the model only through the renal-function ratio RF = CRCL / ",
        "CLcr_std with CLcr_std = 100 mL/min/70 kg = 6 L/h/70 kg (Suzuki 2024 ",
        "Equation 4), which multiplies the renal clearance arm (Equation 5). ",
        "NOTE the normalisation: this is mL/min per 70 kg of body weight, NOT ",
        "the BSA-normalised mL/min/1.73 m^2 that the CRCL register entry ",
        "names as its default, and NOT a raw Cockcroft-Gault value. Body size ",
        "is handled separately and independently by the allometric WT term, ",
        "so substituting a raw or BSA-normalised CLcr here would double-count ",
        "size. Cohort 2.5th / 50th / 97.5th percentiles 6.6 / 52.8 / 148.8 ",
        "mL/min/70 kg (Table 1); the simulations in Figure 3 use 25, 50, 75 ",
        "and 100 mL/min/70 kg. Time-varying within subject as graft function ",
        "recovers."
      ),
      source_name        = "CLcr"
    ),
    POD = list(
      description        = paste0(
        "Post-transplantation day: days elapsed since renal transplantation. ",
        "The paper calls this PTD."
      ),
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-varying within subject. Enters relative bioavailability only, ",
        "as F_PTD = exp(K_PTD * POD / POD_max) with POD_max = 84 days, the ",
        "maximum PTD in the analysis dataset (Suzuki 2024 Equation 6, ",
        "Results Equation 8). NEGATIVE VALUES OCCUR AND ARE GATED: the ",
        "dataset spans 14 days PRE-transplantation to 84 days ",
        "post-transplantation (Results 'PopPK of MPA'), and Appendix S1 ",
        "applies 'IF (PTD.GE.0) THEN FPTD = EXP(KPTD * PTD / 84) ELSE FPTD = ",
        "1', so any pre-transplant record takes F_PTD = 1. That gate is ",
        "reproduced here; do not feed a negative POD into the exponential ",
        "directly. POD_max = 84 is retained as the divisor even when ",
        "simulating beyond 84 days (the paper itself simulates to 90 days in ",
        "Figure 3 and quotes the resulting 2.8-fold rise), so the effect ",
        "extrapolates rather than saturating. Cohort 2.5th / 50th / 97.5th ",
        "percentiles 23 / 40 / 84 days (Table 1)."
      ),
      source_name        = "PTD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    n_observations = "312 plasma MPA concentrations (Suzuki 2024 Results 'PopPK of MPA')",
    age_range      = "29.0-68.9 years (2.5th-97.5th percentile; median 51.0)",
    age_median     = "51.0 years",
    weight_range   = "39.6-77.5 kg (2.5th-97.5th percentile; median 55.9)",
    weight_median  = "55.9 kg",
    height_range   = "1.50-1.80 m (2.5th-97.5th percentile; median 1.66)",
    sex_female_pct = 31.0,
    race_ethnicity = "not reported; single-centre Japanese cohort",
    disease_state  = "adult kidney transplant recipients receiving mycophenolate mofetil immunosuppression",
    renal_function = "creatinine clearance 6.6-148.8 mL/min/70 kg (2.5th-97.5th percentile; median 52.8), i.e. spanning normal to severely impaired graft function",
    dose_range     = "mycophenolate mofetil 500 mg (n = 10), 750 mg (n = 25) or 1000 mg (n = 7) orally twice daily, every 12 h",
    regions        = "Japan (Nagasaki University Hospital, April 2011 to September 2019)",
    co_medication  = "tacrolimus (n = 24, median trough 6.9 ug/L) and prednisone (n = 23, median 10 mg/day) in some patients; neither improved the OFV as a covariate on CL/F or VC/F and neither is in the final model",
    notes          = paste0(
      "Baseline demographics from Suzuki 2024 Table 1 (reported as 2.5th / ",
      "50th / 97.5th percentiles rather than mean +/- SD). Sex 29 male / 13 ",
      "female. Steady-state samples; MPA measured by PETINIA immunoassay ",
      "(Dimension Xpand) with LLOQ 0.2 mg/L and no observations below LLOQ. ",
      "Fitted in NONMEM 7.5.1 with FOCE-I via Wings for NONMEM; internal ",
      "validation by 200-sample bootstrap and prediction-corrected VPC. Only ",
      "TOTAL MPA was modelled: unbound MPA and the MPAG metabolite were ",
      "excluded because the data were unavailable (Discussion, limitations), ",
      "so this model cannot describe the EHC contribution to exposure."
    )
  )

  ini({
    # =========================================================================
    # Structural PK. All values are Suzuki 2024 Table 2 'Final model estimate'
    # and are cross-confirmed by the $THETA block of the NM-TRAN control stream
    # in Appendix S1. Clearances and volumes are APPARENT (per unit
    # bioavailability) and are reported at the standard TBW of 70 kg.
    # =========================================================================
    lcl_nonren <- log(7.56)
    label("Apparent non-renal clearance CL_nonrenal/F at TBW = 70 kg (L/h)")
    # Suzuki 2024 Table 2: 7.56 L/h (%RSE 25.4; bootstrap 95% CI 5.36-13.53).
    # Appendix S1 $THETA '(0,7.56,) ; POP_CLNR L/h, Non-renal clearance'.

    lcl_renal <- log(11.6)
    label("Apparent renal clearance CL_renal/F at TBW = 70 kg and RF = 1 (L/h)")
    # Suzuki 2024 Table 2: 11.6 L/h (%RSE 30.0; bootstrap 95% CI 3.8-18.0).
    # Appendix S1 $THETA '(0,11.6,) ; POP_CLR L/h, Renal clearance'. Enters
    # multiplied by RF (Equation 5), so 11.6 L/h is the renal arm of a patient
    # at the standard CLcr of 100 mL/min/70 kg.

    lvc <- log(104.0)
    label("Apparent central volume of distribution VC/F at TBW = 70 kg (L)")
    # Suzuki 2024 Table 2: 104.0 L (%RSE 18.9; bootstrap 95% CI 68.1-143.0).
    # Appendix S1 $THETA '(0,104.,) ; POP_V2 L'.

    lq <- log(17.3)
    label("Apparent inter-compartmental clearance Q/F at TBW = 70 kg (L/h)")
    # Suzuki 2024 Table 2: 17.3 L/h (%RSE 31.3; bootstrap 95% CI 10.1-35.2).
    # Appendix S1 $THETA '(0,17.3,) ; POP_Q L/h'.

    lvp <- log(169.0)
    label("Apparent peripheral volume of distribution VP/F at TBW = 70 kg (L)")
    # Suzuki 2024 Table 2: 169.0 L (%RSE 89.1; bootstrap 95% CI 68.7-1060.0).
    # Appendix S1 $THETA '(0,169.,) ; POP_V3 L'. The bootstrap distribution is
    # heavily right-skewed (mean 279 vs median 195), i.e. the peripheral volume
    # is the least well determined parameter in the model.

    # -------------------------------------------------------------------------
    # Absorption: sequential zero-order (over D1, after a lag) then first-order.
    # -------------------------------------------------------------------------
    ltabs <- log(0.0453)
    label("Absorption half-life TABS (h); ka = ln(2)/TABS")
    # Suzuki 2024 Table 2: 0.0453 h (%RSE 67.6; bootstrap 95% CI 0.0102-0.1285).
    # Appendix S1 $THETA '(0,0.0453,) ; POP_TABS h, Absorption half-life' and
    # $PK 'KA=LOG(2)/TABS'. Equation 8 prints the same relation as
    # ka (h^-1) = 0.693 / 0.0453 = 15.3 /h. The first-order step is therefore
    # very fast and absorption is rate-limited by the zero-order duration D1.
    # Parameterising on TABS (rather than directly on ka) is deliberate: the
    # paper places the between-subject variability on TABS.

    ltlag <- log(0.283)
    label("Absorption lag time ALAG1 (h)")
    # Suzuki 2024 Table 2: 0.283 h (%RSE 26.4; bootstrap 95% CI 0.145-0.475).
    # Appendix S1 $THETA '(0,0.283,) ; POP_ALAG1 h, Lag time'. Equation 8
    # mislabels the units of this line as '(L)'; Table 2 and the control
    # stream both give h, and a lag time is a time.

    ld1 <- log(1.95)
    label("Duration of the zero-order absorption input D1 (h)")
    # Suzuki 2024 Table 2: 1.95 h (%RSE 15.7; bootstrap 95% CI 1.26-2.27).
    # Appendix S1 $THETA '(0,1.95,) ; POP_D1 h'. Equation 8 mislabels the units
    # of this line as '(L)' as well; Table 2 and the control stream give h.

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F1 at POD = 0, before the MMF-to-MPA molar conversion (unitless)")
    # Suzuki 2024 Table 2: 'Relative bioavailability = 1 fixed'.
    # Appendix S1 $THETA '1. FIX ; POP_F1, Relative bioavailability'.

    # -------------------------------------------------------------------------
    # Covariate effects.
    # -------------------------------------------------------------------------
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of (WT/70) on CL/F and Q/F (unitless)")
    # Suzuki 2024 Methods 'Covariate model': "the allometric exponent (PWR) of
    # F_SIZE was fixed at 0.75 for CL/F (and Q/F)". Appendix S1
    # 'FSIZCL=(NFMCL/NFMCLSTD)**(3/4)'; Equation 8 prints the exponent 3/4 on
    # both CL/F and Q/F.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of (WT/70) on VC/F and VP/F (unitless)")
    # Suzuki 2024 Methods 'Covariate model': "... and at 1 for VC/F (and
    # VP/F)". Appendix S1 'FSIZV=(NFMV/NFMVSTD)**1'; Equation 8 prints VC/F and
    # VP/F scaled by (TBW/TBWstd) to the first power.

    e_pod_fdepot <- 0.956
    label("K_PTD, exponential coefficient of (POD/84) on relative bioavailability (unitless)")
    # Suzuki 2024 Table 2 'Effect of posttransplantation day on relative
    # bioavailability' = 0.956 (%RSE 26.8; bootstrap 95% CI 0.529-1.773).
    # Appendix S1 $THETA '(0,0.956) ; KPTD'. Discussion cross-check: relative
    # bioavailability is "predicted to be 2.8-fold higher at 90 days
    # posttransplantation than on the day of transplantation", and
    # exp(0.956 * 90 / 84) = 2.79.

    # =========================================================================
    # Between-subject variability. Suzuki 2024 Table 2 reports these as "CV%".
    # That CV% is 100 * sqrt(omega^2) -- the SD of the eta -- NOT the exact
    # lognormal CV 100 * sqrt(exp(omega^2) - 1). This is settled, not assumed:
    # the Appendix S1 $OMEGA block gives the variances directly, and
    # 100 * sqrt(omega^2) reproduces all six published CV% values exactly
    # (0.183 -> 42.8, 0.337 -> 58.1, 0.534 -> 73.1, 1.48 -> 121.7,
    # 0.825 -> 90.8, 0.0324 -> 18.0), whereas the lognormal formula does not
    # (it would give 44.8, 63.3, 84.0, 184.2, 113.2, 18.1). The $OMEGA
    # variances are therefore used verbatim below.
    #
    # BSV on Q/F and VP/F was fixed to zero (Table 2 '0 fixed'; Appendix S1
    # '$OMEGA BLOCK(1) 0. FIX ; PPV_Q' and '... PPV_V3'), so those etas are
    # omitted rather than declared with zero variance -- a zero-variance
    # diagonal makes OMEGA singular and breaks the Cholesky sampler.
    #
    # No IIV correlations are reported: Appendix S1 declares eight separate
    # $OMEGA BLOCK(1) statements, i.e. a diagonal OMEGA.
    # =========================================================================
    etalcl ~ 0.183
    # Appendix S1 '$OMEGA BLOCK(1) 0.183 ; PPV_CL'. Table 2 'Apparent
    # clearance 42.8 CV%' = 100 * sqrt(0.183). A SINGLE eta acts on the
    # COMPOSITE apparent clearance, not on the renal and non-renal arms
    # separately: Appendix S1 $PK builds GRPCL = (GRPCLNR + GRPCLR) * FSIZCL
    # and only then applies CL = GRPCL * EXP(PPV_CL).

    etalvc ~ 0.337
    # Appendix S1 '$OMEGA BLOCK(1) 0.337 ; PPV_V2'. Table 2 'Apparent central
    # volume of distribution 58.1 CV%' = 100 * sqrt(0.337).

    etaltabs ~ 0.534
    # Appendix S1 '$OMEGA BLOCK(1) 0.534 ; PPV_TABS'. Table 2 'Absorption
    # half-life 73.1 CV%' = 100 * sqrt(0.534).

    etaltlag ~ 1.48
    # Appendix S1 '$OMEGA BLOCK(1) 1.48 ; PPV_ALAG1'. Table 2 'Lag time
    # 121.7 CV%' = 100 * sqrt(1.48). The largest BSV in the model; the paper
    # attributes it to patients showing a markedly delayed tmax (Figure 1,
    # Discussion paragraph 3).

    etald1 ~ 0.825
    # Appendix S1 '$OMEGA BLOCK(1) 0.825 ; PPV_D1'. Table 2 'Duration of
    # zero-order absorption 90.8 CV%' = 100 * sqrt(0.825).

    etalfdepot ~ 0.0324
    # Appendix S1 '$OMEGA BLOCK(1) 0.0324 ; PPV_F1'. Table 2 'Relative
    # bioavailability 18.0 CV%' = 100 * sqrt(0.0324). NOTE that this eta is
    # ADDITIVE-PROPORTIONAL on the natural scale, not log-additive: Appendix S1
    # $PK uses 'F1=GRPF1*(1+PPV_F1)', not EXP(PPV_F1), and it is applied in that
    # literal form in model() below. The name still carries the 'l' of its
    # parameter 'lfdepot', per the package convention that an IIV term is named
    # eta + the transformed parameter name; the same pattern (a proportional
    # eta named after the transformed parameter) appears in
    # Goulooze_2022_finerenone.R.

    # =========================================================================
    # Residual unexplained variability -- combined proportional + additive.
    # Appendix S1 $ERROR: PROP = CTOTAL * RUV_PROP; ADD = RUV_ADD;
    # SD = SQRT(PROP*PROP + ADD*ADD); Y = CTOTAL + SD*EPS1 with $SIGMA 1 FIX.
    # That is exactly nlmixr2's `prop() + add()` combination, and both THETAs
    # are on the SD scale (the additive term is reported in mg/L in Table 2).
    # =========================================================================
    propSd <- 0.283
    label("Proportional residual error (fraction)")
    # Suzuki 2024 Table 2: 28.3 CV% (%RSE 7.1; bootstrap 95% CI 22.9-30.0).
    # Appendix S1 $THETA '(0,0.283,) ; RUV_PROP'.

    addSd <- 0.0634
    label("Additive residual error (mg/L)")
    # Suzuki 2024 Table 2: 0.0634 mg/L (%RSE 79.0; bootstrap 95% CI
    # 0.0138-0.3053). Appendix S1 $THETA '(0,0.0634,) ; RUV_ADD mg/L'. Well
    # below the 0.2 mg/L assay LLOQ, consistent with a dataset in which no
    # observation fell below LLOQ.
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Derived covariate terms.
    # -----------------------------------------------------------------------
    # Renal function as a ratio to the standard CLcr of 100 mL/min/70 kg
    # (= 6 L/h/70 kg). Suzuki 2024 Equation 4; Appendix S1 'RF=CLCR/100'.
    rf <- CRCL / 100

    # Post-transplantation-day effect on relative bioavailability. Suzuki 2024
    # Equation 6 with PTD_max = 84 days. The (POD > 0) factor reproduces the
    # Appendix S1 gate 'IF (PTD.GE.0) THEN FPTD=EXP(KPTD*PTD/84) ELSE FPTD=1':
    # for a pre-transplant record (POD < 0) the exponent is forced to 0 and
    # F_PTD becomes 1, and at POD = 0 the exponent is 0 either way.
    f_ptd <- exp(e_pod_fdepot * POD * (POD > 0) / 84)

    # Allometric size factors. Suzuki 2024 Equation 1 with the NFM descriptor
    # of Equations 2-3; because Ffat is fixed at 1 for both clearance and
    # volume, NFM = TBW and NFM_std = 70 kg exactly, which is why Equation 8
    # is written directly in terms of TBW/TBW_std.
    fsize_cl <- (WT / 70)^e_wt_cl_q
    fsize_v  <- (WT / 70)^e_wt_vc_vp

    # -----------------------------------------------------------------------
    # 2. Individual PK parameters.
    # -----------------------------------------------------------------------
    # Additive renal + non-renal clearance decomposition, Suzuki 2024
    # Equation 5: CL_total/F = CL_nonrenal/F + CL_renal/F * RF. The single
    # composite eta and the size factor are applied to the SUM, per
    # Appendix S1 $PK.
    cl_nonren <- exp(lcl_nonren)
    cl_renal  <- exp(lcl_renal)
    cl <- (cl_nonren + cl_renal * rf) * fsize_cl * exp(etalcl)

    vc <- exp(lvc + etalvc) * fsize_v
    q  <- exp(lq)  * fsize_cl
    vp <- exp(lvp) * fsize_v

    # Absorption. ka is derived from the absorption half-life, which is where
    # the paper places the BSV (Appendix S1 'KA=LOG(2)/TABS').
    tabs <- exp(ltabs + etaltabs)
    ka   <- log(2) / tabs
    tlag <- exp(ltlag + etaltlag)
    d1   <- exp(ld1 + etald1)

    # Relative bioavailability. Three things multiply together here, exactly as
    # in Appendix S1 'GRPF1=POP_F1*MWMPA*FPTD' followed by
    # 'F1=GRPF1*(1+PPV_F1)':
    #   (a) the fixed F1 of 1,
    #   (b) the MMF -> MPA molecular-weight ratio 320.3 / 433.5 = 0.7389, which
    #       converts an administered mycophenolate mofetil dose in mg into the
    #       equivalent mg of mycophenolic acid (Suzuki 2024 Methods 'PopPK':
    #       "the units of MMF doses were normalized to account for their
    #       molecular weight ratios with MPA (MMF and MPA: 433.5 and 320.3
    #       g/mol)"; MMF is rapidly and essentially completely hydrolysed to
    #       MPA presystemically). Doses supplied to this model are therefore mg
    #       of MMF, not mg of MPA.
    #   (c) the post-transplantation-day factor.
    # The eta is additive-proportional, not exponential.
    fdepot <- exp(lfdepot) * (320.3 / 433.5) * f_ptd * (1 + etalfdepot)

    # -----------------------------------------------------------------------
    # 3. Two-compartment disposition with a first-order depot (ADVAN4 TRANS4).
    # -----------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          cl / vc * central -
                          q / vc * central +
                          q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central -
                          q / vp * peripheral1

    # -----------------------------------------------------------------------
    # 4. Sequential zero- then first-order absorption with a lag. The dose is
    #    delivered into the depot as a zero-order input spread uniformly over
    #    d1, starting tlag after the dose record, and the depot then drains
    #    first-order at ka. Dose records must use rate = -2 so that rxode2
    #    honours the modelled duration.
    # -----------------------------------------------------------------------
    dur(depot) <- d1
    lag(depot) <- tlag
    f(depot)   <- fdepot

    # -----------------------------------------------------------------------
    # 5. Observation. Appendix S1 $ERROR: CTOTAL = A(2)/V2 with S2 = V2.
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
