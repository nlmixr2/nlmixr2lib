Lu_2014_sglt_qsp <- function() {
  description <- paste(
    "QSP. Mechanistic systems pharmacology model of renal glucose",
    "reabsorption by SGLT1 and SGLT2 along the proximal tubules in",
    "humans, with optional competitive inhibition by an SGLT2",
    "inhibitor (calibrated to dapagliflozin; evaluated against",
    "canagliflozin). The proximal convoluted tubules (PCT) are",
    "divided into six sub-segments (PCT1-6, SGLT2-mediated",
    "reabsorption) and the proximal straight tubules into three",
    "(PST1-3, SGLT1-mediated). Filtrate drains into a urinary",
    "bladder. Plasma glucose (GLU, mmol/L) and plasma inhibitor",
    "(CINH, nmol/L) enter as time-varying regressors through",
    "glomerular filtration. Calibrated by hand-tuning in Berkeley",
    "Madonna v8.3.18 against the DeFronzo et al. (2013) urinary",
    "glucose excretion data; evaluated against Polidori et al.",
    "(2013), Mogensen (1971), and Wolf et al. (2009). 23 ODE states;",
    "no fitted IIV or residual error (typical-individual mechanism",
    "model fit to mean per-step data).",
    sep = " "
  )
  reference <- paste(
    "Lu Y, Griffen SC, Boulton DW, Leil TA.",
    "Use of systems pharmacology modeling to elucidate the operating",
    "characteristics of SGLT1 and SGLT2 in renal glucose reabsorption",
    "in humans.",
    "Front Pharmacol. 2014;5:274.",
    "doi:10.3389/fphar.2014.00274.",
    sep = " "
  )
  vignette <- "Lu_2014_sglt_qsp"
  paper_specific_compartments <- c("glu_pct1", "glu_pct2", "glu_pct3", "glu_pct4", "glu_pct5", "glu_pct6", "glu_pst1", "glu_pst2", "glu_pst3", "glu_bladder", "glu_urine", "glu_reabs", "drug_pct1", "drug_pct2", "drug_pct3", "drug_pct4", "drug_pct5", "drug_pct6", "drug_pst1", "drug_pst2", "drug_pst3", "drug_bladder", "drug_urine")

  units <- list(
    time          = "h",
    dosing        = "mmol",
    concentration = "mmol/L"
  )

  covariateData <- list(
    GLU = list(
      description        = "Plasma glucose concentration (time-varying regressor input)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying regressor; supplied at every observation/event",
        "row in the dataset and linearly interpolated between rows",
        "via the linear(GLU, CINH) declaration in model(). Drives the",
        "rate of glucose entry into PCT1 via glomerular filtration:",
        "filtered glucose load = GFR * GLU (mmol/h). The Lu 2014",
        "studies report plasma glucose in mg/dL; convert to the model's",
        "mmol/L units by dividing by 18.02 (1 mg/dL = 0.0555 mmol/L)."
      ),
      source_name        = "Plasma glucose"
    ),
    CINH = list(
      description        = "Plasma SGLT-inhibitor concentration (time-varying regressor input)",
      units              = "nmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying regressor; supplied at every observation/event",
        "row in the dataset and linearly interpolated between rows",
        "via the linear(GLU, CINH) declaration in model(). Drives the",
        "rate of unbound inhibitor entering PCT1 via glomerular",
        "filtration: filtered drug load = GFR * fup * CINH (nmol/h).",
        "Carries total plasma drug concentration (not unbound); the",
        "unbound fraction fup multiplies CINH inside model(). Set CINH",
        "= 0 for baseline (no-inhibitor) simulations. Conventional",
        "clinical reporting is ng/mL; convert by 1 ng/mL = (1000 /",
        "MW_drug) nmol/L (dapagliflozin MW = 409 g/mol so 1 ng/mL =",
        "2.44 nmol/L; canagliflozin MW = 454 g/mol so 1 ng/mL = 2.20",
        "nmol/L)."
      ),
      source_name        = "Plasma SGLT-inhibitor concentration"
    ),
    T2DM = list(
      description        = "Type-2 diabetes mellitus indicator (1 = T2DM, 0 = healthy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy)",
      notes              = paste(
        "Time-fixed at study entry. Multiplicatively shifts the typical",
        "Vmax2 (SGLT2 maximum reabsorption capacity) from the healthy",
        "reference value of 93.5 mmol/h to the diabetic value of 110",
        "mmol/h (+17.6% in T2DM), reflecting the calibrated up-regulation",
        "of SGLT2 expression in the diabetic state described in the Lu",
        "2014 Methods 'Model Parameters and Calibration' section. All",
        "other SGLT-kinetics parameters (Vmax1, Km1, Km2, Ki1, Ki2) are",
        "held common between healthy and T2DM per the paper. The",
        "canonical T2DM column carries 1 = T2DM patient, 0 = healthy",
        "control (inst/references/covariate-columns.md). Healthy",
        "subjects in DeFronzo et al. (2013) carry T2DM = 0; T2DM",
        "subjects in DeFronzo, Polidori, and Wolf carry T2DM = 1;",
        "Mogensen (1971) diabetics are coded as T2DM = 1 for this model",
        "even though Mogensen pre-dates the modern T1DM/T2DM",
        "classification (Lu 2014 does not distinguish them for the",
        "purpose of Vmax2 scaling)."
      ),
      source_name        = "Subject category (healthy vs diabetic)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 93L,
    n_studies      = 4L,
    age_range      = NA_character_,
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Pooled across four published clinical studies that used",
      "stepped hyperglycemic clamp (SHC) or fixed-elevated-glucose",
      "protocols to characterize SGLT-mediated renal glucose",
      "reabsorption: (1) DeFronzo et al. (2013), 12 healthy adults +",
      "12 T2DM, SHC at baseline and after 7 daily doses of 10 mg",
      "dapagliflozin (used for model calibration); (2) Polidori et al.",
      "(2013), 28 T2DM, SHC at baseline and after 8 daily doses of",
      "100 mg canagliflozin; (3) Mogensen (1971), 9 healthy + 10",
      "diabetic, glucose infusion to over 650 mg/dL; (4) Wolf et al.",
      "(2009), 22 T2DM, SHC. Studies 2-4 used for model evaluation."
    ),
    dose_range     = paste(
      "Oral SGLT2 inhibitors at clinical doses fed to the model as",
      "an exogenous plasma-concentration time course (CINH",
      "regressor): dapagliflozin 10 mg QD (DeFronzo et al. 2013);",
      "canagliflozin 100 mg QD (Polidori et al. 2013). The Lu 2014",
      "paper does not fit a PK model for either inhibitor: the",
      "dapagliflozin profile is the interpolated observed mean from",
      "DeFronzo et al. (2013), and the canagliflozin profile is a",
      "two-compartment PK model fitted to Devineni et al. (2013) mean",
      "data, parameters not reported."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Demographics for each study are in the cited primary",
      "publications; the Lu 2014 main text does not reproduce them.",
      "The model fits MEAN per-step data from each study (not",
      "individual-subject data), so per-subject demographics are not",
      "load-bearing for parameter inference. The packaged model is a",
      "typical-individual mechanism model (no IIV, no residual error)",
      "calibrated by hand-tuning in Berkeley Madonna v8.3.18 to the",
      "DeFronzo et al. (2013) cumulative and step-wise UGE data, then",
      "evaluated against three independent data sets. Calibration",
      "agreement is good in the clinically relevant 100-400 mg/dL",
      "plasma glucose range and weaker above 400 mg/dL (paper",
      "Discussion attributes the high-glucose deviation to",
      "unmodelled water reabsorption / hydrodynamic compensation)."
    )
  )

  ini({
    # =====================================================================
    # Anatomy / physiology (Lu 2014 Table 2). All values are literature or
    # assumed inputs to the calibration; they are exposed as estimable
    # typical values so downstream users can refit if needed.
    # =====================================================================

    lvctx  <- log(0.216) ; label("Renal cortex volume VCTX (L)")                                        # Lu 2014 Table 2, Thelwall et al. 2011
    fvptc  <- 0.3        ; label("Proximal tubule volume as fraction of renal cortex VPTC (unitless)")  # Lu 2014 Table 2, Moller and Skriver 1985
    fvpctc <- 0.7        ; label("PCT volume as fraction of proximal tubules VPCTC (unitless)")         # Lu 2014 Table 2, assumed
    lvx    <- log(0.2)   ; label("Urinary bladder volume VX (L)")                                       # Lu 2014 Table 2, Brown et al. 2011

    lgfr <- log(6.65)    ; label("Glomerular filtration rate GFR (L/h)")  # Lu 2014 Table 2 DeFronzo et al. 2013 healthy baseline range 5.66-7.38; midpoint 6.52
    lkx  <- log(0.83)    ; label("Urinary bladder voiding rate KX (L/h)") # Lu 2014 Table 2 DeFronzo et al. 2013 healthy baseline range 0.63-1.20; midpoint 0.92

    # =====================================================================
    # SGLT kinetics (Lu 2014 Table 2). Vmax1, Vmax2 (healthy), Km1, and Km2
    # are the calibrated final values from the DeFronzo et al. (2013) UGE
    # fit. Diabetic-state Vmax2 is reached via the T2DM covariate effect.
    # =====================================================================

    lvmax1 <- log(20.0)  ; label("SGLT1 maximum reabsorption rate Vmax1 (mmol/h)")                                  # Lu 2014 Table 2 calibrated value (PST1-3, SGLT1-mediated)
    lvmax2 <- log(93.5)  ; label("SGLT2 maximum reabsorption rate Vmax2 in healthy subjects (mmol/h)")              # Lu 2014 Table 2 calibrated value (PCT1-6, SGLT2-mediated)
    lkm1   <- log(0.5)   ; label("Glucose Michaelis-Menten constant Km1 for SGLT1 (mmol/L)")                        # Lu 2014 Table 2 calibrated value
    lkm2   <- log(4.0)   ; label("Glucose Michaelis-Menten constant Km2 for SGLT2 (mmol/L)")                        # Lu 2014 Table 2 calibrated value

    # =====================================================================
    # T2DM-vs-healthy adjustment for Vmax2. The paper reports Vmax2 = 110
    # mmol/h in T2DM and 93.5 mmol/h in healthy; the multiplicative shift
    # is 110/93.5 - 1 = 0.176 applied to the healthy reference.
    # =====================================================================

    e_t2dm_vmax2 <- 0.176 ; label("Multiplicative fractional shift on Vmax2 for T2DM = 1 (unitless)")  # Lu 2014 Table 2: Vmax2_T2DM = 110 mmol/h vs Vmax2_healthy = 93.5 mmol/h; ratio - 1 = 0.176

    # =====================================================================
    # Drug-specific parameters (Lu 2014 Table 3). Defaults below are for
    # dapagliflozin (the calibration drug); canagliflozin alternatives are
    # noted in the labels and can be substituted by the user.
    # =====================================================================

    lki1 <- log(400)     ; label("SGLT1 affinity Ki1 of inhibitor (nmol/L; dapagliflozin default = 400, canagliflozin = 200)")  # Lu 2014 Table 3 dapagliflozin from Hummel et al. 2011; canagliflozin from Grempler et al. 2012
    lki2 <- log(0.3)     ; label("SGLT2 affinity Ki2 of inhibitor (nmol/L; dapagliflozin default = 0.3, canagliflozin = 0.6)")  # Lu 2014 Table 3 dapagliflozin calibrated (Hummel et al. 2011 starting point 6 nM); canagliflozin from Grempler et al. 2012
    lfup <- log(0.07)    ; label("Unbound fraction in plasma fup of inhibitor (unitless; dapagliflozin default = 0.07, canagliflozin = 0.01)")  # Lu 2014 Table 3 dapagliflozin (asterisked literature value); canagliflozin from Devineni et al. 2013
  })

  model({
    # ---------------------------------------------------------------------
    # Time-varying regressor inputs. Both columns are required in the data
    # set:
    #   GLU  -- plasma glucose concentration in mmol/L
    #   CINH -- plasma SGLT-inhibitor concentration in nmol/L
    # rxode2 linearly interpolates each between adjacent dataset rows.
    # ---------------------------------------------------------------------
    linear(GLU, CINH)

    # ---------------------------------------------------------------------
    # Derived physiological volumes (L). The total proximal-tubule volume
    # is split equally among 6 PCT and 3 PST sub-segments per the Lu 2014
    # Methods 'Model Structure' paragraph.
    # ---------------------------------------------------------------------
    vctx          <- exp(lvctx)
    vx            <- exp(lvx)
    v_pt          <- fvptc * vctx
    v_pct         <- fvpctc * v_pt
    v_pst         <- (1 - fvpctc) * v_pt
    v_pct_subseg  <- v_pct / 6
    v_pst_subseg  <- v_pst / 3

    # ---------------------------------------------------------------------
    # Filtrate flow rates along the 9 tubular sub-segments (L/h). Per Lu
    # 2014 Table 2, KPCi / KPSj decrement linearly from 0.926 * GFR
    # (outflow of PCT1) to 0.333 * GFR (outflow of PST3) in steps of
    # 0.074 * GFR; this assumes uniform water reabsorption across the 9
    # sub-segments such that 2/3 of filtered water is reabsorbed by the
    # end of the proximal tubule (Koeppen and Stanton 2013).
    # ---------------------------------------------------------------------
    gfr  <- exp(lgfr)
    kx   <- exp(lkx)
    kpc1 <- 0.926 * gfr
    kpc2 <- 0.852 * gfr
    kpc3 <- 0.778 * gfr
    kpc4 <- 0.704 * gfr
    kpc5 <- 0.630 * gfr
    kpc6 <- 0.556 * gfr
    kps1 <- 0.482 * gfr
    kps2 <- 0.408 * gfr
    kps3 <- 0.333 * gfr

    # ---------------------------------------------------------------------
    # SGLT kinetics and drug-specific parameters. Vmax2 is the typical-
    # value healthy Vmax2; the multiplicative T2DM shift gives the
    # individual-state Vmax2 used in the PCT MM kinetics. Vmax is split
    # uniformly across the 6 PCT and 3 PST sub-segments per Lu 2014
    # Methods.
    # ---------------------------------------------------------------------
    vmax1     <- exp(lvmax1)
    vmax2_typ <- exp(lvmax2)
    vmax2     <- vmax2_typ * (1 + e_t2dm_vmax2 * T2DM)
    km1       <- exp(lkm1)
    km2       <- exp(lkm2)
    ki1       <- exp(lki1)
    ki2       <- exp(lki2)
    fup       <- exp(lfup)

    vmax2_sub <- vmax2 / 6
    vmax1_sub <- vmax1 / 3

    # ---------------------------------------------------------------------
    # Luminal concentrations per sub-segment. Glucose amounts are in mmol
    # so cglu_* is in mmol/L; drug amounts are in nmol so cdrug_* is in
    # nmol/L (matching the Ki units).
    # ---------------------------------------------------------------------
    cglu_pct1    <- glu_pct1    / v_pct_subseg
    cglu_pct2    <- glu_pct2    / v_pct_subseg
    cglu_pct3    <- glu_pct3    / v_pct_subseg
    cglu_pct4    <- glu_pct4    / v_pct_subseg
    cglu_pct5    <- glu_pct5    / v_pct_subseg
    cglu_pct6    <- glu_pct6    / v_pct_subseg
    cglu_pst1    <- glu_pst1    / v_pst_subseg
    cglu_pst2    <- glu_pst2    / v_pst_subseg
    cglu_pst3    <- glu_pst3    / v_pst_subseg
    cglu_bladder <- glu_bladder / vx

    cdrug_pct1    <- drug_pct1    / v_pct_subseg
    cdrug_pct2    <- drug_pct2    / v_pct_subseg
    cdrug_pct3    <- drug_pct3    / v_pct_subseg
    cdrug_pct4    <- drug_pct4    / v_pct_subseg
    cdrug_pct5    <- drug_pct5    / v_pct_subseg
    cdrug_pct6    <- drug_pct6    / v_pct_subseg
    cdrug_pst1    <- drug_pst1    / v_pst_subseg
    cdrug_pst2    <- drug_pst2    / v_pst_subseg
    cdrug_pst3    <- drug_pst3    / v_pst_subseg
    cdrug_bladder <- drug_bladder / vx

    # ---------------------------------------------------------------------
    # Per-sub-segment glucose reabsorption rates (Lu 2014 Eq. 2,
    # competitive-inhibition Michaelis-Menten; collapses to Eq. 1 when the
    # luminal inhibitor concentration cdrug = 0).
    #   PCT1-6 (SGLT2-mediated): Vmax2_sub, Km2, Ki2
    #   PST1-3 (SGLT1-mediated): Vmax1_sub, Km1, Ki1
    # Units: mmol/h.
    # ---------------------------------------------------------------------
    r_pct1 <- vmax2_sub * cglu_pct1 / (km2 * (1 + cdrug_pct1 / ki2) + cglu_pct1)
    r_pct2 <- vmax2_sub * cglu_pct2 / (km2 * (1 + cdrug_pct2 / ki2) + cglu_pct2)
    r_pct3 <- vmax2_sub * cglu_pct3 / (km2 * (1 + cdrug_pct3 / ki2) + cglu_pct3)
    r_pct4 <- vmax2_sub * cglu_pct4 / (km2 * (1 + cdrug_pct4 / ki2) + cglu_pct4)
    r_pct5 <- vmax2_sub * cglu_pct5 / (km2 * (1 + cdrug_pct5 / ki2) + cglu_pct5)
    r_pct6 <- vmax2_sub * cglu_pct6 / (km2 * (1 + cdrug_pct6 / ki2) + cglu_pct6)
    r_pst1 <- vmax1_sub * cglu_pst1 / (km1 * (1 + cdrug_pst1 / ki1) + cglu_pst1)
    r_pst2 <- vmax1_sub * cglu_pst2 / (km1 * (1 + cdrug_pst2 / ki1) + cglu_pst2)
    r_pst3 <- vmax1_sub * cglu_pst3 / (km1 * (1 + cdrug_pst3 / ki1) + cglu_pst3)

    # ---------------------------------------------------------------------
    # Glucose mass balance (amounts in mmol). PCT1 receives glomerular
    # filtrate at rate GFR * GLU; each downstream sub-segment receives the
    # outflow of the segment above. The bladder drains at rate KX, feeding
    # the cumulative urine compartment; reabsorbed glucose accumulates in
    # the (sink) glu_reabs compartment. All initial conditions default to
    # zero (rxode2 default); the system equilibrates within seconds of
    # simulation start given the small sub-segment volumes (~0.008 L) and
    # the high filtrate flow (~7 L/h).
    # ---------------------------------------------------------------------
    d/dt(glu_pct1) <- gfr * GLU         - kpc1 * cglu_pct1 - r_pct1
    d/dt(glu_pct2) <- kpc1 * cglu_pct1  - kpc2 * cglu_pct2 - r_pct2
    d/dt(glu_pct3) <- kpc2 * cglu_pct2  - kpc3 * cglu_pct3 - r_pct3
    d/dt(glu_pct4) <- kpc3 * cglu_pct3  - kpc4 * cglu_pct4 - r_pct4
    d/dt(glu_pct5) <- kpc4 * cglu_pct4  - kpc5 * cglu_pct5 - r_pct5
    d/dt(glu_pct6) <- kpc5 * cglu_pct5  - kpc6 * cglu_pct6 - r_pct6
    d/dt(glu_pst1) <- kpc6 * cglu_pct6  - kps1 * cglu_pst1 - r_pst1
    d/dt(glu_pst2) <- kps1 * cglu_pst1  - kps2 * cglu_pst2 - r_pst2
    d/dt(glu_pst3) <- kps2 * cglu_pst2  - kps3 * cglu_pst3 - r_pst3
    d/dt(glu_bladder) <- kps3 * cglu_pst3 - kx * cglu_bladder
    d/dt(glu_urine)   <- kx * cglu_bladder
    d/dt(glu_reabs)   <- r_pct1 + r_pct2 + r_pct3 + r_pct4 + r_pct5 + r_pct6 +
                         r_pst1 + r_pst2 + r_pst3

    # ---------------------------------------------------------------------
    # Inhibitor mass balance (amounts in nmol). PCT1 receives glomerular
    # filtrate of unbound drug at rate GFR * fup * CINH; the inhibitor is
    # not reabsorbed (Lu 2014 Methods: 'similar to glucose but without
    # tubular reabsorption') and flows through the sub-segments to the
    # bladder and out to urine.
    # ---------------------------------------------------------------------
    d/dt(drug_pct1) <- gfr * fup * CINH   - kpc1 * cdrug_pct1
    d/dt(drug_pct2) <- kpc1 * cdrug_pct1  - kpc2 * cdrug_pct2
    d/dt(drug_pct3) <- kpc2 * cdrug_pct2  - kpc3 * cdrug_pct3
    d/dt(drug_pct4) <- kpc3 * cdrug_pct3  - kpc4 * cdrug_pct4
    d/dt(drug_pct5) <- kpc4 * cdrug_pct4  - kpc5 * cdrug_pct5
    d/dt(drug_pct6) <- kpc5 * cdrug_pct5  - kpc6 * cdrug_pct6
    d/dt(drug_pst1) <- kpc6 * cdrug_pct6  - kps1 * cdrug_pst1
    d/dt(drug_pst2) <- kps1 * cdrug_pst1  - kps2 * cdrug_pst2
    d/dt(drug_pst3) <- kps2 * cdrug_pst2  - kps3 * cdrug_pst3
    d/dt(drug_bladder) <- kps3 * cdrug_pst3 - kx * cdrug_bladder
    d/dt(drug_urine)   <- kx * cdrug_bladder

    # ---------------------------------------------------------------------
    # Derived observables for vignette / downstream use:
    #   uge_rate -- instantaneous urinary glucose excretion rate (mmol/h)
    #   r_total  -- total renal glucose reabsorption rate (mmol/h)
    #   r_sglt1  -- SGLT1-mediated reabsorption rate (mmol/h)
    #   r_sglt2  -- SGLT2-mediated reabsorption rate (mmol/h)
    #   oe_sglt1 -- SGLT1 operation efficiency (% of Vmax1, Lu 2014 Eq.
    #               under 'SGLTs operation efficiency')
    #   oe_sglt2 -- SGLT2 operation efficiency (% of Vmax2)
    # The cumulative urine and reabsorbed-glucose amounts are available
    # directly as the ODE states glu_urine and glu_reabs.
    # ---------------------------------------------------------------------
    uge_rate <- kx * cglu_bladder
    r_total  <- r_pct1 + r_pct2 + r_pct3 + r_pct4 + r_pct5 + r_pct6 +
                r_pst1 + r_pst2 + r_pst3
    r_sglt1  <- r_pst1 + r_pst2 + r_pst3
    r_sglt2  <- r_pct1 + r_pct2 + r_pct3 + r_pct4 + r_pct5 + r_pct6
    oe_sglt1 <- 100 * r_sglt1 / vmax1
    oe_sglt2 <- 100 * r_sglt2 / vmax2
  })
}
