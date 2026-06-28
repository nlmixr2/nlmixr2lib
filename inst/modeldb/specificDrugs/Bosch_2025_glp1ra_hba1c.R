Bosch_2025_glp1ra_hba1c <- function() {
  description <- paste(
    "QSP. Integrated glucose-red blood cell-HbA1c (IGRH) sub-model from",
    "the Bosch 2025 4GI-HbA1c systems framework, used to predict",
    "long-term HbA1c response from a time-varying plasma glucose driver",
    "in adults with type 2 diabetes mellitus receiving GLP-1R / GLP-1R +",
    "glucagon receptor agonists (cotadutide, liraglutide). The model is",
    "a 24-state transit chain (12 unglycated red blood cell age cohorts",
    "+ 12 glycated cohorts; NC = 12 transit compartments per Bosch 2025",
    "supplement model code S2) with a glucose-concentration-dependent",
    "shortening of the RBC life span; the HbA1c output is the percentage",
    "glycated fraction of the total RBC pool. All structural parameters",
    "are fixed from the Lledo-Garcia 2013 / Kjellsson 2015 IGRH",
    "publications and held constant during the Bosch 2025 calibration;",
    "only the residual error and the IIV on the RBC life span were",
    "estimated on the cotadutide Ph2a HbA1c dataset (Bosch 2025 Table 2,",
    "third-block 'IGRH model'). Plasma glucose drive is supplied as the",
    "time-varying regressor GLU in mmol/L (linearly interpolated by",
    "rxode2 between dataset rows) and the per-subject baseline glucose",
    "anchor FPG is in mmol/L; both are converted to mg/dL inside model()",
    "to match the published IGRH parameterisation (KG in dL/mg/day,",
    "reference glucose 149 mg/dL = 8.27 mmol/L).",
    sep = " "
  )
  reference <- paste(
    "Bosch R, Petrone M, Arends R, Sijbrands EJG, Hoefman S, Snelder N.",
    "From In Vitro Efficacy to Long-Term HbA1c Response for GLP-1R /",
    "GlucagonR Agonism Using the 4GI-HbA1c Systems Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2025.",
    "doi:10.1002/psp4.70074.",
    "IGRH sub-model structure inherited from Lledo-Garcia R, Kjellsson MC,",
    "Karlsson MO. Br J Clin Pharmacol. 2013;76(2):301-311.",
    "doi:10.1111/bcp.12089.",
    sep = " "
  )
  vignette <- "Bosch_2025_glp1ra_hba1c"
  paper_specific_compartments <- c(
    "hba_n01", "hba_n02", "hba_n03", "hba_n04", "hba_n05", "hba_n06",
    "hba_n07", "hba_n08", "hba_n09", "hba_n10", "hba_n11", "hba_n12",
    "hba_g01", "hba_g02", "hba_g03", "hba_g04", "hba_g05", "hba_g06",
    "hba_g07", "hba_g08", "hba_g09", "hba_g10", "hba_g11", "hba_g12"
  )

  units <- list(
    time          = "day",
    dosing        = "(none; glucose-driven, no exogenous dosing)",
    concentration = "% (HbA1c percentage)"
  )

  covariateData <- list(
    GLU = list(
      description        = "Plasma glucose concentration (time-varying regressor input)",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying regressor; supplied at every observation / event",
        "row in the dataset and linearly interpolated between rows via",
        "the linear(GLU) declaration in model(). Drives the glycation",
        "rate (KG * GLU_mgdl) and the glucose-dependent shortening of",
        "the RBC life span (LS(t) = TVLS * (GLU_mgdl / 149)^gamma_ls).",
        "Bosch 2025 Methods Section 2.3 derives GLU as the predicted",
        "average daily glucose concentration (Cglc,av) from the 4GI",
        "model output; the IGRH model in the source paper consumes",
        "Cglc,av directly. For self-contained simulations, supply any",
        "time-varying glucose trajectory (cotadutide-treated, placebo,",
        "etc.) in mmol/L. The model converts to mg/dL internally as",
        "GLU * 18.02 to match the original Lledo-Garcia 2013 / Bosch",
        "2025 supplement S2 parameterisation."
      ),
      source_name        = "Cglc,av (Bosch 2025 Methods Section 2.3)"
    ),
    FPG = list(
      description        = "Baseline (per-subject) fasting plasma glucose anchor",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject; supplied at study entry. Used to",
        "initialise the transit chain in steady state at the subject's",
        "baseline glycemic state (Bosch 2025 supplement S2 model code:",
        "A_0(2..25) closed-form steady-state expressions driven by AG).",
        "Converted to mg/dL inside model() as FPG * 18.02. Bosch 2025",
        "Ph2a population median FPG = 9.6 mmol/L (SD 2.2); Ph2b median",
        "FPG = 10.1 mmol/L (SD 2.6) per supplement Table S1."
      ),
      source_name        = "AG / IAG / IBGLC (Bosch 2025 supplement S2 $INPUT)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = NA_character_,
    age_median     = "71 years (SD 6.1)",
    weight_range   = NA_character_,
    weight_median  = "93.2 kg (SD 16.6)",
    sex_female_pct = 48.8,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Overweight or obese adults with type 2 diabetes mellitus, baseline",
      "HbA1c 7.87% (SD 1.0), baseline FPG 9.6 mmol/L (SD 2.2), median",
      "duration of T2DM 16.1 years (SD 7.8) (Bosch 2025 supplement Table",
      "S1, Ph2a column). The Ph2a HbA1c dataset on which the IGRH IIV /",
      "residual error were re-estimated is the cotadutide Ph2a study",
      "D5670C00011 / NCT03244800 (Cohort 1 placebo + Cohort 1 cotadutide",
      "n = 39 with weekly HbA1c sampling on Days 7, 15, 49)."
    ),
    dose_range     = paste(
      "Cotadutide once-daily subcutaneous, 50-300 ug titration (Bosch",
      "2025 Table 1, Ph2a Cohort 1: 50 ug for 4 days, 100 for 4 days,",
      "200 for 7 days, 300 for 28 days). Cotadutide PK is NOT modelled",
      "in this IGRH sub-model; the only input is the time-varying",
      "plasma glucose driver GLU produced by an upstream model (4GI in",
      "the source paper) or supplied externally."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Bosch 2025 fits the IGRH HbA1c sub-model sequentially after the",
      "4GI glucose / insulin / GLP-1 / glucagon / GIP module: the 4GI",
      "produces a time series of Cglc,av which is then fed into IGRH as",
      "the AG / CAG inputs. This packaged model captures the IGRH layer",
      "only; the upstream 4GI parameterisation is a separate extraction",
      "(Bosch_2024_cotadutide_qsp). All four IGRH structural parameters",
      "(kg, TVLSP, TVLS, gamma_ls) are FIXED at the published",
      "Lledo-Garcia 2013 / Kjellsson 2015 values; only the residual",
      "error (proportional SD 0.0311) and the IIV on the RBC life span",
      "(omega^2 = 0.00471, ~6.9% CV) are estimated. Bosch 2025 Table 2",
      "third-block reports omega^2 = 0.0145 with 16.5% RSE for the IIV",
      "on LS, but the supplement S2 NONMEM $OMEGA block AND the paper",
      "Results Section 3.1 text (6.9% CV, 18.4% RSE) are mutually",
      "consistent at omega^2 = 0.00471; the Table 2 entry is treated as",
      "a typographical discrepancy and the supplement / text value is",
      "used (see in-file source-trace comment on etals)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # IGRH structural parameters (Bosch 2025 Table 2 'IGRH model' block).
    # All four are FIXED at the Lledo-Garcia 2013 / Kjellsson 2015 values;
    # Bosch 2025 Methods Section 2.3 'In line with Kjellsson et al. [8]'
    # and Results Section 3.1 'Initial model exploration showed that all
    # system parameters could be fixed based on published values'.
    # ---------------------------------------------------------------------

    lkg <- fixed(log(8.37e-6))
    label("Glycation rate constant kg (log dL/mg/day)")
    # Bosch 2025 Table 2: kg = 8.37e-6 dL/mg/day (FIXED, no RSE reported).
    # Supplement S2 NONMEM $THETA(1): (0, 0.008370, 1) FIX; KG = THETA(1)/1000.

    llsp <- fixed(log(8.2))
    label("RBC precursor pool life span LSP (log days)")
    # Bosch 2025 Table 2: TVLSP = 8.2 days (FIXED, no RSE).
    # Supplement S2 NONMEM $THETA(4): (0, 8.2) FIX.

    lls <- fixed(log(91.7))
    label("RBC life span LS at reference glucose (log days)")
    # Bosch 2025 Table 2: TVLS = 91.7 days (FIXED, no RSE).
    # Supplement S2 NONMEM $THETA(3): (0, 91.7) FIX.

    gamma_ls <- fixed(-0.381)
    label("Shape of glucose-concentration-adjusted LS (unitless)")
    # Bosch 2025 Table 2: gammaLS = -0.381 (FIXED, no RSE).
    # Supplement S2 NONMEM $THETA(2): (-0.381) FIX.
    # Negative exponent: LS decreases as glucose rises above 149 mg/dL.

    # ---------------------------------------------------------------------
    # Inter-individual variability on the RBC life span. The published
    # Bosch 2025 IIV on LSP did NOT improve the fit and was dropped
    # (Results Section 3.1: 'IIV on LSP in addition to IIV on LS did not
    # result in a significant drop in OFV (dOFV = 0.001)').
    # ---------------------------------------------------------------------

    etalls ~ 0.00471
    # Bosch 2025 Results Section 3.1: 'only the IIV on the RBC lifespan
    # (6.9% CV (18.4% RSE)) ... was estimated'. Supplement S2 NONMEM
    # $OMEGA block: 0.0047104. Note: Bosch 2025 Table 2 reports omega^2 =
    # 0.0145 with 16.5% RSE for OMEGA LS; this is inconsistent with both
    # the supplement S2 NONMEM code (0.0047104) and the paper text
    # (6.9% CV, which back-solves to omega^2 = log(1+0.069^2) = 0.00475).
    # The supplement S2 value 0.0047104 is treated as authoritative; the
    # Table 2 entry 0.0145 is a typographical discrepancy.

    # ---------------------------------------------------------------------
    # Proportional residual error on HbA1c (%) observations.
    # Bosch 2025 Table 2: residual error = 0.0309 with 9.28% RSE.
    # Supplement S2 NONMEM $THETA(5): (0, 0.0311477); $SIGMA: 1 FIX with
    # W = SQRT(THETA(5)^2 * IPRED^2) => W = THETA(5) * IPRED, a pure
    # proportional residual SD of approximately 3.11% of HbA1c.
    # ---------------------------------------------------------------------

    propSd_Hba1c <- 0.0311
    label("Proportional residual error on HbA1c (fraction)")
    # Bosch 2025 supplement S2 $THETA(5) = 0.0311477; Table 2 value 0.0309
    # is the rounded final estimate. Using 0.0311 from the supplement.
    # Multi-output residual-error suffix form propSd_<output>; the model
    # exposes a single HbA1c output but the canonical lint expects the
    # output-suffixed form for any PD-output residual error other than Cc.
  })

  model({
    # ---------------------------------------------------------------------
    # Time-varying regressor input. GLU is the plasma glucose driver,
    # supplied in mmol/L at every observation / event row in the dataset.
    # rxode2 linearly interpolates GLU between adjacent rows.
    # ---------------------------------------------------------------------
    linear(GLU)

    # ---------------------------------------------------------------------
    # Unit conversion. The Lledo-Garcia 2013 / Bosch 2025 supplement S2
    # IGRH parameterisation uses mg/dL for glucose (KG in dL/mg/day,
    # reference glucose 149 mg/dL). We accept GLU and FPG in mmol/L
    # (matching the 4GI module conventions in supplement Tables S3 and
    # S4) and convert here: 1 mmol/L glucose = 18.016 mg/dL (molecular
    # weight 180.16 g/mol, 1 dL/L = 0.1, => 18.016 mg/dL per mmol/L).
    # ---------------------------------------------------------------------
    mmoll_to_mgdl <- 18.016
    cag           <- GLU * mmoll_to_mgdl           # mg/dL, dynamic
    ag            <- FPG * mmoll_to_mgdl           # mg/dL, baseline
    glu_ref       <- 149                           # mg/dL, IGRH reference

    # ---------------------------------------------------------------------
    # IGRH model parameters (back-transformed). nc = 12 transit
    # compartments per Bosch 2025 supplement S2 model code S2 ('NC = 12').
    # The dynamic transit rate KTR2 depends on the current LS, which in
    # turn depends on the current glucose via the (CAG / 149)^gamma_ls
    # power. PREC is the fraction of new RBCs that escape glycation
    # during their precursor lifetime.
    # ---------------------------------------------------------------------
    nc       <- 12
    kin      <- 1                                  # normalised RBC production rate
    kg       <- exp(lkg)                           # dL/mg/day
    lsp      <- exp(llsp)                          # days
    lsrbc    <- exp(lls + etalls)                  # days, per-subject LS
    agls     <- (cag / glu_ref)^gamma_ls           # dynamic LS multiplier
    ls       <- lsrbc * agls                       # days, current LS
    ktr      <- nc / ls                            # 1/day
    prec     <- exp(-kg * cag * lsp)               # fraction surviving precursor

    # ---------------------------------------------------------------------
    # Steady-state initial conditions computed from the baseline glucose
    # AG (mg/dL) following Bosch 2025 supplement S2 model code A_0
    # expressions. Compute LS and KTR at baseline (AG-based), the
    # precursor survival fraction, and propagate down both transit
    # chains analytically.
    # ---------------------------------------------------------------------
    agls0    <- (ag / glu_ref)^gamma_ls
    ls0      <- lsrbc * agls0
    ktr0     <- nc / ls0
    prec0    <- exp(-kg * ag * lsp)
    denom0   <- ktr0 + kg * ag                     # = KTR + KG*AG at t = 0
    rho0     <- ktr0 / denom0                      # geometric ratio along the non-glyc chain

    # Non-glycated chain: hba_n01 receives KIN*PREC0 and outflows at rate
    # ktr0 + KG*AG; each subsequent age cohort receives the inflow at the
    # outflow rate of the previous, giving the geometric series below.
    ic_n01   <- prec0 * kin / denom0
    ic_n02   <- ic_n01 * rho0
    ic_n03   <- ic_n02 * rho0
    ic_n04   <- ic_n03 * rho0
    ic_n05   <- ic_n04 * rho0
    ic_n06   <- ic_n05 * rho0
    ic_n07   <- ic_n06 * rho0
    ic_n08   <- ic_n07 * rho0
    ic_n09   <- ic_n08 * rho0
    ic_n10   <- ic_n09 * rho0
    ic_n11   <- ic_n10 * rho0
    ic_n12   <- ic_n11 * rho0

    # Glycated chain: hba_g01 receives KG*AG from hba_n01 plus the
    # precursors that were already glycated at birth (KIN*(1-PREC0));
    # each subsequent cohort receives the corresponding non-glycated
    # cohort's KG*AG flux plus the upstream glycated outflow. Bosch
    # 2025 supplement S2 model code A_0(14..25) expressions.
    ic_g01   <- (ic_n01 * kg * ag + kin * (1 - prec0)) / ktr0
    ic_g02   <- (ic_g01 * ktr0 + ic_n02 * kg * ag) / ktr0
    ic_g03   <- (ic_g02 * ktr0 + ic_n03 * kg * ag) / ktr0
    ic_g04   <- (ic_g03 * ktr0 + ic_n04 * kg * ag) / ktr0
    ic_g05   <- (ic_g04 * ktr0 + ic_n05 * kg * ag) / ktr0
    ic_g06   <- (ic_g05 * ktr0 + ic_n06 * kg * ag) / ktr0
    ic_g07   <- (ic_g06 * ktr0 + ic_n07 * kg * ag) / ktr0
    ic_g08   <- (ic_g07 * ktr0 + ic_n08 * kg * ag) / ktr0
    ic_g09   <- (ic_g08 * ktr0 + ic_n09 * kg * ag) / ktr0
    ic_g10   <- (ic_g09 * ktr0 + ic_n10 * kg * ag) / ktr0
    ic_g11   <- (ic_g10 * ktr0 + ic_n11 * kg * ag) / ktr0
    ic_g12   <- (ic_g11 * ktr0 + ic_n12 * kg * ag) / ktr0

    # Initial-condition assignment.
    hba_n01(0) <- ic_n01
    hba_n02(0) <- ic_n02
    hba_n03(0) <- ic_n03
    hba_n04(0) <- ic_n04
    hba_n05(0) <- ic_n05
    hba_n06(0) <- ic_n06
    hba_n07(0) <- ic_n07
    hba_n08(0) <- ic_n08
    hba_n09(0) <- ic_n09
    hba_n10(0) <- ic_n10
    hba_n11(0) <- ic_n11
    hba_n12(0) <- ic_n12

    hba_g01(0) <- ic_g01
    hba_g02(0) <- ic_g02
    hba_g03(0) <- ic_g03
    hba_g04(0) <- ic_g04
    hba_g05(0) <- ic_g05
    hba_g06(0) <- ic_g06
    hba_g07(0) <- ic_g07
    hba_g08(0) <- ic_g08
    hba_g09(0) <- ic_g09
    hba_g10(0) <- ic_g10
    hba_g11(0) <- ic_g11
    hba_g12(0) <- ic_g12

    # ---------------------------------------------------------------------
    # Dynamic glycation flux (mg-glucose / dL / day) and dynamic RBC
    # aging rate (1/day). Bosch 2025 supplement S2 model code S2 $DES
    # block: 'PREC2 = EXP(-KG*CAG*LSP); KTR2 = NC/LS2'.
    # ---------------------------------------------------------------------
    prec_dyn  <- exp(-kg * cag * lsp)
    glyc_flux <- kg * cag                          # 1/day per cohort

    # ---------------------------------------------------------------------
    # Non-glycated transit chain. The youngest age cohort (hba_n01)
    # receives newly-produced RBCs at rate KIN * PREC_DYN; each subsequent
    # cohort receives the upstream cohort's outflow at rate KTR. Each
    # cohort also loses RBCs to glycation at rate KG * CAG (per cohort).
    # ---------------------------------------------------------------------
    d/dt(hba_n01) <- kin * prec_dyn       - (ktr + glyc_flux) * hba_n01
    d/dt(hba_n02) <- ktr * hba_n01        - (ktr + glyc_flux) * hba_n02
    d/dt(hba_n03) <- ktr * hba_n02        - (ktr + glyc_flux) * hba_n03
    d/dt(hba_n04) <- ktr * hba_n03        - (ktr + glyc_flux) * hba_n04
    d/dt(hba_n05) <- ktr * hba_n04        - (ktr + glyc_flux) * hba_n05
    d/dt(hba_n06) <- ktr * hba_n05        - (ktr + glyc_flux) * hba_n06
    d/dt(hba_n07) <- ktr * hba_n06        - (ktr + glyc_flux) * hba_n07
    d/dt(hba_n08) <- ktr * hba_n07        - (ktr + glyc_flux) * hba_n08
    d/dt(hba_n09) <- ktr * hba_n08        - (ktr + glyc_flux) * hba_n09
    d/dt(hba_n10) <- ktr * hba_n09        - (ktr + glyc_flux) * hba_n10
    d/dt(hba_n11) <- ktr * hba_n10        - (ktr + glyc_flux) * hba_n11
    d/dt(hba_n12) <- ktr * hba_n11        - (ktr + glyc_flux) * hba_n12

    # ---------------------------------------------------------------------
    # Glycated transit chain. The youngest glycated cohort receives the
    # non-glycated youngest cohort's glycation flux plus the small fraction
    # of new RBCs born already glycated (KIN * (1 - PREC_DYN)). Each
    # subsequent cohort receives the corresponding non-glycated cohort's
    # glycation flux plus the upstream glycated cohort's outflow at rate
    # KTR. Bosch 2025 supplement S2 model code DADT(14..25).
    # ---------------------------------------------------------------------
    d/dt(hba_g01) <- glyc_flux * hba_n01 + kin * (1 - prec_dyn) + ktr * (0       - hba_g01)
    d/dt(hba_g02) <- glyc_flux * hba_n02                          + ktr * (hba_g01 - hba_g02)
    d/dt(hba_g03) <- glyc_flux * hba_n03                          + ktr * (hba_g02 - hba_g03)
    d/dt(hba_g04) <- glyc_flux * hba_n04                          + ktr * (hba_g03 - hba_g04)
    d/dt(hba_g05) <- glyc_flux * hba_n05                          + ktr * (hba_g04 - hba_g05)
    d/dt(hba_g06) <- glyc_flux * hba_n06                          + ktr * (hba_g05 - hba_g06)
    d/dt(hba_g07) <- glyc_flux * hba_n07                          + ktr * (hba_g06 - hba_g07)
    d/dt(hba_g08) <- glyc_flux * hba_n08                          + ktr * (hba_g07 - hba_g08)
    d/dt(hba_g09) <- glyc_flux * hba_n09                          + ktr * (hba_g08 - hba_g09)
    d/dt(hba_g10) <- glyc_flux * hba_n10                          + ktr * (hba_g09 - hba_g10)
    d/dt(hba_g11) <- glyc_flux * hba_n11                          + ktr * (hba_g10 - hba_g11)
    d/dt(hba_g12) <- glyc_flux * hba_n12                          + ktr * (hba_g11 - hba_g12)

    # ---------------------------------------------------------------------
    # Total non-glycated and glycated RBC pools (sums across age cohorts);
    # HbA1c is the glycated fraction expressed as a percent. Bosch 2025
    # supplement S2 model code $ERROR: NON = A(2)+...+A(13); GLY =
    # A(14)+...+A(25); TOT = GLY+NON; HBA = GLY/TOT*100.
    # ---------------------------------------------------------------------
    rbc_non <- hba_n01 + hba_n02 + hba_n03 + hba_n04 + hba_n05 + hba_n06 +
               hba_n07 + hba_n08 + hba_n09 + hba_n10 + hba_n11 + hba_n12
    rbc_gly <- hba_g01 + hba_g02 + hba_g03 + hba_g04 + hba_g05 + hba_g06 +
               hba_g07 + hba_g08 + hba_g09 + hba_g10 + hba_g11 + hba_g12
    rbc_tot <- rbc_non + rbc_gly
    Hba1c   <- rbc_gly / rbc_tot * 100

    Hba1c ~ prop(propSd_Hba1c)
  })
}
