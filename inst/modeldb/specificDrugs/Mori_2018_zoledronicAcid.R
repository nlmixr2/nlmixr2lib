Mori_2018_zoledronicAcid <- function() {
  description <- "Kinetic-pharmacodynamic (K-PD) PK / bone-turnover-marker / lumbar-spine BMD model for once-yearly intravenous zoledronic acid (ZOL) 5 mg in Japanese patients with primary osteoporosis (ZONE study). A virtual effect-site amount A receives the administered dose and decays first-order at rate KD; the drug-effect signal KD*A enters a sigmoidal Imax factor with Hill coefficient Gamma and half-effect EKD50 that inhibits the zero-order synthesis Kin of the serum bone-resorption marker (tartrate-resistant acid phosphatase 5b, TRACP-5b), which is eliminated first-order at Kout. The observed marker carries a multiplicative disease-progression / supplementation drift (1 + Slope * t + Emax * t / (T50 + t)) capturing the daily oral calcium + vitamin D + magnesium supplementation effect pooled with natural osteoporosis progression (the two effects could not be separated because all subjects received the supplements). Lumbar-spine BMD follows an effect-compartment ODE with rate Ke0 whose target is BMD0 + Scale * (marker - Marker0), where Scale (negative) is the marker-to-BMD coupling. Baseline TRACP-5b (TRACP5B_BL) enters EKD50, Slope, T50, and (active-arm-only) Scale as a power-model covariate centred on the cohort reference 400 mU / dL."
  reference <- paste(
    "Mori Y, Kasai H, Ose A, Serada M, Ishiguro M, Shiraki M, Tanigawara Y (2018).",
    "Modeling and simulation of bone mineral density in Japanese osteoporosis patients",
    "treated with zoledronic acid using tartrate-resistant acid phosphatase 5b, a bone",
    "resorption marker.",
    "Osteoporos Int 29(5):1155-1163.",
    "doi:10.1007/s00198-018-4376-1.",
    sep = " "
  )
  vignette <- "Mori_2018_zoledronicAcid"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "mU/dL"
    # units$concentration documents the serum TRACP-5b biomarker output; this is
    # a K-PD model with no observed plasma drug concentration. Outputs are:
    # TRACP5b (mU / dL, the bone-resorption marker) and BMD (g / cm^2, lumbar
    # spine L2-L4 areal density). The dose unit (mg) is the IV zoledronic-acid
    # 5 mg annual infusion; the effect-site amount A (depot_kpd) carries mg
    # and the KD * A signal enters the inhibition factor in mg / day units
    # which match the published EKD50 of 5.776e-3 mg / day.
  )

  paper_specific_compartments <- c("bmd")

  covariateData <- list(
    ON_TREATMENT = list(
      description        = "Treatment-arm indicator: 1 = subject is in the active zoledronic-acid arm (received yearly IV ZOL 5 mg); 0 = subject is in the placebo arm. Per-subject and time-fixed.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm)",
      notes              = paste(
        "Gates the multiplicative TRACP5B_BL power-effect on Scale (the marker-to-BMD",
        "coupling): the Scale covariate effect is significant only in the ZOL arm",
        "(Mori 2018 Results, Covariate exploration). When ON_TREATMENT = 0 the",
        "Scale-side exponent collapses to 0 (power factor = 1) and the placebo-arm",
        "typical Scale value is recovered exactly; when ON_TREATMENT = 1 the",
        "TRACP5B_BL-centred power exponent -1.112 enters the Scale calculation.",
        "ON_TREATMENT does NOT gate dose administration -- the active-arm dose enters",
        "via dose events on the depot_kpd compartment (5 mg IV at t = 0 and t = 365",
        "days); placebo subjects simply receive no dose events."
      ),
      source_name        = "Treatment group (Mori 2018 Table 1)"
    ),
    TRACP5B_BL = list(
      description        = "Baseline (pre-dose) serum tartrate-resistant acid phosphatase 5b (TRACP-5b) concentration measured by the Osteolinks fragment-absorbed immunocapture enzymatic assay. Time-fixed per subject. Drives (a) the steady-state TRACP-5b pool through Kin = TRACP5B_BL * Kout and the marker initial condition `effect(0) <- TRACP5B_BL`, and (b) a power-model covariate effect on EKD50, Slope, T50 (both arms) and Scale (active-arm only, gated by ON_TREATMENT) centred on the cohort reference 400 mU / dL.",
      units              = "mU/dL",
      type               = "continuous",
      reference_category = "n/a -- power-model standardisation to 400 mU/dL (cohort mean baseline = 401.1 mU/dL; the paper states standardisation to the cohort median but the numeric median is not published)",
      notes              = paste(
        "Mori 2018 Table 1 reports baseline TRACP-5b mean 401.1 +/- 147.9 mU / dL,",
        "range [157, 1240] mU/dL across N = 306 patients. The cohort median is not",
        "published; this file uses 400 mU / dL as the documented standardisation",
        "reference (rounded cohort mean) and notes the assumption in the vignette",
        "Assumptions and deviations. Distinct from a time-varying TRACP-5b output",
        "(which is the modelled state in `effect`, with its observation channel",
        "TRACP5b)."
      ),
      source_name        = "Baseline TRACP-5b (Mori 2018 Table 1, Table 2 covariate rows)"
    ),
    BMD_BL = list(
      description        = "Baseline (pre-dose) lumbar-spine (L2-L4) bone mineral density measured by dual-energy X-ray absorptiometry (DXA, Hologic). Time-fixed per subject. Used as the BMD compartment initial condition `bmd(0) <- BMD_BL` and as the deviation reference inside the BMD effect-compartment ODE.",
      units              = "g/cm^2",
      type               = "continuous",
      reference_category = "n/a -- used as a per-subject anchor; no power-form covariate effect.",
      notes              = paste(
        "Mori 2018 Table 1 reports baseline lumbar-spine BMD mean 0.677 +/- 0.094",
        "g / cm^2, range [0.36, 0.98] across N = 306 patients. The standard DXA",
        "areal-density unit g / cm^2 is preserved (the corresponding T-score range",
        "of -5.51 to -0.25 is reported in Table 1 but not used in the model)."
      ),
      source_name        = "Baseline lumbar BMD (Mori 2018 Table 1)"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male). Mori 2018 reports 94.4% female (289 of 306 patients) and screened sex as a candidate covariate but did not retain it in the final model.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Mori 2018 Methods 'Statistical models' lists sex among the screened patient characteristics; Results 'Covariate exploration' confirms only baseline TRACP-5b survived the forward / backward selection on EKD50, Slope, T50, and Scale. Not retained in the final model."
    ),
    AGE = list(
      description = "Subject age in years at study entry. Mori 2018 cohort mean 72.9 +/- 5.2 years, range [65, 87].",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate via the power-model parameterisation; not retained in the final model."
    ),
    WT = list(
      description = "Body weight in kg at study entry. Mori 2018 cohort mean 52.3 +/- 8.0 kg, range [34.1, 83.6].",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as a continuous covariate via the power-model parameterisation; not retained in the final model."
    ),
    PRIOR_BIO = list(
      description = "Prior bisphosphonate-treatment indicator: 1 = previously treated with a bisphosphonate (with sufficient washout per the ZONE study protocol), 0 = never treated. Mori 2018 cohort 9.5% had prior bisphosphonate use with sufficient washout.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a categorical covariate (modelled in a relative-effect manner per Methods 'Statistical models'); not retained in the final model. PRIOR_BIO is the closest existing canonical (biologic encompasses many therapies); the source paper screens specifically for prior bisphosphonate use. Not used in the final model so the close-but-not-exact canonical is acceptable in covariatesDataExcluded (documentation only)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 306L,
    n_studies      = 1L,
    age_range      = "65-87 years (mean 72.9 +/- 5.2)",
    weight_range   = "34.1-83.6 kg (mean 52.3 +/- 8.0)",
    sex_female_pct = 94.4,
    race_ethnicity = c(Japanese = 100),
    disease_state  = "Primary osteoporosis (Japanese ZONE study; lumbar-spine T-score mean -2.80 +/- 0.80, range [-5.51, -0.25]).",
    dose_range     = "Zoledronic acid 5 mg IV infused over 15 minutes once yearly; placebo arm received no drug. All subjects received daily oral calcium 610 mg + vitamin D 400 IU + magnesium 30 mg supplementation throughout the 2-year study.",
    regions        = "Japan",
    n_observations = c(TRACP5b = 3410L, BMD = 1146L),
    bl_tracp5b_mU_per_dL = "Mean 401.1 +/- 147.9 mU/dL (range [157, 1240])",
    bl_bmd_g_per_cm2     = "Mean 0.677 +/- 0.094 g/cm^2 (range [0.36, 0.98])",
    notes          = paste(
      "ZONE study (Mori 2018 ref [10]) was a 2-year multicentre randomised",
      "double-blind placebo-controlled parallel-group trial; 665 patients were",
      "randomised (333 ZOL, 332 placebo) and 306 patients (145 ZOL, 161 placebo)",
      "who had at least one bone-resorption marker value AND one BMD value were",
      "used in this analysis. Bone-resorption markers TRACP-5b, CTx, and u-NTx",
      "(3410 marker points) and lumbar-spine L2-L4 BMD (1146 BMD points) were",
      "measured under overnight-fasted conditions. The final model uses TRACP-5b",
      "as the marker (selected by statistical significance from the three; CTx",
      "and u-NTx alternative-marker fits are reported in the supplementary tables",
      "but are not extracted here)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters (Mori 2018 Table 2, TRACP-5b final model).
    # Point estimates source-traced to Table 2 rows. SE (Bootstrap, n=200)
    # is documented per-line for transparency but not used in encoding
    # because no IIV / fixed flag is implied by the SE alone.
    # ---------------------------------------------------------------

    # K-PD virtual effect-site kinetics
    lkel <- log(3.719e-3); label("Effect-site first-order disappearance rate KD (1/day)")                # Mori 2018 Table 2: KD = 3.719e-3 /day (SE 2.941e-5)

    # Sigmoidal inhibition of marker synthesis by KD * A
    lhill  <- log(4.583e-1); label("Hill coefficient Gamma (unitless)")                                  # Mori 2018 Table 2: Gamma = 4.583e-1 (SE 2.626e-3)
    lekd50 <- log(5.776e-3); label("Drug-effect signal at half-maximum inhibition EKD50 (mg/day; KD * A units)")  # Mori 2018 Table 2: EKD50 = 5.776e-3 mg/day (SE 3.800e-5)

    # Marker indirect-response turnover (Kin computed from baseline; Kout estimated)
    lkout <- log(4.584e-1); label("TRACP-5b first-order elimination rate Kout (1/day)")                  # Mori 2018 Table 2: Kout = 4.584e-1 /day (SE 3.592e-3)

    # Disease-progression / supplementation drift on observed marker
    # Slope and T50 are positive => log-normal; Emax can be negative => normal.
    lslope_dp <- log(2.688e-4); label("Linear time-effect slope Slope on observed TRACP-5b (1/day)")     # Mori 2018 Table 2: Slope = 2.688e-4 /day (SE 2.877e-6)
    emax_dp   <- -8.266e-2;    label("Relative maximal disease-progression / supplementation effect Emax on observed TRACP-5b (unitless)")  # Mori 2018 Table 2: Emax = -8.266e-2 (SE 2.053e-3)
    lt50      <- log(1.166e2); label("Time at which the disease-progression effect reaches half of Emax T50 (day)")  # Mori 2018 Table 2: T50 = 1.166e2 day (SE 3.091)

    # BMD effect-compartment kinetics
    lke0      <- log(3.802e-3);  label("BMD effect-compartment equilibrium rate Ke0 (1/day)")            # Mori 2018 Table 2: Ke0 = 3.802e-3 /day (SE 6.429e-4)
    scale_bmd <- -2.521e-4;       label("Marker-to-BMD coupling Scale ((g/cm^2)/(mU/dL))")                # Mori 2018 Table 2: Scale = -2.521e-4 (g/cm^2)/(mU/dL) (SE 3.709e-5)

    # ---------------------------------------------------------------
    # Baseline-TRACP-5b power covariate effects (Mori 2018 Table 2).
    # Centred on TRACP5B_BL_REF = 400 mU/dL (cohort reference; see Errata).
    # Scale covariate is significant ONLY in the active arm (ON_TREATMENT = 1).
    # ---------------------------------------------------------------
    e_tracp5b_bl_ekd50 <- -1.534; label("Power exponent for TRACP5B_BL on EKD50 (unitless)")  # Mori 2018 Table 2: TRACP-5b baseline effect on EKD50 = -1.534 (SE 1.699e-2)
    e_tracp5b_bl_slope <- -1.350; label("Power exponent for TRACP5B_BL on Slope (unitless)")  # Mori 2018 Table 2: TRACP-5b baseline effect on Slope = -1.350 (SE 2.093e-2)
    e_tracp5b_bl_t50   <- -1.319; label("Power exponent for TRACP5B_BL on T50 (unitless)")    # Mori 2018 Table 2: TRACP-5b baseline effect on T50 = -1.319 (SE 2.974e-2)
    e_tracp5b_bl_scale <- -1.112; label("Power exponent for TRACP5B_BL on Scale (active arm only, gated by ON_TREATMENT; unitless)")  # Mori 2018 Table 2: TRACP-5b baseline effect on Scale (ZOL arm only) = -1.112 (SE 2.635e-1)

    # ---------------------------------------------------------------
    # Inter-individual variability (Mori 2018 Table 2, variance scale).
    # Positive-valued parameters (EKD50, Slope, T50, Ke0) carry log-normal
    # etas; Emax (signed) and Scale (negative) carry additive normal etas.
    # The Emax / Slope / T50 etas are reported as a 3 x 3 block matrix
    # with the off-diagonals shown in Table 2; EKD50, Ke0 and Scale etas
    # are reported as independent diagonal entries.
    # ---------------------------------------------------------------
    etalekd50 ~ 2.255e-1                                                  # Mori 2018 Table 2: omega^2(EKD50) = 2.255e-1 (SE 4.892e-3)

    # Emax + Slope + T50 block (Mori 2018 Table 2 IIV block):
    #   row 1: omega^2(Emax)   = 8.584e-2 (SE 2.655e-3)
    #   row 2: omega(Emax,Slope) = -7.576e-2 (SE 2.806e-3),
    #          omega^2(Slope) = 1.573e-1 (SE 3.136e-3)
    #   row 3: omega(Emax,T50) = -8.732e-2 (SE 3.911e-3),
    #          omega(Slope,T50) = 8.697e-2 (SE 4.525e-3),
    #          omega^2(T50) = 2.095 (SE 3.662e-2)
    # Variances on the diagonal, off-diagonals from Table 2. The
    # published covariance entries are interpreted on the mixed eta
    # scale (log for Slope / T50, linear for Emax) exactly as reported
    # in the source.
    etaemax_dp + etalslope_dp + etalt50 ~ c(
      8.584e-2,
      -7.576e-2, 1.573e-1,
      -8.732e-2, 8.697e-2, 2.095
    )

    etalke0     ~ 1.406e-3                                                # Mori 2018 Table 2: omega^2(Ke0)   = 1.406e-3 (SE 2.964e-1; high RSE flagged in vignette Errata)
    etascale_bmd ~ 2.361e-8                                               # Mori 2018 Table 2: omega^2(Scale) = 2.361e-8 (SE 1.028e-8; very small, retained per the source)

    # ---------------------------------------------------------------
    # Residual error (Mori 2018 Table 2).
    # Marker: text says "relative error model" but the printed sigma
    # value 35.51 mU / dL is consistent with an additive error on the
    # raw mU / dL scale (35.51 / 401.1 mean baseline ~ 8.85%, matching
    # the 8.9% intra-individual variability the Discussion attributes
    # to the TRACP-5b model). Encoded as additive on the marker scale;
    # the text-vs-value discrepancy is documented in vignette Errata.
    # BMD: text confirms "BMD was modeled with absolute error" =>
    # additive on g / cm^2.
    # ---------------------------------------------------------------
    addSd_TRACP5b <- 35.51;   label("Additive residual SD on TRACP-5b (mU/dL); paper text says 'relative' but the numeric value matches additive (35.51 / 401 = 8.85% ~ 8.9% Discussion CV)")  # Mori 2018 Table 2: sigma (TRACP-5b) = 3.551 x 10 = 35.51 mU/dL (SE 1.204e-1)
    addSd_BMD     <- 2.447e-2; label("Additive residual SD on BMD (g/cm^2)")                          # Mori 2018 Table 2: sigma (BMD) = 2.447e-2 g/cm^2 (SE 1.215e-3)
  })

  model({
    # ---------------------------------------------------------------
    # Centring reference for baseline-TRACP-5b power covariate.
    # The paper says continuous covariates are standardised to their
    # cohort median; the numeric median is not published, so the
    # rounded cohort mean (400 mU/dL; mean 401.1 mU/dL) is used (see
    # covariateData[[TRACP5B_BL]]$notes and vignette Errata).
    # ---------------------------------------------------------------
    TRACP5B_BL_REF <- 400

    # ---------------------------------------------------------------
    # Individual parameters.
    # Positive parameters: typical-value-exponential form
    #   p_i = exp(lp + etalp + cov_effect)
    # Signed parameters (Emax, Scale): typical + additive eta + power covariate
    # The Scale covariate effect is gated by ON_TREATMENT so the placebo
    # arm recovers the typical Scale value exactly (exponent collapses
    # to 0 when ON_TREATMENT = 0).
    # ---------------------------------------------------------------
    log_cov_ratio <- log(TRACP5B_BL / TRACP5B_BL_REF)

    kd    <- exp(lkel)
    hill  <- exp(lhill)
    ekd50 <- exp(lekd50 + etalekd50 + e_tracp5b_bl_ekd50 * log_cov_ratio)
    kout  <- exp(lkout)
    slope_dp  <- exp(lslope_dp + etalslope_dp + e_tracp5b_bl_slope * log_cov_ratio)
    emax_dp_i <- emax_dp + etaemax_dp
    t50       <- exp(lt50 + etalt50 + e_tracp5b_bl_t50 * log_cov_ratio)
    ke0       <- exp(lke0 + etalke0)
    scale_bmd_i <- (scale_bmd + etascale_bmd) *
                   exp(ON_TREATMENT * e_tracp5b_bl_scale * log_cov_ratio)

    # ---------------------------------------------------------------
    # Indirect-response steady-state: Kin = Marker0 * Kout so that the
    # no-drug baseline marker equals TRACP5B_BL (the per-subject
    # baseline TRACP-5b covariate). Marker0 enters as TRACP5B_BL.
    # ---------------------------------------------------------------
    marker0 <- TRACP5B_BL
    kin     <- marker0 * kout

    # ---------------------------------------------------------------
    # K-PD virtual effect-site amount A: dose lands in depot_kpd (mg)
    # and decays first-order at rate KD. There is no measured plasma
    # drug concentration; the drug-effect signal is KD * A (mg / day),
    # which has the same units as EKD50.
    # ---------------------------------------------------------------
    d/dt(depot_kpd) <- -kd * depot_kpd
    drug_signal     <-  kd * depot_kpd

    # ---------------------------------------------------------------
    # Sigmoidal inhibition factor INH = 1 - signal^hill / (EKD50^hill + signal^hill)
    # Equal to 1 when signal = 0 (no drug, full Kin); approaches 0 as
    # signal grows. Multiplies the zero-order synthesis term.
    # ---------------------------------------------------------------
    inh_factor <- 1 - drug_signal^hill / (ekd50^hill + drug_signal^hill)

    # ---------------------------------------------------------------
    # Bone-resorption marker (TRACP-5b) compartment.
    # State variable lives in `effect` (the indirect-response pool);
    # initial condition is the per-subject baseline TRACP5B_BL.
    # ---------------------------------------------------------------
    d/dt(effect) <- kin * inh_factor - kout * effect
    effect(0)   <- marker0

    # ---------------------------------------------------------------
    # Disease-progression / supplementation drift multiplier applied to
    # the observed marker (Mori 2018 Methods, marker-model equation):
    #   Marker(t) = effect * (1 + Slope * t + Emax * t / (T50 + t))
    # At t = 0 the multiplier is 1 so the observation matches the
    # state baseline. Drives the TRACP5b observation only -- the BMD
    # ODE in the published model uses the marker-state deviation
    # (effect - Marker0), not the DP-adjusted observation (see
    # vignette Errata: this is the literal reading of the printed
    # Marker / Marker(t) notation and the schematic in Figure 1).
    # ---------------------------------------------------------------
    dp_multiplier <- 1 + slope_dp * t + emax_dp_i * t / (t50 + t)
    TRACP5b       <- effect * dp_multiplier

    # ---------------------------------------------------------------
    # Lumbar-spine BMD effect-compartment ODE:
    #   dBMD/dt = Ke0 * [Scale * (marker - Marker0) - (BMD - BMD0)]
    # Initial condition is the per-subject baseline BMD_BL. The Scale
    # value Scale (negative) means a lower marker (drug-induced
    # decrease) raises the BMD target. Marker0 here is TRACP5B_BL.
    # ---------------------------------------------------------------
    d/dt(bmd) <- ke0 * (scale_bmd_i * (effect - marker0) - (bmd - BMD_BL))
    bmd(0)    <- BMD_BL

    BMD <- bmd

    # ---------------------------------------------------------------
    # Observation residual error (Mori 2018 Table 2). Multi-output:
    # the TRACP-5b channel uses additive error on the raw marker scale
    # (paper text says "relative error" but the numeric sigma value
    # matches additive; see vignette Errata for the encoding
    # rationale), and the BMD channel uses additive error on g / cm^2
    # (text confirms "absolute error" for BMD).
    #
    # Event-table observation rows must use the observable name on
    # the `cmt` column (`cmt = "TRACP5b"` for marker observations,
    # `cmt = "BMD"` for BMD observations). rxode2 / rxUi auto-allocates
    # virtual compartment slots for the algebraic observables AFTER
    # the ODE-state slots; the algebraic observables (and the ODE
    # states themselves) are returned as columns in the rxSolve output
    # dataframe regardless of which cmt the obs row targets. The model
    # has no ODE states past the auto-injected virtual slots, so no
    # state references are broken by the auto-injection.
    # ---------------------------------------------------------------
    TRACP5b ~ add(addSd_TRACP5b)
    BMD     ~ add(addSd_BMD)
  })
}
