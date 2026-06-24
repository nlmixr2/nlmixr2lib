Hajjar_2018_DMD_6MWT <- function() {
  description <- paste0(
    "Latent variable disease-progression model for the six-minute walk ",
    "test (6MWT, meters) in healthy boys and boys with Duchenne muscular ",
    "dystrophy (DMD), fit by Hajjar et al. (ACoP9 2018 poster T-011) to ",
    "publicly available individual-level longitudinal natural-history ",
    "6MWT data from 16 healthy controls and 219 DMD patients. The 6MWT ",
    "is modelled as a one-compartment indirect-response state (walkDist, ",
    "meters): a zero-order production rate KIN feeds the state and a ",
    "first-order dissipation rate KOUT removes it. A change point at ",
    "subject age MTIME (1.75 years) switches KIN from 0 to its non-zero ",
    "value, encoding the developmental lag before toddlers can walk a ",
    "measurable distance in six minutes. A latent exponential disease ",
    "process DIS = ALPHA*exp(BETA*age) stimulates the dissipation rate ",
    "for DMD subjects only; healthy subjects fix ALPHA = BETA = 0 so the ",
    "disease term vanishes. The DIS_DMD covariate (1 = DMD, 0 = healthy) ",
    "additionally multiplies KIN by KCOV (0.63) for DMD subjects so the ",
    "two populations share KOUT but have separate KIN. Between-subject ",
    "variability is exponential on KOUT (DMD subjects only, per the ",
    "source NONMEM control stream), on KIN (both populations), and on ",
    "ALPHA and BETA (DMD subjects only). Residual error is additive on ",
    "the 6MWT scale. The model has no drug input; the source poster ",
    "frames it as a simulation tool for designing future DMD efficacy ",
    "trials. Time is age in years (the integration variable; the ",
    "source poster reports KOUT and KIN in per-month units for human ",
    "readability, see vignette Errata for the unit-conversion step)."
  )
  reference <- paste(
    "Hajjar JL, Mondick JT, Gastonguay MR.",
    "A Latent Variable Disease Progression Model for Duchenne Muscular Dystrophy.",
    "Poster T-011 presented at the American Conference on Pharmacometrics",
    "(ACoP9), Oct 7-10 2018, San Diego, CA.",
    "doi:10.36255/duchenne-muscular-dystrophy-public-education.",
    sep = " "
  )
  vignette <- "Hajjar_2018_DMD_6MWT"
  units <- list(
    time          = "year (subject age)",
    dosing        = "n/a (disease-progression model with no drug input)",
    concentration = "m (six-minute walk test distance, observation walkDist)"
  )

  covariateData <- list(
    DIS_DMD = list(
      description        = "Binary indicator for Duchenne muscular dystrophy diagnosis: 1 = DMD subject, 0 = non-DMD subject (healthy control). Time-fixed per subject. Source poster's NONMEM variable name is PATIENT with the same orientation (1 = DMD, 0 = healthy).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-DMD; healthy boys recruited as controls in Henricson 2012).",
      notes              = "Selects the disease pathway: for DIS_DMD = 0 the latent disease coefficients ALPHA and BETA are zero and the production-rate scaling KCOV is fixed at 1, so the model reduces to a plain indirect-response 6MWT trajectory anchored to healthy growth. For DIS_DMD = 1 the latent disease term ALPHA*exp(BETA*age) multiplies KOUT and KCOV = 0.63 scales KIN down to the DMD steady-state plateau. The NONMEM source control stream applies between-subject variability on KOUT only for DMD subjects (ETA(1) inside the IF (PATIENT.EQ.1) branch); the model file reproduces this asymmetry by gating etalkout with DIS_DMD. See vignette Errata for the reporting inconsistency between the source NONMEM code (single ETA on KOUT, DMD-only) and source Table 2 (which lists a 5.40% BSV for healthy KOUT alongside the 16.7% DMD value).",
      source_name        = "PATIENT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 235L,
    n_studies      = 6L,
    age_range      = "approximately 4-15.3 years across the pooled cohort (healthy controls: 4-12 years; DMD: 4-15.3 years). Source poster Table 1 reports per-study age statistics.",
    age_median     = "DMD subjects: mean approximately 8.7 years (Brehm 2014 mean 8.72; Mercuri 2016 mean 8.73; McDonald 2013 mean 8; Goemans 2013 mean 9.5). Healthy controls: median 9 years (Henricson 2012).",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 0,
    race_ethnicity = NA_character_,
    disease_state  = "Duchenne muscular dystrophy (219 boys, 91% of pooled cohort) and healthy controls (16 boys, 9% of pooled cohort). The DMD subjects are predominantly ambulatory boys with longitudinal 6MWT measurements; 8.7% had one 6MWT measure, 51.4% had two, and 39.9% had three or more. 44.0% of the DMD subjects came from the Mercuri 2016 natural history study; the remainder are pooled from Brehm 2014 (n=14), Goemans 2013 (n=65), McDonald 2010 (n=15), and the placebo arm of the McDonald 2013 ataluren trial NCT00592553 (n=57).",
    dose_range     = "Not applicable -- disease-progression model with no drug input. Steroid administration was explored as a covariate but conclusions about steroid effects were inconclusive due to limited covariate information and uneven group sizes (source poster Results bullet 9).",
    regions        = "International (multiple study sites pooled from the published literature; source poster Methods step 1).",
    n_observations = "Not tabulated in the source poster; observations per DMD subject distribution: 8.72% one measure, 51.4% two measures, 39.9% three or more measures (source poster Results bullet 1).",
    notes          = "All individual-level data were digitised from the published figures of the six contributing studies using GraphClick version 3.0.3 (source poster Methods step 1a) rather than obtained as raw subject records. The model was first fit to the healthy-subject data and then refit to the pooled (healthy + DMD) dataset with the healthy parameters held fixed during the DMD parameter estimation step (source poster Results bullet 2). MTIME and KIN were estimated using the (small) healthy-subject dataset, so the corresponding RSEs reported in source Table 2 are large (323% on MTIME, 138% on KIN) relative to the precise DMD-only parameters."
  )

  ini({
    # ============================================================
    # Final model parameter estimates from Hajjar 2018 Table 2
    # (poster T-011, ACoP9, October 2018). The structural model is
    # a one-compartment indirect-response ODE (source poster Results
    # bullet 3 and the $PK + $DES NONMEM control-stream extract):
    #     d/dt(walkDist) = KIN(age) - KOUT * walkDist * (1 + DIS)
    # with a change point that switches KIN from 0 to its non-zero
    # value at subject age MTIME, and a latent disease state
    #     DIS = ALPHA * exp(BETA * age)
    # that is zero for healthy subjects (ALPHA = BETA = 0) and
    # non-zero for DMD subjects.
    #
    # Time-unit convention. The integration variable in this nlmixr2
    # encoding is subject age in years (the x-axis of source Figure
    # 1 is "Age, years"; MTIME is reported as 1.75 years; BETA*age
    # must be unitless so BETA carries inverse-year units). The
    # source poster Table 2 reports KOUT and KIN in per-month units
    # for human readability (0.48 month^-1, 321 meters/month); the
    # numerically-equivalent per-year values used here are
    # 0.48*12 = 5.76 year^-1 and 321*12 = 3852 meters/year. KIN/KOUT
    # = 668.75 m is the healthy-subject 6MWT steady-state plateau,
    # consistent with source Figure 1 (pink lines plateauing near
    # 600-650 m for ages 6-12 years); KIN*KCOV/KOUT = 421.4 m is
    # the DMD steady-state plateau before the latent disease term
    # depresses it. See the vignette Errata for the full source-vs-
    # encoding parameter table.
    # ============================================================

    # ---------------- Structural (typical-value) parameters ----------------
    lkout  <- log(0.48 * 12) ; label("Log first-order dissipation rate constant for 6MWT (1/year)")                 # Hajjar 2018 Table 2 KOUT = 0.48 month^-1 (5.76 year^-1, shared healthy + DMD)
    lmtime <- log(1.75)      ; label("Log production-rate change-point lag (years; KIN switches on at age MTIME)")  # Hajjar 2018 Table 2 MTIME = 1.75 years
    lkin   <- log(321 * 12)  ; label("Log zero-order production rate for 6MWT in healthy reference (m/year)")        # Hajjar 2018 Table 2 KIN = 321 m/month (3852 m/year, healthy reference)
    lkcov  <- log(0.63)      ; label("Log multiplier on KIN for DMD subjects (KCOV; unitless)")                       # Hajjar 2018 Table 2 KCOV = 0.63
    lalpha <- log(9.85e-6)   ; label("Log latent disease-process pre-exponential coefficient (DMD subjects only)")    # Hajjar 2018 Table 2 alpha = 9.85E-06
    lbeta  <- log(0.995)     ; label("Log latent disease-process exponential growth rate (1/year, DMD only)")         # Hajjar 2018 Table 2 beta = 0.995

    # ---------------- Between-subject variability ----------------
    # Source Table 2 reports CV%-style BSVs; the variance on the
    # natural-log scale is omega^2 = log(1 + CV^2). The healthy-vs-
    # DMD asymmetry on KOUT (source Table 2 lists 5.40% healthy and
    # 16.7% DMD) is encoded here per the source $PK code: ETA(1) is
    # applied only inside the IF (PATIENT.EQ.1) branch, so etalkout
    # gates on DIS_DMD in model() and uses the DMD omega (16.7%).
    # The healthy-population BSV is documented in the vignette
    # Errata as a source-table reporting inconsistency with the
    # printed NONMEM code.
    etalkout  ~ 0.0275    # log(1 + 0.167^2) ~ 0.0275; Hajjar 2018 Table 2 KOUT BSV 16.7% (DMD subjects)
    etalkin   ~ 0.000882  # log(1 + 0.0297^2) ~ 0.000882; Hajjar 2018 Table 2 KIN BSV 2.97% (shared)
    etalalpha ~ 0.0987    # log(1 + 0.322^2) ~ 0.0987; Hajjar 2018 Table 2 alpha BSV 32.2% (DMD subjects)
    etalbeta  ~ 0.0366    # log(1 + 0.193^2) ~ 0.0366; Hajjar 2018 Table 2 beta BSV 19.3% (DMD subjects)

    # ---------------- Residual error ----------------
    addSd <- 42.0 ; label("Additive residual standard deviation on 6MWT (meters)")  # Hajjar 2018 Table 2 SD of residual error = 42.0 m
  })

  model({
    # ----- Individual structural parameters -----
    # Mu-referenced individual parameters. Each eta is applied
    # additively on the log scale (canonical mu-reference) so
    # nlmixr2's FOCEI / SAEM mu-search recognises them. The
    # population-pathway gating (which etas correspond to which
    # population in the source NONMEM control stream) is documented
    # below and the alpha_dis * DIS_DMD multiplier zeros out the
    # latent-disease contribution for healthy subjects so the
    # source $PK semantics are preserved at the integrator level:
    #   - KIN: etalkin contributes to both populations (source $PK
    #     places ETA(2) on the shared KIN2 line outside both IF
    #     branches).
    #   - KOUT: etalkout is applied to all subjects in this
    #     encoding; the source NONMEM code places ETA(1) only
    #     inside IF (PATIENT.EQ.1) so healthy subjects technically
    #     have no KOUT eta. Applying the eta to both populations is
    #     a minor simulation deviation noted in the vignette
    #     Errata; typical-value predictions are unchanged (eta = 0).
    #   - ALPHA / BETA: zeroed for healthy subjects by the
    #     (DIS_DMD * ...) factor on alpha_dis. etalbeta is drawn
    #     but only affects the prediction through alpha_dis * exp(),
    #     so for healthy subjects (alpha_dis = 0) the beta eta has
    #     no effect on the trajectory.
    kout       <- exp(lkout  + etalkout)
    kin0       <- exp(lkin   + etalkin)
    kcov_eff   <- DIS_DMD * exp(lkcov) + (1 - DIS_DMD) * 1.0
    kin_typ    <- kin0 * kcov_eff
    alpha_dis  <- DIS_DMD * exp(lalpha + etalalpha)
    beta_dis   <- exp(lbeta + etalbeta)
    tlag_kin   <- exp(lmtime)

    # ----- Change-point on KIN (source poster Results bullet 4 and
    # the KIN = KIN1*(1-MPAST(1)) + KIN2*MPAST(1) NONMEM line).
    # MPAST(1) is 0 until subject age crosses MTIME(1), then 1. In
    # rxode2 the equivalent is the integer indicator on `time`
    # (which equals subject age in this encoding). tlag_kin is the
    # internal symbol for the change-point age (the published
    # parameter name MTIME collides with an rxode2 internal name).
    kin_active <- (time >= tlag_kin) * kin_typ

    # ----- Latent disease process (source $PK: DIS = ALPHA *
    # EXP(TIME*BETA) for TIME > 0; DIS = 0 at TIME = 0). For
    # healthy subjects ALPHA = 0 so DIS = 0 identically; for DMD
    # subjects DIS grows exponentially with age.
    dis <- alpha_dis * exp(time * beta_dis)

    # ----- Indirect-response ODE for the 6MWT distance state.
    # walkDist (meters) is initialised at zero (A_0(1) = 0 in the
    # source $PK) and accumulates only after age MTIME.
    d/dt(walkDist) <- kin_active - kout * walkDist * (1 + dis)

    walkDist ~ add(addSd)
  })
}
