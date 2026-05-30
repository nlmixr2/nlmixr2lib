Xiang_2018_baicalein <- function() {
  description <- paste(
    "In vitro (RAW264.7 mouse macrophage cell line).",
    "Semi-mechanism-based cellular pharmacodynamic model for the",
    "anti-inflammatory effect of baicalein on LPS-induced cytokine and",
    "iNOS/NO release in RAW264.7 mouse macrophages. Four indirect-response",
    "states arranged in a TNF-alpha -> {IL-6, iNOS -> NO} cascade:",
    "(1) TNF-alpha indirect response with LPS-stimulated zero-order",
    "production and baicalein's log-linear inhibition f(Bai) =",
    "alpha * log(C_Bai + 1) on the production rate;",
    "(2) IL-6 indirect response with delayed TNF-alpha drive (lag tau1);",
    "(3) iNOS indirect response with delayed TNF-alpha drive (lag tau2)",
    "and elimination held at zero (paper choice for the post-12.5 h plateau);",
    "(4) NO indirect response with iNOS^delta amplification and",
    "elimination held at zero. Baicalein concentration enters as a static",
    "covariate (CONC_BAI_UM); LPS is constant at 1 ug/mL throughout the",
    "experiment and is absorbed into the kinTNF zero-order production",
    "rate (no explicit LPS state). Tau1 and tau2 are encoded via",
    "single-compartment delay states (mean transit time = tau) because",
    "rxode2 does not provide a native delay-differential-equation solver;",
    "see the validation vignette Errata for the impulse-response",
    "implications. Typical-value-only mechanism: no IIV, no residual",
    "error is reported in the source (CV% in Table 1 are RSE / precision",
    "of estimate, not BSV)."
  )

  reference <- paste(
    "Xiang L, Hu Y-F, Wu J-S, Wang L, Huang W-G, Xu C-S, Meng X-L, Wang P.",
    "Semi-Mechanism-Based Pharmacodynamic Model for the Anti-Inflammatory",
    "Effect of Baicalein in LPS-Stimulated RAW264.7 Macrophages.",
    "Front Pharmacol. 2018;9:793. doi:10.3389/fphar.2018.00793."
  )

  vignette <- "Xiang_2018_baicalein"

  units <- list(
    time          = "h",
    dosing        = "uM (static baicalein covariate; no administered events)",
    concentration = "pg/mL (TNF-alpha, IL-6), unitless ratio (iNOS), uM (NO)"
  )

  covariateData <- list(
    CONC_BAI_UM = list(
      description        = paste(
        "Static baicalein concentration in the cell-culture medium",
        "(time-invariant per experimental design)."
      ),
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-invariant per Xiang 2018 Materials and Methods",
        "('Measurement of the Production of TNF-alpha, IL-6 and NO'):",
        "RAW264.7 cells were pretreated with baicalein for 0.5 h",
        "before LPS stimulation and remained exposed to the same",
        "baicalein concentration for the rest of the 24 h experiment.",
        "Tested values: 0 (control), 10, 20, and 40 uM (Materials and",
        "Methods; Figures 3-4). The baicalein-treated PD model in Eq 1",
        "(see model() block) enters this covariate via the log-linear",
        "inhibition f(Bai) = alpha * log(CONC_BAI_UM + 1); the +1",
        "shift encodes a control (CONC_BAI_UM = 0 -> f = 0) without",
        "the undefined ln(0) (see vignette Errata for the choice).",
        "Registered as the canonical CONC_BAI_UM in",
        "inst/references/covariate-columns.md, an in-vitro applied",
        "drug-concentration covariate (the CONC_<DRUG>_<UNITS> family",
        "is distinct from the state-derived plasma concentration Cc",
        "and the CP_<DRUG> plasma-PD-driver family)."
      ),
      source_name        = "C_Bai (Xiang 2018 Eq 2)"
    )
  )

  population <- list(
    species        = "in vitro (RAW264.7 mouse macrophage cell line)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    disease_state  = paste(
      "LPS-induced inflammation in murine RAW264.7 macrophages",
      "(Center of Cellular Resources, Chinese Academy of Sciences,",
      "Shanghai). Cells were cultured in DMEM + 10% FBS + 1%",
      "penicillin-streptomycin at 37 C in 5% CO2 and stimulated with",
      "1 ug/mL LPS from Escherichia coli O55:B5 (Sigma) after a 0.5 h",
      "pretreatment with baicalein at 10, 20, or 40 uM (control = 0 uM)."
    ),
    dose_range     = paste(
      "Static covariate concentrations (no PK dosing events):",
      "baicalein 0, 10, 20, 40 uM and LPS 1 ug/mL throughout the",
      "24.5 h experiment. TNF-alpha and IL-6 measured by ELISA",
      "(MULTI SCIENCES Mouse kits, intra/inter-assay CV 2.3-3.0%",
      "and 2.1-5.3%), iNOS by Western blot (normalized to the t = 0",
      "control as (iNOS/GAPDH) / (iNOS0/GAPDH0)), NO by the Griess",
      "reaction (sodium nitrite calibration). Sampling times after LPS",
      "stimulation: 1, 2, 4, 8, 12, and 24 h; for the PD model the",
      "1.5, 2.5, 4.5, 8.5, and 12.5 h baicalein-treated points were",
      "used (24.5 h excluded per Materials and Methods)."
    ),
    notes          = paste(
      "Estimation in Monolix 2016R1 (Lixoft, Antony, France) with",
      "FOCEELS. Final model retained kinTNF, koutTNF, kinIL6, koutIL6,",
      "tau1, kiniNOS, tau2, kinNO, alpha, and delta as estimated",
      "typical-value population parameters; koutiNOS and koutNO were",
      "fixed to zero to stabilise the post-12.5 h plateau in the",
      "iNOS / NO data (Results, 'Model Simulation'). The reported",
      "CV% column in Table 1 is the precision of the estimate (RSE),",
      "not between-subject variability -- the model is a typical-value",
      "mechanism without etas or a residual-error structure; the per-",
      "output proportional residual placeholders below carry a fixed",
      "10% magnitude so the model is fit-compatible (see vignette",
      "Errata for the rationale)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Estimated typical-value mechanism parameters (Xiang 2018 Table 1).
    # ------------------------------------------------------------------
    lalpha <- log(0.0832)
    label("alpha -- log-linear coefficient of baicalein inhibition on TNF-alpha production (unitless)")
    # Xiang 2018 Table 1, row 1: a = 0.0832, RSE 15%. The Discussion
    # rephrases the magnitude as "TNF-alpha inhibition increases by about
    # 8.32% for each degree rise in baicalein concentration".

    lkin_tnf <- log(3740)
    label("kinTNF -- LPS-stimulated zero-order TNF-alpha production rate (pg/mL/h)")
    # Xiang 2018 Table 1, row 2: kinTNFa = 3.74 x 10^3, RSE 6%.
    # Table 1 lists the unit as "h^-1" but the equation
    # dTNFa/dt = kinTNFa * (1 - f(Bai)) - koutTNFa * TNFa
    # requires kinTNFa to have units of [TNFa]/time = pg/mL/h.
    # Interpreting the Table 1 unit label as a typesetting slip;
    # see vignette Errata.

    lkout_tnf <- log(0.0463)
    label("koutTNF -- first-order TNF-alpha elimination rate (1/h)")
    # Xiang 2018 Table 1, row 3: koutTNFa = 0.0463 h^-1, RSE 32%.

    lkin_il6 <- log(0.353)
    label("kinIL6 -- TNF-alpha-driven IL-6 production rate constant (1/h)")
    # Xiang 2018 Table 1, row 4: kinIL-6 = 0.353 h^-1, RSE 9%.
    # Multiplies the delayed TNF-alpha concentration (Eq 4).

    lkout_il6 <- log(0.143)
    label("koutIL6 -- first-order IL-6 elimination rate (1/h)")
    # Xiang 2018 Table 1, row 5: koutIL-6 = 0.143 h^-1, RSE 23%.

    ltau1 <- log(1.38)
    label("tau1 -- lag time in IL-6 production from TNF-alpha (h)")
    # Xiang 2018 Table 1, row 6: tau1 = 1.38 h, RSE 11%. Encoded
    # via a single-compartment delay state with mean transit time
    # tau1 (see vignette Errata for the DDE-approximation note).

    lkin_inos <- log(0.00169)
    label("kiniNOS -- TNF-alpha-driven iNOS production rate constant (1/h)")
    # Xiang 2018 Table 1, row 7: kiniNOS = 0.00169 h^-1, RSE 9%.
    # Multiplies the delayed TNF-alpha concentration (Eq 6); iNOS is
    # in relative units (normalized to the t = 0 control) so the
    # apparent unit of kiniNOS is (relative iNOS) / (pg/mL * h).

    ltau2 <- log(1.41)
    label("tau2 -- lag time in iNOS production from TNF-alpha (h)")
    # Xiang 2018 Table 1, row 8: tau2 = 1.41 h, RSE 13%. Same
    # single-compartment delay encoding as tau1.

    lkin_no <- log(0.0605)
    label("kinNO -- iNOS-driven NO production rate (uM/h at iNOS^delta = 1)")
    # Xiang 2018 Table 1, row 9: kinNO = 0.0605, RSE 37%. Table 1
    # lists the unit as "h^-1" but Eq 7 (dNO/dt = kinNO * iNOS^delta
    # - koutNO * NO) with iNOS unitless (relative) requires kinNO
    # to have units uM/h. Same Table-1 unit-label slip as kinTNFa;
    # see vignette Errata.

    ldelta <- log(1.35)
    label("delta -- exponential amplification of iNOS for NO production (unitless)")
    # Xiang 2018 Table 1, row 10: delta = 1.35, RSE 9%. NO is
    # amplified from iNOS as iNOS^delta in Eq 7.

    # ------------------------------------------------------------------
    # Fixed parameters (held at zero per source for the post-12.5 h
    # plateau; see Results paragraph 'The observed data of iNOS
    # expression and NO level maintained a proper balance after 12.5 h
    # treatment of baicalein, koutiNOS, and koutNO were fixed to 0 to
    # reduce the model volatility').
    # ------------------------------------------------------------------
    kout_inos <- fixed(0)
    label("koutiNOS -- iNOS elimination rate (fixed to 0 per source)")

    kout_no <- fixed(0)
    label("koutNO -- NO elimination rate (fixed to 0 per source)")

    # ------------------------------------------------------------------
    # Baselines = control (t = 0) measurements per Eqs 3 and 5
    # (TNFa_0 = Control(TNFa); IL-6_0 = Control(IL-6)) and the iNOS
    # Western-blot normalization (iNOS_0 = (iNOS/GAPDH) /
    # (iNOS0/GAPDH0) = 1 at t = 0). The numeric baselines are not
    # tabulated in Table 1; they are approximate digitisations from
    # Figure 4 (control, t = 0). Users supplying their own measured
    # baselines should override these defaults; the figure-derived
    # defaults are documented in the vignette Errata.
    # ------------------------------------------------------------------
    lrbase_tnf <- log(1000)
    label("Baseline TNF-alpha concentration at t = 0 in control wells (pg/mL)")
    # Figure 4A, control curve at t = 0 (digitised).

    lrbase_il6 <- log(100)
    label("Baseline IL-6 concentration at t = 0 in control wells (pg/mL)")
    # Figure 4B, control curve at t = 0 (digitised; IL-6 is near zero
    # before LPS stimulation in RAW264.7 macrophages).

    lrbase_inos <- log(1)
    label("Baseline iNOS expression at t = 0 (relative, definitionally 1)")
    # Materials and Methods 'Western Blotting': iNOS is reported as
    # (iNOS/GAPDH) / (iNOS0/GAPDH0), so iNOS(0) = 1 by definition.

    lrbase_no <- log(1)
    label("Baseline NO concentration at t = 0 in control wells (uM)")
    # Figure 4D, control curve at t = 0 (digitised; NO is near zero
    # before LPS stimulation).

    # ------------------------------------------------------------------
    # Residual-error placeholders. Xiang 2018 Methods describe
    # 'diverse weighting strategies (additive, power, multiplicative)
    # were all based on goodness-of-fit plots, AIC, comparison of the
    # observed and individual predicted concentrations, and the
    # coefficient of variation (CV) in parameter estimations' -- the
    # final-model weighting choice is not stated. Placeholder
    # proportional residual SD of 0.10 (10%) on each of the four
    # outputs so the model is fit-compatible in nlmixr2; see vignette
    # Errata for the rationale.
    # ------------------------------------------------------------------
    propSd_tnf <- fixed(0.10)
    label("Placeholder proportional residual SD on TNF-alpha (unreported in source)")

    propSd_il6 <- fixed(0.10)
    label("Placeholder proportional residual SD on IL-6 (unreported in source)")

    propSd_inos <- fixed(0.10)
    label("Placeholder proportional residual SD on iNOS (unreported in source)")

    propSd_no <- fixed(0.10)
    label("Placeholder proportional residual SD on NO (unreported in source)")
  })

  model({
    # ================================================================
    # 1. Back-transform estimated structural parameters from log scale.
    # ================================================================
    alpha      <- exp(lalpha)
    kin_tnf    <- exp(lkin_tnf)
    kout_tnf   <- exp(lkout_tnf)
    kin_il6    <- exp(lkin_il6)
    kout_il6   <- exp(lkout_il6)
    tau1       <- exp(ltau1)
    kin_inos   <- exp(lkin_inos)
    tau2       <- exp(ltau2)
    kin_no     <- exp(lkin_no)
    delta      <- exp(ldelta)
    rbase_tnf  <- exp(lrbase_tnf)
    rbase_il6  <- exp(lrbase_il6)
    rbase_inos <- exp(lrbase_inos)
    rbase_no   <- exp(lrbase_no)

    # ================================================================
    # 2. Baicalein log-linear inhibition function (Xiang 2018 Eq 2).
    #    Paper writes f(Bai) = alpha * ln(C_Bai); the bare ln is
    #    undefined at C_Bai = 0 (control wells with no baicalein).
    #    Encoded as alpha * log(CONC_BAI_UM + 1) so f(0) = 0 in
    #    control wells. The +1 shift differs from the bare ln by
    #    less than 5% at the three tested concentrations (10, 20,
    #    40 uM); see vignette Errata.
    # ================================================================
    f_bai <- alpha * log(CONC_BAI_UM + 1)

    # ================================================================
    # 3. TNF-alpha indirect response (Xiang 2018 Eq 1):
    #    LPS-stimulated zero-order production kin_tnf scaled by
    #    (1 - f(Bai)) less first-order TNF-alpha elimination.
    # ================================================================
    d/dt(tnf) <- kin_tnf * (1 - f_bai) - kout_tnf * tnf

    # ================================================================
    # 4. Single-compartment delay states for tau1 (TNF -> IL-6) and
    #    tau2 (TNF -> iNOS). The DDE TNFa(t - tau) in Eqs 4 and 6 is
    #    approximated by a one-compartment first-order distribution
    #    with mean transit time tau (impulse response is a single-
    #    exponential ramp, not a Dirac shift). rxode2 does not
    #    provide a native delay-differential-equation solver; see
    #    vignette Errata for the impulse-response implications and
    #    the alternative transit-chain encoding considered.
    # ================================================================
    d/dt(transit1) <- (tnf - transit1) / tau1
    d/dt(transit2) <- (tnf - transit2) / tau2

    # ================================================================
    # 5. IL-6 indirect response (Xiang 2018 Eq 4).
    # ================================================================
    d/dt(il6) <- kin_il6 * transit1 - kout_il6 * il6

    # ================================================================
    # 6. iNOS indirect response (Xiang 2018 Eq 6). koutiNOS held at
    #    0 per source -- the iNOS state grows from the t = 0 baseline
    #    and does not return; the Discussion attributes the absence
    #    of post-12.5 h decay to the (iNOS/GAPDH) / (iNOS0/GAPDH0)
    #    normalisation accumulating measurement bias rather than to
    #    true biological non-degradation.
    # ================================================================
    d/dt(inos) <- kin_inos * transit2 - kout_inos * inos

    # ================================================================
    # 7. NO indirect response (Xiang 2018 Eq 7). iNOS^delta is the
    #    paper's amplifying-effect form and is unitless because iNOS
    #    is unitless; kin_no carries the uM/h units of dNO/dt.
    #    koutNO held at 0 per source (same plateau rationale as
    #    koutiNOS).
    # ================================================================
    d/dt(no) <- kin_no * inos^delta - kout_no * no

    # ================================================================
    # 8. Initial conditions: control (t = 0) measurements. Setting
    #    both delay compartments to bl_tnf at t = 0 makes the delay
    #    states consistent with the baseline TNF-alpha state, so the
    #    delayed signal at t = 0 equals the t = 0 baseline (no
    #    spurious step at t = 0+).
    # ================================================================
    tnf(0)        <- rbase_tnf
    transit1(0) <- rbase_tnf
    transit2(0) <- rbase_tnf
    il6(0)        <- rbase_il6
    inos(0)       <- rbase_inos
    no(0)         <- rbase_no

    # ================================================================
    # 9. Multi-output observations. Each cytokine / mediator is its
    #    own output. The residual-error model is the unreported-source
    #    placeholder (10% proportional on each output) carried in
    #    ini() so the model is nlmixr2-fit-compatible; downstream
    #    users who refit the model should replace propSd_<output>
    #    with whatever residual-error structure their data support.
    # ================================================================
    tnf  ~ prop(propSd_tnf)
    il6  ~ prop(propSd_il6)
    inos ~ prop(propSd_inos)
    no   ~ prop(propSd_no)
  })
}
