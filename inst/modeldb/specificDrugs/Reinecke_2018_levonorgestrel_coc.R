Reinecke_2018_levonorgestrel_coc <- function() {
  description <- paste(
    "Final combined-oral-contraceptive (COC) model of the Reinecke 2018",
    "integrated levonorgestrel (LNG) population PK analysis, for Miranova",
    "(100 ug LNG + 20 ug ethinylestradiol once daily, 21 days on / 7 days",
    "off). The levonorgestrel/SHBG structure is the final iv-oral model:",
    "two-compartment LNG disposition in which only unbound drug is",
    "eliminated or distributed, with the free fraction fuLNG given by the",
    "closed-form solution of reversible binding to SHBG (KDS = 1.82 nmol/L)",
    "and albumin (KDA = 18209 nmol/L, constant ALB = 700000 nmol/L), coupled",
    "to an SHBG indirect-response turnover state. Co-formulated",
    "ethinylestradiol is carried as a three-compartment PK model with",
    "first-order absorption and a lag time, all of its parameters fixed to a",
    "previous ethinylestradiol population analysis, and drives a linear",
    "STIMULATION of SHBG synthesis through its own delay compartment. Two",
    "further ethinylestradiol effects were estimated empirically rather than",
    "mechanistically: unbound LNG clearance is 37% lower than in the",
    "levonorgestrel-only models (175 vs 279 L/h), consistent with",
    "inactivation of CYP3A4 by ethinylestradiol metabolites, and the SHBG",
    "elimination rate is nearly doubled."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S2, 'final COC model' block;",
    "model equations from Supplemental Table S3a. Parameters not listed in",
    "that block are fixed to the final iv-oral model, available as",
    "modellib('Reinecke_2018_levonorgestrel_ivOral'). The fixed",
    "ethinylestradiol disposition parameters come from an earlier",
    "ethinylestradiol / gestodene population analysis cited as reference 38",
    "of the paper."
  )
  vignette <- "Reinecke_2018_levonorgestrel_contraceptives"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel and ethinylestradiol)",
    concentration = "ng/L (total and unbound LNG, and ethinylestradiol, in serum); SHBG in nmol/L"
  )

  compartmentData <- list(
    depot                       = list(analyte = "levonorgestrel",   units = "mg",     specimen = "administration site", verified = TRUE),
    central                     = list(analyte = "levonorgestrel",   units = "mg",     specimen = "serum",               verified = TRUE),
    peripheral1                 = list(analyte = "levonorgestrel",   units = "mg",     specimen = "serum",               verified = TRUE),
    effect                      = list(analyte = "levonorgestrel",   units = "nmol/L", specimen = "not applicable",      verified = TRUE),
    shbg                        = list(analyte = "SHBG",             units = "nmol/L", specimen = "serum",               verified = TRUE),
    depot_ethinylestradiol      = list(analyte = "ethinylestradiol", units = "mg",     specimen = "administration site", verified = TRUE),
    central_ethinylestradiol    = list(analyte = "ethinylestradiol", units = "mg",     specimen = "serum",               verified = TRUE),
    peripheral1_ethinylestradiol = list(analyte = "ethinylestradiol", units = "mg",    specimen = "serum",               verified = TRUE),
    peripheral2_ethinylestradiol = list(analyte = "ethinylestradiol", units = "mg",    specimen = "serum",               verified = TRUE),
    effect_ethinylestradiol     = list(analyte = "ethinylestradiol", units = "nmol/L", specimen = "not applicable",      verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg). Power-function covariate on the absolute bioavailability of levonorgestrel from the combined oral contraceptive; bioavailability decreases with increasing body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, centred on 64 kg, the median of the COC study",
        "population (Table 1, study 94011). The reference weight is pinned",
        "by Table S2 itself: F_COC = 0.751 at median[WGHT] = 64 kg and the",
        "derived values in the same row block (0.791 at 62 kg, 0.813 at",
        "61 kg) are reproduced to three decimal places by",
        "0.751 * (WT/64)^-1.64. No body-weight effect on clearance or on the",
        "SHBG baseline was carried into the COC model."
      ),
      source_name        = "WGHT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 18L,
    n_studies      = 2L,
    studies        = c(
      "Phase 1 study 94011 (combined oral contraceptive Miranova, 100 ug LNG + 20 ug ethinylestradiol once daily, dense LNG and SHBG profiles over 3 treatment cycles, 18 women)",
      "Phase 1 study 92085 (single-dose intravenous 90 ug LNG; included in every model-development step so that absolute bioavailability is estimable)"
    ),
    age_range      = "20-34 years (median 30 years)",
    weight_range   = "51-83 kg (median 64 kg)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women",
    dose_range     = "100 ug levonorgestrel plus 20 ug ethinylestradiol orally once daily in a 28-day cycle (21 days on, 7 days off)",
    co_medication  = "Ethinylestradiol, co-formulated in the same tablet, is modelled explicitly",
    notes          = paste(
      "Demographics from Table 1. Contributing observations: 1166 LNG and",
      "486 SHBG concentrations from study 94011 plus 286 LNG concentrations",
      "from the intravenous study. Ethinylestradiol concentrations were not",
      "fitted because more than 50% were below the lower limit of",
      "quantitation; only the effect of ethinylestradiol on SHBG was",
      "estimated, on top of a fixed ethinylestradiol PK model. No additional",
      "random effect was implemented for ethinylestradiol. Eta and epsilon",
      "shrinkage were below 6%."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Levonorgestrel parameters re-estimated for the combined oral
    # contraceptive. Unbound clearance is 37.3% lower than the 279 L/h of
    # the final iv-oral model, with non-overlapping 95% CIs; the authors
    # attribute this to inactivation of CYP3A4 by ethinylestradiol
    # metabolites (Discussion).
    # ------------------------------------------------------------------
    lfdepot     <- log(0.751)  ; label("Absolute bioavailability F_COC of levonorgestrel at 64 kg (fraction)")            # Table S2, final COC model, F_COC = 0.751 (0.687-0.815) at median[WGHT] = 64 kg
    e_wt_fdepot <- -1.64       ; label("Body-weight power exponent on F_COC, (WT/64)^e_wt_fdepot")                        # Table S2, final COC model, C_WGHT,F,COC = -1.64 (-2.31 to -0.966)
    lka         <- log(2.35)   ; label("Oral first-order absorption rate constant Ka of levonorgestrel (1/h)")            # Table S2, final COC model, Ka_COC = 2.35 (1.77-2.93)
    lcl         <- log(175)    ; label("Clearance of unbound levonorgestrel CL (L/h)")                                    # Table S2, final COC model, CL_COC = 175 (155-195)

    # Fixed to the final iv-oral model (Table S2 footnote h).
    lvc <- fixed(log(22.0))    ; label("Central volume of distribution Vc (L)")                                           # Table S2, final iv-oral model, Vc = 22.0; fixed here (footnote h)
    lvp <- fixed(log(5120))    ; label("Peripheral volume of distribution Vp (L)")                                        # Table S2, final iv-oral model, Vp = 5120; fixed here (footnote h)
    lq  <- fixed(log(649))     ; label("Intercompartmental clearance of unbound levonorgestrel Q (L/h)")                  # Table S2, final iv-oral model, Q = 649; fixed here (footnote h)

    # ------------------------------------------------------------------
    # SHBG turnover. The baseline and the elimination rate were both
    # re-estimated for the COC; KOUT is nearly twice the iv-oral value
    # (0.00615 vs 0.00313), which the authors note is consistent with
    # earlier ethinylestradiol / gestodene modelling.
    # ------------------------------------------------------------------
    lrbase_shbg <- log(62.6)      ; label("Baseline SHBG serum concentration SBL_COC (nmol/L)")                           # Table S2, final COC model, SBL_COC = 62.6 (50.3-74.9)
    lkout_shbg  <- log(0.00615)   ; label("First-order elimination rate constant of SHBG KOUT (1/h)")                     # Table S2, final COC model, KOUT,COC = 0.00615 (0.00431-0.00799)
    ltau        <- fixed(log(13.7))  ; label("Delay time constant tau of the levonorgestrel effect on SHBG (h)")          # Table S2, final iv-oral model, tau = 13.7; fixed here (footnote h)
    lri         <- fixed(log(0.232)) ; label("Linear inhibition factor RI of SHBG synthesis by delayed LNG (L/nmol)")     # Table S2, final iv-oral model, R_I = 0.232; fixed here (footnote h)
    lrs         <- log(17.0)      ; label("Linear stimulation factor RS of SHBG synthesis by delayed ethinylestradiol (L/nmol)")  # Table S2, final COC model, R_S = 17.0 (13.2-20.8)

    # ------------------------------------------------------------------
    # Ethinylestradiol disposition. Every value is fixed to a previous
    # ethinylestradiol population analysis (Table S2 footnote k: "Fixed to
    # value from previous modeling"); only the effect of ethinylestradiol
    # on SHBG (lrs above) was estimated here.
    # ------------------------------------------------------------------
    lfdepot_ethinylestradiol <- fixed(log(0.45))  ; label("Absolute oral bioavailability of ethinylestradiol (fraction)")        # Table S2, final COC model, F_EE = 0.45 (fixed, footnote k)
    ltlag_ethinylestradiol   <- fixed(log(0.39))  ; label("Absorption lag time of ethinylestradiol (h)")                         # Table S2, final COC model, ALAG1_EE = 0.39 (fixed, footnote k)
    lka_ethinylestradiol     <- fixed(log(1.23))  ; label("Oral first-order absorption rate constant of ethinylestradiol (1/h)") # Table S2, final COC model, KA_EE = 1.23 (fixed, footnote k)
    lvc_ethinylestradiol     <- fixed(log(22.5))  ; label("Central volume of distribution of ethinylestradiol (L)")              # Table S2, final COC model, Vc_EE = 22.5 (fixed, footnote k)
    lcl_ethinylestradiol     <- fixed(log(13.5))  ; label("Clearance of ethinylestradiol (L/h)")                                 # Table S2, final COC model, CL_EE = 13.5 (fixed, footnote k)
    lvp_ethinylestradiol     <- fixed(log(211))   ; label("First peripheral volume of distribution of ethinylestradiol (L)")     # Table S2, final COC model, Vp1_EE = 211 (fixed, footnote k)
    lq_ethinylestradiol      <- fixed(log(63.4))  ; label("First intercompartmental clearance of ethinylestradiol (L/h)")        # Table S2, final COC model, Q1_EE = 63.4 (fixed, footnote k)
    lvp2_ethinylestradiol    <- fixed(log(22.5))  ; label("Second peripheral volume of distribution of ethinylestradiol (L)")    # Table S2, final COC model, Vp2_EE = 22.5 (fixed, footnote k)
    lq2_ethinylestradiol     <- fixed(log(79.6))  ; label("Second intercompartmental clearance of ethinylestradiol (L/h)")       # Table S2, final COC model, Q2_EE = 79.6 (fixed, footnote k)
    ltau_ethinylestradiol    <- fixed(log(24))    ; label("Delay time constant of the ethinylestradiol effect on SHBG (h)")      # Table S2, final COC model, tau_EE = 24 (fixed, footnote k)

    # ------------------------------------------------------------------
    # Inter-individual variability on clearance and on the SHBG baseline,
    # with covariance. Diagonals inverted from the reported CV(%) via
    # log(1 + CV^2) per Table S2 footnote d:
    #   CL  : log(1 + 0.270^2) = 0.070365
    #   SBL : log(1 + 0.540^2) = 0.255882
    # implying a correlation of -0.641.
    # ------------------------------------------------------------------
    etalcl + etalrbase_shbg ~ c(0.070365,
                                -0.0860, 0.255882)   # Table S2, final COC model: IIV_CL 27.0% (17.2-34.4), IIV_SBL 54.0% (37.5-68.1), COV_IIV -0.0860 (-0.148 to -0.0245)

    propSd      <- 0.272   ; label("Proportional residual error SD on total levonorgestrel concentration (fraction)")  # Table S2, final COC model, SIGMA_LNG (CV) = 27.2% (24.1-29.9)
    propSd_shbg <- 0.171   ; label("Proportional residual error SD on SHBG concentration (fraction)")                  # Table S2, final COC model, SIGMA_SHBG (CV) = 17.1% (14.1-19.6)
  })

  model({
    # Physicochemical constants for the closed-form free fraction of LNG
    # and for the molar ethinylestradiol signal (Table S3a).
    MWLNG <- 312.5      # Table S3a: molecular weight LNG, g/mol
    MWE   <- 296.4      # Table S3a: molecular weight ethinylestradiol, g/mol
    KDS   <- 1.82       # Table S3a: dissociation constant SHBG, nmol/L
    KDA   <- 18209      # Table S3a: dissociation constant albumin, nmol/L
    ALB   <- 700000     # Table S3a: albumin concentration, nmol/L
    KALB  <- 1 + ALB / KDA

    # Reference body weight for the bioavailability covariate.
    wtRef <- 64

    # Levonorgestrel individual parameters.
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)

    k20 <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    # Ethinylestradiol individual parameters (all fixed).
    ka_ethinylestradiol  <- exp(lka_ethinylestradiol)
    vc_ethinylestradiol  <- exp(lvc_ethinylestradiol)
    cl_ethinylestradiol  <- exp(lcl_ethinylestradiol)
    vp_ethinylestradiol  <- exp(lvp_ethinylestradiol)
    q_ethinylestradiol   <- exp(lq_ethinylestradiol)
    vp2_ethinylestradiol <- exp(lvp2_ethinylestradiol)
    q2_ethinylestradiol  <- exp(lq2_ethinylestradiol)
    tau_ethinylestradiol <- exp(ltau_ethinylestradiol)

    k20e <- cl_ethinylestradiol / vc_ethinylestradiol
    k23e <- q_ethinylestradiol / vc_ethinylestradiol
    k24e <- q2_ethinylestradiol / vc_ethinylestradiol
    k32e <- q_ethinylestradiol / vp_ethinylestradiol
    k42e <- q2_ethinylestradiol / vp2_ethinylestradiol

    # SHBG turnover parameters.
    tau        <- exp(ltau)
    ri         <- exp(lri)
    rs         <- exp(lrs)
    rbase_shbg <- exp(lrbase_shbg + etalrbase_shbg)
    kout_shbg  <- exp(lkout_shbg)
    kin_shbg   <- rbase_shbg * kout_shbg

    # Bioavailabilities.
    fdepot                   <- exp(lfdepot + e_wt_fdepot * log(WT / wtRef))
    fdepot_ethinylestradiol  <- exp(lfdepot_ethinylestradiol)
    tlag_ethinylestradiol    <- exp(ltlag_ethinylestradiol)

    # Free fraction of levonorgestrel; see the ivOral model file for the
    # unit derivation of A3nM.
    A3nM  <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    # Molar ethinylestradiol concentration driving its delay compartment.
    A8nM <- (central_ethinylestradiol / vc_ethinylestradiol) * (1e6 / MWE)

    # Delayed inhibition of SHBG synthesis by LNG and delayed stimulation
    # by ethinylestradiol (Table S3a, "SHBG").
    inh <- ri * effect
    ind <- rs * effect_ethinylestradiol

    # Levonorgestrel ODEs.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * A3nM - (1 / tau) * effect

    # SHBG turnover.
    d/dt(shbg) <- kin_shbg * (1 - inh + ind) - kout_shbg * shbg

    # Ethinylestradiol ODEs.
    d/dt(depot_ethinylestradiol)       <- -ka_ethinylestradiol * depot_ethinylestradiol
    d/dt(central_ethinylestradiol)     <-  ka_ethinylestradiol * depot_ethinylestradiol -
      (k20e + k23e + k24e) * central_ethinylestradiol +
      k32e * peripheral1_ethinylestradiol + k42e * peripheral2_ethinylestradiol
    d/dt(peripheral1_ethinylestradiol) <-  k23e * central_ethinylestradiol - k32e * peripheral1_ethinylestradiol
    d/dt(peripheral2_ethinylestradiol) <-  k24e * central_ethinylestradiol - k42e * peripheral2_ethinylestradiol
    d/dt(effect_ethinylestradiol)      <- (1 / tau_ethinylestradiol) * A8nM -
      (1 / tau_ethinylestradiol) * effect_ethinylestradiol

    f(depot)                  <- fdepot
    f(depot_ethinylestradiol) <- fdepot_ethinylestradiol
    alag(depot_ethinylestradiol) <- tlag_ethinylestradiol
    shbg(0)                   <- rbase_shbg

    Cc                   <- (central / vc) * 1e6
    CcUnbound            <- fuLNG * Cc
    Cc_ethinylestradiol  <- (central_ethinylestradiol / vc_ethinylestradiol) * 1e6

    Cc   ~ prop(propSd)
    shbg ~ prop(propSd_shbg)
  })
}
