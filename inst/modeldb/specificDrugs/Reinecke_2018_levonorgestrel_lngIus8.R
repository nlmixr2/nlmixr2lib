Reinecke_2018_levonorgestrel_lngIus8 <- function() {
  description <- paste(
    "LNG-IUS 8 (Jaydess/Skyla) arm of the Reinecke 2018 final comprehensive",
    "levonorgestrel (LNG) population PK model. In-vivo release from the",
    "intrauterine reservoir uses two processes: a zero-order term c12 and a",
    "time-dependent first-order term c13 * depot / (t1 + t), a form the",
    "authors justify mechanistically as release through the membrane and",
    "through the open ends of the device. Released LNG enters a",
    "two-compartment disposition model in which only unbound drug is",
    "eliminated or distributed, the free fraction fuLNG being the closed-form",
    "solution of reversible binding to SHBG (KDS = 1.82 nmol/L) and to",
    "albumin (KDA = 18209 nmol/L, constant ALB = 700000 nmol/L). SHBG is an",
    "indirect-response turnover state whose synthesis is linearly inhibited",
    "by a delay-compartment-smoothed molar LNG signal. Body weight raises",
    "clearance and lowers the SHBG baseline as power functions centred on",
    "66 kg. Disposition and SHBG parameters are shared verbatim with the",
    "LNG-IUS 20 and LNG-IUS 12 arms because all three devices, plus the",
    "intravenous and oral data, were fitted in one comprehensive NONMEM run."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S2, 'final comprehensive",
    "model' block; model equations from Supplemental Table S3a.",
    "Companion arms: modellib('Reinecke_2018_levonorgestrel_lngIus20') and",
    "modellib('Reinecke_2018_levonorgestrel_lngIus12')."
  )
  vignette <- "Reinecke_2018_levonorgestrel_contraceptives"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel loaded in the intrauterine reservoir)",
    concentration = "ng/L (total and unbound LNG in serum); SHBG in nmol/L; residual device content in mg"
  )

  compartmentData <- list(
    depot       = list(analyte = "levonorgestrel", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "levonorgestrel", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "levonorgestrel", units = "mg", specimen = "serum", verified = TRUE),
    effect      = list(analyte = "levonorgestrel", units = "nmol/L", specimen = "not applicable", verified = TRUE),
    shbg        = list(analyte = "SHBG", units = "nmol/L", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg). Power-function covariate on clearance of unbound levonorgestrel (increasing with weight) and on the SHBG baseline concentration (decreasing with weight).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, centred on 66 kg: the median of the pooled",
        "intrauterine-system data (Table 1 footnote a), the value Table S2",
        "attaches to every F_LNG-IUS row ('median[WGHT]=66kg, IUS data') and",
        "the value the Implant Model section names explicitly. The reported",
        "CL exponent 0.823 is close to the three-quarter power expected from",
        "allometry (Discussion)."
      ),
      source_name        = "WGHT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline (years). Screened as a power-function covariate on levonorgestrel clearance and on the SHBG baseline concentration.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the comprehensive covariate step but not retained; only",
        "body weight reached the P <= .001 forward-inclusion criterion",
        "(Results, 'Comprehensive Model'). No point estimate is reported."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1245L,
    n_studies      = 1L,
    studies        = "Phase 3 study 310442 (levonorgestrel and SHBG concentrations plus residual-content measurements over the full 3 years of intended treatment)",
    age_range      = "18-35 years (median 28 years)",
    weight_range   = "39-137 kg (median 66 kg)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women using a levonorgestrel-releasing intrauterine system for contraception",
    dose_range     = "Single insertion of one LNG-IUS 8 (Jaydess/Skyla) device; nominal release 8 ug/day, indicated for up to 3 years",
    notes          = paste(
      "Demographics and observation counts from Table 1. Contributing",
      "observations: 1457 LNG concentrations, 1612 SHBG concentrations and",
      "763 residual-content measurements. LNG-IUS 8 and LNG-IUS 12 were",
      "both studied in study 310442, over 3 and 5 years of intended",
      "treatment respectively. The loaded reservoir content is not stated",
      "in the text; the y-intercept of the Supplemental Figure S1c",
      "residual-content plot is approximately 13.5 mg, and that value",
      "reproduces the Table 4 in-vivo release rates to within 1%. Nineteen",
      "LNG and 45 SHBG measurements taken after device removal were",
      "excluded across the IUS data before the stochastic model would",
      "converge. Eta and epsilon shrinkage were below 18%."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # In-vivo release from the LNG-IUS 8 reservoir. Table S3a gives the
    # LNG-IUS 12 / LNG-IUS 8 / implant branch as
    #   input1 = C12
    #   input2 = C13 * A(2) / (T1 + T)
    # i.e. a zero-order release plus a time-dependent first-order process.
    # Table S2 reports C12 in ng/h; the 1e-6 conversion to the mg/h scale
    # of the device compartment is applied inside model() so that ini()
    # carries the tabulated number. C13 is dimensionless.
    # ------------------------------------------------------------------
    lc12 <- log(222)     ; label("Zero-order release rate C12 (ng/h)")                          # Table S2, final comprehensive model, C12 LNG-IUS 8 = 222 (219-225) ng/h
    lc13 <- log(0.0257)  ; label("Second release process coefficient C13 (dimensionless)")      # Table S2, final comprehensive model, C13 LNG-IUS 8 = 0.0257 (0.0223-0.0291)
    lt1  <- log(756)     ; label("Time parameter T1 for the drop in release process 2 (h)")     # Table S2, final comprehensive model, T1 LNG-IUS 8 = 756 (540-972)

    # Absolute bioavailability on the logit scale: TF = 4.24 gives
    # F = 0.986, matching the value Table S2 prints in parentheses.
    logitfdepot <- 4.24  ; label("Logit-transformed absolute bioavailability of the LNG-IUS 8 reservoir (logit units)")  # Table S2, final comprehensive model, F LNG-IUS 8 = 4.24 (3.97-4.51), i.e. 0.986 (0.981-0.989)

    # Levonorgestrel disposition, shared across all arms of the
    # comprehensive fit. cl and q are clearances of UNBOUND drug.
    lvc <- log(20.7)  ; label("Central volume of distribution Vc (L)")                          # Table S2, final comprehensive model, Vc = 20.7 (18.7-22.7)
    lcl <- log(243)   ; label("Clearance of unbound levonorgestrel CL (L/h)")                   # Table S2, final comprehensive model, CL = 243 (240-246)
    lvp <- log(4690)  ; label("Peripheral volume of distribution Vp (L)")                       # Table S2, final comprehensive model, Vp = 4690 (3910-5470)
    lq  <- log(600)   ; label("Intercompartmental clearance of unbound levonorgestrel Q (L/h)") # Table S2, final comprehensive model, Q = 600 (534-666)

    # SHBG turnover. The intrauterine baseline applies here; tau, RI and
    # KOUT were fixed to the final iv-oral estimates (Table S2 footnote f).
    lrbase_shbg <- log(51.5)             ; label("Baseline SHBG serum concentration SBL_IUS at 66 kg (nmol/L)")     # Table S2, final comprehensive model, SBL_IUS = 51.5 (50.6-52.4)
    lkout_shbg  <- fixed(log(0.00313))   ; label("First-order elimination rate constant of SHBG KOUT (1/h)")        # Table S2, final comprehensive model, KOUT = 0.00313, fixed to the final iv-oral value (footnote f)
    ltau        <- fixed(log(13.7))      ; label("Delay time constant tau of the levonorgestrel effect on SHBG (h)")# Table S2, final comprehensive model, tau = 13.7, fixed to the final iv-oral value (footnote f)
    lri         <- fixed(log(0.232))     ; label("Linear inhibition factor RI of SHBG synthesis by delayed LNG (L/nmol)") # Table S2, final comprehensive model, R_I = 0.232, fixed to the final iv-oral value (footnote f)

    # Body-weight power exponents.
    e_wt_cl         <-  0.823  ; label("Body-weight power exponent on clearance, (WT/66)^e_wt_cl")                 # Table S2, final comprehensive model, C_WGHT,CL = 0.823 (0.770-0.876)
    e_wt_rbase_shbg <- -0.950  ; label("Body-weight power exponent on the SHBG baseline, (WT/66)^e_wt_rbase_shbg") # Table S2, final comprehensive model, C_WGHT,SBL = -0.950 (-1.04 to -0.863)

    # Inter-individual variability; diagonals inverted from the reported
    # CV(%) via log(1 + CV^2) per Table S2 footnote d, off-diagonal
    # reported untransformed (correlation -0.788).
    etalcl + etalrbase_shbg ~ c(0.039221,
                                -0.0655, 0.176237)   # Table S2, final comprehensive model: IIV_CL 20.0% (17.9-22.0), IIV_SBL 43.9% (42.0-45.7), COV_IIV -0.0655 (-0.0741 to -0.0569)

    # Residual error; the additive residual-content entry in Table S2 is on
    # the variance scale and is square-rooted here (see the vignette
    # Errata for the Table S1 cross-check that pins the scale).
    propSd            <- 0.235          ; label("Proportional residual error SD on total levonorgestrel concentration (fraction)") # Table S2, final comprehensive model, SIGMA_LNG (CV) = 23.5% (21.6-25.3)
    propSd_shbg       <- 0.229          ; label("Proportional residual error SD on SHBG concentration (fraction)")                 # Table S2, final comprehensive model, SIGMA_SHBG (CV) = 22.9% (21.5-24.3)
    addSd_iusResidual <- sqrt(0.0541)   ; label("Additive residual error SD on residual device content (mg)")                      # Table S2, final comprehensive model, SIGMA LNG-IUS 8 res. content = 0.0541 (0.0360-0.0722) as a variance; sqrt = 0.233 mg
  })

  model({
    # Physicochemical constants for the closed-form free fraction of LNG
    # (Table S3a, "Free fraction of LNG").
    MWLNG <- 312.5      # Table S3a: molecular weight LNG, g/mol
    KDS   <- 1.82       # Table S3a: dissociation constant SHBG, nmol/L
    KDA   <- 18209      # Table S3a: dissociation constant albumin, nmol/L
    ALB   <- 700000     # Table S3a: albumin concentration, nmol/L
    KALB  <- 1 + ALB / KDA

    wtRef <- 66

    # Release parameters. c12 converts ng/h to the mg/h scale of the
    # device compartment; c13 is dimensionless.
    c12 <- exp(lc12) * 1e-6
    c13 <- exp(lc13)
    t1  <- exp(lt1)

    vc <- exp(lvc)
    cl <- exp(lcl + e_wt_cl * log(WT / wtRef) + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)

    k20 <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    tau        <- exp(ltau)
    ri         <- exp(lri)
    rbase_shbg <- exp(lrbase_shbg + e_wt_rbase_shbg * log(WT / wtRef) + etalrbase_shbg)
    kout_shbg  <- exp(lkout_shbg)
    kin_shbg   <- rbase_shbg * kout_shbg

    fdepot <- exp(logitfdepot) / (1 + exp(logitfdepot))

    A3nM  <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    inh <- ri * effect

    # LNG-IUS 8 release branch of Table S3a.
    input1 <- c12
    input2 <- c13 * depot / (t1 + t)

    d/dt(depot)       <- -input1 - input2
    d/dt(central)     <-  input1 + input2 - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * A3nM - (1 / tau) * effect
    d/dt(shbg)        <-  kin_shbg * (1 - inh) - kout_shbg * shbg

    f(depot) <- fdepot
    shbg(0)  <- rbase_shbg

    Cc          <- (central / vc) * 1e6
    CcUnbound   <- fuLNG * Cc
    iusResidual <- depot
    releaseRate <- (input1 + input2) * 1000 * 24

    Cc          ~ prop(propSd)
    shbg        ~ prop(propSd_shbg)
    iusResidual ~ add(addSd_iusResidual)
  })
}
