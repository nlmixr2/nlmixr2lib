Reinecke_2018_levonorgestrel_lngIus20 <- function() {
  description <- paste(
    "LNG-IUS 20 (Mirena) arm of the Reinecke 2018 final comprehensive",
    "levonorgestrel (LNG) population PK model. In-vivo release from the",
    "intrauterine reservoir uses two processes: a first-order term",
    "c12 * depot and a time-dependent term c13 * (1 + depot / (t1 + t)).",
    "Released LNG enters a two-compartment disposition model in which only",
    "unbound drug is eliminated or distributed, the free fraction fuLNG",
    "being the closed-form solution of reversible binding to SHBG",
    "(KDS = 1.82 nmol/L) and to albumin (KDA = 18209 nmol/L, constant",
    "ALB = 700000 nmol/L). SHBG is an indirect-response turnover state whose",
    "synthesis is linearly inhibited by a delay-compartment-smoothed molar",
    "LNG signal. Body weight raises clearance and lowers the SHBG baseline as",
    "power functions centred on 66 kg. Disposition and SHBG parameters are",
    "shared verbatim with the LNG-IUS 12 and LNG-IUS 8 arms because all three",
    "devices, plus the intravenous and oral data, were fitted in one",
    "comprehensive NONMEM run."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S2, 'final comprehensive",
    "model' block; model equations from Supplemental Table S3a.",
    "Companion arms: modellib('Reinecke_2018_levonorgestrel_lngIus12') and",
    "modellib('Reinecke_2018_levonorgestrel_lngIus8')."
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
        "the value the Implant Model section names explicitly ('centred",
        "around the median used in the comprehensive model, 66 kg').",
        "Missing weights in study 89532 were set to this median by the",
        "authors. The reported CL exponent 0.823 is close to the",
        "three-quarter power expected from allometry (Discussion)."
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
    n_subjects     = 333L,
    n_studies      = 2L,
    studies        = c(
      "Phase 2 study 308901 (239 women, 3-year data with LNG, SHBG and residual-content measurements)",
      "Phase 3 study 89532 (94 women, 5-year data with LNG concentrations and residual content but no SHBG)"
    ),
    age_range      = "18-40 years (median 28 years across the pooled IUS data, N = 2790)",
    weight_range   = "39-160 kg (median 66 kg across the pooled IUS data, N = 2790)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women using a levonorgestrel-releasing intrauterine system for contraception",
    dose_range     = "Single insertion of one LNG-IUS 20 (Mirena) device; nominal release 20 ug/day, indicated for up to 5 years",
    notes          = paste(
      "Subject counts follow Table 1 read together with the Data section:",
      "Table 1's product labels are offset by one row relative to their",
      "values, so the 94-subject / study-89532 row belongs to LNG-IUS 20",
      "(the Data section states 'LNG-IUS 20. Two studies were included. In",
      "study 89532 ...' and footnote a states 'no age and body weight data",
      "available in study 89532 (LNG-IUS 20)'), not to LNG-IUS 12. The",
      "239 + 94 = 333 total is consistent with Table 2, which reports",
      "N = 323 evaluable subjects for 'LNG-IUS 20 (studies 89532 and",
      "308901)'. Age and weight were not recorded in study 89532 and were",
      "imputed to the IUS-data medians of 28 years and 66 kg. Contributing",
      "observations: 500 + 478 LNG concentrations, 510 SHBG concentrations",
      "and 64 + 87 residual-content measurements. Nineteen LNG and 45 SHBG",
      "measurements taken after device removal were excluded before the",
      "stochastic model would converge. Eta and epsilon shrinkage were",
      "below 18%."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # In-vivo release from the LNG-IUS 20 reservoir. Table S3a gives the
    # LNG-IUS 20 branch as
    #   input1 = C12 * A(2)
    #   input2 = C13 * (1 + A(2) / (T1 + T))
    # i.e. a first-order release plus a time-dependent process. Table S2
    # reports C12 in units of 1e-6 / h and C13 in units of 1e-6 (an amount
    # per hour on the mg scale of the device compartment); both scale
    # factors are applied inside model() so that the ini() values match the
    # tabulated numbers verbatim.
    # ------------------------------------------------------------------
    lc12 <- log(8.57)   ; label("First-order release rate coefficient C12 (x 1e-6 per hour)")            # Table S2, final comprehensive model, C12 LNG-IUS 20 = 8.57 (7.14-1.0 [sic]) x 1e-6/h
    lc13 <- log(303)    ; label("Second release process coefficient C13 (x 1e-6 mg/h)")                  # Table S2, final comprehensive model, C13 LNG-IUS 20 = 303 (254-352) x 1e-6
    lt1  <- log(61.2)   ; label("Time parameter T1 for the drop in release process 2 (h)")               # Table S2, final comprehensive model, T1 LNG-IUS 20 = 61.2 (51.9-70.5)

    # ------------------------------------------------------------------
    # Absolute bioavailability of the released levonorgestrel, estimated on
    # the logit scale: F = exp(TF) / (1 + exp(TF)). TF = 3.51 gives
    # F = 0.971, matching the value Table S2 prints in parentheses.
    # ------------------------------------------------------------------
    logitfdepot <- 3.51 ; label("Logit-transformed absolute bioavailability of the LNG-IUS 20 reservoir (logit units)")  # Table S2, final comprehensive model, F LNG-IUS 20 = 3.51 (3.19-3.83), i.e. 0.971 (0.960-0.979)

    # ------------------------------------------------------------------
    # Levonorgestrel disposition, shared across all arms of the
    # comprehensive fit. cl and q are clearances of UNBOUND drug: the ODEs
    # multiply them by fuLNG.
    # ------------------------------------------------------------------
    lvc <- log(20.7)  ; label("Central volume of distribution Vc (L)")                          # Table S2, final comprehensive model, Vc = 20.7 (18.7-22.7)
    lcl <- log(243)   ; label("Clearance of unbound levonorgestrel CL (L/h)")                   # Table S2, final comprehensive model, CL = 243 (240-246)
    lvp <- log(4690)  ; label("Peripheral volume of distribution Vp (L)")                       # Table S2, final comprehensive model, Vp = 4690 (3910-5470)
    lq  <- log(600)   ; label("Intercompartmental clearance of unbound levonorgestrel Q (L/h)") # Table S2, final comprehensive model, Q = 600 (534-666)

    # ------------------------------------------------------------------
    # SHBG turnover. The comprehensive model estimated the SHBG baseline
    # separately for intrauterine (51.5 nmol/L) and iv/oral (76.4 nmol/L)
    # treatment because the baseline "is known to vary between studies";
    # the intrauterine value applies here. tau, RI and KOUT could not be
    # estimated once the IUS data were added and were fixed to the final
    # iv-oral estimates (Table S2 footnote f).
    # ------------------------------------------------------------------
    lrbase_shbg <- log(51.5)             ; label("Baseline SHBG serum concentration SBL_IUS at 66 kg (nmol/L)")     # Table S2, final comprehensive model, SBL_IUS = 51.5 (50.6-52.4)
    lkout_shbg  <- fixed(log(0.00313))   ; label("First-order elimination rate constant of SHBG KOUT (1/h)")        # Table S2, final comprehensive model, KOUT = 0.00313, fixed to the final iv-oral value (footnote f)
    ltau        <- fixed(log(13.7))      ; label("Delay time constant tau of the levonorgestrel effect on SHBG (h)")# Table S2, final comprehensive model, tau = 13.7, fixed to the final iv-oral value (footnote f)
    lri         <- fixed(log(0.232))     ; label("Linear inhibition factor RI of SHBG synthesis by delayed LNG (L/nmol)") # Table S2, final comprehensive model, R_I = 0.232, fixed to the final iv-oral value (footnote f)

    # Body-weight power exponents.
    e_wt_cl         <-  0.823  ; label("Body-weight power exponent on clearance, (WT/66)^e_wt_cl")            # Table S2, final comprehensive model, C_WGHT,CL = 0.823 (0.770-0.876)
    e_wt_rbase_shbg <- -0.950  ; label("Body-weight power exponent on the SHBG baseline, (WT/66)^e_wt_rbase_shbg") # Table S2, final comprehensive model, C_WGHT,SBL = -0.950 (-1.04 to -0.863)

    # ------------------------------------------------------------------
    # Inter-individual variability on clearance and on the SHBG baseline,
    # with covariance between the two etas. Table S2 reports the diagonals
    # as CV(%) = sqrt(exp(OMEGA^2) - 1) * 100 (footnote d); inverting,
    #   CL  : log(1 + 0.200^2) = 0.039221
    #   SBL : log(1 + 0.439^2) = 0.176237
    # The off-diagonal is reported untransformed and implies a correlation
    # of -0.788.
    # ------------------------------------------------------------------
    etalcl + etalrbase_shbg ~ c(0.039221,
                                -0.0655, 0.176237)   # Table S2, final comprehensive model: IIV_CL 20.0% (17.9-22.0), IIV_SBL 43.9% (42.0-45.7), COV_IIV -0.0655 (-0.0741 to -0.0569)

    # ------------------------------------------------------------------
    # Residual error. The proportional errors are reported as
    # CV(%) = sqrt(SIGMA^2) * 100, so the tabulated percentage is already a
    # standard deviation. The additive residual-content entries in Table S2
    # are on the VARIANCE scale and are square-rooted here; see the file
    # header note in the vignette Errata -- Table S1 reports the same three
    # quantities on the standard-deviation scale (1.43 / 0.435 / 0.234 mg)
    # and sqrt(1.99 / 0.189 / 0.0541) = 1.411 / 0.435 / 0.233 mg reproduces
    # them, which pins the scale of the Table S2 column.
    # ------------------------------------------------------------------
    propSd            <- 0.235          ; label("Proportional residual error SD on total levonorgestrel concentration (fraction)") # Table S2, final comprehensive model, SIGMA_LNG (CV) = 23.5% (21.6-25.3)
    propSd_shbg       <- 0.229          ; label("Proportional residual error SD on SHBG concentration (fraction)")                 # Table S2, final comprehensive model, SIGMA_SHBG (CV) = 22.9% (21.5-24.3)
    addSd_iusResidual <- sqrt(1.99)     ; label("Additive residual error SD on residual device content (mg)")                      # Table S2, final comprehensive model, SIGMA LNG-IUS 20 res. content = 1.99 (1.28-2.70) as a variance; sqrt = 1.411 mg
  })

  model({
    # Physicochemical constants for the closed-form free fraction of LNG
    # (Table S3a, "Free fraction of LNG").
    MWLNG <- 312.5      # Table S3a: molecular weight LNG, g/mol
    KDS   <- 1.82       # Table S3a: dissociation constant SHBG, nmol/L
    KDA   <- 18209      # Table S3a: dissociation constant albumin, nmol/L
    ALB   <- 700000     # Table S3a: albumin concentration, nmol/L
    KALB  <- 1 + ALB / KDA

    # Reference body weight for both power-function covariate effects.
    wtRef <- 66

    # Release parameters, with the 1e-6 scale factors from the Table S2
    # unit column applied here so that ini() carries the tabulated numbers.
    c12 <- exp(lc12) * 1e-6
    c13 <- exp(lc13) * 1e-6
    t1  <- exp(lt1)

    # Individual disposition parameters.
    vc <- exp(lvc)
    cl <- exp(lcl + e_wt_cl * log(WT / wtRef) + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)

    k20 <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    # SHBG turnover parameters.
    tau        <- exp(ltau)
    ri         <- exp(lri)
    rbase_shbg <- exp(lrbase_shbg + e_wt_rbase_shbg * log(WT / wtRef) + etalrbase_shbg)
    kout_shbg  <- exp(lkout_shbg)
    kin_shbg   <- rbase_shbg * kout_shbg

    # Absolute bioavailability of the reservoir (logit scale).
    fdepot <- exp(logitfdepot) / (1 + exp(logitfdepot))

    # Free fraction of levonorgestrel; see the ivOral model file for the
    # unit derivation of A3nM.
    A3nM  <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    # Delayed linear inhibition of SHBG synthesis by levonorgestrel.
    inh <- ri * effect

    # LNG-IUS 20 release branch of Table S3a.
    input1 <- c12 * depot
    input2 <- c13 * (1 + depot / (t1 + t))

    d/dt(depot)       <- -input1 - input2
    d/dt(central)     <-  input1 + input2 - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * A3nM - (1 / tau) * effect
    d/dt(shbg)        <-  kin_shbg * (1 - inh) - kout_shbg * shbg

    f(depot) <- fdepot
    shbg(0)  <- rbase_shbg

    # Observations. iusResidual is the levonorgestrel remaining in the
    # explanted device (mg); releaseRate is the in-vivo release rate in
    # ug/day, the quantity tabulated in Table 4.
    Cc          <- (central / vc) * 1e6
    CcUnbound   <- fuLNG * Cc
    iusResidual <- depot
    releaseRate <- (input1 + input2) * 1000 * 24

    Cc          ~ prop(propSd)
    shbg        ~ prop(propSd_shbg)
    iusResidual ~ add(addSd_iusResidual)
  })
}
