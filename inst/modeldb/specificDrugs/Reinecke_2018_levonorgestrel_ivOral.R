Reinecke_2018_levonorgestrel_ivOral <- function() {
  description <- paste(
    "Final intravenous-oral model of the Reinecke 2018 integrated",
    "levonorgestrel (LNG) population PK analysis: two-compartment LNG",
    "disposition coupled to a sex hormone-binding globulin (SHBG) turnover",
    "model. Only unbound LNG is eliminated from, or distributed out of, the",
    "central compartment, so elimination is cl * fuLNG * central / vc and",
    "distribution is q * fuLNG * central / vc. The free fraction fuLNG is the",
    "closed-form solution of reversible LNG binding to SHBG (KDS = 1.82",
    "nmol/L) and to albumin (KDA = 18209 nmol/L at a constant ALB = 700000",
    "nmol/L). SHBG is an indirect-response turnover state (kin = rbase_shbg *",
    "kout_shbg) whose synthesis is linearly inhibited by a delay-compartment-",
    "smoothed molar LNG signal (delay time constant tau). Oral bioavailability",
    "decreases with body weight as a power function centred on 66 kg. This is",
    "the model the paper fitted to intravenous, progestin-only-pill and",
    "supportive single-dose oral data, and the model whose parameters the",
    "combined-oral-contraceptive model inherits."
  )
  reference <- paste(
    "Reinecke I, Hofmann B, Mesic E, Drenth HJ, Garmann D.",
    "An integrated population pharmacokinetic analysis to characterize",
    "levonorgestrel pharmacokinetics after different administration routes.",
    "J Clin Pharmacol. 2018 Dec;58(12):1639-1654. doi:10.1002/jcph.1288.",
    "Parameter values from Supplemental Table S2, 'final iv-oral model'",
    "block; model equations from Supplemental Table S3a."
  )
  vignette <- "Reinecke_2018_levonorgestrel_contraceptives"
  units <- list(
    time          = "h",
    dosing        = "mg (levonorgestrel)",
    concentration = "ng/L (total and unbound LNG in serum/plasma); SHBG in nmol/L"
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
      description        = "Body weight (kg). Power-function ('allometric') covariate on the absolute oral bioavailability of levonorgestrel; bioavailability decreases with increasing body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight, centred on 66 kg. The reference weight is not",
        "stated in words but is pinned by Supplemental Table S2 itself: the",
        "row labelled 'Absolute oral bioavailability (median[WGHT]=66kg)'",
        "gives 0.770, and the derived values in the same row block (0.825 at",
        "62 kg, 0.819 at 62.4 kg, 0.840 at 61 kg) are reproduced to three",
        "decimal places by 0.770 * (WT/66)^-1.10. Body weight and age were",
        "both screened on clearance and on the SHBG baseline in this model",
        "and neither was significant there; only the effect on oral",
        "bioavailability was retained (Results, 'Intravenous-Oral Model')."
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
        "Not retained: 'The covariates body weight and age did not have a",
        "significant effect on CL or SBL' (Results, 'Intravenous-Oral",
        "Model'). No point estimate is reported, so the effect cannot be",
        "carried into the model."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 3L,
    studies        = c(
      "Phase 1 study 92085 (single-dose intravenous 90 ug LNG and single-dose oral 30 / 90 / 270 ug LNG, crossover in the same 18 women)",
      "Phase 1 study 15687 (progestin-only pill Microlut/Norgeston, 30 ug LNG once daily, 1 treatment cycle, 21 women)",
      "Phase 1 study 90078 (oral 150 ug LNG once daily, 1 treatment cycle, 12 women)"
    ),
    age_range      = "21-45 years (median 35 years across the pooled oral data, N = 51)",
    weight_range   = "51.2-82.9 kg (median 62 kg across the pooled oral data, N = 51)",
    sex_female_pct = 100,
    disease_state  = "Healthy premenopausal women",
    dose_range     = "90 ug LNG single intravenous dose; 30, 90, 150 and 270 ug LNG oral",
    notes          = paste(
      "Demographics from Table 1 and its footnote b. The intravenous and",
      "supportive oral data come from the same 18 subjects in study 92085.",
      "Observations contributing to this step: 286 (IV) + 646 (supportive",
      "oral) + 599 (POP) + 477 (study 90078) LNG concentrations and 621",
      "(POP) + 226 (study 90078) SHBG concentrations. LNG was assayed by",
      "radio-immunoassay except in study 15687 (LC-MS/MS); SHBG by radio-",
      "immunoassay or time-resolved fluoroimmunoassay (LLOQ 9.8 or 10",
      "nmol/L). Both eta and epsilon shrinkage were < 2.5%."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Levonorgestrel disposition. Note that cl and q are clearances of
    # UNBOUND drug: the ODEs multiply them by the free fraction fuLNG, so
    # the total-concentration-based clearance is cl * fuLNG, about
    # 0.015 * 279 = 4.2 L/h (Discussion).
    # ------------------------------------------------------------------
    lka  <- log(2.03)   ; label("Oral first-order absorption rate constant Ka (1/h)")                        # Table S2, final iv-oral model, Ka = 2.03 (1.62-2.44)
    lvc  <- log(22.0)   ; label("Central volume of distribution Vc (L)")                                     # Table S2, final iv-oral model, Vc = 22.0 (19.1-24.9)
    lcl  <- log(279)    ; label("Clearance of unbound levonorgestrel CL (L/h)")                              # Table S2, final iv-oral model, CL = 279 (249-309)
    lvp  <- log(5120)   ; label("Peripheral volume of distribution Vp (L)")                                  # Table S2, final iv-oral model, Vp = 5120 (4040-6200)
    lq   <- log(649)    ; label("Intercompartmental clearance of unbound levonorgestrel Q (L/h)")            # Table S2, final iv-oral model, Q = 649 (556-742)

    # ------------------------------------------------------------------
    # SHBG turnover and the delayed inhibition of SHBG synthesis by LNG.
    # kin is derived inside model() as rbase_shbg * kout_shbg so that the
    # drug-free steady state of the SHBG state equals rbase_shbg.
    # ------------------------------------------------------------------
    lrbase_shbg <- log(70.7)     ; label("Baseline SHBG serum concentration SBL for iv/oral treatment (nmol/L)")   # Table S2, final iv-oral model, SBL_iv,oral = 70.7 (61.7-79.7)
    lkout_shbg  <- log(0.00313)  ; label("First-order elimination rate constant of SHBG KOUT (1/h)")               # Table S2, final iv-oral model, KOUT = 0.00313 (0.00107-0.00519)
    ltau        <- log(13.7)     ; label("Delay time constant tau of the levonorgestrel effect on SHBG (h)")       # Table S2, final iv-oral model, tau = 13.7 (3.17-24.2)
    lri         <- log(0.232)    ; label("Linear inhibition factor RI of SHBG synthesis by delayed LNG (L/nmol)")  # Table S2, final iv-oral model, R_I = 0.232 (0.109-0.355)

    # ------------------------------------------------------------------
    # Absolute oral bioavailability at the reference weight of 66 kg, and
    # its body-weight power exponent.
    # ------------------------------------------------------------------
    lfdepot      <- log(0.770)   ; label("Absolute oral bioavailability F_oral at 66 kg (fraction)")               # Table S2, final iv-oral model, F_oral = 0.770 (0.703-0.837) at median[WGHT] = 66 kg
    e_wt_fdepot  <- -1.10        ; label("Body-weight power exponent on oral bioavailability, (WT/66)^e_wt_fdepot") # Table S2, final iv-oral model, C_WGHT,Foral = -1.10 (-1.54 to -0.655)

    # ------------------------------------------------------------------
    # Inter-individual variability on clearance and on the SHBG baseline,
    # with covariance between the two etas ("included as exponential term
    # with covariance between the 2 terms", Results).
    #
    # Table S2 reports the diagonal entries as coefficients of variation
    # in %, with footnote d giving the transformation the authors used:
    # CV(%) = sqrt(exp(OMEGA^2) - 1) * 100. Inverting it:
    #   CL  : log(1 + 0.258^2) = 0.064442
    #   SBL : log(1 + 0.511^2) = 0.232001
    # The off-diagonal covariance is reported untransformed (-0.0575),
    # which implies a correlation of -0.470.
    # ------------------------------------------------------------------
    etalcl + etalrbase_shbg ~ c(0.064442,
                                -0.0575, 0.232001)   # Table S2, final iv-oral model: IIV_CL 25.8% (17.4-32.3), IIV_SBL 51.1% (26.7-69.6), COV_IIV -0.0575 (-0.110 to -0.00517)

    # ------------------------------------------------------------------
    # Residual error. Both LNG and SHBG carry proportional error. Table S2
    # reports these as CV(%) = sqrt(SIGMA^2) * 100 (footnote d), so the
    # tabulated percentage IS the standard deviation on the proportional
    # scale and is used here directly.
    # ------------------------------------------------------------------
    propSd      <- 0.286   ; label("Proportional residual error SD on total levonorgestrel concentration (fraction)")  # Table S2, final iv-oral model, SIGMA_LNG (CV) = 28.6% (24.2-32.4)
    propSd_shbg <- 0.204   ; label("Proportional residual error SD on SHBG concentration (fraction)")                  # Table S2, final iv-oral model, SIGMA_SHBG (CV) = 20.4% (15.9-24.1)
  })

  model({
    # ------------------------------------------------------------------
    # Physicochemical constants for the closed-form free fraction of LNG
    # (Supplemental Table S3a, "Free fraction of LNG"). KDS and KDA are
    # annotated in the source as unpublished data; ALB is fixed to "an
    # appropriate physiological value" (Discussion: physiological range
    # 530000-830000 nmol/L, i.e. 3.5-5.5 g/dL).
    # ------------------------------------------------------------------
    MWLNG <- 312.5      # Table S3a: molecular weight LNG, g/mol
    KDS   <- 1.82       # Table S3a: dissociation constant SHBG, nmol/L
    KDA   <- 18209      # Table S3a: dissociation constant albumin, nmol/L
    ALB   <- 700000     # Table S3a: albumin concentration, nmol/L
    KALB  <- 1 + ALB / KDA

    # Reference body weight for the oral-bioavailability covariate.
    wtRef <- 66

    # Individual disposition parameters.
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl)
    vp <- exp(lvp)
    q  <- exp(lq)

    k20 <- cl / vc
    k23 <- q / vc
    k32 <- q / vp

    # SHBG turnover parameters.
    tau        <- exp(ltau)
    ri         <- exp(lri)
    rbase_shbg <- exp(lrbase_shbg + etalrbase_shbg)
    kout_shbg  <- exp(lkout_shbg)
    kin_shbg   <- rbase_shbg * kout_shbg

    # Absolute oral bioavailability with the body-weight power function.
    fdepot <- exp(lfdepot + e_wt_fdepot * log(WT / wtRef))

    # ------------------------------------------------------------------
    # Free fraction of levonorgestrel. A3nM is the molar total LNG
    # concentration in the central compartment: central is in mg and vc in
    # L, so (central / vc) is mg/L and the 1e6 / MWLNG factor converts it
    # to nmol/L exactly as written in Table S3a. The 1e-6 guard keeps
    # A3nM strictly positive so the 1 / A3nM factor is finite before the
    # first dose reaches the central compartment.
    # ------------------------------------------------------------------
    A3nM  <- (central / vc) * (1e6 / MWLNG) + 1e-6
    temp1 <- KALB * KDS + shbg - A3nM
    temp2 <- temp1 / (KALB * 2)
    temp3 <- temp2^2 + (KDS / KALB) * A3nM
    fuLNG <- (1 / A3nM) * (-temp2 + sqrt(temp3))

    # Delayed linear inhibition of SHBG synthesis by levonorgestrel.
    inh <- ri * effect

    # ------------------------------------------------------------------
    # ODE system (Table S3a). Only unbound LNG leaves the central
    # compartment, hence the fuLNG factor on both k20 and k23.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (k20 + k23) * fuLNG * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * fuLNG * central - k32 * peripheral1
    d/dt(effect)      <- (1 / tau) * A3nM - (1 / tau) * effect
    d/dt(shbg)        <-  kin_shbg * (1 - inh) - kout_shbg * shbg

    f(depot) <- fdepot
    shbg(0)  <- rbase_shbg

    # ------------------------------------------------------------------
    # Observations. Cc is the total LNG serum/plasma concentration in
    # ng/L: (central [mg] / vc [L]) * 1e6 ng/mg. CcUnbound is reported by
    # the paper (Tables S4, S6, S8) and is derived here without a separate
    # residual error, because the paper attaches residual error to the
    # measured total LNG and SHBG concentrations only.
    # ------------------------------------------------------------------
    Cc        <- (central / vc) * 1e6
    CcUnbound <- fuLNG * Cc

    Cc   ~ prop(propSd)
    shbg ~ prop(propSd_shbg)
  })
}
