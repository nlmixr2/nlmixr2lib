GonzalezSales_2015_testosterone <- function() {
  description <- paste(
    "Stretched-cosine model of the endogenous (baseline) circadian rhythm of",
    "serum testosterone in adult hypogonadal men (Gonzalez-Sales 2015 AAPS J).",
    "T(t) = Base + Amplitude * cos(pi * f(t; tacro, tnadir)) where f is the",
    "piecewise-linear phase function that is 0 at the peak time tacro, 1 at",
    "the nadir time tnadir, and 2 at the next peak; the descending arm lasts",
    "L1 = (tnadir - tacro) mod 24 hours and the ascending arm lasts",
    "L2 = 24 - L1 hours, allowing an asymmetric (stretched) cycle. Baseline",
    "is reduced by 2.40% per decade with age (centred on the pooled median",
    "49.9 y) and elevated by 8.09% during winter and spring (SEMESTER = 1)",
    "relative to summer and fall (SEMESTER = 0). IIV on Base enters via a",
    "Box-Cox-transformed normal eta (lambda = -1.93; Petersson 2009 form).",
    "Pure typical-value + IIV simulation model (no exogenous drug); time is",
    "clock time in hours after midnight."
  )
  reference <- paste(
    "Gonzalez-Sales M, Barriere O, Tremblay PO, Nekka F, Desrochers J,",
    "Tanguay M. Modeling Testosterone Circadian Rhythm in Hypogonadal",
    "Males: Effect of Age and Circannual Variations.",
    "AAPS J. 2016 Jan;18(1):217-227.",
    "doi:10.1208/s12248-015-9841-6 (published online 2015 Nov 9)."
  )
  vignette <- "GonzalezSales_2015_testosterone"
  units <- list(
    time          = "hour (clock time; t = 0 corresponds to midnight)",
    dosing        = "n/a (endogenous baseline model; no exogenous drug)",
    concentration = "ng/dL (serum testosterone)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age (years). Linear-deviation effect on Base centred on the pooled median 49.9 y; -2.40% per 10-year age increase.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Reference 49.9 y is the overall median in Gonzalez-Sales 2015 Table I (range 21-76 y). The continuous covariate enters as `(AGE - 49.9) / 10` so the multiplier on Base reproduces the paper's '-2.40 %/10 years' linear effect.",
      source_name        = "AGE"
    ),
    SEMESTER = list(
      description        = "Paired-season indicator: 1 = winter or spring (December 21 to June 19), 0 = summer or fall (June 20 to December 20). Multiplicative +8.09% effect on Base.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (summer or fall)",
      notes              = "Per Gonzalez-Sales 2015 Methods (Covariate Analysis): 'semester 1 was defined as a dichotomous variable with a value of 0, if summer or fall, and a value of 1, if winter or spring.' Time-fixed per simulation: downstream users would set this based on the calendar date at which the testosterone profile is simulated. Source-paper variable name is `SEMESTER1`; this register uses the unsuffixed `SEMESTER` as the canonical column.",
      source_name        = "SEMESTER1"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 859L,
    n_studies      = 7L,
    age_range      = "21-76 years (overall median 49.9)",
    weight_range   = "60.0-131 kg (overall median 87.0)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 90.6, Black = 5.8, Asian = 1.2, Other = 2.4),
    disease_state  = "Adult hypogonadal men (mean serum testosterone < 300 ng/dL; individual morning serum testosterone <= 350 ng/dL; BMI 18-37 kg/m^2). No other medical conditions.",
    dose_range     = "n/a (endogenous baseline only; no testosterone replacement therapy was administered before the analysed samples)",
    regions        = "Canada (Quebec, Montreal, Toronto), United States (North Carolina, Florida, San Antonio TX), Germany",
    notes          = "Pooled baseline / pre-dose profiles from 7 internal hypogonadism trials at inVentiv Health (4556 testosterone observations). Sampling spanned the 24-hour clock at 2-hour intervals (Table II). 22.8% Hispanic or Latino. Race was tested as a covariate but excluded from selection because of imbalance (90.6% White). Estimation in NONMEM 7.3 via SAEM (with IMP for standard errors); 1000 non-parametric bootstrap replicates all converged."
  )

  ini({
    # Structural typical values -- Gonzalez-Sales 2015 Table III (original-dataset
    # column, fitted parameter values; RSE in parentheses).
    lrbase  <- log(239)            ; label("Baseline (mesor) testosterone Base (ng/dL)")                # Table III: Base = 239 ng/dL; RSE 1.1%
    lamp    <- log(32.1)           ; label("Amplitude of the stretched-cosine oscillation Amplitude (ng/dL)")  # Table III: Amplitude = 32.1 ng/dL; RSE 3.5%
    ltacro  <- log(9 + 22 / 60)    ; label("Clock time of the testosterone peak t_max (h after midnight)")    # Table III: t_max = 9:22 = 9.36667 h; RSE 0.4%
    ltnadir <- log(14 + 2 / 60)    ; label("Clock time of the testosterone nadir t_min (h after midnight)")   # Table III: t_min = 14:02 = 14.03333 h; RSE 0.4%

    # Covariate effects on Base. Encoded so that AGE = 49.9 y (median) and
    # SEMESTER = 0 (summer or fall) reproduce the typical value exp(lrbase).
    e_age_rbase      <- -0.0240    ; label("Linear-deviation effect of AGE on Base per (AGE - 49.9 y) / 10 y; -2.40% per decade")  # Table III: Age on Base = -2.40 %/10 years; RSE 24.0%
    e_semester_rbase <-  0.0809    ; label("Multiplicative effect of SEMESTER (winter or spring) on Base (+8.09%)")              # Table III: Season on Base = 8.09%; RSE 18.8%

    # Box-Cox shape parameter (lambda) for the IIV on Base. Petersson 2009
    # form: eta_BC = (exp(lambda * eta) - 1) / lambda. lambda = -1.93 yields
    # the negative-skewness shape seen in the empirical eta_Base distribution
    # (Gonzalez-Sales 2015 Results: "the additional parameter estimated
    # required for the Box-Cox transformation was close to a value of -2,
    # indicating negative skewness"; numerical value from Table III).
    boxcox_rbase     <- -1.93      ; label("Box-Cox shape parameter lambda for the IIV transformation on Base (dimensionless)")  # Table III: Box-Cox on Base = -1.93; RSE 2.4%

    # IIV variances. Gonzalez-Sales 2015 Methods (Statistical Model): "The
    # magnitude of IIV was expressed as coefficient of variation (CV)
    # calculated as the square root of omega_P^2." That is, CV (as reported
    # in Table III) is the standard deviation of the underlying normally-
    # distributed eta itself, so omega^2 = (CV / 100)^2. This convention
    # differs from the lognormal sqrt(exp(omega^2) - 1) form used in many
    # other popPK papers and is documented in the validation vignette
    # Assumptions section. The paper estimated correlated IIV blocks for
    # (tmax, tmin) and (Amplitude, Base) (Results paragraph 4) but did not
    # publish the off-diagonal covariance values; the packaged model uses a
    # diagonal Omega and surfaces this deviation in the vignette Errata.
    etalrbase  ~ 0.0630    # Table III: eta_Base CV = 25.1%; var = 0.251^2 = 0.06300. Eta enters Base via the Box-Cox transformation above.
    etalamp    ~ 0.2530    # Table III: eta_Amplitude CV = 50.3%; var = 0.503^2 = 0.25301
    etaltacro  ~ 0.01124   # Table III: eta_tmax CV = 10.6%; var = 0.106^2 = 0.011236
    etaltnadir ~ 0.0361    # Table III: eta_tmin CV = 19.0%; var = 0.190^2 = 0.03610

    # Residual error -- proportional only (Gonzalez-Sales 2015 Eq. 4 and
    # Table III). Reported as 13.8% CV; CV equals sigma here per the same
    # Methods convention used for the IIV (sqrt of variance, not the
    # lognormal CV form).
    propSd <- 0.138                ; label("Proportional residual SD on testosterone observation (fraction)")  # Table III: sigma proportional = 13.8% CV; RSE 2.4%
  })

  model({
    # Box-Cox-transformed eta on Base. Petersson 2009 (Eq. 10 of
    # Gonzalez-Sales 2015): eta_BC = (exp(lambda * eta) - 1) / lambda. lambda
    # is `boxcox_rbase` (-1.93 here); the transformed eta is added on the log
    # scale to lrbase, equivalent to multiplying Base by exp(eta_BC) on the
    # natural scale.
    eta_rbase_bc <- (exp(boxcox_rbase * etalrbase) - 1) / boxcox_rbase

    # Covariate-adjusted typical baseline. AGE enters as a linear deviation
    # centred on the pooled median 49.9 y, scaled per decade. SEMESTER enters
    # multiplicatively (no centring; reference is SEMESTER = 0).
    rbase_cov <- exp(lrbase) *
      (1 + e_age_rbase * (AGE - 49.9) / 10) *
      (1 + e_semester_rbase * SEMESTER)

    # Individual structural parameters. tmax / tmin / Amplitude carry
    # standard log-normal IIV; Base carries the Box-Cox-transformed eta.
    rbase_i  <- rbase_cov * exp(eta_rbase_bc)
    amp_i    <- exp(lamp + etalamp)
    tacro_i  <- exp(ltacro + etaltacro)
    tnadir_i <- exp(ltnadir + etaltnadir)

    # Stretched-cosine phase function f(t) (Gonzalez-Sales 2015 Eqs. 5-9):
    #   t24    = (t - tacro_i) mod 24       -- time since the most recent peak (hours)
    #   L1     = (tnadir_i - tacro_i) mod 24 -- length of the descending arm
    #   L2     = 24 - L1                     -- length of the ascending arm
    #   fphase = t24 / L1                    if t24 <= L1   (descent: 0 -> 1)
    #          = 1 + (t24 - L1) / L2         otherwise      (ascent: 1 -> 2)
    # The observation is then Cc = Base + Amplitude * cos(pi * fphase) so
    # that the peak (cos = +1) occurs at fphase = 0 and the nadir (cos = -1)
    # at fphase = 1.
    t24 <- time - tacro_i - 24 * floor((time - tacro_i) / 24)
    L1  <- tnadir_i - tacro_i - 24 * floor((tnadir_i - tacro_i) / 24)
    L2  <- 24 - L1

    # Defensive epsilons in the denominators: under nominal parameter values
    # L1 ~ 4.67 h and L2 ~ 19.33 h, so the epsilons are numerically
    # irrelevant; they only guard against pathological eta draws.
    phase_desc <- t24 / (L1 + 1e-9)
    phase_asc  <- 1 + (t24 - L1) / (L2 + 1e-9)
    fphase     <- (t24 <= L1) * phase_desc + (t24 > L1) * phase_asc

    # Observation: serum testosterone concentration (ng/dL).
    Cc <- rbase_i + amp_i * cos(pi * fphase)

    Cc ~ prop(propSd)
  })
}
