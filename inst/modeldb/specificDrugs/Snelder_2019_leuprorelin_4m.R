Snelder_2019_leuprorelin_4m <- function() {
  description <- paste(
    "Population PK model for the leuprorelin-SR 4-month 15 mg depot",
    "formulation in patients with advanced prostate cancer (study EC402).",
    "Two-compartment disposition with the peripheral volume fixed to the",
    "central volume, fed by three parallel release depots: a fast",
    "first-order route (Ka1), a slow first-order route with a 160 h lag",
    "time (Ka2), and a second slow route whose dose is delivered as a",
    "zero-order process of duration Dsa and released with a rate constant",
    "that ramps linearly with time after the first dose (Ka3_SLP * t). This",
    "is one of the three leuprorelin PK models in Snelder 2019; it supplies",
    "the drug time-course used to externally validate the EC403",
    "testosterone-PSA model, whose parameters were not re-estimated."
  )
  reference <- paste(
    "Snelder N, Drenth HJ, Riber Bergmann K, Wood ND, Hibberd M, Scott G.",
    "Population pharmacokinetic-pharmacodynamic modelling of the",
    "relationship between testosterone and prostate specific antigen in",
    "patients with prostate cancer during treatment with leuprorelin.",
    "Br J Clin Pharmacol. 2019;85(6):1247-1259. doi:10.1111/bcp.13891.",
    "This model is Supplement 1 (Figure S1.1B and Table S1.2).",
    sep = " "
  )
  vignette <- "Snelder_2019_leuprorelin"

  units <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  # All three release depots receive a dose record carrying the FULL nominal
  # dose; f(depot) / f(depot2) / f(depot3) split it. The depot3 record needs
  # rate = -2 so that dur(depot3) supplies the zero-order duration Dsa.
  # Declared explicitly because auto-detection only sees `depot` / `central`.
  dosing <- c("depot", "depot2", "depot3")

  compartmentData <- list(
    depot       = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    depot3      = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "leuprorelin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "leuprorelin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 106L,
    n_studies      = 1L,
    age_range      = "53-85 years (mean 72.0)",
    age_median     = "mean 72.0 years",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = "Advanced prostate cancer requiring chemical castration",
    dose_range     = "Repeated 15 mg leuprorelin-SR 4-month depot injections",
    co_medication  = paste(
      "Cyproterone acetate with the first leuprorelin dose (300 mg",
      "intramuscular, or 200 mg orally daily for one month) to avoid the",
      "testosterone flare"
    ),
    regions        = "not reported",
    notes          = paste(
      "Study EC402, 4-month 15 mg depot arm (Methods 2.1). EC402 was a",
      "randomised, open-label, comparative study of leuprorelin 1M and",
      "leuprorelin-SR 4M in advanced prostate cancer; 106 patients were",
      "enrolled in this arm and 26 in the 3.75 mg 1M arm (packaged",
      "separately as Snelder_2019_leuprorelin_1m). Inclusion required",
      "biopsy-confirmed prostate cancer, PSA >= 1 ng/mL, serum testosterone",
      ">= 150 ng/dL, no prior LHRH-agonist / antiandrogen / oestrogen",
      "treatment, and life expectancy > 9 months. Across both EC402 arms",
      "10.7% of leuprorelin observations were below the 16 pg/mL LOQ (8.9%",
      "in this arm). EC402 was used only for external validation of the",
      "EC403 testosterone-PSA model (Snelder_2019_leuprorelin), not for its",
      "development."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Leuprorelin PK -- Supplement 1, Table S1.2 (study EC402, 4M 15 mg).
    # Structure from Figure S1.1B: fast first-order depot (Fr1), a slow
    # first-order depot with a lag time ((1-Fr1)*Fr2), and a second slow
    # depot with zero-order dose delivery over Dsa and a linearly
    # time-ramping release rate ((1-Fr1)*(1-Fr2), Ka3_SLP * TIME).
    # ------------------------------------------------------------------
    lcl        <- log(20.5)     ; label("Clearance (L/h)")                                                # Table S1.2: CL = 20.5 (RSE 3.73%)
    lvc        <- log(215)      ; label("Central volume of distribution (L)")                             # Table S1.2: Vc = 215 (RSE 5.16%)
    lq         <- log(4.98)     ; label("Inter-compartmental clearance (L/h)")                            # Table S1.2: Q = 4.98 (RSE 4.44%)
    lka_fast   <- log(0.480)    ; label("Fast-release first-order rate constant Ka1 (1/h)")               # Table S1.2: Ka1 = 0.480 (RSE 7.42%)
    lka_slow   <- log(0.00464)  ; label("Lagged slow-release first-order rate constant Ka2 (1/h)")        # Table S1.2: Ka2 = 0.00464 (RSE 5.37%)
    ltlag      <- log(160)      ; label("Lag time on the lagged slow-release depot (h)")                  # Table S1.2: Lag time = 160 (RSE 1.42%)
    lka_slope  <- log(0.757e-6) ; label("Slope of the linearly time-ramping slow-release rate Ka3_SLP (1/h^2)") # Table S1.2: Ka3_SLP = 0.757e-6 (RSE 9.83%)
    ld1        <- fixed(log(2880)); label("Zero-order duration of dose delivery into the ramping slow-release depot Dsa (h)") # Table S1.2: Dsa = 2880 fixed

    # Table S1.2 reports Fr1 and Fr2 on the logit scale; its footnote spells
    # out the back-transforms and gives the three dose fractions:
    #   fast          = expit(Fr1)                      = 0.31
    #   lagged slow   = (1 - expit(Fr1)) * expit(Fr2)   = 0.11
    #   ramping slow  = (1 - expit(Fr1)) * (1 - expit(Fr2)) = 0.58
    logitfrel  <- -0.799        ; label("Logit of the fast-release dose fraction Fr1 (logit units)")      # Table S1.2: Fr1 = -0.799 (RSE 6.88%)
    logitfrel2 <- -1.62         ; label("Logit of the lagged-slow share of the non-fast dose Fr2 (logit units)") # Table S1.2: Fr2 = -1.62 (RSE 21.7%)

    lfdepot    <- fixed(log(1)) ; label("Relative bioavailability RBIO (unitless)")                       # Table S1.2: RBIO = 1

    # Vp was not estimated: Table S1.2 gives "Vp (L) Fixed to Vc", so vp is
    # set equal to vc inside model() rather than carrying its own parameter.

    # ------------------------------------------------------------------
    # Interindividual variability -- Table S1.2 reports four variances and
    # no covariances, so the block is diagonal. The Fr1 / Fr2 etas are
    # additive on the logit scale; the RBIO and Ka3_SLP etas are additive
    # on the log scale (Methods 2.8: exponential random effects).
    # ------------------------------------------------------------------
    # Table S1.2: omega^2 RBIO = 0.0830 (CV 29.4%, RSE 25.8%)
    etalfdepot    ~ 0.0830
    # Table S1.2: omega^2 Ka3_SLP = 0.264 (CV 55.0%, RSE 33.5%)
    etalka_slope  ~ 0.264
    # Table S1.2: omega^2 Fr1 = 0.0625 (CV 25.4%, RSE 26.6%)
    etalogitfrel  ~ 0.0625
    # Table S1.2: omega^2 Fr2 = 4.24 (CV 827%, RSE 30.9%)
    etalogitfrel2 ~ 4.24

    # ------------------------------------------------------------------
    # Residual error -- proportional; the table reports sigma^2.
    # ------------------------------------------------------------------
    propSd <- 0.3768 ; label("Proportional residual error, leuprorelin (fraction)")  # Table S1.2: sigma^2 prop = 0.142 (CV 39.1%) -> sqrt(0.142) = 0.3768
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters
    # ------------------------------------------------------------------
    fdepot    <- exp(lfdepot + etalfdepot)
    cl        <- exp(lcl)
    vc        <- exp(lvc)
    q         <- exp(lq)
    vp        <- vc                      # Table S1.2: "Vp (L) Fixed to Vc"
    ka_fast   <- exp(lka_fast)
    ka_slow   <- exp(lka_slow)
    ka_slope  <- exp(lka_slope + etalka_slope)
    tlag      <- exp(ltlag)
    d1        <- exp(ld1)
    frel      <- expit(logitfrel  + etalogitfrel)
    frel2     <- expit(logitfrel2 + etalogitfrel2)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Ramping slow-depot release rate constant, Figure S1.1B
    # ("Ka3_SLP*TIME"); `max(0, tafd())` is 0 before the first dose.
    k_slow_ramp <- ka_slope * max(0, tafd())

    # ------------------------------------------------------------------
    # 2. ODE system. depot2 is Figure S1.1B's "slow release dose cmt 2"
    #    (lagged first-order) and depot3 its "slow release dose cmt 3"
    #    (zero-order input, ramping release).
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka_fast    * depot
    d/dt(depot2)      <- -ka_slow    * depot2
    d/dt(depot3)      <- -k_slow_ramp * depot3
    d/dt(central)     <-  ka_fast * depot + ka_slow * depot2 +
      k_slow_ramp * depot3 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)     <- fdepot * frel
    f(depot2)    <- fdepot * (1 - frel) * frel2
    f(depot3)    <- fdepot * (1 - frel) * (1 - frel2)
    alag(depot2) <- tlag
    # Zero-order delivery into depot3 ("Dsa: zero-order absorption" in
    # Figure S1.1B); the depot3 dose record must request a modelled
    # duration (rate = -2).
    dur(depot3)  <- d1

    # ------------------------------------------------------------------
    # 3. Observation
    # ------------------------------------------------------------------
    Cc <- central / vc * 1e6             # mg/L -> pg/mL
    Cc ~ prop(propSd)
  })
}
