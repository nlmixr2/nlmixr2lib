Snelder_2019_leuprorelin_1m <- function() {
  description <- paste(
    "Population PK model for the leuprorelin 1-month 3.75 mg depot",
    "formulation in patients with advanced prostate cancer (study EC402).",
    "Two-compartment disposition with the peripheral volume fixed to the",
    "central volume, fed by two parallel release depots: a fast first-order",
    "route (Ka1) and a slow route whose dose is delivered into the depot as",
    "a zero-order process of duration Dsa and released with a rate constant",
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
    "This model is Supplement 1 (Figure S1.1A and Table S1.1).",
    sep = " "
  )
  vignette <- "Snelder_2019_leuprorelin"

  # Doses in mg and volumes in L give mg/L; the *1e6 in model() converts to
  # pg/mL, the unit of the leuprorelin assay (LOQ 16 pg/mL).
  units <- list(time = "h", dosing = "mg", concentration = "pg/mL")

  # Both release depots receive a dose record carrying the FULL nominal dose;
  # f(depot) and f(depot2) split it. The depot2 record additionally needs
  # rate = -2 so that dur(depot2) supplies the zero-order duration Dsa.
  # Declared explicitly because auto-detection only sees `depot` / `central`.
  dosing <- c("depot", "depot2")

  compartmentData <- list(
    depot       = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "leuprorelin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "leuprorelin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "leuprorelin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    age_range      = "53-85 years (mean 71.1)",
    age_median     = "mean 71.1 years",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state  = "Advanced prostate cancer requiring chemical castration",
    dose_range     = "Repeated 3.75 mg leuprorelin 1-month depot injections",
    co_medication  = paste(
      "Cyproterone acetate with the first leuprorelin dose (300 mg",
      "intramuscular, or 200 mg orally daily for one month) to avoid the",
      "testosterone flare"
    ),
    regions        = "not reported",
    notes          = paste(
      "Study EC402, 1-month 3.75 mg depot arm (Methods 2.1). EC402 was a",
      "randomised, open-label, comparative study of leuprorelin 1M and",
      "leuprorelin-SR 4M in advanced prostate cancer; 26 patients were",
      "enrolled in the 3.75 mg 1M arm and 106 in the 15 mg 4M arm (packaged",
      "separately as Snelder_2019_leuprorelin_4m). Inclusion required",
      "biopsy-confirmed prostate cancer, PSA >= 1 ng/mL, serum testosterone",
      ">= 150 ng/dL, no prior LHRH-agonist / antiandrogen / oestrogen",
      "treatment, and life expectancy > 9 months. Across both EC402 arms",
      "10.7% of leuprorelin observations were below the 16 pg/mL LOQ (12.5%",
      "in this arm). EC402 was used only for external validation of the",
      "EC403 testosterone-PSA model (Snelder_2019_leuprorelin), not for its",
      "development."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Leuprorelin PK -- Supplement 1, Table S1.1 (study EC402, 1M 3.75 mg).
    # Structure from Figure S1.1A: fast first-order depot (Fr) plus a slow
    # depot (1 - Fr) with zero-order dose delivery over Dsa and a linearly
    # time-ramping release rate Ka3_SLP * TIME.
    # ------------------------------------------------------------------
    lcl       <- log(21.7)     ; label("Clearance (L/h)")                                                 # Table S1.1: CL = 21.7 (RSE 5.16%)
    lvc       <- log(219)      ; label("Central volume of distribution (L)")                              # Table S1.1: Vc = 219 (RSE 7.53%)
    lq        <- log(4.71)     ; label("Inter-compartmental clearance (L/h)")                             # Table S1.1: Q = 4.71 (RSE 6.41%)
    lka_fast  <- log(0.690)    ; label("Fast-release first-order rate constant Ka1 (1/h)")                # Table S1.1: Ka1 = 0.690 (RSE 14.8%)
    lka_slope <- log(9.83e-6)  ; label("Slope of the linearly time-ramping slow-release rate Ka3_SLP (1/h^2)") # Table S1.1: Ka3_SLP = 9.83e-6 (RSE 15.0%)
    ld1       <- log(803)      ; label("Zero-order duration of dose delivery into the slow-release depot Dsa (h)") # Table S1.1: Dsa = 803 (RSE 0.493%)

    # Unlike the EC403 and EC402-4M models, Table S1.1 reports Fr on the
    # NATURAL scale: its footnote reads "Dose fraction fast absorption:
    # 0.441 / Dose fraction slow absorption compartment: 0.559", and the
    # 95% CI (0.407-0.475) is a fraction interval. Methods 2.8 states that
    # random effects were exponential, so the IIV is log-normal on Fr.
    lfrel     <- log(0.441)    ; label("Fast-release dose fraction Fr (unitless)")                        # Table S1.1: Fr = 0.441 (RSE 3.95%)

    # RBIO scales both release depots; typical value 1, only its IIV was
    # estimated (Table S1.1: "RBIO 1 fixed").
    lfdepot   <- fixed(log(1)) ; label("Relative bioavailability RBIO (unitless)")                        # Table S1.1: RBIO = 1

    # Vp was not estimated: Table S1.1 gives "Vp (L) Fixed to Vc", so vp is
    # set equal to vc inside model() rather than carrying its own parameter.

    # ------------------------------------------------------------------
    # Interindividual variability -- Table S1.1 reports variances for RBIO,
    # Ka3_SLP and Fr plus two of the three possible covariances. Assembled
    # as a 3x3 block in that order, both printed correlations reproduce
    # exactly (0.543 and -0.651) and the matrix is positive definite.
    # ------------------------------------------------------------------
    # Row-by-row source trace for the block below. The comments must live
    # OUTSIDE the `c(...)`: rxode2 5.1.7 rewrites any comment between `c(`
    # and `)` into a bare `;`, after which the model no longer parses.
    #   row 1  omega^2 RBIO     = 0.0535 (CV 23.4%, RSE 28.6%)
    #   row 2  RBIO x Ka3_SLP   -- not reported by Table S1.1, so set to 0;
    #          omega^2 Ka3_SLP  = 0.420 (CV 72.2%, RSE 40.7%)
    #   row 3  Fr x RBIO        = -0.0228 (r -0.651); Ka3_SLP x Fr = 0.0533 (r 0.543);
    #          omega^2 Fr       = 0.0229 (CV 15.2%, RSE 23.9%)
    etalfdepot + etalka_slope + etalfrel ~ c(
      0.0535,
      0,        0.420,
      -0.0228,  0.0533,  0.0229
    )
    # Supplement 1 states "IIV was identified on RBIO, Ka3_SLP, CL and Fr",
    # but Table S1.1 has no omega^2 CL row, so the CL variance is not
    # recoverable from any on-disk source and no eta on CL is carried here.
    # Do not read the absence as "no IIV on CL" -- see the vignette Errata.

    # ------------------------------------------------------------------
    # Residual error -- proportional; the table reports sigma^2.
    # ------------------------------------------------------------------
    propSd <- 0.2805 ; label("Proportional residual error, leuprorelin (fraction)")  # Table S1.1: sigma^2 prop = 0.0787 (CV 28.6%) -> sqrt(0.0787) = 0.2805
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters
    # ------------------------------------------------------------------
    fdepot   <- exp(lfdepot + etalfdepot)
    cl       <- exp(lcl)
    vc       <- exp(lvc)
    q        <- exp(lq)
    vp       <- vc                       # Table S1.1: "Vp (L) Fixed to Vc"
    ka_fast  <- exp(lka_fast)
    ka_slope <- exp(lka_slope + etalka_slope)
    frel     <- exp(lfrel + etalfrel)
    d1       <- exp(ld1)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Slow-depot release rate constant, Figure S1.1A ("Ka3_SLP*TIME"),
    # written as in the EC403 control stream ($DES: K34 = SLP * TAD1 / 1E6).
    # `max(0, tafd())` is 0 before the first dose and does not propagate NA.
    k_slow <- ka_slope * max(0, tafd())

    # ------------------------------------------------------------------
    # 2. ODE system
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka_fast * depot
    d/dt(depot2)      <- -k_slow  * depot2
    d/dt(central)     <-  ka_fast * depot + k_slow * depot2 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot * frel
    f(depot2)   <- fdepot * (1 - frel)
    # Zero-order delivery of the slow-release dose fraction into depot2
    # ("Dsa: zero-order absorption" in Figure S1.1A). The dose record for
    # depot2 must request a modelled duration (rate = -2).
    dur(depot2) <- d1

    # ------------------------------------------------------------------
    # 3. Observation
    # ------------------------------------------------------------------
    Cc <- central / vc * 1e6             # mg/L -> pg/mL
    Cc ~ prop(propSd)
  })
}
