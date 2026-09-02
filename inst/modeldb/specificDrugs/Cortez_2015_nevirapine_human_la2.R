Cortez_2015_nevirapine_human_la2 <- function() {
  description <- paste(
    "Human infant projection of the NVP LA-2 long-acting injectable",
    "nevirapine microparticle suspension for prevention of mother-to-child",
    "HIV transmission in breastfeeding neonates. One-compartment model with",
    "first-order absorption in which the SYSTEMIC parameters (V/F and CL/F)",
    "are the published breastfeeding-infant oral nevirapine population",
    "estimates from two Ugandan studies and the PRESYSTEMIC absorption rate",
    "constant is carried over from the authors' own Sprague-Dawley rat fit",
    "(Cortez_2015_nevirapine_rat_la2), allometrically scaled from the 0.21 kg",
    "rat reference weight with a fixed -0.25 exponent. Deterministic",
    "naive-pooled simulation model; the publication reports neither",
    "between-subject variability nor a residual-error model for the human",
    "projection. IMPORTANT: this encodes the parameters the Methods state.",
    "It reproduces the t = 0 intercept of Cortez 2015 Figure 4C and 4D",
    "(60 mg / 20 L = 3.0 ug/mL) but NOT their decay: those panels fall with",
    "an apparent half-life near 66 days, whereas V/F = 20 L with",
    "CL/F = 0.21 L/h gives the 66 HOURS the same Methods paragraph reports.",
    "See the vignette Assumptions and deviations.",
    sep = " ")
  reference <- paste(
    "Cortez JM Jr, Quintero R, Moss JA, Beliveau M, Smith TJ, Baum MM.",
    "Pharmacokinetics of injectable, long-acting nevirapine for HIV",
    "prophylaxis in breastfeeding infants.",
    "Antimicrob Agents Chemother. 2015;59(1):59-66.",
    "doi:10.1128/AAC.03906-14.",
    "Presystemic absorption carried from the companion rat fit; see",
    "modellib('Cortez_2015_nevirapine_rat_la2').",
    sep = " ")
  vignette <- "Cortez_2015_nevirapine"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot   = list(analyte = "nevirapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nevirapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size covariate on the RAT-DERIVED absorption rate",
        "constant only, referenced to the 0.21 kg mean weight of the",
        "study-2 rats. The human systemic parameters V/F and CL/F are the",
        "published infant population estimates and are NOT weight-scaled,",
        "matching the composition recovered from Cortez 2015 Figure 4A and",
        "4B for the companion LA-1 projection. Simulated infant weights",
        "were 2.0 and 3.9 kg."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "Neonates from birth (WHO prophylaxis window is birth to 4-6 weeks of age)",
    weight_range   = "2.0-3.9 kg",
    disease_state  = paste(
      "HIV-uninfected breastfeeding infants born to HIV-1-infected mothers,",
      "receiving nevirapine as prophylaxis against mother-to-child",
      "transmission."
    ),
    dose_range     = "Single 60 mg subcutaneous injection of NVP LA-2 (Cortez 2015 Figure 4C and 4D)",
    regions        = "Uganda (source of the systemic infant PK estimates)",
    notes          = paste(
      "Cortez 2015 Methods 'Human simulations': the infant systemic",
      "parameters were V/bioavailability (F) = 20.0 liters, CL/F = 0.21",
      "liters h-1 and t1/2 = 66 h, over a newborn weight range of 2.0 to",
      "3.9 kg, and \"These systemic parameters were combined with",
      "presystemic rat parameters, essentially the absorption rate",
      "constant, ka, in order to simulate human exposure.\" The",
      "bioavailable fraction was assumed to be 100% (Cortez 2015 Results",
      "'Human simulation'). The prophylactic target window motivating the",
      "simulation is 0.2 to 3.0 ug/mL (Cortez 2015 Discussion",
      "'Implications for HIV prevention in breastfeeding infants')."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # SYSTEMIC parameters -- human infant population estimates, taken
    # as-is from Cortez 2015 Methods 'Human simulations', and NOT
    # weight-scaled (the composition established by reproducing Figure 4A
    # and 4B in the companion LA-1 projection).
    #
    # The three quantities the paper reports are mutually consistent:
    #   log(2) * 20 L / 0.21 L/h = 66.0 h, the reported t1/2.
    # ------------------------------------------------------------------
    lvc <- log(20)   ; label("Apparent central volume V/F (L), human infant")   # Cortez 2015 Methods 'Human simulations': V/bioavailability (F), 20.0 liters
    lcl <- log(5.04) ; label("Apparent clearance CL/F (L/day), human infant")   # Cortez 2015 Methods 'Human simulations': CL/F = 0.21 L/h; 0.21 * 24 = 5.04 L/day in this model's day time base

    # ------------------------------------------------------------------
    # PRESYSTEMIC parameter -- the rate-limiting release of drug from the
    # microparticles, carried from the study-2 rat fit and allometrically
    # scaled in model() from the 0.21 kg rat reference weight. This is the
    # one quantity the rat study is meant to supply.
    # ------------------------------------------------------------------
    lka <- log(15.4) ; label("Absorption rate constant ka (1/day) at the 0.21 kg rat reference weight")  # Cortez 2015 Table 2 study 2: Ka = 15.4 /day

    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (fraction)")  # Cortez 2015 Results 'Human simulation': "The bioavailable fraction was assumed to be 100% in these simulations."

    # ------------------------------------------------------------------
    # Allometric exponent -- Cortez 2015 Methods 'Structural PK model
    # buildup': "The nominal covariate effect was -0.25 on absorption
    # rate, 1 on volumes, and 0.75 on clearance (19)." Imposed, not
    # estimated.
    # ------------------------------------------------------------------
    e_wt_ka <- fixed(-0.25) ; label("Allometric exponent on ka (unitless)")  # Cortez 2015 Methods 'Structural PK model buildup'

    # ------------------------------------------------------------------
    # No IIV and no residual-error model: Cortez 2015 simulated mean
    # plasma concentrations by a naive-pooled approach and reports
    # neither an omega nor a sigma for the human projection.
    # ------------------------------------------------------------------
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters. Human systemic parameters at their
    #    published values; rat-derived ka allometrically scaled.
    # ------------------------------------------------------------------
    vc <- exp(lvc)
    cl <- exp(lcl)
    ka <- exp(lka) * (WT / 0.21)^e_wt_ka

    # ------------------------------------------------------------------
    # 2. Micro-constant.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system.
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ------------------------------------------------------------------
    # 4. Bioavailability (assumed 100%).
    # ------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 5. Observation. Deterministic simulation; no residual error.
    # ------------------------------------------------------------------
    Cc <- central / vc
  })
}
