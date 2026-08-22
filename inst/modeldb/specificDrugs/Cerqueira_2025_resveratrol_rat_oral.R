Cerqueira_2025_resveratrol_rat_oral <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model with first-order absorption for",
    "resveratrol (3,5,4'-trihydroxy-trans-stilbene) after a single",
    "100 mg/kg oral gavage dose to male Wistar rats. Fitted in Monolix",
    "2020R1 (SAEM) and parameterised in micro-constant form exactly as the",
    "authors report it: first-order absorption ka from the gavage depot, a",
    "single central volume V (L/kg), first-order elimination k10 from",
    "central, and inter-compartmental rate constants k12 and k21. No",
    "bioavailability term was fitted, so the full 100 mg/kg dose enters",
    "the depot and V is an apparent volume V/F; the rate constants remain",
    "true rate constants. Non-compartmental analysis in the same paper put",
    "absolute oral bioavailability at roughly 6%. Because every dose and",
    "volume term in the source is normalised to body weight, the model is",
    "coded on a per-kilogram basis: the dosed amount is ug/kg and the",
    "volume is L/kg, so central/vc lands directly in ng/mL, the units the",
    "assay and Figures 1b and 2b report. The oral arm of the paper; see",
    "Cerqueira_2025_resveratrol_rat_iv for the separately fitted 5 mg/kg",
    "intravenous arm.",
    sep = " ")
  reference <- paste(
    "Cerqueira C, Santos V, Araujo J, Pereira L, Batista F, Soares D,",
    "Azeredo F, Ferreira E.",
    "Development of a Population Pharmacokinetic Model Characterizing the",
    "Tissue Distribution of Resveratrol After Administration by Different",
    "Routes and Doses in Rats.",
    "Nutrients. 2025;17(1):181.",
    "doi:10.3390/nu17010181.",
    sep = " ")
  vignette <- "Cerqueira_2025_resveratrol"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  # Per-kg normalisation: the source reports V in L/kg and dose in mg/kg,
  # so the state amounts are ug/kg and central/vc is ug/L == ng/mL.
  compartmentData <- list(
    depot       = list(analyte = "resveratrol", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "resveratrol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "resveratrol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "rat (Wistar, male)",
    n_subjects     = 6L,
    n_studies      = 1L,
    age_range      = "2-3 months",
    weight_range   = "250-300 g",
    sex_female_pct = 0,
    disease_state  = "Healthy (no disease model)",
    dose_range     = paste(
      "Single 100 mg/kg dose by oral gavage. Resveratrol (purity > 98%)",
      "formulated as a 3 mg/mL solution in 10% DMSO made up to volume",
      "with saline."
    ),
    regions        = "Brazil (Federal University of Bahia, Salvador)",
    notes          = paste(
      "Serial caudal-vein sampling at 0.25, 0.5, 0.75, 1, 1.5, 2, 4, 6, 10",
      "and 24 h post-dose; resveratrol quantified by HPLC-UV at 310 nm",
      "(LLOQ 62.5 ng/mL, calibration range 62.5-5000 ng/mL). Animals were",
      "kept on a 12 h light-dark cycle with ad libitum food and water;",
      "ethics protocol 43/2016, Federal University of Bahia. The mean oral",
      "profile in Figure 1b shows two concentration peaks, at 2 and 6 h,",
      "which the authors attribute to enterohepatic recirculation",
      "(Results 3.2, Discussion); the fitted model has a single first-order",
      "absorption process and therefore does NOT reproduce the second peak.",
      "Plasma protein binding measured separately by ultrafiltration was",
      "79 +/- 5% and linear over 0.5-50 ug/mL (Results 3.4); it is not a",
      "model parameter. See Cerqueira 2025 Methods 2.2, 2.5 and Table 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Cerqueira 2025 Table 1, "Mean Oral (%RSE)"
    # column. Monolix micro-constant parameterisation (ka, V, k10, k12,
    # k21). The paper's k10 is the canonical kel. V is apparent (V/F): no
    # bioavailability term appears in Table 1 and the full 100 mg/kg dose
    # was used in the fit.
    # ------------------------------------------------------------------
    lka  <- log(1.88)  ; label("First-order absorption rate constant ka (log 1/h)")            # Table 1, Mean Oral: ka = 1.88 /h (RSE 25.8%)
    lvc  <- log(22.20) ; label("Apparent central volume of distribution V/F (log L/kg)")       # Table 1, Mean Oral: V = 22.20 L/kg (RSE 29.5%)
    lkel <- log(0.43)  ; label("Elimination rate constant from central k10 (log 1/h)")         # Table 1, Mean Oral: k10 = 0.43 /h (RSE 21.7%)
    lk12 <- log(4.09)  ; label("Rate constant central -> peripheral1 k12 (log 1/h)")           # Table 1, Mean Oral: k12 = 4.09 /h (RSE 23.1%)
    lk21 <- log(0.44)  ; label("Rate constant peripheral1 -> central k21 (log 1/h)")           # Table 1, Mean Oral: k21 = 0.44 /h (RSE 23.1%)

    # ------------------------------------------------------------------
    # IIV -- Cerqueira 2025 Table 1, "BSV Oral (%RSE)" column. Methods 2.5
    # declares log-normal IIV, theta_i = theta * exp(eta), var(eta) =
    # omega^2.
    #
    # SCALE OF THE BSV COLUMN. The paper does not state the scale of the
    # tabulated BSV values. They are read here as the SD of the individual
    # parameter on the NATURAL (untransformed) scale, i.e. an approximate
    # CV of BSV/typical, and converted with omega^2 = log(1 + CV^2).
    # The alternative reading -- BSV = omega, the SD of eta on the log
    # scale -- is untenable, and this arm is what falsifies it: V has
    # typical 22.20 and BSV 5.60, which as an omega is a CV of roughly
    # 5e8 % and cannot be simulated. Under the natural-scale reading every
    # one of the nine tabulated BSV entries across both arms lands in a
    # coherent 11-75% CV band and each BSV tracks the magnitude of its own
    # typical value, which is what a natural-scale SD does and an omega has
    # no reason to do. See the vignette "Assumptions and deviations"
    # section.
    #   ka:   0.49 / 1.88  = 26.06% CV -> log(1 + 0.2606^2) = 0.06572437
    #   V:    5.60 / 22.20 = 25.23% CV -> log(1 + 0.2523^2) = 0.06168871
    #   k10:  0.09 / 0.43  = 20.93% CV -> log(1 + 0.2093^2) = 0.04287505
    #   k12:  0.45 / 4.09  = 11.00% CV -> log(1 + 0.1100^2) = 0.01203270
    #   k21:  0.33 / 0.44  = 75.00% CV -> log(1 + 0.7500^2) = 0.44628710
    # ------------------------------------------------------------------
    etalka  ~ 0.06572437  # Table 1, BSV Oral: ka  = 0.49 (RSE 33%)
    etalvc  ~ 0.06168871  # Table 1, BSV Oral: V   = 5.60 (RSE 13.6%)
    etalkel ~ 0.04287505  # Table 1, BSV Oral: k10 = 0.09 (RSE 26.5%)
    etalk12 ~ 0.01203270  # Table 1, BSV Oral: k12 = 0.45 (RSE 22%)
    etalk21 ~ 0.44628710  # Table 1, BSV Oral: k21 = 0.33 (RSE 29.3%)

    # ------------------------------------------------------------------
    # Residual error -- Methods 2.5: "Residual variability was described
    # by the proportional error model", Y_ij = F_ij * (1 + eps_ij). Table 1
    # reports the proportional coefficient as "b".
    # ------------------------------------------------------------------
    propSd <- 0.41 ; label("Proportional residual error (fraction)")  # Table 1, Mean Oral, Residual variability: b = 0.41 (RSE 12.3%)
  })

  model({
    # 1. Individual parameters (log-normal IIV, Methods 2.5).
    ka  <- exp(lka  + etalka)
    vc  <- exp(lvc  + etalvc)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 2. ODE system -- Cerqueira 2025 Results 3.2, oral equation triple.
    # The source prints the peripheral equation as
    #   dP/dt = +k21*P - k12*C
    # which is a sign transposition: as printed the peripheral compartment
    # is fed by itself and drains into the central compartment, so P grows
    # without bound and no drug ever distributes out of plasma. It also
    # contradicts the same paragraph's own definition of k12 and k21 as
    # "the distribution and redistribution rates" (printed there as a
    # second typo, "k12 and k12"). The standard form (k12*C in, k21*P out)
    # is used here. The depot and central equations are reproduced exactly
    # as printed. See vignette "Assumptions and deviations".
    #
    # The source writes C and P as concentrations and A as an amount while
    # reporting a single volume V; the states here all carry amount
    # (ug/kg), which is the same system scaled by V, is dimensionally
    # consistent, and is what Monolix's (ka, V, k, k12, k21)
    # parameterisation integrates.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-               k12 * central     - k21 * peripheral1

    # 3. Observation. Amount in ug/kg over apparent volume in L/kg gives
    # ug/L, which is ng/mL -- the assay units of Figures 1b and 2b.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
