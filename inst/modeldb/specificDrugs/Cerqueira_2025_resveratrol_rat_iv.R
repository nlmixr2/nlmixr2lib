Cerqueira_2025_resveratrol_rat_iv <- function() {
  description <- paste(
    "Preclinical (rat).",
    "Two-compartment population PK model for resveratrol",
    "(3,5,4'-trihydroxy-trans-stilbene) after a single 5 mg/kg",
    "intravenous bolus into the caudal lateral vein of male Wistar rats.",
    "Fitted in Monolix 2020R1 (SAEM) and parameterised in micro-constant",
    "form exactly as the authors report it: a single central volume V",
    "(L/kg), first-order elimination k10 from central, and",
    "inter-compartmental rate constants k12 and k21 between the central",
    "and peripheral compartments. Because every dose, volume and",
    "clearance term in the source is normalised to body weight, the model",
    "is coded on a per-kilogram basis: the dosed amount is ug/kg and the",
    "volume is L/kg, so central/vc lands directly in ng/mL, the units the",
    "assay and Figures 1a and 2a report. All disposition parameters are",
    "true (non-apparent) values because the dose was given intravenously.",
    "The intravenous arm of the paper; see",
    "Cerqueira_2025_resveratrol_rat_oral for the separately fitted",
    "100 mg/kg oral arm.",
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
      "Single 5 mg/kg intravenous bolus via the caudal lateral vein.",
      "Resveratrol (purity > 98%) formulated as a 3 mg/mL solution in",
      "10% DMSO made up to volume with saline."
    ),
    regions        = "Brazil (Federal University of Bahia, Salvador)",
    notes          = paste(
      "Serial caudal-vein sampling at 0.125, 0.25, 0.5, 1, 1.5, 2, 4, 8",
      "and 24 h post-dose; resveratrol quantified by HPLC-UV at 310 nm",
      "(LLOQ 62.5 ng/mL, calibration range 62.5-5000 ng/mL). Animals were",
      "kept on a 12 h light-dark cycle with ad libitum food and water;",
      "ethics protocol 43/2016, Federal University of Bahia. Plasma",
      "protein binding measured separately by ultrafiltration was",
      "79 +/- 5% and linear over 0.5-50 ug/mL (Results 3.4); it is not a",
      "model parameter. A separate group (n = 12) received 10 mg/kg iv",
      "for the descriptive tissue-distribution assay (Figure 4); no",
      "tissue-distribution model was fitted, so no tissue compartments",
      "appear here. See Cerqueira 2025 Methods 2.2, 2.5 and Table 1."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Cerqueira 2025 Table 1, "Mean i.v. (%RSE)"
    # column. Monolix micro-constant parameterisation (V, k10, k12, k21).
    # The paper's k10 is the canonical kel.
    # ------------------------------------------------------------------
    lvc  <- log(3.60) ; label("Central volume of distribution V (log L/kg)")                  # Table 1, Mean i.v.: V = 3.60 L/kg (RSE 19.9%)
    lkel <- log(0.17) ; label("Elimination rate constant from central k10 (log 1/h)")         # Table 1, Mean i.v.: k10 = 0.17 /h (RSE 23.1%)
    lk12 <- log(1.20) ; label("Rate constant central -> peripheral1 k12 (log 1/h)")           # Table 1, Mean i.v.: k12 = 1.20 /h (RSE 17.4%)
    lk21 <- log(0.26) ; label("Rate constant peripheral1 -> central k21 (log 1/h)")           # Table 1, Mean i.v.: k21 = 0.26 /h (RSE 29.7%)

    # ------------------------------------------------------------------
    # IIV -- Cerqueira 2025 Table 1, "BSV i.v. (%RSE)" column. Methods 2.5
    # declares log-normal IIV, theta_i = theta * exp(eta), var(eta) =
    # omega^2.
    #
    # SCALE OF THE BSV COLUMN. The paper does not state the scale of the
    # tabulated BSV values. They are read here as the SD of the individual
    # parameter on the NATURAL (untransformed) scale, i.e. an approximate
    # CV of BSV/typical, and converted with omega^2 = log(1 + CV^2).
    # The alternative reading -- BSV = omega, the SD of eta on the log
    # scale -- is untenable: it is falsified by the oral sibling arm,
    # where V has typical 22.20 and BSV 5.60, which as an omega is a CV of
    # roughly 5e8 % and cannot be simulated. Under the natural-scale
    # reading every one of the nine tabulated BSV entries across both arms
    # lands in a coherent 11-75% CV band and each BSV tracks the magnitude
    # of its own typical value, which is what a natural-scale SD does and
    # an omega has no reason to do. See the vignette "Assumptions and
    # deviations" section.
    #   V:    0.97 / 3.60 = 26.94% CV -> log(1 + 0.2694^2) = 0.07008590
    #   k10:  0.07 / 0.17 = 41.18% CV -> log(1 + 0.4118^2) = 0.15661921
    #   k12:  0.29 / 1.20 = 24.17% CV -> log(1 + 0.2417^2) = 0.05676096
    #   k21:  0.11 / 0.26 = 42.31% CV -> log(1 + 0.4231^2) = 0.16466160
    # ------------------------------------------------------------------
    etalvc  ~ 0.07008590  # Table 1, BSV i.v.: V   = 0.97 (RSE 29.4%)
    etalkel ~ 0.15661921  # Table 1, BSV i.v.: k10 = 0.07 (RSE 19.1%)
    etalk12 ~ 0.05676096  # Table 1, BSV i.v.: k12 = 0.29 (RSE 24.6%)
    etalk21 ~ 0.16466160  # Table 1, BSV i.v.: k21 = 0.11 (RSE 20.5%)

    # ------------------------------------------------------------------
    # Residual error -- Methods 2.5: "Residual variability was described
    # by the proportional error model", Y_ij = F_ij * (1 + eps_ij). Table 1
    # reports the proportional coefficient as "b".
    # ------------------------------------------------------------------
    propSd <- 0.22 ; label("Proportional residual error (fraction)")  # Table 1, Mean i.v., Residual variability: b = 0.22 (RSE 11.8%)
  })

  model({
    # 1. Individual parameters (log-normal IIV, Methods 2.5).
    vc  <- exp(lvc  + etalvc)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # 2. ODE system -- Cerqueira 2025 Results 3.2, iv equation pair.
    # The source prints the peripheral equation as
    #   dP/dt = +k21*P - k12*C
    # which is a sign transposition: as printed the peripheral compartment
    # is fed by itself and drains into the central compartment, so P grows
    # without bound and no drug ever distributes out of plasma. It also
    # contradicts the same paragraph's own definition of k12 and k21 as
    # "the distribution and redistribution rates, respectively", and it
    # would not be an "open two-compartment model". The standard form
    # (k12*C in, k21*P out) is used here. The central equation is
    # reproduced exactly as printed. See vignette "Assumptions and
    # deviations".
    #
    # The source writes C and P as concentrations while reporting a single
    # volume V; the states here carry amount (ug/kg), which is the same
    # system scaled by V and is what Monolix's (V, k, k12, k21)
    # parameterisation integrates.
    d/dt(central)     <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central     - k21 * peripheral1

    # 3. Observation. Amount in ug/kg over volume in L/kg gives ug/L,
    # which is ng/mL -- the assay units of Figures 1a and 2a.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
