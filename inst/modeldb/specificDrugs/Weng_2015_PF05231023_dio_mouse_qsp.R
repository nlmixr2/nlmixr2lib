Weng_2015_PF05231023_dio_mouse_qsp <- function() {
  description <- paste(
    "Preclinical (mouse, diet-induced obese). QSP. Adipose-target PK,",
    "receptor-occupancy and energy-balance body-weight model for",
    "PF-05231023, a long-acting FGF21 analogue (two FGF21 molecules",
    "conjugated to the CovX-2000 antibody scaffold). PK: SC depot",
    "(absorption half-life 48 h), central (2.4 mL) and adipose",
    "extracellular space (0.8 mL) exchanging by a symmetric flow of",
    "0.0167 mL/h; in every compartment the intact molecule loses its",
    "C-terminus (half-life 1.7 h) to give an N-only metabolite, and both",
    "lose activity by N-terminal clipping (half-life 26.2 h). In adipose,",
    "intact and N-only compete in rapid equilibrium for the KLB/FGFR1c",
    "receptor complex (4.15 nM; KD 0.5 and 94 nM), solved in closed",
    "form. Receptor occupancy drives a sigmoid effect P(RO) that",
    "suppresses the overfeeding response (beta) and adds a",
    "fat-proportional energy expenditure (epsilon); body weight follows",
    "the energy balance of food intake minus resting metabolic rate,",
    "with the vehicle-group body weight as an input curve.",
    "Deterministic: the paper fits mean responses by least squares and",
    "reports no between-animal or residual variability."
  )
  reference <- paste(
    "Weng Y, Chabot JR, Bernardo B, Yan Q, Zhu Y, Brenner MB, Vage C,",
    "Logan A, Calle R, Talukdar S. Pharmacokinetics (PK),",
    "pharmacodynamics (PD) and integrated PK/PD modeling of a novel long",
    "acting FGF21 clinical candidate PF-05231023 in diet-induced obese",
    "and leptin-deficient obese mice. PLoS One. 2015;10(3):e0119104.",
    "doi:10.1371/journal.pone.0119104."
  )
  vignette <- "Weng_2015_PF05231023_mouse_qsp"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "nM"
  )

  # The N-only metabolite (C-terminus clipped, N-terminus intact) is
  # carried alongside the intact molecule in each compartment.
  paper_specific_compartments <- c("depot_nonly", "central_nonly", "adipose_nonly")

  covariateData <- list(
    WT_BASE = list(
      description = "Baseline (day 0) body weight of the treated group",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets the initial condition of the body-weight state, bw(0) =",
        "1000 * WT_BASE g, and the baseline food intake FI(t = 0). The",
        "DIO mice of Fig 7 started at about 0.042-0.043 kg. PK volumes",
        "are fixed for a 45 g mouse and do not scale with WT_BASE."
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "PF-05231023 (intact)", units = "mg", specimen = "administration site", verified = TRUE),
    depot_nonly = list(
      analyte = "PF-05231023 N-only metabolite (C-terminus clipped)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "PF-05231023 (intact)", units = "mg", specimen = "plasma", verified = TRUE),
    central_nonly = list(
      analyte = "PF-05231023 N-only metabolite (C-terminus clipped)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    adipose = list(analyte = "PF-05231023 (intact)", units = "mg", specimen = "tissue", verified = TRUE),
    adipose_nonly = list(
      analyte = "PF-05231023 N-only metabolite (C-terminus clipped)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    bw = list(analyte = "body weight", units = "g", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "mouse (diet-induced obese, 60% high-fat diet)",
    n_subjects = NA_integer_,
    n_studies = NA_integer_,
    age_range = "PK: 9-10 weeks; body weight: about 14 weeks on 60% high-fat diet at the start (male)",
    weight_range = "about 42-43 g at day 0 (Fig 7b and 7c)",
    sex_female_pct = 0,
    disease_state = "Diet-induced obesity (60% high-fat diet)",
    dose_range = "1, 3 or 10 mg/kg IV, or 1 or 10 mg/kg SC, once weekly for three weeks (days 1, 8 and 15)",
    notes = paste(
      "PK in male DIO mice after a single 3 or 10 mg/kg IV or SC dose",
      "(Table 1 and Table 2; n = 3-4 per group). Body-weight study:",
      "male DIO mice, at least 9 per dose group, weighed twice weekly",
      "to day 57 (Fig 4); the model was fitted to all five treated",
      "groups at once (Fig 7). The strain background is not stated."
    )
  )

  ini({
    # PK (Methods, 'PK/PD modeling for BW loss in DIO mouse'; all fixed, none fitted)
    lka <- fixed(log(log(2) / 48)); label("First-order SC absorption rate constant (1/h)") # Methods: 'transfers to the central compartment with a characteristic time (i.e., half-life) of 48 hours'
    thalf_ct <- fixed(1.7); label("C-terminus clipping half-life, intact -> N-only (h)") # Methods: 'DIO mouse-measured half-lives of 1.7 hours and 26.2 hours were assumed for the intact CT and NT'
    thalf_nt <- fixed(26.2); label("N-terminus clipping half-life, loss of activity (h)") # Methods: same sentence
    lvc <- fixed(log(2.4)); label("Central volume (mL)") # Methods: 'central compartment with volume 2.4 ml (based on a 45 g body weight and 7% total Vss ...)'
    lvadipose <- fixed(log(0.8)); label("Adipose extracellular volume (mL)") # Methods: 'volume 0.8 ml (assuming 20 g of fat with 4% extracellular space)'
    lq <- fixed(log(0.0167)); label("Central-adipose exchange flow (mL/h)") # Methods: 'A first-order transport rate of 0.0167 ml/hr was symmetrically applied between these compartments' (= 0.8 mL / 48 h)
    lkp_adipose <- fixed(log(1)); label("Adipose:central equilibrium concentration ratio (unitless)") # Not printed; 'symmetrically applied' exchange implies 1 (Vss = Vc + Vadipose = 7% of 45 g). Not a paper parameter -- see the article for the reading that reproduces Fig 5b/7b/7c
    mw <- fixed(192000); label("Molecular weight of PF-05231023 (g/mol)") # Not printed; read from Fig 5a (10 mg/kg IV, ~1300 nM at t = 0 in 4% of body weight)

    # Receptor binding in adipose (Eq 4)
    rtot <- fixed(4.15); label("Total extracellular KLB/FGFR1c receptor complex in adipose (nM)") # Methods: 'set at 4.15 nM; this corresponds to 10,000 sites per adipocyte'
    kd_intact <- fixed(0.5); label("Dissociation constant of intact PF-05231023, KD,I (nM)") # Methods: 'set to 0.5 nM based on in vitro hASC measurements'
    kd_nonly <- 94; label("Dissociation constant of the N-only pool, KD,N (nM)") # Fig 7a: 'KD,N 94 nM'

    # Transfer function P(RO) (Eq 5) and energy balance (Eqs 6-13)
    delta <- 0.48; label("Receptor occupancy at half-maximal P, delta (unitless)") # Fig 7a: 'delta 0.48 (dimensionless)'
    rho <- 3; label("Sigmoidicity of P(RO), rho (unitless)") # Fig 7a: 'rho 3 (dimensionless)'
    beta <- 1.3; label("Overfeeding exponent, beta (unitless)") # Fig 7a: 'beta 1.3 (dimensionless)'
    eps <- 2.3; label("FGF21 energy expenditure per g body fat at P = 1, epsilon (kJ/(day*g))") # Fig 7a: 'epsilon 2.3 kJ/(day*g FFM)'; Eq 11 and the text apply it to total adipose tissue (TAT)
    rmr <- fixed(0.86); label("Resting metabolic rate (kJ/(day*g fat-free mass))") # Methods: 'measured as 0.86 kJ/(g fat free mass) [23]'
    tat_slope <- fixed(0.46); label("Slope of total adipose tissue on body weight (g/g)") # Eq 6: TAT = 0.46*BW - 9.9 g
    tat_int <- fixed(9.9); label("Intercept of total adipose tissue on body weight (g)") # Eq 6
    ed_fat <- fixed(39.5); label("Energy content of body fat (kJ/g)") # Methods: 'Body fat contains approximately 39.5 kJ/g' [22]
    ed_lean <- fixed(7.6); label("Energy content of lean tissue (kJ/g)") # Methods: 'lean tissue contains 7.6 kJ/g' [22]

    # Vehicle-group body weight BWcontrol(t) = bwc_inf - (bwc_inf - bwc_0)*exp(-kbwc*day).
    # Not printed: the paper smooths the vehicle data; these reproduce the
    # Fig 7b vehicle model curve (IV study), digitised and fitted by the maintainers.
    bwc_0 <- fixed(43.25); label("Vehicle-group body weight at day 0 (g)") # Fig 7b black curve, digitised
    bwc_inf <- fixed(56.83); label("Vehicle-group body-weight asymptote (g)") # Fig 7b black curve, digitised
    kbwc <- fixed(0.01226); label("Vehicle-group body-weight approach rate (1/day)") # Fig 7b black curve, digitised
  })

  model({
    ka <- exp(lka)
    kct <- log(2) / thalf_ct
    knt <- log(2) / thalf_nt
    vc <- exp(lvc)
    vadipose <- exp(lvadipose)
    q <- exp(lq)
    kp_adipose <- exp(lkp_adipose)

    # PK. Doses: IV into central, SC into depot (amounts in mg).
    flux_i <- q * (central / vc - adipose / (vadipose * kp_adipose))
    flux_n <- q * (central_nonly / vc - adipose_nonly / (vadipose * kp_adipose))
    d/dt(depot) <- -(ka + kct + knt) * depot
    d/dt(depot_nonly) <- kct * depot - (ka + knt) * depot_nonly
    d/dt(central) <- ka * depot - (kct + knt) * central - flux_i
    d/dt(central_nonly) <- ka * depot_nonly + kct * central - knt * central_nonly - flux_n
    d/dt(adipose) <- flux_i - (kct + knt) * adipose
    d/dt(adipose_nonly) <- flux_n + kct * adipose - knt * adipose_nonly

    # Concentrations in nM: mg/mL = g/L, / mw -> mol/L
    Cc <- central / vc * 1e9 / mw
    Cc_nonly <- central_nonly / vc * 1e9 / mw
    Cadipose <- adipose / vadipose * 1e9 / mw
    Cadipose_nonly <- adipose_nonly / vadipose * 1e9 / mw

    # Eq 4: rapid-equilibrium competitive binding in adipose. Free receptor
    # is the root in (0, rtot] of R^3 + a R^2 + b R + c = 0 (trigonometric
    # solution of the cubic for two ligands competing for one site).
    cub_a <- kd_intact + kd_nonly + Cadipose + Cadipose_nonly - rtot
    cub_b <- kd_nonly * (Cadipose - rtot) + kd_intact * (Cadipose_nonly - rtot) + kd_intact * kd_nonly
    cub_c <- -kd_intact * kd_nonly * rtot
    cub_arg <- (-2 * cub_a^3 + 9 * cub_a * cub_b - 27 * cub_c) / (2 * sqrt((cub_a^2 - 3 * cub_b)^3))
    cub_arg <- min(max(cub_arg, -1), 1)
    rfree <- -cub_a / 3 + 2 / 3 * sqrt(cub_a^2 - 3 * cub_b) * cos(acos(cub_arg) / 3)
    ro <- (Cadipose / (kd_intact + rfree) + Cadipose_nonly / (kd_nonly + rfree)) * rfree / rtot

    # Eq 5: sigmoid transfer function, normalised so that P(RO = 1) = 1
    peff <- (delta^rho + 1) * ro^rho / (delta^rho + ro^rho)

    # Vehicle-group body weight and food intake (Eq 8; g, g/day, kJ/day)
    day <- t / 24
    ed <- tat_slope * ed_fat + (1 - tat_slope) * ed_lean
    bwc <- bwc_inf - (bwc_inf - bwc_0) * exp(-kbwc * day)
    dbwc <- kbwc * (bwc_inf - bwc_0) * exp(-kbwc * day)
    fic <- rmr * (bwc - (tat_slope * bwc - tat_int)) + dbwc * ed

    # Treated group (Eqs 6, 9-13)
    bw(0) <- 1000 * WT_BASE
    tat <- tat_slope * bw - tat_int
    ffm <- bw - tat
    fic_0 <- rmr * (bwc_0 - (tat_slope * bwc_0 - tat_int)) + kbwc * (bwc_inf - bwc_0) * ed
    fi_0 <- fic_0 * (bwc_0 / (1000 * WT_BASE))^beta
    fi <- fi_0 * peff + fic * (bwc / bw)^beta * (1 - peff)
    de_fgf21 <- eps * peff * tat
    de <- fi - rmr * ffm - de_fgf21
    d/dt(bw) <- de / ed / 24
  })
}
