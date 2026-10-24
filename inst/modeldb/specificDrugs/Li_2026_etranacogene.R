Li_2026_etranacogene <- function() {
  description <- paste(
    "QSP (agent-based gene therapy, deterministic generation-ladder",
    "encoding). Human liver model for the durability of factor IX (FIX)",
    "activity after AAV5 gene therapy with etranacogene dezaparvovec",
    "(AMT-061 / CSL222; Padua R338L FIX transgene) in hemophilia B",
    "(Li 2026). The published model is an agent-based model (ABM) in which",
    "each hepatocyte is an agent that divides, dies and switches ploidy",
    "(diploid vs polyploid); transduced cells carry episomes that are",
    "split between daughters at division and lost at death, and express",
    "FIX from at most EPISOME_CAP episomes per cell. This file encodes the",
    "ABM's expected behaviour exactly (under the paper's even-episome",
    "split) as a deterministic division-count ladder: ND0..ND15 and",
    "NP0..NP15 are the transduced diploid / polyploid hepatocytes grouped",
    "by number of divisions since treatment, WD/WP carry the diluted-",
    "episome tail beyond 15 divisions, and TD/TP are total diploid /",
    "polyploid hepatocytes for normalisation. FIX is produced at quasi-",
    "steady state (half-life ~7 h << dynamics) and scaled to percent of",
    "normal. Deterministic: no IIV or residual error; population",
    "variability comes from sampling the covariates over the Table S2",
    "ranges. See the vignette for the approximations vs the authors' code.",
    sep = " "
  )
  reference <- paste(
    "Li Y, Nandy P, Jordie E, et al.",
    "Development of an Agent-Based Model to Investigate Durability of",
    "Factor IX Activity in Hemophilia B Patients Treated With",
    "Etranacogene Dezaparvovec.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(7):e70286.",
    "doi:10.1002/psp4.70286.",
    "Rate parameters from Supporting Information Table S1; virtual-liver",
    "covariate ranges from Table S2; FIX-activity scaling and the",
    "diploid-fraction-vs-age relation from Materials and Methods.",
    sep = " "
  )
  vignette <- "Li_2026_etranacogene"

  paper_specific_compartments <- c(
    "ND0",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND5",
    "ND6",
    "ND7",
    "ND8",
    "ND9",
    "ND10",
    "ND11",
    "ND12",
    "ND13",
    "ND14",
    "ND15",
    "NP0",
    "NP1",
    "NP2",
    "NP3",
    "NP4",
    "NP5",
    "NP6",
    "NP7",
    "NP8",
    "NP9",
    "NP10",
    "NP11",
    "NP12",
    "NP13",
    "NP14",
    "NP15",
    "WD",
    "WP",
    "TD",
    "TP"
  )

  units <- list(
    time = "year",
    dosing = NA_character_,
    concentration = "FIX activity (% of normal)"
  )

  compartmentData <- list(
    ND0 = list(
      analyte = "transduced diploid hepatocytes after 0 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND1 = list(
      analyte = "transduced diploid hepatocytes after 1 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND2 = list(
      analyte = "transduced diploid hepatocytes after 2 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND3 = list(
      analyte = "transduced diploid hepatocytes after 3 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND4 = list(
      analyte = "transduced diploid hepatocytes after 4 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND5 = list(
      analyte = "transduced diploid hepatocytes after 5 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND6 = list(
      analyte = "transduced diploid hepatocytes after 6 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND7 = list(
      analyte = "transduced diploid hepatocytes after 7 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND8 = list(
      analyte = "transduced diploid hepatocytes after 8 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND9 = list(
      analyte = "transduced diploid hepatocytes after 9 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND10 = list(
      analyte = "transduced diploid hepatocytes after 10 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND11 = list(
      analyte = "transduced diploid hepatocytes after 11 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND12 = list(
      analyte = "transduced diploid hepatocytes after 12 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND13 = list(
      analyte = "transduced diploid hepatocytes after 13 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND14 = list(
      analyte = "transduced diploid hepatocytes after 14 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    ND15 = list(
      analyte = "transduced diploid hepatocytes after 15 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP0 = list(
      analyte = "transduced polyploid hepatocytes after 0 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP1 = list(
      analyte = "transduced polyploid hepatocytes after 1 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP2 = list(
      analyte = "transduced polyploid hepatocytes after 2 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP3 = list(
      analyte = "transduced polyploid hepatocytes after 3 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP4 = list(
      analyte = "transduced polyploid hepatocytes after 4 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP5 = list(
      analyte = "transduced polyploid hepatocytes after 5 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP6 = list(
      analyte = "transduced polyploid hepatocytes after 6 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP7 = list(
      analyte = "transduced polyploid hepatocytes after 7 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP8 = list(
      analyte = "transduced polyploid hepatocytes after 8 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP9 = list(
      analyte = "transduced polyploid hepatocytes after 9 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP10 = list(
      analyte = "transduced polyploid hepatocytes after 10 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP11 = list(
      analyte = "transduced polyploid hepatocytes after 11 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP12 = list(
      analyte = "transduced polyploid hepatocytes after 12 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP13 = list(
      analyte = "transduced polyploid hepatocytes after 13 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP14 = list(
      analyte = "transduced polyploid hepatocytes after 14 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    NP15 = list(
      analyte = "transduced polyploid hepatocytes after 15 division(s) since treatment",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    WD = list(
      analyte = "diluted-episome weight of transduced diploid hepatocytes beyond 15 divisions",
      units = "episome-equivalents per initial hepatocyte",
      specimen = "tissue",
      verified = TRUE
    ),
    WP = list(
      analyte = "diluted-episome weight of transduced polyploid hepatocytes beyond 15 divisions",
      units = "episome-equivalents per initial hepatocyte",
      specimen = "tissue",
      verified = TRUE
    ),
    TD = list(
      analyte = "total diploid hepatocytes (transduced and untransduced)",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    ),
    TP = list(
      analyte = "total polyploid hepatocytes (transduced and untransduced)",
      units = "fraction of initial hepatocyte count",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Patient age at treatment; sets the initial diploid hepatocyte fraction (92 - 0.5125 * AGE) / 100.",
      units = "years",
      type = "continuous",
      source_name = "patient_age"
    ),
    TRANSDUCED_FRAC = list(
      description = "Fraction of hepatocytes transduced by AAV at treatment (Table S2 initially_transduced_cells, 10-70%, as a fraction).",
      units = "fraction",
      type = "continuous",
      source_name = "initially_transduced_cells"
    ),
    EPISOME_MEAN = list(
      description = "Mean episome count per transduced hepatocyte (exponential-distribution mean; Table S2 Init_Episome, log-uniform 10-1000).",
      units = "episomes/hepatocyte",
      type = "continuous",
      source_name = "Init_Episome"
    ),
    EPISOME_CAP = list(
      description = "Maximum number of episomes expressed per hepatocyte (Table S2 max_episome_express, integer 1-5).",
      units = "episomes",
      type = "continuous",
      source_name = "max_episome_express"
    )
  )

  population <- list(
    species = "human",
    disease_state = "Adults with hemophilia B treated with etranacogene dezaparvovec (HOPE-B, NCT03569891; Phase 2b NCT03489291).",
    age_range = "20-80 years (virtual-liver range, Table S2)",
    sex = "male adult livers (Li 2026 Methods)",
    dose_range = "Single IV 2e13 genome copies/kg (one dose level; dose is not an explicit model input).",
    n_subjects = 160000,
    notes = paste(
      "The virtual population of 160,000 livers is built by sampling AGE,",
      "TRANSDUCED_FRAC, EPISOME_MEAN (log-uniform) and EPISOME_CAP over the",
      "Table S2 ranges by Latin hypercube. Rate parameters are fixed from",
      "literature (Heinke 2022 hepatocyte turnover; FIX half-life ~7 h).",
      sep = " "
    )
  )

  ini({
    lprob_div_di <- fixed(log(0.709)); label("Diploid hepatocyte proliferation rate (1/year)") # Table S1 beta2 = 0.709
    lprob_die_di <- fixed(log(0.726)); label("Diploid hepatocyte death rate (1/year)") # Table S1 delta2 = 0.726
    lprob_div_poly <- fixed(log(0.027)); label("Polyploid hepatocyte proliferation rate (1/year)") # Table S1 betap = 0.027 ('[1], tuned')
    lprob_die_poly <- fixed(log(0.019)); label("Polyploid hepatocyte death rate (1/year)") # Table S1 deltap = 0.019 ('[1], tuned')
    lk_di_to_poly <- fixed(log(0.0054)); label("Diploid-to-polyploid transition rate (1/year)") # Table S1 K2p = 0.0054
    lk_poly_to_di <- fixed(log(0.0167)); label("Polyploid-to-diploid transition rate (1/year)") # Table S1 Kp2 = 0.0167
    lkprod_di <- fixed(log(3.15e7)); label("FIX production per expressing episome, diploid (molecules/year)") # Table S1 kprodp diploid = 3.15e7
    lkprod_poly <- fixed(log(6.3e7)); label("FIX production per expressing episome, polyploid (molecules/year)") # Table S1 kprodp polyploid = 6.3e7 (2x diploid)
    lkdeg_fix <- fixed(log(886.16)); label("FIX degradation rate (1/year)") # Table S1 kdeg = 886.16 (half-life ~7 h)
  })

  model({
    prob_div_di <- exp(lprob_div_di)
    prob_die_di <- exp(lprob_die_di)
    prob_div_poly <- exp(lprob_div_poly)
    prob_die_poly <- exp(lprob_die_poly)
    k_di_to_poly <- exp(lk_di_to_poly)
    k_poly_to_di <- exp(lk_poly_to_di)
    kprod_di <- exp(lkprod_di)
    kprod_poly <- exp(lkprod_poly)
    kdeg_fix <- exp(lkdeg_fix)

    # Per-subject constants from the Table S2 covariates.
    diploid_frac <- (92 - 0.5125 * AGE) / 100 # Methods: diploid fraction vs age
    q <- exp(-1 / EPISOME_MEAN) # ceil-exponential episome count: P(n >= k) = q^(k-1)
    Einv <- 1 / (1 - q) # E[n] = mean episomes per transduced cell
    lossD <- prob_div_di + prob_die_di + k_di_to_poly
    lossP <- prob_div_poly + prob_die_poly + k_poly_to_di

    # Expected expressing episomes per living generation-g cell:
    # expr(g) = E[min(EPISOME_CAP, n / 2^g)] = 2^-g (1 - q^(cap 2^g)) / (1 - q).
    expr0 <- (2^(-0)) * (1 - q^(EPISOME_CAP * 1)) / (1 - q)
    expr1 <- (2^(-1)) * (1 - q^(EPISOME_CAP * 2)) / (1 - q)
    expr2 <- (2^(-2)) * (1 - q^(EPISOME_CAP * 4)) / (1 - q)
    expr3 <- (2^(-3)) * (1 - q^(EPISOME_CAP * 8)) / (1 - q)
    expr4 <- (2^(-4)) * (1 - q^(EPISOME_CAP * 16)) / (1 - q)
    expr5 <- (2^(-5)) * (1 - q^(EPISOME_CAP * 32)) / (1 - q)
    expr6 <- (2^(-6)) * (1 - q^(EPISOME_CAP * 64)) / (1 - q)
    expr7 <- (2^(-7)) * (1 - q^(EPISOME_CAP * 128)) / (1 - q)
    expr8 <- (2^(-8)) * (1 - q^(EPISOME_CAP * 256)) / (1 - q)
    expr9 <- (2^(-9)) * (1 - q^(EPISOME_CAP * 512)) / (1 - q)
    expr10 <- (2^(-10)) * (1 - q^(EPISOME_CAP * 1024)) / (1 - q)
    expr11 <- (2^(-11)) * (1 - q^(EPISOME_CAP * 2048)) / (1 - q)
    expr12 <- (2^(-12)) * (1 - q^(EPISOME_CAP * 4096)) / (1 - q)
    expr13 <- (2^(-13)) * (1 - q^(EPISOME_CAP * 8192)) / (1 - q)
    expr14 <- (2^(-14)) * (1 - q^(EPISOME_CAP * 16384)) / (1 - q)
    expr15 <- (2^(-15)) * (1 - q^(EPISOME_CAP * 32768)) / (1 - q)

    # Transduced diploid hepatocytes by division count (generation ladder).
    d/dt(ND0) <- -lossD * ND0 + k_poly_to_di * NP0
    d/dt(ND1) <- 2 * prob_div_di * ND0 - lossD * ND1 + k_poly_to_di * NP1
    d/dt(ND2) <- 2 * prob_div_di * ND1 - lossD * ND2 + k_poly_to_di * NP2
    d/dt(ND3) <- 2 * prob_div_di * ND2 - lossD * ND3 + k_poly_to_di * NP3
    d/dt(ND4) <- 2 * prob_div_di * ND3 - lossD * ND4 + k_poly_to_di * NP4
    d/dt(ND5) <- 2 * prob_div_di * ND4 - lossD * ND5 + k_poly_to_di * NP5
    d/dt(ND6) <- 2 * prob_div_di * ND5 - lossD * ND6 + k_poly_to_di * NP6
    d/dt(ND7) <- 2 * prob_div_di * ND6 - lossD * ND7 + k_poly_to_di * NP7
    d/dt(ND8) <- 2 * prob_div_di * ND7 - lossD * ND8 + k_poly_to_di * NP8
    d/dt(ND9) <- 2 * prob_div_di * ND8 - lossD * ND9 + k_poly_to_di * NP9
    d/dt(ND10) <- 2 * prob_div_di * ND9 - lossD * ND10 + k_poly_to_di * NP10
    d/dt(ND11) <- 2 * prob_div_di * ND10 - lossD * ND11 + k_poly_to_di * NP11
    d/dt(ND12) <- 2 * prob_div_di * ND11 - lossD * ND12 + k_poly_to_di * NP12
    d/dt(ND13) <- 2 * prob_div_di * ND12 - lossD * ND13 + k_poly_to_di * NP13
    d/dt(ND14) <- 2 * prob_div_di * ND13 - lossD * ND14 + k_poly_to_di * NP14
    d/dt(ND15) <- 2 * prob_div_di * ND14 - lossD * ND15 + k_poly_to_di * NP15
    # Transduced polyploid hepatocytes by division count.
    d/dt(NP0) <- -lossP * NP0 + k_di_to_poly * ND0
    d/dt(NP1) <- 2 * prob_div_poly * NP0 - lossP * NP1 + k_di_to_poly * ND1
    d/dt(NP2) <- 2 * prob_div_poly * NP1 - lossP * NP2 + k_di_to_poly * ND2
    d/dt(NP3) <- 2 * prob_div_poly * NP2 - lossP * NP3 + k_di_to_poly * ND3
    d/dt(NP4) <- 2 * prob_div_poly * NP3 - lossP * NP4 + k_di_to_poly * ND4
    d/dt(NP5) <- 2 * prob_div_poly * NP4 - lossP * NP5 + k_di_to_poly * ND5
    d/dt(NP6) <- 2 * prob_div_poly * NP5 - lossP * NP6 + k_di_to_poly * ND6
    d/dt(NP7) <- 2 * prob_div_poly * NP6 - lossP * NP7 + k_di_to_poly * ND7
    d/dt(NP8) <- 2 * prob_div_poly * NP7 - lossP * NP8 + k_di_to_poly * ND8
    d/dt(NP9) <- 2 * prob_div_poly * NP8 - lossP * NP9 + k_di_to_poly * ND9
    d/dt(NP10) <- 2 * prob_div_poly * NP9 - lossP * NP10 + k_di_to_poly * ND10
    d/dt(NP11) <- 2 * prob_div_poly * NP10 - lossP * NP11 + k_di_to_poly * ND11
    d/dt(NP12) <- 2 * prob_div_poly * NP11 - lossP * NP12 + k_di_to_poly * ND12
    d/dt(NP13) <- 2 * prob_div_poly * NP12 - lossP * NP13 + k_di_to_poly * ND13
    d/dt(NP14) <- 2 * prob_div_poly * NP13 - lossP * NP14 + k_di_to_poly * ND14
    d/dt(NP15) <- 2 * prob_div_poly * NP14 - lossP * NP15 + k_di_to_poly * ND15
    # Diluted-episome tail beyond 15 divisions (division conserves the weight).
    d/dt(WD) <- prob_div_di * ND15 * (2^(-15)) - (prob_die_di + k_di_to_poly) * WD + k_poly_to_di * WP
    d/dt(WP) <- prob_div_poly * NP15 * (2^(-15)) - (prob_die_poly + k_poly_to_di) * WP + k_di_to_poly * WD
    # Total hepatocytes (transduced + untransduced) for the normalisation denominator.
    d/dt(TD) <- (prob_div_di - prob_die_di - k_di_to_poly) * TD + k_poly_to_di * TP
    d/dt(TP) <- (prob_div_poly - prob_die_poly - k_poly_to_di) * TP + k_di_to_poly * TD

    ND0(0) <- diploid_frac * TRANSDUCED_FRAC
    NP0(0) <- (1 - diploid_frac) * TRANSDUCED_FRAC
    TD(0) <- diploid_frac
    TP(0) <- 1 - diploid_frac

    expr_epi_di <- ND0 * expr0 + ND1 * expr1 + ND2 * expr2 + ND3 * expr3 + ND4 * expr4 + ND5 * expr5 + ND6 * expr6 + ND7 * expr7 + ND8 * expr8 + ND9 * expr9 + ND10 * expr10 + ND11 * expr11 + ND12 * expr12 + ND13 * expr13 + ND14 * expr14 + ND15 * expr15 + WD * Einv
    expr_epi_poly <- NP0 * expr0 + NP1 * expr1 + NP2 * expr2 + NP3 * expr3 + NP4 * expr4 + NP5 * expr5 + NP6 * expr6 + NP7 * expr7 + NP8 * expr8 + NP9 * expr9 + NP10 * expr10 + NP11 * expr11 + NP12 * expr12 + NP13 * expr13 + NP14 * expr14 + NP15 * expr15 + WP * Einv
    fix_molec <- (kprod_di * expr_epi_di + kprod_poly * expr_epi_poly) / kdeg_fix
    hep_total <- TD + TP
    # Normal FIX molecules per hepatocyte: 87.8 nmol/L * 3.1 L plasma * Avogadro
    # divided by 1561 g liver * 139e6 hepatocytes/g (Methods).
    normal_fix_per_hep <- (87.8e-9 * 3.1 * 6.02214076e23) / (1561 * 139e6)
    # Padua R338L FIX has 6x wild-type specific activity (Methods).
    FIXactivity <- 100 * 6 * (fix_molec / hep_total) / normal_fix_per_hep
    # FIX activity attributable to each ploidy (Figure S1): the diploid share
    # decays over the first ~6 years, the long-lived polyploid share sustains it.
    FIXactivity_di <- 100 * 6 * ((kprod_di * expr_epi_di / kdeg_fix) / hep_total) / normal_fix_per_hep
    FIXactivity_poly <- 100 * 6 * ((kprod_poly * expr_epi_poly / kdeg_fix) / hep_total) / normal_fix_per_hep
  })
}
