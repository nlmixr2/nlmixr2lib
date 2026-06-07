Chatelut_1999_interferon_alfa_2b <- function() {
  description <- "One-compartment population PK model for subcutaneous alpha-2b interferon (Intron A) in adults with chronic hepatitis C virus infection (Chatelut 1999), with sequential zero-order then first-order absorption (a fraction Fz of the bioavailable dose is absorbed at zero-order over duration tk0, the remaining (1 - Fz) is absorbed at first-order rate ka after tk0) and first-order elimination. Apparent oral clearance CL/F is reduced by 63.8% in chronic-haemodialysis patients relative to patients with normal renal function (HEMODIAL = 1 vs 0); apparent central volume of distribution V/F scales linearly with body surface area (BSA). Proportional residual error."
  reference <- "Chatelut E, Rostaing L, Gregoire N, Payen JL, Pujol A, Izopet J, Houin G, Canal P. A pharmacokinetic model for alpha interferon administered subcutaneously. Br J Clin Pharmacol. 1999 Apr;47(4):365-71. doi:10.1046/j.1365-2125.1999.00912.x"
  vignette <- "Chatelut_1999_interferon_alfa_2b"
  units <- list(time = "hour", dosing = "ng", concentration = "pg/mL")

  covariateData <- list(
    HEMODIAL = list(
      description        = "Chronic intermittent haemodialysis treatment-status indicator (1 = chronic haemodialysis patient; 0 = patient with normal renal function).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal renal function)",
      notes              = "Chatelut 1999 Methods/Results uses the symbol DIA (1 if dialysis, 0 otherwise) as a multiplicative covariate on apparent oral clearance: CL/F = theta1 * (1 - theta2 * DIA). The 17 chronic-haemodialysis patients had been on dialysis for more than 4 years; the PK study at the first interferon-alfa-2b dose was performed 8 h after the last dialysis session, and the subsequent dialysis session never occurred before the last blood sample was taken, so for the modelled period the indicator is time-fixed at the subject level. Maps directly to the canonical HEMODIAL covariate (paper symbol DIA).",
      source_name        = "DIA"
    ),
    BSA = list(
      description        = "Body surface area, computed by the DuBois formula from height and weight.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chatelut 1999 Results 'Structural and covariate models' identified BSA as the covariate explaining V/F variability ('V/F was significantly correlated with the body weight, but the correlation was better with the body surface area'). The published structural form V/F = theta * BSA is a linear-proportional scaling (exponent fixed at 1, equivalent to expressing V/F in units of L/m^2). BSA was computed using the DuBois formula (paper Methods, Table 1 footnote). Baseline / time-fixed in this single-dose PK study.",
      source_name        = "BSA"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 27L,
    n_studies       = 1L,
    age_range       = "25-68 years",
    age_median      = "43 years (dialysis), 57 years (normal renal)",
    weight_range    = "41-100 kg",
    weight_median   = "60 kg (dialysis), 74 kg (normal renal)",
    bsa_range       = "1.35-2.26 m^2",
    bsa_median      = "1.69 m^2 (dialysis), 1.75 m^2 (normal renal)",
    sex_female_pct  = 44.4,
    disease_state   = "Chronic hepatitis C virus infection (positive HCV serology, positive HCV viraemia, chronic active hepatitis on liver biopsy). 17 patients on chronic intermittent haemodialysis (more than 4 years on dialysis); 10 patients with normal renal function (median CRCL 79 mL/min, range 48-136).",
    dose_range      = "Single subcutaneous injection of 3,000,000 units (15,000 ng) of alpha-2b interferon (Intron A, Schering Plough) at the start of a 3-times-weekly therapy course scheduled for 1 year. Pharmacokinetic study performed at the time of the first injection. In dialysis patients, the first injection was given 8 h after the last dialysis session.",
    regions         = "France (Toulouse).",
    notes           = "Demographics from Chatelut 1999 Table 1. Plasma concentrations measured by human alpha-interferon ELISA (ENDOGEN; LOQ 4.1 pg/mL, LOD <3 pg/mL); samples drawn pre-dose and at 1, 2, 3, 4, 6, 8, 12, 16, 20, 24, 28, and 32 h after the SC injection. Pre-dose samples were below LOQ in all patients, confirming negligible endogenous interferon interference."
  )

  ini({
    # Structural parameters -- Chatelut 1999 Table 3 ("Pharmacokinetic parameters
    # of alpha interferon"). Typical values are reported as the mean estimate
    # with 95% confidence interval; values below are the mean point estimates.

    lfr_zo <- log(0.24);   label("Fraction of bioavailable dose absorbed via the zero-order arm, Fz (unitless)")  # Table 3: Fz = 0.24 (95% CI 0.16-0.32)
    ltk0   <- log(2.5);    label("Duration of zero-order absorption phase, tk0 (h)")                              # Table 3: tk0 = 2.5 h (95% CI 1.9-3.0)
    lka    <- log(0.18);   label("First-order absorption rate constant for the post-zero-order phase, ka (1/h)")  # Table 3: ka = 0.18 1/h (95% CI 0.14-0.22)
    lvc    <- log(91);     label("Apparent central volume of distribution per m^2 of BSA, V/F (L/m^2)")           # Table 3: V/F = 91 L/m^2 (95% CI 72-110)
    lcl    <- log(36.5);   label("Apparent oral clearance CL/F in patients with normal renal function (L/h)")     # Table 3: CL/F (normal renal) = 36.5 L/h (95% CI 30.3-42.7)

    # Covariate effects.
    # Chatelut 1999 Results 'Structural and covariate models' (paragraph 2):
    # CL/F = theta1 * (1 - theta2 * DIA) with theta2 = 0.638 (95% CI 0.568-0.708).
    # The estimate is reported with a CI -> the coefficient was estimated (not
    # fixed); CL/F in dialysis = 36.5 * (1 - 0.638) = 13.2 L/h matches Table 3.
    e_hemodial_cl <- 0.638; label("Proportional reduction in CL/F for chronic-haemodialysis patients (unitless)")  # Results, paragraph 2: theta2 = 0.638 (95% CI 0.568-0.708)

    # Inter-individual variability. Table 3 reports interindividual variability
    # as %CV with 95% CI. The paper Methods state 'A proportional error model
    # was used for the interpatient variabilities', i.e. theta_i = theta * exp(eta).
    # Convert CV% to log-scale variance: omega^2 = log(1 + CV^2).
    etalfr_zo ~ 0.1034   # Table 3: Fz IIV 33% CV -> log(1 + 0.33^2) = log(1.1089) = 0.1034
    etaltk0   ~ 0.1034   # Table 3: tk0 IIV 33% CV -> same calculation
    etalka    ~ 0.1484   # Table 3: ka IIV 40% CV -> log(1 + 0.40^2) = log(1.1600) = 0.1484
    etalvc    ~ 0.03922  # Table 3: V/F IIV 20% CV -> log(1 + 0.20^2) = log(1.04)   = 0.03922
    etalcl    ~ 0.1349   # Table 3: CL/F IIV 38% CV -> log(1 + 0.38^2) = log(1.1444) = 0.1349

    # Residual error. Methods state 'A proportional error model was used for the
    # interpatient variabilities. Preliminary testing of combination model (i.e.
    # additive plus proportional) for residual variability led to negligible
    # value for the additional term. Then, a proportional model was used for
    # residual variability.' Table 2 reports the final-model residual variability
    # as sigma = 22 %.
    propSd <- 0.22; label("Proportional residual error (fraction)")  # Table 2: sigma = 22% for the final model
  })

  model({
    # Individual PK parameters.
    fr_zo <- exp(lfr_zo + etalfr_zo)
    tk0   <- exp(ltk0 + etaltk0)
    ka    <- exp(lka  + etalka)

    # V/F is linear-proportional to BSA (Chatelut 1999 reports V/F in L/m^2;
    # the structural form V/F = theta * BSA has an implicit exponent of 1).
    vc    <- exp(lvc + etalvc) * BSA

    # CL/F = theta1 * (1 - theta2 * HEMODIAL) (Chatelut 1999 Results, paragraph 2).
    cl    <- exp(lcl + etalcl) * (1 - e_hemodial_cl * HEMODIAL)

    kel <- cl / vc

    # Sequential zero-order then first-order absorption from the SC depot
    # (Chatelut 1999 Figure 1 + Results 'Structural and covariate models',
    # paragraph 1: 'one-compartment model with a zero-order absorption
    # immediately followed by a first-order absorption'). During the
    # interval [0, tk0] the depot empties at the constant zero-order rate
    # `fr_zo * Dose / tk0` (delivering Fz * Dose to central); after tk0 the
    # remaining (1 - Fz) * Dose absorbs at first-order rate ka. This is the
    # rxode2 idiom used in Horita_2018_rifampicin.R for the same structural
    # model.
    mtime(tk_switch) <- tk0

    kzero <- fr_zo * podo(depot) / tk0
    if (tad(depot) > tk0) kzero <- 0.0

    ka_eff <- ka
    if (tad(depot) <= tk0) ka_eff <- 0.0

    d/dt(depot)   <- -ka_eff * depot - kzero
    d/dt(central) <-  ka_eff * depot + kzero - kel * central

    # Concentration in pg/mL: dose in ng / V in L = ng/L = pg/mL (matches the
    # source paper, which reports plasma alpha-interferon concentrations and
    # AUC in pg/mL and pg/mL.h).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
