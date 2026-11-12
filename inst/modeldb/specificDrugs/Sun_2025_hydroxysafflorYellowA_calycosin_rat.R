Sun_2025_hydroxysafflorYellowA_calycosin_rat <- function() {
  description <- paste(
    "Preclinical (rat). Coupled PK-PD model for intravenous co-administration",
    "of hydroxysafflor yellow A (HSYA, 2 mg/kg) and calycosin (CA, 1 mg/kg) in",
    "male Sprague-Dawley rats with middle cerebral artery occlusion",
    "(MCAO; Sun 2025). Each drug follows one-compartment first-order",
    "elimination whose rate is raised by the plasma concentration of the",
    "other drug (dC1/dt = -(k1 + k~1 * C2) * C1, dC2/dt = -(k2 + k~2 * C1) * C2;",
    "Eq. 1), started from fitted initial concentrations C10 and C20 rather",
    "than from a dose. Each drug drives its own effect compartment per",
    "biomarker, and the plasma Caspase-9 (U/L), IL-1beta (pg/mL) and SOD",
    "(U/mL) responses are the sum of two independent sigmoid terms",
    "alpha * E1 + beta * E2 (Eq. 7). Fitted to the mean of 6 rats; no",
    "between-subject variability or residual error was reported. ERRATA: the",
    "printed Table 2 PK parameters do not reproduce the paper's own Figure 7",
    "(predicted HSYA falls to 0.006 ug/mL at 2 h vs 0.65 observed), and the",
    "Table 4 PD 'mean' column equals one rat's estimates (rat 5 for Caspase-9,",
    "rat 6 for IL-1beta and SOD). Values are carried verbatim; see the",
    "vignette before using this model quantitatively."
  )
  reference <- paste(
    "Sun J, Chen Q, Zhuang C, Li X, Yu L, Jin W (2025). Mechanistic insights",
    "into synergistic effects using coupled PK-PD modeling. Sci Rep",
    "15:15631. doi:10.1038/s41598-025-00182-4"
  )
  vignette <- "Sun_2025_hydroxysafflorYellowA_calycosin_rat"
  units <- list(time = "h", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "hydroxysafflor yellow A", units = "ug/mL", specimen = "plasma", verified = TRUE),
    central_calycosin = list(analyte = "calycosin", units = "ug/mL", specimen = "plasma", verified = TRUE),
    effect_casp9 = list(
      analyte = "hydroxysafflor yellow A",
      units = "ug/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect_casp9_calycosin = list(analyte = "calycosin", units = "ug/mL", specimen = "not applicable", verified = TRUE),
    effect_il1b = list(
      analyte = "hydroxysafflor yellow A",
      units = "ug/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect_il1b_calycosin = list(analyte = "calycosin", units = "ug/mL", specimen = "not applicable", verified = TRUE),
    effect_sod = list(
      analyte = "hydroxysafflor yellow A",
      units = "ug/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    effect_sod_calycosin = list(analyte = "calycosin", units = "ug/mL", specimen = "not applicable", verified = TRUE)
  )

  paper_specific_compartments <- c(
    # Six parallel (not chained) effect compartments: each drug has its own
    # effect compartment per biomarker because Sun 2025 estimates a separate
    # ke1 and ke2 for every biomarker (Table 4). The numbered `effect<n>` form
    # would imply a chain, so each carries its biomarker suffix, following
    # Kiriyama_2024_nifedipine_rat.R (`effect_bp`, `effect_hr`).
    "effect_casp9",
    "effect_casp9_calycosin",
    "effect_il1b",
    "effect_il1b_calycosin",
    "effect_sod",
    "effect_sod_calycosin"
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 6,
    n_studies = 1,
    sex = "male",
    weight_range = "260-300 g",
    disease_state = "middle cerebral artery occlusion (MCAO) cerebral ischemia-reperfusion",
    dose_range = "single tail-vein co-injection of HSYA 2 mg/kg + calycosin 1 mg/kg at reperfusion",
    regions = "China",
    notes = paste(
      "Sun 2025 Experiments section: SPF male Sprague-Dawley rats (260-300 g)",
      "from Zhejiang Chinese Medical University, MCAO by the Longa method, n = 6.",
      "Plasma sampled at 0.083, 0.167, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 4, 6, 8,",
      "12 and 24 h; HSYA and calycosin by LC-MS/MS, Caspase-9, IL-1beta and SOD",
      "by ELISA on the same samples. Parameters were estimated in MATLAB with a",
      "genetic algorithm followed by alternating PK / PD optimisation against",
      "the per-time-point mean and SD of the six rats."
    )
  )

  ini({
    # Coupled PK (Sun 2025 Eq. 1; Table 2, mean of six rats; equal to the mean
    # of supplement Table 6). Concentrations in ug/mL, time in h. The model
    # starts from the fitted initial concentrations; it has no dose input.
    lkel <- log(1.570177); label("HSYA first-order elimination rate constant k1 (1/h)") # Table 2, k1, HSYA column
    lkel_calycosin <- log(1.413132); label("Calycosin first-order elimination rate constant k2 (1/h)") # Table 2, k2, CA column
    ks <- 47.406980; label("Enhancement of HSYA elimination per unit calycosin concentration k~1 (mL/ug/h)") # Table 2, k~1, HSYA column
    ks_calycosin <- 0.573247; label("Enhancement of calycosin elimination per unit HSYA concentration k~2 (mL/ug/h)") # Table 2, k~2, CA column
    lc0 <- log(4.441699); label("HSYA initial plasma concentration C10 (ug/mL)") # Table 2, C10, HSYA column
    lc0_calycosin <- log(0.150838); label("Calycosin initial plasma concentration C20 (ug/mL)") # Table 2, C20, CA column

    # Coupled PD (Sun 2025 Eq. 3-7; Table 4). Row order of Table 4 is ke1, ke2,
    # gamma1, EC50_1, gamma2, EC50_2, alpha, beta, eta1, eta2 (confirmed by the
    # labelled supplement Table 8). Suffix _calycosin marks the calycosin arm
    # (index 2 in the paper); the unsuffixed name is the HSYA arm (index 1).
    # The Table 4 'Mean' column equals the Caspase-9 estimates of rat 5 and the
    # IL-1beta and SOD estimates of rat 6 in supplement Table 8.
    # --- Caspase-9 (U/L) ---
    lke0_casp9 <- log(0.785024); label("HSYA effect-compartment rate constant ke1, Caspase-9 (1/h)") # Table 4, ke1, Caspase-9
    lke0_casp9_calycosin <- log(0.050304); label("Calycosin effect-compartment rate constant ke2, Caspase-9 (1/h)") # Table 4, ke2, Caspase-9
    lhill_casp9 <- log(3.108474); label("HSYA Hill coefficient gamma1, Caspase-9 (unitless)") # Table 4, gamma1, Caspase-9
    lec50_casp9 <- log(1.855617); label("HSYA effect-site EC50_1, Caspase-9 (ug/mL)") # Table 4, EC50_1, Caspase-9
    lhill_casp9_calycosin <- log(0.539493); label("Calycosin Hill coefficient gamma2, Caspase-9 (unitless)") # Table 4, gamma2, Caspase-9
    lec50_casp9_calycosin <- log(0.052921); label("Calycosin effect-site EC50_2, Caspase-9 (ug/mL)") # Table 4, EC50_2, Caspase-9
    lemax_casp9 <- log(2272.479702); label("HSYA effect weight alpha, Caspase-9 (U/L)") # Table 4, alpha, Caspase-9
    lemax_casp9_calycosin <- log(1706.856119); label("Calycosin effect weight beta, Caspase-9 (U/L)") # Table 4, beta, Caspase-9
    # --- IL-1beta (pg/mL) ---
    lke0_il1b <- log(0.481205); label("HSYA effect-compartment rate constant ke1, IL-1beta (1/h)") # Table 4, ke1, IL-1beta
    lke0_il1b_calycosin <- log(0.075605); label("Calycosin effect-compartment rate constant ke2, IL-1beta (1/h)") # Table 4, ke2, IL-1beta
    lhill_il1b <- log(5.309134); label("HSYA Hill coefficient gamma1, IL-1beta (unitless)") # Table 4, gamma1, IL-1beta
    lec50_il1b <- log(1.297095); label("HSYA effect-site EC50_1, IL-1beta (ug/mL)") # Table 4, EC50_1, IL-1beta
    lhill_il1b_calycosin <- log(0.512179); label("Calycosin Hill coefficient gamma2, IL-1beta (unitless)") # Table 4, gamma2, IL-1beta
    lec50_il1b_calycosin <- log(0.139198); label("Calycosin effect-site EC50_2, IL-1beta (ug/mL)") # Table 4, EC50_2, IL-1beta
    lemax_il1b <- log(2141.182784); label("HSYA effect weight alpha, IL-1beta (pg/mL)") # Table 4, alpha, IL-1beta
    lemax_il1b_calycosin <- log(1382.785005); label("Calycosin effect weight beta, IL-1beta (pg/mL)") # Table 4, beta, IL-1beta
    # --- SOD (U/mL) ---
    lke0_sod <- log(0.872465); label("HSYA effect-compartment rate constant ke1, SOD (1/h)") # Table 4, ke1, SOD
    lke0_sod_calycosin <- log(0.141701); label("Calycosin effect-compartment rate constant ke2, SOD (1/h)") # Table 4, ke2, SOD
    lhill_sod <- log(5.241084); label("HSYA Hill coefficient gamma1, SOD (unitless)") # Table 4, gamma1, SOD
    lec50_sod <- log(1.846627); label("HSYA effect-site EC50_1, SOD (ug/mL)") # Table 4, EC50_1, SOD
    lhill_sod_calycosin <- log(0.861183); label("Calycosin Hill coefficient gamma2, SOD (unitless)") # Table 4, gamma2, SOD
    lec50_sod_calycosin <- log(0.026937); label("Calycosin effect-site EC50_2, SOD (ug/mL)") # Table 4, EC50_2, SOD
    lemax_sod <- log(866.9796); label("HSYA effect weight alpha, SOD (U/mL)") # Table 4, alpha, SOD
    lemax_sod_calycosin <- log(600.197532); label("Calycosin effect weight beta, SOD (U/mL)") # Table 4, beta, SOD

    # No residual-error model was reported (per-time-point likelihood on the
    # six-rat mean and SD, Eq. 8-15); held at zero for simulation.
    addSd <- fixed(0); label("Additive residual SD, HSYA plasma concentration (ug/mL)") # not reported in Sun 2025
    addSd_calycosin <- fixed(0); label("Additive residual SD, calycosin plasma concentration (ug/mL)") # not reported in Sun 2025
    addSd_casp9 <- fixed(0); label("Additive residual SD, Caspase-9 (U/L)") # not reported in Sun 2025
    addSd_il1b <- fixed(0); label("Additive residual SD, IL-1beta (pg/mL)") # not reported in Sun 2025
    addSd_sod <- fixed(0); label("Additive residual SD, SOD (U/mL)") # not reported in Sun 2025
  })
  model({
    kel <- exp(lkel)
    kel_calycosin <- exp(lkel_calycosin)
    c0 <- exp(lc0)
    c0_calycosin <- exp(lc0_calycosin)

    ke0_casp9 <- exp(lke0_casp9)
    ke0_casp9_calycosin <- exp(lke0_casp9_calycosin)
    hill_casp9 <- exp(lhill_casp9)
    ec50_casp9 <- exp(lec50_casp9)
    hill_casp9_calycosin <- exp(lhill_casp9_calycosin)
    ec50_casp9_calycosin <- exp(lec50_casp9_calycosin)
    emax_casp9 <- exp(lemax_casp9)
    emax_casp9_calycosin <- exp(lemax_casp9_calycosin)

    ke0_il1b <- exp(lke0_il1b)
    ke0_il1b_calycosin <- exp(lke0_il1b_calycosin)
    hill_il1b <- exp(lhill_il1b)
    ec50_il1b <- exp(lec50_il1b)
    hill_il1b_calycosin <- exp(lhill_il1b_calycosin)
    ec50_il1b_calycosin <- exp(lec50_il1b_calycosin)
    emax_il1b <- exp(lemax_il1b)
    emax_il1b_calycosin <- exp(lemax_il1b_calycosin)

    ke0_sod <- exp(lke0_sod)
    ke0_sod_calycosin <- exp(lke0_sod_calycosin)
    hill_sod <- exp(lhill_sod)
    ec50_sod <- exp(lec50_sod)
    hill_sod_calycosin <- exp(lhill_sod_calycosin)
    ec50_sod_calycosin <- exp(lec50_sod_calycosin)
    emax_sod <- exp(lemax_sod)
    emax_sod_calycosin <- exp(lemax_sod_calycosin)

    # States are plasma concentrations (ug/mL); Eq. 1 with C1 = HSYA, C2 = CA.
    central(0) <- c0
    central_calycosin(0) <- c0_calycosin
    d/dt(central) <- -(kel + ks * central_calycosin) * central
    d/dt(central_calycosin) <- -(kel_calycosin + ks_calycosin * central) * central_calycosin

    # Effect compartments (Eq. 3), one pair per biomarker; initial value 0.
    d/dt(effect_casp9) <- ke0_casp9 * (central - effect_casp9)
    d/dt(effect_casp9_calycosin) <- ke0_casp9_calycosin * (central_calycosin - effect_casp9_calycosin)
    d/dt(effect_il1b) <- ke0_il1b * (central - effect_il1b)
    d/dt(effect_il1b_calycosin) <- ke0_il1b_calycosin * (central_calycosin - effect_il1b_calycosin)
    d/dt(effect_sod) <- ke0_sod * (central - effect_sod)
    d/dt(effect_sod_calycosin) <- ke0_sod_calycosin * (central_calycosin - effect_sod_calycosin)

    # Independent sigmoid effects combined linearly (Eq. 4-5): E = alpha * E1 + beta * E2.
    # The text describes fitting Caspase-9 and IL-1beta after E~ = max(E) - E;
    # Figure 8 follows the direct form, which is what is encoded here.
    casp9 <- emax_casp9 * effect_casp9^hill_casp9 / (ec50_casp9^hill_casp9 + effect_casp9^hill_casp9) +
      emax_casp9_calycosin * effect_casp9_calycosin^hill_casp9_calycosin /
        (ec50_casp9_calycosin^hill_casp9_calycosin + effect_casp9_calycosin^hill_casp9_calycosin)
    il1b <- emax_il1b * effect_il1b^hill_il1b / (ec50_il1b^hill_il1b + effect_il1b^hill_il1b) +
      emax_il1b_calycosin * effect_il1b_calycosin^hill_il1b_calycosin /
        (ec50_il1b_calycosin^hill_il1b_calycosin + effect_il1b_calycosin^hill_il1b_calycosin)
    sod <- emax_sod * effect_sod^hill_sod / (ec50_sod^hill_sod + effect_sod^hill_sod) +
      emax_sod_calycosin * effect_sod_calycosin^hill_sod_calycosin /
        (ec50_sod_calycosin^hill_sod_calycosin + effect_sod_calycosin^hill_sod_calycosin)

    Cc <- central
    Cc_calycosin <- central_calycosin

    Cc ~ add(addSd)
    Cc_calycosin ~ add(addSd_calycosin)
    casp9 ~ add(addSd_casp9)
    il1b ~ add(addSd_il1b)
    sod ~ add(addSd_sod)
  })
}
