Yamaguchi_2026_s531011_human <- function() {
  description <- "QSP. Semi-mechanistic PK / receptor-occupancy (RO) model of S-531011, a humanized anti-CCR8 IgG1 monoclonal antibody, in adults with advanced or metastatic solid tumors (Yamaguchi 2026; final human model). Three compartments (central = serum plus rapidly perfused tissues, peripheral = slowly perfused tissues, tumor), each holding free S-531011, free CCR8, the 1:1 complex and the 1:2 (one antibody, two CCR8) complex. S-531011 is dosed intravenously into the central compartment, eliminated from the central compartment only, and distributes to the peripheral and tumor compartments by first-order rate constants; k13 and k31 follow from Vc, the tumor volume, the tumor-to-serum concentration ratio Pdist and the tumor interstitial fraction fint. CCR8 is expressed on Tregs, conventional CD4 T cells and CD8 T cells, its baseline amount in each compartment is computed from cell counts and per-cell expression, it turns over at kdeg (the internalization rate of the complex, also applied to free CCR8 and both complexes), and the in vitro Kd is multiplied by 30 (refined against mouse tumor RO). Drug states are in mg, CCR8 and complex states in nmol. Outputs: serum S-531011 (Cc, ug/mL) and RO (%) in the central, peripheral and tumor compartments; tumor RO is the quantity the paper predicts. The between-subject variability is the variability the authors assumed for their simulations, not an estimate."
  reference <- paste(
    "Yamaguchi D, Nogami W, Sonoda Y, Morihara H, Katsube T.",
    "Semi-Mechanistic Modeling of S-531011, a Humanized Anti-CCR8 Monoclonal",
    "Antibody, for Prediction of CCR8 Receptor Occupancy in Human Tumor Tissues.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(9):e70332.",
    "doi:10.1002/psp4.70332.",
    "Parameter values from Table 1 and Tables S4-S5; equations from Table S1",
    "(ODEs) and Table S2 (derived constants); the authors' SimBiology model",
    "report (Data S1, model 'S-531011_PKRO_human_modified') supplies the",
    "CCR8 expression and cell-count inputs, the molecular weight and the",
    "Kd adjustment factor.",
    sep = " "
  )
  vignette <- "Yamaguchi_2026_s531011"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL (Cc); receptor occupancy in percent (ro_central, ro_peripheral, ro_tumor)"
  )

  # Per-compartment drug, receptor and complex species of the authors'
  # SimBiology model; the tumor-located states and the doubly bound (1:2)
  # complex have no canonical compartment names.
  paper_specific_compartments <- c(
    "drug_tumor",
    "target_tumor",
    "complex_tumor",
    "complex2",
    "complex2_peripheral1",
    "complex2_tumor"
  )

  compartmentData <- list(
    central = list(analyte = "S-531011 (free)", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "S-531011 (free)", units = "mg", specimen = "tissue", verified = TRUE),
    drug_tumor = list(analyte = "S-531011 (free)", units = "mg", specimen = "tumor", verified = TRUE),
    target = list(analyte = "CCR8 (free)", units = "nmol", specimen = "serum", verified = TRUE),
    target_peripheral1 = list(analyte = "CCR8 (free)", units = "nmol", specimen = "tissue", verified = TRUE),
    target_tumor = list(analyte = "CCR8 (free)", units = "nmol", specimen = "tumor", verified = TRUE),
    complex = list(analyte = "S-531011:CCR8 1:1 complex", units = "nmol", specimen = "serum", verified = TRUE),
    complex_peripheral1 = list(
      analyte = "S-531011:CCR8 1:1 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    complex_tumor = list(analyte = "S-531011:CCR8 1:1 complex", units = "nmol", specimen = "tumor", verified = TRUE),
    complex2 = list(analyte = "S-531011:CCR8 1:2 complex", units = "nmol", specimen = "serum", verified = TRUE),
    complex2_peripheral1 = list(
      analyte = "S-531011:CCR8 1:2 complex",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    complex2_tumor = list(analyte = "S-531011:CCR8 1:2 complex", units = "nmol", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 3,
    n_studies = 1,
    age_range = "18 years or older",
    disease_state = paste(
      "Locally advanced (locoregional recurrence not suitable for curative",
      "treatment) or metastatic solid tumors, ECOG performance status 0-1,",
      "without access to, intolerant of, or declining standard therapies",
      "(first-in-human Phase 1b/2 study, Part A-1 monotherapy dose escalation)."
    ),
    dose_range = "24 mg single IV infusion over 2 h (PK fit); 8-1600 mg IV every 3 weeks with 1-h infusions simulated.",
    regions = "Not reported.",
    notes = paste(
      "The human PK parameters (k12, k21, kel, Vc, Vp) were fitted by",
      "unweighted least squares to the MEAN serum concentrations of the three",
      "participants in the 24 mg group of Part A-1 (Phoenix WinNonlin 8.3).",
      "The receptor and binding parameters come from in vitro studies, the",
      "tumor distribution parameters from a CT26.WT tumor-bearing human",
      "CCR8-knock-in mouse study and the literature, and the CCR8 expression",
      "and T-cell counts from published human data. Kd was multiplied by 30",
      "after comparing a mouse version of the model (Yamaguchi_2026_s531011_mouse)",
      "with observed mouse tumor RO. No demographics are reported."
    )
  )

  ini({
    # --- Human two-compartment PK (WinNonlin fit to mean 24 mg data) ---
    lkel <- log(0.00479); label("First-order elimination rate constant of S-531011 from the central compartment (1/h)") # Table 1 kel 4.79e-3 1/h; Table S5 Human kel (CV 20.5%)
    lk12 <- log(0.00695); label("Central-to-peripheral transfer rate constant (1/h)") # Table 1 k12 6.95e-3 1/h; Table S5 Human k12 (CV 12.9%)
    lk21 <- log(0.00575); label("Peripheral-to-central transfer rate constant (1/h)") # Table 1 k21 5.75e-3 1/h; Table S5 Human k21 (CV 44.7%)
    lvc <- log(2.47); label("Volume of the central compartment (L)") # Table 1 Vcen 2.47 L; Table S5 Human Vcen (CV 1.2%)
    lvp <- log(2.99); label("Volume of the peripheral compartment (L)") # Table 1 Vper 2.99 L; Table S5 Human Vper (CV 49.7%)

    # --- Tumor distribution (literature and mouse values) ---
    lvtumor <- fixed(log(0.01)); label("Volume of the tumor compartment (L)") # Table 1 Vtum 0.01 L; Methods 'assumed to be a typical value [16]'
    lpdist <- fixed(log(0.61)); label("Tumor-to-serum concentration ratio of S-531011, Pdist (unitless)") # Results 'mean concentration ratio ... in the tumor to serum ... was 0.61'; Data S1 parameter 29 Pdist 0.61
    lfint <- fixed(log(0.55)); label("Fraction of interstitial space volume in the tumor, fint (unitless)") # Methods 'the value of fint (0.55) for mice was reported [17]'; Data S1 parameter 30 fint 0.55

    # --- Binding and turnover (in vitro) ---
    lkon <- fixed(log(11.85)); label("Association rate constant of S-531011 and CCR8 (L/nmol/h)") # Table 1 kon 11.85; Table S4 kon mean 11.85 L/nmol/hr
    lkd <- fixed(log(0.01862)); label("In vitro equilibrium dissociation constant of S-531011 for CCR8 (nmol/L)") # Table S4 Kd mean 18.62 pmol/L; Data S1 parameter 31 Kd 0.01862 nanomole/liter
    fadj <- fixed(30); label("Multiplier on the in vitro Kd, refined against mouse tumor RO (unitless)") # Results 'optimal value for Kd was ... 558.6 pmol/L, which was 30 times the initial value'; Data S1 parameter 70 Fadj 30
    lthalf_cpx <- fixed(log(1.71)); label("Half-life of cell-surface CCR8-bound S-531011, Thalf,CPX (h)") # Table S4 Thalf,CPX 1.71 hr; Data S1 parameter 47
    mw <- fixed(148653.5); label("Molecular weight of S-531011 (g/mol)") # Data S1 parameter 67 MW 148653.5 gram/mole

    # --- CCR8 expression per cell (molecules/cell) ---
    exp_treg_cen <- fixed(400); label("CCR8 expression per Treg, central compartment (molecules/cell)") # Data S1 parameter 32 exp_Treg_cen 400
    exp_tconv_cen <- fixed(575); label("CCR8 expression per conventional CD4 T cell, central compartment (molecules/cell)") # Data S1 parameter 33 exp_Tconv_cen 575
    exp_cd8_cen <- fixed(0); label("CCR8 expression per CD8 T cell, central compartment (molecules/cell)") # Data S1 parameter 34 exp_CD8_cen 0
    exp_treg_per <- fixed(1803); label("CCR8 expression per Treg, peripheral compartment (molecules/cell)") # Data S1 parameter 35 exp_Treg_per 1803
    exp_tconv_per <- fixed(812); label("CCR8 expression per conventional CD4 T cell, peripheral compartment (molecules/cell)") # Data S1 parameter 36 exp_Tconv_per 812
    exp_cd8_per <- fixed(1274); label("CCR8 expression per CD8 T cell, peripheral compartment (molecules/cell)") # Data S1 parameter 37 exp_CD8_per 1274
    exp_treg_tum <- fixed(3573); label("CCR8 expression per Treg, tumor compartment (molecules/cell)") # Data S1 parameter 38 exp_Treg_tum 3573
    exp_tconv_tum <- fixed(906); label("CCR8 expression per conventional CD4 T cell, tumor compartment (molecules/cell)") # Data S1 parameter 39 exp_Tconv_tum 906
    exp_cd8_tum <- fixed(1790); label("CCR8 expression per CD8 T cell, tumor compartment (molecules/cell)") # Data S1 parameter 40 exp_CD8_tum 1790

    # --- T-cell counts ---
    n_cd4_cen <- fixed(861); label("CD4 T-cell count in the central compartment (cells/uL)") # Data S1 parameter 41 n_CD4_cen 861
    n_cd8_cen <- fixed(505); label("CD8 T-cell count in the central compartment (cells/uL)") # Data S1 parameter 42 n_CD8_cen 505
    f_treg_cd4_cen <- fixed(0.07); label("Fraction of Tregs among CD4 T cells, central compartment (unitless)") # Data S1 parameter 43 f_Treg_CD4_cen 0.07
    n_treg_tum <- fixed(62.31); label("Treg count in the tumor compartment, also used for the peripheral compartment (cells/mm^2)") # Data S1 parameter 44 n_Treg_tum 62.31
    f_treg_cd4_tum <- fixed(0.189); label("Fraction of Tregs among CD4 T cells, tumor compartment (unitless)") # Data S1 parameter 45 f_Treg_CD4_tum 0.189
    r_cd4_cd8_tum <- fixed(1.56); label("Ratio of CD4 to CD8 T cells, tumor compartment (unitless)") # Data S1 parameter 46 r_CD4_CD8_tum 1.56
    funit <- fixed(0.707); label("Tissue volume per unit section area, funit (mm^3/mm^2)") # Data S1 parameter 62 Funit 0.707 millimeter^3/millimeter^2
    dens <- fixed(1); label("Volume conversion of the peripheral and tumor tissue, fvol (uL/mm^3)") # Data S1 parameter 63 DENS 1 microliter/millimeter^3; Methods 'fvol which was assumed to be 1 uL/mm3'

    # --- Between-subject variability assumed for the simulations ---
    # Methods (Simulation): '30% (coefficient of variation; CV) variabilities
    # for kel and Vcen, and 50% variabilities for kon, Pdist, fint, Kd,
    # Thalf,CPX, and Vtum'. Encoded as log-normal: log(1 + 0.30^2) = 0.0862
    # and log(1 + 0.50^2) = 0.2231. Assumed by the authors, not estimated.
    etalkel ~ fixed(0.0862) # Methods Simulation: 30 percent CV on kel
    etalvc ~ fixed(0.0862) # Methods Simulation: 30 percent CV on Vcen
    etalkon ~ fixed(0.2231) # Methods Simulation: 50 percent CV on kon
    etalpdist ~ fixed(0.2231) # Methods Simulation: 50 percent CV on Pdist
    etalfint ~ fixed(0.2231) # Methods Simulation: 50 percent CV on fint
    etalkd ~ fixed(0.2231) # Methods Simulation: 50 percent CV on Kd
    etalthalf_cpx ~ fixed(0.2231) # Methods Simulation: 50 percent CV on Thalf,CPX
    etalvtumor ~ fixed(0.2231) # Methods Simulation: 50 percent CV on Vtum

    # --- Residual error: none in the source (mean-data fit, simulation model) ---
    propSd <- fixed(0); label("Proportional residual error on serum S-531011 (fraction; not reported in the source)") # not reported
  })
  model({
    # Data S1 parameters 21 (Avogadro constant) and 69 (ln2 = 0.693, used for
    # Tdist and kdeg in the authors' code).
    avogadro <- 6.02214076e23
    ln2c <- 0.693

    # Individual parameters
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    vtumor <- exp(lvtumor + etalvtumor)
    pdist <- exp(lpdist + etalpdist)
    fint <- exp(lfint + etalfint)
    kon <- exp(lkon + etalkon)
    kd <- exp(lkd + etalkd)
    thalf_cpx <- exp(lthalf_cpx + etalthalf_cpx)

    # Derived constants (Table S2; Data S1 initial assignments 1-19)
    tdist <- ln2c / k12
    k13 <- (ln2c / tdist) * (pdist / fint) / ((pdist / fint) + vc / (fint * vtumor))
    k31 <- (ln2c / tdist) / (1 + (pdist / fint) * (fint * vtumor) / vc)
    koff <- kd * fadj * kon
    kdeg <- ln2c / thalf_cpx

    n_treg_cen <- f_treg_cd4_cen * n_cd4_cen
    n_tconv_cen <- (1 - f_treg_cd4_cen) * n_cd4_cen
    n_cd4_tum <- n_treg_tum / f_treg_cd4_tum
    n_tconv_tum <- n_cd4_tum - n_treg_tum
    n_cd8_tum <- n_cd4_tum / r_cd4_cd8_tum
    # Peripheral cell counts equal the tumor counts (Methods; Data S1 10-12)
    n_treg_per <- n_treg_tum
    n_tconv_per <- n_tconv_tum
    n_cd8_per <- n_cd8_tum

    # Baseline CCR8 amounts (nmol): molecules/uL x volume (L x 1e6 uL/L)
    # / Avogadro x 1e9 nmol/mol; tissue counts per mm^2 become per uL via
    # 1 / (funit x dens).
    ccr8_cen_int <- (exp_treg_cen * n_treg_cen + exp_tconv_cen * n_tconv_cen + exp_cd8_cen * n_cd8_cen) *
      (vc * 1e6) / avogadro * 1e9
    ccr8_per_int <- (exp_treg_per * n_treg_per + exp_tconv_per * n_tconv_per + exp_cd8_per * n_cd8_per) *
      (1 / (funit * dens)) * (vp * 1e6) / avogadro * 1e9
    ccr8_tum_int <- (exp_treg_tum * n_treg_tum + exp_tconv_tum * n_tconv_tum + exp_cd8_tum * n_cd8_tum) *
      (1 / (funit * dens)) * (vtumor * 1e6) / avogadro * 1e9

    # CCR8 is at steady state before dosing
    ksyn_cen <- kdeg * ccr8_cen_int
    ksyn_per <- kdeg * ccr8_per_int
    ksyn_tum <- kdeg * ccr8_tum_int

    # Free drug in nmol (drug states carry mg; 1 mg = 1e6 / MW nmol)
    nmol_per_mg <- 1e6 / mw
    drug_cen_nmol <- central * nmol_per_mg
    drug_per_nmol <- peripheral1 * nmol_per_mg
    drug_tum_nmol <- drug_tumor * nmol_per_mg

    # Binding fluxes (nmol/h): first site forms the 1:1 complex, a second
    # CCR8 then binds the 1:1 complex to form the 1:2 complex (Table S1)
    bind1_cen <- kon / vc * target * drug_cen_nmol - koff * complex
    bind2_cen <- kon / vc * target * complex - koff * complex2
    bind1_per <- kon / vp * target_peripheral1 * drug_per_nmol - koff * complex_peripheral1
    bind2_per <- kon / vp * target_peripheral1 * complex_peripheral1 - koff * complex2_peripheral1
    bind1_tum <- kon / vtumor * target_tumor * drug_tum_nmol - koff * complex_tumor
    bind2_tum <- kon / vtumor * target_tumor * complex_tumor - koff * complex2_tumor

    target(0) <- ccr8_cen_int
    target_peripheral1(0) <- ccr8_per_int
    target_tumor(0) <- ccr8_tum_int

    d/dt(central) <- -(k12 + k13 + kel) * central + k21 * peripheral1 + k31 * drug_tumor - bind1_cen / nmol_per_mg
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1 - bind1_per / nmol_per_mg
    d/dt(drug_tumor) <- k13 * central - k31 * drug_tumor - bind1_tum / nmol_per_mg

    d/dt(target) <- ksyn_cen - kdeg * target - bind1_cen - bind2_cen
    d/dt(complex) <- bind1_cen - bind2_cen - kdeg * complex
    d/dt(complex2) <- bind2_cen - kdeg * complex2

    d/dt(target_peripheral1) <- ksyn_per - kdeg * target_peripheral1 - bind1_per - bind2_per
    d/dt(complex_peripheral1) <- bind1_per - bind2_per - kdeg * complex_peripheral1
    d/dt(complex2_peripheral1) <- bind2_per - kdeg * complex2_peripheral1

    d/dt(target_tumor) <- ksyn_tum - kdeg * target_tumor - bind1_tum - bind2_tum
    d/dt(complex_tumor) <- bind1_tum - bind2_tum - kdeg * complex_tumor
    d/dt(complex2_tumor) <- bind2_tum - kdeg * complex2_tumor

    # Receptor occupancy (%), Methods equations; Data S1 repeated assignments 1-3
    ro_central <- 100 * (1 - target / ccr8_cen_int)
    ro_peripheral <- 100 * (1 - target_peripheral1 / ccr8_per_int)
    ro_tumor <- 100 * (1 - target_tumor / ccr8_tum_int)

    # Serum S-531011 (mg/L = ug/mL); Data S1 repeated assignment 4
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
