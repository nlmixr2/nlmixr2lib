Yamaguchi_2026_s531011_mouse <- function() {
  description <- "Preclinical (mouse; CT26.WT tumor-bearing human CCR8 knock-in BALB/c). QSP. Semi-mechanistic PK / receptor-occupancy (RO) model of S-531011, a humanized anti-CCR8 IgG1 monoclonal antibody (Yamaguchi 2026; refined mouse model). Same structure as Yamaguchi_2026_s531011_human: central, peripheral and tumor compartments, each holding free S-531011, free CCR8, the 1:1 complex and the 1:2 (one antibody, two CCR8) complex; elimination from the central compartment only; k13 and k31 from Vc, the tumor volume, the tumor-to-serum concentration ratio Pdist and the interstitial fraction fint; CCR8 turnover at kdeg applied to free CCR8 and both complexes; the in vitro Kd multiplied by 30, the factor the authors selected by fitting this model to the observed 24-h tumor RO. CCR8 is expressed on CD4 T cells only and its baseline amount in each compartment is taken directly from Table 1. Drug states are in mg, CCR8 and complex states in nmol. Outputs: serum S-531011 (Cc, ug/mL) and RO (%) in the central, peripheral and tumor compartments. The only between-animal variability is the 50% CV on Kd the authors used for their RO simulations."
  reference <- paste(
    "Yamaguchi D, Nogami W, Sonoda Y, Morihara H, Katsube T.",
    "Semi-Mechanistic Modeling of S-531011, a Humanized Anti-CCR8 Monoclonal",
    "Antibody, for Prediction of CCR8 Receptor Occupancy in Human Tumor Tissues.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(9):e70332.",
    "doi:10.1002/psp4.70332.",
    "Parameter values from Table 1 (Mouse column) and Tables S4-S5; equations",
    "from Table S1 (ODEs), Table S2 (k13, k31, koff, ksyn) and Table S3;",
    "molecular weight and Kd adjustment factor from the authors' SimBiology",
    "model report (Data S1).",
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
    species = "mouse (BALB/c human CCR8 knock-in, CT26.WT tumor-bearing)",
    n_subjects = 4,
    n_studies = 1,
    age_range = "10 weeks at implantation; dosed 8 days after CT26.WT implantation",
    sex_female_pct = 100,
    disease_state = "Subcutaneous CT26.WT colon carcinoma syngeneic tumor.",
    dose_range = "0.05, 0.5 and 5 mg/kg single IV dose (PK fit to the 5 mg/kg group).",
    regions = "Japan (Shionogi).",
    notes = paste(
      "The mouse PK parameters (k12, k21, kel, Vc, Vp) were fitted by",
      "weighted least squares (weight = -2) to the MEAN serum concentrations",
      "of four mice given 5 mg/kg (Phoenix WinNonlin 8.3). A separate study",
      "(10 mice per group, five per day) measured CCR8 RO on Tregs in spleen",
      "and tumor-infiltrating lymphocytes at Day 1 and Day 7; only the Day 1",
      "values were used, adjusted so the 5 mg/kg mean is 100%. The tumor",
      "volume is a literature baseline CT26 tumor volume; the CCR8 amounts are",
      "computed from the mouse CD4 T-cell CCR8 expression (not tabulated) and",
      "the human cell counts."
    )
  )

  ini({
    # --- Mouse two-compartment PK (WinNonlin fit to mean 5 mg/kg data) ---
    lkel <- log(0.0324); label("First-order elimination rate constant of S-531011 from the central compartment (1/h)") # Table 1 Mouse kel 0.0324 1/h; Table S5 Mouse kel (CV 13.2%)
    lk12 <- log(0.136); label("Central-to-peripheral transfer rate constant (1/h)") # Table 1 Mouse k12 0.136 1/h; Table S5 Mouse k12 (CV 38.3%)
    lk21 <- log(0.0907); label("Peripheral-to-central transfer rate constant (1/h)") # Table 1 Mouse k21 0.0907 1/h; Table S5 Mouse k21 (CV 34.2%)
    lvc <- log(1.14e-3); label("Volume of the central compartment (L)") # Table 1 Mouse Vcen 1.14e-3 L; Table S5 Mouse Vcen (CV 12.5%)
    lvp <- log(1.71e-3); label("Volume of the peripheral compartment (L)") # Table 1 Mouse Vper 1.71e-3 L; Table S5 Mouse Vper (CV 17.2%)

    # --- Tumor distribution ---
    lvtumor <- fixed(log(1.70e-4)); label("Volume of the tumor compartment (L)") # Table 1 Mouse Vtum 1.70e-4 L; Methods 'baseline tumor volume in mice was obtained from the literature [23]'
    lpdist <- fixed(log(0.61)); label("Tumor-to-serum concentration ratio of S-531011, Pdist (unitless)") # Results 'mean concentration ratio ... in the tumor to serum ... was 0.61'; Methods mouse Pdist same as human
    lfint <- fixed(log(0.55)); label("Fraction of interstitial space volume in the tumor, fint (unitless)") # Methods 'the value of fint (0.55) for mice was reported [17]'

    # --- Binding and turnover (in vitro; same as human) ---
    lkon <- fixed(log(11.85)); label("Association rate constant of S-531011 and CCR8 (L/nmol/h)") # Table 1 Mouse kon 11.85; Table S4 kon mean 11.85 L/nmol/hr
    lkd <- fixed(log(0.01862)); label("In vitro equilibrium dissociation constant of S-531011 for CCR8 (nmol/L)") # Table S4 Kd mean 18.62 pmol/L; Methods mouse kon and koff same as human
    fadj <- 30; label("Multiplier on the in vitro Kd, selected against observed mouse tumor RO (unitless)") # Results 'optimal value for Kd was ... 558.6 pmol/L, which was 30 times the initial value' (Figure S5)
    lthalf_cpx <- fixed(log(1.71)); label("Half-life of cell-surface CCR8-bound S-531011, Thalf,CPX (h)") # Table S4 Thalf,CPX 1.71 hr; Table 1 Mouse kdeg 0.4053
    mw <- fixed(148653.5); label("Molecular weight of S-531011 (g/mol)") # Data S1 parameter 67 MW 148653.5 gram/mole

    # --- Baseline CCR8 amounts (Table 1, Mouse column) ---
    rbase_target <- fixed(8.476e-7); label("Baseline free CCR8 in the central compartment (nmol)") # Table 1 Mouse CCR8cen 8.476e-7 nmol
    rbase_target_peripheral1 <- fixed(5.791e-8); label("Baseline free CCR8 in the peripheral compartment (nmol)") # Table 1 Mouse CCR8per 5.791e-8 nmol
    rbase_target_tumor <- fixed(4.679e-7); label("Baseline free CCR8 in the tumor compartment (nmol)") # Table 1 Mouse CCR8tum 4.679e-7 nmol

    # --- Between-animal variability for the RO simulations ---
    # Figure S5 and S6 captions: RO 'simulated based on the mouse PK/RO model
    # using optimized Kd with 50%CV (n=500)'. Encoded as log-normal:
    # log(1 + 0.50^2) = 0.2231. Assumed by the authors, not estimated.
    etalkd ~ fixed(0.2231) # Figure S5 caption: 50 percent CV on Kd

    # --- Residual error: none in the source (mean-data fit, simulation model) ---
    propSd <- fixed(0); label("Proportional residual error on serum S-531011 (fraction; not reported in the source)") # not reported
  })
  model({
    # Data S1 parameter 69: ln2 = 0.693, used for Tdist and kdeg
    ln2c <- 0.693

    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc <- exp(lvc)
    vp <- exp(lvp)
    vtumor <- exp(lvtumor)
    pdist <- exp(lpdist)
    fint <- exp(lfint)
    kon <- exp(lkon)
    kd <- exp(lkd + etalkd)
    thalf_cpx <- exp(lthalf_cpx)

    # Derived constants (Table S2, applied to mouse values per Methods)
    tdist <- ln2c / k12
    k13 <- (ln2c / tdist) * (pdist / fint) / ((pdist / fint) + vc / (fint * vtumor))
    k31 <- (ln2c / tdist) / (1 + (pdist / fint) * (fint * vtumor) / vc)
    koff <- kd * fadj * kon
    kdeg <- ln2c / thalf_cpx

    # CCR8 is at steady state before dosing
    ksyn_cen <- kdeg * rbase_target
    ksyn_per <- kdeg * rbase_target_peripheral1
    ksyn_tum <- kdeg * rbase_target_tumor

    # Free drug in nmol (drug states carry mg; 1 mg = 1e6 / MW nmol)
    nmol_per_mg <- 1e6 / mw
    drug_cen_nmol <- central * nmol_per_mg
    drug_per_nmol <- peripheral1 * nmol_per_mg
    drug_tum_nmol <- drug_tumor * nmol_per_mg

    # Binding fluxes (nmol/h), Table S1
    bind1_cen <- kon / vc * target * drug_cen_nmol - koff * complex
    bind2_cen <- kon / vc * target * complex - koff * complex2
    bind1_per <- kon / vp * target_peripheral1 * drug_per_nmol - koff * complex_peripheral1
    bind2_per <- kon / vp * target_peripheral1 * complex_peripheral1 - koff * complex2_peripheral1
    bind1_tum <- kon / vtumor * target_tumor * drug_tum_nmol - koff * complex_tumor
    bind2_tum <- kon / vtumor * target_tumor * complex_tumor - koff * complex2_tumor

    target(0) <- rbase_target
    target_peripheral1(0) <- rbase_target_peripheral1
    target_tumor(0) <- rbase_target_tumor

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

    # Receptor occupancy (%), Methods equations
    ro_central <- 100 * (1 - target / rbase_target)
    ro_peripheral <- 100 * (1 - target_peripheral1 / rbase_target_peripheral1)
    ro_tumor <- 100 * (1 - target_tumor / rbase_target_tumor)

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
