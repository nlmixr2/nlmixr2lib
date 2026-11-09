Choi_2024_cd123nkce_human <- function() {
  description <- "Translational mechanistic synapse-driven PK/PD model of the trifunctional NKp46/CD16a x CD123 natural killer cell engager (CD123-NKCE) in patients with acute myeloid leukemia, used to project the first-in-human starting dose. The monkey model structure is carried over with allometrically scaled linear two-compartment PK, human or AML system parameters (NK and CD123+ blast counts, CD123 density, binding affinities, blast growth and NK turnover), no anti-drug antibodies, and the tolerance feedback replaced by logistic CD123+ cell growth to a maximum capacity. The NKp46-drug-CD123 synapse drives Hill-type CD123+ cell killing and NK-cell redistribution out of blood; extra CD123 in non-blood tissue acts as a target-mediated drug sink (Choi 2024 Appendix S1 section 3 and Table S2). All parameters are fixed; no variability."
  reference <- paste(
    "Choi SL, Valente D, Virone-Oddos A, Mauriac C.",
    "Developing a mechanistic translational PK/PD model for a trifunctional NK cell",
    "engager to predict the first-in-human dose for acute myeloid leukemia.",
    "Clin Transl Sci. 2024;17(1):e13689. doi:10.1111/cts.13689",
    sep = " "
  )
  vignette <- "Choi_2024_cd123nkce"
  # Time is hours. Doses are entered in ug/kg (the paper's dose unit); f()
  # converts them to nmol/kg with the molecular weight, the same conversion
  # the Monolix script applies to its dosing inputs (Table S2
  # 'p = 1/tdur/MW*1e3' over a 0.0001-h input, i.e. an IV bolus). central and
  # peripheral1 hold nmol/kg; every binding state holds nmol/L; the cell
  # states hold cells/uL.
  units <- list(time = "h", dosing = "ug/kg", concentration = "nmol/L")

  # The drug-receptor binary complexes on blood cells and the two circulating
  # cell pools have no canonical compartment name; they are named after the
  # drug-<antigen>-<site> family of drug_cd3_central and the tcell_central
  # cell-pool family.
  paper_specific_compartments <- c(
    "drug_nkp46_central",
    "drug_cd123_central",
    "nkcell_central",
    "cd123cell_central"
  )

  compartmentData <- list(
    central = list(analyte = "CD123-NKCE", units = "nmol/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "CD123-NKCE", units = "nmol/kg", specimen = "tissue", verified = TRUE),
    drug_nkp46_central = list(
      analyte = "CD123-NKCE-NKp46 binary complex",
      units = "nmol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug_cd123_central = list(
      analyte = "CD123-NKCE-CD123 binary complex",
      units = "nmol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    trimer = list(
      analyte = "NKp46-CD123-NKCE-CD123 synapse",
      units = "nmol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    target = list(
      analyte = "CD123 expressed in non-blood tissue",
      units = "nmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex = list(
      analyte = "CD123-NKCE bound to non-blood CD123",
      units = "nmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    nkcell_central = list(analyte = "NKp46+ NK cells", units = "cells/uL", specimen = "whole blood", verified = TRUE),
    cd123cell_central = list(
      analyte = "CD123+ AML cells",
      units = "cells/uL",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0,
    disease_state = "Acute myeloid leukemia (projected first-in-human population; no clinical data fitted)",
    dose_range = "Simulated IV 10-100 ug/kg twice weekly (intrapatient escalation 10/30/100 ug/kg, or a flat 100 ug/kg)",
    notes = paste(
      "Translational projection from the monkey PK/PD model (Choi 2024 Results, 'Translational AML",
      "PK/PD model'). Linear PK allometrically scaled from monkey with exponents 0.81 (CL), 0.57 (Q),",
      "1.04 (V1) and 1.07 (V2). System parameters: baseline NKp46+ NK cells 54 cells/uL (86% of blood",
      "NK cells), CD123+ AML cells 1668 cells/uL, CD123 density 7000 per cell from 24 AML patients",
      "(ref 26), NK half-life 14 days (ref 25), blast doubling time 24 h, maximum CD123+ capacity 95%",
      "above baseline. Drug-specific parameters (kmax, EC50, gamma, kint, SlpNK, kinTP, alpha, beta)",
      "assumed equal to monkey."
    )
  )

  ini({
    # Linear PK, per kg body weight, scaled from monkey (Table 1, AML column)
    lcl <- fixed(log(0.18))
    label("Nonspecific clearance per kg body weight (mL/h/kg)") # Table 1 AML: CL = 0.18 (footnote m: scaled from NHP)
    lq <- fixed(log(0.37))
    label("Intercompartmental clearance per kg body weight (mL/h/kg)") # Table 1 AML: Q = 0.37 (footnote m: scaled from NHP)
    lvc <- fixed(log(0.042))
    label("Central volume of distribution per kg body weight (L/kg)") # Table 1 AML: V1 = 0.042 (footnote m: scaled from NHP)
    lvp <- fixed(log(0.031))
    label("Peripheral volume of distribution per kg body weight (L/kg)") # Table 1 AML: V2 = 0.031 (footnote m: scaled from NHP)
    mw <- fixed(149000)
    label("Molecular weight of CD123-NKCE (g/mol)") # Table S2: MW = 1,49,000 g/mol

    # Target binding (Table 1 AML)
    nkp46_receptors <- fixed(4600)
    label("NKp46 receptors per NK cell (receptors/cell)") # Table 1 AML: RTdE = 4600 (internal data)
    cd123_receptors <- fixed(7000)
    label("CD123 receptors per CD123+ AML cell (receptors/cell)") # Table 1 AML: RTdT = 7000 (footnote e: ref 26)
    kd_nkp46 <- fixed(16.6)
    label("Equilibrium dissociation constant of CD123-NKCE for human NKp46 (nmol/L)") # Table 1 AML: KD1 = 16.6 (internal data)
    kd_cd123 <- fixed(0.404)
    label("Equilibrium dissociation constant of CD123-NKCE for human CD123 (nmol/L)") # Table 1 AML: KD2 = 0.404 (internal data)
    koff_nkp46 <- fixed(0.001)
    label("Dissociation rate constant of CD123-NKCE from NKp46 (1/s)") # Table 1 AML: koff1 = 0.001 1/s (footnote o: ref 15)
    koff_cd123 <- fixed(0.001)
    label("Dissociation rate constant of CD123-NKCE from CD123 (1/s)") # Table 1 AML: koff2 = 0.001 1/s (footnote o: ref 15)
    lkint <- fixed(log(0.15))
    label("Internalisation rate constant of drug-CD123 complexes (1/h)") # Table 1 AML: kint = 0.15 (footnote n: estimated in NHP)

    # Extra CD123 expressed in non-blood tissue (Table 1 AML; Table S2)
    lksyn <- fixed(log(0.21))
    label("Zero-order production rate of non-blood CD123 (nmol/L/h)") # Table 1 AML: kinTP = 0.21 (footnote n: estimated in NHP)
    ltarget_scale <- fixed(log(19.44))
    label("Scale (intercept) of non-blood CD123 versus blood CD123 at baseline (nmol/L)") # Table 1 AML: alpha = 19.44 (footnote n: estimated in NHP)
    ltarget_exp <- fixed(log(0.21))
    label("Power exponent of non-blood CD123 versus blood CD123 at baseline (unitless)") # Table 1 AML: beta = 0.21 (footnote n: estimated in NHP)

    # Synapse-driven killing (Table 1 AML; from the in vitro model)
    lkmax <- fixed(log(0.332))
    label("Maximal CD123+ cell killing rate constant of the synapse (1/h)") # Table 1 AML: kmax = 0.332 (footnote g: in vitro estimate)
    lec50 <- fixed(log(1.46e-9))
    label("Synapse concentration at half-maximal killing rate (nmol/L)") # Table 1 AML: EC50 = 1.46E-9 (footnote g: in vitro estimate)
    lhill <- fixed(log(1.04))
    label("Hill coefficient of synapse-driven killing (unitless)") # Table 1 AML: gamma = 1.04 (footnote g: in vitro estimate)

    # CD123+ AML cell growth (Table 1 AML)
    lkg <- fixed(log(0.0288))
    label("First-order growth rate constant of CD123+ AML cells (1/h)") # Table 1 AML: kg = 0.0288 (footnote i: doubling time 24 h)
    tcmax <- fixed(3252)
    label("Maximum capacity of CD123+ cells in blood (cells/uL)") # Table 1 AML: TCmax = 3252 (footnote q: 95% above baseline)

    # NK-cell redistribution (Table 1 AML)
    lkout <- fixed(log(0.002))
    label("NKp46+ NK cell turnover rate constant (1/h)") # Table 1 AML: kout = 0.002 (footnote k: ref 25, 14-day half-life)
    lslope_nk <- fixed(log(5.72))
    label("Slope of synapse-driven NK-cell redistribution (per 1e9 x synapse in nmol/L)") # Table 1 AML: SlpNK = 5.72 (footnote n: estimated in NHP)

    # Baseline cell counts (Table 1 AML)
    rbase_nk <- fixed(54)
    label("Baseline NKp46+ NK cell count in blood (cells/uL)") # Table 1 AML: EF0 = 54 (internal data)
    rbase_cd123 <- fixed(1668)
    label("Baseline CD123+ AML cell count in blood (cells/uL)") # Table 1 AML: TC0 = 1668 (internal data)
  })

  model({
    # Avogadro's number per nmol (Table S2: NA = 6.022e14 #/nmol)
    navo <- 6.022e14

    # PK; CL and Q in mL/h/kg -> L/h/kg (Table S2 'k01 = (CL/1000)/V')
    cl <- exp(lcl) / 1000
    q <- exp(lq) / 1000
    vc <- exp(lvc)
    vp <- exp(lvp)

    kint <- exp(lkint)
    ksyn <- exp(lksyn)
    target_scale <- exp(ltarget_scale)
    target_exp <- exp(ltarget_exp)
    kmax <- exp(lkmax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    kg <- exp(lkg)
    kout <- exp(lkout)
    slope_nk <- exp(lslope_nk)

    # koff in 1/s -> 1/h; kon = koff / KD (Table S2)
    koff1 <- koff_nkp46 * 3600
    koff2 <- koff_cd123 * 3600
    kon1 <- koff1 / kd_nkp46
    kon2 <- koff2 / kd_cd123

    # Baseline blood CD123 (nmol/L) and non-blood CD123 pool (Table S2 'TR0',
    # 'TP0 = Alpha * TR0^Beta', 'koutTP = kinTP/TP0')
    blood_cd123_0 <- rbase_cd123 * cd123_receptors * 1e6 / navo
    target0 <- target_scale * blood_cd123_0^target_exp
    kdeg_target <- ksyn / target0

    # Free drug and free receptors on circulating cells (nmol/L; Table S2 'E', 'T')
    cfree <- central / vc
    efree <- nkcell_central * nkp46_receptors * 1e6 / navo - drug_nkp46_central - trimer
    tfree <- cd123cell_central * cd123_receptors * 1e6 / navo - drug_cd123_central - trimer

    # Binding flux on free drug (nmol/L/h; Table S2 'ddt_B' binding terms)
    bind <- -kon1 * cfree * efree + koff1 * drug_nkp46_central -
      kon2 * cfree * tfree + koff2 * drug_cd123_central -
      kon2 * cfree * target + koff2 * complex

    # PK (Table S2 'ddt_B', 'ddt_AP'), amounts in nmol/kg
    d/dt(central) <- -cl * cfree - q * cfree + q * peripheral1 / vp + bind * vc
    d/dt(peripheral1) <- q * cfree - q * peripheral1 / vp

    # Binding on circulating cells (Table S2 'ddt_BE', 'ddt_TB', 'ddt_TBE')
    d/dt(drug_nkp46_central) <- kon1 * cfree * efree - koff1 * drug_nkp46_central -
      kon2 * tfree * drug_nkp46_central + koff2 * trimer
    d/dt(drug_cd123_central) <- kon2 * tfree * cfree - koff2 * drug_cd123_central -
      kon1 * drug_cd123_central * efree + koff1 * trimer - kint * drug_cd123_central
    d/dt(trimer) <- kon2 * tfree * drug_nkp46_central - koff2 * trimer +
      kon1 * drug_cd123_central * efree - koff1 * trimer

    # Non-blood CD123 turnover and binding (Table S2 'ddt_TP', 'ddt_TAP'; the
    # script's undefined 'kinttb' in ddt_TAP is kint, per Appendix S1 Eq 17)
    d/dt(target) <- ksyn - kdeg_target * target - kon2 * cfree * target + koff2 * complex
    d/dt(complex) <- kon2 * cfree * target - koff2 * complex - kint * complex
    target(0) <- target0

    # Synapse floored at 0 (Table S2 'TBE2 = max(0, TBE)')
    synapse <- max(0, trimer)
    kkill <- kmax * synapse^hill / (ec50^hill + synapse^hill)

    # CD123+ AML cells: logistic growth to TCmax minus synapse killing (Table S2 'ddt_TA')
    d/dt(cd123cell_central) <- kg * cd123cell_central * (1 - cd123cell_central / tcmax) -
      kkill * cd123cell_central
    cd123cell_central(0) <- rbase_cd123

    # NK-cell redistribution (Table S2 'ddt_EF'; exponent nNK = 1 per Appendix S1 Eq 19)
    d/dt(nkcell_central) <- kout * rbase_nk - kout * (1 + slope_nk * 1e9 * synapse) * nkcell_central
    nkcell_central(0) <- rbase_nk

    # Dose unit conversion ug/kg -> nmol/kg (Table S2 'p = 1/MW*1e3')
    f(central) <- 1e3 / mw

    Cc <- cfree
    cd123_pctcfb <- (cd123cell_central - rbase_cd123) / rbase_cd123 * 100
  })
}
