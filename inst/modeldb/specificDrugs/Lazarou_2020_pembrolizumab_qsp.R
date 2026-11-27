Lazarou_2020_pembrolizumab_qsp <- function() {
  description <- paste(
    "QSP. Minimal Cancer Immunity Cycle (CIC) model calibrated for non-small",
    "cell lung carcinoma (NSCLC) treated with the anti-PD-1 antibody",
    "pembrolizumab. A two-compartment anti-PD-1 PK model with an",
    "injection-site depot drives PD-1 checkpoint blockade in a",
    "tumour-microenvironment (TME) / lymph-node immune model: logistic tumour",
    "growth towards a TME capacity that the tumour itself recruits, neoantigen",
    "release by T-cell killing, antigen uptake by antigen-presenting cells",
    "(APCs) and their migration to the lymph node, APC-driven cytotoxic T-cell",
    "proliferation in the lymph node, T-cell trafficking to the TME limited by",
    "a maximal TME fraction, and PD-1-mediated T-cell exhaustion that the",
    "antibody relieves by receptor occupancy. The baseline TME composition",
    "(leukocyte, cytotoxic T-cell and APC fractions) comes from the Immune",
    "Landscape of Cancer transcriptomic deconvolution of lung adenocarcinoma",
    "and squamous cell biopsies. Nine ODE states (three PK, six",
    "tumour/immune) plus a cumulative antigen pool; deterministic, with no",
    "IIV or residual error. The authors' virtual population draws six",
    "rate constants uniformly between 0.5x and 1x their nominal values,",
    "the initial tumour volume uniformly between 10 and 10000 mm^3, and the",
    "cytotoxic T-cell and APC TME fractions from log-normal distributions;",
    "ini() holds the nominal values.",
    sep = " "
  )
  reference <- paste(
    "Lazarou G, Chelliah V, Small BG, Walker M, van der Graaf PH, Kierzek AM.",
    "Integration of Omics Data Sources to Inform Mechanistic Modeling of",
    "Immune-Oncology Therapies: A Tutorial for Clinical Pharmacologists.",
    "Clin Pharmacol Ther. 2020;107(4):858-870. doi:10.1002/cpt.1786.",
    "Model structure, equations and every parameter value are taken from the",
    "authors' deposited model code (Supporting Information 'Model Code',",
    "CIC_v2.2.6: CSV documentation plus the R and MATLAB exports).",
    sep = " "
  )
  vignette <- "Lazarou_2020_pembrolizumab_qsp"

  # Cancer-immunity-cycle cell populations. `tumor` is the canonical tumour
  # state (here a tumour-cell count); the remaining six are mechanistic
  # pools of this minimal CIC model that do not generalise to a canonical
  # compartment family.
  paper_specific_compartments <- c(
    "capacity",
    "antigens",
    "apc_tme",
    "apc_ln",
    "tc_tme",
    "tc_ln"
  )

  units <- list(
    time = "day",
    dosing = paste(
      "nmol of pembrolizumab, dosed into depot (the model's injection-site",
      "compartment). Convert a mg dose with nmol = mg / 146000 * 1e6; the",
      "paper's 10 mg/kg dose for a 70 kg patient is 700 mg = 4794.521 nmol.",
      sep = " "
    ),
    concentration = paste(
      "ug/mL for the pembrolizumab plasma concentration Cc. Tumour volume",
      "tumor_vol is in mm^3; cell states are cell counts.",
      sep = " "
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pembrolizumab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pembrolizumab", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pembrolizumab", units = "nmol", specimen = "tissue", verified = TRUE),
    tumor = list(analyte = "tumour cells", units = "cells", specimen = "tumor", verified = TRUE),
    capacity = list(analyte = "tumour-microenvironment capacity", units = "cells", specimen = "tumor", verified = TRUE),
    antigens = list(
      analyte = "tumour neoantigen (killed tumour cells)",
      units = "cells",
      specimen = "tumor",
      verified = TRUE
    ),
    apc_tme = list(analyte = "antigen-presenting cells", units = "cells", specimen = "tumor", verified = TRUE),
    apc_ln = list(analyte = "antigen-presenting cells", units = "cells", specimen = "tissue", verified = TRUE),
    tc_tme = list(analyte = "cytotoxic CD8 T cells", units = "cells", specimen = "tumor", verified = TRUE),
    tc_ln = list(
      analyte = "neoantigen-specific cytotoxic T cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUM_VOL = list(
      description = "Tumour volume at the start of the simulation (start of therapy)",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Initial condition only: sets the starting tumour-cell, TME-capacity,",
        "cytotoxic T-cell and APC counts (each its maximal TME percentage of",
        "TUM_VOL / vcell) and is the reference of the percent change from",
        "baseline. Source parameter Initial_tumor_volume, documented in L",
        "(nominal 0.01259253909 L = 12592.53909 mm^3, the day-0 tumour size of",
        "one of the four calibration patients). In the virtual trial it is",
        "drawn uniformly between 0.01 mL and 10 mL (10 to 10000 mm^3;",
        "Assignments.csv A1).",
        sep = " "
      ),
      source_name = "Initial_tumor_volume"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 4L,
    n_studies = 1L,
    disease_state = "Non-small cell lung carcinoma (KEYNOTE-001), anti-PD-1 treated",
    dose_range = "Pembrolizumab 10 mg/kg every 3 weeks (19 doses over 399 days), body weight fixed at 70 kg",
    regions = "Not reported",
    notes = paste(
      "The tumour-growth, tumour-kill and PD-1-checkpoint parameters were",
      "calibrated by hand against the individual tumour-volume trajectories of",
      "four KEYNOTE-001 NSCLC patients treated with pembrolizumab 10 mg/kg Q3W",
      "(three responders, one non-responder). The baseline TME composition is",
      "derived from Immune Landscape of Cancer data on 1156 lung",
      "adenocarcinoma and squamous cell carcinoma biopsies; lymph-node T-cell",
      "numbers come from healthy-donor flow cytometry and T-cell repertoire",
      "data. The paper's virtual trial simulates 489 virtual patients (the",
      "KEYNOTE-001 enrolment) at 70 kg. The model was built for illustration",
      "and is explicitly not qualified for drug-development decisions.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Pembrolizumab PK (Parameters.csv; the authors fitted these to a
    # Simcyp simulation). Two-compartment PK with an injection-site depot
    # that loses the non-absorbed fraction (reactions abs / elIS / el /
    # distr).
    # ------------------------------------------------------------------
    lka <- fixed(log(2)); label("Absorption rate constant from the injection site (1/day)")  # Parameters.csv kabs = 2 1/day
    lfdepot <- fixed(log(0.9815751)); label("Fraction of the injection-site dose absorbed (fraction)")  # Parameters.csv Fabs = 0.9815751
    lkel <- fixed(log(0.07342571)); label("Elimination rate constant from the central compartment (1/day)")  # Parameters.csv kel = 0.07342571 1/day
    lk12 <- fixed(log(0.4)); label("Central to peripheral transfer rate constant (1/day)")  # Parameters.csv k12 = 0.4 1/day
    lk21 <- fixed(log(0.43)); label("Peripheral to central transfer rate constant (1/day)")  # Parameters.csv k21 = 0.43 1/day
    mw <- fixed(146000); label("Molecular weight of pembrolizumab (g/mol)")  # Parameters.csv MWantiPD1 = 146000 g/mol
    vblood <- fixed(2.75); label("Volume converting the central amount to a plasma concentration (L)")  # Parameters.csv Vblood = 2.75 L
    kd_pd1 <- fixed(0.076595741); label("Anti-PD-1 dissociation constant, expressed as a central amount (nmol)")  # Parameters.csv antiPD1_KD = 0.076595741 nano*mol

    # ------------------------------------------------------------------
    # Tumour growth, TME capacity and killing
    # ------------------------------------------------------------------
    kg_tumor <- fixed(0.4); label("Maximal (logistic) tumour growth rate constant (1/day)")  # Parameters.csv mu = 0.4 1/day; virtual trial uniform(0.2, 0.4) (A16)
    kv1 <- fixed(0.01); label("Rate constant of TME capacity recruitment (1/day)")  # Parameters.csv kv1 = 0.01 1/day; virtual trial uniform(0.005, 0.01) (A17)
    kdeath_tumor <- fixed(0.0231); label("Basal tumour-cell death rate constant (1/day)")  # Parameters.csv beta_tumor = 0.0231 1/day; virtual trial uniform(0.01155, 0.0231) (A18)
    kkill <- fixed(0.3); label("Maximal rate constant of tumour killing by cytotoxic T cells (1/day)")  # Parameters.csv kkill = 0.3 1/day; virtual trial uniform(0.15, 0.3) (A19)
    pcc50_kill <- fixed(3.45); label("TME cytotoxic T-cell percentage at half-maximal killing (percent)")  # Parameters.csv Kkill = 3.45 (half saturation of kill function)

    # ------------------------------------------------------------------
    # PD-1 checkpoint (T-cell exhaustion)
    # ------------------------------------------------------------------
    kexh_pd1 <- fixed(300); label("Maximal rate constant of cytotoxic T-cell exhaustion via PD-1 (1/day)")  # Parameters.csv kPD1 = 300 1/day; virtual trial uniform(150, 300) (A20)
    ptum50_pd1 <- fixed(43.5); label("TME tumour-cell percentage at half-maximal PD-1 checkpoint action (percent)")  # Parameters.csv KPD1 = 43.5

    # ------------------------------------------------------------------
    # Baseline TME composition (Immune Landscape of Cancer, LUAD + LUSC)
    # ------------------------------------------------------------------
    leuk_frac <- fixed(0.14); label("Leukocyte fraction of the TME (fraction)")  # Parameters.csv Leuk_frac = 0.14; Supplement p. 1 'fraction of leukocytes in the TME is 0.14'
    tc_frac <- fixed(0.074); label("Relative cytotoxic T-cell fraction of the TME leukocytes (fraction)")  # Parameters.csv Tc_frac = 0.074; virtual trial log-normal (A3: mean_Tc = -2.436, SD_Tc = 0.62335)
    apc_frac <- fixed(0.066); label("Relative APC (dendritic cell + M1 macrophage) fraction of the TME leukocytes (fraction)")  # Parameters.csv APC_frac = 0.066; virtual trial log-normal (A5: mean_APC = -2.58098, SD_APC = 0.69689)
    vcell <- fixed(1e-6); label("Volume of one cell (mm^3/cell)")  # Parameters.csv Vcell = 1E-06 mm^3/cell

    # ------------------------------------------------------------------
    # Antigen-presenting cells
    # ------------------------------------------------------------------
    ksyn_apc_tme <- fixed(108005.973888231); label("Maximal rate of APC development in the TME (cells/mL/day)")  # Parameters.csv ksyn_APC_tme = 108005.973888231 cell/mL/day
    apc_syn_scale <- fixed(5000); label("Constant multiplier of the APC development rate (unitless)")  # Reactions.csv v4 rate law 'ksyn_APC_tme * 5000 * ...', commented 'blood volume = 5 L'
    kdeath_apc_tme <- fixed(0.0466); label("APC death rate constant in the TME (1/day)")  # Parameters.csv beta_APC_tme = 0.0466 1/day
    kdeath_apc_ln <- fixed(0.0466); label("APC death rate constant in the lymph node (1/day)")  # Parameters.csv beta_APC_ln = 0.0466 1/day
    qlymph <- fixed(7.03); label("Lymph flow setting the APC migration rate to the lymph node (L/day)")  # Parameters.csv L = 7.03 L/day; Assignments.csv A14 delta_APC_ln = L / TME
    apc_ln0 <- fixed(1000000); label("Initial number of APCs in the lymph node (cells)")  # Species.csv APC_ln initial state = 1000000 cell

    # ------------------------------------------------------------------
    # Cytotoxic T cells
    # ------------------------------------------------------------------
    kpro_tc_ln <- fixed(1.5); label("Maximal APC-driven cytotoxic T-cell proliferation rate in the lymph node (1/day)")  # Parameters.csv Pro_Tc_ln = 1.5 1/day; virtual trial uniform(0.75, 1.5) (A21)
    kdeath_tc_ln <- fixed(0.0014); label("Cytotoxic T-cell death rate constant in the lymph node (1/day)")  # Parameters.csv beta_Tc_ln = 0.0014 1/day
    kdeath_tc_tme <- fixed(0.0014); label("Cytotoxic T-cell death rate constant in the TME (1/day)")  # Parameters.csv beta_Tc_tme = 0.0014 1/day
    qco <- fixed(4730); label("Cardiac output setting the T-cell trafficking rate to the TME (L/day)")  # Parameters.csv Q = 4730 L/day; Assignments.csv A15 delta_Tc_tme = Q / P
    vplasma <- fixed(5); label("Plasma volume setting the T-cell trafficking rate to the TME (L)")  # Parameters.csv P = 5 L
    tc_ln0 <- fixed(137554.400017073); label("Initial number of neoantigen-specific cytotoxic T cells in the lymph node (cells)")  # Species.csv Tc_ln initial state = 137554.400017073 cell
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters
    # ------------------------------------------------------------------
    ka <- exp(lka)
    fdepot <- exp(lfdepot)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)

    # Maximal TME percentages (Assignments.csv A22-A27). The source's
    # Tc_frac_ratio = Tc_frac / APC_frac and APC_frac_ratio = 1 reduce to
    # these forms; the three percentages sum to 100.
    tumor_max <- 100 * (1 - leuk_frac)
    tc_max <- 100 * leuk_frac * tc_frac / (tc_frac + apc_frac)
    apc_max <- 100 * leuk_frac * apc_frac / (tc_frac + apc_frac)

    # T-cell trafficking rate constant (A15)
    delta_tc_tme <- qco / vplasma

    # ------------------------------------------------------------------
    # 2. Initial conditions: the tumour starts infiltrated to its maximal
    #    leukocyte fractions (A4, A8, A10), with the TME capacity equal to
    #    the tumour-cell count.
    # ------------------------------------------------------------------
    tumor(0) <- tumor_max / 100 * TUM_VOL / vcell
    capacity(0) <- tumor_max / 100 * TUM_VOL / vcell
    tc_tme(0) <- tc_max / 100 * TUM_VOL / vcell
    apc_tme(0) <- apc_max / 100 * TUM_VOL / vcell
    apc_ln(0) <- apc_ln0
    tc_ln(0) <- tc_ln0

    # ------------------------------------------------------------------
    # 3. TME composition (Compartments.csv TME volume; A2, A6, A13)
    # ------------------------------------------------------------------
    ncell_tme <- tumor + apc_tme + tc_tme
    tumor_vol <- ncell_tme * vcell
    pct_tumor <- 100 * tumor / ncell_tme
    pct_apc <- 100 * apc_tme / ncell_tme
    pct_tc <- 100 * tc_tme / ncell_tme
    pct_change <- 100 * (tumor_vol - TUM_VOL) / TUM_VOL

    # APC migration rate constant L / TME (A14); TME converted mm^3 -> L
    delta_apc_ln <- qlymph / (tumor_vol * 1e-6)

    # ------------------------------------------------------------------
    # 4. Rate laws (Reactions.csv)
    # ------------------------------------------------------------------
    # PD-1 occupancy by the antibody, from the central amount (v14)
    occ_pd1 <- central / (central + kd_pd1)

    r_growth <- kg_tumor * tumor * (1 - tumor / capacity)
    r_kill <- kkill * tumor * pct_tc / (pct_tc + pcc50_kill)
    r_death_tumor <- kdeath_tumor * tumor
    r_capacity <- kv1 * tumor * tumor_max / pct_tumor
    # APC development; TME volume converted mm^3 -> mL to match ksyn_apc_tme
    r_apc_syn <- ksyn_apc_tme * apc_syn_scale * (1 - pct_apc / apc_max) * tumor_vol / 1000
    r_apc_mig <- delta_apc_ln * apc_tme * antigens / ncell_tme
    r_tc_pro <- kpro_tc_ln * apc_ln
    r_tc_mig <- delta_tc_tme * tc_ln * (1 - pct_tc / tc_max)
    r_exh <- kexh_pd1 * (1 - occ_pd1) * tc_tme * pct_tumor / (pct_tumor + ptum50_pd1)

    # ------------------------------------------------------------------
    # 5. ODEs
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- fdepot * ka * depot - k12 * central + k21 * peripheral1 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    d/dt(tumor) <- r_growth - r_kill - r_death_tumor
    d/dt(capacity) <- r_capacity
    d/dt(antigens) <- r_kill - r_apc_mig
    d/dt(apc_tme) <- r_apc_syn - kdeath_apc_tme * apc_tme - r_apc_mig
    d/dt(apc_ln) <- r_apc_mig - kdeath_apc_ln * apc_ln
    d/dt(tc_ln) <- r_tc_pro - r_tc_mig - kdeath_tc_ln * tc_ln
    d/dt(tc_tme) <- r_tc_mig - kdeath_tc_tme * tc_tme - r_exh

    # ------------------------------------------------------------------
    # 6. Outputs: plasma concentration (A9; nmol * g/mol / L * 1e-6 =
    #    ug/mL) and tumour volume (mm^3)
    # ------------------------------------------------------------------
    Cc <- central * mw / vblood * 1e-6
  })
}
