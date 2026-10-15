Smith_2022_pah_competitive_metabolism <- function() {
  description <- paste(
    "In vitro (pooled human liver microsomes). Mechanistic competitive-inhibition",
    "metabolism model for the two polycyclic aromatic hydrocarbons (PAHs)",
    "benzo[a]pyrene (BaP) and dibenzo[def,p]chrysene (DBC). Each parent is",
    "cleared by a high-affinity/low-capacity saturable (Michaelis-Menten) enzyme",
    "plus a non-saturable high-capacity/low-affinity linear intrinsic-clearance",
    "arm (the 'Michaelis-Menten clearance' model, Eq. 2); the co-incubated PAH",
    "competitively inhibits the saturable arm by raising its apparent Km (Eq. 4,",
    "the BIC-best competitive model). The two coupled states are the substrate",
    "concentrations in a pooled-human-liver-microsome incubation; the metabolism",
    "rate equations carry the per-mg-microsomal-protein in vitro estimates",
    "(Tables 3-4). Deterministic: the publication reports no IIV and no residual",
    "error. This file is the in vitro metabolism layer only; the paper's",
    "whole-body PBPK interaction model inherits its physiology (compartment",
    "volumes, blood flows, partition coefficients, absorption rates) from Pande",
    "2022 (DBC) and Crowell 2011 (BaP), which are not reproduced here.",
    sep = " "
  )
  reference <- paste(
    "Smith JN, Gaither KA, Pande P. Competitive Metabolism of Polycyclic Aromatic",
    "Hydrocarbons (PAHs): An Assessment Using In Vitro Metabolism and",
    "Physiologically Based Pharmacokinetic (PBPK) Modeling. Int J Environ Res",
    "Public Health. 2022;19(14):8266. doi:10.3390/ijerph19148266. PMCID:",
    "PMC9323266. In vitro Michaelis-Menten-clearance parameters: Table 3;",
    "competitive-inhibition constants (Ki): Table 4; metabolism-rate equations:",
    "Eq. 2 (baseline, no inhibitor) and Eq. 4 (competitive inhibition raising the",
    "apparent Km1 of the saturable arm).",
    sep = " "
  )
  vignette <- "Smith_2022_pah_competitive_metabolism"
  units <- list(
    time = "min",
    dosing = "uM (substrate concentration placed into the incubation at time 0)",
    concentration = "uM"
  )

  # The two parent PAHs are placed into the incubation medium at time 0, so the
  # dosing targets are the substrate states rather than depot or central.
  dosing <- c("bap", "dbc")

  # Paper-mechanistic in vitro states: the two substrate concentrations in the
  # microsomal incubation. Neither is an instance of a canonical PK compartment
  # role, so they are whitelisted here.
  paper_specific_compartments <- c("bap", "dbc")

  compartmentData <- list(
    bap = list(
      analyte = "benzo[a]pyrene (BaP) parent substrate",
      units = "uM",
      specimen = "administration site",
      verified = TRUE
    ),
    dbc = list(
      analyte = "dibenzo[def,p]chrysene (DBC) parent substrate",
      units = "uM",
      specimen = "administration site",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (pooled human liver microsomes)",
    n_subjects = 200L,
    n_studies = 1L,
    age_range = "19-78 years (3 donors 10-18 years)",
    sex_female_pct = 50,
    disease_state = "not applicable (pooled human liver microsomes)",
    notes = paste(
      "Pooled human liver microsomes (Sekisui Xenotech) contributed by 200",
      "donors with an equal male:female ratio, predominantly 19-78 years (3",
      "donors 10-18 years). Incubations (Methods 2.2): 2.0 mg/mL microsomal",
      "protein, 0.1 M phosphate buffer pH 7.4, 3 mM MgCl2, excess (1.5 mM)",
      "NADPH, 37 C. BaP 0.05-2.5 uM incubated 0-30 min; DBC 0.025-1 uM incubated",
      "0-60 min. Competitive-inhibition assays held one PAH at a fixed substrate",
      "concentration (BaP 0.14-0.18 uM; DBC 0.17 uM) and co-incubated the other",
      "PAH across 0.1-10 uM. All kinetic parameters are per mg microsomal",
      "protein. Supermix-10 (a 10-PAH environmental mixture) was measured as an",
      "additional inhibitor (Ki 0.75 uM on BaP, 0.63 uM on DBC; Table 4) but the",
      "authors built no dynamical Supermix-10 model, so it is not encoded here."
    )
  )

  ini({
    # ---- In vitro Michaelis-Menten-clearance parameters (Table 3) ---------
    # Fit to the 'Michaelis-Menten clearance' model (Eq. 2), the BIC-best
    # baseline model for both substrates (Table 2). Vmax1/Km1 is the saturable
    # high-affinity/low-capacity enzyme; Clint2 is the non-saturable linear
    # high-capacity/low-affinity arm that did not saturate over the tested
    # concentrations.
    vmax_bap <- fixed(0.0063)
    label("BaP Vmax1, high-affinity saturable enzyme (nmol/min/mg microsomal protein)") # Table 3, BaP Vmax1 0.0063 (95% CI 0.0044-0.0083)
    km_bap <- fixed(0.088)
    label("BaP Km1, high-affinity saturable enzyme (uM)") # Table 3, BaP Km1 0.088 (0.044-0.15)
    clint2_bap <- fixed(0.0012)
    label("BaP Clint2, non-saturable linear intrinsic clearance (mL/min/mg microsomal protein)") # Table 3, BaP Clint2 0.0012 (1.7e-7 to 0.0028)
    vmax_dbc <- fixed(0.00090)
    label("DBC Vmax1, high-affinity saturable enzyme (nmol/min/mg microsomal protein)") # Table 3, DBC Vmax1 0.00090 (0.00044-0.0023)
    km_dbc <- fixed(0.060)
    label("DBC Km1, high-affinity saturable enzyme (uM)") # Table 3, DBC Km1 0.060 (0.014-0.22)
    clint2_dbc <- fixed(0.0017)
    label("DBC Clint2, non-saturable linear intrinsic clearance (mL/min/mg microsomal protein)") # Table 3, DBC Clint2 0.0017 (6.5e-7 to 0.0028)

    # ---- Competitive-inhibition constants, Ki (Table 4) -------------------
    # Ki modifies the Michaelis-Menten constant Km1 of the saturable arm
    # (Eq. 4); the apparent Km1 becomes Km1 * (1 + [I] / Ki).
    ki_dbc_bap <- fixed(0.44)
    label("Ki for DBC inhibiting BaP metabolism (uM)") # Table 4, substrate BaP / inhibitor DBC Ki 0.44 (0.36-0.54)
    ki_bap_dbc <- fixed(0.061)
    label("Ki for BaP inhibiting DBC metabolism (uM)") # Table 4, substrate DBC / inhibitor BaP Ki 0.061 (0.041-0.12)

    # ---- In vitro incubation condition -----------------------------------
    cmic <- fixed(2.0)
    label("Microsomal protein concentration in the incubation (mg/mL)") # Methods 2.2: microsomes 2.0 mg/mL
  })

  model({
    # Per-mg-microsomal-protein metabolism rates (nmol/min/mg). Eq. 4: the
    # co-incubated PAH competitively raises the apparent Km1 of the saturable
    # arm; the non-saturable Clint2 arm is uninhibited. With the other PAH at
    # zero these collapse to the baseline Eq. 2.
    rate_bap <- vmax_bap * bap / (km_bap * (1 + dbc / ki_dbc_bap) + bap) +
      clint2_bap * bap
    rate_dbc <- vmax_dbc * dbc / (km_dbc * (1 + bap / ki_bap_dbc) + dbc) +
      clint2_dbc * dbc

    # Intrinsic clearance of the first-order (low-concentration) phase,
    # Clint1 = Vmax1/Km1 (Methods 2.5; Table 3: BaP 0.072, DBC 0.015 mL/min/mg).
    clint1_bap <- vmax_bap / km_bap
    clint1_dbc <- vmax_dbc / km_dbc

    # Substrate disappearance in the closed incubation. The per-mg rate
    # (nmol/min/mg) times the microsomal protein concentration (mg/mL) gives a
    # concentration rate (nmol/mL/min = uM/min).
    d/dt(bap) <- -cmic * rate_bap
    d/dt(dbc) <- -cmic * rate_dbc
  })
}
