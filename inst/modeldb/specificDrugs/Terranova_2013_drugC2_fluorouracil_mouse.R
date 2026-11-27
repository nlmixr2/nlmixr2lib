Terranova_2013_drugC2_fluorouracil_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; subcutaneous HT29 human colon adenocarcinoma xenograft, Terranova 2013 Experiment d).",
    "Two-drug tumor-growth-inhibition model for Drug C2 (a polo-like kinase inhibitor coded by the authors; no INN disclosed) combined with 5-fluorouracil (5-FU):",
    "the mouse PK models of both drugs (two-compartment first-order oral absorption for Drug C2, two-compartment IV bolus for 5-fluorouracil; held at their single-agent",
    "estimates) drive a Simeoni exponential-then-linear tumor-growth model with a four-by-four",
    "grid of damaged-cell transit compartments, one axis per drug, plus an interaction term",
    "gamma * Ca * Cb that moves proliferating cells directly into the doubly damaged state",
    "(gamma > 0 synergy, gamma < 0 antagonism). Typical-value fit to mean tumor weight with",
    "no between-animal variability."
  )
  reference <- paste(
    "Terranova N, Germani M, Del Bene F, Magni P. A predictive pharmacokinetic-pharmacodynamic",
    "model of tumor growth kinetics in xenograft mice after administration of anticancer agents",
    "given in combination. Cancer Chemother Pharmacol. 2013;72(2):471-482.",
    "doi:10.1007/s00280-013-2208-8. PK parameters from Supplementary Table 1; tumor and",
    "single-agent drug parameters from Supplementary Table 2; interaction parameter from the",
    "Results; equations from Equations 1-3 and Figure 1 (right panel). Growth model from",
    "Simeoni M, et al. Cancer Res. 2004;64(3):1094-1101. doi:10.1158/0008-5472.CAN-03-2524."
  )
  vignette <- "Terranova_2013_combination_tgi"
  dosing <- c("depot_drugc2", "central_5fu")

  units <- list(
    time = "day",
    dosing = "umol/kg (Drug C2 into depot_drugc2, 5-fluorouracil into central_5fu)",
    concentration = "uM"
  )

  covariateData <- list()

  # Terranova 2013 Equation 2 grid: damaged_cells_a<i>_b<j> is paper state x_ij, a
  # tumor cell that has progressed i stages along the Drug C2 (drug a) chain and j
  # stages along the 5-fluorouracil (drug b) chain; x_00 is cycling_cells.
  paper_specific_compartments <- c(
    "damaged_cells_a0_b1",
    "damaged_cells_a0_b2",
    "damaged_cells_a0_b3",
    "damaged_cells_a1_b0",
    "damaged_cells_a1_b1",
    "damaged_cells_a1_b2",
    "damaged_cells_a1_b3",
    "damaged_cells_a2_b0",
    "damaged_cells_a2_b1",
    "damaged_cells_a2_b2",
    "damaged_cells_a2_b3",
    "damaged_cells_a3_b0",
    "damaged_cells_a3_b1",
    "damaged_cells_a3_b2",
    "damaged_cells_a3_b3"
  )

  compartmentData <- list(
    depot_drugc2 = list(analyte = "Drug C2", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central_drugc2 = list(analyte = "Drug C2", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_drugc2 = list(analyte = "Drug C2", units = "umol/kg", specimen = "tissue", verified = TRUE),
    central_5fu = list(analyte = "5-Fluorouracil", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_5fu = list(analyte = "5-Fluorouracil", units = "umol/kg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a0_b1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a0_b2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a0_b3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a1_b0 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a1_b1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a1_b2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a1_b3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a2_b0 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a2_b1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a2_b2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a2_b3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a3_b0 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a3_b1 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a3_b2 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells_a3_b3 = list(analyte = "drug-damaged tumor cells", units = "g", specimen = "tumor", verified = TRUE)
  )

  population <- list(
    species = "mouse (athymic nude-nu bearing HT29 xenografts)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5 weeks at purchase",
    weight_range = "25-35 g",
    sex_female_pct = 0,
    disease_state = "subcutaneous HT29 human colon adenocarcinoma xenograft, randomized at a palpable tumor of about 100-300 mm^3",
    dose_range = "Drug C2 45 mg/kg (84.5 umol/kg, arm d1) or 60 mg/kg (113 umol/kg, arm d2) orally on days 10, 11, 12, 14, 15 and 16; 5-fluorouracil 50 mg/kg (384 umol/kg) IV on days 9, 13 and 17 (Table 1)",
    regions = "preclinical (Nerviano Medical Sciences, Italy)",
    notes = paste(
      "Terranova 2013 Experiment d: male Balb athymic nude-nu mice, 5 weeks old, 25-35 g.",
      "Control, two single-agent and the combination arms, usually eight animals per",
      "group (exact counts not reported). Time zero is tumor implantation; tumor weight",
      "was computed from caliper length and width as length * width^2 / 2 at unit density.",
      "The seven tumor and single-agent parameters were estimated in Rocchetti 2009 (Experiment 3) from the control and single-agent arms;",
      "the interaction parameter was then estimated on the combination arm with them held",
      "at those values, by weighted least squares on the arm means."
    )
  )

  ini({
    # --- Drug C2 mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    # Table 1 reports no bioavailability for the oral drug, so V1 is read as the
    # apparent volume V/F (bioavailability 1 in the simulation; see the vignette).
    lvc_drugc2 <- fixed(log(2.13)); label("Drug C2 central volume of distribution (L/kg)") # Supplementary Table 1 'Drug C2': V1 = 2.13 l/kg
    lkel_drugc2 <- fixed(log(49.2)); label("Drug C2 elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Drug C2': K10 = 49.2 1/day
    lk12_drugc2 <- fixed(log(141)); label("Drug C2 central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Drug C2': K12 = 141 1/day
    lk21_drugc2 <- fixed(log(10.4)); label("Drug C2 peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Drug C2': K21 = 10.4 1/day
    lka_drugc2 <- fixed(log(18.8)); label("Drug C2 first-order oral absorption rate constant (1/day)") # Supplementary Table 1 'Drug C2': K01 = 18.8 1/day

    # --- 5-Fluorouracil mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_5fu <- fixed(log(0.714)); label("5-Fluorouracil central volume of distribution (L/kg)") # Supplementary Table 1 '5-fluorouracil': V1 = 0.714 l/kg
    lkel_5fu <- fixed(log(151)); label("5-Fluorouracil elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 '5-fluorouracil': K10 = 151 1/day
    lk12_5fu <- fixed(log(5.62)); label("5-Fluorouracil central-to-peripheral rate constant (1/day)") # Supplementary Table 1 '5-fluorouracil': K12 = 5.62 1/day
    lk21_5fu <- fixed(log(2.30)); label("5-Fluorouracil peripheral-to-central rate constant (1/day)") # Supplementary Table 1 '5-fluorouracil': K21 = 2.30 1/day

    # --- Tumor growth and single-agent effects: Supplementary Table 2, Experiment d ---
    ltumorExpGrowth <- log(0.0903); label("Exponential-phase tumor growth rate lambda0 (1/day)") # Supplementary Table 2 Experiment d: lambda0 = 0.0903 1/day
    ltumorLinGrowth <- log(0.0459); label("Linear-phase tumor growth rate lambda1 (g/day)") # Supplementary Table 2 Experiment d: lambda1 = 0.0459 g/day
    # Terranova 2013 keeps the Simeoni 2004 growth equation (Methods, Pharmacodynamic
    # model) without printing psi; Simeoni 2004 sets psi = 20 for a sharp switch.
    psi <- fixed(20); label("Exponential-to-linear growth switch sharpness psi (unitless)") # Simeoni 2004: psi = 20 (not printed in Terranova 2013)
    lrbase_tumor <- log(0.0884); label("Tumor weight at implantation w0 (g)") # Supplementary Table 2 Experiment d: w0 = 0.0884 g
    ldamageTransit_drugc2 <- log(2.46); label("Drug C2 damaged-cell transit rate constant k1a (1/day)") # Supplementary Table 2 Experiment d: k1A = 2.46 1/day
    ldrugSlope_drugc2 <- log(0.241); label("Drug C2 potency k2a (1/(uM*day))") # Supplementary Table 2 Experiment d: k2A = 0.241 1/(uM*day)
    ldamageTransit_5fu <- log(0.612); label("5-Fluorouracil damaged-cell transit rate constant k1b (1/day)") # Supplementary Table 2 Experiment d: k1B = 0.612 1/day
    ldrugSlope_5fu <- log(0.0837); label("5-Fluorouracil potency k2b (1/(uM*day))") # Supplementary Table 2 Experiment d: k2B = 0.0837 1/(uM*day)

    # --- Drug-effect interaction: Results ---
    # The Results print the unit as uM^-2 day^-2; Equation 1 multiplies gamma by
    # Ca * Cb * x00 to give a rate in g/day, so the unit is 1/(uM^2*day). Signed, so
    # it is carried on the linear scale (positive = synergy, negative = antagonism).
    synergy <- 0.131; label("Drug-effect interaction parameter gamma (1/(uM^2*day))") # Results, arm d1: gamma = 0.131 (CV 42%)
    # Arm d2 was fitted independently: gamma = 0.157 1/(uM^2*day) (CV 36%). Simulate it with
    # ini(synergy = 0.157).

    # Weighted least squares with weights 1/y^2 (Methods, Data analysis), i.e. a
    # proportional error whose magnitude is not reported; held at zero.
    propSd <- fixed(0); label("Proportional residual error on tumor weight (fraction)")
  })

  model({
    # --- PK (Supplementary Table 1 micro-constants); doses in umol/kg over V1 in L/kg give uM ---
    vc_drugc2 <- exp(lvc_drugc2)
    kel_drugc2 <- exp(lkel_drugc2)
    k12_drugc2 <- exp(lk12_drugc2)
    k21_drugc2 <- exp(lk21_drugc2)
    ka_drugc2 <- exp(lka_drugc2)
    d/dt(depot_drugc2) <- -ka_drugc2 * depot_drugc2
    d/dt(central_drugc2) <- ka_drugc2 * depot_drugc2 - (kel_drugc2 + k12_drugc2) * central_drugc2 + k21_drugc2 * peripheral1_drugc2
    d/dt(peripheral1_drugc2) <- k12_drugc2 * central_drugc2 - k21_drugc2 * peripheral1_drugc2
    Cc_drugc2 <- central_drugc2 / vc_drugc2

    vc_5fu <- exp(lvc_5fu)
    kel_5fu <- exp(lkel_5fu)
    k12_5fu <- exp(lk12_5fu)
    k21_5fu <- exp(lk21_5fu)
    d/dt(central_5fu) <- -(kel_5fu + k12_5fu) * central_5fu + k21_5fu * peripheral1_5fu
    d/dt(peripheral1_5fu) <- k12_5fu * central_5fu - k21_5fu * peripheral1_5fu
    Cc_5fu <- central_5fu / vc_5fu

    # --- Tumor growth inhibition (Terranova 2013 Equations 1-3, Figure 1 right panel) ---
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    rbase_tumor <- exp(lrbase_tumor)
    damageTransit_drugc2 <- exp(ldamageTransit_drugc2)
    drugSlope_drugc2 <- exp(ldrugSlope_drugc2)
    damageTransit_5fu <- exp(ldamageTransit_5fu)
    drugSlope_5fu <- exp(ldrugSlope_5fu)

    killRate_drugc2 <- drugSlope_drugc2 * Cc_drugc2
    killRate_5fu <- drugSlope_5fu * Cc_5fu
    killRate_interaction <- synergy * Cc_drugc2 * Cc_5fu

    # Equation 3: total tumor weight w(t) is the sum of all 16 states.
    tumor_size <- cycling_cells +
      damaged_cells_a0_b1 +
      damaged_cells_a0_b2 +
      damaged_cells_a0_b3 +
      damaged_cells_a1_b0 +
      damaged_cells_a1_b1 +
      damaged_cells_a1_b2 +
      damaged_cells_a1_b3 +
      damaged_cells_a2_b0 +
      damaged_cells_a2_b1 +
      damaged_cells_a2_b2 +
      damaged_cells_a2_b3 +
      damaged_cells_a3_b0 +
      damaged_cells_a3_b1 +
      damaged_cells_a3_b2 +
      damaged_cells_a3_b3

    # Simeoni growth f_p(w) acting on the proliferating cells.
    growth <- tumorExpGrowth * cycling_cells / (1 + (tumorExpGrowth / tumorLinGrowth * tumor_size)^psi)^(1 / psi)

    # Equation 1; the interaction flux v_11 = gamma * Ca * Cb * x00 leaves x00 and enters x11.
    d/dt(cycling_cells) <- growth - (killRate_drugc2 + killRate_5fu + killRate_interaction) * cycling_cells
    # Equation 2 with the flows of Figure 1 (right panel): a damaged cell hit by the
    # other drug leaves its state, so every off-origin state loses k2 * C of the drug
    # it has not yet been damaged by (see the vignette for this reading).
    d/dt(damaged_cells_a0_b1) <- killRate_5fu * cycling_cells - (killRate_drugc2 + damageTransit_5fu) * damaged_cells_a0_b1
    d/dt(damaged_cells_a0_b2) <- damageTransit_5fu * damaged_cells_a0_b1 - (killRate_drugc2 + damageTransit_5fu) * damaged_cells_a0_b2
    d/dt(damaged_cells_a0_b3) <- damageTransit_5fu * damaged_cells_a0_b2 - (killRate_drugc2 + damageTransit_5fu) * damaged_cells_a0_b3
    d/dt(damaged_cells_a1_b0) <- killRate_drugc2 * cycling_cells - (damageTransit_drugc2 + killRate_5fu) * damaged_cells_a1_b0
    d/dt(damaged_cells_a1_b1) <- killRate_drugc2 * damaged_cells_a0_b1 + killRate_5fu * damaged_cells_a1_b0 + killRate_interaction * cycling_cells - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a1_b1
    d/dt(damaged_cells_a1_b2) <- killRate_drugc2 * damaged_cells_a0_b2 + damageTransit_5fu * damaged_cells_a1_b1 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a1_b2
    d/dt(damaged_cells_a1_b3) <- killRate_drugc2 * damaged_cells_a0_b3 + damageTransit_5fu * damaged_cells_a1_b2 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a1_b3
    d/dt(damaged_cells_a2_b0) <- damageTransit_drugc2 * damaged_cells_a1_b0 - (damageTransit_drugc2 + killRate_5fu) * damaged_cells_a2_b0
    d/dt(damaged_cells_a2_b1) <- damageTransit_drugc2 * damaged_cells_a1_b1 + killRate_5fu * damaged_cells_a2_b0 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a2_b1
    d/dt(damaged_cells_a2_b2) <- damageTransit_drugc2 * damaged_cells_a1_b2 + damageTransit_5fu * damaged_cells_a2_b1 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a2_b2
    d/dt(damaged_cells_a2_b3) <- damageTransit_drugc2 * damaged_cells_a1_b3 + damageTransit_5fu * damaged_cells_a2_b2 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a2_b3
    d/dt(damaged_cells_a3_b0) <- damageTransit_drugc2 * damaged_cells_a2_b0 - (damageTransit_drugc2 + killRate_5fu) * damaged_cells_a3_b0
    d/dt(damaged_cells_a3_b1) <- damageTransit_drugc2 * damaged_cells_a2_b1 + killRate_5fu * damaged_cells_a3_b0 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a3_b1
    d/dt(damaged_cells_a3_b2) <- damageTransit_drugc2 * damaged_cells_a2_b2 + damageTransit_5fu * damaged_cells_a3_b1 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a3_b2
    d/dt(damaged_cells_a3_b3) <- damageTransit_drugc2 * damaged_cells_a2_b3 + damageTransit_5fu * damaged_cells_a3_b2 - (damageTransit_drugc2 + damageTransit_5fu) * damaged_cells_a3_b3

    cycling_cells(0) <- rbase_tumor

    tumor_size ~ prop(propSd)
  })
}
