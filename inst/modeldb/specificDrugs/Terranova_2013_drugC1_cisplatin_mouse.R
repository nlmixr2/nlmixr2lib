Terranova_2013_drugC1_cisplatin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; subcutaneous A2780 human ovarian carcinoma xenograft, Terranova 2013 Experiment b).",
    "Two-drug tumor-growth-inhibition model for Drug C1 (an Aurora kinase inhibitor coded by the authors; no INN disclosed) combined with cisplatin:",
    "the mouse PK models of both drugs (two-compartment IV bolus for Drug C1, two-compartment IV bolus for cisplatin; held at their single-agent",
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
  dosing <- c("central_drugc1", "central_cisplatin")

  units <- list(
    time = "day",
    dosing = "umol/kg (Drug C1 into central_drugc1, cisplatin into central_cisplatin)",
    concentration = "uM"
  )

  covariateData <- list()

  # Terranova 2013 Equation 2 grid: damaged_cells_a<i>_b<j> is paper state x_ij, a
  # tumor cell that has progressed i stages along the Drug C1 (drug a) chain and j
  # stages along the cisplatin (drug b) chain; x_00 is cycling_cells.
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
    central_drugc1 = list(analyte = "Drug C1", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_drugc1 = list(analyte = "Drug C1", units = "umol/kg", specimen = "tissue", verified = TRUE),
    central_cisplatin = list(analyte = "Cisplatin", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_cisplatin = list(analyte = "Cisplatin", units = "umol/kg", specimen = "tissue", verified = TRUE),
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
    species = "mouse (athymic nude-nu bearing A2780 xenografts)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5 weeks at purchase",
    weight_range = "25-35 g",
    sex_female_pct = 0,
    disease_state = "subcutaneous A2780 human ovarian carcinoma xenograft, randomized at a palpable tumor of about 100-300 mm^3",
    dose_range = "arm b1: Drug C1 30 mg/kg (58.7 umol/kg) IV twice daily on days 8-12 then cisplatin 8 mg/kg (26.7 umol/kg) IV on day 13; arm b2 (predicted, not fitted): cisplatin on day 8 then Drug C1 twice daily on days 9-13 (Table 1)",
    regions = "preclinical (Nerviano Medical Sciences, Italy)",
    notes = paste(
      "Terranova 2013 Experiment b: male Balb athymic nude-nu mice, 5 weeks old, 25-35 g.",
      "Control, two single-agent and the combination arms, usually eight animals per",
      "group (exact counts not reported). Time zero is tumor implantation; tumor weight",
      "was computed from caliper length and width as length * width^2 / 2 at unit density.",
      "The seven tumor and single-agent parameters were estimated in this paper by fitting the control, single-agent Drug C1 and single-agent cisplatin arms simultaneously (RMSE 0.26 g);",
      "the interaction parameter was then estimated on the combination arm with them held",
      "at those values, by weighted least squares on the arm means."
    )
  )

  ini({
    # --- Drug C1 mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_drugc1 <- fixed(log(1.420)); label("Drug C1 central volume of distribution (L/kg)") # Supplementary Table 1 'Drug C1': V1 = 1.420 l/kg
    lkel_drugc1 <- fixed(log(28.1)); label("Drug C1 elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Drug C1': K10 = 28.1 1/day
    lk12_drugc1 <- fixed(log(4.94)); label("Drug C1 central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Drug C1': K12 = 4.94 1/day
    lk21_drugc1 <- fixed(log(5.58)); label("Drug C1 peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Drug C1': K21 = 5.58 1/day

    # --- Cisplatin mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_cisplatin <- fixed(log(0.395)); label("Cisplatin central volume of distribution (L/kg)") # Supplementary Table 1 'Cisplatin': V1 = 0.395 l/kg
    lkel_cisplatin <- fixed(log(106)); label("Cisplatin elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Cisplatin': K10 = 106 1/day
    lk12_cisplatin <- fixed(log(109)); label("Cisplatin central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Cisplatin': K12 = 109 1/day
    lk21_cisplatin <- fixed(log(47.3)); label("Cisplatin peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Cisplatin': K21 = 47.3 1/day

    # --- Tumor growth and single-agent effects: Supplementary Table 2, Experiment b ---
    ltumorExpGrowth <- log(0.326); label("Exponential-phase tumor growth rate lambda0 (1/day)") # Supplementary Table 2 Experiment b: lambda0 = 0.326 1/day (7%)
    ltumorLinGrowth <- log(0.286); label("Linear-phase tumor growth rate lambda1 (g/day)") # Supplementary Table 2 Experiment b: lambda1 = 0.286 g/day (8%)
    # Terranova 2013 keeps the Simeoni 2004 growth equation (Methods, Pharmacodynamic
    # model) without printing psi; Simeoni 2004 sets psi = 20 for a sharp switch.
    psi <- fixed(20); label("Exponential-to-linear growth switch sharpness psi (unitless)") # Simeoni 2004: psi = 20 (not printed in Terranova 2013)
    lrbase_tumor <- log(0.0211); label("Tumor weight at implantation w0 (g)") # Supplementary Table 2 Experiment b: w0 = 0.0211 g (19%)
    ldamageTransit_drugc1 <- log(0.608); label("Drug C1 damaged-cell transit rate constant k1a (1/day)") # Supplementary Table 2 Experiment b: k1A = 0.608 1/day (15%)
    ldrugSlope_drugc1 <- log(0.186); label("Drug C1 potency k2a (1/(uM*day))") # Supplementary Table 2 Experiment b: k2A = 0.186 1/(uM*day) (9%)
    ldamageTransit_cisplatin <- log(0.155); label("Cisplatin damaged-cell transit rate constant k1b (1/day)") # Supplementary Table 2 Experiment b: k1B = 0.155 1/day (>100%)
    ldrugSlope_cisplatin <- log(1.89); label("Cisplatin potency k2b (1/(uM*day))") # Supplementary Table 2 Experiment b: k2B = 1.89 1/(uM*day) (48%)

    # --- Drug-effect interaction: Results ---
    # The Results print the unit as uM^-2 day^-2; Equation 1 multiplies gamma by
    # Ca * Cb * x00 to give a rate in g/day, so the unit is 1/(uM^2*day). Signed, so
    # it is carried on the linear scale (positive = synergy, negative = antagonism).
    synergy <- 8.0; label("Drug-effect interaction parameter gamma (1/(uM^2*day))") # Results, arm b1: gamma = 8.0 (CV 15%)

    # Weighted least squares with weights 1/y^2 (Methods, Data analysis), i.e. a
    # proportional error whose magnitude is not reported; held at zero.
    propSd <- fixed(0); label("Proportional residual error on tumor weight (fraction)")
  })

  model({
    # --- PK (Supplementary Table 1 micro-constants); doses in umol/kg over V1 in L/kg give uM ---
    vc_drugc1 <- exp(lvc_drugc1)
    kel_drugc1 <- exp(lkel_drugc1)
    k12_drugc1 <- exp(lk12_drugc1)
    k21_drugc1 <- exp(lk21_drugc1)
    d/dt(central_drugc1) <- -(kel_drugc1 + k12_drugc1) * central_drugc1 + k21_drugc1 * peripheral1_drugc1
    d/dt(peripheral1_drugc1) <- k12_drugc1 * central_drugc1 - k21_drugc1 * peripheral1_drugc1
    Cc_drugc1 <- central_drugc1 / vc_drugc1

    vc_cisplatin <- exp(lvc_cisplatin)
    kel_cisplatin <- exp(lkel_cisplatin)
    k12_cisplatin <- exp(lk12_cisplatin)
    k21_cisplatin <- exp(lk21_cisplatin)
    d/dt(central_cisplatin) <- -(kel_cisplatin + k12_cisplatin) * central_cisplatin + k21_cisplatin * peripheral1_cisplatin
    d/dt(peripheral1_cisplatin) <- k12_cisplatin * central_cisplatin - k21_cisplatin * peripheral1_cisplatin
    Cc_cisplatin <- central_cisplatin / vc_cisplatin

    # --- Tumor growth inhibition (Terranova 2013 Equations 1-3, Figure 1 right panel) ---
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    rbase_tumor <- exp(lrbase_tumor)
    damageTransit_drugc1 <- exp(ldamageTransit_drugc1)
    drugSlope_drugc1 <- exp(ldrugSlope_drugc1)
    damageTransit_cisplatin <- exp(ldamageTransit_cisplatin)
    drugSlope_cisplatin <- exp(ldrugSlope_cisplatin)

    killRate_drugc1 <- drugSlope_drugc1 * Cc_drugc1
    killRate_cisplatin <- drugSlope_cisplatin * Cc_cisplatin
    killRate_interaction <- synergy * Cc_drugc1 * Cc_cisplatin

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
    d/dt(cycling_cells) <- growth - (killRate_drugc1 + killRate_cisplatin + killRate_interaction) * cycling_cells
    # Equation 2 with the flows of Figure 1 (right panel): a damaged cell hit by the
    # other drug leaves its state, so every off-origin state loses k2 * C of the drug
    # it has not yet been damaged by (see the vignette for this reading).
    d/dt(damaged_cells_a0_b1) <- killRate_cisplatin * cycling_cells - (killRate_drugc1 + damageTransit_cisplatin) * damaged_cells_a0_b1
    d/dt(damaged_cells_a0_b2) <- damageTransit_cisplatin * damaged_cells_a0_b1 - (killRate_drugc1 + damageTransit_cisplatin) * damaged_cells_a0_b2
    d/dt(damaged_cells_a0_b3) <- damageTransit_cisplatin * damaged_cells_a0_b2 - (killRate_drugc1 + damageTransit_cisplatin) * damaged_cells_a0_b3
    d/dt(damaged_cells_a1_b0) <- killRate_drugc1 * cycling_cells - (damageTransit_drugc1 + killRate_cisplatin) * damaged_cells_a1_b0
    d/dt(damaged_cells_a1_b1) <- killRate_drugc1 * damaged_cells_a0_b1 + killRate_cisplatin * damaged_cells_a1_b0 + killRate_interaction * cycling_cells - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a1_b1
    d/dt(damaged_cells_a1_b2) <- killRate_drugc1 * damaged_cells_a0_b2 + damageTransit_cisplatin * damaged_cells_a1_b1 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a1_b2
    d/dt(damaged_cells_a1_b3) <- killRate_drugc1 * damaged_cells_a0_b3 + damageTransit_cisplatin * damaged_cells_a1_b2 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a1_b3
    d/dt(damaged_cells_a2_b0) <- damageTransit_drugc1 * damaged_cells_a1_b0 - (damageTransit_drugc1 + killRate_cisplatin) * damaged_cells_a2_b0
    d/dt(damaged_cells_a2_b1) <- damageTransit_drugc1 * damaged_cells_a1_b1 + killRate_cisplatin * damaged_cells_a2_b0 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a2_b1
    d/dt(damaged_cells_a2_b2) <- damageTransit_drugc1 * damaged_cells_a1_b2 + damageTransit_cisplatin * damaged_cells_a2_b1 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a2_b2
    d/dt(damaged_cells_a2_b3) <- damageTransit_drugc1 * damaged_cells_a1_b3 + damageTransit_cisplatin * damaged_cells_a2_b2 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a2_b3
    d/dt(damaged_cells_a3_b0) <- damageTransit_drugc1 * damaged_cells_a2_b0 - (damageTransit_drugc1 + killRate_cisplatin) * damaged_cells_a3_b0
    d/dt(damaged_cells_a3_b1) <- damageTransit_drugc1 * damaged_cells_a2_b1 + killRate_cisplatin * damaged_cells_a3_b0 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a3_b1
    d/dt(damaged_cells_a3_b2) <- damageTransit_drugc1 * damaged_cells_a2_b2 + damageTransit_cisplatin * damaged_cells_a3_b1 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a3_b2
    d/dt(damaged_cells_a3_b3) <- damageTransit_drugc1 * damaged_cells_a2_b3 + damageTransit_cisplatin * damaged_cells_a3_b2 - (damageTransit_drugc1 + damageTransit_cisplatin) * damaged_cells_a3_b3

    cycling_cells(0) <- rbase_tumor

    tumor_size ~ prop(propSd)
  })
}
