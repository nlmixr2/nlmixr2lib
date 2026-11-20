Terranova_2013_drugC5_irinotecan_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; subcutaneous HT29 human colon adenocarcinoma xenograft, Terranova 2013 Experiment f).",
    "Two-drug tumor-growth-inhibition model for Drug C5 (a multikinase inhibitor coded by the authors; no INN disclosed) combined with irinotecan (CPT-11):",
    "the mouse PK models of both drugs (two-compartment IV bolus for Drug C5, two-compartment IV bolus for irinotecan; held at their single-agent",
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
  dosing <- c("central_drugc5", "central_irinotecan")

  units <- list(
    time = "day",
    dosing = "umol/kg (Drug C5 into central_drugc5, irinotecan into central_irinotecan)",
    concentration = "uM"
  )

  covariateData <- list()

  # Terranova 2013 Equation 2 grid: damaged_cells_a<i>_b<j> is paper state x_ij, a
  # tumor cell that has progressed i stages along the Drug C5 (drug a) chain and j
  # stages along the irinotecan (drug b) chain; x_00 is cycling_cells.
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
    central_drugc5 = list(analyte = "Drug C5", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_drugc5 = list(analyte = "Drug C5", units = "umol/kg", specimen = "tissue", verified = TRUE),
    central_irinotecan = list(analyte = "Irinotecan", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_irinotecan = list(analyte = "Irinotecan", units = "umol/kg", specimen = "tissue", verified = TRUE),
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
    dose_range = "Drug C5 10 mg/kg (25.1 umol/kg, arm f1) or 20 mg/kg (50.3 umol/kg, arm f2) IV once daily on days 9-16, given 6 h after irinotecan on shared days; irinotecan 45 mg/kg (72.9 umol/kg) IV on days 8, 12 and 16 (Table 1)",
    regions = "preclinical (Nerviano Medical Sciences, Italy)",
    notes = paste(
      "Terranova 2013 Experiment f: male Balb athymic nude-nu mice, 5 weeks old, 25-35 g.",
      "Control, two single-agent and the combination arms, usually eight animals per",
      "group (exact counts not reported). Time zero is tumor implantation; tumor weight",
      "was computed from caliper length and width as length * width^2 / 2 at unit density.",
      "The seven tumor and single-agent parameters were estimated in this paper from the control, two single-agent Drug C5 (10 and 20 mg/kg) and single-agent irinotecan arms (RMSE 0.19 g);",
      "the interaction parameter was then estimated on the combination arm with them held",
      "at those values, by weighted least squares on the arm means."
    )
  )

  ini({
    # --- Drug C5 mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_drugc5 <- fixed(log(0.577)); label("Drug C5 central volume of distribution (L/kg)") # Supplementary Table 1 'Drug C5': V1 = 0.577 l/kg
    lkel_drugc5 <- fixed(log(89)); label("Drug C5 elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Drug C5': K10 = 89 1/day
    lk12_drugc5 <- fixed(log(13.9)); label("Drug C5 central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Drug C5': K12 = 13.9 1/day
    lk21_drugc5 <- fixed(log(13.7)); label("Drug C5 peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Drug C5': K21 = 13.7 1/day

    # --- Irinotecan mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_irinotecan <- fixed(log(4.85)); label("Irinotecan central volume of distribution (L/kg)") # Supplementary Table 1 'Irinotecan': V1 = 4.85 l/kg
    lkel_irinotecan <- fixed(log(13.3)); label("Irinotecan elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Irinotecan': K10 = 13.3 1/day
    lk12_irinotecan <- fixed(log(0.275)); label("Irinotecan central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Irinotecan': K12 = 0.275 1/day
    lk21_irinotecan <- fixed(log(1.48)); label("Irinotecan peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Irinotecan': K21 = 1.48 1/day

    # --- Tumor growth and single-agent effects: Supplementary Table 2, Experiment f ---
    ltumorExpGrowth <- log(0.231); label("Exponential-phase tumor growth rate lambda0 (1/day)") # Supplementary Table 2 Experiment f: lambda0 = 0.231 1/day (7%)
    ltumorLinGrowth <- log(0.252); label("Linear-phase tumor growth rate lambda1 (g/day)") # Supplementary Table 2 Experiment f: lambda1 = 0.252 g/day (13%)
    # Terranova 2013 keeps the Simeoni 2004 growth equation (Methods, Pharmacodynamic
    # model) without printing psi; Simeoni 2004 sets psi = 20 for a sharp switch.
    psi <- fixed(20); label("Exponential-to-linear growth switch sharpness psi (unitless)") # Simeoni 2004: psi = 20 (not printed in Terranova 2013)
    lrbase_tumor <- log(0.0361); label("Tumor weight at implantation w0 (g)") # Supplementary Table 2 Experiment f: w0 = 0.0361 g (15%)
    ldamageTransit_drugc5 <- log(1.30); label("Drug C5 damaged-cell transit rate constant k1a (1/day)") # Supplementary Table 2 Experiment f: k1A = 1.30 1/day (31%)
    ldrugSlope_drugc5 <- log(0.292); label("Drug C5 potency k2a (1/(uM*day))") # Supplementary Table 2 Experiment f: k2A = 0.292 1/(uM*day) (10%)
    ldamageTransit_irinotecan <- log(0.109); label("Irinotecan damaged-cell transit rate constant k1b (1/day)") # Supplementary Table 2 Experiment f: k1B = 0.109 1/day (36%)
    ldrugSlope_irinotecan <- log(1.65); label("Irinotecan potency k2b (1/(uM*day))") # Supplementary Table 2 Experiment f: k2B = 1.65 1/(uM*day) (9%)

    # --- Drug-effect interaction: Results ---
    # The Results print the unit as uM^-2 day^-2; Equation 1 multiplies gamma by
    # Ca * Cb * x00 to give a rate in g/day, so the unit is 1/(uM^2*day). Signed, so
    # it is carried on the linear scale (positive = synergy, negative = antagonism).
    synergy <- -2.04; label("Drug-effect interaction parameter gamma (1/(uM^2*day))") # Results, arm f1: gamma = -2.04 (CV 18%)
    # Arm f2 was fitted independently: gamma = -1.05 1/(uM^2*day) (CV 30%). Simulate it with
    # ini(synergy = -1.05).

    # Weighted least squares with weights 1/y^2 (Methods, Data analysis), i.e. a
    # proportional error whose magnitude is not reported; held at zero.
    propSd <- fixed(0); label("Proportional residual error on tumor weight (fraction)")
  })

  model({
    # --- PK (Supplementary Table 1 micro-constants); doses in umol/kg over V1 in L/kg give uM ---
    vc_drugc5 <- exp(lvc_drugc5)
    kel_drugc5 <- exp(lkel_drugc5)
    k12_drugc5 <- exp(lk12_drugc5)
    k21_drugc5 <- exp(lk21_drugc5)
    d/dt(central_drugc5) <- -(kel_drugc5 + k12_drugc5) * central_drugc5 + k21_drugc5 * peripheral1_drugc5
    d/dt(peripheral1_drugc5) <- k12_drugc5 * central_drugc5 - k21_drugc5 * peripheral1_drugc5
    Cc_drugc5 <- central_drugc5 / vc_drugc5

    vc_irinotecan <- exp(lvc_irinotecan)
    kel_irinotecan <- exp(lkel_irinotecan)
    k12_irinotecan <- exp(lk12_irinotecan)
    k21_irinotecan <- exp(lk21_irinotecan)
    d/dt(central_irinotecan) <- -(kel_irinotecan + k12_irinotecan) * central_irinotecan + k21_irinotecan * peripheral1_irinotecan
    d/dt(peripheral1_irinotecan) <- k12_irinotecan * central_irinotecan - k21_irinotecan * peripheral1_irinotecan
    Cc_irinotecan <- central_irinotecan / vc_irinotecan

    # --- Tumor growth inhibition (Terranova 2013 Equations 1-3, Figure 1 right panel) ---
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    rbase_tumor <- exp(lrbase_tumor)
    damageTransit_drugc5 <- exp(ldamageTransit_drugc5)
    drugSlope_drugc5 <- exp(ldrugSlope_drugc5)
    damageTransit_irinotecan <- exp(ldamageTransit_irinotecan)
    drugSlope_irinotecan <- exp(ldrugSlope_irinotecan)

    killRate_drugc5 <- drugSlope_drugc5 * Cc_drugc5
    killRate_irinotecan <- drugSlope_irinotecan * Cc_irinotecan
    killRate_interaction <- synergy * Cc_drugc5 * Cc_irinotecan

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
    d/dt(cycling_cells) <- growth - (killRate_drugc5 + killRate_irinotecan + killRate_interaction) * cycling_cells
    # Equation 2 with the flows of Figure 1 (right panel): a damaged cell hit by the
    # other drug leaves its state, so every off-origin state loses k2 * C of the drug
    # it has not yet been damaged by (see the vignette for this reading).
    d/dt(damaged_cells_a0_b1) <- killRate_irinotecan * cycling_cells - (killRate_drugc5 + damageTransit_irinotecan) * damaged_cells_a0_b1
    d/dt(damaged_cells_a0_b2) <- damageTransit_irinotecan * damaged_cells_a0_b1 - (killRate_drugc5 + damageTransit_irinotecan) * damaged_cells_a0_b2
    d/dt(damaged_cells_a0_b3) <- damageTransit_irinotecan * damaged_cells_a0_b2 - (killRate_drugc5 + damageTransit_irinotecan) * damaged_cells_a0_b3
    d/dt(damaged_cells_a1_b0) <- killRate_drugc5 * cycling_cells - (damageTransit_drugc5 + killRate_irinotecan) * damaged_cells_a1_b0
    d/dt(damaged_cells_a1_b1) <- killRate_drugc5 * damaged_cells_a0_b1 + killRate_irinotecan * damaged_cells_a1_b0 + killRate_interaction * cycling_cells - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a1_b1
    d/dt(damaged_cells_a1_b2) <- killRate_drugc5 * damaged_cells_a0_b2 + damageTransit_irinotecan * damaged_cells_a1_b1 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a1_b2
    d/dt(damaged_cells_a1_b3) <- killRate_drugc5 * damaged_cells_a0_b3 + damageTransit_irinotecan * damaged_cells_a1_b2 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a1_b3
    d/dt(damaged_cells_a2_b0) <- damageTransit_drugc5 * damaged_cells_a1_b0 - (damageTransit_drugc5 + killRate_irinotecan) * damaged_cells_a2_b0
    d/dt(damaged_cells_a2_b1) <- damageTransit_drugc5 * damaged_cells_a1_b1 + killRate_irinotecan * damaged_cells_a2_b0 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a2_b1
    d/dt(damaged_cells_a2_b2) <- damageTransit_drugc5 * damaged_cells_a1_b2 + damageTransit_irinotecan * damaged_cells_a2_b1 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a2_b2
    d/dt(damaged_cells_a2_b3) <- damageTransit_drugc5 * damaged_cells_a1_b3 + damageTransit_irinotecan * damaged_cells_a2_b2 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a2_b3
    d/dt(damaged_cells_a3_b0) <- damageTransit_drugc5 * damaged_cells_a2_b0 - (damageTransit_drugc5 + killRate_irinotecan) * damaged_cells_a3_b0
    d/dt(damaged_cells_a3_b1) <- damageTransit_drugc5 * damaged_cells_a2_b1 + killRate_irinotecan * damaged_cells_a3_b0 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a3_b1
    d/dt(damaged_cells_a3_b2) <- damageTransit_drugc5 * damaged_cells_a2_b2 + damageTransit_irinotecan * damaged_cells_a3_b1 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a3_b2
    d/dt(damaged_cells_a3_b3) <- damageTransit_drugc5 * damaged_cells_a2_b3 + damageTransit_irinotecan * damaged_cells_a3_b2 - (damageTransit_drugc5 + damageTransit_irinotecan) * damaged_cells_a3_b3

    cycling_cells(0) <- rbase_tumor

    tumor_size ~ prop(propSd)
  })
}
