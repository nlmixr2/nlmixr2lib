Terranova_2013_drugC4_gemcitabine_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; subcutaneous BxPC3 human pancreatic adenocarcinoma xenograft, Terranova 2013 Experiment e).",
    "Two-drug tumor-growth-inhibition model for Drug C4 (a multikinase inhibitor coded by the authors; no INN disclosed) combined with gemcitabine:",
    "the mouse PK models of both drugs (one-compartment first-order oral absorption for Drug C4, two-compartment IV bolus for gemcitabine; held at their single-agent",
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
  dosing <- c("depot_drugc4", "central_gemcitabine")

  units <- list(
    time = "day",
    dosing = "umol/kg (Drug C4 into depot_drugc4, gemcitabine into central_gemcitabine)",
    concentration = "uM"
  )

  covariateData <- list()

  # Terranova 2013 Equation 2 grid: damaged_cells_a<i>_b<j> is paper state x_ij, a
  # tumor cell that has progressed i stages along the Drug C4 (drug a) chain and j
  # stages along the gemcitabine (drug b) chain; x_00 is cycling_cells.
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
    depot_drugc4 = list(analyte = "Drug C4", units = "umol/kg", specimen = "administration site", verified = TRUE),
    central_drugc4 = list(analyte = "Drug C4", units = "umol/kg", specimen = "plasma", verified = TRUE),
    central_gemcitabine = list(analyte = "Gemcitabine", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_gemcitabine = list(analyte = "Gemcitabine", units = "umol/kg", specimen = "tissue", verified = TRUE),
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
    species = "mouse (athymic nude-nu bearing BxPC3 xenografts)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "5 weeks at purchase",
    weight_range = "25-35 g",
    sex_female_pct = 0,
    disease_state = "subcutaneous BxPC3 human pancreatic adenocarcinoma xenograft, randomized at a palpable tumor of about 100-300 mm^3",
    dose_range = "Drug C4 20 mg/kg (37.5 umol/kg, arm e1) or 40 mg/kg (75 umol/kg, arm e2, predicted) orally twice daily on days 10-12, 14-16 and 18-20; gemcitabine 80 mg/kg (267 umol/kg) IV on days 9, 13 and 17 (Table 1)",
    regions = "preclinical (Nerviano Medical Sciences, Italy)",
    notes = paste(
      "Terranova 2013 Experiment e: male Balb athymic nude-nu mice, 5 weeks old, 25-35 g.",
      "Control, two single-agent and the combination arms, usually eight animals per",
      "group (exact counts not reported). Time zero is tumor implantation; tumor weight",
      "was computed from caliper length and width as length * width^2 / 2 at unit density.",
      "The seven tumor and single-agent parameters were estimated in this paper from the control, single-agent Drug C4 20 mg/kg and single-agent gemcitabine arms (RMSE 0.06 g);",
      "the interaction parameter was then estimated on the combination arm with them held",
      "at those values, by weighted least squares on the arm means."
    )
  )

  ini({
    # --- Drug C4 mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    # Table 1 reports no bioavailability for the oral drug, so V1 is read as the
    # apparent volume V/F (bioavailability 1 in the simulation; see the vignette).
    lvc_drugc4 <- fixed(log(13.3)); label("Drug C4 central volume of distribution (L/kg)") # Supplementary Table 1 'Drug C4': V1 = 13.3 l/kg
    lkel_drugc4 <- fixed(log(3.36)); label("Drug C4 elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Drug C4': K10 = 3.36 1/day
    lka_drugc4 <- fixed(log(206)); label("Drug C4 first-order oral absorption rate constant (1/day)") # Supplementary Table 1 'Drug C4': K01 = 206 1/day

    # --- Gemcitabine mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_gemcitabine <- fixed(log(0.887)); label("Gemcitabine central volume of distribution (L/kg)") # Supplementary Table 1 'Gemcytabine': V1 = 0.887 l/kg
    lkel_gemcitabine <- fixed(log(97.3)); label("Gemcitabine elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Gemcytabine': K10 = 97.3 1/day
    lk12_gemcitabine <- fixed(log(20.4)); label("Gemcitabine central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Gemcytabine': K12 = 20.4 1/day
    lk21_gemcitabine <- fixed(log(45.2)); label("Gemcitabine peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Gemcytabine': K21 = 45.2 1/day

    # --- Tumor growth and single-agent effects: Supplementary Table 2, Experiment e ---
    ltumorExpGrowth <- log(0.153); label("Exponential-phase tumor growth rate lambda0 (1/day)") # Supplementary Table 2 Experiment e: lambda0 = 0.153 1/day (4%)
    ltumorLinGrowth <- log(0.223); label("Linear-phase tumor growth rate lambda1 (g/day)") # Supplementary Table 2 Experiment e: lambda1 = 0.223 g/day (11%)
    # Terranova 2013 keeps the Simeoni 2004 growth equation (Methods, Pharmacodynamic
    # model) without printing psi; Simeoni 2004 sets psi = 20 for a sharp switch.
    psi <- fixed(20); label("Exponential-to-linear growth switch sharpness psi (unitless)") # Simeoni 2004: psi = 20 (not printed in Terranova 2013)
    lrbase_tumor <- log(0.0533); label("Tumor weight at implantation w0 (g)") # Supplementary Table 2 Experiment e: w0 = 0.0533 g (9%)
    ldamageTransit_drugc4 <- log(62.6); label("Drug C4 damaged-cell transit rate constant k1a (1/day)") # Supplementary Table 2 Experiment e: k1A = 62.6 1/day (>100%)
    ldrugSlope_drugc4 <- log(0.0526); label("Drug C4 potency k2a (1/(uM*day))") # Supplementary Table 2 Experiment e: k2A = 0.0526 1/(uM*day) (6%)
    ldamageTransit_gemcitabine <- log(1.49); label("Gemcitabine damaged-cell transit rate constant k1b (1/day)") # Supplementary Table 2 Experiment e: k1B = 1.49 1/day (58%)
    ldrugSlope_gemcitabine <- log(0.104); label("Gemcitabine potency k2b (1/(uM*day))") # Supplementary Table 2 Experiment e: k2B = 0.104 1/(uM*day) (10%)

    # --- Drug-effect interaction: Results ---
    # The Results print the unit as uM^-2 day^-2; Equation 1 multiplies gamma by
    # Ca * Cb * x00 to give a rate in g/day, so the unit is 1/(uM^2*day). Signed, so
    # it is carried on the linear scale (positive = synergy, negative = antagonism).
    synergy <- -0.0455; label("Drug-effect interaction parameter gamma (1/(uM^2*day))") # Results, arm e1: gamma = -0.0455 (CV 28%)

    # Weighted least squares with weights 1/y^2 (Methods, Data analysis), i.e. a
    # proportional error whose magnitude is not reported; held at zero.
    propSd <- fixed(0); label("Proportional residual error on tumor weight (fraction)")
  })

  model({
    # --- PK (Supplementary Table 1 micro-constants); doses in umol/kg over V1 in L/kg give uM ---
    vc_drugc4 <- exp(lvc_drugc4)
    kel_drugc4 <- exp(lkel_drugc4)
    ka_drugc4 <- exp(lka_drugc4)
    d/dt(depot_drugc4) <- -ka_drugc4 * depot_drugc4
    d/dt(central_drugc4) <- ka_drugc4 * depot_drugc4 - kel_drugc4 * central_drugc4
    Cc_drugc4 <- central_drugc4 / vc_drugc4

    vc_gemcitabine <- exp(lvc_gemcitabine)
    kel_gemcitabine <- exp(lkel_gemcitabine)
    k12_gemcitabine <- exp(lk12_gemcitabine)
    k21_gemcitabine <- exp(lk21_gemcitabine)
    d/dt(central_gemcitabine) <- -(kel_gemcitabine + k12_gemcitabine) * central_gemcitabine + k21_gemcitabine * peripheral1_gemcitabine
    d/dt(peripheral1_gemcitabine) <- k12_gemcitabine * central_gemcitabine - k21_gemcitabine * peripheral1_gemcitabine
    Cc_gemcitabine <- central_gemcitabine / vc_gemcitabine

    # --- Tumor growth inhibition (Terranova 2013 Equations 1-3, Figure 1 right panel) ---
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    rbase_tumor <- exp(lrbase_tumor)
    damageTransit_drugc4 <- exp(ldamageTransit_drugc4)
    drugSlope_drugc4 <- exp(ldrugSlope_drugc4)
    damageTransit_gemcitabine <- exp(ldamageTransit_gemcitabine)
    drugSlope_gemcitabine <- exp(ldrugSlope_gemcitabine)

    killRate_drugc4 <- drugSlope_drugc4 * Cc_drugc4
    killRate_gemcitabine <- drugSlope_gemcitabine * Cc_gemcitabine
    killRate_interaction <- synergy * Cc_drugc4 * Cc_gemcitabine

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
    d/dt(cycling_cells) <- growth - (killRate_drugc4 + killRate_gemcitabine + killRate_interaction) * cycling_cells
    # Equation 2 with the flows of Figure 1 (right panel): a damaged cell hit by the
    # other drug leaves its state, so every off-origin state loses k2 * C of the drug
    # it has not yet been damaged by (see the vignette for this reading).
    d/dt(damaged_cells_a0_b1) <- killRate_gemcitabine * cycling_cells - (killRate_drugc4 + damageTransit_gemcitabine) * damaged_cells_a0_b1
    d/dt(damaged_cells_a0_b2) <- damageTransit_gemcitabine * damaged_cells_a0_b1 - (killRate_drugc4 + damageTransit_gemcitabine) * damaged_cells_a0_b2
    d/dt(damaged_cells_a0_b3) <- damageTransit_gemcitabine * damaged_cells_a0_b2 - (killRate_drugc4 + damageTransit_gemcitabine) * damaged_cells_a0_b3
    d/dt(damaged_cells_a1_b0) <- killRate_drugc4 * cycling_cells - (damageTransit_drugc4 + killRate_gemcitabine) * damaged_cells_a1_b0
    d/dt(damaged_cells_a1_b1) <- killRate_drugc4 * damaged_cells_a0_b1 + killRate_gemcitabine * damaged_cells_a1_b0 + killRate_interaction * cycling_cells - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a1_b1
    d/dt(damaged_cells_a1_b2) <- killRate_drugc4 * damaged_cells_a0_b2 + damageTransit_gemcitabine * damaged_cells_a1_b1 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a1_b2
    d/dt(damaged_cells_a1_b3) <- killRate_drugc4 * damaged_cells_a0_b3 + damageTransit_gemcitabine * damaged_cells_a1_b2 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a1_b3
    d/dt(damaged_cells_a2_b0) <- damageTransit_drugc4 * damaged_cells_a1_b0 - (damageTransit_drugc4 + killRate_gemcitabine) * damaged_cells_a2_b0
    d/dt(damaged_cells_a2_b1) <- damageTransit_drugc4 * damaged_cells_a1_b1 + killRate_gemcitabine * damaged_cells_a2_b0 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a2_b1
    d/dt(damaged_cells_a2_b2) <- damageTransit_drugc4 * damaged_cells_a1_b2 + damageTransit_gemcitabine * damaged_cells_a2_b1 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a2_b2
    d/dt(damaged_cells_a2_b3) <- damageTransit_drugc4 * damaged_cells_a1_b3 + damageTransit_gemcitabine * damaged_cells_a2_b2 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a2_b3
    d/dt(damaged_cells_a3_b0) <- damageTransit_drugc4 * damaged_cells_a2_b0 - (damageTransit_drugc4 + killRate_gemcitabine) * damaged_cells_a3_b0
    d/dt(damaged_cells_a3_b1) <- damageTransit_drugc4 * damaged_cells_a2_b1 + killRate_gemcitabine * damaged_cells_a3_b0 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a3_b1
    d/dt(damaged_cells_a3_b2) <- damageTransit_drugc4 * damaged_cells_a2_b2 + damageTransit_gemcitabine * damaged_cells_a3_b1 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a3_b2
    d/dt(damaged_cells_a3_b3) <- damageTransit_drugc4 * damaged_cells_a2_b3 + damageTransit_gemcitabine * damaged_cells_a3_b2 - (damageTransit_drugc4 + damageTransit_gemcitabine) * damaged_cells_a3_b3

    cycling_cells(0) <- rbase_tumor

    tumor_size ~ prop(propSd)
  })
}
