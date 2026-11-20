Terranova_2013_drugC1_gemcitabine_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; subcutaneous BxPC3 human pancreatic adenocarcinoma xenograft, Terranova 2013 Experiment a).",
    "Two-drug tumor-growth-inhibition model for Drug C1 (an Aurora kinase inhibitor coded by the authors; no INN disclosed) combined with gemcitabine:",
    "the mouse PK models of both drugs (two-compartment IV bolus for Drug C1, two-compartment IV bolus for gemcitabine; held at their single-agent",
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
  dosing <- c("central_drugc1", "central_gemcitabine")

  units <- list(
    time = "day",
    dosing = "umol/kg (Drug C1 into central_drugc1, gemcitabine into central_gemcitabine)",
    concentration = "uM"
  )

  covariateData <- list()

  # Terranova 2013 Equation 2 grid: damaged_cells_a<i>_b<j> is paper state x_ij, a
  # tumor cell that has progressed i stages along the Drug C1 (drug a) chain and j
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
    central_drugc1 = list(analyte = "Drug C1", units = "umol/kg", specimen = "plasma", verified = TRUE),
    peripheral1_drugc1 = list(analyte = "Drug C1", units = "umol/kg", specimen = "tissue", verified = TRUE),
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
    weight_range = "17-35 g",
    sex_female_pct = 100,
    disease_state = "subcutaneous BxPC3 human pancreatic adenocarcinoma xenograft, randomized at a palpable tumor of about 100-300 mm^3",
    dose_range = "Drug C1 15 mg/kg (29.4 umol/kg) IV twice daily on days 10-12, 14-16 and 18-20; gemcitabine 80 mg/kg (267 umol/kg) IV on days 9, 13 and 17 (Table 1, arm a1)",
    regions = "preclinical (Nerviano Medical Sciences, Italy)",
    notes = paste(
      "Terranova 2013 Experiment a: female CD1 athymic nude-nu mice, 5 weeks old, 17-35 g.",
      "Control, two single-agent and the combination arms, usually eight animals per",
      "group (exact counts not reported). Time zero is tumor implantation; tumor weight",
      "was computed from caliper length and width as length * width^2 / 2 at unit density.",
      "The seven tumor and single-agent parameters were estimated in Rocchetti 2009 (Experiment 1) from the control and single-agent arms;",
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

    # --- Gemcitabine mouse PK: Supplementary Table 1 (held at the single-agent estimates) ---
    lvc_gemcitabine <- fixed(log(0.887)); label("Gemcitabine central volume of distribution (L/kg)") # Supplementary Table 1 'Gemcytabine': V1 = 0.887 l/kg
    lkel_gemcitabine <- fixed(log(97.3)); label("Gemcitabine elimination rate constant from the central compartment (1/day)") # Supplementary Table 1 'Gemcytabine': K10 = 97.3 1/day
    lk12_gemcitabine <- fixed(log(20.4)); label("Gemcitabine central-to-peripheral rate constant (1/day)") # Supplementary Table 1 'Gemcytabine': K12 = 20.4 1/day
    lk21_gemcitabine <- fixed(log(45.2)); label("Gemcitabine peripheral-to-central rate constant (1/day)") # Supplementary Table 1 'Gemcytabine': K21 = 45.2 1/day

    # --- Tumor growth and single-agent effects: Supplementary Table 2, Experiment a ---
    ltumorExpGrowth <- log(0.149); label("Exponential-phase tumor growth rate lambda0 (1/day)") # Supplementary Table 2 Experiment a: lambda0 = 0.149 1/day
    ltumorLinGrowth <- log(0.203); label("Linear-phase tumor growth rate lambda1 (g/day)") # Supplementary Table 2 Experiment a: lambda1 = 0.203 g/day
    # Terranova 2013 keeps the Simeoni 2004 growth equation (Methods, Pharmacodynamic
    # model) without printing psi; Simeoni 2004 sets psi = 20 for a sharp switch.
    psi <- fixed(20); label("Exponential-to-linear growth switch sharpness psi (unitless)") # Simeoni 2004: psi = 20 (not printed in Terranova 2013)
    lrbase_tumor <- log(0.0566); label("Tumor weight at implantation w0 (g)") # Supplementary Table 2 Experiment a: w0 = 0.0566 g
    ldamageTransit_drugc1 <- log(2.24); label("Drug C1 damaged-cell transit rate constant k1a (1/day)") # Supplementary Table 2 Experiment a: k1A = 2.24 1/day
    ldrugSlope_drugc1 <- log(0.0512); label("Drug C1 potency k2a (1/(uM*day))") # Supplementary Table 2 Experiment a: k2A = 0.0512 1/(uM*day)
    ldamageTransit_gemcitabine <- log(1.68); label("Gemcitabine damaged-cell transit rate constant k1b (1/day)") # Supplementary Table 2 Experiment a: k1B = 1.68 1/day
    ldrugSlope_gemcitabine <- log(0.0984); label("Gemcitabine potency k2b (1/(uM*day))") # Supplementary Table 2 Experiment a: k2B = 0.0984 1/(uM*day)

    # --- Drug-effect interaction: Results ---
    # The Results print the unit as uM^-2 day^-2; Equation 1 multiplies gamma by
    # Ca * Cb * x00 to give a rate in g/day, so the unit is 1/(uM^2*day). Signed, so
    # it is carried on the linear scale (positive = synergy, negative = antagonism).
    synergy <- -0.88; label("Drug-effect interaction parameter gamma (1/(uM^2*day))") # Results, arm a1: gamma = -0.88 (CV 10%)

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
    damageTransit_drugc1 <- exp(ldamageTransit_drugc1)
    drugSlope_drugc1 <- exp(ldrugSlope_drugc1)
    damageTransit_gemcitabine <- exp(ldamageTransit_gemcitabine)
    drugSlope_gemcitabine <- exp(ldrugSlope_gemcitabine)

    killRate_drugc1 <- drugSlope_drugc1 * Cc_drugc1
    killRate_gemcitabine <- drugSlope_gemcitabine * Cc_gemcitabine
    killRate_interaction <- synergy * Cc_drugc1 * Cc_gemcitabine

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
    d/dt(cycling_cells) <- growth - (killRate_drugc1 + killRate_gemcitabine + killRate_interaction) * cycling_cells
    # Equation 2 with the flows of Figure 1 (right panel): a damaged cell hit by the
    # other drug leaves its state, so every off-origin state loses k2 * C of the drug
    # it has not yet been damaged by (see the vignette for this reading).
    d/dt(damaged_cells_a0_b1) <- killRate_gemcitabine * cycling_cells - (killRate_drugc1 + damageTransit_gemcitabine) * damaged_cells_a0_b1
    d/dt(damaged_cells_a0_b2) <- damageTransit_gemcitabine * damaged_cells_a0_b1 - (killRate_drugc1 + damageTransit_gemcitabine) * damaged_cells_a0_b2
    d/dt(damaged_cells_a0_b3) <- damageTransit_gemcitabine * damaged_cells_a0_b2 - (killRate_drugc1 + damageTransit_gemcitabine) * damaged_cells_a0_b3
    d/dt(damaged_cells_a1_b0) <- killRate_drugc1 * cycling_cells - (damageTransit_drugc1 + killRate_gemcitabine) * damaged_cells_a1_b0
    d/dt(damaged_cells_a1_b1) <- killRate_drugc1 * damaged_cells_a0_b1 + killRate_gemcitabine * damaged_cells_a1_b0 + killRate_interaction * cycling_cells - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a1_b1
    d/dt(damaged_cells_a1_b2) <- killRate_drugc1 * damaged_cells_a0_b2 + damageTransit_gemcitabine * damaged_cells_a1_b1 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a1_b2
    d/dt(damaged_cells_a1_b3) <- killRate_drugc1 * damaged_cells_a0_b3 + damageTransit_gemcitabine * damaged_cells_a1_b2 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a1_b3
    d/dt(damaged_cells_a2_b0) <- damageTransit_drugc1 * damaged_cells_a1_b0 - (damageTransit_drugc1 + killRate_gemcitabine) * damaged_cells_a2_b0
    d/dt(damaged_cells_a2_b1) <- damageTransit_drugc1 * damaged_cells_a1_b1 + killRate_gemcitabine * damaged_cells_a2_b0 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a2_b1
    d/dt(damaged_cells_a2_b2) <- damageTransit_drugc1 * damaged_cells_a1_b2 + damageTransit_gemcitabine * damaged_cells_a2_b1 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a2_b2
    d/dt(damaged_cells_a2_b3) <- damageTransit_drugc1 * damaged_cells_a1_b3 + damageTransit_gemcitabine * damaged_cells_a2_b2 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a2_b3
    d/dt(damaged_cells_a3_b0) <- damageTransit_drugc1 * damaged_cells_a2_b0 - (damageTransit_drugc1 + killRate_gemcitabine) * damaged_cells_a3_b0
    d/dt(damaged_cells_a3_b1) <- damageTransit_drugc1 * damaged_cells_a2_b1 + killRate_gemcitabine * damaged_cells_a3_b0 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a3_b1
    d/dt(damaged_cells_a3_b2) <- damageTransit_drugc1 * damaged_cells_a2_b2 + damageTransit_gemcitabine * damaged_cells_a3_b1 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a3_b2
    d/dt(damaged_cells_a3_b3) <- damageTransit_drugc1 * damaged_cells_a2_b3 + damageTransit_gemcitabine * damaged_cells_a3_b2 - (damageTransit_drugc1 + damageTransit_gemcitabine) * damaged_cells_a3_b3

    cycling_cells(0) <- rbase_tumor

    tumor_size ~ prop(propSd)
  })
}
