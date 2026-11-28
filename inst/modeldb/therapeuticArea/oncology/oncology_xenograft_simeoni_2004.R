oncology_xenograft_simeoni_2004 <- function() {
  description <- paste(
    "Oncology tumor growth model in xenograft models (preclinical, mouse):",
    "Simeoni 2004 tumor growth inhibition (TGI) model for subcutaneous human",
    "tumor xenografts in nude mice. Unperturbed growth switches smoothly from",
    "exponential (lambda0) to linear (lambda1); the drug moves proliferating",
    "cells into a three-compartment chain of damaged, non-proliferating cells",
    "at a rate proportional to the plasma concentration Cc (potency k2), and",
    "damaged cells transit to death at rate k1. Parameter values are the",
    "paclitaxel experiment 1 estimates (A2780 ovarian carcinoma). The drug",
    "concentration Cc and the initial tumor weight tumor_vol0 are inputs."
  )
  reference <- paste(
    "Simeoni M, Magni P, Cammia C, De Nicolao G, Croci V, Pesenti E, Germani M,",
    "Poggesi I, Rocchetti M. Predictive pharmacokinetic-pharmacodynamic modeling",
    "of tumor growth kinetics in xenograft models after administration of",
    "anticancer agents. Cancer Res. 2004;64(3):1094-1101.",
    "doi:10.1158/0008-5472.CAN-03-2524"
  )
  vignette <- "oncology_xenograft_simeoni_2004"
  # The drug enters as its plasma concentration Cc, supplied in the data or by
  # a PK model appended upstream; the model itself has no dosing compartment.
  dosing <- character(0)
  units <- list(
    time = "day",
    dosing = "not applicable (drug input is the plasma concentration Cc in ng/mL)",
    concentration = "ng/mL (drug input Cc); tumor weight tumor_vol in g"
  )
  depends <- c("Cc", "tumor_vol0")

  covariateData <- list(
    Cc = list(
      description = "Plasma concentration of the anticancer drug driving tumor-cell damage",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying model input, normally produced by a PK model composed",
        "upstream (Simeoni 2004 used two-compartment PK fitted separately in",
        "WinNonlin, e.g. paclitaxel V1 = 0.81 L/kg, k10 = 0.868 1/h, k12 = 0.0060",
        "1/h, k21 = 0.0838 1/h, Results p. 1098). Must be in ng/mL to match the",
        "units of k2 (Table 2, mL/(ng*day)); 0 for vehicle controls."
      ),
      source_name = "c(t)"
    ),
    tumor_vol0 = list(
      description = "Tumor weight at inoculation (time 0), the initial value of the proliferating-cell state",
      units = "g",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Simeoni 2004 estimates w0 per experiment as a model parameter (Table 2",
        "p. 1098: 0.055 g for paclitaxel experiment 1, 0.033 g for experiment",
        "2); this model takes it from the data instead, so supply the published",
        "w0 to reproduce the paper or a measured/estimated value otherwise.",
        "Tumor weight = length * width^2 / 2 at density 1 mg/mm^3 (Methods",
        "p. 1095), so 1 g corresponds to 1000 mm^3 of caliper volume."
      ),
      source_name = "w0"
    )
  )

  compartmentData <- list(
    cycling_cells = list(analyte = "proliferating tumor cells (x1)", units = "g", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(
      analyte = "drug-damaged non-proliferating tumor cells, stage 1 (x2)",
      units = "g",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells2 = list(
      analyte = "drug-damaged non-proliferating tumor cells, stage 2 (x3)",
      units = "g",
      specimen = "tumor",
      verified = TRUE
    ),
    damaged_cells3 = list(
      analyte = "drug-damaged non-proliferating tumor cells, stage 3 (x4)",
      units = "g",
      specimen = "tumor",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female athymic Hsd nude nu/nu; subcutaneous A2780 human ovarian carcinoma xenograft)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    sex_female_pct = 100,
    age_range = "5-6 weeks",
    weight_range = "20-22 g",
    disease_state = "subcutaneous A2780 human ovarian carcinoma xenograft (tumor fragments implanted in the left flank)",
    dose_range = "paclitaxel 30 mg/kg IV bolus every 4 days for 3 doses (q4d x 3), 8 mice per treated group",
    notes = paste(
      "Simeoni 2004 Methods p. 1094-1095. Parameters were estimated by weighted",
      "nonlinear least squares (WinNonlin 3.1) on group-mean tumor weights of",
      "control and treated arms of two paclitaxel experiments fitted",
      "simultaneously, sharing k1 and k2 and estimating lambda0, lambda1 and w0",
      "per experiment; control-group sizes are not stated. The same structure",
      "was fitted to CPT-11, 5-FU and two Pharmacia candidates (Tables 1, 3-5)."
    )
  )

  # The values for damageTransit (k1), drugSlope (k2), tumorExpGrowth (lambda0)
  # and tumorLinGrowth (lambda1) are the paclitaxel experiment 1 estimates of
  # Table 2 (p. 1098), estimated on the log scale to keep them positive. The
  # lower and upper limits are not from the reference. Table 2 experiment 1 is
  # the parameter set whose curves in Figure 7 start treatment on Day 13 at a
  # tumor weight of about 2 g, although the Methods and the Table 2 and
  # Figure 7 captions place experiment 1's first dose on Day 8; see the
  # vignette. Residual errors are not in the reference.
  ini({
    ldamageTransit <- log(c(0.1, 0.968, 10)) ; label("Transit rate through damage k1 (1/day)") # Table 2 p. 1098: k1 = 0.968 1/day (CV 14.53)
    ldrugSlope <- log(c(0.00001, 0.000629, 0.1)) ; label("Linear drug potency on proliferating cells k2 (mL/(ng*day))") # Table 2 p. 1098: k2 = 6.29e-4 mL/(ng*day) (CV 8.03)
    ltumorExpGrowth <- log(c(0.001, 0.273, 2)) ; label("Exponential tumor growth rate lambda0 (1/day)") # Table 2 p. 1098: lambda0 exp 1 = 0.273 1/day (CV 8.08)
    ltumorLinGrowth <- log(c(0.01, 0.814, 5)) ; label("Linear tumor growth rate lambda1 (g/day)") # Table 2 p. 1098: lambda1 exp 1 = 0.814 g/day (CV 2.34)
    propSd_tumor_vol <- 0.2 ; label("Proportional residual error on tumor weight (fraction)") # not in the reference
    addSd_tumor_vol <- 30 ; label("Additive residual error on tumor weight (g)") # not in the reference
  })
  model({
    damageTransit <- exp(ldamageTransit)
    drugSlope <- exp(ldrugSlope)
    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)

    # tumor_vol0 (w0, g) is provided in the data as the initial tumor weight.
    # It can also be estimated.
    cycling_cells(0) <- tumor_vol0
    # psi = 20 makes the switch from exponential to linear growth sharp
    # (Eq. D and text, p. 1095).
    psi <- 20
    # Total tumor weight w(t) in g (Eq. E, p. 1096).
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    # Cc (ng/mL) is provided in the data or by an appended PK model as the
    # drug concentration; its units set those of k2.
    drugEffectCyclingCells <- drugSlope*Cc
    # Eq. E (p. 1096): x1 = cycling_cells; x2, x3, x4 = damaged_cells1-3.
    d/dt(cycling_cells) <- tumorExpGrowth*cycling_cells/(1 + (tumorExpGrowth/tumorLinGrowth * tumor_vol)^psi)^(1/psi) - drugEffectCyclingCells*cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells*cycling_cells - damageTransit*damaged_cells1
    d/dt(damaged_cells2) <- damageTransit*(damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit*(damaged_cells2 - damaged_cells3)
    tumor_vol ~ prop(propSd_tumor_vol) + add(addSd_tumor_vol)
  })
}

attr(
  oncology_xenograft_simeoni_2004,
  "message"
) <- "You can modify the number of damaged cell compartments in the model using the function updateOncologyXenograftSimeoni2004(model, ncmt)"
oncology_xenograft_simeoni_2004
