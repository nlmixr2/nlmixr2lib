Choi_2021_doxorubicin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; PK in male athymic nude mice, tumor growth in female",
    "SCID mice bearing orthotopic human 143B osteosarcoma xenografts).",
    "Doxorubicin monotherapy PK-PD model (Choi 2021 model A): a",
    "two-compartment linear IV PK model parameterised with micro-constants",
    "drives a Simeoni-type tumor-growth-inhibition model with Koch 2009",
    "exponential-then-linear natural growth and a three-stage damaged-cell",
    "transit chain. Typical-value (naive-pooled SimBiology) fit with an",
    "exponential residual error and no between-animal variability."
  )
  reference <- paste(
    "Choi YH, Zhang C, Liu Z, Tu MJ, Yu AX, Yu AM.",
    "A Novel Integrated Pharmacokinetic-Pharmacodynamic Model to Evaluate",
    "Combination Therapy and Determine In Vivo Synergism.",
    "J Pharmacol Exp Ther. 2021;377(3):305-315. doi:10.1124/jpet.121.000584.",
    "Tumor-growth data and dosing regimen from Jian C, Tu MJ, Ho PY, et al.",
    "Oncotarget. 2017;8(19):30742-30755. doi:10.18632/oncotarget.16372.",
    sep = " "
  )
  vignette <- "Choi_2021_doxorubicin_sorafenib_mouse"

  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "doxorubicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "doxorubicin", units = "mg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (PK: male athymic nude Foxn1nu; tumor growth: female CB17 SCID with orthotopic 143B osteosarcoma xenograft)",
    n_subjects = 20L,
    n_studies = 2L,
    age_range = "7 weeks at purchase (PK mice)",
    weight_range = "approximately 30 g (PK mice)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "PK: tumor-free mice. PD: orthotopic (intratibial) human 143B osteosarcoma cell-line-derived xenograft.",
    dose_range = "PK: single 0.06 mg IV bolus per mouse. PD: 12 ug/mouse IV every other day x4 (cycle 1) then 10 ug/mouse IV every other day x4 (cycle 2), per Jian 2017.",
    regions = "University of California Davis, USA (preclinical)",
    notes = paste(
      "PK: n = 6 male athymic nude mice, serial microsampling to 24 h",
      "(Choi 2021 Methods, PK Studies; Table 1). PD: n = 7 doxorubicin-treated",
      "tumor-bearing mice from Jian 2017 (Choi 2021 Methods; Table 2); natural",
      "growth rates L0 and L1 were estimated from n = 7 vehicle-control mice and",
      "held fixed. The PK and PD were fitted sequentially in SimBiology",
      "(MATLAB 2019a); only the PD parameters w0, k1 and k2 were estimated",
      "against the doxorubicin tumor data."
    )
  )

  ini({
    # --- Doxorubicin PK: Choi 2021 Table 1, 'Dox (iv)' column --------------
    # Rate constants are printed in 1/h; model time is days, so model()
    # multiplies each rate constant by 24.
    lvc  <- log(0.0764); label("Central volume of distribution Vc (L)")                                   # Table 1 Dox (iv): Vc = 0.0764 L (23.4 %CV)
    lvp  <- log(3.52);   label("Peripheral volume of distribution Vp (L)")                                # Table 1 Dox (iv): Vp = 3.52 L (27.1 %CV)
    lk12 <- log(5.86);   label("Central-to-peripheral rate constant k12 (1/h)")                           # Table 1 Dox (iv): k12 = 5.86 1/h (19.1 %CV)
    lk21 <- log(0.127);  label("Peripheral-to-central rate constant k21 (1/h)")                           # Table 1 Dox (iv): k21 = 0.127 1/h (19.1 %CV)
    lkel <- log(1.24);   label("Elimination rate constant from the central compartment ke (1/h)")          # Table 1 Dox (iv): ke = 1.24 1/h (37.1 %CV)

    # --- Tumor growth and inhibition: Choi 2021 Table 2, 'Dox' column -------
    # L0 and L1 were estimated from the vehicle-control arm and held fixed
    # in the doxorubicin fit (Table 2 footnote a).
    ltumorExpGrowth <- fixed(log(0.107)); label("Exponential-phase tumor growth rate L0 (1/day)")           # Table 2 Dox: L0 = 0.107 1/day (fixed; Control estimate, 8.40 %CV)
    ltumorLinGrowth <- fixed(log(0.148)); label("Linear-phase tumor growth rate L1 (cm^3/day)")             # Table 2 Dox: L1 = 0.148 cm^3/day (fixed; Control estimate, 0.809 %CV)
    lrbase_tumor    <- log(0.0386);       label("Initial tumor volume w0 (cm^3)")                           # Table 2 Dox: w0 = 0.0386 cm^3 (10.9 %CV)
    ldamageTransit  <- log(0.130);        label("Damaged-cell transit rate constant k1 (1/day)")            # Table 2 Dox: k1 = 0.130 1/day (6.00 %CV)
    ldrugSlope      <- log(13.8);         label("Doxorubicin potency k2 (L/mg/day)")                        # Table 2 Dox: k2 = 13.8 L/mg/day (8.91 %CV)

    # --- Residual error ------------------------------------------------------
    # Choi 2021 Methods: 'the exponential model showing better goodness of
    # fit was used in final model'. The error-model magnitude is not printed;
    # it is derived as sqrt(MSE) of the reported fit (the MSE of a
    # least-squares fit on log residuals is the residual-variance estimate;
    # Table 2 confirms MSE = SSE / (n - p), e.g. 4.72 / 0.148 = 32 = 35 - 3).
    expSd           <- sqrt(0.290); label("Exponential residual error on doxorubicin plasma concentration (SD, log scale)")  # derived: sqrt(Table 1 Dox (iv) MSE = 0.290)
    expSd_tumor_vol <- sqrt(0.148); label("Exponential residual error on tumor volume (SD, log scale)")                      # derived: sqrt(Table 2 Dox MSE = 0.148)
  })

  model({
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    # Choi 2021 prints the PK rate constants in 1/h; x24 converts to 1/day.
    k12 <- exp(lk12) * 24
    k21 <- exp(lk21) * 24
    kel <- exp(lkel) * 24

    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    rbase_tumor    <- exp(lrbase_tumor)
    damageTransit  <- exp(ldamageTransit)
    drugSlope      <- exp(ldrugSlope)

    # --- Doxorubicin PK (Choi 2021 Equations 1-2; amounts in mg) ------------
    d/dt(central)     <- -k12 * central - kel * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    Cc <- central / vc

    # --- Tumor growth inhibition (Choi 2021 Equations 8-14) -----------------
    # Natural growth follows Koch 2009: 2*L0*L1*x1^2 / ((L1 + 2*L0*x1) * w),
    # which is exponential (rate 2*L0*x1/w) for a small tumor and linear
    # (rate L1*x1/w) for a large one. Drug kills proliferating cells at
    # k2*C and damaged cells pass three transit stages at k1.
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3
    drugEffectCyclingCells <- drugSlope * Cc

    d/dt(cycling_cells)  <- 2 * tumorExpGrowth * tumorLinGrowth * cycling_cells^2 /
      ((tumorLinGrowth + 2 * tumorExpGrowth * cycling_cells) * tumor_vol) -
      drugEffectCyclingCells * cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    cycling_cells(0)  <- rbase_tumor
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0

    Cc ~ lnorm(expSd)
    tumor_vol ~ lnorm(expSd_tumor_vol)
  })
}
