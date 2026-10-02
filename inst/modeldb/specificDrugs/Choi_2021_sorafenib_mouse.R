Choi_2021_sorafenib_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; PK in male athymic nude mice, tumor growth in female",
    "SCID mice bearing orthotopic human 143B osteosarcoma xenografts).",
    "Sorafenib monotherapy PK-PD model (Choi 2021 model B): a two-compartment",
    "linear PK model with first-order oral absorption, parameterised with",
    "micro-constants (disposition from the IV arm, ka from the oral arm),",
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
    depot = list(analyte = "sorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sorafenib", units = "mg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "drug-damaged tumor cells", units = "cm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (PK: male athymic nude Foxn1nu; tumor growth: female CB17 SCID with orthotopic 143B osteosarcoma xenograft)",
    n_subjects = 26L,
    n_studies = 2L,
    age_range = "7 weeks at purchase (PK mice)",
    weight_range = "approximately 30 g (PK mice)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "PK: tumor-free mice. PD: orthotopic (intratibial) human 143B osteosarcoma cell-line-derived xenograft.",
    dose_range = "PK: single 0.02 mg IV bolus or single 0.2 mg oral gavage per mouse. PD: 200 ug/mouse oral gavage on 7 days of cycle 1 then every other day x4 in cycle 2, per Jian 2017.",
    regions = "University of California Davis, USA (preclinical)",
    notes = paste(
      "PK: n = 6 male athymic nude mice per route (IV to 48 h, oral to 96 h;",
      "Choi 2021 Methods, PK Studies; Table 1). The IV data gave Vc, Vp, k12,",
      "k21 and ke; the oral data were then fitted for ka alone with the IV",
      "disposition fixed, so oral bioavailability is implicitly 1. PD: n = 7",
      "sorafenib-treated tumor-bearing mice from Jian 2017 (Choi 2021 Table 2);",
      "natural growth rates L0 and L1 were estimated from n = 7 vehicle-control",
      "mice and held fixed."
    )
  )

  ini({
    # --- Sorafenib PK: Choi 2021 Table 1 -----------------------------------
    # Disposition from the 'Sor (iv)' column; ka from the 'Sor (p.o.)' column
    # (Table 1 footnote a: only ka was estimated for the oral data). Rate
    # constants are printed in 1/h; model time is days, so model()
    # multiplies each rate constant by 24.
    lka  <- log(0.640);  label("First-order oral absorption rate constant ka (1/h)")                      # Table 1 Sor (p.o.): ka = 0.640 1/h (18.8 %CV)
    lvc  <- log(0.0112); label("Central volume of distribution Vc (L)")                                   # Table 1 Sor (iv): Vc = 0.0112 L (9.92 %CV)
    lvp  <- log(0.0120); label("Peripheral volume of distribution Vp (L)")                                # Table 1 Sor (iv): Vp = 0.0120 L (13.7 %CV)
    lk12 <- log(0.663);  label("Central-to-peripheral rate constant k12 (1/h)")                           # Table 1 Sor (iv): k12 = 0.663 1/h (4.06 %CV)
    lk21 <- log(0.620);  label("Peripheral-to-central rate constant k21 (1/h)")                           # Table 1 Sor (iv): k21 = 0.620 1/h (9.19 %CV)
    lkel <- log(0.344);  label("Elimination rate constant from the central compartment ke (1/h)")          # Table 1 Sor (iv): ke = 0.344 1/h (31.4 %CV)

    # --- Tumor growth and inhibition: Choi 2021 Table 2, 'Sor' column -------
    # L0 and L1 were estimated from the vehicle-control arm and held fixed
    # in the sorafenib fit (Table 2 footnote a).
    ltumorExpGrowth <- fixed(log(0.107)); label("Exponential-phase tumor growth rate L0 (1/day)")           # Table 2 Sor: L0 = 0.107 1/day (fixed; Control estimate, 8.40 %CV)
    ltumorLinGrowth <- fixed(log(0.148)); label("Linear-phase tumor growth rate L1 (cm^3/day)")             # Table 2 Sor: L1 = 0.148 cm^3/day (fixed; Control estimate, 0.809 %CV)
    lrbase_tumor    <- log(0.0415);       label("Initial tumor volume w0 (cm^3)")                           # Table 2 Sor: w0 = 0.0415 cm^3 (12.9 %CV)
    ldamageTransit  <- log(17.0);         label("Damaged-cell transit rate constant k1 (1/day)")            # Table 2 Sor: k1 = 17.0 1/day (36.4 %CV)
    ldrugSlope      <- log(0.0267);       label("Sorafenib potency k2 (L/mg/day)")                          # Table 2 Sor: k2 = 0.0267 L/mg/day (35.8 %CV)

    # --- Residual error ------------------------------------------------------
    # Choi 2021 Methods: exponential error model in the final model; the
    # magnitude is not printed and is derived as sqrt(MSE) of the reported
    # fit. The oral-arm MSE is used for the concentration because the
    # tumor-growth study dosed orally (the IV-arm MSE is 0.0741).
    expSd           <- sqrt(0.337); label("Exponential residual error on sorafenib plasma concentration (SD, log scale)")  # derived: sqrt(Table 1 Sor (p.o.) MSE = 0.337)
    expSd_tumor_vol <- sqrt(0.149); label("Exponential residual error on tumor volume (SD, log scale)")                    # derived: sqrt(Table 2 Sor MSE = 0.149)
  })

  model({
    ka  <- exp(lka) * 24
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

    # --- Sorafenib PK (Choi 2021 Equations 3-5; amounts in mg) --------------
    # Oral doses go to depot; IV doses go to central. Bioavailability is
    # implicitly 1 (the oral fit estimated ka only).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - k12 * central - kel * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    Cc <- central / vc

    # --- Tumor growth inhibition (Choi 2021 Equations 8-14) -----------------
    # Natural growth follows Koch 2009: 2*L0*L1*x1^2 / ((L1 + 2*L0*x1) * w).
    # Drug kills proliferating cells at k2*C and damaged cells pass three
    # transit stages at k1.
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
