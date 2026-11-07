Yang_2025_afatinib_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; PK in female BALB/c nude mice, tumor growth in BALB/c",
    "nude mice bearing gefitinib-resistant HCC827 (HCC827-GR) non-small-cell",
    "lung cancer xenografts). Afatinib one-compartment oral PK model with",
    "first-order absorption and elimination driving a Simeoni tumor-growth-",
    "inhibition model (smooth exponential-then-linear natural growth, linear",
    "kill of proliferating cells, three-stage damaged-cell transit chain). The",
    "herbal mixture HangAmDan-B1 (HAD-B1, 400 mg/kg/day) has no PK component",
    "and enters only as a binary synergy factor that multiplies afatinib",
    "potency by (1 + alpha_syn)."
  )
  reference <- paste(
    "Yang SY, Ngo LT, Lee S, Yun HY, Bui TT, Kim DH, Chae JW, Park S.",
    "Model-Based Evaluation of HangAmDan-B1 and Afatinib Combination Therapy",
    "in HCC827 Xenograft Mice with Resistance to Epidermal Growth Factor",
    "Receptor Tyrosine Kinase Inhibitor. Pharmaceuticals (Basel).",
    "2025;18(5):748. doi:10.3390/ph18050748.",
    sep = " "
  )
  vignette <- "Yang_2025_afatinib_mouse"

  units <- list(time = "day", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "afatinib", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "afatinib", units = "mg/kg", specimen = "serum", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list(
    CONMED_HADB1 = list(
      description = "HangAmDan-B1 (HAD-B1) coadministration indicator: 1 = afatinib plus HAD-B1 400 mg/kg/day (combination group), 0 = afatinib monotherapy or vehicle control",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Yang 2025 Equation (1) and Equations (10)-(11): afatinib potency is",
        "KAFT * (1 + alpha_syn * HAD-B1i). HAD-B1 serum constituents were not",
        "detectable, so the indicator carries the whole HAD-B1 contribution and",
        "was estimated at a single HAD-B1 dose (400 mg/kg/day); it has no dose",
        "dimension and no effect on afatinib PK (Section 2.1.1: no significant",
        "difference in afatinib AUC or Cmax with HAD-B1)."
      ),
      source_name = "HAD-B1i"
    )
  )

  population <- list(
    species = "mouse (BALB/c nude; PK in tumor-free females, tumor growth with subcutaneous gefitinib-resistant HCC827 xenograft)",
    n_subjects = 36L,
    n_studies = 2L,
    age_range = "5 weeks at purchase (PK mice)",
    weight_range = "18-20 g at purchase (PK mice)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA,
    disease_state = "PK: tumor-free mice. PD: subcutaneous right-flank xenograft of gefitinib-resistant HCC827 human non-small-cell lung cancer cells (HCC827-GR; 1e6 cells), randomised at 80-100 mm^3.",
    dose_range = "Afatinib 10 mg/kg orally once daily for 30 days, alone or with HAD-B1 400 mg/kg orally once daily; tumor study also had vehicle-control and HAD-B1-alone arms.",
    regions = "Republic of Korea (preclinical)",
    notes = paste(
      "PK: n = 6 female BALB/c nude mice per arm (afatinib alone, afatinib +",
      "HAD-B1), serum sampled 0.5, 2, 4 and 8 h after the day-30 dose (Section",
      "4.2.2). PD: tumor-volume data from Song et al. 2022 (n = 6 mice per",
      "group: vehicle, afatinib, HAD-B1, afatinib + HAD-B1); for modeling, 100",
      "tumor-growth profiles per group were generated from the published group",
      "mean and SD assuming normality (Section 4.3.1). Natural growth was fitted",
      "to the control group first; lambda0 and lambda1 were then fixed and the",
      "treated groups were pooled for the final PK/PD fit (Table 3, 'Treated",
      "Groups'), which is the parameter set encoded here."
    )
  )

  ini({
    # --- Afatinib PK: Yang 2025 Table 2 -----------------------------------
    # Table 2 prints Ka in 1/h and CL/F in L/h/kg; model time is days, so
    # model() multiplies both by 24.
    lka <- log(0.2); label("First-order absorption rate constant Ka (1/h)")          # Table 2: Ka = 0.2 1/h (RSE 39%)
    lcl <- log(24.5); label("Apparent clearance CL/F (L/h/kg)")                      # Table 2: CL/F = 24.5 L/h/kg (RSE 9%)
    lvc <- log(56.1); label("Apparent volume of distribution Vd/F (L/kg)")           # Table 2: Vd/F = 56.1 L/kg (RSE 17%)

    # Table 2 'IIV on CL/F' = 30.9% and 'IIV on Vd/F' = 50.2%, read as CVs:
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.091194 # Table 2: IIV on CL/F = 30.9% -> log(1 + 0.309^2)
    etalvc ~ 0.224745 # Table 2: IIV on Vd/F = 50.2% -> log(1 + 0.502^2)

    propSd <- 0.29; label("Proportional residual error on afatinib serum concentration (fraction)") # Table 2: Prop.ER = 29.0%

    # --- Tumor growth inhibition: Yang 2025 Table 3, 'Treated Groups' -----
    ltumorExpGrowth <- fixed(log(0.156)); label("Exponential tumor growth rate lambda0 (1/day)")            # Table 3: lambda0 = 0.156 1/day (FIX for treated groups, from the control fit)
    ltumorLinGrowth <- fixed(log(34.9)); label("Linear tumor growth rate lambda1 (mm^3/day)")               # Table 3: lambda1 = 34.9 mm^3/day (FIX for treated groups); Section 2.2.1 prose prints 35.5
    psi <- fixed(20); label("Shape factor phi for the exponential-to-linear switch (unitless)")              # Table 3: phi = 20 (FIX); Section 4.3.2: fixed per Simeoni 2004
    lrbase_tumor <- log(32.7); label("Tumor volume at inoculation V0 (mm^3)")                               # Table 3 Treated Groups: V0 = 32.7 mm^3
    ldamageTransit <- log(2.2); label("Damaged-cell transit (cell death) rate constant K1 (1/day)")         # Table 3 Treated Groups: K1 = 2.2 1/day
    ldrugSlope <- log(0.0034); label("Afatinib potency KAFT (mL/ng/day)")                                    # Table 3 Treated Groups: KAFT = 0.0034 mL/ng/day
    e_conmed_hadb1_drugSlope <- 0.45; label("Fractional increase in afatinib potency with HAD-B1 alpha_syn (unitless)") # Table 3 Treated Groups: alpha_syn = 0.45

    propSd_tumor_vol <- 0.215; label("Proportional residual error on tumor volume (fraction)") # Table 3 Treated Groups: Prop.ER = 0.215
  })

  model({
    # Table 2 rate units are per hour; x24 converts to per day.
    ka <- exp(lka) * 24
    cl <- exp(lcl + etalcl) * 24
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    rbase_tumor <- exp(lrbase_tumor)
    damageTransit <- exp(ldamageTransit)
    drugSlope <- exp(ldrugSlope)

    # --- Afatinib PK (Yang 2025 Equations 2-3; amounts in mg/kg) -----------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # mg/kg / (L/kg) = mg/L; x1000 gives ng/mL, the scale of KAFT.
    Cc <- 1000 * central / vc

    # --- Simeoni tumor growth inhibition (Yang 2025 Equations 10-16) -------
    # Equation 10 prints the growth-switch term on the proliferating pool
    # X1(t), not on total tumor volume as in Simeoni 2004; encoded as printed.
    drugEffectCyclingCells <- drugSlope * (1 + e_conmed_hadb1_drugSlope * CONMED_HADB1) * Cc

    d/dt(cycling_cells) <- tumorExpGrowth * cycling_cells /
      (1 + (tumorExpGrowth / tumorLinGrowth * cycling_cells)^psi)^(1 / psi) -
      drugEffectCyclingCells * cycling_cells
    d/dt(damaged_cells1) <- drugEffectCyclingCells * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    cycling_cells(0) <- rbase_tumor
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0

    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    Cc ~ prop(propSd)
    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
