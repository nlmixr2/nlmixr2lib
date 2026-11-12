Yang_2025b_iruplinalkib_mouse_ros1 <- function() {
  description <- paste(
    "Preclinical (mouse; female BALB/c-nu bearing subcutaneous SLC34A2-ROS1",
    "Ba/F3 cell-line-derived xenografts). Iruplinalkib PK-PD tumor growth",
    "inhibition model: a two-compartment PK model with first-order oral",
    "absorption, fitted first to the plasma data of the same mice and then",
    "held fixed, drives a modified Simeoni model in which exponential-then-",
    "linear tumor growth (switch sharpness psi = 20) is opposed by an Emax",
    "kill of proliferating cells (Kmax, KC50) followed by three damaged-cell",
    "transit compartments. Typical-value model."
  )

  reference <- paste(
    "Yang G, Wang Y, Zhao H, Zheng Q, Wang X, Jia L, Xin Q, Ma C, Zhang Y,",
    "Zheng S, Kang X.",
    "Translational Model-Informed Dose Selection for Iruplinalkib, a",
    "Selective Oral ALK/ROS1 Tyrosine Kinase Inhibitor.",
    "Clin Transl Sci. 2025;18(7):e70287.",
    "doi:10.1111/cts.70287"
  )

  vignette <- "Yang_2025b_iruplinalkib"

  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "iruplinalkib", units = "mg/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "iruplinalkib", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "iruplinalkib", units = "mg/kg", specimen = "tissue", verified = TRUE),
    cycling_cells = list(analyte = "proliferating tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells1 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells2 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE),
    damaged_cells3 = list(analyte = "drug-damaged tumor cells", units = "mm^3", specimen = "tumor", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (female BALB/c-nu, SLC34A2-ROS1 Ba/F3 xenograft)",
    n_subjects = 36L,
    n_studies = 1L,
    sex_female_pct = 100,
    disease_state = paste(
      "Subcutaneous (right armpit) xenografts of SLC34A2-ROS1-overexpressing",
      "Ba/F3 cells; treatment started when tumors reached 100-300 mm^3"
    ),
    dose_range = paste(
      "Oral gavage 0 (vehicle), 5, 10 or 20 mg/kg once daily for 14 days;",
      "n = 9 per group"
    ),
    notes = paste(
      "Tumor volume measured twice weekly. Plasma sampled predose and 1, 2,",
      "4, 8 and 24 h after the final (Day 14) dose and pooled for the PK fit",
      "(Table S1). Both PK and PK-PD fits used Phoenix NLME 8.3.4 with",
      "FOCE-ELS. The paper reports no between-animal variability for either",
      "stage and no residual error for the tumor-volume fit."
    )
  )

  ini({
    # Stage 1, mouse PK (Table S1). Estimated first, then 'fixed within the
    # PKPD model to fit the PD parameters' (Methods 'Mouse PKPD Modeling'),
    # hence fixed() here. Apparent oral parameters per kg of body weight.
    lka <- fixed(log(0.366)); label("First-order absorption rate constant Ka (1/h)")          # Table S1: Ka = 0.366 1/h (RSE 8.9%)
    lcl <- fixed(log(2.48)); label("Apparent oral clearance CL/F (L/h/kg)")                   # Table S1: CL/F = 2.48 L/h/kg (RSE 5.1%)
    lvc <- fixed(log(3.59)); label("Apparent central volume V1/F (L/kg)")                     # Table S1: V1/F = 3.59 L/kg (RSE 20%)
    lq <- fixed(log(0.0843)); label("Apparent intercompartmental clearance Q/F (L/h/kg)")     # Table S1: Q/F = 0.0843 L/h/kg (RSE 37%)
    lvp <- fixed(log(1.65)); label("Apparent peripheral volume V2/F (L/kg)")                  # Table S1: V2/F = 1.65 L/kg (RSE 15%)

    # Stage 2, modified Simeoni tumor growth inhibition (Table 1, ROS1
    # positive model; equations Methods S1 Eq1-7).
    ltumorExpGrowth <- log(0.00851); label("First-order (exponential-phase) tumor growth rate lambda0 (1/h)")  # Table 1 ROS1: lambda0 = 0.00851 1/h (RSE 7.2%)
    ltumorLinGrowth <- log(4.45); label("Zero-order (linear-phase) tumor growth rate lambda1 (mm^3/h)")        # Table 1 ROS1: lambda1 = 4.45 mm3/h (RSE 2.1%)
    psi <- fixed(20); label("Exponential-to-linear growth switch sharpness psi (unitless)")                    # Table 1 ROS1: psi = 20 (FIX)
    lkmax <- log(0.0739); label("Maximal tumor cell killing rate Kmax (1/h)")                                 # Table 1 ROS1: Kmax = 0.0739 1/h (RSE 34%)
    lkc50 <- log(1221); label("Plasma concentration giving half of Kmax, KC50 (ng/mL)")                       # Table 1 ROS1: KC50 = 1221 ng/mL (RSE 44%)
    ldamageTransit <- log(0.168); label("Damaged-cell transit rate constant k1 (1/h)")                        # Table 1 ROS1: K1 = 0.168 1/h (RSE 52%)
    lrbase_tumor <- log(221); label("Tumor volume at the start of the model time axis E0 (mm^3)")             # Table 1 ROS1: E0 = 221 mm3 (RSE 5.0%)

    # Residual error. Plasma: proportional (Methods), estimated in stage 1
    # and carried fixed into stage 2. Tumor volume: not reported in any
    # source; encoded fixed(0) (see vignette Errata).
    propSd <- fixed(0.46); label("Proportional residual error on plasma iruplinalkib (fraction)")               # Table S1: Proportional residual = 46% (RSE 5.5%)
    propSd_tumor_vol <- fixed(0); label("Proportional residual error on tumor volume (fraction), not reported")
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    tumorExpGrowth <- exp(ltumorExpGrowth)
    tumorLinGrowth <- exp(ltumorLinGrowth)
    kmax <- exp(lkmax)
    kc50 <- exp(lkc50)
    damageTransit <- exp(ldamageTransit)
    rbase_tumor <- exp(lrbase_tumor)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # (mg/kg) / (L/kg) = mg/L; x1000 -> ng/mL, the units of KC50.
    Cc <- central / vc * 1000

    # Methods S1 Eq6: E = E1 + E2 + E3 + E4.
    tumor_vol <- cycling_cells + damaged_cells1 + damaged_cells2 + damaged_cells3

    # Methods S1 Eq2-5. Growth of the proliferating compartment uses the
    # TOTAL volume in the switch term; Emax kill on the proliferating cells.
    growth <- tumorExpGrowth * cycling_cells /
      (1 + (tumorExpGrowth / tumorLinGrowth * tumor_vol)^psi)^(1 / psi)
    killRate <- kmax * Cc / (kc50 + Cc)

    d/dt(cycling_cells) <- growth - killRate * cycling_cells
    d/dt(damaged_cells1) <- killRate * cycling_cells - damageTransit * damaged_cells1
    d/dt(damaged_cells2) <- damageTransit * (damaged_cells1 - damaged_cells2)
    d/dt(damaged_cells3) <- damageTransit * (damaged_cells2 - damaged_cells3)

    # Methods S1 Eq7.
    cycling_cells(0) <- rbase_tumor
    damaged_cells1(0) <- 0
    damaged_cells2(0) <- 0
    damaged_cells3(0) <- 0

    Cc ~ prop(propSd)
    tumor_vol ~ prop(propSd_tumor_vol)
  })
}
