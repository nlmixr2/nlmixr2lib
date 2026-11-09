Jiang_2025_fanastomig <- function() {
  description <- "Two-compartment population PK model for fanastomig (EMB-02; anti-PD-1 x LAG-3 bispecific antibody) with sigmoid I_max time-dependent clearance and an albumin effect on central volume, plus a direct Emax model of peripheral-blood PD-1 receptor occupancy, in adults with advanced solid tumours (Jiang 2025)"
  reference <- paste(
    "Jiang C, Ren F, Zhang M, Lu Q, Zeng S, Yang G, Zhu Y. Using",
    "pharmacokinetic and pharmacodynamic analysis to optimize the dosing",
    "regimens of fanastomig (EMB-02) in patients with advanced solid tumors.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(5):975-986.",
    "doi:10.1002/psp4.70011.",
    sep = " "
  )
  vignette <- "Jiang_2025_fanastomig"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "fanastomig", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "fanastomig", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect (ALB / 36)^e_alb_vc on the central volume. The reference 36 g/L is the denominator in the final Phoenix model code (Data S1 'stparm(V = tvV *(ALB/36)^dVdALB * exp(nV))'); Table S1 lists a median of 38 g/L for all 47 screened patients, of whom 35 entered the final model. Declared as fcovariate (carried forward) in the model code; the paper does not say whether it varied over time, and Table S1 reports the baseline value.",
      source_name = "ALB"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 35L,
    n_studies = 1L,
    age_range = "37-83 years",
    age_median = "61 years",
    weight_range = "42.1-107 kg",
    weight_median = "72.0 kg",
    sex_female_pct = 59.6,
    race_ethnicity = c(White = 72.3, Asian = 21.3, Black = 4.3, Other = 2.1),
    disease_state = "Advanced solid tumours (first-in-human phase 1, NCT04618393)",
    dose_range = "6, 20, 60, 180, 360, 600 and 900 mg IV QW (Part 1 dose escalation); 60, 180 and 600 mg IV QW (Part 2 dose enrichment)",
    regions = "China, United States, Australia",
    albumin_median = "38.0 g/L (range 23.0-49.6)",
    ada_positive_pct = 95.7,
    notes = "Demographics are Jiang 2025 Table S1 (all 47 patients with PK data). Twelve patients pretreated with a PD-1 inhibitor within 4 months (249 observations, elevated predose anti-PD-(L)1 concentrations) were excluded from model development, leaving 35 patients and 625 observations (Results, first paragraph). PD-1 receptor occupancy was measured by flow cytometry on peripheral-blood CD3+ T cells."
  )

  ini({
    # Structural PK parameters at t = 0 (Jiang 2025 Table 1, final model).
    # The Phoenix model code in Data S1 carries a re-run with slightly
    # different digits (V 5198.95 mL, Cl0 214.70 mL/h, V2 3155.52 mL,
    # Cl2 46.17 mL/h); the Table 1 published estimates are used here.
    lcl <- log(0.215); label("Clearance at t = 0 (L/h)") # Table 1: CL = 0.215 L/h
    lvc <- log(5.20); label("Central volume of distribution at ALB = 36 g/L (L)") # Table 1: Vc = 5.20 L
    lq <- log(0.0463); label("Intercompartmental clearance (L/h)") # Table 1: Q = 0.0463 L/h
    lvp <- log(3.15); label("Peripheral volume of distribution (L)") # Table 1: Vp = 3.15 L

    # Time-dependent clearance, Jiang 2025 Equation 1 and Data S1:
    #   CL(t) = CL0 * exp(Imax * t^Hill / (t^Hill + T50^Hill))
    # Imax is negative; Data S1 'stparm(Imax = tvImax * exp(nImax))' puts a
    # log-normal eta on a negative typical value, so the individual value
    # keeps the typical sign. Stored as log|Imax| with the sign applied in
    # model(), which is the same parameterisation exactly. Maximum reduction
    # 1 - exp(-0.954) = 61.5 %, as stated in the Results.
    lcl_time_max <- log(0.954); label("log|Imax|; magnitude of the log-CL change at t >> T50 (unitless)") # Table 1: Imax = -0.954
    lcl_t50 <- log(180); label("Time at half of the maximal change in CL (h)") # Table 1: T50 = 180 h
    lcl_time_hill <- log(8.49); label("Hill coefficient of the time-on-CL function (unitless)") # Table 1: Hill = 8.49

    # Albumin power effect on the central volume (see covariateData$ALB).
    e_alb_vc <- -0.924; label("Power exponent of ALB/36 on Vc (unitless)") # Table 1: Effect of ALB = -0.924; Data S1 applies it to V (central)

    # PD-1 receptor occupancy, direct Emax model on serum concentration.
    lemax <- log(92.5); label("Maximal peripheral PD-1 receptor occupancy (%)") # digitised by the maintainers from the Figure 3 median curve (92.2 % at 28.5 ug/mL); not printed in the paper
    lec50 <- log(0.084); label("Serum concentration at half-maximal PD-1 receptor occupancy (mg/L)") # Abstract and Results 'Receptor Occupancy': EC50 = 0.084 ug/mL

    # IIV: diagonal, exponential (Table 1 variances; Data S1 ranef diag).
    etalvc ~ 0.035 # Table 1: omega2 Vc = 0.035
    etalcl ~ 0.068 # Table 1: omega2 Cl = 0.068
    etalcl_time_max ~ 0.416 # Table 1: omega2 Imax = 0.416
    etalcl_t50 ~ 0.394 # Table 1: omega2 T50 = 0.394

    # Residual error: proportional, reported as an SD (Data S1 CEps = 0.378).
    propSd <- 0.378; label("Proportional residual error (fraction)") # Table 1: Prop sigma = 0.378
  })

  model({
    vc <- exp(lvc + etalvc) * (ALB / 36)^e_alb_vc
    vp <- exp(lvp)
    q <- exp(lq)

    # Time-dependent CL; t is time since the first dose (h).
    cl_time_max_i <- -exp(lcl_time_max + etalcl_time_max)
    cl_t50 <- exp(lcl_t50 + etalcl_t50)
    cl_time_hill <- exp(lcl_time_hill)
    cl0 <- exp(lcl + etalcl)
    cl <- cl0 * exp(cl_time_max_i * t^cl_time_hill / (t^cl_time_hill + cl_t50^cl_time_hill))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg / volume in L = mg/L (= ug/mL, the paper's reporting unit).
    Cc <- central / vc

    # Peripheral-blood PD-1 receptor occupancy (%), CD3+ T cells. Hyperbolic
    # (no Hill term): the paper's 90 %-of-maximal target 0.756 ug/mL is
    # exactly 9 x EC50. Returned as a derived output; the paper reports no
    # IIV or residual error for this sub-model.
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    RO <- emax * Cc / (ec50 + Cc)

    Cc ~ prop(propSd)
  })
}
