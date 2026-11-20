Li_2026_meropenem <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous meropenem (given as",
    "the 1:1 meropenem-vaborbactam fixed-dose combination by 3-h infusion) in",
    "healthy Chinese adults, with a fixed non-renal clearance plus a sigmoidal",
    "Hill renal clearance driven by time-varying creatinine clearance, and a",
    "power effect of body weight on the peripheral volume (Li 2026)"
  )
  reference <- paste(
    "Li W, Wei Q, Chen M, Liu X, Guo B, Chen Y, Zhang Q, Guan W, Shen L,",
    "Wang X, Li N, Chen Y, Li X, Hu Y, Zhang J, Wu X. Pharmacokinetics,",
    "safety and population pharmacokinetics/pharmacodynamics analysis of",
    "meropenem-vaborbactam in Chinese healthy participants. Front Pharmacol.",
    "2026;17:1887991. doi:10.3389/fphar.2026.1887991"
  )
  vignette <- "Li_2026_meropenem_vaborbactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance in raw mL/min (NOT BSA-normalised); time-varying.",
        "The estimating equation is not stated by the source."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results, 'Population PK': 'time-varying creatinine clearance (CRCL)",
        "was included as a covariate on renal clearance', entering through the",
        "sigmoidal Hill term CLR = CLRE * CRCL^HILL / (CR50^HILL + CRCL^HILL)",
        "with no normalising value. Table 5 labels CR50 in mL/min. The cohort",
        "were healthy adults with eGFR >= 90 mL/min; supplement Table S5",
        "subgroups put CRCL at about 85-170 mL/min (median about 105)."
      ),
      source_name = "CRCL"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Vp as (WT / 65.6)^1.07 (Results equation; Table 5). 65.6 kg",
        "is the Table 2 pooled median (range 50.0-73.2 kg)."
      ),
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 20L,
    n_studies = 1L,
    age_range = "26-43 years",
    age_median = "33.0 years",
    weight_range = "50.0-73.2 kg",
    weight_median = "65.6 kg",
    sex_female_pct = 50,
    race_ethnicity = "100% Asian (Chinese)",
    disease_state = "Healthy adults (eGFR >= 90 mL/min at screening)",
    dose_range = paste(
      "Meropenem-vaborbactam 1 g/1 g or 2 g/2 g (10 per group) as 3-h",
      "intravenous infusions: a single dose on day 1, every 8 h on days 2-4,",
      "and a single final dose on day 5."
    ),
    regions = "China (Phase I Unit, Huashan Hospital, Fudan University)",
    renal_function = "Normal; eGFR >= 90 mL/min was an inclusion criterion.",
    notes = paste(
      "570 meropenem plasma concentrations from 20 subjects (Results,",
      "'Population PK'). Single-centre, open-label phase 1 study",
      "CTR20243815. NONMEM 7.5.0 with PsN 5.2.6. Demographics from Table 2."
    )
  )

  ini({
    # Structural parameters -- Li 2026 Table 5 (final meropenem model) and the
    # Results equations. CL = (CLNR + CLR) * exp(eta_CL), with
    # CLR = CLRE * CRCL^HILL / (CR50^HILL + CRCL^HILL).
    lcl_nonren <- fixed(log(3.94))
    label("Non-renal clearance from the central compartment (L/h)") # Table 5, 'CLNR' = 3.94 L/h, held constant (no RSE)
    lcl_renal_max <- log(10.9)
    label("Maximal renal clearance (L/h)") # Table 5, 'Maximal renal clearance, CLRE' = 10.9 L/h (RSE 9.00%)
    lcrcl50 <- log(75.5)
    label("Creatinine clearance at half-maximal renal clearance (mL/min)") # Table 5, 'CR50' = 75.5 mL/min (RSE 4.67%)
    lhill <- log(3.51)
    label("Hill coefficient of the renal clearance-CRCL relationship (unitless)") # Table 5, 'HILL' = 3.51 (RSE 24.8%)
    lvc <- log(7.60)
    label("Central volume of distribution (L)") # Table 5, 'V C' = 7.60 L (RSE 6.16%)
    lq <- log(6.60)
    label("Intercompartmental clearance (L/h)") # Table 5, 'Q' = 6.60 L/h (RSE 14.3%)
    lvp <- log(4.77)
    label("Peripheral volume of distribution at 65.6 kg (L)") # Table 5, 'V P' = 4.77 L (RSE 8.28%)

    # Covariate effect.
    e_wt_vp <- 1.07
    label("Power exponent of body weight on Vp (unitless)") # Table 5, 'Effect of body weight on V P' = 1.07 (RSE 14.8%)

    # IIV. Table 5 reports %CV; omega^2 = log(CV^2 + 1).
    etalcl ~ 0.0065718 # Table 5, 'Interindividual variation of CL' 8.12 %CV
    etalvp ~ 0.0055794 # Table 5, 'Interindividual variation of V P' 7.48 %CV

    # Residual error (combined additive + proportional).
    propSd <- 0.108
    label("Proportional residual error (fraction)") # Table 5, 'Proportional residual error' = 10.8%
    addSd <- 0.0955
    label("Additive residual error (mg/L)") # Table 5, 'Additive residual error' = 0.0955 mg/L
  })

  model({
    # 1. Renal and non-renal clearance arms.
    hill <- exp(lhill)
    crcl50 <- exp(lcrcl50)
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal_max) * CRCL^hill / (crcl50^hill + CRCL^hill)

    # 2. Individual PK parameters. The eta multiplies total clearance (the
    #    published equation brackets CLNR + CLR before the exponential).
    cl <- (cl_nonren + cl_renal) * exp(etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 65.6)^e_wt_vp

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Zero-order intravenous infusion into central via the
    #    event-table rate or duration.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
