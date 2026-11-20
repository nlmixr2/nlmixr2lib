Li_2026_vaborbactam <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous vaborbactam (given",
    "as the 1:1 meropenem-vaborbactam fixed-dose combination by 3-h infusion)",
    "in healthy Chinese adults, with a fixed non-renal clearance plus a",
    "sigmoidal Hill renal clearance driven by time-varying estimated",
    "glomerular filtration rate, and a power effect of body weight on the",
    "central volume (Li 2026)"
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
    central = list(analyte = "vaborbactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vaborbactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate in raw mL/min (as labelled by",
        "the source; not stated to be BSA-normalised); time-varying. The",
        "estimating equation is not stated by the source."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results, 'Population PK': 'time-varying estimated glomerular",
        "filtration rate (eGFR) was included as a covariate on renal",
        "clearance', entering through the sigmoidal Hill term",
        "CLR = CLRE * eGFR^HILL / (CR50^HILL + eGFR^HILL) with no normalising",
        "value. Table 6 labels CR50 in mL/min. Supplement Table S6 subgroups",
        "put eGFR at about 98-204 mL/min (median about 118). This is a",
        "different renal-function column from the creatinine clearance that",
        "drives the sibling Li_2026_meropenem model."
      ),
      source_name = "eGFR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters Vc as (WT / 65.6)^1.51 (Results equation; Table 6). 65.6 kg",
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
      "588 vaborbactam plasma concentrations from 20 subjects (Results,",
      "'Population PK'). Single-centre, open-label phase 1 study",
      "CTR20243815. NONMEM 7.5.0 with PsN 5.2.6. Demographics from Table 2."
    )
  )

  ini({
    # Structural parameters -- Li 2026 Table 6 (final vaborbactam model) and
    # the Results equations. CL = (CLNR + CLR) * exp(eta_CL), with
    # CLR = CLRE * eGFR^HILL / (CR50^HILL + eGFR^HILL).
    lcl_nonren <- fixed(log(0.164))
    label("Non-renal clearance from the central compartment (L/h)") # Table 6, 'CLNR' = 0.164 L/h, held constant (no RSE)
    lcl_renal_max <- log(10.8)
    label("Maximal renal clearance (L/h)") # Table 6, 'Maximal renal clearance, CLRE' = 10.8 L/h (RSE 5.67%)
    lcrcl50 <- log(73.3)
    label("eGFR at half-maximal renal clearance (mL/min)") # Table 6, 'CR50' = 73.3 mL/min (RSE 3.83%)
    lhill <- log(3.26)
    label("Hill coefficient of the renal clearance-eGFR relationship (unitless)") # Table 6, 'HILL' = 3.26 (RSE 21.6%)
    lvc <- log(10.3)
    label("Central volume of distribution at 65.6 kg (L)") # Table 6, 'V C' = 10.3 L (RSE 4.13%)
    lq <- log(1.64)
    label("Intercompartmental clearance (L/h)") # Table 6, 'Q' = 1.64 L/h (RSE 18.3%)
    lvp <- log(2.53)
    label("Peripheral volume of distribution (L)") # Table 6, 'V P' = 2.53 L (RSE 10.4%)

    # Covariate effect.
    e_wt_vc <- 1.51
    label("Power exponent of body weight on Vc (unitless)") # Table 6, 'Effect of body weight on V C' = 1.51 (RSE 15.9%)

    # IIV. Table 6 reports %CV; omega^2 = log(CV^2 + 1).
    etalcl ~ 0.013138 # Table 6, 'Interindividual variation of CL' 11.5 %CV
    etalvc ~ 0.0058657 # Table 6, 'Interindividual variation of V C' 7.67 %CV

    # Residual error (combined additive + proportional).
    propSd <- 0.131
    label("Proportional residual error (fraction)") # Table 6, 'Proportional residual error' = 13.1%
    addSd <- 0.0689
    label("Additive residual error (mg/L)") # Table 6, 'Additive residual error' = 0.0689 mg/L
  })

  model({
    # 1. Renal and non-renal clearance arms. CRCL carries the source's eGFR.
    hill <- exp(lhill)
    crcl50 <- exp(lcrcl50)
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal_max) * CRCL^hill / (crcl50^hill + CRCL^hill)

    # 2. Individual PK parameters. The eta multiplies total clearance (the
    #    published equation brackets CLNR + CLR before the exponential).
    cl <- (cl_nonren + cl_renal) * exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 65.6)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

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
