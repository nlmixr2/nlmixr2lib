Chua_2023_mirikizumab_phase2 <- function() {
  description <- "Two-compartment population PK model with first-order SC absorption for mirikizumab (anti-IL-23p19 mAb) in adults with moderately to severely active ulcerative colitis (Chua 2023, phase II AMAC model)"
  reference <- "Chua L, Friedrich S, Zhang XC. Mirikizumab Pharmacokinetics in Patients with Moderately to Severely Active Ulcerative Colitis: Results from Phase III LUCENT Studies. Clin Pharmacokinet. 2023;62(10):1479-1491. doi:10.1007/s40262-023-01281-z"
  vignette <- "Chua_2023_mirikizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "mirikizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mirikizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "mirikizumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value. Power effect on Vc only (exponent 0.247), normalised to 75.4 kg (Table 2 footnote d).",
      source_name = "baseline weight"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (ESM Table 2: only time-varying albumin was tested for AMAC). Power effect on CL normalised to 43.2 g/L (Table 2 footnote c).",
      source_name = "albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 233L,
    n_studies = 1L,
    n_observations = 4103L,
    age_range = "17-72 years",
    age_mean = "42.2 years",
    weight_range = "39.6-122.0 kg",
    weight_mean = "74.5 kg",
    bmi_mean = "25.2 kg/m^2 (range 15.0-42.7)",
    albumin_mean = "41.1 g/L (range 24-51)",
    sex_female_pct = 40.8,
    race_ethnicity = c(White = 83.3, Asian = 13.7, Black = 3.0),
    disease_state = "Moderately to severely active ulcerative colitis (total Mayo score 6-12, Mayo endoscopic subscore >= 2) with inadequate response or intolerance to conventional or biologic therapy; 35.2% biologic-naive.",
    dose_range = "Induction 50, 200 or 600 mg IV Q4W x 3 (50 and 200 mg arms exposure-adjusted upward, never above 600 mg; group averages 100 and 250 mg); open-label extension 600 or 1000 mg IV Q4W x 3; maintenance 200 mg SC Q4W or Q12W (lyophilised formulation).",
    regions = "85 sites in 14 countries",
    notes = "AMAC (NCT02589665) phase II trial; Chua 2023 Table 1. BLQ samples (about 4%) excluded. Fit in NONMEM 7.4.2 separately from the phase III LUCENT 1+2 model (see Chua_2023_mirikizumab_phase3)."
  )

  ini({
    # Typical values at the reference subject (75.4 kg, albumin 43.2 g/L).
    # Rates are published per hour; x 24 converts to per day.
    lka <- log(0.00757 * 24); label("SC absorption rate constant (1/day)") # Table 2 AMAC: Ka 0.00757 1/h
    lcl <- log(0.0230 * 24); label("Clearance (L/day)") # Table 2 AMAC: CL 0.0230 L/h
    lvc <- log(3.40); label("Central volume of distribution (L)") # Table 2 AMAC: V1 3.40 L
    lvp <- log(1.58); label("Peripheral volume of distribution (L)") # Table 2 AMAC: V2 1.58 L
    lq <- log(0.00756 * 24); label("Intercompartmental clearance (L/day)") # Table 2 AMAC: Q 0.00756 L/h
    lfdepot <- log(0.422); label("SC bioavailability (fraction)") # Table 2 AMAC: F1 42.2%

    e_wt_vc <- 0.247; label("Power exponent of baseline weight on Vc (unitless)") # Table 2 AMAC and footnote d
    e_alb_cl <- -0.951; label("Power exponent of albumin on CL (unitless)") # Table 2 AMAC and footnote c

    # IIV: Table 2 footnote b, CV% = 100 * sqrt(exp(OMEGA) - 1), so
    # OMEGA = log(1 + CV^2). Covariance CL-V1 printed directly as 0.0426.
    # The AMAC F1 row carries no footnote e (logit domain), so footnote b
    # (log-normal) applies to its IIV.
    etalcl + etalvc ~ c(
      0.094605,
      0.0426, 0.060625
    ) # Table 2 AMAC: IIV CL 31.5%, V1 25.0%, covariance CL-V1 0.0426
    etalvp ~ 0.179932 # Table 2 AMAC: IIV V2 44.4%
    etalfdepot ~ 0.111842 # Table 2 AMAC: IIV F1 34.4%

    propSd <- 0.254; label("Proportional residual error (fraction)") # Table 2 AMAC: 25.4%
    addSd <- 0.0987; label("Additive residual error (ug/mL)") # Table 2 AMAC: 98.7 ng/mL
  })
  model({
    wt_vc <- (WT / 75.4)^e_wt_vc
    alb_cl <- (ALB / 43.2)^e_alb_cl

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * alb_cl
    vc <- exp(lvc + etalvc) * wt_vc
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
