Chua_2023_mirikizumab_phase3 <- function() {
  description <- "Two-compartment population PK model with first-order SC absorption for mirikizumab (anti-IL-23p19 mAb) in adults with moderately to severely active ulcerative colitis (Chua 2023, phase III LUCENT 1+2 final model)"
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
      notes = "Baseline (study-entry) value. Shared power effect on CL and Q (exponent 0.534) and shared power effect on Vc and Vp (exponent 0.525), both normalised to 70.7 kg (Table 2 footnotes c and d).",
      source_name = "baseline weight"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL normalised to 45.0 g/L (Table 2 footnote c). Treated as time-varying: ESM Table 2 lists both baseline and time-varying albumin as tested in LUCENT 1+2, and footnote c writes 'albumin' while explicitly writing 'baseline weight' for weight.",
      source_name = "albumin"
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on the logit of SC bioavailability (Table 2, slope -0.0474). The paper does not print the BMI centring value; 25 kg/m^2 (the rounded LUCENT 1+2 mean, 24.9) is used and reproduces the Figure 3e model curve (F about 59% at BMI 15 and about 19% at BMI 54).",
      source_name = "baseline BMI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1129L,
    n_studies = 2L,
    n_observations = 7578L,
    age_range = "18-79 years",
    age_mean = "42.6 years",
    weight_range = "34.0-152.0 kg",
    weight_mean = "72.6 kg",
    bmi_mean = "24.9 kg/m^2 (range 13.8-53.5)",
    albumin_mean = "42.8 g/L (range 21-54)",
    sex_female_pct = 39.1,
    race_ethnicity = c(
      White = 74.7,
      Asian = 22.2,
      Black = 1.1,
      `American Indian or Alaska Native` = 1.0,
      Multiple = 0.1,
      Missing = 0.9
    ),
    disease_state = "Moderately to severely active ulcerative colitis (modified Mayo score 4-9, Mayo endoscopic subscore >= 2) with inadequate response or intolerance to conventional, biologic or JAK-inhibitor therapy; 43.7% biologic-naive.",
    dose_range = "Induction 300 mg IV (>= 30 min infusion) Q4W x 3 (LUCENT 1); maintenance 200 mg SC Q4W (LUCENT 2); IV 300 mg Q4W rescue / extended induction.",
    regions = "Global (34 countries); 10.7% Japanese.",
    notes = "LUCENT 1 (NCT03518086) and LUCENT 2 (NCT03524092) phase III trials pooled; Chua 2023 Table 1. BLQ samples (about 2%) excluded. Fit in NONMEM 7.4.2."
  )

  ini({
    # Typical values at the reference subject (70.7 kg, albumin 45.0 g/L,
    # BMI 25 kg/m^2). Rates are published per hour; x 24 converts to per day.
    lka <- log(0.00706 * 24); label("SC absorption rate constant (1/day)") # Table 2 LUCENT 1+2: Ka 0.00706 1/h
    lcl <- log(0.022 * 24); label("Clearance (L/day)") # Table 2 LUCENT 1+2: CL 0.022 L/h
    lvc <- log(3.11); label("Central volume of distribution (L)") # Table 2 LUCENT 1+2: V1 3.11 L
    lvp <- log(1.69); label("Peripheral volume of distribution (L)") # Table 2 LUCENT 1+2: V2 1.69 L
    lq <- log(0.0087 * 24); label("Intercompartmental clearance (L/day)") # Table 2 LUCENT 1+2: Q 0.00870 L/h
    logitfdepot <- log(0.476 / (1 - 0.476)); label("Logit of SC bioavailability at BMI 25 kg/m^2 (unitless)") # Table 2 LUCENT 1+2: F1 47.6% (footnote e: estimated on the logit scale)

    e_wt_cl_q <- 0.534; label("Power exponent of baseline weight on CL and Q (unitless)") # Table 2 LUCENT 1+2 and footnote c
    e_wt_vc_vp <- 0.525; label("Power exponent of baseline weight on Vc and Vp (unitless)") # Table 2 LUCENT 1+2 and footnote d
    e_alb_cl <- -0.981; label("Power exponent of albumin on CL (unitless)") # Table 2 LUCENT 1+2 and footnote c
    e_bmi_fdepot <- -0.0474; label("Linear effect of baseline BMI on logit SC bioavailability (per kg/m^2)") # Table 2 LUCENT 1+2

    # IIV: Table 2 footnote b, CV% = 100 * sqrt(exp(OMEGA) - 1), so
    # OMEGA = log(1 + CV^2). Covariance CL-V1 printed directly as 0.0241.
    etalcl + etalvc ~ c(
      0.086729,
      0.0241, 0.032236
    ) # Table 2 LUCENT 1+2: IIV CL 30.1%, V1 18.1%, covariance CL-V1 0.0241
    etalvp ~ 0.301325 # Table 2 LUCENT 1+2: IIV V2 59.3%
    etalogitfdepot ~ 0.413449 # Table 2 LUCENT 1+2: IIV F1 64.3%, footnote e: %IIV = 100 * sqrt(OMEGA), additive on the logit scale

    propSd <- 0.196; label("Proportional residual error (fraction)") # Table 2 LUCENT 1+2: 19.6%
    addSd <- 0.171; label("Additive residual error (ug/mL)") # Table 2 LUCENT 1+2: 171.0 ng/mL
  })
  model({
    wt_cl_q <- (WT / 70.7)^e_wt_cl_q
    wt_vc_vp <- (WT / 70.7)^e_wt_vc_vp
    alb_cl <- (ALB / 45.0)^e_alb_cl

    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * wt_cl_q * alb_cl
    vc <- exp(lvc + etalvc) * wt_vc_vp
    vp <- exp(lvp + etalvp) * wt_vc_vp
    q <- exp(lq) * wt_cl_q

    fdepot <- expit(logitfdepot + e_bmi_fdepot * (BMI - 25) + etalogitfdepot)

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
