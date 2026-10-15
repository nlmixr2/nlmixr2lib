Wang_2022_mycophenolic_acid <- function() {
  description <- paste0(
    "Two-compartment population PK model for mycophenolic acid (MPA) after ",
    "oral mycophenolate mofetil (MMF) dispersible tablets in adult Chinese ",
    "renal transplant recipients on tacrolimus and corticosteroids, sampled ",
    "either in the early post-transplant stage (days 4-8) or in the stable ",
    "state (5.5-10 years after transplantation) (Wang 2022). First-order ",
    "absorption after a lag time and first-order elimination from the ",
    "central compartment; no enterohepatic circulation term. The ",
    "post-transplant stage is the only retained covariate and acts ",
    "exponentially on CL/F and V/F: the stable state lowers CL/F from 23.36 ",
    "to about 10.3 L/h and V/F from 78.07 to about 6.2 L. Doses are mg of ",
    "MMF (no molecular-weight conversion; parameters are apparent with ",
    "respect to the MMF dose) and Cc is MPA in mg/L. IIV is exponential on ",
    "ka, V/F, V2/F, CL/F, Q/F and lag time; residual error is proportional ",
    "(37%)."
  )
  reference <- paste(
    "Wang P, Xie H, Zhang Q, Tian X, Feng Y, Qin Z, Yang J, Shang W, Feng G,",
    "Zhang X. (2022). Population Pharmacokinetics of Mycophenolic Acid in",
    "Renal Transplant Patients: A Comparison of the Early and Stable",
    "Posttransplant Stages. Front Pharmacol 13:859351.",
    "doi:10.3389/fphar.2022.859351",
    sep = " "
  )
  vignette <- "Wang_2022_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Doses are mg of MMF and concentrations mg/L of MPA (Methods; the NCA
  # CL/F of Table 2 equals MMF dose / AUC), so every state holds an amount
  # in MMF-dose-equivalent mg and all volumes and clearances are apparent
  # with respect to the MMF dose.
  compartmentData <- list(
    depot = list(
      analyte = "mycophenolate mofetil (MMF)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "mycophenolic acid (MMF-dose-equivalent mg)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "mycophenolic acid (MMF-dose-equivalent mg)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    POSTTX_STABLE = list(
      description = "Post-transplant stage indicator: 1 = stable state (sampled at least 5 years after renal transplantation), 0 = early post-transplant stage (sampled within the first week after transplantation).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (early post-transplant stage, days 4-8)",
      notes = "Methods: 'posttransplant stages (the early stage as 0 and the stable stage as 1)'. Exponential effect on CL/F and V/F, CL/F = 23.36 x exp(-0.82 x POSTTX_STABLE), V/F = 78.07 x exp(-2.54 x POSTTX_STABLE) (Table 3 dCLdStage, dVdStage). The cohort contains only the two extremes (days 4-8 and 5.5-10 years; Table 1 post-transplant time 4.88 +/- 1.01 and 2499.94 +/- 467.26 days), so the model says nothing about intermediate post-transplant times.",
      source_name = "Stage"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search; not retained (Results; Supplementary Figure S2)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained."
    ),
    WT = list(description = "Body weight.", units = "kg", type = "continuous", notes = "Screened; not retained."),
    ALT = list(
      description = "Alanine aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase.",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    TBILI = list(
      description = "Total bilirubin.",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    WBC = list(
      description = "White blood cell count.",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    HGB = list(description = "Haemoglobin.", units = "g/L", type = "continuous", notes = "Screened; not retained."),
    ALB = list(description = "Serum albumin.", units = "g/L", type = "continuous", notes = "Screened; not retained."),
    CREAT = list(
      description = "Serum creatinine.",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance (as reported, mL/min; not BSA-normalised).",
      units = "mL/min",
      type = "continuous",
      notes = "Screened; not retained. The Discussion notes a possible CrCL-CL correlation that was too small to include once the post-transplant stage was in the model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 99L,
    n_studies = 1L,
    n_observations = "1079 plasma MPA concentrations: 561 from 51 early-stage patients and 518 from 48 stable-state patients (Results). One 12-h steady-state profile per patient: pre-dose plus mainly 0.5, 1, 1.5, 2, 3, 4, 6, 8 and 12 h.",
    age_range = "Early stage 33.39 +/- 8.10 years; stable state 42.29 +/- 9.25 years (mean +/- SD; Table 1). Adults over 18 years.",
    weight_range = "Early stage 64.02 +/- 10.45 kg; stable state 65.31 +/- 10.54 kg (mean +/- SD; Table 1).",
    sex_female_pct = 18.2,
    race_ethnicity = c(Asian = 100),
    disease_state = "Adult renal transplant recipients, sampled either in the early post-transplant stage (days 4-8; post-transplant time 4.88 +/- 1.01 days) or in the stable state (at least 5 years; 2499.94 +/- 467.26 days). Patients with liver dysfunction, combined organ transplantation, rejection or gastrointestinal disease were excluded.",
    dose_range = "Oral MMF dispersible tablets twice daily under fasting conditions: 1.0-3.0 g/day early (1.69 +/- 0.44 g/day) and 0.5-1.5 g/day stable (1.09 +/- 0.29 g/day) (Methods; Table 1).",
    co_medication = "Triple therapy with tacrolimus and corticosteroids (methylprednisolone or prednisone) in all patients; early-stage patients also received pantoprazole (Table 1; Results).",
    renal_function = "Creatinine clearance 62.55 +/- 20.75 mL/min early and 74.51 +/- 19.29 mL/min stable (Table 1).",
    assay = "UPLC with diode-array detection, linear 0.04-40.0 mg/L.",
    regions = "China (First Affiliated Hospital of Zhengzhou University).",
    notes = "Male 43/51 early and 38/48 stable (81/99 overall). Estimation by FOCE-ELS in Phoenix NLME."
  )

  ini({
    lka <- log(1.36); label("Absorption rate constant ka (1/h)") # Table 3 'tvka (1/h)' 1.36
    lvc <- log(78.07); label("Apparent central volume V/F in the early post-transplant stage (L)") # Table 3 'tvV/F (L)' 78.07
    lvp <- log(554.52); label("Apparent peripheral volume V2/F (L)") # Table 3 'tvV2/F (L)' 554.52
    lcl <- log(23.36); label("Apparent clearance CL/F in the early post-transplant stage (L/h)") # Table 3 'tvCL/F (L/h)' 23.36
    lq <- log(29.53); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'tvQ/F (L/h)' 29.53
    ltlag <- log(0.23); label("Absorption lag time (h)") # Table 3 'tvTlag (h)' 0.23

    e_posttx_stable_vc <- -2.54; label("Exponential effect of the stable post-transplant state on V/F (unitless)") # Table 3 'dVdStage' -2.54
    e_posttx_stable_cl <- -0.82; label("Exponential effect of the stable post-transplant state on CL/F (unitless)") # Table 3 'dCLdStage' -0.82

    # IIV: Table 3 reports the omega^2 variances of exponential random
    # effects directly. No covariances are reported, so the matrix is
    # diagonal.
    etalvc ~ 1.03 # Table 3 'omega2 V' 1.03
    etalcl ~ 0.20 # Table 3 'omega2 CL' 0.20
    etalka ~ 0.34 # Table 3 'omega2 Ka' 0.34
    etalvp ~ 1.72 # Table 3 'omega2 V2' 1.72
    etalq ~ 0.98 # Table 3 'omega2 Q' 0.98
    etaltlag ~ 0.74 # Table 3 'omega2 Tlag' 0.74

    propSd <- 0.37; label("Proportional residual error (fraction)") # Table 3 'stdev0' 0.37; Methods proportional form Cobs = Cpred x (1 + eps)
  })

  model({
    # Phoenix NLME exponential categorical covariate form; the stable-state
    # CL/F reproduces the Results value 10.25 L/h (23.36 x exp(-0.82) = 10.29).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + e_posttx_stable_vc * POSTTX_STABLE + etalvc)
    vp <- exp(lvp + etalvp)
    cl <- exp(lcl + e_posttx_stable_cl * POSTTX_STABLE + etalcl)
    q <- exp(lq + etalq)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
