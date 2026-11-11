Alikhani_2025_gentamicin_tobramycin_ca <- function() {
  description <- "One-compartment IV-infusion population PK model for gentamicin and tobramycin (pooled) in hospitalised adults, with CL driven by the inverse of chronological age x serum creatinine and V by drug (Alikhani 2025 model 6.1, the chronological-age comparator to the radiomic biological-age models)."
  reference <- paste(
    "Alikhani R, Horbal SR, Pai MP. A comparison between radiomic biological",
    "age and chronological age in estimating kidney function. Sci Rep.",
    "2025;15:13384. doi:10.1038/s41598-025-98297-1.",
    sep = " "
  )
  vignette <- "Alikhani_2025_gentamicin_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort median 58 years (IQR 22, range 21-93; Table 1). Enters CL through the index 1/(AGE x CREAT).",
      source_name = "Age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort median 0.75 mg/dL (IQR 0.57, range 0.17-5.5; Table 1). Enters CL through the index 1/(AGE x CREAT).",
      source_name = "creat"
    ),
    TRT_TOBRAMYCIN = list(
      description = "Aminoglycoside received: 1 = tobramycin, 0 = gentamicin",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (gentamicin)",
      notes = "Monolix categorical covariate 'drug' with gentamicin as reference (Supplementary Table 4 'Beta_V_drug_tobramycin'). 89 of 156 patients (57.05%) received tobramycin (Table 1).",
      source_name = "drug"
    )
  )

  compartmentData <- list(
    central = list(analyte = "gentamicin or tobramycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 156L,
    n_studies = 1L,
    age_range = "21-93 years",
    age_median = "58 years (IQR 22)",
    weight_range = "43.8-139.3 kg",
    weight_median = "81.0 kg (mean; SD 21.06)",
    height_median = "170.5 cm (mean; SD 10.79)",
    sex_female_pct = 51.92,
    race_ethnicity = "Not reported",
    disease_state = "Hospitalised adults treated with gentamicin (n = 67) or tobramycin (n = 89) who had an abdominal CT scan available for analytic morphomics",
    dose_range = "Not reported; each patient had at least three aminoglycoside doses with measured plasma concentrations",
    regions = "United States (University of Michigan Health System, single centre)",
    renal_function = "Serum creatinine median 0.75 mg/dL (IQR 0.57, range 0.17-5.5)",
    notes = "Retrospective cohort; 156 of 363 adults with CT images met the inclusion criteria (amikacin-treated patients excluded). Demographics from Alikhani 2025 Table 1. Monolix 2024R1, SAEM."
  )

  ini({
    # Supplementary Table 4, column '6.1 CA-based'. Monolix log-normal
    # parameters with covariates on the log scale.
    lcl <- log(5.82); label("Typical CL at 1/(AGE x Scr)/0.03 = 1 (L/h)") # Supp Table 4 (6.1): Cl_pop = 5.82 L/h, SE 0.26, RSE 4.55
    lvc <- log(49.2); label("Typical V for gentamicin (L)") # Supp Table 4 (6.1): V_pop = 49.2 L, SE 3.46, RSE 7.04

    e_age_creat_cl <- 0.19; label("Power exponent on 1/(AGE x Scr)/0.03 for CL (unitless)") # Supp Table 4 (6.1): Beta_Cl_(CA-1/creat) = 0.19, SE 0.062, RSE 32.7
    e_trt_tobramycin_vc <- 0.25; label("Log-scale effect of tobramycin (vs gentamicin) on V (unitless)") # Supp Table 4 (6.1): Beta_V_drug_tobramycin = 0.25, SE 0.086, RSE 33.8

    # Monolix reports omega as the SD of the log-normal random effect;
    # nlmixr2 takes the variance, omega^2.
    etalcl ~ 0.2401 # Supp Table 4 (6.1): Omega_Cl = 0.49 (SD), variance 0.49^2
    etalvc ~ 0.1444 # Supp Table 4 (6.1): Omega_V = 0.38 (SD), variance 0.38^2

    # Monolix combined2 error: y = f + sqrt(a^2 + (b*f)^2) * e.
    addSd <- 0.069; label("Additive residual error SD (mg/L)") # Supp Table 4 (6.1): a = 0.069, SE 0.037, RSE 53.5
    propSd <- 0.4; label("Proportional residual error SD (fraction)") # Supp Table 4 (6.1): b = 0.4, SE 0.014, RSE 3.51
  })
  model({
    # Supplementary Table 3 row 6.1:
    #   log(Cl) = log(Cl_pop) + beta * log((1 / (Age * creat)) / 0.03) + eta_Cl
    cl <- exp(lcl + etalcl) * ((1 / (AGE * CREAT)) / 0.03)^e_age_creat_cl
    vc <- exp(lvc + e_trt_tobramycin_vc * TRT_TOBRAMYCIN + etalvc)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
