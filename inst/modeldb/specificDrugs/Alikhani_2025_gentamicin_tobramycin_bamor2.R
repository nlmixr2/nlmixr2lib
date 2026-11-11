Alikhani_2025_gentamicin_tobramycin_bamor2 <- function() {
  description <- "One-compartment IV-infusion population PK model for gentamicin and tobramycin (pooled) in hospitalised adults, with CL driven by the radiomic biological-age index BAMor_2 (CT psoas muscle area x psoas muscle density x L3 vertebral slab height) divided by serum creatinine, and V by drug (Alikhani 2025 model 25.1)."
  reference <- paste(
    "Alikhani R, Horbal SR, Pai MP. A comparison between radiomic biological",
    "age and chronological age in estimating kidney function. Sci Rep.",
    "2025;15:13384. doi:10.1038/s41598-025-98297-1.",
    sep = " "
  )
  vignette <- "Alikhani_2025_gentamicin_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CT_PSOAS_MUSCLE_AREA = list(
      description = "Total psoas muscle cross-sectional area on abdominal CT at L4",
      units = "cm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at the L4 vertebral level (Alikhani 2025 Results and Supplementary Table 1). The source reports mm^2 (cohort mean 1953.7 mm^2 = 19.537 cm^2, Table 1); supply cm^2. model() multiplies by 100 to recover mm^2, the units in which the published centring constant 4200000 was set.",
      source_name = "total_psoas_expmarea"
    ),
    CT_PSOAS_MUSCLE_DENSITY = list(
      description = "Mean CT attenuation of the psoas muscle at L4",
      units = "HU",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at L4; cohort mean 46.28 HU (SD 11.73, range 0.24-72.64; Table 1).",
      source_name = "avg_psoas_expmmeanhu"
    ),
    CT_VERTEBRAL_SLAB_HEIGHT = list(
      description = "Vertebral body slab height on CT at L3",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at L3; cohort mean 32.78 mm (SD 3.83, range 20.46-42.13; Table 1).",
      source_name = "vbslabheight"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort median 0.75 mg/dL (IQR 0.57, range 0.17-5.5; Table 1). Divides BAMor_2 inside the CL power term.",
      source_name = "Scr"
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
    notes = "Retrospective cohort; 156 of 363 adults with CT images met the inclusion criteria (amikacin-treated patients excluded). Demographics and CT morphomic summaries from Alikhani 2025 Table 1; missing morphomic values were imputed with mice (random forest). Monolix 2024R1, SAEM."
  )

  ini({
    # Supplementary Table 4, column '25.1 BAMor_2-based'. Monolix log-normal
    # parameters with covariates on the log scale.
    lcl <- log(5.72); label("Typical CL at BAMor_2/Scr/4200000 = 1 (L/h)") # Supp Table 4 (25.1): Cl_pop = 5.72 L/h, SE 0.27, RSE 4.65
    lvc <- log(45.87); label("Typical V for gentamicin (L)") # Supp Table 4 (25.1): V_pop = 45.87 L, SE 3.39, RSE 7.39

    e_bamor2_creat_cl <- 0.12; label("Power exponent on BAMor_2/Scr/4200000 for CL (unitless)") # Supp Table 4 (25.1): Beta_Cl_(BAMor/creat) = 0.12, SE 0.046, RSE 38.2
    e_trt_tobramycin_vc <- 0.29; label("Log-scale effect of tobramycin (vs gentamicin) on V (unitless)") # Supp Table 4 (25.1): Beta_V_drug_tobramycin = 0.29, SE 0.087, RSE 29.5

    # Monolix reports omega as the SD of the log-normal random effect;
    # nlmixr2 takes the variance, omega^2.
    etalcl ~ 0.2601 # Supp Table 4 (25.1): Omega_Cl = 0.51 (SD), variance 0.51^2
    etalvc ~ 0.1369 # Supp Table 4 (25.1): Omega_V = 0.37 (SD), variance 0.37^2

    # Monolix combined2 error: y = f + sqrt(a^2 + (b*f)^2) * e.
    addSd <- 0.12; label("Additive residual error SD (mg/L)") # Supp Table 4 (25.1): a = 0.12, SE 0.043, RSE 34.9
    propSd <- 0.38; label("Proportional residual error SD (fraction)") # Supp Table 4 (25.1): b = 0.38, SE 0.014, RSE 3.59
  })
  model({
    # BAMor_2 = psoas muscle area x psoas muscle density x vertebral slab
    # height (Table 2), with the area in the source's mm^2 (cm^2 x 100).
    bamor2 <- (CT_PSOAS_MUSCLE_AREA * 100) * CT_PSOAS_MUSCLE_DENSITY * CT_VERTEBRAL_SLAB_HEIGHT

    # Supplementary Table 3 row 25.1:
    #   log(Cl) = log(Cl_pop) + beta * log(BAMor_2 / creat / 4200000) + eta_Cl
    cl <- exp(lcl + etalcl) * (bamor2 / CREAT / 4200000)^e_bamor2_creat_cl
    vc <- exp(lvc + e_trt_tobramycin_vc * TRT_TOBRAMYCIN + etalvc)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
