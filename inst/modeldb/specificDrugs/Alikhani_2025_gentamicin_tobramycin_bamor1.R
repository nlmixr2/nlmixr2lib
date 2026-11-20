Alikhani_2025_gentamicin_tobramycin_bamor1 <- function() {
  description <- "One-compartment IV-infusion population PK model for gentamicin and tobramycin (pooled) in hospitalised adults, with CL driven by the radiomic biological-age index BAMor_1 (CT psoas muscle area x L3 vertebral bone density / fascia area) divided by serum creatinine, and V by drug (Alikhani 2025 model 19.1)."
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
      notes = "Measured at the L4 vertebral level (Alikhani 2025 Results and Supplementary Table 1). The source reports mm^2 (cohort mean 1953.7 mm^2 = 19.537 cm^2, Table 1); supply cm^2. It enters BAMor_1 only as a ratio to CT_FASCIA_AREA, so the area units cancel.",
      source_name = "total_psoas_expmarea"
    ),
    CT_VERTEBRAL_BMD = list(
      description = "CT trabecular bone density of the L3 vertebral body",
      units = "HU",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean CT attenuation inside a central trabecular bone core sample at L3 (Supplementary Table 2 'bmdhuvbaligned'); cohort mean 182.23 HU (Table 1). Not DXA areal BMD.",
      source_name = "bmdhuvbaligned"
    ),
    CT_FASCIA_AREA = list(
      description = "Area enclosed by the abdominal fascia on CT at L3",
      units = "cm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at L3. The source reports mm^2 (cohort mean 50625 mm^2 = 506.25 cm^2, Table 1); supply cm^2. Enters BAMor_1 only as the denominator of a ratio of two areas, so the units cancel.",
      source_name = "fasciaarea"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort median 0.75 mg/dL (IQR 0.57, range 0.17-5.5; Table 1). Divides BAMor_1 inside the CL power term.",
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
    # Supplementary Table 4, column '19.1 BAMor_1-based'. Monolix log-normal
    # parameters with covariates on the log scale.
    lcl <- log(5.8); label("Typical CL at BAMor_1/Scr/10 = 1 (L/h)") # Supp Table 4 (19.1): Cl_pop = 5.8 L/h, SE 0.26, RSE 4.48
    lvc <- log(49.09); label("Typical V for gentamicin (L)") # Supp Table 4 (19.1): V_pop = 49.09 L, SE 3.27, RSE 6.66

    e_bamor1_creat_cl <- 0.15; label("Power exponent on BAMor_1/Scr/10 for CL (unitless)") # Supp Table 4 (19.1): Beta_Cl_(BAMor/creat) = 0.15, SE 0.047, RSE 30.5
    e_trt_tobramycin_vc <- 0.24; label("Log-scale effect of tobramycin (vs gentamicin) on V (unitless)") # Supp Table 4 (19.1): Beta_V_drug_tobramycin = 0.24, SE 0.083, RSE 34.5

    # Monolix reports omega as the SD of the log-normal random effect;
    # nlmixr2 takes the variance, omega^2.
    etalcl ~ 0.2401 # Supp Table 4 (19.1): Omega_Cl = 0.49 (SD), variance 0.49^2
    etalvc ~ 0.1369 # Supp Table 4 (19.1): Omega_V = 0.37 (SD), variance 0.37^2

    # Monolix combined2 error: y = f + sqrt(a^2 + (b*f)^2) * e.
    addSd <- 0.049; label("Additive residual error SD (mg/L)") # Supp Table 4 (19.1): a = 0.049, SE 0.023, RSE 46.7
    propSd <- 0.4; label("Proportional residual error SD (fraction)") # Supp Table 4 (19.1): b = 0.4, SE 0.013, RSE 3.38
  })
  model({
    # BAMor_1 = psoas muscle area x BMD / fascia area (Table 2). The two areas
    # share units, so the index is in HU whether they are supplied in cm^2 or
    # mm^2.
    bamor1 <- CT_PSOAS_MUSCLE_AREA * CT_VERTEBRAL_BMD / CT_FASCIA_AREA

    # Supplementary Table 3 row 19.1:
    #   log(Cl) = log(Cl_pop) + beta * log(BAMor_1 / creat / 10) + eta_Cl
    cl <- exp(lcl + etalcl) * (bamor1 / CREAT / 10)^e_bamor1_creat_cl
    vc <- exp(lvc + e_trt_tobramycin_vc * TRT_TOBRAMYCIN + etalvc)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
