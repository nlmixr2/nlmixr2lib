Alikhani_2025_gentamicin_tobramycin_bamor3 <- function() {
  description <- "One-compartment IV-infusion population PK model for gentamicin and tobramycin (pooled) in hospitalised adults, with CL driven by the radiomic biological-age index BAMor_3 (CT dorsal muscle group area x L3 vertebral bone density / visceral fat area) divided by serum creatinine and by CT subcutaneous fat density, and V by drug (Alikhani 2025 model 29.2, the paper's best model)."
  reference <- paste(
    "Alikhani R, Horbal SR, Pai MP. A comparison between radiomic biological",
    "age and chronological age in estimating kidney function. Sci Rep.",
    "2025;15:13384. doi:10.1038/s41598-025-98297-1.",
    sep = " "
  )
  vignette <- "Alikhani_2025_gentamicin_tobramycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CT_DORSAL_MUSCLE_AREA = list(
      description = "Dorsal muscle group cross-sectional area on abdominal CT at L2",
      units = "cm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at the L2 vertebral level (Alikhani 2025 Results and Supplementary Table 1). The source reports mm^2 (cohort mean 4030 mm^2 = 40.30 cm^2, Table 1); supply cm^2. It enters BAMor_3 only as a ratio to VISCERAL_ABDOMINAL_FAT, so the area units cancel.",
      source_name = "dmgexpmarea"
    ),
    CT_VERTEBRAL_BMD = list(
      description = "CT trabecular bone density of the L3 vertebral body",
      units = "HU",
      type = "continuous",
      reference_category = NULL,
      notes = "Mean CT attenuation inside a central trabecular bone core sample at L3 (Supplementary Table 2 'bmdhuvbaligned'); cohort mean 182.23 HU (Table 1). Not DXA areal BMD.",
      source_name = "bmdhuvbaligned"
    ),
    VISCERAL_ABDOMINAL_FAT = list(
      description = "Visceral fat area on abdominal CT at L3",
      units = "cm^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at the L3 vertebral level. The source reports mm^2 (cohort mean 12790 mm^2 = 127.90 cm^2, SD 10708.68 mm^2, range 40-49471 mm^2, Table 1); supply cm^2. Enters BAMor_3 only as the denominator of a ratio of two areas, so the units cancel.",
      source_name = "visceralfatarea"
    ),
    CT_SUBCUTANEOUS_FAT_DENSITY = list(
      description = "Median CT attenuation of subcutaneous fat at L3",
      units = "HU",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured at L3; cohort mean -89.03 HU (SD 12.76, range -112.08 to -59.83; Table 1, where the unit is misprinted as mm2 -- Supplementary Table 2 defines it in HU). Enters CL uncentred as exp(e * CT_SUBCUTANEOUS_FAT_DENSITY / 90) (Equation 1 and Supplementary Table 3 'subQhu/90').",
      source_name = "subcutfathu_mean"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort median 0.75 mg/dL (IQR 0.57, range 0.17-5.5; Table 1). Divides BAMor_3 inside the CL power term.",
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
    # Supplementary Table 4, column '29.2 BAMor_3-based'. Monolix log-normal
    # parameters with covariates on the log scale.
    lcl <- log(2.62); label("Typical CL at BAMor_3/(Scr x 300) = 1 and subcutaneous fat density 0 HU (L/h)") # Supp Table 4 (29.2): Cl_pop = 2.62 L/h, SE 1.02, RSE 38.9
    lvc <- log(45.89); label("Typical V for gentamicin (L)") # Supp Table 4 (29.2): V_pop = 45.89 L, SE 3.23, RSE 7.03

    e_bamor3_creat_cl <- 0.11; label("Power exponent on BAMor_3/(Scr x 300) for CL (unitless)") # Supp Table 4 (29.2): Beta_Cl_(BAMor/creat) = 0.11, RSE 29.4
    e_ct_subcutaneous_fat_density_cl <- -0.93; label("Exponential coefficient on subcutaneous fat density/90 for CL (unitless)") # Supp Table 4 (29.2): Beta_Cl_(subcutfathu) = -0.93, SE 0.38, RSE 40.5
    e_trt_tobramycin_vc <- 0.28; label("Log-scale effect of tobramycin (vs gentamicin) on V (unitless)") # Supp Table 4 (29.2): Beta_V_drug_tobramycin = 0.28, SE 0.085, RSE 30.2

    # Monolix reports omega as the SD of the log-normal random effect;
    # nlmixr2 takes the variance, omega^2.
    etalcl ~ 0.2401 # Supp Table 4 (29.2): Omega_Cl = 0.49 (SD), variance 0.49^2
    etalvc ~ 0.1296 # Supp Table 4 (29.2): Omega_V = 0.36 (SD), variance 0.36^2

    # Monolix combined2 error: y = f + sqrt(a^2 + (b*f)^2) * e.
    addSd <- 0.13; label("Additive residual error SD (mg/L)") # Supp Table 4 (29.2): a = 0.13, SE 0.038, RSE 30.5
    propSd <- 0.39; label("Proportional residual error SD (fraction)") # Supp Table 4 (29.2): b = 0.39, SE 0.014, RSE 3.67
  })
  model({
    # BAMor_3 = dorsal muscle group area x BMD / visceral fat area (Table 2).
    # The two areas share units, so the index is in HU whether they are
    # supplied in cm^2 or mm^2.
    bamor3 <- CT_DORSAL_MUSCLE_AREA * CT_VERTEBRAL_BMD / VISCERAL_ABDOMINAL_FAT

    # Equation 1: log(Cl) = log(Cl_pop) + beta * log(BAMor_3 / (Scr * 300))
    #   + beta_subQ * SubcutfatHU / 90 + eta_Cl
    cl <- exp(lcl + etalcl) *
      (bamor3 / (CREAT * 300))^e_bamor3_creat_cl *
      exp(e_ct_subcutaneous_fat_density_cl * CT_SUBCUTANEOUS_FAT_DENSITY / 90)
    vc <- exp(lvc + e_trt_tobramycin_vc * TRT_TOBRAMYCIN + etalvc)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined2()
  })
}
