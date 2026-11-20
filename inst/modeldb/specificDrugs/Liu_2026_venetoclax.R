Liu_2026_venetoclax <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for oral venetoclax in ",
    "adult Chinese patients with hematological malignancies, mostly acute ",
    "myeloid leukemia, treated in routine practice (Liu 2026): first-order ",
    "absorption with ka fixed at 0.16 1/h, linear elimination, CL/F 7.33 L/h ",
    "and V/F 181 L, log-normal between-subject variability on CL/F only, and ",
    "a proportional residual error. No covariate was retained. Developed from ",
    "121 trough concentrations from 76 patients at a single Chongqing centre; ",
    "the companion exposure-response model is Liu_2026_venetoclax_orr."
  )
  reference <- paste0(
    "Liu Z, Liu X, Gong Q, Qin S, Zhu X, Kuan IHS, Mak WY, Xiang X, Jia C, ",
    "Wang Q, Cheng L, Li L. Real-World Pharmacokinetic and Exposure-Response ",
    "Characterization of Venetoclax in Chinese Patients with Hematological ",
    "Malignancies. Drug Des Devel Ther. 2026;20:563629. ",
    "doi:10.2147/DDDT.S563629."
  )
  vignette <- "Liu_2026_venetoclax"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "venetoclax", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Screened in the stepwise covariate search (Methods: demographics plus the
  # blood-routine, liver- and renal-function results of Table 1; forward
  # p < 0.05, backward p < 0.01) and none retained (Results: "there was no
  # statistically covariates included in the final model"). Documentation
  # only; not referenced in model(). Units are those of Table 1.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained (Liu 2026 Methods and Results). Table 1 median 53 (IQR 36-63)."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 60 (IQR 54-65)."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained. 30/76 female (Table 1: 46 male)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 22.13 (IQR 19.31-24.45)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained; no effect on AUC or Cmin in the Figure 3 / Figure S6 forest plots. Table 1 median 14.7."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained; shown in the Figure 3 / Figure S6 forest plots. Table 1 median 16."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained; shown in the Figure 3 / Figure S6 forest plots. Table 1 median 33.5."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 83."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 34.25; most patients were below the normal range (Figure S2)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained; shown in the Figure 3 / Figure S6 forest plots. Table 1 median 15.3."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 4.29."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 66.7. eGFR (median 99.74 mL/min/1.73 m^2) was also screened and appears in the forest plots."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 0.82."
    ),
    HGB = list(
      description = "Hemoglobin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 70."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 21."
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained; shown in the Figure 3 / Figure S6 forest plots. Table 1 median 31."
    ),
    NEUT = list(
      description = "Absolute neutrophil count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained. Table 1 median 0.18."
    ),
    LYMPH_ABS = list(
      description = "Absolute lymphocyte count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained; shown in the Figure 3 / Figure S6 forest plots. Table 1 median 0.44."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 76L,
    n_studies = 1L,
    n_concentrations = 121L,
    age_median = "53 years (IQR 36-63)",
    weight_median = "60 kg (IQR 54-65)",
    bmi_median = "22.13 kg/m^2 (IQR 19.31-24.45)",
    sex_female_pct = 39.5,
    race_ethnicity = c(Asian = 100),
    disease_state = paste0(
      "Adult hematological malignancies: acute myeloid leukemia 63, acute ",
      "lymphoblastic leukemia 3, acute promyelocytic leukemia 1, chronic ",
      "myelomonocytic leukemia 1, myelodysplastic syndrome 7, mixed-phenotype ",
      "acute leukemia 1 (Table 1). Marked cytopenias at the time of sampling ",
      "(median WBC 0.82, neutrophils 0.18 and platelets 31 x 10^9/L; median ",
      "albumin 34.25 g/L)."
    ),
    dose_range = paste0(
      "Oral venetoclax 100, 200, 300 or 400 mg once daily after a ramp-up, ",
      "taken with meals; in the 63 AML patients 7 received 100 mg, 13 ",
      "received 200 mg, 2 received 300 mg and 41 received 400 mg (Table S1). ",
      "Every patient also received one or two chemotherapy agents (mostly ",
      "azacitidine, with homoharringtonine, cytarabine, aclarubicin, ",
      "idarubicin, daunorubicin, vincristine, decitabine or cladribine)."
    ),
    regions = "Single centre, Southwest Hospital, Chongqing, China (July 2021 to March 2024)",
    conmed_exclusions = "Moderate and strong CYP3A inhibitors and rituximab were exclusion criteria.",
    notes = paste0(
      "Prospective, non-interventional real-world therapeutic drug monitoring ",
      "study. Every sample was a pre-dose trough, so ka could not be ",
      "estimated and was fixed at 0.16 1/h from earlier publications. ",
      "Plasma venetoclax by LC-MS/MS, LLOQ 20 ng/mL. Estimation in NONMEM ",
      "7.5 by FOCE-I. One- and two-compartment structures were tested and ",
      "the one-compartment model was kept. Before building this model the ",
      "authors externally evaluated six published venetoclax population PK ",
      "models (Jones 2016, Minocha 2018, Deng 2019, two Samineni 2022 ",
      "analyses and Gong 2023; Table S2) on these data and found that every ",
      "one under-predicted the ",
      "Chinese concentrations, with population-prediction median errors of ",
      "-68% to -94% (Table 2)."
    )
  )

  ini({
    # ----- Structural PK (Liu 2026 Table 3, final model column) -----
    # One-compartment with first-order absorption and linear elimination
    # (Results, 'VNX PK Characteristics in Chinese Patients').
    lka <- fixed(log(0.16)); label("Absorption rate constant (1/h)")  # Table 3, Ka = 0.16 /h FIX; fixed from earlier publications because every sample was a trough
    lcl <- log(7.33); label("Apparent clearance CL/F (L/h)")  # Table 3, CL/F = 7.33 L/h (RSE 10%; bootstrap median 7.26, 95% CI 6.04-8.77)
    lvc <- log(181); label("Apparent volume of distribution V/F (L)")  # Table 3, Vd/F = 181 L (RSE 30%; bootstrap median 183.04, 95% CI 102.01-348.74)

    # ----- Between-subject variability (Liu 2026 Table 3) -----
    # Exponential IIV on CL/F only. Table 3 prints IIV as 53% CV; converted to a
    # log-scale variance with omega^2 = log(1 + CV^2) = 0.2475.
    etalcl ~ log(1 + 0.53^2)  # Table 3, IIV CL/F = 53% CV (RSE 11%, shrinkage 6%)

    # ----- Residual error (Liu 2026 Table 3) -----
    propSd <- 0.32; label("Proportional residual error (fraction)")  # Table 3, sigma prop = 32% (RSE 10%, shrinkage 25%)
  })

  model({
    # Individual PK parameters (no covariates in the final model)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # CL and V are apparent (CL/F, V/F); bioavailability is absorbed into them,
    # so rxode2's default F = 1 applies.

    # Dose in mg and vc in L give mg/L; multiply by 1000 for ng/mL, the unit of
    # the LC-MS/MS assay and of the Cmin / Cave exposure metrics.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
