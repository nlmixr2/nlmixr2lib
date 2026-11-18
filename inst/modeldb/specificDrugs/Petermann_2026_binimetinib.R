Petermann_2026_binimetinib <- function() {
  description <- "Two-compartment population PK model for oral binimetinib at steady state in real-world adults with BRAF V600E/K-mutant unresectable or metastatic melanoma treated with encorafenib plus binimetinib: first-order absorption with ka held constant at 0.76 1/h, linear elimination with no autoinduction term (all samples were drawn at least 25 days after treatment start or a dose change), log-normal IIV on CL/F and Q/F, a proportional residual error, and no retained covariates."
  reference <- paste(
    "Petermann Y. J., Puszkiel A., Khani A., Kramkimel N., Funck-Brentano E.,",
    "Oules B., Bouchereau S., Vidal M., Decleves X., Csajka C., Blanchet B.,",
    "Guidi M. (2026).",
    "Population pharmacokinetics of encorafenib and binimetinib in real-world",
    "patients with BRAF V600E/K-mutant metastatic melanoma.",
    "Cancer Chemotherapy and Pharmacology 96(1):29.",
    "doi:10.1007/s00280-026-04876-y.",
    sep = " "
  )
  vignette <- "Petermann_2026_encorafenib_binimetinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  # Covariates screened in the stepwise covariate analysis (Methods, 'Base and
  # covariate models') and not retained in the final model (Results: 'None of
  # the studied covariates (tested ... on Q and CL for binimetinib) were
  # significantly associated with the IIV'). Documentation only.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and Q/F; not retained. Cohort median 74.00 kg (IQR 64.80-87.00), Table 1."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 24.78 kg/m^2 (IQR 22.32-27.72), Table 1."
    ),
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 62 years (IQR 52-71), Table 1."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Tested; not retained. 20 (29%) female, Table 1."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 39.70 g/L (IQR 36.62-42.42), Table 1."
    ),
    CRP = list(
      description = "C-reactive protein.",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 4.00 mg/L (IQR 2.00-9.55), Table 1; missing in up to 15% of records."
    ),
    LDH = list(
      description = "Lactate dehydrogenase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 226.00 U/L (IQR 187.00-284.00), Table 1."
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status of 1 or worse.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = ECOG 0",
      notes = "ECOG performance status was tested; the paper does not state how the categories were grouped. ECOG 0/1/2/3 = 36/21/1/2 patients, Table 1. Not retained."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic dysfunction by the National Cancer Institute Organ Dysfunction Working Group classification.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = normal hepatic function",
      notes = "Tested; not retained. Normal 40 (58%), mild dysfunction 16 (23%), Table 1; no moderate or severe dysfunction. Interacting comedications and hospitalised versus ambulatory status were also tested and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "binimetinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "binimetinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "binimetinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 69L,
    n_observations = 311L,
    n_studies = 1L,
    age_median = "62 years (IQR 52-71)",
    weight_median = "74.00 kg (IQR 64.80-87.00)",
    sex_female_pct = 29,
    disease_state = "Unresectable or metastatic melanoma with a BRAF V600E or V600K mutation, treated with encorafenib plus binimetinib in routine care.",
    dose_range = "Binimetinib 30-45 mg once or twice daily (median 45 mg); 81% of patients on 45 mg and 19% on 30 mg BID.",
    regions = "France (Cochin University Hospital, Paris, and Ambroise Pare University Hospital, Boulogne-Billancourt), 2018-2024.",
    hepatic_function = "Normal 58%, mild dysfunction 23% (NCI-ODWG); no moderate or severe dysfunction.",
    notes = "Retrospective therapeutic-drug-monitoring data; all samples at steady state (at least 25 days after initiation or a dose change); median 3 (range 1-15) samples per patient, one per occasion; median time after dose 12.00 h (0.25-20.00). Baseline demographics in Petermann 2026 Table 1."
  )

  ini({
    lka <- fixed(log(0.76)); label("First-order absorption rate constant (1/h)") # Results 'Population pharmacokinetic model': ka held constant at 0.76 h-1 for binimetinib, Tmax 1.1 h (Table 2 rounds to 0.7)
    lcl <- log(23); label("Apparent clearance CL/F (L/h)") # Table 2, binimetinib CL = 23 L/h (RSE 4.2%)
    lvc <- log(35); label("Apparent central volume of distribution Vc/F (L)") # Table 2, binimetinib Vc = 35 L (RSE 23%)
    lq <- log(35); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2, binimetinib Q = 35 L/h (RSE 50%)
    lvp <- log(106); label("Apparent peripheral volume of distribution Vp/F (L)") # Table 2, binimetinib Vp = 106 L (RSE 10%)

    etalcl ~ 0.05155 # Table 2 omega CL = 23 CV%; Table 2 footnote CV = sqrt(exp(omega^2) - 1), so omega^2 = log(1 + 0.23^2)
    etalq ~ 1.34088 # Table 2 omega Q = 168 CV%; omega^2 = log(1 + 1.68^2)

    propSd <- 0.43811; label("Proportional residual error (fraction)") # Table 2 sigma prop = 46 CV%; Table 2 footnote applies CV = sqrt(exp(sigma^2) - 1), so sigma = sqrt(log(1 + 0.46^2))
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volume in L give mg/L; x 1000 reports ng/mL as in the paper.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
