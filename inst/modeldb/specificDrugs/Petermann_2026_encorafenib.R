Petermann_2026_encorafenib <- function() {
  description <- "Two-compartment population PK model for oral encorafenib at steady state in real-world adults with BRAF V600E/K-mutant unresectable or metastatic melanoma treated with encorafenib plus binimetinib: first-order absorption with ka held constant at 0.51 1/h, linear elimination with no autoinduction term (all samples were drawn at least 25 days after treatment start or a dose change), log-normal IIV on CL/F, Vc/F and Q/F, a proportional residual error, and no retained covariates."
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
  # the studied covariates (tested on Vc, Q and CL for encorafenib ...) were
  # significantly associated with the IIV'). Documentation only.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F, Vc/F and Q/F; not retained. Cohort median 74.30 kg (IQR 65.00-87.00), Table 1."
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 25.01 kg/m^2 (IQR 22.29-27.74), Table 1."
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
      notes = "Tested; not retained. 21 (32%) female, Table 1."
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 39.70 g/L (IQR 36.40-43.00), Table 1."
    ),
    CRP = list(
      description = "C-reactive protein.",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 4.00 mg/L (IQR 2.00-7.45), Table 1; missing in up to 12% of records."
    ),
    LDH = list(
      description = "Lactate dehydrogenase.",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested; not retained. Cohort median 226.00 U/L (IQR 189.00-271.00), Table 1."
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status of 1 or worse.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = ECOG 0",
      notes = "ECOG performance status was tested; the paper does not state how the categories were grouped. ECOG 0/1/2/3 = 39/23/1/1 patients, Table 1. Not retained."
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic dysfunction by the National Cancer Institute Organ Dysfunction Working Group classification.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = normal hepatic function",
      notes = "Tested; not retained. Normal 43 (65%), mild dysfunction 16 (24%), Table 1."
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description = "Concomitant moderate CYP3A4 inhibitor.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no moderate CYP3A4 inhibitor",
      notes = "Tested; not retained. Only 2 patients (6 samples) took amiodarone and none took a strong CYP3A4 inhibitor (Results, 'Studied population'). Hospitalised versus ambulatory status was also tested and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "encorafenib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "encorafenib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "encorafenib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 66L,
    n_observations = 310L,
    n_studies = 1L,
    age_median = "62 years (IQR 52-71)",
    weight_median = "74.30 kg (IQR 65.00-87.00)",
    sex_female_pct = 32,
    disease_state = "Unresectable or metastatic melanoma with a BRAF V600E or V600K mutation, treated with encorafenib plus binimetinib in routine care.",
    dose_range = "Encorafenib 200-600 mg once daily (median 450 mg); 81% of patients on 450 mg QD and 12% on 300 mg QD.",
    regions = "France (Cochin University Hospital, Paris, and Ambroise Pare University Hospital, Boulogne-Billancourt), 2018-2024.",
    hepatic_function = "Normal 65%, mild dysfunction 24% (NCI-ODWG); no moderate or severe dysfunction.",
    notes = "Retrospective therapeutic-drug-monitoring data; all samples at steady state (at least 25 days after initiation or a dose change); median 3 (range 1-15) samples per patient, one per occasion; median time after dose 15.00 h (0.33-29.00). Baseline demographics in Petermann 2026 Table 1."
  )

  ini({
    lka <- fixed(log(0.51)); label("First-order absorption rate constant (1/h)") # Results 'Population pharmacokinetic model': ka held constant at 0.51 h-1 for encorafenib, Tmax 1.5 h (Table 2 rounds to 0.5)
    lcl <- log(29); label("Apparent clearance CL/F (L/h)") # Table 2, encorafenib CL = 29 L/h (RSE 5.0%)
    lvc <- log(39); label("Apparent central volume of distribution Vc/F (L)") # Table 2, encorafenib Vc = 39 L (RSE 16%)
    lq <- log(1.7); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2, encorafenib Q = 1.7 L/h (RSE 23.8%)
    lvp <- log(26); label("Apparent peripheral volume of distribution Vp/F (L)") # Table 2, encorafenib Vp = 26 L (RSE 12%)

    etalcl ~ 0.03189 # Table 2 omega CL = 18 CV%; Table 2 footnote CV = sqrt(exp(omega^2) - 1), so omega^2 = log(1 + 0.18^2)
    etalvc ~ 0.40819 # Table 2 omega Vc = 71 CV%; omega^2 = log(1 + 0.71^2)
    etalq ~ 0.61324 # Table 2 omega Q = 92 CV%; omega^2 = log(1 + 0.92^2)

    propSd <- 0.44675; label("Proportional residual error (fraction)") # Table 2 sigma prop = 47 CV%; Table 2 footnote applies CV = sqrt(exp(sigma^2) - 1), so sigma = sqrt(log(1 + 0.47^2))
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
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
