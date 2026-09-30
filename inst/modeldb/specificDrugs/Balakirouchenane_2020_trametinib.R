Balakirouchenane_2020_trametinib <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption, an",
    "absorption lag time and linear elimination for oral trametinib in a",
    "real-life cohort of adults with BRAF V600-mutated solid tumours, mostly",
    "metastatic melanoma, co-treated with dabrafenib (Balakirouchenane 2020).",
    "All disposition parameters are apparent oral values. No covariate was",
    "retained. Inter-individual variability on CL/F and Q/F; additive residual",
    "error.",
    sep = " "
  )
  reference <- paste(
    "Balakirouchenane D, Guegan S, Csajka C, Jouinot A, Heidelberger V,",
    "Puszkiel A, Zehou O, Khoudour N, Courlet P, Kramkimel N, Lheure C,",
    "Franck N, Huillard O, Arrondeau J, Vidal M, Goldwasser F, Maubec E,",
    "Dupin N, Aractingi S, Guidi M, Blanchet B.",
    "Population Pharmacokinetics/Pharmacodynamics of Dabrafenib Plus",
    "Trametinib in Patients with BRAF-Mutated Metastatic Melanoma.",
    "Cancers. 2020;12(4):931. doi:10.3390/cancers12040931.",
    "Parameter estimates from Table 3; population from Table 1.",
    sep = " "
  )
  vignette <- "Balakirouchenane_2020_dabrafenib_trametinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "trametinib",
      units = "mg (apparent, i.e. amount/F)",
      specimen = "administration site",
      verified = FALSE
    ),
    central = list(
      analyte = "trametinib",
      units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma",
      verified = FALSE
    ),
    peripheral1 = list(
      analyte = "trametinib",
      units = "mg (apparent, i.e. amount/F)",
      specimen = "plasma",
      verified = FALSE
    )
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested (Methods 4.4.2) but not retained; none of the tested covariates explained trametinib variability (dOFV < 1.426, p > 0.05; Results 2.2.2). Median 73.0 (53.0-166.0) kg (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested but not retained. Median 25.9 (18.3-44.6) kg/m^2 (Table 1)."
    ),
    FFM = list(
      description = "Fat-free mass (Janmahasatian)",
      units = "kg",
      type = "continuous",
      notes = "Tested but not retained. Median 52.7 (34.7-94.4) kg (Table 1)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Tested but not retained. Table 1 prints the trametinib-cohort BSA as '3.1 (2.5-3.0)', an evident typesetting error (the median lies outside its range)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested but not retained. Median 60.0 (20.0-90.0) years (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested but not retained. 27 of 60 patients (45%) female (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained. Median 30.0 (12.0-103.0) U/L (Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested but not retained. Median 28.0 (13.0-116.0) U/L (Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested but not retained. Median 6.3 (1.4-47.0) umol/L (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested but not retained. Median 43.0 (31.0-48.0) g/L (Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Tested but not retained. Median 17.2 (1.0-248.0) mg/L (Table 1)."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Tested explicitly (Results 2.2.2) but not retained. 52 of 60 patients (87%) took a PPI (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60L,
    n_studies = 1L,
    n_observations = "318 trametinib concentrations",
    age_range = "20-90 years (median 60.0)",
    weight_range = "53.0-166.0 kg (median 73.0)",
    sex_female_pct = 45,
    disease_state = paste(
      "Adults with BRAF V600-mutated metastatic solid tumours treated in routine",
      "care with dabrafenib plus trametinib: metastatic melanoma 52 (87%),",
      "other (anaplastic thyroid or non-small-cell lung carcinoma) 8 (13%)",
      "(Table 1).",
      sep = " "
    ),
    dose_range = "Trametinib 2 mg orally once daily (reduced starting doses such as 1 mg allowed at physician discretion), with dabrafenib 150 mg twice daily",
    regions = "France (three Assistance Publique-Hopitaux de Paris hospitals: Cochin, Henri Mondor, Avicenne)",
    notes = paste(
      "Observational multicentre 'real-life' cohort (July 2015 to June 2017);",
      "sparse routine samples at steady state. HPLC-MS/MS assay calibrated",
      "5-50 ng/mL; BLQ samples excluded. Fit in NONMEM 7.4.1 (FOCEI) with PsN",
      "4.8.0; 500-replicate bootstrap, pcVPC, and external validation on 46",
      "samples from 15 Lausanne patients (bias 2%).",
      sep = " "
    )
  )

  ini({
    # ---- Structural parameters (Table 3, final TRA model; apparent values) ----
    lka <- log(0.913); label("First-order absorption rate constant (1/h)") # Table 3 'ka (1/h) 0.913' (RSE 38.2%); bootstrap median 1.04
    ltlag <- log(0.709); label("Absorption lag time (h)") # Table 3 'Tlag (h) 0.709' (RSE 17.1%); bootstrap median 0.728
    lcl <- log(5.83); label("Apparent clearance CL/F (L/h)") # Table 3 'CL/F (L/h) 5.83' (RSE 4.6%); bootstrap median 5.82
    lvc <- log(61.9); label("Apparent central volume V2/F (L)") # Table 3 'V2/F (L) 61.9' (RSE 26.8%); bootstrap median 65.1
    lq <- log(64.9); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F (L/h) 64.9' (RSE 23.4%); bootstrap median 62.6
    lvp <- log(417.0); label("Apparent peripheral volume V3/F (L)") # Table 3 'V3/F (L) 417.0' (RSE 42.2%); bootstrap median 448.2

    # ---- Inter-individual variability (Table 3) ----
    # Reported as CV%; on the exponential-eta scale omega^2 = log(1 + CV^2).
    etalcl ~ 0.083988 # log(1 + 0.296^2); Table 3 'IIV CL/F (%) 29.6'
    etalq ~ 0.496648 # log(1 + 0.802^2); Table 3 'IIV Q (%) 80.2'

    # ---- Residual error (Table 3) ----
    addSd <- 4.14; label("Additive residual error (ng/mL)") # Table 3 'RUV (ng/mL) 4.14' (RSE 6.2%), additive error
  })

  model({
    mgPerLToNgPerMl <- 1000

    ka <- exp(lka)
    tlag <- exp(ltlag)
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

    lag(depot) <- tlag

    Cc <- central / vc * mgPerLToNgPerMl
    Cc ~ add(addSd)
  })
}
