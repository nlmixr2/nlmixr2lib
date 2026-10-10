Courlet_2022_palbociclib <- function() {
  description <- "Final population PK-only model of oral palbociclib in women with advanced breast cancer followed in routine care (Courlet 2022, Supplementary Table S1). Two-compartment model with first-order absorption, an absorption lag time and first-order elimination. Apparent oral clearance takes a separately estimated value when palbociclib is taken under fasting conditions (or with a light meal) together with a proton-pump inhibitor (+56%); PPI use with a meal has no effect. IIV on ka, CL/F and Vc/F; proportional residual error. The jointly estimated PK/PD (neutropenia) model of the same paper is Courlet_2022_palbociclib_anc."
  reference <- paste(
    "Courlet P, Cardoso E, Bandiera C, Stravodimou A, Zurcher JP, Chtioui H, Locatelli I, Decosterd LA,",
    "Darnaud L, Blanchet B, Alexandre J, Wagner AD, Zaman K, Schneider MP, Guidi M, Csajka C. (2022).",
    "Population Pharmacokinetics of Palbociclib and Its Correlation with Clinical Efficacy and Safety in",
    "Patients with Advanced Breast Cancer. Pharmaceutics 14(7):1317. doi:10.3390/pharmaceutics14071317.",
    sep = " "
  )
  vignette <- "Courlet_2022_palbociclib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Unit note. Supplementary Table S1 reports CL/F and Q in L/h, Vc/F and Vp/F
  # in L, ka in 1/h and ALAG in h; doses are in mg (Section 3.1: 75-125 mg per
  # day), so central / vc is in mg/L. Concentrations are reported in ng/mL
  # (Section 2.2: validated range 0.5-500 ng/mL), hence the 1000 factor on Cc.

  compartmentData <- list(
    depot = list(analyte = "palbociclib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "palbociclib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "palbociclib", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor co-administered with the palbociclib dose",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no PPI co-administration)",
      notes = "Per-concentration (time-varying) indicator; 78 of 255 concentrations (31%) were measured under PPI co-administration (Table 1). PPI use alone had no effect on any PK parameter (Section 3.2.1, dOFV > -0.36); it acts only in combination with fasting, through the interaction CONMED_PPI * (1 - FED) that selects the separately estimated clearance lcl_ppifasted.",
      source_name = "PPI"
    ),
    FED = list(
      description = "Palbociclib taken with a meal (1) versus under fasting conditions or with a light meal (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (taken with a meal, as labelled)",
      notes = "Per-concentration indicator. The source flags 'administration under fasting conditions' (54 of 255 concentrations, 21%, Table 1), and Section 3.2.1 groups 'fasting conditions (or with a light meal)' together, so FED = 1 - fasting flag. Fasting alone had no retained effect; it acts only together with a PPI (CONMED_PPI * (1 - FED)).",
      source_name = "fasting conditions"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 65 [IQR 55-75] years.",
      source_name = "Age"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F and Vc/F (Section 2.3.3); not retained. Table 1 median 67 [IQR 61-80] kg.",
      source_name = "Body weight"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained.",
      source_name = "BMI"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 23 [19-28] U/L.",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 20 [15-27] U/L.",
      source_name = "ALT"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 5 [4-7] umol/L.",
      source_name = "BILT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 43 [41-45] g/L.",
      source_name = "Albumin"
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 61 [49-81] U/L.",
      source_name = "ALK"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL/F (Section 2.3.3); not retained. Table 1 median 71 [63-85] umol/L.",
      source_name = "CRT"
    ),
    CRCL = list(
      description = "BSA-normalized renal function (Cockcroft-Gault estimate per the source)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper's 'eGFR, estimated with the Cockroft and Gault formula' was tested on CL/F (Section 2.3.3); not retained. Table 1 reports it as 71 [54-98] in mL/min/1.73 m^2.",
      source_name = "eGFR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    age_range = "median 65 years (IQR 55-75)",
    weight_range = "median 67 kg (IQR 61-80)",
    sex_female_pct = 100,
    disease_state = "Women treated with palbociclib for advanced (metastatic) breast cancer in routine care (OpTAT study, NCT04484064); 70% of ANC records under concomitant fulvestrant.",
    dose_range = "Oral palbociclib 75-125 mg once daily (median 100 mg), 21 days on / 7 days off (28-day cycle).",
    regions = "Switzerland (Lausanne University Hospital); external validation also used 19 patients from Cochin Hospital, Paris, France.",
    observations = "255 plasma concentrations (2.0-159.0 ng/mL; median 5 per patient, range 1-16) from 5 min to 255 h after the last dose, including samples during OFF-treatment periods.",
    notes = "Baseline characteristics from Table 1. Estimated in NONMEM 7.4.3 with FOCEI (ADVAN4 TRANS4). One patient with severe adherence issues was excluded."
  )

  ini({
    # Supplementary Table S1, 'final palbociclib PK-only model'.
    lka <- log(0.73); label("ka: first-order absorption rate constant (1/h)") # Table S1: ka = 0.73 1/h (RSE 34%)
    ltlag <- log(1.9); label("ALAG: absorption lag time (h)") # Table S1: ALAG = 1.9 h (RSE 4%)
    lcl <- log(68); label("CL/F: apparent clearance, fed or without PPI (L/h)") # Table S1: CL = 68 L/h (RSE 5%)
    lcl_ppifasted <- log(106); label("CL/F under fasting conditions with a co-administered PPI (L/h)") # Table S1: CLIPP,no food = 106 L/h (RSE 9%)
    lvc <- log(2730); label("Vc/F: apparent central volume of distribution (L)") # Table S1: Vc = 2730 L (RSE 10%)
    lq <- log(5.1); label("Q/F: apparent intercompartmental clearance (L/h)") # Table S1: Q = 5.1 L/h (RSE 43%)
    lvp <- log(717); label("Vp/F: apparent peripheral volume of distribution (L)") # Table S1: Vp = 717 L (RSE 10%)

    # IIV modelled exponentially (Section 2.3.1); CV% converted with omega^2 = log(CV^2 + 1).
    etalka ~ 0.95072 # Table S1: omega ka = 126% CV -> log(1.26^2 + 1) = 0.95072
    etalcl ~ 0.086178 # Table S1: omega CL = 30% CV -> log(0.30^2 + 1) = 0.086178
    etalvc ~ 0.121864 # Table S1: omega Vc = 36% CV -> log(0.36^2 + 1) = 0.121864

    propSd <- 0.18; label("Proportional residual error (fraction)") # Table S1: proportional residual error = 18% (RSE 19%)
  })

  model({
    # Section 3.2.1: CL/F differs only when palbociclib is taken fasting (or
    # with a light meal) AND with a PPI; PPI use alone has no effect.
    ppi_fasted <- CONMED_PPI * (1 - FED)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp((1 - ppi_fasted) * lcl + ppi_fasted * lcl_ppifasted + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # central / vc is mg/L; 1 mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
