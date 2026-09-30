WerumeusBuning_2021_ceftazidime <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for ceftazidime in",
    "critically ill adults with a proven or suspected Pseudomonas aeruginosa",
    "infection (Werumeus Buning 2021; n = 96 ICU patients, 368 concentrations,",
    "mostly continuous infusion). Clearance is piecewise on continuous",
    "veno-venous hemofiltration (CVVH): patients on CVVH have a single fixed",
    "clearance with no between-subject variability, while patients off CVVH",
    "have a clearance that scales as a power of the CKD-EPI eGFR (centred on",
    "73 mL/min/1.73 m^2) and is multiplied by 1.57 for a hematologic",
    "malignancy and by 1.99 for trauma or head injury. Between-subject",
    "variability is exponential on the off-CVVH clearance and on the volume;",
    "residual error is additive on log-transformed concentrations."
  )
  reference <- paste(
    "Werumeus Buning A, Hodiamont CJ, Lechner NM, Schokkin M, Elbers PWG,",
    "Juffermans NP, Mathot RAA, de Jong MD, van Hest RM (2021).",
    "Population Pharmacokinetics and Probability of Target Attainment of",
    "Different Dosing Regimens of Ceftazidime in Critically Ill Patients with",
    "a Proven or Suspected Pseudomonas aeruginosa Infection.",
    "Antibiotics (Basel) 10(6):612.",
    "doi:10.3390/antibiotics10060612. PMCID PMC8224000.",
    sep = " "
  )
  vignette <- "WerumeusBuning_2021_ceftazidime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Section 4.5: "total serum concentrations were measured" (protein binding
  # approximately 10%) by LC-MS/MS, LLQ 0.5 mg/L.
  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate from the CKD-EPI 2009",
        "creatinine equation (Section 4.4), BSA-normalized."
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on the off-CVVH clearance only, centred on 73, the",
        "cohort median (Appendix C control stream: (CKD/73)**THETA(5), with",
        "THETA(5) = 0.772; Table 2 footnote a). Table 1 gives a baseline",
        "median of 73, range 6-153, with the unit printed as 'mL/min/m 2';",
        "the CKD-EPI equation is normalized to 1.73 m^2 by construction, and",
        "Figure 3 labels the 10th / 50th / 90th percentiles 33 / 73 / 122.",
        "Time-varying in the source dataset (collected 'during treatment';",
        "four missing values were carried forward). The term is switched",
        "off when RRT_CRRT_STATUS = 1."
      ),
      source_name = "CKD"
    ),
    RRT_CRRT_STATUS = list(
      description = paste(
        "1 = patient receiving continuous veno-venous hemofiltration (CVVH);",
        "0 = no CVVH."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects between the two clearance branches of the Appendix C",
        "control stream (IF(CVVH.EQ.0) ... ELSE CL = THETA(4)). On CVVH the",
        "clearance is the fixed 2.9 L/h with no IIV and no eGFR or",
        "comorbidity effect. 20/96 patients (21%) were on CVVH (Table 1).",
        "The dataset carried CVVH as a per-record column, so the flag can",
        "be supplied time-varying."
      ),
      source_name = "CVVH"
    ),
    DIS_HEME_MALIG = list(
      description = "1 = comorbidity hematologic malignancy; 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Multiplies the off-CVVH clearance by 1.57 (THETA(6); COMO category",
        "1 in the Appendix C control stream). 14/96 patients (15%). The",
        "source recorded a single comorbidity category per patient, so",
        "DIS_HEME_MALIG and DIS_TRAUMA are mutually exclusive in the",
        "fitted data; the reference group pools oncologic (solid)",
        "malignancy and 'other'."
      ),
      source_name = "COMO == 1"
    ),
    DIS_TRAUMA = list(
      description = paste(
        "1 = comorbidity acute trauma or head injury; 0 = otherwise."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Multiplies the off-CVVH clearance by 1.99 (THETA(7)). Pools COMO",
        "categories 3 (trauma) and 4 (brain injury) of the Appendix C",
        "control stream; the authors merged them 'due to having the same",
        "underlying mechanism for increasing drug clearance, being the",
        "hyperdynamic state with glomerular hyperfiltration' (Section 2.2).",
        "27/96 patients (28%)."
      ),
      source_name = "COMO == 3 or COMO == 4"
    )
  )

  # Appendix D lists the covariates tested for association with CL and V.
  # None below was retained and no coefficient is reported for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units = "years",
      type = "continuous",
      notes = "Screened (Appendix D); Table 1 median 59, range 20-84. Not retained.",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Appendix D); 38/96 female (Table 1). Not retained.",
      source_name = "sex"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      notes = "Screened (Appendix D); Table 1 median 79, range 44-237. Not retained.",
      source_name = "weight"
    ),
    BMI = list(
      description = "Body mass index.",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Appendix D); Table 1 median 25, range 16-66. Not retained.",
      source_name = "BMI"
    ),
    ALB = list(
      description = "Serum albumin.",
      units = "g/L",
      type = "continuous",
      notes = "Screened (Appendix D). Not retained; no summary is tabulated.",
      source_name = "serum albumin"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 96L,
    n_studies = 1L,
    n_observations = 368L,
    age_range = "20-84 years",
    age_median = "59 years",
    weight_range = "44-237 kg",
    weight_median = "79 kg",
    sex_female_pct = 40,
    disease_state = paste(
      "Critically ill adult ICU patients treated with intravenous",
      "ceftazidime for a proven or suspected clinically relevant",
      "Pseudomonas aeruginosa infection; cystic fibrosis excluded. Median",
      "SOFA 10 (range 4-16, n = 64); 30-day mortality 39%. Primary",
      "infection site pneumonia 39%, meningitis 23%, bloodstream 18%,",
      "abdominal 14%. Comorbidity hematologic malignancy 15%, oncologic",
      "malignancy 13%, trauma or head injury 28%, other 45%. Vasopressors",
      "69%, mechanical ventilation 77%, CVVH 21% (Table 1)."
    ),
    dose_range = paste(
      "Clinical regimens over 2013-2018: intermittent 1 g or 2 g three",
      "times daily, or 3-6 g per 24 h by continuous infusion with or",
      "without a loading dose (a bolus given over several minutes just",
      "before the infusion started). 83% of patients were on continuous",
      "infusion; 65 of 80 continuous-infusion patients received a loading",
      "dose (Table 1)."
    ),
    renal_function = paste(
      "CKD-EPI eGFR median 73 mL/min/1.73 m^2, range 6-153; serum creatinine",
      "median 0.98 mg/dL, range 0.19-7.49 (Table 1)."
    ),
    regions = "The Netherlands (single centre; ICU of Amsterdam University Medical Centre, location AMC).",
    notes = paste(
      "394 samples were drawn from 96 patients (routine TDM plus arterial",
      "blood-gas waste material); 28 (7.1%) were below the 0.5 mg/L LLQ and",
      "were handled with the M5 method, leaving 368 concentrations for the",
      "fit. NONMEM 7.1.2, FOCE with interaction, on log-transformed data.",
      "Final-model condition number 19.81; shrinkage 29% on both etas;",
      "98.2% of 1000 bootstrap runs successful."
    )
  )

  ini({
    # Appendix C control stream ($THETA, $OMEGA) and Table 2 'Final Model'.
    lcl <- log(3.42); label("Clearance off CVVH at eGFR 73, comorbidity 'other' (L/h)") # THETA(2); Table 2 CL non CVVH = 3.42 L/h (RSE 9%)
    lcl_cvvh <- log(2.9); label("Clearance on CVVH (L/h)") # THETA(4); Table 2 CL CVVH = 2.9 L/h (RSE 11%)
    lvc <- log(46.8); label("Volume of distribution (L)") # THETA(3); Table 2 V = 46.8 L (RSE 12%)

    e_crcl_cl <- 0.772; label("Power exponent of CKD-EPI eGFR/73 on off-CVVH CL (unitless)") # THETA(5); Table 2 CKD-EPI = 0.772 (RSE 11%)
    e_heme_malig_cl <- 1.57; label("Multiplicative factor on off-CVVH CL for hematologic malignancy (unitless)") # THETA(6); Table 2 = 1.57 (RSE 17%)
    e_trauma_cl <- 1.99; label("Multiplicative factor on off-CVVH CL for trauma or head injury (unitless)") # THETA(7); Table 2 = 1.99 (RSE 13%)

    # $OMEGA variances; sqrt(exp(omega^2) - 1) reproduces the Table 2 CV%.
    etalcl ~ 0.122 # $OMEGA 1 'IIV CL NON CVVH' = 0.122; Table 2 36.0% CV (applies off CVVH only)
    etalvc ~ 0.721 # $OMEGA 2 'IIV V' = 0.721; Table 2 102.8% CV

    # $ERROR: Y = LOG(F) + THETA(1)*EPS(1) with $SIGMA 1 FIX, i.e. an
    # additive SD of 0.281 on the log scale. Table 2 calls it 'Proportional
    # error'.
    expSd <- 0.281; label("Residual error SD on log-concentration scale") # THETA(1); Table 2 = 0.281 (RSE 12%)
  })
  model({
    # Off-CVVH clearance (Appendix C: THETA(2)*(CKD/73)**THETA(5)*
    # THETA(6)**FLAT1*THETA(7)**FLAT2*EXP(ETA(1))).
    cl_off <- exp(lcl + etalcl) * (CRCL / 73)^e_crcl_cl *
      e_heme_malig_cl^DIS_HEME_MALIG * e_trauma_cl^DIS_TRAUMA
    # On-CVVH clearance (Appendix C: ELSE CL = THETA(4)); no eta.
    cl_on <- exp(lcl_cvvh)
    cl <- cl_off * (1 - RRT_CRRT_STATUS) + cl_on * RRT_CRRT_STATUS

    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
