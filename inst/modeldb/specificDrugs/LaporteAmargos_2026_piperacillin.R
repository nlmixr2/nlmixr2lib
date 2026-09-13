LaporteAmargos_2026_piperacillin <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous piperacillin",
    "(given as piperacillin-tazobactam) in adult hematological patients",
    "with febrile neutropenia enrolled in the BEATLE randomized trial",
    "(Laporte-Amargos 2026), with first-order elimination and time-varying",
    "Cockcroft-Gault creatinine clearance entering clearance as a power",
    "term centered on the cohort median of 99.3 mL/min. Between-subject",
    "variability is carried on clearance and volume of distribution,",
    "between-occasion variability on clearance across the three sampling",
    "occasions, with a combined additive-plus-proportional residual-error",
    "model. Tazobactam concentrations were not measured and are not",
    "described by this model."
  )
  reference <- paste(
    "Laporte-Amargos J, Ulldemolins M, Hernandez-Mitre MP, Roberts JA,",
    "Rigo-Bonnin R, Carmona-Torre F, Huguet M, Puerta-Alcalde P, Arnan M,",
    "del Pozo JL, Torrent A, Garcia-Vidal C, Sureda A, Bergas A,",
    "Sastre-Escola E, Carratala J, Gudiol C (2026).",
    "Population pharmacokinetics and optimized dosing of",
    "piperacillin-tazobactam in hematological patients with febrile",
    "neutropenia.",
    "Antimicrob Agents Chemother 70(1):e01253-25.",
    "doi:10.1128/aac.01253-25.",
    sep = " "
  )
  vignette <- "LaporteAmargos_2026_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "piperacillin", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min, NOT",
        "body-surface-area normalized. The canonical CRCL column carries",
        "mL/min/1.73 m^2 by default; this model uses the raw",
        "un-normalized variant, the same convention as",
        "Delattre_2010_amikacin.R and Chen_2023_nemonoxacin.R. Supplying a",
        "BSA-normalized value would silently rescale the renal-function",
        "term.",
        "TIME-VARYING: Laporte-Amargos 2026 Methods ('Basic model building,",
        "covariate analysis, and model diagnostics') states CrCL was",
        "estimated at each sampling time, and Results ('Population PK",
        "analysis') calls it 'time-varying CrCL'. It was the only covariate",
        "retained, reducing between-subject variability on CL by 18.6%; its",
        "effect on between-occasion variability was negligible because",
        "intra-subject CrCL varied by a median of only 6.6% between",
        "occasions.",
        "The paper also screened CKD-EPI estimated glomerular filtration",
        "rate (an alias of this same canonical column, in",
        "mL/min/1.73 m^2, cohort mean 105.8) but retained the",
        "Cockcroft-Gault estimate; supply the Cockcroft-Gault value.",
        "Cohort median 96.2 mL/min (Table 1, interquartile interval",
        "83.4-125.7); the model's reference value is the median 99.3 mL/min",
        "across sampling occasions quoted in Results. Monte Carlo dosing",
        "simulations were run at fixed CrCL of 60, 90, 120 and 150 mL/min",
        "(Table 3), the range the paper considers representative of this",
        "population. Patients with eGFR < 30 mL/min/1.73 m^2 were excluded,",
        "so the model does not describe moderate-to-severe renal",
        "impairment."
      ),
      source_name        = "CrCL"
    ),
    OCC = list(
      description        = "Integer-valued sampling-occasion indicator for the between-occasion variability on clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Laporte-Amargos 2026 Methods ('Antibiotic dosing, data collection,",
        "and blood sampling') states that plasma samples were collected",
        "during the first 5 days of treatment on THREE occasions, and",
        "Results reports 221 samples from 122 dosing occasions across 44",
        "patients (2.8 occasions per patient on average). Three occasions",
        "are therefore encoded. The supplement (Supplementary material on",
        "Methods) gives the variability model as",
        "theta_ik = theta_p * exp(eta_i) * exp(eta_ik), i.e. the",
        "between-occasion eta is additive with the between-subject eta on",
        "the log scale. Only one between-occasion magnitude was estimated",
        "(Table 2, BOV CL 16.4% CV), so occasions 2 and 3 are fixed equal to",
        "occasion 1. The index is decomposed inside model() into",
        "mutually-exclusive indicators oc1..oc3; a record with OCC outside",
        "1..3 carries no between-occasion variability. For a",
        "single-occasion simulation pass OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at inclusion",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the forward-inclusion / backward-elimination covariate search (Laporte-Amargos 2026 Methods; supplement 'Supplementary material on Methods') but not retained. Cohort mean 55.4 years (SD 10.2, Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. 21 of 44 patients (47.7%) were male (Table 1), i.e. 52.3% female."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort mean 165.6 cm (SD 10.8, Table 1). Height enters the model only indirectly, through the Cockcroft-Gault CrCL that was retained."
    ),
    WT = list(
      description = "Body weight at admission",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened but not retained; the model carries NO allometric size term. Cohort median 70 kg (interquartile interval 62.4-77.6, Table 1). Weight enters only indirectly through the Cockcroft-Gault CrCL."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort mean 3.4 g/dL (SD 0.4, Table 1). The paper's Limitations note that unbound concentrations were derived with an assumed 30% protein binding, considered reasonable because albumin was within the physiological range."
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort mean 5.6 g/dL (SD 0.5, Table 1)."
    ),
    TBILI = list(
      description = "Total plasma bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort median 0.62 mg/dL (interquartile interval 0.47-0.94, Table 1)."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score at the onset of febrile neutropenia",
      units       = "(points)",
      type        = "continuous",
      notes       = "Screened but not retained. Cohort mean 18.0 (SD 3.5, Table 1). Two further baseline severity scores were screened alongside it and are likewise absent from the final model: the Sepsis-Related Organ Failure Assessment score (cohort mean 5, SD 1.6) and the Multinational Association for Supportive Care in Cancer score (20 of 44 patients, 45.5%, in the high-risk stratum with a score below 21). Neither has a canonical covariate column in nlmixr2lib, so they are documented here rather than given placeholder entries."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 1L,
    age_range      = "Adults (>= 18 years by inclusion criterion); mean 55.4 years, SD 10.2",
    age_median     = "Not reported; mean 55.4 years",
    weight_range   = "62.4-77.6 kg (interquartile interval)",
    weight_median  = "70 kg",
    sex_female_pct = 52.3,
    race_ethnicity = "Not reported.",
    disease_state  = paste(
      "Febrile neutropenia (axillary temperature >= 38.0 C with < 500",
      "neutrophils/mm^3, or < 1,000 expected to drop within 24-48 h) in",
      "adults undergoing chemotherapy for acute leukemia or hematopoietic",
      "stem cell transplantation. Underlying malignancy: multiple myeloma",
      "31.8%, lymphoma 29.5%, acute myeloid leukemia or myelodysplastic",
      "syndrome 20.5%, acute lymphoblastic leukemia 6.8%, other 11.4%.",
      "Hematopoietic stem cell transplant was the reason for admission in",
      "90.9%. Only 3 patients (6.8%) were hypotensive at onset and none",
      "required vasoactive drugs; 1 patient (2.3%) was admitted to",
      "intensive care. Patients with eGFR < 30 mL/min/1.73 m^2 by CKD-EPI",
      "were EXCLUDED, so the model does not describe moderate-to-severe",
      "renal impairment. Renal function skewed high: 50% of patients had",
      "CrCL > 90 mL/min and 25% exceeded 120 mL/min, with mean 24-h",
      "diuresis above 2,000 mL, which the authors read as augmented renal",
      "clearance."
    ),
    dose_range     = paste(
      "Piperacillin-tazobactam 4 g / 0.5 g every 6 h for patients with",
      "eGFR > 40 mL/min/1.73 m^2 and every 8 h for eGFR 30-40",
      "mL/min/1.73 m^2. The first dose was always a 30 min infusion;",
      "thereafter 21 patients (47.7%) received 3 h extended infusions and",
      "23 (52.3%) continued with 30 min intermittent infusions."
    ),
    regions        = "Spain (four university hospitals).",
    notes          = paste(
      "Pharmacokinetic substudy of the BEATLE multicenter randomized",
      "controlled trial (extended versus short infusion of beta-lactams in",
      "hematological patients with febrile neutropenia), conducted November",
      "2019 to June 2022. 221 total plasma piperacillin concentrations from",
      "122 dosing occasions; 24 of the 44 patients also underwent intensive",
      "sampling over one dosing occasion. Concentrations were measured by",
      "UHPLC-MS/MS with a lower limit of quantification of 0.58 mg/L and a",
      "measuring interval of 0.58-175 mg/L. Baseline demographics are in",
      "Table 1. The model was fit in Monolix 2024R1 (SAEM); parameter",
      "estimates and a nonparametric bootstrap (n = 1,000) are in Table 2.",
      "TOTAL piperacillin was measured; the paper derives unbound",
      "concentrations for its PK/PD target attainment by assuming 30%",
      "protein binding, i.e. fu = 0.7. This model therefore predicts TOTAL",
      "plasma piperacillin."
    )
  )

  ini({
    # =========================================================================
    # Structural parameters (Laporte-Amargos 2026 Table 2, 'Fixed effects').
    # CL is the typical clearance at the reference CrCL of 99.3 mL/min,
    # because the covariate term in the Results equation is normalized to
    # that median.
    # =========================================================================
    lcl <- log(12.0)
    label("Total body clearance at the reference Cockcroft-Gault creatinine clearance of 99.3 mL/min (L/h)")
    # Table 2 row 'CL (L/h) = 12.0 (6.2%) [16.5%]'; bootstrap median 12.0
    # (2.5-97.5 percentiles 10.8-13.6). Also quoted in the Discussion:
    # 'our population estimate for piperacillin was 12 L/h at a median CrCL
    # of ~100 mL/min'.

    lvc <- log(29.8)
    label("Volume of distribution (L)")
    # Table 2 row 'Vd (L) = 29.8 (10.9%) [18%]'; bootstrap median 29.7
    # (23.4-37.9). Discussion: 'the estimated Vd of piperacillin was
    # approximately 30 L, nearly doubling the values reported in healthy
    # volunteers'. No allometric or other size term was retained.

    # =========================================================================
    # Covariate effect. Results, 'Population PK analysis': the printed
    # equation is CL_i (L/h) = 12 x (CrCL_i / 99.3)^0.64, with the text
    # 'CrCL = CrCL in mL/min estimated with the Cockcroft-Gault equation at
    # each sampling time, normalized to the median CrCL of our patient
    # population (99.3 mL/min)'. The supplement confirms the general form:
    # 'For continuous covariates inclusion, covariates were log-transformed
    # and normalised to their median values', which is exactly a power term
    # on the normalized covariate.
    # =========================================================================
    e_crcl_cl <- 0.64
    label("Power exponent of Cockcroft-Gault creatinine clearance on clearance (unitless)")
    # Table 2 row 'CrCL effect on CL = 0.64 (27.2%)'; bootstrap median 0.62
    # (0.23-0.91). Printed as the exponent in the Results equation.

    # =========================================================================
    # Between-subject variability (Laporte-Amargos 2026 Table 2, 'Random
    # effects'). The supplement states individual parameters are log-normal
    # and BSV / BOV are exponential:
    # theta_ik = theta_p * exp(eta_i) * exp(eta_ik). The Table 2 footnote
    # reads 'BSV, between-subject variability expressed as coefficient of
    # variation (CV %)', so the printed percentages are CV, not omega. For a
    # log-normal parameter CV = sqrt(exp(omega^2) - 1), hence
    # omega^2 = log(1 + CV^2), which is what nlmixr2lib stores.
    # =========================================================================
    etalcl ~ 0.1231411
    # Table 2 row 'BSV CL (CV %) = 36.2 (14.9%)'; bootstrap median 36.1
    # (21.2-50.9). omega^2 = log(1 + 0.362^2) = 0.1231411; sqrt = 0.3509

    etalvc ~ 0.2839218
    # Table 2 row 'BSV Vd (CV %) = 57.3 (15.2%)'; bootstrap median 56.9
    # (35.0-79.5). omega^2 = log(1 + 0.573^2) = 0.2839218; sqrt = 0.5328

    # =========================================================================
    # Between-occasion variability on CL. One magnitude was estimated and
    # three sampling occasions were collected per patient, so occasions 2
    # and 3 are fixed equal to occasion 1.
    # =========================================================================
    etaiov_cl_1 ~ 0.0265407
    # Table 2 row 'BOV CL (CV %) = 16.4 (23.5%)'; bootstrap median 14.1
    # (8.2-23.3). omega^2 = log(1 + 0.164^2) = 0.0265407; sqrt = 0.1629
    etaiov_cl_2 ~ fixed(0.0265407)
    etaiov_cl_3 ~ fixed(0.0265407)

    # =========================================================================
    # Residual variability (Laporte-Amargos 2026 Table 2, 'Residual
    # variability'). The Table 2 footnote defines 'a, constant component of
    # the residual error; b, proportional component of the residual error'.
    # The supplement states 'Additive, proportional or combined (additive +
    # proportional) error models were tested'; the combined model was
    # retained. Monolix parameterizes its combined1 error model as
    # SD = a + b * f and its combined2 model as SD = sqrt(a^2 + (b*f)^2);
    # the paper does not say which was selected, so combined1() is used here
    # -- it is Monolix's first-listed combined model and matches the
    # 'additive + proportional' wording of the supplement. combined1() is
    # required explicitly because rxode2 defaults to the quadrature
    # (combined2) form. See the vignette's Assumptions and deviations
    # section.
    # =========================================================================
    addSd <- 5.5
    label("Constant (additive) component of the residual error (mg/L)")
    # Table 2 row 'a (constant) (mg/L) = 5.5 (19.3%)'; bootstrap median 5.5
    # (0.15-10.5)

    propSd <- 0.20
    label("Proportional component of the residual error (fraction)")
    # Table 2 row 'b (proportional) = 0.20 (11.1%)'; bootstrap median 0.21
    # (0.10-0.26)
  })

  model({
    # --- 1. Occasion indicators and between-occasion variability on CL -----
    # Mutually-exclusive indicators over the three sampling occasions; a
    # record with OCC outside 1..3 zeroes every term and carries no BOV.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3

    # --- 2. Individual parameters -----------------------------------------
    # CL_i = 12 * (CrCL_i / 99.3)^0.64 (Results equation), with log-normal
    # between-subject and between-occasion variability entering additively
    # on the log scale per the supplement's
    # theta_ik = theta_p * exp(eta_i) * exp(eta_ik).
    cl <- exp(lcl + etalcl + iov_cl) * (CRCL / 99.3)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    # --- 3. Micro-constants -----------------------------------------------
    kel <- cl / vc

    # --- 4. ODE system ----------------------------------------------------
    # One compartment with first-order elimination (Results, 'Population PK
    # analysis'). Piperacillin-tazobactam is given intravenously, so doses
    # enter central directly; the source study used 30 min, 3 h and 4 h
    # infusions and the dosing simulations added 24 h continuous infusions.
    d/dt(central) <- -kel * central

    # --- 5. Observation and residual error --------------------------------
    # Dose in mg over volume in L gives mg/L, the unit the paper reports.
    # This is TOTAL plasma piperacillin; the paper's unbound targets assume
    # 30% protein binding, i.e. fCc = 0.7 * Cc.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
