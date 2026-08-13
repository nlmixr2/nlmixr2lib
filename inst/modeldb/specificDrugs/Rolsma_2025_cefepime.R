Rolsma_2025_cefepime <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for",
    "intravenously infused cefepime in children and adults with cystic",
    "fibrosis (Rolsma 2025; 82 participants / 96 enrollments, 368 plasma",
    "concentrations, ages 3 to 54 years). Total clearance is the sum of a",
    "renal arm proportional to lean-body-weight-based creatinine clearance",
    "(2.84 L/h at CLCR,LBW = 77 mL/min) and a constant non-renal arm",
    "(0.928 L/h), giving 3.77 L/h at the cohort reference. The central",
    "volume of distribution scales linearly with fat-free mass (9.95 L at",
    "33 kg; the source reports fat-free mass as 'lean body weight' by the",
    "Janmahasatian formula). Correlated interindividual variability on renal",
    "clearance and volume (correlation 0.67) with a combined proportional",
    "plus additive residual error. CFTR modulator use, CFTR genotype and",
    "CF-related complications were screened and had no substantial effect.",
    sep = " "
  )
  reference <- paste(
    "Rolsma SL, Sokolow A, Patel PC, Sokolow K, Jimenez-Truque N, Fissell WH,",
    "Ryan V, Kirkpatrick CM, Nation RL, Gu K, Teresi M, Fishbane N, Kontos M,",
    "An G, Winokur P, Landersdorfer CB, Creech CB (2025). Population",
    "Pharmacokinetic Modeling of Cefepime, Meropenem, and",
    "Piperacillin-Tazobactam in Patients With Cystic Fibrosis. The Journal of",
    "Infectious Diseases 231(2):e364-e374. doi:10.1093/infdis/jiae451.",
    "Final parameter estimates from Supplemental Table 6A.",
    sep = " "
  )
  vignette <- "Rolsma_2025_betalactams_cysticfibrosis"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated creatinine clearance computed on lean body weight (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source name CLCR,LBW. Supplemental Table 6 abbreviation list:",
        "'CL_CR,LBW, creatinine clearance calculated according to the",
        "Cockcroft-Gault equation using lean body weight (subjects > 12 years",
        "of age) or the Schwartz bedside equation (subjects up to 12 years)'.",
        "The value is a raw mL/min clearance and is NOT BSA-normalized,",
        "matching the raw-mL/min variant of the canonical CRCL column (see",
        "the CRCL entry in inst/references/covariate-columns.md and the",
        "Georges_2009_ceftazidime.R / Delattre_2010_amikacin.R / ",
        "Chen_2023_nemonoxacin.R precedents). Note the two-formula",
        "construction: Cockcroft-Gault above 12 years of age, Schwartz",
        "bedside at or below 12 years, with lean body weight substituted for",
        "total body weight in the Cockcroft-Gault form. The effect on renal",
        "clearance is a through-origin linear normalization to a 77 mL/min",
        "reference, read from the Supplemental Table 6A unit string 'L/h per",
        "77 mL/min CLCR,LBW'; the same table's footnote a ('Sum of CL_R and",
        "CL_NR, for a CL_CR,LBW of 77 mL/min' = 3.77 L/h) confirms the",
        "reference by reproducing 2.84 + 0.928 = 3.768. No separate exponent",
        "parameter is reported, so the relationship is linear rather than a",
        "power model. Renal function was the strongest covariate found for",
        "cefepime (Results 'PopPK Modeling': 'Renal function, expressed as",
        "CLCR adjusted for LBW, best described the elimination clearance of",
        "cefepime and meropenem')."
      ),
      source_name        = "CLCR,LBW"
    ),
    FFM = list(
      description        = "Fat-free mass by the Janmahasatian formula (reported by the source as lean body weight)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source name LBW. The source labels this covariate 'lean body",
        "weight' but computes it with the Janmahasatian equations",
        "(Supplemental Methods: 'lean body weight based on the equations",
        "developed by Janmahasatian et al.'; Supplemental Table 3 footnote a:",
        "'Calculated according to method outlined in: Janmahasatian S,",
        "Duffull SB, Ash S et al. Quantification of lean bodyweight. Clin",
        "Pharmacokinet 2005; 44: 1051-65'). The canonical register",
        "discriminates FFM from LBM by the estimating formula, not by the",
        "paper's label: the FFM entry is defined as the Janmahasatian",
        "quantity, whereas LBM covers the Boer / Hume / James forms.",
        "Accordingly this is registered as FFM with source_name LBW; no value",
        "transformation is required (same quantity, same kg units). The",
        "effect on the central volume is a through-origin linear",
        "normalization to a 33 kg reference, read from the Supplemental Table",
        "6A unit string 'L per 33 kg LBW'; the reference matches the cohort",
        "median lean body weight of 33.45 kg in Supplemental Table 3. No",
        "exponent parameter is reported, so the relationship is linear",
        "(exponent 1) rather than a power model. Cohort lean body weight",
        "34.74 +/- 13.44 kg (median 33.45, range 11.1 to 69.5) per",
        "Supplemental Table 3."
      ),
      source_name        = "LBW"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate size descriptor but not retained in the",
        "final cefepime model, where lean body weight described the volume of",
        "distribution instead (Results 'PopPK Modeling': 'LBW was identified",
        "to affect the volume of distribution of cefepime'). Total body",
        "weight WAS retained in the companion meropenem model",
        "(Rolsma_2025_meropenem.R). Cohort 46.71 +/- 18.68 kg (median 46.86,",
        "range 13.6 to 102.5) per Supplemental Table 3."
      )
    ),
    CFTRMOD = list(
      description = "CFTR modulator use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and explicitly rejected. Results 'PopPK Modeling': 'Other",
        "covariates, such as sample acquisition site, CFTR mutation, use of",
        "CFTR modulators (including highly effective CFTR modulators",
        "ivacaftor and elexacaftor/tezacaftor/ivacaftor), or complications of",
        "CF (eg, diabetes) did not have substantial impact on overall PK.'",
        "The Conclusions restate this as a clinical finding: 'CFTR modulator",
        "therapy neither improved nor impaired the PK/PD of these",
        "antimicrobial agents; therefore, no additional modifications are",
        "needed for patients receiving these highly effective agents.' 35% of",
        "cefepime enrollments reported any modulator use (Supplemental Table",
        "2). No point estimate is reported, so no effect can be encoded."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods 'PopPK Model Development' covariate list) but not",
        "retained; the body-size and renal-function covariates absorbed the",
        "age trend across the 3 to 54 year range (Results 'PopPK Modeling':",
        "'incorporation of renal function and body size effectively described",
        "the PK of cefepime across a wide age range'). Age nevertheless",
        "stratifies the published probability-of-target-attainment analysis",
        "(3-11, 12-17 and > 17 years), which is a simulation stratification",
        "rather than a model covariate. Cohort 17.9 +/- 11.8 years (median",
        "15.0, range 3 to 54) per Supplemental Table 3."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained for cefepime (Results 'Patient",
        "Demographics': 'body surface area and BMI did not differ among",
        "antibiotic groups'). BMI does appear as a band-defining covariate in",
        "the companion meropenem model, where the source reports a separate",
        "non-renal clearance for BMI 45-50 kg/m^2 that it instructs must not",
        "be used for prediction. Cohort 19.71 +/- 4.17 kg/m^2 (median 19.15,",
        "range 13.0 to 35.7) per Supplemental Table 3."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 1.38 +/-",
        "0.37 m^2 (median 1.40, range 0.6 to 2.3) per Supplemental Table 3."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 150.44",
        "+/- 21.73 cm (median 155.70, range 95.0 to 185.0) per Supplemental",
        "Table 3."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. 54% female among",
        "cefepime enrollments (52 of 96) per Supplemental Table 2. Sex",
        "nevertheless enters the model indirectly, as an input to both the",
        "Janmahasatian fat-free-mass equation and the Cockcroft-Gault",
        "creatinine-clearance equation."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 82L,
    n_studies      = 1L,
    n_enrollments  = 96L,
    age_range      = "3 to 54 years",
    age_median     = "15.0 years (mean 17.9 +/- 11.8)",
    weight_range   = "13.6 to 102.5 kg total body weight",
    weight_median  = "46.86 kg total body weight (mean 46.71 +/- 18.68); lean body weight median 33.45 kg (mean 34.74 +/- 13.44)",
    sex_female_pct = 100 * 52 / 96,
    race_ethnicity = c(White = 97, `Black or African American` = 2, `Multi-Racial` = 1, `Hispanic or Latino (ethnicity)` = 3),
    disease_state  = paste(
      "Cystic fibrosis, admitted for a pulmonary exacerbation or for",
      "microbial eradication therapy. 95% of the cefepime group carried at",
      "least one copy of the DF508 mutation. 35% reported any CFTR modulator",
      "use (ivacaftor 16%, ivacaftor/tezacaftor 19%, ivacaftor/lumacaftor 8%,",
      "elexacaftor/tezacaftor/ivacaftor 5%). Staphylococcus aureus and",
      "Pseudomonas aeruginosa were the most common sputum isolates.",
      "Participants with a history of solid-organ or hematologic",
      "transplantation were excluded."
    ),
    renal_function = paste(
      "Estimated creatinine clearance on lean body weight (Cockcroft-Gault",
      "above 12 years of age, Schwartz bedside at or below 12 years). The",
      "structural reference is 77 mL/min per the Supplemental Table 6A unit",
      "string and footnote a. The cohort distribution of CLCR,LBW is not",
      "tabulated in the paper or supplement."
    ),
    dose_range     = paste(
      "Intravenous infusion, standard of care as ordered by the treating",
      "team. All adult doses were 2000 mg (n = 626). Most doses in",
      "participants under 17 years were 45 to 51 mg/kg (n = 1046); 1 dose was",
      "60 mg/kg and the remaining 496 doses were 30 to 44 mg/kg. Infusion",
      "durations in hours were 0.083 (24 doses), 0.167 (12 doses), 0.5 (1788",
      "doses), 3 (344 doses) and 6 (1 dose)."
    ),
    regions        = "United States (Vanderbilt University Medical Center, Nashville TN; University of Iowa Hospital, Iowa City IA).",
    notes          = paste(
      "Opportunistic sampling during hospitalization, January 2018 to March",
      "2020. 368 of the 667 total plasma samples in the study were cefepime.",
      "Assay LC-MS/MS, linear 0.5 to 150 mg/L for cefepime. NONMEM 7.4.3 with",
      "FOCE+I; outliers removed per the FDA popPK guidance (February 2022) at",
      "a base-model weighted residual greater than 5. One- and",
      "two-compartment models were evaluated and the one-compartment model",
      "with linear elimination was selected. Drug administration was modeled",
      "as a zero-order (infusion) process. Model evaluation by",
      "prediction-corrected VPC (Supplemental Figure 2). Baseline demographics",
      "from Supplemental Tables 2 and 3; final parameter estimates from",
      "Supplemental Table 6A. Protein binding of 20% was assumed for cefepime",
      "when converting simulated total concentrations to the unbound",
      "concentrations used in the target-attainment analysis; that factor is",
      "a post-processing step and is not part of this model."
    )
  )

  ini({
    # Structural parameters -- Rolsma 2025 Supplemental Table 6A, 'Population
    # estimate' column. Total clearance is additive across a renal and a
    # non-renal arm; the table reports both arms plus their sum.
    #
    # Supplemental Table 6A verbatim (Parameter / Unit / Population estimate):
    #   CL_R  / L/h per 77 mL/min CLCR,LBW / 2.84   (RSE 10%,  IIV 22%)
    #   CL_NR / L/h                        / 0.928  (RSE 28%)
    #   CL total / L/h                     / 3.77 a
    #   V     / L per 33 kg LBW            / 9.95   (RSE 5.1%, IIV 22%)
    #   CVCP  / -                          / 33%
    #   SDCP  / mg/L                       / 0.241
    #   footnote a: 'Sum of CL_R and CL_NR, for a CL_CR,LBW of 77 mL/min'
    #
    # The unit strings carry the covariate model: 'per 77 mL/min CLCR,LBW' and
    # 'per 33 kg LBW' denote through-origin linear normalization to those
    # reference values. Footnote a is a mass-balance check on that reading:
    # 2.84 + 0.928 = 3.768, which rounds to the published 3.77 L/h.
    lcl_renal  <- log(2.84);  label("Renal clearance at CLCR,LBW = 77 mL/min (L/h)")        # Suppl Table 6A row CL_R  = 2.84 L/h (RSE 10%)
    lcl_nonren <- log(0.928); label("Non-renal clearance (L/h)")                            # Suppl Table 6A row CL_NR = 0.928 L/h (RSE 28%)
    lvc        <- log(9.95);  label("Central volume of distribution at FFM = 33 kg (L)")    # Suppl Table 6A row V     = 9.95 L (RSE 5.1%)

    # Interindividual variability. The Supplemental Methods state the IIV model
    # explicitly ('The inter-individual variability (IIV) of the PK parameters
    # for all antibiotics in this study was estimated using an exponential
    # model'), so the tabulated percentages are CV% on the log-normal scale and
    # convert as omega^2 = log(CV^2 + 1).
    #
    # The 'Covariance' cell of 67% is the CORRELATION between the two etas, not
    # a covariance: read as a covariance it would imply a correlation of
    # 0.67 / 0.047265 = 14.2, which is impossible for a correlation bounded by
    # 1. The companion piperacillin and tazobactam tables (6C, 6D) label the
    # analogous cell 'Correlation (CL, V)', confirming the intended quantity.
    # cov = 0.67 * sqrt(0.047265 * 0.047265) = 0.031668.
    etalcl_renal + etalvc ~ c(0.047265,
                              0.031668, 0.047265)  # Suppl Table 6A: IIV CL_R 22% CV -> log(0.22^2+1) = 0.047265 (shrinkage 20%); IIV V 22% CV -> 0.047265 (shrinkage 44%); correlation 0.67

    # Residual error. Supplemental Methods: 'To describe the residual
    # variability (RV) of each antibiotic, a combined additive and proportional
    # error model was utilized.' CVCP is the proportional component and SDCP
    # the additive component, per the Supplemental Table 6 abbreviation list.
    propSd <- 0.33;  label("Proportional residual error (fraction)")  # Suppl Table 6A row CVCP = 33%    (RSE 12%)
    addSd  <- 0.241; label("Additive residual error (mg/L)")          # Suppl Table 6A row SDCP = 0.241 mg/L (RSE 212%)
  })

  model({
    # 1. Individual PK parameters. Both covariate relationships are
    # through-origin linear normalizations to the cohort reference values
    # printed in the Supplemental Table 6A unit strings (no exponent parameter
    # is reported for either, so neither is a power model).
    #
    # Renal clearance arm, proportional to lean-body-weight-based creatinine
    # clearance. The IIV eta sits on this arm alone, matching the placement of
    # the 22% IIV on the CL_R row of Supplemental Table 6A; the non-renal arm
    # carries no reported IIV.
    cl_renal  <- exp(lcl_renal + etalcl_renal) * (CRCL / 77)
    # Non-renal clearance arm, constant across the cohort (no covariate, no IIV).
    cl_nonren <- exp(lcl_nonren)
    cl        <- cl_renal + cl_nonren

    # Central volume, proportional to fat-free mass (source label: lean body
    # weight, Janmahasatian formula).
    vc <- exp(lvc + etalvc) * (FFM / 33)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment disposition with linear elimination. Cefepime is given
    # by intravenous infusion, so dose records enter the central compartment
    # directly as a zero-order input (Methods: 'Since all antibiotics were
    # administered by intravenous infusions at a constant rate, drug
    # administration was modeled as a zero-order process'); there is no depot.
    d/dt(central) <- -kel * central

    # 4. Observation. Total (not unbound) plasma cefepime in mg/L, matching the
    # LC-MS/MS assay the model was fitted to. Combined proportional plus
    # additive residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
