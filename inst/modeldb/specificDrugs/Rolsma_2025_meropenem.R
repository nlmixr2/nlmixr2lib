Rolsma_2025_meropenem <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for",
    "intravenously infused meropenem in adolescents and adults with cystic",
    "fibrosis (Rolsma 2025; 42 participants / 50 enrollments, 192 plasma",
    "concentrations, ages 13 to 65 years). Total clearance is the sum of a",
    "renal arm proportional to lean-body-weight-based creatinine clearance",
    "(6.44 L/h at CLCR,LBW = 90 mL/min) and a non-renal arm scaling",
    "allometrically with total body weight (4.01 L/h at 61 kg, exponent",
    "0.75), giving 10.5 L/h at the cohort reference. The central volume is",
    "17.5 L with no covariate. Independent interindividual variability on",
    "renal clearance and volume, combined proportional plus additive residual",
    "error. CFTR modulator use, CFTR genotype and CF-related complications",
    "were screened and had no substantial effect.",
    sep = " "
  )
  reference <- paste(
    "Rolsma SL, Sokolow A, Patel PC, Sokolow K, Jimenez-Truque N, Fissell WH,",
    "Ryan V, Kirkpatrick CM, Nation RL, Gu K, Teresi M, Fishbane N, Kontos M,",
    "An G, Winokur P, Landersdorfer CB, Creech CB (2025). Population",
    "Pharmacokinetic Modeling of Cefepime, Meropenem, and",
    "Piperacillin-Tazobactam in Patients With Cystic Fibrosis. The Journal of",
    "Infectious Diseases 231(2):e364-e374. doi:10.1093/infdis/jiae451.",
    "Final parameter estimates from Supplemental Table 6B.",
    sep = " "
  )
  vignette <- "Rolsma_2025_betalactams_cysticfibrosis"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
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
        "Raw mL/min, NOT BSA-normalized, matching the raw-mL/min variant of",
        "the canonical CRCL column (Georges_2009_ceftazidime.R /",
        "Delattre_2010_amikacin.R / Chen_2023_nemonoxacin.R precedents). The",
        "effect on renal clearance is a through-origin linear normalization to",
        "a 90 mL/min reference, read from the Supplemental Table 6B unit",
        "string 'L/h per 90 mL/min CLCR,LBW'; footnote a of the same table",
        "('Sum of CL_R and CL_NR; for the typical CL_CR and TBW in the study",
        "population' = 10.5 L/h) confirms the reference by reproducing",
        "6.44 + 4.01 = 10.45. No exponent parameter is reported, so the",
        "relationship is linear rather than a power model. The 90 mL/min",
        "reference is also the cut point of the published",
        "target-attainment analysis, which contrasts CLCR,LBW below 90 with",
        "at or above 90 mL/min (Table 1, Figure 3): 'increased creatinine",
        "clearance led to reduced PTA' (Abstract)."
      ),
      source_name        = "CLCR,LBW"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source name TBW. Drives allometric scaling of the NON-RENAL",
        "clearance arm only, with a 61 kg reference and an exponent of 0.75,",
        "read from the Supplemental Table 6B unit string 'L/h/61kg^0.75 TBW'.",
        "The 61 kg reference matches the cohort median total body weight of",
        "61.30 kg in Supplemental Table 3, and footnote a describes the total",
        "as evaluated 'for the typical CL_CR and TBW in the study",
        "population'. The 0.75 exponent is printed inside the unit string and",
        "has no row of its own in the table, hence no point estimate, no RSE",
        "and no confidence interval; it is therefore encoded as fixed at the",
        "theoretical allometric value rather than estimated. Results 'PopPK",
        "Modeling': 'total body weight affected the nonrenal component of the",
        "meropenem elimination clearance'; the Discussion adds 'Comparing our",
        "results with a previously published meropenem study, we can confirm",
        "total body weight as an important covariate, but we also identified",
        "CLCR as an important driver of clearance'. Cohort 62.21 +/- 19.51 kg",
        "(median 61.30, range 36.3 to 137.0) per Supplemental Table 3."
      ),
      source_name        = "TBW"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "BMI bands the non-renal clearance arm in the source's fitted model,",
        "but the source instructs that the high-BMI arm must not be used for",
        "prediction, so no BMI effect is encoded here and the packaged model",
        "reproduces the BMI < 45 kg/m^2 arm only. Supplemental Table 6B",
        "reports two non-renal clearances: 'CL_NR for BMI <45 kg/m2' = 4.01",
        "L/h per 61 kg^0.75 TBW (RSE 34%), which is the arm encoded in this",
        "file, and 'CL_NR for BMI 45-50' = 0.16 L/h (RSE 2723%) carrying",
        "footnote b verbatim: 'based on 2 subjects, not for use in",
        "predictions, this parameter was included to avoid bias of the",
        "covariate relationship for BMI<45 kg/m2'. The 2723% relative standard",
        "error confirms the parameter is unidentified. It is a nuisance term",
        "protecting the BMI < 45 relationship from two outlying subjects, not",
        "a predictive covariate effect, and encoding it would make a simulated",
        "BMI 47 patient's non-renal clearance drop about 25-fold. The value is",
        "recorded here and in the validation vignette Errata so nothing from",
        "the source is lost. 96% of meropenem enrollments were adults and the",
        "cohort BMI was 22.21 +/- 6.58 kg/m^2 (median 20.45, range 14.8 to",
        "50.3) per Supplemental Table 3, so only the extreme upper tail of the",
        "cohort fell in the 45-50 band."
      )
    ),
    FFM = list(
      description = "Fat-free mass by the Janmahasatian formula (reported by the source as lean body weight)",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as a size descriptor and NOT retained as a direct covariate",
        "on any meropenem parameter -- total body weight described the",
        "non-renal clearance arm instead (Results 'PopPK Modeling': 'total",
        "body weight affected the nonrenal component of the meropenem",
        "elimination clearance'), and the volume of distribution took no size",
        "covariate at all. Lean body weight nevertheless enters the model",
        "indirectly, as the weight input to the CLCR,LBW covariate above.",
        "Fat-free mass IS a direct covariate in the companion cefepime model",
        "(Rolsma_2025_cefepime.R), where it scales the central volume. Cohort",
        "45.27 +/- 11.00 kg (median 49.15, range 26.7 to 63.5) per",
        "Supplemental Table 3."
      )
    ),
    CFTRMOD = list(
      description = "CFTR modulator use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened and explicitly rejected. Results 'PopPK Modeling': 'Other",
        "covariates, such as sample acquisition site, CFTR mutation, use of",
        "CFTR modulators ..., or complications of CF (eg, diabetes) did not",
        "have substantial impact on overall PK.' 48% of meropenem enrollments",
        "reported any modulator use (Supplemental Table 2). No point estimate",
        "is reported, so no effect can be encoded."
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods 'PopPK Model Development' covariate list) but not",
        "retained. The meropenem cohort is almost entirely adult -- children",
        "under 17 years were 6% of meropenem cases -- so the age range",
        "supporting this model is much narrower than the cefepime one. Cohort",
        "30.9 +/- 11.0 years (median 29.5, range 13 to 65) per Supplemental",
        "Table 3."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 1.68 +/-",
        "0.28 m^2 (median 1.70, range 1.2 to 2.5) per Supplemental Table 3."
      )
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. Cohort 167.01",
        "+/- 10.97 cm (median 167.65, range 145.0 to 188.0) per Supplemental",
        "Table 3."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened (Methods covariate list) but not retained. 48% female among",
        "meropenem enrollments (24 of 50) per Supplemental Table 2. Sex",
        "nevertheless enters the model indirectly, as an input to both the",
        "Janmahasatian fat-free-mass equation and the Cockcroft-Gault",
        "creatinine-clearance equation behind CLCR,LBW."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    n_enrollments  = 50L,
    age_range      = "13 to 65 years",
    age_median     = "29.5 years (mean 30.9 +/- 11.0)",
    weight_range   = "36.3 to 137.0 kg total body weight",
    weight_median  = "61.30 kg total body weight (mean 62.21 +/- 19.51); lean body weight median 49.15 kg (mean 45.27 +/- 11.00)",
    sex_female_pct = 100 * 24 / 50,
    race_ethnicity = c(White = 98, `Multi-Racial` = 2, `Hispanic or Latino (ethnicity)` = 2),
    disease_state  = paste(
      "Cystic fibrosis, admitted for a pulmonary exacerbation or for",
      "microbial eradication therapy. 96% of the meropenem group carried at",
      "least one copy of the DF508 mutation. 48% reported any CFTR modulator",
      "use. Participants receiving meropenem had more CF-related",
      "complications than the cefepime group (CF-related diabetes, hepatic",
      "disease, pancreatic insufficiency, and 10 or more pulmonary",
      "exacerbations in the previous 2 years; Supplemental Table 4).",
      "Staphylococcus aureus and Pseudomonas aeruginosa were the most common",
      "sputum isolates. Participants with a history of solid-organ or",
      "hematologic transplantation were excluded."
    ),
    renal_function = paste(
      "Estimated creatinine clearance on lean body weight (Cockcroft-Gault",
      "above 12 years of age, Schwartz bedside at or below 12 years). The",
      "structural reference is 90 mL/min per the Supplemental Table 6B unit",
      "string, which is also the cut point of the published",
      "target-attainment subgroups (CLCR,LBW below 90 vs at or above 90",
      "mL/min). The cohort distribution of CLCR,LBW is not tabulated in the",
      "paper or supplement."
    ),
    dose_range     = paste(
      "Intravenous infusion, standard of care as ordered by the treating",
      "team. The majority of doses were 2000 mg (n = 918); the remaining 80",
      "doses were 500 to 1800 mg. Most doses in participants under 17 years",
      "were 39 to 41 mg/kg of total body weight (n = 46), the remainder 33",
      "mg/kg (n = 30). Infusion durations in hours were 0.05 (3 doses), 0.5",
      "(85 doses) and 3 (910 doses)."
    ),
    regions        = "United States (Vanderbilt University Medical Center, Nashville TN; University of Iowa Hospital, Iowa City IA).",
    notes          = paste(
      "Opportunistic sampling during hospitalization, January 2018 to March",
      "2020. 192 of the 667 total plasma samples in the study were meropenem.",
      "Assay LC-MS/MS, linear 0.1 to 150 mg/L for meropenem. NONMEM 7.4.3",
      "with FOCE+I; outliers removed per the FDA popPK guidance (February",
      "2022) at a base-model weighted residual greater than 5. One- and",
      "two-compartment models were evaluated and the one-compartment model",
      "with linear elimination was selected. Drug administration was modeled",
      "as a zero-order (infusion) process. Model evaluation by",
      "prediction-corrected VPC (Supplemental Figure 2). Baseline demographics",
      "from Supplemental Tables 2 and 3; final parameter estimates from",
      "Supplemental Table 6B. Protein binding of 2% was assumed for meropenem",
      "when converting simulated total concentrations to the unbound",
      "concentrations used in the target-attainment analysis; that factor is",
      "a post-processing step and is not part of this model. NOTE: the source",
      "also reports a second, explicitly non-predictive non-renal clearance",
      "for the BMI 45-50 kg/m^2 band (0.16 L/h, RSE 2723%, based on 2",
      "subjects) which is documented in covariatesDataExcluded$BMI and in the",
      "vignette Errata but deliberately not encoded."
    )
  )

  ini({
    # Structural parameters -- Rolsma 2025 Supplemental Table 6B, 'Population
    # estimate' column. Total clearance is additive across a renal and a
    # non-renal arm; the table reports both arms plus their sum.
    #
    # Supplemental Table 6B verbatim (Parameter / Unit / Population estimate):
    #   CL_R                    / L/h per 90 mL/min CLCR,LBW / 6.44 (RSE 26%, IIV 21%)
    #   CL_NR for BMI <45 kg/m2 / L/h/61kg^0.75 TBW          / 4.01 (RSE 34%)
    #   CL_T for BMI <45 kg/m2  / L/h                        / 10.5 a
    #   CL_NR for BMI 45-50     / L/h                        / 0.16 b   <- NOT encoded, see below
    #   V                       / L                          / 17.5 (RSE 7.7%, IIV 14%)
    #   CVCP                    / -                          / 41%
    #   SDCP                    / mg/L                       / 0.315
    #   footnote a: 'Sum of CL_R and CL_NR; for the typical CL_CR and TBW in
    #               the study population'
    #   footnote b: 'based on 2 subjects, not for use in predictions, this
    #               parameter was included to avoid bias of the covariate
    #               relationship for BMI<45 kg/m2'
    #
    # The unit strings carry the covariate model: 'per 90 mL/min CLCR,LBW' is a
    # through-origin linear normalization and '/61kg^0.75 TBW' is allometric
    # scaling on a 61 kg reference. Footnote a is a mass-balance check on that
    # reading: 6.44 + 4.01 = 10.45, which rounds to the published 10.5 L/h.
    #
    # The BMI 45-50 non-renal clearance is deliberately NOT encoded, on the
    # source's own instruction that it is 'not for use in predictions'; see
    # covariatesDataExcluded$BMI for the full rationale and the retained value.
    lcl_renal  <- log(6.44); label("Renal clearance at CLCR,LBW = 90 mL/min (L/h)")                 # Suppl Table 6B row CL_R = 6.44 L/h (RSE 26%)
    lcl_nonren <- log(4.01); label("Non-renal clearance at WT = 61 kg, BMI < 45 kg/m^2 (L/h)")      # Suppl Table 6B row 'CL_NR for BMI <45 kg/m2' = 4.01 L/h (RSE 34%)
    lvc        <- log(17.5); label("Central volume of distribution (L)")                            # Suppl Table 6B row V = 17.5 L (RSE 7.7%)

    # Allometric exponent on the non-renal clearance arm. Printed inside the
    # Supplemental Table 6B unit string 'L/h/61kg^0.75 TBW' with no row, no
    # point estimate and no RSE of its own, i.e. fixed a priori at the
    # theoretical allometric value rather than estimated.
    e_wt_cl_nonren <- fixed(0.75); label("Allometric exponent of total body weight on non-renal clearance (unitless)")  # Suppl Table 6B unit string for CL_NR

    # Interindividual variability. The Supplemental Methods state the IIV model
    # explicitly ('The inter-individual variability (IIV) of the PK parameters
    # for all antibiotics in this study was estimated using an exponential
    # model'), so the tabulated percentages are CV% on the log-normal scale and
    # convert as omega^2 = log(CV^2 + 1). Supplemental Table 6B reports no
    # correlation or covariance cell for meropenem, so the two etas are
    # independent (unlike the cefepime, piperacillin and tazobactam models).
    etalcl_renal ~ 0.043155  # Suppl Table 6B: IIV CL_R 21% CV -> log(0.21^2+1) = 0.043155 (shrinkage 17%)
    etalvc       ~ 0.019410  # Suppl Table 6B: IIV V   14% CV -> log(0.14^2+1) = 0.019410 (RSE 124%, shrinkage 57%)

    # Residual error. Supplemental Methods: 'To describe the residual
    # variability (RV) of each antibiotic, a combined additive and proportional
    # error model was utilized.'
    propSd <- 0.41;  label("Proportional residual error (fraction)")  # Suppl Table 6B row CVCP = 41%     (RSE 14%)
    addSd  <- 0.315; label("Additive residual error (mg/L)")          # Suppl Table 6B row SDCP = 0.315 mg/L (RSE 434%)
  })

  model({
    # 1. Individual PK parameters.
    #
    # Renal clearance arm: through-origin linear normalization to the 90 mL/min
    # CLCR,LBW reference. The IIV eta sits on this arm alone, matching the
    # placement of the 21% IIV on the CL_R row of Supplemental Table 6B.
    cl_renal  <- exp(lcl_renal + etalcl_renal) * (CRCL / 90)
    # Non-renal clearance arm: allometric on total body weight, 61 kg
    # reference, exponent fixed at 0.75. No reported IIV. This is the
    # BMI < 45 kg/m^2 arm; the source's separate BMI 45-50 value is
    # deliberately not encoded (see covariatesDataExcluded$BMI).
    cl_nonren <- exp(lcl_nonren) * (WT / 61)^e_wt_cl_nonren
    cl        <- cl_renal + cl_nonren

    # Central volume: no covariate retained for meropenem.
    vc <- exp(lvc + etalvc)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment disposition with linear elimination. Meropenem is
    # given by intravenous infusion, so dose records enter the central
    # compartment directly as a zero-order input (Methods: 'Since all
    # antibiotics were administered by intravenous infusions at a constant
    # rate, drug administration was modeled as a zero-order process'); there is
    # no depot.
    d/dt(central) <- -kel * central

    # 4. Observation. Total (not unbound) plasma meropenem in mg/L, matching
    # the LC-MS/MS assay the model was fitted to. Combined proportional plus
    # additive residual error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
