Rolsma_2026_cefepime <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for",
    "intravenously infused cefepime in critically ill children (Rolsma 2026;",
    "84 patients, 223 plasma concentrations, ages 0.1 to 26.8 years, weights",
    "3.5 to 143 kg, admitted to a neonatal, pediatric or pediatric cardiac",
    "intensive care unit). Clearance is 3.51 L/h and central volume 19.53 L",
    "at the 70 kg / 90.7 mL/min reference, with body weight entering all four",
    "disposition parameters through fixed allometric exponents (0.75 on CL",
    "and Q, 1 on V1 and V2) and creatinine clearance entering clearance as an",
    "estimated power term (CrCl / 90.7)^0.67. Interindividual variability is",
    "diagonal on CL, V1 and V2 (none on Q, which could not be estimated",
    "precisely) with a proportional residual error. ECMO, CRRT, sex, race,",
    "ethnicity and gestational age were screened and not retained.",
    sep = " "
  )
  reference <- paste(
    "Rolsma SL, Blackman M, Beck C, Dbouk A, Morovia F, Glover F, Ess G,",
    "Dohler K, Bridges B, Choi L, Creech CB (2026). Population",
    "Pharmacokinetic Modeling and Target Attainment of Cefepime in Critically",
    "Ill Pediatric Patients. Open Forum Infectious Diseases 13(2):ofag069.",
    "doi:10.1093/ofid/ofag069. PMC12947848. Final parameter estimates from",
    "Supplemental Table 3, column 'Base+Weight+CrCl 2 Cpt Model'",
    "(Obj = 2159.81); model structure from the displayed equation block in",
    "Results, 'Population Pharmacokinetic Analysis'.",
    sep = " "
  )
  vignette <- "Rolsma_2026_cefepime"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: cefepime is given as an intravenous infusion and TOTAL (not
  # unbound) plasma cefepime was assayed by LC-MS/MS at the Vanderbilt Mass
  # Spectrometry Core (Rolsma 2026 Methods, "Quantification of Plasma Drug
  # Concentrations"; Supplementary material text and Supplemental Table 1).
  compartmentData <- list(
    central     = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters all four disposition parameters as the allometric ratio",
        "(WT / 70) raised to a FIXED exponent: 0.75 on the clearances CL and",
        "Q, and 1 on the volumes V1 and V2. Rolsma 2026 Methods states the",
        "scaling explicitly ('starting from modeling weight on PK parameters,",
        "which was allometrically scaled to a factor of 70 kg'), and",
        "Supplemental Table 4 writes the exponents out for all 14 candidate",
        "models. The exponents carry no standard error, confidence interval",
        "or %CV anywhere in Supplemental Table 3, which is the signal that",
        "they were fixed rather than estimated; they are therefore wrapped in",
        "fixed() in ini(). NOTE the 70 kg standardisation is far outside this",
        "cohort: median weight is 13.7 kg (range 3.5-143.0, Table 1), so",
        "3.51 L/h and 19.53 L are extrapolated reference values rather than",
        "typical values for a patient in the study. See the vignette",
        "Assumptions and deviations for the main-text 0.74 / supplement 0.75",
        "discrepancy on the CL exponent.",
        sep = " "
      ),
      source_name        = "wt"
    ),
    CRCL = list(
      description        = "Estimated creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW creatinine clearance in mL/min -- NOT BSA-normalized to",
        "mL/min/1.73 m^2 and NOT weight-normalized. Rolsma 2026 states the",
        "units directly beneath the model equation block ('CrCl_i is",
        "creatinine clearance in mL/min for individual i'), and the cohort",
        "values confirm an absolute scale: median 90.7 mL/min, range",
        "10.8-275.0 in children whose median weight is 13.7 kg (Table 1). A",
        "BSA- or weight-normalized value would not span that range in this",
        "population. The paper does NOT state which estimating equation was",
        "used (no Schwartz / Cockcroft-Gault attribution appears in Methods,",
        "Results or the supplement), so the assay form is unknown; this is",
        "recorded in the vignette Assumptions and deviations. Enters clearance",
        "as the estimated power term (CRCL / 90.7)^0.67. The 90.7 reference is",
        "the cohort MEDIAN (Table 1), consistent with Methods: 'Continuous",
        "covariates, such as gestational age and creatinine clearance were",
        "standardized to the cohort's median values.' Must be strictly",
        "positive -- the covariate enters a power term. Body size is handled",
        "separately by the allometric WT term, so CRCL carries renal function",
        "only and the two are not double-counting size.",
        sep = " "
      ),
      source_name        = "CrCl"
    )
  )

  # Screened in the a priori covariate analysis but NOT retained in the final
  # model, so these are documentation only and are deliberately absent from
  # model(). Rolsma 2026 Results: "We observed that receipt of ECMO, sex, and
  # race did not improve the model fit. The use of the age maturation function
  # with gestational age also did not improve the model fit. Receipt of CRRT,
  # ethnicity, and gestational age did individually improve model fit, but
  # their inclusion in the model with weight and creatinine clearance did not
  # provide evidence of being better than the model with only weight and
  # creatinine clearance." Per-covariate objective function values for every
  # candidate are in Supplemental Table 4 (main cohort) and Supplemental
  # Table 5 (complete-data sensitivity analysis).
  covariatesDataExcluded <- list(
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation treatment-status indicator (1 = receiving ECMO)",
      units       = "(binary)",
      type        = "binary",
      reference_category = "0 (not receiving ECMO)",
      notes       = paste(
        "10 of 84 patients (11.9%, Table 1). Tested as an additive shift on",
        "CL; Supplemental Table 4 gives Obj = 2205.88 versus 2207.18 for the",
        "weight-only model, a negligible improvement. Not retained. The",
        "authors note in the Discussion that they had anticipated an ECMO",
        "effect but that these participants 'only represented a small subset",
        "of our population'. Contrast Valadez_2025_cefepime.R and",
        "Kang_2020_cefpirome.R, which do retain an ECMO term.",
        sep = " "
      ),
      source_name = "ECMO"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy treatment-status indicator (1 = receiving CRRT)",
      units       = "(binary)",
      type        = "binary",
      reference_category = "0 (not receiving CRRT)",
      notes       = paste(
        "9 of 84 patients (10.7%, Table 1). Individually improved fit",
        "(Supplemental Table 4 Obj = 2211.02 as an additive shift on CL) but",
        "was not retained once weight and creatinine clearance were in the",
        "model. Patients receiving intermittent hemodialysis were excluded at",
        "enrollment; CRRT patients were not.",
        sep = " "
      ),
      source_name = "CRRT"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female)",
      units       = "(binary)",
      type        = "binary",
      reference_category = "0 (male)",
      notes       = paste(
        "41 female / 43 male (48.8% female, Table 1). Supplemental Table 4",
        "tested the effect as 'Sex = Male' (i.e. male as the indicator, the",
        "opposite orientation from the canonical SEXF); Obj = 2209.55, not",
        "retained. Because the term was dropped, no coefficient sign needs to",
        "be flipped for the canonical orientation.",
        sep = " "
      ),
      source_name = "Sex"
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      reference_category = "Black or African American (the omitted category in the source's two-indicator parameterisation)",
      notes       = paste(
        "68 of 84 (81.0%, Table 1). Supplemental Table 4 tested race as the",
        "two indicators 'Race = Other' and 'Race = White' simultaneously",
        "(Obj = 2206.80); not retained. One patient had missing race, imputed",
        "with the modal category.",
        sep = " "
      ),
      source_name = "Race"
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator",
      units       = "(binary)",
      type        = "binary",
      reference_category = "Black or African American",
      notes       = paste(
        "The source's 'Other' arm pools Native Hawaiian or other Pacific",
        "Islander (1, 1.2%) with Other race (3, 3.6%), i.e. 4 of 84 (4.8%)",
        "(Table 1). Paired with RACE_WHITE in the Supplemental Table 4",
        "candidate; not retained.",
        sep = " "
      ),
      source_name = "Race"
    ),
    RACE_HISPANIC = list(
      description = "Hispanic or Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      reference_category = "0 (not Hispanic or Latino)",
      notes       = paste(
        "9 of 84 (10.7%, Table 1). Individually improved fit (Supplemental",
        "Table 4 Obj = 2202.48 as 'Eth = Hispanic/Latino'), and even in",
        "combination with weight and creatinine clearance gave the lowest",
        "objective function of any candidate (Obj = 2162.64 vs 2159.81 for the",
        "final model -- note the final model is LOWER, so ethnicity did not",
        "improve on it). Not retained. One patient had missing ethnicity,",
        "imputed with the modal category.",
        sep = " "
      ),
      source_name = "Ethnicity"
    ),
    GA = list(
      description = "Gestational age at birth",
      units       = "days in the source (the canonical GA column is in weeks; divide by 7)",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Median 274 days (39.1 weeks), range 171-280 days, mean 263 (SD 23)",
        "(Table 1). Tested two ways, neither retained: as a continuous power",
        "term (GA / 274) on CL (Supplemental Table 4 Obj = 2211.74), and",
        "inside a Hill maturation function 1 / (1 + (TM50 / (GA/7))^Hill)",
        "(Obj = 2204.48), where the GA/7 converts the source's days to the",
        "weeks the maturation model expects. One patient had missing",
        "gestational age, imputed as full term (37 weeks).",
        sep = " "
      ),
      source_name = "GA"
    ),
    AGE = list(
      description = "Postnatal age",
      units       = "years",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Median 2.61 years, range 0.1-26.8, mean 6.6 (SD 7.2) (Table 1).",
        "Tested as a power term (Age / 2.61) on CL centered at the cohort",
        "median; Supplemental Table 4 Obj = 2191.59 alone and 2166.85 combined",
        "with creatinine clearance, both worse than the final model's 2159.81.",
        "Not retained. Age is not named in the Methods list of a priori",
        "covariates but does appear in the Supplemental Table 4 candidate set.",
        sep = " "
      ),
      source_name = "Age"
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL in the source (the canonical ALB column is g/L; multiply by 10)",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Median 2.9 g/dL, range 1.6-4.8, missing for 8 of 84 (9.5%) (Table 1).",
        "Excluded from the primary covariate analysis for missingness; tested",
        "only in the complete-data sensitivity analysis as (Albumin / 2.9) on",
        "CL (Supplemental Table 5 Obj = 1853.49 vs 1853.52 for weight alone).",
        "Not retained.",
        sep = " "
      ),
      source_name = "Albumin"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Median 39 U/L, range 13-1450, missing for 15 of 84 (17.9%)",
        "(Table 1). Excluded from the primary analysis for missingness; tested",
        "in the sensitivity analysis as (AST / 39) on CL (Supplemental Table 5",
        "Obj = 1846.57 alone, 1818.49 with creatinine clearance, both worse",
        "than the sensitivity final model's 1811.09). Not retained.",
        sep = " "
      ),
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Median 26 U/L, range 6-745, missing for 15 of 84 (17.9%) (Table 1).",
        "Excluded from the primary analysis for missingness; tested in the",
        "sensitivity analysis as (ALT / 26) on CL (Supplemental Table 5",
        "Obj = 1846.04 alone, 1817.81 with creatinine clearance). Not retained.",
        sep = " "
      ),
      source_name = "ALT"
    ),
    TPRO = list(
      description = "Total serum protein",
      units       = "g/dL in the source (the canonical TPRO column is g/L; multiply by 10)",
      type        = "continuous",
      reference_category = NULL,
      notes       = paste(
        "Median 5.3 g/dL, range 3.2-9.0, missing for 15 of 84 (17.0%)",
        "(Table 1). Excluded from the primary analysis for missingness; tested",
        "in the sensitivity analysis as (Total protein / 5.2) on CL",
        "(Supplemental Table 5 Obj = 1847.64 alone, 1817.99 with creatinine",
        "clearance). Not retained. Note the sensitivity analysis centres on",
        "5.2 g/dL while Table 1 reports a cohort median of 5.3 g/dL; the",
        "difference is immaterial because the term was dropped.",
        sep = " "
      ),
      source_name = "Total Protein"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 84L,
    n_studies      = 1L,
    n_samples      = 223L,
    age_range      = "0.1-26.8 years (Table 1). Despite the pediatric framing the cohort reaches adulthood; eligibility required age over 1 month with no upper bound.",
    age_median     = "2.61 years (mean 6.6, SD 7.2; the distribution has a strong right skew)",
    weight_range   = "3.5-143.0 kg (Table 1)",
    weight_median  = "13.7 kg (mean 24.7, SD 25.5)",
    height_range   = "49.0-184.0 cm, median 87.0 (mean 100.0, SD 41.6) (Table 1)",
    sex_female_pct = 48.8,
    race_ethnicity = c(White = 81.0, `Black or African American` = 13.1, `Other race` = 3.6, `Native Hawaiian or other Pacific Islander` = 1.2, Missing = 1.2),
    ethnicity      = "Hispanic or Latino 10.7%, not Hispanic or Latino 88.1%, missing 1.2% (Table 1)",
    disease_state  = "Critically ill children admitted to the Neonatal Intensive Care Unit, Pediatric Critical Care Unit, or Pediatric Cardiac Intensive Care Unit and receiving cefepime as standard of care. 10 of 84 (11.9%) were supported by ECMO and 9 of 84 (10.7%) by CRRT. Patients who were pregnant, receiving hemodialysis or probenecid, or expected to die within 72 hours were excluded.",
    renal_function = "Creatinine clearance median 90.7 mL/min, range 10.8-275.0, mean 105.0 (SD 59.5) (Table 1). The range spans renal impairment through augmented renal clearance, which the Discussion identifies as a common and under-recognised phenomenon in critically ill children.",
    gestational_age = "Median 274 days (39.1 weeks), range 171-280 days, mean 263 (SD 23) (Table 1)",
    dose_range     = "Standard of care as ordered by the treating team; dose median 47 mg/kg, range 11.7-58.4 (mean 44.4, SD 9.00) (Table 1). Median 11 dosing events per patient (range 2-41).",
    regions        = "United States (single centre: Monroe Carell Jr. Children's Hospital at Vanderbilt University Medical Center, Nashville, Tennessee)",
    protein_binding = "Not fitted. Total plasma cefepime was assayed and is what this model outputs as Cc. The probability-of-target-attainment analysis converted total to unbound concentrations using a published protein binding of 20% (unbound fraction 0.80) taken from the literature, not estimated here (Rolsma 2026 Methods, 'Probability of Target Attainment Analysis'). Multiply Cc by 0.80 to obtain the free concentration the paper's fT>MIC targets are evaluated against.",
    sampling       = "Opportunistic residual sampling, prioritising pre-dose, end-of-infusion and 1-2 h post-infusion times. Median 2 concentrations per patient (range 1-10). Of the 223 concentrations, 16.2% were troughs (within 2 h pre-dose), 34.5% peaks (up to 2 h post-dose) and 49.3% random (Supplemental Table 2). Measured concentrations median 60.0 mg/L, range 3.7-201.0. Assay LLOQ 0.05 mg/L; all analysed samples were above it.",
    notes          = "Prospective popPK study of beta-lactam antibiotics with opportunistic sampling, enrolling October 2020 - November 2022. Four cefepime concentrations (1.8%) were excluded as suspected line contamination. Estimated in Monolix 2021R with SAEM. Final model shrinkage: CL 16.2%, V1 63.5%, V2 37.0% (Supplemental Table 6) -- the high V1 shrinkage means individual central volumes are poorly informed by these sparse data, so treat subject-level V1 predictions with caution. A complete-data sensitivity analysis (Supplemental Table 5) selected the same final model."
  )

  ini({
    # ---------------------------------------------------------------------
    # All final estimates are the third column of Rolsma 2026 Supplemental
    # Table 3, headed "Base+Weight+CrCl 2 Cpt Model" with Obj = 2159.81 --
    # the objective function value that Supplemental Table 4 also assigns to
    # the retained weight-plus-creatinine-clearance candidate, confirming the
    # column is the final model. Values are printed there as
    # "Est (SE) [Wald 95% CI]"; only the point estimate is carried here.
    #
    # theta numbering below follows the paper's own labelling in the text
    # under the equation block: "theta_1, theta_3, theta_4, and theta_5 being
    # population PK parameters of CL (L/hour), V1 (L), Q (L/hour), and V2 (L),
    # respectively, and theta_2 being the model parameter for CrCl".
    # ---------------------------------------------------------------------

    lcl <- log(3.51)
    label("Clearance at WT = 70 kg and CRCL = 90.7 mL/min (L/h)")
    # theta_1. Supplemental Table 3: 3.51 (0.29) [2.95, 4.07] L/h.

    lvc <- log(19.53)
    label("Central volume of distribution at WT = 70 kg (L)")
    # theta_3, the paper's V1. Supplemental Table 3: 19.53 (3.83)
    # [12.03, 27.03] L.

    lq <- log(0.70)
    label("Intercompartmental clearance at WT = 70 kg (L/h)")
    # theta_4. Supplemental Table 3: 0.70 (0.05) [0.60, 0.81] L/h.

    lvp <- log(3.56)
    label("Peripheral volume of distribution at WT = 70 kg (L)")
    # theta_5, the paper's V2. Supplemental Table 3: 3.56 (2.34)
    # [-1.03, 8.15] L. The lower confidence limit is negative, i.e. the
    # peripheral volume is estimated imprecisely; this is consistent with the
    # 80% CV and wide interval on its interindividual variability below and
    # with the sparse opportunistic sampling. Carried as published.

    e_crcl_cl <- 0.67
    label("Power exponent on (CRCL / 90.7) for CL (unitless)")
    # theta_2. Supplemental Table 3: 0.67 (0.04) [0.60, 0.75]. Estimated,
    # not fixed -- it is the only exponent in the model that Supplemental
    # Table 3 reports with a standard error and confidence interval, and
    # adding it to the weight model dropped Obj from 2207.18 to 2159.81
    # (Supplemental Table 4).

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on (WT / 70) shared by CL and Q (unitless)")
    # Supplemental Table 4 writes "CLi = Clpop (wt/70)0.75" and
    # "Qi = Qpop (wt/70)0.75" for every one of the 14 candidate models,
    # including the retained weight-plus-creatinine-clearance row
    # (Obj = 2159.81). Fixed, not estimated: Supplemental Table 3 reports no
    # standard error, confidence interval or %CV for any exponent, and
    # Methods describes standard allometric scaling to a 70 kg factor.
    #
    # DISCREPANCY, resolved in favour of 0.75. The main-text equation block
    # in Results prints the CL exponent as 0.74 while printing the Q exponent
    # as 0.75 in the very same block. 0.75 is used here because: (a)
    # Supplemental Table 4 states 0.75 fourteen times, including for this
    # exact final model; (b) Supplemental Table 3 has no estimated-exponent
    # row, so the exponent cannot be a fitted 0.74; and (c) fixing one
    # clearance at 0.74 and the other at 0.75 has no mechanistic reading.
    # The effect of the discrepancy is small (1.6% on CL at the 13.7 kg
    # cohort median, 3% at the 3.5 kg minimum). See the vignette Assumptions
    # and deviations.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on (WT / 70) shared by V1 and V2 (unitless)")
    # Supplemental Table 4: "V1i = V1pop (wt/70)1" and
    # "V2i = V2pop (wt/70)1" for every candidate model; the main-text
    # equation block agrees, printing the exponent 1 on both volumes. Fixed
    # for the same reason as e_wt_cl_q -- no uncertainty is reported for it.

    # ---------------------------------------------------------------------
    # Interindividual variability. Supplemental Table 3 reports these as
    # %CV, and the paper states the random effects follow "a log-normal
    # distribution with means of zero and a covariance matrix that has
    # omega^2_CL, omega^2_V1, and omega^2_V2 along the diagonal". The
    # log-scale variances below are therefore omega^2 = log(CV^2 + 1), the
    # exact inverse of the CV that Monolix (used here, 2021R) reports for a
    # log-normally distributed random effect. Round-trip check:
    #   sqrt(exp(0.1415864) - 1) = 0.39   -> 39%
    #   sqrt(exp(0.0392207) - 1) = 0.20   -> 20%
    #   sqrt(exp(0.4946962) - 1) = 0.80   -> 80%
    #
    # There is NO random effect on Q. Rolsma 2026 Results: "The random
    # effects on all PK parameters except Q were included, since the random
    # effect of Q could not be estimated precisely, and its inclusion did not
    # improve the model fit." The matrix is diagonal -- no correlations are
    # reported in Supplemental Table 3.
    # ---------------------------------------------------------------------

    etalcl ~ 0.1415864
    # omega^2 for CL. Supplemental Table 3: 39 (8) [24, 54] %CV.

    etalvc ~ 0.0392207
    # omega^2 for V1. Supplemental Table 3: 20 (2) [-23, 64] %CV. The printed
    # standard error of 2 is inconsistent with its own confidence interval,
    # which implies about 22 (the interval is centred on 20.5 with a
    # half-width of 43.5, i.e. 43.5 / 1.96 = 22.2); every other row in the
    # table reproduces its interval as estimate +/- 1.96 SE. Only the point
    # estimate is carried into the model, so the typo does not affect the
    # encoding. Note also the high V1 shrinkage (63.5%, Supplemental
    # Table 6).

    etalvp ~ 0.4946962
    # omega^2 for V2. Supplemental Table 3: 80 (108) [-131, 202] %CV --
    # estimated very imprecisely, matching the wide interval on theta_5. This
    # row is also internally inconsistent: 80 +/- 1.96 * 108 is [-131, 292],
    # so the printed lower limit reproduces but the upper limit does not.
    # As with the V1 row, only the point estimate is carried.

    propSd <- 0.42
    label("Proportional residual error (fraction)")
    # Supplemental Table 3: proportional 42 (3) [37, 47] %CV. A proportional
    # error model was selected over additive and combined alternatives on the
    # objective function value (Methods, "Population Pharmacokinetic
    # Analysis").
  })

  model({
    # Covariate reference values, both printed in the Results equation block.
    # 70 kg is the conventional allometric standard; 90.7 mL/min is the
    # cohort median creatinine clearance (Table 1), per the Methods statement
    # that continuous covariates were standardized to the cohort medians.
    wt_ref   <- 70   # kg
    crcl_ref <- 90.7 # mL/min

    # CL_i = theta_1 * (wt_i/70)^0.75 * (CrCl_i/90.7)^theta_2 * exp(eta_CL_i)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl_q * (CRCL / crcl_ref)^e_crcl_cl

    # V1_i = theta_3 * (wt_i/70)^1 * exp(eta_V1_i)
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc_vp

    # Q_i = theta_4 * (wt_i/70)^0.75    -- no random effect, see ini()
    q <- exp(lq) * (WT / wt_ref)^e_wt_cl_q

    # V2_i = theta_5 * (wt_i/70)^1 * exp(eta_V2_i)
    vp <- exp(lvp + etalvp) * (WT / wt_ref)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # TOTAL plasma cefepime, matching the assayed quantity. The unbound
    # fraction of 0.80 used for the fT>MIC targets is a literature constant
    # applied outside the PK model and is recorded in
    # population$protein_binding rather than in ini().
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
