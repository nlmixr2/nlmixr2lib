Zhao_2025_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (HD-MTX, >= 0.5 g/m^2; standard 1.3 g/m^2 given as a 1 h",
    "infusion of one third of the dose followed by an 11 h infusion of the",
    "remainder) in Chinese patients with intracranial germ cell tumors",
    "(Zhao 2025; n = 505 patients, 73.3% children, 5,470 serum",
    "concentrations). Linear first-order elimination from the central",
    "compartment. Clearance is a product of four power-form covariate terms",
    "normalized to the cohort medians -- BSA-normalized eGFR (exponent",
    "0.23), body weight (0.39), total bilirubin (-0.05) and serum albumin",
    "(-0.18) -- and an exponential bleomycin co-medication term that raises",
    "CL 1.08-fold; the central volume carries a body-weight power term",
    "(exponent 0.31). Peripheral volume and intercompartmental clearance",
    "were fixed at the base-model estimates. Exponential between-subject",
    "variability on CL and Vc and a proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Zhao J, Wu R, Zhang S, Lu Q, Wang R, He Y, Zhao Z, Mei S (2025).",
    "Population pharmacokinetic model of high-dose methotrexate in Chinese",
    "patients with intracranial germ cell tumors.",
    "Front Pharmacol 16:1548203. doi:10.3389/fphar.2025.1548203.",
    sep = " "
  )
  vignette <- "Zhao_2025_methotrexate"
  units    <- list(time = "h", dosing = "umol", concentration = "umol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Zhao 2025 Results 3.2 ("A
  # two-compartment model with first-order elimination best described the
  # pharmacokinetics of MTX, with parameters including CL, central compartment
  # volume (Vc), peripheral compartment volume (Vp), and inter-compartmental
  # clearance (Q)") and Methods 2.3 ("Fasting venous blood (2 mL) was collected
  # ... without the addition of anticoagulants. The serum MTX concentrations
  # were assessed by enzyme-multiplied immunoassay"). Amount units are umol
  # because the assay and every reported concentration are in umol/L
  # (Methods 2.3 quantitative range 0.3-2,600 umol/L; Table 1 MTX
  # concentration median 1.2 umol/L, range 0.01-55.80 umol/L).
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "umol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "umol", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "BSA-normalized estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column eGFR. Methods 2.2: 'eGFR was determined via the 2008",
        "Schwartz bedside formula for children and the 2021 CKD-EPI formula",
        "for adults', so the estimating equation differs by age stratum",
        "within the same column. Power-form effect on CL normalized to the",
        "cohort median, per Zhao 2025 Eq. 5:",
        "'CL = 12.88 x (eGFR/102.2)^0.23 x ...'. The 102.2 mL/min/1.73 m^2",
        "reference is the cohort median reported in Table 1 (median 102.20,",
        "range 41.57-446.22). Positive exponent: better glomerular filtration",
        "raises MTX clearance, which is mechanistically expected because",
        "'Approximately 80%-90% of MTX is excreted via the kidneys'",
        "(Introduction / Discussion). eGFR was the single most influential",
        "covariate in the forward step (Table 2 model 2, delta OFV -122.058,",
        "P < 0.01) and the second-largest on backward elimination (Table 2",
        "model 8, delta OFV +162.988, P < 0.001). It was preferred over serum",
        "creatinine, which gave a smaller drop (delta OFV -61.22; Results",
        "3.2). Creatinine clearance was NOT tested at all: 'Due to the lack of",
        "urine data in our study, CLcr could not be accurately estimated in",
        "children. Therefore, it was not included as a tested covariate.'",
        "Cystatin-C-based eGFR was likewise unavailable (Discussion)."
      ),
      source_name        = "eGFR"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BW. Power-form effects on BOTH CL and Vc, each",
        "normalized to the cohort median, per Zhao 2025 Eqs. 5-6:",
        "'CL = ... x (BW/47)^0.39 x ...' and 'Vc = 72.04 x (BW/47)^0.31'.",
        "The 47 kg reference is the cohort median reported in Table 1 (median",
        "47.00, range 14.00-121.00) and restated in Results 3.2 ('Median",
        "values of eGFR, BW, TBIL, and ALB were 102.2 mL/min/1.73 m^2, 47 kg,",
        "15.3 umol/L, and 40.9 g/L, respectively'). Note both exponents are",
        "estimated, not fixed at the allometric 0.75 / 1.0 defaults. Body",
        "surface area was an equally significant alternative on CL (delta OFV",
        "-121.968 for BW vs -116.073 for BSA); Results 3.2: 'Due to",
        "collinearity between the two covariates ... BW outperformed BSA in",
        "-2LL (10,421.315 vs. 10,423.642), AIC (10,443.315 vs. 10,445.642),",
        "and BIC (10,515.992 vs. 10,518.319), therefore, BW was selected as",
        "the final covariate.' BSA is nonetheless the dosing descriptor for",
        "this protocol (1.3 g/m^2), so a simulation must carry BSA for the",
        "dose amount even though only WT enters the model."
      ),
      source_name        = "BW"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column TBIL. Power-form effect on CL normalized to the cohort",
        "median, per Zhao 2025 Eq. 5: 'CL = ... x (TBIL/15.3)^-0.05 x ...'.",
        "The 15.3 umol/L reference is the cohort median reported in Table 1",
        "(median 15.30, range 3.40-78.50) and restated in Results 3.2. SI",
        "units, matching the register canonical. Negative exponent: higher",
        "total bilirubin lowers MTX clearance. Discussion: 'Our study",
        "identified that increased TBIL was associated with decreased MTX CL.",
        "Nakano et al. found that elevated TBIL is an independent risk factor",
        "for delayed MTX elimination, which is consistent with our findings",
        "... Biliary excretion accounts for 10% or less of the administered",
        "dose of MTX.' Added sixth in the forward step (Table 2 model 6,",
        "delta OFV -32.741, P < 0.01; backward removal +29.909, P < 0.001).",
        "Note the printed effect size is small: over the full observed TBIL",
        "range 3.4-78.5 umol/L the multiplier spans only",
        "(3.4/15.3)^-0.05 = 1.08 to (78.5/15.3)^-0.05 = 0.92."
      ),
      source_name        = "TBIL"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column ALB. Power-form effect on CL normalized to the cohort",
        "median, per Zhao 2025 Eq. 5: 'CL = ... x (ALB/40.9)^-0.18'. The",
        "40.9 g/L reference is the cohort median reported in Table 1 (median",
        "40.90, range 30.60-52.30) and restated in Results 3.2. SI units,",
        "matching the register canonical. Negative exponent: LOWER albumin",
        "gives HIGHER MTX clearance, which is the opposite direction to the",
        "Mao 2022 and Pai 2020 models cited in the Discussion. The authors",
        "attribute the reversal to competing mechanisms:",
        "'MTX is about 50% protein-bound, and lower ALB levels may increase",
        "the unbound fraction of MTX, thereby enhancing its clearance and",
        "metabolism ... The increase in third-space fluid and the decrease in",
        "MTX protein binding were two counterbalancing effects on MTX",
        "clearance ... For our patients diagnosed with iGCTs, pleural",
        "effusions were less common. Therefore, the reduction in MTX protein",
        "binding was the predominant factor.' Added last in the forward step",
        "(Table 2 model 7, delta OFV -16.701, P < 0.01) and the weakest term",
        "on backward elimination (+16.701, P < 0.001)."
      ),
      source_name        = "ALB"
    ),
    CONMED_BLEOMYCIN = list(
      description        = "Concomitant bleomycin indicator (1 = co-medicated with bleomycin, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant bleomycin)",
      notes              = paste(
        "Source column BLM. Enters CL as an EXPONENTIAL term, not a power",
        "term, per Zhao 2025 Eq. 5: 'CL = ... x e^BLM x ...' with",
        "'BLM = 0.08 when combined with BLM, otherwise = 0'. Encoded here as",
        "exp(e_bleomycin_cl * CONMED_BLEOMYCIN) so the coefficient 0.08 and",
        "the 0/1 indicator are separated; exp(0.08) = 1.083, matching the",
        "paper's headline claim in the Abstract and Discussion",
        "('Co-administration with bleomycin could increase MTX CL by a factor",
        "of 1.08' / 'MTX CL increases by 1.08-fold when co-administered with",
        "BLM'). Table 1 records bleomycin on 941 (58.67%) of the",
        "concomitant-medication records; Results 3.1 states 'Concurrent use of",
        "BLM was observed in 58.67% of patients'. Treatment schedule (Methods",
        "2.1): 'MTX and vincristine on day 1, BLM on day 2 and cisplatin on",
        "day 3 or 4', so BLM begins about a day after the MTX infusion.",
        "Added fifth in the forward step (Table 2 model 5, delta OFV -73.634,",
        "P < 0.01; backward removal +79.676, P < 0.001). The mechanism is",
        "unexplained -- Discussion: 'the influence of BLM on MTX remains",
        "unexplored. Our study found that BLM slightly accelerates the",
        "excretion of MTX, and the underlying mechanism of this finding",
        "requires further investigation.'"
      ),
      source_name        = "BLM"
    )
  )

  # Covariates the paper collected and screened but did not retain. Results
  # 3.2: 'Other covariates did not show statistically significant effects and
  # were excluded from the model.' Kept here as documentation only -- they are
  # never referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Collected per Methods 2.2 (demographic information) and listed among",
        "the factors reported to influence MTX PK in the Introduction, but not",
        "retained. Table 1 median 14.00 years (range 3.00-48.00); Results 3.1",
        "records 26.7% adults and 73.3% children."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected per Methods 2.2 but not retained. Table 1 reports",
        "'Sex (Male/Female) 357/148'; Results 3.1 states '505 patients (357",
        "females and 148 males)'. The two statements disagree on which count",
        "is male -- see population$notes. Sex is not a model covariate either",
        "way, so the discrepancy does not affect any equation."
      )
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Collected per Methods 2.2 but not retained. Table 1 median 155.00 cm",
        "(range 103.00-194.00). Height is nonetheless an input to the 2008",
        "bedside Schwartz equation used to derive the retained CRCL covariate",
        "in children (Methods 2.2)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Collected per Methods 2.2 but not retained. Table 1 median 19.19",
        "kg/m^2 (range 11.14-36.53)."
      )
    ),
    BSA = list(
      description = "Body surface area (Mosteller equation)",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Significant on CL (delta OFV -116.073, P < 0.01) but deliberately",
        "superseded by body weight because the two are collinear and BW won on",
        "-2LL, AIC and BIC (Results 3.2). Methods 2.2: 'BSA was calculated",
        "using the Mosteller formula'. Table 1 median 1.42 m^2 (range",
        "0.63-2.47). BSA remains the dosing descriptor for the protocol",
        "(1.3 g/m^2, Methods 2.1), so simulations still need it to compute the",
        "dose amount."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and superseded by eGFR. Results 3.2: 'Compared to SCR",
        "(delta OFV = -61.22, P < 0.01), eGFR (delta OFV = -122.058,",
        "P < 0.01) was a more significant predictor of MTX CL'. Table 1 median",
        "62.00 umol/L (range 9.00-147.00)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Tested and not retained. Discussion: 'Consistent with previous",
        "studies ... we did not observe the effect of ALT, AST and ALP on PPK",
        "parameters of MTX.' Table 1 median 26.00 U/L (range 8.00-663.00)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Tested and not retained (same Discussion sentence as AST). Table 1",
        "median 17.00 U/L (range 3.00-530.00)."
      )
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Tested and not retained (same Discussion sentence as AST). Table 1",
        "median 110.00 U/L (range 20.00-593.00)."
      )
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Collected per Methods 2.2 (biological parameters) but not retained.",
        "Table 1 median 212.00 U/L (range 80.00-1,059.00)."
      )
    ),
    CONMED_NSAID = list(
      description = "Concomitant non-steroidal anti-inflammatory drug indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested and not retained. Discussion: 'non-steroidal",
        "anti-inflammatory drugs may reduce MTX CL by inhibiting renal",
        "prostaglandin synthesis, competing with protein binding, and",
        "interfering with organic anion-mediated renal excretion ... However,",
        "our study failed to find the influence of these covariates on model",
        "parameters during model development'. Table 1 records NSAIDs on",
        "1,580 (98.50%) of the concomitant-medication records -- near-universal",
        "exposure leaves almost no contrast to estimate an effect from."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 505L,
    n_studies      = 1L,
    age_range      = "3 to 48 years",
    age_median     = "14 years (26.7% adults, 73.3% children)",
    weight_range   = "14 to 121 kg",
    weight_median  = "47 kg",
    height_range   = "103 to 194 cm (median 155)",
    bmi_range      = "11.14 to 36.53 kg/m^2 (median 19.19)",
    bsa_range      = "0.63 to 2.47 m^2 (median 1.42; Mosteller equation)",
    sex_female_pct = 100 * 148 / 505,
    disease_state  = paste(
      "Intracranial germ cell tumors (iGCTs), a rare CNS neoplasm of",
      "adolescents peaking at ages 12-16 and accounting for about 3%-5% of",
      "paediatric primary CNS tumours. Treatment regimen (Methods 2.1): MTX",
      "and vincristine on day 1, bleomycin on day 2, and cisplatin on day 3",
      "or 4, with calcium folinate rescue. This is the first published popPK",
      "model of methotrexate in iGCTs (Discussion)."
    ),
    renal_function = paste(
      "eGFR median 102.20 mL/min/1.73 m^2 (range 41.57-446.22), computed with",
      "the 2008 bedside Schwartz equation in children and the 2021 CKD-EPI",
      "equation in adults; serum creatinine median 62 umol/L (range 9-147).",
      "Creatinine clearance could not be estimated in children for lack of",
      "urine data and was therefore never tested as a covariate (Results",
      "3.2); cystatin C was unavailable (Discussion)."
    ),
    hepatic_function = paste(
      "Total bilirubin median 15.30 umol/L (range 3.40-78.50); albumin median",
      "40.90 g/L (range 30.60-52.30); AST median 26 U/L (range 8-663); ALT",
      "median 17 U/L (range 3-530); ALP median 110 U/L (range 20-593); total",
      "protein median 62.80 g/L (range 46.90-87.80)."
    ),
    co_medication  = paste(
      "Proportions of the concomitant-medication records in Table 1:",
      "non-steroidal anti-inflammatory drugs 1,580 (98.50%), vancomycin",
      "1,436 (89.53%), bleomycin 941 (58.67%), dexamethasone 751 (46.82%),",
      "cisplatin 627 (39.09%), etoposide 12 (0.75%). Only bleomycin was",
      "retained. Cisplatin was dropped because it is 'administered 2-3 days",
      "after MTX, by which time a large portion of MTX has already been",
      "eliminated and many patients lacked MTX concentrations (with records",
      "available for 53.07% at 48 h, 26.53% at 60 h, and only 1.39% at 72 h)'",
      "(Discussion). Dexamethasone, vancomycin and etoposide were tested and",
      "not retained; they have no canonical covariate column and so are not",
      "listed in covariatesDataExcluded."
    ),
    dose_range     = paste(
      "High-dose intravenous methotrexate >= 0.5 g/m^2 (inclusion criterion,",
      "Methods 2.1), standard dose 1.3 g/m^2: one third of the total dose as",
      "a bolus intravenous infusion over 1 h, the remaining two thirds infused",
      "over 11 h. Calcium folinate 13 mg/m^2 rescue started 12 h after MTX",
      "discontinuation, 5 doses every 6 h, adjusted on serum MTX",
      "concentration."
    ),
    regions        = "China (Beijing Puren Hospital, Beijing).",
    notes          = paste(
      "Retrospective therapeutic-drug-monitoring data from patients",
      "hospitalized between February 2015 and July 2018; 5,470 serum MTX",
      "concentrations from 505 patients (median 1.2 umol/L, range",
      "0.01-55.80). Inclusion required a confirmed iGCT diagnosis, an",
      "intravenous MTX dose >= 0.5 g/m^2, and at least one measured MTX",
      "concentration. Sampling (Methods 2.3): fasting venous blood 24 h after",
      "the completion of each MTX infusion, then every 12 h up to 108 h.",
      "Assay: enzyme-multiplied immunoassay (Siemens SYVA Viva-E),",
      "quantitative range 0.3-2,600 umol/L, though 'concentrations below",
      "0.3 umol/L were still detected' -- consistent with the Table 1 minimum",
      "of 0.01 umol/L, so sub-LLOQ values were retained rather than censored.",
      "Estimation: Phoenix NLME 8.3 with FOCE-ELS; evaluation by GOF plots, a",
      "1,000-run bootstrap, VPC and NPDE from 1,000 Monte Carlo simulations.",
      "Baseline demographics from Table 1; the covariate build from Table 2;",
      "final parameter estimates from Table 3; the individual-parameter",
      "covariate model from Eqs. 5-8.",
      "SEX DISCREPANCY: Table 1 reports 'Sex (Male/Female) 357/148' whereas",
      "Results 3.1 reports '505 patients (357 females and 148 males)'. The two",
      "cannot both be right. sex_female_pct above follows the Table 1 column",
      "header (148 female, 29.3%), which is also the direction consistent with",
      "the well-established male predominance of intracranial germ cell",
      "tumours. Sex is not a covariate in the final model, so the ambiguity",
      "has no effect on any equation."
    )
  )

  ini({
    # Structural parameters -- Zhao 2025 Table 3 'Final model' column, applied
    # in the individual-parameter forms printed as Eqs. 5-8:
    #
    #   Eq. 5  CL = 12.88 x (eGFR/102.2)^0.23 x (BW/47)^0.39 x e^BLM
    #                x (TBIL/15.3)^-0.05 x (ALB/40.9)^-0.18   L/h
    #          BLM = 0.08 when combined with BLM, otherwise = 0
    #   Eq. 6  Vc = 72.04 x (BW/47)^0.31                       L
    #   Eq. 7  Q  = 1.08 (fixed)                               L/h
    #   Eq. 8  Vp = 94.94 (fixed)                              L
    lcl <- log(12.88); label("Clearance at the reference covariate set (L/h)")            # Table 3 final CL = 12.88 L/h (%RSE 0.02, 95% CI 12.49-13.26; bootstrap median 12.85, 95% CI 12.41-13.38); Eq. 5
    lvc <- log(72.04); label("Central volume of distribution at 47 kg body weight (L)")   # Table 3 final Vc = 72.04 L (%RSE 0.01, 95% CI 70.00-74.09; bootstrap median 72.03, 95% CI 69.65-74.60); Eq. 6

    # Q and Vp were held at their base-model estimates. Results 3.2: 'Based on
    # the stable estimates from the base model, Q and Vp were fixed at 1.08 L/h
    # and 94.94 L, respectively. This was due to limited data, which caused
    # significant shrinkage in the estimation of their interindividual
    # variability, necessitating their fixation at zero.' Table 3 marks both
    # 'Fixed' in the base, final and bootstrap columns with no CI.
    lq  <- fixed(log(1.08));  label("Intercompartmental clearance (L/h)")                 # Table 3 'Q (L/h) 1.08 Fixed'; Eq. 7
    lvp <- fixed(log(94.94)); label("Peripheral volume of distribution (L)")              # Table 3 'Vp (L) 94.94 Fixed'; Eq. 8

    # Covariate effects. All four continuous terms are power functions of the
    # covariate divided by its cohort median (Methods 2.5.2: 'Continuous
    # covariates were standardized to their respective median values before
    # analysis'); the bleomycin term is exponential.
    e_crcl_cl      <-  0.23; label("Power exponent of BSA-normalized eGFR on CL (unitless)")      # Table 3 'eGFR on CL' = 0.23 (%RSE 0.14, 95% CI 0.17-0.29; bootstrap 0.23, 95% CI 0.17-0.29); Eq. 5
    e_wt_cl        <-  0.39; label("Power exponent of body weight on CL (unitless)")              # Table 3 'BW on CL' = 0.39 (%RSE 0.09, 95% CI 0.32-0.46; bootstrap 0.39, 95% CI 0.32-0.46); Eq. 5
    e_tbili_cl     <- -0.05; label("Power exponent of total bilirubin on CL (unitless)")          # Table 3 'TBIL on CL' = -0.05 (95% CI -0.07 to -0.02; bootstrap -0.05, same CI); Eq. 5
    e_alb_cl       <- -0.18; label("Power exponent of serum albumin on CL (unitless)")            # Table 3 'ALB on CL' = -0.18 (95% CI -0.31 to -0.04; bootstrap -0.18, same CI); Eq. 5
    e_bleomycin_cl <-  0.08; label("Log-scale effect of concomitant bleomycin on CL (unitless)")  # Table 3 'BLM on CL' = 0.08 (%RSE 0.16, 95% CI 0.05-0.10; bootstrap 0.08, 95% CI 0.05-0.10); Eq. 5 'BLM = 0.08 when combined with BLM, otherwise = 0'; exp(0.08) = 1.083, the 1.08-fold increase quoted in the Abstract
    e_wt_vc        <-  0.31; label("Power exponent of body weight on Vc (unitless)")              # Table 3 'BW on V (L)' = 0.31 (%RSE 0.12, 95% CI 0.24-0.39; bootstrap 0.31, 95% CI 0.24-0.39); Eq. 6

    # Between-subject variability, exponential per Methods 2.5.1 Eq. 1
    # (theta = theta_TV x e^eta, eta ~ N(0, omega^2)). Table 3 labels its two
    # variability rows 'IIV Vc(CV%)' and 'IIV CL(CV%)' and prints 0.32 and 2.98
    # for the final model (1.20 and 5.07 for the base model), with no RSE and no
    # CI on either row. The printed SCALE is ambiguous, and Table 3's headers are
    # demonstrably unreliable about scale: its '(%RSE)' column holds FRACTIONS,
    # not percents, on all nine rows that carry a CI (e.g. Vc 72.04 with CI
    # 70.00-74.09 implies SE 1.043 and RSE 0.0145, printed as '0.01'; TBIL -0.05
    # with CI -0.07 to -0.02 implies RSE -0.255, printed as '-0.27'). So the '%'
    # in 'CV%' cannot be taken at face value either.
    #
    # The scale is settled by a variance-decomposition identity internal to the
    # paper. Covariates can only explain variability the BASE model's IIV
    # actually contained, so
    #     omega^2(base) - omega^2(final)  ~=  var(log covariate contribution),
    # and the right-hand side is computable from the Table 1 covariate ranges and
    # the Eq. 5-6 exponents. Treating the Table 1 ranges as +-3.1 SD of a
    # log-normal (n = 505), the retained covariates inject var(log CL) = 0.0286
    # and var(log Vc) = 0.0116. Against that:
    #
    #   reading                     drop in CL    ratio    drop in Vc   ratio
    #   value is CV in percent        0.00168      0.06      0.00013     0.01
    #   value is CV as a fraction     0.99429     34.77      0.79451    68.32
    #   value is omega^2 x 100        0.02090      0.73      0.00880     0.76
    #   value is omega^2              2.09000     73.10      0.88000    75.67
    #
    # Only 'omega^2 x 100' is the right order of magnitude, and it lands on the
    # SAME ratio for CL and Vc simultaneously (0.73 and 0.76) -- two parameters
    # with different covariate structures agreeing on one scale factor, with the
    # shared shortfall from collinearity between BW/BSA/eGFR and from the
    # range-to-SD conversion. The other three readings miss by 17x to 75x. The
    # VPC (Figure 2) corroborates: simulated 5th-95th spans are ~1.2 decades
    # under this reading against ~1.3 decades read off the figure, whereas the
    # 'CV as a fraction' and 'omega^2' readings give >4 decades.
    #
    # Encoded therefore as omega^2 = printed/100. Final-model CVs are
    # 100*sqrt(exp(omega^2) - 1) = 17.4% on CL and 5.7% on Vc; the base-model
    # values 5.07 and 1.20 correspond to 22.8% and 11.0%, and the covariates
    # absorb the difference, as the base-vs-final ordering requires.
    etalcl ~ 0.0298    # Table 3 'IIV CL(CV%)' final = 2.98 -> omega^2 = 0.0298 (CV 17.4%)
    etalvc ~ 0.0032    # Table 3 'IIV Vc(CV%)' final = 0.32 -> omega^2 = 0.0032 (CV 5.7%)

    # Residual error. Results 3.2: 'Random residual variability was optimally
    # characterized using a proportional error model.' Methods 2.5.1 Eq. 3
    # gives that model as C_obs = C_pred x (1 + epsilon), which is exactly
    # nlmixr2's prop() form.
    propSd <- 0.37; label("Proportional residual error (fraction)")   # Table 3 final 'sigma (proportional)' = 0.37 (%RSE 0.02, 95% CI 0.36-0.38; bootstrap 0.37, same CI); Methods Eq. 3
  })

  model({
    # Individual PK parameters, Eqs. 5-8. Four power terms and one exponential
    # term multiply the typical clearance; a single power term scales Vc.
    # Q and Vp carry neither a covariate nor a random effect.
    cl <- exp(lcl + etalcl) *
      (CRCL / 102.2)^e_crcl_cl *
      (WT / 47)^e_wt_cl *
      exp(e_bleomycin_cl * CONMED_BLEOMYCIN) *
      (TBILI / 15.3)^e_tbili_cl *
      (ALB / 40.9)^e_alb_cl
    vc <- exp(lvc + etalvc) * (WT / 47)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition. Methotrexate is given as an intravenous
    # infusion straight into the central compartment, so there is no
    # absorption phase and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Serum methotrexate concentration in umol/L (amounts in umol, volumes
    # in L), matching the assay units used throughout the paper.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
