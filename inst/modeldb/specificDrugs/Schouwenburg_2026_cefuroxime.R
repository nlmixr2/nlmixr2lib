Schouwenburg_2026_cefuroxime <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenous cefuroxime in critically ill term neonates and children",
    "(0-16.8 years) admitted to a paediatric intensive care unit, pooled",
    "from three Dutch studies (EXPAT Kids, POPSICLE, PERFORM). Clearance",
    "and inter-compartmental clearance carry a priori allometric weight",
    "scaling (fixed exponent 0.75, reference 70 kg); central and peripheral",
    "volumes carry fixed linear weight scaling (exponent 1, reference",
    "70 kg). Clearance additionally carries two estimated power terms,",
    "on bedside-Schwartz creatinine clearance (reference 81.3",
    "mL/min/1.73 m2) and on postnatal age (reference 391 days).",
    "Interindividual variability is on clearance only; residual error is",
    "combined proportional (40.1%) and additive (0.545 mg/L)",
    "(Schouwenburg 2026)."
  )
  reference <- paste(
    "Schouwenburg S, Preijers T, Wosten-van Asperen RM, Hartman SJF,",
    "de Wildt SN, de Hoog M, Koch BCP, Abdulla A, Wildschut ED. Low Target",
    "Attainment of Intravenous Cefuroxime in Critically Ill Term Neonates",
    "and Children: A Pooled Population Pharmacokinetics Study.",
    "Clin Pharmacokinet. 2026;65(1):95-105. doi:10.1007/s40262-025-01577-2.",
    "Open-access supplement (Online Resource 1) retrieved from EuropePMC",
    "PMC12783212 and used for the model-development narrative, the assay",
    "limits, and the Table S1/S2 target-attainment values reproduced in",
    "the validation vignette."
  )
  vignette <- "Schouwenburg_2026_cefuroxime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both entries confirmed against Schouwenburg 2026:
  # the model is a two-compartment model with intravenous administration
  # only (no absorption depot -- Methods and Table 1 describe intravenous
  # bolus/short-infusion dosing exclusively), and all reported
  # concentrations are total cefuroxime in plasma (Methods 'Study Design
  # and Participants': 'datasets reporting on cefuroxime plasma
  # concentrations'; supplement 'Methods of quantification': plasma
  # quantified by UPLC-MS/MS).
  compartmentData <- list(
    central     = list(analyte = "cefuroxime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefuroxime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Current (not birth) body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column BW. Enters all four structural parameters as an",
        "a priori allometric power term normalized to 70 kg",
        "(Schouwenburg 2026 Table 3 equations and Table 3 footnote:",
        "'Body weight is scaled to 70 kg'). Both exponents were FIXED a",
        "priori, not estimated: Methods ('Covariate Relationship Analysis')",
        "states 'Refinements applied to the model were a priori allometric",
        "scaling with a fixed exponent (i.e., 0.75 on clearances, 1 on",
        "distribution volume)', and the supplement's 'Model development'",
        "section repeats 'Allometric scaling was applied a priori and",
        "scaled to 70kg'. Table 3 prints 0.75 inline in the CL and Q",
        "equations and 1.00 inline in the V1 and V2 equations, all four",
        "without an RSE or a SIR confidence interval, confirming they were",
        "not estimated.",
        "The 70 kg anchor is an EXTRAPOLATION for this cohort, whose median",
        "weight is 9.0 kg (Table 2). exp(lcl) = 5.29 L/h is therefore a",
        "70 kg adult-equivalent typical value, not a value any patient in",
        "the dataset had; the authors quote it that way throughout ('5.29",
        "L/h/70 kg') and compare it with the 5.15 L/h/70 kg of Gertler 2018",
        "on the same basis. Observed range 2.8-70.0 kg (Table 2), so",
        "simulations outside that range are extrapolation.",
        "Time-varying in principle; the source datasets recorded weight",
        "during therapy."
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = paste(
        "Body-surface-area-normalized creatinine clearance (estimated",
        "glomerular filtration rate), calculated with the bedside Schwartz",
        "equation."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CRCL, defined in the Table 3 abbreviation list as",
        "'estimated creatinine clearance (mL/min/1.73 m2)'. Methods",
        "('Covariate Relationship Analysis'): 'creatinine clearance (CRCL;",
        "estimated glomerular filtration rate [eGFR] in mL/min/1.73 m2)",
        "... CRCL was calculated using the bedside Schwartz equation'.",
        "Enters clearance as an estimated power term normalized to 81.3",
        "mL/min/1.73 m2, the dataset median (Table 3 footnote: 'Creatinine",
        "clearance and postnatal age are scaled to the dataset median",
        "(81.3 mL/min/1.73 m2 and 391 days)'; Table 2 confirms the pooled",
        "median creatinine clearance is 81.3).",
        "Observed range 25.4-181.5 mL/min/1.73 m2 (Table 2). The paper's",
        "own simulations stratify this covariate into four bands (< 30,",
        "30-80, 80-120, > 120 mL/min/1.73 m2), and its augmented-renal-",
        "clearance cutoffs are 99 mL/min/1.73 m2 below 2 years of age and",
        "140 above; 8.9% of the cohort met those cutoffs.",
        "Time-varying: creatinine was measured repeatedly during therapy.",
        "The power form is singular at CRCL = 0 (clearance goes to zero",
        "with a positive exponent), which is outside the observed range and",
        "outside any physiologically meaningful simulation."
      ),
      source_name        = "CRCL"
    ),
    PNA = list(
      description        = "Postnatal (chronological) age since birth, during therapy.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Schouwenburg 2026 reports postnatal age in DAYS (Table 3",
        "abbreviation list: 'PNA, postnatal age (days)'; the Table 3",
        "clearance equation uses (PNA / 391) with PNA in days; the 391-day",
        "reference is the dataset median per the Table 3 footnote",
        "'Creatinine clearance and postnatal age are scaled to the dataset",
        "median (81.3 mL/min/1.73 m2 and 391 days)', which matches the",
        "Table 2 pooled median postnatal age of 391 days). The canonical",
        "nlmixr2lib PNA column carries MONTHS, so model() reparameterises",
        "the reference as 391 / 30.4375 = 12.846 months; numerator and",
        "denominator carry the same units factor, so the ratio and hence",
        "the exponent are unchanged. This follows the",
        "Schouwenburg_2025_clavulanicAcid, Zhao_2018_omeprazole and",
        "Bardhi_2026_ampicillin_foal precedents. Users supply PNA in",
        "months.",
        "NOTE 1: the Results text describes the retained PNA relationship",
        "as 'an exponential function' ('PNA was found to best describe",
        "maturation of CLcefu as an exponential function'), but the Table 3",
        "equation prints a POWER form, (PNA / 391)^theta_PNA, and the",
        "Methods state that continuous covariates 'were centered on the",
        "median and were evaluated as exponential or power relationships'.",
        "The printed equation is authoritative, so the power form is",
        "encoded. The identical text-versus-equation conflict occurs in the",
        "same group's Schouwenburg_2025_clavulanicAcid and was resolved the",
        "same way.",
        "NOTE 2: the power form is singular at PNA = 0 -- clearance goes to",
        "zero on the day of birth and concentration diverges. Table 2",
        "reports a postnatal-age range of 0-6131 days, so the dataset does",
        "contain day-of-birth records, and the published model cannot be",
        "evaluated at exactly PNA = 0 as written. Use PNA >= 1 day",
        "(0.0329 months) when simulating. The exponent is small (0.0869),",
        "so away from zero the age effect is weak: a 7-day-old is predicted",
        "to have only 30% lower clearance than the 391-day reference at the",
        "same weight and renal function.",
        "NOTE 3: no premature neonates were enrolled (Discussion: 'our",
        "study population does not contain premature neonates'), which is",
        "why a postmenstrual-age Hill maturation function was tested and",
        "rejected in favour of plain postnatal age. The model is therefore",
        "not intended for preterm extrapolation."
      ),
      source_name        = "PNA"
    )
  )

  # Covariates that Schouwenburg 2026 screened on clearance during the
  # stepwise forward-inclusion (P < 0.05) / backward-elimination (P < 0.01)
  # covariate analysis but did NOT retain in the final model (Methods,
  # 'Covariate Relationship Analysis': "The following covariates were
  # tested: sex, postnatal age (PNA), gestational age (GA), postmenstrual
  # age (PMA), serum creatinine (SCR, umol/L; age adjusted), creatinine
  # clearance ..., serum urea (mmol/L), C-reactive protein (mg/L), and
  # study center"). Only PNA and CRCL survived. Documentation only -- no
  # published point estimates exist for these, so they carry no effect in
  # model(). Study centre is a study-design label rather than a patient
  # characteristic and is described in population$notes instead.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Tested on clearance as a categorical covariate (Methods:",
        "'Categorical variables were modeled using a proportional model')",
        "and not retained. Supplementary Fig. S4 plots ETACL against sex",
        "as a boxplot and shows no trend. Pooled cohort was 55.6% female",
        "(Table 2). No point estimate is published, so no effect is",
        "encoded."
      ),
      source_name        = "SEX"
    ),
    GA = list(
      description        = "Gestational age at birth. Time-fixed per subject.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance and not retained; postnatal age alone best",
        "described maturation (Results: 'PNA was found to best describe",
        "maturation of CLcefu as an exponential function'). Pooled median",
        "GA 39.3 weeks, range 26.0-42.0 (Table 2). No point estimate is",
        "published, so no effect is encoded."
      ),
      source_name        = "GA"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age plus postnatal age).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance and not retained. The supplement's 'Model",
        "development' section records that three postmenstrual-age",
        "maturation structures were evaluated -- 'PMA hill, PMA hill (hill",
        "1 fix), PMA hill estimates Gertler et al.' -- and the Discussion",
        "explains the rejection: 'Both a PMA hill equation and PNA were",
        "tested as covariates, whereas a sigmoid function did not improve",
        "model fit or alter population estimates compared with PNA.",
        "Although a PMA hill equation is often applied to describe prenatal",
        "kidney function maturation, our study population does not contain",
        "premature neonates.' No point estimate is published for any of the",
        "PMA structures, so no effect is encoded."
      ),
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine, age-adjusted.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance and not retained; the derived, BSA-normalized",
        "CRCL was retained instead. Methods list 'serum creatinine (SCR,",
        "umol/L; age adjusted)' among the tested covariates, and the",
        "supplement adds that 'A function for age adjusted mean creatinine",
        "was tested on cefuroxime clearance' referencing the Ceriotti 2008",
        "paediatric reference intervals. Pooled median 27 umol/L, range",
        "14-219 (Table 2). No point estimate is published, so no effect is",
        "encoded."
      ),
      source_name        = "SCR"
    ),
    BUN = list(
      description        = paste(
        "Serum urea. Reported as urea (not urea nitrogen) in mmol/L; to",
        "convert to blood-urea-nitrogen mg/dL multiply by 2.80."
      ),
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance and not retained (Methods: 'serum urea",
        "(mmol/L)'). Pooled median 3.5 mmol/L, range 1.4-36.9 (Table 2).",
        "No point estimate is published, so no effect is encoded."
      ),
      source_name        = "UREA"
    ),
    CRP = list(
      description        = "C-reactive protein, an acute-phase inflammatory marker.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested on clearance and not retained (Methods: 'C-reactive",
        "protein (mg/L)'). Pooled median 54.0 mg/L, range 1.3-375.0",
        "(Table 2). No point estimate is published, so no effect is",
        "encoded."
      ),
      source_name        = "CRP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 45,
    n_studies      = 3,
    n_observations = 148,
    age_range      = "0.00-16.80 years (median 1.07 years)",
    pna_range      = "postnatal age 0-6131 days (median 391 days)",
    ga_range       = "gestational age 26.0-42.0 weeks (median 39.3 weeks)",
    weight_range   = "2.8-70.0 kg",
    weight_median  = "9.0 kg",
    height_range   = "47-192 cm (median 68 cm)",
    sex_female_pct = 55.6,
    race_ethnicity = NULL,
    renal_function = paste(
      "bedside-Schwartz creatinine clearance 25.4-181.5 mL/min/1.73 m2",
      "(median 81.3); serum creatinine 14-219 umol/L (median 27);",
      "8.9% (4/45) met the age-adjusted augmented-renal-clearance cutoffs",
      "of 99 mL/min/1.73 m2 below 2 years and 140 above"
    ),
    disease_state  = paste(
      "critically ill term neonates and children admitted to a level 3",
      "paediatric or cardiac intensive care unit and treated with",
      "intravenous cefuroxime for suspected or confirmed bacterial",
      "infection; no premature neonates and no patients on extracorporeal",
      "life support were included"
    ),
    dose_range     = paste(
      "intravenous cefuroxime 65.0-1500.0 mg per administration",
      "(median 225.0 mg; median 25.0 mg/kg/administration), given per the",
      "Dutch Children's Formulary (Kinderformularium) as 70 mg/kg/day q8h",
      "below 1 month of age and 100 mg/kg/day q6h above 1 month"
    ),
    regions        = "The Netherlands (Rotterdam, Utrecht, Nijmegen)",
    notes          = paste(
      "Pooled analysis of two datasets covering three studies",
      "(Schouwenburg 2026 Results and Table 2): (1) EXPAT Kids, a",
      "multicentre observational PK/PD trial of beta-lactam antibiotics",
      "at Erasmus MC Sophia Children's Hospital (Rotterdam) and Wilhelmina",
      "Children's Hospital UMCU (Utrecht), NL9326, n = 31 with 120",
      "concentrations; and (2) the POPSICLE (Radboudumc Nijmegen,",
      "NCT03248349) and PERFORM (Erasmus MC Rotterdam, NCT03502993)",
      "studies, n = 14 combined (PERFORM n = 10, POPSICLE n = 4) with 28",
      "concentrations. Median 3 samples per individual (range 1-8).",
      "Race/ethnicity was not reported.",
      "STUDY CENTRE was screened as a covariate on clearance and not",
      "retained; the supplement also records that 'Separate error models",
      "per study center were tested to account for assay variability and",
      "study procedures' and were likewise not adopted. Supplementary",
      "Fig. S4 shows ETACL by study centre with no trend, and the Fig. 2",
      "visual predictive check is stratified by study.",
      "Modelling was in NONMEM (the Abstract and Methods disagree on the",
      "version -- the Abstract says 7.5, the Methods say 7.4 -- which does",
      "not affect the reported parameter values) with PsN 4.2.0, Pirana",
      "3.0.0 and R 4.2.2.",
      "Concentrations were quantified by a validated multi-analyte",
      "UPLC-MS/MS assay at Erasmus MC with an LLOQ of 1.25 mg/L and a ULOQ",
      "of 50 mg/L (supplement, 'Methods of quantification'); observed",
      "concentrations spanned 0.6-292.0 mg/L. 12 of 148 samples (8.1%)",
      "were below the LLOQ, were RETAINED rather than excluded, and were",
      "handled with Beal's M1 method; the supplement notes 'Different LLOQ",
      "handling methods did not show difference in model fit, leading to",
      "the implementation of M1'. M1 substitutes nothing and simply keeps",
      "the values, so it has no simulation-time counterpart and is not",
      "encoded here.",
      "Free (unbound) concentrations were NOT measured: only total",
      "cefuroxime was available and no protein-binding correction was",
      "applied (Methods: 'Cefuroxime exhibits moderate protein binding up",
      "to 30% in critically ill adults. However, the effects of protein",
      "binding could not be evaluated, as only total concentrations for",
      "cefuroxime were available; owing to limited protein binding, no",
      "correction was performed'). The paper's %T>MIC target attainment is",
      "therefore computed on TOTAL concentrations, which the Discussion",
      "flags as possibly optimistic in hypoalbuminaemic patients.",
      "A one-compartment model was fitted first and rejected in favour of",
      "the two-compartment model (supplement: dOFV = -56.75).",
      "Sparse sampling precluded estimation of interindividual variability",
      "on V1, V2 or Q (Discussion: 'Sparse sampling resulted in the",
      "inability to accurately estimate an IIV for V1, V2, or Q'), so only",
      "clearance carries an eta. No correlation between etas is estimable",
      "with a single eta.",
      "Model evaluation was internal only: goodness-of-fit plots,",
      "visual predictive checks stratified by study, conditional weighted",
      "residuals, and sampling importance resampling (SIR). No external",
      "validation dataset was available. SIR relative standard errors were",
      "below 30% for all structural parameters."
    )
  )

  ini({
    # --- Structural parameters ---
    # Reference subject for every typical value below: WT 70 kg,
    # CRCL 81.3 mL/min/1.73 m2, PNA 391 days -- i.e. all three covariate
    # ratios in the Table 3 equations equal 1.
    lcl <- log(5.29)
    label("Typical clearance TVCL at WT 70 kg, CRCL 81.3, PNA 391 d (L/h)")
    # Schouwenburg 2026 Table 3, 'TVCL (L/h/70kg)': 5.29 (RSE 6.78%);
    # SIR 5.30 (95% CI 4.76-5.92) [RSE 6.6%]. The Abstract and Results
    # both restate it: 'Intravenous cefuroxime clearance was estimated at
    # 5.29 L/h/70 kg', and the Discussion compares it with Gertler 2018's
    # 5.15 L/h/70 kg.

    lvc <- log(5.02)
    label("Typical central volume TVV1 at WT 70 kg (L)")
    # Schouwenburg 2026 Table 3, 'TVV1 (L)': 5.02 (RSE 20.1%);
    # SIR 4.98 (95% CI 3.42-5.92) [RSE 19.7%].

    lvp <- log(12.3)
    label("Typical peripheral volume TVV2 at WT 70 kg (L)")
    # Schouwenburg 2026 Table 3, 'TVV2 (L)': 12.3 (RSE 9.44%);
    # SIR 12.41 (95% CI 10.51-14.28) [RSE 9.3%].

    lq <- log(27.9)
    label("Typical inter-compartmental clearance TVQ at WT 70 kg (L/h)")
    # Schouwenburg 2026 Table 3, 'TVQ (L/h)': 27.9 (RSE 32.9%);
    # SIR 28.23 (95% CI 18.26-44.06) [RSE 29.3%]. This is the least
    # precisely estimated structural parameter, consistent with the
    # Discussion's note that sparse sampling limited the distribution
    # phase.

    # --- Allometric exponents: FIXED a priori, not estimated ---
    # Methods, 'Covariate Relationship Analysis': "Refinements applied to
    # the model were a priori allometric scaling with a fixed exponent
    # (i.e., 0.75 on clearances, 1 on distribution volume)". The
    # supplement's 'Model development' section repeats "Allometric scaling
    # was applied a priori and scaled to 70kg". Table 3 prints all four
    # exponents inline in the parameter equations with no RSE and no SIR
    # confidence interval. One exponent serves both clearances and one
    # serves both volumes, exactly as the single Methods sentence states.
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of body weight on CL and Q (unitless)")
    # Schouwenburg 2026 Table 3 equations:
    #   CLcefu = TVCL * (BW/70)^0.75 * (CRCL/81.3)^theta_CRCL * (PNA/391)^theta_PNA
    #   Qcefu  = TVQ  * (BW/70)^0.75

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of body weight on V1 and V2 (unitless)")
    # Schouwenburg 2026 Table 3 equations:
    #   V1cefu = TVV1 * (BW/70)^1.00
    #   V2cefu = TVV2 * (BW/70)^1.00

    # --- Estimated covariate effects on clearance ---
    e_crcl_cl <- 0.754
    label("Creatinine-clearance power exponent on clearance (unitless; reference 81.3 mL/min/1.73 m2)")
    # Schouwenburg 2026 Table 3, 'theta_CLCR' (under 'Covariate
    # relationships'): 0.754 (RSE 22%); SIR 0.766 (95% CI 0.49-1.02)
    # [RSE 21.6%]. Retained by stepwise forward inclusion (P < 0.05) and
    # backward elimination (P < 0.01).

    e_pna_cl <- 0.0869
    label("Postnatal-age power exponent on clearance (unitless; reference 391 days)")
    # Schouwenburg 2026 Table 3, 'theta_PNA' (under 'Covariate
    # relationships'): 0.0869 (RSE 30.9%); SIR 0.084 (95% CI 0.043-0.130)
    # [RSE 30.0%]. The Results text calls this 'an exponential function'
    # while Table 3 prints the POWER form (PNA/391)^theta_PNA; the printed
    # equation is authoritative (see covariateData$PNA$notes NOTE 1).

    # --- Interindividual variability (clearance only) ---
    etalcl ~ 0.110889
    # Schouwenburg 2026 Table 3, 'Interindividual variability (IIV) /
    # Clearance (%CV)': 33.3 (RSE 26.1%) [shrinkage 10.5%]; SIR 33.8
    # (95% CI 27.9-40.9) [RSE 23.4%]. The SIR interval is reported on the
    # %CV scale (27.9-40.9 brackets the 33.8 %CV point estimate), so the
    # row is unambiguously a CV percentage and not a variance.
    # omega^2 = 0.333^2 = 0.110889, the direct-square convention used by
    # the same author group in Schouwenburg_2025_clavulanicAcid. The
    # strict log-normal conversion log(CV^2 + 1) = 0.10517 (32.4 %CV) is
    # an equally defensible reading; the difference is 3% on the SD scale
    # and immaterial.
    # Sparse sampling precluded IIV on V1, V2 and Q (Discussion), so no
    # eta is encoded on lvc, lvp or lq.

    # --- Residual unexplained variability (combined) ---
    propSd <- 0.401
    label("Proportional residual error (fraction)")
    # Schouwenburg 2026 Table 3, 'Residual variability / Proportional
    # error (%)': 40.1 (RSE 7.9%); SIR 40.6% (95% CI 35.8%-46.1%)
    # [RSE 7.9%]. Results: 'A mixed-error model was used to describe the
    # residual variability'; the supplement's 'Model development' section
    # states 'Residual variability was evaluated using a combined or
    # separate (proportional or additive) error model'.

    addSd <- 0.545
    label("Additive residual error (mg/L)")
    # Schouwenburg 2026 Table 3, 'Residual variability / Additive error
    # (mg/L)': 0.545 (RSE 25.6%); SIR 0.550 (95% CI 0.385-0.796)
    # [RSE 23.3%]. Comparable to the 1.25 mg/L assay LLOQ (supplement,
    # 'Methods of quantification'), as expected for an additive term that
    # absorbs assay noise near the lower limit.
  })

  model({
    # 1. Derived covariate terms
    #    Schouwenburg 2026 Table 3:
    #      CLcefu (L/h) = TVCL * (BW/70)^0.75
    #                          * (CRCL/81.3)^theta_CRCL
    #                          * (PNA/391)^theta_PNA
    #      V1cefu (L)   = TVV1 * (BW/70)^1.00
    #      V2cefu (L)   = TVV2 * (BW/70)^1.00
    #      Qcefu  (L/h) = TVQ  * (BW/70)^0.75
    #    The PNA reference is 391 DAYS in the paper; the canonical PNA
    #    column is in MONTHS, so the reference is converted once here.
    #    The ratio is unit-invariant, so the exponent transfers unchanged.
    pna_ref_months <- 391 / 30.4375
    f_pna <- (PNA / pna_ref_months)^e_pna_cl
    f_crcl <- (CRCL / 81.3)^e_crcl_cl
    allom_cl <- (WT / 70)^e_wt_cl_q
    allom_v <- (WT / 70)^e_wt_vc_vp

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * allom_cl * f_crcl * f_pna
    vc <- exp(lvc) * allom_v
    vp <- exp(lvp) * allom_v
    q <- exp(lq) * allom_cl

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system (two compartments, intravenous administration only --
    #    every dose in the pooled dataset was intravenous, so there is no
    #    absorption depot and no bioavailability term)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error (combined proportional plus additive,
    #    Table 3 'Residual variability')
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
