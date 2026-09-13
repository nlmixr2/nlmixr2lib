Liu_2025_selpercatinib <- function() {
  description <- paste(
    "Two-compartment population PK model for selpercatinib (RETEVMO, a",
    "first-in-class highly selective RET kinase inhibitor approved for",
    "RET-altered lung, thyroid and other solid tumors) in adult,",
    "adolescent and pediatric patients (Liu 2025; N = 830 patients,",
    "8024 plasma concentrations pooled from the phase 1/2 studies",
    "LIBRETTO-001 in patients aged 12 years and older and LIBRETTO-121",
    "in patients aged 6 months to 21 years). Absorption is sequential",
    "zero-order then first-order: the oral dose enters the gut depot over",
    "a zero-order window Dur = 1.09 h and the depot then drains",
    "first-order at ka = 1.47 1/h. Disposition is two-compartment with",
    "first-order elimination; typical apparent values for a 70 kg patient",
    "at the 160 mg reference dose are CL/F = 6.04 L/h, Vc/F = 99.6 L,",
    "Q/F = 29.6 L/h and Vp/F = 91.3 L. Relative bioavailability F1 is",
    "fixed to 1 and carries the only non-weight covariate on absorption:",
    "Asian race raises F1 by 18.3%. Body weight enters allometrically",
    "with the standard fixed exponents referenced to 70 kg (0.75 on CL/F",
    "and Q/F, 1 on Vc/F and Vp/F). The administered dose acts on CL/F as",
    "a time-varying linear covariate centered on 160 mg,",
    "-0.321% per mg, so apparent clearance falls as the dose rises and",
    "exposure is more than dose proportional over the 20-240 mg range",
    "studied. Baseline age, sex, creatinine clearance, liver function",
    "tests and concomitant medication were screened and not retained; an",
    "age effect on CL/F was statistically significant but was removed",
    "during model refinement because it biased predictions in children.",
    "Inter-individual variability is estimated on CL/F (48.8%),",
    "Vc/F (66.1%), ka (63.6%) and Dur (56.2%); residual error is",
    "proportional (25.3%) plus additive (61.5 ng/mL).",
    sep = " "
  )
  reference <- paste(
    "Liu D, van der Walt JS. (2025).",
    "Population pharmacokinetics modeling of selpercatinib to support",
    "posology in pediatric patients with RET-altered metastatic thyroid",
    "cancer or solid tumors.",
    "CPT Pharmacometrics Syst Pharmacol 14(11):1848-1857.",
    "doi:10.1002/psp4.70042",
    sep = " "
  )
  vignette <- "Liu_2025_selpercatinib"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Liu 2025 Figure 1 (schematic view of
  # the final population pharmacokinetic model), whose boxes are 'Gut',
  # 'Central Vc/F' and 'Peripheral Vp/F', and whose caption names the
  # absorption quantities 'Dur, duration of zero-order absorption' and
  # 'ka, first-order absorption rate constant'.
  compartmentData <- list(
    depot       = list(analyte = "selpercatinib", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "selpercatinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "selpercatinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline, time-fixed. Population median 67.0 kg (range",
        "9.6-179; Liu 2025 Table 1 'Weight, kg', Total column), with the",
        "pediatric LIBRETTO-121 cohort at median 48.5 kg (range 9.6-97.7)",
        "and the mostly-adult LIBRETTO-001 cohort at median 67.3 kg",
        "(range 26.8-179).",
        "The allometric reference weight is 70 kg, NOT the cohort median:",
        "Liu 2025 Methods 'Population PK Analysis' states 'body weight was",
        "normalized to 70 kg with exponents fixed at 0.75 (CL/F and Q/F)",
        "or 1 (Vc/F and Vp/F)', and the Table 2 covariate-function block",
        "prints all four terms with an explicit 70 kg denominator.",
        "Weight was carried over from the previously developed adult model",
        "rather than selected by the stepwise covariate search -- Results",
        "'PK Analyses' says 'Body weight effects were predefined and",
        "included in the model on clearances and volumes of distribution'.",
        "Table S2 shows weight was nonetheless the one covariate evaluated",
        "on all four disposition parameters (CL/F, Vc/F, Vp/F, Q/F)."
      ),
      source_name        = "Weight at baseline"
    ),
    DOSE = list(
      description        = "Administered selpercatinib dose per administration",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Register use case (b), time-varying: the current administered",
        "dose amount at the time of the record, in mg per administration",
        "(NOT per day). Liu 2025 Methods 'Population PK Analysis': 'The",
        "relationship between dose and CL/F was modeled as a time-varying",
        "continuous covariate to account for any dose reductions or",
        "increases.'",
        "The centering value is 160 mg, the adult recommended phase 2",
        "dose: the Table 2 covariate function is",
        "CL = theta_CL * (1 + theta_dose * (Dose - 160 mg)) *",
        "(Weight/70 kg)^0.75, and Table S2 footnote (b) states 'Effect of",
        "dose is relative to a 160 mg dose of selpercatinib'. It is a",
        "per-administration amount, not a daily amount, because 160 mg is",
        "the twice-daily unit dose and the recommended pediatric regimens",
        "in Table 3 are likewise stated per administration (40 mg TID,",
        "80/120/160 mg BID).",
        "Observed starting-dose levels span 20 mg QD/BID to 240 mg BID",
        "(Table S3), with 86.5% of patients starting at 160 mg BID.",
        "The effect is LINEAR, not power-form, so it extrapolates poorly:",
        "the multiplier (1 - 0.00321 * (DOSE - 160)) reaches zero at",
        "DOSE = 471.5 mg and turns negative above it. No clamp is imposed",
        "here because the paper imposes none; keep DOSE within the",
        "20-240 mg range the model was fit to."
      ),
      source_name        = "Selpercatinib dose"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian)",
      notes              = paste(
        "Per-subject, time-fixed. Liu 2025 Table S2 footnote (a) defines",
        "the contrast as 'Differences between Asian (Japanese, Chinese,",
        "East-Asian) and non-Asian subjects', i.e. the pooled Asian group",
        "against everyone else, which is exactly the canonical",
        "RACE_ASIAN orientation (1 = Asian, 0 = other).",
        "The retained effect is on relative bioavailability, not on",
        "clearance: Table 2 'Asian race effect on F1 (%)' = 18.3%",
        "(95% CI 10.9/25.8) entering the fractional covariate function",
        "F1 = theta_F1 * (1 + theta_Asian race effect on F1), so Asian",
        "patients have F1 = 1.183 against the reference F1 = 1.00.",
        "Because every disposition parameter in Table 2 is an APPARENT",
        "(/F) value, an 18.3% higher F1 raises exposure by 18.3% at a",
        "given dose. Table S2 shows race was evaluated on both CL/F and F;",
        "the search retained it on F only.",
        "Liu 2025 does not tabulate the race distribution of the analysis",
        "population -- Table 1 reports only sex, age, weight and BSA -- so",
        "the prevalence used in the vignette cohort is an assumption; see",
        "the vignette Errata."
      ),
      source_name        = "Race"
    )
  )

  # Screened in the Liu 2025 stepwise covariate search (Table S2) but NOT
  # retained in the final model, so no coefficient exists to encode. Recorded
  # here to preserve the provenance of the covariate screen without carrying
  # 'declared but not referenced' convention warnings.
  #
  # The age entry is the interesting one: age on CL/F WAS statistically
  # significant and survived backward deletion, then was deliberately removed
  # during model refinement. Liu 2025 Results 'PK Analyses' and Discussion
  # give the numbers and the reasoning, reproduced in that entry's notes.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F (Table S2) and RETAINED by the automated",
        "stepwise covariate search, then removed during model refinement.",
        "Liu 2025 Results 'PK Analyses': 'The inclusion of an age effect",
        "in CL/F resulted in an increase in the typical value of CL/F from",
        "5.94 (95% CI: 5.73-6.15) L/h in the base model to 7.30 (95% CI:",
        "6.73-7.88) L/h in the SCM final backward model. Due to the faster",
        "clearance introduced by the age effect, the model underpredicted",
        "the selpercatinib concentration-time profiles during C1D8 in",
        "patients in study LIBRETTO-121.'",
        "The Discussion reassessed it two ways against the final model --",
        "a linear 'proportional change model' centered on the median age",
        "of 58 years and a segmented 'hockey stick' model with the break",
        "point at 58 years -- which dropped the objective function value",
        "by 43.998 and 49.995 points respectively, yet both 'resulted in",
        "biased population predictions in children'. Neither coefficient",
        "is reported numerically, so neither is encodable even as an",
        "option; the age effect is absent from the final model.",
        "Population median age 58.0 years (range 2.0-92.0), Table 1."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F, Vc/F, Vp/F and Q/F (Table S2); not retained.",
        "403/830 patients (48.6%) were female (Table 1)."
      )
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened on CL/F (Table S2); not retained."
    ),
    ALT = list(
      description = "Baseline alanine transaminase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F as one of the three 'liver function tests' of",
        "Table S2 footnote (c); not retained."
      )
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F as one of the three 'liver function tests' of",
        "Table S2 footnote (c); not retained."
      )
    ),
    BILI = list(
      description = "Baseline total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F as one of the three 'liver function tests' of",
        "Table S2 footnote (c); not retained."
      )
    ),
    CONMED_CYP3A4_INHIBITOR = list(
      description = "Concomitant CYP3A4 inhibitor use during the PK sampling period",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F and F as 'Concomitant medication' (Table S2);",
        "not retained. 274/830 patients (33.0%) used any CYP3A4 inhibitor,",
        "29 (3.5%) a strong one (Table S1). Use was confined to",
        "LIBRETTO-001: no LIBRETTO-121 patient used one, so the pediatric",
        "cohort contributes no information on this covariate",
        "(Methods 'Population PK Analysis')."
      )
    ),
    CONMED_CYP3A4_INDUCER = list(
      description = "Concomitant CYP3A4 inducer use during the PK sampling period",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F and F as 'Concomitant medication' (Table S2);",
        "not retained. Only 6/830 patients (0.4%) used any CYP3A4 inducer",
        "(Table S1), all in LIBRETTO-001."
      )
    ),
    CONMED_ANTACID = list(
      description = "Concomitant acid-reducing agent (proton pump inhibitor or H2 receptor antagonist) use",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on CL/F and F as 'Concomitant medication' (Table S2);",
        "not retained. 367/830 patients (44.2%) used a PPI or an",
        "H2-blocker (Table S1); only 1 of the 27 LIBRETTO-121 patients did",
        "(an H2-blocker)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 830,
    n_observations = 8024,
    n_studies      = 2,
    studies        = paste(
      "LIBRETTO-001 (NCT03157128), phase 1/2 in patients aged 12 years",
      "and older with advanced or metastatic solid tumors: 803 patients,",
      "7723 concentrations. LIBRETTO-121 (NCT03899792), phase 1/2 in",
      "patients aged 6 months to 21 years with an activating RET",
      "alteration and an advanced solid or primary CNS tumor: 27 patients,",
      "301 concentrations. Data cut-off 13 January 2023 for both.",
      sep = " "
    ),
    age_range      = "2-92 years (median 58.0); LIBRETTO-121 median 14.0, range 2.0-20.0",
    age_groups     = paste(
      "6 patients (0.7%) under 12 years, 18 (2.2%) 12 to under 18 years,",
      "806 (97.1%) 18 years and older. The 6 patients under 12 years",
      "contributed only 59 concentrations, which the Discussion flags as",
      "the main limitation of the analysis.",
      sep = " "
    ),
    weight_range   = "9.6-179 kg (median 67.0)",
    bsa_range      = "0.446-2.79 m^2 (median 1.76), Mosteller formula; missing for 22/803 LIBRETTO-001 patients",
    sex_female_pct = 48.6,
    race_ethnicity = paste(
      "Not tabulated in Liu 2025 Table 1. Race entered the covariate",
      "search only as pooled Asian (Japanese, Chinese, East-Asian) versus",
      "non-Asian (Table S2 footnote a), and the Asian effect on F1 was",
      "retained in the final model, but neither the Asian fraction of the",
      "analysis population nor any finer race breakdown is reported.",
      sep = " "
    ),
    disease_state  = paste(
      "RET-altered advanced or metastatic solid tumors, including",
      "medullary thyroid cancer, RET fusion-positive thyroid cancer,",
      "RET fusion-positive non-small-cell lung cancer and primary central",
      "nervous system tumors.",
      sep = " "
    ),
    dose_range     = paste(
      "Oral, continuous 28-day cycles. LIBRETTO-001: 20 mg once daily",
      "through 240 mg twice daily in phase 1 dose escalation, 160 mg twice",
      "daily in phase 2. LIBRETTO-121: body-surface-area-based dosing",
      "starting at the exposure equivalent of adult 160 mg twice daily",
      "(92 mg/m^2 twice daily, capped at 160 mg twice daily). Starting-dose",
      "strata are tabulated in Table S3; 718/830 patients (86.5%) started",
      "at 160 mg twice daily.",
      sep = " "
    ),
    notes          = paste(
      "Bioanalytical lower limit of quantification 1 ng/mL; 133 of the",
      "8157 collected concentrations were below it and were excluded from",
      "the analysis. Missing dosing times were imputed for 2 patients. No",
      "outliers were excluded. Estimation was FOCE with interaction in",
      "NONMEM 7.5.0 with PsN 5.3.0. The model is an update of a previously",
      "developed adult model (reference 10 of Liu 2025), from which the",
      "structural form and the fixed allometric exponents were carried",
      "over; all final estimates encoded here are from this paper's own",
      "Table 2, so no upstream source is required.",
      "Simulations in the paper drew pediatric covariates from the NHANES",
      "DXA database. The authors caution that predictions below 2 years of",
      "age may OVERPREDICT concentrations because selpercatinib is",
      "predominantly metabolized by CYP3A4, which matures around 2 years,",
      "and no patient under 2 years was in the dataset; the model carries",
      "no maturation term. Use it at 2 years and above.",
      sep = " "
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Disposition. Liu 2025 Table 2 'NONMEM parameter estimates for the final
    # model'. Every value is APPARENT (/F) because F1 is fixed to 1 and the
    # drug is dosed orally only: Table 2 lists 'F1 (fraction) 1.00 Fixed' and
    # the Figure 1 caption reads 'F, bioavailability (fixed to 100%)'. The
    # column labels in Table 2 drop the '/F' but Results 'PK Analyses' and
    # the Figure 1 boxes restore it (CL/F, Vc/F, Vp/F, Q/F).
    #
    # These are the values AT the 160 mg reference dose and AT 70 kg: the
    # Table 2 covariate functions collapse to 1 there, since
    # (1 + theta_dose * (160 - 160)) = 1 and (70/70)^p = 1.
    # -----------------------------------------------------------------------
    lcl <- log(6.04)  ; label("Apparent clearance CL/F at 160 mg and 70 kg (L/h)")                        # Liu 2025 Table 2 'CL (L/h)' = 6.04 (95% CI 5.79/6.29)
    lvc <- log(99.6)  ; label("Apparent central volume of distribution Vc/F at 70 kg (L)")                # Liu 2025 Table 2 'Vc (L)' = 99.6 (95% CI 86.3/113)
    lvp <- log(91.3)  ; label("Apparent peripheral volume of distribution Vp/F at 70 kg (L)")             # Liu 2025 Table 2 'Vp (L)' = 91.3 (95% CI 82.9/99.7)
    lq  <- log(29.6)  ; label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)")                 # Liu 2025 Table 2 'Q (L/h)' = 29.6 (95% CI 24.2/35.0)

    # -----------------------------------------------------------------------
    # Absorption. Sequential zero-order then first-order (Abstract, Results
    # 'PK Analyses', Figure 1): the dose is delivered into the gut depot at a
    # constant rate over the window Dur, and the depot then drains into the
    # central compartment first-order at ka.
    #
    # F1 is an estimated-then-fixed anchor rather than an identifiable
    # quantity: with oral-only data every disposition parameter is apparent,
    # so F1 sets the scale and is fixed at 1 ('1.00 Fixed', Table 2). It is
    # kept as a named parameter because the Asian-race covariate acts on it.
    # -----------------------------------------------------------------------
    lka     <- log(1.47)   ; label("First-order absorption rate constant from the depot (1/h)")           # Liu 2025 Table 2 'ka (1/h)' = 1.47 (95% CI 1.25/1.69)
    ld1     <- log(1.09)   ; label("Duration of zero-order absorption into the depot (h)")                # Liu 2025 Table 2 'Dur (h)' = 1.09 (95% CI 1.05/1.13)
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 in non-Asian patients (fraction)")       # Liu 2025 Table 2 'F1 (fraction) 1.00 Fixed'; Figure 1 caption 'F, bioavailability (fixed to 100%)'

    # -----------------------------------------------------------------------
    # Allometry. Both exponents are FIXED at the standard empirical values,
    # not estimated. Liu 2025 Methods 'Population PK Analysis' states it
    # twice: 'Body weight was included following allometric principles with
    # fixed exponents of 0.75 for CL and intercompartmental clearance (Q) and
    # 1.0 for Vc and peripheral volume of distribution (Vp)', and 'body
    # weight was normalized to 70 kg with exponents fixed at 0.75 (CL/F and
    # Q/F) or 1 (Vc/F and Vp/F)'. Neither appears as an estimated row in
    # Table 2 and neither carries a confidence interval, which is the other
    # fixed-parameter signal.
    # -----------------------------------------------------------------------
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F, reference 70 kg (unitless)")     # Liu 2025 Methods 'Population PK Analysis'; Table 2 covariate functions
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on Vc/F and Vp/F, reference 70 kg (unitless)")    # Liu 2025 Methods 'Population PK Analysis'; Table 2 covariate functions

    # -----------------------------------------------------------------------
    # Covariate effects. Both are reported in Table 2 as PERCENTAGES and both
    # enter the paper's fractional (1 + theta * x) covariate functions, so
    # each is entered here as the decimal fraction obtained by dividing the
    # printed percentage by 100.
    #
    # Dose on CL/F: Table 2 prints '-0.321' with the unit '%/mg', i.e.
    # -0.321% per mg, so theta = -0.00321 per mg in
    #   CL = theta_CL * (1 + theta_dose * (Dose - 160 mg)) * (Weight/70)^0.75.
    # The direction is a DECREASE in apparent clearance with increasing dose,
    # hence more-than-dose-proportional exposure: the multiplier is 1.385 at
    # 40 mg, 1.000 at 160 mg and 0.743 at 240 mg. Reading the printed value
    # as a per-mg FRACTION instead (-0.321/mg) is not tenable -- it would
    # drive CL/F negative just 3 mg above the reference dose.
    #
    # Asian race on F1: Table 2 '18.3%' in
    #   F1 = theta_F1 * (1 + theta_Asian race effect on F1),
    # giving F1 = 1.183 for Asian patients against 1.00 for the reference.
    # -----------------------------------------------------------------------
    e_dose_cl          <- -0.00321 ; label("Linear effect of dose on CL/F per mg above 160 mg (fraction/mg)")   # Liu 2025 Table 2 'Dose effect on CL (%/mg)' = -0.321 (95% CI -0.439/-0.203), /100
    e_race_asian_fdepot <- 0.183   ; label("Effect of Asian race on relative bioavailability F1 (fraction)")     # Liu 2025 Table 2 'Asian race effect on F1 (%)' = 18.3 (95% CI 10.9/25.8), /100

    # -----------------------------------------------------------------------
    # Inter-individual variability. Liu 2025 Table 2 column 'IIV', reported
    # for the four parameters the Methods names as carrying IIV terms:
    # 'Interindividual variability terms were included on apparent clearance
    # (CL), apparent central volume of distribution (Vc), first-order
    # absorption rate constant (Ka), and duration of zero-order absorption
    # (D1)'. No IIV is reported for Vp, Q, F1 or any covariate coefficient,
    # and none is invented here.
    #
    # The IIV cells are read as log-normal CV%, so omega^2 = log(CV^2 + 1):
    #   CL/F 48.8% -> log(0.488^2 + 1) = 0.2136135
    #   Vc/F 66.1% -> log(0.661^2 + 1) = 0.3625026
    #   ka   63.6% -> log(0.636^2 + 1) = 0.3396785
    #   Dur  56.2% -> log(0.562^2 + 1) = 0.2744783
    # Table 2 gives no footnote naming the convention and no confidence
    # interval on the IIV column that could settle it by symmetry; see the
    # vignette Errata for the alternative (omega = CV) reading and its
    # numerical consequence. Table 2 reports no off-diagonal covariance
    # terms, so all four etas are independent.
    # -----------------------------------------------------------------------
    etalcl ~ 0.2136135   # Liu 2025 Table 2 'CL (L/h)' IIV = 48.8% (shrinkage 3.5%)
    etalvc ~ 0.3625026   # Liu 2025 Table 2 'Vc (L)' IIV = 66.1% (shrinkage 32.0%)
    etalka ~ 0.3396785   # Liu 2025 Table 2 'ka (1/h)' IIV = 63.6% (shrinkage 57.8%)
    etald1 ~ 0.2744783   # Liu 2025 Table 2 'Dur (h)' IIV = 56.2% (shrinkage 49.2%)

    # -----------------------------------------------------------------------
    # Residual error. Combined additive and proportional, both from Liu 2025
    # Table 2 (epsilon shrinkage 8.7%).
    #
    # UNIT SCALE OF THE ADDITIVE TERM. Table 2 prints 'Additive RUV (mg/L)
    # 61.5', but 61.5 mg/L cannot be a concentration standard deviation for
    # this model: Figure 3 puts the adult C1D8 Cmax at 160 mg twice daily at
    # about 3 mg/L, so 61.5 mg/L would be roughly twenty times the peak
    # concentration. The printed unit tag is wrong and the number is on the
    # analysis-dataset scale of ng/mL, which Figure 2 states outright on its
    # y axis ('Prediction-corrected selpercatinib concentrations (ng/mL)',
    # with the data spanning about 1000-10000 ng/mL) and which matches the
    # 1 ng/mL assay lower limit of quantification. 61.5 ng/mL is 2.5% of the
    # roughly 2500 ng/mL median steady-state concentration -- a plausible
    # additive term. This model is therefore written on the ng/mL scale so
    # the Table 2 value is used verbatim, with the single mg -> ng/mL
    # conversion isolated in the Cc line below.
    #
    # BOTH terms are read as standard deviations rather than variances,
    # because Table 2 tags them with LINEAR scales -- a concentration unit
    # for the additive term and 'fraction' for the proportional one (a
    # variance would carry squared units) -- and because the same table
    # reports its IIV column as CV%, i.e. on the standard-deviation scale
    # throughout. See the vignette Errata for the variance reading.
    # -----------------------------------------------------------------------
    addSd  <- 61.5   ; label("Additive residual error (ng/mL)")                                           # Liu 2025 Table 2 'Additive RUV (mg/L)' = 61.5 (95% CI 37.4/85.6); unit tag read as ng/mL, see Errata
    propSd <- 0.253  ; label("Proportional residual error (fraction)")                                    # Liu 2025 Table 2 'Proportional RUV (fraction)' = 0.253 (95% CI 0.231/0.275)
  })

  model({
    # ---------------------------------------------------------------------
    # 1. Individual parameters, exactly as the Liu 2025 Table 2 covariate
    #    functions are printed:
    #      CL = theta_CL * (1 + theta_dose * (Dose - 160 mg)) * (WT/70)^0.75
    #      Vc = theta_Vc * (WT/70)^1.0
    #      Vp = theta_Vp * (WT/70)^1.0
    #      Q  = theta_Q  * (WT/70)^0.75
    #      F1 = theta_F1 * (1 + theta_Asian)
    #
    #    At WT = 70 kg, DOSE = 160 mg and RACE_ASIAN = 0 every multiplier
    #    collapses to 1 and the parameters reduce exactly to the Table 2
    #    typical values.
    #
    #    DOSE is the per-administration dose in mg and is time-varying (the
    #    paper modelled it that way to absorb dose reductions), so it is read
    #    from a data column rather than from the event table's amt.
    # ---------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (1 + e_dose_cl * (DOSE - 160)) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    q  <- exp(lq)  * (WT / 70)^e_wt_cl
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)
    fdepot <- exp(lfdepot) * (1 + e_race_asian_fdepot * RACE_ASIAN)

    # 2. Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---------------------------------------------------------------------
    # 3. ODE system with sequential zero-order then first-order absorption
    #    (Liu 2025 Figure 1). The dose is delivered into `depot` at a
    #    constant rate over the zero-order window `d1`, and `depot` then
    #    drains into `central` first-order at `ka`. Dose records must
    #    therefore carry rate = -2 so that rxode2 uses the model's
    #    dur(depot); a plain bolus would collapse the zero-order phase and
    #    overstate Cmax.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)   <- fdepot
    dur(depot) <- d1

    # ---------------------------------------------------------------------
    # 4. Observation. Dose is in mg and vc in L, so central / vc is mg/L
    #    (= ug/mL); multiplying by 1000 gives ng/mL, the scale of the
    #    analysis dataset (Figure 2 y axis) and of the 1 ng/mL assay lower
    #    limit of quantification, and the scale the Table 2 additive
    #    residual error is on.
    # ---------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
