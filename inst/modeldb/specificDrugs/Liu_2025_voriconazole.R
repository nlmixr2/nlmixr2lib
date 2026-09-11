Liu_2025_voriconazole <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption for intravenous and oral voriconazole in elderly Chinese inpatients aged 60 years and over (Liu 2025); apparent clearance carries median-normalized power effects of procalcitonin and total bile acids together with an exponential age effect centred at 72 years, so inflammation, cholestasis and advancing age each predict slower clearance. This is the population pharmacokinetic layer of a paper whose headline product is a machine-learning ensemble that consumes the empirical-Bayes CL/F as its most important feature."
  reference <- "Liu R, Ma P, Chen D, Yu M, Xie L, Zhao L, Huang Y, Shang S, Chen Y. A real-time plasma concentration prediction model for voriconazole in elderly patients via machine learning combined with population pharmacokinetics. Drug Des Devel Ther. 2025;19:4021-4034. doi:10.2147/DDDT.S495050"
  vignette <- "Liu_2025_voriconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot   = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "voriconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    PCT = list(
      description        = "Serum procalcitonin concentration",
      units              = "ug/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-normalized power effect on apparent clearance, (PCT / 0.19)^e_pct_cl, per the published",
        "final-model equation on Liu 2025 p. 4024. The divisor 0.19 ug/L is the training-group median of",
        "Liu 2025 Table 1 (IQR 0.13-0.28 ug/L; testing group median also 0.19, IQR 0.12-0.29), so the",
        "printed typical CL/F of 4.35 L/h is the value at the cohort-median procalcitonin -- this is a",
        "median-CENTERED form, unlike the median-scaled-but-uncentered CRP term in the sibling model",
        "Ling_2024_voriconazole.R. The negative exponent encodes the acute-phase suppression of hepatic",
        "cytochrome-P450 activity: a septic rise in procalcitonin lowers voriconazole clearance and raises",
        "exposure. Liu 2025 reports this covariate as ranking fourth by mean absolute SHAP value in the",
        "downstream machine-learning ensemble (Figure 3B).",
        "Whether the column is a single baseline value or a per-sample time-varying measurement is not",
        "stated; the retrospective design drew covariates from the electronic medical record alongside each",
        "therapeutic-drug-monitoring sample, and missing values were median-imputed (Methods, 'Data",
        "Collection and Processing'), which is consistent with a per-sample series.",
        "C-reactive protein, the other inflammation marker, was excluded from the whole analysis because",
        "more than 50% of values were missing (Liu 2025 Discussion, limitation four)."
      ),
      source_name        = "PCT"
    ),
    TBA = list(
      description        = "Total serum bile acids",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Median-normalized power effect on apparent clearance, (TBA / 3.95)^e_tba_cl, per the published",
        "final-model equation on Liu 2025 p. 4024. The divisor 3.95 umol/L is the training-group median of",
        "Liu 2025 Table 1 (IQR 2.50-7.45 umol/L; the testing group median is lower, 3.59, IQR 2.40-4.32).",
        "The negative exponent means that a cholestatic rise in bile acids lowers voriconazole clearance.",
        "Liu 2025 Discussion notes that total bile acids was the only liver-function marker retained by the",
        "stepwise covariate search, ahead of ALT, AST, ALP, GGT and total bilirubin, and cites evidence that",
        "elevated bile acids indicate the hepatic dysfunction that impairs voriconazole metabolism.",
        "Fasting-versus-postprandial status is not stated. Baseline-versus-time-varying status is likewise",
        "not stated; see the PCT notes for the same consideration."
      ),
      source_name        = "TBA"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Exponential effect on apparent clearance centred at 72 years, exp(e_age_cl * (AGE - 72.0)), per",
        "the published final-model equation on Liu 2025 p. 4024. The centring constant 72.0 is the",
        "training-group median age of Liu 2025 Table 1 (IQR 67-78 years; testing group median 73, IQR",
        "68-83). Note the functional form differs from the two laboratory covariates: age enters as an",
        "exponential of the deviation from 72 years, not as a power of a ratio.",
        "Every subject was aged 60 years or over by the inclusion criteria, so the covariate is only",
        "supported over roughly 60-90 years; do not extrapolate below 60. Over the observed interquartile",
        "range 67-78 years the effect spans exp(-0.017 * (78 - 67)) = 0.83, a 17% clearance decrement.",
        "Liu 2025 Discussion attributes the effect to the age-related decline in hepatic function."
      ),
      source_name        = "Age"
    )
  )

  # Liu 2025 does not publish the set of covariates entered into the stepwise
  # search -- the Methods state only that an exploratory correlation analysis
  # against the empirical-Bayes estimates preceded a standard forward-inclusion
  # (dOFV 3.84) / backward-elimination (dOFV 10.83) procedure. The 31 features
  # of Table 1 are the *machine-learning* feature set, not a stated population-PK
  # screening set, so they are deliberately NOT recorded here as screened-and-
  # rejected covariates; doing so would assert a provenance the paper does not
  # support. The Table 1 distributions are reproduced in population$notes.

  population <- list(
    species        = "human",
    n_subjects     = 270L,
    n_studies      = 1L,
    n_observations = 393L,
    age_range      = ">= 60 years (inclusion criterion)",
    age_median     = "72 years (training), 73 years (testing)",
    weight_median  = "58 kg",
    sex_female_pct = 28.3,
    race_ethnicity = c(Chinese = 100),
    disease_state  = "Elderly hospitalized inpatients receiving voriconazole for more than 3 days and undergoing therapeutic drug monitoring. Patients on dialysis, those with concentrations below the lower limit of quantification, and those with incomplete administration records were excluded.",
    dose_range     = "Median daily dose 6.90 mg/kg/day (IQR 6.00-8.00). At the median weight of 58 kg this is 400 mg/day, i.e. the standard 200 mg twice-daily maintenance regimen. 86.0% of training records were intravenous and 14.0% oral (testing 83.5% / 16.5%). Median total treatment duration 6 days (IQR 4-11); median time after dose at sampling 10.57 h (IQR 9.68-11.50), consistent with pre-dose trough sampling on a 12-hourly schedule. Only 5.7% of daily doses exceeded 10 mg/kg/day, the threshold above which the authors cite nonlinear elimination as material.",
    regions        = "Single center: the First Affiliated Hospital of Army Medical University, Chongqing, China.",
    notes          = paste(
      "Retrospective single-center study, March 2022 - December 2023. 393 therapeutic-drug-monitoring",
      "concentrations from 270 patients, measured by LC-MS/MS (Shimadzu LC-30AD with AB Sciex QTRAP 5500).",
      "An additional 48 patients contributing 76 concentrations, enrolled December 2023 - March 2024, form",
      "an external clinical validation group used only for the machine-learning layer.",
      "IMPORTANT: Liu 2025 Table 1 is tabulated per RECORD, not per patient -- the 314 / 79 training and",
      "testing counts are concentrations, and the paper states the randomization was performed at the",
      "sample level rather than the patient level. Baseline statistics below are therefore",
      "record-weighted. Observed voriconazole concentration median 3.40 mg/L (IQR 2.21-5.22) in training",
      "and 3.38 (2.19-5.18) in testing.",
      "Training-group covariate distributions (median and IQR unless stated): age 72 (67-78) years, weight",
      "58 (55-61) kg, BMI 22.66 (21.48-23.40) kg/m2, ALT 17.90 (10.90-30.63) U/L, AST 27.95 (20.10-41.05)",
      "U/L, ALP 85.50 (72.00-115.25) U/L, GGT 54.90 (30.93-107.93) U/L, total bilirubin 10.30 (8.00-14.00)",
      "umol/L, total bile acids 3.95 (2.50-7.45) umol/L, total protein 60.75 (56.10-65.70) g/L, albumin",
      "32.40 +/- 3.56 g/L, eGFR 94.88 (75.95-119.95) mL/min/1.73m2, serum creatinine 67.00 (53.55-80.75)",
      "umol/L, creatinine clearance 68.61 (51.73-88.91) mL/min, platelets 203.90 +/- 117.28 x10^9/L, white",
      "blood cells 7.45 (5.02-11.13) x10^9/L, hemoglobin 99.00 (83.75-116.00) g/L, neutrophils 76.15",
      "(65.35-86.45)%, procalcitonin 0.19 (0.13-0.28) ug/L, interleukin-6 27.34 (20.09-39.31) ng/L.",
      "Co-medications: glucocorticoid 46.2%, meropenem 23.9%, proton-pump inhibitor 55.4%. Concomitant",
      "liver disease 43.6%; malignant hematological disease 16.9%.",
      "Missing covariate values were median-imputed in Python before analysis (Methods, 'Data Collection",
      "and Processing'), so the tabulated distributions are post-imputation.",
      "CYP2C19 genotype was NOT collected. The authors argue this matters less in the elderly, citing their",
      "own earlier finding of a smaller CYP2C19 effect in elderly than in younger adults, and note that the",
      "empirical-Bayes CL/F partly absorbs the missing genotype information (Discussion, limitation one).",
      "This is the principal structural difference from the sibling models Hu_2023_voriconazole.R,",
      "Lin_2018_voriconazole.R and Ling_2024_voriconazole.R, all of which carry an explicit CYP2C19 term.",
      "Final-model estimates, 1000-replicate bootstrap (998 successful, 99.8%) and shrinkage per Table 2;",
      "the published final-model equations are the three display equations on p. 4024."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Reference subject for the typical-value equations: procalcitonin
    # 0.19 ug/L, total bile acids 3.95 umol/L and age 72.0 years -- the
    # training-group medians of Liu 2025 Table 1. All three covariate terms
    # are median-centred, so at the reference the two power ratios and the
    # age exponential are each exactly 1 and CL/F is the printed 4.35 L/h.
    #
    # Every disposition parameter in this model is APPARENT (CL/F, V/F).
    # Liu 2025 pooled intravenous (86% of records) and oral (14%) data and
    # reported CL/F and V/F throughout without estimating a separate
    # bioavailability term, so no lfdepot is encoded here: F is implicitly
    # 1 for both routes. For the intravenous majority CL/F is simply CL.
    # See the vignette Errata for the consequence of that pooling.
    # ---------------------------------------------------------------------

    # Absorption. Estimating ka gave an RSE of 145%, so it was fixed to the
    # literature value; the paper states model parameters and OFV were
    # essentially unchanged after fixing.
    lka <- fixed(log(1.1)); label("Absorption rate constant (1/h)")  # Liu 2025 Methods 'Population Pharmacokinetic Modeling': "Ka was fixed at a value of 1.1 h-1 following reference 18" (Pascual A et al., Clin Infect Dis 2012;55(3):381-390); Table 2 lists "1.1 fixed" in the base model, final model and bootstrap columns. The sibling model Ling_2024_voriconazole.R fixes ka to the same value from the same Pascual 2012 source.

    lcl <- log(4.35); label("Apparent clearance CL/F at the reference covariate values (L/h)")  # Liu 2025 Table 2 final model: CL/F = 4.35 L/h (RSE 4.4%, bootstrap median 4.33, 95% CI 3.98-4.71). Confirmed by the leading constant of the published CL/F equation, p. 4024. Base model 4.16 L/h.

    lvc <- log(140); label("Apparent volume of distribution V/F (L)")  # Liu 2025 Table 2 final model: V/F = 140 L (RSE 12.6%, bootstrap median 139, 95% CI 102-177). Confirmed by the published equation "V/F(L) = 140", p. 4024. Base model 139 L. No covariate was retained on V/F.

    # Covariate effects on CL/F. Note the two laboratory covariates enter as
    # powers of a median-normalized ratio while age enters as an exponential
    # of the deviation from the median -- this asymmetry is exactly as
    # printed in the paper's own equation and is not a transcription slip.
    e_pct_cl <- -0.209; label("Power exponent for PCT on CL/F (unitless)")  # Liu 2025 Table 2 final model, "PCT on CL/F" = -0.209 (RSE 20.1%, bootstrap median -0.210, 95% CI -0.305 to -0.124); appears as the exponent of (PCT/0.19) in the published equation, p. 4024
    e_tba_cl <- -0.158; label("Power exponent for TBA on CL/F (unitless)")  # Liu 2025 Table 2 final model, "TBA on CL/F" = -0.158 (RSE 27.9%, bootstrap median -0.155, 95% CI -0.255 to -0.070); appears as the exponent of (TBA/3.95) in the published equation, p. 4024
    e_age_cl <- -0.017; label("Exponential coefficient for AGE on CL/F, per year above 72 (1/year)")  # Liu 2025 Table 2 final model, "Age on CL/F" = -0.017 (RSE 24.9%, bootstrap median -0.017, 95% CI -0.026 to -0.009); appears as exp[-0.017 * (AGE - 72.0)] in the published equation, p. 4024

    # IIV. Liu 2025 Methods: "An exponential model was applied to explain the
    # inter-individual variability (IIV) of the parameters", and the published
    # equation ends in * exp(eta_CL), confirming CL/F = typical * exp(eta).
    # Table 2 reports the final-model row as "eta CL (%)" = 43.9.
    #
    # That row is read here as the log-scale STANDARD DEVIATION x100, giving
    # a variance of 0.439^2 = 0.192721 -- NOT as a coefficient of variation
    # requiring omega^2 = log(1 + CV^2), which would give 0.4198^2 = 0.176. The
    # two readings differ by only 4.5% in omega, but three independent lines
    # of evidence favour the SD reading and are recorded here because the next
    # reader will otherwise "correct" it:
    #
    # (1) The residual rows of the SAME table settle the table's convention.
    #     "Prop_error (%)" = 28.9 can only be an SD: read as a variance it
    #     would imply a residual CV of sqrt(0.289) = 53.8%, far above the
    #     imprecision of the validated LC-MS/MS assay the paper cites
    #     (reference 17). All rows of one table resolve the same way.
    # (2) Liu 2025 Table 1 publishes the empirical-Bayes CL/F distribution,
    #     median 4.00 with IQR 3.00-5.53 L/h, which reads out a log-scale SD
    #     of log(5.53/3.00) / (2 * qnorm(0.75)) = 0.4534. Reconstructing that
    #     spread from Table 2 -- covariate-explained variance
    #     (0.538^2 - 0.439^2) plus the shrunken eta variance, using NONMEM's
    #     default SD-scale eta-shrinkage of 25.0% from Table 2 -- predicts
    #     sqrt((0.538^2 - 0.439^2) + (0.439 * 0.75)^2) = 0.4529, within 0.1%
    #     of the observed 0.4534. The CV reading predicts 0.4208, 7.2% low.
    #     Honest caveat: under a variance-scale shrinkage convention the two
    #     readings swap (0.4912 vs 0.4584), so this line is corroborating
    #     rather than decisive on its own -- it rests on NONMEM/PsN reporting
    #     eta-shrinkage on the SD scale, which is the default.
    # (3) The same-drug sibling Ling_2024_voriconazole.R -- which also fixes
    #     ka to 1.1 /h from the same Pascual 2012 source -- reads its own
    #     Table 2 "IIV (%)" rows as SDs on the same reasoning.
    etalcl ~ 0.192721  # Liu 2025 Table 2 final model: eta CL = 43.9% (RSE 9.8%, eta-shrinkage 25.0%, bootstrap median 43.9%, 95% CI 35.2-52.2%); var = 0.439^2. Base model 53.8% (shrinkage 21.1%), i.e. the three covariates cut IIV from 53.8% to 43.9%.

    # Residual error. Liu 2025 Methods: "Residual variability was evaluated
    # through a combined additive and proportional error model", and Table 2
    # reports both components, so both are carried.
    propSd <- 0.289; label("Proportional residual error (fraction)")  # Liu 2025 Table 2 final model: Prop_error = 28.9% (RSE 12.3%, eps-shrinkage 19.0%, bootstrap median 27.9%, 95% CI 19.6-34.7%). Base model 31.4%.
    addSd  <- 0.885; label("Additive residual error (mg/L)")  # Liu 2025 Table 2 final model: Add_error = 0.885 mg/L (RSE 9.9%, eps-shrinkage 19.0%, bootstrap median 0.870, 95% CI 0.717-1.066). Base model 0.893 mg/L.
  })

  model({
    # Individual apparent clearance. Reproduces the Liu 2025 published
    # final-model equation verbatim (p. 4024, first display equation):
    #
    #   CL/F (L/h) = 4.35 * (PCT/0.19)^-0.209 * (TBA/3.95)^-0.158
    #                     * e^[-0.017 * (AGE - 72.0)] * e^(eta_CL)
    #
    # The three divisors / centring constants 0.19, 3.95 and 72.0 are the
    # training-group medians of Table 1, so a subject at the cohort median
    # of all three covariates has the printed typical CL/F of 4.35 L/h.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (PCT / 0.19)^e_pct_cl *
      (TBA / 3.95)^e_tba_cl *
      exp(e_age_cl * (AGE - 72.0))

    # Apparent volume of distribution: no covariate was retained, and no IIV
    # was estimated on V/F. Liu 2025 states IIV "was successfully estimated
    # on clearance divided by bioavailability (CL/F)" -- clearance only --
    # and Table 2 reports a single inter-individual variability row. The
    # absence of an etalvc is therefore a faithful transcription of the
    # published model, not an omission.
    vc <- exp(lvc)

    # One-compartment disposition with first-order absorption (NONMEM
    # ADVAN2 TRANS2). Oral doses enter the depot; the intravenous majority
    # bypasses the depot and is dosed directly into central. Because Liu
    # 2025 estimated no bioavailability term, no f(depot) is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # Observation. Amounts in mg over a volume in L give mg/L, the unit in
    # which Liu 2025 reports both the observed concentrations (Table 1) and
    # the additive residual error (Table 2).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
