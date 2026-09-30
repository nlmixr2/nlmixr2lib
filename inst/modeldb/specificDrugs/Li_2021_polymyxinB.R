Li_2021_polymyxinB <- function() {
  description <- "One-compartment intravenous population PK model for polymyxin B in adult renal transplant recipients, most with renal dysfunction (Li 2021). Cockcroft-Gault creatinine clearance is a power covariate on clearance, normalized to the cohort median of 22.2 mL/min with exponent 0.14. Volume of distribution carries no covariate. Exponential IIV on CL and V; proportional residual error. Fit in Phoenix NLME (FOCE-ELS)."
  reference <- paste(
    "Li Y, Deng Y, Zhu ZY, Liu YP, Xu P, Li X, Xie YL, Yao HC, Yang L,",
    "Zhang BK, Zhou YG.",
    "Population pharmacokinetics of polymyxin B and dosage optimization in",
    "renal transplant patients.",
    "Front Pharmacol. 2021;12:727170.",
    "doi:10.3389/fphar.2021.727170. PMCID PMC8424097.",
    sep = " "
  )
  vignette <- "Li_2021_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Li 2021 Methods: polymyxin B1 and B2 were measured in plasma by
  # HPLC-MS/MS and modelled as total polymyxin B in mg/L; the single disposition
  # compartment is the plasma/central compartment of the one-compartment model.
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Li 2021 Methods: 'CrCL was calculated according to the Cockcroft-Gault equation'. Raw mL/min, not normalized to 1.73 m^2 BSA; stored under the canonical CRCL column, which accepts raw mL/min provided the assay form is documented per model (precedent: Yang_2025_polymyxinB.R). Applied as a power covariate on CL, CL = 1.18 * (CrCL/22.2)^0.14 * exp(eta_CL) (Results, 'Development of the PPK Model'). The normalizing constant 22.2 mL/min is printed in that equation and equals the cohort median of Table 1 (22.2, range 4.29-90.7); Methods state 'The median of the covariate was used to normalize the covariate'. The Discussion quotes a slightly different median, 20.89 mL/min (range 4.29-78.84); the printed final-model equation is followed. Eleven of the 50 patients (22%) received continuous renal replacement therapy during treatment; the paper does not describe how CrCL was assigned for them, and CRRT was not screened as a covariate. Time-fixed per subject in this model.",
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Screened (Li 2021 Methods, stepwise forward inclusion p < 0.05 / backward elimination p < 0.01) and not retained; CrCL was the only covariate retained, and none on V."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened and not retained. Li 2021 Discussion: 'a general lack of a significant linear relationship between weight and polymyxin B clearance was found in the renal transplant patient population'."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened and not retained. Table 1 median 34.7 g/L (22.7-49.0)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened and not retained. Table 1 median 20.2 mmol/L (2.65-61.9)."
    ),
    OTHER_TABLE1_LABS = list(
      description = "Remaining Table 1 laboratory parameters screened as candidate covariates",
      units = "(various)",
      type = "continuous",
      notes = "ALT, AST, total bilirubin, direct bilirubin and uric acid were screened per Li 2021 Methods ('The covariates considered for the modeling included age, weight, ALT, AST, TBIL, DBIL, ALB, BUN, UA, and CrCL') and none was retained. Grouped into a single entry because the paper reports the screen collectively rather than per analyte."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = "18-66 years (median 43.5)",
    weight_range = "57.8 +/- 12.4 kg (mean +/- SD, actual body weight)",
    sex_female_pct = 36,
    race_ethnicity = "Chinese (single-centre cohort; race not otherwise reported)",
    disease_state = "Adult renal transplant recipients receiving intravenous polymyxin B for >= 48 h for a (suspected) Gram-negative infection; pulmonary 90%, urinary tract 8%, bloodstream 2%. Renal dysfunction (CrCL <= 80 mL/min) in 92%, 29 of 50 with CrCL < 30 mL/min; 11 (22%) received continuous renal replacement therapy. Hypertension 56%, diabetes mellitus 14%.",
    dose_range = "Polymyxin B (sulfate) IV every 12 h as a 60-120 min infusion for at least 3 days; maintenance 40 mg q12h in 44%, 50 mg q12h in 40%, other doses in 16%; only two patients received a loading dose. Monte Carlo simulations used 2-h infusions of 50 mg loading + 30 or 40 mg q12h, 75 or 100 mg loading + 50 mg q12h, and 150 mg loading + 75 mg q12h.",
    regions = "China (The Second Xiangya Hospital, Central South University, Changsha)",
    renal_function = "Cockcroft-Gault CrCL median 22.2 mL/min (range 4.29-90.7).",
    n_observations = "151 plasma polymyxin B concentrations (0.44-8.15 mg/L); 1-6 samples per patient, drawn 30 min before the sixth dose and at 0, 0.5, 1, 2, 4, 6 and 8 h after the end of that infusion. LLOQ 0.03 mg/L.",
    notes = "Prospective study, ChiCTR1900022231. Baseline demographics in Li 2021 Table 1; sampling distribution in Table 2. Concentrations by validated HPLC-MS/MS (polymyxin B1 + B2). Estimation in Phoenix NLME 8.1 with FOCE-ELS; 1,000-sample bootstrap and pc-VPC reported."
  )

  ini({
    # Structural parameters -- Li 2021 Table 4 (final model), typical values at
    # the cohort median CrCL of 22.2 mL/min.
    lcl <- log(1.18); label("Clearance CL (L/h) at CRCL = 22.2 mL/min") # Li 2021 Table 4: CL = 1.18 L/h (%CV 4.15; bootstrap median 1.17, 95% CI 1.08-1.27)
    lvc <- log(12.09); label("Volume of distribution V (L)") # Li 2021 Table 4: V = 12.09 L (%CV 6.52; bootstrap median 11.98, 95% CI 10.71-13.67)

    # Covariate effect on CL -- Li 2021 Results, final-model equation:
    #   CL (L/h) = 1.18 * (CrCL/22.2)^0.14 * exp(eta_CL)
    # Estimated (it carries a %CV and a bootstrap CI), so not fixed().
    e_crcl_cl <- 0.14; label("Power exponent on (CRCL / 22.2 mL/min) for CL (unitless)") # Li 2021 Table 4: Theta CrCL = 0.14 (%CV 29.35; bootstrap median 0.14, 95% CI 0.05-0.21)

    # Inter-individual variability, exponential (Methods). Table 4 prints the
    # variances rounded to two decimals (omega^2 V = 0.04, omega^2 CL = 0.06)
    # and, in the %CV column, the corresponding CV percentages 20.74 and 24.49.
    # The CV column is sqrt(omega^2) * 100: sqrt(0.06) = 0.2449 exactly, and the
    # Discussion calls 24.49% 'the interindividual variability of CL' in the
    # final model. The variances are therefore taken from the more precise CV
    # column: 0.2074^2 = 0.0430 (rounds to the printed 0.04) and
    # 0.2449^2 = 0.0600 (rounds to the printed 0.06). The IIV-only spread of the
    # day-3 AUC percentiles in Table 6 (P90/P10 = 1.85) independently implies an
    # SD of log(CL) of about 0.24.
    etalcl ~ 0.0600 # Li 2021 Table 4: omega^2 CL = 0.06 (shrinkage 7.16%), %CV 24.49 -> 0.2449^2
    etalvc ~ 0.0430 # Li 2021 Table 4: omega^2 V = 0.04 (shrinkage 40.70%), %CV 20.74 -> 0.2074^2

    # Residual variability -- proportional error selected (Table 3; Results).
    # Phoenix NLME reports the residual epsilon on the SD scale.
    propSd <- 0.17; label("Proportional residual error (fraction)") # Li 2021 Table 4: sigma = 0.17 (bootstrap median 0.17, 95% CI 0.13-0.21)
  })
  model({
    # Individual PK parameters -- Li 2021 Results, 'Development of the PPK Model':
    # CL (L/h) = 1.18*(CrCL/22.2)^0.14*exp(eta CL); V (L) = 12.09*exp(eta V).
    cl <- exp(lcl + etalcl) * (CRCL / 22.2)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment disposition with first-order elimination; polymyxin B is
    # given as an intravenous infusion directly into the central compartment.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
