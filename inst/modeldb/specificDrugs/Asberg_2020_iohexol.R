Asberg_2020_iohexol <- function() {
  description <- "Two-compartment intravenous population PK model for iohexol (Omnipaque 300 mg I/mL) in pediatric and adult patients referred for measured glomerular filtration rate (GFR), developed so that individual iohexol clearance -- i.e. measured GFR -- can be determined from limited sampling within 5 hours even when GFR is below 40 mL/min. All four structural parameters are scaled allometrically by total body weight standardised to 85 kg (the population median), with exponents fixed at 0.75 on CL and Q and 1 on V and Vp; no other covariate was retained (serum creatinine on CL improved population but worsened individual predictions and was rejected). Estimated with the Pmetrics non-parametric adaptive grid (NPAG) algorithm; the source reports the weighted mean, weighted median and 95% CI of the non-parametric distribution but no parametric between-subject variance, so this is a typical-value model (see the vignette Assumptions and deviations). Residual error is the Pmetrics gamma model on the published HPLC-UV assay SD polynomial. Asberg 2020, development cohort n = 176 patients aged 1-82 years, 1131 serum concentrations."
  reference <- "Asberg A, Bjerre A, Almaas R, Luis-Lima S, Robertsen I, Salvador CL, Porrini E, Schwartz GJ, Hartmann A, Bergan S. Measured GFR by Utilizing Population Pharmacokinetic Methods to Determine Iohexol Clearance. Kidney Int Rep. 2020;5(2):189-198. doi:10.1016/j.ekir.2019.11.012. PMC7000849."
  vignette <- "Asberg_2020_iohexol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Iohexol is given as an intravenous bolus into an antecubital cannula and
  # TOTAL serum iohexol was assayed by HPLC-UV (Asberg 2020 Methods,
  # 'Iohexol Administration and Sampling' and 'Bioanalytical Methods').
  compartmentData <- list(
    central = list(analyte = "iohexol", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "iohexol", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Total body weight (Asberg 2020 Results: 'Total body weight (WT) scaling resulted in a somewhat lower individual RMSE (3.0%) compared to body surface area (BSA) scaling (3.3%)'). Enters every structural parameter as the ratio (WT / 85): exponent 0.75 on CL and Q, exponent 1.00 on V and Vp (Table 2 footnote a). 85 kg is the population median body weight of the development cohort (Methods, 'Model Development': body size measures 'were centralized to population median values'). Development-cohort weight by age bin (Table 1, mean +/- SD): 11.4 kg (0-2 yr, n = 2), 45.0 +/- 24.8 (2-21 yr), 78.9 +/- 18.1 (21-60 yr), 87.9 +/- 15.4 (> 60 yr); whole-study mean 74.1 +/- 24.4 kg.",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Plasma creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on CL as CL = CLs * WTc^0.75 * EXP(CLCRE * CREAT) (Asberg 2020 Results, 'Model Development'). It improved the population fit (AIC 6858 -> 6805, population RMSE 52% -> 32%) but raised the individual RMSE 'severalfold to 20%', so it was rejected because the model's purpose is individual (Bayesian) clearance determination. No estimate of CLCRE is printed. Development-cohort plasma creatinine by age bin (Table 1, mean SD): 22, 98 95, 276 168, 249 94 umol/L.",
      source_name = "CREAT"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as the allometric size descriptor in place of total body weight (individual RMSE 3.3% vs 3.0%, AIC 6876 vs 6858) and not chosen (Asberg 2020 Results). Height, BMI and fat-free mass were also screened as size descriptors (Methods). BSA is used outside the model by the clinical Brochner-Mortensen 2-point comparator, f = 0.0032 * BSA^-1.3.",
      source_name = "BSA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 176L,
    n_studies = 3L,
    n_samples = 1131L,
    age_range = "1-82 years across the whole study (Results, 'Study Population'); development cohort by age bin (Table 1): 0-2 yr n = 2, 2-21 yr n = 38, 21-60 yr n = 63, > 60 yr n = 73. The authors consider the model validated only from 5 to 77 years (Discussion).",
    weight_range = "Full range not reported. Development cohort mean by age bin 11.4, 45.0, 78.9 and 87.9 kg (Table 1); population median 85 kg (the allometric reference).",
    sex_female_pct = 24,
    race_ethnicity = "Mainly Caucasian (Results, 'Study Population'); no breakdown reported.",
    disease_state = "Pediatric and adult patients scheduled for a clinical iohexol measured-GFR investigation at Oslo University Hospital-Rikshospitalet (prospective, 2014-2017) plus previously published cohorts from Tenerife and Rochester (refs 17 and 18), spanning normal to severely reduced renal function. Measured GFR 14-149 mL/min across the development and validation cohorts; development-cohort plasma creatinine means 22-276 umol/L by age bin (Table 1).",
    renal_function = "GFR 14-149 mL/min (absolute, not normalised to 1.73 m^2) across the 219 patients; the model was designed to remain accurate below 40 mL/min.",
    dose_range = "Single IV bolus of 5 mL Omnipaque 300 mg I/mL (3235 mg iohexol) flushed with 10 mL saline; children under 2 years received 2 mL (about 1294 mg). The exact administered dose was determined by weighing the syringe.",
    regions = "Norway, Spain (Canary Islands) and the United States",
    sampling = "Development cohort: 1131 concentrations, median 7 (range 2-12) per patient; 44 of 176 patients sampled in the first 2 h (10, 20, 30, 45, 60, 90 min) and the remainder only from 2 h onward, up to 24 h. HPLC-UV, CV < 6%, LLOQ 20 mg/L, linear 20-1100 mg/L.",
    notes = "Validation cohort (n = 43, 395 concentrations) was used for external validation and limited-sampling evaluation; the optimal 4-sample schedule within 5 h was 10 min, 30 min, 2 h and 5 h. NPAG final model: 151 support points after 3209 cycles; population and individual RMSE 52.4% and 3.0%; individual bias -0.0328 mg/L and imprecision 1.81 mg/L; final gamma 1.947."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters: Asberg 2020 Table 2, 'Median' column of the
    # Pmetrics NPAG distribution (weighted by support-point probability),
    # each for an 85 kg patient. Table 2 footnote a: 'To obtain individual
    # values, CLs and Qs should be multiplied by (actual patient weight
    # [kg]/85 [kg])^0.75 and Vs and Vps by (actual patient weight
    # [kg]/85 [kg])^1.00.'
    #
    # The MEDIAN rather than the MEAN column is the typical value: the
    # median parameter set at 85 kg reproduces the prediction-corrected
    # observed median of the supplement's Figure S2 pcVPC (about 273, 183,
    # 143 and 24 mg/L at 0.2, 1, 2 and 24 h) to within a few percent, while
    # the mean parameter set undershoots the 24 h value threefold. See the
    # vignette.
    # ---------------------------------------------------------------------

    lcl <- log(1.55)
    label("Iohexol clearance (= GFR) at WT = 85 kg (L/h)")
    # Table 2: CLs median 1.55 L/h (mean 2.60; 95% CI 2.22-2.97; shrinkage 0.8%)

    lq <- log(5.99)
    label("Intercompartmental clearance at WT = 85 kg (L/h)")
    # Table 2: Qs median 5.99 L/h (mean 9.44; 95% CI 8.03-11.42; shrinkage 4.0%)

    lvc <- log(10.47)
    label("Central volume of distribution at WT = 85 kg (L)")
    # Table 2: Vs median 10.47 L (mean 10.96; 95% CI 10.25-11.68; shrinkage 4.1%)

    lvp <- log(8.02)
    label("Peripheral volume of distribution at WT = 85 kg (L)")
    # Table 2: Vps median 8.02 L (mean 9.44; 95% CI 8.66-10.23; shrinkage 7.3%)

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent of WT on CL and Q (unitless)")
    # Table 2 footnote a: (WT/85)^0.75 on CLs and Qs; stated, not estimated

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent of WT on V and Vp (unitless)")
    # Table 2 footnote a: (WT/85)^1.00 on Vs and Vps; stated, not estimated

    # ---------------------------------------------------------------------
    # Inter-individual variability: NONE is encoded.
    #
    # NPAG estimates a discrete non-parametric joint density (151 support
    # points) rather than an omega matrix, and Table 2 prints no SD, CV% or
    # variance. The 95% CI column brackets the MEAN (2.22-2.97 around mean
    # 2.60 vs median 1.55 for CL) and so is precision of the mean, not
    # between-subject spread. The supplement's Figure S1 marginal
    # support-point distributions are flat, bounded and pile up on the
    # search-grid limits (CL 0.2-10 L/h, V 1-20 L, Vp 1-25 L, Q 0-55 L/h),
    # i.e. visibly not log-normal. The vignette reports two lognormal
    # reconstructions (from mean/median, and from the CI width) for users
    # who need a stochastic prior; neither is encoded. `fixed(0)` etas are
    # deliberately not used because a zero-variance omega breaks rxode2's
    # Cholesky sampler.
    # ---------------------------------------------------------------------

    # ---------------------------------------------------------------------
    # Residual error: the Pmetrics gamma model, error = SD * gamma, with the
    # HPLC-UV assay SD polynomial printed in Methods ('Model Development'):
    #   SD = 0.1523073 + 0.01747435*obs - 0.000003919581*obs^2   (mg/L)
    # and the final fitted gamma = 1.947 (Results). All four are published
    # constants of the final model, hence fixed(). Pmetrics evaluates the
    # polynomial on the OBSERVED concentration; nlmixr2 evaluates it on the
    # prediction (see the vignette Assumptions and deviations).
    # ---------------------------------------------------------------------

    addSd <- fixed(0.1523073)
    label("Assay SD polynomial intercept C0 (mg/L)")
    # Methods, Model Development: SD = 0.1523073 + ...

    propSd <- fixed(0.01747435)
    label("Assay SD polynomial linear coefficient C1 (fraction)")
    # Methods, Model Development: ... + 0.01747435*obs ...

    quadSd <- fixed(-0.000003919581)
    label("Assay SD polynomial quadratic coefficient C2 (L/mg)")
    # Methods, Model Development: ... - 0.000003919581*obs^2

    gammaSd <- fixed(1.947)
    label("Pmetrics gamma multiplier on the assay SD (unitless)")
    # Results, Model Development: 'The final gamma value was 1.947.'
  })

  model({
    wt_ref <- 85 # kg, population median (Table 2 footnote a)

    # Final-model parameterisation (Results): CL = CLs * WTc^0.75,
    # Q = Qs * WTc^0.75, V = Vs * WTc, Vp = Vps * WTc, with WTc = WT / 85.
    cl <- exp(lcl) * (WT / wt_ref)^e_wt_cl_q
    q <- exp(lq) * (WT / wt_ref)^e_wt_cl_q
    vc <- exp(lvc) * (WT / wt_ref)^e_wt_vc_vp
    vp <- exp(lvp) * (WT / wt_ref)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total serum iohexol (mg/L). Iohexol clearance is the measured GFR.
    Cc <- central / vc

    sdCc <- gammaSd * (addSd + propSd * Cc + quadSd * Cc^2)
    Cc ~ add(sdCc)
  })
}
