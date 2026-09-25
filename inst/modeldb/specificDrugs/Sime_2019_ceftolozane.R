Sime_2019_ceftolozane <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for UNBOUND ceftolozane in",
    "twelve critically ill adults without renal dysfunction, admitted to a",
    "quaternary referral intensive care unit in Brisbane, Australia. Fitted",
    "non-parametrically with the NPAG algorithm in Pmetrics 1.5.2 to 133",
    "directly-measured unbound plasma concentrations. Clearance is ADDITIVE in a",
    "covariate-free intercept and an arm linear in measured urinary creatinine",
    "clearance (paper: CL = intercept + slope * CLcr_urinary), and central volume",
    "scales linearly with total body weight (V1 = V * WT/80). The",
    "intercompartmental transfer is parameterised directly as the rate constants",
    "Kcp and Kpc rather than as Q and Vp. Every parameter carries its own",
    "inter-individual variability, taken from the mean and SD of the NPAG",
    "support-point distribution. Residual unexplained variability is carried as",
    "fixed(0) because neither the selected Pmetrics error model nor its assay",
    "error polynomial coefficients were published. Ceftolozane and tazobactam were",
    "fitted in two separate NPAG runs, with different body-weight exponents on",
    "volume, and are supplied as two separate model files; see",
    "modellib('Sime_2019_tazobactam') for the partner component of the fixed 2:1",
    "ceftolozane-tazobactam combination.",
    sep = " "
  )
  reference <- paste(
    "Sime FB, Lassig-Smith M, Starr T, Stuart J, Pandey S, Parker SL, Wallis SC,",
    "Lipman J, Roberts JA. Population pharmacokinetics of unbound ceftolozane and",
    "tazobactam in critically ill patients without renal dysfunction.",
    "Antimicrob Agents Chemother. 2019;63(10):e01265-19.",
    "doi:10.1128/AAC.01265-19. PMCID: PMC6761554.",
    "All structural and variability estimates are Table 2, 'Ceftolozane' block.",
    "The covariate equations are the Results narrative ('CL = intercept + slope *",
    "CL CRurinary' and 'V 1 = V * WT/80'). No supplement was deposited with the",
    "article.",
    sep = " "
  )
  vignette <- "Sime_2019_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The assay measured the UNBOUND plasma concentration
  # directly (Methods, 'Ceftolozane-tazobactam assay': the unbound fraction was
  # isolated by ultracentrifugation before injection), and the full administered
  # milligram amount was dosed into the model, so vc is the apparent volume
  # relating total drug amount to UNBOUND concentration and Cc below is an
  # unbound concentration. No separate fu parameter exists or is needed. The
  # assay validation reports ceftolozane unbound fractions of 90%, 99% and 101%
  # at total concentrations of 160, 20 and 3 mg/liter, so unbound and total are
  # nearly interchangeable for this drug in this cohort -- but the fitted
  # parameters are the unbound ones.
  compartmentData <- list(
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "MEASURED urinary creatinine clearance from a timed urine collection,",
        "body-surface-area normalized to 1.73 m^2. Not a Cockcroft-Gault or",
        "CKD-EPI estimate"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "ASSAY FORM. Methods, 'Clinical data', lists 'renal function markers",
        "(serum creatinine concentration and urinary creatinine clearance)' among",
        "the collected variables, and Table 1 reports 'Urinary creatinine",
        "clearance (ml/min/1.73 m 2)' separately from 'Serum creatinine (umol/",
        "liter)'. The Results call the retained covariate 'CL CRurinary is",
        "measured urinary creatinine clearance'. Do NOT substitute an estimated",
        "creatinine clearance: patients with augmented renal clearance are exactly",
        "the population in which Cockcroft-Gault and measured clearance diverge",
        "most, and the whole point of the paper is the augmented-clearance tail.",
        "FUNCTIONAL FORM. Results: 'ceftolozane and tazobactam clearance linearly",
        "increased with an increase in urinary creatinine clearance. The final",
        "covariate model for clearance of both ceftolozane and tazobactam was",
        "expressed as CL = intercept + slope * CL CRurinary'. Additive, linear,",
        "NOT a power function.",
        "NORMALIZATION BY 100 -- DERIVED, NOT PRINTED. The paper never states the",
        "scale on which CLcr enters, and the printed slope is not usable at face",
        "value in mL/min/1.73 m^2 (6.0 * 107 = 642 L/h). Five independent lines",
        "of evidence fix the divisor at 100, i.e. CL = intercept + slope *",
        "(CRCL/100), so the slope is the renal clearance arm at CRCL = 100.",
        "(1) Table 2 also prints a derived CL of 7.2 L/h for the study",
        "population; solving 0.86 + 6.0 * x = 7.2 gives x = 1.0567.",
        "(2) The SAME x = 1.0571 solves the companion tazobactam block",
        "(6.9 + 17.5 * x = 25.4), which was a separate NPAG fit -- a coincidence",
        "at two different drugs is implausible, and x ~ 1.06 pins the divisor at",
        "100 against the cohort median CLcr of 107 (Table 1) / 108 (Results).",
        "(3) Propagating the tabulated SDs through the same equation reproduces",
        "the tabulated SD of CL: sqrt(0.69^2 + (3.3*1.057)^2) = 3.6 against a",
        "printed 3.2 for ceftolozane, and sqrt(5.6^2 + (6.9*1.057)^2) = 9.2",
        "against a printed 9.4 for tazobactam.",
        "(4) The Discussion prints four simulated steady-state unbound ceftolozane",
        "concentrations from continuous-infusion regimens, at two stated CLcr",
        "values. Rate/CL with CL = 0.86 + 6.0*(CRCL/100) reproduces all four:",
        "4.5 g/24 h combination (= 3 g ceftolozane, 125 mg/h) gives 18.2 and 10.7",
        "mg/liter at CRCL 100 and 180 against a published 19 and 11.2; 9 g/24 h",
        "(= 6 g ceftolozane, 250 mg/h) gives 36.4 and 21.4 against a published 38",
        "and 22.4. All four are within 4.7%, and the published-to-predicted RATIO",
        "is 1.044-1.047 in every case.",
        "(5) That residual is the right SIGN and cannot be obtained any other way.",
        "The published values are MEANS of a 1,000-subject Monte Carlo, and",
        "E[Rate/CL] > Rate/E[CL] strictly by Jensen's inequality, so a correct",
        "divisor must predict BELOW the published mean. Divisors of 107 or 108",
        "(the cohort median, the other obvious candidates) predict 19.3-19.5 and",
        "11.4-11.5 -- ABOVE the published means, which is impossible. Only the",
        "divisor 100 has the admissible sign.",
        "RANGE FITTED. Table 1: median 107, IQR 74-145 mL/min/1.73 m^2. Patients",
        "with renal dysfunction requiring renal replacement therapy were EXCLUDED",
        "(Methods, 'Patients'), so this model carries no information about renal",
        "impairment and must not be used there -- the intercept is the only term",
        "left at CRCL = 0 and it was never informed by a low-clearance patient.",
        "The paper's own target-attainment simulations run the term out to 180",
        "mL/min/1.73 m^2, above the observed IQR but within the augmented-renal",
        "clearance range the study was designed around."
      ),
      source_name = "CL CRurinary"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Results: 'Total body weight (WT) was related to volume of distribution of",
        "the central compartment (V 1) for both ceftolozane (V 1 = V * WT/80) and",
        "tazobactam (V 1 = V * [WT/80] 0.75), where V is the typical value of the",
        "central volume of distribution.' For CEFTOLOZANE the exponent is",
        "structurally 1 -- the equation is printed as a bare proportionality with",
        "no exponent, and Table 2 carries no exponent row and therefore no SD or",
        "CV for it. Encoded as e_wt_vc <- fixed(1). The companion",
        "Sime_2019_tazobactam model uses 0.75 on the same column; the two",
        "exponents genuinely differ because the two analytes were fitted",
        "separately.",
        "REFERENCE 80 kg. The divisor is printed directly in the equation. It is",
        "the rounded cohort median: Table 1 gives 79.5 kg (IQR 64-99), and the",
        "Results restate it as 'body weight (80 kg) of the study population' when",
        "describing the simulation covariate values.",
        "Baseline (time-fixed); the study sampled within a single dosing interval.",
        "Body mass index was screened separately and not retained."
      ),
      source_name = "Wt"
    )
  )

  # Covariates the paper screened but did not retain. Methods, 'Population PK
  # modeling': 'Covariates selected for investigation include serum creatinine,
  # urinary creatinine clearance, body weight, body mass index, albumin
  # concentration, Acute Physiology and Chronic Health Evaluation II (APACHE II)
  # score, and Sequential Organ Failure Assessment (SOFA) score.' Of those, only
  # urinary creatinine clearance and body weight survived forward addition and
  # backward deletion.
  covariatesDataExcluded <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = paste(
        "Screened and not retained; the MEASURED urinary creatinine clearance was",
        "retained instead. Table 1 median 46 umol/L (IQR 39-77) -- note the SI",
        "unit, which is about 0.52 mg/dL at the median and confirms the cohort was",
        "selected for preserved-to-augmented renal function."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = paste(
        "Screened and not retained; total body weight was retained instead.",
        "Table 1 median 28.5 kg/m^2 (IQR 22.1-32.9)."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened and not retained. Table 1 median 25 g/L (IQR 19-28), i.e. the",
        "cohort was uniformly hypoalbuminaemic. The Introduction motivates albumin",
        "at length as a driver of unbound-drug distribution in the critically ill,",
        "but ceftolozane is bound at a low level and the analysis was run on",
        "directly-measured UNBOUND concentrations, which removes the protein-",
        "binding confounder the covariate would otherwise proxy."
      )
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score at ICU admission",
      units = "(score)",
      type = "continuous",
      notes = "Screened as a severity-of-illness covariate and not retained. Table 1 median 19.5 (IQR 16-26)."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "(score)",
      type = "continuous",
      notes = "Screened as a severity-of-illness covariate and not retained. Table 1 median 6 (IQR 3-8)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 12L,
    n_studies = 1L,
    n_concentrations = 133L,
    age_range = "median 56 years (IQR 52-61); enrolment required age >= 18 years",
    weight_range = "median 79.5 kg (IQR 64-99)",
    sex_female_pct = 66.7,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Critically ill adults in a quaternary referral intensive care unit with a",
      "systemic infection known or suspected to be caused by a bacterium",
      "susceptible to ceftolozane-tazobactam. All 12 had a positive culture",
      "(Table 1); lung was the commonest source (75%). Patients were EXCLUDED if",
      "they had renal dysfunction necessitating renal replacement therapy, a known",
      "or suspected cephalosporin allergy, piperacillin-tazobactam in the",
      "preceding 7 days, or were pregnant.",
      sep = " "
    ),
    renal_function = paste(
      "Preserved to augmented, by design. Measured urinary creatinine clearance",
      "median 107 mL/min/1.73 m^2 (IQR 74-145); serum creatinine median 46 umol/L",
      "(IQR 39-77). No patient on renal replacement therapy.",
      sep = " "
    ),
    dose_range = paste(
      "1.5 g or 3.0 g ceftolozane-tazobactam (fixed 2:1 ratio, i.e. 1000/500 mg or",
      "2000/1000 mg) every 8 h as a 1-hour intravenous infusion, at the discretion",
      "of the treating physician.",
      sep = " "
    ),
    regions = "Australia (single centre, Royal Brisbane and Women's Hospital)",
    bmi_range = "median 28.5 kg/m^2 (IQR 22.1-32.9)",
    severity = "APACHE II at admission median 19.5 (IQR 16-26); SOFA median 6 (IQR 3-8)",
    pkpd_target = paste(
      "The paper's Monte Carlo target-attainment analysis used 40%, 60% and 100%",
      "fT>MIC for ceftolozane (Methods, 'Dosing simulations'), against the EUCAST",
      "Pseudomonas aeruginosa MIC distribution with a clinical breakpoint of",
      "4 mg/liter. Those targets are a USE of this model, not part of it.",
      sep = " "
    ),
    notes = paste(
      "Baseline demographics from Table 1. Sampling was intensive but confined to",
      "a single dosing interval: pre-dose, 15 and 45 min after the start of the",
      "1-h infusion, at the end of the post-infusion line flush, at 2, 3, 4, 5, 6",
      "and 7 h, and immediately before the next dose. Unbound plasma was isolated",
      "by ultracentrifugation and assayed by UHPLC-MS/MS, calibrated 1-100",
      "mg/liter for ceftolozane.",
      "ESTIMATION was NONPARAMETRIC: Pmetrics 1.5.2 running the nonparametric",
      "adaptive grid (NPAG) algorithm in R, not NONMEM and not a parametric",
      "maximum-likelihood method. NPAG returns a discrete joint distribution over",
      "support points; Table 2 summarises its MARGINALS as mean, SD and CV. The",
      "joint density -- and therefore every parameter correlation -- is not",
      "recoverable from the publication. See the vignette Assumptions and",
      "deviations for what that costs.",
      "Twelve subjects is a small cohort and the authors say so: 'Another",
      "important limitation of this study is the small sample size, which offers a",
      "limited spread of covariates, limiting broad extrapolation of the study",
      "findings.'",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # STRUCTURAL PARAMETERS -- Sime 2019 Table 2, 'Ceftolozane' block,
    # 'Mean' column. Every value is the MEAN of the NPAG support-point
    # distribution, not a parametric point estimate with a standard error;
    # the accompanying SD column is the spread of that distribution across
    # subjects, i.e. it is inter-individual variability, NOT uncertainty.
    #
    # ENCODING OF THE MEAN. Each tabulated mean is carried below as the
    # MEDIAN of a log-normal marginal (log(x) is the mu of a log-normal
    # whose median is x), matching the two existing Pmetrics/NPAG models in
    # this library, Setiawan_2023_sulbactam.R and
    # Hughes_2024_vancomycin_nonparametric.R. The consequence is explicit:
    # a stochastic cohort simulated from this file has a population mean
    # exp(mu + omega^2/2) = (tabulated mean) * sqrt(1 + CV^2), which exceeds
    # the tabulated mean -- by 28% on the intercept, 14% on the slope, 2% on
    # V, 88% on Kcp and 37% on Kpc. The alternative encoding (shifting mu
    # down so the log-normal MEAN equals the tabulated mean) was rejected
    # because it would break the paper's own published arithmetic: Table 2's
    # derived CL of 7.2 L/h, and the four simulated steady-state
    # concentrations in the Discussion, are all computed FROM these
    # tabulated values as point estimates, and only the median encoding
    # reproduces them under rxode2::zeroRe(). See the vignette Assumptions
    # and deviations.
    # =====================================================================

    # --- Clearance: additive, linear in measured urinary CLcr ------------
    # Results: 'CL = intercept + slope * CL CRurinary'. Mapped onto the
    # registered additive multi-component clearance canonicals lcl_nonren
    # (the intercept, i.e. the clearance remaining at zero measured renal
    # function) and lcl_renal (the slope arm), following the precedent in
    # Kim_2016_tazobactam.R and Bulitta_2011_cefpirome.R. The paper itself
    # says only 'intercept' and 'slope' and makes no mechanistic claim; the
    # non-renal reading is the standard interpretation of that structure and
    # is what the canonical names encode. Note that NO patient in this
    # cohort had impaired renal function, so the intercept is an
    # extrapolation to a region the data never visited.
    lcl_nonren <- log(0.86)
    label("Clearance intercept, i.e. the arm independent of measured renal function (L/h)")
    # Table 2, Ceftolozane / Intercept: mean 0.86 (SD 0.69, CV 80%)
    lcl_renal <- log(6.0)
    label("Renal clearance arm at CRCL = 100 mL/min/1.73 m^2 (L/h)")
    # Table 2, Ceftolozane / Slope: mean 6.0 (SD 3.3, CV 54%). The /100
    # normalization is derived, not printed -- see covariateData$CRCL.
    # Check: 0.86 + 6.0 * (107/100) = 7.28 against Table 2's derived
    # 'CL (liters/h) 7.2', footnoted 'Value calculated for the study
    # population'.

    # --- Central volume ---------------------------------------------------
    lvc <- log(20.4)
    label("Central volume of distribution at WT = 80 kg (L)")
    # Table 2, Ceftolozane / V (liters): mean 20.4 (SD 3.7, CV 18%); also the
    # Abstract, 'V of 20.4 +/- 3.7'.
    e_wt_vc <- fixed(1)
    label("Power exponent of total body weight on central volume, WT/80 (unitless)")
    # Results: 'V 1 = V * WT/80' for ceftolozane -- printed as a bare
    # proportionality with no exponent and no entry in Table 2, so
    # structurally 1 rather than estimated. Contrast the companion
    # tazobactam model, where the same paper prints (WT/80)^0.75.

    # --- Intercompartmental transfer, as rate constants -------------------
    # The paper parameterises distribution directly as Kcp and Kpc rather
    # than as Q and Vp, so those are the primary parameters here and carry
    # the IIV; q and vp are derived inside model(). Mapped onto the
    # canonical micro-constant names k12 / k21, following
    # Setiawan_2023_sulbactam.R (same research group, same Pmetrics
    # parameterisation).
    lk12 <- log(0.46)
    label("Transfer rate constant central -> peripheral1, Kcp (1/h)")
    # Table 2, Ceftolozane / Kcp (h-1): mean 0.46 (SD 0.74, CV 159%)
    lk21 <- log(0.39)
    label("Transfer rate constant peripheral1 -> central, Kpc (1/h)")
    # Table 2, Ceftolozane / Kpc (h-1): mean 0.39 (SD 0.37, CV 94%)

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Sime 2019 Table 2, 'CV (%)' column.
    #
    # SCALE. The printed CV is the ordinary descriptive coefficient of
    # variation of the NPAG support-point distribution, SD/mean: 3.7/20.4 =
    # 18.1% (printed 18), 0.74/0.46 = 161% (printed 159), 0.37/0.39 = 94.9%
    # (printed 94), 0.69/0.86 = 80.2% (printed 80), 3.3/6.0 = 55% (printed
    # 54). Every row reproduces, so the column is a CV and not an omega, and
    # the log-normal variance is the exact conversion omega^2 =
    # log(CV^2 + 1) -- the same convention as
    # Hughes_2024_vancomycin_nonparametric.R.
    #
    # The NPAG joint density is NOT recoverable from the publication, so
    # these marginals are encoded as INDEPENDENT. That is certainly wrong
    # for the intercept/slope pair, which trade off against one another in
    # any linear fit; the paper's own simulated steady-state concentrations
    # have a CV near 29% where independent marginals imply about 49%. See
    # the vignette Assumptions and deviations.
    # =====================================================================
    etalcl_nonren ~ 0.494696 # Table 2 Intercept CV 80%: log(0.80^2 + 1)
    etalcl_renal ~ 0.255882 # Table 2 Slope CV 54%: log(0.54^2 + 1)
    etalvc ~ 0.031886 # Table 2 V CV 18%: log(0.18^2 + 1)
    etalk12 ~ 1.260759 # Table 2 Kcp CV 159%: log(1.59^2 + 1)
    etalk21 ~ 0.633185 # Table 2 Kpc CV 94%: log(0.94^2 + 1)

    # =====================================================================
    # RESIDUAL UNEXPLAINED VARIABILITY is NOT reported anywhere in the
    # paper. Methods, 'Population PK modeling', describes only the menu of
    # forms that was tested -- 'the additive error mode was given by the
    # equation Error = (SD2 + lambda2)0.5, and the multiplicative mode was
    # given by the equation Error = SD * gamma' plus an assay-error
    # polynomial 'Error = C 0 + C 1 * obs, where the coefficients C 0 and
    # C 1 were optimized interactively'. Neither WHICH form was selected nor
    # ANY of lambda, gamma, C0 or C1 is printed, and no supplement was
    # deposited. Carried as fixed(0) rather than invented, following
    # Setiawan_2023_sulbactam.R. See the vignette Errata.
    # =====================================================================
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; not reported in the source)")
    addSd <- fixed(0)
    label("Additive residual SD (mg/L; not reported in the source)")
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters.
    #
    # Clearance (Results): CL = intercept + slope * CLcr_urinary, with the
    # slope arm normalized by 100 mL/min/1.73 m^2 (derivation in
    # covariateData$CRCL):
    #
    #   CL (L/h) = 0.86 + 6.0 * (CRCL / 100)
    #
    # Each arm carries its own eta because Table 2 reports a separate SD
    # and CV for the intercept and for the slope.
    # ------------------------------------------------------------------
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (CRCL / 100)
    cl <- cl_nonren + cl_renal

    # Central volume (Results): V1 = V * WT/80, exponent structurally 1.
    vc <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc

    # Primary distribution rate constants (Table 2 Kcp / Kpc).
    kcp <- exp(lk12 + etalk12)
    kpc <- exp(lk21 + etalk21)

    # ------------------------------------------------------------------
    # 2. Intercompartmental clearance and peripheral volume, derived from
    #    the primary rate constants. The ODEs below MUST be driven by
    #    q / vp rather than by kcp / kpc directly: rxSolve() defaults to
    #    useLinCmt = TRUE and, when the peripheral transfer is written
    #    straight from stored micro-constants, that rewrite can silently
    #    drop peripheral1 and solve a one-compartment model. Routing
    #    through q and vp keeps the closed-form and ODE solvers in
    #    agreement (the vignette asserts this identity). Same treatment as
    #    Setiawan_2023_sulbactam.R.
    # ------------------------------------------------------------------
    q <- kcp * vc
    vp <- q / kpc

    # ------------------------------------------------------------------
    # 3. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 4. Two-compartment disposition. Ceftolozane-tazobactam was given as a
    #    1-hour intravenous infusion into the central compartment; the
    #    paper's simulations also cover 4-hour extended infusions and
    #    24-hour continuous infusions. The infusion is expressed in the
    #    event table (amt + dur or rate), not here.
    # ------------------------------------------------------------------
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ------------------------------------------------------------------
    # 5. Observation. Dose in mg, volumes in L -> mg/L, the unit the paper
    #    reports. Cc is the UNBOUND plasma concentration: the assay
    #    measured unbound drug directly (see compartmentData above), so no
    #    protein-binding conversion is applied anywhere in this file. It is
    #    the quantity the paper's fT>MIC targets are defined on.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
