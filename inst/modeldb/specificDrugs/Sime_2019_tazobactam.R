Sime_2019_tazobactam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for UNBOUND tazobactam in",
    "twelve critically ill adults without renal dysfunction, admitted to a",
    "quaternary referral intensive care unit in Brisbane, Australia. Fitted",
    "non-parametrically with the NPAG algorithm in Pmetrics 1.5.2 to",
    "directly-measured unbound plasma concentrations. Clearance is ADDITIVE in a",
    "covariate-free intercept and an arm linear in measured urinary creatinine",
    "clearance (paper: CL = intercept + slope * CLcr_urinary), and central volume",
    "scales allometrically with total body weight at a fixed 0.75 exponent",
    "(V1 = V * [WT/80]^0.75). The intercompartmental transfer is parameterised",
    "directly as the rate constants Kcp and Kpc rather than as Q and Vp. Every",
    "parameter carries its own inter-individual variability, taken from the mean",
    "and SD of the NPAG support-point distribution; the Kcp distribution is",
    "extremely wide (CV 293%). Residual unexplained variability is carried as",
    "fixed(0) because neither the selected Pmetrics error model nor its assay",
    "error polynomial coefficients were published. Ceftolozane and tazobactam were",
    "fitted in two separate NPAG runs, with different body-weight exponents on",
    "volume, and are supplied as two separate model files; see",
    "modellib('Sime_2019_ceftolozane') for the partner component of the fixed 2:1",
    "ceftolozane-tazobactam combination.",
    sep = " "
  )
  reference <- paste(
    "Sime FB, Lassig-Smith M, Starr T, Stuart J, Pandey S, Parker SL, Wallis SC,",
    "Lipman J, Roberts JA. Population pharmacokinetics of unbound ceftolozane and",
    "tazobactam in critically ill patients without renal dysfunction.",
    "Antimicrob Agents Chemother. 2019;63(10):e01265-19.",
    "doi:10.1128/AAC.01265-19. PMCID: PMC6761554.",
    "All structural and variability estimates are Table 2, 'Tazobactam' block.",
    "The covariate equations are the Results narrative ('CL = intercept + slope *",
    "CL CRurinary' and 'V 1 = V * [WT/80] 0.75'). No supplement was deposited with",
    "the article.",
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
  # assay validation reports tazobactam unbound fractions of 89%, 91% and 92% at
  # total concentrations of 80, 10 and 1.5 mg/liter.
  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
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
        "creatinine clearance.",
        "FUNCTIONAL FORM. Results: 'ceftolozane and tazobactam clearance linearly",
        "increased with an increase in urinary creatinine clearance. The final",
        "covariate model for clearance of both ceftolozane and tazobactam was",
        "expressed as CL = intercept + slope * CL CRurinary'. Additive, linear,",
        "NOT a power function -- the same form, and the same column, as the",
        "companion Sime_2019_ceftolozane model.",
        "NORMALIZATION BY 100 -- DERIVED, NOT PRINTED. The paper never states the",
        "scale on which CLcr enters, and the printed slope is not usable at face",
        "value in mL/min/1.73 m^2 (17.5 * 107 = 1,873 L/h). Solving Table 2's own",
        "derived value, 6.9 + 17.5 * x = 25.4, gives x = 1.0571; the SAME x =",
        "1.0567 solves the separately-fitted ceftolozane block",
        "(0.86 + 6.0 * x = 7.2), which pins the divisor at 100 against the cohort",
        "median CLcr of 107 (Table 1) / 108 (Results). Propagating the tabulated",
        "SDs through the same equation corroborates it:",
        "sqrt(5.6^2 + (6.9*1.057)^2) = 9.2 against a printed SD of 9.4. The full",
        "five-part derivation, including the Jensen's-inequality argument that",
        "excludes divisors of 107 and 108, is written out in the companion file",
        "Sime_2019_ceftolozane.R, whose Discussion simulations supply the",
        "independent check (the paper prints no simulated tazobactam",
        "concentrations, only that all regimens met the tazobactam target).",
        "RANGE FITTED. Table 1: median 107, IQR 74-145 mL/min/1.73 m^2. Patients",
        "with renal dysfunction requiring renal replacement therapy were EXCLUDED",
        "(Methods, 'Patients'), so this model carries no information about renal",
        "impairment and must not be used there -- the intercept is the only term",
        "left at CRCL = 0 and it was never informed by a low-clearance patient.",
        "The paper's own target-attainment simulations run the term out to 180",
        "mL/min/1.73 m^2."
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
        "central volume of distribution.' For TAZOBACTAM the exponent is the",
        "canonical allometric 0.75, printed in the equation itself. Table 2",
        "carries no exponent row and therefore no SD or CV for it, so it is a",
        "structural constant rather than an estimate; encoded as",
        "e_wt_vc <- fixed(0.75).",
        "THE TWO ANALYTES GENUINELY DIFFER. The companion Sime_2019_ceftolozane",
        "model uses a bare proportionality (exponent 1) on the same column,",
        "because the two drugs were fitted in separate NPAG runs. Do not",
        "harmonise the two exponents.",
        "REFERENCE 80 kg. The divisor is printed directly in the equation. It is",
        "the rounded cohort median: Table 1 gives 79.5 kg (IQR 64-99), and the",
        "Results restate it as 'body weight (80 kg) of the study population'.",
        "Baseline (time-fixed); the study sampled within a single dosing interval."
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
        "unit, which is about 0.52 mg/dL at the median."
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
        "cohort was uniformly hypoalbuminaemic. Because the analysis was run on",
        "directly-measured UNBOUND concentrations, the protein-binding confounder",
        "that albumin would otherwise proxy is already removed."
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
    # 133 unbound concentration-time data points entered the population PK
    # analysis (Results); the paper does not split that total between the two
    # analytes, so the same figure is carried in both model files.
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
      "of the treating physician. The TAZOBACTAM dose is therefore 500 or 1000 mg.",
      sep = " "
    ),
    regions = "Australia (single centre, Royal Brisbane and Women's Hospital)",
    bmi_range = "median 28.5 kg/m^2 (IQR 22.1-32.9)",
    severity = "APACHE II at admission median 19.5 (IQR 16-26); SOFA median 6 (IQR 3-8)",
    pkpd_target = paste(
      "For tazobactam the paper used 20% fT>1 mg/liter (Methods, 'Dosing",
      "simulations'), a time-above-a-fixed-minimum-effective-concentration target",
      "rather than an MIC-referenced one. Results: 'for tazobactam, all simulated",
      "dosing regimens had a 100% probability of achieving the recommended",
      "target'. That target is a USE of this model, not part of it.",
      sep = " "
    ),
    notes = paste(
      "Baseline demographics from Table 1. Sampling was intensive but confined to",
      "a single dosing interval: pre-dose, 15 and 45 min after the start of the",
      "1-h infusion, at the end of the post-infusion line flush, at 2, 3, 4, 5, 6",
      "and 7 h, and immediately before the next dose. Unbound plasma was isolated",
      "by ultracentrifugation and assayed by UHPLC-MS/MS, calibrated 0.5-100",
      "mg/liter for tazobactam.",
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
    # STRUCTURAL PARAMETERS -- Sime 2019 Table 2, 'Tazobactam' block, 'Mean'
    # column. Every value is the MEAN of the NPAG support-point
    # distribution, not a parametric point estimate with a standard error;
    # the accompanying SD column is the spread of that distribution across
    # subjects, i.e. it is inter-individual variability, NOT uncertainty.
    #
    # ENCODING OF THE MEAN. Each tabulated mean is carried below as the
    # MEDIAN of a log-normal marginal, matching the two existing
    # Pmetrics/NPAG models in this library, Setiawan_2023_sulbactam.R and
    # Hughes_2024_vancomycin_nonparametric.R, and matching the companion
    # Sime_2019_ceftolozane.R. The consequence is explicit: a stochastic
    # cohort simulated from this file has a population mean
    # (tabulated mean) * sqrt(1 + CV^2), which exceeds the tabulated mean --
    # by 29% on the intercept, 8% on the slope, 5% on V, 5% on Kpc and by a
    # factor of 3.1 on Kcp, whose CV is 293%. See the vignette Assumptions
    # and deviations, which quantifies what that does and does not affect
    # (it does not affect AUC or steady-state concentration, which depend on
    # clearance alone).
    # =====================================================================

    # --- Clearance: additive, linear in measured urinary CLcr ------------
    # Results: 'CL = intercept + slope * CL CRurinary'. Mapped onto the
    # registered additive multi-component clearance canonicals lcl_nonren
    # (the intercept, i.e. the clearance remaining at zero measured renal
    # function) and lcl_renal (the slope arm), following the precedent in
    # Kim_2016_tazobactam.R and Bulitta_2011_cefpirome.R. The paper itself
    # says only 'intercept' and 'slope'. Note that NO patient in this cohort
    # had impaired renal function, so the intercept is an extrapolation to a
    # region the data never visited -- and for tazobactam the intercept is a
    # much larger share of total clearance than for ceftolozane (6.9 of 25.4
    # L/h, 27%, against 0.86 of 7.2 L/h, 12%), which makes that
    # extrapolation correspondingly more load-bearing.
    lcl_nonren <- log(6.9)
    label("Clearance intercept, i.e. the arm independent of measured renal function (L/h)")
    # Table 2, Tazobactam / Intercept: mean 6.9 (SD 5.6, CV 81%)
    lcl_renal <- log(17.5)
    label("Renal clearance arm at CRCL = 100 mL/min/1.73 m^2 (L/h)")
    # Table 2, Tazobactam / Slope: mean 17.5 (SD 6.9, CV 40%). The /100
    # normalization is derived, not printed -- see covariateData$CRCL.
    # Check: 6.9 + 17.5 * (107/100) = 25.6 against Table 2's derived
    # 'CL (liters/h) 25.4', footnoted 'Value calculated for the study
    # population'.

    # --- Central volume ---------------------------------------------------
    lvc <- log(32.4)
    label("Central volume of distribution at WT = 80 kg (L)")
    # Table 2, Tazobactam / V (liters): mean 32.4 (SD 10, CV 31%); also the
    # Abstract, 'V of ... 32.4 +/- 10'.
    e_wt_vc <- fixed(0.75)
    label("Power exponent of total body weight on central volume, WT/80 (unitless)")
    # Results: 'V 1 = V * [WT/80] 0.75' for tazobactam -- the canonical
    # allometric exponent, printed in the equation and absent from Table 2,
    # so structural rather than estimated.

    # --- Intercompartmental transfer, as rate constants -------------------
    # The paper parameterises distribution directly as Kcp and Kpc rather
    # than as Q and Vp, so those are the primary parameters here and carry
    # the IIV; q and vp are derived inside model(). Mapped onto the
    # canonical micro-constant names k12 / k21, following
    # Setiawan_2023_sulbactam.R (same research group, same Pmetrics
    # parameterisation).
    lk12 <- log(2.96)
    label("Transfer rate constant central -> peripheral1, Kcp (1/h)")
    # Table 2, Tazobactam / Kcp (h-1): mean 2.96 (SD 8.69, CV 293%). The
    # Abstract rounds the SD to 8.6; Table 2's 8.69 is used, and it is the
    # value consistent with the printed CV (8.69/2.96 = 294%, printed 293;
    # 8.6/2.96 = 291%).
    lk21 <- log(26.5)
    label("Transfer rate constant peripheral1 -> central, Kpc (1/h)")
    # Table 2, Tazobactam / Kpc (h-1): mean 26.5 (SD 8.4, CV 32%)
    #
    # Kpc >> Kcp means the peripheral compartment equilibrates very fast and
    # holds little drug: at the typical values the derived peripheral volume
    # is vc * Kcp/Kpc = 32.4 * 2.96/26.5 = 3.6 L. Tazobactam disposition is
    # therefore close to one-compartment, which is consistent with Kcp being
    # the least well determined parameter in the table.

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Sime 2019 Table 2, 'CV (%)' column.
    #
    # SCALE. The printed CV is the ordinary descriptive coefficient of
    # variation of the NPAG support-point distribution, SD/mean: 5.6/6.9 =
    # 81.2% (printed 81), 6.9/17.5 = 39.4% (printed 40), 10/32.4 = 30.9%
    # (printed 31), 8.69/2.96 = 294% (printed 293), 8.4/26.5 = 31.7%
    # (printed 32). Every row reproduces, so the column is a CV and not an
    # omega, and the log-normal variance is the exact conversion
    # omega^2 = log(CV^2 + 1) -- the same convention as
    # Hughes_2024_vancomycin_nonparametric.R.
    #
    # The NPAG joint density is NOT recoverable from the publication, so
    # these marginals are encoded as INDEPENDENT. That is certainly wrong
    # for the intercept/slope pair, which trade off against one another in
    # any linear fit. See the vignette Assumptions and deviations.
    # =====================================================================
    etalcl_nonren ~ 0.504465 # Table 2 Intercept CV 81%: log(0.81^2 + 1)
    etalcl_renal ~ 0.148420 # Table 2 Slope CV 40%: log(0.40^2 + 1)
    etalvc ~ 0.091758 # Table 2 V CV 31%: log(0.31^2 + 1)
    etalk12 ~ 2.260189 # Table 2 Kcp CV 293%: log(2.93^2 + 1)
    etalk21 ~ 0.097490 # Table 2 Kpc CV 32%: log(0.32^2 + 1)

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
    #   CL (L/h) = 6.9 + 17.5 * (CRCL / 100)
    #
    # Each arm carries its own eta because Table 2 reports a separate SD
    # and CV for the intercept and for the slope.
    # ------------------------------------------------------------------
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (CRCL / 100)
    cl <- cl_nonren + cl_renal

    # Central volume (Results): V1 = V * (WT/80)^0.75.
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
    #    the quantity the paper's 20% fT>1 mg/liter target is defined on.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
