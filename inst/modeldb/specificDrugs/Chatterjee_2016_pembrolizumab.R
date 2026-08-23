Chatterjee_2016_pembrolizumab <- function() {
  description <- paste(
    "Exposure-response tumor-size (sum of longest diameters, SLD) model for",
    "pembrolizumab in previously treated and treatment-naive advanced",
    "non-small-cell lung cancer (NSCLC), developed by Chatterjee et al.",
    "(Merck) on the KEYNOTE-001 NSCLC expansion cohorts (n = 496 with both",
    "tumor-size and pharmacokinetic data at 2 mg/kg Q3W, 10 mg/kg Q3W and",
    "10 mg/kg Q2W). The structural model splits the baseline SLD into a",
    "treatment-sensitive fraction f that decays first-order at kdeath and a",
    "resistant fraction (1 - f) that grows first-order at kgrowth, giving",
    "the published bi-exponential form",
    "SLD(t) = Baseline * [(1 - f) * exp(kgrowth * t) +",
    "f * exp(-kdeath * max(0, t - delay))].",
    "Pembrolizumab exposure enters as a log-linear (power) effect of the",
    "steady-state 6-week AUC on kdeath, normalized to the population-typical",
    "AUCss-6weeks of 7079 mg*day/L. The final covariate model adds PD-L1",
    "tumor proportion score (four levels) on kdeath and EGFR mutation status",
    "(three levels) on the logit of f. The estimated exposure effect is not",
    "statistically significant (95% CI -0.0784 to 0.47, P = 0.54) and was",
    "retained by the authors only so that the magnitude of any potential",
    "exposure-response relationship could be simulated; the paper's",
    "conclusion is that response is flat over 2-10 mg/kg. There is no PK",
    "input: exposure is supplied per subject as the covariate AUC_PEMBRO,",
    "which the source analysis obtained as dose/CL from the companion",
    "pembrolizumab population-PK model (Ahamadi 2017; available in this",
    "library as Ahamadi_2017_pembrolizumab)."
  )
  reference <- paste(
    "Chatterjee M, Turner DC, Felip E, Lena H, Cappuzzo F, Horn L,",
    "Garon EB, Hui R, Arkenau H-T, Gubens MA, Hellmann MD, Dong D, Li C,",
    "Mayawala K, Freshwater T, Ahamadi M, Stone J, Lubiniecki GM, Zhang J,",
    "Im E, De Alwis DP, Kondic AG, Flotten O.",
    "Systematic evaluation of pembrolizumab dosing in patients with",
    "advanced non-small-cell lung cancer.",
    "Ann Oncol. 2016;27(7):1291-1298. doi:10.1093/annonc/mdw174.",
    "PMID: 27117531.",
    "Structural equation and exposure-effect equation from the main-article",
    "Methods ('tumor size NLME model structure' and 'exposure-efficacy",
    "analysis'); all final parameter values from supplementary Table S6;",
    "covariate parameterization from the Supplementary Methods",
    "('Handling of Covariates').",
    sep = " "
  )
  vignette <- "Chatterjee_2016_pembrolizumab"

  units <- list(
    time          = "day",
    dosing        = "n/a (no PK input; pembrolizumab exposure enters as the per-subject covariate AUC_PEMBRO in mg*day/L)",
    concentration = "mm (the observable `TS` is the RECIST 1.1 sum of the longest diameters of target lesions)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are sub-fractions of the measured tumor
  # diameter (mm), not amounts, so `units` is NA. verified = TRUE: the
  # sensitive / resistant split is stated explicitly in supplementary
  # Figure S1A ("only a fraction (f) of total tumor diameter is accessible
  # and/or sensitive to treatment, which permits the remaining portion
  # (1-f) to undergo unimpeded exponential growth").
  compartmentData <- list(
    growth = list(analyte = "tumor-size", units = NA_character_, specimen = "tumor", verified = TRUE),
    shrink = list(analyte = "tumor-size", units = NA_character_, specimen = "tumor", verified = TRUE)
  )

  covariateData <- list(
    TUM_SLD = list(
      description        = "Observed baseline sum of the longest diameters of target lesions per RECIST 1.1, measured at initial screening.",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "Baseline",
      notes              = paste(
        "The paper's 'Baseline' term: 'the actual measured tumor size (SLD) at initial screening' (Methods, tumor size NLME model structure).",
        "Supplementary Methods, Structural Model Selection: 'Fixing baseline tumor size to observed values was found to improve model stability', so the baseline is a per-subject regressor and is NOT estimated.",
        "It initialises both sub-states, growth(0) = (1 - f) * TUM_SLD and shrink(0) = f * TUM_SLD, so TS(0) = TUM_SLD exactly.",
        "Cohort median 91.70 mm, range 10.40-548.30 mm (supplementary Table S4, N = 505).",
        "Baseline tumor size was also tested as a covariate ON kdeath during the stepwise covariate search (supplementary Table S5) and was removed in the backward-elimination step, so it carries no covariate coefficient in the final model.",
        sep = " "
      )
    ),
    AUC_PEMBRO = list(
      description        = "Per-subject pembrolizumab area under the serum concentration-time curve at steady state over a 6-week interval (AUCss-6weeks).",
      units              = "mg*day/L (equivalently ug*day/mL; the paper prints 'mg/l x day' in Methods and 'ug.day/ml' on the Figure 2 axis -- numerically identical)",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "AUCss-6weeks",
      notes              = paste(
        "Log-linear (power) effect on kdeath: kdeath = TVkdeath * (AUC_PEMBRO / 7079)^theta, per the Methods exposure-effect equation.",
        "The normalizing constant AUCtypical,ss-6weeks = 7079 mg*day/L is the population-typical exposure stated in the Methods and is hard-coded in model() as a structural constant of the published equation.",
        "The source analysis did NOT model pembrolizumab PK here: 'Results from the independent population pharmacokinetics model provided post hoc clearance (CL) estimates, with plasma exposure within the dosing interval at steady state calculated as dose/CL'. That upstream popPK model was 'submitted for publication' when this paper appeared and was published as Ahamadi 2017 (CPT Pharmacometrics Syst Pharmacol 6:49-57), which is packaged in this library as Ahamadi_2017_pembrolizumab.",
        "To reproduce the paper's exposure metric, compute AUC_PEMBRO = (total mg administered over 6 weeks) / CL, with CL the individual clearance in L/day: 2 doses for a Q3W regimen and 3 doses for a Q2W regimen.",
        "AUCss-6weeks was chosen over a time-dependent exposure metric to avoid confounding by early dropout (Methods).",
        sep = " "
      )
    ),
    PDL1_TUM = list(
      description        = "Baseline tumor PD-L1 expression by immunohistochemistry (22C3 clone), reported as the Tumor Proportion Score (TPS): the percentage of tumor cells with membranous PD-L1 staining.",
      units              = "percent (0-100)",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "PD-L1",
      notes              = paste(
        "Carried as the canonical continuous TPS column; model() derives the paper's three non-reference indicator levels from it.",
        "Category boundaries per the Supplementary Methods, Handling of Covariates: TPS >= 50% ('strongly positive'), TPS 1%-49% ('weakly positive'), TPS < 1% ('negative').",
        "TPS 1%-49% is the MOST FREQUENT category (201 of 505, supplementary Table S3) and is therefore the model's reference level, receiving no covariate term.",
        "Subjects whose PD-L1 status was not assignable are carried in the companion column PDL1_TUM_MISSING and must have PDL1_TUM set so that neither derived indicator fires; model() gates both indicators on (1 - PDL1_TUM_MISSING) so any placeholder value is safe.",
        "Assay: prototype assay for enrolment, clinical trial assay (22C3 anti-human PD-1 antibody) for the retrospective TPS used in the model (Methods, treatment and assessments).",
        sep = " "
      )
    ),
    PDL1_TUM_MISSING = list(
      description        = "Binary indicator: 1 = the baseline PD-L1 tumor proportion score could not be assigned for this subject, 0 = a TPS value is available.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (PD-L1 TPS available)",
      source_name        = "PD-L1 'Unknown'",
      notes              = paste(
        "The source analysis retained 'unknown' PD-L1 as a fourth modelled category with its own coefficient (PD-L1_3, supplementary Table S6) rather than imputing it into a measured level; 60 of 505 subjects were unknown (supplementary Table S3).",
        "Supplementary text, Final Exposure-Response Model: 'the only final model parameter that was identified with poor precision was PD-L1_3 on kdeath ... The large uncertainty in this parameter likely reflects that the unknown PD-L1 category includes patients with PD-L1-positive and negative tumors'. RSE 145%; interpret the coefficient with caution.",
        "Mutually exclusive with the derived TPS >= 50% and TPS < 1% indicators.",
        sep = " "
      )
    ),
    TUM_EGFR_MUT = list(
      description        = "Binary indicator of tumor EGFR (epidermal growth factor receptor) mutation status: 1 = EGFR-mutant tumor, 0 = EGFR wild-type tumor.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EGFR wild type, the most frequent category)",
      source_name        = "EGFR status",
      notes              = paste(
        "Additive effect on the logit of the responsive tumor fraction f (EGFR_1, supplementary Table S6).",
        "EGFR wild type is the most frequent category (409 of 505; supplementary Table S3) and is therefore the reference level.",
        "Supplementary text, Covariate Effects: 'EGFR mutation was associated with a lower fraction of the tumor that responds to treatment ... median f was 3.2-fold higher in patients with EGFR wild-type versus mutant tumors'.",
        "Set TUM_EGFR_MUT = 0 when TUM_EGFR_MUT_MISSING = 1; model() gates the mutant term so the unknown stratum receives only its own coefficient.",
        sep = " "
      )
    ),
    TUM_EGFR_MUT_MISSING = list(
      description        = "Binary indicator: 1 = tumor EGFR mutation status was not determined for this subject, 0 = EGFR status (wild type or mutant) is known.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EGFR status known)",
      source_name        = "EGFR 'Unknown'",
      notes              = paste(
        "The source analysis retained 'unknown' EGFR as a third modelled category with its own coefficient on logit(f) (EGFR_2, supplementary Table S6) rather than imputing it into the wild-type reference; 26 of 505 subjects were unknown (supplementary Table S3).",
        "The estimated coefficient is positive (+1.66), i.e. EGFR-unknown subjects had a HIGHER responsive fraction than the wild-type reference; as with the PD-L1 unknown stratum this is a missingness artefact, not a biological effect, and should not be extrapolated.",
        "Mutually exclusive with TUM_EGFR_MUT = 1.",
        sep = " "
      )
    ),
    T_SCAN_TO_DOSE = list(
      description        = "Per-subject delay between the baseline tumor scan (the model's time origin) and the first pembrolizumab dose.",
      units              = "day",
      type               = "continuous",
      reference_category = NULL,
      source_name        = "delay",
      notes              = paste(
        "The paper's 'delay' term in the structural equation: 'the delay between baseline and the first dose' (Methods, tumor size NLME model structure).",
        "Supplementary Methods, Structural Model Selection: an ESTIMATED pharmacologic delay in drug action was tested and dropped (frequency < 5%, high shrinkage), but 'A delay between the baseline scan and the first pembrolizumab dose was retained in the model as a fixed individual-specific parameter' -- i.e. it is per-subject DATA, not an estimated parameter, and so is carried here as a covariate.",
        "Enters the model only through max(0, t - T_SCAN_TO_DOSE): the sensitive sub-state is held at its initial value until the first dose and decays at kdeath thereafter. The resistant sub-state grows from t = 0 regardless.",
        "The source paper does not report the distribution of this delay. Set T_SCAN_TO_DOSE = 0 for a subject dosed on the day of the baseline scan; the vignette uses 0 and states the assumption.",
        sep = " "
      )
    )
  )

  # Covariates that the source paper's stepwise covariate search screened but
  # did NOT retain in the final model. Documented so the provenance of the
  # covariate screen survives without creating declared-but-unreferenced
  # covariates. Sources: main-article Methods ("The patient-specific factors
  # of PD-L1 expression level, smoking history, ECOG performance status,
  # demographics (age, sex, and weight), baseline tumor size, prior treatment,
  # and EGFR mutation status were tested for inclusion"), supplementary
  # Table S5 (forward-addition / backward-elimination log) and the
  # Supplementary Methods, Covariate Effects.
  covariatesDataExcluded <- list(
    SMOKER = list(
      description = "Smoking history (former or current smoker vs nonsmoker).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested during the stepwise covariate search and not selected (supplementary Methods, Covariate Effects: 'covariates such as number of lines of prior treatment, smoking history, baseline tumor size, or ECOG performance status ... were not selected'). Cohort: 382 former/current smokers, 121 nonsmokers, 2 unknown (supplementary Table S3)."
    ),
    WHO_PS = list(
      description = "ECOG (Eastern Cooperative Oncology Group) performance status, 0 or 1.",
      units       = "(integer score)",
      type        = "categorical",
      notes       = "Tested and not selected (supplementary Methods, Covariate Effects). Cohort: ECOG 0 in 177, ECOG 1 in 325, unknown in 3 (supplementary Table S3). Eligibility restricted enrolment to ECOG 0-1."
    ),
    AGE = list(
      description = "Age at baseline.",
      units       = "year",
      type        = "continuous",
      notes       = "Tested as part of the demographics block and not selected (main-article Methods). Cohort median 64 years, range 32-93 (supplementary Table S4)."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as part of the demographics block and not selected (main-article Methods). Cohort: 267 male, 238 female (supplementary Table S3)."
    ),
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested as part of the demographics block and not selected (main-article Methods). Cohort median 70.00 kg, range 35.70-132.00, 1 missing (supplementary Table S4). Note that weight DOES enter the model indirectly, through the mg/kg dose that determines AUC_PEMBRO."
    ),
    PRIOR_THERAPY = list(
      description = "Any prior systemic therapy for advanced NSCLC (yes / no), and the number of prior lines.",
      units       = "(binary / count)",
      type        = "categorical",
      notes       = "Tested and not selected (supplementary Methods, Covariate Effects, explicitly names 'number of lines of prior treatment'). Cohort: 418 previously treated, 87 treatment-naive (supplementary Table S3)."
    ),
    RACE = list(
      description = "Race, grouped as white vs Asian vs other for testing.",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Regrouped for testability and not retained: supplementary Methods, Covariate Effects -- 'because the majority of patients were white and there were <10 patients in 3 of the 6 race categories, races were grouped as white vs Asian vs other to allow for testing'. Cohort: 418 white, 63 Asian, 18 Black or African American, 3 multiracial, 1 American Indian or Alaskan native, 1 Native Hawaiian or other Pacific Islander, 1 unknown (supplementary Table S3)."
    ),
    ALK_TRANSLOC = list(
      description = "ALK gene translocation / rearrangement status.",
      units       = "(binary)",
      type        = "binary",
      notes       = "NOT formally tested: supplementary Methods, Covariate Effects -- 'Because <2% of patients in the dataset harbored an ALK translocation, ALK status was not formally tested as a covariate.' Cohort: 435 wild type, 8 mutation or translocation, 62 unknown (supplementary Table S3)."
    )
  )

  population <- list(
    species         = "human (adults with locally advanced or metastatic NSCLC)",
    n_subjects      = 496L,
    n_studies       = 1L,
    age_range       = "32-93 years (supplementary Table S4, N = 505 with measurable baseline disease)",
    age_median      = "64 years",
    weight_range    = "35.70-132.00 kg",
    weight_median   = "70.00 kg",
    sex_female_pct  = 47.1,
    race_ethnicity  = c(White = 82.8, Asian = 12.5, Black = 3.6, Other = 1.2),
    disease_state   = "locally advanced or metastatic non-small-cell lung cancer, ECOG performance status 0-1, PD-L1 positive by the prototype assay for the final cohort; 83% previously treated",
    dose_range      = "pembrolizumab 2 mg/kg IV Q3W, 10 mg/kg IV Q3W, or 10 mg/kg IV Q2W (not a model input; enters only through AUC_PEMBRO)",
    regions         = "multinational KEYNOTE-001 (NCT01295827), phase Ib, multicenter open-label",
    notes           = paste(
      "Modeling population (main-article Methods, exposure-efficacy analysis): n = 496 patients with both PK data and measurable disease per RECIST v1.1 by central review at baseline --",
      "  6 treatment-naive + 47 previously treated at 2 mg/kg Q3W",
      "  45 treatment-naive + 216 previously treated at 10 mg/kg Q3W",
      "  39 treatment-naive + 143 previously treated at 10 mg/kg Q2W",
      "The covariate summary tables (supplementary Tables S3 and S4) are reported for the slightly larger N = 505 with measurable baseline tumor size; the demographics above are from those tables.",
      "Baseline PD-L1 tumor proportion score distribution (supplementary Table S3): TPS >= 50% 153 (30%), TPS 1%-49% 201 (40%), TPS < 1% 91 (18%), unknown 60 (12%).",
      "Baseline EGFR status (supplementary Table S3): wild type 409, mutant 70, unknown 26.",
      "Baseline sum of longest diameters (supplementary Table S4): median 91.70 mm, range 10.40-548.30 mm.",
      "Population-typical exposure AUCtypical,ss-6weeks = 7079 mg*day/L (main-article Methods).",
      "Tumor assessments by CT or MRI at baseline and every 9 weeks thereafter; SLD by RECIST v1.1 per independent central review.",
      "Estimation used PsN's stepwise covariate modeling (forward inclusion P < 0.01, backward exclusion P < 0.001).",
      "The safety analysis in the same paper (logistic regression of immune-mediated adverse events on AUCss-6weeks, and a time-to-event analysis) reports only P values (P = 0.57 and P = 1.0 respectively) and no coefficient estimates, so it cannot be encoded as a model; see the vignette Errata.",
      sep = "\n"
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters -- typical values at the REFERENCE covariate
    # vector: PD-L1 TPS 1%-49% (most frequent category), EGFR wild type
    # (most frequent category) and AUC_PEMBRO = 7079 mg*day/L.
    # All values from supplementary Table S6 ("Parameter and uncertainty
    # estimates of the final covariate-containing model of tumor growth").
    # ---------------------------------------------------------------------
    lkgrowth <- log(0.00114)
    label("log first-order tumor growth rate constant kgrowth of the resistant fraction (1/day)")   # Table S6: k growth = 0.00114 1/day, RSE 22.7%

    lkdeath <- log(0.00265)
    label("log first-order tumor kill rate constant kdeath of the sensitive fraction (1/day), at the reference covariate vector")  # Table S6: k death = 0.00265 1/day, RSE 21.0%

    # f is estimated on the LOGIT scale: Supplementary Methods, Handling of
    # Covariates -- "For the logit normally distributed f parameter (fraction
    # of the tumor on which killing is occurring), both continuous and
    # categorical covariates were added in a linear fashion to the typical
    # value of the logit transform of the parameter (TVlogit(f))". Table S6
    # reports f on the natural (0, 1) scale, so it is transformed here.
    # Self-consistency check of this reading: expit(logit(0.574) - 1.81) =
    # 0.181, giving a wild-type / mutant f ratio of 3.18, which matches the
    # supplementary text's "median f was 3.2-fold higher in patients with
    # EGFR wild-type versus mutant tumors".
    logitfresp <- qlogis(0.574)
    label("logit of the treatment-sensitive tumor fraction f (unitless), at the reference covariate vector")  # Table S6: f = 0.574, RSE 14.0% (natural scale)

    # ---------------------------------------------------------------------
    # Covariate effects
    # ---------------------------------------------------------------------
    # Continuous covariate parameterization (Supplementary Methods, Handling
    # of Covariates): P* = theta_x * (COV / median)^theta_y, matching the
    # main-article Methods equation
    #   kdeath = TVkdeath * (AUCss-6weeks / AUCtypical,ss-6weeks)^theta.
    e_auc_kdeath <- 0.196
    label("Power exponent of (AUC_PEMBRO / 7079) on kdeath (unitless)")                             # Table S6: AUC ss-6weeks on k death = 0.196, RSE 71.4%; 95% CI -0.0784 to 0.47, P = 0.54 (not significant; retained by the authors for simulation)

    # Categorical covariate parameterization (Supplementary Methods, Handling
    # of Covariates): P* = theta_x for the most frequent category and
    # P* = theta_x * (1 + theta_y) for every other category, so theta_y is a
    # FRACTIONAL deviation from the reference. Reference = TPS 1%-49%.
    # Self-consistency check: kdeath(TPS >= 50%) / kdeath(TPS < 1%) =
    # (1 + 1.74) / (1 - 0.377) = 4.40, matching the supplementary text's
    # "median k death was 4.7-fold higher in strongly PD-L1-positive versus
    # PD-L1-negative patients" (a post-hoc EBE median ratio, so approximate).
    e_pdl1_ge50_kdeath <- 1.74
    label("Fractional deviation of kdeath for PD-L1 TPS >= 50% vs the TPS 1%-49% reference (unitless)")   # Table S6: PD-L1_1 on k death = 1.74, RSE 33%

    e_pdl1_lt1_kdeath <- -0.377
    label("Fractional deviation of kdeath for PD-L1 TPS < 1% vs the TPS 1%-49% reference (unitless)")     # Table S6: PD-L1_2 on k death = -0.377, RSE -48.8%

    e_pdl1_missing_kdeath <- 0.268
    label("Fractional deviation of kdeath for unknown PD-L1 status vs the TPS 1%-49% reference (unitless)")  # Table S6: PD-L1_3 on k death = 0.268, RSE 145% (poorly identified; see supplementary text, Final Exposure-Response Model)

    # Categorical covariates on the logit-normal f enter ADDITIVELY on the
    # logit scale (Supplementary Methods, Handling of Covariates): for the
    # most frequent category TVlogit(f)* = theta_x, for other categories
    # TVlogit(f)* = theta_x + theta_y. Reference = EGFR wild type.
    e_egfr_mut_fresp <- -1.81
    label("Additive deviation on logit(f) for EGFR-mutant vs the EGFR wild-type reference (unitless)")     # Table S6: EGFR_1 on f = -1.81, RSE -28.3%

    e_egfr_missing_fresp <- 1.66
    label("Additive deviation on logit(f) for unknown EGFR status vs the EGFR wild-type reference (unitless)")  # Table S6: EGFR_2 on f = 1.66, RSE 42.3%

    # ---------------------------------------------------------------------
    # Inter-individual variability -- full 3x3 covariance block on kgrowth,
    # kdeath and logit(f). Supplementary Methods, Structural Model Selection:
    # "visual inspection of scatterplots of individual random effects
    # indicated associations between kdeath, f, and kgrowth. Consequently, a
    # full covariance matrix was implemented."
    #
    # Supplementary Table S6 heads this block "Inter-Individual Variability
    # Estimates (Covariance Matrix)" and lists the entries in NONMEM
    # lower-triangular BLOCK(3) order -- omega11, omega21, omega22, omega31,
    # omega32, omega33 -- so the diagonal entries are VARIANCES and the
    # off-diagonal entries are COVARIANCES:
    #   omega^2 k growth        = 1.21    RSE 22.8%, shrinkage 36.9%
    #   omega_xy k death:k growth = -0.33   RSE -85.2%
    #   omega^2 k death         = 1.26    RSE 22.9%, shrinkage 36.2%
    #   omega_xy f:k growth     = -0.814  RSE -65.7%
    #   omega_xy f:k death      = 0.631   RSE 35.8%
    #   omega^2 f               = 2.79    RSE 20.1%, shrinkage 31.9%
    # Implied correlations -0.267 (kgrowth,kdeath), -0.443 (kgrowth,f) and
    # 0.337 (kdeath,f); the matrix is positive definite (eigenvalues 3.37,
    # 1.04, 0.85). kgrowth and kdeath are log-normal and f is logit-normal
    # per the main-article Methods and Supplementary Methods respectively.
    etalkgrowth + etalkdeath + etalogitfresp ~ c(1.21,
                                                 -0.33,  1.26,
                                                 -0.814, 0.631, 2.79)

    # ---------------------------------------------------------------------
    # Residual error -- Supplementary Methods, Structural Model Selection:
    # "An exponential residual error model was found to adequately
    # characterize residual variability." Table S6 reports the VARIANCE.
    expSd <- sqrt(0.0274)
    label("Exponential residual error SD on tumor size (log scale)")                                 # Table S6: Exponential Residual Error Variance = 0.0274, RSE 4.16%, shrinkage 6.8% -> SD = sqrt(0.0274) = 0.1655
  })

  model({
    # -------------------------------------------------------------------
    # 1. Derived covariate indicators
    # -------------------------------------------------------------------
    # PD-L1 TPS category indicators derived from the continuous TPS column
    # (Supplementary Methods, Handling of Covariates: TPS >= 50% strongly
    # positive, 1%-49% weakly positive, < 1% negative). Both are gated on
    # PD-L1 being available so that the unknown stratum receives only its own
    # coefficient; the reference category 1% <= TPS < 50% fires neither.
    pdl1_known <- 1 - PDL1_TUM_MISSING
    pdl1_ge50  <- pdl1_known * (PDL1_TUM >= 50)
    pdl1_lt1   <- pdl1_known * (PDL1_TUM < 1)

    # EGFR mutant indicator, gated on EGFR status being known for the same
    # reason (EGFR wild type is the reference category).
    egfr_mut <- (1 - TUM_EGFR_MUT_MISSING) * TUM_EGFR_MUT

    # -------------------------------------------------------------------
    # 2. Individual tumor-dynamics parameters
    # -------------------------------------------------------------------
    kgrowth <- exp(lkgrowth + etalkgrowth)

    # kdeath: multiplicative fractional-deviation terms for the three
    # non-reference PD-L1 levels (only one can fire per subject), times the
    # power-form exposure effect normalized to the population-typical
    # AUCss-6weeks of 7079 mg*day/L (main-article Methods).
    kdeath <- exp(lkdeath + etalkdeath) *
      (1 + e_pdl1_ge50_kdeath    * pdl1_ge50) *
      (1 + e_pdl1_lt1_kdeath     * pdl1_lt1) *
      (1 + e_pdl1_missing_kdeath * PDL1_TUM_MISSING) *
      (AUC_PEMBRO / 7079)^e_auc_kdeath

    # f: additive covariate terms and the IIV both act on the logit scale.
    logit_fresp <- logitfresp +
      e_egfr_mut_fresp     * egfr_mut +
      e_egfr_missing_fresp * TUM_EGFR_MUT_MISSING +
      etalogitfresp
    fresp <- expit(logit_fresp)

    # -------------------------------------------------------------------
    # 3. ODE system
    # -------------------------------------------------------------------
    # Main-article Methods, tumor size NLME model structure:
    #   Tumor size = Baseline * [(1 - f) * exp(kgrowth * time) +
    #                            f * exp(-kdeath * max(0, time - delay))]
    # Encoded as two sub-states of the baseline diameter:
    #   growth(t) = (1 - f) * TUM_SLD * exp(kgrowth * t)   -- resistant fraction
    #   shrink(t) = f       * TUM_SLD * exp(-kdeath * max(0, t - delay))
    #   TS        = growth + shrink
    # so TS(0) = (1 - f) * TUM_SLD + f * TUM_SLD = TUM_SLD, the observed
    # baseline. The max(0, t - delay) is realised by gating the decay on
    # t >= T_SCAN_TO_DOSE: before the first dose the sensitive sub-state is
    # held at its initial value; the resistant sub-state grows from t = 0.
    dosed <- t >= T_SCAN_TO_DOSE

    d/dt(growth) <-  kgrowth * growth
    d/dt(shrink) <- -kdeath * shrink * dosed

    growth(0) <- (1 - fresp) * TUM_SLD
    shrink(0) <- fresp * TUM_SLD

    # -------------------------------------------------------------------
    # 4. Observation and error
    # -------------------------------------------------------------------
    TS <- growth + shrink
    TS ~ lnorm(expSd)
  })
}
attr(Chatterjee_2016_pembrolizumab, "message") <-
  "Exposure-response tumor-size model for pembrolizumab in advanced NSCLC (Chatterjee 2016; KEYNOTE-001, n = 496). Observable `TS` is the RECIST 1.1 sum of longest diameters in mm. Bi-exponential: a treatment-sensitive fraction f decaying at kdeath plus a resistant fraction (1 - f) growing at kgrowth, both initialised from the observed baseline covariate TUM_SLD. No PK input -- supply exposure per subject as AUC_PEMBRO (mg*day/L), which the source analysis computed as dose/CL from the companion popPK model packaged here as Ahamadi_2017_pembrolizumab. Required covariates: TUM_SLD, AUC_PEMBRO, PDL1_TUM, PDL1_TUM_MISSING, TUM_EGFR_MUT, TUM_EGFR_MUT_MISSING, T_SCAN_TO_DOSE. The exposure effect on kdeath (0.196) is NOT statistically significant (95% CI -0.0784 to 0.47, P = 0.54); the paper's conclusion is that response is flat across 2-10 mg/kg."
Chatterjee_2016_pembrolizumab
