Chatterjee_2017_pembrolizumab_consolidated <- function() {
  description <- paste(
    "Consolidated exposure-response tumor-size (sum of longest diameters,",
    "SLD) model for pembrolizumab in advanced melanoma, developed by",
    "Chatterjee et al. (Merck) on the pooled KEYNOTE-001, -002 and -006",
    "melanoma data (N = 1,366, April 2015 cutoff). The structural model",
    "splits the OBSERVED baseline SLD into a treatment-accessible fraction f",
    "that decays first-order at kdeath after a per-subject delay and a",
    "resistant fraction (1 - f) that grows first-order at kgrowth, giving",
    "the published form",
    "SLD(t) = baseline * [(1 - f) * exp(kgrowth * t) +",
    "f * exp(-kdeath * max(0, t - delay))].",
    "Pembrolizumab exposure enters as a log-linear (power) effect of the",
    "steady-state 6-week AUC on kdeath, normalized to 7079 mg*day/L, with",
    "SEPARATE exponents for ipilimumab-naive and ipilimumab-experienced",
    "patients because dose assignment was confounded with prior-ipilimumab",
    "status in the pooled dataset. The final covariate model adds PD-L1",
    "status (positive / negative / unknown) and baseline tumor size on",
    "kdeath, BRAF mutation status on kgrowth, and baseline tumor size plus",
    "prior ipilimumab treatment on the logit of f. Neither exposure slope is",
    "statistically significant (P = 0.20 and P = 0.25); the authors retained",
    "them only so the magnitude of a potential exposure-response",
    "relationship could be simulated, and the paper's conclusion is that",
    "response is flat over the fivefold 2-10 mg/kg dose range. There is no",
    "PK input: exposure is supplied per subject as the covariate",
    "AUC_PEMBRO, which the source analysis obtained from the companion",
    "pembrolizumab population-PK model (Ahamadi 2017; available in this",
    "library as Ahamadi_2017_pembrolizumab)."
  )
  reference <- paste(
    "Chatterjee MS, Elassaiss-Schaap J, Lindauer A, Turner DC, Sostelly A,",
    "Freshwater T, Mayawala K, Ahamadi M, Stone JA, de Greef R, Kondic AG,",
    "de Alwis DP.",
    "Population pharmacokinetic/pharmacodynamic modeling of tumor size",
    "dynamics in pembrolizumab-treated advanced melanoma.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(1):29-39.",
    "doi:10.1002/psp4.12140. PMID: 27896901. PMCID: PMC5270297.",
    "Structural equation from the main-article Methods ('Consolidated",
    "exposure-response tumor size model'); all final parameter values from",
    "supplementary Table 4C ('Parameter and uncertainty estimates of the",
    "final covariate-containing tumor model'); covariate functional forms,",
    "reference categories and centering constants from the supplementary",
    "NONMEM control stream ('Final tumor size model, consolidated modeling",
    "approach') and the Supplementary Methods ('Consolidated Model Covariate",
    "Parameterization').",
    sep = " "
  )
  vignette <- "Chatterjee_2017_pembrolizumab_melanoma"

  units <- list(
    time = "day",
    dosing = "n/a (no PK input; pembrolizumab exposure enters as the per-subject covariate AUC_PEMBRO in mg*day/L)",
    concentration = "mm (the observable `TS` is the RECIST 1.1 sum of the longest diameters of target lesions)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Both states are sub-fractions of the measured tumor
  # diameter (mm), not amounts, so `units` is NA. verified = TRUE: the
  # accessible / resistant split is stated explicitly in the main-article
  # Methods and drawn in Figure 1a ("the labels 'f' and '1-f' represent
  # proportions of target tumor tissue that are accessible to treatment and
  # undergoing unimpeded exponential growth, respectively").
  compartmentData <- list(
    growth = list(analyte = "tumor-size", units = NA_character_, specimen = "tumor", verified = TRUE),
    shrink = list(analyte = "tumor-size", units = NA_character_, specimen = "tumor", verified = TRUE)
  )

  covariateData <- list(
    TUM_SLD = list(
      description = "Observed baseline sum of the longest diameters of target lesions per RECIST 1.1, measured at screening.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      source_name = "BASE",
      notes = paste(
        "The paper's 'baseline' term: 'baseline is the tumor size at screening' (main-article Methods). It is a per-subject regressor, NOT an estimated parameter, and it initialises both sub-states: growth(0) = (1 - f) * TUM_SLD and shrink(0) = f * TUM_SLD, so TS(0) = TUM_SLD exactly.",
        "It also carries TWO retained covariate effects, both centered at the source median 80.75 mm (supplementary NONMEM control stream):",
        "  (a) power effect on kdeath, (TUM_SLD / 80.75)^-0.186 (KDYINGBASE block; supplementary Table 4C 'Baseline tumor size on kdeath');",
        "  (b) linear effect on the logit of f, -0.00541 * (TUM_SLD - 80.75) (PHIBASE block; supplementary Table 4C 'Baseline Tumor Size on f').",
        "Both coefficients are negative: a larger baseline tumor both shrinks more slowly and has a smaller accessible fraction. Baseline tumor size on f was the FIRST relationship added in the forward search (dOFV 32, P = 1.42e-08; supplementary Table 4B).",
        "The centering constant 80.75 mm is the source control stream's value; the cohort median printed in supplementary Table 5B is 80.60 mm (range 10.00-895.00, N = 1366, none missing). Use 80.75, which is what the published coefficients were estimated against.",
        sep = " "
      )
    ),
    AUC_PEMBRO = list(
      description = "Per-subject pembrolizumab area under the serum concentration-time curve at steady state over a 6-week interval (AUCss-6weeks).",
      units = "mg*day/L (equivalently ug*day/mL; numerically identical)",
      type = "continuous",
      reference_category = NULL,
      source_name = "AUC2",
      notes = paste(
        "Log-linear (power) effect on kdeath with a SEPARATE exponent per prior-ipilimumab stratum (supplementary NONMEM control stream AUCCOV block): kdeath = TVkdeath * (AUC_PEMBRO / 7079)^theta, theta = 0.131 for IPI-naive and 0.100 for IPI-experienced.",
        "The normalizing constant 7079 mg*day/L is annotated in the source control stream as 'the median exposure at the 10 mpk Q3 regime, in the melanoma submission'. It is the same constant the companion NSCLC analysis uses (Chatterjee_2016_pembrolizumab).",
        "The stratified parameterization exists because dose was confounded with prior-ipilimumab status: IPI-naive patients (who respond better) were over-represented at 10 mg/kg and under-represented at 2 mg/kg (main-article Methods and Table 1; supplementary Figure 2).",
        "Neither slope is statistically significant (P = 0.20 IPI-naive, P = 0.25 IPI-experienced; main-article Results) and the RSEs are 77.1% and 87.2%. The authors retained both for simulation only. Do not interpret either as an established exposure-response relationship.",
        "The source analysis did not model pembrolizumab PK here; AUC estimates came from a separate population-PK analysis (Ahamadi 2017, packaged as Ahamadi_2017_pembrolizumab). Records with AUC2 < 0.01 were dropped from the estimation dataset ($DATA IGNORE).",
        sep = " "
      )
    ),
    PRIOR_IPI = list(
      description = "Prior ipilimumab treatment indicator: 1 = previously treated with (or refractory to) ipilimumab, 0 = ipilimumab-naive.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ipilimumab-naive, the most frequent category: 746 of 1366, 54.61%)",
      source_name = "IPIN",
      notes = paste(
        "Two roles in the model. (a) Additive deviation on the logit of f: -0.964 for IPI-experienced versus the IPI-naive reference (supplementary NONMEM control stream PHIIPIN block; supplementary Table 4C 'IPIN_1 on f'). (b) Selector between the two AUC exponents on kdeath (AUCCOV block).",
        "Effect size on f: expit(logit(0.696) - 0.964) = 0.466, so the typical accessible tumor fraction falls from 69.6% in an IPI-naive patient to 46.6% in an IPI-experienced patient at the median baseline tumor size.",
        "Added in the third forward step of the covariate search (dOFV 22, P = 0.000003; supplementary Table 4B).",
        "Cohort split (supplementary Table 5A): IPI-naive 746 (54.61%), prior IPI 620 (45.39%); no unknowns.",
        sep = " "
      )
    ),
    PDL1_TUM_POS = list(
      description = "Baseline tumor PD-L1 positivity indicator by the source study's immunohistochemistry assay and threshold: 1 = PD-L1-positive, 0 = PD-L1-negative.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (PD-L1-positive, the most frequent category: 739 of 1366, 54.1%)",
      source_name = "PDL1 (coded 1 = positive, 0 = negative, 6 = unknown)",
      notes = paste(
        "Fractional-deviation effect on kdeath: kdeath * (1 - 0.614) for PD-L1-negative relative to the PD-L1-positive reference (supplementary NONMEM control stream KDYINGPDL1 block; supplementary Table 4C 'PD-L1_2 on kdeath').",
        "Note the unusual reference: the source's most-frequent-category rule makes PD-L1-POSITIVE the reference, so both PD-L1 coefficients are negative deviations from it.",
        "Set PDL1_TUM_POS to 0 when PDL1_TUM_MISSING = 1; model() gates the negative-stratum indicator on PD-L1 being known, so the unknown stratum receives only its own coefficient.",
        "The source paper dichotomises PD-L1 to positive / negative ('PD-L1 expression level grouped into two categories', supplementary Table 4B footnote a) and does NOT print the assay clone or the percent-staining threshold, so the continuous canonical PDL1_TUM cannot be derived from it without inventing a cut-point. Carry this binary column instead and record the study's assay in covariateData notes when it is known.",
        "PD-L1 on kdeath was the second relationship added in the forward search and the largest single dOFV (37, P = 8.66e-09; supplementary Table 4B). The main-article Discussion: 'Identification of PD-L1 expression as a key determinant of kdeath is consistent with the known mechanism of action for pembrolizumab'.",
        sep = " "
      )
    ),
    PDL1_TUM_MISSING = list(
      description = "Binary indicator: 1 = baseline tumor PD-L1 status could not be assigned for this subject, 0 = a PD-L1 result is available.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (PD-L1 result available)",
      source_name = "PDL1 = 6 ('Unknown')",
      notes = paste(
        "The source analysis retained 'unknown' PD-L1 as a third modelled category with its own coefficient rather than imputing it into a measured level: kdeath * (1 - 0.21) relative to the PD-L1-positive reference (supplementary Table 4C 'PD-L1_1 on kdeath').",
        "402 of 1366 subjects (29.43%) are in this stratum (supplementary Table 5A) -- a much larger share than in the companion NSCLC analysis. As with any missingness coefficient, the estimate mixes PD-L1-positive and PD-L1-negative tumors and should not be extrapolated to a measured subgroup.",
        "Mutually exclusive with the derived PD-L1-negative indicator.",
        sep = " "
      )
    ),
    TUM_BRAF_MUT = list(
      description = "Tumor BRAF (B-type Raf proto-oncogene) mutation indicator: 1 = BRAF-mutant (or translocated) tumor, 0 = BRAF wild-type tumor.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (BRAF wild type, the most frequent category: 985 of 1366, 72.11%)",
      source_name = "BRAFMN (coded 1 = mutation/translocation, 2 = wild type, -99 = unknown)",
      notes = paste(
        "Fractional-deviation effect on kgrowth: kgrowth * (1 + 0.702) for BRAF-mutant relative to the wild-type reference (supplementary NONMEM control stream KGROWBRAFMN block; supplementary Table 4C 'BRAF_1 on kgrowth'). BRAF-mutant melanoma therefore grows about 70% faster in the resistant fraction.",
        "The source stream maps unknown BRAF status (-99) to the wild-type reference multiplier of 1, so set TUM_BRAF_MUT = 0 for an unknown subject; no separate missingness column is needed and only 11 of 1366 subjects (0.81%) are affected (supplementary Table 5A).",
        "Added in the fourth forward step (dOFV 13, P = 0.000309; supplementary Table 4B).",
        sep = " "
      )
    ),
    T_SCAN_TO_DOSE = list(
      description = "Per-subject delay between the baseline tumor scan (the model's time origin) and the onset of drug-driven tumor shrinkage.",
      units = "day",
      type = "continuous",
      reference_category = NULL,
      source_name = "DLYBFRDOSE",
      notes = paste(
        "The paper's 'delay' term in the structural equation, described in the main-article Methods as 'the lag in onset of drug activity for tumor shrinkage, interpreted as the time required for immune system activation'. The source column name DLYBFRDOSE ('delay before dose') and the companion NSCLC analysis instead describe it as the gap between the baseline scan and the first dose; both readings place it in the same place in the equation.",
        "Enters only through max(0, time - T_SCAN_TO_DOSE) (supplementary NONMEM control stream: DYINGTIME = TIME - FLARE - DLYBFRDOSE, floored at 0, with FLARE fixed to 0). Before that time the accessible sub-state is held at its initial value; the resistant sub-state grows from time 0 regardless.",
        "It is per-subject DATA, not an estimated parameter -- the source carries it as an input column and estimates no delay parameter. The paper does not report its distribution. Set T_SCAN_TO_DOSE = 0 for a subject whose shrinkage starts at the baseline scan; the vignette uses 0 and states the assumption.",
        sep = " "
      )
    )
  )

  # Covariates that the source paper's stepwise covariate search screened but
  # did NOT retain in the final model. Documented so the provenance of the
  # covariate screen survives without creating declared-but-unreferenced
  # covariates. Source: main-article Methods, "Covariate exploration" ("the
  # following covariates were tested for inclusion in a stepwise fashion:
  # PD-L1 expression, ECOG status, demographics (age, gender, and weight),
  # baseline tumor size, IPI pretreatment status, and BRAF mutation status")
  # and supplementary Table 4B, which records that no relationship was
  # dropped in backward elimination.
  covariatesDataExcluded <- list(
    WHO_PS = list(
      description = "Baseline ECOG (Eastern Cooperative Oncology Group) performance status, 0 or 1.",
      units = "(integer score)",
      type = "categorical",
      notes = "Tested and not selected (main-article Methods). Cohort: ECOG 0 in 877 (64.2%), ECOG 1 in 486 (35.58%), unknown in 3 (0.22%) (supplementary Table 5A). It WAS retained on the shallow baseline in the paper's earlier mixture model; see Chatterjee_2017_pembrolizumab_mixture."
    ),
    AGE = list(
      description = "Age at baseline.",
      units = "year",
      type = "continuous",
      notes = "Tested as part of the demographics block and not selected (main-article Methods). Cohort median 62.00 years, range 15.00-94.00, none missing (supplementary Table 5B)."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Tested as part of the demographics block and not selected (main-article Methods, listed as 'gender'). Cohort: 833 male (60.98%), 533 female (39.02%) (supplementary Table 5A)."
    ),
    WT = list(
      description = "Body weight at baseline.",
      units = "kg",
      type = "continuous",
      notes = "Tested as part of the demographics block and not selected (main-article Methods). Cohort median 79.30 kg, range 36.20-209.50, none missing (supplementary Table 5B). Weight still enters indirectly through the mg/kg dose that determines AUC_PEMBRO."
    ),
    RACE = list(
      description = "Race.",
      units = "(categorical)",
      type = "categorical",
      notes = "Summarised in supplementary Table 5A but not among the covariates the main-article Methods lists as tested; the cohort is 1330 of 1366 white (97.36%), which leaves no testable contrast. Recorded here so the covariate screen's provenance is complete."
    )
  )

  population <- list(
    species = "human (adults with unresectable or metastatic melanoma)",
    n_subjects = 1366L,
    n_studies = 3L,
    age_range = "15.00-94.00 years (supplementary Table 5B)",
    age_median = "62.00 years",
    weight_range = "36.20-209.50 kg",
    weight_median = "79.30 kg",
    sex_female_pct = 39.02,
    race_ethnicity = c(
      White = 97.36,
      Asian = 1.39,
      Black = 0.51,
      Multiracial = 0.44,
      AmericanIndianOrAlaskaNative = 0.07,
      NativeHawaiianOrOtherPacificIslander = 0.07,
      Unknown = 0.15
    ),
    disease_state = "advanced (unresectable stage III or stage IV) melanoma; 620 of 1366 (45.39%) previously treated with or refractory to ipilimumab",
    dose_range = "pembrolizumab 2 mg/kg IV Q3W, 10 mg/kg IV Q3W, or 10 mg/kg IV Q2W (not a model input; enters only through AUC_PEMBRO)",
    regions = "pooled KEYNOTE-001 (NCT01295827, phase Ib), KEYNOTE-002 (NCT01704287, phase II) and KEYNOTE-006 (NCT01866319, phase III), multinational",
    notes = paste(
      "Pooled modeling dataset (main-article Methods and Table 1): all available melanoma tumor-size data from patients treated with pembrolizumab as of the April 2015 cutoff, restricted to patients with pharmacokinetic data.",
      "Dose imbalance across the prior-ipilimumab strata is the reason two exposure-response slopes were estimated (main-article Table 1): only 64 of the IPI-naive patients (8.58%) were in the 2 mg/kg Q3W group versus 237 of the IPI-experienced patients (38.23%).",
      "Baseline PD-L1 status (supplementary Table 5A): positive 739 (54.1%), negative 225 (16.47%), unknown 402 (29.43%).",
      "Baseline BRAF status (supplementary Table 5A): mutation/translocation 370 (27.09%), wild type 985 (72.11%), unknown 11 (0.81%).",
      "Baseline sum of longest diameters (supplementary Table 5B): median 80.60 mm, range 10.00-895.00 mm, none missing. The model centers its baseline-tumor-size covariate terms at 80.75 mm.",
      "10.8% of patients in the dataset had no post-baseline tumor measurement and therefore contributed no data beyond baseline; unlike the paper's mixture model this structure does not retain them (main-article Discussion).",
      "Tumor size recorded as the sum of the longest dimensions of target lesions by computed tomography or magnetic resonance imaging, per RECIST version 1.1.",
      "Estimation: NONMEM 7.2.0, SAEM (NBURN 3000, NITER 500) followed by importance sampling for the objective function; covariate search by PsN 3.5.3 stepwise covariate modeling.",
      "The paper reports a BASE consolidated model in supplementary Table 4A; per library policy only the FINAL covariate-containing model (supplementary Table 4C) is packaged here.",
      "The visual predictive check in Figure 2e,f censors simulated patients whose tumors grew more than threefold above baseline after at least 6 months on treatment, to mimic study withdrawal. That censoring rule is a post-processing step for the VPC, not part of the model, and is not encoded here.",
      sep = "\n"
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters -- typical values at the REFERENCE covariate
    # vector: PD-L1-positive (most frequent category), BRAF wild type (most
    # frequent category), IPI-naive (most frequent category), baseline tumor
    # size 80.75 mm (the source centering value) and AUC_PEMBRO = 7079
    # mg*day/L. All values from supplementary Table 4C.
    # ---------------------------------------------------------------------
    lkgrowth <- log(0.000512)
    label("log first-order tumor growth rate constant kgrowth of the resistant fraction (1/day)") # Supplementary Table 4C: kgrowth = 0.000512 1/day, |RSE| 12.3%

    lkdeath <- log(0.00609)
    label("log first-order tumor kill rate constant kdeath of the accessible fraction (1/day), at the reference covariate vector") # Supplementary Table 4C: kdeath = 0.00609 1/day, |RSE| 9.15%

    # f is estimated on the LOGIT scale: Supplementary Methods, Consolidated
    # Model Covariate Parameterization -- "For the logit normally distributed
    # f parameter (fraction of the tumor on which killing is occurring), both
    # continuous and categorical covariates were added in a linear fashion to
    # the typical value of the logit transform of the parameter". The source
    # control stream writes exactly TVPHI = LOG(TVFDYING/(1-TVFDYING)).
    # Supplementary Table 4C reports f on the natural (0, 1) scale, so it is
    # transformed here.
    logitfresp <- qlogis(0.696)
    label("logit of the treatment-accessible tumor fraction f (unitless), at the reference covariate vector") # Supplementary Table 4C: f = 0.696, |RSE| 4.51% (natural scale)

    # ---------------------------------------------------------------------
    # Exposure effects on kdeath -- one power exponent per prior-ipilimumab
    # stratum (source control stream AUCCOV block):
    #   kdeath = TVkdeath * (AUC_PEMBRO / 7079)^theta
    # NEITHER is statistically significant; see covariateData[[AUC_PEMBRO]].
    # ---------------------------------------------------------------------
    e_auc_kdeath_ipinaive <- 0.131
    label("Power exponent of (AUC_PEMBRO / 7079) on kdeath in ipilimumab-naive patients (unitless)") # Supplementary Table 4C: 'IPI naive AUC exponent' = 0.131, |RSE| 77.1%; main-article Results P = 0.20 (not significant)

    e_auc_kdeath_ipiexp <- 0.1
    label("Power exponent of (AUC_PEMBRO / 7079) on kdeath in ipilimumab-experienced patients (unitless)") # Supplementary Table 4C: 'IPI experienced exp-rep AUC exponent' = 0.1, |RSE| 87.2%; main-article Results P = 0.25 (not significant)

    # ---------------------------------------------------------------------
    # Categorical covariate effects. Supplementary Methods, Consolidated
    # Model Covariate Parameterization: P* = theta_x for the most frequent
    # category and P* = theta_x * (1 + theta_y) for every other category, so
    # each theta_y below is a FRACTIONAL deviation from the reference.
    # ---------------------------------------------------------------------
    e_pdl1_missing_kdeath <- -0.21
    label("Fractional deviation of kdeath for unknown PD-L1 status versus the PD-L1-positive reference (unitless)") # Supplementary Table 4C: 'PD-L1_1 on kdeath' = -0.21, |RSE| 43.1% (footnote: deviation of PD-L1-unknown from PD-L1-positive)

    e_pdl1_neg_kdeath <- -0.614
    label("Fractional deviation of kdeath for PD-L1-negative versus the PD-L1-positive reference (unitless)") # Supplementary Table 4C: 'PD-L1_2 on kdeath' = -0.614, |RSE| 10.5% (footnote: deviation of PD-L1-negative from PD-L1-positive)

    e_braf_mut_kgrowth <- 0.702
    label("Fractional deviation of kgrowth for BRAF-mutant versus the BRAF wild-type reference (unitless)") # Supplementary Table 4C: 'BRAF_1 on kgrowth' = 0.702, |RSE| 37.9%

    # ---------------------------------------------------------------------
    # Continuous covariate effects. Power form on kdeath (log-log slope,
    # centered at the source median 80.75 mm) and linear-on-logit form on f,
    # per the Supplementary Methods and the source control stream.
    # ---------------------------------------------------------------------
    e_tumsld_kdeath <- -0.186
    label("Log-log slope of kdeath versus baseline tumor size, centered at 80.75 mm (unitless)") # Supplementary Table 4C: 'Baseline tumor size on kdeath' = -0.186, |RSE| 33.5%

    e_tumsld_fresp <- -0.00541
    label("Additive slope on logit(f) per mm of baseline tumor size above 80.75 mm (1/mm)") # Supplementary Table 4C: 'Baseline Tumor Size on f' = -0.00541, |RSE| 22.6%

    e_ipi_fresp <- -0.964
    label("Additive deviation on logit(f) for ipilimumab-experienced versus the ipilimumab-naive reference (unitless)") # Supplementary Table 4C: 'IPIN_1 on f' = -0.964, |RSE| 21.8%

    # ---------------------------------------------------------------------
    # Inter-individual variability -- full 3x3 covariance block on
    # log(kgrowth), log(kdeath) and logit(f). Supplementary Table 4C heads
    # this block "Interindividual variability matrix" and the Table 4A
    # footnote states the entries "represent the variance of the random
    # effects on log(kgrowth), log(kdeath), and logit(f) and their respective
    # covariances (kdeath:kgrowth, f:kgrowth, and f:kdeath)". They are listed
    # in NONMEM lower-triangular BLOCK(3) order, matching the source control
    # stream's $OMEGA BLOCK(3) layout:
    #   omega^2 kgrowth        = 1.64    |RSE| 11.2%, shrinkage 23.5%
    #   omega_xy kdeath:kgrowth = -1.24   |RSE| 10.2%
    #   omega^2 kdeath         = 1.39    |RSE| 10.1%, shrinkage 24.3%
    #   omega_xy f:kgrowth     = -1.54   |RSE| 14.8%
    #   omega_xy f:kdeath      = 1.07    |RSE| 15.3%
    #   omega^2 f              = 5.53    |RSE| 9.76%, shrinkage 24.6%
    # Implied correlations -0.821 (kgrowth,kdeath), -0.511 (kgrowth,f) and
    # 0.386 (kdeath,f); the matrix is positive definite (eigenvalues 6.47,
    # 1.83, 0.254).
    # ---------------------------------------------------------------------
    etalkgrowth + etalkdeath + etalogitfresp ~ c(
      1.64,
      -1.24, 1.39,
      -1.54, 1.07, 5.53
    )

    # ---------------------------------------------------------------------
    # Residual error -- the source control stream writes
    # Y = EXP(LOG(IPRED) + EPS(1)) with $SIGMA 0.0389541, i.e. an exponential
    # (log-normal) residual whose reported value is a VARIANCE.
    # ---------------------------------------------------------------------
    expSd <- sqrt(0.0389)
    label("Exponential residual error SD on tumor size (log scale)") # Supplementary Table 4C: 'Expon. residual' = 0.0389 (variance), RSE 2.3%, shrinkage 3.04% -> SD = sqrt(0.0389) = 0.1972
  })

  model({
    # -------------------------------------------------------------------
    # 1. Derived covariate indicators
    # -------------------------------------------------------------------
    # PD-L1-negative indicator, gated on PD-L1 being known so the unknown
    # stratum receives only its own coefficient. PD-L1-POSITIVE is the
    # reference category and fires neither term.
    pdl1_neg <- (1 - PDL1_TUM_MISSING) * (1 - PDL1_TUM_POS)

    # -------------------------------------------------------------------
    # 2. Individual tumor-dynamics parameters
    # -------------------------------------------------------------------
    kgrowth <- exp(lkgrowth + etalkgrowth) *
      (1 + e_braf_mut_kgrowth * TUM_BRAF_MUT)

    # Exposure exponent selected by prior-ipilimumab status (source control
    # stream AUCCOV block).
    e_auc_kdeath <- e_auc_kdeath_ipinaive * (1 - PRIOR_IPI) +
      e_auc_kdeath_ipiexp * PRIOR_IPI

    kdeath <- exp(lkdeath + etalkdeath) *
      (1 + e_pdl1_missing_kdeath * PDL1_TUM_MISSING) *
      (1 + e_pdl1_neg_kdeath * pdl1_neg) *
      (TUM_SLD / 80.75)^e_tumsld_kdeath *
      (AUC_PEMBRO / 7079)^e_auc_kdeath

    # f: the continuous and categorical covariate terms and the IIV all act
    # on the logit scale.
    logit_fresp <- logitfresp +
      e_tumsld_fresp * (TUM_SLD - 80.75) +
      e_ipi_fresp * PRIOR_IPI +
      etalogitfresp
    fresp <- expit(logit_fresp)

    # -------------------------------------------------------------------
    # 3. ODE system
    # -------------------------------------------------------------------
    # Main-article Methods:
    #   Tumor size = baseline * [(1 - f) * exp(kgrowth * time) +
    #                            f * exp(-kdeath * max(0, time - delay))]
    # Encoded as two sub-states of the observed baseline diameter:
    #   growth(t) = (1 - f) * TUM_SLD * exp(kgrowth * t)   -- resistant
    #   shrink(t) = f       * TUM_SLD * exp(-kdeath * max(0, t - delay))
    #   TS        = growth + shrink
    # so TS(0) = TUM_SLD exactly. The max(0, .) is realised by gating the
    # decay on time >= T_SCAN_TO_DOSE. The source wraps the sum in ABS();
    # both terms are products of positive quantities, so that is a no-op.
    dosed <- time >= T_SCAN_TO_DOSE

    d / dt(growth) <- kgrowth * growth
    d / dt(shrink) <- -kdeath * shrink * dosed

    growth(0) <- (1 - fresp) * TUM_SLD
    shrink(0) <- fresp * TUM_SLD

    # -------------------------------------------------------------------
    # 4. Observation and error
    # -------------------------------------------------------------------
    TS <- growth + shrink
    TS ~ lnorm(expSd)
  })
}
attr(Chatterjee_2017_pembrolizumab_consolidated, "message") <-
  "Consolidated exposure-response tumor-size model for pembrolizumab in advanced melanoma (Chatterjee 2017; pooled KEYNOTE-001/-002/-006, N = 1,366). Observable `TS` is the RECIST 1.1 sum of longest diameters in mm. Bi-exponential: a treatment-accessible fraction f decaying at kdeath after a per-subject delay plus a resistant fraction (1 - f) growing at kgrowth, both initialised from the observed baseline covariate TUM_SLD. No PK input -- supply exposure per subject as AUC_PEMBRO (mg*day/L), which the source computed from the companion popPK model packaged here as Ahamadi_2017_pembrolizumab. Required covariates: TUM_SLD, AUC_PEMBRO, PRIOR_IPI, PDL1_TUM_POS, PDL1_TUM_MISSING, TUM_BRAF_MUT, T_SCAN_TO_DOSE. Note the reference categories: PD-L1-POSITIVE, BRAF wild type and IPI-naive. Neither exposure exponent on kdeath (0.131 IPI-naive, 0.100 IPI-experienced) is statistically significant (P = 0.20 and P = 0.25); the paper's conclusion is that response is flat across the fivefold 2-10 mg/kg dose range."
Chatterjee_2017_pembrolizumab_consolidated
