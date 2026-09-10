Zhang_2025_olanzapine <- function() {
  description <- "One-compartment population PK model for oral olanzapine with first-order absorption in adults with major depressive disorder, built from a routine therapeutic-drug-monitoring database (Zhang 2025). Apparent oral clearance is allometrically scaled on body weight and reduced by 28.9% when paroxetine is co-administered; the absorption rate constant is fixed to a published value. Between-subject variability was retained on CL/F only."
  reference <- paste(
    "Zhang C, Chen L, Duan YY, He SM, Tian YL, Gao Y, Wang DD.",
    "Drug-drug interaction of paroxetine on olanzapine and initial dosage",
    "optimization in patients with major depressive disorder based on population",
    "pharmacokinetics.",
    "Frontiers in Psychiatry. 2025;16:1538996.",
    "doi:10.3389/fpsyt.2025.1538996.",
    "Final model Equations (6) and (7); parameter estimates Table 3.",
    "The fixed absorption rate constant is quoted from the paper's reference 26:",
    "Sun L, Mills R, Sadler BM, Rege B. Population pharmacokinetics of olanzapine and",
    "samidorphan when administered in combination in healthy subjects and patients with",
    "schizophrenia. J Clin Pharmacol. 2021;61(11):1430-1441. doi:10.1002/jcph.1911.",
    "Companion analysis in schizophrenia by the same group:",
    "modellib('Zhang_2024_olanzapine').",
    sep = " "
  )
  vignette <- "Zhang_2025_olanzapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "olanzapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "olanzapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with a 70 kg reference weight, exponent 0.75 on CL/F and 1 on V/F,",
        "both fixed to the Anderson & Holford (2008) canonical values rather than estimated",
        "(Methods, Equation 3, and Results, Equations 6-7). Cohort mean 61.83 kg (SD 11.82),",
        "median 60.00 kg, range 40.00-92.00 kg (Table 1). Treated as a baseline, time-fixed",
        "covariate: the paper draws weight from a retrospective medical record system and",
        "reports a single value per patient."
      ),
      source_name        = "weight"
    ),
    CONMED_PAROXETINE = list(
      description        = "Concomitant paroxetine indicator (1 = patient co-prescribed paroxetine hydrochloride tablets, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant paroxetine)",
      notes              = paste(
        "The paper's `PAR` variable (Results, Equation 6: 'PAR represents paroxetine; when",
        "patients took paroxetine, PAR was 1, otherwise PAR was 0'). 18 of the 72 patients",
        "were on paroxetine hydrochloride tablets (Table 2). This was the only one of the 24",
        "screened concomitant drugs, and the only covariate other than weight, retained by the",
        "stepwise search (Supplementary Table S1: forward inclusion dOFV -8.295, P < 0.05;",
        "backward elimination dOFV +8.295, P < 0.01). The paper attributes the effect to",
        "paroxetine's inhibition of CYP2D6, by which olanzapine is metabolised (Discussion);",
        "it is an empirical exposure shift, not a mechanistically parameterised interaction.",
        "Time-fixed at the analysis baseline in this source."
      ),
      source_name        = "PAR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a potential covariate (Methods, Covariate model) but not retained; cohort mean 48.14 years (SD 20.94), median 52.08, range 16.00-87.90 (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained; 55 women / 17 men (Table 1). Supplementary Table S1 model 26, 'Effect of Gender on CL', dOFV -3.129, P > 0.05. The paper reports 'Gender (men/women)', so a SEXM-oriented source column would need value inversion."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 39.51 g/L (SD 3.41; Table 1)."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 52.72 IU/L (SD 109.29; Table 1)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 37.38 IU/L (SD 48.98; Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 53.88 umol/L (SD 12.18; Table 1)."
    ),
    BUN = list(
      description = "Serum urea",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.48 mmol/L (SD 1.12; Table 1). The paper's Table 1 row is labelled 'Urea (mmol/L)'."
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.66 mmol/L (SD 1.25; Table 1)."
    ),
    TRIG = list(
      description = "Serum triglyceride",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 2.19 mmol/L (SD 1.60; Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 2.43 umol/L (SD 1.55; Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 8.44 umol/L (SD 3.65; Table 1)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 38.16% (SD 3.47; Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 126.81 g/L (SD 13.04; Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72,
    n_studies      = 1,
    age_range      = "16.00-87.90 years",
    age_median     = "52.08 years",
    age_mean       = "48.14 years (SD 20.94)",
    weight_range   = "40.00-92.00 kg",
    weight_median  = "60.00 kg",
    weight_mean    = "61.83 kg (SD 11.82)",
    sex_female_pct = 76.4,
    disease_state  = "Major depressive disorder; inpatients on routine oral olanzapine therapy.",
    dose_range     = "Not reported. Concentrations came from a routine therapeutic-drug-monitoring database and the paper tabulates no administered olanzapine doses or sampling times.",
    regions        = "China (single centre: Xuzhou Oriental Hospital Affiliated to Xuzhou Medical University, Jiangsu)",
    notes          = paste(
      "Retrospective analysis of olanzapine therapeutic-drug-monitoring concentrations",
      "collected between December 2020 and August 2023 (Methods, Data Collection);",
      "1-3 concentration samples per patient (Results, Patient information).",
      "Baseline demographics and laboratory values are Table 1; concomitant medication",
      "counts are Table 2. Estimation used NONMEM 7.",
      "The covariate search screened the demographic and clinical-chemistry variables of",
      "Table 1 -- age, sex, weight, albumin, globulin, alanine transaminase, aspartate",
      "transaminase, creatinine, urea, total protein, total cholesterol, triglyceride,",
      "direct and total bilirubin, hematocrit, hemoglobin, mean corpuscular hemoglobin and",
      "mean corpuscular hemoglobin concentration -- together with the 24 concomitant drugs",
      "of Table 2 (atorvastatin, alprazolam, amlodipine, benzhexol, buspirone, clonazepam,",
      "dexzopiclone, duloxetine, aspirin, escitalopram, irbesartan/hydrochlorothiazide,",
      "levodopa/benserazide, lorazepam, metoprolol, mirtazapine, omeprazole, oxazepam,",
      "paroxetine, propranolol, sertraline, trazodone, valsartan, venlafaxine, zopiclone).",
      "Only body weight and concomitant paroxetine survived (inclusion dOFV > 3.84,",
      "exclusion dOFV > 6.63). The per-drug hypothesis tests are Supplementary Table S1",
      "(base-model OFV 453.159): paroxetine gave the largest single-covariate drop",
      "(dOFV -8.295), and amlodipine was the only other drug significant on its own",
      "(dOFV -5.935) but failed to add to the paroxetine model in the second forward step",
      "(dOFV -2.134, P > 0.05). Globulin, total protein, mean corpuscular hemoglobin and",
      "mean corpuscular hemoglobin concentration have no canonical covariate-column entry",
      "and so are recorded in this note rather than in covariatesDataExcluded.",
      "The model was qualified by goodness-of-fit plots, a VPC, and a bootstrap (Table 3)."
    )
  )

  ini({
    # Structural parameters -- Results, Equations (6) and (7), point estimates
    # also in Table 3. Reference subject: 70 kg, no concomitant paroxetine.
    lcl <- log(19.6); label("Apparent oral clearance CL/F at 70 kg without paroxetine (L/h)")  # Table 3 / Eq. 6: CL/F = 19.6 L/h (SE 7.1%; bootstrap median 19.4 [17.0, 21.9])
    lvc <- log(197);  label("Apparent central volume of distribution V/F at 70 kg (L)")        # Table 3 / Eq. 7: V/F = 197 L (SE 15.3%; bootstrap median 194 [154, 277])

    # Ka was not estimated here. Methods, Modeling: "absorption rate constants
    # [Ka, fixed at 0.861/h (26)]"; Table 3 prints "0.861 (fixed)" with no SE
    # and no bootstrap row. Reference 26 is Sun 2021 (olanzapine + samidorphan).
    lka <- fixed(log(0.861)); label("First-order oral absorption rate constant ka, from Sun 2021 (1/h)")  # Table 3: Ka = 0.861 1/h (fixed), quoted from Sun 2021

    # Allometric exponents. Methods, Equation (3): "F is the allometric
    # coefficient: 0.75 for the CL/F and 1 for the V/F", citing Anderson &
    # Holford (2008). Both are asserted, not estimated -- neither appears in
    # Table 3 and neither carries an SE or a bootstrap interval.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F for WT/70, from Anderson & Holford 2008 (unitless)")  # Methods Eq. 3 and Results Eq. 6 exponent 0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on V/F for WT/70, from Anderson & Holford 2008 (unitless)")   # Methods Eq. 3 and Results Eq. 7 exponent 1

    # Concomitant-paroxetine effect on CL/F, entered in the linear-shift form
    # of Methods Equation (5) for categorical covariates. Results Equation (6)
    # prints the assembled term as (1 - 0.289 * PAR), i.e. a 28.9% reduction in
    # apparent oral clearance; the Results and Discussion state the same result
    # as a clearance ratio of 0.711:1 with vs without paroxetine (Figure 1H),
    # and 1 - 0.289 = 0.711 exactly.
    e_par_cl <- -0.289; label("Proportional change in CL/F with concomitant paroxetine (fraction)")  # Table 3: theta_PAR = -0.289 (SE 30.5%; bootstrap median -0.283 [-0.420, -0.078])

    # Between-subject variability. Methods Equation (1) is Ai = TV(A) * exp(eta_i),
    # i.e. exponential (log-normal) IIV, and defines eta as having "variance
    # omega^2". Table 3's row is labelled `omega_CL/F` (not omega^2), so the
    # printed 0.434 is the STANDARD DEVIATION and the variance is 0.434^2.
    # See the vignette Errata: reading 0.434 as the variance (SD 0.659) drives
    # the model's own target-attainment predictions 20-30 percentage points
    # below every cell of the paper's Table 4, whereas the SD reading
    # reproduces them to within a few points. 0.434 on the log scale is
    # 45.4% CV.
    etalcl ~ 0.434^2  # Table 3: omega_CL/F = 0.434 (SE 11.1%; bootstrap median 0.429 [0.336, 0.535]), read as an SD

    # No IIV on V/F or ka: Table 3 reports a single omega row, for CL/F only.

    # Residual error. Methods Equation (2) is Bi = Ci + Ci*eps_1 + eps_2,
    # a combined proportional-plus-additive model, and defines eps as having
    # "variance sigma^2". Table 3's rows are labelled `sigma_1` and `sigma_2`,
    # so both printed values are STANDARD DEVIATIONS on the same reading used
    # for omega above.
    propSd <- 0.153; label("Proportional residual error (fraction)")  # Table 3: sigma_1 = 0.153 (SE 16.6%; bootstrap median 0.150 [0.061, 0.199])
    addSd  <- 1.005; label("Additive residual error (ng/mL)")         # Table 3: sigma_2 = 1.005 (SE 45.0%; bootstrap median 1.005 [0.306, 2.186])
  })

  model({
    # 1. Individual parameters. Results Equations (6) and (7):
    #      CL/F = 19.6 * (weight/70)^0.75 * (1 - 0.289 * PAR)
    #      V/F  = 197  * (weight/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_par_cl * CONMED_PAROXETINE)
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # 2. Micro-constants
    kel <- cl / vc

    # 3. ODE system -- one-compartment with first-order oral absorption.
    #    F is not separately identifiable from a TDM dataset with oral dosing
    #    only, so cl and vc are the apparent (CL/F, V/F) quantities and no
    #    bioavailability term is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Observation. Doses are in mg and vc is in L, giving mg/L; the factor
    #    of 1000 converts to the ng/mL in which the paper reports olanzapine
    #    concentrations (therapeutic window 20-80 ng/mL, Methods, Simulation).
    Cc <- 1000 * central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
