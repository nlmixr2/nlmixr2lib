Zhang_2024_olanzapine <- function() {
  description <- "One-compartment population PK model for oral olanzapine with first-order absorption in adults with schizophrenia, built from a routine therapeutic-drug-monitoring database (Zhang 2024). Apparent oral clearance is allometrically scaled on body weight and reduced by 39.2% when aripiprazole is co-administered; the absorption rate constant is fixed to a published value. Between-subject variability was retained on CL/F only."
  reference <- paste(
    "Zhang C, Jiang L, Hu K, Chen L, Zhang YJ, Shi HZ, He SM, Chen X, Wang DD.",
    "Effects of Aripiprazole on Olanzapine Population Pharmacokinetics and Initial",
    "Dosage Optimization in Schizophrenia Patients.",
    "Neuropsychiatric Disease and Treatment. 2024;20:479-490.",
    "doi:10.2147/NDT.S455183.",
    "Final model Equations (6) and (7); parameter estimates Table 3.",
    "The fixed absorption rate constant is quoted from the paper's reference 29:",
    "Sun L, Mills R, Sadler BM, Rege B. Population pharmacokinetics of olanzapine and",
    "samidorphan when administered in combination in healthy subjects and patients with",
    "schizophrenia. J Clin Pharmacol. 2021;61(11):1430-1441. doi:10.1002/jcph.1911.",
    sep = " "
  )
  vignette <- "Zhang_2024_olanzapine"
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
        "(Methods, Equation 3, and Results, Equations 6-7). Cohort mean 63.16 kg (SD 10.57;",
        "Table 1). Treated as a baseline, time-fixed covariate: the paper draws weight from a",
        "retrospective medical log and reports a single value per patient."
      ),
      source_name        = "weight"
    ),
    CONMED_ARIPIPRAZOLE = list(
      description        = "Concomitant aripiprazole indicator (1 = patient co-prescribed aripiprazole orally disintegrating tablets, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant aripiprazole)",
      notes              = paste(
        "The paper's `ARI` variable (Results, Equation 6: 'ARI was aripiprazole, when patients",
        "took aripiprazole, ARI was 1, otherwise ARI was 0'). 4 of the 65 patients were on",
        "aripiprazole orally disintegrating tablets (Table 2). This was the only one of the 30",
        "screened concomitant drugs retained by the stepwise search. The paper attributes the",
        "effect to competition for CYP2D6, of which both drugs are substrates (Discussion);",
        "it is an empirical exposure shift, not a mechanistically parameterised interaction.",
        "Time-fixed at the analysis baseline in this source."
      ),
      source_name        = "ARI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a potential covariate (Methods, Covariate model) but not retained; cohort mean 45.92 years (SD 15.79; Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained; 29 women / 36 men (Table 1). The paper reports 'Gender (men/women)', so a SEXM-oriented source column would need value inversion."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 41.53 g/L (SD 3.65; Table 1)."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 23.52 IU/L (SD 22.73; Table 1)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 23.10 IU/L (SD 13.20; Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 66.16 umol/L (SD 13.79; Table 1)."
    ),
    BUN = list(
      description = "Serum urea",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.54 mmol/L (SD 1.39; Table 1). The paper's Table 1 row is labelled 'Urea (mmol/L)'."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 10.56 umol/L (SD 4.59; Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.14 umol/L (SD 1.94; Table 1)."
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.15 mmol/L (SD 0.99; Table 1)."
    ),
    TRIG = list(
      description = "Serum triglyceride",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 1.54 mmol/L (SD 1.08; Table 1)."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 39.30% (SD 4.67; Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 128.66 g/L (SD 16.53; Table 1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 65,
    n_studies      = 1,
    age_mean       = "45.92 years (SD 15.79)",
    weight_mean    = "63.16 kg (SD 10.57)",
    sex_female_pct = 44.6,
    disease_state  = "Schizophrenia; inpatients on routine oral olanzapine therapy.",
    dose_range     = "Not reported. Concentrations came from a routine therapeutic-drug-monitoring database and the paper tabulates no administered olanzapine doses or sampling times.",
    regions        = "China (single centre: Xuzhou Oriental Hospital Affiliated to Xuzhou Medical University, Jiangsu)",
    notes          = paste(
      "Retrospective analysis of olanzapine therapeutic-drug-monitoring concentrations",
      "collected between July 2020 and October 2022 (Methods, Data Collection).",
      "Baseline demographics and laboratory values are Table 1; concomitant medication",
      "counts are Table 2. Estimation used NONMEM 7 with FOCE-I.",
      "The covariate search screened the demographic and clinical-chemistry variables of",
      "Table 1 -- age, sex, weight, albumin, globulin, alanine transaminase, aspartate",
      "transaminase, creatinine, urea, total protein, total cholesterol, triglyceride,",
      "direct and total bilirubin, hematocrit, hemoglobin, mean corpuscular hemoglobin and",
      "mean corpuscular hemoglobin concentration -- together with the 30 concomitant drugs",
      "of Table 2 (alprazolam, amisulpride, aripiprazole, benzhexol, bezafibrate, buspirone,",
      "clonazepam, clopidogrel, docusate sodium, aspirin, escitalopram, finasteride,",
      "irbesartan/hydrochlorothiazide, lamotrigine, lorazepam, metformin, metoprolol,",
      "nimodipine, omeprazole, paliperidone, perospirone, perphenazine, propranolol,",
      "pyrazinamide, risperidone, sertraline, valsartan, zaleplon, ziprasidone, zopiclone).",
      "Only body weight and concomitant aripiprazole survived (inclusion dOFV > 6.63,",
      "exclusion dOFV > 10.8). The per-drug hypothesis tests are Supplementary Table S1,",
      "which is not on disk here. Globulin, total protein, mean corpuscular hemoglobin and",
      "mean corpuscular hemoglobin concentration have no canonical covariate-column entry",
      "and so are recorded in this note rather than in covariatesDataExcluded.",
      "The model was qualified by goodness-of-fit plots, a VPC, and a 1000-replicate",
      "bootstrap (Table 3)."
    )
  )

  ini({
    # Structural parameters -- Results, Equations (6) and (7), point estimates
    # also in Table 3. Reference subject: 70 kg, no concomitant aripiprazole.
    lcl <- log(27.6);  label("Apparent oral clearance CL/F at 70 kg without aripiprazole (L/h)")  # Table 3 / Eq. 6: CL/F = 27.6 L/h (SE 5.6%; bootstrap median 27.3 [24.7, 30.4])
    lvc <- log(854);   label("Apparent central volume of distribution V/F at 70 kg (L)")          # Table 3 / Eq. 7: V/F = 854 L (SE 25.5%; bootstrap median 858 [600, 2530])

    # Ka was not estimated here. Methods, Modeling: "absorption rate constant
    # (Ka, fixed at 0.861/h[29])"; Table 3 prints "0.861 (fixed)" with no SE
    # and no bootstrap row. Reference 29 is Sun 2021 (olanzapine + samidorphan).
    lka <- fixed(log(0.861)); label("First-order oral absorption rate constant ka, from Sun 2021 (1/h)")  # Table 3: Ka = 0.861 1/h (fixed), quoted from Sun 2021

    # Allometric exponents. Methods, Equation (3): "Z represented the allometric
    # coefficient: 0.75 for the CL/F and 1 for the V/F", citing Anderson &
    # Holford (2008). Both are asserted, not estimated -- neither appears in
    # Table 3 and neither carries an SE or a bootstrap interval.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F for WT/70, from Anderson & Holford 2008 (unitless)")  # Methods Eq. 3 and Results Eq. 6 exponent 0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on V/F for WT/70, from Anderson & Holford 2008 (unitless)")   # Methods Eq. 3 and Results Eq. 7 exponent 1

    # Concomitant-aripiprazole effect on CL/F, entered in the linear-shift form
    # of Methods Equation (5) for categorical covariates. Results Equation (6)
    # prints the assembled term as (1 - 0.392 * ARI), i.e. a 39.2% reduction in
    # apparent oral clearance; the Discussion states the same result as a
    # clearance ratio of 0.608:1 with vs without aripiprazole (Figure 1H), and
    # 1 - 0.392 = 0.608 exactly.
    e_ari_cl <- -0.392; label("Proportional change in CL/F with concomitant aripiprazole (fraction)")  # Table 3: theta_ARI = -0.392 (SE 27.8%; bootstrap median -0.377 [-0.535, -0.194])

    # Between-subject variability. Methods Equation (1) is L_i = TV(L) * exp(eta_i),
    # i.e. exponential (log-normal) IIV, and defines eta as having "variance
    # omega^2". Table 3's row is labelled `omega_CL/F` (not omega^2), so the
    # printed 0.316 is the STANDARD DEVIATION and the variance is 0.316^2.
    # See the vignette Errata: reading 0.316 as the variance (SD 0.562) drives
    # the model's own target-attainment predictions 20-30 percentage points
    # below every cell of the paper's Table 4, whereas the SD reading
    # reproduces the once-daily cells to within a few points.
    # 0.316 on the log scale is 32.4% CV.
    etalcl ~ 0.316^2  # Table 3: omega_CL/F = 0.316 (SE 13.4%; bootstrap median 0.321 [0.223, 0.395]), read as an SD

    # No IIV on V/F or ka: Table 3 reports a single omega row, for CL/F only.

    # Residual error. Methods Equation (2) is T_i = U_i + U_i*eps_1 + eps_2,
    # a combined proportional-plus-additive model, and defines eps as having
    # "variance sigma^2". Table 3's rows are labelled `sigma_1` and `sigma_2`,
    # so both printed values are STANDARD DEVIATIONS on the same reading used
    # for omega above.
    propSd <- 0.288;  label("Proportional residual error (fraction)")  # Table 3: sigma_1 = 0.288 (SE 10.1%; bootstrap median 0.278 [0.219, 0.348])
    addSd  <- 3.701;  label("Additive residual error (ng/mL)")         # Table 3: sigma_2 = 3.701 (SE 17.4%; bootstrap median 3.742 [0.703, 4.822])
  })

  model({
    # 1. Individual parameters. Results Equations (6) and (7):
    #      CL/F = 27.6 * (weight/70)^0.75 * (1 - 0.392 * ARI)
    #      V/F  = 854  * (weight/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_ari_cl * CONMED_ARIPIPRAZOLE)
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # 2. Micro-constants
    kel <- cl / vc

    # 3. ODE system -- one-compartment with first-order oral absorption.
    #    F is not separately identifiable from a TDM dataset with oral dosing
    #    only, so CL and vc are the apparent (CL/F, V/F) quantities and no
    #    bioavailability term is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Observation. Doses are in mg and vc is in L, giving mg/L; the factor
    #    of 1000 converts to the ng/mL in which the paper reports olanzapine
    #    concentrations (therapeutic range 20-80 ng/mL, Methods, Simulation).
    Cc <- 1000 * central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
