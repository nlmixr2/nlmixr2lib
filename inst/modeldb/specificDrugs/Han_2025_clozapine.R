Han_2025_clozapine <- function() {
  description <- "One-compartment population PK model for oral clozapine with first-order absorption in adults with schizophrenia, built from routine trough therapeutic-drug-monitoring concentrations (Han 2025). Apparent oral clearance is allometrically scaled on body weight and reduced by 25.4% when zopiclone is co-administered; the absorption rate constant is fixed to a published value. Between-subject variability was retained on CL/F only."
  reference <- paste(
    "Han HH, Zhang Y, Wang J, Tian X, Li Y, He SM, Zhang C, Chen X, Wang DD.",
    "Population pharmacokinetics modelling to predict DDI from zopiclone on",
    "clozapine in schizophrenia patients.",
    "Frontiers in Psychiatry. 2025;16:1664678.",
    "doi:10.3389/fpsyt.2025.1664678.",
    "Final model Equations (6) and (7); parameter estimates Table 3.",
    "The fixed absorption rate constant is quoted from the paper's references 25",
    "and 26: Li LJ, Shang DW, Li WB, et al. Population pharmacokinetics of",
    "clozapine and its primary metabolite norclozapine in Chinese patients with",
    "schizophrenia. Acta Pharmacol Sin. 2012;33:1409-1416. doi:10.1038/aps.2012.71;",
    "see modellib('Li_2012_clozapine'). And Shang DW, Li LJ, Wang XP, et al.",
    "Population pharmacokinetic/pharmacodynamic model of clozapine for",
    "characterizing the relationship between accumulated exposure and PANSS",
    "scores in patients with schizophrenia. Ther Drug Monit. 2014;36:378-386.",
    "doi:10.1097/FTD.0000000000000014.",
    sep = " "
  )
  vignette <- "Han_2025_clozapine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "clozapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clozapine", units = "mg", specimen = "plasma", verified = TRUE)
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
        "(Methods, Equation 3, and Results, Equations 6-7). Cohort mean 70.49 kg (SD 13.53),",
        "median 71.00, range 38.00-120.00 (Table 1). Treated as a baseline, time-fixed",
        "covariate: the paper draws weight from a retrospective medical log and reports a",
        "single value per patient. Note the paper's standard weight of 70 kg is very close to",
        "the cohort median of 71 kg, so the reference subject is close to a typical patient."
      ),
      source_name        = "weight"
    ),
    CONMED_ZOPICLONE = list(
      description        = "Concomitant zopiclone indicator (1 = patient co-prescribed zopiclone tablets, 0 = not)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant zopiclone)",
      notes              = paste(
        "The paper's `ZOP` variable (Results, following Equation 7: 'ZOP denoted zopiclone and",
        "when schizophrenia patients took ZOP, ZOP denoted 1, otherwise ZOP denoted 0').",
        "8 of the 81 patients were on zopiclone tablets (Table 2). This was the only one of the",
        "28 screened concomitant drugs retained by the covariate search. The paper attributes",
        "the effect to competition for CYP3A4, of which both drugs are substrates (Discussion);",
        "it is an empirical exposure shift, not a mechanistically parameterised interaction.",
        "Time-fixed at the analysis baseline in this source."
      ),
      source_name        = "ZOP"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a potential covariate but not retained; cohort mean 49.46 years (SD 11.15), median 50.67, range 20.67-73.11 (Table 1)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained; 44 women / 37 men (Table 1). The paper reports 'Gender (men/women)', so a SEXM-oriented source column would need value inversion."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 39.46 g/L (SD 3.21), median 39.40, range 27.90-47.90 (Table 1)."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 25.47 IU/L (SD 21.61), median 20.00, range 4.00-162.00 (Table 1)."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 21.46 IU/L (SD 12.58), median 19.00, range 9.00-119.00 (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained; cohort mean 61.14 (SD 11.37), median 60.00, range",
        "32.00-96.00 (Table 1). Table 1 tags this row 'Creatinine (mmol/L)', which is a",
        "unit error in the source: 61 mmol/L is roughly a thousand-fold above any survivable",
        "serum creatinine, whereas 61 umol/L is squarely normal and matches the quoted range.",
        "Recorded here in umol/L. The covariate is documentation-only -- it does not enter the",
        "final model -- so the correction cannot affect any prediction."
      )
    ),
    BUN = list(
      description = "Serum urea",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.54 mmol/L (SD 1.39), median 4.33, range 1.82-11.71 (Table 1). The paper's Table 1 row is labelled 'Urea (mmol/L)'."
    ),
    TPRO = list(
      description = "Total protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 65.74 g/L (SD 4.71), median 66.30, range 51.90-76.60 (Table 1)."
    ),
    TCHOL = list(
      description = "Total cholesterol",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 4.08 mmol/L (SD 0.85), median 4.04, range 2.27-6.51 (Table 1)."
    ),
    TRIG = list(
      description = "Serum triglyceride",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 1.58 mmol/L (SD 0.81), median 1.45, range 0.44-5.11 (Table 1)."
    ),
    DBIL = list(
      description = "Direct (conjugated) bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 2.33 (SD 1.18), median 2.00, range 0.50-8.30 (Table 1). Table 1 tags the row 'Direct bilirubin (mmol/L)'; as for creatinine the SI unit is umol/L, and the quoted range is a normal direct-bilirubin range in umol/L."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 7.35 (SD 3.27), median 6.60, range 2.70-21.30 (Table 1). Table 1 tags the row 'Total bilibrubin (mmol/L)' (sic); as for creatinine the SI unit is umol/L."
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 37.98% (SD 3.52), median 37.40, range 31.40-49.20 (Table 1)."
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained; cohort mean 124.93 g/L (SD 14.53), median 124.00, range 21.00-166.00 (Table 1). The lower bound of 21.00 g/L is implausible for a living patient and is very likely a transcription slip in the source; the covariate is documentation-only and does not enter the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 81,
    n_studies      = 1,
    age_mean       = "49.46 years (SD 11.15)",
    age_median     = "50.67 years",
    age_range      = "20.67-73.11 years",
    weight_mean    = "70.49 kg (SD 13.53)",
    weight_median  = "71.00 kg",
    weight_range   = "38.00-120.00 kg",
    sex_female_pct = 54.3,
    disease_state  = "Schizophrenia; inpatients on routine oral clozapine therapy.",
    dose_range     = "Not reported. Concentrations came from routine therapeutic drug monitoring and the paper tabulates no administered clozapine doses, dosing frequencies or sampling times; it states only that dosing 'was mainly based on the instruction' (Methods, Data collection).",
    regions        = "China (single centre: Xuzhou Oriental Hospital Affiliated to Xuzhou Medical University, Jiangsu)",
    notes          = paste(
      "Retrospective analysis of clozapine therapeutic-drug-monitoring concentrations",
      "collected between December 2023 and November 2024 (Methods, Data collection).",
      "All observations are TROUGH concentrations: 'The sample extraction times for plasma",
      "concentrations were before the next administration, which was the value of the trough",
      "concentration.' Clozapine was assayed by homogeneous enzyme immunoassay.",
      "Baseline demographics and laboratory values are Table 1; concomitant medication",
      "counts are Table 2 (28 drugs).",
      "The covariate search screened the demographic and clinical-chemistry variables of",
      "Table 1 -- age, sex, weight, albumin, globulin, alanine transaminase, aspartate",
      "transaminase, creatinine, urea, total protein, total cholesterol, triglyceride,",
      "direct and total bilirubin, hematocrit, hemoglobin, mean corpuscular hemoglobin and",
      "mean corpuscular hemoglobin concentration -- together with the 28 concomitant drugs",
      "of Table 2 (acarbose, alprazolam, amisulpride, amlodipine, aripiprazole, atorvastatin,",
      "bezafibrate, clonazepam, aspirin, glimepiride, lamotrigine, lithium carbonate,",
      "lorazepam, metformin, metoprolol, nifedipine, paliperidone, perphenazine, benzhexol,",
      "propranolol, risperidone oral liquid, risperidone tablets, sertraline, sodium",
      "valproate, sulpiride, valsartan, ziprasidone, zopiclone).",
      "Only body weight and concomitant zopiclone survived; a two-step method was used to",
      "build the covariate model (Methods, Model building). The per-covariate hypothesis",
      "tests are not tabulated in the paper.",
      "Globulin, mean corpuscular hemoglobin and mean corpuscular hemoglobin concentration",
      "have no canonical covariate-column entry and so are recorded in this note rather than",
      "in covariatesDataExcluded; their Table 1 summaries are 26.28 g/L (SD 3.23),",
      "30.17 pg (SD 1.45) and 329.93 g/L (SD 7.89) respectively.",
      "The model was qualified by goodness-of-fit plots, individual plots, a VPC (Figures 1",
      "and 2) and a bootstrap (Table 3)."
    )
  )

  ini({
    # Structural parameters -- Results, Equations (6) and (7), point estimates
    # also in Table 3. Reference subject: 70 kg, no concomitant zopiclone.
    lcl <- log(29.6); label("Apparent oral clearance CL/F at 70 kg without zopiclone (L/h)")  # Table 3 / Eq. 6: CL/F = 29.6 L/h (SE 6.5%; bootstrap median 29.4 [26.1, 33.9])
    lvc <- log(308);  label("Apparent central volume of distribution V/F at 70 kg (L)")       # Table 3 / Eq. 7: V/F = 308 L (SE 14.4%; bootstrap median 309 [230, 421])

    # Ka was not estimated here. Methods, Model building: "CL/F, V/F, and Ka
    # [fixed at 1.3/h (25, 26)] were the main pharmacokinetic parameters";
    # Table 3 prints "1.3 (fixed)" with no SE and no bootstrap row. References
    # 25 and 26 are Li 2012 and Shang 2014, both rich-sampling clozapine popPK
    # analyses in Chinese patients with schizophrenia. The same fixed value is
    # carried by modellib('Li_2012_clozapine').
    lka <- fixed(log(1.3)); label("First-order oral absorption rate constant ka, from Li 2012 / Shang 2014 (1/h)")  # Table 3: Ka = 1.3 1/h (fixed)

    # Allometric exponents. Methods, Equation (3): "W denoted allometric
    # coefficients: 0.75 and 1 for CL/F and V/F, respectively (27)", citing
    # Anderson & Holford (2008). Both are asserted, not estimated -- neither
    # appears in Table 3 and neither carries an SE or a bootstrap interval.
    # Results Equations (6) and (7) print the assembled terms (weight/70)^0.75
    # and (weight/70).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F for WT/70, from Anderson & Holford 2008 (unitless)")  # Methods Eq. 3 and Results Eq. 6 exponent 0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on V/F for WT/70, from Anderson & Holford 2008 (unitless)")   # Methods Eq. 3 and Results Eq. 7 exponent 1

    # Concomitant-zopiclone effect on CL/F, entered in the linear-shift form of
    # Methods Equation (5) for categorical covariates, R_i = TV(R)*(1 + Q*S_i).
    # Results Equation (6) prints the assembled term as (1 - 0.254 * ZOP), i.e.
    # a 25.4% reduction in apparent oral clearance; the Results and Discussion
    # both state the same result in words ("the clozapine clearance of
    # schizophrenia patients was reduced by 25.4%", Figure 3).
    e_zop_cl <- -0.254; label("Proportional change in CL/F with concomitant zopiclone (fraction)")  # Table 3: theta_ZOP = -0.254 (SE 30.8%; bootstrap median -0.241 [-0.408, -0.009])

    # Between-subject variability. Methods Equation (1) is Z_i = TV(Z)*exp(eta_i),
    # i.e. exponential (log-normal) IIV. Table 3's row is labelled `omega_CL/F`
    # (not omega^2), so the printed 0.348 is the STANDARD DEVIATION and the
    # variance is 0.348^2. See the vignette Errata: the SD reading reproduces
    # every cell of the paper's Table 4 safety column to within ~2.5 percentage
    # points, whereas reading 0.348 as the variance (SD 0.590) roughly doubles
    # each one (19-22% against the paper's 9.0-13.0%) and destroys the monotone
    # decline across weight bands that the paper reports. The 11.3% relative
    # standard error is also typical of an omega SD rather than a variance.
    # 0.348 on the log scale is 35.9% CV.
    etalcl ~ 0.348^2  # Table 3: omega_CL/F = 0.348 (SE 11.3%; bootstrap median 0.342 [0.264, 0.422]), read as an SD

    # No IIV on V/F or ka: Table 3 reports a single omega row, for CL/F only.

    # Residual error. Methods Equation (2) is Y_i = X_i + X_i*eps_1, a purely
    # proportional model with no additive term. Table 3's row is labelled
    # `sigma_1` and the footnote calls it "residual variability, proportional
    # error", so the printed value is a STANDARD DEVIATION on the same reading
    # used for omega above -- a 25.7% proportional error.
    propSd <- 0.257; label("Proportional residual error (fraction)")  # Table 3: sigma_1 = 0.257 (SE 6.6%; bootstrap median 0.254 [0.220, 0.289])
  })

  model({
    # 1. Individual parameters. Results Equations (6) and (7):
    #      CL/F = 29.6 * (weight/70)^0.75 * (1 - 0.254 * ZOP)
    #      V/F  = 308  * (weight/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (1 + e_zop_cl * CONMED_ZOPICLONE)
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
    #    of 1000 converts to the ng/mL in which the paper reports clozapine
    #    concentrations (therapeutic range 350-800 ng/mL with a 1000 ng/mL
    #    toxicity threshold, Methods, Dosage simulation).
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd)
  })
}
