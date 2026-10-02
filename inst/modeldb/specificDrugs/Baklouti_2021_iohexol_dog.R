Baklouti_2021_iohexol_dog <- function() {
  description <- paste(
    "Clinical veterinary (dog).",
    "Two-compartment population PK model for the glomerular filtration rate",
    "marker iohexol given as a single 64.7 mg/kg intravenous bolus to 49",
    "client-owned dogs, 29 healthy and 20 with chronic kidney disease (IRIS",
    "criteria). Parameterised per kg body weight in terms of plasma clearance,",
    "intercompartmental clearance and the two volumes of distribution; because",
    "iohexol is eliminated only by glomerular filtration, the individual",
    "clearance is the dog's GFR. Clearance falls exponentially with centred",
    "serum creatinine and is 32% lower in dogs with chronic kidney disease.",
    "Log-normal between-subject variability on clearance and both volumes",
    "(none on Q) and a proportional residual error (Baklouti 2021).",
    sep = " "
  )
  reference <- paste(
    "Baklouti S, Concordet D, Borromeo V, Pocar P, Scarpa P, Cagnardi P.",
    "Population pharmacokinetic model of iohexol in dogs to estimate",
    "glomerular filtration rate and optimize sampling time.",
    "Front Pharmacol. 2021;12:634404. doi:10.3389/fphar.2021.634404.",
    "All parameter values are the final estimates of Table 2; the covariate",
    "equation is the display equation of Results, 'Population",
    "Pharmacokinetics'.",
    sep = " "
  )
  vignette <- "Baklouti_2021_iohexol_dog"

  # Every structural parameter is published per kg body weight (CL and Q in
  # L/min/kg, V1 and V2 in L/kg; Table 2) and the dose is nominal 64.7 mg/kg
  # (Methods, 'Sample Collection and Analysis'). Doses are therefore given in
  # mg/kg and the compartment amounts are carried in mg/kg, so central/vc is
  # mg/L -- numerically the ug/mL of the HPLC assay. Time is in minutes, the
  # unit of Table 2 and of every sampling time.
  units <- list(time = "min", dosing = "mg/kg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "iohexol",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "iohexol",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters log(CL) linearly after centring: ln(Cl) = ln(theta_Cl) +",
        "theta_1 * 1[diseased dogs] + theta_2 * ccreatinine (mg/dL) + eta_Cl",
        "(Results, 'Population Pharmacokinetics'), where 'ccreatinine is the",
        "creatinine centered around an average value'. The centring value is",
        "not printed; 1.47 mg/dL, the cohort mean of Table 1 ('Serum",
        "creatinine (mg/dl) 1.47 +/- 1.99'), is used because the paper calls",
        "it an average, and because it reproduces the paper's one-sample",
        "sampling-time analysis (Table 3) where the median (1.09 mg/dL) does",
        "not; see the vignette. Table 1 range 0.67-14.4 mg/dL. Units are",
        "mg/dL, not umol/L.",
        sep = " "
      ),
      source_name = "creatinine"
    ),
    DIS_RENAL = list(
      description = "Chronic kidney disease indicator (IRIS criteria)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy kidney status, CKD-)",
      notes = paste(
        "Table 1 codes renal status as 'Healthy CKD - n = 29 (Code 0)' and",
        "'Diseased CKD + n = 20 (Code 1)', and the covariate equation uses",
        "1[diseased dogs] = 1 for CKD+ dogs, so DIS_RENAL equals the paper's",
        "code with no inversion. Dogs were classed as healthy or as having",
        "chronic kidney disease by the International Renal Interest Society",
        "(IRIS) guidelines on physical examination, blood count, serum",
        "biochemistry, urinalysis, UPC ratio and ultrasound (Methods,",
        "'Animals'). Chronic kidney disease is the only renal diagnosis in",
        "the cohort, so the pooled DIS_RENAL flag carries it exactly.",
        sep = " "
      ),
      source_name = "kidney status (CKD)"
    )
  )

  # Seven further covariates were screened (Methods, 'Population
  # Pharmacokinetics': 'Nine covariates were tested on individual PK
  # parameters') and not retained. None has a published point estimate.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened and not retained. The disposition parameters are already",
        "per kg, so a per-kg dose carries the dominant weight dependence.",
        "Table 1: 25.8 +/- 9.5 kg, range 3.9-46 kg, median 27.6 kg.",
        sep = " "
      ),
      source_name = "Body weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened and not retained; the Discussion attributes this to age",
        "not being interpretable without breed or body size. Table 1:",
        "5.43 +/- 3.5 years, range 0.4-16 years, median 4.5 years.",
        sep = " "
      ),
      source_name = "Age"
    ),
    UREA = list(
      description = "Serum urea",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened and not retained: the Discussion states that the part of",
        "the clearance variability explained by urea is the same as that",
        "explained by creatinine. This is urea, not blood urea nitrogen, so",
        "the BUN canonical does not apply (urea mg/dL is about 2.14 times BUN",
        "mg/dL). No canonical register entry exists; the name is",
        "documentation only. Table 1: 44.95 +/- 34.81 mg/dL, range",
        "16-181.4 mg/dL, median 36 mg/dL.",
        sep = " "
      ),
      source_name = "Serum urea"
    ),
    USG = list(
      description = "Urine specific gravity",
      units = "(unitless, x1000)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened and not retained. Table 1 prints it on the x1000 scale:",
        "1,036.78 +/- 18.59, range 1,003-1,065, median 1,040. No canonical",
        "register entry exists (USG_CORRECTED is an unrelated",
        "concentration-correction flag); the name is documentation only.",
        sep = " "
      ),
      source_name = "Urine specific gravity (USG)"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Table 1 codes sex as a three-level factor: male n = 19 (code 0),",
        "female n = 6 (code 1), female neutered n = 24 (code 2). SEXF = 1",
        "pools codes 1 and 2 (30 dogs); the neutering split is carried as",
        "NEUTERED below. Screened and not retained.",
        sep = " "
      ),
      source_name = "Sex"
    ),
    NEUTERED = list(
      description = "Surgical neutering / gonadectomy status",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (sexually intact)",
      notes = paste(
        "The third level of the paper's sex factor (Table 1, 'Female",
        "neutered n = 24 (Code 2)'). Castration status of the 19 males is not",
        "reported. Only ever fitted as a level of the sex factor; not",
        "retained.",
        sep = " "
      ),
      source_name = "Sex (code 2)"
    ),
    BREED_PUREBRED = list(
      description = "Purebred (non-mongrel) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (mongrel)",
      notes = paste(
        "Table 1 codes breed as 'Mongrel n = 10 (Code 0)' and 'Other breeds",
        "n = 39 (Code 1)'; the 39 purebred dogs spanned 21 breeds. Screened",
        "and not retained. No canonical register entry exists; the name",
        "follows Cagnardi_2018_cefazolin_dog and is documentation only.",
        sep = " "
      ),
      source_name = "Breed"
    )
  )

  population <- list(
    species = "dog (client-owned Canis lupus familiaris; 10 mongrels and 39 dogs across 21 pure breeds)",
    n_subjects = 49L,
    n_studies = 1L,
    n_observations = 245L,
    age_range = "0.4-16 years",
    age_median = "4.5 years",
    weight_range = "3.9-46 kg",
    weight_median = "27.6 kg",
    sex_female_pct = 61.2,
    disease_state = paste(
      "29 dogs with healthy kidney status (CKD-) and 20 with chronic kidney",
      "disease (CKD+) by IRIS criteria, all scheduled at a veterinary",
      "teaching hospital for various clinical procedures",
      sep = " "
    ),
    dose_range = "64.7 mg/kg (nominal) iohexol as a single 60-s intravenous bolus",
    regions = "Italy (University Veterinary Teaching Hospital, University of Milan)",
    renal_function = "Serum creatinine 1.47 +/- 1.99 mg/dL, range 0.67-14.4 mg/dL, median 1.09 mg/dL; serum urea 44.95 +/- 34.81 mg/dL (Table 1)",
    notes = paste(
      "Demographics from Table 1 and Results, 'Animals and Iohexol",
      "Concentrations'. Five plasma samples per dog at 5, 15, 60, 90 and",
      "180 min (Methods and Supplementary File S1), 245 samples in total;",
      "concentrations ranged 16.4-643.6 ug/mL (HPLC, LOQ 1.80 ug/mL). The",
      "Results list the sampling times as '5, 15, 30, 60, 90, and 180 min',",
      "but 245 = 49 x 5 matches the five-time protocol of the Methods. The",
      "exact dose was measured by weighing the syringe (Supplementary File",
      "S1). Sex: 19 male, 6 female, 24 female neutered, so the female",
      "percentage is 30 / 49.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Disposition. Table 2, 'Synthesis of estimates obtained in the final Pop
    # PK model'. Every value is per kg body weight. theta_Cl is the clearance
    # of a CKD- dog at the centring creatinine value.
    # ------------------------------------------------------------------------
    lcl <- log(0.00212)
    label("Plasma (glomerular) clearance CL (L/min/kg)")
    # Table 2, row 'theta Cl' = 0.00212 L/min/kg (SE 0.00010, RSE 4.68%)

    lvc <- log(0.163)
    label("Central volume of distribution V1 (L/kg)")
    # Table 2, row 'theta V1' = 0.163 L/kg (SE 0.00661, RSE 4.07%)

    lvp <- log(0.058)
    label("Peripheral volume of distribution V2 (L/kg)")
    # Table 2, row 'theta V2' = 0.058 L/kg (SE 0.00387, RSE 6.64%)

    lq <- log(0.0034)
    label("Intercompartmental clearance Q (L/min/kg)")
    # Table 2, row 'theta Q' = 0.0034 L/min/kg (SE 0.00042, RSE 12.21%)

    # ------------------------------------------------------------------------
    # Covariate effects on log(CL). Results display equation:
    #   ln(Cl) = ln(theta_Cl) + theta_1 * 1[diseased dogs]
    #            + theta_2 * ccreatinine (mg/dL) + eta_Cl
    # ------------------------------------------------------------------------
    e_dis_renal_cl <- -0.379
    label("Effect of chronic kidney disease on log(CL) (unitless)")
    # Table 2, row 'theta 1 (diseased dogs)' = -0.379 (SE 0.07002, RSE 18.49%);
    # exp(-0.379) = 0.685, i.e. CL 31.5% lower in CKD+ dogs.

    e_creat_cl <- -0.421
    label("Effect of centred serum creatinine on log(CL) (dL/mg)")
    # Table 2, row 'theta 2 (creatinine)' = -0.421 dl/mg (SE 0.05356, RSE 12.72%)

    # ------------------------------------------------------------------------
    # Between-subject variability. Monolix reports omega as the standard
    # deviation of eta; the Discussion reads Table 2 back as 'the
    # interindividual variability was around 20% on all the parameters' and
    # gives RSEs 'for the standard deviation of the random effects'. The
    # variances below are the squares of the Table 2 values. No IIV on Q
    # (Discussion: 'except intercompartmental clearance that was fixed') and
    # no correlations (Results: 'The correlations between individual
    # parameters were low (< 30%) and were therefore not included').
    # ------------------------------------------------------------------------
    etalcl ~ 0.043264 # Table 2, row 'eta Cl' = 0.208 (SD); 0.208^2
    etalvc ~ 0.061504 # Table 2, row 'eta V1' = 0.248 (SD); 0.248^2
    etalvp ~ 0.039601 # Table 2, row 'eta V2' = 0.199 (SD); 0.199^2

    # ------------------------------------------------------------------------
    # Residual error. Results: 'The variability error model was best
    # described by a proportional error ... The value of the residual error
    # was 6.17%.'
    # ------------------------------------------------------------------------
    propSd <- 0.0617
    label("Proportional residual SD (fraction)")
    # Table 2, row 'Residual error epsilon' = 0.0617
  })

  model({
    # ---- Individual parameters (Results covariate equation) -----------------
    # Creatinine is centred on 1.47 mg/dL, the Table 1 cohort mean (the paper
    # says only 'centered around an average value'; see covariateData).
    cl <- exp(lcl + e_dis_renal_cl * DIS_RENAL + e_creat_cl * (CREAT - 1.47) + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    # ---- Micro-constants ----------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Two-compartment disposition, intravenous bolus into central --------
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Observation --------------------------------------------------------
    # central is in mg/kg and vc in L/kg, so central/vc is mg/L == ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
