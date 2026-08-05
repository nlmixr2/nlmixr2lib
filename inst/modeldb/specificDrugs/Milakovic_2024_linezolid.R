Milakovic_2024_linezolid <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination for",
    "intravenous linezolid in nine critically ill adults with",
    "COVID-19-associated acute respiratory distress syndrome (CARDS)",
    "supported by veno-venous extracorporeal membrane oxygenation (vv ECMO),",
    "who received a higher-than-standard 600 mg dose as a 30-min infusion",
    "every 8 h. Between-subject variability is exponential on both clearance",
    "and volume of distribution, estimated as a correlated 2x2 block with a",
    "strong negative CL-Vd covariance; residual variability is proportional.",
    "No covariate was retained: the automated covariate search found none",
    "significant, which the authors attribute to the small, deliberately",
    "homogeneous sample. The model was used for Monte Carlo probability of",
    "target attainment (PTA) and cumulative fraction of response (CFR)",
    "analyses comparing 600 mg every 8 h against the standard 600 mg every",
    "12 h."
  )
  reference <- paste(
    "Milakovic D, Kovacevic T, Kovacevic P, Barisic V, Avram S, Dragic S,",
    "Zlojutro B, Momcicevic D, Miljkovic B, Vucicevic K. (2024).",
    "Population Pharmacokinetic Model of Linezolid and Probability of Target",
    "Attainment in Patients with COVID-19-Associated Acute Respiratory",
    "Distress Syndrome on Veno-Venous Extracorporeal Membrane",
    "Oxygenation-A Step toward Correct Dosing.",
    "Pharmaceutics 16(2):253. doi:10.3390/pharmaceutics16020253"
  )
  vignette <- "Milakovic_2024_linezolid"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # No covariate was retained in the final model (Milakovic 2024 Results:
  # "none of the covariates showed a significant effect on the PK
  # parameters"). The candidate covariates collected and screened are
  # documented below so the provenance of the covariate search is preserved;
  # none of them is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age at ICU admission",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened in the automated covariate search but not retained.",
        "Median 40 years (range 30-62; IQR 30.4-59.8), Milakovic 2024",
        "Table 1."
      ),
      source_name        = "Age"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Collected (Milakovic 2024 Methods 'Data': body weight and height",
        "were recorded) and screened but not retained. Body weight itself",
        "is not tabulated in the paper; only the derived BMI is reported in",
        "Table 1, so no cohort weight distribution can be recovered from the",
        "publication."
      ),
      source_name        = "body weight"
    ),
    BMI = list(
      description        = "Body mass index at ICU admission",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Median 27.7 kg/m^2 (range 23.5-39.2;",
        "IQR 23.7-38.28), Milakovic 2024 Table 1. This is the only body-size",
        "descriptor tabulated in the paper."
      ),
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened but not retained. Milakovic 2024 Table 1 reports 5 male",
        "(55.6%) of the 9 analysed patients, i.e., 4 female (44.4%). The",
        "paper reports the male count, so SEXF = 1 - SEXM."
      ),
      source_name        = "Male (number)"
    ),
    CREAT = list(
      description        = "Plasma creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Reported both at ICU admission (median",
        "70 umol/L, range 53-349) and on the PK sampling day (median 55",
        "umol/L, range 34-201), Milakovic 2024 Table 1. Patients on renal",
        "replacement therapy were excluded by protocol."
      ),
      source_name        = "Plasma creatinine concentration"
    ),
    BUN = list(
      description        = "Serum urea concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Reported as urea (not blood urea",
        "nitrogen): median 5.5 mmol/L (range 4.5-36.3) at admission and 7.2",
        "mmol/L (range 4.5-19.9) on the PK sampling day, Milakovic 2024",
        "Table 1. Urea in mmol/L converts to BUN in mg/dL by multiplying by",
        "2.80."
      ),
      source_name        = "Urea"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Median 32 g/L (range 30-45) at admission",
        "and 38 g/L (range 30-46) on the PK sampling day, Milakovic 2024",
        "Table 1. The Discussion notes that all patients had normal albumin",
        "levels, in contrast with a published ECMO case report where",
        "hypoalbuminaemia was invoked to explain a very low linezolid trough."
      ),
      source_name        = "Albumin"
    ),
    TPRO = list(
      description        = "Total serum protein concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Median 60 g/L (range 53-66) at admission",
        "and 51 g/L (range 47-65) on the PK sampling day, Milakovic 2024",
        "Table 1."
      ),
      source_name        = "Total protein"
    ),
    TBILI = list(
      description        = "Total serum bilirubin concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Median 9.9 umol/L (range 3.9-25.7) at",
        "admission and 37.4 umol/L (range 13.9-156.3) on the PK sampling",
        "day, Milakovic 2024 Table 1."
      ),
      source_name        = "Bilirubin"
    ),
    DDIMER = list(
      description        = "Plasma D-dimer concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Median 1.1 mg/L (range 0.54-49.32),",
        "Milakovic 2024 Table 1. Reported in mg/L rather than the register's",
        "ng/mL reference unit; 1 mg/L = 1000 ng/mL."
      ),
      source_name        = "D-dimer"
    ),
    CRP = list(
      description        = "C-reactive protein concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Listed among the collected clinical variables screened as covariates",
        "(Milakovic 2024 Methods 'Data') but not tabulated in Table 1 and",
        "not retained in the final model."
      ),
      source_name        = "C-reactive protein"
    ),
    AST = list(
      description        = "Aspartate aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Collected as part of the 'liver enzymes' panel screened as",
        "covariates (Milakovic 2024 Methods 'Data'); not tabulated in Table 1",
        "and not retained."
      ),
      source_name        = "liver enzymes"
    ),
    ALT = list(
      description        = "Alanine aminotransferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Collected as part of the 'liver enzymes' panel screened as",
        "covariates (Milakovic 2024 Methods 'Data'); not tabulated in Table 1",
        "and not retained."
      ),
      source_name        = "liver enzymes"
    ),
    SAPS_II = list(
      description        = "Simplified Acute Physiology Score II on ICU admission day",
      units              = "points",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Median 21 points (range 13-36; IQR",
        "13.2-35.2), Milakovic 2024 Table 1."
      ),
      source_name        = "SAPS II"
    ),
    URINE_FLOW = list(
      description        = "Urine output over 24 h",
      units              = "mL/24 h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Collected as '24 h urine output' and screened as a covariate",
        "(Milakovic 2024 Methods 'Data'); not tabulated in Table 1 and not",
        "retained. Reported per 24 h rather than as the register's",
        "instantaneous mL/h rate."
      ),
      source_name        = "24 h urine output"
    ),
    ECMO_PUMP_SPEED = list(
      description        = "ECMO centrifugal pump rotational speed on the PK sampling day",
      units              = "RPM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate whose test statistic the paper reports:",
        "including ECMO centrifugal pump speed in the base model dropped the",
        "objective function value by only 2.04 units, short of the 3.84",
        "required for one degree of freedom at p < 0.05, so it was not",
        "retained (Milakovic 2024 Results). Median 3500 RPM (range 2900-4300;",
        "IQR 2920-4240), Table 1."
      ),
      source_name        = "ECMO centrifugal pump speed"
    ),
    T_ECMO = list(
      description        = "Time on ECMO support at the PK sampling day",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained. Reported as 'ECMO duration', median 6",
        "days (range 1.75-17; IQR 1.86-16.25), Milakovic 2024 Table 1."
      ),
      source_name        = "ECMO duration"
    )
  )

  compartmentData <- list(
    central = list(
      analyte  = "linezolid",
      units    = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 9L,
    n_studies      = 1L,
    n_pk_observations = 53L,
    age_range      = "30-62 years (median 40; IQR 30.4-59.8)",
    age_median     = "40 years",
    weight_range   = paste(
      "not reported; body weight was collected but only the derived BMI is",
      "tabulated (median 27.7 kg/m^2, range 23.5-39.2)"
    ),
    sex_female_pct = 44.4,
    race_ethnicity = "not reported",
    disease_state  = paste(
      "COVID-19-associated acute respiratory distress syndrome (CARDS)",
      "requiring veno-venous extracorporeal membrane oxygenation (vv ECMO).",
      "All patients were RT-PCR-confirmed SARS-CoV-2 positive and met a lung",
      "injury (Murray) score >= 3 or uncompensated hypercapnic acidosis with",
      "pH < 7.2; all 9 had a P/F ratio < 100 and a median Murray score of 3",
      "(range 3-3.8). Median SOFA on the PK sampling day was 10 (range 7-20).",
      "Patients under 18 years, pregnant, allergic to linezolid, or who had",
      "undergone therapeutic plasma exchange within 24 h or renal replacement",
      "therapy were excluded."
    ),
    dose_range     = "600 mg linezolid intravenously as a 30-min infusion every 8 h",
    regions        = paste(
      "Bosnia and Herzegovina (Medical Intensive Care Unit, University",
      "Clinical Centre of the Republic of Srpska, Banja Luka)"
    ),
    notes          = paste(
      "Prospective, observational, single-centre PK study conducted between",
      "1 January and 31 December 2021. Eleven patients were sampled; two were",
      "excluded from the PopPK analysis (one on continuous veno-venous",
      "haemodialysis combined with a CytoSorb device, one who received a blood",
      "transfusion during the sampling window), leaving 9 patients and 53",
      "steady-state serum concentrations. Sampling was rich: pre-dose and 30,",
      "60, 120, 240 and 360 min after the start of infusion, after at least",
      "six doses had been given. Linezolid was assayed by homogeneous enzyme",
      "immunoassay (ARK Linezolid Assay on a Beckman Coulter DxC 700 AU),",
      "measuring range 0.75-30 mg/L. ECMO support used a Maquet Cardiohelp",
      "system with a Quadrox-D oxygenator. Estimation was in NONMEM 7.5 with",
      "PsN 5.2.6 for the VPC (1000 simulations) and 1000 bootstrap runs.",
      "Additional clinical variables tabulated in Table 1 but with no",
      "canonical covariate-register entry, and therefore not listed in",
      "covariatesDataExcluded, are: SOFA score, Murray lung-injury score,",
      "P/F ratio, platelet count, ECMO blood flow rate (median 4.6 L/min),",
      "ECMO sweep gas flow (median 5 L/min) and 24-h fluid balance."
    )
  )

  ini({
    # ---- Structural parameters (Milakovic 2024 Table 3) ----
    # One-compartment model with first-order elimination; linezolid is given
    # intravenously, so there is no absorption or bioavailability term.
    lvc <- log(41.1) ; label("Volume of distribution (L)")   # Milakovic 2024 Table 3 (Volume of distribution 41.1 L, bootstrap 95% CI 30.57-51.77)
    lcl <- log(5.9)  ; label("Clearance (L/h)")              # Milakovic 2024 Table 3 (Clearance 5.9 L/h, bootstrap 95% CI 4.98-7.10)

    # ---- Between-subject variability (Milakovic 2024 Table 3) ----
    # Table 3 reports the interindividual variabilities as percentages
    # (36.3% on Vd, 24.8% on CL) and the Vd-CL covariance on the variance
    # scale (-0.0901). The percentages are omega on the log scale, NOT
    # log-normal CV%: back-calculating the individual CL and Vd implied by
    # Table 2 (CL_i = 1800 / AUC24_i and V_i = t_half_i * CL_i / ln 2)
    # reproduces geometric means of 5.92 L/h and 41.1 L, and population
    # log-scale SDs of 0.249 and 0.358 -- i.e., the reported 24.8% and 36.3%
    # are sd(log CL) and sd(log Vd) directly. Hence omega^2 = 0.248^2 and
    # 0.363^2 with no CV-to-variance conversion.
    #
    # The same back-calculation gives a Vd-CL covariance of -0.0883
    # (correlation -0.991), matching the reported -0.0901. Because Table 3's
    # rounding makes |r| = 1.0008 > 1, the published block as printed is
    # marginally NOT positive definite (determinant -1.37e-05), and rxode2's
    # Cholesky sampler cannot decompose it. Scaling the off-diagonal by 0.99
    # restores positive definiteness and lands on correlation -0.991, which is
    # the value the Table 2 individual estimates actually imply. The reported
    # variances are kept exactly.
    # Block entries below, in order:
    #   0.131769       = 0.363^2, IIV in Vd 36.3% (bootstrap 95% CI 16.11-44.96%)
    #   0.99 * -0.0901 = Vd-CL covariance -0.0901, off-diagonal nudged by 0.99
    #   0.061504       = 0.248^2, IIV in CL 24.8% (bootstrap 95% CI 13.03-31.25%)
    etalvc + etalcl ~ c(
      0.131769,
      0.99 * -0.0901, 0.061504
    )

    # ---- Residual variability (Milakovic 2024 Table 3) ----
    # Table 3's "Proportional residual error 0.114" is the proportional
    # residual SD (11.4%), not a variance: comparing the observed peak and
    # trough concentrations in Table 2 against the model-implied individual
    # predictions built from the same table's AUC24 and half-life gives a
    # proportional residual SD of about 0.076, consistent with 11.4% and
    # incompatible with the 33.8% that reading 0.114 as a variance implies.
    propSd <- 0.114 ; label("Proportional residual error (fraction)")  # Milakovic 2024 Table 3 (bootstrap 95% CI 0.079-0.148)
  })

  model({
    # ---- Individual parameters ----
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # ---- Micro-constant ----
    kel <- cl / vc

    # ---- ODE system: one compartment, first-order elimination ----
    # Doses are 30-min intravenous infusions into central (use rate or dur
    # in the event table).
    d/dt(central) <- -kel * central

    # ---- Observation and error ----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
