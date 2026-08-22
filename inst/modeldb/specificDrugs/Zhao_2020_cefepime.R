Zhao_2020_cefepime <- function() {
  description <- "One-compartment IV population PK model for cefepime in 85 Chinese neonates and young infants with postmenstrual age below 48 weeks (Zhao 2020); body-weight allometric scaling referenced to 3.352 kg (fixed exponents 0.75 on CL, 1 on Vc), a power effect of postmenstrual age on CL referenced to 40 weeks, and an inverse power effect of serum creatinine on CL referenced to 28.5 umol/L."
  reference <- paste(
    "Zhao Y, Yao BF, Kou C, Xu HY, Tang BH, Wu YE, Hao GX, Zhang XP,",
    "Zhao W. Developmental population pharmacokinetics and dosing",
    "optimization of cefepime in neonates and young infants.",
    "Front Pharmacol. 2020;11:14. doi:10.3389/fphar.2020.00014.",
    "Used as 'pharmacometric model 3' (sensitivity analysis C) in",
    "Gotta V, Csajka C, Glauser A, Berger C, Pfister M, Paioni P.",
    "Risk of potentially neurotoxic exposure in infants under high-dose",
    "cefepime treatment - a pharmacometric simulation study.",
    "Pharmaceutics. 2025;17(5):544. doi:10.3390/pharmaceutics17050544",
    "(Supplemental Methods S1).",
    sep = " "
  )
  vignette <- "Gotta_2025_cefepime_infant_neurotoxicity"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Current body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling referenced to 3.352 kg (the cohort median",
        "current weight, reported as 3353 g in Zhao 2020 Table 1 and",
        "written as 3,352 g in the Table 3 model expressions) on both CL",
        "(exponent 0.75) and Vc (exponent 1). Zhao 2020 Results,",
        "'Covariate Analysis': the allometric size approach was applied",
        "a priori with coefficients fixed at 0.75 for CL and 1 for V.",
        "Zhao 2020 Table 3 writes weight in grams (CW/3,352); this model",
        "uses the canonical WT in kg with the equivalent reference",
        "3.352 kg. Cohort range 950-4350 g (0.95-4.35 kg)."
      ),
      source_name        = "CW"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL: (PMA_weeks / 40)^1.16 per Zhao 2020 Table 3",
        "(F_age = (PMA/40)^theta_3, theta_3 = 1.16). The source paper",
        "reports PMA in weeks; this model converts the canonical PAGE",
        "(months) back to weeks via pma_wk = PAGE * 4.35 so the Table 3",
        "constant applies unchanged. This is the same conversion used in",
        "Shoji_2016_cefepime.R. Cohort PMA range 30.6-45.1 weeks",
        "(median 40.1); the study inclusion criterion was PMA < 48 weeks,",
        "so extrapolation much beyond 48 weeks is outside the",
        "model-development range."
      ),
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Inverse power effect on CL. Zhao 2020 Table 3 writes the renal",
        "function factor as RF = 1/(CREA/28.5)^theta_4 with",
        "theta_4 = 0.218, i.e. (CREAT / 28.5)^-0.218. The reference",
        "28.5 umol/L is the cohort median (Table 1; mean 34.3, range",
        "11.5-92.4).",
        "SIGN ERRATUM: Gotta 2025 Supplemental Methods S1 reprints this",
        "term as (Scr/28.5)^0.218 -- dropping the reciprocal and thereby",
        "inverting the direction of the effect. The negative exponent",
        "encoded here is the primary's, and is confirmed by the Zhao 2020",
        "abstract: 'The clearance of cefepime increased with current",
        "weight and decreased with increased serum creatinine",
        "concentration'. A positive exponent would make clearance rise",
        "with creatinine, which is backwards for a renally eliminated",
        "cephalosporin."
      ),
      source_name        = "CREA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 85L,
    n_studies      = 1L,
    age_range      = "Postnatal age 1-25 days (mean 7.58, SD 3.83, median 8)",
    age_median     = "8 days postnatal age",
    weight_range   = "0.95-4.35 kg (mean 3.21, SD 0.678, median 3.353)",
    weight_median  = "3.353 kg",
    disease_state  = "Neonates and young infants treated with cefepime for neonatal infection",
    dose_range     = "30 mg/kg every 12 h intravenously; actual doses 30-190 mg/dose (mean 106, SD 31.8), 25.2-53.9 mg/kg/dose (mean 33.3, SD 8.31)",
    regions        = "China (Beijing Obstetrics and Gynecology Hospital and Shandong Provincial Qianfoshan Hospital)",
    ga_range       = "28.0-41.6 weeks (mean 38.1, SD 2.80, median 39.0)",
    pma_range      = "30.6-45.1 weeks (mean 39.2, SD 3.35, median 40.1)",
    creat_range    = "11.5-92.4 umol/L (mean 34.3, SD 17.1, median 28.5)",
    n_observations = 100L,
    notes          = paste(
      "Baseline demographics per Zhao 2020 Table 1. 85 neonates enrolled",
      "2017-2018 contributing 100 plasma concentrations under an",
      "opportunistic sampling design. Inclusion criteria: PMA < 48 weeks,",
      "regular cefepime treatment, written parental consent. Exclusion",
      "criteria: concomitant other antibiotics, expected survival shorter",
      "than the treatment cycle, or investigator-determined unsuitability.",
      "Final model fit in NONMEM with an exponential inter-individual",
      "variability model and a proportional residual-error model;",
      "evaluated by goodness-of-fit plots, a bootstrap (Table 3 medians",
      "and 5th-95th percentiles) and NPDE, plus an external validation.",
      "Birth weight mean 3092 g (SD 620, range 980-4210)."
    )
  )

  ini({
    # Structural PK -- Zhao 2020 Table 3 final estimates (full dataset).
    #   V  (L)   = theta_1 * (CW / 3352)
    #   CL (L/h) = theta_2 * (CW / 3352)^0.75 * F_age * RF
    lvc <- log(2.07);  label("Typical Vc at WT = 3.352 kg (L)")                                                # Zhao 2020 Table 3: theta_1 = 2.07, RSE 8.40% (bootstrap 2.06 [1.79-2.46])
    lcl <- log(0.589); label("Typical CL at WT = 3.352 kg, PMA = 40 weeks and CREAT = 28.5 umol/L (L/h)")      # Zhao 2020 Table 3: theta_2 = 0.589, RSE 6.20% (bootstrap 0.586 [0.530-0.649])

    # Allometric exponents. Zhao 2020 Results, 'Covariate Analysis': "The
    # allometric size approach was used by incorporating a priori the current
    # weight into the basic model (allometric coefficients of 0.75 for CL, 1
    # for V)" -- fixed a priori, not estimated, and correspondingly absent
    # from the Table 3 theta list.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)") # Zhao 2020 Results, Covariate Analysis (fixed a priori)
    e_wt_vc <- fixed(1.00); label("Allometric exponent on Vc (unitless)") # Zhao 2020 Results, Covariate Analysis (fixed a priori)

    # Covariate effects on CL (Zhao 2020 Table 3).
    e_page_cl  <- 1.16; label("Power exponent on (PMA / 40 weeks) for CL (unitless)")   # Zhao 2020 Table 3: theta_3 = 1.16, RSE 49.5% (bootstrap 1.21 [0.283-2.042])
    # NEGATIVE by construction: Table 3 gives RF = 1/(CREA/28.5)^theta_4 with
    # theta_4 = +0.218, so the exponent applied to the ratio is -0.218. See
    # covariateData$CREAT$notes for the Gotta 2025 sign erratum.
    e_creat_cl <- -0.218; label("Power exponent on (CREAT / 28.5 umol/L) for CL (unitless)") # Zhao 2020 Table 3: theta_4 = 0.218 inside RF = 1/(CREA/28.5)^theta_4, RSE 45.4% (bootstrap 0.238 [0.068-0.363])

    # Reference (centring) covariate values -- structural constants of the
    # published equations, not estimated quantities.
    wt_ref    <- fixed(3.352); label("Reference current body weight (kg)")     # Zhao 2020 Table 3 denominator (3,352 g); Table 1 cohort median 3353 g
    page_ref  <- fixed(40);    label("Reference postmenstrual age (weeks)")    # Zhao 2020 Table 3: F_age = (PMA/40)^theta_3
    creat_ref <- fixed(28.5);  label("Reference serum creatinine (umol/L)")    # Zhao 2020 Table 3: RF = 1/(CREA/28.5)^theta_4; Table 1 cohort median 28.5

    # Inter-individual variability -- Zhao 2020 Table 3 "Inter-individual
    # variability (%)": CL 15.3, V 26.8, under the exponential IIV model
    # declared in Results ("The exponential model best described the
    # inter-individual variability"). Encoded with the standard NONMEM
    # convention that the reported percentage is omega x 100 (the
    # approximate CV), so omega^2 = (pct/100)^2:
    #   CL: 0.153^2 = 0.023409
    #   V : 0.268^2 = 0.071824
    # The alternative lognormal reading omega^2 = log(1 + CV^2) would give
    # 0.023139 and 0.069353 -- 1.2% and 3.6% smaller in variance, which is
    # immaterial against these parameters' own uncertainty (RSE 87.6% and
    # 56.1%, bootstrap 5th-95th 3.88-24.2% and 7.75-35.3%). The paper prints
    # neither an explicit exponential-IIV term in its equations nor an omega
    # covariance block, so the two readings cannot be separated arithmetically;
    # the choice is documented rather than hidden.
    etalcl ~ 0.023409 # Zhao 2020 Table 3: IIV CL 15.3%, RSE 87.6% (bootstrap 15.5 [3.88-24.2])
    etalvc ~ 0.071824 # Zhao 2020 Table 3: IIV V 26.8%, RSE 56.1% (bootstrap 23.8 [7.75-35.3])

    # Residual error -- Zhao 2020 Table 3 "Residual variability (%) 36.6",
    # Results: "Residual variability was best described using a proportional
    # model."
    propSd <- 0.366; label("Proportional residual error (fraction)") # Zhao 2020 Table 3: 36.6%, RSE 20.6% (bootstrap 35.5 [28.4-47.7])
  })

  model({
    # ----- Derived covariate terms -----
    # Convert canonical PAGE (months) back to source-paper PMA (weeks) so the
    # Zhao 2020 Table 3 constants apply unchanged.
    pma_wk <- PAGE * 4.35

    # Allometric size on CL and Vc.
    wt_cl_factor <- (WT / wt_ref)^e_wt_cl
    wt_vc_factor <- (WT / wt_ref)^e_wt_vc

    # F_age = (PMA / 40)^1.16 (Zhao 2020 Table 3).
    page_cl_factor <- (pma_wk / page_ref)^e_page_cl

    # RF = 1 / (CREAT / 28.5)^0.218 = (CREAT / 28.5)^-0.218 (Zhao 2020
    # Table 3). Clearance falls as serum creatinine rises.
    creat_cl_factor <- (CREAT / creat_ref)^e_creat_cl

    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * wt_cl_factor * page_cl_factor * creat_cl_factor
    vc <- exp(lvc + etalvc) * wt_vc_factor

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system (one-compartment IV) -----
    d/dt(central) <- -kel * central

    # ----- Observation -----
    # Dose in mg, vc in L -> central / vc has units mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
