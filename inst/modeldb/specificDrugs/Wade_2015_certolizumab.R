Wade_2015_certolizumab <- function() {
  description <- "One-compartment population PK model with first-order SC absorption and an additive baseline concentration for certolizumab pegol in adults with Crohn's disease (Wade 2015)"
  reference <- "Wade JR, Parker G, Kosutic G, Feagen BG, Sandborn WJ, Laveille C, Oliver R. Population pharmacokinetic analysis of certolizumab pegol in patients with Crohn's disease. J Clin Pharmacol. 2015;55(8):866-874. doi:10.1002/jcph.491"
  vignette <- "Wade_2015_certolizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    ADA_POS = list(
      description        = "Anti-certolizumab antibody positivity (per CZP concentration, time-varying)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (anti-CZP antibody negative)",
      notes              = "Wade 2015 codes this time-varying at the per-concentration level (ATB present/absent per CZP concentration, not per subject). A concentration is classified ADA-positive if measured anti-CZP antibody > 2.4 units/mL; negative otherwise. Source column label 'ATB'; renamed to the canonical ADA_POS per covariate-columns.md.",
      source_name        = "ATB"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation effect on both CL/F and V/F: (1 + theta * (BSA - 1.76)). Reference 1.76 m^2 is the typical BSA used in Table 3 footnotes. The Wade 2015 paper does not state the BSA computation formula; assume the formula used in the source studies' case-report forms.",
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Non-linear piecewise-linear effect on CL/F with break at 41 g/L (population median). Form: (1 + theta7 * (ALB - 41)) for ALB <= 41 and (1 + theta8 * (ALB - 41)) for ALB > 41 per Table 3 footnotes. Baseline value (time-fixed).",
      source_name        = "ALB"
    ),
    CRP = list(
      description        = "Baseline C-reactive protein concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Non-linear piecewise-linear effect on CL/F with break at 8 mg/L (population median). Form: (1 + theta10 * (CRP - 8)) for CRP <= 8 and (1 + theta11 * (CRP - 8)) for CRP > 8 per Table 3 footnotes. Standard CRP assay (baseline; Crohn's disease cohort with median CRP well above the hs-CRP sensitivity range).",
      source_name        = "CRP"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White, Indian Asian, Hispanic, or Other (pooled reference per Wade 2015 Table 3 footnotes)",
      notes              = "Multiplicative fractional effect on V/F: (1 + theta14 * RACE_BLACK). Source NONMEM coding: IF(RACE.EQ.2) VRACE = (1 + theta14). N = 29 Black subjects in the model-development cohort (Table 2).",
      source_name        = "RACE"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator (Wade 2015: excludes Japanese and Indian Asian)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White, Indian Asian, Hispanic, or Other (pooled reference per Wade 2015 Table 3 footnotes)",
      notes              = "Multiplicative fractional effect on V/F: (1 + theta15 * RACE_ASIAN). Source NONMEM coding: IF(RACE.EQ.3) VRACE = (1 + theta15). In Wade 2015 coding Japanese (RACE.EQ.8) and Indian Asian (RACE.EQ.4) are separate categories -- the RACE_ASIAN column here captures only RACE.EQ.3 (N = 11 subjects in Table 2). Indian Asian is folded into the reference group; Japanese has its own indicator RACE_JAPANESE.",
      source_name        = "RACE"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese race/ethnicity indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = White, Indian Asian, Hispanic, or Other (pooled reference per Wade 2015 Table 3 footnotes)",
      notes              = "Multiplicative fractional effect on V/F: (1 + theta13 * RACE_JAPANESE). Source NONMEM coding: IF(RACE.EQ.8) VRACE = (1 + theta13). N = 89 Japanese subjects in the model-development cohort (Table 2); the largest non-reference race category.",
      source_name        = "RACE"
    )
  )

  population <- list(
    n_subjects     = 2157L,
    n_studies      = 9L,
    age_range      = "16-80 years",
    age_median     = "35 years",
    weight_range   = "31-151 kg",
    weight_median  = "65 kg",
    bmi_range      = "13-56 kg/m^2",
    bmi_median     = "23 kg/m^2",
    bsa_range      = "1.2-2.7 m^2",
    bsa_median     = "1.8 m^2",
    sex_female_pct = 55.5,
    race_ethnicity = c(White = 91.1, Black = 1.3, Asian = 0.5, Indian = 0.8, Hispanic = 0.2, Japanese = 4.1, Other = 1.9),
    disease_state  = "Moderately to severely active Crohn's disease (active or maintenance treatment populations pooled)",
    dose_range     = "100, 200, or 400 mg SC every 2 weeks (Q2W) or every 4 weeks (Q4W); most data from 400 mg SC Q2W / Q4W regimens",
    regions        = "Multinational (9 pooled studies: C87005, C87031, C87032, C87037, C87042, C87043, C87047, C87048, C87085).",
    baseline_covariates = list(
      crcl_mean_ml_min    = 112,
      crp_mean_mg_L       = 20,
      crp_median_mg_L     = 8,
      albumin_mean_g_L    = 41,
      albumin_median_g_L  = 41,
      lymphocytes_median_10e9_L = 1.5,
      cdai_median         = 290
    ),
    ada_positive_pct_concentrations = 2.0,
    ada_positive_pct_subjects       = 6.4,
    immunosuppressant_use_pct       = 41.2,
    notes          = "Baseline demographics from Wade 2015 Table 2. Population CZP concentration data: 13,561 total observations across 2157 subjects (median 6 samples per subject, range 1-17). Anti-CZP antibodies detected in 139 subjects (6.4%) contributing 270 ADA-positive concentrations (2.0% of all observations). Reference covariate values for typical-subject predictions: BSA = 1.76 m^2, ALB = 41 g/L, CRP = 8 mg/L, ADA-negative, non-Japanese / non-Black / non-Asian (i.e., White, Indian Asian, Hispanic, or Other)."
  )

  ini({
    # Structural parameters (Wade 2015 Table 3, final population PK model; typical
    # ADA-negative subject at the reference covariate values BSA = 1.76 m^2,
    # ALB = 41 g/L, CRP = 8 mg/L, non-Japanese / non-Black / non-Asian).
    lka       <- log(0.200); label("First-order SC absorption rate constant (ka, 1/day)")    # Wade 2015 Table 3 (theta1)
    lcl       <- log(0.685); label("Apparent clearance when ADA-negative (CL/F, L/day)")     # Wade 2015 Table 3 (theta2)
    lvc       <- log(7.61);  label("Apparent volume of distribution (V/F, L)")               # Wade 2015 Table 3 (theta3)
    lbaseline <- log(1.23);  label("Additive baseline concentration (ug/mL)")                # Wade 2015 Table 3 (theta4)

    # ADA-positive fractional-change factor on CL/F (Table 3 footnote:
    # CL/F = ((1 - ATB) * theta2 + ATB * theta5) * ...). theta5 / theta2 gives the
    # relative clearance of ADA-positive vs ADA-negative subjects (2.74 / 0.685 = 4.0x).
    e_ada_cl <- 2.74 / 0.685 - 1;  label("ADA-positive fractional increase in CL/F (unitless)")  # Wade 2015 Table 3 (theta5 = 2.74 L/day; theta2 = 0.685 L/day)

    # Continuous covariate effects (Wade 2015 Table 3 footnotes).
    # Albumin piecewise-linear on CL/F with break at 41 g/L.
    e_alb_lo_cl <- -0.0598;  label("Albumin effect on CL/F for ALB <= 41 g/L (per g/L, centered at 41)")  # Wade 2015 Table 3 (theta7)
    e_alb_hi_cl <- -0.0177;  label("Albumin effect on CL/F for ALB >  41 g/L (per g/L, centered at 41)")  # Wade 2015 Table 3 (theta8)
    # CRP piecewise-linear on CL/F with break at 8 mg/L.
    e_crp_lo_cl <-  0.0205;   label("CRP effect on CL/F for CRP <= 8 mg/L (per mg/L, centered at 8)")     # Wade 2015 Table 3 (theta10)
    e_crp_hi_cl <-  0.000561; label("CRP effect on CL/F for CRP >  8 mg/L (per mg/L, centered at 8)")     # Wade 2015 Table 3 (theta11)
    # BSA linear-deviation effects on CL/F and V/F with reference 1.76 m^2.
    e_bsa_cl    <-  0.715;    label("BSA linear-deviation effect on CL/F (per m^2, centered at 1.76)")     # Wade 2015 Table 3 (theta9)
    e_bsa_vc    <-  0.656;    label("BSA linear-deviation effect on V/F (per m^2, centered at 1.76)")      # Wade 2015 Table 3 (theta12)
    # Race fractional-change effects on V/F (reference = White / Indian Asian / Hispanic / Other).
    e_japanese_vc <- -0.250;  label("Japanese fractional change in V/F vs. pooled reference")               # Wade 2015 Table 3 (theta13)
    e_black_vc    <-  0.265;  label("Black fractional change in V/F vs. pooled reference")                  # Wade 2015 Table 3 (theta14)
    e_asian_vc    <-  0.415;  label("Asian fractional change in V/F vs. pooled reference (excludes Japanese / Indian Asian)")  # Wade 2015 Table 3 (theta15)

    # Inter-individual variability. Wade 2015 reported IIV as %CV; the internal
    # variance on the log scale is omega^2 = log(1 + CV^2). The paper reports
    # two stratified IIVs for CL/F (ADA-negative 27.5%, ADA-positive 83.6%) and
    # for baseline (ADA-negative 95.1%, ADA-positive 72.3%). nlmixr2 expresses a
    # single eta per parameter; we use the ADA-negative values below because
    # they describe the reference-population variability and the ADA-positive
    # subjects are a minority (6.4% of patients, 2.0% of concentrations). The
    # ADA-positive IIV is documented here for the validation vignette.
    etalka       ~ log(1 + 0.501^2)  # Wade 2015 Table 3: IIV KA CV = 50.1%
    etalcl       ~ log(1 + 0.275^2)  # Wade 2015 Table 3: IIV CL/F (ADA-) CV = 27.5% (ADA+ CV = 83.6%)
    etalvc       ~ log(1 + 0.167^2)  # Wade 2015 Table 3: IIV V/F CV = 16.7%
    etalbaseline ~ log(1 + 0.951^2)  # Wade 2015 Table 3: IIV baseline (ADA-) CV = 95.1% (ADA+ CV = 72.3%)

    # Residual error. Wade 2015 used a proportional error model described as
    # "additive on the Ln scale" with theta6 = 34.6%; this maps to a
    # proportional error in linear concentration space for nlmixr2.
    propSd <- 0.346;  label("Proportional residual error (fraction)")   # Wade 2015 Table 3 (theta6)
  })

  model({
    # 1. Covariate-effect terms.
    # Piecewise-linear ALB effect on CL/F with break at 41 g/L.
    alb_cl <- 1 + e_alb_lo_cl * (ALB - 41) * (ALB <= 41) +
                  e_alb_hi_cl * (ALB - 41) * (ALB >  41)
    # Piecewise-linear CRP effect on CL/F with break at 8 mg/L.
    crp_cl <- 1 + e_crp_lo_cl * (CRP - 8) * (CRP <= 8) +
                  e_crp_hi_cl * (CRP - 8) * (CRP >  8)
    # BSA linear-deviation effects (reference 1.76 m^2).
    bsa_cl <- 1 + e_bsa_cl * (BSA - 1.76)
    bsa_vc <- 1 + e_bsa_vc * (BSA - 1.76)
    # Race fractional-change multiplier on V/F (at most one of the three
    # indicators is 1 per subject; all zero -> VRACE = 1 for the reference).
    race_vc <- 1 + e_japanese_vc * RACE_JAPANESE +
                   e_black_vc    * RACE_BLACK    +
                   e_asian_vc    * RACE_ASIAN
    # ADA fractional-change multiplier on CL/F (time-varying ADA_POS).
    ada_cl <- 1 + e_ada_cl * ADA_POS

    # 2. Individual PK parameters.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * ada_cl * alb_cl * crp_cl * bsa_cl
    vc <- exp(lvc + etalvc) * race_vc * bsa_vc
    baseline <- exp(lbaseline + etalbaseline)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. ODE system (SC dose into depot; apparent parameters).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observation: predicted concentration = central/V + baseline.
    # Dose in mg / V in L -> mg/L = ug/mL.
    Cc <- central / vc + baseline
    Cc ~ prop(propSd)
  })
}
