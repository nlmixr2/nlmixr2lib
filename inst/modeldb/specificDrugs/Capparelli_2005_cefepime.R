Capparelli_2005_cefepime <- function() {
  description <- "One-compartment population PK model for cefepime in preterm and term neonates (Capparelli 2005); additive renal-plus-non-renal CL on serum creatinine, additive Vc step for PCA < 30 weeks."
  reference <- "Capparelli E, Hochwald C, Rasmussen M, Parham A, Bradley J, Moya F. Population pharmacokinetics of cefepime in the neonate. Antimicrob Agents Chemother. 2005;49(7):2760-2766. doi:10.1128/aac.49.7.2760-2766.2005"
  vignette <- "Capparelli_2005_cefepime"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-kg parameterization in Capparelli 2005 Methods: 'parameters were scaled by subject weight before evaluation of other potential covariates'. CL and Vc are computed per kg in ini() and multiplied by WT in model().",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the renal arm of CL via cl_renal = exp(lcl_renal)/CREAT (mL/min/kg). Capparelli 2005 Table 1 cohort range 0.3-1.5 mg/dL, mean 0.8. Source column name 'SCr' (enzymatic at San Diego, modified Jaffe at Houston; both bilirubin-corrected).",
      source_name        = "SCr"
    ),
    PAGE = list(
      description        = "Postmenstrual (postconceptional) age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Capparelli 2005 defines postconceptional age (PCA) as GA at birth (weeks) + postnatal age (weeks); this matches the modern postmenstrual age (PMA) definition. Canonical PAGE is in months; the PCA < 30 weeks threshold from Table 3 maps to PAGE < 30/4.345 months ~= 6.904 months. The covariate enters model() as the derived binary indicator pca30 = (PAGE < 30/4.345).",
      source_name        = "PCA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 2L,
    age_range      = "postnatal 1-62 days",
    age_median     = "postnatal age 14.5 days (mean +/- SD 14.7 +/- 14.5)",
    weight_range   = "0.58-4.70 kg",
    weight_median  = "1.91 kg (mean +/- SD)",
    sex_female_pct = 48,
    race_ethnicity = "Not reported in source paper",
    disease_state  = "Premature and term neonates in NICUs (suspected or documented infection; some subjects also studied with a single study dose)",
    dose_range     = "50 mg/kg IV (30-min infusion); subset received steady-state Q12H dosing",
    regions        = "United States (San Diego, CA and Houston, TX)",
    notes          = "Capparelli 2005 Table 1 and Table 2. Gestational age at birth 30.5 +/- 5.3 weeks (range 22.1-42.3); 42 preterm (<36 weeks) plus 12 term; serum creatinine 0.8 +/- 0.3 mg/dL (range 0.3-1.5). 55 enrolled, 54 used in the final fit (one outlier with massive transmural fluid loss excluded). 35 of 54 from San Diego, 19 from Houston.",
    ga_range       = "22.1-42.3 weeks at birth"
  )

  ini({
    # Structural CL [mL/min/kg] = 0.26 + 0.59/SCr (Capparelli 2005 Table 3).
    # Encoded as the canonical multi-component CL pair:
    #   cl_nonren = exp(lcl_nonren)             (mL/min/kg)
    #   cl_renal  = exp(lcl_renal) / CREAT      (mL/min/kg, CREAT in mg/dL)
    lcl_nonren <- log(0.26); label("Non-renal cefepime clearance per kg (mL/min/kg)")        # Table 3
    lcl_renal  <- log(0.59); label("Renal CL coefficient: cl_renal [mL/min/kg] = exp(lcl_renal)/CREAT[mg/dL]")  # Table 3

    # Structural Vc [L/kg] = 0.385 + 0.122 * PCA30 (Capparelli 2005 Table 3).
    # PCA30 = 1 if postconceptional age < 30 weeks, else 0.
    lvc        <- log(0.385); label("Central volume of distribution per kg for PCA >= 30 weeks (L/kg)")       # Table 3
    e_pca30_vc <- 0.122;      label("Additive shift on Vc for PCA < 30 weeks (L/kg; binary indicator)")        # Table 3

    # Inter-individual variability: Capparelli 2005 used an "exponential normal distribution
    # with the full variance-covariance matrix" but Table 3 reports only the diagonal CV%.
    # omega^2 = log(CV^2 + 1) per the IIV convention; off-diagonal correlation is not in the
    # paper and is encoded as zero here (see vignette Assumptions and deviations).
    # CL: CV = 25%  -> omega^2 = log(1 + 0.25^2) = 0.06062
    # Vc: CV = 29%  -> omega^2 = log(1 + 0.29^2) = 0.08075
    etalcl ~ 0.06062  # Table 3 (CV 25%)
    etalvc ~ 0.08075  # Table 3 (CV 29%)

    # Residual error: proportional 13% (Capparelli 2005 Table 3).
    propSd <- 0.13; label("Proportional residual error (fraction)")  # Table 3
  })
  model({
    # Derived binary indicator for postconceptional age < 30 weeks.
    # PAGE is the canonical postmenstrual age in months; 30 weeks = 30 / 4.345 months ~= 6.904.
    pca30 <- ifelse(PAGE < 30 / 4.345, 1, 0)

    # Per-kg typical PK parameters (Capparelli 2005 Table 3).
    cl_nonren <- exp(lcl_nonren)              # mL/min/kg
    cl_renal  <- exp(lcl_renal) / CREAT       # mL/min/kg (CREAT in mg/dL)
    cl_per_kg <- cl_nonren + cl_renal         # mL/min/kg
    vc_per_kg <- exp(lvc) + e_pca30_vc * pca30 # L/kg

    # Individual PK: exponential IIV on the typical-value per-kg parameters; multiply by WT
    # and convert clearance from mL/min/kg to L/h (factor 60/1000 = 0.06).
    cl <- cl_per_kg * WT * 60 / 1000 * exp(etalcl)  # L/h
    vc <- vc_per_kg * WT * exp(etalvc)              # L

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Concentration: dose mg / volume L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
