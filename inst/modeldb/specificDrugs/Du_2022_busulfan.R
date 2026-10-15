Du_2022_busulfan <- function() {
  description <- "One-compartment IV popPK model for busulfan in Chinese pediatric hematopoietic stem cell transplantation recipients, with normal-fat-mass allometric scaling (estimated fat fraction, exponent fixed at 3/4) and a postmenstrual-age sigmoid maturation function on CL, and linear fat-free-mass scaling on V (Du 2022)."
  reference <- "Du X, Huang C, Xue L, Jiao Z, Zhu M, Li J, Lu J, Xiao P, Zhou X, Mao C, Zhu Z, Dong J, Liu X, Chen Z, Zhang S, Ding Y, Hu S, Miao L. The Correlation Between Busulfan Exposure and Clinical Outcomes in Chinese Pediatric Patients: A Population Pharmacokinetic Study. Front Pharmacol. 2022;13:905879. doi:10.3389/fphar.2022.905879"
  vignette <- "Du_2022_busulfan"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "busulfan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at baseline. Combined with FFM to build the clearance-specific normal fat mass, NFM = FFM + Ffat * (WT - FFM) (Supplementary Equation S4). The standard subject is a 70 kg, 176 cm adult male (Methods, Equations 1-2).",
      source_name = "ABW"
    ),
    FFM = list(
      description = "Fat-free mass (Janmahasatian 2005 semi-mechanistic equation, derived from body weight, height and sex)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed at baseline. Supplementary Equation S3: FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT), HT in metres, with WHSmax = 42.92 and WHS50 = 30.93 kg/m^2 for males and 37.99 and 35.98 kg/m^2 for females. The paper applies this adult equation to children without a paediatric correction. FFMSTD = 56.1 kg is the male equation at 70 kg and 1.76 m (Table 1 and 2 footnotes).",
      source_name = "FFM"
    ),
    AGE = list(
      description = "Postnatal age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Combined with GA inside model() to form postmenstrual age, PMA (weeks) = AGE * 52 + GA, exactly as written in Supplementary Equation S6 (52 weeks per year).",
      source_name = "Age"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Second component of postmenstrual age (Supplementary Equation S6, where it is called PW). Actual gestational ages were recorded for this cohort, mean 39.82 weeks (range 32-40; Supplementary Table S3).",
      source_name = "PW"
    )
  )

  # Screened in the stepwise covariate search (Supplementary Table S5 and Results)
  # but not retained in the final model; documentation only.
  covariatesDataExcluded <- list(
    BUN = list(
      description = "Urea nitrogen",
      units = "mmol/L",
      type = "continuous",
      notes = "Passed forward inclusion on CL (dOFV -23.64) and survived backward elimination, but its positive coefficient (0.244) was judged physiologically implausible and it was removed from the final model (Results; Supplementary Table S5 rows 2 and 26).",
      source_name = "UREA"
    ),
    GGT = list(
      description = "Gamma-glutamyl transpeptidase",
      units = "U/L",
      type = "continuous",
      notes = "Passed forward inclusion on V (dOFV -15.84) and survived backward elimination, but its coefficient (-0.1) was judged physiologically implausible and it was removed (Results; Supplementary Table S5 rows 7 and 23)."
    )
  )
  # Also screened and not retained (Results; Discussion; Supplementary Table
  # S5): malignant vs nonmalignant disease, GSTA1 / GSTP1 genotypes, daily
  # doses of fludarabine, phenytoin and metronidazole, and the remaining
  # hematological and biochemical indicators of Supplementary Table S3.

  population <- list(
    species = "human",
    n_subjects = 128,
    n_studies = 1,
    n_centers = 1,
    n_observations = 467,
    age_range = "0.6-17.0 years (mean 6.11)",
    weight_range = "7.5-96.5 kg (mean 23.99)",
    height_range = "67-185 cm (mean 115.11)",
    ga_range = "32-40 weeks (mean 39.82)",
    pma_range = "70.6-926 weeks (mean 368.32)",
    sex_female_pct = 29.7,
    race_ethnicity = "Chinese",
    disease_state = "Pediatric patients receiving IV busulfan as part of conditioning before hematopoietic stem cell transplantation; 70% malignant (AML, ALL, MDS, MPN, JMML) and 30% nonmalignant disease (Wiskott-Aldrich syndrome, thalassemia, severe aplastic anemia)",
    dose_range = "0.8, 1.0 or 1.2 mg/kg per dose (actual or adjusted body weight), 2 h IV infusion four times daily for 2-4 days (8-16 doses)",
    regions = "China (Children's Hospital of Soochow University, Suzhou)",
    notes = "Patients enrolled July 2018 to February 2021 (Methods). Sampling at 0, 2 and 4 h after the end of the first 2 h infusion, plus a pre-dose sample before the fifth dose in 61 patients. Assay LLOQ 0.1 ug/mL. Demographics from Supplementary Table S3; dose strata from Supplementary Table S2. Between-occasion variability on CL and V was tested and not retained (Supplementary Table S5 rows 20-21)."
  )

  ini({
    # Structural parameters standardised to a 70 kg, 176 cm adult male
    # (Table 2 final model = Table 1 Model III).
    lcl <- log(7.71);  label("Typical clearance for a 70 kg, 176 cm adult male (L/h)")              # Table 2 CL_STD = 7.71 L/h
    lvc <- log(42.4);  label("Typical volume of distribution for a 70 kg, 176 cm adult male (L)")   # Table 2 V_STD = 42.4 L

    # Fraction of fat mass contributing to the clearance NFM (Supplementary
    # Equation S4). The corresponding fraction for V was zero (Results), so V
    # scales on FFM alone and has no parameter here.
    ffat_cl <- 0.692;  label("Fraction of fat mass contributing to NFM on CL (unitless)")           # Table 2 Ffat_CL = 0.692

    # Allometric exponent on CL fixed to the theory-based value in Model III.
    e_nfm_cl <- fixed(0.75);  label("Allometric exponent of NFM on CL (unitless)")                  # Methods; Table 1 Model III k1 = 0.75 (fixed)

    # Sigmoid maturation of CL on postmenstrual age (Supplementary Equation S5).
    tm50_mat <- 31.0;  label("Postmenstrual age at which CL reaches 50% of the adult value (weeks)") # Table 2 TM50 = 31.0 weeks
    hill_mat <- 2.03;  label("Hill coefficient of the CL maturation function (unitless)")           # Table 2 HILL = 2.03

    # Between-subject variability (exponential model, Supplementary Equation S1).
    # Table 2 BSV_CL 0.234 and BSV_V 0.240 (Results: 23.4% and 24.0%) are read
    # as the SD of eta, so omega^2 = 0.234^2 and 0.240^2.
    etalcl ~ 0.054756  # Table 2 BSV_CL 0.234 -> 0.234^2
    etalvc ~ 0.0576    # Table 2 BSV_V 0.240 -> 0.240^2

    # Combined residual error, Supplementary Equation S2:
    # Y = Con + sqrt(Con^2 * thetaPROP^2 + thetaADD^2) * eps, i.e. variances add.
    propSd <- 0.130;  label("Proportional residual error (fraction)")                               # Table 2 RUV_PROP = 0.130
    addSd  <- 0.048;  label("Additive residual error (mg/L)")                                       # Table 2 RUV_ADD = 0.048 mg/L
  })

  model({
    # Normal fat mass for CL (Supplementary Equation S4). The standard NFM is the
    # same formula evaluated for the 70 kg adult male with FFMSTD = 56.1 kg
    # (Methods: 'NFMSTD was calculated using Supplementary Equations S3, S4').
    nfm_cl     <- FFM + ffat_cl * (WT - FFM)
    nfm_std_cl <- 56.1 + ffat_cl * (70 - 56.1)

    # Postmenstrual age (Supplementary Equation S6) and maturation (S5).
    pma  <- AGE * 52 + GA
    fmat <- 1 / (1 + (pma / tm50_mat)^(-hill_mat))

    # Equation 1 (CL) and Equation 2 (V; exponent 1, Ffat_V = 0 so NFM = FFM).
    cl <- exp(lcl + etalcl) * (nfm_cl / nfm_std_cl)^e_nfm_cl * fmat
    vc <- exp(lvc + etalvc) * (FFM / 56.1)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
