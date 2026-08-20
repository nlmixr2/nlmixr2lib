Li_2021_ganciclovir <- function() {
  description <- paste(
    "One-compartment population PK model for IV ganciclovir in critically ill",
    "children (Li 2021), with a kidney-function power effect and near-linear",
    "body-weight scaling on clearance.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Li 2021 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Li S, Shu C, Wu S, Xu H, Wang Y. Population Pharmacokinetics and Dose",
    "Optimization of Ganciclovir in Critically Ill Children.",
    "Front Pharmacol. 2020;11:614164. doi:10.3389/fphar.2020.614164.",
    "Reported with a 2021 publication year in the Yang 2023 repository (Table 3).",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central = list(analyte = "ganciclovir", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling referenced to the cohort median 12.0 kg. Both exponents are",
        "estimated rather than fixed at canonical allometric values (1.02 on CL,",
        "0.80 on Vc). Cohort median 12.0 kg, range 2.5-55.0 kg (Yang 2023 Table 2).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (Gao formula, pediatric)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the model through the derived kidney-function ratio",
        "KF = eGFR / 120, i.e. CRCL normalized to 120 mL/min/1.73 m^2",
        "(Yang 2023 Table 3 footnote: 'KF: kidney function, KF = eGFR/(120",
        "mL/min/1.73 m^2)'). The Yang 2023 Discussion notes that this study used",
        "the Gao formula to estimate eGFR, which is best suited to children with",
        "moderate renal failure -- one of the reasons the Yang 2023 authors give",
        "for this model's simulated profiles differing from the others.",
        sep = " "
      ),
      source_name        = "eGFR (via KF)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 104L,
    n_studies      = 1L,
    n_observations = 138L,
    age_median     = "2.46 years (range 0.10-12.83); mean 3.06 (SD 2.99) years",
    weight_median  = "12.0 kg (range 2.5-55.0); mean 13.7 (SD 8.3) kg",
    sex_female_pct = 48.1,
    race_ethnicity = "Chinese (single-centre Chinese cohort; not further specified).",
    disease_state  = "Critically ill pediatric patients.",
    dose_range     = "IV ganciclovir 5 mg/kg/12 h.",
    regions        = "China (retrospective study).",
    bioassay       = "HPLC, LLOQ 0.1 ug/mL.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Sampling strategy not",
      "reported; the dataset is very sparse (138 samples from 104 subjects), which",
      "the Yang 2023 authors cite as a likely limitation of this model's",
      "performance. Only 11.5% of patients were older than six years.",
      "Covariates tested: weight, sex, age, height, BSA, BUN, SCR, uric acid,",
      "total bilirubin, ALT, AST and kidney function; weight and KF were retained",
      "(Yang 2023 Table 4). Critical illness characterises the whole cohort and is",
      "not a covariate in this model.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Li et al. (2021) row. Reference subject:
    # WT = 12.0 kg, eGFR = 120 mL/min/1.73 m^2 (i.e. KF = 1). Clearance in L/h,
    # volume in L. One-compartment IV model; no absorption parameters.
    lcl <- log(5.23) ; label("Clearance at WT = 12.0 kg and KF = 1 (CL, L/h)")     # Yang 2023 Table 3 (Li 2021): CL = 5.23 * KF^0.92 * (BW/12.0)^1.02
    lvc <- log(11.35); label("Central volume of distribution at WT = 12.0 kg (Vc, L)") # Yang 2023 Table 3 (Li 2021): Vc = 11.35 * (BW/12.0)^0.80

    # Covariate effects. All three exponents are non-canonical estimated values.
    e_crcl_cl <- 0.92; label("Power exponent of kidney function KF = CRCL/120 on CL (unitless)") # Yang 2023 Table 3 (Li 2021): KF^0.92
    e_wt_cl   <- 1.02; label("Power exponent of WT on CL (unitless; reference 12.0 kg)")         # Yang 2023 Table 3 (Li 2021): (BW/12.0)^1.02
    e_wt_vc   <- 0.80; label("Power exponent of WT on Vc (unitless; reference 12.0 kg)")         # Yang 2023 Table 3 (Li 2021): (BW/12.0)^0.80

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2.
    etalcl ~ 0.016641  # Yang 2023 Table 3 (Li 2021): BSV CL = 12.9% -> 0.129^2
    etalvc ~ 0.432964  # Yang 2023 Table 3 (Li 2021): BSV Vc = 65.8% -> 0.658^2

    # Residual unexplained variability: proportional.
    propSd <- 0.0823; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Li 2021): 8.23% proportional error
  })

  model({
    # Derived kidney-function ratio KF, normalized to 120 mL/min/1.73 m^2.
    kf <- CRCL / 120

    cl <- exp(lcl + etalcl) * kf^e_crcl_cl * (WT / 12.0)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 12.0)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
