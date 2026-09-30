deVelde_2020_imipenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 26 critically",
    "ill adults treated with imipenem-cilastatin in a Geneva intensive care",
    "unit (de Velde 2020, parametric NONMEM arm), parameterised in",
    "micro-constants rather than clearances: an elimination rate constant,",
    "two distribution rate constants and a central volume. The elimination",
    "rate constant scales as a power of absolute (BSA-unadjusted) CKD-EPI",
    "eGFR in mL/min, entered as a time-varying covariate. Between-subject",
    "variability is exponential on the elimination rate constant only, and",
    "residual error is additive on log-transformed concentrations.",
    "The same data were also fitted non-parametrically in Pmetrics; that arm",
    "is the sibling model deVelde_2020_imipenem_nonparametric.",
    sep = " "
  )
  reference <- paste(
    "de Velde F, de Winter BCM, Neely MN, Yamada WM, Koch BCP,",
    "Harbarth S, von Dach E, van Gelder T, Huttner A, Mouton JW, on behalf",
    "of COMBACTE-NET consortium.",
    "Population pharmacokinetics of imipenem in critically ill patients: a",
    "parametric and nonparametric model converge on CKD-EPI estimated",
    "glomerular filtration rate as an impactful covariate.",
    "Clin Pharmacokinet. 2020;59(7):885-898. doi:10.1007/s40262-020-00859-1",
    sep = " "
  )
  vignette <- "deVelde_2020_imipenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Glomerular filtration rate estimated by the CKD-EPI equation and",
        "then DE-NORMALISED to an absolute per-patient rate in mL/min by",
        "multiplying the BSA-normalised value by the patient's body surface",
        "area (de Velde 2020 Sect. 2.5 and Table 1 footnote: 'CKD-EPI-abs",
        "absolute CKD-EPI (i.e. CKD-EPI multiplied by BSA)')."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "UNIT HAZARD. This column is an ABSOLUTE rate in mL/min, NOT the",
        "mL/min/1.73 m^2 that a CKD-EPI calculator returns by default and",
        "that the CRCL canonical otherwise carries. The reference of 119",
        "mL/min is the cohort median CKD-EPI-abs at inclusion (Table 1:",
        "119, IQR 110-139); the cohort median BSA-normalised CKD-EPI was 116",
        "mL/min/1.73 m^2 at a median BSA of 1.89 m^2. Enters the elimination",
        "rate constant as (CRCL/119)^0.655 (Eq. 9). TIME-VARYING in the",
        "source analysis: a median of three creatinine samples per patient",
        "(Table 1) gave per-sample eGFR values, carried NOCB in NONMEM",
        "(Sect. 2.5). Observed CKD-EPI-abs range in the VPC 18-190 mL/min",
        "(Sect. 3.3). Four renal-function estimators (Cockcroft-Gault,",
        "MDRD, CKD-EPI, Jelliffe; the last three also in absolute form)",
        "were screened head to head and CKD-EPI-abs gave the lowest OFV",
        "(Sect. 3.3). A cap on Ke above an eGFR of 150, 120 or 90 was",
        "tested and not retained."
      ),
      source_name = "CKD-EPI-abs"
    )
  )

  # Screened during covariate model building and not retained (Sect. 2.5 and
  # 3.3: 'None of the tested measures of body weight improved the model as a
  # covariate on Ke (dOFV < 3.84) during the univariate analysis').
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as TBW on Ke (power, exponent fixed -0.25 or estimated), not retained (Sect. 3.3). Cohort median 75 kg, IQR 66-85 (Table 1)."
    ),
    IBW = list(
      description = "Ideal body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on Ke, not retained (Sect. 3.3). Cohort median 70 kg, IQR 59-73 (Table 1)."
    ),
    LBW = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on Ke, not retained (Sect. 3.3). Cohort median 58 kg, IQR 46-64 (Table 1). Not a registered canonical in inst/references/covariate-columns.md; recorded here as documentation only, since it is never referenced in model()."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 1L,
    age_median = "51 years (IQR 39-54); inclusion 18-60 years",
    weight_median = "75 kg (IQR 66-85)",
    sex_female_pct = 30.8,
    race_ethnicity = NULL,
    disease_state = paste(
      "Critically ill adults with suspected or documented severe bacterial",
      "infection (lower respiratory tract 62%, intra-abdominal 15%,",
      "bloodstream 12%, other 12%); APACHE II median 22 (IQR 17-27); no",
      "continuous renal replacement therapy. Exclusion: Cockcroft-Gault",
      "eGFR < 60 mL/min, BMI < 18 or > 30 kg/m^2, pregnancy."
    ),
    dose_range = paste(
      "Imipenem/cilastatin 500 mg/500 mg four times daily as a 30-min",
      "intermittent IV infusion (Sect. 2.1)."
    ),
    regions = "Switzerland (Geneva University Hospitals ICU)",
    renal_function = "CKD-EPI 116 mL/min/1.73 m^2 (IQR 104-124); CKD-EPI-abs 119 mL/min (IQR 110-139); Cockcroft-Gault 146 mL/min (IQR 123-170) at inclusion (Table 1).",
    n_concentrations = 125L,
    notes = paste(
      "Data from the last 27 of 54 imipenem-treated patients of a",
      "prospective ICU cohort (2010-2013) whose exact dosing and sampling",
      "times were known; one excluded for missing height (Sect. 2.1).",
      "138 samples (peak, intermediate, trough on days 1, 2, 3, 4 and 6);",
      "13 (9.4%) below the 0.5 mg/L LOQ were excluded, leaving 125",
      "concentrations after 84 doses (Sect. 3.2). HPLC-UV assay, linear",
      "0.5-80 mg/L. FOCE-I in NONMEM 7.2 on log-transformed data.",
      "Covariates: forward inclusion p < 0.05 (dOFV 3.84), backward",
      "elimination p < 0.001 (dOFV 10.83), tested only on Ke (the only",
      "parameter with BSV). Eta shrinkage 14%."
    )
  )

  ini({
    # ===== Structural PK -- de Velde 2020 Table 2, NONMEM 'Final model'
    # column. The model is parameterised in MICRO-CONSTANTS (Ke, Kcp, Kpc)
    # plus the central volume; clearance is a derived quantity
    # (CL = Vc x Ke, Discussion). BSV on Vc, Kcp and Kpc 'did not
    # significantly improve the model' (Sect. 3.3), so those carry no eta. =====
    lkel <- log(0.637); label("Elimination rate constant at CRCL = 119 mL/min (1/h)")  # Table 2 NONMEM: Ke = 0.637 1/h (bootstrap median 0.634, 95% CI 0.543-0.805)
    lk12 <- log(0.166); label("Central-to-peripheral rate constant Kcp (1/h)")         # Table 2 NONMEM: Kcp = 0.166 1/h (bootstrap 95% CI 0.092-0.436)
    lk21 <- log(0.195); label("Peripheral-to-central rate constant Kpc (1/h)")         # Table 2 NONMEM: Kpc = 0.195 1/h (bootstrap 95% CI 0.079-0.604)
    lvc  <- log(29.6);  label("Central volume of distribution (L)")                    # Table 2 NONMEM: Vc = 29.6 L (bootstrap 95% CI 22.9-34.4)

    # ===== Covariate effect -- Eq. 9:
    #   Ke_i = 0.637 x (CKD-EPI-abs_i / 119)^0.655 x e^eta
    # The covariate acts on the rate constant, not on clearance. =====
    e_crcl_kel <- 0.655; label("Power exponent on (CRCL/119) for kel (unitless)")  # Table 2 NONMEM: Ke(cov) = 0.655 (bootstrap 95% CI 0.474-1.184); Eq. 9

    # ===== Between-subject variability -- Sect. 3.3 text below Eq. 8:
    # 'variance omega^2 (estimated from the data as 0.0354)'. Table 2 prints
    # the equivalent CV 19.0% from Eq. 2, CV = sqrt(exp(omega^2) - 1). =====
    etalkel ~ 0.0354  # Sect. 3.3: omega^2 = 0.0354 (Table 2 CV 19.0%)

    # ===== Residual error -- Eq. 1 on log-transformed data:
    #   log(OBS) = log(IPRED) + sqrt(error^2) x eps, eps fixed to 1 (variance)
    # so the estimated 'Exponential error' is the SD on the log scale. =====
    expSd <- 0.348; label("Residual SD on the log scale (unitless)")  # Table 2 NONMEM: Exponential error = 0.348 (bootstrap 95% CI 0.281-0.413)
  })

  model({
    # ----- Individual PK parameters (Eq. 9) -----
    kel <- exp(lkel + etalkel) * (CRCL / 119)^e_crcl_kel
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    vc  <- exp(lvc)

    # ----- ODE system -----
    # Imipenem-cilastatin given as a 30-min IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
