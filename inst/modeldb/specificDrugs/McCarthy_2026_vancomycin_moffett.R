McCarthy_2026_vancomycin_moffett <- function() {
  description <- paste(
    "One-compartment IV population PK model for vancomycin in paediatric cardiovascular",
    "intensive care unit (CVICU) patients: the Moffett 2019 paediatric cardiac-surgery structure",
    "re-estimated by McCarthy 2026 on 219 paediatric CVICU patients. Clearance scales allometrically",
    "with (WT/70)^0.75, with a power of bedside-Schwartz creatinine clearance (CRCL/84) computed",
    "inside the model from height and serum creatinine, and with a postmenstrual-age maturation",
    "1 / (1 + (PMA/50)^theta); volume scales linearly with WT/70. CL, V, the maturation and CRCL",
    "exponents, IIV and the proportional residual error were re-estimated.",
    sep = " "
  )
  reference <- paste(
    "McCarthy MG, Keizer RJ, Hughes JH. Evaluation and Improvement of Specialized Vancomycin",
    "Pharmacokinetic Models for Pediatric Cardiovascular Intensive Care Unit and Pediatric Oncology",
    "Patients. Clin Pharmacokinet. 2026;65(4):527. doi:10.1007/s40262-026-01623-7.",
    "(Table 5 'Moffett 2019 (refit)'; Supplementary file 2 'Moffett 2019 (refit)' NONMEM control",
    "stream; Table S2 structure.) Structural model from Moffett BS, Resendiz K, Morris J,",
    "Akcan-Arikan A, Checchia PA. Population pharmacokinetics of vancomycin in the pediatric",
    "cardiac surgical population. J Pediatr Pharmacol Ther. 2019;24:107-116.",
    "doi:10.5863/1551-6776-24.2.107.",
    sep = " "
  )
  vignette <- "McCarthy_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "(WT/70)^0.75 on CL and WT/70 on V in the control stream. McCarthy 2026 Table 1,",
        "Study 2 CVICU: median 9 kg (range 3-91).",
        sep = " "
      ),
      source_name = "WT"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only in the creatinine clearance computed inside the model,",
        "CRCL = 0.413 * HT / CR (control stream, commented 'Schwartz equation'; mL/min/1.73 m^2",
        "with HT in cm and CR in mg/dL). Table 1, Study 2 CVICU: median 73 cm (range 51-187).",
        sep = " "
      ),
      source_name = "HT"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Denominator of the in-model creatinine clearance CRCL = 0.413 * HT / CR. Table 1,",
        "Study 2 CVICU: median 0.31 mg/dL (range 0.06-2.44).",
        sep = " "
      ),
      source_name = "CR"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Maturation FMAT = 1 / (1 + (PMA/50)^THETA(3)) with PMA in weeks (Table 5 footnote",
        "'PMA post-menstrual age (wk)'). Where gestational age was unknown McCarthy 2026 imputed",
        "PMA as chronological age plus 40 weeks. Table 1, Study 2 CVICU: median 99 weeks",
        "(range 46-970).",
        sep = " "
      ),
      source_name = "PMA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 219L,
    n_studies = 1L,
    age_range = "0.1-17.8 years",
    age_median = "1.1 years",
    weight_range = "3-91 kg",
    weight_median = "9 kg",
    sex_female_pct = 46.1,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Paediatric patients (< 20 years) tagged as cardiovascular intensive care unit (CVICU)",
      "patients and treated with intermittent intravenous vancomycin. Patients on haemodialysis,",
      "CRRT, ECMO or a ventricular assist device, or also tagged NICU or PICU, were excluded.",
      sep = " "
    ),
    dose_range = "Routine-care intermittent IV vancomycin; median dose 125 mg (range 25-2000), median interval 7.9 h, median infusion 1 h",
    regions = "United States (11 sites)",
    renal_function = "Serum creatinine median 0.31 mg/dL (range 0.06-2.44)",
    pma_range = "46-970 weeks",
    n_concentrations = 1136L,
    notes = paste(
      "McCarthy 2026 Table 1, Study 2 (CVICU only; the matched non-CVICU controls were not used",
      "in the refit). 391 treatment courses, each treated as a separate individual in NONMEM",
      "(Methods, Pharmacokinetic Analysis). 101/219 female. Routine therapeutic drug monitoring",
      "records from the InsightRX Nova dosing software, January 2018 to April 2025.",
      sep = " "
    )
  )

  ini({
    lcl <- log(7.4986); label("Clearance at WT 70 kg, CRCL 84 mL/min/1.73 m^2, before maturation (L/h)")
    # Table 5 refit thetaCL 7.50 (SE 0.03); control stream $THETA(1) 7.4986
    lvc <- log(70.4378); label("Volume of distribution at WT 70 kg (L)")
    # Table 5 refit thetaVD 70.44 (SE 0.08); control stream $THETA(2) 70.4378

    e_wt_cl <- fixed(0.75); label("Allometric exponent on WT/70 for CL (unitless)")
    # control stream TVCL '(WT/70)**0.75'; Table S2
    e_page_cl <- -0.169922; label("Exponent on PAGE/50 weeks in the maturation 1/(1+(PAGE/50)^e) for CL (unitless)")
    # Table 5 refit thetaPMA -0.17 (SE 0.37); control stream $THETA(3) -0.169922
    tmat50 <- fixed(50); label("Postmenstrual age at 50% clearance maturation (weeks)")
    # control stream FMAT '(PMA/50)'; Table S2 'PMA/50'
    e_crcl_cl <- 0.768862; label("Exponent on CRCL/84 for CL (unitless)")
    # Table 5 refit thetaCRCL 0.77 (SE 0.08); control stream $THETA(4) 0.768862

    # IIV: $OMEGA variances; Table 5 prints their square roots (omega).
    # sqrt(0.103392) = 0.322 -> Table 5 refit omegaCL 0.32
    etalcl ~ 0.103392
    # sqrt(0.123482) = 0.351 -> Table 5 refit omegaVD 0.35
    etalvc ~ 0.123482

    # Residual error: W = SQRT(IPRED**2 * PROP**2 + ADD**2), Y = IPRED + W*ERR(1),
    # $SIGMA 1 FIX, so PROP and ADD are standard deviations.
    propSd <- 0.238867; label("Proportional residual error (fraction)")
    # Table 5 refit thetaPROP 0.24 (SE 0.03); control stream $THETA(5) 0.238867
    addSd <- fixed(0.0001); label("Additive residual error (mg/L)")
    # control stream $THETA(6) 0.0001 FIX (a numerical floor; not in Table 5)
  })
  model({
    # Creatinine clearance (mL/min/1.73 m^2) computed in the control stream as
    # CRCL = 0.413 * HT/CR.
    crcl <- 0.413 * HT / CREAT

    fmat <- 1 / (1 + (PAGE / tmat50)^e_page_cl)

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (crcl / 84)^e_crcl_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)

    kel <- cl / vc

    # The control stream's second state, DADT(2) = A(1)/V1, is an AUC
    # accumulator with no effect on the PK; it is omitted here.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
