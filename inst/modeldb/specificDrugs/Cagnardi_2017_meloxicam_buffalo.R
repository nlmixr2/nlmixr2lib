Cagnardi_2017_meloxicam_buffalo <- function() {
  description <- paste(
    "Preclinical (Mediterranean buffalo calf).",
    "Two-compartment intravenous pharmacokinetic model for meloxicam in",
    "male Mediterranean buffalo calves (Bubalus bubalis) given a single",
    "0.5 mg/kg intravenous bolus immediately after surgical orchiectomy",
    "(Cagnardi 2017). The authors fitted each of the six treated calves",
    "separately in WinNonlin (classical two-compartment model, weights",
    "1/Yhat^2) and published only the mean and SD of the individual",
    "estimates (Table 2). The packaged typical values are the Table 2 means",
    "of the central volume V1, body clearance ClB and peripheral volume V2,",
    "with the intercompartmental clearance taken as mean K12 x mean V1.",
    "Volumes and clearances are per kg body weight, so the dose is given in",
    "ug/kg and central/vc is ug/mL. No between-animal variability or",
    "residual-error magnitude was estimated by the two-stage analysis, so",
    "none is carried and the proportional residual SD is fixed at zero.",
    sep = " "
  )
  reference <- paste(
    "Cagnardi P, Guccione J, Villa R, D'Andrea L, Di Loria A, Ferrante MC,",
    "Borriello G, Zicarelli L, Ciaramella P.",
    "Clinical efficacy and pharmacokinetics of meloxicam in Mediterranean",
    "buffalo calves (Bubalus bubalis).",
    "PLoS One. 2017;12(10):e0187252.",
    "doi:10.1371/journal.pone.0187252.",
    sep = " "
  )
  vignette <- "Cagnardi_2017_meloxicam_buffalo"

  # Table 2 reports every volume in mL/kg and the clearance in mL/h/kg, so the
  # dose is carried per kg in ug/kg (0.5 mg/kg = 500 ug/kg). central/vc is then
  # ug/mL, the unit of the HPLC-UV assay and of Figure 2.
  units <- list(time = "h", dosing = "ug/kg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(
      analyte = "meloxicam",
      units = "ug/kg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "meloxicam",
      units = "ug/kg",
      specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "Mediterranean buffalo (Bubalus bubalis), male calves",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "3-4 months (Methods; the abstract says 4-5 months)",
    age_median = "3.65 months",
    weight_range = "127-135 kg",
    weight_median = "130 kg",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy calves (ASA physical status class I) undergoing elective",
      "surgical orchiectomy under dexmedetomidine sedation and lidocaine",
      "local anaesthesia, reversed with atipamezole",
      sep = " "
    ),
    dose_range = paste(
      "Meloxicam 0.5 mg/kg single intravenous bolus (Metacam 20 mg/mL)",
      "through a jugular catheter immediately after surgery",
      sep = " "
    ),
    regions = "Italy (Caserta district, Campania)",
    co_medication = paste(
      "Dexmedetomidine 5 ug/kg IV, lidocaine 2% into each spermatic cord and",
      "the incision line, atipamezole 0.05 mg/kg IV, oxytetracycline IM and",
      "topical",
      sep = " "
    ),
    notes = paste(
      "Ten male calves were enrolled; six formed the treated group whose",
      "serum meloxicam concentrations were modelled and four were untreated",
      "controls used only for the clinical comparisons. Serum was sampled at",
      "5 and 10 min and at 15, 30 and 45 min and 1, 1.5, 2, 4, 6, 8, 12, 24,",
      "36, 48, 72 and 96 h after surgery (Methods; Supporting Information S1",
      "for the mean profile). Meloxicam was assayed by HPLC-UV (LOQ",
      "0.01 ug/mL).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Disposition. Table 2, 'Mean and S.D. of meloxicam pharmacokinetic
    # parameters after IV administration at 0.5 mg/kg in 6 buffalo calves'.
    # Each row is the arithmetic mean of six individual WinNonlin two-
    # compartment fits, not a population estimate. Means of individual
    # estimates are not mutually consistent (for example mean ClB / mean V1 =
    # 0.108 1/h against mean K10 = 0.13 1/h), so three of them are taken
    # directly and the fourth is derived:
    #   V1, ClB and V2 are the Table 2 means. Together they reproduce the
    #   Table 2 Vdss (V1 + V2 = 142.31 mL/kg, exact), AUC within 4%
    #   (Dose / ClB = 114.2 vs 118.81 h*ug/mL) and MRT within 2%
    #   (Vdss / ClB = 32.5 vs 32.89 h).
    #   Q is mean K12 x mean V1 (4.18 x 40.51). The alternative mean K21 x
    #   mean V2 = 178.2 mL/h/kg changes no concentration at the sampling
    #   times by more than 3%; see the vignette.
    # The K10-based alternative (CL = 0.13 x 40.51 = 5.27 mL/h/kg) would put
    # AUC 20% below the Table 2 mean and is not used.
    # ------------------------------------------------------------------------
    lvc <- fixed(log(40.51))
    label("Central volume of distribution V1 (mL/kg)")
    # Table 2, row 'V1' = 40.51 mL/kg (SD 15.16)

    lcl <- fixed(log(4.38))
    label("Body clearance ClB (mL/h/kg)")
    # Table 2, row 'ClB' = 4.38 mL/h/kg (SD 0.95)

    lvp <- fixed(log(101.8))
    label("Peripheral volume of distribution V2 (mL/kg)")
    # Table 2, row 'V2' = 101.8 mL/kg (SD 67.31)

    lq <- fixed(log(169.33))
    label("Intercompartmental clearance (mL/h/kg)")
    # Derived: Table 2 'K12' 4.18 1/h x 'V1' 40.51 mL/kg = 169.33 mL/h/kg

    # ------------------------------------------------------------------------
    # Residual error. Methods 'Pharmacokinetic analysis': 'All data points
    # were weighted by the inverse square of the fitted value', i.e. a
    # proportional error structure for each individual least-squares fit. No
    # residual variance was reported, so the proportional SD is fixed at 0.
    # ------------------------------------------------------------------------
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; not reported)")
  })

  model({
    vc <- exp(lvc)
    cl <- exp(lcl)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Classical two-compartment model with intravenous bolus into central
    # (Methods 'Pharmacokinetic analysis').
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ug/kg over mL/kg gives ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
