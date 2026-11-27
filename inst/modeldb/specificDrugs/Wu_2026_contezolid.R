Wu_2026_contezolid <- function() {
  description <- paste(
    "Population PK model for oral contezolid in plasma and pulmonary",
    "epithelial lining fluid (ELF) of 15 Chinese adults with pneumonia",
    "receiving 800 mg every 12 h (Wu 2026). First-order absorption with a",
    "lag time from a depot into a central compartment that exchanges",
    "drug with an ELF compartment through a directional pair of",
    "first-order rate constants; linear elimination from central.",
    "Bioavailability is 0.557 fasted and 1.83-fold higher fed, with",
    "inter-occasion variability; between-subject variability on the",
    "absorption rate constant and on the central-to-ELF rate constant.",
    "Proportional residual error differs between the two recruiting",
    "hospitals for plasma. The published albumin effect on the",
    "absorption rate constant is not encoded because the albumin scale",
    "is not reported; the absorption rate constant is its value at the",
    "albumin reference."
  )
  reference <- paste(
    "Wu H, Cheng Z, Ding F, Hua J, Wang H, Meng X, Cao Y, Chen Y, Chen Y,",
    "Guo B, Chen T, Zhang J, Zhang J. Evaluation of the clinical efficacy,",
    "safety, and permeability in pulmonary epithelial lining fluid of",
    "contezolid, a novel oxazolidinone drug, in adult patients with",
    "pneumonia. Antimicrob Agents Chemother. 2026;70(4):e01033-25.",
    "doi:10.1128/aac.01033-25"
  )
  vignette <- "Wu_2026_contezolid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Wu 2026 Table S1 estimates one proportional plasma residual variance per
  # recruiting hospital; the canonical single propSd cannot hold both.
  paper_specific_residual_sds <- c("propSd_zhongshan", "propSd_firsthosp")

  compartmentData <- list(
    depot = list(analyte = "contezolid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "contezolid", units = "mg", specimen = "plasma", verified = TRUE),
    elf = list(analyte = "contezolid", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  covariateData <- list(
    FED = list(
      description = "1 = dose taken after a meal, 0 = dose taken fasted",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Multiplies bioavailability by 1.83^FED (supplement Results",
        "equation 'BIO = 0.557 x 1.83^FASTING x exp(eta3 + kappa)';",
        "Table S1 row 'BIO(FASTING)' = 1.83). The source column is named",
        "FASTING but the supplement codes it 'value=0: fast, value=1:",
        "food', so FED = FASTING with no inversion. In the cohort 13 of 15",
        "patients took the sampled doses fasted (required before",
        "bronchoscopy) and 2 took them after a meal. Fed bioavailability",
        "is 0.557 x 1.83 = 1.02."
      ),
      source_name = "FASTING"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for inter-occasion variability on bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Two occasions, decomposed inside model() into oc1 / oc2 that",
        "multiplex the IOV etas on log-bioavailability (the supplement's",
        "kappa term, Table S1 'pi2 BIO' = 0.066). The paper does not say",
        "how an occasion was delimited. Blood and BAL sampling both",
        "followed 'the fifth or sixth dose', so OCC = 1 is read as the",
        "fifth-dose interval and OCC = 2 as the sixth-dose interval.",
        "Records with any other OCC value (e.g. 0) carry no IOV. For a",
        "single-occasion simulation pass OCC = 1."
      ),
      source_name = "OCC"
    ),
    STUDY_SHANGHAI_FIRST = list(
      description = "1 = patient recruited at Shanghai First People's Hospital, 0 = Zhongshan Hospital (Fudan University)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the plasma proportional residual SD (supplement",
        "Methods: 'site is an indicator variable (0 for Zhongshan",
        "Hospital, 1 for First Hospital)'). No structural effect. ELF",
        "residual error is estimated for Zhongshan only (Table S1), so",
        "the ELF residual does not depend on this column."
      ),
      source_name = "site"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Baseline serum albumin",
      units = "unknown scale (see notes)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Retained in the published final model as a power effect on the",
        "absorption rate constant, Ka = 0.48 x (ALB / 0.93)^8.02 x",
        "exp(eta1) (supplement Results equation; Table S1 'Ka(ALB)' =",
        "8.02, bootstrap 6.89, RSE 39.3%). The reference 0.93 is not",
        "albumin in g/L or g/dL, and albumin is not summarised anywhere",
        "in the paper, so the scale of ALB in the authors' dataset cannot",
        "be recovered. The same group's nemonoxacin analysis",
        "(Chen_2023_nemonoxacin.R) reported albumin as measured value /",
        "midpoint of the normal range, on which 0.93 would be a plausible",
        "median, but Wu 2026 does not say so. The effect is therefore not",
        "encoded; the model's Ka is the value at ALB = 0.93. Ka only",
        "shapes absorption: steady-state AUC, the ELF penetration ratio",
        "and the food effect are unaffected."
      ),
      source_name = "ALB"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 1L,
    n_observations = "75 plasma concentrations (15 patients) and 11 ELF concentrations (11 patients)",
    age_range = "36 to 73 years",
    age_median = "58 years (mean 55)",
    weight_range = "47 to 105 kg",
    weight_median = "65 kg (mean 69)",
    height_range = "1.52 to 1.84 m",
    bmi_range = "17.3 to 29.3 kg/m^2 (median 25.4)",
    sex_female_pct = 100 * 6 / 15,
    race_ethnicity = c(Asian = 100),
    renal_function = "eGFR 79 to 117 mL/min (median 96.0, mean 97.3)",
    disease_state = paste(
      "Hospitalised adults with clinically confirmed or suspected",
      "gram-positive bacterial pneumonia: 14 community-acquired, 1",
      "hospital-acquired. Mechanically ventilated, vasopressor-dependent",
      "and PSI class V patients were excluded."
    ),
    dose_range = paste(
      "Contezolid tablets 800 mg orally every 12 h for 7 to 14 days.",
      "Serial plasma samples after the fifth or sixth dose (group 1: 1, 2,",
      "4, 6 h; group 2: 2, 4, 6, 12 h); bronchoalveolar lavage at 2, 4 or",
      "6 h after the fifth or sixth dose in 13 fasted patients, 11 with",
      "evaluable ELF (urea-corrected)."
    ),
    regions = "China (Zhongshan Hospital, Fudan University, and Shanghai First People's Hospital, Shanghai).",
    notes = paste(
      "Prospective single-arm open-label study. Demographics from Table",
      "1; structural model, covariate equations and ODEs from the",
      "supplement (Fig. S1 and the Results equation); estimates from",
      "Table S1. NONMEM 7.5.6, FOCE with interaction."
    )
  )

  ini({
    # Structural parameters: Table S1 'Estimate based on original dataset'
    # and the supplement Results equation.
    lka <- log(0.48); label("Absorption rate constant at the albumin reference (1/h)") # Table S1 'Ka (1/h)' = 0.48 (bootstrap 0.54, RSE 53.8%); equation Ka = 0.48 x (ALB/0.93)^8.02 x exp(eta1)
    ltlag <- log(0.083); label("Absorption lag time (h)") # Table S1 'TLAG (h)' = 0.083 (bootstrap 0.081, RSE 145.0%)
    lcl <- log(9.42); label("Clearance before bioavailability correction (L/h)") # Table S1 'CL' = 9.42 (bootstrap 9.43, RSE 12.6%); equation 'CL (L/h) = 9.42'
    lvc <- log(13.11); label("Central volume before bioavailability correction (L)") # Table S1 'Vc (L)' = 13.11 (bootstrap 12.53, RSE 27.0%)
    lvelf <- log(12.82); label("ELF volume before bioavailability correction (L)") # Table S1 'VELF (L)' = 12.82 (bootstrap 21.21, RSE 32.5%)
    lk_central_elf <- log(0.99); label("Central-to-ELF transfer rate constant (1/h)") # Table S1 'KCELF (1/h)' = 0.99 (bootstrap 1.00, RSE 2.5%)
    lk_elf_central <- log(2.99); label("ELF-to-central transfer rate constant (1/h)") # Table S1 'KELFC (1/h)' = 2.99 (bootstrap 2.98, RSE 4.7%)
    lfdepot <- log(0.557); label("Oral bioavailability, fasted (fraction)") # supplement Results equation 'BIO = 0.557 x 1.83^FASTING x exp(eta3 + kappa)'; no Table S1 row
    e_fed_fdepot <- 1.83; label("Multiplicative food effect on bioavailability, base of the FED power (unitless)") # Table S1 'BIO(FASTING)' = 1.83 (bootstrap 1.90, RSE 19.6%)

    # Between-subject variability (variances; Table S1 omega2 rows).
    etalka ~ 0.839 # Table S1 'omega2 Ka' = 0.839 (shrinkage 14.8%)
    etalk_central_elf ~ 0.47 # Table S1 'omega2 KCELF' = 0.47 (shrinkage 13.8%)
    # The BIO equation carries eta3 but Table S1 reports no omega2 for it;
    # held at zero variance so only the kappa (IOV) term varies BIO.
    etalfdepot ~ fixed(0)

    # Inter-occasion variability on bioavailability (kappa), two occasions;
    # occasion 2 repeats the occasion-1 variance (NONMEM BLOCK SAME).
    etaiov_fdepot_1 ~ 0.066 # Table S1 'pi2 BIO' = 0.066 (shrinkage 46.0%)
    etaiov_fdepot_2 ~ fixed(0.066)

    # Residual error, proportional; Table S1 reports variances.
    propSd_zhongshan <- 0.5099; label("Proportional residual SD, plasma, Zhongshan Hospital (fraction)") # Table S1 'sigma2 Zhongshan_Plasma' = 0.26; sqrt(0.26) = 0.5099
    propSd_firsthosp <- 0.3606; label("Proportional residual SD, plasma, Shanghai First People's Hospital (fraction)") # Table S1 'sigma2 First_hospital_Plasma' = 0.13; sqrt(0.13) = 0.3606
    propSd_Celf <- 0.001; label("Proportional residual SD, ELF (fraction)") # Table S1 'sigma2 Zhongshan_ELF' = 1.00E-06 (shrinkage 99.7%); sqrt(1e-6) = 0.001
  })
  model({
    # Inter-occasion indicators (OCC outside 1-2 gives no IOV).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp(lcl)
    vc <- exp(lvc)
    velf <- exp(lvelf)
    k_central_elf <- exp(lk_central_elf + etalk_central_elf)
    k_elf_central <- exp(lk_elf_central)

    # BIO = 0.557 x 1.83^FASTING x exp(eta3 + kappa).
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) * e_fed_fdepot^FED

    # The supplement divides CL, Vc and VELF by BIO and also sets the
    # initial depot amount to Dose x BIO. Bioavailability is applied once
    # here, through f(depot): CL/BIO over Vc/BIO leaves the elimination
    # rate constant cl / vc unchanged, and the concentrations equal
    # depot-input x BIO over the uncorrected volumes.
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k_central_elf * central + k_elf_central * elf
    d/dt(elf) <- k_central_elf * central - k_elf_central * elf

    f(depot) <- fdepot
    alag(depot) <- tlag

    Cc <- central / vc
    Celf <- elf / velf

    propSd_i <- propSd_zhongshan * (1 - STUDY_SHANGHAI_FIRST) + propSd_firsthosp * STUDY_SHANGHAI_FIRST
    Cc ~ prop(propSd_i)
    Celf ~ prop(propSd_Celf)
  })
}
