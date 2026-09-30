Cojutti_2020_meropenem <- function() {
  description <- "One-compartment IV population PK model for continuous-infusion meropenem in 61 adult febrile neutropenic patients with hematologic malignancies (Cojutti 2020). Pmetrics NPAG non-parametric fit; clearance is a linear function of CKD-EPI creatinine clearance (CL = theta1 + theta2 * CRCL, both terms with their own between-subject distribution) and volume carries no covariate. Age, height, weight and sex were screened but not retained."
  reference <- "Cojutti PG, Candoni A, Lazzarotto D, Fili C, Zannier M, Fanin R, Pea F. Population Pharmacokinetics of Continuous-Infusion Meropenem in Febrile Neutropenic Patients with Hematologic Malignancies: Dosing Strategies for Optimizing Empirical Treatment against Enterobacterales and P. aeruginosa. Pharmaceutics. 2020;12(9):785. doi:10.3390/pharmaceutics12090785"
  vignette <- "Cojutti_2020_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the CKD-EPI equation, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cojutti 2020 Methods 2.1: 'Serum creatinine was collected at each TDM assessment, and CLCR was",
        "estimated by means of the Chronic Kidney Disease Epidemiology (CKD-EPI) formula'. Time-varying (one",
        "value per TDM assessment). Enters clearance UNCENTRED as the additive linear slope",
        "CL = theta1 + theta2 * CLCR (Results 3.2, Equation 1), so theta1 is the clearance at CLCR = 0.",
        "Table 1: median 107.3 mL/min/1.73 m^2 (IQR 96.1-123.6); Figure 3 histogram spans roughly 10-170.",
        "Source-paper alias: 'CLCR'."
      ),
      source_name = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and V by linear regression on the MAP Bayesian estimates plus forward/backward elimination (Methods 2.2); not retained. Table 1: median 55 years (IQR 54-60)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened (Methods 2.2); not retained. Not summarised in Table 1."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 2.2); not retained. Table 1: median 77 kg (IQR 63-85). Volume and clearance are therefore absolute, not weight-scaled."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened as 'gender' (Methods 2.2); not retained. Table 1: 37 male / 24 female."
    )
  )

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 61L,
    n_studies = 1L,
    n_concentrations = 178L,
    age_median = "55 years (IQR 54-60)",
    weight_median = "77 kg (IQR 63-85)",
    sex_female_pct = 39.3,
    race_ethnicity = "Not reported (single-centre Italian cohort)",
    disease_state = paste(
      "Febrile neutropenia in adults with hematologic malignancies: acute myeloid leukemia 57.4%, lymphoma",
      "19.7%, acute lymphocytic leukemia 18.0%, multiple myeloma 4.9%. Clinically documented infection in",
      "49.2% (pneumonia 18.0%, bloodstream infection 16.4%); 7 documented Gram-negative isolates. 91.8% cured."
    ),
    renal_function = paste(
      "CKD-EPI CLCR median 107.3 mL/min/1.73 m^2 (IQR 96.1-123.6); 22.9% had augmented renal clearance",
      "(CLCR >= 130 mL/min/1.73 m^2)."
    ),
    dose_range = paste(
      "1 g loading dose over 30 min, then 1 g q8h by continuous infusion over 8 h (CLCR >= 60) or 0.5 g q6h",
      "by continuous infusion over 6 h (CLCR < 60), subsequently adjusted by real-time TDM to a steady-state",
      "target of 8-16 mg/L. Median dose 1 g q8h CI; median therapy 9 days (IQR 7-12.3)."
    ),
    regions = "Italy (Santa Maria della Misericordia University-Hospital of Udine)",
    notes = paste(
      "Demographics from Cojutti 2020 Table 1. 61 of the 100 patients of a prospective monocentric",
      "interventional TDM study (the other 39 excluded for inadequate sampling). 178 steady-state",
      "concentrations (median 3 TDM assessments per patient, from day 2-3 then every 48-72 h), HPLC,",
      "limit of detection 0.5 mg/L. Fitted with the Pmetrics NPAG algorithm (Pmetrics 1.5.0, R 3.4.4)."
    )
  )

  ini({
    # Structural parameters: Cojutti 2020 Table 2, 'Mean' column of the
    # Pmetrics NPAG non-parametric population distribution. Table 2 also
    # reports a 'Median' column (theta1 0.20, theta2 0.13, V 20.00). The mean
    # is used as the typical value, as in the sibling Pmetrics extraction
    # Braune_2018_meropenem.R; it is also the column the paper's own
    # summary sentence uses ('CL = 13.04 (4.85) L/h and V = 21.88 (5.85) L')
    # and the one that reproduces that mean clearance at the cohort's median
    # CLCR (0.27 + 0.12 * 107.3 = 13.15 L/h; the median column gives 14.15).
    lcl <- log(0.27)
    label("Clearance intercept theta1 at CRCL = 0 (L/h)") # Table 2: theta1 mean 0.27, SD 0.13, CV 48.53%, median 0.20
    e_crcl_cl <- 0.12
    label("Additive slope theta2 of CL on CRCL (L/h per mL/min/1.73 m^2)") # Table 2: theta2 mean 0.12, SD 0.03, CV 27.44%, median 0.13
    lvc <- log(21.88)
    label("Volume of distribution V (L)") # Table 2: V mean 21.88, SD 5.85, CV 26.71%, median 20.00

    # Interindividual variability. NPAG estimates a discrete joint density,
    # summarised in Table 2 by mean, SD and CV% for every estimated
    # parameter -- including the CRCL slope theta2, which therefore carries
    # its own random effect. Each CV% is carried into a log-normal random
    # effect with omega^2 = log(CV^2 + 1); correlations of the joint density
    # are not reported, so the etas are independent.
    #   theta1 : 48.53% CV -> log(0.4853^2 + 1) = 0.211397
    #   theta2 : 27.44% CV -> log(0.2744^2 + 1) = 0.072604
    #   V      : 26.71% CV -> log(0.2671^2 + 1) = 0.068940
    etalcl ~ 0.211397 # Table 2 theta1, CV 48.53%
    etae_crcl_cl ~ 0.072604 # Table 2 theta2, CV 27.44%
    etalvc ~ 0.068940 # Table 2 V, CV 26.71%

    # Residual error. Methods 2.2: 'A first-order polynomial relationship
    # between drug concentrations and the standard deviation of the
    # observations was used (C0 = 0.224, C1 = 0.060). Extra process noise
    # was captured with a gamma (G) model (G = 5).' Pmetrics' gamma model
    # multiplies the assay SD by gamma: SD = G * (C0 + C1 * C). The additive
    # and proportional parts are summed LINEARLY, which is nlmixr2's
    # combined1() form.
    #   addSd  = 5 * 0.224 = 1.12 mg/L
    #   propSd = 5 * 0.060 = 0.30
    addSd <- 1.12
    label("Additive residual SD, gamma * C0 (mg/L)") # Methods 2.2: C0 = 0.224, G = 5
    propSd <- 0.30
    label("Proportional residual SD, gamma * C1 (fraction)") # Methods 2.2: C1 = 0.060, G = 5
  })
  model({
    # Cojutti 2020 Results 3.2, Equation 1: CLi = theta1 + theta2 * CLCRi,
    # with each term drawn from its own NPAG marginal (Table 2).
    crcl_cl_i <- e_crcl_cl * exp(etae_crcl_cl)
    cl <- exp(lcl + etalcl) + crcl_cl_i * CRCL
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # One-compartment model with zero-order input and first-order
    # elimination (Methods 2.2).
    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
