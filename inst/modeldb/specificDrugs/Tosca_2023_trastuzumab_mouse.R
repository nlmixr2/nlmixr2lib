Tosca_2023_trastuzumab_mouse <- function() {
  description <- paste(
    "Preclinical (mouse; tumor-bearing SCID mice with KPL-4 and nude mice with",
    "SKOV3 xenografts). Two-compartment PK model with first-order absorption",
    "for trastuzumab after IV and intraperitoneal administration, fitted by",
    "a naive-pooled approach to literature radiolabel biodistribution data",
    "converted to plasma concentrations. Proportional residual error and no",
    "between-animal variability; volumes and clearances are whole-mouse values."
  )
  reference <- paste(
    "Tosca EM, Borella E, Piana C, Bouchene S, Merlino G, Fiascarelli A,",
    "Mazzei P, Magni P. Model-based prediction of",
    "effective target exposure for MEN1611 in combination with trastuzumab in",
    "HER2-positive advanced or metastatic breast cancer patients.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1626-1639.",
    "doi:10.1002/psp4.12910. Parameter values from Supplementary Table S3.",
    sep = " "
  )
  vignette <- "Tosca_2023_men1611_trastuzumab"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "trastuzumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "trastuzumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "trastuzumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "mouse (SCID bearing KPL-4 tumors; nude bearing SKOV3 tumors)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    weight_range = "22.5 g assumed typical body weight for the dose conversion",
    disease_state = "tumor-bearing mice (KPL-4 breast cancer and SKOV3 ovarian cancer xenografts)",
    dose_range = "Study 1: single IV 0.1, 1.4 or 17 mg/kg unlabelled trastuzumab plus 111In- and 125I-trastuzumab (3.43, 32.68 and 383.68 ug total). Study 2: single IP 86Y-trastuzumab (17.65 ug).",
    regions = "preclinical (literature data)",
    notes = paste(
      "TZB-PK study 1 (blood %ID/g over 168 h, three replicates per dose and",
      "radioisotope) and TZB-PK study 2 (%IA/g digitised over 3 days) were",
      "taken from the literature and converted to plasma ng/mL as described in",
      "Tosca 2023 Supplementary Material S1 (blood density 1 g/mL, no",
      "blood-to-plasma correction). Animal numbers are not reported."
    )
  )

  ini({
    # Tosca 2023 Supplementary Table S3, 'Estimates of TZB PK model parameters
    # in TB mice'. Volumes and clearances are whole-mouse values (doses were
    # converted to ug for a 22.5 g mouse; Material S1).
    lka <- log(0.151); label("Intraperitoneal absorption rate constant (1/h)") # Table S3: Kabs = 0.151 1/h (7% CV)
    lcl <- log(1.53e-5); label("Clearance (L/h)") # Table S3: CL = 1.53E-5 L/h (7% CV)
    lvc <- log(2.45e-3); label("Central volume of distribution (L)") # Table S3: V1 = 2.45E-3 L (3% CV)
    lq <- log(2.58e-4); label("Intercompartmental clearance (L/h)") # Table S3: Q = 2.58e-4 L/h (10% CV)
    lvp <- log(2.22e-3); label("Peripheral volume of distribution (L)") # Table S3: V2 = 2.22e-3 L (7% CV)
    lfdepot <- log(0.5557); label("Intraperitoneal bioavailability (fraction)") # Table S3: F = 0.5557 (4% CV)

    propSd <- 0.095; label("Proportional residual error (fraction)") # Table S3: b = 0.095 (9% CV)
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intraperitoneal doses go to depot; IV doses go to central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    # central is in mg and vc in L, so central / vc is mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
