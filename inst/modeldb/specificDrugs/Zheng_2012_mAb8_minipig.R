Zheng_2012_mAb8_minipig <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment IV-only population",
    "PK model for the non-disclosed humanized IgG antibody mAb8",
    "(Zheng 2012; pI 8.7) with linear central clearance in Gottingen",
    "minipigs (Table 1 population mean parameter estimates). No SC arm",
    "was run for mAb8, so no Ka or F is defined."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (mAb8 row);",
    "minipig study details from Table 5; pI from Table 4."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "mAb8", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb8", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 5L,
    n_studies      = 1L,
    weight_range   = "IV: 9.7 kg mean (Table 5); SC arm not conducted",
    sex_female_pct = 0,
    disease_state  = "Healthy Gottingen minipig; non-disclosed humanized IgG1 or stabilized IgG4 (approx 150 kDa) mAb preclinical PK (IV only).",
    dose_range     = "10 mg/kg IV bolus (n = 5); no SC arm (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5; pI = 8.7 (Table 4). Non-disclosed",
      "humanized IgG (IgG1 or stabilized IgG4) around 150 kDa. Male",
      "cohort. IV-only study; Ka, F and Tmax are 'not available' in Table 1.",
      "IIV and residual error were reported in prose (log-normal IIV;",
      "proportional plus additive residual error) but no variance magnitudes",
      "were tabulated, so this model encodes typical population values",
      "only; see the companion vignette Errata for details."
    )
  )

  ini({
    lcl <- log(2.45) ; label("Linear clearance CL (mL/day/kg)")                            # Table 1: CL = 2.45 +/- 0.12 mL/day/kg
    lvc <- log(36.0) ; label("Central volume of distribution Vc (mL/kg)")                  # Table 1: Vc = 36.0 +/- 1.1 mL/kg
    lvp <- log(36.4) ; label("Peripheral volume of distribution Vp (mL/kg)")               # Table 1: Vp = 36.4 +/- 1.8 mL/kg
    lq  <- log(27.6) ; label("Inter-compartmental clearance Q (mL/day/kg)")                # Table 1: Q  = 27.6 +/- 2.8 mL/day/kg
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV-only 2-compartment model; no absorption depot.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- 1000 * central / vc
  })
}
