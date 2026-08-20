Zheng_2012_mAb6 <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment population PK model",
    "for the non-disclosed humanized IgG antibody mAb6 (Zheng 2012;",
    "pI 9.2) with first-order SC absorption and linear central clearance",
    "in Gottingen minipigs (Table 1 population mean parameter estimates)."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (mAb6 row);",
    "minipig study details from Table 5; pI from Table 4."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "mAb6", units = NA_character_, specimen = "administration site", verified = FALSE),
    central     = list(analyte = "mAb6", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb6", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 10L,
    n_studies      = 1L,
    weight_range   = "IV: 8.9 kg mean; SC: 8.6 kg mean (Table 5)",
    sex_female_pct = 100,
    disease_state  = "Healthy Gottingen minipig; non-disclosed humanized IgG1 or stabilized IgG4 (approx 150 kDa) mAb preclinical PK.",
    dose_range     = "10 mg/kg IV bolus (n = 5); 120 mg SC inguinal (n = 5) (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5; pI = 9.2 (Table 4). Non-disclosed",
      "humanized IgG (IgG1 or stabilized IgG4) around 150 kDa. IIV and",
      "residual error were reported in prose (log-normal IIV; proportional",
      "plus additive residual error) but no variance magnitudes were",
      "tabulated, so this model encodes typical population values only;",
      "see the companion vignette Errata for details."
    )
  )

  ini({
    lka     <- log(0.530) ; label("First-order SC absorption rate Ka (1/day)")            # Table 1: Ka = 0.530 +/- 0.090 /day
    lcl     <- log(8.06)  ; label("Linear clearance CL (mL/day/kg)")                       # Table 1: CL = 8.06 +/- 0.81 mL/day/kg
    lvc     <- log(58.5)  ; label("Central volume of distribution Vc (mL/kg)")             # Table 1: Vc = 58.5 +/- 3.5 mL/kg
    lvp     <- log(26.0)  ; label("Peripheral volume of distribution Vp (mL/kg)")          # Table 1: Vp = 26.0 +/- 4.9 mL/kg
    lq      <- log(10.2)  ; label("Inter-compartmental clearance Q (mL/day/kg)")           # Table 1: Q  = 10.2 +/- 3.1 mL/day/kg
    lfdepot <- log(0.689) ; label("SC bioavailability / fraction absorbed F (unitless)")    # Table 1: F  = 68.9 +/- 5.5 %
  })

  model({
    ka     <- exp(lka)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- 1000 * central / vc
  })
}
