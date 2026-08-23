Zheng_2012_mAb3_minipig <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment population PK model",
    "for the non-disclosed humanized IgG antibody mAb3 (Zheng 2012;",
    "pI 9.1) with first-order SC absorption and linear central clearance",
    "in Gottingen minipigs (Table 1 population mean parameter estimates)."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (mAb3 row);",
    "minipig study details from Table 5; pI from Table 4."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "mAb3", units = NA_character_, specimen = "administration site", verified = FALSE),
    central     = list(analyte = "mAb3", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb3", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 10L,
    n_studies      = 1L,
    weight_range   = "IV: 8.2 kg mean; SC: 8.3 kg mean (Table 5)",
    sex_female_pct = 100,
    disease_state  = "Healthy Gottingen minipig; non-disclosed humanized IgG1 or stabilized IgG4 (approx 150 kDa) mAb preclinical PK.",
    dose_range     = "5 mg/kg IV bolus (n = 5); 5 mg/kg SC scapular (n = 5) (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5; pI = 9.1 (Table 4). Non-disclosed",
      "humanized IgG (IgG1 or stabilized IgG4) around 150 kDa. IIV and",
      "residual error were reported in prose (log-normal IIV; proportional",
      "plus additive residual error) but no variance magnitudes were",
      "tabulated, so this model encodes typical population values only;",
      "see the companion vignette Errata for details."
    )
  )

  ini({
    lka     <- log(0.635) ; label("First-order SC absorption rate Ka (1/day)")            # Table 1: Ka = 0.635 +/- 0.244 /day
    lcl     <- log(4.48)  ; label("Linear clearance CL (mL/day/kg)")                       # Table 1: CL = 4.48 +/- 0.75 mL/day/kg
    lvc     <- log(36.6)  ; label("Central volume of distribution Vc (mL/kg)")             # Table 1: Vc = 36.6 +/- 13.5 mL/kg
    lvp     <- log(70.9)  ; label("Peripheral volume of distribution Vp (mL/kg)")          # Table 1: Vp = 70.9 +/- 12.9 mL/kg
    lq      <- log(72.5)  ; label("Inter-compartmental clearance Q (mL/day/kg)")           # Table 1: Q  = 72.5 +/- 18.0 mL/day/kg
    lfdepot <- log(0.658) ; label("SC bioavailability / fraction absorbed F (unitless)")    # Table 1: F  = 65.8 +/- 14.4 %
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
