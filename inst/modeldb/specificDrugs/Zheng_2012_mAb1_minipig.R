Zheng_2012_mAb1_minipig <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment population PK model",
    "for the non-disclosed humanized IgG antibody mAb1 (Zheng 2012;",
    "pI 6.1) with first-order SC absorption and linear central clearance",
    "in Gottingen minipigs (Table 1 population mean parameter estimates)."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (mAb1 row);",
    "minipig study details from Table 5; pI from Table 4."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "mAb1", units = NA_character_, specimen = "administration site", verified = FALSE),
    central     = list(analyte = "mAb1", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "mAb1", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 9L,
    n_studies      = 1L,
    weight_range   = "IV: 9.4 kg mean; SC: 9.9 kg mean (Table 5)",
    sex_female_pct = 100,
    disease_state  = "Healthy Gottingen minipig; non-disclosed humanized IgG1 or stabilized IgG4 (approx 150 kDa) mAb preclinical PK.",
    dose_range     = "5 mg/kg IV bolus (n = 4); 5 mg/kg SC scapular (n = 5) (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5; pI = 6.1 (Table 4). Non-disclosed",
      "humanized IgG (IgG1 or stabilized IgG4) around 150 kDa. IIV and",
      "residual error were reported in prose (log-normal IIV; proportional",
      "plus additive residual error) but no variance magnitudes were",
      "tabulated, so this model encodes typical population values only",
      "(no eta / residual error blocks); see the companion vignette Errata",
      "for details."
    )
  )

  ini({
    lka     <- log(1.02)  ; label("First-order SC absorption rate Ka (1/day)")            # Table 1: Ka = 1.02 +/- 0.33 /day
    lcl     <- log(2.83)  ; label("Linear clearance CL (mL/day/kg)")                       # Table 1: CL = 2.83 +/- 0.60 mL/day/kg
    lvc     <- log(51.9)  ; label("Central volume of distribution Vc (mL/kg)")             # Table 1: Vc = 51.9 +/- 2.2 mL/kg
    lvp     <- log(54.6)  ; label("Peripheral volume of distribution Vp (mL/kg)")          # Table 1: Vp = 54.6 +/- 9.9 mL/kg
    lq      <- log(128)   ; label("Inter-compartmental clearance Q (mL/day/kg)")           # Table 1: Q  = 128 +/- 25 mL/day/kg
    lfdepot <- log(0.975) ; label("SC bioavailability / fraction absorbed F (unitless)")    # Table 1: F  = 97.5 +/- 15.1 %
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

    # Paper Eq. 1: 2-compartment linear model with first-order SC absorption.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Cc [ug/mL] = 1000 * central [mg/kg] / vc [mL/kg]  (i.e. central / (vc/1000))
    Cc <- 1000 * central / vc
  })
}
