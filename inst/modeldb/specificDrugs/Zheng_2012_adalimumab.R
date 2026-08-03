Zheng_2012_adalimumab <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment population PK model",
    "for adalimumab with first-order SC absorption and linear central",
    "clearance in Gottingen minipigs, per Zheng 2012 minipig IV+SC study",
    "(Table 1 population mean parameter estimates)."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (adalimumab",
    "row); minipig study details from Table 5."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "adalimumab", units = NA_character_, specimen = "administration site", verified = FALSE),
    central     = list(analyte = "adalimumab", units = NA_character_, specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "adalimumab", units = NA_character_, specimen = "plasma", verified = FALSE)
  )

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 6L,
    n_studies      = 1L,
    weight_range   = "IV: 10.2 kg mean; SC: 10.0 kg mean (Table 5)",
    sex_female_pct = 100,
    disease_state  = "Healthy Gottingen minipig; adalimumab preclinical PK.",
    dose_range     = "40 mg IV bolus (n = 3); 40 mg SC inguinal (n = 3) (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5. Adalimumab formulated at 50 mg/mL",
      "(commercial Abbott product). Concentration-time profiles showed",
      "evidence of anti-therapeutic antibody (ATA) formation; ATA-affected",
      "periods were excluded from the PK fit per paper Methods and",
      "Discussion. IIV and residual error were reported in prose",
      "(log-normal IIV; proportional plus additive residual error) but no",
      "variance magnitudes were tabulated, so this model encodes typical",
      "population values only (no eta / residual error blocks); see the",
      "companion vignette Errata for details."
    )
  )

  ini({
    # Structural population means (Zheng 2012 Table 1, adalimumab row)
    lka     <- log(4.61)  ; label("First-order SC absorption rate Ka (1/day)")            # Table 1: Ka = 4.61 +/- 0.65 /day
    lcl     <- log(3.48)  ; label("Linear clearance CL (mL/day/kg)")                       # Table 1: CL = 3.48 +/- 0.49 mL/day/kg
    lvc     <- log(55.4)  ; label("Central volume of distribution Vc (mL/kg)")             # Table 1: Vc = 55.4 +/- 7.8 mL/kg
    lvp     <- log(13.5)  ; label("Peripheral volume of distribution Vp (mL/kg)")          # Table 1: Vp = 13.5 +/- 5.1 mL/kg
    lq      <- log(24.0)  ; label("Inter-compartmental clearance Q (mL/day/kg)")           # Table 1: Q  = 24.0 +/- 18.2 mL/day/kg
    lfdepot <- log(0.829) ; label("SC bioavailability / fraction absorbed F (unitless)")    # Table 1: F  = 82.9 +/- 8.3 %
  })

  model({
    ka     <- exp(lka)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    fdepot <- exp(lfdepot)

    # Rate constants: mL/day/kg over mL/kg -> 1/day
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Paper Eq. 1: dAabs/dt = -Ka*Aabs; Vc*dCc/dt = -Q(Cc-Cp) - CL*Cc + F*Ka*Aabs; Vp*dCp/dt = +Q(Cc-Cp)
    # Encoded in mass (mg/kg) space with F applied via f(depot).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # Concentration units: dose (mg/kg) / (Vc mL/kg / 1000 mL/L) = mg/L = ug/mL
    Cc <- 1000 * central / vc
  })
}
