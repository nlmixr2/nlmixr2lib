Zheng_2012_mAb2 <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment population PK model",
    "for the non-disclosed humanized IgG antibody mAb2 (Zheng 2012;",
    "pI 8.7; unusually fast systemic clearance in minipig) with first-order",
    "SC absorption and linear central clearance in Gottingen minipigs",
    "(Table 1 population mean parameter estimates)."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (mAb2 row);",
    "minipig study details from Table 5; pI from Table 4."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 10L,
    n_studies      = 1L,
    weight_range   = "IV: 9.5 kg mean; SC: 9.7 kg mean (Table 5)",
    sex_female_pct = 100,
    disease_state  = "Healthy Gottingen minipig; non-disclosed humanized IgG1 or stabilized IgG4 (approx 150 kDa) mAb preclinical PK.",
    dose_range     = "5 mg/kg IV bolus (n = 5); 5 mg/kg SC scapular (n = 5) (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5; pI = 8.7 (Table 4). Non-disclosed",
      "humanized IgG (IgG1 or stabilized IgG4) around 150 kDa. mAb2 showed",
      "an unusually fast clearance in minipig (36.4 mL/day/kg) despite",
      "normal CL in human and cyno; excluded from the allometric-scaling",
      "analysis (Table 2). IIV and residual error were reported in prose",
      "(log-normal IIV; proportional plus additive residual error) but no",
      "variance magnitudes were tabulated, so this model encodes typical",
      "population values only; see the companion vignette Errata for",
      "details."
    )
  )

  ini({
    lka     <- log(0.316) ; label("First-order SC absorption rate Ka (1/day)")            # Table 1: Ka = 0.316 +/- 0.042 /day
    lcl     <- log(36.4)  ; label("Linear clearance CL (mL/day/kg)")                       # Table 1: CL = 36.4 +/- 2.2 mL/day/kg
    lvc     <- log(56.4)  ; label("Central volume of distribution Vc (mL/kg)")             # Table 1: Vc = 56.4 +/- 4.1 mL/kg
    lvp     <- log(134)   ; label("Peripheral volume of distribution Vp (mL/kg)")          # Table 1: Vp = 134 +/- 16 mL/kg
    lq      <- log(24.2)  ; label("Inter-compartmental clearance Q (mL/day/kg)")           # Table 1: Q  = 24.2 +/- 8.9 mL/day/kg
    lfdepot <- log(0.799) ; label("SC bioavailability / fraction absorbed F (unitless)")    # Table 1: F  = 79.9 +/- 9.4 %
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
