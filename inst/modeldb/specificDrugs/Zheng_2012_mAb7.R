Zheng_2012_mAb7 <- function() {
  description <- paste(
    "Preclinical (Gottingen minipig). Two-compartment population PK model",
    "with first-order SC absorption plus parallel linear and saturable",
    "(Michaelis-Menten) central elimination for the non-disclosed humanized",
    "IgG antibody mAb7 (Zheng 2012; pI 8.9) in Gottingen minipigs (Table 1",
    "population mean parameter estimates)."
  )
  reference <- paste(
    "Zheng Y, Tesar DB, Benincosa L, et al. Minipig as a potential",
    "translatable model for monoclonal antibody pharmacokinetics after",
    "intravenous and subcutaneous administration. mAbs. 2012;4(2):243-255.",
    "doi:10.4161/mabs.4.2.19387. Parameter values from Table 1 (mAb7 row,",
    "including nonlinear Vmax / Km from footnote g); minipig study details",
    "from Table 5; pI from Table 4."
  )
  vignette <- "Zheng_2012_minipig_mab"
  units <- list(time = "day", dosing = "mg/kg", concentration = "ug/mL")

  covariateData <- list()

  population <- list(
    species        = "Gottingen minipig",
    n_subjects     = 10L,
    n_studies      = 1L,
    weight_range   = "IV: 8.4 kg mean; SC: 8.2 kg mean (Table 5)",
    sex_female_pct = 100,
    disease_state  = "Healthy Gottingen minipig; non-disclosed humanized IgG1 or stabilized IgG4 (approx 150 kDa) mAb preclinical PK.",
    dose_range     = "9 mg/kg IV bolus (n = 5); 108 mg SC inguinal (n = 5) (Table 5)",
    regions        = "Contract research organizations (Denmark, UK, USA)",
    notes          = paste(
      "Zheng 2012 Tables 1 and 5; pI = 8.9 (Table 4). Non-disclosed",
      "humanized IgG (IgG1 or stabilized IgG4) around 150 kDa. mAb7",
      "exhibited nonlinear pharmacokinetics; per paper Eq. 2 the effective",
      "clearance is the sum of a linear component (CL_lin) and a saturable",
      "Michaelis-Menten arm (Vmax / (Km + Cc)). Table 1 tabulates only the",
      "linear CL for mAb7; the nonlinear Vmax and Km are given in the",
      "Table 1 footnote (g). T1/2 beta is not reported for mAb7 (Table 1",
      "footnote h) because half-life is not a well-defined constant under",
      "the nonlinear model. IIV and residual error were reported in prose",
      "(log-normal IIV; proportional plus additive residual error) but no",
      "variance magnitudes were tabulated, so this model encodes typical",
      "population values only; see the companion vignette Errata for",
      "details."
    )
  )

  ini({
    lka     <- log(0.828) ; label("First-order SC absorption rate Ka (1/day)")            # Table 1: Ka = 0.828 +/- 0.157 /day
    lcl     <- log(5.06)  ; label("Linear clearance component CL_lin (mL/day/kg)")         # Table 1: CL_lin = 5.06 +/- 0.46 mL/day/kg (footnote g: nonlinear PK; only linear CL tabulated)
    lvc     <- log(61.6)  ; label("Central volume of distribution Vc (mL/kg)")             # Table 1: Vc = 61.6 +/- 2.5 mL/kg
    lvp     <- log(37.4)  ; label("Peripheral volume of distribution Vp (mL/kg)")          # Table 1: Vp = 37.4 +/- 3.7 mL/kg
    lq      <- log(17.6)  ; label("Inter-compartmental clearance Q (mL/day/kg)")           # Table 1: Q  = 17.6 +/- 3.0 mL/day/kg
    lfdepot <- log(0.854) ; label("SC bioavailability / fraction absorbed F (unitless)")    # Table 1: F  = 85.4 +/- 2.6 %
    lvmax   <- log(165)   ; label("Michaelis-Menten maximum elimination rate Vmax (ug/day/kg)")  # Table 1 footnote g: Vmax = 165 ug/day/kg
    lkm     <- log(8.27)  ; label("Michaelis-Menten half-saturation concentration Km (ug/mL)")   # Table 1 footnote g: Km = 8.27 ug/mL
  })

  model({
    ka     <- exp(lka)
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    vp     <- exp(lvp)
    q      <- exp(lq)
    fdepot <- exp(lfdepot)
    vmax   <- exp(lvmax)   # ug/day/kg
    km     <- exp(lkm)     # ug/mL

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Concentration in ug/mL from mg/kg central amount and mL/kg volume.
    Cc <- 1000 * central / vc

    # Paper Eq. 2 nonlinear elimination in concentration space:
    #   CL_eff * Cc = CL * Cc + Vmax * Cc / (Km + Cc)  [ug/day/kg]
    # Convert the MM term from ug/day/kg to mg/day/kg for the amount
    # balance on central [mg/kg].
    mm_rate_mg <- vmax * Cc / (km + Cc) / 1000

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - mm_rate_mg
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot
  })
}
