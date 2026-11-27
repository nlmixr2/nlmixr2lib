Gomeni_2019_methylphenidate_ir <- function() {
  description <- paste(
    "One-compartment PK model with first-order absorption, an absorption lag",
    "time and first-order elimination for immediate-release methylphenidate",
    "(Ritalin IR tablets) in healthy adults (Gomeni 2019, step 1 of the",
    "IVIVC analysis). Fitted to the mean concentration-time profile after",
    "two 20 mg tablets given 4 h apart, so there is no between-subject",
    "variability; additive residual error. Its elimination rate and volume",
    "are the disposition fixed in the convolution IVIVC model",
    "(Gomeni_2019_methylphenidate_ivivc)."
  )
  reference <- paste(
    "Gomeni R, Fang L, Bressolle-Gomeni F, Spencer TJ, Faraone SV, Babiskin A.",
    "A General Framework for Assessing In vitro/In vivo Correlation as a Tool",
    "for Maximizing the Benefit-Risk Ratio of a Treatment Using a",
    "Convolution-Based Modeling Approach. CPT Pharmacometrics Syst Pharmacol.",
    "2019;8(2):97-106. doi:10.1002/psp4.12378. Parameter values from Table S1",
    "(kel to three decimals from the NONMEM control stream)."
  )
  vignette <- "Gomeni_2019_methylphenidate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "methylphenidate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "methylphenidate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    studies = paste(
      "Wang Y et al. Biopharm Drug Dispos 2004;25:91-98: Ritalin IR reference",
      "arm (two 20 mg tablets 4 h apart) of a four-period single-dose",
      "crossover in healthy volunteers after at least a 10-hour fast."
    ),
    disease_state = "Healthy volunteers",
    dose_range = "Two 20 mg IR tablets given 4 h apart (40 mg total)",
    regions = "Not reported",
    notes = paste(
      "Fitted to the published mean profile digitised from Wang 2004. Age,",
      "weight and sex are not reported in Gomeni 2019. All Table S1 RSEs are",
      "printed as <0.01%."
    )
  )

  ini({
    lka <- log(2.39); label("First-order absorption rate constant ka (1/h)") # Table S1: ka = 2.39 1/h
    ltlag <- log(0.63); label("Absorption lag time (h)") # Table S1: Lag = 0.63 h
    lkel <- log(0.254); label("Elimination rate constant kel (1/h)") # Table S1: kel = 0.25 1/h; the control stream fixes KEL = 0.254 from this fit
    lvc <- log(1200); label("Volume of distribution V (L)") # Table S1: V = 1.20 in kL (dose mg, concentration ng/mL; control stream S1 = 1.2), i.e. 1200 L
    addSd <- 0.60; label("Additive residual error (ng/mL)") # Table S1: err_add = 0.60
  })

  model({
    ka <- exp(lka)
    tlag <- exp(ltlag)
    kel <- exp(lkel)
    vc <- exp(lvc)

    d/dt(depot) <- -ka * depot
    lag(depot) <- tlag
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V in L give mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
