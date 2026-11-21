Yin_2024_erlotinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral erlotinib in adults with",
    "EGFR-mutant non-small cell lung cancer (Yin 2024): first-order",
    "absorption with a lag time and first-order elimination, fitted to",
    "intensive steady-state profiles from two food-interaction studies",
    "pooled with sparse real-world trough samples from the START-TKI study.",
    "Bioavailability is fixed to 1 with an IIV on F that is shared by all",
    "apparent clearances and volumes. No covariate was retained. The",
    "paper's companion tumor-dynamics models are",
    "Yin_2024_erlotinib_tumor and Yin_2024_erlotinib_tumor_ctdna.",
    sep = " "
  )
  reference <- paste(
    "Yin A, Veerman GDM, van Hasselt JGC, Steendam CMJ, Dubbink HJ,",
    "Guchelaar HJ, Friberg LE, Dingemans AMC, Mathijssen RHJ, Moes DJAR.",
    "Quantitative modeling of tumor dynamics and development of drug",
    "resistance in non-small cell lung cancer patients treated with",
    "erlotinib. CPT Pharmacometrics Syst Pharmacol. 2024;13(4):612-623. doi:10.1002/psp4.13105",
    sep = " "
  )
  vignette <- "Yin_2024_erlotinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  compartmentData <- list(
    depot = list(analyte = "erlotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "erlotinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "erlotinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 45,
    n_studies = 3,
    age_range = "35-78 years",
    age_median = "63 years (intensive PK cohort); 66 years (START-TKI cohort)",
    weight_range = "46.1-109 kg",
    weight_median = "74 kg (intensive PK cohort); 69.5 kg (START-TKI cohort)",
    sex_female_pct = 61.7,
    disease_state = "Non-small cell lung cancer with an activating EGFR mutation, treated with erlotinib",
    dose_range = "50-150 mg orally once daily (150 mg in 86% of the intensive cohort and all START-TKI patients)",
    regions = "The Netherlands (Erasmus MC Cancer Institute and the START-TKI multicentre real-world study)",
    notes = paste(
      "Yin 2024 Table 1. 29 patients contributed intensive steady-state",
      "profiles (13 samples each, predose to 24 h) from the control arms of",
      "two food-interaction studies; 18 START-TKI patients contributed",
      "146 sparse samples (2 of the 18 had no PK data). Sex percentage is",
      "pooled over all 47 patients (29 of 47 female)."
    )
  )

  ini({
    lka <- log(1.61)
    label("Absorption rate constant ka (1/h)") # Table 2 'Ka (/h)' = 1.61 (RSE 23%)
    lcl <- log(4.10)
    label("Apparent clearance CL/F (L/h)") # Table 2 'CL/F (L/h)' = 4.10 (RSE 5%)
    lvc <- log(142)
    label("Apparent central volume Vc/F (L)") # Table 2 'Vc/F (L)' = 142 (RSE 7%)
    lvp <- log(2420)
    label("Apparent peripheral volume Vp/F (L)") # Table 2 'Vp/F (L)' = 2420 (RSE 12%)
    lq <- log(0.548)
    label("Apparent inter-compartmental clearance Q/F (L/h)") # Table 2 'Q/F (L/h)' = 0.548 (RSE 40%)
    ltlag <- log(0.400)
    label("Absorption lag time (h)") # Table 2 'Tlag (h)' = 0.400 (RSE 5%)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability F (unitless)") # Table 2 'F' = 1 fixed; Methods 'the typical value of F was fixed to 1'

    # Table 2 reports IIV as CV%; converted to log-normal variances as
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.024350 # Table 2 IIV CL/F 15.7% CV -> log(1 + 0.157^2)
    etalvc ~ 0.040383 # Table 2 IIV Vc/F 20.3% CV -> log(1 + 0.203^2)
    etalq ~ 1.564260 # Table 2 IIV Q/F 194.4% CV -> log(1 + 1.944^2)
    etalka ~ 0.936103 # Table 2 IIV Ka 124.5% CV -> log(1 + 1.245^2)
    etalfdepot ~ 0.026222 # Table 2 IIV F 16.3% CV -> log(1 + 0.163^2)

    propSd <- 0.154
    label("Proportional residual error (fraction)") # Table 2 'Prop. err. (CV%)' = 15.4 (RSE 6%)
    addSd <- 44.5
    label("Additive residual error (ng/mL)") # Table 2 'Add. err. (SD, ng/mL)' = 44.5 (RSE 25%)
  })
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq + etalq)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag
    f(depot) <- fdepot

    # Dose in mg and volume in L give mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
