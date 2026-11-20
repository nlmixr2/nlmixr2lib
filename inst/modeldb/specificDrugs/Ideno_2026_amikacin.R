Ideno_2026_amikacin <- function() {
  description <- "Two-compartment population PK model for intravenous amikacin without covariates, fitted by Ideno 2026 to the publicly available University of Strathclyde amikacin therapeutic-drug-monitoring dataset (Thomson 2021) as the base model for the SCOUT individual-optimal-dose covariate screen (Part 2). Exponential IIV on CL, V1 and V2 (none on Q) and a combined proportional + additive residual error."
  reference <- paste(
    "Ideno Y, Kasai H, Tanigawara Y.",
    "SCOUT: An Exploratory Approach to Scouting Dose-Relevant Covariates.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(4):e70235.",
    "doi:10.1002/psp4.70235.",
    "Parameter estimates from Data S2 (population PK model for amikacin, Part 2).",
    "Source data: Thomson A, Robb F, Siebinga H. Amikacin Data for Population",
    "Pharmacokinetic Modelling. University of Strathclyde; 2021.",
    "doi:10.15129/d5b80e9f-2f49-45d6-9432-c23470f4e851.",
    sep = " "
  )
  vignette <- "Ideno_2026_scout"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds. The SCOUT control stream (Data S1, Part 2)
  # doses amikacin into CENTRAL and divides by V1 to give mg/L. The sampling
  # matrix of the Strathclyde dataset is not described in Ideno 2026.
  compartmentData <- list(
    central = list(analyte = "amikacin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "amikacin", units = "mg", specimen = "tissue", verified = FALSE)
  )

  # No covariates were included in the population model ("No covariate
  # modeling was performed -- no covariates were included in the model from
  # the outset", Ideno 2026 Methods, Part 2). SEX, AGE, WT, IBW, AJBW, HT,
  # SCREAT and CRCL were carried in the SCOUT dataset only to be plotted
  # against the estimated optimal doses (Figure 3).
  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "approximately 18-90 years (read from the x-axis of Ideno 2026 Figure 3)",
    weight_range = "adjusted body weight approximately 35-93 kg (Ideno 2026 Figure 3)",
    sex_female_pct = NA_real_,
    disease_state = "Adult patients receiving intravenous amikacin with therapeutic drug monitoring (publicly available University of Strathclyde dataset, Thomson 2021)",
    dose_range = "Clinical amikacin doses from the source dataset; SCOUT target was Cmax 60 mg/L after a 30-min IV infusion",
    renal_function = "Creatinine clearance approximately 20-185 mL/min (Ideno 2026 Figure 3)",
    notes = "Ideno 2026 does not tabulate the size or demographics of the Strathclyde cohort; the ranges above are read from the axes of Figure 3 (one point per patient, both sexes represented). The number of subjects is not stated in the text."
  )

  ini({
    lcl <- log(3.15); label("Clearance CL (L/h)")                          # Data S2 table: CL (L/hr) = 3.15 (RSE 0.0740)
    lvc <- log(19.3); label("Central volume of distribution V1 (L)")       # Data S2 table: V1 (L) = 19.3 (RSE 0.0404)
    lq  <- log(0.692); label("Intercompartmental clearance Q (L/h)")       # Data S2 table: Q (L/hr) = 0.692 (RSE 0.159)
    lvp <- log(7.65); label("Peripheral volume of distribution V2 (L)")    # Data S2 table: V2 (L) = 7.65 (RSE 0.506)

    # Data S2 reports Omega-CL/V1/V2 as SDs of an exponential IIV model;
    # nlmixr2 takes the variance, so each value below is SD^2.
    etalcl ~ 0.431649  # Data S2 table: Omega-CL SD = 0.657 (shrinkage 0.112); variance 0.657^2
    etalvc ~ 0.0961    # Data S2 table: Omega-V1 SD = 0.310 (shrinkage 0.239); variance 0.310^2
    etalvp ~ 0.651249  # Data S2 table: Omega-V2 SD = 0.807 (shrinkage 0.770); variance 0.807^2

    propSd <- 0.170; label("Proportional residual error (fraction)")  # Data S2 table: RUV prop = 0.170 (RSE 0.0600)
    addSd  <- 0.994; label("Additive residual error (mg/L)")          # Data S2 table: RUV add = 0.994 (RSE 0.175)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Data S1 Part 2: DADT(1) = -(CL+Q)*CC + Q*CP; DADT(2) = Q*(CC - CP)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
