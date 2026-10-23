Zhou_2016_ceftazidime <- function() {
  description <- paste(
    "Two-compartment IV population PK model for ceftazidime in adults, as",
    "tabulated by Zhou 2016 (Supplementary Table S4) for adult-to-pediatric",
    "extrapolation of renally cleared drugs. The model was developed by",
    "Li 2012 (ICAAC conference abstract) in healthy volunteers and adults with",
    "complicated intra-abdominal infection; Zhou 2016 removed its age and",
    "disease effects before simulation, leaving creatinine clearance as the",
    "only covariate, acting on clearance as a power function of CRCL/100.",
    "Correlated IIV on CL and V1; combined proportional and additive",
    "residual error."
  )
  reference <- paste(
    "Zhou W, Johnson TN, Xu H, Cheung SYA, Bui KH, Li J, Al-Huniti N, Zhou D.",
    "Predictive Performance of Physiologically Based Pharmacokinetic and",
    "Population Pharmacokinetic Modeling of Renally Cleared Drugs in Children.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(9):475-483.",
    "doi:10.1002/psp4.12101. Parameter values are Supplementary Table S4",
    "('Summary of POPPK models', ceftazidime row). Original model: Li J,",
    "Knebel W, Riggs M, Zhou D, Nichols WW, Das S. Population pharmacokinetic",
    "modeling of ceftazidime and avibactam in healthy subjects and in patients",
    "with complicated intra-abdominal infection (cIAI). 52nd Interscience",
    "Conference on Antimicrobial Agents and Chemotherapy (ICAAC), San",
    "Francisco, CA, 2012 (Zhou 2016 reference 31).",
    sep = " "
  )
  vignette <- "Zhou_2016_ceftazidime_avibactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (Zhou 2016 Supplementary Table S4, 'CRCL effect on CL (CRCL/100)^theta')",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw mL/min, NOT BSA-normalized: Table S4 prints no BSA term and centres the effect at 100 mL/min. The estimating equation used in the original adult model is not stated in Zhou 2016; for its pediatric simulations Zhou 2016 Methods ('Predictions in pediatrics') computed CRCL with the Schwartz equation for neonates and infants and the Rhodin equation from 2 years of age.",
      source_name = "CRCL"
    )
  )

  # Covariates of the original Li 2012 model that Zhou 2016 removed before
  # simulating ('Age and disease effects removed for simulations', Table S4
  # comment column). Their functional forms and estimates are not printed.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant in the original adult model; removed by Zhou 2016 ('The age effect was fixed to one if the age was a significant covariate in the adult PopPK model', Methods, 'PopPK models in adults'). Form and estimate not reported.",
      source_name = "Age"
    ),
    DIS_CIAI = list(
      description = "Complicated intra-abdominal infection indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = healthy volunteer",
      notes = "The original model pooled healthy volunteers and cIAI patients (Li 2012 title) and carried a disease effect, which Zhou 2016 removed. Which parameters it acted on and its estimate are not reported.",
      source_name = "disease"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "Healthy adult volunteers and adults with complicated intra-abdominal infection (original Li 2012 model). As used by Zhou 2016, the disease effect is removed, so the model describes the reference population of the original fit.",
    age_range = "Adults; the original model's age range is not reported. The adult validation study Zhou 2016 simulated (Supplementary Table S1) enrolled 16 healthy males aged 22-46 years.",
    dose_range = "Adult validation study: ceftazidime 2000 mg as a 2-hour IV infusion (Zhou 2016 Supplementary Table S1). Pediatric evaluation: 50 mg/kg (maximum 2000 mg) as a 2-hour IV infusion in children aged 3 months to 18 years.",
    renal_function = "Creatinine clearance enters clearance with a reference of 100 mL/min; the cohort's CRCL range is not reported.",
    notes = "Zhou 2016 is a predictive-performance study: it did not fit this model but adapted an adult model from the literature (Li 2012 ICAAC abstract, not a full publication) and used it to predict clearance in children. Every value here is Zhou 2016's transcription of that model, with age and disease effects removed. N, demographics and the parameter uncertainties of the original fit are not reported in Zhou 2016. If the original analysis is later published in full, it should be extracted under its own author name."
  )

  ini({
    # Structural parameters: Zhou 2016 Supplementary Table S4, ceftazidime row
    lcl <- log(7.4); label("Clearance at CRCL = 100 mL/min (L/h)") # Table S4: CL 7.4 (L/h)
    lvc <- log(13.3); label("Central volume of distribution (L)") # Table S4: V1 13.3 (L)
    lq <- log(20.5); label("Inter-compartmental clearance (L/h)") # Table S4: Q 20.5 (L/h)
    lvp <- log(9.0); label("Peripheral volume of distribution (L)") # Table S4: V2 9.0 (L)

    e_crcl_cl <- 0.4; label("Exponent of (CRCL/100) on CL (unitless)") # Table S4: 'CRCL effect on CL (CRCL/100)^theta' = 0.4

    # IIV. Table S4 'IIV (% CV)' is the log-scale omega SD x 100 (see the
    # meropenem cross-check in Zhou_2016_avibactam.R). 'Omega block between
    # CL-V1 47.5%' is read as the CL-V1 correlation coefficient r = 0.475,
    # so cov = r * omega_CL * omega_V1.
    etalcl + etalvc ~ c(
      0.358^2,
      0.475 * 0.358 * 0.702,
      0.702^2
    ) # Table S4: IIV CL 35.8, V1 70.2 (% CV); 'Omega block between CL-V1 47.5%'

    # Residual error: Table S4 'P err 41.3% A err 79.2 (ng/mL)'
    propSd <- 0.413; label("Proportional residual error (fraction)") # Table S4: P err 41.3%
    addSd <- 0.0792; label("Additive residual error (mg/L)") # Table S4: A err 79.2 ng/mL = 0.0792 mg/L
  })
  model({
    cl <- exp(lcl + etalcl) * (CRCL / 100)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> total plasma concentration in mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
