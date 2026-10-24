Zhou_2016_avibactam <- function() {
  description <- paste(
    "Two-compartment IV population PK model for avibactam in adults, as",
    "tabulated by Zhou 2016 (Supplementary Table S4) for adult-to-pediatric",
    "extrapolation of renally cleared drugs. The model was developed by",
    "Li 2012 (ICAAC conference abstract) in healthy volunteers and adults with",
    "complicated intra-abdominal infection; Zhou 2016 removed its age and",
    "disease effects before simulation, leaving creatinine clearance as the",
    "only covariate, acting on clearance as a power function of CRCL/100",
    "whose exponent switches at 50 mL/min. Correlated IIV on CL and V1;",
    "combined proportional and additive residual error."
  )
  reference <- paste(
    "Zhou W, Johnson TN, Xu H, Cheung SYA, Bui KH, Li J, Al-Huniti N, Zhou D.",
    "Predictive Performance of Physiologically Based Pharmacokinetic and",
    "Population Pharmacokinetic Modeling of Renally Cleared Drugs in Children.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(9):475-483.",
    "doi:10.1002/psp4.12101. Parameter values are Supplementary Table S4",
    "('Summary of POPPK models', avibactam row). Original model: Li J,",
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
    central = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "avibactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (Zhou 2016 Supplementary Table S4, 'CRCL effect on CL (CRCL/100)^theta')",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw mL/min, NOT BSA-normalized: Table S4 prints no BSA term and centres the effect at 100 mL/min. The estimating equation used in the original adult model is not stated in Zhou 2016; for its pediatric simulations Zhou 2016 Methods ('Predictions in pediatrics') computed CRCL with the Schwartz equation for neonates and infants and the Rhodin equation from 2 years of age. The exponent on (CRCL/100) is 1.1 below 50 mL/min and 0.5 above it, applied literally as printed; the two arms do not meet at 50 mL/min (see the ini() comment on e_crcl_cl_lt50).",
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
    dose_range = "Adult validation study: avibactam 500 mg as a 2-hour IV infusion (Zhou 2016 Supplementary Table S1). Pediatric evaluation: 12.5 mg/kg (maximum 500 mg) as a 2-hour IV infusion in children aged 3 months to 18 years.",
    renal_function = "Creatinine clearance enters clearance with a reference of 100 mL/min; the cohort's CRCL range is not reported.",
    notes = "Zhou 2016 is a predictive-performance study: it did not fit this model but adapted an adult model from the literature (Li 2012 ICAAC abstract, not a full publication) and used it to predict clearance in children. Every value here is Zhou 2016's transcription of that model, with age and disease effects removed. N, demographics and the parameter uncertainties of the original fit are not reported in Zhou 2016. If the original analysis is later published in full, it should be extracted under its own author name."
  )

  ini({
    # Structural parameters: Zhou 2016 Supplementary Table S4, avibactam row
    lcl <- log(11.1); label("Clearance at CRCL = 100 mL/min (L/h)") # Table S4: CL 11.1 (L/h)
    lvc <- log(12.9); label("Central volume of distribution (L)") # Table S4: V1 12.9 (L)
    lq <- log(6.1); label("Inter-compartmental clearance (L/h)") # Table S4: Q 6.1 (L/h)
    lvp <- log(7.1); label("Peripheral volume of distribution (L)") # Table S4: V2 7.1 (L)

    # Table S4: 'CRCL effect on CL (CRCL/100)^theta' = '1.1 (CRCL<50) 0.5
    # (CRCL>50)'. Encoded literally. With a common reference of 100 mL/min
    # the two arms do not meet at the 50 mL/min switch: the factor is
    # 0.5^1.1 = 0.467 just below it and 0.5^0.5 = 0.707 just above.
    e_crcl_cl_lt50 <- 1.1; label("Exponent of (CRCL/100) on CL below 50 mL/min (unitless)") # Table S4: theta = 1.1 (CRCL<50)
    e_crcl_cl_ge50 <- 0.5; label("Exponent of (CRCL/100) on CL at or above 50 mL/min (unitless)") # Table S4: theta = 0.5 (CRCL>50)

    # IIV. Table S4 'IIV (% CV)' is the log-scale omega SD x 100: the
    # meropenem row of the same table (34.3/37.8/53.9/31.9) equals the square
    # roots of the Li 2006 meropenem omega variances (0.118/0.143/0.290/0.102).
    # 'Omega block between CL-V1 36.7%' is read as the CL-V1 correlation
    # coefficient r = 0.367, so cov = r * omega_CL * omega_V1.
    etalcl + etalvc ~ c(
      0.367^2,
      0.367 * 0.367 * 0.524,
      0.524^2
    ) # Table S4: IIV CL 36.7, V1 52.4 (% CV); 'Omega block between CL-V1 36.7%'

    # Residual error: Table S4 'P err 26.5% A err 14.5 (ng/mL)'
    propSd <- 0.265; label("Proportional residual error (fraction)") # Table S4: P err 26.5%
    addSd <- 0.0145; label("Additive residual error (mg/L)") # Table S4: A err 14.5 ng/mL = 0.0145 mg/L
  })
  model({
    # Renal function: power of CRCL/100 with an exponent that switches at
    # 50 mL/min (Table S4)
    e_crcl_cl <- e_crcl_cl_lt50 * (CRCL < 50) + e_crcl_cl_ge50 * (CRCL >= 50)

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
