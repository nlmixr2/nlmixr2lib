Delattre_2010_amikacin <- function() {
  description <- "Two-compartment IV population PK model for amikacin in critically ill adult patients with severe sepsis or septic shock during the first 24 hours of antibiotic treatment (Delattre 2010)"
  reference <- "Delattre IK, Musuamba FT, Nyberg J, Taccone FS, Laterre PF, Verbeeck RK, Jacobs F, Wallemacq PE. Population pharmacokinetic modeling and optimal sampling strategy for Bayesian estimation of amikacin exposure in critically ill septic patients. Ther Drug Monit. 2010;32(6):749-756. doi:10.1097/FTD.0b013e3181f675c2"
  vignette <- "Delattre_2010_amikacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 55.5 mL/min (population median, Delattre 2010 Table 1). The effect is applied to CL as a divisive normalization (CRCL / 55.5) inside an additive linear covariate term: TVCL = exp(lcl) + e_crcl_cl * (CRCL / 55.5).",
      source_name        = "CLCR"
    )
  )

  population <- list(
    n_subjects     = 88L,
    n_studies      = 1L,
    age_range      = "22-89 years",
    age_median     = "65 years",
    weight_range   = "38-125 kg",
    weight_median  = "70 kg",
    sex_female_pct = 35,
    race_ethnicity = "Not reported (Belgian university-hospital ICU population)",
    disease_state  = "Severe sepsis or septic shock (29% septic shock, 71% severe sepsis); ICU patients during first 24 hours of antibiotic treatment",
    dose_range     = "25 mg/kg amikacin IV infusion over 30 minutes (single first dose), combined with broad-spectrum beta-lactam (piperacillin, ceftazidime, cefepime, or meropenem)",
    regions        = "Belgium (4 university-hospital ICUs: Cliniques universitaires St-Luc, Hopital Erasme, Universitair Ziekenhuis Brussel, Clinique St-Pierre)",
    apache_ii      = "median 20 (range 6-45)",
    sofa           = "median 8 (range 1-19)",
    renal_function = "Cockcroft-Gault creatinine clearance median 55.5 mL/min (range 12.3-408.3); raw mL/min, not BSA-normalized",
    notes          = "Baseline demographics per Delattre 2010 Table 1. 88 adults enrolled over 20 months in 4 Belgian ICUs. Mechanical ventilation 52%, catecholamine administration 53%. Pregnant / breastfeeding women, dialysis patients, and those with prior allergic reaction to study drugs were excluded. The study captured the early phase of the septic process (first 24 h of treatment), distinguishing it from prior aminoglycoside ICU PK analyses conducted at steady state."
  )

  ini({
    # Structural parameters at the typical septic patient (median CRCL 55.5 mL/min);
    # Delattre 2010 Table 2 final-model column.
    lcl <- log(0.77); label("Non-renal / baseline component of CL (intercept) (L/h)") # Delattre 2010 Table 2: CL = 0.77 L/h (intercept of additive linear CL ~ CRCL model)
    lvc <- log(19.2); label("Central volume of distribution V1 (L)")                  # Delattre 2010 Table 2: V1 = 19.2 L
    lvp <- log(9.38); label("Peripheral volume of distribution V2 (L)")               # Delattre 2010 Table 2: V2 = 9.38 L
    lq  <- log(4.38); label("Intercompartmental clearance Q (L/h)")                   # Delattre 2010 Table 2: Q = 4.38 L/h

    # Covariate effect: additive linear CL ~ CRCL with divisive normalization to
    # the population median CRCL (55.5 mL/min). Per the operator-resolved sidecar
    # (request-001 / response-001), the paper's "covariates centered to their
    # median" wording is interpreted as divisive normalization (cov = CRCL / 55.5)
    # because subtractive centering yields physically absurd CL values across the
    # observed CRCL range whereas divisive normalization reconciles with the
    # structural-model CL of 2.21 L/h at typical patient (0.77 + 1.42 = 2.19 L/h).
    e_crcl_cl <- 1.42; label("Renal CL slope per (CRCL / 55.5) (L/h)") # Delattre 2010 Table 2: theta_CL-CLCR = 1.42 (slope of additive linear CL ~ CRCL covariate term)

    # Inter-individual variability (Delattre 2010 Table 2 final-model CV%);
    # omega^2 = log(CV^2 + 1) for log-normal etas.
    etalcl ~ 0.30045 # log(0.592^2 + 1); 59.2% CV on CL
    etalvc ~ 0.14226 # log(0.391^2 + 1); 39.1% CV on V1
    etalvp ~ 0.17403 # log(0.436^2 + 1); 43.6% CV on V2
    etalq  ~ 0.02751 # log(0.167^2 + 1); 16.7% CV on Q

    # Combined residual error (Delattre 2010 Table 2 final-model row).
    propSd <- 0.268; label("Proportional residual error (fraction)") # Delattre 2010 Table 2: proportional error CV = 26.8%
    addSd  <- 1.03;  label("Additive residual error (mg/L)")         # Delattre 2010 Table 2: additive error SD = 1.03 mg/L
  })
  model({
    # Individual PK parameters. CL has an additive linear covariate term on CRCL
    # with divisive normalization (TVCL = intercept + slope * CRCL / 55.5), wrapped
    # in exp(etalcl) to give log-normal IIV.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL / 55.5)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
