Jeon_2014_piperacillin <- function() {
  description <- "Two-compartment IV population PK model for piperacillin in 50 Korean adult burn-ICU patients receiving piperacillin-tazobactam 4.5 g (4 g piperacillin + 0.5 g tazobactam) every 8 h as a 30-min infusion (Jeon 2014)"
  reference   <- "Jeon S, Han S, Lee J, Hong T, Paek J, Woo H, Yim DS. Population pharmacokinetic analysis of piperacillin in burn patients. Antimicrob Agents Chemother. 2014;58(7):3744-3751. doi:10.1128/AAC.02089-13"
  vignette    <- "Jeon_2014_piperacillin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 132 mL/min (cohort mean per Jeon 2014 Table 1; rounded from 132.1 mL/min in the paper text). The effect is applied to CL as a multiplicative renal scaling on the structural intercept: TVCL = exp(lcl) * (CRCL / 132) + e_dai_cl * DAY_AFTER_INJURY.",
      source_name        = "CLCR"
    ),
    DIS_SEPSIS = list(
      description        = "Active sepsis indicator at the PK-sampling window",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column 'sepsis'. 1 = active sepsis at PK sampling (12 of 50 patients per Table 1), 0 = no sepsis. Time-fixed at the PK-sampling window in this paper. Jeon 2014 Methods do not name a formal sepsis diagnostic-criteria reference; the indicator reflects the attending clinician's diagnosis. Used as an additive shift on central volume: TVV1 = exp(lvc) + e_sepsis_vc * DIS_SEPSIS, so septic patients have V1 = 25.3 + 14.8 = 40.1 L versus 25.3 L in non-septic patients (capillary leakage / interstitial edema).",
      source_name        = "sepsis"
    ),
    DAY_AFTER_INJURY = list(
      description        = "Days elapsed since burn injury at the start of PK sampling",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column DAI. Continuous integer days since the burn injury (cohort mean 12.8 days, range 2-68 per Jeon 2014 Table 1). Paper-specific covariate (NOT a registered canonical in inst/references/covariate-columns.md) per the operator-resolved sidecar response (request-001 / response-001, 2026-06-07): treated as a paper-specific column inside this single model only. checkModelConventions() will note the absence of a canonical entry, which is expected. Used as an additive linear effect on CL: TVCL = exp(lcl) * (CRCL / 132) + e_dai_cl * DAY_AFTER_INJURY with e_dai_cl = -0.0874 L/h per day (CL decreases as the hypermetabolic phase resolves). A future continuous-days-since-burn covariate canonical (suggested name POSTBURN_DAYS per the Han 2013 fluconazole sibling-task sidecar) may motivate promoting this concept when a second burn-cohort extraction reuses it."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "20-83 years",
    age_median     = "50.14 years (mean per Table 1)",
    weight_range   = "50-90 kg",
    weight_median  = "66.9 kg (mean per Table 1)",
    sex_female_pct = 20,
    race_ethnicity = "Not reported explicitly; single Korean burn ICU cohort",
    disease_state  = "Adult burn patients in a burn intensive care unit; total body surface area burned mean 34.56% (range 1-81%); serum albumin mean 2.58 g/dL (range 1.6-3.5); 12/50 with active sepsis; 16/50 with clinical edema (puffy face and pitting leg edema); 5/50 on continuous renal replacement therapy",
    dose_range     = "Piperacillin-tazobactam 4.5 g (4 g piperacillin + 0.5 g tazobactam, 8:1 ratio) IV infusion over 30 minutes every 8 h. PK sampling performed at steady state after 5 or more doses",
    regions        = "South Korea (single center: Hangang Sacred Heart Hospital, Hallym University Medical Center, Burn Intensive Care Unit)",
    renal_function = "Cockcroft-Gault creatinine clearance: mean 132.1 mL/min, range 39-231.4 mL/min (raw mL/min, not BSA-normalized)",
    notes          = "Baseline demographics per Jeon 2014 Table 1. 50 adults admitted to the Burn Intensive Care Unit between November 2011 and August 2012. Days after burn injury mean 12.8 (range 2-68). Patients excluded if pregnant, breastfeeding, < 18 years old, or allergic to penicillin."
  )

  ini({
    # Structural parameters at the typical patient (cohort mean CRCL = 132 mL/min,
    # DAY_AFTER_INJURY = 0, no sepsis); Jeon 2014 Table 2 final-model column.
    lcl <- log(16.6);  label("Typical CL at CRCL = 132 mL/min and DAY_AFTER_INJURY = 0 (L/h)") # Jeon 2014 Table 2 theta1 = 16.6 L/h
    lvc <- log(25.3);  label("Typical central volume V1 in non-septic patients (L)")          # Jeon 2014 Table 2 theta2 = 25.3 L
    lvp <- log(16.1);  label("Typical peripheral volume V2 (L)")                              # Jeon 2014 Table 2 theta3 = V2 = 16.1 L
    lq  <- log(0.636); label("Typical intercompartmental clearance Q (L/h)")                  # Jeon 2014 Table 2 theta4 = Q = 0.636 L/h

    # Covariate effects (Jeon 2014 Table 2 final-model column). CRCL enters as
    # multiplicative renal scaling (CRCL/132) without an estimated slope -- the
    # reference 132 mL/min is the cohort mean per Table 1. Sepsis and
    # DAY_AFTER_INJURY enter as additive shifts inside the structural arithmetic.
    e_dai_cl    <- -0.0874; label("Linear DAY_AFTER_INJURY effect on CL (L/h per day)")  # Jeon 2014 Table 2 theta5 = -0.0874 L/h
    e_sepsis_vc <- 14.8;    label("Additive sepsis effect on central volume V1 (L)")     # Jeon 2014 Table 2 theta6 = 14.8 L

    # Inter-individual variability. Paper reports CV; convert to log-normal
    # variance via omega^2 = log(CV^2 + 1). IIV on CL and V1 are correlated
    # (rho = 0.434) per Jeon 2014 Table 2; encoded as a 2x2 block. No IIV on
    # V2 (Jeon 2014 Results: omega V2 not estimated, fixed at 0).
    # var_CL  = log(0.354^2 + 1) = 0.118092
    # var_V1  = log(0.424^2 + 1) = 0.165390
    # cov     = 0.434 * sqrt(0.118092 * 0.165390) = 0.060651
    # var_Q   = log(0.903^2 + 1) = 0.596277
    etalcl + etalvc ~ c(0.118092,
                        0.060651, 0.165390)  # Jeon 2014 Table 2: 35.4 CV CL, 42.4 CV V1, rho 0.434
    etalq           ~ 0.596277               # Jeon 2014 Table 2: 90.3 CV Q

    # Combined residual error (Jeon 2014 Table 2 final-model row).
    propSd <- 0.185; label("Proportional residual error (fraction)")                       # Jeon 2014 Table 2 sigma proportional = 18.5%
    addSd  <- 0.359; label("Additive residual error (mg/L)")                               # Jeon 2014 Table 2 sigma additive = 0.359 mg/L
  })

  model({
    # Individual PK parameters. CL combines multiplicative CRCL scaling on the
    # intercept and an additive DAY_AFTER_INJURY shift, wrapped in exp(etalcl)
    # for log-normal IIV. V1 has an additive sepsis shift wrapped in exp(etalvc).
    # No IIV on V2.
    cl <- (exp(lcl) * (CRCL / 132) + e_dai_cl * DAY_AFTER_INJURY) * exp(etalcl)
    vc <- (exp(lvc) + e_sepsis_vc * DIS_SEPSIS) * exp(etalvc)
    vp <- exp(lvp)
    q  <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
