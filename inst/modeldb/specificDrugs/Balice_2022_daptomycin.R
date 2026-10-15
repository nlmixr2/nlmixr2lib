Balice_2022_daptomycin <- function() {
  description <- paste(
    "One-compartment population PK model with first-order elimination for",
    "intravenous daptomycin (30-min infusion, 4-12 mg/kg once daily) in",
    "hospitalised adults with severe Gram-positive infections, built from",
    "routine therapeutic drug monitoring peak and trough samples (Pisa,",
    "Italy; DAPTOLIN study). Clearance is additive-linear in Cockcroft-Gault",
    "creatinine clearance centred at 63.35 mL/min with an additive",
    "(normal, L/h) between-subject random effect; the volume of distribution",
    "is additive-linear in female sex centred at the model-building cohort's",
    "female fraction of 0.309. Residual error is piecewise: an",
    "intercept-plus-slope standard deviation (additive + proportional summed",
    "linearly) for individual predictions below 80 mg/L and a purely",
    "proportional standard deviation at or above 80 mg/L."
  )
  reference <- paste(
    "Balice G, Passino C, Bongiorni MG, Segreti L, Russo A, Lastella M,",
    "Luci G, Falcone M, Di Paolo A. Daptomycin Population Pharmacokinetics",
    "in Patients Affected by Severe Gram-Positive Infections: An Update.",
    "Antibiotics (Basel). 2022;11(7):914. doi:10.3390/antibiotics11070914."
  )
  vignette <- "Balice_2022_daptomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The between-subject random effect on clearance is additive on the linear
  # (L/h) scale: Balice 2022 Eq 1 prints 'Cl (l/h) = theta1 + theta5 x
  # (ClCr - 63.35)/100 + eta'. It therefore acts on the linear typical value
  # rather than on the log-scale intercept lcl_int, so it is declared here
  # (precedent: Lin 2021 vancomycin).
  paper_specific_etas <- c("etacl")

  # The residual standard deviation switches form at an individual
  # prediction of 80 mg/L (Results, after Table 3). propSdHigh is the
  # 'Error Slope for Higher iPRED' (theta7) that replaces the
  # additive + proportional pair above the switch.
  paper_specific_residual_sds <- c("propSdHigh")

  compartmentData <- list(
    central = list(analyte = "daptomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source symbol ClCr (CrCl in the text). Methods section 4.1: CrCl was",
        "'calculated through the Cockroft-Gault equation', which yields raw",
        "mL/min. Table 1 heads the row 'CrCl (mL/min/1.73m2)', which conflicts",
        "with the stated estimating equation; the Methods equation is taken as",
        "authoritative and the value is treated as raw Cockcroft-Gault",
        "mL/min (the body-weight convention inside the equation is not",
        "stated). Stored under the canonical CRCL column, which accepts raw",
        "mL/min with the assay form documented per model (precedents",
        "Delattre 2010 amikacin, Wu 2024 daptomycin). Enters CL additively as",
        "theta5 x (CRCL - 63.35) / 100 L/h (Eq 1); 63.35 mL/min is the",
        "centring value printed in Eq 1 (not identified by the paper as a",
        "cohort statistic; the model-building cohort mean is 74.6 +/- 39.7",
        "mL/min, Table 1). Treated as time-fixed per subject.",
        collapse = " "
      ),
      source_name = "ClCr"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column Sex, 'coded as 0 for males and 1 for females'",
        "(Results 2.3), identical to the canonical SEXF coding, so no value",
        "transformation is needed. Enters V additively as",
        "- theta6 x (SEXF - 0.309) L (Eq 2) with theta6 = -2.524, so females",
        "have a LARGER volume (12.67 L) than males (10.15 L). The centring",
        "value 0.309 equals the female fraction of the model-building cohort",
        "(29 of 94, Table 1).",
        collapse = " "
      ),
      source_name = "Sex"
    )
  )

  # Covariates Balice 2022 screened but did not retain (Methods 4.4 GAM /
  # lasso screen; Discussion). Documentation only -- none enters model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Considered for subset balancing and screened on V; the lasso",
        "selected sex rather than weight on V (Discussion: 'we expected to",
        "find an effect of weight, instead of sex, on Vd'). Model-building",
        "cohort 72.6 +/- 10.9 kg (Table 1). Weight still matters for",
        "mg/kg dosing because the dose in mg is weight x mg/kg.",
        collapse = " "
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened on V and CL; no model including albumin improved the OFV",
        "(Discussion). Table 1 reports 3.2 +/- 0.6 under a 'mg/dL' label,",
        "which is only plausible as g/dL (32 g/L).",
        collapse = " "
      )
    ),
    TPRO = list(
      description = "Serum total protein",
      units = "g/L",
      type = "continuous",
      notes = paste(
        "Screened on V and CL; not retained (Discussion). Table 1 reports",
        "6.7 +/- 1.0 under a 'mg/dL' label, which is only plausible as g/dL",
        "(67 g/L).",
        collapse = " "
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Not retained; the Discussion notes age entered a post-hoc GAM screen for CL on the full dataset only. Model-building cohort 65.7 +/- 13.2 years (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 94L,
    n_studies = 1L,
    n_observations = "424 plasma concentrations across all 156 enrolled patients (peak and trough TDM samples)",
    age_range = "adults >= 18 years; model-building cohort mean 65.7 +/- 13.2 years",
    weight_range = "model-building cohort mean 72.6 +/- 10.9 kg",
    sex_female_pct = 30.9,
    race_ethnicity = "not reported (single-centre Italian cohort)",
    disease_state = "Hospitalised adults (medical or surgical ward) with a severe Gram-positive infection receiving daptomycin with routine therapeutic drug monitoring",
    dose_range = "Daptomycin 4-12 mg/kg once daily (mean 6.8 +/- 1.6 mg/kg; 491.2 +/- 115.8 mg/day), usually as a 30-min IV infusion",
    regions = "Italy (Pisa University Hospital)",
    renal_function = "Cockcroft-Gault creatinine clearance 74.6 +/- 39.7 mL/min in the model-building cohort (Table 1)",
    sampling = "Peak 1 h after the start of infusion and trough about 23-23.5 h later; plasma by HPLC-UV; 5 BLQ troughs excluded",
    notes = paste(
      "DAPTOLIN observational retrospective study. 134 patients enrolled in a",
      "first round were split by stratified sampling into a model-building",
      "subset (n = 94, used to estimate the final model in NONMEM 7.2) and an",
      "external-validation subset (n = 40); 22 further patients enrolled in a",
      "second round (156 in total) were used only for exploratory covariate",
      "insights. Demographics above are the model-building subset (Table 1).",
      collapse = " "
    )
  )

  ini({
    # Final model, Balice 2022 Table 3 and Eqs 1-2.
    lcl_int <- log(0.636); label("Intercept of the additive-linear CL equation: typical CL at CRCL = 63.35 mL/min (L/h)") # Table 3 theta1 Cl = 0.636 L/h (RSE 6%)
    lvc <- log(10.925); label("Typical volume of distribution at SEXF = 0.309 (L)") # Table 3 theta2 V = 10.925 L (RSE 4%)

    e_crcl_cl <- 0.109; label("Additive CL change per 100 mL/min of CRCL (L/h)") # Table 3 theta5 kCrCl = 0.109 (RSE 164%); Eq 1
    e_sexf_vc <- -2.524; label("Sex coefficient on V, entering as -theta6 x (SEXF - 0.309) (L)") # Table 3 theta6 kSex = -2.524 (RSE 37%); Eq 2

    # Additive (L/h) between-subject variability on CL (Eq 1 '+ eta').
    # Table 3 'IIV Cl' = 0.027 is the NONMEM OMEGA variance, (L/h)^2
    # (SD 0.164 L/h, about 25% of the typical CL); the abstract's '2.7% IIV'
    # restates the variance as a percentage. The variance reading reproduces
    # the paper's Monte-Carlo Table 4 -- see the vignette.
    etacl ~ 0.027 # Table 3 'IIV Cl' eta = 0.027 (RSE 30%)

    # Residual error (Results after Table 3): SD = theta3 + theta4 x IPRED for
    # IPRED < 80 mg/L and SD = theta7 x IPRED at or above 80 mg/L.
    addSd <- 3.805; label("Additive residual SD component, IPRED < 80 mg/L (mg/L)") # Table 3 theta3 = 3.805 mg/L (RSE 42%)
    propSd <- 0.296; label("Proportional residual SD component, IPRED < 80 mg/L (fraction)") # Table 3 theta4 = 0.296 (RSE 14%)
    propSdHigh <- -0.546; label("Proportional residual SD, IPRED >= 80 mg/L (fraction; sign immaterial)") # Table 3 theta7 'Error Slope for Higher iPRED' = -0.546 (RSE 49%)
  })
  model({
    crcl_ref <- 63.35 # mL/min, Eq 1 centring value
    sexf_ref <- 0.309 # Eq 2 centring value (29/94 female, Table 1)

    # Eq 1: additive-linear CL with an additive random effect (L/h)
    cl <- exp(lcl_int) + e_crcl_cl * (CRCL - crcl_ref) / 100 + etacl
    # Eq 2: V = theta2 - theta6 x (Sex - 0.309)
    vc <- exp(lvc) - e_sexf_vc * (SEXF - sexf_ref)

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L: central / vc is mg/L.
    Cc <- central / vc

    # Piecewise residual SD: intercept-plus-slope below 80 mg/L, purely
    # proportional (theta7, applied as its absolute value) at or above.
    sdCc <- (addSd + propSd * Cc) * (Cc < 80) + abs(propSdHigh) * Cc * (Cc >= 80)
    Cc ~ add(sdCc)
  })
}
