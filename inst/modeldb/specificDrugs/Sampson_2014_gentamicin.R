Sampson_2014_gentamicin <- function() {
  description <- "One-compartment IV population PK model of gentamicin in term neonates with hypoxic-ischaemic encephalopathy undergoing whole-body hypothermia, as reported (model originally developed by Frymoyer 2013; this model file reproduces the parameter values stated by Sampson 2014 during the model's external predictive-performance evaluation). Allometric body-weight scaling on CL (fixed exponent 0.75) and linear body-weight scaling on V (exponent 1) referenced to a 3.3 kg neonate, with a power effect of serum creatinine on CL (exponent -0.566); inter-individual variability on CL only and proportional residual error."
  reference <- paste(
    "Sampson MR, Frymoyer A, Rattray B, Cotten CM, Smith B, Capparelli E,",
    "Bonifacio SL, Cohen-Wolkowiez M.",
    "Predictive performance of a gentamicin population pharmacokinetic model",
    "in neonates receiving full-body hypothermia.",
    "Ther Drug Monit. 2014;36(5):584-589.",
    "doi:10.1097/FTD.0000000000000056.",
    "Model originally developed in",
    "Frymoyer A, Lee S, Bonifacio SL, Meng L, Lucas SS, Guglielmo BJ, Sun Y,",
    "Verotta D.",
    "Every 36-h gentamicin dosing in neonates with hypoxic-ischemic",
    "encephalopathy receiving hypothermia.",
    "J Perinatol. 2013;33(10):778-782. doi:10.1038/jp.2013.59",
    "(PMID 23553582); the present file reproduces the parameter values stated",
    "in Sampson 2014 Methods (page 585).",
    sep = " "
  )
  vignette <- "Sampson_2014_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (birth weight in the source cohorts).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL (fixed exponent 0.75) and linear scaling on",
        "V (exponent 1), each referenced to a 3.3 kg neonate. The 3.3 kg",
        "reference is the median birth weight in the Frymoyer 2013",
        "model-development cohort and is reused verbatim in Sampson 2014",
        "(page 585). Cohort birth-weight range across the two Sampson 2014",
        "validation datasets: 1.9-4.6 kg (median 3.3)."
      ),
      source_name        = "BW"
    ),
    CREAT = list(
      description        = "Serum creatinine.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL with exponent -0.566 in the form",
        "CL = TVCL * (BW/3.3)^0.75 * (1/SCR)^0.566 (Sampson 2014 page 585;",
        "Frymoyer 2013 final model). The equation has no explicit SCR",
        "reference value -- the implicit anchor is 1 mg/dL, so the SCR factor",
        "is identically 1 when CREAT = 1 mg/dL. Time-varying per measurement.",
        "Cohort SCR range across the two Sampson 2014 validation datasets:",
        "0.4-1.9 mg/dL (median 1.0). SCR was measured by the Jaffe method in",
        "the Frymoyer 2013 model-development cohort and the Sampson 2014",
        "Validation A cohort; the Validation B cohort used enzymatic methods",
        "before January 2008 and Jaffe methods afterward (Sampson 2014 page",
        "584; analytical-method differences are discussed as a possible",
        "driver of the under-prediction observed in Validation B)."
      ),
      source_name        = "SCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 29L,
    n_studies      = 1L,
    n_observations = 47L,
    ga_range       = "36-42 weeks (median 40) at birth",
    pna_range      = "1-4 days (median 2) at the first PK sample",
    weight_range   = "approx. 1.9-4.6 kg birth weight (median 3.3); model-development cohort range reported via the Frymoyer 2013 publication and reproduced in Sampson 2014",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Term newborns with hypoxic-ischaemic encephalopathy treated with",
      "whole-body hypothermia (target core temperature 33.5 degC) per the",
      "standard NICHD/ICEH protocol (cooling for 72 h initiated within 6 h",
      "of birth, then 8 h of rewarming to normothermia)."
    ),
    dose_range     = paste(
      "Gentamicin IV given per institutional dosing during therapeutic drug",
      "monitoring. The Frymoyer 2013 model-development cohort received",
      "gentamicin as part of routine care at UCSF; the Sampson 2014",
      "validation cohorts used 5 mg/kg q36h (Validation A, UCSF) and 3.5-4",
      "mg/kg q24h or q36h (Validation B, Duke)."
    ),
    regions        = "United States (UCSF; Frymoyer 2013 model-development cohort).",
    notes          = paste(
      "This entry reproduces the published Frymoyer 2013 one-compartment IV",
      "popPK model as cited verbatim in Sampson 2014 Methods 'Pharmacokinetic",
      "Analysis' (page 585):",
      "CL(L/h) = 0.111 * (BW(kg)/3.3)^0.75 * (1/SCR(mg/dL))^0.566;",
      "V(L) = 1.56 * (BW(kg)/3.3); IIV(CL) = 16.1% CV;",
      "proportional residual error = 16.2%. The Sampson 2014 paper does NOT",
      "re-estimate the model -- it evaluates the model's predictive",
      "performance on two external datasets (Validation A from UCSF,",
      "Validation B from Duke) and reports adequate prediction in",
      "Validation A but under-prediction in Validation B. Estimation in the",
      "original Frymoyer 2013 publication: NONMEM 7.2. Companion library",
      "entries for the same Frymoyer 2013 model (extracted from the",
      "original Frymoyer 2013 paper) may also exist; this file is the",
      "Sampson 2014 representation."
    ),
    validation_cohort = paste(
      "Sampson 2014 evaluated the model on two external retrospective",
      "cohorts: Validation A = 18 UCSF neonates (33 samples), Validation B",
      "= 23 Duke neonates (43 samples, of which 76% within 96 h of birth).",
      "Median AFE 1.1 in Validation A (acceptable) and 0.6 in Validation B",
      "(under-prediction); Validation B was hypothesised to be driven by",
      "differences in serum-creatinine assay (enzymatic vs Jaffe) and by",
      "dataset structure (model built on initial peak + trough only,",
      "Validation B included post-therapy samples)."
    )
  )

  ini({
    # Structural population PK parameters reported in Sampson 2014 Methods
    # 'Pharmacokinetic Analysis' (page 585), where the published Frymoyer
    # 2013 model is reproduced verbatim:
    #   CL(L/h) = 0.111 * (BW(kg)/3.3)^0.75 * (1/SCR(mg/dL))^0.566
    #   V(L)    = 1.56  * (BW(kg)/3.3)
    # Typical values are reported at the cohort reference (BW = 3.3 kg,
    # SCR = 1 mg/dL).
    lcl <- log(0.111); label("Clearance at reference BW = 3.3 kg, SCR = 1 mg/dL (L/h)")  # Sampson 2014 page 585: CL = 0.111 L/h (Frymoyer 2013 final model)
    lvc <- log(1.56);  label("Volume of distribution at reference BW = 3.3 kg (L)")       # Sampson 2014 page 585: V  = 1.56  L  (Frymoyer 2013 final model)

    # Allometric body-weight exponents. Sampson 2014 page 585 quotes the
    # Frymoyer 2013 model using the canonical theoretical-allometry
    # exponent of 0.75 on CL (no SE / RSE reported -- fixed at the West
    # canonical value) and a linear (exponent 1) body-weight effect on V.
    e_wt_cl <- fixed(0.75); label("Allometric body-weight exponent on CL (unitless)")  # Sampson 2014 page 585: (BW/3.3)^0.75 (theoretical-allometry exponent, fixed)
    e_wt_vc <- fixed(1);    label("Linear body-weight exponent on V (unitless)")        # Sampson 2014 page 585: V = 1.56 * (BW/3.3) implies exponent 1 (fixed)

    # Power-form serum-creatinine effect on CL. Sampson 2014 page 585
    # writes the covariate as (1/SCR(mg/dL))^0.566; the implicit reference
    # SCR is 1 mg/dL, so the encoded form is (1 / CREAT)^e_creat_cl with
    # e_creat_cl = 0.566.
    e_creat_cl <- 0.566; label("Power exponent of (1/CREAT_mg_per_dL) on CL (unitless)")  # Sampson 2014 page 585: (1/SCR)^0.566 (Frymoyer 2013 final model)

    # Inter-individual variability on CL only (Sampson 2014 page 585 and
    # Discussion: 'the model only included inter-individual variability
    # for one of two primary PK parameters (CL)'). Reported as 16.1% CV;
    # converted to log-scale variance via omega^2 = log(1 + CV^2)
    # = log(1 + 0.161^2) = 0.025591.
    etalcl ~ 0.025591  # Sampson 2014 page 585: IIV CL = 16.1% CV

    # Residual error (proportional). Reported as 16.2% CV.
    propSd <- 0.162; label("Proportional residual error (fraction)")  # Sampson 2014 page 585: proportional residual error = 16.2%
  })

  model({
    # Individual PK parameters (Sampson 2014 page 585 model equations).
    # CREAT enters as (1/CREAT)^e_creat_cl per the source equation
    # (1/SCR(mg/dL))^0.566.
    cl <- exp(lcl + etalcl) * (WT / 3.3)^e_wt_cl * (1 / CREAT)^e_creat_cl
    vc <- exp(lvc)          * (WT / 3.3)^e_wt_vc

    kel <- cl / vc

    # One-compartment IV PK. Gentamicin is administered as IV doses
    # (typically as a 30-min infusion in the source TDM cohorts); the
    # library model does not hard-code an infusion duration so users
    # supply rate or dur per dose in their event table.
    d/dt(central) <- -kel * central

    # Plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
