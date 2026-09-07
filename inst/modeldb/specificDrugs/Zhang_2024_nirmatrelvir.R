Zhang_2024_nirmatrelvir <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination (NONMEM ADVAN2 TRANS2) for oral nirmatrelvir",
    "coadministered with ritonavir 100 mg (Paxlovid) in Chinese inpatients",
    "with mild-to-moderate COVID-19 (Zhang 2024; single-centre retrospective",
    "trough-dominated sparse-sampling study, 130 samples from 129 patients,",
    "median age 76 years). Apparent clearance carries a priori allometric",
    "body-weight scaling with the exponent fixed at 0.75 on a 70 kg reference,",
    "and an estimated power effect of creatinine clearance referenced to the",
    "cohort median 52.9 mL/min. Apparent volume of distribution carries",
    "allometric body-weight scaling with the exponent fixed at 1, and both the",
    "apparent volume (39 L) and the absorption rate constant (0.8 1/h) were",
    "fixed to literature values because the trough-dominated sparse sampling",
    "could not support estimating them. Interindividual variability on",
    "apparent clearance only; combined additive plus proportional residual",
    "error with both magnitudes fixed at 10 percent."
  )
  reference <- paste(
    "Zhang R, Fan J, Han L, Mao J, Sun L, Yu Y, Fan W, Xie J, Lin B, Lin N",
    "(2024). Population Pharmacokinetics and Dosing Regimen Analysis of",
    "Nirmatrelvir in Chinese Patients with COVID-19 Infection.",
    "Drug Design, Development and Therapy 18:5515-5525.",
    "doi:10.2147/DDDT.S479561. PMCID PMC11622681.",
    sep = " "
  )
  vignette <- "Zhang_2024_nirmatrelvir"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: the two states follow directly from
  # the paper's ADVAN2 TRANS2 oral one-compartment structure (Results, Model
  # Development), and the assayed matrix is plasma (Methods, Sample
  # Processing: "100 uL plasma sample"; Model Development: "PopPK of
  # nirmatrelvir in plasma concentration-time data"). The Abstract's
  # "serum concentrations" wording is inconsistent with the Methods; the
  # Methods matrix is used.
  compartmentData <- list(
    depot   = list(analyte = "nirmatrelvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nirmatrelvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "A priori allometric scaling, not estimated. Zhang 2024 Methods,",
        "Model Development: 'a priori allometric scaling was employed, with",
        "volume terms multiplied by (WT/70) and clearance terms multiplied by",
        "(WT/70)^0.75'. Restated in the Discussion: 'allometric scaling in",
        "CL/F ((WT/70)^0.75) and V/F (WT/70)'. The final-model equation in",
        "Results, Model Development carries the same 70 kg reference:",
        "CL/F = 3.41 x (WT/70)^0.75 x (CRCL/52.9)^theta x exp(eta).",
        "Cohort weight 61.2 +/- 9.3 kg, median 61.1, range 37.5-96.0 kg",
        "(Table 1); the 70 kg reference is therefore above the cohort median.",
        "Table 4 simulates 40, 65, 90 and 115 kg."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CrCl. The only covariate retained after stepwise",
        "forward inclusion and backward elimination (Results, Model",
        "Development: 'the model retained only CrCl (objective function",
        "value = 37.275; chi2 test, df = 1, p < 0.001) as a covariate for",
        "CL/F'). Entered as a centred power term on CL/F with the reference",
        "at the cohort median: CL/F = 3.41 x (WT/70)^0.75 x",
        "(CRCL/52.9)^0.429 x exp(eta). Table 1 reports Creatinine Clearance",
        "in mL/min (56.7 +/- 33.1; median 52.9, range 4.8-289.2), so the",
        "value is the raw un-normalized variant of the canonical CRCL column",
        "-- see the CRCL entry in inst/references/covariate-columns.md and",
        "the Delattre_2010_amikacin.R / Chen_2023_nemonoxacin.R precedents.",
        "The paper does not state which equation was used to estimate CrCl",
        "(no Cockcroft-Gault / MDRD / CKD-EPI attribution appears anywhere in",
        "the Methods, Results or table footnotes), so the estimating equation",
        "is unknown; the raw mL/min units are the only assay fact reported.",
        "Because the model also carries an allometric weight term, supplying",
        "a BSA-normalized value here would silently rescale the renal term.",
        "Tables 3 simulates CrCl 15, 45, 70 and 100 mL/min."
      ),
      source_name        = "CrCl"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 129L,
    n_studies      = 1L,
    n_observations = 130L,
    age_range      = "18.0-97.0 years",
    age_median     = "76.0 years",
    age_mean       = "73.2 years (SD 14.7)",
    weight_range   = "37.5-96.0 kg",
    weight_median  = "61.1 kg",
    weight_mean    = "61.2 kg (SD 9.3)",
    sex_female_pct = 49.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "mild to moderate COVID-19 infection",
    renal_function = paste(
      "Creatinine clearance 56.7 +/- 33.1 mL/min, median 52.9,",
      "range 4.8-289.2 mL/min (Table 1); the cohort spans normal renal",
      "function through severe renal impairment."
    ),
    co_medication  = paste(
      "Ritonavir 100 mg twice daily as the fixed nirmatrelvir/ritonavir",
      "combination. Patients who had used cytochrome P450 or P-glycoprotein",
      "inhibitors or inducers within the prior week were excluded. Other",
      "co-medication was not collected (Results, Demographic Parameters:",
      "'Information related to N/R co-administration not collected')."
    ),
    dose_range     = "300 mg nirmatrelvir / 100 mg ritonavir orally twice daily for 5 days",
    regions        = "China (Changxing People's Hospital, Zhejiang)",
    sampling       = paste(
      "Sparse sampling; the sampling times 'mainly represent valley",
      "concentrations' (Results, Model Development), which is why V/F and Ka",
      "were fixed rather than estimated."
    ),
    notes          = paste(
      "Single-centre retrospective PK trial run December 2022 to June 2023.",
      "129 patients (65 male, 64 female) contributing 130 plasma samples.",
      "Baseline demographics in Table 1. Adults aged 18 years or older with",
      "mild-to-moderate COVID-19; hepatitis B / hepatitis C / HIV positivity",
      "were exclusion criteria."
    )
  )

  ini({
    # Structural parameters, referenced to 70 kg body weight and
    # 52.9 mL/min creatinine clearance.
    lcl <- log(3.41)
    label("Apparent clearance CL/F (L/h) at 70 kg and CrCl 52.9 mL/min")           # Table 2, "CL/F, L/h" = 3.41 (RSE 4.00%; bootstrap median 3.36, 2.5-97.5th 3.16-3.92)

    lvc <- fixed(log(39))
    label("Apparent central volume of distribution V/F (L) at 70 kg")              # Table 2, "V/F, L" = 39.00 (fixed); Results, Model Development: "The V/F and Ka values were set at 39 L and 0.8 h-1"

    lka <- fixed(log(0.8))
    label("First-order absorption rate constant (1/h)")                            # Table 2, "Ka, h-1" = 0.80 (fixed); Results, Model Development: "The V/F and Ka values were set at 39 L and 0.8 h-1"

    # Allometric exponents, fixed a priori rather than estimated.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")                 # Methods, Model Development: "clearance terms multiplied by (WT/70)^0.75"; Results equation "CL/F = 3.41 x (WT/70)^0.75 x (CRCL/52.9)^theta x exp(eta)"

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V/F (unitless)")                  # Methods, Model Development: "volume terms multiplied by (WT/70)"; Discussion: "allometric scaling in CL/F ((WT/70)^0.75) and V/F (WT/70)"

    # Covariate effect.
    e_crcl_cl <- 0.429
    label("Power exponent of creatinine clearance on CL/F, referenced to 52.9 mL/min (unitless)")  # Table 2, "Effect of CrCl on CL/F" = 0.429 (RSE 28.00%; bootstrap median 0.40, 2.5-97.5th 0.21-0.81); reference 52.9 from the Results equation "(CRCL/52.9)^theta"

    # Interindividual variability. Table 2 reports IIV on CL/F under the
    # heading "Interindividual variability, CV%" as 52.30 (RSE 7.00%;
    # bootstrap median 52.90, 2.5-97.5th 44.70-60.80), restated in the text
    # as "reduced the IIV from 58.9% to 52.3%". Converted to a log-normal
    # variance the standard way: omega^2 = log(CV^2 + 1)
    #                                    = log(0.523^2 + 1) = 0.241792.
    etalcl ~ 0.241792                                                              # Table 2, "Interindividual variability, CV%" / "CL/F, L/h" = 52.30

    # Residual error. Results, Model Development: "Residual errors were
    # described by additive plus constant coefficient of variation models,
    # with values set at 10%". Table 2 lists epsilon 1 and epsilon 2 both as
    # "10.00 (fixed)" under the heading "Residual variability, CV%".
    propSd <- fixed(0.1)
    label("Proportional residual error (fraction)")                                # Table 2, "Residual variability, CV%" epsilon 1 = 10.00 (fixed)

    addSd <- fixed(0.1)
    label("Additive residual error (ug/mL)")                                       # Table 2, "Residual variability, CV%" epsilon 2 = 10.00 (fixed); read as sigma = 0.1 in concentration units -- see the vignette Assumptions and deviations section
  })

  model({
    # Individual parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CRCL / 52.9)^e_crcl_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    # Micro-constant.
    kel <- cl / vc

    # ODE system (NONMEM ADVAN2 TRANS2).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation and error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
