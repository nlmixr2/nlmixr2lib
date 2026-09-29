Lin_2020_glasdegib_healthy <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption for",
    "oral glasdegib in 49 healthy adult volunteers (studies B1371010 and",
    "B1371014), reported by Lin 2020 as an exploratory comparison with the",
    "cancer-patient model (see Lin_2020_glasdegib). Allometric body weight",
    "scaling (assumed exponents 0.75 on CL/F and Q/F and 1 on Vc/F and",
    "Vp/F, 70 kg reference, as in the patient model). Exponential IIV on",
    "CL/F, Vc/F and ka, with IIV on Vp/F and Q/F fixed at 3.2%;",
    "proportional residual error. The paper's empirical linear food effect",
    "on ka is not quantified in the source and is omitted.",
    sep = " "
  )
  reference <- paste(
    "Lin S, Shaik N, Martinelli G, Wagner AJ, Cortes J, Ruiz-Garcia A.",
    "Population Pharmacokinetics of Glasdegib in Patients With Advanced",
    "Hematologic Malignancies and Solid Tumors.",
    "J Clin Pharmacol. 2020;60(5):605-616. doi:10.1002/jcph.1556.",
    "Healthy-volunteer model: Table 3 'Healthy Volunteer Model' column.",
    sep = " "
  )
  vignette <- "Lin_2020_glasdegib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 3 footnote b says only 'allometric body weight scaling'.",
        "The exponents (0.75 on clearances, 1 on volumes) and the 70 kg",
        "reference are NOT stated for this model; they are assumed equal",
        "to the patient model, where they were fixed at these values."
      ),
      source_name = "BWT"
    )
  )

  covariatesDataExcluded <- list(
    FED = list(
      description = "Fed-state indicator at dosing, 1 = fed, 0 = fasted",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Table 3 footnote b: 'Food was empirically added as a covariate on",
        "ka by linear function.' The coefficient and the reference food",
        "state are not reported, so the effect cannot be encoded; ka here",
        "is the typical value in the (unstated) reference food state."
      )
    )
  )

  compartmentData <- list(
    depot = list(analyte = "glasdegib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "glasdegib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "glasdegib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 2L,
    n_observations = "1937 plasma glasdegib concentrations",
    disease_state = "Healthy adult volunteers (noncancer subjects)",
    dose_range = paste(
      "Oral glasdegib; the only regimen the paper describes is the",
      "B1371010 ketoconazole drug-interaction study (200 mg single dose",
      "alone and on day 4 of ketoconazole 400 mg QD). Doses in B1371014",
      "are not given."
    ),
    notes = paste(
      "Demographics are not reported for this cohort. The paper does not",
      "say whether the ketoconazole-period data were included or how they",
      "were handled; no CYP3A inhibitor covariate is reported. FOCEI",
      "estimation. OFV 83.523 (Table 3)."
    )
  )

  ini({
    # Lin 2020 Table 3 'Healthy Volunteer Model' column: estimate (RSE%).
    lka <- log(0.79); label("First-order absorption rate constant ka, reference food state (1/h)") # Table 3 HV 'ka, hour-1' = 0.79 (8.9%)
    lcl <- log(10.1); label("Apparent clearance CL/F at 70 kg (L/h)") # Table 3 HV 'CL/F, L/h' = 10.1 (4.1%)
    lvc <- log(112); label("Apparent central volume Vc/F at 70 kg (L)") # Table 3 HV 'Vc/F, L' = 112 (4.1%)
    lvp <- log(53.5); label("Apparent peripheral volume Vp/F at 70 kg (L)") # Table 3 HV 'Vp/F, L' = 53.5 (6.9%)
    lq <- log(1.92); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table 3 HV 'Q/F, L/h' = 1.92 (11.4%)

    # Exponents and reference weight not reported for this model; assumed
    # equal to the patient model (Results 'Base Model').
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # assumed; Table 3 footnote b 'allometric body weight scaling'
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # assumed; Table 3 footnote b 'allometric body weight scaling'

    # IIV: CV% converted with omega^2 = log(1 + CV^2).
    etalcl ~ 0.0376965 # Table 3 HV IIV 'CL/F' = 19.6% (5.8%)
    etalvc ~ 0.0120274 # Table 3 HV IIV 'Vc/F' = 11.0% (52.0%)
    etalvp ~ fixed(0.00102348) # Table 3 HV IIV 'Vp/F' = 3.2%, not estimated
    etalq ~ fixed(0.00102348) # Table 3 HV IIV 'Q/F' = 3.2%, not estimated
    etalka ~ 0.25088 # Table 3 HV IIV 'ka' = 53.4% (10.6%)

    propSd <- 0.586; label("Proportional residual error (fraction)") # Table 3 HV 'Residual error, proportional error' = 58.6%; footnote b 'proportional error model of all subjects'
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = ug/mL; x 1000 gives ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
