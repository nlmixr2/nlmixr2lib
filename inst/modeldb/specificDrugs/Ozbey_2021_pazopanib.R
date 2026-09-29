Ozbey_2021_pazopanib <- function() {
  description <- paste(
    "One-compartment population PK model for oral pazopanib in adult and paediatric cancer patients",
    "(therapeutic-drug-monitoring cohort with soft-tissue or Ewing sarcoma, pooled with a phase I/II",
    "pazopanib + temozolomide glioblastoma trial), parameterised on apparent clearance and volume with",
    "first-order absorption and linear elimination. Aspartate aminotransferase enters V/F as a power",
    "model centred at the cohort median of 36.5 U/L. Exponential inter-individual variability on ka,",
    "V/F and CL/F, inter-occasion variability on V/F and CL/F, and a combined additive-plus-proportional",
    "residual error. Developed to estimate an individual AUC from a single randomly timed sample and to",
    "define a target AUC of 750 mg*h/L equivalent to the 20.5 mg/L trough target.",
    sep = " "
  )
  reference <- paste(
    "Ozbey AC, Combarel D, Poinsignon V, Lovera C, Saada E, Mir O, Paci A.",
    "Population Pharmacokinetic Analysis of Pazopanib in Patients and Determination of Target AUC.",
    "Pharmaceuticals (Basel). 2021;14(9):927. doi:10.3390/ph14090927.",
    sep = " "
  )
  vignette <- "Ozbey_2021_pazopanib"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "pazopanib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pazopanib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AST = list(
      description = "Serum aspartate aminotransferase (ASAT)",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Ozbey 2021 Table 4: min 22.0, Q1 35.4, median 36.5, Q3 37.0, max 233 UI/L. Power model on V/F",
        "centred at the median 36.5 UI/L (Equation 4: log(V/Fi) = log(Vpop/F) + beta_V/F *",
        "log(ASATi / 36.5) + eta_Vi), with beta_V/F = -0.838 (Table 1). It was the only covariate",
        "retained by the AIC-based stepwise search, and lowered the V/F IIV from 39.8% to 24.8%",
        "(Results). Missing values were imputed with the patient's own median, else the population",
        "median (Methods 4.1); the paper does not say whether ASAT was handled as time-varying.",
        sep = " "
      ),
      source_name = "ASAT"
    ),
    OCC = list(
      description = "Occasion index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Decomposed inside model() into indicators oc1..oc6 that select the per-occasion IOV etas on",
        "CL/F and V/F (Methods Equations 1-2, eta_Vocc and eta_Clocc). Ozbey 2021 does not state how",
        "occasions were delimited in the Monolix dataset. The trial patients had full profiles on",
        "day 1 of cycles 1, 2 and 4 (three occasions); the monitoring patients had 1-6 samples each.",
        "Six occasions are encoded, which covers the largest per-patient sample count. Records with OCC",
        "outside 1..6 carry no IOV. For a single-occasion simulation pass OCC = 1. See the vignette",
        "Assumptions and deviations.",
        sep = " "
      ),
      source_name = "OCC"
    )
  )

  # Covariates Ozbey 2021 screened in the stepwise search but did not retain.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Table 4: 3.00-87.0 years, median 50.5. Tested on all parameters, not retained (Methods 4.3)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Table 4: 44 male / 29 female. Tested on all parameters, not retained (Methods 4.3)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Table 4: 44.0-180 umol/L, median 65.6. Tested on all parameters, not retained (Methods 4.3)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Table 4: 34.0-63.0 g/L, median 39.5. Tested on all parameters, not retained (Methods 4.3)."
    ),
    ALT = list(
      description = "Serum alanine aminotransferase (ALAT)",
      units = "U/L",
      type = "continuous",
      notes = "Table 4: 12.0-370 UI/L, median 35.0. Tested on all parameters, not retained (Methods 4.3)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 73L,
    n_studies = 2L,
    n_concentrations = 406L,
    age_range = "3-87 years (3 children)",
    age_median = "50.5 years",
    weight_range = "Not reported",
    sex_female_pct = 39.7,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "58 patients from routine pazopanib therapeutic drug monitoring at Gustave Roussy (2012-2018;",
      "55 adults and 3 children treated for soft-tissue sarcoma and Ewing sarcoma), pooled with 15 adults",
      "with glioblastoma from a phase I/II trial of pazopanib plus temozolomide (NCT02331498).",
      sep = " "
    ),
    hepatic_function = "ASAT median 36.5 UI/L (range 22.0-233); ALAT median 35.0 UI/L (range 12.0-370) (Table 4).",
    dose_range = paste(
      "Oral pazopanib 200-800 mg once daily. 52 patients were dosed fasted and 15 with food; the fed",
      "state was unknown for 6. Food was not tested as a covariate.",
      sep = " "
    ),
    regions = "France (Gustave Roussy, Villejuif; Centre Antoine Lacassagne, Nice)",
    notes = paste(
      "406 concentrations: 126 monitoring samples (1-6 per patient, mean 2.2, drawn at a mean 24.6 h",
      "post-dose) and 280 samples from 36 full profiles (0, 0.5, 1, 2, 4, 6, 8, 24 h on day 1 of cycles",
      "1, 2 and 4). HPLC-UV assay, LLOQ 1 mg/L. Fitted in Monolix 2018R1 by SAEM; evaluated by",
      "goodness-of-fit, VPC (1000 replicates), NPDE, 100-sample bootstrap and leave-one-out",
      "cross-validation in 10 steady-state trial patients (Methods 4.1-4.3).",
      sep = " "
    )
  )

  ini({
    # Structural parameters (Table 1, 'Value (RSE%)'); typical values at the
    # median ASAT of 36.5 UI/L.
    lka <- log(0.976); label("Absorption rate constant (1/h)")
    # Table 1 row 'ka (h-1)' 0.976 (RSE 12.3%); Results text 'ka = 0.976 h-1'.
    lvc <- log(22.3); label("Apparent volume of distribution V/F at ASAT 36.5 U/L (L)")
    # Table 1 row 'V/F (L)' 22.3 (RSE 9.25%); Results text 'V/F = 22.3 L'.
    lcl <- log(0.458); label("Apparent clearance CL/F (L/h)")
    # Table 1 row 'Cl/F (L h-1)' 0.458 (RSE 9.73%); Results text 'Cl/F = 0.458 L h-1'.

    # Covariate effect: power model on V/F (Methods Equation 4).
    e_ast_vc <- -0.838; label("Power exponent on (AST / 36.5 U/L) for V/F (unitless)")
    # Table 1 'Covariate (RSE%)' column 'ASAT: -0.838 (29)'; Methods text 'beta V/F = -0.838 (29%)'.

    # Inter-individual variability. The Table 1 footnote says IIV is
    # 'expressed in omega^2', but the printed values are Monolix omega (SD)
    # values: the Results text reads the V/F entry 0.248 directly as '24.8%',
    # and the Figure 4 parameter-distribution curves only reproduce with the
    # SD reading (the ka density peaks at about 2.0 per h^-1; the log-normal
    # modal density exp(s^2 / 2) / (0.976 * s * sqrt(2 * pi)) is 1.98 at
    # s = 0.211 and 0.99 at s = sqrt(0.211)). Variances below are therefore
    # the squared Table 1 values.
    etalka ~ 0.044521
    # Table 1 row 'ka' IIV 0.211 (RSE 44.6%), SD scale; 0.211^2 = 0.044521.
    etalvc ~ 0.061504
    # Table 1 row 'V/F' IIV 0.248 (RSE 46.9%), SD scale; 0.248^2 = 0.061504.
    etalcl ~ 0.509796
    # Table 1 row 'Cl/F' IIV 0.714 (RSE 13.3%), SD scale; 0.714^2 = 0.509796.

    # Inter-occasion variability on V/F and CL/F (Methods Equations 1-2),
    # same SD scale as the IIV. Occasion-indicator expansion; occasions after
    # the first share the occasion-1 variance and are fixed(). The number of
    # occasions is not reported; six are provided (see the OCC notes).
    etaiov_vc_1 ~ 0.147456
    etaiov_vc_2 ~ fixed(0.147456)
    etaiov_vc_3 ~ fixed(0.147456)
    etaiov_vc_4 ~ fixed(0.147456)
    etaiov_vc_5 ~ fixed(0.147456)
    etaiov_vc_6 ~ fixed(0.147456)
    # Table 1 row 'V/F' IOV 0.384 (RSE 18%), SD scale; 0.384^2 = 0.147456.
    etaiov_cl_1 ~ 0.137641
    etaiov_cl_2 ~ fixed(0.137641)
    etaiov_cl_3 ~ fixed(0.137641)
    etaiov_cl_4 ~ fixed(0.137641)
    etaiov_cl_5 ~ fixed(0.137641)
    etaiov_cl_6 ~ fixed(0.137641)
    # Table 1 row 'Cl/F' IOV 0.371 (RSE 12.4%), SD scale; 0.371^2 = 0.137641.

    # Residual error: 'The residual variability is described using a combined
    # model' (Methods 4.3). Encoded as the Monolix 2018 default combined form
    # (combined1, SD = a + b * f).
    addSd <- 4.9; label("Additive residual error a (mg/L)")
    # Table 1 row 'Constant residual error' 4.9 (RSE 8.07%).
    propSd <- 0.05; label("Proportional residual error b (fraction)")
    # Table 1 row 'Proportional residual error' 0.05 (RSE 26.2%).
  })

  model({
    # Inter-occasion variability terms, selected by the occasion index OCC.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5 + oc6 * etaiov_vc_6
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # Individual parameters (Methods Equations 1, 2 and 4).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc + iov_vc) * (AST / 36.5)^e_ast_vc
    cl <- exp(lcl + etalcl + iov_cl)

    kel <- cl / vc

    # One-compartment model with first-order absorption and linear elimination.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg over V/F in L gives mg/L (the paper's ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
