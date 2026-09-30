MerinoBohorquez_2021_piperacillin <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous piperacillin",
    "(given as piperacillin-tazobactam 4/0.5 g) in hospitalized,",
    "non-critically ill adults with Enterobacteriaceae bloodstream infection",
    "in Seville, Spain. Fitted non-parametrically with the NPAG algorithm in",
    "Pmetrics. Clearance is ADDITIVE in a covariate-free intercept and an arm",
    "linear in Cockcroft-Gault creatinine clearance expressed in L/h (paper:",
    "CL = Intercept + Slope x ClCr); the central volume carries no covariate.",
    "Every parameter carries inter-individual variability reconstructed as a",
    "log-normal from the mean and SD of the NPAG support-point distribution.",
    "Residual error is the published Pmetrics assay-error polynomial",
    "SD = 0.4388 + 0.027 x C (combined1), carried with the unreported gamma",
    "multiplier at 1. Merino-Bohorquez 2021, n = 27 subjects, 102 samples.",
    sep = " "
  )
  reference <- paste(
    "Merino-Bohorquez V, Docobo-Perez F, Valiente-Mendez A, Delgado-Valverde M,",
    "Camean M, Hope WW, Pascual A, Rodriguez-Bano J. Population",
    "Pharmacokinetics of Piperacillin in Non-Critically Ill Patients with",
    "Bacteremia Caused by Enterobacteriaceae. Antibiotics (Basel).",
    "2021;10(4):348. doi:10.3390/antibiotics10040348. PMCID: PMC8064303.",
    "Structural and variability estimates are Table 2; the structural",
    "equation is Equation (1); the residual-error polynomial is Methods 4.4.",
    "No supplementary material beyond the article figures was deposited.",
    sep = " "
  )
  vignette <- "MerinoBohorquez_2021_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The dose is the piperacillin component of piperacillin-tazobactam
  # (4 g of the 4/0.5 g vial) in mg, and the HPLC assay measured TOTAL
  # piperacillin in serum (Methods 4.3), so Cc is a total concentration.
  # The paper applies an unbound fraction of 0.7 only inside its target-
  # attainment simulations (Methods 4.5); that factor is not part of the
  # fitted model and lives in the vignette.
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance estimated with the Cockcroft-Gault equation,",
        "raw mL/min, NOT body-surface-area normalized"
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "ASSAY FORM. Methods 4.1: 'creatinine clearance was calculated daily",
        "using the Cockcroft-Gault equation'; Table 1 reports 'CrCl in mL/min,",
        "median (range) 50.7 (45.3-255.3)'. Supply a Cockcroft-Gault value in",
        "mL/min; the model converts it to L/h internally (x 60 / 1000).",
        "SCALE ON WHICH IT ENTERS -- L/h. Table 2 prints the covariate equation",
        "as 'CL = Intercept + slope x creatinine clearance (L/h)' and Equation",
        "(1) as dX/dt = R(1) - ((Intercept + Slope x ClCr)/Vc) x X1. Reading",
        "ClCr in L/h gives a typical CL of 4.556 + 1.353 x 3.04 = 8.67 L/h at the",
        "median 50.7 mL/min, consistent with the observed steady-state",
        "concentrations of Figure 1 (up to ~250 mg/L after 4 g q8h) and with the",
        "Table 3 target-attainment grid (reproduced in the vignette). Reading",
        "ClCr in mL/min would give 73 L/h, which is incompatible with both.",
        "TIME-FIXED in practice: all samples were drawn within one steady-state",
        "dosing interval (Methods 4.3), so one value per subject applies.",
        "RANGE CAVEAT: Table 1's minimum of 45.3 mL/min conflicts with Results",
        "2.1, which says three patients had CrCl < 20 mL/min/1.73 m^2 and",
        "received q12h dosing; the Discussion adds that 'few' patients had",
        "severe renal impairment. The model was nevertheless used by the",
        "authors for simulations down to kidney failure (< 15 mL/min/1.73 m^2)."
      ),
      source_name = "ClCr"
    )
  )

  # Covariates the paper screened but did not retain (Methods 4.4:
  # 'Creatinine clearance, weight, age, sex, and body mass index (BMI) were
  # explored as covariates'; Results 2.2: 'No other covariates improved the
  # final model').
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 4.4) and not retained. Cohort weight distribution not reported."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods 4.4) and not retained. Table 1: median 76.5 years (range 48-86)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods 4.4) and not retained. Table 1: 17 of 27 (62.96%) male."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened (Methods 4.4) and not retained. Table 1: BMI >= 25 in 19 patients (79.1%)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 27L,
    n_studies = 1L,
    n_samples = 102L,
    age_range = "48-86 years",
    age_median = "76.5 years",
    sex_female_pct = 37.04,
    disease_state = paste(
      "Hospitalized, non-critically ill (not admitted to intensive care)",
      "adults with monomicrobial Enterobacteriaceae bloodstream infection",
      "treated with piperacillin-tazobactam monotherapy started within 12 h of",
      "blood culture. Sources: urinary tract 18 (66.7%), biliary tract 7",
      "(25.9%), other intra-abdominal 2 (7.4%). Organisms: E. coli 15,",
      "K. oxytoca 6, K. pneumoniae 3, E. aerogenes 2, E. cloacae 1; 3 ESBL",
      "producers. Median Pitt score 2 (0-5), SOFA 3 (0-8), Charlson 2.5 (0-8).",
      "Neutropenic patients were excluded."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance median 50.7 mL/min (Table 1 range",
      "45.3-255.3, but Results 2.1 reports three patients with CrCl < 20",
      "mL/min/1.73 m^2 who received q12h dosing)."
    ),
    dose_range = paste(
      "Piperacillin-tazobactam 4/0.5 g by 4 h extended infusion q8h, either",
      "alone or preceded by a 4/0.5 g 30 min loading infusion (local protocol",
      "changed during the study); 4/0.5 g q12h in patients with CrCl < 20",
      "mL/min/1.73 m^2."
    ),
    regions = "Spain (Hospital Universitario Virgen Macarena, Seville), October 2012 - February 2015",
    notes = paste(
      "Prospective single-centre study. Four serum samples per patient at",
      "steady state, 1, 4, 6 and 8 h after the start of an infusion within a",
      "single dosing interval; none below the 1 mg/L LLOQ. Total piperacillin",
      "by HPLC (range 1-500 mg/L). Demographics from Table 1."
    )
  )

  ini({
    # =====================================================================
    # STRUCTURAL PARAMETERS -- Table 2, 'Mean' column of the Pmetrics NPAG
    # support-point distribution. Table 2 also prints an SD column (the
    # spread of that distribution across subjects, i.e. inter-individual
    # variability, not uncertainty) and a Median column (noted per line).
    # Following the other Pmetrics/NPAG models in this library
    # (Sime_2019_tazobactam.R, Tsai_2023_ceftriaxone.R), each MEAN is carried
    # as the typical value, i.e. the median of the log-normal marginal.
    # =====================================================================

    # --- Clearance: CL = Intercept + Slope x ClCr (L/h) -------------------
    # Equation (1) and Table 2. Mapped onto the registered additive
    # clearance-arm canonicals lcl_nonren (intercept: clearance at zero
    # creatinine clearance) and lcl_renal (slope arm), as in
    # Sime_2019_tazobactam.R.
    lcl_nonren <- log(4.556)
    label("Clearance intercept, the arm independent of creatinine clearance (L/h)")
    # Table 2 'Intercept (L/h)': mean 4.556, SD 5.035, median 3.503
    lcl_renal <- log(1.353)
    label("Slope of clearance on Cockcroft-Gault creatinine clearance, both in L/h (unitless)")
    # Table 2 'Slope': mean 1.353, SD 1.032, median 1.39

    # --- Volume -----------------------------------------------------------
    lvc <- log(30.68)
    label("Central volume of distribution Vc (L)")
    # Table 2 'Volume of distribution, Vc (L)': mean 30.68, SD 23.349, median 20.039

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY -- Table 2 SD column, as CV = SD/mean,
    # converted with omega^2 = log(CV^2 + 1). A parametric approximation of
    # the non-parametric NPAG density; the marginals are encoded as
    # independent because the joint (support-point) distribution is not
    # published.
    #   Intercept: 5.035 / 4.556 = 1.10514 -> log(1.10514^2 + 1) = 0.798104
    #   Slope:     1.032 / 1.353 = 0.76275 -> log(0.76275^2 + 1) = 0.458555
    #   Vc:       23.349 / 30.68 = 0.76105 -> log(0.76105^2 + 1) = 0.456916
    # =====================================================================
    etalcl_nonren ~ 0.798104 # Table 2 Intercept: SD 5.035 / mean 4.556
    etalcl_renal ~ 0.458555 # Table 2 Slope: SD 1.032 / mean 1.353
    etalvc ~ 0.456916 # Table 2 Vc: SD 23.349 / mean 30.68

    # =====================================================================
    # RESIDUAL ERROR -- Methods 4.4: 'The data were weighted by the inverse
    # of the estimated assay variance ... given by SD (mg/L) = gamma x
    # (0.4388 + 0.027 x C)'. The assay polynomial coefficients are fixed
    # inputs to Pmetrics (derived from quality-control samples), so they are
    # wrapped in fixed(). The multiplier gamma was estimated but its final
    # value is not reported anywhere in the paper; it is carried here as 1
    # (the unscaled assay polynomial), as in Tsai_2023_ceftriaxone.R. The
    # Pmetrics polynomial ADDS the two terms linearly, hence combined1().
    # =====================================================================
    addSd <- fixed(0.4388)
    label("Additive residual SD, assay polynomial C0 (mg/L)")
    # Methods 4.4: C0 = 0.4388 mg/L (gamma unreported, taken as 1)
    propSd <- fixed(0.027)
    label("Proportional residual SD, assay polynomial C1 (fraction)")
    # Methods 4.4: C1 = 0.027 (gamma unreported, taken as 1)
  })

  model({
    # Cockcroft-Gault creatinine clearance, mL/min -> L/h (the scale on which
    # Table 2 / Equation (1) apply the slope).
    crcl_lh <- CRCL * 60 / 1000

    # Individual parameters. Each clearance arm carries its own eta because
    # Table 2 reports a separate SD for the intercept and for the slope.
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl_renal <- exp(lcl_renal + etalcl_renal) * crcl_lh
    cl <- cl_nonren + cl_renal
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # Equation (1): dX1/dt = R(1) - ((Intercept + Slope x ClCr)/Vc) x X1.
    # The intravenous infusion R(1) enters central via the event table.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
