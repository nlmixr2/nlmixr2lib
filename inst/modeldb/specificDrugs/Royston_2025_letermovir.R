Royston_2025_letermovir <- function() {
  description <- "One-compartment population PK model with first-order absorption and absorption lag time for oral letermovir in adult allogeneic hematopoietic cell transplant recipients, with a cyclosporine effect on apparent oral clearance and inter-occasion variability on CL/F"
  reference <- paste(
    "Royston L, Kunz C, Tonoli D, Lescuyer P, Neofytos D, Gotta V (2025).",
    "Population pharmacokinetic analysis of letermovir in adult hematopoietic",
    "cell transplant recipients. Antimicrob Agents Chemother 69(10):e00697-25.",
    "doi:10.1128/aac.00697-25.",
    "Ka, Tlag and the IIV on Ka and V/F were fixed from the phase III model of",
    "Prohn et al. (2021) CPT Pharmacometrics Syst Pharmacol 10:255-267,",
    "doi:10.1002/psp4.12593, as reproduced in Royston 2025 Table 1.",
    sep = " "
  )
  vignette <- "Royston_2025_letermovir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    CONMED_CSA = list(
      description        = "Concomitant cyclosporine (ciclosporin) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant cyclosporine)",
      notes              = paste(
        "Royston 2025 Table 1: fCSA = -0.38 (RSE 69.5%) on log CL/F.",
        "Six of the 40 patients received cyclosporine and were dosed at 240 mg",
        "once daily rather than 480 mg once daily. The effect is multiplicative",
        "on the log scale: CL/F is multiplied by exp(-0.38) = 0.684, i.e. a 32%",
        "reduction, which is exactly the reduction quoted in the Results text",
        "('reduced by 32% with cyclosporin'). A linear 1 - 0.38 form would give",
        "a 38% reduction and is therefore excluded. Note that the Discussion",
        "separately restates the effect as '-38%' ('estimated in this analysis",
        "to -38% with however large RSE'); that is the raw coefficient, not the",
        "back-transformed effect, and it is the Results' 32% that identifies the",
        "functional form. Monolix, used for this fit, parameterises a",
        "categorical covariate on a log-normal parameter exactly this way.",
        "Cyclosporine inhibits OATP1B1-mediated hepatic uptake of letermovir.",
        sep = " "
      ),
      source_name        = "CSA"
    ),
    OCC = list(
      description        = "Integer-valued occasion index for inter-occasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Royston 2025 Methods: 'each patient visit was treated as an occasion at",
        "steady state'. IOV on CL/F was estimated at SD 0.62 (RSE 7%) on the log",
        "scale (Table 1), which exceeds the IIV of 0.58 and is therefore the",
        "dominant variance component of this model. The paper does not report the",
        "number of visits per patient; values 1..8 are multiplexed inside model()",
        "via binary indicators oc1..oc8, a cap chosen to cover the observed",
        "sampling density (296 observations across 40 patients, ~7.4 per patient).",
        "The cap is an implementation choice, not a paper value. Records with",
        "OCC = 0 or any value outside 1..8 zero every indicator and yield the",
        "typical-value CL/F with no IOV applied. Extending beyond eight occasions",
        "requires adding further etaiov_cl_<n> blocks.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates screened by Royston 2025 in univariable analysis but NOT retained
  # in the final model reported in Table 1. The paper reports the direction of
  # each association in prose only; no regression coefficients are tabulated for
  # any of them, so none can be implemented. Documented here to preserve the
  # provenance of the covariate screen.
  #
  # The list below covers every screened covariate that maps onto an existing
  # canonical column in inst/references/covariate-columns.md. Royston 2025 also
  # screened the following, which have no canonical column and are recorded here
  # in prose rather than as list entries so that no new canonical is implied for
  # an effect that cannot be implemented: vomiting, nausea and mild-to-moderate
  # diarrhea (<10x daily) -- all associated with CL/F in the reported direction
  # only; acute graft-versus-host disease (the register carries only the
  # organ-specific AGVHD_SKIN / AGVHD_LIVER / AGVHD_INTESTINE canonicals, while
  # the paper reports aGvHD generically); posaconazole and pantoprazole
  # coadministration; days since letermovir start; documented infectious
  # complications; and measured cyclosporin blood concentration in the six
  # cyclosporine-treated patients.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Higher weight associated with increased CL/F (allometric scaling) in univariable analysis; no coefficient reported, not retained in the final model.",
      source_name        = "weight"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Younger age associated with increased CL/F in univariable analysis; no coefficient reported, not retained in the final model.",
      source_name        = "age"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Higher serum albumin associated with increased CL/F in univariable analysis (albumin-mediated hepatic uptake of OATP substrates); no coefficient reported, not retained in the final model.",
      source_name        = "albumin"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened; no association with CL/F found. Not retained.",
      source_name        = "gender"
    ),
    CRCL = list(
      description        = "Creatinine clearance / glomerular filtration rate",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as glomerular filtration rate; no association with CL/F found. Not retained.",
      source_name        = "glomerular filtration rate"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened (reported as ALAT); no association with CL/F found. Not retained.",
      source_name        = "ALAT"
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened (reported as ASAT); no association with CL/F found. Not retained.",
      source_name        = "ASAT"
    ),
    CRP = list(
      description        = "C-reactive protein",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened; no association with CL/F found. Not retained.",
      source_name        = "C-reactive protein"
    ),
    WBC = list(
      description        = "White blood cell (leukocyte) count",
      units              = "10^9/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as absolute leukocyte count; no association with CL/F found. Not retained.",
      source_name        = "leukocyte count"
    ),
    NEUT = list(
      description        = "Absolute neutrophil count",
      units              = "cells/mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as absolute neutrophil count; no association with CL/F found. Not retained. The paper reports no unit for this covariate; the canonical unit is given here.",
      source_name        = "absolute neutrophil count"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as one of the liver enzymes/function tests (ALAT, ASAT, and bilirubin); no association with CL/F found. Not retained. The paper reports no unit for this covariate.",
      source_name        = "bilirubin"
    ),
    CONMED_STEROID = list(
      description        = "Systemic corticosteroid (prednisone) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no systemic corticosteroid)",
      notes              = "Screened as prednisone use and associated with DECREASED CL/F in univariable analysis; no coefficient reported, not retained in the final model. The paper's Discussion groups steroid use with acute GvHD as a marker of later time after HCT.",
      source_name        = "prednisone"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "letermovir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species       = "human",
    n_subjects    = 40,
    n_studies     = 1,
    disease_state = "adult CMV-seropositive allogeneic hematopoietic cell transplant recipients receiving primary CMV prophylaxis",
    dose_range    = "480 mg orally once daily, or 240 mg orally once daily when coadministered with cyclosporine",
    regions       = "Switzerland (University Hospitals of Geneva)",
    notes         = paste(
      "Post hoc pharmacometric analysis of a prospective observational therapeutic",
      "drug monitoring study conducted between 1 March 2020 and 20 April 2021",
      "(Royston 2022, Antimicrob Agents Chemother 66:e0065722,",
      "doi:10.1128/aac.00657-22). 296 oral letermovir plasma concentrations were",
      "analysed: 217 trough samples (24 +/- 2 h post-dose) and 79 non-trough",
      "samples from 40 patients, of whom 6 received cyclosporine. The model was",
      "fitted in Monolix 2023 on a time-after-dose scale with each patient visit",
      "treated as an occasion at steady state. This short-form publication does",
      "not tabulate baseline demographics (age, weight, sex distribution); those",
      "are reported in the parent study, Royston 2022. Observed median (IQR)",
      "trough concentrations were 260 (123-518) ng/mL at 480 mg/day and 418",
      "(235-945) ng/mL at 240 mg/day with cyclosporine.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Royston 2025 Table 1, "Adjusted model estimate
    # from one-compartmental model fitted to real-world data" column.

    # Fixed from the Prohn 2021 phase III model; Table 1 marks both "(fixed)".
    lka <- fixed(log(0.15))
    label("Absorption rate constant (1/h)")  # Table 1, Ka row: 0.15 (fixed)
    ltlag <- fixed(log(0.674))
    label("Absorption lag time (h)")  # Table 1, Tlag row: 0.674 (fixed)

    # Estimated on the real-world data.
    lcl <- log(26.5)
    label("Apparent oral clearance CL/F (L/h)")  # Table 1, CL/F row: 26.5 (RSE 12.6%)
    lvc <- log(115.2)
    label("Apparent oral volume of distribution V/F (L)")  # Table 1, V/F row: 115.2 (RSE 12%)

    # Cyclosporine effect on log CL/F. Multiplicative on the log scale:
    # exp(-0.38) = 0.684, i.e. the 32% reduction quoted in the Results text.
    e_csa_cl <- -0.38
    label("Effect of concomitant cyclosporine on log CL/F")  # Table 1, fCSA row: -0.38 (RSE 69.5%)

    # IIV. Table 1 footnote a: variability is "presented as standard deviation
    # of log-transformed parameters", so each variance below is that SD squared.
    etalcl ~ 0.58^2  # Table 1, CL/F row [IIV 0.58 (RSE 17%)]
    etalka ~ fixed(0.72^2)  # Table 1, Ka row [IIV 0.72]; carried from Prohn 2021
    etalvc ~ fixed(0.23^2)  # Table 1, V/F row [IIV 0.23]; carried from Prohn 2021

    # IOV on CL/F, SD 0.62 on the log scale, shared across all occasions.
    # Occasion 1 carries the estimated variance; occasions 2-8 repeat it as a
    # fixed value (the equivalent of NONMEM's $OMEGA BLOCK(1) SAME).
    etaiov_cl_1 ~ 0.62^2  # Table 1, CL/F row [IOV 0.62 (RSE 7%)]
    etaiov_cl_2 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance
    etaiov_cl_3 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance
    etaiov_cl_4 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance
    etaiov_cl_5 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance
    etaiov_cl_6 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance
    etaiov_cl_7 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance
    etaiov_cl_8 ~ fixed(0.62^2)  # Table 1, CL/F row [IOV 0.62]; shared variance

    # Residual error. Table 1 reports "29%/-" in the "Residual error (%/ug/L)"
    # row: a 29% proportional component and no additive component.
    propSd <- 0.29
    label("Proportional residual error (fraction)")  # Table 1, Residual error row: 29%
  })

  model({
    # Decompose the integer-valued OCC column into binary occasion indicators
    # for IOV multiplexing on CL/F. OCC = 1..8 selects the matching per-occasion
    # eta; OCC = 0 or any value outside 1..8 zeros every indicator and leaves
    # CL/F at its typical value.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8

    # Individual parameters. CL/F and V/F are apparent oral parameters: the
    # adjusted model was fitted to oral data only and does not separate F.
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl + iov_cl) * exp(e_csa_cl * CONMED_CSA)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
