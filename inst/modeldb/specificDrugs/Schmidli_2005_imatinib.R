Schmidli_2005_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption of ",
    "duration D1 into the central compartment and first-order elimination ",
    "for oral imatinib 400 mg once daily in adults with chronic-phase ",
    "chronic myeloid leukaemia, from the phase III IRIS study (Schmidli ",
    "2005). Both CL/F and Vc/F carry an additive shift between the day-1 ",
    "occasion and the day-29-or-later occasion, plus power-form effects of ",
    "total body weight (reference 80 kg), hemoglobin (reference 13 g/dL) ",
    "and white blood cell count (reference 16 x 10^9/L). Inter-individual ",
    "variability is estimated on CL/F and Vc/F as a correlated 2x2 block; ",
    "residual error is combined proportional plus additive. TRANSCRIBED ",
    "FROM A SECONDARY SOURCE: the parameter values come from Table 1 of ",
    "Yang 2025, an external evaluation of 15 published imatinib population ",
    "PK models, not from the primary publication. Re-extract from Schmidli ",
    "2005 when that paper is obtained."
  )
  reference <- paste0(
    "Schmidli H, Peng B, Riviere GJ, Capdeville R, Hensley M, Gathmann I, ",
    "Bolton AE, Racine-Poon A. Population pharmacokinetics of imatinib ",
    "mesylate in patients with chronic-phase chronic myeloid leukaemia: ",
    "results of a phase III study. Br J Clin Pharmacol. ",
    "2005;60(1):35-44. doi:10.1111/j.1365-2125.2005.02372.x. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Schmidli et al. ",
    "(2005)' and Table 1 footnote b. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    OCC = list(
      description        = "Occasion indicator: 1 = the day-1 (first-dose) occasion, 2 = the day-29-or-later occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (day 1, first dose)",
      notes              = paste0(
        "Yang 2025 Table 1 footnote b defines the source column as ",
        "'Occasion (OCC) with OCC = 0 for day 1 and OCC = 1 (used in this ",
        "current validation study) for day >= 29'. The canonical OCC column ",
        "is 1-based (see inst/references/covariate-columns.md), so the ",
        "paper's OCC = 0 maps to OCC = 1 here and the paper's OCC = 1 maps ",
        "to OCC = 2. The model decomposes OCC into the binary indicator ",
        "occ_late = (OCC == 2), which reproduces the published equation ",
        "term for term. This is a structural typical-value shift, not an ",
        "inter-occasion random effect: no eta is attached to the occasion. ",
        "It captures the well-described fall in imatinib CL/F between the ",
        "first dose and steady state. Yang 2025 evaluated the model at ",
        "OCC = 1 in the paper's coding (day >= 29), i.e. OCC = 2 here, ",
        "which is the appropriate setting for therapeutic drug monitoring ",
        "data. Same column and coding as Demetri_2009_imatinib.R and ",
        "Petain_2008_imatinib.R, which Yang 2025 annotates with the same ",
        "footnote b."
      ),
      source_name        = "OCC"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as (WT/80)^0.301 and Vc/F as (WT/80)^0.405 (Yang 2025 ",
        "Table 1). The reference 80 kg is the centring constant printed ",
        "inside the covariate term. Note that neither exponent is the ",
        "allometric 0.75 / 1.0 pair, so Yang 2025 additionally evaluated a ",
        "variant of this model in which standard allometric scaling ",
        "REPLACED the published weight terms; that variant is Yang 2025's ",
        "own modification and is NOT encoded here. The allometrically ",
        "scaled form of this model was among the two best performers in ",
        "Yang 2025 (median prediction error -8.82%, median absolute ",
        "prediction error 39.7%, Table 3)."
      ),
      source_name        = "TBW"
    ),
    HGB = list(
      description        = "Hemoglobin concentration",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as (HGB/13)^0.897 and Vc/F as (HGB/13)^0.676 (Yang ",
        "2025 Table 1). Yang 2025 Table 1 abbreviation list gives the unit ",
        "explicitly: 'HG hemoglobin (g/dL)'. The reference 13 g/dL is the ",
        "centring constant printed inside the covariate term and is a ",
        "plausible central value for a chronic-phase CML cohort. Data ",
        "assemblers holding hemoglobin in g/L must divide by 10 before ",
        "using this model. The register permits either unit for HGB and ",
        "requires the per-model unit to be declared here, which it is."
      ),
      source_name        = "HG"
    ),
    WBC = list(
      description        = "White blood cell count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as (WBC/16)^-0.105 and Vc/F as (WBC/16)^-0.07 (Yang ",
        "2025 Table 1). Yang 2025 Table 1 abbreviation list: 'WBC white ",
        "blood cell count (10^9/L)'. The reference 16 x 10^9/L is the ",
        "centring constant printed inside the covariate term; it is high ",
        "relative to a healthy population because the IRIS cohort was ",
        "sampled in newly diagnosed chronic-phase CML, where leukocytosis ",
        "is the defining laboratory abnormality. Both exponents are small ",
        "and negative."
      ),
      source_name        = "WBC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 371L,
    n_studies      = 1L,
    n_observations = "1930 imatinib plasma concentrations (Yang 2025 Table 1)",
    age_range      = "18-70 years",
    disease_state  = "Adults with chronic-phase chronic myeloid leukaemia (CML)",
    dose_range     = "Oral imatinib 400 mg total daily dose",
    regions        = "Multi-country international phase III study",
    bioanalytical  = "LC-MS/MS, limit of quantification 25 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "The largest of the three international-trial-derived models in Yang ",
      "2025 Table 1 apart from Gotta 2014. Demographic detail beyond the ",
      "row above (weight range, sex split, race) is not reported by the ",
      "secondary source and must be read from the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Schmidli row) -----
    # Typical values are for the REFERENCE subject: day-1 occasion
    # (OCC == 1), WT = 80 kg, HGB = 13 g/dL, WBC = 16 x 10^9/L, at which
    # every covariate factor equals 1.
    lcl <- log(13.8); label("Apparent oral clearance CL/F at the reference subject on the day-1 occasion (L/h)")  # Yang 2025 Table 1: CL/F = (13.8 + (-3.18) x OCC) x ...
    lvc <- log(252); label("Apparent central volume Vc/F at the reference subject on the day-1 occasion (L)")  # Yang 2025 Table 1: Vc/F = (252 + (-7.82) x OCC) x ...
    ld1 <- log(1.5); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 1.5

    # ----- Occasion shifts (ADDITIVE, in the parameter's own units) -----
    # These are additive offsets inside the parentheses of the published
    # equation, not multiplicative factors, so they are carried on the
    # linear scale and added before the covariate factors are applied.
    e_occ_cl <- -3.18; label("Additive shift in CL/F on the day-29-or-later occasion (L/h)")  # Yang 2025 Table 1: (13.8 + (-3.18) x OCC)
    e_occ_vc <- -7.82; label("Additive shift in Vc/F on the day-29-or-later occasion (L)")  # Yang 2025 Table 1: (252 + (-7.82) x OCC)

    # ----- Covariate exponents -----
    e_wt_cl <- 0.301; label("Power exponent of (WT / 80) on CL/F (unitless)")  # Yang 2025 Table 1
    e_hgb_cl <- 0.897; label("Power exponent of (HGB / 13) on CL/F (unitless)")  # Yang 2025 Table 1
    e_wbc_cl <- -0.105; label("Power exponent of (WBC / 16) on CL/F (unitless)")  # Yang 2025 Table 1
    e_wt_vc <- 0.405; label("Power exponent of (WT / 80) on Vc/F (unitless)")  # Yang 2025 Table 1
    e_hgb_vc <- 0.676; label("Power exponent of (HGB / 13) on Vc/F (unitless)")  # Yang 2025 Table 1
    e_wbc_vc <- -0.07; label("Power exponent of (WBC / 16) on Vc/F (unitless)")  # Yang 2025 Table 1

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports 'CV%(CL): 31.9%', 'CV%(Vc): 31.4%' and
    # 'Cov(eta_CL, eta_Vc): 0.071'. Following the convention used throughout
    # this transcription (and the shipped Jiang_2023_imatinib.R precedent),
    # the tabulated CV% is taken as omega, the standard deviation of the
    # log-scale random effect, so the variance is (CV/100)^2. The resulting
    # correlation is 0.071 / (0.319 * 0.314) = 0.709, which is admissible.
    etalcl + etalvc ~ c(0.101761,
                        0.071, 0.098596)  # Yang 2025 Table 1: omega^2 = 0.319^2 and 0.314^2; Cov = 0.071

    # ----- Residual unexplained variability -----
    # The additive term is tabulated as 0.249 mg/L; Cc is reported in ng/mL,
    # so 0.249 mg/L x 1000 = 249 ng/mL.
    propSd <- 0.26; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 26%
    addSd <- 249; label("Additive residual error (ng/mL)")  # Yang 2025 Table 1: Add 0.249 mg/L = 249 ng/mL
  })

  model({
    # ----- 1. Occasion indicator -----
    # Canonical OCC is 1-based; occ_late is 1 on the day-29-or-later
    # occasion and 0 on the day-1 occasion, matching the paper's own 0/1
    # OCC coding term for term.
    occ_late <- (OCC == 2)

    # ----- 2. Covariate factors -----
    cl_cov <- (WT / 80)^e_wt_cl * (HGB / 13)^e_hgb_cl * (WBC / 16)^e_wbc_cl
    vc_cov <- (WT / 80)^e_wt_vc * (HGB / 13)^e_hgb_vc * (WBC / 16)^e_wbc_vc

    # ----- 3. Individual parameters -----
    # The occasion shift is additive on the typical value and is applied
    # BEFORE the multiplicative covariate factors and the exponential eta,
    # exactly as the published equation is bracketed.
    cl <- (exp(lcl) + e_occ_cl * occ_late) * cl_cov * exp(etalcl)
    vc <- (exp(lvc) + e_occ_vc * occ_late) * vc_cov * exp(etalvc)
    d1 <- exp(ld1)

    # ----- 4. Micro-constants -----
    kel <- cl / vc

    # ----- 5. ODE system -----
    # One compartment, zero-order input. Dose records must carry rate = -2
    # so rxode2 uses the modelled duration dur(central) = d1.
    d/dt(central) <- -kel * central
    dur(central) <- d1

    # ----- 6. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
