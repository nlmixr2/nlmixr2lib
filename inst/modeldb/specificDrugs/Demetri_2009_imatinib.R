Demetri_2009_imatinib <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption of ",
    "duration D1 into the central compartment and first-order elimination ",
    "for oral imatinib in adults with unresectable or metastatic ",
    "gastrointestinal stromal tumor (Demetri 2009). CL/F and Vc/F both ",
    "carry power-form effects of serum albumin (reference 38.3 g/L) and ",
    "white blood cell count (reference 7 x 10^9/L), and Vc/F additionally ",
    "carries an additive shift between the day-1 occasion and the ",
    "day-29-or-later occasion. Inter-individual variability is estimated ",
    "on CL/F and Vc/F as a correlated 2x2 block; residual error is ",
    "combined proportional plus a very small additive term. TRANSCRIBED ",
    "FROM A SECONDARY SOURCE: the parameter values come from Table 1 of ",
    "Yang 2025, an external evaluation of 15 published imatinib population ",
    "PK models, not from the primary publication. Re-extract from Demetri ",
    "2009 when that paper is obtained; note in particular the unusual ",
    "repetition of the CL/F covariate exponents on Vc/F, flagged in the ",
    "covariate notes below."
  )
  reference <- paste0(
    "Demetri GD, Wang Y, Wehrle E, Racine A, Nikolova Z, Blanke CD, ",
    "Joensuu H, von Mehren M. Imatinib plasma levels are correlated with ",
    "clinical benefit in patients with unresectable/metastatic ",
    "gastrointestinal stromal tumors. J Clin Oncol. ",
    "2009;27(19):3141-3147. doi:10.1200/JCO.2008.20.4818. ",
    "PARAMETER SOURCE (secondary): Yang T, Rasmussen ASB, Weimann A, ",
    "Thastrup M, Rank CU, Als-Nielsen B, Malmros J, Wik HS, Lohi O, ",
    "Overgaard U, Johannsdottir IMR, Vaitkeviciene G, Dalhoff K, ",
    "Schmiegelow K, Lund TM. Published population pharmacokinetic models ",
    "of imatinib perform poorly on TDM data from pediatric patients. ",
    "Target Oncol. 2025;20(5):871-886. Table 1, row 'Demetri et al. ",
    "(2009)' and Table 1 footnote b. doi:10.1007/s11523-025-01172-2."
  )
  vignette <- "Yang_2025_imatinib_external_evaluation"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "imatinib", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as (ALB/38.3)^1.66 and Vc/F as (ALB/38.3)^1.66 (Yang ",
        "2025 Table 1). Yang 2025 Table 1 abbreviation list: 'ALB albumin ",
        "(g/L)', matching the canonical ALB unit. The reference 38.3 g/L ",
        "is the centring constant printed inside the covariate term. ",
        "TRANSCRIPTION CAVEAT: Yang 2025 Table 1 prints the SAME two ",
        "exponents (1.66 for albumin and -0.418 for white blood cell ",
        "count) on both CL/F and Vc/F. This was checked against the ",
        "original PDF with 'pdftotext -layout' rather than the trimmed ",
        "markdown, and the duplication is genuinely what the secondary ",
        "source prints -- it is not an artifact of table flattening. It is ",
        "nonetheless unusual for a covariate to act with an identical ",
        "exponent on a clearance and on a volume, so this is the highest ",
        "priority item to confirm when the primary is obtained. The ",
        "consequence if the primary differs is confined to Vc/F: the ",
        "steady-state trough prediction that drives imatinib TDM depends ",
        "mainly on CL/F, so the risk to the model's principal use is ",
        "limited."
      ),
      source_name        = "ALB"
    ),
    WBC = list(
      description        = "White blood cell count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters CL/F as (WBC/7)^-0.418 and Vc/F as (WBC/7)^-0.418 (Yang ",
        "2025 Table 1). Yang 2025 Table 1 abbreviation list: 'WBC white ",
        "blood cell count (10^9/L)'. The reference 7 x 10^9/L is the ",
        "centring constant printed inside the covariate term; unlike the ",
        "CML-derived models in the same table it sits in the normal adult ",
        "range, which is consistent with a solid-tumor (GIST) cohort. See ",
        "the ALB note above for the transcription caveat about the exponent ",
        "being repeated on Vc/F."
      ),
      source_name        = "WBC"
    ),
    OCC = list(
      description        = "Occasion indicator: 1 = the day-1 (first-dose) occasion, 2 = the day-29-or-later occasion",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (day 1, first dose)",
      notes              = paste0(
        "Yang 2025 Table 1 attaches footnote b to this model's Vc/F term ",
        "'(168 + 58.5 x OCC)'. Footnote b reads 'Occasion (OCC) with ",
        "OCC = 0 for day 1 and OCC = 1 (used in this current validation ",
        "study) for day >= 29'. The canonical OCC column is 1-based, so ",
        "the paper's OCC = 0 maps to OCC = 1 here and the paper's OCC = 1 ",
        "maps to OCC = 2. The effect is ADDITIVE and applies to Vc/F only, ",
        "not to CL/F: the apparent central volume rises from 168 L on the ",
        "first-dose occasion to 226.5 L at steady state. No eta is ",
        "attached to the occasion, so this is a structural typical-value ",
        "shift, not inter-occasion variability. Same column and coding as ",
        "Schmidli_2005_imatinib.R and Petain_2008_imatinib.R."
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 73L,
    n_studies      = 1L,
    n_observations = "not reported in Yang 2025 Table 1",
    age_range      = "25-79 years",
    disease_state  = "Adults with unresectable or metastatic gastrointestinal stromal tumor (GIST)",
    dose_range     = "Oral imatinib 400 or 600 mg total daily dose",
    regions        = "Europe and USA",
    bioanalytical  = "LC-MS, limit of quantification 4 ng/mL (Yang 2025 Table 1)",
    notes          = paste0(
      "This model was one of only four whose ORIGINAL (unscaled) form met ",
      "Yang 2025's bias criterion of a median prediction error within ",
      "+/- 15% on their external dataset (-13.71%, Table 3); its precision ",
      "was nonetheless poor (median absolute prediction error 47.36%). ",
      "Demographic detail beyond the row above (weight range, sex split, ",
      "race) is not reported by the secondary source and must be read from ",
      "the primary publication."
    )
  )

  ini({
    # ----- Structural parameters (Yang 2025 Table 1, Demetri row) -----
    # Typical values are for the REFERENCE subject: ALB = 38.3 g/L,
    # WBC = 7 x 10^9/L, on the day-1 occasion, at which every covariate
    # factor equals 1.
    lcl <- log(8.18); label("Apparent oral clearance CL/F at the reference subject (L/h)")  # Yang 2025 Table 1: CL/F = 8.18 x (ALB/38.3)^1.66 x (WBC/7)^-0.418
    lvc <- log(168); label("Apparent central volume Vc/F at the reference subject on the day-1 occasion (L)")  # Yang 2025 Table 1: Vc/F = (168 + 58.5 x OCC) x (ALB/38.3)^1.66 x (WBC/7)^-0.418
    ld1 <- log(1.69); label("Zero-order absorption duration D1 (h)")  # Yang 2025 Table 1: D1 = 1.69

    # ----- Occasion shift on Vc/F (ADDITIVE, in litres) -----
    e_occ_vc <- 58.5; label("Additive shift in Vc/F on the day-29-or-later occasion (L)")  # Yang 2025 Table 1: (168 + 58.5 x OCC)

    # ----- Covariate exponents -----
    # See the covariateData notes: the same two exponents are printed on
    # both CL/F and Vc/F in the secondary source. They are carried as four
    # separate parameters (rather than two shared ones) so that a future
    # correction against the primary can change the Vc/F pair without
    # touching the CL/F pair.
    e_alb_cl <- 1.66; label("Power exponent of (ALB / 38.3) on CL/F (unitless)")  # Yang 2025 Table 1
    e_wbc_cl <- -0.418; label("Power exponent of (WBC / 7) on CL/F (unitless)")  # Yang 2025 Table 1
    e_alb_vc <- 1.66; label("Power exponent of (ALB / 38.3) on Vc/F (unitless)")  # Yang 2025 Table 1
    e_wbc_vc <- -0.418; label("Power exponent of (WBC / 7) on Vc/F (unitless)")  # Yang 2025 Table 1

    # ----- Inter-individual variability -----
    # Yang 2025 Table 1 reports 'CV%(CL): 34.6%', 'CV%(Vc): 35.7%' and
    # 'Cov(eta_CL, eta_Vc): 0.119'. The tabulated CV% is taken as omega
    # (the log-scale standard deviation), so the variance is (CV/100)^2.
    # This row is what settles that convention for the whole
    # transcription: with omega = CV the implied correlation is
    # 0.119 / (0.346 * 0.357) = 0.963, which is admissible, whereas the
    # exact log-normal conversion omega^2 = log(1 + CV^2) would give
    # 0.119 / sqrt(0.1128 * 0.1198) = 1.024 > 1, which is impossible.
    etalcl + etalvc ~ c(0.119716,
                        0.119, 0.127449)  # Yang 2025 Table 1: omega^2 = 0.346^2 and 0.357^2; Cov = 0.119

    # ----- Residual unexplained variability -----
    # The additive term is tabulated as 0.004 mg/L; Cc is reported in
    # ng/mL, so 0.004 mg/L x 1000 = 4 ng/mL. It is far below the assay
    # limit of quantification (4 ng/mL) and is effectively negligible
    # against the 34.9% proportional term, but it is retained so the
    # published error model is reproduced exactly.
    propSd <- 0.349; label("Proportional residual error (fraction)")  # Yang 2025 Table 1: Prop 34.9%
    addSd <- 4; label("Additive residual error (ng/mL)")  # Yang 2025 Table 1: Add 0.004 mg/L = 4 ng/mL
  })

  model({
    # ----- 1. Occasion indicator -----
    occ_late <- (OCC == 2)

    # ----- 2. Individual parameters -----
    # The occasion shift on Vc/F is additive on the typical value and is
    # applied BEFORE the multiplicative covariate factors and the
    # exponential eta, exactly as the published equation is bracketed.
    cl <- exp(lcl + etalcl) * (ALB / 38.3)^e_alb_cl * (WBC / 7)^e_wbc_cl
    vc <- (exp(lvc) + e_occ_vc * occ_late) *
      (ALB / 38.3)^e_alb_vc * (WBC / 7)^e_wbc_vc * exp(etalvc)
    d1 <- exp(ld1)

    # ----- 3. Micro-constants -----
    kel <- cl / vc

    # ----- 4. ODE system -----
    # One compartment, zero-order input. Dose records must carry rate = -2
    # so rxode2 uses the modelled duration dur(central) = d1.
    d/dt(central) <- -kel * central
    dur(central) <- d1

    # ----- 5. Observation and error -----
    # central is mg and vc is L, so central/vc is mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
