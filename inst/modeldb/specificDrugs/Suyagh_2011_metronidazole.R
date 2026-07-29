Suyagh_2011_metronidazole <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous metronidazole in",
    "32 preterm neonates receiving treatment of or prophylaxis against",
    "necrotising enterocolitis, with dried-blood-spot HPLC sampling",
    "(Suyagh 2011). Clearance is described by an allometric 3/4-power",
    "scaling on body weight (reference 1.0 kg) and a linear postmenstrual-",
    "age maturation term centred at 30 weeks; volume of distribution is",
    "proportional to body weight. The publication is open access only at",
    "the abstract level, so inter-individual variability on CL and V and",
    "the residual-error magnitude are FIXED at 0 here; users running",
    "stochastic VPCs must supply their own variability terms (see the",
    "validation vignette's Errata section)."
  )
  reference <- paste(
    "Suyagh M, Collier PS, Millership JS, Iheagwaram G, Millar M,",
    "Halliday HL, McElnay JC.",
    "Metronidazole population pharmacokinetics in preterm neonates",
    "using dried blood-spot sampling.",
    "Pediatrics. 2011;127(2):e367-e374.",
    "doi:10.1542/peds.2010-0807. PMID 21220396.",
    sep = " "
  )
  vignette <- "Suyagh_2011_metronidazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying body weight (kg). Used as the allometric size",
        "descriptor on CL with the canonical 3/4-power exponent and",
        "reference 1.0 kg, and as the linear size descriptor on V",
        "(V = 0.726 * WT). Preterm-neonate cohort; baseline weights are",
        "not enumerated in the abstract, so a virtual cohort spanning",
        "0.7-2.0 kg (typical preterm range at the abstract's PMA range",
        "of 25-32 weeks) is used in the validation vignette."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = paste(
        "Postmenstrual age in months. Canonical units for PAGE are",
        "months (PAGE = GA_weeks / 4.35 + postnatal_months); the source",
        "paper expresses postmenstrual age in weeks (PMA_weeks). The",
        "conversion PMA_weeks = PAGE * 4.35 is applied inside model() so",
        "that the paper's centring at 30 weeks (PAGE * 4.35 - 30) is",
        "preserved."
      ),
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the linear PMA maturation term on CL.",
        "Abstract reports a half-life trajectory from approximately 40 h",
        "at 25 weeks' PMA (PAGE 5.747 months) to 19 h at 32 weeks' PMA",
        "(PAGE 7.356 months). Source column 'PMA' in the paper is in",
        "weeks; the conversion factor 4.35 weeks per month matches the",
        "PAGE canonical (see inst/references/covariate-columns.md)."
      ),
      source_name        = "PMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    n_observations = 203L,
    age_range      = paste(
      "Postmenstrual age range approximately 25-32 weeks (PAGE",
      "approximately 5.7-7.4 months) per the abstract's half-life",
      "trajectory; the abstract does not tabulate exact PMA quantiles."
    ),
    weight_range   = "Not enumerated in the abstract (preterm neonate cohort).",
    disease_state  = paste(
      "Preterm neonates receiving intravenous metronidazole for the",
      "treatment of or prophylaxis against necrotising enterocolitis."
    ),
    dose_range     = paste(
      "Intravenous metronidazole; per-administration doses and",
      "dosing-interval distribution are not enumerated in the abstract.",
      "Conclusion section notes the paper suggests a weight-and-PMA-",
      "based dosing scheme for preterm neonates."
    ),
    regions        = "United Kingdom (Royal Maternity Hospital, Belfast; single-centre).",
    notes          = paste(
      "Sampling: 203 dried blood spots on filter paper, analysed by",
      "HPLC. Modelling: nonlinear mixed-effect modelling (specific",
      "NONMEM control stream not given in the abstract). The package",
      "extraction is sourced from the open-access abstract; the full",
      "text is paywalled. Random-effect variances and residual-error",
      "magnitudes are not reported in the abstract and are encoded as",
      "fixed(0) here (see ini() comments and the vignette Errata",
      "section)."
    )
  )

  ini({
    # =========================================================================
    # Structural PK parameters -- Suyagh 2011 abstract Results paragraph.
    # Reference subject for the typical values: WT = 1.00 kg, PMA = 30 weeks
    # (PAGE = 30 / 4.35 = 6.897 months). CL units L/h, V units L, WT units kg,
    # PMA units weeks (paper). Final-model equations:
    #   CL = 0.0247 * (WT / 1.00)^0.75 * (1 + 0.107 * (PMA - 30))
    #   V  = 0.726 * WT
    # =========================================================================
    lcl <- log(0.0247)
    label("Clearance for the 1.0 kg, 30-week PMA reference (CL, L/h)")
    # Suyagh 2011 abstract: CL typical value 0.0247 L/h.

    lvc <- log(0.726)
    label("Per-kg volume of distribution (V/WT, L/kg)")
    # Suyagh 2011 abstract: V = 0.726 * WT, so the typical V/WT slope is
    # 0.726 L/kg.

    # Allometric exponent on CL. The abstract reports the value as 0.75 with
    # no uncertainty and the canonical 3/4-power form, consistent with a
    # theoretical (fixed) Anderson-Holford allometric exponent rather than
    # an estimated one. Encoded as fixed() here; should the operator obtain
    # the full text and find the exponent was actually estimated, drop the
    # fixed() wrapper and supply the point estimate.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL (unitless; fixed at the canonical 3/4-power value)")
    # Suyagh 2011 abstract: (WT / 1.00)^0.75.

    # PMA maturation slope. The abstract reports postmenstrual age as a
    # "significant covariate" on CL with the linear-deviation form
    # (1 + 0.107 * (PMA - 30)). Significance testing implies the slope was
    # estimated (not fixed).
    e_pma_cl <- 0.107
    label("PMA maturation slope on CL (per week PMA, centred at 30 weeks)")
    # Suyagh 2011 abstract: (1 + 0.107 * (PMA - 30)).

    # =========================================================================
    # Inter-individual variability -- NOT REPORTED in the Suyagh 2011 abstract.
    # The paper used "nonlinear mixed-effect modelling" so an OMEGA structure
    # exists in the full text, but its magnitude is paywalled. Per project
    # convention for abstract-only extractions, IIV on CL and V is encoded as
    # fixed(0). The vignette Errata section flags this so downstream users do
    # not interpret the zero variances as the published estimate.
    # =========================================================================
    etalcl ~ fixed(0)
    etalvc ~ fixed(0)

    # =========================================================================
    # Residual error -- NOT REPORTED in the Suyagh 2011 abstract. The HPLC
    # assay specifications and the residual-error form are paywalled with the
    # full text; only the structural-model equations appear in the abstract.
    # Encoded as fixed(0) per project convention; the vignette Errata section
    # documents this as a transcription gap rather than a paper-derived value.
    # =========================================================================
    propSd <- fixed(0)
    label("Proportional residual SD (fraction; FIXED at 0 -- not reported in the abstract-only source)")
  })

  model({
    # PMA conversion: canonical PAGE is in months; the paper's covariate
    # equation uses postmenstrual age in weeks. PMA_weeks = PAGE * 4.35.
    pma_weeks <- PAGE * 4.35

    # Individual PK parameters with the allometric and PMA-maturation effects
    # on CL, and the linear-weight scaling on V.
    cl <- exp(lcl + etalcl) * (WT / 1.0)^e_wt_cl * (1 + e_pma_cl * (pma_weeks - 30))
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    # One-compartment IV disposition: central compartment receives the dose
    # directly (IV) and is eliminated first-order.
    d/dt(central) <- -kel * central

    # Concentration in mg/L: dose in mg, vc in L, central in mg.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
