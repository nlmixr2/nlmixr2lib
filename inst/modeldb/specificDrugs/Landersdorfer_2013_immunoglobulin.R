Landersdorfer_2013_immunoglobulin <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G (IVIg and SCIg) in primary immunodeficiency (Landersdorfer 2013)"
  reference   <- "Landersdorfer CB, Bexon M, Edelman J, Rojavin M, Kirkpatrick CM, Lu J, et al. Pharmacokinetic modeling and simulation of biweekly subcutaneous immunoglobulin dosing in primary immunodeficiency. Postgrad Med. 2013;125(6):53-61. doi:10.3810/pgm.2013.11.2712 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 40)"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Body weight was the only covariate tested (van der Zeeuw 2026 Table 3) but was NOT retained on any parameter in the final model. Landersdorfer 2013 is absent from every 'BW on CL' and 'BW on Vc' list in van der Zeeuw 2026 sections 3.2.1.3 and 3.2.1.4, and its Table 4 row carries bare parameter values with no allometric term.",
      source_name = "BW"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "immunoglobulin G", units = "g", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 151L,
    n_studies      = 4L,
    age_range      = "3-81 years (per-study medians 18.0-32.0 years)",
    age_median     = "NCT00168025: 25.0; NCT00322556: 23.0; NCT00419341: 32.0; NCT00542997: 18.0 years",
    weight_range   = "13.0-135.0 kg (per-study medians 53.5-66.5 kg)",
    weight_median  = "NCT00168025: 66.5; NCT00322556: 62.0; NCT00419341: 66.0; NCT00542997: 53.5 kg",
    sex_female_pct = round(100 * (34 + 29 + 27 + 16) / (46 + 34 + 26 + 29 + 22 + 27 + 35 + 16), 1),
    race_ethnicity = "Not reported",
    disease_state  = "Primary immunodeficiency (PID) on immunoglobulin replacement therapy",
    dose_range     = "IVIg 200-888 mg/kg once every 3 to 4 weeks (NCT00168025, NCT00322556); SCIg 54-406 and 72-262 mg/kg weekly (NCT00419341, NCT00542997)",
    regions        = "United States, Europe",
    notes          = "Pooled analysis of four clinical trials (NCT00168025, NCT00322556, NCT00419341, NCT00542997). Demographics and dosing from van der Zeeuw 2026 Tables 1 and 2. Baseline IgG was not reported. IIV and residual error were not reported anywhere in the secondary source; see the ini() notes and the vignette Errata."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Landersdorfer
    # et al. (2013) [40]'. No covariate terms appear on this row: the model
    # carries bare typical values in L and L/day.
    lcl     <- log(0.142); label("Clearance (L/day)")                                   # van der Zeeuw 2026 Table 4: CL = 0.142
    lvc     <- log(3.94);  label("Central volume of distribution (L)")                  # van der Zeeuw 2026 Table 4: Vc = 3.94
    lq      <- log(0.252); label("Intercompartmental clearance (L/day)")                # van der Zeeuw 2026 Table 4: Q = 0.252
    lvp     <- log(4.18);  label("Peripheral volume of distribution (L)")               # van der Zeeuw 2026 Table 4: Vp = 4.18
    lka     <- log(0.439); label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 Table 4: Ka = 0.439
    lfdepot <- log(0.660); label("Subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 = 66.0%

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'In several studies,
    # the endogenous IgG concentration was fixed to 4 g/L [40, 41, 44, 54]' --
    # reference 40 is this model. Held constant, so wrapped in fixed().
    bl_igg  <- fixed(4); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 section 3.2.1.6, held constant at 4 g/L

    # Residual unexplained variability. van der Zeeuw 2026 Table 4 prints 'NR'
    # in BOTH residual-error columns for this row, and section 3.2.1.7 states
    # 'The remaining study did not report the error structure [40]'. Set to
    # zero rather than invented; see the vignette Errata.
    propSd  <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported in source)")  # van der Zeeuw 2026 Table 4: Prop = NR

    # Inter-individual variability is 'NR' in every IIV column of this row, so
    # no eta terms are declared at all. A zero-variance eta would make OMEGA
    # singular; a typical-value-only model is the faithful encoding.
  })
  model({
    cl     <- exp(lcl)
    vc     <- exp(lvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    ka     <- exp(lka)
    fdepot <- exp(lfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # `depot` receives subcutaneous doses; `central` receives intravenous
    # doses directly. The states hold EXOGENOUS (therapeutic) IgG only;
    # endogenous IgG enters at the observation step via `bl_igg`.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    f(depot)          <-  fdepot

    # Observed total plasma IgG = exogenous concentration + endogenous
    # baseline. Dose in g, volume in L -> g/L; bl_igg in g/L.
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
