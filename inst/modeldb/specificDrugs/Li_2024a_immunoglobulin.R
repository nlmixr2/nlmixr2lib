Li_2024a_immunoglobulin <- function() {
  description <- "One-compartment population PK model for intravenous polyclonal immunoglobulin G in multifocal motor neuropathy, scaled on lean body mass (Li 2024, Frontiers in Neurology)"
  reference <- "Li Z, Roepcke S, Franke R, Yel L. Dose, exposure, and treatment regimen of intravenous immunoglobulin G in multifocal motor neuropathy. Front Neurol. 2024;15:1478419. doi:10.3389/fneur.2024.1478419 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 39)"
  vignette <- "vanderZeeuw_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 56.54 kg, with an estimated allometric exponent of 2.23 on the volume of distribution -- far outside the physiological 0.75-1.0 range and the most extreme exponent in van der Zeeuw 2026 (section 3.2.2.3 reports it without comment). The exponent is transcribed as printed; users extrapolating outside the observed LBM range should expect implausible volumes. This model scales the volume on LBM only; no covariate is reported on elimination (section 3.2.2.4: 'Li et al. did not report allometric scaling on elimination').",
      source_name = "LBM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested but not retained: van der Zeeuw 2026 Table 3 lists 'Age, sex, BW, BMI, LBM, creatinine clearance' as covariates tested and only 'LBM on Vc' in the final model.",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Tested but not retained (van der Zeeuw 2026 Table 3). IgG is not renally eliminated, so a null result is expected.",
      source_name = "Creatinine clearance"
    )
  )

  compartmentData <- list(
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 44L,
    n_studies = 1L,
    age_range = "31.0-72.0 years",
    age_median = "52.0 years",
    weight_range = "56.3-107.0 kg",
    weight_median = "83.5 kg",
    sex_female_pct = round(100 * 12 / 44, 1),
    race_ethnicity = "Not reported",
    disease_state = "Multifocal motor neuropathy (MMN)",
    dose_range = "IVIg 400-2000 mg/kg once every 2-4 weeks",
    regions = "United States, Canada, Europe",
    notes = "Single study, NCT00666263 (van der Zeeuw 2026 Table 1). Treatment-naive (endogenous) IgG median 20.2 g/L, range 11.6-37.0. This model and Li_2024b_immunoglobulin.R were fitted to the SAME patient cohort; van der Zeeuw 2026 Figure 2 notes 'both models from Li et al. overlap as they report the same popPK model', with the Annals of Clinical and Translational Neurology version (Li_2024b) including additional data and adding a grip-strength pharmacodynamic layer. Only trough concentrations were available, which is why a one-compartment model was used (section 3.2.2.3)."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Li et al. (2024)
    # (Front Neurol.) [39]'. This model is parameterised with an elimination
    # RATE CONSTANT rather than a clearance -- Table 4 prints 'Kel = 0.05784'
    # in the CL column -- so lkel is used instead of lcl.
    lkel    <- log(0.05784); label("First-order elimination rate constant (1/day)")     # van der Zeeuw 2026 Table 4: Kel = 0.05784
    lvc     <- log(6.59);    label("Volume of distribution at LBM 56.54 kg (L)")        # van der Zeeuw 2026 Table 4: Vc = 6.59 (LBM/56.54)^2.23

    # Allometric exponent on lean body mass -- estimated.
    e_lbm_vc <- 2.23; label("Allometric exponent on Vc (unitless)")                     # van der Zeeuw 2026 Table 4 and section 3.2.2.3: (LBM/56.54)^2.23

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.2.5: 'Li et al. estimated
    # a typical endogenous IgG value (CBASE) with associated IIV [39].' The
    # estimate itself is NOT reported in Table 4 or anywhere in the review, so
    # the cohort's observed treatment-naive median from Table 1 is used.
    # Wrapped in fixed() because it is an observed summary statistic
    # substituting for an unreported parameter; see the vignette Errata.
    bl_igg  <- fixed(20.2); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 Table 1: treatment-naive/endogenous IgG median 20.2 g/L (11.6-37.0); the estimated CBASE is not reported

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2). Table 4 prints 'NR' in the CL column of
    # the IIV block for this row, so no eta is declared on elimination.
    etalvc ~ 0.031886  # 18.0% CV; van der Zeeuw 2026 Table 4 IIV 'V = 18.0'

    # Residual error. van der Zeeuw 2026 Table 4 prints a VARIANCE in the
    # proportional column, so it is entered here as its square root.
    propSd <- 0.095131; label("Proportional residual error (fraction)")                 # van der Zeeuw 2026 Table 4: Prop sigma^2 = 0.00905; sqrt = 0.095131
  })
  model({
    kel <- exp(lkel)
    vc  <- exp(lvc + etalvc) * (LBM / 56.54)^e_lbm_vc

    # Intravenous administration only: doses go directly to `central`.
    # The state holds EXOGENOUS (therapeutic) IgG only.
    d/dt(central) <- -kel * central

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
