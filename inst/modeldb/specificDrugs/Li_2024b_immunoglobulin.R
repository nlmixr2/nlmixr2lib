Li_2024b_immunoglobulin <- function() {
  description <- "One-compartment population PK model with a grip-strength indirect-response pharmacodynamic layer for intravenous polyclonal immunoglobulin G in multifocal motor neuropathy (Li 2024, Annals of Clinical and Translational Neurology)"
  reference   <- "Li Z, Roepcke S, Franke R, Yel L. Dose-exposure-efficacy response of intravenous immunoglobulin G 10% in multifocal motor neuropathy. Ann Clin Transl Neurol. 2024;11(8):1977-87. doi:10.1002/acn3.52107 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 23), section 3.3.1, and Supplementary Equations S1-S2"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  paper_specific_compartments <- c("gs")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 56.54 kg, with an estimated allometric exponent of 2.17 on the volume of distribution -- far outside the physiological 0.75-1.0 range. Transcribed as printed; users extrapolating outside the observed LBM range should expect implausible volumes. No covariate is reported on elimination.",
      source_name        = "LBM"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested but not retained: van der Zeeuw 2026 Table 3 lists 'Age, sex, BW, BMI, LBM, creatinine clearance' as covariates tested and only 'LBM on Vc' in the final model. Covariate relationships were also tested on the pharmacodynamic baseline G_BASE and on the residual-error term, and none were retained (section 3.3.1).",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested but not retained (van der Zeeuw 2026 Table 3). IgG is not renally eliminated, so a null result is expected.",
      source_name = "Creatinine clearance"
    )
  )

  compartmentData <- list(
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    gs      = list(analyte = "grip strength (latent pharmacodynamic state)", units = "kg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 1L,
    age_range      = "31.0-72.0 years",
    age_median     = "52.0 years",
    weight_range   = "56.3-107.0 kg",
    weight_median  = "83.5 kg",
    sex_female_pct = round(100 * 12 / 44, 1),
    race_ethnicity = "Not reported",
    disease_state  = "Multifocal motor neuropathy (MMN)",
    dose_range     = "IVIg 400-2000 mg/kg once every 2-4 weeks",
    regions        = "United States, Canada, Europe",
    notes          = "Single study, NCT00666263 (van der Zeeuw 2026 Table 1). Treatment-naive (endogenous) IgG median 20.2 g/L, range 11.6-37.0. Same patient cohort as Li_2024a_immunoglobulin.R; this version includes additional data and adds the grip-strength pharmacodynamic layer. Grip strength was measured by dynamometer. A change of 4 kg was taken as the minimal clinically meaningful improvement; Monte Carlo simulation showed >=70% of patients reaching that threshold at IVIg doses >=1600 mg/kg/month, and dose-splitting did not affect grip strength above 1000 mg/kg (section 3.3.1)."
  )

  ini({
    # ---- Pharmacokinetics -------------------------------------------------
    # van der Zeeuw 2026 Table 4, row 'Li et al. (2024) (Ann Clin Transl
    # Neurol.) [23]'. Parameterised with an elimination RATE CONSTANT rather
    # than a clearance, so lkel is used instead of lcl.
    lkel     <- log(0.05832); label("First-order elimination rate constant (1/day)")    # van der Zeeuw 2026 Table 4: Kel = 0.05832
    lvc      <- log(6.48);    label("Volume of distribution at LBM 56.54 kg (L)")       # van der Zeeuw 2026 Table 4: Vc = 6.48 (LBM/56.54)^2.17
    e_lbm_vc <- 2.17;         label("Allometric exponent on Vc (unitless)")             # van der Zeeuw 2026 Table 4: (LBM/56.54)^2.17

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.2.5 reports that a typical
    # endogenous IgG value was estimated, but no point estimate appears in
    # Table 4 or elsewhere in the review, so the cohort's observed
    # treatment-naive median from Table 1 is used. This value also defines the
    # pharmacodynamic driver DRV (IgG in excess of baseline); see model().
    bl_igg   <- fixed(20.2); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 Table 1: treatment-naive/endogenous IgG median 20.2 g/L (11.6-37.0)

    # ---- Pharmacodynamics: grip strength ----------------------------------
    # Indirect-response model in which IgG in excess of baseline INHIBITS the
    # deterioration (loss) of grip strength. van der Zeeuw 2026 Equation 1:
    #   dGS/dt = MNT - (1 - Imax*DRV/(C50 + DRV)) * DTR * GS
    # MNT is the production rate of grip strength (canonical `kin`), DTR the
    # deterioration rate (canonical `kout`), C50 the DRV giving half-maximal
    # inhibition (canonical `ec50`), and DRV the IgG concentration in excess
    # of the baseline IgG concentration.
    #
    # G_BASE is the reported estimate, so the model is parameterised on it and
    # MNT is derived as kin = rbase * kout -- the steady-state identity that
    # falls out of Equation 1 at DRV = 0. Supplementary Equation S1 prints this
    # identity as 'Gbase = MNT/DRV'; DRV there is a typo for DTR, since
    # MNT/DRV is dimensionally inconsistent and contradicts Equation 1. See
    # the vignette Errata.
    lrbase   <- log(10.6); label("Baseline grip strength in the absence of treatment (kg)")  # van der Zeeuw 2026 section 3.3.1: G_BASE = 10.6 kg

    # Deterioration rate, reported as 0.023 per HOUR. Converted to the model's
    # day time base by multiplying by 24 -- the same hour-to-day conversion the
    # review applies to PK estimates in section 2.3.
    lkout    <- log(0.552); label("Grip-strength deterioration rate constant (1/day)")  # van der Zeeuw 2026 section 3.3.1: DTR = 0.023 1/hour; x24 = 0.552 1/day

    lec50    <- log(9.41);  label("DRV giving half-maximal inhibition of grip-strength loss (g/L)")  # van der Zeeuw 2026 section 3.3.1: C50 = 9.41 mg/mL (= 9.41 g/L)

    # Maximum inhibitory effect, reported on the logit scale as LIMAX = 0.433.
    # Supplementary Equation S2: Imax = exp(LIMAX) / (1 + exp(LIMAX)) =
    # 0.6066. The review's own back-transform in section 3.3.1 gives 0.61,
    # confirming the value (its stated unit 'mg/mL' there is a slip -- Imax is
    # a dimensionless fraction).
    imax     <- 0.6066;    label("Maximum fractional inhibition of grip-strength loss (unitless)")  # van der Zeeuw 2026 section 3.3.1 and Supplementary Equation S2: LIMAX = 0.433 -> Imax = 0.6066

    # ---- Random effects ---------------------------------------------------
    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2). Table 4 prints 'NR' in the CL column of
    # the IIV block, so no eta is declared on elimination.
    etalvc   ~ 0.033652  # 18.5% CV; van der Zeeuw 2026 Table 4 IIV 'V = 18.5'
    etalrbase ~ 0.561570 # 86.8% CV; van der Zeeuw 2026 section 3.3.1 IIV on G_BASE

    # Residual error, one term per output.
    propSd <- 0.095131; label("Proportional residual error on total IgG (fraction)")  # van der Zeeuw 2026 Table 4: Prop sigma^2 = 0.00905; sqrt = 0.095131
    propSd_gs <- 0.0602;   label("Proportional residual error on grip strength (fraction)")  # van der Zeeuw 2026 section 3.3.1: proportional error term 0.0602 (shrinkage 7.0%)
    # NOTE: the source also estimates inter-individual variability ON the
    # grip-strength residual-error term (24.5% CV, shrinkage 0.8%). nlmixr2
    # has no direct encoding for IIV on sigma, so that layer is omitted; see
    # the vignette Errata.
  })
  model({
    kel   <- exp(lkel)
    vc    <- exp(lvc + etalvc) * (LBM / 56.54)^e_lbm_vc
    rbase <- exp(lrbase + etalrbase)
    kout  <- exp(lkout)
    ec50  <- exp(lec50)

    # Steady-state identity from Equation 1 evaluated at DRV = 0.
    kin <- rbase * kout

    # Intravenous administration only: doses go directly to `central`.
    # The state holds EXOGENOUS (therapeutic) IgG only.
    d/dt(central) <- -kel * central

    # DRV: the IgG concentration IN EXCESS of the baseline IgG concentration.
    # Since Cc = central/vc + bl_igg, the excess is exactly central/vc.
    drv <- central / vc

    # Grip strength: IgG inhibits the deterioration (kout) arm.
    d/dt(gs) <- kin - (1 - imax * drv / (ec50 + drv)) * kout * gs
    gs(0)    <- rbase

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg

    Cc ~ prop(propSd)
    gs ~ prop(propSd_gs)
  })
}
