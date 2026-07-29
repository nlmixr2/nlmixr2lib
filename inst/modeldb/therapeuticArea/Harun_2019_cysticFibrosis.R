Harun_2019_cysticFibrosis <- function() {
  description <- "Non-linear mixed-effects disease-progression model of forced expiratory volume in 1 second (FEV1) percent predicted versus age in children with cystic fibrosis (Harun 2019). The model describes the typical sigmoid-Emax decline of FEV1% predicted from a baseline at age 5 years to an asymptote, with covariate effects of BMI z-score and severe air trapping at age 5 on the baseline, and time-varying hospitalisation due to pulmonary exacerbation on the maximum drop and the half-effect age."
  reference <- paste(
    "Harun SN, Wainwright CE, Grimwood K, Hennig S; Australasian Cystic Fibrosis",
    "Bronchoalveolar Lavage (ACFBAL) study group. Aspergillus and progression of",
    "lung disease in children with cystic fibrosis. Thorax 2019;74(2):125-131.",
    "doi:10.1136/thoraxjnl-2018-211550.",
    sep = " "
  )
  vignette <- "Harun_2019_cysticFibrosis"
  units <- list(time = "year", dosing = "n/a (no dosing)", concentration = "% predicted (FEV1)")

  covariateData <- list(
    BMIZ = list(
      description        = "Body-mass-index z-score (age- and sex-standardised) at the time of the FEV1% predicted measurement",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying per-visit z-score; reference value is 0 (population mean for age and sex). Linear-deviation effect on baseline FEV1% predicted: e_bmi_baseline * (BMIZ - 0). Source paper does not state which growth-reference standard was used to compute the z-score; ACFBAL paediatric CF cohorts conventionally use the WHO 2007 Growth Reference for school-aged children and adolescents.",
      source_name        = "BMI"
    ),
    AIR_TRAP_5Y = list(
      description        = "Severe air trapping (Brody-II HRCT component score > 0) on the chest HRCT scan performed at age 5 years",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no air trapping at age 5)",
      notes              = "Time-fixed per subject; captures the single end-of-study HRCT performed at age 5 in the ACFBAL study. Coefficient e_at_baseline applies only to baseline FEV1% predicted at age 5 (not to dmax or t50max).",
      source_name        = "ATS5C"
    ),
    HOSPRA = list(
      description        = "Hospitalisation due to a pulmonary exacerbation at the time of the FEV1% predicted measurement",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not hospitalised at this visit)",
      notes              = "Time-varying per-visit indicator; affects both the maximum drop in FEV1% predicted (e_hpe_dmax) and the age at which 50% of the maximum drop occurs (e_hpe_t50max). Hospitalised visits accelerate both the magnitude and the onset of FEV1% decline.",
      source_name        = "HOSPRA"
    )
  )

  population <- list(
    n_subjects     = 79L,
    n_studies      = 2L,
    age_range      = "4.8-14.4 years",
    age_median     = "9.1 years",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 46.8,
    race_ethnicity = NA_character_,
    disease_state  = "Children with classic cystic fibrosis (two CFTR mutations, sweat chloride > 60 mmol/L, pancreatic insufficiency, or meconium ileus). Subset of the 156 ACFBAL participants who consented to longitudinal follow-up via the Australian Cystic Fibrosis Data Registry between ages 5 and 14 years.",
    dose_range     = "Not applicable -- disease-progression model with no drug dosing.",
    regions        = "Australia (multi-site: New South Wales, Northern Territory, Queensland, South Australia, Victoria, Australian Capital Territory) and New Zealand",
    notes          = "n=79 from the ACFBAL longitudinal sub-cohort; 2651 FEV1% predicted observations (median 33 per child) collected over a median of 8.0 years of follow-up. The age-5-years baseline FEV1% measurements were obtained from the ACFBAL study; subsequent FEV1% measurements at ages 6-14 years were obtained from the Australian Cystic Fibrosis Data Registry. Population characteristics from Table 1 of the source paper (n=79 column) and Table 2."
  )

  ini({
    # Structural parameters -- final non-linear progression model, Table E3 of the
    # Harun 2019 supplementary material (PDF page 12 of supplementary-material-1).
    lfev1pp_baseline <- log(99.7)        ; label("Typical baseline FEV1% predicted at age 5 (% predicted)")  # Table E3 final, also NMTRAN $THETA(1)
    ldmax            <- fixed(log(40))   ; label("Maximum lifetime change in FEV1% predicted (% predicted) -- FIXED to literature value")  # Table E3 final, $THETA(2) FIX; literature value from Harun 2016 systematic review (PMID 26597232)
    lt50max          <- log(8.38)        ; label("Age at which 50% of the maximum FEV1% change occurs (years)")  # Table E3 final, $THETA(3)
    lhill            <- log(3.08)        ; label("Hill coefficient on age (unitless)")  # Table E3 final, $THETA(4)

    # Covariate effects -- linear-deviation form (1 + e * cov) per Harun 2019
    # supplementary text page 13 ("Structural parameters after inclusion of the
    # influencing factors"). Numerical values from Table E3 final-model column;
    # NMTRAN $THETA(6)..$THETA(9) carry the same values to additional precision.
    e_bmi_baseline   <- 0.0382           ; label("Fractional change in baseline FEV1% per BMI z-score unit deviation from 0")  # Table E3 (0.038); $THETA(6) = 0.0382
    e_at_baseline    <- -0.0417          ; label("Fractional change in baseline FEV1% for severe air trapping (AIR_TRAP_5Y = 1) at age 5")  # Table E3 (-0.04); $THETA(9) = -0.0417
    e_hpe_dmax       <- -0.22            ; label("Fractional change in maximum FEV1% drop when hospitalised due to a pulmonary exacerbation (HOSPRA = 1)")  # Table E3 (-0.22); $THETA(7)
    e_hpe_t50max     <- -0.235           ; label("Fractional change in age of 50% maximum FEV1% drop when hospitalised due to a pulmonary exacerbation (HOSPRA = 1)")  # Table E3 (-0.24); $THETA(8) = -0.235

    # Inter-individual variability -- Harun 2019 supplementary equations E2 and
    # E3 (PDF page 9). The paper uses an exponential stochastic model (Eq E2)
    # for t50max and a proportional stochastic model (Eq E3) for baseline FEV1%,
    # dmax, and the Hill coefficient. The OMEGA matrix is applied to the etas
    # in their respective forms in model() below. Variance values are the
    # NMTRAN $OMEGA entries on PDF page 20; bootstrap-median equivalents in
    # Table E3 are within rounding.
    etalhill                             ~ 0.272                       # $OMEGA(1,1); proportional IIV on Hill (BSV ~ sqrt(0.272) = 52.2%)
    etalfev1pp_baseline                  ~ 0.0172                      # $OMEGA(2,2); proportional IIV on baseline FEV1% (BSV ~ sqrt(0.0172) = 13.1%)
    etaldmax + etalt50max                ~ c(0.405, 0.15, 0.09)        # $OMEGA BLOCK(2): var(dmax), cov(dmax,t50max), var(t50max). Proportional IIV on dmax (BSV ~ sqrt(0.405) = 63.6%); exponential IIV on t50max (BSV ~ sqrt(0.09) = 30%)

    # Residual error -- additive on FEV1% predicted, $SIGMA via THETA(5) in NMTRAN
    addSd <- 9.32 ; label("Additive residual error on FEV1% predicted (% predicted)")  # Table E3 final (9.32 [9.27; 8.64, 9.96])
  })

  model({
    # The independent variable for this disease-progression model is the
    # subject's chronological age in years. The model has no dosing events; the
    # `t` variable inside model() is taken to be age in years (the dataset's
    # `time` column carries age in years). The NMTRAN source converts the input
    # column from hours to years via TIME1 = ((TIME/24)/365); in the nlmixr2
    # data convention the operator carries age directly in the time column.
    age_yr <- t

    # Individual structural parameters. Each parameter applies its IIV and
    # covariate effects in the same multiplicative form as the NMTRAN $PRED
    # block (PDF pages 19-20 of the supplementary material).
    fev1pp_baseline <- exp(lfev1pp_baseline) *
      (1 + etalfev1pp_baseline) *
      (1 + e_at_baseline * AIR_TRAP_5Y) *
      (1 + e_bmi_baseline * (BMIZ - 0))
    dmax <- exp(ldmax) *
      (1 + etaldmax) *
      (1 + e_hpe_dmax * HOSPRA)
    t50max <- exp(lt50max + etalt50max) *
      (1 + e_hpe_t50max * HOSPRA)
    hill <- exp(lhill) *
      (1 + etalhill)

    # Sigmoid-Emax disease-progression model (Equation 1 of Harun 2019; the
    # right-hand-side recurs as Equation E2 of the supplementary material).
    fev1pp <- fev1pp_baseline -
      (dmax * age_yr^hill) / (t50max^hill + age_yr^hill)

    fev1pp ~ add(addSd)
  })
}
