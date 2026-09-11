Hartmann_2026_nintedanib_fvcz <- function() {
  description <- paste0(
    "Pediatric exposure-response model for the FVC Z-score under ",
    "nintedanib in children and adolescents 6 to less than 18 years of age ",
    "with clinically significant fibrosing interstitial lung disease ",
    "(Hartmann 2026; the phase 3 InPedILD trial and its open-label ",
    "extension InPedILD-ON). The FVC Z-score is a single state starting at ",
    "an estimated baseline and changing at a constant annual rate: a ",
    "linear placebo disease-progression slope, offset for pediatric ",
    "patients, plus an Emax disease-modifying drug effect driven by the ",
    "individual steady-state nintedanib trough concentration. ",
    "Inter-individual variability sits on baseline, on slope and on the ",
    "residual-error magnitude; residual error is additive. Estimated in ",
    "NONMEM with the NWPRI frequentist-prior functionality using the adult ",
    "FVC Z-score exposure-response meta-model as prior; slope, Emax, EC50 ",
    "and the inter-individual variability on slope were supported by that ",
    "prior, while baseline, residual error, the variability terms on ",
    "baseline and on residual error, and the pediatric slope offset were ",
    "estimated from the pediatric data alone. There is no PK layer: the ",
    "exposure metric is supplied as the CTROUGH data column. Sibling model ",
    "for the percent-predicted form of the same endpoint: ",
    "modellib('Hartmann_2026_nintedanib_fvcpp')."
  )
  reference <- paste(
    "Hartmann S, Chan Kwong A, Ribbing J, Gahlemann M, Korell J.",
    "Population Pharmacokinetics and Exposure-Response Model-Based",
    "Bayesian Extrapolation of FVC-Based Efficacy Endpoints From Adults to",
    "Pediatric Patients Receiving Nintedanib.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70135.",
    "doi:10.1002/psp4.70135. PMCID PMC12823301.",
    "Parameter values are the final estimates in Table 3; the model",
    "structure is Section 3.5 together with the parallel FVC percent",
    "predicted NONMEM control stream reproduced in Data S1, which Section",
    "3.5 states the Z-score model mirrors.",
    "The exposure metric derives from the companion pediatric population",
    "pharmacokinetic model in the same paper; see",
    "modellib('Hartmann_2026_nintedanib').",
    sep = " "
  )
  vignette <- "Hartmann_2026_nintedanib"

  # fvcz is a first sighting of this endpoint as an ODE state in the
  # library, so it is declared paper-specific rather than promoted to a
  # canonical compartment. A second FVC Z-score model is the trigger to
  # promote it alongside the register's existing fev1pp entry.
  paper_specific_compartments <- c("fvcz")

  units <- list(
    time          = "year",
    dosing        = "(no dose events; nintedanib exposure enters through the CTROUGH covariate column, in nM)",
    concentration = "Z-score (FVC expressed as standard deviations from the age-, sex- and height-standardised reference mean; the modelled state is a dimensionless standardised score rather than a drug concentration, so the dosing-versus-concentration dimensional check is not applicable and the dosing string is parenthesised to skip it)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    fvcz = list(analyte = "forced vital capacity, standardised Z-score", units = "Z-score", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual model-predicted steady-state nintedanib plasma trough concentration, the exposure driver of the Emax disease-modifying effect.",
      units              = "nM",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma concentration, at STEADY STATE, and an",
        "EMPIRICAL-BAYES prediction rather than an observed trough.",
        "Hartmann 2026 Section 2.2: Individual PK model predicted",
        "nintedanib trough concentrations at steady-state (Ctrough,ss) were",
        "used as exposure metrics for ER models. Predicted Ctrough,ss",
        "values were obtained using empirical Bayes estimates from the",
        "popPK model based on the pediatric patients included in the model.",
        "Section 2.4 adds that the values fed to the ER models were derived",
        "for all patients from the recorded dosing and demographic",
        "information, taking dose adjustments and treatment interruptions",
        "into account. The column is therefore TIME-VARYING within a",
        "subject: it steps whenever the dose level changes and drops to 0",
        "during a treatment interruption and for every placebo subject.",
        "Reproduce it with modellib('Hartmann_2026_nintedanib') solved to",
        "steady state on the patient's current weight-band dose.",
        "UNITS ARE LOAD-BEARING AND ARE nM, not ng/mL. Table 3 gives the",
        "EC50 unit as nM, and the paper reports every nintedanib",
        "concentration in nM (Figure 1 caption adult reference geometric",
        "means Cmax,ss 33 nM, Ctrough,ss 20 nM, Cav,ss 26 nM,",
        "AUCtau,ss 316 nM h).",
        "Enters through a SATURABLE Emax form emax * CTROUGH / (ec50 +",
        "CTROUGH), the same form the parallel FVC percent predicted control",
        "stream codes as DREFF = EMAX*CPREPRED/(EC50+CPREPRED).",
        "The Z-score EC50 of 8.12 nM is estimated with POOR precision",
        "(RSE 102%, Table 3) but is close to the 8.05 nM of the",
        "percent-predicted sibling, and both sit below the adult reference",
        "Ctrough,ss geometric mean of 20 nM."
      ),
      source_name        = "CPREPRED"
    ),
    CHILD = list(
      description        = "Pediatric-versus-adult indicator (1 = pediatric patient, under 18 years of age; 0 = adult). Additive offset on the annual rate of change in the FVC Z-score.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult), the level the adult prior describes.",
      notes              = paste(
        "Time-fixed per subject. Effect on the annual rate of change:",
        "slope = slope_placebo + e_child_slope_placebo * CHILD, i.e. an",
        "ADDITIVE offset in Z-score per year, mirroring the parallel FVC",
        "percent predicted control stream (Data S1) TVPLSL = PLSLCOV +",
        "TVPLSL.",
        "Age cutoff: under 18 years. Hartmann 2026 Table S3 footnote b",
        "lists the three forms in which pediatric age was tested on slope,",
        "the third being dichotomous paediatric patients versus adults;",
        "Table 3 labels the retained parameter Pediatric change in slope,",
        "and Section 3.5 describes the model as a linear placebo model with",
        "a separate annual rate of decline for pediatric patients.",
        "EVERY patient in this analysis is pediatric, so CHILD = 1",
        "reproduces the paper's pediatric predictions; the value is kept as",
        "a covariate rather than folded into the slope so that setting",
        "CHILD = 0 recovers the adult slope of -0.308 Z-score per year for",
        "the adult-versus-pediatric comparison the paper draws in Figures 2",
        "to 4."
      ),
      source_name        = "STUDYN337 (the InPedILD study indicator, which in the pooled analysis data set is the pediatric-patient flag)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age, screened on baseline, slope, Emax and the inter-individual variability on slope",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened in the stepwise covariate search (Hartmann 2026",
        "Table S3) as continuous with a cut-off at 18 years, and as two",
        "dichotomous contrasts. Only the dichotomous",
        "pediatric-versus-adult contrast on slope was retained; that one is",
        "carried as the CHILD covariate."
      )
    ),
    DIS_SSC_ILD = list(
      description = "Systemic-sclerosis-associated ILD indicator, carried by the adult model as a study effect on baseline and slope",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Hartmann 2026 Section 2.4: the adult models retained a study",
        "effect for SENSCIS, an SSc-ILD population, on slope and on",
        "baseline. It was NOT retained in the pediatric model because a",
        "diagnosis with SSc-ILD was not sufficiently prevalent in the",
        "pediatric trial population (9 of 53 patients, Table 1). Carried as",
        "documentation only."
      )
    ),
    RACE_JAPANESE = list(
      description = "Japanese-heritage race indicator, carried by the adult model on baseline and on slope",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Hartmann 2026 Section 2.4 lists Chinese, Korean, Indian or",
        "Japanese ethnicity on baseline and Japanese ethnicity on slope",
        "among the covariates identified for the ADULT models. None was",
        "retained in the pediatric model because the pediatric trial",
        "population did not have any patients from the respective",
        "ethnicities. Carried as documentation only; the same applies to",
        "RACE_CHINESE, RACE_KOREAN and RACE_INDIAN."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 53L,
    n_studies      = 2L,
    n_observations = "517 FVC observations",
    age_range      = "6 to less than 18 years (trial eligibility); cohort mean 13.2 years (SD 3.27)",
    age_median     = "mean 13.2 years (SD 3.27); no median reported",
    weight_range   = "13.5 kg lower eligibility bound; cohort mean 43.1 kg (SD 18.3)",
    weight_median  = "mean 43.1 kg (SD 18.3); no median reported",
    sex_female_pct = 60.4,
    race_ethnicity = "Not tabulated for the 53-patient FVC set. The overlapping 44-patient PK set was 77% Caucasian, 9.1% Black, 6.8% American Indian/Alaska Native, 4.5% Other Asian and 2.3% missing (Hartmann 2026 Table S4). No patient was Chinese, Korean, Indian or Japanese.",
    disease_state  = "Clinically significant fibrosing interstitial lung disease of mixed aetiology: pediatric autoimmune ILD 32%, surfactant protein deficiency 26%, other ILDs 25%, toxic/radiation/drug-induced pneumonitis 11%, chronic hypersensitivity pneumonia 3.8%, post-HSCT fibrosis 1.9%. Systemic sclerosis-associated ILD in 17%.",
    dose_range     = "Oral nintedanib twice daily, dosed by body-weight bin (50, 75, 100 or 150 mg BID; Hartmann 2026 Table S1), or matching placebo. InPedILD randomised 2:1 nintedanib to placebo over 24 weeks; the extension is open-label active treatment.",
    regions        = "Multinational; the InPedILD trial and the InPedILD-ON open-label extension",
    baseline_endpoint = "FVC Z-score at baseline: mean -3.48 (SD 1.82) overall; -2.97 (SD 1.35) placebo and -3.99 (SD 1.73) active in the 6 to less than 12 year group, -3.27 (SD 2.18) placebo and -3.37 (SD 1.84) active in the 12 to less than 18 year group (Hartmann 2026 Table 1).",
    treatment_arm_breakdown = c(`6_to_lt12_placebo_n` = 4, `6_to_lt12_active_n` = 13, `12_to_lt18_placebo_n` = 9, `12_to_lt18_active_n` = 27),
    notes          = paste0(
      "Exposure-response analysis set: all 53 patients enrolled in ",
      "InPedILD and its open-label extension InPedILD-ON. Baseline ",
      "characteristics are Hartmann 2026 Table 1. The FVC Z-score is the ",
      "second of the two derived response measures modelled, alongside FVC ",
      "percent predicted; separate exposure-response models were developed ",
      "for each endpoint on the same 517 observations. The Z-score is the ",
      "more informative of the two in a growing population because it is ",
      "standardised for age, sex and height rather than for a single ",
      "predicted value. Because the pediatric data set is small the model ",
      "was estimated with the adult FVC Z-score exposure-response ",
      "meta-model as a frequentist prior through the NONMEM NWPRI ",
      "functionality, after an external validation step. The trial was NOT ",
      "powered to estimate efficacy."
    )
  )

  # Implementation notes (see the vignette section 'Assumptions and
  # deviations' for the full justification of each item):
  # * Time is in YEARS. The slope, the pediatric slope offset and Emax
  #   are all rates per year (Table 3 unit column), and they are
  #   integrated directly as the derivative of the state, so the
  #   integration variable has to be years for the units to close.
  # * No FVC Z-score control stream is included in the supplement; Data
  #   S1 carries only the pediatric popPK and the FVC percent predicted
  #   streams. The structure encoded here is the one Section 3.5 states
  #   in words -- Similarly to the model for FVC %predicted, the final
  #   model for FVC Z-score was a linear placebo model with a separate
  #   annual rate of decline for pediatric patients, and an Emax model
  #   describing the disease-modifying effect of nintedanib on FVC
  #   Z-score. IIV terms were supported on RUV, baseline and slope. RUV
  #   was described by an additive model -- mapped onto the sibling
  #   stream equation for equation. Two arithmetic identities confirm
  #   the mapping: -0.308 + 0.133 = -0.175 against the Table 3 footnote
  #   a value of -0.174 Z-score per year for a pediatric patient quoted
  #   again in Section 3.5, and the Discussion contrast of -0.174 per
  #   year in children against -0.309 per year in adults.
  # * The ONE structural difference from the percent-predicted sibling
  #   is the scale of the baseline eta, and it is forced by the sign of
  #   the endpoint. Table 3's unit column gives IIV baseline as an SD of
  #   1.84 on the Z-score scale, where Table 2 gives its baseline IIV as
  #   a CV of 0.401. A Z-score baseline of -3.49 cannot carry a
  #   lognormal eta, so the baseline eta here is ADDITIVE and the
  #   typical value is on the natural scale (rbase, not lrbase). Reading
  #   1.84 as a CV instead would be a category error and would put the
  #   whole cohort on the wrong side of zero.
  # * The eta on slope is likewise ADDITIVE, in Z-score per year
  #   (Table 3 unit column 'SD'), exactly as in the sibling model. With
  #   a typical pediatric slope of about -0.175 per year and a
  #   between-subject SD of 0.379 per year, an individual slope may take
  #   either sign.
  # * Residual error is additive with a per-subject magnitude, mirroring
  #   the sibling stream's Y = IPRED + EPS(1)*EXP(ETA(3)). rxode2's
  #   add() helper takes a plain variable name, so the eta-scaled
  #   magnitude is materialised into addSd_fvcz_i first.
  # * A placebo patient, or a patient during a treatment interruption,
  #   carries CTROUGH = 0, which makes the Emax term exactly 0 and
  #   leaves the linear disease-progression slope.
  ini({
    # ----- Baseline (estimated from the pediatric data, without the adult prior) -----
    # Bare, NOT log-transformed: the typical baseline Z-score is negative.
    rbase <- -3.49; label("Typical baseline FVC Z-score (Z-score)")  # Table 3 row 'Baseline' = -3.49 Z-score (RSE 7.29%). Section 3.5: baseline, RUV and IIV on baseline and RUV were estimated independently from the prior. Table 1 gives the observed cohort mean baseline as -3.48 (SD 1.82)

    # ----- Linear placebo disease-progression slope (supported by the adult prior) -----
    slope_placebo <- -0.308; label("Typical adult annual rate of change in FVC Z-score under placebo (Z-score/year)")  # Table 3 row 'Slope' = -0.308 Z-score/year (RSE 4.71%), footnote b 'Supported by adult priors'. The Discussion quotes the adult rate as -0.309 Z-score/year

    # ----- Pediatric offset on the slope (released from the adult prior) -----
    e_child_slope_placebo <- 0.133; label("Additive change in the annual rate of change in FVC Z-score for a pediatric patient (Z-score/year)")  # Table 3 row 'Pediatric change in slope' = 0.133 Z-score/year (RSE 49.7%), footnote a 'Translating into a slope of -0.174 Z-score/year for a pediatric patient'; -0.308 + 0.133 = -0.175. Section 3.5: with the covariate describing the difference in slope for pediatric patients, the pediatric slope was released from the adult prior

    # ----- Emax drug effect (supported by the adult prior) -----
    emax <- 0.292; label("Maximum improvement in the annual rate of change in FVC Z-score over placebo (Z-score/year)")  # Table 3 row 'Rate of change in FVC Z-score at maximum drug effect (Emax)' = 0.292 Z-score/year (RSE 30.7%), footnote b 'Supported by adult priors'. Section 3.5 repeats: the Emax was an improvement in the rate of decline of 0.292 Z-score/year (difference from placebo)
    ec50 <- 8.12; label("Steady-state nintedanib trough concentration producing half the maximum effect on FVC Z-score (nM)")  # Table 3 row 'EC50' = 8.12 nM (RSE 102%), footnote b 'Supported by adult priors'. Section 3.5 repeats: The EC50 estimated for the final FVC Z-score model was 8.12 nM. The 102% RSE is the least precisely estimated parameter in either exposure-response model

    # ----- Inter-individual variability -----
    # Table 3 states the scale of each term in its own unit column: CV for
    # the exponential eta on the residual magnitude, SD on the endpoint's
    # own scale for the two additive etas. Variances are the squares.
    etarbase         ~ 3.3856     # Table 3 row 'IIV baseline' SD = 1.84 Z-score (RSE 9.78%, shrinkage 0%); 1.84^2 = 3.3856. ADDITIVE eta on the Z-score scale, not a CV: a negative typical baseline forbids a lognormal eta. Estimated without the adult prior
    etaslope_placebo ~ 0.143641   # Table 3 row 'IIV slope' SD = 0.379 Z-score/year (RSE 2.18%, shrinkage 17.3%), footnote b 'Supported by adult priors'; 0.379^2 = 0.143641. ADDITIVE eta in Z-score/year
    etaaddSd_fvcz    ~ 0.215296   # Table 3 row 'IIV RUV' CV = 0.464 (RSE 12.7%, shrinkage 5.72%); 0.464^2 = 0.215296. Estimated without the adult prior

    # ----- Residual unexplained variability -----
    addSd_fvcz <- 0.276; label("Typical additive residual SD on the FVC Z-score (Z-score)")  # Table 3 row 'Add. RUV' = 0.276 Z-score (RSE 7.74%, shrinkage 2.14%). Estimated without the adult prior
  })
  model({
    # ----- Individual baseline -----
    # Additive eta on the Z-score scale; see the implementation notes.
    base_fvcz <- rbase + etarbase

    # ----- Individual annual rate of change -----
    # Additive eta in Z-score per year, so an individual slope may be
    # positive.
    slope_i <- slope_placebo + e_child_slope_placebo * CHILD + etaslope_placebo

    # ----- Emax disease-modifying drug effect -----
    # A placebo patient, or a patient during a treatment interruption,
    # carries CTROUGH = 0 and therefore contributes exactly zero drug
    # effect.
    dreff <- emax * CTROUGH / (ec50 + CTROUGH)

    # ----- Disease-progression ODE, integrated in years -----
    fvcz(0)    <- base_fvcz
    d/dt(fvcz) <- slope_i + dreff

    # ----- Residual error: additive, with per-subject magnitude -----
    addSd_fvcz_i <- addSd_fvcz * exp(etaaddSd_fvcz)
    fvcz ~ add(addSd_fvcz_i)
  })
}
