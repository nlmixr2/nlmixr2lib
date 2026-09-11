Hartmann_2026_nintedanib_fvcpp <- function() {
  description <- paste0(
    "Pediatric exposure-response model for FVC percent predicted under ",
    "nintedanib in children and adolescents 6 to less than 18 years of age ",
    "with clinically significant fibrosing interstitial lung disease ",
    "(Hartmann 2026; the phase 3 InPedILD trial and its open-label ",
    "extension InPedILD-ON). FVC percent predicted is a single state ",
    "starting at an estimated baseline and changing at a constant annual ",
    "rate: a linear placebo disease-progression slope, offset for ",
    "pediatric patients, plus an Emax disease-modifying drug effect driven ",
    "by the individual steady-state nintedanib trough concentration. ",
    "Inter-individual variability sits on baseline, on slope and on the ",
    "residual-error magnitude; residual error is additive. Estimated in ",
    "NONMEM with the NWPRI frequentist-prior functionality using the adult ",
    "FVC percent predicted exposure-response meta-model as prior; slope, ",
    "Emax, EC50 and the inter-individual variability on slope were ",
    "supported by that prior, while baseline, residual error, the ",
    "variability terms on baseline and on residual error, and the ",
    "pediatric slope offset were estimated from the pediatric data alone. ",
    "There is no PK layer: the exposure metric is supplied as the CTROUGH ",
    "data column. Sibling model for the Z-score form of the same endpoint: ",
    "modellib('Hartmann_2026_nintedanib_fvcz')."
  )
  reference <- paste(
    "Hartmann S, Chan Kwong A, Ribbing J, Gahlemann M, Korell J.",
    "Population Pharmacokinetics and Exposure-Response Model-Based",
    "Bayesian Extrapolation of FVC-Based Efficacy Endpoints From Adults to",
    "Pediatric Patients Receiving Nintedanib.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70135.",
    "doi:10.1002/psp4.70135. PMCID PMC12823301.",
    "Parameter values are the final estimates in Table 2; the model",
    "equations are the final pediatric FVC percent predicted NONMEM",
    "control stream reproduced in Data S1.",
    "The exposure metric derives from the companion pediatric population",
    "pharmacokinetic model in the same paper; see",
    "modellib('Hartmann_2026_nintedanib').",
    sep = " "
  )
  vignette <- "Hartmann_2026_nintedanib"

  # fvcpp is a first sighting of this endpoint as an ODE state in the
  # library, so it is declared paper-specific rather than promoted to a
  # canonical compartment. The register already carries fev1pp (FEV1
  # percent predicted, Harun 2019 cystic fibrosis); a second FVC percent
  # predicted model is the trigger to promote fvcpp alongside it.
  paper_specific_compartments <- c("fvcpp")

  units <- list(
    time          = "year",
    dosing        = "(no dose events; nintedanib exposure enters through the CTROUGH covariate column, in nM)",
    concentration = "% predicted (FVC percent predicted; the modelled state is a spirometric endpoint expressed as a percentage of a reference-equation predicted value, not a drug concentration, so the dosing-versus-concentration dimensional check is not applicable and the dosing string is parenthesised to skip it)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    fvcpp = list(analyte = "forced vital capacity, percent of predicted", units = "% predicted", specimen = "not applicable", verified = TRUE)
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
        "Perfect adherence to the randomized dose of nintedanib and",
        "steady-state conditions was assumed for the predictions; however,",
        "subjects' current dosage regimen was accounted for.",
        "Section 2.4 adds that the values fed to the ER models were derived",
        "for all patients from the recorded dosing and demographic",
        "information, taking dose adjustments and treatment interruptions",
        "into account. The column is therefore TIME-VARYING within a",
        "subject: it steps whenever the dose level changes and drops to 0",
        "during a treatment interruption and for every placebo subject.",
        "Reproduce it with modellib('Hartmann_2026_nintedanib') solved to",
        "steady state on the patient's current weight-band dose.",
        "UNITS ARE LOAD-BEARING AND ARE nM, not ng/mL. Table 2 gives the",
        "EC50 unit as nM, and the paper reports every nintedanib",
        "concentration in nM (Figure 1 caption adult reference geometric",
        "means Cmax,ss 33 nM, Ctrough,ss 20 nM, Cav,ss 26 nM,",
        "AUCtau,ss 316 nM h). Supplying ng/mL instead would inflate the",
        "driver by the 1.853 ng/mL per nM factor and push the model",
        "further onto the Emax plateau.",
        "Enters through a SATURABLE Emax form emax * CTROUGH / (ec50 +",
        "CTROUGH), not linearly and not on a log scale, per the final",
        "control stream DREFF = EMAX*CPREPRED/(EC50+CPREPRED).",
        "For scale: the estimated EC50 is 8.05 nM against an adult",
        "reference Ctrough,ss geometric mean of 20 nM at the approved",
        "150 mg twice-daily dose, so the typical adult and the typical",
        "pediatric patient both sit ABOVE the EC50, on the rising shoulder",
        "of the curve; the paper states that most pediatric patients",
        "exceeded the EC50 and some reached the EC80."
      ),
      source_name        = "CPREPRED (the control stream $INPUT column; CPREPREDN is its dose-normalised companion used only for visual predictive checks)"
    ),
    CHILD = list(
      description        = "Pediatric-versus-adult indicator (1 = pediatric patient, under 18 years of age; 0 = adult). Additive offset on the annual rate of change in FVC percent predicted.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult), the level the adult prior describes.",
      notes              = paste(
        "Time-fixed per subject. Effect on the annual rate of change:",
        "slope = slope_placebo + e_child_slope_placebo * CHILD, i.e. an",
        "ADDITIVE offset in percent predicted per year, from the final",
        "control stream (Data S1) TVPLSL = PLSLCOV + TVPLSL with",
        "PLSLCOV = PLSLSTUDYN337 = THETA(5).",
        "Age cutoff: under 18 years. Hartmann 2026 Table S3 footnote b",
        "lists the three forms in which pediatric age was tested on slope,",
        "the third being dichotomous paediatric patients versus adults;",
        "Table 2 labels the retained parameter Pediatric change in slope,",
        "and Section 3.4 describes it as a covariate effect describing the",
        "change in pediatric annual rate of decline.",
        "The control stream spells the same covariate PLSLSTUDYN337 after",
        "the study number of InPedILD (11990337) because the PsN stepwise",
        "covariate search names a covariate after the data column that",
        "carries it, and in the pooled adult-plus-pediatric analysis data",
        "set that column flags exactly the pediatric patients. The paper's",
        "own text, table label and covariate-screening table all describe",
        "the retained effect as pediatric versus adult, so the canonical",
        "pediatric indicator is used here.",
        "EVERY patient in this analysis is pediatric, so CHILD = 1",
        "reproduces the paper's pediatric predictions; the value is kept as",
        "a covariate rather than folded into the slope so that setting",
        "CHILD = 0 recovers the adult slope of -4.74 percent predicted per",
        "year for the adult-versus-pediatric comparison the paper draws in",
        "Figures 2 to 4."
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
        "carried as the CHILD covariate. No age effect on baseline, on",
        "Emax or on the variability terms survived."
      )
    ),
    DIS_SSC_ILD = list(
      description = "Systemic-sclerosis-associated ILD indicator, carried by the adult model as a study effect on baseline and slope",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Hartmann 2026 Section 2.4: the adult FVC percent predicted model",
        "retained a study effect for SENSCIS, an SSc-ILD population, on",
        "slope and on baseline. It was NOT retained in the pediatric model",
        "because a diagnosis with SSc-ILD was not sufficiently prevalent in",
        "the pediatric trial population (9 of 53 patients, Table 1).",
        "Carried as documentation only; the parameter is not in the",
        "pediatric model and no pediatric estimate exists."
      )
    ),
    RACE_JAPANESE = list(
      description = "Japanese-heritage race indicator, carried by the adult model on baseline and on slope",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Hartmann 2026 Section 2.4 lists Chinese, Korean, Indian or",
        "Japanese ethnicity on baseline and Japanese ethnicity on slope",
        "among the covariates identified for the ADULT FVC percent",
        "predicted model. None was retained in the pediatric model because",
        "the pediatric trial population did not have any patients from the",
        "respective ethnicities. Carried as documentation only; the same",
        "applies to RACE_CHINESE, RACE_KOREAN and RACE_INDIAN, which are",
        "not repeated as separate entries here."
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
    baseline_endpoint = "FVC percent predicted at baseline: mean 59.3% (SD 21.2) overall; 64.7% (SD 17.2) placebo and 52.9% (SD 20.6) active in the 6 to less than 12 year group, 62.1% (SD 25.6) placebo and 60.7% (SD 20.9) active in the 12 to less than 18 year group (Hartmann 2026 Table 1).",
    treatment_arm_breakdown = c(`6_to_lt12_placebo_n` = 4, `6_to_lt12_active_n` = 13, `12_to_lt18_placebo_n` = 9, `12_to_lt18_active_n` = 27),
    notes          = paste0(
      "Exposure-response analysis set: all 53 patients enrolled in ",
      "InPedILD and its open-label extension InPedILD-ON, nine more than ",
      "contributed to the companion population PK model. Baseline ",
      "characteristics are Hartmann 2026 Table 1. FVC was measured at ",
      "screening, at the week 0 baseline visit and at several later visits ",
      "of variable total duration; screening observations were excluded ",
      "from the fit (control stream IGNORE(PLCBL.LT.0)). Because the ",
      "pediatric data set is small the model was estimated with the adult ",
      "FVC percent predicted exposure-response meta-model as a frequentist ",
      "prior through the NONMEM NWPRI functionality, after an external ",
      "validation step. The trial was NOT powered to estimate efficacy."
    )
  )

  # Implementation notes (see the vignette section 'Assumptions and
  # deviations' for the full justification of each item):
  # * Time is in YEARS. The slope, the pediatric slope offset and Emax
  #   are all rates per year (Table 2 unit column), and the control
  #   stream integrates them directly as DADT(1) = PLSL + DREFF, so the
  #   integration variable has to be years for the units to close.
  # * Parameter values are the FINAL estimates in Hartmann 2026 Table 2,
  #   NOT the $THETA and $OMEGA values printed in the Data S1 control
  #   stream. Those are initial estimates: every prior-supported one is
  #   numerically identical to its $THETAP / $OMEGAP prior entry
  #   (-4.78463, 4.33916, 8.2449, 31.1633), which is the tell. Contrast
  #   the companion popPK control stream, whose $THETA block does carry
  #   the finals. Each line below quotes the Table 2 value and the
  #   control-stream initial value it supersedes.
  # * Two arithmetic identities in the source confirm the additive
  #   parameterisation of the pediatric slope offset:
  #   -4.74 + 2.28 = -2.46, against the Table 2 footnote b value of
  #   -2.45 percent predicted per year for a pediatric patient quoted
  #   again in Section 3.4; and the Discussion contrast of -2.45 per
  #   year in children against -4.78 per year in adults.
  # * Table 2 reports the variance terms in a mixture of scales, marked
  #   in its own unit column: IIV RUV and IIV baseline are given as CV
  #   (exponential etas) and IIV slope as an SD in percent predicted per
  #   year (an additive eta). The additive form on slope is load-bearing
  #   and is what the control stream codes
  #   (PLSL = TVPLSL + ETA(1)*TVETAPLSL): with a typical slope of about
  #   -2.5 per year and a between-subject SD of 5.55 per year, an
  #   individual pediatric slope can take either sign, which is what the
  #   paper's per-patient panels show. A lognormal eta could not do
  #   that.
  # * Residual error is additive with a per-subject magnitude:
  #   Y = IPRED + EPS(1)*EXP(ETA(3)) in the control stream. rxode2's
  #   add() helper takes a plain variable name, so the eta-scaled
  #   magnitude is materialised into addSd_fvcpp_i first, following the
  #   idiom already used by Ezzati_2014_dexmedetomidine_piglet.R and
  #   Tang_2023_tenecteplase.R.
  # * The control stream carries a commented-out combined
  #   proportional-plus-additive alternative
  #   (;Y = IPRED + (EPS(1)*IPRED + EPS(2))*EXP(ETA(3)) and
  #   ;$SIGMA 0.000882475 ; 1. RUV prop). It is commented out in the
  #   final run and Table 2 reports only an additive residual, so the
  #   additive form alone is encoded.
  # * There is NO drug effect on the placebo arm and no separate placebo
  #   model: a placebo patient simply has CTROUGH = 0, which makes the
  #   Emax term exactly 0 and leaves the linear disease-progression
  #   slope. This reproduces the paper's difference-from-placebo
  #   figures directly as the difference between a CTROUGH > 0 and a
  #   CTROUGH = 0 solve.
  ini({
    # ----- Baseline (estimated from the pediatric data, without the adult prior) -----
    lrbase <- log(54.9); label("Log typical baseline FVC percent predicted (% predicted)")  # Table 2 row 'Baseline' = 54.9 % (RSE 5.54%); Data S1 $THETA 4 initial (0,54.8049). Section 3.4: Baseline, RUV and IIV on baseline and RUV were estimated independently from the prior

    # ----- Linear placebo disease-progression slope (supported by the adult prior) -----
    slope_placebo <- -4.74; label("Typical adult annual rate of change in FVC percent predicted under placebo (% predicted/year)")  # Table 2 row 'Slope' = -4.74 %/year (RSE 4.24%), footnote a 'Supported by adult priors'; Data S1 $THETA 1 initial -4.78463, identical to the $THETAP adult prior. The Discussion quotes the adult rate as -4.78%/year

    # ----- Pediatric offset on the slope (released from the adult prior) -----
    e_child_slope_placebo <- 2.28; label("Additive change in the annual rate of change in FVC percent predicted for a pediatric patient (% predicted/year)")  # Table 2 row 'Pediatric change in slope' = 2.28 %/year (RSE 41.2%), footnote b 'Translating into a slope of -2.45FVC%predicted/year for a pediatric patient'; -4.74 + 2.28 = -2.46. Data S1 $THETA 5 initial (-1.00,4.7832,20). Section 3.4: with the covariate describing the difference in slope for pediatric patients, the pediatric slope was released from the adult prior

    # ----- Emax drug effect (supported by the adult prior) -----
    emax <- 4.17; label("Maximum improvement in the annual rate of change in FVC percent predicted over placebo (% predicted/year)")  # Table 2 row 'Rate of change in FVC %predicted at maximum drug effect (Emax)' = 4.17 %/year (RSE 12.6%), footnote a 'Supported by adult priors'; Data S1 $THETA 2 initial 4.33916, identical to the $THETAP adult prior
    ec50 <- 8.05; label("Steady-state nintedanib trough concentration producing half the maximum effect on FVC percent predicted (nM)")  # Table 2 row 'EC50' = 8.05 nM (RSE 28.8%), footnote a 'Supported by adult priors'; Data S1 $THETA 3 initial (0,8.2449), essentially the $THETAP adult prior 8.24491. Section 3.4 repeats: The EC50 estimated in the final model was 8.05 nM

    # ----- Inter-individual variability -----
    # Table 2 states the scale of each term in its own unit column: CV for
    # the two exponential etas, SD in %/year for the additive eta on
    # slope. Variances below are the squares of those Table 2 values.
    etalrbase        ~ 0.160801   # Table 2 row 'IIV baseline' CV = 0.401 (RSE 9.85%, shrinkage 0%); 0.401^2 = 0.160801. Data S1 $OMEGA 2 initial 0.161468. Estimated without the adult prior
    etaslope_placebo ~ 30.8025    # Table 2 row 'IIV slope' SD = 5.55 %/year (RSE 2.15%, shrinkage 26.1%), footnote a 'Supported by adult priors'; 5.55^2 = 30.8025. Data S1 $OMEGA 1 initial 30.8275 against the $OMEGAP adult prior 31.1633. ADDITIVE eta in %/year, not a CV
    etaaddSd_fvcpp   ~ 0.219024   # Table 2 row 'IIV RUV' CV = 0.468 (RSE 12.7%, shrinkage 5.48%); 0.468^2 = 0.219024. Data S1 $OMEGA 3 initial 0.223249. Estimated without the adult prior

    # ----- Residual unexplained variability -----
    addSd_fvcpp <- 3.26; label("Typical additive residual SD on FVC percent predicted (% predicted)")  # Table 2 row 'Add. RUV' = 3.26 % (RSE 7.77%, shrinkage 2.16%); Data S1 $SIGMA 2 initial 10.5292, whose square root is 3.24. Estimated without the adult prior
  })
  model({
    # ----- Individual baseline -----
    # Control stream: TVBASE = THETA(4); BASE = TVBASE*EXP(ETA(2));
    # A_0(1) = BASE.
    base_fvcpp <- exp(lrbase + etalrbase)

    # ----- Individual annual rate of change -----
    # Control stream: TVPLSL = THETA(1); TVPLSL = PLSLCOV + TVPLSL;
    # PLSL = TVPLSL + ETA(1)*TVETAPLSL. The eta is ADDITIVE, in
    # % predicted per year, so an individual slope may be positive.
    slope_i <- slope_placebo + e_child_slope_placebo * CHILD + etaslope_placebo

    # ----- Emax disease-modifying drug effect -----
    # Control stream: DREFF = EMAX*CPREPRED/(EC50+CPREPRED). A placebo
    # patient, or a patient during a treatment interruption, carries
    # CTROUGH = 0 and therefore contributes exactly zero drug effect.
    dreff <- emax * CTROUGH / (ec50 + CTROUGH)

    # ----- Disease-progression ODE -----
    # Control stream: DADT(1) = PLSL + DREFF, integrated in years.
    fvcpp(0)   <- base_fvcpp
    d/dt(fvcpp) <- slope_i + dreff

    # ----- Residual error: additive, with per-subject magnitude -----
    # Control stream: Y = IPRED + (EPS(1))*EXP(ETA(3)).
    addSd_fvcpp_i <- addSd_fvcpp * exp(etaaddSd_fvcpp)
    fvcpp ~ add(addSd_fvcpp_i)
  })
}
