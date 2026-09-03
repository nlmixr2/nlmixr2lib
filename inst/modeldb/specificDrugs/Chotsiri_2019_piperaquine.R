# Population pharmacokinetic-pharmacodynamic model for piperaquine given as
# dihydroartemisinin-piperaquine for seasonal malaria chemoprevention (SMC) in
# young children in Burkina Faso (Chotsiri 2019, Nat Commun 10(1):480;
# doi:10.1038/s41467-019-08297-9; parent trial NCT00941785).
#
# NOTE: this is NOT the same paper as Chotsiri_2019_lumefantrine.R, which is
# Chotsiri 2019 Clin Pharmacol Ther 106(6):1299 (doi:10.1002/cpt.1531).

Chotsiri_2019_piperaquine <- function() {
  description <- paste(
    "Population PK-PD model for oral piperaquine given as",
    "dihydroartemisinin-piperaquine (DHA-PQ) for seasonal malaria",
    "chemoprevention in 179 children aged 2.33-58.1 months in Burkina Faso",
    "(Chotsiri 2019 Nat Commun; NCT00941785). Two transit-absorption",
    "compartments followed by three-compartment disposition, with the",
    "absorption rate constant set equal to the transit rate constant",
    "(Savic 2007) so ktr = ka = 3 / MTT. Estimated with a NONMEM $PRIOR",
    "frequentist prior taken from the dense-sampling paediatric treatment",
    "study of Tarning 2012 in the same region. Fixed allometric scaling of",
    "all clearances (power 3/4) and volumes (power 1) on body weight",
    "centred at 18.0 kg, and a study-level relative bioavailability of",
    "0.726 that carries the unexplained lower exposure of this cohort",
    "relative to the prior study. Capillary (finger-prick) and venous",
    "plasma were fitted simultaneously through a population conversion",
    "factor: the disposition parameters are scaled to the capillary matrix",
    "and the venous prediction Cc is 0.380 x the capillary prediction Ccap.",
    "Coupled interval-censored time-to-event model for new P. falciparum",
    "infection: constant baseline hazard (6.28 infections/year) multiplied",
    "by a sigmoid Emax protective effect of venous piperaquine",
    "(IC50 = 3.66 ng/mL, gamma = 1.79, Emax fixed at 1) and by a",
    "dihydroartemisinin effect that sets the hazard to zero from 6 days",
    "before the first dose of each treatment round to 24 h after its last",
    "dose, exposed as the survivor function `sur`. Companion piperaquine",
    "model: modellib('Hoglund_2017_piperaquine').",
    sep = " "
  )
  reference <- paste(
    "Chotsiri P, Zongo I, Milligan P, Compaore YD, Some AF, Chandramohan D,",
    "et al. (2019). Optimal dosing of dihydroartemisinin-piperaquine for",
    "seasonal malaria chemoprevention in young children. Nature",
    "Communications 10(1):480. doi:10.1038/s41467-019-08297-9.",
    "The frequentist prior for the pharmacokinetic structure is Tarning J,",
    "Zongo I, Some FA, Rouamba N, Parikh S, Rosenthal PJ, et al. (2012).",
    "Population pharmacokinetics and pharmacodynamics of piperaquine in",
    "children with uncomplicated falciparum malaria. Clinical Pharmacology",
    "and Therapeutics 91(3):497-505. doi:10.1038/clpt.2011.254.",
    sep = " "
  )
  vignette <- "Chotsiri_2019_piperaquine"
  units <- list(time = "h", dosing = "mg (piperaquine base)", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = FALSE),
    transit1    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = FALSE),
    transit2    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = FALSE),
    cumhaz      = list(analyte = "Cumulative hazard of new P. falciparum infection", units = NA_character_, specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject at admission. Chotsiri 2019 Table 1 reports",
        "median (range) body weight 11.0 (4.20-18.3) kg in the 179 children",
        "of the PKPD group. Applied as a fixed allometric function on every",
        "clearance (exponent 0.75) and every volume (exponent 1.00)",
        "parameter, centred at 18.0 kg (Methods, 'Pharmacokinetic",
        "analysis': 'Body weight (BW) was introduced as an allometric",
        "function for all clearance (exponent of n = 0.75) and volume",
        "(exponent of n = 1.00) parameters, centralized to 18 kg of body",
        "weight according to the typical patient in the prior model',",
        "theta_i = theta * exp(eta_i) * (BW_i / 18.0)^n, Eq. 1). The 18.0 kg",
        "centring value comes from the typical child of the Tarning 2012",
        "prior model, not from this cohort's median, and Table 2 footnote c",
        "confirms the population estimates are 'for a typical child of",
        "18.0 kg body weight'. Body weight also determines the weight-band",
        "tablet count in the WHO 2010 and WHO 2015 dosing regimens",
        "(Supplementary Table 2).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CONMED_DHA = list(
      description        = "Dihydroartemisinin parasite-clearance coverage window (1 = covered)",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Time-varying record-level indicator, 1 while the dihydroartemisinin",
        "component of a DHA-PQ treatment round is deemed to eliminate any",
        "emerging blood-stage infection and 0 otherwise. Chotsiri 2019",
        "Eq. 4 defines DHA_EFF = 0 'from 6 days before the first dose to",
        "24 hours after the last dose' of each three-day treatment round and",
        "DHA_EFF = 1 elsewhere, so CONMED_DHA = 1 - DHA_EFF and the hazard",
        "term in model() is dhaeff = 1 - CONMED_DHA. With doses on days 0, 1",
        "and 2 of a round, the window spans t = -144 h to t = +72 h relative",
        "to that round's first dose (9 days). The six-day lead-in is not a",
        "drug-exposure window: it is the interval over which a blood-stage",
        "infection emerging from the liver would still be below the",
        "microscopy detection limit when the treatment arrives, and would",
        "therefore be cleared before it could be observed (Methods,",
        "'Pharmacodynamic analysis'). Build the column in the event table",
        "from the planned round start times; it is not derivable inside",
        "model() because it depends on FUTURE doses.",
        sep = " "
      ),
      source_name        = "DHA_EFF (reverse-coded)"
    )
  )

  # Screened by Chotsiri 2019 and not retained in the final model. The
  # stepwise addition (p < 0.05) / elimination (p < 0.01) search covered
  # parasitaemia, gender, age and nutritional status on the pharmacokinetic
  # parameters (Methods, 'Pharmacokinetic analysis') and body weight, gender
  # and age on the baseline hazard (Methods, 'Pharmacodynamic analysis');
  # 'No other covariates were significant' (Results) and 'No covariate
  # effects on the pharmacodynamic parameters were found' (Results). None is
  # referenced in model(), so they are documented here rather than in
  # covariateData.
  covariatesDataExcluded <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units       = "months",
      type        = "continuous",
      notes       = paste(
        "Driver of the enzyme-maturation factor evaluated on elimination",
        "clearance, MF = PMA^HILL / (TM50^HILL + PMA^HILL) (Chotsiri 2019",
        "Eq. 2, Anderson & Holford form). The Results state that including",
        "the maturation effect 'improved the model fit significantly, but",
        "resulted in an unrealistic and poorly estimated enzyme maturation",
        "half-life (TM50) of less than one month. Therefore, the maturation",
        "effect was omitted in the final pharmacokinetic model.' No TM50 or",
        "HILL estimate is published, so the rejected term is omitted",
        "entirely rather than carried with a null coefficient.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Chronological age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Chotsiri 2019 Table 1 reports median (range) age 32.1 (2.33-58.1)",
        "months in the PKPD group; the trial enrolled children aged 3-59",
        "months. Screened on the pharmacokinetic parameters and on the",
        "baseline hazard; not retained.",
        sep = " "
      )
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on the pharmacokinetic parameters and on the baseline hazard; not retained. Chotsiri 2019 Table 1 reports 93/179 (51.9%) male in the PKPD group, i.e. 48.1% female."
    ),
    PARA = list(
      description = "Plasmodium falciparum parasitaemia",
      units       = "parasites/uL",
      type        = "continuous",
      notes       = paste(
        "Screened as a pharmacokinetic covariate and not retained.",
        "Parasitaemia at the time a recurrent infection is detected does",
        "enter the pharmacodynamic ANALYSIS, but as data rather than as a",
        "model covariate: Chotsiri 2019 Eqs. 8-9 back-extrapolate the",
        "censoring interval of the time-to-event likelihood from PAR_OBS,",
        "I_Start = ln(PAR_OBS / 1e4) / K_Growth,slow and I_End =",
        "ln(PAR_OBS / 1e5) / K_Growth,fast, with a 5-fold (slow) and",
        "10-fold (fast) parasite multiplication every 48 h. That interval",
        "censoring shapes the estimation likelihood, not the hazard, so it",
        "has no counterpart in a simulation model. Chotsiri 2019 Table 1",
        "reports median (range) parasitaemia at detection of 48,926",
        "(64-1,496,212) parasites/uL in the PKPD group.",
        sep = " "
      )
    ),
    MAL_NOURISH = list(
      description = "Nutritional status indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on the pharmacokinetic parameters ('All other covariates (parasitaemia, gender, age, and nutritional status) were investigated by a stepwise addition ... approach'); not retained. Children with severe malnutrition were excluded from the parent trial."
    ),
    BODYTEMP = list(
      description = "Axillary body temperature at admission",
      units       = "degC",
      type        = "continuous",
      notes       = "Reported in Chotsiri 2019 Table 1 (median 36.7 degC, range 35.0-39.3 in the PKPD group) but not listed among the covariates tested; not in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 179L,
    n_subjects_pk  = 179L,
    n_subjects_pd  = 741L,
    n_studies      = 1L,
    age_range      = "2.33-58.1 months; median 32.1 months (Table 1, PKPD group). Trial inclusion criterion 3-59 months.",
    weight_range   = "4.20-18.3 kg; median 11.0 kg (Table 1, PKPD group). Allometric centring value 18.0 kg, taken from the typical child of the Tarning 2012 prior model.",
    sex_female_pct = 48.1,
    disease_state  = paste(
      "Apparently healthy children receiving seasonal malaria",
      "chemoprevention during the high-transmission season. Malaria at",
      "enrolment was not an exclusion criterion: 71/179 (39.6%) of the PKPD",
      "group had a positive blood smear at admission, and those children",
      "were treated with artemether-lumefantrine instead of receiving that",
      "month's SMC dose. Children with severe malnutrition or a chronic",
      "condition requiring hospitalisation were excluded.",
      sep = " "
    ),
    dose_range     = paste(
      "Three-day courses of a fixed oral dihydroartemisinin-piperaquine",
      "combination (Duocotexin, 40 mg dihydroartemisinin + 320 mg",
      "piperaquine tetra-phosphate per tablet), once monthly for three",
      "consecutive months, weight-band dosed and rounded to the nearest",
      "quarter tablet, targeting 2 mg/kg/day dihydroartemisinin and",
      "18 mg/kg/day piperaquine phosphate. Table 1 reports a total monthly",
      "piperaquine BASE dose of 29.2 (18.0-39.0) mg/kg in the PKPD group.",
      "Doses were supervised and given on an empty stomach. The",
      "translational simulations use the paediatric tablet (20 mg",
      "dihydroartemisinin + 160 mg piperaquine tetra-phosphate = 85.72 mg",
      "piperaquine base) under the WHO 2010 and WHO 2015 weight bands",
      "(Supplementary Table 2).",
      sep = " "
    ),
    regions        = "Burkina Faso (three rural health facilities in the district of Lena, 40-50 km from Bobo-Dioulasso)",
    notes          = paste(
      "Parent trial NCT00941785, an SMC efficacy trial comparing",
      "sulfadoxine-pyrimethamine + amodiaquine against DHA-PQ; only the",
      "DHA-PQ arm contributed pharmacokinetic data. 466 capillary",
      "(finger-prick) and 71 venous plasma piperaquine concentrations were",
      "assayed by LC-MS/MS (LLOQ 1.50 ng/mL); 1 capillary and 3 venous",
      "samples were below the LLOQ and were omitted. Sampling was sparse",
      "(~3-4 samples per child), so the structural model was stabilised with",
      "a NONMEM $PRIOR frequentist prior from Tarning 2012. Eta shrinkage in",
      "the final model was high: CL/F 54%, Vc/F 46%, Q1/F 74%, Vp2/F 65%,",
      "MTT 76%, F 21%, and epsilon shrinkage 26.8%. The time-to-event model",
      "was developed on the 179 PKPD children and externally validated on",
      "the 562 children of the PD group (741 in total); 432/741 (58.2%) had",
      "a malaria episode during the 4-month study period, with a median time",
      "to malaria of 90 days.",
      sep = " "
    )
  )

  ini({
    # ---- Absorption -------------------------------------------------
    # 'Observed capillary plasma piperaquine concentration-time data were
    # described accurately by a two-transit-compartment absorption model,
    # followed by a three compartment distribution model as reported
    # previously' (Results, 'Pharmacokinetic model'). With n = 2 transit
    # compartments the mean transit time from the dosing compartment to the
    # central compartment spans n + 1 = 3 first-order transfers, so
    # ktr = ka = 3 / MTT (Savic 2007) -- the same parameterisation the
    # sibling Chotsiri_2019_lumefantrine.R uses. That convention is
    # confirmed here by the paper's own secondary parameters: it reproduces
    # the published median Tmax of 2.85 h (3.0 h simulated for the typical
    # 11.0-kg child) whereas ktr = 2 / MTT gives ~4 h.
    lmtt <- log(1.37)  ; label("Mean absorption transit time MTT (h)")                                        # Chotsiri 2019 Table 2: MTT = 1.37 h (bootstrap 95% CI 0.506-1.93; %RSE 26.9)

    # ---- Disposition ------------------------------------------------
    # Three-compartment disposition. Table 2 population estimates are for a
    # typical child of 18.0 kg body weight (Table 2 footnote c).
    #
    # Table 2's unit column is transposed for three rows: it prints
    # 'V_P1/F (L h^-1)', 'Q_2/F (L)' and 'V_P2/F (L h^-1)'. The magnitudes
    # and the disposition structure identify the intended units
    # unambiguously (Vp1 = 274 L, Q2 = 10.8 L/h, Vp2 = 3,490 L), and the
    # resulting terminal half-life of 21.6 days reproduces the published
    # 21.3 days (Supplementary Table 1). See vignette Errata.
    lcl  <- log(7.36)  ; label("Apparent elimination clearance CL/F (L/h, typical 18.0-kg child)")            # Chotsiri 2019 Table 2: CL/F = 7.36 L/h (bootstrap 95% CI 7.52-7.84; %RSE 1.04)
    lvc  <- log(314)   ; label("Apparent central volume of distribution Vc/F (L, typical 18.0-kg child)")     # Chotsiri 2019 Table 2: VC/F = 314 L (bootstrap 95% CI 282-356; %RSE 5.80)
    lq   <- log(9.78)  ; label("Apparent intercompartmental clearance to peripheral1, Q1/F (L/h, typical 18.0-kg child)")  # Chotsiri 2019 Table 2: Q_1/F = 9.78 L/h (bootstrap 95% CI 6.89-12.7; %RSE 15.1)
    lvp  <- log(274)   ; label("Apparent first peripheral volume of distribution Vp1/F (L, typical 18.0-kg child)")        # Chotsiri 2019 Table 2: V_P1/F = 274 (bootstrap 95% CI 266-284; %RSE 1.69)
    lq2  <- log(10.8)  ; label("Apparent intercompartmental clearance to peripheral2, Q2/F (L/h, typical 18.0-kg child)")  # Chotsiri 2019 Table 2: Q_2/F = 10.8 (bootstrap 95% CI 10.5-11.1; %RSE 1.30)
    lvp2 <- log(3490)  ; label("Apparent second peripheral volume of distribution Vp2/F (L, typical 18.0-kg child)")       # Chotsiri 2019 Table 2: V_P2/F = 3,490 (bootstrap 95% CI 3,410-3,580; %RSE 1.30)

    # ---- Relative bioavailability -----------------------------------
    # 'The present study showed a significantly lower exposure to
    # piperaquine compared to the prior model. This was, therefore,
    # corrected for by implementing a categorical covariate on the relative
    # bioavailability' (Results, 'Pharmacokinetic model'). The prior study
    # is the reference category (F = 1) and this SMC cohort is F = 0.726,
    # so a model of THIS population uses 0.726. The paper could not explain
    # the difference; a high-fat-meal effect on piperaquine absorption is
    # offered as one possibility (Discussion).
    lfdepot <- log(0.726) ; label("Relative bioavailability F of this SMC cohort, referenced to the Tarning 2012 prior study (unitless)")  # Chotsiri 2019 Table 2 'Covariates', Relative F = 0.726 (bootstrap 95% CI 0.675-0.781; %RSE 3.71)

    # ---- Allometric scaling -----------------------------------------
    # Methods, 'Pharmacokinetic analysis': 'Body weight (BW) was introduced
    # as an allometric function for all clearance (exponent of n = 0.75) and
    # volume (exponent of n = 1.00) parameters, centralized to 18 kg of body
    # weight according to the typical patient in the prior model (Eq. 1)'.
    # Fixed, not estimated, and reported without uncertainty.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on CL/F, Q1/F and Q2/F (unitless)")     # Chotsiri 2019 Methods, 'Pharmacokinetic analysis': clearance exponent n = 0.75
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on Vc/F, Vp1/F and Vp2/F (unitless)")   # Chotsiri 2019 Methods, 'Pharmacokinetic analysis': volume exponent n = 1.00

    # ---- Capillary / venous matrix conversion -----------------------
    # 'Venous and capillary piperaquine plasma concentrations were modelled
    # simultaneously by estimating a venous-capillary conversion factor at a
    # population level' (Results). The paper's parameter is
    # Conversion_CAP-VEN = 0.380, i.e. venous = 0.380 x capillary, so
    # capillary concentrations are 163% higher than venous ones in the same
    # child (Discussion). The library canonical `cfcap` runs the other way
    # (`Ccap = Cc * cfcap`, capillary:venous), so it takes the reciprocal,
    # 1 / 0.380 = 2.6316. Three independent published pairs confirm the
    # direction and the value: IC50 3.66 ng/mL venous = 9.62 ng/mL
    # capillary (9.62 / 3.66 = 2.63), MIC 12.9 ng/mL venous = 33.9 ng/mL
    # capillary (2.63), and the Supplementary Figure 1 regression of venous
    # on capillary concentrations with slope 0.4269 (95% CI 0.3843-0.4695),
    # which brackets 0.380.
    cfcap <- 1 / 0.380 ; label("Capillary:venous conversion factor for piperaquine plasma concentrations (unitless)")  # Chotsiri 2019 Table 2: Conversion_CAP-VEN = 0.380 (bootstrap 95% CI 0.313-0.450; %RSE 8.99), reported as venous/capillary -- inverted here to match the library's capillary:venous `cfcap`

    # ---- Pharmacodynamics: time to new P. falciparum infection ------
    # Constant baseline hazard (Methods, 'Pharmacodynamic analysis': 'A
    # constant baseline hazard of malaria infection during the
    # high-transmission season was assumed'), reported per year and
    # converted to the model's hourly time base:
    # 6.28 / (365.25 * 24) = 7.1642e-4 /h.
    lbase <- log(6.28 / (365.25 * 24)) ; label("Log baseline hazard of new P. falciparum infection (1/h)")    # Chotsiri 2019 Table 2: BASE = 6.28 per year (bootstrap 95% CI 5.13-11.2; %RSE 9.35). The Discussion quotes a different CI (5.50-14.2) for the same point estimate -- see vignette Errata.

    # Sigmoid Emax protective effect of piperaquine on the hazard
    # (Chotsiri 2019 Eq. 3):
    #   PQ_EFF = 1 - EMAX * CP(t)^gamma / (IC50^gamma + CP(t)^gamma)
    # CP(t) is the predicted VENOUS plasma concentration (Methods,
    # 'Pharmacodynamic analysis'), which is `Cc` in this file.
    lic50 <- log(3.66) ; label("Venous piperaquine concentration giving a 50% reduction of the infection hazard (ng/mL)")  # Chotsiri 2019 Table 2: IC50 = 3.66 ng/mL venous (bootstrap 95% CI 2.09-5.40; %RSE 15.1); equivalent to 9.62 ng/mL capillary (Discussion)
    lhill <- log(1.79) ; label("Shape factor gamma of the sigmoid drug effect on the infection hazard (unitless)")         # Chotsiri 2019 Table 2: gamma = 1.79 (bootstrap 95% CI 1.12-2.45; %RSE 12.5)

    # EMAX is named in Eq. 3 but is absent from Table 2 and from every other
    # table in the paper and its supplement. It is fixed to 1, which the
    # paper proves for itself: Table 2 footnote a defines IC50 as the
    # 'piperaquine venous plasma concentrations associated with a reduction
    # of the baseline hazard by 50%', and in Eq. 3 the hazard at CP = IC50
    # is reduced by EMAX/2, so that definition holds only for EMAX = 1.
    # EMAX = 1 also lets the hazard reach zero at saturating concentrations,
    # which is what the DHA effect of Eq. 4 does explicitly. The sibling
    # Chotsiri_2019_lumefantrine.R makes the same assumption. Recorded in
    # vignette Errata.
    emax <- fixed(1) ; label("Maximum fractional reduction of the infection hazard by piperaquine (unitless)")  # Not tabulated by Chotsiri 2019 -- fixed to 1; proven by the Table 2 footnote a definition of IC50, see vignette Errata

    # ---- Inter-individual variability -------------------------------
    # Table 2's 'Inter-individual variability (%CV)' block prints the
    # variance with the %CV in parentheses, and footnote c gives the
    # conversion 100 * sqrt(exp(omega^2) - 1). Each printed pair satisfies
    # that identity exactly (e.g. 100 * sqrt(exp(0.574) - 1) = 88.1%), so
    # the tabulated numbers are VARIANCES and enter ini() unchanged. IIV was
    # estimated on MTT, CL/F, Vc/F, Q2/F, Vp2/F and F, but not on Q1/F or
    # Vp1/F.
    etalmtt    ~ 0.574    # Chotsiri 2019 Table 2: IIV MTT = 0.574 (88.1% CV; bootstrap 95% CI 0.440-0.827; %RSE 9.35); eta shrinkage 76%
    etalcl     ~ 0.0438   # Chotsiri 2019 Table 2: IIV CL/F = 0.0438 (21.2% CV; bootstrap 95% CI 0.0362-0.0540; %RSE 5.30); eta shrinkage 54%
    etalvc     ~ 0.665    # Chotsiri 2019 Table 2: IIV VC/F = 0.665 (97.2% CV; bootstrap 95% CI 0.0825-1.14; %RSE 18.8); eta shrinkage 46%
    etalq2     ~ 0.0478   # Chotsiri 2019 Table 2: IIV Q_2/F = 0.0478 (22.1% CV; bootstrap 95% CI 0.0444-0.0531; %RSE 2.37)
    etalvp2    ~ 0.0486   # Chotsiri 2019 Table 2: IIV V_P2/F = 0.0486 (22.3% CV; bootstrap 95% CI 0.00000486-0.283; %RSE 65.0); eta shrinkage 65%
    etalfdepot ~ 0.114    # Chotsiri 2019 Table 2: IIV F = 0.114 (34.7% CV; bootstrap 95% CI 0.0805-0.164; %RSE 9.40); eta shrinkage 21%

    # ---- Residual error ---------------------------------------------
    # 'The unexplained residual variability was modelled separately for
    # capillary and venous plasma concentrations and implemented as
    # proportional error models on the log-transformed concentrations,
    # essentially equivalent to proportional errors on an arithmetic scale'
    # (Methods, 'Pharmacokinetic analysis') -- so `prop()` in nlmixr2's
    # linear concentration space. Table 2 footnote a defines both sigma
    # entries as the 'variance of proportional residual error', so the SDs
    # are the square roots: sqrt(0.666) = 0.816 venous and
    # sqrt(0.305) = 0.552 capillary. The venous residual is the larger of
    # the two because only 71 venous samples were available and the
    # capillary-to-venous conversion factor carries no IIV, so all
    # between-subject variability in that ratio lands in the venous
    # residual. See vignette Errata for the variance-vs-SD reading.
    propSd      <- sqrt(0.666) ; label("Proportional residual SD for venous plasma piperaquine (unitless)")     # Chotsiri 2019 Table 2: sigma_VP = 0.666 (bootstrap 95% CI 0.489-0.797; %RSE 11.9), reported as a VARIANCE -> SD = sqrt(0.666)
    propSd_Ccap <- sqrt(0.305) ; label("Proportional residual SD for capillary plasma piperaquine (unitless)")  # Chotsiri 2019 Table 2: sigma_CP = 0.305 (bootstrap 95% CI 0.256-0.346; %RSE 7.76), reported as a VARIANCE -> SD = sqrt(0.305)

    # Placeholder additive residual on the survivor-function output. The
    # source model is an interval-censored time-to-event likelihood, not an
    # additively-observed continuous endpoint; `sur` is exposed as a derived
    # output so the protective-efficacy simulations of the vignette can be
    # run, and this residual exists only so the nlmixr2 likelihood machinery
    # accepts a further endpoint. Same convention as
    # Chotsiri_2019_lumefantrine.R and Lindauer_2017_lacosamide_seizure.R.
    addSd <- 0.001 ; label("Placeholder additive residual error on the malaria-free survival output `sur` (unitless); not from the source -- see vignette Assumptions")
  })

  model({
    # ---- Individual parameters --------------------------------------
    # Chotsiri 2019 Eq. 1: theta_i = theta * exp(eta_i,theta) *
    # (BW_i / 18.0)^n, with n = 0.75 on clearances and n = 1.00 on volumes.
    mtt <- exp(lmtt + etalmtt)
    cl  <- exp(lcl + etalcl)   * (WT / 18.0)^e_wt_cl
    vc  <- exp(lvc + etalvc)   * (WT / 18.0)^e_wt_vc
    q   <- exp(lq)             * (WT / 18.0)^e_wt_cl
    vp  <- exp(lvp)            * (WT / 18.0)^e_wt_vc
    q2  <- exp(lq2 + etalq2)   * (WT / 18.0)^e_wt_cl
    vp2 <- exp(lvp2 + etalvp2) * (WT / 18.0)^e_wt_vc

    # Transit-absorption rate constant. Two transit compartments plus the
    # dosing compartment give three first-order transfers between the dosing
    # site and the central compartment, so MTT = 3 / ktr and ka = ktr
    # (Savic 2007).
    ktr <- 3 / mtt
    ka  <- ktr

    # Three-compartment disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Relative bioavailability of this cohort (0.726 of the prior study).
    f(depot) <- exp(lfdepot + etalfdepot)

    # ---- ODE system -------------------------------------------------
    # Dose enters `depot`, the first compartment of the transit chain, in mg
    # of piperaquine BASE (Chotsiri 2019 Table 1 tabulates the administered
    # dose as 'Total monthly dose of piperaquine base (mg kg^-1)').
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ka  * transit2 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- Observables ------------------------------------------------
    # The disposition parameters are apparent values scaled to the
    # CAPILLARY matrix (the prior model and 466 of the 537 concentrations in
    # this study are capillary finger-prick samples), so 1000 * central / vc
    # is the capillary prediction in ng/mL (mg / L * 1000). The venous
    # prediction is the paper's Conversion_CAP-VEN times that, written here
    # in the library's canonical `Ccap = Cc * cfcap` direction with
    # cfcap = 1 / 0.380. `Cc` is therefore the VENOUS concentration, which
    # is the matrix the pharmacodynamic IC50 is defined on.
    Cc   <- 1000 * central / vc / cfcap
    Ccap <- Cc * cfcap

    # ---- Time to new P. falciparum infection ------------------------
    # Chotsiri 2019 Eqs. 3-6:
    #   PQ_EFF  = 1 - EMAX * CP^gamma / (IC50^gamma + CP^gamma)
    #   DHA_EFF = 0 from 6 days before the first dose of a treatment round
    #             to 24 h after its last dose; 1 otherwise
    #   Hz(t)   = theta_BASE * PQ_EFF * DHA_EFF
    #   S(t)    = exp(-integral_0^t Hz(t) dt)
    # The interval-censoring machinery of Eqs. 8-11 belongs to the
    # estimation likelihood (it turns an observed parasitaemia into a
    # censoring window for the event time) and has no counterpart in a
    # forward simulation, so it is not encoded here -- see vignette Errata.
    hill   <- exp(lhill)
    ic50   <- exp(lic50)
    pqeff  <- 1 - emax * Cc^hill / (ic50^hill + Cc^hill)
    dhaeff <- 1 - CONMED_DHA
    hazard <- exp(lbase) * pqeff * dhaeff

    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)

    # Separate proportional residual errors for the venous and capillary
    # matrices, plus the placeholder residual that exposes the survivor
    # function as a further endpoint (see ini()).
    Cc   ~ prop(propSd)
    Ccap ~ prop(propSd_Ccap)
    sur  ~ add(addSd)
  })
}
