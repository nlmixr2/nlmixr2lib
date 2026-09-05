Wallender_2021_piperaquine_malaria <- function() {
  description <- "Repeated-time-to-event exposure-response model for incident Plasmodium falciparum malaria in Ugandan children receiving intermittent preventive treatment with dihydroartemisinin-piperaquine, from Wallender 2021. The hazard of incident malaria is a constant (exponential) baseline hazard multiplied by a calendar high-transmission-period factor and by a sigmoidal inhibitory effect of the time-varying plasma piperaquine concentration with a maximum effect of complete hazard suppression: h(t) = h0 * transmission * (1 - Cc^gamma / (EC50^gamma + Cc^gamma)) * exp(eta). The three-compartment piperaquine population PK model of the same paper is embedded so the hazard is driven by simulated concentration, matching the source's sequential fit in which the exposure-response model used model-derived individual PK parameters. The cumulative hazard is carried as an ODE state and exposed alongside the survivor function and the fractional protective efficacy that Figure 5B plots. A piperaquine concentration of 15.4 ng/mL reduces the malaria hazard by 95% and was the target concentration for the paper's dosing-regimen simulations. Sister model files from the same paper: modellib('Wallender_2021_piperaquine') (population PK) and modellib('Wallender_2021_piperaquine_qtc') (Bazett-corrected QT interval)."
  reference <- paste(
    "Wallender E, Ali AM, Hughes E, Kakuru A, Jagannathan P, Muhindo MK,",
    "Opira B, Whalen M, Huang L, Duvalsaint M, Legac J, Kajubi R, Aweeka F,",
    "Dorsey G, Kamya MR, Rosenthal PJ, Savic RM.",
    "Identifying an optimal dihydroartemisinin-piperaquine dosing regimen",
    "for malaria prevention in young Ugandan children.",
    "Nat Commun. 2021;12(1):6714.",
    "doi:10.1038/s41467-021-27051-8. PMC8602248.",
    "Open Access under CC BY 4.0.",
    "The hazard equation is Eq. 3 and the final estimates are in Table 2",
    "under 'Pharmacodynamic parameters'; the embedded PK parameters are the",
    "'Pharmacokinetic parameters' rows of the same table.",
    "Sister model files from the same paper:",
    "modellib('Wallender_2021_piperaquine') and",
    "modellib('Wallender_2021_piperaquine_qtc').",
    sep = " "
  )
  vignette <- "Wallender_2021_piperaquine"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Acts on the EMBEDDED PK model only (allometric exponent 0.75 on all clearances and 1 on all volumes, referenced to the study-median 8.6 kg), and thereby on the malaria hazard indirectly through piperaquine concentration. Body weight was screened as a direct predictor of the malaria hazard and was not retained.",
      source_name        = "WEIGHT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age).",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months. Acts on the EMBEDDED PK model only, through the Emax maturation term PMA / (PMA + 96) on clearance. Age was screened as a direct covariate on the malaria hazard and was not retained.",
      source_name        = "Post menstrual age"
    ),
    WAZ = list(
      description        = "Weight-for-age z-score against the WHO Child Growth Standards.",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Acts on the EMBEDDED PK model only, through relative oral bioavailability, and therefore reaches the malaria hazard only through reduced piperaquine exposure. This is the paper's mechanism for the higher predicted malaria incidence in malnourished children: 'Sex, IPT arm, maternal IPT regimen, WAZ, WHZ, and HAZ were not associated with the hazard of incident malaria' -- the association runs entirely through exposure.",
      source_name        = "WAZ"
    ),
    SELFADMIN = list(
      description        = "Self-administered dosing-occasion indicator: 1 = the DP course was taken at home without direct observation, 0 = every daily dose of the course was directly observed.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (directly observed therapy)",
      notes              = "Per dosing occasion. Acts on the EMBEDDED PK model only, multiplying relative oral bioavailability by 0.397. This is the covariate the paper's adherence scenarios manipulate: full adherence uses SELFADMIN = 0 (the bioavailability of the directly observed group) and 1/3 adherence uses SELFADMIN = 1 (the bioavailability observed for non-directly-observed therapy).",
      source_name        = "Self-administered DP"
    ),
    OCC = list(
      description        = "Integer-valued occasion index identifying the DP treatment course, for between-occasion variability on relative oral bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "One occasion per three-day DP course; values 1 to 25 are multiplexed inside model() onto etaiov_fdepot_1 .. etaiov_fdepot_25 sharing one variance. OCC = 0, or any value outside 1..25, yields the occasion-free typical value. See modellib('Wallender_2021_piperaquine') for the full annotation.",
      source_name        = "OCC"
    ),
    TRANSM_HIGH_2015 = list(
      description        = "2015 high-malaria-transmission calendar-period indicator: 1 = the record falls in the 2015 high-transmission period, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled low-transmission periods across all study years)",
      notes              = "TIME-VARYING. High transmission is defined in Methods as 1 March to 31 August annually in Tororo, Uganda. The three annual indicators are mutually exclusive one-hot columns; all three set to 0 selects the low-transmission reference, whose hazard multiplier is 1. The paper included the transmission period as a time-varying categorical covariate coded 0 for low-intensity seasons and 1-3 for the three high-transmission seasons (Peer Review File, response to reviewer 3 comment on line 153); the authors would have preferred a parametric surge or cosine function but could not use one because the periods were fixed to calendar time rather than to age and the peaks differed by year for ecological reasons.",
      source_name        = "Transmission period 2015"
    ),
    TRANSM_HIGH_2016 = list(
      description        = "2016 high-malaria-transmission calendar-period indicator: 1 = the record falls in the 2016 high-transmission period, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled low-transmission periods across all study years)",
      notes              = "Time-varying; same calendar definition and reference as TRANSM_HIGH_2015. The three multipliers span 1.29 to 7.83 because district-wide indoor residual spraying changed insecticide between years (bendiocarb through 2015, pirimiphos-methyl from 2016), so the annual peaks are not exchangeable."
    ),
    TRANSM_HIGH_2017 = list(
      description        = "2017 high-malaria-transmission calendar-period indicator: 1 = the record falls in the 2017 high-transmission period, 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pooled low-transmission periods across all study years)",
      notes              = "Time-varying; same calendar definition and reference as TRANSM_HIGH_2015."
    )
  )

  # Covariates the source screened on the malaria hazard but did not retain,
  # or retained only in an exploratory model. Documentation only.
  covariatesDataExcluded <- list(
    SES_MATERNAL = list(
      description        = "Maternal socioeconomic status, a propensity score summarising household property and income, valued between -1 and 3.",
      units              = "unitless (propensity score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "SIGNIFICANT but NOT RETAINED. In univariate analysis each 1-unit increase in maternal socioeconomic status was associated with a 26.2% decreased risk of malaria (delta-OFV -7.21), but 'when we incorporated SES into the full PK/PD model we encountered unacceptable model instability and confidence intervals could not be reliably acquired by bootstrap, so maternal SES was not included in the final model'. No point estimate is tabulated -- only the 26.2% univariate risk reduction -- so no coefficient can be encoded. The authors note it reduced intraindividual variability but did not modify the key exposure-response relationships or the baseline hazard, and recommend an externally validated socioeconomic measure in future studies (Discussion)."
    ),
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened on the malaria hazard; explicitly reported as not associated with it."
    ),
    HAZ = list(
      description        = "Height-for-age z-score (stunting axis).",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a time-varying covariate on the malaria hazard; explicitly reported as not associated with it. It does reach the hazard indirectly in the PK model, where it was significant on bioavailability but was not the retained malnutrition index -- see modellib('Wallender_2021_piperaquine')."
    ),
    WHZ = list(
      description        = "Weight-for-height z-score (wasting axis).",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a time-varying covariate on the malaria hazard; explicitly reported as not associated with it."
    ),
    IPT_MATERNAL = list(
      description        = "Maternal intermittent preventive treatment regimen received during pregnancy: sulfadoxine-pyrimethamine every 8 weeks, dihydroartemisinin-piperaquine every 8 weeks, or dihydroartemisinin-piperaquine every 4 weeks.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "sulfadoxine-pyrimethamine every 8 weeks",
      notes              = "Screened on the malaria hazard and explicitly reported as not associated with it; also screened on the PK parameters, where 'the maternal chemoprevention regimen was not associated with PK exposure'. Table 1 gives the distribution. Not encoded because no coefficient is reported."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    cumhaz      = list(analyte = "Cumulative hazard of incident malaria", units = NA_character_, specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 280L,
    n_studies      = 1L,
    age_range      = "2 to 36 months (dosed from 8 to 104 weeks of age; malaria surveillance to 156 weeks)",
    weight_median  = "8.6 kg (study median, used as the allometric reference)",
    sex_female_pct = 48.9,
    disease_state  = "Healthy Ugandan infants and toddlers receiving intermittent preventive treatment against Plasmodium falciparum malaria in a high-transmission setting with district-wide indoor residual spraying.",
    dose_range     = "Weight-band dosed dihydroartemisinin-piperaquine once daily for three consecutive days per course, every 4 weeks (n = 96) or every 12 weeks (n = 184) from 8 to 104 weeks of age (Supplementary Table 1).",
    regions        = "Tororo District, Uganda",
    n_observations = 326L,
    biomarkers     = "Repeated time to incident malaria, defined as a temperature above 38 degrees C or a history of fever in the last 24 h together with a positive thick blood smear, more than 14 days after any prior malaria episode. 326 events in 280 children.",
    notes          = "Randomised controlled trial NCT02163447. Malaria incidence was 0.017 episodes per person-year in the every-4-weeks arm versus 0.322 in the every-12-weeks arm (incidence rate ratio 0.05, 95% CI 0.012-0.16), a protective efficacy of 95% (95% CI 84-99%). Exponential, Weibull and Gompertz baseline distributions were tested and the exponential (constant) baseline fit best. Confidence intervals are bootstrap-derived (n = 1000). Two structural limitations the authors state: there was no placebo arm, so the baseline hazard rests on the every-12-weeks arm adjusted for piperaquine exposure and may be underestimated; and the parasiticidal contribution of dihydroartemisinin was not quantified, so cumulative survival was assumed to return to 100% whenever DP was given. Estimation in NONMEM 7.4."
  )

  ini({
    # ==================================================================
    # EMBEDDED PIPERAQUINE PK -- identical to
    # modellib('Wallender_2021_piperaquine'), Table 2 "Pharmacokinetic
    # parameters". See that file for the full annotation, including the
    # note that the Peer Review File carries a SUPERSEDED earlier
    # revision of Table 2 (clearance 435 L/d with a power function of
    # age) that must not be used.
    #
    # The capillary-to-venous conversion of the PK model is a MEASUREMENT
    # model for finger-prick samples and is deliberately absent here: the
    # hazard is driven by plasma concentration, and the paper's
    # exposure-response covariate is "time-varying plasma PPQ
    # concentration as defined by model-derived individual PK parameters"
    # (Methods).
    # ==================================================================
    lcl  <- log(867)   ; label("Apparent oral clearance at the reference weight, extrapolated to full enzyme maturation (CL/F, L/day)")     # Table 2 Clearance 867 L/d
    lvc  <- log(592)   ; label("Apparent central volume of distribution at the reference weight (Vc/F, L)")                                 # Table 2 Volume of central compartment 592 L
    lq   <- log(511)   ; label("Apparent intercompartmental clearance to the first peripheral compartment at the reference weight (Q1/F, L/day)")   # Table 2 Intercompartmental clearance 1 = 511 L/d
    lvp  <- log(7240)  ; label("Apparent first peripheral volume of distribution at the reference weight (Vp1/F, L)")                       # Table 2 Volume of peripheral compartment 1 = 7240 L
    lq2  <- log(671)   ; label("Apparent intercompartmental clearance to the second peripheral compartment at the reference weight (Q2/F, L/day)")  # Table 2 Intercompartmental clearance 2 = 671 L/d
    lvp2 <- log(1060)  ; label("Apparent second peripheral volume of distribution at the reference weight (Vp2/F, L)")                      # Table 2 Volume of peripheral compartment 2 = 1060 L
    lmtt <- log(0.045) ; label("Mean transit time of the absorption chain (MTT, day)")                                                      # Table 2 Absorption transit time 0.045 d; ktr = (n+1)/MTT with n = 2

    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on all clearance parameters, referenced to 8.6 kg (unitless)")       # Methods, "Population PK model"
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of body weight on all volume parameters, referenced to 8.6 kg (unitless)")          # Methods, "Population PK model"

    ltm50_cl <- log(96) ; label("Postmenstrual age at which the maturation of piperaquine clearance is half-maximal (weeks)")               # Table 2 theta Post menstrual age EC50 = 96 weeks

    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability at the reference covariates and on a directly observed dosing occasion (F, unitless)")  # Table 2 Relative bioavailability (F) = 1, pre-specified anchor

    e_waz_fdepot       <- 0.113 ; label("Fractional change in relative oral bioavailability per unit of weight-for-age z-score above the cohort median of -0.5 (unitless)")  # Table 2 theta Weight for age z-score = 0.113
    e_selfadmin_fdepot <- 0.397 ; label("Multiplicative effect of a self-administered rather than directly observed dosing occasion on relative oral bioavailability (unitless)")  # Table 2 theta Self-administered DP = 0.397

    # Table 2 Clearance IIV 27.1%; log(0.271^2 + 1) = 0.0708694
    etalcl  ~ 0.0708694
    # Table 2 Volume of central compartment IIV 32.8%; log(0.328^2 + 1) = 0.1021811
    etalvc  ~ 0.1021811
    # Table 2 Absorption transit time IIV 43.2%; log(0.432^2 + 1) = 0.1711123
    etalmtt ~ 0.1711123

    # Table 2 theta Between occasion variability 66.9%; log(0.669^2 + 1) = 0.3698801.
    # Expanded into 25 per-occasion etas sharing one variance because rxode2
    # has no NONMEM-style "| occ" level; see the sister PK model file.
    etaiov_fdepot_1  ~ 0.3698801
    etaiov_fdepot_2  ~ fixed(0.3698801)
    etaiov_fdepot_3  ~ fixed(0.3698801)
    etaiov_fdepot_4  ~ fixed(0.3698801)
    etaiov_fdepot_5  ~ fixed(0.3698801)
    etaiov_fdepot_6  ~ fixed(0.3698801)
    etaiov_fdepot_7  ~ fixed(0.3698801)
    etaiov_fdepot_8  ~ fixed(0.3698801)
    etaiov_fdepot_9  ~ fixed(0.3698801)
    etaiov_fdepot_10 ~ fixed(0.3698801)
    etaiov_fdepot_11 ~ fixed(0.3698801)
    etaiov_fdepot_12 ~ fixed(0.3698801)
    etaiov_fdepot_13 ~ fixed(0.3698801)
    etaiov_fdepot_14 ~ fixed(0.3698801)
    etaiov_fdepot_15 ~ fixed(0.3698801)
    etaiov_fdepot_16 ~ fixed(0.3698801)
    etaiov_fdepot_17 ~ fixed(0.3698801)
    etaiov_fdepot_18 ~ fixed(0.3698801)
    etaiov_fdepot_19 ~ fixed(0.3698801)
    etaiov_fdepot_20 ~ fixed(0.3698801)
    etaiov_fdepot_21 ~ fixed(0.3698801)
    etaiov_fdepot_22 ~ fixed(0.3698801)
    etaiov_fdepot_23 ~ fixed(0.3698801)
    etaiov_fdepot_24 ~ fixed(0.3698801)
    etaiov_fdepot_25 ~ fixed(0.3698801)

    # ==================================================================
    # MALARIA HAZARD -- Eq. 3 and Table 2 "Pharmacodynamic parameters".
    #
    #   h(t) = theta_Baseline/1000 * theta_TransmissionPeriod
    #            * (1 - [PPQ]^gamma / (theta_EC50^gamma + [PPQ]^gamma))
    #            * exp(eta)
    #
    # SOURCE CONFLICT, ADJUDICATED. Table 2 footnote e abbreviates the
    # survival function WITHOUT the leading "1 -". Eq. 3 in the body
    # carries it. Eq. 3 is correct: without the "1 -" the hazard would
    # RISE with drug concentration, inverting the paper's entire result.
    # Per the standing policy (printed equation outranks a table
    # footnote), Eq. 3 is used, and the reading is confirmed
    # arithmetically -- 15.4^3.13 / (6.00^3.13 + 15.4^3.13) = 0.9503,
    # reproducing the paper's headline claim that 15.4 ng/mL reduces the
    # malaria hazard by 95%.
    #
    # The maximum effect is complete suppression: "Higher PPQ
    # concentrations reduced the risk of malaria with a maximum effect of
    # 100% reduction in hazard", and Eq. 3 accordingly carries no Imax
    # term, so none is introduced here.
    # ==================================================================
    # Table 2 tabulates the baseline hazard already divided by 1000
    # ("Baseline hazard/1000 (per day) = 0.402"), and footnote e restates
    # the /1000 inside the survival function, so the hazard on the
    # natural scale is 0.402/1000 = 4.02e-4 per day = 0.147 per year.
    lhaz_base <- log(0.402 / 1000) ; label("Baseline hazard of incident malaria during low-transmission periods (1/day)")  # Table 2 Baseline hazard/1000 (per day) = 0.402 (17.2% RSE, 95% CI 0.268-0.536)

    e_transm_high_2015_haz_base <- 1.29 ; label("Multiplicative adjustment of the baseline malaria hazard during the 2015 high-transmission period (unitless)")  # Table 2 theta Transmission period 2015 = 1.29 (33.9% RSE, 95% CI 0.65-2.40)
    e_transm_high_2016_haz_base <- 5.20 ; label("Multiplicative adjustment of the baseline malaria hazard during the 2016 high-transmission period (unitless)")  # Table 2 theta Transmission period 2016 = 5.20 (16.1% RSE, 95% CI 3.92-7.09)
    e_transm_high_2017_haz_base <- 7.83 ; label("Multiplicative adjustment of the baseline malaria hazard during the 2017 high-transmission period (unitless)")  # Table 2 theta Transmission period 2017 = 7.83 (16.7% RSE, 95% CI 5.75-10.96)

    lec50 <- log(6.00) ; label("Piperaquine concentration producing a half-maximal reduction of the incident-malaria hazard (EC50, ng/mL)")  # Table 2 theta PPQ EC50 = 6.00 ng/mL (14.4% RSE, 95% CI 4.25-7.58)
    lhill <- log(3.13) ; label("Hill coefficient of the piperaquine hazard-reduction sigmoid (unitless)")                                    # Table 2 theta PPQ gamma = 3.13 (24.3% RSE, 95% CI 1.96-4.89)

    # Table 2 baseline-hazard interindividual variability 69.6% (66.1% RSE,
    # 95% CI 43.3-146); log(0.696^2 + 1) = 0.3950214
    etalhaz_base ~ 0.3950214

    # The source fits this layer with a repeated-time-to-event likelihood on
    # event records, so there is no observation-error model to translate.
    # This placeholder additive residual is attached to the survivor-function
    # output so the nlmixr2 likelihood machinery accepts the model for forward
    # simulation. It is NOT from the source -- see the vignette Assumptions
    # and deviations section.
    addSd <- fixed(0.001) ; label("Additive residual standard deviation on the survivor-function output sur (unitless); not paper-derived")
  })

  model({
    # ---- embedded PK ------------------------------------------------
    allom_cl <- (WT / 8.6)^e_wt_cl
    allom_v  <- (WT / 8.6)^e_wt_vc
    mat_cl  <- PAGE / (PAGE + exp(ltm50_cl))

    iov_fdepot <-
      (OCC ==  1) * etaiov_fdepot_1  + (OCC ==  2) * etaiov_fdepot_2  +
      (OCC ==  3) * etaiov_fdepot_3  + (OCC ==  4) * etaiov_fdepot_4  +
      (OCC ==  5) * etaiov_fdepot_5  + (OCC ==  6) * etaiov_fdepot_6  +
      (OCC ==  7) * etaiov_fdepot_7  + (OCC ==  8) * etaiov_fdepot_8  +
      (OCC ==  9) * etaiov_fdepot_9  + (OCC == 10) * etaiov_fdepot_10 +
      (OCC == 11) * etaiov_fdepot_11 + (OCC == 12) * etaiov_fdepot_12 +
      (OCC == 13) * etaiov_fdepot_13 + (OCC == 14) * etaiov_fdepot_14 +
      (OCC == 15) * etaiov_fdepot_15 + (OCC == 16) * etaiov_fdepot_16 +
      (OCC == 17) * etaiov_fdepot_17 + (OCC == 18) * etaiov_fdepot_18 +
      (OCC == 19) * etaiov_fdepot_19 + (OCC == 20) * etaiov_fdepot_20 +
      (OCC == 21) * etaiov_fdepot_21 + (OCC == 22) * etaiov_fdepot_22 +
      (OCC == 23) * etaiov_fdepot_23 + (OCC == 24) * etaiov_fdepot_24 +
      (OCC == 25) * etaiov_fdepot_25

    cl  <- exp(lcl + etalcl) * allom_cl * mat_cl
    vc  <- exp(lvc + etalvc) * allom_v
    q   <- exp(lq)   * allom_cl
    vp  <- exp(lvp)  * allom_v
    q2  <- exp(lq2)  * allom_cl
    vp2 <- exp(lvp2) * allom_v

    mtt <- exp(lmtt + etalmtt)

    fdepot <- exp(lfdepot + iov_fdepot) *
      (1 + e_waz_fdepot * (WAZ - (-0.5))) *
      e_selfadmin_fdepot^SELFADMIN

    ktr <- 3 / mtt
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    f(depot) <- fdepot

    # Plasma piperaquine concentration in ng/mL (amounts in mg, volumes
    # in L, so the factor of 1000 converts mg/L to ng/mL).
    Cc <- 1000 * central / vc

    # ---- malaria hazard ---------------------------------------------
    # One-hot transmission-period adjustment. All three indicators at 0
    # selects the pooled low-transmission reference, whose multiplier is 1.
    transm <- 1 +
      (e_transm_high_2015_haz_base - 1) * TRANSM_HIGH_2015 +
      (e_transm_high_2016_haz_base - 1) * TRANSM_HIGH_2016 +
      (e_transm_high_2017_haz_base - 1) * TRANSM_HIGH_2017

    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Numerical guard only: a solver undershoot to a small negative
    # concentration would make the fractional power NaN.
    cppq <- max(Cc, 0)

    # Fractional reduction of the malaria hazard attributable to
    # piperaquine. This is the quantity plotted in Figure 5B as
    # "protective efficacy", and it reaches 0.95 at 15.4 ng/mL.
    protection <- cppq^hill / (ec50^hill + cppq^hill)

    hazard <- exp(lhaz_base + etalhaz_base) * transm * (1 - protection)

    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0

    # Probability of remaining free of incident malaria since the start
    # of the integration window. The source assumed cumulative survival
    # returned to 100% each time a DP course was given, to stand in for
    # the unquantified parasiticidal effect of dihydroartemisinin; that
    # reset is a property of the simulation design rather than of the
    # hazard, so it is applied in the vignette per treatment course
    # rather than hard-coded here.
    sur <- exp(-cumhaz)

    sur ~ add(addSd)
  })
}
