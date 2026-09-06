Wallender_2021_piperaquine <- function() {
  description <- "Population pharmacokinetic model for piperaquine (PPQ) given as intermittent preventive treatment with dihydroartemisinin-piperaquine (DP) to Ugandan children from 2 to 36 months of age, from Wallender 2021. Three-compartment disposition fed by a pre-specified transit absorption chain of two transit compartments. Clearances are allometrically scaled with a fixed exponent of 0.75 and volumes with a fixed exponent of 1, both normalised to the study-median body weight of 8.6 kg, and clearance additionally carries an Emax maturation function of postmenstrual age representing CYP3A4 ontogeny. Relative oral bioavailability falls linearly with weight-for-age z-score (11.3% per z-score unit, centred on the cohort median of -0.5), is multiplied by 0.397 on self-administered rather than directly observed dosing occasions, and carries between-occasion variability of 66.9%. Venous and capillary finger-prick samples were pooled into one joint fit: capillary concentrations are predicted from the venous prediction through the published log-linear conversion Ccap = Cven^0.922 and selected per record by SAMPLE_CAPILLARY, so a single pooled proportional residual serves both matrices exactly as fitted. Sister model files from the same paper: modellib('Wallender_2021_piperaquine_malaria') (incident-malaria repeated-time-to-event hazard) and modellib('Wallender_2021_piperaquine_qtc') (Bazett-corrected QT interval)."
  reference <- paste(
    "Wallender E, Ali AM, Hughes E, Kakuru A, Jagannathan P, Muhindo MK,",
    "Opira B, Whalen M, Huang L, Duvalsaint M, Legac J, Kajubi R, Aweeka F,",
    "Dorsey G, Kamya MR, Rosenthal PJ, Savic RM.",
    "Identifying an optimal dihydroartemisinin-piperaquine dosing regimen",
    "for malaria prevention in young Ugandan children.",
    "Nat Commun. 2021;12(1):6714.",
    "doi:10.1038/s41467-021-27051-8. PMC8602248.",
    "Open Access under CC BY 4.0.",
    "Final parameter estimates are in Table 2; the structural equations are",
    "Eq. 1 (clearance) and Eq. 2 (bioavailability) with the capillary-to-venous",
    "conversion in Table 2 footnote c.",
    "The Peer Review File (41467_2021_27051_MOESM2_ESM.pdf) carries a SUPERSEDED",
    "earlier revision of Table 2 in which clearance is 435 L/d with a power",
    "function of age; those values must not be used.",
    "Sister model files from the same paper:",
    "modellib('Wallender_2021_piperaquine_malaria') and",
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
      notes              = "Time-varying; children were followed monthly from 8 to 156 weeks of age. Allometrically scaled a priori (not estimated) with exponent 0.75 on all clearance parameters and 1 on all volume parameters, normalised to the study-median weight of 8.6 kg (Methods, 'Population PK model'). Body weight was explicitly REMOVED from the covariate screen because the allometry was applied a priori (Peer Review File, response to reviewer comment on line 140).",
      source_name        = "WEIGHT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age).",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "WEEKS, not the register-default months -- Eq. 1 and the maturation half-time of 96 weeks are only meaningful on the week scale. Time-varying. Drives the Emax maturation function PMA / (PMA + 96) on clearance, representing CYP3A4 ontogeny during infancy (piperaquine is primarily metabolised by CYP3A4). Observed range approximately 42 to 196 weeks: gestational age at birth had a median of 39.9 weeks (2.5-97.5 percentile 33.2-41.4) in the every-12-weeks arm and 39.0 weeks (33.0-41.8) in the every-4-weeks arm (Table 1), and postnatal age spanned 8 to 156 weeks.",
      source_name        = "Post menstrual age"
    ),
    WAZ = list(
      description        = "Weight-for-age z-score against the WHO Child Growth Standards.",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TIME-VARYING -- the paper explicitly uses the z-score measured at the time of DP dosing, not a baseline value (Methods lists 'time-varying WAZ' among the screened covariates, and the reviewer request that drove this is answered in the Peer Review File). Enters relative oral bioavailability as the linear-deviation term (1 + 0.113 * (WAZ - (-0.5))), centred on the cohort median WAZ of -0.5 given verbatim in Eq. 2; each 1 z-score DECREASE lowers F by 11.3%. Table 1 reports the WAZ at 8 weeks of age as a median of -0.22 (2.5-97.5 percentile -2.92 to 1.36) in the every-12-weeks arm and -0.31 (-3.75 to 1.21) in the every-4-weeks arm. WAZ was selected over HAZ and WHZ because it gave the largest drop in objective function, the largest effect size, and the best visual predictive check; see covariatesDataExcluded for the two that were screened but not retained.",
      source_name        = "WAZ"
    ),
    SELFADMIN = list(
      description        = "Self-administered dosing-occasion indicator: 1 = the DP course was taken at home without direct observation, 0 = every daily dose of the course was directly observed.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (directly observed therapy)",
      notes              = "PER DOSING OCCASION and time-varying within a subject. Trial design: the first of the three daily DP doses was given in the study clinic and the remaining two were handed to the guardian for home administration, so routine IPT courses have SELFADMIN = 1; all three doses were directly observed for the 32 children in the intensive PK substudy, so those courses have SELFADMIN = 0. The multiplier 0.397 means apparent bioavailability was 60% lower on unobserved courses. The authors read this as an ADHERENCE effect rather than a true absorption effect: 'the magnitude of the association with self-administration is most consistent with a low adherence effect' (Discussion). A mixture model was tested to see whether distinct high- and low-bioavailability subpopulations existed and none could be identified.",
      source_name        = "Self-administered DP"
    ),
    OCC = list(
      description        = "Integer-valued occasion index identifying the DP treatment course, for between-occasion variability on relative oral bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "One occasion per three-day DP course. Between-occasion variability on F was estimated at 66.9% (5.9% RSE, 95% CI 62.3-71.3%) in Table 2 and is the largest single variance component of the PK model. rxode2 has no NONMEM-style '| occ' variability level, so occasions 1 to 25 are multiplexed inside model() by binary indicators onto etaiov_fdepot_1 .. etaiov_fdepot_25 sharing one variance (the equivalent of NONMEM $OMEGA BLOCK(1) SAME). The cap of 25 covers the trial exactly: the every-4-weeks arm received a DP course every 4 weeks from 8 to 104 weeks of age, which is 25 courses. OCC = 0, or any value outside 1..25, zeroes every indicator and yields the occasion-free typical value. Simulating more than 25 courses requires adding further etaiov_fdepot_<n> terms.",
      source_name        = "OCC"
    ),
    SAMPLE_CAPILLARY = list(
      description        = "Per-observation sampling-site indicator: 1 = capillary finger-prick sample, 0 = venous plasma sample.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (venous plasma)",
      notes              = "PER OBSERVATION. The authors first built a model for venous plasma PPQ and then added capillary concentrations to form a joint model, linked by the log-linear relationship ln([PPQ]capillary) = 0.922 * ln([PPQ]venous) (Table 2 footnote c and Supplementary Figure 2, from 70 children with paired samples). The relationship did not vary with age over the study period. This indicator selects which prediction the residual error is applied to; unlike the other registered uses of SAMPLE_CAPILLARY it does NOT switch the residual magnitude, because separate capillary and venous residual errors were tested and the data did not support them (Peer Review File, response to reviewer 3 comment on line 133). Routine study visits used capillary sampling except at 40 and 104 weeks of age; the intensive substudy used venipuncture for the first 24 h and finger-prick thereafter.",
      source_name        = "sample matrix"
    )
  )

  # Anthropometric malnutrition indices the source screened on relative oral
  # bioavailability and reports effect sizes for, but did NOT retain in the
  # final model. Documentation only -- checkModelConventions() does not require
  # these to be referenced in model(). See the vignette Assumptions section.
  covariatesDataExcluded <- list(
    HAZ = list(
      description        = "Height-for-age z-score (stunting axis).",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a time-varying covariate on relative oral bioavailability and found significant (delta-OFV -5.93), with each 1 z-score decrease lowering bioavailability by 4.7%. Not retained: WAZ gave the greatest statistical significance, the greatest effect size and the best visual predictive check, and only one malnutrition index was carried into the final model. Table 1 median height-for-age z-score at 8 weeks: 0.03 (2.5-97.5 percentile -3.27 to 2.06) every 12 weeks, -0.39 (-4.38 to 1.18) every 4 weeks."
    ),
    WHZ = list(
      description        = "Weight-for-height z-score (wasting axis).",
      units              = "unitless (z-score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a time-varying covariate on relative oral bioavailability and found significant (delta-OFV -6.47), with each 1 z-score decrease lowering bioavailability by 4.4%. Not retained, for the same reason as HAZ."
    ),
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened in the stepwise covariate model building (Methods, 'Population PK model'); not retained. Table 1 reports 50% female in the every-12-weeks arm and 47% in the every-4-weeks arm."
    ),
    MUAC = list(
      description        = "Mid-upper arm circumference.",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "NOT COLLECTED in this trial and therefore not screenable, although the authors identify it as a malnutrition marker linked to antimalarial exposure elsewhere in the literature (Peer Review File, response to the reviewer question on MUAC; Discussion). Recorded here so the gap is visible rather than silent."
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "piperaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "piperaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 280L,
    n_studies      = 1L,
    age_range      = "2 to 36 months (dosed from 8 to 104 weeks of age; PK sampled from 12 to 104 weeks)",
    weight_median  = "8.6 kg (study median, used as the allometric reference)",
    sex_female_pct = 48.9,
    disease_state  = "Healthy Ugandan infants and toddlers born to mothers enrolled in a trial of intermittent preventive treatment in pregnancy, receiving intermittent preventive treatment against Plasmodium falciparum malaria in a high-transmission setting.",
    dose_range     = "Weight-band dosed dihydroartemisinin-piperaquine once daily for three consecutive days per course: <6 kg 10/80 mg, 6-<11 kg 20/160 mg, 11-<15 kg 30/240 mg, 15-<20 kg 40/320 mg (manufacturer package insert; Supplementary Table 1). Courses were given every 4 weeks (n = 96) or every 12 weeks (n = 184) from 8 to 104 weeks of age.",
    regions        = "Tororo District, Uganda",
    n_observations = 4573L,
    notes          = "Randomised controlled trial NCT02163447; participant characteristics in Table 1. 4573 PPQ concentrations, of which 578 (12.6%) were below the 0.5 ng/mL lower limit of quantification and were handled by the M6 (BLQ/2) method with BLQ values set to 0.25 ng/mL. All 280 children contributed sparse samples (median 12 per child, range 1-20); a substudy of 32 children contributed intensive profiles at 32 and 104 weeks of age (median 31 per child, range 16-33). None of the children randomised to DP every 4 weeks were born to mothers who received sulfadoxine-pyrimethamine in pregnancy, by trial design. Median birth weight 3000 g (2.5-97.5 percentile 1932-3807) and 2965 g (1694-3688) in the two arms; 11.4% and 14.6% had low birth weight and 7.6% and 12.5% were preterm. Estimation in NONMEM 7.4."
  )

  ini({
    # ------------------------------------------------------------------
    # STRUCTURAL PK -- Wallender 2021 Table 2, "Pharmacokinetic
    # parameters". Every disposition parameter is an APPARENT ORAL value
    # (X/F) for the reference child: body weight 8.6 kg (the study
    # median) and fully mature enzyme activity.
    #
    # WHICH TABLE IS FINAL. The Peer Review File carries an EARLIER
    # revision of Table 2 in which clearance is 435 L/d and age enters as
    # a power function (Age / 57 weeks)^0.269. Reviewer 3 asked for a
    # biologically plausible maturation function; the authors refit with
    # the Emax form and the PUBLISHED Table 2 -- transcribed here --
    # carries the refit estimates. The peer-review numbers are
    # superseded. The two agree where they should: the Emax maturation
    # evaluated at the cohort median postmenstrual age of about 97 weeks
    # gives 867 * 97 / (97 + 96) = 436 L/d, matching the superseded
    # median-centred estimate of 435 L/d.
    # ------------------------------------------------------------------
    lcl  <- log(867)   ; label("Apparent oral clearance at the reference weight, extrapolated to full enzyme maturation (CL/F, L/day)")     # Table 2 Clearance 867 L/d (10.2% RSE, 95% CI 727-1064)
    lvc  <- log(592)   ; label("Apparent central volume of distribution at the reference weight (Vc/F, L)")                                 # Table 2 Volume of central compartment 592 L (10.1% RSE, 95% CI 484-717)
    lq   <- log(511)   ; label("Apparent intercompartmental clearance to the first peripheral compartment at the reference weight (Q1/F, L/day)")  # Table 2 Intercompartmental clearance 1 = 511 L/d (12.4% RSE, 95% CI 401-645)
    lvp  <- log(7240)  ; label("Apparent first peripheral volume of distribution at the reference weight (Vp1/F, L)")                       # Table 2 Volume of peripheral compartment 1 = 7240 L (8.2% RSE, 95% CI 6210-8440)
    lq2  <- log(671)   ; label("Apparent intercompartmental clearance to the second peripheral compartment at the reference weight (Q2/F, L/day)") # Table 2 Intercompartmental clearance 2 = 671 L/d (16.8% RSE, 95% CI 461-934)
    lvp2 <- log(1060)  ; label("Apparent second peripheral volume of distribution at the reference weight (Vp2/F, L)")                      # Table 2 Volume of peripheral compartment 2 = 1060 L (18.2% RSE, 95% CI 731-1510)

    # Mean transit time of the whole absorption chain. Two transit
    # compartments were PRE-SPECIFIED (Table 2 footnote b; the Peer
    # Review File tabulates "Absorption Compartments = 2"), and no
    # separate first-order absorption rate constant is reported, so the
    # chain empties directly into the central compartment. Encoded in
    # the Savic 2007 parameterisation ktr = (n + 1) / MTT with n = 2,
    # i.e. three first-order steps between the dose record and the
    # central compartment, giving ktr = 3 / 0.045 = 66.7 /day. The paper
    # never writes the absorption ODEs; this reading was ratified by the
    # operator (sidecar oare_PMC8602248 request-001 q4, option A) and the
    # sensitivity of the alternative reading is recorded in the vignette.
    lmtt <- log(0.045) ; label("Mean transit time of the absorption chain (MTT, day)")                                                      # Table 2 Absorption transit time 0.045 d (9.1% RSE, 95% CI 0.034-0.048)

    # Allometric exponents applied A PRIORI, not estimated: "Clearance
    # and volume parameters were allometrically scaled for bodyweight a
    # priori by normalizing the child's weight to the median weight of
    # the study population (8.6 kg) and raising to the power of 0.75 for
    # all clearance parameters and to the power of 1 for all volume PK
    # parameters" (Methods, "Population PK model").
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on all clearance parameters, referenced to 8.6 kg (unitless)")       # Methods, "Population PK model"; also Table 2 footnote a
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of body weight on all volume parameters, referenced to 8.6 kg (unitless)")          # Methods, "Population PK model"

    # Emax maturation of clearance on postmenstrual age (Eq. 1 and
    # Table 2 footnote a): CL scales by PMA / (PMA + tm50_cl). The
    # paper labels this parameter "theta Post menstrual age EC50" and
    # defines it as "the PMA when maturation reaches 50%", which is the
    # canonical maturation half-time.
    ltm50_cl <- log(96) ; label("Postmenstrual age at which the maturation of piperaquine clearance is half-maximal (weeks)")               # Table 2 theta Post menstrual age EC50 = 96 weeks (15.4% RSE, 95% CI 73-130)

    # ------------------------------------------------------------------
    # RELATIVE ORAL BIOAVAILABILITY -- Eq. 2 and Table 2 footnote d:
    #   F = 1 * theta_SelfAdministeredDP * (1 + theta_WAZ * (WAZ - (-0.5)))
    #         * exp(BOV)
    # The population value theta_F was ASSUMED to be 100% ("the
    # population bioavailability, which was assumed to be 100%") and is
    # tabulated as "Relative bioavailability (F) = 1" under the
    # pre-specified footnote, so it is an anchor rather than an estimate.
    # ------------------------------------------------------------------
    lfdepot <- fixed(log(1)) ; label("Relative oral bioavailability at the reference covariates and on a directly observed dosing occasion (F, unitless)")  # Table 2 Relative bioavailability (F) = 1, pre-specified anchor; Eq. 2 theta_F

    e_waz_fdepot       <- 0.113 ; label("Fractional change in relative oral bioavailability per unit of weight-for-age z-score above the cohort median of -0.5 (unitless)")  # Table 2 theta Weight for age z-score = 0.113 (18.9% RSE, 95% CI 0.061-0.137)
    e_selfadmin_fdepot <- 0.397 ; label("Multiplicative effect of a self-administered rather than directly observed dosing occasion on relative oral bioavailability (unitless)")  # Table 2 theta Self-administered DP = 0.397 (7.8% RSE, 95% CI 0.344-0.465); 1 - 0.397 = 60% lower

    # ------------------------------------------------------------------
    # CAPILLARY-TO-VENOUS CONVERSION -- Table 2 footnote c:
    #   ln([PPQ]capillary) = theta_slope * ln([PPQ]venous)
    # i.e. Ccap = Cven^0.922. Determined by linear regression on paired
    # samples from 70 children (Supplementary Figure 2); Table 2 reports
    # no RSE or confidence interval for it. The exponent is applied to a
    # concentration expressed in ng/mL, which is the scale on which the
    # regression was performed.
    # ------------------------------------------------------------------
    e_sample_capillary_cc <- 0.922 ; label("Exponent of the log-linear capillary-to-venous piperaquine concentration conversion (unitless)")  # Table 2 theta Capillary to venous conversion = 0.922 (no uncertainty reported)

    # ------------------------------------------------------------------
    # INTERINDIVIDUAL VARIABILITY -- Table 2 "Interindividual
    # variability, %" column, reported as %CV. Converted to the
    # log-normal variance with omega^2 = log(CV^2 + 1).
    # ------------------------------------------------------------------
    # Table 2 Clearance IIV 27.1% (15.1% RSE, 95% CI 22.0-31.0); log(0.271^2 + 1) = 0.0708694
    etalcl  ~ 0.0708694
    # Table 2 Volume of central compartment IIV 32.8% (56.3% RSE, 95% CI 11.8-48.5); log(0.328^2 + 1) = 0.1021811
    etalvc  ~ 0.1021811
    # Table 2 Absorption transit time IIV 43.2% (42.2% RSE, 95% CI 22.4-60.1); log(0.432^2 + 1) = 0.1711123
    etalmtt ~ 0.1711123

    # ------------------------------------------------------------------
    # BETWEEN-OCCASION VARIABILITY ON BIOAVAILABILITY -- Table 2
    # "theta Between occasion variability = 66.9% (5.9% RSE, 95% CI
    # 62.3-71.3%)", entered in Eq. 2 as exp(BOV). This is the largest
    # variance component in the model. log(0.669^2 + 1) = 0.3698801.
    #
    # rxode2 has no NONMEM-style "| occ" variability level, so the single
    # published BOV term is expanded into 25 per-occasion etas selected
    # by binary indicators on the OCC column, all sharing one variance
    # (the equivalent of NONMEM $OMEGA BLOCK(1) SAME). 25 covers the
    # every-4-weeks arm exactly: one DP course every 4 weeks from 8 to
    # 104 weeks of age. Indicator multiplexing makes rxode2 report "some
    # etas defaulted to non-mu referenced" on every solve; that warning
    # concerns estimation-time optimisation only and is inert for
    # simulation.
    # ------------------------------------------------------------------
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

    # ------------------------------------------------------------------
    # RESIDUAL ERROR -- a single proportional error shared by venous and
    # capillary records. Separate residuals for the two matrices were
    # tested and "the data did not support separate residual errors"
    # (Peer Review File, response to reviewer 3 comment on line 133).
    # ------------------------------------------------------------------
    propSd <- 0.446 ; label("Proportional residual error (fraction)")  # Table 2 Proportional error 44.6% (1.7% RSE, 95% CI 43.0-46.0%)
  })

  model({
    # ---- 1. Derived covariate terms ---------------------------------
    # Allometry on the study-median weight of 8.6 kg (Methods).
    allom_cl <- (WT / 8.6)^e_wt_cl
    allom_v  <- (WT / 8.6)^e_wt_vc

    # Eq. 1 maturation term: PMA / (PMA + theta_PMA_EC50), postmenstrual
    # age in weeks. Approaches 1 as maturation completes, so exp(lcl) is
    # the fully mature clearance.
    mat_cl <- PAGE / (PAGE + exp(ltm50_cl))

    # Between-occasion variability on bioavailability. OCC outside
    # 1..25 zeroes every indicator and gives the occasion-free value.
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

    # ---- 2. Individual parameters -----------------------------------
    cl  <- exp(lcl + etalcl) * allom_cl * mat_cl
    vc  <- exp(lvc + etalvc) * allom_v
    q   <- exp(lq)   * allom_cl
    vp  <- exp(lvp)  * allom_v
    q2  <- exp(lq2)  * allom_cl
    vp2 <- exp(lvp2) * allom_v

    mtt <- exp(lmtt + etalmtt)

    # Eq. 2. The self-administration multiplier is raised to the
    # indicator so that a directly observed occasion (SELFADMIN = 0) is
    # unmodified.
    fdepot <- exp(lfdepot + iov_fdepot) *
      (1 + e_waz_fdepot * (WAZ - (-0.5))) *
      e_selfadmin_fdepot^SELFADMIN

    # ---- 3. Micro-constants -----------------------------------------
    # Savic 2007 transit parameterisation: three first-order steps of
    # rate ktr carry the dose from the depot through two transit
    # compartments into the central compartment, so MTT = 3 / ktr.
    ktr <- 3 / mtt

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 4. ODE system ----------------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # ---- 5. Bioavailability -----------------------------------------
    f(depot) <- fdepot

    # ---- 6. Observation and residual error --------------------------
    # Amounts are in mg and volumes in L, so central / vc is mg/L; the
    # factor of 1000 converts to the ng/mL scale on which every
    # piperaquine concentration in the paper is reported and on which
    # the capillary regression and the malaria EC50 are defined.
    Cven <- 1000 * central / vc

    # Table 2 footnote c: ln(Ccap) = 0.922 * ln(Cven). The max() is a
    # numerical guard only -- a solver undershoot to a small negative
    # value would otherwise make the fractional power NaN.
    Ccap <- max(Cven, 0)^e_sample_capillary_cc

    # Prediction for the record as actually drawn. The joint fit pooled
    # venous and capillary observations into one DV column under one
    # proportional residual, so a single endpoint switched by the
    # per-record sampling-site indicator reproduces the published
    # structure exactly.
    Cc <- Cven * (1 - SAMPLE_CAPILLARY) + Ccap * SAMPLE_CAPILLARY

    Cc ~ prop(propSd)
  })
}
