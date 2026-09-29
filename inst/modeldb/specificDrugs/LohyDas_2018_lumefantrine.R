LohyDas_2018_lumefantrine <- function() {
  description <- paste(
    "Population PK model of lumefantrine (LF) in 22 Rwandese pregnant women in",
    "their second or third trimester treated for uncomplicated Plasmodium",
    "falciparum malaria with the standard fixed-dose oral artemether-lumefantrine",
    "combination (80/480 mg twice daily for 3 days). Absorption is described by a",
    "5-transit-compartment chain (n = 5 fixed) followed by a two-compartment",
    "disposition; the mean transit time is 4.04 h and the terminal half-life is",
    "about 4 days. Relative bioavailability F is fixed at 1 and carries a large",
    "between-occasion variability whose distribution is Box-Cox transformed",
    "(shape -0.605); mean transit time carries both between-subject and",
    "between-occasion variability. Allometric body-weight scaling is applied to",
    "all clearances (exponent 0.75) and volumes (exponent 1.0), centered on the",
    "cohort-median 59 kg. Concentrations and doses are on a molar basis",
    "(nmol/L, nmol), matching the paper's log-molar estimation scale. The",
    "companion artemether-dihydroartemisinin model from the same paper is",
    "LohyDas_2018_artemether; the authors modelled the two drugs separately",
    "because a simultaneous fit was unstable."
  )
  reference <- paste(
    "Lohy Das J, Rulisa S, de Vries PJ, Mens PF, Kaligirwa N, Agaba S,",
    "Tarning J, Karlsson MO, Dorlo TPC (2018).",
    "Population pharmacokinetics of artemether, dihydroartemisinin, and",
    "lumefantrine in Rwandese pregnant women treated for uncomplicated",
    "Plasmodium falciparum malaria.",
    "Antimicrobial Agents and Chemotherapy 62(10):e00518-18.",
    "doi:10.1128/AAC.00518-18.",
    sep = " "
  )
  vignette <- "LohyDas_2018_artemether_lumefantrine"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Units are derived from the `units` block above.
  compartmentData <- list(
    depot = list(analyte = "lumefantrine", units = "nmol", specimen = "administration site", verified = FALSE),
    transit1 = list(analyte = "lumefantrine", units = "nmol", specimen = "administration site", verified = FALSE),
    transit2 = list(analyte = "lumefantrine", units = "nmol", specimen = "administration site", verified = FALSE),
    transit3 = list(analyte = "lumefantrine", units = "nmol", specimen = "administration site", verified = FALSE),
    transit4 = list(analyte = "lumefantrine", units = "nmol", specimen = "administration site", verified = FALSE),
    transit5 = list(analyte = "lumefantrine", units = "nmol", specimen = "administration site", verified = FALSE),
    central = list(analyte = "lumefantrine", units = "nmol", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "lumefantrine", units = "nmol", specimen = "tissue", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at enrollment",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed at enrollment. Lohy Das 2018 Results ('Covariates'):",
        "'Considering biological plausibility and previous population PK",
        "reports, total body weight, centered to median body weight, was",
        "implemented allometrically on clearances (raised to the power of",
        "0.75) and volumes of distribution (raised to the power of 1.0).'",
        "Methods gives the centering value explicitly: 'All size descriptors",
        "were scaled to their respective medians (i.e., total body weight",
        "[59 kg], ...)'. Ideal body weight, fat-free mass and normal fat mass",
        "were also explored allometrically (Methods Eqs. 7-9) but none",
        "improved the fit, so only total body weight is retained here.",
        "Cohort weight range 40.0-65.0 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    OCC = list(
      description = "Dose occasion index, 1 to 6 across the six-dose artemether-lumefantrine regimen",
      units = "(count)",
      type = "categorical",
      reference_category = "1 (the first dose, at hour 0)",
      notes = paste(
        "Integer occasion column taking value k on the interval starting at",
        "the k-th dose. Doses are given at 0, 8, 24, 36, 48 and 60 h (Lohy Das",
        "2018 Methods, 'Study design': 'Enrolled patients were prescribed 4",
        "tablets of a fixed oral combination of ARM and LF twice daily under",
        "supervision for 3 days (at 0 h [initial dose] and 8, 24, 36, 48, and",
        "60 h)'), so OCC = 1 for 0 <= t < 8 h, 2 for 8 <= t < 24 h, 3 for",
        "24 <= t < 36 h, 4 for 36 <= t < 48 h, 5 for 48 <= t < 60 h and 6 for",
        "t >= 60 h. Methods, 'Population pharmacokinetics' defines the",
        "occasion: 'dosing occasion (OCC; i.e., each dose given was considered",
        "single dosing occasion)'. The column serves purely as the grouping",
        "for the between-occasion variability on relative bioavailability and",
        "on mean transit time; unlike the sibling Ding_2026_artemether model",
        "there is NO fixed time-dependent clearance term on OCC in this paper",
        "(dosing occasion on F was selected in the forward step but was 'not",
        "maintained in the backward elimination step', Results 'Covariates').",
        "For simulations of the alternative 5-day regimen (10 doses, Results",
        "'Model-based simulations') the paper does not state how occasions",
        "beyond 6 were handled; holding OCC at 6 thereafter reuses the",
        "occasion-6 random effects and is the encoding used in the vignette."
      ),
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    EGA = list(
      description = "Estimated gestational age at enrollment",
      units = "weeks",
      type = "continuous",
      notes = paste(
        "Screened on relative bioavailability (described linearly and with a",
        "spline) and additionally as a categorical trimester-2-versus-3",
        "contrast, and selected during forward selection (P <= 0.05), but not",
        "retained in backward elimination (P <= 0.01): Lohy Das 2018 Results",
        "'Covariates'. The paper reports no point estimate for the effect, so",
        "it cannot be encoded. Cohort EGA range 15.7-39.0 weeks (Table 1)."
      )
    ),
    PARASITEMIA = list(
      description = "Asexual Plasmodium falciparum parasite density (baseline and time-varying)",
      units = "parasites/uL",
      type = "continuous",
      notes = paste(
        "Baseline and time-varying parasite density were screened on mean",
        "transit time and on relative bioavailability. Observed parasitemia",
        "density on MTT was selected during forward selection but was 'not",
        "maintained in the backward elimination step' (Lohy Das 2018 Results",
        "'Covariates'). No point estimate is published. Cohort baseline range",
        "3,060-160,000 parasites/uL (Table 1)."
      )
    ),
    TEMP = list(
      description = "Body temperature at enrollment",
      units = "degrees Celsius",
      type = "continuous",
      notes = paste(
        "Listed among the covariates considered for exploration (Lohy Das 2018",
        "Methods, 'Population pharmacokinetics') but not reported as selected",
        "at any step and carrying no published point estimate. Cohort range",
        "34.9-38.6 degrees Celsius (Table 1)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 22L,
    n_studies = 1L,
    age_range = "18-39 years (Table 1 median 26)",
    weight_range = "40.0-65.0 kg (Table 1 median 59.0)",
    weight_typical = "59 kg (population median; allometric reference, Methods)",
    sex_female_pct = 100,
    race_ethnicity = "East African (Rwandese)",
    disease_state = paste(
      "Uncomplicated Plasmodium falciparum malaria, confirmed by light",
      "microscopy, during pregnancy. Eleven patients were in the second",
      "trimester (EGA range 15.7-27.6 weeks) and eleven in the third",
      "(EGA range 28.3-39.0 weeks); median EGA 27.9 weeks. Median baseline",
      "parasitemia 24,970 parasites/uL (range 3,060-160,000). Median body",
      "mass index 21.9 kg/m2 (range 15.6-25.4). Eligibility: pregnant women",
      "older than 18 years with microscopy-confirmed uncomplicated P.",
      "falciparum malaria; pregnancy confirmed by hCG urine test and",
      "gestational age estimated by ultrasound. All patients cleared",
      "microscopic parasitemia by day 3; one recrudescence (4.55%)."
    ),
    dose_range = paste(
      "Fixed oral combination of artemether and lumefantrine (Coartem,",
      "Novartis; 20 mg artemether / 120 mg lumefantrine per tablet), 4 tablets",
      "(80 mg artemether / 480 mg lumefantrine) twice daily for 3 days, given",
      "at 0, 8, 24, 36, 48 and 60 h under supervision with a glass of milk",
      "and/or a small cake (a fatty meal) to enhance absorption. Each 480 mg",
      "lumefantrine dose is 907,475 nmol using the lumefantrine molar mass",
      "528.94 g/mol."
    ),
    sampling = paste(
      "Venous plasma. A pre-dose sample, then samples at 2 and 4 h after each",
      "dose and immediately before doses 2-6 (troughs), with additional",
      "samples at 0.25 h after doses 1 and 2, at 6, 8 and 12 h after the last",
      "dose, and a scheduled day 7 sample. 363 lumefantrine samples were",
      "analyzed (9 excluded for hemolysis). Lumefantrine was measured by",
      "automated solid-phase extraction with LC-MS/MS; LLOQ 24.86 ng/mL",
      "(47 nmol/L). Fewer than 8% of lumefantrine observations were below the",
      "LLOQ and these were excluded from the fit (FOCEI estimation)."
    ),
    regions = paste(
      "Rwanda (obstetrics and gynecology ward, Rwamagana district hospital,",
      "eastern province; mesoendemic transmission). Study conducted June 2007",
      "to July 2009; Rwanda National Ethics Committee study IRB 00001497. The",
      "PK study was nested within a pharmacovigilance study of ACT use in",
      "pregnancy."
    ),
    notes = paste(
      "NONMEM 7.3, FOCEI. Because the lumefantrine data were sparse in the",
      "elimination phase, frequentist informative priors from a previous",
      "Ugandan pregnant/non-pregnant lumefantrine study (Kloprogge 2013, the",
      "paper's reference 15) were applied to the lumefantrine parameter",
      "estimates, which is what allowed a peripheral compartment to be",
      "identified; the prior for Q/F was recalculated to represent pregnant",
      "women. The values encoded here are the paper's FINAL posterior",
      "estimates in Table 2, not the priors. Parameter imprecision was",
      "obtained by sampling importance resampling (SIR). Eta shrinkage was",
      "25.8% and 11.1% for the two between-subject random effects and 14.5-39.4%",
      "(MTT) and 32.5-65.3% (F) for the between-occasion random effects;",
      "epsilon shrinkage was 19.9%."
    )
  )

  ini({
    # Structural PK parameters, Lohy Das 2018 Table 2 (Lumefantrine block).
    # Concentrations were modelled on a natural-log MOLAR scale (Methods,
    # 'Population pharmacokinetics': 'The molar units of LF, ARM, and DHA
    # concentration were transformed to their natural logarithms for this
    # modeling analysis'), so doses are supplied in nmol and Cc is nmol/L.
    # Typical values are apparent (CL/F, V/F) at the cohort median WT = 59 kg.

    lfdepot <- fixed(log(1))
    label("Reference relative oral bioavailability of lumefantrine, F (unitless) (the source paper fixes F at 1)")  # Lohy Das 2018 Table 2: F = '1 fixed'; Methods: 'The typical relative F was implemented as a fixed parameter for the parent analyte, i.e., LF and ARM (100% relative bioavailability).'

    lmtt <- log(4.04)
    label("Mean transit time of the 5-compartment transit-absorption chain, MTT (h)")  # Lohy Das 2018 Table 2: MTT = 4.04 h (%RSE 5.16, 90% CI 3.71-4.41)

    lcl <- log(4.49)
    label("Apparent lumefantrine elimination clearance, CL/F at WT = 59 kg (L/h)")  # Lohy Das 2018 Table 2: CL/F = 4.49 L/h (%RSE 6.59, 90% CI 4.18-5.17)

    lvc <- log(139)
    label("Apparent lumefantrine central volume of distribution, Vc/F at WT = 59 kg (L)")  # Lohy Das 2018 Table 2: Vc/F = 139 L (%RSE 6.77, 90% CI 119-149)

    lq <- log(0.924)
    label("Apparent lumefantrine intercompartmental clearance, Q/F at WT = 59 kg (L/h)")  # Lohy Das 2018 Table 2: Q/F = 0.924 L/h (%RSE 13.3, 90% CI 0.770-1.21)

    lvp <- log(111)
    label("Apparent lumefantrine peripheral volume of distribution, Vp/F at WT = 59 kg (L)")  # Lohy Das 2018 Table 2: Vp/F = 111 L (%RSE 8.69, 90% CI 96.5-129)

    # Box-Cox shape for the relative-bioavailability random effect. Applied in
    # model() in the Petersson (2009) form documented for the canonical
    # boxcox_<param> family. Estimated (not fixed) by the source paper.
    boxcox_lfdepot <- -0.605
    label("Box-Cox shape parameter (lambda) for the relative-bioavailability random-effect distribution, unitless")  # Lohy Das 2018 Table 2: 'Box-Cox shape parameter for BSV on F' = -0.605 (%RSE 34.9, 90% CI -0.590 to -0.180)

    # Allometric exponents. Lohy Das 2018 Methods: 'All size descriptors were
    # scaled to their respective medians (i.e., total body weight [59 kg] ...)
    # on PK parameters using allometric power exponents of 0.75 for clearances
    # (CL/F, Q/F, CL_ARM/F, and CL_DHA/F) and 1 for volumes of distribution
    # (Vc/F, Vp/F, V_ARM/F, and V_DHA/F).' The exponents are structural choices
    # imposed by the authors rather than estimated quantities (no RSE or CI is
    # given for either), so both are encoded as fixed.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on clearances (CL/F, Q/F), unitless")  # Lohy Das 2018 Methods: allometric power exponent 0.75 for clearances
    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on volumes (Vc/F, Vp/F), unitless")  # Lohy Das 2018 Methods: allometric power exponent 1 for volumes of distribution

    # Between-subject variability. Lohy Das 2018 Table 2 footnote a states the
    # scale explicitly: 'Coefficient of variation (CV) for BSV and BOV was
    # calculated as 100 x (variance)^1/2', i.e. the printed %CV is 100 times
    # the eta-scale SD, so variance = (CV/100)^2. BSV was NOT estimated for
    # CL/F, Q/F or Vp/F (Results: these 'were estimated with poor precision
    # (>50% residual standard error)'), so those parameters carry no eta.
    etalmtt ~ 1.7424   # Table 2 row 'MTT (h)', BSV = 132% CV (%RSE 37.9, 90% CI 72.6-178); variance = 1.32^2
    etalvc ~ 0.2372    # Table 2 row 'Vc/F (liters)', BSV = 48.7% CV (%RSE 56.8, 90% CI 17.8-77.8); variance = 0.487^2

    # Between-occasion variability, one eta per dosing occasion (six doses).
    # Results: 'The addition of between-occasion variability (BOV) was
    # significant for mean transit time (MTT) (change in objective function
    # value [dOFV] = -15.1) and relative bioavailability (F) (dOFV = -124.2).'
    # A single BOV variance is shared across occasions in NONMEM, so only the
    # first slot is estimable here and the remaining five are fixed to it.
    etaiov_fdepot_1 ~ 2.0736
    # Table 2 row 'F' (Lumefantrine), BOV = 144% CV (%RSE 19.7, 90% CI 106-189); variance = 1.44^2
    etaiov_fdepot_2 ~ fixed(2.0736)
    etaiov_fdepot_3 ~ fixed(2.0736)
    etaiov_fdepot_4 ~ fixed(2.0736)
    etaiov_fdepot_5 ~ fixed(2.0736)
    etaiov_fdepot_6 ~ fixed(2.0736)

    etaiov_mtt_1 ~ 0.2116
    # Table 2 row 'MTT (h)' second line, BOV = 46.0% CV (%RSE 43.6, 90% CI 13.1-64.5); variance = 0.46^2
    etaiov_mtt_2 ~ fixed(0.2116)
    etaiov_mtt_3 ~ fixed(0.2116)
    etaiov_mtt_4 ~ fixed(0.2116)
    etaiov_mtt_5 ~ fixed(0.2116)
    etaiov_mtt_6 ~ fixed(0.2116)

    # Residual error. Lohy Das 2018 Methods, 'Population pharmacokinetics':
    # 'The unexplained residual error was estimated using an additive error
    # model on the logarithmic scale for all drugs, which equates to an
    # exponential error model on an arithmetic scale.' By the standing
    # nlmixr2lib convention an additive-on-log-scale residual maps to a
    # proportional residual in linear space, with propSd equal to the
    # log-scale SD. Table 2 reports it as a percentage.
    propSd <- 0.487
    label("Proportional residual SD for lumefantrine plasma concentration (SD on the log scale, approximately CV in linear space)")  # Lohy Das 2018 Table 2: RUV = 48.7% (%RSE 4.82, 90% CI 45.8-53.5)
  })

  model({
    # Dose-occasion indicators. OCC is 1-6 across the six-dose regimen; see the
    # covariateData notes for the mapping from time to occasion.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 +
      oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4 +
      oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6

    # Individual PK parameters. Allometric body-weight scaling on both apparent
    # clearances (exponent 0.75) and both apparent volumes (exponent 1),
    # centered on the cohort median WT = 59 kg (Lohy Das 2018 Methods; Results
    # 'Covariates'). CL/F, Q/F and Vp/F carry no between-subject random effect.
    cl <- exp(lcl) * (WT / 59)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 59)^e_wt_vc
    q <- exp(lq) * (WT / 59)^e_wt_cl
    vp <- exp(lvp) * (WT / 59)^e_wt_vc

    # Mean transit time carries both between-subject and between-occasion
    # random effects. Transit-chain rate constant in the Savic (2007)
    # parameterisation, which the paper states explicitly in the Fig. 1A
    # legend: 'ktr, transit absorption rate constant [ktr = (n + 1)/mean
    # transit time]'. With n = 5 transit compartments (Results: 'an absorption
    # model consisting of 5 first-order transit compartments delivering the
    # absorbed amount to the central compartment') the chain depot -> transit1
    # ... -> transit5 -> central has n + 1 = 6 equal-rate transfers, so
    # ktr = 6 / MTT.
    mtt <- exp(lmtt + etalmtt + iov_mtt)
    ktr <- 6 / mtt

    # Disposition micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODE system: oral depot -> 5 transit compartments -> two-compartment
    # lumefantrine disposition (Lohy Das 2018 Fig. 1A).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(transit4) <- ktr * transit3 - ktr * transit4
    d/dt(transit5) <- ktr * transit4 - ktr * transit5
    d/dt(central) <- ktr * transit5 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Relative bioavailability. The typical value is fixed at 1; the
    # occasion-level random effect is Box-Cox transformed before exponentiation
    # (Methods Eq. 3 and Table 2 'Box-Cox shape parameter'), written out here
    # in the Petersson (2009) form because rxode2's boxCox() attaches to the
    # residual-error model and cannot transform an eta. Note that with a
    # negative shape the distribution is bounded above (F < exp(1/0.605) = 5.2)
    # but has a long lower tail, so occasional occasions absorb very little.
    phi_fdepot <- (exp(iov_fdepot)^boxcox_lfdepot - 1) / boxcox_lfdepot
    f(depot) <- exp(lfdepot + phi_fdepot)

    # Plasma concentration in nmol/L (dose in nmol, volume in L).
    Cc <- central / vc

    # NONMEM additive-on-log-scale residual maps to nlmixr2 proportional.
    Cc ~ prop(propSd)
  })
}
