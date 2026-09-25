LohyDas_2018_artemether <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model of artemether (ARM) and its",
    "active metabolite dihydroartemisinin (DHA) in 22 Rwandese pregnant women",
    "in their second or third trimester treated for uncomplicated Plasmodium",
    "falciparum malaria with the standard fixed-dose oral",
    "artemether-lumefantrine combination (80/480 mg twice daily for 3 days).",
    "Artemether absorption is described by a 2-transit-compartment chain",
    "(n = 2 fixed, mean transit time 0.738 h) followed by one-compartment",
    "artemether disposition; complete and irreversible in-vivo conversion of",
    "artemether to dihydroartemisinin is assumed, and dihydroartemisinin",
    "disposition is one-compartment. Artemether shows time-dependent PK: an",
    "enzyme turnover compartment (believed to represent CYP3A4) is driven by an",
    "Emax function of the artemether plasma concentration and multiplies the",
    "pre-induced artemether clearance, giving a 1.43-fold clearance increase by",
    "the sixth dosing occasion with an enzyme half-life of 30.4 h. Allometric",
    "body-weight scaling is applied to all clearances (exponent 0.75) and",
    "volumes (exponent 1.0), centered on the cohort-median 59 kg.",
    "Concentrations and doses are on a molar basis (nmol/L, nmol), matching the",
    "paper's log-molar estimation scale and the molar EC50. The companion",
    "lumefantrine model from the same paper is LohyDas_2018_lumefantrine; the",
    "authors modelled the two drugs separately because a simultaneous fit was",
    "unstable."
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
  # biological matrix. Units are derived from the `units` block above. The
  # enzyme state is dimensionless: it is normalised to 1 at baseline.
  compartmentData <- list(
    depot = list(analyte = "artemether", units = "nmol", specimen = "administration site", verified = FALSE),
    transit1 = list(analyte = "artemether", units = "nmol", specimen = "administration site", verified = FALSE),
    transit2 = list(analyte = "artemether", units = "nmol", specimen = "administration site", verified = FALSE),
    central = list(analyte = "artemether", units = "nmol", specimen = "plasma", verified = FALSE),
    central_dihydroart = list(analyte = "dihydroartemisinin", units = "nmol", specimen = "plasma", verified = FALSE),
    enzyme = list(
      analyte = "inducible metabolising enzyme (believed CYP3A4)",
      units = "(fraction of baseline)",
      specimen = "tissue",
      verified = FALSE
    )
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
        "improved the fit. Results: 'The final PK models for ARM-DHA and LF",
        "therefore incorporated only body weight implemented allometrically on",
        "clearances and volumes of distribution.' Cohort weight range",
        "40.0-65.0 kg (Table 1)."
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
        "single dosing occasion)'. Here OCC is purely the grouping for the",
        "between-occasion variability on relative bioavailability and on mean",
        "transit time (Results: 'The addition of BOV on both MTT (dOFV =",
        "-12.6) and F (dOFV = -9.39) yielded significant improvements in the",
        "fit'). The time dependence of artemether clearance is NOT an",
        "occasion-indexed covariate in this model: it is generated",
        "mechanistically by the enzyme turnover compartment, which is driven by",
        "the artemether concentration itself and therefore extrapolates to",
        "dosing regimens other than the studied six doses."
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
        "Screened as a covariate on the artemether and dihydroartemisinin PK",
        "parameters and found not significant: Lohy Das 2018 Results,",
        "'Covariates': 'As for ARM-DHA, none of the covariates tested",
        "(time-varying parasitemia density and EGA) had any effect on the PK",
        "parameters except for baseline parasitemia density on MTT in the",
        "forward step, but again this parameter was not maintained during the",
        "backward elimination step.' No point estimate is published. Cohort",
        "EGA range 15.7-39.0 weeks (Table 1)."
      )
    ),
    PARASITEMIA = list(
      description = "Asexual Plasmodium falciparum parasite density (baseline and time-varying)",
      units = "parasites/uL",
      type = "continuous",
      notes = paste(
        "Baseline parasitemia density on mean transit time was selected in the",
        "forward step but not retained in backward elimination (P <= 0.01):",
        "Lohy Das 2018 Results, 'Covariates'. No point estimate is published.",
        "Cohort baseline range 3,060-160,000 parasites/uL (Table 1)."
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
      "and/or a small cake (a fatty meal) to enhance absorption. Each 80 mg",
      "artemether dose is 268,124 nmol using the artemether molar mass",
      "298.37 g/mol."
    ),
    sampling = paste(
      "Venous plasma. A pre-dose sample, then samples at 2 and 4 h after each",
      "dose and immediately before doses 2-6 (troughs), with additional",
      "samples at 0.25 h after doses 1 and 2, at 6, 8 and 12 h after the last",
      "dose, and a scheduled day 7 sample. 387 samples were analyzed for",
      "artemether and dihydroartemisinin (9 excluded for hemolysis).",
      "Concentrations were measured by LC-MS/MS; LLOQ 1.43 ng/mL for both",
      "analytes, i.e. 4.79 nmol/L (artemether) and 5.03 nmol/L",
      "(dihydroartemisinin). 24% (artemether) and 28% (dihydroartemisinin) of",
      "observations were below the LLOQ; all were retained and handled with",
      "the likelihood-based M3 method under Laplacian estimation."
    ),
    regions = paste(
      "Rwanda (obstetrics and gynecology ward, Rwamagana district hospital,",
      "eastern province; mesoendemic transmission). Study conducted June 2007",
      "to July 2009; Rwanda National Ethics Committee study IRB 00001497. The",
      "PK study was nested within a pharmacovigilance study of ACT use in",
      "pregnancy."
    ),
    notes = paste(
      "NONMEM 7.3, Laplacian estimation with the M3 method for the censored",
      "artemether and dihydroartemisinin observations. Parameter imprecision",
      "was obtained by sampling importance resampling (SIR). Artemether and",
      "dihydroartemisinin were fitted simultaneously assuming complete and",
      "irreversible in-vivo conversion. The enzyme turnover model follows",
      "Hassan et al. and Smythe et al. (the paper's references 37 and 59).",
      "Two alternative autoinduction structures (inducible plus non-inducible",
      "parallel elimination pathways, in either order of which pathway forms",
      "dihydroartemisinin) were explored and did not improve the fit. Attempts",
      "to fit artemether-dihydroartemisinin and lumefantrine simultaneously,",
      "to capture the correlation in absorption between the two co-formulated",
      "drugs, were unstable and were abandoned, which is why this model and",
      "LohyDas_2018_lumefantrine are separate."
    )
  )

  ini({
    # Structural PK parameters, Lohy Das 2018 Table 2 (Artemether and
    # Dihydroartemisinin blocks). Concentrations were modelled on a natural-log
    # MOLAR scale (Methods, 'Population pharmacokinetics': 'The molar units of
    # LF, ARM, and DHA concentration were transformed to their natural
    # logarithms for this modeling analysis'; 'ARM and DHA, expressed as molar
    # concentrations, were characterized simultaneously assuming complete and
    # irreversible in vivo conversion of ARM into DHA'), so doses are supplied
    # in nmol and concentrations are nmol/L. Working in molar units is not
    # optional here: the autoinduction EC50 is reported in nM and is compared
    # directly against the artemether concentration. Typical values are
    # apparent (CL/F, V/F) at the cohort median WT = 59 kg.

    lfdepot <- fixed(log(1))
    label("Reference relative oral bioavailability of artemether, F (unitless) (the source paper fixes F at 1)")  # Lohy Das 2018 Table 2: F = '1 fixed'; Methods: 'The typical relative F was implemented as a fixed parameter for the parent analyte, i.e., LF and ARM (100% relative bioavailability).'

    lmtt <- log(0.738)
    label("Mean transit time of the 2-compartment transit-absorption chain, MTT (h)")  # Lohy Das 2018 Table 2: MTT = 0.738 h (%RSE 12.5, 90% CI 0.569-0.840); Discussion cross-check 'the mean absorption times differed substantially between the two drugs (45 min versus 4 h)'

    lcl <- log(467)
    label("Apparent PRE-INDUCED artemether elimination clearance, CL_ARM/F at WT = 59 kg (L/h); multiplied by the relative enzyme amount to give the time-varying clearance")  # Lohy Das 2018 Table 2: CL_ARM/F = 467 L/h (%RSE 17.9, 90% CI 298-508); Abstract 'The typical oral clearance, which started at 467 liters/h, increased 1.43-fold at the end of treatment.'

    lvc <- log(3000)
    label("Apparent artemether central volume of distribution, V_ARM/F at WT = 59 kg (L)")  # Lohy Das 2018 Table 2: V_ARM/F = 3,000 L (%RSE 14.1, 90% CI 2,050-3,180)

    lcl_dihydroart <- log(611)
    label("Apparent dihydroartemisinin elimination clearance, CL_DHA/F at WT = 59 kg (L/h)")  # Lohy Das 2018 Table 2: CL_DHA/F = 611 L/h (%RSE 15.4, 90% CI 486-782)

    lvc_dihydroart <- log(137)
    label("Apparent dihydroartemisinin central volume of distribution, V_DHA/F at WT = 59 kg (L)")  # Lohy Das 2018 Table 2: V_DHA/F = 137 L (%RSE 38.9, 90% CI 99.8-251)

    # Enzyme turnover / autoinduction parameters (Lohy Das 2018 Eqs. 1, 4, 5, 6
    # and Table 2 Artemether block). Eq. 1: EFF = Emax * CP / (EC50 + CP),
    # where CP is the artemether plasma concentration; Eq. 4:
    # dA_ENZ/dt = K_ENZ * (1 + EFF) - K_ENZ * A_ENZ; Eq. 5:
    # (CL_ARM/F)_induced = (CL_ARM/F)_preinduced * A_ENZ; Eq. 6:
    # K_ENZ = ln(2) / t_(1/2)ENZ.
    lemax <- log(0.986)
    label("Maximal autoinduction effect on the enzyme production rate, Emax (unitless)")  # Lohy Das 2018 Table 2: Emax = 0.986 (%RSE 22.8, 90% CI 0.623-1.42). Table 2 prints the unit as 'h-1', but Eq. 4 adds Emax-scaled EFF to 1 inside K_ENZ * (1 + EFF), so EFF and hence Emax must be dimensionless; the printed unit is treated as a typographical error (see the vignette Assumptions and deviations section).

    lec50 <- log(9.37)
    label("Artemether concentration producing half-maximal autoinduction, EC50 (nmol/L)")  # Lohy Das 2018 Table 2: EC50 = 9.37 nM (%RSE 25.4, 90% CI 6.16-14.4); Results 'the EC50 was 9.37 nM'

    lkenz <- log(log(2) / 30.4)
    label("First-order enzyme degradation rate constant, K_ENZ (1/h)")  # Lohy Das 2018 Table 2: TIME_ENZ (enzyme half-life) = 30.4 h (%RSE 42.1, 90% CI 7.59-41.9). The paper estimates the half-life and defines K_ENZ from it in Eq. 6: K_ENZ = ln(2)/t_(1/2)ENZ = 0.6931/30.4 = 0.02280 1/h. K_ENZ (not the half-life) is the quantity that appears in the enzyme ODE.

    # Allometric exponents. Lohy Das 2018 Methods: 'All size descriptors were
    # scaled to their respective medians (i.e., total body weight [59 kg] ...)
    # on PK parameters using allometric power exponents of 0.75 for clearances
    # (CL/F, Q/F, CL_ARM/F, and CL_DHA/F) and 1 for volumes of distribution
    # (Vc/F, Vp/F, V_ARM/F, and V_DHA/F).' The exponents are structural choices
    # imposed by the authors rather than estimated quantities (no RSE or CI is
    # given for either), so both are encoded as fixed.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on clearances (CL_ARM/F, CL_DHA/F), unitless")  # Lohy Das 2018 Methods: allometric power exponent 0.75 for clearances
    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on volumes (V_ARM/F, V_DHA/F), unitless")  # Lohy Das 2018 Methods: allometric power exponent 1 for volumes of distribution

    # Between-subject variability. Lohy Das 2018 Table 2 footnote a states the
    # scale explicitly: 'Coefficient of variation (CV) for BSV and BOV was
    # calculated as 100 x (variance)^1/2', i.e. the printed %CV is 100 times
    # the eta-scale SD, so variance = (CV/100)^2. Results: 'BSV was estimated
    # for all structural parameters except for the Emax of the autoinduction
    # effect of ARM, due to instability of the model' -- Table 2 additionally
    # prints no BSV for EC50 or TIME_ENZ, and the table governs here.
    etalfdepot ~ 0.3318   # Table 2 row 'F' (Artemether), BSV = 57.6% CV (%RSE 36.8, 90% CI 43.2-78.8); variance = 0.576^2
    etalmtt ~ 1.2100      # Table 2 row 'MTT (h)' (Artemether), BSV = 110% CV (%RSE 32.5, 90% CI 86.2-143.2); variance = 1.10^2
    etalcl ~ 0.0778       # Table 2 row 'CL_ARM/F', BSV = 27.9% CV (%RSE 44.1, 90% CI 21.5-43.5); variance = 0.279^2
    etalvc ~ 0.0420       # Table 2 row 'V_ARM/F', BSV = 20.5% CV (%RSE 43.1, 90% CI 15.4-31.0); variance = 0.205^2
    etalcl_dihydroart ~ 0.0428   # Table 2 row 'CL_DHA/F', BSV = 20.7% CV (%RSE 50.2, 90% CI 12.9-29.9); variance = 0.207^2
    etalvc_dihydroart ~ 0.1640   # Table 2 row 'V_DHA/F', BSV = 40.5% CV (%RSE 48.9, 90% CI 17.5-51.5); variance = 0.405^2

    # Between-occasion variability, one eta per dosing occasion (six doses).
    # A single BOV variance is shared across occasions in NONMEM, so only the
    # first slot is estimable here and the remaining five are fixed to it.
    etaiov_fdepot_1 ~ 0.2323
    # Table 2 row 'F' (Artemether) second line, BOV = 48.2% CV (%RSE 35.6, 90% CI 39.1-63.9); variance = 0.482^2
    etaiov_fdepot_2 ~ fixed(0.2323)
    etaiov_fdepot_3 ~ fixed(0.2323)
    etaiov_fdepot_4 ~ fixed(0.2323)
    etaiov_fdepot_5 ~ fixed(0.2323)
    etaiov_fdepot_6 ~ fixed(0.2323)

    etaiov_mtt_1 ~ 0.2830
    # Table 2 row 'MTT (h)' (Artemether) second line, BOV = 53.2% CV (%RSE 21.6, 90% CI 54.1-97.9); variance = 0.532^2
    etaiov_mtt_2 ~ fixed(0.2830)
    etaiov_mtt_3 ~ fixed(0.2830)
    etaiov_mtt_4 ~ fixed(0.2830)
    etaiov_mtt_5 ~ fixed(0.2830)
    etaiov_mtt_6 ~ fixed(0.2830)

    # Residual error. Lohy Das 2018 Methods, 'Population pharmacokinetics':
    # 'The unexplained residual error was estimated using an additive error
    # model on the logarithmic scale for all drugs, which equates to an
    # exponential error model on an arithmetic scale. In the case of ARM and
    # DHA, a separate additive error model was used for each analyte.' By the
    # standing nlmixr2lib convention an additive-on-log-scale residual maps to
    # a proportional residual in linear space, with propSd equal to the
    # log-scale SD. Table 2 reports each as a percentage. Both are large,
    # consistent with the erratic absorption and the 24-28% of observations
    # below the limit of quantification.
    propSd <- 0.984
    label("Proportional residual SD for artemether plasma concentration (SD on the log scale, approximately CV in linear space)")  # Lohy Das 2018 Table 2: RUV (Artemether) = 98.4% (%RSE 5.54, 90% CI 92.0-108)
    propSd_dihydroart <- 1.13
    label("Proportional residual SD for dihydroartemisinin plasma concentration (SD on the log scale, approximately CV in linear space)")  # Lohy Das 2018 Table 2: RUV (Dihydroartemisinin) = 113% (%RSE 6.01, 90% CI 109-129)
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
    # 'Covariates'). `cl` is the PRE-INDUCED artemether clearance; the
    # time-varying induced clearance is formed below from the enzyme state.
    cl <- exp(lcl + etalcl) * (WT / 59)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 59)^e_wt_vc
    cl_dihydroart <- exp(lcl_dihydroart + etalcl_dihydroart) * (WT / 59)^e_wt_cl
    vc_dihydroart <- exp(lvc_dihydroart + etalvc_dihydroart) * (WT / 59)^e_wt_vc

    emax <- exp(lemax)
    ec50 <- exp(lec50)
    kenz <- exp(lkenz)

    # Mean transit time carries both between-subject and between-occasion
    # random effects. Transit-chain rate constant in the Savic (2007)
    # parameterisation stated in the Lohy Das 2018 Fig. 1A legend, 'ktr,
    # transit absorption rate constant [ktr = (n + 1)/mean transit time]'. With
    # n = 2 transit compartments (Results: 'An absorption model consisting of 2
    # transit compartments was superior to other explored absorption models')
    # the chain depot -> transit1 -> transit2 -> central has n + 1 = 3
    # equal-rate transfers, so ktr = 3 / MTT.
    mtt <- exp(lmtt + etalmtt + iov_mtt)
    ktr <- 3 / mtt

    # Plasma concentrations in nmol/L (dose in nmol, volume in L). Cc must be
    # defined before it is used in the autoinduction term below.
    Cc <- central / vc
    Cc_dihydroart <- central_dihydroart / vc_dihydroart

    # Autoinduction. Eq. 1: the artemether plasma concentration stimulates the
    # enzyme production rate through an Emax relationship. Eq. 5: the relative
    # enzyme amount multiplies the pre-induced artemether clearance. EC50 is in
    # nmol/L, matching Cc.
    eff <- emax * Cc / (ec50 + Cc)
    cl_induced <- cl * enzyme

    # ODE system: oral depot -> 2 transit compartments -> artemether central ->
    # dihydroartemisinin central -> elimination (Lohy Das 2018 Fig. 1B), plus
    # the enzyme turnover compartment. Under the paper's assumption of complete
    # and irreversible in-vivo conversion, ALL artemether clearance is
    # metabolic conversion, so the flux leaving artemether central enters
    # dihydroartemisinin central unchanged. Both analytes are tracked on a
    # molar basis, so no molecular-weight factor is applied at the conversion
    # step (the paper modelled on the same molar basis).
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(central) <- ktr * transit2 - (cl_induced / vc) * central
    d/dt(central_dihydroart) <- (cl_induced / vc) * central -
      (cl_dihydroart / vc_dihydroart) * central_dihydroart

    # Eq. 4: enzyme turnover. Methods: 'The enzyme concentration was
    # initialized at 1 in order to normalize it to unity at baseline; i.e., the
    # zero-order production rate of the enzyme was set to K_ENZ.' The state is
    # therefore dimensionless and equals the enzyme amount relative to its
    # pre-induced baseline.
    d/dt(enzyme) <- kenz * (1 + eff) - kenz * enzyme
    enzyme(0) <- 1

    # Relative bioavailability on the depot compartment. lfdepot is fixed at
    # log(1) (F = 1 fixed in Table 2); the between-subject and between-occasion
    # random effects carry the variability.
    f(depot) <- exp(lfdepot + etalfdepot + iov_fdepot)

    # NONMEM additive-on-log-scale residual maps to nlmixr2 proportional, one
    # per analyte.
    Cc ~ prop(propSd)
    Cc_dihydroart ~ prop(propSd_dihydroart)
  })
}
