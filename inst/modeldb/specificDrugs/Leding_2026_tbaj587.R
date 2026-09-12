Leding_2026_tbaj587 <- function() {
  description <- "Joint parent + metabolite population PK model for the second-in-class diarylquinoline (DARQ) antitubercular TBAJ-587 and its two main metabolites M2 and M3 in healthy adult volunteers after single oral doses of 25-800 mg, followed to Day 126. TBAJ-587 has Savic 2007 analytical transit-compartment absorption (non-integer NN = 2.36 transit compartments feeding an absorption compartment that empties at ka) and three-compartment disposition; M3 has three-compartment and M2 two-compartment disposition, each formed in parallel from the parent central compartment. Both relative fractions metabolised are fixed at 1, so every metabolite clearance and volume is apparent relative to an unidentifiable true fraction metabolised, and each metabolite receives a formation flux equal to the parent's whole elimination flux. Dose (a continuous covariate referenced to 200 mg) acts as a power function on parent clearance, on both fractions metabolised and on M2 clearance, and as an exponential function on ka; a high-calorie high-fat meal raises relative bioavailability and mean transit time and lowers both fractions metabolised. Residual variability is additive on the natural-log scale (lognormal) and estimated separately for each of the three analytes. Because no intravenous data were collected, all parent clearances and volumes are apparent oral values."
  reference <- paste(
    "Leding A. A. M., Bruinenberg P., Conradie A., Nedelman J.,",
    "Lombardi A., Hickman D., Simonsson U. S. H. (2026).",
    "Population pharmacokinetics of TBAJ-587 and its main metabolites -",
    "Evaluation of different loading dose strategies and early dose",
    "selection.",
    "British Journal of Clinical Pharmacology 92(4):1058-1068.",
    "doi:10.1002/bcp.70333.",
    "Structure and all covariate / error equations transcribed from the",
    "final NONMEM control stream in Supporting Information Code S1;",
    "final parameter values from Table 1 (the Code S1 $THETA / $OMEGA",
    "records for the metabolite sub-models are initial estimates, not",
    "final ones).",
    sep = " "
  )
  vignette <- "Leding_2026_tbaj587"

  # Modelling was performed on natural-log-transformed MOLAR concentration
  # data: Code S1 $INPUT declares 'AMT ; dose in nmol' and
  # 'DV ; concentration in ln(nmol/L)'. The amount unit is therefore nmol
  # and the concentration unit nmol/L. The separate mg dose that drives the
  # dose covariates is carried by the DOSE_TBAJ587_MG covariate column
  # (Code S1 $INPUT 'DOSE ; dose in mg'), because the paper does not report
  # the molecular weights needed to convert between the two.
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(
      analyte = "TBAJ-587", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "TBAJ-587", units = "nmol", specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "TBAJ-587", units = "nmol", specimen = "plasma", verified = TRUE
    ),
    peripheral2 = list(
      analyte = "TBAJ-587", units = "nmol", specimen = "plasma", verified = TRUE
    ),
    central_m3 = list(
      analyte = "TBAJ-587 metabolite M3", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_m3 = list(
      analyte = "TBAJ-587 metabolite M3", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral2_m3 = list(
      analyte = "TBAJ-587 metabolite M3", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    central_m2 = list(
      analyte = "TBAJ-587 metabolite M2", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_m2 = list(
      analyte = "TBAJ-587 metabolite M2", units = "nmol",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    DOSE_TBAJ587_MG = list(
      description        = "Administered TBAJ-587 oral dose in mg, referenced to 200 mg",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Dose was treated as a CONTINUOUS covariate in stepwise covariate",
        "modelling (Leding 2026 Methods 2.2: 'Dose was treated as a",
        "continuous covariate and evaluated on absorption (ka and MTT),",
        "relative bioavailability (F), metabolite formation (fm,M2 and",
        "fm,M3) and elimination parameters'), referenced to 200 mg.",
        "Five effects were retained, in the exact forms given by Code S1:",
        "power on parent CL, CLDOSE = (DOSE/200)**THETA(12);",
        "exponential on ka, KADOSE = EXP(THETA(13)*(DOSE - 200));",
        "power on fm,M3, FM3DOSE = (DOSE/200)**THETA(22);",
        "power on fm,M2, FM2DOSE = (DOSE/200)**THETA(29);",
        "and power on CL_M2, CLM2DOSE = (DOSE/200)**THETA(31).",
        "This column is NOT the dosing amount: Code S1 doses AMT in nmol",
        "while this covariate is the mg dose that labels the record, so",
        "both are needed and the paper reports no molecular weight to",
        "convert between them. Studied levels were 25, 50, 100, 200, 400",
        "and 800 mg (Leding 2026 Methods 2.1). Set it on every record of a",
        "subject to that subject's mg dose level; for multiple-dose",
        "simulation it is the daily mg dose, and during a loading-dose",
        "period it takes the loading mg dose (Leding 2026 Methods 2.4).",
        "Member of the DOSE_<drug>_<units> canonical family."
      ),
      source_name        = "DOSE"
    ),
    FED_HIGHFAT = list(
      description        = "High-calorie high-fat meal taken with the dose (1 = fed, 0 = fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted; the six single-ascending-dose cohorts)",
      notes              = paste(
        "Leding 2026 Methods 2.1: 'Nine subjects in the food-effect cohort",
        "received a high-calorie and high-fat meal together with a single",
        "200 mg oral dose of TBAJ-587', so FED_HIGHFAT rather than the",
        "generic FED carries the high-fat semantic. Code S1 names the",
        "column FED ('food state, fed = 1, fasted = 0') with the same",
        "orientation as the canonical. Three retained effects, each the",
        "linear-deviation form P_fed = P_fasted * (1 + Covariate_Fed)",
        "given in Table 1 footnote d: relative bioavailability",
        "F * (1 + 0.688), mean transit time MTT * (1 + 0.958), and both",
        "fractions metabolised, fm,M3 * (1 - 0.479) and",
        "fm,M2 * (1 - 0.547). The effect is BETWEEN-subject here, not a",
        "crossover: the food-effect cohort is a separate group of nine",
        "subjects dosed only at 200 mg, which is why Leding 2026",
        "Discussion cautions that fed simulations assume the",
        "fasting-derived dose dependence of apparent oral clearance also",
        "holds under fed conditions."
      ),
      source_name        = "FED"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    age_range      = "18-64 years",
    bmi_range      = "15.5-32.0 kg/m^2",
    weight_range   = "minimum 50.0 kg (no upper bound or median reported)",
    disease_state  = "healthy volunteers",
    dose_range     = "single oral doses of 25, 50, 100, 200, 400 or 800 mg (n = 6, 4, 5, 8, 4 and 6 respectively), plus a nine-subject food-effect cohort dosed 200 mg with a high-calorie high-fat meal",
    trial          = "NCT04890535 Part 1; partially blinded, placebo-controlled, randomized single ascending dose with food-effect cohort",
    regions        = "Netherlands (single site; ethics approval ID NL73973.056.20)",
    notes          = paste(
      "Leding 2026 Methods 2.1 and Results 3.1. Participants were male or",
      "female of non-childbearing potential; the sex split is not",
      "reported, so sex_female_pct is deliberately absent rather than",
      "guessed. 1929 observations were available for each of the three",
      "analytes, collected to Day 126; 29 of 42 subjects (69%) were",
      "sampled beyond Day 21 and 23 of those 29 (79%) beyond Day 112.",
      "Below-quantification-limit observations (LLOQ 1 ng/mL) made up 2%,",
      "46% and 28% of records for TBAJ-587, M2 and M3 respectively; all",
      "parent BQL records were omitted, and metabolite BQL records were",
      "omitted except the first and last in the decreasing and ascending",
      "curve, which were set to LLOQ/2 (Leding 2026 Results 3.1)."
    )
  )

  ini({
    # ---- TBAJ-587 (parent) sub-model -------------------------------------
    # Table 1 footnote a: because no intravenous data were available, all
    # clearances and volumes are apparent oral values. Only the three
    # by-fiat anchors (F, fm,M3, fm,M2) are wrapped in fixed(); every other
    # value carries an RSE and a SIR confidence interval in Table 1 and so
    # was estimated. Code S1 marks the parent block FIX because the final
    # model estimated the metabolite sub-models with the parent sub-model
    # held at its own final estimates (Leding 2026 Methods 2.2), which is
    # an estimation-sequence device rather than a structurally fixed value.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability in the fasted state (unitless)")  # Table 1 'F fasted' = 1 FIX; Code S1 $THETA(2) '(1) FIX ; 2 BIO'

    lka <- log(0.0866)
    label("Absorption rate constant from the absorption compartment at a 200 mg dose (1/h)")  # Table 1 'ka (h-1)' = 0.0866, RSE 3.5%, CI 0.0818-0.0930; Code S1 $THETA(1)

    lcl <- log(5.40)
    label("Apparent oral clearance of TBAJ-587, fasted, at a 200 mg dose (L/h)")  # Table 1 'CL TBAJ-587 (L/h)' = 5.40, RSE 8.4%, CI 4.70-6.24; Code S1 $THETA(3)

    lvc <- log(88.2)
    label("Apparent central volume of distribution of TBAJ-587 (L)")  # Table 1 'V C,TBAJ-587 (L)' = 88.2, RSE 6.8%, CI 76.7-100; Code S1 $THETA(4)

    lq <- log(31.1)
    label("Apparent inter-compartmental clearance to the first TBAJ-587 peripheral compartment (L/h)")  # Table 1 'Q P1,TBAJ-587 (L/h)' = 31.1, RSE 6.4%, CI 28.9-33.7; Code S1 $THETA(5)

    lvp <- log(1910)
    label("Apparent first peripheral volume of distribution of TBAJ-587 (L)")  # Table 1 'V P1,TBAJ-587 (L)' = 1910, RSE 8.7%, CI 1640-2180; Code S1 $THETA(6)

    lq2 <- log(31.5)
    label("Apparent inter-compartmental clearance to the second TBAJ-587 peripheral compartment (L/h)")  # Table 1 'Q P2,TBAJ-587 (L/h)' = 31.5, RSE 4.0%, CI 29.2-33.5; Code S1 $THETA(7)

    lvp2 <- log(22500)
    label("Apparent second peripheral volume of distribution of TBAJ-587 (L)")  # Table 1 'V P2,TBAJ-587 (L)' = 22 500, RSE 3.9%, CI 20 700-24 900; Code S1 $THETA(8)

    lmtt <- log(0.837)
    label("Mean absorption transit time in the fasted state (h)")  # Table 1 'MTT fasted (h)' = 0.837, RSE 4.5%, CI 0.763-0.912; Code S1 $THETA(9)

    lnn <- log(2.36)
    label("Number of absorption transit compartments (unitless, non-integer)")  # Table 1 'NN' = 2.36, RSE 8%, CI 2.14-2.64; Code S1 $THETA(10)

    e_fed_highfat_fdepot <- 0.688
    label("Fractional change in relative bioavailability with a high-fat meal (unitless)")  # Table 1 'Fed state on F fasted' = 0.688, RSE 18.6%, CI 0.392-0.965; applied as F*(1 + 0.688) per Table 1 footnote d; Code S1 $THETA(11)

    e_dose_cl <- 0.298
    label("Power exponent on (dose/200 mg) for apparent oral clearance of TBAJ-587 (unitless)")  # Table 1 'Dose covariate on CL TBAJ-587' = 0.298, RSE 28.8%, CI 0.144-0.466; power form per Table 1 footnote e; Code S1 $THETA(12)

    e_dose_ka <- -0.000462
    label("Exponential coefficient on (dose - 200 mg) for the absorption rate constant (1/mg)")  # Code S1 $THETA(13) '(-0.008, -0.000462, 0.008) FIX ; KADOSE1' -- Table 1 'Dose covariate on ka' prints the same value rounded to -0.00046 (RSE 32.0%, CI -0.000669 to -0.000277); exponential form per Table 1 footnote f

    e_fed_highfat_mtt <- 0.958
    label("Fractional change in mean absorption transit time with a high-fat meal (unitless)")  # Table 1 'Fed state on MTT fasted' = 0.958, RSE 14.0%, CI 0.643-1.29; applied as MTT*(1 + 0.958) per Table 1 footnote d; Code S1 $THETA(14)

    # ---- M3 sub-model ----------------------------------------------------
    # Table 1 footnote a: for the metabolites, clearances and volumes are
    # additionally relative to the unknown fraction metabolised.
    lfm_m3 <- fixed(log(1))
    label("Relative fraction of TBAJ-587 metabolised to M3 at a 200 mg fasted dose (unitless)")  # Table 1 'f m,M3' = 1 FIX; Code S1 $THETA(15) '1 FIX ; 15 FM3'

    lcl_m3 <- log(18.5)
    label("Apparent clearance of M3 (L/h)")  # Table 1 'CL M3 (L/h)' = 18.5, RSE 7.8%, CI 16.7-20.6 (Code S1 $THETA(16) 20.121 is the INITIAL estimate)

    lvc_m3 <- log(816)
    label("Apparent central volume of distribution of M3 (L)")  # Table 1 'V C,M3 (L)' = 816, RSE 6.0%, CI 739-906 (Code S1 $THETA(17) 820.718 is the INITIAL estimate)

    lvp_m3 <- log(1630)
    label("Apparent first peripheral volume of distribution of M3 (L)")  # Table 1 'V P1,M3 (L)' = 1630, RSE 22.1%, CI 1170-2330 (Code S1 $THETA(18) V6 = 937.257 is the INITIAL estimate)

    lq_m3 <- log(9.09)
    label("Apparent inter-compartmental clearance to the first M3 peripheral compartment (L/h)")  # Table 1 'Q P1,M3 (L/h)' = 9.09, RSE 17.9%, CI 5.97-11.8 (Code S1 $THETA(19) Q3 = 17.4191 is the INITIAL estimate)

    lvp2_m3 <- log(3010)
    label("Apparent second peripheral volume of distribution of M3 (L)")  # Table 1 'V P2,M3 (L)' = 3010, RSE 6.5%, CI 2740-3340 (Code S1 $THETA(20) V7 = 2914.39 is the INITIAL estimate)

    lq2_m3 <- log(99.6)
    label("Apparent inter-compartmental clearance to the second M3 peripheral compartment (L/h)")  # Table 1 'Q P2,M3 (L)' = 99.6, RSE 6.4%, CI 90.9-108 (Code S1 $THETA(21) Q4 = 93.6415 is the INITIAL estimate)

    e_dose_fm_m3 <- -0.418
    label("Power exponent on (dose/200 mg) for the relative fraction metabolised to M3 (unitless)")  # Table 1 'Dose covariate on f m,M3' = -0.418, RSE 14.3%, CI -0.489 to -0.346; power form per Table 1 footnote e (Code S1 $THETA(22) -0.408881 is the INITIAL estimate)

    e_fed_highfat_fm_m3 <- -0.479
    label("Fractional change in the relative fraction metabolised to M3 with a high-fat meal (unitless)")  # Table 1 'Fed state on f m,M3' = -0.479, RSE 24.8%, CI -0.595 to -0.341, i.e. 47.1% lower fed (Results 3.2); applied as fm,M3*(1 - 0.479) per Table 1 footnote d (Code S1 $THETA(23) -0.472581 is the INITIAL estimate)

    # ---- M2 sub-model ----------------------------------------------------
    lfm_m2 <- fixed(log(1))
    label("Relative fraction of TBAJ-587 metabolised to M2 at a 200 mg fasted dose (unitless)")  # Table 1 'f m,M2' = 1 FIX; Code S1 $THETA(24) '1 FIX ; 24 FM2'

    lcl_m2 <- log(34.1)
    label("Apparent clearance of M2 at a 200 mg dose (L/h)")  # Table 1 'CL M2 (L/h)' = 34.1, RSE 6.1%, CI 30.9-37.5 (Code S1 $THETA(25) 34.8785 is the INITIAL estimate)

    lvc_m2 <- log(247)
    label("Apparent central volume of distribution of M2 (L)")  # Table 1 'V C,M2 (L)' = 247, RSE 8.5%, CI 216-288 (Code S1 $THETA(26) 246.255 is the INITIAL estimate)

    lvp_m2 <- log(12200)
    label("Apparent first peripheral volume of distribution of M2 (L)")  # Table 1 'V P1,M2 (L)' = 12 200, RSE 6.5%, CI 10 800-13 700 (Code S1 $THETA(27) V9 = 12650.2 is the INITIAL estimate)

    lq_m2 <- log(174)
    label("Apparent inter-compartmental clearance to the first M2 peripheral compartment (L/h)")  # Table 1 'Q P1,M2 (L/h)' = 174, RSE 5.0%, CI 131-189 (Code S1 $THETA(28) Q5 = 181.293 is the INITIAL estimate)

    e_dose_fm_m2 <- -0.373
    label("Power exponent on (dose/200 mg) for the relative fraction metabolised to M2 (unitless)")  # Table 1 'Dose covariate on f m,M2' = -0.373, RSE 16.5%, CI -0.451 to -0.294; power form per Table 1 footnote e (Code S1 $THETA(29) -0.377743 is the INITIAL estimate)

    e_fed_highfat_fm_m2 <- -0.547
    label("Fractional change in the relative fraction metabolised to M2 with a high-fat meal (unitless)")  # Table 1 'Fed state on f m,M2' = -0.547, RSE 20.8%, CI -0.657 to -0.0424, i.e. 54.7% lower fed (Results 3.2); applied as fm,M2*(1 - 0.547) per Table 1 footnote d (Code S1 $THETA(30) -0.543959 is the INITIAL estimate)

    e_dose_cl_m2 <- 0.146
    label("Power exponent on (dose/200 mg) for apparent M2 clearance (unitless)")  # Table 1 'Dose covariate on CL M2' = 0.146, RSE 31.2%, CI 0.0666-0.217; power form per Table 1 footnote e (Code S1 $THETA(31) 0.144014 is the INITIAL estimate)

    # ---- Interindividual variability -------------------------------------
    # Exponential IIV throughout, P_i = P * exp(eta_i) (Leding 2026
    # Equation 1), so each value below is an omega^2 on the log scale.
    # Parent omegas are taken from the Code S1 $OMEGA block, which carries
    # them to more digits than Table 1 and reproduces every Table 1 CV%
    # exactly under the paper's own footnote-c formula
    # CV = sqrt(exp(omega^2) - 1)*100%. Metabolite omegas are NOT in Code S1
    # at their final values (that block holds initial estimates), so they
    # are back-transformed from the Table 1 CV% column by inverting the same
    # formula: omega^2 = log(1 + CV^2).
    etalka     ~ 0.0354   # Code S1 $OMEGA(1) 'KA' = 0.0354 FIX; reproduces Table 1 ka CV 19.0%
    etalfdepot ~ 0.162    # Code S1 $OMEGA(2) 'BIO' = 0.162 FIX; reproduces Table 1 'F fasted' CV 41.9%
    etalcl     ~ 0.164    # Code S1 $OMEGA(3) 'CL' = 0.164 FIX; reproduces Table 1 CL TBAJ-587 CV 42.2%
    etalvc     ~ 0.209    # Code S1 $OMEGA(4) 'V2' = 0.209 FIX; reproduces Table 1 V C,TBAJ-587 CV 48.2%
    etalvp     ~ 0.138    # Code S1 $OMEGA(5) 'V3' = 0.138 FIX; reproduces Table 1 V P1,TBAJ-587 CV 38.5%
    etalvp2    ~ 0.111    # Code S1 $OMEGA(6) 'V4' = 0.111 FIX; reproduces Table 1 V P2,TBAJ-587 CV 34.3%
    etalmtt    ~ 0.129    # Code S1 $OMEGA(7) 'MTT' = 0.129 FIX; reproduces Table 1 MTT fasted CV 37.1%
    etalcl_m3  ~ 0.048532 # log(1 + 0.223^2) from Table 1 CL M3 IIV CV 22.3% (RSE 30.0%, CI 12.9-32.8)
    etalvc_m3  ~ 0.094605 # log(1 + 0.315^2) from Table 1 V C,M3 IIV CV 31.5% (RSE 18.1%, CI 23.5-43.1)
    etalfm_m2  ~ 0.015504 # log(1 + 0.125^2) from Table 1 f m,M2 IIV CV 12.5% (RSE 15.2%, CI 9.15-16.8)
    etalvc_m2  ~ 0.311905 # log(1 + 0.605^2) from Table 1 V C,M2 IIV CV 60.5% (RSE 12.9%, CI 46.5-78.4)
    etalvp_m2  ~ 0.219156 # log(1 + 0.495^2) from Table 1 V P1,M2 IIV CV 49.5% (RSE 11.2%, CI 38.7-61.4)

    # ---- Residual variability --------------------------------------------
    # Code S1 $ERROR sets IPRED = LOG(C) and Y = IPRED + EPS, i.e. an
    # additive error on the natural-log scale, which is exactly nlmixr2's
    # lnorm() residual. Table 1 reports the $SIGMA VARIANCE (its parent
    # value 0.0337 is identical to the Code S1 $SIGMA record), so each
    # expSd below is sqrt(variance).
    expSd <- 0.183576
    label("TBAJ-587 residual SD, additive on the natural-log scale (log units)")  # sqrt(0.0337); Table 1 'Additive residual error on logarithmic scale' = 0.0337, RSE 5.0%, CI 0.0315-0.0363; Code S1 $SIGMA 0.0337 FIX

    expSd_m3 <- 0.195448
    label("M3 residual SD, additive on the natural-log scale (log units)")  # sqrt(0.0382); Table 1 'M3 additive residual error on logarithmic scale' = 0.0382, RSE 6.6%, CI 0.0356-0.0412 (Code S1 $SIGMA 0.0372996 is the INITIAL estimate)

    expSd_m2 <- 0.163401
    label("M2 residual SD, additive on the natural-log scale (log units)")  # sqrt(0.0267); Table 1 'M2 additive residual error on logarithmic scale' = 0.0267, RSE 8.3%, CI 0.0246-0.0293 (Code S1 $SIGMA 0.0265703 is the INITIAL estimate)
  })

  model({
    # 1. Dose-level covariate factors, referenced to 200 mg. Code S1 $PK
    #    builds each of these as a named block; the power form is
    #    (DOSE/200)**THETA and the ka form is EXP(THETA*(DOSE - 200)).
    #    DOSE_TBAJ587_MG is the mg dose label of the record, NOT the nmol
    #    dosing amount.
    dosefac_cl    <- (DOSE_TBAJ587_MG / 200)^e_dose_cl
    dosefac_ka    <- exp(e_dose_ka * (DOSE_TBAJ587_MG - 200))
    dosefac_cl_m2 <- (DOSE_TBAJ587_MG / 200)^e_dose_cl_m2
    dosefac_fm_m3 <- (DOSE_TBAJ587_MG / 200)^e_dose_fm_m3
    dosefac_fm_m2 <- (DOSE_TBAJ587_MG / 200)^e_dose_fm_m2

    # 2. High-fat-meal factors, each the linear-deviation form
    #    P_fed = P_fasted * (1 + Covariate_Fed) of Table 1 footnote d
    #    (Code S1 'IF(FED.EQ.1) XFED = (1 + THETA(i))').
    fedfac_fdepot <- 1 + e_fed_highfat_fdepot * FED_HIGHFAT
    fedfac_mtt    <- 1 + e_fed_highfat_mtt    * FED_HIGHFAT
    fedfac_fm_m3  <- 1 + e_fed_highfat_fm_m3  * FED_HIGHFAT
    fedfac_fm_m2  <- 1 + e_fed_highfat_fm_m2  * FED_HIGHFAT

    # 3. Individual parameters. Code S1 applies the covariate factor to the
    #    typical value and the eta on top of that (TVKA = KACOV*TVKA then
    #    KA = TVKA*EXP(ETA(1))), which is the multiplicative order used
    #    here. Parameters with no Table 1 IIV row carry no eta: Q and Q2 for
    #    the parent, every M3 parameter except CL_M3 and V_C,M3, and CL_M2,
    #    Q_P1,M2 and f m,M3.
    ka     <- exp(lka + etalka) * dosefac_ka
    fdepot <- exp(lfdepot + etalfdepot) * fedfac_fdepot
    cl     <- exp(lcl + etalcl) * dosefac_cl
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp + etalvp)
    q2     <- exp(lq2)
    vp2    <- exp(lvp2 + etalvp2)
    mtt    <- exp(lmtt + etalmtt) * fedfac_mtt
    nn     <- exp(lnn)

    fm_m3  <- exp(lfm_m3) * dosefac_fm_m3 * fedfac_fm_m3
    cl_m3  <- exp(lcl_m3 + etalcl_m3)
    vc_m3  <- exp(lvc_m3 + etalvc_m3)
    vp_m3  <- exp(lvp_m3)
    q_m3   <- exp(lq_m3)
    vp2_m3 <- exp(lvp2_m3)
    q2_m3  <- exp(lq2_m3)

    fm_m2  <- exp(lfm_m2 + etalfm_m2) * dosefac_fm_m2 * fedfac_fm_m2
    cl_m2  <- exp(lcl_m2) * dosefac_cl_m2
    vc_m2  <- exp(lvc_m2 + etalvc_m2)
    vp_m2  <- exp(lvp_m2 + etalvp_m2)
    q_m2   <- exp(lq_m2)

    # 4. Savic 2007 analytical transit-compartment absorption input. Code S1
    #    computes KTR = (NN+1)/MTT -- note the NN+1, which only the control
    #    stream settles -- and normalises the gamma density with the
    #    STIRLING approximation to log(NN!),
    #    LNFAC = LOG(2.5066) + (NN+0.5)*LOG(NN) - NN,
    #    rather than the exact lgamma(NN+1). At NN = 2.36 Stirling
    #    understates log(NN!) by 0.0351, so this input rate integrates to
    #    1.0357 * fdepot * dose rather than exactly fdepot * dose. The
    #    approximation is reproduced verbatim because it is the expression
    #    against which every Table 1 value (CL especially) was estimated;
    #    see the vignette Assumptions and deviations section, which
    #    measures the factor. The two 0.00001 offsets are Code S1's own
    #    guards against log(0) at the dosing instant.
    ktr   <- (nn + 1) / mtt
    lnfac <- log(2.5066) + (nn + 0.5) * log(nn) - nn
    trin  <- exp(log(fdepot * podo(depot) + 0.00001) + log(ktr) +
                   nn * log(ktr * tad(depot) + 0.00001) -
                   ktr * tad(depot) - lnfac)

    # 5. ODE system, one line per Code S1 $DES equation.
    #    A(1) = depot (Code S1 COMP(ABSP)), A(2)-A(4) = TBAJ-587 central
    #    and two peripherals, A(5)-A(7) = M3 central and two peripherals,
    #    A(8)-A(9) = M2 central and one peripheral.
    #
    #    PARALLEL FORMATION WITH BOTH FRACTIONS ANCHORED AT 1. Code S1 sets
    #    K25 = CL/V2*FM3 and K28 = CL/V2*FM2 as inflows to the two
    #    metabolite central compartments while DADT(2) retains its full
    #    -K20*A(2) = -CL/V2*A(2) elimination term. The parent's total
    #    elimination is therefore exactly CL/V2 and each metabolite
    #    receives a formation flux equal to fm times that same whole flux,
    #    so the fm values anchor at 1 instead of summing to 1 and the
    #    metabolite amounts are not mass-conserved against the parent.
    #    That is the paper's stated design: Methods 2.2 develops the
    #    metabolites "assuming parallel formation ... and the assumption of
    #    f m, the fraction of TBAJ-587 cleared to metabolite, equalling
    #    one", and Table 1 footnote a records that metabolite clearances
    #    and volumes are "additionally relative to unknown fractions
    #    metabolized". Do not add a compensating loss term to the parent:
    #    it would change every apparent metabolite parameter in Table 1.
    d/dt(depot)          <- trin - ka * depot

    d/dt(central)        <- ka * depot -
      cl * central / vc -
      q  * central / vc + q  * peripheral1 / vp -
      q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1)    <- q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2)    <- q2 * central / vc - q2 * peripheral2 / vp2

    d/dt(central_m3)     <- fm_m3 * cl * central / vc -
      cl_m3 * central_m3 / vc_m3 -
      q_m3  * central_m3 / vc_m3 + q_m3  * peripheral1_m3 / vp_m3 -
      q2_m3 * central_m3 / vc_m3 + q2_m3 * peripheral2_m3 / vp2_m3
    d/dt(peripheral1_m3) <- q_m3  * central_m3 / vc_m3 - q_m3  * peripheral1_m3 / vp_m3
    d/dt(peripheral2_m3) <- q2_m3 * central_m3 / vc_m3 - q2_m3 * peripheral2_m3 / vp2_m3

    d/dt(central_m2)     <- fm_m2 * cl * central / vc -
      cl_m2 * central_m2 / vc_m2 -
      q_m2  * central_m2 / vc_m2 + q_m2 * peripheral1_m2 / vp_m2
    d/dt(peripheral1_m2) <- q_m2  * central_m2 / vc_m2 - q_m2 * peripheral1_m2 / vp_m2

    # 6. Code S1 'F1 = 0 ; The amount is explicitly used in differential
    #    equation describing the absorption process'. Relative
    #    bioavailability enters through trin instead, via podo(depot).
    f(depot) <- 0

    # 7. Observed plasma concentrations: amount (nmol) / volume (L).
    Cc    <- central    / vc
    Cc_m3 <- central_m3 / vc_m3
    Cc_m2 <- central_m2 / vc_m2

    Cc    ~ lnorm(expSd)
    Cc_m3 ~ lnorm(expSd_m3)
    Cc_m2 ~ lnorm(expSd_m2)
  })
}
