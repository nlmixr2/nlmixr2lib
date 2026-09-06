Sager_2023_sotrovimab <- function() {
  description <- paste0(
    "Two-compartment population PK model with linear elimination for the ",
    "anti-SARS-CoV-2 spike-protein monoclonal antibody sotrovimab (VIR-7831 ",
    "/ GSK4182136) after a single intravenous or intramuscular dose in ",
    "non-hospitalized adults and adolescents with mild-to-moderate COVID-19 ",
    "at high risk of progression, plus healthy volunteers (Sager 2023, ",
    "n = 1984 across five studies: COMET-ICE, COMET-TAIL, COMET-PEAK, ",
    "BLAZE-4 and a Japanese/Caucasian healthy-volunteer study). ",
    "Intramuscular absorption is the authors' 'sigmoid absorption model': a ",
    "zero-order input of rate R1 = 130 mg/h into a depot compartment ",
    "followed by first-order absorption KA into the central compartment, ",
    "which together produce an S-shaped cumulative-absorption profile. ",
    "Intravenous doses enter the central compartment directly. Covariate ",
    "effects retained after forward selection / backward elimination are ",
    "body weight on CL and on the peripheral volume V3, sex on the ",
    "intramuscular bioavailability F_IM and on KA, and body mass index on ",
    "KA; the authors report that none of these reach clinical relevance ",
    "(all geometric-mean ratios within the prespecified 0.5-2.0 bounds). ",
    "Interindividual variability is a full 5x5 covariance block on CL, V2, ",
    "V3, KA and the logit of F_IM, and residual variability is combined ",
    "additive plus constant-coefficient-of-variation. The companion ",
    "exposure-response model for progression of COVID-19 through day 29 is ",
    "Sager_2023_sotrovimab_progression, which consumes the 168 h ",
    "concentration this model predicts."
  )
  reference <- paste(
    "Sager JE, El-Zailik A, Passarell J, Roepcke S, Li X, Aldinger M,",
    "Nader A, Skingsley A, Alexander EL, Yeh WW, Mogalian E, Garner C,",
    "Peppercorn A, Shapiro AE, Reyes M.",
    "Population pharmacokinetics and exposure-response analysis of a single",
    "dose of sotrovimab in the early treatment of patients with mild to",
    "moderate COVID-19.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(6):853-864.",
    "doi:10.1002/psp4.12958.",
    "Final parameter estimates are Table 1; the analysis-population",
    "demographics and the derived-parameter summary are Tables S5 and S6",
    "of the Supporting Information, and the NONMEM control stream is",
    "MODEL CODE S1 of the same document.",
    "Companion exposure-response model:",
    "modellib('Sager_2023_sotrovimab_progression').",
    sep = " "
  )
  vignette <- "Sager_2023_sotrovimab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters as a power model on CL and on the peripheral volume V3, ",
        "each normalised to 83.6 kg -- the median body weight of the PopPK ",
        "analysis population (Sager 2023 Results, 'PopPK demographics and ",
        "baseline characteristics', and the Table 1 row labels ",
        "'... in participants of 83.6 kg'). Note that the exponents are ",
        "estimated (0.494 on CL, 0.757 on V3), not fixed at the allometric ",
        "0.75 / 1. Body weight was not retained on V2 or Q."
      ),
      source_name        = "body weight"
    ),
    BMI = list(
      description        = "Body mass index at baseline.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters as a power model on KA normalised to 30.41 kg/m^2 (Sager ",
        "2023 Table 1 row label 'Absorption rate in male participants with ",
        "BMI of 30.41 kg/m2'); the population median BMI is reported as ",
        "30.4 kg/m^2 in the Results text. The exponent is negative ",
        "(-0.711): higher BMI slows intramuscular absorption. Relevant only ",
        "to intramuscular dosing."
      ),
      source_name        = "body mass index (BMI)"
    ),
    SEXF = list(
      description        = "Sex indicator; 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; every Table 1 typical value is quoted for male participants)",
      notes              = paste0(
        "Two distinct effects, both intramuscular-only. (1) On KA as a ",
        "proportional shift: KA_female = KA_male * (1 - 0.323). (2) On F_IM ",
        "as an ADDITIVE shift of -0.449 on the LOGIT scale, per the Sager ",
        "2023 Table 1 footnote b transformation. The paper's own numbers ",
        "pin this down: logit(0.582) - 0.449 = -0.1180, and ",
        "expit(-0.1180) = 0.4705, matching the 0.471 female bioavailability ",
        "quoted in the Results text. The analysis population was 44.9% ",
        "male, i.e. 55.1% female."
      ),
      source_name        = "sex"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Sager 2023 Methods, 'PopPK model development') but not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin concentration at baseline.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened but not retained. Albumin was unavailable in the healthy-volunteer study, so the screen for this factor was restricted to the four COVID-19 studies (Sager 2023 Methods)."
    ),
    DIS_COVID19 = list(
      description = "Disease-state indicator; 1 = patient with COVID-19, 0 = healthy volunteer.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'disease state (healthy vs. COVID-19)' but not retained. The analysis set contained 1891 patients with COVID-19 and 38 healthy volunteers."
    ),
    SARS_VLOAD = list(
      description = "Baseline SARS-CoV-2 viral load.",
      units       = "log10 copies/mL",
      type        = "continuous",
      notes       = "Screened but not retained in the PopPK model. It is also screened, and likewise not retained, in the companion exposure-response model."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "sotrovimab", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "sotrovimab", units = "mg",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "sotrovimab", units = "mg",
      specimen = "serum", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1984L,
    n_studies      = 5L,
    n_observations = "11,772 serum sotrovimab concentrations retained of 14,269 contributed (2497 excluded for missing/duplicate sample times, below-quantification predose and postdose samples, measurable predose concentrations, CWRES outliers beyond +/-5, or non-physiologic end-of-infusion concentrations); Sager 2023 Results and Table S3",
    age_median     = "49 years (range 15-96)",
    weight_median  = "83.6 kg (range 44.0-183.0)",
    bmi_median     = "30.4 kg/m^2 (range 15.9-71.1)",
    sex_female_pct = 55.1,
    race_ethnicity = c(White = 88.4, `Black/African American` = 6.2,
                       Asian = 4.1, Other = 0.7, Missing = 0.6),
    disease_state  = "non-hospitalized mild-to-moderate COVID-19 at high risk of progression to hospitalization or death (1891 participants), plus 38 healthy volunteers of Japanese or Caucasian descent",
    renal_function = "1415 normal, 401 mild, 65 moderate, 5 severe impairment (MDRD-based eGFR categories)",
    hepatic_function = "1393 normal, 487 mild, 5 moderate impairment (NCI criteria)",
    dose_range     = "single dose: 500 mg intravenous, or 250 mg or 500 mg intramuscular",
    regions        = "COMET-ICE, COMET-TAIL, COMET-PEAK and BLAZE-4 were multinational; COMET-TAIL enrolment was 85% from Florida, USA",
    notes          = paste0(
      "Five studies: COMET-ICE (NCT04545060), COMET-TAIL (NCT04913675), ",
      "COMET-PEAK (NCT04779879), BLAZE-4 (NCT04634409) and a healthy-",
      "volunteer study in participants of Japanese or Caucasian descent ",
      "(NCT04988152). In BLAZE-4 sotrovimab 500 mg i.v. was co-administered ",
      "with bamlanivimab. Demographics are summarised in Sager 2023 Results ",
      "('PopPK demographics and baseline characteristics') and Table S5; ",
      "the per-study participant and sample counts are in Table S4. Three ",
      "population sizes appear in the source and are not interchangeable: ",
      "1984 participants CONTRIBUTED concentrations, Table S5 tabulates ",
      "demographics for the 1929 in the analysis set (COMET-ICE 503, ",
      "COMET-TAIL 945, COMET-PEAK 348, BLAZE-4 95, GSK 217653 38), and ",
      "Table S6 summarises derived PK parameters for 1927 of them (739 of ",
      "whom received an intramuscular dose and so have a bioavailability ",
      "estimate). Route split across the analysis set: 1190 (61.7%) 500 mg ",
      "i.v., 478 (24.8%) 500 mg i.m., 261 (13.5%) 250 mg i.m. Table S6 ",
      "also gives the derived-parameter distributions the vignette gates ",
      "against: geometric-mean CL 0.095 L/day (35.0 %CV), V2 3.287 L, V3 ",
      "4.429 L, Vss 7.881 L, terminal half-life 60.67 days (median 61.2), ",
      "and intramuscular bioavailability median 0.497."
    )
  )

  ini({
    # ==================================================================
    # All values are the final population-mean estimates in Sager 2023
    # Table 1 ("Final population PK parameter estimates and covariate
    # effects"), corroborated by the Results narrative under "Final
    # model".
    #
    # TIME UNIT. Table 1 is internally mixed: CL and Q are reported in
    # L/day while KA and R1 are reported per hour. This file uses hours
    # throughout (units$time = "h"), because the absorption process and
    # the exposure metric that drives the companion exposure-response
    # model (concentration at 168 h) are both expressed in hours. CL and
    # Q are therefore divided by 24 at the trace site so the published
    # L/day number stays visible.
    #
    # CONCENTRATION UNIT. Doses are in mg and volumes in L, so
    # central/vc is mg/L = ug/mL, matching the assay units (lower limit
    # of quantification 0.1 ug/mL; Sager 2023 Methods, "PK sampling and
    # assay"). No conversion factor is needed.
    # ==================================================================

    # ----- Disposition -----
    lcl <- log(0.0960 / 24) ; label("Elimination clearance at 83.6 kg body weight (CL, L/h)")     # Sager 2023 Table 1: 0.0960 L/day (RSE 1.33%); /24 -> 0.00400 L/h
    lvc <- log(3.33)        ; label("Central volume of distribution (V2, L)")                     # Sager 2023 Table 1: 3.33 L (RSE 2.20%)
    lq  <- log(0.667 / 24)  ; label("Distribution clearance (Q, L/h)")                            # Sager 2023 Table 1: 0.667 L/day (RSE 1.49%); /24 -> 0.0277917 L/h
    lvp <- log(4.51)        ; label("Peripheral volume of distribution at 83.6 kg body weight (V3, L)") # Sager 2023 Table 1: 4.51 L (RSE 1.32%)

    # ----- Intramuscular absorption -----
    # Table 1 prints the KA units as "L/h", which is a typographical
    # error: KA is the first-order absorption rate constant of the depot
    # and must be 1/h. The Results text repeats the same slip ("the rate
    # of absorption following i.m. dosing was 0.00643 L/h"). This is
    # CONFIRMED, not inferred: the NONMEM control stream in the
    # Supporting Information annotates the same THETA as
    # ";--th5- KA: Absorption Rate in Male Participants with BMI of
    # 30.41 m^2 (1/h)" and uses it as `DADT(1) = -KA * A(1)`
    # (MODEL CODE S1, $THETA and $DES). Read as 1/h the value gives an
    # absorption half-life of log(2)/0.00643 = 108 h (4.5 days), the
    # expected order for intramuscular monoclonal-antibody absorption.
    lka <- log(0.00643)     ; label("First-order absorption rate constant from the i.m. depot in male participants with BMI 30.41 kg/m2 (KA, 1/h)") # Sager 2023 Table 1: 0.00643 (RSE 4.67%)
    lr1 <- log(130)         ; label("Zero-order input rate R1 into the i.m. depot (mg/h)")        # Sager 2023 Table 1: 130 mg/h (RSE 6.69%)

    # F_IM is estimated on the logit scale (Sager 2023 Table 1 footnote
    # b defines LFIM = log(F_IM / (1 - F_IM)) + sexf * (-0.449)). The
    # footnote's back-transformation is printed as
    # "F_IM = exp(LFIM)/(1 - exp(LFIM))", which is a sign typo -- that
    # expression is negative for any LFIM > 0 and cannot be a
    # bioavailability. The intended inverse logit is
    # exp(LFIM)/(1 + exp(LFIM)), and this too is CONFIRMED rather than
    # inferred: MODEL CODE S1 of the Supporting Information writes
    # `TVLFIM = log (TVFIM / (1-TVFIM)) + COV4`, `LFIM = TVLFIM + ETA(3)`
    # and `FIM = EXP(LFIM) / (1 + EXP(LFIM))`. The paper's own numbers
    # agree: the male estimate 0.582 with the female shift -0.449 maps
    # onto 0.4705, matching the 0.471 female value quoted in the Results.
    logitfdepot <- log(0.582 / (1 - 0.582)) ; label("Logit of intramuscular bioavailability F_IM in male participants (fraction)") # Sager 2023 Table 1: F_IM = 0.582 in males (RSE 2.64%); logit(0.582) = 0.33096

    # ----- Covariate effects -----
    e_wt_cl          <- 0.494  ; label("Power of body weight on CL, normalised to 83.6 kg (unitless)")             # Sager 2023 Table 1, "Power of body weight effect" under CL (RSE 7.18%)
    e_wt_vp          <- 0.757  ; label("Power of body weight on the peripheral volume V3, normalised to 83.6 kg (unitless)") # Sager 2023 Table 1, "Power of body weight effect" under V3 (RSE 6.08%)
    e_bmi_ka         <- -0.711 ; label("Power of body mass index on KA, normalised to 30.41 kg/m2 (unitless)")     # Sager 2023 Table 1, "Power of BMI effect" under KA (RSE 20.0%)
    e_sexf_ka        <- -0.323 ; label("Proportional shift in KA for female sex vs male reference (fraction)")     # Sager 2023 Table 1, "Proportional shift in female participants" under KA (RSE 11.9%)
    e_sexf_logitfdepot <- -0.449 ; label("Additive shift on the logit of F_IM for female sex vs male reference (unitless logit)") # Sager 2023 Table 1, "Shift in female participants, on logit scale" under F_IM (RSE 16.4%)

    # ----- Interindividual variability: full 5x5 covariance block -----
    # Sager 2023 Table 1 reports the DIAGONAL elements as %CV and the
    # OFF-DIAGONAL elements directly as covariances (the ten "cov(IIV in
    # A, IIV in B)" rows), each with the implied correlation coefficient
    # in footnotes d-m.
    #
    # The %CV convention is sqrt(exp(omega^2) - 1) * 100, NOT omega*100.
    # This is over-determined and verified: recomputing each printed
    # correlation as cov / (omega_a * omega_b) using omegas derived by
    # that formula reproduces all six correlations among CL, V2, V3 and
    # KA to three decimals (0.714 -> 0.7130, -0.122 -> -0.1220,
    # -0.318 -> -0.3175, 0.631 -> 0.6324, 0.227 -> 0.2272,
    # -0.115 -> -0.1148). The alternative omega = %CV/100 convention
    # misses every one of them by a factor of ~1.03-1.09.
    #
    #   CL:  log(1 + 0.382^2) = 0.136211
    #   V2:  log(1 + 0.572^2) = 0.283059
    #   V3:  log(1 + 0.294^2) = 0.082903
    #   KA:  log(1 + 0.554^2) = 0.267670
    #
    # F_IM's diagonal is NOT directly recoverable from its printed
    # "42.9 %CV", because that figure is a delta-method transformation
    # of the logit-scale omega onto the bioavailability scale
    # (Table 1 footnote a: 100 * (1 - 0.582) * omega). It is instead
    # back-solved from the four printed covariance/correlation pairs
    # that involve F_IM, which agree tightly:
    #   from cov(F_IM, CL): 0.180 / 0.480 / 0.369068 = 1.01607
    #   from cov(F_IM, V2): 0.317 / 0.588 / 0.532033 = 1.01331
    #   from cov(KA, F_IM): 0.303 / 0.577 / 0.517368 = 1.01500
    #   from cov(V3, F_IM): 0.00744 / 0.0255 / 0.287928 = 1.01332
    # Mean omega = 1.01443, so the variance is 1.029064. That value
    # reproduces all four printed correlations to within 0.0008
    # (0.4808 vs 0.480, 0.5874 vs 0.588, 0.0255 vs 0.0255, 0.5773 vs
    # 0.577). Footnote a's rounded "1.03" would instead need omega =
    # 42.9/41.8 = 1.02632, which misses each printed correlation by
    # ~1.2%; the covariance block is the more reliable route and is the
    # one used here. See the vignette's Assumptions and deviations.
    #
    # The resulting 5x5 matrix is positive definite (smallest eigenvalue
    # 6.02e-05).
    etalcl + etalvc + etalvp + etalka + etalogitfdepot ~
      c(0.136211,
        0.140000, 0.283059,
        0.067200, 0.034800, 0.082903,
       -0.023300, -0.087400, -0.017100, 0.267670,
        0.180000, 0.317000, 0.007440, 0.303000, 1.029064)

    # ----- Residual variability: combined additive + CCV -----
    # Table 1 reports the two components as VARIANCES (0.0175 and
    # 0.0312); footnote n gives the SD as
    # sqrt(0.0175 * F^2 + 0.0312), i.e. the square roots are the
    # proportional and additive SDs. Confirmed by the footnote's own
    # bracketed range: at F = 0.1 ug/mL the %CV is 177.1 and at
    # F = 4000 ug/mL it is 13.2, matching the printed "177-13.2 %CV
    # F [0.1-4000]". Also consistent with the base-model narrative
    # ("%CV < 14% for concentrations > 10 ug/mL"): at F = 10 the %CV is
    # 13.3.
    propSd <- sqrt(0.0175) ; label("Proportional (CCV) residual error (fraction)")  # Sager 2023 Table 1: CCV component 0.0175 (RSE 0.693%); sqrt = 0.132288
    addSd  <- sqrt(0.0312) ; label("Additive residual error (ug/mL)")               # Sager 2023 Table 1: additive component 0.0312 (RSE 10.3%); sqrt = 0.176635
  })

  model({
    # 1. Individual parameters. Body weight is normalised to the
    #    population median 83.6 kg and BMI to 30.41 kg/m2, the reference
    #    values named in the Table 1 row labels.
    #    Q carries no random effect: MODEL CODE S1 writes `Q = TVQ`,
    #    with no EXP(ETA(.)) factor, and Table 1 reports its magnitude of
    #    variability as "NE" (not estimated).
    cl <- exp(lcl + etalcl) * (WT / 83.6)^e_wt_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp + etalvp) * (WT / 83.6)^e_wt_vp

    # KA carries a BMI power effect and a proportional female shift. The
    # shift is multiplicative-on-one, not exponential: Table 1 labels it
    # "Proportional shift in female participants" and MODEL CODE S1
    # writes `COV5 = (1+THETA(12)*SEXF)` with `TVKA = THETA(5)*COV3*COV5`.
    ka <- exp(lka + etalka) * (BMI / 30.41)^e_bmi_ka * (1 + e_sexf_ka * SEXF)

    # F_IM: additive sex shift on the logit scale, then inverse logit.
    fdepot <- expit(logitfdepot + e_sexf_logitfdepot * SEXF + etalogitfdepot)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment disposition with a depot for the i.m. route.
    #    Intravenous doses are given directly into `central`, so the
    #    depot is simply unused for those records.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. The authors' "sigmoid absorption model" for the i.m. route: a
    #    zero-order input into the depot at rate R1 followed by
    #    first-order absorption KA out of it. rxode2 reads rate() at
    #    solve time and spreads the (bioavailability-scaled) dose over
    #    amt * fdepot / r1 hours, matching NONMEM's modelled-RATE
    #    semantics. For the 500 mg i.m. dose in a typical male that is
    #    500 * 0.582 / 130 = 2.24 h. Zero-order delivery into a
    #    compartment that empties first-order is what makes the
    #    cumulative absorption S-shaped rather than exponential.
    f(depot)    <- fdepot
    rate(depot) <- exp(lr1)

    # 5. Observation. mg / L = ug/mL, the assay units.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
