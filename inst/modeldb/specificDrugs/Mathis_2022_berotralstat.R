Mathis_2022_berotralstat <- function() {
  description <- "Three-compartment population PK model with first-order absorption, an absorption lag time and linear elimination for the oral plasma kallikrein inhibitor berotralstat in healthy adults and adults and adolescents with hereditary angioedema (Mathis 2022). Distribution is parameterized directly in the first-order rate constants K23, K32, K24 and K42, as the authors did. Clearance and central volume carry power effects of body weight centered on 75.2 kg, and relative bioavailability rises with dose as a power of dose over 300 mg. Between-subject variability on clearance, central volume, ka and the lag time plus a between-study random effect on clearance are the control-stream $OMEGA values, which the paper prints only as initial estimates; the final variances were not published."
  reference <- "Mathis A, Sale M, Cornpropst M, Sheridan WP, Ma SC. Population pharmacokinetic modeling and simulations of berotralstat for prophylactic treatment of attacks of hereditary angioedema. Clin Transl Sci. 2022;15(4):1027-1035. doi:10.1111/cts.13233. Fixed-effect estimates from Table 1; model code and random-effect values from Supporting Information Text S2."
  vignette <- "Mathis_2022_berotralstat"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "berotralstat", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "berotralstat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "berotralstat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "berotralstat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on CL and V2, centered on 75.2 kg (Mathis 2022 Table 1 footnote: 'Clearance and volume were estimated for a 75.2 kg subject (the median weight across all studies)'; Supporting Information Text S2: 'CWTKG =WTKG/75.2', 'TVCL= THETA(1)*CWTKG**THETA(15)...', 'TVV2=THETA(2)*CWTKG**THETA(13)...'). Observed range 40.1-150 kg (Results, Data summary).",
      source_name = "WTKG"
    ),
    DOSE_BEROTRALSTAT_MG = list(
      description = "Administered berotralstat dose on the dose record",
      units = "mg",
      type = "continuous",
      reference_category = "300 mg (the median dose; relative bioavailability = 1)",
      notes = "Enters relative bioavailability as (DOSE_BEROTRALSTAT_MG / 300)^0.497 (Mathis 2022 Equation 1; Supporting Information Text S2: 'CDOSE = 300', 'F1=1*(DOSE/CDOSE)**THETA(11)*(1+FOODONF)', with the food term fixed to 0 in the final model). Set it to the same value as the dose record's amt. The paper's studies record API-in-capsule doses as the HCl salt and blend-in-capsule doses as the free base (Supporting Information Table S1 note); the commercial 150 mg blend-in-capsule dose is free base.",
      source_name = "DOSE"
    )
  )

  # Screened during forward addition / backward elimination but NOT in the
  # final model (Mathis 2022 Methods, Covariate model development; Results,
  # Final model and model evaluation). The Text S2 control stream keeps each of
  # them as a 0 FIX THETA.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested on CL and V2 (centered on 36 years, Text S2 'CAGE =AGE/36'); retained after forward addition, then removed from V2 at backward elimination and dropped from CL from the final model. THETA(17) V~AGE and THETA(19) CL~AGE are 0 FIX in Text S2."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Pre-specified effect on CL; not retained at forward addition."
    ),
    RACE_WHITE = list(
      description = "White race indicator (the paper tested 'race other than White')",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL and V2 as (1 + ARACE*THETA) with ARACE = 1 for any race other than White; removed at backward elimination. THETA(18) and THETA(20) are 0 FIX in Text S2."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (paper column EGFR)",
      units = "mL/min/1.73m^2",
      type = "continuous",
      notes = "Pre-specified effect on CL (centered on 93.6, Text S2 'CEGFR =EGFR/93.6'); not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Pre-specified effect on CL (centered on 19 U/L in Text S2); not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Pre-specified effect on CL (centered on 18 U/L in Text S2); not retained."
    ),
    TBILI = list(
      description = "Total bilirubin (paper column BILI)",
      units = "umol/L",
      type = "continuous",
      notes = "Statistically significant on CL (centered on 0.5 mg/dL = 8.55 umol/L, Text S2 'CBILI =BILI/0.5', the source column being in mg/dL) but removed because it was not clinically significant on the Forest plots and its direction (higher bilirubin, higher clearance) was physiologically implausible. THETA(16) is 0 FIX in Text S2."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Exploratory effect on CL (centered on 4.4 g/dL in Text S2) prompted by a post hoc ETA-versus-covariate trend; not retained."
    ),
    FED = list(
      description = "Fed-state indicator (any food, regular or high-fat)",
      units = "(binary)",
      type = "binary",
      notes = "Tested on ka and on relative bioavailability for the API-in-capsule formulation only; statistically significant but removed because more than 90 percent of the values were imputed and well-controlled studies showed no food effect. THETA(12) and THETA(14) are 0 FIX in Text S2."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 771,
    n_studies = 13,
    n_observations = 10437,
    age_range = "12-74 years",
    age_mean = "38.1 years (SD 13.9); median 37 years (Supporting Information Table S3)",
    weight_range = "40.1-150 kg",
    weight_median = "75.1 kg (mean 77.7 kg, SD 17.2; Supporting Information Table S3); the model is centered on 75.2 kg",
    sex_female_pct = 46.3,
    race_ethnicity = c(
      White = 83.3,
      Asian = 8.7,
      Black = 4.5,
      Other = 2.7,
      `Native American` = 0.4,
      `Pacific Islander` = 0.4
    ),
    disease_state = "Healthy adults (including 7 with severe renal impairment and 18 with Child-Pugh A-C hepatic impairment) and patients with hereditary angioedema (400 of 771, studies BCX7353-203, -204, -301 and -302); 16 subjects were 17 years or younger",
    dose_range = "Oral single doses of 10-900 mg and once-daily repeat doses of 110-450 mg (Methods, Clinical trial subjects); the commercial dose is 150 mg once daily",
    formulation = "API-in-capsule (dose expressed as the HCl salt; studies BCX7353-101, -102, -103, -105, -203) and blend-in-capsule (dose expressed as the free base; the commercial formulation); bioequivalent in study BCX7353-103",
    regions = "International; study BCX7353-301 (APeX-J) enrolled Japanese patients and Part 3 of BCX7353-101 enrolled Japanese healthy subjects",
    notes = "Thirteen phase 1-3 studies pooled (Supporting Information Table S1); demographics in Supporting Information Table S3. Assay LLOQ 1.00 ng/mL; BLQ records were handled with the M3 likelihood in the control stream (Text S2). Final estimates by SAEM in NONMEM 7.3."
  )

  ini({
    # Structural parameters - typical values for a 75.2 kg subject receiving
    # 300 mg (every covariate term equal to 1). Mathis 2022 Table 1.
    lcl <- log(47.3); label("Clearance CL (L/h)")                     # Table 1 'Clearance (L/h)' = 47.3 (RSE 2.4%)
    lvc <- log(1650); label("Central volume of distribution V2 (L)")  # Table 1 'Volume (L)' = 1650 (RSE 2.49%)

    # Distribution is parameterized directly in first-order rate constants, as
    # the authors did (Text S2: 'K23 = THETA(6)' ... 'K42 = THETA(10)'), not in
    # intercompartmental clearances. NONMEM ADVAN12 numbers the compartments
    # 1 = depot, 2 = central, 3 = peripheral1, 4 = peripheral2, so K23 / K32 map
    # to the canonical k12 / k21 and K24 / K42 map to k13 / k31. None carries a
    # weight effect or a random effect.
    lk12 <- log(0.0812);  label("Rate constant, central to peripheral1, K23 (1/h)")   # Table 1 'K23 (1/h)' = 0.0812 (RSE 2.32%)
    lk21 <- log(0.0309);  label("Rate constant, peripheral1 to central, K32 (1/h)")   # Table 1 'K32 (1/h)' = 0.0309 (RSE 2.52%)
    lk13 <- log(0.00281); label("Rate constant, central to peripheral2, K24 (1/h)")   # Table 1 'K24 (1/h)' = 0.00281 (RSE 17.7%)
    lk31 <- log(0.00136); label("Rate constant, peripheral2 to central, K42 (1/h)")   # Table 1 'K42 (1/h)' = 0.00136 (RSE 25.4%)

    # Absorption
    lka   <- log(1.12);  label("Absorption rate constant ka (1/h)")   # Table 1 'Ka (1/h)' = 1.12 (RSE 0.0891%)
    ltlag <- log(0.468); label("Absorption lag time ALAG1 (h)")       # Table 1 'Absorption lag time (h)' = 0.468 (RSE 2.15%)

    # Covariate effects (power models; Text S2 $PK)
    e_dose_f <- 0.497; label("Power exponent on dose relative to 300 mg for relative bioavailability (unitless)")  # Table 1 'Bioavailability as a function of dose' = 0.497 (RSE 5.51%); Equation 1 THETA(11)
    e_wt_vc  <- 1.00;  label("Power exponent on body weight relative to 75.2 kg for central volume (unitless)")  # Table 1 'Volume as a function of weight' = 1.00 (RSE 8.28%)
    e_wt_cl  <- 0.480; label("Power exponent on body weight relative to 75.2 kg for clearance (unitless)")       # Table 1 'Clearance as a function of weight' = 0.480 (RSE 14.3%)

    # Between-subject variability: exponential etas on CL, V2, KA and ALAG1
    # (Text S2: 'CL=TVCL*EXP(ETA(1)+ETA(5))', 'V2=TVV2*EXP(ETA(2))',
    # 'KA=TVKA*EXP(ETA(3))', 'ALAG1 = THETA(8)*EXP(ETA(4))'). The paper does not
    # report the final OMEGA estimates. The values below are the $OMEGA BLOCK(4)
    # INITIAL estimates of the Text S2 final-model control stream, transcribed
    # verbatim (log-scale variances and covariances; the CL-ALAG1 covariance is
    # 0 in the source). Its $THETA initials sit within a few percent of the
    # Table 1 finals, so these are near-converged values, not the finals.
    etalcl + etalvc + etalka + etaltlag ~ c(
      0.141588,
      0.0667178, 0.109354,
      -0.0532533, 0.061502, 0.701395,
      0, -0.0432653, -0.190778, 0.245715
    ) # Text S2 $OMEGA BLOCK(4) initial estimates: ETA(1) CLEARANCE, ETA(2) VOLUME, ETA(3) KA, ETA(4) ALAG1

    # Between-study random effect on CL (Text S2 '$LEVEL STDY=(5[1])', ETA(5)).
    # Also an initial estimate; the final value is unreported. When simulated it
    # is drawn per subject; set it to 0 to reproduce the paper's Table 2, whose
    # geometric CVs are consistent with between-subject variability only.
    eta_study_lcl ~ 0.0817099 # Text S2 $OMEGA initial estimate 0.0817099 'ETA(5) BETWEEN STUDY ON CLEARANCE'

    # Residual error. Text S2 $ERROR: 'W=SQRT(ADDERR**2+PROPERR**2*IPRE**2)',
    # 'Y=IPRED+W*ERR(1)' with $SIGMA 1 FIX, so THETA(3) and THETA(4) are the
    # additive and proportional standard deviations of a combined (variance-sum)
    # error model.
    addSd  <- 0.483; label("Additive residual error (ng/mL)")       # Table 1 'Additive error (ng/ml)' = 0.483 (RSE 3.75%)
    propSd <- 0.286; label("Proportional residual error (fraction)") # Table 1 'Proportional error' = 0.286 (RSE 0.404%)
  })

  model({
    # Covariate model (Text S2 $PK, final-model terms only; the bilirubin, age,
    # race and food THETAs are 0 FIX there and drop out).
    # ETA(1) (between subject) and ETA(5) (between study) both act on CL.
    clbase <- exp(lcl + etalcl)
    cl <- clbase * exp(eta_study_lcl) * (WT / 75.2)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 75.2)^e_wt_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Relative bioavailability rises with dose (Equation 1: F1 = (Dose/300)^THETA(11)).
    f(depot) <- (DOSE_BEROTRALSTAT_MG / 300)^e_dose_f
    alag(depot) <- tlag

    # Doses in mg and volumes in L give mg/L; x1000 converts to ng/mL, the
    # source's 'S2 = V2/1000' scaling.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
