Zou_2022_pretomanid <- function() {
  description <- "One-compartment population PK model with one transit compartment and first-order absorption for pretomanid dispersible (pediatric) and marketed tablets given with a high-fat meal to healthy adults, with dose-dependent relative bioavailability and inter-occasion variability on absorption"
  reference <- "Zou Y, Nedelman J, Lombardi A, Pappas F, Karlsson MO, Svensson EM. Characterizing Absorption Properties of Dispersible Pretomanid Tablets Using Population Pharmacokinetic Modelling. Clin Pharmacokinet. 2022;61(11):1585-1593. doi:10.1007/s40262-022-01163-w"
  vignette <- "Zou_2022_pretomanid"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")
  # Oral doses enter the transit compartment, not `depot`: `depot` here is the
  # source control stream's ABS compartment, which sits downstream of the
  # transit compartment and empties into `central` at KA (ESM Section C
  # $MODEL: COMP(TRANS1,DEFDOSE), K13 = KTR, K32 = KA). Same layout as
  # Salinger_2019_pretomanid, the starting model.
  dosing <- "transit1"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling on CL (exponent 0.75) and V (exponent 1), both fixed and normalised to a 55 kg reference participant (Table 2 footnote a), carried over from the Salinger 2019 starting model. The 55 kg reference lies below the enrolled fed-panel range (64.4-117 kg).",
      source_name = "WT"
    ),
    DOSE_PRETOMANID_MG = list(
      description = "Administered pretomanid oral dose level",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-function effect on relative bioavailability normalised to 200 mg: F = (DOSE/200)^0.0822 (Table 2; ESM control stream TVLF1 = LOG(DOSE/200) * THETA(6)). Evaluated on the dose record. Studied levels were 10, 50 and 200 mg; the 200 mg dispersible-tablet dose was given as 4 x 50 mg tablets.",
      source_name = "DOSE"
    ),
    FORM_PRETOMANID_DT = list(
      description = "Pretomanid pediatric dispersible-tablet formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (dispersible tablet) is the model reference for KA; 0 = marketed 200 mg tablet",
      notes = "The source control stream forms MF = 1 - DTF and multiplies KA by 1.65 for the marketed formulation. Formulation was tested on bioavailability (ratio 1.00, 90% CI 0.87-1.14) but not retained.",
      source_name = "DTF"
    ),
    OCC = list(
      description = "Dosing-period occasion (1-4) for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "A new occasion starts at each of the four single doses of the crossover, which were separated by a 7-day washout (Methods 'Pharmacokinetic Analysis'). The control stream carries OCC1-OCC4 indicator columns; here they are derived from OCC inside model(). Values outside 1-4 switch inter-occasion variability off.",
      source_name = "OCC1, OCC2, OCC3, OCC4"
    )
  )

  covariatesDataExcluded <- list(
    FED = list(
      description = "Fed-versus-fasted state at the dose record",
      units = "(binary)",
      type = "binary",
      notes = "Only the fed panel (FDA standard high-fat, high-calorie breakfast) was modelled; every record in the analysis was fed, so no food effect is estimable from this model. The fasted panel was analysed by NCA only (ESM Tables D2, D4, D5)."
    )
  )

  compartmentData <- list(
    transit1 = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = TRUE),
    depot = list(analyte = "pretomanid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pretomanid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    n_observations = 1377L,
    age_range = "23-50 years",
    age_median = "39.0 years",
    weight_range = "64.4-117 kg",
    weight_median = "76.2 kg",
    height_median = "170 cm (range 155-193)",
    bmi_median = "27.0 kg/m^2 (range 22.6-31.5)",
    sex_female_pct = 33.3,
    disease_state = "healthy adult volunteers",
    dose_range = "single oral doses of 10, 50 and 200 mg (4 x 50 mg) dispersible tablet and 200 mg marketed tablet, each after an FDA standard high-fat, high-calorie breakfast, in a four-period crossover with 7-day washouts",
    regions = "United States (San Antonio, Texas)",
    notes = "Phase 1 relative-bioavailability and food-effect study CL-011 (NCT04309656). The model uses the fed panel only (24 of 48 enrolled participants; Table 1 'Fed' column). Observations below the 1 ng/mL quantification limit (n = 39, 2.8%) were excluded."
  )

  ini({
    # ---- Absorption ------------------------------------------------------
    # One transit compartment (rate KTR = 1/MTT) feeds a first-order
    # absorption compartment (rate KA); ESM Section C $MODEL / $PK.
    lka <- log(0.396)
    label("First-order absorption rate constant for the dispersible tablet (1/h)") # Table 2 'Absorption rate of DTF (KA, h-1)' = 0.396 (RSE 14%); ESM THETA(3)
    e_form_mf_ka <- 1.65
    label("Fold change in KA for the marketed formulation relative to the dispersible tablet (unitless)") # Table 2 'Proportional effect (theta3) of MF on KA' = 1.65 (RSE 26%); ESM THETA(5)
    lmtt <- log(1.13)
    label("Mean transit time through the transit compartment (h)") # Table 2 'Mean transit time (MTT, h)' = 1.13 (RSE 23%); ESM THETA(4)

    # ---- Disposition -----------------------------------------------------
    lcl <- log(2.81)
    label("Apparent oral clearance for a 55 kg reference participant (L/h)") # Table 2 'Apparent clearance (CL/F, L/h)' = 2.81 (RSE 14%); ESM THETA(1)
    lvc <- log(68.0)
    label("Apparent volume of distribution for a 55 kg reference participant (L)") # Table 2 'Apparent volume of distribution (Vd/F, L)' = 68.0 (RSE 3.7%); ESM THETA(2)
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT/55) for CL (unitless)") # Table 2 'Weight scaling (theta1) on CL' = 0.75 (fixed)
    e_wt_vc <- fixed(1)
    label("Allometric exponent on (WT/55) for V (unitless)") # Table 2 'Weight scaling (theta2) on Vd' = 1 (fixed)

    # ---- Relative bioavailability ----------------------------------------
    lfdepot <- fixed(log(1))
    label("Relative bioavailability for a 200 mg dose (unitless)") # Table 2 'Bioavailability (F)' = 1 (fixed)
    e_dose_fdepot <- 0.0822
    label("Power exponent on (DOSE/200) for relative bioavailability (unitless)") # Table 2 'Dose effect (theta4) on F' = 0.0822 (RSE 19%); ESM THETA(6)

    # ---- Interindividual variability -------------------------------------
    # ESM $OMEGA BLOCK(2) over CL and V; Table 2 prints the CVs as
    # sqrt(variance) x 100 (footnote b): sqrt(0.0622) = 24.9%,
    # sqrt(0.00751) = 8.67%. The covariance is printed only in the ESM.
    etalcl + etalvc ~ c(
      0.0622,
      0.0111, 0.00751
    ) # ESM $OMEGA BLOCK(2) 0.0622 / 0.0111 / 0.00751; Table 2 IIV CL 24.9% (RSE 32%), Vd 8.67% (RSE 34%)

    # ---- Interoccasion variability ---------------------------------------
    # ESM $OMEGA BLOCK(3) on F1, KA and MTT with zero off-diagonals, then
    # three BLOCK(3) SAME blocks: one draw per occasion with the same
    # variances, so occasions 2-4 are fixed to the occasion-1 values.
    etaiov_fdepot_1 ~ 0.00561 # ESM $OMEGA 'IOV.F1' = 0.00561; Table 2 IOV F 7.49% (RSE 13%)
    etaiov_fdepot_2 ~ fixed(0.00561) # same variance as occasion 1
    etaiov_fdepot_3 ~ fixed(0.00561) # same variance as occasion 1
    etaiov_fdepot_4 ~ fixed(0.00561) # same variance as occasion 1
    etaiov_ka_1 ~ 0.284 # ESM $OMEGA 'IOV.KA' = 0.284; Table 2 IOV KA 53.3% (RSE 26%)
    etaiov_ka_2 ~ fixed(0.284) # same variance as occasion 1
    etaiov_ka_3 ~ fixed(0.284) # same variance as occasion 1
    etaiov_ka_4 ~ fixed(0.284) # same variance as occasion 1
    etaiov_mtt_1 ~ 0.77 # ESM $OMEGA 'IOV.MATT' = 0.77; Table 2 IOV MTT 87.8% (RSE 20%)
    etaiov_mtt_2 ~ fixed(0.77) # same variance as occasion 1
    etaiov_mtt_3 ~ fixed(0.77) # same variance as occasion 1
    etaiov_mtt_4 ~ fixed(0.77) # same variance as occasion 1

    # ---- Residual error --------------------------------------------------
    # ESM $ERROR: Y = IPRED + ERRT * (EPS(1) * IPRED + EPS(2)), diagonal
    # $SIGMA, so the proportional and additive variances add (combined2).
    propSd <- 0.0889
    label("Proportional residual SD from 10 h after the dose (fraction)") # Table 2 'Proportional error' = 8.89% (RSE 5.2%); ESM $SIGMA 0.0079, sqrt = 0.0889
    addSd <- 0.401
    label("Additive residual SD from 10 h after the dose (ng/mL)") # Table 2 'Additive error (ng/mL)' = 0.401 (RSE 3.6%); ESM $SIGMA 0.161, sqrt = 0.401
    e_tad_early_ruv <- 3.46
    label("Multiplicative factor on both residual SDs for samples taken less than 10 h after the dose (unitless)") # Table 2 'Time-varying error term (theta5)' = 3.46 (RSE 6.6%); ESM THETA(7) ERRT
  })

  model({
    # ---- 1. Derived covariate terms --------------------------------------
    mf <- 1 - FORM_PRETOMANID_DT
    doseRatio <- DOSE_PRETOMANID_MG / 200

    # ---- 2. Inter-occasion random effects --------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4
    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 + oc4 * etaiov_mtt_4

    # ---- 3. Individual parameters ----------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 55)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 55)^e_wt_vc
    ka <- exp(lka + iov_ka) * e_form_mf_ka^mf
    mtt <- exp(lmtt + iov_mtt)
    fdepot <- exp(lfdepot + iov_fdepot) * doseRatio^e_dose_fdepot

    # ---- 4. Micro-constants ----------------------------------------------
    # One transit compartment, so KTR = 1/MTT (ESM KTR = 1/MTT, K13 = KTR).
    # Mean absorption time = 1/KA + MTT.
    ktr <- 1 / mtt
    kel <- cl / vc

    # ---- 5. ODE system ---------------------------------------------------
    d/dt(transit1) <- -ktr * transit1
    d/dt(depot) <- ktr * transit1 - ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(transit1) <- fdepot

    # ---- 6. Observation and residual error -------------------------------
    # central is in mg and vc in L, so central/vc is mg/L; x 1000 gives
    # ng/mL, the units of the additive residual SD.
    Cc <- 1000 * central / vc

    # ESM $ERROR: ERRT = THETA(7) when TAD < 10 h, else 1; it scales both
    # residual components. tad() is the time since the most recent dose.
    tadNow <- tad()
    errt <- 1 + (e_tad_early_ruv - 1) * (tadNow < 10)
    propSdCc <- propSd * errt
    addSdCc <- addSd * errt
    Cc ~ add(addSdCc) + prop(propSdCc) + combined2()
  })
}
