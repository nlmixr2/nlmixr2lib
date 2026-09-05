Chandasana_2024b_dolutegravir <- function() {
  description <- "One-compartment oral population PK model with first-order absorption, allometric weight scaling and a postmenstrual-age enzyme-maturation function on clearance for dolutegravir in infants, children and adolescents with HIV-1, applied without re-estimation to the ABC/DTG/3TC fixed-dose combination (dispersible tablet and tablet) in IMPAACT 2019 (Chandasana 2024)"
  reference <- paste(
    "Chandasana H, van Dijkman SC, Mehta R, Bush M, Rabie H, Flynn P,",
    "Cressey TR, Acosta EP, Brooks KM; for the IMPAACT 2019 Study Team.",
    "Population pharmacokinetic modeling of abacavir/dolutegravir/lamivudine",
    "to support a fixed-dose combination in children with HIV-1.",
    "Infect Dis Ther. 2024;13(8):1877-1891. doi:10.1007/s40121-024-01008-y.",
    "The dolutegravir model itself was developed in Chandasana H, Thapar M,",
    "Hayes S, Baker M, Gibb DM, Turkova A, et al. Population pharmacokinetic",
    "modeling of dolutegravir to optimize pediatric dosing in HIV-1-infected",
    "infants, children, and adolescents. Clin Pharmacokinet. 2023;62(10):1445-1459",
    "(Chandasana 2024 reference 15) and is reproduced in Chandasana 2024 Table 2;",
    "Chandasana 2024 applied it to IMPAACT 2019 without re-estimation",
    "(NONMEM MAXEVAL = 0 external validation).",
    "That primary paper is now itself extracted as",
    "modellib('Chandasana_2023_dolutegravir'), and the parameter values here",
    "were corrected against it: Chandasana 2024 Table 2 reproduces the",
    "formulation effects as ratios and the residual errors as variances, and",
    "defers the random-effects covariance matrix to the primary.",
    sep = " "
  )
  vignette <- "Chandasana_2024_abacavir_dolutegravir_lamivudine_pediatric"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power effect on CL/F (exponent 0.455) and V/F (exponent 0.556) with reference weight 70 kg (Chandasana 2024 Table 2). Both exponents were estimated, not fixed at the canonical 0.75 / 1.0. Weight range 3.9-91.0 kg in the pooled model-development population; 8.15-39.30 kg in the IMPAACT 2019 external-validation cohort.",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the sigmoidal enzyme-maturation function on CL/F, FMAT = PMA^HILL / (PMA^HILL + TM50^HILL), with TM50 = 52.2 postmenstrual weeks (equivalently 12 weeks postnatal age) and Hill = 3.43, both FIXED (Chandasana 2024 Table 2 and its footnote, citing Anderson and Holford 2009). The source model is written in postmenstrual WEEKS with PMA (weeks) = PNA (years) * 52 + 40; the canonical PAGE column is in months, so model() converts months to weeks with 1 month = 4.348125 weeks before evaluating FMAT. FMAT saturates to 1 in older children, so the maturation term is only material in infants.",
      source_name        = "PMA (weeks)"
    ),
    FORM_DTG_DT = list(
      description        = "Dolutegravir dispersible-tablet / granule formulation indicator (1 = dispersible tablet or granules for oral suspension, 0 = film-coated tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (film-coated tablet, FCT)",
      notes              = "The source model accounts for differences in absorption rate and bioavailability across dolutegravir formulations (Chandasana 2024 'DTG Pediatric PopPK Model'). Both reported effects are RATIOS to the film-coated tablet: Chandasana 2024 Table 2's row 'KA ~ DT and granules = 2.04' is the multiplier THETA(17) of the primary source's control stream (KAFORM = THETA(17)**FORMK), so the dispersible-tablet Ka is 0.854 x 2.04 = 1.74 1/h, which Chandasana 2023 Table S1 tabulates as 1.74 (95% CI 1.20-2.28); and relative bioavailability is 1.00 for the fasted film-coated tablet versus 1.53 for the fasted dispersible tablet / granules. The dispersible tablet and granules share a single set of estimates because a healthy-adult study found them bioequivalent at the same dose, so one indicator covers both. Named in the FORM_<drug>_<formulation> family because the pooled DT + granule contrast is dolutegravir-formulation-development specific; promote to a general FORM_DISPERSIBLE_TABLET canonical if a second paper ratifies the same encoding.",
      source_name        = "formulation"
    ),
    FED = list(
      description        = "Dose administered without regard to food (1 = dosed without regard to food, 0 = dosed fasted)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Chandasana 2024 Table 2 reports relative bioavailability 1.00 for the fasted film-coated tablet and 1.10 for the film-coated tablet given without regard to food. The source category is 'without regard to food' (i.e. dosing was not required to be fasted) rather than a controlled fed challenge with a defined meal, so the general FED indicator applies rather than FED_HIGHFAT. The effect is NOT formulation-specific: the primary source's control stream forms F1 = 1 * THETA(4)**SFLAG * THETA(9)**FFLAG, whose own comment block enumerates fasted FCT = 1, fasted DT = THETA(9), fed FCT = THETA(4) and fed DT = THETA(9) * THETA(4), and Chandasana 2023 Table S1 accordingly tabulates F for the dispersible tablet without regard to food as 1.68 = 1.10 x 1.53. Chandasana 2024 Table 2 lists a fed estimate on the film-coated-tablet row only because that row is the reference level.",
      source_name        = "food"
    ),
    STUDY_ODYSSEY = list(
      description        = "ODYSSEY study indicator (1 = record from the ODYSSEY trial, 0 = record from IMPAACT P1093)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (IMPAACT P1093)",
      notes              = "The source model is a pooled analysis of two clinical studies (Chandasana 2024 'DTG Pediatric PopPK Model') with study-specific residual error. Chandasana 2024 Table 2's additive rows print the NONMEM $SIGMA VARIANCES 0.00164 (P1093) and 0.0900 (ODYSSEY); nlmixr2's add() takes the standard-deviation scale, so the values used here are their square roots, which Chandasana 2023 Table 2 prints explicitly as SD = 0.0405 ug/mL and SD = 0.300 ug/mL. The residual error is therefore proportional 28.6% plus additive 0.0405 ug/mL for P1093, and proportional 11.1% plus additive 0.300 ug/mL for ODYSSEY. Both sets are retained here and selected inside model() by this indicator. IMPAACT 2019 records are neither study; set STUDY_ODYSSEY = 0 to use the P1093 residual error, which is the pediatric dolutegravir single-entity study most comparable to IMPAACT 2019.",
      source_name        = "study"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion variability on CL/F and Ka",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Chandasana 2024 Table 2 reports a single inter-occasion variability magnitude per parameter (CL/F 33.9% CV, Ka 91.7% CV) shared across occasions, i.e. the NONMEM $OMEGA BLOCK(1) SAME idiom. Two occasions are encoded here (the minimum that makes IOV operational); set OCC = 1 for every record to reproduce the single steady-state occasion simulated in Chandasana 2024 Table 4. Users needing more occasions extend the oc<k> / etaiov_*_<k> pattern with additional ~ fix(...) slots.",
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 239L,
    n_studies      = 2L,
    age_range      = "2 months to 18 years",
    weight_range   = "3.9-91.0 kg",
    disease_state  = "Infants, children and adolescents living with HIV-1 receiving oral dolutegravir",
    dose_range     = "Oral dolutegravir; in the IMPAACT 2019 confirmatory simulations 15 mg (>=6 to <10 kg), 20 mg (>=10 to <14 kg), 25 mg (>=14 to <20 kg) and 30 mg (>=20 to <25 kg) once daily as the ABC/DTG/3TC dispersible tablet, and 50 mg once daily as the ABC/DTG/3TC tablet (>=25 to <40 kg)",
    notes          = "Model-development population: a pooled analysis of two clinical studies (IMPAACT P1093 and ODYSSEY) reported in Chandasana et al. Clin Pharmacokinet 2023;62(10):1445-1459 (Chandasana 2024 reference 15), summarised in Chandasana 2024 'DTG Pediatric PopPK Model' and Table 2. External-validation cohort: IMPAACT 2019 (NCT03760458), an international phase I/II open-label study in children <12 years living with HIV-1 enrolled into five weight bands (>=6 to <10, >=10 to <14, >=14 to <20, >=20 to <25 and >=25 to <40 kg); 55 participants contributed 598 dolutegravir intensive and sparse PK samples, median (min-max) baseline age 6.0 (1.00-11.0) years and weight 17.00 (8.15-39.30) kg, 45.5% female, 67% Black and 31% Asian (Chandasana 2024 Results). The existing model was applied to the IMPAACT 2019 data with NONMEM MAXEVAL = 0, i.e. no parameter was re-estimated. Predefined exposure targets for dose confirmation were a geometric-mean C24 of 0.697-2.26 ug/mL and a geometric-mean AUC0-24 of 37-134 ug*h/mL (Chandasana 2024 Methods)."
  )

  ini({
    # Structural parameters at the reference body weight of 70 kg, fully mature
    # (FMAT = 1), fasted film-coated tablet. All values from Chandasana 2024
    # Table 2 (previously reported DTG pediatric PopPK parameter estimates;
    # primary source Chandasana 2023 Clin Pharmacokinet, reference 15).
    lka          <- log(0.854); label("Absorption rate constant Ka, film-coated tablet (1/h)")            # Chandasana 2024 Table 2: KA, FCT = 0.854 1/h (%RSE 11.2)
    e_form_dt_ka <- log(2.04);  label("log ratio of the dispersible-tablet / granule absorption rate constant to the film-coated-tablet value (unitless; exp = 2.04)")  # Chandasana 2024 Table 2 row "KA ~ DT and granules" = 2.04 (%RSE 15.7) is the RATIO THETA(17) of the primary source's control stream (KAFORM = THETA(17)**FORMK), not an absolute Ka; the dispersible-tablet Ka is 0.854 x 2.04 = 1.74 1/h, as Chandasana 2023 Table S1 tabulates (1.74, 95% CI 1.20-2.28).
    lcl     <- log(1.03);  label("Apparent oral clearance CL/F at WT = 70 kg and full maturation (L/h)")  # Chandasana 2024 Table 2: CL/F = 1.03 L/h (%RSE 2.31)
    lvc     <- log(13.6);  label("Apparent volume of distribution V/F at WT = 70 kg (L)")                 # Chandasana 2024 Table 2: V/F = 13.6 L (%RSE 2.42)

    # Allometric weight exponents; both estimated (not fixed at 0.75 / 1.0).
    e_wt_cl <- 0.455; label("Power exponent of body weight on CL/F, reference 70 kg (unitless)")  # Chandasana 2024 Table 2: CL/F~(WT/70) = 0.455 (%RSE 4.15)
    e_wt_vc <- 0.556; label("Power exponent of body weight on V/F, reference 70 kg (unitless)")   # Chandasana 2024 Table 2: V/F~(WT/70)  = 0.556 (%RSE 3.87)

    # Sigmoidal enzyme-maturation function on CL/F, in postmenstrual WEEKS:
    #   FMAT = PMA^hill_mat / (PMA^hill_mat + tmat50^hill_mat)
    # Both parameters reported as FIX in Chandasana 2024 Table 2.
    tmat50   <- fixed(52.2); label("Postmenstrual age at half-maximal CL/F maturation TM50 (weeks)")   # Chandasana 2024 Table 2: maturation half time TM50 = 52.2 PMA weeks FIX (12 weeks postnatal age)
    hill_mat <- fixed(3.43); label("Hill coefficient of the CL/F enzyme-maturation function (unitless)") # Chandasana 2024 Table 2: Hill coefficient = 3.43 FIX

    # Formulation and food effects on relative bioavailability, encoded on the
    # log scale so model() uses exp(effect * indicator). Reference is the fasted
    # film-coated tablet (relative bioavailability 1.00, Chandasana 2024 Table 2).
    e_food_fdepot    <- log(1.10); label("log relative bioavailability of the film-coated tablet given without regard to food (unitless; exp = 1.10)") # Chandasana 2024 Table 2: F, without regard to food FCT = 1.10 (%RSE 3.03)
    e_form_dt_fdepot <- log(1.53); label("log relative bioavailability of the fasted dispersible tablet / granules (unitless; exp = 1.53)")            # Chandasana 2024 Table 2: F, fasted DT/granules = 1.53 (%RSE 3.26)

    # IIV - a full 3x3 NONMEM $OMEGA BLOCK(3) on log-CL/F, log-V/F and log-Ka.
    # Chandasana 2024 Table 2 reports only the CV% diagonal and defers the
    # covariance matrix to its reference [15] ("Further details about covariance
    # matrix and full model can be found in the reference [15]"). That primary
    # source is now on disk, so the variances and covariances below are taken
    # directly from Chandasana 2023 Table 2 rather than back-solved from the
    # rounded CV% column. Every value reproduces its printed CV% and correlation.
    etalcl + etalvc + etalka ~ c(0.0863,
                                 0.0499, 0.0698,
                                 0.0953, 0.138,  0.762)
    # Chandasana 2023 Table 2: omega^2 CL = 0.0863 (CV 29.4%); Covar eta_CL,eta_V = 0.0499 (R = 0.643);
    #   omega^2 V = 0.0698 (CV 26.4%); Covar eta_CL,eta_KA = 0.0953 (R = 0.372);
    #   Covar eta_V,eta_KA = 0.138 (R = 0.598); omega^2 KA = 0.762 (CV 107%).
    #   Chandasana 2023 ESM Supplementary Text 1 control stream: $OMEGA BLOCK(3).

    # IOV on CL/F and Ka - a single shared magnitude across occasions (NONMEM
    # $OMEGA BLOCK(1) SAME); occasions after the first are fixed to the first
    # occasion's variance.
    etaiov_lcl_1 ~ 0.10878        # Chandasana 2024 Table 2: IOV-CL/F = 33.9% CV; log(1 + 0.339^2) = 0.10878
    etaiov_lcl_2 ~ fix(0.10878)   # Chandasana 2024 Table 2: same IOV-CL/F magnitude on every occasion
    etaiov_lka_1 ~ 0.61025        # Chandasana 2024 Table 2: IOV-KA   = 91.7% CV; log(1 + 0.917^2) = 0.61025
    etaiov_lka_2 ~ fix(0.61025)   # Chandasana 2024 Table 2: same IOV-KA magnitude on every occasion

    # Residual error - combined additive + proportional, study-specific.
    # Both pooled studies' estimates are retained and selected inside model()
    # by STUDY_ODYSSEY.
    # The additive rows of Chandasana 2024 Table 2 carry NONMEM $SIGMA
    # VARIANCES; nlmixr2's add() takes the standard-deviation scale. Chandasana
    # 2023 Table 2 prints both, giving SD = sqrt(0.00164) = 0.0405 ug/mL for
    # P1093 and SD = sqrt(0.0900) = 0.300 ug/mL for ODYSSEY, and its control
    # stream confirms the error model is Y = F + F*ERR(prop) + ERR(add).
    propSd          <- 0.286;   label("Proportional residual error, IMPAACT P1093 (fraction)")  # Chandasana 2024 Table 2: proportional error study P1093 = 28.6% (variance 0.0818)
    addSd           <- 0.0405;  label("Additive residual error, IMPAACT P1093 (ug/mL)")         # Chandasana 2023 Table 2: additive error P1093, variance 0.00164, SD = 0.0405 ug/mL
    propSd_odyssey  <- 0.111;   label("Proportional residual error, ODYSSEY (fraction)")        # Chandasana 2024 Table 2 continued: proportional error ODYSSEY = 11.1% (variance 0.0123)
    addSd_odyssey   <- 0.300;   label("Additive residual error, ODYSSEY (ug/mL)")               # Chandasana 2023 Table 2: additive error ODYSSEY, variance 0.0900, SD = 0.300 ug/mL
  })

  model({
    # Inter-occasion variability, multiplexed by the occasion indicator.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2
    iov_ka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2

    # Enzyme maturation on CL/F. The source model is written in postmenstrual
    # WEEKS (Chandasana 2024 Table 2 footnote: PMA (weeks) = PNA (years) * 52 +
    # 40); the canonical PAGE column is postmenstrual age in MONTHS, so convert
    # with 1 month = 30.4375 / 7 = 4.348125 weeks, matching the PAGE schema in
    # inst/references/covariate-columns.md.
    pma_weeks <- PAGE * 4.348125
    fmat <- pma_weeks^hill_mat / (pma_weeks^hill_mat + tmat50^hill_mat)

    # Formulation- and food-dependent relative bioavailability. Reference is the
    # fasted film-coated tablet. The two ratios multiply and the food effect is
    # NOT formulation-specific: the primary source's control stream forms
    # F1 = 1 * THETA(4)**SFLAG * THETA(9)**FFLAG, and Chandasana 2023 Table S1
    # tabulates F for the dispersible tablet without regard to food as the
    # product 1.68 = 1.10 x 1.53. Chandasana 2024 Table 2 lists a fed estimate
    # on the film-coated-tablet row only because that row is the reference level.
    frel <- exp(e_food_fdepot * FED + e_form_dt_fdepot * FORM_DTG_DT)

    # Individual PK parameters. Reference body weight 70 kg, full maturation.
    # The formulation effect on Ka is a multiplicative ratio, giving
    # 0.854 x 2.04 = 1.74 1/h for the dispersible tablet / granules.
    ka <- exp(lka + e_form_dt_ka * FORM_DTG_DT + etalka + iov_ka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * fmat
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    f(depot) <- frel

    # Study-specific residual error magnitudes.
    propSd_i <- propSd * (1 - STUDY_ODYSSEY) + propSd_odyssey * STUDY_ODYSSEY
    addSd_i  <- addSd  * (1 - STUDY_ODYSSEY) + addSd_odyssey  * STUDY_ODYSSEY

    Cc <- central / vc
    Cc ~ add(addSd_i) + prop(propSd_i)
  })
}
