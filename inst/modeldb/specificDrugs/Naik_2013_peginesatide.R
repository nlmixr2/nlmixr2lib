Naik_2013_peginesatide <- function() {
  description <- "Two-compartment population PK/PD model for peginesatide in adult chronic kidney disease (CKD) patients (Naik 2013). PK: first-order subcutaneous absorption with saturable Michaelis-Menten elimination and fixed inter-compartmental clearance. PD: modified precursor-dependent lifespan indirect-response (LIDR) model of hemoglobin (1 progenitor compartment + 7 red-blood-cell aging compartments) with a peginesatide Emax stimulation on progenitor production and an empirical exponential downward-drift factor on the progenitor-to-RBC transit."
  reference <- "Naik H, Tsai MC, Fiedler-Kelly J, Qiu P, Vakilynejad M. A Population Pharmacokinetic and Pharmacodynamic Analysis of Peginesatide in Patients with Chronic Kidney Disease on Dialysis. PLoS ONE. 2013;8(6):e66422. doi:10.1371/journal.pone.0066422"
  vignette <- "Naik_2013_peginesatide"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL", hemoglobin = "g/dL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to convert the paper-reported per-kg volumes (V2 in mL/kg, V3 in mL/kg) and per-kg inter-compartmental clearance (Q in mL/kg/hr) to absolute volumes (L) and flow (L/hr) inside model(). Naik 2013 PK-population mean 79.4 kg (Table 3).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation covariate on Vc (mL/kg per year, reference 59 yr; Naik 2013 eq 14) and exponential-log covariate on the empirical Hgb drift factor CF (per year, reference 58 yr; Naik 2013 eq 17). Paper-reported median 58.4 yr (PK population) / 57.9 yr (PK-PD population) (Table 3).",
      source_name        = "AGE"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on Vc (reference 26 kg/m^2; Naik 2013 eq 14). Paper-reported median BMI 26 kg/m^2 used as the centering value in the discussion.",
      source_name        = "BMI"
    ),
    TBILI = list(
      description        = "Total serum bilirubin at baseline",
      units              = "umol/L (paper-reported as 'g/L', interpreted here as umol/L since the values mean 9.1, range 2-38 are inconsistent with g/L but match clinical-PK total-bilirubin ranges in umol/L)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation covariate on Vc (mL/kg per umol/L, reference 9; Naik 2013 eq 14). See vignette Errata for the units interpretation. Source column TBIL renamed to canonical TBILI.",
      source_name        = "TBIL"
    ),
    ALP = list(
      description        = "Alkaline phosphatase at baseline",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on the Michaelis-Menten Km (reference 87 U/L; Naik 2013 eq 15).",
      source_name        = "ALP"
    ),
    CREAT = list(
      description        = "Serum creatinine at baseline",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear-deviation covariate on Ka, applied only in non-dialysis subjects via the (1 - HEMODIAL) gate (reference 3.3 mg/dL; Naik 2013 eq 13). Source column CR renamed to canonical CREAT.",
      source_name        = "CR"
    ),
    HEMODIAL = list(
      description        = "Hemodialysis-treatment-status indicator at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on hemodialysis)",
      notes              = "1 = subject on hemodialysis (89.6% of PK population per Table 2). Used to gate the CREAT effect on Ka: the CREAT slope applies only when HEMODIAL = 0, matching the paper PDIA indicator (PDIA = 1 - HEMODIAL). Source column PDIA renamed to canonical HEMODIAL with value inversion (HEMODIAL = 1 - PDIA).",
      source_name        = "PDIA"
    ),
    RACE_HISPANIC = list(
      description        = "Hispanic / Latino ethnicity indicator at baseline",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Hispanic)",
      notes              = "1 = Hispanic (14.7% of PK population per Table 2). Used as an additive shift on Ka (Naik 2013 eq 13): the paper-reported ETHN column codes 1 = non-Hispanic, so RACE_HISPANIC = 1 - ETHN and the paper term 0.00811 * (1 - ETHN) maps to 0.00811 * RACE_HISPANIC.",
      source_name        = "ETHN"
    ),
    ESAD = list(
      description        = "Previous erythropoiesis-stimulating-agent (ESA) dose at baseline",
      units              = "units/week (epoetin-equivalent activity per week)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Recorded prior ESA dose for hemodialysis subjects switching to peginesatide. Reference 7996 units/week (paper-reported median; Naik 2013 eq 16). When ESAD = 0 the ESADF indicator gate disables the covariate effect on HgbBL (matching paper text 'no effect of ESAD was incorporated for subjects whose ESAD dose information was not available'). Paper-reported PK-PD population mean 11030 units/week (SD 11866, range 954-105910). New entry proposed alongside the Naik 2013 peginesatide extraction; first-of-class covariate for ESA-switch popPK/PD models.",
      source_name        = "ESAD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 672L,
    n_studies      = 5L,
    age_range      = "21.0-93.0 years (PK population, Table 3)",
    age_median     = "58.4 years (PK population mean, SD 14.6)",
    weight_range   = "38.0-187.5 kg (PK population, Table 3)",
    weight_median  = "79.4 kg (PK population mean, SD 21.6)",
    sex_female_pct = 39.3,
    race_ethnicity = c(White = 57.1, Black = 37.2, Asian = 3.6, Other = 2.1, Hispanic = 14.7),
    disease_state  = "Chronic kidney disease (CKD) on or not on hemodialysis. PK cohort N = 672 (89.6% on dialysis, 10.4% non-dialysis); PK-PD cohort N = 517 hemodialysis subjects only. 2,665 plasma peginesatide concentrations and 18,857 hemoglobin observations.",
    dose_range     = "Peginesatide 0.025-0.16 mg/kg or fixed doses 3-16 mg, IV or subcutaneous, every 4 weeks (some cohorts every 2 weeks)",
    regions        = "Not stated explicitly; Takeda-sponsored phase 2 and phase 3 trials",
    notes          = "Trials AFX01-02, AFX01-03 (NCT00228449), AFX01-04 (NCT00228436), AFX01-07, and AFX01-14 (NCT00597584 EMERALD 1). 63 subjects had hemoglobin samples excluded around transfusions, phlebotomy, GI bleeding, or trauma. Peginesatide was approved as OMONTYS for treatment of CKD-associated anemia in adults on dialysis but was withdrawn from the market by the manufacturer prior to publication of this analysis."
  )

  ini({
    # === PK structural parameters (Naik 2013 Table 6, reference covariates) ===
    # All volumes and inter-compartmental clearance are per kg body weight
    # in the source paper; conversion to absolute units happens in model().
    lka      <- log(0.00865)     ; label("Absorption rate constant base (Ka, 1/hr) for non-Hispanic dialysis subject")  # Naik 2013 Table 6
    lfdepot  <- log(0.498)       ; label("Subcutaneous bioavailability (F1, fraction)")                                  # Naik 2013 Table 6
    lvc      <- log(35.6)        ; label("Central volume of distribution (V2, mL/kg) at reference BMI 26, age 59, TBILI 9")  # Naik 2013 Table 6
    lvp      <- log(7.44)        ; label("Peripheral volume of distribution (V3, mL/kg)")                                 # Naik 2013 Table 6
    lq       <- fixed(log(5.23)) ; label("Inter-compartmental clearance (Q, mL/kg/hr) -- fixed for model stability")      # Naik 2013 Table 6 (Fixed)
    lvmax    <- log(45.3)        ; label("Maximum rate of elimination from central (Vmax, ng/mL/hr)")                     # Naik 2013 Table 6
    lkm      <- log(1880)        ; label("Concentration giving 50% of Vmax (Km, ng/mL) at reference ALP 87")              # Naik 2013 Table 6

    # PK covariate effects (Naik 2013 eq 13, 14, 15)
    e_creat_ka         <-  0.000784 ; label("Linear slope of (CREAT-3.3) on Ka in non-dialysis subjects ((1/hr) per mg/dL)")  # Naik 2013 Table 6 / eq 13
    e_race_hispanic_ka <-  0.00811  ; label("Additive shift on Ka for Hispanic vs non-Hispanic (1/hr)")                       # Naik 2013 Table 6 / eq 13
    e_bmi_vc           <- -0.491    ; label("Power exponent of (BMI/26) on Vc (unitless)")                                    # Naik 2013 Table 6 / eq 14
    e_age_vc           <- -0.125    ; label("Linear slope of (AGE-59) on Vc (mL/kg per year)")                                # Naik 2013 Table 6 / eq 14 (paper-reported units 'L/yr' interpreted as mL/kg/yr; see vignette Errata)
    e_tbili_vc         <-  0.477    ; label("Linear slope of (TBILI-9) on Vc (mL/kg per umol/L)")                              # Naik 2013 Table 6 / eq 14 (paper-reported units 'L/(g/L)' interpreted as mL/kg/(umol/L); see vignette Errata)
    e_alp_km           <- -0.194    ; label("Power exponent of (ALP/87) on Km (unitless)")                                    # Naik 2013 Table 6 / eq 15

    # === PK IIV (Naik 2013 Table 6) ===
    # NONMEM v^2 reported directly as omega^2 (variance on log-transformed parameters).
    # Ka and Vc are correlated: 'Cov between V2 and Ka = -0.0928'.
    etalka + etalvc ~ c(0.197,
                       -0.0928, 0.101)   # var_lka, cov(lka, lvc), var_lvc
    etalkm          ~ 0.0589              # diagonal

    # === PK residual error (Naik 2013 Table 6) ===
    # 'Ratio of proportional to additive residual variability' interpreted as the
    # proportional-component variance sigma2_prop = 0.0218 (CV ~ 14.8%);
    # 's2 (additive component)' = 81.8 (ng/mL)^2 so additive SD ~ 9.04 ng/mL.
    propSd <- 0.1476 ; label("Proportional residual error on peginesatide Cc (fraction)")  # Naik 2013 Table 6 (sqrt of 0.0218)
    addSd  <- 9.0443 ; label("Additive residual error on peginesatide Cc (ng/mL)")          # Naik 2013 Table 6 (sqrt of 81.8)

    # === PD structural parameters (Naik 2013 Table 7) ===
    lhgbbl <- log(11.5)    ; label("Baseline hemoglobin (HgbBL, g/dL) for ESA-naive subject (ESAD = 0)")          # Naik 2013 Table 7
    lec50  <- log(401)     ; label("Concentration for 50% of maximum stimulatory effect (EC50, ng/mL)")           # Naik 2013 Table 7
    lemax  <- log(0.542)   ; label("Maximum stimulatory effect of peginesatide on K0 (Emax, unitless)")           # Naik 2013 Table 7
    lmtt   <- log(1640)    ; label("Mean red-blood-cell lifespan (MTT, hours ~ 68.3 days)")                       # Naik 2013 Table 7
    lmtp   <- log(462)     ; label("Mean transit time for progenitor cells (MTP, hours ~ 19.3 days)")             # Naik 2013 Table 7
    lrsa   <- log(0.153)   ; label("Residual prior-ESA effective concentration (RSA, ng/mL; paper symbol RSA)")   # Naik 2013 Table 7
    lcf    <- log(2.75e-4) ; label("Correction factor for empirical Hgb drift (CF, 1/hr) at reference age 58")    # Naik 2013 Table 7

    # PD covariate effects (Naik 2013 eq 16, 17)
    e_esad_lhgbbl <- -4.49e-7 ; label("Log-scale slope of (ESAD-7996) on HgbBL ((1/(units/week)))")  # Naik 2013 Table 7 / eq 16
    e_age_lcf     <- -0.00314 ; label("Log-scale slope of (AGE-58) on CF (1/year)")                  # Naik 2013 Table 7 / eq 17

    # === PD IIV (Naik 2013 Table 7) ===
    # Diagonal Omega (no cross-correlations reported); v^2 values are variances on
    # log-transformed parameters. EC50 and CF have very high IIV (CV 298.7% and 325.6%).
    etalrsa   ~ 0.0130
    etalec50  ~ 8.92      # IIV 298.7% CV
    etalhgbbl ~ 0.00485
    etalcf    ~ 10.6      # IIV 325.6% CV

    # === PD residual error (Naik 2013 Table 7) ===
    # 's2 (additive component) = 0.00478', so additive SD ~ 0.0691 g/dL
    # (Table 7 footer 'Residual Variability 0.07 SD' matches).
    addSd_Hgb <- 0.0691 ; label("Additive residual error on hemoglobin (g/dL)")  # Naik 2013 Table 7 (sqrt of 0.00478)
  })

  model({
    # =================================================================
    # PK section
    # =================================================================

    # Individual Ka (Naik 2013 eq 13): linear / additive covariates on the bare
    # parameter, then multiplicative IIV on the typical value. The paper's PDIA
    # indicator equals 1 for non-dialysis subjects, so the CREAT slope is gated
    # by (1 - HEMODIAL). The paper's ETHN indicator equals 1 for non-Hispanic;
    # the (1 - ETHN) shift maps to + e_race_hispanic_ka * RACE_HISPANIC.
    ka_tv <- exp(lka) +
             e_creat_ka * (CREAT - 3.3) * (1 - HEMODIAL) +
             e_race_hispanic_ka * RACE_HISPANIC
    ka    <- ka_tv * exp(etalka)

    # Individual Vc (Naik 2013 eq 14): power on BMI, linear-deviation additions
    # on AGE and TBILI. Result in mL/kg, then converted to L via WT scaling.
    v2_per_kg <- exp(lvc) * (BMI / 26)^e_bmi_vc +
                 e_age_vc   * (AGE   - 59) +
                 e_tbili_vc * (TBILI - 9)
    v2 <- v2_per_kg * exp(etalvc)        # mL/kg
    vc <- v2 * WT / 1000                 # L

    vp <- exp(lvp) * WT / 1000           # L (V3 has no covariates)
    q  <- exp(lq)  * WT / 1000           # L/hr (Q fixed in mL/kg/hr)

    # Individual Km and Vmax (Naik 2013 eq 15; Vmax has no covariate or IIV).
    km   <- exp(lkm + etalkm) * (ALP / 87)^e_alp_km   # ng/mL
    vmax <- exp(lvmax)                                 # ng/mL/hr = ug/L/hr

    # Peginesatide concentration: dose in ug, central in ug, vc in L
    # so Cc = central / vc has units ug/L = ng/mL (no further scaling).
    Cc <- central / vc

    # Two-compartment ODE with Michaelis-Menten elimination from central.
    # Vmax * Cc / (Km + Cc) is a concentration rate (ng/mL/hr = ug/L/hr);
    # multiplying by vc (L) converts to a mass rate (ug/hr) for d/dt(central).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          vmax * Cc / (km + Cc) * vc -
                          q * (central / vc - peripheral1 / vp)
    d/dt(peripheral1) <-  q * (central / vc - peripheral1 / vp)

    # Subcutaneous bioavailability F1 applied to depot dosing; IV doses entered
    # directly into central use the default f(central) = 1.
    f(depot) <- exp(lfdepot)

    # =================================================================
    # PD section: precursor-dependent LIDR (Naik 2013 eq 1-7)
    # =================================================================
    # Compartment map (uses the canonical 'precursor<n>' chain prefix):
    #   precursor1            -- bone-marrow progenitor cell pool (paper PRC, A(1))
    #   precursor2..precursor8 -- 7 mature red-blood-cell aging compartments
    #                             (paper Hgb1..Hgb7, A(2)..A(NRBC+1)) discretising the
    #                             RBC lifespan into NRBC = 7 first-order transit stages
    # Hemoglobin is the sum of the 7 RBC stages (Naik 2013 eq 5).

    # Individual PD parameters with covariate effects.
    # HgbBL (Naik 2013 eq 16): exponential covariate on ESAD, gated by ESADF.
    esadf <- (ESAD > 0)
    hgbbl <- exp(lhgbbl + etalhgbbl + e_esad_lhgbbl * (ESAD - 7996) * esadf)

    # CF (Naik 2013 eq 17): exponential covariate on (AGE - 58).
    cf <- exp(lcf + etalcf + e_age_lcf * (AGE - 58))

    # Other PD parameters (Naik 2013 Table 7; no covariate effects in the final model).
    ec50 <- exp(lec50 + etalec50)
    emax <- exp(lemax)
    mtt  <- exp(lmtt)
    mtp  <- exp(lmtp)
    rsa  <- exp(lrsa + etalrsa)

    # LIDR rate constants. K1 = 1/MTP (progenitor first-order transit); KT = NRBC/MTT
    # with NRBC = 7 aging stages so each stage has mean residence MTT/7. K0 is
    # calibrated from the baseline hemoglobin (Naik 2013 eq 8): the (1 + Emax * RSA /
    # (EC50 + RSA)) factor captures the residual prior-ESA stimulation at baseline.
    k1      <- 1 / mtp
    kt      <- 7 / mtt
    k0_stim <- 1 + emax * rsa / (ec50 + rsa)
    k0      <- hgbbl / mtt * k0_stim                # g/dL/hr

    # Peginesatide stimulation of progenitor production (paper eq 1 STM term).
    stm <- 1 + emax * Cc / (ec50 + Cc)

    # Empirical trial-period downward-drift factor on the progenitor -> RBC
    # transit. Naik 2013 Methods describes INT as the "exponential function to
    # empirically account for the downward shift in hemoglobin levels during
    # trial"; implemented here as exp(-CF * t). See vignette Assumptions and
    # deviations for the OCR-degraded equation 1 interpretation.
    int_t <- exp(-cf * t)

    # ODEs.
    d/dt(precursor1) <- k0 * stm - k1 * int_t * precursor1
    d/dt(precursor2) <- k1 * int_t * precursor1 - kt * precursor2
    d/dt(precursor3) <- kt * (precursor2 - precursor3)
    d/dt(precursor4) <- kt * (precursor3 - precursor4)
    d/dt(precursor5) <- kt * (precursor4 - precursor5)
    d/dt(precursor6) <- kt * (precursor5 - precursor6)
    d/dt(precursor7) <- kt * (precursor6 - precursor7)
    d/dt(precursor8) <- kt * (precursor7 - precursor8)

    # Steady-state initial conditions consistent with K0 calibration (Naik 2013
    # eq 6-7). With INT(0) = 1 and Cc(0) = 0 the production-clearance balance
    # gives PRC(0) = K0 * MTP and RBC_j(0) = HgbBL / 7 for j = 1..7.
    precursor1(0) <- k0 * mtp
    precursor2(0) <- hgbbl / 7
    precursor3(0) <- hgbbl / 7
    precursor4(0) <- hgbbl / 7
    precursor5(0) <- hgbbl / 7
    precursor6(0) <- hgbbl / 7
    precursor7(0) <- hgbbl / 7
    precursor8(0) <- hgbbl / 7

    # Hemoglobin observation = sum of 7 RBC aging compartments (paper eq 5).
    Hgb <- precursor2 + precursor3 + precursor4 + precursor5 +
           precursor6 + precursor7 + precursor8

    # Combined additive + proportional error on peginesatide Cc;
    # additive-only error on Hgb (paper-equivalent proportional on the
    # untransformed observation).
    Cc  ~ add(addSd) + prop(propSd)
    Hgb ~ add(addSd_Hgb)
  })
}
