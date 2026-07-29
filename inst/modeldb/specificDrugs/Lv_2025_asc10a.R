Lv_2025_asc10a <- function() {
  description <- paste(
    "Population PK model of ASC10-A (also called NHC,",
    "beta-D-N4-hydroxycytidine), the active nucleoside metabolite of",
    "the oral double prodrug ASC10, in Chinese healthy adult volunteers",
    "(Lv 2025). ASC10 is under clinical development for COVID-19 and",
    "shares an active moiety with molnupiravir. Data pooled 57 subjects",
    "with 1,634 ASC10-A plasma concentrations from a Phase I trial",
    "(NCT05523141) with a multiple-ascending-dose part (50-800 mg BID x",
    "5 days + one AM dose on day 6) and a fasted/fed food-effect part.",
    "The final model is a two-compartment disposition with first-order",
    "elimination and a two-transit-compartment absorption chain (rate",
    "constant KTR). Body weight enters as an estimated power on CL/F",
    "(exponent 0.903) and a fixed power on Vc/F (exponent 1.0),",
    "normalised to the median 61.4 kg. Food status enters KTR only:",
    "fed KTR is scaled by theta_KTR-food = 0.474, i.e., ~52.6%",
    "reduction relative to fasted, with no effect on relative",
    "bioavailability (F is fixed at 1). Residual error is combined",
    "additive (2.3 ng/mL, fixed) + proportional (SD = 0.41).",
    sep = " "
  )
  reference <- paste(
    "Lv D, Li S, Li Y, Lin M, Zhai Y, Wu M, Qiu Y, Zhao Q, Liu J.",
    "Population Pharmacokinetic Modeling Analysis of ASC10, a Novel",
    "Antiviral Agent Targeted COVID-19, in Chinese Healthy Subjects.",
    "Drug Des Devel Ther. 2025;19:7391-7402.",
    "doi:10.2147/DDDT.S517282.",
    sep = " "
  )
  vignette <- "Lv_2025_asc10a"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at enrolment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at enrolment. Enters CL/F through an estimated power",
        "(theta_CL-WT = 0.903) and Vc/F through a fixed unit power",
        "(theta_Vc-WT = 1 FIX), both normalised to the population median",
        "of 61.4 kg (Lv 2025 Results paragraph 1 and Table 2). BMI was",
        "screened as an alternative body-size covariate but was not",
        "retained in the final model."
      ),
      source_name        = "WT"
    ),
    FED = list(
      description        = "Fed vs fasted status at dosing",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 = fasted",
      notes              = paste(
        "1 = fed (postprandial), 0 = fasted. Enters KTR via the",
        "categorical multiplicative form P_i = theta_pop * theta_cov^cov",
        "(Lv 2025 equation 2), so KTR_fed = 7.02 * 0.474 = 3.33 h^-1",
        "and KTR_fasted = 7.02 h^-1. Food effect delays absorption",
        "(reduces peak, prolongs Tmax) without altering relative",
        "bioavailability (F is fixed at 1); AUCss was reported to be",
        "essentially unchanged between fasted and fed."
      ),
      source_name        = "FED"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at enrolment",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a potential covariate but not retained (Lv 2025",
        "Results 'PopPK Model Development'). Range 18-44 y (median 31)."
      )
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 = male",
      notes              = paste(
        "Screened as a potential covariate but not retained (Lv 2025",
        "Results 'PopPK Model Development'). Cohort was 42.1% female."
      )
    ),
    BMI = list(
      description        = "Body mass index at enrolment",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Tested as an alternative body-size covariate for CL/F and",
        "Vc/F but did not improve model fit relative to WT (Lv 2025",
        "Results 'PopPK Model Development'). Range 19-25 kg/m^2."
      )
    ),
    RACE_HAN = list(
      description        = "Han Chinese ethnicity indicator",
      units              = "(binary)",
      type               = "categorical",
      reference_category = "0 = non-Han",
      notes              = paste(
        "Screened as a potential covariate but not retained (Lv 2025",
        "Results 'PopPK Model Development'). Cohort was 96.5% Han."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_studies      = 1L,
    n_observations = 1634L,
    age_range      = "18-44 years",
    age_median     = "31 years",
    weight_range   = "47.1-77.9 kg",
    weight_median  = "61.4 kg (population median; used as normalisation reference for allometric WT covariate)",
    bmi_range      = "19-25 kg/m^2",
    sex_female_pct = 42.1,
    race_ethnicity = c(Han = 96.5, Other = 3.5),
    disease_state  = "Healthy Chinese adult volunteers (no COVID-19 or other disease at enrolment).",
    dose_range     = paste(
      "MAD (n = 45): oral ASC10 tablets 50, 100, 200, 400, 600, or",
      "800 mg twice daily for 5 consecutive days plus a single morning",
      "dose on day 6. FE (n = 12): single 800 mg oral dose in fasted",
      "and fed states with >= 7 day washout (crossover)."
    ),
    regions        = "Mainland China (Hangzhou; single-centre)",
    trial_registration = "ClinicalTrials.gov NCT05523141",
    bql_handling   = paste(
      "22.6% of observations were below the quantification limit",
      "(LLOQ = 10.0 ng/mL for ASC10-A by UPLC-MS/MS). Handled via the",
      "M3 likelihood-based method with BQL data included as categorical",
      "observations; population parameters estimated with the Laplacian",
      "algorithm in NONMEM 7.4 (Lv 2025 Methods 'PopPK Model",
      "Development' and Results 'Demographics and Data Summary')."
    ),
    notes          = paste(
      "Modelled analyte is the active nucleoside metabolite ASC10-A",
      "(equivalent to NHC, beta-D-N4-hydroxycytidine). ASC10 is a",
      "double prodrug that is rapidly biotransformed to molnupiravir",
      "and then to ASC10-A, or transformed directly to ASC10-A. Doses",
      "in the event table are the parent ASC10 doses in mg (F applied",
      "at the depot); the model does not resolve the parent-to-active",
      "conversion mechanistically. Compare with Bihorel 2022 and Gouda",
      "2023 popPK models for molnupiravir/NHC."
    )
  )

  ini({
    # =========================================================================
    # Structural PK: two-compartment disposition with first-order elimination
    # and a two-transit-compartment absorption chain. Population typical
    # values from Lv 2025 Table 2 ('Final Model Estimate' column), all
    # apparent (F-adjusted) since F is fixed at 1.
    # =========================================================================
    lcl     <- log(79.7);  label("Apparent ASC10-A clearance CL/F at the reference 61.4 kg (L/h)")             # Lv 2025 Table 2 row 'CL (L/h) = 79.7 (RSE 2.9%)'
    lvc     <- log(139);   label("Apparent ASC10-A central volume Vc/F at the reference 61.4 kg (L)")          # Lv 2025 Table 2 row 'V_C (L) = 139 (RSE 3.1%)'
    lq      <- log(1.24);  label("Apparent inter-compartmental clearance Q/F (L/h)")                            # Lv 2025 Table 2 row 'Q (L/h) = 1.24 (RSE 9.3%)'
    lvp     <- log(31.6);  label("Apparent ASC10-A peripheral volume Vp/F (L)")                                 # Lv 2025 Table 2 row 'V_P (L) = 31.6 (RSE 11.3%)'
    lktr    <- log(7.02);  label("Transit rate constant KTR for fasted state (1/h)")                            # Lv 2025 Table 2 row 'KTR (1/h) = 7.02 (RSE 19.5%)' (typical value = fasted; food scaling applied below)
    lfdepot <- fixed(log(1)); label("Reference relative oral bioavailability F (fraction; fixed at 100%)")     # Lv 2025 Table 2 row 'F = 1 FIX'

    # =========================================================================
    # Covariate effect parameters. Continuous covariates enter as
    # theta_pop * (cov/cov_median)^theta_cov (Lv 2025 equation 1);
    # categorical covariates as theta_pop * theta_cov^cov (Lv 2025
    # equation 2, with cov in {0, 1}).
    # =========================================================================
    e_wt_cl   <- 0.903;         label("Power exponent for body weight on ASC10-A CL/F (unitless)")              # Lv 2025 Table 2 row 'theta_CL-WT = 0.903 (RSE 21.2%)'; CL_i = 79.7 * (WT/61.4)^0.903
    e_wt_vc   <- fixed(1);      label("Power exponent for body weight on ASC10-A Vc/F (unitless; fixed)")       # Lv 2025 Table 2 row 'theta_Vc-WT = 1 FIX'; Vc_i = 139 * (WT/61.4)^1
    r_ktr_fed <- 0.474;         label("Fed/fasted multiplier on KTR (unitless; KTR_fed = KTR_fasted * this)")   # Lv 2025 Table 2 row 'theta_KTR-food = 0.474 (RSE 16.8%)'; KTR_fed = 7.02 * 0.474 = 3.33 h^-1

    # =========================================================================
    # Between-subject variability (IIV). Reported as omega (square root of
    # interindividual variance for exponential IIV) with %CV interpretation.
    # Converted to log-normal omega^2 via omega^2 = log(1 + CV^2),
    # following the convention used in Bihorel 2022 (cited by this paper)
    # and Hien 2017. No IIV was reported on Q, Vp, or F (only IIV on F is
    # implicit in the fixed F = 1 anchoring).
    # =========================================================================
    etalcl  ~ 0.02884  # Lv 2025 Table 2 IIV on CL = 17.1% CV (RSE 10.6%); omega^2 = log(1 + 0.171^2) = 0.02884
    etalvc  ~ 0.01055  # Lv 2025 Table 2 IIV on Vc = 10.3% CV (RSE 22.2%); omega^2 = log(1 + 0.103^2) = 0.01055
    etalktr ~ 0.09693  # Lv 2025 Table 2 IIV on KTR = 31.9% CV (RSE 9.5%); omega^2 = log(1 + 0.319^2) = 0.09693

    # =========================================================================
    # Residual variability: combined additive + proportional. Proportional
    # SD is reported on the linear plasma concentration scale (fraction);
    # additive SD is in ng/mL and was fixed during estimation.
    # =========================================================================
    propSd <- 0.41;        label("Proportional residual error on plasma ASC10-A concentration (fraction)")      # Lv 2025 Table 2 row 'sigma_prop = 0.41 (RSE 5.2%)'
    addSd  <- fixed(2.3);  label("Additive residual error on plasma ASC10-A concentration (ng/mL; fixed)")      # Lv 2025 Table 2 row 'sigma_add (ng/mL) = 2.3 FIX'
  })

  model({
    # =========================================================================
    # 1. Individual PK parameters. Body weight scales CL/F (estimated
    #    exponent) and Vc/F (fixed unit exponent), both normalised to
    #    the population median of 61.4 kg. Food status scales KTR via
    #    the multiplicative categorical form (Lv 2025 equation 2).
    # =========================================================================
    cl     <- exp(lcl  + etalcl) * (WT / 61.4)^e_wt_cl
    vc     <- exp(lvc  + etalvc) * (WT / 61.4)^e_wt_vc
    q      <- exp(lq)
    vp     <- exp(lvp)
    ktr    <- exp(lktr + etalktr) * r_ktr_fed^FED
    fdepot <- exp(lfdepot)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Plasma ASC10-A concentration (ng/mL) from central amount (mg) and
    # apparent central volume (L): amt(mg) / V(L) = mg/L = 1000 * ng/mL.
    Cc <- 1000 * central / vc

    # =========================================================================
    # 2. Two-transit-compartment absorption chain into two-compartment
    #    disposition. Dose enters depot; depot -> transit1 -> transit2 ->
    #    central, all three transitions at the individual rate ktr
    #    (Lv 2025 Results 'PopPK Model Development': the final model has
    #    two transit compartments and first-order elimination). The
    #    dashed-arrow schematic in Figure 1 corresponds to depot ->
    #    transit1 -> transit2 -> central.
    # =========================================================================
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Bioavailability applied at the dosing compartment.
    f(depot) <- fdepot

    # =========================================================================
    # 3. Observation and residual error (combined additive + proportional).
    # =========================================================================
    Cc ~ add(addSd) + prop(propSd)
  })
}
