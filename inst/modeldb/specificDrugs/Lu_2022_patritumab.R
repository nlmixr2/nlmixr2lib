Lu_2022_patritumab <- function() {
  description <- "Joint two-analyte population PK model for patritumab deruxtecan (HER3-DXd, an anti-HER3 antibody-drug conjugate) in adults with HER3-expressing solid tumors (Lu 2022). DXd-conjugated antibody (intact ADC) is described by a 2-compartment model with parallel linear and Michaelis-Menten clearance. Released unconjugated DXd (MAAA-1181a, exatecan-derivative payload) is described by a 1-compartment model with linear clearance and a first-order, time-dependent release rate driven by the level of DXd-conjugated antibody in the central compartment, scaled by the molecular-weight ratio MW_DXd/MW_DXdAb and a payload-to-intact-drug ratio PIR modulated by a cycle-1-vs-later (factor1) and a within-cycle exponential (factor2) modifier."
  reference <- "Lu Y, Tang N, Hayashi L, Lin Z, Sasaki R, Asakawa M, Liu X, Sahasranaman S, Yamamoto N, Nakagawa K, Janne PA, Schmid P. Population Pharmacokinetics of Patritumab Deruxtecan in Patients With Solid Tumors. J Clin Pharmacol. 2023;63(1):77-88. doi:10.1002/jcph.2137. PMID 36053771."
  vignette <- "Lu_2022_patritumab"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Power effects on CLlin (exponent 0.911) and Vc (exponent 0.654) of DXd-conjugated antibody, and on Krel (exponent -0.467) of unconjugated DXd. Reference 60 kg per Lu 2022 simulation reference patient ('a male with weight of 60 kg, albumin level of 39 g/L, NSCLC, and normal hepatic function', Methods: Simulations of Covariate Effects); the Lu 2022 Table S2 covariate-function header states the reference is the population median, which is 58.6 kg per Table S4 distribution of continuous baseline factors. The 60 kg value is preferred here because it matches the explicit simulation reference in the paper text and is used as the baseline for the published forest-plot exposures.",
      source_name        = "Weight"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value (SI units; not g/dL). Power effect on CLlin of DXd-conjugated antibody (exponent -0.795); reference 39 g/L per Lu 2022 simulation reference patient (Methods: Simulations of Covariate Effects) and population median (Table S4).",
      source_name        = "Albumin level"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) per the canonical SEXF column. Lu 2022 also uses male as the reference category and reports the female effect.",
      notes              = "Time-fixed. Multiplicative fractional effect e_female_cl on CLlin of DXd-conjugated antibody applied as (1 + e_female_cl * SEXF) so SEXF = 1 yields the paper-reported 0.865 multiplier (e_female_cl = -0.135) and SEXF = 0 yields multiplier 1 (male reference, unchanged). Lu 2022 Table 2 reports 'Female-CLlin = 0.865' (relative-to-reference female multiplier; Table S2 row Sex states the reference is male and the covariate function is PV = PVref * theta^i for indicator i).",
      source_name        = "Sex"
    ),
    TUMTP_BC = list(
      description        = "Breast-cancer tumor-type indicator, 1 = breast cancer, 0 = NSCLC or colorectal cancer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = NSCLC or colorectal cancer pooled. NSCLC is the implicit reference category in Lu 2022; colorectal cancer was tested and found insignificant relative to NSCLC and is therefore pooled into the reference (Lu 2022 Results: 'the colorectal cancer effect was insignificant relative to the reference NSCLC').",
      notes              = "Time-fixed. Multiplicative fractional effect e_bc_cl on CLlin of DXd-conjugated antibody applied as (1 + e_bc_cl * TUMTP_BC); TUMTP_BC = 1 yields the paper-reported 0.811 multiplier (e_bc_cl = -0.189) and TUMTP_BC = 0 yields multiplier 1. Lu 2022 Table 2 reports 'Breast cancer-CLlin = 0.811'. CRC was retained at the covariate-evaluation stage and dropped on backward elimination.",
      source_name        = "Tumor type (breast cancer level)"
    ),
    HEPIMP_MILD = list(
      description        = "Mild hepatic impairment indicator (NCI ODWG group 2)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, NCI ODWG group 1; mutually exclusive with HEPIMP_MOD_OR_MISSING).",
      notes              = "Time-fixed. Multiplicative fractional effect e_hepmild_cl_dxd on CLDXd of unconjugated DXd applied as (1 + e_hepmild_cl_dxd * HEPIMP_MILD); HEPIMP_MILD = 1 yields the paper-reported 0.706 multiplier (e_hepmild_cl_dxd = -0.294). Lu 2022 Table 3 reports 'Hepatic impairment mild-CLDXd = 0.706'. Categorization per the National Cancer Institute Organ Dysfunction Working Group: mild = total bilirubin <= ULN with AST > ULN, OR total bilirubin > 1.0 to 1.5 x ULN with any AST.",
      source_name        = "Hepatic function (mild impairment)"
    ),
    HEPIMP_MOD_OR_MISSING = list(
      description        = "Composite moderate-hepatic-impairment-or-data-missing indicator (paper-specific)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function or mild impairment; mutually exclusive with HEPIMP_MILD when both indicators are zero, the patient is in the normal-function reference).",
      notes              = "Time-fixed. Multiplicative fractional effect e_hepmod_cl_dxd on CLDXd of unconjugated DXd applied as (1 + e_hepmod_cl_dxd * HEPIMP_MOD_OR_MISSING); HEPIMP_MOD_OR_MISSING = 1 yields the paper-reported 0.532 multiplier (e_hepmod_cl_dxd = -0.468). Lu 2022 Table 3 reports 'Hepatic impairment moderate/data missing-CLDXd = 0.532' as a single pooled coefficient because the moderate-impairment subgroup (n = 6 per Lu 2022 Table S5) and the missing/unknown subgroup (n = 6 per Table S5) were both individually too small to estimate as separate effects. Lu 2022 Discussion notes that 'the effect of moderate hepatic function was not estimated, as the relevant data were too limited (n = 6) for the evaluation to be reliable using the population PK approach' - the reported 0.532 effect is therefore not a moderate-only estimate but a moderate-plus-missing composite.",
      source_name        = "Hepatic function (moderate impairment OR data missing)"
    ),
    CYCLE = list(
      description        = "Treatment cycle number, 1 = first dosing cycle, 2 = second cycle, ...",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = "Required for the unconjugated-DXd submodel only. CYCLE = 1 in the cycle-1 reference; cycle 2 and later receive a multiplicative scaling of the payload-to-intact-drug ratio PIR by factor1 = theta = 0.648 (Lu 2022 Eqs. 5 and 7). Time-varying across a multi-cycle treatment course; supply CYCLE >= 1 at every observation row in the data set, incrementing at the start of each new dosing cycle. Does not affect DXd-conjugated antibody PK.",
      source_name        = "Cycle (derived from time and dosing schedule)"
    )
  )

  population <- list(
    n_subjects     = 425L,
    n_studies      = 3L,
    n_observations_dxdab = 6986L,
    n_observations_dxd   = 7152L,
    age_range      = "29-83 years (median 60)",
    age_median     = "60 years",
    weight_range   = "32.4-114.1 kg (median 58.6)",
    weight_median  = "58.6 kg (population median per Lu 2022 Table S4); the simulation reference patient is rounded to 60 kg",
    sex_female_pct = 75.3,
    race_ethnicity = c(White = 35.3, Asian = 57.6, Black = 2.8, Other_Unknown = 4.2),
    disease_state  = "HER3-expressing advanced or metastatic solid tumors. 50.8% (216/425) NSCLC, 39.8% (169/425) breast cancer (HR+/HER2- and triple-negative), 9.4% (40/425) colorectal cancer.",
    dose_range     = "1.6 to 8.0 mg/kg IV every 3 weeks (Q3W) for the majority of patients; a small cohort (n = 12) received 4.2 mg/kg every 2 weeks (Q2W) for cycles 1-3 with a switch to 6.4 mg/kg Q3W from cycle 4 onwards (J101 cohort 7). 5.6 mg/kg Q3W was the dosage selected for further development (55% of patients in the analysis dataset).",
    regions        = "Multi-regional. Three studies: J101 (NCT02980341, dose escalation in Japan with global dose expansion in HER3-positive metastatic breast cancer), U102 (NCT03260491, NSCLC, US/EU/Asia), U202 (NCT04479436, advanced/metastatic colorectal cancer, global).",
    hepatic_function = c(Normal = 71.8, Mild_impairment = 25.4, Moderate_impairment = 1.4, Missing_unknown = 1.4),
    renal_function   = c(Normal = 37.4, Mild_impairment = 41.6, Moderate_impairment = 20.7, Severe_impairment = 0.2),
    reference_subject = "Male, weight 60 kg, albumin 39 g/L, NSCLC, normal hepatic function (cycle 1 for the unconjugated-DXd submodel) per Lu 2022 Methods: Simulations of Covariate Effects.",
    notes          = "NONMEM 7.3 with FOCEI on the Metworx platform. PK observations beyond day 200 (10 + 12 trough samples for DXd-conjugated antibody and unconjugated DXd respectively) were excluded because they shortened the run time by >= 40% without affecting modelling results. CrCL, ECOG, race, race-country, baseline tumor SOD, and drug formulation were tested as candidate covariates and not retained in the final models. Anti-drug antibody effect was not evaluated because immunogenicity was rare (2.3% positive). The two analytes were fit sequentially: DXd-conjugated antibody first, then unconjugated DXd with the antibody parameters fixed to their final-model estimates."
  )

  ini({
    # ============================================================
    # DXd-conjugated antibody (intact ADC) -- Lu 2022 Table 2
    # ============================================================
    # Reference subject for typical-value parameters: male (SEXF = 0),
    # weight 60 kg, albumin 39 g/L, NSCLC (TUMTP_BC = 0). Time unit kept
    # as day to match the paper's L/d / ug/d tabulation.
    lcl   <- log(0.342);   label("Linear clearance of DXd-conjugated antibody CLlin at reference covariates (L/day)")  # Lu 2022 Table 2: CLlin = 0.342 L/d
    lvc   <- log(2.943);   label("Central volume of distribution of DXd-conjugated antibody Vc at reference (L)")     # Lu 2022 Table 2: Vc = 2.943 L
    lvmax <- log(4.3522);  label("Maximum nonlinear clearance Vmax of DXd-conjugated antibody (mg/day)")              # Lu 2022 Table 2: Vmax = 4352.2 ug/d = 4.3522 mg/d
    lkm   <- log(0.5759);  label("Michaelis-Menten constant Km of DXd-conjugated antibody (mg/L; equivalent to 575.9 ng/mL)") # Lu 2022 Table 2: Km = 575.9 ng/mL = 0.5759 mg/L
    lq    <- log(0.443);   label("Intercompartmental clearance Q of DXd-conjugated antibody (L/day)")                 # Lu 2022 Table 2: Q = 0.443 L/d
    lvp   <- log(5.369);   label("Peripheral volume of distribution Vp of DXd-conjugated antibody (L)")               # Lu 2022 Table 2: Vp = 5.369 L

    # Covariate effects on DXd-conjugated antibody (Lu 2022 Table 2 and
    # Table S2 covariate-function header). Continuous covariates use
    # power-of-ratio scaling (PV = PVref * (cov / cov_ref)^theta);
    # categorical covariates use fractional multipliers (1 + theta * indicator)
    # so that e_female_cl = -0.135 reproduces the paper-reported 0.865
    # multiplier when SEXF = 1, and e_bc_cl = -0.189 reproduces the
    # paper-reported 0.811 multiplier when TUMTP_BC = 1.
    e_wt_cl     <-  0.911;  label("Power exponent of WT on CLlin (unitless)")                                       # Lu 2022 Table 2: Weight-CLlin = 0.911
    e_alb_cl    <- -0.795;  label("Power exponent of ALB on CLlin (unitless)")                                      # Lu 2022 Table 2: Albumin-CLlin = -0.795
    e_female_cl <- -0.135;  label("Fractional multiplicative effect of female sex on CLlin (unitless)")             # Lu 2022 Table 2: Female-CLlin = 0.865 -> multiplier (1 - 0.135) = 0.865
    e_bc_cl     <- -0.189;  label("Fractional multiplicative effect of breast-cancer tumor type on CLlin (unitless)") # Lu 2022 Table 2: Breast cancer-CLlin = 0.811 -> multiplier (1 - 0.189) = 0.811
    e_wt_vc     <-  0.654;  label("Power exponent of WT on Vc (unitless)")                                          # Lu 2022 Table 2: Weight-Vc = 0.654

    # IIV on DXd-conjugated antibody. Lu 2022 Table 2 'Variance of CLlin
    # IIV', 'Variance of Vc IIV', and 'Covariance of Vc and CLlin'
    # report omega^2 / cov directly on the log scale (no CV<->omega^2
    # conversion needed).
    etalcl + etalvc ~ c(0.215,
                        0.030, 0.023)  # Lu 2022 Table 2: Var(CLlin) = 0.215, Cov(Vc, CLlin) = 0.030, Var(Vc) = 0.023

    # ============================================================
    # Unconjugated DXd (MAAA-1181a, payload) -- Lu 2022 Table 3
    # ============================================================
    # Time units in Table 3 are HOURS (Krel in 1/h, CLDXd in L/h, beta in
    # 1/h). Converted to per-day for consistency with the antibody
    # parameters (multiply rates / clearances by 24).
    lkrel    <- log(0.030 * 24);   label("First-order DXd release-rate constant Krel from intact ADC (1/day; converted from 0.030 1/h)")  # Lu 2022 Table 3: Krel = 0.030 1/h
    lcl_dxd   <- log(5.713 * 24);   label("Linear clearance of unconjugated DXd CLDXd (L/day; converted from 5.713 L/h)")                  # Lu 2022 Table 3: CLDXd = 5.713 L/h
    ltheta   <- log(0.648);        label("Factor 1 multiplier on PIR for cycles >= 2 (unitless)")                                          # Lu 2022 Table 3: theta = 0.648
    lbeta    <- log(0.180 * 24);   label("Exponential decline rate constant beta of factor 2 within each dosing cycle (1/day; converted from 0.180 1/h)") # Lu 2022 Table 3: beta = 0.180 1/h

    # Covariate effects on unconjugated DXd (Lu 2022 Table 3).
    e_wt_krel         <- -0.467;  label("Power exponent of WT on Krel (unitless)")                                                                          # Lu 2022 Table 3: Weight-Krel = -0.467
    e_hepmild_cl_dxd  <- -0.294;  label("Fractional multiplicative effect of mild hepatic impairment on CLDXd (unitless)")                                  # Lu 2022 Table 3: Hep mild-CLDXd = 0.706 -> multiplier (1 - 0.294) = 0.706
    e_hepmod_cl_dxd   <- -0.468;  label("Fractional multiplicative effect of moderate-or-data-missing hepatic impairment on CLDXd (unitless)")              # Lu 2022 Table 3: Hep moderate/missing-CLDXd = 0.532 -> multiplier (1 - 0.468) = 0.532

    # IIV on unconjugated DXd. Lu 2022 Table 3 reports omega^2 / cov on
    # the log scale directly (Var(Krel IIV) = 0.194, Var(CLDXd IIV) = 0.360,
    # Cov(Krel, CLDXd) = 0.142).
    etalkrel + etalcl_dxd ~ c(0.194,
                             0.142, 0.360)   # Lu 2022 Table 3: Var(Krel)=0.194, Cov(Krel,CLDXd)=0.142, Var(CLDXd)=0.360

    # Residual error for both analytes. Lu 2022 Table 2 (DXd-conjugated
    # antibody): proportional error 0.236 (fraction), additive error fixed
    # to 100 ng/mL = 0.1 ug/mL (the LLOQ for the DXd-conjugated antibody
    # ELISA). Lu 2022 Table 3 (unconjugated DXd): proportional error 0.392
    # (fraction), additive error fixed to 0.01 ng/mL (the LLOQ for the
    # unconjugated-DXd LC-MS assay).
    addSd     <- fix(0.1);   label("Additive residual error on DXd-conjugated antibody Cc (ug/mL; equivalent to 100 ng/mL = LLOQ, fixed)") # Lu 2022 Table 2: Additive error 100 ng/mL Fixed
    propSd    <- 0.236;      label("Proportional residual error on DXd-conjugated antibody Cc (fraction)")                                  # Lu 2022 Table 2: Proportional error 0.236
    addSd_dxd     <- fix(0.01);  label("Additive residual error on unconjugated DXd Cc_dxd (ng/mL; equivalent to LLOQ for DXd, fixed)")         # Lu 2022 Table 3: Additive error 0.01 ng/mL Fixed
    propSd_dxd    <- 0.392;      label("Proportional residual error on unconjugated DXd Cc_dxd (fraction)")                                     # Lu 2022 Table 3: Proportional error 0.392
  })

  model({
    # ============================================================
    # Constants for the DXd release flux (Lu 2022 Eqs. 3, 5, 6, 7)
    # ============================================================
    # Molecular weights from Lu 2022 Methods text immediately following
    # Eq. 3: 'MWDXd 493.5 g/mol' and 'MWDXdAb approx. 150 000 g/mol'.
    # Payload-to-intact-drug ratio PIR0 = 8 (Lu 2022 Methods after Eq. 3
    # and Discussion: 'a drug-to-antibody ratio of approximately 8').
    # Factor-2 floor alpha = 0.125 (Lu 2022 Eq. 6 footnote: 'alpha = 0.125
    # ... resulting in a modified PIR between 1 and 8').
    mw_dxd     <- 493.5
    mw_dxdab   <- 150000
    pir0       <- 8
    alpha_pir  <- 0.125

    # ============================================================
    # Individual PK parameters for DXd-conjugated antibody
    # ============================================================
    # Reference covariates: WT 60 kg, ALB 39 g/L, SEXF = 0 (male),
    # TUMTP_BC = 0 (NSCLC reference). Continuous covariates enter as
    # power-of-ratio; categorical covariates enter as fractional
    # multipliers so SEXF = 0 and TUMTP_BC = 0 yield multiplier 1.
    cl <- exp(lcl + etalcl) *
      (WT  / 60)^e_wt_cl *
      (ALB / 39)^e_alb_cl *
      (1 + e_female_cl * SEXF) *
      (1 + e_bc_cl     * TUMTP_BC)

    vc <- exp(lvc + etalvc) * (WT / 60)^e_wt_vc

    vmax <- exp(lvmax)
    km   <- exp(lkm)
    q    <- exp(lq)
    vp   <- exp(lvp)

    # ============================================================
    # Individual PK parameters for unconjugated DXd
    # ============================================================
    krel  <- exp(lkrel  + etalkrel) * (WT / 60)^e_wt_krel
    cl_dxd <- exp(lcl_dxd + etalcl_dxd) *
      (1 + e_hepmild_cl_dxd * HEPIMP_MILD) *
      (1 + e_hepmod_cl_dxd  * HEPIMP_MOD_OR_MISSING)
    theta_factor1 <- exp(ltheta)
    beta_pir      <- exp(lbeta)

    # Apparent volume of distribution for unconjugated DXd. Lu 2022
    # Figure 1 caption identifies VDXd as a model parameter but neither
    # Table 3 (final), Table S7 (base), the main text, nor the supplement
    # report a numerical value for it. Fixed here to 1 L so that
    # CLDXd / V_DXd has time-units 1/day and the unconjugated-DXd ODE
    # mass balance is well-defined; the resulting DXd concentrations are
    # determined predominantly by the formation flux (release-rate-limited
    # regime) so the choice of VDXd has only a transient (~1/CLDXd) effect
    # on the simulated profile shape. Documented as a deviation in the
    # validation vignette's 'Errata' and 'Assumptions and deviations'
    # sections.
    vdxd <- 1.0

    # ============================================================
    # Time-dependent PIR modulation (Lu 2022 Eqs. 5-7)
    # ============================================================
    # Factor 1 differentiates cycle 1 (factor1 = 1) from later cycles
    # (factor1 = theta = 0.648). CYCLE is supplied per-observation by
    # the user (Lu 2022 Eq. 5).
    factor1 <- 1 + (theta_factor1 - 1) * (CYCLE > 1)
    # Factor 2 is an exponential decay of PIR within each dosing cycle,
    # bounded between alpha (= 0.125) and 1 (Lu 2022 Eq. 6). tad() returns
    # time after the most recent dose (days); beta_pir has units 1/day.
    factor2 <- alpha_pir + (1 - alpha_pir) * exp(-beta_pir * tad())
    pir <- pir0 * factor1 * factor2

    # ============================================================
    # ODE system: 2-cmt DXd-conjugated antibody + 1-cmt DXd
    # ============================================================
    # Antibody: parallel linear (CLlin) and Michaelis-Menten (CLnonlin =
    # Vmax / (Km + C)) clearance from the central compartment, plus
    # 2-compartment distribution (Lu 2022 Eqs. 1-2 and Figure 1).
    # central in mg, vc in L -> C = central / vc in mg/L = ug/mL.
    cdxdab <- central / vc
    cl_nonlin <- vmax / (km + cdxdab)

    d/dt(central)     <- -cl       * cdxdab -
                          cl_nonlin * cdxdab -
                          q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    # Unconjugated DXd: first-order release driven by the AMOUNT of
    # DXd-conjugated antibody in the central compartment, scaled by PIR
    # (with factor 1 and factor 2 modulation) and the MW ratio
    # MW_DXd / MW_DXdAb (Lu 2022 Eq. 3); linear clearance via CLDXd from
    # an apparent V_DXd compartment (V_DXd = 1 L; see deviation note above).
    rrelease <- krel * central * pir * (mw_dxd / mw_dxdab)
    d/dt(central_dxd) <- rrelease - cl_dxd / vdxd * central_dxd

    # ============================================================
    # Observations and error model
    # ============================================================
    # DXd-conjugated antibody concentration in ug/mL (= mg/L).
    Cc <- cdxdab
    # Unconjugated DXd concentration in ng/mL. central_dxd is in mg, vdxd
    # in L -> mg/L = ug/mL; multiply by 1000 to obtain ng/mL (the unit
    # used by Lu 2022 to report unconjugated-DXd concentrations).
    Cc_dxd <- (central_dxd / vdxd) * 1000

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_dxd ~ add(addSd_dxd) + prop(propSd_dxd)
  })
}
