Lin_2021_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in adult Chinese ICU patients, with additive-linear covariate equations on CL (dopamine co-medication, serum creatinine, severe-burn status, total body weight, CRRT status) and on V (serum creatinine, age, total body weight) and additive-normal between-subject variability, fitted by the EM algorithm in Kinetica (Lin 2021)"
  reference <- "Lin Z, Chen DY, Zhu YW, Jiang ZL, Cui K, Zhang S, Chen LH. Population pharmacokinetic modeling and clinical application of vancomycin in Chinese patients hospitalized in intensive care units. Sci Rep. 2021;11(1):2670. doi:10.1038/s41598-021-82312-2"
  vignette <- "Lin_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # The between-subject variability is additive-normal on the linear
  # parameter scale (Lin 2021 Methods: 'beta_i = Z_i beta + eta_j', eta
  # normal with mean 0), so the etas act on the linear typical values cl and
  # vc rather than on the log-scale intercepts lcl_int / lvc_int.
  paper_specific_etas <- c("etacl", "etavc")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column TBW. Enters both covariate equations as an additive-linear term with no centring: +0.067 L/h per kg on CL and +0.16 L per kg on V (Lin 2021 Tables 3 and 4). Cohort median 65.0 kg (range 40.0-90.6; Table 1, n = 374).",
      source_name = "TBW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column Age. Enters the V equation only, additive-linear with no centring: -0.07 L per year (Lin 2021 Table 4). Cohort median 62 years (range 18-93; Table 1).",
      source_name = "Age"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Source column Cr, in umol/L (Lin 2021 Methods, creatine-oxidase method on a Beckmann CX5). Enters both covariate equations additive-linear with no centring: -0.007 L/h per umol/L on CL and -0.04 L per umol/L on V (Tables 3 and 4). Cohort median 71.0 umol/L (range 28.0-581.0; Table 1). Because the CL equation is additive-linear it crosses zero at high creatinine combined with low weight (e.g. about 494 umol/L at 40 kg without CRRT, 435 umol/L at 40 kg on CRRT); see the vignette.",
      source_name = "Cr"
    ),
    CONMED_DOPA = list(
      description = "Concomitant dopamine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Source column DA, coded 1 = dopamine not used, 2 = dopamine used (Lin 2021 Methods covariate list). Recoded as CONMED_DOPA = DA - 1 (1 = dopamine, 0 = none); model() rebuilds the source code as CONMED_DOPA + 1 so the printed coefficient +0.036 L/h on DA applies unchanged. Dopamine therefore raises CL by 0.036 L/h, matching the paper's statement that CL is 'accelerated in the presence of rescue drugs (DA)'. The fraction of the cohort on dopamine is not reported.",
      source_name = "DA"
    ),
    DIS_BURN_RECENT = list(
      description = "Severe acute-phase burn indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Source column Burn-S, coded 1 = burn degree > 50%, acute convalescence; 2 = burn degree < 50%, 20 days after burn (Lin 2021 Methods covariate list). Table 1 tabulates Burn-S as 'Burn' 32 (8.6%) versus 'No burn' 342 (91.4%), so patients without a burn carry the code 2 with the resolved / minor burns. Recoded as DIS_BURN_RECENT = 2 - BurnS (1 = severe burn in the acute phase, 0 = otherwise); model() rebuilds the source code as 2 - DIS_BURN_RECENT so the printed coefficient -0.43 L/h on Burn-S applies unchanged. A severe acute burn therefore raises CL by 0.43 L/h, matching the abstract ('CL ... raised under ... severe burn status') and the Discussion ('burn status increased the CL'). Recency threshold: acute convalescence versus 20 days after burn; severity threshold: 50% burn degree.",
      source_name = "Burn-S"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy status indicator",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Source column CRRT-S, coded 1 = treated by CRRT, 2 = not (Lin 2021 Methods covariate list). Recoded as RRT_CRRT_STATUS = 2 - CRRTS (1 = CRRT, 0 = none); model() rebuilds the source code as 2 - RRT_CRRT_STATUS so the printed coefficient +0.41 L/h on CRRT-S applies unchanged. CRRT therefore LOWERS CL by 0.41 L/h, matching the abstract ('reduced under ... continuous renal replacement therapy status') and the Discussion, which notes that this contradicts earlier reports. 87/374 (23.3%) on CRRT (Table 1). CRRT modality and effluent flow are not reported; trough samples in CRRT patients were drawn after the CRRT session.",
      source_name = "CRRT-S"
    )
  )

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 294L,
    n_studies = 1L,
    age_range = "18-93 years",
    age_median = "62 years",
    weight_range = "40.0-90.6 kg",
    weight_median = "65.0 kg",
    sex_female_pct = 32.1,
    race_ethnicity = "Chinese (single-centre cohort, Taizhou Hospital of Zhejiang Province)",
    disease_state = "Adult ICU patients (> 18 years) treated with vancomycin for suspected or documented MRSA / Gram-positive infection (blood 55.6%, pulmonary 18.2%, skin / soft tissue 11.2%, abdominal 6.4%, other 8.6%). Excluded: < 7 days of vancomycin, children, pregnancy, end-stage renal disease, terminal tumour, hepatic failure.",
    dose_range = "Vancomycin IV 0.5 g qd, 0.5 g q12h, 0.5 g q8h, 0.5 g q6h, 1 g qd or 1 g q12h, chosen by clinicians on condition and renal function. Infusion duration not reported.",
    regions = "China (Taizhou Hospital of Zhejiang Province, Linhai)",
    renal_function = "Serum creatinine median 71.0 umol/L (range 28.0-581.0); 87/374 (23.3%) on CRRT",
    burn_status = "32/374 (8.6%) with burn",
    sampling = "837 trough samples (0.5 h before the next dose) and 156 'approach peak' samples (1 h after the end of infusion), mostly at steady state after the 6th dose; 993 concentrations in total. Serum vancomycin by HPLC (calibration 0.4-100 ug/mL).",
    notes = "Table 1 summarises all 374 patients; 294 were used for model building and 80 (randomly selected by Kinetica) for external validation. A further 92 patients whose steady-state trough was outside 10-20 ug/mL formed the Bayesian dose-adjustment application set. Enrolled July 2015 - December 2017. Kinetica 4.4.1, expectation-maximisation (EM) algorithm, one-compartment model."
  )

  ini({
    # Final covariate equations (Lin 2021 Table 3, Table 4, Results text and
    # Table 7; theta numbering from Table 6):
    #   CL = 0.78 + 0.036*DA - 0.007*Cr - 0.43*Burn-S + 0.067*TBW + 0.41*CRRT-S
    #   V  = 46.47 - 0.04*Cr - 0.07*Age + 0.16*TBW
    # with the source codings DA, Burn-S, CRRT-S each taking the values 1 or 2
    # (see covariateData). Table 6 prints the coefficient magnitudes
    # (theta3 = 0.007, theta4 = 0.43, theta8 = 0.04, theta9 = 0.07); the signs
    # are taken from the printed equations, and each sign agrees with the
    # direction the abstract and Discussion state.
    lcl_int <- log(0.78); label("Intercept of the additive-linear CL equation (L/h)") # Table 6 theta1 = 0.78 (RSEE 33.46%); Table 3 final equation
    lvc_int <- log(46.47); label("Intercept of the additive-linear V equation (L)") # Table 6 theta7 = 46.47 (RSEE 36.21%); Table 4 final equation

    e_conmed_dopa_cl <- 0.036; label("CL slope on the dopamine code DA (1 = none, 2 = dopamine) (L/h)") # Table 6 theta2 = 0.036 (RSEE 36.15%); Table 3 '+0.036*DA'
    e_creat_cl <- -0.007; label("CL slope on serum creatinine (L/h per umol/L)") # Table 6 theta3 = 0.007 (RSEE 32.32%); Table 3 '-0.007*Cr'
    e_dis_burn_recent_cl <- -0.43; label("CL slope on the burn code Burn-S (1 = severe acute burn, 2 = otherwise) (L/h)") # Table 6 theta4 = 0.43 (RSEE 36.51%); Table 3 '-0.43*Burn-S'
    e_wt_cl <- 0.067; label("CL slope on total body weight (L/h per kg)") # Table 6 theta5 = 0.067 (RSEE 22.15%); Table 3 '+0.067*TBW'
    e_rrt_crrt_status_cl <- 0.41; label("CL slope on the CRRT code CRRT-S (1 = CRRT, 2 = none) (L/h)") # Table 6 theta6 = 0.41 (RSEE 35.65%); Table 3 '+0.41*CRRT-S'

    e_creat_vc <- -0.04; label("V slope on serum creatinine (L per umol/L)") # Table 6 theta8 = 0.04 (RSEE 32.43%); Table 4 '-0.04*Cr'
    e_age_vc <- -0.07; label("V slope on age (L per year)") # Table 6 theta9 = 0.07 (RSEE 38.17%); Table 4 '-0.07*Age'
    e_wt_vc <- 0.16; label("V slope on total body weight (L per kg)") # Table 6 theta10 = 0.16 (RSEE 32.63%); Table 4 '+0.16*TBW'

    # Between-subject variability, additive-normal on the linear scale
    # (Methods: 'beta_i = Z_i beta + eta_j'). Table 6 reports the final-model
    # population coefficient of variation relative to the population mean
    # parameter (CL 3.16 L/h, V 60.71 L), so the eta SD is CV * mean:
    #   SD(CL) = 0.4172 * 3.16  = 1.3184 L/h -> variance 1.7381
    #   SD(V)  = 0.3514 * 60.71 = 21.333 L   -> variance 455.11
    # No CL-V covariance is reported; the etas are diagonal.
    etacl ~ 1.7381 # Table 6 final model '%CV of CL' = 41.72%, times CL = 3.16 L/h, squared
    etavc ~ 455.11 # Table 6 final model '%CV of V' = 35.14%, times V = 60.71 L, squared

    # Residual error: the Kinetica EM fit estimated a residual variance but
    # the paper does not report its value or its weighting scheme. Fixed at
    # zero so that simulations return individual predictions; see the
    # vignette.
    addSd <- fixed(0); label("Additive residual error (ug/mL); not reported") # not reported in Lin 2021
  })
  model({
    # Rebuild the source 1/2 codings from the canonical 0/1 indicators so the
    # printed coefficients apply unchanged.
    DA_code <- CONMED_DOPA + 1
    BURNS_code <- 2 - DIS_BURN_RECENT
    CRRTS_code <- 2 - RRT_CRRT_STATUS

    tvcl <- exp(lcl_int) +
      e_conmed_dopa_cl * DA_code +
      e_creat_cl * CREAT +
      e_dis_burn_recent_cl * BURNS_code +
      e_wt_cl * WT +
      e_rrt_crrt_status_cl * CRRTS_code
    tvvc <- exp(lvc_int) + e_creat_vc * CREAT + e_age_vc * AGE + e_wt_vc * WT

    cl <- tvcl + etacl
    vc <- tvvc + etavc

    kel <- cl / vc
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L: central / vc is mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
