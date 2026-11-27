Traynard_2020_remifentanil <- function() {
  description <- "Three-compartment population PK model for intravenous remifentanil in healthy adults, with power effects of age on CL, Q2, Q3 and V3, of lean body mass on CL and V1, and of sex on V3 (Traynard 2020 MonolixSuite tutorial, final run PK_16, fitted to the Minto 1997 data)."
  reference <- paste(
    "Traynard P, Ayral G, Twarogowska M, Chauvin J. Efficient Pharmacokinetic",
    "Modeling Workflow With the MonolixSuite: A Case Study of Remifentanil.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(4):198-210.",
    "doi:10.1002/psp4.12500. Final estimates from the authors' deposited",
    "Monolix run PK_16 (Supporting Information, Supinfo 2).",
    sep = " "
  )
  vignette <- "Traynard_2020_remifentanil"
  units <- list(
    time = "min",
    dosing = "ug",
    concentration = "ng/mL"
  )
  # Supinfo 1 'Data set': TIME in minutes, AMT in ug, RATE in ug/min, DV in
  # ug/L. With volumes in L, central / vc is ug/L = ng/mL and clearances are
  # in L/min.

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect centered at 36.6299 years (Monolix logtAGE = log(AGE/36.6299); main text, covariate section; PK_16.mlxtran [COVARIATE]). 36.6299 is the observation-weighted mean age of the data set. Effects on CL, Q2 (to peripheral1), Q3 (to peripheral2) and V3 (peripheral2).",
      source_name = "AGE"
    ),
    LBM = list(
      description = "Lean body mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect centered at 55.1554 kg (Monolix logtLBM = log(LBM/55.1554); PK_16.mlxtran [COVARIATE]), the observation-weighted mean LBM. Effects on CL and V1. The tutorial does not state the LBM formula; the data set is the Minto 1997 remifentanil data, whose LBM column matches the James formula from weight, height and sex (cohort range 36.5-75.6 kg in the deposited data set).",
      source_name = "LBM"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female) is the paper's reference; see notes",
      notes = "The tutorial data set codes SEX as 1 = male, 0 = female (Supinfo 1, 'Data set'), and Monolix estimates beta_V3_SEX_1 for SEX = 1 with SEX = 0 (female) as the reference. Under the canonical SEXF (1 = female) the effect is applied in model() as exp(e_sex_vp2 * (1 - SEXF)), so a female keeps the typical V3 and a male has V3 x exp(0.271) (about 1.31-fold).",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    central = list(analyte = "remifentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remifentanil", units = "ug", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "remifentanil", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 1L,
    n_observations = 1992L,
    age_range = "20-85 years",
    lbm_range = "36.5-75.6 kg",
    sex_female_pct = 41.5,
    disease_state = "Healthy adult volunteers",
    dose_range = "Constant-rate intravenous infusion of 1-8 ug/kg/min lasting 4-20 min (single infusion per subject)",
    notes = "The tutorial re-analyses the publicly available Minto 1997 remifentanil data set (main text ref. 7: Minto CF et al. Anesthesiology 1997;86:10-23; also distributed as nlme::Remifentanil): dense concentration measurements in 65 healthy adults after a single constant-rate i.v. infusion. Age, LBM and sex counts (27 female, 38 male) and the observation count were read from the deposited data file PK_reminfentanil.csv (Supinfo 2). Estimated with SAEM in Monolix 2019R1."
  )

  ini({
    # Final estimates from the deposited Monolix run PK_16
    # (final_workflow/PK_16/populationParameters.txt, Supinfo 2). PK_16 is
    # the final model in the main text ('Final model overview'). Typical
    # values are for a female (SEX = 0) at AGE 36.6299 years and
    # LBM 55.1554 kg. Time in min, volumes in L, clearances in L/min.
    lcl <- log(2.56943); label("Clearance CL (L/min)") # PK_16 populationParameters.txt Cl_pop = 2.56943
    lvc <- log(4.85867); label("Central volume V1 (L)") # PK_16 V1_pop = 4.85867
    lq <- log(0.0825548); label("Intercompartmental clearance Q2, central-peripheral1 (L/min)") # PK_16 Q2_pop = 0.0825548
    lvp <- log(7.59029); label("Peripheral volume V2, peripheral1 (L)") # PK_16 V2_pop = 7.59029
    lq2 <- log(1.92566); label("Intercompartmental clearance Q3, central-peripheral2 (L/min)") # PK_16 Q3_pop = 1.92566
    lvp2 <- log(8.20430); label("Peripheral volume V3, peripheral2 (L)") # PK_16 V3_pop = 8.20430

    # Covariate coefficients on the log-transformed parameter (Monolix
    # 'coefficient'); for logtAGE / logtLBM these are power exponents.
    e_age_cl <- -0.320405; label("Power exponent of AGE/36.6299 on CL (unitless)") # PK_16 beta_Cl_logtAGE = -0.320405
    e_lbm_cl <- 0.410657; label("Power exponent of LBM/55.1554 on CL (unitless)") # PK_16 beta_Cl_logtLBM = 0.410657
    e_lbm_vc <- 1.26150; label("Power exponent of LBM/55.1554 on V1 (unitless)") # PK_16 beta_V1_logtLBM = 1.26150
    e_age_q <- -0.950579; label("Power exponent of AGE/36.6299 on Q2 (unitless)") # PK_16 beta_Q2_logtAGE = -0.950579
    e_age_q2 <- -0.869000; label("Power exponent of AGE/36.6299 on Q3 (unitless)") # PK_16 beta_Q3_logtAGE = -0.869000
    e_age_vp2 <- -0.462430; label("Power exponent of AGE/36.6299 on V3 (unitless)") # PK_16 beta_V3_logtAGE = -0.462430
    e_sex_vp2 <- 0.271176; label("Log-scale effect of male sex on V3 (unitless; applied as exp(e_sex_vp2 * (1 - SEXF)))") # PK_16 beta_V3_SEX_1 = 0.271176

    # IIV: Monolix omega is the SD of the normal random effect on the log
    # scale; variances are omega^2 and covariances r * omega_i * omega_j.
    # Block 1 (Cl, Q3, V3): omega_Cl 0.142692, omega_Q3 0.412362,
    # omega_V3 0.313035, corr_Q3_Cl 0.675947, corr_V3_Cl 0.709079,
    # corr_V3_Q3 0.821229 (PK_16 populationParameters.txt).
    etalcl + etalq2 + etalvp2 ~ c(
      0.0203610,
      0.0397732, 0.170042,
      0.0316729, 0.106007, 0.0979912
    )
    # Block 2 (Q2, V2): omega_Q2 0.385790, omega_V2 0.853822,
    # corr_V2_Q2 0.650297 (PK_16 populationParameters.txt).
    etalq + etalvp ~ c(
      0.148834,
      0.214205, 0.729012
    )
    etalvc ~ 0.0700334 # PK_16 omega_V1 = 0.264638; variance = 0.264638^2

    propSd <- 0.112461; label("Proportional residual error (fraction)") # PK_16 b = 0.112461 (errorModel = proportional(b))
  })
  model({
    sex_male <- 1 - SEXF

    cl <- exp(lcl + etalcl) * (AGE / 36.6299)^e_age_cl * (LBM / 55.1554)^e_lbm_cl
    vc <- exp(lvc + etalvc) * (LBM / 55.1554)^e_lbm_vc
    q <- exp(lq + etalq) * (AGE / 36.6299)^e_age_q
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2 + etalq2) * (AGE / 36.6299)^e_age_q2
    vp2 <- exp(lvp2 + etalvp2) * (AGE / 36.6299)^e_age_vp2 * exp(e_sex_vp2 * sex_male)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Monolix library model infusion_3cpt_ClV1Q2V2Q3V3: i.v. infusion into
    # the central compartment, linear elimination, two peripheral
    # compartments.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
