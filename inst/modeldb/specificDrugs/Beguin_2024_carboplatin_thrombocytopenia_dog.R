Beguin_2024_carboplatin_thrombocytopenia_dog <- function() {
  description <- paste(
    "Veterinary (client-owned dog).",
    "Sigmoidal Emax exposure-toxicity model for the fractional fall in",
    "thrombocyte count 14 days after a single dose of carboplatin, in",
    "a SECOND cohort of 14 client-owned dogs that took no part in",
    "fitting the PK model. Beguin 2024 fitted two separate Hill",
    "regressions to the same endpoint (Table 2): one driven by the",
    "model-predicted AUC0-inf, which fitted better (r2 = 0.73), and",
    "one driven by the administered dose in mg/kg (r2 = 0.57). Both",
    "are carried here, as thromboRedAuc and thromboRedDose, so the",
    "comparison the paper draws can be reproduced. The one-compartment",
    "PK layer that generates the exposure is reproduced verbatim from",
    "Beguin 2024 Table 1 and is identical to the companion model",
    "Beguin_2024_carboplatin_dog; it is carried inline so this model",
    "runs end to end. Dose in mg/kg is supplied as the covariate",
    "column DOSE_CARBOPLATIN_MGKG rather than read from the event",
    "table, because the two Hill arms are alternative predictors of",
    "one endpoint and the dose arm must stay available even when the",
    "PK layer is not solved. The Hill regressions were fitted by",
    "non-linear least squares in GraphPad Prism, so they carry no",
    "residual-error model and no random effects - only the PK layer",
    "is stochastic."
  )
  reference <- paste(
    "Beguin J, Mahfoudhi S, Uzel M, Rostang A, Ibish C, Ferran AA,",
    "Pelligand L, Hulin A, Kohlhauer M. (2024).",
    "Population pharmacokinetics modelling for clinical dose",
    "adjustment of carboplatin in dogs.",
    "BMC Veterinary Research 20:575.",
    "doi:10.1186/s12917-024-04404-1.",
    sep = " "
  )
  vignette <- "Beguin_2024_carboplatin_dog"

  # Bookkeeping state that integrates Cc so the exposure-response layer can
  # read AUC off the solve, reproducing the AUC output the authors added by
  # hand to their Monolix model. Same idiom as auc_central in
  # Assmus_2025_benznidazole_qpcr.R; not a biological compartment.
  paper_specific_compartments <- c("auc_central")


  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ug/L"
  )

  compartmentData <- list(
    central     = list(analyte = "carboplatin", units = "ug", specimen = "plasma", verified = TRUE),
    auc_central = list(analyte = "carboplatin", units = "ug*h/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Plasma creatinine concentration on the day of carboplatin administration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Uncentred; enters clearance as CREAT^-0.25 per Beguin 2024",
        "Eq. (3) and Table 1. In this second, toxicity cohort the",
        "mean was 10.9 +/- 3.6 mg/L, appreciably higher than the",
        "7.32 +/- 1.86 mg/L of the PK cohort."
      ),
      source_name        = "plasma creatinine concentration"
    ),
    NEUTERED = list(
      description        = "Neutering status (1 = surgically neutered, i.e. spayed female or castrated male; 0 = sexually intact)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sexually intact)",
      notes              = paste(
        "Beguin 2024 'sterilization status'. In the toxicity cohort",
        "only 2 of 14 dogs were neutered, both female; all 9 males",
        "were intact."
      ),
      source_name        = "neutered status"
    ),
    DOSE_CARBOPLATIN_MGKG = list(
      description        = "Administered carboplatin dose per kilogram of body weight",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The explanatory variable X of the dose arm of Beguin 2024",
        "Eq. (4) and the second column of Table 2. Carried as a",
        "covariate rather than derived from the event table so that",
        "the dose arm stays evaluable independently of the PK",
        "solve, and because the paper's dose arm was fitted to the",
        "PRESCRIBED dose, not to an amount reconstructed from the",
        "infusion record. Toxicity cohort mean 11.5 +/- 2.4 mg/kg,",
        "corresponding to 306.1 +/- 22.6 mg/m2."
      ),
      source_name        = "Dose (mg/kg)"
    ),
    OCC = list(
      description        = "Occasion index; retained so the PK layer matches the companion model",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Every dog in the toxicity cohort contributed a single",
        "carboplatin administration, so OCC = 1 throughout this",
        "cohort. The five-occasion structure is kept so that the PK",
        "layer is bit-identical to Beguin_2024_carboplatin_dog."
      ),
      source_name        = "occasion"
    )
  )

  population <- list(
    species      = "dog (client-owned, mixed breeds)",
    n_subjects   = 14L,
    n_studies    = 1L,
    age_range    = "mean 9.34 +/- 2.91 years",
    weight_range = "mean 22.43 +/- 11.28 kg",
    disease_state = "Solid tumours treated with carboplatin; a separate cohort from the 16 dogs used to fit the PK model",
    dose_range   = "mean 306.1 +/- 22.6 mg/m2, equal to 11.5 +/- 2.4 mg/kg, single administration",
    regions      = "France (National Veterinary School of Alfort and Oniris VetAgroBio)",
    renal_function = "Mean plasma creatinine 10.9 +/- 3.6 mg/L",
    notes        = paste(
      "Recruited January 2022 to June 2023. 17 dogs enrolled, 3",
      "excluded for loss to follow-up, 14 analysed. Thrombocyte and",
      "neutrophil counts were taken on the day of treatment and again",
      "at day 14, the expected nadir. Mean observed reduction was",
      "66 +/- 16% for thrombocytes and 65 +/- 25% for neutrophils.",
      "AUC0-inf was predicted per dog in Simulx as the median of 1000",
      "Monte-Carlo replicates. Neither AUC0-inf nor dose correlated",
      "with the NEUTROPHIL reduction, so no neutrophil arm is",
      "reported by the paper and none is encoded here. Predicted",
      "AUC0-inf median 3342 [3121-4017]; see the vignette Errata for",
      "the unit label. See Beguin 2024 Results, 'Toxicity prediction'."
    )
  )

  ini({
    # --- PK layer, reproduced verbatim from Beguin 2024 Table 1 so that this
    # model generates its own exposure. Identical to the companion model
    # Beguin_2024_carboplatin_dog; see that file for the full source trace.
    lvc <- log(4.05)
    label("Volume of distribution (L/kg)")                                  # Table 1 V = 4.05 L/kg [RSE 6.04%]
    lcl <- log(6.9)
    label("Clearance (L/h/kg)")                                             # Table 1 Cl = 6.9 L/h/kg [RSE 24.7%]

    e_creat_cl <- -0.25
    label("Exponent of plasma creatinine on clearance (unitless)")          # Table 1 beta_creatinine = -0.25 [RSE 50.3%]
    e_neutered_cl <- -0.22
    label("Log-scale effect of neutering on clearance (unitless)")          # Table 1 beta_neutered = -0.22 [RSE 34.2%]

    etalvc ~ fixed(0.01)                                                    # Table 1 omega_V = 0.1 SD, held constant by the authors per footnote a; 0.1^2
    etalcl ~ 0.005476                                                       # Table 1 omega_Cl = 0.074 SD [RSE 53.1%]; 0.074^2

    etaiov_vc_1 ~ fixed(0.01)                                               # Table 1 gamma_V = 0.1 SD, held constant by the authors per footnote a; 0.1^2
    etaiov_vc_2 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_vc_3 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_vc_4 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_vc_5 ~ fixed(0.01)                                               # shared with occasion 1
    etaiov_cl_1 ~ 0.0121                                                    # Table 1 gamma_Cl = 0.11 SD [RSE 28.3%]; 0.11^2
    etaiov_cl_2 ~ fixed(0.0121)                                             # shared with occasion 1
    etaiov_cl_3 ~ fixed(0.0121)                                             # shared with occasion 1
    etaiov_cl_4 ~ fixed(0.0121)                                             # shared with occasion 1
    etaiov_cl_5 ~ fixed(0.0121)                                             # shared with occasion 1

    addSd <- 12.28
    label("Additive residual error (ug/L)")                                 # Table 1 a = 12.28 [RSE 17.2%]
    propSd <- 0.23
    label("Proportional residual error (fraction)")                         # Table 1 b = 0.23 [RSE 12.2%]

    # --- Exposure-driven arm of the toxicity model (Beguin 2024 Eq. 4 with
    # X = AUC0-inf; Table 2, first estimate column). Fitted by non-linear
    # least squares in GraphPad Prism 10.2.3 to 14 dogs; r2 = 0.73.
    lemax_auc <- log(0.99)
    label("Maximum fractional thrombocyte reduction, exposure arm (fraction)")   # Table 2 Emax = 0.99 [95% CI 0.76-1]
    leauc50 <- log(2679)
    label("AUC0-inf giving half the maximum thrombocyte reduction (ug*h/L)")     # Table 2 E50 = 2679 [95% CI 2221-2886]
    lhill_auc <- log(3.35)
    label("Hill coefficient, exposure arm (unitless)")                           # Table 2 Hill coefficient = 3.35 [95% CI 2.14-7.44]

    # --- Dose-driven arm of the toxicity model (Beguin 2024 Eq. 4 with
    # X = dose in mg/kg; Table 2, second estimate column). The paper reports
    # this as the inferior comparator, r2 = 0.57.
    lemax_dose <- log(0.76)
    label("Maximum fractional thrombocyte reduction, dose arm (fraction)")       # Table 2 Emax = 0.76 [95% CI 0.68-1]
    led50_dose <- log(9.0)
    label("Dose giving half the maximum thrombocyte reduction (mg/kg)")          # Table 2 E50 = 9.0 [95% CI 7.64-9.32]
    lhill_dose <- log(17.80)
    label("Hill coefficient, dose arm (unitless)")                               # Table 2 Hill coefficient = 17.80 [95% CI 3.47-47.14]
  })

  model({
    # 1. PK layer (Beguin 2024 Eq. 2, 3 and Table 1).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)

    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 +
      oc4 * etaiov_vc_4 + oc5 * etaiov_vc_5
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    vc <- exp(lvc + etalvc + iov_vc)
    cl <- exp(lcl + etalcl + iov_cl) * CREAT^e_creat_cl *
      exp(e_neutered_cl * NEUTERED)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    Cc <- central / vc

    # 2. Cumulative exposure. Carboplatin disposition is linear and
    #    mono-exponential with a half-life under an hour, so auc_central
    #    reaches AUC0-inf for practical purposes within the first day and
    #    is the X of the exposure arm below.
    d/dt(auc_central) <- Cc

    # 3. Toxicity arms (Beguin 2024 Eq. 4,
    #    Y = Emax * X^n / (E50^n + X^n)). Y is the FRACTION of thrombocyte
    #    reduction between day 0 and day 14, bounded 0 to 1. Both arms
    #    predict the same endpoint from different explanatory variables and
    #    are alternatives, not additive contributions.
    emax_auc <- exp(lemax_auc)
    eauc50 <- exp(leauc50)
    hill_auc <- exp(lhill_auc)
    emax_dose <- exp(lemax_dose)
    ed50_dose <- exp(led50_dose)
    hill_dose <- exp(lhill_dose)

    thromboRedAuc <- emax_auc * auc_central^hill_auc /
      (eauc50^hill_auc + auc_central^hill_auc)
    thromboRedDose <- emax_dose * DOSE_CARBOPLATIN_MGKG^hill_dose /
      (ed50_dose^hill_dose + DOSE_CARBOPLATIN_MGKG^hill_dose)

    Cc ~ add(addSd) + prop(propSd)
  })
}
