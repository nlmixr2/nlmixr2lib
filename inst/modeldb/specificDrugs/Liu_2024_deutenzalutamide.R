Liu_2024_deutenzalutamide <- function() {
  description <- "Two-compartment population PK model with first-order absorption for deutenzalutamide (deuterated enzalutamide, HC-1119) in patients with metastatic castration-resistant prostate cancer (Liu 2024)"
  reference <- "Liu Y, He Y, Qi X, Li X, Zhou Y, Chen Y, Wang Z, Zheng L. Population Pharmacokinetics Modeling and Simulation of Deutenzalutamide, A Novel Androgen Receptor Antagonist, in Patients With Metastatic Castration-Resistant Prostate Cancer. Clin Pharmacol Drug Dev. 2024;13(12):1291-1300. doi:10.1002/cpdd.1477"
  vignette <- "Liu_2024_deutenzalutamide"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot       = list(analyte = "deutenzalutamide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "deutenzalutamide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "deutenzalutamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. The only covariate retained in the final model.",
        "Enters V2/F as the power form of Liu 2024 Equation 3,",
        "theta_i = theta_TV * (COV_i / COV_bar)^theta_2, with COV_bar = 64.0 kg",
        "(the cohort median, Table S3, and the 'typical value of body weight (64.0 kg)'",
        "named in Liu 2024 Results / Simulation). Liu 2024 Methods prose calls Equation 3",
        "the 'exponential mode', but the printed equation is a power (ratio-exponent) model;",
        "the printed equation governs.",
        sep = " "
      ),
      source_name        = "BW"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Baseline age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 005: dOFV = -3.484, short of the 3.84 forward-inclusion criterion) and in the general covariate screen; not retained in the final model."
    ),
    HT = list(
      description = "Baseline height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 006: dOFV = -3.573, short of 3.84); not retained."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 008: dOFV = -0.052) and discussed in Liu 2024 Discussion as a renal-function marker; not retained. Cohort values were all within the normal range, so the final model should not be extrapolated to renal impairment."
    ),
    WBC = list(
      description = "Baseline leukocyte count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 009: dOFV = -0.036); not retained."
    ),
    NEUT = list(
      description = "Baseline neutrophil count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 010: dOFV = -0.046); not retained."
    ),
    PLT = list(
      description = "Baseline platelet count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 011: dOFV = -0.059); not retained."
    ),
    TBIL = list(
      description = "Baseline total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 012: dOFV = -1.419); not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 013: dOFV = 0); not retained. Cohort values were within the normal range, so the final model should not be extrapolated to hepatic impairment."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 014: dOFV = -0.040); not retained."
    ),
    HGB = list(
      description = "Baseline hemoglobin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened on ka (Table S4, run 047: dOFV = -0.124); not retained."
    ),
    ECOG = list(
      description = "Eastern Cooperative Oncology Group performance status grade",
      units       = "(grade)",
      type        = "categorical",
      notes       = "Listed in Liu 2024 Methods (Covariate Model) among the screened covariates; not retained in the final model. Cohort was ECOG 0 (n = 5) or 1 (n = 19) only (Table S3)."
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 24,
    n_studies       = 1,
    n_observations  = 543,
    age_range       = "55-80 years",
    age_median      = "69.5 years",
    weight_range    = "46-81 kg",
    weight_median   = "64 kg",
    height_range    = "149-178 cm",
    height_median   = "165 cm",
    sex_female_pct  = 0,
    disease_state   = "metastatic castration-resistant prostate cancer (mCRPC)",
    dose_range      = "40 mg (n = 3), 80 mg (n = 9), 160 mg (n = 9), or 200 mg (n = 3) oral soft capsules once daily for 12 weeks",
    regions         = "China",
    performance_status = "ECOG 0 in 5 subjects (20.8%), ECOG 1 in 19 subjects (79.2%)",
    renal_function  = "Serum creatinine 0.41-1.56 mg/dL (median 0.81); all subjects within the normal range",
    hepatic_function = "Total bilirubin 0.33-1.18 mg/dL (median 0.60); ALT 9-34 U/L (median 19); all subjects within the normal range",
    notes           = paste(
      "Open-label 3 + 3 dose-escalation Phase Ia trial NCT03774056 (Liu 2024 Methods,",
      "Study Design and Population). Baseline demographics from Liu 2024 Table S3.",
      "4% of deutenzalutamide concentrations were below the 40 ng/mL lower limit of",
      "quantification and were handled with the Beal M1 method. Because renal and",
      "hepatic laboratory values were all within the normal range, Liu 2024 Discussion",
      "cautions that the model should not be extrapolated to hepatic or renal impairment.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Liu 2024 Table 1, "Final model" estimates column.
    # V1/F and CLd/F carry no RSE, no bootstrap median and no bootstrap CI in
    # Table 1, and the table marks each "(fixed)".
    lka <- log(1.32)         ; label("Absorption rate constant (1/h)")                                        # Liu 2024 Table 1 (RSE 16.3%)
    lcl <- log(0.166)        ; label("Apparent clearance (L/h)")                                            # Liu 2024 Table 1 (RSE 3.9%)
    lvc <- fixed(log(17.5))  ; label("Apparent central volume of distribution (L)")                         # Liu 2024 Table 1, marked "(fixed)"
    lvp <- log(58.1)         ; label("Apparent peripheral volume of distribution at 64 kg body weight (L)") # Liu 2024 Table 1 (RSE 8.3%)
    lq  <- fixed(log(12.5))  ; label("Apparent intercompartmental clearance (L/h)")                        # Liu 2024 Table 1, marked "(fixed)"

    # Covariate effect. Power (ratio-exponent) form of Liu 2024 Equation 3,
    # theta_i = theta_TV * (COV_i / COV_bar)^theta_2, with COV_bar = 64.0 kg.
    e_wt_vp <- 1.83          ; label("Body-weight exponent on V2/F (unitless)")                                   # Liu 2024 Table 1, "BW on V2/F" (RSE 25.6%)

    # Inter-individual variability. Liu 2024 Equation 1 is exponential,
    # P_i = P_hat * exp(eta_i). Table 1 labels these rows "omega", but the
    # tabulated numbers are the NONMEM OMEGA *variances*: read as standard
    # deviations they would imply 2.7% CV on CL/F and 2.1% proportional
    # residual error, both below the assay's own reported imprecision
    # (intra-assay 0.6-4.7%, inter-assay 3.4-4.9%, Liu 2024 Bioanalytical
    # Methods). Read as variances they give the plausible values below.
    etalcl ~ 0.0274  # Liu 2024 Table 1, omega CL/F  (RSE 41.6%, shrinkage 2.40%); sqrt(exp(0.0274) - 1) ~= 16.7% CV
    etalvp ~ 0.0446  # Liu 2024 Table 1, omega V2/F  (RSE 38.1%, shrinkage 6.40%); sqrt(exp(0.0446) - 1) ~= 21.3% CV
    etalka ~ 0.440   # Liu 2024 Table 1, omega ka    (RSE 40.2%, shrinkage 4.80%); sqrt(exp(0.440) - 1) ~= 74.3% CV

    # Residual error. Liu 2024 Population Pharmacokinetic Modeling: "the
    # residual variability was fitted by proportional model", i.e. only the
    # eps_1 term of Equation 2 was retained. Table 1 tabulates sigma^2.
    propSd <- sqrt(0.0208); label("Proportional residual error (fraction)")  # Liu 2024 Table 1, sigma prop (RSE 2.50%, shrinkage 5.80%); sqrt(0.0208) ~= 0.144
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    # Body weight acts only on the peripheral volume (Liu 2024 Results:
    # "Covariates with significant effects on pharmacokinetic parameters were
    # selected by body weight on V2").
    vp <- exp(lvp + etalvp) * (WT / 64)^e_wt_vp
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L = ug/mL
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
