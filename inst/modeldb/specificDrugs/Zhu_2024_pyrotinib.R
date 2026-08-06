Zhu_2024_pyrotinib <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for oral pyrotinib in Chinese patients with HER2-positive breast cancer, with a serum total protein effect on apparent clearance (Zhu 2024)"
  reference <- paste(
    "Zhu Y, Xu Y, Zhao H, Qie H, Gao X, Gao J, Feng Z, Bai J, Feng R, Wang M.",
    "A validated UPLC-MS/MS method for quantification of pyrotinib and population",
    "pharmacokinetic study of pyrotinib in HER2-positive breast cancer patients.",
    "Front Pharmacol. 2024;15:1432944. doi:10.3389/fphar.2024.1432944.",
    "Ka was not estimated here; it was fixed at the value reported by the earlier",
    "pyrotinib population PK model of Wen et al. (2021), as transcribed in",
    "Zhu 2024 Table 4 and Section 4.3.2.",
    sep = " "
  )
  vignette <- "Zhu_2024_pyrotinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "pyrotinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pyrotinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    TPRO = list(
      description        = "Serum total protein (albumin plus globulins)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model (Zhu 2024 Table 5:",
        "TP on CL/F survived backward elimination at p < 0.001; height on CL/F and",
        "on V/F did not). Cohort median 67.2 g/L, range 49.0-80.7 g/L (Table 3).",
        "Entered as the paper's Equation 7 exponential form theta1 * theta2^(cov/cov_median)",
        "with theta2 = exp(0.376), which is NOT median-centred: the covariate factor is",
        "exp(0.376) = 1.46 at the median rather than 1. See the vignette Errata.",
        "Baseline value used throughout (Zhu 2024 collected baseline laboratory",
        "measurements from medical records; no time-varying TP was modelled)."
      ),
      source_name        = "TP"
    )
  )

  # Screened during covariate model building (Zhu 2024 Sections 3.2.2, 4.2.2 and
  # 4.3.2) but NOT retained in the final model. Documentation only -- none of
  # these appear in model(). Beyond the canonical names listed here the paper
  # also screened occurrence of diarrhoea, combined montmorillonite powder,
  # combined loperamide capsules, and oestrogen / progesterone receptor status;
  # those have no canonical covariate-register name and are described in the
  # vignette narrative instead of being given invented names here.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened; not retained. Median 53 years (range 34-68), Table 3."
    ),
    HT = list(
      description = "Height", units = "cm", type = "continuous",
      notes = paste(
        "Entered the model during forward inclusion on both V/F and CL/F",
        "(Table 5, steps 2 and 3) but was dropped at backward elimination",
        "(steps 4 and 5, p > 0.001). Median 158 cm (range 150-166)."
      )
    ),
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      notes = "Screened; not retained. Median 61.5 kg (range 52.5-86), Table 3."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = "Correlated with WT; WT was carried forward instead (Section 4.3.2)."
    ),
    SOD = list(
      description = "Serum sodium", units = "mmol/L", type = "continuous",
      notes = "Screened; not retained. Median 138.6 mmol/L (range 123.2-145.2)."
    ),
    ALB = list(
      description = "Serum albumin", units = "g/L", type = "continuous",
      notes = paste(
        "Screened; not retained. Median 41.4 g/L. A component of TPRO,",
        "which is the covariate the final model retained."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Correlated with ALT (R > 0.7); neither retained (Section 4.3.2)."
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = "Correlated with AST (R > 0.7); neither retained (Section 4.3.2)."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      notes = paste(
        "Carried forward as the representative bilirubin measure (DBIL and IBIL",
        "were correlated with it) but not retained in the final model."
      )
    ),
    CREAT = list(
      description = "Serum creatinine", units = "umol/L", type = "continuous",
      notes = "Screened; not retained. Median 60.5 umol/L (range 30.8-419.0)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50,
    n_studies      = 1,
    n_observations = 158,
    age_range      = "34-68 years",
    age_median     = "53 years",
    weight_range   = "52.5-86 kg",
    weight_median  = "61.5 kg",
    height_range   = "150-166 cm",
    height_median  = "158 cm",
    sex_female_pct = "not reported (HER2-positive breast cancer cohort)",
    race_ethnicity = "not reported (single-centre Chinese cohort)",
    disease_state  = "HER2-positive advanced or metastatic breast cancer",
    dose_range     = "240, 320 or 400 mg pyrotinib maleate tablets orally once daily, taken 30 min after a meal",
    regions        = "China (single centre: The Fourth Hospital of Hebei Medical University, Shijiazhuang)",
    total_protein_range  = "49.0-80.7 g/L",
    total_protein_median = "67.2 g/L",
    notes          = paste(
      "Real-world therapeutic-drug-monitoring cohort recruited November 2020 to",
      "November 2023; opportunistic (sparse) blood sampling. Doses were dynamic:",
      "some patients had dose reductions for intolerable adverse reactions, so the",
      "modelling dataset spans three dose strengths. Baseline demographics are",
      "Zhu 2024 Table 3. NOTE: the mean +/- SD column of Table 3 is scrambled in the",
      "published article (e.g. height reads 73.4 +/- 24.4 cm and weight",
      "163.7 +/- 23.0 kg, which are transposed and physiologically impossible), so",
      "only the median and range values from that table -- which are internally",
      "consistent -- are transcribed here. The TP median of 67.2 g/L is independently",
      "confirmed by the final-model covariate equation.",
      "Concomitant medication reported: montmorillonite powder 9/50 (18%),",
      "loperamide capsules 15/50 (30%); diarrhoea occurred in 24/50 (48%).",
      "ER positive in 28/50 (56%), PR positive in 25/50 (50%)."
    )
  )

  ini({
    # Structural parameters -- Zhu 2024 Table 7 (final model) and the final-model
    # equations in Section 4.3.2 (page 9):
    #   CL/F (L/h) = 88.8 * e^{(TP / 67.2) * 0.376}
    #   V/F (L)    = 3940
    #   KA (1/h)   = 0.357 FIXED
    lka <- fixed(log(0.357)); label("Absorption rate constant (Ka, 1/h), taken from the earlier pyrotinib population PK model of Wen 2021")  # Zhu 2024 Table 7 / Section 4.3.2 final-model equation; sensitivity analysis in Table 6
    lcl <- log(88.8); label("Apparent clearance scale factor theta1 (L/h); CL/F = theta1 * exp(e_tpro_cl * TPRO / 67.2)")  # Zhu 2024 Table 7 CL/F = 88.8 (RSE 15.1%)
    lvc <- log(3940); label("Apparent central volume of distribution (V/F, L)")  # Zhu 2024 Table 7 V/F = 3,940 (RSE 25.8%)

    # Covariate effect -- exponential form (Zhu 2024 Equation 7,
    # theta_i = theta1 * theta2^(cov_i / cov_median), with theta2 = exp(0.376)).
    e_tpro_cl <- 0.376; label("Serum total protein effect on CL/F, exponential scale (unitless)")  # Zhu 2024 Table 7 "TP-CL/F" = 0.376 (RSE 18.0%)

    # Inter-individual variability. Zhu 2024 reports BSV as "%CV" in Table 7
    # (51.5% on CL/F, 96.5% on V/F) and as the same numbers on the fraction scale
    # in the Table 6 sensitivity analysis (0.515, 0.965), i.e. the reported
    # quantity is the omega standard deviation, not the variance. Variances below
    # are omega^2.
    etalcl ~ 0.265225  # 0.515^2; Zhu 2024 Table 7 BSV-CL 51.5 %CV (shrinkage 26.5%)
    etalvc ~ 0.931225  # 0.965^2; Zhu 2024 Table 7 BSV-V 96.5 %CV (shrinkage 23.2%)

    # Residual error -- proportional only (Zhu 2024 Section 4.3.1, Equation 2)
    propSd <- 0.270; label("Proportional residual error (fraction)")  # Zhu 2024 Table 7 ERR-1 27.0 %CV (shrinkage 16.5%)
  })
  model({
    ka <- exp(lka)

    # Apparent clearance. Zhu 2024 prints the final model as
    #   CL/F = 88.8 * e^{(TP / 67.2) * 0.376}
    # which instantiates the paper's own Equation 7 and is deliberately NOT
    # median-centred: at the cohort median TP of 67.2 g/L the covariate factor is
    # exp(0.376) = 1.457, giving CL/F = 129 L/h. That reconciles with the base
    # model (127 L/h, Table 4) and with the previously published pyrotinib
    # estimates (127 L/h, Wen 2021; 141 L/h, product label). See vignette Errata.
    cl <- exp(lcl + etalcl) * exp(e_tpro_cl * TPRO / 67.2)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg / volume in L gives mg/L; x 1000 converts to ng/mL, the scale on
    # which the UPLC-MS/MS assay was validated (1-1,000 ng/mL, Zhu 2024 Section 4.1.2)
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
