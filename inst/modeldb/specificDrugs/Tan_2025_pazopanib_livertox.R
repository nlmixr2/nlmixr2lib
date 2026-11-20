Tan_2025_pazopanib_livertox <- function() {
  description <- paste(
    "Parametric time-to-first-event model for CTCAE grade >= 2 liver toxicity (ALT, AST or alkaline",
    "phosphatase elevation) in real-world adults with metastatic renal cell carcinoma or soft tissue",
    "sarcoma on oral pazopanib (Tan 2025). Gompertz baseline hazard h(t) = lambda * exp(shape * t)",
    "with a negative shape (hazard falling over time), multiplied by exp(1.21) = 3.35 while the",
    "predicted steady-state pazopanib trough exceeds 34 mg/L. The trough is supplied as the",
    "time-varying data column CTROUGH, which the source derived from the companion population PK",
    "model (Tan_2025_pazopanib). Outputs the hazard, the cumulative hazard and the toxicity-free",
    "survival probability `sur`; follow-up was censored at 365 days.",
    sep = " "
  )
  reference <- paste(
    "Tan Z, Voller S, Yin A, Rieborn A, Gelderblom H, van der Hulle T, Knibbe CAJ, Moes DJAR.",
    "Model-Informed Dose Optimization of Pazopanib in Real-World Patients with Cancer.",
    "Clin Pharmacokinet. 2025;64(5):715-728. doi:10.1007/s40262-025-01504-5.",
    sep = " "
  )
  vignette <- "Tan_2025_pazopanib"
  units <- list(
    time = "day",
    dosing = "n/a (no dosing events; pazopanib exposure enters as the data column CTROUGH)",
    concentration = "probability (the model output `sur` is a toxicity-free survival probability)"
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "cumulative hazard of CTCAE grade >= 2 liver toxicity",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Predicted steady-state pazopanib trough plasma concentration on the current day",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "TOTAL pazopanib, time-varying. The source computed it per patient and day from the",
        "individual (post hoc) parameters of the companion population PK model Tan_2025_pazopanib",
        "('the trough concentration on the day of the toxicity event or the day of the censoring',",
        "Methods). Enters the hazard as a step: the hazard is multiplied by exp(1.21) = 3.35 when",
        "CTROUGH > 34 mg/L and by 1 otherwise (supplement control stream",
        "`IF(trough.GT.34) RF = EXP(THETA(3))`). The cohort median predicted Cmin,ss was 27 mg/L",
        "(IQR 22-35, Table 1).",
        sep = " "
      ),
      source_name = "trough"
    )
  )

  # Covariates tested on the base hazard and not retained (Results, 'Exposure-Liver Toxicity
  # Analysis'); only the most significant covariate (Cmin,ss) was kept.
  covariatesDataExcluded <- list(
    DOSE = list(
      description = "Pazopanib dose on the day before the toxicity event, and initial dose",
      units = "mg",
      type = "continuous",
      notes = "Both tested on the base hazard and not retained; correlated with the exposure metrics."
    ),
    CMAX = list(
      description = "Predicted pazopanib Cmax on the event / censoring day",
      units = "mg/L",
      type = "continuous",
      notes = "Tested as the exposure metric; Cmin,ss gave the larger drop in OFV (8.03) and was selected."
    ),
    AUC = list(
      description = "Predicted pazopanib AUC over 24 h on the event / censoring day",
      units = "mg*h/L",
      type = "continuous",
      notes = "Tested as the exposure metric; not selected."
    ),
    AGE = list(description = "Age", units = "years", type = "continuous", notes = "Tested, not retained."),
    WT = list(description = "Body weight", units = "kg", type = "continuous", notes = "Tested, not retained."),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Tested, not retained."
    ),
    TUMTP_SARC = list(
      description = "Sarcoma tumor-type indicator (1 = soft tissue sarcoma, 0 = metastatic renal cell carcinoma)",
      units = "(binary)",
      type = "binary",
      notes = "Tested, not retained; the Kaplan-Meier curves for mRCC and STS did not differ (Supplementary Fig. S6b)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 135L,
    n_studies = 1L,
    n_events = 27L,
    age_range = "21-89 years (median 67)",
    weight_range = "45-130 kg (median 79)",
    sex_female_pct = 37.0,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Adults with metastatic renal cell carcinoma (n = 96) or soft tissue sarcoma (n = 39) on",
      "pazopanib at Leiden University Medical Center, 2014-2023. 27 of 135 (20%) had a CTCAE v5.0",
      "grade >= 2 ALT, AST or alkaline phosphatase elevation within 1 year of starting pazopanib",
      "(Table 1: 11 grade 2, 16 grade >= 3).",
      sep = " "
    ),
    dose_range = "Oral pazopanib once daily, fasted; initial dose median 800 mg (range 200-800).",
    regions = "The Netherlands (Leiden University Medical Center)",
    notes = paste(
      "Time to first event from the start of pazopanib; censored at the last laboratory record",
      "before the end of treatment or at 365 days. NONMEM 7.4.4, FO with LIKELIHOOD (Methods).",
      sep = " "
    )
  )

  ini({
    llam_haz <- log(0.0021); label("Gompertz baseline hazard at treatment start, lambda (1/day)") # Table 3 row 'Lamda (1/day)' 0.0021 (RSE 33%); control stream THETA(1)
    shape_haz <- -0.012; label("Gompertz shape, the exponential time trend of the hazard (1/day)") # Table 3 row 'SHAPE (1/day)' -0.012 (RSE 28%); control stream THETA(2)
    e_ctrough_haz <- 1.21; label("Log hazard ratio while CTROUGH is above the 34 mg/L cutoff (unitless)") # Table 3 row 'Cmin,ss coefficient (coef)' 1.21 (RSE 32%); exp(1.21) = 3.35; control stream THETA(3)
    ctrough_cut <- fixed(34); label("Trough cutoff above which the hazard is raised (mg/L)") # Table 3 footnote and supplement control stream 'IF(trough.GT.34)'; chosen by a grid over 20.5-43 mg/L (Supplementary Table S1)

    # The likelihood is the survival / event density ($ESTIMATION LIKE), with
    # no residual error. A fixed placeholder lets the model compile.
    addSd_sur <- fixed(0.001); label("Placeholder additive residual SD on the survival probability (unitless); not from the source")
  })
  model({
    lam_haz <- exp(llam_haz)

    # Supplement control stream: RF = EXP(THETA(3)) if trough > 34, else 1.
    rf_ctrough <- exp(e_ctrough_haz * (CTROUGH > ctrough_cut))

    # Gompertz hazard, DADT(1) = LAM * RF * EXP(SHP * T).
    hazard <- lam_haz * rf_ctrough * exp(shape_haz * t)
    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    sur <- exp(-cumhaz)

    sur ~ add(addSd_sur)
  })
}
