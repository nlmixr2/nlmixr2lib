Tan_2025_pazopanib_tumor_sts <- function() {
  description <- paste(
    "Tumor size model for the sum of longest diameters (SLD) of target lesions in real-world adults",
    "with soft tissue sarcoma on oral pazopanib (Tan 2025). Claret-type dynamics with",
    "exponential growth and a first-order decay whose rate falls exponentially with time (acquired",
    "resistance): dSLD/dt = p * SLD - kkill * exp(-res * t) * SLD, started from the observed",
    "baseline SLD. A two-class mixture separates primary-sensitive patients from primary-resistant",
    "patients (13.4%) whose decay rate is zero; the class is supplied as the binary covariate",
    "MIX_PRIM_RESIST. No exposure or dose dependence was found on any rate. Exponential IIV on the",
    "decay rate only (IIV on growth was fixed to zero) and a combined proportional-plus-additive",
    "residual error. The metastatic renal cell carcinoma fit is Tan_2025_pazopanib_tumor_mrcc.",
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
    dosing = "n/a (no dosing events; the model has no exposure dependence)",
    concentration = "mm (sum of longest diameters of target lesions)"
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "sum of longest diameters of target lesions",
      units = "mm",
      specimen = "tumor",
      verified = TRUE
    )
  )

  covariateData <- list(
    TUM_SLD = list(
      description = "Baseline sum of longest diameters of target lesions (RECIST 1.1)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject; sets the initial condition tumor_size(0) (supplement control",
        "stream `A_0(3) = BASE`). The observed baseline was used as an input rather than estimated",
        "(Methods), because the STS tumor data were sparse. STS cohort median 88 mm (IQR 64-148, Table 1).",
        sep = " "
      ),
      source_name = "BASE"
    ),
    MIX_PRIM_RESIST = list(
      description = paste(
        "Mixture-model class indicator: 1 = primary-resistant subpopulation (no drug-induced decay,",
        "kkill = 0); 0 = primary-sensitive subpopulation",
        sep = " "
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (primary-sensitive subpopulation; the majority class)",
      notes = paste(
        "Not a measured covariate: the latent class of the NONMEM $MIXTURE block. MIXNUM = 1 is",
        "the sensitive class with probability THETA(7) = 0.866; MIXNUM = 2 is the resistant class",
        "with KD = THETA(2) = 0 FIX (supplement control stream). MIX_PRIM_RESIST =",
        "as.integer(MIXNUM == 2). For population simulation draw MIX_PRIM_RESIST ~ Bernoulli(0.27)",
        "per subject (Table 4 row 'Subpopulation of primary resistance (%)' 13.4%, RSE 8%).",
        sep = " "
      ),
      source_name = "MIXNUM (NONMEM $MIXTURE class)"
    )
  )

  # Exposure metrics tested on the decay rate (linear, exponential and Emax) and not retained.
  covariatesDataExcluded <- list(
    CTROUGH = list(
      description = "Predicted daily pazopanib trough concentration",
      units = "mg/L",
      type = "continuous",
      notes = "Tested on KD; not significant (Results, 'Exposure-Tumor Size Dynamics Analysis')."
    ),
    CMAX = list(
      description = "Predicted daily pazopanib Cmax",
      units = "mg/L",
      type = "continuous",
      notes = "Tested on KD; not significant."
    ),
    DOSE = list(
      description = "Pazopanib dose",
      units = "mg",
      type = "continuous",
      notes = "Tested on KD; not significant."
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    n_observations = 81L,
    age_range = "21-82 years (median 51), whole STS cohort",
    weight_range = "45-124 kg (median 73), whole STS cohort",
    sex_female_pct = 59.0,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Adults with soft tissue sarcoma on pazopanib at Leiden University Medical",
      "Center, 2014-2023, with a baseline SLD and at least one on-treatment CT scan.",
      sep = " "
    ),
    dose_range = "Oral pazopanib once daily, fasted; initial dose median 800 mg (range 200-800).",
    regions = "The Netherlands (Leiden University Medical Center)",
    notes = paste(
      "111 of 135 patients had tumor size data (403 SLD records); 81 SLD records came from the",
      "STS patients, median 2 per patient (Table 1). The paper does not split the 111 by tumor",
      "type. The demographic ranges and female share (Table 1: 16 male / 23 female of 39) are for",
      "the whole STS PK cohort. NONMEM 7.4.4, FOCE-I.",
      sep = " "
    )
  )

  ini({
    lp <- log(0.0086); label("Exponential tumor growth rate KG (1/day)") # Table 4 STS row 'Growth rate KG (day-1)' 0.0086 (RSE 22%); control stream THETA(1)
    lkkill <- log(0.008); label("First-order tumor decay rate KD in the primary-sensitive class (1/day)") # Table 4 STS row 'Decay rate KD (day-1)' 0.008 (RSE 24%); control stream THETA(3)
    lres <- log(0.0003); label("Acquired resistance rate lambda, the exponential decline of KD over time (1/day)") # Table 4 STS row 'Acquired resistant rate (day-1)' 0.0003 (RSE 59%); control stream THETA(6)

    # Table 4 reports IIV as CV% = sqrt(omega): the control stream $OMEGA
    # final 0.017 (KD) gives 13.0%, matching the table. IIV on KG and on lambda
    # was 0 FIX and is omitted.
    etalkkill ~ 0.017 # Table 4 STS row 'IIV_KD' 13% (RSE 33%); control stream $OMEGA 0.017

    propSd <- 0.166; label("Proportional residual error (fraction)") # Table 4 STS row 'Proportional error (%)' 16.6% (RSE 18%); control stream THETA(4)
    addSd <- 8.1; label("Additive residual error (mm)") # Table 4 STS row 'Additive error (mm)' 8.1 (RSE 71%); control stream THETA(5)
  })
  model({
    p <- exp(lp)
    # Primary-resistant class: KD = THETA(2) = 0 FIX (control stream $PK, MIXNUM = 2).
    kkill <- exp(lkkill + etalkkill) * (1 - MIX_PRIM_RESIST)
    res <- exp(lres)

    # Figure 2 Eq. 7 (sensitive) and Eq. 8 (resistant, kkill = 0).
    d/dt(tumor_size) <- p * tumor_size - kkill * exp(-res * t) * tumor_size
    tumor_size(0) <- TUM_SLD

    # Control stream W = SQRT(THETA(4)**2 * IPRED**2 + THETA(5)**2): combined2.
    tumor_size ~ add(addSd) + prop(propSd) + combined2()
  })
}
