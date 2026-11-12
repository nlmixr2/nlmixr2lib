SanchezVidaurre_2025_denosumab_phase1 <- function() {
  description <- paste(
    "One-compartment population PK model with first-order subcutaneous absorption and parallel linear",
    "and Michaelis-Menten elimination for denosumab (MB09 biosimilar and EU- and US-sourced Xgeva,",
    "pooled) after a single 35 mg subcutaneous dose in healthy adult men (phase I study",
    "MB09-A-01-19). Body weight enters linear clearance as a power term centred on 70 kg; log-normal",
    "inter-individual variability on ka, V/F and CL/F; proportional residual error. This is the",
    "phase-I-only fit of Sanchez-Vidaurre 2025 (Table 3, the abstract's headline estimates) used for",
    "the model-based bioequivalence assessment; the companion pooled phase I + phase III final model",
    "is SanchezVidaurre_2025_denosumab."
  )
  reference <- paste(
    "Sanchez-Vidaurre S, Paravisini A, Queiruga-Parada J. Population PK Modeling of Denosumab",
    "Biosimilar MB09 and Reference Denosumab to Establish PK Similarity. Pharmaceutics.",
    "2025;17(9):1146. doi:10.3390/pharmaceutics17091146"
  )
  vignette <- "SanchezVidaurre_2025_denosumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters linear CL/F as a power term (WT/70)^BW~CL (Sanchez-Vidaurre 2025 Table 3 row 'BW~CL';",
        "Methods 2.3.3 'the continuous covariates were centered to 70 kg'). Body weight was included",
        "as a structural covariate on CL/F only (Results 3.2.1)."
      ),
      source_name = "body weight"
    )
  )

  covariatesDataExcluded <- list(
    TRT_MB09 = list(
      description = "MB09 biosimilar treatment-arm indicator (vs EU- or US-sourced Xgeva)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (EU-Xgeva or US-Xgeva)",
      notes = paste(
        "Treatment (source TRT = 1 MB09, 2 EU-Xgeva, 3 US-Xgeva; Table 1) was tested on CL/F, V/F and",
        "KA in three pairwise comparisons and gave no significant improvement in fit (Results 3.2.2,",
        "Supplementary Figure S1), so it is not in the model and all three products share the same",
        "parameters."
      ),
      source_name = "TRT"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "denosumab",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "denosumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 255L,
    n_observations = 4256L,
    n_studies = 1L,
    age_range = "28-55 years (254 of 255 aged 28 to <55)",
    age_mean = "39.5 years (SD 6.90)",
    weight_mean = "82.97 kg (SD 8.492)",
    height_mean = "178.66 cm (SD 6.227)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "Healthy adult male volunteers",
    dose_range = "Single 35 mg subcutaneous dose (MB09, EU-Xgeva or US-Xgeva; 85 per arm)",
    notes = paste(
      "Phase I study MB09-A-01-19 (NCT05299073): randomised, double-blind, three-arm, parallel,",
      "single-dose study (Sanchez-Vidaurre 2025 Table 1 and Table 2). Rich sampling from pre-dose to",
      "day 253. Serum denosumab by validated Meso Scale Discovery electrochemiluminescence assay.",
      "NONMEM 7.4.3, ADVAN13, FOCE with interaction, MU-referenced log-normal IIV. Four of the 255",
      "subjects had ADA-positive samples; ADA was not tested as a covariate."
    )
  )

  ini({
    # Sanchez-Vidaurre 2025 Table 3 ('PopPK parameters of PK model phase I study').
    lka <- log(0.406)
    label("First-order subcutaneous absorption rate constant (1/day)") # Table 3, row 'KA (1/day)' = 0.406 [3.92% RSE, 95% CI 0.375-0.437]
    lvc <- log(9.33)
    label("Apparent volume of distribution V/F (L)") # Table 3, row 'V/F (L)' = 9.33 [1.21% RSE, 95% CI 9.11-9.55]
    lcl <- log(0.123)
    label("Apparent linear clearance CL/F at 70 kg (L/day)") # Table 3, row 'CL/F (L/day)' = 0.123 [3.55% RSE, 95% CI 0.114-0.132]
    e_wt_cl <- 1.32
    label("Power exponent on (WT/70) for CL/F (unitless)") # Table 3, row 'BW~CL' = 1.32 [16.2% RSE, 95% CI 0.901-1.74]

    # Michaelis-Menten arm. Table 3 prints Km in 'ng/mL' and Vm in 'ng/day';
    # the values are encoded in ug/mL (= mg/L) and mg/day. See model() and the
    # vignette Errata for the arithmetic that establishes the corrected units.
    lkm <- log(0.124)
    label("Michaelis-Menten constant Km (ug/mL)") # Table 3, row 'Km (ng/mL)' = 0.124 [9.27% RSE, 95% CI 0.101-0.147]; printed unit corrected to ug/mL
    lvmax <- log(0.139)
    label("Maximum Michaelis-Menten elimination rate Vm (mg/day)") # Table 3, row 'Vm (ng/day)' = 0.139 [5.28% RSE, 95% CI 0.125-0.153]; printed unit corrected to mg/day

    # IIV is reported as a percentage for exponential (log-normal) random
    # effects (Methods 2.3.2, Eq 2); variances are written as omega^2 = log(1 + CV^2).
    etalka ~ log(1 + 0.61^2) # Table 3, KA inter-individual variability 61%
    etalvc ~ log(1 + 0.13^2) # Table 3, V/F inter-individual variability 13%
    etalcl ~ log(1 + 0.37^2) # Table 3, CL/F inter-individual variability 37%

    propSd <- 0.17
    label("Proportional residual error (fraction)") # Table 3, row 'Residual variability' = 17%; proportional model (Results 3.2.1)
  })

  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    km <- exp(lkm)
    vmax <- exp(lvmax)

    # Amounts in mg and volume in L, so Cc is in mg/L (= ug/mL). The
    # saturable arm vmax * Cc / (km + Cc) is an amount rate in mg/day.
    # UNIT CORRECTION on Km and Vm. Read literally (Km 0.124 ng/mL, Vm 0.139
    # ng/day) the Michaelis-Menten arm is either negligible or, with Vm taken
    # as mg/day, a constant 0.139 mg/day sink that drives the 35 mg profile
    # negative. With Km in ug/mL and Vm in mg/day the typical 83 kg subject
    # (Table 2 phase I mean) reaches Cmax about 3,170 ng/mL and AUCinf about
    # 141,600 day*ng/mL, matching the Supplementary Table 2 geometric means
    # (Cmax 3,020-3,188; AUCinf 136,115-148,433); the linear arm alone would
    # give an AUCinf of about 227,000.
    Cc <- central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cc - vmax * Cc / (km + Cc)

    Cc ~ prop(propSd)
  })
}
