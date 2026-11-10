SanchezVidaurre_2025_denosumab <- function() {
  description <- paste(
    "One-compartment population PK model with first-order subcutaneous absorption and parallel linear",
    "and Michaelis-Menten elimination for denosumab (MB09 biosimilar and reference Xgeva / Prolia,",
    "pooled), fitted to a phase I study in healthy adult men (single 35 mg SC dose) and a phase III",
    "study in postmenopausal women with osteoporosis (60 mg SC every 6 months). Body weight enters",
    "linear clearance as a power term centred on 70 kg and study enters V/F; log-normal",
    "inter-individual variability on ka, V/F and CL/F; proportional residual error. Final model of",
    "Sanchez-Vidaurre 2025 (Table 4). The study effect on V/F (-0.534) is encoded as printed with",
    "the phase I study as reference, but the paper does not state its form or coding and this",
    "reading does not reproduce the paper's own phase III figures; see the vignette before using",
    "STUDY_MB09C0119 = 1. The phase-I-only fit is SanchezVidaurre_2025_denosumab_phase1."
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
        "Enters linear CL/F as a power term (WT/70)^BW~CL (Sanchez-Vidaurre 2025 Table 4 row 'BW~CL';",
        "Methods 2.3.3 'the continuous covariates were centered to 70 kg'). Body weight was included",
        "on CL/F as a structural covariate (Results 3.2.3)."
      ),
      source_name = "body weight"
    ),
    STUDY_MB09C0119 = list(
      description = "Phase III study MB09-C-01-19 (SIMBA) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase I study MB09-A-01-19, healthy men)",
      notes = paste(
        "1 = postmenopausal woman with osteoporosis from the phase III study MB09-C-01-19, 0 =",
        "healthy man from the phase I study MB09-A-01-19. Enters V/F as exp(-0.534 * STUDY_MB09C0119)",
        "(Table 4 row 'Study~V'). The paper prints neither the functional form nor the study coding",
        "nor the reference study; the exponential MU-referenced form with the phase I study as",
        "reference was chosen because all parameters were MU-referenced (Results 3.2.1) and the",
        "pooled V/F of 9.00 L sits next to the phase-I-only 9.33 L. Under this reading phase III V/F",
        "is 5.27 L, which over-predicts the peak of the paper's own phase III simulation (Figure 3)",
        "by about 1.7-fold; with the indicator at 0 the peaks of Figures 3 and 4 are reproduced,",
        "though their terminal phase falls between the two readings.",
        "The question has been sent to the corresponding author."
      ),
      source_name = "Study"
    )
  )

  covariatesDataExcluded <- list(
    TRT_MB09 = list(
      description = "MB09 biosimilar treatment-arm indicator (vs reference denosumab)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (EU-Xgeva, US-Xgeva or EU-Prolia)",
      notes = paste(
        "Treatment (source TRT; Table 1) was tested on CL/F, V/F and KA and gave no significant",
        "improvement in fit (Results 3.2.2, Supplementary Figure S1), so it is not in the model."
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
    n_subjects = 810L,
    n_observations = 7282L,
    n_studies = 2L,
    age_mean = "39.5 years (SD 6.90) phase I; 65.8 years (SD 5.94) phase III",
    weight_mean = "82.97 kg (SD 8.492) phase I; 63.196 kg (SD 8.787) phase III",
    sex_female_pct = 68.5,
    race_ethnicity = c(White = 99.5, `American Indian or Alaska Native` = 0.5),
    disease_state = paste(
      "Healthy adult men (phase I MB09-A-01-19, n = 255) and postmenopausal women with osteoporosis",
      "(phase III MB09-C-01-19, n = 555)"
    ),
    dose_range = paste(
      "Single 35 mg SC dose (phase I; MB09, EU-Xgeva or US-Xgeva) or 60 mg SC every 6 months",
      "(phase III; MB09 or EU-Prolia)"
    ),
    notes = paste(
      "Pooled from Sanchez-Vidaurre 2025 Tables 1-2 and Results 3.1: 4,256 phase I samples (rich, to",
      "day 253) and 3,026 phase III samples (sparse: day 1 pre-dose, day 11, months 1, 3, 6 and 12).",
      "Female percentage = 555/810. Serum denosumab by validated Meso Scale Discovery",
      "electrochemiluminescence assay. NONMEM 7.4.3, ADVAN13, FOCE with interaction, MU-referenced",
      "log-normal IIV. Ten of 810 subjects had ADA-positive samples; ADA was not tested as a covariate."
    )
  )

  ini({
    # Sanchez-Vidaurre 2025 Table 4 ('PopPK parameters of final PK model, both studies').
    lka <- log(0.349)
    label("First-order subcutaneous absorption rate constant (1/day)") # Table 4, row 'KA (1/day)' = 0.349 [2.7% RSE, 95% CI 0.330-0.368]
    lvc <- log(9.00)
    label("Apparent volume of distribution V/F in the phase I study (L)") # Table 4, row 'V/F (L)' = 9.00 [0.8% RSE, 95% CI 8.86-9.14]
    lcl <- log(0.143)
    label("Apparent linear clearance CL/F at 70 kg (L/day)") # Table 4, row 'CL/F (L/day)' = 0.143 [2.6% RSE, 95% CI 0.136-0.150]
    e_wt_cl <- 1.07
    label("Power exponent on (WT/70) for CL/F (unitless)") # Table 4, row 'BW~CL' = 1.07 [9.8% RSE, 95% CI 0.864-1.28]
    e_study_mb09c0119_vc <- -0.534
    label("Effect of the phase III study on V/F, log-scale (unitless)") # Table 4, row 'Study~V' = -0.534 [33.7% RSE, 95% CI -0.887 to -0.181]; form and coding unprinted, see covariateData

    # Michaelis-Menten arm. Table 4 prints Km in 'ng/mL' and Vm in 'ng/day';
    # the values are encoded in ug/mL (= mg/L) and mg/day, as in the phase I
    # model (SanchezVidaurre_2025_denosumab_phase1, where the corrected units
    # reproduce Supplementary Table 2).
    lkm <- log(0.162)
    label("Michaelis-Menten constant Km (ug/mL)") # Table 4, row 'Km (ng/mL)' = 0.162 [18.8% RSE, 95% CI 0.102-0.222]; printed unit corrected to ug/mL
    lvmax <- log(0.128)
    label("Maximum Michaelis-Menten elimination rate Vm (mg/day)") # Table 4, row 'Vm (ng/day)' = 0.128 [7.6% RSE, 95% CI 0.109-0.147]; printed unit corrected to mg/day

    # IIV is reported as a percentage for exponential (log-normal) random
    # effects (Methods 2.3.2, Eq 2); variances are written as omega^2 = log(1 + CV^2).
    etalka ~ log(1 + 0.60^2) # Table 4, KA inter-individual variability 60%
    etalvc ~ log(1 + 0.157^2) # Table 4, V/F inter-individual variability 15.7%
    etalcl ~ log(1 + 0.346^2) # Table 4, CL/F inter-individual variability 34.6%

    propSd <- 0.195
    label("Proportional residual error (fraction)") # Table 4, row 'Residual variability' = 19.5%; proportional model (Methods 2.3.2)
  })

  model({
    ka <- exp(lka + etalka)
    vc <- exp(lvc + e_study_mb09c0119_vc * STUDY_MB09C0119 + etalvc)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    km <- exp(lkm)
    vmax <- exp(lvmax)

    # Amounts in mg and volume in L, so Cc is in mg/L (= ug/mL); the
    # saturable arm vmax * Cc / (km + Cc) is an amount rate in mg/day.
    Cc <- central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * Cc - vmax * Cc / (km + Cc)

    Cc ~ prop(propSd)
  })
}
