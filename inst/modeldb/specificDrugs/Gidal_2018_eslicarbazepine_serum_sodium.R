Gidal_2018_eslicarbazepine_serum_sodium <- function() {
  description <- paste0(
    "Linear exposure-response model for the change from baseline in SERUM ",
    "SODIUM concentration as a function of eslicarbazepine exposure in ",
    "adults with focal-onset seizures taking adjunctive eslicarbazepine ",
    "acetate (ESL) (Gidal 2018, 3,354 serum sodium measurements from ",
    "1,128 patients in the phase 3 trials 2093-301, 2093-302 and ",
    "2093-304). The model is Na_ij - BaselineNa_i = ",
    "-0.0000041 * AUC0-24_ij (Gidal 2018 Appendix S1 Eq. E-6, ",
    "Table S-6): reductions in serum sodium are strictly proportional to ",
    "eslicarbazepine exposure, with no intercept, no covariates and no ",
    "saturation. Linear, power and exponential structural forms were ",
    "screened and the linear one was selected. The relationship is weak ",
    "by design and the paper says so: a 400 mg increase in ESL dose ",
    "raises AUC0-24 by about 165 ug*h/mL and is predicted to lower serum ",
    "sodium by only 0.68 mmol/L, roughly one fifteenth of the 10 mmol/L ",
    "decrease the phase 3 trials required before calling a change ",
    "clinically meaningful. Interindividual variability on the slope is ",
    "very large (134.54 %CV), which is how the model accommodates the ",
    "minority of outlier patients who did show clinically significant ",
    "hyponatraemia. There is no PK layer and no ODE: exposure enters as ",
    "the static per-patient column AUC_ESL, an empirical-Bayes ",
    "prediction from modellib('Gidal_2018_eslicarbazepine'). Minimum ",
    "objective function 9,122.097."
  )
  reference <- paste(
    "Gidal BE, Jacobson MP, Ben-Menachem E, Carreno M, Blum D,",
    "Soares-da-Silva P, Falcao A, Rocha F, Moreira J, Grinnell T,",
    "Ludwig E, Fiedler-Kelly J, Passarell J, Sunkaraneni S.",
    "Exposure-safety and efficacy response relationships and population",
    "pharmacokinetics of eslicarbazepine acetate.",
    "Acta Neurol Scand. 2018;138(3):203-211. doi:10.1111/ane.12950.",
    "Parameter table and equation are in Appendix S1 (supporting",
    "information), Table S-6 and Equation E-6.",
    "Exposure metric produced by modellib('Gidal_2018_eslicarbazepine').",
    sep = " "
  )
  vignette <- "Gidal_2018_eslicarbazepine_exposure_response"
  units <- list(
    time = "n/a (exposure-response regression on repeated visit measurements; the model carries no time term)",
    dosing = "n/a (no dose events; exposure enters as the covariate AUC_ESL)",
    concentration = "dsod (change from baseline in serum sodium, mmol/L; the absolute value sod is derived from the SOD baseline)"
  )

  covariateData <- list(
    AUC_ESL = list(
      description = paste(
        "Individual predicted eslicarbazepine area under the plasma",
        "concentration-time curve over the 24-hour dosing interval at the",
        "visit in question."
      ),
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Empirical-Bayes prediction from the population PK model of the",
        "same paper; compute it as dose / CL/F with",
        "modellib('Gidal_2018_eslicarbazepine'). This model is the one",
        "that PINS THE UNITS for the whole Gidal_2018_eslicarbazepine_*",
        "family: Table S-6 prints the slope units explicitly as",
        "(mmol/L)/(ng x h/mL), and the main-text worked example -- a",
        "400 mg dose increase raises AUC0-24 by 165 ug*h/mL and lowers",
        "serum sodium by 0.68 mmol/L -- reproduces exactly on that scale",
        "(165,000 x 0.0000041 = 0.6765). It also independently confirms",
        "the PK model's CL/F: 400 mg / 2.43 L/h = 164.6 ug*h/mL.",
        "Time-varying across visits in principle, since the model is",
        "fitted to repeated measurements, although a patient on a stable",
        "maintenance dose has an essentially constant value. Set to 0 for",
        "placebo patients, which makes the predicted change exactly 0."
      ),
      source_name = "AUC_0-24_ij (Eq. E-6)"
    ),
    SOD = list(
      description = "Serum sodium concentration at baseline.",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-patient baseline value, used ONLY to convert the model's",
        "native change-from-baseline prediction into an absolute serum",
        "sodium concentration (sod = SOD + dsod). It carries no estimated",
        "coefficient: Eq. E-6 is a pure change-from-baseline equation with",
        "no intercept and no baseline-dependent term, so setting SOD does",
        "not alter dsod. Median baseline in the analysis population was",
        "141 mmol/L, range 121-156; the lowest pretreatment value",
        "(121 mmol/L) was in a patient randomised to ESL 400 mg once",
        "daily. mmol/L and mEq/L are numerically identical for the",
        "monovalent sodium ion, and Gidal 2018 uses mEq/L in the main text",
        "and mmol/L in Table S-6."
      ),
      source_name = "BaselineNa_i (Eq. E-6)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1128L,
    n_studies = 3L,
    n_observations = "3,354 serum sodium measurements",
    age_median = "37 years",
    disease_state = "adults with focal-onset seizures, at least four in the 4 weeks before screening despite 1-3 concomitant antiepileptic drugs",
    dose_range = "eslicarbazepine acetate 400, 800 or 1,200 mg orally once daily, or placebo",
    regions = "Western Europe, North America, Latin America and Rest of World",
    baseline_serum_sodium = "median 141 mmol/L, range 121-156 (Gidal 2018 Results)",
    notes = paste(
      "Serum sodium showed no time trend in the placebo group. Some ESL",
      "patients, particularly in the 1,200 mg group, showed a downward",
      "trend; the paper treats these as outliers and reports that levels",
      "stabilised or recovered after about 8 weeks of exposure. A pooled",
      "analysis of the same three trials found the proportion of patients",
      "with plasma sodium at or below 125 mmol/L, and the proportion with",
      "hyponatraemia reported as a treatment-emergent adverse event, to be",
      "dose-related. A visual predictive check using 1,000 simulated",
      "datasets supported the model. Table S-6 is labelled a BASE",
      "STRUCTURAL model rather than a final covariate model -- no",
      "covariate was retained on the slope."
    )
  )

  ini({
    # ==================================================================
    # Gidal 2018 Appendix S1 Table S-6 and Equation E-6:
    #
    #   Na_ij - BaselineNa_i = -0.0000041 * AUC_0-24_ij
    #
    # Two estimated rows only: the slope and the residual SD. There is no
    # intercept term, so a patient with zero exposure is predicted to have
    # no change from baseline. No covariate was retained.
    # ==================================================================

    # ----- Exposure slope -----
    e_auc_esl_dsod <- -0.0000041 ; label("Change in serum sodium per 1 ng*h/mL of eslicarbazepine AUC0-24 ((mmol/L)/(ng*h/mL))")  # Table S-6, 'Slope [(mmol/L)/(ng x h/mL)]' -0.0000041, 8.7% SEM; Eq. E-6. Check: 165,000 ng*h/mL x 0.0000041 = 0.68 mmol/L, the main-text worked example for a 400 mg dose increase

    # ----- Interindividual variability on the slope -----
    # Table S-6 reports IIV as 134.54 %CV; the variance on the log scale
    # is omega^2 = log(CV^2 + 1) = log(1.3454^2 + 1) = 1.0666. The IIV is
    # carried on the log of the slope MAGNITUDE so that the sign of the
    # effect is preserved for every simulated subject, which matches the
    # source's exponential-IIV parameterisation and keeps the direction
    # of the sodium change physiologically consistent.
    etae_auc_esl_dsod ~ 1.0666  # Table S-6, slope IIV 134.54 %CV, 8.9% SEM: log(1.3454^2 + 1) = 1.0666. Very large, which is how the model accommodates the outlier patients described in the Results

    # ----- Residual error -----
    addSd_dsod <- 2.08 ; label("Additive residual SD on serum sodium (mmol/L)")  # Table S-6, 'Residual error (SD in mmol/L)' 2.08, 5.2% SEM
  })

  model({
    # ----- Individual slope (Gidal 2018 Eq. E-6) -----
    # The slope is negative, so the exponential IIV is applied to its
    # magnitude and the sign is reattached; every subject therefore has a
    # non-positive slope, as the source's structural model requires.
    slope_dsod <- e_auc_esl_dsod * exp(etae_auc_esl_dsod)

    # ----- Predicted change from baseline -----
    dsod <- slope_dsod * AUC_ESL

    # ----- Absolute serum sodium, for convenience -----
    # Eq. E-6 is a change-from-baseline equation; adding the measured
    # baseline recovers the concentration scale of the observed data.
    sod <- SOD + dsod

    # ----- Observation -----
    dsod ~ add(addSd_dsod)
  })
}
