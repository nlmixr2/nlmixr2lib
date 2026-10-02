Jiang_2021_ivosidenib_QTcF <- function() {
  description <- paste(
    "Linear concentration-QTc model relating the change from baseline in",
    "the Fridericia-corrected QT interval (DeltaQTcF, ms) to the plasma",
    "ivosidenib concentration, pooled across the phase 1 studies",
    "AG120-C-001 (IDH1-mutant hematologic malignancies), AG120-C-002",
    "(IDH1-mutant solid tumors) and AG120-C-004 (healthy volunteers)",
    "(Jiang 2021). Direct effect with no hysteresis:",
    "DeltaQTcF = e0 + slope * CP_IVOSIDENIB_NGML. Typical-value model",
    "only: the paper reports the slope (0.00258 ms per ng/mL) and the",
    "predicted DeltaQTcF at the 500 mg once-daily geometric-mean Cmax",
    "(17.2 ms at 6551 ng/mL), from which the intercept is back-solved;",
    "the between-subject variances of the intercept and slope, the",
    "residual error, and the intercept covariate coefficients are not",
    "reported. PD-only model: the ivosidenib concentration is supplied",
    "as a time-varying covariate, for example from a simulation of",
    "Jiang_2021_ivosidenib."
  )
  reference <- paste(
    "Jiang X, Wada R, Poland B, Kleijn HJ, Fan B, Liu G, Liu H, Kapsalis S,",
    "Yang H, Le K. Population pharmacokinetic and exposure-response analyses",
    "of ivosidenib in patients with IDH1-mutant advanced hematologic",
    "malignancies. Clin Transl Sci. 2021;14(3):942-953.",
    "doi:10.1111/cts.12959"
  )
  vignette <- "Jiang_2021_ivosidenib"
  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external ivosidenib plasma-concentration covariate)",
    concentration = "(observation QTcF is the CHANGE FROM BASELINE in the Fridericia-corrected QT interval, DeltaQTcF, in ms; driving covariate CP_IVOSIDENIB_NGML is in ng/mL)"
  )

  covariateData <- list(
    CP_IVOSIDENIB_NGML = list(
      description = "Instantaneous plasma ivosidenib concentration at the time of each time-matched ECG, supplied as a time-varying covariate",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-matched observed plasma concentrations in the source analysis",
        "(1203 triplicate ECGs with time-matched PK samples from 171",
        "participants; Jiang 2021 Results 'QTc analysis'). No hysteresis",
        "was found, so the concentration enters directly. The slope is",
        "reported in ms per ng/mL, so no in-model unit conversion is",
        "needed. Reference value: the geometric-mean Cmax at 500 mg once",
        "daily is 6551 ng/mL. Set to 0 before the first dose."
      ),
      source_name = "Plasma ivosidenib concentration"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Retained on the intercept of the final C-QTc model ('Increasing age ... associated with increased DeltaQTcF', Jiang 2021 Results), but the coefficient is not reported."
    ),
    QTC_BL = list(
      description = "Baseline QTcF",
      units = "ms",
      type = "continuous",
      notes = paste(
        "Retained on the intercept (lower baseline QTcF associated with",
        "increased DeltaQTcF); coefficient not reported. The final model",
        "also retained lower serum calcium and magnesium and lower use of",
        "medications with known QT-prolongation risk on the intercept",
        "(Jiang 2021 Results 'QTc analysis'), none with a reported",
        "coefficient or coding, so none can be encoded."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 171L,
    n_studies = 3L,
    disease_state = paste(
      "Pooled: patients with IDH1-mutant advanced hematologic malignancies",
      "(AG120-C-001), patients with IDH1-mutant advanced solid tumors",
      "(AG120-C-002), and healthy volunteers (AG120-C-004)"
    ),
    dose_range = "Ivosidenib 100 mg twice daily to 1200 mg once daily orally (AG120-C-001 escalation and 500 mg once-daily expansion); per-study regimens of AG120-C-002 and AG120-C-004 not described in Jiang 2021",
    notes = paste(
      "1203 triplicate ECG measurements with time-matched plasma",
      "concentrations from 171 participants (Jiang 2021 Results 'QTc",
      "analysis'). Demographics of the pooled C-QTc dataset are not",
      "tabulated in the paper."
    )
  )

  ini({
    e0 <- 0.30; label("Intercept of DeltaQTcF at zero concentration (ms)") # Not printed; back-solved as 17.2 - 0.00258 * 6551 = 0.30 ms from the Results 'QTc analysis' prediction of 17.2 ms (90% CI 14.7-19.7) at the 6551 ng/mL geometric-mean Cmax
    slope <- 0.00258; label("Linear slope of DeltaQTcF on plasma ivosidenib concentration (ms per ng/mL)") # Results 'QTc analysis': 'DeltaQTcF was predicted to increase with plasma ivosidenib concentration at 0.00258 msec/(ng/ml)'
    addSd <- fixed(0); label("Additive residual error on DeltaQTcF (ms); not reported, set to 0") # Residual error not reported in Jiang 2021
  })

  model({
    QTcF <- e0 + slope * CP_IVOSIDENIB_NGML
    QTcF ~ add(addSd)
  })
}
