Ojara_2023_paclitaxel_tgi <- function() {
  description <- paste(
    "Tumor growth inhibition (TGI) model for paclitaxel/platinum-treated",
    "advanced non-small cell lung cancer (CEPAC-TDM phase III study, 365",
    "patients, paclitaxel q3w for up to six cycles). Linear (zero-order)",
    "tumor growth plus first-order tumor decay driven by the per-cycle",
    "paclitaxel AUC (AUC_PTX, AUC from start to end of a 21-day cycle after",
    "the single day-1 dose), with an exponential decline of the drug effect",
    "over time (resistance): d/dt(tumor_size) = kgl - kkill *",
    "exp(-res * t) * AUC_PTX * tumor_size, tumor_size(0) = observed",
    "baseline. Includes the exponential parametric time-to-event dropout",
    "model whose hazard falls with tumor shrinkage relative to baseline.",
    "AUC_PTX must be supplied as a time-varying step covariate (0 once",
    "treatment has ended) and t = 0 must be the first paclitaxel dose."
  )
  reference <- paste(
    "Ojara FW, Henrich A, Frances N, Nassar YM, Huisinga W, Hartung N,",
    "Geiger K, Holdenrieder S, Joerger M, Kloft C. A prognostic baseline",
    "blood biomarker and tumor growth kinetics integrated model in",
    "paclitaxel/platinum treated advanced non-small cell lung cancer",
    "patients. CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1714-1725.",
    "doi:10.1002/psp4.12937"
  )
  vignette <- "Ojara_2023_paclitaxel_nsclc"
  units <- list(
    time = "h",
    dosing = "none (paclitaxel exposure enters as the AUC_PTX covariate)",
    concentration = "umol/L (paclitaxel, via AUC_PTX in umol*h/L)"
  )

  covariateData <- list(
    TUMSZ = list(
      description = "Observed baseline tumor size: sum of diameters of the RECIST 1.1 target lesions at cycle 1 day 1. Sets the initial condition of the tumor_size state and is the reference BASE of the dropout hazard covariate (BASE - tumor_size) / BASE.",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Per-subject, time-fixed. The paper reports tumor size in cm (median 8.3 cm across all 365 patients, Results; 8.40 cm, range 1.00-38.3 cm, in the 221-patient biomarker subset, Table S1); multiply by 10 to obtain mm. The paper estimates no IIV on the baseline, so the observed value is used directly.",
      source_name = "BSL"
    ),
    AUC_PTX = list(
      description = "Paclitaxel AUCcycle: area under the total plasma paclitaxel concentration-time curve from the start to the end of a 21-day cycle, from the single paclitaxel dose given on day 1 of that cycle (3-h intravenous infusion). Drives the first-order tumor decay rate kkill * exp(-res * t) * AUC_PTX.",
      units = "umol*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying step covariate: constant within each cycle, updated at each cycle start, and 0 once the last cycle (at most six) has ended. Use last-observation-carried-forward interpolation. Ojara 2023 derived it from the individual paclitaxel PK model of Henrich et al. 2017 (J Pharmacol Exp Ther 362:347-358) for the PK-guided arm and by stochastic simulation (multiple imputation, 50 replicates) for the BSA-guided arm. Median 18.4 umol*h/L (range 3.82-33.3) across cycles in all patients (Results). 1 umol*h/L of paclitaxel = 0.854 mg*h/L (molecular weight 853.9 g/mol).",
      source_name = "AUCcycle"
    )
  )

  compartmentData <- list(
    tumor_size = list(
      analyte = "tumor (sum of target-lesion diameters)",
      units = "mm",
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative hazard of dropout from tumor assessment",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 365L,
    n_studies = 1L,
    age_range = "Adults; 41-78 years (median 64) in the 221-patient biomarker subset (Table S1); the paper does not tabulate the full 365-patient cohort, which is described in Joerger et al. 2016 (Ann Oncol 27:1895-1902).",
    weight_range = "BSA 1.34-2.49 m^2 (median 1.86) in the 221-patient biomarker subset (Table S1).",
    sex_female_pct = 33.5,
    race_ethnicity = "Not reported.",
    disease_state = "Newly diagnosed advanced (stage IIIB/IV) non-small cell lung cancer, first-line paclitaxel plus carboplatin (target AUC 6 mg*min/mL) or cisplatin (80 mg/m^2) every 3 weeks for up to six cycles.",
    dose_range = "Paclitaxel 3-h IV infusion q3w: 200 mg/m^2 in the BSA-guided arm (182 patients, 734 cycles); PK- and neutropenia-guided doses in the PK-guided arm (183 patients, 720 cycles).",
    regions = "Switzerland and Germany (CEPAC-TDM, NCT01326767).",
    notes = "1353 tumor size observations (RECIST 1.1 sum of diameters) at baseline, day 1 of cycles 3 and 5, end of treatment and every 3 months during follow-up; data after week 30 excluded. TGI parameters for both arms combined (Table 1) were estimated by multiple imputation of paclitaxel exposure in the BSA-guided arm. Sex percentage is from the 221-patient biomarker subset (Table S1). The dropout model was fit to both arms combined (Table S4, right column)."
  )

  ini({
    lkgl <- log(1.03e-3); label("Linear (zero-order) tumor growth rate Kg (mm/h)") # Table 1: Kg = 1.03 um/h (RSE 43.2%) = 1.03e-3 mm/h
    lkkill <- log(2.30e-5); label("Paclitaxel AUCcycle-driven tumor decay rate constant at t = 0, beta02_pac (1/(umol*h/L)/h)") # Table 1: beta02_pac = 2.30e-5 (1/(umol/L*h))/h (RSE 12.8%)
    lres <- log(8.75e-4); label("Rate constant for the exponential decline of the drug effect over time, lambda (1/h)") # Table 1: lambda = 8.75e-4 1/h (RSE 15.4%)

    # Table 1 gives IIV as CV% of an exponential model; the variances below
    # are (CV/100)^2. The paper does not state the CV formula; Figure 3c
    # (best-case parameter percentiles) is reproduced with this reading and
    # not with omega^2 = log(1 + CV^2). See the vignette.
    etalkgl ~ 1.3689 # Table 1: IIV Kg CV 117% -> 1.17^2
    etalkkill ~ 0.744769 # Table 1: IIV beta02_pac CV 86.3% -> 0.863^2
    etalres ~ 0.290521 # Table 1: IIV lambda CV 53.9% -> 0.539^2

    expSd <- 0.179; label("Exponential residual error on tumor size (log-scale SD)") # Table 1: residual variability, exponential model, CV 17.9%

    llam_haz <- log(0.000473); label("Baseline hazard of dropout from tumor assessment, lambda0 (1/h)") # Table S4, right column (both arms): lambda0 = 0.000473 1/h (RSE 9.2%)
    e_tsshrink_haz <- -1.88; label("Effect of relative tumor shrinkage (BASE - TS)/BASE on the log dropout hazard (unitless)") # Table S4, right column (both arms): beta = -1.88 (1/RelTS) (RSE 16.9%); Appendix S1 equation 6
  })

  model({
    kgl <- exp(lkgl + etalkgl)
    kkill <- exp(lkkill + etalkkill)
    res <- exp(lres + etalres)

    # Equations 1-4: linear growth minus first-order decay; the decay rate per
    # unit AUCcycle declines exponentially from the first dose (t = 0).
    tumor_size(0) <- TUMSZ
    d/dt(tumor_size) <- kgl - kkill * exp(-res * t) * AUC_PTX * tumor_size

    # Appendix S1 equation 6: exponential dropout hazard driven by tumor
    # shrinkage relative to baseline; shrinkage lowers the hazard.
    lam_haz <- exp(llam_haz)
    ts_shrink <- (TUMSZ - tumor_size) / TUMSZ
    hazard <- lam_haz * exp(e_tsshrink_haz * ts_shrink)
    d/dt(cumhaz) <- hazard
    sur_drop <- exp(-cumhaz)

    tumor_size ~ lnorm(expSd)
  })
}
