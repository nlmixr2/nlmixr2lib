Park_2023_pregabalin <- function() {
  description <- "One-compartment population PK model with first-order absorption, an absorption lag time and first-order elimination for a once-daily sustained-release pregabalin tablet (YHD1119, 75 and 150 mg) in Korean adults with near-normal renal function or moderate renal impairment (Park 2023). Apparent clearance is proportional to time-varying MDRD eGFR up to an estimated breakpoint of 91.2 mL/min/1.73 m^2, above which it is flat. Between-subject variability on CL/F and V/F, inter-occasion variability on the lag time, and a proportional residual error."
  reference <- "Park M, Choi S, Han S, Shin W, Kim A, Han S, Kim B, Lim Y, Yoo H. Pharmacokinetic properties of a new sustained-release pregabalin tablet in subjects with reduced renal function. Transl Clin Pharmacol. 2023;31(4):226-237. doi:10.12793/tcp.2023.31.e20"
  vignette <- "Park_2023_pregabalin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the Modification of Diet in Renal Disease (MDRD) equation, BSA-normalized (mL/min/1.73 m^2). Enters CL/F as the through-origin linear ratio min(CRCL, crcl_hinge) / 73.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Park 2023 Results 'Population PK analysis': 'eGFR values were",
        "incorporated as time-variant covariates of CL/F in a linear",
        "relationship centered on the weighted mean of eGFR up to the",
        "breakpoint.' The weighted mean itself is not reported; the",
        "reference value 73 was back-solved by the maintainers from the",
        "steady-state curves of Figure 4 (see the vignette). Cohort medians",
        "91.5 (79.0-117.0) for near-normal renal function and 45.5",
        "(36.0-59.0) for moderate renal impairment (Table 1, screening",
        "values)."
      ),
      source_name = "eGFR"
    ),
    OCC = list(
      description = "Integer-valued occasion (study period) indicator for inter-occasion variability on the absorption lag time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Park 2023 Methods: 'Inter-occasional variability (IOV) was allowed",
        "only for absorption parameters'; the final model carries IOV on",
        "Tlag only (Results). The near-normal renal function cohort had two",
        "periods (OCC = 1: 75 mg tablet; OCC = 2: 150 mg tablet after a",
        "four-day washout); the moderate renal impairment cohort had one",
        "(OCC = 1). Decomposed inside model() into oc1 / oc2 indicators",
        "multiplexing the IOV etas (Chen_2023_nemonoxacin.R precedent). For",
        "single-occasion simulation pass OCC = 1."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "pregabalin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "pregabalin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "21-67 years",
    age_median = "38 years (near-normal renal function) and 60 years (moderate renal impairment)",
    weight_range = "46.6-97.0 kg",
    weight_median = "70.6 kg (near-normal renal function) and 73.5 kg (moderate renal impairment)",
    sex_female_pct = 6.25,
    race_ethnicity = c(Asian = 100),
    disease_state = "Adults with near-normal renal function (MDRD eGFR > 60 mL/min/1.73 m^2, n = 8) or moderate renal impairment (30-60 mL/min/1.73 m^2, n = 8); no other clinically significant disease.",
    dose_range = "Single oral dose of the YHD1119 sustained-release tablet after a standard meal: 75 mg then 150 mg (four-day washout) in the near-normal cohort; 75 mg only in the moderate renal impairment cohort.",
    regions = "South Korea (Seoul St. Mary's Hospital and CHA Bundang Medical Center)",
    renal_function = "MDRD eGFR median 91.5 (79.0-117.0) and 45.5 (36.0-59.0) mL/min/1.73 m^2 by cohort (Table 1)",
    notes = "Two-centre open-label study NCT05012436. Plasma sampled pre-dose and 1-72 h post-dose (14 samples per period); LC-MS/MS, LLOQ 20 ng/mL, BQL data handled with the M4 method. Fitted with SAEM in MonolixSuite 2023R1. Demographics from Table 1."
  )

  ini({
    # Table 4 fixed effects (Monolix final estimates, all estimated).
    ltlag <- log(0.889); label("Absorption lag time Tlag (h)") # Table 4 'Tlag, h' = 0.889 (RSE 16.2%)
    lka <- log(0.347); label("First-order absorption rate constant ka (1/h)") # Table 4 'ka, h-1' = 0.347 (RSE 8.02%)
    lvc <- log(61.7); label("Apparent volume of distribution V/F (L)") # Table 4 'Vd/F, L' = 61.7 (RSE 5.24%)
    lcl <- log(5.55); label("Apparent clearance CL/F at the reference eGFR of 73 mL/min/1.73 m^2 (L/h)") # Table 4 'CL/F, L/h' = 5.55 (RSE 6.67%)
    lcrcl_hinge <- log(91.2); label("eGFR breakpoint above which CL/F no longer rises (mL/min/1.73 m^2)") # Table 4 'Breakpoint' = 91.2 (RSE 4.38%); Results text gives 92.10

    # Monolix reports omega as the SD of the log-normal eta (Table 4 CV% =
    # sqrt(exp(omega^2) - 1)), so the nlmixr2 variances are omega^2.
    etalcl ~ 0.058564 # Table 4 BSV 'CL/F' omega = 0.242 (CV 24.6%)
    etalvc ~ 0.021609 # Table 4 BSV 'Vd/F' omega = 0.147 (CV 14.8%)

    # Inter-occasion variability on Tlag, one eta per period with the variance
    # shared across occasions.
    etaiov_tlag_1 ~ 0.527076 # Table 4 IOV 'Tlag' omega = 0.726 (CV 83.3%)
    etaiov_tlag_2 ~ fix(0.527076) # same IOV variance as occasion 1

    propSd <- 0.286; label("Proportional residual error (fraction)") # Table 4 'Proportional error' = 0.286 (RSE 5.16%)
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_tlag <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2

    # CL/F proportional to eGFR up to the breakpoint, flat above it. The
    # centring value 73 mL/min/1.73 m^2 stands in for the unreported
    # 'weighted mean of eGFR up to the breakpoint'; it is back-solved from the
    # Figure 4 steady-state median curves (72.3-73.7 across the three curves).
    crcl_hinge <- exp(lcrcl_hinge)
    cl <- exp(lcl + etalcl) * min(CRCL, crcl_hinge) / 73
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)
    tlag <- exp(ltlag + iov_tlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # central (mg) / vc (L) = mg/L; x 1000 gives ng/mL, the assay units.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
