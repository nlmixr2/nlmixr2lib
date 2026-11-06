Fouad_2024_cefiderocol <- function() {
  description <- "Two-compartment population PK model for cefiderocol in critically ill adults on continuous venovenous haemodiafiltration (CVVHDF), fitted nonparametrically with the Pmetrics NPAG algorithm and parameterised by clearance, central volume and the intercompartmental micro-rate constants k12 and k21"
  reference <- paste(
    "Fouad A, Kobic E, Nicolasora NP, Thompson Bastin ML, Adams PM, Shen Y,",
    "Fratoni AJ, Ye X, Kuti JL, Nicolau DP, Asempa TE.",
    "Validation of cefiderocol package insert dosing recommendation for",
    "patients receiving continuous renal replacement therapy: a prospective",
    "multicenter pharmacokinetic study.",
    "Open Forum Infect Dis. 2024;11(10):ofae451.",
    "doi:10.1093/ofid/ofae451",
    sep = " "
  )
  vignette <- "Fouad_2024_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The FINAL model is the covariate-free base model (see
  # covariatesDataExcluded for the screen). FU is the one covariate the model
  # consumes, and it acts only on the free-concentration observable that drives
  # the paper's %fT > MIC analysis.
  covariateData <- list(
    FU = list(
      description = paste(
        "Fraction of cefiderocol unbound in plasma. Measured per patient by",
        "ultrafiltration (Centrifree) of a prefilter plasma sample drawn 3 h",
        "after the start of the final infusion, in triplicate, as the",
        "ultrafiltrate concentration divided by the total plasma",
        "concentration. Enters ONLY the free-concentration observable Cfree,",
        "which is what the paper's %fT > MIC analysis is computed on; it does",
        "not scale clearance, volume or any transfer constant."
      ),
      units = "fraction (unitless; 0.84 = 84% unbound)",
      type = "continuous",
      source_name = "free (protein unbound) fraction",
      notes = paste(
        "Cohort mean 0.84 +/- 0.1 (range 0.69-0.98), Results 'Cefiderocol",
        "Exposures'. Per-patient values are not published. The Discussion",
        "notes this is far above the 40-60% binding (FU 0.4-0.6) reported in",
        "non-critically ill patients, which the authors attribute to critical",
        "illness and hypoalbuminaemia (baseline albumin 3.0 +/- 0.6 g/dL)."
      )
    )
  )

  # Covariates that Fouad 2024 screened and did NOT retain. Each was regressed
  # against the individual Bayesian CL and Vc estimates (GraphPad Prism,
  # Supplementary Figure S1). CL correlated significantly with all three, but
  # 'these covariates did not improve AIC from the base model', so the base
  # model was kept as final (Results 'Pharmacokinetics').
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened against CL and Vc (Supplementary Figure S1) and not",
        "retained. Cohort 113.5 +/- 57.3 kg, range 62.4-265 (Table 1). The",
        "Discussion attributes the null AIC result to effluent flow rate",
        "already being prescribed per kg (effluent dose x body weight)."
      )
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "Total CRRT effluent flow rate",
      units = "mL/h",
      type = "continuous",
      source_name = "Effluent flow rate (L/h)",
      notes = paste(
        "Screened against CL and Vc (Supplementary Figure S1) and not",
        "retained. Reported in L/h; multiply by 1000 on ingestion. Cohort",
        "3.34 +/- 0.9 L/h, range 2.08-5.1 (Table 1); effluent dose",
        "32.7 +/- 10 mL/kg/h. The label dose bands (<= 2, 2.1-3, 3.1-4,",
        ">= 4.1 L/h) are applied to this column in the vignette's",
        "simulations, but the model itself does not use it."
      )
    ),
    QEFF = list(
      description = "Transmembrane (CRRT filter) clearance, CLTM",
      units = "L/h",
      type = "continuous",
      source_name = "CLTM",
      notes = paste(
        "Screened against CL and Vc (Supplementary Figure S1) and not",
        "retained. Computed per patient as CLTM = SA * Quf * Qb / (Qb + Qrep)",
        "from the measured saturation coefficient SA (mean 0.8, range",
        "0.34-1.1), ultrafiltrate rate Quf, blood flow Qb and prefilter",
        "replacement fluid rate Qrep (Methods 'Transmembrane Clearance').",
        "Cohort 2.3 +/- 0.6 L/h, range 1-3.5 (Results)."
      )
    )
  )

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = FALSE)
  )

  population <- list(
    species = "human",
    n_subjects = 14,
    n_studies = 1,
    age_range = "33-77 years",
    age_median = "55.4 years (mean, SD 14.9)",
    weight_range = "62.4-265 kg",
    weight_median = "113.5 kg (mean, SD 57.3)",
    bmi = "37.6 +/- 16.2 kg/m^2 (range 21-73)",
    sex_female_pct = 35.7,
    disease_state = paste(
      "critically ill adults in the intensive care unit receiving continuous",
      "renal replacement therapy (all on CVVHDF, Prismaflex), given",
      "cefiderocol for actual or suspected infection (5 of 14) or for PK",
      "sampling only (9 of 14). APACHE II 19 +/- 6.24 (range 10-32)"
    ),
    dose_range = paste(
      "cefiderocol 1.5 or 2 g as a 3 h intravenous infusion: 2 g q6h (n = 5),",
      "2 g q8h (n = 1), 1.5 g q8h (n = 3), 2 g q12h (n = 4), 1.5 g q12h",
      "(n = 1). Five patients were dosed per the label's effluent-flow bands;",
      "nine had adjustments for obesity, residual renal function or",
      "fluctuating effluent flow"
    ),
    renal_function = paste(
      "on CVVHDF; effluent flow 3.34 +/- 0.9 L/h (range 2.08-5.1), effluent",
      "dose 32.7 +/- 10 mL/kg/h. Filters M100 (n = 9), HF1400 (n = 2) and",
      "ST100 (n = 3). Some patients had residual urine output (residual",
      "creatinine clearance 0.5-0.7 L/h in two patients, Supplementary",
      "Table S3)"
    ),
    regions = "United States (3 sites: Hartford CT, Phoenix AZ, Lexington KY)",
    notes = paste(
      "Demographics from Fouad 2024 Table 1. Enrolled August 2022 to",
      "December 2023 (NCT05373615). Prefilter samples at 0, 1.5, 3, 4, 6, 8",
      "and 12 h (12 h only for q12h regimens) after the start of a dose given",
      "after at least 2 prior doses, i.e. at steady state. Fitted in Pmetrics",
      "(NPAG). The two-compartment model was selected on AIC 625.3 against",
      "656.9 for one compartment; fitted gamma 1.6; shrinkage 0.27-5.1%."
    )
  )

  ini({
    # Structural parameters: Fouad 2024 Table 2 'Median' column. NPAG
    # estimates a nonparametric joint density, and Table 2 summarises it by
    # mean +/- SD and median. The median is encoded as the lognormal median
    # exp(mu), and the SD column becomes the omega below. The two columns are
    # mutually consistent under that reading: the lognormal mean
    # median * exp(omega^2 / 2) reproduces the printed mean for CL (3.58 vs
    # 3.5), k12 (3.85 vs 3.9) and k21 (2.12 vs 2.2), and for Vc it gives 9.4
    # against 10.7. Encoding the means as medians instead would overshoot the
    # printed means by exp(omega^2 / 2) (up to 27% for Vc). The vignette
    # compares both columns against the paper's AUCdaily.
    lcl <- log(3.3); label("Total body clearance including CRRT (L/h)") # Table 2: CL median 3.3 (mean 3.5 +/- 1.46)
    lvc <- log(7.4); label("Central volume of distribution (L)") # Table 2: Vc median 7.4 (mean 10.7 +/- 8.4)
    lk12 <- log(3.5); label("Central-to-peripheral transfer rate constant (1/h)") # Table 2: k12 median 3.5 (mean 3.9 +/- 1.8)
    lk21 <- log(1.5); label("Peripheral-to-central transfer rate constant (1/h)") # Table 2: k21 median 1.5 (mean 2.2 +/- 2.2)

    # Interindividual variability. NPAG places every parameter in the joint
    # density, so Table 2's SD column is the between-patient SD on the natural
    # scale. Converted to lognormal variances with omega^2 = log(1 + CV^2),
    # CV = SD / mean:
    #   CL   CV = 1.46 / 3.5 = 41.7%  -> omega^2 = 0.16043
    #   Vc   CV = 8.4 / 10.7 = 78.5%  -> omega^2 = 0.48013
    #   k12  CV = 1.8 / 3.9  = 46.2%  -> omega^2 = 0.19311
    #   k21  CV = 2.2 / 2.2  = 100%   -> omega^2 = 0.69315
    # No covariance matrix is published, so the etas are independent.
    etalcl ~ 0.16043; label("IIV on clearance (log-scale variance)")
    etalvc ~ 0.48013; label("IIV on central volume (log-scale variance)")
    etalk12 ~ 0.19311; label("IIV on k12 (log-scale variance)")
    etalk21 ~ 0.69315; label("IIV on k21 (log-scale variance)")

    # Residual error. Methods 'Pharmacokinetic Analyses' gives the Pmetrics
    # assay polynomial SD = gamma * (0.006 + 0.06 * C), and Results report the
    # fitted gamma as 1.6, so the residual SD is 1.6 * (0.006 + 0.06 * Cc) =
    # 0.0096 + 0.096 * Cc. Pmetrics sums the two parts linearly, which is
    # nlmixr2's combined1() structure. Not fixed(): gamma was estimated.
    addSd <- 0.0096; label("Additive residual SD (mg/L); gamma 1.6 x C0 0.006") # Methods: C0 = 0.006; Results: gamma = 1.6
    propSd <- 0.096; label("Proportional residual SD (fraction); gamma 1.6 x C1 0.06") # Methods: C1 = 0.06; Results: gamma = 1.6
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)

    # Fouad 2024 parameterises the model as CL, Vc, k12 and k21 (Table 2).
    # Q = k12 * Vc and Vp = Q / k21 are algebraic identities of those
    # micro-constants (Q = 25.9 L/h, Vp = 17.3 L at the medians). Defining
    # them lets rxode2's linear-compartment detection recognise a
    # TWO-compartment system; with only cl and vc defined it is detected as
    # one-compartment and the peripheral ODE is dropped.
    q <- k12 * vc
    vp <- q / k21
    kel <- cl / vc

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and vc in L, so central / vc is mg/L.
    Cc <- central / vc

    # Free (unbound) concentration, the quantity the paper's %fT > MIC
    # analysis uses. Only total prefilter plasma cefiderocol was fitted, so the
    # residual error attaches to Cc.
    Cfree <- FU * Cc

    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
