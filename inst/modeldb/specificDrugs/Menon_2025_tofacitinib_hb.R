Menon_2025_tofacitinib_hb <- function() {
  description <- paste(
    "Indirect-response model of blood hemoglobin in adults with active",
    "psoriatic arthritis given placebo or tofacitinib 5 or 10 mg twice daily",
    "in OPAL Broaden and OPAL Beyond (Menon 2025). Zero-order production and",
    "first-order loss of hemoglobin, with the loss stimulated linearly by",
    "the individual steady-state average tofacitinib concentration",
    "(EFF = 1 + CAV / C'). CAV is supplied from a population PK model such",
    "as Xie_2019_tofacitinib. Exponential IIV on the baseline and on the",
    "drug slope 1/C'; additive residual error on the log scale."
  )
  reference <- paste(
    "Menon S, Shoji S, Tsuchiwata S, Fallon L, Kanik K. Exposure-Response",
    "Analysis of Tofacitinib in Active Psoriatic Arthritis: Results from Two",
    "Phase 3 Studies. J Clin Pharmacol. 2025;65(3):369-377.",
    "doi:10.1002/jcph.6147. PMID: 39453735. PMCID: PMC11867917."
  )
  vignette <- "Menon_2025_tofacitinib"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "hb (blood hemoglobin, g/dL)"
  )

  compartmentData <- list(
    hb = list(analyte = "hemoglobin", units = "g/dL", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state tofacitinib plasma concentration over the twice-daily dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Patient-specific Cavg predicted from the tofacitinib psoriatic",
        "arthritis population PK model of Xie 2019 (Menon 2025 Methods;",
        "modellib('Xie_2019_tofacitinib')). Applied as a constant from the",
        "first dose, so steady state is assumed from day 0 (tofacitinib",
        "half-life is a few hours against a hemoglobin turnover half-life of",
        "about 50 days). 0 for placebo and for pre-dose records. Median",
        "(range) 17.6 (6.67-47.2) and 36.1 (15.1-88.1) ng/mL at 5 and 10 mg",
        "twice daily (Results)."
      ),
      source_name = "Cavg"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 685L,
    n_studies = 2L,
    age_range = "adults >= 18 years (median 50 years)",
    weight_range = "median 83 kg",
    sex_female_pct = 54.9,
    disease_state = paste(
      "Active psoriatic arthritis (CASPAR criteria, >= 3 swollen and >= 3",
      "tender joints, active plaque psoriasis) with inadequate response to",
      "conventional synthetic DMARDs (OPAL Broaden) or to TNF inhibitors",
      "(OPAL Beyond), on a stable background csDMARD"
    ),
    dose_range = "Placebo (up to Month 3) or tofacitinib 5 or 10 mg orally twice daily (up to Month 6 or Month 12)",
    regions = "Multinational (OPAL Broaden NCT01877668, OPAL Beyond NCT01882439)",
    baseline_hemoglobin = "median 13.8 g/dL",
    notes = paste(
      "Menon 2025 Results: 685 patients (309 male, 376 female), median",
      "weight 83 kg, age 50 years and baseline hemoglobin 13.8 g/dL. Fitted",
      "in NONMEM 7.3."
    )
  )

  ini({
    lrbase <- log(13.8); label("Baseline hemoglobin R0 (g/dL)")  # Table 3, 'R0 (g/dL)' = 13.8 (%RSE 0.343)
    lkout <- log(0.014); label("First-order loss rate constant of hemoglobin (1/day)")  # Table 3, 'k out (day-1)' = 0.014 (%RSE 0.52)
    lcprime <- log(1670); label("Cavg scale C' of the stimulation of hemoglobin loss, EFF = 1 + CAV/C' (ng/mL)")  # Table 3, 'C' (ng/mL)' = 1670 (%RSE 1.71)

    # Table 3 prints the variances to one significant figure (0.008, 0.980)
    # and its footnote gives %CV = 100 * sqrt(omega^2) to three (8.81, 99.0);
    # the variances below are the squares of the footnote values.
    etalrbase ~ 0.0881^2 # Table 3, 'IIV for R0 (variance)' = 0.008; footnote %CV = 8.81
    # The paper puts the eta on the slope 1/C'. An exponential eta on 1/C'
    # is an exponential eta of opposite sign on C', and a zero-mean normal
    # eta is symmetric, so carrying it on lcprime is the same distribution.
    etalcprime ~ 0.99^2 # Table 3, 'IIV for 1/C' (variance)' = 0.980; footnote %CV = 99.0

    # log(Y) = log(R) + theta_residual * eps (Methods): additive on the log
    # scale, which is rxode2's lnorm().
    expSd <- 0.0381; label("Residual SD on the log scale (log g/dL)")  # Table 3, 'RV (CV)' = 3.81% (%RSE 2.11)
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout)
    cprime <- exp(lcprime + etalcprime)
    kin <- rbase * kout

    eff <- 1 + CAV / cprime

    hb(0) <- rbase
    d / dt(hb) <- kin - kout * hb * eff

    hb ~ lnorm(expSd)
  })
}
