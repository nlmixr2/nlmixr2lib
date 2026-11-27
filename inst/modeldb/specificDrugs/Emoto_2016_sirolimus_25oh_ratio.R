Emoto_2016_sirolimus_25oh_ratio <- function() {
  description <- paste(
    "Simple Emax regression of the whole-blood 25-hydroxysirolimus to",
    "sirolimus concentration ratio on postmenstrual age in sirolimus-treated",
    "patients younger than 4 years with complicated vascular anomalies",
    "(Emoto 2016 Figure 3A). The ratio rises from zero at PMA 0 towards",
    "emax with half-maximum at tm50 weeks: mpr_25oh = emax * PMA / (tm50 +",
    "PMA). Fitted by unweighted nonlinear least squares in GraphPad Prism",
    "(no random effects), so the only variability is an additive residual",
    "error. The paper prints only R^2 = 0.34; emax, tm50 and the residual SD",
    "were digitised by the maintainers from the Figure 3A regression line",
    "and 95% prediction band. A CYP3A-activity marker that parallels the",
    "clearance maturation in modellib('Emoto_2016_sirolimus'); no dosing,",
    "volume or concentration-time course."
  )
  reference <- paste(
    "Emoto C, Fukuda T, Mizuno T, Schniedewind B, Christians U, Adams DM,",
    "Vinks AA.",
    "Characterizing the Developmental Trajectory of Sirolimus Clearance in",
    "Neonates and Infants.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(8):411-417.",
    "doi:10.1002/psp4.12096.",
    sep = " "
  )
  vignette <- "Emoto_2016_sirolimus"
  units <- list(
    time = "week (unused by the model; the vignette uses postmenstrual age in weeks as the time axis)",
    dosing = "n/a (no dosing; metabolite-ratio regression)",
    concentration = "unitless (25-hydroxysirolimus / sirolimus whole-blood concentration ratio, observation mpr_25oh)"
  )

  covariateData <- list(
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Supplied in WEEKS, the x-axis of Emoto 2016 Figure 3A. The cohort",
        "had no premature infants, so PMA = postnatal age (weeks) + 40",
        "(Methods). Data in Figure 3A span roughly 55-240 weeks."
      ),
      source_name = "PMA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = "postmenstrual age about 55-240 weeks (patients younger than 4 years)",
    weight_range = "not reported for this subset",
    disease_state = paste(
      "Complicated vascular anomalies in the concentration-controlled",
      "sirolimus phase II trial NCT00975819."
    ),
    dose_range = "Oral sirolimus titrated to a trough of 10-15 ng/mL.",
    regions = "United States (Cincinnati Children's Hospital Medical Center)",
    notes = paste(
      "About 29 metabolite-ratio observations are plotted in Figure 3A; the",
      "number of patients they came from is not reported. Metabolites",
      "assayed by LC/MS-MS (LLOQ 0.05 ng/mL). The 25-OH / sirolimus ratio",
      "correlated with the Bayesian sirolimus clearance estimate (r = 0.66,",
      "P < 0.0001)."
    )
  )

  ini({
    # Emoto 2016 Methods 'Statistical analyses': simple Emax model of the
    # ratio against age, nonlinear regression in GraphPad Prism. No
    # parameter values are printed; values digitised by the maintainers
    # from the Figure 3A regression line (158 points along the curve,
    # least-squares fit of the Emax form; curve passes through the origin).
    lemax <- log(0.184)
    label("Maximum 25-hydroxysirolimus / sirolimus concentration ratio (unitless)") # Figure 3A regression line, digitised: emax = 0.184
    ltm50 <- log(129)
    label("Postmenstrual age at half-maximal ratio (weeks)") # Figure 3A regression line, digitised: tm50 = 129 weeks

    # Additive residual SD from the Figure 3A 95% prediction band: at PMA 0
    # the fitted line and its standard error are both zero, so the band
    # half-width there is t(0.975, n - 2) * s = 0.059 with n = 29 points
    # (dotted band extrapolated to PMA 0 from its first three dots).
    addSd <- 0.029
    label("Additive residual SD of the metabolite / parent ratio (unitless)") # Figure 3A prediction band at PMA 0 = 0.059; s = 0.059 / qt(0.975, 27) = 0.029
  })

  model({
    emax <- exp(lemax)
    tm50 <- exp(ltm50)
    mpr_25oh <- emax * PAGE / (tm50 + PAGE)
    mpr_25oh ~ add(addSd)
  })
}
