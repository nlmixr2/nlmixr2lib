Wen_2016_enrofloxacin_MIC0p01 <- function() {
  description <- "In vitro (Pasteurella multocida bovine-lung isolate, enrofloxacin MIC = 0.01 ug/mL). Additive inhibitory sigmoidal Emax pharmacodynamic model for the rate of growth/decline of the bacterial population under constant exposure to enrofloxacin. Wen 2016 Eq (1) parameterises the rate of change of log10(CFU/mL) as E(C) = E0 - Emax * C^H / (EC50^H + C^H), with C the enrofloxacin concentration in the broth. The packaged model encodes the bacterial density bact (linear CFU/mL) with d/dt(bact) = ln(10) * E(C) * bact so a constant drug exposure yields the linear log10-CFU-vs-time slope that Wen 2016 fit in Phoenix NLME. This is the most-susceptible isolate (concentration-dependent PD per the paper Discussion): Emax = 1.64, EC50 = 0.11 ug/mL (10x MIC), Hill = 0.97, E0 = 0.32 log10(CFU/mL)/h. The enrofloxacin concentration is an external time-varying covariate Cenrofloxacin (no PK component); the experimental design holds it constant for 24 h at 0, 0.5, 0.75, 1, 2, 3, 5, or 10 multiples of the isolate MIC. Random effects are NOT included -- the paper added 'the experiment was added as a random effect' (Methods, Fitting PD model) but did not report the variance, so the packaged model is intended for typical-value / parametric-uncertainty simulation matching the 1,000-Monte-Carlo simulations in Figs 2 and 3. The packaged additive residual SD on log10(CFU/mL) is FIXED at 0 because Wen 2016 did not report a density-scale residual SD (their residual was on the growth-rate regression)."
  reference <- paste(
    "Wen X, Gehring R, Stallbaumer A, Riviere JE, Volkova VV. (2016).",
    "Limitations of MIC as sole metric of pharmacodynamic response across the",
    "range of antimicrobial susceptibilities within a single bacterial species.",
    "Scientific Reports 6:37907.",
    "doi:10.1038/srep37907.",
    sep = " "
  )
  vignette <- "Wen_2016_enrofloxacin_Pmultocida"
  units <- list(
    time = "hour",
    dosing = "ug/mL (enrofloxacin in broth)",
    concentration = "log10 CFU/mL (observation); ug/mL (drug covariate)"
  )

  # Cenrofloxacin is the static enrofloxacin concentration in the broth,
  # an experimentally-controlled in-vitro input (NOT a host PK output).
  # Declared in depends so the canonical-covariate register does not apply
  # (the register is for human pop-PK covariates; this in-vitro PD model
  # documents the input in covariateData for provenance).
  depends <- c("Cenrofloxacin")

  # bact is a paper-mechanistic bacterial-population compartment (linear
  # CFU/mL) and is not in the canonical compartment register.
  paper_specific_compartments <- c("bact")

  covariateData <- list(
    Cenrofloxacin = list(
      description        = "Enrofloxacin concentration in the broth",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Static enrofloxacin concentration in BHI broth held constant for the 24 h time-kill experiment. Wen 2016 Methods (Time-kill experiments paragraph): seven concentrations per isolate normalised as 0 (control), 0.5, 0.75, 1, 2, 3, 5, and 10 multiples of the isolate's MIC. For this isolate (MIC = 0.01 ug/mL), the studied concentrations were 0, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.05, and 0.1 ug/mL. In-vitro experimental input -- not in inst/references/covariate-columns.md.",
      source_name        = "Enrofloxacin concentration (Wen 2016 Methods, Time-kill experiments paragraph)"
    )
  )

  population <- list(
    species             = "in vitro (Pasteurella multocida, bovine-lung isolate)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Pasteurella multocida (bovine-lung isolate from Kansas State Veterinary Diagnostic Laboratory); enrofloxacin MIC = 0.01 ug/mL (broth microdilution per CLSI; confirmed by E-test)",
    system              = "Static time-kill experiments at 37 C in brain heart infusion (BHI) broth with shaking (200 rpm); 20 mL cultures per concentration",
    medium              = "BHI broth (Remel)",
    duration            = "24 h with sampling at 0, 1, 2, 3, 4, 5, 6, 7, 8, 12, and 24 h",
    starting_inoculum   = "5 x 10^5 CFU/mL (1:20 dilution from a 30-min recovery sub-culture)",
    mic                 = "0.01 ug/mL (broth microdilution per CLSI; E-test confirmed)",
    isolate_source      = "Bovine lung",
    regimens            = "Static enrofloxacin concentrations equal to 0 (control), 0.5, 0.75, 1, 2, 3, 5, and 10 multiples of the isolate MIC = 0.01 ug/mL (i.e. 0, 0.005, 0.0075, 0.01, 0.02, 0.03, 0.05, 0.1 ug/mL); two independent experiments per concentration",
    notes               = "Most-susceptible isolate in Wen 2016 (Bovine lung, MIC = 0.01 ug/mL). The PD against this isolate exhibits the expected concentration-dependent characteristics (Discussion): EC50 (0.11 ug/mL) is 10x the MIC, Emax is 1.64, and the inhibitory effect continues to increase at concentrations several multiples of the MIC (Fig. 3a,b). Two further isolates (bovine nasal swabs, MIC = 1.5 and 2.0 ug/mL) are packaged separately as Wen_2016_enrofloxacin_MIC1p5 and Wen_2016_enrofloxacin_MIC2p0."
  )

  ini({
    # =============================================================
    # Wen 2016 Eq (1) parameters (Bovine-lung isolate, MIC = 0.01 ug/mL)
    # =============================================================
    # E(C) = E0 - Emax * C^H / (EC50^H + C^H) is the rate of change
    # (slope) of log10(CFU/mL) under constant enrofloxacin exposure C
    # (Methods, Fitting PD model paragraph). Parameter values come from
    # Wen 2016 Table 1, row "Bovine lung MIC 0.01 ug/mL".
    le0 <- log(0.32)
    label("Bacterial growth rate in absence of drug E0 (log10(CFU/mL) per h)")  # Wen 2016 Table 1, Bovine lung MIC=0.01 ug/mL, E0 = 0.32 (CV 8.86%, CI 0.26-0.38)
    lemax <- log(1.64)
    label("Maximum drug-induced inhibition of growth rate Emax (log10(CFU/mL) per h)")  # Wen 2016 Table 1, Bovine lung, Emax = 1.64 (CV 8.47%, CI 1.34-1.93)
    lec50 <- log(0.11)
    label("Enrofloxacin concentration giving 50% of Emax (EC50; ug/mL)")  # Wen 2016 Table 1, Bovine lung, EC50 = 0.11 (CV 11.49%, CI 0.08-0.13)
    lhill <- log(0.97)
    label("Hill coefficient H (unitless)")  # Wen 2016 Table 1, Bovine lung, H = 0.97 (CV 11.97%, CI 0.73-1.22)

    # =============================================================
    # Initial inoculum
    # =============================================================
    # Wen 2016 Methods, Time-kill experiments paragraph: "The starting
    # bacterial densities in the experiments were 10^5 to 10^6 CFU per
    # mL". The Monte-Carlo simulations in Figs 2 and 3 used a starting
    # density of 5 x 10^5 CFU/mL (Fig. 2 / Fig. 3 captions).
    log10_cfu0 <- log10(5e5)
    label("log10 initial inoculum (log10 CFU/mL; default 5e5 CFU/mL from Fig. 2/3 captions)")  # Wen 2016 Fig. 2 + Fig. 3 captions, "starting density 5 x 10^5 CFU/mL"

    # =============================================================
    # Residual error
    # =============================================================
    # Wen 2016 did NOT report a residual SD on the log10(CFU/mL)
    # observations -- the regression residual was on the rates E(C),
    # whose magnitude is captured by the per-parameter CV% reported
    # in Table 1. The packaged model holds the density-scale residual
    # SD at zero for deterministic typical-value simulation; see the
    # vignette Assumptions and deviations section for details.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (FIXED 0; not reported in Wen 2016)")  # Wen 2016 did not report a density-scale residual SD
  })

  model({
    e0   <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Wen 2016 Eq (1): rate of change of log10(CFU/mL) under constant
    # enrofloxacin exposure Cenrofloxacin.
    growth_rate <- e0 - emax * Cenrofloxacin^hill / (ec50^hill + Cenrofloxacin^hill)

    # Linear-space first-order growth on CFU/mL. d(log10 N)/dt =
    # growth_rate implies d(N)/dt = ln(10) * growth_rate * N at constant
    # drug exposure, so the slope of the log10-density vs. time line is
    # exactly growth_rate.
    d/dt(bact) <- log(10) * growth_rate * bact
    bact(0) <- 10 ^ log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor (matches the
    # Bulitta 2009 / 2010 / Yadav 2017 in-vitro PD convention so the
    # log10 stays finite if bact is driven below 1 CFU/mL).
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
