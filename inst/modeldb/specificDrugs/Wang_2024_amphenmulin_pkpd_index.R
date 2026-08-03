Wang_2024_amphenmulin_pkpd_index <- function() {
  description <- "In vitro (Mycoplasma gallisepticum strain S6). Inhibitory sigmoid Emax PK/PD-index model for the anti-mycoplasma effect of amphenmulin, a novel pleuromutilin derivative, against M. gallisepticum strain S6 in Wang 2024's in-vitro dynamic model. Wang 2024 Materials and Methods (Integration and modeling of pharmacokinetics/pharmacodynamics) parameterises the effect over a 24 h dosing interval as E = E0 - (E0 - Emax) * Ce^N / (EC50^N + Ce^N), where E is the SIGNED change in M. gallisepticum counts over 24 h of cultivation in log10(CFU/mL) (NEGATIVE = bacterial reduction), E0 is the corresponding change in the untreated control, Emax is the maximal (most negative) achievable change, and Ce is a PK/PD index. The packaged model uses the paper's best-correlating index AUC24h/MIC (R = 0.9657, versus 0.8995 for Cmax/MIC), formed as the per-interval covariate AUC_AMPH divided by the parameter mic. Parameters from Wang 2024 Table 3, AUC24h/MIC column: Emax = -2.4214 log10(CFU/mL), E0 = -0.3845 log10(CFU/mL), EC50 = 1199.4720 h, Hill N = 3.1997. Note the sign convention differs from the companion Wang_2024_amphenmulin_killrate model, whose Table 2 parameters are kill RATES with positive meaning killing; here the effect is a signed change with negative meaning killing. Substituting Table 3 into the equation at the paper's headline target of AUC24h/MIC = 904.05 h returns E = -0.97 log10(CFU/mL) against the stated 1 log10 reduction, confirming the sign convention and the equation form to within the rounding of the four-significant-figure parameters. The bacterial density bact (linear CFU/mL) is integrated as d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E across each 24 h interval, reproducing the paper's per-interval model exactly at every 24 h boundary. There is no PK compartment: exposure enters as the externally supplied per-interval AUC_AMPH, because Wang 2024 derived the AUC of the in-vitro dynamic system non-compartmentally in Phoenix. Wang 2024 reports no between-subject variability and no residual error magnitude, so no eta parameters are present and addSd is held at zero for typical-value simulation."
  reference <- paste(
    "Wang W, Yu J, Ji X, Xia X, Ding H. (2024).",
    "Pharmacokinetic/pharmacodynamic integration of amphenmulin:",
    "a novel pleuromutilin derivative against Mycoplasma gallisepticum.",
    "Microbiology Spectrum 12(2):e03675-23.",
    "doi:10.1128/spectrum.03675-23.",
    sep = " "
  )
  vignette <- "Wang_2024_amphenmulin"
  units <- list(
    time = "hour",
    dosing = "h*ug/mL (amphenmulin AUC0-24h per dosing interval, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_AMPH")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Mycoplasma gallisepticum strain S6", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    AUC_AMPH = list(
      description        = "Amphenmulin area under the medium concentration-time curve over the current 24 h dosing interval",
      units              = "h*ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-24-h-interval amphenmulin exposure in the reaction chamber of the in-vitro dynamic model, supplied as a piecewise-constant time-varying covariate (one value per 24 h interval). Wang 2024 derived these non-compartmentally in Phoenix and reports the ranges in the Results section (In vitro pharmacokinetics and the effects on M. gallisepticum): across the 0.1-1.5 ug/mL regimens the AUC was 0.34-6.19 h*ug/mL over 0-24 h, 0.57-7.69 h*ug/mL over 24-48 h and 0.62-8.23 h*ug/mL over 48-72 h. Set to 0 for the untreated control so the sigmoid term vanishes and the predicted 24 h change equals E0. Because the paper published no structural PK model for the apparatus, this covariate is the model's only route for drug exposure.",
      source_name        = "AUC24h (Wang 2024 Materials and Methods, Integration and modeling of pharmacokinetics/pharmacodynamics; ranges in Results, In vitro pharmacokinetics and the effects on M. gallisepticum)"
    )
  )

  population <- list(
    species        = "in vitro (Mycoplasma gallisepticum strain S6 culture)",
    n_subjects     = 3L,
    n_studies      = 1L,
    organism       = "Mycoplasma gallisepticum standard strain S6 (National Center for Veterinary Culture Collection, Beijing); amphenmulin MIC = 0.0039 ug/mL by broth microdilution at a 10^7 CFU/mL inoculum",
    system         = "In-vitro dynamic model: a reservoir chamber of fresh drug-free medium feeding a three-necked reaction flask holding 300 mL of medium in the external compartment and a 10 mL semipermeable-membrane internal compartment carrying the 10^7 CFU/mL culture, draining to a waste chamber under peristaltic-pump control at 1.63 mL/min. The flow rate was set so that amphenmulin elimination matched its 2.13 h intravenous half-life in chickens",
    disease_state  = "Not applicable (in-vitro culture); starting inoculum 10^7 CFU/mL",
    dose_range     = "Peak concentrations of 0.1, 0.2, 0.5, 0.8, 1.0 and 1.5 ug/mL, dosed every 24 h for three doses",
    sampling       = "Drug: reaction-chamber samples at 0.083, 0.167, 0.25, 0.50, 0.75, 1, 2, 3, 4, 6, 9, 12, 24, 24.08, 25, 30, 36, 48.08, 54, 60 and 72 h. Counts: internal-compartment cultures at 0, 3, 6, 9, 12, 24, 30, 36, 48, 54, 60 and 72 h. Each sample in triplicate",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "Wang 2024 Table 3 also reports the same model fitted against Cmax/MIC (Emax = -3.3571, E0 = -0.3575, EC50 = 441.4603, Hill N = 1.3842, R = 0.8995); AUC24h/MIC is the better correlate (R = 0.9657) and is the parameterisation packaged here. The paper's headline targets for a 1 log10 CFU/mL reduction are AUC24h/MIC = 904.05 h and Cmax/MIC = 190.11. A reduction of at least 3 log10 CFU/mL was defined as bactericidal; regimens at or above 0.8 ug/mL reached that level over the three-day course."
  )

  ini({
    # =================================================================
    # Wang 2024 Table 3, AUC24h/MIC column
    # =================================================================
    # Wang 2024 Materials and Methods (Integration and modeling of
    # pharmacokinetics/pharmacodynamics):
    #   E = E0 - (E0 - Emax) * Ce^N / (EC50^N + Ce^N)
    # which is algebraically E0 + (Emax - E0) * Ce^N / (EC50^N + Ce^N),
    # i.e. E runs from E0 at zero exposure to Emax at saturating
    # exposure. E is the SIGNED change in log10(CFU/mL) over a 24 h
    # interval, so both Table 3 values are NEGATIVE and neither can be
    # log-transformed; both stay on the natural scale.
    #
    # Numerical confirmation of the sign convention and equation form:
    # substituting the four values below at the paper's stated target
    # of AUC24h/MIC = 904.05 h returns E = -0.97 log10(CFU/mL) against
    # the stated 1 log10 reduction, and the Cmax/MIC parameterisation
    # at its stated target of 190.11 returns E = -1.07. See the
    # vignette for the reproduction.
    e0 <- -0.3845
    label("Change in mycoplasma count in the untreated control over 24 h E0 (log10 CFU/mL)")  # Wang 2024 Table 3, AUC24h/MIC column, E0 = -0.3845 Log10 CFU/mL
    emax <- -2.4214
    label("Maximal change in mycoplasma count over 24 h Emax (log10 CFU/mL; negative = reduction)")  # Wang 2024 Table 3, AUC24h/MIC column, Emax = -2.4214 Log10 CFU/mL
    lec50 <- log(1199.4720)
    label("AUC24h/MIC producing 50% of the maximal effect EC50 (h)")  # Wang 2024 Table 3, AUC24h/MIC column, EC50 = 1,199.4720
    lhill <- log(3.1997)
    label("Hill coefficient N defining the slope of the fitted curve (unitless)")  # Wang 2024 Table 3, AUC24h/MIC column, Hill's slope = 3.1997

    # =================================================================
    # Minimum inhibitory concentration
    # =================================================================
    # The PK/PD index driving the sigmoid is AUC24h/MIC, so the MIC is
    # required numerically. FIXED because it is a measured property of
    # the challenge strain, not an estimated parameter. Change it to
    # apply the model to an isolate with a different susceptibility.
    mic <- fixed(0.0039)
    label("Amphenmulin MIC against M. gallisepticum S6 (ug/mL; measured by broth microdilution)")  # Wang 2024 Results, Susceptibility assay: MIC = 0.0039 ug/mL by broth microdilution at a 10^7 CFU/mL inoculum

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(7)
    label("log10 starting mycoplasma density in the internal compartment (log10 CFU/mL)")  # Wang 2024 In vitro dynamic model: 10 mL of M. gallisepticum culture at 10^7 CFU/mL

    # =================================================================
    # Residual error
    # =================================================================
    # Wang 2024 reported only the correlation coefficient of the
    # sigmoid fit (R = 0.9657 for AUC24h/MIC; Table 3) and gave no
    # residual standard deviation, so the residual SD is held at zero
    # for deterministic typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Wang 2024)")  # Wang 2024 reports R only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index: amphenmulin AUC over the current 24 h dosing
    # interval divided by the MIC (Wang 2024 Integration and modeling
    # of pharmacokinetics/pharmacodynamics). Units: h.
    aucmic <- AUC_AMPH / mic

    # Wang 2024 inhibitory sigmoid Emax equation. eff_log10 is the
    # SIGNED change in log10(CFU/mL) accrued over the current 24 h
    # interval (negative = bacterial reduction); it equals e0 at zero
    # exposure and approaches emax at saturating exposure.
    eff_log10 <- e0 - (e0 - emax) * aucmic^hill / (ec50^hill + aucmic^hill)

    # Wang 2024 fitted only the per-interval change and counted
    # mycoplasma at 24 h boundaries. Spreading the interval change
    # uniformly across the interval makes log10(bact) change by
    # exactly eff_log10 over each 24 h interval, so the trajectory
    # matches the paper's model at every counted time point:
    #   d(log10 N)/dt = eff_log10 / 24
    # => d(N)/dt      = ln(10) * (eff_log10 / 24) * N
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL.
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
