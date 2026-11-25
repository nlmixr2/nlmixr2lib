Lei_2018_tildipirosin <- function() {
  description <- paste0(
    "Veterinary (pig). Ex vivo inhibitory sigmoid Emax PK/PD-integration model for the ",
    "antibacterial effect of tildipirosin against Pasteurella multocida (clinical isolate PM04) ",
    "in serum from pigs dosed 4 mg/kg intramuscularly. The 24 h change in log10 CFU/mL is ",
    "E = e0 - (e0 - emax) * C^N / (C^N + EC50^N), with C the PK/PD index AUC24h/MIC (h), formed ",
    "in the model as the externally supplied per-24-h serum exposure AUC_TILDIPIROSIN divided by the ",
    "serum MIC of PM04 (0.5 ug/mL, a fixed parameter that can be changed to re-target the model ",
    "to another isolate). e0 = +2.54 log10 CFU/mL (drug-free 24 h growth, Lei 2018 Table 4 ",
    "'Emax'), EC50 = 23.52 h and N = 3.48 are from Table 4. The maximal-effect asymptote emax is ",
    "NOT printed in the paper and was recovered by the maintainers by refitting the eight ",
    "observed points digitized from Figure 5 with the Table 4 values held fixed (emax = -5.62); ",
    "back-solving the paper's own bactericidal target (AUC24h/MIC = 29.13 h at E = -3) gives ",
    "-5.63. The bacterial density bact (CFU/mL, starting at 10^6) is integrated as ",
    "d/dt(bact) = ln(10) * (E / 24) * bact so that log10(bact) changes by exactly E over each ",
    "24 h window. There is no PK component: Lei 2018 analysed the pig plasma PK ",
    "non-compartmentally only. No between-subject variability or residual error was reported, ",
    "so there are no etas and addSd is fixed at 0 (typical-value simulation)."
  )
  reference <- paste(
    "Lei Z, Liu Q, Qi Y, Yang B, Khaliq H, Xiong J, Moku GK, Ahmed S, Li K, Zhang H,",
    "Zhang W, Cao J, He Q. Optimal Regimens and Cutoff Evaluation of Tildipirosin Against",
    "Pasteurella multocida. Front Pharmacol. 2018;9:765.",
    "doi:10.3389/fphar.2018.00765.",
    "Sigmoid Emax parameters from Table 4; serum MIC from Table 1; maximal-effect asymptote",
    "digitized from Figure 5.",
    sep = " "
  )
  vignette <- "Lei_2018_tildipirosin"

  units <- list(
    time = "h",
    dosing = "h*ug/mL (tildipirosin serum AUC over a 24 h interval, supplied as a covariate)",
    concentration = "log10 CFU/mL (observation)"
  )

  depends <- c("AUC_TILDIPIROSIN")
  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    bact = list(analyte = "Pasteurella multocida strain PM04", units = "CFU/mL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AUC_TILDIPIROSIN = list(
      description = "Tildipirosin serum area under the concentration-time curve over the current 24 h interval (AUC24h)",
      units = "h*ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Numerator of the AUC24h/MIC index that drives the sigmoid Emax model (Lei 2018 eq 1, ",
        "where C is the AUC24/MIC ratio). Lei 2018 Table 2 reports AUC24h = 3.94 +/- 0.62 ",
        "h*ug/mL after a single 4 mg/kg intramuscular dose and 7.10 +/- 0.91 h*ug/mL after 4 ",
        "mg/kg intravenously; Table 3's AUC24h/MIC of 7.88 h is the intramuscular AUC24h divided ",
        "by the serum MIC of 0.5 ug/mL. Set to 0 for drug-free control records so the sigmoid ",
        "term vanishes and the predicted 24 h change reduces to e0."
      ),
      source_name = "AUC24 h -- Lei 2018 Table 2 (row 'AUC24 h'), Table 3 and Table 4 (index 'AUC24h/MIC')"
    )
  )

  population <- list(
    species = "pig",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "4-5 weeks",
    weight_range = "15-20 kg",
    sex_female_pct = NA_real_,
    disease_state = "healthy pigs (PK and serum donors); ex vivo serum infection with Pasteurella multocida",
    dose_range = "single 4 mg/kg tildipirosin intravenously, then (after an approximately 2-week washout) a single 4 mg/kg intramuscularly",
    regions = "China (Hubei)",
    organism = paste0(
      "Pasteurella multocida PM04, a toxA-positive clinical isolate from Anhui province chosen ",
      "from 112 pig isolates (collected 2013-2016; tildipirosin MIC50 0.5 ug/mL, MIC90 2 ug/mL) ",
      "because its broth MIC (2 ug/mL in TSB) equals the MIC90; MIC 0.5 ug/mL and MBC 1 ug/mL ",
      "in pig serum (Table 1)"
    ),
    notes = paste0(
      "Eight healthy pigs of both sexes (sex split not reported). Plasma sampled at 0.25, 0.5, ",
      "1, 2, 4, 6, 8, 10, 12, 24, 48 and 72 h and assayed by HPLC; PK analysed ",
      "non-compartmentally in Phoenix WinNonlin 6.1. Unbound fraction in plasma 78.00 +/- 5.23%. ",
      "Ex vivo time-kill: serum collected before and 0.25, 0.5, 1, 2, 4, 12 and 24 h after the ",
      "4 mg/kg intramuscular dose was inoculated with PM04 at 10^6 CFU/mL and incubated for ",
      "24 h; the sigmoid Emax model was fitted in WinNonlin."
    )
  )

  ini({
    # =================================================================
    # Lei 2018 Table 4 -- inhibitory sigmoid Emax PK/PD integration
    # =================================================================
    # Lei 2018 eq 1 is printed as E = Emax * (1 - C^N / (C^N + EC50^N)),
    # with "Emax" defined as the log10 change over 24 h in the drug-free
    # control. As printed it can only run from Emax (+2.54) down to 0,
    # yet Table 4 reports index values for E = -3 and E = -4 and the
    # fitted line in Figure 5 descends to about -5.3 at AUC24h/MIC = 50 h.
    # The fit therefore used the standard inhibitory sigmoid form with a
    # negative maximal-effect asymptote (the Aliabadi and Lees form the
    # paper cites):
    #   E = e0 - (e0 - emax) * C^N / (C^N + EC50^N)
    # with e0 = the paper's "Emax" (zero-exposure growth) and emax = the
    # maximal drug effect, which the paper does not print.
    #
    # E is a SIGNED change in log10(CFU/mL), so e0 and emax stay on the
    # natural scale.
    e0 <- 2.54
    label("Change in P. multocida count over 24 h in drug-free serum (log10 CFU/mL)") # Lei 2018 Table 4, 'Emax' = 2.54 +/- 0.37 LgCFU/ml (defined in Methods eq 1 text as the control 24 h change)
    emax <- -5.62
    label("Maximal change in P. multocida count over 24 h (log10 CFU/mL; negative = reduction)") # Not printed in Lei 2018. Digitized by the maintainers: unweighted refit of the 8 observed points in Figure 5 with e0, EC50 and N fixed at Table 4 values gives -5.62; back-solving the Table 4 bactericidal target (29.13 h at E = -3) gives -5.63
    lec50 <- log(23.52)
    label("AUC24h/MIC producing half the maximal effect, EC50 (h)") # Lei 2018 Table 4, EC50 = 23.52 +/- 1.80 h
    lhill <- log(3.48)
    label("Hill coefficient N, steepness of the AUC24h/MIC effect curve (unitless)") # Lei 2018 Table 4, N = 3.48 +/- 0.81

    # =================================================================
    # Minimum inhibitory concentration of the challenge isolate
    # =================================================================
    # Measured property of PM04, not an estimated parameter. The SERUM
    # MIC is used: the ex vivo model was built in serum, and dividing
    # the Table 2 intramuscular AUC24h (3.94) by 0.5 reproduces the
    # Table 3 index (7.88 h) exactly. Change it to apply the model to an
    # isolate of different susceptibility.
    mic <- fixed(0.5)
    label("Tildipirosin MIC against P. multocida PM04 in pig serum (ug/mL)") # Lei 2018 Table 1, 'Serum' row, MIC = 0.5 ug/ml (2 ug/ml in TSB)

    # =================================================================
    # Starting bacterial density
    # =================================================================
    log10_cfu0 <- fixed(6)
    label("log10 starting P. multocida density in the ex vivo serum incubation (log10 CFU/mL)") # Lei 2018 Methods 'Time-Killing Curves in Vitro and ex Vivo': 'The bacteria (10^6 CFU/ml) were co-incubated'

    # =================================================================
    # Residual error
    # =================================================================
    # Lei 2018 reports only R^2 = 0.988 for the fit (Discussion), no
    # residual SD; held at zero for typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Lei 2018)") # Lei 2018 reports R^2 only, no residual SD
  })

  model({
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # PK/PD index AUC24h/MIC (h), Lei 2018 eq 1 driver C.
    ce <- AUC_TILDIPIROSIN / mic

    # Inhibitory sigmoid Emax (see the ini() note on why the printed
    # eq 1 cannot be used literally). eff_log10 is the SIGNED change in
    # log10(CFU/mL) over a 24 h window: e0 at zero exposure, approaching
    # emax at saturating exposure.
    eff_log10 <- e0 - (e0 - emax) * ce^hill / (ce^hill + ec50^hill)

    # Spreading the 24 h change uniformly over the window makes
    # log10(bact) change by exactly eff_log10 across each 24 h, matching
    # the paper's endpoint (bacterial count at 0 and 24 h):
    #   d(log10 N)/dt = eff_log10 / 24  =>  dN/dt = ln(10) * (eff_log10 / 24) * N
    d/dt(bact) <- log(10) * (eff_log10 / 24) * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL.
    log_cfu <- log10(bact + 1)
    log_cfu ~ add(addSd)
  })
}
