Yang_2026_APTM_timekill <- function() {
  description <- "In vitro (Mycoplasma gallisepticum strain S6, ATCC 15302). Semi-mechanistic time-kill pharmacodynamic model for APTM (14-O-[(4-amino-6-hydroxy-pyrimidine-2-yl) thioacetyl] mutilin), a novel semi-synthetic pleuromutilin derivative, against M. gallisepticum. Yang 2026 Materials and methods (Pharmacokinetic, pharmacodynamic, and statistical analysis) fits the static time-kill curves to the printed ODE dN/dt = kgrowth * N - (Emax * C^gamma / (EC50^gamma + C^gamma)) * N, i.e. a net first-order growth rate reduced by a sigmoid Emax kill-rate term acting on the same state; both kgrowth and Emax are first-order rate constants in 1/h and the drug effect is a rate constant, not a log10-per-hour kill rate (contrast the sibling Wang_2024_amphenmulin_killrate, whose printed equation is a rate on the log10 scale). Parameters from Yang 2026 Table 1: kgrowth = 0.0578 1/h, Emax = 0.2932 1/h, EC50 = 0.0766 ug/mL, Hill gamma = 4.12. APTM exposure is static: the time-kill assay used a single drug addition to a sealed macrodilution tube with no medium exchange and Yang 2026 reports no degradation rate, so the `aptm` state holds the bath concentration in ug/mL and is integrated as d/dt(aptm) = 0, following the in-vitro convention of HernandezLozano_2025_apramycin_invitro. The printed equation has no carrying-capacity term, so the drug-free control grows without bound; over the 48 h horizon of the source experiment this is immaterial and reproduces the paper's own descriptive statement (kgrowth * 48 / ln(10) = 1.20 log10 CFU/mL against the reported growth of 'approximately 1 log10CFU/mL over 48 hours'). Yang 2026 reports neither between-subject variability nor a residual error magnitude for this fit, so no eta parameters are present and addSd is FIXED at 0 for deterministic typical-value simulation. Sibling models from the same paper: Yang_2026_APTM_aucmic and Yang_2026_APTM_cmaxmic (the in vivo chicken inhibitory sigmoid Emax PK/PD-index fits)."
  reference <- paste(
    "Yang W, Ding H, Ma X, Lv T, Wang L. (2026).",
    "Pharmacokinetic/pharmacodynamic relationship of a novel pleuromutilin",
    "derivative APTM against Mycoplasma gallisepticum.",
    "Poultry Science 105:106560.",
    "doi:10.1016/j.psj.2026.106560. PMCID: PMC12919259.",
    "Model equation: Materials and methods, 'Pharmacokinetic, pharmacodynamic,",
    "and statistical analysis'. Parameter estimates: Table 1.",
    "MIC values and the descriptive time-kill behaviour: Results,",
    "'In vitro susceptibility of MG to APTM' and 'In Vitro time-kill kinetics'.",
    sep = " "
  )
  vignette <- "Yang_2026_APTM"
  units <- list(
    time = "h",
    dosing = "ug/mL (APTM bath concentration placed directly into the static `aptm` state at time 0)",
    concentration = "log10 CFU/mL (observation)"
  )

  # The dosing target is neither `depot` nor `central`: this is a static in
  # vitro system, so the "dose" is the bath concentration placed directly into
  # the `aptm` state at time 0.
  dosing <- c("aptm")

  paper_specific_compartments <- c("aptm", "bact")

  compartmentData <- list(
    aptm = list(analyte = "APTM (14-O-[(4-amino-6-hydroxy-pyrimidine-2-yl) thioacetyl] mutilin)", units = "ug/mL", specimen = "administration site", verified = TRUE),
    bact = list(analyte = "Mycoplasma gallisepticum strain S6 (ATCC 15302)", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "in vitro (Mycoplasma gallisepticum strain S6 culture)",
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    organism       = "Mycoplasma gallisepticum standard strain S6 (ATCC 15302; China Institute of Veterinary Drug Control), cultured at 37 C in a humidified 5% CO2 atmosphere in M. gallisepticum basal medium supplemented with 10% porcine serum, 2% penicillin, 0.013% reduced NADH and L-cysteine. APTM MIC = 0.03125 ug/mL by broth microdilution and 0.0625 ug/mL by broth macrodilution, constant across initial inoculum densities of 10^5, 10^6 and 10^7 CFU/mL",
    system         = "Static time-kill assay by the broth macrodilution method in 2 mL tubes, incubated at 37 C, with viable counts by the drop plate method after ten-fold serial dilution",
    medium         = "M. gallisepticum basal medium (Qingdao Hope Biological Technology) with 10% porcine serum, 2% penicillin, 0.013% NADH and L-cysteine",
    temperature    = "37 C",
    duration       = "48 h",
    starting_inoculum = "approximately 10^5 CFU/mL",
    limit_of_detection = "10 CFU/mL (drop plate method)",
    concentration_range = "0 (drug-free growth control) to 16 x MIC. Because the time-kill assay used the macrodilution method, the relevant MIC for converting the x MIC multiples to absolute concentrations is the macrodilution value 0.0625 ug/mL, giving 0.0625 to 1 ug/mL",
    disease_state  = "not applicable (in vitro)",
    sampling       = "Viable counts at 0, 4, 8, 12, 24, 36 and 48 h; each experiment performed in triplicate",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "Ethics approval 2025C037 (Animal Ethics Committee of South China Agricultural University). APTM was supplied by Shandong Qilu KingPhar Pharmaceutical Co., Ltd.; the reference standard (96.5% purity) was used for the in vitro work. Yang 2026 Results reports that concentrations at and above 2 x MIC reduced the count to the 10 CFU/mL limit of detection within 24 h, but the Table 1 parameter set cannot reproduce that: the maximum attainable net rate is kgrowth - Emax = -0.2354 1/h, which caps the 24 h decline at 2.45 log10 CFU/mL, so a 10^5 CFU/mL inoculum cannot reach 10 CFU/mL (a 4 log10 drop) in 24 h at any concentration. The fitted Emax therefore under-predicts the maximal killing that the raw curves display. This is a property of the published fit, not of the encoding; see the vignette Errata."
  )

  ini({
    # =================================================================
    # Yang 2026 Table 1 -- in vitro time-kill model
    # =================================================================
    # Yang 2026 Materials and methods, "Pharmacokinetic, pharmacodynamic,
    # and statistical analysis", printed equation:
    #
    #   dN/dt = kgrowth * N - (Emax * C^gamma / (EC50^gamma + C^gamma)) * N
    #
    # kgrowth and Emax are both first-order rate constants (1/h), so the
    # kill term is a rate constant multiplying N rather than a log10-scale
    # decline rate. All four parameters are positive and are carried on
    # the natural-log scale for the positivity constraint.
    lkgrow <- log(0.0578)
    label("Log net growth rate constant of M. gallisepticum in the absence of drug kgrowth (1/h)")  # Yang 2026 Table 1, K growth = 0.0578 hr-1
    lemax <- log(0.2932)
    label("Log maximum killing rate constant produced by APTM Emax (1/h)")  # Yang 2026 Table 1, E max = 0.2932 hr-1
    lec50 <- log(0.0766)
    label("Log APTM concentration producing 50% of the maximum killing effect EC50 (ug/mL)")  # Yang 2026 Table 1, EC 50 = 0.0766 ug/mL
    lhill <- log(4.12)
    label("Log Hill coefficient gamma defining the steepness of the concentration-effect curve (unitless)")  # Yang 2026 Table 1, gamma = 4.12

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(5)
    label("log10 starting mycoplasma density in the time-kill tube (log10 CFU/mL)")  # Yang 2026 Materials and methods, "In vitro susceptibility testing and time-kill curve experiment": initial inoculum of approximately 10^5 CFU/mL

    # =================================================================
    # Residual error
    # =================================================================
    # Yang 2026 Table 1 reports the four point estimates only -- no
    # standard errors, no RSEs and no residual standard deviation -- so
    # the residual SD is held at zero for deterministic typical-value
    # simulation. See the vignette Assumptions and deviations section.
    addSd <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (0; not reported in Yang 2026)")  # Yang 2026 Table 1 reports point estimates only, with no uncertainty and no residual SD
  })

  model({
    kgrow <- exp(lkgrow)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Static in vitro exposure. The time-kill assay was a single drug
    # addition to a sealed macrodilution tube with no medium exchange, and
    # Yang 2026 reports no APTM degradation rate, so the concentration is
    # held constant. A bolus of `amt` into `aptm` sets the bath
    # concentration directly in ug/mL. Give this state a first-order
    # elimination to simulate a dynamic-exposure variant of the experiment.
    d/dt(aptm) <- 0

    # Sigmoid Emax kill-rate constant (1/h) at the current bath
    # concentration, from the printed equation.
    kkill <- emax * aptm^hill / (ec50^hill + aptm^hill)

    # Yang 2026 printed ODE, encoded verbatim. There is no
    # carrying-capacity term in the source equation, so the drug-free
    # control grows exponentially at kgrow.
    d/dt(bact) <- kgrow * bact - kkill * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL (the in-vitro PD convention
    # used by Chen_2023_tilmicosin and Wang_2024_amphenmulin_killrate).
    Cc <- log10(bact + 1)
    Cc ~ add(addSd)
  })
}
