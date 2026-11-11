Gumbo_2026_ethambutol_hfs <- function() {
  description <- paste(
    "In vitro (hollow fiber system model of Mycobacterium kansasii, HFS-MKC,",
    "intracellular THP-1 infection, reference strain ATCC 12478).",
    "Exposure-response model for ethambutol monotherapy (Gumbo 2026): on each",
    "sampling day (4, 7, 14 and 21) the bacterial burden is an inhibitory",
    "sigmoid Emax function of the AUC0-24/MIC exposure, fitted separately for",
    "the agar log10 CFU/mL readout and for the MGIT time-to-positivity readout",
    "(Table 2), and the ethambutol-resistant subpopulation (log10 CFU/mL on",
    "agar with 3x MIC ethambutol) is a quadratic function of AUC0-24/MIC (the",
    "U-shaped curve of the antibiotic resistance arrow of time model). The",
    "model is algebraic: it has no ODE, the day-specific parameter set is",
    "selected from the simulation time in days, and predictions are defined",
    "only at the four sampling days. Exposure enters as the absolute daily AUC",
    "covariate AUC_ETHAMBUTOL divided by the strain MIC (8 mg/L). The",
    "resistant-subpopulation coefficients were digitised from Figure 4C",
    "because the paper does not tabulate them; like every fit here they hold",
    "only over the tested exposure range, AUC0-24/MIC 0 to 2.07, and the",
    "quadratic rises without bound beyond it. No residual error was",
    "reported, so all residual SDs are fixed at 0."
  )
  reference <- paste(
    "Gumbo T, Boorgula GD, Srivastava S.",
    "Preclinical pharmacokinetics/pharmacodynamics studies defining the role of",
    "ethambutol in Mycobacterium kansasii lung disease. Antimicrob Agents",
    "Chemother. 2026;70(2):e01447-25. doi:10.1128/aac.01447-25.",
    "Sigmoid Emax parameters from Table 2; resistant-subpopulation quadratics",
    "digitised from Figure 4C; strain MIC from Results."
  )
  vignette <- "Gumbo_2026_ethambutol"
  units <- list(
    time = "day",
    dosing = "mg*h/L (ethambutol AUC0-24, supplied as a covariate)",
    concentration = "log10 CFU/mL (log_cfu, log_cfu_res); days (ttp)"
  )

  depends <- c("AUC_ETHAMBUTOL")

  # ttp (MGIT time-to-positivity) and log_cfu_res (ethambutol-resistant
  # subpopulation) are algebraic observation outputs of this paper.
  paper_specific_compartments <- c("ttp", "log_cfu_res")

  covariateData <- list(
    AUC_ETHAMBUTOL = list(
      description = "Ethambutol area under the concentration-time curve over one 24 h dosing interval (AUC0-24) in the hollow-fiber central compartment",
      units = "mg*h/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Gumbo 2026 expresses exposure as AUC0-24/MIC. The ratio is split",
        "into the absolute AUC0-24 and the model parameter mic = 8 mg/L",
        "(the ATCC 12478 MIC) because that MIC reproduces the paper's own",
        "exposures: the HFS-MKC AUC0-24 values plotted in Figure 2B (about",
        "1.0, 2.2, 4, 11.3, 14 and 16.6 mg*h/L for regimens R1-R6) equal 8 x",
        "the AUC0-24/MIC values printed in the Figure 3 panel titles (0.13,",
        "0.28, 0.5, 1.41, 1.74, 2.07). Set to 0 for non-treated controls."
      ),
      source_name = "AUC0-24 (AUC0-24/MIC x MIC); Figure 2B, Figure 3 panel titles, Table 2"
    )
  )

  population <- list(
    species = "in vitro (hollow fiber system, M. kansasii ATCC 12478 in THP-1 monocytes)",
    n_subjects = 8L,
    n_studies = 1L,
    organism = "Mycobacterium kansasii complex reference strain ATCC 12478; ethambutol MIC 8 mg/L (broth microdilution); mutation frequency on 3x MIC agar 4.23 +/- 1.54 x 10^-2 (4.23% prevalence)",
    system = paste(
      "Eight HFS-MKC units (cellulosic cartridges, FiberCell Systems): six",
      "ethambutol regimens (R1-R6) and two non-treated controls. THP-1 cells",
      "infected for 4 h at a multiplicity of infection of 1:1 were washed and",
      "20 mL inoculated into the peripheral compartment (2.73 x 10^6",
      "cells/mL); circulating medium RPMI-1640 with 2% FBS. Ethambutol was",
      "infused once daily over 1 h, with the dilution rate set to mimic the",
      "unbound concentration-time profile at the centre of tuberculosis lung",
      "cavities. Central-compartment PK samples pre-dose and 2, 5, 8, 12, 18",
      "and 23.5 h post-infusion; peripheral-compartment bacterial burden",
      "(agar CFU, agar with 24 mg/L ethambutol, MGIT time-to-positivity) on",
      "days 4, 7, 14 and 21."
    ),
    disease_state = "in vitro intracellular and extracellular M. kansasii infection model",
    dose_range = "Six once-daily ethambutol exposures, AUC0-24/MIC 0.13, 0.28, 0.5, 1.41, 1.74 and 2.07 (Cmax/MIC 0.01-0.21), for 21 days",
    regions = "USA (University of Texas at Tyler)",
    notes = paste(
      "Day-0 bacterial burden (B0) was about 5.3 log10 CFU/mL (Figure 3),",
      "with time-to-positivity about 3.2 days (Figure 3A); no exposure killed",
      "below B0. THP-1 cells survived to day 7, so the model represents both",
      "intracellular and extracellular bacilli. AUC0-24/MIC had a better",
      "corrected AIC than Cmax/MIC or %TMIC on every sampling day. The EC80",
      "averaged over the eight day-by-readout fits, AUC0-24/MIC 5.47 (95% CI",
      "1.17-9.77), is the target exposure used in the Monte Carlo",
      "experiments (see Gumbo_2026_ethambutol)."
    )
  )

  ini({
    # Ethambutol MIC of the HFS-MKC strain; design constant, not estimated.
    mic <- fixed(8)
    label("Ethambutol MIC of M. kansasii ATCC 12478 (mg/L)") # Gumbo 2026 Results, 'The ethambutol MIC of the MKC standard laboratory isolate (American Type Culture Collection ATCC#12478) was 8 mg/L'

    # ----- Gumbo 2026 Table 2, CFU/mL readout --------------------------
    # log_cfu = e0 - emax * X^hill / (X^hill + ec50^hill), X = AUC0-24/MIC
    le0_cfu_d4 <- log(7.45)
    label("Day 4 drug-free bacterial burden Econ (log10 CFU/mL)") # Table 2, CFU/mL readout, Econ, Day 4 = 7.45 +/- 0.22
    lemax_cfu_d4 <- log(0.62)
    label("Day 4 maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Table 2, CFU/mL readout, Emax, Day 4 = 0.62 +/- 1.09
    lhill_cfu_d4 <- log(1.60)
    label("Day 4 Hill coefficient H, CFU readout (unitless)") # Table 2, CFU/mL readout, H, Day 4 = 1.60 +/- 4.36
    lec50_cfu_d4 <- log(3.84)
    label("Day 4 AUC0-24/MIC at half-maximal effect EC50, CFU readout (h)") # Table 2, CFU/mL readout, EC50, Day 4 = 3.84 +/- 11.42

    le0_cfu_d7 <- log(7.62)
    label("Day 7 drug-free bacterial burden Econ (log10 CFU/mL)") # Table 2, CFU/mL readout, Econ, Day 7 = 7.62 +/- 0.04
    lemax_cfu_d7 <- log(0.64)
    label("Day 7 maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Table 2, CFU/mL readout, Emax, Day 7 = 0.64 +/- 0.09
    lhill_cfu_d7 <- log(40.22)
    label("Day 7 Hill coefficient H, CFU readout (unitless)") # Table 2, CFU/mL readout, H, Day 7 = 40.22 +/- 231.9
    lec50_cfu_d7 <- log(1.81)
    label("Day 7 AUC0-24/MIC at half-maximal effect EC50, CFU readout (h)") # Table 2, CFU/mL readout, EC50, Day 7 = 1.81 +/- 0.43

    le0_cfu_d14 <- log(7.71)
    label("Day 14 drug-free bacterial burden Econ (log10 CFU/mL)") # Table 2, CFU/mL readout, Econ, Day 14 = 7.71 +/- 0.07
    lemax_cfu_d14 <- log(1.42)
    label("Day 14 maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Table 2, CFU/mL readout, Emax, Day 14 = 1.42 +/- 1.55
    lhill_cfu_d14 <- log(1.0)
    label("Day 14 Hill coefficient H, CFU readout (unitless)") # Table 2, CFU/mL readout, H, Day 14 = 1.0 +/- 1.5
    lec50_cfu_d14 <- log(0.97)
    label("Day 14 AUC0-24/MIC at half-maximal effect EC50, CFU readout (h)") # Table 2, CFU/mL readout, EC50, Day 14 = 0.97 +/- 2.05

    le0_cfu_d21 <- log(8.14)
    label("Day 21 drug-free bacterial burden Econ (log10 CFU/mL)") # Table 2, CFU/mL readout, Econ, Day 21 = 8.14 +/- 0.11
    lemax_cfu_d21 <- log(1.81)
    label("Day 21 maximal reduction in bacterial burden Emax (log10 CFU/mL)") # Table 2, CFU/mL readout, Emax, Day 21 = 1.81 +/- 0.17
    lhill_cfu_d21 <- log(0.98)
    label("Day 21 Hill coefficient H, CFU readout (unitless)") # Table 2, CFU/mL readout, H, Day 21 = 0.98 +/- 0.49
    lec50_cfu_d21 <- log(0.09)
    label("Day 21 AUC0-24/MIC at half-maximal effect EC50, CFU readout (h)") # Table 2, CFU/mL readout, EC50, Day 21 = 0.09 +/- 0.03

    # ----- Gumbo 2026 Table 2, time-to-positivity readout ---------------
    # Longer time-to-positivity means fewer bacilli, so drug effect ADDS:
    # ttp = e0 + emax * X^hill / (X^hill + ec50^hill)
    le0_ttp_d4 <- log(1.23)
    label("Day 4 drug-free time-to-positivity Econ (day)") # Table 2, time-to-positivity readout, Econ, Day 4 = 1.23 +/- 0.02
    lemax_ttp_d4 <- log(0.58)
    label("Day 4 maximal increase in time-to-positivity Emax (day)") # Table 2, time-to-positivity readout, Emax, Day 4 = 0.58 +/- 0.03
    lhill_ttp_d4 <- log(2.34)
    label("Day 4 Hill coefficient H, TTP readout (unitless)") # Table 2, time-to-positivity readout, H, Day 4 = 2.34 +/- 0.33
    lec50_ttp_d4 <- log(0.45)
    label("Day 4 AUC0-24/MIC at half-maximal effect EC50, TTP readout (h)") # Table 2, time-to-positivity readout, EC50, Day 4 = 0.45 +/- 0.02

    le0_ttp_d7 <- log(0.86)
    label("Day 7 drug-free time-to-positivity Econ (day)") # Table 2, time-to-positivity readout, Econ, Day 7 = 0.86 +/- 0.11
    lemax_ttp_d7 <- log(0.83)
    label("Day 7 maximal increase in time-to-positivity Emax (day)") # Table 2, time-to-positivity readout, Emax, Day 7 = 0.83 +/- 2.49
    lhill_ttp_d7 <- log(1.0)
    label("Day 7 Hill coefficient H, TTP readout (unitless)") # Table 2, time-to-positivity readout, H, Day 7 = 1.0 +/- 0.71
    lec50_ttp_d7 <- log(3.84)
    label("Day 7 AUC0-24/MIC at half-maximal effect EC50, TTP readout (h)") # Table 2, time-to-positivity readout, EC50, Day 7 = 3.84 +/- 21.44

    le0_ttp_d14 <- log(0.70)
    label("Day 14 drug-free time-to-positivity Econ (day)") # Table 2, time-to-positivity readout, Econ, Day 14 = 0.70 +/- 0.18
    lemax_ttp_d14 <- log(1.04)
    label("Day 14 maximal increase in time-to-positivity Emax (day)") # Table 2, time-to-positivity readout, Emax, Day 14 = 1.04 +/- 0.94
    lhill_ttp_d14 <- log(1.0)
    label("Day 14 Hill coefficient H, TTP readout (unitless)") # Table 2, time-to-positivity readout, H, Day 14 = 1.0 +/- 0.25
    lec50_ttp_d14 <- log(2.10)
    label("Day 14 AUC0-24/MIC at half-maximal effect EC50, TTP readout (h)") # Table 2, time-to-positivity readout, EC50, Day 14 = 2.10 +/- 5.16

    le0_ttp_d21 <- log(0.48)
    label("Day 21 drug-free time-to-positivity Econ (day)") # Table 2, time-to-positivity readout, Econ, Day 21 = 0.48 +/- 0.08
    lemax_ttp_d21 <- log(0.45)
    label("Day 21 maximal increase in time-to-positivity Emax (day)") # Table 2, time-to-positivity readout, Emax, Day 21 = 0.45 +/- 0.17
    lhill_ttp_d21 <- log(1.0)
    label("Day 21 Hill coefficient H, TTP readout (unitless)") # Table 2, time-to-positivity readout, H, Day 21 = 1.0 +/- 0.1
    lec50_ttp_d21 <- log(0.98)
    label("Day 21 AUC0-24/MIC at half-maximal effect EC50, TTP readout (h)") # Table 2, time-to-positivity readout, EC50, Day 21 = 0.98 +/- 1.40

    # ----- Ethambutol-resistant subpopulation (Figure 4C) ---------------
    # log_cfu_res = res0 + res1 * X + res2 * X^2, X = AUC0-24/MIC.
    # Not tabulated in Gumbo 2026: coefficients obtained by the maintainers
    # from a quadratic least-squares fit to the pixels of each day's
    # regression curve in Figure 4C (axes calibrated on the printed ticks).
    res0_d4 <- 6.30
    label("Day 4 resistant subpopulation at zero exposure (log10 CFU/mL)") # Digitised by the maintainers from Gumbo 2026 Figure 4C, Day 4 curve
    res1_d4 <- -1.42
    label("Day 4 linear coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 4 curve
    res2_d4 <- 0.67
    label("Day 4 quadratic coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per squared AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 4 curve

    res0_d7 <- 6.82
    label("Day 7 resistant subpopulation at zero exposure (log10 CFU/mL)") # Digitised from Gumbo 2026 Figure 4C, Day 7 curve
    res1_d7 <- -2.29
    label("Day 7 linear coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 7 curve
    res2_d7 <- 1.08
    label("Day 7 quadratic coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per squared AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 7 curve

    res0_d14 <- 8.02
    label("Day 14 resistant subpopulation at zero exposure (log10 CFU/mL)") # Digitised from Gumbo 2026 Figure 4C, Day 14 curve
    res1_d14 <- -2.81
    label("Day 14 linear coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 14 curve
    res2_d14 <- 1.27
    label("Day 14 quadratic coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per squared AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 14 curve

    res0_d21 <- 7.15
    label("Day 21 resistant subpopulation at zero exposure (log10 CFU/mL)") # Digitised from Gumbo 2026 Figure 4C, Day 21 curve
    res1_d21 <- -3.77
    label("Day 21 linear coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 21 curve
    res2_d21 <- 1.79
    label("Day 21 quadratic coefficient of the resistant-subpopulation quadratic (log10 CFU/mL per squared AUC0-24/MIC unit)") # Digitised from Gumbo 2026 Figure 4C, Day 21 curve

    # ----- Residual error ----------------------------------------------
    # Table 2 reports standard errors of the estimates and r2 only; no
    # residual SD, so each endpoint is deterministic.
    addSd_log_cfu <- fixed(0)
    label("Additive residual SD, total log10 CFU/mL (0; not reported)") # Gumbo 2026 Table 2 gives r2 only, no residual SD
    addSd_ttp <- fixed(0)
    label("Additive residual SD, time-to-positivity (0; not reported)") # Gumbo 2026 Table 2 gives r2 only, no residual SD
    addSd_log_cfu_res <- fixed(0)
    label("Additive residual SD, resistant log10 CFU/mL (0; not reported)") # Gumbo 2026 Figure 4C, no residual SD reported
  })

  model({
    # PK/PD exposure index AUC0-24/MIC (h).
    aucmic <- AUC_ETHAMBUTOL / mic

    # Each sampling day has its own independent fit. The day is chosen from
    # the simulation time (days) with midpoints between sampling days as
    # cut points; predictions are meaningful only at t = 4, 7, 14 and 21.
    is_d4 <- (t < 5.5)
    is_d7 <- (t >= 5.5) * (t < 10.5)
    is_d14 <- (t >= 10.5) * (t < 17.5)
    is_d21 <- (t >= 17.5)

    e0_cfu <- is_d4 * exp(le0_cfu_d4) + is_d7 * exp(le0_cfu_d7) + is_d14 * exp(le0_cfu_d14) + is_d21 * exp(le0_cfu_d21)
    emax_cfu <- is_d4 * exp(lemax_cfu_d4) + is_d7 * exp(lemax_cfu_d7) + is_d14 * exp(lemax_cfu_d14) + is_d21 * exp(lemax_cfu_d21)
    hill_cfu <- is_d4 * exp(lhill_cfu_d4) + is_d7 * exp(lhill_cfu_d7) + is_d14 * exp(lhill_cfu_d14) + is_d21 * exp(lhill_cfu_d21)
    ec50_cfu <- is_d4 * exp(lec50_cfu_d4) + is_d7 * exp(lec50_cfu_d7) + is_d14 * exp(lec50_cfu_d14) + is_d21 * exp(lec50_cfu_d21)

    e0_ttp <- is_d4 * exp(le0_ttp_d4) + is_d7 * exp(le0_ttp_d7) + is_d14 * exp(le0_ttp_d14) + is_d21 * exp(le0_ttp_d21)
    emax_ttp <- is_d4 * exp(lemax_ttp_d4) + is_d7 * exp(lemax_ttp_d7) + is_d14 * exp(lemax_ttp_d14) + is_d21 * exp(lemax_ttp_d21)
    hill_ttp <- is_d4 * exp(lhill_ttp_d4) + is_d7 * exp(lhill_ttp_d7) + is_d14 * exp(lhill_ttp_d14) + is_d21 * exp(lhill_ttp_d21)
    ec50_ttp <- is_d4 * exp(lec50_ttp_d4) + is_d7 * exp(lec50_ttp_d7) + is_d14 * exp(lec50_ttp_d14) + is_d21 * exp(lec50_ttp_d21)

    res0 <- is_d4 * res0_d4 + is_d7 * res0_d7 + is_d14 * res0_d14 + is_d21 * res0_d21
    res1 <- is_d4 * res1_d4 + is_d7 * res1_d7 + is_d14 * res1_d14 + is_d21 * res1_d21
    res2 <- is_d4 * res2_d4 + is_d7 * res2_d7 + is_d14 * res2_d14 + is_d21 * res2_d21

    # Inhibitory sigmoid Emax on total burden (Gumbo 2026 Table 2).
    log_cfu <- e0_cfu - emax_cfu * aucmic^hill_cfu / (aucmic^hill_cfu + ec50_cfu^hill_cfu)
    # Time-to-positivity rises as the burden falls (Gumbo 2026 Figure 4A).
    ttp <- e0_ttp + emax_ttp * aucmic^hill_ttp / (aucmic^hill_ttp + ec50_ttp^hill_ttp)
    # Quadratic resistant subpopulation (Gumbo 2026 Figure 4C); defined only
    # over the tested range AUC0-24/MIC 0-2.07.
    log_cfu_res <- res0 + res1 * aucmic + res2 * aucmic^2

    log_cfu ~ add(addSd_log_cfu)
    ttp ~ add(addSd_ttp)
    log_cfu_res ~ add(addSd_log_cfu_res)
  })
}
