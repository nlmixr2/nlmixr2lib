Wang_2024_amphenmulin_killrate <- function() {
  description <- "In vitro (Mycoplasma gallisepticum strain S6). Concentration-dependent kill-rate PK/PD model for amphenmulin, a novel pleuromutilin derivative, coupling the one-compartment first-order-elimination PK of Wang 2024's in-vitro dynamic system to the sigmoid Emax kill-rate relationship fitted to the static time-kill curves. Wang 2024 Materials and Methods (Static time-kill curves) parameterises the kill rate as E = E0 + Emax * Ce^N / (EC50^N + Ce^N), where E is the kill rate in log10(CFU/mL) per hour (POSITIVE = killing), E0 is the corresponding rate of change in the untreated control (NEGATIVE, i.e. net growth), Ce is the amphenmulin concentration and Emax is an INCREMENT above E0 rather than the asymptote, so the saturating kill rate is E0 + Emax. Note that Wang 2024's Results text calls Emax 'the maximum kill rate'; the printed equation, which matches the Phoenix Sigmoid Emax model the authors used, is taken as authoritative. Parameters are Wang 2024 Table 2, row 0-24 h, which the paper identifies as the optimal fit (R = 0.9936) and plots in Figure 4b: Emax = 0.1261 1/h, EC50 = 0.0325 ug/mL, E0 = -0.0093 1/h, Hill N = 0.9238. The bacterial density bact (linear CFU/mL) is integrated as d/dt(bact) = -ln(10) * E * bact so that log10(bact) falls at exactly E log10(CFU/mL) per hour, reproducing the paper's fitted kill rate exactly at any constant concentration. The PK is the paper's in-vitro dynamic model (Materials and Methods, In vitro dynamic model): first-order elimination C = C0 * exp(-k*t) with k set from the chicken intravenous half-life T1/2Ke = 2.13 h, giving k = 0.3254 1/h, realised in the apparatus as a 1.63 mL/min flow through the 300 mL reaction chamber (1.63*60/300 = 0.326 1/h). IMPORTANT PROVENANCE NOTE: Wang 2024 fitted the kill-rate parameters to STATIC concentrations and ran the dynamic apparatus as a separate experiment; the two were not fitted jointly. This model composes them, which is the paper's stated PK/PD integration written as an ODE system. Setting the elimination rate to zero recovers the static experiment and reproduces Table 2 exactly. Wang 2024 reports no between-subject variability and no residual error magnitude for either component, so no eta parameters are present and both residual SDs are held at zero for typical-value simulation."
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
    dosing = "ug (amphenmulin into the 300 mL reaction chamber; amt/300 is the resulting concentration in ug/mL)",
    concentration = "ug/mL (Cc); log10 CFU/mL (logCfu)"
  )

  paper_specific_compartments <- c("bact")

  compartmentData <- list(
    central = list(analyte = "amphenmulin", units = "ug", specimen = "not applicable", verified = TRUE),
    bact    = list(analyte = "Mycoplasma gallisepticum strain S6", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "in vitro (Mycoplasma gallisepticum strain S6 culture)",
    n_subjects     = 3L,
    n_studies      = 1L,
    organism       = "Mycoplasma gallisepticum standard strain S6 (National Center for Veterinary Culture Collection, Beijing), cultured at 37 C in M. gallisepticum artificial medium base supplemented with swine serum, NADH and cysteine. Amphenmulin MIC = 0.0039 ug/mL by broth microdilution and 0.0078 ug/mL by agar dilution; MIC99 = 0.0077 ug/mL, MBC = 0.0156 ug/mL, MPC = 0.0500 ug/mL at a 10^9 CFU/mL inoculum, giving a mutant selection window of 0.0077-0.05 ug/mL",
    system         = "Static time-kill curves in 10 mL Celine bottles (7 mL blank medium + 0.2 mL amphenmulin solution + 0.8 mL mycoplasma suspension) for the kill-rate parameters; the PK component is the separate in-vitro dynamic apparatus, a reservoir chamber of fresh drug-free medium feeding a three-necked reaction flask holding 300 mL of medium in the external compartment and a 10 mL semipermeable-membrane internal compartment carrying the culture, draining to a waste chamber under peristaltic-pump control at 1.63 mL/min",
    disease_state  = "Not applicable (in-vitro culture); starting inoculum 10^7 CFU/mL",
    dose_range     = "Static time-kill: 1, 2, 4, 6, 8, 16, 32 and 64 times MIC (0.0039 to 0.2496 ug/mL) plus a drug-free growth control. In-vitro dynamic apparatus: peak concentrations of 0.1, 0.2, 0.5, 0.8, 1.0 and 1.5 ug/mL every 24 h for three doses",
    sampling       = "Static: 100 uL sampled at 0, 6, 9, 12, 24, 36 and 48 h for colony counting, lower counting limit 100 CFU/mL, each experiment repeated three times",
    regions        = "China (South China Agricultural University, Guangzhou)",
    notes          = "Wang 2024 Table 2 reports the sigmoid Emax kill-rate fit over seven time windows (0-24, 0-36, 0-48, 6-24, 6-36, 6-48 and 12-48 h) with R from 0.9649 to 0.9936. This model packages the 0-24 h window, which the paper identifies as the optimal fit and plots in Figure 4b; the other six parameterisations are tabulated in the vignette. The kill-rate model has no carrying-capacity term, so the drug-free control grows without bound; over the 48-72 h horizon of the source experiments this is immaterial (the control gains about 0.22 log10 CFU/mL per 24 h)."
  )

  ini({
    # =================================================================
    # Wang 2024 Table 2, row 0-24 h -- sigmoid Emax kill rate
    # =================================================================
    # Wang 2024 Materials and Methods (Static time-kill curves):
    #   E = E0 + Emax * Ce^N / (EC50^N + Ce^N)
    # E is the kill rate in log10(CFU/mL) per hour, positive = killing.
    # Emax is an INCREMENT above E0 in this printed form (the Phoenix
    # Sigmoid Emax parameterisation the authors used), so the
    # saturating kill rate is E0 + Emax = 0.1168 1/h, not 0.1261 1/h.
    # The Results text loosely calls Emax "the maximum kill rate"; the
    # printed equation governs. See the vignette Errata.
    #
    # E0 stays on the natural scale because it is NEGATIVE: the
    # untreated control grows, so its kill rate is below zero and the
    # parameter is not sign-constrained.
    e0 <- -0.0093
    label("Kill rate of the untreated control E0 (log10 CFU/mL per h; negative = net growth)")  # Wang 2024 Table 2, row 0-24 h, E0 = -0.0093 1/h
    lemax <- log(0.1261)
    label("Maximum increment in kill rate above the control Emax (log10 CFU/mL per h)")  # Wang 2024 Table 2, row 0-24 h, Emax = 0.1261 1/h
    lec50 <- log(0.0325)
    label("Amphenmulin concentration producing half the maximal kill-rate increment EC50 (ug/mL)")  # Wang 2024 Table 2, row 0-24 h, EC50 = 0.0325 ug/mL
    lhill <- log(0.9238)
    label("Hill coefficient N defining the slope of the kill-rate curve (unitless)")  # Wang 2024 Table 2, row 0-24 h, Hill's slope = 0.9238

    # =================================================================
    # In-vitro dynamic-system PK
    # =================================================================
    # Wang 2024 Materials and Methods (In vitro dynamic model):
    # C = C0 * exp(-k*t), with the apparatus flow rate set from the
    # chicken intravenous elimination half-life T1/2Ke = 2.13 h
    # (Table 1), so k = log(2)/2.13 = 0.3254 1/h. The apparatus
    # realises this as 1.63 mL/min through the 300 mL reaction
    # chamber: 1.63*60/300 = 0.326 1/h. Both are design constants of
    # the experiment rather than estimated quantities.
    lkel <- fixed(log(0.3254))
    label("First-order elimination rate constant of the in-vitro dynamic system (1/h)")  # Wang 2024 In vitro dynamic model; k = log(2)/T1/2Ke with T1/2Ke = 2.13 h from Table 1
    lvc <- fixed(log(300))
    label("Reaction-chamber volume (mL)")  # Wang 2024 In vitro dynamic model: 300 mL of medium in the external compartment of the reaction chamber

    # =================================================================
    # Starting bacterial density
    # =================================================================
    # Experimental design input, not an estimated parameter.
    log10_cfu0 <- fixed(7)
    label("log10 starting mycoplasma density (log10 CFU/mL)")  # Wang 2024 Static time-kill curves: final mycoplasma count of 10^7 CFU/mL

    # =================================================================
    # Residual error
    # =================================================================
    # Wang 2024 reported only the correlation coefficient of the
    # sigmoid fit (R = 0.9936 for the 0-24 h window) and no residual
    # standard deviation for either the concentrations or the counts,
    # so both residual SDs are held at zero for deterministic
    # typical-value simulation.
    addSd <- fixed(0)
    label("Additive residual SD on amphenmulin concentration (ug/mL; not reported in Wang 2024)")  # Wang 2024 reports R only, no residual SD
    addSd_logCfu <- fixed(0)
    label("Additive residual SD on log10 CFU/mL (not reported in Wang 2024)")  # Wang 2024 reports R only, no residual SD
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    kel <- exp(lkel)
    vc <- exp(lvc)

    # In-vitro dynamic-system PK: first-order elimination only. A bolus
    # of `amt` ug into `central` gives a peak concentration of
    # amt/300 ug/mL, so the paper's six regimens (0.1, 0.2, 0.5, 0.8,
    # 1.0 and 1.5 ug/mL) correspond to amt = 30, 60, 150, 240, 300 and
    # 450 ug.
    d/dt(central) <- -kel * central
    Cc <- central / vc

    # Wang 2024 sigmoid Emax kill rate. kill_log10 is the rate of
    # decline of log10(CFU/mL) in log10 units per hour; it is negative
    # at zero concentration (the control grows) and rises towards
    # e0 + emax at saturating concentration.
    kill_log10 <- e0 + emax * Cc^hill / (ec50^hill + Cc^hill)

    # d(log10 N)/dt = -kill_log10  =>  d(N)/dt = -ln(10)*kill_log10*N
    d/dt(bact) <- -log(10) * kill_log10 * bact
    bact(0) <- 10^log10_cfu0

    # log10 CFU/mL observation with a 1-CFU/mL floor so the log10 stays
    # finite if bact is driven below 1 CFU/mL (the same in-vitro PD
    # convention used by Chen_2023_tilmicosin).
    logCfu <- log10(bact + 1)

    Cc ~ add(addSd)
    logCfu ~ add(addSd_logCfu)
  })
}
