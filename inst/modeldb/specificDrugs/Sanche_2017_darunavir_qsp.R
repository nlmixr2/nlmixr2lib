Sanche_2017_darunavir_qsp <- function() {
  description <- "QSP. Deterministic (mean-field) HIV-1 viral-dynamics model for ritonavir-boosted darunavir monotherapy (Sanche 2017). A two-compartment first-order-absorption population PK model (Molto 2013) drives the unaffected fraction U = 1 / (1 + (C / (IC50 * kp * kl))^m) in lymph-node CD4+ T cells, with darunavir's 1% lymph-node penetration (kl = 100); ritonavir acts only through darunavir's PK. One (wild-type) strain follows the supplement's boxed delay differential equation with a 1-day eclipse period; no resistance mutation is modelled, as in the paper. Expected-value form of Sanche_2017_darunavir_qsp_tauleap. VL is the viral load (copies/mL). Do not dose with ss = 1: rxode2 would also move the infected-cell state to its steady-state solution, which is not a valid start for a delay equation. Simulation only: no residual error."
  reference <- paste(
    "Sanche S, Sheehan N, Mesplede T, Wainberg MA, Li J, Nekka F. (2017).",
    "A mathematical model to predict HIV virological failure and elucidate",
    "the role of lymph node drug penetration. CPT Pharmacometrics Syst",
    "Pharmacol 6(7):469-476. doi:10.1002/psp4.12200.",
    "Viral-dynamics framework, per-mutation resistance parameters and the",
    "setpoint viral-load distribution from Rosenbloom DIS, Hill AL, Rabi SA,",
    "Siliciano RF, Nowak MA. (2012). Antiretroviral dynamics determines HIV",
    "evolution and predicts therapy outcome. Nat Med 18(9):1378-1385.",
    "doi:10.1038/nm.2892.",
    sep = " "
  )
  vignette <- "Sanche_2017_hiv_virological_failure"

  # Productively infected CD4+ T cells (wild type) have no canonical compartment analogue.
  paper_specific_compartments <- c("inf_wt")

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L for drug concentrations; productively infected CD4+ T cells (count) for the inf_* states; copies/mL for VL"
  )

  compartmentData <- list(
    depot = list(analyte = "darunavir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "darunavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "darunavir", units = "mg", specimen = "tissue", verified = TRUE),
    inf_wt = list(
      analyte = "productively infected CD4+ T cells, wild-type strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    HIV_VLOAD = list(
      description = "Pre-treatment (setpoint) plasma HIV-1 RNA viral load",
      units = "copies/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "The untreated viral setpoint, which fixes the CD4+ T-cell supply rate lambda = setpoint_cells * dy * R0 / (R0 - 1) with setpoint_cells = HIV_VLOAD * 1500 (Rosenbloom 2012 Suppl. section 1.2; 3 L plasma, 2 RNA copies per virion, as in the authors' MATLAB code). For a treatment-naive start (start_suppressed = 0) the trial also starts at this setpoint; for a suppressed start it only sets lambda. Sanche 2017 samples it from Rosenbloom 2012 Suppl. Fig. 2a (log10 bins of width 0.5 from 1.5 to 6.5; the authors' example input file uses the bin midpoints).",
      source_name = "viral setpoint"
    )
  )

  population <- list(
    species = "human",
    disease_state = "Virtual HIV-1-infected adults, virologically suppressed on established therapy, switched to once-daily darunavir/ritonavir 800/100 mg monotherapy for 48 weeks (compared with the MONOI/MONET-type darunavir/ritonavir monotherapy arm of Antinori 2015, 14% virological failure).",
    dose_range = "darunavir/ritonavir 800/100 mg orally once daily",
    regions = "Not applicable (virtual patients; parameters from published population PK models, in vitro assays and viral-kinetic studies)",
    notes = "No subject-level data were fitted: every value is taken from an independent study. The paper simulates groups of 300 virtual patients sharing an adherence level (10-100%), a probability of taking a dose after a missed dose (PIM, 5-100% in 5% steps, weighted per adherence level by Table S1) and an infection phase (acute or chronic; 90% chronic assumed for the clinical comparison), with 1,200-12,000 virtual patients per adherence level. Setpoint viral loads follow the Rosenbloom 2012 Suppl. Fig. 2a distribution."
  )

  ini({
    # ---- Darunavir PK: two-compartment, first-order absorption (Molto 2013) ----
    lka <- fixed(log(0.3)); label("Darunavir absorption rate constant ka (1/h)")  # Supplement 'Summary of parameter values': ka = 0.3 h^-1
    lcl <- fixed(log(10.7)); label("Darunavir apparent clearance CL/F (L/h)")  # Supplement: CL = 10.7 L/h
    lvc <- fixed(log(36.9)); label("Darunavir apparent central volume V1/F (L)")  # Supplement: V1 = 36.9 L
    lvp <- fixed(log(1130)); label("Darunavir apparent peripheral volume V2/F (L)")  # Supplement: V2 = 1130 L
    lq <- fixed(log(21.3)); label("Darunavir apparent intercompartmental clearance Q/F (L/h)")  # Supplement: Q = 21.3 L/h
    
    # ---- Drug effect in lymph-node CD4+ T cells ----
    lic50_drv <- fixed(-6.28); label("Darunavir log IC50 without plasma protein, wild type (log mg/L)")  # Table S2: mu = -6.28
    kp_drv <- fixed(14); label("Darunavir IC50 fold-change for 100% plasma serum, kp (unitless)")  # Table S2: kp = 14
    kl_drv <- fixed(100); label("Darunavir PBMC-to-lymph-node-cell concentration ratio, kl (unitless)")  # Table S2: kl = 100.0 (1% penetration)
    hill_drv <- fixed(3.55); label("Darunavir Hill slope m (unitless)")  # Rosenbloom 2012 Suppl. Table 1: DRV/r slope 3.55
    
    # ---- Within-host viral dynamics ----
    lr_chronic <- fixed(-0.78); label("Mean log viral expansion rate, chronic phase (log 1/day)")  # Supplement: chronic ln(r) ~ N(-0.78, 0.52^2) (Ioannidis 2000)
    lrsd_chronic <- fixed(0.52); label("SD of log viral expansion rate, chronic phase (log 1/day)")  # Supplement: sigma = 0.52
    lr_acute <- fixed(-0.07); label("Mean log viral expansion rate, acute phase (log 1/day)")  # Supplement: acute ln(r) ~ N(-0.07, 0.44^2) (Ribeiro 2010)
    lrsd_acute <- fixed(0.44); label("SD of log viral expansion rate, acute phase (log 1/day)")  # Supplement: sigma = 0.44
    acute_phase <- fixed(0); label("Infection-phase switch: 0 = chronic, 1 = acute (set per simulated group)")  # Methods 'Simulations': groups share an infection phase
    tau_ecl <- fixed(1); label("Eclipse period tau (day)")  # Supplement 'Eclipse Period': tau = 1 day
    dy <- fixed(1); label("Death rate of productively infected cells d_y (1/day)")  # Main text Methods: 1 per day (Markowitz 2003)
    atot <- fixed(3000); label("Total latent-reservoir activation rate sum A_i (cells/day)")  # Main text Eq. (4) and Supplement: 3,000 cells per day
    cells_per_copy <- fixed(1500); label("Infected cells per plasma HIV-1 RNA copy/mL (cells per copies/mL)")  # Authors' MATLAB code: viral load = sum(y) / 1500 (3 L plasma, 2 RNA copies per virion)
    start_suppressed <- fixed(0); label("Trial-start switch: 0 = untreated at setpoint, 1 = suppressed on established therapy")  # Rosenbloom 2012 Suppl. section 1.3: suppression vs maintenance trials
    
    # ---- Between-subject variability (all values taken from the cited sources, none estimated here) ----
    etalcl ~ fixed(0.101124)  # Supplement: darunavir IIV_CL = 0.318, variance 0.318^2
    etalvc ~ fixed(0.346921)  # Supplement: darunavir IIV_V1 = 0.589, variance 0.589^2
    etalvp ~ fixed(0.674041)  # Supplement: darunavir IIV_V2 = 0.821, variance 0.821^2
    etalic50_drv ~ fixed(0.09)  # Table S2: darunavir IC50 log-normal sigma = 0.3, variance 0.3^2
    etalr ~ fixed(1)  # standard-normal deviate; scaled in model() by the phase-specific SD lrsd_chronic / lrsd_acute
  })
  model({
    # ---- Darunavir PK (ritonavir-boosted; time in days, source rates per hour) ----
    ka <- exp(lka) * 24
    cl <- exp(lcl + etalcl) * 24
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq) * 24

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    Cc <- central / vc

    # ---- Per-subject viral and drug-susceptibility parameters ----
    lic50i_drv <- lic50_drv + etalic50_drv
    # viral expansion rate r (1/day): log-normal, phase-specific mean and SD; R0 = (1 + r / dy) * exp(r * tau)
    rgrowth <- exp((1 - acute_phase) * lr_chronic + acute_phase * lr_acute + etalr * ((1 - acute_phase) * lrsd_chronic + acute_phase * lrsd_acute))
    r0 <- (1 + rgrowth / dy) * exp(rgrowth * tau_ecl)
    fit_wt <- 1
    f_wt <- 1
    # setpoint (cells) and CD4+ T-cell supply rate lambda (cells/day) from the pre-treatment setpoint viral load
    setpoint <- HIV_VLOAD * cells_per_copy
    lambda <- setpoint * dy * r0 / (r0 - 1)
    y0tot <- (1 - start_suppressed) * setpoint + start_suppressed * atot / dy

    # ---- Delay terms: state and concentration one eclipse period (tau) earlier ----
    # delay() returns the initial condition before t = 0. Before t = tau the lagged
    # concentration is 0 for a treatment-naive start and equals the current one for a
    # start on established once-daily therapy, whose profile repeats every 24 h.
    clag_drv <- delay(central, tau_ecl) / vc
    c_drv_lag <- (t >= tau_ecl) * clag_drv + (t < tau_ecl) * start_suppressed * Cc
    ylag_wt <- delay(inf_wt, tau_ecl)

    # ---- Drug effect at the infection time t - tau ----
    if (c_drv_lag > 0) { ufr_wt <- 1 / (1 + (c_drv_lag / (exp(lic50i_drv) * kp_drv * kl_drv))^hill_drv) } else { ufr_wt <- 1 }
    rx_wt <- r0 * fit_wt * ufr_wt

    # ---- Productively infected cells (supplement boxed delay equation + de novo mutation transfer) ----
    yden <- lambda + dy * (rx_wt * ylag_wt)
    ninf_wt <- dy * lambda * rx_wt * ylag_wt / yden
    d/dt(inf_wt) <- atot * f_wt - dy * inf_wt + ninf_wt
    inf_wt(0) <- y0tot * f_wt

    # ---- Outputs ----
    ytot <- inf_wt
    VL <- ytot / cells_per_copy
  })
}
