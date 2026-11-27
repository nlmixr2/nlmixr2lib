Sanche_2017_efavirenz_tenofovir_emtricitabine_qsp <- function() {
  description <- "QSP. Deterministic (mean-field) HIV-1 viral-dynamics and resistance model for efavirenz monotherapy or efavirenz + tenofovir disoproxil fumarate + emtricitabine (Sanche 2017). Population PK models (efavirenz one-compartment; tenofovir and emtricitabine two-compartment, first-order absorption; dose the regimen's drugs into depot, depot_tfv and depot_ftc) drive the fraction of single-cycle infections unaffected in lymph-node CD4+ T cells, U = 1 / (1 + (C / (rho * IC50 * kp * kl))^(m (1 + sigma))), combined across drugs by the Jilek 2012 degree of independence between Bliss independence and Loewe additivity. Sixteen viral strains (wild type, K103N, Y181C, G190S, M184V, K65R and their ten double mutants) follow the supplement's boxed delay differential equation dy_i/dt = A_i - dy y_i + dy lambda R_i(t - tau) y_i(t - tau) / (lambda + dy sum_j R_j(t - tau) y_j(t - tau)) with a 1-day eclipse period, plus the de novo mutation transfer of the main-text Eq. (3). This is the expected-value form of the authors' as-run stochastic tau-leap (Sanche_2017_efavirenz_tenofovir_emtricitabine_qsp_tauleap); because mutant cells are carried as fractional expected counts it predicts resistance emergence whenever any resistant strain can grow, so resistance probabilities need the stochastic file. VL is the viral load (copies/mL) and fmut the mutant fraction. Do not dose with ss = 1: rxode2 would also move the infected-cell states to its steady-state solution, which is not a valid start for a delay equation. Simulation only: no residual error."
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
  dosing <- c("depot", "depot_tfv", "depot_ftc")

  # Productively infected CD4+ T cells per viral strain (wild type, five single
  # mutants, ten double mutants) have no canonical compartment analogue.
  paper_specific_compartments <- c(
    "inf_wt",
    "inf_k103n",
    "inf_y181c",
    "inf_g190s",
    "inf_m184v",
    "inf_k65r",
    "inf_k103n_y181c",
    "inf_k103n_g190s",
    "inf_k103n_m184v",
    "inf_k103n_k65r",
    "inf_y181c_g190s",
    "inf_y181c_m184v",
    "inf_y181c_k65r",
    "inf_g190s_m184v",
    "inf_g190s_k65r",
    "inf_m184v_k65r"
  )

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "mg/L for drug concentrations; productively infected CD4+ T cells (count) for the inf_* states; copies/mL for VL"
  )

  compartmentData <- list(
    depot = list(analyte = "efavirenz", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "efavirenz", units = "mg", specimen = "plasma", verified = TRUE),
    depot_tfv = list(
      analyte = "tenofovir",
      units = "mg (tenofovir equivalents)",
      specimen = "administration site",
      verified = TRUE
    ),
    central_tfv = list(analyte = "tenofovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_tfv = list(analyte = "tenofovir", units = "mg", specimen = "tissue", verified = TRUE),
    depot_ftc = list(analyte = "emtricitabine", units = "mg", specimen = "administration site", verified = TRUE),
    central_ftc = list(analyte = "emtricitabine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ftc = list(analyte = "emtricitabine", units = "mg", specimen = "tissue", verified = TRUE),
    inf_wt = list(
      analyte = "productively infected CD4+ T cells, wild-type strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_k103n = list(
      analyte = "productively infected CD4+ T cells, K103N strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_y181c = list(
      analyte = "productively infected CD4+ T cells, Y181C strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_g190s = list(
      analyte = "productively infected CD4+ T cells, G190S strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_m184v = list(
      analyte = "productively infected CD4+ T cells, M184V strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_k65r = list(
      analyte = "productively infected CD4+ T cells, K65R strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_k103n_y181c = list(
      analyte = "productively infected CD4+ T cells, K103N+Y181C strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_k103n_g190s = list(
      analyte = "productively infected CD4+ T cells, K103N+G190S strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_k103n_m184v = list(
      analyte = "productively infected CD4+ T cells, K103N+M184V strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_k103n_k65r = list(
      analyte = "productively infected CD4+ T cells, K103N+K65R strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_y181c_g190s = list(
      analyte = "productively infected CD4+ T cells, Y181C+G190S strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_y181c_m184v = list(
      analyte = "productively infected CD4+ T cells, Y181C+M184V strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_y181c_k65r = list(
      analyte = "productively infected CD4+ T cells, Y181C+K65R strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_g190s_m184v = list(
      analyte = "productively infected CD4+ T cells, G190S+M184V strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_g190s_k65r = list(
      analyte = "productively infected CD4+ T cells, G190S+K65R strain",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    inf_m184v_k65r = list(
      analyte = "productively infected CD4+ T cells, M184V+K65R strain",
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
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Raw (not BSA-normalised) creatinine clearance in mL/min, used only for emtricitabine clearance as the power term (CRCL / 104.5)^0.278. The supplement states that creatinine clearance 'in liters per hour followed a log-N(mu, sigma^2) distribution with mu = 104.5 and sigma = 0.34'; 104.5 L/h is not physiological, so the maintainers read it as a median of 104.5 mL/min with a log-scale SD of 0.34, and centre the power model on that median.",
      source_name = "creatinine clearance"
    )
  )

  population <- list(
    species = "human",
    disease_state = "Virtual HIV-1-infected adults: (a) treatment-naive patients at their viral setpoint starting once-daily efavirenz 600 mg monotherapy for 48 weeks, and (b) patients on once-daily efavirenz 600 mg + tenofovir disoproxil fumarate 300 mg + emtricitabine 200 mg followed for 12-18 months (compared with the pharmacy-refill adherence cohort of Gordon 2015).",
    dose_range = "efavirenz 600 mg, tenofovir disoproxil fumarate 300 mg and emtricitabine 200 mg, orally once daily",
    regions = "Not applicable (virtual patients; parameters from published population PK models, in vitro assays and viral-kinetic studies)",
    notes = "No subject-level data were fitted: every value is taken from an independent study. The paper simulates groups of 300 virtual patients sharing an adherence level (10-100%), a probability of taking a dose after a missed dose (PIM, 5-100% in 5% steps, weighted per adherence level by Table S1) and an infection phase (acute or chronic; 90% chronic assumed for the clinical comparison), with 1,200-12,000 virtual patients per adherence level. Setpoint viral loads follow the Rosenbloom 2012 Suppl. Fig. 2a distribution."
  )

  ini({
    # ---- Efavirenz PK: one-compartment, first-order absorption (Sanchez 2011 base model) ----
    lka <- fixed(log(0.3)); label("Efavirenz absorption rate constant ka (1/h)")  # Supplement 'Summary of parameter values': ka = 0.3 h^-1
    lcl <- fixed(log(9.61)); label("Efavirenz apparent clearance CL/F (L/h)")  # Supplement: CL = 9.61 L/h
    lvc <- fixed(log(291)); label("Efavirenz apparent volume V/F (L)")  # Supplement: V = 291 L
    
    # ---- Tenofovir PK: two-compartment, first-order absorption (Baheti 2011) ----
    lka_tfv <- fixed(log(0.822)); label("Tenofovir absorption rate constant ka (1/h)")  # Supplement: ka = 0.822 h^-1
    lcl_tfv <- fixed(log(42.1)); label("Tenofovir apparent clearance CL/F (L/h)")  # Supplement: CL = 42.1 L/h
    lvc_tfv <- fixed(log(224)); label("Tenofovir apparent central volume V1/F (L)")  # Supplement: V1 = 224 L
    lvp_tfv <- fixed(log(512)); label("Tenofovir apparent peripheral volume V2/F (L)")  # Supplement: V2 = 512 L
    lq_tfv <- fixed(log(175)); label("Tenofovir apparent intercompartmental clearance Q/F (L/h)")  # Supplement: Q = 175 L/h
    
    # ---- Emtricitabine PK: two-compartment, first-order absorption, CrCL on CL ----
    lka_ftc <- fixed(log(0.53)); label("Emtricitabine absorption rate constant ka (1/h)")  # Supplement: ka = 0.53 h^-1
    lcl_ftc <- fixed(log(15.1)); label("Emtricitabine apparent clearance CL/F at CrCL 104.5 mL/min (L/h)")  # Supplement: CL = 15.1 L/h
    lvc_ftc <- fixed(log(42.3)); label("Emtricitabine apparent central volume V1/F (L)")  # Supplement: V1 = 42.3 L
    lvp_ftc <- fixed(log(55.4)); label("Emtricitabine apparent peripheral volume V2/F (L)")  # Supplement: V2 = 55.4 L
    lq_ftc <- fixed(log(5.75)); label("Emtricitabine apparent intercompartmental clearance Q/F (L/h)")  # Supplement: Q = 5.75 L/h
    e_crcl_cl_ftc <- fixed(0.278); label("Power exponent of CrCL on emtricitabine CL (unitless)")  # Supplement: theta_CLCR = 0.278
    
    # ---- Drug effect in lymph-node CD4+ T cells (Supplement Table S2, Hill slopes from Rosenbloom 2012 Suppl. Table 1) ----
    lic50_efv <- fixed(-7.49); label("Efavirenz log IC50 without plasma protein, wild type (log mg/L)")  # Table S2: mu = -7.49
    lic50_tfv <- fixed(-1.79); label("Tenofovir log IC50 without plasma protein, wild type (log mg/L)")  # Table S2: mu = -1.79
    lic50_ftc <- fixed(-4.98); label("Emtricitabine log IC50 without plasma protein, wild type (log mg/L)")  # Table S2: mu = -4.98
    kp_efv <- fixed(27); label("Efavirenz IC50 fold-change for 100% plasma serum, kp (unitless)")  # Table S2: kp = 27
    kp_tfv <- fixed(1); label("Tenofovir IC50 fold-change for 100% plasma serum, kp (unitless)")  # Table S2: kp = 1
    kp_ftc <- fixed(1); label("Emtricitabine IC50 fold-change for 100% plasma serum, kp (unitless)")  # Table S2: kp = 1
    kl_efv <- fixed(16.7); label("Efavirenz PBMC-to-lymph-node-cell concentration ratio, kl (unitless)")  # Table S2: kl = 16.7
    kl_tfv <- fixed(5.0); label("Tenofovir PBMC-to-lymph-node-cell concentration ratio, kl (unitless)")  # Table S2: kl = 5.0
    kl_ftc <- fixed(2.9); label("Emtricitabine PBMC-to-lymph-node-cell concentration ratio, kl (unitless)")  # Table S2: kl = 2.9
    hill_efv <- fixed(1.69); label("Efavirenz Hill slope m (unitless)")  # Rosenbloom 2012 Suppl. Table 1: EFV slope 1.69
    hill_tfv <- fixed(0.97); label("Tenofovir Hill slope m (unitless)")  # Rosenbloom 2012 Suppl. Table 1: TDF slope 0.97
    hill_ftc <- fixed(1.20); label("Emtricitabine Hill slope m (unitless)")  # Rosenbloom 2012 Suppl. Table 1: FTC slope 1.20
    di_efv_tfv <- fixed(0.46); label("Degree of independence, efavirenz + tenofovir (unitless)")  # Supplement: DI = 0.46 (Jilek 2012)
    di_efv_ftc <- fixed(0.71); label("Degree of independence, efavirenz + emtricitabine (unitless)")  # Supplement: DI = 0.71 (Jilek 2012)
    di_tfv_ftc <- fixed(0.33); label("Degree of independence, tenofovir + emtricitabine (unitless)")  # Supplement: DI = 0.33 (Jilek 2012)
    
    # ---- Resistance mutations (Rosenbloom 2012 Suppl. Table 3: cost s, mutation probability u, IC50 fold-change rho, slope change sigma) ----
    s_k103n <- fixed(0.3); label("Fitness cost of K103N (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K103N s = 0.3
    u_k103n <- fixed(1.5e-06); label("Probability of de novo K103N per infection (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K103N u = 1.5e-06
    rho_k103n_efv <- fixed(85); label("EFV IC50 fold-change for K103N (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K103N / EFV rho = 85
    sigma_k103n_efv <- fixed(-0.17); label("EFV Hill-slope change for K103N (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K103N / EFV sigma = -0.17
    s_y181c <- fixed(0.26); label("Fitness cost of Y181C (unitless)")  # Rosenbloom 2012 Suppl. Table 3: Y181C s = 0.26
    u_y181c <- fixed(1.1e-05); label("Probability of de novo Y181C per infection (unitless)")  # Rosenbloom 2012 Suppl. Table 3: Y181C u = 1.1e-05
    rho_y181c_efv <- fixed(2.6); label("EFV IC50 fold-change for Y181C (unitless)")  # Rosenbloom 2012 Suppl. Table 3: Y181C / EFV rho = 2.6
    sigma_y181c_efv <- fixed(-0.11); label("EFV Hill-slope change for Y181C (unitless)")  # Rosenbloom 2012 Suppl. Table 3: Y181C / EFV sigma = -0.11
    s_g190s <- fixed(0.79); label("Fitness cost of G190S (unitless)")  # Rosenbloom 2012 Suppl. Table 3: G190S s = 0.79
    u_g190s <- fixed(2.2e-05); label("Probability of de novo G190S per infection (unitless)")  # Rosenbloom 2012 Suppl. Table 3: G190S u = 2.2e-05
    rho_g190s_efv <- fixed(70); label("EFV IC50 fold-change for G190S (unitless)")  # Rosenbloom 2012 Suppl. Table 3: G190S / EFV rho = 70
    sigma_g190s_efv <- fixed(-0.40); label("EFV Hill-slope change for G190S (unitless)")  # Rosenbloom 2012 Suppl. Table 3: G190S / EFV sigma = -0.40
    s_m184v <- fixed(0.46); label("Fitness cost of M184V (unitless)")  # Rosenbloom 2012 Suppl. Table 3: M184V s = 0.46
    u_m184v <- fixed(1.1e-05); label("Probability of de novo M184V per infection (unitless)")  # Rosenbloom 2012 Suppl. Table 3: M184V u = 1.1e-05
    rho_m184v_ftc <- fixed(1186); label("FTC IC50 fold-change for M184V (unitless)")  # Rosenbloom 2012 Suppl. Table 3: M184V / FTC rho = 1186
    sigma_m184v_ftc <- fixed(-0.49); label("FTC Hill-slope change for M184V (unitless)")  # Rosenbloom 2012 Suppl. Table 3: M184V / FTC sigma = -0.49
    rho_m184v_tfv <- fixed(3); label("TFV IC50 fold-change for M184V (unitless)")  # Rosenbloom 2012 Suppl. Table 3: M184V / TDF rho = 3
    sigma_m184v_tfv <- fixed(-0.27); label("TFV Hill-slope change for M184V (unitless)")  # Rosenbloom 2012 Suppl. Table 3: M184V / TDF sigma = -0.27
    s_k65r <- fixed(0.41); label("Fitness cost of K65R (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K65R s = 0.41
    u_k65r <- fixed(1.1e-05); label("Probability of de novo K65R per infection (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K65R u = 1.1e-05
    rho_k65r_ftc <- fixed(29); label("FTC IC50 fold-change for K65R (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K65R / FTC rho = 29
    sigma_k65r_ftc <- fixed(-0.04); label("FTC Hill-slope change for K65R (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K65R / FTC sigma = -0.04
    rho_k65r_tfv <- fixed(43); label("TFV IC50 fold-change for K65R (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K65R / TDF rho = 43
    sigma_k65r_tfv <- fixed(0.00); label("TFV Hill-slope change for K65R (unitless)")  # Rosenbloom 2012 Suppl. Table 3: K65R / TDF sigma = 0.00
    
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
    etalcl ~ fixed(0.190969)  # Supplement: efavirenz IIV_CL = 0.437 (SD of log CL), variance 0.437^2
    etalvc ~ fixed(1.0201)  # Supplement: efavirenz IIV_V = 1.01, variance 1.01^2
    etalcl_tfv ~ fixed(0.128881)  # Supplement: tenofovir IIV_CL = 0.359, variance 0.359^2
    etalvc_tfv ~ fixed(0.630436)  # Supplement: tenofovir IIV_V1 = 0.794, variance 0.794^2
    etalvp_tfv ~ fixed(0.173889)  # Supplement: tenofovir IIV_V2 = 0.417, variance 0.417^2
    etalcl_ftc ~ fixed(0.030276)  # Supplement: emtricitabine IIV_CL = 0.174, variance 0.174^2
    etalic50_efv ~ fixed(0.1681)  # Table S2: efavirenz IC50 log-normal sigma = 0.41, variance 0.41^2
    etalic50_tfv ~ fixed(0.6724)  # Table S2: tenofovir IC50 log-normal sigma = 0.82, variance 0.82^2
    etalic50_ftc ~ fixed(0.9216)  # Table S2: emtricitabine IC50 log-normal sigma = 0.96, variance 0.96^2
    etalr ~ fixed(1)  # standard-normal deviate; scaled in model() by the phase-specific SD lrsd_chronic / lrsd_acute
  })
  model({
    # ---- Pharmacokinetics (time in days; the source rate constants are per hour) ----
    ka <- exp(lka) * 24
    cl <- exp(lcl + etalcl) * 24
    vc <- exp(lvc + etalvc)
    ka_tfv <- exp(lka_tfv) * 24
    cl_tfv <- exp(lcl_tfv + etalcl_tfv) * 24
    vc_tfv <- exp(lvc_tfv + etalvc_tfv)
    vp_tfv <- exp(lvp_tfv + etalvp_tfv)
    q_tfv <- exp(lq_tfv) * 24
    ka_ftc <- exp(lka_ftc) * 24
    cl_ftc <- exp(lcl_ftc + etalcl_ftc) * (CRCL / 104.5)^e_crcl_cl_ftc * 24
    vc_ftc <- exp(lvc_ftc)
    vp_ftc <- exp(lvp_ftc)
    q_ftc <- exp(lq_ftc) * 24

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl / vc * central
    d/dt(depot_tfv) <- -ka_tfv * depot_tfv
    d/dt(central_tfv) <- ka_tfv * depot_tfv - (cl_tfv + q_tfv) / vc_tfv * central_tfv + q_tfv / vp_tfv * peripheral1_tfv
    d/dt(peripheral1_tfv) <- q_tfv / vc_tfv * central_tfv - q_tfv / vp_tfv * peripheral1_tfv
    d/dt(depot_ftc) <- -ka_ftc * depot_ftc
    d/dt(central_ftc) <- ka_ftc * depot_ftc - (cl_ftc + q_ftc) / vc_ftc * central_ftc + q_ftc / vp_ftc * peripheral1_ftc
    d/dt(peripheral1_ftc) <- q_ftc / vc_ftc * central_ftc - q_ftc / vp_ftc * peripheral1_ftc

    Cc <- central / vc
    Cc_tfv <- central_tfv / vc_tfv
    Cc_ftc <- central_ftc / vc_ftc

    # ---- Per-subject viral and drug-susceptibility parameters ----
    lic50i_efv <- lic50_efv + etalic50_efv
    lic50i_tfv <- lic50_tfv + etalic50_tfv
    lic50i_ftc <- lic50_ftc + etalic50_ftc
    # viral expansion rate r (1/day): log-normal, phase-specific mean and SD; R0 = (1 + r / dy) * exp(r * tau)
    rgrowth <- exp((1 - acute_phase) * lr_chronic + acute_phase * lr_acute + etalr * ((1 - acute_phase) * lrsd_chronic + acute_phase * lrsd_acute))
    r0 <- (1 + rgrowth / dy) * exp(rgrowth * tau_ecl)
    # strain replicative fitness (1 - s); double mutants: midpoint of min and product of the single-mutant fitnesses
    fit_wt <- 1
    fit_k103n <- (1 - s_k103n)
    fit_y181c <- (1 - s_y181c)
    fit_g190s <- (1 - s_g190s)
    fit_m184v <- (1 - s_m184v)
    fit_k65r <- (1 - s_k65r)
    fit_k103n_y181c <- (min((1 - s_k103n), (1 - s_y181c)) + (1 - s_k103n) * (1 - s_y181c)) / 2
    fit_k103n_g190s <- (min((1 - s_k103n), (1 - s_g190s)) + (1 - s_k103n) * (1 - s_g190s)) / 2
    fit_k103n_m184v <- (min((1 - s_k103n), (1 - s_m184v)) + (1 - s_k103n) * (1 - s_m184v)) / 2
    fit_k103n_k65r <- (min((1 - s_k103n), (1 - s_k65r)) + (1 - s_k103n) * (1 - s_k65r)) / 2
    fit_y181c_g190s <- (min((1 - s_y181c), (1 - s_g190s)) + (1 - s_y181c) * (1 - s_g190s)) / 2
    fit_y181c_m184v <- (min((1 - s_y181c), (1 - s_m184v)) + (1 - s_y181c) * (1 - s_m184v)) / 2
    fit_y181c_k65r <- (min((1 - s_y181c), (1 - s_k65r)) + (1 - s_y181c) * (1 - s_k65r)) / 2
    fit_g190s_m184v <- (min((1 - s_g190s), (1 - s_m184v)) + (1 - s_g190s) * (1 - s_m184v)) / 2
    fit_g190s_k65r <- (min((1 - s_g190s), (1 - s_k65r)) + (1 - s_g190s) * (1 - s_k65r)) / 2
    fit_m184v_k65r <- (min((1 - s_m184v), (1 - s_k65r)) + (1 - s_m184v) * (1 - s_k65r)) / 2
    # expected strain composition of the infected-cell population in the absence of drug
    w_wt <- 1
    w_k103n <- u_k103n / s_k103n
    w_y181c <- u_y181c / s_y181c
    w_g190s <- u_g190s / s_g190s
    w_m184v <- u_m184v / s_m184v
    w_k65r <- u_k65r / s_k65r
    w_k103n_y181c <- u_k103n * u_y181c / (1 - fit_k103n_y181c)
    w_k103n_g190s <- u_k103n * u_g190s / (1 - fit_k103n_g190s)
    w_k103n_m184v <- u_k103n * u_m184v / (1 - fit_k103n_m184v)
    w_k103n_k65r <- u_k103n * u_k65r / (1 - fit_k103n_k65r)
    w_y181c_g190s <- u_y181c * u_g190s / (1 - fit_y181c_g190s)
    w_y181c_m184v <- u_y181c * u_m184v / (1 - fit_y181c_m184v)
    w_y181c_k65r <- u_y181c * u_k65r / (1 - fit_y181c_k65r)
    w_g190s_m184v <- u_g190s * u_m184v / (1 - fit_g190s_m184v)
    w_g190s_k65r <- u_g190s * u_k65r / (1 - fit_g190s_k65r)
    w_m184v_k65r <- u_m184v * u_k65r / (1 - fit_m184v_k65r)
    wtot <- w_wt + w_k103n + w_y181c + w_g190s + w_m184v + w_k65r + w_k103n_y181c + w_k103n_g190s + w_k103n_m184v + w_k103n_k65r + w_y181c_g190s + w_y181c_m184v + w_y181c_k65r + w_g190s_m184v + w_g190s_k65r + w_m184v_k65r
    f_wt <- w_wt / wtot
    f_k103n <- w_k103n / wtot
    f_y181c <- w_y181c / wtot
    f_g190s <- w_g190s / wtot
    f_m184v <- w_m184v / wtot
    f_k65r <- w_k65r / wtot
    f_k103n_y181c <- w_k103n_y181c / wtot
    f_k103n_g190s <- w_k103n_g190s / wtot
    f_k103n_m184v <- w_k103n_m184v / wtot
    f_k103n_k65r <- w_k103n_k65r / wtot
    f_y181c_g190s <- w_y181c_g190s / wtot
    f_y181c_m184v <- w_y181c_m184v / wtot
    f_y181c_k65r <- w_y181c_k65r / wtot
    f_g190s_m184v <- w_g190s_m184v / wtot
    f_g190s_k65r <- w_g190s_k65r / wtot
    f_m184v_k65r <- w_m184v_k65r / wtot
    # setpoint (cells) and CD4+ T-cell supply rate lambda (cells/day) from the pre-treatment setpoint viral load
    setpoint <- HIV_VLOAD * cells_per_copy
    lambda <- setpoint * dy * r0 / (r0 - 1)
    y0tot <- (1 - start_suppressed) * setpoint + start_suppressed * atot / dy

    # ---- Delay terms: state and concentration one eclipse period (tau) earlier ----
    # delay() returns the initial condition before t = 0. Before t = tau the lagged
    # concentration is 0 for a treatment-naive start and equals the current one for a
    # start on established once-daily therapy, whose profile repeats every 24 h.
    clag_efv <- delay(central, tau_ecl) / vc
    clag_tfv <- delay(central_tfv, tau_ecl) / vc_tfv
    clag_ftc <- delay(central_ftc, tau_ecl) / vc_ftc
    c_efv_lag <- (t >= tau_ecl) * clag_efv + (t < tau_ecl) * start_suppressed * Cc
    c_tfv_lag <- (t >= tau_ecl) * clag_tfv + (t < tau_ecl) * start_suppressed * Cc_tfv
    c_ftc_lag <- (t >= tau_ecl) * clag_ftc + (t < tau_ecl) * start_suppressed * Cc_ftc
    ylag_wt <- delay(inf_wt, tau_ecl)
    ylag_k103n <- delay(inf_k103n, tau_ecl)
    ylag_y181c <- delay(inf_y181c, tau_ecl)
    ylag_g190s <- delay(inf_g190s, tau_ecl)
    ylag_m184v <- delay(inf_m184v, tau_ecl)
    ylag_k65r <- delay(inf_k65r, tau_ecl)
    ylag_k103n_y181c <- delay(inf_k103n_y181c, tau_ecl)
    ylag_k103n_g190s <- delay(inf_k103n_g190s, tau_ecl)
    ylag_k103n_m184v <- delay(inf_k103n_m184v, tau_ecl)
    ylag_k103n_k65r <- delay(inf_k103n_k65r, tau_ecl)
    ylag_y181c_g190s <- delay(inf_y181c_g190s, tau_ecl)
    ylag_y181c_m184v <- delay(inf_y181c_m184v, tau_ecl)
    ylag_y181c_k65r <- delay(inf_y181c_k65r, tau_ecl)
    ylag_g190s_m184v <- delay(inf_g190s_m184v, tau_ecl)
    ylag_g190s_k65r <- delay(inf_g190s_k65r, tau_ecl)
    ylag_m184v_k65r <- delay(inf_m184v_k65r, tau_ecl)

    # ---- Drug effect at the infection time t - tau ----
    # strain wt
    n1 <- hill_efv * (1 + 0)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_wt <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_wt <- expit(-uc)
    }
    rx_wt <- r0 * fit_wt * ufr_wt
    # strain k103n
    n1 <- hill_efv * (1 + (sigma_k103n_efv) / 1)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_k103n_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_k103n <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_k103n <- expit(-uc)
    }
    rx_k103n <- r0 * fit_k103n * ufr_k103n
    # strain y181c
    n1 <- hill_efv * (1 + (sigma_y181c_efv) / 1)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_y181c <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_y181c <- expit(-uc)
    }
    rx_y181c <- r0 * fit_y181c * ufr_y181c
    # strain g190s
    n1 <- hill_efv * (1 + (sigma_g190s_efv) / 1)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_g190s <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_g190s <- expit(-uc)
    }
    rx_g190s <- r0 * fit_g190s * ufr_g190s
    # strain m184v
    n1 <- hill_efv * (1 + 0)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 1)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 1)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_m184v <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_m184v <- expit(-uc)
    }
    rx_m184v <- r0 * fit_m184v * ufr_m184v
    # strain k65r
    n1 <- hill_efv * (1 + 0)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 1)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 1)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_k65r <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_k65r <- expit(-uc)
    }
    rx_k65r <- r0 * fit_k65r * ufr_k65r
    # strain k103n_y181c
    n1 <- hill_efv * (1 + (sigma_k103n_efv + sigma_y181c_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_k103n_efv * rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_k103n_y181c <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_k103n_y181c <- expit(-uc)
    }
    rx_k103n_y181c <- r0 * fit_k103n_y181c * ufr_k103n_y181c
    # strain k103n_g190s
    n1 <- hill_efv * (1 + (sigma_k103n_efv + sigma_g190s_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_k103n_efv * rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_k103n_g190s <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_k103n_g190s <- expit(-uc)
    }
    rx_k103n_g190s <- r0 * fit_k103n_g190s * ufr_k103n_g190s
    # strain k103n_m184v
    n1 <- hill_efv * (1 + (sigma_k103n_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_k103n_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_k103n_m184v <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_k103n_m184v <- expit(-uc)
    }
    rx_k103n_m184v <- r0 * fit_k103n_m184v * ufr_k103n_m184v
    # strain k103n_k65r
    n1 <- hill_efv * (1 + (sigma_k103n_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_k103n_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_k103n_k65r <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_k103n_k65r <- expit(-uc)
    }
    rx_k103n_k65r <- r0 * fit_k103n_k65r * ufr_k103n_k65r
    # strain y181c_g190s
    n1 <- hill_efv * (1 + (sigma_y181c_efv + sigma_g190s_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_y181c_efv * rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_y181c_g190s <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_y181c_g190s <- expit(-uc)
    }
    rx_y181c_g190s <- r0 * fit_y181c_g190s * ufr_y181c_g190s
    # strain y181c_m184v
    n1 <- hill_efv * (1 + (sigma_y181c_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_y181c_m184v <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_y181c_m184v <- expit(-uc)
    }
    rx_y181c_m184v <- r0 * fit_y181c_m184v * ufr_y181c_m184v
    # strain y181c_k65r
    n1 <- hill_efv * (1 + (sigma_y181c_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_y181c_k65r <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_y181c_k65r <- expit(-uc)
    }
    rx_y181c_k65r <- r0 * fit_y181c_k65r * ufr_y181c_k65r
    # strain g190s_m184v
    n1 <- hill_efv * (1 + (sigma_g190s_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_g190s_m184v <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_g190s_m184v <- expit(-uc)
    }
    rx_g190s_m184v <- r0 * fit_g190s_m184v * ufr_g190s_m184v
    # strain g190s_k65r
    n1 <- hill_efv * (1 + (sigma_g190s_efv) / 2)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_g190s_k65r <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_g190s_k65r <- expit(-uc)
    }
    rx_g190s_k65r <- r0 * fit_g190s_k65r * ufr_g190s_k65r
    # strain m184v_k65r
    n1 <- hill_efv * (1 + 0)
    if (c_efv_lag > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_lag) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv + sigma_k65r_tfv) / 2)
    if (c_tfv_lag > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_lag) - log(rho_m184v_tfv * rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc + sigma_k65r_ftc) / 2)
    if (c_ftc_lag > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_lag) - log(rho_m184v_ftc * rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
    if (p1 + p2 + p3 < 0.5) {
      ufr_m184v_k65r <- 1
    } else {
      # Bliss independence (product of single-drug unaffected fractions), on the log-odds scale
      sb <- 0
      if (p1 > 0.5) { if (z1 > 30) { sb <- sb + z1 } else { sb <- sb + log1p(exp(z1)) } }
      if (p2 > 0.5) { if (z2 > 30) { sb <- sb + z2 } else { sb <- sb + log1p(exp(z2)) } }
      if (p3 > 0.5) { if (z3 > 30) { sb <- sb + z3 } else { sb <- sb + log1p(exp(z3)) } }
      if (sb > 30) { ub <- sb } else { ub <- log(expm1(sb)) }
      # Loewe additivity: root of sum_j exp((z_j - x) / n_j) = 1 (Newton from the largest single-drug root)
      xl <- -1e10
      if (p1 > 0.5) { xl <- max(xl, z1) }
      if (p2 > 0.5) { xl <- max(xl, z2) }
      if (p3 > 0.5) { xl <- max(xl, z3) }
      dxl <- 1
      itl <- 0
      while (abs(dxl) > 1e-10 * (1 + abs(xl)) && itl < 100) {
        e1 <- 0; e2 <- 0; e3 <- 0
        if (p1 > 0.5) { e1 <- exp((z1 - xl) / n1) }
        if (p2 > 0.5) { e2 <- exp((z2 - xl) / n2) }
        if (p3 > 0.5) { e3 <- exp((z3 - xl) / n3) }
        dxl <- (e1 + e2 + e3 - 1) / (e1 / n1 + e2 / n2 + e3 / n3)
        xl <- xl + dxl
        itl <- itl + 1
      }
      # degree-of-independence weights from single-drug affected fractions (Jilek 2012 DI pairs)
      f1 <- 0; f2 <- 0; f3 <- 0
      if (p1 > 0.5) { f1 <- expit(z1) }
      if (p2 > 0.5) { f2 <- expit(z2) }
      if (p3 > 0.5) { f3 <- expit(z3) }
      wsum <- f1 * f2 + f1 * f3 + f2 * f3
      if (wsum > 0) { dicomb <- (f1 * f2 * di_efv_tfv + f1 * f3 * di_efv_ftc + f2 * f3 * di_tfv_ftc) / wsum } else { dicomb <- 0 }
      uc <- xl + dicomb * (ub - xl)
      ufr_m184v_k65r <- expit(-uc)
    }
    rx_m184v_k65r <- r0 * fit_m184v_k65r * ufr_m184v_k65r

    # ---- Productively infected cells (supplement boxed delay equation + de novo mutation transfer) ----
    yden <- lambda + dy * (rx_wt * ylag_wt + rx_k103n * ylag_k103n + rx_y181c * ylag_y181c + rx_g190s * ylag_g190s + rx_m184v * ylag_m184v + rx_k65r * ylag_k65r + rx_k103n_y181c * ylag_k103n_y181c + rx_k103n_g190s * ylag_k103n_g190s + rx_k103n_m184v * ylag_k103n_m184v + rx_k103n_k65r * ylag_k103n_k65r + rx_y181c_g190s * ylag_y181c_g190s + rx_y181c_m184v * ylag_y181c_m184v + rx_y181c_k65r * ylag_y181c_k65r + rx_g190s_m184v * ylag_g190s_m184v + rx_g190s_k65r * ylag_g190s_k65r + rx_m184v_k65r * ylag_m184v_k65r)
    ninf_wt <- dy * lambda * rx_wt * ylag_wt / yden
    ninf_k103n <- dy * lambda * rx_k103n * ylag_k103n / yden
    ninf_y181c <- dy * lambda * rx_y181c * ylag_y181c / yden
    ninf_g190s <- dy * lambda * rx_g190s * ylag_g190s / yden
    ninf_m184v <- dy * lambda * rx_m184v * ylag_m184v / yden
    ninf_k65r <- dy * lambda * rx_k65r * ylag_k65r / yden
    ninf_k103n_y181c <- dy * lambda * rx_k103n_y181c * ylag_k103n_y181c / yden
    ninf_k103n_g190s <- dy * lambda * rx_k103n_g190s * ylag_k103n_g190s / yden
    ninf_k103n_m184v <- dy * lambda * rx_k103n_m184v * ylag_k103n_m184v / yden
    ninf_k103n_k65r <- dy * lambda * rx_k103n_k65r * ylag_k103n_k65r / yden
    ninf_y181c_g190s <- dy * lambda * rx_y181c_g190s * ylag_y181c_g190s / yden
    ninf_y181c_m184v <- dy * lambda * rx_y181c_m184v * ylag_y181c_m184v / yden
    ninf_y181c_k65r <- dy * lambda * rx_y181c_k65r * ylag_y181c_k65r / yden
    ninf_g190s_m184v <- dy * lambda * rx_g190s_m184v * ylag_g190s_m184v / yden
    ninf_g190s_k65r <- dy * lambda * rx_g190s_k65r * ylag_g190s_k65r / yden
    ninf_m184v_k65r <- dy * lambda * rx_m184v_k65r * ylag_m184v_k65r / yden
    d/dt(inf_wt) <- atot * f_wt - dy * inf_wt + ninf_wt * (1 - (u_k103n + u_y181c + u_g190s + u_m184v + u_k65r))
    d/dt(inf_k103n) <- atot * f_k103n - dy * inf_k103n + ninf_k103n * (1 - (u_y181c + u_g190s + u_m184v + u_k65r)) + u_k103n * ninf_wt
    d/dt(inf_y181c) <- atot * f_y181c - dy * inf_y181c + ninf_y181c * (1 - (u_k103n + u_g190s + u_m184v + u_k65r)) + u_y181c * ninf_wt
    d/dt(inf_g190s) <- atot * f_g190s - dy * inf_g190s + ninf_g190s * (1 - (u_k103n + u_y181c + u_m184v + u_k65r)) + u_g190s * ninf_wt
    d/dt(inf_m184v) <- atot * f_m184v - dy * inf_m184v + ninf_m184v * (1 - (u_k103n + u_y181c + u_g190s + u_k65r)) + u_m184v * ninf_wt
    d/dt(inf_k65r) <- atot * f_k65r - dy * inf_k65r + ninf_k65r * (1 - (u_k103n + u_y181c + u_g190s + u_m184v)) + u_k65r * ninf_wt
    d/dt(inf_k103n_y181c) <- atot * f_k103n_y181c - dy * inf_k103n_y181c + ninf_k103n_y181c * (1 - (0)) + u_y181c * ninf_k103n + u_k103n * ninf_y181c
    d/dt(inf_k103n_g190s) <- atot * f_k103n_g190s - dy * inf_k103n_g190s + ninf_k103n_g190s * (1 - (0)) + u_g190s * ninf_k103n + u_k103n * ninf_g190s
    d/dt(inf_k103n_m184v) <- atot * f_k103n_m184v - dy * inf_k103n_m184v + ninf_k103n_m184v * (1 - (0)) + u_m184v * ninf_k103n + u_k103n * ninf_m184v
    d/dt(inf_k103n_k65r) <- atot * f_k103n_k65r - dy * inf_k103n_k65r + ninf_k103n_k65r * (1 - (0)) + u_k65r * ninf_k103n + u_k103n * ninf_k65r
    d/dt(inf_y181c_g190s) <- atot * f_y181c_g190s - dy * inf_y181c_g190s + ninf_y181c_g190s * (1 - (0)) + u_g190s * ninf_y181c + u_y181c * ninf_g190s
    d/dt(inf_y181c_m184v) <- atot * f_y181c_m184v - dy * inf_y181c_m184v + ninf_y181c_m184v * (1 - (0)) + u_m184v * ninf_y181c + u_y181c * ninf_m184v
    d/dt(inf_y181c_k65r) <- atot * f_y181c_k65r - dy * inf_y181c_k65r + ninf_y181c_k65r * (1 - (0)) + u_k65r * ninf_y181c + u_y181c * ninf_k65r
    d/dt(inf_g190s_m184v) <- atot * f_g190s_m184v - dy * inf_g190s_m184v + ninf_g190s_m184v * (1 - (0)) + u_m184v * ninf_g190s + u_g190s * ninf_m184v
    d/dt(inf_g190s_k65r) <- atot * f_g190s_k65r - dy * inf_g190s_k65r + ninf_g190s_k65r * (1 - (0)) + u_k65r * ninf_g190s + u_g190s * ninf_k65r
    d/dt(inf_m184v_k65r) <- atot * f_m184v_k65r - dy * inf_m184v_k65r + ninf_m184v_k65r * (1 - (0)) + u_k65r * ninf_m184v + u_m184v * ninf_k65r
    inf_wt(0) <- y0tot * f_wt
    inf_k103n(0) <- y0tot * f_k103n
    inf_y181c(0) <- y0tot * f_y181c
    inf_g190s(0) <- y0tot * f_g190s
    inf_m184v(0) <- y0tot * f_m184v
    inf_k65r(0) <- y0tot * f_k65r
    inf_k103n_y181c(0) <- y0tot * f_k103n_y181c
    inf_k103n_g190s(0) <- y0tot * f_k103n_g190s
    inf_k103n_m184v(0) <- y0tot * f_k103n_m184v
    inf_k103n_k65r(0) <- y0tot * f_k103n_k65r
    inf_y181c_g190s(0) <- y0tot * f_y181c_g190s
    inf_y181c_m184v(0) <- y0tot * f_y181c_m184v
    inf_y181c_k65r(0) <- y0tot * f_y181c_k65r
    inf_g190s_m184v(0) <- y0tot * f_g190s_m184v
    inf_g190s_k65r(0) <- y0tot * f_g190s_k65r
    inf_m184v_k65r(0) <- y0tot * f_m184v_k65r

    # ---- Outputs ----
    ytot <- inf_wt + inf_k103n + inf_y181c + inf_g190s + inf_m184v + inf_k65r + inf_k103n_y181c + inf_k103n_g190s + inf_k103n_m184v + inf_k103n_k65r + inf_y181c_g190s + inf_y181c_m184v + inf_y181c_k65r + inf_g190s_m184v + inf_g190s_k65r + inf_m184v_k65r
    VL <- ytot / cells_per_copy
    fmut <- 1 - inf_wt / ytot
  })
}
