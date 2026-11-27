Sanche_2017_efavirenz_tenofovir_emtricitabine_qsp_tauleap <- function() {
  description <- "QSP. Stochastic tau-leap HIV-1 viral-dynamics and resistance model for efavirenz monotherapy or efavirenz + tenofovir disoproxil fumarate + emtricitabine, as run by Sanche 2017 (Poisson counts of reservoir activations, infected-cell deaths, new infections and de novo mutations each leap, 0.1-h step, 1-day eclipse period). Same PK, drug-effect and sixteen-strain structure as Sanche_2017_efavirenz_tenofovir_emtricitabine_qsp, whose delay ODE is this process's expected value. Random numbers come from an MRG32k3a generator carried in six model states (mrg1..mrg6) and seeded per subject through mrg_seed1..mrg_seed6 (independent L'Ecuyer-CMRG streams, e.g. from parallel::nextRNGStream); rxode2's own random draws cannot be used because they are not drawn while the dynamics are solved. The model schedules its own leaps, so it adds one output row per leap; solve with rxSolve(..., maxExtra = 0L, method = 'lsoda') and set tend to the trial length (with the default liblsoda the step count accumulates across the model-scheduled intervals until maxsteps stops the solve). Counts with a mean above 100 use a normal approximation. Simulation only."
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

  # Productively infected CD4+ T cells per viral strain (wild type, five single
  # mutants, ten double mutants) have no canonical compartment analogue.
  # mrg1..mrg6 hold the MRG32k3a random-number generator and kleap the leap
  # counter; their derivatives are zero and only in-model replace() events move them.
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
    "inf_m184v_k65r",
    "mrg1",
    "mrg2",
    "mrg3",
    "mrg4",
    "mrg5",
    "mrg6",
    "kleap"
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
    ),
    mrg1 = list(
      analyte = "MRG32k3a generator state (not a drug or cell amount)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    mrg2 = list(
      analyte = "MRG32k3a generator state (not a drug or cell amount)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    mrg3 = list(
      analyte = "MRG32k3a generator state (not a drug or cell amount)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    mrg4 = list(
      analyte = "MRG32k3a generator state (not a drug or cell amount)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    mrg5 = list(
      analyte = "MRG32k3a generator state (not a drug or cell amount)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    mrg6 = list(
      analyte = "MRG32k3a generator state (not a drug or cell amount)",
      units = "unitless",
      specimen = "not applicable",
      verified = TRUE
    ),
    kleap = list(
      analyte = "tau-leap counter (not a drug or cell amount)",
      units = "count",
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
    dt_leap <- fixed(1 / 240); label("Tau-leap step (day)")  # Supplement 'Algorithm': Delta t = 0.1 hour
    tend <- fixed(336); label("Last time at which leaps and delayed arrivals are scheduled (day)")  # Main text: 48-week trials
    mrg_seed1 <- fixed(12345); label("MRG32k3a generator seed 1 (set per subject)")
    mrg_seed2 <- fixed(12345); label("MRG32k3a generator seed 2 (set per subject)")
    mrg_seed3 <- fixed(12345); label("MRG32k3a generator seed 3 (set per subject)")
    mrg_seed4 <- fixed(12345); label("MRG32k3a generator seed 4 (set per subject)")
    mrg_seed5 <- fixed(12345); label("MRG32k3a generator seed 5 (set per subject)")
    mrg_seed6 <- fixed(12345); label("MRG32k3a generator seed 6 (set per subject)")
    
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

    # ---- Stochastic state: infected cells per strain, MRG32k3a generator, leap counter ----
    d/dt(inf_wt) <- 0
    d/dt(inf_k103n) <- 0
    d/dt(inf_y181c) <- 0
    d/dt(inf_g190s) <- 0
    d/dt(inf_m184v) <- 0
    d/dt(inf_k65r) <- 0
    d/dt(inf_k103n_y181c) <- 0
    d/dt(inf_k103n_g190s) <- 0
    d/dt(inf_k103n_m184v) <- 0
    d/dt(inf_k103n_k65r) <- 0
    d/dt(inf_y181c_g190s) <- 0
    d/dt(inf_y181c_m184v) <- 0
    d/dt(inf_y181c_k65r) <- 0
    d/dt(inf_g190s_m184v) <- 0
    d/dt(inf_g190s_k65r) <- 0
    d/dt(inf_m184v_k65r) <- 0
    d/dt(mrg1) <- 0
    d/dt(mrg2) <- 0
    d/dt(mrg3) <- 0
    d/dt(mrg4) <- 0
    d/dt(mrg5) <- 0
    d/dt(mrg6) <- 0
    d/dt(kleap) <- 0
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
    mrg1(0) <- mrg_seed1
    mrg2(0) <- mrg_seed2
    mrg3(0) <- mrg_seed3
    mrg4(0) <- mrg_seed4
    mrg5(0) <- mrg_seed5
    mrg6(0) <- mrg_seed6
    kleap(0) <- 0
    # leap k happens at t = k / nld exactly, so leap times coincide exactly with whole-day
    # dosing records instead of drifting from them by rounding error
    nld <- floor(1 / dt_leap + 0.5)
    ntau <- floor(tau_ecl * nld + 0.5)

    # ---- Drug effect at the current time ----
    c_efv_cur <- Cc
    c_tfv_cur <- Cc_tfv
    c_ftc_cur <- Cc_ftc
    # strain wt
    n1 <- hill_efv * (1 + 0)
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_k103n_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 1)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 1)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 1)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 1)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_k103n_efv * rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_k103n_efv * rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_k103n_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_k103n_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_y181c_efv * rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + 0)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(1) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + 0)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(1) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_y181c_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_m184v_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_m184v_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(rho_g190s_efv) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_k65r_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_k65r_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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
    if (c_efv_cur > 0) { p1 <- 1; z1 <- n1 * (log(c_efv_cur) - log(1) - lic50i_efv - log(kp_efv) - log(kl_efv)) } else { p1 <- 0; z1 <- 0 }
    n2 <- hill_tfv * (1 + (sigma_m184v_tfv + sigma_k65r_tfv) / 2)
    if (c_tfv_cur > 0) { p2 <- 1; z2 <- n2 * (log(c_tfv_cur) - log(rho_m184v_tfv * rho_k65r_tfv) - lic50i_tfv - log(kp_tfv) - log(kl_tfv)) } else { p2 <- 0; z2 <- 0 }
    n3 <- hill_ftc * (1 + (sigma_m184v_ftc + sigma_k65r_ftc) / 2)
    if (c_ftc_cur > 0) { p3 <- 1; z3 <- n3 * (log(c_ftc_cur) - log(rho_m184v_ftc * rho_k65r_ftc) - lic50i_ftc - log(kp_ftc) - log(kl_ftc)) } else { p3 <- 0; z3 <- 0 }
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

    # ---- Outputs ----
    ytot <- inf_wt + inf_k103n + inf_y181c + inf_g190s + inf_m184v + inf_k65r + inf_k103n_y181c + inf_k103n_g190s + inf_k103n_m184v + inf_k103n_k65r + inf_y181c_g190s + inf_y181c_m184v + inf_y181c_k65r + inf_g190s_m184v + inf_g190s_k65r + inf_m184v_k65r
    VL <- ytot / cells_per_copy
    fmut <- 1 - inf_wt / ytot

    # ---- Tau-leap step (Supplement 'Simulating viral kinetics'), run once per leap time ----
    # The model schedules its own leaps (obs()); every count is drawn from the
    # in-model generator held in mrg1..mrg6. A leap at time t updates the state from
    # y(t) to y(t + dt_leap): reservoir activations and deaths are applied at t, and the
    # new infections drawn at t are delivered tau later (between the leaps at t + tau and
    # t + tau + dt_leap), when those cells start producing virus.
    if (t >= kleap / nld - 1e-9 && (kleap + 1) / nld <= tend + 1e-9) {
      r10 <- mrg1; r11 <- mrg2; r12 <- mrg3; r20 <- mrg4; r21 <- mrg5; r22 <- mrg6
      ycur_wt <- inf_wt
      ycur_k103n <- inf_k103n
      ycur_y181c <- inf_y181c
      ycur_g190s <- inf_g190s
      ycur_m184v <- inf_m184v
      ycur_k65r <- inf_k65r
      ycur_k103n_y181c <- inf_k103n_y181c
      ycur_k103n_g190s <- inf_k103n_g190s
      ycur_k103n_m184v <- inf_k103n_m184v
      ycur_k103n_k65r <- inf_k103n_k65r
      ycur_y181c_g190s <- inf_y181c_g190s
      ycur_y181c_m184v <- inf_y181c_m184v
      ycur_y181c_k65r <- inf_y181c_k65r
      ycur_g190s_m184v <- inf_g190s_m184v
      ycur_g190s_k65r <- inf_g190s_k65r
      ycur_m184v_k65r <- inf_m184v_k65r
      npre <- 0
      if (kleap < 0.5) {
        # first leap: Poisson initial composition, then the infections generated during the
        # eclipse period before t = 0 by that constant pre-trial state
        mu <- y0tot * f_wt
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_wt <- nn
        mu <- y0tot * f_k103n
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_k103n <- nn
        mu <- y0tot * f_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_y181c <- nn
        mu <- y0tot * f_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_g190s <- nn
        mu <- y0tot * f_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_m184v <- nn
        mu <- y0tot * f_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_k65r <- nn
        mu <- y0tot * f_k103n_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_k103n_y181c <- nn
        mu <- y0tot * f_k103n_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_k103n_g190s <- nn
        mu <- y0tot * f_k103n_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_k103n_m184v <- nn
        mu <- y0tot * f_k103n_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_k103n_k65r <- nn
        mu <- y0tot * f_y181c_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_y181c_g190s <- nn
        mu <- y0tot * f_y181c_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_y181c_m184v <- nn
        mu <- y0tot * f_y181c_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_y181c_k65r <- nn
        mu <- y0tot * f_g190s_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_g190s_m184v <- nn
        mu <- y0tot * f_g190s_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_g190s_k65r <- nn
        mu <- y0tot * f_m184v_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ycur_m184v_k65r <- nn
        npre <- ntau
      }
      # latent-reservoir activations, then deaths of productively infected cells
      mu <- atot * f_wt * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_wt <- nn
      mu <- atot * f_k103n * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_k103n <- nn
      mu <- atot * f_y181c * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_y181c <- nn
      mu <- atot * f_g190s * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_g190s <- nn
      mu <- atot * f_m184v * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_m184v <- nn
      mu <- atot * f_k65r * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_k65r <- nn
      mu <- atot * f_k103n_y181c * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_k103n_y181c <- nn
      mu <- atot * f_k103n_g190s * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_k103n_g190s <- nn
      mu <- atot * f_k103n_m184v * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_k103n_m184v <- nn
      mu <- atot * f_k103n_k65r * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_k103n_k65r <- nn
      mu <- atot * f_y181c_g190s * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_y181c_g190s <- nn
      mu <- atot * f_y181c_m184v * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_y181c_m184v <- nn
      mu <- atot * f_y181c_k65r * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_y181c_k65r <- nn
      mu <- atot * f_g190s_m184v * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_g190s_m184v <- nn
      mu <- atot * f_g190s_k65r * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_g190s_k65r <- nn
      mu <- atot * f_m184v_k65r * dt_leap
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      nres_m184v_k65r <- nn
      mu <- ycur_wt * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_wt <- nn
      ndth_wt <- min(ndth_wt, ycur_wt)
      mu <- ycur_k103n * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_k103n <- nn
      ndth_k103n <- min(ndth_k103n, ycur_k103n)
      mu <- ycur_y181c * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_y181c <- nn
      ndth_y181c <- min(ndth_y181c, ycur_y181c)
      mu <- ycur_g190s * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_g190s <- nn
      ndth_g190s <- min(ndth_g190s, ycur_g190s)
      mu <- ycur_m184v * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_m184v <- nn
      ndth_m184v <- min(ndth_m184v, ycur_m184v)
      mu <- ycur_k65r * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_k65r <- nn
      ndth_k65r <- min(ndth_k65r, ycur_k65r)
      mu <- ycur_k103n_y181c * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_k103n_y181c <- nn
      ndth_k103n_y181c <- min(ndth_k103n_y181c, ycur_k103n_y181c)
      mu <- ycur_k103n_g190s * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_k103n_g190s <- nn
      ndth_k103n_g190s <- min(ndth_k103n_g190s, ycur_k103n_g190s)
      mu <- ycur_k103n_m184v * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_k103n_m184v <- nn
      ndth_k103n_m184v <- min(ndth_k103n_m184v, ycur_k103n_m184v)
      mu <- ycur_k103n_k65r * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_k103n_k65r <- nn
      ndth_k103n_k65r <- min(ndth_k103n_k65r, ycur_k103n_k65r)
      mu <- ycur_y181c_g190s * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_y181c_g190s <- nn
      ndth_y181c_g190s <- min(ndth_y181c_g190s, ycur_y181c_g190s)
      mu <- ycur_y181c_m184v * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_y181c_m184v <- nn
      ndth_y181c_m184v <- min(ndth_y181c_m184v, ycur_y181c_m184v)
      mu <- ycur_y181c_k65r * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_y181c_k65r <- nn
      ndth_y181c_k65r <- min(ndth_y181c_k65r, ycur_y181c_k65r)
      mu <- ycur_g190s_m184v * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_g190s_m184v <- nn
      ndth_g190s_m184v <- min(ndth_g190s_m184v, ycur_g190s_m184v)
      mu <- ycur_g190s_k65r * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_g190s_k65r <- nn
      ndth_g190s_k65r <- min(ndth_g190s_k65r, ycur_g190s_k65r)
      mu <- ycur_m184v_k65r * (1 - exp(-dy * dt_leap))
      p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
      p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
      if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
      nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
      ndth_m184v_k65r <- nn
      ndth_m184v_k65r <- min(ndth_m184v_k65r, ycur_m184v_k65r)
      # new productive infections and de novo mutations: npre pre-trial passes, then this leap
      ipass <- 0
      while (ipass <= npre) {
        if (ipass < npre) {
          rsrc_wt <- (1 - start_suppressed) * r0 * fit_wt + start_suppressed * rx_wt
          rsrc_k103n <- (1 - start_suppressed) * r0 * fit_k103n + start_suppressed * rx_k103n
          rsrc_y181c <- (1 - start_suppressed) * r0 * fit_y181c + start_suppressed * rx_y181c
          rsrc_g190s <- (1 - start_suppressed) * r0 * fit_g190s + start_suppressed * rx_g190s
          rsrc_m184v <- (1 - start_suppressed) * r0 * fit_m184v + start_suppressed * rx_m184v
          rsrc_k65r <- (1 - start_suppressed) * r0 * fit_k65r + start_suppressed * rx_k65r
          rsrc_k103n_y181c <- (1 - start_suppressed) * r0 * fit_k103n_y181c + start_suppressed * rx_k103n_y181c
          rsrc_k103n_g190s <- (1 - start_suppressed) * r0 * fit_k103n_g190s + start_suppressed * rx_k103n_g190s
          rsrc_k103n_m184v <- (1 - start_suppressed) * r0 * fit_k103n_m184v + start_suppressed * rx_k103n_m184v
          rsrc_k103n_k65r <- (1 - start_suppressed) * r0 * fit_k103n_k65r + start_suppressed * rx_k103n_k65r
          rsrc_y181c_g190s <- (1 - start_suppressed) * r0 * fit_y181c_g190s + start_suppressed * rx_y181c_g190s
          rsrc_y181c_m184v <- (1 - start_suppressed) * r0 * fit_y181c_m184v + start_suppressed * rx_y181c_m184v
          rsrc_y181c_k65r <- (1 - start_suppressed) * r0 * fit_y181c_k65r + start_suppressed * rx_y181c_k65r
          rsrc_g190s_m184v <- (1 - start_suppressed) * r0 * fit_g190s_m184v + start_suppressed * rx_g190s_m184v
          rsrc_g190s_k65r <- (1 - start_suppressed) * r0 * fit_g190s_k65r + start_suppressed * rx_g190s_k65r
          rsrc_m184v_k65r <- (1 - start_suppressed) * r0 * fit_m184v_k65r + start_suppressed * rx_m184v_k65r
          tarr <- (ipass + 0.5) / nld
        } else {
          rsrc_wt <- rx_wt
          rsrc_k103n <- rx_k103n
          rsrc_y181c <- rx_y181c
          rsrc_g190s <- rx_g190s
          rsrc_m184v <- rx_m184v
          rsrc_k65r <- rx_k65r
          rsrc_k103n_y181c <- rx_k103n_y181c
          rsrc_k103n_g190s <- rx_k103n_g190s
          rsrc_k103n_m184v <- rx_k103n_m184v
          rsrc_k103n_k65r <- rx_k103n_k65r
          rsrc_y181c_g190s <- rx_y181c_g190s
          rsrc_y181c_m184v <- rx_y181c_m184v
          rsrc_y181c_k65r <- rx_y181c_k65r
          rsrc_g190s_m184v <- rx_g190s_m184v
          rsrc_g190s_k65r <- rx_g190s_k65r
          rsrc_m184v_k65r <- rx_m184v_k65r
          tarr <- (kleap + ntau + 0.5) / nld
        }
        sden <- lambda + dy * (rsrc_wt * ycur_wt + rsrc_k103n * ycur_k103n + rsrc_y181c * ycur_y181c + rsrc_g190s * ycur_g190s + rsrc_m184v * ycur_m184v + rsrc_k65r * ycur_k65r + rsrc_k103n_y181c * ycur_k103n_y181c + rsrc_k103n_g190s * ycur_k103n_g190s + rsrc_k103n_m184v * ycur_k103n_m184v + rsrc_k103n_k65r * ycur_k103n_k65r + rsrc_y181c_g190s * ycur_y181c_g190s + rsrc_y181c_m184v * ycur_y181c_m184v + rsrc_y181c_k65r * ycur_y181c_k65r + rsrc_g190s_m184v * ycur_g190s_m184v + rsrc_g190s_k65r * ycur_g190s_k65r + rsrc_m184v_k65r * ycur_m184v_k65r)
        mu <- dt_leap * dy * lambda * rsrc_wt * ycur_wt / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_wt <- nn
        mu <- dt_leap * dy * lambda * rsrc_k103n * ycur_k103n / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_k103n <- nn
        mu <- dt_leap * dy * lambda * rsrc_y181c * ycur_y181c / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_y181c <- nn
        mu <- dt_leap * dy * lambda * rsrc_g190s * ycur_g190s / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_g190s <- nn
        mu <- dt_leap * dy * lambda * rsrc_m184v * ycur_m184v / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_m184v <- nn
        mu <- dt_leap * dy * lambda * rsrc_k65r * ycur_k65r / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_k65r <- nn
        mu <- dt_leap * dy * lambda * rsrc_k103n_y181c * ycur_k103n_y181c / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_k103n_y181c <- nn
        mu <- dt_leap * dy * lambda * rsrc_k103n_g190s * ycur_k103n_g190s / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_k103n_g190s <- nn
        mu <- dt_leap * dy * lambda * rsrc_k103n_m184v * ycur_k103n_m184v / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_k103n_m184v <- nn
        mu <- dt_leap * dy * lambda * rsrc_k103n_k65r * ycur_k103n_k65r / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_k103n_k65r <- nn
        mu <- dt_leap * dy * lambda * rsrc_y181c_g190s * ycur_y181c_g190s / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_y181c_g190s <- nn
        mu <- dt_leap * dy * lambda * rsrc_y181c_m184v * ycur_y181c_m184v / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_y181c_m184v <- nn
        mu <- dt_leap * dy * lambda * rsrc_y181c_k65r * ycur_y181c_k65r / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_y181c_k65r <- nn
        mu <- dt_leap * dy * lambda * rsrc_g190s_m184v * ycur_g190s_m184v / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_g190s_m184v <- nn
        mu <- dt_leap * dy * lambda * rsrc_g190s_k65r * ycur_g190s_k65r / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_g190s_k65r <- nn
        mu <- dt_leap * dy * lambda * rsrc_m184v_k65r * ycur_m184v_k65r / sden
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        ninf_m184v_k65r <- nn
        mu <- ninf_wt * u_k103n
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut1 <- nn
        mu <- ninf_wt * u_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut2 <- nn
        mu <- ninf_wt * u_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut3 <- nn
        mu <- ninf_wt * u_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut4 <- nn
        mu <- ninf_wt * u_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut5 <- nn
        mu <- ninf_k103n * u_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut6 <- nn
        mu <- ninf_k103n * u_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut7 <- nn
        mu <- ninf_k103n * u_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut8 <- nn
        mu <- ninf_k103n * u_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut9 <- nn
        mu <- ninf_y181c * u_k103n
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut10 <- nn
        mu <- ninf_y181c * u_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut11 <- nn
        mu <- ninf_y181c * u_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut12 <- nn
        mu <- ninf_y181c * u_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut13 <- nn
        mu <- ninf_g190s * u_k103n
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut14 <- nn
        mu <- ninf_g190s * u_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut15 <- nn
        mu <- ninf_g190s * u_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut16 <- nn
        mu <- ninf_g190s * u_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut17 <- nn
        mu <- ninf_m184v * u_k103n
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut18 <- nn
        mu <- ninf_m184v * u_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut19 <- nn
        mu <- ninf_m184v * u_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut20 <- nn
        mu <- ninf_m184v * u_k65r
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut21 <- nn
        mu <- ninf_k65r * u_k103n
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut22 <- nn
        mu <- ninf_k65r * u_y181c
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut23 <- nn
        mu <- ninf_k65r * u_g190s
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut24 <- nn
        mu <- ninf_k65r * u_m184v
        p1 <- 1403580 * r11 - 810728 * r10; p1 <- p1 - 4294967087 * floor(p1 / 4294967087); r10 <- r11; r11 <- r12; r12 <- p1
        p2 <- 527612 * r22 - 1370589 * r20; p2 <- p2 - 4294944443 * floor(p2 / 4294944443); r20 <- r21; r21 <- r22; r22 <- p2
        if (p1 > p2) { u <- (p1 - p2) * 2.328306549295727688e-10 } else { u <- (p1 - p2 + 4294967087) * 2.328306549295727688e-10 }
        nn <- 0; if (mu > 100) { nn <- max(0, floor(mu + sqrt(mu) * probit(u) + 0.5)) } else { if (mu > 0) { pk <- exp(-mu); cdf <- pk; while (u > cdf && nn < 1000) { nn <- nn + 1; pk <- pk * mu / nn; cdf <- cdf + pk } } }
        nmut25 <- nn
        narr_wt <- ninf_wt - min(nmut1 + nmut2 + nmut3 + nmut4 + nmut5, ninf_wt) + 0
        narr_k103n <- ninf_k103n - min(nmut6 + nmut7 + nmut8 + nmut9, ninf_k103n) + nmut1
        narr_y181c <- ninf_y181c - min(nmut10 + nmut11 + nmut12 + nmut13, ninf_y181c) + nmut2
        narr_g190s <- ninf_g190s - min(nmut14 + nmut15 + nmut16 + nmut17, ninf_g190s) + nmut3
        narr_m184v <- ninf_m184v - min(nmut18 + nmut19 + nmut20 + nmut21, ninf_m184v) + nmut4
        narr_k65r <- ninf_k65r - min(nmut22 + nmut23 + nmut24 + nmut25, ninf_k65r) + nmut5
        narr_k103n_y181c <- ninf_k103n_y181c - 0 + nmut6 + nmut10
        narr_k103n_g190s <- ninf_k103n_g190s - 0 + nmut7 + nmut14
        narr_k103n_m184v <- ninf_k103n_m184v - 0 + nmut8 + nmut18
        narr_k103n_k65r <- ninf_k103n_k65r - 0 + nmut9 + nmut22
        narr_y181c_g190s <- ninf_y181c_g190s - 0 + nmut11 + nmut15
        narr_y181c_m184v <- ninf_y181c_m184v - 0 + nmut12 + nmut19
        narr_y181c_k65r <- ninf_y181c_k65r - 0 + nmut13 + nmut23
        narr_g190s_m184v <- ninf_g190s_m184v - 0 + nmut16 + nmut20
        narr_g190s_k65r <- ninf_g190s_k65r - 0 + nmut17 + nmut24
        narr_m184v_k65r <- ninf_m184v_k65r - 0 + nmut21 + nmut25
        if (narr_wt > 0 && tarr <= tend) { evid_(tarr, 1, narr_wt, 9) }
        if (narr_k103n > 0 && tarr <= tend) { evid_(tarr, 1, narr_k103n, 10) }
        if (narr_y181c > 0 && tarr <= tend) { evid_(tarr, 1, narr_y181c, 11) }
        if (narr_g190s > 0 && tarr <= tend) { evid_(tarr, 1, narr_g190s, 12) }
        if (narr_m184v > 0 && tarr <= tend) { evid_(tarr, 1, narr_m184v, 13) }
        if (narr_k65r > 0 && tarr <= tend) { evid_(tarr, 1, narr_k65r, 14) }
        if (narr_k103n_y181c > 0 && tarr <= tend) { evid_(tarr, 1, narr_k103n_y181c, 15) }
        if (narr_k103n_g190s > 0 && tarr <= tend) { evid_(tarr, 1, narr_k103n_g190s, 16) }
        if (narr_k103n_m184v > 0 && tarr <= tend) { evid_(tarr, 1, narr_k103n_m184v, 17) }
        if (narr_k103n_k65r > 0 && tarr <= tend) { evid_(tarr, 1, narr_k103n_k65r, 18) }
        if (narr_y181c_g190s > 0 && tarr <= tend) { evid_(tarr, 1, narr_y181c_g190s, 19) }
        if (narr_y181c_m184v > 0 && tarr <= tend) { evid_(tarr, 1, narr_y181c_m184v, 20) }
        if (narr_y181c_k65r > 0 && tarr <= tend) { evid_(tarr, 1, narr_y181c_k65r, 21) }
        if (narr_g190s_m184v > 0 && tarr <= tend) { evid_(tarr, 1, narr_g190s_m184v, 22) }
        if (narr_g190s_k65r > 0 && tarr <= tend) { evid_(tarr, 1, narr_g190s_k65r, 23) }
        if (narr_m184v_k65r > 0 && tarr <= tend) { evid_(tarr, 1, narr_m184v_k65r, 24) }
        ipass <- ipass + 1
      }
      ynew_wt <- ycur_wt + nres_wt - ndth_wt
      ynew_k103n <- ycur_k103n + nres_k103n - ndth_k103n
      ynew_y181c <- ycur_y181c + nres_y181c - ndth_y181c
      ynew_g190s <- ycur_g190s + nres_g190s - ndth_g190s
      ynew_m184v <- ycur_m184v + nres_m184v - ndth_m184v
      ynew_k65r <- ycur_k65r + nres_k65r - ndth_k65r
      ynew_k103n_y181c <- ycur_k103n_y181c + nres_k103n_y181c - ndth_k103n_y181c
      ynew_k103n_g190s <- ycur_k103n_g190s + nres_k103n_g190s - ndth_k103n_g190s
      ynew_k103n_m184v <- ycur_k103n_m184v + nres_k103n_m184v - ndth_k103n_m184v
      ynew_k103n_k65r <- ycur_k103n_k65r + nres_k103n_k65r - ndth_k103n_k65r
      ynew_y181c_g190s <- ycur_y181c_g190s + nres_y181c_g190s - ndth_y181c_g190s
      ynew_y181c_m184v <- ycur_y181c_m184v + nres_y181c_m184v - ndth_y181c_m184v
      ynew_y181c_k65r <- ycur_y181c_k65r + nres_y181c_k65r - ndth_y181c_k65r
      ynew_g190s_m184v <- ycur_g190s_m184v + nres_g190s_m184v - ndth_g190s_m184v
      ynew_g190s_k65r <- ycur_g190s_k65r + nres_g190s_k65r - ndth_g190s_k65r
      ynew_m184v_k65r <- ycur_m184v_k65r + nres_m184v_k65r - ndth_m184v_k65r
      if (ynew_wt != inf_wt) { replace(ynew_wt, inf_wt) }
      if (ynew_k103n != inf_k103n) { replace(ynew_k103n, inf_k103n) }
      if (ynew_y181c != inf_y181c) { replace(ynew_y181c, inf_y181c) }
      if (ynew_g190s != inf_g190s) { replace(ynew_g190s, inf_g190s) }
      if (ynew_m184v != inf_m184v) { replace(ynew_m184v, inf_m184v) }
      if (ynew_k65r != inf_k65r) { replace(ynew_k65r, inf_k65r) }
      if (ynew_k103n_y181c != inf_k103n_y181c) { replace(ynew_k103n_y181c, inf_k103n_y181c) }
      if (ynew_k103n_g190s != inf_k103n_g190s) { replace(ynew_k103n_g190s, inf_k103n_g190s) }
      if (ynew_k103n_m184v != inf_k103n_m184v) { replace(ynew_k103n_m184v, inf_k103n_m184v) }
      if (ynew_k103n_k65r != inf_k103n_k65r) { replace(ynew_k103n_k65r, inf_k103n_k65r) }
      if (ynew_y181c_g190s != inf_y181c_g190s) { replace(ynew_y181c_g190s, inf_y181c_g190s) }
      if (ynew_y181c_m184v != inf_y181c_m184v) { replace(ynew_y181c_m184v, inf_y181c_m184v) }
      if (ynew_y181c_k65r != inf_y181c_k65r) { replace(ynew_y181c_k65r, inf_y181c_k65r) }
      if (ynew_g190s_m184v != inf_g190s_m184v) { replace(ynew_g190s_m184v, inf_g190s_m184v) }
      if (ynew_g190s_k65r != inf_g190s_k65r) { replace(ynew_g190s_k65r, inf_g190s_k65r) }
      if (ynew_m184v_k65r != inf_m184v_k65r) { replace(ynew_m184v_k65r, inf_m184v_k65r) }
      replace(r10, mrg1)
      replace(r11, mrg2)
      replace(r12, mrg3)
      replace(r20, mrg4)
      replace(r21, mrg5)
      replace(r22, mrg6)
      replace(kleap + 1, kleap)
      obs((kleap + 1) / nld - t)
    }
  })
}
