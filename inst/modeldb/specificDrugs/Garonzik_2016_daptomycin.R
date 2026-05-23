Garonzik_2016_daptomycin <- function() {
  description <- "In vitro (Staphylococcus aureus USA300, methicillin-resistant CA-MRSA reference strain). Mechanism-based mathematical pharmacodynamic (MBM) model of daptomycin time-kill activity in supplemented Mueller-Hinton broth with 0%, 10%, 30%, 50%, or 70% v/v heat-inactivated human serum. The bacterial population is split into three subpopulations of decreasing daptomycin susceptibility (susceptible, intermediate, resistant), each described by two states (state 1 vegetative, state 2 replicating; six bacterial compartments total). Replication of state 2 cells back into state 1 is gated by a successful-replication probability (REP = 2 x Plateau, with Plateau saturating at a maximum CFU/mL CFUm), and the vegetative-to-replicating transition k12 is modulated by an exponential lag-phase term (Eq 3) and a saturable carrying-capacity term (Eq 7) parameterised by Imax_k12 and IC50_k12. Daptomycin acts on each subpopulation via two mechanisms: a Hill-type stimulation of the probability of death (STI; reduces successful replication via IREP = 1 - STI) and a Hill-type direct killing of bacteria (Kill); the relative balance of the two is the dominant pharmacodynamic feature, with SC50 (0.05 mg/L) much lower than KC50 (4.8 mg/L). The intermediate and resistant subpopulations share the same SC50 and KC50 but have reduced Smax and Kmax (Smax_r and Kmax_r fixed to 0) and the resistant subpopulation has a slower vegetative-to-replicating transition (FR_K12r = 0.0442). Protein binding by human serum is encoded as an 'active fraction' factive(HS) multiplying the total static daptomycin concentration to give an effective drug concentration DAP_EF; the active fraction takes five experimental levels (factive = 1 at 0% HS, then 0.346, 0.284, 0.239, 0.252 at 10%, 30%, 50%, 70% HS). The model is in-vitro PD only -- there is no human PK component; daptomycin is dosed once at t = 0 into the dap compartment and is chemically stable in the medium for the 24-h experiment. Random effects (eta) are NOT present: the paper reports replicate-level experimental fits with additive plus small-count Poisson residual error on log10 CFU/mL."
  reference <- "Garonzik SM, Lenhard JR, Forrest A, Holden PN, Bulitta JB, Tsuji BT. Defining the active fraction of daptomycin against methicillin-resistant Staphylococcus aureus (MRSA) using a pharmacokinetic and pharmacodynamic approach. PLoS ONE. 2016;11(6):e0156131. doi:10.1371/journal.pone.0156131."
  vignette <- "Garonzik_2016_daptomycin"
  units <- list(
    time          = "hour",
    dosing        = "mg/L (static initial daptomycin concentration in the broth)",
    concentration = "log10 CFU/mL (observation); mg/L (daptomycin state)"
  )

  covariateData <- list(
    HS = list(
      description        = "Human serum percentage (v/v) supplementing Mueller-Hinton broth in the in-vitro time-kill experiment",
      units              = "% v/v",
      type               = "categorical",
      reference_category = "0 (no human serum; factive forced to 1)",
      notes              = paste(
        "Five experimental levels were studied: HS in {0, 10, 30, 50, 70}.",
        "Drives the active-fraction factive(HS) that multiplies the total",
        "static daptomycin concentration to give the effective concentration",
        "DAP_EF (Garonzik 2016 Eq 2). factive estimates from Table 2:",
        "factive(10) = 0.346, factive(30) = 0.284, factive(50) = 0.239,",
        "factive(70) = 0.252; factive(0) = 1 by construction (no protein",
        "binding in the absence of serum). HS values outside this discrete",
        "set make factive = 0 inside the model (none of the indicator",
        "expressions fire) -- a deliberately conspicuous failure rather",
        "than silent interpolation; if a user wants intermediate HS levels",
        "they must add the corresponding factive value or interpolate",
        "outside the model. In-vitro experimental indicator -- not in",
        "inst/references/covariate-columns.md (the canonical register is",
        "for human pop-PK covariates and does not apply to this in-vitro",
        "PD model)."
      ),
      source_name        = "% Human Serum (paper Methods + Table 2)"
    )
  )

  population <- list(
    species             = "in vitro (Staphylococcus aureus USA300 NRS384/FRP3757, methicillin-resistant CA-MRSA reference strain)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Staphylococcus aureus USA300 (NRS 384, FRP3757), daptomycin MIC = 0.5 mg/L. MBM was fit to USA300 only -- selected as the most common pulsed-field gel electrophoresis type in the USA.",
    second_strain       = "Time-kill data were also collected for VISA Mu50 (NRS 4, HIP5836), daptomycin MIC = 1.0 mg/L, but the MBM was NOT fit to Mu50 (Hill-only PD parameters for Mu50 in Table 1; this model file reproduces the MBM in Table 2 for USA300).",
    system              = "Static time-kill experiments, 24-hour exposure, dense sampling at 0, 1, 2, 4, 8, 24 h.",
    medium              = "Mueller-Hinton broth supplemented with calcium and magnesium (12.5 mg/L; supplemented Mueller-Hinton broth, SMHB) plus heat-inactivated human serum at the specified HS percentage. Calcium concentration in each batch titrated to physiologic conditions (1.1-1.3 mmol/L).",
    temperature         = "35 C (standard bacteriology incubation for S. aureus).",
    duration            = "24 h.",
    starting_inoculum   = "approximately 10^6 CFU/mL (model estimate Log10CFU0 = 6.22, Table 2).",
    concentration_range = "Daptomycin 0, 0.125 (USA300 only), 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128 mg/L (and 256 mg/L for USA300 per Fig 2 legend); the MBM in Table 2 was fit to these static-concentration time-kill data.",
    hs_options          = "0%, 10%, 30%, 50%, 70% v/v human serum / SMHB.",
    notes               = paste(
      "In-vitro pharmacodynamic study; no human or animal subjects. Limit",
      "of detection 10^2 CFU/mL (5 colonies per agar plate from an",
      "undiluted sample). Bactericidal activity defined as >=3.0 log10",
      "CFU/mL decrease from initial inoculum. Residual variability was",
      "modelled in NONMEM VI (FOCE) using additive plus Poisson error",
      "models per Bulitta et al. (paper reference 11); the single",
      "additive residual SD on log10 CFU/mL (epsilon_CFU = 0.558) is",
      "carried into this model file via addSd, with the small-count",
      "Poisson term and the additional additive component for counts < 5",
      "(both fixed to 1.00 and 0.250 respectively in Table 2) omitted",
      "as a deliberate simplification documented in the vignette",
      "Assumptions and deviations. The published r^2 was 0.95 for the",
      "USA300 MBM fit (Results, Mechanism based modeling)."
    )
  )

  ini({
    # ---- Active daptomycin fraction at each HS level (Garonzik 2016 Table 2) ----
    # Each factive is bounded in (0, 1); log-transform keeps positivity and is
    # consistent with the rest of the parameter set. factive(0%) = 1 by
    # construction (no serum, no protein binding) and is not estimated.
    lfact_hs10 <- log(0.346); label("Active daptomycin fraction at 10% HS (factive10, unitless)")  # Garonzik 2016 Table 2 (factive 10%)
    lfact_hs30 <- log(0.284); label("Active daptomycin fraction at 30% HS (factive30, unitless)")  # Garonzik 2016 Table 2 (factive 30%)
    lfact_hs50 <- log(0.239); label("Active daptomycin fraction at 50% HS (factive50, unitless)")  # Garonzik 2016 Table 2 (factive 50%)
    lfact_hs70 <- log(0.252); label("Active daptomycin fraction at 70% HS (factive70, unitless)")  # Garonzik 2016 Table 2 (factive 70%)

    # ---- Bacterial inoculum and subpopulation fractions ----
    # The paper parameterises CFU0, FR_I, FR_r on log10 scale (Table 2). These
    # parameters are kept on the same scale to make the round-trip from the
    # published values exact.
    cfu0_log10 <- 6.22;  label("Initial bacterial inoculum (log10 CFU/mL)")                          # Garonzik 2016 Table 2 (Log10CFUo)
    fri_log10  <- -3.65; label("Intermediate population as a fraction of initial inoculum (log10)") # Garonzik 2016 Table 2 (Log10 FR_I)
    frr_log10  <- -5.67; label("Resistant population as a fraction of initial inoculum (log10)")    # Garonzik 2016 Table 2 (Log10 FR_r)

    # ---- Bacterial life-cycle structural parameters ----
    lmtt_lag    <- log(75.5);                label("Mean transit time for lag phase (MTT_lag, h)")                       # Garonzik 2016 Table 2 (MTT_lag)
    beta_lag    <- fixed(10.0);              label("Sigmoidicity constant for lag-phase function (beta, unitless; FIXED)") # Garonzik 2016 Table 2 (beta; reported with 0% SE = FIXED)
    lmtt_k12    <- log(20.2);                label("Mean transit time S1 -> S2 (MTT_K12, h; k12 = 1/MTT_K12)")          # Garonzik 2016 Table 2 (MTT_K12)
    ic50_k12_log10 <- 7.81;                  label("CFU/mL for 50% inhibition of k12 transition (log10 IC50_K12)")        # Garonzik 2016 Table 2 (Log10IC50_K12)
    imax_k12    <- fixed(0.99);              label("Maximum inhibition of k12 at high CFU/mL (Imax_k12, unitless; FIXED)") # Garonzik 2016 Table 2 (IMAX_K12; FIXED)
    lk21        <- fixed(log(50.0));         label("S2 -> S1 transition rate (k21, 1/h; FIXED -- doubling fast)")        # Garonzik 2016 Table 2 (K21; FIXED at 50)
    cfum_log10  <- 9.20;                     label("CFU/mL at which probability of successful replication is 50% (log10 CFU_M)") # Garonzik 2016 Table 2 (Log10 CFU_M)

    # ---- Subpopulation k12 ratios (relative to susceptible) ----
    fr_k12i     <- fixed(1.00);              label("Ratio k12_intermediate / k12_susceptible (FR_K12i, unitless; FIXED)") # Garonzik 2016 Table 2 (FR_K12i; FIXED -- estimated close to 1)
    fr_k12r     <- 0.0442;                   label("Ratio k12_resistant / k12_susceptible (FR_K12r, unitless)")           # Garonzik 2016 Table 2 (FR_K12r)

    # ---- Stimulation of probability of death (per subpopulation) ----
    smax_s      <- fixed(0.99);              label("Maximum stimulation of death probability, susceptible (Smax_s, unitless; FIXED)") # Garonzik 2016 Table 2 footnote: "Smax_s was estimated to be very close to 1 so was fixed to 0.99"
    smax_i      <- 0.515;                    label("Maximum stimulation of death probability, intermediate (Smax_i, unitless)")        # Garonzik 2016 Table 2 (Smax_i)
    smax_r      <- fixed(0);                 label("Maximum stimulation of death probability, resistant (Smax_r, unitless; FIXED)")    # Garonzik 2016 Table 2 footnote: "Smax_r was estimated to be close to zero so was fixed at zero"
    lsc50       <- log(0.0468);              label("Effective daptomycin concentration for 50% Smax (SC50, mg/L)")                     # Garonzik 2016 Table 2 (SC50_s); shared across subpopulations

    # ---- Direct killing (per subpopulation) ----
    lkmax_s     <- log(14.0);                label("Maximum direct killing rate constant, susceptible (Kmax_s, 1/h)")    # Garonzik 2016 Table 2 (Kmax_s)
    lkmax_i     <- log(1.45);                label("Maximum direct killing rate constant, intermediate (Kmax_i, 1/h)")   # Garonzik 2016 Table 2 (Kmax_i)
    kmax_r      <- fixed(0);                 label("Maximum direct killing rate constant, resistant (Kmax_r, 1/h; FIXED)") # Garonzik 2016 Table 2 footnote: "Kmax_r was estimated close to zero and was thus fixed to zero"
    lkc50       <- log(4.81);                label("Effective daptomycin concentration for 50% Kmax (KC50, mg/L)")        # Garonzik 2016 Table 2 (KC50_s); shared across subpopulations

    # ---- Residual error ----
    # The paper reports three residual variability components on log10 CFU/mL:
    # epsilon_CFU = 0.558 (additive on log10 scale, the dominant term),
    # epsilon_Pois = 1.00 (FIXED, Poisson error scaling at low counts), and
    # epsilon_Add = 0.250 (FIXED, additional additive for counts < 5). The
    # nlmixr2 model file carries only the dominant additive epsilon_CFU as
    # addSd; the small-count Poisson and additional-additive components are
    # omitted as a deliberate simplification (see vignette Assumptions and
    # deviations). For simulation use this approximates well above the limit
    # of detection (~log10 = 2), but understates noise as CFU approaches 0.
    addSd <- 0.558; label("Additive residual SD on log10 CFU/mL (epsilon_CFU)")                                              # Garonzik 2016 Table 2 (epsilon_CFU; epsilon_Pois and epsilon_Add omitted)
  })

  model({
    # ---- Per-experiment effective parameters ----
    fact_hs10 <- exp(lfact_hs10)
    fact_hs30 <- exp(lfact_hs30)
    fact_hs50 <- exp(lfact_hs50)
    fact_hs70 <- exp(lfact_hs70)

    cfu0    <- 10 ^ cfu0_log10
    fri     <- 10 ^ fri_log10
    frr     <- 10 ^ frr_log10
    cfu_m   <- 10 ^ cfum_log10
    ic50_k12 <- 10 ^ ic50_k12_log10

    klag    <- 1 / exp(lmtt_lag)
    k12_s   <- 1 / exp(lmtt_k12)
    k12_i   <- fr_k12i * k12_s
    k12_r   <- fr_k12r * k12_s
    k21     <- exp(lk21)

    sc50    <- exp(lsc50)
    kmax_s  <- exp(lkmax_s)
    kmax_i  <- exp(lkmax_i)
    kc50    <- exp(lkc50)

    # ---- Active fraction selection by experimental HS level ----
    # Boolean indicator multipliers (rxode2 evaluates booleans as 0/1). HS = 0
    # gives factive = 1 (no serum, no protein binding); the four estimated
    # factive values cover HS = 10, 30, 50, 70. Any HS not in {0,10,30,50,70}
    # yields factive = 0 (deliberately conspicuous, not interpolated).
    i_hs0  <- HS == 0
    i_hs10 <- HS == 10
    i_hs30 <- HS == 30
    i_hs50 <- HS == 50
    i_hs70 <- HS == 70
    factive <- i_hs0  * 1 +
               i_hs10 * fact_hs10 +
               i_hs30 * fact_hs30 +
               i_hs50 * fact_hs50 +
               i_hs70 * fact_hs70

    # ---- Effective daptomycin concentration (Garonzik 2016 Eq 2) ----
    dap_eff <- factive * dap

    # ---- Total CFU and saturable growth modulators (Eqs 4-7) ----
    cfu_total <- s1 + s2 + i1 + i2 + r1 + r2
    plateau   <- 1 - cfu_total / (cfu_total + cfu_m)        # Eq 5
    rep_term  <- 2 * plateau                                 # Eq 6

    lag_term  <- 1 - exp(-(klag * t) ^ beta_lag)             # Eq 3
    inhib_growth <- 1 - imax_k12 * cfu_total / (cfu_total + ic50_k12)
    k12es_s   <- lag_term * k12_s * inhib_growth             # Eq 7 (susceptible)
    k12es_i   <- lag_term * k12_i * inhib_growth             # Eq 7 (intermediate)
    k12es_r   <- lag_term * k12_r * inhib_growth             # Eq 7 (resistant)

    # ---- Stimulation of probability of death (Eqs 8-9) and direct killing (Eq 10) ----
    # Shared SC50 and KC50 across subpopulations per Garonzik 2016 Table 2
    # (only Smax_s/Smax_i/Smax_r and Kmax_s/Kmax_i/Kmax_r vary across the
    # three subpopulations; SC50 and KC50 are reported as single estimates).
    sti_s   <- smax_s * dap_eff / (dap_eff + sc50)
    sti_i   <- smax_i * dap_eff / (dap_eff + sc50)
    sti_r   <- smax_r * dap_eff / (dap_eff + sc50)
    irep_s  <- 1 - sti_s
    irep_i  <- 1 - sti_i
    irep_r  <- 1 - sti_r
    kill_s  <- kmax_s * dap_eff / (dap_eff + kc50)
    kill_i  <- kmax_i * dap_eff / (dap_eff + kc50)
    kill_r  <- kmax_r * dap_eff / (dap_eff + kc50)

    # ---- ODE system ----
    # Daptomycin solution concentration: static (no chemical degradation in
    # the broth over 24 h per paper Methods).
    d/dt(dap) <- 0

    # Susceptible subpopulation (Garonzik 2016 Eqs 11-12).
    d/dt(s1) <- rep_term * k21 * s2 * irep_s - k12es_s * s1 - kill_s * s1
    d/dt(s2) <- -k21 * s2 + k12es_s * s1 - kill_s * s2

    # Intermediate subpopulation (same structure with subpopulation-specific
    # Smax_i, Kmax_i, k12_i parameters).
    d/dt(i1) <- rep_term * k21 * i2 * irep_i - k12es_i * i1 - kill_i * i1
    d/dt(i2) <- -k21 * i2 + k12es_i * i1 - kill_i * i2

    # Resistant subpopulation (Smax_r = Kmax_r = 0, k12_r much slower).
    d/dt(r1) <- rep_term * k21 * r2 * irep_r - k12es_r * r1 - kill_r * r1
    d/dt(r2) <- -k21 * r2 + k12es_r * r1 - kill_r * r2

    # Initial conditions. All bacteria start in state 1 (vegetative); state 2
    # is empty. Subpopulations are partitioned by FR_I and FR_R relative to
    # the total inoculum CFU0; the susceptible IC absorbs the remainder.
    s1(0) <- cfu0 * (1 - fri - frr)
    i1(0) <- cfu0 * fri
    r1(0) <- cfu0 * frr

    # ---- Observation ----
    # Observation Cc is log10 of total CFU/mL across the six bacterial
    # compartments; a 1e-6 floor protects log10 when all states are driven
    # to zero. The Cc name is the nlmixr2lib single-output convention; the
    # underlying quantity is log10 CFU/mL, not a drug concentration -- see
    # units$concentration. Additive residual SD = 0.558 on log10 CFU/mL
    # (Garonzik 2016 Table 2 epsilon_CFU).
    Cc <- log10(cfu_total + 1e-6)
    Cc ~ add(addSd)
  })
}
