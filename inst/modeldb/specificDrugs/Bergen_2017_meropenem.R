Bergen_2017_meropenem <- function() {
  description <- "In vitro (hollow-fiber infection model). Mechanism-based PK/PD (life-cycle growth) model of meropenem bacterial killing and resistance against Pseudomonas aeruginosa 1280 (meropenem MIC 0.25 mg/L) across simulated critically ill patient renal-function profiles (augmented renal clearance, normal, and impaired). The bacterial population is split into three pre-existing subpopulations of decreasing meropenem susceptibility (susceptible, intermediate, resistant), each described by two states (state 1 preparing for replication, state 2 immediately before replication; six bacterial compartments total). Meropenem acts via inhibition of successful bacterial replication (a Hill-type Inh_Rep function per subpopulation; no direct killing term). The intermediate and resistant subpopulations have higher IC50_Rep and steeper or shallower Hill coefficients than the susceptible subpopulation; the susceptible subpopulation has Imax_Rep and Hill fixed to 1. Meropenem disposition in the HFIM is a fixed-half-life first-order decline parameterised from the upstream popPK model (Mattioli 2016, reference 20 in the source paper); the default half-life is 1.1 h (normal renal function); 0.6 h (augmented renal clearance) and 4.0 h (impaired renal function) are obtained by overriding thalf_mem at simulation time. No patient covariates and no random effects: this is the typical-value MBM fit (Bergen 2017 Table 3) to the simultaneous P. aeruginosa 1280 HFIM data across the three renal-function scenarios and four dosing regimens (2, 1, or 0.5 g q8h plus 1 g q12h for impaired)."
  reference <- "Bergen PJ, Bulitta JB, Kirkpatrick CMJ, Rogers KE, McGregor MJ, Wallis SC, Paterson DL, Nation RL, Lipman J, Roberts JA, Landersdorfer CB. Substantial impact of altered pharmacokinetics in critically ill patients on the antibacterial effects of meropenem evaluated via the dynamic hollow-fiber infection model. Antimicrob Agents Chemother. 2017;61(5):e02642-16. doi:10.1128/AAC.02642-16. Model differential equations (Eqs 1-5) and final parameter estimates (Table 3) are in the main text Materials and Methods + Discussion; HFIM dosing scenarios and concentration summaries are Table 4. Meropenem PK profiles were simulated from the upstream popPK model in reference 20 (Mattioli 2016, AAC; not packaged here)."
  vignette <- "Bergen_2017_meropenem"
  units <- list(
    time          = "hour",
    dosing        = "mg/L (concentration dosed directly into the cmem state, as in the in-vitro HFIM)",
    concentration = "log10 CFU/mL (observation); mg/L (cmem state)"
  )

  # No patient covariates: this is an in vitro mechanism-based model. The
  # meropenem exposure is a state variable (cmem) dosed by the user, not a
  # covariate column. Renal-function scenarios are switched at simulation time
  # by overriding thalf_mem (1.1 h normal, 0.6 h augmented renal clearance,
  # 4.0 h impaired) and choosing dosing amounts to hit the published Cmax/Cmin
  # (Bergen 2017 Table 4).
  covariateData <- list()

  population <- list(
    species          = "in vitro (Pseudomonas aeruginosa 1280, clinical isolate from a critically ill patient with soft-tissue infection; meropenem MIC = 0.25 mg/L)",
    n_subjects       = 1L,
    n_studies        = 1L,
    disease_state    = "Critically ill ICU patient pharmacokinetic scenarios simulated against a susceptible P. aeruginosa isolate",
    model_system     = "Dynamic hollow-fiber in vitro infection model (HFIM) at 36 C, cation-adjusted Mueller-Hinton broth (CAMHB), 10-day duration with periodic sampling for viable counts",
    initial_inoculum = "~10^7.5 CFU/mL (paper Methods); model estimate Log10CFU0 = 6.97 (Table 3) -- the model places ~1 log10 below the observed inoculum because the three subpopulations partition mutation-frequency-seeded cells away from the dominant susceptible pool, so the susceptible CFU_S0 absorbs the remainder of CFU0 and is not directly the observed total inoculum.",
    dose_range       = "Meropenem 2, 1, or 0.5 g q8h as 30-min IV infusions, plus 1 g q12h (impaired renal function only), with concentration-time profiles simulated from the upstream popPK model (Mattioli 2016) for three renal-function scenarios: ARC (CLcr ~285 mL/min, meropenem CL = 34.0 L/h, t_1/2 = 0.6 h), normal (CLcr 120 mL/min, CL = 16.3 L/h, t_1/2 = 1.1 h), and impaired (CLcr ~10 mL/min, CL = 4.1 L/h, t_1/2 = 4.0 h)",
    notes            = paste(
      "Mechanism-based MBM (S-ADAPT, importance sampling pmethod=4) fit simultaneously to all renal functions and dosing regimens.",
      "Single P. aeruginosa isolate (1280); a separate growth-control arm per renal function was included.",
      "Meropenem protein binding = 2% (free fraction = 0.98); fraction excreted unchanged in urine = 79% (Methods).",
      "Resistance was screened on antibiotic-containing agar at 5x and 10x MIC (1.25 and 2.5 mg/L) but the MBM's CFU_I and CFU_R subpopulations are mechanistic (mutation-frequency-seeded) and do NOT directly reflect bacterial counts on those plates (Discussion).",
      "Between-curve variability was set to 15% CV during estimation; no covariate effects were included."
    )
  )

  ini({
    # --- Bacterial life-cycle structural parameters (Table 3) ---
    log10cfu0   <- 6.97;  label("Initial inoculum (log10 CFU/mL)")                          # Bergen 2017 Table 3 (Log CFU0; SE 1.5%)
    log10cfumax <- 9.98;  label("Maximum population size (log10 CFU/mL)")                   # Bergen 2017 Table 3 (Log CFUmax; SE 0.9%)
    lk21        <- fixed(log(50.0)); label("Log replication rate constant k21 (1/h; FIXED)") # Bergen 2017 Table 3 (k21 = 50, fixed: "Bacterial replication was assumed to be fast")

    # Mean generation time per subpopulation (minutes); growth rate k12 = 60 / MGT (1/h).
    # Table 3 calls these rows "k12,X^-1" with units of minutes -- they are mean
    # generation times, not rate constants. The corresponding k12 (1/h) is
    # derived in model().
    mgt_s <- 49.3; label("Mean generation time, susceptible (min)")        # Bergen 2017 Table 3 (k12,S^-1 = 49.3 min; SE 9.5%)
    mgt_i <- 683;  label("Mean generation time, intermediate (min)")       # Bergen 2017 Table 3 (k12,I^-1 = 683 min; SE 19%)
    mgt_r <- 78.5; label("Mean generation time, resistant (min)")          # Bergen 2017 Table 3 (k12,R^-1 = 78.5 min; SE 10.4%)

    # Log10 mutation frequencies seeding the two pre-existing
    # less-susceptible subpopulations. Note these are MECHANISTIC estimates,
    # not the observed agar-plate log10 mutation frequencies in Table 1
    # (Discussion: "the estimated subpopulations did not directly reflect
    # bacterial counts on meropenem-containing agar plates at 5x and 10x MIC").
    log10mf_i <- -3.66; label("Log10 mutation frequency, intermediate subpopulation seed (LogMF_I, log10 unitless)") # Bergen 2017 Table 3 (LogMF_I = -3.66; SE 5.2%)
    log10mf_r <- -6.28; label("Log10 mutation frequency, resistant subpopulation seed (LogMF_R, log10 unitless)")    # Bergen 2017 Table 3 (LogMF_R = -6.28; SE 6.9%)

    # --- Inhibition of successful replication by meropenem (Table 3) ---
    # Imax_Rep_S was estimated at 0.999 and fixed at 1.0 (footnote b);
    # Hill_S was estimated at 1.03 and fixed at 1.0 (footnote c).
    imax_rep_s <- fixed(1.0); label("Maximum inhibition of replication, susceptible (Imax_Rep_S, unitless; FIXED)") # Bergen 2017 Table 3 (footnote b: estimated 0.999, fixed to 1.0)
    imax_rep_i <- 0.673;      label("Maximum inhibition of replication, intermediate (Imax_Rep_I, unitless)")        # Bergen 2017 Table 3 (Imax_Rep_I = 0.673; SE 56.3%)
    imax_rep_r <- 0.956;      label("Maximum inhibition of replication, resistant (Imax_Rep_R, unitless)")           # Bergen 2017 Table 3 (Imax_Rep_R = 0.956; SE 25.0%)

    ic50_rep_s <- 0.648; label("Meropenem concn for 50% Imax_Rep, susceptible (mg/L)")  # Bergen 2017 Table 3 (IC50_Rep_S = 0.648; SE 39.8%)
    ic50_rep_i <- 2.96;  label("Meropenem concn for 50% Imax_Rep, intermediate (mg/L)") # Bergen 2017 Table 3 (IC50_Rep_I = 2.96; SE 45.3%)
    ic50_rep_r <- 6.09;  label("Meropenem concn for 50% Imax_Rep, resistant (mg/L)")    # Bergen 2017 Table 3 (IC50_Rep_R = 6.09; SE 23.3%)

    hill_s <- fixed(1.0); label("Hill coefficient, susceptible (unitless; FIXED)") # Bergen 2017 Table 3 (footnote c: estimated 1.03, fixed to 1.0)
    hill_i <- 7.14;       label("Hill coefficient, intermediate (unitless)")        # Bergen 2017 Table 3 (Hill_I = 7.14; SE 23.9%)
    hill_r <- 2.19;       label("Hill coefficient, resistant (unitless)")           # Bergen 2017 Table 3 (Hill_R = 2.19; SE 62.3%)

    # --- Residual error (Table 3) ---
    addSd <- 0.493; label("Additive residual SD on log10 scale (log10 CFU/mL)") # Bergen 2017 Table 3 (SD_CFU = 0.493)

    # --- Simulated meropenem disposition in the HFIM ---
    # Half-life of the simulated meropenem concentration-time profile in the
    # HFIM. Three scenarios per the paper's Table 4: 0.6 h (augmented renal
    # clearance, CLcr ~285 mL/min), 1.1 h (normal renal function, CLcr 120
    # mL/min), 4.0 h (impaired renal function, CLcr ~10 mL/min). Default is the
    # normal renal function scenario; override at simulation time for ARC or
    # impaired-renal-function scenarios. Fixed because it is not an MBM
    # estimate -- it is an input simulated from the upstream popPK model
    # (Mattioli 2016, paper reference 20) and reported in Table 4.
    thalf_mem <- fixed(1.1); label("Simulated meropenem half-life in HFIM (h; default = normal renal function)") # Bergen 2017 Table 4 (t_1/2 = 1.1 h at CLcr 120 mL/min); not an MBM estimate
  })

  model({
    # 0. Back-transform the log replication rate constant.
    k21 <- exp(lk21)

    # 1. Carrying capacity and total inoculum (log10 -> linear, CFU/mL).
    cfumax <- 10^log10cfumax
    cfu0   <- 10^log10cfu0

    # 2. First-order growth rate constants from mean generation time (Table 3
    #    footer convention: k12,X^-1 reported in minutes -> k12 = 60/MGT in 1/h).
    k12s <- 60 / mgt_s
    k12i <- 60 / mgt_i
    k12r <- 60 / mgt_r

    # 3. Meropenem elimination rate constant from the simulated HFIM half-life.
    kel_mem <- log(2) / thalf_mem    # 1/h

    # 4. Inhibition of successful replication per subpopulation (Hill / Imax;
    #    Bergen 2017 Eq 4). At cmem = 0, inh_rep = 0 and irep = 1 (no inhibition);
    #    at cmem >> IC50, irep approaches 1 - Imax_Rep (= 0 for susceptible).
    #    inh_rep = 0.5 corresponds to net bacterial stasis on this subpopulation
    #    (Methods: "An Inh_Rep_S value of 0.50 resulted in net stasis, whereas an
    #    Inh_Rep_S value of > 0.50 led to bacterial killing due to the
    #    elimination of bacteria that replicate unsuccessfully").
    inh_rep_s <- imax_rep_s * cmem^hill_s / (ic50_rep_s^hill_s + cmem^hill_s)
    inh_rep_i <- imax_rep_i * cmem^hill_i / (ic50_rep_i^hill_i + cmem^hill_i)
    inh_rep_r <- imax_rep_r * cmem^hill_r / (ic50_rep_r^hill_r + cmem^hill_r)
    irep_s <- 1 - inh_rep_s
    irep_i <- 1 - inh_rep_i
    irep_r <- 1 - inh_rep_r

    # 5. Total viable count and replication factor (Bergen 2017 Eq 1 + Eq 3).
    #    REP = 2 * (1 - CFUall / CFUmax) approaches 2 at low density
    #    (100% probability of successful replication) and 1 at saturation
    #    (50% probability). At no-drug steady state, REP settles to 1
    #    (so CFUall stationary = 0.5 * CFUmax = 10^(9.98 - 0.301) = 10^9.68).
    cfu_all <- bact_susceptible1 + bact_susceptible2 +
               bact_intermediate1 + bact_intermediate2 +
               bact_resistant1 + bact_resistant2
    rep_factor <- 2 * (1 - cfu_all / cfumax)

    # 6. Life-cycle growth model: two states per subpopulation (Bergen 2017
    #    Eq 2 + Eq 5). State 1 = preparing for replication; state 2 = immediately
    #    before replication. State-2 cells leave at k21 (attempt to replicate);
    #    the fraction (REP * irep) / 2 of attempts are successful and produce
    #    two daughter cells back into state 1; the rest die. Drug acts only on
    #    the replication arm (no direct killing term in this paper).
    d/dt(bact_susceptible1)  <- rep_factor * irep_s * k21 * bact_susceptible2  - k12s * bact_susceptible1
    d/dt(bact_susceptible2)  <- -k21 * bact_susceptible2  + k12s * bact_susceptible1
    d/dt(bact_intermediate1) <- rep_factor * irep_i * k21 * bact_intermediate2 - k12i * bact_intermediate1
    d/dt(bact_intermediate2) <- -k21 * bact_intermediate2 + k12i * bact_intermediate1
    d/dt(bact_resistant1)    <- rep_factor * irep_r * k21 * bact_resistant2    - k12r * bact_resistant1
    d/dt(bact_resistant2)    <- -k21 * bact_resistant2    + k12r * bact_resistant1

    # 7. Meropenem concentration in the HFIM (mg/L): dosed by the user as a
    #    concentration, declines with the simulated half-life. To replicate a
    #    clinical regimen (e.g., 2 g q8h as a 0.5-h infusion at normal renal
    #    function), the vignette converts the target Cmax to an infusion rate
    #    using the steady-state IV-infusion formula and doses cmem at that
    #    rate; see vignettes/articles/Bergen_2017_meropenem.Rmd.
    d/dt(cmem) <- -kel_mem * cmem

    # 8. Initial conditions (Methods + Table 3). The total inoculum CFU0 is
    #    partitioned by mutation frequency: the intermediate and resistant
    #    subpopulations are seeded at 10^logMF_I and 10^logMF_R of CFU0, and
    #    the susceptible subpopulation absorbs the remainder. All bacteria
    #    start in state 1; state-2 initial conditions are 0 (balanced-growth
    #    state-2 fraction is ~k12/k21 = ~2-3% for susceptible, < 1% for
    #    intermediate and resistant -- a brief start-up transient).
    bact_susceptible1(0)  <- cfu0 * (1 - 10^log10mf_i - 10^log10mf_r)
    bact_intermediate1(0) <- cfu0 * 10^log10mf_i
    bact_resistant1(0)    <- cfu0 * 10^log10mf_r

    # 9. Derived output: less-susceptible viable count (intermediate plus
    #    resistant subpopulations across both life-cycle states). These do
    #    not directly correspond to colony counts on antibiotic-containing
    #    agar plates (5x/10x MIC) because the MBM's subpopulations are
    #    mechanism-fit, not plate-defined (Discussion).
    cfu_less_susc <- bact_intermediate1 + bact_intermediate2 +
                     bact_resistant1 + bact_resistant2

    # 10. Observation: total viable count on the log10 scale; an additive 1e-6
    #     floor protects log10 when all states are driven near zero by sustained
    #     bactericidal exposure. Cc is the nlmixr2lib single-output convention;
    #     the underlying quantity is log10 CFU/mL, not a drug concentration --
    #     see units$concentration. Additive residual SD on log10 scale matches
    #     Table 3 SD_CFU.
    Cc <- log10(cfu_all + 1e-6)
    Cc ~ add(addSd)
  })
}
