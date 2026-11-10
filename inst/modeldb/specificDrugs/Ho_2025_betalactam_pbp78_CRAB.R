Ho_2025_betalactam_pbp78_CRAB <- function() {
  description <- "In vitro (Acinetobacter baumannii HUMC1 wild type and the isogenic HUMC1 delta-PBP7/8 (pbpG) knockout). Mechanism-based pharmacodynamic model of static-concentration time-kill against carbapenem-resistant A. baumannii by five antibiotics -- meropenem, imipenem, aztreonam, ceftazidime/avibactam and polymyxin B -- as monotherapy and in beta-lactam-plus-polymyxin-B combinations. Two co-existing bacterial subpopulations (a majority susceptible population and a resistant population seeded at the estimated mutant frequency) each follow the Bulitta life-cycle growth structure: a growing state 1 that transits to a replicative state 2 at the mean-generation-time rate and divides back to two state-1 cells at the fast (not rate-limiting) rate KDiv, gated toward the carrying capacity by a logistic replication factor. Each drug kills through its own Hill function with separate maximal killing rate constants for the susceptible and resistant subpopulations and a shared KC50. PBP7/8 knockout enters as a covariate on KC50 (the strain indicator BACT: 0 = wild type, 1 = delta-PBP7/8 knockout), lowering the effective KC50 of every drug through outer-membrane permeabilisation. This is the 'pharmacodynamic evaluation of adjuvant targets' (PEAT) model used to quantify PBP7/8 as a beta-lactam adjuvant target."
  reference <- paste(
    "Ho BM, Jin J, Sanborn J, Nguyen TD, Singh N, Cheng C, Nasief NN,",
    "Carlino-MacDonald U, Russo TA, Smith NM. Pharmacodynamic Evaluation of",
    "Adjuvant Targets: Low Molecular Weight PBP7/8 Effects on beta-Lactam",
    "Activity Against Carbapenem-Resistant Acinetobacter Baumannii.",
    "Pharmaceuticals (Basel). 2025;18(6):918. doi:10.3390/ph18060918. PMID:",
    "40573313. PMCID: PMC12196172. Final-model parameter estimates: Table 1.",
    "Structural model equations: supplementary Figure S5 ('Model Equations',",
    "equations 1-21). Life-cycle growth / subpopulation framework and the",
    "fast division rate KDiv (fixed at 50 /h, not rate-limiting): Bulitta JB,",
    "Ly NS, Yang JC, Forrest A, Jusko WJ, Tsuji BT, Antimicrob Agents",
    "Chemother 2009;53:46-56 (doi:10.1128/AAC.00489-08), reference 29 of the",
    "source paper and extracted in this package as Bulitta_2009_ceftazidime.R.",
    sep = " "
  )
  vignette <- "Ho_2025_betalactam_pbp78_CRAB"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # Static-concentration time-kill (SCTK): the "dose" is the fixed bath
  # concentration placed directly into each drug state at time 0, so the dosing
  # targets are neither `depot` nor `central`.
  dosing <- c("cmem", "cipm", "catm", "ccef", "cpmb")

  # `cmem` is the canonical meropenem bath-concentration state. `cipm` follows
  # the established imipenem usage (Yadav_2017_imipenem_tobramycin_*), and
  # `cpmb` the established polymyxin B usage (Mahadevan_2026_polymyxinB_*);
  # `catm` (aztreonam) and `ccef` (ceftazidime/avibactam) extend the same
  # `c<drug>` bath-state family and are declared here. The bacterial states use
  # the canonical single-axis `bact_<phenotype><lifecycle-digit>` scheme and
  # need no declaration.
  paper_specific_compartments <- c("cipm", "catm", "ccef", "cpmb")

  compartmentData <- list(
    cmem = list(analyte = "meropenem", units = "mg/L", specimen = "administration site", verified = TRUE),
    cipm = list(analyte = "imipenem", units = "mg/L", specimen = "administration site", verified = TRUE),
    catm = list(analyte = "aztreonam", units = "mg/L", specimen = "administration site", verified = TRUE),
    ccef = list(
      analyte = "ceftazidime/avibactam (concentration of the ceftazidime component)",
      units = "mg/L",
      specimen = "administration site",
      verified = TRUE
    ),
    cpmb = list(analyte = "polymyxin B", units = "mg/L", specimen = "administration site", verified = TRUE),
    bact_susceptible1 = list(
      analyte = "Acinetobacter baumannii HUMC1, susceptible subpopulation, growing life-cycle state 1",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_susceptible2 = list(
      analyte = "Acinetobacter baumannii HUMC1, susceptible subpopulation, replicative life-cycle state 2",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_resistant1 = list(
      analyte = "Acinetobacter baumannii HUMC1, resistant subpopulation, growing life-cycle state 1",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    bact_resistant2 = list(
      analyte = "Acinetobacter baumannii HUMC1, resistant subpopulation, replicative life-cycle state 2",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    BACT = list(
      description = "Bacterial-strain indicator for the isogenic HUMC1 pair, entering as a covariate on the KC50 of every drug (the PBP7/8-knockout effect). 0 = HUMC1 wild type, 1 = HUMC1 delta-PBP7/8 (pbpG) knockout.",
      units = "(categorical)",
      type = "categorical",
      source_name = "PBP7/8 knockout status",
      reference_category = 0,
      notes = paste(
        "0 = HUMC1 wild type (MIC meropenem 32 mg/L, MIC polymyxin B 2 mg/L); 1 = HUMC1 delta-PBP7/8, the isogenic low-molecular-weight PBP7/8 (pbpG) knockout.",
        "The knockout enters each drug's effective KC50 multiplicatively as kc50_eff = kc50 * exp(e_bact_kc50_<drug> * BACT), so BACT must be coded 0/1 with 0 as the wild-type reference; the beta is negative for every drug, i.e. knockout lowers the KC50 (greater potency).",
        "The paper modelled both isolates simultaneously and tested the PBP7/8 knockout as a covariate effect on KC50 (Methods, 'Mechanism-Based Mathematical Model'; Results, 'Mechanism-Based Model of PBP7/8 Effects').",
        "The off-target effects of the knockout were excluded with the complemented derivative HUMC1 delta-pbpG/pNLAC1[pbpG] (reference 16 of the source paper, Russo 2023, doi:10.1128/aac.01033-22)."
      )
    )
  )

  population <- list(
    species = "in vitro (Acinetobacter baumannii HUMC1 wild type and isogenic delta-PBP7/8 knockout)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "Acinetobacter baumannii HUMC1 (carbapenem-resistant clinical isolate) and its isogenic HUMC1 delta-PBP7/8 (pbpG) knockout",
    system = "24 h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton broth (25 mg/L Ca2+, 12.5 mg/L Mg2+), 37 C under constant shaking; viable counts at 0, 1, 2, 4, 6, 8 and 24 h",
    temperature = "37 C",
    duration = "24 h",
    starting_inoculum = "~10^7 CFU/mL target; the model estimated 10^6.71 CFU/mL",
    limit_of_detection = "10^2 CFU/mL",
    mic_values = c(
      `HUMC1 wild type, polymyxin B` = "2 mg/L",
      `HUMC1 wild type, meropenem` = "32 mg/L"
    ),
    dose_range = paste(
      "static bath concentrations: meropenem 25 / 50 / 100 mg/L, imipenem 7.7 / 15.3 / 30.6 mg/L,",
      "aztreonam 33 / 66 / 132 mg/L, ceftazidime/avibactam 37-5 / 73-9.9 / 146-19.8 mg/L,",
      "polymyxin B 0.75 / 1.5 / 3 / 6 mg/L; monotherapy and each beta-lactam plus polymyxin B.",
      "Concentrations chosen as 0.5x / 1x / 2x the predicted free Cmax (beta-lactams) or free steady-state average (polymyxin B) from Monte Carlo simulation of published population-PK studies."
    ),
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "16 monotherapy arms and 48 beta-lactam-plus-polymyxin-B combination arms per isolate (Results, 'Static Concentration Time-Kill Experiments').",
      "Mechanism-based model fit in Monolix 2023R1 to the log10-transformed total-population counts of both isolates simultaneously; counts below the 10^2 CFU/mL limit were handled as censored (M3) data.",
      "Interexperimental variability was fixed at a 5% coefficient of variation on every model parameter (assumed, not estimated). Residual variability was additive on the log10 scale, which corresponds to proportional error on the natural scale; its magnitude is not reported in Table 1.",
      "The average KC50 decline from PBP7/8 knockout was 59.7% across the beta-lactams and 39.0% for polymyxin B.",
      "A prospective 7-day hollow-fiber infection model (HFIM) validated the model under dynamic dosing (Figures 2 and 3)."
    )
  )

  ini({
    # ---- Bacterial growth and subpopulation structure (Table 1) ------------
    log10cfumax <- 8.20; label("Maximum carrying-capacity population size (log10 CFU/mL)")                 # Table 1: LOGCFUmax = 8.20 (0.91%)
    log10inoc <- 6.71; label("Starting inoculum (log10 CFU/mL)")                                          # Table 1: LOGINOC = 6.71 (0.38%)
    log10mutr <- 5.03; label("Negative log10 population mutation frequency (unitless); MUTR = 10^-log10mutr") # Table 1: LOGMUTR = 5.03 (0.21%); eq S5-5: MUTR = 10^(-1*LOGMUTR)
    mgts <- 75.3; label("Mean generation time of the susceptible (S) subpopulation (min)")                # Table 1: MGTS = 75.3 (0.51%)
    mgtr <- 46.4; label("Mean generation time of the resistant (R) subpopulation (min)")                  # Table 1: MGTR = 46.4 (1.06%)

    # KDiv, the fast state-2 -> state-1 division rate, is used in equations
    # S5-14/16/18/20 but is not reported in Table 1. It is a fixed structural
    # constant of the Bulitta 2009 life-cycle framework (reference 29), where
    # the fast doubling step is fixed to 50 /h as not rate-limiting; the same
    # convention is used by Bulitta_2009_ceftazidime.R and the Mahadevan_2026
    # siblings in this package. Dynamics are governed by the slow generation
    # step (60/MGT), so the model is insensitive to the exact KDiv value.
    lkdiv <- fixed(log(50)); label("Log fast division (state 2 -> state 1) rate constant KDiv (1/h; not rate-limiting)") # Bulitta 2009 (ref 29), not reported in Ho 2025

    # ---- Maximal killing rate constants, susceptible subpopulation (Table 1) ----
    kmax_mer_s <- 2.09; label("Maximal killing rate of the S subpopulation by meropenem (1/h)")           # Table 1: Kmax,S,MER = 2.09 (2.96%)
    kmax_imi_s <- 1.56; label("Maximal killing rate of the S subpopulation by imipenem (1/h)")            # Table 1: Kmax,S,IMI = 1.56 (2.54%)
    kmax_atm_s <- 1.12; label("Maximal killing rate of the S subpopulation by aztreonam (1/h)")           # Table 1: Kmax,S,ATM = 1.12 (0.98%)
    kmax_cef_s <- 1.52; label("Maximal killing rate of the S subpopulation by ceftazidime/avibactam (1/h)") # Table 1: Kmax,S,CEF = 1.52 (3.12%)
    kmax_pmb_s <- 2.98; label("Maximal killing rate of the S subpopulation by polymyxin B (1/h)")         # Table 1: Kmax,S,PMB = 2.98 (1.11%)

    # ---- Maximal killing rate constants, resistant subpopulation (Table 1) ----
    kmax_mer_r <- 1.13; label("Maximal killing rate of the R subpopulation by meropenem (1/h)")           # Table 1: Kmax,R,MER = 1.13 (2.01%)
    kmax_imi_r <- 0.83; label("Maximal killing rate of the R subpopulation by imipenem (1/h)")            # Table 1: Kmax,R,IMI = 0.83 (1.90%)
    kmax_atm_r <- 0.76; label("Maximal killing rate of the R subpopulation by aztreonam (1/h)")           # Table 1: Kmax,R,ATM = 0.76 (1.60%)
    kmax_cef_r <- 0.89; label("Maximal killing rate of the R subpopulation by ceftazidime/avibactam (1/h)") # Table 1: Kmax,R,CEF = 0.89 (1.97%)
    kmax_pmb_r <- 1.12; label("Maximal killing rate of the R subpopulation by polymyxin B (1/h)")         # Table 1: Kmax,R,PMB = 1.12 (0.84%)

    # ---- KC50: concentration for 50% of maximal killing (Table 1) ----------
    # Shared between the S and R subpopulations (equations S5-10 and S5-11 both
    # use KC50,Drug); only Kmax differs by subpopulation.
    kc50_mer <- 47.6; label("Meropenem concentration for 50% of maximal killing (mg/L)")                  # Table 1: KC50,MER = 47.6 (9.42%)
    kc50_imi <- 16.8; label("Imipenem concentration for 50% of maximal killing (mg/L)")                   # Table 1: KC50,IMI = 16.8 (10.6%)
    kc50_atm <- 78.6; label("Aztreonam concentration for 50% of maximal killing (mg/L)")                  # Table 1: KC50,ATM = 78.6 (14.4%)
    kc50_cef <- 76.0; label("Ceftazidime/avibactam concentration for 50% of maximal killing (mg/L)")      # Table 1: KC50,CEF = 76.0 (10.8%)
    kc50_pmb <- 1.00; label("Polymyxin B concentration for 50% of maximal killing (mg/L)")                # Table 1: KC50,PMB = 1.00 (0.34%)

    # ---- Effect of PBP7/8 knockout (BACT) on each drug KC50 (Table 1) -------
    # kc50_eff = kc50 * exp(beta * BACT); all beta negative (knockout lowers KC50).
    e_bact_kc50_mer <- -0.80; label("Log fractional change in meropenem KC50 from PBP7/8 knockout (unitless)")           # Table 1: betaKC50,MER = -0.80 (13.0%); calc. KO KC50 = 21.4
    e_bact_kc50_imi <- -1.10; label("Log fractional change in imipenem KC50 from PBP7/8 knockout (unitless)")            # Table 1: betaKC50,IMI = -1.10 (10.4%); calc. KO KC50 listed as 6.6 (see vignette errata)
    e_bact_kc50_atm <- -1.02; label("Log fractional change in aztreonam KC50 from PBP7/8 knockout (unitless)")           # Table 1: betaKC50,ATM = -1.02 (19.1%); calc. KO KC50 = 28.3
    e_bact_kc50_cef <- -0.89; label("Log fractional change in ceftazidime/avibactam KC50 from PBP7/8 knockout (unitless)") # Table 1: betaKC50,CEF = -0.89 (17.2%); calc. KO KC50 = 31.2
    e_bact_kc50_pmb <- -0.50; label("Log fractional change in polymyxin B KC50 from PBP7/8 knockout (unitless)")         # Table 1: betaKC50,PMB = -0.50 (5.25%); calc. KO KC50 = 0.61

    # ---- Hill shape coefficients (Table 1) ---------------------------------
    hill_mer <- 0.91; label("Hill coefficient of meropenem killing (unitless)")                           # Table 1: HillMER = 0.91 (1.46%)
    hill_imi <- 0.88; label("Hill coefficient of imipenem killing (unitless)")                            # Table 1: HillIMI = 0.88 (1.19%)
    hill_atm <- 0.51; label("Hill coefficient of aztreonam killing (unitless)")                           # Table 1: HillATM = 0.51 (4.08%)
    hill_cef <- 0.90; label("Hill coefficient of ceftazidime/avibactam killing (unitless)")               # Table 1: HillCER [sic, ceftazidime/avibactam] = 0.90 (1.99%)
    hill_pmb <- 0.93; label("Hill coefficient of polymyxin B killing (unitless)")                         # Table 1: HillPMB = 0.93 (0.35%)

    # ---- Residual variability ----------------------------------------------
    # Additive on the log10 scale (= proportional on the natural scale); the
    # paper states the error model but does not report its magnitude in Table 1.
    addSd <- fixed(0); label("Additive residual SD on the log10 bacterial count (log10 CFU/mL); magnitude not reported") # Methods: additive residual on log10 data; value not in Table 1
  })

  model({
    # ---- 1. Back-transforms and derived rate constants (equations S5-1,2,5) ----
    kdiv <- exp(lkdiv)
    cfumax <- 10^log10cfumax
    inoc <- 10^log10inoc
    mutr <- 10^(-log10mutr)
    kgr_s <- 60 / mgts
    kgr_r <- 60 / mgtr

    # ---- 2. Initial subpopulation sizes (equations S5-6,7) -----------------
    ic_r <- inoc * mutr
    ic_s <- inoc - ic_r

    # ---- 3. Effective KC50 under PBP7/8 knockout (covariate on KC50) --------
    kc50_mer_eff <- kc50_mer * exp(e_bact_kc50_mer * BACT)
    kc50_imi_eff <- kc50_imi * exp(e_bact_kc50_imi * BACT)
    kc50_atm_eff <- kc50_atm * exp(e_bact_kc50_atm * BACT)
    kc50_cef_eff <- kc50_cef * exp(e_bact_kc50_cef * BACT)
    kc50_pmb_eff <- kc50_pmb * exp(e_bact_kc50_pmb * BACT)

    # ---- 4. Per-drug Hill killing of each subpopulation (equations S5-10,11) ----
    kill_s <- kmax_mer_s * cmem^hill_mer / (kc50_mer_eff^hill_mer + cmem^hill_mer) +
      kmax_imi_s * cipm^hill_imi / (kc50_imi_eff^hill_imi + cipm^hill_imi) +
      kmax_atm_s * catm^hill_atm / (kc50_atm_eff^hill_atm + catm^hill_atm) +
      kmax_cef_s * ccef^hill_cef / (kc50_cef_eff^hill_cef + ccef^hill_cef) +
      kmax_pmb_s * cpmb^hill_pmb / (kc50_pmb_eff^hill_pmb + cpmb^hill_pmb)
    kill_r <- kmax_mer_r * cmem^hill_mer / (kc50_mer_eff^hill_mer + cmem^hill_mer) +
      kmax_imi_r * cipm^hill_imi / (kc50_imi_eff^hill_imi + cipm^hill_imi) +
      kmax_atm_r * catm^hill_atm / (kc50_atm_eff^hill_atm + catm^hill_atm) +
      kmax_cef_r * ccef^hill_cef / (kc50_cef_eff^hill_cef + ccef^hill_cef) +
      kmax_pmb_r * cpmb^hill_pmb / (kc50_pmb_eff^hill_pmb + cpmb^hill_pmb)

    # ---- 5. Total viable count and the replication plateau (equations S5-12,13) ----
    CFUtot <- bact_susceptible1 + bact_susceptible2 + bact_resistant1 + bact_resistant2
    repf <- 2 * (1 - CFUtot / (cfumax + CFUtot))

    # ---- 6. Bacterial life-cycle system (equations S5-14,16,18,20) ----------
    d/dt(bact_susceptible1) <-
      repf * (kdiv * bact_susceptible2) - kgr_s * bact_susceptible1 - kill_s * bact_susceptible1
    d/dt(bact_susceptible2) <-
      kgr_s * bact_susceptible1 - kdiv * bact_susceptible2 - kill_s * bact_susceptible2
    d/dt(bact_resistant1) <-
      repf * (kdiv * bact_resistant2) - kgr_r * bact_resistant1 - kill_r * bact_resistant1
    d/dt(bact_resistant2) <-
      kgr_r * bact_resistant1 - kdiv * bact_resistant2 - kill_r * bact_resistant2

    # ---- 7. Static drug exposures ------------------------------------------
    d/dt(cmem) <- 0
    d/dt(cipm) <- 0
    d/dt(catm) <- 0
    d/dt(ccef) <- 0
    d/dt(cpmb) <- 0

    # ---- 8. Initial conditions (equations S5-15,17,19,21) -------------------
    bact_susceptible1(0) <- ic_s
    bact_resistant1(0) <- ic_r

    # ---- 9. Observation -----------------------------------------------------
    # The 1e-6 CFU/mL regularisation keeps log_cfu finite when a regimen drives
    # the population toward numerical extinction (the solver can undershoot to a
    # small negative value, and log10 of that is NaN). It is far below the
    # 10^2 CFU/mL limit of detection, so it never perturbs an observable count.
    # Same device as Mahadevan_2026_polymyxinB_meropenem_fosfomycin_*.R.
    log_cfu <- log10(CFUtot + 1e-6)
    log_cfu ~ add(addSd)
  })
}
