Landersdorfer_2018_imipenem_tobramycin <- function() {
  description <- "In vitro (carbapenem-resistant Acinetobacter baumannii FADDI-AB034). Mechanism-based pharmacodynamic model for an imipenem-plus-tobramycin combination in a 7-day hollow-fiber infection model (HFIM). Three pre-existing bacterial subpopulations differing in imipenem and tobramycin susceptibility (population 1 IPM-S/TOB-S, population 2 IPM-R/TOB-I, population 3 IPM-I/TOB-R) each follow a Bulitta two-state life-cycle growth model (S1 -> S2 -> 2*S1 with replication rate k21 fixed and per-subpopulation mean generation time 1/k12). Logistic carrying-capacity attenuation applies to the replicating step. Imipenem kills with a sigmoidal Hill function (Hill = 3 fixed) and tobramycin with an Emax function; mechanistic synergy is encoded as a discrete 70-fold reduction of the imipenem KC50 against population 3 when the tobramycin concentration meets or exceeds 1.15 mg/L (i.e. tobramycin permeabilizing the outer bacterial membrane toward imipenem). Imipenem and tobramycin concentrations are external time-varying inputs (covariates CONC_IPM_MGL and CONC_TOB_MGL); the model contains no human PK component."
  reference <- paste(
    "Landersdorfer CB, Yadav R, Rogers KE, Kim TH, Shin BS, Boyce JD, Nation RL, Bulitta JB. (2018).",
    "Combating carbapenem-resistant Acinetobacter baumannii by an optimized imipenem-plus-tobramycin",
    "dosage regimen: prospective validation via hollow-fiber infection and mathematical modeling.",
    "Antimicrobial Agents and Chemotherapy 62(4):e02053-17.",
    "doi:10.1128/AAC.02053-17.",
    sep = " "
  )
  vignette <- "Landersdorfer_2018_imipenem_tobramycin"
  units <- list(time = "hour", dosing = "mg/L (drug input concentration)", concentration = "log10 CFU/mL (observation); mg/L (drug covariates)")

  covariateData <- list(
    CONC_IPM_MGL = list(
      description        = "Unbound imipenem concentration in the hollow-fiber growth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally; the HFIM experiment used continuous infusion targeting the 5th-percentile (7.6 mg/L), median (13.4 mg/L), and 95th-percentile (23.3 mg/L) unbound concentrations expected from imipenem 4 g/day continuous infusion in critically ill patients. In-vitro applied drug-concentration covariate, registered as the canonical CONC_IPM_MGL in inst/references/covariate-columns.md.",
      source_name        = "Imipenem concentration (paper Methods and Fig. 1 legend)"
    ),
    CONC_TOB_MGL = list(
      description        = "Unbound tobramycin concentration in the hollow-fiber growth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally; the HFIM experiment simulated the two-compartment unbound tobramycin profile produced by 7 mg/kg q24h 0.5-h infusions (observed peak 12.3 mg/L at 1.2 h, trough 1.37 mg/L at 23 h; pump-flow rate switched at 5 h each day to mimic alpha/beta phases). In-vitro applied drug-concentration covariate, registered as the canonical CONC_TOB_MGL in inst/references/covariate-columns.md.",
      source_name        = "Tobramycin concentration (paper Methods, Fig. S1 reference)"
    )
  )

  population <- list(
    species             = "in vitro (Acinetobacter baumannii, carbapenem-resistant clinical isolate FADDI-AB034)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Acinetobacter baumannii FADDI-AB034 (clinical CRAB isolate; imipenem MIC 32 mg/L, tobramycin MIC 2 mg/L)",
    system              = "Hollow-fiber dynamic in vitro infection model (HFIM); 7-day exposure; total and resistant-subpopulation viable counts on antibiotic-free and antibiotic-containing agar (imipenem 1.75x and 3x MIC, tobramycin 3x and 5x MIC)",
    medium              = "Standard cation-adjusted Mueller-Hinton broth circulating through the hollow-fiber cartridge",
    temperature         = "37 C",
    duration            = "168 h (7 days); dense sampling within the first 24 h and at 47, 71, 95, 119, 143, and 168 h",
    inoculum            = "~10^7.2 CFU/mL (mirrors bacterial densities in severe clinical infections)",
    mic_values          = c(imipenem = "32 mg/L", tobramycin = "2 mg/L"),
    regimens            = "Imipenem monotherapy at constant 7.6, 13.4, and 23.3 mg/L unbound (corresponding to 5th, 50th, and 95th percentile steady-state concentrations from 4 g/day continuous infusion with a 1 g loading dose); tobramycin monotherapy at 7 mg/kg q24h (0.5-h infusions); each imipenem level combined with the same tobramycin regimen; antibiotic-free growth controls",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. The MBM was fit in S-ADAPT to viable-count data; the coefficient of correlation for observed vs individual (population) fitted log10 viable counts was 0.995 (0.968). Random effects (eta) are NOT estimated in the source: the paper reports population mean parameter estimates with relative standard errors only, on a single bacterial isolate; consequently the packaged model contains no etas and is intended for typical-value simulation only. See Landersdorfer 2018 Methods (page 632 = e02053-17 p. 2) and Table 1."
  )

  ini({
    # =============================================================
    # Bacterial life-cycle (Bulitta two-state LCGM)
    # =============================================================
    # k21 = replication rate constant (S2 -> 2*S1 doubling step);
    # 1/k12 = mean generation time for the slow S1 -> S2 transition.
    # The replication rate constant is shared across subpopulations
    # and fixed at 50 /h per Table 1 footnote (a sensitivity-tested
    # value, not an estimate).
    lk21 <- fixed(log(50))
    label("Replication rate constant k21 (1/h; S2 -> 2*S1 doubling; FIXED)")  # Landersdorfer 2018 Table 1

    # Mean generation time per subpopulation (minutes). Inverse
    # gives the slow S1 -> S2 rate constant k12_p (1/min initially
    # -- converted to 1/h inside model()). Table 1 footnote (c)
    # states the IPM-I/TOB-R generation time was held equal to the
    # IPM-S/TOB-S value; the file follows the same constraint by
    # using mgt_ss for both populations 1 and 3 in model().
    mgt_ss <- 33.7
    label("Mean generation time IPM-S/TOB-S subpop (population 1; min)")  # Landersdorfer 2018 Table 1 (also shared by population 3 per footnote c)
    mgt_ri <- 63.2
    label("Mean generation time IPM-R/TOB-I subpop (population 2; min)")  # Landersdorfer 2018 Table 1

    # Initial total inoculum and carrying capacity (log10 CFU/mL).
    logcfu0   <- 7.29
    label("Initial total inoculum (log10 CFU/mL)")  # Landersdorfer 2018 Table 1 (log CFU0)
    logcfumax <- 10.1
    label("Maximum population size (log10 CFU/mL)")  # Landersdorfer 2018 Table 1 (logCFUmax)

    # =============================================================
    # Initial subpopulation fractions (log10 mutation frequencies)
    # =============================================================
    # The IPM-R/TOB-I subpopulation (population 2) has initial
    # fraction 10^log10_mut_ipm of the total; the IPM-I/TOB-R
    # subpopulation (population 3) has 10^log10_mut_tob; the
    # remainder is IPM-S/TOB-S (population 1).
    log10_mut_ipm <- -5.59
    label("log10 mutation frequency at t = 0 for IPM-R/TOB-I (population 2)")  # Landersdorfer 2018 Table 1 (log MUT,IPM)
    log10_mut_tob <- -3.22
    label("log10 mutation frequency at t = 0 for IPM-I/TOB-R (population 3)")  # Landersdorfer 2018 Table 1 (log MUT,TOB)

    # =============================================================
    # Imipenem killing (Hill model)
    # =============================================================
    lkmax_ipm <- log(1.74)
    label("Maximum imipenem killing rate Kmax,IPM (1/h)")  # Landersdorfer 2018 Table 1
    lhill_ipm <- fixed(log(3.0))
    label("Imipenem Hill coefficient (unitless; FIXED per sensitivity analysis)")  # Landersdorfer 2018 Table 1 footnote (d)

    # Imipenem KC50 per subpopulation. Population 3 has two
    # values: the monotherapy KC50 (without sufficient tobramycin)
    # and the combination KC50 (when CONC_TOB_MGL >= tob_cut). model() will
    # switch between them at the threshold.
    lkc50_ss_ipm        <- log(0.175)
    label("Imipenem KC50 against IPM-S/TOB-S (population 1; mg/L)")  # Landersdorfer 2018 Table 1
    lkc50_ri_ipm        <- log(645)
    label("Imipenem KC50 against IPM-R/TOB-I (population 2; mg/L)")  # Landersdorfer 2018 Table 1
    lkc50_ir_ipm_mono   <- log(112)
    label("Imipenem KC50 against IPM-I/TOB-R (population 3; monotherapy or CONC_TOB_MGL < tob_cut; mg/L)")  # Landersdorfer 2018 Table 1
    lkc50_ir_ipm_combo  <- log(1.60)
    label("Imipenem KC50 against IPM-I/TOB-R (population 3; with CONC_TOB_MGL >= tob_cut; mg/L)")  # Landersdorfer 2018 Table 1

    # Tobramycin concentration triggering mechanistic synergy
    # (outer-membrane permeabilization shifts KC50,IR,IPM 70-fold).
    tob_cut <- 1.15
    label("Tobramycin threshold for mechanistic synergy (mg/L)")  # Landersdorfer 2018 Table 1 (TOB_cut)

    # =============================================================
    # Tobramycin killing (Emax model; Hill = 1 implicit)
    # =============================================================
    lkmax_tob_ss <- log(4.69)
    label("Maximum tobramycin killing rate against IPM-S/TOB-S (population 1; 1/h; also applied to population 3 per Table 1 footnote e)")  # Landersdorfer 2018 Table 1
    lkmax_tob_ri <- log(0.992)
    label("Maximum tobramycin killing rate against IPM-R/TOB-I (population 2; 1/h)")  # Landersdorfer 2018 Table 1

    lkc50_ss_tob <- log(0.156)
    label("Tobramycin KC50 against IPM-S/TOB-S (population 1; mg/L)")  # Landersdorfer 2018 Table 1
    lkc50_ri_tob <- log(0.316)
    label("Tobramycin KC50 against IPM-R/TOB-I (population 2; mg/L)")  # Landersdorfer 2018 Table 1
    lkc50_ir_tob <- log(27.7)
    label("Tobramycin KC50 against IPM-I/TOB-R (population 3; mg/L)")  # Landersdorfer 2018 Table 1

    # =============================================================
    # Residual error
    # =============================================================
    # Table 1: "SD of residual error on log10 scale" = 0.304. The
    # packaged observation Cc is log10(total viable count) (per
    # the Wicha 2017 convention for HFIM PD models, where Cc
    # carries log10 CFU/mL rather than a drug concentration); the
    # additive residual SD enters in log10 units directly.
    addSd <- 0.304
    label("Additive residual SD on log10 total viable count (log10 CFU/mL)")  # Landersdorfer 2018 Table 1 (SD_CFU)
  })

  model({
    # Effective parameters (exponentiation of log-transformed inputs).
    # No etas are estimated in the source MBM; see population$notes.
    k21       <- exp(lk21)
    kmax_ipm  <- exp(lkmax_ipm)
    hill_ipm  <- exp(lhill_ipm)
    kc50_ss_ipm       <- exp(lkc50_ss_ipm)
    kc50_ri_ipm       <- exp(lkc50_ri_ipm)
    kc50_ir_ipm_mono  <- exp(lkc50_ir_ipm_mono)
    kc50_ir_ipm_combo <- exp(lkc50_ir_ipm_combo)
    kmax_tob_ss <- exp(lkmax_tob_ss)
    kmax_tob_ri <- exp(lkmax_tob_ri)
    kc50_ss_tob <- exp(lkc50_ss_tob)
    kc50_ri_tob <- exp(lkc50_ri_tob)
    kc50_ir_tob <- exp(lkc50_ir_tob)

    # Per-subpopulation k12 (slow S1 -> S2 transition, 1/h). Mean
    # generation times are in minutes per Table 1, so divide
    # 60 min/h by the minutes-per-generation to get 1/h. Per Table
    # 1 footnote (c) populations 1 and 3 share mgt_ss.
    k12_ss_base <- 60 / mgt_ss
    k12_ri_base <- 60 / mgt_ri
    k12_ir_base <- 60 / mgt_ss

    # Logistic carrying-capacity attenuation on the S1 -> S2
    # transition (linear-scale form matching Wicha 2017 and the
    # standard Bulitta LCGM). Scaling the slow transition rather
    # than the doubling step makes the no-drug asymptote land at
    # CFUmax = 10^logcfumax (with k21*S2 the only growth term, the
    # system stops growing when k12_eff -> 0 i.e. cfu_total ->
    # CFUmax). The exact attenuation site is not stated in the
    # paper main text and the Fig S2 schematic was not on disk --
    # see the vignette Assumptions and deviations section.
    cfu_total <- s1_ss + s2_ss + s1_ri + s2_ri + s1_ir + s2_ir
    growth_attn <- 1 - cfu_total / (10 ^ logcfumax)
    k12_ss <- k12_ss_base * growth_attn
    k12_ri <- k12_ri_base * growth_attn
    k12_ir <- k12_ir_base * growth_attn

    # =============================================================
    # Mechanistic synergy switch for population 3 (IPM-I/TOB-R)
    # =============================================================
    # Paper text (Results): "the KC50,IR,IPM was 112 mg/L in the
    # absence of and 1.60 mg/L in the presence of at least 1.15
    # mg/L tobramycin". Discrete switch at tob_cut.
    kc50_ir_ipm_eff <- kc50_ir_ipm_mono
    if (CONC_TOB_MGL >= tob_cut) kc50_ir_ipm_eff <- kc50_ir_ipm_combo

    # =============================================================
    # Per-subpopulation killing rates (1/h)
    # =============================================================
    # Imipenem: sigmoidal Hill (Emax fixed at Kmax,IPM, shared
    # across subpopulations per Table 1).
    ipm_hp <- CONC_IPM_MGL ^ hill_ipm
    kill_ss_ipm <- kmax_ipm * ipm_hp / (kc50_ss_ipm ^ hill_ipm + ipm_hp)
    kill_ri_ipm <- kmax_ipm * ipm_hp / (kc50_ri_ipm ^ hill_ipm + ipm_hp)
    kill_ir_ipm <- kmax_ipm * ipm_hp / (kc50_ir_ipm_eff ^ hill_ipm + ipm_hp)

    # Tobramycin: Emax (Hill = 1; no Hill exponent is tabulated).
    # Kmax for population 3 (IR) is held equal to population 1 (SS)
    # per Table 1 footnote (e).
    kill_ss_tob <- kmax_tob_ss * CONC_TOB_MGL / (kc50_ss_tob + CONC_TOB_MGL)
    kill_ri_tob <- kmax_tob_ri * CONC_TOB_MGL / (kc50_ri_tob + CONC_TOB_MGL)
    kill_ir_tob <- kmax_tob_ss * CONC_TOB_MGL / (kc50_ir_tob + CONC_TOB_MGL)

    # Total per-state killing rate per subpopulation (additive).
    kill_ss <- kill_ss_ipm + kill_ss_tob
    kill_ri <- kill_ri_ipm + kill_ri_tob
    kill_ir <- kill_ir_ipm + kill_ir_tob

    # =============================================================
    # ODE system -- Bulitta two-state LCGM per subpopulation
    # =============================================================
    # For each population p in {ss, ri, ir}:
    #   d/dt(S1_p) = -k12_p_eff * S1_p
    #                + 2 * k21 * S2_p          (doubling)
    #                - kill_p * S1_p
    #   d/dt(S2_p) =  k12_p_eff * S1_p - k21 * S2_p
    #                - kill_p * S2_p
    # where k12_p_eff = k12_p_base * growth_attn carries the
    # carrying-capacity logistic limit. Killing acts on both states
    # (standard assumption; the supplement Fig. S2 / S4 schematic
    # was not on disk -- see the vignette Assumptions and
    # deviations section).

    d/dt(s1_ss) <- -k12_ss * s1_ss + 2 * k21 * s2_ss - kill_ss * s1_ss
    d/dt(s2_ss) <-  k12_ss * s1_ss - k21 * s2_ss - kill_ss * s2_ss

    d/dt(s1_ri) <- -k12_ri * s1_ri + 2 * k21 * s2_ri - kill_ri * s1_ri
    d/dt(s2_ri) <-  k12_ri * s1_ri - k21 * s2_ri - kill_ri * s2_ri

    d/dt(s1_ir) <- -k12_ir * s1_ir + 2 * k21 * s2_ir - kill_ir * s1_ir
    d/dt(s2_ir) <-  k12_ir * s1_ir - k21 * s2_ir - kill_ir * s2_ir

    # =============================================================
    # Initial conditions
    # =============================================================
    # Total inoculum 10^logcfu0 is partitioned across the three
    # subpopulations using the mutation frequencies, then split
    # between S1 and S2 at the LCGM pseudo-steady-state ratio
    # (k12*S1 = k21*S2 implies S2/total = k12/(k12+k21)). This
    # avoids the otherwise large initial S1 -> S2 transient that a
    # naive "all in S1" initialisation would impose during the
    # first ~1/k21 ~= 1.2 min and would distort the very dense
    # early sampling of the HFIM.
    total0 <- 10 ^ logcfu0
    frac_ri <- 10 ^ log10_mut_ipm
    frac_ir <- 10 ^ log10_mut_tob
    frac_ss <- 1 - frac_ri - frac_ir

    f_s2_ss <- k12_ss_base / (k12_ss_base + k21)
    f_s2_ri <- k12_ri_base / (k12_ri_base + k21)
    f_s2_ir <- k12_ir_base / (k12_ir_base + k21)

    s1_ss(0) <- total0 * frac_ss * (1 - f_s2_ss)
    s2_ss(0) <- total0 * frac_ss * f_s2_ss
    s1_ri(0) <- total0 * frac_ri * (1 - f_s2_ri)
    s2_ri(0) <- total0 * frac_ri * f_s2_ri
    s1_ir(0) <- total0 * frac_ir * (1 - f_s2_ir)
    s2_ir(0) <- total0 * frac_ir * f_s2_ir

    # =============================================================
    # Observation: log10 total viable count (CFU/mL)
    # =============================================================
    # 1 CFU/mL floor mirrors the experimental limit of counting
    # (paper: counts below 1.0 log10 CFU/mL plotted as zero).
    # Variable named Cc per nlmixr2lib single-output convention;
    # values are log10 CFU/mL (not a drug concentration). Matches
    # Wicha 2017 pharmacodynamics/ pattern.
    cfu_obs_floor <- cfu_total + 1
    Cc <- log10(cfu_obs_floor)
    Cc ~ add(addSd)
  })
}
