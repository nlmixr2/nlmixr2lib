Yadav_2017_imipenem_tobramycin <- function() {
  description <- "Preclinical (mouse, neutropenic murine thigh infection model; clinical Pseudomonas aeruginosa isolate FADDI-PA088). Mechanism-based pharmacodynamic model for an imipenem-plus-tobramycin combination with humanized dosing schemes. Two pre-existing bacterial subpopulations differing in imipenem and tobramycin susceptibility (population 1 IPM-S/TOB-R, population 2 IPM-I/TOB-R) each follow a Bulitta two-state life-cycle growth model (S1 -> S2 -> 2*S1 with replication rate k21 fixed; mean generation time 1/k12 shared across the two subpopulations). A plateau factor PLAT = 1 - CFUall/(CFUall + CFUmax) attenuates the replication step so the no-drug viable count plateaus at CFUmax. Imipenem kills with a sigmoidal Hill function whose effective KC50 is multiplied by the mechanistic-synergy factor OM_effect = 1 - Imax,OM,TOB * Ctob / (Ctob + IC50,OM,TOB) (i.e. tobramycin permeabilizing the outer bacterial membrane toward imipenem). Both subpopulations are tobramycin-resistant, so tobramycin has no direct killing term in this model; its only role is to lower the effective imipenem KC50 via OM_effect. Imipenem and tobramycin concentrations are external time-varying inputs (covariates Cipm and Ctob); the model contains no rodent PK component (the paper imported the murine one-compartment PK of imipenem and tobramycin from external references and the PK parameter values are not reported in the source on disk)."
  reference <- paste(
    "Yadav R, Bulitta JB, Wang J, Nation RL, Landersdorfer CB. (2017).",
    "Evaluation of pharmacokinetic/pharmacodynamic model-based optimized",
    "combination regimens against multidrug-resistant Pseudomonas aeruginosa",
    "in a murine thigh infection model by using humanized dosing schemes.",
    "Antimicrobial Agents and Chemotherapy 61(12):e01268-17.",
    "doi:10.1128/AAC.01268-17.",
    sep = " "
  )
  vignette <- "Yadav_2017_imipenem_tobramycin"
  units <- list(time = "hour", dosing = "mg/L (drug input concentration)", concentration = "log10 CFU/thigh (observation); mg/L (drug covariates)")

  # Cipm / Ctob are the time-varying imipenem / tobramycin plasma
  # unbound concentrations supplied externally from a PK driver (the
  # paper used the published murine one-compartment PK of each drug,
  # equations 1-3 with parameter values imported from references 23
  # [Katsube 2008, imipenem] and 24 [Moffie 1993, tobramycin]; those
  # numerical PK values are not reported in the present paper on
  # disk). Declared in `depends` so the canonical-covariate register
  # check does not apply (they are documented in covariateData below
  # for provenance).
  depends <- c("Cipm", "Ctob")

  # bact_* state names are paper-mechanistic life-cycle states for the
  # two-subpopulation Bulitta growth model (Fig. 3) and are explicitly
  # registered for the convention checker via the regex pattern.
  paper_specific_compartment_pattern <- "^bact_"

  covariateData <- list(
    Cipm = list(
      description        = "Unbound imipenem plasma concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally; the murine experiment delivered humanized exposure via imipenem 60 or 77 mg/kg s.c. every 2 h (total daily doses 720 or 924 mg/kg), targeting the plasma unbound concentration profiles of 4 or 5 g/day continuous infusion (with a 1 g loading dose) in critically ill patients. The driver PK is the one-compartment s.c. model of equations 1-3 with ka, ke, V/F and fu from the murine imipenem PK reference (Katsube 2008) -- numerical PK parameter values not reported in the present paper on disk. Preclinical experimental input -- not in inst/references/covariate-columns.md (the canonical register is for human pop-PK covariates and does not apply to this in-vivo murine PD model).",
      source_name        = "Imipenem concentration (paper Methods, Eqs 1-3; Fig. 1 legend)"
    ),
    Ctob = list(
      description        = "Unbound tobramycin plasma concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate supplied externally; the murine experiment delivered humanized exposure via tobramycin s.c. fractionated doses (33.3 percent at 0 h; 16.65 percent at 4, 8, and 12 h; 8.37 percent at 16 and 20 h; total daily 73 mg/kg) targeting the plasma unbound profile of 7 mg/kg q24h as a 0.5 h infusion in humans. The driver PK is the one-compartment s.c. model of equations 1-3 with ka, ke, V/F and fu from the murine tobramycin PK reference (Moffie 1993) -- numerical PK parameter values not reported in the present paper on disk. Preclinical experimental input -- not in inst/references/covariate-columns.md.",
      source_name        = "Tobramycin concentration (paper Methods, Eqs 1-3; Fig. 1 legend)"
    )
  )

  population <- list(
    species             = "mouse (Swiss, 7-week-old males, 25-30 g; neutropenic via cyclophosphamide 150 mg/kg i.p. 4 days pre-infection plus 100 mg/kg 1 day pre-infection)",
    n_subjects          = NA_integer_,
    n_studies           = 2L,
    organism            = "Pseudomonas aeruginosa FADDI-PA088 (carbapenem- and aminoglycoside-resistant clinical isolate; MIC imipenem 16 mg/L, MIC tobramycin 32 mg/L)",
    system              = "Neutropenic murine thigh infection model; two mice (four thighs) per dose-regimen and time-point combination; viable counts on antibiotic-free and antibiotic-containing (3x MIC) CAMHA agar at 2, 6, and 24 h post-treatment-initiation",
    inoculum            = "approximately 10^5 CFU/thigh (50 uL of 2x10^6 CFU/mL bacterial suspension injected into each posterior thigh muscle under isoflurane anaesthesia 2 h before treatment start)",
    mic_values          = c(imipenem = "16 mg/L", tobramycin = "32 mg/L"),
    duration            = "24 h (sampling at 0, 2, 6, and 24 h post-treatment-initiation)",
    regimens            = "Untreated control; tobramycin monotherapy (humanized 7 mg/kg q24h 0.5-h infusion); imipenem monotherapy at humanized 4 or 5 g/day continuous infusion with a 1 g loading dose (60 or 77 mg/kg s.c. every 2 h in mice); imipenem + tobramycin combinations using the same humanized regimens",
    notes               = "Experiment 1 mimicked imipenem 4 g/day plus tobramycin (Log CFU0 4.93); experiment 2 mimicked imipenem 5 g/day plus tobramycin (Log CFU0 4.78). The packaged ini() uses experiment 1's Log CFU0 = 4.93 as the default; switch to 4.78 to reproduce experiment 2. The MBM was fit in S-ADAPT using importance sampling (pmethod = 4); the coefficient of correlation for observed versus population-fitted log10 viable counts was at least 0.98. Random effects (eta) are NOT estimated in the source -- the paper reports population mean parameter estimates with relative standard errors only (between-curves variability was fixed to a final small CV), so the packaged model contains no etas and is intended for typical-value simulation only. See Yadav 2017 Methods (page 8 = e01268-17 p. 8) and Table 1."
  )

  ini({
    # =============================================================
    # Bacterial life-cycle (Bulitta two-state LCGM)
    # =============================================================
    # k21 = replication rate constant (S2 -> 2*S1 doubling step;
    # 1/k12 = mean generation time for the slow S1 -> S2 transition.
    # The replication rate constant is shared across subpopulations
    # and fixed at 50 /h per Table 1 ("Replication rate constant
    # (h-1) k21 = 50 (fixed)").
    lk21 <- fixed(log(50))
    label("Replication rate constant k21 (1/h; S2 -> 2*S1 doubling; FIXED)")  # Yadav 2017 Table 1

    # Mean generation time (minutes) for the slow S1 -> S2 transition.
    # Table 1 reports k12,SR^-1 = 142 (min) and k12,IR^-1 = 142 (min)
    # -- identical for the two subpopulations -- with a single SE of
    # 9.21 percent, so a single parameter is encoded.
    mgt <- 142
    label("Mean generation time (min; shared across IPM-S/TOB-R and IPM-I/TOB-R)")  # Yadav 2017 Table 1 (k12,SR^-1 = k12,IR^-1)

    # Initial total inoculum and carrying capacity (log10 CFU/thigh).
    # Experiment 1 (IPM 4 g/day + TOB) Log CFU0 = 4.93; experiment 2
    # (IPM 5 g/day + TOB) Log CFU0 = 4.78. The packaged default uses
    # experiment 1's value; see population$notes.
    logcfu0   <- 4.93
    label("Initial total inoculum (log10 CFU/thigh; experiment 1 default)")  # Yadav 2017 Table 1 (Log CFU0)
    logcfumax <- 8.00
    label("Maximum population size (log10 CFU/thigh)")  # Yadav 2017 Table 1 (LogCFUmax)

    # =============================================================
    # Initial subpopulation fractions (log10 mutation frequencies)
    # =============================================================
    # The IPM-I/TOB-R subpopulation (population 2) has initial
    # fraction 10^log10_mut_ipm of the total; the remainder is the
    # IPM-S/TOB-R subpopulation (population 1). Only one mutation
    # frequency is estimated in this two-population model (for IR).
    log10_mut_ipm <- -3.09
    label("log10 mutation frequency at t = 0 for IPM-I/TOB-R (population 2)")  # Yadav 2017 Table 1 (Log MUT,IPM)

    # =============================================================
    # Imipenem killing (sigmoidal Hill model)
    # =============================================================
    lkmax_ipm <- log(3.38)
    label("Maximum imipenem killing rate Kmax,IPM (1/h)")  # Yadav 2017 Table 1
    lhill_ipm <- log(2.05)
    label("Imipenem Hill coefficient (unitless)")  # Yadav 2017 Table 1 (Hill IPM)

    # Imipenem KC50 per subpopulation.
    lkc50_sr_ipm <- log(0.202)
    label("Imipenem KC50 against IPM-S/TOB-R (population 1; mg/L)")  # Yadav 2017 Table 1 (KC50,SR,IPM)
    lkc50_ir_ipm <- log(96.6)
    label("Imipenem KC50 against IPM-I/TOB-R (population 2; mg/L)")  # Yadav 2017 Table 1 (KC50,IR,IPM)

    # =============================================================
    # Mechanistic synergy via outer-membrane permeabilization
    # =============================================================
    # OM_effect = 1 - Imax,OM,TOB * Ctob / (Ctob + IC50,OM,TOB)
    # (paper equation 7). Effective KC50,IPM = OM_effect * KC50,IPM
    # inside the Hill denominator (paper equations 5-6). At maximal
    # tobramycin the floor of OM_effect is 1 - Imax,OM,TOB = 0.353,
    # giving a maximum 1/(1 - Imax,OM,TOB) = 2.83-fold reduction in
    # KC50,IPM. The 95 percent CI for Imax,OM,TOB is 0.44 - 0.816
    # (Table 1 footnote c); the point estimate 0.647 is used here.
    imax_om_tob  <- 0.647
    label("Maximum fractional reduction of KC50,IPM by tobramycin (unitless; OM permeabilization)")  # Yadav 2017 Table 1 (Imax,OM,TOB)
    lic50_om_tob <- log(2.99)
    label("Tobramycin concentration for half-maximal OM permeabilization (mg/L)")  # Yadav 2017 Table 1 (IC50,OM,TOB)

    # =============================================================
    # Residual error
    # =============================================================
    # Table 1: "SD of residual error on log10 scale" = 0.114. The
    # packaged observation Cc is log10(total viable count) (per the
    # in-vivo PD-model convention shared with Landersdorfer 2018);
    # the additive residual SD enters in log10 units directly.
    addSd <- 0.114
    label("Additive residual SD on log10 total viable count (log10 CFU/thigh)")  # Yadav 2017 Table 1 (SD_CFU)
  })

  model({
    # Effective parameters (exponentiation of log-transformed inputs).
    # No etas are estimated in the source MBM; see population$notes.
    k21          <- exp(lk21)
    kmax_ipm     <- exp(lkmax_ipm)
    hill_ipm     <- exp(lhill_ipm)
    kc50_sr_ipm  <- exp(lkc50_sr_ipm)
    kc50_ir_ipm  <- exp(lkc50_ir_ipm)
    ic50_om_tob  <- exp(lic50_om_tob)

    # Slow S1 -> S2 transition rate (1/h). MGT is in minutes per
    # Table 1, so divide 60 min/h by the minutes-per-generation.
    # Shared across the two subpopulations per Table 1.
    k12 <- 60 / mgt

    # =============================================================
    # Plateau factor (logistic carrying-capacity attenuation)
    # =============================================================
    # Paper Methods: "The plateau factor (PLAT) is defined as
    # 1 - [CFUall / (CFUall + CFUmax)]". Applied to the doubling
    # step 2 * k21 * S2 (paper equation 5). At no-drug steady state
    # PLAT = 0.5, which makes the total viable count plateau at
    # exactly CFUmax (algebra: 2*k21*S2*PLAT = k21*S2 implies
    # PLAT = 0.5 -> CFUall = CFUmax).
    cfu_total <- bact_susceptible_resistant1 + bact_susceptible_resistant2 + bact_intermediate_resistant1 + bact_intermediate_resistant2
    cfumax    <- 10 ^ logcfumax
    plat      <- 1 - cfu_total / (cfu_total + cfumax)

    # =============================================================
    # Mechanistic synergy (outer-membrane permeabilization)
    # =============================================================
    # Paper equation 7: OM_effect = 1 - Imax,OM,TOB * Ctob / (Ctob + IC50,OM,TOB).
    # Applied to BOTH subpopulations (paper Results: "Mechanistic
    # synergy was expressed as a decrease in the imipenem
    # concentration causing 50 percent killing (KC50,IPM) of both
    # bacterial populations with increasing tobramycin
    # concentrations").
    om_effect <- 1 - imax_om_tob * Ctob / (Ctob + ic50_om_tob)

    # =============================================================
    # Per-subpopulation imipenem killing rate (1/h)
    # =============================================================
    # Paper equations 5-6: kill = Kmax,IPM * Cipm^Hill /
    #   (Cipm^Hill + (OM_effect * KC50,IPM)^Hill). Tobramycin has no
    # direct killing term -- both populations are tobramycin-resistant
    # (TOB-R), so tobramycin's only contribution is the OM_effect
    # multiplier on KC50,IPM.
    ipm_hp        <- Cipm ^ hill_ipm
    kc50_sr_eff_h <- (om_effect * kc50_sr_ipm) ^ hill_ipm
    kc50_ir_eff_h <- (om_effect * kc50_ir_ipm) ^ hill_ipm
    kill_sr       <- kmax_ipm * ipm_hp / (ipm_hp + kc50_sr_eff_h)
    kill_ir       <- kmax_ipm * ipm_hp / (ipm_hp + kc50_ir_eff_h)

    # =============================================================
    # ODE system -- Bulitta two-state LCGM per subpopulation
    # =============================================================
    # For each population p in {sr, ir} (paper Fig. 3, equations 5-6):
    #   d/dt(S1_p) = 2 * PLAT * k21 * S2_p
    #                - k12 * S1_p
    #                - kill_p * S1_p
    #   d/dt(S2_p) = k12 * S1_p
    #                - k21 * S2_p
    #                - kill_p * S2_p
    # Killing acts on both states (paper Methods: "Kmax,IPM, KC50,...,
    # and HillIPM affected both states ... of the population").

    d/dt(bact_susceptible_resistant1) <- 2 * plat * k21 * bact_susceptible_resistant2 - k12 * bact_susceptible_resistant1 - kill_sr * bact_susceptible_resistant1
    d/dt(bact_susceptible_resistant2) <- k12 * bact_susceptible_resistant1 - k21 * bact_susceptible_resistant2 - kill_sr * bact_susceptible_resistant2

    d/dt(bact_intermediate_resistant1) <- 2 * plat * k21 * bact_intermediate_resistant2 - k12 * bact_intermediate_resistant1 - kill_ir * bact_intermediate_resistant1
    d/dt(bact_intermediate_resistant2) <- k12 * bact_intermediate_resistant1 - k21 * bact_intermediate_resistant2 - kill_ir * bact_intermediate_resistant2

    # =============================================================
    # Initial conditions
    # =============================================================
    # Total inoculum 10^logcfu0 is partitioned between the two
    # subpopulations using the IR mutation frequency, then split
    # between S1 and S2 at the LCGM pseudo-steady-state ratio
    # (k12 * S1 = k21 * S2 implies S2/total = k12 / (k12 + k21)).
    # This avoids the otherwise large initial S1 -> S2 transient
    # that a naive "all in S1" initialisation would impose during
    # the first ~1/k21 minute. Paper Methods cite the
    # previously-described implementation (refs 40, 63 in Yadav 2017).
    total0  <- 10 ^ logcfu0
    frac_ir <- 10 ^ log10_mut_ipm
    frac_sr <- 1 - frac_ir
    f_s2    <- k12 / (k12 + k21)

    bact_susceptible_resistant1(0)  <- total0 * frac_sr * (1 - f_s2)
    bact_susceptible_resistant2(0)  <- total0 * frac_sr * f_s2
    bact_intermediate_resistant1(0) <- total0 * frac_ir * (1 - f_s2)
    bact_intermediate_resistant2(0) <- total0 * frac_ir * f_s2

    # =============================================================
    # Observation: log10 total viable count (CFU/thigh)
    # =============================================================
    # A 1-CFU floor mirrors the experimental limit of counting and
    # avoids log10(0) if the model drives the total to zero in
    # simulation. Variable named Cc per nlmixr2lib single-output
    # convention; values are log10 CFU/thigh (not a drug
    # concentration), matching the Landersdorfer 2018 in-vitro
    # imipenem-tobramycin PD model in pharmacodynamics/.
    cfu_obs_floor <- cfu_total + 1
    Cc <- log10(cfu_obs_floor)
    Cc ~ add(addSd)
  })
}
