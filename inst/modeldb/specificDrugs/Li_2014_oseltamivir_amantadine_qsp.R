Li_2014_oseltamivir_amantadine_qsp <- function() {
  description <- "QSP. Full mechanistic model (the '+IM' model) of influenza A virus infection in an average adult, with innate (IFN-alpha) and cell-mediated adaptive (CD8+ cytotoxic T lymphocyte) immune responses, driven by oral oseltamivir phosphate (OP) and oral amantadine (AMT). Each drug has a one-compartment PK module with a first-order dissolution delay and first-order absorption; OP is converted first-order to oseltamivir carboxylate (OC). The viral life-cycle module (uninfected, IFN-refractory, eclipse-phase and virus-producing lung epithelial cells plus free virus) is inhibited at the replication step by AMT (1/(1 + C_AMT/IC50_AMT)) and at the virus-release step by OC (1/(1 + C_OC/IC50_OC)). Infected lung antigen-presenting cells (APC) migrate to the lymph node after a 24-h delay, activate CD8+ T cells, and 15% of the lymph-node T-cell count reaches the lung after a further 12-h delay, where it kills infected cells and induces IFN-alpha secretion; IFN-alpha drives NK-mediated killing of infected cells and an antiviral refractory state in uninfected cells. A delay differential equation model (rxode2 delay()); solve with method = 'dop853' and a large maxsteps. Deterministic typical-value model: no between-subject variability and no residual error. Time 0 is the initial virus exposure; the paper starts treatment 96 h later. Rate constant K37 (APC migration to the lymph node) is not tabulated by Li 2014 and was back-solved by the maintainers from the paper's own simulated curves (Figure 3a; see the vignette). The companion target-cell-limited model without immune responses is modellib('Li_2014_oseltamivir_amantadine_noimmune_qsp')."
  reference <- paste(
    "Li Z, Zhou H, Lu Y, Colatsky T. (2014). A Critical Role for Immune System",
    "Response in Mediating Anti-influenza Drug Synergies Assessed by Mechanistic",
    "Modeling. CPT Pharmacometrics Syst Pharmacol 3(9):e135.",
    "doi:10.1038/psp.2014.32. Equations from Supplementary Table S1, parameter",
    "values from Supplementary Table S2 (+IM column), initial values from",
    "Supplementary Table S3, and the authors' deSolve implementation of the",
    "immune-response model (Supplementary file delaymodelfun.R)."
  )
  vignette <- "Li_2014_oseltamivir_amantadine_influenza"

  # Oral oseltamivir phosphate is dosed into `depot` and oral amantadine into
  # `depot_amantadine`.
  dosing <- c("depot", "depot_amantadine")

  # Viral life-cycle and immune-cell states with no canonical analogue in
  # inst/references/compartment-names.md. `target`, `infected`, `virus` and
  # `ifna` are canonical; the amantadine PK states carry a sibling-drug suffix
  # that is not registered as a metabolite suffix.
  paper_specific_compartments <- c(
    "depot_amantadine",
    "transit1_amantadine",
    "central_amantadine",
    "refractory",
    "eclipse",
    "apc_naive",
    "apc_infected",
    "apc_activated",
    "ctl"
  )

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL for the plasma drug concentrations Cc, Cc_oselcarb and Cc_amantadine; TCID50/mL for virus; pg/mL for ifna; cells for the epithelial and immune-cell states"
  )

  # verified = TRUE -- each state checked against Li 2014 Supplementary
  # Table S1 (equations), Table S3 (initial values and units) and Figure 1.
  compartmentData <- list(
    depot = list(analyte = "oseltamivir phosphate", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "oseltamivir phosphate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oseltamivir phosphate", units = "mg", specimen = "plasma", verified = TRUE),
    central_oselcarb = list(analyte = "oseltamivir carboxylate", units = "mg", specimen = "plasma", verified = TRUE),
    depot_amantadine = list(analyte = "amantadine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1_amantadine = list(analyte = "amantadine", units = "mg", specimen = "administration site", verified = TRUE),
    central_amantadine = list(analyte = "amantadine", units = "mg", specimen = "plasma", verified = TRUE),
    target = list(
      analyte = "uninfected susceptible lung epithelial cells (E1)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    refractory = list(
      analyte = "IFN-induced refractory (protected) lung epithelial cells (E2)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    eclipse = list(
      analyte = "lung epithelial cells carrying influenza A virus particles, eclipse phase (E3)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    infected = list(
      analyte = "lung epithelial cells with replicated influenza A virus, virus-producing (E4)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus = list(analyte = "free influenza A virus", units = "TCID50/mL", specimen = "tissue", verified = TRUE),
    apc_naive = list(
      analyte = "inactivated antigen-presenting cells in the lung (pulmonary dendritic cells)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    apc_infected = list(
      analyte = "virus-infected antigen-presenting cells in the lung",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    apc_activated = list(
      analyte = "activated antigen-presenting cells in the lymph node",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    ctl = list(
      analyte = "virus-specific CD8+ cytotoxic T lymphocytes in the lymph node",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    ifna = list(analyte = "interferon-alpha (IFNA1)", units = "pg/mL", specimen = "tissue", verified = TRUE)
  )

  # A typical-value model for an "average" adult; no subject-level covariates.
  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = "Typical ('average') adult calibrated against published summary data; no individual-level data",
    n_studies = 7L,
    disease_state = "Acute uncomplicated influenza A virus infection (experimental and natural); healthy volunteers for the PK modules",
    dose_range = "Oseltamivir phosphate 75-150 mg twice daily and amantadine 100 mg twice daily in the calibration data; 0-400 mg twice daily of each in the combination simulations",
    notes = paste(
      "Li 2014 calibrated the model to published summary data rather than fitting individual data.",
      "PK modules: oseltamivir loading-dose (150 + 75 mg) multidose study in Thai volunteers (Wattanagoon 2009),",
      "150 mg single dose in Japanese and Caucasian volunteers (Jhee 2008), amantadine 100 mg once daily",
      "(Reuman 1989) and 200 mg single dose averaged over young and elderly adults (Hayden 1985) -- Figure 2.",
      "Immune modules: nasal-lavage viral titre and IFN-alpha after experimental influenza A infection",
      "(Hayden 1998) -- Figure 3. Treatment calibration: oseltamivir 75 and 150 mg twice daily versus",
      "placebo (Treanor 2000) and amantadine versus placebo (Knight 1970) -- Figure 4.",
      "Parameters were estimated with a genetic algorithm (GAPSM) within literature-based ranges",
      "(Supplementary Table S2)."
    )
  )

  ini({
    # ==================================================================
    # PK modules -- Li 2014 Supplementary Table S2 (shared by the +IM and
    # -IM models). Clearances and volumes are tabulated in mL/h and mL for
    # a 70-kg adult and are kept in those units; model() converts the
    # mg/mL ratio of a dose in mg to ng/mL.
    # ==================================================================
    lktr <- log(6.412)
    label("OP dissolution (delay) rate constant from depot to transit1, K1a (1/h)") # Table S2: K1a = 6.412
    lka <- log(0.3652)
    label("OP absorption rate constant from transit1 to central, K1b (1/h)") # Table S2: K1b = 0.3652
    lkm <- log(1.768)
    label("OP to OC first-order conversion rate constant, K2 (1/h)") # Table S2: K2 = 1.768
    lcl <- log(78060)
    label("Apparent OP clearance, CLOP (mL/h)") # Table S2: CLOP = 78060
    lvc <- log(212600)
    label("Apparent OP volume of distribution, VDOP (mL)") # Table S2: VDOP = 212600
    lcl_oselcarb <- log(26460)
    label("Apparent OC clearance, CLOC (mL/h)") # Table S2: CLOC = 26460
    lvc_oselcarb <- log(115000)
    label("Apparent OC volume of distribution, VDOC (mL)") # Table S2: VDOC = 115000
    lktr_amantadine <- log(62.69)
    label("AMT dissolution (delay) rate constant from depot to transit1, K19 (1/h)") # Table S2: K19 = 62.69
    lka_amantadine <- log(3.356)
    label("AMT absorption rate constant from transit1 to central, K20 (1/h)") # Table S2: K20 = 3.356
    lcl_amantadine <- log(12850)
    label("Apparent AMT clearance, CLAMT (mL/h)") # Table S2: CLAMT = 12850
    lvc_amantadine <- log(325900)
    label("Apparent AMT volume of distribution, VDAMT (mL)") # Table S2: VDAMT = 325900

    # ==================================================================
    # Drug effects -- Li 2014 Supplementary Table S2 (+IM value, the first
    # of each '+IM; -IM' pair) and Table S1 footnotes b and d. The paper
    # parameterises each inhibition by 1/IC50 (K17 for OC, K22 for AMT);
    # the IC50 is encoded here as log(1/K).
    # ==================================================================
    lic50_oselcarb <- log(1 / 0.0445)
    label("In vivo plasma OC IC50 for inhibition of virus release, 1/K17 (ng/mL)") # Table S2: K17 = 0.0445 (+IM)
    lic50_amantadine <- log(1 / 0.07987)
    label("In vivo plasma AMT IC50 for inhibition of viral replication, 1/K22 (ng/mL)") # Table S2: K22 = 0.07987 (+IM)

    # ==================================================================
    # Virus life-cycle module -- Li 2014 Supplementary Table S2, +IM value.
    # ==================================================================
    lbeta <- log(1.04e-5)
    label("Epithelial cell infection rate constant, K11 (1/(TCID50/mL)/h)") # Table S2: K11 = 1.04*10-5 (+IM)
    lk_eclipse <- log(0.06191)
    label("Eclipse-phase transition rate constant to virus-producing cells, K13 (1/h)") # Table S2: K13 = 0.06191 (+IM)
    ldelta <- log(0.4196)
    label("Death rate constant of virus-producing cells, K14 (1/h)") # Table S2: K14 = 0.4196 (+IM)
    lp_virus <- log(0.00231)
    label("Virus release rate per virus-producing cell, K16 ((TCID50/mL)/cell/h)") # Table S2: K16 = 0.00231 (+IM)
    lc_virus <- log(0.6007)
    label("Free-virus elimination rate constant, K18 (1/h)") # Table S2: K18 = 0.6007 (+IM)

    # ==================================================================
    # Lung and lymphatic immune-response modules -- Li 2014 Supplementary
    # Table S2.
    # ==================================================================
    lk_ifn_ctl <- log(0.01439)
    label("CTL-induced IFN-alpha secretion rate, K24 ((pg/mL)/cell/h)") # Table S2: K24 = 0.01439
    lkdeg_ifn <- log(1.856)
    label("IFN-alpha decay rate constant, K25 (1/h)") # Table S2: K25 = 1.856
    lk_protect <- log(0.0003995)
    label("IFN-induced protection (refractory-state) rate constant, K26 (1/(pg/mL)/h)") # Table S2: K26 = 0.0003995
    lk_revert <- log(8.066e-5)
    label("Reversion rate constant from the refractory state, K27 (1/h)") # Table S2: K27 = 8.066*10-5
    lk_nk <- log(1.614e-7)
    label("IFN-induced NK killing rate constant of infected cells, K28 (1/(pg/mL)/h)") # Table S2: K28 = 1.614*10-7
    lk_ctl_act <- log(0.1135)
    label("Activation rate constant of CD8+ T cells by lymph-node APC, K32 (1/h)") # Table S2: K32 = 0.1135
    lk_ctl_kill <- log(5.002e-6)
    label("CTL killing rate constant of infected cells, K34 (1/cell/h)") # Table S2: K34 = 5.002*10-6
    # K37 appears in the Table S1 equation for activated lymph-node APC but
    # has no row in Table S2. Back-solved by the maintainers from the
    # paper's simulated viral-titre and IFN-alpha curves of Figure 3a
    # (least squares on the drawn model line), and cross-checked against
    # the IFN-alpha second-to-first peak ratio of Figure 3b (0.18 vs ~0.19)
    # and the +IM curves of Figure 4. Only the product K32*K37 is
    # identifiable from any model output.
    lk_apc_mig <- log(0.111)
    label("Rate constant of lymph-node APC activation by migrated infected lung APC, K37 (1/h)") # back-solved from Figure 3a; not tabulated
    lkdeg_apc_act <- log(0.07744)
    label("Loss rate constant of activated APC in the lymph node, K42 (1/h)") # Table S2: K42 = 0.07744
    lkdeg_ctl <- log(0.04555)
    label("CTL death rate constant, K45 (1/h)") # Table S2: K45 = 0.04555
    lk_ifn_inf <- log(5.45e-14)
    label("IFN-alpha secretion rate from virus-producing epithelial cells, K47 ((pg/mL)/cell/h)") # Table S2: K47 = 5.45*10-14
    lk_apc_inf <- log(0.004715)
    label("APC infection rate constant, K48 (1/(TCID50/mL)/h)") # Table S2: K48 = 0.004715
    lkdeg_apc_inf <- log(0.09598)
    label("Death rate constant of infected APC in the lung, K53 (1/h)") # Table S2: K53 = 0.09598
    lk_ifn_apc <- log(0.01177)
    label("IFN-alpha secretion rate from infected APC, K57 ((pg/mL)/cell/h)") # Table S2: K57 = 0.01177
    lk_apc_turn <- log(6.459e-7)
    label("Turnover rate constant of uninfected lung APC, K59 (1/h)") # Table S2: K59 = 6.459*10-7

    # ==================================================================
    # Assumed constants -- Li 2014 Methods 'Simulation of immune responses
    # and therapeutic effects' and Table S1 footnotes c and f.
    # ==================================================================
    f_ctl_lung <- fixed(0.15)
    label("CTL count in the lung as a fraction of the lymph-node count (fraction)") # Methods; Table S1 footnote c
    tau_ctl_mig <- fixed(12)
    label("Delay of CTL migration from lymph node to lung (h)") # Methods; Table S1 footnote c
    tau_apc_mig <- fixed(24)
    label("Delay of infected-APC migration from lung to lymph node (h)") # Methods; Table S1 footnote f
    apc0 <- fixed(1e5)
    label("Balance (and initial) level of inactivated lung APC (cells)") # Methods; Table S1 'K59*(1E5 - A)'; Table S3: A = 1*105
    target0 <- fixed(4e8)
    label("Initial uninfected lung epithelial cells (cells)") # Table S3: EpitheliumLung = 4*108
    # Initial free-virus titre. Li 2014 fitted it separately for each
    # simulated data set (Table S3: 1.474*10-5 for Figure 3; 4.136*10-7,
    # 1.849*10-7, 2.486*10-6 for the Figure 4a placebo, 75-mg and 150-mg
    # groups; 5.541*10-9 for the Figure 4b placebo group and Figure 4c;
    # 1.005*10-10 for the Figure 4b amantadine group) and fixed it at 10
    # for the drug-combination simulations.
    virus0 <- fixed(10)
    label("Initial free influenza A virus titre at exposure (TCID50/mL)") # Methods 'Calculation of synergy': initial titre of 10
  })

  model({
    ktr <- exp(lktr)
    ka <- exp(lka)
    km <- exp(lkm)
    cl <- exp(lcl)
    vc <- exp(lvc)
    cl_oselcarb <- exp(lcl_oselcarb)
    vc_oselcarb <- exp(lvc_oselcarb)
    ktr_amantadine <- exp(lktr_amantadine)
    ka_amantadine <- exp(lka_amantadine)
    cl_amantadine <- exp(lcl_amantadine)
    vc_amantadine <- exp(lvc_amantadine)
    ic50_oselcarb <- exp(lic50_oselcarb)
    ic50_amantadine <- exp(lic50_amantadine)
    beta <- exp(lbeta)
    k_eclipse <- exp(lk_eclipse)
    delta <- exp(ldelta)
    p_virus <- exp(lp_virus)
    c_virus <- exp(lc_virus)
    k_ifn_ctl <- exp(lk_ifn_ctl)
    kdeg_ifn <- exp(lkdeg_ifn)
    k_protect <- exp(lk_protect)
    k_revert <- exp(lk_revert)
    k_nk <- exp(lk_nk)
    k_ctl_act <- exp(lk_ctl_act)
    k_ctl_kill <- exp(lk_ctl_kill)
    k_apc_mig <- exp(lk_apc_mig)
    kdeg_apc_act <- exp(lkdeg_apc_act)
    kdeg_ctl <- exp(lkdeg_ctl)
    k_ifn_inf <- exp(lk_ifn_inf)
    k_apc_inf <- exp(lk_apc_inf)
    kdeg_apc_inf <- exp(lkdeg_apc_inf)
    k_ifn_apc <- exp(lk_ifn_apc)
    k_apc_turn <- exp(lk_apc_turn)

    # ------------------------------------------------------------------
    # PK modules -- Table S1. 'OP Oral' -> 'OP Delay' -> 'OP Central' ->
    # 'OC Central', and 'AMT Oral' -> 'AMT Delay' -> 'AMT Central'. The OP
    # to OC conversion carries no molecular-weight correction, as printed.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ktr * depot
    d/dt(transit1) <- ktr * depot - ka * transit1
    d/dt(central) <- ka * transit1 - km * central - cl / vc * central
    d/dt(central_oselcarb) <- km * central - cl_oselcarb / vc_oselcarb * central_oselcarb
    d/dt(depot_amantadine) <- -ktr_amantadine * depot_amantadine
    d/dt(transit1_amantadine) <- ktr_amantadine * depot_amantadine - ka_amantadine * transit1_amantadine
    d/dt(central_amantadine) <- ka_amantadine * transit1_amantadine - cl_amantadine / vc_amantadine * central_amantadine

    # Amounts in mg over volumes in mL give mg/mL; x 1e6 gives ng/mL, the
    # units of the 1/IC50 parameters K17 and K22.
    Cc <- central / vc * 1e6
    Cc_oselcarb <- central_oselcarb / vc_oselcarb * 1e6
    Cc_amantadine <- central_amantadine / vc_amantadine * 1e6

    # Table S1 footnotes b and d: inhibitory saturation models with a
    # maximum effect of 1.
    inh_replication <- 1 / (1 + Cc_amantadine / ic50_amantadine)
    inh_release <- 1 / (1 + Cc_oselcarb / ic50_oselcarb)

    # Delayed immune signals (Table S1 footnotes c and f). Before the delay
    # has elapsed the history is the zero initial condition, matching the
    # authors' deSolve code (lag set to 0 for t < 12 h and t < 24 h).
    ctl_lung <- f_ctl_lung * delay(ctl, tau_ctl_mig)
    apc_migrated <- delay(apc_infected, tau_apc_mig)

    # ------------------------------------------------------------------
    # Virus life-cycle module -- Table S1 (model with immune systems).
    # E1 = target, E2 = refractory, E3 = eclipse, E4 = infected.
    # ------------------------------------------------------------------
    d/dt(target) <- -beta * virus * target + k_revert * refractory - k_protect * ifna * target
    d/dt(refractory) <- -k_revert * refractory + k_protect * ifna * target
    d/dt(eclipse) <- beta * virus * target - inh_replication * k_eclipse * eclipse -
      k_ctl_kill * eclipse * ctl_lung - k_nk * ifna * eclipse
    d/dt(infected) <- inh_replication * k_eclipse * eclipse - delta * infected -
      k_ctl_kill * infected * ctl_lung - k_nk * ifna * infected
    d/dt(virus) <- inh_release * p_virus * infected - c_virus * virus

    # ------------------------------------------------------------------
    # Lung and lymphatic immune-response modules -- Table S1.
    # ------------------------------------------------------------------
    d/dt(apc_naive) <- k_apc_turn * (apc0 - apc_naive) - k_apc_inf * virus * apc_naive
    d/dt(apc_infected) <- k_apc_inf * virus * apc_naive - kdeg_apc_inf * apc_infected
    d/dt(apc_activated) <- k_apc_mig * apc_migrated - kdeg_apc_act * apc_activated
    d/dt(ctl) <- k_ctl_act * apc_activated - kdeg_ctl * ctl
    d/dt(ifna) <- k_ifn_ctl * ctl_lung + k_ifn_inf * infected + k_ifn_apc * apc_infected - kdeg_ifn * ifna

    target(0) <- target0
    apc_naive(0) <- apc0
    virus(0) <- virus0

    # Li 2014 plots titres as log10(titre + 1): 'a pseudo count of 1 is
    # added to each viral titer data point ... for the purpose of log
    # transformation' (Table S3 note).
    log10_viral_load <- log10(virus + 1)
    # Fraction of epithelial cells not infected (uninfected plus
    # refractory), the dashed curve of Figure 4c.
    frac_uninfected <- (target + refractory) / target0
  })
}
