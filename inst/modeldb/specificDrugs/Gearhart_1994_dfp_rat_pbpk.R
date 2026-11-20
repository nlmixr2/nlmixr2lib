Gearhart_1994_dfp_rat_pbpk <- function() {
  description <- paste(
    "Preclinical (rat). PBPK (whole-body, flow-limited, seven tissues plus",
    "venous and arterial blood) with mechanistic acetylcholinesterase (AChE)",
    "inhibition. Diisopropylfluorophosphate (DFP) disposition and the",
    "resulting cholinesterase inhibition, developed by Gearhart et al. 1994",
    "for organophosphate risk assessment. Brain, liver, kidney and richly",
    "perfused tissue are flow-limited with saturable A-esterase (DFPase)",
    "hydrolysis; fat, slowly perfused tissue and diaphragm are flow-limited",
    "storage; venous and arterial blood are explicit compartments that also",
    "hydrolyse DFP (blood A-esterase strongly limits the DFP reaching target",
    "tissue). The lung is an algebraic steady-state gas exchanger coupling",
    "cardiac output and alveolar ventilation through the blood:air partition,",
    "so inhalation enters as an air concentration; iv and sc/im routes enter",
    "as dose records. AChE in brain and in plasma is tracked as a",
    "free/inhibited/aged pool (Fig 2): free enzyme is irreversibly inhibited",
    "by local DFP at a bimolecular rate, inhibited enzyme either regenerates",
    "or ages, and free enzyme turns over by zero-order synthesis and",
    "first-order degradation; the readout is activity as a fraction of",
    "control. Deterministic typical-value model: the publication reports no",
    "IIV and no residual error. The covalent DFP-esterase sink on DFP mass",
    "balance is negligible (authors' Discussion) and is not fed back. The",
    "sc/im absorption rate constant (needed for the rat repeated-sc dosing",
    "data, Figs 6-7) was visually optimised by the authors but not tabulated;",
    "it is a placeholder here (see vignette Errata).",
    sep = " "
  )
  reference <- paste(
    "Gearhart JM, Jepson GW, Clewell HJ 3rd, Andersen ME, Conolly RB.",
    "Physiologically based pharmacokinetic model for the inhibition of",
    "acetylcholinesterase by organophosphate esters.",
    "Environ Health Perspect. 1994 Dec;102 Suppl 11:51-60.",
    "doi:10.1289/ehp.94102s1151. PMID:7737042; PMCID:PMC1566752.",
    "Metabolic constants Table 1; AChE-inhibition parameters Table 3;",
    "physiological parameters Table 4; partition coefficients Table 5;",
    "mass-balance and enzyme ODEs Appendix II.",
    sep = " "
  )
  vignette <- "Gearhart_1994_dfp"

  paper_specific_compartments <- c(
    "a_rich",
    "a_diaphragm",
    "a_hydrolyzed",
    "a_exhaled",
    "ache_free_brain",
    "ache_inh_brain",
    "ache_aged_brain",
    "ache_free_plasma",
    "ache_inh_plasma",
    "ache_aged_plasma"
  )

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L",
    amount = "mg",
    weight = "kg"
  )

  compartmentData <- list(
    depot = list(
      analyte = "Diisopropylfluorophosphate",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    a_venous = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "whole blood", verified = TRUE),
    a_arterial = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "whole blood", verified = TRUE),
    a_brain = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_kidney = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_rich = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_slow = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_diaphragm = list(analyte = "Diisopropylfluorophosphate", units = "mg", specimen = "tissue", verified = TRUE),
    a_hydrolyzed = list(
      analyte = "Diisopropylfluorophosphate, cumulative amount hydrolysed by A-esterase (DFPase)",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_exhaled = list(
      analyte = "Diisopropylfluorophosphate, cumulative amount exhaled",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    ache_free_brain = list(
      analyte = "Acetylcholinesterase, free (active), fraction of baseline activity",
      units = "fraction",
      specimen = "tissue",
      verified = TRUE
    ),
    ache_inh_brain = list(
      analyte = "Acetylcholinesterase, inhibited (reversible), fraction of baseline activity",
      units = "fraction",
      specimen = "tissue",
      verified = TRUE
    ),
    ache_aged_brain = list(
      analyte = "Acetylcholinesterase, aged (permanently inhibited), fraction of baseline activity",
      units = "fraction",
      specimen = "tissue",
      verified = TRUE
    ),
    ache_free_plasma = list(
      analyte = "Acetylcholinesterase, free (active), fraction of baseline activity",
      units = "fraction",
      specimen = "plasma",
      verified = TRUE
    ),
    ache_inh_plasma = list(
      analyte = "Acetylcholinesterase, inhibited (reversible), fraction of baseline activity",
      units = "fraction",
      specimen = "plasma",
      verified = TRUE
    ),
    ache_aged_plasma = list(
      analyte = "Acetylcholinesterase, aged (permanently inhibited), fraction of baseline activity",
      units = "fraction",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (male Fischer-344 for in vitro parameters; male Wistar for the repeated-dose validation)",
    n_subjects = 0L,
    n_studies = 0L,
    disease_state = "healthy",
    dose_range = paste(
      "Validation data: repeated sc DFP in male Wistar rats, first dose",
      "1.1 mg/kg then 0.7 mg/kg every second day for 22 days (Traina &",
      "Serpietri 1984; Michalek et al. 1982). Plasma and brain AChE activity",
      "were simulated (Figures 6-7).",
      sep = " "
    ),
    notes = paste(
      "Deterministic PBPK/PD, not a population fit (n_subjects = 0). DFP",
      "tissue:blood partition coefficients, A-esterase (DFPase) hydrolysis",
      "Vmax/Km and the AChE bimolecular inhibition rate were measured in vitro",
      "in male Fischer-344 rat tissue homogenates (200-250 g); the rat Km was",
      "used directly (Methods; Table 1). Physiological parameters are standard",
      "allometric rat values (Table 4). AChE synthesis, degradation,",
      "regeneration and aging rate constants are from the cited literature",
      "(Table 3 footnotes).",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Physiological parameters -- Table 4, Rat column. Deterministic model;
    # every value is a literature constant, so all are fixed().
    # =====================================================================
    bw <- fixed(0.22)
    label("Body weight (kg)") # Table 4, Body weight, rat = 0.22 kg
    qc <- fixed(4.57)
    label("Cardiac output (L/h)") # Table 4, Cardiac output, rat = 4.57 l/hr
    qp <- fixed(4.57)
    label("Alveolar ventilation (L/h)") # Table 4, Alveolar ventilation, rat = 4.57 l/hr

    # ---- tissue volumes, fraction of body weight (Table 4) --------------
    fv_brain <- fixed(0.0116)
    label("Brain volume, fraction of body weight") # Table 4, Brain, rat = 1.16%
    fv_liver <- fixed(0.0400)
    label("Liver volume, fraction of body weight") # Table 4, Liver, rat = 4.00%
    fv_kidney <- fixed(0.0073)
    label("Kidney volume, fraction of body weight") # Table 4, Kidney, rat = 0.73%
    fv_rich <- fixed(0.0200)
    label("Richly perfused volume, fraction of body weight") # Table 4, Richly perfused, rat = 2.00%
    fv_fat <- fixed(0.0700)
    label("Fat volume, fraction of body weight") # Table 4, Fat, rat = 7.00%
    fv_slow <- fixed(0.6866)
    label("Slowly perfused volume, fraction of body weight") # Table 4, Slowly perfused, rat = 68.66%
    fv_diaphragm <- fixed(0.00030)
    label("Diaphragm volume, fraction of body weight") # Table 4, Diaphragm, rat = 0.030%
    fv_venous <- fixed(0.0400)
    label("Venous blood volume, fraction of body weight") # Table 4, Venous blood, rat = 4.00%
    fv_arterial <- fixed(0.0200)
    label("Arterial blood volume, fraction of body weight") # Table 4, Arterial blood, rat = 2.00%

    # ---- tissue blood flows, fraction of cardiac output (Table 4) -------
    fq_brain <- fixed(0.0300)
    label("Brain blood flow, fraction of cardiac output") # Table 4, Brain, rat = 3.00%
    fq_liver <- fixed(0.2500)
    label("Liver blood flow, fraction of cardiac output") # Table 4, Liver, rat = 25.00%
    fq_kidney <- fixed(0.2000)
    label("Kidney blood flow, fraction of cardiac output") # Table 4, Kidney, rat = 20.00%
    fq_rich <- fixed(0.2796)
    label("Richly perfused blood flow, fraction of cardiac output") # Table 4, Richly perfused, rat = 27.96%
    fq_fat <- fixed(0.0900)
    label("Fat blood flow, fraction of cardiac output") # Table 4, Fat, rat = 9.00%
    fq_slow <- fixed(0.1440)
    label("Slowly perfused blood flow, fraction of cardiac output") # Table 4, Slowly perfused, rat = 14.40%
    fq_diaphragm <- fixed(0.0060)
    label("Diaphragm blood flow, fraction of cardiac output") # Table 4, Diaphragm, rat = 0.60%

    # =====================================================================
    # Partition coefficients -- Table 5, DFP column (species-independent).
    # =====================================================================
    pc_blood_air <- fixed(12.57)
    label("Blood:air partition coefficient (unitless)") # Table 5, DFP blood-air = 12.57
    pc_brain <- fixed(0.67)
    label("Brain:blood partition coefficient (unitless)") # Table 5, DFP brain-blood = 0.67
    pc_liver <- fixed(1.53)
    label("Liver:blood partition coefficient (unitless)") # Table 5, DFP liver-blood = 1.53
    pc_kidney <- fixed(1.63)
    label("Kidney:blood partition coefficient (unitless)") # Table 5, DFP kidney-blood = 1.63
    pc_rich <- fixed(0.67)
    label("Richly perfused:blood partition coefficient (unitless)") # Table 5, DFP richly perfused blood = 0.67
    pc_fat <- fixed(17.6)
    label("Fat:blood partition coefficient (unitless)") # Table 5, DFP fat-blood = 17.6
    pc_slow <- fixed(0.77)
    label("Slowly perfused:blood partition coefficient (unitless)") # Table 5, DFP slowly perfused blood = 0.77
    pc_diaphragm <- fixed(0.77)
    label("Diaphragm:blood partition coefficient (unitless)") # Table 5, DFP diaphragm-blood = 0.77

    # =====================================================================
    # A-esterase (DFPase) hydrolysis, Michaelis-Menten -- Table 1, Rat.
    # Vmax in mg/h, Km in mg/L.
    # =====================================================================
    vmax_brain <- fixed(9.18)
    label("DFPase Vmax, brain (mg/h)") # Table 1, Brain Vmax, rat = 9.18
    km_brain <- fixed(439.80)
    label("DFPase Km, brain (mg/L)") # Table 1, Brain Km, rat = 439.80
    vmax_liver <- fixed(1380.00)
    label("DFPase Vmax, liver (mg/h)") # Table 1, Liver Vmax, rat = 1380.00
    km_liver <- fixed(237.36)
    label("DFPase Km, liver (mg/L)") # Table 1, Liver Km, rat = 237.36
    vmax_kidney <- fixed(103.32)
    label("DFPase Vmax, kidney (mg/h)") # Table 1, Kidney Vmax, rat = 103.32
    km_kidney <- fixed(134.3)
    label("DFPase Km, kidney (mg/L)") # Table 1, Kidney Km, rat = 134.3
    vmax_rich <- fixed(11.04)
    label("DFPase Vmax, richly perfused (mg/h)") # Table 1, Richly perfused Vmax, rat = 11.04
    km_rich <- fixed(50.89)
    label("DFPase Km, richly perfused (mg/L)") # Table 1, Richly perfused Km, rat = 50.89
    vmax_venous <- fixed(320.18)
    label("DFPase Vmax, venous blood (mg/h)") # Table 1, Venous blood Vmax, rat = 320.18
    km_venous <- fixed(198.70)
    label("DFPase Km, venous blood (mg/L)") # Table 1, Venous blood Km, rat = 198.70
    vmax_arterial <- fixed(106.73)
    label("DFPase Vmax, arterial blood (mg/h)") # Table 1, Arterial blood Vmax, rat = 106.73
    km_arterial <- fixed(198.70)
    label("DFPase Km, arterial blood (mg/L)") # Table 1, Arterial blood Km, rat = 198.70

    mw <- fixed(184.15)
    label("Molecular weight of DFP (g/mol)") # DFP C6H14FO3P; Appendix I uses 184 mg/mmole for the uM conversion

    # =====================================================================
    # AChE inhibition kinetics -- Table 3, Rat column.
    # =====================================================================
    ki <- fixed(14.16)
    label("AChE bimolecular inhibition rate constant (1/uM/h)") # Table 3, Enzyme inhibition, rat = 14.16 (= Table 1 AChE)
    kreg <- fixed(0.016)
    label("AChE regeneration rate constant, inhibited -> free (1/h)") # Table 3, Enzyme regeneration, rat = 0.016
    kage <- fixed(0.263)
    label("AChE aging rate constant, inhibited -> aged (1/h)") # Table 3, Enzyme aging, rat = 0.263
    kdeg_brain <- fixed(0.010)
    label("AChE first-order degradation rate constant, brain (1/h)") # Table 3, Enzyme degradation Brain, rat = 0.010
    kdeg_plasma <- fixed(0.100)
    label("AChE first-order degradation rate constant, plasma (1/h)") # Table 3, Enzyme degradation Plasma, rat = 0.100

    # =====================================================================
    # Extravascular (sc/im) absorption. NOT REPORTED numerically by Gearhart
    # 1994: the sc absorption rate was visually optimised (Methods: "the rate
    # of absorption of DFP into blood following a single sc injection ...
    # Simulation of the initial fall in plasma AChE activity was particularly
    # sensitive to this parameter"), but never tabulated. Placeholder ka =
    # 1/h; dose into depot for the repeated-sc protocol. See vignette Errata.
    # =====================================================================
    lka <- fixed(0.0)
    label("First-order sc/im absorption rate constant, log (log 1/h) -- NOT REPORTED") # Methods: sc absorption visually optimised, not tabulated; placeholder ka = 1/h
  })

  model({
    # -------------------------------------------------------------------
    # 1. Scaled volumes (L) and blood flows (L/h). Density 1 kg/L.
    # -------------------------------------------------------------------
    v_brain <- fv_brain * bw
    v_liver <- fv_liver * bw
    v_kidney <- fv_kidney * bw
    v_rich <- fv_rich * bw
    v_fat <- fv_fat * bw
    v_slow <- fv_slow * bw
    v_diaphragm <- fv_diaphragm * bw
    v_venous <- fv_venous * bw
    v_arterial <- fv_arterial * bw

    q_brain <- fq_brain * qc
    q_liver <- fq_liver * qc
    q_kidney <- fq_kidney * qc
    q_rich <- fq_rich * qc
    q_fat <- fq_fat * qc
    q_slow <- fq_slow * qc
    q_diaphragm <- fq_diaphragm * qc

    ka <- exp(lka)

    # -------------------------------------------------------------------
    # 2. Tissue concentrations (mg/L) and venous-effluent (capillary)
    #    concentrations cv = c / partition.
    # -------------------------------------------------------------------
    c_brain <- a_brain / v_brain
    cv_brain <- c_brain / pc_brain
    c_liver <- a_liver / v_liver
    cv_liver <- c_liver / pc_liver
    c_kidney <- a_kidney / v_kidney
    cv_kidney <- c_kidney / pc_kidney
    c_rich <- a_rich / v_rich
    cv_rich <- c_rich / pc_rich
    c_fat <- a_fat / v_fat
    cv_fat <- c_fat / pc_fat
    c_slow <- a_slow / v_slow
    cv_slow <- c_slow / pc_slow
    c_diaphragm <- a_diaphragm / v_diaphragm
    cv_diaphragm <- c_diaphragm / pc_diaphragm

    c_venous <- a_venous / v_venous
    c_arterial <- a_arterial / v_arterial

    # Lung is an algebraic steady-state gas exchanger: mixed venous blood
    # equilibrates with alveolar air across the blood:air partition and
    # leaves as the blood entering the arterial compartment. No inhalation
    # here (c_inh = 0), so this is a small pulmonary exhalation loss.
    c_inh <- 0
    c_lung_exit <- (qc * c_venous + qp * c_inh) / (qc + qp / pc_blood_air)

    # -------------------------------------------------------------------
    # 3. Michaelis-Menten DFPase hydrolysis (mg/h). Tissue terms act on the
    #    capillary concentration cv (Appendix II, CVBr); blood terms on the
    #    compartment concentration.
    # -------------------------------------------------------------------
    r_brain <- vmax_brain * cv_brain / (km_brain + cv_brain)
    r_liver <- vmax_liver * cv_liver / (km_liver + cv_liver)
    r_kidney <- vmax_kidney * cv_kidney / (km_kidney + cv_kidney)
    r_rich <- vmax_rich * cv_rich / (km_rich + cv_rich)
    r_venous <- vmax_venous * c_venous / (km_venous + c_venous)
    r_arterial <- vmax_arterial * c_arterial / (km_arterial + c_arterial)

    # -------------------------------------------------------------------
    # 4. DFP mass-balance ODE system (mg), Appendix II. The covalent
    #    DFP-esterase binding sink is negligible (Discussion) and omitted.
    # -------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(a_venous) <- (q_brain * cv_brain + q_liver * cv_liver +
      q_kidney * cv_kidney + q_rich * cv_rich + q_fat * cv_fat +
      q_slow * cv_slow + q_diaphragm * cv_diaphragm) -
      qc * c_venous - r_venous + ka * depot
    d/dt(a_arterial) <- qc * (c_lung_exit - c_arterial) - r_arterial
    d/dt(a_brain) <- q_brain * (c_arterial - cv_brain) - r_brain
    d/dt(a_liver) <- q_liver * (c_arterial - cv_liver) - r_liver
    d/dt(a_kidney) <- q_kidney * (c_arterial - cv_kidney) - r_kidney
    d/dt(a_rich) <- q_rich * (c_arterial - cv_rich) - r_rich
    d/dt(a_fat) <- q_fat * (c_arterial - cv_fat)
    d/dt(a_slow) <- q_slow * (c_arterial - cv_slow)
    d/dt(a_diaphragm) <- q_diaphragm * (c_arterial - cv_diaphragm)

    # Cumulative-amount accumulators for mass balance: total hydrolysed by
    # A-esterase, and total lost by pulmonary exhalation. At any time,
    # dose = (sum of DFP compartments) + a_hydrolyzed + a_exhaled.
    d/dt(a_hydrolyzed) <- r_brain + r_liver + r_kidney + r_rich + r_venous + r_arterial
    d/dt(a_exhaled) <- qc * (c_venous - c_lung_exit)

    # -------------------------------------------------------------------
    # 5. AChE free/inhibited/aged dynamics (fraction of baseline activity),
    #    Figure 2 / Appendix II. DFP driving concentrations converted to uM.
    #    At baseline (DFP = 0) synthesis balances degradation, so the free
    #    pool holds at 1: d(free)/dt = kdeg*(1 - free) encodes zero-order
    #    synthesis (kdeg * baseline) minus first-order degradation.
    # -------------------------------------------------------------------
    dfp_brain_uM <- c_brain * 1000 / mw
    dfp_plasma_uM <- c_venous * 1000 / mw

    d/dt(ache_free_brain) <- kdeg_brain * (1 - ache_free_brain) -
      ki * dfp_brain_uM * ache_free_brain + kreg * ache_inh_brain
    d/dt(ache_inh_brain) <- ki * dfp_brain_uM * ache_free_brain -
      kreg * ache_inh_brain - kage * ache_inh_brain
    d/dt(ache_aged_brain) <- kage * ache_inh_brain

    d/dt(ache_free_plasma) <- kdeg_plasma * (1 - ache_free_plasma) -
      ki * dfp_plasma_uM * ache_free_plasma + kreg * ache_inh_plasma
    d/dt(ache_inh_plasma) <- ki * dfp_plasma_uM * ache_free_plasma -
      kreg * ache_inh_plasma - kage * ache_inh_plasma
    d/dt(ache_aged_plasma) <- kage * ache_inh_plasma

    ache_free_brain(0) <- 1
    ache_free_plasma(0) <- 1

    # -------------------------------------------------------------------
    # 6. Outputs. Cc is venous plasma DFP; Cbrain is brain DFP. AChE
    #    activities are fractions of control (Figures 6-7).
    # -------------------------------------------------------------------
    Cc <- c_venous
    Cbrain <- c_brain
    ache_brain <- ache_free_brain
    ache_plasma <- ache_free_plasma
  })
}
