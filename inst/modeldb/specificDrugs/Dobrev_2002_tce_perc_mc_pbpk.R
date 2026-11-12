Dobrev_2002_tce_perc_mc_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited, five tissue groups per chemical;",
    "interactive ternary-mixture inhalation model recoded from the ACSL /",
    "Berkeley Madonna listing in the Appendix). Human toxicokinetics of a",
    "co-inhaled mixture of three volatile organic solvents: trichloroethylene",
    "(TCE), tetrachloroethylene / perchloroethylene (PERC) and 1,1,1-",
    "trichloroethane / methylchloroform (MC). Each solvent has its own set of",
    "four flow-limited well-mixed tissue compartments (fat, richly perfused,",
    "slowly perfused, liver); arterial blood is an algebraic steady-state lung",
    "equation and mixed venous blood is the flow-weighted average of the tissue",
    "effluents, so inhalation enters as a square-wave air concentration (ppm)",
    "rather than a dose record and no blood compartment is carried as a state.",
    "Hepatic metabolism is saturable (Michaelis-Menten) oxidation; the three",
    "single-chemical oxidations are coupled through competitive inhibition in",
    "the shared liver so that co-exposure to PERC and MC suppresses TCE",
    "oxidation and raises TCE blood levels. TCE additionally carries a",
    "first-order glutathione(GSH)-conjugation branch producing the nephrotoxic",
    "S-dichlorovinylcysteine (DCVC) metabolites. Deterministic typical-value",
    "model (no IIV, no residual error): the authors integrated previously",
    "validated single-chemical human PBPK models with rat-derived interaction",
    "constants and evaluated it against published human kinetic data rather",
    "than fitting it. Outputs are peak arterial TCE/PERC/MC blood concentration",
    "and the cumulative oxidative and DCVC metabolite amounts used as",
    "biomarkers of exposure.",
    sep = " "
  )
  reference <- paste(
    "Dobrev ID, Andersen ME, Yang RSH. In silico toxicology: simulating",
    "interaction thresholds for human exposure to mixtures of trichloroethylene,",
    "tetrachloroethylene, and 1,1,1-trichloroethane.",
    "Environ Health Perspect. 2002;110(10):1031-1039.",
    "doi:10.1289/ehp.021101031. PMCID: PMC1241030.",
    "Physiologic parameters Table 1; partition coefficients and biochemical",
    "constants Table 2; mass-balance differential equations in the Appendix",
    "(Equations 1-7). Single-chemical human models from Fisher et al. 1998 (TCE),",
    "Reitz et al. 1996 (PERC) and Reitz et al. 1988 (MC); interaction and GSH",
    "constants from Dobrev et al. 2001 and Clewell et al. 2000.",
    sep = " "
  )
  vignette <- "Dobrev_2002_tce_perc_mc"

  units <- list(
    time = "h",
    dosing = paste(
      "not a dose record: inhalation enters as the square-wave air-concentration",
      "parameters conc_air_ppm_tce / _perc / _mc (ppm), active while",
      "time < inhale_time"
    ),
    concentration = "mg/L",
    amount = "mg",
    weight = "kg"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "The paper performs every human simulation at a standard body weight of",
        "70 kg ('a standard body weight of 70 kg', Individual model calibration).",
        "Carried as a covariate so the model can be re-scaled: all tissue volumes",
        "are linear in body weight and cardiac output, alveolar ventilation and",
        "the metabolic Vmax are allometric at exponent 0.75 (Krishnan et al.",
        "1994a). The apparent affinity constants Km and the inhibition constants",
        "Ki are enzyme-specific and are NOT body-weight scaled.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  # Every ODE state is a flow-limited tissue AMOUNT (mg) for one of the three
  # solvents, plus cumulative metabolite accumulators. Arterial and mixed
  # venous blood are algebraic (Appendix Equations 6-7), so blood is not a
  # state. `a_metab_*` integrate the saturable oxidative rate (the paper's
  # 'total TCE metabolites' biomarker is a_metab_tce); `a_dcvc_tce` integrates
  # the first-order GSH/DCVC branch (TCE only).
  compartmentData <- list(
    a_fat_tce = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_rapid_tce = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_slow_tce = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_liver_tce = list(analyte = "Trichloroethylene", units = "mg", specimen = "tissue", verified = TRUE),
    a_fat_perc = list(
      analyte = "Tetrachloroethylene (perchloroethylene)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_rapid_perc = list(
      analyte = "Tetrachloroethylene (perchloroethylene)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_slow_perc = list(
      analyte = "Tetrachloroethylene (perchloroethylene)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_liver_perc = list(
      analyte = "Tetrachloroethylene (perchloroethylene)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_fat_mc = list(
      analyte = "1,1,1-trichloroethane (methylchloroform)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_rapid_mc = list(
      analyte = "1,1,1-trichloroethane (methylchloroform)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_slow_mc = list(
      analyte = "1,1,1-trichloroethane (methylchloroform)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_liver_mc = list(
      analyte = "1,1,1-trichloroethane (methylchloroform)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    a_metab_tce = list(
      analyte = "Trichloroethylene, cumulative oxidative (CYP450) metabolites",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_metab_perc = list(
      analyte = "Tetrachloroethylene, cumulative oxidative (CYP450) metabolites",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_metab_mc = list(
      analyte = "1,1,1-trichloroethane, cumulative oxidative (CYP450) metabolites",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_dcvc_tce = list(
      analyte = "Trichloroethylene, cumulative GSH-conjugation (S-DCVC) metabolites",
      units = "mg",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # All sixteen ODE states are paper-mechanistic PBPK compartments outside the
  # canonical register (per-chemical tissue amounts and metabolite integrators).
  paper_specific_compartments <- c(
    "a_fat_tce",
    "a_rapid_tce",
    "a_slow_tce",
    "a_liver_tce",
    "a_fat_perc",
    "a_rapid_perc",
    "a_slow_perc",
    "a_liver_perc",
    "a_fat_mc",
    "a_rapid_mc",
    "a_slow_mc",
    "a_liver_mc",
    "a_metab_tce",
    "a_metab_perc",
    "a_metab_mc",
    "a_dcvc_tce"
  )

  population <- list(
    species = "human",
    n_subjects = 0L,
    n_studies = 0L,
    age_range = "adult (single 70 kg reference adult; no age term in the model)",
    weight_range = "70 kg (standard body weight; Individual model calibration)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "healthy reference adult (occupational-exposure scenarios)",
    dose_range = paste(
      "Constant 6-hour occupational threshold-limit-value (TLV) inhalation",
      "exposures: TCE 50 ppm with PERC and MC co-exposure sweeps 0-2,000 ppm",
      "(Figures 6-8); interaction-threshold tables at TCE 50 ppm, MC 0-350 ppm,",
      "PERC 0-25 ppm (Table 3). No oral or dermal route.",
      sep = " "
    ),
    regions = "United States (occupational / environmental risk assessment)",
    notes = paste(
      "This is a forward in-silico risk-assessment simulation, not a fit, so",
      "n_subjects = 0. The human mixture model integrates previously developed",
      "and validated single-chemical human PBPK models (Fisher et al. 1998 TCE;",
      "Reitz et al. 1996 PERC; Reitz et al. 1988 MC) linked through competitive",
      "inhibition of hepatic oxidation. Interaction (Ki) and GSH (K) constants",
      "are rat-derived (Dobrev et al. 2001) and human DCVC parameterisation is",
      "from Clewell et al. 2000. Physiologic parameters are human reference",
      "values (Table 1). The paper also reports a rat model, but that uses a",
      "different three-compartment lumped structure from the upstream Dobrev et",
      "al. 2001 paper (not reproduced here); only the fully specified human",
      "five-compartment model is implemented.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Physiologic parameters -- Table 1 (human column). Cardiac output and
    # alveolar ventilation are reported per kg but are applied with
    # allometric body-weight scaling at exponent 0.75 (Krishnan et al.
    # 1994a); the per-kg figure is therefore the coefficient on BW^0.75.
    # =====================================================================
    qc_per_kg <- fixed(15.0)
    label("Cardiac output coefficient (L/h/kg^0.75)") # Table 1, 'Cardiac output (L/hr/kg)' 15.0
    qp_per_kg <- fixed(15.0)
    label("Alveolar ventilation coefficient (L/h/kg^0.75)") # Table 1, 'Alveolar ventilation (L/hr/kg)' 15.0
    e_wt_flow <- fixed(0.75)
    label("Allometric exponent on cardiac output, ventilation and Vmax (unitless)") # Approach and Methods: 'Allometric scaling was applied to rates of metabolism (Vmax, K), cardiac output (QC), and pulmonary ventilation (QP)' (Krishnan et al. 1994a)

    # ---- blood flows, fraction of cardiac output (Table 1, human) -------
    fq_liver <- fixed(0.25)
    label("Blood flow to liver, fraction of cardiac output (unitless)") # Table 1, 'Blood flow to liver' 0.25
    fq_rapid <- fixed(0.449)
    label("Blood flow to richly perfused tissue, fraction of cardiac output (unitless)") # Table 1, 'Blood flow to richly perfused tissue' 0.449
    fq_slow <- fixed(0.249)
    label("Blood flow to slowly perfused tissue, fraction of cardiac output (unitless)") # Table 1, 'Blood flow to slowly perfused tissue' 0.249
    fq_fat <- fixed(0.052)
    label("Blood flow to fat, fraction of cardiac output (unitless)") # Table 1, 'Blood flow to fat' 0.052

    # ---- tissue volumes, fraction of body weight (Table 1, human) -------
    fv_slow <- fixed(0.651)
    label("Slowly perfused tissue volume, fraction of body weight (unitless)") # Table 1, 'Fraction slowly perfused' 0.651
    fv_rapid <- fixed(0.118)
    label("Richly perfused tissue volume, fraction of body weight (unitless)") # Table 1, 'Fraction richly perfused' 0.118
    fv_fat <- fixed(0.214)
    label("Fat volume, fraction of body weight (unitless)") # Table 1, 'Fraction fat' 0.214
    fv_liver <- fixed(0.026)
    label("Liver volume, fraction of body weight (unitless)") # Table 1, 'Fraction liver' 0.026

    # =====================================================================
    # Partition coefficients -- Table 2 (human columns).
    # =====================================================================
    # ---- trichloroethylene (TCE), Fisher et al. 1998 --------------------
    pc_blood_air_tce <- fixed(11.2)
    label("TCE blood/air partition coefficient (unitless)") # Table 2, TCE human Blood/air 11.2
    pc_fat_tce <- fixed(52.3)
    label("TCE fat/blood partition coefficient (unitless)") # Table 2, TCE human Fat/blood 52.3
    pc_liver_tce <- fixed(4.9)
    label("TCE liver/blood partition coefficient (unitless)") # Table 2, TCE human Liver/blood 4.9
    pc_rapid_tce <- fixed(1.1)
    label("TCE richly perfused tissue/blood partition coefficient (unitless)") # Table 2, TCE human Richly perfused tissue/blood 1.1
    pc_slow_tce <- fixed(1.4)
    label("TCE slowly perfused tissue/blood partition coefficient (unitless)") # Table 2, TCE human Slowly perfused tissue/blood 1.4

    # ---- tetrachloroethylene (PERC), Reitz et al. 1996 ------------------
    pc_blood_air_perc <- fixed(10.3)
    label("PERC blood/air partition coefficient (unitless)") # Table 2, PERC human Blood/air 10.3
    pc_fat_perc <- fixed(119.1)
    label("PERC fat/blood partition coefficient (unitless)") # Table 2, PERC human Fat/blood 119.1
    pc_liver_perc <- fixed(5.9)
    label("PERC liver/blood partition coefficient (unitless)") # Table 2, PERC human Liver/blood 5.9
    pc_rapid_perc <- fixed(5.9)
    label("PERC richly perfused tissue/blood partition coefficient (unitless)") # Table 2, PERC human Richly perfused tissue/blood 5.9
    pc_slow_perc <- fixed(3.1)
    label("PERC slowly perfused tissue/blood partition coefficient (unitless)") # Table 2, PERC human Slowly perfused tissue/blood 3.1

    # ---- 1,1,1-trichloroethane (MC), Reitz et al. 1988 -----------------
    pc_blood_air_mc <- fixed(2.5)
    label("MC blood/air partition coefficient (unitless)") # Table 2, MC human Blood/air 2.5
    pc_fat_mc <- fixed(104.0)
    label("MC fat/blood partition coefficient (unitless)") # Table 2, MC human Fat/blood 104.0
    pc_liver_mc <- fixed(3.4)
    label("MC liver/blood partition coefficient (unitless)") # Table 2, MC human Liver/blood 3.4
    pc_rapid_mc <- fixed(3.4)
    label("MC richly perfused tissue/blood partition coefficient (unitless)") # Table 2, MC human Richly perfused tissue/blood 3.4
    pc_slow_mc <- fixed(0.6)
    label("MC slowly perfused tissue/blood partition coefficient (unitless)") # Table 2, MC human Slowly perfused tissue/blood 0.6

    # =====================================================================
    # Biochemical constants -- Table 2 (human columns). Vmax is per kg and
    # is allometrically scaled (BW^0.75); Km and Ki are enzyme-specific and
    # unscaled. Ki is the competitive-inhibition constant of each chemical
    # acting as an inhibitor of the shared hepatic oxidation.
    # =====================================================================
    vmax_per_kg_tce <- fixed(4.0)
    label("TCE maximum oxidative metabolic rate coefficient (mg/h/kg^0.75)") # Table 2, TCE human Vmax 4.0 mg/hr/kg
    km_tce <- fixed(1.80)
    label("TCE Michaelis-Menten affinity constant (mg/L)") # Table 2, TCE human Km 1.80
    vmax_per_kg_perc <- fixed(1.42)
    label("PERC maximum oxidative metabolic rate coefficient (mg/h/kg^0.75)") # Table 2, PERC human Vmax 1.42 mg/hr/kg
    km_perc <- fixed(4.66)
    label("PERC Michaelis-Menten affinity constant (mg/L)") # Table 2, PERC human Km 4.66
    vmax_per_kg_mc <- fixed(0.42)
    label("MC maximum oxidative metabolic rate coefficient (mg/h/kg^0.75)") # Table 2, MC human Vmax 0.42 mg/hr/kg
    km_mc <- fixed(5.75)
    label("MC Michaelis-Menten affinity constant (mg/L)") # Table 2, MC human Km 5.75

    ki_tce <- fixed(100000)
    label("TCE inhibition constant as an inhibitor of PERC/MC oxidation (mg/L)") # Methods: 'only inhibition of TCE metabolism was implemented ... (i.e., Ki = 100,000 for PERC and MC)' -- TCE's inhibitory effect on the others is disabled
    ki_perc <- fixed(6.3)
    label("PERC inhibition constant as an inhibitor of TCE oxidation (mg/L)") # Table 2, PERC Ki 6.3 (Dobrev et al. 2001 rat estimate)
    ki_mc <- fixed(4.2)
    label("MC inhibition constant as an inhibitor of TCE oxidation (mg/L)") # Table 2, MC Ki 4.2 (Dobrev et al. 2001 rat estimate)

    kgsh_tce <- fixed(0.015)
    label("TCE first-order GSH/DCVC conjugation rate constant (1/h)") # Table 2, TCE K 0.015 (Clewell et al. 2000)

    # =====================================================================
    # Molecular weights and molar gas volume for the ppm <-> mg/L
    # conversion. Not printed in the paper; standard chemical constants and
    # the 24450 mL/mol ideal-gas molar volume (25 C, 1 atm) used by the
    # ACSL / Berkeley Madonna toolkit (see vignette Errata).
    # =====================================================================
    mw_tce <- fixed(131.39)
    label("Molecular weight of trichloroethylene, C2HCl3 (g/mol)") # standard constant; not printed in paper (see Errata)
    mw_perc <- fixed(165.83)
    label("Molecular weight of tetrachloroethylene, C2Cl4 (g/mol)") # standard constant; not printed in paper (see Errata)
    mw_mc <- fixed(133.40)
    label("Molecular weight of 1,1,1-trichloroethane, C2H3Cl3 (g/mol)") # standard constant; not printed in paper (see Errata)
    molar_vol <- fixed(24450)
    label("Ideal-gas molar volume for the ppm <-> mg/L conversion (mL/mol)") # standard ACSL/Berkeley Madonna constant (25 C, 1 atm); not printed in paper (see Errata)

    # =====================================================================
    # Exposure parameters. The published simulations are constant 6-hour
    # inhalations; conc_air_ppm_* set the three air concentrations (ppm) and
    # inhale_time the window length. Defaults give a single 6-hour TCE-only
    # 50 ppm exposure; override per scenario through rxSolve(params = ...).
    # =====================================================================
    conc_air_ppm_tce <- fixed(50)
    label("Inhaled TCE air concentration during the exposure window (ppm)") # Results: 'constant 6-hr occupational TLV exposure to TCE (50 ppm)'
    conc_air_ppm_perc <- fixed(0)
    label("Inhaled PERC air concentration during the exposure window (ppm)") # co-exposure sweep variable (Figures 6-8; Table 3)
    conc_air_ppm_mc <- fixed(0)
    label("Inhaled MC air concentration during the exposure window (ppm)") # co-exposure sweep variable (Figures 6-8; Table 3)
    inhale_time <- fixed(6)
    label("Duration of the inhalation exposure window (h)") # Results: 'constant 6-hr ... exposure'
    inhale_interval <- fixed(100000)
    label("Period between repeated inhalation exposures (h); a large value gives a single window") # single-window default (not a repeated protocol in this paper)

    # =====================================================================
    # No IIV and no residual error: the authors integrated validated
    # deterministic single-chemical PBPK models and evaluated the mixture
    # model against published human kinetic data 'without further
    # statistical optimization' (Individual model calibration).
    # =====================================================================
  })

  model({
    # -------------------------------------------------------------------
    # 1. Scaled volumes (L), blood flows (L/h) and metabolic capacities.
    # -------------------------------------------------------------------
    v_fat <- fv_fat * WT
    v_rapid <- fv_rapid * WT
    v_slow <- fv_slow * WT
    v_liver <- fv_liver * WT

    qc <- qc_per_kg * WT^e_wt_flow
    qp <- qp_per_kg * WT^e_wt_flow
    q_liver <- fq_liver * qc
    q_rapid <- fq_rapid * qc
    q_slow <- fq_slow * qc
    q_fat <- fq_fat * qc

    vmax_tce <- vmax_per_kg_tce * WT^e_wt_flow
    vmax_perc <- vmax_per_kg_perc * WT^e_wt_flow
    vmax_mc <- vmax_per_kg_mc * WT^e_wt_flow

    # -------------------------------------------------------------------
    # 2. Inhalation square wave (idiom from Boone_2025_vinylchloride_pbpk).
    #    MOD(time, interval) is rendered as time - interval*floor(.); the
    #    Berkeley Madonna IF/THEN/ELSE air switch becomes the complement of
    #    the comparison, which rxode2 evaluates to 1 or 0. ppm -> mg/L uses
    #    Cair = ppm * MW / molar_vol (Appendix; standard gas-law conversion).
    # -------------------------------------------------------------------
    time_mod <- time - inhale_interval * floor(time / inhale_interval)
    air_on <- 1 - (time_mod >= inhale_time)
    cair_tce <- conc_air_ppm_tce * mw_tce / molar_vol * air_on
    cair_perc <- conc_air_ppm_perc * mw_perc / molar_vol * air_on
    cair_mc <- conc_air_ppm_mc * mw_mc / molar_vol * air_on

    # -------------------------------------------------------------------
    # 3. Tissue concentrations (mg/L) and venous-effluent concentrations
    #    (tissue concentration / partition coefficient, Appendix Eq 1).
    # -------------------------------------------------------------------
    c_liver_tce <- a_liver_tce / v_liver
    cv_liver_tce <- c_liver_tce / pc_liver_tce
    cv_rapid_tce <- (a_rapid_tce / v_rapid) / pc_rapid_tce
    cv_slow_tce <- (a_slow_tce / v_slow) / pc_slow_tce
    cv_fat_tce <- (a_fat_tce / v_fat) / pc_fat_tce

    c_liver_perc <- a_liver_perc / v_liver
    cv_liver_perc <- c_liver_perc / pc_liver_perc
    cv_rapid_perc <- (a_rapid_perc / v_rapid) / pc_rapid_perc
    cv_slow_perc <- (a_slow_perc / v_slow) / pc_slow_perc
    cv_fat_perc <- (a_fat_perc / v_fat) / pc_fat_perc

    c_liver_mc <- a_liver_mc / v_liver
    cv_liver_mc <- c_liver_mc / pc_liver_mc
    cv_rapid_mc <- (a_rapid_mc / v_rapid) / pc_rapid_mc
    cv_slow_mc <- (a_slow_mc / v_slow) / pc_slow_mc
    cv_fat_mc <- (a_fat_mc / v_fat) / pc_fat_mc

    # -------------------------------------------------------------------
    # 4. Mixed venous (Appendix Eq 7, flow-weighted tissue effluents) and
    #    arterial (Appendix Eq 6, steady-state lung) blood. Both algebraic;
    #    blood is not a compartment.
    # -------------------------------------------------------------------
    cv_tce <- (q_liver * cv_liver_tce + q_rapid * cv_rapid_tce +
      q_slow * cv_slow_tce + q_fat * cv_fat_tce) / qc
    ca_tce <- (qc * cv_tce + qp * cair_tce) / (qc + qp / pc_blood_air_tce)

    cv_perc <- (q_liver * cv_liver_perc + q_rapid * cv_rapid_perc +
      q_slow * cv_slow_perc + q_fat * cv_fat_perc) / qc
    ca_perc <- (qc * cv_perc + qp * cair_perc) / (qc + qp / pc_blood_air_perc)

    cv_mc <- (q_liver * cv_liver_mc + q_rapid * cv_rapid_mc +
      q_slow * cv_slow_mc + q_fat * cv_fat_mc) / qc
    ca_mc <- (qc * cv_mc + qp * cair_mc) / (qc + qp / pc_blood_air_mc)

    # -------------------------------------------------------------------
    # 5. Saturable hepatic oxidation with competitive inhibition
    #    (Appendix Eq 3). The apparent Km of each chemical is inflated by
    #    the other two inhibitors' effluent liver concentrations divided by
    #    their Ki. TCE is inhibited by PERC and MC; the paper disables TCE's
    #    inhibition of PERC and MC with ki_tce = 100,000.
    # -------------------------------------------------------------------
    ram_tce <- vmax_tce * cv_liver_tce /
      (km_tce * (1 + cv_liver_perc / ki_perc + cv_liver_mc / ki_mc) + cv_liver_tce)
    ram_perc <- vmax_perc * cv_liver_perc /
      (km_perc * (1 + cv_liver_tce / ki_tce + cv_liver_mc / ki_mc) + cv_liver_perc)
    ram_mc <- vmax_mc * cv_liver_mc /
      (km_mc * (1 + cv_liver_tce / ki_tce + cv_liver_perc / ki_perc) + cv_liver_mc)

    # First-order GSH/DCVC branch for TCE (Appendix Eq 4): production =
    # K * liver volume * TCE liver concentration = kgsh_tce * a_liver_tce.
    rgsh_tce <- kgsh_tce * v_liver * c_liver_tce

    # -------------------------------------------------------------------
    # 6. Mass-balance ODEs (Appendix Eqs 1-5). Liver carries metabolism;
    #    the other tissues are pure flow-limited uptake.
    # -------------------------------------------------------------------
    d/dt(a_fat_tce) <- q_fat * (ca_tce - cv_fat_tce)
    d/dt(a_rapid_tce) <- q_rapid * (ca_tce - cv_rapid_tce)
    d/dt(a_slow_tce) <- q_slow * (ca_tce - cv_slow_tce)
    d/dt(a_liver_tce) <- q_liver * (ca_tce - cv_liver_tce) - ram_tce - rgsh_tce

    d/dt(a_fat_perc) <- q_fat * (ca_perc - cv_fat_perc)
    d/dt(a_rapid_perc) <- q_rapid * (ca_perc - cv_rapid_perc)
    d/dt(a_slow_perc) <- q_slow * (ca_perc - cv_slow_perc)
    d/dt(a_liver_perc) <- q_liver * (ca_perc - cv_liver_perc) - ram_perc

    d/dt(a_fat_mc) <- q_fat * (ca_mc - cv_fat_mc)
    d/dt(a_rapid_mc) <- q_rapid * (ca_mc - cv_rapid_mc)
    d/dt(a_slow_mc) <- q_slow * (ca_mc - cv_slow_mc)
    d/dt(a_liver_mc) <- q_liver * (ca_mc - cv_liver_mc) - ram_mc

    # Cumulative metabolite accumulators (biomarkers of exposure).
    d/dt(a_metab_tce) <- ram_tce
    d/dt(a_metab_perc) <- ram_perc
    d/dt(a_metab_mc) <- ram_mc
    d/dt(a_dcvc_tce) <- rgsh_tce

    # -------------------------------------------------------------------
    # 7. Outputs. Cc is arterial TCE blood concentration -- the 'peak TCE
    #    blood level' biomarker the paper reports (end of the 6-hour
    #    exposure). Arterial PERC/MC and the liver concentrations are also
    #    exposed for inspection.
    # -------------------------------------------------------------------
    Cc <- ca_tce
    Cc_perc <- ca_perc
    Cc_mc <- ca_mc
    Cliver_tce <- c_liver_tce
    Cliver_perc <- c_liver_perc
    Cliver_mc <- c_liver_mc
  })
}
