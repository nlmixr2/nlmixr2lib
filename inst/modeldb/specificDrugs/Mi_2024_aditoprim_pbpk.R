Mi_2024_aditoprim_pbpk <- function() {
  description <- paste(
    "Veterinary (pig). PBPK (whole-body, seven-compartment, Berkeley",
    "Madonna 8.3.23.0) for the diaminopyrimidine aditoprim (ADP) after",
    "intramuscular injection in swine, developed to optimise the dosage",
    "regimen against Streptococcus suis and to predict the withdrawal",
    "interval for edible tissues (Mi et al. 2024, Front Pharmacol",
    "15:1378034). A single well-mixed blood pool, liver, kidney, muscle,",
    "fat and a lumped rest-of-body compartment are perfusion-limited and",
    "well stirred; liver, kidney, muscle and fat are carried as separate",
    "states because they are the edible tissues against which the maximum",
    "residue limits are set. Intramuscular absorption is a two-step",
    "process: Frac = 92% of the injected dose is immediately available at",
    "the injection site and is absorbed into the venous mixing pool at",
    "first-order rate Kim, while the remaining 8% is released from a slow",
    "depot at first-order rate Kdiss into that same injection-site pool.",
    "Elimination is renal (from kidney) plus hepatic metabolism (from",
    "liver). Plasma protein binding rescales the free plasma concentration",
    "only -- it does not gate tissue perfusion, which is driven by TOTAL",
    "arterial concentration in the published code. The seven parameters",
    "that Mi 2024 Table 3 carried through the 1000-animal Monte Carlo",
    "pop-PBPK analysis carry between-animal variability here; every other",
    "parameter is a fixed physiological or fitted chemical-specific",
    "constant. The model was hand-calibrated in Berkeley Madonna against",
    "digitised literature data, so no standard errors and no residual-error",
    "model are reported -- the propSd term is a placeholder. Several",
    "Table 2 / supplement-code discrepancies are reproduced as coded; see",
    "the vignette Errata."
  )
  reference <- paste(
    "Mi K, Sun L, Zhang L, Tang A, Tian X, Hou Y, Sun L, Huang L.",
    "A physiologically based pharmacokinetic/pharmacodynamic model to",
    "determine dosage regimens and withdrawal intervals of aditoprim",
    "against Streptococcus suis. Front Pharmacol. 2024;15:1378034.",
    "doi:10.3389/fphar.2024.1378034.",
    "Model equations transcribed from Supplementary Material section 3.1",
    "(Berkeley Madonna code) and Supplementary Table 3; final parameter",
    "values from Table 2; Monte Carlo distributions from Table 3.",
    sep = " "
  )
  vignette <- "Mi_2024_aditoprim"

  # `metabolized` is the cumulative amount cleared by hepatic metabolism
  # (`AML` in the Supplementary Material code, fed by `RML = KML*CVL*BW`).
  # It is a terminal bookkeeping sink with no feedback into the dynamics;
  # it exists so the published mass-balance identity (Supplementary
  # Material 3.1, `; Mass balance`) can be checked. There is no canonical
  # compartment for a metabolic elimination sink in
  # inst/references/compartment-names.md, so it is declared here rather
  # than silently extending the register. The sibling model
  # Mi_2023_cefquinome_pbpk.R declares its biliary sink `bile` the same
  # way.
  paper_specific_compartments <- c("metabolized")

  # Time in hours, doses in mg (the paper prescribes mg/kg; the vignette
  # multiplies by WT), concentrations in ug/mL. States hold amounts in mg
  # and volumes are in L, so amount/volume is mg/L, which is numerically
  # identical to the ug/mL the paper reports for plasma and to the ug/g
  # used for the maximum residue limits in liver and kidney.
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mi 2024 Table 2 fixes BW = 30 kg, and the Supplementary Material",
        "code sets `BW=30`. Every organ volume and blood flow is a fixed",
        "fraction of BW, and the renal and hepatic clearances are per-kg",
        "constants, so WT is carried here as a covariate exactly as BW",
        "enters the published code; the model is therefore linear in dose",
        "and scales isometrically with WT. Note that Mi 2024 section 2.1",
        "says the physiological parameters were taken from Upton (2008)",
        "for 25 kg pigs while Table 2 and the code both use BW = 30 kg;",
        "see vignette Errata. The underlying residue-depletion studies",
        "used market-age swine."
      ),
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "aditoprim", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "aditoprim", units = "mg", specimen = "administration site", verified = TRUE),
    blood       = list(analyte = "aditoprim", units = "mg", specimen = "whole blood", verified = TRUE),
    liver       = list(analyte = "aditoprim", units = "mg", specimen = "tissue", verified = TRUE),
    kidney      = list(analyte = "aditoprim", units = "mg", specimen = "tissue", verified = TRUE),
    muscle      = list(analyte = "aditoprim", units = "mg", specimen = "tissue", verified = TRUE),
    adipose     = list(analyte = "aditoprim", units = "mg", specimen = "tissue", verified = TRUE),
    other       = list(analyte = "aditoprim", units = "mg", specimen = "tissue", verified = TRUE),
    urine       = list(analyte = "aditoprim", units = "mg", specimen = "urine", verified = TRUE),
    metabolized = list(analyte = "aditoprim", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "pig (swine)",
    n_subjects     = 92L,
    n_studies      = 4L,
    age_range      = "grower / finisher swine",
    weight_range   = "model reference BW = 30 kg; Upton (2008) physiology is quoted for 25 kg pigs",
    disease_state  = "healthy",
    dose_range     = "5 and 10 mg/kg intramuscular aditoprim injection in the calibration and validation datasets; 0, 2.5, 5, 10, 12.5, 15, 20 and 25 mg/kg once or twice daily were simulated",
    regions        = "China",
    notes          = paste(
      "No animals were dosed for this modelling paper. The structural",
      "model was calibrated and validated against four previously",
      "published swine datasets (Mi 2024 Table 1), with graphical data",
      "digitised using WebPlotDigitizer 3.10. Calibration: Wang et al.",
      "(2022) plasma (n = 6, single 5 mg/kg IM) and Wang (2020) liver,",
      "kidney, muscle and fat (n = 40, 10 mg/kg every 12 h for 14 doses).",
      "Validation: Qu et al. (2022) plasma (n = 6, single 5 mg/kg IM) and",
      "Wang (2016) liver, kidney, muscle and fat (n = 40, 5 mg/kg every",
      "24 h for 7 doses). The tissue concentrations behind the two",
      "residue-depletion datasets are reproduced in Supplementary Tables",
      "1 and 2. All data were from healthy animals. The pop-PBPK layer",
      "simulates 1000 virtual animals by Monte Carlo sampling of the",
      "seven influential parameters in Table 3; no individual-level",
      "fitting was performed and no between-animal variance was estimated",
      "from data -- the coefficients of variation are literature defaults",
      "as described in Mi 2024 section 2.5."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Intramuscular absorption (Mi 2024 Table 2, "Absorption rate
    # constant"; all three "Model fitting"). Note that Frac is the FAST
    # fraction here: Supplementary Material 3.1 has
    # `Rpenim = Rinputim*(Frac)` feeding the injection-site pool and
    # `Rppgim = Rinputim*(1-Frac)` feeding the slow depot, and Mi 2024
    # section 2.3 says "Drugs were distributed into fast absorption
    # (Frac*doseim) and slow absorption [(1-Frac)*doseim]". This is the
    # OPPOSITE sense to the sibling Mi_2023_cefquinome_pbpk.R, where
    # Frac = 0.1 is the slow fraction.
    # ---------------------------------------------------------------
    lka <- log(1.3)
    label("Absorption rate constant Kim, injection site to venous mixing pool (1/h)")   # Mi 2024 Table 2 (Kim = 1.3 /h, Model fitting); Suppl 3.1 `Kim = 1.3`

    lkdiss <- log(0.0118)
    label("Dissolution rate constant Kdiss, slow depot to injection site (1/h)")        # Mi 2024 Table 2 (Kdiss = 0.0118 /h, Model fitting); Suppl 3.1 has the pre-optimisation `Kdiss = 0.0115`

    lfrac <- log(0.92)
    label("Fraction of the intramuscular dose entering the fast injection-site pool (unitless)") # Mi 2024 Table 2 (Frac = 0.92, Model fitting); Suppl 3.1 has the pre-optimisation `Frac = 0.89`

    # ---------------------------------------------------------------
    # Tissue-to-blood partition coefficients (Mi 2024 Table 2, all
    # "Calculated/optimized"). Initial values were computed by the AUC
    # method as AUCtissue/AUCplasma (Mi 2024 section 2.3, after Wang
    # et al. 2016) and then optimised against the residual data. The
    # Supplementary Material code carries the pre-optimisation values
    # (PL 4, PK 5, PM 0.75, PF 1, POT 0.26); Table 2 is the final record
    # and is what is encoded here. See vignette Errata.
    # ---------------------------------------------------------------
    lkp_liver <- log(5.249)
    label("Liver-to-blood partition coefficient PL (unitless)")                         # Mi 2024 Table 2 (PL = 5.249, Calculated/optimized)

    lkp_kidney <- log(6)
    label("Kidney-to-blood partition coefficient PK (unitless)")                        # Mi 2024 Table 2 (PK = 6, Calculated/optimized)

    lkp_muscle <- log(0.79)
    label("Muscle-to-blood partition coefficient PM (unitless)")                        # Mi 2024 Table 2 (PM = 0.79, Calculated/optimized)

    lkp_adipose <- log(1.1)
    label("Fat-to-blood partition coefficient PF (unitless)")                           # Mi 2024 Table 2 (PF = 1.1, Calculated/optimized)

    lkp_other <- log(0.18)
    label("Rest-of-body-to-blood partition coefficient PT (unitless)")                  # Mi 2024 Table 2 (PT = 0.18, Calculated/optimized)

    # ---------------------------------------------------------------
    # Plasma protein binding (Mi 2024 Table 2, "Model fitting"). Table 2
    # reports PB = 0.82 as the BOUND fraction ("Percentage of plasma
    # protein binding") and the code uses it as `CAfree = CA*(1-PB)`, so
    # the unbound fraction is 1 - 0.82 = 0.18. Carried on the bound scale
    # because that is the scale Mi 2024 Table 3 gives the Monte Carlo
    # lognormal on; `fu` is derived in model(). PB does NOT gate tissue
    # perfusion in the published code -- the tissue ODEs are driven by
    # total arterial concentration -- so it rescales the free plasma
    # output (and hence the PD driver and fAUC/MIC) only.
    # ---------------------------------------------------------------
    lfb <- log(0.82)
    label("Fraction of aditoprim bound to plasma protein PB (unitless)")                # Mi 2024 Table 2 (PB = 0.82, Model fitting); Suppl 3.1 has the pre-optimisation `PB=0.87`, `CAfree = CA*(1-PB)`

    # ---------------------------------------------------------------
    # Elimination (Mi 2024 Table 2, both "Model fitting"). Both are per-kg
    # clearances and are multiplied by WT in model().
    # ---------------------------------------------------------------
    lcl_renal <- log(0.1)
    label("Renal clearance KurineC (L/h/kg), applied to kidney venous concentration")   # Mi 2024 Table 2 (KurineC = 0.1 L/h/kg); Suppl 3.1 has the pre-optimisation `KurineC =0.11`

    lcl_nonren <- log(0.01)
    label("Hepatic metabolic clearance KML (L/h/kg), applied to liver venous concentration") # Mi 2024 Table 2 (KML = 0.01 L/h/kg); Suppl 3.1 `KML=0.01`

    # ---------------------------------------------------------------
    # Cardiac output. Fixed swine physiology (Upton 2008), but kept in
    # ini() so it is visible alongside the other scaling constants. Every
    # other physiological constant is written directly in model().
    # ---------------------------------------------------------------
    qcc <- fixed(4.944)
    label("Cardiac output QCC (L/h/kg)")                                                # Mi 2024 Table 2 (QCC = 4.944 L/h/kg, Upton 2008); Suppl 3.1 `QCAR=4.944`

    # ---------------------------------------------------------------
    # Between-animal variability, from the Monte Carlo distributions of
    # Mi 2024 Table 3. These are NOT estimated variances -- they are the
    # default coefficients of variation the authors assumed (Mi 2024
    # section 2.5), so all are fixed().
    #
    # Every Table 3 row is lognormal and is parameterised so that its
    # ARITHMETIC MEAN equals the Table 2 point estimate. The log-scale
    # variance is therefore log(1 + CV^2), and all seven pairs of printed
    # 2.5 / 97.5 percentile bounds are reproduced exactly by that
    # parameterisation (worked arithmetic in the vignette source-trace
    # table). As in the sibling Mi_2023_cefquinome_pbpk.R, the Table 2
    # value is used here as the typical value (i.e. the MEDIAN) so that a
    # typical-value simulation reproduces the deterministic Table 2 model
    # and Figures 2, 3 and 6 exactly; that places the median about 1-7%
    # below the Table 3 arithmetic mean. See vignette Errata.
    #
    # Mi 2024 section 2.5 says the default CVs were 30% for
    # chemical-specific and 40% for physiological parameters, but Table 3
    # prints 40% for the six chemical-specific rows and 10% for Frac, and
    # the Table 3 bounds confirm the printed CVs. Table 3 is followed.
    # ---------------------------------------------------------------
    etalkp_liver ~ fixed(0.1484200)
    # Mi 2024 Table 3 (PL, Lognormal, CV 0.40): log(1 + 0.40^2) = 0.1484200
    etalkp_kidney ~ fixed(0.1484200)
    # Mi 2024 Table 3 (PK, Lognormal, CV 0.40)
    etalkp_muscle ~ fixed(0.1484200)
    # Mi 2024 Table 3 (PM, Lognormal, CV 0.40)
    etalkp_adipose ~ fixed(0.1484200)
    # Mi 2024 Table 3 (PF, Lognormal, CV 0.40)
    etalcl_renal ~ fixed(0.1484200)
    # Mi 2024 Table 3 (KurineC, Lognormal, CV 0.40)
    etalfb ~ fixed(0.1484200)
    # Mi 2024 Table 3 (PB, Lognormal, CV 0.40)
    etalfrac ~ fixed(0.00995033)
    # Mi 2024 Table 3 (Frac, Lognormal, CV 0.10): log(1 + 0.10^2) = 0.00995033

    # ---------------------------------------------------------------
    # Mi 2024 calibrated the model by hand in Berkeley Madonna 8.3.23.0
    # against digitised literature data and reports no residual-error
    # model, no standard errors and no objective function. nlmixr2 model
    # definitions require a residual-error term, so propSd below is a
    # fixed placeholder for syntactic completeness only and must NOT be
    # read as an estimate. Same convention as
    # Mi_2023_cefquinome_pbpk and Kang_2023_artesunate_hamster_pbpk.
    # ---------------------------------------------------------------
    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                         # not reported in Mi 2024; placeholder only
  })

  model({
    # =================================================================
    # Swine physiology (Mi 2024 Table 2, Upton 2008; Supplementary
    # Material 3.1 {Physiological parameters}). Blood flows are fractions
    # of cardiac output; volumes are fractions of body weight.
    # =================================================================
    q_co <- qcc * WT                                  # L/h; Suppl 3.1 `QC=QCAR*BW`

    # Fractional flows. QFC is reported ONLY in the Supplementary
    # Material code -- Mi 2024 Table 2 omits both the fat blood flow and
    # the fat volume even though fat is a modelled, calibrated and
    # validated output tissue (Figures 2D, 3) and PF is one of the seven
    # Monte Carlo parameters. The rest-of-body fraction is computed as
    # the complement exactly as the code does. Note this evaluates to
    # 0.1278, not the 0.3055 printed in Mi 2024 Table 2; the printed
    # value is the complement taken WITHOUT the fat compartment and is
    # incompatible with flow balance once fat is present. See vignette
    # Errata.
    qlc_frac <- 0.3053                                # Mi 2024 Table 2 (QLC = 0.3053)
    qkc_frac <- 0.1398                                # Mi 2024 Table 2 (QKC = 0.1398)
    qmc_frac <- 0.2524                                # Mi 2024 Table 2 (QMC = 0.2524)
    qfc_frac <- 0.1747                                # Suppl 3.1 `QFC=0.1747`; absent from Table 2
    qoc_frac <- 1 - qlc_frac - qkc_frac - qmc_frac - qfc_frac
    # Suppl 3.1 `QOT=QC-(QL+QK+QM+QF)`

    q_liver <- qlc_frac * q_co
    q_kidney <- qkc_frac * q_co
    q_muscle <- qmc_frac * q_co
    q_adipose <- qfc_frac * q_co
    q_other <- qoc_frac * q_co

    # Organ volumes (L). VFC is likewise code-only. The blood pool is a
    # SINGLE well-mixed compartment of VBC = 0.06 of body weight, which
    # is exactly the sum of the arterial (VartC 0.016) and venous
    # (VvenC 0.044) fractions printed in Mi 2024 Table 2. The rest-of-body
    # volume is the complement, and evaluates to 0.2066 rather than the
    # 0.232 printed in Table 2. Table 2 also lists a lung (VLUC 0.01,
    # QLUC 1) that appears in the Figure 1A schematic but has no state
    # and no equation in either Supplementary Table 3 or the code; it is
    # therefore absent here. See vignette Errata.
    v_liver <- 0.0294 * WT                            # Mi 2024 Table 2 (VLC = 0.0294)
    v_kidney <- 0.004 * WT                            # Mi 2024 Table 2 (VKC = 0.004)
    v_muscle <- 0.4 * WT                              # Mi 2024 Table 2 (VMC = 0.4)
    v_adipose <- 0.3 * WT                             # Suppl 3.1 `VFC=0.3`; absent from Table 2
    v_blood <- 0.06 * WT                              # Suppl 3.1 `VBC=0.06` = Table 2 VartC 0.016 + VvenC 0.044
    v_other <- (1 - 0.0294 - 0.004 - 0.4 - 0.3 - 0.06) * WT
    # Suppl 3.1 `VOT=BW-(VL+VK+VM+VF+Vblood)`

    # =================================================================
    # Individual chemical-specific parameters
    # =================================================================
    ka <- exp(lka)
    kdiss <- exp(lkdiss)

    # Mi 2024 section 2.5: the Monte Carlo draws were confined to the
    # 2.5 / 97.5 percentile bounds printed in Table 3 ("The distributions
    # within the lower bound (2.5%) and upper bound (97.5%) were
    # described by the model code"). Those bounds are reproduced here so
    # a population simulation cannot wander outside the published
    # support -- this matters most for Frac and PB, which are fractions
    # and whose Table 3 upper bounds are truncated at 1.00 and 0.99.
    frac <- min(max(exp(lfrac + etalfrac), 0.75), 1)          # Mi 2024 Table 3 (Frac bounds 0.75, 1.00)
    pb <- min(max(exp(lfb + etalfb), 0.36), 0.99)             # Mi 2024 Table 3 (PB bounds 0.36, 0.99)
    fu <- 1 - pb                                              # Suppl 3.1 `CAfree = CA*(1-PB)`

    kp_liver <- min(max(exp(lkp_liver + etalkp_liver), 2.29), 10.37)      # Mi 2024 Table 3 (PL bounds)
    kp_kidney <- min(max(exp(lkp_kidney + etalkp_kidney), 2.62), 11.85)   # Mi 2024 Table 3 (PK bounds)
    kp_muscle <- min(max(exp(lkp_muscle + etalkp_muscle), 0.34), 1.56)    # Mi 2024 Table 3 (PM bounds)
    kp_adipose <- min(max(exp(lkp_adipose + etalkp_adipose), 0.48), 2.17) # Mi 2024 Table 3 (PF bounds)
    kp_other <- exp(lkp_other)                                # not carried through Monte Carlo

    # Per-kg clearances scaled to the individual (L/h).
    cl_renal <- min(max(exp(lcl_renal + etalcl_renal), 0.04), 0.20) * WT  # Mi 2024 Table 3 (KurineC bounds); Suppl 3.1 `Kurine = KurineC * BW`
    cl_nonren <- exp(lcl_nonren) * WT                         # Suppl 3.1 `RML=KML*CVL*BW`

    # =================================================================
    # Concentrations (ug/mL == mg/L). States hold amounts in mg.
    # =================================================================
    # Emergent venous concentrations leaving each perfusion-limited
    # organ, Suppl 3.1 `CVL = AL/(VL * PL)` and siblings.
    cv_liver <- liver / (v_liver * kp_liver)
    cv_kidney <- kidney / (v_kidney * kp_kidney)
    cv_muscle <- muscle / (v_muscle * kp_muscle)
    cv_adipose <- adipose / (v_adipose * kp_adipose)
    cv_other <- other / (v_other * kp_other)

    # Arterial (== whole-pool) blood concentration, Suppl 3.1
    # `CA = AA/Vblood`.
    c_arterial <- blood / v_blood

    # Absorption flux (mg/h), Suppl 3.1 `Rim = Kim*Amtsiteim`.
    r_absorb <- ka * depot

    # Mixed venous concentration. Unlike the sibling Mi 2023 model this
    # is ALGEBRAIC, not a state: the published code computes it as the
    # flow-weighted mean of the organ outflows plus the absorbed
    # intramuscular drug, Suppl 3.1
    # `CV = ((QL*CVL+QK*CVK+QM*CVM+QOT*CVOT+QF*CVF+Rim)/QC)`. All of the
    # blood volume (0.06 of body weight) sits in the single arterial
    # pool. Supplementary Table 3 instead prints venous blood as a
    # dynamic state and omits the absorption input entirely; the code is
    # the runnable and self-consistent record and is what is encoded.
    # See vignette Errata.
    c_venous <-
      (q_liver * cv_liver +
         q_kidney * cv_kidney +
         q_muscle * cv_muscle +
         q_adipose * cv_adipose +
         q_other * cv_other +
         r_absorb) / q_co

    # Elimination fluxes (mg/h).
    r_renal <- cl_renal * cv_kidney                   # Suppl 3.1 `Rurine = Kurine * CVK`
    r_metab <- cl_nonren * cv_liver                   # Suppl 3.1 `RML=KML*CVL*BW`

    # =================================================================
    # ODE system (Supplementary Material 3.1; Supplementary Table 3)
    # =================================================================
    # Intramuscular absorption. `depot` is the fast injection-site pool
    # (Suppl 3.1 `Amtsiteim`) and `depot2` the slow-release depot
    # (Suppl 3.1 `DOSEppgim`). The dose is split between them by the
    # bioavailability terms below rather than by an in-ODE input rate.
    d/dt(depot2) <- -kdiss * depot2
    d/dt(depot) <- kdiss * depot2 - r_absorb

    # Blood pool, Suppl 3.1 `RA = QC*(CV-CA)`.
    d/dt(blood) <- q_co * (c_venous - c_arterial)

    # Perfusion-limited, well-stirred organs. Note these are driven by
    # TOTAL arterial concentration, not the unbound concentration --
    # Suppl 3.1 `RL=QL*(CA-CVL)` etc. use `CA`, and `CAfree` is computed
    # but never fed back into the disposition. This differs from the
    # sibling Mi_2023_cefquinome_pbpk.R, where tissues see `CAfree`.
    d/dt(liver) <- q_liver * (c_arterial - cv_liver) - r_metab
    d/dt(kidney) <- q_kidney * (c_arterial - cv_kidney) - r_renal
    d/dt(muscle) <- q_muscle * (c_arterial - cv_muscle)
    d/dt(adipose) <- q_adipose * (c_arterial - cv_adipose)
    d/dt(other) <- q_other * (c_arterial - cv_other)

    # Terminal elimination sinks, carried so the published mass balance
    # can be checked, Suppl 3.1 `d/dt(Aurine)` and `d/dt(AML)`.
    d/dt(urine) <- r_renal
    d/dt(metabolized) <- r_metab

    # Split the intramuscular dose: Frac is immediately available at the
    # injection site, 1 - Frac goes to the slow-release depot. Suppl 3.1
    # `Rpenim = Rinputim*(Frac)` and `Rppgim = Rinputim*(1-Frac)`. Dose
    # the same amount to both compartments and let these fractions
    # apportion it.
    f(depot) <- frac
    f(depot2) <- 1 - frac

    # =================================================================
    # Observations
    # =================================================================
    # Mi 2024 compares model output against plasma, liver, kidney, muscle
    # and fat. The code makes no blood-to-plasma correction, so the mixed
    # venous concentration is the paper's plasma prediction; it is also
    # the quantity the code accumulates for AUC (`d/dt(AUCCV) = CV`).
    Cc <- c_venous
    # Total tissue concentrations, the quantities compared against the
    # maximum residue limits (liver 0.303 ug/g, kidney 0.084 ug/g).
    Cliver <- liver / v_liver
    Ckidney <- kidney / v_kidney
    Cmuscle <- muscle / v_muscle
    Cadipose <- adipose / v_adipose
    # Unbound arterial concentration, Suppl 3.1 `CAfree = CA*(1-PB)`.
    # This is the driver of the fAUC/MIC index Mi 2024 selected as the
    # best PK/PD parameter, and the concentration that replaces Cstastic
    # in the integrated PBPK/PD model.
    Cfree <- c_arterial * fu

    Cc ~ prop(propSd)
  })
}
