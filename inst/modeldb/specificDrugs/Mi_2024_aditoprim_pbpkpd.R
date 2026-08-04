Mi_2024_aditoprim_pbpkpd <- function() {
  description <- paste(
    "Veterinary (pig) + in vitro (Streptococcus suis ATCC 49619).",
    "Integrated PBPK/PD model for the diaminopyrimidine aditoprim (ADP)",
    "after intramuscular injection in swine (Mi et al. 2024, Front",
    "Pharmacol 15:1378034). The disposition layer is the seven-compartment",
    "whole-body PBPK model of modellib('Mi_2024_aditoprim_pbpk'); the",
    "effect layer is a two-subpopulation semi-mechanistic",
    "pharmacodynamic model fitted in Monolix 2018R1 to in vitro",
    "time-killing curves against S. suis ATCC 49619 (MIC 0.5 ug/mL).",
    "Susceptible (S) and resistant (R) subpopulations share one logistic",
    "growth term capped at the stationary-phase count Bmax and are killed",
    "by a common sigmoid Emax function that differs only in potency",
    "(EC50_R = 1.63 mg/L is 2.4-fold EC50_S = 0.685 mg/L), so resistant",
    "organisms are selected under sub-optimal exposure. The starting",
    "inoculum is split by the fixed mutation frequency MF = 1e-3. The two",
    "layers are coupled one-way: the UNBOUND plasma concentration",
    "predicted by the PBPK model replaces the static broth concentration",
    "Cstastic of the time-kill experiment, and the bacteria do not feed",
    "back on the drug. The unbound (not total) concentration is the",
    "driver -- that is what the published code computes as CAfree and is",
    "the only reading that reproduces all four dose levels of Mi 2024",
    "Figure 6; see the vignette. Used to show that 20 mg/kg every 12 h",
    "for 3 days eradicates both subpopulations while 5 mg/kg selects for",
    "resistance. Neither layer reports between-animal or between-replicate",
    "random effects on the PD parameters, so the PD block is",
    "typical-value only and propSd is a placeholder."
  )
  reference <- paste(
    "Mi K, Sun L, Zhang L, Tang A, Tian X, Hou Y, Sun L, Huang L.",
    "A physiologically based pharmacokinetic/pharmacodynamic model to",
    "determine dosage regimens and withdrawal intervals of aditoprim",
    "against Streptococcus suis. Front Pharmacol. 2024;15:1378034.",
    "doi:10.3389/fphar.2024.1378034.",
    "PBPK equations transcribed from Supplementary Material section 3.1",
    "(Berkeley Madonna code) and Supplementary Table 3, with final values",
    "from Table 2 and Monte Carlo distributions from Table 3; PD equations",
    "from Supplementary Material section 3.2 (Monolix code) and Eq 2-3,",
    "with values from Table 4.",
    sep = " "
  )
  vignette <- "Mi_2024_aditoprim"

  # `metabolized` is the cumulative amount cleared by hepatic metabolism
  # (`AML` in the Supplementary Material code); a terminal bookkeeping
  # sink carried so the published mass balance can be checked. See
  # Mi_2024_aditoprim_pbpk.R for the full note. `S` and `R` are the
  # canonical susceptible / resistant bacterial subpopulations and need
  # no declaration.
  paper_specific_compartments <- c("metabolized")

  # Time in hours, doses in mg, concentrations in ug/mL. Drug states hold
  # amounts in mg and volumes are in L, so amount/volume is mg/L == ug/mL.
  # The two bacterial states are CONCENTRATIONS (CFU/mL), not amounts --
  # they live in the in vitro broth system of the time-kill experiment,
  # not in a swine tissue, and are never dosed.
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Mi 2024 Table 2 fixes BW = 30 kg, and the Supplementary Material",
        "code sets `BW=30`. Every organ volume and blood flow is a fixed",
        "fraction of BW and both clearances are per-kg constants, so the",
        "PBPK layer scales isometrically with WT. The PD layer has no",
        "body-weight dependence -- it was fitted to in vitro broth",
        "cultures. See Mi_2024_aditoprim_pbpk.R for the 25 kg / 30 kg",
        "discrepancy note."
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
    metabolized = list(analyte = "aditoprim", units = "mg", specimen = "not applicable", verified = TRUE),
    S           = list(analyte = "Streptococcus suis ATCC 49619 (susceptible subpopulation)", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    R           = list(analyte = "Streptococcus suis ATCC 49619 (resistant subpopulation)", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "pig (swine) + in vitro (Streptococcus suis ATCC 49619)",
    n_subjects     = 92L,
    n_studies      = 4L,
    age_range      = "grower / finisher swine",
    weight_range   = "model reference BW = 30 kg",
    disease_state  = "healthy swine; Streptococcus suis infection is simulated, not observed in vivo",
    dose_range     = "0, 2.5, 5, 10, 15, 20 and 25 mg/kg intramuscular once or twice daily were simulated; 5, 12.5, 15 and 20 mg/kg every 12 h for 3 consecutive days are shown in Figure 6",
    regions        = "China",
    notes          = paste(
      "The PBPK layer carries the swine population of",
      "modellib('Mi_2024_aditoprim_pbpk') -- four previously published",
      "swine datasets (Mi 2024 Table 1) digitised with WebPlotDigitizer.",
      "The PD layer is INDEPENDENT of that population: it was fitted by",
      "non-linear mixed effects in Monolix 2018R1 to in vitro",
      "time-killing curves of ADP against S. suis ATCC 49619 performed in",
      "triplicate (Mi 2024 section 2.6), with the MIC determined per CLSI",
      "2018 as 0.5 ug/mL. There is therefore no single cohort of animals",
      "in which the coupled model was observed; the integration is an in",
      "silico projection performed in Mlxplore 2018a (Mi 2024",
      "section 2.7). No between-animal or between-replicate variance was",
      "reported for any PD parameter -- Table 4 gives relative standard",
      "errors of the estimates only, which are parameter uncertainty and",
      "NOT between-subject variability, so no PD etas are encoded."
    )
  )

  ini({
    # =================================================================
    # PBPK layer -- identical to Mi_2024_aditoprim_pbpk.R. See that file
    # for the full provenance notes on each value (Table 2 final values
    # versus the pre-optimisation values in the supplement code, and the
    # Table 3 Monte Carlo parameterisation).
    # =================================================================
    lka <- log(1.3)
    label("Absorption rate constant Kim, injection site to venous mixing pool (1/h)")   # Mi 2024 Table 2 (Kim = 1.3 /h, Model fitting)

    lkdiss <- log(0.0118)
    label("Dissolution rate constant Kdiss, slow depot to injection site (1/h)")        # Mi 2024 Table 2 (Kdiss = 0.0118 /h, Model fitting)

    lfrac <- log(0.92)
    label("Fraction of the intramuscular dose entering the fast injection-site pool (unitless)") # Mi 2024 Table 2 (Frac = 0.92, Model fitting)

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

    lfb <- log(0.82)
    label("Fraction of aditoprim bound to plasma protein PB (unitless)")                # Mi 2024 Table 2 (PB = 0.82, Model fitting); Suppl 3.1 `CAfree = CA*(1-PB)`

    lcl_renal <- log(0.1)
    label("Renal clearance KurineC (L/h/kg), applied to kidney venous concentration")   # Mi 2024 Table 2 (KurineC = 0.1 L/h/kg)

    lcl_nonren <- log(0.01)
    label("Hepatic metabolic clearance KML (L/h/kg), applied to liver venous concentration") # Mi 2024 Table 2 (KML = 0.01 L/h/kg)

    qcc <- fixed(4.944)
    label("Cardiac output QCC (L/h/kg)")                                                # Mi 2024 Table 2 (QCC = 4.944 L/h/kg, Upton 2008)

    # =================================================================
    # Semi-mechanistic PD layer (Mi 2024 Table 4; Supplementary Material
    # section 3.2 Monolix code). All values are the final Monolix
    # estimates; the bracketed percentages in Table 4 are relative
    # standard errors of the estimate (parameter uncertainty), not
    # between-subject variability, so none of them becomes an eta.
    # =================================================================
    lkgrow <- log(0.833)
    label("Net natural bacterial growth rate constant Kgrowth (1/h)")                   # Mi 2024 Table 4 (Kgrowth = 0.833 1/h, R.S.E. 15.6%); Suppl 3.2 `kg`

    lbmax <- log(4.61e7)
    label("Stationary-phase bacterial count Bmax (CFU/mL)")                             # Mi 2024 Table 4 (Bmax = 4.61e7 CFU/mL, no R.S.E. reported); Suppl 3.2 `Bmax`

    lemax <- log(1.45)
    label("Maximum aditoprim kill rate constant Emax (1/h)")                            # Mi 2024 Table 4 (Emax = 1.45 1/h, R.S.E. 4.47%); Suppl 3.2 `Emax`

    lec50_s <- log(0.685)
    label("Aditoprim concentration giving 50% of Emax in the susceptible subpopulation EC50_S (mg/L)") # Mi 2024 Table 4 (EC50_S = 0.685 mg/L, R.S.E. 8.12%); Suppl 3.2 `EC50_E`

    lec50_r <- log(1.63)
    label("Aditoprim concentration giving 50% of Emax in the resistant subpopulation EC50_R (mg/L)")   # Mi 2024 Table 4 (EC50_R = 1.63 mg/L, R.S.E. 4.2%); Suppl 3.2 `EC50_R`

    hill <- 2.5
    label("Sigmoid (Hill) coefficient gamma of the antibacterial effect (unitless)")    # Mi 2024 Table 4 (gamma = 2.5, R.S.E. 6.5%); Suppl 3.2 `gama`

    # MF is reported as a FIXED parameter in Mi 2024 Table 4. It is not
    # used as a rate anywhere in the Monolix code; it is the fraction of
    # the starting inoculum assigned to the resistant subpopulation. The
    # code hard-codes the resulting split as `E_0 =99900` and `R_0=100`,
    # which is exactly 1e5 CFU/mL partitioned by MF = 1e-3. Both are
    # carried here as parameters so the inoculum can be varied.
    mf <- fixed(1e-3)
    label("Mutation frequency MF, resistant fraction of the starting inoculum (unitless)") # Mi 2024 Table 4 (MF = 10^-3, fixed); Suppl 3.2 `E_0 =99900`, `R_0=100`

    log10inoc <- fixed(5)
    label("Total starting inoculum (log10 CFU/mL)")                                     # Suppl 3.2 `E_0 =99900` + `R_0=100` = 1e5 CFU/mL

    mic <- fixed(0.5)
    label("Minimum inhibitory concentration of aditoprim against S. suis ATCC 49619 (ug/mL)") # Mi 2024 section 3.4 (MIC = 0.5 ug/mL, CLSI 2018)

    # =================================================================
    # Between-animal variability on the PBPK layer, from the Monte Carlo
    # distributions of Mi 2024 Table 3. Assumed CVs, not estimates, so
    # all fixed(). See Mi_2024_aditoprim_pbpk.R for the derivation.
    # =================================================================
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

    propSd <- fixed(0.10)
    label("Proportional residual error placeholder (fraction)")                         # not reported in Mi 2024; placeholder only
  })

  model({
    # =================================================================
    # Swine physiology (Mi 2024 Table 2, Upton 2008; Supplementary
    # Material 3.1). QFC and VFC are reported ONLY in the supplement
    # code -- Table 2 omits both even though fat is a modelled output
    # tissue. The rest-of-body flow and volume are complements, as the
    # code computes them. See Mi_2024_aditoprim_pbpk.R and the vignette
    # Errata for the discrepancies against the printed Table 2 values.
    # =================================================================
    q_co <- qcc * WT                                  # Suppl 3.1 `QC=QCAR*BW`

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

    v_liver <- 0.0294 * WT                            # Mi 2024 Table 2 (VLC = 0.0294)
    v_kidney <- 0.004 * WT                            # Mi 2024 Table 2 (VKC = 0.004)
    v_muscle <- 0.4 * WT                              # Mi 2024 Table 2 (VMC = 0.4)
    v_adipose <- 0.3 * WT                             # Suppl 3.1 `VFC=0.3`; absent from Table 2
    v_blood <- 0.06 * WT                              # Suppl 3.1 `VBC=0.06` = Table 2 VartC 0.016 + VvenC 0.044
    v_other <- (1 - 0.0294 - 0.004 - 0.4 - 0.3 - 0.06) * WT
    # Suppl 3.1 `VOT=BW-(VL+VK+VM+VF+Vblood)`

    # =================================================================
    # Individual chemical-specific parameters. The Monte Carlo draws are
    # confined to the 2.5 / 97.5 percentile bounds printed in Mi 2024
    # Table 3, per section 2.5.
    # =================================================================
    ka <- exp(lka)
    kdiss <- exp(lkdiss)

    frac <- min(max(exp(lfrac + etalfrac), 0.75), 1)          # Mi 2024 Table 3 (Frac bounds 0.75, 1.00)
    pb <- min(max(exp(lfb + etalfb), 0.36), 0.99)             # Mi 2024 Table 3 (PB bounds 0.36, 0.99)
    fu <- 1 - pb                                              # Suppl 3.1 `CAfree = CA*(1-PB)`

    kp_liver <- min(max(exp(lkp_liver + etalkp_liver), 2.29), 10.37)      # Mi 2024 Table 3 (PL bounds)
    kp_kidney <- min(max(exp(lkp_kidney + etalkp_kidney), 2.62), 11.85)   # Mi 2024 Table 3 (PK bounds)
    kp_muscle <- min(max(exp(lkp_muscle + etalkp_muscle), 0.34), 1.56)    # Mi 2024 Table 3 (PM bounds)
    kp_adipose <- min(max(exp(lkp_adipose + etalkp_adipose), 0.48), 2.17) # Mi 2024 Table 3 (PF bounds)
    kp_other <- exp(lkp_other)                                # not carried through Monte Carlo

    cl_renal <- min(max(exp(lcl_renal + etalcl_renal), 0.04), 0.20) * WT  # Mi 2024 Table 3 (KurineC bounds)
    cl_nonren <- exp(lcl_nonren) * WT                         # Suppl 3.1 `RML=KML*CVL*BW`

    # =================================================================
    # Concentrations (ug/mL == mg/L)
    # =================================================================
    cv_liver <- liver / (v_liver * kp_liver)
    cv_kidney <- kidney / (v_kidney * kp_kidney)
    cv_muscle <- muscle / (v_muscle * kp_muscle)
    cv_adipose <- adipose / (v_adipose * kp_adipose)
    cv_other <- other / (v_other * kp_other)

    c_arterial <- blood / v_blood                     # Suppl 3.1 `CA = AA/Vblood`
    r_absorb <- ka * depot                            # Suppl 3.1 `Rim = Kim*Amtsiteim`

    # Mixed venous concentration -- algebraic, not a state. Suppl 3.1
    # `CV = ((QL*CVL+QK*CVK+QM*CVM+QOT*CVOT+QF*CVF+Rim)/QC)`.
    c_venous <-
      (q_liver * cv_liver +
         q_kidney * cv_kidney +
         q_muscle * cv_muscle +
         q_adipose * cv_adipose +
         q_other * cv_other +
         r_absorb) / q_co

    r_renal <- cl_renal * cv_kidney                   # Suppl 3.1 `Rurine = Kurine * CVK`
    r_metab <- cl_nonren * cv_liver                   # Suppl 3.1 `RML=KML*CVL*BW`

    # =================================================================
    # PBPK ODE system (Supplementary Material 3.1)
    # =================================================================
    d/dt(depot2) <- -kdiss * depot2
    d/dt(depot) <- kdiss * depot2 - r_absorb
    d/dt(blood) <- q_co * (c_venous - c_arterial)     # Suppl 3.1 `RA = QC*(CV-CA)`

    # Driven by TOTAL arterial concentration, as coded.
    d/dt(liver) <- q_liver * (c_arterial - cv_liver) - r_metab
    d/dt(kidney) <- q_kidney * (c_arterial - cv_kidney) - r_renal
    d/dt(muscle) <- q_muscle * (c_arterial - cv_muscle)
    d/dt(adipose) <- q_adipose * (c_arterial - cv_adipose)
    d/dt(other) <- q_other * (c_arterial - cv_other)

    d/dt(urine) <- r_renal
    d/dt(metabolized) <- r_metab

    f(depot) <- frac                                  # Suppl 3.1 `Rpenim = Rinputim*(Frac)`
    f(depot2) <- 1 - frac                             # Suppl 3.1 `Rppgim = Rinputim*(1-Frac)`

    # =================================================================
    # Semi-mechanistic PD layer (Mi 2024 Eq 2-3; Suppl 3.2)
    # =================================================================
    # PD driver. Mi 2024 section 2.7: "the predicted dynamic
    # concentrations from the PBPK model replace the static ADP
    # concentration (Cstastic) in the semi-mechanistic PD model". The
    # UNBOUND concentration is the one that does so: the published code
    # computes `CAfree = CA*(1-PB)` and never uses it anywhere in the
    # disposition, and Mi 2024 section 3.5.1 selects fAUC/MIC (r2 = 0.99)
    # as the best PK/PD index. Driving the bacteria with the total
    # concentration instead eradicates S. suis at EVERY simulated dose
    # including 5 mg/kg, contradicting Figure 6 and the paper's central
    # conclusion; the unbound driver reproduces all four Figure 6 dose
    # levels. See the vignette for that comparison.
    cdrive <- c_arterial * fu

    kgrow <- exp(lkgrow)
    bmax <- exp(lbmax)
    emax <- exp(lemax)
    ec50_s <- exp(lec50_s)
    ec50_r <- exp(lec50_r)

    # Shared logistic growth term. Mi 2024 Eq 2-3 as typeset drop the
    # leading "1 -" of this factor, which would make growth vanish at low
    # density and peak at carrying capacity -- the reverse of a
    # stationary-phase limit. The Monolix code in Supplementary Material
    # 3.2 prints it correctly as `kg*(1-(E+R)/Bmax)`, and that is what is
    # encoded. See vignette Errata.
    logistic_growth <- kgrow * (1 - (S + R) / bmax)

    # Sigmoid Emax killing, common Emax and gamma, subpopulation-specific
    # potency. Suppl 3.2 `((Emax*Cc^gama)/(Cc^gama+EC50_E^gama))`.
    kill_s <- emax * cdrive^hill / (cdrive^hill + ec50_s^hill)
    kill_r <- emax * cdrive^hill / (cdrive^hill + ec50_r^hill)

    d/dt(S) <- logistic_growth * S - kill_s * S       # Mi 2024 Eq 2; Suppl 3.2 `ddt_E`
    d/dt(R) <- logistic_growth * R - kill_r * R       # Mi 2024 Eq 3; Suppl 3.2 `ddt_R`

    # Starting inoculum split by the mutation frequency. Suppl 3.2
    # `E_0 =99900` and `R_0=100`, i.e. 1e5 CFU/mL split by MF = 1e-3.
    S(0) <- 10^log10inoc * (1 - mf)
    R(0) <- 10^log10inoc * mf

    # =================================================================
    # Observations
    # =================================================================
    Cc <- c_venous
    Cliver <- liver / v_liver
    Ckidney <- kidney / v_kidney
    Cmuscle <- muscle / v_muscle
    Cadipose <- adipose / v_adipose
    Cfree <- cdrive                                   # Suppl 3.1 `CAfree = CA*(1-PB)`
    # Free-concentration-to-MIC ratio. Mi 2024 section 3.5.1 selected
    # fAUC/MIC (r2 = 0.99) over %fT>MIC (r2 = 0.98) as the PK/PD index
    # best related to the bacteriological effect. Integrating this over
    # 24 h gives fAUC24/MIC; the fraction of the dosing interval for
    # which it exceeds 1 gives %fT>MIC.
    Cfree_over_mic <- cdrive / mic

    # Total bacterial count and its log10, the quantity plotted in
    # Mi 2024 Figures 4 and 6. Suppl 3.2 `B=E+R`, `output =B`.
    bacteria <- S + R
    log_cfu <- log10(bacteria)

    Cc ~ prop(propSd)
  })
}
