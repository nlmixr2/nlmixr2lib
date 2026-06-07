Bulitta_2010_colistin_URMC1 <- function() {
  description <- "In vitro (Pseudomonas aeruginosa URMC1; S-ADAPT fit). Mechanism-based population pharmacodynamic model for colistin static time-kill experiments against the clinical P. aeruginosa isolate URMC1, used in Bulitta 2010 alongside PAO1 and URMC2 to test whether the structural model generalises beyond the reference strain. Identical structure to Bulitta_2010_colistin_PAO1: three pre-existing bacterial subpopulations (S / I / R) with different second-order colistin killing rate constants; a lag compartment that holds susceptible cells at t = 0 and transfers them into the replicating susceptible compartment at rate klag; a saturable Michaelis-Menten-style growth function parameterised by POPmax and per-subpopulation low-density growth half-lives; first-order natural death at kd = 0.3 /h; receptor-occupancy submodel for Mg2+/Ca2+ competition followed by a Hill-10 mapping into effective colistin at the target site; and a signal-molecule compartment tracking CFUALL whose Hill-1 inhibition of replication and killing produces the observed inoculum effect. Drug concentration Ccolistin (mg/L) and cation concentration Ccations (umol/L) are external time-varying covariates from the data. Random effects (eta) are NOT included -- URMC1 is reported in the S-ADAPT column of Table 1 without IIV."
  reference <- paste(
    "Bulitta JB, Yang JC, Yohonn L, Ly NS, Brown SV, D'Hondt RE, Jusko WJ, Forrest A, Tsuji BT. (2010).",
    "Attenuation of colistin bactericidal activity by high inoculum of Pseudomonas aeruginosa",
    "characterized by a new mechanism-based population pharmacodynamic model.",
    "Antimicrobial Agents and Chemotherapy 54(5):2051-2062.",
    "doi:10.1128/AAC.00881-09.",
    sep = " "
  )
  vignette <- "Bulitta_2010_colistin"
  units <- list(time = "hour", dosing = "mg/L (colistin in broth)", concentration = "log10 CFU/mL (observation); mg/L (drug covariate); umol/L (cation covariate)")

  depends <- c("Ccolistin", "Ccations")
  paper_specific_compartments <- c("bact_slag", "bact_s", "bact_i", "bact_r", "signal")

  covariateData <- list(
    Ccolistin = list(
      description        = "Colistin concentration in growth medium",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Static concentration in supplemented LB broth held constant over the 24 h experiment. URMC1 was studied at 9 colistin concentrations up to 64 mg/L (64x the LB-broth MIC of 1 mg/L). In-vitro experimental input -- not in inst/references/covariate-columns.md.",
      source_name        = "Ccolistin (Bulitta 2010 Methods, Time-kill experiments)"
    ),
    Ccations = list(
      description        = "Sum of Mg2+ and Ca2+ molar concentration in the growth medium",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Default 1138 umol/L recovers the published LB-broth supplemented condition (0.514 mmol/L Mg2+ + 0.624 mmol/L Ca2+ per Bulitta 2010 Table 1 footnote f).",
      source_name        = "Ccations (Bulitta 2010 Methods + Table 1 footnote f)"
    )
  )

  population <- list(
    species             = "in vitro (Pseudomonas aeruginosa, URMC1 clinical isolate)",
    n_subjects          = NA_integer_,
    n_studies           = 1L,
    organism            = "Pseudomonas aeruginosa URMC1 (clinical isolate from the University of Rochester Medical Center; LB-broth MIC 1.0 mg/L, MHB MIC 0.5 mg/L; four-fold more colistin-susceptible than PAO1)",
    system              = "Static time-kill experiments at 37 C in supplemented cation-adjusted LB broth (25 mg/L Ca2+, 12.5 mg/L Mg2+); 20 mL cultures in constant shaking water bath",
    medium              = "Supplemented LB broth (cation-adjusted)",
    duration            = "24 h with sampling at 0, 0.5, 1, 2, 4, 6, 8, and 24 h",
    inoculum_range      = "10^6, 10^8, and 10^9 CFU/mL",
    mic_values          = c(colistin_LB_broth = "1.0 mg/L", colistin_MHB = "0.5 mg/L"),
    regimens            = "9 colistin concentrations up to 64 mg/L (64x MIC in LB broth) at each initial inoculum; antibiotic-free growth controls; viable counts in duplicate.",
    notes               = "In-vitro pharmacodynamic study; no human or animal subjects. URMC1 is the first of two clinical P. aeruginosa isolates used to externally qualify the structural model developed against PAO1; the same equations are re-fit in S-ADAPT with strain-specific parameter values (Table 1, Strain URMC1 column). Random effects are not reported for the URMC strains. The companion model files Bulitta_2010_colistin_PAO1 (NONMEM primary, 11 colistin concentrations, 10^4-10^9 inocula) and Bulitta_2010_colistin_URMC2 (a different clinical isolate) cover the full strain set described in Bulitta 2010."
  )

  ini({
    # =============================================================
    # Bacterial growth -- URMC1 / S-ADAPT (Bulitta 2010 Table 1)
    # =============================================================
    t12_kg_low_s <- 17.9
    label("Susceptible-population low-density growth half-life (min; t1/2(kg,low CFU)_S)")  # Bulitta 2010 Table 1, Strain URMC1
    t12_kg_low_i <- 21.7
    label("Intermediate-population low-density growth half-life (min; t1/2(kg,low CFU)_I)")  # Bulitta 2010 Table 1, URMC1
    t12_kg_low_r <- 86.8
    label("Least-susceptible-population low-density growth half-life (min; t1/2(kg,low CFU)_R)")  # Bulitta 2010 Table 1, URMC1

    kd_nat <- fixed(0.3)
    label("First-order natural death rate constant (1/h; kd; FIXED to Meagher 2004 ref 45)")  # Bulitta 2010 Table 1 footnote (a)

    t12_klag <- 1.36
    label("Growth-lag half-life (h; t1/2(klag); klag = ln(2)/t12_klag)")  # Bulitta 2010 Table 1, URMC1

    log10_popmax <- 9.68
    label("log10 maximum population size (log10 CFU/mL; POPmax without signal inhibition)")  # Bulitta 2010 Table 1, URMC1

    # =============================================================
    # Initial inoculum and subpopulation fractions
    # =============================================================
    log10_cfu0 <- 6.35
    label("log10 initial inoculum (default = 10^6-class value for URMC1)")  # Bulitta 2010 Table 1, Log10 CFUo,6, URMC1

    log10_fr_i <- -3.97
    label("log10 intermediate-population fraction of initial inoculum (unitless; Log10 FrI)")  # Bulitta 2010 Table 1, URMC1
    log10_fr_r <- -7.16
    label("log10 least-susceptible-population fraction of initial inoculum (unitless; Log10 FrR)")  # Bulitta 2010 Table 1, URMC1

    # =============================================================
    # Signal molecules
    # =============================================================
    log10_fr_sig <- -1.96
    label("log10 initial signal-molecule concentration relative to inoculum (ml/CFU; Log10 FrSig)")  # Bulitta 2010 Table 1, URMC1
    log10_ic50_sig <- 6.27
    label("log10 signal-molecule concentration for 50% maximal inhibition (Log10 IC50)")  # Bulitta 2010 Table 1, URMC1
    t12_kdeg <- 1.30
    label("Signal-molecule degradation half-life (h; t1/2(kdeg))")  # Bulitta 2010 Table 1, URMC1
    imax_rep <- 0.569
    label("Maximal fractional inhibition of bacterial replication by signal molecules (ImaxRep, unitless)")  # Bulitta 2010 Table 1, URMC1
    imax_kill <- 0.999
    label("Maximal fractional inhibition of bacterial killing by signal molecules (ImaxKill, unitless)")  # Bulitta 2010 Table 1, URMC1

    # =============================================================
    # Receptor occupancy and bacterial killing
    # =============================================================
    ec50_rec <- 0.316
    label("Fraction of receptors not occupied by Mg2+/Ca2+ giving 50% effective colistin (EC50, unitless)")  # Bulitta 2010 Table 1, URMC1
    hill_rec <- fixed(10)
    label("Receptor-occupancy Hill coefficient (unitless; gamma; FIXED at 10 per Table 1 footnote e)")  # Bulitta 2010 Table 1 footnote (e)
    kdiss_cation <- fixed(200)
    label("Receptor dissociation constant for Mg2+/Ca2+ (umol/L; KdCations; FIXED at 200)")  # Bulitta 2010 Methods, Table 1 footnote (g)
    kdiss_colistin <- fixed(0.3)
    label("Receptor dissociation constant for colistin (umol/L; KdColistin; FIXED at 0.3)")  # Bulitta 2010 Methods, Table 1 footnote (g)
    mw_colistin <- fixed(1.163)
    label("Mean colistin A+B molar mass (mg/umol = g/mmol; equivalent to 1163 g/mol; FIXED)")  # Bulitta 2010 Methods, paragraph after Eq. 1

    lk2s <- log(7.88)
    label("Susceptible-population second-order killing rate constant (L/(mg*h); k2S)")  # Bulitta 2010 Table 1, URMC1
    lk2i <- log(0.627)
    label("Intermediate-population second-order killing rate constant (L/(mg*h); k2I)")  # Bulitta 2010 Table 1, URMC1
    lk2r <- log(0.00573)
    label("Least-susceptible-population second-order killing rate constant (L/(mg*h); k2R)")  # Bulitta 2010 Table 1, URMC1

    # =============================================================
    # Residual error
    # =============================================================
    # The paper reports the SD of additive residual error on log10
    # scale in Table 1 footnote (h) as 0.478 from S-ADAPT (well-
    # comparable to 0.474 from NONMEM; PAO1 NONMEM was the primary
    # fit and the S-ADAPT value carries over to URMC1 and URMC2).
    addSd <- 0.478
    label("Additive residual SD on log10 total viable count (log10 CFU/mL; SD_CFU)")  # Bulitta 2010 Table 1 footnote (h), S-ADAPT
  })

  model({
    # ---- exponentiate log-transformed structural parameters ----
    k2s <- exp(lk2s)
    k2i <- exp(lk2i)
    k2r <- exp(lk2r)

    # ---- rate-constant derivations from half-lives ----
    kg_low_s <- log(2) * 60 / t12_kg_low_s
    kg_low_i <- log(2) * 60 / t12_kg_low_i
    kg_low_r <- log(2) * 60 / t12_kg_low_r
    klag <- log(2) / t12_klag
    kdeg <- log(2) / t12_kdeg

    # ---- linear-scale population parameters ----
    popmax   <- 10 ^ log10_popmax
    cfu0     <- 10 ^ log10_cfu0
    fr_i     <- 10 ^ log10_fr_i
    fr_r     <- 10 ^ log10_fr_r
    fr_sig   <- 10 ^ log10_fr_sig
    ic50_sig <- 10 ^ log10_ic50_sig

    # ---- growth-function back-derivation (CFUm and VGmax_X) ----
    cfum    <- popmax * kd_nat / (kg_low_s - kd_nat)
    vgmax_s <- kg_low_s * cfum
    vgmax_i <- kg_low_i * cfum
    vgmax_r <- kg_low_r * cfum

    # ---- total viable bacteria (Eq. 3) ----
    cfu_all <- bact_slag + bact_s + bact_i + bact_r

    # ---- receptor occupancy (Eq. 1) ----
    ccolistin_um <- Ccolistin / mw_colistin
    frcations    <- Ccations / (kdiss_cation + Ccations + (kdiss_cation / kdiss_colistin) * ccolistin_um)

    # ---- effective colistin concentration at the target site (Eq. 2) ----
    not_occ       <- 1 - frcations
    not_occ_h     <- not_occ ^ hill_rec
    ec50_h        <- ec50_rec ^ hill_rec
    ccolistin_eff <- (not_occ_h / (ec50_h + not_occ_h)) * Ccolistin

    # ---- signal-molecule inhibition (Eqs. 4, 5) ----
    inh_factor <- signal / (ic50_sig + signal)
    inh_kill   <- 1 - imax_kill * inh_factor
    inh_rep    <- 1 - imax_rep * inh_factor

    # ---- per-subpopulation growth and killing rates ----
    growth_factor <- 1 / (cfum + cfu_all)
    growth_s <- inh_rep * vgmax_s * growth_factor
    growth_i <- inh_rep * vgmax_i * growth_factor
    growth_r <- inh_rep * vgmax_r * growth_factor

    kill_s <- inh_kill * k2s * ccolistin_eff
    kill_i <- inh_kill * k2i * ccolistin_eff
    kill_r <- inh_kill * k2r * ccolistin_eff

    # ---- ODE system (Eqs. 6-10) ----
    d/dt(bact_slag) <- (-klag - kill_s) * bact_slag
    d/dt(bact_s)    <- (growth_s - kd_nat - kill_s) * bact_s + klag * bact_slag
    d/dt(bact_i)    <- (growth_i - kd_nat - kill_i) * bact_i
    d/dt(bact_r)    <- (growth_r - kd_nat - kill_r) * bact_r
    d/dt(signal)    <- kdeg * (cfu_all - signal)

    # ---- initial conditions ----
    cfu_i0    <- fr_i * cfu0
    cfu_r0    <- fr_r * cfu0
    cfu_slag0 <- cfu0 - cfu_i0 - cfu_r0

    bact_slag(0) <- cfu_slag0
    bact_s(0)    <- 0
    bact_i(0)    <- cfu_i0
    bact_r(0)    <- cfu_r0
    signal(0)    <- fr_sig * cfu0

    # ---- observation (Eq. 11) ----
    cfu_obs <- cfu_all + 1
    Cc      <- log10(cfu_obs)
    Cc      ~ add(addSd)
  })
}
