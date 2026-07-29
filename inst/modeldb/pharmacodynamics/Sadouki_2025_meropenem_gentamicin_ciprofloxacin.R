Sadouki_2025_meropenem_gentamicin_ciprofloxacin <- function() {
  description <- "In-vitro static-time-kill pharmacodynamic model for two- and three-way combinations of meropenem, gentamicin, and ciprofloxacin against Escherichia coli NCTC 12,241. Logistic bacterial growth (knet, Bmax) is killed by an Emax-Hill function for each drug; emergence of regrowth is captured by a time-decay term parameterised by BETA (loss of effect) and TAU (time-shape). Meropenem chemical degradation in CAMHB at 37.5 C is embedded as a first-order decay of the meropenem solution concentration. Drug-drug interactions are encoded as: a fixed -1 categorical shift on BETA whenever a 2- or 3-way combination is present (so the regrowth term reverses sign and effect is sustained), proportional reductions of -0.353 and -0.576 in ciprofloxacin IC50 in the presence of meropenem and gentamicin respectively (synergy on potency), and concentration-dependent Emax shifts of BETA for gentamicin and ciprofloxacin. The model is in-vitro PD only -- there is no human or animal PK component; bacterial counts (CFU/mL) are observed on log scale."
  reference <- "Sadouki Z, Wey EQ, Read L, Bayliss M, Noel A, Balakrishnan I, McHugh TD, Kloprogge F. Pharmacodynamic interactions among meropenem ciprofloxacin and gentamicin in an in-vitro model. Sci Rep. 2025 Nov 24;15:45244. doi:10.1038/s41598-025-29354-y."
  vignette <- "Sadouki_2025_meropenem_gentamicin_ciprofloxacin"
  units <- list(time = "hour", dosing = "mg/L", concentration = "mg/L")

  covariateData <- list(
    CONMED_MER = list(
      description        = "Indicator that meropenem is present in the experimental regimen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (meropenem absent)",
      notes              = "1 if the regimen includes meropenem, 0 otherwise. Drives the combination categorical effect on BETA and the proportional MER-presence effect on Cip IC50 (Sadouki 2025 Table 1, Drug interactions section). In-vitro experimental indicator -- not in inst/references/covariate-columns.md (the canonical register is for human pop-PK covariates and does not apply to this in-vitro PD model).",
      source_name        = "Mer regimen flag (paper Methods)"
    ),
    CONMED_GEN = list(
      description        = "Indicator that gentamicin is present in the experimental regimen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (gentamicin absent)",
      notes              = "1 if the regimen includes gentamicin, 0 otherwise. Drives the combination categorical effect on BETA and the proportional GEN-presence effect on Cip IC50 (Sadouki 2025 Table 1).",
      source_name        = "Gen regimen flag (paper Methods)"
    ),
    CONMED_CIP = list(
      description        = "Indicator that ciprofloxacin is present in the experimental regimen",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ciprofloxacin absent)",
      notes              = "1 if the regimen includes ciprofloxacin, 0 otherwise. Used together with CONMED_MER and CONMED_GEN to compute the combination indicator that triggers the -1 categorical shift on BETA (Sadouki 2025 Table 1).",
      source_name        = "Cip regimen flag (paper Methods)"
    ),
    CONMED_GEN_CC = list(
      description        = "Gentamicin concentration applied during the time-kill experiment (constant over the 24-h window)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Static concentration in CAMHB. The paper studied multiples of the MIC of 1 mg/L (range 0.25-16x MIC = 0.25-16 mg/L). Set to 0 when gentamicin is absent.",
      source_name        = "Gen concentration (paper Methods)"
    ),
    CONMED_CIP_CC = list(
      description        = "Ciprofloxacin concentration applied during the time-kill experiment (constant over the 24-h window)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Static concentration in CAMHB. The paper studied multiples of the MIC of 0.015 mg/L (range 0.25-16x MIC = 0.00375-0.24 mg/L). Set to 0 when ciprofloxacin is absent.",
      source_name        = "Cip concentration (paper Methods)"
    ),
    LOWINOC = list(
      description        = "Indicator of low-inoculum experiment",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (10^5 CFU/mL standard inoculum)",
      notes              = "1 if the starting inoculum was 10^3 CFU/mL, 0 if 10^5 CFU/mL (the modal experimental setup). Drives the additive shifts on B0 and Bmax (Sadouki 2025 Table 1, Inoculum effect rows).",
      source_name        = "Inoculum size (paper Methods)"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 1L,
    organism       = "Escherichia coli NCTC 12,241 (susceptible laboratory reference strain)",
    system         = "Static time-kill in 96-well plates, total volume 200 uL, biological duplicates with technical triplicates",
    medium         = "Cation-adjusted Mueller-Hinton broth (CAMHB)",
    temperature    = "37.5 C",
    duration       = "24 h, with hourly sampling for the first 8 h and at 24 h",
    mic_values     = c(meropenem = "0.03 mg/L", ciprofloxacin = "0.015 mg/L", gentamicin = "1 mg/L"),
    concentration_range = "0.25 to 16 x MIC for each drug",
    inoculum_options = c("10^3 CFU/mL (low)", "10^5 CFU/mL (standard)"),
    regimens        = "Seven antibiotic-containing regimens: monotherapy (Mer, Gen, Cip), two-way combinations (Mer+Gen, Mer+Cip, Gen+Cip), three-way combination (Mer+Gen+Cip), plus antibiotic-free growth controls",
    notes           = "In-vitro pharmacodynamic study; no human or animal subjects. Random effects (eta) in the structural model represent variability *between experimental replicates*, not between-subject IIV. CFU were enumerated on Mueller-Hinton agar (and on 2x and 8x MIC supplemented agar at 24 h to detect resistant subpopulations). See Sadouki 2025 Methods (page 2) and Figure 1 for the experimental design."
  )

  ini({
    # ---- Growth model (logistic with carrying capacity) ----
    lknet <- log(1.35)
    label("Net bacterial growth rate (knet, 1/h)")  # Sadouki 2025 Table 1, Growth model parameters
    b0    <- 5.29
    label("Baseline bacterial concentration (B0, log10 CFU/mL)")  # Sadouki 2025 Table 1
    bmax  <- 10
    label("Maximum carrying capacity (Bmax, log10 CFU/mL)")  # Sadouki 2025 Table 1

    # ---- Inoculum effect (low vs reference; reference = 10^5 CFU/mL) ----
    cat_b0_lowinoc   <- -0.326
    label("Categorical low-inoculum (10^3 vs 10^5) effect on B0 (log10 CFU/mL)")  # Sadouki 2025 Table 1, Inoculum effect rows
    cat_bmax_lowinoc <- -0.110
    label("Categorical low-inoculum (10^3 vs 10^5) effect on Bmax (log10 CFU/mL)")  # Sadouki 2025 Table 1

    # ---- Meropenem -- bacterial killing ----
    lemax_mer <- log(4.18)
    label("Meropenem maximum effect (Emax_Mer, 1/h)")  # Sadouki 2025 Table 1, Meropenem killing
    lic50_mer <- log(0.0781)
    label("Meropenem half-maximum-effect concentration (IC50_Mer, mg/L)")  # Sadouki 2025 Table 1
    lhill_mer <- log(2.76)
    label("Meropenem Hill exponent (hill_Mer, unitless)")  # Sadouki 2025 Table 1

    # Meropenem chemical degradation in CAMHB at 37 C.
    # Paper reports 10% loss in 8 h => kdeg = -log(0.9)/8 = 0.01317 1/h.
    # The numeric kdeg itself is not tabulated in Table 1; the rate is
    # implied by the 10%-in-8-h Suppl. Fig. 3 result quoted in the
    # Discussion (page 7).
    lkdeg_mer <- fixed(log(0.01317))
    label("Meropenem chemical degradation rate (kdeg_Mer, 1/h; FIXED, derived from Sadouki 2025 Discussion: 10% loss in 8 h)")  # Sadouki 2025 page 7 + Suppl. Fig. 3 (derived; not in Table 1)

    # ---- Meropenem -- regrowth dynamics ----
    beta_mer <- 0.922
    label("Meropenem max loss of antimicrobial effect (BETA_Mer, unitless)")  # Sadouki 2025 Table 1, Meropenem regrowth
    ltau_mer <- fixed(log(0.570))
    label("Meropenem time-shape of effect loss (TAU_Mer, 1/h; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    coef_taumer_mer <- 0.00155
    label("Proportional Mer concentration effect on TAU_Mer (per mg/L Mer)")  # Sadouki 2025 Table 1, Proportional MER-on-TAU row

    # ---- Gentamicin -- bacterial killing ----
    lemax_gen <- log(5.47)
    label("Gentamicin maximum effect (Emax_Gen, 1/h)")  # Sadouki 2025 Table 1, Gentamicin killing
    lic50_gen <- log(1.12)
    label("Gentamicin half-maximum-effect concentration (IC50_Gen, mg/L)")  # Sadouki 2025 Table 1
    lhill_gen <- log(3.63)
    label("Gentamicin Hill exponent (hill_Gen, unitless)")  # Sadouki 2025 Table 1

    # ---- Gentamicin -- regrowth dynamics ----
    beta_gen <- fixed(0.829)
    label("Gentamicin max loss of antimicrobial effect (BETA_Gen, unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    ltau_gen <- fixed(log(0.517))
    label("Gentamicin time-shape of effect loss (TAU_Gen, 1/h; FIXED)")  # Sadouki 2025 Table 1 (FIXED)

    # ---- Gentamicin -- concentration-dependent shift of BETA_Gen (Emax-on-BETA) ----
    emax_betagen  <- fixed(-2.97)
    label("Max additive shift of Gen concentration on BETA_Gen (unitless; FIXED)")  # Sadouki 2025 Table 1, Emax-on-BETA-Gen rows (FIXED)
    lic50_betagen <- fixed(log(5.72))
    label("IC50 of Gen-concentration effect on BETA_Gen (mg/L; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    lhill_betagen <- fixed(log(20))
    label("Hill exponent of Gen-concentration effect on BETA_Gen (unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)

    # ---- Ciprofloxacin -- bacterial killing ----
    lemax_cip <- log(4.55)
    label("Ciprofloxacin maximum effect (Emax_Cip, 1/h)")  # Sadouki 2025 Table 1, Ciprofloxacin killing
    lic50_cip <- log(0.0106)
    label("Ciprofloxacin half-maximum-effect concentration (IC50_Cip, mg/L)")  # Sadouki 2025 Table 1
    lhill_cip <- log(3.58)
    label("Ciprofloxacin Hill exponent (hill_Cip, unitless)")  # Sadouki 2025 Table 1

    # ---- Ciprofloxacin -- regrowth dynamics ----
    beta_cip <- fixed(0.674)
    label("Ciprofloxacin max loss of antimicrobial effect (BETA_Cip, unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    ltau_cip <- fixed(log(0.359))
    label("Ciprofloxacin time-shape of effect loss (TAU_Cip, 1/h; FIXED)")  # Sadouki 2025 Table 1 (FIXED)

    # ---- Ciprofloxacin -- concentration-dependent shift of BETA_Cip ----
    emax_betacip  <- fixed(-4)
    label("Max additive shift of Cip concentration on BETA_Cip (unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    lic50_betacip <- fixed(log(0.017))
    label("IC50 of Cip-concentration effect on BETA_Cip (mg/L; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    lhill_betacip <- fixed(log(20))
    label("Hill exponent of Cip-concentration effect on BETA_Cip (unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)

    # ---- Drug-drug interactions ----
    combo_beta <- fixed(-1)
    label("Categorical 2- or 3-way combination shift on BETA (additive, unitless; FIXED)")  # Sadouki 2025 Table 1, Drug interactions (FIXED)
    mer_on_ic50cip <- fixed(-0.353)
    label("Proportional Mer-presence effect on IC50_Cip (unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)
    gen_on_ic50cip <- fixed(-0.576)
    label("Proportional Gen-presence effect on IC50_Cip (unitless; FIXED)")  # Sadouki 2025 Table 1 (FIXED)

    # ---- Between-experiment random effects ----
    # Sadouki 2025 Table 1 BSV column. For log-normal parameters BSV is %CV;
    # internal variance is omega^2 = log(1 + (CV/100)^2). For BETA the table
    # footnote says values are additive eta variance estimates (linear scale).
    etalknet     ~ 0.01649  # 12.9% CV: log(1 + 0.129^2) = 0.01649; Sadouki 2025 Table 1 (Growth knet BSV)
    etalemax_mer ~ 0.00294  # 5.43% CV: log(1 + 0.0543^2) = 0.00294; Sadouki 2025 Table 1
    etalic50_mer ~ 0.40988  # 72.3% CV: log(1 + 0.723^2) = 0.40988; Sadouki 2025 Table 1
    etalhill_mer ~ 0.06738  # 26.4% CV: log(1 + 0.264^2) = 0.06738; Sadouki 2025 Table 1
    etabeta_mer  ~ 1.58     # additive eta variance on BETA_Mer (linear scale); Sadouki 2025 Table 1, footnote
    etalemax_gen ~ 0.06035  # 25.5% CV: log(1 + 0.255^2) = 0.06035; Sadouki 2025 Table 1 (Gen Emax BSV)
    etalhill_gen ~ 0.01624  # 12.8% CV: log(1 + 0.128^2) = 0.01624; Sadouki 2025 Table 1
    etalemax_cip ~ 0.12805  # 37.6% CV: log(1 + 0.376^2) = 0.12805; Sadouki 2025 Table 1
    etalic50_cip ~ 0.00733  # 8.58% CV: log(1 + 0.0858^2) = 0.00733; Sadouki 2025 Table 1
    etalhill_cip ~ 0.01451  # 12.1% CV: log(1 + 0.121^2) = 0.01451; Sadouki 2025 Table 1

    # ---- Residual error ----
    # 'Additive on logarithmic data' per Sadouki 2025 Table 1. Implemented as
    # an additive error on log_cfu = log(bacteria) (natural log of CFU/mL).
    addSd_log_cfu <- 0.864
    label("Additive residual SD on log(bacteria) (log CFU/mL)")  # Sadouki 2025 Table 1, Residual variability
  })

  model({
    # ---- Per-experiment effective parameters ----
    knet     <- exp(lknet + etalknet)
    emax_mer <- exp(lemax_mer + etalemax_mer)
    ic50_mer <- exp(lic50_mer + etalic50_mer)
    hill_mer <- exp(lhill_mer + etalhill_mer)
    kdeg_mer <- exp(lkdeg_mer)
    tau_mer  <- exp(ltau_mer)
    emax_gen <- exp(lemax_gen + etalemax_gen)
    ic50_gen <- exp(lic50_gen)
    hill_gen <- exp(lhill_gen + etalhill_gen)
    tau_gen  <- exp(ltau_gen)
    ic50_betagen <- exp(lic50_betagen)
    hill_betagen <- exp(lhill_betagen)
    emax_cip <- exp(lemax_cip + etalemax_cip)
    ic50_cip_base <- exp(lic50_cip + etalic50_cip)
    hill_cip <- exp(lhill_cip + etalhill_cip)
    tau_cip  <- exp(ltau_cip)
    ic50_betacip <- exp(lic50_betacip)
    hill_betacip <- exp(lhill_betacip)

    # ---- Inoculum-adjusted growth parameters ----
    # B0 and Bmax shift by additive categorical effects when LOWINOC = 1
    # (10^3 CFU/mL inoculum). Reference category is 10^5 CFU/mL.
    b0_eff   <- b0   + cat_b0_lowinoc   * LOWINOC
    bmax_eff <- bmax + cat_bmax_lowinoc * LOWINOC

    # ---- Combination-presence indicator (1 if 2- or 3-way regimen) ----
    n_drugs   <- CONMED_MER + CONMED_GEN + CONMED_CIP
    i_combo   <- n_drugs >= 2

    # ---- BETA modifiers per drug ----
    # Mer BETA is shifted by combo only.
    beta_mer_eff <- beta_mer + etabeta_mer + combo_beta * i_combo

    # Gen BETA is shifted by gentamicin concentration (Emax on BETA) and combo.
    beta_gen_eff <- beta_gen +
      emax_betagen * (CONMED_GEN_CC ^ hill_betagen) /
        (ic50_betagen ^ hill_betagen + CONMED_GEN_CC ^ hill_betagen) +
      combo_beta * i_combo

    # Cip BETA is shifted by ciprofloxacin concentration and combo.
    beta_cip_eff <- beta_cip +
      emax_betacip * (CONMED_CIP_CC ^ hill_betacip) /
        (ic50_betacip ^ hill_betacip + CONMED_CIP_CC ^ hill_betacip) +
      combo_beta * i_combo

    # ---- Mer TAU shifted by Mer concentration (proportional) ----
    tau_mer_eff <- tau_mer * (1 + coef_taumer_mer * mer)

    # ---- Cip IC50 shifted by Mer or Gen presence (proportional, multiplicative) ----
    ic50_cip_eff <- ic50_cip_base *
      (1 + mer_on_ic50cip * CONMED_MER + gen_on_ic50cip * CONMED_GEN)

    # ---- Drug effect terms (Hill * (1 - BETA*(1 - exp(-t*TAU)))) ----
    # Each Effect_drug is the per-drug killing rate (1/h); summed below.
    eff_mer <- emax_mer * (mer ^ hill_mer) /
      (ic50_mer ^ hill_mer + mer ^ hill_mer) *
      (1 - beta_mer_eff * (1 - exp(-t * tau_mer_eff)))

    eff_gen <- emax_gen * (CONMED_GEN_CC ^ hill_gen) /
      (ic50_gen ^ hill_gen + CONMED_GEN_CC ^ hill_gen) *
      (1 - beta_gen_eff * (1 - exp(-t * tau_gen)))

    eff_cip <- emax_cip * (CONMED_CIP_CC ^ hill_cip) /
      (ic50_cip_eff ^ hill_cip + CONMED_CIP_CC ^ hill_cip) *
      (1 - beta_cip_eff * (1 - exp(-t * tau_cip)))

    # ---- ODE system ----
    # State `mer` -- meropenem solution concentration (mg/L), first-order chemical
    #                degradation in CAMHB at 37 C.
    # State `bacteria` -- E. coli concentration (CFU/mL on linear scale), logistic
    #                     growth with summed Hill killing terms.
    d/dt(mer)      <- -kdeg_mer * mer
    d/dt(bacteria) <- (knet * (1 - bacteria / (10 ^ bmax_eff)) -
                         (eff_mer + eff_gen + eff_cip)) * bacteria

    # ---- Initial conditions ----
    bacteria(0) <- 10 ^ b0_eff

    # ---- Observation: log(bacteria) (natural log of CFU/mL) with additive error
    log_cfu <- log(bacteria)
    log_cfu ~ add(addSd_log_cfu)
  })
}
