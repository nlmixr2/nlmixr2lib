Lacroix_2025_polymyxinB_AB121D0 <- function() {
  description <- "In vitro (Acinetobacter baumannii AB121-D0; multidrug-resistant clinical isolate recovered from a patient BEFORE colistin treatment; polymyxin B MIC 1 mg/L in CAMHB). Semi-mechanistic time-kill PK/PD model of polymyxin B (PMB) in cation-adjusted Mueller-Hinton broth with or without 1% porcine-stomach mucin. Two bacterial subpopulations (bact_s susceptible, bact_r resistant) share one logistic net growth rate; PMB killing follows a non-sigmoidal Emax function of the UNBOUND PMB concentration for both subpopulations. Mucin acts twice: (i) it binds PMB, so the unbound fraction fu follows a saturable sigmoidal Emax function of the TOTAL PMB concentration (Eq 1), and (ii) beyond binding it increases the resistant-subpopulation Emax 1.4-fold (categorical mucin covariate on Emax_R, dOFV = 152). Companion model: Lacroix_2025_polymyxinB_AB122D12 (post-colistin isolate; sigmoidal Emax, mucin acting on Knet and EC50 instead)."
  reference <- paste(
    "Lacroix M, Moreau J, Zampaloni C, Bissantz C, Mirfendereski H, Shirvani H,",
    "Marchand S, Couet W, Chauzy A. (2025).",
    "In vitro pharmacokinetic/pharmacodynamic modeling of the effect of mucin on",
    "polymyxin B activity against Acinetobacter baumannii.",
    "Antimicrobial Agents and Chemotherapy 69(5):e01535-24.",
    "doi:10.1128/aac.01535-24.",
    "Structural model: supplemental text (AAC01535-24-S0002) Eqs S1-S6 and",
    "supplemental Fig. S2 (left panel, AB121-D0), which prints the final kill-rate",
    "and covariate equations. Unbound-fraction model: main-text Eq 1.",
    "Parameter estimates: main-text Table 1 (fu model) and Table 3, columns",
    "'AB121-D0 without mucin' and 'AB121-D0 with mucin'.",
    sep = " "
  )
  vignette <- "Lacroix_2025_polymyxinB_mucin"
  units <- list(
    time          = "h",
    dosing        = "none (PMB is held at a static total concentration supplied by the CONC_PMB_MGL covariate; the TK tubes receive no dosing events)",
    concentration = "log10 CFU/mL (Cc, total viable count of bact_s + bact_r)"
  )

  # Bacterial subpopulations of the hetero-resistance (S/R) model. Matches the
  # registered "^bact_" paper-specific compartment pattern and the sibling
  # A. baumannii / polymyxin models (Cheah_2016_polymyxin_*).
  paper_specific_compartments <- c("bact_s", "bact_r")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Bacterial subpopulations are PD states, not a
  # biological matrix, so specimen is "not applicable" (see conventions.R).
  compartmentData <- list(
    bact_s = list(analyte = "Acinetobacter baumannii AB121-D0 polymyxin-B-susceptible subpopulation", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r = list(analyte = "Acinetobacter baumannii AB121-D0 polymyxin-B-resistant subpopulation", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    CONC_PMB_MGL = list(
      description = "Static TOTAL polymyxin B concentration applied to the time-kill tube. Drives the unbound-fraction model (Eq 1) and hence the unbound concentration Cu = CONC_PMB_MGL * fu that produces the killing effect.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- AB121-D0 was studied at total PMB concentrations of 0.5 to 128 mg/L (Materials and Methods, TK experiments); 0 mg/L is the drug-free growth control.",
      notes = "PMB was assumed stable over the 30-h experiment (Materials and Methods, Pharmacodynamic modeling), so the concentration is time-invariant and carried as a covariate rather than an ODE state. In CAMHB without mucin the unbound and total concentrations coincide (fu = 1).",
      source_name = "C_T"
    ),
    MUCIN_PRESENT = list(
      description = "Binary indicator of 1% (w/v) porcine-stomach mucin supplementation of the cation-adjusted Mueller-Hinton broth used for the time-kill experiment. 1 = CAMHB + 1% mucin; 0 = CAMHB alone (control).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = CAMHB without mucin",
      notes = "Enters the model twice. (i) Binding: when MUCIN_PRESENT = 1 the unbound fraction follows the sigmoidal Emax function of Eq 1; when 0 there is no mucin to bind PMB and fu is 1. (ii) Residual activity effect: main-text Eq 4 selects Emax_R = emax_r_camhb when 0 and emax_r_mucin when 1 (supplemental Fig. S2, left panel).",
      source_name = "Mucin"
    )
  )

  population <- list(
    species = "in vitro (Acinetobacter baumannii AB121-D0)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    organism = "A. baumannii AB121-D0, a multidrug-resistant clinical isolate recovered from a patient before colistin treatment (the pre-treatment member of the AB121-D0 / AB122-D12 isogenic-pair). Polymyxin B MIC 1 mg/L in CAMHB and 64 mg/L in CAMHB + 1% mucin, i.e. an unbound MIC (MICu) of 1 vs 13 mg/L, a 13-fold MICu increase (Table 2).",
    model_system = "Static time-kill (TK) experiments in glass tubes, 35 degrees C with shaking at 150-170 rpm, in commercial CAMHB with or without 1% (w/v) porcine-stomach mucin. Starting inoculum approximately 1e6 CFU/mL. Viable counts sampled at 0, 1, 3, 6, 10, 24 and 30 h, plated on Mueller-Hinton agar supplemented with 1% activated charcoal to prevent PMB carry-over. Experiments in triplicate.",
    dose_range = "Total PMB 0.5 to 128 mg/L (plus a drug-free growth control); with 1% mucin these correspond to unbound PMB concentrations of about 0.03 to 13 mg/L after correction by fu.",
    lloq = "400 CFU/mL, i.e. 2.6 log10 CFU/mL (CLSI); observations below it were handled by Beal's M3 method during estimation.",
    notes = paste(
      "Estimated in NONMEM 7.4.2 with the LAPLACIAN algorithm on log10-transformed bacterial counts;",
      "uncertainty by sampling importance resampling (SIR).",
      "The mucin-free and mucin-containing TK data were fit simultaneously in one data set.",
      "The only between-experiment random effects in the published model are the FIXED-variance additive etas on the four unbound-fraction parameters,",
      "which propagate the uncertainty of the Eq 1 binding model into each TK experiment (Materials and Methods, Eqs 2-3).",
      "The PD parameters themselves carry no IIV.",
      "Whole-genome sequencing of regrowth isolates found pmrABC mutations in 9 of 11 samples with a raised MIC, but mutations appeared randomly across replicates and only partially explained regrowth."
    )
  )

  ini({
    # ===============================================================
    # PMB binding to mucin -- sigmoidal Emax unbound-fraction model
    # (main-text Eq 1, parameters from Table 1, pooled data set)
    # ===============================================================
    fu_e0 <- 0.061
    label("Minimum PMB unbound fraction in CAMHB + 1% mucin (unitless; E0)")  # Lacroix 2025 Table 1
    fu_imax <- 0.56
    label("Maximum increase in PMB unbound fraction above E0 (unitless; Imax)")  # Lacroix 2025 Table 1
    fu_ic50 <- 114.37
    label("Total PMB concentration giving 50% of the maximum fu increase (mg/L; IC50)")  # Lacroix 2025 Table 1
    fu_hill <- 1.97
    label("Hill coefficient of the unbound-fraction model (unitless; epsilon)")  # Lacroix 2025 Table 1

    # Between-TK-experiment additive random effects on the Eq 1 parameters.
    # Variances are FIXED to the estimation variances of Table 1 so that the
    # uncertainty in fu propagates into each simulated experiment
    # (Materials and Methods, Eqs 2-3: "the variance omega^2_E0 fixed to the
    # value obtained for each parameter estimate of equation 1 (Table 1)").
    etafu_e0 ~ fixed(0.0004)    # Lacroix 2025 Table 1, variance column for E0
    etafu_imax ~ fixed(0.0009)  # Lacroix 2025 Table 1, variance column for Imax
    etafu_ic50 ~ fixed(196)     # Lacroix 2025 Table 1, variance column for IC50
    etafu_hill ~ fixed(0.252)   # Lacroix 2025 Table 1, variance column for epsilon

    # ===============================================================
    # Bacterial growth sub-model (Eqs S1-S2)
    # ===============================================================
    log10_inoc <- 6.0
    label("log10 initial total inoculum (log10 CFU/mL; INOC)")  # Lacroix 2025 Table 3, AB121-D0
    log10_cfumax <- 8.47
    label("log10 maximum bacterial population size / carrying capacity (log10 CFU/mL; BMAX)")  # Lacroix 2025 Table 3, AB121-D0
    log10_fmut <- -4.72
    label("log10 proportion of resistant bacteria in the initial inoculum (log10; MUT)")  # Lacroix 2025 Table 3, AB121-D0
    knet <- 2.83
    label("Apparent net growth rate constant, shared by both subpopulations (1/h; Knet)")  # Lacroix 2025 Table 3, AB121-D0 (no mucin effect for this strain)

    # ===============================================================
    # PMB killing sub-model -- Emax (non-sigmoidal), Eqs S5-S6
    # ===============================================================
    emax_s <- 21.1
    label("Maximum PMB kill rate constant, susceptible subpopulation (1/h; Emax_S)")  # Lacroix 2025 Table 3, AB121-D0
    emax_r_camhb <- 2.84
    label("Maximum PMB kill rate constant, resistant subpopulation, without mucin (1/h; Emax_R)")  # Lacroix 2025 Table 3, AB121-D0 without mucin
    emax_r_mucin <- 3.90
    label("Maximum PMB kill rate constant, resistant subpopulation, with 1% mucin (1/h; Emax_R)")  # Lacroix 2025 Table 3, AB121-D0 with mucin
    ec50 <- 0.40
    label("Unbound PMB concentration giving 50% of Emax, both subpopulations (mg/L; EC50)")  # Lacroix 2025 Table 3, AB121-D0

    # ===============================================================
    # Residual error
    # ===============================================================
    addSd <- 0.35
    label("Additive residual SD on the log10 total bacterial count (log10 CFU/mL; sigma)")  # Lacroix 2025 Table 3, AB121-D0
  })

  model({
    # --- Unbound PMB concentration -------------------------------------
    # Additive between-experiment deviations on the Eq 1 parameters (Eq 3).
    e0i <- fu_e0 + etafu_e0
    imaxi <- fu_imax + etafu_imax
    ic50i <- fu_ic50 + etafu_ic50
    epsi <- fu_hill + etafu_hill

    ct <- CONC_PMB_MGL
    ct_e <- ct^epsi
    ic50_e <- ic50i^epsi
    # Eq 2: fu_i = E0_i + Imax_i * C^eps_i / (IC50_i^eps_i + C^eps_i).
    fu_mucin <- e0i + imaxi * ct_e / ((ic50_e) + (ct_e))
    # Without mucin there is nothing to bind PMB, so the unbound and total
    # concentrations coincide (Results, "MICs in the absence of mucin
    # corresponded directly to MICu values").
    fu <- (1 - MUCIN_PRESENT) + (fu_mucin * MUCIN_PRESENT)
    cu <- ct * fu

    # --- Mucin effect on the PD parameters (main-text Eq 4) -------------
    # For AB121-D0 mucin acts on Emax_R only (supplemental Fig. S2, left panel).
    emax_r <- (emax_r_camhb * (1 - MUCIN_PRESENT)) + (emax_r_mucin * MUCIN_PRESENT)

    # --- Kill rates (Eqs S5-S6): Emax, no Hill coefficient for this strain
    kill_s <- emax_s * cu / ((ec50) + (cu))
    kill_r <- emax_r * cu / ((ec50) + (cu))

    # --- Bacterial dynamics (Eqs S3-S4) --------------------------------
    cfumax <- 10^log10_cfumax
    btot <- bact_s + bact_r
    logistic <- 1 - btot / cfumax

    d/dt(bact_s) <- knet * logistic * bact_s - kill_s * bact_s
    d/dt(bact_r) <- knet * logistic * bact_r - kill_r * bact_r

    # Initial conditions: MUT is the log10 PROPORTION of the inoculum that is
    # resistant, so the inoculum is split between the two subpopulations.
    inoc0 <- 10^log10_inoc
    fmut <- 10^log10_fmut
    bact_s(0) <- inoc0 * (1 - fmut)
    bact_r(0) <- inoc0 * fmut

    # --- Observation ----------------------------------------------------
    # +1 CFU/mL floors the count so log10 stays finite when PMB sterilises the
    # tube; 1 CFU/mL is 2.6 log10 below the 400 CFU/mL limit of quantification.
    cfu_obs <- bact_s + bact_r + 1
    Cc <- log10(cfu_obs)
    Cc ~ add(addSd)
  })
}
