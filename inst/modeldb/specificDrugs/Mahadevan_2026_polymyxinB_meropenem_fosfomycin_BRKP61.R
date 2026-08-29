Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP61 <- function() {
  description <- "In vitro (Klebsiella pneumoniae BRKP61, blaKPC-2 + ompK35/ompK36 porin mutations; PMB MIC <0.5, MEM MIC 128, FOF MIC 256 mg/L). Mechanism-based PK/PD (Bulitta life-cycle growth) model of static-concentration time-kill by polymyxin B, meropenem and fosfomycin as monotherapy and in double / triple combination against a carbapenem-resistant clinical isolate. Two pre-existing bacterial subpopulations, each split over the two-state replication life cycle: PMB-intermediate / MEM-resistant / FOF-resistant (the majority population) and PMB-resistant / MEM-intermediate / FOF-intermediate (seeded at the estimated mutant frequency). Each drug kills through its own Hill function, and polymyxin B additionally lowers the effective meropenem and fosfomycin KC50 through mechanistic synergy (outer-membrane disruption). One of six isolate-specific files from this paper"
  reference <- "Mahadevan R, Garcia E, Sharma R, Qiu H, Elsheikh A, Parambi R, Abboud CS, Pasteran F, Ramirez MS, Kaye KS, Bonomo RA, Rao GG. A mechanism-based pharmacokinetic/pharmacodynamic analysis of polymyxin B-based combination therapy against carbapenem-resistant Klebsiella pneumoniae isolates with diverse phenotypic and genotypic resistance mechanisms. Antimicrob Agents Chemother. 2026 Feb;70(2):e00782-25. doi:10.1128/aac.00782-25. PMCID: PMC12888851. Structural model: Materials and Methods, 'Mechanism based PK/PD model development', equations 1-7. Parameter estimates: Table 2, column BRKP61. Isolate resistance genes and MICs: Table 1 and supplemental Table S2. Model-predicted triple-combination time-kill profiles: Figure 3, row BRKP61"
  vignette <- "Mahadevan_2026_polymyxinB_combination_CRKP"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # The "dose" in a static-concentration time-kill (SCTK) assay is the bath
  # concentration placed directly into each drug state at time 0, so the dosing
  # targets are neither `depot` nor `central`.
  dosing <- c("cpmb", "cmem", "cfof")

  # `cmem` is the canonical meropenem bath-concentration state. `cpmb` and
  # `cfof` follow the same documented `c<drug>` bath-state family (cmem, ccip,
  # ctob) but are not yet registered canonicals, so they are declared here --
  # the same route Kroemer_2024_ceftazidime_avibactam_fosfomycin_hfim.R takes
  # for its fosfomycin state.
  #
  # The bacterial states use the canonical N-drug subpopulation scheme
  # `bact_<pheno>..._<pheno><digit>` (compartment-names.md, "Bacterial
  # subpopulation states"): one spelled-out phenotype token per drug in the
  # model's drug order -- here PMB, MEM, FOF -- and a trailing digit for the
  # Bulitta / Wicha life-cycle state (1 = vegetative, 2 = replicating). These
  # six files are the founding three-drug example of that scheme;
  # `bacterialSubpopRegex` in R/conventions.R was widened from at-most-two
  # tokens to any number in the same commit, so the names are canonical and
  # need no paper-specific declaration.
  paper_specific_compartments <- c("cpmb", "cfof")

  compartmentData <- list(
    cpmb = list(analyte = "polymyxin B", units = "mg/L", specimen = "administration site", verified = TRUE),
    cmem = list(analyte = "meropenem", units = "mg/L", specimen = "administration site", verified = TRUE),
    cfof = list(analyte = "fosfomycin", units = "mg/L", specimen = "administration site", verified = TRUE),
    bact_intermediate_resistant_resistant1 = list(analyte = "Klebsiella pneumoniae BRKP61, polymyxin B-intermediate / meropenem-resistant / fosfomycin-resistant subpopulation, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_intermediate_resistant_resistant2 = list(analyte = "Klebsiella pneumoniae BRKP61, polymyxin B-intermediate / meropenem-resistant / fosfomycin-resistant subpopulation, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_intermediate_intermediate1 = list(analyte = "Klebsiella pneumoniae BRKP61, polymyxin B-resistant / meropenem-intermediate / fosfomycin-intermediate subpopulation, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_intermediate_intermediate2 = list(analyte = "Klebsiella pneumoniae BRKP61, polymyxin B-resistant / meropenem-intermediate / fosfomycin-intermediate subpopulation, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Klebsiella pneumoniae BRKP61, carbapenem-resistant clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    organism = "Klebsiella pneumoniae BRKP61, obtained from an individual patient at Instituto Dante Pazzanese de Cardiologia, Sao Paulo, Brazil",
    resistance_genes = "blaKPC-2; functional mgrB; non-functional ompK35 (not amplified) and ompK36 (stop codon at 2nd amino acid); ompK37 comparable in size with wild type (Table 1, Table S2)",
    mic_values = c(
      `polymyxin B` = "<0.5 mg/L (intermediate, CLSI)",
      meropenem = "128 mg/L (resistant, CLSI)",
      fosfomycin = "256 mg/L (resistant, EUCAST E. coli breakpoint)"
    ),
    system = "24 h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton broth (25.0 mg/L Ca2+, 12.5 mg/L Mg2+) supplemented with 25 mg/L glucose-6-phosphate; viable counts at 0, 1, 2, 4, 6, 8 and 24 h",
    temperature = "37 C",
    duration = "24 h",
    starting_inoculum = "~10^6 CFU/mL; the model estimated 10^6.11 CFU/mL for this isolate",
    limit_of_quantification = "2 log10 CFU/mL; observations below it were handled with the Beal M3 method",
    dose_range = "polymyxin B 0.5-64 mg/L, meropenem 10-120 mg/L, fosfomycin 75-500 mg/L, as monotherapy and in double and triple combinations",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Mechanism-based PK/PD model (MBM) fit in S-ADAPT 1.57 via SADAPT-TRAN with the importance-sampling algorithm (pmethod = 4); observations below the LLOQ were fit with the Beal M3 method.",
      "The correlation coefficient between observed and model-predicted log10 CFU/mL was 0.79 for this isolate (Figure S2).",
      "Mechanistic synergy of polymyxin B on meropenem and fosfomycin was retained for this isolate; it was also retained for BRKP76 and KP0016-1 and was not supported for KP0052-1, BRKP67 or BRKP28.",
      "Companion files: Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP76, _KP00161, _KP00521, _BRKP67, _BRKP28."
    )
  )

  ini({
    # ---- Bacterial growth and subpopulations (Table 2, column BRKP61) -------
    log10cfu0 <- 6.11; label("Bacterial initial inoculum (log10 CFU/mL)")                                 # Table 2: LogCFU0 = 6.11 (%CV 0.64)
    log10cfumax <- 9.33; label("Maximum population size (log10 CFU/mL)")                                  # Table 2: LogCFUmax = 9.33 (%CV 1.87)
    lk21 <- fixed(log(50)); label("Log replication rate constant, state 2 to state 1 (1/h)")              # Table 2: k21 = 50 (fixed); Methods: "k21 was fixed to 50 h-1"

    mgt_irr <- 83.3; label("Mean generation time, PMB-I / MEM-R / FOF-R subpopulation (min)")             # Table 2: MGT, PMB_I/MEM_R/FOF_R = 83.3 (%CV 3.73)
    mgt_rii <- 13.5; label("Mean generation time, PMB-R / MEM-I / FOF-I subpopulation (min)")             # Table 2: MGT, PMB_R/MEM_I/FOF_I = 13.5 (%CV 3.73)

    log10mf_rii <- -5.10; label("Log10 mutant frequency at baseline of the PMB-R / MEM-I / FOF-I subpopulation (unitless fraction)") # Table 2: log10 MF, PMB_R/MEM_I/FOF_I = -5.10 (%CV 1.91)

    # ---- Drug effect of polymyxin B (Table 2) ------------------------------
    kmax_pmb <- 13.42; label("Maximum killing rate constant of polymyxin B (1/h)")                        # Results text: "other isolates (6.61-13.42 h-1)", BRKP61 being the maximum; Table 2 rounds the same estimate to 13.4 (%CV 2.76)
    kc50_pmb_i <- 0.82; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-intermediate subpopulation (mg/L)")  # Table 2: KC50,PMB,I = 0.82 (%CV 3.53)
    kc50_pmb_r <- 41.03; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-resistant subpopulation (mg/L)")    # Table 2: KC50,PMB,R = 41.03 (%CV 5.83)
    hill_pmb <- 0.64; label("Hill coefficient of polymyxin B killing (unitless)")                         # Table 2: Hill coefficient of polymyxin B = 0.64 (%CV 3.83)

    # ---- Drug effect of meropenem (Table 2) --------------------------------
    kmax_mem <- 2.13; label("Maximum killing rate constant of meropenem (1/h)")                           # Table 2: KmaxMEM = 2.13 (%CV 4.39)
    kc50_mem_i <- 24.7; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-intermediate subpopulation (mg/L)")    # Table 2: KC50,MEM,I = 24.7 (%CV 6.50)
    kc50_mem_r <- 61.4; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-resistant subpopulation (mg/L)")       # Table 2: KC50,MEM,R = 61.4 (%CV 4.30)
    hill_mem <- 1.85; label("Hill coefficient of meropenem killing (unitless)")                           # Table 2: Hill coefficient of meropenem = 1.85 (%CV 6.54)

    # ---- Drug effect of fosfomycin (Table 2) -------------------------------
    kmax_fof <- 3.74; label("Maximum killing rate constant of fosfomycin (1/h)")                          # Table 2: KmaxFOF = 3.74 (%CV 7.62)
    kc50_fof_i <- 42.5; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-intermediate subpopulation (mg/L)")   # Table 2: KC50,FOF,I = 42.5 (%CV 7.57)
    kc50_fof_r <- 44.71; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-resistant subpopulation (mg/L)")     # Table 2: KC50,FOF,R = 44.71 (%CV 4.26)
    hill_fof <- 0.61; label("Hill coefficient of fosfomycin killing (unitless)")                          # Table 2: Hill coefficient of fosfomycin = 0.61 (%CV 9.59)

    # ---- Mechanistic synergy of polymyxin B on meropenem (Table 2) ---------
    imax_mem_pmb <- 0.84; label("Maximum fractional decrease of the meropenem KC50 by polymyxin B via outer-membrane disruption (unitless)") # Table 2: Imax,M,PMB = 0.84 (%CV 7.32); Results: "Imax,M,PMB values of 0.84, 0.83, and 0.88 for BRKP61, BRKP76, and KP0016-1"
    ic50_mem_pmb <- 0.51; label("Polymyxin B concentration causing 50% of Imax,M,PMB (mg/L)")             # Table 2: IC50,M,PMB = 0.51 (%CV 6.92)
    hill_syn_mem <- 0.88; label("Hill coefficient for the mechanistic synergy on meropenem (unitless)")   # Table 2: Hill coefficient for the mechanistic synergy (meropenem block) = 0.88 (%CV 6.50)

    # ---- Mechanistic synergy of polymyxin B on fosfomycin (Table 2) --------
    imax_fof_pmb <- 0.99; label("Maximum fractional decrease of the fosfomycin KC50 by polymyxin B via outer-membrane disruption (unitless)") # Table 2: Imax,F,PMB = 0.99 (%CV 9.89); Results: "the Imax,F,PMB values were 0.99, 0.81, and 0.89". The Abstract's "81%-98%" range appears to round the same estimate down; see the vignette Errata
    ic50_fof_pmb <- 0.64; label("Polymyxin B concentration causing 50% of Imax,F,PMB (mg/L)")             # Table 2: IC50,F,PMB = 0.64 (%CV 9.35)
    hill_syn_fof <- 0.72; label("Hill coefficient for the mechanistic synergy on fosfomycin (unitless)")  # Table 2: Hill coefficient for the mechanistic synergy (fosfomycin block) = 0.72 (%CV 7.22)

    # ---- Residual variability (Table 2, Variability model) -----------------
    addSd <- 0.39; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)")              # Table 2: Additive residual variability = 0.39 (%CV 7.07)
  })

  model({
    # ---- 1. Back-transforms and derived rate constants ---------------------
    k21 <- exp(lk21)
    cfu0 <- 10^log10cfu0
    cfumax <- 10^log10cfumax
    # Table 2 reports the mean generation time in minutes; k12 is its inverse
    # (Methods: "k12,IRR, the inverse of mean replication time from state 1 to
    # state 2"), so 60 converts min to 1/h.
    k12_irr <- 60 / mgt_irr
    k12_rii <- 60 / mgt_rii
    # The mutant subpopulation carries fraction 10^log10MF of the inoculum; the
    # majority subpopulation carries the remainder, so the two initial
    # conditions sum to CFU0 (equation 2: IC_CFU_IRR,1 = CFU0 * MF_IRR).
    mf_rii <- 10^log10mf_rii
    mf_irr <- 1 - mf_rii

    # ---- 2. Total viable count and the replication plateau (equation 1) ----
    CFUtot <- bact_intermediate_resistant_resistant1 +
      bact_intermediate_resistant_resistant2 +
      bact_resistant_intermediate_intermediate1 +
      bact_resistant_intermediate_intermediate2
    # REP = 2 * (1 - CFUtot / (CFUMAX + CFUtot)) (equation 3 footnote). At
    # CFUtot = CFUMAX the plateau factor is 0.5, so REP = 1 and replication
    # exactly replaces the parent cell: the population settles at CFUMAX.
    plat <- 1 - CFUtot / (cfumax + CFUtot)

    # ---- 3. Mechanistic synergy of polymyxin B (equation 5) ----------------
    # Mechanistic_synergy = 1 - IMAX * CPMB^h / (CPMB^h + IC50^h). It multiplies
    # the meropenem and fosfomycin KC50 inside the Hill denominator (equations 6
    # and 7). The Methods sentence introducing equations 6-7 lists all four
    # KC50 terms (KC50,MEM,I, KC50,FOF,I, KC50,MEM,R, KC50,FOF,R) as enhanced,
    # and the Table 2 row labels read "KC50,MEM,I and/or KC50,MEM,R", so the
    # synergy is applied to BOTH subpopulations for both drugs.
    syn_mem <- 1 - imax_mem_pmb * cpmb^hill_syn_mem /
      (cpmb^hill_syn_mem + ic50_mem_pmb^hill_syn_mem)
    syn_fof <- 1 - imax_fof_pmb * cpmb^hill_syn_fof /
      (cpmb^hill_syn_fof + ic50_fof_pmb^hill_syn_fof)

    # ---- 4. Per-subpopulation killing (equations 4, 6, 7) ------------------
    # Subpopulation 1 is PMB-intermediate / MEM-resistant / FOF-resistant, so it
    # takes KC50,PMB,I, KC50,MEM,R and KC50,FOF,R; subpopulation 2 is the
    # mirror image (Figure 2B and the Table 2 row labels for BRKP61).
    kill_irr <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_i^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + (kc50_mem_r * syn_mem)^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + (kc50_fof_r * syn_fof)^hill_fof)
    kill_rii <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_r^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + (kc50_mem_i * syn_mem)^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + (kc50_fof_i * syn_fof)^hill_fof)

    # ---- 5. Bacterial life-cycle system (equations 2 and 3) ----------------
    d/dt(bact_intermediate_resistant_resistant1) <-
      2 * plat * k21 * bact_intermediate_resistant_resistant2 -
      k12_irr * bact_intermediate_resistant_resistant1 -
      kill_irr * bact_intermediate_resistant_resistant1
    d/dt(bact_intermediate_resistant_resistant2) <-
      -k21 * bact_intermediate_resistant_resistant2 +
      k12_irr * bact_intermediate_resistant_resistant1 -
      kill_irr * bact_intermediate_resistant_resistant2

    d/dt(bact_resistant_intermediate_intermediate1) <-
      2 * plat * k21 * bact_resistant_intermediate_intermediate2 -
      k12_rii * bact_resistant_intermediate_intermediate1 -
      kill_rii * bact_resistant_intermediate_intermediate1
    d/dt(bact_resistant_intermediate_intermediate2) <-
      -k21 * bact_resistant_intermediate_intermediate2 +
      k12_rii * bact_resistant_intermediate_intermediate1 -
      kill_rii * bact_resistant_intermediate_intermediate2

    # ---- 6. Static drug exposures ------------------------------------------
    # Static-concentration time-kill: the bath concentrations are constant. A
    # user simulating a clinical regimen overrides these with elimination terms
    # or drives them from an external population PK model.
    d/dt(cpmb) <- 0
    d/dt(cmem) <- 0
    d/dt(cfof) <- 0

    # ---- 7. Initial conditions (equations 2 and 3) --------------------------
    # All bacteria start in the vegetative state 1 (IC_CFU,2 = 0).
    bact_intermediate_resistant_resistant1(0) <- cfu0 * mf_irr
    bact_resistant_intermediate_intermediate1(0) <- cfu0 * mf_rii

    # ---- 8. Observation -----------------------------------------------------
    # The 1e-6 CFU/mL regularisation keeps Cc finite when a regimen drives the
    # population to numerical extinction (the solver can undershoot to a small
    # negative value, and log10 of that is NaN). It is 8 orders of magnitude
    # below the 2 log10 CFU/mL LLOQ, so it never perturbs an observable count.
    # Same device as HernandezLozano_2025_apramycin_invitro.R.
    Cc <- log10(CFUtot + 1e-6)
    Cc ~ add(addSd)
  })
}
