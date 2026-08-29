Mahadevan_2026_polymyxinB_meropenem_fosfomycin_KP00161 <- function() {
  description <- "In vitro (Klebsiella pneumoniae KP0016-1, non-carbapenemase producer with ompK35/ompK36/ompK37 porin mutations; PMB MIC <0.5, MEM MIC 64, FOF MIC 64 mg/L). Mechanism-based PK/PD (Bulitta life-cycle growth) model of static-concentration time-kill by polymyxin B, meropenem and fosfomycin as monotherapy and in double / triple combination against a carbapenem-resistant clinical isolate. Two pre-existing bacterial subpopulations, each split over the two-state replication life cycle: intermediate to all three drugs (the majority population) and resistant to all three (seeded at the estimated mutant frequency). Each drug kills through its own Hill function, and polymyxin B additionally lowers the effective meropenem and fosfomycin KC50 through mechanistic synergy (outer-membrane disruption). This isolate had the lowest estimated meropenem and fosfomycin Kmax of the six. One of six isolate-specific files from this paper"
  reference <- "Mahadevan R, Garcia E, Sharma R, Qiu H, Elsheikh A, Parambi R, Abboud CS, Pasteran F, Ramirez MS, Kaye KS, Bonomo RA, Rao GG. A mechanism-based pharmacokinetic/pharmacodynamic analysis of polymyxin B-based combination therapy against carbapenem-resistant Klebsiella pneumoniae isolates with diverse phenotypic and genotypic resistance mechanisms. Antimicrob Agents Chemother. 2026 Feb;70(2):e00782-25. doi:10.1128/aac.00782-25. PMCID: PMC12888851. Structural model: Materials and Methods, 'Mechanism based PK/PD model development', equations 1-7. Parameter estimates: Table 2, column KP0016-1. Isolate resistance genes and MICs: Table 1 and supplemental Table S2. Model-predicted triple-combination time-kill profiles: Figure 3, row KP0016-1"
  vignette <- "Mahadevan_2026_polymyxinB_combination_CRKP"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  # The "dose" in a static-concentration time-kill (SCTK) assay is the bath
  # concentration placed directly into each drug state at time 0, so the dosing
  # targets are neither `depot` nor `central`.
  dosing <- c("cpmb", "cmem", "cfof")

  # `cmem` is the canonical meropenem bath-concentration state; `cpmb` / `cfof`
  # follow the same `c<drug>` family but are not yet registered canonicals, so
  # only those two are declared here. The bacterial states use the canonical
  # `bact_<pheno>..._<pheno><digit>` N-drug scheme (tokens in PMB, MEM, FOF
  # order; trailing digit is the Bulitta life-cycle state) and need no
  # declaration. See the BRKP61 sibling file for the full note.
  paper_specific_compartments <- c("cpmb", "cfof")

  compartmentData <- list(
    cpmb = list(analyte = "polymyxin B", units = "mg/L", specimen = "administration site", verified = TRUE),
    cmem = list(analyte = "meropenem", units = "mg/L", specimen = "administration site", verified = TRUE),
    cfof = list(analyte = "fosfomycin", units = "mg/L", specimen = "administration site", verified = TRUE),
    bact_intermediate_intermediate_intermediate1 = list(analyte = "Klebsiella pneumoniae KP0016-1, polymyxin B-intermediate / meropenem-intermediate / fosfomycin-intermediate subpopulation, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_intermediate_intermediate_intermediate2 = list(analyte = "Klebsiella pneumoniae KP0016-1, polymyxin B-intermediate / meropenem-intermediate / fosfomycin-intermediate subpopulation, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_resistant1 = list(analyte = "Klebsiella pneumoniae KP0016-1, subpopulation resistant to all three drugs, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_resistant2 = list(analyte = "Klebsiella pneumoniae KP0016-1, subpopulation resistant to all three drugs, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Klebsiella pneumoniae KP0016-1, carbapenem-resistant clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    organism = "Klebsiella pneumoniae KP0016-1, obtained from an individual patient at Siriraj Hospital, Bangkok, Thailand",
    resistance_genes = "no carbapenemase detected (blaKPC, blaGES, blaNDM, blaVIM, blaIMP, blaOXA-40 and blaOXA-48 all not amplified); functional mgrB (comparable in size with wild type); non-functional ompK35 (terminated at 91 amino acids), ompK36 (not amplified) and ompK37 (not amplified) -- three porin mutations (Table 1, Table S2)",
    mic_values = c(
      `polymyxin B` = "<0.5 mg/L (intermediate, CLSI)",
      meropenem = "64 mg/L (resistant, CLSI)",
      fosfomycin = "64 mg/L (resistant, EUCAST E. coli breakpoint)"
    ),
    system = "24 h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton broth (25.0 mg/L Ca2+, 12.5 mg/L Mg2+) supplemented with 25 mg/L glucose-6-phosphate; viable counts at 0, 1, 2, 4, 6, 8 and 24 h",
    temperature = "37 C",
    duration = "24 h",
    starting_inoculum = "~10^6 CFU/mL; the model estimated 10^6.0 CFU/mL for this isolate",
    limit_of_quantification = "2 log10 CFU/mL; observations below it were handled with the Beal M3 method",
    dose_range = "polymyxin B 0.5-64 mg/L, meropenem 10-120 mg/L, fosfomycin 75-500 mg/L, as monotherapy and in double and triple combinations",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Mechanism-based PK/PD model (MBM) fit in S-ADAPT 1.57 via SADAPT-TRAN with the importance-sampling algorithm (pmethod = 4); observations below the LLOQ were fit with the Beal M3 method.",
      "The correlation coefficient between observed and model-predicted log10 CFU/mL was 0.84 for this isolate, the highest of the six (Figure S2).",
      "This is the only non-carbapenemase-producing isolate of the six; it carries three porin mutations.",
      "Mechanistic synergy of polymyxin B on meropenem and fosfomycin was retained for this isolate; it was also retained for BRKP61 and BRKP76 and was not supported for KP0052-1, BRKP67 or BRKP28.",
      "Companion files: Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP61, _BRKP76, _KP00521, _BRKP67, _BRKP28."
    )
  )

  ini({
    # ---- Bacterial growth and subpopulations (Table 2, column KP0016-1) ----
    log10cfu0 <- 6.0; label("Bacterial initial inoculum (log10 CFU/mL)")                                  # Table 2: LogCFU0 = 6.0 (%CV 1.04)
    log10cfumax <- 8.43; label("Maximum population size (log10 CFU/mL)")                                  # Table 2: LogCFUmax = 8.43 (%CV 2.51)
    lk21 <- fixed(log(50)); label("Log replication rate constant, state 2 to state 1 (1/h)")              # Table 2: k21 = 50 (fixed); Methods: "k21 was fixed to 50 h-1"

    mgt_iii <- 98.3; label("Mean generation time, PMB-I / MEM-I / FOF-I subpopulation (min)")             # Table 2: MGT, PMB_I/MEM_I/FOF_I = 98.3 (%CV 5.07)
    mgt_rrr <- 69.26; label("Mean generation time, PMB-R / MEM-R / FOF-R subpopulation (min)")            # Table 2: MGT, PMB_R/MEM_R/FOF_R = 69.26 (%CV 2.31)

    log10mf_rrr <- -4.29; label("Log10 mutant frequency at baseline of the PMB-R / MEM-R / FOF-R subpopulation (unitless fraction)") # Table 2: log10 MF, PMB_R/MEM_R/FOF_R = -4.29 (%CV 3.14)

    # ---- Drug effect of polymyxin B (Table 2) ------------------------------
    kmax_pmb <- 10.9; label("Maximum killing rate constant of polymyxin B (1/h)")                         # Table 2: KmaxPMB = 10.9 (%CV 9.31)
    kc50_pmb_i <- 6.56; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-intermediate subpopulation (mg/L)")  # Table 2: KC50,PMB,I = 6.56 (%CV 13.7)
    kc50_pmb_r <- 61.16; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-resistant subpopulation (mg/L)")    # Table 2: KC50,PMB,R = 61.16 (%CV 7.32)
    hill_pmb <- 1.09; label("Hill coefficient of polymyxin B killing (unitless)")                         # Table 2: Hill coefficient of polymyxin B = 1.09 (%CV 4.14)

    # ---- Drug effect of meropenem (Table 2) --------------------------------
    kmax_mem <- 1.56; label("Maximum killing rate constant of meropenem (1/h)")                           # Table 2: KmaxMEM = 1.56 (%CV 12.9); Results: "the lowest killing rate constant for both meropenem (1.56 h-1)"
    kc50_mem_i <- 112; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-intermediate subpopulation (mg/L)")     # Table 2: KC50,MEM,I = 112 (%CV 9.39)
    kc50_mem_r <- 424; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-resistant subpopulation (mg/L)")        # Table 2: KC50,MEM,R = 424 (%CV 2.83)
    hill_mem <- 1.99; label("Hill coefficient of meropenem killing (unitless)")                           # Table 2: Hill coefficient of meropenem = 1.99 (%CV 9.72)

    # ---- Drug effect of fosfomycin (Table 2) -------------------------------
    kmax_fof <- 0.64; label("Maximum killing rate constant of fosfomycin (1/h)")                          # Table 2: KmaxFOF = 0.64 (%CV 8.00); Results: "and fosfomycin (0.64 h-1)"
    kc50_fof_i <- 20.2; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-intermediate subpopulation (mg/L)")   # Table 2: KC50,FOF,I = 20.2 (%CV 7.54)
    kc50_fof_r <- 448; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-resistant subpopulation (mg/L)")       # Table 2: KC50,FOF,R = 448 (%CV 7.89)
    hill_fof <- 0.79; label("Hill coefficient of fosfomycin killing (unitless)")                          # Table 2: Hill coefficient of fosfomycin = 0.79 (%CV 22.5)

    # ---- Mechanistic synergy of polymyxin B on meropenem (Table 2) ---------
    imax_mem_pmb <- 0.88; label("Maximum fractional decrease of the meropenem KC50 by polymyxin B via outer-membrane disruption (unitless)") # Table 2: Imax,M,PMB = 0.88 (%CV 18.7); Results: "Imax,M,PMB values of 0.84, 0.83, and 0.88 for BRKP61, BRKP76, and KP0016-1"
    ic50_mem_pmb <- 0.50; label("Polymyxin B concentration causing 50% of Imax,M,PMB (mg/L)")             # Table 2: IC50,M,PMB = 0.50 (%CV 26.65)
    hill_syn_mem <- 0.63; label("Hill coefficient for the mechanistic synergy on meropenem (unitless)")   # Table 2: Hill coefficient for the mechanistic synergy (meropenem block) = 0.63 (%CV 14.7)

    # ---- Mechanistic synergy of polymyxin B on fosfomycin (Table 2) --------
    imax_fof_pmb <- 0.89; label("Maximum fractional decrease of the fosfomycin KC50 by polymyxin B via outer-membrane disruption (unitless)") # Table 2: Imax,F,PMB = 0.89 (%CV 13.6); Results: "the Imax,F,PMB values were 0.99, 0.81, and 0.89"
    ic50_fof_pmb <- 0.59; label("Polymyxin B concentration causing 50% of Imax,F,PMB (mg/L)")             # Table 2: IC50,F,PMB = 0.59 (%CV 15.7)
    hill_syn_fof <- 0.56; label("Hill coefficient for the mechanistic synergy on fosfomycin (unitless)")  # Table 2: Hill coefficient for the mechanistic synergy (fosfomycin block) = 0.56 (%CV 12.3)

    # ---- Residual variability (Table 2, Variability model) -----------------
    addSd <- 0.45; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)")              # Table 2: Additive residual variability = 0.45 (%CV 5.67)
  })

  model({
    # ---- 1. Back-transforms and derived rate constants ---------------------
    k21 <- exp(lk21)
    cfu0 <- 10^log10cfu0
    cfumax <- 10^log10cfumax
    k12_iii <- 60 / mgt_iii
    k12_rrr <- 60 / mgt_rrr
    mf_rrr <- 10^log10mf_rrr
    mf_iii <- 1 - mf_rrr

    # ---- 2. Total viable count and the replication plateau (equation 1) ----
    CFUtot <- bact_intermediate_intermediate_intermediate1 +
      bact_intermediate_intermediate_intermediate2 +
      bact_resistant_resistant_resistant1 +
      bact_resistant_resistant_resistant2
    plat <- 1 - CFUtot / (cfumax + CFUtot)

    # ---- 3. Mechanistic synergy of polymyxin B (equation 5) ----------------
    syn_mem <- 1 - imax_mem_pmb * cpmb^hill_syn_mem /
      (cpmb^hill_syn_mem + ic50_mem_pmb^hill_syn_mem)
    syn_fof <- 1 - imax_fof_pmb * cpmb^hill_syn_fof /
      (cpmb^hill_syn_fof + ic50_fof_pmb^hill_syn_fof)

    # ---- 4. Per-subpopulation killing (equations 4, 6, 7) ------------------
    kill_iii <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_i^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + (kc50_mem_i * syn_mem)^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + (kc50_fof_i * syn_fof)^hill_fof)
    kill_rrr <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_r^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + (kc50_mem_r * syn_mem)^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + (kc50_fof_r * syn_fof)^hill_fof)

    # ---- 5. Bacterial life-cycle system (equations 2 and 3) ----------------
    d/dt(bact_intermediate_intermediate_intermediate1) <-
      2 * plat * k21 * bact_intermediate_intermediate_intermediate2 -
      k12_iii * bact_intermediate_intermediate_intermediate1 -
      kill_iii * bact_intermediate_intermediate_intermediate1
    d/dt(bact_intermediate_intermediate_intermediate2) <-
      -k21 * bact_intermediate_intermediate_intermediate2 +
      k12_iii * bact_intermediate_intermediate_intermediate1 -
      kill_iii * bact_intermediate_intermediate_intermediate2

    d/dt(bact_resistant_resistant_resistant1) <-
      2 * plat * k21 * bact_resistant_resistant_resistant2 -
      k12_rrr * bact_resistant_resistant_resistant1 -
      kill_rrr * bact_resistant_resistant_resistant1
    d/dt(bact_resistant_resistant_resistant2) <-
      -k21 * bact_resistant_resistant_resistant2 +
      k12_rrr * bact_resistant_resistant_resistant1 -
      kill_rrr * bact_resistant_resistant_resistant2

    # ---- 6. Static drug exposures ------------------------------------------
    d/dt(cpmb) <- 0
    d/dt(cmem) <- 0
    d/dt(cfof) <- 0

    # ---- 7. Initial conditions (equations 2 and 3) --------------------------
    bact_intermediate_intermediate_intermediate1(0) <- cfu0 * mf_iii
    bact_resistant_resistant_resistant1(0) <- cfu0 * mf_rrr

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
