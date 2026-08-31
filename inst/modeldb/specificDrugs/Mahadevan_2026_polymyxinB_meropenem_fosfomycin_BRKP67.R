Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP67 <- function() {
  description <- "In vitro (Klebsiella pneumoniae BRKP67, blaKPC-2 + non-functional MgrB (ISKpn13 insertion) + ompK35/ompK36 porin mutations; PMB MIC 8, MEM MIC 64, FOF MIC 32 mg/L). Mechanism-based PK/PD (Bulitta life-cycle growth) model of static-concentration time-kill by polymyxin B, meropenem and fosfomycin as monotherapy and in double / triple combination against a polymyxin B-resistant carbapenem-resistant clinical isolate. Two pre-existing bacterial subpopulations, each split over the two-state replication life cycle: intermediate to all three drugs (the majority population) and resistant to all three (seeded at the estimated mutant frequency). Each drug kills through its own Hill function. Adding mechanistic synergy did not improve the fit for this isolate, so no synergy term is present. One of six isolate-specific files from this paper"
  reference <- "Mahadevan R, Garcia E, Sharma R, Qiu H, Elsheikh A, Parambi R, Abboud CS, Pasteran F, Ramirez MS, Kaye KS, Bonomo RA, Rao GG. A mechanism-based pharmacokinetic/pharmacodynamic analysis of polymyxin B-based combination therapy against carbapenem-resistant Klebsiella pneumoniae isolates with diverse phenotypic and genotypic resistance mechanisms. Antimicrob Agents Chemother. 2026 Feb;70(2):e00782-25. doi:10.1128/aac.00782-25. PMCID: PMC12888851. Structural model: Materials and Methods, 'Mechanism based PK/PD model development', equations 1-4. Parameter estimates: Table 2, column BRKP67 (the mechanistic-synergy rows are '-' for this isolate). Isolate resistance genes and MICs: Table 1 and supplemental Table S2. Model-predicted triple-combination time-kill profiles: Figure 3, row BRKP67"
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
    bact_intermediate_intermediate_intermediate1 = list(analyte = "Klebsiella pneumoniae BRKP67, polymyxin B-intermediate / meropenem-intermediate / fosfomycin-intermediate subpopulation, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_intermediate_intermediate_intermediate2 = list(analyte = "Klebsiella pneumoniae BRKP67, polymyxin B-intermediate / meropenem-intermediate / fosfomycin-intermediate subpopulation, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_resistant1 = list(analyte = "Klebsiella pneumoniae BRKP67, subpopulation resistant to all three drugs, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_resistant2 = list(analyte = "Klebsiella pneumoniae BRKP67, subpopulation resistant to all three drugs, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Klebsiella pneumoniae BRKP67, carbapenem-resistant clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    organism = "Klebsiella pneumoniae BRKP67, obtained from an individual patient at Instituto Dante Pazzanese de Cardiologia, Sao Paulo, Brazil",
    resistance_genes = "blaKPC-2; non-functional mgrB (insertion of ISKpn13, a 1148 bp IS5-like element); non-functional ompK35 (terminated at 144 amino acids) and ompK36 (few amino-acid substitutions); ompK37 comparable in size with wild type (Table 1, Table S2)",
    mic_values = c(
      `polymyxin B` = "8 mg/L (resistant, CLSI)",
      meropenem = "64 mg/L (resistant, CLSI)",
      fosfomycin = "32 mg/L (resistant, EUCAST E. coli breakpoint)"
    ),
    system = "24 h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton broth (25.0 mg/L Ca2+, 12.5 mg/L Mg2+) supplemented with 25 mg/L glucose-6-phosphate; viable counts at 0, 1, 2, 4, 6, 8 and 24 h",
    temperature = "37 C",
    duration = "24 h",
    starting_inoculum = "~10^6 CFU/mL; the model estimated 10^5.93 CFU/mL for this isolate",
    limit_of_quantification = "2 log10 CFU/mL; observations below it were handled with the Beal M3 method",
    dose_range = "polymyxin B 0.5-64 mg/L, meropenem 10-120 mg/L, fosfomycin 75-500 mg/L, as monotherapy and in double and triple combinations",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Mechanism-based PK/PD model (MBM) fit in S-ADAPT 1.57 via SADAPT-TRAN with the importance-sampling algorithm (pmethod = 4); observations below the LLOQ were fit with the Beal M3 method.",
      "The correlation coefficient between observed and model-predicted log10 CFU/mL was 0.73 for this isolate (Figure S2).",
      "Mechanistic synergy of polymyxin B on meropenem and fosfomycin was NOT supported for this isolate, so no synergy term is present; it was also unsupported for KP0052-1 and BRKP28, and was retained for BRKP61, BRKP76 and KP0016-1.",
      "Companion files: Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP61, _BRKP76, _KP00161, _KP00521, _BRKP28."
    )
  )

  ini({
    # ---- Bacterial growth and subpopulations (Table 2, column BRKP67) -------
    log10cfu0 <- 5.93; label("Bacterial initial inoculum (log10 CFU/mL)")                                 # Table 2: LogCFU0 = 5.93 (%CV 1.64)
    log10cfumax <- 9.44; label("Maximum population size (log10 CFU/mL)")                                  # Table 2: LogCFUmax = 9.44 (%CV 2.02)
    lk21 <- fixed(log(50)); label("Log replication rate constant, state 2 to state 1 (1/h)")              # Table 2: k21 = 50 (fixed); Methods: "k21 was fixed to 50 h-1"

    mgt_iii <- 32.3; label("Mean generation time, PMB-I / MEM-I / FOF-I subpopulation (min)")             # Table 2: MGT, PMB_I/MEM_I/FOF_I = 32.3 (%CV 12.5)
    mgt_rrr <- 31.9; label("Mean generation time, PMB-R / MEM-R / FOF-R subpopulation (min)")             # Table 2: MGT, PMB_R/MEM_R/FOF_R = 31.9 (%CV 5.08)

    log10mf_rrr <- -4.80; label("Log10 mutant frequency at baseline of the PMB-R / MEM-R / FOF-R subpopulation (unitless fraction)") # Table 2: log10 MF, PMB_R/MEM_R/FOF_R = -4.80 (%CV 1.64)

    # ---- Drug effect of polymyxin B (Table 2) ------------------------------
    kmax_pmb <- 9.08; label("Maximum killing rate constant of polymyxin B (1/h)")                         # Table 2: KmaxPMB = 9.08 (%CV 6.68); Results: "BRKP67 ... showed a slightly higher KmaxPMB of 9.08 h-1"
    kc50_pmb_i <- 2.77; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-intermediate subpopulation (mg/L)")  # Table 2: KC50,PMB,I = 2.77 (%CV 6.59)
    kc50_pmb_r <- 96.78; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-resistant subpopulation (mg/L)")    # Table 2: KC50,PMB,R = 96.78 (%CV 7.07)
    hill_pmb <- 0.82; label("Hill coefficient of polymyxin B killing (unitless)")                         # Table 2: Hill coefficient of polymyxin B = 0.82 (%CV 10.16)

    # ---- Drug effect of meropenem (Table 2) --------------------------------
    kmax_mem <- 4.16; label("Maximum killing rate constant of meropenem (1/h)")                           # Table 2: KmaxMEM = 4.16 (%CV 7.09)
    kc50_mem_i <- 8.91; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-intermediate subpopulation (mg/L)")    # Table 2: KC50,MEM,I = 8.91 (%CV 25.8)
    kc50_mem_r <- 227; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-resistant subpopulation (mg/L)")        # Table 2: KC50,MEM,R = 227 (%CV 12.3)
    hill_mem <- 1.17; label("Hill coefficient of meropenem killing (unitless)")                           # Table 2: Hill coefficient of meropenem = 1.17 (%CV 8.20)

    # ---- Drug effect of fosfomycin (Table 2) -------------------------------
    kmax_fof <- 3.45; label("Maximum killing rate constant of fosfomycin (1/h)")                          # Table 2: KmaxFOF = 3.45 (%CV 11.7)
    kc50_fof_i <- 20.4; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-intermediate subpopulation (mg/L)")   # Table 2: KC50,FOF,I = 20.4 (%CV 22.2)
    kc50_fof_r <- 1000; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-resistant subpopulation (mg/L)")      # Table 2: KC50,FOF,R = 1000 (%CV 7.38)
    hill_fof <- 0.36; label("Hill coefficient of fosfomycin killing (unitless)")                          # Table 2: Hill coefficient of fosfomycin = 0.36 (%CV 9.64)

    # ---- Residual variability (Table 2, Variability model) -----------------
    addSd <- 0.52; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)")              # Table 2: Additive residual variability = 0.52 (%CV 6.85)
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

    # ---- 3. Per-subpopulation killing (equation 4) -------------------------
    # No mechanistic-synergy multiplier: adding it did not improve the fit for
    # this isolate, so the Table 2 Imax / IC50 / Hill synergy rows are "-".
    kill_iii <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_i^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + kc50_mem_i^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + kc50_fof_i^hill_fof)
    kill_rrr <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_r^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + kc50_mem_r^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + kc50_fof_r^hill_fof)

    # ---- 4. Bacterial life-cycle system (equations 2 and 3) ----------------
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

    # ---- 5. Static drug exposures ------------------------------------------
    d/dt(cpmb) <- 0
    d/dt(cmem) <- 0
    d/dt(cfof) <- 0

    # ---- 6. Initial conditions (equations 2 and 3) --------------------------
    bact_intermediate_intermediate_intermediate1(0) <- cfu0 * mf_iii
    bact_resistant_resistant_resistant1(0) <- cfu0 * mf_rrr

    # ---- 7. Observation -----------------------------------------------------
    # The 1e-6 CFU/mL regularisation keeps Cc finite when a regimen drives the
    # population to numerical extinction (the solver can undershoot to a small
    # negative value, and log10 of that is NaN). It is 8 orders of magnitude
    # below the 2 log10 CFU/mL LLOQ, so it never perturbs an observable count.
    # Same device as HernandezLozano_2025_apramycin_invitro.R.
    Cc <- log10(CFUtot + 1e-6)
    Cc ~ add(addSd)
  })
}
