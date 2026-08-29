Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP28 <- function() {
  description <- "In vitro (Klebsiella pneumoniae BRKP28, blaKPC-2 + non-functional MgrB (premature stop codon) + ompK35/ompK36 porin mutations; PMB MIC >128, MEM MIC 256, FOF MIC 128 mg/L). Mechanism-based PK/PD (Bulitta life-cycle growth) model of static-concentration time-kill by polymyxin B, meropenem and fosfomycin as monotherapy and in double / triple combination against an isolate resistant to all three drugs. Two pre-existing bacterial subpopulations, each split over the two-state replication life cycle: polymyxin B-resistant / meropenem-resistant / fosfomycin-intermediate (the majority population) and resistant to all three (seeded at the estimated mutant frequency). Each drug kills through its own Hill function. Adding mechanistic synergy did not improve the fit for this isolate, so no synergy term is present; the estimated maximum polymyxin B killing rate constant (3.61 1/h) is the lowest of the six isolates. One of six isolate-specific files from this paper"
  reference <- "Mahadevan R, Garcia E, Sharma R, Qiu H, Elsheikh A, Parambi R, Abboud CS, Pasteran F, Ramirez MS, Kaye KS, Bonomo RA, Rao GG. A mechanism-based pharmacokinetic/pharmacodynamic analysis of polymyxin B-based combination therapy against carbapenem-resistant Klebsiella pneumoniae isolates with diverse phenotypic and genotypic resistance mechanisms. Antimicrob Agents Chemother. 2026 Feb;70(2):e00782-25. doi:10.1128/aac.00782-25. PMCID: PMC12888851. Structural model: Materials and Methods, 'Mechanism based PK/PD model development', equations 1-4. Parameter estimates: Table 2, column BRKP28 (the KC50,PMB,I, KC50,MEM,I and mechanistic-synergy rows are '-' for this isolate). Isolate resistance genes and MICs: Table 1 and supplemental Table S2. Model-predicted triple-combination time-kill profiles: Figure 3, row BRKP28"
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
    bact_resistant_resistant_intermediate1 = list(analyte = "Klebsiella pneumoniae BRKP28, polymyxin B-resistant / meropenem-resistant / fosfomycin-intermediate subpopulation, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_intermediate2 = list(analyte = "Klebsiella pneumoniae BRKP28, polymyxin B-resistant / meropenem-resistant / fosfomycin-intermediate subpopulation, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_resistant1 = list(analyte = "Klebsiella pneumoniae BRKP28, subpopulation resistant to all three drugs, vegetative state 1", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_resistant_resistant_resistant2 = list(analyte = "Klebsiella pneumoniae BRKP28, subpopulation resistant to all three drugs, replicating state 2", units = "CFU/mL", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "in vitro (Klebsiella pneumoniae BRKP28, carbapenem-resistant clinical isolate)",
    n_subjects = 1L,
    n_studies = 1L,
    organism = "Klebsiella pneumoniae BRKP28, obtained from an individual patient at Instituto Dante Pazzanese de Cardiologia, Sao Paulo, Brazil",
    resistance_genes = "blaKPC-2; non-functional mgrB (premature stop codon); non-functional ompK35 (terminated at 144 amino acids) and ompK36 (not amplified); ompK37 comparable in size with wild type (Table 1, Table S2)",
    mic_values = c(
      `polymyxin B` = ">128 mg/L (resistant, CLSI)",
      meropenem = "256 mg/L (resistant, CLSI)",
      fosfomycin = "128 mg/L (resistant, EUCAST E. coli breakpoint)"
    ),
    system = "24 h static-concentration time-kill (SCTK) assay in cation-adjusted Mueller-Hinton broth (25.0 mg/L Ca2+, 12.5 mg/L Mg2+) supplemented with 25 mg/L glucose-6-phosphate; viable counts at 0, 1, 2, 4, 6, 8 and 24 h",
    temperature = "37 C",
    duration = "24 h",
    starting_inoculum = "~10^6 CFU/mL; the model estimated 10^6.32 CFU/mL for this isolate",
    limit_of_quantification = "2 log10 CFU/mL; observations below it were handled with the Beal M3 method",
    dose_range = "polymyxin B 0.5-64 mg/L, meropenem 10-120 mg/L, fosfomycin 75-500 mg/L, as monotherapy and in double and triple combinations",
    disease_state = "not applicable (in vitro)",
    notes = paste(
      "Mechanism-based PK/PD model (MBM) fit in S-ADAPT 1.57 via SADAPT-TRAN with the importance-sampling algorithm (pmethod = 4); observations below the LLOQ were fit with the Beal M3 method.",
      "The correlation coefficient between observed and model-predicted log10 CFU/mL was 0.91 for this isolate (Figure S2).",
      "This isolate is resistant to all three drugs and was the least responsive of the six: it has the lowest estimated KmaxPMB (3.61 1/h) and the largest 0-24 h AUC_CFU (50.3 vs 40.1 log10 CFU*h/mL for BRKP67). The Discussion attributes the residual non-response to additional uncharacterized resistance mechanisms.",
      "Mechanistic synergy of polymyxin B on meropenem and fosfomycin was NOT supported for this isolate, so no synergy term is present; it was also unsupported for KP0052-1 and BRKP67, and was retained for BRKP61, BRKP76 and KP0016-1.",
      "Companion files: Mahadevan_2026_polymyxinB_meropenem_fosfomycin_BRKP61, _BRKP76, _KP00161, _KP00521, _BRKP67."
    )
  )

  ini({
    # ---- Bacterial growth and subpopulations (Table 2, column BRKP28) -------
    log10cfu0 <- 6.32; label("Bacterial initial inoculum (log10 CFU/mL)")                                 # Table 2: LogCFU0 = 6.32 (%CV 0.45)
    log10cfumax <- 9.35; label("Maximum population size (log10 CFU/mL)")                                  # Table 2: LogCFUmax = 9.35 (%CV 0.39)
    lk21 <- fixed(log(50)); label("Log replication rate constant, state 2 to state 1 (1/h)")              # Table 2: k21 = 50 (fixed); Methods: "k21 was fixed to 50 h-1"

    mgt_rri <- 44.25; label("Mean generation time, PMB-R / MEM-R / FOF-I subpopulation (min)")            # Table 2: MGT, PMB_R/MEM_R/FOF_I = 44.25 (%CV 6.50)
    mgt_rrr <- 19.1; label("Mean generation time, PMB-R / MEM-R / FOF-R subpopulation (min)")             # Table 2: MGT, PMB_R/MEM_R/FOF_R = 19.1 (%CV 7.39)

    log10mf_rrr <- -6.55; label("Log10 mutant frequency at baseline of the PMB-R / MEM-R / FOF-R subpopulation (unitless fraction)") # Table 2: log10 MF, PMB_R/MEM_R/FOF_R = -6.55 (%CV 15.4)

    # ---- Drug effect of polymyxin B (Table 2) ------------------------------
    # Both subpopulations are polymyxin B-resistant, so a single KC50,PMB,R
    # applies to both; Table 2 reports "-" for KC50,PMB,I.
    kmax_pmb <- 3.61; label("Maximum killing rate constant of polymyxin B (1/h)")                         # Table 2: KmaxPMB = 3.61 (%CV 13.4); Abstract: "BRKP28 had a reduced maximum killing rate constant for PMB (3.61 h-1)"
    kc50_pmb_r <- 38.82; label("Polymyxin B concentration causing 50% of KmaxPMB in the PMB-resistant subpopulations (mg/L)")   # Table 2: KC50,PMB,R = 38.82 (%CV 13.35)
    hill_pmb <- 1.44; label("Hill coefficient of polymyxin B killing (unitless)")                         # Table 2: Hill coefficient of polymyxin B = 1.44 (%CV 20.82)

    # ---- Drug effect of meropenem (Table 2) --------------------------------
    # Both subpopulations are meropenem-resistant; Table 2 reports "-" for
    # KC50,MEM,I.
    kmax_mem <- 3.70; label("Maximum killing rate constant of meropenem (1/h)")                           # Table 2: KmaxMEM = 3.70 (%CV 16.5); Results: "BRKP28 (MIC: 256 mg/L) ... demonstrated correspondingly lower Kmax values for meropenem (... 3.70 h-1)"
    kc50_mem_r <- 228; label("Meropenem concentration causing 50% of KmaxMEM in the MEM-resistant subpopulations (mg/L)")       # Table 2: KC50,MEM,R = 228 (%CV 11.1)
    hill_mem <- 1.43; label("Hill coefficient of meropenem killing (unitless)")                           # Table 2: Hill coefficient of meropenem = 1.43 (%CV 12.2)

    # ---- Drug effect of fosfomycin (Table 2) -------------------------------
    kmax_fof <- 3.42; label("Maximum killing rate constant of fosfomycin (1/h)")                          # Table 2: KmaxFOF = 3.42 (%CV 7.37)
    kc50_fof_i <- 22.0; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-intermediate subpopulation (mg/L)")   # Table 2: KC50,FOF,I = 22.0 (%CV 20.6)
    kc50_fof_r <- 1342; label("Fosfomycin concentration causing 50% of KmaxFOF in the FOF-resistant subpopulation (mg/L)")      # Table 2: KC50,FOF,R = 1342 (%CV 15.9)
    hill_fof <- 0.35; label("Hill coefficient of fosfomycin killing (unitless)")                          # Table 2: Hill coefficient of fosfomycin = 0.35 (%CV 23.73)

    # ---- Residual variability (Table 2, Variability model) -----------------
    addSd <- 0.25; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)")              # Table 2: Additive residual variability = 0.25 (%CV 5.02)
  })

  model({
    # ---- 1. Back-transforms and derived rate constants ---------------------
    k21 <- exp(lk21)
    cfu0 <- 10^log10cfu0
    cfumax <- 10^log10cfumax
    k12_rri <- 60 / mgt_rri
    k12_rrr <- 60 / mgt_rrr
    mf_rrr <- 10^log10mf_rrr
    mf_rri <- 1 - mf_rrr

    # ---- 2. Total viable count and the replication plateau (equation 1) ----
    CFUtot <- bact_resistant_resistant_intermediate1 +
      bact_resistant_resistant_intermediate2 +
      bact_resistant_resistant_resistant1 +
      bact_resistant_resistant_resistant2
    plat <- 1 - CFUtot / (cfumax + CFUtot)

    # ---- 3. Per-subpopulation killing (equation 4) -------------------------
    # No mechanistic-synergy multiplier (Results). The two subpopulations share
    # the polymyxin B and meropenem KC50 because both are resistant to those
    # drugs; they differ only in the fosfomycin KC50.
    kill_rri <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_r^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + kc50_mem_r^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + kc50_fof_i^hill_fof)
    kill_rrr <- kmax_pmb * cpmb^hill_pmb / (cpmb^hill_pmb + kc50_pmb_r^hill_pmb) +
      kmax_mem * cmem^hill_mem / (cmem^hill_mem + kc50_mem_r^hill_mem) +
      kmax_fof * cfof^hill_fof / (cfof^hill_fof + kc50_fof_r^hill_fof)

    # ---- 4. Bacterial life-cycle system (equations 2 and 3) ----------------
    d/dt(bact_resistant_resistant_intermediate1) <-
      2 * plat * k21 * bact_resistant_resistant_intermediate2 -
      k12_rri * bact_resistant_resistant_intermediate1 -
      kill_rri * bact_resistant_resistant_intermediate1
    d/dt(bact_resistant_resistant_intermediate2) <-
      -k21 * bact_resistant_resistant_intermediate2 +
      k12_rri * bact_resistant_resistant_intermediate1 -
      kill_rri * bact_resistant_resistant_intermediate2

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
    bact_resistant_resistant_intermediate1(0) <- cfu0 * mf_rri
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
