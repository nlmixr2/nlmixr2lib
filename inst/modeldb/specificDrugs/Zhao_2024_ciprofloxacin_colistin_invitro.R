Zhao_2024_ciprofloxacin_colistin_invitro <- function() {
  description <- "In vitro (Escherichia coli; one clinical urinary isolate C47 plus MG1655 wild type LM347 and its isogenic gyrA1/marR mutants LM378 and LM421). Semi-mechanistic PK/PD model of static-concentration time-kill for ciprofloxacin and colistin alone and in combination. Five pre-existing subpopulations (the 2x2 ciprofloxacin/colistin susceptible-resistant grid plus a pre-existing resting subpopulation present only in C47), each cycling between a growing S state, a density-driven resting R state and a transient ciprofloxacin-induced non-colony-forming Nc state. Killing is additive across drugs (subpopulation synergy) with concentration-independent interaction factors shifting the ciprofloxacin and colistin EC50 in combination. Colistin loss to labware is described by an Emax free-fraction model at time zero plus saturable binding kinetics during the experiment; ciprofloxacin is held at its nominal concentration. Sibling models: Zhao_2024_ciprofloxacin_colistin_plasma, Zhao_2024_ciprofloxacin_colistin_kidney."
  reference <- "Zhao C, Kristoffersson AN, Khan DD, Lagerback P, Lustig U, Cao S, Annerstedt C, Cars O, Andersson DI, Hughes D, Nielsen EI, Friberg LE. Quantifying combined effects of colistin and ciprofloxacin against Escherichia coli in an in silico pharmacokinetic-pharmacodynamic model. Sci Rep. 2024 May 22;14(1):11603. doi:10.1038/s41598-024-61518-0. Structural equations: main text Eqs 1-17. Parameter estimates: main text Table 1 (PK/PD model) and Supplementary Table S1 (colistin binding model). Every value here was cross-checked against the authors' deposited final NONMEM control stream (Supplementary zip, Supplementary/run422b_clean.mod), which carries the final estimates in its $THETA/$SIGMA records."
  vignette <- "Zhao_2024_ciprofloxacin_colistin"
  units <- list(time = "h", dosing = "CFU/mL", concentration = "mg/L")

  # Fifteen states of the deposited control stream ($MODEL of run422b_clean.mod).
  # Subpopulation indices follow Fig. 1 of the paper exactly:
  #   1 = ciprofloxacin susceptible, colistin susceptible
  #   2 = ciprofloxacin susceptible, colistin resistant
  #   3 = ciprofloxacin resistant,   colistin susceptible
  #   4 = ciprofloxacin resistant,   colistin resistant
  # Within each, S = growing, R = resting/non-growing, Nc = transiently
  # non-colony-forming (filamented under ciprofloxacin). bact_pr is the
  # pre-existing resting ("persister") subpopulation, present only in C47.
  # cst_unbound / cst_bound are the free and labware-bound colistin in the tube.
  paper_specific_compartments <- c(
    "bact_s1", "bact_r1", "bact_nc1",
    "bact_s2", "bact_r2", "bact_nc2",
    "bact_s3", "bact_r3", "bact_nc3",
    "bact_s4", "bact_r4", "bact_nc4",
    "bact_pr", "cst_unbound", "cst_bound"
  )

  # Every state is a concentration (CFU/mL or mg/L) in the 2 mL time-kill
  # tube, so no state carries a separate distribution volume. `specimen` is
  # "not applicable" throughout because the controlled vocabulary covers
  # clinical matrices only; the actual matrix is cation-adjusted Mueller
  # Hinton II broth (and, for cst_bound, the polypropylene tube wall).
  compartmentData <- list(
    bact_s1 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-susceptible subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r1 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-susceptible subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc1 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-susceptible subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s2 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-resistant subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r2 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-resistant subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc2 = list(analyte = "Escherichia coli, ciprofloxacin-susceptible and colistin-resistant subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s3 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-susceptible subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r3 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-susceptible subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc3 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-susceptible subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_s4 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-resistant subpopulation, growing (S) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_r4 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-resistant subpopulation, resting (R) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_nc4 = list(analyte = "Escherichia coli, ciprofloxacin-resistant and colistin-resistant subpopulation, non-colony-forming (Nc) state", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    bact_pr = list(analyte = "Escherichia coli, pre-existing resting (persister) subpopulation, present only in strain C47", units = "CFU/mL", specimen = "not applicable", verified = TRUE),
    cst_unbound = list(analyte = "colistin, unbound in the broth", units = "mg/L", specimen = "not applicable", verified = TRUE),
    cst_bound = list(analyte = "colistin, adsorbed to the polypropylene tube wall", units = "mg/L", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    BACT = list(
      description = "Bacterial strain indicator. The strain selects the colistin PD parameters that Table 1 reports per strain, and the fraction of pre-existing resting bacteria in the inoculum.",
      units = "(categorical)",
      type = "categorical",
      source_name = "FSTRAIN",
      notes = paste(
        "Level coding is taken verbatim from the FSTRAIN column of the authors' NONMEM dataset (CIPCST_clean.csv) and the IF(FSTRAIN.EQ.nn) blocks of run422b_clean.mod:",
        "47 = C47, a clinical isolate from a patient with uncomplicated urinary tract infection (MIC_CIP 0.047 mg/L, MIC_CST 0.75 mg/L); the only strain with a pre-existing resting subpopulation.",
        "347 = LM347, the E. coli MG1655 wild-type laboratory strain (MIC_CIP 0.023 mg/L, MIC_CST 0.5 mg/L).",
        "378 = LM378, isogenic mutant carrying gyrA1 S83L (MIC_CIP 0.38 mg/L, MIC_CST 0.5 mg/L).",
        "421 = LM421, isogenic mutant carrying gyrA1 S83L plus a marR knockout (MIC_CIP 1.0 mg/L, MIC_CST 0.5 mg/L).",
        "LM378 and LM421 share one set of colistin PD parameters (Table 1 column 'LM378, LM421'). There is no reference level: every strain selects its own parameter set, and the model has no reference-versus-test contrast."
      )
    )
  )

  population <- list(
    species = "in vitro (Escherichia coli; four strains)",
    n_subjects = 303L,
    n_studies = 1L,
    disease_state = "Four E. coli strains of differing ciprofloxacin susceptibility. C47 is a clinical isolate from a patient with an uncomplicated urinary tract infection (MIC_CIP 0.047 mg/L, MIC_CST 0.75 mg/L). LM347 is the MG1655 wild type (MIC_CIP 0.023 mg/L, MIC_CST 0.5 mg/L). LM378 (gyrA1 S83L) and LM421 (gyrA1 S83L plus marR knockout) are isogenic laboratory mutants with MIC_CIP 0.38 and 1.0 mg/L respectively and MIC_CST 0.5 mg/L for both. MICs were determined by macrodilution in cation-adjusted Mueller Hinton II broth.",
    model_system = "Static-concentration time-kill experiments in cation-adjusted Mueller Hinton II (CAMHII) broth at 37 degrees C in 2 mL Falcon polypropylene tubes, sampled for viable count on CAMHII agar. Colistin monodrug: 0.125-16 x MIC plus growth control, sampled at 0, 1, 2, 4, 6, 9, 12, 24 and 32 h. Combination: ciprofloxacin 0.125-0.75 x MIC with colistin 0.125-0.375 x MIC (concentrations and sampling times chosen by optimal design in PopED), plus 2 x MIC monodrug controls and a growth control, sampled at 0, 0.75, 2, 4, 6, 8, 24 and 32 h. Ciprofloxacin monodrug curves were re-used from the earlier dataset of Khan 2015 (ref 19 of the paper). Experiments were run in at least duplicate, in different laboratories, by several technicians, over several years.",
    initial_inoculum = "~10^6 CFU/mL (log-phase culture diluted 1:100 twice from an overnight culture)",
    dose_range = "Ciprofloxacin 0.125-2 x MIC; colistin 0.125-16 x MIC (nominal 0.0625-8 mg/L)",
    notes = paste(
      "303 time-kill curves in total: 129 ciprofloxacin (1976 CFU counts, 16 experimental days), 102 colistin (1701 counts, 15 days) and 72 combination (1028 counts, 8 days) - Supplementary Table S2.",
      "Estimated in NONMEM 7.5.0 with the LAPLACIAN method and ADVAN13/TOL9; parameter uncertainty by sampling importance resampling in PsN 5.3.0. Data were transformed to log10 CFU/mL (transform-both-sides) and the M3 method handled the 1 log10 CFU/mL limit of detection.",
      "All ciprofloxacin monodrug parameters were FIXED to the values reported by Khan 2015 (ref 19), which were estimated on a larger E. coli dataset that included these four strains; only the colistin, combination, Bmax and residual parameters were estimated here.",
      "Emax of colistin on the colistin-susceptible subpopulations was fixed to 50 /h (a ~1 min minimum half-life). A sensitivity analysis over 30, 50, 100 and 150 found 50 to be the lowest value that did not compromise the fit while stabilising the model.",
      "Subpopulation 4 (resistant to both drugs) never exceeded 1 CFU/mL in any investigated scenario but is retained for completeness, as the authors did.",
      "Bacteria in the Nc state do not form colonies and are therefore excluded from the observed count, but they DO count towards the total density that drives the S-to-R transfer (Eq. 10). The model keeps the two totals separate.",
      "Sibling models applying this PD model to simulated human exposures: Zhao_2024_ciprofloxacin_colistin_plasma (bloodstream infection) and Zhao_2024_ciprofloxacin_colistin_kidney (pyelonephritis)."
    )
  )

  ini({
    # ---- Experimental design inputs -------------------------------------
    # These are data columns in the authors' NONMEM dataset (AMT, CIPT, CSTT,
    # CIPMIC), not estimated parameters. They are exposed here as fixed ini()
    # entries so a user can set them per simulated time-kill curve with
    # rxSolve(params = c(...)). Defaults reproduce a C47 combination curve.
    inoc        <- fixed(1e6);   label("Start inoculum (CFU/mL); experimental design input, AMT column") # Methods, Time-kill experiments: start inoculum around 10^6 CFU/mL
    cipNominal  <- fixed(0.047); label("Nominal (constant) ciprofloxacin concentration (mg/L); experimental design input, CIPT column") # run422b_clean.mod $DES: CIPC = CIPT, i.e. ciprofloxacin is held at its nominal concentration
    cstNominal  <- fixed(0.75);  label("Nominal colistin concentration added at time zero (mg/L); experimental design input, CSTT column") # run422b_clean.mod $PK: A_0(13) = CSTT*FF
    micCip      <- fixed(0.047); label("Ciprofloxacin MIC of the strain (mg/L); experimental design input, CIPMIC column") # Methods: MIC_CIP 0.047 (C47), 0.023 (LM347), 0.38 (LM378), 1.0 (LM421) mg/L

    # ---- Bacterial growth model (Table 1, "Bacterial growth model") -----
    kgs  <- fixed(1.69);  label("Growth rate constant of the ciprofloxacin-susceptible subpopulations 1 and 2 (1/h)") # Table 1: k_g1,2 = 1.69 FIX (fixed to Khan 2015); run422b_clean.mod THETA(1)
    kgr  <- fixed(0.316); label("Growth rate constant of the ciprofloxacin-resistant subpopulations 3 and 4 (1/h)")  # Table 1: k_g3,4 = 0.316 FIX (fitness cost); run422b_clean.mod THETA(16)
    kd   <- fixed(0.179); label("Natural death rate constant (1/h)")                                                  # Table 1: k_d = 0.179 FIX; run422b_clean.mod $PK KD = 0.179
    lbmax <- log(1.56e9); label("Log of the maximum bacterial density in the system (log CFU/mL)")                     # Table 1: Bmax = 1.56e9 CFU/mL (11% RSE); run422b_clean.mod THETA(2) = 1.55754 x 10^9. Carried on the log scale per the library naming convention and back-transformed in model()

    # ---- Ciprofloxacin monodrug model (Table 1; all FIXED to Khan 2015) --
    emaxCip     <- fixed(6.75);  label("Emax of ciprofloxacin on k_CIP,1,2 and k_CIP,3,4 (1/h)")                      # Table 1: E_maxCIP = 6.75 FIX; run422b_clean.mod THETA(18)
    ec50CipSlp  <- fixed(1.38);  label("Multiplier of the ciprofloxacin EC50-MIC power relation for subpopulations 1 and 2 (mg/L)") # Table 1: EC50_CIP,1,2 = 1.38 x MIC^0.996 FIX; run422b_clean.mod THETA(19)
    ec50CipPwr  <- fixed(0.996); label("Exponent of the ciprofloxacin EC50-MIC power relation for subpopulations 1 and 2")          # Table 1: EC50_CIP,1,2 = 1.38 x MIC^0.996 FIX; run422b_clean.mod THETA(20)
    ec50CipR    <- fixed(1.53);  label("EC50 of ciprofloxacin on the ciprofloxacin-resistant subpopulations 3 and 4 (mg/L)")        # Table 1: EC50_CIP,3,4 = 1.53 FIX; run422b_clean.mod THETA(21)
    hillCip     <- fixed(1.60);  label("Hill factor of ciprofloxacin on k_CIP,1,2 and k_CIP,3,4")                     # Table 1: gamma_CIP = 1.60 FIX; run422b_clean.mod THETA(22)
    fmutCip     <- fixed(1.23e-6); label("Fraction of ciprofloxacin-resistant bacteria in the inoculum")               # Table 1: f_MUTCIP = 1.23 per 10^6 FIX; run422b_clean.mod THETA(17) x 10^-6
    fmutPrC47   <- fixed(3.91e-3); label("Fraction of pre-existing resting (PR) bacteria in the inoculum, strain C47") # Table 1: f_MUTP = 3910 per 10^6 FIX for C47, 0 for the lab strains; run422b_clean.mod THETA(28) x 10^-6
    ksncMax     <- fixed(4.69);  label("Emax of the S-to-Nc transfer rate constant k_SNc (1/h)")                      # Table 1: k_SNc,max = 4.69 FIX; run422b_clean.mod THETA(23)
    ksnc50      <- fixed(0.183); label("EC50 of the S-to-Nc transfer rate constant, on the C_CIP/EC50_CIP ratio scale (tr50)") # Table 1: t_r50 = 0.183 FIX; run422b_clean.mod THETA(24)
    hillSnc     <- fixed(20);    label("Hill factor of the S-to-Nc transfer rate constant k_SNc")                     # Table 1: gamma_SNc = 20 FIX; run422b_clean.mod THETA(25)
    kncs        <- fixed(0.449); label("Scaling factor of the Nc-to-S transfer rate constant k_NcS (1/h)")            # Table 1: k_NcS = 0.449 FIX; run422b_clean.mod THETA(26), sfNpS
    tnc         <- fixed(3.2);   label("Time after which transfer from S to Nc is switched off (h)")                  # Table 1: t_Nc = 3.2 FIX; run422b_clean.mod THETA(27) via MTIME(2)

    # ---- Colistin monodrug model (Table 1, "CST monodrug model") --------
    emaxCstS      <- fixed(50);    label("Emax of colistin on the colistin-susceptible subpopulations 1 and 3 (1/h)")   # Table 1: E_maxCST,1,3 = 50 FIX (fixed to stabilise the model; sensitivity analysis over 30/50/100/150); run422b_clean.mod THETA(6)
    emaxCstRC47   <- 1.44;         label("Emax of colistin on the colistin-resistant subpopulations 2 and 4, strain C47 (1/h)")        # Table 1: E_maxCST,2,4 = 1.44 (1.1% RSE) for C47; run422b_clean.mod THETA(7) = 1.43838
    emaxCstRLab   <- 2.44;         label("Emax of colistin on the colistin-resistant subpopulations 2 and 4, strains LM347/LM378/LM421 (1/h)") # Table 1: E_maxCST,2,4 = 2.44 (1.1% RSE); run422b_clean.mod THETA(8) = 2.43987, shared with E_maxCST,PR (Table 1 footnote 4)
    ec50Cst       <- 0.110;        label("EC50 of colistin for k_CST,1,3 and k_CST,2,4, shared across strains (mg/L)")  # Table 1: EC50_CST = 0.110 (1.4% RSE); run422b_clean.mod THETA(9) = 109.81 x 0.001
    emaxCstPr     <- 2.44;         label("Emax of colistin on the pre-existing resting subpopulation, strain C47 (1/h)") # Table 1: E_maxCST,PR = 2.44 (1.1% RSE); run422b_clean.mod EMAXCST3 = THETA(8), the same estimate as emaxCstRLab (Table 1 footnote 4)
    ec50CstPrFrac <- 1.2807;       label("Fractional increase of the colistin EC50 for the pre-existing resting subpopulation; EC50_CST,PR = EC50_CST x (1 + this)") # run422b_clean.mod THETA(15): EC50CST3 = EC50CST*(1+THETA(15)); 0.110 x (1 + 1.2807) = 0.251 mg/L = Table 1 EC50_CST,PR (16% RSE)
    hillCstSC47   <- 5.26;         label("Hill factor of colistin on subpopulations 1 and 3, strain C47")               # Table 1: gamma_CST,1,3 = 5.26 (1.6% RSE) for C47; run422b_clean.mod THETA(10) = 5.25557
    hillCstS347   <- 3.55;         label("Hill factor of colistin on subpopulations 1 and 3, strain LM347")             # Table 1: gamma_CST,1,3 = 3.55 (2.0% RSE) for LM347; run422b_clean.mod THETA(11) = 3.55309
    hillCstSLm    <- fixed(20);    label("Hill factor of colistin on subpopulations 1 and 3, strains LM378 and LM421")  # Table 1: gamma_CST,1,3 = 20 FIX for LM378/LM421; run422b_clean.mod THETA(12)
    hillCstRLab   <- 0.270;        label("Hill factor of colistin on subpopulations 2 and 4, strains LM347/LM378/LM421") # Table 1: gamma_CST,2,4 = 0.270 (5.7% RSE); run422b_clean.mod THETA(13) = 0.269542. For C47, gamma_CST,2,4 equals gamma_CST,1,3 (Table 1 footnote 4)
    fmutCst       <- 3.80e-6;      label("Fraction of colistin-resistant bacteria in the inoculum")                     # Table 1: f_MUTCST = 3.80 per 10^6 (1.7% RSE); run422b_clean.mod THETA(14) = 3.80213 x 10^-6

    # ---- Combination model (Table 1, "CIP and CST combination model") ---
    intCstC47 <- 0.541; label("Factor on the colistin EC50 for k_CST,1,3 and k_CST,PR in combination, strains C47/LM378/LM421") # Table 1: I_CST = 0.541 (1.3% RSE) for C47 and, per footnote 4, for LM378/LM421; run422b_clean.mod THETA(29) = 0.541187
    intCst347 <- 1.19;  label("Factor on the colistin EC50 for k_CST,1,3 in combination, strain LM347")                         # Table 1: I_CST = 1.19 (3.2% RSE) for LM347; run422b_clean.mod THETA(31) = 1.18776
    intCip    <- 1.60;  label("Factor on the ciprofloxacin EC50 for k_CIP,1,2 in combination, shared across strains")           # Table 1: I_CIP = 1.60 (0.75% RSE); run422b_clean.mod THETA(30) = 1.59999

    # ---- Colistin labware binding model (Supplementary Table S1) --------
    fuMin   <- fixed(0.328);  label("Minimum free fraction of colistin at the start of the experiment")             # Table S1: fu_min = 0.328 (10.7% RSE); run422b_clean.mod LOGBASE = -0.719, BASE = exp(-0.719)/(1+exp(-0.719)) = 0.3276
    fuMax   <- fixed(1);      label("Maximum free fraction of colistin at the start of the experiment")             # Table S1: fu_max = 1 FIX (footnote 2: estimating it did not improve the fit)
    fuC50   <- fixed(0.522);  label("Nominal colistin concentration giving half the maximum free fraction (mg/L)")  # Table S1: fu_c50 = 0.522 (18.4% RSE); run422b_clean.mod F50 = 0.522
    hillFu  <- fixed(1);      label("Sigmoid factor of the free-fraction model")                                    # run422b_clean.mod LDIL = 1; Table S1 does not list gamma_fu, consistent with it being fixed to 1 (a plain Emax function)
    kbinMax <- fixed(0.189);  label("Maximum labware binding rate constant of colistin (1/h)")                      # Table S1: k_binMAX = 0.189 (5.4% RSE); run422b_clean.mod THETA(4)
    b50     <- fixed(0.0347); label("Bound colistin concentration giving half the maximum binding rate (mg/L)")     # Table S1: B50 = 0.0347 (13.5% RSE); run422b_clean.mod THETA(5)
    kunb    <- fixed(0.122);  label("Rate constant of colistin dissociation from labware (1/h)")                    # Table S1: k_unb = 0.122 (9.5% RSE); run422b_clean.mod THETA(3)
    hillBin <- fixed(1);      label("Sigmoid factor of the labware binding rate model")                             # run422b_clean.mod $PK HILL = 1

    # ---- Residual error --------------------------------------------------
    # run422b_clean.mod uses a two-level additive residual on the log10 scale:
    # EPS(1) (variance 1.66489, shared within an L2 group = one tube at one
    # sampling time) plus one of EPS(2)-EPS(5) (variance 0.0205917, BLOCK(1)
    # SAME, selected by the FL2 replicate-plate index). rxode2 has no analogue
    # of NONMEM's L2 data item, so the two components are combined into the
    # single total SD that the control stream itself forms in $ERROR as
    # W = SQRT(SIGMA(1) + SIGMA(2)) and uses for the M3 likelihood.
    # See the vignette Assumptions section.
    addSd <- 1.2983; label("Additive residual SD on the log10 bacterial count (log10 CFU/mL)") # run422b_clean.mod $ERROR: W = SQRT(SIGMA(1) + SIGMA(2)) = sqrt(1.66489 + 0.0205917) = 1.2983. Table 1 reports the two components separately as 1.29 (= sqrt(1.66489)) and 0.0144; see the vignette Errata for the replicate-component discrepancy
  })

  model({
    # Guard against log10() of a non-positive solver excursion.
    eps <- 1e-30
    bmax <- exp(lbmax)

    # ---- 1. Strain-specific parameter selection (BACT) ------------------
    # run422b_clean.mod encodes these as IF(FSTRAIN.EQ.nn) blocks in $PK.
    # C47 shares one Hill factor across the colistin-susceptible, colistin-
    # resistant and pre-existing-resting effects (Table 1 footnote 4).
    if (BACT == 47) {
      emaxCstR <- emaxCstRC47
      hillCstS <- hillCstSC47
      hillCstR <- hillCstSC47
      hillCstPr <- hillCstSC47
      emaxCstPrStrain <- emaxCstPr
      ec50CstPr <- ec50Cst * (1 + ec50CstPrFrac)
      fmutPr <- fmutPrC47
      intCst <- intCstC47
    } else if (BACT == 347) {
      emaxCstR <- emaxCstRLab
      hillCstS <- hillCstS347
      hillCstR <- hillCstRLab
      hillCstPr <- hillCstS347
      emaxCstPrStrain <- 0
      # Inert: the laboratory strains have no pre-existing resting
      # subpopulation (E_maxCST,PR = 0 and f_MUTP = 0), so k_CST,PR is
      # identically zero. EC50_CST is used here rather than 0 only to keep the
      # denominator of the k_CST,PR Emax expression strictly positive.
      ec50CstPr <- ec50Cst
      fmutPr <- 0
      intCst <- intCst347
    } else {
      # LM378 and LM421 share one parameter set (Table 1 column "LM378, LM421").
      emaxCstR <- emaxCstRLab
      hillCstS <- hillCstSLm
      hillCstR <- hillCstRLab
      hillCstPr <- hillCstSLm
      emaxCstPrStrain <- 0
      # Inert, as for LM347 above.
      ec50CstPr <- ec50Cst
      fmutPr <- 0
      intCst <- intCstC47
    }

    # ---- 2. Ciprofloxacin EC50 from the strain MIC (Table 1) -------------
    ec50Cip <- ec50CipSlp * micCip^ec50CipPwr

    # ---- 3. Interaction factors (Eq. 17) --------------------------------
    # Both factors are concentration-independent over the studied range and
    # apply only when BOTH drugs are present, and only to the susceptible
    # subpopulations k_CIP,1,2 and k_CST,1,3 (and k_CST,PR).
    # run422b_clean.mod initialises ICST/ICIP to 1 and overrides them inside
    # IF(FSTRAIN.EQ.nn .AND. CSTT.GT.0 .AND. CIPT.GT.0).
    if (cstNominal > 0 && cipNominal > 0) {
      iCst <- intCst
      iCip <- intCip
    } else {
      iCst <- 1
      iCip <- 1
    }

    # ---- 4. Colistin free fraction at time zero (Eq. 13) ----------------
    fu0 <- fuMin + (fuMax - fuMin) * cstNominal^hillFu / (fuC50^hillFu + cstNominal^hillFu)

    # ---- 5. Drug concentrations driving the effect ----------------------
    # Ciprofloxacin is held at its nominal concentration throughout the
    # experiment ($DES: CIPC = CIPT); colistin acts through its unbound
    # concentration in the tube, which declines by labware binding.
    cipConc <- cipNominal
    cstConc <- cst_unbound

    # ---- 6. Ciprofloxacin killing (Eq. 12) ------------------------------
    # I_CIP shifts the EC50 of the ciprofloxacin-susceptible subpopulations
    # only; the resistant subpopulations keep EC50_CIP,3,4 unchanged.
    # Guarded exactly as $DES guards KKCIP1/KKCIP2 with IF(CIPT>0).
    if (cipNominal > 0) {
      kCipS <- emaxCip * cipConc^hillCip / ((ec50Cip * iCip)^hillCip + cipConc^hillCip)
      kCipR <- emaxCip * cipConc^hillCip / (ec50CipR^hillCip + cipConc^hillCip)
    } else {
      kCipS <- 0
      kCipR <- 0
    }

    # ---- 7. Colistin killing (Eq. 16) -----------------------------------
    # I_CST shifts the EC50 of the colistin-susceptible subpopulations and of
    # the pre-existing resting subpopulation only.
    # Guarded exactly as $DES guards KKCST1/KKCST2/KKCST3 with IF(CSTT>0).
    # The guard is load-bearing, not cosmetic: for the three laboratory strains
    # both E_maxCST,PR and EC50_CST,PR are zero, so without it the k_CST,PR term
    # would evaluate 0 * 0^gamma / (0^gamma + 0^gamma) = 0/0 = NaN in a colistin-
    # free (ciprofloxacin monodrug or growth control) experiment.
    if (cstNominal > 0) {
      kCstS <- emaxCstS * cstConc^hillCstS / ((ec50Cst * iCst)^hillCstS + cstConc^hillCstS)
      kCstR <- emaxCstR * cstConc^hillCstR / (ec50Cst^hillCstR + cstConc^hillCstR)
      kCstPr <- emaxCstPrStrain * cstConc^hillCstPr /
        ((ec50CstPr * iCst)^hillCstPr + cstConc^hillCstPr)
    } else {
      kCstS <- 0
      kCstR <- 0
      kCstPr <- 0
    }

    # ---- 8. Additive combined killing per subpopulation -----------------
    # A simple sum of the single-drug rate constants (subpopulation synergy);
    # the interaction enters only through the shifted EC50 values above.
    kTot1 <- kCipS + kCstS
    kTot2 <- kCipS + kCstR
    kTot3 <- kCipR + kCstS
    kTot4 <- kCipR + kCstR
    # The pre-existing resting subpopulation is resistant to ciprofloxacin
    # alone but regains ciprofloxacin susceptibility when colistin is
    # co-administered (Drug combination effect model, and $DES KKTOT5).
    if (cstNominal > 0 && cipNominal > 0) {
      kTotPr <- kCipS + kCstPr
    } else {
      kTotPr <- kCstPr
    }

    # ---- 9. S-to-Nc and Nc-to-S transfer (CIP-driven) -------------------
    # Guarded exactly as in $DES: both rate constants are identically zero
    # without ciprofloxacin (k_NcS would otherwise divide by zero).
    if (cipConc > 0) {
      ksncS <- ksncMax * (cipConc / ec50Cip)^hillSnc /
        ((cipConc / ec50Cip)^hillSnc + ksnc50^hillSnc)
      ksncR <- ksncMax * (cipConc / ec50CipR)^hillSnc /
        ((cipConc / ec50CipR)^hillSnc + ksnc50^hillSnc)
      kncsS <- kncs * ec50Cip / cipConc
      kncsR <- kncs * ec50CipR / cipConc
    } else {
      ksncS <- 0
      ksncR <- 0
      kncsS <- 0
      kncsR <- 0
    }
    # MTIME(2) shuts the S-to-Nc transfer off after t_Nc hours.
    if (t <= tnc) {
      flagNc <- 1
    } else {
      flagNc <- 0
    }

    # ---- 10. Density-dependent S-to-R transfer (Eqs. 9 and 10) ----------
    # B_TOT includes the non-colony-forming states; the observed count below
    # does not.
    btot <- bact_s1 + bact_r1 + bact_nc1 +
      bact_s2 + bact_r2 + bact_nc2 +
      bact_s3 + bact_r3 + bact_nc3 +
      bact_s4 + bact_r4 + bact_nc4 + bact_pr
    ksrS <- btot / bmax * (kgs - kd)
    ksrR <- btot / bmax * (kgr - kd)

    # ---- 11. Growth cut-off below one bacterium in the tube -------------
    # $DES: growth of a subpopulation stops once fewer than 0.5 bacteria
    # remain in the 2 mL culture, preventing regrowth from a fractional cell.
    vol <- 2
    if (bact_s1 * vol < 0.5) {
      grow1 <- 0
    } else {
      grow1 <- 1
    }
    if (bact_s2 * vol < 0.5) {
      grow2 <- 0
    } else {
      grow2 <- 1
    }
    if (bact_s3 * vol < 0.5) {
      grow3 <- 0
    } else {
      grow3 <- 1
    }
    if (bact_s4 * vol < 0.5) {
      grow4 <- 0
    } else {
      grow4 <- 1
    }

    # ---- 12. Bacterial ODEs (Eqs. 6-8 and 11) ---------------------------
    d/dt(bact_s1) <- grow1 * kgs * bact_s1 - (kd + kTot1) * bact_s1 -
      ksrS * bact_s1 + kncsS * bact_nc1 - ksncS * bact_s1 * flagNc
    d/dt(bact_r1) <- -kd * bact_r1 + ksrS * bact_s1
    d/dt(bact_nc1) <- ksncS * bact_s1 * flagNc - kncsS * bact_nc1 - (kd + kTot1) * bact_nc1

    d/dt(bact_s2) <- grow2 * kgs * bact_s2 - (kd + kTot2) * bact_s2 -
      ksrS * bact_s2 + kncsS * bact_nc2 - ksncS * bact_s2 * flagNc
    d/dt(bact_r2) <- -kd * bact_r2 + ksrS * bact_s2
    d/dt(bact_nc2) <- ksncS * bact_s2 * flagNc - kncsS * bact_nc2 - (kd + kTot2) * bact_nc2

    d/dt(bact_s3) <- grow3 * kgr * bact_s3 - (kd + kTot3) * bact_s3 -
      ksrR * bact_s3 + kncsR * bact_nc3 - ksncR * bact_s3 * flagNc
    d/dt(bact_r3) <- -kd * bact_r3 + ksrR * bact_s3
    d/dt(bact_nc3) <- ksncR * bact_s3 * flagNc - kncsR * bact_nc3 - (kd + kTot3) * bact_nc3

    d/dt(bact_s4) <- grow4 * kgr * bact_s4 - (kd + kTot4) * bact_s4 -
      ksrR * bact_s4 + kncsR * bact_nc4 - ksncR * bact_s4 * flagNc
    d/dt(bact_r4) <- -kd * bact_r4 + ksrR * bact_s4
    d/dt(bact_nc4) <- ksncR * bact_s4 * flagNc - kncsR * bact_nc4 - (kd + kTot4) * bact_nc4

    # Pre-existing resting bacteria neither multiply nor move between states.
    d/dt(bact_pr) <- -(kd + kTotPr) * bact_pr

    # ---- 13. Colistin labware binding (Eqs. 14 and 15) ------------------
    kbin <- kbinMax - kbinMax * cst_bound^hillBin / (b50^hillBin + cst_bound^hillBin)
    d/dt(cst_unbound) <- -cst_unbound * kbin + kunb * cst_bound
    d/dt(cst_bound) <- cst_unbound * kbin - kunb * cst_bound

    # ---- 14. Initial conditions (Eqs. 1-5) ------------------------------
    bact_s1(0) <- inoc * (1 - fmutCip) * (1 - fmutCst) * (1 - fmutPr)
    bact_s2(0) <- inoc * (1 - fmutCip) * fmutCst * (1 - fmutPr)
    bact_s3(0) <- inoc * fmutCip * (1 - fmutCst) * (1 - fmutPr)
    bact_s4(0) <- inoc * fmutCip * fmutCst * (1 - fmutPr)
    bact_pr(0) <- inoc * fmutPr
    cst_unbound(0) <- cstNominal * fu0
    cst_bound(0) <- 0

    # ---- 15. Output (log10 CFU/mL) --------------------------------------
    # $ERROR ATOT sums the S and R states of all four subpopulations plus the
    # pre-existing resting subpopulation. The Nc states are excluded: those
    # bacteria are viable but do not form colonies, so they are invisible to
    # the plate count that the model is fitted to.
    cfuObserved <- bact_s1 + bact_r1 + bact_s2 + bact_r2 +
      bact_s3 + bact_r3 + bact_s4 + bact_r4 + bact_pr
    # max() rather than a bare "+ eps": once a subpopulation has been driven to
    # ~1e-70 CFU/mL the solver can return a small negative excursion, which
    # would make log10() undefined.
    Cc <- log10(max(cfuObserved, eps))
    Cc ~ add(addSd)
  })
}
