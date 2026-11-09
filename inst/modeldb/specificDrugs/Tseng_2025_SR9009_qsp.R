Tseng_2025_SR9009_qsp <- function() {
  description <- "QSP. Theoretical (in silico) lung-cell circadian clock model coupled to SARS-CoV-2 entry and replication, with the REV-ERB agonist SR9009 as a timed perturbation (Tseng 2025). Twelve clock ODEs (Per/Cry, Clock/Bmal1, Dbp, Ror and Rev-Erb mRNAs and proteins, plus the PER:CRY and CLOCK:BMAL1 complexes) with Hill-type transcriptional activation and repression and first-order degradation give a free-running 24 h oscillator. Four infection ODEs couple to it: ACE2 (induced by CLOCK:BMAL1, repressed by PER:CRY) gates entry of free virus into the endosome, endosomal virus cycles with replicated virus, and interferon-stimulated genes (ISG; repressed by CLOCK:BMAL1, induced by PER:CRY) repress replication. SR9009 adds a zero-order term to Rev-Erb transcription while a dose is on; a cortisol pulse adds to Per/Cry transcription and blocks Rev-Erb transcription. Inoculation, SR9009 and cortisol exposures are given as events: a dose of 1 into vfree inoculates, and a dose of 1 into sr9009 or cortisol switches the exposure on until a replacement (evid = 5, amt = 0) event switches it off. All levels are in arbitrary units and time in hours. Deterministic: no between-subject variability and no residual error."
  reference <- paste(
    "Tseng Y-Y. (2025). A theoretical systems chronopharmacology approach for",
    "COVID-19: Modeling circadian regulation of lung infection and potential",
    "precision therapies. CPT Pharmacometrics Syst Pharmacol 14(2):340-350.",
    "doi:10.1002/psp4.13277. Equations from Data S1, parameter values from",
    "Table S1, and initial conditions, inoculum, SR9009 and cortisol switch",
    "values from the SBML model deposited as File S1."
  )
  vignette <- "Tseng_2025_SR9009"

  # Paper-mechanistic states with no canonical analogue in
  # inst/references/compartment-names.md: the clock mRNAs and proteins, the
  # PER:CRY and CLOCK:BMAL1 complexes, the infection layer (ACE2, free,
  # endosomal and replicated virus, ISG), and two exposure switches that are
  # states only so that SR9009 and cortisol pulses can be given as events
  # (the d/dt = 0 switch pattern of Bouhaddou_2020_iadademstat_invitro.R).
  paper_specific_compartments <- c(
    "cb_mrna",
    "pc_mrna",
    "dbp_mrna",
    "ror_mrna",
    "rev_mrna",
    "cb_prot",
    "pc_prot",
    "rev_prot",
    "ror_prot",
    "dbp_prot",
    "cbc",
    "pcc",
    "ace2",
    "vfree",
    "vendo",
    "vrep",
    "isg",
    "sr9009",
    "cortisol"
  )

  dosing <- c("vfree", "sr9009", "cortisol")

  units <- list(
    time = "h",
    dosing = "unitless (exposure switch level; 1 = on) for sr9009 and cortisol; a.u. for vfree",
    concentration = "a.u. (relative level) for every state"
  )

  compartmentData <- list(
    cb_mrna = list(analyte = "Clock/Bmal1 mRNA (cb_RNA)", units = "a.u.", specimen = "tissue", verified = TRUE),
    pc_mrna = list(analyte = "Per/Cry mRNA (pc_RNA)", units = "a.u.", specimen = "tissue", verified = TRUE),
    dbp_mrna = list(analyte = "Dbp mRNA (dbp_RNA)", units = "a.u.", specimen = "tissue", verified = TRUE),
    ror_mrna = list(analyte = "Ror mRNA (ror_RNA)", units = "a.u.", specimen = "tissue", verified = TRUE),
    rev_mrna = list(analyte = "Rev-Erb mRNA (rev_RNA)", units = "a.u.", specimen = "tissue", verified = TRUE),
    cb_prot = list(analyte = "CLOCK/BMAL1 protein (CB)", units = "a.u.", specimen = "tissue", verified = TRUE),
    pc_prot = list(analyte = "PER/CRY protein (PC)", units = "a.u.", specimen = "tissue", verified = TRUE),
    rev_prot = list(analyte = "REV-ERB protein (REV)", units = "a.u.", specimen = "tissue", verified = TRUE),
    ror_prot = list(analyte = "ROR protein (ROR)", units = "a.u.", specimen = "tissue", verified = TRUE),
    dbp_prot = list(analyte = "DBP protein (DBP)", units = "a.u.", specimen = "tissue", verified = TRUE),
    cbc = list(analyte = "CLOCK:BMAL1 complex (CBC)", units = "a.u.", specimen = "tissue", verified = TRUE),
    pcc = list(analyte = "PER:CRY complex (PCC)", units = "a.u.", specimen = "tissue", verified = TRUE),
    ace2 = list(analyte = "ACE2 receptor protein", units = "a.u.", specimen = "tissue", verified = TRUE),
    vfree = list(
      analyte = "free SARS-CoV-2 inoculum outside the cell (Vfree)",
      units = "a.u.",
      specimen = "not applicable",
      verified = TRUE
    ),
    vendo = list(
      analyte = "SARS-CoV-2 in the endosome (Vendo)",
      units = "a.u.",
      specimen = "endosome",
      verified = TRUE
    ),
    vrep = list(analyte = "replicated SARS-CoV-2 (Vrep)", units = "a.u.", specimen = "tissue", verified = TRUE),
    isg = list(
      analyte = "interferon-stimulated gene products (ISG)",
      units = "a.u.",
      specimen = "tissue",
      verified = TRUE
    ),
    sr9009 = list(
      analyte = "SR9009 exposure switch (Ka_rev_SR9009 on/off)",
      units = "unitless (1 = on)",
      specimen = "not applicable",
      verified = TRUE
    ),
    cortisol = list(
      analyte = "cortisol pulse switch (kact_pc_cortisol / Krep_rev_cortisol on/off)",
      units = "unitless (1 = on)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species = "human (in silico lung epithelial cell)",
    n_subjects = NA_integer_,
    n_studies = 0L,
    disease_state = "SARS-CoV-2 infection of lung cells (theoretical model; no subjects)",
    dose_range = "SR9009 as a 1 h switch-on of Ka_rev_SR9009 = 0.5 once daily for 5 days, started at circadian times 0-24 on the day after inoculation; cortisol as single 1 h pulses at 4 h intervals over one circadian day",
    notes = "Theoretical systems-pharmacology model built in CellDesigner 4.4.2 and simulated with SOSlib; no subject data were fitted. mRNA degradation rates were taken from experimental data or published models and the remaining parameters and initial conditions were tuned by hand so that the clock has a 24 h period and the mRNA phases fall inside the experimental windows of Figure 2b. The infection layer was qualitatively matched to the BMAL1 / ACE2 / ISG findings in Calu-3 lung epithelial cells (Zhuang 2021, the paper's ref. 9). SR9009 has no PK model: its exposure is the on/off level of a transcription term."
  )

  ini({
    # ==================================================================
    # Clock transcription (steps 1-5) -- Data S1 Eqs. 1-5, Table S1.
    # ==================================================================
    ktx_cb <- 0.55
    label("Clock/Bmal1 (cb) maximal transcription rate k_cb (a.u./h)") # Table S1 k_cb = 0.55
    ktx_dbp <- 1.24
    label("Dbp maximal transcription rate k_dbp (a.u./h)") # Table S1 k_dbp = 1.24
    ktx_pc <- 3.9
    label("Per/Cry (pc) maximal transcription rate k_pc (a.u./h)") # Table S1 k_pc = 3.9
    ktx_rev <- 1.15
    label("Rev-Erb maximal transcription rate k_rev (a.u./h)") # Table S1 k_rev = 1.15
    ktx_ror <- 0.29
    label("Ror maximal transcription rate k_ror (a.u./h)") # Table S1 k_ror = 0.29

    # ==================================================================
    # Translation and complex formation (steps 6-12) -- Eqs. 6-12.
    # ==================================================================
    ktl_cb <- 0.6
    label("CLOCK/BMAL1 translation rate constant k_CB (1/h)") # Table S1 k_CB = 0.6
    ktl_dbp <- 0.7
    label("DBP translation rate constant k_DBP (1/h)") # Table S1 k_DBP = 0.7
    ktl_pc <- 0.5
    label("PER/CRY translation rate constant k_PC (1/h)") # Table S1 k_PC = 0.5
    ktl_rev <- 0.6
    label("REV-ERB translation rate constant k_REV (1/h)") # Table S1 k_REV = 0.6
    ktl_ror <- 0.3
    label("ROR translation rate constant k_ROR (1/h)") # Table S1 k_ROR = 0.3
    kf_cbc <- 0.3
    label("CLOCK:BMAL1 complex formation rate constant k_CBC (1/h)") # Table S1 k_CBC = 0.3
    kf_pcc <- 0.5
    label("PER:CRY complex formation rate constant k_PCC (1/h)") # Table S1 k_PCC = 0.5

    # ==================================================================
    # Transcriptional regulation thresholds -- Table S1 'Regulation
    # thresholds'. ka_* activate, ki_* repress.
    # ==================================================================
    ka_cb_ror <- 0.85
    label("Activation threshold of cb transcription by ROR Ka_cb_ROR (a.u.)") # Table S1 Ka_cb_ROR = 0.85
    ka_dbp_cbc <- 1
    label("Activation threshold of Dbp transcription by CLOCK:BMAL1 Ka_dbp_CBC (a.u.)") # Table S1 Ka_dbp_CBC = 1
    ka_pc_cbc <- 0.3
    label("Activation threshold of pc transcription by CLOCK:BMAL1 Ka_pc_CBC (a.u.)") # Table S1 Ka_pc_CBC = 0.3
    ka_pc_dbp <- 0.5
    label("Activation threshold of pc transcription by DBP Ka_pc_DBP (a.u.)") # Table S1 Ka_pc_DBP = 0.5
    ka_pc_ror <- 1.85
    label("Activation threshold of pc transcription by ROR Ka_pc_ROR (a.u.)") # Table S1 Ka_pc_ROR = 1.85
    ka_rev_cbc <- 0.85
    label("Activation threshold of Rev-Erb transcription by CLOCK:BMAL1 Ka_rev_CBC (a.u.)") # Table S1 Ka_rev_CBC = 0.85
    ka_rev_dbp <- 0.1
    label("Activation threshold of Rev-Erb transcription by DBP Ka_rev_DBP (a.u.)") # Table S1 Ka_rev_DBP = 0.1
    ka_ror_dbp <- 0.33
    label("Activation threshold of Ror transcription by DBP Ka_ror_DBP (a.u.)") # Table S1 Ka_ror_DBP = 0.33
    ki_cb_rev <- 0.6
    label("Repression threshold of cb transcription by REV-ERB Ki_cb_REV (a.u.)") # Table S1 Ki_cb_REV = 0.6
    ki_dbp_pcc <- 0.6
    label("Repression threshold of Dbp transcription by PER:CRY Ki_dbp_PCC (a.u.)") # Table S1 Ki_dbp_PCC = 0.6
    ki_pc_pcc <- 1
    label("Repression threshold of pc transcription by PER:CRY Ki_pc_PCC (a.u.)") # Table S1 Ki_pc_PCC = 1
    ki_pc_rev <- 4
    label("Repression threshold of pc transcription by REV-ERB Ki_pc_REV (a.u.)") # Table S1 Ki_pc_REV = 4
    ki_rev_pcc <- 0.8
    label("Repression threshold of Rev-Erb transcription by PER:CRY Ki_rev_PCC (a.u.)") # Table S1 Ki_rev_PCC = 0.8

    # ==================================================================
    # Hill coefficients of the clock regulations -- Table S1.
    # ==================================================================
    hill_cb_rev <- 4
    label("Hill coefficient, cb repression by REV-ERB H_cb_REV (unitless)") # Table S1 H_cb_REV = 4
    hill_cb_ror <- 2
    label("Hill coefficient, cb activation by ROR H_cb_ROR (unitless)") # Table S1 H_cb_ROR = 2
    hill_dbp_cbc <- 4
    label("Hill coefficient, Dbp activation by CLOCK:BMAL1 H_dbp_CBC (unitless)") # Table S1 H_dbp_CBC = 4
    hill_dbp_pcc <- 2
    label("Hill coefficient, Dbp repression by PER:CRY H_dbp_PCC (unitless)") # Table S1 H_dbp_PCC = 2
    hill_pc_cbc <- 2
    label("Hill coefficient, pc activation by CLOCK:BMAL1 H_pc_CBC (unitless)") # Table S1 H_pc_CBC = 2
    hill_pc_dbp <- 2
    label("Hill coefficient, pc activation by DBP H_pc_DBP (unitless)") # Table S1 H_pc_DBP = 2
    hill_pc_pcc <- 4
    label("Hill coefficient, pc repression by PER:CRY H_pc_PCC (unitless)") # Table S1 H_pc_PCC = 4
    hill_pc_rev <- 4
    label("Hill coefficient, pc repression by REV-ERB H_pc_REV (unitless)") # Table S1 H_pc_REV = 4
    hill_pc_ror <- 4
    label("Hill coefficient, pc activation by ROR H_pc_ROR (unitless)") # Table S1 H_pc_ROR = 4
    hill_rev_cbc <- 4
    label("Hill coefficient, Rev-Erb activation by CLOCK:BMAL1 H_rev_CBC (unitless)") # Table S1 H_rev_CBC = 4
    hill_rev_dbp <- 2
    label("Hill coefficient, Rev-Erb activation by DBP H_rev_DBP (unitless)") # Table S1 H_rev_DBP = 2
    hill_rev_pcc <- 2
    label("Hill coefficient, Rev-Erb repression by PER:CRY H_rev_PCC (unitless)") # Table S1 H_rev_PCC = 2
    hill_ror_dbp <- 2
    label("Hill coefficient, Ror activation by DBP H_ror_DBP (unitless)") # Table S1 H_ror_DBP = 2

    # ==================================================================
    # First-order degradation (steps 13-24) -- Table S1. The mRNA
    # degradation rates were set from experimental data or published
    # models rather than tuned (Methods, Model construction).
    # ==================================================================
    kdeg_cb_mrna <- fixed(0.37)
    label("cb mRNA degradation rate constant kd_cb (1/h)") # Table S1 kd_cb = 0.37
    kdeg_dbp_mrna <- fixed(0.36)
    label("Dbp mRNA degradation rate constant kd_dbp (1/h)") # Table S1 kd_dbp = 0.36
    kdeg_pc_mrna <- fixed(0.3)
    label("pc mRNA degradation rate constant kd_pc (1/h)") # Table S1 kd_pc = 0.3
    kdeg_rev_mrna <- fixed(0.55)
    label("Rev-Erb mRNA degradation rate constant kd_rev (1/h)") # Table S1 kd_rev = 0.55
    kdeg_ror_mrna <- fixed(0.31)
    label("Ror mRNA degradation rate constant kd_ror (1/h)") # Table S1 kd_ror = 0.31
    kdeg_cb_prot <- 0.2
    label("CLOCK/BMAL1 protein degradation rate constant kd_CB (1/h)") # Table S1 kd_CB = 0.2
    kdeg_dbp_prot <- 0.3
    label("DBP protein degradation rate constant kd_DBP (1/h)") # Table S1 kd_DBP = 0.3
    kdeg_pc_prot <- 0.3
    label("PER/CRY protein degradation rate constant kd_PC (1/h)") # Table S1 kd_PC = 0.3
    kdeg_rev_prot <- 0.25
    label("REV-ERB protein degradation rate constant kd_REV (1/h)") # Table S1 kd_REV = 0.25
    kdeg_ror_prot <- 0.18
    label("ROR protein degradation rate constant kd_ROR (1/h)") # Table S1 kd_ROR = 0.18
    kdeg_cbc <- 0.27
    label("CLOCK:BMAL1 complex degradation rate constant kd_CBC (1/h)") # Table S1 kd_CBC = 0.27
    kdeg_pcc <- 0.3
    label("PER:CRY complex degradation rate constant kd_PCC (1/h)") # Table S1 kd_PCC = 0.3

    # ==================================================================
    # ACE2, ISG and virus (steps 25-32) -- Data S1 Eqs. 13-16, Table S1.
    # ==================================================================
    ktl_ace2 <- 0.25
    label("ACE2 translation rate constant k_ACE2 (1/h)") # Table S1 k_ACE2 = 0.25
    k_isg <- 0.3
    label("ISG production rate constant from precursor k_ISG (1/h)") # Table S1 k_ISG = 0.3
    k_vendo <- 0.5
    label("Return rate constant of replicated virus to the endosome k_Vendo (1/h)") # Table S1 k_Vendo = 0.5
    k_vrep <- 5
    label("Maximal virus replication rate constant k_Vrep (1/h)") # Table S1 k_Vrep = 5
    kin_vfree <- 0.3
    label("Maximal virus entry rate constant kin_Vfree (1/h)") # File S1 SBML global parameter kin_Vfree = 0.3 (not in Table S1)
    ka_ace2_cbc <- 0.3
    label("Activation threshold of ACE2 by CLOCK:BMAL1 Ka_ACE2_CBC (a.u.)") # Table S1 Ka_ACE2_CB = 0.3
    ka_isg_pcc <- 0.01
    label("Activation threshold of ISG by PER:CRY Ka_ISG_PCC (a.u.)") # Table S1 Ka_ISG_PC = 0.01
    ka_vin_ace2 <- 0.3
    label("Activation threshold of virus entry by ACE2 Ka_Vin_ACE2 (a.u.)") # Table S1 Ka_Vin_ACE2 = 0.3
    ki_ace2_pcc <- 2.4
    label("Repression threshold of ACE2 by PER:CRY Ki_ACE2_PCC (a.u.)") # Table S1 Ki_ACE2_PCC = 2.4
    ki_isg_cbc <- 1
    label("Repression threshold of ISG by CLOCK:BMAL1 Ki_ISG_CBC (a.u.)") # Table S1 Ki_ISG_CBC = 1
    ki_vrep_isg <- 1
    label("Repression threshold of virus replication by ISG Ki_Vrep_ISG (a.u.)") # Table S1 Ki_Vrep_ISG = 1
    hill_ace2_cbc <- 4
    label("Hill coefficient, ACE2 activation by CLOCK:BMAL1 H_ACE2_CBC (unitless)") # Table S1 H_ACE2_CBC = 4
    hill_ace2_pcc <- 2.01
    label("Hill coefficient, ACE2 repression by PER:CRY H_ACE2_PCC (unitless)") # Table S1 H_ACE2_PCC = 2.01
    hill_isg_cbc <- 2
    label("Hill coefficient, ISG repression by CLOCK:BMAL1 H_ISG_CBC (unitless)") # Table S1 H_ISG_CBC = 2
    hill_isg_pcc <- 2
    label("Hill coefficient, ISG activation by PER:CRY H_ISG_PCC (unitless)") # Table S1 H_ISG_PCC = 2
    hill_vin_ace2 <- 2
    label("Hill coefficient, virus entry activation by ACE2 H_Vin_ACE2 (unitless)") # Table S1 H_Vin_ACE2 = 2
    hill_vrep_isg <- 2
    label("Hill coefficient, virus replication repression by ISG H_Vrep_ISG (unitless)") # Table S1 H_Vrep_ISG = 2
    kdeg_ace2 <- 0.3
    label("ACE2 degradation rate constant kd_ACE2 (1/h)") # Table S1 kd_ACE2 = 0.3
    kdeg_isg <- 0.3
    label("ISG degradation rate constant kd_ISG (1/h)") # Table S1 kd_ISG = 0.3
    kdeg_vendo <- 3.164
    label("Endosomal virus degradation rate constant kd_Vendo (1/h)") # Table S1 kd_Vendo = 3.164
    ace2_mrna <- fixed(1)
    label("Constant Ace2 mRNA level ace2_RNA (a.u.)") # File S1 SBML species ace2_RNA initialAmount = 1, boundaryCondition = true
    isg_pre <- fixed(1)
    label("Constant ISG precursor level ISGpre (a.u.)") # File S1 SBML species ISGpre initialAmount = 1, boundaryCondition = true

    # ==================================================================
    # Exposure effects while a switch state is 1 -- File S1 events.
    # ==================================================================
    ka_rev_sr9009 <- 0.5
    label("Zero-order Rev-Erb transcription added by SR9009 Ka_rev_SR9009 (a.u./h)") # File S1 Ka_rev_SR9009_on = 0.5; Figure 5b right axis pulses to 0.5
    kact_pc_cortisol <- 10
    label("Zero-order Per/Cry transcription added by a cortisol pulse kact_pc_cortisol (a.u./h)") # Methods 'from 0 to 10'; File S1 kact_pc_cortisol_on = 10
    krep_rev_cortisol <- 0
    label("Rev-Erb transcription multiplier during a cortisol pulse Krep_rev_cortisol (unitless)") # Methods 'altered from 1 to 0'; File S1 Krep_rev_cortisol_on = 0
  })

  model({
    # Hill activation x^h / (x^h + Ka^h) and repression Ki^h / (x^h + Ki^h),
    # written out per regulation as in Data S1.
    act_cb_ror <- ror_prot^hill_cb_ror / (ror_prot^hill_cb_ror + ka_cb_ror^hill_cb_ror)
    rep_cb_rev <- ki_cb_rev^hill_cb_rev / (rev_prot^hill_cb_rev + ki_cb_rev^hill_cb_rev)

    act_pc_ror <- ror_prot^hill_pc_ror / (ror_prot^hill_pc_ror + ka_pc_ror^hill_pc_ror)
    rep_pc_pcc <- ki_pc_pcc^hill_pc_pcc / (pcc^hill_pc_pcc + ki_pc_pcc^hill_pc_pcc)
    act_pc_dbp <- dbp_prot^hill_pc_dbp / (dbp_prot^hill_pc_dbp + ka_pc_dbp^hill_pc_dbp)
    act_pc_cbc <- cbc^hill_pc_cbc / (cbc^hill_pc_cbc + ka_pc_cbc^hill_pc_cbc)
    rep_pc_rev <- ki_pc_rev^hill_pc_rev / (rev_prot^hill_pc_rev + ki_pc_rev^hill_pc_rev)

    act_dbp_cbc <- cbc^hill_dbp_cbc / (cbc^hill_dbp_cbc + ka_dbp_cbc^hill_dbp_cbc)
    rep_dbp_pcc <- ki_dbp_pcc^hill_dbp_pcc / (pcc^hill_dbp_pcc + ki_dbp_pcc^hill_dbp_pcc)

    act_ror_dbp <- dbp_prot^hill_ror_dbp / (dbp_prot^hill_ror_dbp + ka_ror_dbp^hill_ror_dbp)

    act_rev_cbc <- cbc^hill_rev_cbc / (cbc^hill_rev_cbc + ka_rev_cbc^hill_rev_cbc)
    rep_rev_pcc <- ki_rev_pcc^hill_rev_pcc / (pcc^hill_rev_pcc + ki_rev_pcc^hill_rev_pcc)
    act_rev_dbp <- dbp_prot^hill_rev_dbp / (dbp_prot^hill_rev_dbp + ka_rev_dbp^hill_rev_dbp)

    # Exposure terms. The switches hold their level between events (d/dt = 0);
    # a cortisol pulse moves Krep_rev_cortisol from 1 to krep_rev_cortisol.
    d/dt(sr9009) <- 0
    d/dt(cortisol) <- 0
    sr9009_rev_tx <- ka_rev_sr9009 * sr9009
    cortisol_pc_tx <- kact_pc_cortisol * cortisol
    cortisol_rev_rep <- 1 - cortisol * (1 - krep_rev_cortisol)

    # Clock mRNAs, Data S1 Eqs. 1-5.
    d/dt(cb_mrna) <- ktx_cb * act_cb_ror * rep_cb_rev - kdeg_cb_mrna * cb_mrna
    d/dt(pc_mrna) <- ktx_pc * act_pc_ror * rep_pc_pcc * act_pc_dbp * act_pc_cbc * rep_pc_rev +
      cortisol_pc_tx - kdeg_pc_mrna * pc_mrna
    d/dt(dbp_mrna) <- ktx_dbp * act_dbp_cbc * rep_dbp_pcc - kdeg_dbp_mrna * dbp_mrna
    d/dt(ror_mrna) <- ktx_ror * act_ror_dbp - kdeg_ror_mrna * ror_mrna
    d/dt(rev_mrna) <- ktx_rev * act_rev_cbc * rep_rev_pcc * act_rev_dbp * cortisol_rev_rep +
      sr9009_rev_tx - kdeg_rev_mrna * rev_mrna

    # Clock proteins and complexes, Data S1 Eqs. 6-12.
    d/dt(cb_prot) <- ktl_cb * cb_mrna - kf_cbc * cb_prot - kdeg_cb_prot * cb_prot
    d/dt(pc_prot) <- ktl_pc * pc_mrna - kf_pcc * pc_prot - kdeg_pc_prot * pc_prot
    d/dt(rev_prot) <- ktl_rev * rev_mrna - kdeg_rev_prot * rev_prot
    d/dt(ror_prot) <- ktl_ror * ror_mrna - kdeg_ror_prot * ror_prot
    d/dt(dbp_prot) <- ktl_dbp * dbp_mrna - kdeg_dbp_prot * dbp_prot
    d/dt(cbc) <- kf_cbc * cb_prot - kdeg_cbc * cbc
    d/dt(pcc) <- kf_pcc * pc_prot - kdeg_pcc * pcc

    # Infection layer, Data S1 Eqs. 13-16. Free virus is an input that entry
    # does not deplete: Data S1 adds four ODEs and uses [Vfree] only as a
    # driver, and Figure 5a shows a sustained infection.
    act_ace2_cbc <- cbc^hill_ace2_cbc / (cbc^hill_ace2_cbc + ka_ace2_cbc^hill_ace2_cbc)
    rep_ace2_pcc <- ki_ace2_pcc^hill_ace2_pcc / (pcc^hill_ace2_pcc + ki_ace2_pcc^hill_ace2_pcc)
    act_vin_ace2 <- ace2^hill_vin_ace2 / (ace2^hill_vin_ace2 + ka_vin_ace2^hill_vin_ace2)
    rep_vrep_isg <- ki_vrep_isg^hill_vrep_isg / (isg^hill_vrep_isg + ki_vrep_isg^hill_vrep_isg)
    rep_isg_cbc <- ki_isg_cbc^hill_isg_cbc / (cbc^hill_isg_cbc + ki_isg_cbc^hill_isg_cbc)
    act_isg_pcc <- pcc^hill_isg_pcc / (pcc^hill_isg_pcc + ka_isg_pcc^hill_isg_pcc)

    virus_entry <- vfree * kin_vfree * act_vin_ace2
    virus_replication <- vendo * k_vrep * rep_vrep_isg

    d/dt(ace2) <- ace2_mrna * ktl_ace2 * rep_ace2_pcc * act_ace2_cbc - kdeg_ace2 * ace2
    d/dt(vfree) <- 0
    d/dt(vendo) <- virus_entry + k_vendo * vrep - virus_replication - kdeg_vendo * vendo
    d/dt(vrep) <- virus_replication - k_vendo * vrep
    d/dt(isg) <- isg_pre * k_isg * rep_isg_cbc * act_isg_pcc - kdeg_isg * isg

    # Initial conditions: File S1 SBML species initialAmount (a point on the
    # free-running limit cycle; virus states start empty).
    cb_mrna(0) <- 0.8379059396
    pc_mrna(0) <- 0.3302851508
    dbp_mrna(0) <- 0.0337391156
    ror_mrna(0) <- 0.498005573
    rev_mrna(0) <- 0.038319912
    cb_prot(0) <- 0.7789134637
    pc_prot(0) <- 0.248584977
    rev_prot(0) <- 0.1932572563
    ror_prot(0) <- 1.194949787
    dbp_prot(0) <- 0.190787271
    cbc(0) <- 0.477004444
    pcc(0) <- 0.591785229
    ace2(0) <- 0.622120294388117
    isg(0) <- 0.762080921516905

    # Viral load plotted in Figure 5 and Figure S8 (File S1 assignment rule
    # total_virus = Vendo + Vrep).
    total_virus <- vendo + vrep
  })
}
