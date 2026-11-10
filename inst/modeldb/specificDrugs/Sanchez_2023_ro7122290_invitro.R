Sanchez_2023_ro7122290_invitro <- function() {
  description <- "QSP. In vitro (MKN-45 tumor cell + PBMC + FAP-expressing fibroblast co-culture). FAP-4-1BBL (RO7122290) combined with the CEA x CD3 T-cell bispecific cibisatamab: Model 1 (cibisatamab-driven 4-1BB expression per CD8+ T cell through a 7-compartment transit chain, nonlinear mixed-effects fit), the FAP-4-1BBL / FAP / 4-1BB binding and trimeric-complex model in the culture well, Model 3 (sigmoid Emax fold increase in tumor cell killing versus the time-averaged trimeric complexes per T cell) and Model 2 (generalized logistic decrease of cibisatamab monotherapy killing with the proportion of FAP-positive fibroblasts) (Sanchez 2023)."
  reference <- "Sanchez J, Claus C, Albrecht R, Gaillard BC, Marinho J, McIntyre C, Tanos T, Boehnke A, Friberg LE, Jonsson S, Frances N. A model-based approach leveraging in vitro data to support dose selection from the outset: A framework for bispecific antibodies in immuno-oncology. CPT Pharmacometrics Syst Pharmacol. 2023;12(11):1804-1818. doi:10.1002/psp4.13065"
  vignette <- "Sanchez_2023_ro7122290_cibisatamab"
  units <- list(
    time = "h",
    dosing = "nM",
    concentration = "nM",
    expr_41bb = "receptors/cell",
    trimer_per_tcell = "trimeric complexes per T cell",
    fold_kill = "fold increase over cibisatamab monotherapy",
    kill_mono_pct = "% of cibisatamab monotherapy killing without FAP-positive fibroblasts"
  )

  # The free-drug pool, the receptor pools and their complexes are
  # concentrations in the culture well; the transit chain and expr_41bb are
  # receptor counts per CD8+ T cell; auc_trimer integrates trimer_per_tcell so
  # that Model 3 can read the time-averaged value.
  paper_specific_compartments <- c(
    "expr_41bb",
    "target_41bb",
    "target_fap",
    "complex_fap",
    "complex_41bb",
    "medium",
    "auc_trimer"
  )

  covariateData <- list(
    CONC_CIBISATAMAB_NM = list(
      description = "Static cibisatamab concentration applied to the co-culture well",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the sigmoid Emax Signal of Model 1 (Equation 2). The Monolix regressor 'CEA' of the deposited Model 1 code. Experiment 1 tested 0, 0.2, 2, 20, 200 and 2000 nM; Experiments 2 and 3 used 2 nM and 0.2 nM.",
      source_name = "CEA"
    )
  )

  compartmentData <- list(
    transit1 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit4 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit5 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit6 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit7 = list(
      analyte = "4-1BB induction signal",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    expr_41bb = list(
      analyte = "4-1BB receptor",
      units = "receptors/cell",
      specimen = "not applicable",
      verified = TRUE
    ),
    target_41bb = list(analyte = "free 4-1BB receptor", units = "nM", specimen = "not applicable", verified = TRUE),
    medium = list(analyte = "free FAP-4-1BBL (RO7122290)", units = "nM", specimen = "not applicable", verified = TRUE),
    target_fap = list(analyte = "free FAP receptor", units = "nM", specimen = "not applicable", verified = TRUE),
    complex_fap = list(analyte = "FAP-4-1BBL-FAP dimer", units = "nM", specimen = "not applicable", verified = TRUE),
    complex_41bb = list(analyte = "FAP-4-1BBL-4-1BB dimer", units = "nM", specimen = "not applicable", verified = TRUE),
    trimer = list(
      analyte = "FAP-FAP-4-1BBL-4-1BB trimeric complex",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_trimer = list(
      analyte = "cumulative trimeric complexes per T cell",
      units = "trimeric complexes per T cell * h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (MKN-45 gastric tumor cells, human PBMC, FAP-expressing fibroblasts)",
    n_subjects = NA_integer_,
    n_studies = 3L,
    disease_state = "CEACAM5-expressing gastric cancer cell line co-cultured with PBMC effector cells (effector:target 1:1) and fibroblasts",
    dose_range = "Cibisatamab 0.2-2000 nM; FAP-4-1BBL 0.001-10000 nM (static applied concentrations)",
    notes = paste(
      "Experiment 1 (Model 1): 10,000 MKN-45 cells, 10,000 PBMC and 10,000 NIH-3T3 fibroblasts at 92,000 FAP/fibroblast in 200 uL;",
      "cibisatamab 0, 0.2, 2, 20, 200, 2000 nM, or 20 nM + 1 nM FAP-4-1BBL; 4-1BB per CD8+ T cell by flow cytometry at 4, 21, 45 and 69 h,",
      "median of three technical replicates.",
      "Experiment 2 (Model 3): 5,000 MKN-45-NLR cells, 5,000 PBMC, 10,000 fibroblasts at 92,000 (high) or 12,700 (low) FAP/fibroblast in 200 uL;",
      "2 nM cibisatamab plus FAP-4-1BBL 0-10,000 nM; IncuCyte imaging over 138 h.",
      "Experiment 3 (Model 2): 5,000 MKN-45-NLR cells, 5,000 PBMC, 10,000 hTERT fibroblasts of which 0-100% expressed 6,800 FAP/fibroblast,",
      "240 uL; 0.2 nM cibisatamab. Supplementary Materials and Methods."
    )
  )

  ini({
    # Model 1 -- 4-1BB expression time course (Table 2; Table S1 for the extra decimals).
    lemax <- log(537.682); label("Maximum Signal for cibisatamab, SignalMAX (receptors/cell)") # Table S1 'Emax 537.682 (12.91)'; Table 2 538
    lec50 <- log(16.004); label("Cibisatamab concentration giving half-maximal Signal (nM)") # Table S1 '16.004 (19.172)'; Table 2 16
    lhill <- log(0.210); label("Hill coefficient of Signal versus cibisatamab concentration (unitless)") # Table 2 and Table S1 0.210
    lsignal2 <- log(381.92); label("Signal in the 20 nM cibisatamab + 1 nM FAP-4-1BBL condition (receptors/cell)") # Table S1 '381.92 (28.11)'; Table 2 382
    lt_stop <- log(16.18); label("Stop time after which Signal and kin are zero (h)") # Table S1 '16.18 (18.85)'; Table 2 16.2
    lktr <- log(0.2801); label("Transit rate constant of the 7-compartment chain (1/h)") # Table S1 '0.2801 (5.2)'; Table 2 0.280
    lkin <- log(2.952); label("Basal 4-1BB synthesis rate (receptors/cell/h)") # Table S1 '2.952 (21.49)'; Table 2 2.95
    lkout <- log(0.2832); label("4-1BB elimination rate constant (1/h)") # Table S1 '0.2832 (19)'; Table 2 0.283
    cibi_signal2 <- fixed(20); label("Cibisatamab concentration of the combination condition that defines Signal2 (nM)") # Table 2 row Signal2; Supplement Experiment 1

    # Monolix lognormal random effects; Table 2 reports their standard deviations.
    etalemax ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors
    etalec50 ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors
    etalhill ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors
    etalsignal2 ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors
    etalt_stop ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors
    etalktr ~ 0.0081 # Table 2 SD 0.09 (31.7% RSE)
    etalkin ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors
    etalkout ~ fixed(0.01) # Table 2 SD 0.1, held constant by the authors

    addSd <- 5.36; label("Additive residual error on 4-1BB receptors per CD8+ T cell (receptors/cell)") # Table 2 'a 5.36 (8.6)'

    # FAP-4-1BBL binding (Table S1).
    kd_fap <- fixed(0.7); label("Equilibrium dissociation constant of FAP-4-1BBL for FAP (nM)") # Table S1 KD,FAP 0.7 (in vitro binding data)
    koff_fap <- fixed(0.1); label("Dissociation rate constant of FAP-4-1BBL from FAP (1/h)") # Table S1 koff,FAP 0.1 (assumed)
    kd_41bb <- fixed(0.2); label("Equilibrium dissociation constant of FAP-4-1BBL for 4-1BB (nM)") # Table S1 KD,4-1BB 0.2 (in vitro binding data)
    koff_41bb <- fixed(0.1); label("Dissociation rate constant of FAP-4-1BBL from 4-1BB (1/h)") # Table S1 koff,4-1BB 0.1 (assumed)

    # Culture-well composition (Supplementary Materials and Methods, Experiment 2).
    vol_well <- fixed(200e-6); label("Culture volume per well (L)") # Supplement Experiments 1-2: 200 uL per well
    n_fib <- fixed(10000); label("Fibroblasts per well (cells)") # Supplement Experiments 1-3: 10,000 fibroblasts
    fap_per_fib <- fixed(92000); label("FAP receptors per FAP-positive fibroblast (receptors/cell)") # Supplement Experiment 1: NIH-3T3 clone 19, 92,000 FAP/fibroblast
    fap_pos_pct <- fixed(100); label("Percentage of fibroblasts expressing FAP (%)") # Supplement Experiments 1-2: all fibroblasts FAP-transfected
    n_tumor <- fixed(5000); label("Tumor cells per well (cells)") # Supplement Experiment 2: 5,000 MKN-45-NLR cells
    et_ratio <- fixed(1); label("Effector (PBMC) to tumor cell ratio (unitless)") # Supplement Experiment 2: E:T 1:1
    tcell_pct <- fixed(40); label("Percentage of effector cells that are T cells (%)") # deposited RxODE code s04 'Tcell_percentage_in_tumor = 40'
    cd8_frac <- fixed(0.4); label("Fraction of T cells carrying the modelled 4-1BB (CD8+) (unitless)") # deposited RxODE code s04 factor 0.4 in the T-cell count

    # Model 3 -- fold increase in killing versus trimeric complexes per T cell (Table 3).
    kill_max <- 4.37; label("Maximum fold increase in tumor cell killing over cibisatamab monotherapy (unitless)") # Table 3 'KillingMAX,combination 4.37 (3.7)'
    trimer50 <- 3.91e-3; label("Trimeric complexes per T cell giving half-maximal fold increase (complexes/cell)") # Table 3 'Trimerics50 3.91e-3 (16.5)'
    hill_kill <- 1.16; label("Hill coefficient of fold increase versus trimeric complexes (unitless)") # Table 3 'HillKill,Combination 1.16 (17.4)'
    kill_base <- 1.04; label("Fold killing without trimeric complex formation (unitless)") # Table 3 'BasalKilling,monotherapy 1.04 (8.2)'

    # Model 2 -- cibisatamab monotherapy killing versus FAP-positive fibroblasts (Table 3).
    kill_nofap <- 101.9; label("Cibisatamab monotherapy killing without FAP-positive fibroblasts (%)") # Table 3 'BasalNo FAP+ 101.9 (3.6)'
    kill_maxdec <- 67.4; label("Maximum decrease in monotherapy killing due to FAP-positive fibroblasts (%)") # Table 3 'Maximum decrease 67.4 (7.8)'
    hill_fap <- 0.86; label("Sigmoid coefficient of the generalized logistic model (unitless)") # Table 3 'SigmoidCoefficient 0.86 (16.8)'
    logit_fap50 <- -1.1; label("Logit FAP-positive proportion at half-maximal decrease (logit units)") # Table 3 'Midpoint -1.1 (21.3)'
  })
  model({
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)
    signal2 <- exp(lsignal2 + etalsignal2)
    t_stop <- exp(lt_stop + etalt_stop)
    ktr <- exp(lktr + etalktr)
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout)
    kon_fap <- koff_fap / kd_fap
    kon_41bb <- koff_41bb / kd_41bb

    # Equation 1: Signal and kin switch off at t_stop.
    stepf <- 0.5 * (tanh(1e15 * (t_stop - t)) + 1)

    # Equation 2 (monotherapy). With FAP-4-1BBL present the deposited code
    # multiplies the monotherapy Signal by an added value, capped at emax; the
    # added value is Signal2 over the monotherapy Signal at the 20 nM
    # combination condition, so the combination arm reproduces Equation 3
    # (Signal = Signal2) at the tested condition.
    sig_mono <- emax * CONC_CIBISATAMAB_NM^hill / (CONC_CIBISATAMAB_NM^hill + ec50^hill)
    sig_ref <- emax * cibi_signal2^hill / (cibi_signal2^hill + ec50^hill)
    added_value <- signal2 / sig_ref
    drug_tot <- medium + complex_fap + complex_41bb + trimer
    sig <- sig_mono
    if (drug_tot > 0) {
      sig <- min(sig_mono * added_value, emax)
    }
    signal <- sig * stepf
    # Equation 4; the deposited Monolix code keeps basal synthesis running in
    # the cibisatamab-free control wells.
    kin_t <- kin
    if (CONC_CIBISATAMAB_NM > 0) {
      kin_t <- kin * stepf
    }

    # T cells per well and the receptors-per-cell to nM conversion (deposited code s04).
    n_tcell <- et_ratio * n_tumor * tcell_pct / 100 * cd8_frac
    conv <- n_tcell / 6.022e23 * 1e9 / vol_well
    # Equation 22 with the plate volume.
    fap0 <- n_fib * fap_pos_pct / 100 * fap_per_fib / 6.022e23 * 1e9 / vol_well

    # Equations 5-8: transit chain and 4-1BB turnover per CD8+ T cell.
    d/dt(transit1) <- ktr * (signal - transit1)
    d/dt(transit2) <- ktr * (transit1 - transit2)
    d/dt(transit3) <- ktr * (transit2 - transit3)
    d/dt(transit4) <- ktr * (transit3 - transit4)
    d/dt(transit5) <- ktr * (transit4 - transit5)
    d/dt(transit6) <- ktr * (transit5 - transit6)
    d/dt(transit7) <- ktr * (transit6 - transit7)
    d/dt(expr_41bb) <- ktr * transit7 - kout * expr_41bb + kin_t
    expr_41bb(0) <- kin / kout

    # Equations 16-21 in the deposited-code form: total 4-1BB in nM follows the
    # per-cell turnover, binding moves it between the free and bound pools.
    d/dt(target_41bb) <- (ktr * transit7 - kout * expr_41bb + kin_t) * conv +
      koff_41bb * (complex_41bb + trimer) - kon_41bb * target_41bb * (medium + complex_fap)
    target_41bb(0) <- kin / kout * conv
    d/dt(medium) <- koff_fap * complex_fap + koff_41bb * complex_41bb -
      kon_fap * target_fap * medium - kon_41bb * target_41bb * medium
    d/dt(complex_fap) <- kon_fap * target_fap * medium + koff_41bb * trimer -
      kon_41bb * target_41bb * complex_fap - koff_fap * complex_fap
    d/dt(complex_41bb) <- kon_41bb * target_41bb * medium + koff_fap * trimer -
      kon_fap * target_fap * complex_41bb - koff_41bb * complex_41bb
    d/dt(trimer) <- kon_41bb * target_41bb * complex_fap + kon_fap * target_fap * complex_41bb -
      (koff_fap + koff_41bb) * trimer
    d/dt(target_fap) <- koff_fap * (trimer + complex_fap) - kon_fap * target_fap * (medium + complex_41bb)
    target_fap(0) <- fap0

    # Equations 23-26.
    ro_fap <- (complex_fap + trimer) / (target_fap + complex_fap + trimer) * 100
    ro_41bb <- (complex_41bb + trimer) / (target_41bb + complex_41bb + trimer) * 100
    trimer_per_tcell <- trimer / conv

    # Model 3 (Equation 27) on the trimeric complexes per T cell averaged from 0 to t.
    d/dt(auc_trimer) <- trimer_per_tcell
    trimer_avg <- trimer_per_tcell
    if (t > 0) {
      trimer_avg <- auc_trimer / t
    }
    fold_kill <- (kill_max - kill_base) * trimer_avg^hill_kill /
      (trimer_avg^hill_kill + trimer50^hill_kill) + kill_base

    # Model 2 (Equation 9); 0% and 100% FAP-positive map to logit -7 and 7.
    p_fap <- fap_pos_pct / 100
    logit_fap <- -7
    if (p_fap >= 1) {
      logit_fap <- 7
    } else if (p_fap > 0) {
      logit_fap <- log(p_fap / (1 - p_fap))
    }
    kill_mono_pct <- kill_nofap - kill_maxdec / (1 + exp(-hill_fap * (logit_fap - logit_fap50)))
    kill_abs <- fold_kill * kill_mono_pct / 100

    expr_41bb ~ add(addSd)
  })
}
