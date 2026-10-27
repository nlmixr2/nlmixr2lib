Gall_2023_fluorouracil_mouse <- function() {
  description <- paste(
    "QSP. Preclinical (mouse, C57BL/6J small-intestinal epithelium). Gall 2023",
    "dynamic model of 5-fluorouracil (5-FU)-induced intestinal toxicity fitted",
    "to mouse data. A concentration-based two-compartment 5-FU plasma PK model",
    "(fitted to Simcyp mouse predictions) feeds a reduced intracellular",
    "metabolic cascade (FUMP -> FUDP -> FUTP, with FdUMP and FdUTP in rapid",
    "equilibrium with FUDP). FUTP drives, through a Hill function, cell-cycle",
    "arrest and apoptosis of proliferative crypt cells in a four-compartment",
    "epithelial model (stem, proliferative progenitors, non-proliferative",
    "crypt and villus cells) with a villus-to-crypt proliferation feedback and",
    "arrested/apoptotic cell compartments in crypt and villus. The epithelial",
    "state also drives plasma citrulline (production proportional to villus",
    "cell number squared, with first-order clearance), the probability of",
    "diarrhea (a logistic function of fractional villus-cell loss) and the",
    "z-scores of six p53-pathway genes (instantaneous linear functions of the",
    "arrested + apoptotic cell fraction). Deterministic typical-value model",
    "(no IIV, no residual error). Time is in hours; the plasma and metabolite",
    "states are concentrations, so an intravenous dose in mg/kg is converted",
    "to its instantaneous plasma-concentration increment through f(central).",
    "Equations S1-S24 of Gall 2023."
  )
  reference <- paste(
    "Gall L, Jardi F, Lammens L, Pinero J, Souza TM, Rodrigues D, Jennen DGJ,",
    "de Kok TM, Coyle L, Chung SW, Ferreira S, Jo H, Beattie KA, Kelly C,",
    "Duckworth CA, Pritchard DM, Pin C. (2023).",
    "A dynamic model of the intestinal epithelium integrates multiple sources",
    "of preclinical data and enables clinical translation of drug-induced",
    "toxicity. CPT Pharmacometrics Syst Pharmacol 12(10):1511-1528.",
    "doi:10.1002/psp4.13029.",
    "Numeric values are from the authors' MATLAB code deposited at BioModels",
    "(MODEL2212120003), which generated the published figures and agrees with",
    "Gall 2023 Tables 1 and 2 to the precision printed there."
  )
  vignette <- "Gall_2023_intestinal_epithelium_toxicity"

  # Mechanistic states of the Gall 2023 epithelial model (metabolites and
  # epithelial cell compartments) with no canonical analogue in
  # inst/references/compartment-names.md.
  paper_specific_compartments <- c(
    "fump",
    "fudp",
    "futp",
    "stem",
    "prolif",
    "nonprolif",
    "villus",
    "carr",
    "capo",
    "varr",
    "vapo",
    "citr"
  )

  units <- list(
    time = "h",
    dosing = "mg/kg (intravenous; converted to a plasma-concentration increment through f(central))",
    concentration = "ng/mL for 5-FU plasma, metabolites and citrulline; cell counts (two-dimensional crypt-villus strip) for epithelial states"
  )

  compartmentData <- list(
    central = list(analyte = "5-fluorouracil", units = "ng/mL", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "5-fluorouracil", units = "ng/mL", specimen = "tissue", verified = TRUE),
    fump = list(analyte = "fluorouridine monophosphate (FUMP)", units = "ng/mL", specimen = "tissue", verified = TRUE),
    fudp = list(analyte = "fluorouridine diphosphate (FUDP)", units = "ng/mL", specimen = "tissue", verified = TRUE),
    futp = list(analyte = "fluorouridine triphosphate (FUTP)", units = "ng/mL", specimen = "tissue", verified = TRUE),
    stem = list(analyte = "stem cells", units = "cells", specimen = "tissue", verified = TRUE),
    prolif = list(
      analyte = "proliferative crypt progenitor cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    nonprolif = list(
      analyte = "non-proliferative (differentiated) crypt cells",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    villus = list(analyte = "healthy villus epithelial cells", units = "cells", specimen = "tissue", verified = TRUE),
    carr = list(analyte = "drug-arrested crypt cells", units = "cells", specimen = "tissue", verified = TRUE),
    capo = list(analyte = "apoptotic crypt cells", units = "cells", specimen = "tissue", verified = TRUE),
    varr = list(analyte = "drug-arrested villus cells", units = "cells", specimen = "tissue", verified = TRUE),
    vapo = list(analyte = "apoptotic villus cells", units = "cells", specimen = "tissue", verified = TRUE),
    citr = list(analyte = "citrulline", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "mouse (C57BL/6J; small-intestinal epithelium)",
    n_subjects = "Literature group-mean data; no individual animals were modelled (deterministic typical-value fit).",
    n_studies = 1L,
    sex_female_pct = 0,
    disease_state = "Healthy small-intestinal epithelium challenged with 5-FU",
    dose_range = "20 or 50 mg/kg intravenously every 12 h for 4 days, then a 2-day recovery period (Jardi et al. 2022, the source of the mouse cell-count, citrulline, diarrhea and transcriptomic data)",
    regions = "Janssen Pharmaceutica, Belgium (preclinical)",
    notes = paste(
      "The epithelial, citrulline, diarrhea and gene-transcription parameters",
      "were fitted to C57BL/6J mouse jejunum data from Jardi et al. (2022,",
      "Arch Toxicol), in which mice received 20 or 50 mg/kg 5-FU twice daily",
      "for 4 days. 5-FU plasma PK and the metabolite rate constants were",
      "fitted to Simcyp mouse predictions and to measured mouse metabolite",
      "concentrations. Ki-67+ cells were taken to represent all proliferative",
      "crypt cells (stem plus proliferative progenitors). No inter-individual",
      "variability and no residual error: the model is deterministic, used to",
      "reproduce group-mean trajectories."
    )
  )

  ini({
    # ==================================================================
    # 5-FU plasma PK and metabolism -- Gall 2023 Table 2, '5-FU PK model'
    # (posterior means). The plasma and metabolite states are concentrations
    # (ng/mL); a mg/kg intravenous dose is scaled to its ng/mL plasma
    # increment through f(central) = kadj. ke_fuplasma is the combined plasma
    # loss 'Ke FUplasma + K FUMP'; the non-metabolic elimination used in
    # model() is kel = ke_fuplasma - kfump.
    # ==================================================================
    kadj <- 1123
    label("Dose-to-concentration adjustment constant Ka (ng x kg/mg/mL)") # Table 2, Ka = 1123 +/- 122.9
    ke_fuplasma <- 5.625
    label("Combined 5-FU plasma loss rate Ke FUplasma + K FUMP (1/h)") # Table 2, KeFUplasma + KFUMP = 5.625 +/- 0.3919
    kper <- 2.160
    label("5-FU plasma-to-peripheral transfer rate Kper (1/h)") # Table 2, Kper = 2.160 +/- 0.1707
    kfromper <- 1.516
    label("5-FU peripheral-to-plasma transfer rate Kfromper (1/h)") # Table 2, Kfromper = 1.516 +/- 0.0385
    kfump <- 1.555
    label("5-FU conversion rate into FUMP KFUMP (1/h)") # Table 2, KFUMP = 1.555 +/- 0.2071
    kfudp <- 0.2973
    label("FUMP conversion rate into FUDP KFUDP (1/h)") # Table 2, KFUDP = 0.2973 +/- 0.0169
    kfutp <- 0.8478
    label("FUDP conversion rate into FUTP KFUTP (1/h)") # Table 2, KFUTP = 0.8478 +/- 0.1069
    kefutp <- 1.010
    label("FUTP clearance rate KeFUTP (1/h)") # Table 2, KeFUTP = 1.010 +/- 0.1548
    kfdump <- 0.1452
    label("Rapid-equilibrium constant for FdUMP KFdUMP (unitless)") # Table 2, KFdUMP = 0.1452 +/- 0.01580
    kfdutp <- 3.154
    label("Rapid-equilibrium constant for FdUTP KFdUTP (unitless)") # Table 2, KFdUTP = 3.154 +/- 0.6003

    # ==================================================================
    # Epithelial cell dynamics -- Gall 2023 Table 1, 'Mouse 5-FU' column.
    # Proliferation rates and the apoptotic detachment rate are literature
    # values (fixed); steady-state cell numbers are direct data estimates.
    # ==================================================================
    lambda_s <- fixed(log(2) / 21)
    label("Stem-cell proliferation rate lambda_s (1/h)") # Table 1: ln(2)/21, literature value [54]
    lambda_p <- fixed(log(2) / 12)
    label("Proliferative-cell proliferation rate lambda_P (1/h)") # Table 1: ln(2)/12, literature value [9]
    tau_apo <- fixed(1 / 12)
    label("Apoptotic-cell detachment rate tau_Apo (1/h)") # Table 1: 1/12, literature value [55]
    gam <- 3.41394
    label("Villus-to-crypt proliferation feedback strength gamma (unitless)") # Table 1: 3.41 +/- 0.0793
    mm <- 1.468566
    label("Crypt-to-villus migration shape coefficient m (unitless)") # Table 1: 1.469 +/- 0.0418
    sss <- fixed(2)
    label("Steady-state stem-cell number s_ss (cells)") # Table 1: 2, literature value [8]
    pss <- fixed(18.3)
    label("Steady-state proliferative-progenitor number p_ss (cells)") # Table 1: 18.3, direct data estimate (Ki67+ minus s_ss)
    dss <- fixed(5.3)
    label("Steady-state non-proliferative crypt-cell number d_ss (cells)") # Table 1: 5.3, direct data estimate (crypt minus p_ss minus s_ss)
    vss <- 106
    label("Steady-state villus-cell number v_ss (cells)") # Table 1: 106 +/- 15.1

    # ==================================================================
    # Drug-effect parameters -- Gall 2023 Table 1, 'Mouse 5-FU'. kkill_fdutp
    # is 0 (FdUTP is modelled but assigned no cytotoxicity), so c50_fdutp is
    # inert and fixed at the authors' placeholder of 1.
    # ==================================================================
    kinhib <- 0.793992
    label("5-FU proliferation-inhibition (arrest) strength k_inhib (unitless)") # Table 1: 0.7940 +/- 0.0818
    kkill <- 0.55509
    label("5-FU cell-death (apoptosis) strength k_kill (unitless)") # Table 1: 0.5551 +/- 0.0660
    kkill_fdutp <- fixed(0)
    label("FdUTP cell-death strength (unitless)") # authors' code: k_kill_FdUTP = 0 (FdUTP assigned no cytotoxicity)
    c50 <- 1317.88
    label("FUTP concentration for half-maximal effect C50 (ng/mL)") # Table 1: 1318 +/- 49.48
    c50_fdutp <- fixed(1)
    label("FdUTP concentration for half-maximal effect (ng/mL)") # authors' code placeholder; inert because kkill_fdutp = 0
    hill <- fixed(5)
    label("Hill coefficient of the FUTP effect (unitless)") # Table 1: 5 (fixed; fitted value >> 5)

    # ==================================================================
    # Plasma citrulline -- Gall 2023 Table 1. kcitru is fitted; the clearance
    # Cl_Citru is derived at steady state in model() as kcitru v_ss^2/citru_ss.
    # ==================================================================
    citss <- 8237.6362
    label("Steady-state plasma citrulline concentration citru_ss (ng/mL)") # Table 1: 8237.6 +/- 1321
    kcitru <- 0.02325
    label("Citrulline production rate k_Citru (ng/mL/h/cells^2)") # Table 1: 0.0233 +/- 0.00222

    # ==================================================================
    # Probability of diarrhea -- Gall 2023 Table 1 logistic regression on
    # fractional villus-cell loss. beta_0 is fixed to give ~0 placebo
    # probability.
    # ==================================================================
    b0 <- fixed(3.5)
    label("Diarrhea logistic intercept beta_0 (unitless)") # Table 1: 3.5 (fixed to achieve placebo probability ~0)
    b1 <- -5.1298
    label("Diarrhea logistic slope beta_1 (unitless)") # Table 1: -5.1298 +/- 0.8558
    loss50 <- 0.069511
    label("Villus-cell-loss value for half-maximal diarrhea odds Loss50 (unitless)") # Table 1: 0.06951 +/- 0.0421

    # ==================================================================
    # p53-pathway gene transcription -- Gall 2023 Table 1. Each z-score is an
    # instantaneous linear function of the arrested + apoptotic cell fraction.
    # ==================================================================
    zss_p53 <- 9.842
    label("Steady-state z-score of p53 (unitless)") # Table 1: 9.842 +/- 0.06588
    k_p53 <- 4.5405
    label("Slope of the p53 z-score change (unitless)") # authors' code k_TP53 = 4.5405 (Table 1 prints 5.540 +/- 1.041; see vignette Errata)
    zss_p21 <- 12.3035
    label("Steady-state z-score of p21 (unitless)") # Table 1: 12.30 +/- 0.06203
    k_p21 <- 9.7568
    label("Slope of the p21 z-score change (unitless)") # Table 1: 9.757 +/- 0.9799
    zss_mdm2 <- 8.207
    label("Steady-state z-score of mdm2 (unitless)") # Table 1: 8.207 +/- 0.07607
    k_mdm2 <- 9.7394
    label("Slope of the mdm2 z-score change (unitless)") # Table 1: 9.739 +/- 1.201
    zss_bax <- 7.8517
    label("Steady-state z-score of bax (unitless)") # Table 1: 7.852 +/- 0.09150
    k_bax <- 12.4458
    label("Slope of the bax z-score change (unitless)") # Table 1: 12.45 +/- 1.445
    zss_pidd1 <- 6.71
    label("Steady-state z-score of pidd1 (unitless)") # Table 1: 6.710 +/- 0.09470
    k_pidd1 <- 11.8412
    label("Slope of the pidd1 z-score change (unitless)") # Table 1: 11.84 +/- 1.496
    zss_sesn2 <- 7.3233
    label("Steady-state z-score of sesn2 (unitless)") # Table 1: 7.323 +/- 0.06949
    k_sesn2 <- 10.667
    label("Slope of the sesn2 z-score change (unitless)") # Table 1: 10.67 +/- 1.098
  })
  model({
    # --- Steady-state-derived epithelial transfer and shedding rates
    # (Gall 2023 Table 1 'Steady-state estimation' rows). ----------------
    delta_sp <- lambda_s
    delta_pd <- lambda_p * (1 + gam / 2) + delta_sp * sss / pss
    delta_dv <- delta_pd * pss / dss
    tau_v <- delta_dv * dss / vss
    clcitru <- kcitru * vss^2 / citss

    # --- 5-FU plasma PK and intracellular metabolism (Eqs S1-S7). The
    # non-metabolic plasma elimination is the combined rate minus KFUMP. ---
    kel <- ke_fuplasma - kfump
    d/dt(central) <- -(kel + kper + kfump) * central + kfromper * peripheral1
    d/dt(peripheral1) <- kper * central - kfromper * peripheral1
    d/dt(fump) <- kfump * central - kfudp * fump
    d/dt(fudp) <- kfudp * fump - kfutp * fudp
    d/dt(futp) <- kfutp * fudp - kefutp * futp
    f(central) <- kadj

    # FdUMP and FdUTP in rapid equilibrium with FUDP (Eqs S6-S7); both are
    # inert outputs here (kkill_fdutp = 0).
    fdump <- kfdump * fudp
    fdutp <- kfdutp * fudp

    # --- Drug effects on proliferative cells (Eqs S19-S20). --------------
    emax_futp <- futp^hill / (futp^hill + c50^hill)
    emax_fdutp <- (kfdutp * futp)^hill / ((kfdutp * futp)^hill + c50_fdutp^hill)
    arrest <- kinhib * emax_futp
    kill <- kkill_fdutp * emax_fdutp + kkill * emax_futp

    # --- Epithelial cell dynamics (Eqs S10-S18). ------------------------
    feedback <- gam * vss / (vss + villus)
    migration <- ((stem + prolif) / (sss + pss))^mm
    d/dt(stem) <- lambda_s * stem * (1 - arrest - kill) - delta_sp * stem * stem / sss
    d/dt(prolif) <- delta_sp * stem * stem / sss +
      lambda_p * prolif * (1 + feedback - arrest - kill) - delta_pd * prolif
    d/dt(nonprolif) <- delta_pd * prolif - delta_dv * nonprolif * migration
    d/dt(villus) <- delta_dv * nonprolif * migration - tau_v * villus
    d/dt(carr) <- arrest * (lambda_s * stem + lambda_p * prolif) - delta_dv * carr * migration
    d/dt(capo) <- kill * (lambda_s * stem + lambda_p * prolif) -
      delta_dv * capo * migration - tau_apo * capo
    d/dt(varr) <- delta_dv * carr * migration - tau_v * varr
    d/dt(vapo) <- delta_dv * capo * migration - tau_v * vapo - tau_apo * vapo

    # --- Plasma citrulline (Eq S21). ------------------------------------
    d/dt(citr) <- kcitru * villus^2 - clcitru * citr

    stem(0) <- sss
    prolif(0) <- pss
    nonprolif(0) <- dss
    villus(0) <- vss
    citr(0) <- citss

    # --- Clinical-readout observables. ----------------------------------
    Cc <- central
    vtot <- villus + varr + vapo
    cellloss <- max(0, 1 - vtot / vss)
    # Probability of diarrhea (Eqs S22-S23), a logistic function of villus loss.
    prob_diarrhea <- 1 - expit(b0 + b1 * cellloss / (cellloss + loss50))
    # p53-pathway gene z-scores (Eq S24): instantaneous linear functions of
    # the arrested + apoptotic cell fraction.
    apo_arr_frac <- (carr + capo + varr + vapo) /
      (stem + prolif + nonprolif + villus + carr + capo + varr + vapo)
    z_p53 <- zss_p53 + k_p53 * apo_arr_frac
    z_p21 <- zss_p21 + k_p21 * apo_arr_frac
    z_mdm2 <- zss_mdm2 + k_mdm2 * apo_arr_frac
    z_bax <- zss_bax + k_bax * apo_arr_frac
    z_pidd1 <- zss_pidd1 + k_pidd1 * apo_arr_frac
    z_sesn2 <- zss_sesn2 + k_sesn2 * apo_arr_frac
  })
}
