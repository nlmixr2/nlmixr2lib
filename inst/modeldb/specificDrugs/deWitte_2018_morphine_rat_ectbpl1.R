deWitte_2018_morphine_rat_ectbpl1 <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). de Witte 2018 model EC-TBPL1, the",
    "selected combined effect-compartment / target-binding model on",
    "plasma concentrations (Table 6): plasma morphine equilibrates with",
    "an effect compartment (keo), where it binds its target with kon =",
    "koff/KD and koff (Supplement S1 eqs 13-15), and the EEG delta-band",
    "(0.5-4.5 Hz) amplitude rises linearly with target occupancy. Total",
    "target Rtot is fixed at 1 nM. GF120918 (elacridar) co-infusion lowers",
    "the apparent KD.",
    "Plasma morphine follows a three-compartment model (Supplement S1 Table",
    "S1) and drives a brain extracellular-fluid (ECF) compartment with",
    "passive diffusion, saturable active influx and first-order efflux",
    "(Supplement S1 eq 9, Table S2); the plasma and ECF layers were fitted",
    "sequentially and then held at their individual estimates for the EEG",
    "fit. The EEG baseline changes linearly with experiment time, E0 + slope",
    "* t, where t = 0 is the start of EEG recording (the 10-min morphine",
    "infusion starts about 60 min later). The saturable-influx C50 is not",
    "reported in de Witte 2018 and is fixed to the Groenendaal 2007 estimate",
    "for the same dataset (9.92 ng/mL = 34.8 nM).",
    sep = " "
  )
  reference <- paste(
    "de Witte WEA, Rottschafer V, Danhof M, van der Graaf PH, Peletier LA,",
    "de Lange ECM. Modelling the delay between pharmacokinetics and EEG",
    "effects of morphine in rats: binding kinetic versus effect compartment",
    "models. J Pharmacokinet Pharmacodyn. 2018;45(4):621-635.",
    "doi:10.1007/s10928-018-9593-x. The saturable brain-influx C50, which",
    "de Witte 2018 does not report, is fixed to the estimate for the same",
    "dataset in Groenendaal D et al. Br J Pharmacol. 2007;151(5):701-712,",
    "doi:10.1038/sj.bjp.0707257 (Table 3).",
    sep = " "
  )
  vignette <- "deWitte_2018_morphine_rat_eeg"

  # `lkdiff_bbb`, `lkefflux_bbb`, `lNstarMax` and `lc50` follow
  # Groenendaal_2007_morphine_brain_rat.R, the earlier model of the same
  # rat dataset and the same brain-ECF structure.
  paper_specific_etas <- c("etalkdiff_bbb", "etalNstarMax")
  paper_specific_residual_sds <- c("propSd_Cbrain_ecf", "propSd_eeg")

  units <- list(
    time = "min",
    dosing = "nmol",
    concentration = "nM"
  )

  covariateData <- list(
    CONMED_ELACRIDAR = list(
      description = "Indicator for co-infusion of the P-glycoprotein inhibitor GF120918 (elacridar): 1 = GF120918 arm, 0 = vehicle arm.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (vehicle co-infusion, no elacridar)",
      notes = paste(
        "de Witte 2018 Table 1: 9 of the 29 rats given 4 mg/kg morphine",
        "received GF120918 as a continuous infusion; the other dose levels",
        "were vehicle only. Regimen (Groenendaal 2007, the source",
        "experiment): 6 mg/kg 1-min bolus then a continuous infusion, started",
        "120 min before morphine and maintained throughout, so the indicator",
        "is time-fixed per animal. GF120918 had no retained effect on plasma",
        "PK (Table S1) or on brain ECF distribution (Table S2 selected",
        "column); in this model it acts on the apparent equilibrium dissociation constant KD (Table 6, separately estimated KD-Pgp)."
      ),
      source_name = "Pgp inhibitor (GF120918, Table 1 '+' / '-')"
    )
  )

  compartmentData <- list(
    central = list(analyte = "morphine", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "morphine", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "morphine", units = "nmol", specimen = "plasma", verified = TRUE),
    # Supplement S1 eq 9 is written for the ECF CONCENTRATION
    # d(A_ECF/V_ECF)/dt and V_ECF is never identified, so the state holds
    # nM rather than an amount (same as Groenendaal_2007_morphine_brain_rat.R).
    brain_ecf = list(analyte = "morphine", units = "nM (concentration state)", specimen = "brain ISF", verified = TRUE),
    # Effect-compartment morphine as a plasma-equivalent CONCENTRATION.
    effect = list(
      analyte = "morphine",
      units = "nM (concentration state)",
      specimen = "not applicable",
      verified = TRUE
    ),
    # Receptor-bound morphine (drug-target complex RL) as a CONCENTRATION;
    # S1 eqs 13-17 are written for A_RL/V_RL. The bound drug is not
    # removed from the driving compartment.
    complex = list(
      analyte = "morphine-opioid receptor complex",
      units = "nM (concentration state)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "rat (male Wistar)",
    n_subjects = 59L,
    n_studies = 1L,
    age_range = "adult (not reported)",
    weight_range = "not reported in de Witte 2018; 250-350 g in the source experiment (Groenendaal 2007)",
    sex_female_pct = 0,
    race_ethnicity = NA,
    disease_state = "healthy, chronically instrumented rats (EEG electrodes, microdialysis probe in a subset)",
    dose_range = "0 (vehicle), 4, 10 or 40 mg/kg morphine as a 10-min intravenous infusion; 4 mg/kg with and without GF120918",
    regions = "Leiden, The Netherlands (single laboratory)",
    notes = paste(
      "Reanalysis of the Groenendaal 2007 experiments (de Witte 2018 Table",
      "1): 59 dosed animals (0 mg/kg N = 5; 4 mg/kg N = 20 vehicle + 9",
      "GF120918; 10 mg/kg N = 11; 40 mg/kg N = 14). Plasma PK in 54",
      "animals, brain ECF (microdialysis) in 29, EEG in 54. The EEG",
      "endpoint is the delta-band (0.5-4.5 Hz) amplitude in uV, recorded",
      "every minute and averaged over 3-min intervals. LLOQ 88 nM (plasma)",
      "and 1.75 nM (ECF). EC-TBPL1 was fitted to all EEG animals."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Plasma PK -- Supplement S1 Table S1, column '3-cmp model 4 IIV
    # parameters' (OFV 9205), the selected plasma model. Whole-animal
    # volumes (L) and clearances (L/min); doses in nmol give nM. S1 eqs
    # 1-8 name the inter-compartmental clearances Q1 (to 'Tissue', V2)
    # and Q2 (to 'Tissue2', V3).
    # ------------------------------------------------------------------
    lcl <- log(0.028); label("Plasma clearance CL (L/min)") # Table S1: CL = 0.028 L/min (CV 22%)
    lvc <- log(0.17); label("Central (plasma) volume Vc (L)") # Table S1: Vc = 0.17 L (CV 49%)
    lq <- log(0.019); label("Inter-compartmental clearance to peripheral1, Q1 (L/min)") # Table S1: Q1 = 0.019 L/min (CV 66%)
    lvp <- log(1.3); label("First peripheral volume V2 (L)") # Table S1: V2 = 1.3 L (CV 48%)
    lq2 <- log(0.031); label("Inter-compartmental clearance to peripheral2, Q2 (L/min)") # Table S1: Q2 = 0.031 L/min (CV 30%)
    lvp2 <- log(0.36); label("Second peripheral volume V3 (L)") # Table S1: V3 = 0.36 L (CV 32%)

    # ------------------------------------------------------------------
    # Brain ECF -- Supplement S1 eq 9 and Table S2, column 'IIV on kdiff,
    # Nmax' (OFV -1126), the selected ECF model (S1 text: the GF120918
    # effect did not lower the OFV enough to be retained).
    # ------------------------------------------------------------------
    lkdiff_bbb <- log(0.0025); label("Passive diffusion rate constant between plasma and brain ECF, kdiff (1/min)") # Table S2: kdiff = 0.0025 /min (CV 17%)
    lkefflux_bbb <- log(0.020); label("First-order efflux rate constant from brain ECF, keff (1/min)") # Table S2: keff = 0.020 /min (CV 11%)
    lNstarMax <- log(2.6); label("Maximal saturable active influx rate into brain ECF, Nmax (nM/min)") # Table S2: Nmax = 2.6 nM/min (CV 21%)
    # NOT REPORTED in de Witte 2018 (no table, text, figure or supplement
    # value). Fixed to Groenendaal 2007 Table 3, C50 = 9.92 ng/mL (same rat
    # dataset, same ECF structure, which S1 calls 'the original one
    # compartment ECF model'), converted with the morphine free-base molar
    # mass 285.34 g/mol: 9.92 * 1000 / 285.34 = 34.77 nM.
    lc50 <- fixed(log(9.92 * 1000 / 285.34)); label("Plasma concentration giving half-maximal saturable influx, C50 (nM)") # Groenendaal 2007 Table 3 (non-paper provenance; see vignette Errata)

    # ------------------------------------------------------------------
    # EEG PD -- Table 6, column 'EC-TBPL1 selected model' (OFV 44790.9,
    # AIC 44808.9). Supplement S1 eqs 13-15.
    # ------------------------------------------------------------------
    lke0 <- log(0.0327); label("Effect-compartment equilibration rate constant, keo (1/min)") # Table 6: keo = 0.0327 /min (CV 17%)
    lkoff <- log(0.0275); label("Drug-target dissociation rate constant, koff (1/min)") # Table 6: koff = 0.0275 /min (CV 14%); the Discussion's '0.275' is a typo
    lkd <- log(1520); label("Equilibrium dissociation constant, KD (nM)") # Table 6: KD = 1520 nM (CV 34%)
    # Table 6 estimates a separate KD-Pgp for the GF120918 arm; koff and keo are shared.
    e_conmed_elacridar_kd <- log(296 / 1520); label("Log-ratio of KD with GF120918 to KD without (unitless)") # Table 6: KD-Pgp = 296 nM (CV 93%)
    lrtot <- fixed(log(1)); label("Total target concentration, Rtot (nM)") # Table 6: Rtot = 1 FIX nM
    le0 <- log(45.0); label("Baseline EEG delta amplitude at t = 0, E0 (uV)") # Table 6: E0 = 45.0 uV (CV 4%)
    lemax <- log(31.8); label("EEG amplitude increase at full target occupancy, Emax (uV)") # Table 6: Emax = 31.8 uV (CV 11%)
    drift <- -0.0276; label("Linear change of the EEG amplitude per minute of experiment time, slope (uV/min)") # Table 6: slope = -0.0276 uV/min (CV 15%)

    # ------------------------------------------------------------------
    # Inter-animal variability, log-normal (paper eqs 1-2); the table rows
    # are omega^2 variances.
    # ------------------------------------------------------------------
    etalcl ~ 0.33 # Table S1: omega^2 CL = 0.33 (CV 25%)
    etalvc ~ 0.63 # Table S1: omega^2 Vc = 0.63 (CV 35%)
    etalq ~ 0.62 # Table S1: omega^2 Q1 = 0.62 (CV 33%)
    etalvp ~ 0.66 # Table S1: omega^2 V2 = 0.66 (CV 28%)
    etalkdiff_bbb ~ 0.36 # Table S2: omega^2 kdiff = 0.36 (CV 39%)
    etalNstarMax ~ 0.42 # Table S2: omega^2 Nmax = 0.42 (CV 55%); omega^2 keff = 0 FIX, so keff has no eta
    etale0 ~ 0.111 # Table 6: omega^2 E0 = 0.111 (CV 20%)

    # ------------------------------------------------------------------
    # Residual error, proportional (paper eq 3); tables give sigma^2.
    # Table S1 prints the plasma sigma^2 as '0.79 (14)'. It is read as
    # 0.079: the neighbouring columns print 0.0766 and 0.071, and the
    # Figure S2 observed-versus-IPRED scatter is far tighter than a 0.89
    # proportional SD (see vignette Errata). The additive plasma term is
    # 0 FIX in this column.
    # ------------------------------------------------------------------
    propSd <- sqrt(0.079); label("Proportional residual SD, plasma morphine (fraction)") # Table S1: sigma^2 prop printed '0.79 (14)', read as 0.079
    propSd_Cbrain_ecf <- sqrt(0.11); label("Proportional residual SD, brain ECF morphine (fraction)") # Table S2: sigma^2 prop = 0.11 (CV 18%)
    propSd_eeg <- sqrt(0.057); label("Proportional residual SD, EEG delta amplitude (fraction)") # Table 6: sigma^2 prop = 0.057 (CV 7%)

  })

  model({
    # ---- Plasma PK (Supplement S1 eqs 1-8) ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc

    # ---- Brain ECF (Supplement S1 eq 9) ----
    # The paper assumes ECF exchange does not change plasma concentrations,
    # so nothing feeds back into d/dt(central).
    kdiff_bbb <- exp(lkdiff_bbb + etalkdiff_bbb)
    kefflux_bbb <- exp(lkefflux_bbb)
    NstarMax <- exp(lNstarMax + etalNstarMax)
    c50 <- exp(lc50)

    d/dt(brain_ecf) <- kdiff_bbb * (Cc - brain_ecf) +
      NstarMax * Cc / (c50 + Cc) -
      kefflux_bbb * brain_ecf

    Cbrain_ecf <- brain_ecf

    # ---- Effect compartment (Supplement S1 eq 13) ----
    # S1 eq 13 is printed as keo * (Ceff - Cp), which diverges from zero
    # initial conditions; the sign is reversed here to the standard
    # keo * (Cp - Ceff) of the Figure 3 caption (see vignette Errata).
    ke0 <- exp(lke0)
    d/dt(effect) <- ke0 * (Cc - effect)

    # ---- Target binding (Supplement S1 eq 14) ----
    #   dRL/dt = kon * C * (Rtot - RL) - koff * RL,  kon = koff / KD
    # Rtot is fixed at 1 nM, far below the drug concentrations, so drug
    # depletion by binding is negligible and occupancy RL / Rtot does not
    # depend on the Rtot value.
    koff <- exp(lkoff)
    kd <- exp(lkd + e_conmed_elacridar_kd * CONMED_ELACRIDAR)
    kon <- koff / kd
    rtot <- exp(lrtot)

    d/dt(complex) <- kon * effect * (rtot - complex) - koff * complex

    # ---- EEG amplitude, linear in occupancy (Supplement S1 eq 15) ----
    e0 <- exp(le0 + etale0)
    emax <- exp(lemax)
    eeg <- e0 + drift * t + emax * complex / rtot

    Cc ~ prop(propSd)
    Cbrain_ecf ~ prop(propSd_Cbrain_ecf)
    eeg ~ prop(propSd_eeg)

  })
}
