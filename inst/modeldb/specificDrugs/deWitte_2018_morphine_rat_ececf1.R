deWitte_2018_morphine_rat_ececf1 <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). de Witte 2018 model ECECF1 (Table",
    "7): brain-ECF morphine equilibrates with an effect compartment",
    "(keo, Figure 4b) whose concentration drives the EEG delta-band",
    "(0.5-4.5 Hz) amplitude through a sigmoidal Emax model. GF120918",
    "(elacridar) co-infusion raises Emax.",
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
        "column); in this model it acts on the maximal EEG effect Emax (Table 7, separately estimated Emax-Pgp)."
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
    # Effect-compartment morphine as an ECF-equivalent CONCENTRATION.
    effect = list(
      analyte = "morphine",
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
      "and 1.75 nM (ECF). The ECF-driven models (Table 7) were fitted to the reduced dataset of animals with brain ECF data; the reduced-dataset animal count is not printed."
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
    # EEG PD -- Table 7, column 'ECECF1 effect compartment' (OFV 26255.1,
    # AIC 26273.1). Figure 4b.
    # ------------------------------------------------------------------
    lke0 <- log(0.161); label("Effect-compartment equilibration rate constant, keo (1/min)") # Table 7: keo = 0.161 /min (CV 40%)
    le0 <- log(49.1); label("Baseline EEG delta amplitude at t = 0, E0 (uV)") # Table 7: E0 = 49.1 uV (CV 6%)
    lemax <- log(24.9); label("Maximal morphine-induced increase in EEG amplitude, Emax (uV)") # Table 7: Emax = 24.9 uV (CV 19%)
    e_conmed_elacridar_emax <- log(43.2 / 24.9); label("Log-ratio of Emax with GF120918 to Emax without (unitless)") # Table 7: Emax-Pgp = 43.2 uV (CV 15%)
    lec50 <- log(182); label("Effect-compartment concentration giving half-maximal effect, EC50 (nM)") # Table 7: EC50 = 182 nM (CV 26%)
    lhill <- log(2.02); label("Hill coefficient, NH (unitless)") # Table 7: NH = 2.02 (CV 43%)
    drift <- -0.0373; label("Linear change of the EEG amplitude per minute of experiment time, slope (uV/min)") # Table 7: slope = -0.0373 uV/min (CV 19%)

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
    etale0 ~ 0.0696 # Table 7: omega^2 E0 = 0.0696 (CV 26%)

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
    propSd_eeg <- sqrt(0.0593); label("Proportional residual SD, EEG delta amplitude (fraction)") # Table 7: sigma^2 prop = 0.0593 (CV 10%)

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

    # ---- Effect compartment fed by brain ECF (Figure 4b, keo in and out) ----
    ke0 <- exp(lke0)
    d/dt(effect) <- ke0 * (brain_ecf - effect)

    # ---- EEG amplitude, sigmoidal Emax (Figure 4b) ----
    e0 <- exp(le0 + etale0)
    emax <- exp(lemax + e_conmed_elacridar_emax * CONMED_ELACRIDAR)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    cdrv <- max(effect, 0)
    eeg <- e0 + drift * t + emax * cdrv^hill / (ec50^hill + cdrv^hill)

    Cc ~ prop(propSd)
    Cbrain_ecf ~ prop(propSd_Cbrain_ecf)
    eeg ~ prop(propSd_eeg)

  })
}
