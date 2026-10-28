Thorsted_2023_neutrophil_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). Semi-mechanistic model of neutrophil-mediated",
    "phagocytosis and digestion of Acinetobacter baumannii in a murine",
    "pneumonia model, with cyclophosphamide-induced neutropenia (Thorsted",
    "2023), in vivo arm. Intraperitoneal cyclophosphamide (depot -> central,",
    "one compartment; rat parameters allometrically scaled to a 25 g mouse",
    "and fixed) kills neutrophil progenitors through a linear slope in a",
    "Friberg semi-physiological myelosuppression chain (prol, transit1-3,",
    "circ) with (circ0 / circ)^gamma feedback. Circulating neutrophils",
    "migrate into the lung (neut_lung) at the maturation rate ktr, and",
    "bacterial deposition triggers a surge in that influx with amplitude",
    "NAMP, centre NT0 and width NSW. Bacteria are dosed in CFU/g lung into the",
    "growing state bact_susceptible (no lag in vivo), move to the resting",
    "state bact_resting at kSR = (S + R) * (kgrowth - kdeath) / Bmax, and are",
    "phagocytosed by lung neutrophils into bact_phagocytosed and digested",
    "from there, with the in vitro Emax-on-neutrophils and phagocytosed-to-",
    "neutrophil-ratio structure, a lower in vivo kN,max and no time-related",
    "loss. An additional in vivo phagocytic capacity (kM,kill, represented",
    "by a macrophage) moves bacteria into bact_phagocytosed_mac and digests",
    "them at the same rate; it scales with circ / circ0 and decays with time",
    "since inoculation (kM,loss). Outputs: log10 CFU/g lung (log_cfu), log10",
    "circulating neutrophils per mL (log_anc) and log10 lung neutrophils per",
    "lung (log_neut_lung). The in vitro arm of the same simultaneous fit is",
    "Thorsted_2023_neutrophil_invitro."
  )
  reference <- paste(
    "Thorsted A, Pham AD, Friberg LE, Nielsen EI. Model-based assessment of",
    "neutrophil-mediated phagocytosis and digestion of bacteria across in",
    "vitro and in vivo studies. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(12):1972-1987. doi:10.1002/psp4.13046. PMID: 37700716.",
    "PMCID: PMC10725272. Parameter estimates: Table 2, 'Final model' column.",
    "Equations: Materials and Methods, Equations 1-14. Model code: Data S2",
    "(NONMEM control stream); digitized data: Table S1. Myelosuppression",
    "structure and the rat cyclophosphamide parameters scaled here: Friberg",
    "LE, Sandstrom M, Karlsson MO. Invest New Drugs. 2010;28:744-753.",
    sep = " "
  )
  vignette <- "Thorsted_2023_neutrophil_bacterial_killing"
  units <- list(
    time = "h",
    dosing = "mg (cyclophosphamide, intraperitoneal); CFU/g lung (bacterial inoculum)",
    concentration = "log10 CFU/g lung, log10 cells/mL and log10 cells/lung (observations); mg/mL (cyclophosphamide)"
  )

  covariateData <- list()

  # Bacterial life-cycle states split one phenotype along a growth axis, not
  # along the drug-resistance axis of the registered ^bact_(susceptible|...)
  # regex (see Saporta_2026_meropenem). neut_lung is the lung neutrophil pool
  # that the paper adds downstream of the Friberg circulating compartment.
  paper_specific_compartments <- c(
    "bact_susceptible",
    "bact_resting",
    "bact_phagocytosed",
    "bact_phagocytosed_mac",
    "neut_lung",
    "log_anc",
    "log_neut_lung"
  )

  compartmentData <- list(
    depot = list(analyte = "cyclophosphamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cyclophosphamide", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(
      analyte = "neutrophil progenitor cells",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "maturing neutrophils, transit 1",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "maturing neutrophils, transit 2",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "maturing neutrophils, transit 3",
      units = "cells/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(analyte = "circulating neutrophils", units = "cells/mL", specimen = "whole blood", verified = TRUE),
    neut_lung = list(analyte = "lung neutrophils", units = "cells/g lung", specimen = "tissue", verified = TRUE),
    bact_susceptible = list(
      analyte = "Acinetobacter baumannii, growing (S) state",
      units = "CFU/g lung",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_resting = list(
      analyte = "Acinetobacter baumannii, resting (R) state",
      units = "CFU/g lung",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_phagocytosed = list(
      analyte = "Acinetobacter baumannii phagocytosed by neutrophils and not yet digested (PN)",
      units = "CFU/g lung",
      specimen = "tissue",
      verified = TRUE
    ),
    bact_phagocytosed_mac = list(
      analyte = "Acinetobacter baumannii taken up by the additional in vivo phagocytic capacity and not yet digested (PM)",
      units = "CFU/g lung",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (C57BL/6 and Swiss-Webster)",
    n_subjects = 21L,
    n_studies = 5L,
    age_range = "6-12 weeks (where reported)",
    weight_median = "25 g assumed for allometric scaling",
    disease_state = "Acinetobacter baumannii pneumonia (four studies) and uninfected steady-state neutrophil tissue data (one study); one pneumonia study pretreated with cyclophosphamide",
    dose_range = paste(
      "Cyclophosphamide intraperitoneally 96 and 24 h before inoculation:",
      "50/50, 100/50 or 200/150 mg/kg, or vehicle (Mice-2). Bacterial",
      "inocula 5.71 x 10^6 to 5.71 x 10^8 CFU/g lung."
    ),
    notes = paste(
      "Five published mouse studies were digitized (Table 1, Mice-1 to Mice-5): 21 in vivo experiments with 56 CFU, 5 blood-neutrophil and 16 lung-neutrophil observations (Results). n_subjects counts experiments, the unit that NONMEM treated as an individual.",
      "Lung observations were converted to per-gram units with an approximate C57BL/6 lung weight of 0.175 g (Materials and Methods).",
      "No between-experiment variability was estimated ($OMEGA 0 FIX in Data S2); the model carries typical values plus residual error only."
    )
  )

  ini({
    # ---- Bacterial dynamics, in vivo (Table 2, Final model) ----------------
    lkgrowth <- log(0.502); label("Bacterial first-order growth rate constant, in vivo (1/h)") # Table 2: kgrowth (in vivo) = 0.502 1/h, RSE 18%
    lkdeath <- fixed(log(0.179)); label("Bacterial first-order natural death rate constant (1/h)") # Table 2: kdeath = 0.179 1/h, Fixed (shared with in vitro)
    lbmax <- fixed(log(2.5e9)); label("Maximum system bacterial capacity (CFU/g lung)") # Table 2: Bmax = 2.5 x 10^9, Fixed (Data S2: BMAX shared by in vitro and in vivo)

    # ---- Neutrophil-mediated killing, in vivo (Table 2, Final model) -------
    lkn_max <- log(0.820); label("Neutrophil first-order maximum phagocytosis and digestion rate constant, in vivo (1/h)") # Table 2: kN,max (in vivo) = 0.820 1/h, RSE 33%
    kn_loss <- fixed(0); label("Time-related decrease in neutrophil phagocytosis and digestion, in vivo (1/h)") # Table 2: kN,loss (in vivo) = 0, Fixed
    # N50 and P/N50 are shared with the in vitro arm (Data S2 uses the same
    # THETAs for every study).
    ln50_dig <- log(6.55e6); label("Neutrophil concentration giving 50% of the maximum digestion rate (cells/g lung)") # Table 2: N50,dig = 6.55 x 10^6 N/mL, RSE 1.4%
    ln50_phag <- log(1.19e6); label("Neutrophil concentration giving 50% of the maximum phagocytosis rate (cells/g lung)") # Table 2: N50,phag = 1.19 x 10^6 N/mL, RSE 4.3%
    lpn50 <- log(3.41); label("Phagocytosed-bacteria-to-neutrophil ratio halving the phagocytosis and digestion rates (CFU/cell)") # Table 2: P/N50 = 3.41 CFU/N, RSE 22%
    lkm_kill <- log(0.287); label("First-order rate constant of the additional in vivo phagocytic capacity (1/h)") # Table 2: kM,kill = 0.287 1/h, RSE 18%
    lkm_loss <- log(0.0111); label("Time-related decrease in the additional in vivo phagocytic capacity (1/h)") # Table 2: kM,loss = 0.0111 1/h, RSE 45%

    # ---- Cyclophosphamide PK (Table 2: 'Scaled', not estimated) ------------
    # Rat values from Friberg 2010 scaled to a 25 g mouse from a 295 g rat
    # with exponents 0.75 (CL), 1 (V) and -0.25 (ka) (Materials and Methods;
    # Data S2 $PK writes the same expressions).
    lka <- fixed(log(11.3 * (0.025 / 0.295)^(-0.25))); label("Cyclophosphamide intraperitoneal absorption rate constant (1/h)") # Table 2: ka = 20.9 1/h (Scaled); Data S2: KA = 11.3 * (0.025/0.295)**(-0.25)
    lcl <- fixed(log(1.98 * (0.025 / 0.295)^0.75)); label("Cyclophosphamide clearance (mL/h)") # Table 2: CL = 0.311 mL/h (Scaled); Data S2: CL = 1.98 * (0.025/0.295)**0.75
    lvc <- fixed(log(1.8 * (0.025 / 0.295))); label("Cyclophosphamide volume of distribution (mL)") # Table 2: V = 0.153 mL (Scaled); Data S2: VC = 1.8 * (0.025/0.295)**1

    # ---- Neutrophil dynamics (Table 2, Final model) ------------------------
    lmtt <- fixed(log(52.8 * (0.025 / 0.295)^0.136)); label("Mean transit time of neutrophil bone-marrow maturation (h)") # Table 2: MTT = 37.7 h (Scaled); Data S2: MTT = 52.8 * (0.025/0.295)**0.136
    lgamma <- fixed(log(0.149)); label("Feedback exponent on (circ0 / circ) (unitless)") # Table 2: gamma = 0.149, Fixed (rat estimate)
    lslope <- log(1.26); label("Linear cyclophosphamide effect on neutrophil progenitors (mL/mg)") # Table 2: DSLP = 1.26 mL/mg, RSE 11%
    lcirc0 <- log(1.69e6); label("Baseline circulating neutrophil count (cells/mL)") # Table 2: Ncirc,T0 = 1.69 x 10^6 N/mL, RSE 0.7%
    lnamp <- log(44.2); label("Amplitude of the post-inoculation surge in lung neutrophil influx (unitless)") # Table 2: NAMP = 44.2, RSE 24%
    lnt0 <- log(50.5); label("Time after inoculation of the lung neutrophil influx surge peak (h)") # Table 2: NT0 = 50.5 h, RSE 3.0%
    lnsw <- log(28.2); label("Width of the lung neutrophil influx surge (h)") # Table 2: NSW = 28.2 h, RSE 14%
    lung_weight <- fixed(0.175); label("Lung weight converting lung neutrophils per gram to per lung (g)") # Data S2: THETA(15) NLNG = 0.175 FIX; Materials and Methods: approximate C57BL/6 lung weight 0.175 g

    # ---- Residual error (log10 transform-both-sides, additive) -------------
    addSd_log_cfu <- 0.600; label("Additive residual SD on log10 in vivo CFU (log10 CFU/g lung)") # Table 2: ERR (in vivo CFU) = 0.600, RSE 9.4%
    # One SIGMA serves blood and lung neutrophils (Data S2: SIGMA(3,3) for
    # LOCA >= 5). nlmixr2 needs one residual parameter per endpoint, so the
    # single estimate is carried twice with the same value.
    addSd_log_anc <- 0.161; label("Additive residual SD on log10 circulating neutrophils (log10 cells/mL)") # Table 2: ERR (in vivo neutrophils) = 0.161, RSE 10%
    addSd_log_neut_lung <- 0.161; label("Additive residual SD on log10 lung neutrophils (log10 cells/lung; same estimate as addSd_log_anc)") # Table 2: ERR (in vivo neutrophils) = 0.161, single shared SIGMA
  })

  model({
    kgrowth <- exp(lkgrowth)
    kdeath <- exp(lkdeath)
    bmax <- exp(lbmax)
    kn_max <- exp(lkn_max)
    n50_dig <- exp(ln50_dig)
    n50_phag <- exp(ln50_phag)
    pn50 <- exp(lpn50)
    km_kill <- exp(lkm_kill)
    km_loss <- exp(lkm_loss)

    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc

    mtt <- exp(lmtt)
    ktr <- 4 / mtt
    gamma <- exp(lgamma)
    slope <- exp(lslope)
    circ0 <- exp(lcirc0)
    namp <- exp(lnamp)
    nt0 <- exp(lnt0)
    nsw <- exp(lnsw)

    # Time since bacterial inoculation. Data S2 measures the surge and the
    # kM,loss decay from the inoculation time (T - EXTRT, with NT0 + EXTRT);
    # here that is the time since the first bacterial dose into
    # bact_susceptible. Without an inoculation the surge is off (Data S2:
    # NEUSUR = 0 when INOC = 0).
    tinoc <- tafd(bact_susceptible)
    inoculated <- 0
    if (!is.na(tinoc)) inoculated <- 1
    tsince <- 0
    if (inoculated == 1) tsince <- tinoc

    # Equations 7-8: cyclophosphamide PK.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    Cc <- central / vc

    # Equations 9-13: Friberg myelosuppression chain, all initialized at circ0.
    deff <- slope * Cc
    d/dt(prol) <- -ktr * prol + ktr * prol * (1 - deff) * (circ0 / circ)^gamma
    d/dt(transit1) <- -ktr * transit1 + ktr * prol
    d/dt(transit2) <- -ktr * transit2 + ktr * transit1
    d/dt(transit3) <- -ktr * transit3 + ktr * transit2
    d/dt(circ) <- -ktr * circ + ktr * transit3

    # Equation 14: lung neutrophils with the post-inoculation influx surge.
    surge <- inoculated * namp / (((tsince - nt0) / nsw)^6 + 1)
    d/dt(neut_lung) <- -ktr * neut_lung + ktr * (1 + surge) * circ

    prol(0) <- circ0
    transit1(0) <- circ0
    transit2(0) <- circ0
    transit3(0) <- circ0
    circ(0) <- circ0
    neut_lung(0) <- circ0

    # Equations 4-6 with the lung neutrophil concentration as N.
    pn_ratio <- 0
    if (neut_lung > 0) pn_ratio <- bact_phagocytosed / neut_lung
    kn_common <- kn_max * exp(-kn_loss * tsince) * (1 - pn_ratio / (pn_ratio + pn50))
    kn_dig <- kn_common * neut_lung / (neut_lung + n50_dig)
    kn_phag <- kn_common * neut_lung / (neut_lung + n50_phag)

    # Additional in vivo phagocytic capacity: one rate constant for uptake and
    # digestion, reduced in proportion to circ / circ0 by cyclophosphamide and
    # decaying with time since inoculation (Results; Data S2 KPGM = KDGM).
    km <- km_kill * (circ / circ0) * exp(-km_loss * tsince)

    bact_ext <- bact_susceptible + bact_resting
    ksr <- bact_ext * (kgrowth - kdeath) / bmax

    d/dt(bact_susceptible) <- -kdeath * bact_susceptible - kn_phag * bact_susceptible - km * bact_susceptible - ksr * bact_susceptible + kgrowth * bact_susceptible
    d/dt(bact_resting) <- -kdeath * bact_resting - kn_phag * bact_resting - km * bact_resting + ksr * bact_susceptible
    d/dt(bact_phagocytosed) <- kn_phag * bact_ext - kn_dig * bact_phagocytosed
    d/dt(bact_phagocytosed_mac) <- km * bact_ext - km * bact_phagocytosed_mac

    # Observations (Data S2 $ERROR): lung CFU counts every viable bacterium,
    # including those inside phagocytes; a non-positive total is floored at
    # 1 CFU/g before the log10.
    cfu_total <- bact_ext + bact_phagocytosed + bact_phagocytosed_mac
    if (cfu_total <= 0) cfu_total <- 1
    log_cfu <- log10(cfu_total)
    log_anc <- log10(circ)
    log_neut_lung <- log10(neut_lung * lung_weight)

    log_cfu ~ add(addSd_log_cfu)
    log_anc ~ add(addSd_log_anc)
    log_neut_lung ~ add(addSd_log_neut_lung)
  })
}
