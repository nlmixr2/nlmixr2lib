Shimizu_2021_lusutrombopag <- function() {
  description <- paste(
    "QSP. Thrombopoiesis and platelet life-cycle model with thrombopoietin",
    "(TPO) target-mediated disposition, driven by the oral TPO-receptor",
    "agonist lusutrombopag, in healthy adults and in thrombocytopenic",
    "patients with chronic liver disease (CLD). 27 dividing megakaryocyte",
    "progenitor compartments (26 doublings), 5 megakaryocyte maturation",
    "compartments, 9 platelet aging compartments, a 3-state TPO binding",
    "model and a 3-compartment first-order-absorption lusutrombopag PK",
    "model. TPO and unbound lusutrombopag add on one shared Emax term that",
    "scales the progenitor division rate; splenic sequestration sets the",
    "plasma fraction of the total platelet pool.",
    sep = " "
  )
  reference <- paste(
    "Shimizu R, Katsube T, Wajima T. Quantitative systems pharmacology model",
    "of thrombopoiesis and platelet life-cycle, and its application to",
    "thrombocytopenia based on chronic liver disease. CPT Pharmacometrics",
    "Syst Pharmacol. 2021;10(5):489-499. doi:10.1002/psp4.12623.",
    "Model equations from Supplementary Text S1/S2 and the deposited MATLAB",
    "model code (Supplementary Model Code). Lusutrombopag PK from Katsube T",
    "et al. Clin Pharmacokinet. 2016;55:1423-1433; TPO binding kinetics from",
    "Jin F, Krzyzanski W. AAPS PharmSci. 2004;6:E9.",
    sep = " "
  )
  vignette <- "Shimizu_2021_lusutrombopag"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DIS_HEALTHY = list(
      description = "Population indicator: 1 = healthy adult, 0 = thrombocytopenic patient with chronic liver disease",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the population-specific parameter set of Shimizu 2021.",
        "DIS_HEALTHY = 1 uses the healthy-subject baselines (TPO0 1.4 pM,",
        "PLT0 20 x 10^4/uL) with splenic sequestration fixed at 1/3 and the",
        "platelets-per-megakaryocyte count derived from the steady state",
        "(Supplementary Text S1; deposited Healthy/ scripts).",
        "DIS_HEALTHY = 0 uses the CLD baselines (TPO0 0.78 pM, PLT0",
        "4 x 10^4/uL) with 2500 platelets per megakaryocyte and splenic",
        "sequestration derived from the steady state (deposited CLD/",
        "scripts; 78.6% at typical values). The complement cohort is",
        "Japanese patients with chronic liver disease and a platelet count",
        "below 5 x 10^4/uL (lusutrombopag phase II).",
        sep = " "
      ),
      source_name = "none (the paper describes the two populations as separate parameter sets, Table 1)"
    )
  )

  # Model-local state names. `depot`, `central`, `peripheral1`,
  # `peripheral2` and `precursor1..27` are canonical. The megakaryocyte and
  # platelet aging chains reuse the `mk<N>` / `plt<N>` model-local pattern of
  # Cao_2025_ferricCarboxymaltose_rat; the three thrombopoietin states are
  # specific to the Jin & Krzyzanski 2004 TPO binding model.
  paper_specific_compartments <- c("tpo", "tpo_ns", "tpo_complex")
  paper_specific_compartment_pattern <- "^(mk|plt)[0-9]+$"
  # The two baseline etas act on a population-selected typical value
  # (lrbase_plt_healthy or lrbase_plt_cld), so their names carry the shared
  # stem rather than either ini() parameter.
  paper_specific_etas <- c("etalrbase_plt", "etalrbase_tpo")

  compartmentData <- list(
    depot = list(analyte = "lusutrombopag", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lusutrombopag", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lusutrombopag", units = "mg", specimen = "tissue", verified = TRUE),
    peripheral2 = list(analyte = "lusutrombopag", units = "mg", specimen = "tissue", verified = TRUE),
    tpo = list(analyte = "thrombopoietin (free)", units = "pM", specimen = "plasma", verified = TRUE),
    tpo_ns = list(
      analyte = "thrombopoietin (nonspecifically bound)",
      units = "pM",
      specimen = "tissue",
      verified = TRUE
    ),
    tpo_complex = list(
      analyte = "thrombopoietin-c-Mpl receptor complex on platelets",
      units = "pM",
      specimen = "blood cell",
      verified = TRUE
    ),
    precursor1 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor2 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor3 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor4 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor5 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor6 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor7 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor8 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor9 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor10 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor11 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor12 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor13 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor14 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor15 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor16 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor17 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor18 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor19 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor20 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor21 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor22 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor23 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor24 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor25 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor26 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    precursor27 = list(
      analyte = "megakaryocyte progenitor cells (bone marrow)",
      units = "cells",
      specimen = "tissue",
      verified = TRUE
    ),
    mk1 = list(analyte = "megakaryocytes (bone marrow)", units = "cells", specimen = "tissue", verified = TRUE),
    mk2 = list(analyte = "megakaryocytes (bone marrow)", units = "cells", specimen = "tissue", verified = TRUE),
    mk3 = list(analyte = "megakaryocytes (bone marrow)", units = "cells", specimen = "tissue", verified = TRUE),
    mk4 = list(analyte = "megakaryocytes (bone marrow)", units = "cells", specimen = "tissue", verified = TRUE),
    mk5 = list(analyte = "megakaryocytes (bone marrow)", units = "cells", specimen = "tissue", verified = TRUE),
    plt1 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt2 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt3 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt4 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt5 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt6 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt7 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt8 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    ),
    plt9 = list(
      analyte = "platelets (whole body, circulating plus splenic pool)",
      units = "platelets",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 4,
    age_range = "adults (not tabulated in the paper)",
    disease_state = paste(
      "Healthy adults (Japanese and non-Japanese) and Japanese patients with",
      "chronic liver disease and thrombocytopenia",
      sep = " "
    ),
    dose_range = paste(
      "Lusutrombopag oral: single 1, 2, 4, 10, 25 or 50 mg (healthy Japanese);",
      "1 mg once daily for 14 days (healthy non-Japanese); 2 mg once daily for",
      "14 days (healthy Japanese); 3 mg once daily for 7 days (Japanese CLD",
      "patients, phase II)",
      sep = " "
    ),
    regions = "Japan; non-Japanese healthy-subject study (Table S1)",
    notes = paste(
      "The model was NOT fitted to data: the structure and parameters were",
      "assembled from physiology and the literature and the simulations were",
      "compared against observed platelet counts from the lusutrombopag",
      "clinical studies summarised in Supplementary Table S1 (phase I",
      "single- and multiple-dose studies in healthy subjects, and the phase",
      "II study in CLD patients). For the CLD simulation (Figure 4), 200",
      "virtual patients were generated by resampling patient demographics",
      "from the clinical data with 40% interindividual variability on PLT0,",
      "TPO0, PP and kout. Subject counts per study are in Table S1 and are",
      "not needed by the model.",
      sep = " "
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Lusutrombopag PK (Table 1; carried from Katsube 2016, ref 13).
    # Table 1 prints two-significant-figure roundings; the deposited
    # initial_values_constants_lusu.m holds the values actually simulated,
    # used here.
    # -------------------------------------------------------------------
    lka <- fixed(log(7.176)); label("First-order absorption rate constant, ka (1/day)") # Table 1 ka = 7.2 /day; deposited code ka = 7.176
    lkel <- fixed(log(1.358)); label("First-order elimination rate constant, ke (1/day)") # Table 1 ke = 1.4 /day; deposited code ke = 1.358
    lk12 <- fixed(log(1.2816)); label("Central to peripheral 1 rate constant, k12 (1/day)") # Table 1 k12 = 1.3 /day; deposited code 1.2816
    lk21 <- fixed(log(2.088)); label("Peripheral 1 to central rate constant, k21 (1/day)") # Table 1 k21 = 2.1 /day; deposited code 2.088
    lk13 <- fixed(log(0.03408)); label("Central to peripheral 2 rate constant, k13 (1/day)") # Table 1 k13 = 0.034 /day; deposited code 0.03408
    lk31 <- fixed(log(0.14088)); label("Peripheral 2 to central rate constant, k31 (1/day)") # Table 1 k31 = 0.14 /day; deposited code 0.14088
    lvc <- fixed(log(13.7)); label("Central volume of distribution, Vc (L)") # Table 1 Vc = 13.7 L; deposited code V = 13.7

    # -------------------------------------------------------------------
    # Emax stimulation of progenitor division (Equation 12 / Text S1)
    # -------------------------------------------------------------------
    lemax <- fixed(log(4.52)); label("Maximum effect via the TPO receptor, Emax (unitless)") # Table 1 Emax = 4.52 (from Katsube 2019 PK/PD, ref 14)
    lec50 <- fixed(log(183)); label("Lusutrombopag concentration at 50% of Emax, EC50Lusu (ng/mL)") # Table 1 EC50Lusu = 183 ng/mL (in vitro CD34+ cells, ref 25)
    # Equation 11: EC50TPO = (Emax - ETPO,ss) * TPO0 / ETPO,ss with ETPO,ss = 1
    # and the healthy TPO0 = 1.4 pM, i.e. (4.52 - 1) * 1.4 = 4.928 pM. The
    # deposited code hard-codes 4.928 for both populations.
    lec50_tpo <- fixed(log(4.928)); label("TPO concentration at 50% of Emax, EC50TPO (pM)") # Table 1 EC50TPO = 4.9 pM (footnote b, Equation 11); deposited code 4.928

    # -------------------------------------------------------------------
    # Thrombopoietin binding model (Jin & Krzyzanski 2004, ref 10)
    # -------------------------------------------------------------------
    kon <- fixed(1.32); label("TPO binding rate to platelet c-Mpl receptor, kon (1/pM/day)") # Table 1 kon = 1.3; deposited code kon = 1.32
    koff <- fixed(60); label("TPO dissociation rate from the receptor, koff (1/day)") # Table 1 koff = 60; deposited code koff = 60
    knf <- fixed(1.152); label("TPO nonspecific-binding return rate, knf (1/day)") # Table 1 knf = 1.2; deposited code knf = 1.152
    kfn <- fixed(3.12); label("TPO nonspecific-binding rate, kfn (1/day)") # Table 1 kfn = 3.1; deposited code kfn = 3.12
    # kint is named in Text S1 (the TPO-receptor complex internalization rate)
    # but no value is printed in the paper; the deposited healthy script sets
    # kint = 2.4 and the CLD script kint = 2.4 * 0.78 / 1.4 (scaled by TPO0).
    kint <- fixed(2.4); label("TPO-receptor complex internalization rate at healthy TPO0, kint (1/day)") # deposited Healthy/initial_values_constants_for_platelet.m kint = 2.4
    # The deposited CLD script scales Rp,0 by PLT0: mpl = 164 * PLT0 / 20.
    rp0 <- fixed(164); label("TPO receptor concentration on platelets at healthy PLT0, Rp,0 (pM)") # Table 1 Rp,0 = 164 pM; deposited code mpl = 164

    # -------------------------------------------------------------------
    # Baselines (Table 1; 'in house data')
    # -------------------------------------------------------------------
    lrbase_tpo_healthy <- fixed(log(1.4)); label("Baseline plasma TPO, healthy subjects, TPO0 (pM)") # Table 1 TPO0 healthy = 1.4 pM
    lrbase_tpo_cld <- fixed(log(0.78)); label("Baseline plasma TPO, CLD patients, TPO0 (pM)") # Table 1 TPO0 CLD = 0.78 pM
    lrbase_plt_healthy <- fixed(log(20)); label("Baseline plasma platelet count, healthy subjects, PLT0 (x10^4/uL)") # Table 1 PLT0 healthy = 20 x 10,000/uL
    lrbase_plt_cld <- fixed(log(4)); label("Baseline plasma platelet count, CLD patients, PLT0 (x10^4/uL)") # Table 1 PLT0 CLD = 4 x 10,000/uL

    # -------------------------------------------------------------------
    # Platelet production
    # -------------------------------------------------------------------
    lkout <- fixed(log(1)); label("Cell division / maturation / aging rate constant without TPO effect, kout (1/day)") # Table 1 kout = 1 /day
    # Used for CLD patients only; for healthy subjects PP is derived from the
    # steady state (Text S1: PP = PLT0 / %SPS / (2^26 * 9 / blood volume)).
    lpp <- fixed(log(2500)); label("Platelets produced per megakaryocyte in CLD patients, PP (count)") # Table 1 PP = 2,500; deposited CLD code PP = 2500
    # Text S1: '%SPS ... (fixed as 1/3 in healthy subjects)'; deposited
    # healthy code spleen = 2/3 is the complementary plasma fraction.
    sps_healthy <- fixed(1 / 3); label("Splenic platelet sequestration fraction in healthy subjects, %SPS (fraction)") # Text S1 %SPS = 1/3; Discussion 33.3%

    # -------------------------------------------------------------------
    # Interindividual variability used for the CLD simulation (Figure 4):
    # 'IIV for PLT0, TPO0, PP, and kout were set at 40% as arbitrary values'.
    # Encoded as log-normal with omega^2 = log(1 + 0.40^2) = 0.1484.
    # -------------------------------------------------------------------
    etalrbase_plt ~ fixed(0.1484) # Methods 'IIV ... set at 40%'; omega^2 = log(1 + 0.4^2)
    etalrbase_tpo ~ fixed(0.1484) # Methods 'IIV ... set at 40%'; omega^2 = log(1 + 0.4^2)
    etalpp ~ fixed(0.1484) # Methods 'IIV ... set at 40%'; omega^2 = log(1 + 0.4^2)
    etalkout ~ fixed(0.1484) # Methods 'IIV ... set at 40%'; omega^2 = log(1 + 0.4^2)
  })

  model({
    # Structural constants (Methods, 'Assumptions for platelet model
    # development'): 26 progenitor doublings, 9 one-day platelet aging
    # compartments, and a 5 L blood volume expressed in uL * 10^4 so that
    # whole-body platelet counts divide to the x10^4/uL reporting unit
    # (deposited code divisor 50000000000 = 5e10).
    nprolif <- 26
    nplt_cmt <- 9
    blood_volume <- 5e10

    # 1. Population-selected typical values with IIV
    rbase_tpo <- exp(DIS_HEALTHY * lrbase_tpo_healthy + (1 - DIS_HEALTHY) * lrbase_tpo_cld + etalrbase_tpo)
    rbase_plt <- exp(DIS_HEALTHY * lrbase_plt_healthy + (1 - DIS_HEALTHY) * lrbase_plt_cld + etalrbase_plt)
    kout <- exp(lkout + etalkout)
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    ec50_tpo <- exp(lec50_tpo)

    # Deposited CLD script: kint = 2.4 * TPO0 / 1.4 and mpl = 164 * PLT0 / 20,
    # i.e. both scale with the population baseline relative to the healthy
    # typical value (identity for healthy subjects).
    kint_i <- kint * rbase_tpo / exp(lrbase_tpo_healthy)
    rp0_i <- rp0 * rbase_plt / exp(lrbase_plt_healthy)

    # Baseline progenitor division rate relative to kout (Text S2): the TPO
    # Emax term at TPO0. Equals 1 for healthy subjects by construction of
    # EC50TPO (Equation 11) and 0.618 for CLD patients (Methods: 'kout1 was
    # restricted to be 0.62/day or higher').
    kout1_0 <- emax * rbase_tpo / (ec50_tpo + rbase_tpo)

    # Platelets per megakaryocyte and plasma (non-sequestered) fraction.
    # Healthy: %SPS fixed at 1/3 and PP derived so PLT(0) = PLT0.
    # CLD: PP = 2500 and %SPS derived so PLT(0) = PLT0 (Equation 10).
    pp_healthy <- rbase_plt * blood_volume / ((1 - sps_healthy) * 2^nprolif * nplt_cmt * kout1_0)
    pp_cld <- exp(lpp + etalpp)
    pp <- DIS_HEALTHY * pp_healthy + (1 - DIS_HEALTHY) * pp_cld
    fplasma_cld <- rbase_plt * blood_volume / (pp_cld * 2^nprolif * nplt_cmt * kout1_0)
    fplasma <- DIS_HEALTHY * (1 - sps_healthy) + (1 - DIS_HEALTHY) * fplasma_cld
    sps <- 1 - fplasma

    # 2. Lusutrombopag PK (Text S1 X45-X48)
    ka <- exp(lka)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    vc <- exp(lvc)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2 - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc

    # 3. Plasma platelet count (Equation 9): whole-body platelets over the
    # blood volume, times the non-sequestered fraction (1 - %SPS).
    plt_total <- plt1 + plt2 + plt3 + plt4 + plt5 + plt6 + plt7 + plt8 + plt9
    PLT <- plt_total / blood_volume * fplasma

    # 4. Progenitor division rate (Equation 12 / Text S1). The unbound
    # fraction fuLusu multiplies both the concentration and EC50Lusu and
    # cancels; the deposited code omits it. kout1 is floored at its baseline
    # value (Methods: 'restricted to be the initial value ... or higher').
    kout1_raw <- kout * emax * (tpo / (ec50_tpo + tpo) + Cc / (ec50 + Cc))
    kout1_floor <- kout * kout1_0
    kout1 <- max(kout1_raw, kout1_floor)

    # 5. Thrombopoietin binding model (Text S1 X42-X44). The receptor pool
    # scales with the plasma platelet count relative to its baseline.
    d/dt(tpo) <- kint_i * (kon * rbase_tpo * rp0_i / (koff + kint_i)) - kon * tpo * rp0_i * (PLT / rbase_plt) + koff * tpo_complex - kfn * tpo + knf * tpo_ns
    d/dt(tpo_ns) <- kfn * tpo - knf * tpo_ns
    d/dt(tpo_complex) <- kon * tpo * rp0_i * (PLT / rbase_plt) - (koff + kint_i) * tpo_complex

    tpo(0) <- rbase_tpo
    tpo_ns(0) <- kfn / knf * rbase_tpo
    tpo_complex(0) <- kon * rbase_tpo * rp0_i / (koff + kint_i)

    # 6. Megakaryocyte progenitor proliferation (Text S1 X1-X27): each
    # compartment divides once per 1/kout1, doubling the cell count.
    d/dt(precursor1) <- kout1 - kout1 * precursor1
    d/dt(precursor2) <- 2 * kout1 * precursor1 - kout1 * precursor2
    d/dt(precursor3) <- 2 * kout1 * precursor2 - kout1 * precursor3
    d/dt(precursor4) <- 2 * kout1 * precursor3 - kout1 * precursor4
    d/dt(precursor5) <- 2 * kout1 * precursor4 - kout1 * precursor5
    d/dt(precursor6) <- 2 * kout1 * precursor5 - kout1 * precursor6
    d/dt(precursor7) <- 2 * kout1 * precursor6 - kout1 * precursor7
    d/dt(precursor8) <- 2 * kout1 * precursor7 - kout1 * precursor8
    d/dt(precursor9) <- 2 * kout1 * precursor8 - kout1 * precursor9
    d/dt(precursor10) <- 2 * kout1 * precursor9 - kout1 * precursor10
    d/dt(precursor11) <- 2 * kout1 * precursor10 - kout1 * precursor11
    d/dt(precursor12) <- 2 * kout1 * precursor11 - kout1 * precursor12
    d/dt(precursor13) <- 2 * kout1 * precursor12 - kout1 * precursor13
    d/dt(precursor14) <- 2 * kout1 * precursor13 - kout1 * precursor14
    d/dt(precursor15) <- 2 * kout1 * precursor14 - kout1 * precursor15
    d/dt(precursor16) <- 2 * kout1 * precursor15 - kout1 * precursor16
    d/dt(precursor17) <- 2 * kout1 * precursor16 - kout1 * precursor17
    d/dt(precursor18) <- 2 * kout1 * precursor17 - kout1 * precursor18
    d/dt(precursor19) <- 2 * kout1 * precursor18 - kout1 * precursor19
    d/dt(precursor20) <- 2 * kout1 * precursor19 - kout1 * precursor20
    d/dt(precursor21) <- 2 * kout1 * precursor20 - kout1 * precursor21
    d/dt(precursor22) <- 2 * kout1 * precursor21 - kout1 * precursor22
    d/dt(precursor23) <- 2 * kout1 * precursor22 - kout1 * precursor23
    d/dt(precursor24) <- 2 * kout1 * precursor23 - kout1 * precursor24
    d/dt(precursor25) <- 2 * kout1 * precursor24 - kout1 * precursor25
    d/dt(precursor26) <- 2 * kout1 * precursor25 - kout1 * precursor26
    d/dt(precursor27) <- 2 * kout1 * precursor26 - kout1 * precursor27

    # Megakaryocyte maturation and marrow reservoir (Text S1 X28-X32): the
    # thrombopoietin-stimulated rate feeds MK1, then a fixed kout transit.
    d/dt(mk1) <- kout1 * precursor27 - kout * mk1
    d/dt(mk2) <- kout * mk1 - kout * mk2
    d/dt(mk3) <- kout * mk2 - kout * mk3
    d/dt(mk4) <- kout * mk3 - kout * mk4
    d/dt(mk5) <- kout * mk4 - kout * mk5

    # Platelet life-cycle (Text S1 X33-X41): pp platelets shed per
    # megakaryocyte, then nine one-day aging compartments.
    d/dt(plt1) <- pp * kout * mk5 - kout * plt1
    d/dt(plt2) <- kout * plt1 - kout * plt2
    d/dt(plt3) <- kout * plt2 - kout * plt3
    d/dt(plt4) <- kout * plt3 - kout * plt4
    d/dt(plt5) <- kout * plt4 - kout * plt5
    d/dt(plt6) <- kout * plt5 - kout * plt6
    d/dt(plt7) <- kout * plt6 - kout * plt7
    d/dt(plt8) <- kout * plt7 - kout * plt8
    d/dt(plt9) <- kout * plt8 - kout * plt9

    # Initial conditions (Text S2): the system starts at its own steady state.
    precursor1(0) <- 2^0
    precursor2(0) <- 2^1
    precursor3(0) <- 2^2
    precursor4(0) <- 2^3
    precursor5(0) <- 2^4
    precursor6(0) <- 2^5
    precursor7(0) <- 2^6
    precursor8(0) <- 2^7
    precursor9(0) <- 2^8
    precursor10(0) <- 2^9
    precursor11(0) <- 2^10
    precursor12(0) <- 2^11
    precursor13(0) <- 2^12
    precursor14(0) <- 2^13
    precursor15(0) <- 2^14
    precursor16(0) <- 2^15
    precursor17(0) <- 2^16
    precursor18(0) <- 2^17
    precursor19(0) <- 2^18
    precursor20(0) <- 2^19
    precursor21(0) <- 2^20
    precursor22(0) <- 2^21
    precursor23(0) <- 2^22
    precursor24(0) <- 2^23
    precursor25(0) <- 2^24
    precursor26(0) <- 2^25
    precursor27(0) <- 2^26
    mk1(0) <- 2^nprolif * kout1_0
    mk2(0) <- 2^nprolif * kout1_0
    mk3(0) <- 2^nprolif * kout1_0
    mk4(0) <- 2^nprolif * kout1_0
    mk5(0) <- 2^nprolif * kout1_0
    plt1(0) <- pp * 2^nprolif * kout1_0
    plt2(0) <- pp * 2^nprolif * kout1_0
    plt3(0) <- pp * 2^nprolif * kout1_0
    plt4(0) <- pp * 2^nprolif * kout1_0
    plt5(0) <- pp * 2^nprolif * kout1_0
    plt6(0) <- pp * 2^nprolif * kout1_0
    plt7(0) <- pp * 2^nprolif * kout1_0
    plt8(0) <- pp * 2^nprolif * kout1_0
    plt9(0) <- pp * 2^nprolif * kout1_0
  })
}
