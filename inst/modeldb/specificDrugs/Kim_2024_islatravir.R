Kim_2024_islatravir <- function() {
  description <- paste(
    "QSP (semi-mechanistic PK-PD with HIV-1 viral dynamics). Islatravir",
    "(ISL) for HIV-1 pre-exposure prophylaxis after oral tablets or a",
    "subcutaneous (subdermal) implant, from Kim et al. 2024. Plasma ISL is a",
    "three-compartment model (central plasma volume fixed at 3.5 L) fed",
    "either by first-order oral absorption or by an implant whose first-order",
    "release rate decays exponentially from insertion into a subcutaneous",
    "compartment that also has its own first-order loss; a loading-dose",
    "dependent fraction of the implant load is placed directly in the",
    "subcutaneous compartment at insertion. Intracellular ISL-triphosphate",
    "(ISL-TP) in peripheral blood mononuclear cells is driven by plasma ISL",
    "with a dose- and route-specific uptake/conversion rate and first-order",
    "loss. ISL-TP inhibits reverse transcription through an Emax model",
    "(IC50 = 429.707 nM) inside an established eight-state HIV-1 viral",
    "dynamics model (T cells and macrophages, uninfected / early infected /",
    "productively infected, plus infectious and non-infectious virus) that",
    "starts at its untreated chronic-infection steady state; the observation",
    "is the log10 change from baseline in total plasma virus. Also reports",
    "the prophylactic efficacy that the current ISL-TP level would give if",
    "held constant. Fitted to published mean concentration-time profiles",
    "by least squares: typical-value model with no IIV and no residual",
    "error."
  )
  reference <- paste(
    "Kim S, Zhang L, von Kleist M. Modeling of HIV-1 prophylactic efficacy",
    "and toxicity with islatravir shows non-superiority for oral dosing, but",
    "promise as a subcutaneous implant. CPT Pharmacometrics Syst Pharmacol.",
    "2024;13(10):1693-1705. doi:10.1002/psp4.13212. PMCID PMC11494919.",
    "Structure from Methods Equations (1)-(8) and Figure 1a; PK parameters",
    "from Table 1; dose-specific parameters from Supplementary Table S2 and",
    "the Figure S6 legend; IC50 from Results; prophylactic-efficacy",
    "propensities from Supplementary Text S1 Eqs. (S1)-(S6). Unit",
    "conventions and the viral-dynamics constants are from the authors'",
    "code deposit https://github.com/KleistLab/ISL_PK-PD-PrEP (GPL-3.0,",
    "commit abc6d6f), whose viral-dynamics model is that of Duwal S,",
    "Schuette C, von Kleist M. PLoS One. 2012;7(7):e40382.",
    "doi:10.1371/journal.pone.0040382."
  )
  vignette <- "Kim_2024_islatravir"

  # Doses and the three dosing-site states are in mg. Plasma ISL and the two
  # peripheral states are CONCENTRATION states in nM (Kim 2024 Eq. 4 writes
  # dC/dt for a concentration, with the oral input divided by VC). The
  # intracellular ISL-TP state is also in nM; Cpbmc_isltp converts it to the
  # clinically reported pmol per 10^6 cells. The viral-dynamics states are
  # cell and virion counts.
  units <- list(time = "h", dosing = "mg", concentration = "nM")

  # Three dosing targets: depot (oral tablet), implant (the implant's load)
  # and sc (the fraction of the load placed directly in the subcutaneous
  # tissue at insertion).
  dosing <- c("depot", "implant", "sc")

  paper_specific_compartments <- c(
    "implant",
    "sc",
    "pbmc_isltp",
    "tcell_uninfected",
    "tcell_early",
    "tcell_infected",
    "macro_uninfected",
    "macro_early",
    "macro_infected",
    "virus_noninfectious"
  )

  covariateData <- list(
    DOSE_ISLATRAVIR_MG = list(
      description = paste(
        "Nominal islatravir amount of the subject's regimen, in mg: the oral",
        "dose per administration for tablet regimens, or the implant loading",
        "dose for implant regimens."
      ),
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Selects the dose-specific parameters of Kim 2024 Supplementary",
        "Table S2. Oral (ROUTE_ORAL = 1): the plasma-to-PBMC uptake/conversion",
        "rate kC->E is looked up for the ten fitted doses (0.5, 1, 2, 5, 10,",
        "15, 30, 100, 200, 400 mg) and the three doses the paper extrapolated",
        "in the Figure S6 legend (0.1, 0.25, 0.75 mg); any other oral dose",
        "uses the Figure S6 power law of the weekly steady-state plasma",
        "maximum. Implant (ROUTE_ORAL = 0): kC->E and the release rate",
        "kImpl->SC are defined only for the 48, 52, 54, 56 and 62 mg loads,",
        "and the initial subcutaneous fraction pi is 4.12% for loads of",
        "54 mg or less and 5.83% above; any other implant load returns NA",
        "for kC->E and kImpl->SC. A value of 0 (an untreated run) gives",
        "kC->E = 0. Constant within a subject; carry it on every record,",
        "placed after the event columns."
      ),
      source_name = "ld (implant loading dose); dose in mg (oral)"
    ),
    ROUTE_ORAL = list(
      description = paste(
        "1 = oral islatravir tablet regimen, 0 = subcutaneous islatravir",
        "implant regimen."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (subcutaneous implant)",
      notes = paste(
        "Selects which half of Supplementary Table S2 supplies kC->E",
        "(oral vs implant columns). The route-specific absorption",
        "parameters need no switch: ka acts only on depot and kSC->C /",
        "kSC->0 only on sc, so the dose record's cmt carries the rest of",
        "the route. Not a data column in the source, which fitted each",
        "route and dose as its own scenario."
      ),
      source_name = NA_character_
    )
  )

  compartmentData <- list(
    depot = list(analyte = "islatravir", units = "mg", specimen = "administration site", verified = TRUE),
    implant = list(analyte = "islatravir", units = "mg", specimen = "administration site", verified = TRUE),
    sc = list(analyte = "islatravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "islatravir", units = "nM", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "islatravir", units = "nM", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "islatravir", units = "nM", specimen = "plasma", verified = TRUE),
    pbmc_isltp = list(
      analyte = "islatravir triphosphate",
      units = "nM",
      specimen = "blood cell",
      verified = TRUE
    ),
    tcell_uninfected = list(
      analyte = "uninfected CD4+ T cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_early = list(
      analyte = "early infected CD4+ T cells (pre-integration)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    tcell_infected = list(
      analyte = "productively infected CD4+ T cells",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    macro_uninfected = list(
      analyte = "uninfected macrophages",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    macro_early = list(
      analyte = "early infected macrophages (pre-integration)",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    macro_infected = list(
      analyte = "productively infected macrophages",
      units = "cells",
      specimen = "not applicable",
      verified = TRUE
    ),
    virus = list(analyte = "infectious HIV-1 virions", units = "virions", specimen = "plasma", verified = TRUE),
    virus_noninfectious = list(
      analyte = "non-infectious HIV-1 virions",
      units = "virions",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 114,
    n_studies = 5,
    age_range = "18-60 years (mean 45.4)",
    sex_female_pct = 27,
    disease_state = paste(
      "Adults without HIV (oral single doses 5-400 mg and all implant",
      "cohorts) and treatment-naive adults with HIV-1 (oral single doses",
      "0.5-30 mg, also the viral-load data)"
    ),
    dose_range = paste(
      "Oral single doses 0.5-400 mg; subdermal implants loaded with 48, 52,",
      "54, 56 or 62 mg, removed at week 12"
    ),
    regions = "Not reported (phase 1/2 trials sponsored by Merck)",
    notes = paste(
      "Model fitted to digitised published MEAN profiles, not individual",
      "data (Supplementary Table S1): 10 oral plasma ISL and 11 oral PBMC",
      "ISL-TP datasets (Schuermann 2020, Matthews 2021 CTS), implant plasma",
      "ISL and ISL-TP at 54 and 62 mg (Matthews 2021 Nat Med) and implant",
      "ISL-TP at 48, 52 and 56 mg (Matthews 2022 JAIDS). About 83 of 114",
      "individuals were male; mean BMI 26.1 (17.1-32) kg/m2; demographics",
      "were only partially available. IC50 fitted to mean viral-load",
      "declines after single oral doses of 0.5, 1, 2, 10 and 30 mg in 30",
      "treatment-naive men with HIV-1 (6 per dose)."
    )
  )

  ini({
    # ---- Plasma ISL PK (Table 1; rate constants in 1/h) -------------------
    lka <- log(0.4232); label("Oral first-order absorption rate constant ka (1/h)") # Table 1, ka = 0.4232 (oral)
    lksc_c <- log(0.0151); label("Subcutaneous-to-plasma transfer rate constant kSC->C (1/h)") # Table 1, kSC->C = 0.0151 (implant)
    lksc_0 <- log(0.1521); label("Subcutaneous loss rate constant kSC->0 (1/h)") # Table 1, kSC->0 = 0.1521 (implant)
    lkel <- log(6.8999); label("Elimination rate constant from plasma kC->0 (1/h)") # Table 1, kC->0 = 6.8999
    lk12 <- log(1.8494); label("Plasma to peripheral I rate constant kC->P1 (1/h)") # Table 1, kC->P1 = 1.8494
    lk21 <- log(0.0023); label("Peripheral I to plasma rate constant kP1->C (1/h)") # Table 1, kP1->C = 0.0023
    lk13 <- log(5.2726); label("Plasma to peripheral II rate constant kC->P2 (1/h)") # Table 1, kC->P2 = 5.2726
    lk31 <- log(0.0237); label("Peripheral II to plasma rate constant kP2->C (1/h)") # Table 1, kP2->C = 0.0237
    lvc <- fixed(log(3.5)); label("Central (plasma) volume VC (L)") # Table 1, VC = 3.5 L (fixed to physiological plasma volume)
    lkel_pbmc <- log(0.0088); label("Elimination rate constant of intracellular ISL-TP kE->0 (1/h)") # Table 1, kE->0 = 0.0088

    # ---- Implant release (Supplementary Table S2) ---------------------------
    lkr <- log(2.18e-4); label("Decay rate constant of the implant release rate kR (1/h)") # Table S2, kR = 2.18e-4 (global)
    lkimpl_48 <- log(1.92e-4); label("Implant release rate constant, 48 mg load (1/h)") # Table S2, kImpl->SC(48) = 1.92e-4
    lkimpl_52 <- log(2.77e-4); label("Implant release rate constant, 52 mg load (1/h)") # Table S2, kImpl->SC(52) = 2.77e-4
    lkimpl_54 <- log(3.07e-4); label("Implant release rate constant, 54 mg load (1/h)") # Table S2, kImpl->SC(54) = 3.07e-4
    lkimpl_56 <- log(3.45e-4); label("Implant release rate constant, 56 mg load (1/h)") # Table S2, kImpl->SC(56) = 3.45e-4
    lkimpl_62 <- log(4.60e-4); label("Implant release rate constant, 62 mg load (1/h)") # Table S2, kImpl->SC(62) = 4.60e-4
    lfsc_le54 <- log(0.0412); label("Fraction of the implant load placed in sc at insertion, loads <= 54 mg (fraction)") # Table S2, pi48 = pi52 = pi54 = 4.12 %
    lfsc_gt54 <- log(0.0583); label("Fraction of the implant load placed in sc at insertion, loads > 54 mg (fraction)") # Table S2, pi56 = pi62 = 5.83 %

    # ---- Plasma-to-PBMC uptake and conversion kC->E (1/h) ------------------
    lkce_oral_0p1 <- log(56.64); label("kC->E for 0.1 mg oral (1/h)") # Figure S6 legend, kC->E(0.1) = 56.64 (extrapolated by the authors)
    lkce_oral_0p25 <- log(40.32); label("kC->E for 0.25 mg oral (1/h)") # Figure S6 legend, kC->E(0.25) = 40.32 (extrapolated by the authors)
    lkce_oral_0p5 <- log(37.69); label("kC->E for 0.5 mg oral (1/h)") # Table S2, kC->E (oral) 0.5 mg = 37.69
    lkce_oral_0p75 <- log(26.82); label("kC->E for 0.75 mg oral (1/h)") # Figure S6 legend, kC->E(0.75) = 26.82 (extrapolated by the authors)
    lkce_oral_1 <- log(28.17); label("kC->E for 1 mg oral (1/h)") # Table S2, kC->E (oral) 1 mg = 28.17
    lkce_oral_2 <- log(16.56); label("kC->E for 2 mg oral (1/h)") # Table S2, kC->E (oral) 2 mg = 16.56
    lkce_oral_5 <- log(12.59); label("kC->E for 5 mg oral (1/h)") # Table S2, kC->E (oral) 5 mg = 12.59
    lkce_oral_10 <- log(17.46); label("kC->E for 10 mg oral (1/h)") # Table S2, kC->E (oral) 10 mg = 17.46
    lkce_oral_15 <- log(13.81); label("kC->E for 15 mg oral (1/h)") # Table S2, kC->E (oral) 15 mg = 13.81
    lkce_oral_30 <- log(16.95); label("kC->E for 30 mg oral (1/h)") # Table S2, kC->E (oral) 30 mg = 16.95
    lkce_oral_100 <- log(8.83); label("kC->E for 100 mg oral (1/h)") # Table S2, kC->E (oral) 100 mg = 08.83
    lkce_oral_200 <- log(10.57); label("kC->E for 200 mg oral (1/h)") # Table S2, kC->E (oral) 200 mg = 10.57
    lkce_oral_400 <- log(9.10); label("kC->E for 400 mg oral (1/h)") # Table S2, kC->E (oral) 400 mg = 09.10
    lkce_implant_48 <- log(25.30); label("kC->E for the 48 mg implant (1/h)") # Table S2, kC->E (implant) 48 mg = 25.30
    lkce_implant_52 <- log(32.55); label("kC->E for the 52 mg implant (1/h)") # Table S2, kC->E (implant) 52 mg = 32.55
    lkce_implant_54 <- log(19.35); label("kC->E for the 54 mg implant (1/h)") # Table S2, kC->E (implant) 54 mg = 19.35
    lkce_implant_56 <- log(41.98); label("kC->E for the 56 mg implant (1/h)") # Table S2, kC->E (implant) 56 mg = 41.98
    lkce_implant_62 <- log(36.87); label("kC->E for the 62 mg implant (1/h)") # Table S2, kC->E (implant) 62 mg = 36.87
    kce_s6_int <- 4.4; label("Figure S6 power law: log kC->E at a 1 nM steady-state plasma maximum (log 1/h)") # Figure S6 caption, kC->E = x^-0.371 * e^4.4
    kce_s6_exp <- -0.371; label("Figure S6 power law: exponent on the weekly steady-state plasma maximum (unitless)") # Figure S6 caption, kC->E = x^-0.371 * e^4.4

    # ---- Pharmacodynamics ---------------------------------------------------
    lic50 <- log(429.707); label("Intracellular ISL-TP concentration for 50% inhibition of reverse transcription IC50 (nM)") # Results 'Pharmacodynamics', IC50 = 429.707 nM

    # ---- HIV-1 viral dynamics (rates in 1/day, as in the source) -----------
    # Taken unchanged from the literature model (not estimated by Kim 2024);
    # values are those of params_vd in code/steady_state.py of the authors'
    # deposit, which generated Figure 2 and the IC50 fit.
    lambda_t <- fixed(2e9); label("Production of uninfected T cells lambdaT (cells/day)") # deposit steady_state.py 'lambdaT': 2e9; Duwal 2012 Table 2
    delta_t <- fixed(0.02); label("Death rate of uninfected T cells deltaT (1/day)") # deposit 'deltaT': 0.02; Duwal 2012 Table 2
    delta_t1 <- fixed(0.5); label("Death rate of early infected T cells deltaT1 (1/day)") # deposit 'deltaT1': 0.5
    delta_t2 <- fixed(1.4); label("Death rate of productively infected T cells deltaT2 (1/day)") # deposit 'deltaT2': 1.4
    delta_pic_t <- fixed(0.35); label("Intracellular destruction of the pre-integration complex in T cells deltaPIC,T (1/day)") # deposit 'deltaPICT': 0.35; Duwal 2012 Table 2
    k_t <- fixed(0.35); label("Integration (early to productive) rate in T cells kT (1/day)") # deposit 'kT': 0.35; Duwal 2012 Table 2
    beta_t <- fixed(8e-12); label("Infection rate of T cells without drug betaT (1/(virion day))") # deposit 'betaT0': 8e-12; Duwal 2012 Table 2
    n_t <- fixed(670); label("Infectious virions released per productively infected T cell NT (virions/(cell day))") # deposit 'NT': 670 = 0.67 x 1000 (Duwal 2012 Table 2)
    n_t_noninf <- fixed(330); label("Non-infectious virions released per productively infected T cell (virions/(cell day))") # deposit 'NTdiff': 330 = 1000 - 670
    lambda_m <- fixed(6.9e7); label("Production of uninfected macrophages lambdaM (cells/day)") # deposit 'lambdaM': 6.9e7; Duwal 2012 Table 2
    delta_m <- fixed(0.0069); label("Death rate of uninfected macrophages deltaM (1/day)") # deposit 'deltaM': 0.0069; Duwal 2012 Table 2
    delta_m1 <- fixed(0.0069); label("Death rate of early infected macrophages deltaM1 (1/day)") # deposit 'deltaM1': 0.0069; Duwal 2012 Table 2
    delta_m2 <- fixed(0.09); label("Death rate of productively infected macrophages deltaM2 (1/day)") # deposit 'deltaM2': 0.09; Duwal 2012 Table 2
    delta_pic_m <- fixed(0.0035); label("Intracellular destruction of the pre-integration complex in macrophages deltaPIC,M (1/day)") # deposit 'deltaPICM': 0.0035; Duwal 2012 Table 2
    k_m <- fixed(0.07); label("Integration (early to productive) rate in macrophages kM (1/day)") # deposit 'kM': 0.07; Duwal 2012 Table 2
    beta_m <- fixed(1e-13); label("Infection rate of macrophages without drug betaM (1/(virion day))") # deposit 'betaM0': 1e-13
    n_m <- fixed(67); label("Infectious virions released per productively infected macrophage NM (virions/(cell day))") # deposit 'NM': 67 = 0.67 x 100 (Duwal 2012 Table 2)
    n_m_noninf <- fixed(33); label("Non-infectious virions released per productively infected macrophage (virions/(cell day))") # deposit 'NMdiff': 33 = 100 - 67
    cl_virus <- fixed(23); label("Clearance rate of free virus in infected individuals CL (1/day)") # deposit 'CL': 23; Duwal 2012 Table 2 CL(infected)
    cl_virus_naive <- fixed(2.3); label("Clearance rate of free virus at exposure, prophylaxis setting (1/day)") # deposit pgs_ISL.py CL = 2.3; Duwal 2012 Table 2 CL(naive)
  })

  model({
    # ---- 1. Unit conversions -------------------------------------------------
    # ISL molar mass 293.258 g/mol (deposit pgs_ISL.py M_ISL), so 1 mg =
    # 1e6/293.258 nmol. PBMC ISL-TP: 1 pmol per 10^6 cells = 5555.6 nM
    # intracellular (deposit pk_optimization_bootstrap.py converts the
    # reported pmol/10^6 cells by x 5555.6, i.e. 0.18 pL per PBMC).
    nmol_per_mg <- 1e6 / 293.258
    nm_per_pmol_per_mcell <- 5555.6

    # ---- 2. Individual parameters (no IIV) ------------------------------------
    ka <- exp(lka)
    ksc_c <- exp(lksc_c)
    ksc_0 <- exp(lksc_0)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k13 <- exp(lk13)
    k31 <- exp(lk31)
    vc <- exp(lvc)
    kel_pbmc <- exp(lkel_pbmc)
    kr <- exp(lkr)
    ic50 <- exp(lic50)

    # ---- 3. Dose- and route-specific parameters (Table S2) --------------------
    # Exact-match lookup on the nominal amount; the tolerance only absorbs
    # floating-point representation of e.g. 0.1 mg.
    is_0p1 <- abs(DOSE_ISLATRAVIR_MG - 0.1) < 1e-6
    is_0p25 <- abs(DOSE_ISLATRAVIR_MG - 0.25) < 1e-6
    is_0p5 <- abs(DOSE_ISLATRAVIR_MG - 0.5) < 1e-6
    is_0p75 <- abs(DOSE_ISLATRAVIR_MG - 0.75) < 1e-6
    is_1 <- abs(DOSE_ISLATRAVIR_MG - 1) < 1e-6
    is_2 <- abs(DOSE_ISLATRAVIR_MG - 2) < 1e-6
    is_5 <- abs(DOSE_ISLATRAVIR_MG - 5) < 1e-6
    is_10 <- abs(DOSE_ISLATRAVIR_MG - 10) < 1e-6
    is_15 <- abs(DOSE_ISLATRAVIR_MG - 15) < 1e-6
    is_30 <- abs(DOSE_ISLATRAVIR_MG - 30) < 1e-6
    is_48 <- abs(DOSE_ISLATRAVIR_MG - 48) < 1e-6
    is_52 <- abs(DOSE_ISLATRAVIR_MG - 52) < 1e-6
    is_54 <- abs(DOSE_ISLATRAVIR_MG - 54) < 1e-6
    is_56 <- abs(DOSE_ISLATRAVIR_MG - 56) < 1e-6
    is_62 <- abs(DOSE_ISLATRAVIR_MG - 62) < 1e-6
    is_100 <- abs(DOSE_ISLATRAVIR_MG - 100) < 1e-6
    is_200 <- abs(DOSE_ISLATRAVIR_MG - 200) < 1e-6
    is_400 <- abs(DOSE_ISLATRAVIR_MG - 400) < 1e-6

    # Oral kC->E: tabulated for 13 doses. Any other oral dose uses the
    # Figure S6 power law kC->E = x^-0.371 * e^4.4, where x is the
    # steady-state plasma maximum (nM) of the same dose given once weekly.
    # Plasma PK is linear, so x = cmax_ss_weekly_per_mg * dose. 26.62 nM/mg
    # is the weekly steady-state Cmax this model gives per mg (26.618), and
    # the same value back-solves from each of the three kC->E values the
    # paper read off the power law (0.1, 0.25, 0.75 mg; Figure S6 legend).
    cmax_ss_weekly_per_mg <- 26.62
    oral_tab <- is_0p1 + is_0p25 + is_0p5 + is_0p75 + is_1 + is_2 + is_5 +
      is_10 + is_15 + is_30 + is_100 + is_200 + is_400
    kce_oral_tab <- is_0p1 * exp(lkce_oral_0p1) + is_0p25 * exp(lkce_oral_0p25) +
      is_0p5 * exp(lkce_oral_0p5) + is_0p75 * exp(lkce_oral_0p75) +
      is_1 * exp(lkce_oral_1) + is_2 * exp(lkce_oral_2) +
      is_5 * exp(lkce_oral_5) + is_10 * exp(lkce_oral_10) +
      is_15 * exp(lkce_oral_15) + is_30 * exp(lkce_oral_30) +
      is_100 * exp(lkce_oral_100) + is_200 * exp(lkce_oral_200) +
      is_400 * exp(lkce_oral_400)
    kce_oral_s6 <- 0
    if (DOSE_ISLATRAVIR_MG > 0) {
      kce_oral_s6 <- exp(kce_s6_int) *
        (cmax_ss_weekly_per_mg * DOSE_ISLATRAVIR_MG)^kce_s6_exp
    }
    kce_oral <- oral_tab * kce_oral_tab + (1 - oral_tab) * kce_oral_s6

    # Implant: kC->E and the release rate exist only for the five studied
    # loads; pi follows the paper's <= 54 mg / > 54 mg grouping.
    implant_tab <- is_48 + is_52 + is_54 + is_56 + is_62
    kce_implant <- is_48 * exp(lkce_implant_48) + is_52 * exp(lkce_implant_52) +
      is_54 * exp(lkce_implant_54) + is_56 * exp(lkce_implant_56) +
      is_62 * exp(lkce_implant_62)
    kimpl <- is_48 * exp(lkimpl_48) + is_52 * exp(lkimpl_52) +
      is_54 * exp(lkimpl_54) + is_56 * exp(lkimpl_56) +
      is_62 * exp(lkimpl_62)
    if (ROUTE_ORAL == 0 && implant_tab == 0) {
      kce_implant <- NA
      kimpl <- NA
    }
    fsc <- exp(lfsc_le54)
    if (DOSE_ISLATRAVIR_MG > 54) fsc <- exp(lfsc_gt54)

    kce <- ROUTE_ORAL * kce_oral + (1 - ROUTE_ORAL) * kce_implant

    # ---- 4. PK ODEs (Eqs. 1-7) --------------------------------------------
    # Eq. 1: the implant's first-order release rate decays as exp(-kR t)
    # from insertion. tad(implant) restarts the clock at every insertion;
    # removal is a replacement event (evid = 5, amt = 0) on implant, which
    # stops release as in the source (release set to 0 at week 12).
    t_impl <- tad(implant)
    if (is.na(t_impl)) t_impl <- 0
    krel <- kimpl * exp(-kr * t_impl)
    if (ROUTE_ORAL == 1) krel <- 0

    d/dt(depot) <- -ka * depot
    d/dt(implant) <- -krel * implant
    d/dt(sc) <- krel * implant - (ksc_c + ksc_0) * sc
    # Eq. 4: plasma concentration (nM). The oral input is divided by VC; the
    # subcutaneous input is not (Eq. 4 and the deposit write kSC->C * SC
    # with SC in nmol), so the implant parameters were fitted on that scale.
    d/dt(central) <- ka * depot * nmol_per_mg / vc +
      ksc_c * sc * nmol_per_mg +
      k21 * peripheral1 + k31 * peripheral2 -
      (k12 + k13 + kel) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    # Eq. 7: intracellular ISL-TP (nM).
    d/dt(pbmc_isltp) <- kce * central - kel_pbmc * pbmc_isltp

    # Initial fraction pi of the implant load starts in sc, the rest in the
    # implant: dose the full load to both implant and sc at insertion.
    f(implant) <- 1 - fsc
    f(sc) <- fsc

    # ---- 5. HIV-1 viral dynamics (rates in 1/day; divided by 24) ----------
    # Eq. 8: direct effect of ISL-TP on reverse transcription.
    eta_rt <- pbmc_isltp / (ic50 + pbmc_isltp)
    beta_t_drug <- beta_t * (1 - eta_rt)
    beta_m_drug <- beta_m * (1 - eta_rt)
    cl_t_drug <- (1 + eta_rt) * beta_t
    cl_m_drug <- (1 + eta_rt) * beta_m

    # Untreated steady state (the authors' baseline), in closed form: with
    # no drug, dT1 = dT2 = dM1 = dM2 = dV = 0 reduce to a quadratic in V.
    a_t <- delta_t1 + k_t + delta_pic_t
    a_m <- delta_m1 + k_m + delta_pic_m
    g_t <- beta_t * (delta_t1 + k_t) / a_t
    g_m <- beta_m * (delta_m1 + k_m) / a_m
    r_t <- beta_t * (n_t * k_t / (a_t * delta_t2) - 2)
    r_m <- beta_m * (n_m * k_m / (a_m * delta_m2) - 2)
    qa <- cl_virus * g_t * g_m
    qb <- cl_virus * (delta_t * g_m + delta_m * g_t) - r_t * lambda_t * g_m -
      r_m * lambda_m * g_t
    qc <- cl_virus * delta_t * delta_m - r_t * lambda_t * delta_m -
      r_m * lambda_m * delta_t
    v_ss <- (-qb + sqrt(qb^2 - 4 * qa * qc)) / (2 * qa)
    tu_ss <- lambda_t / (delta_t + g_t * v_ss)
    mu_ss <- lambda_m / (delta_m + g_m * v_ss)
    t1_ss <- beta_t * v_ss * tu_ss / a_t
    m1_ss <- beta_m * v_ss * mu_ss / a_m
    t2_ss <- k_t * t1_ss / delta_t2
    m2_ss <- k_m * m1_ss / delta_m2
    vn_ss <- (n_t_noninf * t2_ss + n_m_noninf * m2_ss) / cl_virus

    tcell_uninfected(0) <- tu_ss
    tcell_early(0) <- t1_ss
    tcell_infected(0) <- t2_ss
    macro_uninfected(0) <- mu_ss
    macro_early(0) <- m1_ss
    macro_infected(0) <- m2_ss
    virus(0) <- v_ss
    virus_noninfectious(0) <- vn_ss

    d/dt(tcell_uninfected) <- (lambda_t - delta_t * tcell_uninfected -
      beta_t_drug * virus * tcell_uninfected + delta_pic_t * tcell_early) / 24
    d/dt(macro_uninfected) <- (lambda_m - delta_m * macro_uninfected -
      beta_m_drug * virus * macro_uninfected + delta_pic_m * macro_early) / 24
    d/dt(tcell_early) <- (beta_t_drug * virus * tcell_uninfected -
      a_t * tcell_early) / 24
    d/dt(macro_early) <- (beta_m_drug * virus * macro_uninfected -
      a_m * macro_early) / 24
    d/dt(tcell_infected) <- (k_t * tcell_early - delta_t2 * tcell_infected) / 24
    d/dt(macro_infected) <- (k_m * macro_early - delta_m2 * macro_infected) / 24
    d/dt(virus) <- (n_m * macro_infected + n_t * tcell_infected -
      virus * (cl_virus + (beta_t_drug + cl_t_drug) * tcell_uninfected +
        (cl_m_drug + beta_m_drug) * macro_uninfected)) / 24
    d/dt(virus_noninfectious) <- (n_t_noninf * tcell_infected +
      n_m_noninf * macro_infected - cl_virus * virus_noninfectious) / 24

    # ---- 6. Prophylactic efficacy at a constant ISL-TP level (Text S1) -----
    # Extinction probabilities of the branching process with propensities
    # (S1)-(S6) at uninfected T-cell steady state, solved as the quadratic
    # of the authors' getExtinctionProbability(). This is the efficacy the
    # CURRENT intracellular level would give if held constant (the basis of
    # EC85 / EC90 in Figure S4), not the time-resolved efficacy of Eq. 9.
    tu_naive <- lambda_t / delta_t
    pe_a2 <- delta_pic_t + delta_t1
    pe_b1 <- pe_a2 / (pe_a2 + k_t)
    pe_b3 <- k_t / (pe_a2 + k_t)
    pe_g <- n_t / (delta_t2 + n_t)
    pe_ec <- 1 - pe_g
    # drug-free reference
    pe_a1_0 <- cl_virus_naive + beta_t * tu_naive
    pe_a4_0 <- beta_t * tu_naive
    pe_al3_0 <- pe_a4_0 / (pe_a1_0 + pe_a4_0)
    pe_ea_0 <- pe_al3_0 * pe_b3 * pe_g
    pe_eb_0 <- pe_g - 1 - pe_ea_0
    pe_t2_0 <- (-pe_eb_0 - sqrt(pe_eb_0^2 - 4 * pe_ea_0 * pe_ec)) / (2 * pe_ea_0)
    pe_none <- (1 - pe_al3_0) + pe_al3_0 * (pe_b1 + pe_b3 * pe_t2_0)
    # with the current ISL-TP level
    pe_a1 <- cl_virus_naive + (1 + eta_rt) * beta_t * tu_naive
    pe_a4 <- (1 - eta_rt) * beta_t * tu_naive
    pe_al3 <- pe_a4 / (pe_a1 + pe_a4)
    pe_ea <- pe_al3 * pe_b3 * pe_g
    pe_eb <- pe_g - 1 - pe_ea
    pe_t2 <- (-pe_eb - sqrt(pe_eb^2 - 4 * pe_ea * pe_ec)) / (2 * pe_ea)
    pe_drug <- (1 - pe_al3) + pe_al3 * (pe_b1 + pe_b3 * pe_t2)
    phi_const <- 1 - (1 - pe_drug) / (1 - pe_none)

    # ---- 7. Observations ---------------------------------------------------
    Cc <- central
    Cpbmc_isltp <- pbmc_isltp / nm_per_pmol_per_mcell
    vl_total0 <- v_ss + vn_ss
    d_viral_load <- log10((virus + virus_noninfectious) / vl_total0)
  })
}
