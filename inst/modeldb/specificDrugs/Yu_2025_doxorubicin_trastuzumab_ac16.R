Yu_2025_doxorubicin_trastuzumab_ac16 <- function() {
  description <- "In vitro (AC16 human cardiomyocyte cell line). QSP / quantitative systems toxicology model of doxorubicin (DOX), trastuzumab (TmAb) and DOX + TmAb cardiotoxicity in AC16 cells: an apoptosis-signalling protein-dynamics layer (caspase-9 driven by a TmAb-stimulated upstream compartment and by DOX, feeding a three-compartment caspase-3 transit chain) drives a cardiomyocyte-death layer (a Hill function of the caspase-3 rise from baseline, delayed through two transit compartments) acting on an exponentially growing relative-viability state whose growth rate DOX also inhibits linearly. DOX degrades first-order in the culture medium; the TmAb concentration is held constant over the experiment. Observed outputs are relative caspase-9 activity (cas92), relative caspase-3 activity (cas33) and relative cell viability (viability), each normalised to its value at time 0. Random effects are inter-treatment-group variability (one 'subject' = one treatment group), on the DOX and TmAb caspase-9 slopes (protein-dynamics fit) and on the DOX growth-inhibition slope (viability fit). Sibling models from the same paper: Yu_2025_doxorubicin_trastuzumab_bnp_ac16 (in vitro BNP release) and Yu_2025_doxorubicin_trastuzumab_bnp_human (clinical BNP scale-up)."
  reference <- "Yu S, Mody H, Vaidya TR, Kagan L, Ait-Oudhia S. Mitigating Trastuzumab-Doxorubicin Cardiotoxicity With Multiscale Quantitative Systems Toxicology and PBPK-Toxicodynamic Predictive Modeling Framework. CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1625-1636. doi:10.1002/psp4.70087"
  vignette <- "Yu_2025_trastuzumab_doxorubicin_cardiotoxicity"

  units <- list(
    time = "h",
    dosing = "nM (doxorubicin concentration applied to the culture medium at t = 0)",
    concentration = "nM (medium doxorubicin concentration); caspase activity and cell viability are dimensionless ratios to their t = 0 values"
  )

  # The caspase / transit / viability states are the paper's own Figure 1
  # compartments (Cas9S, Cas91, Cas92, Cas31-Cas33, Tr1, Tr2); conc_dox
  # follows the in-vitro medium-concentration precedent of
  # Mody_2023_doxorubicin_dexrazoxane_jimt1.R from the same group.
  paper_specific_compartments <- c(
    "conc_dox",
    "cas9s",
    "cas91",
    "cas92",
    "cas31",
    "cas32",
    "cas33",
    "tr1",
    "tr2"
  )

  compartmentData <- list(
    conc_dox = list(analyte = "doxorubicin", units = "nM", specimen = "administration site", verified = TRUE),
    cas9s = list(
      analyte = "hypothetical caspase-9 stimulatory signal (TmAb-driven)",
      units = "ratio to t = 0",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas91 = list(
      analyte = "caspase-9 activation transit",
      units = "ratio to t = 0",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas92 = list(analyte = "active caspase-9", units = "ratio to t = 0", specimen = "not applicable", verified = TRUE),
    cas31 = list(
      analyte = "caspase-3 activation transit 1",
      units = "ratio to t = 0",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas32 = list(
      analyte = "caspase-3 activation transit 2",
      units = "ratio to t = 0",
      specimen = "not applicable",
      verified = TRUE
    ),
    cas33 = list(analyte = "active caspase-3", units = "ratio to t = 0", specimen = "not applicable", verified = TRUE),
    tr1 = list(analyte = "cell-death signal transit 1", units = "1/h", specimen = "not applicable", verified = TRUE),
    tr2 = list(analyte = "cell-death signal transit 2", units = "1/h", specimen = "not applicable", verified = TRUE),
    viability = list(
      analyte = "AC16 cell viability",
      units = "ratio to t = 0",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONC_TRASTUZUMAB_NM = list(
      description = "Static trastuzumab concentration applied to the AC16 culture medium",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Constant for the whole experiment (the paper models no TmAb loss in this layer). 0 for control and DOX-only wells. Experimental levels: 25,000 nM (25 uM) and 12,500 nM (12.5 uM, with 50 nM DOX) in the caspase experiment; 13,330 / 23,000 / 28,380 nM (13.3 / 23 / 28.4 uM) in the viability experiment (Supplement model code, REG blocks).",
      source_name = "CTmab"
    )
  )

  population <- list(
    species = "in vitro (AC16 human cardiomyocyte cell line)",
    n_subjects = NA_integer_,
    n_studies = 2L,
    disease_state = "Immortalised HER2-positive human cardiomyocytes (AC16, MilliporeSigma) used as an in vitro model of anthracycline / anti-HER2 cardiotoxicity.",
    dose_range = "Caspase experiment: control, DOX 100 nM, TmAb 25 uM, DOX 50 nM + TmAb 12.5 uM, sampled at 0, 3, 6, 12, 24, 48 and 72 h. Viability experiment: control, DOX 1.1 / 10.5 / 25 nM, TmAb 13.3 / 23 / 28.4 uM and all nine DOX x TmAb combinations, over 96 h.",
    cell_line = "AC16 (MilliporeSigma-EMD Millipore); HER2-positive immortalised human cardiomyocytes.",
    culture = "DMEM with 10% FBS, 1% sodium bicarbonate, 1% MEM non-essential amino acids and 1% penicillin/streptomycin; 37 C, humidified 5% CO2; seeded at 10^4 cells/well (100 uL) in 96-well plates for the CCK-8 viability assay.",
    notes = "Caspase-9 and -3 activity measured by MILLIPLEX / Luminex MAGPIX; viability by CCK-8 absorbance at 450 nm; at least triplicate wells per condition. Fitted in Monolix 2021R1 with inter-treatment-group (not inter-individual) random effects. Each output is normalised to its own t = 0 value, so relative viability grows ~29-fold over 96 h in untreated wells."
  )

  ini({
    # ---------------------------------------------------------------
    # DOX degradation in the culture medium (Eq 1)
    # ---------------------------------------------------------------
    lkdeg_dox <- fixed(log(0.022)); label("First-order degradation rate constant of DOX in medium kdeg (1/h)") # Table 1 'k deg' 0.022 (fix), from Vaidya 2021 (ref 27); Supplement code Kdeg = 0.022

    # ---------------------------------------------------------------
    # Protein-dynamics (QST) layer -- Table 1, cellular protein dynamic model
    # ---------------------------------------------------------------
    c1 <- 1.68; label("Power coefficient on DOX concentration for caspase-9 stimulation c1 (unitless)") # Table 1 'c1' 1.68 (RSE 3%)
    lsc1 <- log(0.046e-2); label("Slope for the DOX effect on caspase-9 activation Sc1 (nM^-c1)") # Table 1 'Sc1 (10-2 nM-1)' 0.046, i.e. 0.046 x 10^-2 (RSE 21%)
    lkcas9s <- log(0.093e-3); label("Turnover rate constant of the TmAb-driven caspase-9 stimulatory compartment kcas9s (1/h)") # Table 1 'k cas9s (10-3 h-1)' 0.093, i.e. 0.093 x 10^-3 (RSE 67%)
    lsc2 <- log(0.057); label("Slope for the TmAb effect on caspase-9 activation Sc2 (nM^-c2)") # Table 1 'Sc2 (nM-1)' 0.057 (RSE 104%)
    lkcas9 <- log(0.083); label("Caspase-9 turnover rate constant kcas9 (1/h)") # Table 1 'k cas9 (h-1)' 0.083 (RSE 15%)
    lkcas3 <- log(0.87); label("Caspase-3 turnover rate constant kcas3 (1/h)") # Table 1 'k cas3 (h-1)' 0.87 (RSE 9%)
    c2 <- 0.99; label("Power coefficient on TmAb concentration for caspase-9 stimulation c2 (unitless)") # Table 1 'c2' 0.99 (RSE 0.3%)
    c3 <- 4.48; label("Power coefficient of caspase-9 on caspase-3 activation c3 (unitless)") # Table 1 'c3' 4.48 (RSE 32%)

    # ---------------------------------------------------------------
    # Cardiomyocyte-death layer -- Table 1, cardiomyocytes-death model
    # ---------------------------------------------------------------
    lsmax_cas3 <- log(0.23); label("Maximal stimulation of cell death through caspase-3 activation SmaxCas3 (1/h)") # Table 1 'Smaxcas3' 0.23 (RSE 59%)
    c4 <- 1.5; label("Hill coefficient for caspase-3-driven cell death c4 (unitless)") # Table 1 'c4' 1.5 (RSE 34%)
    lsc50_cas3 <- log(2.56); label("Caspase-3 rise from baseline giving half-maximal death stimulation SC50Cas3 (ratio units)") # Table 1 'SC50cas3' 2.56 (RSE 66%)
    lktr_death <- log(0.023); label("Transit rate constant for the stimulation of AC16 cell death kTr (1/h)") # Table 1 'k tr (h-1)' 0.023 (RSE 39%)
    lkgrow <- log(0.035); label("First-order growth rate constant of AC16 cells kg (1/h)") # Table 1 'k g (h-1)' 0.035 (RSE 2%)
    lkd2 <- log(0.022); label("Slope for DOX inhibition of the cell growth rate kd2 (1/nM)") # Table 1 'kd2 (nM-1)' 0.022 (RSE 31%)

    # ---------------------------------------------------------------
    # Inter-treatment-group variability. Monolix reports omega as the SD
    # of eta; the paper's %CV column equals sqrt(exp(omega^2) - 1)
    # (0.049 -> 5%, 0.99 -> 129%, 0.69 -> 78%), confirming omega = SD.
    # ---------------------------------------------------------------
    etalsc1 ~ 0.049^2 # Table 1 'IIVSc1' omega 0.049 (5% CV)
    etalsc2 ~ 0.99^2 # Table 1 'IIVSc2' omega 0.99 (129% CV)
    etalkd2 ~ 0.69^2 # Table 1 'IIVkd2' omega 0.69 (78% CV)

    # ---------------------------------------------------------------
    # Residual error (Monolix proportional, y = f + b*f*e)
    # ---------------------------------------------------------------
    propSd_cas92 <- 0.18; label("Proportional residual error on relative caspase-9 activity (fraction)") # Table 1 protein model 'Ba' 0.18 (caspase-9)
    propSd_cas33 <- 0.43; label("Proportional residual error on relative caspase-3 activity (fraction)") # Table 1 protein model 'Bb' 0.43 (caspase-3)
    propSd_viability <- 0.36; label("Proportional residual error on relative AC16 cell viability (fraction)") # Table 1 death model 'Ba' 0.36 (cell viability)
  })

  model({
    kdeg_dox <- exp(lkdeg_dox)
    sc1 <- exp(lsc1 + etalsc1)
    kcas9s <- exp(lkcas9s)
    sc2 <- exp(lsc2 + etalsc2)
    kcas9 <- exp(lkcas9)
    kcas3 <- exp(lkcas3)
    smax_cas3 <- exp(lsmax_cas3)
    sc50_cas3 <- exp(lsc50_cas3)
    ktr_death <- exp(lktr_death)
    kgrow <- exp(lkgrow)
    kd2 <- exp(lkd2 + etalkd2)

    # Eq 1: the applied DOX concentration enters as a bolus "dose" (nM)
    # into conc_dox at t = 0, so the state IS the medium concentration.
    d/dt(conc_dox) <- -kdeg_dox * conc_dox

    # Eq 2: TmAb-driven upstream stimulatory compartment (constant TmAb)
    d/dt(cas9s) <- kcas9s * (1 + sc2 * CONC_TRASTUZUMAB_NM^c2) - kcas9s * cas9s

    # Eqs 3-4: caspase-9; DOX and TmAb act multiplicatively (Cas91 x Cas9S)
    d/dt(cas91) <- kcas9 * (1 + sc1 * conc_dox^c1) - kcas9 * cas91 * cas9s
    d/dt(cas92) <- kcas9 * cas91 * cas9s - kcas9 * cas92

    # Eqs 5-7: caspase-3 transit chain driven by caspase-9 to the power c3.
    # Printed Eq 5 multiplies the input by Cas9S; the deposited Monolix code
    # (Supplement) does not. The code form reproduces the published Figure 2
    # fits (TmAb-only caspase-3 peak ~4, not ~11), so it is used here.
    d/dt(cas31) <- kcas3 * cas92^c3 - kcas3 * cas31
    d/dt(cas32) <- kcas3 * cas31 - kcas3 * cas32
    d/dt(cas33) <- kcas3 * cas32 - kcas3 * cas33

    # Eq 8: caspase-3-driven death-rate signal (Hill on the rise from
    # baseline; the rise is floored at 0 so the fractional power stays real
    # if caspase-3 dips below baseline)
    dcas3 <- max(cas33 - 1, 0)
    kd1 <- smax_cas3 * dcas3^c4 / (dcas3^c4 + sc50_cas3^c4)

    # Eqs 9-10: death-signal transit chain
    d/dt(tr1) <- ktr_death * kd1 - ktr_death * tr1
    d/dt(tr2) <- ktr_death * tr1 - ktr_death * tr2

    # Eq 11: exponential growth, DOX-inhibited, minus caspase-3-driven death
    d/dt(viability) <- (1 - kd2 * conc_dox) * kgrow * viability - tr2 * viability

    cas9s(0) <- 1
    cas91(0) <- 1
    cas92(0) <- 1
    cas31(0) <- 1
    cas32(0) <- 1
    cas33(0) <- 1
    tr1(0) <- 0
    tr2(0) <- 0
    viability(0) <- 1

    cas92 ~ prop(propSd_cas92)
    cas33 ~ prop(propSd_cas33)
    viability ~ prop(propSd_viability)
  })
}
