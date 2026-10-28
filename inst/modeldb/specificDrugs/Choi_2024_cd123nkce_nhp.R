Choi_2024_cd123nkce_nhp <- function() {
  description <- "Preclinical (nonhuman primate). Mechanistic synapse-driven PK/PD model of the trifunctional NKp46/CD16a x CD123 natural killer cell engager (CD123-NKCE) in healthy monkeys after IV and SC dosing. Two-compartment linear PK with first-order SC absorption, target-mediated disposition by sequential binding to NKp46 on circulating NK cells and CD123 on circulating CD123+ cells (forming the NKp46-drug-CD123 synapse) plus binding to an extra, turned-over CD123 pool in non-blood tissue whose baseline scales with blood CD123; internalisation of drug-CD123 complexes; an onset-time switch for anti-drug antibodies that adds a first-order elimination and stops further drug input. The synapse drives Hill-type CD123+ cell killing (with a moderator-mediated tolerance feedback on CD123+ cell regrowth) and a linear-slope increase of NK-cell redistribution out of blood (Choi 2024 Appendix S1 Eqs 9-21)."
  reference <- paste(
    "Choi SL, Valente D, Virone-Oddos A, Mauriac C.",
    "Developing a mechanistic translational PK/PD model for a trifunctional NK cell",
    "engager to predict the first-in-human dose for acute myeloid leukemia.",
    "Clin Transl Sci. 2024;17(1):e13689. doi:10.1111/cts.13689",
    sep = " "
  )
  vignette <- "Choi_2024_cd123nkce"
  # Time is hours. Doses are entered in ug/kg (the paper's dose unit); f()
  # converts them to nmol/kg with the molecular weight, the same conversion
  # the Monolix script applies to its dosing inputs (Table S2
  # 'p = 1/MW*1e3'). depot, central and peripheral1 therefore hold nmol/kg;
  # every binding state holds nmol/L; the cell states hold cells/uL.
  units <- list(time = "h", dosing = "ug/kg", concentration = "nmol/L")

  # The drug-receptor binary complexes on blood cells and the two circulating
  # cell pools have no canonical compartment name; they are named after the
  # drug-<antigen>-<site> family of drug_cd3_central and the tcell_central
  # cell-pool family.
  paper_specific_compartments <- c(
    "drug_nkp46_central",
    "drug_cd123_central",
    "nkcell_central",
    "cd123cell_central"
  )

  compartmentData <- list(
    depot = list(analyte = "CD123-NKCE", units = "nmol/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "CD123-NKCE", units = "nmol/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "CD123-NKCE", units = "nmol/kg", specimen = "tissue", verified = TRUE),
    drug_nkp46_central = list(
      analyte = "CD123-NKCE-NKp46 binary complex",
      units = "nmol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    drug_cd123_central = list(
      analyte = "CD123-NKCE-CD123 binary complex",
      units = "nmol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    trimer = list(
      analyte = "NKp46-CD123-NKCE-CD123 synapse",
      units = "nmol/L",
      specimen = "whole blood",
      verified = TRUE
    ),
    target = list(
      analyte = "CD123 expressed in non-blood tissue",
      units = "nmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    complex = list(
      analyte = "CD123-NKCE bound to non-blood CD123",
      units = "nmol/L",
      specimen = "tissue",
      verified = TRUE
    ),
    nkcell_central = list(analyte = "NKp46+ NK cells", units = "cells/uL", specimen = "whole blood", verified = TRUE),
    cd123cell_central = list(analyte = "CD123+ cells", units = "cells/uL", specimen = "whole blood", verified = TRUE),
    moderator1 = list(
      analyte = "tolerance moderator of CD123+ cell regrowth",
      units = "cells/uL",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "nonhuman primate (monkey)",
    n_subjects = 53,
    n_studies = 5,
    sex_female_pct = 41.5,
    disease_state = "Healthy",
    dose_range = "IV 1-h infusion 0.5-3000 ug/kg (single, weekly or twice weekly); SC 3.6-3000 ug/kg (single or twice weekly)",
    notes = paste(
      "Five studies (Table S1): study 1, 3000 ug/kg IV single dose (3 females); study 2, 0.5, 3",
      "and 3000 ug/kg IV and 1000 ug/kg SC single dose (9 males); study 3, 100 and 3000 ug/kg IV",
      "weekly x4 (4 males, 4 females); study 4, 3.6 ug/kg SC single dose (3 males); study 5,",
      "3000 ug/kg SC or 100 and 3000 ug/kg IV twice weekly x9 (15 males, 15 females). PK and PD",
      "from 23 animals were fitted; the 30 twice-weekly animals contributed PK only and their PD",
      "was used for external validation. Baseline NKp46+ NK cells median 247 [144-1110] and",
      "CD123+ cells median 28 [1.59-103] cells/uL (Table 1). Monolix 2020 (SAEM), BLQ data by M3."
    )
  )

  ini({
    # Linear PK, per kg body weight (Table 1, in vivo NHP column; estimated)
    lcl <- log(0.32)
    label("Nonspecific clearance per kg body weight (mL/h/kg)") # Table 1 NHP: CL = 0.32 (RSE 7.32%)
    lq <- log(1.35)
    label("Intercompartmental clearance per kg body weight (mL/h/kg)") # Table 1 NHP: Q = 1.35 (RSE 4.24%)
    lvc <- log(0.040)
    label("Central volume of distribution per kg body weight (L/kg)") # Table 1 NHP: V1 = 0.040 (RSE 4.23%)
    lvp <- log(0.020)
    label("Peripheral volume of distribution per kg body weight (L/kg)") # Table 1 NHP: V2 = 0.020 (RSE 14.3%)
    lka <- log(0.04)
    label("First-order SC absorption rate constant (1/h)") # Table 1 NHP: ka = 0.04 (RSE 9.12%)
    logitfdepot <- log(0.74 / (1 - 0.74))
    label("SC bioavailability (fraction)") # Table 1 NHP: F = 74.0% (RSE 6.87%)
    mw <- fixed(149000)
    label("Molecular weight of CD123-NKCE (g/mol)") # Table S2: MW = 1,49,000 g/mol

    # Anti-drug antibodies (Table 1 NHP; estimated)
    ltada <- log(333)
    label("Onset time of emergent anti-drug antibodies (h)") # Table 1 NHP: TADA = 333 (RSE 4.11%)
    lkada <- log(1.94)
    label("CD123-NKCE elimination rate constant due to anti-drug antibodies (1/h)") # Table 1 NHP: KADA = 1.94 (RSE 84.9%)

    # Target binding (Table 1 NHP)
    nkp46_receptors <- fixed(4600)
    label("NKp46 receptors per NK cell (receptors/cell)") # Table 1 NHP: RTdE = 4600 (internal data)
    cd123_receptors <- fixed(25000)
    label("CD123 receptors per CD123+ cell (receptors/cell)") # Table 1 NHP: RTdT = 25,000 (footnote d: ref 12)
    kd_nkp46 <- fixed(29)
    label("Equilibrium dissociation constant of CD123-NKCE for monkey NKp46 (nmol/L)") # Table 1 NHP: KD1 = 29 (internal data)
    kd_cd123 <- fixed(1.23)
    label("Equilibrium dissociation constant of CD123-NKCE for monkey CD123 (nmol/L)") # Table 1 NHP: KD2 = 1.23 (internal data)
    koff_nkp46 <- fixed(0.001)
    label("Dissociation rate constant of CD123-NKCE from NKp46 (1/s)") # Table 1 NHP: koff1 = 0.001 1/s (footnote o: ref 15)
    koff_cd123 <- fixed(0.001)
    label("Dissociation rate constant of CD123-NKCE from CD123 (1/s)") # Table 1 NHP: koff2 = 0.001 1/s (footnote o: ref 15)
    lkint <- log(0.15)
    label("Internalisation rate constant of drug-CD123 complexes (1/h)") # Table 1 NHP: kint = 0.15 (RSE 1.76%)

    # Extra CD123 expressed in non-blood tissue (Table 1 NHP; Appendix S1 Eqs 16-18)
    lksyn <- log(0.21)
    label("Zero-order production rate of non-blood CD123 (nmol/L/h)") # Table 1 NHP: kinTP = 0.21 (RSE 15.4%)
    ltarget_scale <- log(19.44)
    label("Scale (intercept) of non-blood CD123 versus blood CD123 at baseline (nmol/L)") # Table 1 NHP: alpha = 19.44 (RSE 16.3%)
    ltarget_exp <- log(0.21)
    label("Power exponent of non-blood CD123 versus blood CD123 at baseline (unitless)") # Table 1 NHP: beta = 0.21 (RSE 0.842%)

    # Synapse-driven killing (Table 1 NHP; carried from the in vitro model)
    lkmax <- fixed(log(0.332))
    label("Maximal CD123+ cell killing rate constant of the synapse (1/h)") # Table 1 NHP: kmax = 0.332 (from in vitro model, footnote g)
    lec50 <- fixed(log(1.46e-9))
    label("Synapse concentration at half-maximal killing rate (nmol/L)") # Table 1 NHP: EC50 = 1.46E-9 (from in vitro model, footnote g)
    lhill <- fixed(log(1.04))
    label("Hill coefficient of synapse-driven killing (unitless)") # Table 1 NHP: gamma = 1.04 (from in vitro model, footnote g)

    # CD123+ cell growth with tolerance feedback (Table 1 NHP; Appendix S1 Eqs 20-21)
    lkg <- log(0.31)
    label("First-order growth rate constant of CD123+ cells (1/h)") # Table 1 NHP: kg = 0.31 (RSE 50%)
    lktol <- log(0.0032)
    label("Rate constant of the tolerance moderator of CD123+ cell regrowth (1/h)") # Table 1 NHP: ktol = 0.0032 (RSE 55.9%)
    ldelta_tol <- log(0.18)
    label("Power coefficient of the tolerance feedback on CD123+ cell growth (unitless)") # Table 1 NHP: delta = 0.18 (RSE 32.4%)

    # NK-cell redistribution (Table 1 NHP; Appendix S1 Eq 19)
    lkout <- fixed(log(0.2))
    label("NKp46+ NK cell turnover rate constant (1/h)") # Table 1 NHP: kout = 0.2 (footnote j: calibrated)
    lslope_nk <- log(5.72)
    label("Slope of synapse-driven NK-cell redistribution (per 1e9 x synapse in nmol/L)") # Table 1 NHP: SlpNK = 5.72 (RSE 0.951%)

    # Baseline cell counts. The paper mapped each animal's observed baselines as
    # regressors; the cohort medians are carried here and can be overridden per
    # subject through rxSolve(params = ).
    rbase_nk <- fixed(247)
    label("Baseline NKp46+ NK cell count in blood (cells/uL)") # Table 1 NHP: EF0 median = 247 [144-1110] (regressor)
    rbase_cd123 <- fixed(28)
    label("Baseline CD123+ cell count in blood (cells/uL)") # Table 1 NHP: TC0 median = 28 [1.59-103] (regressor)

    # IIV (Table 1 NHP, IIV column: SD of the random effect, footnote a;
    # variance = SD^2; logit-normal for F, footnote p; no correlations reported)
    etalcl ~ 0.04 # Table 1 NHP: IIV CL = 0.20 (RSE 31.2%)
    etalvc ~ 0.0169 # Table 1 NHP: IIV V1 = 0.130 (RSE 31.3%)
    etalvp ~ 0.0676 # Table 1 NHP: IIV V2 = 0.26 (RSE 57.8%)
    etalogitfdepot ~ 0.0784 # Table 1 NHP: IIV F = 0.28 (RSE 55.2%), logit scale
    etaltada ~ 0.04 # Table 1 NHP: IIV TADA = 0.20 (RSE 15.1%)
    etalkada ~ 17.7241 # Table 1 NHP: IIV KADA = 4.21 (RSE 17.7%)
    etalksyn ~ 0.3249 # Table 1 NHP: IIV kinTP = 0.57 (RSE 22.6%)
    etaltarget_scale ~ 0.3721 # Table 1 NHP: IIV alpha = 0.61 (RSE 19.8%)
    etalkg ~ 0.6889 # Table 1 NHP: IIV kg = 0.83 (RSE 26.8%)
    etalktol ~ 2.6569 # Table 1 NHP: IIV ktol = 1.63 (RSE 43.6%)
    etaldelta_tol ~ 0.2809 # Table 1 NHP: IIV delta = 0.53 (RSE 35.5%)
  })

  model({
    # Avogadro's number per nmol (Table S2: NA = 6.022e14 #/nmol)
    navo <- 6.022e14

    # PK; CL and Q in mL/h/kg -> L/h/kg (Table S2 'k01 = (CL/1000)/V')
    cl <- exp(lcl + etalcl) / 1000
    q <- exp(lq) / 1000
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka)
    fdepot <- expit(logitfdepot + etalogitfdepot)
    tada <- exp(ltada + etaltada)
    kada <- exp(lkada + etalkada)

    kint <- exp(lkint)
    ksyn <- exp(lksyn + etalksyn)
    target_scale <- exp(ltarget_scale + etaltarget_scale)
    target_exp <- exp(ltarget_exp)
    kmax <- exp(lkmax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    kg <- exp(lkg + etalkg)
    ktol <- exp(lktol + etalktol)
    delta_tol <- exp(ldelta_tol + etaldelta_tol)
    kout <- exp(lkout)
    slope_nk <- exp(lslope_nk)

    # koff in 1/s -> 1/h; kon = koff / KD
    koff1 <- koff_nkp46 * 3600
    koff2 <- koff_cd123 * 3600
    kon1 <- koff1 / kd_nkp46
    kon2 <- koff2 / kd_cd123

    # Anti-drug antibody switch (Appendix S1 Eq 14): after TADA, KADA adds a
    # first-order elimination and the IV and SC inputs to central stop.
    ada_on <- 0
    if (t > tada) ada_on <- 1

    # Baseline blood CD123 (nmol/L) and non-blood CD123 pool (Appendix S1 Eq 18)
    blood_cd123_0 <- rbase_cd123 * cd123_receptors * 1e6 / navo
    target0 <- target_scale * blood_cd123_0^target_exp
    kdeg_target <- ksyn / target0

    # Free drug and free receptors on circulating cells (nmol/L)
    cfree <- central / vc
    efree <- nkcell_central * nkp46_receptors * 1e6 / navo - drug_nkp46_central - trimer
    tfree <- cd123cell_central * cd123_receptors * 1e6 / navo - drug_cd123_central - trimer

    # Binding flux on free drug (nmol/L/h; Appendix S1 Eq 14 binding terms)
    bind <- -kon1 * cfree * efree + koff1 * drug_nkp46_central -
      kon2 * cfree * tfree + koff2 * drug_cd123_central -
      kon2 * cfree * target + koff2 * complex

    # PK (Appendix S1 Eqs 11-14), amounts in nmol/kg
    d/dt(depot) <- -ka * depot
    d/dt(central) <- (1 - ada_on) * ka * depot - cl * cfree - q * cfree + q * peripheral1 / vp +
      bind * vc - ada_on * kada * central
    d/dt(peripheral1) <- q * cfree - q * peripheral1 / vp

    # Binding on circulating cells (Appendix S1 Eqs 2, 4 and 15)
    d/dt(drug_nkp46_central) <- kon1 * cfree * efree - koff1 * drug_nkp46_central -
      kon2 * tfree * drug_nkp46_central + koff2 * trimer
    d/dt(drug_cd123_central) <- kon2 * tfree * cfree - koff2 * drug_cd123_central -
      kon1 * drug_cd123_central * efree + koff1 * trimer - kint * drug_cd123_central
    d/dt(trimer) <- kon2 * tfree * drug_nkp46_central - koff2 * trimer +
      kon1 * drug_cd123_central * efree - koff1 * trimer

    # Non-blood CD123 turnover and binding (Appendix S1 Eqs 16-17)
    d/dt(target) <- ksyn - kdeg_target * target - kon2 * cfree * target + koff2 * complex
    d/dt(complex) <- kon2 * cfree * target - koff2 * complex - kint * complex
    target(0) <- target0

    # Synapse floored at 0 as in Table S2 'TBE2 = max(0, TBE)'
    synapse <- max(0, trimer)
    kkill <- kmax * synapse^hill / (ec50^hill + synapse^hill)

    # NK-cell redistribution (Appendix S1 Eq 19)
    d/dt(nkcell_central) <- kout * rbase_nk - kout * (1 + slope_nk * 1e9 * synapse) * nkcell_central
    nkcell_central(0) <- rbase_nk

    # CD123+ cells with moderator-mediated tolerance (Appendix S1 Eqs 20-21)
    d/dt(cd123cell_central) <- kg * (1 - (moderator1 / rbase_cd123)^delta_tol) * cd123cell_central -
      kkill * cd123cell_central
    d/dt(moderator1) <- ktol * (cd123cell_central - moderator1)
    cd123cell_central(0) <- rbase_cd123
    moderator1(0) <- rbase_cd123

    # Dose unit conversion ug/kg -> nmol/kg (Table S2 'p = F/MW*1e3'); no
    # further input once anti-drug antibodies emerge
    f(depot) <- fdepot * 1e3 / mw
    f(central) <- (1 - ada_on) * 1e3 / mw

    Cc <- cfree
    cd123_pctcfb <- (cd123cell_central - rbase_cd123) / rbase_cd123 * 100
  })
}
