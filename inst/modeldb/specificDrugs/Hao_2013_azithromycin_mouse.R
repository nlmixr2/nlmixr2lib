Hao_2013_azithromycin_mouse <- function() {
  description <- paste(
    "Preclinical (mouse, male CD1). Mechanism-based PK-PD model of",
    "azithromycin (AZI) inhibiting lipopolysaccharide (LPS)-induced",
    "depressive-like behaviour (Hao 2013). AZI (intragastric) and LPS",
    "(intraperitoneal) each follow one-compartment disposition with",
    "first-order absorption, with condition-specific PK: AZI PK differs in",
    "LPS-challenged mice (LIPOPOLYSACCHARIDE_CHALLENGE) and LPS PK differs in AZI-treated",
    "mice (CONMED_AZITHROMYCIN). Circulating LPS stimulates (Smax / SC50)",
    "and AZI inhibits (Imax / IC50) production of total proinflammatory",
    "cytokines (PCs), which drive plasma kynurenine (KYN), which in turn",
    "drives forced-swimming-test (FST) and tail-suspension-test (TST)",
    "immobility time through a series of turnover (transduction) models.",
    "The KYN and FST / TST production terms are baseline-normalised",
    "(kin * (1 + X / X0) with kin = X0 * kout / 2), a maintainers'",
    "reconstruction of Eqs 6-7, whose printed form kin * (1 + X) cannot",
    "hold the reported baselines. Naive-pooled ADAPT II fit: typical",
    "values only, no between-animal variability or residual error."
  )
  reference <- paste(
    "Hao K, Qi Q, Hao H, Wang G, Chen Y, Liang Y, Xie L. (2013). The",
    "pharmacokinetic-pharmacodynamic model of azithromycin for",
    "lipopolysaccharide-induced depressive-like behavior in mice.",
    "PLoS ONE 8(1):e54981. doi:10.1371/journal.pone.0054981.",
    "Model equations from Eqs 1-7 and the ADAPT II code of Appendix S1;",
    "parameter values from Tables 1-3."
  )
  vignette <- "Hao_2013_azithromycin_mouse"

  # LPS (the challenge agent) carries its own absorption and central pools,
  # named as the drug-disposition canonicals with an `_lps` qualifier. The
  # four PD states are the paper's named turnover pools: total
  # proinflammatory cytokines, plasma kynurenine, and the immobility times
  # of the two behavioural tests. None has a canonical analogue.
  paper_specific_compartments <- c(
    "depot_lps",
    "central_lps",
    "pcs",
    "kyn",
    "fst",
    "tst"
  )

  dosing <- c("depot", "depot_lps")

  units <- list(
    time = "h",
    dosing = "mg/kg",
    concentration = paste(
      "ug/mL for Cc (azithromycin); EU/mL for Cc_lps (LPS); pg/mL for pcs;",
      "ng/mL for kyn; s for fst and tst"
    )
  )

  covariateData <- list(
    LIPOPOLYSACCHARIDE_CHALLENGE = list(
      description = paste(
        "1 = mouse challenged with intraperitoneal lipopolysaccharide",
        "(LPS); 0 = not challenged."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal mouse; Table 1 'Normal Mice' AZI PK)",
      notes = paste(
        "Selects the azithromycin PK parameter set: 0 = Table 1 'Normal",
        "Mice', 1 = Table 1 'LPS-Challenging Mice' (ka 0.85 1/h, ke",
        "0.042 1/h, V/F 28.5 L/kg). Per-animal, time-fixed. It does not",
        "itself add LPS: the LPS dose is given as an event into depot_lps.",
        "Set it to 1 exactly when an LPS dose is in the event table."
      ),
      source_name = "group (normal vs LPS-challenging mice)"
    ),
    CONMED_AZITHROMYCIN = list(
      description = "1 = azithromycin given to the mouse; 0 = no azithromycin.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal mouse; Table 2 'Normal Mice' LPS PK)",
      notes = paste(
        "Selects the LPS PK parameter set: 0 = Table 2 'Normal Mice', 1 =",
        "Table 2 'Azithromycin-Treating Mice' (ka_LPS 1.63 1/h, ke_LPS",
        "0.088 1/h, V_LPS/F 15.5 L/kg). Per-animal, time-fixed; in the",
        "study AZI 100 mg/kg was given intragastrically 20 min before LPS.",
        "Set it to 1 exactly when an azithromycin dose is in the event table."
      ),
      source_name = "group (normal vs azithromycin-treating mice)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "azithromycin",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "azithromycin",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    depot_lps = list(
      analyte = "lipopolysaccharide (E. coli serotype O127:B8)",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_lps = list(
      analyte = "lipopolysaccharide (E. coli serotype O127:B8)",
      units = "mg/kg",
      specimen = "plasma",
      verified = TRUE
    ),
    pcs = list(
      analyte = "total proinflammatory cytokines (IL-1beta + IL-6 + IFN-gamma + TNF-alpha)",
      units = "pg/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    kyn = list(
      analyte = "kynurenine",
      units = "ng/mL",
      specimen = "plasma",
      verified = TRUE
    ),
    fst = list(
      analyte = "forced swimming test immobility time",
      units = "s",
      specimen = "not applicable",
      verified = TRUE
    ),
    tst = list(
      analyte = "tail suspension test immobility time",
      units = "s",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (CD1, male)",
    n_subjects = "8 mice per time point per group, sacrificed at each time point (naive pooled)",
    n_studies = 1L,
    age_range = "10-14 weeks",
    weight_range = "35-40 g",
    sex_female_pct = 0,
    disease_state = "LPS-induced depressive-like behaviour (sickness / inflammation model)",
    dose_range = paste(
      "Single intragastric azithromycin 100 mg/kg given 20 min before a",
      "single intraperitoneal LPS 0.8 mg/kg"
    ),
    regions = "China (China Pharmaceutical University, Nanjing)",
    notes = paste(
      "Groups: control, LPS alone, AZI + LPS, and AZI alone (AZI-alone PD",
      "not shown in the paper). One set of animals underwent the",
      "behavioural tests and a second set gave plasma for AZI, LPS, PCs",
      "and KYN (Methods, Study Design)."
    )
  )

  ini({
    # Azithromycin one-compartment oral PK (Eqs 1-2), Table 1.
    lka <- log(0.78)
    label("AZI absorption rate constant ka, normal mice (1/h)") # Table 1 ka Normal Mice = 0.78
    lkel <- log(0.040)
    label("AZI elimination rate constant ke, normal mice (1/h)") # Table 1 ke Normal Mice = 0.040
    lvc <- log(42.1)
    label("AZI apparent volume V_AZI/F_AZI, normal mice (L/kg)") # Table 1 V/F Normal Mice = 42.1
    lka_chall <- log(0.85)
    label("AZI absorption rate constant ka, LPS-challenged mice (1/h)") # Table 1 ka LPS-Challenging Mice = 0.85
    lkel_chall <- log(0.042)
    label("AZI elimination rate constant ke, LPS-challenged mice (1/h)") # Table 1 ke LPS-Challenging Mice = 0.042
    lvc_chall <- log(28.5)
    label("AZI apparent volume V_AZI/F_AZI, LPS-challenged mice (L/kg)") # Table 1 V/F LPS-Challenging Mice = 28.5

    # LPS one-compartment intraperitoneal PK (Eqs 3-4), Table 2.
    lka_lps <- log(1.51)
    label("LPS absorption rate constant ka_LPS, normal mice (1/h)") # Table 2 ka_LPS Normal Mice = 1.51
    lkel_lps <- log(0.082)
    label("LPS elimination rate constant ke_LPS, normal mice (1/h)") # Table 2 ke_LPS Normal Mice = 0.082
    lvc_lps <- log(16.3)
    label("LPS apparent volume V_LPS/F_LPS, normal mice (L/kg)") # Table 2 V_LPS/F_LPS Normal Mice = 16.3
    lka_lps_azi <- log(1.63)
    label("LPS absorption rate constant ka_LPS, AZI-treated mice (1/h)") # Table 2 ka_LPS Azithromycin-Treating Mice = 1.63
    lkel_lps_azi <- log(0.088)
    label("LPS elimination rate constant ke_LPS, AZI-treated mice (1/h)") # Table 2 ke_LPS Azithromycin-Treating Mice = 0.088
    lvc_lps_azi <- log(15.5)
    label("LPS apparent volume V_LPS/F_LPS, AZI-treated mice (L/kg)") # Table 2 V_LPS/F_LPS Azithromycin-Treating Mice = 15.5
    eu_per_ng_lps <- fixed(3.4)
    label("LPS endotoxin activity per mass (EU/ng)") # Methods (after Eq 4): 'per nanogram LPS equal to 3.4 Eu approximately'

    # Proinflammatory cytokines (Eq 5), Table 3.
    lrbase_pcs <- log(49.8)
    label("Baseline total proinflammatory cytokines PCs0 (pg/mL)") # Table 3 PCs,0 = 49.8
    lkout_pcs <- log(0.234)
    label("PCs first-order elimination rate constant kout_PCs (1/h)") # Table 3 kout_PCs = 0.234
    lsmax <- log(27.5)
    label("Maximal LPS stimulation of PCs production Smax (unitless)") # Table 3 Smax = 27.5
    lsc50 <- log(103.5)
    label("LPS concentration giving half-maximal stimulation SC50 (EU/mL)") # Table 3 SC50 = 103.5
    limax <- log(0.83)
    label("Maximal AZI inhibition of PCs production Imax (fraction)") # Table 3 Imax = 0.83
    lic50 <- log(6.3)
    label("AZI concentration giving half-maximal inhibition IC50 (ug/mL)") # Table 3 IC50 = 6.3

    # Kynurenine (Eq 6), Table 3.
    lrbase_kyn <- log(104.3)
    label("Baseline plasma kynurenine KYN0 (ng/mL)") # Table 3 KYN0 = 104.3
    lkout_kyn <- log(0.167)
    label("KYN first-order elimination rate constant kout_KYN (1/h)") # Table 3 kout_KYN = 0.167

    # Behavioural immobility times (Eq 7), Table 3.
    lrbase_fst <- log(51.3)
    label("Baseline forced swimming test immobility time FST0 (s)") # Table 3 FST0 = 51.3
    lkout_fst <- log(0.253)
    label("FST hypothetical elimination rate constant kout_FST (1/h)") # Table 3 kout_FST = 0.253
    lrbase_tst <- log(108.2)
    label("Baseline tail suspension test immobility time TST0 (s)") # Table 3 TST0 = 108.2
    lkout_tst <- log(0.143)
    label("TST hypothetical elimination rate constant kout_TST (1/h)") # Table 3 kout_TST = 0.143
  })

  model({
    # Condition-specific PK parameter sets (Tables 1 and 2).
    ka <- exp(lka * (1 - LIPOPOLYSACCHARIDE_CHALLENGE) + lka_chall * LIPOPOLYSACCHARIDE_CHALLENGE)
    kel <- exp(lkel * (1 - LIPOPOLYSACCHARIDE_CHALLENGE) + lkel_chall * LIPOPOLYSACCHARIDE_CHALLENGE)
    vc <- exp(lvc * (1 - LIPOPOLYSACCHARIDE_CHALLENGE) + lvc_chall * LIPOPOLYSACCHARIDE_CHALLENGE)
    ka_lps <- exp(
      lka_lps * (1 - CONMED_AZITHROMYCIN) + lka_lps_azi * CONMED_AZITHROMYCIN
    )
    kel_lps <- exp(
      lkel_lps * (1 - CONMED_AZITHROMYCIN) + lkel_lps_azi * CONMED_AZITHROMYCIN
    )
    vc_lps <- exp(
      lvc_lps * (1 - CONMED_AZITHROMYCIN) + lvc_lps_azi * CONMED_AZITHROMYCIN
    )

    smax <- exp(lsmax)
    sc50 <- exp(lsc50)
    imax <- exp(limax)
    ic50 <- exp(lic50)

    rbase_pcs <- exp(lrbase_pcs)
    kout_pcs <- exp(lkout_pcs)
    kin_pcs <- rbase_pcs * kout_pcs
    rbase_kyn <- exp(lrbase_kyn)
    kout_kyn <- exp(lkout_kyn)
    rbase_fst <- exp(lrbase_fst)
    kout_fst <- exp(lkout_fst)
    rbase_tst <- exp(lrbase_tst)
    kout_tst <- exp(lkout_tst)
    # Eqs 6-7 drive production with (1 + X / X0); halving kin keeps X0 a
    # true steady state (2 * kin = X0 * kout at X = X0). See the vignette
    # for why the printed (1 + X) form cannot be used.
    kin_kyn <- rbase_kyn * kout_kyn / 2
    kin_fst <- rbase_fst * kout_fst / 2
    kin_tst <- rbase_tst * kout_tst / 2

    # Concentrations: AZI mg/kg / (L/kg) = mg/L = ug/mL; LPS mg/kg converted
    # to EU/kg (1e6 ng/mg * 3.4 EU/ng) and L/kg to mL/kg (x 1000).
    Cc <- central / vc
    Cc_lps <- central_lps * 1e6 * eu_per_ng_lps / (vc_lps * 1000)

    stim_lps <- 1 + smax * Cc_lps / (Cc_lps + sc50)
    inh_azi <- 1 - imax * Cc / (Cc + ic50)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    d/dt(depot_lps) <- -ka_lps * depot_lps
    d/dt(central_lps) <- ka_lps * depot_lps - kel_lps * central_lps
    d/dt(pcs) <- kin_pcs * stim_lps * inh_azi - kout_pcs * pcs
    d/dt(kyn) <- kin_kyn * (1 + pcs / rbase_pcs) - kout_kyn * kyn
    d/dt(fst) <- kin_fst * (1 + kyn / rbase_kyn) - kout_fst * fst
    d/dt(tst) <- kin_tst * (1 + kyn / rbase_kyn) - kout_tst * tst

    pcs(0) <- rbase_pcs
    kyn(0) <- rbase_kyn
    fst(0) <- rbase_fst
    tst(0) <- rbase_tst
  })
}
