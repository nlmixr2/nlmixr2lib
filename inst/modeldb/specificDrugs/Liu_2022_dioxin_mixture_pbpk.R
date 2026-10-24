Liu_2022_dioxin_mixture_pbpk <- function() {
  description <- paste(
    "PBPK/QSP (whole-body, lifetime human; original implementation in GNU",
    "Octave). Binary mixture of 2,3,7,8-TCDD and 2,3,4,7,8-PeCDF, a",
    "two-congener instance of the dioxin-like compound (DLC) mixture PBPK",
    "framework of Liu et al. 2022 (Toxics). It demonstrates the paper's",
    "central result (its Figure 3A): two congeners sharing the hepatic aryl",
    "hydrocarbon receptor (AHR) and CYP1A2 cross-induce CYP1A2, which",
    "accelerates sequestration and metabolism so that each congener's",
    "extrahepatic tissue burden in the mixture is lower than in single",
    "exposure, while the liver burden may rise. Each congener has its own",
    "GI tract, fat/rest-of-body/liver blood and diffusion-limited tissue",
    "subcompartments, and its own AHR and CYP1A2 complexes; the free CYP1A2",
    "protein, its mRNA and total AHR are single shared pools that couple the",
    "two congeners. Blood is quasi-steady-state (algebraic). Body weight and",
    "tissue volume/flow fractions are polynomial functions of age, so the",
    "model is a forward lifetime simulation of constant daily dietary intake",
    "with no dose records. Deterministic; no IIV and no residual error. The",
    "framework extends to an arbitrary number of congeners; this file fixes",
    "two. Parameter values are the authors' published Octave code",
    "(https://github.com/pulsatility/2022-DLC-Mixture-PBPK), which the paper",
    "Section 2.6 designates as the model source."
  )
  reference <- paste(
    "Liu R, Zacharewski TR, Conolly RB, Zhang Q. A Physiologically Based",
    "Pharmacokinetic (PBPK) Modeling Framework for Mixtures of Dioxin-like",
    "Compounds. Toxics. 2022;10(11):700. doi:10.3390/toxics10110700.",
    "Physiological and biochemical parameters from the authors' Octave code",
    "DLC_human_cmd.m / DLC_human_ode.m; chemical-specific parameters from",
    "DLC_human_parameters.xlsx (rows '2,3,7,8-TCDD' and '2,3,4,7,8-PeCDF';",
    "paper Tables S1-S3; https://github.com/pulsatility/2022-DLC-Mixture-PBPK)."
  )
  vignette <- "Liu_2022_dioxin_mixture_pbpk"

  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nM",
    amount = "nmol",
    weight = "kg"
  )

  paper_specific_compartments <- c(
    "gi_tcdd",
    "fat_blood_tcdd",
    "fat_tcdd",
    "rest_blood_tcdd",
    "rest_tcdd",
    "liver_blood_tcdd",
    "liver_tcdd",
    "dioxin_ahr_tcdd",
    "dioxin_cyp1a2_tcdd",
    "gi_pecdf",
    "fat_blood_pecdf",
    "fat_pecdf",
    "rest_blood_pecdf",
    "rest_pecdf",
    "liver_blood_pecdf",
    "liver_pecdf",
    "dioxin_ahr_pecdf",
    "dioxin_cyp1a2_pecdf",
    "cyp1a2_mrna",
    "cyp1a2"
  )

  compartmentData <- list(
    gi_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "administration site", verified = TRUE),
    fat_blood_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "whole blood", verified = TRUE),
    fat_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "tissue", verified = TRUE),
    rest_blood_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "whole blood", verified = TRUE),
    rest_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "tissue", verified = TRUE),
    liver_blood_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver_tcdd = list(analyte = "TCDD", units = "nmol", specimen = "tissue", verified = TRUE),
    dioxin_ahr_tcdd = list(analyte = "TCDD-AHR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    dioxin_cyp1a2_tcdd = list(analyte = "TCDD-CYP1A2 complex", units = "nmol", specimen = "tissue", verified = TRUE),
    gi_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "administration site", verified = TRUE),
    fat_blood_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "whole blood", verified = TRUE),
    fat_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "tissue", verified = TRUE),
    rest_blood_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "whole blood", verified = TRUE),
    rest_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "tissue", verified = TRUE),
    liver_blood_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver_pecdf = list(analyte = "2,3,4,7,8-PeCDF", units = "nmol", specimen = "tissue", verified = TRUE),
    dioxin_ahr_pecdf = list(analyte = "PeCDF-AHR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    dioxin_cyp1a2_pecdf = list(analyte = "PeCDF-CYP1A2 complex", units = "nmol", specimen = "tissue", verified = TRUE),
    cyp1a2_mrna = list(analyte = "CYP1A2 mRNA", units = "nM", specimen = "tissue", verified = TRUE),
    cyp1a2 = list(analyte = "CYP1A2 protein", units = "nM", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    age_range = "0-70 years (lifetime forward simulation)",
    weight_range = "age-dependent (birth 3.65 kg to adult ~64 kg, female growth curve)",
    sex_female_pct = 100,
    disease_state = paste(
      "Healthy. No clinical study underlies the model; it is a forward",
      "simulation of lifetime dietary intake of a TCDD + 2,3,4,7,8-PeCDF",
      "mixture. Results in the paper are presented for women unless stated."
    ),
    dose_range = paste(
      "Constant daily dietary intake per congener. The paper runs an EPA RfD",
      "(0.0007 ng/kg bw/day) and PoD (0.02 ng/kg bw/day) scenario; both",
      "congener intakes default to the PoD dose."
    ),
    regions = "exposure scenario (EPA RfD/PoD); dietary-intake composition from a 1994 Dutch survey",
    notes = paste(
      "Two-congener instance of the DLC mixture framework. 2,3,4,7,8-PeCDF",
      "is the representative second congener of the paper's Figure 3A. The",
      "framework scales to any number of congeners; the full 11-congener",
      "parameter set is given in the vignette."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Exposure, per congener (ng/kg bw/day).
    # -------------------------------------------------------------------
    MSTOT_tcdd <- fixed(0.02)
    label("Daily oral dietary intake of TCDD (ng/kg bw/day)") # DLC_human_cmd.m PoD = 0.02
    MSTOT_pecdf <- fixed(0.02)
    label("Daily oral dietary intake of 2,3,4,7,8-PeCDF (ng/kg bw/day)") # DLC_human_cmd.m PoD = 0.02

    # -------------------------------------------------------------------
    # Chemical-specific parameters, TCDD (xlsx row '2,3,7,8-TCDD').
    # -------------------------------------------------------------------
    PL_tcdd <- fixed(6.0)
    label("TCDD liver:blood partition coefficient (unitless)") # xlsx 'PL' = 6
    PF_tcdd <- fixed(100.0)
    label("TCDD fat:blood partition coefficient (unitless)") # xlsx 'PF' = 100
    PRB_tcdd <- fixed(1.5)
    label("TCDD rest-of-body:blood partition coefficient (unitless)") # xlsx 'PRB' = 1.5
    KDAHR_tcdd <- fixed(0.1)
    label("TCDD-AHR dissociation constant Kd,AHR (nM)") # xlsx 'KDAHR' = 0.1
    kelim_tcdd <- fixed(0.0011)
    label("TCDD CYP1A2-mediated hepatic elimination rate constant (1/h)") # xlsx 'kelim' = 0.0011
    KABS_tcdd <- fixed(0.06)
    label("TCDD GI-tract absorption rate constant (1/h)") # xlsx 'KABS' = 0.06
    MW_tcdd <- fixed(321.976)
    label("TCDD molecular weight (g/mol)") # xlsx 'MW' = 321.976

    # -------------------------------------------------------------------
    # Chemical-specific parameters, 2,3,4,7,8-PeCDF (xlsx row 1).
    # -------------------------------------------------------------------
    PL_pecdf <- fixed(27.551020)
    label("PeCDF liver:blood partition coefficient (unitless)") # xlsx 'PL' = 27.55102
    PF_pecdf <- fixed(136.03239)
    label("PeCDF fat:blood partition coefficient (unitless)") # xlsx 'PF' = 136.03239
    PRB_pecdf <- fixed(3.1666667)
    label("PeCDF rest-of-body:blood partition coefficient (unitless)") # xlsx 'PRB' = 3.1666667
    KDAHR_pecdf <- fixed(0.1492537)
    label("PeCDF-AHR dissociation constant Kd,AHR (nM)") # xlsx 'KDAHR' = 0.1492537
    kelim_pecdf <- fixed(0.0006888889)
    label("PeCDF CYP1A2-mediated hepatic elimination rate constant (1/h)") # xlsx 'kelim' = 0.0006888889
    KABS_pecdf <- fixed(0.06061856)
    label("PeCDF GI-tract absorption rate constant (1/h)") # xlsx 'KABS' = 0.06061856
    MW_pecdf <- fixed(340.422)
    label("PeCDF molecular weight (g/mol)") # xlsx 'MW' = 340.422

    # -------------------------------------------------------------------
    # Shared chemical parameters (equal across all 11 congeners in the xlsx).
    # -------------------------------------------------------------------
    PALC <- fixed(0.35)
    label("Liver permeability-area coefficient, fraction of liver blood flow (unitless)") # xlsx 'PALC' = 0.35
    PAFC <- fixed(0.12)
    label("Fat permeability-area coefficient, fraction of fat blood flow (unitless)") # xlsx 'PAFC' = 0.12
    PARBC <- fixed(0.03)
    label("Rest-of-body permeability-area coefficient, fraction of RB blood flow (unitless)") # xlsx 'PARBC' = 0.03
    Fu <- fixed(1.0)
    label("Fraction of congener unbound in tissue blood (unitless)") # xlsx 'Fu' = 1
    kbAHR <- fixed(6.0)
    label("Congener-AHR dissociation rate constant kb,AHR (1/h)") # paper Section 2.3; xlsx 'kbAHR' = 6
    KD1A2 <- fixed(40.0)
    label("Congener-CYP1A2 dissociation constant Kd,1A2 (nM)") # xlsx 'KD1A2' = 40
    kb1A2 <- fixed(6.0)
    label("Congener-CYP1A2 dissociation rate constant kb,1A2 (1/h)") # xlsx 'kb1A2' = 6
    KST <- fixed(0.01)
    label("GI-tract transit (loss) rate constant (1/h)") # xlsx 'KST' = 0.01
    a <- fixed(0.7)
    label("Fraction of absorbed congener entering blood via lymph (unitless)") # xlsx 'a' = 0.7
    CLURI <- fixed(4.17e-08)
    label("Slow first-order (urinary/non-CYP1A2) blood clearance (L/h)") # xlsx 'CLURI' = 4.17e-8

    # -------------------------------------------------------------------
    # Shared biochemical parameters, hepatic AHR / CYP1A2 (DLC_human_cmd.m).
    # -------------------------------------------------------------------
    AHRtot <- fixed(0.35)
    label("Total hepatic AHR concentration (nM)") # g.AHRtot = 0.35
    CYP1A2_1BASAL <- fixed(1600.0)
    label("Basal (and initial) free CYP1A2 protein concentration (nM)") # g.CYP1A2_1BASAL = 1600
    CYP1A2_mRNA_1BASAL <- fixed(160.0)
    label("Basal (and initial) CYP1A2 mRNA concentration (nM)") # g.CYP1A2_mRNA_1BASAL = 160
    CYP1A2_1EC50 <- fixed(130.0)
    label("Dioxin-AHR concentration for half-maximal CYP1A2 induction EC50 (nM)") # g.CYP1A2_1EC50 = 130
    CYP1A2_1EMAX <- fixed(9300.0)
    label("Maximal fold induction of CYP1A2 transcription Emax (unitless)") # g.CYP1A2_1EMAX = 9300
    hill_n <- fixed(0.6)
    label("Hill coefficient for CYP1A2 transcriptional induction (unitless)") # DLC_human_ode.m exponent 0.6
    ktranscription_1A2 <- fixed(16.0)
    label("CYP1A2 mRNA transcription rate constant (nM/h)") # g.ktranscription_1A2 = 16
    kdegCYP1A2_mRNA <- fixed(0.1)
    label("CYP1A2 mRNA degradation rate constant (1/h)") # g.kdegCYP1A2_mRNA = 0.1
    ktranslation_1A2 <- fixed(1.25)
    label("CYP1A2 protein translation rate constant (1/h)") # g.ktranslation_1A2 = 1.25
    kdegCYP1A2 <- fixed(0.125)
    label("CYP1A2 protein degradation rate constant (1/h)") # g.kdegCYP1A2 = 0.125
  })

  model({
    # =================================================================
    # Fixed physiology (DLC_human_cmd.m; Emond human, Wang et al. 1997).
    # =================================================================
    QCC <- 15.36
    QLC <- 0.26
    QFC <- 0.05
    QRBC <- 1 - QFC - QLC
    VLBC <- 0.266
    VFBC <- 0.050
    VRBBC <- 0.030
    ORAL_DURATION <- 24

    # =================================================================
    # Time-varying anatomy (female; DLC_human_ode.m).
    # =================================================================
    age <- time / 24 / 365
    BW <- 0.0006 * age^3 - 0.0912 * age^2 + 4.3200 * age + 3.6520
    BWg <- BW * 1000
    VLC <- 3.59e-2 - 4.76e-7 * BWg + 8.50e-12 * BWg^2 - 5.45e-17 * BWg^3
    VFC <- -6.36e-20 * BWg^4 + 1.12e-14 * BWg^3 - 5.8e-10 * BWg^2 + 1.2e-5 * BWg + 5.91e-2
    VRBC <- (0.91 - (VLBC * VLC + VFBC * VFC + VLC + VFC)) / (1 + VRBBC)

    VL <- VLC * BW
    VF <- VFC * BW
    VRB <- VRBC * BW
    VLB <- VLBC * VL
    VFB <- VFBC * VF
    VRBB <- VRBBC * VRB

    QC <- QCC * BW^0.75
    QF <- QFC * QC
    QL <- QLC * QC
    QRB <- QRBC * QC

    PAF <- PAFC * QF
    PARB <- PARBC * QRB
    PAL <- PALC * QL

    # Association rate constants (1/nM/h). kb,AHR, kb,1A2, Kd,1A2 shared, so
    # kf1A2 is common; kfAHR differs through each congener's Kd,AHR.
    kfAHR_tcdd <- kbAHR / KDAHR_tcdd
    kfAHR_pecdf <- kbAHR / KDAHR_pecdf
    kf1A2 <- kb1A2 / KD1A2
    b <- 1 - a

    ODR_tcdd <- (MSTOT_tcdd / MW_tcdd) * BW / ORAL_DURATION
    ODR_pecdf <- (MSTOT_pecdf / MW_pecdf) * BW / ORAL_DURATION

    # =================================================================
    # Quasi-steady-state blood (algebraic), per congener.
    # =================================================================
    CA_tcdd <- (QF * (fat_blood_tcdd / VFB) + QRB * (rest_blood_tcdd / VRBB) +
      QL * (liver_blood_tcdd / VLB) + KABS_tcdd * gi_tcdd * a) / (QC + CLURI)
    CA_pecdf <- (QF * (fat_blood_pecdf / VFB) + QRB * (rest_blood_pecdf / VRBB) +
      QL * (liver_blood_pecdf / VLB) + KABS_pecdf * gi_pecdf * a) / (QC + CLURI)

    # Liver-tissue free concentrations and complex concentrations (nM).
    CFLL_tcdd <- liver_tcdd / VL / PL_tcdd * Fu
    CFLL_pecdf <- liver_pecdf / VL / PL_pecdf * Fu
    Dioxin_AHR_tcdd <- dioxin_ahr_tcdd / VL
    Dioxin_AHR_pecdf <- dioxin_ahr_pecdf / VL
    Dioxin_CYP1A2_tcdd <- dioxin_cyp1a2_tcdd / VL
    Dioxin_CYP1A2_pecdf <- dioxin_cyp1a2_pecdf / VL

    # Shared pools: total AHR-bound and total CYP1A2-bound congener (nM).
    sum_Dioxin_AHR <- Dioxin_AHR_tcdd + Dioxin_AHR_pecdf
    sum_Dioxin_CYP1A2 <- Dioxin_CYP1A2_tcdd + Dioxin_CYP1A2_pecdf

    # =================================================================
    # ODEs (amounts in nmol). DLC_human_ode.m, M = 2.
    # =================================================================
    d/dt(gi_tcdd) <- -(KST + KABS_tcdd) * gi_tcdd + ODR_tcdd
    d/dt(gi_pecdf) <- -(KST + KABS_pecdf) * gi_pecdf + ODR_pecdf

    d/dt(fat_blood_tcdd) <- QF * (CA_tcdd - fat_blood_tcdd / VFB) -
      PAF * (fat_blood_tcdd / VFB * Fu - (fat_tcdd / VF) / PF_tcdd * Fu)
    d/dt(fat_tcdd) <- PAF * (fat_blood_tcdd / VFB * Fu - (fat_tcdd / VF) / PF_tcdd * Fu)
    d/dt(fat_blood_pecdf) <- QF * (CA_pecdf - fat_blood_pecdf / VFB) -
      PAF * (fat_blood_pecdf / VFB * Fu - (fat_pecdf / VF) / PF_pecdf * Fu)
    d/dt(fat_pecdf) <- PAF * (fat_blood_pecdf / VFB * Fu - (fat_pecdf / VF) / PF_pecdf * Fu)

    d/dt(rest_blood_tcdd) <- QRB * (CA_tcdd - rest_blood_tcdd / VRBB) -
      PARB * (rest_blood_tcdd / VRBB * Fu - rest_tcdd / VRB / PRB_tcdd * Fu)
    d/dt(rest_tcdd) <- PARB * (rest_blood_tcdd / VRBB * Fu - rest_tcdd / VRB / PRB_tcdd * Fu)
    d/dt(rest_blood_pecdf) <- QRB * (CA_pecdf - rest_blood_pecdf / VRBB) -
      PARB * (rest_blood_pecdf / VRBB * Fu - rest_pecdf / VRB / PRB_pecdf * Fu)
    d/dt(rest_pecdf) <- PARB * (rest_blood_pecdf / VRBB * Fu - rest_pecdf / VRB / PRB_pecdf * Fu)

    d/dt(liver_blood_tcdd) <- QL * (CA_tcdd - liver_blood_tcdd / VLB) -
      PAL * (liver_blood_tcdd / VLB * Fu - CFLL_tcdd) + KABS_tcdd * gi_tcdd * b
    d/dt(liver_blood_pecdf) <- QL * (CA_pecdf - liver_blood_pecdf / VLB) -
      PAL * (liver_blood_pecdf / VLB * Fu - CFLL_pecdf) + KABS_pecdf * gi_pecdf * b

    d/dt(liver_tcdd) <- PAL * (liver_blood_tcdd / VLB * Fu - CFLL_tcdd) +
      (-kelim_tcdd * CFLL_tcdd * (cyp1a2 + sum_Dioxin_CYP1A2 - CYP1A2_1BASAL) / CYP1A2_1BASAL -
        kfAHR_tcdd * CFLL_tcdd * (AHRtot - sum_Dioxin_AHR) + kbAHR * Dioxin_AHR_tcdd -
        kf1A2 * CFLL_tcdd * cyp1a2 + kb1A2 * Dioxin_CYP1A2_tcdd) * VL
    d/dt(liver_pecdf) <- PAL * (liver_blood_pecdf / VLB * Fu - CFLL_pecdf) +
      (-kelim_pecdf * CFLL_pecdf * (cyp1a2 + sum_Dioxin_CYP1A2 - CYP1A2_1BASAL) / CYP1A2_1BASAL -
        kfAHR_pecdf * CFLL_pecdf * (AHRtot - sum_Dioxin_AHR) + kbAHR * Dioxin_AHR_pecdf -
        kf1A2 * CFLL_pecdf * cyp1a2 + kb1A2 * Dioxin_CYP1A2_pecdf) * VL

    d/dt(dioxin_ahr_tcdd) <- (kfAHR_tcdd * CFLL_tcdd * (AHRtot - sum_Dioxin_AHR) - kbAHR * Dioxin_AHR_tcdd) * VL
    d/dt(dioxin_ahr_pecdf) <- (kfAHR_pecdf * CFLL_pecdf * (AHRtot - sum_Dioxin_AHR) - kbAHR * Dioxin_AHR_pecdf) * VL
    d/dt(dioxin_cyp1a2_tcdd) <- (kf1A2 * CFLL_tcdd * cyp1a2 - kb1A2 * Dioxin_CYP1A2_tcdd) * VL
    d/dt(dioxin_cyp1a2_pecdf) <- (kf1A2 * CFLL_pecdf * cyp1a2 - kb1A2 * Dioxin_CYP1A2_pecdf) * VL

    d/dt(cyp1a2_mrna) <- ktranscription_1A2 *
      (1 + CYP1A2_1EMAX * sum_Dioxin_AHR^hill_n / (CYP1A2_1EC50^hill_n + sum_Dioxin_AHR^hill_n)) -
      kdegCYP1A2_mRNA * cyp1a2_mrna
    d/dt(cyp1a2) <- ktranslation_1A2 * cyp1a2_mrna - kdegCYP1A2 * cyp1a2 -
      (kf1A2 * CFLL_tcdd * cyp1a2 + kf1A2 * CFLL_pecdf * cyp1a2) +
      (kb1A2 * Dioxin_CYP1A2_tcdd + kb1A2 * Dioxin_CYP1A2_pecdf)

    cyp1a2(0) <- CYP1A2_1BASAL
    cyp1a2_mrna(0) <- CYP1A2_mRNA_1BASAL

    # =================================================================
    # Observations (concentrations, nM).
    # =================================================================
    Cc_tcdd <- CA_tcdd
    Cc_pecdf <- CA_pecdf
    Cfat_tcdd <- fat_tcdd / VF
    Cfat_pecdf <- fat_pecdf / VF
    Crest_tcdd <- rest_tcdd / VRB
    Crest_pecdf <- rest_pecdf / VRB
    Cliver_total_tcdd <- (liver_tcdd + dioxin_ahr_tcdd + dioxin_cyp1a2_tcdd) / VL
    Cliver_total_pecdf <- (liver_pecdf + dioxin_ahr_pecdf + dioxin_cyp1a2_pecdf) / VL
    Ccyp1a2 <- cyp1a2
  })
}
