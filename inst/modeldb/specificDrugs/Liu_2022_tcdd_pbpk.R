Liu_2022_tcdd_pbpk <- function() {
  description <- paste(
    "PBPK/QSP (whole-body, lifetime human; original implementation in GNU",
    "Octave). 2,3,7,8-tetrachlorodibenzo-p-dioxin (TCDD) in an adult woman",
    "(Liu et al. 2022, Toxics). Single-congener instance of the dioxin-like",
    "compound (DLC) mixture PBPK framework, adapted from the Emond et al.",
    "human TCDD model and used as the framework's validation case (its",
    "Figure 2 reproduces the Emond model to R-squared = 1). Blood is a",
    "quasi-steady-state algebraic compartment; fat, rest-of-body and liver",
    "are each split into a blood subcompartment and a diffusion-limited",
    "tissue subcompartment. In the liver tissue the free congener binds",
    "reversibly to the aryl hydrocarbon receptor (AHR) and to CYP1A2;",
    "AHR-liganded congener induces CYP1A2 transcription through a Hill",
    "function (mRNA then protein), and the induced CYP1A2 both sequesters",
    "and metabolises the congener. A slow first-order blood clearance",
    "(CLURI) accounts for non-CYP1A2 elimination. Body weight and tissue",
    "volume/flow fractions are polynomial functions of age, so the model is",
    "run as a forward lifetime simulation of constant daily dietary intake",
    "with no dose records. Deterministic; no IIV and no residual error.",
    "Parameter values are the authors' published Octave code",
    "(https://github.com/pulsatility/2022-DLC-Mixture-PBPK), which the paper",
    "Section 2.6 designates as the model source."
  )
  reference <- paste(
    "Liu R, Zacharewski TR, Conolly RB, Zhang Q. A Physiologically Based",
    "Pharmacokinetic (PBPK) Modeling Framework for Mixtures of Dioxin-like",
    "Compounds. Toxics. 2022;10(11):700. doi:10.3390/toxics10110700.",
    "Physiological and biochemical parameters from the authors' Octave code",
    "DLC_human_cmd.m / DLC_human_ode.m and chemical-specific TCDD parameters",
    "from DLC_human_parameters.xlsx",
    "(https://github.com/pulsatility/2022-DLC-Mixture-PBPK; paper Tables",
    "S1-S3). Structural parent model: Emond C et al. Environ Health Perspect",
    "2005;113(12):1666-1668 and Toxicol Sci 2004;80(1):115-133."
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
    "gi",
    "lymph",
    "portal",
    "fat_blood",
    "fat",
    "rest_blood",
    "rest",
    "liver_blood",
    "dioxin_ahr",
    "dioxin_cyp1a2",
    "cyp1a2_mrna",
    "cyp1a2"
  )

  compartmentData <- list(
    gi = list(analyte = "TCDD", units = "nmol", specimen = "administration site", verified = TRUE),
    lymph = list(analyte = "TCDD", units = "nmol", specimen = "lymph", verified = TRUE),
    portal = list(analyte = "TCDD", units = "nmol", specimen = "not applicable", verified = TRUE),
    urine = list(analyte = "TCDD", units = "nmol", specimen = "urine", verified = TRUE),
    fat_blood = list(analyte = "TCDD", units = "nmol", specimen = "whole blood", verified = TRUE),
    fat = list(analyte = "TCDD", units = "nmol", specimen = "tissue", verified = TRUE),
    rest_blood = list(analyte = "TCDD", units = "nmol", specimen = "whole blood", verified = TRUE),
    rest = list(analyte = "TCDD", units = "nmol", specimen = "tissue", verified = TRUE),
    liver_blood = list(analyte = "TCDD", units = "nmol", specimen = "whole blood", verified = TRUE),
    liver = list(analyte = "TCDD", units = "nmol", specimen = "tissue", verified = TRUE),
    dioxin_ahr = list(analyte = "TCDD-AHR complex", units = "nmol", specimen = "tissue", verified = TRUE),
    dioxin_cyp1a2 = list(analyte = "TCDD-CYP1A2 complex", units = "nmol", specimen = "tissue", verified = TRUE),
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
      "simulation of lifetime dietary TCDD intake. Results in the paper are",
      "presented for women unless otherwise stated."
    ),
    dose_range = paste(
      "Constant daily dietary intake. The paper runs an EPA reference dose",
      "(RfD, 0.0007 ng/kg bw/day) and point-of-departure dose (PoD,",
      "0.02 ng/kg bw/day) lifetime scenario; MSTOT defaults to the PoD dose."
    ),
    regions = "exposure scenario (EPA RfD/PoD); dietary-intake composition from a 1994 Dutch survey",
    notes = paste(
      "Single-congener (TCDD) instance of the DLC mixture framework. The",
      "Monte Carlo population analysis (paper Figure 5) samples total AHR,",
      "CYP1A2 transcription/basal levels, CYP1A2 Emax, the body-weight",
      "polynomial and dietary intake; those distributions are not encoded as",
      "random effects here (deterministic typical-value model)."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Exposure. Constant daily dietary intake, oral (ng/kg bw/day).
    # -------------------------------------------------------------------
    MSTOT <- fixed(0.02)
    label("Daily oral dietary intake of TCDD (ng/kg bw/day)") # DLC_human_cmd.m g.MSTOT PoD = 0.02 (RfD = 7e-4); paper Section 3.1

    # -------------------------------------------------------------------
    # Chemical-specific parameters, TCDD.
    # Source: DLC_human_parameters.xlsx row '2,3,7,8-TCDD' (paper Table
    # S3); unchanged from the Emond TCDD model. All FIXED.
    # -------------------------------------------------------------------
    PL <- fixed(6.0)
    label("Liver tissue:blood partition coefficient (unitless)") # xlsx 'PL' TCDD = 6
    PF <- fixed(100.0)
    label("Fat tissue:blood partition coefficient (unitless)") # xlsx 'PF' TCDD = 100
    PRB <- fixed(1.5)
    label("Rest-of-body tissue:blood partition coefficient (unitless)") # xlsx 'PRB' TCDD = 1.5
    PALC <- fixed(0.35)
    label("Liver permeability-area coefficient, fraction of liver blood flow (unitless)") # xlsx 'PALC' = 0.35
    PAFC <- fixed(0.12)
    label("Fat permeability-area coefficient, fraction of fat blood flow (unitless)") # xlsx 'PAFC' = 0.12
    PARBC <- fixed(0.03)
    label("Rest-of-body permeability-area coefficient, fraction of RB blood flow (unitless)") # xlsx 'PARBC' = 0.03
    Fu <- fixed(1.0)
    label("Fraction of congener unbound in tissue blood (unitless)") # xlsx 'Fu' = 1
    KDAHR <- fixed(0.1)
    label("TCDD-AHR dissociation constant Kd,AHR (nM)") # xlsx 'KDAHR' TCDD = 0.1
    kbAHR <- fixed(6.0)
    label("TCDD-AHR dissociation rate constant kb,AHR (1/h)") # paper Section 2.3 'we set it to 6/h'; xlsx 'kbAHR' = 6
    KD1A2 <- fixed(40.0)
    label("TCDD-CYP1A2 dissociation constant Kd,1A2 (nM)") # xlsx 'KD1A2' = 40
    kb1A2 <- fixed(6.0)
    label("TCDD-CYP1A2 dissociation rate constant kb,1A2 (1/h)") # xlsx 'kb1A2' = 6
    kelim <- fixed(0.0011)
    label("CYP1A2-mediated hepatic elimination rate constant (1/h)") # xlsx 'kelim' TCDD = 0.0011
    KST <- fixed(0.01)
    label("GI-tract transit (loss) rate constant (1/h)") # xlsx 'KST' = 0.01
    KABS <- fixed(0.06)
    label("GI-tract absorption rate constant (1/h)") # xlsx 'KABS' TCDD = 0.06
    a <- fixed(0.7)
    label("Fraction of absorbed congener entering blood via lymph (unitless)") # xlsx 'a' = 0.7 (portal fraction b = 1-a)
    CLURI <- fixed(4.17e-08)
    label("Slow first-order (urinary/non-CYP1A2) blood clearance (L/h)") # xlsx 'CLURI' = 4.17e-8
    MW <- fixed(321.976)
    label("TCDD molecular weight (g/mol)") # xlsx 'MW' TCDD = 321.976

    # -------------------------------------------------------------------
    # Biochemical parameters, hepatic AHR / CYP1A2 (shared).
    # Source: DLC_human_cmd.m (translation/degradation/transcription of
    # CYP1A2 adjusted from Emond so basal and maximally induced
    # steady-state levels match; paper Section 2.3). All FIXED.
    # -------------------------------------------------------------------
    AHRtot <- fixed(0.35)
    label("Total hepatic AHR concentration (nM)") # DLC_human_cmd.m g.AHRtot = 0.35
    CYP1A2_1BASAL <- fixed(1600.0)
    label("Basal (and initial) free CYP1A2 protein concentration (nM)") # DLC_human_cmd.m g.CYP1A2_1BASAL = 1600
    CYP1A2_mRNA_1BASAL <- fixed(160.0)
    label("Basal (and initial) CYP1A2 mRNA concentration (nM)") # DLC_human_cmd.m g.CYP1A2_mRNA_1BASAL = 160
    CYP1A2_1EC50 <- fixed(130.0)
    label("Dioxin-AHR concentration for half-maximal CYP1A2 induction EC50 (nM)") # DLC_human_cmd.m g.CYP1A2_1EC50 = 130
    CYP1A2_1EMAX <- fixed(9300.0)
    label("Maximal fold induction of CYP1A2 transcription Emax (unitless)") # DLC_human_cmd.m g.CYP1A2_1EMAX = 9300
    hill_n <- fixed(0.6)
    label("Hill coefficient for CYP1A2 transcriptional induction (unitless)") # DLC_human_ode.m exponent 0.6 on sum(Dioxin_AHR)
    ktranscription_1A2 <- fixed(16.0)
    label("CYP1A2 mRNA transcription rate constant (nM/h)") # DLC_human_cmd.m g.ktranscription_1A2 = 16
    kdegCYP1A2_mRNA <- fixed(0.1)
    label("CYP1A2 mRNA degradation rate constant (1/h)") # DLC_human_cmd.m g.kdegCYP1A2_mRNA = 0.1
    ktranslation_1A2 <- fixed(1.25)
    label("CYP1A2 protein translation rate constant (1/h)") # DLC_human_cmd.m g.ktranslation_1A2 = 1.25
    kdegCYP1A2 <- fixed(0.125)
    label("CYP1A2 protein degradation rate constant (1/h)") # DLC_human_cmd.m g.kdegCYP1A2 = 0.125
  })

  model({
    # =================================================================
    # Fixed physiology (DLC_human_cmd.m; Emond human, Wang et al. 1997).
    # =================================================================
    QCC <- 15.36 # cardiac-output constant (L/h/kg^0.75)
    QLC <- 0.26 # liver fraction of cardiac output
    QFC <- 0.05 # fat fraction of cardiac output
    QRBC <- 1 - QFC - QLC # rest-of-body fraction of cardiac output
    VLBC <- 0.266 # liver blood volume fraction of liver volume
    VFBC <- 0.050 # fat blood volume fraction of fat volume
    VRBBC <- 0.030 # RB blood volume fraction of RB volume
    ORAL_DURATION <- 24 # daily intake spread over 24 h (continuous)

    # =================================================================
    # Time-varying anatomy. age in years; BW and volume/flow fractions
    # are polynomial functions of age (female; DLC_human_ode.m).
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

    PAF <- PAFC * QF # fat permeability-area cross product (L/h)
    PARB <- PARBC * QRB # rest-of-body permeability-area cross product (L/h)
    PAL <- PALC * QL # liver permeability-area cross product (L/h)

    kfAHR <- kbAHR / KDAHR # AHR association rate constant (1/nM/h)
    kf1A2 <- kb1A2 / KD1A2 # CYP1A2 association rate constant (1/nM/h)
    b <- 1 - a # portal fraction of absorbed congener

    # =================================================================
    # Continuous dietary input (nmol/h), from ng/kg/day and time-varying BW.
    # =================================================================
    ORAL_DOSE_RATE <- (MSTOT / MW) * BW / ORAL_DURATION

    # =================================================================
    # Blood is a quasi-steady-state algebraic compartment (no ODE).
    # =================================================================
    CA <- (QF * (fat_blood / VFB) + QRB * (rest_blood / VRBB) +
      QL * (liver_blood / VLB) + KABS * gi * a) / (QC + CLURI)

    # Liver-tissue concentrations
    CFLL <- liver / VL / PL * Fu # free + nonspecific-bound free concentration (nM)
    Dioxin_CYP1A2 <- dioxin_cyp1a2 / VL
    Dioxin_AHR <- dioxin_ahr / VL

    # =================================================================
    # ODEs (amounts in nmol unless noted). DLC_human_ode.m, M = 1.
    # =================================================================
    d/dt(gi) <- -(KST + KABS) * gi + ORAL_DOSE_RATE
    d/dt(lymph) <- KABS * gi * a # cumulative lymph-absorbed amount
    d/dt(portal) <- KABS * gi * b # cumulative portal-absorbed amount
    d/dt(urine) <- CLURI * CA # cumulative slow blood clearance

    d/dt(fat_blood) <- QF * (CA - fat_blood / VFB) -
      PAF * (fat_blood / VFB * Fu - (fat / VF) / PF * Fu)
    d/dt(fat) <- PAF * (fat_blood / VFB * Fu - (fat / VF) / PF * Fu)

    d/dt(rest_blood) <- QRB * (CA - rest_blood / VRBB) -
      PARB * (rest_blood / VRBB * Fu - rest / VRB / PRB * Fu)
    d/dt(rest) <- PARB * (rest_blood / VRBB * Fu - rest / VRB / PRB * Fu)

    d/dt(liver_blood) <- QL * (CA - liver_blood / VLB) -
      PAL * (liver_blood / VLB * Fu - CFLL) + KABS * gi * b
    d/dt(liver) <- PAL * (liver_blood / VLB * Fu - CFLL) +
      (-kelim * CFLL * (cyp1a2 + Dioxin_CYP1A2 - CYP1A2_1BASAL) / CYP1A2_1BASAL -
        kfAHR * CFLL * (AHRtot - Dioxin_AHR) + kbAHR * Dioxin_AHR -
        kf1A2 * CFLL * cyp1a2 + kb1A2 * Dioxin_CYP1A2) * VL

    d/dt(dioxin_ahr) <- (kfAHR * CFLL * (AHRtot - Dioxin_AHR) - kbAHR * Dioxin_AHR) * VL
    d/dt(dioxin_cyp1a2) <- (kf1A2 * CFLL * cyp1a2 - kb1A2 * Dioxin_CYP1A2) * VL

    d/dt(cyp1a2_mrna) <- ktranscription_1A2 *
      (1 + CYP1A2_1EMAX * Dioxin_AHR^hill_n / (CYP1A2_1EC50^hill_n + Dioxin_AHR^hill_n)) -
      kdegCYP1A2_mRNA * cyp1a2_mrna
    d/dt(cyp1a2) <- ktranslation_1A2 * cyp1a2_mrna - kdegCYP1A2 * cyp1a2 -
      kf1A2 * CFLL * cyp1a2 + kb1A2 * Dioxin_CYP1A2

    # Initial conditions: CYP1A2 and its mRNA start at basal.
    cyp1a2(0) <- CYP1A2_1BASAL
    cyp1a2_mrna(0) <- CYP1A2_mRNA_1BASAL

    # =================================================================
    # Observations (concentrations, nM).
    # =================================================================
    Cc <- CA # arterial/venous blood (CA)
    Cfat <- fat / VF # CF
    Crest <- rest / VRB # CRB
    Cliver_free <- liver / VL / PL # CLfree
    Cliver_total <- (liver + dioxin_ahr + dioxin_cyp1a2) / VL # CL
    Ccyp1a2 <- cyp1a2 # free CYP1A2 protein
  })
}
