Wang_2022_daidzein_equol_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited) for the dietary isoflavone daidzein",
    "and its gut-microbial metabolite S-equol in the adult human. Seven",
    "flow-limited tissue groups (blood, liver, fat, rapidly perfused,",
    "slowly perfused, small-intestine tissue) plus a small-intestine lumen",
    "and a large-intestine lumen; the large-intestine lumen is the microbiota",
    "compartment where daidzein is converted by capacity-limited",
    "(Michaelis-Menten) gut-microbial metabolism to dihydrodaidzein (DHD),",
    "S-equol (only in S-equol producers) and O-desmethylangolensin (O-DMA).",
    "Small-intestine and hepatic phase-II glucuronidation/sulfation of",
    "daidzein, and hepatic glucuronidation/sulfation of S-equol, are also",
    "capacity-limited. A coupled S-equol sub-model (large-intestine lumen,",
    "liver, fat, rapidly/slowly perfused, blood) receives the microbially",
    "formed S-equol via a first-order transfer from lumen to liver. All",
    "kinetic constants come from in vitro anaerobic human fecal incubations",
    "(microbial steps) and pooled human liver S9 incubations (S-equol",
    "conjugation), scaled to the whole body. The default parameter set is the",
    "S-equol PRODUCER; the nonproducer is the same structure with the S-equol",
    "sub-model inactive and different microbial DHD/O-DMA kinetics (see the",
    "vignette). Deterministic: the publication reports no IIV and no residual",
    "error. Daidzein enters as an oral dose into the small-intestine lumen.",
    sep = " "
  )
  reference <- paste(
    "Wang Q, Spenkelink B, Boonpawa R, Rietjens IMCM. Use of Physiologically",
    "Based Pharmacokinetic Modeling to Predict Human Gut Microbial Conversion",
    "of Daidzein to S-Equol. J Agric Food Chem. 2022 Jan 19;70(2):343-352.",
    "doi:10.1021/acs.jafc.1c03950. PMCID: PMC8759082. Microbial kinetic",
    "constants Table 2; S-equol conjugation kinetics Table 3; the complete",
    "Berkeley Madonna ODE listing with all physiological, partition and",
    "scaling parameters is Supporting Information 2 (PBPK model code).",
    sep = " "
  )
  vignette <- "Wang_2022_daidzein_equol"
  units <- list(
    time = "h",
    dosing = paste(
      "umol (daidzein, into the small-intestine lumen; an oral dose of",
      "D mg/kg is amt = D * BW * 1000 / 254.23 umol)"
    ),
    concentration = "umol/L (= uM; blood-to-plasma ratio taken as 1)"
  )

  # Daidzein enters the small-intestine lumen as an oral dose (the deposited
  # code sets Init ASILuDAI = ODOSEumol). All other states start at zero.
  dosing <- c(
    "si_lumen" # oral daidzein into the small-intestine lumen
  )

  # Every ODE state, in amount units (umol) unless noted. The metabolite
  # `a_*` states are cumulative amounts formed (integrators); they carry no
  # volume and no concentration and exist so the mass balance closes and the
  # conjugate output can be read (e.g. S-equol urinary excretion, Table 4).
  compartmentData <- list(
    si_lumen = list(
      analyte = "Daidzein remaining in the small-intestine lumen (unabsorbed)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    si_tissue = list(
      analyte = "Daidzein in small-intestine tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    a_si_dai7g = list(
      analyte = "Daidzein-7-O-glucuronide formed in small-intestine tissue (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_si_dai4ig = list(
      analyte = "Daidzein-4'-O-glucuronide formed in small-intestine tissue (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_si_dais = list(
      analyte = "Daidzein sulfate formed in small-intestine tissue (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    li_lumen = list(
      analyte = "Daidzein in the large-intestine lumen (microbiota compartment)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_dhd = list(
      analyte = "Dihydrodaidzein (DHD) formed by gut microbiota (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_equol_formed = list(
      analyte = "S-equol formed by gut microbiota in the large-intestine lumen (cumulative source term for the S-equol sub-model)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_odma = list(
      analyte = "O-desmethylangolensin (O-DMA) formed by gut microbiota (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    liver = list(
      analyte = "Daidzein in liver",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    a_liv_dai7g = list(
      analyte = "Daidzein-7-O-glucuronide formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_liv_dai4ig = list(
      analyte = "Daidzein-4'-O-glucuronide formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_liv_dais = list(
      analyte = "Daidzein sulfate formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    fat = list(
      analyte = "Daidzein in fat tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    rapid = list(
      analyte = "Daidzein in rapidly perfused tissue (heart, lung, brain, etc.)",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    slow = list(
      analyte = "Daidzein in slowly perfused tissue (skin, muscle, bone, etc.)",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    blood = list(
      analyte = "Daidzein in blood",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    auc_dai = list(
      analyte = "Time integral of daidzein blood AMOUNT (SI2 AUC'=AB). Equals the paper's reported daidzein AUC(0-4h) of 0.60 at 1 mg/kg; the paper labels that value uM*h/L but it is an amount integral in umol*h (see the vignette Errata)",
      units = "umol*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    li_lumen_equol = list(
      analyte = "S-equol in the large-intestine lumen",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    liver_equol = list(
      analyte = "S-equol in liver",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    a_liv_equg1 = list(
      analyte = "S-equol glucuronide-1 formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_liv_equg2 = list(
      analyte = "S-equol glucuronide-2 formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    a_liv_equs = list(
      analyte = "S-equol sulfate formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    fat_equol = list(
      analyte = "S-equol in fat tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    rapid_equol = list(
      analyte = "S-equol in rapidly perfused tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    slow_equol = list(
      analyte = "S-equol in slowly perfused tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    blood_equol = list(
      analyte = "S-equol in blood",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    ),
    auc_equol = list(
      analyte = "Time integral of S-equol blood AMOUNT (SI2 AUCEQU'=ABEQU). Equals the paper's reported S-equol AUC(0-4h) of 2.02 nmol*h at 1 mg/kg; the paper labels that value nmol*h/L but it is an amount integral in nmol*h (see the vignette Errata)",
      units = "umol*h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  # All ODE states are paper-mechanistic PBPK compartments not in the
  # canonical register (`liver` is canonical; the rest are declared here).
  paper_specific_compartments <- c(
    "si_lumen",
    "si_tissue",
    "a_si_dai7g",
    "a_si_dai4ig",
    "a_si_dais",
    "li_lumen",
    "a_dhd",
    "a_equol_formed",
    "a_odma",
    "a_liv_dai7g",
    "a_liv_dai4ig",
    "a_liv_dais",
    "fat",
    "rapid",
    "slow",
    "blood",
    "auc_dai",
    "li_lumen_equol",
    "liver_equol",
    "a_liv_equg1",
    "a_liv_equg2",
    "a_liv_equs",
    "fat_equol",
    "rapid_equol",
    "slow_equol",
    "blood_equol",
    "auc_equol"
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy adults",
    age_range = "adult",
    weight_median = "70 kg (reference body weight)",
    dose_range = paste(
      "Model predictions made for oral daidzein 0.09-3.34 mg/kg bw (literature",
      "comparison, Figure 4); Figure 5 profiles simulated at 1 mg/kg;",
      "sensitivity analysis at 1, 10 and 100 mg/kg"
    ),
    notes = paste(
      "NOT a population fit -- a deterministic in vitro-in silico PBPK.",
      "Physiological parameters are human reference values (Brown et al.",
      "1997); tissue/blood partition coefficients were computed by the QPPR",
      "of DeJongh et al. 1997. Gut-microbial kinetic constants (Vmax, Km for",
      "DHD, S-equol, O-DMA formation) were derived in vitro from pooled",
      "anaerobic fecal incubations of 6 S-equol producers and 9 nonproducers",
      "(15 volunteers screened, Table 1). S-equol glucuronidation/sulfation",
      "kinetics were derived from pooled human liver S9 fractions (25 donors,",
      "mixed gender; Table 3). Daidzein phase-II (SI and liver) kinetics and",
      "the SI/liver S9 protein yields were taken from the literature (Islam",
      "et al. 2014; Cubitt et al. 2011). The default ini() values are the",
      "S-equol PRODUCER parameter set (the complete model, with the S-equol",
      "sub-model active). The nonproducer is reproduced by setting",
      "VmaxLIEQUc = 0 (no microbial S-equol formation) and switching the",
      "microbial DHD/O-DMA constants to the nonproducer column of Table 2",
      "(VmaxLIDHDc = 0.008, KmLIDHD = 2.55, VmaxLIODMAc = 0.0007,",
      "KmLIODMA = 5.12)."
    )
  )

  ini({
    # =====================================================================
    # Physiological parameters -- SI2 code block 'Physiological parameters'
    # (Brown et al. 1997, reference 1 in SI2)
    # =====================================================================
    BW <- fixed(70)
    label("Body weight (kg)") # SI2, BW = 70 {Kg}
    VSIc <- fixed(0.0091)
    label("Small-intestine tissue volume as fraction of body weight (unitless)") # SI2, VSIc
    VLc <- fixed(0.0257)
    label("Liver tissue volume as fraction of body weight (unitless)") # SI2, VLc
    VRc <- fixed(0.0442)
    label("Rapidly perfused tissue volume as fraction of body weight (unitless)") # SI2, VRc
    VSc <- fixed(0.5428)
    label("Slowly perfused tissue volume as fraction of body weight (unitless)") # SI2, VSc
    VFc <- fixed(0.2142)
    label("Fat tissue volume as fraction of body weight (unitless)") # SI2, VFc
    VBc <- fixed(0.0790)
    label("Blood volume as fraction of body weight (unitless)") # SI2, VBc
    VMB <- fixed(0.0140)
    label("Gastrointestinal-tract (feces) contents as fraction of body weight (unitless)") # SI2, VMB (14 mL feces/kg bw)

    # ---- blood flows (SI2 'Blood flow rates', Brown et al. 1997) ---------
    QC <- fixed(347.9)
    label("Cardiac output (L/h)") # SI2, QC = 347.9 {L/h} (= 15 * BW^0.74)
    QSIc <- fixed(0.09)
    label("Fraction of cardiac output to small intestine (unitless)") # SI2, QSIc
    QLc <- fixed(0.137)
    label("Fraction of cardiac output to liver (unitless)") # SI2, QLc (0.227 - QSIc)
    QRc <- fixed(0.473)
    label("Fraction of cardiac output to rapidly perfused tissue (unitless)") # SI2, QRc (0.7 - QSIc - QLc)
    QSc <- fixed(0.248)
    label("Fraction of cardiac output to slowly perfused tissue (unitless)") # SI2, QSc (0.3 - QFc)
    QFc <- fixed(0.052)
    label("Fraction of cardiac output to fat (unitless)") # SI2, QFc

    # =====================================================================
    # Partition coefficients -- SI2 'Physicochemical parameters'
    # (QPPR of DeJongh et al. 1997, reference 2 in SI2)
    # =====================================================================
    PIDAI <- fixed(1.29)
    label("Daidzein intestine/blood partition coefficient (unitless)") # SI2, PIDAI
    PLDAI <- fixed(1.29)
    label("Daidzein liver/blood partition coefficient (unitless)") # SI2, PLDAI
    PRDAI <- fixed(1.29)
    label("Daidzein rapidly perfused tissue/blood partition coefficient (unitless)") # SI2, PRDAI
    PSDAI <- fixed(0.56)
    label("Daidzein slowly perfused tissue/blood partition coefficient (unitless)") # SI2, PSDAI
    PFDAI <- fixed(39.9)
    label("Daidzein fat/blood partition coefficient (unitless)") # SI2, PFDAI
    PLEQU <- fixed(1.83)
    label("S-equol liver/blood partition coefficient (unitless)") # SI2, PLEQU
    PREQU <- fixed(1.83)
    label("S-equol rapidly perfused tissue/blood partition coefficient (unitless)") # SI2, PREQU
    PSEQU <- fixed(0.65)
    label("S-equol slowly perfused tissue/blood partition coefficient (unitless)") # SI2, PSEQU
    PFEQU <- fixed(77.2)
    label("S-equol fat/blood partition coefficient (unitless)") # SI2, PFEQU

    # =====================================================================
    # Absorption/transfer rate constants -- SI2 'absorption/transfer rates'
    # =====================================================================
    Ka <- fixed(0.46)
    label("Absorption rate of daidzein from small-intestine lumen to tissue (1/h)") # SI2, Ka (reference 3: Steensma et al. 2004)
    Kb <- fixed(4.56)
    label("Transfer rate of daidzein from large-intestine lumen to liver (1/h)") # SI2, Kb (reference 3)
    Ksl <- fixed(1.16)
    label("Transfer rate of daidzein from small- to large-intestine lumen (to feces route) (1/h)") # SI2, Ksl (reference 4: Kimura & Higaki 2002)
    Kll <- fixed(4.56)
    label("Transfer rate of S-equol from large-intestine lumen to liver (1/h)") # SI2, Kll (reference 3)

    # =====================================================================
    # Scaling factors for in vitro-to-in vivo Vmax scaling -- SI2
    # =====================================================================
    S9SI <- fixed(38.6)
    label("Small-intestine S9 protein yield (mg S9 protein/g intestine)") # SI2, S9SI (reference 5: Cubitt et al. 2011)
    VLS9 <- fixed(143)
    label("Liver S9 protein yield (mg S9 protein/g liver)") # SI2, VLS9 (reference 7; 108 cytosolic + 35 microsomal, Cubitt et al. 2011)

    # =====================================================================
    # Small-intestine daidzein phase-II kinetics -- SI2 (Islam et al. 2014)
    # unscaled Vmax {nmol/min/mg S9 protein}; Km {umol/L}
    # =====================================================================
    VmaxSIDAI7Gc <- fixed(0.2)
    label("Unscaled Vmax, daidzein-7-O-glucuronide by small intestine (nmol/min/mg S9)") # SI2, VmaxSIDAI7Gc (reference 6)
    VmaxSIDAI4iGc <- fixed(0.1)
    label("Unscaled Vmax, daidzein-4'-O-glucuronide by small intestine (nmol/min/mg S9)") # SI2, VmaxSIDAI4iGc (reference 6)
    VmaxSIDAISc <- fixed(0.02)
    label("Unscaled Vmax, daidzein sulfate by small intestine (nmol/min/mg S9)") # SI2, VmaxSIDAISc (reference 6)
    KmSIDAI7G <- fixed(2.7)
    label("Km, daidzein-7-O-glucuronide by small intestine (umol/L)") # SI2, KmSIDAI7G (reference 6)
    KmSIDAI4iG <- fixed(2.9)
    label("Km, daidzein-4'-O-glucuronide by small intestine (umol/L)") # SI2, KmSIDAI4iG (reference 6)
    KmSIDAIS <- fixed(0.35)
    label("Km, daidzein sulfate by small intestine (umol/L)") # SI2, KmSIDAIS (reference 6)

    # =====================================================================
    # Large-intestine (microbial) daidzein kinetics -- Table 2, S-equol
    # PRODUCER column; unscaled Vmax {umol/h/g feces}; Km {umol/L}
    # =====================================================================
    VmaxLIDHDc <- fixed(0.024)
    label("Unscaled Vmax, DHD formation by gut microbiota (umol/h/g feces)") # Table 2, producers DHD Vmax 0.024; SI2 VmaxLIDHDc
    VmaxLIEQUc <- fixed(0.009)
    label("Unscaled Vmax, S-equol formation by gut microbiota (umol/h/g feces)") # Table 2, producers S-equol Vmax 0.009; SI2 VmaxLIEQUc. Set to 0 for nonproducers.
    VmaxLIODMAc <- fixed(0.001)
    label("Unscaled Vmax, O-DMA formation by gut microbiota (umol/h/g feces)") # Table 2, producers O-DMA Vmax 0.001; SI2 VmaxLIODMAc
    KmLIDHD <- fixed(6.239)
    label("Km, DHD formation by gut microbiota (umol/L)") # SI2, KmLIDHD (Table 2 producers DHD Km 6.24)
    KmLIEQU <- fixed(7.243)
    label("Km, S-equol formation by gut microbiota (umol/L)") # SI2, KmLIEQU (Table 2 producers S-equol Km 7.24)
    KmLIODMA <- fixed(18.070)
    label("Km, O-DMA formation by gut microbiota (umol/L)") # SI2, KmLIODMA (Table 2 producers O-DMA Km 18.07)

    # =====================================================================
    # Liver daidzein phase-II kinetics -- SI2 (Islam et al. 2014)
    # unscaled Vmax {nmol/min/mg S9 protein}; Km {umol/L}
    # =====================================================================
    VmaxLDAI7Gc <- fixed(1.0)
    label("Unscaled Vmax, daidzein-7-O-glucuronide by liver (nmol/min/mg S9)") # SI2, VmaxLDAI7Gc (reference 6)
    VmaxLDAI4iGc <- fixed(0.2)
    label("Unscaled Vmax, daidzein-4'-O-glucuronide by liver (nmol/min/mg S9)") # SI2, VmaxLDAI4iGc (reference 6)
    VmaxLDAISc <- fixed(0.02)
    label("Unscaled Vmax, daidzein sulfate by liver (nmol/min/mg S9)") # SI2, VmaxLDAISc (reference 6)
    KmLDAI7G <- fixed(18.9)
    label("Km, daidzein-7-O-glucuronide by liver (umol/L)") # SI2, KmLDAI7G (reference 6)
    KmLDAI4iG <- fixed(72.1)
    label("Km, daidzein-4'-O-glucuronide by liver (umol/L)") # SI2, KmLDAI4iG (reference 6)
    KmLDAIS <- fixed(0.82)
    label("Km, daidzein sulfate by liver (umol/L)") # SI2, KmLDAIS (reference 6)

    # =====================================================================
    # Liver S-equol phase-II kinetics -- Table 3; from human liver S9
    # incubations. unscaled Vmax {nmol/min/mg S9 protein}; Km {umol/L}
    # =====================================================================
    VmaxLEQUG1c <- fixed(4.62)
    label("Unscaled Vmax, S-equol glucuronide-1 by liver (nmol/min/mg S9)") # Table 3, S-equol glucuronide-1 Vmax 4.62; SI2 VmaxLEQUG1c
    VmaxLEQUG2c <- fixed(0.61)
    label("Unscaled Vmax, S-equol glucuronide-2 by liver (nmol/min/mg S9)") # Table 3, S-equol glucuronide-2 Vmax 0.61; SI2 VmaxLEQUG2c
    VmaxLEQUSc <- fixed(9.24)
    label("Unscaled Vmax, S-equol sulfate by liver (nmol/min/mg S9)") # Table 3, S-equol sulfate Vmax 9.24; SI2 VmaxLEQUSc
    KmLEQUG1 <- fixed(20.28)
    label("Km, S-equol glucuronide-1 by liver (umol/L)") # Table 3, S-equol glucuronide-1 Km 20.28; SI2 KmLEQUG1
    KmLEQUG2 <- fixed(29.39)
    label("Km, S-equol glucuronide-2 by liver (umol/L)") # Table 3, S-equol glucuronide-2 Km 29.39; SI2 KmLEQUG2
    KmLEQUS <- fixed(6.50)
    label("Km, S-equol sulfate by liver (umol/L)") # Table 3, S-equol sulfate Km 6.50; SI2 KmLEQUS
  })

  model({
    # =====================================================================
    # Scaling calculations (SI2)
    # =====================================================================
    # ---- tissue volumes (L) ---------------------------------------------
    VSI <- VSIc * BW
    VL <- VLc * BW
    VR <- VRc * BW
    VS <- VSc * BW
    VF <- VFc * BW
    VB <- VBc * BW
    # ---- gram-of-tissue-per-kg factors used in the S9 Vmax scaling ------
    SI <- VSIc * 1000
    L <- VLc * 1000
    # ---- blood flows (L/h) ----------------------------------------------
    QSI <- QSIc * QC
    QL <- QLc * QC
    QR <- QRc * QC
    QS <- QSc * QC
    QF <- QFc * QC

    # ---- scaled Vmax values (umol/h) ------------------------------------
    # Small-intestine S9: nmol/min/mg -> umol/h over the whole SI S9 pool.
    VmaxSIDAI7G <- VmaxSIDAI7Gc / 1000 * 60 * S9SI * SI * BW
    VmaxSIDAI4iG <- VmaxSIDAI4iGc / 1000 * 60 * S9SI * SI * BW
    VmaxSIDAIS <- VmaxSIDAISc / 1000 * 60 * S9SI * SI * BW
    # Large-intestine microbiota: umol/h/g feces -> umol/h over fecal mass.
    VmaxLIDHD <- VmaxLIDHDc * 1000 * VMB * BW
    VmaxLIEQU <- VmaxLIEQUc * 1000 * VMB * BW
    VmaxLIODMA <- VmaxLIODMAc * 1000 * VMB * BW
    # Liver S9: nmol/min/mg -> umol/h over the whole liver S9 pool.
    VmaxLDAI7G <- VmaxLDAI7Gc / 1000 * 60 * VLS9 * L * BW
    VmaxLDAI4iG <- VmaxLDAI4iGc / 1000 * 60 * VLS9 * L * BW
    VmaxLDAIS <- VmaxLDAISc / 1000 * 60 * VLS9 * L * BW
    VmaxLEQUG1 <- VmaxLEQUG1c / 1000 * 60 * VLS9 * L * BW
    VmaxLEQUG2 <- VmaxLEQUG2c / 1000 * 60 * VLS9 * L * BW
    VmaxLEQUS <- VmaxLEQUSc / 1000 * 60 * VLS9 * L * BW

    # =====================================================================
    # Concentrations (umol/L). CV* are venous-equilibrium concentrations
    # leaving each tissue (tissue concentration / partition coefficient).
    # =====================================================================
    CSIDAI <- si_tissue / VSI
    CVSIDAI <- CSIDAI / PIDAI
    CLIDAI <- li_lumen / (VMB * BW)
    CLDAI <- liver / VL
    CVLDAI <- CLDAI / PLDAI
    CF <- fat / VF
    CVF <- CF / PFDAI
    CR <- rapid / VR
    CVR <- CR / PRDAI
    CS <- slow / VS
    CVS <- CS / PSDAI
    CB <- blood / VB

    CLEQU <- liver_equol / VL
    CVLEQU <- CLEQU / PLEQU
    CFEQU <- fat_equol / VF
    CVFEQU <- CFEQU / PFEQU
    CREQU <- rapid_equol / VR
    CVREQU <- CREQU / PREQU
    CSEQU <- slow_equol / VS
    CVSEQU <- CSEQU / PSEQU
    CBEQU <- blood_equol / VB

    # =====================================================================
    # Capacity-limited (Michaelis-Menten) formation rates (umol/h)
    # =====================================================================
    rSI7G <- VmaxSIDAI7G * CVSIDAI / (KmSIDAI7G + CVSIDAI)
    rSI4iG <- VmaxSIDAI4iG * CVSIDAI / (KmSIDAI4iG + CVSIDAI)
    rSIS <- VmaxSIDAIS * CVSIDAI / (KmSIDAIS + CVSIDAI)
    rLIDHD <- VmaxLIDHD * CLIDAI / (KmLIDHD + CLIDAI)
    rLIEQU <- VmaxLIEQU * CLIDAI / (KmLIEQU + CLIDAI)
    rLIODMA <- VmaxLIODMA * CLIDAI / (KmLIODMA + CLIDAI)
    rLD7G <- VmaxLDAI7G * CVLDAI / (KmLDAI7G + CVLDAI)
    rLD4iG <- VmaxLDAI4iG * CVLDAI / (KmLDAI4iG + CVLDAI)
    rLDS <- VmaxLDAIS * CVLDAI / (KmLDAIS + CVLDAI)
    rEQUG1 <- VmaxLEQUG1 * CVLEQU / (KmLEQUG1 + CVLEQU)
    rEQUG2 <- VmaxLEQUG2 * CVLEQU / (KmLEQUG2 + CVLEQU)
    rEQUS <- VmaxLEQUS * CVLEQU / (KmLEQUS + CVLEQU)

    # =====================================================================
    # Daidzein main-model ODEs (SI2 'Main model calculations/dynamics')
    # =====================================================================
    d/dt(si_lumen) <- -Ka * si_lumen - Ksl * si_lumen
    d/dt(si_tissue) <- Ka * si_lumen + QSI * (CB - CVSIDAI) - rSI7G - rSI4iG - rSIS
    d/dt(a_si_dai7g) <- rSI7G
    d/dt(a_si_dai4ig) <- rSI4iG
    d/dt(a_si_dais) <- rSIS
    d/dt(li_lumen) <- Ksl * si_lumen - rLIDHD - rLIEQU - rLIODMA - Kb * li_lumen
    d/dt(a_dhd) <- rLIDHD
    d/dt(a_equol_formed) <- rLIEQU
    d/dt(a_odma) <- rLIODMA
    d/dt(liver) <- QL * CB + QSI * CVSIDAI - (QL + QSI) * CVLDAI -
      rLD7G - rLD4iG - rLDS + Kb * li_lumen
    d/dt(a_liv_dai7g) <- rLD7G
    d/dt(a_liv_dai4ig) <- rLD4iG
    d/dt(a_liv_dais) <- rLDS
    d/dt(fat) <- QF * (CB - CVF)
    d/dt(rapid) <- QR * (CB - CVR)
    d/dt(slow) <- QS * (CB - CVS)
    d/dt(blood) <- (QL + QSI) * CVLDAI + QF * CVF + QS * CVS + QR * CVR - QC * CB
    d/dt(auc_dai) <- blood

    # =====================================================================
    # S-equol sub-model ODEs (SI2 'Sub-model calculations/dynamics: equol')
    # =====================================================================
    d/dt(li_lumen_equol) <- rLIEQU - Kll * li_lumen_equol
    d/dt(liver_equol) <- Kll * li_lumen_equol + (QL + QSI) * CBEQU -
      (QSI + QL) * CVLEQU - rEQUG1 - rEQUG2 - rEQUS
    d/dt(a_liv_equg1) <- rEQUG1
    d/dt(a_liv_equg2) <- rEQUG2
    d/dt(a_liv_equs) <- rEQUS
    d/dt(fat_equol) <- QF * (CBEQU - CVFEQU)
    d/dt(rapid_equol) <- QR * (CBEQU - CVREQU)
    d/dt(slow_equol) <- QS * (CBEQU - CVSEQU)
    d/dt(blood_equol) <- (QL + QSI) * CVLEQU + QF * CVFEQU +
      QR * CVREQU + QS * CVSEQU - QC * CBEQU
    d/dt(auc_equol) <- blood_equol

    # =====================================================================
    # Observations. Blood-to-plasma ratio taken as 1, so blood == plasma.
    # =====================================================================
    Cc <- CB
    Cc_equol <- CBEQU
    # No residual-error or IIV model: Wang 2022 is a deterministic
    # in vitro-in silico PBPK and reports neither.
  })
}
