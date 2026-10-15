deBruijn_2022_bileacid_bosentan_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, flow-limited) for bile-acid homeostasis and",
    "bosentan-induced cholestasis in the adult human (de Bruijn 2022). A",
    "lumped bile-acid pool represented by glycochenodeoxycholic acid (GCDCA)",
    "circulates enterohepatically through eight flow-limited compartments",
    "(gall bladder, liver, intestinal lumen, intestinal tissue, fat, rapidly",
    "perfused, slowly perfused, blood). GCDCA is synthesised de novo in the",
    "liver (Ks), actively effluxed from liver to bile canaliculi by BSEP",
    "(Michaelis-Menten, Vmax scaled from an in vitro vesicular assay via an",
    "absolute BSEP-abundance IVIVE), split 50/50 between the common bile duct",
    "(direct to intestinal lumen) and gall-bladder storage, reabsorbed from",
    "the lumen (ka), and excreted faecally (Kf = Ks). The gall bladder",
    "contracts at three simulated daytime meals (08:00, 12:00, 16:00),",
    "emptying its contents into the intestinal lumen. Two coupled drug",
    "sub-models (each five flow-limited compartments plus liver metabolism)",
    "describe the BSEP inhibitor bosentan and its active desmethyl metabolite",
    "RO 47-8634; their free intrahepatic concentrations non-competitively",
    "inhibit BSEP-mediated BA efflux (modulation factor",
    "1 + CVLbos/Kibos + CVLdes/KiDES), so the systemic BA pool rises under",
    "bosentan treatment. With no bosentan dose the model reduces to healthy",
    "BA homeostasis. Deterministic: the publication reports no IIV and no",
    "residual error. Interindividual variability is explored two ways -- a",
    "log-normal BSEP abundance (aBSEP) and an empirical total-pool scaling",
    "factor (sens) -- both left as overridable fixed parameters; see the",
    "validation vignette. All but one parameter were derived experimentally.",
    sep = " "
  )
  reference <- paste(
    "de Bruijn VMP, Rietjens IMCM, Bouwmeester H. Population pharmacokinetic",
    "model to generate mechanistic insights in bile acid homeostasis and",
    "drug-induced cholestasis. Arch Toxicol. 2022 Sep;96(9):2541-2558.",
    "doi:10.1007/s00204-022-03345-8. PMCID: PMC9352636. The complete Berkeley",
    "Madonna ODE listing with all physiological, partition, scaling and",
    "kinetic parameters is Supplementary file II (PBK model code); the",
    "physicochemical properties are Table 1.",
    sep = " "
  )
  vignette <- "deBruijn_2022_bile_acid_cholestasis"
  units <- list(
    time = "h",
    dosing = paste(
      "umol (bosentan, into bos_stomach; an oral dose of D mg is",
      "amt = D * 1000 / 551.6 umol, with f(bos_stomach) = Fa). The gall-bladder",
      "bile-acid content is an initial condition, not a dose."
    ),
    concentration = "umol/L (= uM)"
  )

  # Every ODE state is in amount units (umol). The `*_formed` and `*_bile`
  # drug states are cumulative integrators (no volume, no concentration): they
  # exist so the mass balance closes and the metabolite source term can be
  # read. All states are paper-mechanistic PBK compartments outside the
  # canonical register, so each is declared in paper_specific_compartments.
  compartmentData <- list(
    ba_gallbladder = list(
      analyte = "GCDCA stored in the gall bladder",
      units = "umol",
      specimen = "bile",
      verified = TRUE
    ),
    ba_liver = list(analyte = "GCDCA in liver", units = "umol", specimen = "tissue", verified = TRUE),
    ba_lumen = list(
      analyte = "GCDCA in the intestinal lumen",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    ba_intestine = list(
      analyte = "GCDCA in intestinal tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    ba_fat = list(analyte = "GCDCA in fat tissue", units = "umol", specimen = "tissue", verified = TRUE),
    ba_rapid = list(
      analyte = "GCDCA in rapidly perfused tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    ba_slow = list(
      analyte = "GCDCA in slowly perfused tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    ba_blood = list(analyte = "GCDCA in blood", units = "umol", specimen = "whole blood", verified = TRUE),
    bos_stomach = list(
      analyte = "Bosentan remaining in the stomach (absorption site)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    bos_liver = list(analyte = "Bosentan in liver", units = "umol", specimen = "tissue", verified = TRUE),
    bos_oh_formed = list(
      analyte = "Hydroxy-bosentan RO 48-5033 formed in liver (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    bos_des_formed = list(
      analyte = "Desmethyl-bosentan RO 47-8634 formed in liver (cumulative source term for the metabolite sub-model)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    bos_bile = list(
      analyte = "Bosentan excreted via bile (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    bos_fat = list(analyte = "Bosentan in fat tissue", units = "umol", specimen = "tissue", verified = TRUE),
    bos_rapid = list(
      analyte = "Bosentan in rapidly perfused tissue (lumped with intestine)",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    bos_slow = list(
      analyte = "Bosentan in slowly perfused tissue (lumped with gall bladder + lumen)",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    bos_blood = list(analyte = "Bosentan in blood", units = "umol", specimen = "whole blood", verified = TRUE),
    des_liver = list(
      analyte = "Desmethyl-bosentan RO 47-8634 in liver",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    des_bile = list(
      analyte = "Desmethyl-bosentan RO 47-8634 excreted via bile (cumulative)",
      units = "umol",
      specimen = "not applicable",
      verified = TRUE
    ),
    des_fat = list(
      analyte = "Desmethyl-bosentan RO 47-8634 in fat tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    des_rapid = list(
      analyte = "Desmethyl-bosentan RO 47-8634 in rapidly perfused tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    des_slow = list(
      analyte = "Desmethyl-bosentan RO 47-8634 in slowly perfused tissue",
      units = "umol",
      specimen = "tissue",
      verified = TRUE
    ),
    des_blood = list(
      analyte = "Desmethyl-bosentan RO 47-8634 in blood",
      units = "umol",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  paper_specific_compartments <- c(
    "ba_gallbladder",
    "ba_liver",
    "ba_lumen",
    "ba_intestine",
    "ba_fat",
    "ba_rapid",
    "ba_slow",
    "ba_blood",
    "bos_stomach",
    "bos_liver",
    "bos_oh_formed",
    "bos_des_formed",
    "bos_bile",
    "bos_fat",
    "bos_rapid",
    "bos_slow",
    "bos_blood",
    "des_liver",
    "des_bile",
    "des_fat",
    "des_rapid",
    "des_slow",
    "des_blood"
  )

  # Bosentan enters the stomach (absorption site) as an oral dose; the
  # gall-bladder bile-acid content is set as an initial condition in model().
  dosing <- c(
    "bos_stomach"
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy adults; bosentan-treated and bosentan-induced cholestasis explored",
    age_range = "adult",
    weight_median = "70 kg (reference body weight)",
    dose_range = "Bosentan 500 mg orally twice a day (08:00 and 20:00); three daytime meals at 08:00, 12:00, 16:00 drive gall-bladder contraction.",
    notes = paste(
      "NOT a population fit -- a deterministic in vitro / in silico PBK model.",
      "Physiological volumes and blood flows are human reference values (Brown",
      "et al. 1997; gall-bladder volume Van Erpecum et al. 1992). GCDCA, bosentan",
      "and RO 47-8634 tissue:blood partition coefficients were computed by the",
      "QPPR of Rodgers & Rowland 2006 (QIVIVE toolbox, Punt et al. 2020). BSEP",
      "Vmax/Km are from an Sf9 vesicular transport assay with physiological",
      "cholesterol (Kis et al. 2009); bosentan metabolism Vmax/Km and",
      "non-saturable clearance from human liver microsomes (Sato et al. 2018);",
      "BSEP inhibition constants Ki from Fattinger et al. 2001 (taurocholate,",
      "assumed equal for GCDCA). GCDCA ka and bosentan/metabolite absorption +",
      "biliary-excretion rate constants were fitted to in vivo data (Hepner &",
      "Demers 1977; Ponz de Leon et al. 1978; Weber et al. 1999). Two",
      "interindividual-variability scenarios are left as overridable fixed",
      "parameters: (1) BSEP abundance aBSEP is log-normal with meanlog -0.26,",
      "sdlog 0.403 (Burt et al. 2016 meta-analysis of Caucasian hepatic",
      "transporter abundances, truncated at +/- 3 SD, i.e. aBSEP in 0.23-2.58",
      "pmol/10^6 hepatocytes); the deterministic reference uses aBSEP = 0.839.",
      "(2) The empirical pool-scaling factor sens (0.5 / 1 / 1.5) multiplies the",
      "gall-bladder initial content, de novo synthesis, faecal excretion and",
      "fasting plasma concentration to span the reported between-subject range",
      "in total BA pool size (reference 3020 umol gall-bladder content; total",
      "pool ~3079 umol at sens = 1)."
    )
  )

  ini({
    # =====================================================================
    # Physiological parameters -- Supplementary file II (Suppl. II),
    # fractions of body weight and of cardiac output (Brown et al. 1997;
    # gall bladder Van Erpecum et al. 1992). All fixed.
    # =====================================================================
    BW <- fixed(70)
    label("Body weight (kg)") # Suppl. II, BW = 70 {Kg}
    VFc <- fixed(0.214)
    label("Fat tissue volume fraction of BW (unitless)") # Suppl. II, VFc
    VLc <- fixed(0.026)
    label("Liver tissue volume fraction of BW (unitless)") # Suppl. II, VLc
    VRc <- fixed(0.054)
    label("Rapidly perfused tissue volume fraction of BW (unitless)") # Suppl. II, VRc
    VSc <- fixed(0.6033)
    label("Slowly perfused tissue volume fraction of BW (unitless)") # Suppl. II, VSc
    VBc <- fixed(0.079)
    label("Blood volume fraction of BW (unitless)") # Suppl. II, VBc
    VIc <- fixed(0.009)
    label("Intestinal tissue volume fraction of BW (unitless)") # Suppl. II, VIc
    VGc <- fixed(0.0007)
    label("Gall-bladder volume fraction of BW (unitless)") # Suppl. II, VGc (Van Erpecum 1992)
    VLuc <- fixed(0.014)
    label("Intestinal-lumen volume fraction of BW (unitless)") # Suppl. II, VLuc

    QFc <- fixed(0.052)
    label("Fraction of cardiac output to fat (unitless)") # Suppl. II, QFc
    QLc <- fixed(0.046)
    label("Fraction of cardiac output to liver, excl. portal vein (unitless)") # Suppl. II, QLc
    QSc <- fixed(0.248)
    label("Fraction of cardiac output to slowly perfused tissue (unitless)") # Suppl. II, QSc
    QRc <- fixed(0.473)
    label("Fraction of cardiac output to rapidly perfused tissue (unitless)") # Suppl. II, QRc
    QIc <- fixed(0.181)
    label("Fraction of cardiac output to intestines (unitless)") # Suppl. II, QIc

    # =====================================================================
    # GCDCA physicochemical parameters -- Suppl. II. Partition coefficients
    # are tissue:plasma numerators (QPPR of Rodgers & Rowland 2006) divided
    # by the blood:plasma ratio RGCDCA in model(); stored here as numerators
    # to match the deposited code (e.g. PF = 0.05 / RGCDCA).
    # =====================================================================
    RGCDCA <- fixed(0.55)
    label("GCDCA blood:plasma ratio (1 - Hct, assumption)") # Suppl. II, RGCDCA; Table 1
    PFnum <- fixed(0.05)
    label("GCDCA fat:plasma partition numerator (unitless)") # Suppl. II, PF = 0.05 / RGCDCA
    PLnum <- fixed(0.09)
    label("GCDCA liver:plasma partition numerator (unitless)") # Suppl. II, PL = 0.09 / RGCDCA
    PRnum <- fixed(0.125)
    label("GCDCA rapidly perfused:plasma partition numerator (unitless)") # Suppl. II, PR = 0.125 / RGCDCA
    PSnum <- fixed(0.19)
    label("GCDCA slowly perfused:plasma partition numerator (unitless)") # Suppl. II, PS = 0.19 / RGCDCA
    PGnum <- fixed(0.16)
    label("GCDCA gut:plasma partition numerator (unitless)") # Suppl. II, PG = 0.16 / RGCDCA

    # =====================================================================
    # GCDCA kinetic parameters -- Suppl. II.
    # =====================================================================
    Ka <- fixed(1.047)
    label("GCDCA absorption rate, intestinal lumen to liver (1/h)") # Suppl. II, ka (fitted; Hepner 1977, De Leon 1978)
    KsPerH <- fixed(46.8)
    label("GCDCA de novo hepatic synthesis rate (umol/h, before sens scaling)") # Suppl. II, Ks = 0.78 * 60 (Kullak-Ublick 2004)
    VmaxBSEPc <- fixed(5.848)
    label("BSEP Vmax for GCDCA efflux (umol/min/mg BSEP)") # Suppl. II, VmaxBSEPc (Kis 2009, +cholesterol)
    KmBSEP <- fixed(4.3)
    label("BSEP Michaelis constant for GCDCA efflux (umol/L)") # Suppl. II, KmBSEP (Kis 2009)
    aBSEP <- fixed(0.839)
    label("BSEP protein abundance (pmol/10^6 hepatocytes; reference individual)") # Suppl. II, deterministic aBSEP = 0.839 (Burt 2016)
    MWBSEP <- fixed(140000)
    label("BSEP molecular weight (g/mol)") # Suppl. II, 140 kDa
    Hep <- fixed(99)
    label("Hepatocellularity (10^6 hepatocytes/g liver)") # Suppl. II, Hep (Barter 2007)
    liverGPerKg <- fixed(20)
    label("Liver weight per kg body weight (g/kg)") # Suppl. II, WL = 20 * BW (Soars 2002)
    QIbfrac <- fixed(0.5)
    label("Fraction of biliary GCDCA going directly to intestinal lumen (unitless)") # Suppl. II, QIb (Molino 1986)
    CBfs0 <- fixed(2.4)
    label("Fasting systemic plasma GCDCA concentration (umol/L, before sens scaling)") # Suppl. II, CBfs (Garcia-Canaveras 2012)
    Gdose0 <- fixed(3020)
    label("Gall-bladder GCDCA content at start (umol, before sens scaling)") # Suppl. II, Gdose = 3020 (Sips 2018)
    sens <- fixed(1)
    label("Empirical total-BA-pool scaling factor (unitless; 0.5 / 1 / 1.5 scenarios)") # Suppl. II, sens (reference individual = 1)

    # =====================================================================
    # BSEP inhibition constants -- Suppl. II (Fattinger 2001; taurocholate
    # values assumed equal for GCDCA).
    # =====================================================================
    Kibos <- fixed(12)
    label("BSEP inhibition constant for bosentan (umol/L)") # Suppl. II, Kibos (Fattinger 2001)
    KiDES <- fixed(8.5)
    label("BSEP inhibition constant for RO 47-8634 (umol/L)") # Suppl. II, KiDES (Fattinger 2001)

    # =====================================================================
    # Bosentan physicochemical parameters -- Suppl. II (partition numerators
    # over the blood:plasma ratio Rbos).
    # =====================================================================
    Rbos <- fixed(0.6)
    label("Bosentan blood:plasma ratio (unitless)") # Suppl. II, Rbos; Table 1 (EMA 2004, Meyer 1996)
    PFbosNum <- fixed(0.05)
    label("Bosentan fat:plasma partition numerator (unitless)") # Suppl. II, PFbos = 0.05 / Rbos
    PLbosNum <- fixed(0.11)
    label("Bosentan liver:plasma partition numerator (unitless)") # Suppl. II, PLbos = 0.11 / Rbos
    PRbosNum <- fixed(0.14)
    label("Bosentan rapidly perfused:plasma partition numerator (unitless)") # Suppl. II, PRbos = 0.14 / Rbos
    PSbosNum <- fixed(0.21)
    label("Bosentan slowly perfused:plasma partition numerator (unitless)") # Suppl. II, PSbos = 0.21 / Rbos

    # =====================================================================
    # RO 47-8634 (desmethyl bosentan) physicochemical parameters -- Suppl. II.
    # =====================================================================
    RDES <- fixed(0.55)
    label("RO 47-8634 blood:plasma ratio (1 - Hct, assumption)") # Suppl. II, RDES; Table 1
    PFdesNum <- fixed(0.06)
    label("RO 47-8634 fat:plasma partition numerator (unitless)") # Suppl. II, PFDES = 0.06 / RDES
    PLdesNum <- fixed(0.15)
    label("RO 47-8634 liver:plasma partition numerator (unitless)") # Suppl. II, PLDES = 0.15 / RDES
    PRdesNum <- fixed(0.18)
    label("RO 47-8634 rapidly perfused:plasma partition numerator (unitless)") # Suppl. II, PRDES = 0.18 / RDES
    PSdesNum <- fixed(0.30)
    label("RO 47-8634 slowly perfused:plasma partition numerator (unitless)") # Suppl. II, PSDES = 0.30 / RDES

    # =====================================================================
    # Bosentan / metabolite kinetic parameters -- Suppl. II. Absorption and
    # biliary rate constants fitted to Weber 1999; metabolism from Sato 2018
    # human liver microsomes; scaled with MPPGL and liver weight in model().
    # =====================================================================
    kabos <- fixed(0.130)
    label("Bosentan absorption rate, stomach to liver (1/h)") # Suppl. II, kabos (fitted; Weber 1999)
    kbilebos <- fixed(23.660)
    label("Bosentan biliary excretion rate (1/h)") # Suppl. II, kbilebos (fitted; Weber 1999)
    kbileDES <- fixed(133.924)
    label("RO 47-8634 biliary excretion rate (1/h)") # Suppl. II, kbileDES (fitted; Weber 1999)
    Fa <- fixed(0.5)
    label("Bosentan fraction absorbed (unitless)") # Suppl. II, Fa (Weber 1996)
    MPPGL <- fixed(32)
    label("Microsomal protein per gram of liver (mg/g)") # Suppl. II, MPPGL (Barter 2007)
    VmaxOHc <- fixed(16.4)
    label("Unscaled Vmax, RO 48-5033 (hydroxy-bosentan) formation (pmol/min/mg microsomal protein)") # Suppl. II, VmaxOHc (Sato 2018)
    VmaxDESc <- fixed(7.53)
    label("Unscaled Vmax, RO 47-8634 (desmethyl-bosentan) formation (pmol/min/mg microsomal protein)") # Suppl. II, VmaxDESc (Sato 2018)
    KmOH <- fixed(6.4)
    label("Km, RO 48-5033 formation (umol/L)") # Suppl. II, KmOH (Sato 2018)
    KmDES <- fixed(4.8)
    label("Km, RO 47-8634 formation (umol/L)") # Suppl. II, KmDES (Sato 2018)
    CLOHc <- fixed(0.158)
    label("Unscaled non-saturable clearance, RO 48-5033 pathway (uL/min/mg microsomal protein)") # Suppl. II, CLOHc (Sato 2018)
    CLDESc <- fixed(0.273)
    label("Unscaled non-saturable clearance, RO 47-8634 pathway (uL/min/mg microsomal protein)") # Suppl. II, CLDESc (Sato 2018)
  })

  model({
    # =====================================================================
    # Scaled volumes (L) and blood flows (L/h). QC = 15 * BW^0.74 (Brown
    # et al. 1997). Bosentan / metabolite lump intestine into rapidly
    # perfused and gall bladder + lumen into slowly perfused (Suppl. II).
    # =====================================================================
    VF <- VFc * BW
    VL <- VLc * BW
    VR <- VRc * BW
    VS <- VSc * BW
    VB <- VBc * BW
    VI <- VIc * BW
    VLu <- VLuc * BW

    QC <- 15 * BW^0.74
    QF <- QFc * QC
    QL <- QLc * QC
    QS <- QSc * QC
    QR <- QRc * QC
    QI <- QIc * QC

    # bosentan / metabolite lumped tissue groups
    VRbos <- (VRc + VIc) * BW
    VSbos <- (VSc + VGc + VLuc) * BW
    QRbos <- (QRc + QIc) * QC

    # =====================================================================
    # Partition coefficients (numerator / blood:plasma ratio).
    # =====================================================================
    PF <- PFnum / RGCDCA
    PL <- PLnum / RGCDCA
    PR <- PRnum / RGCDCA
    PS <- PSnum / RGCDCA
    PG <- PGnum / RGCDCA

    PFbos <- PFbosNum / Rbos
    PLbos <- PLbosNum / Rbos
    PRbos <- PRbosNum / Rbos
    PSbos <- PSbosNum / Rbos

    PFdes <- PFdesNum / RDES
    PLdes <- PLdesNum / RDES
    PRdes <- PRdesNum / RDES
    PSdes <- PSdesNum / RDES

    # =====================================================================
    # Sens-scaled pool quantities and BSEP Vmax IVIVE.
    # =====================================================================
    Ks <- KsPerH * sens # de novo synthesis (umol/h)
    Kf <- Ks # faecal excretion equals de novo synthesis (mass balance)
    CBfs <- CBfs0 * sens # fasting plasma GCDCA (umol/L)
    Gdose <- Gdose0 * sens # gall-bladder GCDCA content at start (umol)

    WL <- liverGPerKg * BW # liver weight (g)
    # scaling factor (mg BSEP/entire liver); 60 min/h, 1e-9 pg->mg (Suppl. II Eq. 2)
    SF <- aBSEP * MWBSEP * Hep * WL * 60 * 1e-9
    VmaxBSEP <- VmaxBSEPc * SF # umol/h/entire liver

    # bosentan-metabolite liver Vmax (pmol/min/mg -> umol/h/entire liver) and
    # non-saturable clearance (uL/min/mg -> L/h/entire liver).
    VmaxOH <- VmaxOHc * MPPGL * WL * 60 * 1e-6
    VmaxDES <- VmaxDESc * MPPGL * WL * 60 * 1e-6
    CLOH <- CLOHc * MPPGL * WL * 60 * 1e-6
    CLDES <- CLDESc * MPPGL * WL * 60 * 1e-6

    # =====================================================================
    # Concentrations. CV* are venous-equilibrium concentrations leaving a
    # tissue (tissue concentration / partition coefficient).
    # =====================================================================
    CL <- ba_liver / VL
    CVL <- CL / PL
    CI <- ba_intestine / VI
    CVI <- CI / PG
    CF <- ba_fat / VF
    CVF <- CF / PF
    CR <- ba_rapid / VR
    CVR <- CR / PR
    CS <- ba_slow / VS
    CVS <- CS / PS
    CB <- ba_blood / VB

    CLbos <- bos_liver / VL
    CVLbos <- CLbos / PLbos
    CFbos <- bos_fat / VF
    CVFbos <- CFbos / PFbos
    CRbos <- bos_rapid / VRbos
    CVRbos <- CRbos / PRbos
    CSbos <- bos_slow / VSbos
    CVSbos <- CSbos / PSbos
    CBbos <- bos_blood / VB

    CLdes <- des_liver / VL
    CVLdes <- CLdes / PLdes
    CFdes <- des_fat / VF
    CVFdes <- CFdes / PFdes
    CRdes <- des_rapid / VRbos
    CVRdes <- CRdes / PRdes
    CSdes <- des_slow / VSbos
    CVSdes <- CSdes / PSdes
    CBdes <- des_blood / VB

    # =====================================================================
    # BSEP-mediated GCDCA efflux with non-competitive inhibition by free
    # intrahepatic bosentan and RO 47-8634 (Suppl. II; main-text Eq. 7).
    # With no bosentan dose CVLbos = CVLdes = 0 so VmaxBSEPapp = VmaxBSEP.
    # =====================================================================
    VmaxBSEPapp <- VmaxBSEP / (1 + CVLbos / Kibos + CVLdes / KiDES)
    bsepEfflux <- VmaxBSEPapp * CVL / (KmBSEP + CVL)
    QGbfrac <- 1 - QIbfrac # fraction of biliary GCDCA stored in gall bladder

    # =====================================================================
    # Gall-bladder contraction at three simulated daytime meals (08:00,
    # 12:00, 16:00). The deposited model empties the entire gall-bladder
    # content instantaneously (a Dirac pulse); rxode2 integrates continuous
    # ODEs, so the maintainers reproduce the full emptying with a short
    # high-rate window (~0.25 h, rate 30/h -> >99.9% emptied) gated to the
    # three daytime meal times via the time-of-day tod. This matches the
    # paper's stated meal schedule and overnight fasting; see the vignette.
    # Downstream users must solve on a grid fine enough to resolve the
    # 0.25 h emptying window.
    # =====================================================================
    tod <- t - 24 * floor(t / 24)
    meal1 <- (tod >= 8) * (tod < 8.25)
    meal2 <- (tod >= 12) * (tod < 12.25)
    meal3 <- (tod >= 16) * (tod < 16.25)
    atMeal <- meal1 + meal2 + meal3
    gbEmpty <- atMeal * 30 * ba_gallbladder # umol/h, maintainer approximation of the instantaneous pulse

    # =====================================================================
    # Bile-acid (GCDCA) sub-model ODEs (Suppl. II).
    # =====================================================================
    d/dt(ba_gallbladder) <- -gbEmpty + bsepEfflux * QGbfrac
    d/dt(ba_lumen) <- gbEmpty + bsepEfflux * QIbfrac - Kf - Ka * ba_lumen
    d/dt(ba_liver) <- QL * (CB - CVL) - bsepEfflux + Ks + Ka * ba_lumen
    d/dt(ba_intestine) <- QI * (CB - CVI)
    d/dt(ba_fat) <- QF * (CB - CVF)
    d/dt(ba_rapid) <- QR * (CB - CVR)
    d/dt(ba_slow) <- QS * (CB - CVS)
    d/dt(ba_blood) <- QF * CVF + QL * CVL + QS * CVS + QR * CVR + QI * CVI -
      (QF + QL + QS + QR + QI) * CB
    ba_gallbladder(0) <- Gdose

    # =====================================================================
    # Bosentan sub-model ODEs (Suppl. II). rMOH / rMDES are the hydroxy- and
    # desmethyl-metabolite formation rates; rMDES is the source term for the
    # RO 47-8634 sub-model. rbilebos is biliary excretion.
    # =====================================================================
    rMOH <- VmaxOH * CVLbos / (KmOH + CVLbos) + CLOH * CVLbos
    rMDES <- VmaxDES * CVLbos / (KmDES + CVLbos) + CLDES * CVLbos
    rbilebos <- kbilebos * bos_liver

    d/dt(bos_stomach) <- -kabos * bos_stomach
    d/dt(bos_liver) <- kabos * bos_stomach + QL * (CBbos - CVLbos) -
      rMOH - rMDES - rbilebos
    d/dt(bos_oh_formed) <- rMOH
    d/dt(bos_des_formed) <- rMDES
    d/dt(bos_bile) <- rbilebos
    d/dt(bos_fat) <- QF * (CBbos - CVFbos)
    d/dt(bos_rapid) <- QRbos * (CBbos - CVRbos)
    d/dt(bos_slow) <- QS * (CBbos - CVSbos)
    d/dt(bos_blood) <- QL * CVLbos + QF * CVFbos + QS * CVSbos + QRbos * CVRbos -
      (QL + QF + QS + QRbos) * CBbos
    f(bos_stomach) <- Fa

    # =====================================================================
    # RO 47-8634 (desmethyl bosentan) sub-model ODEs (Suppl. II). Formation
    # from bosentan enters the liver (+rMDES); rbileDES is biliary excretion.
    # =====================================================================
    rbileDES <- kbileDES * des_liver

    d/dt(des_liver) <- QL * (CBdes - CVLdes) + rMDES - rbileDES
    d/dt(des_bile) <- rbileDES
    d/dt(des_fat) <- QF * (CBdes - CVFdes)
    d/dt(des_rapid) <- QRbos * (CBdes - CVRdes)
    d/dt(des_slow) <- QS * (CBdes - CVSdes)
    d/dt(des_blood) <- QL * CVLdes + QF * CVFdes + QS * CVSdes + QRbos * CVRdes -
      (QL + QF + QS + QRbos) * CBdes

    # =====================================================================
    # Observations (umol/L = uM). Whole-blood concentrations converted to
    # plasma via the blood:plasma ratio; the bile-acid output adds the
    # fasting systemic plasma concentration CBfs (Suppl. II CBtot).
    # Deterministic model: de Bruijn 2022 reports no residual error / IIV.
    # =====================================================================
    Cc <- CB / RGCDCA + CBfs # systemic plasma GCDCA
    Cc_bosentan <- CBbos / Rbos # plasma bosentan
    Cc_desmethyl <- CBdes / RDES # plasma RO 47-8634
  })
}
