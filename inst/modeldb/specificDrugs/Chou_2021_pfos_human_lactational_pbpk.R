Chou_2021_pfos_human_lactational_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, lactational; original implementation in R/mrgsolve).",
    "Perfluorooctane sulfonate (PFOS) in the lactating woman and her",
    "breastfed neonate (Chou and Lin 2021, Environ Health Perspect). The",
    "mother is a two-compartment gastrointestinal tract feeding the liver,",
    "plus plasma, liver, fat, mammary gland, rest of body, a milk",
    "compartment and a three-subcompartment kidney (kidney blood, proximal",
    "tubule cells, filtrate) with glomerular filtration, saturable",
    "basolateral (Oat1/Oat3) and apical (Oatp1a1) reabsorption, passive",
    "diffusion and first-order efflux. PFOS diffuses from the mammary gland",
    "into milk and is transferred to the neonate by breastfeeding; the",
    "neonatal submodel is gut, plasma, liver, rest of body and a",
    "three-subcompartment kidney. Maternal and neonatal physiology change",
    "with postnatal week/month (Loccisano et al. 2013; Yang et al. 2019).",
    "PFOS is not metabolised; only the unbound fraction exchanges.",
    "Deterministic, no random effects. Parameter values are the calibrated",
    "(as-run) set from the authors' openly published mrgsolve code and fit",
    "object (https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Human/HMod.R",
    "object LHumanPBPK and LFit_H.rds)."
  )
  reference <- paste(
    "Chou WC, Lin Z. Development of a Gestational and Lactational",
    "Physiologically Based Pharmacokinetic (PBPK) Model for",
    "Perfluorooctane Sulfonate (PFOS) in Rats and Humans and Its",
    "Implications in the Derivation of Health-Based Toxicity Values.",
    "Environ Health Perspect. 2021;129(3):037004. doi:10.1289/EHP7671.",
    "Structure and lactational growth equations from the Methods",
    "(Loccisano et al. 2013; Yang et al. 2019) and the authors' published",
    "code (ModFit/Human/HMod.R, object LHumanPBPK); calibrated parameters",
    "from LFit_H.rds (as-run) and Table 2 (human, lactating column)",
    "(https://github.com/KSUICCM/PFOS-Ges-Lac)."
  )
  vignette <- "Chou_2021_pfos_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c(
    "kidney_blood",
    "ptc",
    "filtrate",
    "fat",
    "mammary",
    "rest",
    "milk",
    "feces",
    "plasma_neo",
    "liver_neo",
    "rest_neo",
    "kidney_blood_neo",
    "ptc_neo",
    "filtrate_neo",
    "gi_neo",
    "urine_neo",
    "feces_neo",
    "auc_plasma",
    "auc_plasma_neo"
  )

  compartmentData <- list(
    stomach = list(analyte = "PFOS", units = "mg", specimen = "administration site", verified = TRUE),
    intestine = list(analyte = "PFOS", units = "mg", specimen = "administration site", verified = TRUE),
    plasma = list(analyte = "PFOS", units = "mg", specimen = "plasma", verified = TRUE),
    liver = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_blood = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    ptc = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    filtrate = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    fat = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    mammary = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    rest = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    milk = list(analyte = "PFOS", units = "mg", specimen = "milk", verified = TRUE),
    urine = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    feces = list(analyte = "PFOS", units = "mg", specimen = "faeces", verified = TRUE),
    plasma_neo = list(analyte = "PFOS", units = "mg", specimen = "plasma", verified = TRUE),
    liver_neo = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    rest_neo = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_blood_neo = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    ptc_neo = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    filtrate_neo = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    gi_neo = list(analyte = "PFOS", units = "mg", specimen = "administration site", verified = TRUE),
    urine_neo = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    feces_neo = list(analyte = "PFOS", units = "mg", specimen = "faeces", verified = TRUE),
    auc_plasma = list(analyte = "PFOS", units = "mg/L*h", specimen = "not applicable", verified = TRUE),
    auc_plasma_neo = list(analyte = "PFOS", units = "mg/L*h", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Maternal body weight at delivery excluding the fetus, BW0. The",
        "authors' code uses BW0 = 68.05 kg. Maternal weight then declines",
        "over postpartum weeks (BW = 0.0014*WK^2 - 0.1227*WK + BW0) and",
        "plasma, fat and mammary volumes and blood flows follow postpartum",
        "growth equations; glomerular filtration and milk production have",
        "their own postnatal-time equations. The breastfed neonate grows by",
        "its own postnatal-month body-weight equation independent of WT."
      ),
      source_name = "BW0"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 6L,
    age_range = "lactating women and breastfed neonates, birth to ~6 months postpartum",
    weight_range = "BW0 68.05 kg (authors' code default)",
    sex_female_pct = 100,
    disease_state = paste(
      "Lactating women and breastfed neonates. Calibrated against maternal",
      "plasma, neonatal plasma and breast-milk PFOS biomonitoring data",
      "(Karrman et al. 2007; von Ehrenstein et al. 2009; Fromme et al.",
      "2010; Lee et al. 2018) and evaluated against further cohorts."
    ),
    dose_range = paste(
      "Chronic dietary intake. Population-specific estimates 0.35-2.3",
      "ng/kg/day (Loccisano et al. 2013), applied to the stomach."
    ),
    regions = "worldwide human biomonitoring",
    notes = paste(
      "Calibrated (as-run) parameters are the sensitivity-selected subset",
      "in LFit_H.rds: maternal Free, apical RAFapi and the",
      "mammary-to-milk permeability PAMilkC. All other values are the",
      "code's $PARAM defaults. Milk feeding stops after month 6",
      "(KMilk = 0). A single neonate is modelled."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific and milk-transfer parameters for PFOS in the
    # lactating woman and breastfed neonate. As-run values: the code's
    # LHumanPBPK $PARAM defaults with the sensitivity-selected subset
    # overridden by LFit_H.rds (calibrated, log-scale). Source of record:
    # authors' published code (https://github.com/KSUICCM/PFOS-Ges-Lac,
    # ModFit/Human/HMod.R, LFit_H.rds). Every value is FIXED; the model is
    # deterministic.
    # ---------------------------------------------------------------------
    lfu <- fixed(log(0.0371512)); label("Free fraction of PFOS in maternal plasma (unitless)")  # LFit_H calibrated 'Free' = 0.0371512 (Table 2 human lactating 0.037)

    lk0     <- fixed(log(1.000));   label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # HMod LHumanPBPK 'K0C' = 1
    lkabs   <- fixed(log(2.120));   label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # HMod LHumanPBPK 'KabsC' = 2.12
    lkunabs <- fixed(log(7.05e-5)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # HMod LHumanPBPK 'KunabsC' = 7.05e-5
    lge     <- fixed(log(3.510));   label("Gastric emptying rate constant GEC (1/h per BW^0.25)")              # HMod LHumanPBPK 'GEC' = 3.510

    lkp_liver   <- fixed(log(2.03)); label("Maternal liver-to-plasma partition coefficient PL (unitless)")           # HMod LHumanPBPK 'PL' = 2.03
    lkp_kidney  <- fixed(log(1.26)); label("Maternal kidney-to-plasma partition coefficient PK (unitless)")          # HMod LHumanPBPK 'PK' = 1.26
    lkp_fat     <- fixed(log(0.13)); label("Maternal fat-to-plasma partition coefficient PF (unitless)")             # HMod LHumanPBPK 'PF' = 0.13
    lkp_mammary <- fixed(log(0.16)); label("Maternal mammary-to-plasma partition coefficient PM (unitless)")         # HMod LHumanPBPK 'PM' = 0.16
    lkp_rest    <- fixed(log(0.20)); label("Maternal rest-of-body-to-plasma partition coefficient PRest (unitless)") # HMod LHumanPBPK 'PRest' = 0.20

    lpmilk_m <- fixed(log(1.9));       label("Milk-to-mammary-gland partition coefficient PMilkM (unitless)")   # HMod LHumanPBPK 'PMilkM' = 1.9
    lpamilk  <- fixed(log(0.0028404)); label("Mammary-to-milk permeability-area product PAMilkC (L/h per kg)")  # LFit_H calibrated 'PAMilkC' = 0.0028404 (Table 2 human lactating 0.0028)

    lkbile  <- fixed(log(1.3e-4)); label("Maternal biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # HMod LHumanPBPK 'KbileC' = 1.3e-4
    lkurine <- fixed(log(0.096));  label("Maternal urinary elimination rate constant KurineC (1/h per BW^-0.25)") # HMod LHumanPBPK 'KurineC' = 0.096

    lvmax_baso_invitro   <- fixed(log(479));   label("Maternal in vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)")  # HMod LHumanPBPK 'Vmax_baso_invitro' = 479
    lvmax_apical_invitro <- fixed(log(51803)); label("Maternal in vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")         # HMod LHumanPBPK 'Vmax_apical_invitro' = 51803
    lkm_baso             <- fixed(log(20.1));  label("Maternal Michaelis constant, basolateral transporters (mg/L)")         # HMod LHumanPBPK 'Km_baso' = 20.1
    lkm_apical           <- fixed(log(64.4));  label("Maternal Michaelis constant, apical transporters (mg/L)")              # HMod LHumanPBPK 'Km_apical' = 64.4
    lrafbaso             <- fixed(log(1));         label("Relative activity factor, basolateral transporters (unitless)")    # HMod LHumanPBPK 'RAFbaso' = 1
    lrafapi              <- fixed(log(0.525075));  label("Relative activity factor, apical transporters (unitless)")         # LFit_H calibrated 'RAFapi' = 0.525075 (Table 2 human lactating 0.525)
    lkdif                <- fixed(log(0.001));  label("Maternal diffusion rate, kidney blood to proximal tubule cells (L/h)") # HMod LHumanPBPK 'Kdif' = 0.001
    lkefflux             <- fixed(log(0.150));  label("Maternal efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)")  # HMod LHumanPBPK 'KeffluxC' = 0.150

    # ---- Breastfed-neonate chemical parameters ----
    lfu_neo       <- fixed(log(0.014));  label("Free fraction of PFOS in neonatal plasma (unitless)")             # HMod LHumanPBPK 'Free_neo' = 0.014
    lkp_liver_neo <- fixed(log(2.03));   label("Neonatal liver-to-plasma partition coefficient PL_neo (unitless)") # HMod LHumanPBPK 'PL_neo' = 2.03
    lkp_rest_neo  <- fixed(log(0.20));   label("Neonatal rest-of-body-to-plasma partition coefficient PRest_neo (unitless)") # HMod LHumanPBPK 'PRest_neo' = 0.20
    lkabs_neo     <- fixed(log(2.120));  label("Neonatal small-intestine absorption rate constant KabsC_neo (1/h per BW^-0.25)") # HMod LHumanPBPK 'KabsC_neo' = 2.12
    lkbile_neo    <- fixed(log(1.3e-4)); label("Neonatal biliary elimination rate constant KbileC_neo (1/h per BW^-0.25)")       # HMod LHumanPBPK 'KbileC_neo' = 1.3e-4
    lkurine_neo   <- fixed(log(0.001));  label("Neonatal urinary elimination rate constant KurineC_neo (1/h per BW^-0.25)")      # HMod LHumanPBPK 'KurineC_neo' = 0.001
    lkdif_neo     <- fixed(log(0.001));  label("Neonatal diffusion rate, kidney blood to proximal tubule cells (L/h)")           # HMod LHumanPBPK 'Kdif_neo' = 0.001
    lkefflux_neo  <- fixed(log(0.150));  label("Neonatal efflux rate constant, tubule cells to plasma KeffluxC_neo (1/h per BW^-0.25)") # HMod LHumanPBPK 'KeffluxC_neo' = 0.150

    # Residual error placeholder: no residual-error model is reported for
    # this deterministic PBPK. FIXED, carries no information. See Errata.
    propSd <- fixed(0.30); label("Proportional residual error placeholder, plasma (fraction)")  # not reported in Chou 2021; placeholder only
  })

  model({
    # =====================================================================
    # 0. Back-transforms
    # =====================================================================
    fu       <- exp(lfu)
    K0C      <- exp(lk0)
    KabsC    <- exp(lkabs)
    KunabsC  <- exp(lkunabs)
    GEC      <- exp(lge)
    PL       <- exp(lkp_liver)
    PK       <- exp(lkp_kidney)
    PF       <- exp(lkp_fat)
    PM       <- exp(lkp_mammary)
    PRest    <- exp(lkp_rest)
    PMilkM   <- exp(lpmilk_m)
    PAMilkC  <- exp(lpamilk)
    KbileC   <- exp(lkbile)
    KurineC  <- exp(lkurine)
    Vmax_baso_invitro   <- exp(lvmax_baso_invitro)
    Vmax_apical_invitro <- exp(lvmax_apical_invitro)
    Km_baso   <- exp(lkm_baso)
    Km_apical <- exp(lkm_apical)
    RAFbaso   <- exp(lrafbaso)
    RAFapi    <- exp(lrafapi)
    Kdif      <- exp(lkdif)
    KeffluxC  <- exp(lkefflux)
    fu_neo       <- exp(lfu_neo)
    PL_neo       <- exp(lkp_liver_neo)
    PRest_neo    <- exp(lkp_rest_neo)
    KabsC_neo    <- exp(lkabs_neo)
    KbileC_neo   <- exp(lkbile_neo)
    KurineC_neo  <- exp(lkurine_neo)
    Kdif_neo     <- exp(lkdif_neo)
    KeffluxC_neo <- exp(lkefflux_neo)

    # =====================================================================
    # 1. Fixed physiology (HMod LHumanPBPK $PARAM). A single neonate.
    # =====================================================================
    BW0     <- WT        # maternal delivery weight excluding fetus (default 68.05 kg)
    Htc     <- 0.44      # HMod LHumanPBPK Htc
    QCC     <- 16.4      # HMod LHumanPBPK QCC (L/h/kg^0.75)
    QLC     <- 0.25      # HMod LHumanPBPK QLC
    QKC     <- 0.141     # HMod LHumanPBPK QKC
    QMC     <- 0.027     # HMod LHumanPBPK QMC
    QFC     <- 0.052     # HMod LHumanPBPK QFC
    VLC     <- 0.026     # HMod LHumanPBPK VLC
    VKC     <- 0.004     # HMod LHumanPBPK VKC
    VMC     <- 0.0062    # HMod LHumanPBPK VMC
    VMilk   <- 0.25      # HMod LHumanPBPK residual milk volume (L)
    VFilC   <- 8.4e-4    # HMod LHumanPBPK VFilC (L/kg BW)
    VPTCC   <- 1.35e-4   # HMod LHumanPBPK VPTCC (L/kg kidney)
    protein <- 2.0e-6    # HMod LHumanPBPK protein (mg protein/PTC)
    MW      <- 500.126   # HMod LHumanPBPK MW (g/mol)
    Vmax_baso_invitro_neo   <- 479   # HMod LHumanPBPK neonatal Vmax_baso_invitro_neo
    Vmax_apical_invitro_neo <- 51803 # HMod LHumanPBPK neonatal Vmax_apical_invitro_neo
    Km_baso_neo   <- 20.1  # HMod LHumanPBPK neonatal Km_baso_neo
    Km_apical_neo <- 64.4  # HMod LHumanPBPK neonatal Km_apical_neo

    # =====================================================================
    # 2. Postnatal time and maternal growth equations
    # =====================================================================
    PND <- t / 24
    WK  <- PND / 7
    Mon <- PND / 30
    Yr  <- Mon / 12
    PMA <- 40 + WK

    BW      <- 0.0014 * WK^2 - 0.1227 * WK + BW0
    VPlasC  <- (1e-4) * WK^2 - 0.0024 * WK + 0.0469
    VFC     <- (-7e-4) * WK + 0.3026
    VFC_0   <- 0.3026
    VFC_F   <- (-7e-4) * 48 + 0.3026
    KMilk_0 <- (1e-9) * PND^3 - (1e-6) * PND^2 + 0.0002 * PND + 0.0211
    if (Mon <= 6) {
      KMilk <- KMilk_0
    } else {
      KMilk <- 0
    }
    VMT <- (-3e-6) * PND^2 + 0.0006 * PND + 0.0059
    GFR <- 0.0259 * WK^2 - 0.4369 * WK + 7.7972

    VPlas <- VPlasC * BW
    VL    <- VLC * BW0
    VK    <- VKC * BW0
    VF    <- VFC * BW
    VF_0  <- VFC_0 * BW
    VM    <- (1 + (VMT / VMC)) * (VMC * BW0)
    VM_0  <- (VMC + 0.05) * BW0
    MK    <- VKC * BW * 1000
    VKb   <- VK * 0.16
    VFil  <- VFilC * BW
    VPTC  <- MK * VPTCC
    VRest <- 0.93 * BW - (VL + VK + VM + VF + VPlas + VMilk)

    # ---- Maternal blood flows ----
    QC_0 <- QCC * BW0^0.75 * (1 - Htc)
    QC_P <- QC_0 * (1 + (QFC * ((0.4 / VFC_F) - 1)) + (QMC * ((0.05 / VMC) - 1)))
    QF_0 <- QFC * QC_0 * (0.4 / VFC_F)
    QF   <- QF_0 * (VF / VF_0)
    QM_0 <- QMC * QC_0 * (0.05 / VMC)
    QM   <- QM_0 * (VM / VM_0)
    QL   <- QLC * QC_0
    QK   <- QKC * QC_0
    QC   <- QC_P + (QF - QF_0) + (QM - QM_0)
    QRest <- QC - (QL + QK + QM + QF)

    # ---- Maternal scaled kinetic parameters ----
    Kabs   <- KabsC * BW^(-0.25)
    Kunabs <- KunabsC * BW^(-0.25)
    Kbile  <- KbileC * BW^(-0.25)
    Kurine <- KurineC * BW^(-0.25)
    GE     <- GEC * BW^(-0.25)
    K0     <- K0C * BW^(-0.25)
    PTC          <- VKC * 6e7 * 1000
    Vmax_basoC   <- Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC <- Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso    <- Vmax_basoC * BW^0.75
    Vmax_apical  <- Vmax_apicalC * BW^0.75
    Kefflux      <- KeffluxC * BW^(-0.25)

    # =====================================================================
    # 3. Neonatal physiology (Yang et al. 2019)
    # =====================================================================
    BW_neo <- (-0.0183) * Mon^2 + 0.7732 * Mon + 3.7711
    if (WK > 27) {
      Htc_neo <- (0.744 * WK + 24.656) / 100
    } else {
      Htc_neo <- 0.37
    }
    VPlasC_neo <- (4e-6) * Yr^3 - 0.0003 * Yr^2 + 0.0059 * Yr + 0.05
    VPlas_neo  <- VPlasC_neo * BW_neo
    VL_neo     <- (0.036 * BW_neo * 1000 + 17.851) / 1000
    VK_neo     <- ((0.0034 * BW_neo * 1000 + 24.36) + (0.0036 * BW_neo * 1000 + 2.207)) / 1000
    MK_neo     <- VKC * BW_neo * 1000
    VPTC_neo   <- MK_neo * VPTCC
    VKb_neo    <- VK_neo * 0.16
    VFil_neo   <- VFilC * BW_neo
    VRest_neo  <- 0.93 * BW_neo - (VL_neo + VPlas_neo + VK_neo + VPTC_neo + VFil_neo)
    QC_neo     <- (214.81 * BW_neo + 76.57) * 0.06
    QK_neo     <- 0.221 * VK_neo * 1000 + 0.572
    QL_neo     <- 0.026 * QC_neo / 0.06
    QRest_neo  <- QC_neo - (QL_neo + QK_neo)
    Kabs_neo   <- KabsC_neo * BW_neo^(-0.25)

    Vmax_basoC_neo   <- Vmax_baso_invitro_neo * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC_neo <- Vmax_apical_invitro_neo * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso_neo    <- Vmax_basoC_neo * BW_neo^0.75
    Vmax_apical_neo  <- Vmax_apicalC_neo * BW_neo^0.75
    Kefflux_neo <- KeffluxC_neo * BW_neo^(-0.25)
    PAMilk      <- PAMilkC * BW_neo^0.75
    Kurine_neo  <- KurineC_neo * BW_neo^(-0.25)
    Kbile_neo   <- KbileC_neo * BW_neo^(-0.25)
    GFR_neo     <- 0.00108 * exp(0.1328 * PMA)

    # =====================================================================
    # 4. Concentrations (mg/L)
    # =====================================================================
    CVL    <- liver / (VL * PL)
    CVK    <- kidney_blood / VKb
    CVM    <- mammary / (VM * PM)
    CVF    <- fat / (VF * PF)
    CVRest <- rest / ((VRest + 1e-7) * PRest)
    CPlas_free <- plasma / VPlas
    CPlas      <- CPlas_free / fu
    CKb        <- kidney_blood / VKb
    CPTC       <- ptc / VPTC
    CFil       <- filtrate / (VFil + 1e-7)
    CMilk      <- milk / VMilk
    CKid       <- CVK * PK
    CM         <- mammary / VM
    CF         <- fat / VF
    CRest      <- rest / VRest
    CL         <- liver / VL

    CVL_neo    <- liver_neo / (VL_neo * PL_neo)
    CVK_neo    <- kidney_blood_neo / VKb_neo
    CVRest_neo <- rest_neo / ((VRest_neo + 1e-7) * PRest_neo)
    CPlas_free_neo <- plasma_neo / (VPlas_neo + 1e-7)
    CPlas_neo  <- CPlas_free_neo / fu_neo
    CKb_neo    <- kidney_blood_neo / VKb_neo
    CPTC_neo   <- ptc_neo / VPTC_neo
    CFil_neo   <- filtrate_neo / VFil_neo
    CL_neo     <- liver_neo / (VL_neo + 1e-7)

    # =====================================================================
    # 5. Fluxes (mg/h)
    # =====================================================================
    RA_baso   <- Vmax_baso * (kidney_blood / VKb) / (Km_baso + kidney_blood / VKb)
    RA_apical <- Vmax_apical * (filtrate / (VFil + 1e-7)) / (Km_apical + filtrate / (VFil + 1e-7))
    Rdif      <- Kdif * (CKb - CPTC)
    RAefflux  <- Kefflux * ptc
    RCI       <- CPlas * GFR * fu
    Rtrans    <- KMilk * CMilk

    RA_baso_neo   <- Vmax_baso_neo * (kidney_blood_neo / VKb_neo) / (Km_baso_neo + kidney_blood_neo / VKb_neo)
    RA_apical_neo <- Vmax_apical_neo * (filtrate_neo / VFil_neo) / (Km_apical_neo + filtrate_neo / VFil_neo)
    Rdif_neo      <- Kdif_neo * (CKb_neo - CPTC_neo)
    RAefflux_neo  <- Kefflux_neo * ptc_neo
    RCI_neo       <- CPlas_neo * GFR_neo * fu_neo

    # =====================================================================
    # 6. Maternal ODEs
    # =====================================================================
    d/dt(stomach)   <- -K0 * stomach - GE * stomach
    d/dt(intestine) <- GE * stomach - Kabs * intestine - Kunabs * intestine
    d/dt(liver)     <- QL * (CPlas - CVL) * fu - Kbile * liver + Kabs * intestine + K0 * stomach
    d/dt(fat)       <- QF * (CPlas - CVF) * fu
    d/dt(mammary)   <- QM * (CPlas - CVM) * fu - PAMilk * (CVM * fu - CMilk / PMilkM)
    d/dt(milk)      <- PAMilk * (CVM * fu - CMilk / PMilkM) - Rtrans
    d/dt(rest)      <- QRest * (CPlas - CVRest) * fu
    d/dt(kidney_blood) <- QK * (CPlas - CVK) * fu - RCI - Rdif - RA_baso
    d/dt(ptc)          <- Rdif + RA_apical + RA_baso - RAefflux
    d/dt(filtrate)     <- RCI - RA_apical - Kurine * filtrate
    d/dt(plasma)    <- QF * CVF * fu + QM * CVM * fu + QRest * CVRest * fu +
                       QK * CVK * fu + QL * CVL * fu - QC * CPlas * fu + RAefflux
    d/dt(urine)     <- Kurine * filtrate
    d/dt(feces)     <- Kbile * liver + Kunabs * intestine

    # =====================================================================
    # 7. Neonatal ODEs (gut fed by suckled milk Rtrans)
    # =====================================================================
    d/dt(gi_neo)           <- Rtrans - Kabs_neo * gi_neo
    d/dt(liver_neo)        <- QL_neo * (CPlas_neo - CVL_neo) * fu_neo + Kabs_neo * gi_neo - Kbile_neo * liver_neo * fu_neo
    d/dt(rest_neo)         <- QRest_neo * (CPlas_neo - CVRest_neo) * fu_neo
    d/dt(kidney_blood_neo) <- QK_neo * (CPlas_neo - CVK_neo) * fu_neo - RCI_neo - Rdif_neo - RA_baso_neo
    d/dt(ptc_neo)          <- Rdif_neo + RA_apical_neo + RA_baso_neo - RAefflux_neo
    d/dt(filtrate_neo)     <- RCI_neo - RA_apical_neo - filtrate_neo * Kurine_neo
    d/dt(plasma_neo)       <- QRest_neo * CVRest_neo * fu_neo + QL_neo * CVL_neo * fu_neo +
                              QK_neo * CVK_neo * fu_neo - QC_neo * CPlas_neo * fu_neo + RAefflux_neo
    d/dt(urine_neo)        <- filtrate_neo * Kurine_neo
    d/dt(feces_neo)        <- Kbile_neo * liver_neo * fu_neo

    # =====================================================================
    # 8. Exposure accumulators
    # =====================================================================
    d/dt(auc_plasma)     <- CPlas
    d/dt(auc_plasma_neo) <- CPlas_neo

    # =====================================================================
    # 9. Observations
    # =====================================================================
    Cc        <- CPlas
    Cliver    <- CL
    Ckidney   <- CKid
    Cmilk     <- CMilk
    Cc_neo    <- CPlas_neo
    Cliver_neo <- CL_neo

    Cc ~ prop(propSd)
  })
}
