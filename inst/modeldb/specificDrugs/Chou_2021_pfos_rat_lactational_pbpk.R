Chou_2021_pfos_rat_lactational_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, lactational; original implementation in R/mrgsolve).",
    "Perfluorooctane sulfonate (PFOS) in the lactating Sprague-Dawley dam",
    "and her nursing litter (Chou and Lin 2021, Environ Health Perspect).",
    "The dam is a two-compartment gastrointestinal tract feeding the liver,",
    "plus plasma, liver, fat, mammary gland, rest of body, a milk",
    "compartment and a three-subcompartment kidney (kidney blood, proximal",
    "tubule cells, filtrate) with glomerular filtration, saturable",
    "basolateral (Oat1/Oat3) and apical (Oatp1a1) reabsorption, passive",
    "diffusion and first-order efflux. PFOS diffuses from the mammary",
    "gland into milk and is suckled by the whole litter (N = 8); the",
    "nursing-pup submodel is gut, plasma, liver, rest of body and a",
    "three-subcompartment kidney, with pup urine returned to the dam gut",
    "(coprophagy). Dam and pup physiology grow with postnatal day",
    "(Loccisano et al. 2012; Mirfazaelian and Fisher 2007). PFOS is not",
    "metabolised; only the unbound fraction exchanges. Deterministic, no",
    "random effects. Parameter values are the calibrated (as-run) set from",
    "the authors' openly published mrgsolve code and fit object",
    "(https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Rat/RMod.R object",
    "LRatPBPK and LFit_R.rds)."
  )
  reference <- paste(
    "Chou WC, Lin Z. Development of a Gestational and Lactational",
    "Physiologically Based Pharmacokinetic (PBPK) Model for",
    "Perfluorooctane Sulfonate (PFOS) in Rats and Humans and Its",
    "Implications in the Derivation of Health-Based Toxicity Values.",
    "Environ Health Perspect. 2021;129(3):037004. doi:10.1289/EHP7671.",
    "Structure and lactational growth equations from the Methods",
    "(Loccisano et al. 2012; Mirfazaelian and Fisher 2007) and the",
    "authors' published code (ModFit/Rat/RMod.R, object LRatPBPK);",
    "calibrated parameters from LFit_R.rds (as-run) and Table 2 (rat,",
    "lactating column) (https://github.com/KSUICCM/PFOS-Ges-Lac)."
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
    "plasma_pup",
    "liver_pup",
    "rest_pup",
    "kidney_blood_pup",
    "ptc_pup",
    "filtrate_pup",
    "gi_pup",
    "urine_pup",
    "feces_pup",
    "auc_plasma",
    "auc_plasma_pup"
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
    plasma_pup = list(analyte = "PFOS", units = "mg", specimen = "plasma", verified = TRUE),
    liver_pup = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    rest_pup = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    kidney_blood_pup = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    ptc_pup = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    filtrate_pup = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    gi_pup = list(analyte = "PFOS", units = "mg", specimen = "administration site", verified = TRUE),
    urine_pup = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    feces_pup = list(analyte = "PFOS", units = "mg", specimen = "faeces", verified = TRUE),
    auc_plasma = list(analyte = "PFOS", units = "mg/L*h", specimen = "not applicable", verified = TRUE),
    auc_plasma_pup = list(analyte = "PFOS", units = "mg/L*h", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Dam body weight at the start of lactation (PND0 = GD22), BW0. The",
        "authors' code uses BW0 = 0.247 kg (Shirley 1984). Dam weight then",
        "grows linearly with postnatal day (BW = 0.0021*PND + BW0) and all",
        "dam tissue volumes and blood flows follow postnatal-day growth",
        "equations. The nursing pups grow by their own generalised",
        "Michaelis-Menten body-weight curve (Mirfazaelian and Fisher 2007)",
        "independent of WT. The authors dose against a lactational",
        "reference weight of 0.242 kg when converting a mg/kg/day gavage",
        "dose to a mg amount; the vignette follows that convention."
      ),
      source_name = "BW0"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = NA_integer_,
    n_studies = 4L,
    age_range = "lactating dam and nursing litter, postnatal day 0 to 21",
    weight_range = "BW0 0.247 kg; lactational dosing reference 0.242 kg",
    sex_female_pct = 100,
    disease_state = paste(
      "Lactating Sprague-Dawley dam and nursing pups. Calibrated against",
      "lactational rat toxicokinetic data (Chang et al. 2009; Luebker et",
      "al. 2005a, 2005b)."
    ),
    dose_range = "Oral gavage 0.1-2 mg/kg per day during lactation (study-specific).",
    regions = NA_character_,
    notes = paste(
      "Calibrated (as-run) parameters are the sensitivity-selected subset",
      "in LFit_R.rds: dam Free, PL, PRest, milk suckling KMilk0, apical",
      "Vmax_apical_invitro and pup PL_pup. All other values are the code's",
      "$PARAM defaults. The litter is one pup compartment scaled by N = 8.",
      "Pup urine is returned to the dam gut (coprophagy), as in the code."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific and milk-transfer parameters for PFOS in the
    # lactating rat and nursing litter. As-run values: the code's LRatPBPK
    # $PARAM defaults with the sensitivity-selected subset overridden by
    # LFit_R.rds (calibrated, log-scale). Source of record: authors'
    # published code (https://github.com/KSUICCM/PFOS-Ges-Lac,
    # ModFit/Rat/RMod.R, LFit_R.rds). Every value is FIXED; the model is
    # deterministic.
    # ---------------------------------------------------------------------
    lfu <- fixed(log(0.00373148)); label("Free fraction of PFOS in dam plasma (unitless)")  # LFit_R calibrated 'Free' = 0.00373148 (Table 2 rat lactating 0.0037)

    lk0     <- fixed(log(1.000));   label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # RMod LRatPBPK 'K0C' = 1
    lkabs   <- fixed(log(2.12));    label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # RMod LRatPBPK 'KabsC' = 2.12
    lkunabs <- fixed(log(7.05e-5)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # RMod LRatPBPK 'KunabsC' = 7.05e-5
    lge     <- fixed(log(1.4));     label("Gastric emptying rate constant GEC (1/h per BW^0.25)")              # RMod LRatPBPK 'GEC' = 1.4

    lkp_liver   <- fixed(log(2.72303));   label("Dam liver-to-plasma partition coefficient PL (unitless)")           # LFit_R calibrated 'PL' = 2.72303 (Table 2 rat lactating 2.72)
    lkp_kidney  <- fixed(log(0.80));      label("Dam kidney-to-plasma partition coefficient PK (unitless)")          # RMod LRatPBPK 'PK' = 0.80
    lkp_fat     <- fixed(log(0.13));      label("Dam fat-to-plasma partition coefficient PF (unitless)")             # RMod LRatPBPK 'PF' = 0.13
    lkp_mammary <- fixed(log(0.16));      label("Dam mammary-to-plasma partition coefficient PM (unitless)")         # RMod LRatPBPK 'PM' = 0.16
    lkp_rest    <- fixed(log(0.0327343)); label("Dam rest-of-body-to-plasma partition coefficient PRest (unitless)") # LFit_R calibrated 'PRest' = 0.0327343 (Table 2 rat lactating 0.03)

    lpmilk_m  <- fixed(log(1.9));  label("Milk-to-mammary-gland partition coefficient PMilkM (unitless)")   # RMod LRatPBPK 'PMilkM' = 1.9
    lpamilk   <- fixed(log(0.5));  label("Mammary-to-milk permeability-area product PAMilkC (L/h per kg)")  # RMod LRatPBPK 'PAMilkC' = 0.5
    lkmilk0   <- fixed(log(0.278959)); label("Milk suckling/production rate constant KMilk0 (L/h per kg)")  # LFit_R calibrated 'KMilk0' = 0.278959 (Table 2 rat lactating 0.28)

    lkbile  <- fixed(log(0.0026)); label("Dam biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # RMod LRatPBPK 'KbileC' = 0.0026
    lkurine <- fixed(log(1.60));   label("Dam urinary elimination rate constant KurineC (1/h per BW^-0.25)") # RMod LRatPBPK 'KurineC' = 1.60
    lgfr    <- fixed(log(41.04));  label("Glomerular filtration rate constant GFRC (L/h per kg kidney)")     # RMod LRatPBPK 'GFRC' = 41.04

    lvmax_baso_invitro   <- fixed(log(393.45));  label("Dam in vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)")  # RMod LRatPBPK 'Vmax_baso_invitro' = 393.45
    lvmax_apical_invitro <- fixed(log(4140.65)); label("Dam in vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")         # LFit_R calibrated 'Vmax_apical_invitro' = 4140.65 (Table 2 rat lactating 4141)
    lkm_baso             <- fixed(log(27.2));    label("Dam Michaelis constant, basolateral transporters (mg/L)")         # RMod LRatPBPK 'Km_baso' = 27.2
    lkm_apical           <- fixed(log(278));     label("Dam Michaelis constant, apical transporters (mg/L)")              # RMod LRatPBPK 'Km_apical' = 278
    lrafbaso             <- fixed(log(1.90));    label("Relative activity factor, basolateral transporters (unitless)")   # RMod LRatPBPK 'RAFbaso' = 1.90
    lrafapi              <- fixed(log(4.15));    label("Relative activity factor, apical transporters (unitless)")        # RMod LRatPBPK 'RAFapi' = 4.15
    lkdif                <- fixed(log(5.1e-4));  label("Dam diffusion rate, kidney blood to proximal tubule cells (L/h)") # RMod LRatPBPK 'Kdif' = 5.1e-4
    lkefflux             <- fixed(log(2.09));    label("Dam efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)")  # RMod LRatPBPK 'KeffluxC' = 2.09

    # ---- Nursing-pup chemical parameters ----
    lfu_pup       <- fixed(log(0.022));   label("Free fraction of PFOS in pup plasma (unitless)")            # RMod LRatPBPK 'Free_pup' = 0.022
    lkp_liver_pup <- fixed(log(2.54638)); label("Pup liver-to-plasma partition coefficient PL_pup (unitless)") # LFit_R calibrated 'PL_pup' = 2.54638 (Table 2 rat lactating 2.55)
    lkp_kidney_pup <- fixed(log(0.80));   label("Pup kidney-to-plasma partition coefficient PK_pup (unitless)") # RMod LRatPBPK 'PK_pup' = 0.80
    lkp_rest_pup  <- fixed(log(0.22));    label("Pup rest-of-body-to-plasma partition coefficient PRest_pup (unitless)") # RMod LRatPBPK 'PRest_pup' = 0.22
    lvmax_apical_invitro_pup <- fixed(log(1808)); label("Pup in vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")  # RMod LRatPBPK 'Vmax_apical_invitro_p' = 1808 (uncalibrated pup default)
    lkdif_pup <- fixed(log(0.001)); label("Pup diffusion rate, kidney blood to proximal tubule cells (L/h)")  # RMod LRatPBPK 'Kdif_pup_0' = 0.001

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
    KMilk0   <- exp(lkmilk0)
    KbileC   <- exp(lkbile)
    KurineC  <- exp(lkurine)
    GFRC     <- exp(lgfr)
    Vmax_baso_invitro   <- exp(lvmax_baso_invitro)
    Vmax_apical_invitro <- exp(lvmax_apical_invitro)
    Km_baso   <- exp(lkm_baso)
    Km_apical <- exp(lkm_apical)
    RAFbaso   <- exp(lrafbaso)
    RAFapi    <- exp(lrafapi)
    Kdif      <- exp(lkdif)
    KeffluxC  <- exp(lkefflux)
    fu_pup    <- exp(lfu_pup)
    PL_pup    <- exp(lkp_liver_pup)
    PK_pup    <- exp(lkp_kidney_pup)
    PRest_pup <- exp(lkp_rest_pup)
    Vmax_apical_invitro_p <- exp(lvmax_apical_invitro_pup)
    Kdif_pup  <- exp(lkdif_pup)

    # =====================================================================
    # 1. Fixed physiology and growth-curve constants (RMod LRatPBPK $PARAM)
    # =====================================================================
    BW0     <- WT        # dam PND0 weight (default 0.247 kg)
    Htc     <- 0.46      # RMod LRatPBPK Htc
    QLC0    <- 0.196     # RMod LRatPBPK dam QLC0
    QKC0    <- 0.1161    # RMod LRatPBPK dam QKC0
    QMC0    <- 0.0887    # RMod LRatPBPK dam QMC0
    QFC     <- 0.07      # RMod LRatPBPK dam QFC
    QCl0    <- 28.645    # RMod LRatPBPK dam QCl0 (L/h/kg)
    VMilk   <- 0.002     # RMod LRatPBPK milk volume (L)
    VMC0    <- 0.049     # RMod LRatPBPK dam VMC0
    VLC0    <- 0.046     # RMod LRatPBPK dam VLC0
    VKC0    <- 0.011     # RMod LRatPBPK dam VKC0
    VPlasC  <- 0.0466    # RMod LRatPBPK dam VPlasC
    VFilC   <- 8.4e-4    # RMod LRatPBPK dam VFilC (L/kg BW)
    VPTCC   <- 1.35e-4   # RMod LRatPBPK dam VPTCC (L/kg kidney)
    protein <- 2.0e-6    # RMod LRatPBPK protein (mg protein/PTC)
    MW      <- 500.126   # RMod LRatPBPK MW (g/mol)
    N       <- 8         # RMod LRatPBPK number of pups (litter)
    KabsC_pup  <- 2.12   # RMod LRatPBPK pup KabsC_pup
    KbileC_pup <- 0.0026 # RMod LRatPBPK pup KbileC_pup
    KurineC_pup <- 1.6   # RMod LRatPBPK pup KurineC_pup
    KeffluxC_p <- 2.09   # RMod LRatPBPK pup KeffluxC_p
    Vmax_baso_invitro_p <- 393.45 # RMod LRatPBPK pup Vmax_baso_invitro_p
    Km_baso_p  <- 27.2   # RMod LRatPBPK pup Km_baso_p
    Km_apical_p <- 278   # RMod LRatPBPK pup Km_apical_p
    # Pup body-weight / liver / kidney generalised Michaelis-Menten growth
    Wt0     <- 0.00731   # RMod LRatPBPK pup Wt0 (kg, one pup at PND0)
    Kg      <- 63.21     # RMod LRatPBPK pup K (half-maximal, days)
    gg      <- 2.01      # RMod LRatPBPK pup g (Hill)
    Wtmax   <- 0.52      # RMod LRatPBPK pup Wtmax (kg)
    Wt_LIV0 <- 2.96e-4   # RMod LRatPBPK pup Wt_LIV0
    K_LIV   <- 43.49     # RMod LRatPBPK pup K_LIV
    g_LIV   <- 2.76      # RMod LRatPBPK pup g_LIV
    Wtmax_LIV <- 1.54e-2 # RMod LRatPBPK pup Wtmax_LIV
    Wt_KID0 <- 6.01e-5   # RMod LRatPBPK pup Wt_KID0
    K_KID   <- 50.83     # RMod LRatPBPK pup K_KID
    g_KID   <- 1.78      # RMod LRatPBPK pup g_KID
    Wtmax_KID <- 3.8e-3  # RMod LRatPBPK pup Wtmax_KID
    QCmax   <- 8.72      # RMod LRatPBPK pup QCmax (L/h)
    BW50    <- 0.189     # RMod LRatPBPK pup BW50 (kg)

    # =====================================================================
    # 2. Postnatal time and dam growth equations (Loccisano et al. 2012)
    # =====================================================================
    PND <- t / 24

    if (PND > 0) {
      BW  <- 0.0021 * PND + BW0
      QCl <- 0.0123 * PND^3 - 0.4059 * PND^2 + 3.9661 * PND + QCl0
      QMC <- (-0.0001) * PND^2 + 0.0051 * PND + QMC0
      QLC <- 0.0079 * PND + 0.2208
      QKC <- (-0.0014) * PND + QKC0
      VMC <- (-1e-5) * PND^3 + 0.0004 * PND^2 + 0.0027 * PND + 0.049
      VLC <- 0.0384 * PND^(0.096)
      VKC <- 0.0685 * PND^(0.0514)
      KMilkC <- (-7e-6) * PND^3 + 0.0003 * PND^2 - 0.0032 * PND + KMilk0
    } else {
      BW  <- BW0
      QCl <- QCl0
      QMC <- QMC0
      QLC <- QLC0
      QKC <- QKC0
      VMC <- VMC0
      VLC <- VLC0
      VKC <- VKC0
      KMilkC <- KMilk0
    }
    if (PND > 16) {
      VFC <- 0.07
    } else {
      VFC <- (-0.0012) * PND^2 + 0.0162 * PND + 0.1245
    }

    # ---- Dam blood flows and volumes ----
    QC1   <- QCl * BW
    QC    <- QC1 * (1 - Htc)
    QL    <- QLC * QC
    QM    <- QMC * QC
    QK    <- QKC * QC
    QF    <- QFC * QC
    QRest <- QC - (QL + QK + QM + QF)

    VPlas <- VPlasC * BW
    VM    <- VMC * BW
    VL    <- VLC * BW
    VK    <- VKC * BW
    VF    <- VFC * BW
    MK    <- VKC * BW * 1000
    VPTC  <- MK * VPTCC
    VKb   <- VK * 0.16
    VFil  <- VFilC * BW
    VRest <- 0.93 * BW - (VL + VK + VF + VM + VPlas)

    # ---- Dam scaled kinetic parameters ----
    GFR    <- GFRC * (VKC / 1000)
    GE     <- GEC * BW^(-0.25)
    K0     <- K0C * BW^(-0.25)
    Kbile  <- KbileC * BW^(-0.25)
    Kurine <- KurineC * BW^(-0.25)
    Kabs   <- KabsC * BW^(-0.25)
    Kunabs <- KunabsC * BW^(-0.25)
    PTC          <- VKC * 6e7 * 1000
    Vmax_basoC   <- Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC <- Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso    <- Vmax_basoC * BW^0.75
    Vmax_apical  <- Vmax_apicalC * BW^0.75
    Kefflux      <- KeffluxC * BW^(-0.25)

    # =====================================================================
    # 3. Pup growth (Mirfazaelian and Fisher 2007)
    # =====================================================================
    BW_pup   <- (Wt0 * Kg^gg + Wtmax * PND^gg) / (Kg^gg + PND^gg)
    BW_pup_W <- BW_pup * N
    QC_pup_i <- (QCmax * BW_pup) / (BW50 + BW_pup)
    QC_pup_W <- QC_pup_i * N
    QLC_pup  <- (4e-6) * PND^3 - 0.0005 * PND^2 + 0.0175 * PND + 0.0571
    QKC_pup  <- (-4e-5) * PND^2 + 0.0033 * PND + 0.0126
    Htc_pup  <- (-0.0136) * PND + 0.3814
    VPlasC_pup <- (-6e-6) * PND^3 + (1e-4) * PND^2 + 0.0006 * PND + 0.0454
    KMilk    <- KMilkC * BW_pup_W
    PAMilk   <- PAMilkC * (BW_pup^0.75) * N
    VL_pup   <- (Wt_LIV0 * K_LIV^g_LIV + Wtmax_LIV * PND^g_LIV) / (K_LIV^g_LIV + PND^g_LIV)
    VL_pup_W <- VL_pup * N
    VK_pup   <- (Wt_KID0 * K_KID^g_KID + Wtmax_KID * PND^g_KID) / (K_KID^g_KID + PND^g_KID)
    VK_pup_W <- VK_pup * N
    MK_pup   <- VKC * BW_pup * 1000
    VPTC_pup_W <- MK_pup * VPTCC * N
    VKb_pup_W  <- VK_pup_W * 0.16 * N
    VFil_pup_W <- VFilC * BW_pup_W
    VPlas_pup_W <- VPlasC_pup * BW_pup_W
    VRest_pup_W <- 0.92 * BW_pup_W - (VL_pup_W + VPlas_pup_W + VK_pup_W + VPTC_pup_W + VFil_pup_W)

    GFR_pup_W    <- GFRC * (MK_pup / 1000) * N
    Kbile_pup_W  <- KbileC_pup * BW_pup^(-0.25) * N
    Kurine_pup_W <- KurineC_pup * BW_pup^(-0.25) * N
    Kabs_pup_W   <- KabsC_pup * BW^(-0.25) * N

    QPlas_pup <- (1 - Htc_pup)
    QC_pup    <- QC_pup_W * QPlas_pup
    QL_pup    <- QC_pup * QLC_pup
    QK_pup    <- QC_pup * QKC_pup
    QRest_pup <- QC_pup - (QL_pup + QK_pup)

    Vmax_basoC_p   <- Vmax_baso_invitro_p * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC_p <- Vmax_apical_invitro_p * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso_pup   <- Vmax_basoC_p * BW_pup^0.75
    Vmax_apical_pup <- Vmax_apicalC_p * BW_pup^0.75
    Kefflux_pup <- KeffluxC_p * BW_pup^(-0.25)
    Km_baso_pup   <- Km_baso_p
    Km_apical_pup <- Km_apical_p

    # =====================================================================
    # 4. Concentrations (mg/L)
    # =====================================================================
    CPlas_free <- plasma / VPlas
    CPlas      <- CPlas_free / fu
    CKb        <- kidney_blood / VKb
    CVK        <- CKb
    CKid       <- CKb * PK
    CPTC       <- ptc / VPTC
    CFil       <- filtrate / (VFil + 1e-7)
    CMilk      <- milk / VMilk
    CL         <- liver / VL
    CVL        <- CL / PL
    CM         <- mammary / VM
    CVM        <- CM / PM
    CF         <- fat / VF
    CVF        <- CF / PF
    CRest      <- rest / VRest
    CVRest     <- CRest / PRest

    CPlas_free_pup <- plasma_pup / (VPlas_pup_W + 1e-7)
    CPlas_pup      <- CPlas_free_pup / fu_pup
    CL_pup         <- liver_pup / (VL_pup_W + 1e-7)
    CVL_pup        <- CL_pup / PL_pup
    CKb_pup        <- kidney_blood_pup / (VKb_pup_W + 1e-7)
    CVK_pup        <- CKb_pup
    CPTC_pup       <- ptc_pup / (VPTC_pup_W + 1e-7)
    CFil_pup       <- filtrate_pup / (VFil_pup_W + 1e-7)
    CRest_pup      <- rest_pup / (VRest_pup_W + 1e-7)
    CVRest_pup     <- CRest_pup / PRest_pup

    # =====================================================================
    # 5. Fluxes (mg/h)
    # =====================================================================
    RA_baso   <- Vmax_baso * (kidney_blood / VKb) / (Km_baso + kidney_blood / VKb)
    RA_apical <- Vmax_apical * (filtrate / (VFil + 1e-7)) / (Km_apical + filtrate / (VFil + 1e-7))
    Rdif      <- Kdif * (CKb - CPTC)
    RAefflux  <- Kefflux * ptc
    RCI       <- CPlas * GFR * fu
    Rtrans    <- KMilk * CMilk

    RAbaso_pup   <- Vmax_baso_pup * (kidney_blood_pup / (VKb_pup_W + 1e-7)) / (Km_baso_pup + kidney_blood_pup / (VKb_pup_W + 1e-7))
    RAapical_pup <- Vmax_apical_pup * (filtrate_pup / (VFil_pup_W + 1e-7)) / (Km_apical_pup + filtrate_pup / (VFil_pup_W + 1e-7))
    Rdif_pup     <- Kdif_pup * (CKb_pup - CPTC_pup)
    RAefflux_pup <- Kefflux_pup * ptc_pup
    RCI_pup      <- CPlas_pup * GFR_pup_W * fu_pup

    # =====================================================================
    # 6. Dam ODEs (dam gut receives pup urine via coprophagy)
    # =====================================================================
    d/dt(stomach)   <- -K0 * stomach - GE * stomach + filtrate_pup * Kurine_pup_W
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
    # 7. Pup ODEs (gut fed by suckled milk Rtrans)
    # =====================================================================
    d/dt(gi_pup)           <- Rtrans - Kabs_pup_W * gi_pup
    d/dt(liver_pup)        <- QL_pup * (CPlas_pup - CVL_pup) * fu_pup + Kabs_pup_W * gi_pup - Kbile_pup_W * liver_pup * fu_pup
    d/dt(rest_pup)         <- QRest_pup * (CPlas_pup - CVRest_pup) * fu_pup
    d/dt(kidney_blood_pup) <- QK_pup * (CPlas_pup - CVK_pup) * fu_pup - RCI_pup - Rdif_pup - RAbaso_pup
    d/dt(ptc_pup)          <- Rdif_pup + RAapical_pup + RAbaso_pup - RAefflux_pup
    d/dt(filtrate_pup)     <- RCI_pup - RAapical_pup - filtrate_pup * Kurine_pup_W
    d/dt(plasma_pup)       <- QRest_pup * CVRest_pup * fu_pup + QL_pup * CVL_pup * fu_pup +
                              QK_pup * CVK_pup * fu_pup - QC_pup * CPlas_pup * fu_pup + RAefflux_pup
    d/dt(urine_pup)        <- filtrate_pup * Kurine_pup_W
    d/dt(feces_pup)        <- Kbile_pup_W * liver_pup * fu_pup

    # =====================================================================
    # 8. Exposure accumulators
    # =====================================================================
    d/dt(auc_plasma)     <- CPlas
    d/dt(auc_plasma_pup) <- CPlas_pup

    # =====================================================================
    # 9. Observations
    # =====================================================================
    Cc        <- CPlas
    Cliver    <- CL
    Ckidney   <- CKid
    Cmilk     <- CMilk
    Cc_pup    <- CPlas_pup
    Cliver_pup <- CL_pup

    Cc ~ prop(propSd)
  })
}
