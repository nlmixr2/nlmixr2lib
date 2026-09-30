Chou_2021_pfos_rat_gestational_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, gestational; original implementation in",
    "R/mrgsolve). Perfluorooctane sulfonate (PFOS) in the pregnant",
    "Sprague-Dawley rat and her litter (Chou and Lin 2021, Environ Health",
    "Perspect). Extends the pre-pregnancy maternal model with placenta and",
    "fetal compartments and time-varying gestational physiology (growth",
    "of fat, mammary gland, placenta, plasma, amniotic fluid and the",
    "fetus, and the corresponding blood-flow changes, all functions of",
    "gestational day). The maternal side is a two-compartment",
    "gastrointestinal tract feeding the liver, plus plasma, liver, fat,",
    "mammary gland, rest of body, placenta and a three-subcompartment",
    "kidney (kidney blood, proximal tubule cells, filtrate) with",
    "glomerular filtration, saturable basolateral (Oat1/Oat3) and apical",
    "(Oatp1a1) reabsorption, passive diffusion and first-order efflux.",
    "The fetal submodel (whole litter of 8) is plasma, liver, rest of",
    "body and amniotic fluid, connected to the dam by bidirectional",
    "placental and amniotic-fluid diffusion. Elimination is maternal",
    "urinary and faecal. PFOS is not metabolised; only the unbound",
    "fraction exchanges. Deterministic, no random effects. Parameter",
    "values are the calibrated (as-run) set from the authors' openly",
    "published mrgsolve code and fit object",
    "(https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Rat/RMod.R object",
    "GRatPBPK and GFit_R.rds)."
  )
  reference <- paste(
    "Chou WC, Lin Z. Development of a Gestational and Lactational",
    "Physiologically Based Pharmacokinetic (PBPK) Model for",
    "Perfluorooctane Sulfonate (PFOS) in Rats and Humans and Its",
    "Implications in the Derivation of Health-Based Toxicity Values.",
    "Environ Health Perspect. 2021;129(3):037004. doi:10.1289/EHP7671.",
    "Structure and gestational growth equations from the Methods and the",
    "authors' published code (ModFit/Rat/RMod.R, object GRatPBPK);",
    "calibrated parameters from GFit_R.rds (as-run, sensitivity-selected",
    "subset) and Table 2 (rat, pregnant column)",
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
    "placenta",
    "plasma_fet",
    "liver_fet",
    "rest_fet",
    "amniotic",
    "feces",
    "auc_plasma",
    "auc_plasma_fet"
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
    placenta = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    plasma_fet = list(analyte = "PFOS", units = "mg", specimen = "plasma", verified = TRUE),
    liver_fet = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    rest_fet = list(analyte = "PFOS", units = "mg", specimen = "tissue", verified = TRUE),
    amniotic = list(analyte = "PFOS", units = "mg", specimen = "not applicable", verified = TRUE),
    urine = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    feces = list(analyte = "PFOS", units = "mg", specimen = "faeces", verified = TRUE),
    auc_plasma = list(analyte = "PFOS", units = "mg/L*h", specimen = "not applicable", verified = TRUE),
    auc_plasma_fet = list(analyte = "PFOS", units = "mg/L*h", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Pre-pregnancy maternal body weight (BW0) that drives the",
        "gestational growth equations. The authors' code uses BW0 = 0.185",
        "kg. Gestational maternal weight BW_P is computed internally as",
        "BW0 plus the growing fat, mammary, placenta, fetal and amniotic",
        "volumes. Note the authors dose against a gestational reference",
        "weight of 0.225 kg (Loccisano et al. 2012) when converting a",
        "mg/kg/day gavage dose to a mg amount; the vignette follows that",
        "convention."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = NA_integer_,
    n_studies = 4L,
    age_range = "pregnant dam, gestational day 0 to 21 (litter of 8)",
    weight_range = "BW0 0.185 kg; gestational dosing reference 0.225 kg",
    sex_female_pct = 100,
    disease_state = paste(
      "Pregnant Sprague-Dawley rat. Calibrated against gestational rat",
      "toxicokinetic data (Thibodeaux et al. 2003; Chang et al. 2009;",
      "Luebker et al. 2005a, 2005b)."
    ),
    dose_range = "Oral gavage 0.1-10 mg/kg per day during gestation (study-specific).",
    regions = NA_character_,
    notes = paste(
      "Calibrated (as-run) parameters are the sensitivity-selected subset",
      "in GFit_R.rds: Free, PL, PL_Fet, PRest_Fet, KbileC,",
      "Vmax_baso_invitro, Km_apical, Ktrans1C and Ktrans3C. All other",
      "values are the code's $PARAM defaults. The litter is treated as one",
      "fetal compartment (N = 8). See the vignette Errata for the Table 2",
      "vs code Km_baso/Km_apical labelling note."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific and placental-transfer parameters for PFOS in the
    # pregnant rat. As-run values: the code's GRatPBPK $PARAM defaults with
    # the sensitivity-selected subset overridden by GFit_R.rds (calibrated,
    # log-scale). Source of record: authors' published code
    # (https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Rat/RMod.R,
    # GFit_R.rds). Every value is FIXED; the model is deterministic.
    # ---------------------------------------------------------------------
    lfu <- fixed(log(0.0190056)); label("Free fraction of PFOS in maternal plasma (unitless)")  # GFit_R calibrated 'Free' = 0.0190056 (Table 2 rat pregnant 0.019)

    lk0     <- fixed(log(1.000));   label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # RMod GRatPBPK 'K0C' = 1
    lkabs   <- fixed(log(2.12));    label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # RMod GRatPBPK 'KabsC' = 2.12
    lkunabs <- fixed(log(7.05e-5)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # RMod GRatPBPK 'KunabsC' = 7.05e-5
    lge     <- fixed(log(1.4));     label("Gastric emptying rate constant GEC (1/h per BW^0.25)")              # RMod GRatPBPK 'GEC' = 1.4

    lkp_liver    <- fixed(log(3.16684)); label("Liver-to-plasma partition coefficient PL (unitless)")            # GFit_R calibrated 'PL' = 3.16684 (Table 2 rat pregnant 3.17)
    lkp_kidney   <- fixed(log(0.80));    label("Kidney-to-plasma partition coefficient PK (unitless)")           # RMod GRatPBPK 'PK' = 0.80
    lkp_fat      <- fixed(log(0.13));    label("Fat-to-plasma partition coefficient PF (unitless)")              # RMod GRatPBPK 'PF' = 0.13
    lkp_mammary  <- fixed(log(0.16));    label("Mammary-to-plasma partition coefficient PM (unitless)")          # RMod GRatPBPK 'PM' = 0.16
    lkp_rest     <- fixed(log(0.22));    label("Rest-of-body-to-plasma partition coefficient PRest (unitless)")  # RMod GRatPBPK 'PRest' pregnant = 0.22
    lkp_placenta <- fixed(log(0.41));    label("Placenta-to-plasma partition coefficient PPla (unitless)")       # RMod GRatPBPK 'PPla' = 0.41

    lkbile  <- fixed(log(0.00683402)); label("Biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # GFit_R calibrated 'KbileC' = 0.00683402 (Table 2 rat pregnant 0.007)
    lkurine <- fixed(log(1.60));       label("Urinary elimination rate constant KurineC (1/h per BW^-0.25)") # RMod GRatPBPK 'KurineC' = 1.60
    lgfr    <- fixed(log(41.04));      label("Glomerular filtration rate constant GFRC (L/h per kg kidney)")  # RMod GRatPBPK 'GFRC' = 41.04

    lvmax_baso_invitro   <- fixed(log(220.535)); label("In vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)")  # GFit_R calibrated 'Vmax_baso_invitro' = 220.535 (Table 2 rat pregnant 221)
    lvmax_apical_invitro <- fixed(log(1808));    label("In vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")         # RMod GRatPBPK 'Vmax_apical_invitro' = 1808
    lkm_baso             <- fixed(log(27.2));    label("Michaelis constant, basolateral transporters (mg/L)")         # RMod GRatPBPK 'Km_baso' = 27.2 (code default; not calibrated)
    lkm_apical           <- fixed(log(19.8859)); label("Michaelis constant, apical transporters (mg/L)")              # GFit_R calibrated code parameter 'Km_apical' = 19.8859 (see Errata re Table 2 labelling)
    lrafbaso             <- fixed(log(1.90));    label("Relative activity factor, basolateral transporters (unitless)") # RMod GRatPBPK 'RAFbaso' = 1.90
    lrafapi              <- fixed(log(4.15));    label("Relative activity factor, apical transporters (unitless)")      # RMod GRatPBPK 'RAFapi' = 4.15
    lkdif                <- fixed(log(5.1e-4));  label("Diffusion rate, kidney blood to proximal tubule cells (L/h)")   # RMod GRatPBPK 'Kdif' = 5.1e-4
    lkefflux             <- fixed(log(2.09));    label("Efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)")  # RMod GRatPBPK 'KeffluxC' = 2.09

    lktrans1 <- fixed(log(1.26751)); label("Mother-to-fetus placental transfer rate Ktrans1C (L/h per kg^0.75)")  # GFit_R calibrated 'Ktrans1C' = 1.26751 (Table 2 rat pregnant 1.27)
    lktrans2 <- fixed(log(1.00));    label("Fetus-to-mother placental transfer rate Ktrans2C (L/h per kg^0.75)")  # RMod GRatPBPK 'Ktrans2C' = 1
    lktrans3 <- fixed(log(0.229298)); label("Fetus-to-amniotic-fluid transfer rate Ktrans3C (L/h per kg^0.75)")  # GFit_R calibrated 'Ktrans3C' = 0.229298 (Table 2 rat pregnant 0.23)
    lktrans4 <- fixed(log(0.001));   label("Amniotic-fluid-to-fetus transfer rate Ktrans4C (L/h per kg^0.75)")   # RMod GRatPBPK 'Ktrans4C' = 0.001

    lfree_fet     <- fixed(log(0.022));    label("Free fraction of PFOS in fetal plasma (unitless)")             # RMod GRatPBPK 'Free_Fet' = 0.022
    lkp_liver_fet <- fixed(log(1.30378));  label("Fetal liver-to-plasma partition coefficient PL_Fet (unitless)") # GFit_R calibrated 'PL_Fet' = 1.30378 (Table 2 1.30)
    lkp_rest_fet  <- fixed(log(0.109815)); label("Fetal rest-of-body-to-plasma partition coefficient PRest_Fet (unitless)") # GFit_R calibrated 'PRest_Fet' = 0.109815 (Table 2 0.11)

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
    PPla     <- exp(lkp_placenta)
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
    Ktrans1C  <- exp(lktrans1)
    Ktrans2C  <- exp(lktrans2)
    Ktrans3C  <- exp(lktrans3)
    Ktrans4C  <- exp(lktrans4)
    fu_fet    <- exp(lfree_fet)
    PL_Fet    <- exp(lkp_liver_fet)
    PRest_Fet <- exp(lkp_rest_fet)

    # =====================================================================
    # 1. Fixed physiology (RMod GRatPBPK $PARAM). N = whole litter.
    # =====================================================================
    Htc     <- 0.46      # RMod rat Htc
    QLC     <- 0.183     # RMod rat QLC
    QKC     <- 0.141     # RMod rat QKC
    QMC     <- 0.002     # RMod rat QMC
    QFC     <- 0.07      # RMod rat QFC
    VLC     <- 0.035     # RMod rat VLC
    VKC     <- 0.0084    # RMod rat VKC
    VMC     <- 0.01      # RMod rat VMC
    VFC     <- 0.07      # RMod rat VFC
    VPlasC  <- 0.0466    # RMod rat VPlasC
    VFilC   <- 8.4e-4    # RMod rat VFilC (L/kg BW)
    VPTCC   <- 1.35e-4   # RMod rat VPTCC (L/kg kidney)
    protein <- 2.0e-6    # RMod rat protein (mg protein/PTC)
    MW      <- 500.126   # RMod rat MW (g/mol)
    N          <- 8      # RMod GRatPBPK number of fetuses (litter)
    VPlasC_Fet <- 0.047  # RMod GRatPBPK fetal plasma fraction
    QLC_Fet    <- 0.061  # RMod GRatPBPK fetal liver blood-flow fraction

    # =====================================================================
    # 2. Gestational time (BW0 = WT drives growth). GD in days, GA in weeks
    # =====================================================================
    GD <- t / 24
    GA <- t / 168

    # ---- Placenta volume (L), Loccisano et al. 2012 growth/decay ----
    if (GD <= 6) {
      VPla <- 1e-7
    } else if (GD <= 10) {
      VPla <- (N * (8 * (GD - 6))) / 1.0e6
    } else {
      VPla <- ((N * (32 * exp(-0.23 * (GD - 10)))) + (40 * (exp(0.28 * (GD - 10)) - 1))) / 1.0e6
    }

    # ---- Plasma flow to placenta (L/h) ----
    Qdec1 <- 0
    Qdec2 <- 0
    Qcap  <- 0
    if (GD > 6 && GD <= 10) {
      Qdec1 <- 0.55 * (GD - 6)
    }
    if (GD > 10 && GD <= 12) {
      Qdec2 <- 2.2 * exp(-0.23 * (GD - 10))
    }
    if (GD > 12) {
      Qcap <- (0.1207 * (GD - 12))^4.36
    }
    Qdec  <- Qdec1 + Qdec2
    QPla1 <- (N * (0.02 * Qdec + Qcap)) / 24
    QPla  <- QPla1 * (1 - Htc)

    # ---- Maternal tissue volumes (L); VM_P/VF_P expand with GD ----
    VM   <- VMC * WT
    if (GD > 3) {
      VM_P <- VM * (1 + 0.2 * GD)
    } else {
      VM_P <- VM
    }
    VF   <- VFC * WT
    VF_P <- VF * (1 + 0.0182 * GD)
    VL   <- VLC * WT
    VK   <- VKC * WT
    VPlas <- VPlasC * WT
    MK   <- VKC * WT * 1000
    VPTC <- MK * VPTCC
    VKb  <- VK * 0.16
    VFil <- VFilC * WT

    # ---- Fetal growth (whole litter) ----
    VFet_1 <- (0.1089 + (16 * exp(-exp(5.515 - 0.2565 * GD)))) / 1000
    VFet   <- VFet_1 * N
    if (GD >= 10) {
      VAmX <- (-4e-6) * GD^3 + 0.0002 * GD^2 - 0.0023 * GD + 0.0099
    } else {
      VAmX <- 1e-7
    }
    VAm    <- VAmX * N
    VPlas_Fet <- VPlasC_Fet * VFet
    if (GD < 17) {
      VLC_Fet <- 0
    } else {
      VLC_Fet <- (-0.0013) * GD^3 + 0.0731 * GD^2 - 1.375 * GD + 8.6997
    }
    VL_Fet    <- VLC_Fet * VFet
    VRest_Fet <- 0.93 * VFet - VPlas_Fet - VL_Fet

    # ---- Fetal blood flows ----
    QC_Fet    <- QPla / (1 + 20000 * exp(-0.55 * GD))
    QL_Fet    <- QC_Fet * QLC_Fet
    QRest_Fet <- QC_Fet - QL_Fet

    # ---- Maternal weight and rest-of-body volume during pregnancy ----
    BW_P  <- WT + (VF_P - VF) + (VM_P - VM) + VPla + VFet + VAm
    VRest <- 0.93 * BW_P - (VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm)

    # ---- Maternal cardiac output and plasma flows ----
    QCl_P <- 24.56 - 0.1323 * GD
    QC_P  <- QCl_P * BW_P * (1 - Htc)
    QCl   <- 24.56
    QC    <- QCl * WT * (1 - Htc)
    QF    <- QFC * QC
    QF_P  <- QF * (VF_P / VF)
    QM    <- QMC * QC
    QM_P  <- QM * (VM_P / VM)
    QL    <- QLC * QC
    QK    <- QKC * QC
    QK_P  <- (-0.0014) * GD^2 + 0.0308 * GD + 0.449
    QRest <- QC_P - (QK_P + QL + QM_P + QF_P + QPla)

    # =====================================================================
    # 3. Scaled kinetic parameters (allometry on BW_P)
    # =====================================================================
    GFR     <- GFRC * (MK / 1000)
    GE      <- GEC * BW_P^(-0.25)
    K0      <- K0C * BW_P^(-0.25)
    Kbile   <- KbileC * BW_P^(-0.25)
    Kurine  <- KurineC * BW_P^(-0.25)
    Kabs    <- KabsC * BW_P^(-0.25)
    Kunabs  <- KunabsC * BW_P^(-0.25)
    PTC          <- VKC * 6e7 * 1000
    Vmax_basoC   <- Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC <- Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso    <- Vmax_basoC * BW_P^0.75
    Vmax_apical  <- Vmax_apicalC * BW_P^0.75
    Kefflux      <- KeffluxC * BW_P^(-0.25)
    Ktrans_1 <- Ktrans1C * (VFet_1^0.75 * N)
    Ktrans_2 <- Ktrans2C * (VFet_1^0.75 * N)
    Ktrans_3 <- Ktrans3C * (VFet_1^0.75 * N)
    Ktrans_4 <- Ktrans4C * (VFet_1^0.75 * N)

    # =====================================================================
    # 4. Concentrations (mg/L). Maternal fat/mammary use the GD0 volume, as
    #    in the authors' code (CM = AM/VM, CF = AF/VF).
    # =====================================================================
    CPlas_free <- plasma / VPlas
    CPlas      <- CPlas_free / fu
    CL         <- liver / VL
    CVL        <- CL / PL
    CKb        <- kidney_blood / VKb
    CVK        <- CKb
    CKid       <- CKb * PK
    CM         <- mammary / VM
    CVM        <- CM / PM
    CF         <- fat / VF
    CVF        <- CF / PF
    CAm        <- amniotic / (VAm + 1e-7)
    CRest      <- rest / VRest
    CVRest     <- CRest / PRest
    CPTC       <- ptc / VPTC
    CFil       <- filtrate / VFil
    CPla       <- placenta / VPla
    CVPla      <- CPla / PPla

    CPlas_Fet_free <- plasma_fet / (VPlas_Fet + 1e-7)
    CPlas_Fet      <- CPlas_Fet_free / fu_fet
    CRest_Fet      <- rest_fet / (VRest_Fet + 1e-7)
    CVRest_Fet     <- CRest_Fet / PRest_Fet
    CL_Fet         <- liver_fet / (VL_Fet + 1e-7)
    CVL_Fet        <- CL_Fet / PL_Fet

    # =====================================================================
    # 5. Fluxes (mg/h)
    # =====================================================================
    RA_baso   <- Vmax_baso * (kidney_blood / VKb) / (Km_baso + kidney_blood / VKb)
    RA_apical <- Vmax_apical * (filtrate / VFil) / (Km_apical + filtrate / VFil)
    Rdif      <- Kdif * (CKb - CPTC)
    RAefflux  <- Kefflux * ptc
    RCI       <- CPlas * GFR * fu
    Rtrans_1  <- Ktrans_1 * CVPla * fu
    Rtrans_2  <- Ktrans_2 * CPlas_Fet * fu
    Rtrans_3  <- Ktrans_3 * CVRest_Fet * fu_fet
    Rtrans_4  <- Ktrans_4 * CAm

    # =====================================================================
    # 6. Maternal ODEs
    # =====================================================================
    d/dt(stomach)   <- -K0 * stomach - GE * stomach
    d/dt(intestine) <- GE * stomach - Kabs * intestine - Kunabs * intestine
    d/dt(liver)     <- QL * (CPlas - CVL) * fu - Kbile * liver + Kabs * intestine + K0 * stomach
    d/dt(fat)       <- QF_P * (CPlas - CVF) * fu
    d/dt(mammary)   <- QM_P * (CPlas - CVM) * fu
    d/dt(rest)      <- QRest * (CPlas - CVRest) * fu
    d/dt(kidney_blood) <- QK_P * (CPlas - CVK) * fu - RCI - Rdif - RA_baso
    d/dt(ptc)          <- Rdif + RA_apical + RA_baso - RAefflux
    d/dt(filtrate)     <- RCI - RA_apical - Kurine * filtrate
    d/dt(placenta)  <- QPla * (CPlas - CVPla) * fu + Rtrans_2 - Rtrans_1
    d/dt(plasma)    <- QRest * CVRest * fu + QK_P * CVK * fu + QL * CVL * fu +
                       QM_P * CVM * fu + QF_P * CVF * fu + QPla * CVPla * fu -
                       QC_P * CPlas * fu + RAefflux
    d/dt(urine)     <- Kurine * filtrate
    d/dt(feces)     <- Kbile * liver + Kunabs * intestine

    # =====================================================================
    # 7. Fetal ODEs
    # =====================================================================
    d/dt(liver_fet)  <- QL_Fet * (CPlas_Fet - CVL_Fet) * fu_fet
    d/dt(rest_fet)   <- QRest_Fet * (CPlas_Fet - CVRest_Fet) * fu_fet - Rtrans_3 + Rtrans_4
    d/dt(plasma_fet) <- QRest_Fet * CVRest_Fet * fu_fet + QL_Fet * CVL_Fet * fu_fet -
                        QC_Fet * CPlas_Fet * fu_fet + Rtrans_1 - Rtrans_2
    d/dt(amniotic)   <- Rtrans_3 - Rtrans_4

    # =====================================================================
    # 8. Exposure accumulators
    # =====================================================================
    d/dt(auc_plasma)     <- CPlas
    d/dt(auc_plasma_fet) <- CPlas_Fet

    # =====================================================================
    # 9. Observations
    # =====================================================================
    Cc         <- CPlas
    Cliver     <- CL
    Ckidney    <- CKid
    Cplacenta  <- CPla
    Cc_fet     <- CPlas_Fet
    Cliver_fet <- CL_Fet

    Cc ~ prop(propSd)
  })
}
