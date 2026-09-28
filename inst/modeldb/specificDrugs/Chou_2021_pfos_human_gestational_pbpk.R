Chou_2021_pfos_human_gestational_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, gestational; original implementation in",
    "R/mrgsolve). Perfluorooctane sulfonate (PFOS) in the pregnant woman",
    "and her fetus (Chou and Lin 2021, Environ Health Perspect). Extends",
    "the pre-pregnancy maternal model with placenta and fetal compartments",
    "and time-varying gestational physiology (plasma volume, haematocrit,",
    "glomerular filtration, fat, mammary gland, placenta, amniotic fluid",
    "and fetal growth, all functions of gestational age) taken from",
    "Kapraun et al. 2019 and Loccisano et al. 2013. The maternal side is a",
    "two-compartment gastrointestinal tract feeding the liver, plus",
    "plasma, liver, fat, mammary gland, rest of body, placenta and a",
    "three-subcompartment kidney (kidney blood, proximal tubule cells,",
    "filtrate) with glomerular filtration, saturable basolateral",
    "(Oat1/Oat3) and apical (Oatp1a1) reabsorption, passive diffusion and",
    "first-order efflux. The fetal submodel is plasma, liver, rest of body",
    "and amniotic fluid, connected to the mother by bidirectional",
    "placental and amniotic-fluid diffusion. Elimination is maternal",
    "urinary and faecal. PFOS is not metabolised; only the unbound",
    "fraction exchanges. Deterministic, no random effects. Parameter",
    "values are the calibrated (as-run) set from the authors' openly",
    "published mrgsolve code and fit object",
    "(https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Human/HMod.R object",
    "GHumanPBPK and GFit_H.rds)."
  )
  reference <- paste(
    "Chou WC, Lin Z. Development of a Gestational and Lactational",
    "Physiologically Based Pharmacokinetic (PBPK) Model for",
    "Perfluorooctane Sulfonate (PFOS) in Rats and Humans and Its",
    "Implications in the Derivation of Health-Based Toxicity Values.",
    "Environ Health Perspect. 2021;129(3):037004. doi:10.1289/EHP7671.",
    "Structure and gestational growth equations from the Methods (Kapraun",
    "et al. 2019; Loccisano et al. 2013) and the authors' published code",
    "(ModFit/Human/HMod.R, object GHumanPBPK); calibrated parameters from",
    "GFit_H.rds (as-run) and Table 2 (human, pregnant column)",
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
        "Pre-pregnancy maternal body weight (BW0) that anchors the",
        "gestational growth equations. The authors' code uses BW0 = 60 kg",
        "(Yoon et al. 2011). Gestational maternal weight BW_P is computed",
        "internally as BW0 plus the growing fat, mammary, placenta, fetal",
        "and amniotic volumes. Plasma volume, haematocrit and glomerular",
        "filtration follow their own gestational-age growth equations",
        "rather than scaling with WT."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 15L,
    age_range = "pregnant women, gestational age 0 to 40 weeks",
    weight_range = "BW0 60 kg (authors' code default)",
    sex_female_pct = 100,
    disease_state = paste(
      "Pregnant women. Calibrated against maternal-plasma, cord-blood and",
      "placental/fetal PFOS biomonitoring data (Inoue et al. 2004; Fei et",
      "al. 2007; Kato et al. 2014; Pan et al. 2017; Mamsen et al. 2019;",
      "Midasch et al. 2007) and evaluated against further cohorts."
    ),
    dose_range = paste(
      "Chronic dietary intake. Population-specific estimates 1.2-3.7",
      "ng/kg/day (Loccisano et al. 2013), applied to the stomach."
    ),
    regions = "worldwide human biomonitoring",
    notes = paste(
      "Calibrated (as-run) parameters are the sensitivity-selected subset",
      "in GFit_H.rds: Free, PPla, Free_Fet, PL_Fet, PRest_Fet, KeffluxC,",
      "Km_apical, Ktrans1C and Ktrans2C. All other values are the code's",
      "$PARAM defaults. A single fetus is modelled. See the vignette",
      "Errata for the Table 2 vs code Km labelling note."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific and placental-transfer parameters for PFOS in the
    # pregnant woman. As-run values: the code's GHumanPBPK $PARAM defaults
    # with the sensitivity-selected subset overridden by GFit_H.rds
    # (calibrated, log-scale). Source of record: authors' published code
    # (https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Human/HMod.R,
    # GFit_H.rds). Every value is FIXED; the model is deterministic.
    # ---------------------------------------------------------------------
    lfu <- fixed(log(0.00269267)); label("Free fraction of PFOS in maternal plasma (unitless)")  # GFit_H calibrated 'Free' = 0.00269267 (Table 2 human pregnant 0.0027)

    lk0     <- fixed(log(1.000));   label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # HMod GHumanPBPK 'K0C' = 1
    lkabs   <- fixed(log(2.120));   label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # HMod GHumanPBPK 'KabsC' = 2.12
    lkunabs <- fixed(log(7.05e-5)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # HMod GHumanPBPK 'KunabsC' = 7.05e-5
    lge     <- fixed(log(3.510));   label("Gastric emptying rate constant GEC (1/h per BW^0.25)")              # HMod GHumanPBPK 'GEC' = 3.510

    lkp_liver    <- fixed(log(2.03)); label("Liver-to-plasma partition coefficient PL (unitless)")            # HMod GHumanPBPK 'PL' = 2.03
    lkp_kidney   <- fixed(log(1.26)); label("Kidney-to-plasma partition coefficient PK (unitless)")           # HMod GHumanPBPK 'PK' = 1.26
    lkp_fat      <- fixed(log(0.13)); label("Fat-to-plasma partition coefficient PF (unitless)")              # HMod GHumanPBPK 'PF' = 0.13
    lkp_mammary  <- fixed(log(0.16)); label("Mammary-to-plasma partition coefficient PM (unitless)")          # HMod GHumanPBPK 'PM' = 0.16
    lkp_rest     <- fixed(log(0.20)); label("Rest-of-body-to-plasma partition coefficient PRest (unitless)")  # HMod GHumanPBPK 'PRest' = 0.20
    lkp_placenta <- fixed(log(0.131036)); label("Placenta-to-plasma partition coefficient PPla (unitless)")   # GFit_H calibrated 'PPla' = 0.131036 (Table 2 human pregnant 0.13)

    lkbile  <- fixed(log(1.3e-4)); label("Biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # HMod GHumanPBPK 'KbileC' = 1.3e-4
    lkurine <- fixed(log(0.096));  label("Urinary elimination rate constant KurineC (1/h per BW^-0.25)") # HMod GHumanPBPK 'KurineC' = 0.096

    lvmax_baso_invitro   <- fixed(log(479));      label("In vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)")  # HMod GHumanPBPK 'Vmax_baso_invitro' = 479
    lvmax_apical_invitro <- fixed(log(51803));    label("In vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")         # HMod GHumanPBPK 'Vmax_apical_invitro' = 51803
    lkm_baso             <- fixed(log(20.1));     label("Michaelis constant, basolateral transporters (mg/L)")         # HMod GHumanPBPK 'Km_baso' = 20.1 (code default; not calibrated)
    lkm_apical           <- fixed(log(247.817));  label("Michaelis constant, apical transporters (mg/L)")              # GFit_H calibrated code parameter 'Km_apical' = 247.817 (Table 2 human pregnant 248; see Errata re labelling)
    lrafbaso             <- fixed(log(1));        label("Relative activity factor, basolateral transporters (unitless)") # HMod GHumanPBPK 'RAFbaso' = 1
    lrafapi              <- fixed(log(0.001));    label("Relative activity factor, apical transporters (unitless)")      # HMod GHumanPBPK 'RAFapi' = 0.001
    lkdif                <- fixed(log(0.001));    label("Diffusion rate, kidney blood to proximal tubule cells (L/h)")   # HMod GHumanPBPK 'Kdif' = 0.001
    lkefflux             <- fixed(log(0.0147148)); label("Efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)")  # GFit_H calibrated 'KeffluxC' = 0.0147148 (Table 2 human pregnant 0.015)

    lktrans1 <- fixed(log(0.788627)); label("Mother-to-fetus placental transfer rate Ktrans1C (L/h per kg^0.75)")  # GFit_H calibrated 'Ktrans1C' = 0.788627 (Table 2 human pregnant 0.79)
    lktrans2 <- fixed(log(1.11518));  label("Fetus-to-mother placental transfer rate Ktrans2C (L/h per kg^0.75)")  # GFit_H calibrated 'Ktrans2C' = 1.11518 (Table 2 human pregnant 1.12)
    lktrans3 <- fixed(log(0.006));    label("Fetus-to-amniotic-fluid transfer rate Ktrans3C (L/h per kg^0.75)")   # HMod GHumanPBPK 'Ktrans3C' = 0.006
    lktrans4 <- fixed(log(0.001));    label("Amniotic-fluid-to-fetus transfer rate Ktrans4C (L/h per kg^0.75)")   # HMod GHumanPBPK 'Ktrans4C' = 0.001

    lfree_fet     <- fixed(log(0.00382031)); label("Free fraction of PFOS in fetal plasma (unitless)")            # GFit_H calibrated 'Free_Fet' = 0.00382031 (Table 2 human 0.0038)
    lkp_liver_fet <- fixed(log(0.579717));   label("Fetal liver-to-plasma partition coefficient PL_Fet (unitless)") # GFit_H calibrated 'PL_Fet' = 0.579717 (Table 2 0.58)
    lkp_rest_fet  <- fixed(log(2.31331));    label("Fetal rest-of-body-to-plasma partition coefficient PRest_Fet (unitless)") # GFit_H calibrated 'PRest_Fet' = 2.31331 (Table 2 2.3)

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
    # 1. Fixed physiology (HMod GHumanPBPK $PARAM). A single fetus.
    # =====================================================================
    QCC     <- 16.4      # HMod human QCC (L/h/kg^0.75)
    QLC     <- 0.25      # HMod human QLC
    QKC     <- 0.141     # HMod human QKC
    QMC     <- 0.027     # HMod human QMC
    QFC     <- 0.052     # HMod human QFC
    VLC     <- 0.026     # HMod human VLC
    VKC     <- 0.004     # HMod human VKC
    VMC     <- 0.0062    # HMod human VMC
    VFC     <- 0.214     # HMod human VFC
    VFilC   <- 8.4e-4    # HMod human VFilC (L/kg BW)
    VPTCC   <- 1.35e-4   # HMod human VPTCC (L/kg kidney)
    protein <- 2.0e-6    # HMod human protein (mg protein/PTC)
    MW      <- 500.126   # HMod human MW (g/mol)
    VPlasC_Fet <- 0.0428 # HMod GHumanPBPK fetal plasma fraction

    # =====================================================================
    # 2. Gestational time and growth equations (Kapraun et al. 2019)
    # =====================================================================
    GD <- t / 24
    GA <- GD / 7

    VPlas <- (1.2406 / (1 + exp(-0.31338 * (GA - 17.813)))) + 2.4958
    Htc   <- (39.192 - 0.10562 * GA - (7.1045e-4) * GA^2) / 100
    GFR   <- (113.73 + 3.5784 * GA - 0.067272 * GA^2) * 0.06
    VAm   <- (822.34 / (1 + exp(-0.26988 * (GA - 20.150)))) / 1000
    if (GA < 2) {
      VPla <- 0
    } else {
      VPla <- (-1.7646 * GA + 0.91775 * GA^2 - 0.011543 * GA^3) / 1000
    }
    VF_P <- (1 / 0.95) * (17.067 + 0.14937 * GA)
    VF_0 <- WT * VFC
    VM_0 <- WT * VMC
    VM_P <- WT * (VMC + (0.0065 * exp(-7.444868 * exp(-0.000678 * (GD * 24)))))
    VL   <- VLC * WT
    VK   <- VKC * WT
    MK   <- VKC * WT * 1000
    VKb  <- VK * 0.16
    VFil <- VFilC * WT
    VPTC <- MK * VPTCC

    # ---- Fetal growth (single fetus) ----
    Htc_Fet   <- (4.5061 * GA - 0.18487 * GA^2 + 0.0026766 * GA^3) / 100
    VFet      <- (0.0018282 * exp(15.12691 * (1 - exp(-0.077577 * GA)))) / 1000
    VPlas_Fet <- VPlasC_Fet * VFet
    VL_Fet    <- (0.0075 * exp(10.68 * (1 - exp(-0.062 * GA)))) / 1050
    VRest_Fet <- 0.93 * VFet - VPlas_Fet - VL_Fet
    QCC_Fet   <- 54
    QC_Fet    <- QCC_Fet * VPlas_Fet * (1 - Htc_Fet)
    QL_Fet    <- (6.5 / 54) * (1 - 26.5 / 75) * QC_Fet
    QRest_Fet <- QC_Fet - QL_Fet

    # ---- Maternal weight and rest-of-body volume ----
    BW_P  <- WT + (VF_P - VF_0) + (VM_P - VM_0) + VPla + VFet + VAm
    VRest <- 0.93 * BW_P - (VL + VK + VM_P + VF_P + VPlas + VPla + VFet + VAm)

    # ---- Maternal cardiac output and plasma flows ----
    QC_0 <- QCC * WT^0.75 * (1 - Htc)
    QC   <- QC_0 + 3.2512 * GA + 0.15947 * GA^2 - 0.0047059 * GA^3
    QF_P <- (0.01) * (8.5 + (-0.0175) * GA) * QC
    QF   <- QFC * QC_0
    QK_P <- (0.01) * (17 + (-0.01) * GA) * QC
    QK   <- QKC * QC_0
    QL_P <- (0.01) * (27 + (-0.175) * GA) * QC
    QL   <- QLC * QC_0
    if (GA < 3.6) {
      QPla <- 0
    } else {
      QPla <- (0.00022) * (GA - 3.6) * (0.4 + 0.29 * GA) * QC
    }
    QM   <- QMC * QC_0
    QM_P <- QM * (VM_P / VM_0)
    QRest <- QC - (QK_P + QL_P + QM_P + QF_P + QPla)

    # =====================================================================
    # 3. Scaled kinetic parameters (allometry on BW_P)
    # =====================================================================
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
    Ktrans_1 <- Ktrans1C * (VFet^0.75)
    Ktrans_2 <- Ktrans2C * (VFet^0.75)
    Ktrans_3 <- Ktrans3C * (VFet^0.75)
    Ktrans_4 <- Ktrans4C * (VFet^0.75)

    # =====================================================================
    # 4. Concentrations (mg/L). Maternal fat/mammary use the expanded
    #    gestational volumes (CM = AM/VM_P, CF = AF/VF_P), per the code.
    # =====================================================================
    CPlas_free <- plasma / VPlas
    CPlas      <- CPlas_free / fu
    CL         <- liver / VL
    CVL        <- CL / PL
    CKb        <- kidney_blood / VKb
    CVK        <- CKb
    CKid       <- CKb * PK
    CM         <- mammary / VM_P
    CVM        <- CM / PM
    CF         <- fat / VF_P
    CVF        <- CF / PF
    CAm        <- amniotic / (VAm + 1e-7)
    CRest      <- rest / VRest
    CVRest     <- CRest / PRest
    CPTC       <- ptc / VPTC
    CFil       <- filtrate / VFil
    CPla       <- placenta / (VPla + 1e-7)
    CVPla      <- CPla / PPla

    CPlas_Fet_free <- plasma_fet / (VPlas_Fet + 1e-7)
    CPlas_Fet      <- CPlas_Fet_free / fu_fet
    CRest_Fet      <- rest_fet / (VRest_Fet + 1e-7)
    CVRest_Fet     <- CRest_Fet / PRest_Fet
    CL_Fet         <- liver_fet / (VL_Fet + 1e-7)
    CVL_Fet        <- CL_Fet / PL_Fet

    # =====================================================================
    # 5. Fluxes (mg/h). Maternal Free multiplies the placenta->fetus term
    #    Rtrans_1; the fetus->mother term Rtrans_2 uses the fetal free
    #    fraction, per the authors' code.
    # =====================================================================
    RA_baso   <- Vmax_baso * (kidney_blood / VKb) / (Km_baso + kidney_blood / VKb)
    RA_apical <- Vmax_apical * (filtrate / VFil) / (Km_apical + filtrate / VFil)
    Rdif      <- Kdif * (CKb - CPTC)
    RAefflux  <- Kefflux * ptc
    RCI       <- CPlas * GFR * fu
    Rtrans_1  <- Ktrans_1 * CVPla * fu
    Rtrans_2  <- Ktrans_2 * CPlas_Fet * fu_fet
    Rtrans_3  <- Ktrans_3 * CVRest_Fet * fu_fet
    Rtrans_4  <- Ktrans_4 * CAm

    # =====================================================================
    # 6. Maternal ODEs
    # =====================================================================
    d/dt(stomach)   <- -K0 * stomach - GE * stomach
    d/dt(intestine) <- GE * stomach - Kabs * intestine - Kunabs * intestine
    d/dt(liver)     <- QL_P * (CPlas - CVL) * fu - Kbile * liver + Kabs * intestine + K0 * stomach
    d/dt(fat)       <- QF_P * (CPlas - CVF) * fu
    d/dt(mammary)   <- QM_P * (CPlas - CVM) * fu
    d/dt(rest)      <- QRest * (CPlas - CVRest) * fu
    d/dt(kidney_blood) <- QK_P * (CPlas - CVK) * fu - RCI - Rdif - RA_baso
    d/dt(ptc)          <- Rdif + RA_apical + RA_baso - RAefflux
    d/dt(filtrate)     <- RCI - RA_apical - Kurine * filtrate
    d/dt(placenta)  <- QPla * (CPlas - CVPla) * fu + Rtrans_2 - Rtrans_1
    d/dt(plasma)    <- QRest * CVRest * fu + QK_P * CVK * fu + QL_P * CVL * fu +
                       QM_P * CVM * fu + QF_P * CVF * fu + QPla * CVPla * fu -
                       QC * CPlas * fu + RAefflux
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
