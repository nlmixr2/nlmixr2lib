Zhang_2024_f53b_human_prepregnancy_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, non-pregnant adult female; original",
    "implementation in R/mrgsolve). Chlorinated polyfluorinated ether",
    "sulfonic acid (F-53B, 6:2 Cl-PFESA) in a non-pregnant woman (Zhang",
    "et al. 2024, Environ Sci Technol). This is the pre-conception",
    "companion to Zhang_2024_f53b_human_pbpk: the paper simulates a",
    "constant dietary intake from birth to age 30 years with this model",
    "and uses the resulting body burden as the initial condition of the",
    "38-week gestational run. Structure is the maternal half of the",
    "gestational model with time-invariant physiology and no placental",
    "or fetal compartments: two-compartment gastrointestinal tract",
    "(stomach, small intestine) feeding the liver by the portal vein,",
    "plus plasma, liver, fat, mammary gland, rest of body and a",
    "three-subcompartment kidney (kidney blood, proximal tubule cells,",
    "filtrate) with glomerular filtration, saturable basolateral",
    "(Oat1/Oat3) and apical (Oatp1a1) reabsorption, passive diffusion",
    "and first-order efflux back to plasma. Elimination is urinary and",
    "faecal (biliary plus unabsorbed). Only the plasma-unbound fraction",
    "exchanges with tissues. Deterministic, no random effects.",
    "IMPORTANT PROVENANCE: the chemical-specific parameters of this",
    "sub-model are NOT tabulated in the paper's Supporting Information,",
    "which gives only pregnant-mouse and pregnant-human columns. They",
    "are taken from the authors' own openly published model code",
    "(https://github.com/choulab210, Code/Model_Humans.R, object",
    "PreGHumanPBPK), which the paper cites as the source of record for",
    "the model. See the vignette Errata."
  )
  reference <- paste(
    "Zhang J, Li SP, Li QQ, Zhang YT, Dong GH, Canchola A, Zeng X, Chou",
    "WC. Development of a Physiologically Based Pharmacokinetic (PBPK)",
    "Model for F-53B in Pregnant Mice and Its Extrapolation to Humans.",
    "Environ Sci Technol. 2024;58(43):18928-18939.",
    "doi:10.1021/acs.est.4c05405. Structure from Supporting Information",
    "Sections S4.1-S4.4; physiological parameters from Table S4",
    "(prepregnant human column); chemical-specific parameters from the",
    "authors' published code repository, Model_Humans.R / PreGHumanPBPK",
    "(https://github.com/choulab210)."
  )
  vignette <- "Zhang_2024_f53b_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c(
    "kidney_blood", "ptc", "filtrate", "fat", "mammary", "rest", "feces"
  )

  compartmentData <- list(
    stomach      = list(analyte = "F-53B", units = "mg", specimen = "administration site",      verified = TRUE),
    intestine    = list(analyte = "F-53B", units = "mg", specimen = "administration site",       verified = TRUE),
    plasma       = list(analyte = "F-53B", units = "mg", specimen = "plasma",                verified = TRUE),
    liver        = list(analyte = "F-53B", units = "mg", specimen = "tissue",                 verified = TRUE),
    kidney_blood = list(analyte = "F-53B", units = "mg", specimen = "tissue",          verified = TRUE),
    ptc          = list(analyte = "F-53B", units = "mg", specimen = "tissue", verified = TRUE),
    filtrate     = list(analyte = "F-53B", units = "mg", specimen = "urine",      verified = TRUE),
    fat          = list(analyte = "F-53B", units = "mg", specimen = "tissue",        verified = TRUE),
    mammary      = list(analyte = "F-53B", units = "mg", specimen = "tissue",         verified = TRUE),
    rest         = list(analyte = "F-53B", units = "mg", specimen = "tissue",          verified = TRUE),
    urine        = list(analyte = "F-53B", units = "mg", specimen = "urine",                 verified = TRUE),
    feces        = list(analyte = "F-53B", units = "mg", specimen = "faeces",                 verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Body weight of the non-pregnant adult woman. Zhang 2024 Table",
        "S4 gives 54 kg for the prepregnant human, which is the default",
        "carried in ini(). Every tissue volume and blood flow is a",
        "fraction of WT, and WT drives the allometric BW^-0.25",
        "rate-constant and BW^0.75 Vmax scaling. Because the companion",
        "gestational model uses a 60 kg pregnant reference (Table S4),",
        "a simulation that chains the two should carry the body-weight",
        "change across the handover explicitly."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = "birth to 30 years (simulated pre-conception exposure window)",
    weight_range   = "54 kg (Table S4, prepregnant human)",
    sex_female_pct = 100,
    disease_state  = paste(
      "Healthy non-pregnant adult woman. No clinical study underlies",
      "this sub-model; it is a forward simulation used to accumulate a",
      "pre-conception F-53B body burden from chronic dietary intake."
    ),
    dose_range     = paste(
      "Chronic dietary intake. Published Chinese estimated daily",
      "intakes range from 0.067 to 1.87 ng/kg/day (Table S8); the paper",
      "simulates 0.122 ng/kg/day from birth to age 30 years."
    ),
    regions        = "China (exposure scenario)",
    notes          = paste(
      "Chemical-specific parameters come from the authors' published",
      "code (Model_Humans.R, PreGHumanPBPK) rather than the Supporting",
      "Information, which tabulates only the pregnant-mouse and",
      "pregnant-human parameter sets. Where a code value can be checked",
      "against the paper's own allometric equation",
      "K_human = K_mouse * (54 / 0.025)^-0.25 = K_mouse * 0.146637,",
      "KabsC, KunabsC and KeffluxC reproduce exactly; K0C, KbileC,",
      "KurineC and Kdif do not, having been carried directly from the",
      "authors' earlier PFOS model (Chou and Lin 2019). Each ini() entry",
      "records which case applies."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific parameters for F-53B in the non-pregnant woman.
    # Source: authors' published code, Code/Model_Humans.R, object
    # PreGHumanPBPK $PARAM block (https://github.com/choulab210) --
    # these are NOT in the paper's Supporting Information. Every value is
    # FIXED; the model is deterministic.
    #
    # Allometric cross-check performed during extraction against the
    # paper's own Eq. 1 using the mouse finals of Table S7 and
    # (54 / 0.025)^-0.25 = 0.146637 is recorded per line below.
    # ---------------------------------------------------------------------

    lfu <- fixed(log(0.067)); label("Free fraction of F-53B in plasma (unitless)")  # PreGHumanPBPK Free = 0.067; equals the Table S7 human/mouse calibrated value

    lk0     <- fixed(log(0.827));    label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # PreGHumanPBPK K0C = 0.827 (allometric check gives 0.794; value carried from Chou and Lin 2019)
    lkabs   <- fixed(log(0.356));    label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # PreGHumanPBPK KabsC = 0.356 = 2.430 * 0.146637, allometric check exact
    lkunabs <- fixed(log(7.921e-5)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # PreGHumanPBPK KunabsC = 7.921e-05 = 5.400e-04 * 0.146637, allometric check exact
    lge     <- fixed(log(3.510));    label("Gastric emptying rate constant GEC (1/h per BW^-0.25)")             # PreGHumanPBPK GEC = 3.510; matches Table S4 human GEC

    lkp_liver   <- fixed(log(2.1));  label("Liver-to-plasma partition coefficient PL (unitless)")           # PreGHumanPBPK PL = 2.1; matches Table S7 human
    lkp_kidney  <- fixed(log(1.26)); label("Kidney-to-plasma partition coefficient PK (unitless)")          # PreGHumanPBPK PK = 1.26; matches Table S7 human
    lkp_fat     <- fixed(log(0.27)); label("Fat-to-plasma partition coefficient PF (unitless)")             # PreGHumanPBPK PF = 0.27 (Table S7 human prints 0.29)
    lkp_mammary <- fixed(log(0.16)); label("Mammary-to-plasma partition coefficient PM (unitless)")         # PreGHumanPBPK PM = 0.16; matches Table S7 human
    lkp_rest    <- fixed(log(0.43)); label("Rest-of-body-to-plasma partition coefficient PRest (unitless)") # PreGHumanPBPK PRest = 0.43; matches Table S7 human

    lkbile  <- fixed(log(1.35e-6)); label("Biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # PreGHumanPBPK KbileC = 1.35e-06 (allometric check gives 1.466e-06)
    lkurine <- fixed(log(0.005));   label("Urinary elimination rate constant KurineC (1/h per BW^-0.25)") # PreGHumanPBPK KurineC = 0.005 (allometric check gives 2.93e-03)
    lgfr    <- fixed(log(27.28));   label("Glomerular filtration rate constant GFRC (L/h per kg kidney)")  # PreGHumanPBPK GFRC = 27.28; matches Table S4 prepregnant human

    lvmax_baso_invitro   <- fixed(log(479));    label("In vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)")  # PreGHumanPBPK Vmax_baso_invitro = 479; matches Table S7 human
    lvmax_apical_invitro <- fixed(log(51803));  label("In vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")         # PreGHumanPBPK Vmax_apical_invitro = 51803; matches Table S7 human
    lkm_baso             <- fixed(log(64.40));  label("Michaelis constant, basolateral transporters (mg/L)")         # PreGHumanPBPK Km_baso = 64.40; matches Table S7 human
    lkm_apical           <- fixed(log(20.1));   label("Michaelis constant, apical transporters (mg/L)")              # PreGHumanPBPK Km_apical = 20.1; matches Table S7 human
    lrafbaso             <- fixed(log(1));      label("Relative activity factor, basolateral transporters (unitless)")  # PreGHumanPBPK RAFbaso = 1; matches Table S7 human
    lrafapi              <- fixed(log(0.001));  label("Relative activity factor, apical transporters (unitless)")       # PreGHumanPBPK RAFapi = 0.001; matches Table S7 human
    lkdif                <- fixed(log(1.46e-6)); label("Diffusion rate, kidney blood to proximal tubule cells (L/h)")   # PreGHumanPBPK Kdif = 1.46e-06 (allometric check gives 7.92e-06)
    lkefflux             <- fixed(log(0.821));  label("Efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)")  # PreGHumanPBPK KeffluxC = 0.821 = 5.600 * 0.146637, allometric check exact

    # Residual error placeholder: the paper reports no residual-error
    # model for this forward-simulation sub-model. FIXED, carries no
    # information from the paper. See the vignette Errata.
    propSd <- fixed(0.30); label("Proportional residual error placeholder, plasma (fraction)")  # not reported in Zhang 2024; placeholder only
  })

  model({
    # =====================================================================
    # 0. Back-transforms
    # =====================================================================
    fu      <- exp(lfu)
    K0C     <- exp(lk0)
    KabsC   <- exp(lkabs)
    KunabsC <- exp(lkunabs)
    GEC     <- exp(lge)
    PL      <- exp(lkp_liver)
    PK      <- exp(lkp_kidney)
    PF      <- exp(lkp_fat)
    PM      <- exp(lkp_mammary)
    PRest   <- exp(lkp_rest)
    KbileC  <- exp(lkbile)
    KurineC <- exp(lkurine)
    GFRC    <- exp(lgfr)
    Vmax_baso_invitro   <- exp(lvmax_baso_invitro)
    Vmax_apical_invitro <- exp(lvmax_apical_invitro)
    Km_baso   <- exp(lkm_baso)
    Km_apical <- exp(lkm_apical)
    RAFbaso   <- exp(lrafbaso)
    RAFapi    <- exp(lrafapi)
    Kdif      <- exp(lkdif)
    KeffluxC  <- exp(lkefflux)

    # =====================================================================
    # 1. Fixed physiology, non-pregnant adult female
    #    (Zhang 2024 Table S4, "Human / prepregnant"). Physiology is
    #    time-invariant here: there are no gestational growth equations.
    # =====================================================================
    Htc    <- 0.44     # Zhang 2024 Table S4, prepregnant human Htc
    QCC    <- 16.4     # Zhang 2024 Table S4, human QCC (L/h/kg^0.75)
    QLC    <- 0.25     # Zhang 2024 Table S4, human QLC
    QKC    <- 0.175    # Zhang 2024 Table S4, human QKC
    QMC    <- 0.027    # Zhang 2024 Table S4, human QMC
    QFC    <- 0.052    # Zhang 2024 Table S4, human QFC
    VLC    <- 0.026    # Zhang 2024 Table S4, human VLC
    VKC    <- 0.004    # Zhang 2024 Table S4, human VKC
    VMC    <- 0.0062   # Zhang 2024 Table S4, human VMC
    VFC    <- 0.214    # Zhang 2024 Table S4, human VFC
    VPlasC <- 0.0428   # Zhang 2024 Table S4, human VPlasC
    VFilC  <- 0.00084  # Zhang 2024 Table S4, human VFilC (L/kg BW)
    VPTCC  <- 1.35e-4  # Zhang 2024 Table S4, VPTCC (L/kg kidney)
    protein <- 2.0e-6  # Zhang 2024 Table S4, protein (mg protein/PTC)
    MW     <- 570.67   # Zhang 2024 Table S1, F-53B molecular weight (g/mol)

    # =====================================================================
    # 2. Volumes (L) and plasma flows (L/h)
    # =====================================================================
    VL    <- VLC * WT
    VK    <- VKC * WT
    VM    <- VMC * WT
    VF    <- VFC * WT
    VPlas <- VPlasC * WT
    VFil  <- VFilC * WT
    VKb   <- VK * 0.16
    VPTC  <- VK * VPTCC
    MK    <- VKC * WT * 1000   # kidney mass (g)
    VRest <- 0.93 * WT - VL - VK - VM - VF - VPlas

    QC    <- QCC * WT^0.75 * (1 - Htc)
    QK    <- QKC * QC
    QL    <- QLC * QC
    QM    <- QMC * QC
    QF    <- QFC * QC
    QRest <- QC - QK - QL - QM - QF

    # =====================================================================
    # 3. Allometrically scaled kinetic parameters
    # =====================================================================
    GFR    <- GFRC * (MK / 1000)
    GE     <- GEC * WT^(-0.25)
    K0     <- K0C * WT^(-0.25)
    Kbile  <- KbileC * WT^(-0.25)
    Kurine <- KurineC * WT^(-0.25)
    Kabs   <- KabsC * WT^(-0.25)
    Kunabs <- KunabsC * WT^(-0.25)

    PTC          <- VKC * 6e7 * 1000
    Vmax_basoC   <- Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC <- Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso    <- Vmax_basoC * WT^0.75
    Vmax_apical  <- Vmax_apicalC * WT^0.75
    Kefflux      <- KeffluxC * WT^(-0.25)

    # =====================================================================
    # 4. Concentrations
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
    CRest      <- rest / VRest
    CVRest     <- CRest / PRest
    CPTC       <- ptc / VPTC
    CFil       <- filtrate / VFil

    # =====================================================================
    # 5. Renal and biliary fluxes (mg/h). Michaelis-Menten terms written
    #    with the state expression INLINE (see the mouse model comment).
    # =====================================================================
    RA_baso   <- Vmax_baso * (kidney_blood / VKb) / (Km_baso + kidney_blood / VKb)
    RA_apical <- Vmax_apical * (filtrate / VFil) / (Km_apical + filtrate / VFil)
    Rdif      <- Kdif * (CKb - CPTC)
    RAefflux  <- Kefflux * ptc
    RCI       <- CPlas * GFR * fu

    # =====================================================================
    # 6. Mass-balance ODEs
    # =====================================================================
    d/dt(stomach)   <- -K0 * stomach - GE * stomach
    d/dt(intestine) <- GE * stomach - Kabs * intestine - Kunabs * intestine

    d/dt(liver)   <- QL * (CPlas - CVL) * fu - Kbile * liver +
                     Kabs * intestine + K0 * stomach
    d/dt(fat)     <- QF * (CPlas - CVF) * fu
    d/dt(mammary) <- QM * (CPlas - CVM) * fu
    d/dt(rest)    <- QRest * (CPlas - CVRest) * fu

    d/dt(kidney_blood) <- QK * (CPlas - CVK) * fu - RCI - Rdif - RA_baso
    d/dt(ptc)          <- Rdif + RA_apical + RA_baso - RAefflux
    d/dt(filtrate)     <- RCI - RA_apical - Kurine * filtrate

    d/dt(plasma) <- QRest * CVRest * fu + QK * CVK * fu + QL * CVL * fu +
                    QM * CVM * fu + QF * CVF * fu - QC * CPlas * fu + RAefflux

    d/dt(urine) <- Kurine * filtrate
    d/dt(feces) <- Kbile * liver + Kunabs * intestine

    # =====================================================================
    # 7. Observations
    # =====================================================================
    Cc       <- CPlas
    Cliver   <- CL
    Ckidney  <- CKid
    Cfat     <- CF
    Cmammary <- CM
    Crest    <- CRest

    Cc ~ prop(propSd)
  })
}
