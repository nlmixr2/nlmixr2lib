Chou_2021_pfos_human_prepregnant_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, non-pregnant adult woman; original implementation",
    "in R/mrgsolve). Perfluorooctane sulfonate (PFOS) in a pre-conception",
    "woman (Chou and Lin 2021, Environ Health Perspect). This is the",
    "pre-pregnancy sub-model of the gestational/lactational PFOS PBPK",
    "model: the paper runs a constant dietary intake from birth to age 30",
    "years with this model and uses the resulting steady-state body burden",
    "as the initial condition of the gestational run. Structure is a",
    "two-compartment gastrointestinal tract (stomach, small intestine)",
    "feeding the liver by the portal vein, plus plasma, liver, fat,",
    "mammary gland, rest of body and a three-subcompartment kidney (kidney",
    "blood, proximal tubule cells, filtrate) with glomerular filtration,",
    "saturable basolateral (Oat1/Oat3) and apical (Oatp1a1)",
    "transporter-mediated renal reabsorption, passive diffusion and",
    "first-order efflux back to plasma. Elimination is urinary and faecal",
    "(biliary plus unabsorbed). PFOS is not metabolised. Only the",
    "plasma-unbound fraction exchanges with tissues. Deterministic, no",
    "random effects. Model code is the authors' openly published mrgsolve",
    "source (https://github.com/KSUICCM/PFOS-Ges-Lac, ModFit/Human/HMod.R,",
    "object PreGHumanPBPK)."
  )
  reference <- paste(
    "Chou WC, Lin Z. Development of a Gestational and Lactational",
    "Physiologically Based Pharmacokinetic (PBPK) Model for",
    "Perfluorooctane Sulfonate (PFOS) in Rats and Humans and Its",
    "Implications in the Derivation of Health-Based Toxicity Values.",
    "Environ Health Perspect. 2021;129(3):037004. doi:10.1289/EHP7671.",
    "Chemical parameters from Table 2 (human, before pregnancy);",
    "physiological parameters and structure from the authors' published",
    "code, ModFit/Human/HMod.R / PreGHumanPBPK",
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
    "feces"
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
    urine = list(analyte = "PFOS", units = "mg", specimen = "urine", verified = TRUE),
    feces = list(analyte = "PFOS", units = "mg", specimen = "faeces", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Body weight of the non-pregnant adult woman. The authors' code",
        "holds BW = 60 kg (Yoon et al. 2011) for the prepregnant burn-in,",
        "which is the default carried in the vignette. Every tissue volume",
        "and blood flow is a fraction of WT, and WT drives the allometric",
        "BW^-0.25 rate-constant and BW^0.75 Vmax scaling. The paper's text",
        "describes an age-dependent BW from birth to 30 years (its Eq. 24);",
        "the published prepregnant code simplifies this to a constant BW,",
        "which is reproduced here."
      ),
      source_name = "BW"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    age_range = "birth to 30 years (simulated pre-conception exposure window)",
    weight_range = "60 kg (authors' code default)",
    sex_female_pct = 100,
    disease_state = paste(
      "Healthy non-pregnant adult woman. No clinical study underlies this",
      "sub-model; it is a forward simulation used to accumulate a",
      "pre-conception PFOS body burden from chronic dietary intake."
    ),
    dose_range = paste(
      "Chronic dietary intake. The paper simulates 0.19-4.4 ng/kg/day PFOS",
      "(Loccisano et al. 2013) from birth to age 30 years."
    ),
    regions = "worldwide human biomonitoring (exposure scenario)",
    notes = paste(
      "Chemical-specific parameters (Table 2, human, before-pregnancy",
      "column) were optimised in the authors' adult PFOS model (Chou and",
      "Lin 2019). Physiological fractions come from the authors' published",
      "mrgsolve code (HMod.R / PreGHumanPBPK)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific parameters for PFOS in the non-pregnant woman.
    # Source: Chou and Lin 2021 Table 2 (human, before-pregnancy column),
    # cross-checked against the authors' published code
    # ModFit/Human/HMod.R, object PreGHumanPBPK $PARAM block
    # (https://github.com/KSUICCM/PFOS-Ges-Lac). Every value is FIXED; the
    # model is deterministic.
    # ---------------------------------------------------------------------
    lfu <- fixed(log(0.014)); label("Free fraction of PFOS in plasma (unitless)")  # Table 2 human before-pregnancy 'Free' = 0.014

    lk0     <- fixed(log(1.000));   label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # Table 2 'K0C' = 1
    lkabs   <- fixed(log(2.120));   label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # Table 2 'KabsC' = 2.12
    lkunabs <- fixed(log(7.05e-5)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # Table 2 'KunabsC' = 7.05e-5
    lge     <- fixed(log(3.510));   label("Gastric emptying rate constant GEC (1/h per BW^0.25)")              # HMod.R PreGHumanPBPK 'GEC' = 3.510 (Yang et al. 2014)

    lkp_liver   <- fixed(log(2.03)); label("Liver-to-plasma partition coefficient PL (unitless)")           # Table 2 human 'PL' = 2.03
    lkp_kidney  <- fixed(log(1.26)); label("Kidney-to-plasma partition coefficient PK (unitless)")          # Table 2 human 'PK' = 1.26
    lkp_fat     <- fixed(log(0.13)); label("Fat-to-plasma partition coefficient PF (unitless)")             # Table 2 human 'PF' = 0.13
    lkp_mammary <- fixed(log(0.16)); label("Mammary-to-plasma partition coefficient PM (unitless)")         # Table 2 human 'PM' = 0.16
    lkp_rest    <- fixed(log(0.20)); label("Rest-of-body-to-plasma partition coefficient PRest (unitless)") # Table 2 human 'PRest' = 0.20

    lkbile  <- fixed(log(1.3e-4)); label("Biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # Table 2 human 'KbileC' = 1.3e-4
    lkurine <- fixed(log(0.096));  label("Urinary elimination rate constant KurineC (1/h per BW^-0.25)") # Table 2 human 'KurineC' = 0.096
    lgfr    <- fixed(log(27.28));  label("Glomerular filtration rate constant GFRC (L/h per kg kidney)")  # HMod.R PreGHumanPBPK 'GFRC' = 27.28 (Corley 2005)

    lvmax_baso_invitro   <- fixed(log(479));   label("In vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)")  # Table 2 human 'Vmax_baso_invitro' = 479
    lvmax_apical_invitro <- fixed(log(51803)); label("In vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")         # Table 2 human 'Vmax_apical_invitro' = 51803
    lkm_baso             <- fixed(log(20.1));  label("Michaelis constant, basolateral transporters (mg/L)")         # Table 2 human 'Km_baso' = 20.1
    lkm_apical           <- fixed(log(64.4));  label("Michaelis constant, apical transporters (mg/L)")              # Table 2 human 'Km_apical' before-pregnancy = 64.4
    lrafbaso             <- fixed(log(1));      label("Relative activity factor, basolateral transporters (unitless)") # Table 2 human 'RAFbaso' = 1
    lrafapi              <- fixed(log(0.001));  label("Relative activity factor, apical transporters (unitless)")      # Table 2 human 'RAFapi' before-pregnancy = 0.001
    lkdif                <- fixed(log(0.001));  label("Diffusion rate, kidney blood to proximal tubule cells (L/h)")   # Table 2 human 'Kdif' = 0.001
    lkefflux             <- fixed(log(0.150));  label("Efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)")  # Table 2 human 'KeffluxC' before-pregnancy = 0.150

    # Residual error placeholder: the paper reports no residual-error
    # model for this deterministic PBPK. FIXED, carries no information from
    # the paper. See the vignette Errata.
    propSd <- fixed(0.30); label("Proportional residual error placeholder, plasma (fraction)")  # not reported in Chou 2021; placeholder only
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
    # 1. Fixed physiology, non-pregnant adult woman
    #    (HMod.R / PreGHumanPBPK $PARAM). Time-invariant here.
    # =====================================================================
    Htc     <- 0.44      # HMod.R human Htc (ICRP 89, 2003)
    QCC     <- 16.4      # HMod.R human QCC (L/h/kg^0.75)
    QLC     <- 0.25      # HMod.R human QLC
    QKC     <- 0.141     # HMod.R human QKC
    QMC     <- 0.027     # HMod.R human QMC
    QFC     <- 0.052     # HMod.R human QFC
    VLC     <- 0.026     # HMod.R human VLC
    VKC     <- 0.004     # HMod.R human VKC
    VMC     <- 0.0062    # HMod.R human VMC
    VFC     <- 0.214     # HMod.R human VFC
    VPlasC  <- 0.0428    # HMod.R human VPlasC
    VFilC   <- 8.4e-4    # HMod.R human VFilC (L/kg BW)
    VPTCC   <- 1.35e-4   # HMod.R human VPTCC (L/kg kidney)
    protein <- 2.0e-6    # HMod.R human protein (mg protein/PTC)
    MW      <- 500.126   # HMod.R human MW, PFOS molecular weight (g/mol)

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
    # 4. Concentrations (mg/L)
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
    # 5. Renal fluxes (mg/h). Michaelis-Menten written with the state
    #    expression inline for numerical robustness.
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
