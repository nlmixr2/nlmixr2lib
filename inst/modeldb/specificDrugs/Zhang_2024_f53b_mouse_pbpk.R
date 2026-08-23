Zhang_2024_f53b_mouse_pbpk <- function() {
  description <- paste(
    "Preclinical (mouse). PBPK (whole-body, gestational; original",
    "implementation in R/mrgsolve). Chlorinated polyfluorinated ether",
    "sulfonic acid (F-53B, 6:2 Cl-PFESA) in pregnant C57BL/6J mice dosed",
    "80 ug/kg orally or intravenously on gestation day 13 (Zhang et al.",
    "2024, Environ Sci Technol). Maternal model: two-compartment",
    "gastrointestinal tract (stomach, small intestine) feeding the liver",
    "by the portal vein, plus plasma, liver, fat, mammary gland,",
    "placenta, rest-of-body and a three-subcompartment kidney (kidney",
    "blood, proximal tubule cells, filtrate). Renal handling combines",
    "glomerular filtration with saturable basolateral (Oat1/Oat3) and",
    "apical (Oatp1a1) reabsorption described by Michaelis-Menten",
    "kinetics scaled from in vitro Vmax by relative activity factors,",
    "passive diffusion between kidney blood and proximal tubule cells,",
    "and first-order efflux from the tubule cells back to plasma.",
    "Elimination is urinary (from the filtrate) and faecal (biliary plus",
    "unabsorbed intestinal drug). All flow-limited tissues exchange only",
    "the plasma-unbound fraction. The fetal submodel (fetal plasma,",
    "brain, liver and rest-of-body plus amniotic fluid) is a separate",
    "circulation reached only by bidirectional placental diffusion, with",
    "a second bidirectional exchange between the fetal rest-of-body and",
    "the amniotic fluid. Maternal and fetal physiology are time-varying",
    "functions of gestational day (O'Flaherty growth equations), so",
    "model time is hours since GD0 and the GD13 dose is given at t =",
    "312 h. Intended for typical-value simulation; the paper reports no",
    "between-subject random effects."
  )
  reference <- paste(
    "Zhang J, Li SP, Li QQ, Zhang YT, Dong GH, Canchola A, Zeng X, Chou",
    "WC. Development of a Physiologically Based Pharmacokinetic (PBPK)",
    "Model for F-53B in Pregnant Mice and Its Extrapolation to Humans.",
    "Environ Sci Technol. 2024;58(43):18928-18939.",
    "doi:10.1021/acs.est.4c05405. Model equations and parameters from",
    "Supporting Information Sections S4-S5 and Tables S3-S5, S7."
  )
  vignette <- "Zhang_2024_f53b_pbpk"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # States not covered by the canonical compartment register. These are the
  # gestational-PBPK states of Zhang 2024 Figure 1 (maternal kidney
  # subcompartments, gestational tissues, and the fetal submodel). The
  # `_fet` suffix marks a fetal-circulation state. `trans_mf` / `trans_fm`
  # are cumulative placental-transfer accumulators carried only so the
  # vignette can run an exact maternal/fetal mass-balance gate; they do not
  # feed back into the dynamics.
  paper_specific_compartments <- c(
    "kidney_blood", "ptc", "filtrate", "fat", "mammary", "rest",
    "placenta", "feces", "plasma_fet", "liver_fet", "brain_fet",
    "rest_fet", "amniotic", "trans_mf", "trans_fm"
  )

  compartmentData <- list(
    stomach      = list(analyte = "F-53B", units = "mg", specimen = "administration site",     verified = TRUE),
    intestine    = list(analyte = "F-53B", units = "mg", specimen = "administration site",      verified = TRUE),
    plasma       = list(analyte = "F-53B", units = "mg", specimen = "plasma",      verified = TRUE),
    liver        = list(analyte = "F-53B", units = "mg", specimen = "tissue",       verified = TRUE),
    kidney_blood = list(analyte = "F-53B", units = "mg", specimen = "tissue",         verified = TRUE),
    ptc          = list(analyte = "F-53B", units = "mg", specimen = "tissue", verified = TRUE),
    filtrate     = list(analyte = "F-53B", units = "mg", specimen = "urine",     verified = TRUE),
    fat          = list(analyte = "F-53B", units = "mg", specimen = "tissue",         verified = TRUE),
    mammary      = list(analyte = "F-53B", units = "mg", specimen = "tissue",        verified = TRUE),
    rest         = list(analyte = "F-53B", units = "mg", specimen = "tissue", verified = TRUE),
    placenta     = list(analyte = "F-53B", units = "mg", specimen = "tissue",             verified = TRUE),
    urine        = list(analyte = "F-53B", units = "mg", specimen = "urine",                verified = TRUE),
    feces        = list(analyte = "F-53B", units = "mg", specimen = "faeces",                verified = TRUE),
    plasma_fet   = list(analyte = "F-53B", units = "mg", specimen = "plasma",         verified = TRUE),
    liver_fet    = list(analyte = "F-53B", units = "mg", specimen = "tissue",          verified = TRUE),
    brain_fet    = list(analyte = "F-53B", units = "mg", specimen = "tissue",          verified = TRUE),
    rest_fet     = list(analyte = "F-53B", units = "mg", specimen = "tissue",   verified = TRUE),
    amniotic     = list(analyte = "F-53B", units = "mg", specimen = "tissue",       verified = TRUE),
    trans_mf     = list(analyte = "F-53B", units = "mg", specimen = "not applicable", verified = TRUE),
    trans_fm     = list(analyte = "F-53B", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Maternal pre-gestational body weight. Zhang 2024 Table S4 gives",
        "0.025 kg for the pregnant C57BL/6J dams (in-house TK study);",
        "that is the default carried in ini(). Every maternal tissue",
        "volume and blood flow is a fraction of WT (Table S4), and the",
        "gestational body weight BW_P = WT + the fat, mammary, placental,",
        "fetal and amniotic increments (Table S5) drives all allometric",
        "BW^-0.25 rate-constant and BW^0.75 Vmax scaling. WT is the most",
        "influential parameter in the paper's own local sensitivity",
        "analysis (Table S13: normalized sensitivity coefficient -0.89 to",
        "-0.90 across every maternal and fetal AUC) and is the one",
        "physiological parameter varied in the Monte Carlo analysis",
        "(Table S9, normal, CV 30%)."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "mouse (C57BL/6J, pregnant)",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "adult, gestation day 13 at dosing",
    weight_range   = "0.025 kg (Table S4)",
    sex_female_pct = 100,
    disease_state  = paste(
      "Healthy pregnant C57BL/6J dams. Eighty females and males were",
      "mated 1:1; pregnant dams were split into an oral group (n = 20)",
      "and an intravenous group (n = 20). Maternal plasma, liver,",
      "placenta, urine and faeces and fetal brain and liver were",
      "sampled from GD13 to GD17, with amniotic fluid, fat, brain,",
      "heart, spleen, kidney, small intestine and stomach added at GD17.",
      "Mean litter size was 8 fetuses. Approved by the ethics committee",
      "of Sun Yat-sen University (SYSUIACUC-2022-001602)."
    ),
    dose_range     = paste(
      "Single 80 ug/kg body weight dose of F-53B on gestation day 13,",
      "given either orally or by lateral tail vein injection. At the",
      "0.025 kg reference weight that is 0.002 mg."
    ),
    regions        = "Sun Yat-sen University, Guangzhou, China",
    notes          = paste(
      "F-53B was measured by UPLC-MS/MS. Plasma protein binding was",
      "98.70% (oral) and 99.54% (IV) by ultrafiltration (Table S11);",
      "the model's calibrated free fraction is 0.067. The PBPK model was",
      "calibrated against maternal plasma, maternal liver, placenta,",
      "fetal brain and fetal liver and evaluated against maternal fat,",
      "maternal rest-of-body, amniotic fluid and fetal-liver data held",
      "out of calibration (Table 1). Adjusted R-squared was 0.68 with",
      "more than 72% of predictions within 2-fold of observation. The",
      "structural model was inherited from the authors' PFOS PBPK models",
      "(Chou and Lin 2019, 2021). No between-subject random effects were",
      "estimated; population variability was addressed only by a Monte",
      "Carlo analysis of the human model (Table S9)."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Chemical-specific parameters for F-53B in the pregnant mouse.
    # Zhang 2024 Supporting Information Table S7, "Mice / Pregnant" column.
    # Values footnoted "a" in Table S7 were calibrated in this study; the
    # remainder are literature values carried in from the cited PFOS/PFOA
    # models. Every value is FIXED here: the paper reports a deterministic
    # calibrated PBPK model with no estimated random effects, so nothing in
    # this file is an nlmixr2 estimation target.
    #
    # Cross-check performed during extraction: every Table S7 row footnoted
    # "b" in the Human column reproduces exactly as
    # K_human = K_mouse * (60 / 0.025)^-0.25 = K_mouse * 0.14288, which
    # confirms the Mice column holds the FINAL calibrated values (the
    # authors' GitHub $PARAM block holds an earlier, pre-calibration
    # iteration -- see the vignette Errata).
    # ---------------------------------------------------------------------

    # Plasma protein binding
    lfu     <- fixed(log(0.067)); label("Free fraction of F-53B in maternal plasma (unitless)")  # Zhang 2024 Table S7, Free = 0.067 (calibrated)
    lfu_fet <- fixed(log(0.076)); label("Free fraction of F-53B in fetal plasma (unitless)")     # Zhang 2024 Table S7, Free_Fet = 0.076 (calibrated)

    # Absorption / gastric emptying
    lk0     <- fixed(log(5.42));     label("Stomach absorption rate constant K0C (1/h per BW^-0.25)")           # Zhang 2024 Table S7, K0C = 5.42 (calibrated)
    lkabs   <- fixed(log(2.430));    label("Small-intestine absorption rate constant KabsC (1/h per BW^-0.25)") # Zhang 2024 Table S7, KabsC = 2.430
    lkunabs <- fixed(log(5.400e-4)); label("Unabsorbed-to-faeces rate constant KunabsC (1/h per BW^-0.25)")     # Zhang 2024 Table S7, KunabsC = 5.400E-04
    lge     <- fixed(log(0.54));     label("Gastric emptying rate constant GEC (1/h per BW^-0.25)")             # Zhang 2024 Table S4, mouse GEC = 0.54

    # Tissue-to-plasma partition coefficients
    lkp_liver     <- fixed(log(2.10)); label("Liver-to-plasma partition coefficient PL (unitless)")            # Zhang 2024 Table S7, PL = 2.10
    lkp_kidney    <- fixed(log(0.80)); label("Kidney-to-plasma partition coefficient PK (unitless)")           # Zhang 2024 Table S7, PK = 0.80
    lkp_fat       <- fixed(log(0.29)); label("Fat-to-plasma partition coefficient PF (unitless)")              # Zhang 2024 Table S7, PF = 0.29 (calibrated)
    lkp_mammary   <- fixed(log(0.16)); label("Mammary-to-plasma partition coefficient PM (unitless)")          # Zhang 2024 Table S7, PM = 0.16
    lkp_placenta  <- fixed(log(0.08)); label("Placenta-to-plasma partition coefficient PPla (unitless)")       # Zhang 2024 Table S7, PPla = 0.08 (calibrated)
    lkp_rest      <- fixed(log(0.43)); label("Rest-of-body-to-plasma partition coefficient PRest (unitless)")  # Zhang 2024 Table S7, PRest = 0.43 (calibrated)
    lkp_liver_fet <- fixed(log(2.10)); label("Fetal liver-to-plasma partition coefficient PL_Fet (unitless)")  # Zhang 2024 Table S7, PL_Fet = 2.10 (assumed same as dam)
    lkp_brain_fet <- fixed(log(1.55)); label("Fetal brain-to-plasma partition coefficient PB_Fet (unitless)")  # Zhang 2024 Table S7, PB_Fet = 1.55
    lkp_rest_fet  <- fixed(log(0.22)); label("Fetal rest-of-body-to-plasma partition coefficient PRest_Fet (unitless)") # Zhang 2024 Table S7, PRest_Fet = 0.22

    # Placental and amniotic-fluid bidirectional transfer
    lktrans1 <- fixed(log(1.39));  label("Dam-to-fetus placental transfer constant Ktrans1C (L/h per kg^0.75)")      # Zhang 2024 Table S7, Ktrans1C = 1.39 (calibrated)
    lktrans2 <- fixed(log(0.60));  label("Fetus-to-dam placental transfer constant Ktrans2C (L/h per kg^0.75)")      # Zhang 2024 Table S7, Ktrans2C = 0.60 (calibrated)
    lktrans3 <- fixed(log(0.230)); label("Fetus-to-amniotic-fluid transfer constant Ktrans3C (L/h per kg^0.75)")     # Zhang 2024 Table S7, Ktrans3C = 0.230
    lktrans4 <- fixed(log(0.001)); label("Amniotic-fluid-to-fetus transfer constant Ktrans4C (L/h per kg^0.75)")     # Zhang 2024 Table S7, Ktrans4C = 0.001

    # Elimination
    lkbile  <- fixed(log(0.00001)); label("Biliary elimination rate constant KbileC (1/h per BW^-0.25)")  # Zhang 2024 Table S7, KbileC = 0.00001 (calibrated)
    lkurine <- fixed(log(0.02));    label("Urinary elimination rate constant KurineC (1/h per BW^-0.25)") # Zhang 2024 Table S7, KurineC = 0.02 (calibrated)
    lgfr    <- fixed(log(59));      label("Glomerular filtration rate constant GFRC (L/h per kg kidney)") # Zhang 2024 Table S4, mouse GFRC = 59

    # Renal reabsorption (in vitro Vmax scaled by relative activity factors)
    lvmax_baso_invitro   <- fixed(log(393.450));  label("In vitro Vmax, basolateral Oat1/Oat3 (pmol/mg protein/min)") # Zhang 2024 Table S7, Vmax_baso_invitro = 393.450
    lvmax_apical_invitro <- fixed(log(1632.576)); label("In vitro Vmax, apical Oatp1a1 (pmol/mg protein/min)")        # Zhang 2024 Table S7, Vmax_apical_invitro = 1632.576 (calibrated)
    lkm_baso             <- fixed(log(50.157));   label("Michaelis constant, basolateral transporters (mg/L)")        # Zhang 2024 Table S7, Km_baso = 50.157 (calibrated)
    lkm_apical           <- fixed(log(140.987));  label("Michaelis constant, apical transporters (mg/L)")             # Zhang 2024 Table S7, Km_apical = 140.987 (calibrated)
    lrafbaso             <- fixed(log(3.990));    label("Relative activity factor, basolateral transporters (unitless)") # Zhang 2024 Table S7, RAFbaso = 3.990
    lrafapi              <- fixed(log(2.810));    label("Relative activity factor, apical transporters (unitless)")      # Zhang 2024 Table S7, RAFapi = 2.810
    lkdif                <- fixed(log(5.4e-5));   label("Diffusion rate, kidney blood to proximal tubule cells (L/h)")   # Zhang 2024 Table S7, Kdif = 5.4E-05 (calibrated)
    lkefflux             <- fixed(log(5.600));    label("Efflux rate constant, tubule cells to plasma KeffluxC (1/h per BW^-0.25)") # Zhang 2024 Table S7, KeffluxC = 5.600

    # Litter size. Not a subject-level covariate: it is a study constant of
    # the in-house TK experiment and enters the placental, fetal and
    # amniotic growth equations of Table S5 as the whole-litter multiplier.
    n_fetus <- fixed(8); label("Number of fetuses per litter (count)")  # Zhang 2024 Table S4, mouse N = 8 (in-house experiment)

    # Residual error. Zhang 2024 fit the PBPK model deterministically by
    # least squares against tissue time courses and reports no residual-
    # error model (no sigma, no CV%, no weighting scheme is tabulated).
    # These are FIXED placeholders present only so the model is a complete
    # nlmixr2 object; they carry no information from the paper and must not
    # be interpreted as published estimates. See the vignette Errata.
    propSd            <- fixed(0.30); label("Proportional residual error placeholder, maternal plasma (fraction)")  # not reported in Zhang 2024; placeholder only
    propSd_Cliver     <- fixed(0.30); label("Proportional residual error placeholder, maternal liver (fraction)")   # not reported in Zhang 2024; placeholder only
    propSd_Cplacenta  <- fixed(0.30); label("Proportional residual error placeholder, placenta (fraction)")         # not reported in Zhang 2024; placeholder only
    propSd_Cbrain_fet <- fixed(0.30); label("Proportional residual error placeholder, fetal brain (fraction)")      # not reported in Zhang 2024; placeholder only
    propSd_Cliver_fet <- fixed(0.30); label("Proportional residual error placeholder, fetal liver (fraction)")      # not reported in Zhang 2024; placeholder only
  })

  model({
    # =====================================================================
    # 0. Back-transforms
    # =====================================================================
    fu     <- exp(lfu)
    fu_fet <- exp(lfu_fet)
    K0C    <- exp(lk0)
    KabsC  <- exp(lkabs)
    KunabsC <- exp(lkunabs)
    GEC    <- exp(lge)
    PL     <- exp(lkp_liver)
    PK     <- exp(lkp_kidney)
    PF     <- exp(lkp_fat)
    PM     <- exp(lkp_mammary)
    PPla   <- exp(lkp_placenta)
    PRest  <- exp(lkp_rest)
    PL_Fet <- exp(lkp_liver_fet)
    PB_Fet <- exp(lkp_brain_fet)
    PRest_Fet <- exp(lkp_rest_fet)
    Ktrans1C <- exp(lktrans1)
    Ktrans2C <- exp(lktrans2)
    Ktrans3C <- exp(lktrans3)
    Ktrans4C <- exp(lktrans4)
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
    N         <- n_fetus

    # =====================================================================
    # 1. Time base
    #
    # Model time is HOURS SINCE GD0, because every maternal and fetal
    # volume and blood flow below is an explicit function of gestational
    # day. The GD13 dose is therefore administered at t = 312 h and the
    # 96 h observation window runs to t = 408 h (GD17). This matches the
    # authors' own data files, whose Time column is hours since GD0.
    # =====================================================================
    GD <- time / 24  # Zhang 2024 SI code: GD = TIME/24

    # =====================================================================
    # 2. Fixed maternal physiology (Zhang 2024 Table S4, "Mice / Pregnant")
    #    Volumes are fractions of body weight; flows are fractions of
    #    cardiac output. QCC is a plasma-flow scalar (the (1 - Htc) factor
    #    below converts whole-blood cardiac output to plasma flow).
    # =====================================================================
    Htc    <- 0.48    # Zhang 2024 Table S4, mouse Htc
    QCC    <- 16.5    # Zhang 2024 Table S4, mouse QCC (L/h/kg^0.75)
    QLC    <- 0.161   # Zhang 2024 Table S4, mouse QLC
    QKC    <- 0.091   # Zhang 2024 Table S4, mouse QKC
    QMC    <- 0.002   # Zhang 2024 Table S4, mouse QMC
    QFC    <- 0.07    # Zhang 2024 Table S4, mouse QFC
    VLC    <- 0.055   # Zhang 2024 Table S4, mouse VLC
    VKC    <- 0.017   # Zhang 2024 Table S4, mouse VKC
    VMC    <- 0.01    # Zhang 2024 Table S4, mouse VMC (rat value used)
    VFC    <- 0.068   # Zhang 2024 Table S4, mouse VFC
    VPlasC <- 0.049   # Zhang 2024 Table S4, mouse VPlasC
    VFilC  <- 0.0017  # Zhang 2024 Table S4, mouse VFilC (L/kg BW)
    VPTCC  <- 1.35e-4 # Zhang 2024 Table S4, VPTCC (L/g kidney)
    protein <- 2.0e-6 # Zhang 2024 Table S4, protein (mg protein/PTC)
    MW     <- 570.67  # Zhang 2024 Table S1, F-53B molecular weight (g/mol)

    # Fetal physiological constants (Zhang 2024 Table S4, fetal block)
    VPlasC_Fet <- 0.049  # Zhang 2024 Table S4, VPlasC_Fet (assumed same as dam)
    QLC_Fet    <- 0.161  # Zhang 2024 Table S4, QLC_Fet (assumed same as dam)
    QBC_Fet    <- 0.1055 # Zhang 2024 Table S4, QBC_Fet

    # Rat-to-mouse body-weight correction factor applied to the mammary,
    # fat, amniotic-fluid and fetal-brain growth equations, which were
    # published for the rat (Zhang 2024 Table S5 footnote d: "Correction of
    # rat parameters by weight"; SI code uses 30/300).
    RATMOUSE <- 30.0 / 300.0

    # =====================================================================
    # 3. Maternal gestational growth equations (Zhang 2024 Table S5)
    # =====================================================================
    # Placenta volume for the whole litter (L). Segments are exclusive:
    # nothing before GD6, decidual growth to GD9.25, then decidual
    # regression plus chorioallantoic capillary growth. The 1e6 divisor
    # converts the published mm^3 to L. A 1e-10 floor keeps CPla finite
    # before implantation (the authors' code uses the same guard).
    VPla <- 1e-10
    if (GD > 6 && GD <= 9.25) {
      VPla <- N * 8 * (GD - 6) / 1e6
    }
    if (GD > 9.25) {
      VPla <- (N * 32 * exp(-0.23 * (GD - 9.25)) +
               N * 40 * (exp(0.28 * (GD - 9.25)) - 1)) / 1e6
    }

    # Plasma flow to the placenta (L/h). Qdec1 (yolk-sac placenta), Qdec2
    # (yolk-sac regression) and Qcap (chorioallantoic placenta) are
    # mutually exclusive gestational segments.
    Qdec1 <- 0
    Qdec2 <- 0
    Qcap  <- 0
    if (GD > 5.5 && GD <= 9.25) {
      Qdec1 <- 0.55 * (GD - 5.5)
    }
    if (GD > 9.25 && GD <= 11) {
      Qdec2 <- 2.2 * exp(-0.23 * (GD - 9.25))
    }
    if (GD > 11) {
      Qcap <- (0.1207 * (GD - 11))^4.36
    }
    Qdec <- Qdec1 + Qdec2
    QPla <- (N * (0.02 * Qdec + Qcap) / 24) * (1 - Htc)

    # Maternal tissue volumes (L). Mammary and fat expand with gestation.
    VM0   <- VMC * WT
    VM    <- VM0 * (1 + 0.27 * GD * RATMOUSE)
    VF0   <- VFC * WT
    VF    <- VF0 * (1 + 0.0165 * GD * RATMOUSE)
    VL    <- VLC * WT
    VK    <- VKC * WT
    VPlas <- VPlasC * WT
    VFil  <- VFilC * WT
    MK    <- VKC * WT * 1000          # kidney mass (g)
    VPTC  <- MK * VPTCC               # proximal tubule cell volume (L)
    VKb   <- VK * 0.16                # kidney blood volume (L)

    # =====================================================================
    # 4. Fetal gestational growth equations (Zhang 2024 Table S5, fetus)
    # =====================================================================
    # Volume of one fetus (kg == L). The three O'Flaherty segments are
    # cumulative: the GD8.6-15.8 segment adds to the frozen GD8.6 value of
    # the first, and the final segment interpolates linearly from the
    # GD15.8 value (273.38 mm^3) to the 1250 mm^3 term size at GD19.
    VFet_1 <- (0.12 * GD)^4.53 / 1e6
    if (GD > 8.6 && GD <= 15.8) {
      VFet_1 <- ((0.12 * 8.6)^4.53 + (1.2 * (GD - 8.6))^2.6) / 1e6
    }
    if (GD > 15.8) {
      VFet_1 <- (273.38 + ((GD - 15.8) / (19 - 15.8)) * (1250 - 273.38)) / 1e6
    }
    VFet <- VFet_1 * N

    # Amniotic fluid (L), whole litter; zero before GD8.
    VAmX <- 1e-10
    if (GD >= 8) {
      VAmX <- ((-4e-6) * GD^3 + 0.0002 * GD^2 - 0.0023 * GD + 0.0099) * RATMOUSE
    }
    VAm <- VAmX * N

    VPlas_Fet <- VPlasC_Fet * VFet
    VL_Fet    <- (0.406 / (1 + exp((14.716 - GD) / 0.907))) / 1000
    VB_Fet    <- (4.191 * exp(-exp(2.554 - 0.06726 * GD)) / 1000) * N * RATMOUSE
    VRest_Fet <- 0.93 * VFet - VPlas_Fet - VL_Fet - VB_Fet
    # As published, the Table S5 growth equations make the fetal
    # rest-of-body volume NEGATIVE before about GD3.5, because the fetal
    # liver and brain equations are non-zero while whole-fetus volume is
    # still negligible. A negative volume makes the fetal concentrations
    # diverge. Floor it at 1e-9 L. The window is physiologically inert:
    # placental plasma flow is zero until GD5.5 and the placental
    # transfer clearances scale as VFet_1^0.75, so no F-53B can reach the
    # fetus while the floor is binding (and this model is dosed at
    # GD13). See the vignette Errata.
    if (VRest_Fet < 1e-9) {
      VRest_Fet <- 1e-9
    }

    # Fetal plasma flows (L/h)
    QC_Fet    <- QPla / (1 + 20000 * exp(-0.55 * GD))
    QL_Fet    <- QC_Fet * QLC_Fet
    QB_Fet    <- QC_Fet * QBC_Fet
    QRest_Fet <- QC_Fet - QL_Fet - QB_Fet

    # =====================================================================
    # 5. Gestational maternal body weight, remaining volume and flows
    # =====================================================================
    BW_P  <- WT + (VF - VF0) + (VM - VM0) + VPla + VFet + VAm
    VRest <- 0.93 * BW_P - (VL + VK + VM + VF + VPlas + VPla + VFet + VAm)

    QC    <- QCC * WT^0.75 * (1 - Htc)
    QC_P  <- QC + (N * (Qdec + Qcap) / 24) * (1 - Htc)
    QF    <- QFC * QC
    QF_P  <- QF * (VF / VF0)
    QM    <- QMC * QC
    QM_P  <- QM * (VM / VM0)
    QL    <- QLC * QC
    QK    <- QKC * QC
    QRest <- QC_P - (QK + QL + QM_P + QF_P + QPla)

    # =====================================================================
    # 6. Allometrically scaled kinetic parameters
    #    Rate constants scale as BW_P^-0.25 and Vmax as BW_P^0.75, both on
    #    the GESTATIONAL body weight (Zhang 2024 SI code).
    # =====================================================================
    GFR    <- GFRC * (MK / 1000)   # L/h; GFRC is per kg kidney
    GE     <- GEC * BW_P^(-0.25)
    K0     <- K0C * BW_P^(-0.25)
    Kbile  <- KbileC * BW_P^(-0.25)
    Kurine <- KurineC * BW_P^(-0.25)
    Kabs   <- KabsC * BW_P^(-0.25)
    Kunabs <- KunabsC * BW_P^(-0.25)

    # In vitro to in vivo extrapolation of the renal transporter capacity.
    # PTC is the proximal tubule cell count per kg body weight (60 million
    # cells per gram kidney); the 60 converts per-minute to per-hour and
    # MW/1e12 * 1000 converts pmol to mg.
    PTC          <- VKC * 6e7 * 1000
    Vmax_basoC   <- Vmax_baso_invitro * RAFbaso * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_apicalC <- Vmax_apical_invitro * RAFapi * PTC * protein * 60 * (MW / 1e12) * 1000
    Vmax_baso    <- Vmax_basoC * BW_P^0.75
    Vmax_apical  <- Vmax_apicalC * BW_P^0.75
    Kefflux      <- KeffluxC * BW_P^(-0.25)

    # Placental / amniotic transfer clearances (L/h), scaled by the
    # 0.75 power of single-fetus volume times litter size.
    Ktrans_1 <- Ktrans1C * VFet_1^0.75 * N
    Ktrans_2 <- Ktrans2C * VFet_1^0.75 * N
    Ktrans_3 <- Ktrans3C * VFet_1^0.75 * N
    Ktrans_4 <- Ktrans4C * VFet_1^0.75 * N

    # =====================================================================
    # 7. Concentrations
    #    The plasma state holds the amount associated with the unbound
    #    pool; the measurable TOTAL plasma concentration is recovered by
    #    dividing by the free fraction (Zhang 2024 SI code:
    #    CPlas = CPlas_free / Free).
    # =====================================================================
    CPlas_free <- plasma / VPlas
    CPlas      <- CPlas_free / fu
    CL         <- liver / VL
    CVL        <- CL / PL
    CKb        <- kidney_blood / VKb
    CVK        <- CKb
    CKid       <- CVK * PK
    CM         <- mammary / VM
    CVM        <- CM / PM
    CF         <- fat / VF
    CVF        <- CF / PF
    CRest      <- rest / VRest
    CVRest     <- CRest / PRest
    CPTC       <- ptc / VPTC
    CFil       <- filtrate / VFil
    CPla       <- placenta / (VPla + 1e-10)
    CVPla      <- CPla / PPla
    CAm        <- amniotic / (VAm + 1e-10)

    CPlas_Fet_free <- plasma_fet / (VPlas_Fet + 1e-10)
    CPlas_Fet      <- CPlas_Fet_free / fu_fet
    CRest_Fet      <- rest_fet / (VRest_Fet + 1e-10)
    CVRest_Fet     <- CRest_Fet / PRest_Fet
    CL_Fet         <- liver_fet / (VL_Fet + 1e-10)
    CVL_Fet        <- CL_Fet / PL_Fet
    CB_Fet         <- brain_fet / (VB_Fet + 1e-10)
    CVB_Fet        <- CB_Fet / PB_Fet

    # =====================================================================
    # 8. Renal, biliary and placental fluxes (mg/h)
    #    The two Michaelis-Menten terms are written with the state
    #    expression INLINE rather than through the named CKb / CFil
    #    intermediates: a named intermediate that reads an ODE state can
    #    silently evaluate to zero inside a d/dt() in the nlmixr2 model
    #    function form, which would delete the saturable reabsorption
    #    without any error.
    # =====================================================================
    RA_baso   <- Vmax_baso * (kidney_blood / VKb) / (Km_baso + kidney_blood / VKb)
    RA_apical <- Vmax_apical * (filtrate / VFil) / (Km_apical + filtrate / VFil)
    Rdif      <- Kdif * (CKb - CPTC)
    RAefflux  <- Kefflux * ptc
    RCI       <- CPlas * GFR * fu

    # Bidirectional placental and amniotic transfer. Rtrans_2 uses the
    # MATERNAL free fraction as published (Zhang 2024 eq S17 and SI code).
    Rtrans_1 <- Ktrans_1 * CVPla * fu
    Rtrans_2 <- Ktrans_2 * CPlas_Fet * fu
    Rtrans_3 <- Ktrans_3 * CVRest_Fet * fu_fet
    Rtrans_4 <- Ktrans_4 * CAm

    # =====================================================================
    # 9. Mass-balance ODEs
    # =====================================================================
    # Gastrointestinal tract (Zhang 2024 eqs S3-S5)
    d/dt(stomach)   <- -K0 * stomach - GE * stomach
    d/dt(intestine) <- GE * stomach - Kabs * intestine - Kunabs * intestine

    # Maternal flow-limited tissues (Zhang 2024 eqs S14-S15). The liver
    # additionally receives portal absorption from stomach and intestine
    # and loses drug to bile.
    d/dt(liver)   <- QL * (CPlas - CVL) * fu - Kbile * liver +
                     Kabs * intestine + K0 * stomach
    d/dt(fat)     <- QF_P * (CPlas - CVF) * fu
    d/dt(mammary) <- QM_P * (CPlas - CVM) * fu
    d/dt(rest)    <- QRest * (CPlas - CVRest) * fu

    # Kidney subcompartments (Zhang 2024 eqs S6-S11)
    d/dt(kidney_blood) <- QK * (CPlas - CVK) * fu - RCI - Rdif - RA_baso
    d/dt(ptc)          <- Rdif + RA_apical + RA_baso - RAefflux
    d/dt(filtrate)     <- RCI - RA_apical - Kurine * filtrate

    # Placenta: perfused from maternal plasma and exchanging with the
    # fetal circulation (Zhang 2024 eqs S16-S17)
    d/dt(placenta) <- QPla * (CPlas - CVPla) * fu + Rtrans_2 - Rtrans_1

    # Maternal plasma. Not printed explicitly in the SI equation list; it
    # is the venous-return closure of the flow terms above plus the
    # proximal-tubule efflux, exactly as in the authors' SI code.
    d/dt(plasma) <- QRest * CVRest * fu + QK * CVK * fu + QL * CVL * fu +
                    QM_P * CVM * fu + QF_P * CVF * fu + QPla * CVPla * fu -
                    QC_P * CPlas * fu + RAefflux

    # Excreta (Zhang 2024 eqs S12-S13)
    d/dt(urine) <- Kurine * filtrate
    d/dt(feces) <- Kbile * liver + Kunabs * intestine

    # Fetal circulation (Zhang 2024 eqs S16-S19)
    d/dt(liver_fet) <- QL_Fet * (CPlas_Fet - CVL_Fet) * fu_fet
    d/dt(brain_fet) <- QB_Fet * (CPlas_Fet - CVB_Fet) * fu_fet
    d/dt(rest_fet)  <- QRest_Fet * (CPlas_Fet - CVRest_Fet) * fu_fet -
                       Rtrans_3 + Rtrans_4
    d/dt(plasma_fet) <- QRest_Fet * CVRest_Fet * fu_fet +
                        QL_Fet * CVL_Fet * fu_fet +
                        QB_Fet * CVB_Fet * fu_fet -
                        QC_Fet * CPlas_Fet * fu_fet + Rtrans_1 - Rtrans_2
    d/dt(amniotic)  <- Rtrans_3 - Rtrans_4

    # Cumulative placental transfer, carried for the mass-balance gate only
    d/dt(trans_mf) <- Rtrans_1
    d/dt(trans_fm) <- Rtrans_2

    # =====================================================================
    # 10. Observations
    # =====================================================================
    Cc          <- CPlas      # maternal plasma, total concentration (mg/L)
    Cliver      <- CL
    Cplacenta   <- CPla
    Cbrain_fet  <- CB_Fet
    Cliver_fet  <- CL_Fet
    Cplasma_fet <- CPlas_Fet
    Camniotic   <- CAm
    Ckidney     <- CKid
    Cfat        <- CF
    Cmammary    <- CM
    Crest       <- CRest
    Crest_fet   <- CRest_Fet

    Cc         ~ prop(propSd)
    Cliver     ~ prop(propSd_Cliver)
    Cplacenta  ~ prop(propSd_Cplacenta)
    Cbrain_fet ~ prop(propSd_Cbrain_fet)
    Cliver_fet ~ prop(propSd_Cliver_fet)
  })
}
