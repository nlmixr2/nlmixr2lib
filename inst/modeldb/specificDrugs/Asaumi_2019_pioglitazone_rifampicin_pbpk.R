Asaumi_2019_pioglitazone_rifampicin_pbpk <- function() {
  description <- paste(
    "PBPK (five-unit tandem liver, segregated-flow intestine for rifampicin;",
    "Napp). Rifampicin-pioglitazone drug-drug interaction in healthy adults,",
    "used to estimate the rifampicin CYP2C8 induction parameter. The rifampicin",
    "perpetrator model (saturable OATP1B uptake, UGT autoinduction) is coupled",
    "to a perfusion-limited pioglitazone victim model whose hepatic metabolism",
    "is 83.6% CYP2C8 and 16.4% CYP3A; rifampicin induces hepatic CYP2C8 (Emax",
    "2.55, fitted here) and CYP3A (Emax 4.57) by turnover models (EC50,u",
    "0.0639 uM). Deterministic: no IIV or residual error published."
  )
  reference <- paste(
    "Asaumi R, Menzel K, Lee W, Nunoya K, Imawaka H, Kusuhara H, Sugiyama Y.",
    "Expanded Physiologically-Based Pharmacokinetic Model of Rifampicin for",
    "Predicting Interactions With Drugs and an Endogenous Biomarker via Complex",
    "Mechanisms Including Organic Anion Transporting Polypeptide 1B Induction.",
    "CPT Pharmacometrics Syst Pharmacol. 2019;8:845.",
    "doi:10.1002/psp4.12457 (PMC6875706).",
    "Model code, ODEs and Tables S1-S4 are the paper's Supporting Information",
    "(PSP4-8-845-s003..s008).",
    ""
  )
  vignette <- "Asaumi_2019_rifampicin_ddi_pbpk"

  # State names follow the sibling Sugiyama-laboratory PBPK models
  # (Toshimoto_2017_irinotecan_pbpk.R, Tsuchitani_2024_telmisartan_pbpk.R):
  # `is_liver<n>` = hepatic extracellular (sinusoidal) unit n,
  # `int_liver<n>` = hepatocyte unit n of the five-unit tandem liver,
  # `serosa` = gut serosa, `intestine_ent` / `intestine_muc` = enterocyte /
  # mucosal blood of the segregated-flow intestine, `ehc<n>` = biliary
  # transit chain. The victim drug carries the bare names and rifampicin
  # the registered `_rif` suffix. The induction states hold relative
  # enzyme / transporter amounts per hepatocyte unit (`<enzyme>_liver<n>`)
  # or in the enterocyte (`<enzyme>_gut`); none is registered, so all are
  # declared paper-specific.
  paper_specific_compartments <- c(
    "liver1",
    "liver2",
    "liver3",
    "liver4",
    "liver5",
    "is_liver1_rif",
    "is_liver2_rif",
    "is_liver3_rif",
    "is_liver4_rif",
    "is_liver5_rif",
    "int_liver1_rif",
    "int_liver2_rif",
    "int_liver3_rif",
    "int_liver4_rif",
    "int_liver5_rif",
    "serosa_rif",
    "intestine_ent_rif",
    "intestine_muc_rif",
    "enzyme_ugt_liver1",
    "enzyme_ugt_liver2",
    "enzyme_ugt_liver3",
    "enzyme_ugt_liver4",
    "enzyme_ugt_liver5",
    "enzyme_ugt_gut",
    "enzyme_3a4_liver1",
    "enzyme_3a4_liver2",
    "enzyme_3a4_liver3",
    "enzyme_3a4_liver4",
    "enzyme_3a4_liver5",
    "enzyme_2c8_liver1",
    "enzyme_2c8_liver2",
    "enzyme_2c8_liver3",
    "enzyme_2c8_liver4",
    "enzyme_2c8_liver5"
  )

  # Rifampicin: doses in umol (MW 822.94 g/mol, 600 mg = 729 umol), states
  # in umol/L. Pioglitazone: doses in ug, states in ug/L = ng/mL.
  units <- list(
    time = "h",
    dosing = "umol (rifampicin); ug (pioglitazone)",
    concentration = "ng/mL (Cc, pioglitazone); umol/L (Cc_rif)"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Every blood flow, tissue volume and clearance of Table S1A / Table 1",
        "is per kg and is multiplied by WT, exactly as the Supplementary Model",
        "Code multiplies by BW. Table S1 footnote: body weight 'was in",
        "accordance with the original reported value or assumed to be 70 kg'."
      ),
      source_name = "BW"
    )
  )

  compartmentData <- list(
    central = list(analyte = "pioglitazone", units = "ng/mL", specimen = "whole blood", verified = TRUE),
    liver1 = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    liver2 = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    liver3 = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    liver4 = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    liver5 = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "pioglitazone", units = "ng/mL", specimen = "tissue", verified = TRUE),
    gut_lumen = list(analyte = "pioglitazone", units = "ug", specimen = "administration site", verified = TRUE),
    central_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "whole blood", verified = TRUE),
    is_liver1_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver2_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver3_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver4_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    is_liver5_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver1_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver2_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver3_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver4_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    int_liver5_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    serosa_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    muscle_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    skin_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    adipose_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    gut_lumen_rif = list(analyte = "rifampicin", units = "umol", specimen = "administration site", verified = TRUE),
    intestine_ent_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "tissue", verified = TRUE),
    intestine_muc_rif = list(analyte = "rifampicin", units = "umol/L", specimen = "whole blood", verified = TRUE),
    enzyme_ugt_liver1 = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_liver2 = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_liver3 = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_liver4 = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_liver5 = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_gut = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver1 = list(
      analyte = "CYP3A4 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver2 = list(
      analyte = "CYP3A4 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver3 = list(
      analyte = "CYP3A4 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver4 = list(
      analyte = "CYP3A4 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver5 = list(
      analyte = "CYP3A4 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c8_liver1 = list(
      analyte = "CYP2C8 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c8_liver2 = list(
      analyte = "CYP2C8 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c8_liver3 = list(
      analyte = "CYP2C8 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c8_liver4 = list(
      analyte = "CYP2C8 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c8_liver5 = list(
      analyte = "CYP2C8 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy volunteers",
    weight_median = "70 kg (assumed; Table S1 footnote)",
    dose_range = "rifampicin 600 mg PO once daily for 6 days; pioglitazone 30 mg single PO",
    notes = paste(
      "Jaakkola 2006 (ref. 16): pioglitazone 30 mg before and after rifampicin",
      "600 mg QD for 6 days; mean blood profiles of healthy volunteers used to",
      "estimate the pioglitazone parameters and the CYP2C8 Emax."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Physiological constants -- Supplementary Table S1A (per kg body
    # weight; multiplied by WT in model()). Identical in all five programs
    # of the Supplementary Model Code.
    # ------------------------------------------------------------------
    q_liver <- fixed(1.24)
    label("Hepatic blood flow (L/h/kg)") # Table S1A Liver 1.24; code Qh = 1.24 * BW
    q_muscle <- fixed(0.642)
    label("Muscle blood flow (L/h/kg)") # Table S1A Muscle 0.642; code Qm
    q_skin <- fixed(0.257)
    label("Skin blood flow (L/h/kg)") # Table S1A Skin 0.257; code Qs
    q_adipose <- fixed(0.223)
    label("Adipose blood flow (L/h/kg)") # Table S1A Adipose 0.223; code Qa
    q_serosa <- fixed(0.274)
    label("Intestinal serosa blood flow (L/h/kg)") # Table S1A Serosa 0.274; code Qser
    q_villi <- fixed(0.257)
    label("Intestinal mucosal (villous) blood flow (L/h/kg)") # Table S1A Mucosal blood 0.257; code Qvilli
    v_blood <- fixed(0.0743)
    label("Blood (central) volume (L/kg)") # Table S1A Blood 0.0743; code Vb = 0.0743 * BW
    v_hc <- fixed(0.0174)
    label("Hepatocyte volume, all five units (L/kg)") # Table S1A Hepatocytes 0.0174; code Vh
    v_he <- fixed(0.0067)
    label("Hepatic extracellular volume, all five units (L/kg)") # Table S1A Hepatic extracellular space 0.0067; code Vi
    v_muscle <- fixed(0.429)
    label("Muscle volume (L/kg)") # Table S1A Muscle 0.429; code Vm
    v_skin <- fixed(0.111)
    label("Skin volume (L/kg)") # Table S1A Skin 0.111; code Vs
    v_adipose <- fixed(0.143)
    label("Adipose volume (L/kg)") # Table S1A Adipose 0.143; code Va
    v_serosa <- fixed(0.00893)
    label("Intestinal serosa volume (L/kg)") # Table S1A Serosa 0.00893; code Vser
    v_mucb <- fixed(0.00099)
    label("Intestinal mucosal blood volume (L/kg)") # Table S1A Mucosal blood 0.00099; code Vmucb
    v_ent <- fixed(0.00739)
    label("Enterocyte volume (L/kg)") # Table S1A Enterocytes 0.00739 (Simcyp v14.1 default); code Vent

    # ------------------------------------------------------------------
    # Rifampicin (perpetrator) PBPK -- main-text Table 1 'Rifampicin'
    # column and its footnote; all carried from Asaumi 2018 (ref. 5) and
    # held fixed in this paper. Concentrations are in umol/L (the Km,u
    # and EC50,u of Table 1 are in uM), so rifampicin doses are in umol.
    # ------------------------------------------------------------------
    ltlag_rif <- fixed(log(0.255))
    label("Rifampicin lag time in intestinal absorption (h)") # Table 1 Rifampicin Tlag 0.255
    lkin_rif <- fixed(log(4.0))
    label("Rifampicin lumen-to-enterocyte absorption rate constant kin (1/h)") # Table 1 Rifampicin kin 4.0
    fa_rif <- fixed(1)
    label("Rifampicin fraction absorbed from gut lumen Fa (unitless)") # Table 1 footnote: Fa (1); code FaRif = 1
    fg_rif <- fixed(0.943)
    label("Rifampicin intestinal availability Fg (unitless)") # Table 1 footnote: Fg (0.943); Table 1 FaFg 0.943
    fb_rif <- fixed(0.0778)
    label("Rifampicin unbound fraction in blood (unitless)") # Table 1 Rifampicin fB 0.0778
    fh_rif <- fixed(0.0814)
    label("Rifampicin unbound fraction in hepatocytes (unitless)") # Table 1 Rifampicin fH 0.0814
    fent_rif <- fixed(0.115)
    label("Rifampicin unbound fraction in enterocytes fE (unitless)") # Table 1 footnote: fE (0.115)
    sfkp_rif <- fixed(6.65)
    label("Rifampicin common scaling factor on in silico Kp values (unitless)") # Table 1 Rifampicin SF Kp 6.65
    lkp_muscle_rif <- fixed(log(0.0947))
    label("Rifampicin muscle-to-blood partition coefficient (unitless)") # Table 1 Rifampicin Kp,muscle 0.0947
    lkp_skin_rif <- fixed(log(0.326))
    label("Rifampicin skin-to-blood partition coefficient (unitless)") # Table 1 Rifampicin Kp,skin 0.326
    lkp_adipose_rif <- fixed(log(0.0629))
    label("Rifampicin adipose-to-blood partition coefficient (unitless)") # Table 1 Rifampicin Kp,adipose 0.0629
    lkp_serosa_rif <- fixed(log(0.200))
    label("Rifampicin serosa-to-blood partition coefficient (unitless)") # Table 1 Rifampicin Kp,serosa 0.200
    fbclint_all_rif <- fixed(0.251)
    label("Rifampicin unbound-blood overall hepatic intrinsic clearance fB*CLint,all (L/h/kg)") # Table 1 Rifampicin fB CLint,all 0.251
    rdif_rif <- fixed(0.129)
    label("Rifampicin Rdif = PSdif,inf / PSact,inf (unitless)") # Table 1 Rifampicin Rdif 0.129
    beta_rif <- fixed(0.2)
    label("Rifampicin beta, rate-determining fraction of hepatic clearance (unitless)") # Table 1 Rifampicin beta 0.2
    gamma_rif <- fixed(0.778)
    label("Rifampicin gamma = PSdif,inf / PSdif,eff (unitless)") # Table 1 Rifampicin gamma 0.778
    km_oatp_rif <- fixed(0.177)
    label("Rifampicin unbound Michaelis-Menten constant for hepatic uptake Km,u (umol/L)") # Table 1 footnote: Km,u for hepatic uptake (0.177 uM)
    ps_dif_ent_rif <- fixed(0.14)
    label("Rifampicin passive diffusion clearance on the enterocyte basolateral membrane PSdif,E (L/h/kg)") # Table 1 footnote: PSdif,E (0.14 L/hour/kg)
    fm_ugt <- fixed(0.759)
    label("Rifampicin fraction metabolised by UGT, liver and enterocytes (unitless)") # Table 1 fm in liver 0.759 (UGT); footnote: fm via UGT in enterocytes (0.759)
    lcl_renal_rif <- fixed(log(0.011))
    label("Rifampicin renal clearance (L/h/kg)") # Table 1 Rifampicin CLrenal 0.011
    ec50u_rif <- fixed(0.0639)
    label("Rifampicin unbound EC50 for induction, shared by all enzymes and OATP1B (umol/L)") # Table 1 footnote: EC50,u (0.0639 uM); Methods: OATP1B EC50,u set equal to CYP3A
    emax_ugt <- fixed(1.34)
    label("Rifampicin maximum induction effect on UGT (unitless)") # Table 1 footnote: Emax for UGT autoinduction (1.34)
    emax_3a4 <- fixed(4.57)
    label("Rifampicin maximum induction effect on CYP3A, liver and enterocytes (unitless)") # Table 1 footnote: Emax for CYP3A (4.57)
    emax_2c8 <- 2.55
    label("Rifampicin maximum induction effect on CYP2C8 (unitless)") # Table S2B Emax for CYP2C8 2.55 +/- 0.32 (estimated in this paper, pioglitazone data)
    kdeg_ugt_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic UGT (1/h)") # Table S1B UGT hepatocytes 0.0158 (assumed equal to CYP3A4)
    kdeg_ugt_gut <- fixed(0.0288)
    label("Degradation rate constant of enterocyte UGT (1/h)") # Table S1B UGT enterocytes 0.0288 (assumed equal to CYP3A4)
    kdeg_3a4_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic CYP3A4 (1/h)") # Table S1B CYP3A4 hepatocytes 0.0158
    kdeg_2c8_liver <- fixed(0.0301)
    label("Degradation rate constant of hepatic CYP2C8 (1/h)") # Table S1B CYP2C8 hepatocytes 0.0301

    # ------------------------------------------------------------------
    # Pioglitazone (victim; CYP2C8 probe). Estimated in this paper: Table
    # S2A 'Pioglitazone' column (Table 1 marks these 'a Initial value').
    # The remainder is Table 1 'Pioglitazone' column. fB (0.03) is listed
    # in Table 1 but the model uses fB*CLint,met directly, so it is not a
    # separate parameter here.
    # ------------------------------------------------------------------
    ltlag <- log(0.87)
    label("Pioglitazone lag time in intestinal absorption (h)") # Table S2A Pioglitazone Tlag 0.87 +/- 1.50
    lka <- log(3.9)
    label("Pioglitazone absorption rate constant (1/h)") # Table S2A Pioglitazone ka 3.9 +/- 21.1
    sfkp <- 1.6
    label("Pioglitazone common scaling factor on in silico Kp values (unitless)") # Table S2A Pioglitazone SF Kp 1.6 +/- 0.15
    fbclint_met <- 0.050
    label("Pioglitazone unbound-blood hepatic metabolic intrinsic clearance fB*CLint,met (L/h/kg)") # Table S2A Pioglitazone fB CLint,met 0.050 +/- 0.003
    fafg <- fixed(0.854)
    label("Pioglitazone intestinal availability FaFg (unitless)") # Table 1 Pioglitazone FaFg 0.854 (from bioavailability 0.83)
    lkp_liver <- fixed(log(0.324))
    label("Pioglitazone liver-to-blood partition coefficient (unitless)") # Table 1 Pioglitazone Kp,liver 0.324 (in silico)
    lkp_muscle <- fixed(log(0.162))
    label("Pioglitazone muscle-to-blood partition coefficient (unitless)") # Table 1 Pioglitazone Kp,muscle 0.162
    lkp_skin <- fixed(log(0.585))
    label("Pioglitazone skin-to-blood partition coefficient (unitless)") # Table 1 Pioglitazone Kp,skin 0.585
    lkp_adipose <- fixed(log(0.818))
    label("Pioglitazone adipose-to-blood partition coefficient (unitless)") # Table 1 Pioglitazone Kp,adipose 0.818
    fm_cyp2c8 <- fixed(0.836)
    label("Pioglitazone fraction of hepatic metabolism by CYP2C8 (unitless)") # Table 1 Pioglitazone fm 0.836 (CYP2C8); Supplementary Text Eq. 1 with alpha 2.99
    fm_cyp3a4 <- fixed(0.164)
    label("Pioglitazone fraction of hepatic metabolism by CYP3A (unitless)") # Table 1 Pioglitazone fm 0.164 (CYP3A) = 1 - fm,CYP2C8
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology scaled to body weight (Supplementary Model Code,
    #    'preliminary calculating formula' column).
    # ------------------------------------------------------------------
    qh <- q_liver * WT
    qm <- q_muscle * WT
    qs <- q_skin * WT
    qa <- q_adipose * WT
    qser <- q_serosa * WT
    qvilli <- q_villi * WT
    qhart <- qh - qser - qvilli
    vb <- v_blood * WT
    vh <- v_hc * WT
    vi <- v_he * WT
    vm <- v_muscle * WT
    vs <- v_skin * WT
    va <- v_adipose * WT
    vser <- v_serosa * WT
    vmucb <- v_mucb * WT
    vent <- v_ent * WT

    # ------------------------------------------------------------------
    # 2. Rifampicin PBPK (Supplementary Text 'Case 1: rifampicin'; code
    #    states y51-y68). Hepatic clearance decomposed by the extended
    #    clearance concept ('Other equations'), with saturable OATP1B
    #    uptake; UGT autoinduction scales hepatic and enterocyte
    #    metabolism.
    # ------------------------------------------------------------------
    kin_rif <- exp(lkin_rif)
    cl_int_all_rif <- fbclint_all_rif * WT / fb_rif
    vmax_oatp_rif <- 1 / (1 + rdif_rif) * cl_int_all_rif / beta_rif * km_oatp_rif
    ps_dif_inf_rif <- rdif_rif / (1 + rdif_rif) * cl_int_all_rif / beta_rif
    ps_dif_eff_rif <- rdif_rif / (1 + rdif_rif) * cl_int_all_rif / beta_rif / gamma_rif
    cl_int_met_rif <- rdif_rif / (1 + rdif_rif) * cl_int_all_rif / (1 - beta_rif) / gamma_rif
    ps_ent_rif <- ps_dif_ent_rif * WT
    q_gut_rif <- fent_rif * ps_ent_rif * qvilli / (qvilli + fb_rif * ps_ent_rif)
    cl_met_ent_rif <- (q_gut_rif * (1 / fg_rif - 1) - (1 - fa_rif) * fent_rif * ps_ent_rif) / fent_rif
    clr_rif <- exp(lcl_renal_rif) * WT
    kp_muscle_rif <- exp(lkp_muscle_rif)
    kp_skin_rif <- exp(lkp_skin_rif)
    kp_adipose_rif <- exp(lkp_adipose_rif)
    kp_serosa_rif <- exp(lkp_serosa_rif)
    upt1_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * is_liver1_rif) + ps_dif_inf_rif) * is_liver1_rif
    eff1_rif <- fh_rif * ps_dif_eff_rif * int_liver1_rif
    met1_rif <- fh_rif * cl_int_met_rif * (1 + fm_ugt * (enzyme_ugt_liver1 - 1)) * int_liver1_rif
    upt2_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * is_liver2_rif) + ps_dif_inf_rif) * is_liver2_rif
    eff2_rif <- fh_rif * ps_dif_eff_rif * int_liver2_rif
    met2_rif <- fh_rif * cl_int_met_rif * (1 + fm_ugt * (enzyme_ugt_liver2 - 1)) * int_liver2_rif
    upt3_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * is_liver3_rif) + ps_dif_inf_rif) * is_liver3_rif
    eff3_rif <- fh_rif * ps_dif_eff_rif * int_liver3_rif
    met3_rif <- fh_rif * cl_int_met_rif * (1 + fm_ugt * (enzyme_ugt_liver3 - 1)) * int_liver3_rif
    upt4_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * is_liver4_rif) + ps_dif_inf_rif) * is_liver4_rif
    eff4_rif <- fh_rif * ps_dif_eff_rif * int_liver4_rif
    met4_rif <- fh_rif * cl_int_met_rif * (1 + fm_ugt * (enzyme_ugt_liver4 - 1)) * int_liver4_rif
    upt5_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * is_liver5_rif) + ps_dif_inf_rif) * is_liver5_rif
    eff5_rif <- fh_rif * ps_dif_eff_rif * int_liver5_rif
    met5_rif <- fh_rif * cl_int_met_rif * (1 + fm_ugt * (enzyme_ugt_liver5 - 1)) * int_liver5_rif

    d/dt(central_rif) <- (qh * is_liver5_rif +
      qm * (muscle_rif / (sfkp_rif * kp_muscle_rif) - central_rif) +
      qs * (skin_rif / (sfkp_rif * kp_skin_rif) - central_rif) +
      qa * (adipose_rif / (sfkp_rif * kp_adipose_rif) - central_rif) -
      (qhart + qser + qvilli) * central_rif - clr_rif * central_rif) / vb # code y51
    d/dt(is_liver1_rif) <- (qhart * central_rif + qvilli * intestine_muc_rif +
      qser * serosa_rif / (sfkp_rif * kp_serosa_rif) - qh * is_liver1_rif +
      (eff1_rif - upt1_rif) / 5) / (vi / 5) # code y52
    d/dt(int_liver1_rif) <- (upt1_rif - eff1_rif - met1_rif) / vh # code y53
    d/dt(is_liver2_rif) <- (qh * (is_liver1_rif - is_liver2_rif) + (eff2_rif - upt2_rif) / 5) / (vi / 5) # code y54
    d/dt(int_liver2_rif) <- (upt2_rif - eff2_rif - met2_rif) / vh # code y55
    d/dt(is_liver3_rif) <- (qh * (is_liver2_rif - is_liver3_rif) + (eff3_rif - upt3_rif) / 5) / (vi / 5) # code y56
    d/dt(int_liver3_rif) <- (upt3_rif - eff3_rif - met3_rif) / vh # code y57
    d/dt(is_liver4_rif) <- (qh * (is_liver3_rif - is_liver4_rif) + (eff4_rif - upt4_rif) / 5) / (vi / 5) # code y58
    d/dt(int_liver4_rif) <- (upt4_rif - eff4_rif - met4_rif) / vh # code y59
    d/dt(is_liver5_rif) <- (qh * (is_liver4_rif - is_liver5_rif) + (eff5_rif - upt5_rif) / 5) / (vi / 5) # code y60
    d/dt(int_liver5_rif) <- (upt5_rif - eff5_rif - met5_rif) / vh # code y61
    d/dt(serosa_rif) <- qser * (central_rif - serosa_rif / (sfkp_rif * kp_serosa_rif)) / vser # code y62
    d/dt(muscle_rif) <- qm * (central_rif - muscle_rif / (sfkp_rif * kp_muscle_rif)) / vm # code y63
    d/dt(skin_rif) <- qs * (central_rif - skin_rif / (sfkp_rif * kp_skin_rif)) / vs # code y64
    d/dt(adipose_rif) <- qa * (central_rif - adipose_rif / (sfkp_rif * kp_adipose_rif)) / va # code y65
    # Segregated-flow intestine: lumen amount (umol), enterocyte and
    # mucosal-blood concentrations (umol/L).
    d/dt(gut_lumen_rif) <- -kin_rif / fa_rif * gut_lumen_rif + fent_rif * ps_ent_rif * intestine_ent_rif # code y66
    d/dt(intestine_ent_rif) <- (kin_rif * gut_lumen_rif + fb_rif * ps_ent_rif * intestine_muc_rif -
      fent_rif * (2 * ps_ent_rif + cl_met_ent_rif * (1 + fm_ugt * (enzyme_ugt_gut - 1))) * intestine_ent_rif) / vent # code y67
    d/dt(intestine_muc_rif) <- (qvilli * (central_rif - intestine_muc_rif) + fent_rif * ps_ent_rif * intestine_ent_rif -
      fb_rif * ps_ent_rif * intestine_muc_rif) / vmucb # code y68

    # ------------------------------------------------------------------
    # 3. Induction turnover (Supplementary Text Eq. 2): each relative
    #    enzyme / transporter amount starts at 1 and is driven by the
    #    unbound rifampicin concentration in its own hepatocyte unit (or in
    #    the enterocyte for the gut isoforms).
    # ------------------------------------------------------------------
    d/dt(enzyme_ugt_liver1) <- kdeg_ugt_liver * (1 + emax_ugt * fh_rif * int_liver1_rif / (fh_rif * int_liver1_rif + ec50u_rif) - enzyme_ugt_liver1) # code y81
    d/dt(enzyme_ugt_liver2) <- kdeg_ugt_liver * (1 + emax_ugt * fh_rif * int_liver2_rif / (fh_rif * int_liver2_rif + ec50u_rif) - enzyme_ugt_liver2) # code y82
    d/dt(enzyme_ugt_liver3) <- kdeg_ugt_liver * (1 + emax_ugt * fh_rif * int_liver3_rif / (fh_rif * int_liver3_rif + ec50u_rif) - enzyme_ugt_liver3) # code y83
    d/dt(enzyme_ugt_liver4) <- kdeg_ugt_liver * (1 + emax_ugt * fh_rif * int_liver4_rif / (fh_rif * int_liver4_rif + ec50u_rif) - enzyme_ugt_liver4) # code y84
    d/dt(enzyme_ugt_liver5) <- kdeg_ugt_liver * (1 + emax_ugt * fh_rif * int_liver5_rif / (fh_rif * int_liver5_rif + ec50u_rif) - enzyme_ugt_liver5) # code y85
    d/dt(enzyme_ugt_gut) <- kdeg_ugt_gut * (1 + emax_ugt * fent_rif * intestine_ent_rif / (fent_rif * intestine_ent_rif + ec50u_rif) - enzyme_ugt_gut) # code y86
    d/dt(enzyme_3a4_liver1) <- kdeg_3a4_liver * (1 + emax_3a4 * fh_rif * int_liver1_rif / (fh_rif * int_liver1_rif + ec50u_rif) - enzyme_3a4_liver1) # code y101
    d/dt(enzyme_3a4_liver2) <- kdeg_3a4_liver * (1 + emax_3a4 * fh_rif * int_liver2_rif / (fh_rif * int_liver2_rif + ec50u_rif) - enzyme_3a4_liver2) # code y102
    d/dt(enzyme_3a4_liver3) <- kdeg_3a4_liver * (1 + emax_3a4 * fh_rif * int_liver3_rif / (fh_rif * int_liver3_rif + ec50u_rif) - enzyme_3a4_liver3) # code y103
    d/dt(enzyme_3a4_liver4) <- kdeg_3a4_liver * (1 + emax_3a4 * fh_rif * int_liver4_rif / (fh_rif * int_liver4_rif + ec50u_rif) - enzyme_3a4_liver4) # code y104
    d/dt(enzyme_3a4_liver5) <- kdeg_3a4_liver * (1 + emax_3a4 * fh_rif * int_liver5_rif / (fh_rif * int_liver5_rif + ec50u_rif) - enzyme_3a4_liver5) # code y105
    d/dt(enzyme_2c8_liver1) <- kdeg_2c8_liver * (1 + emax_2c8 * fh_rif * int_liver1_rif / (fh_rif * int_liver1_rif + ec50u_rif) - enzyme_2c8_liver1) # code y121
    d/dt(enzyme_2c8_liver2) <- kdeg_2c8_liver * (1 + emax_2c8 * fh_rif * int_liver2_rif / (fh_rif * int_liver2_rif + ec50u_rif) - enzyme_2c8_liver2) # code y122
    d/dt(enzyme_2c8_liver3) <- kdeg_2c8_liver * (1 + emax_2c8 * fh_rif * int_liver3_rif / (fh_rif * int_liver3_rif + ec50u_rif) - enzyme_2c8_liver3) # code y123
    d/dt(enzyme_2c8_liver4) <- kdeg_2c8_liver * (1 + emax_2c8 * fh_rif * int_liver4_rif / (fh_rif * int_liver4_rif + ec50u_rif) - enzyme_2c8_liver4) # code y124
    d/dt(enzyme_2c8_liver5) <- kdeg_2c8_liver * (1 + emax_2c8 * fh_rif * int_liver5_rif / (fh_rif * int_liver5_rif + ec50u_rif) - enzyme_2c8_liver5) # code y125

    enzyme_ugt_liver1(0) <- 1
    enzyme_ugt_liver2(0) <- 1
    enzyme_ugt_liver3(0) <- 1
    enzyme_ugt_liver4(0) <- 1
    enzyme_ugt_liver5(0) <- 1
    enzyme_ugt_gut(0) <- 1
    enzyme_3a4_liver1(0) <- 1
    enzyme_3a4_liver2(0) <- 1
    enzyme_3a4_liver3(0) <- 1
    enzyme_3a4_liver4(0) <- 1
    enzyme_3a4_liver5(0) <- 1
    enzyme_2c8_liver1(0) <- 1
    enzyme_2c8_liver2(0) <- 1
    enzyme_2c8_liver3(0) <- 1
    enzyme_2c8_liver4(0) <- 1
    enzyme_2c8_liver5(0) <- 1

    # ------------------------------------------------------------------
    # 4. Pioglitazone PBPK (Supplementary Text 'Case 2'; code states
    #    y1-y10). Perfusion-limited five-unit liver (each unit holds the
    #    whole-liver concentration, volume (Vh + Vi) / 5); hepatic
    #    CYP2C8 and CYP3A are induced by rifampicin (Eq. 5).
    #    The gut lumen is emptied at ka * FaFg and the same flux enters
    #    liver unit 1, exactly as printed in the Supplementary Text and
    #    the Supplementary Model Code (so FaFg scales the absorption rate
    #    and the whole dose is absorbed).
    # ------------------------------------------------------------------
    ka <- exp(lka)
    clr <- 0 # Table 1 Pioglitazone CLrenal 0 (urinary excretion not detected)
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_adipose <- exp(lkp_adipose)
    kp_liver <- exp(lkp_liver)
    vliv <- 0.2 * (vh + vi)
    kpl <- sfkp * kp_liver
    fbcl_met <- fbclint_met * WT
    ind1 <- 1 + fm_cyp3a4 * (enzyme_3a4_liver1 - 1) + fm_cyp2c8 * (enzyme_2c8_liver1 - 1)
    ind2 <- 1 + fm_cyp3a4 * (enzyme_3a4_liver2 - 1) + fm_cyp2c8 * (enzyme_2c8_liver2 - 1)
    ind3 <- 1 + fm_cyp3a4 * (enzyme_3a4_liver3 - 1) + fm_cyp2c8 * (enzyme_2c8_liver3 - 1)
    ind4 <- 1 + fm_cyp3a4 * (enzyme_3a4_liver4 - 1) + fm_cyp2c8 * (enzyme_2c8_liver4 - 1)
    ind5 <- 1 + fm_cyp3a4 * (enzyme_3a4_liver5 - 1) + fm_cyp2c8 * (enzyme_2c8_liver5 - 1)

    d/dt(central) <- (qh * liver5 / kpl - qh * central +
      qm * (muscle / (sfkp * kp_muscle) - central) +
      qs * (skin / (sfkp * kp_skin) - central) +
      qa * (adipose / (sfkp * kp_adipose) - central) - clr * central) / vb # code y1
    d/dt(liver1) <- (qh * central - qh * liver1 / kpl - fbcl_met / 5 * liver1 / kpl * ind1 +
      ka * fafg * gut_lumen) / vliv # code y2
    d/dt(liver2) <- (qh * (liver1 - liver2) - fbcl_met / 5 * liver2 * ind2) / kpl / vliv # code y3
    d/dt(liver3) <- (qh * (liver2 - liver3) - fbcl_met / 5 * liver3 * ind3) / kpl / vliv # code y4
    d/dt(liver4) <- (qh * (liver3 - liver4) - fbcl_met / 5 * liver4 * ind4) / kpl / vliv # code y5
    d/dt(liver5) <- (qh * (liver4 - liver5) - fbcl_met / 5 * liver5 * ind5) / kpl / vliv # code y6
    d/dt(muscle) <- qm * (central - muscle / (sfkp * kp_muscle)) / vm # code y7
    d/dt(skin) <- qs * (central - skin / (sfkp * kp_skin)) / vs # code y8
    d/dt(adipose) <- qa * (central - adipose / (sfkp * kp_adipose)) / va # code y9
    d/dt(gut_lumen) <- -ka * fafg * gut_lumen # code y10

    # Oral victim doses go into gut_lumen after the lag time.
    alag(gut_lumen) <- exp(ltlag)

    # Rifampicin oral doses go into gut_lumen_rif (umol) after the Table 1
    # lag time. Intravenous rifampicin (glibenclamide study) is dosed into
    # central_rif, a CONCENTRATION state, so the amount is divided by the
    # blood volume on the way in; give the infusion duration with `dur` in
    # the event table.
    alag(gut_lumen_rif) <- exp(ltlag_rif)
    f(central_rif) <- 1 / vb

    # ------------------------------------------------------------------
    # 5. Observations: every state is already a blood concentration. The
    #    paper fits and reports BLOOD concentrations. No residual-error
    #    model is published (Napp nonlinear least squares, weight = square
    #    root of the value), so the model is deterministic.
    # ------------------------------------------------------------------
    Cc <- central
    Cc_rif <- central_rif
  })
}
