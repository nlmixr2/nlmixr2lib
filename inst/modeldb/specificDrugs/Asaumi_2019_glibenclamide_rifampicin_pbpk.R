Asaumi_2019_glibenclamide_rifampicin_pbpk <- function() {
  description <- paste(
    "PBPK (five-unit tandem liver, extended clearance concept, segregated-flow",
    "intestine for rifampicin; Napp). Rifampicin-glibenclamide (glyburide)",
    "drug-drug interaction in healthy adults, with intravenous and repeated",
    "oral rifampicin. The rifampicin perpetrator model (saturable OATP1B",
    "uptake, UGT autoinduction) is coupled to a glibenclamide victim model with",
    "a Qgut model of intestinal CYP3A first-pass; rifampicin competitively",
    "inhibits OATP1B uptake (in vivo Ki,u 0.13 uM) and induces OATP1B (Emax",
    "2.32), hepatic CYP2C9 (Emax 2.41) and hepatic and enterocyte CYP3A (Emax",
    "4.57) by turnover models (EC50,u 0.0639 uM). Glibenclamide hepatic",
    "metabolism is 85% CYP2C9 and 15% CYP3A at the beta = 0.2 used here.",
    "Deterministic: no IIV or residual error published."
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
    "Glibenclamide fm by beta from Asaumi R et al. CPT Pharmacometrics Syst Pharmacol. 2018;7:186-196 (doi:10.1002/psp4.12275)."
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
    "is_liver1",
    "is_liver2",
    "is_liver3",
    "is_liver4",
    "is_liver5",
    "int_liver1",
    "int_liver2",
    "int_liver3",
    "int_liver4",
    "int_liver5",
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
    "enzyme_3a4_gut",
    "enzyme_2c9_liver1",
    "enzyme_2c9_liver2",
    "enzyme_2c9_liver3",
    "enzyme_2c9_liver4",
    "enzyme_2c9_liver5",
    "oatp1b_liver1",
    "oatp1b_liver2",
    "oatp1b_liver3",
    "oatp1b_liver4",
    "oatp1b_liver5"
  )

  # Rifampicin: doses in umol (MW 822.94 g/mol, 600 mg = 729 umol), states
  # in umol/L. Glibenclamide: doses in ug, states in ug/L = ng/mL.
  units <- list(
    time = "h",
    dosing = "umol (rifampicin); ug (glibenclamide)",
    concentration = "ng/mL (Cc, glibenclamide); umol/L (Cc_rif)"
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
    central = list(analyte = "glibenclamide", units = "ng/mL", specimen = "whole blood", verified = TRUE),
    is_liver1 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver2 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver3 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver4 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver5 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver1 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "glibenclamide", units = "ng/mL", specimen = "tissue", verified = TRUE),
    portal = list(analyte = "glibenclamide", units = "ng/mL", specimen = "whole blood", verified = TRUE),
    gut_lumen = list(analyte = "glibenclamide", units = "ug", specimen = "administration site", verified = TRUE),
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
    enzyme_3a4_gut = list(
      analyte = "CYP3A4 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c9_liver1 = list(
      analyte = "CYP2C9 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c9_liver2 = list(
      analyte = "CYP2C9 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c9_liver3 = list(
      analyte = "CYP2C9 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c9_liver4 = list(
      analyte = "CYP2C9 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_2c9_liver5 = list(
      analyte = "CYP2C9 relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    oatp1b_liver1 = list(
      analyte = "OATP1B relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    oatp1b_liver2 = list(
      analyte = "OATP1B relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    oatp1b_liver3 = list(
      analyte = "OATP1B relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    oatp1b_liver4 = list(
      analyte = "OATP1B relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    oatp1b_liver5 = list(
      analyte = "OATP1B relative amount (induced / control)",
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
    dose_range = "rifampicin 600 mg single IV, or 600 mg PO once daily for 6-7 days (+/- IV on day 7); glibenclamide 1.25 mg single PO",
    notes = paste(
      "Clinical DDI AUC ratios of Zheng 2009 (ref. 19): glibenclamide with a",
      "single IV 600 mg rifampicin, with IV rifampicin after 6 days of oral",
      "600 mg QD, and 48 h after 7 days of oral 600 mg QD. Observed profiles",
      "were not published; healthy volunteers."
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
    q_portal <- fixed(0.531)
    label("Portal vein blood flow (L/h/kg)") # Table S1A Portal vein 0.531; code Qpv (glibenclamide program only)
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
    v_portal <- fixed(0.001)
    label("Portal vein volume (L/kg)") # Table S1A Portal vein 0.001; code Vpv (glibenclamide program only)

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
    emax_2c9 <- fixed(2.41)
    label("Rifampicin maximum induction effect on CYP2C9 (unitless)") # Table 1 footnote: Emax for CYP2C9 (2.41)
    emax_oatp <- 2.32
    label("Rifampicin maximum induction effect on OATP1B (unitless)") # Table S2B Emax for OATP1B 2.32 +/- 0.18 at pravastatin beta = 0.2 (2.26 / 2.23 at beta = 0.5 / 0.8); main text 2.3
    kdeg_ugt_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic UGT (1/h)") # Table S1B UGT hepatocytes 0.0158 (assumed equal to CYP3A4)
    kdeg_ugt_gut <- fixed(0.0288)
    label("Degradation rate constant of enterocyte UGT (1/h)") # Table S1B UGT enterocytes 0.0288 (assumed equal to CYP3A4)
    kdeg_3a4_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic CYP3A4 (1/h)") # Table S1B CYP3A4 hepatocytes 0.0158
    kdeg_3a4_gut <- fixed(0.0288)
    label("Degradation rate constant of enterocyte CYP3A4 (1/h)") # Table S1B CYP3A4 enterocytes 0.0288
    kdeg_2c9_liver <- fixed(0.00666)
    label("Degradation rate constant of hepatic CYP2C9 (1/h)") # Table S1B CYP2C9 hepatocytes 0.00666
    kdeg_oatp_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic OATP1B (1/h)") # Table S1B OATP1B hepatocytes 0.0158 (assumed equal to CYP3A4)

    # ------------------------------------------------------------------
    # Glibenclamide (victim) -- Table 1 'Glibenclamide' column and its
    # footnote; carried from Asaumi 2018 (ref. 5) and held fixed here.
    # ------------------------------------------------------------------
    ltlag <- fixed(log(0.773))
    label("Glibenclamide lag time in intestinal absorption (h)") # Table 1 Glibenclamide Tlag 0.773
    lka <- fixed(log(0.445))
    label("Glibenclamide absorption rate constant (1/h)") # Table 1 Glibenclamide ka 0.445
    fa <- fixed(1)
    label("Glibenclamide fraction absorbed Fa (unitless)") # Table 1 footnote: Fa (1)
    sfkp <- fixed(0.57)
    label("Glibenclamide common scaling factor on in silico Kp values (unitless)") # Table 1 Glibenclamide SF Kp 0.57
    lkp_muscle <- fixed(log(0.104))
    label("Glibenclamide muscle-to-blood partition coefficient (unitless)") # Table 1 Glibenclamide Kp,muscle 0.104
    lkp_skin <- fixed(log(0.447))
    label("Glibenclamide skin-to-blood partition coefficient (unitless)") # Table 1 Glibenclamide Kp,skin 0.447
    lkp_adipose <- fixed(log(0.0795))
    label("Glibenclamide adipose-to-blood partition coefficient (unitless)") # Table 1 Glibenclamide Kp,adipose 0.0795
    fbclint_all <- fixed(0.123)
    label("Glibenclamide unbound-blood overall hepatic intrinsic clearance fB*CLint,all (L/h/kg)") # Table 1 Glibenclamide fB CLint,all 0.123
    fb <- fixed(0.000774)
    label("Glibenclamide unbound fraction in blood (unitless)") # Table 1 Glibenclamide fB 0.000774
    fh <- fixed(0.0221)
    label("Glibenclamide unbound fraction in hepatocytes (unitless)") # Table 1 Glibenclamide fH 0.0221
    rdif <- fixed(0.246)
    label("Glibenclamide Rdif = PSdif,inf / PSact,inf (unitless)") # Table 1 Glibenclamide Rdif 0.246
    beta <- fixed(0.2)
    label("Glibenclamide beta, rate-determining fraction of hepatic clearance (unitless)") # Table 1 Glibenclamide beta 0.2/0.5/0.8 (sensitivity values); 0.2 used here
    gamma <- fixed(0.24)
    label("Glibenclamide gamma = PSdif,inf / PSdif,eff (unitless)") # Table 1 Glibenclamide gamma 0.24
    f_oatp1b <- fixed(1)
    label("Glibenclamide fraction of active hepatic uptake mediated by OATP1B (unitless)") # Methods: 'f OATP1B values were set to be unity'
    fm_cyp2c9 <- fixed(0.85)
    label("Glibenclamide fraction of hepatic metabolism by CYP2C9 (unitless)") # Table 1 fm 0.85-1 (CYP2C9); Asaumi 2018 Discussion: 0.85 / 0.94 / 1.00 at beta 0.2 / 0.5 / 0.8
    fm_cyp3a4 <- fixed(0.15)
    label("Glibenclamide fraction of hepatic metabolism by CYP3A (unitless)") # Table 1 fm 0-0.15 (CYP3A); Asaumi 2018 Discussion: 0.15 / 0.06 / 0.0 at beta 0.2 / 0.5 / 0.8
    fecl_int_ent <- fixed(0.0131)
    label("Glibenclamide unbound enterocyte intrinsic clearance fE*CLint,E (L/h/kg)") # Table 1 footnote: fE CLint,E (0.0131 L/hour/kg)
    cl_perm_ent <- fixed(0.103)
    label("Glibenclamide enterocyte permeability clearance for the Qgut model (L/h/kg)") # Table 1 footnote: clearance permeability in enterocytes (0.103 L/hour/kg)
    frac_ent_cyp3a4 <- fixed(1)
    label("Glibenclamide fraction of enterocyte metabolism by CYP3A (unitless)") # Table 1 footnote: fm via CYP3A pathway in enterocytes (1)
    ki_oatp_rif <- fixed(0.13)
    label("Rifampicin unbound Ki for OATP1B-mediated glibenclamide uptake (umol/L)") # Methods: 0.13 uM = 0.44 x 0.19 / 0.65
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
    qpv <- q_portal * WT
    vb <- v_blood * WT
    vh <- v_hc * WT
    vi <- v_he * WT
    vm <- v_muscle * WT
    vs <- v_skin * WT
    va <- v_adipose * WT
    vser <- v_serosa * WT
    vmucb <- v_mucb * WT
    vent <- v_ent * WT
    vpv <- v_portal * WT

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
    d/dt(enzyme_3a4_gut) <- kdeg_3a4_gut * (1 + emax_3a4 * fent_rif * intestine_ent_rif / (fent_rif * intestine_ent_rif + ec50u_rif) - enzyme_3a4_gut) # code y106
    d/dt(enzyme_2c9_liver1) <- kdeg_2c9_liver * (1 + emax_2c9 * fh_rif * int_liver1_rif / (fh_rif * int_liver1_rif + ec50u_rif) - enzyme_2c9_liver1) # code y111
    d/dt(enzyme_2c9_liver2) <- kdeg_2c9_liver * (1 + emax_2c9 * fh_rif * int_liver2_rif / (fh_rif * int_liver2_rif + ec50u_rif) - enzyme_2c9_liver2) # code y112
    d/dt(enzyme_2c9_liver3) <- kdeg_2c9_liver * (1 + emax_2c9 * fh_rif * int_liver3_rif / (fh_rif * int_liver3_rif + ec50u_rif) - enzyme_2c9_liver3) # code y113
    d/dt(enzyme_2c9_liver4) <- kdeg_2c9_liver * (1 + emax_2c9 * fh_rif * int_liver4_rif / (fh_rif * int_liver4_rif + ec50u_rif) - enzyme_2c9_liver4) # code y114
    d/dt(enzyme_2c9_liver5) <- kdeg_2c9_liver * (1 + emax_2c9 * fh_rif * int_liver5_rif / (fh_rif * int_liver5_rif + ec50u_rif) - enzyme_2c9_liver5) # code y115
    d/dt(oatp1b_liver1) <- kdeg_oatp_liver * (1 + emax_oatp * fh_rif * int_liver1_rif / (fh_rif * int_liver1_rif + ec50u_rif) - oatp1b_liver1) # code y141
    d/dt(oatp1b_liver2) <- kdeg_oatp_liver * (1 + emax_oatp * fh_rif * int_liver2_rif / (fh_rif * int_liver2_rif + ec50u_rif) - oatp1b_liver2) # code y142
    d/dt(oatp1b_liver3) <- kdeg_oatp_liver * (1 + emax_oatp * fh_rif * int_liver3_rif / (fh_rif * int_liver3_rif + ec50u_rif) - oatp1b_liver3) # code y143
    d/dt(oatp1b_liver4) <- kdeg_oatp_liver * (1 + emax_oatp * fh_rif * int_liver4_rif / (fh_rif * int_liver4_rif + ec50u_rif) - oatp1b_liver4) # code y144
    d/dt(oatp1b_liver5) <- kdeg_oatp_liver * (1 + emax_oatp * fh_rif * int_liver5_rif / (fh_rif * int_liver5_rif + ec50u_rif) - oatp1b_liver5) # code y145

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
    enzyme_3a4_gut(0) <- 1
    enzyme_2c9_liver1(0) <- 1
    enzyme_2c9_liver2(0) <- 1
    enzyme_2c9_liver3(0) <- 1
    enzyme_2c9_liver4(0) <- 1
    enzyme_2c9_liver5(0) <- 1
    oatp1b_liver1(0) <- 1
    oatp1b_liver2(0) <- 1
    oatp1b_liver3(0) <- 1
    oatp1b_liver4(0) <- 1
    oatp1b_liver5(0) <- 1

    # ------------------------------------------------------------------
    # 4. Glibenclamide PBPK (Supplementary Text 'Case 4'; code states
    #    y1-y16). Oral absorption through a portal-vein compartment with a
    #    Qgut model of CYP3A intestinal first-pass (induced by
    #    enterocyte rifampicin); OATP1B uptake is induced and
    #    competitively inhibited; hepatic CYP2C9 and CYP3A are induced.
    # ------------------------------------------------------------------
    ka <- exp(lka)
    clr <- 0 # Table 1 Glibenclamide CLrenal 0
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_adipose <- exp(lkp_adipose)
    cl_perm <- cl_perm_ent * WT
    q_gut <- qvilli * cl_perm / (qvilli + cl_perm)
    fecl_ent <- fecl_int_ent * WT
    cl_int_all <- fbclint_all * WT / fb
    ps_act_inf <- 1 / (1 + rdif) * cl_int_all / beta
    ps_dif_inf <- rdif / (1 + rdif) * cl_int_all / beta
    ps_dif_eff <- rdif / (1 + rdif) * cl_int_all / beta / gamma
    cl_int_h <- rdif / (1 + rdif) * cl_int_all / (1 - beta) / gamma
    upt1 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver1 - 1)) / (1 + fb_rif * is_liver1_rif / ki_oatp_rif) + ps_dif_inf) * is_liver1
    eff1 <- fh * ps_dif_eff * int_liver1
    elim1 <- fh * cl_int_h * (1 + fm_cyp3a4 * (enzyme_3a4_liver1 - 1) + fm_cyp2c9 * (enzyme_2c9_liver1 - 1)) * int_liver1
    upt2 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver2 - 1)) / (1 + fb_rif * is_liver2_rif / ki_oatp_rif) + ps_dif_inf) * is_liver2
    eff2 <- fh * ps_dif_eff * int_liver2
    elim2 <- fh * cl_int_h * (1 + fm_cyp3a4 * (enzyme_3a4_liver2 - 1) + fm_cyp2c9 * (enzyme_2c9_liver2 - 1)) * int_liver2
    upt3 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver3 - 1)) / (1 + fb_rif * is_liver3_rif / ki_oatp_rif) + ps_dif_inf) * is_liver3
    eff3 <- fh * ps_dif_eff * int_liver3
    elim3 <- fh * cl_int_h * (1 + fm_cyp3a4 * (enzyme_3a4_liver3 - 1) + fm_cyp2c9 * (enzyme_2c9_liver3 - 1)) * int_liver3
    upt4 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver4 - 1)) / (1 + fb_rif * is_liver4_rif / ki_oatp_rif) + ps_dif_inf) * is_liver4
    eff4 <- fh * ps_dif_eff * int_liver4
    elim4 <- fh * cl_int_h * (1 + fm_cyp3a4 * (enzyme_3a4_liver4 - 1) + fm_cyp2c9 * (enzyme_2c9_liver4 - 1)) * int_liver4
    upt5 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver5 - 1)) / (1 + fb_rif * is_liver5_rif / ki_oatp_rif) + ps_dif_inf) * is_liver5
    eff5 <- fh * ps_dif_eff * int_liver5
    elim5 <- fh * cl_int_h * (1 + fm_cyp3a4 * (enzyme_3a4_liver5 - 1) + fm_cyp2c9 * (enzyme_2c9_liver5 - 1)) * int_liver5

    d/dt(central) <- (qh * is_liver5 - (qh - qpv) * central - qpv * central +
      qm * (muscle / (sfkp * kp_muscle) - central) +
      qs * (skin / (sfkp * kp_skin) - central) +
      qa * (adipose / (sfkp * kp_adipose) - central) - clr * central) / vb # code y1
    d/dt(is_liver1) <- ((qh - qpv) * central + qpv * portal - qh * is_liver1 + (eff1 - upt1) / 5) / (vi / 5) # code y2
    d/dt(int_liver1) <- (upt1 - eff1 - elim1) / vh # code y3
    d/dt(is_liver2) <- (qh * (is_liver1 - is_liver2) + (eff2 - upt2) / 5) / (vi / 5) # code y4
    d/dt(int_liver2) <- (upt2 - eff2 - elim2) / vh # code y5
    d/dt(is_liver3) <- (qh * (is_liver2 - is_liver3) + (eff3 - upt3) / 5) / (vi / 5) # code y6
    d/dt(int_liver3) <- (upt3 - eff3 - elim3) / vh # code y7
    d/dt(is_liver4) <- (qh * (is_liver3 - is_liver4) + (eff4 - upt4) / 5) / (vi / 5) # code y8
    d/dt(int_liver4) <- (upt4 - eff4 - elim4) / vh # code y9
    d/dt(is_liver5) <- (qh * (is_liver4 - is_liver5) + (eff5 - upt5) / 5) / (vi / 5) # code y10
    d/dt(int_liver5) <- (upt5 - eff5 - elim5) / vh # code y11
    d/dt(muscle) <- qm * (central - muscle / (sfkp * kp_muscle)) / vm # code y12
    d/dt(skin) <- qs * (central - skin / (sfkp * kp_skin)) / vs # code y13
    d/dt(adipose) <- qa * (central - adipose / (sfkp * kp_adipose)) / va # code y14
    d/dt(portal) <- (qpv * (central - portal) +
      ka * q_gut / (q_gut + fecl_ent * (1 + frac_ent_cyp3a4 * (enzyme_3a4_gut - 1))) * gut_lumen) / vpv # code y15
    d/dt(gut_lumen) <- -ka / fa * gut_lumen # code y16

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
