Asaumi_2019_coproporphyrin_I_rifampicin_pbpk <- function() {
  description <- paste(
    "PBPK (five-unit tandem liver, extended clearance concept, segregated-flow",
    "intestine for rifampicin; Napp). Coproporphyrin I (CP-I), an endogenous",
    "OATP1B biomarker, under single and repeated oral rifampicin in healthy",
    "adults. CP-I is synthesised in the hepatocytes at a rate solved for a",
    "steady-state baseline, taken up by OATP1B, excreted into bile by MRP2 with",
    "three-compartment enterohepatic circulation and cleared renally. The",
    "rifampicin perpetrator model competitively inhibits OATP1B (in vivo Ki,u",
    "0.106 uM) and MRP2 (Ki,u 0.87 uM) and induces OATP1B by a turnover model",
    "(Emax 2.32, EC50,u 0.0639 uM, kdeg 0.0158 1/h). CP-I beta fixed at 0.2",
    "(parameters for 0.5 / 0.8 in the file comments). Deterministic: no IIV or",
    "residual error published."
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
    "CP-I Ki,u per beta and Mw from Yoshikado T et al. CPT Pharmacometrics Syst Pharmacol. 2018;7:739-747 (doi:10.1002/psp4.12348)."
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
    "ehc1",
    "ehc2",
    "ehc3",
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
    "oatp1b_liver1",
    "oatp1b_liver2",
    "oatp1b_liver3",
    "oatp1b_liver4",
    "oatp1b_liver5"
  )

  # Rifampicin: doses in umol (MW 822.94 g/mol, 600 mg = 729 umol), states
  # in umol/L. CP-I is endogenous (no dose); its states are ng/mL (ng in
  # the biliary and lumen amount states) and Cc is reported in nmol/L.
  units <- list(time = "h", dosing = "umol (rifampicin)", concentration = "nmol/L (Cc, CP-I); umol/L (Cc_rif)")

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
    central = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "whole blood", verified = TRUE),
    is_liver1 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver2 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver3 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver4 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    is_liver5 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver1 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "coproporphyrin I", units = "ng/mL", specimen = "tissue", verified = TRUE),
    ehc1 = list(analyte = "coproporphyrin I", units = "ng", specimen = "bile", verified = TRUE),
    ehc2 = list(analyte = "coproporphyrin I", units = "ng", specimen = "bile", verified = TRUE),
    ehc3 = list(analyte = "coproporphyrin I", units = "ng", specimen = "bile", verified = TRUE),
    gut_lumen = list(analyte = "coproporphyrin I", units = "ng", specimen = "administration site", verified = TRUE),
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
    dose_range = "rifampicin 300 or 600 mg single PO, or 600 mg PO once daily for 6 days",
    notes = paste(
      "CP-I profiles from Takehara 2018 (300 / 600 mg single dose; basis of the",
      "Yoshikado 2018 fit), Lai 2016 (600 mg single) and Kunze 2018 (600 mg",
      "single and 6-day repeated). Mean blood profiles of healthy volunteers;",
      "the predose CP-I concentration differs between studies (bl_cpi)."
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
    emax_oatp <- 2.32
    label("Rifampicin maximum induction effect on OATP1B (unitless)") # Table S2B Emax for OATP1B 2.32 +/- 0.18 at pravastatin beta = 0.2 (2.26 / 2.23 at beta = 0.5 / 0.8); main text 2.3
    kdeg_ugt_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic UGT (1/h)") # Table S1B UGT hepatocytes 0.0158 (assumed equal to CYP3A4)
    kdeg_ugt_gut <- fixed(0.0288)
    label("Degradation rate constant of enterocyte UGT (1/h)") # Table S1B UGT enterocytes 0.0288 (assumed equal to CYP3A4)
    kdeg_oatp_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic OATP1B (1/h)") # Table S1B OATP1B hepatocytes 0.0158 (assumed equal to CYP3A4)

    # ------------------------------------------------------------------
    # Coproporphyrin I (endogenous victim) -- Table 1 'Coproporphyrin I'
    # column; values are the beta = 0.2 set (Table 1 lists FaFg and
    # fB CLint,all as 0.318/0.309/0.256 and 0.453/0.402/0.369 for
    # beta = 0.2/0.5/0.8). All were optimised in Yoshikado 2018 (ref. 29,
    # Table 2) and are held fixed here.
    # ------------------------------------------------------------------
    lka <- fixed(log(3.0))
    label("CP-I reabsorption rate constant from the gut lumen (1/h)") # Table 1 Coproporphyrin I ka 3.0
    fafg <- fixed(0.318)
    label("CP-I intestinal availability of reabsorbed CP-I FaFg (unitless)") # Table 1 Coproporphyrin I FaFg 0.318 (beta 0.2)
    fbclint_all <- fixed(0.453)
    label("CP-I unbound-blood overall hepatic intrinsic clearance fB*CLint,all (L/h/kg)") # Table 1 Coproporphyrin I fB CLint,all 0.453 (beta 0.2); = 43.1 x 0.0105 (Yoshikado 2018 Table 2 C)
    lkbile <- fixed(log(5.2))
    label("CP-I transit rate constant in enterohepatic circulation (1/h)") # Table 1 Coproporphyrin I kbile 5.2
    sfkp <- fixed(1)
    label("CP-I common scaling factor on in silico Kp values (unitless)") # Table 1 Coproporphyrin I SF Kp 1
    lkp_muscle <- fixed(log(0.103))
    label("CP-I muscle-to-blood partition coefficient (unitless)") # Table 1 Coproporphyrin I Kp,muscle 0.103
    lkp_skin <- fixed(log(0.442))
    label("CP-I skin-to-blood partition coefficient (unitless)") # Table 1 Coproporphyrin I Kp,skin 0.442
    lkp_adipose <- fixed(log(0.079))
    label("CP-I adipose-to-blood partition coefficient (unitless)") # Table 1 Coproporphyrin I Kp,adipose 0.079
    fb <- fixed(0.0105)
    label("CP-I unbound fraction in blood (unitless)") # Table 1 Coproporphyrin I fB 0.0105
    fh <- fixed(0.0697)
    label("CP-I unbound fraction in hepatocytes (unitless)") # Table 1 Coproporphyrin I fH 0.0697
    rdif <- fixed(0.035)
    label("CP-I Rdif = PSdif,inf / PSact,inf (unitless)") # Table 1 Coproporphyrin I Rdif 0.035
    beta <- fixed(0.2)
    label("CP-I beta, rate-determining fraction of hepatic clearance (unitless)") # Table 1 Coproporphyrin I beta 0.2/0.5/0.8 (sensitivity values); 0.2 used here
    gamma <- fixed(0.020)
    label("CP-I gamma = PSdif,inf / PSdif,eff (unitless)") # Table 1 Coproporphyrin I gamma 0.020
    fbile <- fixed(0.84)
    label("CP-I fraction of hepatocyte elimination that is biliary excretion (unitless)") # Table 1 Coproporphyrin I fbile 0.84
    f_oatp1b <- fixed(1)
    label("CP-I fraction of active hepatic uptake mediated by OATP1B (unitless)") # Methods: 'f OATP1B values were set to be unity'
    lcl_renal <- fixed(log(0.0421))
    label("CP-I renal clearance (L/h/kg)") # Table 1 Coproporphyrin I CLrenal 0.0421
    ki_oatp_rif <- fixed(0.106)
    label("Rifampicin unbound Ki for OATP1B-mediated CP-I uptake (umol/L)") # Methods '~0.1 uM' (ref. 29); Yoshikado 2018 Table 2 C (beta 0.2) 0.106 (0.0965 / 0.0824 at beta 0.5 / 0.8)
    ki_mrp2_rif <- fixed(0.87)
    label("Rifampicin unbound Ki for MRP2-mediated biliary excretion of CP-I (umol/L)") # Methods: Ki,u for MRP2 (0.87 uM, ref. 29)
    mw_cpi <- fixed(654)
    label("CP-I molecular weight (g/mol)") # Supplementary Model Code y19 = y1 / Mw * 1000; Yoshikado 2018 Table S1 Mw 654
    bl_cpi <- fixed(0.49)
    label("Pre-dose (baseline) blood CP-I concentration (nmol/L)") # Figure 6 legend: predose CP-I set to 0.49 nM (b, Takehara), 0.58 nM (c, Lai), 0.45 nM (d, Kunze)
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
    oatp1b_liver1(0) <- 1
    oatp1b_liver2(0) <- 1
    oatp1b_liver3(0) <- 1
    oatp1b_liver4(0) <- 1
    oatp1b_liver5(0) <- 1

    # ------------------------------------------------------------------
    # 4. Coproporphyrin I PBPK (Supplementary Text 'Case 3'; code states y1-y18).
    #    OATP1B uptake is induced (oatp1b_liver<n>) and competitively
    #    inhibited by unbound extracellular rifampicin (Eq. 3).
    # ------------------------------------------------------------------
    ka <- exp(lka)
    kbile <- exp(lkbile)
    clr <- exp(lcl_renal) * WT
    kp_muscle <- exp(lkp_muscle)
    kp_skin <- exp(lkp_skin)
    kp_adipose <- exp(lkp_adipose)
    cl_int_all <- fbclint_all * WT / fb
    ps_act_inf <- 1 / (1 + rdif) * cl_int_all / beta
    ps_dif_inf <- rdif / (1 + rdif) * cl_int_all / beta
    ps_dif_eff <- rdif / (1 + rdif) * cl_int_all / beta / gamma
    cl_int_h <- rdif / (1 + rdif) * cl_int_all / (1 - beta) / gamma
    cl_int_bile <- cl_int_h * fbile
    cl_int_met <- cl_int_h * (1 - fbile)
    upt1 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver1 - 1)) / (1 + fb_rif * is_liver1_rif / ki_oatp_rif) + ps_dif_inf) * is_liver1
    eff1 <- fh * ps_dif_eff * int_liver1
    mrp2_1 <- 1 / (1 + fh_rif * int_liver1_rif / ki_mrp2_rif)
    bile1 <- fh * cl_int_bile * mrp2_1 * int_liver1
    elim1 <- bile1 + fh * cl_int_met * int_liver1
    upt2 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver2 - 1)) / (1 + fb_rif * is_liver2_rif / ki_oatp_rif) + ps_dif_inf) * is_liver2
    eff2 <- fh * ps_dif_eff * int_liver2
    mrp2_2 <- 1 / (1 + fh_rif * int_liver2_rif / ki_mrp2_rif)
    bile2 <- fh * cl_int_bile * mrp2_2 * int_liver2
    elim2 <- bile2 + fh * cl_int_met * int_liver2
    upt3 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver3 - 1)) / (1 + fb_rif * is_liver3_rif / ki_oatp_rif) + ps_dif_inf) * is_liver3
    eff3 <- fh * ps_dif_eff * int_liver3
    mrp2_3 <- 1 / (1 + fh_rif * int_liver3_rif / ki_mrp2_rif)
    bile3 <- fh * cl_int_bile * mrp2_3 * int_liver3
    elim3 <- bile3 + fh * cl_int_met * int_liver3
    upt4 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver4 - 1)) / (1 + fb_rif * is_liver4_rif / ki_oatp_rif) + ps_dif_inf) * is_liver4
    eff4 <- fh * ps_dif_eff * int_liver4
    mrp2_4 <- 1 / (1 + fh_rif * int_liver4_rif / ki_mrp2_rif)
    bile4 <- fh * cl_int_bile * mrp2_4 * int_liver4
    elim4 <- bile4 + fh * cl_int_met * int_liver4
    upt5 <- fb * (ps_act_inf * (1 + f_oatp1b * (oatp1b_liver5 - 1)) / (1 + fb_rif * is_liver5_rif / ki_oatp_rif) + ps_dif_inf) * is_liver5
    eff5 <- fh * ps_dif_eff * int_liver5
    mrp2_5 <- 1 / (1 + fh_rif * int_liver5_rif / ki_mrp2_rif)
    bile5 <- fh * cl_int_bile * mrp2_5 * int_liver5
    elim5 <- bile5 + fh * cl_int_met * int_liver5

    # ------------------------------------------------------------------
    # 4b. Steady-state biosynthesis and initial conditions, transcribed
    #    from the CP-I program of the Supplementary Model Code (the
    #    'preliminary calculating formula' column). The synthesis rate
    #    vsyn (ng/h, total over the five hepatocyte units) is solved so
    #    that the drug-free system sits at the baseline blood
    #    concentration c1; the states are ng/mL internally.
    # ------------------------------------------------------------------
    c1 <- bl_cpi * mw_cpi / 1000
    ps_inf <- ps_act_inf + ps_dif_inf
    ss_a <- fb * ps_inf / fh / (ps_dif_eff + cl_int_h)
    ss_b <- (qh + 0.2 * fb * ps_inf) / qh - 0.2 * fh * ps_dif_eff * ss_a / qh
    ss_x <- 0.2 * fh * cl_int_bile
    ss_y <- -5 * qh * ss_x / (0.2 * fh * ps_dif_eff) + 4 * ss_a * ss_x + 3 * ss_a * ss_x * ss_b +
      2 * ss_a * ss_x * ss_b^2 + ss_a * ss_x * ss_b^3 -
      (qh + 0.2 * fb * ps_inf) / fafg * (1 + ss_b + ss_b^2 + ss_b^3) +
      0.2 * fh * ps_dif_eff / fafg * (-qh / (0.2 * fh * ps_dif_eff) + ss_a + ss_a * ss_b + ss_a * ss_b^2 + ss_a * ss_b^3)
    ss_z <- (qh + clr) / qh * c1 * ((qh + 0.2 * fb * ps_inf) / fafg * ss_b^4 - qh^2 / fafg / (qh + clr) -
      0.2 * fh * ps_dif_eff / fafg * ss_a * ss_b^4 - ss_a * ss_x * (1 + ss_b + ss_b^2 + ss_b^3 + ss_b^4))
    vsyn <- -qh * fh * (ps_dif_eff + cl_int_h) / (0.2 * fh * ps_dif_eff) * ss_z / ss_y
    x18_0 <- (vsyn - clr * c1) / (ka / fafg / fbile - ka)
    x15_0 <- ka / fafg * x18_0 / kbile
    ss_in <- qh * c1 + ka * x18_0 + vsyn / 5
    kph_0 <- (vsyn * qh / ss_in + fb * ps_inf) / (fh * (ps_dif_eff + cl_int_h) + vsyn / 5 * fh * cl_int_h / ss_in)
    c2_0 <- ss_in / (qh + fh * cl_int_h / 5 * kph_0)
    c4_0 <- (qh * c2_0 + vsyn / 5) / (qh + fh * cl_int_h / 5 * kph_0)
    c6_0 <- (qh * c4_0 + vsyn / 5) / (qh + fh * cl_int_h / 5 * kph_0)
    c8_0 <- (qh * c6_0 + vsyn / 5) / (qh + fh * cl_int_h / 5 * kph_0)
    c10_0 <- (qh * c8_0 + vsyn / 5) / (qh + fh * cl_int_h / 5 * kph_0)
    central(0) <- c1
    is_liver1(0) <- c2_0
    int_liver1(0) <- c2_0 * kph_0
    is_liver2(0) <- c4_0
    int_liver2(0) <- c4_0 * kph_0
    is_liver3(0) <- c6_0
    int_liver3(0) <- c6_0 * kph_0
    is_liver4(0) <- c8_0
    int_liver4(0) <- c8_0 * kph_0
    is_liver5(0) <- c10_0
    int_liver5(0) <- c10_0 * kph_0
    muscle(0) <- c1 * kp_muscle * sfkp
    skin(0) <- c1 * kp_skin * sfkp
    adipose(0) <- c1 * kp_adipose * sfkp
    ehc1(0) <- x15_0
    ehc2(0) <- x15_0
    ehc3(0) <- x15_0
    gut_lumen(0) <- x18_0

    d/dt(central) <- (qh * (is_liver5 - central) +
      qm * (muscle / (sfkp * kp_muscle) - central) +
      qs * (skin / (sfkp * kp_skin) - central) +
      qa * (adipose / (sfkp * kp_adipose) - central) - clr * central) / vb # code y1
    d/dt(is_liver1) <- (qh * central - qh * is_liver1 + (eff1 - upt1) / 5 + ka * gut_lumen) / (vi / 5) # code y2
    d/dt(int_liver1) <- (vsyn + upt1 - eff1 - elim1) / vh # code y3
    d/dt(is_liver2) <- (qh * (is_liver1 - is_liver2) + (eff2 - upt2) / 5) / (vi / 5) # code y4
    d/dt(int_liver2) <- (vsyn + upt2 - eff2 - elim2) / vh # code y5
    d/dt(is_liver3) <- (qh * (is_liver2 - is_liver3) + (eff3 - upt3) / 5) / (vi / 5) # code y6
    d/dt(int_liver3) <- (vsyn + upt3 - eff3 - elim3) / vh # code y7
    d/dt(is_liver4) <- (qh * (is_liver3 - is_liver4) + (eff4 - upt4) / 5) / (vi / 5) # code y8
    d/dt(int_liver4) <- (vsyn + upt4 - eff4 - elim4) / vh # code y9
    d/dt(is_liver5) <- (qh * (is_liver4 - is_liver5) + (eff5 - upt5) / 5) / (vi / 5) # code y10
    d/dt(int_liver5) <- (vsyn + upt5 - eff5 - elim5) / vh # code y11
    d/dt(muscle) <- qm * (central - muscle / (sfkp * kp_muscle)) / vm # code y12
    d/dt(skin) <- qs * (central - skin / (sfkp * kp_skin)) / vs # code y13
    d/dt(adipose) <- qa * (central - adipose / (sfkp * kp_adipose)) / va # code y14
    # Enterohepatic circulation: three biliary transit compartments and
    # the gut lumen (amounts), reabsorbed into hepatic extracellular unit 1.
    d/dt(ehc1) <- (bile1 + bile2 + bile3 + bile4 + bile5) / 5 - kbile * ehc1 # code y15
    d/dt(ehc2) <- kbile * (ehc1 - ehc2) # code y16
    d/dt(ehc3) <- kbile * (ehc2 - ehc3) # code y17
    d/dt(gut_lumen) <- kbile * ehc3 - ka / fafg * gut_lumen # code y18

    # Rifampicin oral doses go into gut_lumen_rif (umol) after the Table 1
    # lag time. Intravenous rifampicin (glibenclamide study) is dosed into
    # central_rif, a CONCENTRATION state, so the amount is divided by the
    # blood volume on the way in; give the infusion duration with `dur` in
    # the event table.
    alag(gut_lumen_rif) <- exp(ltlag_rif)
    f(central_rif) <- 1 / vb

    # ------------------------------------------------------------------
    # 5. Observations. CP-I states are ng/mL; the paper reports blood CP-I
    #    in nmol/L (code y19 = y1 / Mw * 1000). No residual-error model is
    #    published, so the model is deterministic.
    # ------------------------------------------------------------------
    Cc <- central / mw_cpi * 1000
    Cc_rif <- central_rif

    # Relative transporter activities in hepatocyte unit 1 (Figure 2c):
    # OATP1B = induction x competitive inhibition, MRP2 = inhibition only.
    oatp1b_activity <- oatp1b_liver1 / (1 + fb_rif * is_liver1_rif / ki_oatp_rif)
    mrp2_activity <- mrp2_1
  })
}
