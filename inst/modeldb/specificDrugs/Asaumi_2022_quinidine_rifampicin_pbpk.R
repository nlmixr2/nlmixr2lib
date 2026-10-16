Asaumi_2022_quinidine_rifampicin_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body; five-unit tandem liver, three-segment segregated-flow",
    "intestine, nephron-resolved kidney; Napp). Rifampicin-quinidine",
    "P-glycoprotein (P-gp) and CYP3A drug-drug interaction in healthy adults.",
    "The rifampicin perpetrator model (saturable OATP1B hepatic uptake, UGT",
    "autoinduction) induces P-gp (Emax 4.01) and CYP3A (Emax 4.57) by",
    "turnover models (EC50,u 0.0639 uM) and inhibits P-gp competitively",
    "(Ki,u 402 ng/mL) in the gut enterocytes and the renal proximal tubule.",
    "The quinidine victim model has a stomach, a PERFUSION-limited liver with",
    "saturable CYP3A metabolism (Km,u 1996 ng/mL), saturable intestinal CYP3A",
    "metabolism and saturable P-gp efflux (Km,u 409 ng/mL) into the gut lumen",
    "and the proximal-tubule lumen; it has no enterohepatic circulation. The",
    "study-specific quinidine parameter sets (Maeda 2011 default; Twum-Barima",
    "1981 IV and PO and Damkier 1999 via STUDY_ indicators) are switched",
    "inside the model. Deterministic: no IIV or residual error published."
  )
  reference <- paste(
    "Asaumi R, Nunoya K, Yamaura Y, Taskar KS, Sugiyama Y. Robust",
    "physiologically based pharmacokinetic model of rifampicin for predicting",
    "drug-drug interactions via P-glycoprotein induction and inhibition in",
    "the intestine, liver, and kidney. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11:919-933. doi:10.1002/psp4.12807 (PMC9286720).",
    "ODEs: Supplementary Material S1 ('Case 1: rifampicin', 'Case 4:",
    "quinidine', which lists only its differences from 'Case 2: digoxin').",
    "Physiology: Table S1; fitted parameters: Table S2."
  )
  vignette <- "Asaumi_2022_rifampicin_pgp_ddi_pbpk"

  # State-name conventions as in Asaumi_2022_talinolol_rifampicin_pbpk.R.
  # The quinidine liver is perfusion-limited, so each `liver<n>` unit is one
  # well-stirred compartment of volume (V_HE + V_HC) / 5 holding a tissue
  # concentration, rather than a separate `is_liver` / `int_liver` pair.
  # `enzyme_3a4_*` and `pgp_*` hold the relative CYP3A / P-gp amount
  # (induced / control).
  paper_specific_compartments <- c(
    "liver1",
    "liver2",
    "liver3",
    "liver4",
    "liver5",
    "serosa",
    "duodenum_muc",
    "jejunum_muc",
    "ileum_muc",
    "duodenum_ent",
    "jejunum_ent",
    "ileum_ent",
    "duodenum_lumen",
    "jejunum_lumen",
    "ileum_lumen",
    "a_feces",
    "glom_b",
    "glom_u",
    "pt_b",
    "pt_c",
    "pt_u",
    "dt_b",
    "dt_c",
    "dt_u",
    "cd_b",
    "cd_c",
    "cd_u",
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
    "duodenum_muc_rif",
    "jejunum_muc_rif",
    "ileum_muc_rif",
    "duodenum_ent_rif",
    "jejunum_ent_rif",
    "ileum_ent_rif",
    "duodenum_lumen_rif",
    "jejunum_lumen_rif",
    "ileum_lumen_rif",
    "enzyme_ugt_liver1",
    "enzyme_ugt_liver2",
    "enzyme_ugt_liver3",
    "enzyme_ugt_liver4",
    "enzyme_ugt_liver5",
    "enzyme_ugt_duodenum",
    "enzyme_ugt_jejunum",
    "enzyme_ugt_ileum",
    "enzyme_3a4_liver1",
    "enzyme_3a4_liver2",
    "enzyme_3a4_liver3",
    "enzyme_3a4_liver4",
    "enzyme_3a4_liver5",
    "enzyme_3a4_duodenum",
    "enzyme_3a4_jejunum",
    "enzyme_3a4_ileum",
    "pgp_duodenum",
    "pgp_jejunum",
    "pgp_ileum",
    "pgp_kidney"
  )

  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Every blood flow, tissue volume and clearance of Table S1 and Table 1",
        "is per kg body weight and is multiplied by WT. Table S1A footnote:",
        "body weight 'was in accordance with the original reported value or",
        "assumed to be 70 kg'."
      ),
      source_name = "body weight"
    ),
    STUDY_TWUM1981_IV = list(
      description = "1 = Twum-Barima 1981 intravenous quinidine arm (350 mg); 0 = Maeda 2011.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table S2B fits SF_Kp and fB*CLint,all per study; selects the Twum 1981 infusion set.",
      source_name = "Table S2B column 'Twum et al. (1981) (infusion)'"
    ),
    STUDY_TWUM1981_PO = list(
      description = "1 = Twum-Barima 1981 oral quinidine arm (350 mg); 0 = Maeda 2011.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table S2B fits SF_Kp and fB*CLint,all per study; selects the Twum 1981 oral set.",
      source_name = "Table S2B column 'Twum et al. (1981) (p.o.)'"
    ),
    STUDY_DAMKIER1999 = list(
      description = "1 = Damkier 1999 oral quinidine study (166 mg); 0 = Maeda 2011.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table S2B fits SF_Kp and fB*CLint,all per study; selects the Damkier 1999 set.",
      source_name = "Table S2B column 'Damkier et al. (1999)'"
    )
  )

  compartmentData <- list(
    central = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    stomach = list(analyte = "quinidine", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum_lumen = list(analyte = "quinidine", units = "ug", specimen = "administration site", verified = TRUE),
    jejunum_lumen = list(analyte = "quinidine", units = "ug", specimen = "administration site", verified = TRUE),
    ileum_lumen = list(analyte = "quinidine", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum_ent = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum_ent = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    ileum_ent = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    duodenum_muc = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    jejunum_muc = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    ileum_muc = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    a_feces = list(analyte = "quinidine", units = "ug", specimen = "faeces", verified = TRUE),
    serosa = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver1 = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver2 = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver3 = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver4 = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    liver5 = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    glom_b = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    glom_u = list(analyte = "quinidine", units = "ug", specimen = "urine", verified = TRUE),
    pt_b = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    pt_c = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    pt_u = list(analyte = "quinidine", units = "ug", specimen = "urine", verified = TRUE),
    dt_b = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    dt_c = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    dt_u = list(analyte = "quinidine", units = "ug", specimen = "urine", verified = TRUE),
    cd_b = list(analyte = "quinidine", units = "ug", specimen = "whole blood", verified = TRUE),
    cd_c = list(analyte = "quinidine", units = "ug", specimen = "tissue", verified = TRUE),
    cd_u = list(analyte = "quinidine", units = "ug", specimen = "urine", verified = TRUE),
    a_urine = list(analyte = "quinidine", units = "ug", specimen = "urine", verified = TRUE),
    central_rif = list(analyte = "rifampicin", units = "ug", specimen = "whole blood", verified = TRUE),
    duodenum_lumen_rif = list(analyte = "rifampicin", units = "ug", specimen = "administration site", verified = TRUE),
    jejunum_lumen_rif = list(analyte = "rifampicin", units = "ug", specimen = "administration site", verified = TRUE),
    ileum_lumen_rif = list(analyte = "rifampicin", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum_ent_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum_ent_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    ileum_ent_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    duodenum_muc_rif = list(analyte = "rifampicin", units = "ug", specimen = "whole blood", verified = TRUE),
    jejunum_muc_rif = list(analyte = "rifampicin", units = "ug", specimen = "whole blood", verified = TRUE),
    ileum_muc_rif = list(analyte = "rifampicin", units = "ug", specimen = "whole blood", verified = TRUE),
    serosa_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    muscle_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    skin_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    adipose_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver1_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver2_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver3_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver4_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver5_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver1_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver2_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver3_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver4_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver5_rif = list(analyte = "rifampicin", units = "ug", specimen = "tissue", verified = TRUE),
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
    enzyme_ugt_duodenum = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_jejunum = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_ugt_ileum = list(
      analyte = "UGT relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver1 = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver2 = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver3 = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver4 = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_liver5 = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_duodenum = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_jejunum = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    enzyme_3a4_ileum = list(
      analyte = "CYP3A relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_duodenum = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_jejunum = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_ileum = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_kidney = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    disease_state = "healthy volunteers",
    weight_median = "70 kg (assumed unless reported; Table S1A footnote)",
    dose_range = paste(
      "rifampicin 600 mg PO once daily for 7 or 11 days; quinidine 350 mg IV",
      "(30-min infusion) or PO 24 h after the last rifampicin dose (Twum-Barima",
      "1981), 166 mg PO simultaneous with rifampicin (Damkier 1999), or",
      "0.1 to 100 mg PO alone (Maeda 2011)"
    ),
    notes = paste(
      "Mean blood profiles digitised from the literature. Maeda 2011 (ref. 8):",
      "single oral quinidine microdoses to 100 mg for the nonlinear (P-gp and",
      "CYP3A saturation) fit; Twum-Barima 1981 (ref. 9) and Damkier 1999",
      "(ref. 10): quinidine before and after rifampicin treatment."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Physiology, per kg body weight (Table S1A-E). Shared with rifampicin.
    # ------------------------------------------------------------------
    q_liver <- fixed(1.24)
    label("Hepatic blood flow (L/h/kg)") # Table S1A Liver (total) 1.24
    q_hart <- fixed(0.310)
    label("Hepatic arterial blood flow (L/h/kg)") # Table S1A Hepatic arterial 0.310
    q_serosa <- fixed(0.744)
    label("Intestinal serosa blood flow (L/h/kg)") # Table S1A Serosa 0.744
    q_muc_duo <- fixed(0.0164)
    label("Duodenum mucosal blood flow (L/h/kg)") # Table S1A Mucosal blood, duodenum part 0.0164
    q_muc_jej <- fixed(0.0900)
    label("Jejunum mucosal blood flow (L/h/kg)") # Table S1A Mucosal blood, jejunum part 0.0900
    q_muc_ile <- fixed(0.0796)
    label("Ileum mucosal blood flow (L/h/kg)") # Table S1A Mucosal blood, ileum part 0.0796
    q_muscle <- fixed(0.642)
    label("Muscle blood flow (L/h/kg)") # Table S1A Muscle 0.642
    q_skin <- fixed(0.257)
    label("Skin blood flow (L/h/kg)") # Table S1A Skin 0.257
    q_adipose <- fixed(0.223)
    label("Adipose blood flow (L/h/kg)") # Table S1A Adipose 0.223
    v_blood <- fixed(0.0743)
    label("Blood (central) volume (L/kg)") # Table S1A Blood 0.0743
    v_hc <- fixed(0.0174)
    label("Hepatocyte volume, all five units (L/kg)") # Table S1A Hepatocytes 0.0174
    v_he <- fixed(0.0067)
    label("Hepatic extracellular volume, all five units (L/kg)") # Table S1A Hepatic extracellular space 0.0067
    v_muscle <- fixed(0.429)
    label("Muscle volume (L/kg)") # Table S1A Muscle 0.429
    v_skin <- fixed(0.111)
    label("Skin volume (L/kg)") # Table S1A Skin 0.111
    v_adipose <- fixed(0.143)
    label("Adipose volume (L/kg)") # Table S1A Adipose 0.143
    v_serosa <- fixed(0.00702)
    label("Intestinal serosa volume (L/kg)") # Table S1A Serosa 0.00702
    v_muc_duo <- fixed(0.000113)
    label("Duodenum mucosal blood volume (L/kg)") # Table S1A Mucosal blood, duodenum part 0.000113
    v_muc_jej <- fixed(0.000400)
    label("Jejunum mucosal blood volume (L/kg)") # Table S1A Mucosal blood, jejunum part 0.000400
    v_muc_ile <- fixed(0.000477)
    label("Ileum mucosal blood volume (L/kg)") # Table S1A Mucosal blood, ileum part 0.000477
    v_ent_duo <- fixed(0.000842)
    label("Duodenum enterocyte volume (L/kg)") # Table S1A Enterocytes, duodenum part 0.000842
    v_ent_jej <- fixed(0.00299)
    label("Jejunum enterocyte volume (L/kg)") # Table S1A Enterocytes, jejunum part 0.00299
    v_ent_ile <- fixed(0.00356)
    label("Ileum enterocyte volume (L/kg)") # Table S1A Enterocytes, ileum part 0.00356
    v_lumen_duo <- fixed(0.000935)
    label("Duodenum lumen volume (L/kg)") # Table S1A Gut lumen, duodenum part 0.000935
    v_lumen_jej <- fixed(0.00331)
    label("Jejunum lumen volume (L/kg)") # Table S1A Gut lumen, jejunum part 0.00331
    v_lumen_ile <- fixed(0.00395)
    label("Ileum lumen volume (L/kg)") # Table S1A Gut lumen, ileum part 0.00395
    k_feces_duo <- fixed(5.13)
    label("Duodenum lumen transit rate constant (1/h)") # Table S1B k_feces duodenum part 5.13
    k_feces_jej <- fixed(1.47)
    label("Jejunum lumen transit rate constant (1/h)") # Table S1B k_feces jejunum part 1.47
    k_feces_ile <- fixed(1.21)
    label("Ileum lumen transit rate constant to faeces (1/h)") # Table S1B k_feces ileum part 1.21
    fsa_duo <- fixed(0.076)
    label("Duodenum fraction of intestinal surface area (fraction)") # Table S1C surface area 7.6%; Suppl. S1 x 0.076
    fsa_jej <- fixed(0.394)
    label("Jejunum fraction of intestinal surface area (fraction)") # Table S1C surface area 39.4%; Suppl. S1 x 0.394
    fsa_ile <- fixed(0.530)
    label("Ileum fraction of intestinal surface area (fraction)") # Table S1C surface area 53.0%; Suppl. S1 x 0.530
    fcyp_duo <- fixed(0.138)
    label("Duodenum fraction of intestinal CYP3A metabolism (fraction)") # Suppl. S1 Case 4: CL_met,CYP3A,duodenum x 0.138
    fcyp_jej <- fixed(0.545)
    label("Jejunum fraction of intestinal CYP3A metabolism (fraction)") # Suppl. S1 Case 4: CL_met,CYP3A,jejunum x 0.545
    fcyp_ile <- fixed(0.317)
    label("Ileum fraction of intestinal CYP3A metabolism (fraction)") # Suppl. S1 Case 4: CL_met,CYP3A,ileum x 0.317
    q_renal <- fixed(1.06)
    label("Renal blood flow (L/h/kg)") # Table S1D Kidney 1.06
    gfr <- fixed(0.107)
    label("Glomerular filtration rate (L/h/kg)") # Table S1D Glomerulus urinary lumen 0.107
    q_glom_pt <- fixed(0.956)
    label("Glomerulus-to-proximal vessel blood flow (L/h/kg)") # Table S1D Glomerulus blood vessels 0.956 (outlet)
    q_pt_dt <- fixed(0.991)
    label("Proximal-to-distal vessel blood flow (L/h/kg)") # Table S1D Proximal tubule blood vessels 0.991 (outlet)
    q_dt_cd <- fixed(1.03)
    label("Distal-to-collecting-duct vessel blood flow (L/h/kg)") # Table S1D Distal tubule blood vessels 1.03 (outlet)
    qu_pt_dt <- fixed(0.0370)
    label("Proximal-to-distal tubular fluid flow (L/h/kg)") # Table S1D Proximal tubule urinary lumen 0.0370 (outlet)
    qu_dt_cd <- fixed(0.00994)
    label("Distal-to-collecting-duct tubular fluid flow (L/h/kg)") # Table S1D Distal tubule urinary lumen 0.00994 (outlet)
    qu_cd <- fixed(0.000833)
    label("Urine flow out of the collecting duct (L/h/kg)") # Table S1D Collecting duct urinary lumen 0.000833 (outlet)
    v_glom_b <- fixed(0.0000465)
    label("Glomerulus blood-vessel volume (L/kg)") # Table S1D glomerulus 0.0000465
    v_glom_u <- fixed(0.0000713)
    label("Glomerulus urinary-lumen volume (L/kg)") # Table S1D glomerulus 0.0000713
    v_pt_b <- fixed(0.000474)
    label("Proximal tubule blood-vessel volume (L/kg)") # Table S1D proximal tubule 0.000474
    v_pt_c <- fixed(0.000395)
    label("Proximal tubule cell volume (L/kg)") # Table S1D proximal tubule 0.000395
    v_pt_u <- fixed(0.000727)
    label("Proximal tubule urinary-lumen volume (L/kg)") # Table S1D proximal tubule 0.000727
    v_dt_b <- fixed(0.000102)
    label("Distal tubule blood-vessel volume (L/kg)") # Table S1D distal tubule 0.000102
    v_dt_c <- fixed(0.0000852)
    label("Distal tubule cell volume (L/kg)") # Table S1D distal tubule 0.0000852
    v_dt_u <- fixed(0.000157)
    label("Distal tubule urinary-lumen volume (L/kg)") # Table S1D distal tubule 0.000157
    v_cd_b <- fixed(0.000577)
    label("Collecting duct blood-vessel volume (L/kg)") # Table S1D collecting duct 0.000577
    v_cd_c <- fixed(0.000480)
    label("Collecting duct cell volume (L/kg)") # Table S1D collecting duct 0.000480
    v_cd_u <- fixed(0.000884)
    label("Collecting duct urinary-lumen volume (L/kg)") # Table S1D collecting duct 0.000884
    kdeg_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic UGT, CYP3A and P-gp (1/h)") # Table S1E hepatic k_deg 0.0158
    kdeg_gut <- fixed(0.0288)
    label("Degradation rate constant of intestinal UGT, CYP3A and P-gp (1/h)") # Table S1E intestinal k_deg 0.0288
    kdeg_pgp_kidney <- fixed(0.0158)
    label("Degradation rate constant of renal P-gp (1/h)") # Table S1E renal P-gp k_deg 0.0158 (assumed = hepatic CYP3A4)

    # ------------------------------------------------------------------
    # Rifampicin (perpetrator). Table 1 'Rifampicin' column and note.
    # ------------------------------------------------------------------
    ltlag_rif <- fixed(log(0.255))
    label("Rifampicin lag time in intestinal absorption (h)") # Table 1 Rifampicin T_lag 0.255
    fb_rif <- fixed(0.0778)
    label("Rifampicin unbound fraction in blood (fraction)") # Table 1 Rifampicin f_B 0.0778
    fh_rif <- fixed(0.0814)
    label("Rifampicin unbound fraction in hepatocytes (fraction)") # Table 1 Rifampicin f_hepatocytes 0.0814
    fent_rif <- fixed(0.115)
    label("Rifampicin unbound fraction in enterocytes (fraction)") # Table 1 Rifampicin f_enterocytes 0.115
    sfkp_rif <- fixed(6.65)
    label("Rifampicin common scaling factor on in silico Kp values (unitless)") # Table 1 Rifampicin SF_Kp 6.65
    lkp_muscle_rif <- fixed(log(0.0947))
    label("Rifampicin muscle-to-blood partition coefficient (unitless)") # Table 1 Rifampicin K_p,muscle 0.0947
    lkp_skin_rif <- fixed(log(0.326))
    label("Rifampicin skin-to-blood partition coefficient (unitless)") # Table 1 Rifampicin K_p,skin 0.326
    lkp_adipose_rif <- fixed(log(0.0629))
    label("Rifampicin adipose-to-blood partition coefficient (unitless)") # Table 1 Rifampicin K_p,adipose 0.0629
    lkp_serosa_rif <- fixed(log(0.200))
    label("Rifampicin serosa-to-blood partition coefficient (unitless)") # Table 1 Rifampicin K_p,serosa 0.200
    fbclint_all_rif <- fixed(0.251)
    label("Rifampicin unbound-blood overall hepatic intrinsic clearance fB*CLint,all (L/h/kg)") # Table 1 Rifampicin f_B CL_int,all 0.251
    rdif_rif <- fixed(0.129)
    label("Rifampicin Rdif,inf,liver (unitless)") # Table 1 Rifampicin R_dif,inf,liver 0.129
    beta_rif <- fixed(0.2)
    label("Rifampicin hepatic beta (unitless)") # Table 1 Rifampicin beta_liver 0.2
    gamma_rif <- fixed(0.778)
    label("Rifampicin gamma_inf,liver (unitless)") # Table 1 Rifampicin gamma_inf,liver 0.778
    fm_ugt <- fixed(0.759)
    label("Rifampicin fraction metabolised by UGT (fraction)") # Table 1 Rifampicin f_m 0.759 (UGT)
    km_oatp_rif <- fixed(146)
    label("Rifampicin unbound Km for OATP1B uptake (ng/mL)") # Table 1 Rifampicin K_m,u 146 (OATP1B)
    ps_dif_gut_rif <- fixed(0.08)
    label("Rifampicin passive gut-lumen-to-enterocyte clearance, whole intestine (L/h/kg)") # Table 1 Rifampicin PS_dif,gut lumen to enterocytes 0.08
    gamma_ent_rif <- fixed(0.778)
    label("Rifampicin enterocyte gamma_inf = gamma_eff (unitless)") # Table 1 Rifampicin gamma_inf,enterocytes 0.778
    cl_met_ent_rif <- fixed(0.005)
    label("Rifampicin enterocyte metabolic intrinsic clearance, whole intestine (L/h/kg)") # Table 1 Rifampicin CL_int,met,enterocytes 0.005
    lcl_renal_rif <- fixed(log(0.011))
    label("Rifampicin renal clearance (L/h/kg)") # Table 1 note: rifampicin CL_renal (0.011 L/h/kg)
    ec50u_rif <- fixed(0.0639)
    label("Rifampicin unbound EC50 for induction (umol/L)") # Table 1 note EC50,u 63.9 nmol/L
    emax_ugt <- fixed(1.34)
    label("Rifampicin maximum induction effect on UGT (unitless)") # Table 1 note: E_max,UGT (1.34)
    emax_cyp3a <- fixed(4.57)
    label("Rifampicin maximum induction effect on CYP3A (unitless)") # Table 1 note: E_max,CYP3A (4.57)
    emax_pgp <- fixed(4.01)
    label("Rifampicin maximum induction effect on P-gp (unitless)") # Table S2C E_max for P-gp 4.01
    ki_pgp_rif <- fixed(402)
    label("Rifampicin unbound inhibition constant for P-gp (ng/mL)") # Table S2C K_i for P-gp 402 ng/mL (0.488 umol/L)

    # ------------------------------------------------------------------
    # Quinidine (victim). Table 1 'Quinidine' column and note; SF_Kp and
    # fB*CLint,all are study-specific (Table S2B).
    # ------------------------------------------------------------------
    lkstomach <- log(5)
    label("Quinidine stomach-to-duodenum transit rate constant (1/h)") # Table 1 Quinidine k_stomach 5 (initial); Table S2B fitted 6.0
    fb <- fixed(0.218)
    label("Quinidine unbound fraction in blood (fraction)") # Table 1 Quinidine f_B 0.218
    rb <- fixed(0.919)
    label("Quinidine blood-to-plasma concentration ratio (ratio)") # Table 1 Quinidine R_B 0.919
    fent <- fixed(0.215)
    label("Quinidine unbound fraction in enterocytes (fraction)") # Table 1 Quinidine f_enterocytes 0.215
    frenal <- fixed(0.856)
    label("Quinidine unbound fraction in renal cells (fraction)") # Table 1 Quinidine f_renal cell 0.856
    sfkp <- 6.1
    label("Quinidine common scaling factor on in silico Kp values, Maeda 2011 (unitless)") # Table S2B Maeda (2011) SF_Kp 6.1
    sfkp_twiv <- 5.5
    label("Quinidine SF_Kp, Twum 1981 infusion (unitless)") # Table S2B Twum (1981) (infusion) SF_Kp 5.5
    sfkp_twpo <- 1.2
    label("Quinidine SF_Kp, Twum 1981 oral (unitless)") # Table S2B Twum (1981) (p.o.) SF_Kp 1.2
    sfkp_dam <- 6.2
    label("Quinidine SF_Kp, Damkier 1999 (unitless)") # Table S2B Damkier (1999) SF_Kp 6.2
    lkp_liver <- fixed(log(0.409))
    label("Quinidine liver-to-blood partition coefficient (unitless)") # Table 1 Quinidine K_p,liver 0.409
    lkp_muscle <- fixed(log(0.434))
    label("Quinidine muscle-to-blood partition coefficient (unitless)") # Table 1 Quinidine K_p,muscle 0.434
    lkp_skin <- fixed(log(0.589))
    label("Quinidine skin-to-blood partition coefficient (unitless)") # Table 1 Quinidine K_p,skin 0.589
    lkp_adipose <- fixed(log(0.229))
    label("Quinidine adipose-to-blood partition coefficient (unitless)") # Table 1 Quinidine K_p,adipose 0.229
    lkp_serosa <- fixed(log(0.523))
    label("Quinidine serosa-to-blood partition coefficient (unitless)") # Table 1 Quinidine K_p,serosa 0.523
    fbclint_all <- 0.35
    label("Quinidine unbound-blood hepatic intrinsic clearance, Maeda 2011 (L/h/kg)") # Table S2B Maeda (2011) f_B CL_int 0.35
    fbclint_twiv <- 0.17
    label("Quinidine fB*CLint, Twum 1981 infusion (L/h/kg)") # Table S2B Twum (1981) (infusion) f_B CL_int 0.17
    fbclint_twpo <- 0.12
    label("Quinidine fB*CLint, Twum 1981 oral (L/h/kg)") # Table S2B Twum (1981) (p.o.) f_B CL_int 0.12
    fbclint_dam <- 0.14
    label("Quinidine fB*CLint, Damkier 1999 (L/h/kg)") # Table S2B Damkier (1999) f_B CL_int 0.14
    fm_cyp3a4 <- fixed(0.962)
    label("Quinidine fraction of hepatic metabolism by CYP3A (fraction)") # Table 1 Quinidine f_m,liver 0.962 (CYP3A); enterocyte f_m,CYP3A 1 (all gut metabolism is CYP3A)
    ps_dif_gut <- fixed(0.011)
    label("Quinidine passive gut-lumen-to-enterocyte clearance, whole intestine (L/h/kg)") # Table 1 Quinidine PS_dif,gut lumen to enterocytes 0.011 (footnote h)
    cl_met_ent <- fixed(0.008)
    label("Quinidine enterocyte metabolic intrinsic clearance, whole intestine (L/h/kg)") # Table 1 Quinidine CL_int,met,enterocytes 0.008 (footnote h)
    rdif_eff_duo <- fixed(1.38)
    label("Quinidine Rdif,eff in the duodenum (unitless)") # Table 1 Quinidine R_dif,eff,duodenum 1.38
    rdif_eff_jej <- fixed(0.385)
    label("Quinidine Rdif,eff in the jejunum (unitless)") # Table 1 Quinidine R_dif,eff,jejunum 0.385
    rdif_eff_ile <- fixed(0.315)
    label("Quinidine Rdif,eff in the ileum (unitless)") # Table 1 Quinidine R_dif,eff,ileum 0.315
    gamma_ent <- fixed(1)
    label("Quinidine enterocyte gamma_inf = gamma_eff (unitless)") # Table 1 Quinidine gamma_inf,enterocytes 1
    fpgp_ent <- fixed(1)
    label("Quinidine P-gp fraction of enterocyte apical efflux (fraction)") # Table 1 Quinidine f_P-gp,enterocytes 1
    rdif_inf_kidney <- fixed(0.02)
    label("Quinidine Rdif,inf,kidney (unitless)") # Table 1 Quinidine R_dif,inf,kidney 0.02 (footnote k)
    rdif_eff_kidney <- fixed(0.3)
    label("Quinidine Rdif,eff,kidney (unitless)") # Table 1 Quinidine R_dif,eff,kidney 0.3 (footnote k)
    fpgp_kidney <- fixed(1)
    label("Quinidine P-gp fraction of proximal-tubule apical efflux (fraction)") # Table 1 Quinidine f_P-gp,kidney 1
    ps_dif_pc_v <- fixed(0.0104)
    label("Quinidine passive proximal cell-to-vessel clearance (L/h/kg)") # Table 1 Quinidine PS_dif,proximal cells to vessels 0.0104
    ps_dif_pc_t <- fixed(0.0784)
    label("Quinidine passive proximal cell-to-tubule clearance (L/h/kg)") # Table 1 Quinidine PS_dif,proximal cells to tubule 0.0784
    ps_dif_dc <- fixed(0.0027)
    label("Quinidine passive distal cell-to-vessel = cell-to-tubule clearance (L/h/kg)") # Table 1 Quinidine PS_dif,distal cells to vessels 0.0027 (footnote m)
    ps_dif_cc <- fixed(0.000579)
    label("Quinidine passive collecting-duct cell-to-vessel = cell-to-tubule clearance (L/h/kg)") # Table 1 Quinidine PS_dif,collecting duct cells to vessels 0.000579 (footnote n)
    sf_renal <- fixed(0.1)
    label("Quinidine scaling factor on renal passive permeability (unitless)") # Table 1 note: quinidine SF_renal permeability (0.1)
    gamma_inf_pc <- fixed(7.22)
    label("Quinidine gamma_inf in proximal cells (unitless)") # Table 1 note: quinidine gamma_inf,proximal cells (7.22)
    gamma_inf_dc <- fixed(5.96)
    label("Quinidine gamma_inf in distal cells (unitless)") # Table 1 note: quinidine gamma_inf,distal cells (5.96)
    gamma_inf_cc <- fixed(7.22)
    label("Quinidine gamma_inf in collecting-duct cells (unitless)") # Table 1 note: quinidine gamma_inf,collecting duct cells (7.22)
    gamma_eff_pc <- fixed(6.41)
    label("Quinidine gamma_eff in proximal cells (unitless)") # Table 1 note: quinidine gamma_eff,proximal cells (6.41)
    gamma_eff_dc <- fixed(4.35)
    label("Quinidine gamma_eff in distal cells (unitless)") # Table 1 note: quinidine gamma_eff,distal cells (4.35)
    gamma_eff_cc <- fixed(2.50)
    label("Quinidine gamma_eff in collecting-duct cells (unitless)") # Table 1 note: quinidine gamma_eff,collecting duct cells (2.50)
    km_pgp <- fixed(409)
    label("Quinidine unbound Km for P-gp efflux (ng/mL)") # Table 1 Quinidine K_m,u 409 (P-gp)
    km_cyp3a <- fixed(1996)
    label("Quinidine unbound Km for CYP3A metabolism (ng/mL)") # Table 1 Quinidine K_m,u 1996 (CYP3A) = 409 * 4.88 (Suppl. S1 Case 4)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Physiology scaled to body weight.
    # ------------------------------------------------------------------
    qh <- q_liver * WT
    qhart <- q_hart * WT
    qser <- q_serosa * WT
    qmd <- q_muc_duo * WT
    qmj <- q_muc_jej * WT
    qmi <- q_muc_ile * WT
    qm <- q_muscle * WT
    qs <- q_skin * WT
    qa <- q_adipose * WT
    vb <- v_blood * WT
    vhc5 <- v_hc * WT / 5
    vhe5 <- v_he * WT / 5
    vliv5 <- (v_hc + v_he) * WT / 5
    vm <- v_muscle * WT
    vs <- v_skin * WT
    va <- v_adipose * WT
    vser <- v_serosa * WT
    vmd <- v_muc_duo * WT
    vmj <- v_muc_jej * WT
    vmi <- v_muc_ile * WT
    ved <- v_ent_duo * WT
    vej <- v_ent_jej * WT
    vei <- v_ent_ile * WT
    vld <- v_lumen_duo * WT
    vlj <- v_lumen_jej * WT
    vli <- v_lumen_ile * WT
    qr <- q_renal * WT
    qgfr <- gfr * WT
    qr_gp <- q_glom_pt * WT
    qr_pd <- q_pt_dt * WT
    qr_dc <- q_dt_cd * WT
    qu_pd <- qu_pt_dt * WT
    qu_dc <- qu_dt_cd * WT
    qu_out <- qu_cd * WT
    ec50_rif <- ec50u_rif * 822.94

    # ------------------------------------------------------------------
    # 2. Rifampicin (Supplementary Material S1 'Case 1: rifampicin').
    # ------------------------------------------------------------------
    cl_int_all_rif <- fbclint_all_rif * WT / fb_rif
    vmax_oatp_rif <- 1 / (1 + rdif_rif) * cl_int_all_rif / beta_rif * km_oatp_rif
    ps_dif_inf_rif <- rdif_rif / (1 + rdif_rif) * cl_int_all_rif / beta_rif
    ps_dif_eff_rif <- rdif_rif / (1 + rdif_rif) * cl_int_all_rif / beta_rif / gamma_rif
    cl_met_rif <- rdif_rif / (1 + rdif_rif) * cl_int_all_rif / (1 - beta_rif) / gamma_rif
    clr_rif <- exp(lcl_renal_rif) * WT
    kp_muscle_rif <- sfkp_rif * exp(lkp_muscle_rif)
    kp_skin_rif <- sfkp_rif * exp(lkp_skin_rif)
    kp_adipose_rif <- sfkp_rif * exp(lkp_adipose_rif)
    kp_serosa_rif <- sfkp_rif * exp(lkp_serosa_rif)
    psd_rif <- ps_dif_gut_rif * WT * fsa_duo
    psj_rif <- ps_dif_gut_rif * WT * fsa_jej
    psi_rif <- ps_dif_gut_rif * WT * fsa_ile
    cmd_rif <- cl_met_ent_rif * WT * fsa_duo
    cmj_rif <- cl_met_ent_rif * WT * fsa_jej
    cmi_rif <- cl_met_ent_rif * WT * fsa_ile

    c_b_rif <- central_rif / vb
    c_he1_rif <- is_liver1_rif / vhe5
    c_he2_rif <- is_liver2_rif / vhe5
    c_he3_rif <- is_liver3_rif / vhe5
    c_he4_rif <- is_liver4_rif / vhe5
    c_he5_rif <- is_liver5_rif / vhe5
    c_hc1_rif <- int_liver1_rif / vhc5
    c_hc2_rif <- int_liver2_rif / vhc5
    c_hc3_rif <- int_liver3_rif / vhc5
    c_hc4_rif <- int_liver4_rif / vhc5
    c_hc5_rif <- int_liver5_rif / vhc5
    c_ser_rif <- serosa_rif / vser
    c_m_rif <- muscle_rif / vm
    c_s_rif <- skin_rif / vs
    c_a_rif <- adipose_rif / va
    c_md_rif <- duodenum_muc_rif / vmd
    c_mj_rif <- jejunum_muc_rif / vmj
    c_mi_rif <- ileum_muc_rif / vmi
    c_ed_rif <- duodenum_ent_rif / ved
    c_ej_rif <- jejunum_ent_rif / vej
    c_ei_rif <- ileum_ent_rif / vei
    c_ld_rif <- duodenum_lumen_rif / vld
    c_lj_rif <- jejunum_lumen_rif / vlj
    c_li_rif <- ileum_lumen_rif / vli

    upt1_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * c_he1_rif) + ps_dif_inf_rif) / 5 * c_he1_rif
    upt2_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * c_he2_rif) + ps_dif_inf_rif) / 5 * c_he2_rif
    upt3_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * c_he3_rif) + ps_dif_inf_rif) / 5 * c_he3_rif
    upt4_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * c_he4_rif) + ps_dif_inf_rif) / 5 * c_he4_rif
    upt5_rif <- fb_rif * (vmax_oatp_rif / (km_oatp_rif + fb_rif * c_he5_rif) + ps_dif_inf_rif) / 5 * c_he5_rif
    eff1_rif <- fh_rif * ps_dif_eff_rif / 5 * c_hc1_rif
    eff2_rif <- fh_rif * ps_dif_eff_rif / 5 * c_hc2_rif
    eff3_rif <- fh_rif * ps_dif_eff_rif / 5 * c_hc3_rif
    eff4_rif <- fh_rif * ps_dif_eff_rif / 5 * c_hc4_rif
    eff5_rif <- fh_rif * ps_dif_eff_rif / 5 * c_hc5_rif

    d/dt(central_rif) <- qh * (c_he5_rif - c_b_rif) +
      qm * (c_m_rif / kp_muscle_rif - c_b_rif) +
      qs * (c_s_rif / kp_skin_rif - c_b_rif) +
      qa * (c_a_rif / kp_adipose_rif - c_b_rif) - clr_rif * c_b_rif
    d/dt(is_liver1_rif) <- qhart * c_b_rif + qmd * c_md_rif + qmj * c_mj_rif + qmi * c_mi_rif +
      qser * c_ser_rif / kp_serosa_rif - qh * c_he1_rif + eff1_rif - upt1_rif
    d/dt(is_liver2_rif) <- qh * (c_he1_rif - c_he2_rif) + eff2_rif - upt2_rif
    d/dt(is_liver3_rif) <- qh * (c_he2_rif - c_he3_rif) + eff3_rif - upt3_rif
    d/dt(is_liver4_rif) <- qh * (c_he3_rif - c_he4_rif) + eff4_rif - upt4_rif
    d/dt(is_liver5_rif) <- qh * (c_he4_rif - c_he5_rif) + eff5_rif - upt5_rif
    d/dt(int_liver1_rif) <- upt1_rif - eff1_rif - fh_rif * cl_met_rif * (1 + fm_ugt * (enzyme_ugt_liver1 - 1)) / 5 * c_hc1_rif
    d/dt(int_liver2_rif) <- upt2_rif - eff2_rif - fh_rif * cl_met_rif * (1 + fm_ugt * (enzyme_ugt_liver2 - 1)) / 5 * c_hc2_rif
    d/dt(int_liver3_rif) <- upt3_rif - eff3_rif - fh_rif * cl_met_rif * (1 + fm_ugt * (enzyme_ugt_liver3 - 1)) / 5 * c_hc3_rif
    d/dt(int_liver4_rif) <- upt4_rif - eff4_rif - fh_rif * cl_met_rif * (1 + fm_ugt * (enzyme_ugt_liver4 - 1)) / 5 * c_hc4_rif
    d/dt(int_liver5_rif) <- upt5_rif - eff5_rif - fh_rif * cl_met_rif * (1 + fm_ugt * (enzyme_ugt_liver5 - 1)) / 5 * c_hc5_rif
    d/dt(muscle_rif) <- qm * (c_b_rif - c_m_rif / kp_muscle_rif)
    d/dt(skin_rif) <- qs * (c_b_rif - c_s_rif / kp_skin_rif)
    d/dt(adipose_rif) <- qa * (c_b_rif - c_a_rif / kp_adipose_rif)
    d/dt(serosa_rif) <- qser * (c_b_rif - c_ser_rif / kp_serosa_rif)
    d/dt(duodenum_muc_rif) <- qmd * (c_b_rif - c_md_rif) + fent_rif * psd_rif / gamma_ent_rif * c_ed_rif - fb_rif * psd_rif * c_md_rif
    d/dt(jejunum_muc_rif) <- qmj * (c_b_rif - c_mj_rif) + fent_rif * psj_rif / gamma_ent_rif * c_ej_rif - fb_rif * psj_rif * c_mj_rif
    d/dt(ileum_muc_rif) <- qmi * (c_b_rif - c_mi_rif) + fent_rif * psi_rif / gamma_ent_rif * c_ei_rif - fb_rif * psi_rif * c_mi_rif
    d/dt(duodenum_ent_rif) <- psd_rif * c_ld_rif + fb_rif * psd_rif * c_md_rif -
      fent_rif * (2 * psd_rif / gamma_ent_rif + cmd_rif * (1 + fm_ugt * (enzyme_ugt_duodenum - 1))) * c_ed_rif
    d/dt(jejunum_ent_rif) <- psj_rif * c_lj_rif + fb_rif * psj_rif * c_mj_rif -
      fent_rif * (2 * psj_rif / gamma_ent_rif + cmj_rif * (1 + fm_ugt * (enzyme_ugt_jejunum - 1))) * c_ej_rif
    d/dt(ileum_ent_rif) <- psi_rif * c_li_rif + fb_rif * psi_rif * c_mi_rif -
      fent_rif * (2 * psi_rif / gamma_ent_rif + cmi_rif * (1 + fm_ugt * (enzyme_ugt_ileum - 1))) * c_ei_rif
    d/dt(duodenum_lumen_rif) <- -psd_rif * c_ld_rif - k_feces_duo * duodenum_lumen_rif + fent_rif * psd_rif / gamma_ent_rif * c_ed_rif
    d/dt(jejunum_lumen_rif) <- k_feces_duo * duodenum_lumen_rif - psj_rif * c_lj_rif - k_feces_jej * jejunum_lumen_rif +
      fent_rif * psj_rif / gamma_ent_rif * c_ej_rif
    d/dt(ileum_lumen_rif) <- k_feces_jej * jejunum_lumen_rif - psi_rif * c_li_rif - k_feces_ile * ileum_lumen_rif +
      fent_rif * psi_rif / gamma_ent_rif * c_ei_rif

    # ------------------------------------------------------------------
    # 3. Induction turnover (Supplementary Material S1 Eq. 1): UGT (gut +
    #    liver), CYP3A (gut + liver) and P-gp (gut + kidney).
    # ------------------------------------------------------------------
    uh1 <- fh_rif * c_hc1_rif
    uh2 <- fh_rif * c_hc2_rif
    uh3 <- fh_rif * c_hc3_rif
    uh4 <- fh_rif * c_hc4_rif
    uh5 <- fh_rif * c_hc5_rif
    ued <- fent_rif * c_ed_rif
    uej <- fent_rif * c_ej_rif
    uei <- fent_rif * c_ei_rif
    ub <- fb_rif * c_b_rif
    d/dt(enzyme_ugt_liver1) <- kdeg_liver * (1 + emax_ugt * uh1 / (ec50_rif + uh1) - enzyme_ugt_liver1)
    d/dt(enzyme_ugt_liver2) <- kdeg_liver * (1 + emax_ugt * uh2 / (ec50_rif + uh2) - enzyme_ugt_liver2)
    d/dt(enzyme_ugt_liver3) <- kdeg_liver * (1 + emax_ugt * uh3 / (ec50_rif + uh3) - enzyme_ugt_liver3)
    d/dt(enzyme_ugt_liver4) <- kdeg_liver * (1 + emax_ugt * uh4 / (ec50_rif + uh4) - enzyme_ugt_liver4)
    d/dt(enzyme_ugt_liver5) <- kdeg_liver * (1 + emax_ugt * uh5 / (ec50_rif + uh5) - enzyme_ugt_liver5)
    d/dt(enzyme_ugt_duodenum) <- kdeg_gut * (1 + emax_ugt * ued / (ec50_rif + ued) - enzyme_ugt_duodenum)
    d/dt(enzyme_ugt_jejunum) <- kdeg_gut * (1 + emax_ugt * uej / (ec50_rif + uej) - enzyme_ugt_jejunum)
    d/dt(enzyme_ugt_ileum) <- kdeg_gut * (1 + emax_ugt * uei / (ec50_rif + uei) - enzyme_ugt_ileum)
    d/dt(enzyme_3a4_liver1) <- kdeg_liver * (1 + emax_cyp3a * uh1 / (ec50_rif + uh1) - enzyme_3a4_liver1)
    d/dt(enzyme_3a4_liver2) <- kdeg_liver * (1 + emax_cyp3a * uh2 / (ec50_rif + uh2) - enzyme_3a4_liver2)
    d/dt(enzyme_3a4_liver3) <- kdeg_liver * (1 + emax_cyp3a * uh3 / (ec50_rif + uh3) - enzyme_3a4_liver3)
    d/dt(enzyme_3a4_liver4) <- kdeg_liver * (1 + emax_cyp3a * uh4 / (ec50_rif + uh4) - enzyme_3a4_liver4)
    d/dt(enzyme_3a4_liver5) <- kdeg_liver * (1 + emax_cyp3a * uh5 / (ec50_rif + uh5) - enzyme_3a4_liver5)
    d/dt(enzyme_3a4_duodenum) <- kdeg_gut * (1 + emax_cyp3a * ued / (ec50_rif + ued) - enzyme_3a4_duodenum)
    d/dt(enzyme_3a4_jejunum) <- kdeg_gut * (1 + emax_cyp3a * uej / (ec50_rif + uej) - enzyme_3a4_jejunum)
    d/dt(enzyme_3a4_ileum) <- kdeg_gut * (1 + emax_cyp3a * uei / (ec50_rif + uei) - enzyme_3a4_ileum)
    d/dt(pgp_duodenum) <- kdeg_gut * (1 + emax_pgp * ued / (ec50_rif + ued) - pgp_duodenum)
    d/dt(pgp_jejunum) <- kdeg_gut * (1 + emax_pgp * uej / (ec50_rif + uej) - pgp_jejunum)
    d/dt(pgp_ileum) <- kdeg_gut * (1 + emax_pgp * uei / (ec50_rif + uei) - pgp_ileum)
    d/dt(pgp_kidney) <- kdeg_pgp_kidney * (1 + emax_pgp * ub / (ec50_rif + ub) - pgp_kidney)
    enzyme_ugt_liver1(0) <- 1
    enzyme_ugt_liver2(0) <- 1
    enzyme_ugt_liver3(0) <- 1
    enzyme_ugt_liver4(0) <- 1
    enzyme_ugt_liver5(0) <- 1
    enzyme_ugt_duodenum(0) <- 1
    enzyme_ugt_jejunum(0) <- 1
    enzyme_ugt_ileum(0) <- 1
    enzyme_3a4_liver1(0) <- 1
    enzyme_3a4_liver2(0) <- 1
    enzyme_3a4_liver3(0) <- 1
    enzyme_3a4_liver4(0) <- 1
    enzyme_3a4_liver5(0) <- 1
    enzyme_3a4_duodenum(0) <- 1
    enzyme_3a4_jejunum(0) <- 1
    enzyme_3a4_ileum(0) <- 1
    pgp_duodenum(0) <- 1
    pgp_jejunum(0) <- 1
    pgp_ileum(0) <- 1
    pgp_kidney(0) <- 1

    inh_d <- 1 + ued / ki_pgp_rif
    inh_j <- 1 + uej / ki_pgp_rif
    inh_i <- 1 + uei / ki_pgp_rif
    inh_k <- 1 + ub / ki_pgp_rif
    pgp_activity_duodenum <- pgp_duodenum / inh_d
    pgp_activity_jejunum <- pgp_jejunum / inh_j
    pgp_activity_ileum <- pgp_ileum / inh_i
    pgp_activity_kidney <- pgp_kidney / inh_k

    # ------------------------------------------------------------------
    # 4. Quinidine (Supplementary Material S1 'Case 4: quinidine').
    #    Perfusion-limited liver with saturable, induced CYP3A metabolism;
    #    saturable intestinal CYP3A and P-gp; no enterohepatic circulation.
    # ------------------------------------------------------------------
    sfkp_study <- sfkp +
      STUDY_TWUM1981_IV * (sfkp_twiv - sfkp) +
      STUDY_TWUM1981_PO * (sfkp_twpo - sfkp) +
      STUDY_DAMKIER1999 * (sfkp_dam - sfkp)
    fbclint_study <- fbclint_all +
      STUDY_TWUM1981_IV * (fbclint_twiv - fbclint_all) +
      STUDY_TWUM1981_PO * (fbclint_twpo - fbclint_all) +
      STUDY_DAMKIER1999 * (fbclint_dam - fbclint_all)
    kstomach <- exp(lkstomach)
    cl_int_liver <- fbclint_study * WT / fb
    cl_met_cyp3a <- cl_int_liver * fm_cyp3a4
    cl_met_noncyp3a <- cl_int_liver * (1 - fm_cyp3a4)
    kpl <- sfkp_study * exp(lkp_liver)
    kp_muscle <- sfkp_study * exp(lkp_muscle)
    kp_skin <- sfkp_study * exp(lkp_skin)
    kp_adipose <- sfkp_study * exp(lkp_adipose)
    kp_serosa <- sfkp_study * exp(lkp_serosa)

    psd <- ps_dif_gut * WT * fsa_duo
    psj <- ps_dif_gut * WT * fsa_jej
    psi <- ps_dif_gut * WT * fsa_ile
    pgpd <- psd / gamma_ent / rdif_eff_duo * fpgp_ent
    pgpj <- psj / gamma_ent / rdif_eff_jej * fpgp_ent
    pgpi <- psi / gamma_ent / rdif_eff_ile * fpgp_ent
    npgpd <- psd / gamma_ent / rdif_eff_duo * (1 - fpgp_ent)
    npgpj <- psj / gamma_ent / rdif_eff_jej * (1 - fpgp_ent)
    npgpi <- psi / gamma_ent / rdif_eff_ile * (1 - fpgp_ent)
    # Intestinal CYP3A metabolic clearance per segment (CYP3A share fractions
    # 0.138 / 0.545 / 0.317). The enterocyte f_m,CYP3A is 1 (Table 1), so all
    # intestinal metabolism is the induced CYP3A term with no non-CYP3A part.
    cm3d <- cl_met_ent * fcyp_duo
    cm3j <- cl_met_ent * fcyp_jej
    cm3i <- cl_met_ent * fcyp_ile

    ps_pcv <- ps_dif_pc_v * sf_renal * WT
    ps_vpc <- ps_pcv * gamma_inf_pc
    ps_act_vpc <- ps_pcv * gamma_inf_pc / rdif_inf_kidney
    ps_pct <- ps_dif_pc_t * sf_renal * WT
    ps_tpc <- ps_pct * gamma_eff_pc
    pgp_pc <- ps_pct / rdif_eff_kidney * fpgp_kidney
    npgp_pc <- ps_pct / rdif_eff_kidney * (1 - fpgp_kidney)
    ps_dcv <- ps_dif_dc * sf_renal * WT
    ps_vdc <- ps_dcv * gamma_inf_dc
    ps_dct <- ps_dif_dc * sf_renal * WT
    ps_tdc <- ps_dct * gamma_eff_dc
    ps_ccv <- ps_dif_cc * sf_renal * WT
    ps_vcc <- ps_ccv * gamma_inf_cc
    ps_cct <- ps_dif_cc * sf_renal * WT
    ps_tcc <- ps_cct * gamma_eff_cc
    fp <- fb * rb

    c_b <- central / vb
    c_l1 <- liver1 / vliv5
    c_l2 <- liver2 / vliv5
    c_l3 <- liver3 / vliv5
    c_l4 <- liver4 / vliv5
    c_l5 <- liver5 / vliv5
    c_ser <- serosa / vser
    c_m <- muscle / vm
    c_s <- skin / vs
    c_a <- adipose / va
    c_md <- duodenum_muc / vmd
    c_mj <- jejunum_muc / vmj
    c_mi <- ileum_muc / vmi
    c_ed <- duodenum_ent / ved
    c_ej <- jejunum_ent / vej
    c_ei <- ileum_ent / vei
    c_ld <- duodenum_lumen / vld
    c_lj <- jejunum_lumen / vlj
    c_li <- ileum_lumen / vli
    c_glom_b <- glom_b / (v_glom_b * WT)
    c_glom_u <- glom_u / (v_glom_u * WT)
    c_pt_b <- pt_b / (v_pt_b * WT)
    c_pt_c <- pt_c / (v_pt_c * WT)
    c_pt_u <- pt_u / (v_pt_u * WT)
    c_dt_b <- dt_b / (v_dt_b * WT)
    c_dt_c <- dt_c / (v_dt_c * WT)
    c_dt_u <- dt_u / (v_dt_u * WT)
    c_cd_b <- cd_b / (v_cd_b * WT)
    c_cd_c <- cd_c / (v_cd_c * WT)
    c_cd_u <- cd_u / (v_cd_u * WT)

    # Hepatic saturable CYP3A metabolism (induced by rifampicin) plus
    # non-CYP3A metabolism, per perfusion-limited unit. The unbound hepatic
    # concentration is fB * C_liver / (SF_Kp * Kp_liver).
    full1 <- fb * c_l1 / kpl
    full2 <- fb * c_l2 / kpl
    full3 <- fb * c_l3 / kpl
    full4 <- fb * c_l4 / kpl
    full5 <- fb * c_l5 / kpl
    met1 <- (cl_met_cyp3a / 5 * km_cyp3a * enzyme_3a4_liver1 / (km_cyp3a + full1) + cl_met_noncyp3a / 5) * full1
    met2 <- (cl_met_cyp3a / 5 * km_cyp3a * enzyme_3a4_liver2 / (km_cyp3a + full2) + cl_met_noncyp3a / 5) * full2
    met3 <- (cl_met_cyp3a / 5 * km_cyp3a * enzyme_3a4_liver3 / (km_cyp3a + full3) + cl_met_noncyp3a / 5) * full3
    met4 <- (cl_met_cyp3a / 5 * km_cyp3a * enzyme_3a4_liver4 / (km_cyp3a + full4) + cl_met_noncyp3a / 5) * full4
    met5 <- (cl_met_cyp3a / 5 * km_cyp3a * enzyme_3a4_liver5 / (km_cyp3a + full5) + cl_met_noncyp3a / 5) * full5

    # Enterocyte apical efflux: passive + saturable induced/inhibited P-gp;
    # enterocyte metabolism: saturable induced CYP3A + non-CYP3A.
    effd <- fent * (psd / gamma_ent + pgpd * km_pgp * pgp_duodenum / (km_pgp * inh_d + fent * c_ed) + npgpd) * c_ed
    effj <- fent * (psj / gamma_ent + pgpj * km_pgp * pgp_jejunum / (km_pgp * inh_j + fent * c_ej) + npgpj) * c_ej
    effi <- fent * (psi / gamma_ent + pgpi * km_pgp * pgp_ileum / (km_pgp * inh_i + fent * c_ei) + npgpi) * c_ei
    metd <- fent * cm3d * km_cyp3a * enzyme_3a4_duodenum / (km_cyp3a + fent * c_ed) * c_ed
    metj <- fent * cm3j * km_cyp3a * enzyme_3a4_jejunum / (km_cyp3a + fent * c_ej) * c_ej
    meti <- fent * cm3i * km_cyp3a * enzyme_3a4_ileum / (km_cyp3a + fent * c_ei) * c_ei

    eff_pc <- frenal * (pgp_pc * km_pgp * pgp_kidney / (km_pgp * inh_k + frenal * c_pt_c) + npgp_pc + ps_pct) * c_pt_c

    d/dt(central) <- qh * (c_l5 / kpl - c_b) +
      qm * (c_m / kp_muscle - c_b) +
      qs * (c_s / kp_skin - c_b) +
      qa * (c_a / kp_adipose - c_b) +
      (qr - qu_out) * c_cd_b - qr * c_b
    d/dt(liver1) <- qhart * c_b + qmd * c_md + qmj * c_mj + qmi * c_mi +
      qser * c_ser / kp_serosa - qh * c_l1 / kpl - met1
    d/dt(liver2) <- qh * (c_l1 - c_l2) / kpl - met2
    d/dt(liver3) <- qh * (c_l2 - c_l3) / kpl - met3
    d/dt(liver4) <- qh * (c_l3 - c_l4) / kpl - met4
    d/dt(liver5) <- qh * (c_l4 - c_l5) / kpl - met5
    d/dt(muscle) <- qm * (c_b - c_m / kp_muscle)
    d/dt(skin) <- qs * (c_b - c_s / kp_skin)
    d/dt(adipose) <- qa * (c_b - c_a / kp_adipose)
    d/dt(serosa) <- qser * (c_b - c_ser / kp_serosa)
    d/dt(duodenum_muc) <- qmd * (c_b - c_md) + fent * psd / gamma_ent * c_ed - fb * psd * c_md
    d/dt(jejunum_muc) <- qmj * (c_b - c_mj) + fent * psj / gamma_ent * c_ej - fb * psj * c_mj
    d/dt(ileum_muc) <- qmi * (c_b - c_mi) + fent * psi / gamma_ent * c_ei - fb * psi * c_mi
    d/dt(duodenum_ent) <- psd * c_ld + fb * psd * c_md - fent * psd / gamma_ent * c_ed - effd - metd
    d/dt(jejunum_ent) <- psj * c_lj + fb * psj * c_mj - fent * psj / gamma_ent * c_ej - effj - metj
    d/dt(ileum_ent) <- psi * c_li + fb * psi * c_mi - fent * psi / gamma_ent * c_ei - effi - meti
    d/dt(stomach) <- -kstomach * stomach
    d/dt(duodenum_lumen) <- kstomach * stomach - psd * c_ld - k_feces_duo * duodenum_lumen + effd
    d/dt(jejunum_lumen) <- k_feces_duo * duodenum_lumen - psj * c_lj - k_feces_jej * jejunum_lumen + effj
    d/dt(ileum_lumen) <- k_feces_jej * jejunum_lumen - psi * c_li - k_feces_ile * ileum_lumen + effi
    d/dt(a_feces) <- k_feces_ile * ileum_lumen

    d/dt(glom_b) <- qr * c_b - qr_gp * c_glom_b - fp * qgfr * c_glom_b / rb
    d/dt(glom_u) <- fp * qgfr * c_glom_b / rb - qgfr * c_glom_u
    d/dt(pt_b) <- qr_gp * c_glom_b - qr_pd * c_pt_b + frenal * ps_pcv * c_pt_c - fb * (ps_act_vpc + ps_vpc) * c_pt_b
    d/dt(pt_c) <- -frenal * ps_pcv * c_pt_c + fb * (ps_act_vpc + ps_vpc) * c_pt_b - eff_pc + ps_tpc * c_pt_u
    d/dt(pt_u) <- qgfr * c_glom_u - qu_pd * c_pt_u + eff_pc - ps_tpc * c_pt_u
    d/dt(dt_b) <- qr_pd * c_pt_b - qr_dc * c_dt_b + frenal * ps_dcv * c_dt_c - fb * ps_vdc * c_dt_b
    d/dt(dt_c) <- -frenal * ps_dcv * c_dt_c + fb * ps_vdc * c_dt_b - frenal * ps_dct * c_dt_c + ps_tdc * c_dt_u
    d/dt(dt_u) <- qu_pd * c_pt_u - qu_dc * c_dt_u + frenal * ps_dct * c_dt_c - ps_tdc * c_dt_u
    d/dt(cd_b) <- qr_dc * c_dt_b - (qr - qu_out) * c_cd_b + frenal * ps_ccv * c_cd_c - fb * ps_vcc * c_cd_b
    d/dt(cd_c) <- -frenal * ps_ccv * c_cd_c + fb * ps_vcc * c_cd_b - frenal * ps_cct * c_cd_c + ps_tcc * c_cd_u
    d/dt(cd_u) <- qu_dc * c_dt_u - qu_out * c_cd_u + frenal * ps_cct * c_cd_c - ps_tcc * c_cd_u
    d/dt(a_urine) <- qu_out * c_cd_u

    # Rifampicin oral doses go to duodenum_lumen_rif after the Table 1 lag;
    # oral quinidine to `stomach` (Table 1 quinidine T_lag 0, no lag);
    # intravenous quinidine to `central` (give the 30-min infusion with `dur`).
    alag(duodenum_lumen_rif) <- exp(ltlag_rif)

    # ------------------------------------------------------------------
    # 5. Observations: blood concentrations (ng/mL). Deterministic.
    # ------------------------------------------------------------------
    Cc <- c_b
    Cc_rif <- c_b_rif
  })
}
