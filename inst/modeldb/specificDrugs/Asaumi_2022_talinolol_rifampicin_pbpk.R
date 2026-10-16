Asaumi_2022_talinolol_rifampicin_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body; five-unit tandem liver, three-segment segregated-flow",
    "intestine, nephron-resolved kidney; Napp). Rifampicin-talinolol",
    "P-glycoprotein (P-gp) drug-drug interaction in healthy adults. The",
    "rifampicin perpetrator model (saturable OATP1B hepatic uptake, UGT",
    "autoinduction, duodenum / jejunum / ileum gut) induces P-gp by a",
    "turnover model (Emax 4.01, EC50,u 0.0639 uM) and inhibits it",
    "competitively (Ki,u 402 ng/mL) in the enterocytes of each intestinal",
    "segment, in each hepatocyte unit and in the renal proximal tubule. The",
    "talinolol victim model has a stomach, intestinal active uptake,",
    "saturable P-gp efflux (Km,u 4071 ng/mL) into the gut lumen, the bile",
    "and the proximal-tubule lumen, a permeability-limited liver with",
    "enterohepatic circulation, and perfusion-limited muscle, skin, adipose",
    "and serosa. The study-specific talinolol parameter sets (Westphal 2000",
    "default, Tubic 2006 via STUDY_TUBIC2006) are switched inside the model.",
    "Deterministic: no IIV or residual error published."
  )
  reference <- paste(
    "Asaumi R, Nunoya K, Yamaura Y, Taskar KS, Sugiyama Y. Robust",
    "physiologically based pharmacokinetic model of rifampicin for predicting",
    "drug-drug interactions via P-glycoprotein induction and inhibition in",
    "the intestine, liver, and kidney. CPT Pharmacometrics Syst Pharmacol.",
    "2022;11:919-933. doi:10.1002/psp4.12807 (PMC9286720).",
    "ODEs: Supplementary Material S1 ('Case 1: rifampicin', 'Case 3:",
    "talinolol', which lists only its differences from 'Case 2: digoxin').",
    "Physiology: Table S1; fitted parameters: Table S2."
  )
  vignette <- "Asaumi_2022_rifampicin_pgp_ddi_pbpk"

  # State names follow the Sugiyama-laboratory PBPK models already in the
  # package (Asaumi_2019_*_rifampicin_pbpk.R, Tsuchitani_2024_telmisartan_pbpk.R,
  # Nishiyama_2019_metformin_pbpk.R): `is_liver<n>` = hepatic extracellular
  # unit n and `int_liver<n>` = hepatocyte unit n of the five-unit tandem
  # liver; `<segment>_lumen` / `<segment>_ent` / `<segment>_muc` = gut lumen,
  # enterocytes and mucosal blood of the duodenum, jejunum and ileum;
  # `ehc<n>` = biliary transit chain; nephron spaces are `glom_*`, `pt_*`,
  # `dt_*`, `cd_*` (glomerulus, proximal tubule, distal tubule, collecting
  # duct) with `_b` (blood vessel), `_c` (tubular cells) and `_u` (urinary
  # lumen). Talinolol carries the bare names and rifampicin the `_rif`
  # suffix. `enzyme_ugt_*` and `pgp_*` hold the relative UGT / P-gp amount
  # (induced / control) of Supplementary Material S1 Eq. 1.
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
    "ehc1",
    "ehc2",
    "ehc3",
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
    "pgp_liver1",
    "pgp_liver2",
    "pgp_liver3",
    "pgp_liver4",
    "pgp_liver5",
    "pgp_duodenum",
    "pgp_jejunum",
    "pgp_ileum",
    "pgp_kidney"
  )

  # Every state is an AMOUNT in ug (talinolol and rifampicin) or a unitless
  # relative enzyme / transporter amount. Volumes are in L, so amount /
  # volume is ug/L = ng/mL, the unit of the printed Km,u and Ki values.
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
    STUDY_TUBIC2006 = list(
      description = paste(
        "1 = simulate the Tubic 2006 single-dose talinolol study (25 to 400 mg",
        "PO; nonlinear absorption); 0 = Westphal 2000 (the rifampicin DDI study)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Table S2B fits k_stomach, SF_Kp and fB*CLint,all separately to each",
        "clinical study; this indicator selects the Tubic 2006 set."
      ),
      source_name = "Table S2B column 'Tubic et al. (2006)'"
    )
  )

  compartmentData <- list(
    central = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    stomach = list(analyte = "talinolol", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum_lumen = list(analyte = "talinolol", units = "ug", specimen = "administration site", verified = TRUE),
    jejunum_lumen = list(analyte = "talinolol", units = "ug", specimen = "administration site", verified = TRUE),
    ileum_lumen = list(analyte = "talinolol", units = "ug", specimen = "administration site", verified = TRUE),
    duodenum_ent = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    jejunum_ent = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    ileum_ent = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    duodenum_muc = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    jejunum_muc = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    ileum_muc = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    a_feces = list(analyte = "talinolol", units = "ug", specimen = "faeces", verified = TRUE),
    serosa = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    muscle = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    skin = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    adipose = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver1 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver2 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver3 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver4 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    is_liver5 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver1 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver2 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver3 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver4 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    int_liver5 = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    ehc1 = list(analyte = "talinolol", units = "ug", specimen = "bile", verified = TRUE),
    ehc2 = list(analyte = "talinolol", units = "ug", specimen = "bile", verified = TRUE),
    ehc3 = list(analyte = "talinolol", units = "ug", specimen = "bile", verified = TRUE),
    glom_b = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    glom_u = list(analyte = "talinolol", units = "ug", specimen = "urine", verified = TRUE),
    pt_b = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    pt_c = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    pt_u = list(analyte = "talinolol", units = "ug", specimen = "urine", verified = TRUE),
    dt_b = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    dt_c = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    dt_u = list(analyte = "talinolol", units = "ug", specimen = "urine", verified = TRUE),
    cd_b = list(analyte = "talinolol", units = "ug", specimen = "whole blood", verified = TRUE),
    cd_c = list(analyte = "talinolol", units = "ug", specimen = "tissue", verified = TRUE),
    cd_u = list(analyte = "talinolol", units = "ug", specimen = "urine", verified = TRUE),
    a_urine = list(analyte = "talinolol", units = "ug", specimen = "urine", verified = TRUE),
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
    pgp_liver1 = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_liver2 = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_liver3 = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_liver4 = list(
      analyte = "P-gp relative amount (induced / control)",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    pgp_liver5 = list(
      analyte = "P-gp relative amount (induced / control)",
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
      "rifampicin 600 mg PO once daily for 7 or 9 days; talinolol 30 mg IV",
      "(30-min infusion) or 100 mg PO 13 h after the last rifampicin dose",
      "(Westphal 2000); talinolol 25 to 400 mg PO alone (Tubic 2006)"
    ),
    notes = paste(
      "Mean blood profiles digitised from the literature, not individual data.",
      "Westphal 2000 (ref. 10): talinolol IV and PO before and after",
      "rifampicin; Tubic 2006 (ref. 23): single oral talinolol doses of 25 to",
      "400 mg used for the nonlinear (P-gp saturation) fit."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Physiology, per kg body weight (Table S1A, S1B, S1C, S1D). Shared by
    # rifampicin and talinolol.
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
    k_bile <- fixed(0.605)
    label("Transit rate constant in the enterohepatic circulation (1/h)") # Table S1B k_bile 0.605
    k_feces_duo <- fixed(5.13)
    label("Duodenum lumen transit rate constant (1/h)") # Table S1B k_feces duodenum part 5.13
    k_feces_jej <- fixed(1.47)
    label("Jejunum lumen transit rate constant (1/h)") # Table S1B k_feces jejunum part 1.47
    k_feces_ile <- fixed(1.21)
    label("Ileum lumen transit rate constant to faeces (1/h)") # Table S1B k_feces ileum part 1.21
    fsa_duo <- fixed(0.076)
    label("Duodenum fraction of intestinal surface area (fraction)") # Table S1C surface area 7.6%; Suppl. S1 'Other equations' x 0.076
    fsa_jej <- fixed(0.394)
    label("Jejunum fraction of intestinal surface area (fraction)") # Table S1C surface area 39.4%; Suppl. S1 x 0.394
    fsa_ile <- fixed(0.530)
    label("Ileum fraction of intestinal surface area (fraction)") # Table S1C surface area 53.0%; Suppl. S1 x 0.530
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
    label("Glomerulus blood-vessel volume (L/kg)") # Table S1D glomerulus 0.0000465 (footnotes a,b: blood vessels)
    v_glom_u <- fixed(0.0000713)
    label("Glomerulus urinary-lumen volume (L/kg)") # Table S1D glomerulus 0.0000713 (footnotes a,e: urinary lumen)
    v_pt_b <- fixed(0.000474)
    label("Proximal tubule blood-vessel volume (L/kg)") # Table S1D proximal tubule 0.000474 (a,b)
    v_pt_c <- fixed(0.000395)
    label("Proximal tubule cell volume (L/kg)") # Table S1D proximal tubule 0.000395 (a,f: cells)
    v_pt_u <- fixed(0.000727)
    label("Proximal tubule urinary-lumen volume (L/kg)") # Table S1D proximal tubule 0.000727 (a,e)
    v_dt_b <- fixed(0.000102)
    label("Distal tubule blood-vessel volume (L/kg)") # Table S1D distal tubule 0.000102 (a,b)
    v_dt_c <- fixed(0.0000852)
    label("Distal tubule cell volume (L/kg)") # Table S1D distal tubule 0.0000852 (a,f)
    v_dt_u <- fixed(0.000157)
    label("Distal tubule urinary-lumen volume (L/kg)") # Table S1D distal tubule 0.000157 (a,e)
    v_cd_b <- fixed(0.000577)
    label("Collecting duct blood-vessel volume (L/kg)") # Table S1D collecting duct 0.000577 (a,b)
    v_cd_c <- fixed(0.000480)
    label("Collecting duct cell volume (L/kg)") # Table S1D collecting duct 0.000480 (a,f)
    v_cd_u <- fixed(0.000884)
    label("Collecting duct urinary-lumen volume (L/kg)") # Table S1D collecting duct 0.000884 (a,e)
    kdeg_liver <- fixed(0.0158)
    label("Degradation rate constant of hepatic UGT and P-gp (1/h)") # Table S1E hepatic k_deg 0.0158 (UGT, P-gp assumed = CYP3A4)
    kdeg_gut <- fixed(0.0288)
    label("Degradation rate constant of intestinal UGT and P-gp (1/h)") # Table S1E intestinal k_deg 0.0288 (UGT, P-gp assumed = CYP3A4)
    kdeg_pgp_kidney <- fixed(0.0158)
    label("Degradation rate constant of renal P-gp (1/h)") # Table S1E renal P-gp k_deg 0.0158 (assumed = hepatic CYP3A4)

    # ------------------------------------------------------------------
    # Rifampicin (perpetrator). Table 1 'Rifampicin' column and its note
    # (carried from Asaumi 2019, ref. 15), with the P-gp Emax and Ki
    # estimated in this paper (Table S2C).
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
    label("Rifampicin Rdif,inf,liver = PSdif,inf / PSact,inf (unitless)") # Table 1 Rifampicin R_dif,inf,liver 0.129
    beta_rif <- fixed(0.2)
    label("Rifampicin hepatic beta (unitless)") # Table 1 Rifampicin beta_liver 0.2
    gamma_rif <- fixed(0.778)
    label("Rifampicin gamma_inf,liver = PSdif,inf / PSdif,eff (unitless)") # Table 1 Rifampicin gamma_inf,liver 0.778
    fm_ugt <- fixed(0.759)
    label("Rifampicin fraction metabolised by UGT, liver and enterocytes (fraction)") # Table 1 Rifampicin f_m,liver and f_m,enterocytes 0.759 (UGT)
    km_oatp_rif <- fixed(146)
    label("Rifampicin unbound Km for OATP1B hepatic uptake (ng/mL)") # Table 1 Rifampicin K_m,u 146 (OATP1B)
    ps_dif_gut_rif <- fixed(0.08)
    label("Rifampicin passive gut-lumen-to-enterocyte clearance, whole intestine (L/h/kg)") # Table 1 Rifampicin PS_dif,gut lumen to enterocytes 0.08 (footnote g)
    gamma_ent_rif <- fixed(0.778)
    label("Rifampicin enterocyte gamma_inf = gamma_eff (unitless)") # Table 1 Rifampicin gamma_inf,enterocytes 0.778 (footnote j: = gamma_eff,enterocytes)
    cl_met_ent_rif <- fixed(0.005)
    label("Rifampicin enterocyte metabolic intrinsic clearance, whole intestine (L/h/kg)") # Table 1 Rifampicin CL_int,met,enterocytes 0.005 (footnote g)
    lcl_renal_rif <- fixed(log(0.011))
    label("Rifampicin renal clearance (L/h/kg)") # Table 1 note: rifampicin CL_renal (0.011 L/h/kg)
    ec50u_rif <- fixed(0.0639)
    label("Rifampicin unbound EC50 for induction, all enzymes and P-gp (umol/L)") # Table 1 note EC50,u 63.9 nmol/L; Table S2C 0.0639 umol/L
    emax_ugt <- fixed(1.34)
    label("Rifampicin maximum induction effect on UGT (unitless)") # Table 1 note: E_max,UGT (1.34, autoinduction)
    emax_pgp <- 4.01
    label("Rifampicin maximum induction effect on P-gp (unitless)") # Table S2C E_max for P-gp 4.01 +/- 1.31
    ki_pgp_rif <- 402
    label("Rifampicin unbound inhibition constant for P-gp (ng/mL)") # Table S2C K_i for P-gp 402 +/- 352 ng/mL (0.488 umol/L)

    # ------------------------------------------------------------------
    # Talinolol (victim). Table 1 'Talinolol' column and note; the three
    # study-fitted parameters are Table S2B (Westphal 2000 is the default,
    # Tubic 2006 is selected by STUDY_TUBIC2006).
    # ------------------------------------------------------------------
    lkstomach <- log(0.27)
    label("Talinolol stomach-to-duodenum transit rate constant, Westphal 2000 (1/h)") # Table S2B Westphal (2000) k_stomach 0.27 +/- 0.06
    lkstomach_tubic2006 <- log(0.59)
    label("Talinolol stomach-to-duodenum transit rate constant, Tubic 2006 (1/h)") # Table S2B Tubic (2006) k_stomach 0.59 +/- 0.11
    sfkp <- 2.3
    label("Talinolol common scaling factor on in silico Kp values, Westphal 2000 (unitless)") # Table S2B Westphal (2000) SF_Kp 2.3 +/- 0.5
    sfkp_tubic2006 <- 1.7
    label("Talinolol common scaling factor on in silico Kp values, Tubic 2006 (unitless)") # Table S2B Tubic (2006) SF_Kp 1.7 +/- 0.9
    fbclint_all <- 0.64
    label("Talinolol unbound-blood overall hepatic intrinsic clearance, Westphal 2000 (L/h/kg)") # Table S2B Westphal (2000) f_B CL_int,all 0.64 +/- 0.05
    fbclint_all_tubic2006 <- 0.81
    label("Talinolol unbound-blood overall hepatic intrinsic clearance, Tubic 2006 (L/h/kg)") # Table S2B Tubic (2006) f_B CL_int,all 0.81 +/- 0.07
    fb <- fixed(0.478)
    label("Talinolol unbound fraction in blood (fraction)") # Table 1 Talinolol f_B 0.478
    rb <- fixed(0.941)
    label("Talinolol blood-to-plasma concentration ratio (ratio)") # Table 1 Talinolol R_B 0.941
    fh <- fixed(0.893)
    label("Talinolol unbound fraction in hepatocytes (fraction)") # Table 1 Talinolol f_hepatocytes 0.893
    fent <- fixed(0.724)
    label("Talinolol unbound fraction in enterocytes (fraction)") # Table 1 Talinolol f_enterocytes 0.724
    frenal <- fixed(0.881)
    label("Talinolol unbound fraction in renal cells (fraction)") # Table 1 Talinolol f_renal cell 0.881
    lkp_muscle <- fixed(log(0.861))
    label("Talinolol muscle-to-blood partition coefficient (unitless)") # Table 1 Talinolol K_p,muscle 0.861
    lkp_skin <- fixed(log(0.877))
    label("Talinolol skin-to-blood partition coefficient (unitless)") # Table 1 Talinolol K_p,skin 0.877
    lkp_adipose <- fixed(log(0.383))
    label("Talinolol adipose-to-blood partition coefficient (unitless)") # Table 1 Talinolol K_p,adipose 0.383
    lkp_serosa <- fixed(log(0.908))
    label("Talinolol serosa-to-blood partition coefficient (unitless)") # Table 1 Talinolol K_p,serosa 0.908
    rdif_liver <- fixed(1)
    label("Talinolol Rdif,inf,liver = PSdif,inf / PSact,inf (unitless)") # Table 1 Talinolol R_dif,inf,liver 1 (footnote e: assumption)
    beta_liver <- fixed(0.5)
    label("Talinolol hepatic beta (unitless)") # Table 1 Talinolol beta_liver 0.2/0.5/0.8; Table S2B fitted with 0.5
    gamma_liver <- fixed(4.15)
    label("Talinolol gamma_inf,liver = PSdif,inf / PSdif,eff (unitless)") # Table 1 Talinolol gamma_inf,liver 4.15
    f_bile <- fixed(0.914)
    label("Talinolol biliary fraction CLint,bile / (CLint,bile + CLint,met) (fraction)") # Table 1 Talinolol f_bile 0.914
    fpgp_liver <- fixed(1)
    label("Talinolol P-gp fraction of biliary efflux (fraction)") # Table 1 Talinolol f_P-gp,liver 1 (footnote e: assumption)
    ps_dif_gut <- fixed(0.016)
    label("Talinolol passive gut-lumen-to-enterocyte clearance, whole intestine (L/h/kg)") # Table 1 Talinolol PS_dif,gut lumen to enterocytes 0.016 (footnote h)
    rdif_inf_gut <- fixed(0.5)
    label("Talinolol Rdif,inf,gut lumen to enterocytes = PSdif / PSact (unitless)") # Table 1 note: talinolol R_dif,inf,gut lumen to enterocytes (0.5)
    rdif_eff_duo <- fixed(1.05)
    label("Talinolol Rdif,eff in the duodenum = PSdif,eff / PS_P-gp (unitless)") # Table 1 Talinolol R_dif,eff,duodenum 1.05
    rdif_eff_jej <- fixed(0.292)
    label("Talinolol Rdif,eff in the jejunum = PSdif,eff / PS_P-gp (unitless)") # Table 1 Talinolol R_dif,eff,jejunum 0.292
    rdif_eff_ile <- fixed(0.239)
    label("Talinolol Rdif,eff in the ileum = PSdif,eff / PS_P-gp (unitless)") # Table 1 Talinolol R_dif,eff,ileum 0.239
    gamma_ent <- fixed(1)
    label("Talinolol enterocyte gamma_inf = gamma_eff (unitless)") # Table 1 Talinolol gamma_inf,enterocytes 1 (footnotes e, j)
    fpgp_ent <- fixed(1)
    label("Talinolol P-gp fraction of enterocyte apical efflux (fraction)") # Table 1 Talinolol f_P-gp,enterocytes 1 (footnote e)
    rdif_inf_kidney <- fixed(0.012)
    label("Talinolol Rdif,inf,kidney = PSdif,inf / PSact,inf, basolateral (unitless)") # Table 1 Talinolol R_dif,inf,kidney 0.012 (footnote k)
    rdif_eff_kidney <- fixed(1)
    label("Talinolol Rdif,eff,kidney = PSdif,eff / PS_P-gp, apical (unitless)") # Table 1 Talinolol R_dif,eff,kidney 1 (footnote k)
    fpgp_kidney <- fixed(1)
    label("Talinolol P-gp fraction of proximal-tubule apical efflux (fraction)") # Table 1 Talinolol f_P-gp,kidney 1 (footnote e)
    ps_dif_pc_v <- fixed(0.000617)
    label("Talinolol passive proximal cell-to-vessel clearance (L/h/kg)") # Table 1 Talinolol PS_dif,proximal cells to vessels 0.000617
    ps_dif_pc_t <- fixed(0.00464)
    label("Talinolol passive proximal cell-to-tubule clearance (L/h/kg)") # Table 1 Talinolol PS_dif,proximal cells to tubule 0.00464
    ps_dif_dc <- fixed(0.00016)
    label("Talinolol passive distal cell-to-vessel = cell-to-tubule clearance (L/h/kg)") # Table 1 Talinolol PS_dif,distal cells to vessels 0.00016 (footnote m)
    ps_dif_cc <- fixed(0.0000343)
    label("Talinolol passive collecting-duct cell-to-vessel = cell-to-tubule clearance (L/h/kg)") # Table 1 Talinolol PS_dif,collecting duct cells to vessels 0.0000343 (footnote n)
    sf_renal <- fixed(1)
    label("Talinolol scaling factor on renal passive permeability (unitless)") # Suppl. S1 talinolol 'Other equations' multiply by SF_renal permeability; Table 1 note prints a value only for quinidine
    gamma_inf_pc <- fixed(8.95)
    label("Talinolol gamma_inf in proximal cells (unitless)") # Table 1 note: talinolol gamma_inf,proximal cells (8.95)
    gamma_inf_dc <- fixed(7.01)
    label("Talinolol gamma_inf in distal cells (unitless)") # Table 1 note: talinolol gamma_inf,distal cells (7.01)
    gamma_inf_cc <- fixed(8.95)
    label("Talinolol gamma_inf in collecting-duct cells (unitless)") # Table 1 note: talinolol gamma_inf,collecting duct cells (8.95)
    gamma_eff_pc <- fixed(8.04)
    label("Talinolol gamma_eff in proximal cells (unitless)") # Table 1 note: talinolol gamma_eff,proximal cells (8.04)
    gamma_eff_dc <- fixed(5.07)
    label("Talinolol gamma_eff in distal cells (unitless)") # Table 1 note: talinolol gamma_eff,distal cells (5.07)
    gamma_eff_cc <- fixed(2.72)
    label("Talinolol gamma_eff in collecting-duct cells (unitless)") # Table 1 note: talinolol gamma_eff,collecting duct cells (2.72)
    km_pgp <- fixed(4071)
    label("Talinolol unbound Km for P-gp efflux, intestine / liver / kidney (ng/mL)") # Table 1 Talinolol K_m,u 4071 (P-gp); Results 11.2 umol/L
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

    # Rifampicin MW 822.94 g/mol converts the molar EC50,u to ng/mL
    # (Table S2C prints Ki as 402 ng/mL = 0.488 umol/L, the same ratio).
    ec50_rif <- ec50u_rif * 822.94

    # ------------------------------------------------------------------
    # 2. Rifampicin (Supplementary Material S1 'Case 1: rifampicin').
    #    Extended-clearance decomposition ('Other equations') with
    #    saturable OATP1B uptake; UGT autoinduction scales hepatic and
    #    enterocyte metabolism.
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
    # Segregated-flow intestine. Lumen-to-enterocyte and mucosal-blood-to-
    # enterocyte clearances are equal; the reverse directions divide by
    # gamma_inf (basolateral) or gamma_eff (apical).
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
    # 3. Induction turnover (Supplementary Material S1 Eq. 1), driven by the
    #    unbound rifampicin concentration in each hepatocyte unit, in each
    #    segment's enterocytes, or in blood for the kidney.
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
    d/dt(pgp_liver1) <- kdeg_liver * (1 + emax_pgp * uh1 / (ec50_rif + uh1) - pgp_liver1)
    d/dt(pgp_liver2) <- kdeg_liver * (1 + emax_pgp * uh2 / (ec50_rif + uh2) - pgp_liver2)
    d/dt(pgp_liver3) <- kdeg_liver * (1 + emax_pgp * uh3 / (ec50_rif + uh3) - pgp_liver3)
    d/dt(pgp_liver4) <- kdeg_liver * (1 + emax_pgp * uh4 / (ec50_rif + uh4) - pgp_liver4)
    d/dt(pgp_liver5) <- kdeg_liver * (1 + emax_pgp * uh5 / (ec50_rif + uh5) - pgp_liver5)
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
    pgp_liver1(0) <- 1
    pgp_liver2(0) <- 1
    pgp_liver3(0) <- 1
    pgp_liver4(0) <- 1
    pgp_liver5(0) <- 1
    pgp_duodenum(0) <- 1
    pgp_jejunum(0) <- 1
    pgp_ileum(0) <- 1
    pgp_kidney(0) <- 1

    # Competitive P-gp inhibition denominators (1 + unbound rifampicin / Ki):
    # hepatocyte unit, enterocytes of each segment, and blood as the
    # surrogate for the proximal tubular cells.
    inh_l1 <- 1 + uh1 / ki_pgp_rif
    inh_l2 <- 1 + uh2 / ki_pgp_rif
    inh_l3 <- 1 + uh3 / ki_pgp_rif
    inh_l4 <- 1 + uh4 / ki_pgp_rif
    inh_l5 <- 1 + uh5 / ki_pgp_rif
    inh_d <- 1 + ued / ki_pgp_rif
    inh_j <- 1 + uej / ki_pgp_rif
    inh_i <- 1 + uei / ki_pgp_rif
    inh_k <- 1 + ub / ki_pgp_rif

    # Relative P-gp activity (induced amount / inhibition), 1 = no
    # rifampicin; Figure 2c-d plots these for liver unit 1, the three gut
    # segments and the kidney.
    pgp_activity_liver1 <- pgp_liver1 / inh_l1
    pgp_activity_duodenum <- pgp_duodenum / inh_d
    pgp_activity_jejunum <- pgp_jejunum / inh_j
    pgp_activity_ileum <- pgp_ileum / inh_i
    pgp_activity_kidney <- pgp_kidney / inh_k

    # ------------------------------------------------------------------
    # 4. Talinolol (Supplementary Material S1 'Case 3: talinolol' and the
    #    unchanged parts of 'Case 2: digoxin').
    # ------------------------------------------------------------------
    kstomach <- exp(lkstomach * (1 - STUDY_TUBIC2006) + lkstomach_tubic2006 * STUDY_TUBIC2006)
    sfkp_study <- sfkp * (1 - STUDY_TUBIC2006) + sfkp_tubic2006 * STUDY_TUBIC2006
    fbclint_study <- fbclint_all * (1 - STUDY_TUBIC2006) + fbclint_all_tubic2006 * STUDY_TUBIC2006

    cl_int_all <- fbclint_study * WT / fb
    ps_act_inf <- 1 / (1 + rdif_liver) * cl_int_all / beta_liver
    ps_dif_inf <- rdif_liver / (1 + rdif_liver) * cl_int_all / beta_liver
    ps_dif_eff <- rdif_liver / (1 + rdif_liver) * cl_int_all / beta_liver / gamma_liver
    cl_hep <- rdif_liver / (1 + rdif_liver) * cl_int_all / (1 - beta_liver) / gamma_liver
    cl_bile_pgp <- cl_hep * f_bile * fpgp_liver
    cl_bile_nonpgp <- cl_hep * f_bile * (1 - fpgp_liver)
    cl_met <- cl_hep * (1 - f_bile)
    kp_muscle <- sfkp_study * exp(lkp_muscle)
    kp_skin <- sfkp_study * exp(lkp_skin)
    kp_adipose <- sfkp_study * exp(lkp_adipose)
    kp_serosa <- sfkp_study * exp(lkp_serosa)

    # Intestine: passive clearance per segment scales with surface area;
    # active uptake is PSdif / Rdif,inf,gut (twice passive); apical P-gp
    # efflux is PSdif / gamma_eff / Rdif,eff,segment.
    psd <- ps_dif_gut * WT * fsa_duo
    psj <- ps_dif_gut * WT * fsa_jej
    psi <- ps_dif_gut * WT * fsa_ile
    pgpd <- psd / gamma_ent / rdif_eff_duo * fpgp_ent
    pgpj <- psj / gamma_ent / rdif_eff_jej * fpgp_ent
    pgpi <- psi / gamma_ent / rdif_eff_ile * fpgp_ent
    npgpd <- psd / gamma_ent / rdif_eff_duo * (1 - fpgp_ent)
    npgpj <- psj / gamma_ent / rdif_eff_jej * (1 - fpgp_ent)
    npgpi <- psi / gamma_ent / rdif_eff_ile * (1 - fpgp_ent)

    # Kidney (Suppl. S1 'Other equations'): basolateral and apical passive
    # clearances of each tubular segment, the basolateral uptake transporter
    # and the apical P-gp of the proximal cells.
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
    c_he1 <- is_liver1 / vhe5
    c_he2 <- is_liver2 / vhe5
    c_he3 <- is_liver3 / vhe5
    c_he4 <- is_liver4 / vhe5
    c_he5 <- is_liver5 / vhe5
    c_hc1 <- int_liver1 / vhc5
    c_hc2 <- int_liver2 / vhc5
    c_hc3 <- int_liver3 / vhc5
    c_hc4 <- int_liver4 / vhc5
    c_hc5 <- int_liver5 / vhc5
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

    # Hepatic uptake / efflux / metabolism and saturable, induced and
    # inhibited P-gp biliary excretion of each unit.
    upt1 <- fb * (ps_act_inf + ps_dif_inf) / 5 * c_he1
    upt2 <- fb * (ps_act_inf + ps_dif_inf) / 5 * c_he2
    upt3 <- fb * (ps_act_inf + ps_dif_inf) / 5 * c_he3
    upt4 <- fb * (ps_act_inf + ps_dif_inf) / 5 * c_he4
    upt5 <- fb * (ps_act_inf + ps_dif_inf) / 5 * c_he5
    eff1 <- fh * ps_dif_eff / 5 * c_hc1
    eff2 <- fh * ps_dif_eff / 5 * c_hc2
    eff3 <- fh * ps_dif_eff / 5 * c_hc3
    eff4 <- fh * ps_dif_eff / 5 * c_hc4
    eff5 <- fh * ps_dif_eff / 5 * c_hc5
    bile1 <- fh * (cl_bile_pgp / 5 * km_pgp * pgp_liver1 / (km_pgp * inh_l1 + fh * c_hc1) + cl_bile_nonpgp / 5) * c_hc1
    bile2 <- fh * (cl_bile_pgp / 5 * km_pgp * pgp_liver2 / (km_pgp * inh_l2 + fh * c_hc2) + cl_bile_nonpgp / 5) * c_hc2
    bile3 <- fh * (cl_bile_pgp / 5 * km_pgp * pgp_liver3 / (km_pgp * inh_l3 + fh * c_hc3) + cl_bile_nonpgp / 5) * c_hc3
    bile4 <- fh * (cl_bile_pgp / 5 * km_pgp * pgp_liver4 / (km_pgp * inh_l4 + fh * c_hc4) + cl_bile_nonpgp / 5) * c_hc4
    bile5 <- fh * (cl_bile_pgp / 5 * km_pgp * pgp_liver5 / (km_pgp * inh_l5 + fh * c_hc5) + cl_bile_nonpgp / 5) * c_hc5

    # Enterocyte apical efflux to the lumen of each segment: passive plus
    # saturable, induced and inhibited P-gp plus non-P-gp transport.
    effd <- fent * (psd / gamma_ent + pgpd * km_pgp * pgp_duodenum / (km_pgp * inh_d + fent * c_ed) + npgpd) * c_ed
    effj <- fent * (psj / gamma_ent + pgpj * km_pgp * pgp_jejunum / (km_pgp * inh_j + fent * c_ej) + npgpj) * c_ej
    effi <- fent * (psi / gamma_ent + pgpi * km_pgp * pgp_ileum / (km_pgp * inh_i + fent * c_ei) + npgpi) * c_ei

    # Proximal-cell apical efflux to the tubular lumen.
    eff_pc <- frenal * (pgp_pc * km_pgp * pgp_kidney / (km_pgp * inh_k + frenal * c_pt_c) + npgp_pc + ps_pct) * c_pt_c

    d/dt(central) <- qh * (c_he5 - c_b) +
      qm * (c_m / kp_muscle - c_b) +
      qs * (c_s / kp_skin - c_b) +
      qa * (c_a / kp_adipose - c_b) +
      (qr - qu_out) * c_cd_b - qr * c_b
    d/dt(is_liver1) <- qhart * c_b + qmd * c_md + qmj * c_mj + qmi * c_mi +
      qser * c_ser / kp_serosa - qh * c_he1 + eff1 - upt1
    d/dt(is_liver2) <- qh * (c_he1 - c_he2) + eff2 - upt2
    d/dt(is_liver3) <- qh * (c_he2 - c_he3) + eff3 - upt3
    d/dt(is_liver4) <- qh * (c_he3 - c_he4) + eff4 - upt4
    d/dt(is_liver5) <- qh * (c_he4 - c_he5) + eff5 - upt5
    d/dt(int_liver1) <- upt1 - eff1 - fh * cl_met / 5 * c_hc1 - bile1
    d/dt(int_liver2) <- upt2 - eff2 - fh * cl_met / 5 * c_hc2 - bile2
    d/dt(int_liver3) <- upt3 - eff3 - fh * cl_met / 5 * c_hc3 - bile3
    d/dt(int_liver4) <- upt4 - eff4 - fh * cl_met / 5 * c_hc4 - bile4
    d/dt(int_liver5) <- upt5 - eff5 - fh * cl_met / 5 * c_hc5 - bile5
    d/dt(ehc1) <- bile1 + bile2 + bile3 + bile4 + bile5 - k_bile * ehc1
    d/dt(ehc2) <- k_bile * (ehc1 - ehc2)
    d/dt(ehc3) <- k_bile * (ehc2 - ehc3)
    d/dt(muscle) <- qm * (c_b - c_m / kp_muscle)
    d/dt(skin) <- qs * (c_b - c_s / kp_skin)
    d/dt(adipose) <- qa * (c_b - c_a / kp_adipose)
    d/dt(serosa) <- qser * (c_b - c_ser / kp_serosa)
    d/dt(duodenum_muc) <- qmd * (c_b - c_md) + fent * psd / gamma_ent * c_ed - fb * psd * c_md
    d/dt(jejunum_muc) <- qmj * (c_b - c_mj) + fent * psj / gamma_ent * c_ej - fb * psj * c_mj
    d/dt(ileum_muc) <- qmi * (c_b - c_mi) + fent * psi / gamma_ent * c_ei - fb * psi * c_mi
    d/dt(duodenum_ent) <- (psd / rdif_inf_gut + psd) * c_ld + fb * psd * c_md - fent * psd / gamma_ent * c_ed - effd
    d/dt(jejunum_ent) <- (psj / rdif_inf_gut + psj) * c_lj + fb * psj * c_mj - fent * psj / gamma_ent * c_ej - effj
    d/dt(ileum_ent) <- (psi / rdif_inf_gut + psi) * c_li + fb * psi * c_mi - fent * psi / gamma_ent * c_ei - effi
    d/dt(stomach) <- -kstomach * stomach
    d/dt(duodenum_lumen) <- k_bile * ehc3 + kstomach * stomach -
      (psd / rdif_inf_gut + psd) * c_ld - k_feces_duo * duodenum_lumen + effd
    d/dt(jejunum_lumen) <- k_feces_duo * duodenum_lumen -
      (psj / rdif_inf_gut + psj) * c_lj - k_feces_jej * jejunum_lumen + effj
    d/dt(ileum_lumen) <- k_feces_jej * jejunum_lumen -
      (psi / rdif_inf_gut + psi) * c_li - k_feces_ile * ileum_lumen + effi
    d/dt(a_feces) <- k_feces_ile * ileum_lumen

    # Kidney: glomerulus, then proximal tubule, distal tubule and collecting
    # duct, each with blood vessel, cells and urinary lumen.
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
    # oral talinolol to `stomach`; intravenous talinolol to `central` (give
    # the 30-min infusion with `dur`).
    alag(duodenum_lumen_rif) <- exp(ltlag_rif)

    # ------------------------------------------------------------------
    # 5. Observations: blood concentrations (ng/mL), as fitted and reported
    #    in the paper. No residual-error model is published (Napp nonlinear
    #    least squares), so the model is deterministic.
    # ------------------------------------------------------------------
    Cc <- c_b
    Cc_rif <- c_b_rif
  })
}
