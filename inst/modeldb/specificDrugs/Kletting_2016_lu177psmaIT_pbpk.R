Kletting_2016_lu177psmaIT_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, SAAM II / Simulink). Physiologically based",
    "pharmacokinetic model of the PSMA-targeting radioligand",
    "[177Lu]Lu-PSMA I&T in men with metastasizing prostate cancer, built to",
    "find the peptide amount and injected activity that maximise the",
    "tumour-to-kidney and tumour-to-salivary-gland biologically effective",
    "dose ratio. Labelled and unlabelled peptide are carried in two parallel",
    "circulations with identical physiology, coupled by competition for the",
    "same finite pool of PSMA binding sites and by the physical decay of",
    "177Lu, which converts a labelled molecule into an unlabelled one in",
    "place. Nine PSMA-positive regions (two delineated tumour lesions, a",
    "lumped tumour remainder, salivary glands, liver, spleen, GI tract plus",
    "pancreas, prostate and kidneys) carry vascular, interstitial,",
    "receptor-bound and internalised sub-compartments; the kidneys add",
    "glomerular filtration, tubular excretion and an unspecific",
    "proximal-tubular-cell pool. Eight PSMA-negative tissues carry vascular",
    "and interstitial spaces, the brain a vascular space only, and arteries,",
    "veins and a serum-protein-bound pool complete the circulation.",
    "Ten parameters were fitted per patient to planar gamma-camera",
    "time-activity curves; the shipped values are the paper's typical",
    "patient of Figure 3 (median normal-tissue parameters, tumour 1 = median",
    "and tumour 2 = mean tumour parameters of Table 2). Amounts are nmol and",
    "time is minutes.",
    sep = " "
  )
  reference <- paste(
    "Kletting P, Schuchardt C, Kulkarni HR, Shahinfar M, Singh A, Glatting G,",
    "Baum RP, Beer AJ. Investigating the Effect of Ligand Amount and Injected",
    "Therapeutic Activity: A Simulation Study for 177Lu-Labeled",
    "PSMA-Targeting Peptides. PLoS One. 2016;11(9):e0162303.",
    "doi:10.1371/journal.pone.0162303.",
    sep = " "
  )
  vignette <- "Kletting_2016_lu177psmaIT_pbpk"
  dosing <- c("pven_lab", "pven_unl")
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/L")

  paper_specific_compartment_pattern <- "^(pv|pint|rp|pintern|pintra|part|pven|pprp|urine|cleared)_"

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Scales ten ICRP Reference Man organ volumes, each written in S1 Table A",
        "as coefficient * BW/71, and sets the rest-of-body volume V_REST =",
        "V_BW - sum of organ volumes (1 mL = 1 g). S1 Table A lists BW as",
        "measured, but the per-patient values are not published in the paper",
        "or its supplements. 71 kg, the reference individual the BW/71",
        "normalisation encodes, reproduces Reference Man exactly."
      ),
      source_name = "BW"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "% (volume fraction x 100)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as the fraction H = HCT/100 (S1 Table A carries H as unity).",
        "Sets the total serum volume V_P = 2.8*(1-H)*BSA and hence the total",
        "serum flow and every vascular sub-volume, and also v_TU,v = 0.05*(1-H),",
        "f_PRO = 0.18*(1-H), V_PRO,v = 0.004*(1-H)*V_PRO,total and V_SAL,v =",
        "0.03*(1-H)*V_SAL,total. Listed as measured in S1 Table A; per-patient",
        "values are not published."
      ),
      source_name = "H"
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured per patient, Table 1 (P1-P5: 2.0, 1.9, 1.8, 2.1, 2.0 m^2).",
        "Used directly rather than recomputed with the Du Bois formula BSA =",
        "0.007184*BH^0.725*BW^0.425 that S1 Table A also lists, because body",
        "height and weight are not published."
      ),
      source_name = "BSA"
    ),
    TER_MAG3 = list(
      description = "Tubular extraction rate measured by 99mTc-MAG3 renography",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured per patient, Table 1 (P1-P5: 198, 201, 136, 252, 176 mL/min).",
        "S1 Table A derives the 51Cr-EDTA glomerular filtration rate as GFR =",
        "TER/3 * 20/15, which is then scaled to the molecular size of PSMA I&T",
        "by the sieving-coefficient ratio phi = 0.66."
      ),
      source_name = "TER"
    ),
    ORGVOL_KIDNEY = list(
      description = "Measured total kidney volume",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CT-measured per patient, left plus right kidney, Table 1 (P1-P5: 321,",
        "311, 394, 268, 296 mL). Split into vascular (5.5 percent), interstitial",
        "(15 percent) and intracellular-tubular (two thirds of the remainder)",
        "sub-volumes."
      ),
      source_name = "Kidneys (measured volume)"
    ),
    ORGVOL_SALGLAND = list(
      description = "Measured total salivary gland volume (left plus right parotid)",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CT-measured per patient, left plus right parotid gland, Table 1",
        "(P1-P5: 54, 21, 17, 52, 29 mL). The parotids are taken as",
        "representative of all salivary glands."
      ),
      source_name = "Parotid glands (measured volume)"
    ),
    ORGVOL_TUMOR1 = list(
      description = "Measured volume of explicitly modelled tumour lesion 1",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CT-measured per patient, Table 1 (P1-P5: 0.5, 1, 2, 4, 1.5 mL)."
      ),
      source_name = "Tumor 1 (measured volume)"
    ),
    ORGVOL_TUMOR2 = list(
      description = "Measured volume of explicitly modelled tumour lesion 2",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "CT-measured per patient, Table 1 (P1-P5: 1, 34, 13, 3, 1 mL)."
      ),
      source_name = "Tumor 2 (measured volume)"
    ),
    ORGVOL_TUMORREST = list(
      description = "Assumed volume of the lumped remaining-tumour-lesion compartment",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Not measured. S1 Table A defines the tumour-remainder volumes through",
        "the binding-site number, V = R_TU,Rest,0 / [R_TU,Rest,0], with",
        "[R_TU,Rest,0] = 266 nmol/L and R_TU,Rest,0 = 13, 2.6, 13, 2.6, 13 nmol",
        "for P1-P5 (footnote c: a 10 mL or 50 mL additional tumour volume).",
        "The quotient gives 48.87 or 9.774 mL."
      ),
      source_name = "R_TU,Rest,0 / [R_TU,Rest,0]"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Reported per patient in Table 1 (53-78 years) as a cohort descriptor; not used in any model equation.",
      source_name = "Age"
    )
  )

  compartmentData <- list(
    pv_tumor1_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumor2_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumorrest_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_salgland_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_liver_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_spleen_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_gi_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_prostate_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_kidney_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_muscle_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_lungs_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_bone_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_redmarrow_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_skin_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_adipose_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_heart_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_rest_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_brain_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pint_tumor1_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumor2_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumorrest_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_salgland_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_liver_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_spleen_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_gi_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_prostate_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_kidney_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_muscle_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_lungs_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_bone_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_redmarrow_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_skin_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_adipose_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_heart_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_rest_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_tumor1_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor1_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumor2_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor2_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumorrest_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumorrest_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_salgland_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_salgland_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_liver_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_liver_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_spleen_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_spleen_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_gi_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_gi_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_prostate_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_prostate_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_kidney_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_kidney_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintra_kidney_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    part_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pven_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pprp_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    urine_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "urine", verified = TRUE),
    cleared_lab = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "not applicable", verified = TRUE),
    pv_tumor1_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumor2_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumorrest_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_salgland_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_liver_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_spleen_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_gi_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_prostate_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_kidney_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_muscle_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_lungs_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_bone_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_redmarrow_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_skin_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_adipose_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_heart_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_rest_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_brain_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pint_tumor1_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumor2_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumorrest_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_salgland_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_liver_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_spleen_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_gi_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_prostate_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_kidney_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_muscle_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_lungs_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_bone_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_redmarrow_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_skin_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_adipose_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_heart_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_rest_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_tumor1_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor1_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumor2_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor2_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumorrest_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumorrest_unl = list(
      analyte = "PSMA I&T (unlabelled)",
      units = "nmol",
      specimen = "tumor",
      verified = TRUE
    ),
    rp_salgland_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_salgland_unl = list(
      analyte = "PSMA I&T (unlabelled)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    rp_liver_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_liver_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_spleen_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_spleen_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_gi_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_gi_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_prostate_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_prostate_unl = list(
      analyte = "PSMA I&T (unlabelled)",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    rp_kidney_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_kidney_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintra_kidney_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    part_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pven_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pprp_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    urine_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "urine", verified = TRUE),
    cleared_unl = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 5,
    n_studies = 1,
    age_range = "53-78 years",
    age_mean = "66 years",
    sex_female_pct = 0,
    disease_state = "metastasizing prostate cancer (PSMA-positive on 68Ga-PSMA PET/CT) treated with 177Lu-PSMA I&T",
    dose_range = paste(
      "one therapy cycle per patient: 74-302 nmol total PSMA I&T (mean 142",
      "nmol, 235 +/- 155 ug) labelled with 5.4-6.0 GBq 177Lu, 10 min",
      "intravenous infusion, with lysine/arginine co-infusion from 0.5 h",
      "before; the simulation study spans 2 to 8192 nmol",
      sep = " "
    ),
    regions = "Germany (Zentralklinik Bad Berka)",
    notes = paste(
      "Baseline measurements in Table 1: age, body surface area,",
      "99mTc-MAG3 tubular extraction rate, CT volumes of the parotid glands,",
      "kidneys and two tumour lesions (bone and lymph node), injected amount",
      "and activity. Planar whole-body gamma-camera images at 0.5 h, 2 h and",
      "1-5 days post-infusion (4-5 time points per patient); patient 1",
      "contributes his second cycle, patients 2-5 their first.",
      sep = " "
    )
  )

  ini({
    # ---- Binding, internalisation, protein binding, decay (S1 Table A) ----
    kon        <- fixed(0.046)    ; label("Association rate onto PSMA (L/nmol/min)")           # S1 Table A; Methods 'PBPK Model validation and selection', k_on = 0.046 l/nmol/min
    kd         <- fixed(1)        ; label("PSMA dissociation constant (nmol/L)")               # Model 1: k_off = 0.046 /min (Methods), K_D = k_off/k_on = 1 nmol/L (S1 Table A lists 1 or 8); Model 1 had Akaike weight > 0.9 in every patient (Results)
    kint       <- fixed(0.001)    ; label("Internalisation rate, all PSMA-positive regions (1/min)") # Methods: lambda_i,int = 0.001 /min for all organs and tumour; S1 Table A
    kpr        <- fixed(4.7e-4)   ; label("Binding rate of peptide to serum protein (1/min)")  # S1 Table A, k_PR = 4.7e-4 /min
    lambda_phy <- fixed(7.15e-5)  ; label("Physical decay rate of 177Lu (1/min)")              # S1 Table A, lambda_phy = 7.15e-5 /min

    # ---- Release rates (Table 2, typical patient of Figure 3) -------------
    krel_kidney   <- fixed(2.9e-4)  ; label("Release rate, kidneys and other PSMA-positive normal tissue (1/min)") # Table 2 Median 2.9 (x1e-4 /min); S1 Table A sets lambda_L, lambda_S and lambda_NT,release = lambda_K,release
    krel_salgland <- fixed(3.7e-4)  ; label("Release rate, salivary glands (1/min)")               # Table 2 Median 3.7 (x1e-4 /min)
    krel_tumor1   <- fixed(2.0e-4)  ; label("Release rate, tumour lesion 1 (1/min)")               # Table 2 Median 2.0 (x1e-4 /min): the median tumour of Figure 3
    krel_tumor2   <- fixed(1.76e-4) ; label("Release rate, tumour lesion 2 (1/min)")              # Table 2 Mean 1.76 (x1e-4 /min): the mean tumour of Figure 3

    # ---- Binding-site densities (nmol/L) ----------------------------------
    rdens_kidney    <- fixed(24)    ; label("PSMA binding-site density, kidneys (nmol/L)")         # Table 2 Median [R_K,0] = 24
    rdens_salgland  <- fixed(42)    ; label("PSMA binding-site density, salivary glands (nmol/L)")  # Table 2 Median [R_SAL,0] = 42
    rdens_tumor1    <- fixed(110)   ; label("PSMA binding-site density, tumour lesion 1 (nmol/L)")  # Table 2 Median [R_TU,0] = 110: the median tumour of Figure 3
    rdens_tumor2    <- fixed(380)   ; label("PSMA binding-site density, tumour lesion 2 (nmol/L)")  # Table 2 Mean [R_TU,0] = 380: the mean tumour of Figure 3
    rdens_tumorrest <- fixed(266)   ; label("PSMA binding-site density, tumour remainder (nmol/L)") # S1 Table A, [R_TU,Rest,0] = 266
    rdens_prostate  <- fixed(26.6)  ; label("PSMA binding-site density, prostate (nmol/L)")        # S1 Table A: [R_PRO,0] = [R_TU,Rest,0] * 0.1
    rdens_liver     <- fixed(1.33)  ; label("PSMA binding-site density, liver (nmol/L)")           # S1 Table A: [R_L,0] = [R_PRO,0] * 0.05
    rdens_spleen    <- fixed(0.532) ; label("PSMA binding-site density, spleen (nmol/L)")          # S1 Table A: [R_S,0] = [R_PRO,0] * 0.02
    rdens_gi        <- fixed(1.596) ; label("PSMA binding-site density, GI plus pancreas (nmol/L)") # S1 Table A: [R_GI,0] = [R_PRO,0] * 0.06

    # ---- Perfusion and permeability densities (mL/min/g) ------------------
    fperf_tumor     <- fixed(0.5)   ; label("Serum flow density, tumour (mL/min/g)")                # Methods: f_TU = 0.5 ml/g/min assumed; S1 Table A
    fperf_prostate  <- fixed(0.18)  ; label("Serum flow density, prostate, before the (1-H) factor (mL/min/g)") # S1 Table A: f_PRO = 0.18*(1-H)
    fperf_salgland  <- fixed(0.16)  ; label("Serum flow density, salivary glands (mL/min/g)")       # Table 2 Median f_SAL = 0.16
    kps_tumor       <- fixed(0.6)   ; label("Permeability-surface density, tumour (mL/min/g)")      # S1 Table A, k_TU = 0.6
    kps_muscle      <- fixed(0.02)  ; label("Permeability-surface density, muscle-like tissue (mL/min/g)") # S1 Table A, k_MUS = 0.02; k_GI, k_SKIN, k_ADI, k_HRT, k_BONE, k_REST are set equal to it
    kps_prostate    <- fixed(0.1)   ; label("Permeability-surface density, prostate (mL/min/g)")    # S1 Table A, k_PRO = 0.1
    kps_liver       <- fixed(2)     ; label("Permeability-surface density, liver-like tissue (mL/min/g)") # S1 Table A: k_L = k_MUS * 100; k_S, k_LU, k_SAL, k_RM are set equal to it

    # ---- Renal handling ---------------------------------------------------
    phi_gfr <- fixed(0.66) ; label("Ratio of sieving coefficients, PSMA I&T to 51Cr-EDTA (unitless)") # S1 Table A, phi = 0.66
    f_ex    <- fixed(0.96) ; label("Tubular excretion fraction of the filtered load (unitless)")      # S1 Table A, f_ex = 0.96

    # ---- Tumour ROI background correction (observation only) --------------
    bgcorr_tumor1 <- fixed(0) ; label("Muscle background fraction in the tumour 1 ROI, paper symbol x_tu,1 (unitless)") # Fitted per patient (S2 File) to absorb ROI background; not summarised in Table 2; 0 gives the lesion content used for dosimetry
    bgcorr_tumor2 <- fixed(0) ; label("Muscle background fraction in the tumour 2 ROI, paper symbol x_tu,2 (unitless)") # Fitted per patient (S2 File) to absorb ROI background; not summarised in Table 2; 0 gives the lesion content used for dosimetry
  })

  model({
    # =====================================================================
    # 1. Anatomy and physiology. All volumes in L, all flows in L/min, all
    #    amounts in nmol, time in min. Because the paper reports density
    #    parameters in mL/min/g and assumes rho = 1 g/mL (S1 Table A), a flow
    #    or PS product in L/min is numerically the density times the organ
    #    volume in L, with no conversion factor.
    # =====================================================================
    hct_frac <- HCT / 100                      # S1 Table A carries H as a fraction
    vserum   <- 2.8 * (1 - hct_frac) * BSA     # V_P, S1 Table A
    fserum   <- vserum * 1.23                  # F, S1 Table A and footnote b (6500/5300 mL/min per mL)

    # ---- Total organ volumes (S1 Table A) ----------------------------------
    vtot_kidney     <- ORGVOL_KIDNEY / 1000    # measured, Table 1
    vtot_salgland   <- ORGVOL_SALGLAND / 1000  # measured, Table 1
    vtot_tumor1     <- ORGVOL_TUMOR1 / 1000    # measured, Table 1
    vtot_tumor2     <- ORGVOL_TUMOR2 / 1000    # measured, Table 1
    vtot_tumorrest  <- ORGVOL_TUMORREST / 1000 # S1 Table A footnote c
    vtot_liver      <- 1.91                    # S1 Table A
    vtot_spleen     <- 0.183                   # S1 Table A
    vtot_prostate   <- 0.016  * WT / 71        # S1 Table A
    vtot_lungs      <- 1      * WT / 71        # S1 Table A
    vtot_muscle     <- 30.078 * WT / 71        # S1 Table A
    vtot_gi         <- (0.385 + 0.548 + 0.104 + 0.15) * WT / 71  # GI plus pancreas, S1 Table A
    vtot_skin       <- 3.408  * WT / 71        # S1 Table A
    vtot_adipose    <- 13.465 * WT / 71        # S1 Table A
    vtot_redmarrow  <- 1.1    * WT / 71        # S1 Table A
    vtot_bone       <- 10.165 * WT / 71 - vtot_redmarrow  # S1 Table A
    vtot_heart      <- 0.341  * WT / 71        # S1 Table A
    vtot_brain      <- 1.45   * WT / 71        # S1 Table A
    # V_REST,total = V_BW - sum over i of V_i,total, i = all organs except
    # tumour (S1 Table A), with V_BW = BW at 1 mL = 1 g.
    vtot_rest <- WT -
      vtot_salgland -
      vtot_liver -
      vtot_spleen -
      vtot_gi -
      vtot_prostate -
      vtot_kidney -
      vtot_muscle -
      vtot_lungs -
      vtot_bone -
      vtot_redmarrow -
      vtot_skin -
      vtot_adipose -
      vtot_heart -
      vtot_brain

    # ---- Vascular (serum) sub-volumes (S1 Table A) --------------------------
    vv_lungs     <- 0.105 * vserum
    vv_muscle    <- 0.14  * vserum
    vv_gi        <- 0.076 * vserum
    vv_skin      <- 0.03  * vserum
    vv_adipose   <- 0.05  * vserum
    vv_redmarrow <- 0.04  * vserum
    vv_bone      <- 0.07  * vserum - vv_redmarrow
    vv_heart     <- 0.01  * vserum
    vv_brain     <- 0.012 * vserum
    vv_art       <- 0.105 * vserum            # 0.06 arterial plus half the cardiac serum
    vv_ven       <- 0.225 * vserum            # 0.18 venous plus half the cardiac serum
    vv_liver     <- 0.085 * vtot_liver
    vv_spleen    <- 0.12  * vtot_spleen
    vv_kidney    <- 0.055 * vtot_kidney
    vv_prostate  <- 0.004 * (1 - hct_frac) * vtot_prostate
    vv_salgland  <- 0.03  * (1 - hct_frac) * vtot_salgland
    vv_tumor1    <- 0.05  * (1 - hct_frac) * vtot_tumor1     # v_TU,v
    vv_tumor2    <- 0.05  * (1 - hct_frac) * vtot_tumor2     # v_TU,v
    vv_tumorrest <- 0.05  * (1 - hct_frac) * vtot_tumorrest  # v_TU,v
    # V_REST,v = V_P - sum over i of V_i,v, i = all organs except tumour
    # (S1 Table A); the arterial and venous pools are part of V_P and are
    # subtracted with the organs.
    vv_rest <- vserum - vv_art - vv_ven -
      vv_salgland -
      vv_liver -
      vv_spleen -
      vv_gi -
      vv_prostate -
      vv_kidney -
      vv_muscle -
      vv_lungs -
      vv_bone -
      vv_redmarrow -
      vv_skin -
      vv_adipose -
      vv_heart -
      vv_brain

    # ---- Interstitial sub-volumes (S1 Table A) ------------------------------
    vi_liver     <- 0.2  * vtot_liver
    vi_spleen    <- 0.2  * vtot_spleen
    vi_kidney    <- 0.15 * vtot_kidney
    vi_prostate  <- 0.25 * vtot_prostate
    vi_salgland  <- 0.23 * vtot_salgland
    vi_tumor1    <- 0.38 * vtot_tumor1        # v_TU,int
    vi_tumor2    <- 0.38 * vtot_tumor2        # v_TU,int
    vi_tumorrest <- 0.38 * vtot_tumorrest     # v_TU,int
    vi_lungs      <- vv_lungs      * 5.5         # alpha_LUNGS
    vi_muscle     <- vv_muscle     * 5.9         # alpha_MUSCLE
    vi_gi         <- vv_gi         * 8.8         # alpha_GI
    vi_skin       <- vv_skin       * 8.9         # alpha_SKIN
    vi_adipose    <- vv_adipose    * 15.5        # alpha_ADIPOSE
    vi_redmarrow  <- vv_redmarrow  * 3.7         # alpha_REDMARROW
    vi_bone       <- vv_bone       * 8.4         # alpha_BONE
    vi_heart      <- vv_heart      * 3.7         # alpha_HEART
    vi_rest       <- vv_rest       * 4.1         # alpha_REST
    # Intracellular kidney: two thirds of what is neither vascular nor interstitial
    vintra_kidney <- (vtot_kidney - vi_kidney - vv_kidney) * 2 / 3  # S1 Table A footnote d

    # ---- Serum flows (S1 Table A) -------------------------------------------
    fq_liver      <- 0.065  * fserum
    fq_spleen     <- 0.03   * fserum
    fq_kidney     <- 0.19   * fserum
    fq_muscle     <- 0.17   * fserum
    fq_gi         <- 0.16   * fserum
    fq_skin       <- 0.05   * fserum
    fq_adipose    <- 0.05   * fserum
    fq_redmarrow  <- 0.03   * fserum
    fq_bone       <- 0.05   * fserum
    fq_heart      <- 0.04   * fserum
    fq_brain      <- 0.12   * fserum
    fq_prostate  <- fperf_prostate * (1 - hct_frac) * vtot_prostate  # f_PRO * V_PRO,total
    fq_salgland  <- fperf_salgland * vtot_salgland                   # f_SAL * V_SAL,total (f_SAL fitted)
    fq_tumor1    <- fperf_tumor * vtot_tumor1                        # f_TU * V_TU,total
    fq_tumor2    <- fperf_tumor * vtot_tumor2
    fq_tumorrest <- fperf_tumor * vtot_tumorrest
    # The tabulated fractions of F account for 95.5 % of systemic flow; the
    # remainder, less prostate and salivary glands, perfuses the rest compartment.
    fq_rest <- fserum - fq_prostate - fq_salgland -
      fq_liver -
      fq_spleen -
      fq_kidney -
      fq_muscle -
      fq_gi -
      fq_skin -
      fq_adipose -
      fq_redmarrow -
      fq_bone -
      fq_heart -
      fq_brain
    # F_TOTAL = F + F_TU,1 + F_TU,2 + F_TU,REST (S1 Table A): the cardiac
    # serum output that returns through the veins and lungs, tumours included.
    fq_total <- fserum + fq_tumor1 + fq_tumor2 + fq_tumorrest

    # ---- Permeability-surface-area products (S1 Table A) --------------------
    ps_tumor1     <- kps_tumor    * vtot_tumor1
    ps_tumor2     <- kps_tumor    * vtot_tumor2
    ps_tumorrest  <- kps_tumor    * vtot_tumorrest
    ps_salgland   <- kps_liver    * vtot_salgland
    ps_liver      <- kps_liver    * vtot_liver
    ps_spleen     <- kps_liver    * vtot_spleen
    ps_lungs      <- kps_liver    * vtot_lungs
    ps_redmarrow  <- kps_liver    * vtot_redmarrow
    ps_prostate   <- kps_prostate * vtot_prostate
    ps_gi         <- kps_muscle   * vtot_gi
    ps_muscle     <- kps_muscle   * vtot_muscle
    ps_skin       <- kps_muscle   * vtot_skin
    ps_adipose    <- kps_muscle   * vtot_adipose
    ps_heart      <- kps_muscle   * vtot_heart
    ps_bone       <- kps_muscle   * vtot_bone
    ps_rest       <- kps_muscle   * vtot_rest
    # Brain PS = 0 (S1 Eq 4 note); the brain has no interstitial state.

    # ---- Renal filtration and excretion (S1 Table A) ------------------------
    gfr    <- (TER_MAG3 / 1000) / 3 * 20 / 15  # 51Cr-EDTA GFR from measured TER
    f_fil  <- gfr * phi_gfr                    # filtration scaled for molecular size
    f_exc  <- f_fil * f_ex                     # excreted fraction of the filtered load
    f_reab <- f_fil - f_exc                    # returned to the vasculature via the cells

    # ---- Total binding sites, R_i,0 = [R_i,0] * V_i,total (S1 Table A) ------
    r0_tumor1     <- rdens_tumor1     * vtot_tumor1
    r0_tumor2     <- rdens_tumor2     * vtot_tumor2
    r0_tumorrest  <- rdens_tumorrest  * vtot_tumorrest
    r0_salgland   <- rdens_salgland   * vtot_salgland
    r0_liver      <- rdens_liver      * vtot_liver
    r0_spleen     <- rdens_spleen     * vtot_spleen
    r0_gi         <- rdens_gi         * vtot_gi
    r0_prostate   <- rdens_prostate   * vtot_prostate
    r0_kidney     <- rdens_kidney     * vtot_kidney

    # ---- Rate constants ---------------------------------------------------
    koff <- kd * kon                           # S1 Table A: k_off = K_D * k_on
    krel_tumorrest  <- (krel_tumor1 + krel_tumor2) / 2  # S1 Table A: lambda_TU,Rest,release
    krel_liver      <- krel_kidney
    krel_spleen     <- krel_kidney
    krel_gi         <- krel_kidney
    krel_prostate   <- krel_kidney

    # =====================================================================
    # 2. ODE system. S1 Eqs 1-13. Each state exists twice: _lab is
    #    radiolabelled [177Lu]Lu-PSMA I&T and _unl is unlabelled peptide.
    #    Physical decay (lambda_phy) is the only coupling between the two
    #    circulations: it removes labelled peptide from a state and creates
    #    unlabelled peptide in the SAME state.
    # =====================================================================

    # ---------------------------------------------------------------------
    # Radiolabelled circulation
    # ---------------------------------------------------------------------

    # Free binding sites, S1 Eq 1: R_0,i = R_i + RP_i + RP*_i.
    # The pool is shared, so labelled and unlabelled peptide compete.
    rfree_tumor1     <- r0_tumor1 - rp_tumor1_lab - rp_tumor1_unl
    rfree_tumor2     <- r0_tumor2 - rp_tumor2_lab - rp_tumor2_unl
    rfree_tumorrest  <- r0_tumorrest - rp_tumorrest_lab - rp_tumorrest_unl
    rfree_salgland   <- r0_salgland - rp_salgland_lab - rp_salgland_unl
    rfree_liver      <- r0_liver - rp_liver_lab - rp_liver_unl
    rfree_spleen     <- r0_spleen - rp_spleen_lab - rp_spleen_unl
    rfree_gi         <- r0_gi - rp_gi_lab - rp_gi_unl
    rfree_prostate   <- r0_prostate - rp_prostate_lab - rp_prostate_unl
    rfree_kidney     <- r0_kidney - rp_kidney_lab - rp_kidney_unl

    # Arteries, S1 Eq 8. The printed outflow term multiplies P_i,v where P_ART is
    # meant, and the printed inflow uses F; F_TOTAL is used for both so that the
    # arterial outflow to every organ (tumours included) balances the inflow.
    d/dt(part_lab) <- fq_total * (pv_lungs_lab / vv_lungs - part_lab / vv_art) - lambda_phy * part_lab

    # Veins, S1 Eq 7, plus the outflow to the lungs that the printed equation
    # omits (F_TOTAL, as in S1 Eq 5). Spleen and GI drain into the liver rather than
    # into the veins, so the liver returns F_L + F_S + F_GI (S1 compartment diagrams B2).
    d/dt(pven_lab) <- - kpr * pven_lab - fq_total * pven_lab / vv_ven +
      (fq_liver + fq_spleen + fq_gi) * pv_liver_lab / vv_liver +
      fq_tumor1 * pv_tumor1_lab / vv_tumor1 +
      fq_tumor2 * pv_tumor2_lab / vv_tumor2 +
      fq_tumorrest * pv_tumorrest_lab / vv_tumorrest +
      fq_salgland * pv_salgland_lab / vv_salgland +
      fq_prostate * pv_prostate_lab / vv_prostate +
      fq_kidney * pv_kidney_lab / vv_kidney +
      fq_muscle * pv_muscle_lab / vv_muscle +
      fq_bone * pv_bone_lab / vv_bone +
      fq_redmarrow * pv_redmarrow_lab / vv_redmarrow +
      fq_skin * pv_skin_lab / vv_skin +
      fq_adipose * pv_adipose_lab / vv_adipose +
      fq_heart * pv_heart_lab / vv_heart +
      fq_rest * pv_rest_lab / vv_rest +
      fq_brain * pv_brain_lab / vv_brain - lambda_phy * pven_lab

    # Peptide bound to serum protein, S1 Eq 13
    d/dt(pprp_lab) <- kpr * pven_lab - lambda_phy * pprp_lab

    # Free peptide, vascular spaces
    d/dt(pv_tumor1_lab) <- ps_tumor1 * (pint_tumor1_lab / vi_tumor1 - pv_tumor1_lab / vv_tumor1) +
      fq_tumor1 * (part_lab / vv_art - pv_tumor1_lab / vv_tumor1) - lambda_phy * pv_tumor1_lab
    d/dt(pv_tumor2_lab) <- ps_tumor2 * (pint_tumor2_lab / vi_tumor2 - pv_tumor2_lab / vv_tumor2) +
      fq_tumor2 * (part_lab / vv_art - pv_tumor2_lab / vv_tumor2) - lambda_phy * pv_tumor2_lab
    d/dt(pv_tumorrest_lab) <- ps_tumorrest * (pint_tumorrest_lab / vi_tumorrest - pv_tumorrest_lab / vv_tumorrest) +
      fq_tumorrest * (part_lab / vv_art - pv_tumorrest_lab / vv_tumorrest) - lambda_phy * pv_tumorrest_lab
    d/dt(pv_salgland_lab) <- ps_salgland * (pint_salgland_lab / vi_salgland - pv_salgland_lab / vv_salgland) +
      fq_salgland * (part_lab / vv_art - pv_salgland_lab / vv_salgland) - lambda_phy * pv_salgland_lab
    # Liver, S1 Eq 4 with the triple inflow of S1 compartment diagrams B2: arterial,
    # splenic and gastrointestinal.
    d/dt(pv_liver_lab) <- ps_liver * (pint_liver_lab / vi_liver - pv_liver_lab / vv_liver) +
      fq_liver * part_lab / vv_art + fq_spleen * pv_spleen_lab / vv_spleen +
      fq_gi * pv_gi_lab / vv_gi -
      (fq_liver + fq_spleen + fq_gi) * pv_liver_lab / vv_liver - lambda_phy * pv_liver_lab
    d/dt(pv_spleen_lab) <- ps_spleen * (pint_spleen_lab / vi_spleen - pv_spleen_lab / vv_spleen) +
      fq_spleen * (part_lab / vv_art - pv_spleen_lab / vv_spleen) - lambda_phy * pv_spleen_lab
    d/dt(pv_gi_lab) <- ps_gi * (pint_gi_lab / vi_gi - pv_gi_lab / vv_gi) +
      fq_gi * (part_lab / vv_art - pv_gi_lab / vv_gi) - lambda_phy * pv_gi_lab
    d/dt(pv_prostate_lab) <- ps_prostate * (pint_prostate_lab / vi_prostate - pv_prostate_lab / vv_prostate) +
      fq_prostate * (part_lab / vv_art - pv_prostate_lab / vv_prostate) - lambda_phy * pv_prostate_lab
    # Kidneys, S1 Eq 6: filtration leaves the vasculature and the
    # intratubular cells return the non-excreted fraction.
    d/dt(pv_kidney_lab) <- - pv_kidney_lab / vv_kidney * (f_fil + fq_kidney) +
      fq_kidney * part_lab / vv_art + f_reab * pintra_kidney_lab / vintra_kidney - lambda_phy * pv_kidney_lab
    d/dt(pv_muscle_lab) <- ps_muscle * (pint_muscle_lab / vi_muscle - pv_muscle_lab / vv_muscle) +
      fq_muscle * (part_lab / vv_art - pv_muscle_lab / vv_muscle) - lambda_phy * pv_muscle_lab
    # Lungs, S1 Eq 5: perfused from the veins at F_TOTAL, drained to the arteries.
    d/dt(pv_lungs_lab) <- ps_lungs * (pint_lungs_lab / vi_lungs - pv_lungs_lab / vv_lungs) +
      fq_total * (pven_lab / vv_ven - pv_lungs_lab / vv_lungs) - lambda_phy * pv_lungs_lab
    d/dt(pv_bone_lab) <- ps_bone * (pint_bone_lab / vi_bone - pv_bone_lab / vv_bone) +
      fq_bone * (part_lab / vv_art - pv_bone_lab / vv_bone) - lambda_phy * pv_bone_lab
    d/dt(pv_redmarrow_lab) <- ps_redmarrow * (pint_redmarrow_lab / vi_redmarrow - pv_redmarrow_lab / vv_redmarrow) +
      fq_redmarrow * (part_lab / vv_art - pv_redmarrow_lab / vv_redmarrow) - lambda_phy * pv_redmarrow_lab
    d/dt(pv_skin_lab) <- ps_skin * (pint_skin_lab / vi_skin - pv_skin_lab / vv_skin) +
      fq_skin * (part_lab / vv_art - pv_skin_lab / vv_skin) - lambda_phy * pv_skin_lab
    d/dt(pv_adipose_lab) <- ps_adipose * (pint_adipose_lab / vi_adipose - pv_adipose_lab / vv_adipose) +
      fq_adipose * (part_lab / vv_art - pv_adipose_lab / vv_adipose) - lambda_phy * pv_adipose_lab
    d/dt(pv_heart_lab) <- ps_heart * (pint_heart_lab / vi_heart - pv_heart_lab / vv_heart) +
      fq_heart * (part_lab / vv_art - pv_heart_lab / vv_heart) - lambda_phy * pv_heart_lab
    d/dt(pv_rest_lab) <- ps_rest * (pint_rest_lab / vi_rest - pv_rest_lab / vv_rest) +
      fq_rest * (part_lab / vv_art - pv_rest_lab / vv_rest) - lambda_phy * pv_rest_lab
    # Brain, S1 Eq 4 with PS = 0: no interstitial exchange.
    d/dt(pv_brain_lab) <- fq_brain * (part_lab / vv_art - pv_brain_lab / vv_brain) - lambda_phy * pv_brain_lab

    # Free peptide, interstitial spaces
    d/dt(pint_tumor1_lab) <- - kon * pint_tumor1_lab * rfree_tumor1 / vi_tumor1 + koff * rp_tumor1_lab +
      ps_tumor1 * (pv_tumor1_lab / vv_tumor1 - pint_tumor1_lab / vi_tumor1) - lambda_phy * pint_tumor1_lab
    d/dt(pint_tumor2_lab) <- - kon * pint_tumor2_lab * rfree_tumor2 / vi_tumor2 + koff * rp_tumor2_lab +
      ps_tumor2 * (pv_tumor2_lab / vv_tumor2 - pint_tumor2_lab / vi_tumor2) - lambda_phy * pint_tumor2_lab
    d/dt(pint_tumorrest_lab) <- - kon * pint_tumorrest_lab * rfree_tumorrest / vi_tumorrest + koff * rp_tumorrest_lab +
      ps_tumorrest * (pv_tumorrest_lab / vv_tumorrest - pint_tumorrest_lab / vi_tumorrest) - lambda_phy * pint_tumorrest_lab
    d/dt(pint_salgland_lab) <- - kon * pint_salgland_lab * rfree_salgland / vi_salgland + koff * rp_salgland_lab +
      ps_salgland * (pv_salgland_lab / vv_salgland - pint_salgland_lab / vi_salgland) - lambda_phy * pint_salgland_lab
    d/dt(pint_liver_lab) <- - kon * pint_liver_lab * rfree_liver / vi_liver + koff * rp_liver_lab +
      ps_liver * (pv_liver_lab / vv_liver - pint_liver_lab / vi_liver) - lambda_phy * pint_liver_lab
    d/dt(pint_spleen_lab) <- - kon * pint_spleen_lab * rfree_spleen / vi_spleen + koff * rp_spleen_lab +
      ps_spleen * (pv_spleen_lab / vv_spleen - pint_spleen_lab / vi_spleen) - lambda_phy * pint_spleen_lab
    d/dt(pint_gi_lab) <- - kon * pint_gi_lab * rfree_gi / vi_gi + koff * rp_gi_lab +
      ps_gi * (pv_gi_lab / vv_gi - pint_gi_lab / vi_gi) - lambda_phy * pint_gi_lab
    d/dt(pint_prostate_lab) <- - kon * pint_prostate_lab * rfree_prostate / vi_prostate + koff * rp_prostate_lab +
      ps_prostate * (pv_prostate_lab / vv_prostate - pint_prostate_lab / vi_prostate) - lambda_phy * pint_prostate_lab
    # Kidneys, S1 Eq 9: fed by glomerular filtration, not by PS.
    d/dt(pint_kidney_lab) <- - kon * pint_kidney_lab * rfree_kidney / vi_kidney +
      koff * rp_kidney_lab +
      f_fil * (pv_kidney_lab / vv_kidney - pint_kidney_lab / vi_kidney) - lambda_phy * pint_kidney_lab
    d/dt(pint_muscle_lab) <- ps_muscle * (pv_muscle_lab / vv_muscle - pint_muscle_lab / vi_muscle) - lambda_phy * pint_muscle_lab
    d/dt(pint_lungs_lab) <- ps_lungs * (pv_lungs_lab / vv_lungs - pint_lungs_lab / vi_lungs) - lambda_phy * pint_lungs_lab
    d/dt(pint_bone_lab) <- ps_bone * (pv_bone_lab / vv_bone - pint_bone_lab / vi_bone) - lambda_phy * pint_bone_lab
    d/dt(pint_redmarrow_lab) <- ps_redmarrow * (pv_redmarrow_lab / vv_redmarrow - pint_redmarrow_lab / vi_redmarrow) - lambda_phy * pint_redmarrow_lab
    d/dt(pint_skin_lab) <- ps_skin * (pv_skin_lab / vv_skin - pint_skin_lab / vi_skin) - lambda_phy * pint_skin_lab
    d/dt(pint_adipose_lab) <- ps_adipose * (pv_adipose_lab / vv_adipose - pint_adipose_lab / vi_adipose) - lambda_phy * pint_adipose_lab
    d/dt(pint_heart_lab) <- ps_heart * (pv_heart_lab / vv_heart - pint_heart_lab / vi_heart) - lambda_phy * pint_heart_lab
    d/dt(pint_rest_lab) <- ps_rest * (pv_rest_lab / vv_rest - pint_rest_lab / vi_rest) - lambda_phy * pint_rest_lab

    # Unspecific uptake into the proximal tubular cells, S1 Eq 12
    d/dt(pintra_kidney_lab) <- f_reab * (pint_kidney_lab / vi_kidney - pintra_kidney_lab / vintra_kidney) - lambda_phy * pintra_kidney_lab

    # Cumulative urinary excretion (bookkeeping state, not in the paper; it
    # collects the excreted fraction f_ex of the filtered load so that mass
    # balance can be checked).
    d/dt(urine_lab) <- f_exc * pint_kidney_lab / vi_kidney - lambda_phy * urine_lab

    # Cumulative peptide lost by release from the internalised pool. S1 Eq 2
    # gives release no destination compartment, so released peptide and 177Lu
    # leave the body. Bookkeeping state, not in the paper; it closes the mass
    # balance exactly.
    d/dt(cleared_lab) <-
      krel_tumor1 * pintern_tumor1_lab + krel_tumor2 * pintern_tumor2_lab +
      krel_tumorrest * pintern_tumorrest_lab + krel_salgland * pintern_salgland_lab +
      krel_liver * pintern_liver_lab + krel_spleen * pintern_spleen_lab +
      krel_gi * pintern_gi_lab + krel_prostate * pintern_prostate_lab +
      krel_kidney * pintern_kidney_lab - lambda_phy * cleared_lab

    # Peptide bound to PSMA on the cell surface, S1 Eq 3
    d/dt(rp_tumor1_lab) <- kon * pint_tumor1_lab * rfree_tumor1 / vi_tumor1 -
      (koff + kint) * rp_tumor1_lab - lambda_phy * rp_tumor1_lab
    d/dt(rp_tumor2_lab) <- kon * pint_tumor2_lab * rfree_tumor2 / vi_tumor2 -
      (koff + kint) * rp_tumor2_lab - lambda_phy * rp_tumor2_lab
    d/dt(rp_tumorrest_lab) <- kon * pint_tumorrest_lab * rfree_tumorrest / vi_tumorrest -
      (koff + kint) * rp_tumorrest_lab - lambda_phy * rp_tumorrest_lab
    d/dt(rp_salgland_lab) <- kon * pint_salgland_lab * rfree_salgland / vi_salgland -
      (koff + kint) * rp_salgland_lab - lambda_phy * rp_salgland_lab
    d/dt(rp_liver_lab) <- kon * pint_liver_lab * rfree_liver / vi_liver -
      (koff + kint) * rp_liver_lab - lambda_phy * rp_liver_lab
    d/dt(rp_spleen_lab) <- kon * pint_spleen_lab * rfree_spleen / vi_spleen -
      (koff + kint) * rp_spleen_lab - lambda_phy * rp_spleen_lab
    d/dt(rp_gi_lab) <- kon * pint_gi_lab * rfree_gi / vi_gi -
      (koff + kint) * rp_gi_lab - lambda_phy * rp_gi_lab
    d/dt(rp_prostate_lab) <- kon * pint_prostate_lab * rfree_prostate / vi_prostate -
      (koff + kint) * rp_prostate_lab - lambda_phy * rp_prostate_lab
    d/dt(rp_kidney_lab) <- kon * pint_kidney_lab * rfree_kidney / vi_kidney -
      (koff + kint) * rp_kidney_lab - lambda_phy * rp_kidney_lab

    # Internalised peptide, S1 Eq 2. Released peptide and free 177Lu are
    # assumed cleared from the body directly, so release is a pure sink.
    d/dt(pintern_tumor1_lab) <- kint * rp_tumor1_lab - krel_tumor1 * pintern_tumor1_lab - lambda_phy * pintern_tumor1_lab
    d/dt(pintern_tumor2_lab) <- kint * rp_tumor2_lab - krel_tumor2 * pintern_tumor2_lab - lambda_phy * pintern_tumor2_lab
    d/dt(pintern_tumorrest_lab) <- kint * rp_tumorrest_lab - krel_tumorrest * pintern_tumorrest_lab - lambda_phy * pintern_tumorrest_lab
    d/dt(pintern_salgland_lab) <- kint * rp_salgland_lab - krel_salgland * pintern_salgland_lab - lambda_phy * pintern_salgland_lab
    d/dt(pintern_liver_lab) <- kint * rp_liver_lab - krel_liver * pintern_liver_lab - lambda_phy * pintern_liver_lab
    d/dt(pintern_spleen_lab) <- kint * rp_spleen_lab - krel_spleen * pintern_spleen_lab - lambda_phy * pintern_spleen_lab
    d/dt(pintern_gi_lab) <- kint * rp_gi_lab - krel_gi * pintern_gi_lab - lambda_phy * pintern_gi_lab
    d/dt(pintern_prostate_lab) <- kint * rp_prostate_lab - krel_prostate * pintern_prostate_lab - lambda_phy * pintern_prostate_lab
    d/dt(pintern_kidney_lab) <- kint * rp_kidney_lab - krel_kidney * pintern_kidney_lab - lambda_phy * pintern_kidney_lab

    # ---------------------------------------------------------------------
    # Unlabelled circulation
    # ---------------------------------------------------------------------

    # Arteries, S1 Eq 8. The printed outflow term multiplies P_i,v where P_ART is
    # meant, and the printed inflow uses F; F_TOTAL is used for both so that the
    # arterial outflow to every organ (tumours included) balances the inflow.
    d/dt(part_unl) <- fq_total * (pv_lungs_unl / vv_lungs - part_unl / vv_art) + lambda_phy * part_lab

    # Veins, S1 Eq 7, plus the outflow to the lungs that the printed equation
    # omits (F_TOTAL, as in S1 Eq 5). Spleen and GI drain into the liver rather than
    # into the veins, so the liver returns F_L + F_S + F_GI (S1 compartment diagrams B2).
    d/dt(pven_unl) <- - kpr * pven_unl - fq_total * pven_unl / vv_ven +
      (fq_liver + fq_spleen + fq_gi) * pv_liver_unl / vv_liver +
      fq_tumor1 * pv_tumor1_unl / vv_tumor1 +
      fq_tumor2 * pv_tumor2_unl / vv_tumor2 +
      fq_tumorrest * pv_tumorrest_unl / vv_tumorrest +
      fq_salgland * pv_salgland_unl / vv_salgland +
      fq_prostate * pv_prostate_unl / vv_prostate +
      fq_kidney * pv_kidney_unl / vv_kidney +
      fq_muscle * pv_muscle_unl / vv_muscle +
      fq_bone * pv_bone_unl / vv_bone +
      fq_redmarrow * pv_redmarrow_unl / vv_redmarrow +
      fq_skin * pv_skin_unl / vv_skin +
      fq_adipose * pv_adipose_unl / vv_adipose +
      fq_heart * pv_heart_unl / vv_heart +
      fq_rest * pv_rest_unl / vv_rest +
      fq_brain * pv_brain_unl / vv_brain + lambda_phy * pven_lab

    # Peptide bound to serum protein, S1 Eq 13
    d/dt(pprp_unl) <- kpr * pven_unl + lambda_phy * pprp_lab

    # Free peptide, vascular spaces
    d/dt(pv_tumor1_unl) <- ps_tumor1 * (pint_tumor1_unl / vi_tumor1 - pv_tumor1_unl / vv_tumor1) +
      fq_tumor1 * (part_unl / vv_art - pv_tumor1_unl / vv_tumor1) + lambda_phy * pv_tumor1_lab
    d/dt(pv_tumor2_unl) <- ps_tumor2 * (pint_tumor2_unl / vi_tumor2 - pv_tumor2_unl / vv_tumor2) +
      fq_tumor2 * (part_unl / vv_art - pv_tumor2_unl / vv_tumor2) + lambda_phy * pv_tumor2_lab
    d/dt(pv_tumorrest_unl) <- ps_tumorrest * (pint_tumorrest_unl / vi_tumorrest - pv_tumorrest_unl / vv_tumorrest) +
      fq_tumorrest * (part_unl / vv_art - pv_tumorrest_unl / vv_tumorrest) + lambda_phy * pv_tumorrest_lab
    d/dt(pv_salgland_unl) <- ps_salgland * (pint_salgland_unl / vi_salgland - pv_salgland_unl / vv_salgland) +
      fq_salgland * (part_unl / vv_art - pv_salgland_unl / vv_salgland) + lambda_phy * pv_salgland_lab
    # Liver, S1 Eq 4 with the triple inflow of S1 compartment diagrams B2: arterial,
    # splenic and gastrointestinal.
    d/dt(pv_liver_unl) <- ps_liver * (pint_liver_unl / vi_liver - pv_liver_unl / vv_liver) +
      fq_liver * part_unl / vv_art + fq_spleen * pv_spleen_unl / vv_spleen +
      fq_gi * pv_gi_unl / vv_gi -
      (fq_liver + fq_spleen + fq_gi) * pv_liver_unl / vv_liver + lambda_phy * pv_liver_lab
    d/dt(pv_spleen_unl) <- ps_spleen * (pint_spleen_unl / vi_spleen - pv_spleen_unl / vv_spleen) +
      fq_spleen * (part_unl / vv_art - pv_spleen_unl / vv_spleen) + lambda_phy * pv_spleen_lab
    d/dt(pv_gi_unl) <- ps_gi * (pint_gi_unl / vi_gi - pv_gi_unl / vv_gi) +
      fq_gi * (part_unl / vv_art - pv_gi_unl / vv_gi) + lambda_phy * pv_gi_lab
    d/dt(pv_prostate_unl) <- ps_prostate * (pint_prostate_unl / vi_prostate - pv_prostate_unl / vv_prostate) +
      fq_prostate * (part_unl / vv_art - pv_prostate_unl / vv_prostate) + lambda_phy * pv_prostate_lab
    # Kidneys, S1 Eq 6: filtration leaves the vasculature and the
    # intratubular cells return the non-excreted fraction.
    d/dt(pv_kidney_unl) <- - pv_kidney_unl / vv_kidney * (f_fil + fq_kidney) +
      fq_kidney * part_unl / vv_art + f_reab * pintra_kidney_unl / vintra_kidney + lambda_phy * pv_kidney_lab
    d/dt(pv_muscle_unl) <- ps_muscle * (pint_muscle_unl / vi_muscle - pv_muscle_unl / vv_muscle) +
      fq_muscle * (part_unl / vv_art - pv_muscle_unl / vv_muscle) + lambda_phy * pv_muscle_lab
    # Lungs, S1 Eq 5: perfused from the veins at F_TOTAL, drained to the arteries.
    d/dt(pv_lungs_unl) <- ps_lungs * (pint_lungs_unl / vi_lungs - pv_lungs_unl / vv_lungs) +
      fq_total * (pven_unl / vv_ven - pv_lungs_unl / vv_lungs) + lambda_phy * pv_lungs_lab
    d/dt(pv_bone_unl) <- ps_bone * (pint_bone_unl / vi_bone - pv_bone_unl / vv_bone) +
      fq_bone * (part_unl / vv_art - pv_bone_unl / vv_bone) + lambda_phy * pv_bone_lab
    d/dt(pv_redmarrow_unl) <- ps_redmarrow * (pint_redmarrow_unl / vi_redmarrow - pv_redmarrow_unl / vv_redmarrow) +
      fq_redmarrow * (part_unl / vv_art - pv_redmarrow_unl / vv_redmarrow) + lambda_phy * pv_redmarrow_lab
    d/dt(pv_skin_unl) <- ps_skin * (pint_skin_unl / vi_skin - pv_skin_unl / vv_skin) +
      fq_skin * (part_unl / vv_art - pv_skin_unl / vv_skin) + lambda_phy * pv_skin_lab
    d/dt(pv_adipose_unl) <- ps_adipose * (pint_adipose_unl / vi_adipose - pv_adipose_unl / vv_adipose) +
      fq_adipose * (part_unl / vv_art - pv_adipose_unl / vv_adipose) + lambda_phy * pv_adipose_lab
    d/dt(pv_heart_unl) <- ps_heart * (pint_heart_unl / vi_heart - pv_heart_unl / vv_heart) +
      fq_heart * (part_unl / vv_art - pv_heart_unl / vv_heart) + lambda_phy * pv_heart_lab
    d/dt(pv_rest_unl) <- ps_rest * (pint_rest_unl / vi_rest - pv_rest_unl / vv_rest) +
      fq_rest * (part_unl / vv_art - pv_rest_unl / vv_rest) + lambda_phy * pv_rest_lab
    # Brain, S1 Eq 4 with PS = 0: no interstitial exchange.
    d/dt(pv_brain_unl) <- fq_brain * (part_unl / vv_art - pv_brain_unl / vv_brain) + lambda_phy * pv_brain_lab

    # Free peptide, interstitial spaces
    d/dt(pint_tumor1_unl) <- - kon * pint_tumor1_unl * rfree_tumor1 / vi_tumor1 + koff * rp_tumor1_unl +
      ps_tumor1 * (pv_tumor1_unl / vv_tumor1 - pint_tumor1_unl / vi_tumor1) + lambda_phy * pint_tumor1_lab
    d/dt(pint_tumor2_unl) <- - kon * pint_tumor2_unl * rfree_tumor2 / vi_tumor2 + koff * rp_tumor2_unl +
      ps_tumor2 * (pv_tumor2_unl / vv_tumor2 - pint_tumor2_unl / vi_tumor2) + lambda_phy * pint_tumor2_lab
    d/dt(pint_tumorrest_unl) <- - kon * pint_tumorrest_unl * rfree_tumorrest / vi_tumorrest + koff * rp_tumorrest_unl +
      ps_tumorrest * (pv_tumorrest_unl / vv_tumorrest - pint_tumorrest_unl / vi_tumorrest) + lambda_phy * pint_tumorrest_lab
    d/dt(pint_salgland_unl) <- - kon * pint_salgland_unl * rfree_salgland / vi_salgland + koff * rp_salgland_unl +
      ps_salgland * (pv_salgland_unl / vv_salgland - pint_salgland_unl / vi_salgland) + lambda_phy * pint_salgland_lab
    d/dt(pint_liver_unl) <- - kon * pint_liver_unl * rfree_liver / vi_liver + koff * rp_liver_unl +
      ps_liver * (pv_liver_unl / vv_liver - pint_liver_unl / vi_liver) + lambda_phy * pint_liver_lab
    d/dt(pint_spleen_unl) <- - kon * pint_spleen_unl * rfree_spleen / vi_spleen + koff * rp_spleen_unl +
      ps_spleen * (pv_spleen_unl / vv_spleen - pint_spleen_unl / vi_spleen) + lambda_phy * pint_spleen_lab
    d/dt(pint_gi_unl) <- - kon * pint_gi_unl * rfree_gi / vi_gi + koff * rp_gi_unl +
      ps_gi * (pv_gi_unl / vv_gi - pint_gi_unl / vi_gi) + lambda_phy * pint_gi_lab
    d/dt(pint_prostate_unl) <- - kon * pint_prostate_unl * rfree_prostate / vi_prostate + koff * rp_prostate_unl +
      ps_prostate * (pv_prostate_unl / vv_prostate - pint_prostate_unl / vi_prostate) + lambda_phy * pint_prostate_lab
    # Kidneys, S1 Eq 9: fed by glomerular filtration, not by PS.
    d/dt(pint_kidney_unl) <- - kon * pint_kidney_unl * rfree_kidney / vi_kidney +
      koff * rp_kidney_unl +
      f_fil * (pv_kidney_unl / vv_kidney - pint_kidney_unl / vi_kidney) + lambda_phy * pint_kidney_lab
    d/dt(pint_muscle_unl) <- ps_muscle * (pv_muscle_unl / vv_muscle - pint_muscle_unl / vi_muscle) + lambda_phy * pint_muscle_lab
    d/dt(pint_lungs_unl) <- ps_lungs * (pv_lungs_unl / vv_lungs - pint_lungs_unl / vi_lungs) + lambda_phy * pint_lungs_lab
    d/dt(pint_bone_unl) <- ps_bone * (pv_bone_unl / vv_bone - pint_bone_unl / vi_bone) + lambda_phy * pint_bone_lab
    d/dt(pint_redmarrow_unl) <- ps_redmarrow * (pv_redmarrow_unl / vv_redmarrow - pint_redmarrow_unl / vi_redmarrow) + lambda_phy * pint_redmarrow_lab
    d/dt(pint_skin_unl) <- ps_skin * (pv_skin_unl / vv_skin - pint_skin_unl / vi_skin) + lambda_phy * pint_skin_lab
    d/dt(pint_adipose_unl) <- ps_adipose * (pv_adipose_unl / vv_adipose - pint_adipose_unl / vi_adipose) + lambda_phy * pint_adipose_lab
    d/dt(pint_heart_unl) <- ps_heart * (pv_heart_unl / vv_heart - pint_heart_unl / vi_heart) + lambda_phy * pint_heart_lab
    d/dt(pint_rest_unl) <- ps_rest * (pv_rest_unl / vv_rest - pint_rest_unl / vi_rest) + lambda_phy * pint_rest_lab

    # Unspecific uptake into the proximal tubular cells, S1 Eq 12
    d/dt(pintra_kidney_unl) <- f_reab * (pint_kidney_unl / vi_kidney - pintra_kidney_unl / vintra_kidney) + lambda_phy * pintra_kidney_lab

    # Cumulative urinary excretion (bookkeeping state, not in the paper; it
    # collects the excreted fraction f_ex of the filtered load so that mass
    # balance can be checked).
    d/dt(urine_unl) <- f_exc * pint_kidney_unl / vi_kidney + lambda_phy * urine_lab

    # Cumulative peptide lost by release from the internalised pool. S1 Eq 2
    # gives release no destination compartment, so released peptide and 177Lu
    # leave the body. Bookkeeping state, not in the paper; it closes the mass
    # balance exactly.
    d/dt(cleared_unl) <-
      krel_tumor1 * pintern_tumor1_unl + krel_tumor2 * pintern_tumor2_unl +
      krel_tumorrest * pintern_tumorrest_unl + krel_salgland * pintern_salgland_unl +
      krel_liver * pintern_liver_unl + krel_spleen * pintern_spleen_unl +
      krel_gi * pintern_gi_unl + krel_prostate * pintern_prostate_unl +
      krel_kidney * pintern_kidney_unl + lambda_phy * cleared_lab

    # Peptide bound to PSMA on the cell surface, S1 Eq 3
    d/dt(rp_tumor1_unl) <- kon * pint_tumor1_unl * rfree_tumor1 / vi_tumor1 -
      (koff + kint) * rp_tumor1_unl + lambda_phy * rp_tumor1_lab
    d/dt(rp_tumor2_unl) <- kon * pint_tumor2_unl * rfree_tumor2 / vi_tumor2 -
      (koff + kint) * rp_tumor2_unl + lambda_phy * rp_tumor2_lab
    d/dt(rp_tumorrest_unl) <- kon * pint_tumorrest_unl * rfree_tumorrest / vi_tumorrest -
      (koff + kint) * rp_tumorrest_unl + lambda_phy * rp_tumorrest_lab
    d/dt(rp_salgland_unl) <- kon * pint_salgland_unl * rfree_salgland / vi_salgland -
      (koff + kint) * rp_salgland_unl + lambda_phy * rp_salgland_lab
    d/dt(rp_liver_unl) <- kon * pint_liver_unl * rfree_liver / vi_liver -
      (koff + kint) * rp_liver_unl + lambda_phy * rp_liver_lab
    d/dt(rp_spleen_unl) <- kon * pint_spleen_unl * rfree_spleen / vi_spleen -
      (koff + kint) * rp_spleen_unl + lambda_phy * rp_spleen_lab
    d/dt(rp_gi_unl) <- kon * pint_gi_unl * rfree_gi / vi_gi -
      (koff + kint) * rp_gi_unl + lambda_phy * rp_gi_lab
    d/dt(rp_prostate_unl) <- kon * pint_prostate_unl * rfree_prostate / vi_prostate -
      (koff + kint) * rp_prostate_unl + lambda_phy * rp_prostate_lab
    d/dt(rp_kidney_unl) <- kon * pint_kidney_unl * rfree_kidney / vi_kidney -
      (koff + kint) * rp_kidney_unl + lambda_phy * rp_kidney_lab

    # Internalised peptide, S1 Eq 2. Released peptide and free 177Lu are
    # assumed cleared from the body directly, so release is a pure sink.
    d/dt(pintern_tumor1_unl) <- kint * rp_tumor1_unl - krel_tumor1 * pintern_tumor1_unl + lambda_phy * pintern_tumor1_lab
    d/dt(pintern_tumor2_unl) <- kint * rp_tumor2_unl - krel_tumor2 * pintern_tumor2_unl + lambda_phy * pintern_tumor2_lab
    d/dt(pintern_tumorrest_unl) <- kint * rp_tumorrest_unl - krel_tumorrest * pintern_tumorrest_unl + lambda_phy * pintern_tumorrest_lab
    d/dt(pintern_salgland_unl) <- kint * rp_salgland_unl - krel_salgland * pintern_salgland_unl + lambda_phy * pintern_salgland_lab
    d/dt(pintern_liver_unl) <- kint * rp_liver_unl - krel_liver * pintern_liver_unl + lambda_phy * pintern_liver_lab
    d/dt(pintern_spleen_unl) <- kint * rp_spleen_unl - krel_spleen * pintern_spleen_unl + lambda_phy * pintern_spleen_lab
    d/dt(pintern_gi_unl) <- kint * rp_gi_unl - krel_gi * pintern_gi_unl + lambda_phy * pintern_gi_lab
    d/dt(pintern_prostate_unl) <- kint * rp_prostate_unl - krel_prostate * pintern_prostate_unl + lambda_phy * pintern_prostate_lab
    d/dt(pintern_kidney_unl) <- kint * rp_kidney_unl - krel_kidney * pintern_kidney_unl + lambda_phy * pintern_kidney_lab

    # =====================================================================
    # 3. Observed quantities. S1 'Assigning data to compartments' gives the
    #    labelled-peptide amount seen by the gamma camera in each region of
    #    interest: the organ self-content plus that organ's share of the
    #    serum-protein-bound pool, and for the tumour ROIs a fitted fraction
    #    of the labelled peptide in muscle as background correction.
    # =====================================================================
    vprp_denom <- vserum + vv_tumor1 + vv_tumor2 + vv_tumorrest
    prp_muscle <- pv_muscle_lab + pint_muscle_lab + pprp_lab * vv_muscle / vprp_denom

    # Data_K, including the unspecific tubular-cell pool P*_K,intra
    Pkidney <- pv_kidney_lab + pint_kidney_lab + rp_kidney_lab + pintern_kidney_lab +
      pintra_kidney_lab + pprp_lab * vv_kidney / vprp_denom
    # Data_SAL
    Psalgland <- pv_salgland_lab + pint_salgland_lab + rp_salgland_lab + pintern_salgland_lab +
      pprp_lab * vv_salgland / vprp_denom
    # Data_TU,1 (x_tu,1 = bgcorr_tumor1)
    Ptumor1 <- pv_tumor1_lab + pint_tumor1_lab + rp_tumor1_lab + pintern_tumor1_lab +
      pprp_lab * vv_tumor1 / vprp_denom + bgcorr_tumor1 * prp_muscle
    # Data_TU,2 (x_tu,2 = bgcorr_tumor2)
    Ptumor2 <- pv_tumor2_lab + pint_tumor2_lab + rp_tumor2_lab + pintern_tumor2_lab +
      pprp_lab * vv_tumor2 / vprp_denom + bgcorr_tumor2 * prp_muscle
    # The lumped remaining lesions follow the same form, with no ROI background
    # correction because they are not delineated on the gamma camera.
    Ptumorrest <- pv_tumorrest_lab + pint_tumorrest_lab + rp_tumorrest_lab +
      pintern_tumorrest_lab + pprp_lab * vv_tumorrest / vprp_denom
    Ptumortotal <- Ptumor1 + Ptumor2 + Ptumorrest

    # Total body: all labelled ('hot') compartments; and the venous serum
    # concentration.
    Pbody <- pven_lab + part_lab + pprp_lab + pintra_kidney_lab +
      pv_tumor1_lab + pv_tumor2_lab + pv_tumorrest_lab +
      pv_salgland_lab + pv_liver_lab + pv_spleen_lab +
      pv_gi_lab + pv_prostate_lab + pv_kidney_lab +
      pv_muscle_lab + pv_lungs_lab + pv_bone_lab +
      pv_redmarrow_lab + pv_skin_lab + pv_adipose_lab +
      pv_heart_lab + pv_rest_lab + pv_brain_lab +
      pint_tumor1_lab + pint_tumor2_lab + pint_tumorrest_lab +
      pint_salgland_lab + pint_liver_lab + pint_spleen_lab +
      pint_gi_lab + pint_prostate_lab + pint_kidney_lab +
      pint_muscle_lab + pint_lungs_lab + pint_bone_lab +
      pint_redmarrow_lab + pint_skin_lab + pint_adipose_lab +
      pint_heart_lab + pint_rest_lab + rp_tumor1_lab +
      rp_tumor2_lab + rp_tumorrest_lab + rp_salgland_lab +
      rp_liver_lab + rp_spleen_lab + rp_gi_lab +
      rp_prostate_lab + rp_kidney_lab + pintern_tumor1_lab +
      pintern_tumor2_lab + pintern_tumorrest_lab + pintern_salgland_lab +
      pintern_liver_lab + pintern_spleen_lab + pintern_gi_lab +
      pintern_prostate_lab + pintern_kidney_lab
    Cc <- pven_lab / vv_ven
    # Fraction of PSMA sites occupied (labelled plus unlabelled) in the two lesions.
    occ_tumor1 <- (rp_tumor1_lab + rp_tumor1_unl) / r0_tumor1
    occ_tumor2 <- (rp_tumor2_lab + rp_tumor2_unl) / r0_tumor2
  })
}
