Golzaryan_2025_lu177psmaIT_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 21 compartments, SimBiology). Personalized",
    "physiologically based pharmacokinetic model for the PSMA-targeted",
    "radioligand [177Lu]Lu-PSMA I&T in men with metastatic",
    "castration-resistant prostate cancer, built to compare metronomic",
    "(multiple lower-dose) against conventional single-dose",
    "radiopharmaceutical therapy. Two parallel circulations are carried:",
    "one for radiolabelled peptide and one for unlabelled peptide, coupled",
    "only by the physical decay of 177Lu, which converts labelled species",
    "into their unlabelled counterparts in every compartment. The two",
    "species compete for the same finite pool of PSMA binding sites. Nine",
    "PSMA-positive organs (two explicit tumour lesions, a lumped",
    "tumour-remainder lesion, salivary glands, liver, spleen, GI tract",
    "plus pancreas, prostate and kidneys) carry vascular, interstitial,",
    "receptor-bound and internalised sub-compartments; the kidneys add a",
    "fifth intratubular sub-compartment together with glomerular",
    "filtration and tubular excretion. Eight PSMA-negative organs carry",
    "vascular and interstitial sub-compartments, the brain carries a",
    "vascular sub-compartment only (PS = 0), and arteries, veins and a",
    "serum-protein-bound pool complete the 21. Transcapillary exchange is",
    "by permeability-surface-area product, binding is second order onto a",
    "conserved receptor pool, and internalised peptide is released",
    "irreversibly. Ten parameters (four binding-site densities, three",
    "release rates, two tumour-ROI background corrections and salivary",
    "gland perfusion) were fitted per patient to gamma-camera",
    "time-activity curves; the shipped values are patient 2, the paper's",
    "worked example. Amounts are nmol and time is minutes.",
    sep = " "
  )
  reference <- paste(
    "Golzaryan A, Soltani M, Moradi Kashkooli F, Saboury B, Rahmim A.",
    "Personalized metronomic radiopharmaceutical therapy through injection",
    "profile optimization via physiologically based pharmacokinetic (PBPK)",
    "modeling. Sci Rep. 2025;15(1):4046. doi:10.1038/s41598-025-86159-9.",
    "Patient measurements, the gamma-camera time-activity curves the model",
    "was fitted to, and the k_off model selection are from the upstream",
    "primary Kletting P, Thieme A, Eberhardt N, et al. Modeling and",
    "Predicting Tumor Response in Radioligand Therapy. PLoS One.",
    "2016;11(9):e0162303. doi:10.1371/journal.pone.0162303.",
    sep = " "
  )
  vignette <- "Golzaryan_2025_lu177psmaIT_pbpk"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/L")

  paper_specific_compartment_pattern <- "^(pv|pint|rp|pintern|pintra|part|pven|pprp|urine|cleared)_"

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales ten ICRP-23 Reference Man organ volumes, each written in Table S2",
        "as coefficient * BW/71. NOT published for these five patients: absent",
        "from the paper, its supplement, the upstream Kletting 2016 primary and",
        "that paper's own S1 File. The shipped demonstration value is 71 kg, the",
        "reference individual the paper's own BW/71 normalisation encodes (its",
        "cited ref 22 is ICRP Publication 23); assumed, not paper-derived."
      ),
      source_name        = "BW"
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "% (volume fraction x 100)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as the fraction H = HCT/100. Sets total serum volume V_P =",
        "2.8*(1-H)*BSA (Table S2) and hence total serum flow and every vascular",
        "sub-volume, and also v_TU,v = 0.05*(1-H), f_PRO = 0.18*(1-H) and V_SAL,v",
        "= 0.03*(1-H)*V_SAL,total. Declared measured in Table S2 but not",
        "published for these five patients; the shipped 45 % is assumed, not",
        "paper-derived."
      ),
      source_name        = "H"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured per patient, Table S1 (P1-P5: 2.0, 1.9, 1.8, 2.1, 2.0 m^2).",
        "Used directly rather than via the Du Bois formula BSA =",
        "0.007184*BH^0.725*BW^0.425 that Table S2 also lists, because BSA is the",
        "quantity the paper actually tabulates."
      ),
      source_name        = "BSA"
    ),
    TER_MAG3 = list(
      description        = "Tubular extraction rate measured by 99mTc-MAG3 renography",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured per patient, Table S1 (P1-P5: 198, 201, 136, 252, 176 mL/min).",
        "Table S2 derives the 51Cr-EDTA glomerular filtration rate from it as GFR",
        "= TER/3 * 20/15, which is then scaled to peptide molecular size by phi =",
        "0.66."
      ),
      source_name        = "TER"
    ),
    ORGVOL_KIDNEY = list(
      description        = "Measured total kidney volume",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured per patient, Table S1 (P1-P5: 321, 311, 394, 268, 296 mL).",
        "Split into vascular (5.5 %), interstitial (15 %) and intratubular (two",
        "thirds of the remainder) sub-volumes."
      ),
      source_name        = "Kidney (measured volume)"
    ),
    ORGVOL_SALGLAND = list(
      description        = "Measured total salivary gland volume (left plus right parotid)",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured per patient, Table S1 (P1-P5: 54, 21, 17, 52, 29 mL)."
      ),
      source_name        = "Salivary Glands (measured volume)"
    ),
    ORGVOL_TUMOR1 = list(
      description        = "Measured volume of explicitly modelled tumour lesion 1",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured per patient, Table S1 (P1-P5: 0.5, 1, 2, 4, 1.5 mL). Patient 4",
        "is 4 mL, not the 5 mL printed in Table S1; Table S4 tabulates S-values",
        "for exactly the 13 distinct volumes the model evaluates and has a row",
        "for 4 mL and none for 5 mL, and the upstream Kletting 2016 Table 1 that",
        "Table S1 reproduces also gives 4 mL."
      ),
      source_name        = "Tumor 1 (measured volume)"
    ),
    ORGVOL_TUMOR2 = list(
      description        = "Measured volume of explicitly modelled tumour lesion 2",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Measured per patient, Table S1 (P1-P5: 1, 34, 13, 3, 1 mL)."
      ),
      source_name        = "Tumor 2 (measured volume)"
    ),
    ORGVOL_TUMORREST = list(
      description        = "Assumed volume of the lumped remaining-tumour-lesion compartment",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Table S2 footnote c assumes either a 10 mL or a 50 mL additional tumour",
        "volume at 266 nmol/L binding-site density, giving the tabulated total",
        "binding sites R_TU,Rest,0 of 2.6 and 13 nmol respectively. Back-solving",
        "those against P1-P5 (13, 2.6, 13, 2.6, 13 nmol) gives 50, 10, 50, 10, 50",
        "mL."
      ),
      source_name        = "R_TU,Rest,0 / [R_TU,Rest,0]"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported per patient in Table S1 (54-78 years) as a cohort descriptor.",
        "Not used in any model equation; declared for provenance only."
      ),
      source_name        = "Age"
    )
  )

  compartmentData <- list(
    pv_tumor1_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumor2_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumorrest_lab           = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_salgland_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_liver_lab               = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_spleen_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_gi_lab                  = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_prostate_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_kidney_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_muscle_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_lungs_lab               = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_bone_lab                = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_redmarrow_lab           = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_skin_lab                = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_adipose_lab             = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_heart_lab               = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_rest_lab                = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pv_brain_lab               = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pint_tumor1_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumor2_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumorrest_lab         = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_salgland_lab          = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_liver_lab             = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_spleen_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_gi_lab                = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_prostate_lab          = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_kidney_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_muscle_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_lungs_lab             = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_bone_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_redmarrow_lab         = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_skin_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_adipose_lab           = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_heart_lab             = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_rest_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_tumor1_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor1_lab         = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumor2_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor2_lab         = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumorrest_lab           = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumorrest_lab      = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_salgland_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_salgland_lab       = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_liver_lab               = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_liver_lab          = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_spleen_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_spleen_lab         = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_gi_lab                  = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_gi_lab             = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_prostate_lab            = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_prostate_lab       = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_kidney_lab              = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_kidney_lab         = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    pintra_kidney_lab          = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "tissue", verified = TRUE),
    part_lab                   = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pven_lab                   = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    pprp_lab                   = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "serum", verified = TRUE),
    urine_lab                  = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "urine", verified = TRUE),
    cleared_lab                = list(analyte = "[177Lu]Lu-PSMA I&T", units = "nmol", specimen = "not applicable", verified = TRUE),
    pv_tumor1_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumor2_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_tumorrest_unl           = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_salgland_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_liver_unl               = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_spleen_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_gi_unl                  = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_prostate_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_kidney_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_muscle_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_lungs_unl               = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_bone_unl                = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_redmarrow_unl           = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_skin_unl                = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_adipose_unl             = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_heart_unl               = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_rest_unl                = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_brain_unl               = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pint_tumor1_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumor2_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_tumorrest_unl         = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_salgland_unl          = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_liver_unl             = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_spleen_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_gi_unl                = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_prostate_unl          = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_kidney_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_muscle_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_lungs_unl             = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_bone_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_redmarrow_unl         = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_skin_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_adipose_unl           = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_heart_unl             = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_rest_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_tumor1_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor1_unl         = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumor2_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumor2_unl         = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_tumorrest_unl           = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_tumorrest_unl      = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_salgland_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_salgland_unl       = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_liver_unl               = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_liver_unl          = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_spleen_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_spleen_unl         = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_gi_unl                  = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_gi_unl             = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_prostate_unl            = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_prostate_unl       = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_kidney_unl              = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_kidney_unl         = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintra_kidney_unl          = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    part_unl                   = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pven_unl                   = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pprp_unl                   = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    urine_unl                  = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "urine", verified = TRUE),
    cleared_unl                = list(analyte = "PSMA I&T (unlabelled)", units = "nmol", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 5,
    n_studies      = 1,
    age_range      = "53-78 years",
    age_median     = "69 years",
    sex_female_pct = 0,
    disease_state  = "metastatic castration-resistant prostate cancer (mCRPC)",
    dose_range     = paste(
      "clinically administered cycle: 74-302 nmol total peptide with 5.4-6.0 GBq",
      "of [177Lu]Lu-PSMA I&T as a 10 min intravenous infusion; the in-silico",
      "metronomic exploration spans 1-22 GBq, 2^5 to 2^10 nmol, 2-6 injections",
      "and 12-36 h between injections",
      sep = " "
    ),
    notes          = paste(
      "Baseline measurements are Table S1 of the supplement (age, body surface",
      "area, 99mTc-MAG3 tubular extraction rate, salivary gland / kidney / tumour",
      "1 / tumour 2 volumes, injected peptide amount and activity), reproduced",
      "from Table 1 of the upstream Kletting 2016 PLoS One primary. Body weight",
      "and hematocrit are declared measured in Table S2 but were never published",
      "for these patients; see covariateData notes. The individual parameters",
      "were fitted to gamma-camera time-activity curves acquired at 0.5 h, 2 h",
      "and 1-5 days post-infusion during one therapy cycle (Tables S5-S9).",
      sep = " "
    )
  )

  ini({
    # ---- Binding and turnover (Table S2) ---------------------------------
    kon        <- fixed(0.046)    ; label("Association rate onto PSMA (L/nmol/min)")           # Table S2, ref (4)
    kd         <- fixed(1)        ; label("PSMA dissociation constant (nmol/L)")               # Table S2 lists 1 or 8; Kletting 2016 selected the K_D = 1 model (Akaike weight > 0.9 in every patient)
    kint       <- fixed(0.001)    ; label("Internalisation rate, all organs (1/min)")          # Table S2, ref (12)
    kpr        <- fixed(4.7e-4)   ; label("Binding rate of peptide to serum protein (1/min)")  # Table S2, ref (12)
    lambda_phy <- fixed(7.15e-5)  ; label("Physical decay rate of 177Lu (1/min)")              # Table S2

    # ---- Fitted release rates -- patient 2 (Table S6) ---------------------
    krel_kidney   <- fixed(0.00029) ; label("Release rate, kidneys and all normal tissue (1/min)")  # Table S6; Table S2 sets lambda_L,release = lambda_S,release = lambda_NT,release = lambda_K,release
    krel_salgland <- fixed(0.00042) ; label("Release rate, salivary glands (1/min)")                # Table S6
    krel_tumor    <- fixed(0.00015) ; label("Release rate, all tumour compartments (1/min)")        # Table S6; Table S2 sets lambda_TU,Rest,release to the mean of the two lesion rates, which is this single fitted value

    # ---- Binding-site densities (nmol/L) ----------------------------------
    rdens_kidney    <- fixed(14)    ; label("PSMA binding-site density, kidneys (nmol/L)")        # Table S6 (fitted)
    rdens_salgland  <- fixed(38)    ; label("PSMA binding-site density, salivary glands (nmol/L)")# Table S6 (fitted)
    rdens_tumor1    <- fixed(57)    ; label("PSMA binding-site density, tumour lesion 1 (nmol/L)")# Table S6 (fitted)
    rdens_tumor2    <- fixed(19)    ; label("PSMA binding-site density, tumour lesion 2 (nmol/L)")# Table S6 (fitted)
    rdens_tumorrest <- fixed(266)   ; label("PSMA binding-site density, tumour remainder (nmol/L)")# Table S2, ref (11)
    rdens_prostate  <- fixed(26.6)  ; label("PSMA binding-site density, prostate (nmol/L)")       # Table S2: [R_PRO,0] = [R_TU,Rest,0] * 0.1, ref (25)
    rdens_liver     <- fixed(1.33)  ; label("PSMA binding-site density, liver (nmol/L)")          # Table S2: [R_L,0] = [R_PRO,0] * 0.05, ref (19)
    rdens_spleen    <- fixed(0.532) ; label("PSMA binding-site density, spleen (nmol/L)")         # Table S2: [R_S,0] = [R_PRO,0] * 0.02, ref (19)
    rdens_gi        <- fixed(1.596) ; label("PSMA binding-site density, GI plus pancreas (nmol/L)")# Table S2: [R_GI,0] = [R_PRO,0] * 0.06, ref (19)

    # ---- Perfusion and permeability densities (mL/min/g) ------------------
    fperf_tumor     <- fixed(0.5)   ; label("Serum flow density, tumour (mL/min/g)")              # Table S2, refs (9, 10)
    fperf_prostate  <- fixed(0.18)  ; label("Serum flow density, prostate, before the (1-H) factor (mL/min/g)") # Table S2, ref (8)
    fperf_salgland  <- fixed(0.074) ; label("Serum flow density, salivary glands (mL/min/g)")     # Table S6 (fitted)
    kps_tumor       <- fixed(0.6)   ; label("Permeability-surface density, tumour (mL/min/g)")    # Table S2, ref (9)
    kps_muscle      <- fixed(0.02)  ; label("Permeability-surface density, muscle-like tissue (mL/min/g)") # Table S2, ref (18); GI, skin, adipose, heart, bone and rest are assumed equal to muscle
    kps_prostate    <- fixed(0.1)   ; label("Permeability-surface density, prostate (mL/min/g)")  # Table S2, ref (8)
    kps_liver       <- fixed(2)     ; label("Permeability-surface density, liver-like tissue (mL/min/g)")  # Table S2: k_L = k_MUS * 100, ref (18); spleen, lungs, salivary glands and red marrow share it

    # ---- Renal handling ---------------------------------------------------
    phi_gfr <- fixed(0.66) ; label("Ratio of sieving coefficients, PSMA I&T to 51Cr-EDTA (unitless)") # Table S2, ref (15)
    f_ex    <- fixed(0.96) ; label("Tubular excretion fraction of the filtered load (unitless)")      # Table S2, ref (17)

    # ---- Tumour ROI background correction -- patient 2 (Table S6) ---------
    bgcorr_tumor1 <- fixed(0.0064) ; label("Muscle background fraction in the tumour 1 ROI, paper symbol c1 (unitless)") # Table S6 (fitted)
    bgcorr_tumor2 <- fixed(0.0013) ; label("Muscle background fraction in the tumour 2 ROI, paper symbol c2 (unitless)") # Table S6 (fitted)
  })

  model({
    # =====================================================================
    # 1. Anatomy and physiology. All volumes in L, all flows in L/min, all
    #    amounts in nmol, time in min. Because the paper reports density
    #    parameters in mL/min/g and assumes rho = 1 g/mL (Table S2), a flow
    #    or PS product in L/min is numerically the density times the organ
    #    volume in L, with no conversion factor.
    # =====================================================================
    hct_frac <- HCT / 100                      # Table S2 carries H as a fraction
    vserum   <- 2.8 * (1 - hct_frac) * BSA     # V_P, Table S2 (7)
    fserum   <- vserum * 1.23                  # F, Table S2 (6) and footnote b

    # ---- Total organ volumes (Table S2) ----------------------------------
    vtot_kidney     <- ORGVOL_KIDNEY / 1000    # measured, Table S1
    vtot_salgland   <- ORGVOL_SALGLAND / 1000  # measured, Table S1
    vtot_tumor1     <- ORGVOL_TUMOR1 / 1000    # measured, Table S1
    vtot_tumor2     <- ORGVOL_TUMOR2 / 1000    # measured, Table S1
    vtot_tumorrest  <- ORGVOL_TUMORREST / 1000 # Table S2 footnote c
    vtot_liver      <- 1.91                    # Table S2, ref (13)
    vtot_spleen     <- 0.183                   # Table S2
    vtot_prostate   <- 0.016  * WT / 71        # Table S2, ref (22) ICRP-23
    vtot_lungs      <- 1      * WT / 71        # Table S2, ref (22)
    vtot_muscle     <- 30.078 * WT / 71        # Table S2, ref (22)
    vtot_gi         <- (0.385 + 0.548 + 0.104 + 0.15) * WT / 71  # GI plus pancreas, Table S2, ref (22)
    vtot_skin       <- 3.408  * WT / 71        # Table S2, ref (22)
    vtot_adipose    <- 13.465 * WT / 71        # Table S2
    vtot_redmarrow  <- 1.1    * WT / 71        # Table S2, ref (22)
    vtot_bone       <- 10.165 * WT / 71 - vtot_redmarrow  # Table S2, ref (22)
    vtot_heart      <- 0.341  * WT / 71        # Table S2, ref (22)
    vtot_brain      <- 1.45   * WT / 71        # Table S2, ref (22)
    # V_REST,total = V_BW - sum of all organ volumes except the tumours (Table S2;
    # the printed formula is an undecodable embedded image, the wording is not).
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

    # ---- Vascular (serum) sub-volumes (Table S2) --------------------------
    vv_lungs     <- 0.105 * vserum            # ref (6)
    vv_muscle    <- 0.14  * vserum            # ref (6)
    vv_gi        <- 0.076 * vserum            # ref (6)
    vv_skin      <- 0.03  * vserum            # ref (6)
    vv_adipose   <- 0.05  * vserum            # ref (6)
    vv_redmarrow <- 0.04  * vserum            # ref (6)
    vv_bone      <- 0.07  * vserum - vv_redmarrow  # ref (6)
    vv_heart     <- 0.01  * vserum            # ref (6)
    vv_brain     <- 0.012 * vserum            # ref (6)
    vv_art       <- 0.105 * vserum            # 0.06 arterial plus half the cardiac serum, ref (6)
    vv_ven       <- 0.225 * vserum            # 0.18 venous plus half the cardiac serum, ref (6)
    vv_liver     <- 0.085 * vtot_liver        # ref (14)
    vv_spleen    <- 0.12  * vtot_spleen       # ref (14)
    vv_kidney    <- 0.055 * vtot_kidney       # ref (14)
    vv_prostate  <- 0.004 * (1 - hct_frac) * vtot_prostate   # ref (8)
    vv_salgland  <- 0.03  * (1 - hct_frac) * vtot_salgland   # ref (23)
    vv_tumor1    <- 0.05  * (1 - hct_frac) * vtot_tumor1     # v_TU,v, ref (9)
    vv_tumor2    <- 0.05  * (1 - hct_frac) * vtot_tumor2     # v_TU,v, ref (9)
    vv_tumorrest <- 0.05  * (1 - hct_frac) * vtot_tumorrest  # v_TU,v, ref (9)
    # V_REST,v is the serum left over once every other serum pool is assigned.
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

    # ---- Interstitial sub-volumes (Table S2) ------------------------------
    vi_liver     <- 0.2  * vtot_liver         # ref (14)
    vi_spleen    <- 0.2  * vtot_spleen        # ref (14)
    vi_kidney    <- 0.15 * vtot_kidney        # ref (14)
    vi_prostate  <- 0.25 * vtot_prostate      # ref (8)
    vi_salgland  <- 0.23 * vtot_salgland      # ref (23)
    vi_tumor1    <- 0.38 * vtot_tumor1        # v_TU,int, ref (8)
    vi_tumor2    <- 0.38 * vtot_tumor2        # v_TU,int, ref (8)
    vi_tumorrest <- 0.38 * vtot_tumorrest     # v_TU,int, ref (8)
    vi_lungs      <- vv_lungs      * 5.5         # alpha_LUNGS, ref (14)
    vi_muscle     <- vv_muscle     * 5.9         # alpha_MUSCLE, ref (14)
    vi_gi         <- vv_gi         * 8.8         # alpha_GI, ref (14)
    vi_skin       <- vv_skin       * 8.9         # alpha_SKIN, ref (14)
    vi_adipose    <- vv_adipose    * 15.5        # alpha_ADIPOSE, ref (14)
    vi_redmarrow  <- vv_redmarrow  * 3.7         # alpha_REDMARROW, ref (14)
    vi_bone       <- vv_bone       * 8.4         # alpha_BONE, ref (14)
    vi_heart      <- vv_heart      * 3.7         # alpha_HEART, ref (14)
    vi_rest       <- vv_rest       * 4.1         # alpha_REST, ref (14)
    # Intracellular kidney: two thirds of what is neither vascular nor interstitial
    vintra_kidney <- (vtot_kidney - vi_kidney - vv_kidney) * 2 / 3  # Table S2 footnote d

    # ---- Serum flows (Table S2) -------------------------------------------
    fq_liver      <- 0.065  * fserum        # ref (6)
    fq_spleen     <- 0.03   * fserum        # ref (6)
    fq_kidney     <- 0.19   * fserum        # ref (6)
    fq_muscle     <- 0.17   * fserum        # ref (6)
    fq_gi         <- 0.16   * fserum        # ref (6)
    fq_skin       <- 0.05   * fserum        # ref (6)
    fq_adipose    <- 0.05   * fserum        # ref (6)
    fq_redmarrow  <- 0.03   * fserum        # ref (6)
    fq_bone       <- 0.05   * fserum        # ref (6)
    fq_heart      <- 0.04   * fserum        # ref (6)
    fq_brain      <- 0.12   * fserum        # ref (6)
    fq_prostate  <- fperf_prostate * (1 - hct_frac) * vtot_prostate  # f_PRO * V_PRO,total, ref (8)
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
    # Total cardiac serum flow through the lungs, tumours included.
    fq_total <- fserum + fq_tumor1 + fq_tumor2 + fq_tumorrest

    # ---- Permeability-surface-area products (Table S2) --------------------
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
    # Brain PS = 0 (equation S4 note); the brain has no interstitial state.

    # ---- Renal filtration and excretion (Table S2) ------------------------
    gfr    <- (TER_MAG3 / 1000) / 3 * 20 / 15  # 51Cr-EDTA GFR from measured TER, ref (16)
    f_fil  <- gfr * phi_gfr                    # filtration scaled for molecular size
    f_exc  <- f_fil * f_ex                     # excreted fraction of the filtered load
    f_reab <- f_fil - f_exc                    # returned to the vasculature via the cells

    # ---- Total binding sites, R_i,0 = [R_i,0] * V_i,total (Table S2) ------
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
    koff <- kd * kon                           # Table S2: k_off = K_D * k_on
    krel_tumor1     <- krel_tumor
    krel_tumor2     <- krel_tumor
    krel_tumorrest  <- krel_tumor
    krel_liver      <- krel_kidney
    krel_spleen     <- krel_kidney
    krel_gi         <- krel_kidney
    krel_prostate   <- krel_kidney

    # =====================================================================
    # 2. ODE system. Equations S1-S13. Each state exists twice: _lab is
    #    radiolabelled [177Lu]Lu-PSMA I&T and _unl is unlabelled peptide.
    #    Physical decay (lambda_phy) is the only coupling between the two
    #    circulations: it removes labelled peptide from a state and creates
    #    unlabelled peptide in the SAME state.
    # =====================================================================

    # ---------------------------------------------------------------------
    # Radiolabelled circulation
    # ---------------------------------------------------------------------

    # Free binding sites, equation S1: R_0,i = R_i + RP_i + RP*_i.
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

    # Arteries, equation S8
    d/dt(part_lab) <- fq_total * (pv_lungs_lab / vv_lungs - part_lab / vv_art) - lambda_phy * part_lab

    # Veins, equation S7. Spleen and GI drain into the liver rather than
    # into the veins, so the liver returns F_L + F_S + F_GI (Fig. S1 legend).
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

    # Peptide bound to serum protein, equation S13
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
    # Liver, equation S4 with the triple inflow of Fig. S1a: arterial,
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
    # Kidneys, equation S6: filtration leaves the vasculature and the
    # intratubular cells return the non-excreted fraction.
    d/dt(pv_kidney_lab) <- - pv_kidney_lab / vv_kidney * (f_fil + fq_kidney) +
      fq_kidney * part_lab / vv_art + f_reab * pintra_kidney_lab / vintra_kidney - lambda_phy * pv_kidney_lab
    d/dt(pv_muscle_lab) <- ps_muscle * (pint_muscle_lab / vi_muscle - pv_muscle_lab / vv_muscle) +
      fq_muscle * (part_lab / vv_art - pv_muscle_lab / vv_muscle) - lambda_phy * pv_muscle_lab
    # Lungs, equation S5: perfused from the veins, drained to the arteries.
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
    # Brain, equation S4 with PS = 0: no interstitial exchange.
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
    # Kidneys, equation S9: fed by glomerular filtration, not by PS.
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

    # Unspecific uptake into the proximal tubular cells, equation S12
    d/dt(pintra_kidney_lab) <- f_reab * (pint_kidney_lab / vi_kidney - pintra_kidney_lab / vintra_kidney) - lambda_phy * pintra_kidney_lab

    # Cumulative urinary excretion (bookkeeping state, not in the paper; it
    # collects the excreted fraction f_ex of the filtered load so that mass
    # balance can be checked).
    d/dt(urine_lab) <- f_exc * pint_kidney_lab / vi_kidney - lambda_phy * urine_lab

    # Cumulative peptide lost by release from the internalised pool, which the
    # paper assumes is cleared from the body directly (Mechanism and Structure
    # of PBPK Model). Bookkeeping state, not in the paper; it closes the mass
    # balance exactly.
    d/dt(cleared_lab) <-
      krel_tumor1 * pintern_tumor1_lab + krel_tumor2 * pintern_tumor2_lab +
      krel_tumorrest * pintern_tumorrest_lab + krel_salgland * pintern_salgland_lab +
      krel_liver * pintern_liver_lab + krel_spleen * pintern_spleen_lab +
      krel_gi * pintern_gi_lab + krel_prostate * pintern_prostate_lab +
      krel_kidney * pintern_kidney_lab - lambda_phy * cleared_lab

    # Peptide bound to PSMA on the cell surface, equation S3
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

    # Internalised peptide, equation S2. Released peptide and free 177Lu are
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

    # Arteries, equation S8
    d/dt(part_unl) <- fq_total * (pv_lungs_unl / vv_lungs - part_unl / vv_art) + lambda_phy * part_lab

    # Veins, equation S7. Spleen and GI drain into the liver rather than
    # into the veins, so the liver returns F_L + F_S + F_GI (Fig. S1 legend).
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

    # Peptide bound to serum protein, equation S13
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
    # Liver, equation S4 with the triple inflow of Fig. S1a: arterial,
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
    # Kidneys, equation S6: filtration leaves the vasculature and the
    # intratubular cells return the non-excreted fraction.
    d/dt(pv_kidney_unl) <- - pv_kidney_unl / vv_kidney * (f_fil + fq_kidney) +
      fq_kidney * part_unl / vv_art + f_reab * pintra_kidney_unl / vintra_kidney + lambda_phy * pv_kidney_lab
    d/dt(pv_muscle_unl) <- ps_muscle * (pint_muscle_unl / vi_muscle - pv_muscle_unl / vv_muscle) +
      fq_muscle * (part_unl / vv_art - pv_muscle_unl / vv_muscle) + lambda_phy * pv_muscle_lab
    # Lungs, equation S5: perfused from the veins, drained to the arteries.
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
    # Brain, equation S4 with PS = 0: no interstitial exchange.
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
    # Kidneys, equation S9: fed by glomerular filtration, not by PS.
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

    # Unspecific uptake into the proximal tubular cells, equation S12
    d/dt(pintra_kidney_unl) <- f_reab * (pint_kidney_unl / vi_kidney - pintra_kidney_unl / vintra_kidney) + lambda_phy * pintra_kidney_lab

    # Cumulative urinary excretion (bookkeeping state, not in the paper; it
    # collects the excreted fraction f_ex of the filtered load so that mass
    # balance can be checked).
    d/dt(urine_unl) <- f_exc * pint_kidney_unl / vi_kidney + lambda_phy * urine_lab

    # Cumulative peptide lost by release from the internalised pool, which the
    # paper assumes is cleared from the body directly (Mechanism and Structure
    # of PBPK Model). Bookkeeping state, not in the paper; it closes the mass
    # balance exactly.
    d/dt(cleared_unl) <-
      krel_tumor1 * pintern_tumor1_unl + krel_tumor2 * pintern_tumor2_unl +
      krel_tumorrest * pintern_tumorrest_unl + krel_salgland * pintern_salgland_unl +
      krel_liver * pintern_liver_unl + krel_spleen * pintern_spleen_unl +
      krel_gi * pintern_gi_unl + krel_prostate * pintern_prostate_unl +
      krel_kidney * pintern_kidney_unl + lambda_phy * cleared_lab

    # Peptide bound to PSMA on the cell surface, equation S3
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

    # Internalised peptide, equation S2. Released peptide and free 177Lu are
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
    # 3. Observed quantities. Equations S20-S23 give the labelled-peptide
    #    amount seen by the gamma camera in each region of interest: the
    #    organ self-content plus that organ's share of the
    #    serum-protein-bound pool, and for the tumour ROIs a background
    #    correction for labelled peptide in muscle.
    # =====================================================================
    vprp_denom <- vserum + vv_tumor1 + vv_tumor2 + vv_tumorrest
    prp_muscle <- pv_muscle_lab + pint_muscle_lab + pprp_lab * vv_muscle / vprp_denom

    # Equation S20. Note that the printed equation does not include the
    # intratubular pool pintra_kidney_lab; it is transcribed as printed.
    Pkidney <- pv_kidney_lab + pint_kidney_lab + rp_kidney_lab + pintern_kidney_lab +
      pprp_lab * vv_kidney / vprp_denom
    # Equation S21
    Psalgland <- pv_salgland_lab + pint_salgland_lab + rp_salgland_lab + pintern_salgland_lab +
      pprp_lab * vv_salgland / vprp_denom
    # Equation S22
    Ptumor1 <- pv_tumor1_lab + pint_tumor1_lab + rp_tumor1_lab + pintern_tumor1_lab +
      pprp_lab * vv_tumor1 / vprp_denom + bgcorr_tumor1 * prp_muscle
    # Equation S23
    Ptumor2 <- pv_tumor2_lab + pint_tumor2_lab + rp_tumor2_lab + pintern_tumor2_lab +
      pprp_lab * vv_tumor2 / vprp_denom + bgcorr_tumor2 * prp_muscle
    # The lumped remaining lesions follow the same form, with no ROI background
    # correction because they are not delineated on the gamma camera.
    Ptumorrest <- pv_tumorrest_lab + pint_tumorrest_lab + rp_tumorrest_lab +
      pintern_tumorrest_lab + pprp_lab * vv_tumorrest / vprp_denom
    Ptumortotal <- Ptumor1 + Ptumor2 + Ptumorrest

    # Whole-body retained labelled peptide, and the venous serum concentration.
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
    # Receptor occupancy in the two delineated lesions (Figs. 6 and 7).
    occ_tumor1 <- (rp_tumor1_lab + rp_tumor1_unl) / r0_tumor1
    occ_tumor2 <- (rp_tumor2_lab + rp_tumor2_unl) / r0_tumor2
  })
}
