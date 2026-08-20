Budiansah_2025_dotatate_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, NONMEM 7.5.1). Seventeen-region whole-body",
    "physiologically-based pharmacokinetic model for the somatostatin-receptor",
    "(sst2) targeted peptide DOTATATE, fitted within a non-linear mixed-effects",
    "framework to [111In]In-DOTA-TATE biokinetics in eight patients (four",
    "neuroendocrine tumours, four meningiomas) undergoing pre-therapeutic",
    "dosimetry for [90Y]Y-DOTA-TATE peptide receptor radionuclide therapy.",
    "114 ODE states: two parallel circulations, one carrying total peptide",
    "(labelled plus unlabelled, in nmol) and one carrying labelled peptide as a",
    "fraction of the injected activity, coupled through a shared, conserved pool",
    "of sst2 binding sites. Ten sst2-positive regions (tumour, liver, spleen,",
    "kidneys, GI plus pancreas, muscle, red marrow, prostate or uterus, adrenals",
    "and a lumped rest of body) each carry vascular, interstitial,",
    "receptor-bound and internalised sub-compartments; the kidneys add a fifth",
    "intracellular pool fed by glomerular filtration, five sst2-negative regions",
    "(lungs, skin, adipose, heart, bone) carry vascular and interstitial pools,",
    "the brain is vascular-only because its permeability-surface product is set",
    "to zero, and arteries, veins and an irreversibly serum-protein-bound pool",
    "complete the structure. Spleen and GI drain into the liver through the",
    "portal circulation. Five sst2 densities, two internalisation rates, two",
    "degradation rates and the serum-protein binding rate were estimated with",
    "exponential interindividual variability; every physiological parameter is",
    "fixed to a measured, allometrically scaled or literature value. Six",
    "observations (serum, whole body, kidneys, spleen, liver, tumour) each carry",
    "their own proportional residual error.",
    sep = " "
  )
  reference <- paste(
    "Budiansah I, Hardiansyah D, Riana A, Pawiro SA, Beer AJ, Glatting G.",
    "Accuracy and precision analyses of single-time-point dosimetry utilising",
    "physiologically-based pharmacokinetic modelling and non-linear",
    "mixed-effects modelling. EJNMMI Phys. 2025;12(1):26.",
    "doi:10.1186/s40658-025-00726-7.",
    "Structural equations and fixed physiological parameters from the upstream",
    "framework paper: Kletting P, Kull T, Maass C, Malik N, Luster M, Beer AJ,",
    "Glatting G. Optimized peptide amount and activity for 90Y-labeled DOTATATE",
    "therapy. J Nucl Med. 2016;57(4):503-8. doi:10.2967/jnumed.115.164699.",
    sep = " "
  )
  vignette <- "Budiansah_2025_dotatate_pbpk"
  units <- list(time = "min", dosing = "nmol", concentration = "fraction of injected activity (serum: per L)")
  # The intravenous infusion enters the venous pool of each circulation, not a
  # state named depot or central, so buildModelDb()'s two-name heuristic cannot
  # infer it and the registry would otherwise report no dosing compartment.
  dosing <- c("pven_tot", "pven_lab")

  covariateData <- list(
    WT = list(
      description = "Total body weight, used to allometrically scale every organ volume that was not measured directly (Table S1, all `x * BW / 71` entries) and to set the total body volume via the 1 mL = 1 g assumption.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table S1 lists body weight as 'individually measured', but the",
        "per-patient values are not tabulated in this paper or in the upstream",
        "Kletting 2016 framework paper (whose Table 1 reports body surface area,",
        "not weight). Downstream users must supply it. The reference individual",
        "is the 71 kg ICRP-23 adult that every `x * BW / 71` scaling is written",
        "against. Baseline value, held constant per subject.",
        sep = " "
      ),
      source_name = "BW"
    ),
    BSA = list(
      description = "Body surface area, which sets the total body serum volume V_P (Table S1: 2.8 * (1 - H) * BSA in men and 2.4 * (1 - H) * BSA in women).",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reported per patient in Kletting 2016 Table 1 (cohort range",
        "1.57-2.05 m^2). V_P in turn sets total serum flow, every vascular",
        "sub-volume expressed as a fraction of V_P, and the arterial and venous",
        "pools, so BSA propagates into essentially every flow term. Baseline",
        "value, held constant per subject.",
        sep = " "
      ),
      source_name = "BSA"
    ),
    HCT = list(
      description = "Hematocrit, which converts whole-blood quantities to serum quantities in V_P and in the perfusion of the tumour, adrenals and prostate or uterus.",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table S1 lists hematocrit as 'individually measured' and carries it as",
        "a fraction ('unity'); it is stored here in the canonical percent and",
        "divided by 100 inside model(). Per-patient values are not tabulated in",
        "any on-disk source, so downstream users must supply it. Baseline value,",
        "held constant per subject.",
        sep = " "
      ),
      source_name = "H"
    ),
    SEXF = list(
      description = "Sex, 1 = female. Switches the sex-specific reproductive organ (uterus in women, prostate in men) and the coefficient of the total body serum volume.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "Enters in five places, all from Table S1: the V_P coefficient",
        "(2.8 male, 2.4 female); the reproductive organ volume (0.016 * BW / 71",
        "for prostate, 0.080 * BW / 71 for uterus); its vascular fraction (0.04",
        "prostate, 0.07 uterus); its interstitial fraction (0.25 prostate, 0.5",
        "uterus); its serum flow density (0.18 prostate, 1 uterus); its",
        "permeability-surface density (0.1 prostate, 0.2 uterus); and its sst2",
        "density relative to the kidneys (0.26 prostate, 0.092 uterus). The",
        "cohort was three women and five men in the upstream Kletting 2016",
        "Table 1.",
        sep = " "
      ),
      source_name = "Sex"
    ),
    CRCL = list(
      description = "Glomerular filtration rate measured with 51Cr-EDTA, which sets the filtration flow F_fil = GFR * phi and hence the entire renal handling of the peptide.",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Kletting 2016 Table 1 reports the measured 51Cr-EDTA GFR in absolute",
        "L/min (cohort range 0.028-0.13 L/min). It is carried here in the",
        "canonical BSA-normalised mL/min/1.73 m^2 and converted back to an",
        "absolute L/min inside model() as CRCL * BSA / 1.73 / 1000, because the",
        "filtration flow in the ODEs is an absolute volumetric flow rather than",
        "a normalised one. Baseline value, held constant per subject.",
        sep = " "
      ),
      source_name = "GFR"
    ),
    ORGVOL_KIDNEY = list(
      description = "Measured total volume of both kidneys, delineated on CT (Table S1: 'individually measured').",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Split into vascular (5.5 %), interstitial (15 %) and intracellular",
        "(two thirds of the remainder, Table S1 footnote c) sub-volumes, and",
        "multiplies the fitted sst2 density to give the renal receptor number.",
        "Cohort range 125-233 mL in the upstream Kletting 2016 Table 1.",
        sep = " "
      ),
      source_name = "Kidneys (measured volume)"
    ),
    ORGVOL_LIVER = list(
      description = "Measured total liver volume, delineated on CT (Table S1: 'individually measured').",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Split into vascular (8.5 %) and interstitial (20 %) sub-volumes and",
        "multiplies the fitted hepatic sst2 density. Cohort range 1500-4876 mL",
        "in the upstream Kletting 2016 Table 1.",
        sep = " "
      ),
      source_name = "Liver (measured volume)"
    ),
    ORGVOL_SPLEEN = list(
      description = "Measured total spleen volume, delineated on CT (Table S1: 'individually measured').",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Split into vascular (12 %) and interstitial (20 %) sub-volumes and",
        "multiplies the fitted splenic sst2 density. Cohort range 110-320 mL in",
        "the upstream Kletting 2016 Table 1; one patient in the Budiansah 2025",
        "cohort had a splenectomy and contributed no spleen data.",
        sep = " "
      ),
      source_name = "Spleen (measured volume)"
    ),
    TUM_VOL = list(
      description = "Total volume of the delineated tumour lesions, summed over the two largest lesions, measured on CT.",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Methods, 'Biokinetic data': the regions of interest of all visible",
        "lesions were drawn, the two largest were delineated and quantified, and",
        "the sum of their biokinetic data was used. Reported in mL in the",
        "upstream Kletting 2016 Table 1 (cohort range 2.5-134 mL over the eight",
        "patients retained here) and carried in the canonical mm^3, so multiply",
        "the published mL by 1000. Sets the tumour vascular and interstitial",
        "sub-volumes, its perfusion, its permeability-surface product and its",
        "receptor number. Baseline value, held constant per subject.",
        sep = " "
      ),
      source_name = "Tumor 1 + Tumor 2 (measured volume)"
    ),
    TUMTP_NET = list(
      description = "Tumour-type indicator, 1 = neuroendocrine tumour, 0 = meningioma.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = meningioma",
      notes = paste(
        "The tumour is the only region whose physiology is tumour-type",
        "specific. Table S1 gives four switched values: the interstitial",
        "fraction (0.3 NET, 0.23 meningioma), the vascular fraction (0.1 NET,",
        "0.11 meningioma), the serum flow density (1.0 NET, 0.9 meningioma",
        "mL/min/g) and the permeability-surface density (0.2 NET, 0.31",
        "meningioma mL/min/g). The perfusion values are the setting II prior",
        "knowledge that the paper fixed rather than fitted (Study workflow, and",
        "Table S2). The cohort was four NET and four meningioma patients.",
        sep = " "
      ),
      source_name = "Disease"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 8,
    n_studies = 1,
    disease_state = "Four patients with metastasised neuroendocrine tumours and four with meningiomas, all scheduled for two to three cycles of [90Y]Y-DOTA-TATE peptide receptor radionuclide therapy.",
    dose_range = "75 +/- 10 nmol DOTATATE labelled with 140 +/- 14 MBq 111In, given as a 51 +/- 8 min intravenous infusion",
    co_medication = "Arginine and lysine (1000 mL, 2.5 % infusion) over 2 h, starting 0.5 h before the [111In]In-DOTA-TATE administration, for renal protection.",
    regions = "Germany (Ulm University)",
    sampling = "Serum at 5 and 15 min, 0.5, 1, 2 and 4 h, and 1, 2 and 3 d; planar whole-body scintigraphy at 2, 4, 24, 48 and 72 h after injection.",
    notes = paste(
      "The eight patients are a subset of the nine-patient cohort of the",
      "upstream Kletting 2016 framework paper, whose Table 1 is the only",
      "on-disk source of the per-patient body surface areas, measured GFRs,",
      "organ volumes, tumour volumes and sexes. Body weight and hematocrit are",
      "described as individually measured but are not tabulated anywhere on",
      "disk. [111In]In-DOTA-TATE biokinetics were used as a surrogate for",
      "[90Y]Y-DOTA-TATE, so the model is fitted on 111In decay and simulated on",
      "90Y decay by changing lambda_phy.",
      sep = " "
    )
  )

  compartmentData <- list(
    pv_tumor_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pv_liver_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_spleen_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_gi_tot             = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_muscle_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_redmarrow_tot      = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_reprod_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_adrenals_tot       = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_rest_tot           = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_kidney_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_lungs_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_skin_tot           = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_adipose_tot        = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_heart_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_bone_tot           = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pv_brain_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pint_tumor_tot        = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pint_liver_tot        = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_spleen_tot       = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_gi_tot           = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_muscle_tot       = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_redmarrow_tot    = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_reprod_tot       = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_adrenals_tot     = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_rest_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_kidney_tot       = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_lungs_tot        = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_skin_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_adipose_tot      = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_heart_tot        = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pint_bone_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_tumor_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    rp_liver_tot          = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_spleen_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_gi_tot             = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_muscle_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_redmarrow_tot      = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_reprod_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_adrenals_tot       = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_rest_tot           = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    rp_kidney_tot         = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_tumor_tot     = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tumor", verified = TRUE),
    pintern_liver_tot     = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_spleen_tot    = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_gi_tot        = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_muscle_tot    = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_redmarrow_tot = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_reprod_tot    = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_adrenals_tot  = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_rest_tot      = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintern_kidney_tot    = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    pintra_kidney_tot     = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "tissue", verified = TRUE),
    part_tot              = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pven_tot              = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    pprp_tot              = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "serum", verified = TRUE),
    urine_tot             = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "urine", verified = TRUE),
    cleared_tot           = list(analyte = "DOTATATE (labelled plus unlabelled)", units = "nmol", specimen = "not applicable", verified = TRUE),
    pv_tumor_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tumor", verified = TRUE),
    pv_liver_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_spleen_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_gi_lab             = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_muscle_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_redmarrow_lab      = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_reprod_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_adrenals_lab       = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_rest_lab           = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_kidney_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_lungs_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_skin_lab           = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_adipose_lab        = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_heart_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_bone_lab           = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pv_brain_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pint_tumor_lab        = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tumor", verified = TRUE),
    pint_liver_lab        = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_spleen_lab       = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_gi_lab           = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_muscle_lab       = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_redmarrow_lab    = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_reprod_lab       = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_adrenals_lab     = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_rest_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_kidney_lab       = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_lungs_lab        = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_skin_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_adipose_lab      = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_heart_lab        = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pint_bone_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_tumor_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tumor", verified = TRUE),
    rp_liver_lab          = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_spleen_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_gi_lab             = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_muscle_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_redmarrow_lab      = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_reprod_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_adrenals_lab       = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_rest_lab           = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    rp_kidney_lab         = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_tumor_lab     = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tumor", verified = TRUE),
    pintern_liver_lab     = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_spleen_lab    = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_gi_lab        = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_muscle_lab    = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_redmarrow_lab = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_reprod_lab    = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_adrenals_lab  = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_rest_lab      = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintern_kidney_lab    = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    pintra_kidney_lab     = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "tissue", verified = TRUE),
    part_lab              = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pven_lab              = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    pprp_lab              = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "serum", verified = TRUE),
    urine_lab             = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "urine", verified = TRUE),
    cleared_lab           = list(analyte = "[111In]In-DOTA-TATE", units = "fraction of injected activity", specimen = "not applicable", verified = TRUE)
  )


  ini({
    # =====================================================================
    # Estimated parameters, Table 1 of Budiansah 2025 (parameter setting II,
    # the set selected by an Akaike weight of nearly 100 %). Every one carries
    # exponential interindividual variability; the paper reports the random
    # effect as a coefficient of variation in percent, converted here to a
    # log-scale variance as omega^2 = log(CV^2 + 1).
    # =====================================================================
    lrdens_kidney <- log(4.55)   ; label("Log sst2 receptor density, kidneys (nmol/L)")            # Table 1, [R_K,0] = 4.55 (15 % RSE)
    lrdens_liver  <- log(0.94)   ; label("Log sst2 receptor density, liver (nmol/L)")              # Table 1, [R_L,0] = 0.94 (15 % RSE)
    lrdens_spleen <- log(7.37)   ; label("Log sst2 receptor density, spleen (nmol/L)")             # Table 1, [R_S,0] = 7.37 (19 % RSE)
    lrdens_tumor  <- log(15.9)   ; label("Log sst2 receptor density, tumour (nmol/L)")             # Table 1, [R_TU,0] = 15.9 (30 % RSE)
    lrdens_rest   <- log(0.36)   ; label("Log sst2 receptor density, rest of body (nmol/L)")       # Table 1, [R_Rest,0] = 0.36 (30 % RSE)
    lkpr          <- log(5.50e-4); label("Log binding rate of peptide to serum protein (1/min)")   # Table 1, k_on,Alb = 5.50e-4 (18 % RSE)
    lkint_nt      <- log(2.36e-3); label("Log sst2 internalisation rate, normal tissue (1/min)")   # Table 1, lambda_int,K = 2.36e-3 (11 % RSE); Table S1 sets lambda_int,L = lambda_int,S = lambda_int,NT = lambda_int,K
    lkint_tumor   <- log(1.39e-3); label("Log sst2 internalisation rate, tumour (1/min)")          # Table 1, lambda_int,TU = 1.39e-3 (14 % RSE)
    lkdeg_nt      <- log(1.01e-4); label("Log degradation rate, normal tissue (1/min)")            # Table 1, lambda_deg,NT = 1.01e-4 (10 % RSE); Table S1 sets lambda_deg,L = lambda_deg,S = lambda_deg,K = lambda_deg,NT
    lkdeg_tumor   <- log(1.15e-4); label("Log degradation rate, tumour (1/min)")                   # Table 1, lambda_deg,TU = 1.15e-4 (24 % RSE)

    etalrdens_kidney ~ 0.184403  ; label("IIV on log sst2 density, kidneys")        # Table 1 RE 45 % CV: log(0.45^2 + 1)
    etalrdens_liver  ~ 0.169658  ; label("IIV on log sst2 density, liver")          # Table 1 RE 43 % CV
    etalrdens_spleen ~ 0.239332  ; label("IIV on log sst2 density, spleen")         # Table 1 RE 52 % CV
    etalrdens_tumor  ~ 0.703147  ; label("IIV on log sst2 density, tumour")         # Table 1 RE 101 % CV
    etalrdens_rest   ~ 0.663152  ; label("IIV on log sst2 density, rest of body")   # Table 1 RE 97 % CV
    etalkpr          ~ 0.239332  ; label("IIV on log serum-protein binding rate")   # Table 1 RE 52 % CV
    etalkint_nt      ~ 0.075478  ; label("IIV on log internalisation, normal")      # Table 1 RE 28 % CV
    etalkint_tumor   ~ 0.128305  ; label("IIV on log internalisation, tumour")      # Table 1 RE 37 % CV
    etalkdeg_nt      ~ 0.035464  ; label("IIV on log degradation, normal")          # Table 1 RE 19 % CV
    etalkdeg_tumor   ~ 0.121864  ; label("IIV on log degradation, tumour")          # Table 1 RE 36 % CV

    # =====================================================================
    # Fixed drug-dependent parameters, Table S1 (adapted from the upstream
    # Kletting 2016 supplemental Table 1).
    # =====================================================================
    koff       <- fixed(0.04)    ; label("sst2 dissociation rate (1/min)")                       # Table S1
    kdiss      <- fixed(0.4)     ; label("sst2 dissociation constant (nmol/L)")                  # Table S1, K_D; k_on = k_off / K_D
    lambda_phy <- fixed(1.71e-4) ; label("Physical decay rate of the radionuclide (1/min)")      # Table S1, 111In; set to 1.80e-4 to simulate 90Y
    phi_gfr    <- fixed(0.66)    ; label("Ratio of sieving coefficients, DOTATATE to 51Cr-EDTA") # Table S1, phi
    f_ex_frac  <- fixed(0.98)    ; label("Excreted fraction of the filtered load")               # Table S1, f_ex; the remaining 2 % re-enters the serum

    # ---- Tumour-type-specific physiology (Table S1) ----------------------
    vfrac_tu_int_net <- fixed(0.3)  ; label("Tumour interstitial volume fraction, NET")           # Table S1
    vfrac_tu_int_men <- fixed(0.23) ; label("Tumour interstitial volume fraction, meningioma")    # Table S1
    vfrac_tu_v_net   <- fixed(0.1)  ; label("Tumour vascular volume fraction, NET")               # Table S1
    vfrac_tu_v_men   <- fixed(0.11) ; label("Tumour vascular volume fraction, meningioma")        # Table S1
    fperf_tu_net     <- fixed(1.0)  ; label("Tumour serum flow density, NET (mL/min/g)")          # Table S1; setting II prior knowledge
    fperf_tu_men     <- fixed(0.9)  ; label("Tumour serum flow density, meningioma (mL/min/g)")   # Table S1; setting II prior knowledge
    kps_tu_net       <- fixed(0.2)  ; label("Tumour permeability-surface density, NET (mL/min/g)")        # Table S1
    kps_tu_men       <- fixed(0.31) ; label("Tumour permeability-surface density, meningioma (mL/min/g)") # Table S1

    # ---- Permeability-surface and perfusion densities (Table S1) ---------
    kps_muscle   <- fixed(0.02) ; label("Permeability-surface density, muscle-like tissue (mL/min/g)")  # Table S1, k_MUS; GI, skin, adipose, red marrow, heart, bone and rest are set equal to it
    kps_liver    <- fixed(2)    ; label("Permeability-surface density, liver-like tissue (mL/min/g)")   # Table S1, k_L = k_MUS * 100; spleen, lungs and adrenals share it
    kps_prost    <- fixed(0.1)  ; label("Permeability-surface density, prostate (mL/min/g)")            # Table S1
    kps_uterus   <- fixed(0.2)  ; label("Permeability-surface density, uterus (mL/min/g)")              # Table S1
    fperf_prost  <- fixed(0.18) ; label("Serum flow density, prostate (mL/min/g)")                      # Table S1, f_PRO
    fperf_uterus <- fixed(1)    ; label("Serum flow density, uterus (mL/min/g)")                        # Table S1, f_UT
    fperf_adren  <- fixed(6)    ; label("Serum flow density, adrenals (mL/min/g)")                      # Table S1, f_AD

    # ---- sst2 densities of the non-delineated organs, relative to the
    #      kidneys (Table S1; regional sst2 expression of Boy et al.) ------
    rratio_prost  <- fixed(0.26)   ; label("sst2 density, prostate, relative to kidneys")    # Table S1
    rratio_uterus <- fixed(0.092)  ; label("sst2 density, uterus, relative to kidneys")      # Table S1
    rratio_adren  <- fixed(1.65)   ; label("sst2 density, adrenals, relative to kidneys")    # Table S1
    rratio_muscle <- fixed(0.0056) ; label("sst2 density, muscle, relative to kidneys")      # Table S1
    rratio_gi     <- fixed(0.16)   ; label("sst2 density, GI plus pancreas, relative to kidneys")  # Table S1
    rratio_rm     <- fixed(0.028)  ; label("sst2 density, red marrow, relative to kidneys")  # Table S1

    # =====================================================================
    # Residual error. Equation 2 is a purely proportional model, and setting
    # II gives each observed region its own error parameter a, reported in
    # Table 1 as a percent fractional standard deviation.
    # =====================================================================
    propSd_Aserum      <- 0.23  ; label("Proportional residual error, blood serum")   # Table 1, a = 23 % (8 % RSE)
    propSd_Abody  <- 0.088 ; label("Proportional residual error, whole body")    # Table 1, a = 8.8 % (14 % RSE)
    propSd_Akidney     <- 0.053 ; label("Proportional residual error, kidneys")       # Table 1, a = 5.3 % (11 % RSE)
    propSd_Aspleen     <- 0.067 ; label("Proportional residual error, spleen")        # Table 1, a = 6.7 % (14 % RSE)
    propSd_Aliver      <- 0.10  ; label("Proportional residual error, liver")         # Table 1, a = 10 % (9 % RSE)
    propSd_Atumor      <- 0.060 ; label("Proportional residual error, tumour")        # Table 1, a = 6.0 % (10 % RSE)
  })

  model({
    # =====================================================================
    # 1. Individual parameters. Equation 1 of Budiansah 2025 is an
    #    exponential interindividual error model, P_ki = theta_k * exp(eta_ki).
    # =====================================================================
    rdens_kidney <- exp(lrdens_kidney + etalrdens_kidney)
    rdens_liver  <- exp(lrdens_liver + etalrdens_liver)
    rdens_spleen <- exp(lrdens_spleen + etalrdens_spleen)
    rdens_tumor  <- exp(lrdens_tumor + etalrdens_tumor)
    rdens_rest   <- exp(lrdens_rest + etalrdens_rest)
    kpr          <- exp(lkpr + etalkpr)
    kint_nt      <- exp(lkint_nt + etalkint_nt)
    kint_tumor   <- exp(lkint_tumor + etalkint_tumor)
    kdeg_nt      <- exp(lkdeg_nt + etalkdeg_nt)
    kdeg_tumor   <- exp(lkdeg_tumor + etalkdeg_tumor)

    kon <- koff / kdiss                        # Table S1, k_on = k_off / K_D (L/nmol/min)

    # =====================================================================
    # 2. Anatomy and physiology, all from Table S1. Volumes are in L, flows
    #    in L/min, amounts in nmol and time in min. Table S1 reports the
    #    perfusion and permeability densities in mL/min/g and assumes
    #    1 mL = 1 g, so a density in mL/min/g multiplied by a volume in L is
    #    already a flow in L/min and needs no conversion factor.
    # =====================================================================
    hct    <- HCT / 100                                        # Table S1 carries H as a fraction
    vserum <- (2.8 * (1 - SEXF) + 2.4 * SEXF) * (1 - hct) * BSA  # V_P, Table S1
    fserum <- 1.23 * vserum                                    # F, Table S1 footnote b

    # Kletting 2016 Table 1 reports the 51Cr-EDTA GFR in absolute L/min; CRCL
    # is carried BSA-normalised, so undo the normalisation here.
    gfr_abs <- CRCL * BSA / 1.73 / 1000
    f_fil   <- gfr_abs * phi_gfr               # Table S1, F_fil = GFR * phi
    f_exc   <- f_fil * f_ex_frac               # Table S1, F_ex = F_fil * f_ex
    f_reab  <- f_fil - f_exc                   # the 2 % that returns to the serum via the tubular cells

    # ---- Tumour-type switches (Table S1) ---------------------------------
    vfrac_tu_int <- vfrac_tu_int_net * TUMTP_NET + vfrac_tu_int_men * (1 - TUMTP_NET)
    vfrac_tu_v   <- vfrac_tu_v_net * TUMTP_NET + vfrac_tu_v_men * (1 - TUMTP_NET)
    fperf_tumor  <- fperf_tu_net * TUMTP_NET + fperf_tu_men * (1 - TUMTP_NET)
    kps_tumor    <- kps_tu_net * TUMTP_NET + kps_tu_men * (1 - TUMTP_NET)

    # ---- Sex switches for the reproductive organ (Table S1) --------------
    vcoef_reprod  <- 0.016 * (1 - SEXF) + 0.080 * SEXF   # V_PRO,total / V_UT,total coefficient
    vfrac_reprod_v   <- 0.04 * (1 - SEXF) + 0.07 * SEXF  # vascular fraction
    vfrac_reprod_int <- 0.25 * (1 - SEXF) + 0.5 * SEXF   # interstitial fraction
    fperf_reprod  <- fperf_prost * (1 - SEXF) + fperf_uterus * SEXF
    kps_reprod_d  <- kps_prost * (1 - SEXF) + kps_uterus * SEXF
    rratio_reprod <- rratio_prost * (1 - SEXF) + rratio_uterus * SEXF

    # ---- Total region volumes (L), Table S1 ------------------------------
    vtot_kidney    <- ORGVOL_KIDNEY / 1000     # measured
    vtot_liver     <- ORGVOL_LIVER / 1000      # measured
    vtot_spleen    <- ORGVOL_SPLEEN / 1000     # measured
    vtot_tumor     <- TUM_VOL / 1e6            # measured, canonical mm^3 -> L
    vtot_reprod    <- vcoef_reprod * WT / 71
    vtot_lungs     <- 1      * WT / 71
    vtot_adrenals  <- 0.014  * WT / 71
    vtot_muscle    <- 30.078 * WT / 71
    vtot_gi        <- (0.385 + 0.548 + 0.104 + 0.105) * WT / 71   # GI plus pancreas
    vtot_skin      <- 3.408  * WT / 71
    vtot_adipose   <- 13.465 * WT / 71
    vtot_redmarrow <- 1.1    * WT / 71
    vtot_bone      <- 10.165 * WT / 71 - vtot_redmarrow
    vtot_heart     <- 0.341  * WT / 71
    vtot_brain     <- 1.45   * WT / 71
    # V_REST,total = V_BW - sum over all regions except the tumour, with
    # V_BW taken from body weight through the 1 mL = 1 g assumption.
    vtot_rest <- WT - vtot_kidney - vtot_liver - vtot_spleen - vtot_reprod -
      vtot_lungs - vtot_adrenals - vtot_muscle - vtot_gi - vtot_skin -
      vtot_adipose - vtot_redmarrow - vtot_bone - vtot_heart - vtot_brain

    # ---- Vascular (serum) sub-volumes (L), Table S1 ----------------------
    vv_lungs     <- 0.105 * vserum
    vv_muscle    <- 0.14  * vserum
    vv_gi        <- 0.076 * vserum
    vv_skin      <- 0.03  * vserum
    vv_adipose   <- 0.05  * vserum
    vv_redmarrow <- 0.04  * vserum
    vv_bone      <- 0.07  * vserum - vv_redmarrow
    vv_heart     <- 0.01  * vserum
    vv_brain     <- 0.012 * vserum
    vv_art       <- 0.06 * vserum + 0.045 * vserum   # arterial serum plus half the cardiac serum
    vv_ven       <- 0.18 * vserum + 0.045 * vserum   # venous serum plus half the cardiac serum
    vv_liver     <- 0.085 * vtot_liver
    vv_spleen    <- 0.12  * vtot_spleen
    vv_kidney    <- 0.055 * vtot_kidney
    vv_tumor     <- vfrac_tu_v * vtot_tumor
    vv_reprod    <- vfrac_reprod_v * (1 - hct) * vtot_reprod
    vv_adrenals  <- 0.03 * (1 - hct) * vtot_adrenals
    # V_REST,v is the serum left once every other serum pool is assigned.
    vv_rest <- vserum - vv_art - vv_ven - vv_lungs - vv_muscle - vv_gi -
      vv_skin - vv_adipose - vv_redmarrow - vv_bone - vv_heart - vv_brain -
      vv_liver - vv_spleen - vv_kidney - vv_reprod - vv_adrenals

    # ---- Interstitial sub-volumes (L), Table S1 --------------------------
    vi_liver     <- 0.2  * vtot_liver
    vi_spleen    <- 0.2  * vtot_spleen
    vi_kidney    <- 0.15 * vtot_kidney
    vi_tumor     <- vfrac_tu_int * vtot_tumor
    vi_reprod    <- vfrac_reprod_int * vtot_reprod
    vi_adrenals  <- 0.24 * vtot_adrenals              # the salivary gland value is used
    vi_lungs     <- vv_lungs     * 5.5                # alpha_LU
    vi_muscle    <- vv_muscle    * 5.9                # alpha_MUS
    vi_gi        <- vv_gi        * 8.8                # alpha_GI
    vi_skin      <- vv_skin      * 8.9                # alpha_SKIN
    vi_adipose   <- vv_adipose   * 15.5               # alpha_ADI
    vi_redmarrow <- vv_redmarrow * 3.7                # alpha_RM
    vi_bone      <- vv_bone      * 9.3                # alpha_BONE
    vi_heart     <- vv_heart     * 3.7                # alpha_HRT
    vi_rest      <- vv_rest      * 3.7                # alpha_REST
    # Two thirds of the remaining renal volume is proximal tubular cells
    # (Table S1 footnote c).
    vintra_kidney <- (vtot_kidney - vi_kidney - vv_kidney) * 2 / 3

    # ---- Serum flows (L/min), Table S1 -----------------------------------
    fq_liver     <- 0.065 * fserum      # hepatic arterial flow only
    fq_spleen    <- 0.03  * fserum
    fq_kidney    <- 0.19  * fserum
    fq_muscle    <- 0.17  * fserum
    fq_gi        <- 0.16  * fserum
    fq_skin      <- 0.05  * fserum
    fq_adipose   <- 0.05  * fserum
    fq_redmarrow <- 0.03  * fserum
    fq_bone      <- 0.05  * fserum
    fq_heart     <- 0.04  * fserum
    fq_brain     <- 0.12  * fserum
    fq_reprod    <- fperf_reprod * (1 - hct) * vtot_reprod
    fq_adrenals  <- fperf_adren  * (1 - hct) * vtot_adrenals
    fq_tumor     <- fperf_tumor  * (1 - hct) * vtot_tumor
    # F_REST = F - sum over all regions except the tumour.
    fq_rest <- fserum - fq_liver - fq_spleen - fq_kidney - fq_muscle - fq_gi -
      fq_skin - fq_adipose - fq_redmarrow - fq_bone - fq_heart - fq_brain -
      fq_reprod - fq_adrenals
    # F_TOTAL = F + F_TU (Table S1): the tumour is perfused in addition to the
    # systemic flow, so the lungs, arteries and veins all carry F_TOTAL. Table
    # S1 also prints F_LU = F, which would leave the tumour flow unbalanced;
    # F_TOTAL is used because mass balance requires it.
    fq_total <- fserum + fq_tumor

    # ---- Permeability-surface area products (L/min), Table S1 ------------
    ps_tumor     <- kps_tumor    * vtot_tumor
    ps_liver     <- kps_liver    * vtot_liver
    ps_spleen    <- kps_liver    * vtot_spleen      # k_S = k_L, similar capillary structure
    ps_lungs     <- kps_liver    * vtot_lungs       # k_LU = k_MUS * 100
    ps_adrenals  <- kps_liver    * vtot_adrenals    # k_AD = k_MUS * 100
    ps_reprod    <- kps_reprod_d * vtot_reprod
    ps_muscle    <- kps_muscle   * vtot_muscle
    ps_gi        <- kps_muscle   * vtot_gi
    ps_skin      <- kps_muscle   * vtot_skin
    ps_adipose   <- kps_muscle   * vtot_adipose
    ps_redmarrow <- kps_muscle   * vtot_redmarrow
    ps_bone      <- kps_muscle   * vtot_bone
    ps_heart     <- kps_muscle   * vtot_heart
    ps_rest      <- kps_muscle   * vtot_rest
    # The kidneys exchange by glomerular filtration rather than by a
    # permeability-surface product, and the brain has PS = 0 (equation S4).

    # ---- Total sst2 receptor numbers (nmol), Table S1 --------------------
    r0_kidney    <- rdens_kidney * vtot_kidney
    r0_liver     <- rdens_liver  * vtot_liver
    r0_spleen    <- rdens_spleen * vtot_spleen
    r0_tumor     <- rdens_tumor  * vtot_tumor
    r0_rest      <- rdens_rest   * vtot_rest
    r0_reprod    <- rdens_kidney * rratio_reprod * vtot_reprod
    r0_adrenals  <- rdens_kidney * rratio_adren  * vtot_adrenals
    r0_muscle    <- rdens_kidney * rratio_muscle * vtot_muscle
    r0_gi        <- rdens_kidney * rratio_gi     * vtot_gi
    r0_redmarrow <- rdens_kidney * rratio_rm     * vtot_redmarrow

    # ---- Free receptors, equation S1 -------------------------------------
    # The receptor pool is conserved and is shared by labelled and unlabelled
    # peptide, so the free number is the total less the TOTAL bound.
    rfree_kidney    <- r0_kidney    - rp_kidney_tot
    rfree_liver     <- r0_liver     - rp_liver_tot
    rfree_spleen    <- r0_spleen    - rp_spleen_tot
    rfree_tumor     <- r0_tumor     - rp_tumor_tot
    rfree_rest      <- r0_rest      - rp_rest_tot
    rfree_reprod    <- r0_reprod    - rp_reprod_tot
    rfree_adrenals  <- r0_adrenals  - rp_adrenals_tot
    rfree_muscle    <- r0_muscle    - rp_muscle_tot
    rfree_gi        <- r0_gi        - rp_gi_tot
    rfree_redmarrow <- r0_redmarrow - rp_redmarrow_tot

    # =====================================================================
    # 3. The two circulations. The paper writes one system for unlabelled and
    #    one for labelled peptide (equations S2-S13), coupled by physical
    #    decay, which converts labelled into unlabelled peptide, and by the
    #    shared receptor pool. Because labelled and unlabelled peptide are
    #    assumed kinetically identical, that pair is written here in the
    #    exactly equivalent TOTAL / LABELLED basis: summing the two printed
    #    equations cancels every decay term, so the total-peptide states obey
    #    the same equations with no decay, while the labelled states keep
    #    -lambda_phy and lose their decay-gain term. The unlabelled states of
    #    the paper are recovered as P_unl = P_tot - P_lab * n_inj_labelled.
    #    The change of basis is exact, and it removes a scale separation of
    #    about 1e3 between the two printed systems: 140 MBq of 111In is only
    #    0.082 nmol of the 75 nmol of DOTATATE injected (140 MBq over the
    #    1716 MBq/nmol molar activity of carrier-free 111In). The total states
    #    carry nmol; the labelled states are dosed with 1 and therefore carry
    #    the fraction of injected activity directly. Figure S1 plots the same
    #    curves as a labelled peptide amount, which is this fraction times the
    #    0.082 nmol of labelled peptide injected.
    # =====================================================================
    # ---- total peptide (nmol) --------------------------------------

    # Free peptide, vascular. Equation S4: transcapillary exchange plus
    # arterial delivery and venous drainage.
    d/dt(pv_tumor_tot) <- ps_tumor * (pint_tumor_tot / vi_tumor - pv_tumor_tot / vv_tumor) +
      fq_tumor * (part_tot / vv_art - pv_tumor_tot / vv_tumor)
    d/dt(pv_liver_tot) <- ps_liver * (pint_liver_tot / vi_liver - pv_liver_tot / vv_liver) +
      fq_liver * part_tot / vv_art + fq_spleen * pv_spleen_tot / vv_spleen +
      fq_gi * pv_gi_tot / vv_gi -
      (fq_liver + fq_spleen + fq_gi) * pv_liver_tot / vv_liver
    d/dt(pv_spleen_tot) <- ps_spleen * (pint_spleen_tot / vi_spleen - pv_spleen_tot / vv_spleen) +
      fq_spleen * (part_tot / vv_art - pv_spleen_tot / vv_spleen)
    d/dt(pv_gi_tot) <- ps_gi * (pint_gi_tot / vi_gi - pv_gi_tot / vv_gi) +
      fq_gi * (part_tot / vv_art - pv_gi_tot / vv_gi)
    d/dt(pv_muscle_tot) <- ps_muscle * (pint_muscle_tot / vi_muscle - pv_muscle_tot / vv_muscle) +
      fq_muscle * (part_tot / vv_art - pv_muscle_tot / vv_muscle)
    d/dt(pv_redmarrow_tot) <- ps_redmarrow * (pint_redmarrow_tot / vi_redmarrow - pv_redmarrow_tot / vv_redmarrow) +
      fq_redmarrow * (part_tot / vv_art - pv_redmarrow_tot / vv_redmarrow)
    d/dt(pv_reprod_tot) <- ps_reprod * (pint_reprod_tot / vi_reprod - pv_reprod_tot / vv_reprod) +
      fq_reprod * (part_tot / vv_art - pv_reprod_tot / vv_reprod)
    d/dt(pv_adrenals_tot) <- ps_adrenals * (pint_adrenals_tot / vi_adrenals - pv_adrenals_tot / vv_adrenals) +
      fq_adrenals * (part_tot / vv_art - pv_adrenals_tot / vv_adrenals)
    d/dt(pv_rest_tot) <- ps_rest * (pint_rest_tot / vi_rest - pv_rest_tot / vv_rest) +
      fq_rest * (part_tot / vv_art - pv_rest_tot / vv_rest)
    d/dt(pv_lungs_tot) <- ps_lungs * (pint_lungs_tot / vi_lungs - pv_lungs_tot / vv_lungs) +
      fq_total * (pven_tot / vv_ven - pv_lungs_tot / vv_lungs)
    d/dt(pv_skin_tot) <- ps_skin * (pint_skin_tot / vi_skin - pv_skin_tot / vv_skin) +
      fq_skin * (part_tot / vv_art - pv_skin_tot / vv_skin)
    d/dt(pv_adipose_tot) <- ps_adipose * (pint_adipose_tot / vi_adipose - pv_adipose_tot / vv_adipose) +
      fq_adipose * (part_tot / vv_art - pv_adipose_tot / vv_adipose)
    d/dt(pv_heart_tot) <- ps_heart * (pint_heart_tot / vi_heart - pv_heart_tot / vv_heart) +
      fq_heart * (part_tot / vv_art - pv_heart_tot / vv_heart)
    d/dt(pv_bone_tot) <- ps_bone * (pint_bone_tot / vi_bone - pv_bone_tot / vv_bone) +
      fq_bone * (part_tot / vv_art - pv_bone_tot / vv_bone)
    d/dt(pv_kidney_tot) <- - pv_kidney_tot / vv_kidney * (f_fil + fq_kidney) +
      fq_kidney * part_tot / vv_art +
      f_reab * pintra_kidney_tot / vintra_kidney
    d/dt(pv_brain_tot) <- fq_brain * (part_tot / vv_art - pv_brain_tot / vv_brain)

    # Free peptide, interstitial. Equation S11 where sst2 is expressed,
    # equation S10 where it is not, equation S9 for the kidneys.
    d/dt(pint_tumor_tot) <- - kon * pint_tumor_tot * rfree_tumor / vi_tumor + koff * rp_tumor_tot +
      ps_tumor * (pv_tumor_tot / vv_tumor - pint_tumor_tot / vi_tumor)
    d/dt(pint_liver_tot) <- - kon * pint_liver_tot * rfree_liver / vi_liver + koff * rp_liver_tot +
      ps_liver * (pv_liver_tot / vv_liver - pint_liver_tot / vi_liver)
    d/dt(pint_spleen_tot) <- - kon * pint_spleen_tot * rfree_spleen / vi_spleen + koff * rp_spleen_tot +
      ps_spleen * (pv_spleen_tot / vv_spleen - pint_spleen_tot / vi_spleen)
    d/dt(pint_gi_tot) <- - kon * pint_gi_tot * rfree_gi / vi_gi + koff * rp_gi_tot +
      ps_gi * (pv_gi_tot / vv_gi - pint_gi_tot / vi_gi)
    d/dt(pint_muscle_tot) <- - kon * pint_muscle_tot * rfree_muscle / vi_muscle + koff * rp_muscle_tot +
      ps_muscle * (pv_muscle_tot / vv_muscle - pint_muscle_tot / vi_muscle)
    d/dt(pint_redmarrow_tot) <- - kon * pint_redmarrow_tot * rfree_redmarrow / vi_redmarrow + koff * rp_redmarrow_tot +
      ps_redmarrow * (pv_redmarrow_tot / vv_redmarrow - pint_redmarrow_tot / vi_redmarrow)
    d/dt(pint_reprod_tot) <- - kon * pint_reprod_tot * rfree_reprod / vi_reprod + koff * rp_reprod_tot +
      ps_reprod * (pv_reprod_tot / vv_reprod - pint_reprod_tot / vi_reprod)
    d/dt(pint_adrenals_tot) <- - kon * pint_adrenals_tot * rfree_adrenals / vi_adrenals + koff * rp_adrenals_tot +
      ps_adrenals * (pv_adrenals_tot / vv_adrenals - pint_adrenals_tot / vi_adrenals)
    d/dt(pint_rest_tot) <- - kon * pint_rest_tot * rfree_rest / vi_rest + koff * rp_rest_tot +
      ps_rest * (pv_rest_tot / vv_rest - pint_rest_tot / vi_rest)
    d/dt(pint_kidney_tot) <- - kon * pint_kidney_tot * rfree_kidney / vi_kidney + koff * rp_kidney_tot +
      f_fil * (pv_kidney_tot / vv_kidney - pint_kidney_tot / vi_kidney)
    d/dt(pint_lungs_tot) <- ps_lungs * (pv_lungs_tot / vv_lungs - pint_lungs_tot / vi_lungs)
    d/dt(pint_skin_tot) <- ps_skin * (pv_skin_tot / vv_skin - pint_skin_tot / vi_skin)
    d/dt(pint_adipose_tot) <- ps_adipose * (pv_adipose_tot / vv_adipose - pint_adipose_tot / vi_adipose)
    d/dt(pint_heart_tot) <- ps_heart * (pv_heart_tot / vv_heart - pint_heart_tot / vi_heart)
    d/dt(pint_bone_tot) <- ps_bone * (pv_bone_tot / vv_bone - pint_bone_tot / vi_bone)

    # Peptide bound to sst2 on the cell surface, equation S3.
    d/dt(rp_tumor_tot) <- kon * pint_tumor_tot * rfree_tumor / vi_tumor -
      (koff + kint_tumor) * rp_tumor_tot
    d/dt(rp_liver_tot) <- kon * pint_liver_tot * rfree_liver / vi_liver -
      (koff + kint_nt) * rp_liver_tot
    d/dt(rp_spleen_tot) <- kon * pint_spleen_tot * rfree_spleen / vi_spleen -
      (koff + kint_nt) * rp_spleen_tot
    d/dt(rp_gi_tot) <- kon * pint_gi_tot * rfree_gi / vi_gi -
      (koff + kint_nt) * rp_gi_tot
    d/dt(rp_muscle_tot) <- kon * pint_muscle_tot * rfree_muscle / vi_muscle -
      (koff + kint_nt) * rp_muscle_tot
    d/dt(rp_redmarrow_tot) <- kon * pint_redmarrow_tot * rfree_redmarrow / vi_redmarrow -
      (koff + kint_nt) * rp_redmarrow_tot
    d/dt(rp_reprod_tot) <- kon * pint_reprod_tot * rfree_reprod / vi_reprod -
      (koff + kint_nt) * rp_reprod_tot
    d/dt(rp_adrenals_tot) <- kon * pint_adrenals_tot * rfree_adrenals / vi_adrenals -
      (koff + kint_nt) * rp_adrenals_tot
    d/dt(rp_rest_tot) <- kon * pint_rest_tot * rfree_rest / vi_rest -
      (koff + kint_nt) * rp_rest_tot
    d/dt(rp_kidney_tot) <- kon * pint_kidney_tot * rfree_kidney / vi_kidney -
      (koff + kint_nt) * rp_kidney_tot

    # Internalised peptide, equation S2. Degraded peptide and released free
    # radionuclide are assumed cleared from the body directly.
    d/dt(pintern_tumor_tot) <- kint_tumor * rp_tumor_tot - kdeg_tumor * pintern_tumor_tot
    d/dt(pintern_liver_tot) <- kint_nt * rp_liver_tot - kdeg_nt * pintern_liver_tot
    d/dt(pintern_spleen_tot) <- kint_nt * rp_spleen_tot - kdeg_nt * pintern_spleen_tot
    d/dt(pintern_gi_tot) <- kint_nt * rp_gi_tot - kdeg_nt * pintern_gi_tot
    d/dt(pintern_muscle_tot) <- kint_nt * rp_muscle_tot - kdeg_nt * pintern_muscle_tot
    d/dt(pintern_redmarrow_tot) <- kint_nt * rp_redmarrow_tot - kdeg_nt * pintern_redmarrow_tot
    d/dt(pintern_reprod_tot) <- kint_nt * rp_reprod_tot - kdeg_nt * pintern_reprod_tot
    d/dt(pintern_adrenals_tot) <- kint_nt * rp_adrenals_tot - kdeg_nt * pintern_adrenals_tot
    d/dt(pintern_rest_tot) <- kint_nt * rp_rest_tot - kdeg_nt * pintern_rest_tot
    d/dt(pintern_kidney_tot) <- kint_nt * rp_kidney_tot - kdeg_nt * pintern_kidney_tot

    # Unspecific uptake into the proximal tubular cells, equation S12.
    d/dt(pintra_kidney_tot) <- f_reab * (pint_kidney_tot / vi_kidney - pintra_kidney_tot / vintra_kidney)

    # Arteries (S8), veins (S7) and the serum-protein-bound pool (S13).
    d/dt(part_tot) <- fq_total * (pv_lungs_tot / vv_lungs - part_tot / vv_art)
    d/dt(pven_tot) <- - kpr * pven_tot - fq_total * pven_tot / vv_ven +
      fq_tumor * pv_tumor_tot / vv_tumor + fq_muscle * pv_muscle_tot / vv_muscle +
      fq_redmarrow * pv_redmarrow_tot / vv_redmarrow + fq_reprod * pv_reprod_tot / vv_reprod +
      fq_adrenals * pv_adrenals_tot / vv_adrenals + fq_rest * pv_rest_tot / vv_rest +
      fq_kidney * pv_kidney_tot / vv_kidney + fq_skin * pv_skin_tot / vv_skin +
      fq_adipose * pv_adipose_tot / vv_adipose + fq_heart * pv_heart_tot / vv_heart +
      fq_bone * pv_bone_tot / vv_bone + (fq_liver + fq_spleen + fq_gi) * pv_liver_tot / vv_liver +
      fq_brain * pv_brain_tot / vv_brain
    d/dt(pprp_tot) <- kpr * pven_tot

    # Bookkeeping states (not in the paper): cumulative urinary excretion and
    # cumulative loss by degradation, so that mass balance can be checked.
    d/dt(urine_tot) <- f_exc * pint_kidney_tot / vi_kidney
    d/dt(cleared_tot) <-
      kdeg_tumor * pintern_tumor_tot + kdeg_nt * pintern_liver_tot +
      kdeg_nt * pintern_spleen_tot + kdeg_nt * pintern_gi_tot +
      kdeg_nt * pintern_muscle_tot + kdeg_nt * pintern_redmarrow_tot +
      kdeg_nt * pintern_reprod_tot + kdeg_nt * pintern_adrenals_tot +
      kdeg_nt * pintern_rest_tot + kdeg_nt * pintern_kidney_tot

    # ---- labelled peptide (fraction of injected activity) ----------

    # Free peptide, vascular. Equation S4: transcapillary exchange plus
    # arterial delivery and venous drainage.
    d/dt(pv_tumor_lab) <- ps_tumor * (pint_tumor_lab / vi_tumor - pv_tumor_lab / vv_tumor) +
      fq_tumor * (part_lab / vv_art - pv_tumor_lab / vv_tumor) - lambda_phy * pv_tumor_lab
    d/dt(pv_liver_lab) <- ps_liver * (pint_liver_lab / vi_liver - pv_liver_lab / vv_liver) +
      fq_liver * part_lab / vv_art + fq_spleen * pv_spleen_lab / vv_spleen +
      fq_gi * pv_gi_lab / vv_gi -
      (fq_liver + fq_spleen + fq_gi) * pv_liver_lab / vv_liver - lambda_phy * pv_liver_lab
    d/dt(pv_spleen_lab) <- ps_spleen * (pint_spleen_lab / vi_spleen - pv_spleen_lab / vv_spleen) +
      fq_spleen * (part_lab / vv_art - pv_spleen_lab / vv_spleen) - lambda_phy * pv_spleen_lab
    d/dt(pv_gi_lab) <- ps_gi * (pint_gi_lab / vi_gi - pv_gi_lab / vv_gi) +
      fq_gi * (part_lab / vv_art - pv_gi_lab / vv_gi) - lambda_phy * pv_gi_lab
    d/dt(pv_muscle_lab) <- ps_muscle * (pint_muscle_lab / vi_muscle - pv_muscle_lab / vv_muscle) +
      fq_muscle * (part_lab / vv_art - pv_muscle_lab / vv_muscle) - lambda_phy * pv_muscle_lab
    d/dt(pv_redmarrow_lab) <- ps_redmarrow * (pint_redmarrow_lab / vi_redmarrow - pv_redmarrow_lab / vv_redmarrow) +
      fq_redmarrow * (part_lab / vv_art - pv_redmarrow_lab / vv_redmarrow) - lambda_phy * pv_redmarrow_lab
    d/dt(pv_reprod_lab) <- ps_reprod * (pint_reprod_lab / vi_reprod - pv_reprod_lab / vv_reprod) +
      fq_reprod * (part_lab / vv_art - pv_reprod_lab / vv_reprod) - lambda_phy * pv_reprod_lab
    d/dt(pv_adrenals_lab) <- ps_adrenals * (pint_adrenals_lab / vi_adrenals - pv_adrenals_lab / vv_adrenals) +
      fq_adrenals * (part_lab / vv_art - pv_adrenals_lab / vv_adrenals) - lambda_phy * pv_adrenals_lab
    d/dt(pv_rest_lab) <- ps_rest * (pint_rest_lab / vi_rest - pv_rest_lab / vv_rest) +
      fq_rest * (part_lab / vv_art - pv_rest_lab / vv_rest) - lambda_phy * pv_rest_lab
    d/dt(pv_lungs_lab) <- ps_lungs * (pint_lungs_lab / vi_lungs - pv_lungs_lab / vv_lungs) +
      fq_total * (pven_lab / vv_ven - pv_lungs_lab / vv_lungs) - lambda_phy * pv_lungs_lab
    d/dt(pv_skin_lab) <- ps_skin * (pint_skin_lab / vi_skin - pv_skin_lab / vv_skin) +
      fq_skin * (part_lab / vv_art - pv_skin_lab / vv_skin) - lambda_phy * pv_skin_lab
    d/dt(pv_adipose_lab) <- ps_adipose * (pint_adipose_lab / vi_adipose - pv_adipose_lab / vv_adipose) +
      fq_adipose * (part_lab / vv_art - pv_adipose_lab / vv_adipose) - lambda_phy * pv_adipose_lab
    d/dt(pv_heart_lab) <- ps_heart * (pint_heart_lab / vi_heart - pv_heart_lab / vv_heart) +
      fq_heart * (part_lab / vv_art - pv_heart_lab / vv_heart) - lambda_phy * pv_heart_lab
    d/dt(pv_bone_lab) <- ps_bone * (pint_bone_lab / vi_bone - pv_bone_lab / vv_bone) +
      fq_bone * (part_lab / vv_art - pv_bone_lab / vv_bone) - lambda_phy * pv_bone_lab
    d/dt(pv_kidney_lab) <- - pv_kidney_lab / vv_kidney * (f_fil + fq_kidney) +
      fq_kidney * part_lab / vv_art +
      f_reab * pintra_kidney_lab / vintra_kidney - lambda_phy * pv_kidney_lab
    d/dt(pv_brain_lab) <- fq_brain * (part_lab / vv_art - pv_brain_lab / vv_brain) - lambda_phy * pv_brain_lab

    # Free peptide, interstitial. Equation S11 where sst2 is expressed,
    # equation S10 where it is not, equation S9 for the kidneys.
    d/dt(pint_tumor_lab) <- - kon * pint_tumor_lab * rfree_tumor / vi_tumor + koff * rp_tumor_lab +
      ps_tumor * (pv_tumor_lab / vv_tumor - pint_tumor_lab / vi_tumor) - lambda_phy * pint_tumor_lab
    d/dt(pint_liver_lab) <- - kon * pint_liver_lab * rfree_liver / vi_liver + koff * rp_liver_lab +
      ps_liver * (pv_liver_lab / vv_liver - pint_liver_lab / vi_liver) - lambda_phy * pint_liver_lab
    d/dt(pint_spleen_lab) <- - kon * pint_spleen_lab * rfree_spleen / vi_spleen + koff * rp_spleen_lab +
      ps_spleen * (pv_spleen_lab / vv_spleen - pint_spleen_lab / vi_spleen) - lambda_phy * pint_spleen_lab
    d/dt(pint_gi_lab) <- - kon * pint_gi_lab * rfree_gi / vi_gi + koff * rp_gi_lab +
      ps_gi * (pv_gi_lab / vv_gi - pint_gi_lab / vi_gi) - lambda_phy * pint_gi_lab
    d/dt(pint_muscle_lab) <- - kon * pint_muscle_lab * rfree_muscle / vi_muscle + koff * rp_muscle_lab +
      ps_muscle * (pv_muscle_lab / vv_muscle - pint_muscle_lab / vi_muscle) - lambda_phy * pint_muscle_lab
    d/dt(pint_redmarrow_lab) <- - kon * pint_redmarrow_lab * rfree_redmarrow / vi_redmarrow + koff * rp_redmarrow_lab +
      ps_redmarrow * (pv_redmarrow_lab / vv_redmarrow - pint_redmarrow_lab / vi_redmarrow) - lambda_phy * pint_redmarrow_lab
    d/dt(pint_reprod_lab) <- - kon * pint_reprod_lab * rfree_reprod / vi_reprod + koff * rp_reprod_lab +
      ps_reprod * (pv_reprod_lab / vv_reprod - pint_reprod_lab / vi_reprod) - lambda_phy * pint_reprod_lab
    d/dt(pint_adrenals_lab) <- - kon * pint_adrenals_lab * rfree_adrenals / vi_adrenals + koff * rp_adrenals_lab +
      ps_adrenals * (pv_adrenals_lab / vv_adrenals - pint_adrenals_lab / vi_adrenals) - lambda_phy * pint_adrenals_lab
    d/dt(pint_rest_lab) <- - kon * pint_rest_lab * rfree_rest / vi_rest + koff * rp_rest_lab +
      ps_rest * (pv_rest_lab / vv_rest - pint_rest_lab / vi_rest) - lambda_phy * pint_rest_lab
    d/dt(pint_kidney_lab) <- - kon * pint_kidney_lab * rfree_kidney / vi_kidney + koff * rp_kidney_lab +
      f_fil * (pv_kidney_lab / vv_kidney - pint_kidney_lab / vi_kidney) - lambda_phy * pint_kidney_lab
    d/dt(pint_lungs_lab) <- ps_lungs * (pv_lungs_lab / vv_lungs - pint_lungs_lab / vi_lungs) - lambda_phy * pint_lungs_lab
    d/dt(pint_skin_lab) <- ps_skin * (pv_skin_lab / vv_skin - pint_skin_lab / vi_skin) - lambda_phy * pint_skin_lab
    d/dt(pint_adipose_lab) <- ps_adipose * (pv_adipose_lab / vv_adipose - pint_adipose_lab / vi_adipose) - lambda_phy * pint_adipose_lab
    d/dt(pint_heart_lab) <- ps_heart * (pv_heart_lab / vv_heart - pint_heart_lab / vi_heart) - lambda_phy * pint_heart_lab
    d/dt(pint_bone_lab) <- ps_bone * (pv_bone_lab / vv_bone - pint_bone_lab / vi_bone) - lambda_phy * pint_bone_lab

    # Peptide bound to sst2 on the cell surface, equation S3.
    d/dt(rp_tumor_lab) <- kon * pint_tumor_lab * rfree_tumor / vi_tumor -
      (koff + kint_tumor) * rp_tumor_lab - lambda_phy * rp_tumor_lab
    d/dt(rp_liver_lab) <- kon * pint_liver_lab * rfree_liver / vi_liver -
      (koff + kint_nt) * rp_liver_lab - lambda_phy * rp_liver_lab
    d/dt(rp_spleen_lab) <- kon * pint_spleen_lab * rfree_spleen / vi_spleen -
      (koff + kint_nt) * rp_spleen_lab - lambda_phy * rp_spleen_lab
    d/dt(rp_gi_lab) <- kon * pint_gi_lab * rfree_gi / vi_gi -
      (koff + kint_nt) * rp_gi_lab - lambda_phy * rp_gi_lab
    d/dt(rp_muscle_lab) <- kon * pint_muscle_lab * rfree_muscle / vi_muscle -
      (koff + kint_nt) * rp_muscle_lab - lambda_phy * rp_muscle_lab
    d/dt(rp_redmarrow_lab) <- kon * pint_redmarrow_lab * rfree_redmarrow / vi_redmarrow -
      (koff + kint_nt) * rp_redmarrow_lab - lambda_phy * rp_redmarrow_lab
    d/dt(rp_reprod_lab) <- kon * pint_reprod_lab * rfree_reprod / vi_reprod -
      (koff + kint_nt) * rp_reprod_lab - lambda_phy * rp_reprod_lab
    d/dt(rp_adrenals_lab) <- kon * pint_adrenals_lab * rfree_adrenals / vi_adrenals -
      (koff + kint_nt) * rp_adrenals_lab - lambda_phy * rp_adrenals_lab
    d/dt(rp_rest_lab) <- kon * pint_rest_lab * rfree_rest / vi_rest -
      (koff + kint_nt) * rp_rest_lab - lambda_phy * rp_rest_lab
    d/dt(rp_kidney_lab) <- kon * pint_kidney_lab * rfree_kidney / vi_kidney -
      (koff + kint_nt) * rp_kidney_lab - lambda_phy * rp_kidney_lab

    # Internalised peptide, equation S2. Degraded peptide and released free
    # radionuclide are assumed cleared from the body directly.
    d/dt(pintern_tumor_lab) <- kint_tumor * rp_tumor_lab - kdeg_tumor * pintern_tumor_lab - lambda_phy * pintern_tumor_lab
    d/dt(pintern_liver_lab) <- kint_nt * rp_liver_lab - kdeg_nt * pintern_liver_lab - lambda_phy * pintern_liver_lab
    d/dt(pintern_spleen_lab) <- kint_nt * rp_spleen_lab - kdeg_nt * pintern_spleen_lab - lambda_phy * pintern_spleen_lab
    d/dt(pintern_gi_lab) <- kint_nt * rp_gi_lab - kdeg_nt * pintern_gi_lab - lambda_phy * pintern_gi_lab
    d/dt(pintern_muscle_lab) <- kint_nt * rp_muscle_lab - kdeg_nt * pintern_muscle_lab - lambda_phy * pintern_muscle_lab
    d/dt(pintern_redmarrow_lab) <- kint_nt * rp_redmarrow_lab - kdeg_nt * pintern_redmarrow_lab - lambda_phy * pintern_redmarrow_lab
    d/dt(pintern_reprod_lab) <- kint_nt * rp_reprod_lab - kdeg_nt * pintern_reprod_lab - lambda_phy * pintern_reprod_lab
    d/dt(pintern_adrenals_lab) <- kint_nt * rp_adrenals_lab - kdeg_nt * pintern_adrenals_lab - lambda_phy * pintern_adrenals_lab
    d/dt(pintern_rest_lab) <- kint_nt * rp_rest_lab - kdeg_nt * pintern_rest_lab - lambda_phy * pintern_rest_lab
    d/dt(pintern_kidney_lab) <- kint_nt * rp_kidney_lab - kdeg_nt * pintern_kidney_lab - lambda_phy * pintern_kidney_lab

    # Unspecific uptake into the proximal tubular cells, equation S12.
    d/dt(pintra_kidney_lab) <- f_reab * (pint_kidney_lab / vi_kidney - pintra_kidney_lab / vintra_kidney) - lambda_phy * pintra_kidney_lab

    # Arteries (S8), veins (S7) and the serum-protein-bound pool (S13).
    d/dt(part_lab) <- fq_total * (pv_lungs_lab / vv_lungs - part_lab / vv_art) - lambda_phy * part_lab
    d/dt(pven_lab) <- - kpr * pven_lab - fq_total * pven_lab / vv_ven +
      fq_tumor * pv_tumor_lab / vv_tumor + fq_muscle * pv_muscle_lab / vv_muscle +
      fq_redmarrow * pv_redmarrow_lab / vv_redmarrow + fq_reprod * pv_reprod_lab / vv_reprod +
      fq_adrenals * pv_adrenals_lab / vv_adrenals + fq_rest * pv_rest_lab / vv_rest +
      fq_kidney * pv_kidney_lab / vv_kidney + fq_skin * pv_skin_lab / vv_skin +
      fq_adipose * pv_adipose_lab / vv_adipose + fq_heart * pv_heart_lab / vv_heart +
      fq_bone * pv_bone_lab / vv_bone + (fq_liver + fq_spleen + fq_gi) * pv_liver_lab / vv_liver +
      fq_brain * pv_brain_lab / vv_brain - lambda_phy * pven_lab
    d/dt(pprp_lab) <- kpr * pven_lab - lambda_phy * pprp_lab

    # Bookkeeping states (not in the paper): cumulative urinary excretion and
    # cumulative loss by degradation, so that mass balance can be checked.
    d/dt(urine_lab) <- f_exc * pint_kidney_lab / vi_kidney - lambda_phy * urine_lab
    d/dt(cleared_lab) <-
      kdeg_tumor * pintern_tumor_lab + kdeg_nt * pintern_liver_lab +
      kdeg_nt * pintern_spleen_lab + kdeg_nt * pintern_gi_lab +
      kdeg_nt * pintern_muscle_lab + kdeg_nt * pintern_redmarrow_lab +
      kdeg_nt * pintern_reprod_lab + kdeg_nt * pintern_adrenals_lab +
      kdeg_nt * pintern_rest_lab + kdeg_nt * pintern_kidney_lab - lambda_phy * cleared_lab

    # =====================================================================
    # 4. Observed regions. The paper fitted six time-activity data sets
    #    (Table 1 residual errors): blood serum, whole body, kidneys,
    #    spleen, liver and tumour. An organ region of interest contains
    #    every sub-compartment of that organ; the upstream Kletting 2016
    #    supplement states this explicitly for the red marrow, whose time-
    #    integrated activity coefficient was obtained by 'integrating all
    #    compartments describing the red marrow including red marrow serum,
    #    interstitial space, bound and internalized peptide'. Because the
    #    labelled states are dosed with 1, each output is already a
    #    fraction of the injected activity.
    # =====================================================================
    Akidney <- pv_kidney_lab + pint_kidney_lab + rp_kidney_lab + pintern_kidney_lab + pintra_kidney_lab
    Aliver <- pv_liver_lab + pint_liver_lab + rp_liver_lab + pintern_liver_lab
    Aspleen <- pv_spleen_lab + pint_spleen_lab + rp_spleen_lab + pintern_spleen_lab
    Atumor <- pv_tumor_lab + pint_tumor_lab + rp_tumor_lab + pintern_tumor_lab

    # Serum: an activity CONCENTRATION, not an amount. The measured datum is
    # the activity concentration of a serum aliquot, so the prediction is the
    # labelled peptide dissolved anywhere in serum -- the arterial, venous and
    # every regional vascular sub-volume, plus the irreversibly
    # serum-protein-bound pool of equation S13 -- divided by the total body
    # serum volume V_P. Its units are therefore fraction of injected activity
    # per litre, so its time integral is min/L rather than the min that
    # Table S3 prints for all six regions alike. Dividing by V_P is what
    # reproduces the published serum rTIACs (see the vignette).
    Aserum_amt <-
      part_lab + pven_lab + pprp_lab +
      pv_tumor_lab + pv_liver_lab + pv_spleen_lab +
      pv_gi_lab + pv_muscle_lab + pv_redmarrow_lab +
      pv_reprod_lab + pv_adrenals_lab + pv_rest_lab +
      pv_kidney_lab + pv_lungs_lab + pv_skin_lab +
      pv_adipose_lab + pv_heart_lab + pv_bone_lab +
      pv_brain_lab
    Aserum <- Aserum_amt / vserum

    # Whole body: everything still in the body, i.e. every labelled state
    # except the two cumulative-loss bookkeeping states.
    Abody <-
      part_lab + pven_lab + pprp_lab +
      pintra_kidney_lab + pv_tumor_lab + pv_liver_lab +
      pv_spleen_lab + pv_gi_lab + pv_muscle_lab +
      pv_redmarrow_lab + pv_reprod_lab + pv_adrenals_lab +
      pv_rest_lab + pv_kidney_lab + pv_lungs_lab +
      pv_skin_lab + pv_adipose_lab + pv_heart_lab +
      pv_bone_lab + pv_brain_lab + pint_tumor_lab +
      pint_liver_lab + pint_spleen_lab + pint_gi_lab +
      pint_muscle_lab + pint_redmarrow_lab + pint_reprod_lab +
      pint_adrenals_lab + pint_rest_lab + pint_kidney_lab +
      pint_lungs_lab + pint_skin_lab + pint_adipose_lab +
      pint_heart_lab + pint_bone_lab + rp_tumor_lab +
      rp_liver_lab + rp_spleen_lab + rp_gi_lab +
      rp_muscle_lab + rp_redmarrow_lab + rp_reprod_lab +
      rp_adrenals_lab + rp_rest_lab + rp_kidney_lab +
      pintern_tumor_lab + pintern_liver_lab + pintern_spleen_lab +
      pintern_gi_lab + pintern_muscle_lab + pintern_redmarrow_lab +
      pintern_reprod_lab + pintern_adrenals_lab + pintern_rest_lab +
      pintern_kidney_lab

    # Total peptide still in the body, used for mass-balance checking and
    # for receptor-occupancy diagnostics.
    occ_tumor  <- rp_tumor_tot / r0_tumor
    occ_kidney <- rp_kidney_tot / r0_kidney

    # =====================================================================
    # 5. Residual error. Equation 2 is purely proportional, and setting II
    #    gives each observed region its own error parameter.
    # =====================================================================
    Aserum ~ prop(propSd_Aserum)
    Abody ~ prop(propSd_Abody)
    Akidney ~ prop(propSd_Akidney)
    Aspleen ~ prop(propSd_Aspleen)
    Aliver ~ prop(propSd_Aliver)
    Atumor ~ prop(propSd_Atumor)
  })
}
