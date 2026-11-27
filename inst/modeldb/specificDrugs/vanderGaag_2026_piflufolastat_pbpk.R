vanderGaag_2026_piflufolastat_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 74 ODEs, hand-coded deSolve). Physiologically based",
    "pharmacokinetic model of the PSMA PET tracer [18F]F-piflufolastat",
    "([18F]DCFPyL) in men with metastatic prostate cancer, rebuilt from the",
    "Begum 2018 [68Ga]Ga-PSMA-11 / [177Lu]Lu-PSMA-617 model. Eleven",
    "PSMA-expressing tissues (two delineated tumour lesions, a lumped",
    "remaining-tumour compartment, prostate, kidneys, liver, spleen, GI tract,",
    "parotid, submandibular and lacrimal glands) carry vascular, interstitial,",
    "receptor-bound, internalised and free-receptor states; in the kidney the",
    "interstitial state is the tubular lumen, fed by glomerular filtration",
    "(eGFR times a sieving coefficient) and drained by urinary excretion of a",
    "fixed fraction of the filtered load. Eight PSMA-negative tissues (lungs,",
    "red marrow, muscle, skin, adipose, heart, bone, rest) carry vascular and",
    "interstitial states, the brain is vascular only, and arteries and veins",
    "complete the circulation. Binding is second order onto a conserved",
    "receptor pool; internalised ligand is released irreversibly. All",
    "parameters are fixed literature or physiological values; only the tumour",
    "PSMA receptor density was fitted per patient to the 120 min PET image",
    "(paper Table S2), and the shipped tumour densities are the typical",
    "Table S1 values. Covariates are body weight, height, age, haematocrit,",
    "eGFR and the measured tumour volumes. Amounts are nmol, time is minutes;",
    "tissue outputs are decay-corrected activity concentrations in kBq/mL.",
    sep = " "
  )
  reference <- paste(
    "van der Gaag S, Yaqub M, Oprea-Lager DE, Bartelink IH, et al.",
    "Leveraging [18F]F-Piflufolastat PET for PBPK Prediction of Tissue",
    "Distribution and Therapy Optimization in Prostate Cancer.",
    "Clin Transl Sci. 2026;19:e70583. doi:10.1111/cts.70583.",
    "Drug-specific constants from main-text Table 1 ([18F]F-Piflufolastat",
    "column); equations from Supporting Information Part A (Suppl. Eq 1-15);",
    "physiological and remaining drug-specific parameters from Table S1;",
    "patient covariates from Table 2; fitted tumour receptor densities from",
    "Table S2. Where Table S1 is internally inconsistent or garbled, the",
    "authors' deposited code (github.com/ibartelink/Target-Site-PBPK-in-oncology)",
    "settles it; see the vignette Errata.",
    sep = " "
  )
  vignette <- "vanderGaag_2026_piflufolastat_pbpk"
  units <- list(time = "min", dosing = "nmol", concentration = "nmol/L")
  dosing <- "venous"

  paper_specific_compartment_pattern <- c(
    "^(vasc|int|bound|intern|rfree)_",
    "^(venous|arterial|lumen_kidney|released)$"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured per patient (Table 2, 64-119 kg). Enters the Du Bois body",
        "surface area (and hence serum volume and cardiac output), the",
        "BW/71-scaled organ volumes, the liver (0.023*BW) and spleen",
        "(0.003*BW) volumes, and the rest-of-body volume (Table S1)."
      ),
      source_name = "BW"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured per patient (Table 2). Enters only the Du Bois body surface area (Table S1).",
      source_name = "BH"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured per patient (Table 2). Sets the kidney serum-flow density",
        "F_k,age = 4.3 - 0.026*Age mL/min/g (Table S1)."
      ),
      source_name = "Age"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "% (volume fraction x 100)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters as the fraction H = HCT/100. Sets the serum volume",
        "2.8*(1-H)*BSA and hence total flow and every serum-scaled vascular",
        "volume, and the (1-H) factors on tumour, prostate and salivary-gland",
        "vascular volumes and on prostate and kidney flow (Table S1). Declared",
        "measured in Table S1 but not published per patient; the vignette uses",
        "the Table S4 prostate-cancer population mean of 43.3 %."
      ),
      source_name = "H"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, absolute (not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "eGFR in absolute mL/min (Table S1 lists eGFR in mL/min and Table S4",
        "gives the population mean 83 mL/min); NOT normalised to 1.73 m^2,",
        "as the register's CRCL aliases GFR / AGFR also allow. Filtration flow",
        "F_fil = eGFR * sieving coefficient, converted to L/min. Not published",
        "per patient; the vignette uses the Table S4 mean."
      ),
      source_name = "eGFR"
    ),
    ORGVOL_TUMOR1 = list(
      description = "Measured volume of delineated tumour lesion 1",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Measured per patient by PET/CT delineation (Table 2, cm^3 = mL). Must be > 0.",
      source_name = "Tumor lesion 1"
    ),
    ORGVOL_TUMOR2 = list(
      description = "Measured volume of delineated tumour lesion 2",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured per patient (Table 2). Must be > 0 (the vascular volume",
        "divides the ODEs); for a patient with a single lesion use a",
        "negligible volume such as 1e-6 mL."
      ),
      source_name = "Tumor lesion 2"
    ),
    ORGVOL_TUMORREST = list(
      description = "Measured volume of all remaining tumour lesions (lumped)",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured per patient (Table 2, 'Tumor REST'). Used directly as",
        "V_tot,turest; Table S1 lists a measured-to-actual ratio X_V = 0.64",
        "but no equation applies it, and neither does the authors' code.",
        "Must be > 0."
      ),
      source_name = "Tumor REST"
    )
  )

  # Every state holds an amount in nmol: free, bound or internalised ligand,
  # or (rfree_*) unoccupied PSMA receptors.
  .lig <- "[18F]F-piflufolastat"
  .psma <- c(
    "tumor1",
    "tumor2",
    "tumorrest",
    "prostate",
    "liver",
    "spleen",
    "gi",
    "parotid",
    "submandibular",
    "lacrimal",
    "kidney"
  )
  .npsma <- c("lungs", "redmarrow", "muscle", "skin", "adipose", "heart", "bone", "rest")
  .cmt <- c(
    "venous",
    "arterial",
    paste0("vasc_", c(.psma, .npsma, "brain")),
    paste0("int_", c(.psma[-11], .npsma)),
    "lumen_kidney",
    paste0("bound_", .psma),
    paste0("intern_", .psma),
    paste0("rfree_", .psma),
    "urine",
    "released"
  )
  .analyte <- ifelse(
    grepl("^rfree_", .cmt),
    "PSMA (free receptor)",
    ifelse(
      grepl("^bound_", .cmt),
      paste(.lig, "bound to PSMA"),
      ifelse(grepl("^intern_", .cmt), paste(.lig, "(internalised)"), .lig)
    )
  )
  .specimen <- ifelse(
    .cmt %in% c("venous", "arterial") | grepl("^vasc_", .cmt),
    "serum",
    ifelse(
      .cmt == "urine",
      "urine",
      ifelse(.cmt == "released", "not applicable", ifelse(grepl("tumor", .cmt), "tumor", "tissue"))
    )
  )
  compartmentData <- stats::setNames(
    Map(
      function(a, s) list(analyte = a, units = "nmol", specimen = s, verified = TRUE),
      .analyte,
      .specimen,
      USE.NAMES = FALSE
    ),
    .cmt
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 3L,
    age_range = "56-86 years",
    weight_range = "64-119 kg",
    height_range = "170-191 cm",
    sex_female_pct = 0,
    disease_state = "metastatic prostate cancer (PSMA-positive)",
    dose_range = paste(
      "245-319 MBq [18F]F-piflufolastat as an intravenous bolus (molar dose =",
      "activity / specific activity, 0.045 GBq/nmol); two patients also",
      "received a first therapeutic cycle of [177Lu]Lu-PSMA-617"
    ),
    regions = "Netherlands (Amsterdam UMC)",
    notes = paste(
      "Three verification cohorts (Table 2): internal, 8 patients with dynamic",
      "0-120 min PET (from Jansen et al.); external, 5 patients with a single",
      "total-body PET/CT at 120 min; therapeutic, 2 patients with diagnostic",
      "[18F]F-piflufolastat PET plus 24 h and 7 day [177Lu]Lu-PSMA-617 SPECT.",
      "This is a bottom-up PBPK: no parameter was estimated across the",
      "cohort, there is no between-subject variability term, and only the",
      "tumour receptor density was fitted per patient (Table S2)."
    )
  )

  ini({
    # ---- Drug-specific binding and renal handling (Table 1, piflufolastat) ----
    kon <- fixed(0.046); label("PSMA association rate constant (L/nmol/min)") # Table 1, k on = 0.046 L/nmol/min
    kd <- fixed(0.49); label("PSMA dissociation constant KD (nmol/L)") # Table 1, K d = 0.49 nmol/L
    phi_kid <- fixed(0.75); label("Sieving coefficient of the tracer at the glomerulus (unitless)") # Table 1, Sieving of tracer = 0.75 (Table S1 prints Phi kid = 1; see Errata)
    fraex <- fixed(0.96); label("Excreted fraction of the filtered load (unitless)") # Table 1, Fraex = 0.96 (Table S1 prints Fra ex = 1; see Errata)
    aspec <- fixed(0.045); label("Specific activity (GBq/nmol)") # Table 1, A spec = 0.045 GBq/nmol

    # ---- Internalisation and release (Table S1 section 1a) --------------------
    kint <- fixed(0.001); label("Internalisation rate constant, all PSMA tissues (1/min)") # Table S1 1a, Labda int,tu = 0.001; every other tissue = Labda int,tu
    krel_tumor <- fixed(1.4e-4); label("Release rate of internalised ligand, tumours (1/min)") # Table S1 1a, Labda rel,tu = 1.4e-4; rest tumour equal
    krel_kidney <- fixed(2.3e-4); label("Release rate of internalised ligand, kidney and normal tissue (1/min)") # Table S1 1a, Labda rel,kid = 2.3e-4; nt, prostate, liver, spleen, GI equal
    krel_parotid <- fixed(3.7e-4); label("Release rate of internalised ligand, salivary and lacrimal glands (1/min)") # Table S1 1a, Labda rel,par = 0.00037; sm, lg equal

    # ---- Tumour perfusion and receptor density (Table S1 section 1b) ----------
    fdens_tumor1 <- fixed(0.1351538); label("Serum flow density, tumour lesion 1 (mL/min/g)") # Table S1, F dens,tu1 = 0.1351538; Results 'fixed at 0.135'
    fdens_tumor2 <- fixed(0.2676154); label("Serum flow density, tumour lesion 2 (mL/min/g)") # Table S1, F dens,tu2 = 0.2676154
    xf_tumorrest <- fixed(0.53); label("Ratio of actual to assumed flow density, rest tumour (unitless)") # Table S1, X F = 0.53
    rdens_tumor1 <- fixed(44.93846); label("PSMA receptor density, tumour lesion 1 (nmol/L)") # Table S1, R dens,tu1 = 44.93846
    rdens_tumor2 <- fixed(47.38462); label("PSMA receptor density, tumour lesion 2 (nmol/L)") # Table S1, R dens,tu2 = 47.38462
    xr_tumorrest <- fixed(1.4); label("Ratio of assumed to actual receptor density, rest tumour (unitless)") # Table S1, X R = 1.4
    rdens_kidney <- fixed(18); label("PSMA receptor density, kidneys (nmol/L)") # Table S1, R dens,kid = 18
    rdens_parotid <- fixed(42); label("PSMA receptor density, parotid, submandibular and lacrimal glands (nmol/L)") # Table S1, R dens,par = 42; sm, lg equal
  })

  model({
    # Volumes in L, flows and PS products in L/min, amounts in nmol, time in
    # min. Densities in mL/min/g times a volume in L give L/min directly
    # (tissue density 1 g/mL).
    hfrac <- HCT / 100
    bsa <- 0.007184 * HT^0.725 * WT^0.425 # Du Bois, Table S1

    # ---- Total tissue volumes (Table S1) -------------------------------------
    vtot_ven <- 0.0452 * WT
    vtot_art <- 0.0224 * WT
    vtot_lungs <- WT / 71
    vtot_tumor1 <- ORGVOL_TUMOR1 / 1000
    vtot_tumor2 <- ORGVOL_TUMOR2 / 1000
    vtot_tumorrest <- ORGVOL_TUMORREST / 1000
    vtot_prostate <- 0.016 * WT / 71 # Table S1 prints (0.016+BW)/71; see Errata
    vtot_kidney <- 0.3
    vtot_liver <- 0.023 * WT
    vtot_spleen <- 0.003 * WT
    vtot_gi <- (0.395 + 0.548 + 0.104 + 0.15) * WT / 71
    vtot_parotid <- 0.032
    vtot_submandibular <- 0.0095
    vtot_lacrimal <- 0.00068
    vtot_redmarrow <- 1.1 * WT / 71
    vtot_muscle <- 30.078 * WT / 71
    vtot_skin <- 3.408 * WT / 71
    vtot_brain <- 1.45 * WT / 71
    vtot_adipose <- 13.465 * WT / 71
    vtot_heart <- 0.341 * WT / 71
    vtot_bone <- 10.165 * WT / 71 - vtot_redmarrow
    vtot_rest <- WT - (vtot_ven + vtot_art + vtot_lungs + vtot_tumor1 + vtot_tumor2 +
      vtot_tumorrest + vtot_prostate + vtot_kidney + vtot_liver + vtot_spleen +
      vtot_gi + vtot_parotid + vtot_submandibular + vtot_lacrimal + vtot_redmarrow +
      vtot_muscle + vtot_skin + vtot_brain + vtot_adipose + vtot_heart + vtot_bone)

    # ---- Serum and vascular volumes (Table S1) -------------------------------
    vser <- 2.8 * (1 - hfrac) * bsa
    vser_ven <- 0.18 * vser + 0.045 * vser
    vser_art <- 0.06 * vser + 0.045 * vser
    fvasc_tumor <- 0.05 * (1 - hfrac)
    vv_lungs <- vtot_lungs * 0.055
    vv_tumor1 <- vtot_tumor1 * fvasc_tumor
    vv_tumor2 <- vtot_tumor2 * fvasc_tumor
    vv_tumorrest <- vtot_tumorrest * fvasc_tumor
    vv_prostate <- 0.004 * (1 - hfrac) * vtot_prostate
    vv_kidney <- vtot_kidney * 0.055
    vv_liver <- vtot_liver * 0.085
    vv_spleen <- vtot_spleen * 0.12
    vv_gi <- 0.076 * vser
    vv_parotid <- 0.03 * (1 - hfrac) * vtot_parotid
    vv_submandibular <- 0.03 * (1 - hfrac) * vtot_submandibular
    vv_lacrimal <- 0.03 * (1 - hfrac) * vtot_lacrimal
    vv_redmarrow <- 0.04 * vser
    vv_muscle <- 0.14 * vser
    vv_skin <- 0.03 * vser
    vv_brain <- 0.012 * vser
    vv_adipose <- 0.05 * vser
    vv_heart <- 0.01 * vser
    vv_bone <- 0.07 * vser - vv_redmarrow
    vv_rest <- vser - (vser_ven + vser_art + vv_lungs + vv_tumor1 + vv_tumor2 +
      vv_tumorrest + vv_prostate + vv_kidney + vv_liver + vv_spleen + vv_gi +
      vv_parotid + vv_submandibular + vv_lacrimal + vv_redmarrow + vv_muscle +
      vv_skin + vv_brain + vv_adipose + vv_heart + vv_bone)

    # ---- Interstitial volumes (Table S1); the kidney lumen is Fra_int,kid ----
    vi_lungs <- vtot_lungs * 0.3
    vi_tumor1 <- vtot_tumor1 * 0.38
    vi_tumor2 <- vtot_tumor2 * 0.38
    vi_tumorrest <- vtot_tumorrest * 0.38
    vi_prostate <- 0.25 * vtot_prostate
    vlumen_kidney <- vtot_kidney * 0.15
    vi_liver <- vtot_liver * 0.2
    vi_spleen <- vtot_spleen * 0.2
    vi_gi <- vtot_gi * 0.1739634
    vi_parotid <- 0.23 * vtot_parotid
    vi_submandibular <- 0.23 * vtot_submandibular
    vi_lacrimal <- 0.23 * vtot_lacrimal
    vi_redmarrow <- vv_redmarrow * 3.7
    vi_muscle <- vv_muscle * 5.9
    vi_skin <- vv_skin * 8.9
    vi_adipose <- vv_adipose * 15.5
    vi_heart <- vv_heart * 3.7
    vi_bone <- vv_bone * 8.4
    vi_rest <- vtot_rest * 0.1712696

    # ---- Serum flows (Table S1) ----------------------------------------------
    ftot <- vser * 1.23
    fdens_tumorrest <- xf_tumorrest * (fdens_tumor1 + fdens_tumor2) / 2
    fdens_prostate <- 0.18 * (1 - hfrac)
    fq_tumor1 <- fdens_tumor1 * vtot_tumor1
    fq_tumor2 <- fdens_tumor2 * vtot_tumor2
    fq_tumorrest <- fdens_tumorrest * vtot_tumorrest # Table S1 prints F_dens,pro here; see Errata
    fq_prostate <- fdens_prostate * vtot_prostate
    fq_kidney <- (4.3 - 0.026 * AGE) * vtot_kidney * (1 - hfrac)
    fq_liver <- 0.065 * ftot
    fq_spleen <- 0.03 * ftot
    fq_gi <- 0.16 * ftot
    fq_parotid <- 0.16 * vtot_parotid
    fq_submandibular <- 0.16 * vtot_submandibular
    fq_lacrimal <- 0.16 * vtot_lacrimal
    fq_redmarrow <- 0.03 * ftot
    fq_muscle <- 0.17 * ftot
    fq_skin <- 0.05 * ftot
    fq_brain <- 0.12 * ftot
    fq_adipose <- 0.05 * ftot
    fq_heart <- 0.04 * ftot
    fq_bone <- 0.05 * ftot
    fq_rest <- ftot - (fq_tumor1 + fq_tumor2 + fq_tumorrest + fq_prostate +
      fq_kidney + fq_liver + fq_spleen + fq_gi + fq_parotid + fq_submandibular +
      fq_lacrimal + fq_redmarrow + fq_muscle + fq_skin + fq_brain + fq_adipose +
      fq_heart + fq_bone)
    fq_portal <- fq_liver + fq_spleen + fq_gi # hepatic venous outflow

    # ---- Renal filtration (Table S1, Table 1) --------------------------------
    f_fil <- CRCL / 1000 * phi_kid
    f_ex <- f_fil * fraex
    f_ret <- f_fil - f_ex

    # ---- Permeability-surface-area products (Table S1) -----------------------
    ps_lungs <- 2 * vtot_lungs # K_lu = K_mu * 100
    ps_tumor1 <- 0.6 * vtot_tumor1
    ps_tumor2 <- 0.6 * vtot_tumor2
    ps_tumorrest <- 0.6 * vtot_tumorrest
    ps_prostate <- 0.1 * vtot_prostate # K_pro listed twice (0.1 and 0.02); see Errata
    ps_liver <- 2 * vtot_liver
    ps_spleen <- 2 * vtot_spleen
    ps_gi <- 0.02 * vtot_gi
    ps_parotid <- 0.4 * vtot_parotid
    ps_submandibular <- 0.4 * vtot_submandibular
    ps_lacrimal <- 0.4 * vtot_lacrimal
    ps_redmarrow <- 2 * vtot_redmarrow
    ps_muscle <- 0.02 * vtot_muscle
    ps_skin <- 0.02 * vtot_skin
    ps_adipose <- 0.02 * vtot_adipose
    ps_heart <- 0.02 * vtot_heart
    ps_bone <- 0.02 * vtot_bone
    ps_rest <- 0.02 * vtot_rest # K_rest not tabulated; muscle-like value, see Errata

    # ---- Receptor densities and initial free receptors (Table S1) ------------
    koff <- kd * kon
    rdens_tumorrest <- xr_tumorrest * (rdens_tumor1 + rdens_tumor2) / 2
    rdens_prostate <- rdens_tumorrest * 0.1
    rdens_liver <- rdens_prostate * 0.05
    rdens_spleen <- rdens_kidney * 0.2
    rdens_gi <- rdens_prostate * 0.06
    rfree_tumor1(0) <- rdens_tumor1 * vtot_tumor1
    rfree_tumor2(0) <- rdens_tumor2 * vtot_tumor2
    rfree_tumorrest(0) <- rdens_tumorrest * vtot_tumorrest
    rfree_prostate(0) <- rdens_prostate * vtot_prostate
    rfree_kidney(0) <- rdens_kidney * vtot_kidney
    rfree_liver(0) <- rdens_liver * vtot_liver
    rfree_spleen(0) <- rdens_spleen * vtot_spleen
    rfree_gi(0) <- rdens_gi * vtot_gi
    rfree_parotid(0) <- rdens_parotid * vtot_parotid
    rfree_submandibular(0) <- rdens_parotid * vtot_submandibular
    rfree_lacrimal(0) <- rdens_parotid * vtot_lacrimal

    cart <- arterial / vser_art

    # ---- Veins and arteries (Suppl. Eq 1-2) ----------------------------------
    d/dt(venous) <- fq_tumor1 * vasc_tumor1 / vv_tumor1 + fq_tumor2 * vasc_tumor2 / vv_tumor2 +
      fq_tumorrest * vasc_tumorrest / vv_tumorrest + fq_prostate * vasc_prostate / vv_prostate +
      fq_kidney * vasc_kidney / vv_kidney + fq_portal * vasc_liver / vv_liver +
      fq_parotid * vasc_parotid / vv_parotid + fq_submandibular * vasc_submandibular / vv_submandibular +
      fq_lacrimal * vasc_lacrimal / vv_lacrimal + fq_redmarrow * vasc_redmarrow / vv_redmarrow +
      fq_muscle * vasc_muscle / vv_muscle + fq_skin * vasc_skin / vv_skin +
      fq_brain * vasc_brain / vv_brain + fq_adipose * vasc_adipose / vv_adipose +
      fq_heart * vasc_heart / vv_heart + fq_bone * vasc_bone / vv_bone +
      fq_rest * vasc_rest / vv_rest - ftot * venous / vser_ven
    d/dt(arterial) <- ftot * vasc_lungs / vv_lungs - ftot * cart

    # ---- Lungs (Suppl. Eq 4, 9) ----------------------------------------------
    d/dt(vasc_lungs) <- ftot * (venous / vser_ven - vasc_lungs / vv_lungs) +
      ps_lungs * (int_lungs / vi_lungs - vasc_lungs / vv_lungs)
    d/dt(int_lungs) <- ps_lungs * (vasc_lungs / vv_lungs - int_lungs / vi_lungs)

    # ---- PSMA-expressing tissues (Suppl. Eq 5, 7, 8, 11, 13, 14) -------------
    d/dt(vasc_tumor1) <- fq_tumor1 * (cart - vasc_tumor1 / vv_tumor1) +
      ps_tumor1 * (int_tumor1 / vi_tumor1 - vasc_tumor1 / vv_tumor1)
    d/dt(int_tumor1) <- ps_tumor1 * (vasc_tumor1 / vv_tumor1 - int_tumor1 / vi_tumor1) -
      kon * int_tumor1 * rfree_tumor1 / vi_tumor1 + koff * bound_tumor1
    d/dt(bound_tumor1) <- kon * int_tumor1 * rfree_tumor1 / vi_tumor1 - (koff + kint) * bound_tumor1
    d/dt(intern_tumor1) <- kint * bound_tumor1 - krel_tumor * intern_tumor1
    d/dt(rfree_tumor1) <- (koff + kint) * bound_tumor1 - kon * int_tumor1 * rfree_tumor1 / vi_tumor1

    d/dt(vasc_tumor2) <- fq_tumor2 * (cart - vasc_tumor2 / vv_tumor2) +
      ps_tumor2 * (int_tumor2 / vi_tumor2 - vasc_tumor2 / vv_tumor2)
    d/dt(int_tumor2) <- ps_tumor2 * (vasc_tumor2 / vv_tumor2 - int_tumor2 / vi_tumor2) -
      kon * int_tumor2 * rfree_tumor2 / vi_tumor2 + koff * bound_tumor2
    d/dt(bound_tumor2) <- kon * int_tumor2 * rfree_tumor2 / vi_tumor2 - (koff + kint) * bound_tumor2
    d/dt(intern_tumor2) <- kint * bound_tumor2 - krel_tumor * intern_tumor2
    d/dt(rfree_tumor2) <- (koff + kint) * bound_tumor2 - kon * int_tumor2 * rfree_tumor2 / vi_tumor2

    d/dt(vasc_tumorrest) <- fq_tumorrest * (cart - vasc_tumorrest / vv_tumorrest) +
      ps_tumorrest * (int_tumorrest / vi_tumorrest - vasc_tumorrest / vv_tumorrest)
    d/dt(int_tumorrest) <- ps_tumorrest * (vasc_tumorrest / vv_tumorrest - int_tumorrest / vi_tumorrest) -
      kon * int_tumorrest * rfree_tumorrest / vi_tumorrest + koff * bound_tumorrest
    d/dt(bound_tumorrest) <- kon * int_tumorrest * rfree_tumorrest / vi_tumorrest - (koff + kint) * bound_tumorrest
    d/dt(intern_tumorrest) <- kint * bound_tumorrest - krel_tumor * intern_tumorrest
    d/dt(rfree_tumorrest) <- (koff + kint) * bound_tumorrest - kon * int_tumorrest * rfree_tumorrest / vi_tumorrest

    d/dt(vasc_prostate) <- fq_prostate * (cart - vasc_prostate / vv_prostate) +
      ps_prostate * (int_prostate / vi_prostate - vasc_prostate / vv_prostate)
    d/dt(int_prostate) <- ps_prostate * (vasc_prostate / vv_prostate - int_prostate / vi_prostate) -
      kon * int_prostate * rfree_prostate / vi_prostate + koff * bound_prostate
    d/dt(bound_prostate) <- kon * int_prostate * rfree_prostate / vi_prostate - (koff + kint) * bound_prostate
    d/dt(intern_prostate) <- kint * bound_prostate - krel_kidney * intern_prostate
    d/dt(rfree_prostate) <- (koff + kint) * bound_prostate - kon * int_prostate * rfree_prostate / vi_prostate

    # Liver receives arterial, splenic and GI inflow (Suppl. Eq 5)
    d/dt(vasc_liver) <- fq_liver * cart + fq_spleen * vasc_spleen / vv_spleen +
      fq_gi * vasc_gi / vv_gi - fq_portal * vasc_liver / vv_liver +
      ps_liver * (int_liver / vi_liver - vasc_liver / vv_liver)
    d/dt(int_liver) <- ps_liver * (vasc_liver / vv_liver - int_liver / vi_liver) -
      kon * int_liver * rfree_liver / vi_liver + koff * bound_liver
    d/dt(bound_liver) <- kon * int_liver * rfree_liver / vi_liver - (koff + kint) * bound_liver
    d/dt(intern_liver) <- kint * bound_liver - krel_kidney * intern_liver
    d/dt(rfree_liver) <- (koff + kint) * bound_liver - kon * int_liver * rfree_liver / vi_liver

    d/dt(vasc_spleen) <- fq_spleen * (cart - vasc_spleen / vv_spleen) +
      ps_spleen * (int_spleen / vi_spleen - vasc_spleen / vv_spleen)
    d/dt(int_spleen) <- ps_spleen * (vasc_spleen / vv_spleen - int_spleen / vi_spleen) -
      kon * int_spleen * rfree_spleen / vi_spleen + koff * bound_spleen
    d/dt(bound_spleen) <- kon * int_spleen * rfree_spleen / vi_spleen - (koff + kint) * bound_spleen
    d/dt(intern_spleen) <- kint * bound_spleen - krel_kidney * intern_spleen
    d/dt(rfree_spleen) <- (koff + kint) * bound_spleen - kon * int_spleen * rfree_spleen / vi_spleen

    d/dt(vasc_gi) <- fq_gi * (cart - vasc_gi / vv_gi) +
      ps_gi * (int_gi / vi_gi - vasc_gi / vv_gi)
    d/dt(int_gi) <- ps_gi * (vasc_gi / vv_gi - int_gi / vi_gi) -
      kon * int_gi * rfree_gi / vi_gi + koff * bound_gi
    d/dt(bound_gi) <- kon * int_gi * rfree_gi / vi_gi - (koff + kint) * bound_gi
    d/dt(intern_gi) <- kint * bound_gi - krel_kidney * intern_gi
    d/dt(rfree_gi) <- (koff + kint) * bound_gi - kon * int_gi * rfree_gi / vi_gi

    d/dt(vasc_parotid) <- fq_parotid * (cart - vasc_parotid / vv_parotid) +
      ps_parotid * (int_parotid / vi_parotid - vasc_parotid / vv_parotid)
    d/dt(int_parotid) <- ps_parotid * (vasc_parotid / vv_parotid - int_parotid / vi_parotid) -
      kon * int_parotid * rfree_parotid / vi_parotid + koff * bound_parotid
    d/dt(bound_parotid) <- kon * int_parotid * rfree_parotid / vi_parotid - (koff + kint) * bound_parotid
    d/dt(intern_parotid) <- kint * bound_parotid - krel_parotid * intern_parotid
    d/dt(rfree_parotid) <- (koff + kint) * bound_parotid - kon * int_parotid * rfree_parotid / vi_parotid

    d/dt(vasc_submandibular) <- fq_submandibular * (cart - vasc_submandibular / vv_submandibular) +
      ps_submandibular * (int_submandibular / vi_submandibular - vasc_submandibular / vv_submandibular)
    d/dt(int_submandibular) <- ps_submandibular * (vasc_submandibular / vv_submandibular - int_submandibular / vi_submandibular) -
      kon * int_submandibular * rfree_submandibular / vi_submandibular + koff * bound_submandibular
    d/dt(bound_submandibular) <- kon * int_submandibular * rfree_submandibular / vi_submandibular - (koff + kint) * bound_submandibular
    d/dt(intern_submandibular) <- kint * bound_submandibular - krel_parotid * intern_submandibular
    d/dt(rfree_submandibular) <- (koff + kint) * bound_submandibular - kon * int_submandibular * rfree_submandibular / vi_submandibular

    d/dt(vasc_lacrimal) <- fq_lacrimal * (cart - vasc_lacrimal / vv_lacrimal) +
      ps_lacrimal * (int_lacrimal / vi_lacrimal - vasc_lacrimal / vv_lacrimal)
    d/dt(int_lacrimal) <- ps_lacrimal * (vasc_lacrimal / vv_lacrimal - int_lacrimal / vi_lacrimal) -
      kon * int_lacrimal * rfree_lacrimal / vi_lacrimal + koff * bound_lacrimal
    d/dt(bound_lacrimal) <- kon * int_lacrimal * rfree_lacrimal / vi_lacrimal - (koff + kint) * bound_lacrimal
    d/dt(intern_lacrimal) <- kint * bound_lacrimal - krel_parotid * intern_lacrimal
    d/dt(rfree_lacrimal) <- (koff + kint) * bound_lacrimal - kon * int_lacrimal * rfree_lacrimal / vi_lacrimal

    # ---- Kidney (Suppl. Eq 3, 10, 12, 13, 15) --------------------------------
    # The lumen loses F_ex to urine and returns F_ret to the kidney vasculature.
    # Suppl. Eq 10 prints only the F_ex term; F_ret is subtracted here so that
    # the F_ret return in Suppl. Eq 3 conserves mass (see Errata).
    d/dt(vasc_kidney) <- fq_kidney * (cart - vasc_kidney / vv_kidney) -
      f_fil * vasc_kidney / vv_kidney + f_ret * lumen_kidney / vlumen_kidney
    d/dt(lumen_kidney) <- f_fil * vasc_kidney / vv_kidney - (f_ex + f_ret) * lumen_kidney / vlumen_kidney -
      kon * lumen_kidney * rfree_kidney / vlumen_kidney + koff * bound_kidney
    d/dt(bound_kidney) <- kon * lumen_kidney * rfree_kidney / vlumen_kidney - (koff + kint) * bound_kidney
    d/dt(intern_kidney) <- kint * bound_kidney - krel_kidney * intern_kidney
    d/dt(rfree_kidney) <- (koff + kint) * bound_kidney - kon * lumen_kidney * rfree_kidney / vlumen_kidney

    # ---- PSMA-negative tissues (Suppl. Eq 6, 7, 9) ---------------------------
    d/dt(vasc_redmarrow) <- fq_redmarrow * (cart - vasc_redmarrow / vv_redmarrow) +
      ps_redmarrow * (int_redmarrow / vi_redmarrow - vasc_redmarrow / vv_redmarrow)
    d/dt(int_redmarrow) <- ps_redmarrow * (vasc_redmarrow / vv_redmarrow - int_redmarrow / vi_redmarrow)
    d/dt(vasc_muscle) <- fq_muscle * (cart - vasc_muscle / vv_muscle) +
      ps_muscle * (int_muscle / vi_muscle - vasc_muscle / vv_muscle)
    d/dt(int_muscle) <- ps_muscle * (vasc_muscle / vv_muscle - int_muscle / vi_muscle)
    d/dt(vasc_skin) <- fq_skin * (cart - vasc_skin / vv_skin) +
      ps_skin * (int_skin / vi_skin - vasc_skin / vv_skin)
    d/dt(int_skin) <- ps_skin * (vasc_skin / vv_skin - int_skin / vi_skin)
    d/dt(vasc_brain) <- fq_brain * (cart - vasc_brain / vv_brain)
    d/dt(vasc_adipose) <- fq_adipose * (cart - vasc_adipose / vv_adipose) +
      ps_adipose * (int_adipose / vi_adipose - vasc_adipose / vv_adipose)
    d/dt(int_adipose) <- ps_adipose * (vasc_adipose / vv_adipose - int_adipose / vi_adipose)
    d/dt(vasc_heart) <- fq_heart * (cart - vasc_heart / vv_heart) +
      ps_heart * (int_heart / vi_heart - vasc_heart / vv_heart)
    d/dt(int_heart) <- ps_heart * (vasc_heart / vv_heart - int_heart / vi_heart)
    d/dt(vasc_bone) <- fq_bone * (cart - vasc_bone / vv_bone) +
      ps_bone * (int_bone / vi_bone - vasc_bone / vv_bone)
    d/dt(int_bone) <- ps_bone * (vasc_bone / vv_bone - int_bone / vi_bone)
    d/dt(vasc_rest) <- fq_rest * (cart - vasc_rest / vv_rest) +
      ps_rest * (int_rest / vi_rest - vasc_rest / vv_rest)
    d/dt(int_rest) <- ps_rest * (vasc_rest / vv_rest - int_rest / vi_rest)

    # ---- Cumulative elimination (bookkeeping; closes the mass balance) -------
    d/dt(urine) <- f_ex * lumen_kidney / vlumen_kidney
    d/dt(released) <- krel_tumor * (intern_tumor1 + intern_tumor2 + intern_tumorrest) +
      krel_kidney * (intern_prostate + intern_liver + intern_spleen + intern_gi + intern_kidney) +
      krel_parotid * (intern_parotid + intern_submandibular + intern_lacrimal)

    # ---- Outputs -------------------------------------------------------------
    # Arterial serum concentration (nmol/L) and decay-corrected tissue activity
    # concentrations (kBq/mL): nmol/L * aspec (GBq/nmol) * 1e6 kBq/GBq / 1000 mL/L.
    Cc <- cart
    kbq <- aspec * 1000
    act_blood <- cart * kbq
    act_tumor1 <- (vasc_tumor1 + int_tumor1 + bound_tumor1 + intern_tumor1) / vtot_tumor1 * kbq
    act_tumor2 <- (vasc_tumor2 + int_tumor2 + bound_tumor2 + intern_tumor2) / vtot_tumor2 * kbq
    act_lungs <- (vasc_lungs + int_lungs) / vtot_lungs * kbq
    act_muscle <- (vasc_muscle + int_muscle) / vtot_muscle * kbq
    act_kidney <- (vasc_kidney + lumen_kidney + bound_kidney + intern_kidney) / vtot_kidney * kbq
    act_liver <- (vasc_liver + int_liver + bound_liver + intern_liver) / vtot_liver * kbq
    act_spleen <- (vasc_spleen + int_spleen + bound_spleen + intern_spleen) / vtot_spleen * kbq
    act_parotid <- (vasc_parotid + int_parotid + bound_parotid + intern_parotid) / vtot_parotid * kbq
  })
}
