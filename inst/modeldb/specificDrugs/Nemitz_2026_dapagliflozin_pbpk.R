Nemitz_2026_dapagliflozin_pbpk <- function() {
  description <- paste(
    "PBPK/PD (whole-body, 31 ODEs, SBML/libroadrunner).",
    "Dapagliflozin and its major UGT1A9 metabolite dapagliflozin-3-O-",
    "glucuronide (D3OG) in type 2 diabetes mellitus, built by Nemitz,",
    "Elias and Koenig (2026) from 28 curated clinical studies. The",
    "whole-body circulation links gut, liver, kidney and lung submodels",
    "through venous, arterial, portal and hepatic-venous plasma pools,",
    "with a lumped rest-of-body compartment. Oral drug dissolves from the",
    "dose compartment into the gut lumen, 84% is absorbed to the portal",
    "vein and 16% is excreted in feces. Dapagliflozin is taken into liver",
    "and kidney tissue by carrier-mediated transport and glucuronidated to",
    "D3OG by UGT1A9 (irreversible Michaelis-Menten, renal activity ten-fold",
    "the hepatic rate per litre of tissue); parent and metabolite are then",
    "excreted in urine by first-order renal processes. The pharmacodynamic",
    "layer links kidney plasma dapagliflozin to urinary glucose excretion",
    "through an Emax reduction of the renal threshold for glucose (RTG),",
    "with glucose excreted only while fasting plasma glucose exceeds RTG.",
    "Three scenario scalars reproduce the paper's variability analyses:",
    "relative renal function (KDIGO categories) scales GFR together with",
    "renal transport, metabolism and excretion; relative liver function",
    "sets the cirrhotic parenchymal-volume loss and portal/arterial",
    "shunting (Child-Turcotte-Pugh classes); and a fed-state indicator cuts",
    "the intestinal absorption rate to 30% of fasted. The model is",
    "deterministic - the authors fit a typical individual and report no",
    "between-subject variability or residual error. Two other dapagliflozin",
    "models are available: modellib('vanderWalt_2013_dapagliflozin') for",
    "the parent/D3OG population PK model and",
    "modellib('Kobuchi_2025_dapagliflozin')."
  )
  reference <- paste(
    "Nemitz N, Elias M, Koenig M. A Physiologically Based Pharmacokinetic",
    "and Pharmacodynamic (PBPK/PD) Model of Dapagliflozin in Type 2",
    "Diabetes Mellitus: The Effect of Dosing, Hepatorenal Impairment, and",
    "Food. Pharmaceutics. 2026;18(3):287.",
    "doi:10.3390/pharmaceutics18030287.",
    "Submodel rate laws and ODEs from Supplementary Materials Section S4",
    "(Intestine, Liver, Kidney and Pharmacodynamics Models); optimized",
    "pharmacokinetic parameters from Supplementary Table S2 and",
    "pharmacodynamic parameters from Supplementary Table S3; fractional",
    "organ volumes, blood flows, partition coefficient and the renal,",
    "hepatic and prandial scaling factors from Methods Section 2.2; fixed",
    "Michaelis constants from Methods Section 2.3 and Supplementary S4.",
    "The whole-body circulation ODEs are not written out in the",
    "Supplementary Materials and were taken from the paper's own archived",
    "SBML model, cited in Methods Section 2.2 and Supplementary S4 as the",
    "definitive model source: Koenig M. matthiaskoenig/dapagliflozin-model",
    "v0.9.8, Zenodo. doi:10.5281/zenodo.18011516 (file",
    "models/dapagliflozin_body_flat.md). See the vignette Errata for the",
    "pharmacodynamic parameter divergence between Table S3 and that",
    "archive.",
    sep = " "
  )
  vignette <- "Nemitz_2026_dapagliflozin_pbpk"
  units <- list(time = "min", dosing = "mg", concentration = "umol/L")

  # Both routes the paper simulates: oral (all but one curated study) and the
  # single 80 microgram intravenous dose of Boulton 2013. `depot_iv` is outside
  # the auto-detected depot/central set, so it is declared explicitly.
  dosing <- c("depot", "depot_iv")

  # States are CONCENTRATIONS (mmol/L) throughout the circulation, following
  # the source SBML, so doses must be given into `depot` (oral, mg) or
  # `depot_iv` (intravenous, mg) - never directly into a plasma or tissue
  # compartment.
  paper_specific_compartments <- c(
    "depot_iv", "gut_lumen", "gut_plasma", "gut_plasma_d3og",
    "portal", "portal_d3og", "liver_plasma", "liver_plasma_d3og",
    "liver_d3og", "hepatic_vein", "hepatic_vein_d3og",
    "kidney_plasma", "kidney_plasma_d3og", "kidney_d3og",
    "lung_plasma", "lung_plasma_d3og", "rest", "rest_plasma",
    "rest_plasma_d3og", "venous", "venous_d3og", "arterial",
    "arterial_d3og", "feces", "urine_d3og", "glc_urine"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scales every absolute organ volume and, through cardiac output",
        "(COBW = 1.548 mL/s/kg), every blood flow: Methods Section 2.2",
        "states that 'absolute organ volumes and blood flows were",
        "calculated by scaling the corresponding fractional values with",
        "body weight'. The reference individual is 75 kg (source SBML",
        "parameter BW). Methods Section 2.5 states that study-specific mean",
        "bodyweight was applied where reported."
      ),
      source_name        = "BW"
    ),
    FPG = list(
      description        = "Fasting plasma glucose concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the entire pharmacodynamic layer. Model Assumptions",
        "(Section 2.3) state that diurnal variation in plasma glucose was",
        "not modelled and a constant fasting plasma glucose was assumed:",
        "study-specific values were used when reported, otherwise 5 mM for",
        "healthy subjects and 7.5 mM for both type 1 and type 2 diabetes",
        "mellitus. FPG enters twice - it raises the baseline renal",
        "threshold for glucose through RTG_m_fpg above the healthy",
        "reference of 5 mM, and it is the concentration filtered at the",
        "glomerulus, so glucose is excreted only while FPG exceeds RTG.",
        "In the source SBML this is the state KI__glc_ext, which has a",
        "zero derivative and is therefore exactly a covariate; carrying it",
        "as a data column (rather than a state) additionally allows a",
        "time-varying glucose profile to be supplied."
      ),
      source_name        = "KI__glc_ext"
    ),
    RENALFUNC_REL = list(
      description        = "Relative renal function as a fraction of normal",
      units              = "(dimensionless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "1 = normal renal function. Methods Section 2.2 maps KDIGO",
        "categories onto this scalar: normal (eGFR >= 90 mL/min) = 1.00,",
        "mild (GFR 50-89) = 0.69, moderate (GFR 30-49) = 0.32 and severe",
        "(GFR <= 30) = 0.19. It multiplies the glomerular filtration rate",
        "AND the renal transport, metabolism and excretion terms, so it",
        "reduces urinary glucose excretion through two distinct routes:",
        "less glucose filtered, and less dapagliflozin delivered into",
        "kidney tissue to inhibit reabsorption. This dual action is why the",
        "paper reports that renal impairment cut urinary glucose excretion",
        "by 40-60% while barely changing parent plasma exposure.",
        "Source name f_renal_function (SBML KI__f_renal_function)."
      ),
      source_name        = "f_renal_function"
    ),
    HEPFUNC_REL = list(
      description        = "Relative liver function as a fraction of normal",
      units              = "(dimensionless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "1 = normal liver function. The source parameterises hepatic",
        "impairment with the OPPOSITE orientation, as a cirrhosis severity",
        "f_cirrhosis on [0, 0.95], so this column is the complement:",
        "f_cirrhosis = 1 - HEPFUNC_REL. Methods Section 2.2 maps",
        "Child-Turcotte-Pugh classes onto f_cirrhosis as A (mild, 5-6",
        "points) = 0.40, B (moderate, 7-9) = 0.70 and C (severe, 10-15) =",
        "0.80, which are HEPFUNC_REL values of 0.60, 0.30 and 0.20; normal",
        "liver function is f_cirrhosis = 0, i.e. HEPFUNC_REL = 1.",
        "f_cirrhosis acts twice in the source model, as the fraction of",
        "parenchymal liver volume lost (shrinking the metabolising tissue",
        "volume Vli_tissue) and as the fraction of portal and hepatic-",
        "arterial blood shunted past the liver. The canonical column is",
        "used with the fraction-of-normal orientation because that is the",
        "register's definition; take care with the inversion when",
        "reproducing the paper's figures, which are labelled by",
        "f_cirrhosis."
      ),
      source_name        = "f_cirrhosis"
    ),
    FED = list(
      description        = "Fed state at the time of dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste(
        "1 = dose taken in the fed state. Methods Section 2.2 implements",
        "the food effect as the intestinal absorption scaling factor",
        "f_absorption, fasted = 1.00 and fed = 0.30, so the model computes",
        "f_absorption = 1 - 0.7 * FED. Section 2.3 gives the reasoning: a",
        "standard high-fat meal lengthens gastric emptying from under 30",
        "min to about 100 min or more, and the absorption rate was assumed",
        "to change in inverse proportion, giving a fed-state factor of 0.3.",
        "Because the factor scales absorption RATE and not the absorbed",
        "FRACTION, it lowers Cmax while leaving AUC essentially unchanged,",
        "which is the food effect the paper reproduces (Cmax down 30-50%).",
        "Note that the paper's Figure 6 scans f_absorption continuously",
        "over 0.1-10.0; to reproduce that scan, set f_absorption directly",
        "rather than through this binary indicator."
      ),
      source_name        = "f_absorption"
    )
  )

  # Plasma and tissue states are CONCENTRATIONS (mmol/L) as in the source
  # SBML; only the two dosing compartments and the cumulative excretion
  # pools are amounts.
  compartmentData <- list(
    depot                = list(analyte = "dapagliflozin", units = "mg",     specimen = "administration site", verified = TRUE),
    depot_iv             = list(analyte = "dapagliflozin", units = "mg",     specimen = "administration site", verified = TRUE),
    gut_lumen            = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "administration site", verified = TRUE),
    feces                = list(analyte = "dapagliflozin", units = "mmol",   specimen = "faeces",              verified = TRUE),
    gut_plasma           = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    gut_plasma_d3og      = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    portal               = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    portal_d3og          = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    liver_plasma         = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    liver_plasma_d3og    = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    liver                = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "tissue",              verified = TRUE),
    liver_d3og           = list(analyte = "D3OG",          units = "mmol/L", specimen = "tissue",              verified = TRUE),
    hepatic_vein         = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    hepatic_vein_d3og    = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    kidney_plasma        = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    kidney_plasma_d3og   = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    kidney               = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "tissue",              verified = TRUE),
    kidney_d3og          = list(analyte = "D3OG",          units = "mmol/L", specimen = "tissue",              verified = TRUE),
    urine                = list(analyte = "dapagliflozin", units = "mmol",   specimen = "urine",               verified = TRUE),
    urine_d3og           = list(analyte = "D3OG",          units = "mmol",   specimen = "urine",               verified = TRUE),
    lung_plasma          = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    lung_plasma_d3og     = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    lung                 = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "tissue",              verified = TRUE),
    rest_plasma          = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    rest_plasma_d3og     = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    rest                 = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "tissue",              verified = TRUE),
    venous               = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    venous_d3og          = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    arterial             = list(analyte = "dapagliflozin", units = "mmol/L", specimen = "plasma",              verified = TRUE),
    arterial_d3og        = list(analyte = "D3OG",          units = "mmol/L", specimen = "plasma",              verified = TRUE),
    glc_urine            = list(analyte = "glucose",       units = "mmol",   specimen = "urine",               verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = 28L,
    age_range      = "adults; not tabulated per study in the source",
    weight_range   = paste(
      "study-specific mean bodyweight where reported (Methods Section",
      "2.5); 75 kg reference individual otherwise"
    ),
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Healthy volunteers, patients with type 1 or type 2 diabetes",
      "mellitus, and subjects with renal or hepatic impairment. Pediatric",
      "and animal studies were excluded (Methods Section 2.1)."
    ),
    dose_range     = paste(
      "Oral 0.001-500 mg single and multiple dose across the curated",
      "studies (Table 1); 10 mg is the standard maintenance dose. A single",
      "80 microgram intravenous dose (Boulton 2013) informs the",
      "distribution parameters."
    ),
    regions        = "not reported; the curated studies are international",
    notes          = paste(
      "Deterministic typical-individual model. Methods Section 2.2 states",
      "explicitly that 'all simulations were performed deterministically",
      "using the optimized parameter set representing the typical (mean)",
      "individual' and that inter-individual variability was NOT included,",
      "so no IIV and no residual-error model are reported and none are",
      "invented here. Data were curated into PK-DB (study identifiers in",
      "Table 1) from a PubMed and PKPDAI search screening 190 records down",
      "to 28 studies. Parameters were fitted sequentially - pharmacokinetic",
      "first (Table S2), then pharmacodynamic (Table S3) - by weighted",
      "least squares over a subset of single-dose data from healthy",
      "subjects, type 2 diabetics and subjects with renal impairment, with",
      "100 multi-start local optimisation runs. The authors note that",
      "renal threshold for glucose measurements were available from only",
      "one study, fecal excretion from one study, and urinary glucose",
      "excretion under fed conditions from none, so the pharmacodynamic",
      "layer rests on a narrow evidence base."
    )
  )

  ini({
    # ================= Reference physiology ==============================
    # Methods Section 2.2 gives the fractional volumes and blood flows;
    # values below are the source SBML parameters, which carry the extra
    # blood-pool fractions the manuscript text does not tabulate.
    covbw   <- fixed(1.548);  label("Cardiac output per unit body weight (mL/s/kg)")  # SBML COBW
    fcardio <- fixed(1);      label("Cardiac function scaling factor (unitless)")     # SBML f_cardiac_function, normal = 1
    hct     <- fixed(0.51);   label("Hematocrit (fraction)")                          # SBML HCT
    fblood  <- fixed(0.02);   label("Blood fraction of organ volume (unitless)")      # SBML Fblood

    fvar <- fixed(0.0257);  label("Arterial blood fractional volume (L/kg)")        # SBML FVar
    fvgu <- fixed(0.0171);  label("Gut fractional volume (L/kg)")                   # Methods 2.2 'FV gu = 1.71%'
    fvhv <- fixed(0.001);   label("Hepatic venous fractional volume (L/kg)")        # SBML FVhv
    fvki <- fixed(0.0044);  label("Kidney fractional volume (L/kg)")                # Methods 2.2 'FV ki = 0.44%'
    fvli <- fixed(0.021);   label("Liver fractional volume (L/kg)")                 # Methods 2.2 'FV li = 2.10%'
    fvlu <- fixed(0.0076);  label("Lung fractional volume (L/kg)")                  # Methods 2.2 'FV lu = 0.76%'
    fvpo <- fixed(0.001);   label("Portal fractional volume (L/kg)")                # SBML FVpo
    fvve <- fixed(0.0514);  label("Venous blood fractional volume (L/kg)")          # SBML FVve

    fqgu <- fixed(0.18);   label("Gut fractional blood flow (unitless)")                     # Methods 2.2 'FQ gu = 18.00%'
    fqh  <- fixed(0.215);  label("Hepatic venous-side fractional blood flow (unitless)")     # Methods 2.2 'FQ h = 21.50%'
    fqki <- fixed(0.19);   label("Kidney fractional blood flow (unitless)")                  # Methods 2.2 'FQ ki = 19.00%'
    fqlu <- fixed(1);      label("Lung fractional blood flow (unitless)")                    # Methods 2.2 'FQ lu = 100%'

    # ================= Molecular weights and conversions =================
    mr_dap    <- fixed(408.873); label("Molecular weight of dapagliflozin (g/mol)")   # SBML Mr_dap
    mr_glc    <- fixed(180);     label("Molecular weight of glucose (g/mol)")         # SBML KI__Mr_glc
    cf_mg_g   <- fixed(1000);    label("Conversion factor, mg per g")                 # SBML KI__cf_mg_per_g
    cf_ml_l   <- fixed(1000);    label("Conversion factor, mL per L")                 # SBML KI__cf_ml_per_l

    # ================= Distribution ======================================
    # Supplementary Table S2 (optimized). Tissue partitioning was assumed
    # identical across all tissues for dapagliflozin and none was assumed
    # for D3OG (Methods Section 2.2).
    kp_dap      <- fixed(25.517); label("Tissue-to-plasma partition coefficient, dapagliflozin (unitless)")  # Table S2 Kp_dap = 25.517
    ftissue_dap <- fixed(0.01);   label("Tissue distribution rate, dapagliflozin (L/min)")                   # Table S2 ftissue_dap = 0.01

    # ================= Absorption ========================================
    ka_dis_dap <- fixed(0.84842);  label("Dissolution rate of dapagliflozin in the GI tract (1/hr)")  # Table S2 GU__Ka_dis_dap = 0.84842
    dapabs_k   <- fixed(0.05946);  label("Absorption rate of dapagliflozin in the GI tract (1/min)")  # Table S2 GU__DAPABS_k = 0.05946
    f_dap_abs  <- fixed(0.84);     label("Fraction of the dissolved dose absorbed (unitless)")        # Figure 1B / Supplementary S2: 84% absorbed, 16% excreted in feces

    # ================= Hepatic and renal transport =======================
    # Methods Section 2.3: transport of dapagliflozin and D3OG in liver and
    # kidney was assumed fast and reversible relative to metabolic
    # conversion and therefore NOT rate-limiting, so the Vmax values were
    # not fitted; they are held at the source SBML values. The Michaelis
    # constants are literature values from the FDA review cited in
    # Supplementary S4 (OAT3).
    dapim_vmax    <- fixed(10);    label("Vmax, dapagliflozin uptake into liver and kidney (mmol/min/L)")  # SBML LI__DAPIM_Vmax / KI__DAPIM_Vmax; not rate-limiting per Methods 2.3
    dapim_km_dap  <- fixed(0.033); label("Km, dapagliflozin uptake, OAT3 (mmol/L)")                        # Supplementary S4 liver model: 'Km = 33 uM'
    d3gex_vmax    <- fixed(10);    label("Vmax, D3OG export from liver (mmol/min/L)")                      # SBML LI__D3GEX_Vmax; not rate-limiting per Methods 2.3
    d3gex_km_d3og <- fixed(0.115); label("Km, D3OG export from liver, OAT3 (mmol/L)")                      # Supplementary S4 liver model: 'Km = 115 uM'
    d3gim_vmax    <- fixed(10);    label("Vmax, D3OG uptake into kidney (mmol/min/L)")                     # SBML KI__D3GIM_Vmax; not rate-limiting per Methods 2.3
    d3gim_km_d3og <- fixed(0.033); label("Km, D3OG uptake into kidney (mmol/L)")                           # SBML KI__D3GIM_Km_d3g

    # ================= Metabolism (UGT1A9) ===============================
    # Methods Section 2.3: irreversible Michaelis-Menten; the Michaelis
    # constant was FIXED to the kidney-microsome value and only the maximum
    # rates were estimated.
    dap2d3g_vmax   <- fixed(0.01992); label("Vmax, UGT1A9 conversion of dapagliflozin to D3OG in liver (mmol/min/L)")  # Table S2 DAP2D3G_Vmax = 0.01992
    dap2d3g_km_dap <- fixed(0.479);   label("Km, UGT1A9 conversion of dapagliflozin to D3OG (mmol/L)")                 # Methods 2.3: 'Km = 479 uM' from kidney microsomes, fixed not fitted
    f_dap2d3g_ki   <- fixed(10);      label("Renal UGT1A9 activity relative to liver (unitless)")                      # Table S2 KI__f_DAP2D3G = 10.0 (at its upper bound of 10)
    f_ugt1a9       <- fixed(1);       label("UGT1A9 activity scaling factor (unitless)")                               # SBML f_ugt1a9, normal = 1

    # ================= Renal excretion ===================================
    dapex_k <- fixed(0.01815); label("Urinary excretion rate constant, dapagliflozin (1/min)")  # Table S2 KI__DAPEX_k = 0.01815
    d3gex_k <- fixed(0.45036); label("Urinary excretion rate constant, D3OG (1/min)")           # Table S2 KI__D3GEX_k = 0.45036

    # ================= Pharmacodynamics ==================================
    # Supplementary Table S3 (optimized). These are the published values;
    # the archived SBML v0.9.8 carries a DIFFERENT pharmacodynamic set. The
    # paper's own Figure S19 shows the renal threshold for glucose starting
    # at 8.0 mM in a healthy subject, which identifies Table S3 (rtg_base =
    # 8.00) as the set that generated the published figures. See the
    # vignette Errata.
    gfr_healthy       <- fixed(100);     label("Glomerular filtration rate at normal renal function (mL/min)")  # SBML KI__GFR_healthy
    fpg_healthy       <- fixed(5);       label("Reference healthy fasting plasma glucose (mmol/L)")             # SBML KI__fpg_healthy; Methods 2.3 assumes 5 mM for healthy subjects
    rtg_base          <- fixed(8);       label("Baseline renal threshold for glucose at the healthy reference FPG (mmol/L)")  # Table S3 KI__RTG_base = 8.00
    rtg_m_fpg         <- fixed(1.2533);  label("Slope of the renal threshold for glucose on FPG (unitless)")    # Table S3 KI__RTG_m_fpg = 1.2533
    rtg_e50           <- fixed(6.49e-6); label("Kidney plasma dapagliflozin giving half-maximal RTG reduction (mmol/L)")  # Table S3 KI__RTG_E50 = 6.49e-6
    rtg_gamma         <- fixed(1.036);   label("Hill coefficient for the RTG reduction (unitless)")             # Table S3 KI__RTG_gamma = 1.036 (tabulated with unit 'mM'; a Hill coefficient is dimensionless - see vignette Errata)
    rtg_max_inhib     <- fixed(0.70673); label("Maximum fractional reduction of RTG (unitless)")                # Table S3 KI__RTG_max_inhibition = 0.70673

    # ================= Intravenous input =================================
    ti_dap <- fixed(10); label("Nominal intravenous injection time (s)")  # SBML ti_dap; drives a first-order release with half-life ti_dap

    # ================= Residual error ====================================
    # Methods Section 2.2 states that no population variability was
    # included and all simulations were deterministic, so no residual-error
    # model is reported; fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported by the source)")  # Methods 2.2: deterministic typical-individual simulations
  })

  model({
    # ================= Scenario scalars ==================================
    # Renal function multiplies GFR and every renal transport, metabolism
    # and excretion term (Methods Section 2.2).
    f_renal <- RENALFUNC_REL

    # The source parameterises hepatic impairment as cirrhosis SEVERITY,
    # the complement of the canonical fraction-of-normal column.
    f_cirrhosis    <- 1 - HEPFUNC_REL
    f_shunts       <- f_cirrhosis
    f_tissue_loss  <- f_cirrhosis

    # Fed state cuts the absorption rate to 30% of fasted (Methods 2.3).
    f_absorption <- 1 - 0.7 * FED

    # ================= Volumes ===========================================
    fvre <- 1 - (fvgu + fvki + fvli + fvlu + fvve + fvar)

    vgu <- WT * fvgu
    vki <- WT * fvki
    vli <- WT * fvli
    vlu <- WT * fvlu
    vre <- WT * fvre

    # Arterial and venous blood pools are corrected for the blood held
    # within the organ volumes.
    var_  <- WT * fvar - (fvar / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vve   <- WT * fvve - (fvve / (fvar + fvve)) * WT * fblood * (1 - fvve - fvar)
    vpo   <- (1 - hct) * (WT * fvpo - (fvpo / (fvar + fvve + fvpo + fvhv)) * WT * fblood * (1 - (fvar + fvve + fvpo + fvhv)))
    vhv   <- (1 - hct) * (WT * fvhv - (fvhv / (fvar + fvve + fvpo + fvhv)) * WT * fblood * (1 - (fvar + fvve + fvpo + fvhv)))

    vgu_plasma <- vgu * fblood * (1 - hct)
    vki_plasma <- vki * fblood * (1 - hct)
    vli_plasma <- vli * fblood * (1 - hct)
    vlu_plasma <- vlu * fblood * (1 - hct)
    vre_plasma <- vre * fblood * (1 - hct)

    vki_tissue <- vki * (1 - fblood)
    vlu_tissue <- vlu * (1 - fblood)
    vre_tissue <- vre * (1 - fblood)
    # Cirrhosis destroys functional parenchyma, so the metabolising liver
    # tissue volume shrinks with severity.
    vli_tissue <- vli * (1 - f_tissue_loss) * (1 - fblood)

    # ================= Blood flows =======================================
    co <- fcardio * WT * covbw
    qc <- (co / 1000) * 60

    fqre <- 1 - (fqki + fqh)

    qgu <- qc * fqgu
    qh  <- qc * fqh
    qki <- qc * fqki
    qlu <- qc * fqlu
    qre <- qc * fqre
    qha <- qh - qgu
    qpo <- qgu

    # ================= Absorption ========================================
    dissolution <- (ka_dis_dap / 60) * depot / mr_dap
    absorption  <- f_absorption * dapabs_k * vgu * gut_lumen
    dapabs      <- f_dap_abs * absorption
    dapexc      <- (1 - f_dap_abs) * absorption

    # ================= Intravenous input =================================
    ki_dap <- (0.693 / ti_dap) * 60
    iv_dap <- ki_dap * depot_iv / mr_dap

    # ================= Hepatic processes =================================
    # Reversible Michaelis-Menten uptake (OAT3) and export, and
    # irreversible UGT1A9 glucuronidation (Supplementary S4 liver model).
    li_dapim <- (dapim_vmax / dapim_km_dap) * vli_tissue *
      (liver_plasma - liver) /
      (1 + liver_plasma / dapim_km_dap + liver / dapim_km_dap)

    li_dap2d3g <- f_ugt1a9 * dap2d3g_vmax * vli_tissue *
      liver / (liver + dap2d3g_km_dap)

    li_d3gex <- (d3gex_vmax / d3gex_km_d3og) * vli_tissue *
      (liver_d3og - liver_plasma_d3og) /
      (1 + liver_d3og / d3gex_km_d3og + liver_plasma_d3og / d3gex_km_d3og)

    # ================= Renal processes ===================================
    # Supplementary S4 kidney model. Note that the manuscript's printed
    # DAPEX and D3GEX equations have their rate constants transposed; the
    # correct pairing (each species with its own rate constant) is used
    # here - see the vignette Errata.
    ki_dapim <- (f_renal * dapim_vmax / dapim_km_dap) * vki_tissue *
      (kidney_plasma - kidney) /
      (1 + kidney_plasma / dapim_km_dap + kidney / dapim_km_dap)

    ki_d3gim <- (f_renal * d3gim_vmax / d3gim_km_d3og) * vki_tissue *
      (kidney_plasma_d3og - kidney_d3og) /
      (1 + kidney_plasma_d3og / d3gim_km_d3og + kidney_d3og / d3gim_km_d3og)

    ki_dap2d3g <- f_renal * f_ugt1a9 * f_dap2d3g_ki * dap2d3g_vmax *
      vki_tissue * kidney / (kidney + dap2d3g_km_dap)

    ki_dapex <- f_renal * dapex_k * vki_tissue * kidney_plasma
    ki_d3gex <- f_renal * d3gex_k * vki_tissue * kidney_plasma_d3og

    # ================= Perfusion-limited tissue distribution =============
    # Lung and rest-of-body take up dapagliflozin only; Methods Section 2.2
    # states that no tissue partitioning was assumed for D3OG.
    transport_lu <- ftissue_dap * (lung_plasma * kp_dap - lung)
    transport_re <- ftissue_dap * (rest_plasma * kp_dap - rest)

    # ================= Convective flows ==================================
    flow_ar_gu <- qgu * arterial
    flow_ar_ki <- qki * arterial
    flow_ar_re <- qre * arterial
    flow_gu_po <- qgu * gut_plasma
    flow_hv_ve <- qh * hepatic_vein
    flow_ki_ve <- qki * kidney_plasma
    flow_lu_ar <- qlu * lung_plasma
    flow_re_ve <- qre * rest_plasma
    flow_ve_lu <- qlu * venous

    # Cirrhotic shunting diverts portal and hepatic-arterial blood straight
    # into the hepatic vein, bypassing the metabolising liver.
    flow_arli_hv <- f_shunts * qha * arterial
    flow_arli_li <- (1 - f_shunts) * qha * arterial
    flow_li_hv   <- (1 - f_shunts) * (qpo + qha) * liver_plasma
    flow_po_hv   <- f_shunts * qpo * portal
    flow_po_li   <- (1 - f_shunts) * qpo * portal

    flow_ar_gu_d3og <- qgu * arterial_d3og
    flow_ar_ki_d3og <- qki * arterial_d3og
    flow_ar_re_d3og <- qre * arterial_d3og
    flow_gu_po_d3og <- qgu * gut_plasma_d3og
    flow_hv_ve_d3og <- qh * hepatic_vein_d3og
    flow_ki_ve_d3og <- qki * kidney_plasma_d3og
    flow_lu_ar_d3og <- qlu * lung_plasma_d3og
    flow_re_ve_d3og <- qre * rest_plasma_d3og
    flow_ve_lu_d3og <- qlu * venous_d3og

    flow_arli_hv_d3og <- f_shunts * qha * arterial_d3og
    flow_arli_li_d3og <- (1 - f_shunts) * qha * arterial_d3og
    flow_li_hv_d3og   <- (1 - f_shunts) * (qpo + qha) * liver_plasma_d3og
    flow_po_hv_d3og   <- f_shunts * qpo * portal_d3og
    flow_po_li_d3og   <- (1 - f_shunts) * qpo * portal_d3og

    # ================= Pharmacodynamics ==================================
    # Supplementary S4 pharmacodynamics model. Renal function scales the
    # filtered glucose load as well as drug delivery to the kidney.
    gfr <- f_renal * gfr_healthy

    rtg_fpg   <- rtg_base + rtg_m_fpg * (FPG - fpg_healthy)
    rtg_delta <- rtg_fpg * rtg_max_inhib

    # The Emax term is asymmetric in the source: the numerator is driven by
    # kidney PLASMA dapagliflozin and the denominator by kidney TISSUE
    # dapagliflozin. Both the Supplementary S4 equation and the archived
    # SBML agree on this, so it is reproduced faithfully; because transport
    # into kidney tissue is fast and non-rate-limiting the two are nearly
    # equal in practice. See the vignette Errata.
    #
    # Both bases are clamped at zero before being raised to the
    # non-integer Hill exponent. Both states start at exactly 0, and a
    # solver step that takes either infinitesimally negative makes
    # `(-1e-20)^1.036` NaN, which poisons every state from the first step
    # onward - this is what happens for moderate and severe cirrhosis
    # without the guard. Concentrations are physically non-negative, so the
    # clamp changes nothing mathematically; it only removes the numerical
    # failure mode.
    rtg <- rtg_fpg - rtg_delta * max(kidney_plasma, 0)^rtg_gamma /
      (rtg_e50^rtg_gamma + max(kidney, 0)^rtg_gamma)

    # Glucose is excreted only above the renal threshold. Written as a
    # positive part rather than a branch so the right-hand side stays
    # continuous across the threshold crossings.
    glcex <- (gfr / cf_ml_l) * max(FPG - rtg, 0)

    # ================= ODEs ==============================================
    # Dosing compartments (mg). Oral drug dissolves into the gut lumen;
    # the intravenous compartment empties by a fast first-order release.
    d/dt(depot)    <- -dissolution * mr_dap
    d/dt(depot_iv) <- -iv_dap * mr_dap

    # Intestine (Supplementary S4 intestine model).
    d/dt(gut_lumen) <- -dapabs / vgu - dapexc / vgu + dissolution / vgu
    d/dt(feces)     <- dapexc

    # Gut plasma receives the absorbed fraction and drains to the portal vein.
    d/dt(gut_plasma)      <- flow_ar_gu / vgu_plasma - flow_gu_po / vgu_plasma +
      dapabs / vgu_plasma
    d/dt(gut_plasma_d3og) <- flow_ar_gu_d3og / vgu_plasma - flow_gu_po_d3og / vgu_plasma

    d/dt(portal)      <- -flow_po_li / vpo - flow_po_hv / vpo + flow_gu_po / vpo
    d/dt(portal_d3og) <- -flow_po_li_d3og / vpo - flow_po_hv_d3og / vpo +
      flow_gu_po_d3og / vpo

    # Liver (Supplementary S4 liver model).
    d/dt(liver_plasma)      <- flow_arli_li / vli_plasma + flow_po_li / vli_plasma -
      flow_li_hv / vli_plasma - li_dapim / vli_plasma
    d/dt(liver_plasma_d3og) <- flow_arli_li_d3og / vli_plasma +
      flow_po_li_d3og / vli_plasma - flow_li_hv_d3og / vli_plasma +
      li_d3gex / vli_plasma

    d/dt(liver)      <- li_dapim / vli_tissue - li_dap2d3g / vli_tissue
    d/dt(liver_d3og) <- li_dap2d3g / vli_tissue - li_d3gex / vli_tissue

    d/dt(hepatic_vein)      <- flow_arli_hv / vhv + flow_po_hv / vhv +
      flow_li_hv / vhv - flow_hv_ve / vhv
    d/dt(hepatic_vein_d3og) <- flow_arli_hv_d3og / vhv + flow_po_hv_d3og / vhv +
      flow_li_hv_d3og / vhv - flow_hv_ve_d3og / vhv

    # Kidney (Supplementary S4 kidney model).
    d/dt(kidney_plasma)      <- flow_ar_ki / vki_plasma - flow_ki_ve / vki_plasma -
      ki_dapim / vki_plasma - ki_dapex / vki_plasma
    d/dt(kidney_plasma_d3og) <- flow_ar_ki_d3og / vki_plasma -
      flow_ki_ve_d3og / vki_plasma - ki_d3gim / vki_plasma - ki_d3gex / vki_plasma

    d/dt(kidney)      <- ki_dapim / vki_tissue - ki_dap2d3g / vki_tissue
    d/dt(kidney_d3og) <- ki_d3gim / vki_tissue + ki_dap2d3g / vki_tissue

    d/dt(urine)      <- ki_dapex
    d/dt(urine_d3og) <- ki_d3gex

    # Lung sits between the venous and arterial pools.
    d/dt(lung_plasma)      <- -transport_lu / vlu_plasma + flow_ve_lu / vlu_plasma -
      flow_lu_ar / vlu_plasma
    d/dt(lung_plasma_d3og) <- flow_ve_lu_d3og / vlu_plasma - flow_lu_ar_d3og / vlu_plasma
    d/dt(lung)             <- transport_lu / vlu_tissue

    # Rest of body.
    d/dt(rest_plasma)      <- -transport_re / vre_plasma + flow_ar_re / vre_plasma -
      flow_re_ve / vre_plasma
    d/dt(rest_plasma_d3og) <- flow_ar_re_d3og / vre_plasma - flow_re_ve_d3og / vre_plasma
    d/dt(rest)             <- transport_re / vre_tissue

    # Systemic circulation.
    d/dt(venous)      <- iv_dap / vve + flow_ki_ve / vve + flow_hv_ve / vve -
      flow_ve_lu / vve + flow_re_ve / vve
    d/dt(venous_d3og) <- flow_ki_ve_d3og / vve + flow_hv_ve_d3og / vve -
      flow_ve_lu_d3og / vve + flow_re_ve_d3og / vve

    d/dt(arterial)      <- -flow_ar_ki / var_ - flow_arli_li / var_ -
      flow_arli_hv / var_ + flow_lu_ar / var_ - flow_ar_gu / var_ - flow_ar_re / var_
    d/dt(arterial_d3og) <- -flow_ar_ki_d3og / var_ - flow_arli_li_d3og / var_ -
      flow_arli_hv_d3og / var_ + flow_lu_ar_d3og / var_ - flow_ar_gu_d3og / var_ -
      flow_ar_re_d3og / var_

    # Urinary glucose.
    d/dt(glc_urine) <- glcex

    # ================= Observations ======================================
    # Venous plasma concentrations in umol/L, matching the paper's figures.
    Cc      <- 1000 * venous
    Cc_d3og <- 1000 * venous_d3og

    # Cumulative urinary and fecal recovery (mmol) and the pharmacodynamic
    # readouts: renal threshold for glucose (mmol/L) and cumulative urinary
    # glucose excretion (g).
    Aurine_dap  <- urine
    Aurine_d3og <- urine_d3og
    Afeces_dap  <- feces
    RTG         <- rtg
    UGE         <- glc_urine * mr_glc / cf_mg_g

    Cc ~ prop(propSd)
  })
}
