Granda_2024_tenofovir_pbpk <- function() {
  description <- paste(
    "PBPK (mechanistic kidney, 35 states).",
    "Individualised prediction of tenofovir renal clearance from",
    "per-subject biomarker measurements of kidney blood flow and",
    "OAT1/3-mediated tubular secretory capacity, on top of",
    "tracer-measured GFR (Granda 2024, Clin Transl Sci).",
    "The kidney is resolved into 11 longitudinal subsegments (proximal",
    "tubule S1-S3, descending and ascending loop of Henle, distal",
    "tubule, and five collecting-duct subsegments), each with a tubular",
    "lumen, a tubular epithelial cell and a peritubular blood subspace",
    "(33 states), plus a systemic blood compartment and a bladder.",
    "Mechanisms: unbound glomerular filtration, OAT1/3-mediated active",
    "secretion in the proximal tubule only, pH-dependent bidirectional",
    "passive diffusion along the whole nephron, and CKD tubular-flow",
    "adaptation. The structural equations, segment volumes, surface",
    "areas, tubular pH profile and adaptive flow factors are the",
    "Huang & Isoherranen 2018 / 2020 framework; Granda 2024 replaces the",
    "framework population defaults for kidney blood flow and secretory",
    "clearance with per-subject measured / fitted values.",
    "Predicts each subject's tenofovir renal clearance from their",
    "measured GFR and kidney blood flow plus the kynurenic-acid-derived",
    "secretory capacity, scaled by 0.033. Reproduces the Predicted TFV",
    "CLr column of Table 4.",
    "Deterministic typical-value model: the paper reports no IIV and no",
    "residual error.",
    sep = " "
  )
  reference <- paste(
    "Granda ML, Huang W, Yeung CK, Isoherranen N, Kestenbaum B.",
    "Predicting complex kidney drug handling using a physiologically-based",
    "pharmacokinetic model informed by biomarker-estimated secretory",
    "clearance and blood flow. Clin Transl Sci. 2024;17(1):e13678.",
    "doi:10.1111/cts.13678.",
    "Structural framework and physiological constants:",
    "Huang W, Isoherranen N. Development of a dynamic physiologically based",
    "mechanistic kidney model to predict renal clearance.",
    "CPT Pharmacometrics Syst Pharmacol. 2018;7(9):593-602.",
    "doi:10.1002/psp4.12321;",
    "Huang W, Isoherranen N. Novel mechanistic PBPK model to predict renal",
    "clearance in varying stages of CKD by incorporating tubular adaptation",
    "and dynamic passive reabsorption.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(10):571-583.",
    "doi:10.1002/psp4.12553 (adaptive tubular flows, segment volumes and",
    "surface areas, and the distributed MATLAB / Simulink implementation",
    "PSP4-9-571-s002 from which the unit conversions are taken).",
    sep = " "
  )
  vignette <- "Granda_2024_kidney_pbpk"

  # The 33 kidney states are the paper-mechanistic longitudinal
  # subsegments of the Huang & Isoherranen mechanistic kidney model:
  # <subsegment>_<lumen|cell|blood>. None corresponds to a canonical
  # compartment role; `central` and `urine` are canonical.
  paper_specific_compartment_pattern <-
    "^(pt1|pt2|pt3|lohdesc|lohasc|distal|cd1|cd2|cd3|cd4|cd5)_(lumen|cell|blood)$"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Iohexol renal clearance (measured GFR by an exogenous filtration",
        "tracer), RAW and NOT BSA-normalized. Table 3 tabulates the",
        "absolute per-subject value in mL/min (range 20-146); Table 1",
        "separately reports the BSA-normalized cohort summary of",
        "76 +/- 30 mL/min/1.73 m^2. The mechanistic kidney model consumes",
        "the ABSOLUTE whole-organ filtration flow, so the raw Table 3",
        "column is the correct input."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters twice: as the glomerular filtration flow qgfr = CRCL * 0.06",
        "(L/h) and as the tubular-adaptation driver pr = CRCL / 120, where",
        "120 mL/min is the framework healthy-adult anchor. Seven subjects",
        "have CRCL above 120, so pr exceeds 1 for them; that is the",
        "framework behaviour and is not clamped."
      ),
      source_name        = "Iohexol CLr"
    ),
    KBF = list(
      description        = paste(
        "Kidney blood flow, estimated per subject as the measured renal",
        "clearance of isovalerylglycine -- an endogenous solute chosen for",
        "its rapid clearance (more than fourfold GFR), low protein binding",
        "and high extraction ratio. Range 76-1692 mL/min (Table 3)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as the peritubular vascular flow qk = KBF * 0.06 (L/h),",
        "REPLACING the framework population default Q_kidney = 60 * pr L/h.",
        "Substituting this per-subject measurement is the paper's central",
        "methodological contribution. The model enforces the physiological",
        "constraint that secretory clearance cannot exceed KBF, so a",
        "subject whose measured biomarker clearance exceeds their KBF",
        "cannot be fitted (participants 38, 41 and 45)."
      ),
      source_name        = "Isovalerylglycine CLr"
    )
  )

  compartmentData <- list(
    pt1_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, proximal tubule S1
    pt1_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, proximal tubule S1
    pt1_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, proximal tubule S1
    pt2_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, proximal tubule S2
    pt2_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, proximal tubule S2
    pt2_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, proximal tubule S2
    pt3_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, proximal tubule S3
    pt3_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, proximal tubule S3
    pt3_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, proximal tubule S3
    lohdesc_lumen    = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, descending loop of Henle
    lohdesc_cell     = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, descending loop of Henle
    lohdesc_blood    = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, descending loop of Henle
    lohasc_lumen     = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, ascending loop of Henle
    lohasc_cell      = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, ascending loop of Henle
    lohasc_blood     = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, ascending loop of Henle
    distal_lumen     = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, distal tubule
    distal_cell      = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, distal tubule
    distal_blood     = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, distal tubule
    cd1_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, connecting tubule
    cd1_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, connecting tubule
    cd1_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, connecting tubule
    cd2_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, initial collecting duct
    cd2_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, initial collecting duct
    cd2_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, initial collecting duct
    cd3_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, cortical collecting duct
    cd3_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, cortical collecting duct
    cd3_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, cortical collecting duct
    cd4_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, medullary collecting duct
    cd4_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, medullary collecting duct
    cd4_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, medullary collecting duct
    cd5_lumen        = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE),  # tubular fluid, papillary collecting duct
    cd5_cell         = list(analyte = "Tenofovir", units = "mg", specimen = "tissue", verified = TRUE),  # epithelial cell, papillary collecting duct
    cd5_blood        = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # peritubular blood, papillary collecting duct
    central          = list(analyte = "Tenofovir", units = "mg", specimen = "whole blood", verified = TRUE),  # systemic venous blood (V_ven0 = 42 L)
    urine            = list(analyte = "Tenofovir", units = "mg", specimen = "urine", verified = TRUE)  # cumulative amount excreted into bladder urine
  )

  population <- list(
    species        = "human",
    n_subjects     = 27,
    n_studies      = 1,
    age_mean       = "54 years (SD 15)",
    sex_female_pct = 37,
    race_ethnicity = c(Black = 41, White = 56, Other = 3.7),
    bmi_mean       = "29.0 kg/m^2 (SD 5.7)",
    disease_state  = paste(
      "Adult outpatients spanning the full range of kidney function",
      "(CKD stages 1-5, not on renal replacement therapy). iGFR",
      "76 +/- 30 mL/min/1.73 m^2: 6 subjects (22%) below 45, 4 (15%)",
      "45-60, 9 (33%) 60-90, 8 (30%) above 90. Exclusions: nephrotic",
      "syndrome, cirrhosis, renal replacement therapy."
    ),
    renal_function = paste(
      "Per-subject iohexol CLr 20-146 mL/min and isovalerylglycine CLr",
      "76-1692 mL/min (Table 3); kynurenic acid CLr 226 +/- 138 mL/min",
      "and creatinine clearance 105 +/- 50 mL/min (Table 1)."
    ),
    dose_range     = paste(
      "Single oral tenofovir alafenamide 50 mg plus single oral",
      "oseltamivir 30 mg, with 13 plasma samples over 600 min and a",
      "concurrent 10-h timed urine collection."
    ),
    regions        = "United States (University of Washington, Seattle)",
    notes          = paste(
      "Ancillary study to PROCLAIM (Proximal Tubular Clearance of Renal",
      "Medications); 54 subjects completed the primary visit and 27",
      "completed the second visit analysed here. Baseline demographics",
      "are Table 1; per-subject model inputs are Table 3."
    )
  )

  ini({
    ## ---- The one fitted parameter -------------------------------------
    ## Unbound intrinsic secretory clearance of kynurenic acid, TOTAL
    ## across the three proximal-tubule subsegments (the model divides by
    ## 3 internally). Granda fitted this per subject over a 10-1000 L/h
    ## grid; the typical value here is the median of the 27 optimized
    ## values in Table 3. Supply per-subject values via rxSolve(params=).
    lclintsec <- log(220)
    label("Unbound intrinsic secretory clearance, total over the 3 proximal-tubule subsegments (L/h)")  # Granda 2024 Table 3, median of the 27 optimized values

    ## ---- Drug-specific in vitro inputs (Granda 2024 Table 2) ----------
    fup  <- fixed(0.993)
    label("Unbound fraction in plasma")  # Granda 2024 Table 2 (footnote b) and Methods ("set as 99.3% ... for all subjects based on literature")
    bp   <- fixed(0.55)
    label("Blood-to-plasma concentration ratio")  # Granda 2024 Table 2 (footnote e)
    papp <- fixed(0.37)
    label("Apparent passive permeability (10^-6 cm/s)")  # Granda 2024 Table 2 (footnote e)
    secscalar <- fixed(0.033)
    label("Biomarker-to-drug secretory-clearance scalar (unitless)")  # Granda 2024 Results, Mechanistic prediction section: "We calculated optimized scalars of 0.033 for tenofovir"

    ## ---- Ionization ----------------------------------------------------
    ## NOT REPORTED by Granda 2024 (Table 2 lists MW, LogP, fu,p, B/P and
    ## permeability but no pKa), nor by Huang 2018/2020 or Chang 2023. The
    ## value below is MODEL-IDENTIFIED, not transcribed: reproducing the
    ## paper's own Tables 3 and 4 discriminates decisively between the two
    ## candidate readings (fully-ionized acid: 24/25 subjects within 5% for
    ## kynurenic acid and 27/27 within 5% for both drugs; neutral/
    ## non-ionizing: 13/25 and 1/27). Any pka_acid at or below ~4.5 gives
    ## numerically IDENTICAL predictions, because the un-ionized fraction is
    ## then below 0.1% across the tubular pH range 6.5-7.2 and the passive
    ## diffusion terms vanish; the identified quantity is therefore
    ## "un-ionized fraction ~ 0", not a specific pKa. pka_base = 1 is the
    ## framework idiom for "no basic moiety" (Huang 2020 MATLAB comment:
    ## "a very low number for the pKa_base --> no ionization").
    pka_acid <- fixed(2)
    label("Acid pKa (model-identified as fully ionized; see comment)")  # NOT reported in any on-disk source; identified by reproduction of Tables 3-4
    pka_base <- fixed(1)
    label("Base pKa (framework sentinel for no basic moiety)")  # Huang & Isoherranen 2020 MATLAB pKa_b = 1

    ## ---- Framework system parameters -----------------------------------
    lvc <- fixed(log(42))
    label("Central (venous blood) volume of distribution (L)")  # Huang & Isoherranen 2020 MATLAB V_ven0 = 42 L
    fu_invitro <- fixed(1)
    label("Un-ionized fraction in the transwell donor chamber (unitless)")  # Huang & Isoherranen 2020 MATLAB In_vitro_unionization = 1 (framework default; not reported by Granda 2024, and provably immaterial here because the passive terms vanish)
    clh <- fixed(0)
    label("Hepatic (non-renal) blood clearance (L/h)")  # Huang & Isoherranen 2020 MATLAB CLh = 0; Granda 2024 simulates renal clearance only
    clreabs_api <- fixed(0)
    label("Active apical reabsorption clearance (L/h)")  # Huang & Isoherranen 2020 MATLAB CL_api_reabs = 0; no active reabsorption reported for this compound
    clreabs_bsl <- fixed(0)
    label("Active basolateral reabsorption clearance (L/h)")  # Huang & Isoherranen 2020 MATLAB CL_bsl_reabs = 0
    clint_met <- fixed(0)
    label("Renal cellular metabolic intrinsic clearance (L/h)")  # Huang & Isoherranen 2020 MATLAB CL_kidney_intrinsic = 0; no renal metabolism reported for this compound

    ## ---- Residual error --------------------------------------------------
    ## Granda 2024 reports no residual-error or IIV estimates: the model is
    ## a deterministic per-subject forward predictor whose inputs are all
    ## measured. Fixed at zero rather than invented.
    propSd <- fixed(0)
    label("Proportional residual error (fraction; not estimated by the paper)")  # not reported; fixed at 0
  })

  model({
    ## =====================================================================
    ## 1. Kidney-function scaling. `pr` is the fraction of healthy kidney
    ##    function; 120 mL/min is the healthy-adult GFR anchor of the
    ##    framework (Huang & Isoherranen 2020 MATLAB: Proportional_Reduction
    ##    = GFR_test / 120). 1 mL/min = 0.06 L/h (framework `Unit = 7.2/120`).
    ##    CRCL is the RAW iohexol renal clearance (mL/min) and KBF the RAW
    ##    isovalerylglycine renal clearance (mL/min), both per subject.
    ## =====================================================================
    pr   <- CRCL / 120
    qgfr <- CRCL * 0.06   # glomerular filtration flow (L/h)
    qk   <- KBF * 0.06    # peritubular (kidney) blood flow (L/h)

    ## =====================================================================
    ## 2. Adaptive CKD tubular flows (L/h):
    ##      Q_i = Q_i,healthy * (1 - (1 - pr) * (1 - AF_i))
    ##    Healthy flows (mL/min) and adaptation factors AF_i are the
    ##    "For Adaptive Model" block of the Huang & Isoherranen 2020
    ##    distributed MATLAB code (PSP4-9-571-s002).
    ## =====================================================================
    q1       <-   120 * 0.06 * (1 - (1 - pr) * (1 - 0))  # proximal tubule S1
    q2       <-    94 * 0.06 * (1 - (1 - pr) * (1 - 0.009922))  # proximal tubule S2
    q3       <-    68 * 0.06 * (1 - (1 - pr) * (1 - 0.01753))  # proximal tubule S3
    q4       <-    43 * 0.06 * (1 - (1 - pr) * (1 - 0.03848))  # descending loop of Henle
    q5       <-    24 * 0.06 * (1 - (1 - pr) * (1 - 0.09769))  # ascending loop of Henle
    q6       <-    24 * 0.06 * (1 - (1 - pr) * (1 - 0.09769))  # distal tubule
    q7       <-    11 * 0.06 * (1 - (1 - pr) * (1 - 0.2606))  # connecting tubule
    q8       <-     9 * 0.06 * (1 - (1 - pr) * (1 - 0.3119))  # initial collecting duct
    q9       <-     7 * 0.06 * (1 - (1 - pr) * (1 - 0.3735))  # cortical collecting duct
    q10      <-     5 * 0.06 * (1 - (1 - pr) * (1 - 0.4428))  # medullary collecting duct
    q11      <-     3 * 0.06 * (1 - (1 - pr) * (1 - 0.5114))  # papillary collecting duct
    q_urine  <-     1 * 0.06 * (1 - (1 - pr) * (1 - 0.5611))  # urine formation

    ## =====================================================================
    ## 3. Segment volumes (L) and luminal surface areas (dm^2), each scaled
    ##    proportionally with residual kidney function. Huang & Isoherranen
    ##    2020 MATLAB V_PT/V_HT/V_DT/V_CDT and S_PT/S_HT/S_DT/S_CD. The
    ##    lumen, cell and blood subspaces of a segment share its volume.
    ## =====================================================================
    v_pt     <- 0.0305 * pr
    v_loh    <- 0.0027 * pr
    v_distal <- 0.0194 * pr
    v_cd     <- 0.0090 * pr
    s_pt     <- 6107 * pr
    s_loh    <-   61 * pr
    s_distal <-  156 * pr
    s_cd     <-  6.7 * pr

    ## =====================================================================
    ## 4. Passive diffusion clearances (L/h). `papp` is the in vitro apparent
    ##    permeability in 10^-6 cm/s; dividing by the in vitro un-ionized
    ##    donor fraction and by the 1.5 microvilli factor gives the intrinsic
    ##    permeability, and 0.00036 converts 10^-6 cm/s to dm/h so that
    ##    dm^2 * dm/h = L/h. Basolateral surface area of the proximal tubule
    ##    is 1/30 of the apical (brush-border) area.
    ## =====================================================================
    pint  <- papp / fu_invitro / 1.5
    pdmhr <- 0.00036 * pint
    clpd_pt_api <- s_pt * pdmhr
    clpd_pt_bsl <- s_pt / 30 * pdmhr
    clpd_loh    <- s_loh * pdmhr
    clpd_distal <- s_distal * pdmhr
    clpd_cd     <- s_cd * pdmhr

    ## =====================================================================
    ## 5. Un-ionized fractions along the tubular pH gradient (pH 7.2 -> 6.5),
    ##    in the epithelial cell (pH 7.2) and in blood (pH 7.4). Product of
    ##    the acid and base Henderson-Hasselbalch terms, exactly the
    ##    alpha / beta / gamma quantities of the framework MATLAB code.
    ## =====================================================================
    fnu_t1  <- 1 / (1 + 10^(7.2 - pka_acid)) * 1 / (1 + 10^(pka_base - 7.2))
    fnu_t2  <- 1 / (1 + 10^(7.1 - pka_acid)) * 1 / (1 + 10^(pka_base - 7.1))
    fnu_t3  <- 1 / (1 + 10^(7 - pka_acid)) * 1 / (1 + 10^(pka_base - 7))
    fnu_t4  <- 1 / (1 + 10^(7 - pka_acid)) * 1 / (1 + 10^(pka_base - 7))
    fnu_t5  <- 1 / (1 + 10^(7 - pka_acid)) * 1 / (1 + 10^(pka_base - 7))
    fnu_t6  <- 1 / (1 + 10^(6.9 - pka_acid)) * 1 / (1 + 10^(pka_base - 6.9))
    fnu_t7  <- 1 / (1 + 10^(6.8 - pka_acid)) * 1 / (1 + 10^(pka_base - 6.8))
    fnu_t8  <- 1 / (1 + 10^(6.7 - pka_acid)) * 1 / (1 + 10^(pka_base - 6.7))
    fnu_t9  <- 1 / (1 + 10^(6.6 - pka_acid)) * 1 / (1 + 10^(pka_base - 6.6))
    fnu_t10 <- 1 / (1 + 10^(6.5 - pka_acid)) * 1 / (1 + 10^(pka_base - 6.5))
    fnu_t11 <- 1 / (1 + 10^(6.5 - pka_acid)) * 1 / (1 + 10^(pka_base - 6.5))
    fnu_cell  <- 1 / (1 + 10^(7.2 - pka_acid)) * 1 / (1 + 10^(pka_base - 7.2))
    fnu_blood <- 1 / (1 + 10^(7.4 - pka_acid)) * 1 / (1 + 10^(pka_base - 7.4))

    ## =====================================================================
    ## 6. Binding, secretion and the vascular inlet.
    ##    `fupbp` is the unbound fraction referenced to BLOOD; the framework
    ##    sets the intracellular unbound fraction equal to the plasma one
    ##    (MATLAB `fu_k = fu0`). `clsec` is the PER-SUBSEGMENT unbound
    ##    intrinsic secretory clearance: the tabulated value is the total
    ##    across the three proximal-tubule subsegments, and the same value is
    ##    applied to basolateral uptake and to apical efflux.
    ## =====================================================================
    vc     <- exp(lvc)
    fupbp  <- fup / bp
    fucell <- fup
    clsec  <- exp(lclintsec) * secscalar / 3

    cb      <- central / vc                      # blood concentration (mg/L)
    cbowman <- fupbp * cb                        # Bowman's capsule filtrate
    cpb0    <- (1 - fupbp * qgfr / qk) * cb      # post-filtration peritubular inlet

    ## ---- section concentrations (mg/L) ----------------------------------
    c_pt1_lumen <- pt1_lumen / v_pt
    c_pt1_cell <- pt1_cell / v_pt
    c_pt1_blood <- pt1_blood / v_pt
    c_pt2_lumen <- pt2_lumen / v_pt
    c_pt2_cell <- pt2_cell / v_pt
    c_pt2_blood <- pt2_blood / v_pt
    c_pt3_lumen <- pt3_lumen / v_pt
    c_pt3_cell <- pt3_cell / v_pt
    c_pt3_blood <- pt3_blood / v_pt
    c_lohdesc_lumen <- lohdesc_lumen / v_loh
    c_lohdesc_cell <- lohdesc_cell / v_loh
    c_lohdesc_blood <- lohdesc_blood / v_loh
    c_lohasc_lumen <- lohasc_lumen / v_loh
    c_lohasc_cell <- lohasc_cell / v_loh
    c_lohasc_blood <- lohasc_blood / v_loh
    c_distal_lumen <- distal_lumen / v_distal
    c_distal_cell <- distal_cell / v_distal
    c_distal_blood <- distal_blood / v_distal
    c_cd1_lumen <- cd1_lumen / v_cd
    c_cd1_cell <- cd1_cell / v_cd
    c_cd1_blood <- cd1_blood / v_cd
    c_cd2_lumen <- cd2_lumen / v_cd
    c_cd2_cell <- cd2_cell / v_cd
    c_cd2_blood <- cd2_blood / v_cd
    c_cd3_lumen <- cd3_lumen / v_cd
    c_cd3_cell <- cd3_cell / v_cd
    c_cd3_blood <- cd3_blood / v_cd
    c_cd4_lumen <- cd4_lumen / v_cd
    c_cd4_cell <- cd4_cell / v_cd
    c_cd4_blood <- cd4_blood / v_cd
    c_cd5_lumen <- cd5_lumen / v_cd
    c_cd5_cell <- cd5_cell / v_cd
    c_cd5_blood <- cd5_blood / v_cd

    ## =====================================================================
    ## 7. The 33 kidney ODEs: 11 longitudinal subsegments x (lumen, cell,
    ##    peritubular blood). Active secretion (`clsec`) and active
    ##    reabsorption act only in the three proximal-tubule subsegments.
    ## =====================================================================
    ## proximal tubule S1 (subsegment 1 of 11).
    d/dt(pt1_lumen) <- q1 * cbowman - q2 * c_pt1_lumen +
      clsec * c_pt1_cell * fucell - clreabs_api * c_pt1_lumen +
      clpd_pt_api * (c_pt1_cell * fucell * fnu_cell - c_pt1_lumen * fnu_t1)
    d/dt(pt1_cell) <- clsec * c_pt1_blood * fupbp - clsec * c_pt1_cell * fucell +
      clreabs_api * c_pt1_lumen - clreabs_bsl * c_pt1_cell * fucell +
      clpd_pt_bsl * (c_pt1_blood * fupbp * fnu_blood - c_pt1_cell * fucell * fnu_cell) +
      clpd_pt_api * (c_pt1_lumen * fnu_t1 - c_pt1_cell * fucell * fnu_cell) -
      clint_met * c_pt1_cell * fucell
    d/dt(pt1_blood) <- qk * (cpb0 - c_pt1_blood) -
      clsec * c_pt1_blood * fupbp + clreabs_bsl * c_pt1_cell * fucell +
      clpd_pt_bsl * (c_pt1_cell * fucell * fnu_cell - c_pt1_blood * fupbp * fnu_blood)

    ## proximal tubule S2 (subsegment 2 of 11).
    d/dt(pt2_lumen) <- q2 * c_pt1_lumen - q3 * c_pt2_lumen +
      clsec * c_pt2_cell * fucell - clreabs_api * c_pt2_lumen +
      clpd_pt_api * (c_pt2_cell * fucell * fnu_cell - c_pt2_lumen * fnu_t2)
    d/dt(pt2_cell) <- clsec * c_pt2_blood * fupbp - clsec * c_pt2_cell * fucell +
      clreabs_api * c_pt2_lumen - clreabs_bsl * c_pt2_cell * fucell +
      clpd_pt_bsl * (c_pt2_blood * fupbp * fnu_blood - c_pt2_cell * fucell * fnu_cell) +
      clpd_pt_api * (c_pt2_lumen * fnu_t2 - c_pt2_cell * fucell * fnu_cell) -
      clint_met * c_pt2_cell * fucell
    d/dt(pt2_blood) <- qk * (c_pt1_blood - c_pt2_blood) -
      clsec * c_pt2_blood * fupbp + clreabs_bsl * c_pt2_cell * fucell +
      clpd_pt_bsl * (c_pt2_cell * fucell * fnu_cell - c_pt2_blood * fupbp * fnu_blood)

    ## proximal tubule S3 (subsegment 3 of 11).
    d/dt(pt3_lumen) <- q3 * c_pt2_lumen - q4 * c_pt3_lumen +
      clsec * c_pt3_cell * fucell - clreabs_api * c_pt3_lumen +
      clpd_pt_api * (c_pt3_cell * fucell * fnu_cell - c_pt3_lumen * fnu_t3)
    d/dt(pt3_cell) <- clsec * c_pt3_blood * fupbp - clsec * c_pt3_cell * fucell +
      clreabs_api * c_pt3_lumen - clreabs_bsl * c_pt3_cell * fucell +
      clpd_pt_bsl * (c_pt3_blood * fupbp * fnu_blood - c_pt3_cell * fucell * fnu_cell) +
      clpd_pt_api * (c_pt3_lumen * fnu_t3 - c_pt3_cell * fucell * fnu_cell) -
      clint_met * c_pt3_cell * fucell
    d/dt(pt3_blood) <- qk * (c_pt2_blood - c_pt3_blood) -
      clsec * c_pt3_blood * fupbp + clreabs_bsl * c_pt3_cell * fucell +
      clpd_pt_bsl * (c_pt3_cell * fucell * fnu_cell - c_pt3_blood * fupbp * fnu_blood)

    ## descending loop of Henle (subsegment 4 of 11).
    d/dt(lohdesc_lumen) <- q4 * c_pt3_lumen - q5 * c_lohdesc_lumen +
      clpd_loh * (c_lohdesc_cell * fucell * fnu_cell - c_lohdesc_lumen * fnu_t4)
    d/dt(lohdesc_cell) <- clpd_loh * (c_lohdesc_blood * fupbp * fnu_blood + c_lohdesc_lumen * fnu_t4 -
      2 * c_lohdesc_cell * fucell * fnu_cell)
    d/dt(lohdesc_blood) <- qk * (c_pt3_blood - c_lohdesc_blood) +
      clpd_loh * (c_lohdesc_cell * fucell * fnu_cell - c_lohdesc_blood * fupbp * fnu_blood)

    ## ascending loop of Henle (subsegment 5 of 11).
    d/dt(lohasc_lumen) <- q5 * c_lohdesc_lumen - q6 * c_lohasc_lumen +
      clpd_loh * (c_lohasc_cell * fucell * fnu_cell - c_lohasc_lumen * fnu_t5)
    d/dt(lohasc_cell) <- clpd_loh * (c_lohasc_blood * fupbp * fnu_blood + c_lohasc_lumen * fnu_t5 -
      2 * c_lohasc_cell * fucell * fnu_cell)
    d/dt(lohasc_blood) <- qk * (c_lohdesc_blood - c_lohasc_blood) +
      clpd_loh * (c_lohasc_cell * fucell * fnu_cell - c_lohasc_blood * fupbp * fnu_blood)

    ## distal tubule (subsegment 6 of 11).
    d/dt(distal_lumen) <- q6 * c_lohasc_lumen - q7 * c_distal_lumen +
      clpd_distal * (c_distal_cell * fucell * fnu_cell - c_distal_lumen * fnu_t6)
    d/dt(distal_cell) <- clpd_distal * (c_distal_blood * fupbp * fnu_blood + c_distal_lumen * fnu_t6 -
      2 * c_distal_cell * fucell * fnu_cell)
    d/dt(distal_blood) <- qk * (c_lohasc_blood - c_distal_blood) +
      clpd_distal * (c_distal_cell * fucell * fnu_cell - c_distal_blood * fupbp * fnu_blood)

    ## connecting tubule (subsegment 7 of 11).
    d/dt(cd1_lumen) <- q7 * c_distal_lumen - q8 * c_cd1_lumen +
      clpd_cd * (c_cd1_cell * fucell * fnu_cell - c_cd1_lumen * fnu_t7)
    d/dt(cd1_cell) <- clpd_cd * (c_cd1_blood * fupbp * fnu_blood + c_cd1_lumen * fnu_t7 -
      2 * c_cd1_cell * fucell * fnu_cell)
    d/dt(cd1_blood) <- qk * (c_distal_blood - c_cd1_blood) +
      clpd_cd * (c_cd1_cell * fucell * fnu_cell - c_cd1_blood * fupbp * fnu_blood)

    ## initial collecting duct (subsegment 8 of 11).
    d/dt(cd2_lumen) <- q8 * c_cd1_lumen - q9 * c_cd2_lumen +
      clpd_cd * (c_cd2_cell * fucell * fnu_cell - c_cd2_lumen * fnu_t8)
    d/dt(cd2_cell) <- clpd_cd * (c_cd2_blood * fupbp * fnu_blood + c_cd2_lumen * fnu_t8 -
      2 * c_cd2_cell * fucell * fnu_cell)
    d/dt(cd2_blood) <- qk * (c_cd1_blood - c_cd2_blood) +
      clpd_cd * (c_cd2_cell * fucell * fnu_cell - c_cd2_blood * fupbp * fnu_blood)

    ## cortical collecting duct (subsegment 9 of 11).
    d/dt(cd3_lumen) <- q9 * c_cd2_lumen - q10 * c_cd3_lumen +
      clpd_cd * (c_cd3_cell * fucell * fnu_cell - c_cd3_lumen * fnu_t9)
    d/dt(cd3_cell) <- clpd_cd * (c_cd3_blood * fupbp * fnu_blood + c_cd3_lumen * fnu_t9 -
      2 * c_cd3_cell * fucell * fnu_cell)
    d/dt(cd3_blood) <- qk * (c_cd2_blood - c_cd3_blood) +
      clpd_cd * (c_cd3_cell * fucell * fnu_cell - c_cd3_blood * fupbp * fnu_blood)

    ## medullary collecting duct (subsegment 10 of 11).
    d/dt(cd4_lumen) <- q10 * c_cd3_lumen - q11 * c_cd4_lumen +
      clpd_cd * (c_cd4_cell * fucell * fnu_cell - c_cd4_lumen * fnu_t10)
    d/dt(cd4_cell) <- clpd_cd * (c_cd4_blood * fupbp * fnu_blood + c_cd4_lumen * fnu_t10 -
      2 * c_cd4_cell * fucell * fnu_cell)
    d/dt(cd4_blood) <- qk * (c_cd3_blood - c_cd4_blood) +
      clpd_cd * (c_cd4_cell * fucell * fnu_cell - c_cd4_blood * fupbp * fnu_blood)

    ## papillary collecting duct (subsegment 11 of 11).
    d/dt(cd5_lumen) <- q11 * c_cd4_lumen - q_urine * c_cd5_lumen +
      clpd_cd * (c_cd5_cell * fucell * fnu_cell - c_cd5_lumen * fnu_t11)
    d/dt(cd5_cell) <- clpd_cd * (c_cd5_blood * fupbp * fnu_blood + c_cd5_lumen * fnu_t11 -
      2 * c_cd5_cell * fucell * fnu_cell)
    d/dt(cd5_blood) <- qk * (c_cd4_blood - c_cd5_blood) +
      clpd_cd * (c_cd5_cell * fucell * fnu_cell - c_cd5_blood * fupbp * fnu_blood)

    ## =====================================================================
    ## 8. Systemic compartment and bladder.
    ## =====================================================================
    d/dt(central) <- qk * (c_cd5_blood - cb) - clh * cb / bp
    d/dt(urine)   <- q_urine * c_cd5_lumen

    ## =====================================================================
    ## 9. Outputs. `Cc` is the plasma concentration; `CLr` is the renal
    ##    clearance in mL/min (urinary excretion rate / plasma concentration,
    ##    times 1000/60 to convert L/h to mL/min), which is the quantity the
    ##    paper tabulates in Tables 3 and 4.
    ## =====================================================================
    Cc  <- cb / bp
    CLr <- 1000 / 60 * q_urine * c_cd5_lumen / Cc
    Cc ~ prop(propSd)
  })
}
