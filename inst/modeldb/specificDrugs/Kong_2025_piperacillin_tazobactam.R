Kong_2025_piperacillin_tazobactam <- function() {
  description <- paste(
    "Joint four-compartment population PK model for piperacillin and",
    "tazobactam in 20 adults with end-stage kidney disease receiving",
    "thrice-weekly high-flux intermittent haemodialysis (Kong 2025).",
    "The model is an extension of the authors' general-purpose",
    "piperacillin/tazobactam model: every structural parameter, the",
    "clearance maturation and decline functions, the serum-creatinine",
    "effect and all inter-individual variances are FIXED to the",
    "general-purpose values, and only four things are estimated from the",
    "ESKD data -- the residual endogenous clearance fraction (theta_ESKD,",
    "with residual diuresis as a covariate), the piperacillin dialyser",
    "extraction ratio, the tazobactam dialyser extraction ratio (three",
    "levels by vascular access type), and the residual error.",
    "Dialyser clearance is added to the central compartment of each drug",
    "as extraction ratio times blood flow rate, gated on by a",
    "time-varying during-session indicator. Both drugs are fitted",
    "simultaneously via the NONMEM L2 data item. Piperacillin uses the",
    "unsuffixed canonical compartment / parameter set; tazobactam carries",
    "the sibling-drug suffix _taz throughout."
  )
  reference <- paste(
    "Kong D, Koomen JV, Vanommeslaeghe F, Delanghe S, Van Biesen W,",
    "Colin PJ, Eloot S. A Population Pharmacokinetic Analysis for",
    "Piperacillin/Tazobactam in Patients with End-Stage Kidney Disease",
    "Undergoing Intermittent Haemodialysis: Extension of a",
    "General-Purpose Model. Clin Pharmacokinet. 2025.",
    "doi:10.1007/s40262-025-01527-y.",
    "Structural parameters, maturation / decline functions,",
    "serum-creatinine effect and inter-individual variances are fixed to",
    "the general-purpose parent model: Kong D, et al. A pooled",
    "pharmacokinetic analysis for piperacillin/tazobactam across",
    "different patient populations: from premature infants to the",
    "elderly. Clin Pharmacokinet. 2025;64(1):107-126. The parent-model",
    "values used here are transcribed from the final NONMEM control",
    "stream reproduced verbatim in the Kong 2025 Electronic",
    "Supplementary Material, so nothing is imported from the parent",
    "publication itself.",
    sep = " "
  )
  vignette <- "Kong_2025_piperacillin_tazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central         = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    central_taz     = list(analyte = "tazobactam",   units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_taz = list(analyte = "tazobactam",   units = "mg", specimen = "plasma", verified = TRUE)
  )

  dosing <- c("central", "central_taz")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scalar on every structural parameter, normalised to a",
        "70 kg reference individual (Kong 2025 Methods 2.3; ESM $PK",
        "'FSIZE = WT/70'). Conventional fixed West-Brown-Enquist exponents",
        "are used: 0.75 on the clearances and 1.00 on the volumes. The",
        "inter-compartmental clearances instead use COMPARTMENTAL",
        "allometry, scaling on the individual peripheral volume rather than",
        "on weight directly: Q = theta_Q * (V2_i / theta_V2)^0.75, where",
        "V2_i already carries both the weight scalar and eta_V2 (ESM $PK",
        "'V2SZ = V2/REFV2'; 'TVQ2 = EXP(THETA(4)) * (V2SZ**0.75)').",
        "Time-fixed per subject in the source analysis. Cohort median",
        "66.3 kg, range 27.0-93.8 kg (Table 1)."
      ),
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives the sigmoidal maturation and decline functions on the",
        "endogenous clearance of both drugs, and the postmenstrual-age",
        "standardisation of serum creatinine (Eq. 1). Kong 2025 computes",
        "postmenstrual age by adding an assumed gestational age of 40",
        "weeks to the patient's chronological age (Methods 2.4) and",
        "expresses it in YEARS in the model code (ESM $PK 'RPMA =",
        "35+(40/52)'). The canonical PAGE column carries MONTHS, so",
        "model() converts with PAGE / 12 before use -- the same",
        "reparameterisation pattern documented for Llanos-Paez 2017 and",
        "Zhao 2018. The reference individual is a 35-year-old (35 + 40/52",
        "= 35.769 postmenstrual years = 429.2 months). Cohort median",
        "chronological age 71.5 years, range 20-84 (Table 1); the",
        "maturation ratio is therefore ~1 throughout this all-adult",
        "cohort and only the decline term moves."
      ),
      source_name        = "PMA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters the endogenous clearance through the exponential",
        "correction F_SCR = exp(-theta_SCR * (CREAT - CREAT_std)) of Kong",
        "2025 Eq. 2, where CREAT_std is the postmenstrual-age-standardised",
        "serum creatinine of Eq. 1 and theta_SCR = 0.346 dL/mg is fixed",
        "from the parent general-purpose model. Units are mg/dL, pinned by",
        "the dL/mg units of theta_SCR (Table 3).",
        "IMPORTANT: Eq. 2 sets F_SCR = 1 for patients on intermittent",
        "haemodialysis, so CREAT is IGNORED whenever",
        "RRT_HEMODIAL_STATUS = 1 and only matters for the",
        "normal-kidney-function comparator arm the paper simulates. The",
        "source dataset encodes this with a CR = 0 sentinel that the",
        "control stream maps to CREAT_std (ESM $PK 'SCR = CR; IF(CR.EQ.0)",
        "SCR = STCR'), making F_SCR = 1 by construction; model() here",
        "implements the published Eq. 2 gate on RRT_HEMODIAL_STATUS",
        "directly instead of relying on the sentinel."
      ),
      source_name        = "CR"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Chronic intermittent-haemodialysis (end-stage kidney disease) status",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on chronic intermittent haemodialysis)",
      notes              = paste(
        "Subject-level switch between the two branches of Kong 2025",
        "Eq. 2 and Eq. 4. When 1, F_SCR = 1 and the estimated residual",
        "endogenous clearance fraction F_ESKD applies; when 0, F_SCR is",
        "the serum-creatinine correction and F_ESKD = 1, which recovers the",
        "parent general-purpose model for a patient with normal kidney",
        "function. Every subject in the fitted dataset has value 1 (the",
        "cohort is ESKD-only), so the 0 branch is not gated by any",
        "parameter estimated here -- it exists because Kong 2025 uses",
        "exactly that branch to simulate the normal-kidney-function",
        "comparator arm of Table 4 and Figs. 3-4.",
        "Distinct from RRT_HEMODIAL_ACTIVE, which is the time-varying",
        "within-session gate on the dialyser clearance."
      ),
      source_name        = "(implicit; the fitted dataset is ESKD-only)"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Haemodialysis-session-active indicator (time-varying)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic period)",
      notes              = paste(
        "Time-varying per-session gate: 1 while blood is circulating",
        "through the dialyser, 0 in the interdialytic interval (ESM $PK",
        "indicator IND_DIA). Multiplies the dialyser clearance arm",
        "CL_DIA = ER * blood flow rate of Eq. 5, so total clearance is the",
        "endogenous clearance alone between sessions and endogenous plus",
        "dialyser clearance during a session (Methods 2.4). The same",
        "indicator also selects the post-dialyser observable: sampled",
        "blood exiting the dialyser carries C * (1 - ER), and equals the",
        "inlet concentration when no session is running (ESM $ERROR).",
        "Paired with BFR exactly as in Liesenfeld_2013_dabigatran.R.",
        "Simulations in the source paper assume a thrice-weekly schedule",
        "with 4-hour sessions (Methods 2.5)."
      ),
      source_name        = "IND_DIA"
    ),
    BFR = list(
      description        = "Blood flow rate through the dialyser",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Recorded dialyser blood flow rate; enters the dialyser clearance",
        "of Eq. 5 directly as CL_DIA = ER * BFR after conversion to L/h",
        "(ESM $PK 'FLOW = IFLOW*60/1000'). Meaningful only while",
        "RRT_HEMODIAL_ACTIVE = 1; the value is sentinel in the",
        "interdialytic period because the arm is gated off. Cohort median",
        "225 mL/min, range 160-350 (Table 1).",
        "For sessions using a SINGLE-needle arteriovenous fistula, Kong",
        "2025 supplied a TIME-AVERAGED blood flow rate to correct for",
        "recirculation rather than the nominal pump setting (Methods 2.4),",
        "so a user reproducing the analysis should apply the same",
        "correction to BFR when VASCACC_AVF1N = 1. Blood flow rate was",
        "also screened as a continuous covariate on the logit-scale",
        "extraction ratio (Eq. 7) and was not retained (Methods 2.4;",
        "Table 2)."
      ),
      source_name        = "IFLOW"
    ),
    URINE_VOL_24H = list(
      description        = "24-hour residual diuresis (residual urine output)",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used as the binary preserved-diuresis gate (URINE_VOL_24H > 100)",
        "selecting between the two estimated values of theta_ESKD: 0.214",
        "in patients WITHOUT residual diuresis and 0.370 in patients WITH",
        "residual diuresis (Table 3). Retained in the final model after a",
        "re-evaluation that followed the merge of the two drug-specific",
        "theta_ESKD parameters into one (Table 2 model 11, dOFV -12.37,",
        "p < 0.001); the drug-specific versions of the same effect had",
        "both been dropped in backward elimination (models 7 and 8).",
        "PROVENANCE CAVEAT: Kong 2025 reports residual diuresis only as a",
        "binary present / absent status (11/20 present, Table 1) with no",
        "urine volumes, and the source column UVOL0123 is gated as",
        "IF(UVOL0123.NE.0) in the ESM. The 100 mL/24h threshold used here",
        "is the anuria cutoff carried by the canonical URINE_VOL_24H",
        "register entry and the Huppe_2023_fosfomycin.R precedent, NOT a",
        "value stated by Kong 2025. Any positive value above 100 mL/24h",
        "reproduces the paper's 'residual diuresis present' branch."
      ),
      source_name        = "UVOL0123"
    ),
    VASCACC_AVF1N = list(
      description        = "Single-needle arteriovenous fistula vascular access",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tunnelled dialysis catheter, shared with VASCACC_AVF2N)",
      notes              = paste(
        "Selects the tazobactam dialyser extraction ratio: 73.9%",
        "(95% CI 65.3-81.0) for a single-needle arteriovenous fistula",
        "versus the 80.1% tunnelled-dialysis-catheter reference (Table 3).",
        "Together with VASCACC_AVF2N this pair decomposes the paper's",
        "three-level access classification (source column CVC0SN1DN2,",
        "0 = TDC / 1 = single-needle AVF / 2 = double-needle AVF) into two",
        "binary indicators with TDC as the implicit reference, matching the",
        "NONMEM default branch (ESM $PK 'LGT_TAZ = THETA(18);",
        "IF(CVC0SN1DN2.EQ.1) LGT_TAZ = THETA(19)').",
        "This effect applies to TAZOBACTAM ONLY. The analogous access-type",
        "effect on the piperacillin extraction ratio was added in forward",
        "inclusion (Table 2 model 5, dOFV -8.02, p = 0.0183) but removed",
        "in the more conservative backward-elimination step (model 9,",
        "dOFV +11.26, p < 0.001), so piperacillin keeps a single 64.0%",
        "extraction ratio. Cohort split 6 (32%) TDC / 4 (21%) AVF 1N /",
        "9 (47%) AVF 2N (Table 1)."
      ),
      source_name        = "CVC0SN1DN2 (level 1)"
    ),
    VASCACC_AVF2N = list(
      description        = "Double-needle arteriovenous fistula vascular access",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (tunnelled dialysis catheter, shared with VASCACC_AVF1N)",
      notes              = paste(
        "Selects the tazobactam dialyser extraction ratio: 73.5%",
        "(95% CI 65.0-80.5) for a double-needle arteriovenous fistula",
        "versus the 80.1% tunnelled-dialysis-catheter reference (Table 3;",
        "ESM $PK 'IF(CVC0SN1DN2.EQ.2) LGT_TAZ = THETA(20)').",
        "Mutually exclusive with VASCACC_AVF1N -- setting both to 1 is",
        "invalid. Kong 2025 judged the access-type differences",
        "statistically significant but clinically limited (Discussion)."
      ),
      source_name        = "CVC0SN1DN2 (level 2)"
    )
  )

  covariatesDataExcluded <- list(
    DIALYSIS_MODE_HDF = list(
      description        = "Dialysis modality (haemodialysis vs pre- or post-dilution haemodiafiltration)",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "haemodialysis",
      notes              = paste(
        "Screened as a covariate on the logit-scale dialyser extraction",
        "ratio in the stepwise covariate search (Kong 2025 Methods 2.4)",
        "and not retained in the final model. Cohort split 14 (74%) HD /",
        "1 (5%) HDF-pre / 4 (21%) HDF-post (Table 1). The Discussion notes",
        "that post-dilution haemodiafiltration is exactly the modality in",
        "which a post-dialyser blood-flow correction would have mattered,",
        "and that the modality was not a significant covariate here.",
        "Documented for provenance only; not referenced in model()."
      )
    ),
    DIALYZER_PRIMEVOL = list(
      description        = "Priming volume of the dialyser blood compartment",
      units              = "mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a covariate on the logit-scale dialyser extraction",
        "ratio (Kong 2025 Methods 2.4) and not retained. Cohort median",
        "115 mL, range 83-115 (Table 1). The source dataset carries it as",
        "KNVOL and the control stream computes V_DIA = KNVOL/1000 in the",
        "$PK block, but V_DIA is never used in the final $DES or $ERROR --",
        "a vestige of an explored dialyser-compartment parameterisation",
        "that the authors report they did not adopt (Discussion: 'we did",
        "not consider implementation of these approaches').",
        "Documented for provenance only; not referenced in model()."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    age_range      = "20-84 years (median 71.5)",
    age_median     = "71.5 years",
    weight_range   = "27.0-93.8 kg (median 66.3)",
    weight_median  = "66.3 kg",
    sex_female_pct = 35,
    disease_state  = paste(
      "End-stage kidney disease on thrice-weekly high-flux intermittent",
      "haemodialysis, with a documented or presumed infection requiring",
      "piperacillin/tazobactam. Residual diuresis present in 11/20 (55%)",
      "and absent in 9/20 (45%). Patients under 18 years, with known",
      "beta-lactam / beta-lactamase-inhibitor hypersensitivity, or",
      "pregnant were excluded (Kong 2025 Methods 2.1)."
    ),
    renal_function = paste(
      "End-stage kidney disease; the model sets the serum-creatinine",
      "correction F_SCR to 1 for these patients (Eq. 2) because serum",
      "creatinine is itself removed by the dialyser, and instead estimates",
      "a residual endogenous clearance fraction theta_ESKD. Endogenous",
      "clearance was estimated to be 78.6% lower (95% CI 66.3-86.4%)",
      "without residual diuresis and 63.0% lower (95% CI 49.5-73.0%) with",
      "residual diuresis, relative to a patient with normal kidney",
      "function (Kong 2025 Results 3.4)."
    ),
    dose_range     = paste(
      "Piperacillin 0.60-4.00 g (median 3.08) with tazobactam",
      "0.075-0.50 g (median 0.385) per dose, administered intravenously",
      "in the clinical fixed 8:1 ratio; infusion duration 0.25-24.50 h",
      "(median 2.85). Regimens were chosen by the treating physician, not",
      "by protocol. The actual delivered dose was corrected for drug",
      "remaining in the used infusion lines (Kong 2025 Methods 2.1,",
      "Table 1)."
    ),
    dialysis       = paste(
      "Thrice-weekly high-flux intermittent haemodialysis; session",
      "duration 89.4-333.4 min (median 240.0), dialyser blood flow rate",
      "160-350 mL/min (median 225), priming volume 83-115 mL",
      "(median 115). Modality 14 (74%) haemodialysis, 1 (5%) pre-dilution",
      "and 4 (21%) post-dilution haemodiafiltration. Vascular access",
      "6 (32%) tunnelled dialysis catheter, 4 (21%) single-needle and",
      "9 (47%) double-needle arteriovenous fistula. Membranes were",
      "polysulfone (n = 14), polyacrylonitrile (3), polyamix (1),",
      "cellulosic triacetate (1) and polyarylethersulfone (1), surface",
      "area 1.3-2 m^2. One patient had no intradialytic sampling, so the",
      "dialysis characteristics are reported for n = 19 (Kong 2025",
      "Results 3.1, Table 1)."
    ),
    regions        = "Belgium (single centre, Ghent University Hospital)",
    notes          = paste(
      "Monocentric prospective observational study, ClinicalTrials.gov",
      "NCT03909698, completed 31 March 2022. 195 blood samples were",
      "collected with both piperacillin and tazobactam quantified in each;",
      "4 piperacillin (2.1%) and 3 tazobactam (1.5%) concentrations were",
      "excluded for bioanalytical anomalies or missing sample-type",
      "information. Intradialytic samples were drawn SIMULTANEOUSLY from",
      "the arterial line entering the dialyser and the venous line",
      "exiting it, which is what makes the extraction ratio identifiable.",
      "Lower limits of quantification 0.5 mg/L (piperacillin) and",
      "0.25 mg/L (tazobactam). Estimation used NONMEM 7.5 FOCE-I",
      "(Kong 2025 Methods 2.1-2.2, 2.6, Results 3.1-3.2).",
      "The Monte Carlo simulations of Table 4 and Figs. 3-4 drew 1000",
      "virtual patients, sampling sex from the observed 35% female / 65%",
      "male split and continuous covariates from normal distributions",
      "truncated to the observed range; sex itself is not a covariate in",
      "the model."
    )
  )

  ini({
    # =====================================================================
    # Structural parameters -- ALL FIXED to the general-purpose parent
    # model. Values are the log-scale $THETA entries transcribed verbatim
    # from the final NONMEM control stream in the Kong 2025 ESM, each of
    # which back-transforms to the corresponding "FIX" row of Table 3.
    # Reference individual: 70 kg, 35 years old (35 + 40/52 postmenstrual
    # years), standardised serum creatinine.
    # =====================================================================
    lvc <- fixed(2.34); label("Piperacillin central volume V1 at 70 kg (L)")                    # ESM $THETA(1) 2.34 FIX -> exp = 10.38 L; Table 3 theta_V1 piperacillin = 10.4 L/70 kg FIX
    lcl <- fixed(2.36); label("Piperacillin endogenous clearance CL at the reference patient (L/h)") # ESM $THETA(2) 2.36 FIX -> exp = 10.59 L/h; Table 3 theta_CL piperacillin = 10.6 L/h/70 kg FIX
    lvp <- fixed(2.45); label("Piperacillin peripheral volume V2 at 70 kg (L)")                 # ESM $THETA(3) 2.45 FIX -> exp = 11.59 L; Table 3 theta_V2 piperacillin = 11.6 L/70 kg FIX
    lq  <- fixed(2.72); label("Piperacillin inter-compartmental clearance Q2 at 70 kg (L/h)")   # ESM $THETA(4) 2.72 FIX -> exp = 15.18 L/h; Table 3 theta_Q2 piperacillin = 15.2 L/h/70 kg FIX

    lvc_taz <- fixed(2.35); label("Tazobactam central volume V3 at 70 kg (L)")                  # ESM $THETA(5) 2.35 FIX -> exp = 10.49 L; Table 3 theta_V1 tazobactam = 10.5 L/70 kg FIX
    lcl_taz <- fixed(2.26); label("Tazobactam endogenous clearance CL at the reference patient (L/h)") # ESM $THETA(6) 2.26 FIX -> exp = 9.583 L/h; Table 3 theta_CL tazobactam = 9.58 L/h/70 kg FIX
    lvp_taz <- fixed(2.62); label("Tazobactam peripheral volume V4 at 70 kg (L)")               # ESM $THETA(7) 2.62 FIX -> exp = 13.74 L; Table 3 theta_V2 tazobactam = 13.7 L/70 kg FIX
    lq_taz  <- fixed(2.82); label("Tazobactam inter-compartmental clearance Q4 at 70 kg (L/h)") # ESM $THETA(8) 2.82 FIX -> exp = 16.78 L/h; Table 3 theta_Q2 tazobactam = 16.8 L/h/70 kg FIX

    # =====================================================================
    # Allometric exponents -- conventional fixed West-Brown-Enquist values,
    # not estimated (Kong 2025 Methods 2.3: "the conventional fixed
    # West-Brown-Enquist exponents of 0.75 for clearance and 1.00 for
    # volume of distribution ... but for the intercompartmental clearance
    # (Q2), compartmental allometry was assumed, where Q2 is scaled by the
    # individual estimated size of the peripheral volume of distribution
    # (V2) with a fixed West-Brown-Enquist exponent of 0.75").
    # =====================================================================
    e_wt_cl <- fixed(0.75); label("Body-weight allometric exponent on endogenous clearance (unitless)") # Methods 2.3: fixed West-Brown-Enquist 0.75 on clearance; ESM $PK 'FSIZE**0.75' on TVCL1 / TVCL3
    e_wt_vc <- fixed(1.0);  label("Body-weight allometric exponent on all volumes of distribution (unitless)") # Methods 2.3: fixed West-Brown-Enquist 1.00 on volume; ESM $PK 'TVV1 = EXP(THETA(1)) * FSIZE' (exponent 1 implicit)
    e_wt_q  <- fixed(0.75); label("Compartmental-allometry exponent of peripheral volume on inter-compartmental clearance (unitless)") # Methods 2.3: compartmental allometry, fixed exponent 0.75; ESM $PK 'TVQ2 = EXP(THETA(4)) * (V2SZ**0.75)'

    # =====================================================================
    # Clearance maturation and decline with postmenstrual age -- FIXED to
    # the general-purpose parent model. Both are sigmoidal in
    # postmenstrual AGE IN YEARS and are normalised by their value at the
    # 35-year reference (35 + 40/52 = 35.769 y), so the ratio is ~1 across
    # the all-adult ESKD cohort. The shape factors are shared between the
    # two drugs; the decline half-age is drug-specific.
    # =====================================================================
    lhill_mat <- fixed(1.21);   label("Log maturation Hill / shape factor for endogenous clearance (log unitless)")         # ESM $THETA(9) 1.21 FIX -> exp = 3.353; Table 3 gamma_1 = 3.35 FIX
    ltmat50   <- fixed(0.0358); label("Log postmenstrual age at 50% of mature endogenous clearance (log years)") # ESM $THETA(10) 0.0358 FIX -> exp = 1.0364 years = 54.1 weeks; Table 3 MAT50 = 54.2 weeks FIX (the ESM code value is used because it is the implemented one; the 0.2% difference is a weeks-per-year rounding and is immaterial in an all-adult cohort)
    lhill_dec <- fixed(0.653);  label("Log decline Hill / shape factor for endogenous clearance (log unitless)")            # ESM $THETA(11) 0.653 FIX -> exp = 1.921; Table 3 gamma_2 = 1.92 FIX
    ltdec50     <- fixed(4.49); label("Log postmenstrual age at 50% decline of piperacillin clearance (log years)") # ESM $THETA(12) 4.49 FIX -> exp = 89.12 years; Table 3 DEC50 piperacillin = 89.1 years FIX
    ltdec50_taz <- fixed(4.12); label("Log postmenstrual age at 50% decline of tazobactam clearance (log years)")   # ESM $THETA(13) 4.12 FIX -> exp = 61.56 years; Table 3 DEC50 tazobactam = 61.6 years FIX

    # =====================================================================
    # Serum-creatinine effect on endogenous clearance -- FIXED to the
    # general-purpose parent model (Kong 2025: "The influence of serum
    # creatinine is based on the model parameter theta_SCR fixed to the
    # value estimated in the development of the general-purpose model").
    # Applies only when the patient is NOT on intermittent haemodialysis
    # (Eq. 2).
    # =====================================================================
    le_creat_cl <- fixed(-1.06); label("Log serum-creatinine effect coefficient on endogenous clearance (log dL/mg)") # ESM $THETA(14) -1.06 FIX -> exp = 0.3465 dL/mg; Table 3 theta_SCR = 0.346 dL/mg FIX

    # =====================================================================
    # ESTIMATED from the ESKD data: residual endogenous clearance fraction
    # (Eq. 4 F_ESKD), shared between piperacillin and tazobactam after the
    # two drug-specific parameters were merged without loss of fit
    # (Table 2 model 10, dOFV +0.92, p = 0.343). Residual diuresis
    # selects between the two values.
    # =====================================================================
    lf_eskd      <- log(0.214); label("Log residual endogenous clearance fraction in ESKD without residual diuresis (fraction of normal-kidney-function CL)") # ESM $THETA(15) estimated; Table 3 theta_ESKD without residual diuresis = 0.214 (95% CI 0.136-0.337) -> log(0.214) = -1.542
    lf_eskd_diur <- log(0.370); label("Log residual endogenous clearance fraction in ESKD with residual diuresis (fraction of normal-kidney-function CL)")    # ESM $THETA(16) estimated; Table 3 theta_ESKD with residual diuresis = 0.370 (95% CI 0.270-0.505) -> log(0.370) = -0.994

    # =====================================================================
    # ESTIMATED from the ESKD data: dialyser extraction ratios, on the
    # logit scale so that ER stays inside (0, 1) (Eq. 6). Piperacillin has
    # a single ER; the tazobactam ER is a three-level function of vascular
    # access type with tunnelled dialysis catheter as the reference.
    # =====================================================================
    logitedia <- log(0.640 / (1 - 0.640)); label("Logit piperacillin dialyser extraction ratio (logit fraction)")                              # ESM $THETA(17) estimated; Table 3 ER piperacillin = 64.0% (95% CI 57.3-70.2) -> logit = 0.575
    logitedia_taz       <- log(0.801 / (1 - 0.801)); label("Logit tazobactam dialyser extraction ratio, tunnelled dialysis catheter (reference)") # ESM $THETA(18) estimated; Table 3 ER of TDC tazobactam = 80.1% (95% CI 75.4-84.1) -> logit = 1.393
    logitedia_taz_avf1n <- log(0.739 / (1 - 0.739)); label("Logit tazobactam dialyser extraction ratio, single-needle arteriovenous fistula (logit fraction)")     # ESM $THETA(19) estimated; Table 3 ER of AVF 1N tazobactam = 73.9% (95% CI 65.3-81.0) -> logit = 1.041
    logitedia_taz_avf2n <- log(0.735 / (1 - 0.735)); label("Logit tazobactam dialyser extraction ratio, double-needle arteriovenous fistula (logit fraction)")     # ESM $THETA(20) estimated; Table 3 ER of AVF 2N tazobactam = 73.5% (95% CI 65.0-80.5) -> logit = 1.020

    # =====================================================================
    # Inter-individual variability -- ALL FIXED to the general-purpose
    # parent model ($OMEGA ... FIX). Table 3 reports these as CV% computed
    # as sqrt(exp(omega) - 1) x 100% (Table 3 footnote b), so the omega
    # values below are variances on the log scale and reproduce the
    # published CV% exactly.
    #
    # Three of the five etas are SHARED across the two drugs, which is the
    # parent model's finding that these parameters are correlated between
    # compounds (Methods 2.3): eta(1) on both central volumes, eta(4) on
    # both peripheral volumes and eta(5) on both inter-compartmental
    # clearances. Only clearance has drug-specific etas. The sharing is
    # implemented by using the SAME eta in both drugs' expressions inside
    # model(), which is why only one variance is declared for each shared
    # eta.
    # =====================================================================
    etalvc     ~ fixed(0.167); label("IIV on central volume, SHARED between piperacillin and tazobactam") # ESM $OMEGA 0.167 FIX '; V1 & V3'; Table 3 IIV of V1 = 42.6% FIX -> sqrt(exp(0.167)-1) = 42.6%
    etalcl     ~ fixed(0.171); label("IIV on piperacillin endogenous clearance")                          # ESM $OMEGA 0.171 FIX '; CL1'; Table 3 IIV of CL piperacillin = 43.2% FIX -> sqrt(exp(0.171)-1) = 43.2%
    etalcl_taz ~ fixed(0.159); label("IIV on tazobactam endogenous clearance")                            # ESM $OMEGA 0.159 FIX '; CL3'; Table 3 IIV of CL tazobactam = 41.5% FIX -> sqrt(exp(0.159)-1) = 41.5%
    etalvp     ~ fixed(0.548); label("IIV on peripheral volume, SHARED between piperacillin and tazobactam") # ESM $OMEGA 0.548 FIX '; V2 & V4'; Table 3 IIV of V2 = 85.4% FIX -> sqrt(exp(0.548)-1) = 85.4%
    etalq      ~ fixed(0.358); label("IIV on inter-compartmental clearance, SHARED between piperacillin and tazobactam") # ESM $OMEGA 0.358 FIX '; Q2 & Q4'; Table 3 IIV of Q2 = 65.6% FIX -> sqrt(exp(0.358)-1) = 65.6%

    # =====================================================================
    # Residual error -- ESTIMATED. Combined proportional plus additive on
    # the linear scale: the ESM builds SD = sqrt(IPRED^2 * PROP^2 + ADD^2)
    # with $SIGMA fixed to 1, which is exactly nlmixr2's default
    # "combined2" add() + prop() form. Removing either additive term
    # worsened the fit (Table 2 models 12-13), so both are retained.
    #
    # The cross-drug residual correlation rho = 95.0% (95% CI 91.4-97.1;
    # ESM $ERROR 'TRE = (EPS(1)*RHO + EPS(2)*SQRT(1-RHO*RHO)) * TSD') is
    # NOT expressible in rxode2 / nlmixr2, which has no cross-endpoint
    # residual correlation for a combined error model. The residuals are
    # therefore declared independently per drug and the correlation is
    # documented in the vignette's Assumptions and deviations section --
    # the same treatment Soto_2014_ampicillin_sulbactam.R applies to its
    # own rho = 0.946.
    # =====================================================================
    propSd     <- 0.412; label("Piperacillin proportional residual SD (fraction)") # ESM $THETA(21) estimated (PPROP = EXP(THETA(21))); Table 3 SD_prop piperacillin = 41.2% (95% CI 0.340-0.499)
    addSd      <- 0.588; label("Piperacillin additive residual SD (mg/L)")         # ESM $THETA(22) estimated (PADD = EXP(THETA(22))); Table 3 SD_add piperacillin = 0.588 mg/L (95% CI 0.346-1.00)
    propSd_taz <- 0.423; label("Tazobactam proportional residual SD (fraction)")   # ESM $THETA(23) estimated (TPROP = EXP(THETA(23))); Table 3 SD_prop tazobactam = 42.3% (95% CI 0.339-0.526)
    addSd_taz  <- 0.214; label("Tazobactam additive residual SD (mg/L)")           # ESM $THETA(24) estimated (TADD = EXP(THETA(24))); Table 3 SD_add tazobactam = 0.214 mg/L (95% CI 0.140-0.327)
  })

  model({
    # ------------------------------------------------------------------
    # 1a. Postmenstrual age in YEARS. The canonical PAGE column carries
    #     months; Kong 2025 works in years throughout (ESM $PK
    #     'RPMA = 35+(40/52)').
    # ------------------------------------------------------------------
    pma  <- PAGE / 12
    rpma <- 35 + 40 / 52

    # ------------------------------------------------------------------
    # 1b. Allometric size factor (ESM $PK 'FSIZE = WT/70').
    # ------------------------------------------------------------------
    fsize <- WT / 70

    # ------------------------------------------------------------------
    # 1c. Sigmoidal maturation of endogenous clearance, normalised to the
    #     35-year reference individual (ESM $PK FMAT / RFMAT).
    # ------------------------------------------------------------------
    hill_mat <- exp(lhill_mat)
    tmat50   <- exp(ltmat50)
    fmat     <- pma^hill_mat  / (pma^hill_mat  + tmat50^hill_mat)
    fmat_ref <- rpma^hill_mat / (rpma^hill_mat + tmat50^hill_mat)

    # ------------------------------------------------------------------
    # 1d. Sigmoidal decline of endogenous clearance with advancing age,
    #     normalised to the same reference (ESM $PK PFD / PRFD for
    #     piperacillin and TFD / TRFD for tazobactam). Note the form is
    #     T50^g / (PMA^g + T50^g), i.e. DECREASING in age, and the shape
    #     factor is shared while the half-age is drug-specific.
    # ------------------------------------------------------------------
    hill_dec     <- exp(lhill_dec)
    tdec50       <- exp(ltdec50)
    tdec50_taz   <- exp(ltdec50_taz)
    fdec         <- tdec50^hill_dec     / (pma^hill_dec  + tdec50^hill_dec)
    fdec_ref     <- tdec50^hill_dec     / (rpma^hill_dec + tdec50^hill_dec)
    fdec_taz     <- tdec50_taz^hill_dec / (pma^hill_dec  + tdec50_taz^hill_dec)
    fdec_taz_ref <- tdec50_taz^hill_dec / (rpma^hill_dec + tdec50_taz^hill_dec)

    # ------------------------------------------------------------------
    # 1e. Postmenstrual-age-standardised serum creatinine (Kong 2025
    #     Eq. 1 / Eq. 3, ESM $PK STCR):
    #       SCR_std = exp[1.42 - (1.17 + 0.203 ln(PMA/100)) / sqrt(PMA/100)]
    #     with PMA in years. The hardcoded 1.42 / 1.17 / 0.203 / 100 are
    #     the published constants of Eq. 1, not fitted parameters.
    # ------------------------------------------------------------------
    creat_std <- exp(1.42 - (1.17 + 0.203 * log(pma / 100)) / sqrt(pma / 100))

    # ------------------------------------------------------------------
    # 1f. Kong 2025 Eq. 2 -- the serum-creatinine correction on endogenous
    #     clearance is the exponential term for patients NOT on
    #     intermittent haemodialysis and exactly 1 for patients on it.
    # ------------------------------------------------------------------
    e_creat_cl <- exp(le_creat_cl)
    f_scr <- (1 - RRT_HEMODIAL_STATUS) * exp(-e_creat_cl * (CREAT - creat_std)) +
      RRT_HEMODIAL_STATUS

    # ------------------------------------------------------------------
    # 1g. Kong 2025 Eq. 4 F_ESKD -- the residual endogenous clearance
    #     fraction, which applies only to ESKD patients on intermittent
    #     haemodialysis and reduces to 1 for a normal-kidney-function
    #     patient. Residual diuresis selects between the two estimated
    #     values (ESM $PK 'CLCKD = EXP(THETA(15)); IF(UVOL0123.NE.0)
    #     CLCKD = EXP(THETA(16))').
    # ------------------------------------------------------------------
    diuresis  <- (URINE_VOL_24H > 100)
    f_eskd_hd <- exp(lf_eskd) * (1 - diuresis) + exp(lf_eskd_diur) * diuresis
    f_eskd    <- (1 - RRT_HEMODIAL_STATUS) + RRT_HEMODIAL_STATUS * f_eskd_hd

    # ------------------------------------------------------------------
    # 2. Individual PK parameters, Kong 2025 Eq. 4:
    #      CL = CL_normal x F_SCR x F_ESKD x F_MAT x F_DEC x F_SIZE x e^eta
    #    Central and peripheral volumes carry the linear weight scalar;
    #    the inter-compartmental clearance uses compartmental allometry on
    #    the INDIVIDUAL peripheral volume (which already carries etalvp).
    #    etalvc, etalvp and etalq are each used TWICE -- once per drug --
    #    implementing the parent model's shared random effects.
    # ------------------------------------------------------------------
    vc     <- exp(lvc     + etalvc) * fsize^e_wt_vc
    vc_taz <- exp(lvc_taz + etalvc) * fsize^e_wt_vc

    cl <- exp(lcl + etalcl) * fsize^e_wt_cl *
      (fmat / fmat_ref) * (fdec / fdec_ref) * f_scr * f_eskd
    cl_taz <- exp(lcl_taz + etalcl_taz) * fsize^e_wt_cl *
      (fmat / fmat_ref) * (fdec_taz / fdec_taz_ref) * f_scr * f_eskd

    vp     <- exp(lvp     + etalvp) * fsize^e_wt_vc
    vp_taz <- exp(lvp_taz + etalvp) * fsize^e_wt_vc

    q     <- exp(lq     + etalq) * (vp     / exp(lvp))^e_wt_q
    q_taz <- exp(lq_taz + etalq) * (vp_taz / exp(lvp_taz))^e_wt_q

    # ------------------------------------------------------------------
    # 3a. Dialyser clearance, Kong 2025 Eq. 5 and 6:
    #       ER      = e^theta_lgtER / (1 + e^theta_lgtER)
    #       CL_DIA  = ER x blood flow rate (L/h)
    #     gated on by the time-varying during-session indicator. The
    #     tazobactam logit is selected by vascular access type, with
    #     tunnelled dialysis catheter as the reference (both indicators 0).
    # ------------------------------------------------------------------
    logitedia_taz_i <- logitedia_taz +
      (logitedia_taz_avf1n - logitedia_taz) * VASCACC_AVF1N +
      (logitedia_taz_avf2n - logitedia_taz) * VASCACC_AVF2N

    edia     <- exp(logitedia)       / (1 + exp(logitedia))
    edia_taz <- exp(logitedia_taz_i) / (1 + exp(logitedia_taz_i))

    bfr_lph <- BFR * 60 / 1000  # ESM $PK 'FLOW = IFLOW*60/1000' -- mL/min to L/h

    cl_dialysis     <- edia     * bfr_lph * RRT_HEMODIAL_ACTIVE
    cl_dialysis_taz <- edia_taz * bfr_lph * RRT_HEMODIAL_ACTIVE

    # ------------------------------------------------------------------
    # 3b. Micro-constants (ESM $PK k10 / k12 / k21 / kdia_PIP and the
    #     tazobactam counterparts k30 / k34 / k43 / kdia_TAZ).
    # ------------------------------------------------------------------
    kel  <- cl / vc
    k12  <- q  / vc
    k21  <- q  / vp
    kdia <- cl_dialysis / vc

    kel_taz  <- cl_taz / vc_taz
    k12_taz  <- q_taz  / vc_taz
    k21_taz  <- q_taz  / vp_taz
    kdia_taz <- cl_dialysis_taz / vc_taz

    # ------------------------------------------------------------------
    # 4. Two-compartment IV disposition per drug (ESM $DES). The two drugs
    #    share random effects but do not interconvert, and each is dosed
    #    into its own central compartment in the fixed clinical 8:1 ratio.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12 + kdia) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central                - k21 * peripheral1

    d/dt(central_taz)     <- -(kel_taz + k12_taz + kdia_taz) * central_taz + k21_taz * peripheral1_taz
    d/dt(peripheral1_taz) <-   k12_taz * central_taz                        - k21_taz * peripheral1_taz

    # ------------------------------------------------------------------
    # 5. Observations. Dose in mg and volumes in L give mg/L = ug/mL.
    #
    #    Cc / Cc_taz are the PRE-dialyser (dialyser-inlet, arterial-line)
    #    plasma concentrations and are the model's declared endpoints,
    #    matching ESM $ERROR CPRE_PIP / CPRE_TAZ.
    #
    #    Cc_post / Cc_post_taz are the POST-dialyser (dialyser-outlet,
    #    venous-line) concentrations of ESM $ERROR CPOST_PIP / CPOST_TAZ:
    #    blood leaving the dialyser has lost the extracted fraction, so the
    #    outlet concentration is C x (1 - ER) while a session is running
    #    and equals the inlet concentration otherwise. The source dataset
    #    used a PRE1POST2 column to route each observation to the matching
    #    prediction; here both are returned as derived variables so a user
    #    can compare against either sampling site.
    # ------------------------------------------------------------------
    Cc     <- central     / vc
    Cc_taz <- central_taz / vc_taz

    Cc_post     <- Cc     * (1 - edia     * RRT_HEMODIAL_ACTIVE)
    Cc_post_taz <- Cc_taz * (1 - edia_taz * RRT_HEMODIAL_ACTIVE)

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_taz ~ add(addSd_taz) + prop(propSd_taz)
  })
}
