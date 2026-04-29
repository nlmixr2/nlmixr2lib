Ide_2020_elotuzumab <- function() {
  description <- "Two-compartment population PK model for elotuzumab (anti-SLAMF7 humanized IgG1) in Japanese and non-Japanese patients with multiple myeloma (Ide 2020); parallel linear and Michaelis-Menten elimination from the central compartment plus second-order target-mediated elimination from the peripheral compartment driven by a non-renewable target pool, with time-varying serum M protein on Vmax."
  reference   <- "Ide T, Roy A, Imai Y, Vezina HE. Model-Based Determination of Elotuzumab Pharmacokinetics in Japanese Patients With Multiple Myeloma Incorporating Time-Varying M Protein. J Clin Pharmacol. 2021;61(1):64-73. doi:10.1002/jcph.1698"
  vignette    <- "Ide_2020_elotuzumab"
  units       <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL, VC, Q (fixed exponent 0.75), and VP. Reference 75 kg (Ide 2020 supplement S2 'Supplemental Information for Table 2' covariate-equation reference patient).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 65 years (Ide 2020 supplement S2 reference patient).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) — matches Ide 2020 reference patient (supplement S2: 'Male').",
      notes              = "Exponential effects on both CL (CLSEX = 1.06) and VC (VCSEX = 0.808). The source NONMEM column SEXF already follows the canonical 1 = female / 0 = male convention.",
      source_name        = "SEXF"
    ),
    RACE_ASIAN = list(
      description        = "Indicator for Asian race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian) — matches Ide 2020 reference patient (supplement S2: 'non-Asian').",
      notes              = "Exponential effects on both CL (CLRACE = 0.897) and VC (VCRACE = 0.853). Decomposed from the source NONMEM column RACEN via ASIAN = as.integer(RACEN == 3).",
      source_name        = "RACEN"
    ),
    CRCL = list(
      description        = "Baseline estimated glomerular filtration rate (eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 100 mL/min/1.73 m^2 (Ide 2020 Methods cite 'estimated glomerular filtration rate (eGFR)'; Table 1 reports eGFR in mL/min/1.73 m^2; supplement S2 lists the reference as 'GFR = 100 mL/min' which omits the /1.73 m^2 qualifier — see vignette Errata). Source NONMEM column GFR; stored under the canonical CRCL.",
      source_name        = "GFR"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 200 U/L (Ide 2020 supplement S2 reference patient). The source NONMEM column LDH_Z carries imputation sentinel -99 for missing values; the control stream replaces missing with the population median 194 U/L (close to but not equal to the 200 U/L reference; both come from Ide 2020 supplement 7 NONMEM control stream).",
      source_name        = "LDH_Z"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL with reference 3.5 g/dL (Ide 2020 supplement S2 reference patient).",
      source_name        = "ALB"
    ),
    B2M = list(
      description        = "Baseline serum beta-2-microglobulin",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in Ide 2020 only as two binary thresholded indicators: B2M >= 2.0 mg/L (CLB2MICG>2 = 1.11, VCB2MICG>2 = 1.05) and B2M >= 3.5 mg/L (CLB2MICG>3.5 = 1.01, VCB2MICG>3.5 = 1.07). The thresholds and effects act simultaneously: a subject with B2M >= 3.5 mg/L receives both factors. The source NONMEM column B2MICG_Z is reported in mg/dL (Table 1: 0.04-3.47 with median 0.330; thresholds 0.20 and 0.35 mg/dL); supplement S2 abbreviation list states the thresholds in mg/L (2.0 and 3.5 mg/L). 1 mg/dL = 10 mg/L, so the two unit conventions describe the same cutpoints. To match the canonical B2M unit (mg/L) and the supplement's abbreviation list, supply B2M in mg/L. Missing-value sentinel -99 in the source dataset is imputed to the population median 0.33 mg/dL = 3.3 mg/L by the control stream (which is below the 3.5 mg/L threshold).",
      source_name        = "B2MICG_Z"
    ),
    HEPIMP = list(
      description        = "Baseline hepatic-impairment indicator (NCI ODWG mild or worse vs. normal)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (normal hepatic function, NCI ODWG group 1) — matches Ide 2020 reference patient (supplement S2: 'Normal Liver Function').",
      notes              = "Exponential effect on CL (CLHEP = 0.91). Decomposed from the source NONMEM column HEP via HEPIMP = as.integer(HEP > 0.5). Per Ide 2020 Methods, hepatic impairment was specified using NCI Organ Dysfunction Working Group criteria.",
      source_name        = "HEP"
    ),
    ECOG_GE1 = list(
      description        = "Baseline ECOG performance-status indicator (1 if ECOG >= 1, else 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG = 0) — matches Ide 2020 reference patient (supplement S2: 'ECOG = 0').",
      notes              = "Paired with ECOG_GE2 to retain the three-level ordinal effect Ide 2020 reports for ECOG = 0 / 1 / >=2. Exponential effect on CL (CLECOG>0 = 1.03). Decomposed from the source NONMEM column ECOG101 via ECOG_GE1 = as.integer(ECOG101 >= 1).",
      source_name        = "ECOG101"
    ),
    ECOG_GE2 = list(
      description        = "Baseline ECOG performance-status indicator (1 if ECOG >= 2, else 0)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ECOG <= 1) — supplied alongside ECOG_GE1 to encode the three-level ECOG = 0 / 1 / >=2 ordinal effect; Ide 2020 reference patient has both indicators = 0 (ECOG = 0).",
      notes              = "Exponential effect on CL (CLECOG>1 = 1.15). Decomposed from the source NONMEM column ECOG101 via ECOG_GE2 = as.integer(ECOG101 >= 2). Always paired with ECOG_GE1 in the Ide 2020 model.",
      source_name        = "ECOG101"
    ),
    LINE_1L = list(
      description        = "First-line-therapy indicator (1 = treatment-naive, 0 = >=1 prior line of therapy)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (>=1 prior line of therapy) — matches Ide 2020 reference patient (supplement S2: 'with >= 1 prior line of therapy'). The Ide 2020 cohort enriched the prior analysis with newly-diagnosed previously-untreated Japanese patients (NCT02272803), so the LINE_1L = 1 group is a meaningful subset of the data.",
      notes              = "Exponential effects on both CL (CLLINE=0 = 0.921) and Vmax of the Michaelis-Menten elimination (VMAXLINE=0 = 1.01). Decomposed from the source NONMEM column LINEF via LINE_1L = as.integer(LINEF == 0); LINEF = 0 in the source dataset means 'no prior line of therapy', i.e., first-line / treatment-naive.",
      source_name        = "LINEF"
    ),
    MCPROT = list(
      description        = "Time-varying serum M (monoclonal) protein concentration",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying multiple-myeloma tumor-burden marker. Enters the model un-log-transformed via exp(0.277 * MCPROT) on Vmax of the Michaelis-Menten target-mediated elimination from the central compartment. The supplement S2 reference for VMAX,REF = 12.2 ug/mL/day is MCPROT = 0 g/dL (the value at which the exp(*MCPROT) term is unity); the figure-1 covariate-effect plot uses MCPROT = 2.0 g/dL as the reference patient. Time-varying values must be supplied at every PK observation time; Ide 2020 Methods specify linear interpolation between observed M-protein measurements with last-observation-carried-forward beyond the last sample. The source NONMEM column TMCPROT_Z (with imputation sentinel -99 replaced by population median 2.1 g/dL) is the time-varying counterpart of the baseline-only MCPROT_Z column carried for descriptive statistics.",
      source_name        = "TMCPROT"
    ),
    COMBO_LEN_DEX = list(
      description        = "Lenalidomide + dexamethasone combination-therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no Ld coadministration) — note: this differs from the Ide 2020 reference patient, which has Ld coadministration (supplement S2: 'with lenalidomide/dexamethasone co-administration'). The model encodes the effects so that COMBO_LEN_DEX = 1 yields factor 1 (paper's reference) and COMBO_LEN_DEX = 0 activates the (LENDEX-1)-shifted exponential coefficients exp(theta * (COMBO_LEN_DEX - 1)) for both CL and KINT.",
      notes              = "Multiplicative effects: CLLd = 0.74 on nonspecific CL (CL is 35% higher in the Ld-free arm relative to the Ld arm) and KINTLd = 10.1 on the second-order target-mediated elimination rate from the peripheral compartment (KINT is ~10x lower in the Ld-free arm; lenalidomide-activated NK cells are hypothesized to drive the elevated target-binding-mediated elimination in the Ld arm). Decomposed from the source NONMEM column LENDEX, derived in the control stream as `LENDEX = as.integer(STUDY != 204011)` because study 204011 was the Ld-free elotuzumab-monotherapy cohort.",
      source_name        = "LENDEX"
    )
  )

  population <- list(
    n_subjects     = 420L,
    n_observations = 8125L,
    n_studies      = 5L,
    age_range      = "37-88 years (median 67; mean 66.1, SD 9.66)",
    age_median     = "67 years",
    weight_range   = "33.8-150 kg (median 74.0; mean 74.3, SD 17.6)",
    weight_median  = "74.0 kg",
    sex_female_pct = 42.1,
    race_ethnicity = c(`non-Asian` = 81.0, Asian = 19.0),
    ethnicity_japanese_pct = 18.3,
    ecog_distribution = "ECOG 0 51.7%, ECOG 1 41.0%, ECOG 2 7.4%",
    renal_function = "Baseline eGFR mean 73.5 (SD 23.5) mL/min/1.73 m^2; range 4.58-124",
    hepatic_function = "Normal 91.4%, mild-or-worse impairment 8.6%",
    coadministration = "Lenalidomide+dexamethasone (Ld) coadministration: 92.6% (n=389); none: 7.4% (n=31). Ld-free arm came from study 204011 (NCT01441973 Cohort 1), the elotuzumab-monotherapy cohort.",
    line_of_therapy = "0 prior lines (treatment-naive): 16.4% (n=69); 1 prior line: 39.3% (n=165); 2 or more prior lines: 44.3% (n=186)",
    mprotein_baseline = "Baseline serum M-protein mean 2.39 g/dL (SD 1.59); median 2.10 g/dL (range 0.00-7.70)",
    b2m_baseline   = "Baseline serum beta-2-microglobulin mean 0.440 mg/dL (SD 0.375); median 0.330 mg/dL (range 0.040-3.47); 2.4% missing",
    albumin_baseline = "Baseline serum albumin mean 3.78 g/dL (SD 0.581); median 3.80 g/dL (range 1.90-5.00)",
    ldh_baseline   = "Baseline serum LDH mean 239 U/L (SD 146); median 194 U/L (range 54.0-1900); 6.4% missing",
    disease_state  = "Multiple myeloma (mix of newly diagnosed previously untreated and relapsed/refractory). 16.4% treatment-naive (NCT02272803, Japanese newly-diagnosed; pivotal addition over the prior PPK analysis); 83.6% with >=1 prior line of therapy.",
    dose_range     = "10 or 20 mg/kg IV infusion. Standard regimen: 10 mg/kg weekly for cycles 1-2, then 10 mg/kg every 2 weeks from cycle 3 (continuing thereafter, or switching to 20 mg/kg every 4 weeks beginning cycle 19 in NCT02272803). NCT01441973 Cohort 1 used 20 mg/kg every 4 weeks from cycle 2 as monotherapy.",
    regions        = "Five clinical studies pooled: 2 phase 1 (NCT01241292, NCT01393964), 2 phase 2 (NCT01441973, NCT02272803), and 1 phase 3 (NCT01239797 / ELOQUENT-2). NCT02272803 is Japan-only; the others are global.",
    notes          = "Baseline demographics from Ide 2020 Table 1 (n = 420 patients, 8125 elotuzumab serum concentrations). Of 420 patients, 77 (18.3%) were Japanese. Bioanalytical assay: quantitative ELISA, LLOQ 190 ng/mL. Studies pooled via supplement S1 (Table S1 in PMID_32656777_supplement_5_trimmed.md): NCT02272803 (Japanese newly-diagnosed; 10 mg/kg Q2W cycles 3-18, 20 mg/kg Q4W from cycle 19), NCT01239797 (ELOQUENT-2 phase 3, 10 mg/kg Q2W from cycle 3), NCT01241292 (phase 1, 10 or 20 mg/kg Q2W from cycle 3), NCT01393964 (phase 1, 10 mg/kg Q2W from cycle 4), NCT01441973 (phase 2, monotherapy: Cohort 1 20 mg/kg Q4W or Cohort 2 10 mg/kg Q2W)."
  )

  ini({
    # Structural parameters - reference values for the supplement S2 reference
    # patient (WT = 75 kg, GFR = 100 mL/min/1.73 m^2, LDH = 200 U/L, ALB = 3.5
    # g/dL, AGE = 65 years, Male, non-Asian, normal liver function, with
    # lenalidomide/dexamethasone co-administration, ECOG = 0, B2M < 2.0 mg/L,
    # >= 1 prior line of therapy, MCPROT = 0 g/dL for VMAX scaling).
    lcl    <- log(0.0806);    label("Nonspecific clearance CL_REF (L/day)")                                            # Ide 2020 Table 2: CLREF = 0.0806 L/day
    lvc    <- log(3.94);      label("Central volume of distribution VC_REF (L)")                                       # Ide 2020 Table 2: VCREF = 3.94 L
    lq     <- log(0.515);     label("Intercompartmental clearance Q_REF (L/day)")                                      # Ide 2020 Table 2: QREF = 0.515 L/day
    lvp    <- log(2.01);      label("Peripheral volume of distribution VP_REF (L)")                                    # Ide 2020 Table 2: VPREF = 2.01 L

    # Target-mediated elimination - parallel Michaelis-Menten on central
    # (Vmax/Km) and second-order on peripheral (KINT * peripheral_amount *
    # target_concentration). Initial target concentration = RMAX (per Ide 2020
    # supplement S2: "RMAX = initial target concentration in the peripheral
    # compartment" and the NONMEM control-stream initial-condition `A_0(3) =
    # RMAX`). Target pool is non-renewable: it depletes monotonically over
    # treatment time as drug binds.
    lvmax  <- log(12.2);      label("Maximum Michaelis-Menten elimination rate VMAX_REF (ug/mL/day) at MCPROT = 0")    # Ide 2020 Table 2: VMAX,REF = 12.2 ug/mL/day; supplement S2 footnote: VMAX,REF reference includes MCPROT = 0 g/dL with >= 1 prior line of therapy
    lkm    <- log(298);       label("Michaelis-Menten constant KM (ug/mL)")                                            # Ide 2020 Table 2: KM = 298 ug/mL
    lrmax  <- log(832);       label("Initial target concentration in peripheral compartment RMAX (ug/mL)")             # Ide 2020 Table 2: RMAX = 832 ug/mL
    lkint  <- log(0.207e-3);  label("Second-order target-mediated elimination rate KINT_REF (1/day/(ug/mL))")          # Ide 2020 Table 2: KINT = 0.207 x 10^-3 /day/(ug/mL); supplement S2 footnote: reference includes Ld co-administration; the source NONMEM control stream rescales as KINT = exp(MU)/1000

    # Covariate effects on CL (Ide 2020 Table 2). Power exponents on continuous
    # covariates (multiplicative as (cov/ref)^theta). Exponential coefficients
    # on categorical / binary covariates (multiplicative as exp(theta * X) =
    # factor^X where Table 2 reports factor = exp(theta) directly).
    e_wt_cl         <-  1.33;          label("Power exponent of WT on CL (unitless)")                                  # Ide 2020 Table 2: CLWT = 1.33
    e_combo_len_dex_cl   <-  log(0.74);     label("Exponential coefficient of (COMBO_LEN_DEX - 1) on CL (unitless)")             # Ide 2020 Table 2: CLLd = 0.74; encoded as exp(log(0.74) * (COMBO_LEN_DEX - 1)) so COMBO_LEN_DEX = 1 (Ld+) gives factor 1 and COMBO_LEN_DEX = 0 (Ld-) gives factor 1.35 (paper Discussion: 35% higher CL in Ld-free arm)
    e_sexf_cl       <-  log(1.06);     label("Exponential coefficient of SEXF on CL (unitless)")                       # Ide 2020 Table 2: CLSEX = 1.06
    e_race_asian_cl <-  log(0.897);    label("Exponential coefficient of RACE_ASIAN on CL (unitless)")                 # Ide 2020 Table 2: CLRACE = 0.897
    e_age_cl        <-  0.179;         label("Power exponent of AGE on CL (unitless)")                                 # Ide 2020 Table 2: CLAGE = 0.179
    e_crcl_cl       <-  0.121;         label("Power exponent of CRCL (eGFR) on CL (unitless)")                         # Ide 2020 Table 2: CLeGFR = 0.121
    e_ldh_cl        <-  0.0816;        label("Power exponent of LDH on CL (unitless)")                                 # Ide 2020 Table 2: CLLDH = 0.0816
    e_hepimp_cl     <-  log(0.91);     label("Exponential coefficient of HEPIMP on CL (unitless)")                     # Ide 2020 Table 2: CLHEP = 0.91
    e_alb_cl        <- -0.346;         label("Power exponent of ALB on CL (unitless)")                                 # Ide 2020 Table 2: CLALB = -0.346
    e_ecog_ge1_cl   <-  log(1.03);     label("Exponential coefficient of ECOG_GE1 on CL (unitless)")                   # Ide 2020 Table 2: CLECOG>0 = 1.03
    e_ecog_ge2_cl   <-  log(1.15);     label("Exponential coefficient of ECOG_GE2 on CL (unitless)")                   # Ide 2020 Table 2: CLECOG>1 = 1.15
    e_b2m_ge2_cl    <-  log(1.11);     label("Exponential coefficient of (B2M >= 2 mg/L) on CL (unitless)")            # Ide 2020 Table 2: CLB2MICG>=0.20 = 1.11 (mg/dL threshold; equivalent to >= 2 mg/L)
    e_b2m_ge35_cl   <-  log(1.01);     label("Exponential coefficient of (B2M >= 3.5 mg/L) on CL (unitless)")          # Ide 2020 Table 2: CLB2MICG>=0.35 = 1.01 (mg/dL threshold; equivalent to >= 3.5 mg/L)
    e_line_1l_cl    <-  log(0.921);    label("Exponential coefficient of LINE_1L on CL (unitless)")                    # Ide 2020 Table 2: CLLINE=0 = 0.921

    # Covariate effects on VC (Ide 2020 Table 2).
    e_wt_vc         <-  0.348;         label("Power exponent of WT on VC (unitless)")                                  # Ide 2020 Table 2: VCWT = 0.348
    e_sexf_vc       <-  log(0.808);    label("Exponential coefficient of SEXF on VC (unitless)")                       # Ide 2020 Table 2: VCSEX = 0.808
    e_race_asian_vc <-  log(0.853);    label("Exponential coefficient of RACE_ASIAN on VC (unitless)")                 # Ide 2020 Table 2: VCRACE = 0.853
    e_b2m_ge2_vc    <-  log(1.05);     label("Exponential coefficient of (B2M >= 2 mg/L) on VC (unitless)")            # Ide 2020 Table 2: VCB2MICG>=0.20 = 1.05
    e_b2m_ge35_vc   <-  log(1.07);     label("Exponential coefficient of (B2M >= 3.5 mg/L) on VC (unitless)")          # Ide 2020 Table 2: VCB2MICG>=0.35 = 1.07

    # Covariate effects on Q and VP - body-weight allometry only.
    # Q: exponent fixed at allometric value 0.75 because the estimated value
    # was imprecise (Ide 2020 Table 2 footnote c).
    e_wt_q          <- fix(0.75);      label("Power exponent of WT on Q (unitless; fixed allometric value)")           # Ide 2020 Table 2: QWT = 0.75 FIXED
    e_wt_vp         <-  0.623;         label("Power exponent of WT on VP (unitless)")                                  # Ide 2020 Table 2: VPWT = 0.623

    # Covariate effects on VMAX and KINT (Ide 2020 Table 2). MCPROT enters
    # un-log-transformed (i.e., exp(theta * MCPROT) directly, NOT (MCPROT/ref)^theta).
    e_mcprot_vmax   <-  0.277;         label("Exponential coefficient of MCPROT on VMAX (1/(g/dL))")                   # Ide 2020 Table 2: VMAXMCPROT = 0.277
    e_line_1l_vmax  <-  log(1.01);     label("Exponential coefficient of LINE_1L on VMAX (unitless)")                  # Ide 2020 Table 2: VMAXLINE=0 = 1.01
    e_combo_len_dex_kint <-  log(10.1);     label("Exponential coefficient of (COMBO_LEN_DEX - 1) on KINT (unitless)")           # Ide 2020 Table 2: KINTLd = 10.1; encoded as exp(log(10.1) * (COMBO_LEN_DEX - 1)) so COMBO_LEN_DEX = 1 (Ld+) gives factor 1 and COMBO_LEN_DEX = 0 (Ld-) gives factor 1/10.1 (paper Discussion: ~10-fold lower KINT in Ld-free arm)

    # Inter-individual variability (Ide 2020 Table 2 omega^2 values).
    # All etas are independent log-normal in the source NONMEM control stream
    # (no $OMEGA BLOCK; etas appear in successive `eta(1)..eta(8)` slots on
    # individual MU equations). VMAX has effectively zero IIV (omega^2 = 1e-4
    # FIXED, retained per the IMPMAP method requirement; supplement S2
    # footnote d).
    etalcl   ~ 0.156      # Ide 2020 Table 2: omega^2_CL    = 0.156
    etalvc   ~ 0.0355     # Ide 2020 Table 2: omega^2_VC    = 0.0355
    etalq    ~ 0.427      # Ide 2020 Table 2: omega^2_Q     = 0.427
    etalvp   ~ 0.137      # Ide 2020 Table 2: omega^2_VP    = 0.137
    etalrmax ~ 0.193      # Ide 2020 Table 2: omega^2_Rmax  = 0.193
    etalkint ~ 1.84       # Ide 2020 Table 2: omega^2_KINT  = 1.84
    etalvmax ~ fix(0.0001)  # Ide 2020 Table 2: omega^2_VMAX = 0.0001 FIXED (set near zero per IMPMAP requirement)
    etalkm   ~ 0.392      # Ide 2020 Table 2: omega^2_KM    = 0.392

    # Residual error - saturable concentration-dependent log-scale additive
    # error model. Ide 2020 Table 2 reports four parameters: SDL = 2.78 (SD at
    # low concentrations), SDH = 0.0984 (SD at high concentrations), SD50 =
    # 5.56 (concentration at which SD = (SDL+SDH)/2), and SDphase1,2 = 0.707
    # (study-specific multiplier on residual SD for non-204004 studies). The
    # library model uses the saturable W(Cc) form W = SDL - (SDL-SDH)*Cc /
    # (SD50+Cc) directly; the study multiplier (THETA(16)^STOTHER) and the
    # IIV on residual magnitude (omega^2_EPS = 0.164) are intentionally
    # omitted (deviations documented in the vignette Assumptions section).
    sdL  <- 2.78;    label("Log-scale residual SD at low concentrations (unitless)")                                   # Ide 2020 Table 2: SDL  = 2.78
    sdH  <- 0.0984;  label("Log-scale residual SD at high concentrations (unitless)")                                  # Ide 2020 Table 2: SDH  = 0.0984
    sd50 <- 5.56;    label("Concentration at which residual SD is (SDL+SDH)/2 (ug/mL)")                                # Ide 2020 Table 2: SD50 = 5.56
  })

  model({
    # ---- Derived binary indicators from the canonical continuous B2M -------
    # Ide 2020 supplement 7 NONMEM control stream (PMID_32656777_supplement_7_trimmed.md):
    #   B2MICG1 = 0; IF(B2MICG.GE.0.2) B2MICG1 = 1   ; threshold 0.2 mg/dL = 2.0 mg/L
    #   B2MICG2 = 0; IF(B2MICG.GE.0.35) B2MICG2 = 1  ; threshold 0.35 mg/dL = 3.5 mg/L
    # The two indicators act simultaneously (a subject with B2M >= 3.5 mg/L
    # has both B2M_GE2 = 1 and B2M_GE35 = 1).
    B2M_GE2  <- (B2M >= 2.0)
    B2M_GE35 <- (B2M >= 3.5)

    # ---- Individual structural parameters with covariate adjustments ------
    # Reference patient covariate values are listed in supplement S2; effects
    # match the NONMEM control-stream MU equations (PMID_32656777_supplement_7).
    cl <- exp(lcl + etalcl) *
      (WT   / 75)^e_wt_cl *
      (AGE  / 65)^e_age_cl *
      (CRCL / 100)^e_crcl_cl *
      (LDH  / 200)^e_ldh_cl *
      (ALB  / 3.5)^e_alb_cl *
      exp(e_combo_len_dex_cl   * (COMBO_LEN_DEX - 1)) *
      exp(e_sexf_cl       * SEXF) *
      exp(e_race_asian_cl * RACE_ASIAN) *
      exp(e_hepimp_cl     * HEPIMP) *
      exp(e_ecog_ge1_cl   * ECOG_GE1) *
      exp(e_ecog_ge2_cl   * ECOG_GE2) *
      exp(e_b2m_ge2_cl    * B2M_GE2) *
      exp(e_b2m_ge35_cl   * B2M_GE35) *
      exp(e_line_1l_cl    * LINE_1L)

    vc <- exp(lvc + etalvc) *
      (WT / 75)^e_wt_vc *
      exp(e_sexf_vc       * SEXF) *
      exp(e_race_asian_vc * RACE_ASIAN) *
      exp(e_b2m_ge2_vc    * B2M_GE2) *
      exp(e_b2m_ge35_vc   * B2M_GE35)

    q  <- exp(lq  + etalq)  * (WT / 75)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp

    # Vmax: MCPROT enters un-log-transformed; LINE_1L enters as an exponential
    # coefficient on the un-transformed indicator.
    vmax <- exp(lvmax + etalvmax) *
      exp(e_mcprot_vmax  * MCPROT) *
      exp(e_line_1l_vmax * LINE_1L)

    km   <- exp(lkm   + etalkm)
    rmax <- exp(lrmax + etalrmax)

    # KINT: reference is COMBO_LEN_DEX = 1 (paper's "with Ld coadministration"
    # reference). Encoded as exp(log(10.1) * (COMBO_LEN_DEX - 1)) so COMBO_LEN_DEX = 1
    # gives factor 1 and COMBO_LEN_DEX = 0 gives factor 1/10.1.
    kint <- exp(lkint + etalkint) * exp(e_combo_len_dex_kint * (COMBO_LEN_DEX - 1))

    # ---- Micro-constants for the linear two-compartment skeleton ----------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ---- Concentrations ---------------------------------------------------
    # central is in mg (amount); peripheral1 is in mg (amount); target is the
    # peripheral target *concentration* in ug/mL (Ide 2020 supplement S2:
    # "Division by VP is required as A2 is amount while A3 is concentration").
    Cc <- central / vc
    Cp <- peripheral1 / vp

    # ---- ODE system (Ide 2020 supplement 6 + 7 NONMEM $DES block) ---------
    # Central (A1, mg): linear elimination -kel * A1, distribution to/from
    # peripheral, and Michaelis-Menten elimination -VMAX * A1 / (Cc + KM)
    # (the M-M term is a rate of mass elimination scaled by the central
    # amount; equivalent to -VMAX * vc * Cc / (Cc + KM)).
    # Peripheral drug amount (A2, mg): distribution and second-order target
    # elimination -KINT * A2 * target_concentration.
    # Target concentration in peripheral (A3, ug/mL): non-renewable
    # consumption only, -KINT * (A2/VP) * A3 = -KINT * Cp * target.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - vmax * central / (Cc + km)
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1 - kint * peripheral1 * target
    d/dt(target)      <- -kint * Cp * target

    # Initial peripheral target concentration = RMAX (non-renewable target
    # pool; depletes monotonically over treatment time).
    target(0) <- rmax

    # ---- Observation and saturable log-scale residual error ---------------
    # Ide 2020 supplement 7 $ERROR block:
    #   TY  = A1/VC
    #   LTY = LOG(TY)  (with floor LTY = -2.5 if TY <= 0)
    #   W   = SDL - (SDL - SDH) * TY / (SD50 + TY)         ; saturable SD
    #   Y   = LTY + W * EPS(1)                              ; log-additive
    # Y = LTY + W * eps with eps ~ N(0,1) is equivalent to lnorm(W) on
    # linear-space Cc, with W concentration-dependent. Floor is omitted in
    # the library model because rxode2 simulations of a positive-state PK
    # model do not produce non-positive concentrations.
    W <- sdL - (sdL - sdH) * Cc / (sd50 + Cc)
    Cc ~ lnorm(W)
  })
}
