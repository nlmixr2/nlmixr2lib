Ding_2026_vancomycin <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for vancomycin in adult cardiac-surgery patients,",
    "parameterised on clearance and distribution volume. Fat-free mass is computed inside the model by",
    "the Janmahasatian height/weight/sex equations and drives theory-based allometric scaling of both",
    "parameters (exponent 0.75 on CL, 1 on V) against a 56.1 kg standard. Clearance carries three",
    "power-model covariates - serum creatinine, cystatin C and NT-proBNP - so that both renal function",
    "and cardiac stress reduce vancomycin elimination; distribution volume carries neutrophil count,",
    "cystatin C and age. Exponential inter-individual and inter-occasion variability on CL and V, and a",
    "combined additive-plus-proportional residual error whose magnitudes switch between the two",
    "bioanalytical assays used (CMIA and EMIT).",
    sep = " "
  )
  reference <- paste(
    "Ding Y, Xue L, Qin Q, Huang H, Chen Y, Shen H, Hu Y, Chen Y, Miao L, Shen Z.",
    "Exploring the effects of renal and cardiac functions on the pharmacokinetics of vancomycin in",
    "patients undergoing cardiac surgery: a population pharmacokinetic analysis.",
    "Int J Clin Pharm. 2026;48:808-818. doi:10.1007/s11096-025-02077-w.",
    "Structural equations, allometric scaling and covariate forms are taken from the Online Resource",
    "(Equations S1-S9) and from the NONMEM control stream deposited at the end of that Online Resource",
    "(\"Code for individual vancomycin dosing on the basis of target concentration approach (TCA)",
    "theory\"), which reproduces the final model's $PK and $ERROR blocks verbatim.",
    sep = " "
  )
  vignette <- "Ding_2026_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 65 kg (2.5-97.5 percentiles 41-101). Does not enter CL or V",
        "directly. It is an input to the Janmahasatian fat-free-mass equation (Online Resource",
        "Equation S6), which is what the allometric terms scale on. The control stream reads it as",
        "the `Weight` column and assigns `WT = Weight ;kg`.",
        sep = " "
      ),
      source_name        = "Weight"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 169 cm (2.5-97.5 percentiles 150-181). The control stream divides",
        "by 100 to obtain metres before squaring (`HTM = HTCM/100`), which fixes the unit as cm.",
        "Enters only through fat-free mass.",
        sep = " "
      ),
      source_name        = "Hight"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Ding 2026 Table 1: 257 male (80.3%) / 63 female (19.7%). The control stream uses the",
        "OPPOSITE polarity, `M1F0 = Sex` with 1 = male: `IF (M1F0.EQ.0) THEN` selects the female",
        "Janmahasatian constants (WHSmax 37.99, WHS50 35.98) and the ELSE branch selects the male",
        "constants (WHSmax 42.92, WHS50 30.93). Converted to the canonical SEXF (1 = female) via",
        "SEXF = 1 - M1F0, so SEXF selects the female branch. Sex was tested as a categorical covariate",
        "on CL and V and rejected (dOFV 0.138, p = 0.71 and 1.096, p = 0.295; Online Resource Results),",
        "so it acts only through fat-free mass.",
        sep = " "
      ),
      source_name        = "Sex"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 65 umol/L (2.5-97.5 percentiles 36-226.6). Power model on CL",
        "centred at 80 umol/L: `FScr_CL = (Scr/80)**FScr_CL1` in the control stream, exponent -0.458",
        "(Table 3). The centring constant 80 umol/L is a rounded normal-range value, not the cohort",
        "median (65). Ding 2026 deliberately used serum creatinine rather than an eGFR formula because",
        "eGFR equations are unreliable in advanced heart failure (Discussion).",
        sep = " "
      ),
      source_name        = "Serum_creatinine"
    ),
    CYSC = list(
      description        = "Serum cystatin C",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 1.23 mg/L (2.5-97.5 percentiles 0.77-4.01). Enters BOTH CL and V as",
        "separate power models, each centred at 1.5 mg/L: `FCysC_CL = (CysC/1.5)**FCysC_CL1`",
        "(exponent -0.650) and `FCysc_V = (CysC/1.5)**FCysc_V1` (exponent -0.294). A sigmoid form",
        "(Online Resource Equation S4) was also tested on CL and gave no meaningful improvement while",
        "inflating the typical CL to 19.9 L/h, so the power model was retained.",
        sep = " "
      ),
      source_name        = "Cystatin_C"
    ),
    NTPROBNP = list(
      description        = "Serum N-terminal pro-B-type natriuretic peptide",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 1.095 ng/mL (2.5-97.5 percentiles 0.117-19.12). Power model on CL",
        "centred at 10 ng/mL: `FBNP_CL = (BNP/10)**FBNP_CL1`, exponent -0.0823 (Table 3). NOTE the",
        "unit: this paper reports NT-proBNP in ng/mL, whereas the common clinical reporting unit is",
        "pg/mL (1 ng/mL = 1000 pg/mL); a cohort median of 1.095 ng/mL is 1095 pg/mL. The assay's upper",
        "limit of detection is 35 ng/mL and results reported as > 35 ng/mL were set to 35 ng/mL before",
        "modelling (Results, Model development). This is the paper's headline covariate: it is the",
        "index of cardiac stress through which reduced cardiac function lowers vancomycin clearance.",
        sep = " "
      ),
      source_name        = "NT_proBNP"
    ),
    NEUT = list(
      description        = "Absolute neutrophil count",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 7.01 x 10^9/L (2.5-97.5 percentiles 2.06-19.76). Power model on V",
        "centred at 8.0 x 10^9/L: `FNEUT_V = (NEUT/8.0)**FNEUT_V1`, exponent -0.128 (Table 3). The",
        "control stream annotates the column `NEUT = Neutrophil ; 10^9/L`, which fixes the unit; the",
        "NEUT canonical's default unit is cells/mm^3 (numerically 1000x this value), so a user",
        "supplying cells/mm^3 must divide by 1000 before use with this model. Ding 2026 is the first",
        "report of neutrophil count as a covariate on vancomycin V (Discussion).",
        sep = " "
      ),
      source_name        = "Neutrophil"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Ding 2026 Table 1: median 56.9 years (2.5-97.5 percentiles 24.9-77.4); 41.3% of the cohort",
        "were elderly. Power model on V centred at 45 years: `FAge_V = (AGEY/45)**FAge_V1`, exponent",
        "+0.768 (Table 3), i.e. V rises with age. 45 years is the age of the paper's virtual standard",
        "patient, not the cohort median. Age was also tested as a categorical covariate on CL",
        "(< 60 vs >= 60 years) and rejected (dOFV 0.77, p = 0.38).",
        sep = " "
      ),
      source_name        = "Age"
    ),
    OCC = list(
      description        = "Occasion index for inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Decomposed inside model() into the indicators occ1..occ4 that select the per-occasion IOV",
        "etas on CL and V. Ding 2026 reports the two IOV magnitudes (Table 3, IOV_CL and IOV_V) and",
        "confirms exactly two IOV variances via the 2-degree-of-freedom OFV drop for Model 8 (Table",
        "2), but never states how many occasions were defined in the estimation dataset, nor how an",
        "occasion was delimited. Four occasions are encoded here, covering the cohort's median of",
        "three samples per patient; a user with more sampling occasions can extend the pattern by",
        "adding further etaiov_cl_<k> / etaiov_vc_<k> slots at the same fixed variances. Records",
        "outside 1-4 receive no IOV because all four indicators evaluate to zero. See the vignette",
        "Errata.",
        sep = " "
      ),
      source_name        = "OCC"
    ),
    ASSAY_CMIA = list(
      description        = paste(
        "Bioanalytical assay indicator: 1 = chemiluminescent microparticle immunoassay (CMIA,",
        "Abbott Architect i1000), 0 = enzyme-multiplied immunoassay technique (EMIT, Viva-E)",
        sep = " "
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (EMIT)",
      notes              = paste(
        "Per-observation indicator. Switches BOTH residual-error magnitudes. CMIA: proportional 0.198",
        "and additive 0.898 mg/L; EMIT: proportional FIXED at 0 and additive 1.83 mg/L (Table 3;",
        "control stream `$ERROR` block branching on `TYPE = Concentration_Method`, with the comments",
        "`IF (TYPE.EQ.0) THEN; CMIA i1000` and `ELSE; EMIT Viva-E`). The source column is therefore",
        "coded 0 = CMIA / 1 = EMIT and is INVERTED here to the assay-name-positive canonical",
        "(ASSAY_CMIA = 1 - Concentration_Method). CMIA carried 1062 of 1120 concentrations (94.8%);",
        "EMIT carried 58 (5.2%). Assay LOQs differ: 3.0 mg/L with a 3.0-100 mg/L linear range for",
        "CMIA, 2.0 mg/L with a 2.0-50 mg/L linear range for EMIT (Methods, Sampling and analysis).",
        sep = " "
      ),
      source_name        = "Concentration_Method"
    )
  )

  # Covariates Ding 2026 screened and explicitly rejected. Documented here so the
  # provenance of the covariate search is preserved without triggering a
  # "declared but not referenced" convention warning.
  covariatesDataExcluded <- list(
    LVEF = list(
      description = "Left ventricular ejection fraction (echocardiographic), as a fraction",
      units       = "(fraction)",
      type        = "continuous",
      notes       = paste(
        "Ding 2026 Table 1: median 0.56 (2.5-97.5 percentiles 0.29-0.69). Tested on CL as a",
        "continuous covariate (dOFV 1.074, df = 1, p = 0.3) and as a categorical covariate at",
        "cut-offs 0.4 (dOFV 0.68, p = 0.410) and 0.3 (no decrease in OFV). Not retained; NT-proBNP",
        "is the cardiac-function covariate that survived. The Discussion attributes the null result",
        "to insufficient data rather than to absence of an effect, and contrasts it with Shimamoto",
        "et al., who did find lower CL below LVEF 0.6.",
        sep = " "
      )
    ),
    CRRT = list(
      description = "Continuous renal replacement therapy indicator (1 = receiving CRRT)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Ding 2026 Table 1: 26 of 320 patients (8.1%). Tested as a categorical covariate on CL and",
        "rejected (dOFV 1.403, df = 1, p = 0.236; Online Resource Results, Covariate analysis).",
        sep = " "
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 320L,
    n_studies        = 1L,
    n_concentrations = 1120L,
    age_range        = "24.9-77.4 years (2.5-97.5 percentiles)",
    age_median       = "56.9 years",
    weight_range     = "41-101 kg (2.5-97.5 percentiles)",
    weight_median    = "65 kg",
    height_range     = "150-181 cm (2.5-97.5 percentiles)",
    height_median    = "169 cm",
    sex_female_pct   = 19.7,
    race_ethnicity   = "Not reported (single-centre Chinese cohort)",
    disease_state    = paste(
      "Adults (>= 18 years) hospitalised for cardiac surgery and treated with intravenous vancomycin,",
      "either for infective endocarditis (95 patients, 29.7%, who began vancomycin before surgery) or",
      "for post-surgical infection - pneumonia (62.2%), septicaemia (2.5%), mediastinal infection",
      "(1.6%) - or as prophylaxis for procedures implanting prosthetic material. Patients with tumours",
      "were excluded. 26 patients (8.1%) received continuous renal replacement therapy.",
      sep = " "
    ),
    renal_function   = paste(
      "Serum creatinine median 65 umol/L (2.5-97.5 percentiles 36-226.6); cystatin C median 1.23 mg/L",
      "(0.77-4.01). Renal function spans normal to severely impaired.",
      sep = " "
    ),
    cardiac_function = paste(
      "NT-proBNP median 1.095 ng/mL (2.5-97.5 percentiles 0.117-19.12); left ventricular ejection",
      "fraction median 0.56 (0.29-0.69).",
      sep = " "
    ),
    dose_range       = paste(
      "Intravenous vancomycin per routine clinical practice, initial dose per the label and adjusted",
      "on concentration where necessary (adjusted in 83 patients, 25.9%). The paper's dosage-strategy",
      "simulations span 250-1000 mg at 8-, 12- and 24-hour intervals as 2-hour infusions.",
      sep = " "
    ),
    regions          = "China (single centre: The First Affiliated Hospital of Soochow University, Suzhou)",
    notes            = paste(
      "Retrospective analysis of routine therapeutic-drug-monitoring records (Ding 2026 Table 1).",
      "1120 concentrations from 320 patients, median 3 samples per patient; concentration median",
      "13.77 mg/L (range 1.30-60.64). 542 samples (48.4%) were troughs, 566 (50.5%) were drawn between",
      "05:00 and 08:00, and 12 (1.1%) were drawn 1-41 min after the end of an infusion. Eleven CMIA",
      "concentrations (0.98%) fell below the 3.0 mg/L LOQ and were retained as continuous data.",
      "Fitted in NONMEM 7.6.0 with FOCE-I via Wings for NONMEM 760; evaluated by goodness-of-fit",
      "plots, visual predictive check and nonparametric bootstrap (Table 3).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters. Typical values apply to a subject at the standard
    # fat-free mass (56.1 kg) and at every covariate reference value
    # (Scr 80 umol/L, CysC 1.5 mg/L, NT-proBNP 10 ng/mL, NEUT 8.0 x 10^9/L,
    # age 45 years).
    # ------------------------------------------------------------------------
    lcl <- log(3.22); label("Clearance at reference covariates (L/h)")
    # Table 3 row "CL (L/h)" Original = 3.22 (bootstrap median 3.22, 2.5-97.5%
    # 2.98-3.50, RSE 3.9%); control stream $THETA "(0,3.22) ;CL_STD".
    lvc <- log(88.3); label("Distribution volume at reference covariates (L)")
    # Table 3 row "V (L)" Original = 88.3 (bootstrap median 87.9, 2.5-97.5%
    # 73.8-108, RSE 10.1%); control stream $THETA "(0,88.3) ;V_STD".

    # ------------------------------------------------------------------------
    # Theory-based allometric exponents on fat-free mass. Not estimated: Online
    # Resource Equations S8 and S9 write them as the theoretical 3/4 and 1, and
    # the control stream hardcodes them ("FSIZ_CL=(FFM/56.1)**0.75",
    # "FSIZ_V=FFM/56.1") rather than carrying them as $THETA records. Encoded
    # with fixed() accordingly.
    # ------------------------------------------------------------------------
    e_ffm_cl <- fixed(0.75); label("Allometric exponent on (FFM / 56.1) for CL (unitless)")
    # Online Resource Equation S8; control stream $PK "FSIZ_CL=(FFM/56.1)**0.75".
    e_ffm_vc <- fixed(1);    label("Allometric exponent on (FFM / 56.1) for V (unitless)")
    # Online Resource Equation S9 (V scales linearly on NFM); control stream $PK
    # "FSIZ_V=FFM/56.1".

    # ------------------------------------------------------------------------
    # Covariate effects. All are power-model exponents (Online Resource
    # Equation S3, P = P_STD * (Cov / Mean_Cov)^FCov_P). The centring constants
    # live in model() because the control stream hardcodes them there.
    # ------------------------------------------------------------------------
    e_creat_cl    <- -0.458;  label("Power exponent on (CREAT / 80 umol/L) for CL (unitless)")
    # Table 3 row "FScr_CL" Original = -0.458 (bootstrap -0.546 to -0.366, RSE
    # 10.1%); control stream $THETA "-0.458 ;FScr_CL1".
    e_cysc_cl     <- -0.650;  label("Power exponent on (CYSC / 1.5 mg/L) for CL (unitless)")
    # Table 3 row "FCysc_CL" Original = -0.650 (bootstrap -0.765 to -0.528, RSE
    # 9.6%); control stream $THETA "-0.65 ;FCysC_CL1".
    e_ntprobnp_cl <- -0.0823; label("Power exponent on (NTPROBNP / 10 ng/mL) for CL (unitless)")
    # Table 3 row "FNT-proBNP _CL" Original = -0.0823 (bootstrap -0.115 to
    # -0.0507, RSE 20.1%); control stream $THETA "-0.0823 ;FBNP_CL1".
    e_neut_vc     <- -0.128;  label("Power exponent on (NEUT / 8.0 x 10^9/L) for V (unitless)")
    # Table 3 row "FNEUT_V" Original = -0.128 (bootstrap -0.225 to -0.0380, RSE
    # 36.2%); control stream $THETA "-0.128 ;FNEUT_V1".
    e_cysc_vc     <- -0.294;  label("Power exponent on (CYSC / 1.5 mg/L) for V (unitless)")
    # Table 3 row "FCysc_V" Original = -0.294 (bootstrap -0.455 to -0.195, RSE
    # 22.1%); control stream $THETA "-0.294 ;FCysc_V1".
    e_age_vc      <-  0.768;  label("Power exponent on (AGE / 45 years) for V (unitless)")
    # Table 3 row "FAge_V" Original = 0.768 (bootstrap 0.350 to 1.14, RSE
    # 27.2%); control stream $THETA "0.768 ;FAge_V1".

    # ------------------------------------------------------------------------
    # Inter-individual variability, exponential (Online Resource Equation S1).
    # The control stream $OMEGA carries the VARIANCES 0.0675 and 0.239, whose
    # square roots (0.2598 and 0.4889) are exactly the Table 3 rows "IIV_CL"
    # 0.260 and "IIV_V" 0.489 and the Abstract's "26.0%" and "48.9%". Table 3
    # therefore reports apparent CVs (i.e. omega on the SD scale), and the
    # variances below are taken from the control stream directly rather than
    # back-transformed. Diagonal: the stream declares two separate $OMEGA
    # records, not a BLOCK, so CL and V IIV are uncorrelated.
    # ------------------------------------------------------------------------
    etalcl ~ 0.0675
    # Control stream $OMEGA "0.0675 ;BSV_CL"; sqrt = 0.2598 -> Table 3 IIV_CL 0.260.
    etalvc ~ 0.239
    # Control stream $OMEGA "0.239 ;BSV_V"; sqrt = 0.4889 -> Table 3 IIV_V 0.489.

    # ------------------------------------------------------------------------
    # Inter-occasion variability, exponential (Online Resource Equation S1, the
    # kappa term). Table 3 reports IOV_CL = 0.089 and IOV_V = 0.136 on the same
    # SD scale as the IIV rows just above (Abstract: "the corresponding
    # interoccasion variabilities were 8.9% and 13.6%"), so the variances are
    # 0.089^2 = 0.007921 and 0.136^2 = 0.018496. Adding IOV to Model 7 dropped
    # the OFV by 9.68 on 2 df (Table 2, Model 8), confirming exactly two IOV
    # variances.
    #
    # Implemented by the occasion-indicator expansion rather than rxode2's
    # `~ var | OCC` syntax, which parses but cannot be simulated from an rxUi.
    # This is also the more faithful translation of a NONMEM $PK IOV block.
    # Occasions after the first repeat the same variance, the analogue of
    # NONMEM's $OMEGA BLOCK(1) SAME, and are therefore fixed().
    #
    # The paper does not report how many occasions were defined; four are
    # provided here, which covers the cohort's median of three samples per
    # patient. See the vignette Errata.
    # ------------------------------------------------------------------------
    etaiov_cl_1 ~ 0.007921
    etaiov_cl_2 ~ fixed(0.007921)
    etaiov_cl_3 ~ fixed(0.007921)
    etaiov_cl_4 ~ fixed(0.007921)
    # Table 3 row "IOV_CL" Original = 0.089 (bootstrap 0.001-0.141, RSE 34.6%).
    etaiov_vc_1 ~ 0.018496
    etaiov_vc_2 ~ fixed(0.018496)
    etaiov_vc_3 ~ fixed(0.018496)
    etaiov_vc_4 ~ fixed(0.018496)
    # Table 3 row "IOV_V" Original = 0.136 (bootstrap 0.001-0.265, RSE 59.4%).

    # ------------------------------------------------------------------------
    # Residual error, combined additive + proportional on the variance scale
    # (Online Resource Equation S2:
    #   Y = Con + sqrt(Con^2 * theta_PROP^2 + theta_ADD^2) * eps, eps ~ N(0,1)),
    # which is exactly nlmixr2's `add() + prop()`. $SIGMA is 1 FIX, so the
    # thetas below are the residual SDs themselves. Magnitudes differ by
    # bioanalytical assay and are selected per observation in model() by
    # ASSAY_CMIA.
    # ------------------------------------------------------------------------
    propSdCmia <- 0.198; label("Proportional residual error, CMIA assay (fraction)")
    # Table 3 row "RUV_PROP1" Original = 0.198 (bootstrap 0.149-0.235, RSE
    # 11.5%); control stream $THETA "(0,0.198) ;RUV_PROP1". Results text:
    # "The proportional residual error and additional residual error of CMIA
    # were 19.8% and 0.898 mg/L".
    addSdCmia  <- 0.898; label("Additive residual error, CMIA assay (mg/L)")
    # Table 3 row "RUV_ADD1 (mg/L)" Original = 0.898 (bootstrap 0.306-1.20, RSE
    # 25.0%); control stream $THETA "(0,0.898) ;RUV_ADD1".
    propSdEmit <- fixed(0); label("Proportional residual error, EMIT assay (fraction)")
    # Control stream $THETA "0 FIX ;RUV_PROP2". Results text: "The 95% CI lower
    # bound of the proportional residual error for EMIT was close to 0;
    # therefore, it was fixed at 0." Table 3 has no RUV_PROP2 row, consistent
    # with the parameter being fixed rather than estimated.
    addSdEmit  <- 1.83;  label("Additive residual error, EMIT assay (mg/L)")
    # Table 3 row "RUV_ADD2 (mg/L)" Original = 1.83 (bootstrap 0.718-2.46, RSE
    # 24.0%); control stream $THETA "(0,1.83) ;RUV_ADD2". Results text: "The
    # additional residual error of EMIT was 1.83 mg/L".
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Derived covariate terms.
    #
    # Fat-free mass by the Janmahasatian height-squared form, Online Resource
    # Equation S6, transcribed from the control stream $PK block:
    #   HTM  = HTCM/100
    #   FFM  = WHSMAX * HTM**2 * WT / (WHS50 * HTM**2 + WT)
    # with WHSMAX 42.92 kg/m^2 and WHS50 30.93 kg/m^2 for males, and 37.99 and
    # 35.98 for females. The stream selects the branch on M1F0 (1 = male); with
    # the canonical SEXF = 1 - M1F0, SEXF interpolates to the female constants.
    #
    # Normal fat mass (Equation S7, NFM = FFM + Ffat * (WT - FFM)) collapses to
    # FFM because Ffat was fixed at 0 for both CL and V: its 95% bootstrap CIs
    # spanned zero (-0.457 to 1.981 for CL, -0.862 to 1.47 for V), so fat mass
    # explained no variability and FFM alone was the size descriptor (Results,
    # Model development). Accordingly the control stream scales on FFM/56.1 with
    # no Ffat term, and so does this model.
    #
    # 56.1 kg is the standard fat-free mass, i.e. Equation S6 evaluated for the
    # paper's standard male subject of 70 kg and 176 cm:
    #   42.92 * 1.76^2 * 70 / (30.93 * 1.76^2 + 70) = 56.13.
    # ----------------------------------------------------------------------
    htm    <- HT / 100
    whsmax <- 42.92 + SEXF * (37.99 - 42.92)
    whs50  <- 30.93 + SEXF * (35.98 - 30.93)
    ffm    <- whsmax * htm^2 * WT / (whs50 * htm^2 + WT)

    # Inter-occasion variability terms, selected by the occasion index OCC.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 +
      occ3 * etaiov_cl_3 + occ4 * etaiov_cl_4
    iov_vc <- occ1 * etaiov_vc_1 + occ2 * etaiov_vc_2 +
      occ3 * etaiov_vc_3 + occ4 * etaiov_vc_4

    # ----------------------------------------------------------------------
    # 2. Individual PK parameters, transcribed from the control stream $PK:
    #   GRPCL = CL_STD * FSIZ_CL * FScr_CL * FCysC_CL * FBNP_CL
    #   GRPV  = V_STD  * FSIZ_V  * FNEUT_V * FCysc_V  * FAge_V
    # with the covariate factors written as power models about their hardcoded
    # centring constants:
    #   FScr_CL  = (Scr/80)**FScr_CL1        FNEUT_V = (NEUT/8.0)**FNEUT_V1
    #   FCysC_CL = (CysC/1.5)**FCysC_CL1     FCysc_V = (CysC/1.5)**FCysc_V1
    #   FBNP_CL  = (BNP/10)**FBNP_CL1        FAge_V  = (AGEY/45)**FAge_V1
    # ----------------------------------------------------------------------
    cl <- exp(lcl + etalcl + iov_cl) *
      (ffm / 56.1)^e_ffm_cl *
      (CREAT / 80)^e_creat_cl *
      (CYSC / 1.5)^e_cysc_cl *
      (NTPROBNP / 10)^e_ntprobnp_cl
    vc <- exp(lvc + etalvc + iov_vc) *
      (ffm / 56.1)^e_ffm_vc *
      (NEUT / 8.0)^e_neut_vc *
      (CYSC / 1.5)^e_cysc_vc *
      (AGE / 45)^e_age_vc

    # 3. Micro-constant.
    kel <- cl / vc

    # ----------------------------------------------------------------------
    # 4. ODE system. The source uses the closed-form $SUBROUTINE ADVAN1 TRANS2
    #    (one compartment, IV, parameterised on CL and V); the equivalent ODE
    #    is written out explicitly for house consistency.
    # ----------------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ----------------------------------------------------------------------
    # 5. Observation and error. S1 = V in the stream, so dose in mg over volume
    #    in L gives mg/L. The $ERROR block branches on the assay to pick the
    #    residual magnitudes; the equivalent linear switch on ASSAY_CMIA is
    #    used here (ASSAY_CMIA = 1 for CMIA, 0 for EMIT).
    # ----------------------------------------------------------------------
    Cc <- central / vc
    propSd <- propSdCmia * ASSAY_CMIA + propSdEmit * (1 - ASSAY_CMIA)
    addSd  <- addSdCmia  * ASSAY_CMIA + addSdEmit  * (1 - ASSAY_CMIA)
    Cc ~ add(addSd) + prop(propSd)
  })
}
