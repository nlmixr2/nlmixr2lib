Liang_2026_remimazolam <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for remimazolam",
    "(administered as remimazolam tosilate) in healthy Chinese adults,",
    "developed on pooled arterial-plasma data from the HR7056-Ia single",
    "ascending dose and HR7056-Ib loading-plus-maintenance-infusion Phase I",
    "studies. Remimazolam is an ultra-short-acting benzodiazepine that",
    "carboxylesterase 1 (CES1) hydrolyses to the inactive, renally excreted",
    "acid metabolite CNS7054, so no active metabolite is tracked. Body",
    "weight enters every clearance and volume term as theory-based",
    "allometric scaling with exponents fixed to 0.75 and 1, which is what",
    "lets the adult model be extrapolated to children and adolescents; no",
    "other covariate survived stepwise selection. The paper's purpose was",
    "pediatric dose selection by adult exposure matching, and the",
    "recommended regimens were approved by China's Center for Drug",
    "Evaluation. A companion whole-body PBPK model was built in PK-Sim v12",
    "and is NOT reproduced here -- the paper publishes only a six-row",
    "drug-specific parameter table for it and takes all system physiology,",
    "organ partitioning and CES1 tissue expression from the platform's",
    "built-in libraries, so no ODE system is recoverable from any on-disk",
    "source. See the vignette Errata."
  )
  reference <- paste(
    "Liang QY, Hu HH, Djebli N, Huang YY, Jiang H (2026). Dose",
    "Recommendation of Remimazolam Tosilate for General Anesthesia in",
    "Children and Adolescents: Synergistic Combination of PopPK and PBPK",
    "Approaches. Pharmaceutics 18(3):315.",
    "doi:10.3390/pharmaceutics18030315.",
    "Structural model and final parameter estimates from Table 4;",
    "allometric scaling from Equation 3; residual-error and random-effects",
    "structure, the 62 kg allometric reference weight and the ng/mL",
    "concentration scaling from the final NONMEM control stream in",
    "Supplementary Material Section S2; covariate screen from Table S3;",
    "baseline demographics from Table S1; validation targets from Tables 5",
    "and 7."
  )
  vignette <- "Liang_2026_remimazolam"
  units <- list(time = "min", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model, entered as",
        "theory-based allometric scaling on all three volumes and all three",
        "clearances (Equation 3) with the exponents fixed rather than",
        "estimated. Adding it to the base model dropped the OFV by 31.557",
        "(Discussion). Reference weight is 62 kg: the final control stream",
        "(Supplementary Material S2) writes every term as",
        "(WEIGHTBL/62)**THETA and labels its thetas 'L/min/62 kg', whereas",
        "the Table 4 row headers read '63 kg' -- see the vignette Errata.",
        "Analysis-population weight was median 62.8 kg, range 50.2-83.8 kg",
        "(Table S1). The paper's pediatric extrapolation drives this term",
        "down to a cohort median of 6.6 kg (Table S4), so the model is",
        "routinely evaluated far below the weight range it was fit on."
      ),
      source_name        = "WEIGHTBL"
    )
  )

  # Covariates the paper screened and rejected. None of these is referenced in
  # model(); they are recorded so the provenance of the covariate screen
  # survives. Table S3 lists which parameter each was tested on, and shows
  # that the full model (V2+sex, CL+BUN, V1+ALB) lost all three terms during
  # backward elimination, leaving allometric weight alone in the final model.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "The strongest candidate covariate: tested on V2, it entered the",
        "full model first (dOFV -10.742, p = 0.001047) but was then removed",
        "during backward elimination because re-adding it raised the OFV by",
        "only 10.742, short of the 10.83 retention threshold (Table S3).",
        "Also screened on CL, Q2, Q3, V1 and V3. Cohort was 15 female",
        "(21.1%) of 71 (Table S1)."
      ),
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on CL, Q2, Q3, V1, V2 and V3 (Table S3 note); never",
        "reached the forward-inclusion threshold. Cohort median 27.0 years,",
        "range 18.0-51.0 (Table S1). Age is the driver of the PBPK model's",
        "CES1 ontogeny function instead (see the FCES1 entry)."
      ),
      source_name        = "AGE"
    ),
    HT = list(
      description        = "Body height at baseline.",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on all six disposition parameters (Table S3 note) and not",
        "retained. Correlated with weight and BSA at R > 0.8, so the",
        "screening rules kept only one member of that group (Methods 2.6.2).",
        "Cohort median 168.0 cm, range 155.0-184.0 (Table S1)."
      ),
      source_name        = "HEIGHTBL"
    ),
    BMI = list(
      description        = "Body mass index at baseline.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Listed among the assessed covariates in Methods 2.6.2 and derived",
        "via Equation 7, but it does not appear in the Table S3",
        "parameter-by-parameter screen, so it was dropped at the |R| >= 0.3",
        "correlation pre-screen. Cohort median 22.3 kg/m^2, range",
        "18.4-25.7 (Table S1)."
      ),
      source_name        = "BMIBL"
    ),
    BSA = list(
      description        = "Body surface area at baseline.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Assessed per Methods 2.6.2 and carried in the control stream's",
        "$INPUT as BSABL, but excluded as a body-size descriptor because it",
        "correlates with height and weight at R > 0.8 (Methods 2.6.2).",
        "Still required upstream of the analysis as an input to the MDRD",
        "eGFR calculation (Equation 6)."
      ),
      source_name        = "BSABL"
    ),
    ALT = list(
      description        = "Alanine aminotransferase.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on all six disposition parameters (Table S3 note) and not",
        "retained. Cohort median 13.0 U/L, range 5.0-33.0 (Table S1) -- a",
        "healthy-volunteer range with no hepatic impairment, so the screen",
        "had little signal to find."
      ),
      source_name        = "ALTBL"
    ),
    AST = list(
      description        = "Aspartate aminotransferase.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on all six disposition parameters (Table S3 note) and not",
        "retained; correlated with ALT at R > 0.8 (Methods 2.6.2). Cohort",
        "median 19.0 U/L, range 11.0-32.0 (Table S1)."
      ),
      source_name        = "ASTBL"
    ),
    TBILI = list(
      description        = "Total serum bilirubin.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on all six disposition parameters (Table S3 note) and not",
        "retained. Cohort median 11.2 umol/L, range 5.3-26.2 (Table S1);",
        "reported in SI units, so no mg/dL conversion is needed."
      ),
      source_name        = "TBILBL"
    ),
    ALB = list(
      description        = "Serum albumin.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered the full model on V1 (dOFV -4.612, p = 0.031753) and was",
        "the first term removed during backward elimination (Table S3).",
        "Also screened on CL, Q2, Q3, V2 and V3. Mechanistically plausible",
        "as remimazolam is ~92% bound to albumin (Methods 2.5.1), but the",
        "effect did not survive the p < 0.001 retention criterion. No",
        "albumin summary is tabulated in Table S1 despite ALB being listed",
        "among the assessed covariates."
      ),
      source_name        = "ALBBL"
    ),
    ALP = list(
      description        = "Alkaline phosphatase.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on all six disposition parameters (Table S3 note) and not",
        "retained. No summary tabulated in Table S1."
      ),
      source_name        = "ALPBL"
    ),
    CREAT = list(
      description        = "Serum creatinine.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Assessed per Methods 2.6.2 but excluded from the",
        "parameter-by-parameter screen because it correlates with eGFR and",
        "sex at R > 0.8; eGFR was carried forward instead. Cohort median",
        "72.0 umol/L, range 46.0-100.0 (Table S1). Used upstream as the",
        "input to the MDRD eGFR calculation (Equation 6)."
      ),
      source_name        = "CRBL"
    ),
    BUN = list(
      description        = "Blood urea nitrogen.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Entered the full model on CL (dOFV -8.241, p = 0.004096) and was",
        "removed during backward elimination (Table S3). Also screened on",
        "Q2, Q3. Cohort median 4.5 mmol/L, range 1.9-8.4 (Table S1)."
      ),
      source_name        = "BUNBL"
    ),
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate from the four-variable MDRD",
        "equation (Equation 6)."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on CL, Q2 and Q3 (Table S3 note) and not retained --",
        "expected, since remimazolam is cleared by CES1 hydrolysis and only",
        "the inactive metabolite is renally excreted (Introduction). Two",
        "cautions if this column is ever reused. (1) Units are absolute",
        "mL/min, not the canonical mL/min/1.73 m^2: Equation 6 multiplies",
        "the BSA-normalised MDRD result by BSA/1.73. (2) The tabulated",
        "values are not reproducible from Equation 6 as printed -- the",
        "median 62.8 kg, 168 cm, 27-year-old male with creatinine 72",
        "umol/L evaluates to about 124 mL/min, against a tabulated cohort",
        "median of 48.3 mL/min (range 27.9-107.7, Table S1), which is",
        "implausibly low for healthy young adults. Because the covariate",
        "was rejected, neither issue affects any simulation from this model."
      ),
      source_name        = "EGFRBL"
    ),
    FCES1 = list(
      description        = paste(
        "Fractional hepatic CES1 protein abundance relative to the adult",
        "maximum, derived from age by the sigmoidal ontogeny function of",
        "Boberg et al."
      ),
      units              = "(fraction of adult)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on CL only (Methods 2.6.2) and not retained: it changed",
        "the OFV negligibly and did not improve diagnostics, which is",
        "expected because the analysis population was entirely adult and",
        "therefore CES1-mature (Discussion). Computed as",
        "FCES1 = (Adultmax - Fbirth) / (Age50^n + Age^n) * Age^n + Fbirth",
        "with Fbirth = 0.20, Adultmax = 1, Age50 = 1.10 years and n = 0.56",
        "(Equation 2 and Supplementary Material S1), so it reaches about",
        "0.9 of adult abundance by 3 years of age. Cohort values were",
        "median 0.886, range 0.862-0.916 (Table S1) -- a 5% spread, which",
        "is why adult data cannot inform it. Retained in the PBPK model,",
        "where age-dependent CES1 abundance is the mechanism driving the",
        "lower exposures predicted below about 30 kg. Not a register",
        "canonical: it is a derived age function rather than an observed",
        "column, and it is excluded here, so no new covariate name is",
        "being claimed."
      ),
      source_name        = "FCES1"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "remimazolam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 71,
    n_studies      = 2,
    age_range      = "18.0-51.0 years",
    age_median     = "27.0 years",
    weight_range   = "50.2-83.8 kg",
    weight_median  = "62.8 kg",
    sex_female_pct = 21.1,
    race_ethnicity = c(Asian = 100),
    disease_state  = "healthy volunteers",
    dose_range     = paste(
      "0.007-0.32 mg/kg as a 1 min intravenous injection (HR7056-Ia,",
      "eleven single ascending dose groups); 0.29 mg/kg loading dose",
      "infused over 1 min plus 1.08 mg/kg/h maintenance infusion for 2 h",
      "(HR7056-Ib, two-period crossover). Doses are expressed as the free",
      "base."
    ),
    regions        = "China",
    notes          = paste(
      "Baseline demographics from Table S1; study designs from Table 1.",
      "Pooled from the single-centre Phase I studies HR7056-Ia (63",
      "subjects, 1197 planned samples) and HR7056-Ib (8 subjects, 416",
      "planned samples), contributing 1439 arterial and 518 venous",
      "concentrations in total. The PopPK model was fit to the ARTERIAL",
      "plasma concentrations only (Methods 2.2); the venous data were used",
      "for the companion PBPK model. Marked arterio-venous concentration",
      "gradients were observed, with venous peaks well below the",
      "contemporaneous arterial Cmax (Results 3.1), so this model predicts",
      "arterial concentrations and will overpredict a venous sample around",
      "the peak. HR7056-Ib subjects also received 0.5 mg flumazenil or",
      "saline 1 h 55 min into dosing; flumazenil is a GABA-A receptor",
      "antagonist with no reported effect on remimazolam disposition and",
      "the model carries no term for it. Assay range 0.5-1000 ng/mL with",
      "an LLOQ of 0.5 ng/mL (Methods 2.3), so observed concentrations at",
      "the higher doses exceeded the calibration range and were presumably",
      "assayed on dilution. No pediatric data were used: the",
      "3-18-year-old predictions the paper reports come from extrapolating",
      "this adult model through the allometric weight term (Methods 2.6.4)."
    )
  )

  ini({
    # Typical values are the estimates of the final model, Table 4 'Final
    # Model / Estimates' column, each at the 62 kg allometric reference
    # weight. The control stream's $THETA block (Supplementary Material S2)
    # holds the INITIAL estimates (1.03, 2.1, 11.1, 1.51, 19.7, 0.264), which
    # differ slightly and are not used here.
    lcl  <- log(1.03)  ; label("Clearance CL at 62 kg (L/min)")                                                        # Table 4: CL 1.03 (RSE 1.9%), 95% CI 0.992-1.068
    lvc  <- log(2.08)  ; label("Central volume of distribution V1 at 62 kg (L)")                                       # Table 4: V1 2.08 (RSE 2.6%), 95% CI 1.97-2.19
    lvp  <- log(10.9)  ; label("First peripheral volume of distribution V2 at 62 kg (L)")                              # Table 4: V2 10.9 (RSE 3.8%), 95% CI 10.1-11.7
    lq   <- log(1.49)  ; label("Intercompartmental clearance Q2 (central <-> peripheral1) at 62 kg (L/min)")           # Table 4: Q2 1.49 (RSE 4.1%), 95% CI 1.37-1.61
    lvp2 <- log(19.7)  ; label("Second peripheral volume of distribution V3 at 62 kg (L)")                             # Table 4: V3 19.7 (RSE 4.1%), 95% CI 18.1-21.3
    lq2  <- log(0.266) ; label("Intercompartmental clearance Q3 (central <-> peripheral2) at 62 kg (L/min)")           # Table 4: Q3 0.266 (RSE 6.1%), 95% CI 0.234-0.298

    # Theory-based allometric exponents. All six were fixed rather than
    # estimated (Methods 2.6.1 after Equation 3; Table 4 footnote), and carry
    # the FIX flag as THETA(7)-THETA(12) of the control stream.
    e_wt_cl  <- fixed(0.75) ; label("Allometric exponent of body weight on CL (unitless)")               # Equation 3: 0.75 for CL
    e_wt_vc  <- fixed(1)    ; label("Allometric exponent of body weight on V1 (unitless)")               # Equation 3: 1.0 for V1
    e_wt_vp  <- fixed(1)    ; label("Allometric exponent of body weight on V2 (unitless)")               # Equation 3: 1.0 for V2
    e_wt_q   <- fixed(0.75) ; label("Allometric exponent of body weight on Q2 (unitless)")               # Equation 3: 0.75 for Q2
    e_wt_vp2 <- fixed(1)    ; label("Allometric exponent of body weight on V3 (unitless)")               # Equation 3: 1.0 for V3
    e_wt_q2  <- fixed(0.75) ; label("Allometric exponent of body weight on Q3 (unitless)")               # Equation 3: 0.75 for Q3

    # Exponential IIV (Equation 4) on all six disposition parameters, with two
    # 2x2 OMEGA blocks matching $OMEGA BLOCK(2) twice then two diagonal
    # elements in the control stream. Values are the Table 4 final-model
    # variances and covariances.
    #
    # Block 1, CL with V1. Table 4 rows 'omega^2 (CL)' = 0.0203,
    # 'omega (CL): omega (V1)' = 0.00948 and 'omega^2 (V1)' = 0.00546. The
    # paper flags IIV on V1 as negligible and imprecise (RSE 99.5%, 95% CI
    # spanning zero) but retained the block because eta-CL and eta-V1
    # correlate and the block improved the OFV (Results 3.3). The implied
    # correlation is 0.90 and the block is positive definite
    # (determinant 2.10e-5), so it is encoded as published.
    etalcl + etalvc ~ c(0.0203,
                        0.00948, 0.00546)

    # Block 2, V2 with Q2. Table 4 rows 'omega^2 (V2)' = 0.0509,
    # 'omega (V2): omega (Q2)' = 0.0555 and 'omega^2 (Q2)' = 0.107; implied
    # correlation 0.75. Note that Equation 4's prose list of parameters
    # carrying IIV omits Q2, but Table 4 reports its variance with good
    # precision (RSE 18.3%, shrinkage 5%) and the control stream declares
    # ETA(4) on Q2, so the prose list is incomplete.
    etalvp + etalq ~ c(0.0509,
                       0.0555, 0.107)

    etalvp2 ~ 0.072                                                                                      # Table 4: omega^2 (V3) = 0.072 (RSE 17.5%, shrinkage 7%)
    etalq2  ~ 0.0777                                                                                     # Table 4: omega^2 (Q3) = 0.0777 (RSE 27.9%, shrinkage 12.2%)

    # Residual error. The paper's Equation 5 prints a plain additive error on
    # the concentration scale, but the control stream shows the fit was run on
    # LOG-TRANSFORMED concentrations: its $PROBLEM description reads
    # 'logDV+ADD', $ERROR sets IPRED = LOG(F) and returns
    # Y = IPRED * (1 + ERR(1)) + ERR(2), and $SIGMA fixes the proportional
    # term to 0 leaving only the additive one. Additive on the log scale is
    # nlmixr2's `Cc ~ lnorm(expSd)`, i.e. approximately proportional error in
    # linear space. Encoding Equation 5 literally would give an additive SD of
    # 0.127 ng/mL against concentrations up to 5790 ng/mL, which the
    # goodness-of-fit plots (Figure S3, drawn on log axes) contradict.
    expSd <- sqrt(0.0162) ; label("Log-scale residual SD (~CV 12.8%)")                                   # Table 4: sigma^2 (ADD) = 0.0162 (RSE 8.5%, shrinkage 11%)
  })

  model({
    # Allometric reference weight. Table 4's row headers read '63 kg', but the
    # final control stream centres every term on WEIGHTBL/62 and labels its
    # own thetas 'L/min/62 kg'; the analysis-population median was 62.8 kg
    # (Table S1). The control stream is authoritative for the centring
    # constant, so 62 is used here -- see the vignette Errata. The difference
    # is a 1.3% scale factor on the volumes and 1.0% on the clearances.
    wt_ref <- 62

    # Individual disposition parameters, Equation 3 combined with Equation 4.
    cl  <- exp(lcl  + etalcl)  * (WT / wt_ref)^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * (WT / wt_ref)^e_wt_vc
    vp  <- exp(lvp  + etalvp)  * (WT / wt_ref)^e_wt_vp
    q   <- exp(lq   + etalq)   * (WT / wt_ref)^e_wt_q
    vp2 <- exp(lvp2 + etalvp2) * (WT / wt_ref)^e_wt_vp2
    q2  <- exp(lq2  + etalq2)  * (WT / wt_ref)^e_wt_q2

    # Three-compartment linear disposition with intravenous input, written out
    # explicitly. This is the ODE form of the source model's ADVAN11 TRANS4
    # (CL / V1 / Q2 / V2 / Q3 / V3) parameterisation. Dose intravenously into
    # 'central'.
    d/dt(central)     <- -cl * central / vc -
                          q  * central / vc + q  * peripheral1 / vp -
                          q2 * central / vc + q2 * peripheral2 / vp2
    d/dt(peripheral1) <-  q  * central / vc - q  * peripheral1 / vp
    d/dt(peripheral2) <-  q2 * central / vc - q2 * peripheral2 / vp2

    # Dose in mg with vc in L gives mg/L; multiply by 1000 for the ng/mL of
    # the published concentrations and of the 0.5-1000 ng/mL assay. This
    # reproduces the control stream's S1 = V1/1000 scaling.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
