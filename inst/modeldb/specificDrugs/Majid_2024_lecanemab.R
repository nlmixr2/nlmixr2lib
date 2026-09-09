Majid_2024_lecanemab <- function() {
  description <- paste(
    "Population pharmacokinetics of the anti-amyloid-beta protofibril",
    "monoclonal antibody lecanemab (Leqembi) in 1619 subjects with early",
    "Alzheimer's disease (mild cognitive impairment due to AD or mild AD",
    "dementia), from 21,929 serum concentrations pooled across two phase I",
    "studies (101, 104), the phase II study 201 Core and open-label",
    "extension, and the phase III Clarity AD study 301 Core and open-label",
    "extension. Linear two-compartment model with first-order elimination",
    "from the central compartment following a 1 h intravenous infusion,",
    "parameterized for CL, V1, V2 and Q. Covariate effects: body weight",
    "(power, reference 72 kg) and albumin (power, reference 43 g/L) on CL;",
    "female sex and sample-level ADA-positive status as multiplicative",
    "ratios on CL; body weight (power), female sex and Japanese",
    "race/ethnicity as ratios on V1; Japanese race/ethnicity as a ratio on",
    "V2. Q carries no covariate and no IIV. A manufacturing-process",
    "comparability factor F is applied to the intravenous dose: F is fixed",
    "at 1 for Process A and estimated at 0.904 for Process B, and the",
    "between-subject variability on F applies to Process B records ONLY",
    "(supplement Text S1 $PK codes it inside an IF (FORM.EQ.1) branch).",
    "CL and V1 IIV are correlated (R = 0.144). None of the retained",
    "covariates shifted steady-state AUC or Cmax outside the 0.8-1.25",
    "acceptance interval, and age was not a significant covariate. The",
    "typical terminal half-life is 14.5 days. Companion exposure-response",
    "model: Majid_2024_lecanemab_ariae.",
    sep = " "
  )
  reference <- paste(
    "Majid O, Cao Y, Willis BA, Hayato S, Takenaka O, Lalovic B,",
    "Sreerama Reddy SH, Penner N, Reyderman L, Yasuda S, Hussein Z (2024).",
    "Population pharmacokinetics and exposure-response analyses of safety",
    "(ARIA-E and isolated ARIA-H) of lecanemab in subjects with early",
    "Alzheimer's disease.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(12):2111-2123.",
    "doi:10.1002/psp4.13224.",
    sep = " "
  )
  vignette <- "Majid_2024_lecanemab"

  units <- list(
    time          = "h",
    dosing        = "mg (intravenous infusion into central over 60 +/- 10 minutes in every contributing study; weight-based mg/kg regimens, so a 10 mg/kg dose for a 72 kg subject is 720 mg)",
    concentration = "ug/mL (equivalently mg/L; dose in mg divided by a volume in L)"
  )

  compartmentData <- list(
    central     = list(analyte = "lecanemab", units = "mg", specimen = "serum", verified = FALSE),
    peripheral1 = list(analyte = "lecanemab", units = "mg", specimen = "serum", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power function normalized to 72 kg on BOTH CL and V1,",
        "with separately estimated exponents (0.353 and 0.513) rather than",
        "the fixed 0.75 / 1 allometric pair. The 72 kg normalization",
        "constant is the population median, confirmed by Table S2 (PK",
        "analysis set median 72 kg, range 37.7-130.5) and named as the",
        "reference subject weight in the Figure 1 forest-plot caption. The",
        "Figure 1 test categories 49 and 99 kg are the 5th and 95th",
        "percentiles of the PK analysis set. Because lecanemab is dosed in",
        "mg/kg, dose rises linearly with weight while CL rises only as",
        "WT^0.353, so heavier subjects have higher exposure; the paper's",
        "conclusion is nonetheless that the net effect stays inside the",
        "0.8-1.25 acceptance interval."
      ),
      source_name        = "WGT (supplement Text S1 $INPUT); BW (printed CL and V1 equations)"
    ),
    ALB = list(
      description        = "Serum albumin concentration at the time of the PK sample.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a power function normalized to 43 g/L on CL with a",
        "NEGATIVE exponent (-0.374): clearance DECLINES as albumin rises.",
        "The paper's Discussion states the sign in words ('Lecanemab",
        "clearance was found to decline with increasing albumin levels with",
        "an exponent of 0.374'), which reads as a positive number only",
        "because the direction is carried by the verb; Table 1 and the",
        "printed CL equation both give -0.374 and are authoritative. The",
        "mechanistic rationale offered is FcRn-mediated recycling, which",
        "handles IgG and albumin jointly, so a higher albumin indicates more",
        "FcRn and slower lecanemab elimination. The 43 g/L normalization",
        "constant is the population median (Table S2: median 43, range",
        "35-54); Figure 1 tests 39 and 48 g/L as the 5th and 95th",
        "percentiles. Register units are g/L (SI), which is exactly what",
        "this paper reports -- no conversion needed."
      ),
      source_name        = "ALB"
    ),
    SEXF = list(
      description        = "Female sex indicator; 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "The source column SEXN is already coded 0 = male / 1 = female (the",
        "footnote under the printed equations states 'SEX, 0 (male) or 1",
        "(female)'), so SEXF equals it directly with no 1 - x inversion and",
        "no sign change on the coefficients. Applied as a ratio raised to",
        "the indicator on both CL (0.791, i.e. 20.9% lower in females) and",
        "V1 (0.868, i.e. 13.2% lower in females); the two percentages are",
        "quoted in exactly this form in the Discussion, which confirms that",
        "males are the reference. PK analysis set: 800 females (49.4%), 819",
        "males (50.6%) (Table S2)."
      ),
      source_name        = "SEXN"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positive status; 1 = ADA-positive, 0 = ADA-negative.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = paste(
        "TIME-VARYING at the SAMPLE level, not a per-subject flag: the",
        "Methods specify 'ADA status at sample level (positive or",
        "negative)', so a subject who seroconverts mid-study switches this",
        "column from 0 to 1 at that point and the 1.13 ratio (13% higher CL)",
        "applies only to the ADA-positive records. Both ADA-negative",
        "conclusive and ADA-negative inconclusive were pooled as negative,",
        "and all PK observations with a MISSING ADA status were assumed",
        "negative, so a downstream user reconstructing the analysis dataset",
        "should impute missing to 0 rather than dropping the record. PK",
        "analysis set: 1225/21,929 observations (5.6%) ADA-positive (Table",
        "S2). ADA titer at the time of the sample was screened separately",
        "and not retained; see covariatesDataExcluded."
      ),
      source_name        = "ADA"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese race/ethnicity indicator; 1 = Japanese, 0 = non-Japanese.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = paste(
        "Japanese heritage is broken out on its own, NOT pooled with other",
        "Asian groups: supplement Text S1 $PK derives the indicator as",
        "'RACE=0; IF (RACEN.EQ.3.1) RACE=1', and Figure S3 gives the RACEN",
        "code list in which 3.1 = Japanese while 3.2 = Chinese, 3.3 =",
        "Korean and 3.8 = Other Asian are each distinct levels that all map",
        "to 0 here. RACE_JAPANESE is therefore the correct canonical rather",
        "than RACE_ASIAN or RACE_NEAS (the North East Asian composite, which",
        "would wrongly sweep in the Chinese and Korean subjects). Applied as",
        "a ratio raised to the indicator on V1 (0.920) and on V2 (0.671).",
        "The V1 effect is the ONE structural difference between this model",
        "and the previously published lecanemab PK model, and adding it",
        "dropped the objective function from 173,949.983 to 169,260.786. PK",
        "analysis set: 138 Japanese subjects (8.5%) (Table S2)."
      ),
      source_name        = "RACEN == 3.1 (raw); RACE (derived indicator in $PK); JPN (printed equations)"
    ),
    FORM_LEC_PROCESSB = list(
      description        = "Lecanemab manufacturing Process B drug-product indicator; 1 = Process B, 0 = Process A.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Process A, the earlier drug product; F fixed at 1)",
      notes              = paste(
        "Per-DOSE-RECORD indicator, not a per-subject flag: subjects in the",
        "Study 201 open-label extension switched from Process A to Process B",
        "part-way through (United States, Canada and South Korea from May",
        "2020; Japan from July 2020), so a single subject can contribute",
        "records under both processes. All Study 201 Core and all Study 101",
        "/ 104 lecanemab is Process A; all Study 301 Core lecanemab is",
        "Process B. Relative bioavailability is 0.904 for Process B, i.e.",
        "9.6% lower exposure. NOTE the asymmetry in the between-subject",
        "variability: supplement Text S1 $PK codes 'F1=1; IF (FORM.EQ.1)",
        "F1=THETA(5)*EXP(ETA(4))', so ETA(4) enters ONLY on Process B",
        "records and Process A bioavailability is exactly 1 with no",
        "variability -- the model() block reproduces that branch literally",
        "rather than putting the eta on a shared anchor. Eta shrinkage on",
        "this term is high (60.2%). PK analysis set: 8595 observations",
        "(39.2%) Process A, 13,334 (60.8%) Process B (Table S2). Set to 0",
        "to simulate the earlier drug product and to 1 for the commercial",
        "material used in Clarity AD."
      ),
      source_name        = "FORM"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units       = "year",
      type        = "continuous",
      notes       = paste(
        "Screened on CL, V1 and V2 and retained on none of them. This is a",
        "headline published NULL result rather than a reporting gap: the",
        "Abstract states 'Importantly, age, a well-recognized risk factor",
        "for AD, was not found to significantly affect lecanemab PK', and",
        "the Discussion repeats it. No point estimate exists, so the correct",
        "encoding is the omission of the term rather than a fixed(0)",
        "coefficient. PK analysis set median 72 years, range 50-93 (Table",
        "S2) -- a narrow, uniformly elderly cohort, which is part of why the",
        "null is unsurprising."
      )
    ),
    ADA_TITER = list(
      description = "Anti-drug antibody titer at the time of the PK sample.",
      units       = "(titer)",
      type        = "continuous",
      notes       = paste(
        "Screened on CL as a graded alternative to the binary ADA_POS",
        "indicator and not retained; the binary sample-level status is what",
        "the final model carries. PK analysis set median titer 16, range",
        "1-50,000 (Table S2). Carried as column NWTIT in the supplement Text",
        "S1 $INPUT list."
      )
    ),
    DOSE = list(
      description = "Administered lecanemab dose level.",
      units       = "mg/kg",
      type        = "continuous",
      notes       = paste(
        "Screened on CL as a test for dose-nonlinearity across the 0.3-15",
        "mg/kg range and not retained, supporting the linear two-compartment",
        "structure. The Results note a possible nonlinearity at the very low",
        "(<1 mg/kg) Study 101 doses as one candidate explanation for the",
        "slight CWRES trend in that study, but no nonlinear term entered the",
        "final model."
      )
    ),
    RACE_ASIAN = list(
      description = "Non-Japanese Asian race indicator (Chinese, Korean, or Other Asian).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race/ethnicity was screened as a whole on CL, V1 and V2, and only",
        "the Japanese level survived backward elimination. The Chinese (n =",
        "6, 0.4%), Korean (n = 54, 3.3%) and Other Asian (n = 21, 1.3%)",
        "levels are therefore explicitly NOT folded into the retained",
        "RACE_JAPANESE indicator; see that entry. Listed here so a",
        "downstream user does not mistake the absence of a general Asian",
        "term for a transcription gap."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1619L,
    n_studies      = 4L,
    n_observations = "21,929 serum lecanemab concentrations retained: 653 (3.0%) from Study 101, 395 (1.8%) from Study 104, 7991 (36.4%) from Study 201 Core + OLE, 12,890 (58.8%) from Study 301 Core + OLE. 614 samples were excluded (337 BLQ or BLQ with time-after-dose over 2000 h, 61 missing sampling time, 107 with CWRES > 5, 46 above 600 ug/mL as within-subject outliers, 30 with time-after-dose over 2000 h)",
    age_range      = "median 72 years, range 50-93 (Table S2)",
    weight_range   = "median 72 kg, range 37.7-130.5 (Table S2)",
    sex_female_pct = 49.4,
    race_ethnicity = c(
      White = 80.7, Japanese = 8.5, Korean = 3.3,
      `Black/African American` = 3.0,
      `Asian excluding Chinese/Japanese/Korean` = 1.3,
      Chinese = 0.4,
      `American Indian/Alaskan/Other/Missing` = 2.8
    ),
    disease_state  = "early Alzheimer's disease -- mild cognitive impairment due to AD or mild AD dementia, with confirmed amyloid-beta pathology (studies 201 and 301); studies 101 and 104 also contributed subjects with mild-to-moderate AD",
    dose_range     = "0.3, 1, 3, 10 and 15 mg/kg as single intravenous doses; 2.5, 5 and 10 mg/kg bi-weekly; 0.3, 1, 3, 5 and 10 mg/kg monthly. Every infusion ran over 60 +/- 10 minutes. The approved regimen is 10 mg/kg bi-weekly (1113 of 1619 subjects)",
    regions        = "multicentre international; United States, Canada, Japan, South Korea and Europe across studies 101, 104, 201 and 301 (Clarity AD, NCT03887455)",
    albumin_range  = "median 43 g/L, range 35-54 (Table S2)",
    ada_status     = "sample level: 20,703 observations (94.4%) ADA-negative, 1225 (5.6%) ADA-positive; median ADA titer 16, range 1-50,000 (Table S2)",
    drug_product   = "8595 observations (39.2%) Process A, 13,334 (60.8%) Process B (Table S2)",
    notes          = paste(
      "The base structural model was inherited from the previously",
      "published lecanemab population PK analysis and re-estimated on this",
      "larger pooled dataset; every parameter in ini() is a FINAL estimate",
      "from Table 1 of this paper, not an inherited or initial value. The",
      "supplement Text S1 $THETA / $OMEGA / $SIGMA blocks contain INITIAL",
      "estimates only and were deliberately not used as a source. Model",
      "evaluation used goodness-of-fit plots, prediction-corrected VPCs",
      "stratified by study, and non-parametric bootstrap; bootstrap medians",
      "in Table 1 agree with the point estimates throughout. Precision was",
      "high (%RSE <= 4.23% for core parameters, <= 11.8% for covariate",
      "effects). Eta shrinkage: CL 6.02%, V1 29.3%, V2 25.8%, F 60.2%."
    )
  )

  ini({
    # ==================================================================
    # All values are FINAL estimates from Majid 2024 Table 1 ("Population
    # PK parameters and bootstrap CIs for the final lecanemab covariate
    # model"), cross-checked against the four covariate equations printed
    # immediately below that table and against supplement Text S1 $PK.
    #
    # The printed equations are:
    #   CL = 0.0154 * (BW/72)^0.353 * (ALB/43)^-0.374 * 0.791^SEX * 1.13^ADA
    #   V1 = 3.24   * (BW/72)^0.513 * 0.868^SEX * 0.920^JPN
    #   V2 = 2.00   * 0.671^JPN
    #   F  = 1      * 0.904^FORM
    # with SEX 0/1 = male/female, ADA 0/1 = negative/positive,
    # JPN 0/1 = non-Japanese/Japanese, FORM 0/1 = Process A/Process B.
    #
    # CLOSED-FORM CHECK on the structural block: the two-compartment
    # terminal half-life implied by these four values is
    #   k10 = 0.0154/3.24, k12 = 0.00718/3.24, k21 = 0.00718/2.00,
    #   beta = 0.5*(k10+k12+k21 - sqrt((k10+k12+k21)^2 - 4*k21*k10))
    # giving log(2)/beta = 348 h = 14.5 days, reproducing the Discussion's
    # stated "approximately 14.5 days" exactly. Supplement Text S1 $PK
    # computes T12 by this same formula.
    # ==================================================================

    # ----- Structural parameters (Table 1, "PK parameters" block) -----
    lcl <- log(0.0154) ; label("Clearance for the reference subject (CL, L/h)")                       # Majid 2024 Table 1: CL 0.0154 L/h, %RSE 1.60, bootstrap median 0.0154 (95% CI 0.0147-0.0160)
    lvc <- log(3.24)   ; label("Central volume of distribution for the reference subject (V1, L)")    # Majid 2024 Table 1: V1 3.24 L, %RSE 0.799, bootstrap median 3.24 (95% CI 3.18-3.30)
    lvp <- log(2.00)   ; label("Peripheral volume of distribution for the reference subject (V2, L)") # Majid 2024 Table 1: V2 2.00 L, %RSE 4.09, bootstrap median 2.02 (95% CI 1.83-2.21)
    lq  <- log(0.00718); label("Intercompartmental clearance (Q, L/h)")                               # Majid 2024 Table 1: Q 0.00718 L/h, %RSE 4.23, bootstrap median 0.00701 (95% CI 0.00155-0.0125); no covariate and no IIV (Text S1 $PK: "TVQ=THETA(4); Q=TVQ")

    # Process A bioavailability is a STRUCTURAL ANCHOR, not an estimate:
    # supplement Text S1 $PK opens the branch with a literal "F1=1" and no
    # THETA, so it is fixed here rather than estimated.
    lfcentral <- fixed(log(1)); label("Relative bioavailability of the intravenous dose for manufacturing Process A (F, unitless)")  # Majid 2024 supplement Text S1 $PK: "F1=1" (reference process; not an estimated parameter)

    # ----- Covariate effects on CL (Table 1, "Covariate effects" block) -----
    e_wt_cl     <-  0.353; label("Power exponent on (WT/72 kg) for CL (unitless)")                            # Majid 2024 Table 1 "Weight ~ CL (exponent)" 0.353, %RSE 10.5, bootstrap median 0.344 (95% CI 0.250-0.460)
    e_alb_cl    <- -0.374; label("Power exponent on (ALB/43 g/L) for CL; NEGATIVE, so CL falls as albumin rises (unitless)")  # Majid 2024 Table 1 "Albumin ~ CL (exponent)" -0.374, %RSE 9.71, bootstrap median -0.372 (95% CI -0.482 to -0.261)
    e_female_cl <-  0.791; label("Multiplicative CL ratio for females vs the male reference, applied as ratio^SEXF (unitless)")         # Majid 2024 Table 1 "Females ~ CL (ratio to males)" 0.791, %RSE 2.17, bootstrap median 0.791 (95% CI 0.758-0.824); Discussion: "20.9% [lower] for CL"
    e_ada_cl    <-  1.13 ; label("Multiplicative CL ratio for ADA-positive vs ADA-negative samples, applied as ratio^ADA_POS (unitless)")  # Majid 2024 Table 1 "ADApositive ~ CL (ratio to ADAnegative)" 1.13, %RSE 0.860, bootstrap median 1.13 (95% CI 1.09-1.17)

    # ----- Covariate effects on V1 -----
    e_wt_vc       <- 0.513; label("Power exponent on (WT/72 kg) for V1 (unitless)")                             # Majid 2024 Table 1 "Weight ~ V1 (exponent)" 0.513, %RSE 5.01, bootstrap median 0.514 (95% CI 0.469-0.558)
    e_female_vc   <- 0.868; label("Multiplicative V1 ratio for females vs the male reference, applied as ratio^SEXF (unitless)")  # Majid 2024 Table 1 "Females ~ V1 (ratio to males)" 0.868, %RSE 1.04, bootstrap median 0.868 (95% CI 0.853-0.884); Discussion: "13.2% [lower] for V1"
    e_japanese_vc <- 0.920; label("Multiplicative V1 ratio for Japanese vs the non-Japanese reference, applied as ratio^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1 "Japanese ethnicity ~ V1 (ratio to non-Japanese)" 0.920, %RSE 1.58, bootstrap median 0.920 (95% CI 0.896-0.945); the one term new relative to the previously published model

    # ----- Covariate effect on V2 -----
    e_japanese_vp <- 0.671; label("Multiplicative V2 ratio for Japanese vs the non-Japanese reference, applied as ratio^RACE_JAPANESE (unitless)")  # Majid 2024 Table 1 "Japanese ethnicity ~ V2 (ratio to non-Japanese)" 0.671, %RSE 11.8, bootstrap median 0.665 (95% CI 0.475-0.835)

    # ----- Manufacturing-process comparability on the intravenous dose -----
    e_processb_f <- 0.904; label("Relative bioavailability ratio for Process B vs the Process A reference (unitless)")  # Majid 2024 Table 1 "F (comparability) for process B" 0.904, %RSE 0.750, bootstrap median 0.904 (95% CI 0.890-0.918); Discussion: 9.6% lower exposure

    # ==================================================================
    # Inter-individual variability. Table 1 reports these as "%CV", and
    # the table footnote DEFINES that column as "CV%, square root of
    # variance x 100" -- i.e. the printed numbers are omega SDs on the
    # log scale multiplied by 100, NOT log-normal coefficients of
    # variation. The usual omega^2 = log(CV^2 + 1) back-transform is
    # therefore WRONG here and is deliberately not applied; each variance
    # below is simply (printed %CV / 100)^2.
    #
    # CL and V1 share an OMEGA BLOCK(2) (supplement Text S1
    # "$OMEGA BLOCK(2)"); the off-diagonal is reconstructed from the
    # printed correlation R = 0.144 as R * omega_CL * omega_V1.
    # ==================================================================
    etalcl + etalvc ~ c(0.121801,
                        0.006131, 0.014884)   # Majid 2024 Table 1: CL 34.9% -> 0.349^2 = 0.121801; correlation CL~V1 R = 0.144 -> 0.144*0.349*0.122 = 0.006131; V1 12.2% -> 0.122^2 = 0.014884
    etalvp ~ 0.894916                         # Majid 2024 Table 1: V2 94.6% -> 0.946^2 = 0.894916 (eta shrinkage 25.8%)
    etalfcentral ~ 0.007242                   # Majid 2024 Table 1: F 8.51% -> 0.0851^2 = 0.007242 (eta shrinkage 60.2%); NOTE this eta acts on Process B records ONLY -- see model()

    # ----- Residual variability (Table 1, "Residual variability" block) -----
    propSd <- 0.210; label("Proportional residual error (fraction)")   # Majid 2024 Table 1 "Proportional (%CV)" 21.0, %RSE 1.23, bootstrap median 21.0 (95% CI 20.4-21.5)
    addSd  <- 1.12 ; label("Additive residual error (ug/mL)")          # Majid 2024 Table 1 "Additive (SD; ug/mL)" 1.12, %RSE 16.8, bootstrap median 1.12 (95% CI <0-1.62)
  })

  model({
    # ----- Individual parameters ------------------------------------
    # Each line reproduces one printed covariate equation from the block
    # below Majid 2024 Table 1, with the binary covariates entering as
    # ratio^indicator exactly as printed (identical to the supplement
    # Text S1 $PK "THETA(n)**COV" forms).
    cl <- exp(lcl + etalcl) * (WT / 72)^e_wt_cl * (ALB / 43)^e_alb_cl *
      e_female_cl^SEXF * e_ada_cl^ADA_POS
    vc <- exp(lvc + etalvc) * (WT / 72)^e_wt_vc *
      e_female_vc^SEXF * e_japanese_vc^RACE_JAPANESE
    vp <- exp(lvp + etalvp) * e_japanese_vp^RACE_JAPANESE
    q  <- exp(lq)

    # ----- Micro-constants ------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- Disposition ----------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # ----- Manufacturing-process comparability ----------------------
    # Text S1 $PK is an IF branch, not a shared multiplier:
    #     F1 = 1
    #     IF (FORM.EQ.1) F1 = THETA(5)*EXP(ETA(4))
    # so Process A bioavailability is exactly 1 with NO between-subject
    # variability, and both the 0.904 ratio and its eta apply to Process
    # B records only. The explicit two-branch form below reproduces that;
    # it is algebraically identical to the paper's printed compact form
    # F = 1 * 0.904^FORM once the eta is carried along. The Process B
    # branch is assigned on its own simple line so etalfcentral stays
    # mu-referenced for estimation.
    fprocessb  <- exp(lfcentral + etalfcentral) * e_processb_f
    f(central) <- (1 - FORM_LEC_PROCESSB) + FORM_LEC_PROCESSB * fprocessb

    # ----- Observation ----------------------------------------------
    # Dose in mg over a volume in L gives mg/L = ug/mL, matching the
    # IP/LC-MS/MS assay range of 0.5-150 ug/mL and the additive residual
    # SD of 1.12 ug/mL.
    Cc <- central / vc
    # Text S1 $ERROR is "W = F + 0.01; Y = W + W*ERR(1) + ERR(2)": a
    # combined proportional-plus-additive error in which the proportional
    # term multiplies (prediction + 0.01 ug/mL). The 0.01 ug/mL offset is
    # a numerical guard sitting 50-fold below the 0.5 ug/mL assay LLOQ
    # and is not reproduced here; see the vignette Assumptions section.
    Cc ~ prop(propSd) + add(addSd)
  })
}
