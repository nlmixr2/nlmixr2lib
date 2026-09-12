Yang_2025_osimertinib <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for osimertinib and its",
    "active metabolite AZ5104 in 2,196 patients with EGFR-mutated advanced",
    "non-small cell lung cancer (NSCLC), pooled across six studies (AURA,",
    "AURA2, AURA3, FLAURA, ADAURA, and FLAURA2) (Yang 2025). This is the",
    "FLAURA2 update of the earlier Brown 2017 model: it retains the same",
    "published structure (first-order oral absorption into a one-compartment",
    "parent disposition, feeding a one-compartment metabolite in series, with",
    "first-order elimination from both) and re-estimates it on the pooled",
    "dataset that adds the FLAURA2 osimertinib-plus-chemotherapy arm. The",
    "fraction of parent clearance appearing as AZ5104 is fixed at 0.25.",
    "Retained covariates are baseline body weight and baseline serum albumin",
    "on parent CL/F and V/F, baseline body weight, baseline serum albumin and",
    "race (Japanese and Asian-other, each vs a White reference) on AZ5104",
    "CL/F, and baseline serum albumin on AZ5104 V/F. Chinese and Other race",
    "indicators were carried in the final control stream but their",
    "coefficients were fixed to zero after backward elimination found them",
    "non-significant. Adding chemotherapy was tested on both parent and",
    "metabolite clearance and was not significant, i.e. the model found no",
    "pemetrexed-platinum drug interaction. No covariate had a clinically",
    "meaningful effect. Concentrations are in nM, matching the source.",
    sep = " "
  )
  reference <- paste(
    "Yang J, Olabode D, Sawant-Basak A, Baldry R, Vishwanathan K, Bachina S,",
    "Todd A, Ghiorghiu D, Rukazenkov Y, Zhou D, Shahraz A. Population",
    "Pharmacokinetics and Exposure-Response Analysis of First-Line",
    "Osimertinib Plus Chemotherapy in Patients with EGFR-Mutated Advanced",
    "NSCLC. Clin Pharmacol Ther. 2025;118(5):1110-1120.",
    "doi:10.1002/cpt.3759",
    sep = " "
  )
  vignette <- "Yang_2025_osimertinib"
  units <- list(time = "h", dosing = "mg", concentration = "nM")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are carried in mg of parent-mass equivalent
  # throughout, exactly as in the source NONMEM control stream: the
  # central -> metabolite transfer is 1:1 in amount (no molar
  # stoichiometric correction), and each analyte's concentration is then
  # formed with its OWN molecular weight via the control stream's scaling
  # factors S2 = V1 * 499.61 / 1e6 and S3 = VM1 * 485.59 / 1e6. Verified
  # against Table S5: this convention reproduces the observed AZ5104
  # steady-state AUC to within 1%, whereas a molar 1:1 transfer is low by
  # about 3%.
  compartmentData <- list(
    depot          = list(analyte = "osimertinib", units = "mg", specimen = "administration site", verified = TRUE),
    central        = list(analyte = "osimertinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_az5104 = list(analyte = "AZ5104", units = "mg (parent-mass equivalent)", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_value    = "62 kg (Yang 2025 Table 1 overall median; the reference used in Eqs 2-4 and in the Figure 3 forest-plot reference patient).",
      notes              = "Power (allometric-form) effect on parent CL/F (exponent 0.36), parent V/F (0.64) and AZ5104 CL/F (0.74). The source control stream maps a missing weight (coded -999) onto a covariate factor of 1, i.e. the typical value; this implementation has no missingness code and expects an observed weight.",
      source_name        = "WT (Yang 2025 supplementary NONMEM control stream $INPUT)"
    ),
    ALB = list(
      description        = "Serum albumin at baseline.",
      units              = "g/L",
      type               = "continuous",
      reference_value    = "40 g/L (Yang 2025 Table 1 overall median; the reference used in Eqs 2-5).",
      notes              = "Power effect on parent CL/F (exponent 0.67), parent V/F (1.57), AZ5104 CL/F (0.72) and AZ5104 V/F (-0.65). The AZ5104 V/F exponent is NEGATIVE: Table 2 prints its magnitude (0.65) only, but Eq 5 and $THETA(13) = -0.652465 both carry the sign. Reported in g/L (SI); US-convention g/dL values must be multiplied by 10.",
      source_name        = "BALB (Yang 2025 supplementary NONMEM control stream $INPUT)"
    ),
    RACE_ASIAN_OTH = list(
      description        = "Asian-other race indicator (1 = Asian heritage other than Chinese or Japanese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper-defined race reference category is White).",
      notes              = "Linear additive effect (1 + 0.161274 * RACE_ASIAN_OTH) on AZ5104 CL/F, i.e. 16.1% higher metabolite clearance than a White patient. Corresponds to ETHL = 1 in the source control stream. Yang 2025 carries White (reference), Asian-other, Chinese, Japanese and Other as five mutually exclusive levels; the dominant reference cohort is White, not Chinese.",
      source_name        = "ETHL = 1 (Yang 2025 supplementary NONMEM control stream); 'Asian (excluding Chinese and Japanese)' (Table 1); 'Asian (NonCHN or nonJPN) on CLm/F' (Table 2)"
    ),
    RACE_JAPANESE = list(
      description        = "Japanese-heritage race indicator (1 = Japanese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper-defined race reference category is White).",
      notes              = "Linear additive effect (1 + 0.181966 * RACE_JAPANESE) on AZ5104 CL/F, i.e. 18.2% higher metabolite clearance than a White patient. Corresponds to ETHL = 3 in the source control stream.",
      source_name        = "ETHL = 3 (Yang 2025 supplementary NONMEM control stream); 'JPN on CLm/F' (Table 2)"
    ),
    RACE_CHINESE = list(
      description        = "Chinese-heritage race indicator (1 = Chinese, 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper-defined race reference category is White).",
      notes              = "Carried in the final model structure but with its coefficient FIXED TO ZERO, so it has no effect on AZ5104 CL/F. Backward elimination (Table S2) found that removing Chinese from CLm/F was 'Not significant', and the source control stream accordingly holds $THETA(15) at '0 FIX' for ETHL = 2. Retained here so the covariate screen is auditable rather than silently dropped.",
      source_name        = "ETHL = 2 (Yang 2025 supplementary NONMEM control stream $THETA(15), 0 FIX)"
    ),
    RACE_OTHER = list(
      description        = "Race-category 'Other' indicator (1 = Hispanic/Latino, Native American, Native Alaskan/Inuit, Native Hawaiian/Pacific Islander, African, African-American, African-Caribbean, or missing; 0 = otherwise).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the paper-defined race reference category is White).",
      notes              = "Carried in the final model structure but with its coefficient FIXED TO ZERO, so it has no effect on AZ5104 CL/F. Backward elimination (Table S2) found that removing other ethnicity from CLm/F was 'Not significant', and the source control stream accordingly holds $THETA(17) at '0 FIX' for ETHL = 99. The composite membership is defined in the Yang 2025 Table 1 footnote a.",
      source_name        = "ETHL = 99 (Yang 2025 supplementary NONMEM control stream $THETA(17), 0 FIX)"
    )
  )

  # Covariates that Yang 2025 evaluated but did not retain. Recorded so the
  # covariate screen stays auditable; none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.", units = "years", type = "continuous",
      notes = "Median 62.0 (range 25.0-91.0) overall (Table 1). Not retained in the final PopPK model."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline.", units = "mL/min", type = "continuous",
      notes = "Forward inclusion on parent CL/F was 'Not significant' (Table S2)."
    ),
    RENAL_IMPAIRMENT = list(
      description = "Grouped renal impairment status (normal / mild / moderate-severe).", units = "(category)", type = "categorical",
      notes = "Forward inclusion on parent CL/F and AZ5104 CL/F was 'Not significant' (Table S2); the paper concludes no renal dose adjustment is needed."
    ),
    HEPATIC_IMPAIRMENT = list(
      description = "Grouped hepatic impairment status (normal / at least mild).", units = "(category)", type = "categorical",
      notes = "Forward inclusion on parent CL/F and AZ5104 CL/F was 'Not significant' (Table S2); the paper concludes no hepatic dose adjustment is needed."
    ),
    CHEMO_COMBINATION = list(
      description = "Co-administration of pemetrexed plus platinum chemotherapy.", units = "(binary)", type = "binary",
      notes = "The pre-specified covariate of interest for this analysis. Forward inclusion on both parent CL/F and AZ5104 CL/F was 'Not significant' (Table S2), i.e. no osimertinib-chemotherapy drug interaction; this is the paper's central PK finding."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2196,
    n_studies      = 6,
    n_samples      = 52338,
    studies        = "AURA (including its extension phase), AURA2, AURA3, FLAURA, ADAURA, FLAURA2 (NCT04035486)",
    disease_state  = "EGFR-mutated locally advanced or metastatic non-small cell lung cancer (NSCLC)",
    age_range      = "25.0-91.0 years (overall median 62.0; mean 61.5, SD 10.7)",
    weight_range   = "29.0-128 kg (overall median 62.0; mean 63.5, SD 14.1)",
    albumin_range  = "17.0-53.3 g/L (overall median 40.0; mean 39.7, SD 4.92)",
    crcl_range     = "20.5-196 mL/min (overall median 82.1; mean 85.1, SD 26.8)",
    bmi_range      = "12.9-42.6 kg/m2 (overall median 23.4)",
    sex_female_pct = 64.5,
    race_ethnicity = list(
      White                                  = 26.6,
      `Asian (excluding Chinese and Japanese)` = 22.9,
      Chinese                                = 22.8,
      Japanese                               = 17.5,
      Other                                  = 10.2
    ),
    line_of_therapy = list(`First-line` = 38.5, `Second-line` = 21.3, `Third-line and later` = 25.4, Adjuvant = 14.8),
    who_ps          = list(`0` = 40.8, `1` = 59.2),
    dose_range      = "Osimertinib 80 mg once daily (reduction to 40 mg once daily permitted for toxicity). In the FLAURA2 combination arm, given with pemetrexed 500 mg/m2 plus either cisplatin 75 mg/m2 or carboplatin AUC 5 mg/mL/min on Day 1 of 21-day cycles for 4 cycles, then pemetrexed 500 mg/m2 maintenance every 3 weeks.",
    notes           = "Demographics reproduced from Yang 2025 Table 1 (Overall, N = 2,196). Of 52,338 plasma samples, 3,636 (7.0%) were removed and 3,250 (6.2%) below-quantification-limit samples were handled by the M1 method. Race percentages are of the overall pooled population; the reference (typical) patient used for the Figure 3 forest plot is a White patient of 62 kg with albumin 40 g/L."
  )

  ini({
    # ---- Structural parameters -------------------------------------------
    # Typical values for the reference patient: White, 62 kg, albumin 40 g/L.
    # Final estimates are taken at full precision from the $THETA block of the
    # supplementary NONMEM control stream ($PROBLEM FINAL MODEL); each rounds
    # to the value printed in Table 2 and in Eqs 1-5.
    lka           <- log(0.260946);  label("First-order oral absorption rate constant, Ka (1/h)")                       # Yang 2025 $THETA(3) = 0.260946; Table 2 Ka = 0.26; Eq 1
    lcl           <- log(14.3452);   label("Apparent total parent clearance, CLptot/F, reference patient (L/h)")        # Yang 2025 $THETA(1) = 14.3452; Table 2 CLptot/F = 14.35; Eq 2
    lvc           <- log(1151.11);   label("Apparent parent volume of distribution, Vp/F, reference patient (L)")       # Yang 2025 $THETA(4) = 1151.11; Table 2 Vp/F = 1,151; Eq 3
    lcl_az5104    <- log(32.0484);   label("Apparent AZ5104 clearance, CLm/F, reference patient (L/h)")                 # Yang 2025 $THETA(2) = 32.0484; Table 2 CLm/F = 32.05; Eq 4
    lvc_az5104    <- log(151.014);   label("Apparent AZ5104 volume of distribution, Vm/F, reference patient (L)")       # Yang 2025 $THETA(6) = 151.014; Table 2 Vm/F = 151; Eq 5

    # Fraction of parent clearance appearing as AZ5104. Fixed, not estimated:
    # fm, Vm/F and CLm/F are mutually confounded in a joint parent-metabolite
    # model fitted to plasma data alone.
    fm            <- fixed(0.25);    label("Fraction of parent clearance forming AZ5104 (unitless)")                    # Yang 2025 $THETA(5) = 0.25 FIX; Table 2 Fm = 0.25 (Fixed); Methods "Fm was fixed at 25%"

    # ---- Continuous-covariate power effects ------------------------------
    # Form: (WT / 62)^exponent and (ALB / 40)^exponent, per Eqs 2-5.
    e_wt_cl            <- 0.359774;  label("Power exponent for body weight on parent CL/F (unitless)")                  # Yang 2025 $THETA(7) = 0.359774; Table 2 WT on CLptot/F = 0.36; Eq 2
    e_alb_cl           <- 0.674479;  label("Power exponent for albumin on parent CL/F (unitless)")                      # Yang 2025 $THETA(8) = 0.674479; Table 2 ALB on CLptot/F = 0.67; Eq 2
    e_wt_vc            <- 0.640227;  label("Power exponent for body weight on parent V/F (unitless)")                   # Yang 2025 $THETA(9) = 0.640227; Table 2 WT on Vp/F = 0.64; Eq 3
    e_alb_vc           <- 1.56611;   label("Power exponent for albumin on parent V/F (unitless)")                       # Yang 2025 $THETA(10) = 1.56611; Table 2 ALB on Vp/F = 1.57; Eq 3
    e_wt_cl_az5104     <- 0.740376;  label("Power exponent for body weight on AZ5104 CL/F (unitless)")                  # Yang 2025 $THETA(11) = 0.740376; Table 2 WT on CLm/F = 0.74; Eq 4
    e_alb_cl_az5104    <- 0.722460;  label("Power exponent for albumin on AZ5104 CL/F (unitless)")                      # Yang 2025 $THETA(12) = 0.722460; Table 2 ALB on CLm/F = 0.72; Eq 4
    # NEGATIVE exponent: Table 2 prints only the magnitude (0.65), but Eq 5
    # writes (ALB/40)^-0.65 and $THETA(13) = -0.652465 carries the sign.
    e_alb_vc_az5104    <- -0.652465; label("Power exponent for albumin on AZ5104 V/F (unitless)")                       # Yang 2025 $THETA(13) = -0.652465; Eq 5 exponent -0.65; Table 2 ALB on Vm/F prints 0.65 unsigned

    # ---- Categorical-covariate linear effects ----------------------------
    # Form: (1 + coefficient * indicator) on AZ5104 CL/F, per Eq 4 and the
    # control stream's CLM1ETHL block. White (ETHL = 0) is the reference.
    e_race_asian_oth_cl_az5104 <- 0.161274;  label("Linear coefficient for Asian-other vs White on AZ5104 CL/F (unitless)")  # Yang 2025 $THETA(14) = 0.161274; Table 2 = 0.16; Eq 4 factor 1.16
    e_race_japanese_cl_az5104  <- 0.181966;  label("Linear coefficient for Japanese vs White on AZ5104 CL/F (unitless)")     # Yang 2025 $THETA(16) = 0.181966; Table 2 = 0.18; Eq 4 factor 1.18
    # Fixed to zero after backward elimination judged them non-significant
    # (Table S2); kept so the covariate screen stays visible.
    e_race_chinese_cl_az5104   <- fixed(0);  label("Linear coefficient for Chinese vs White on AZ5104 CL/F (unitless; no retained effect)")  # Yang 2025 $THETA(15) = 0 FIX; Table S2 "Removing Chinese from CLm/F: Not significant"
    e_race_other_cl_az5104     <- fixed(0);  label("Linear coefficient for Other vs White on AZ5104 CL/F (unitless; no retained effect)")    # Yang 2025 $THETA(17) = 0 FIX; Table S2 "Removing other ethnicity from CLm/F: Not significant"

    # ---- Between-subject variability -------------------------------------
    # Table 2's "Between subject variability" column is on the VARIANCE scale.
    # Confirmed three ways: (i) the values equal the control stream's $OMEGA
    # entries exactly; (ii) back-transforming each as a variance via
    # sqrt(exp(omega2) - 1) reproduces every published %CV to within 1
    # percentage point (46/52/101/85/82% for CLptot/F, CLm/F, Ka, Vp/F, Vm/F),
    # whereas reading them as SDs gives 19/24/80/58/55%; (iii) the off-diagonal
    # is a $OMEGA BLOCK(2) covariance, not a correlation coefficient -- it
    # implies a CLptot/F ~ CLm/F correlation of 0.88, which Table 2's
    # "Correlation" row reports unconverted.
    etalcl + etalcl_az5104 ~ c(0.191276,
                               0.186815, 0.237930)                                # Yang 2025 $OMEGA BLOCK(2); Table 2 IIV on CLptot/F = 0.19, Correlation = 0.19, IIV on CLm/F = 0.24
    etalka        ~ 0.699828                                                      # Yang 2025 $OMEGA 3; Table 2 IIV on Ka = 0.70 (%CV 101)
    etalvc        ~ 0.542208                                                      # Yang 2025 $OMEGA 4; Table 2 IIV on Vp/F = 0.54 (%CV 85)
    etalvc_az5104 ~ 0.512184                                                      # Yang 2025 $OMEGA 5; Table 2 IIV on Vm/F = 0.51 (%CV 82)

    # ---- Residual unexplained variability ---------------------------------
    # Combined additive plus proportional, estimated separately per analyte.
    # The control stream's $ERROR builds W = sqrt(THETA_add^2 + (THETA_prop *
    # IPRED)^2) with $SIGMA 1 FIX, so each THETA is a standard deviation on the
    # nM concentration scale and maps directly onto add()/prop().
    propSd         <- 0.219108;  label("Proportional residual error on osimertinib (fraction)")   # Yang 2025 $THETA(19) = 0.219108; Table 2 parent proportional component = 0.22
    addSd          <- 29.5202;   label("Additive residual error on osimertinib (nM)")             # Yang 2025 $THETA(18) = 29.5202; Table 2 parent additive component = 29.52 nM
    propSd_az5104  <- 0.231740;  label("Proportional residual error on AZ5104 (fraction)")        # Yang 2025 $THETA(21) = 0.231740; Table 2 metabolite proportional component = 0.23
    addSd_az5104   <- 0.358422;  label("Additive residual error on AZ5104 (nM)")                  # Yang 2025 $THETA(20) = 0.358422; Table 2 metabolite additive component = 0.36 nM
  })

  model({
    # Molecular weights used by the source control stream's scaling factors
    # S2 = V1 * 499.61 / 1e6 (osimertinib) and S3 = VM1 * 485.59 / 1e6
    # (AZ5104, N-desmethyl osimertinib). They convert an amount in mg and a
    # volume in L into a concentration in nM, which is the scale on which
    # Yang 2025 reports every concentration, the assay range (16-8,010 nM
    # osimertinib; 1.65-824 nM AZ5104), the additive residual errors and the
    # AUCss quartiles.
    mw_parent <- 499.61
    mw_az5104 <- 485.59

    # Reference covariate values (Yang 2025 Table 1 overall medians; the
    # denominators written into Eqs 2-5 and the control stream).
    ref_wt  <- 62
    ref_alb <- 40

    # ---- Individual parameters (Eqs 1-5) ---------------------------------
    ka <- exp(lka + etalka)

    cl <- exp(lcl + etalcl) *
      (WT / ref_wt)^e_wt_cl *
      (ALB / ref_alb)^e_alb_cl

    vc <- exp(lvc + etalvc) *
      (WT / ref_wt)^e_wt_vc *
      (ALB / ref_alb)^e_alb_vc

    cl_az5104 <- exp(lcl_az5104 + etalcl_az5104) *
      (WT / ref_wt)^e_wt_cl_az5104 *
      (ALB / ref_alb)^e_alb_cl_az5104 *
      (1 + e_race_asian_oth_cl_az5104 * RACE_ASIAN_OTH) *
      (1 + e_race_japanese_cl_az5104  * RACE_JAPANESE) *
      (1 + e_race_chinese_cl_az5104   * RACE_CHINESE) *
      (1 + e_race_other_cl_az5104     * RACE_OTHER)

    vc_az5104 <- exp(lvc_az5104 + etalvc_az5104) *
      (ALB / ref_alb)^e_alb_vc_az5104

    # ---- ODE system ------------------------------------------------------
    # Reproduces the control stream's ADVAN7 rate constants exactly:
    #   K12 = KA, K20 = CL * (1 - FM) / V1, K23 = CL * FM / V1,
    #   K30 = CLM1 / VM1.
    # Total efflux from central is therefore cl/vc (the K20 + K23 sum), of
    # which the fraction fm is routed to the metabolite. The transfer is 1:1
    # in amount with NO molar stoichiometric correction -- that is what the
    # source fitted, and it is what reproduces the published AZ5104 exposure
    # (see compartmentData above).
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot - (cl / vc) * central
    d/dt(central_az5104) <-  fm * (cl / vc) * central -
      (cl_az5104 / vc_az5104) * central_az5104

    # ---- Observations (nM) -----------------------------------------------
    Cc        <- central        / (vc        * mw_parent / 1e6)
    Cc_az5104 <- central_az5104 / (vc_az5104 * mw_az5104 / 1e6)

    Cc        ~ prop(propSd)        + add(addSd)
    Cc_az5104 ~ prop(propSd_az5104) + add(addSd_az5104)
  })
}
