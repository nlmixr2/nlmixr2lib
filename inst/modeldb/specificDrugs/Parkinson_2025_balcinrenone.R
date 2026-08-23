Parkinson_2025_balcinrenone <- function() {
  description <- paste0(
    "Two-compartment population pharmacokinetic model for balcinrenone ",
    "(AZD9977), a selective mineralocorticoid receptor modulator, pooled ",
    "across six clinical studies (184 participants, 1882 observations) in ",
    "healthy participants, participants with renal impairment, and patients ",
    "with heart failure and chronic kidney disease receiving the ",
    "immediate-release capsule. Absorption is sequential zero-order input ",
    "into the depot (D1 = 0.631 h fasted) followed by first-order transfer ",
    "(KA = 0.381 1/h). Apparent clearance (9.39 L/h in patients) carries a ",
    "power effect of baseline CKD-EPI eGFR ((eGFR/57.68)^0.435) and a ",
    "study-type effect (2.06-fold higher CL/F in the single-dose phase 1 ",
    "studies than in the phase 1b/2b heart-failure/CKD patient studies). ",
    "Food state acts on three parameters: relative bioavailability ",
    "(+13.4% fed), zero-order duration (+52.1% fed) and intercompartmental ",
    "clearance (-79.0% fed). The absorption rate constant decreases with ",
    "dose ((dose/150 mg)^-0.262), consistent with solubility-limited ",
    "absorption. Body weight enters through allometric scaling with ",
    "exponents fixed at 0.75 (CL/F, Q/F) and 1 (Vc/F, Vp/F) referenced to ",
    "70 kg. Interindividual variability is exponential on all seven ",
    "structural parameters and residual variability is combined ",
    "proportional (43.7% CV) plus additive (5.05 nmol/L)."
  )
  reference <- paste(
    "Parkinson J, Astrand M, Melin J, Ericsson H. (2025). Population",
    "Pharmacokinetic Analysis of Balcinrenone in Healthy Participants and",
    "Participants with Heart Failure and Chronic Kidney Disease.",
    "Clin Pharmacokinet. doi:10.1007/s40262-025-01572-7. PMCID: PMC12618412.",
    "Structural encoding follows the final NONMEM control stream reproduced",
    "in the Supplementary Materials ('Final NONMEM Model'); all point",
    "estimates are the final estimates of Table 1 (the control stream's",
    "$THETA / $OMEGA / $SIGMA blocks hold initial values only).",
    sep = " "
  )
  vignette <- "Parkinson_2025_balcinrenone"

  # Dose amounts are in nmol because the NONMEM $ERROR block computes
  # CP = A(2)/S2 with S2 = VC (L) against observations in nmol/L (Supplementary
  # Materials 'Summary of Bioanalytical Methods': assay QC concentrations 30,
  # 500 and 8000 nmol/L; Methods 2.2: LLOQ 2 nmol/L for NCT04595370 and
  # 10 nmol/L elsewhere). Milligram doses convert with the balcinrenone molar
  # mass 399.4 g/mol; see covariateData$DOSE_BALCINRENONE_MG.
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  compartmentData <- list(
    depot = list(
      analyte = "balcinrenone", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "balcinrenone", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "balcinrenone", units = "nmol",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Baseline body weight. Allometric scaling referenced to 70 kg with ",
        "exponents FIXED at 0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F ",
        "(supplement 'Final NONMEM Model': ALLOCL = (BWT/70)**0.75, ",
        "ALLOV = (BWT/70)**1). Body weight was NOT selected by the stepwise ",
        "covariate search (supplement Table S4 does not test it); Results ",
        "3.2 states it was added afterwards by allometry with fixed ",
        "exponents and the SCM search then repeated. Cohort body weight ",
        "83.8 (15.1) kg mean (SD), median 83.0 kg, range 47.0-128.2 kg ",
        "(supplement Table S2)."
      ),
      source_name        = "BWT"
    ),
    CRCL = list(
      description        = paste0(
        "Baseline estimated glomerular filtration rate, CKD-EPI equation ",
        "(Levey 2009), BSA-normalised"
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed baseline value (the source column is BEGFR). Enters CL/F ",
        "as the power ratio (CRCL/57.68)^0.435; 57.68 is the pooled-dataset ",
        "median baseline eGFR taken from the supplement's final control ",
        "stream, which Table 1's footnote rounds to 57.7 (supplement ",
        "Table S2 reports the pooled median as 57.7 mL/min/1.73 m^2). ",
        "Cohort eGFR 63.5 (25.0) mean (SD), median 57.7, range 21-127.7 ",
        "mL/min/1.73 m^2 (supplement Table S2). Discussion cautions that ",
        "the power function was only validated over the observed eGFR range ",
        "and may be biased as eGFR approaches 0. Only one renal-function ",
        "column exists in this analysis, so the canonical CRCL is used ",
        "rather than the Wahlby-style CRCL_BASE (which is reserved for ",
        "papers that split a time-varying CRCL from its baseline)."
      ),
      source_name        = "BEGFR"
    ),
    FED = list(
      description        = "1 = dose taken in the fed state, 0 = dose taken fasted",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste0(
        "Dose-record level. The source column is FOOD and is coded with the ",
        "OPPOSITE polarity: the supplement's final control stream sets the ",
        "multiplier to 1 when FOOD = 1 and to (1 + THETA) when FOOD = 0, so ",
        "FED = 1 - FOOD. The published reference state is unchanged by the ",
        "recoding -- fasted remains the reference (all multipliers equal 1), ",
        "which matches the canonical FED reference category and the paper's ",
        "reference participant ('received balcinrenone in a fasted state', ",
        "Methods 2.6). Six independent confirmations that FOOD = 1 is the ",
        "fasted state: Table 1 labels the D1 estimate 0.631 h 'fasted ",
        "state'; the Table 1 footnote states 'Fed state F1 = 1 + F1~Food'; ",
        "Results 3.2 reports longer D1, higher F1 and lower Q/F in the fed ",
        "state, matching the signs +0.521, +0.134 and -0.790; and the ",
        "forest plot's fed-vs-fasted AUCss ratio 1.13 equals 1 + 0.134. ",
        "Studies NCT03804645 and NCT04798222 used a high-fat, high-calorie ",
        "breakfast, but the phase 1b/2b patient records were assumed fed ",
        "without a controlled meal (Discussion), so the general FED ",
        "indicator applies rather than FED_HIGHFAT."
      ),
      source_name        = "FOOD (inverted: FED = 1 - FOOD)"
    ),
    STUDY_BALCINRENONE_PHASE1 = list(
      description        = paste0(
        "1 = participant from one of the four single-dose phase 1 studies ",
        "(NCT03843060, NCT03804645, NCT04469907, NCT04798222; healthy ",
        "participants and participants with renal impairment); 0 = ",
        "participant from the multiple-dose phase 1b/2b studies in patients ",
        "with heart failure and chronic kidney disease (NCT03682497, ",
        "NCT04595370)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (phase 1b/2b patient with heart failure and CKD)",
      notes              = paste0(
        "Subject-level, time-fixed. The source column is PTSFLAG and is ",
        "coded with the OPPOSITE polarity (PTSFLAG = 1 for patients, the ",
        "most common group at 135 of 189 participants), so ",
        "STUDY_BALCINRENONE_PHASE1 = 1 - PTSFLAG. The reference category is ",
        "preserved: the supplement's control stream sets CLPTSFLAG = 1 when ",
        "PTSFLAG = 1, and Methods 2.6 defines the reference participant as ",
        "'a patient with HF and CKD'. Phase 1 participants have 2.06-fold ",
        "higher CL/F (= 1 + 1.06), corresponding to 0.48-fold AUCss ",
        "(Results 3.3). Discussion notes the covariate is empirical -- it ",
        "persists after adjusting for eGFR and body weight and may reflect ",
        "single- vs repeated-dose design, less accurate dose/sample timing ",
        "in the patient studies, or unrecorded food state."
      ),
      source_name        = "PTSFLAG (inverted: STUDY_BALCINRENONE_PHASE1 = 1 - PTSFLAG)"
    ),
    DOSE_BALCINRENONE_MG = list(
      description        = "Balcinrenone dose administered on the current dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Per-dose-record covariate driving the absorption-rate effect ",
        "KA = TVKA * (DOSE/150)^-0.262 (supplement 'Final NONMEM Model' ",
        "KADOSE block; Table 1 footnote 'Dose KA = (dose/150)**KA~Dose'). ",
        "A dedicated mg column is required because the model's amounts are ",
        "in nmol while the published relationship is calibrated in mg; ",
        "convert with the balcinrenone molar mass 399.4 g/mol ",
        "(C20H18FN3O5; PubChem CID 118599727), i.e. ",
        "amt [nmol] = DOSE_BALCINRENONE_MG * 1e6 / 399.4. The molar mass is ",
        "NOT reported in Parkinson 2025 and is not paper-derived. Doses ",
        "contributing to the analysis were 15, 50, 100, 150, 200 and 300 mg ",
        "(supplement Table S1); the reference 150 mg is the renal-impairment ",
        "study dose."
      ),
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 184L,
    n_studies      = 6L,
    age_range      = "20-90 years",
    age_median     = "69 years",
    weight_range   = "47.0-128.2 kg",
    weight_median  = "83.0 kg",
    sex_female_pct = 25.9,
    race_ethnicity = c(White = 79.37, Black = 11.11, Asian = 8.99, Other = 0.53),
    disease_state  = paste0(
      "pooled: healthy participants (NCT03843060, NCT03804645, ",
      "NCT04798222), participants with renal impairment without heart ",
      "failure (NCT04469907), and patients with heart failure (HFmrEF or ",
      "HFpEF) and chronic kidney disease (NCT03682497, NCT04595370)"
    ),
    renal_function = paste0(
      "baseline CKD-EPI eGFR 63.5 (25.0) mL/min/1.73 m^2 mean (SD), median ",
      "57.7, range 21-127.7 (supplement Table S2); the phase 2b cohort ",
      "median was 49.1 and the healthy bioavailability-study cohorts ",
      "99.6-106.6"
    ),
    co_medication  = paste0(
      "dapagliflozin co-administered in 107 of 189 participants (56.6%); ",
      "tested in the SCM and not retained in the final model (supplement ",
      "Table S4)"
    ),
    dose_range     = paste0(
      "15-300 mg oral immediate-release capsule; single doses of 50, 100, ",
      "150 and 300 mg in the phase 1 studies and once-daily 15, 50, 100, ",
      "150 or 200 mg for 12-28 days in the phase 1b/2b studies ",
      "(supplement Table S1)"
    ),
    n_observations = "1882 plasma concentrations used from 2031 available (17% below the limit of quantification, handled by Beal's M3 method in the final model)",
    notes          = paste0(
      "Demographics are the pooled 'Total' column of supplement Table S2, ",
      "which summarises all 189 participants evaluable for PK; 184 ",
      "participants and 1882 observations entered the final analysis after ",
      "exclusion of samples with unreliable dosing/sampling times ",
      "(Results 3.1 and supplement 'Summary of Excluded Pharmacokinetic ",
      "Samples'). Ethnicity was 90.5% not Hispanic or Latino. Estimation ",
      "used NONMEM 7.4 importance sampling; the covariance step provided ",
      "parameter uncertainty and a 1000-sample bootstrap the reported ",
      "bootstrap intervals."
    )
  )

  ini({
    # ========================================================================
    # Structural parameters -- Parkinson 2025 Table 1, 'Structural model
    # parameters'. Typical values apply to the reference participant: a
    # phase 1b/2b patient with heart failure and CKD (STUDY_BALCINRENONE_PHASE1
    # = 0), 70 kg, baseline eGFR 57.68 mL/min/1.73 m^2, dosed fasted
    # (FED = 0) with a 150 mg dose.
    # ========================================================================
    lcl <- log(9.39)
    label("Apparent clearance CL/F (L/h)")
    # Table 1: CL/F = 9.39 L/h (95% CI 8.21-10.6; bootstrap median 9.38,
    # bootstrap RSE 6.98%).

    lvc <- log(23.4)
    label("Apparent central volume Vc/F (L)")
    # Table 1: Vc/F = 23.4 L (95% CI 17.2-29.6; bootstrap median 23.56,
    # bootstrap RSE 14.82%).

    lq <- log(9.53)
    label("Apparent intercompartmental clearance Q/F (L/h)")
    # Table 1: Q/F = 9.53 L/h (95% CI 7.25-11.8; bootstrap median 9.30,
    # bootstrap RSE 17.13%).

    lvp <- log(48.8)
    label("Apparent peripheral volume Vp/F (L)")
    # Table 1: Vp/F = 48.8 L (95% CI 34.1-63.5; bootstrap median 49.03,
    # bootstrap RSE 17.64%).

    lka <- log(0.381)
    label("First-order absorption rate constant KA (1/h) at a 150 mg dose")
    # Table 1: KA = 0.381 1/h (95% CI 0.339-0.423; bootstrap median 0.380,
    # bootstrap RSE 4.80%). Reported at the 150 mg reference dose; see
    # e_dose_balcinrenone_mg_ka for the dose dependence.

    ld1 <- log(0.631)
    label("Zero-order input duration into the depot D1 (h), fasted state")
    # Table 1: D1 = 0.631 h (95% CI 0.404-0.858; bootstrap median 0.643,
    # bootstrap RSE 19.29%), described as 'Zero-order duration, fasted
    # state'. The dose is delivered into `depot` at a constant rate over
    # [0, D1] and then transfers first-order at KA to `central`
    # (Results 3.2; supplement 'Final NONMEM Model').

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F1 (unitless) in the fasted state")
    # Supplement 'Final NONMEM Model' $THETA line 6: `1 FIX ; 6. F1`. F1 is
    # a fixed anchor of 1 in the fasted reference state; only the fed-state
    # deviation e_fed_fdepot was estimated. F1 is a RELATIVE bioavailability
    # -- the absolute oral bioavailability of balcinrenone is 52% per the
    # ADME study (Introduction; Lindmark 2023 doi:10.1124/dmd.122.001240) --
    # so CL/F and the volumes are apparent parameters throughout.

    # ========================================================================
    # Allometric body-weight scaling -- Results 3.2 and supplement 'Final
    # NONMEM Model'. Both exponents were fixed at their canonical values
    # rather than estimated; body weight was not selected by the SCM search.
    # ========================================================================
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on (WT/70) shared by CL/F and Q/F (unitless)")
    # Supplement 'Final NONMEM Model': ALLOCL = (BWT/70)**0.75, applied to
    # both TVCL and TVQ. Results 3.2: 'allometric scaling with the fixed
    # exponents (0.75 for clearance and 1.0 for volume of distribution)'.

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on (WT/70) shared by Vc/F and Vp/F (unitless)")
    # Supplement 'Final NONMEM Model': ALLOV = (BWT/70)**1, applied to both
    # TVVC and TVVP.

    # ========================================================================
    # Covariate effects -- Parkinson 2025 Table 1, 'Covariate effect
    # parameters'. Every food-state and study-type effect enters as the
    # multiplicative form (1 + theta * indicator), following the
    # (1 + THETA(n)) branches of the supplement's control stream.
    # ========================================================================
    e_crcl_cl <- 0.435
    label("Power exponent on (baseline eGFR / 57.68) for CL/F (unitless)")
    # Table 1 CL/F~eGFR = 0.435 (95% CI 0.232-0.638; bootstrap median 0.433,
    # bootstrap RSE 29.00%). Table 1 footnote:
    # CL/F = TVCL/F * ((BEGFR/57.7)**CL/F~eGFR). The centring constant is
    # 57.68 in the supplement's final control stream (CLBEGFR block), which
    # Table 1 rounds to 57.7; the difference changes CL/F by 0.02%.

    e_study_balcinrenone_phase1_cl <- 1.06
    label("Fractional increase in CL/F for phase 1 study participants (unitless)")
    # Table 1 CL/F~Study = 1.06 (95% CI 0.685-1.44; bootstrap median 1.058,
    # bootstrap RSE 19.04%). Supplement control stream: CLPTSFLAG = 1 when
    # PTSFLAG = 1 (patients, the reference) and (1 + THETA(12)) when
    # PTSFLAG = 0 (phase 1). 1 + 1.06 = 2.06 reproduces the 2.06-fold
    # (95% CI 1.68-2.46) CL/F difference of Results 3.3, i.e. 0.48-fold
    # AUCss.

    e_fed_fdepot <- 0.134
    label("Fractional increase in relative bioavailability F1 in the fed state (unitless)")
    # Table 1 F1~Food = 0.134 (95% CI 0.0593-0.208; bootstrap median 0.136,
    # bootstrap RSE 37.88%). Table 1 footnote: 'Fed state F1 = 1 + F1~Food'.
    # 1 + 0.134 = 1.134 reproduces the 1.13-fold (95% CI 1.06-1.21) fed-vs-
    # fasted AUCss ratio of Results 3.3.

    e_fed_d1 <- 0.521
    label("Fractional increase in the zero-order duration D1 in the fed state (unitless)")
    # Table 1 D1~Food = 0.521 (95% CI 0.251-0.791; bootstrap median 0.528,
    # bootstrap RSE 75.47%). Results 3.2: 'longer duration in the fed state'.
    # Fed D1 = 0.631 * 1.521 = 0.960 h.

    e_fed_q <- -0.790
    label("Fractional change in Q/F in the fed state (unitless)")
    # Table 1 Q/F~Food = -0.790 (95% CI -0.817 to -0.764; bootstrap median
    # -0.792, bootstrap RSE 4.68%). Results 3.2: 'lower intercompartmental
    # clearance in the fed state'. Fed Q/F = 9.53 * 0.210 = 2.00 L/h; this
    # is what produces the lower fed-state Ctrough reported in Results 3.3.

    e_dose_balcinrenone_mg_ka <- -0.262
    label("Power exponent on (dose / 150 mg) for KA (unitless)")
    # Table 1 KA~Dose = -0.262 (95% CI -0.331 to -0.192; bootstrap median
    # -0.261, bootstrap RSE 23.15%). Table 1 footnote:
    # 'Dose KA = (dose/150)**KA~Dose'. Negative exponent = faster absorption
    # at lower doses, consistent with the solubility-limited absorption
    # described in the Introduction and Discussion.

    # ========================================================================
    # Between-subject variability -- Parkinson 2025 Table 1. Values are
    # log-scale VARIANCES for exponential (log-normal) IIV; Table 1's
    # footnote gives CV% of log-normal omegas = sqrt(exp(estimate) - 1) * 100,
    # which each bracketed CV% below reproduces. All omegas are diagonal
    # (the supplement's $OMEGA block declares no BLOCK structure). Results
    # 3.2: 'As the importance sampling method was used for estimation in the
    # final model, interindividual variability was included on all
    # parameters.'
    # ========================================================================
    etalcl ~ 0.148
    # Table 1 IIV-CL = 0.148 [CV% = 40.0] (95% CI 0.0943-0.202; shrinkage 21.1%).

    etalvc ~ 0.352
    # Table 1 IIV-Vc = 0.352 [CV% = 65.0] (95% CI 0.0555-0.649; shrinkage 48.9%).

    etalq ~ 0.525
    # Table 1 IIV-Q = 0.525 [CV% = 83.1] (95% CI 0.249-0.802; shrinkage 47.5%).

    etalvp ~ 1.21
    # Table 1 IIV-Vp = 1.21 [CV% = 153] (95% CI 0.685-1.73; shrinkage 41.0%).

    etalka ~ 0.00508
    # Table 1 IIV-KA = 0.00508 [CV% = 7.13] (95% CI -0.0187 to 0.0289;
    # shrinkage 85.3%). The confidence interval spans zero and shrinkage is
    # 85%, so this variance is poorly informed; it is retained because the
    # paper retained it.

    etalfdepot ~ 0.0281
    # Table 1 IIV-F1 = 0.0281 [CV% = 16.9] (95% CI -0.00178 to 0.0581;
    # shrinkage 69.2%).

    etald1 ~ 1.26
    # Table 1 IIV-D1 = 1.26 [CV% = 159] (95% CI 0.455-2.07; shrinkage 43.9%).
    # Discussion: 'the PK of balcinrenone was highly variable, especially in
    # the absorption phase, as judged by the large estimate of
    # between-subject variability for the zero-order infusion into the
    # dosing compartment'.

    # ========================================================================
    # Residual variability -- Parkinson 2025 Table 1, 'Residual variability'.
    # Table 1 reports VARIANCES; its footnotes give CV% of sigma =
    # sqrt(estimate) * 100 for the proportional term and SD = sqrt(estimate)
    # for the additive term. Supplement $ERROR:
    # Y = IPRED + IPRED*EPS(1) + EPS(2), i.e. combined proportional plus
    # additive on the linear scale.
    # ========================================================================
    propSd <- sqrt(0.191)
    label("Proportional residual SD (fraction)")
    # Table 1 Proportional = 0.191 [CV% = 43.7] (95% CI 0.173-0.209;
    # shrinkage 10.8%). sqrt(0.191) = 0.437.

    addSd <- sqrt(25.5)
    label("Additive residual SD (nmol/L)")
    # Table 1 Additive = 25.5 [SD = 5.05] (95% CI 2.95-48.0; shrinkage
    # 10.8%). sqrt(25.5) = 5.05 nmol/L, close to the 2 and 10 nmol/L assay
    # LLOQs of Methods 2.2.
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual parameters. Order of operations follows the supplement's
    #    final control stream exactly: the typical value is scaled first by
    #    allometry, then by the covariate multipliers, and the exponential
    #    eta is applied last.
    #
    #    Both binary covariates are recoded to the canonical polarity
    #    (FED = 1 - FOOD, STUDY_BALCINRENONE_PHASE1 = 1 - PTSFLAG) so the
    #    reference level -- a fasted phase 1b/2b heart-failure/CKD patient --
    #    is the all-zero record, exactly as in the control stream.
    # ----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (WT / 70)^e_wt_cl_q *
      (CRCL / 57.68)^e_crcl_cl *
      (1 + e_study_balcinrenone_phase1_cl * STUDY_BALCINRENONE_PHASE1)

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp

    q <- exp(lq + etalq) *
      (WT / 70)^e_wt_cl_q *
      (1 + e_fed_q * FED)

    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    ka <- exp(lka + etalka) *
      (DOSE_BALCINRENONE_MG / 150)^e_dose_balcinrenone_mg_ka

    fdepot <- exp(lfdepot + etalfdepot) * (1 + e_fed_fdepot * FED)

    d1 <- exp(ld1 + etald1) * (1 + e_fed_d1 * FED)

    # ----------------------------------------------------------------------
    # 2. Micro-constants. The control stream uses ADVAN5 with
    #    K12 = KA (depot -> central), K20 = CL/VC, K23 = Q/VC and
    #    K32 = Q/VP; the names below follow the nlmixr2lib convention where
    #    k12 / k21 are the central <-> peripheral1 rates.
    # ----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----------------------------------------------------------------------
    # 3. Two-compartment disposition with sequential zero-order then
    #    first-order absorption (Results 3.2). The dose enters `depot` at a
    #    constant rate over the window [0, D1] and `depot` then drains
    #    first-order at KA into `central`. Dose records must carry
    #    rate = -2 so rxode2 uses the model's dur(depot); a plain bolus
    #    collapses the zero-order phase and biases Cmax and Tmax.
    # ----------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    dur(depot) <- d1
    f(depot) <- fdepot

    # ----------------------------------------------------------------------
    # 4. Observation. Amounts are in nmol and vc in L, so central / vc is
    #    nmol/L, matching the bioanalytical assay units (supplement 'Summary
    #    of Bioanalytical Methods'). This reproduces the control stream's
    #    S2 = VC and CP = A(2)/S2.
    # ----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
