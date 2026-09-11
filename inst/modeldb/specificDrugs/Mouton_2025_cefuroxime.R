Mouton_2025_cefuroxime <- function() {
  description <- "Two-compartment population PK model for intravenous cefuroxime in critically ill adults admitted to the intensive care unit. Clearance is directly proportional to absolute (NOT BSA-indexed) MDRD estimated glomerular filtration rate, normalised to 93 mL/min, and carries no non-renal term; both volumes and the intercompartmental clearance scale allometrically on Janmahasatian fat-free mass derived inside the model from body weight, height and sex, normalised to a 58.2 kg reference (exponent 1 fixed on the volumes, 0.75 fixed on Q). Note that clearance itself carries no allometric term in the final model. Fitted by NONMEM 7.5.0 FOCE+I to TOTAL (protein-bound plus unbound) plasma concentrations with a proportional residual error and interindividual variability on CL and V1 only. The model additionally returns Cu, the unbound cefuroxime concentration, from the saturable Thonnings albumin-binding relationship the authors applied to derive their probability-of-target-attainment results; the unbound fraction is Cu / Cc. Mouton 2025, n = 20 patients, 105 analysed plasma samples over a single dosing interval."
  reference <- "Mouton JWA, Machiels JD, Pistorius AMA, ter Heine R, Frenzel T, Jager NGL, Schouten JA, Janssen PKC, Aarnoutse RE, Bruggemann RJ. Population pharmacokinetics and optimized dosing of cefuroxime in critically ill patients. Br J Clin Pharmacol. 2025;91(9):2755-2761. doi:10.1002/bcp.70144. PMC12381626. Structural model, final parameter estimates and the unbound-concentration relationship are all taken from the Supporting Information (Table S2 and the two NONMEM control streams); the main article reports no parameter table."
  vignette <- "Mouton_2025_cefuroxime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefuroxime was given as an intravenous bolus over 5 min
  # and TOTAL (protein-bound plus unbound) plasma cefuroxime was assayed by a
  # validated UPLC-MS/MS method with an LLOQ of 0.1 mg/L (Methods 2.5). The
  # model therefore predicts total, not unbound, concentration; Cu in model()
  # is derived from the total by an external binding relationship.
  compartmentData <- list(
    central     = list(analyte = "cefuroxime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefuroxime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the MDRD equation, re-expressed for standardised creatinine assays, reported as an ABSOLUTE (not BSA-indexed) value in mL/min",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The ONLY covariate retained in the final model, and only on",
        "clearance. Enters as the bare linear ratio (CRCL / 93) -- the",
        "supplement's Model development equation (4),",
        "CL = (eGFR / eGFR_mean) * theta_CL, transcribed in the estimation",
        "control stream $PK block as TVCL = THETA(1) * MDRDABS / 93. The",
        "relationship is a strict proportionality with unit slope, so",
        "e_crcl_cl is fixed at 1 rather than estimated, and the model has NO",
        "non-renal clearance intercept: the supplement's alternative",
        "equation (5), CL = theta_CL_nonrenal + (eGFR / eGFR_mean) *",
        "theta_CL_renal, was tested during model development but is not the",
        "final form. A consequence is that predicted clearance goes to zero",
        "as CRCL goes to zero, which is an extrapolation the 24-168 mL/min",
        "observed range does not support; see the vignette Assumptions and",
        "deviations.",
        "SIZE NORMALISATION: this is the RAW un-normalised mL/min variant of",
        "the CRCL canonical, NOT the mL/min/1.73 m^2 default. The paper is",
        "explicit and repeated about this ('the absolute (not indexed by BSA)",
        "values', Methods 2.4; Table 1 row 'MDRD absolute, mL/min'). Supplying",
        "a BSA-indexed value would silently rescale renal clearance. The",
        "supplement gives the conversion used: eGFR_absolute =",
        "eGFR_normalized * BSA / 1.73.",
        "The 93 mL/min denominator is the cohort mean absolute MDRD, and is",
        "confirmed independently by the Table S2 row label",
        "'CL (L/h/93mL/min)'. It is NOT the Table 1 median, which is 90",
        "[60-117.5] mL/min (Table 1 reports medians, the normalisation uses",
        "the mean).",
        "ASSAY CHOICE: MDRD absolute, CKD-EPI absolute and a measured 24-h",
        "urine creatinine clearance were all tested as the renal-function",
        "estimator and did not differ significantly in model performance.",
        "MDRD was chosen for the final model because it had the best fit and",
        "is the most commonly used parameter in clinical practice",
        "(Supporting Information, Final Model). CKD-EPI was 2.77 objective",
        "function units worse and changed the clearance estimate by 0.6%, so",
        "the authors state the results are also applicable to CKD-EPI",
        "(Discussion). A user may therefore supply a CKD-EPI absolute value",
        "in this column with negligible bias.",
        "Time-fixed: demographic and biochemical data were obtained once, on",
        "the pharmacokinetic sampling day (Methods 2.4). Must be strictly",
        "positive."
      ),
      source_name        = "MDRDABS"
    ),
    WT = list(
      description        = "Total body weight on the pharmacokinetic sampling day",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Not used directly as a covariate. WT enters only through the",
        "Janmahasatian fat-free mass derived inside model() from WT, HT and",
        "SEXF; FFM then scales V1, V2 and Q. Body weight was included a",
        "priori (not by covariate selection) because fat-free mass best",
        "describes the influence of weight for a hydrophilic drug such as",
        "cefuroxime (Methods 2.6). Table 1: median 85 [75-100] kg, range",
        "55-120 kg; the Results paragraph additionally gives mean 85.2",
        "(SD 18.0) kg. Must be strictly positive."
      ),
      source_name        = "WT"
    ),
    HT = list(
      description        = "Body height on the pharmacokinetic sampling day",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Not used directly as a covariate; enters only through the",
        "Janmahasatian fat-free mass derived in model(). The control stream",
        "computes BMI = WT / (HT/100)**2, so HT is supplied in CENTIMETRES",
        "and the model divides by 100 internally -- supplying metres would",
        "inflate BMI by a factor of 10^4 and collapse FFM. Height is not",
        "tabulated in Table 1; the simulated typical patient used 172 cm",
        "(Methods 2.7). Must be strictly positive."
      ),
      source_name        = "HT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Not used directly as a covariate; enters only through the",
        "sex-specific coefficients of the Janmahasatian fat-free mass",
        "equation. VALUE TRANSFORMATION: the source control stream codes",
        "SEX = 1 for male and SEX = 0 for female (stated verbatim in the",
        "$PK comment '; Male: SEX=1 Female: SEX=0'), which is the inverse of",
        "the SEXF canonical, so SEXF = 1 - SEX. The substitution is exact",
        "and was verified algebraically rather than assumed: the control",
        "stream writes the denominator as",
        "(8780 - SEX*2100) + ((244 - SEX*28) * BMI), which under SEX =",
        "1 - SEXF becomes (6680 + 2100*SEXF) + ((216 + 28*SEXF) * BMI) --",
        "i.e. 6680 + 216*BMI for males (SEXF = 0) and 8780 + 244*BMI for",
        "females (SEXF = 1), the published Janmahasatian coefficients. The",
        "model is written in the sex-branch-interpolation form used by",
        "Hughes_2024_vancomycin_parametric.R for the identical equation.",
        "No effect-coefficient sign is involved, because sex enters a",
        "structural body-composition formula rather than a fitted covariate",
        "effect. Table 1: 10 of 20 patients (50%) were male."
      ),
      source_name        = "SEX"
    )
  )

  # Covariates the source SCREENED but did not retain in the final model.
  # Documentation only -- none is referenced in model(). All were assessed by
  # visual inspection of empirical Bayes estimates of the PK parameters
  # against the covariate, with forward selection (p < 0.05) and backward
  # elimination (p > 0.001) for any relationship that looked present
  # (Methods 2.6 and Supporting Information, Model development). The
  # Supporting Information "Final Model" paragraph records the outcome: "No
  # other covariates were identified during the visual inspection".
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened, not retained. Table 1: median 69 [65-75] years, range 29-86; Results gives mean 66 (SD 14.0) years."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened as a covariate in its own right and not retained. BMI is",
        "nevertheless computed inside model() as the intermediate bmi_i,",
        "because the Janmahasatian fat-free-mass equation is parameterised",
        "on it; that internal use is a body-composition formula, not a",
        "covariate effect, which is why BMI is listed here rather than in",
        "covariateData. Table 1: median 26.1 [24.8-32.8] kg/m^2, range",
        "20.7-46.9."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a covariate on the PK parameters and not retained.",
        "Table 1: median 24.0 [19.75-26.25] g/L, range 5-31 -- markedly",
        "hypoalbuminaemic, as expected in an ICU cohort. Albumin does NOT",
        "enter the unbound-concentration relationship implemented in",
        "model() either: the Thonnings binding constants are fixed at the",
        "healthy-population values. Supplementary Figure S1 explores scaling",
        "the binding capacity by the albumin ratio 23/45 and reports that",
        "halving albumin raises the free fraction by only about 0.12 around",
        "the 8 mg/L MIC cutoff, which the authors judge 'too small to be of",
        "practical importance'. See the vignette Assumptions and deviations.",
        "The ALB column was present in the estimation dataset ($INPUT) and",
        "output ($TABLE) but is unused by the final $PK block."
      )
    ),
    WBC = list(
      description = "Leukocyte count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened, not retained. Not tabulated in Table 1."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score at ICU admission",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened, not retained. Table 1: median 21 [13-24], range 8-34. Present in the estimation dataset as APA2."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score on the day of pharmacokinetic sampling",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened, not retained. Table 1: median 8 [6-10], range 3-17. Present in the estimation dataset as SOFA."
    ),
    DOSE = list(
      description = "Administered cefuroxime dose level (750 mg or 1500 mg q8h)",
      units       = "mg",
      type        = "continuous",
      notes       = paste(
        "Screened as a covariate and not retained: 'In this study, dose did",
        "not significantly impact PK parameters' (Discussion). 16 of 20",
        "patients received 750 mg q8h for selective digestive",
        "decontamination and 4 received 1500 mg q8h for empirical sepsis",
        "therapy. The authors therefore treat cefuroxime PK as linear and",
        "extrapolate to the higher simulated doses on that basis."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20,
    n_studies      = 1,
    age_range      = "29-86 years",
    age_median     = "69 years",
    weight_range   = "55-120 kg",
    weight_median  = "85 kg",
    sex_female_pct = 50,
    disease_state  = "critically ill adults admitted to the intensive care unit, receiving cefuroxime as standard care either for selective digestive decontamination or as empirical antibiotic therapy for sepsis",
    dose_range     = "750 mg or 1500 mg intravenously three times daily, each administered as a 5-min bolus",
    regions        = "Netherlands (single centre: Radboud University Medical Center, Nijmegen)",
    renal_function = "MDRD absolute median 90 [60-117.5] mL/min, range 24-168; CKD-EPI absolute median 93.5 [60-109.8] mL/min, range 22-142; measured 24-h creatinine clearance median 113 [72-150.8] mL/min, range 27-240. Augmented renal clearance was present in part of the cohort and every ARC patient had traumatic brain injury, subarachnoid haemorrhage or burns (Discussion).",
    notes          = paste(
      "Baseline demographics are Table 1 of the main article; the Results",
      "'Patient characteristics and sampling' paragraph adds the means.",
      "Single-centre prospective observational PK study, ClinicalTrials.gov",
      "NCT04470973. Sampling was pre-dose and at 0.5, 1, 3, 5 and 8 h after",
      "an intravenous dose within a SINGLE dosing interval, in the first 72 h",
      "of therapy (Methods 2.4); 120 samples were drawn and 15 excluded",
      "(drawn before or during administration, or with no recorded sampling",
      "time), leaving 105 for analysis. Observed total concentrations ranged",
      "0.2-112.8 mg/L. The sample size was set by stochastic simulation and",
      "re-estimation to estimate clearance with bias and imprecision below",
      "15% (Methods 2.2). All patients had a central venous or arterial",
      "catheter and were older than 18 years."
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # SOURCE OF THE FINAL ESTIMATES.
    #
    # The main article is a Short Communication and contains NO parameter
    # table -- Methods 2.6 defers to "Supporting Information". Two independent
    # sources in the supplement agree, and both are used below:
    #
    #   (a) Table S2 "Final model parameters including covariates", which
    #       prints the estimates rounded to 2-3 significant figures together
    #       with RSEs and a 1000-sample bootstrap mean and 95% CI; and
    #   (b) the second NONMEM control stream, headed "NONMEM control stream
    #       simulations" (based on run052), whose $THETA / $OMEGA / $SIGMA
    #       blocks carry the final estimates at full precision -- this is the
    #       stream that generated the paper's Figures 1 and 2.
    #
    # The values below are taken from (b) because it is more precise, and
    # every one is cross-checked against (a) in the trailing comment. This is
    # NOT the "initial estimates" trap: the FIRST control stream in the
    # supplement (based on run026) is the estimation run, and its $THETA
    # block holds bounded starting values -- (0, 9.0, 20) for CL, (0, 10.0,
    # 20) for V1, (0, 13.0, 50) for V2, (0, 18.0, 50) for Q, with $OMEGA 0.1 /
    # 0.1 and $SIGMA 0.1. Those starting values are deliberately NOT used
    # here; they disagree with Table S2, which is the tell.
    # ------------------------------------------------------------------------

    # Structural parameters. The reference subject is one with an absolute
    # MDRD eGFR of 93 mL/min and a fat-free mass of 58.2 kg.
    lcl <- log(7.95); label("Clearance at CRCL = 93 mL/min (L/h)")
    # Simulation control stream $THETA 7.95 ; CL renal.
    # Table S2 "CL (L/h/93mL/min)" = 8.0, RSE 5.4%, bootstrap mean 8.0
    # (95% CI 7.1-8.8).
    lvc <- log(6.26); label("Central volume V1 at FFM = 58.2 kg (L)")
    # Simulation control stream $THETA 6.26 ; V1.
    # Table S2 "V1 (L)" = 6.3, RSE 17.9%, bootstrap mean 6.4 (95% CI 4.2-9.3).
    lvp <- log(14.3); label("Peripheral volume V2 at FFM = 58.2 kg (L)")
    # Simulation control stream $THETA 14.3 ; V2.
    # Table S2 "V2 (L)" = 14.3, RSE 6.1%, bootstrap mean 14.3
    # (95% CI 12.1-16.5).
    lq <- log(26.5); label("Intercompartmental clearance Q at FFM = 58.2 kg (L/h)")
    # Simulation control stream $THETA 26.5 ; Q.
    # Table S2 "Q (L/h)" = 26.5, RSE 25.7%, bootstrap mean 26.4
    # (95% CI 12.3-40.1). This is the least well identified parameter in the
    # model; the bootstrap CI spans more than a threefold range.

    # Covariate effects. Both are structurally fixed, not estimated: no
    # exponent for either appears in Table S2, and both are written as
    # literal constants in the control stream.
    e_crcl_cl <- fixed(1); label("Exponent of (CRCL / 93) on CL (unitless; a value of 1 makes CL strictly proportional to renal function)")
    # Estimation and simulation control streams, $PK:
    #   TVCL = THETA(1) * MDRDABS / 93
    # i.e. a bare linear ratio. Written here as a power with the exponent
    # fixed at 1 so that the "fixed, not estimated" status is explicit and
    # the term is overridable; (CRCL/93)^1 is identically (CRCL/93).
    e_ffm_vc <- fixed(1); label("Allometric exponent of (FFM / 58.2) on V1 (unitless)")
    # Estimation and simulation control streams, $PK: ALLOV = (FFM/58.2),
    # TVV1 = THETA(2) * ALLOV. Also Supporting Information, Model development
    # equation (1): TVV1 = (FFM/58.2) * theta_V1. Methods 2.6 states the
    # exponent 1.0 for volumes was imposed a priori, citing Anderson &
    # Holford; it was not estimated.
    e_ffm_vp <- fixed(1); label("Allometric exponent of (FFM / 58.2) on V2 (unitless)")
    # Same ALLOV term: TVV2 = THETA(3) * ALLOV; Model development equation (2).
    e_ffm_q <- fixed(0.75); label("Allometric exponent of (FFM / 58.2) on Q (unitless)")
    # Estimation and simulation control streams, $PK: ALLOCL = (FFM/58.2)**0.75,
    # TVQ = THETA(4) * ALLOCL. Also Model development equation (3):
    # TVQ = (FFM/58.2)^0.75 * theta_Q. Methods 2.6: exponent 0.75 for
    # clearances, imposed a priori.
    #
    # NOTE that ALLOCL is applied ONLY to Q. Despite the Methods 2.6 sentence
    # "FFM was incorporated using an allometric exponent of 0.75 for
    # clearances", the final $PK block scales CL by renal function alone
    # (TVCL = THETA(1)*MDRDABS/93, with no ALLOCL factor). The control stream
    # is the operative source here and is unambiguous; there is no e_ffm_cl
    # in this model. See the vignette Assumptions and deviations.

    # ------------------------------------------------------------------------
    # Interindividual variability. NONMEM $OMEGA values are VARIANCES on the
    # log scale (exponential IIV, CL = TVCL*EXP(ETA(1))), which the Table S2
    # footnote confirms by giving the back-transform it used:
    # "%CV = 100*sqrt(exp(omega^2) - 1)". Both round-trip exactly:
    #   0.054 -> 100*sqrt(exp(0.054)-1) = 23.6% -> Table S2 "24"
    #   0.373 -> 100*sqrt(exp(0.373)-1) = 67.2% -> Table S2 "67"
    # so the omega scale is settled by arithmetic, not assumed.
    #
    # The supplement's Model development section states directly that
    # "Inter-individual variability was assumed to be log-normally
    # distributed", and the Final Model section that IIV was estimated on
    # clearance and the central volume only.
    # ------------------------------------------------------------------------
    etalcl ~ 0.054
    # Simulation control stream $OMEGA 0.054 ; IIV CL. Table S2 "IIV on CL
    # (% CV)" = 24, shrinkage 3.9%, bootstrap mean 23 (95% CI 11-33).
    etalvc ~ 0.373
    # Simulation control stream $OMEGA 0.373 ; IIV V1. Table S2 "IIV on V1
    # (% CV)" = 67, shrinkage 17.4%, bootstrap mean 59 (95% CI 32-96). The
    # paper flags this in the Discussion: "Our results suggest a high
    # interindividual variability of the CL and V1 parameters, as expected in
    # the critically ill patient population."
    #
    # ETA(3) on V2 and ETA(4) on Q are declared in both control streams as
    # "0 FIX" and are correspondingly absent from Table S2. They are OMITTED
    # here rather than written as `~ fixed(0)`: a zero diagonal makes the
    # omega matrix singular, and rxode2 then fails to Cholesky-decompose it
    # when simulating a cohort. Omitting the eta is numerically identical
    # (V2 and Q take their typical values for every subject) and is the same
    # choice made by Assmus_2025_benznidazole_mouse.R.

    # Residual error. Proportional only -- the supplement's Final Model
    # section states "The residual error of the final model was best
    # described with a proportional error model", and the $ERROR block is
    # Y = IPRED*(1 + ERR(1)) with a single $SIGMA. NONMEM $SIGMA is a
    # VARIANCE, so the nlmixr2 SD is its square root:
    # sqrt(0.0483) = 0.21977, and 100*sqrt(0.0483) = 22.0% reproduces the
    # Table S2 "Proportional (% CV)" row of 22 exactly.
    propSd <- 0.2198; label("Proportional residual error (fraction)")
    # Simulation control stream $SIGMA 0.0483 ; PROP ERR. Table S2
    # "Proportional (% CV)" = 22, RSE 13.8%, bootstrap mean 21
    # (95% CI 14-27).
  })

  model({
    # 1. Derived body-composition terms, transcribed verbatim from the $PK
    #    block of both control streams:
    #      BMI = WT/(HT/100)**2
    #      FFM = (9270*WT) / ((8780 - SEX*2100) + ((244 - SEX*28)*BMI))
    #    with SEX = 1 for male and SEX = 0 for female. Under the canonical
    #    SEXF = 1 - SEX this is the standard Janmahasatian et al. (Clin
    #    Pharmacokinet 2005;44:1051-1065) pair -- the control stream cites
    #    PMID 16176118 inline -- written here in the sex-branch-interpolation
    #    form so that SEXF selects the female branch:
    #      male   (SEXF = 0): FFM = 9270*WT / (6680 + 216*BMI)
    #      female (SEXF = 1): FFM = 9270*WT / (8780 + 244*BMI)
    bmi_i    <- WT / (HT / 100)^2
    ffm_male <- 9270 * WT / (6680 + 216 * bmi_i)
    ffm_fem  <- 9270 * WT / (8780 + 244 * bmi_i)
    ffm_i    <- ffm_male + SEXF * (ffm_fem - ffm_male)

    # 2. Individual PK parameters. CL scales on renal function only (no
    #    allometric term, no non-renal intercept), the two volumes and Q on
    #    fat-free mass only. Only CL and V1 carry an eta.
    #
    #    The 58.2 kg FFM reference is hardcoded in the control stream, whose
    #    comment describes it as "FFMref=58.2 (male; HT=1.80; WT=70)". That
    #    parenthesis does not reproduce: the Janmahasatian equation above
    #    returns 57.19 kg for a 70 kg, 180 cm male, and 58.2 kg corresponds
    #    to about 184 cm. The discrepancy is in the source's comment, not in
    #    its arithmetic -- 58.2 is the literal divisor the model was fitted
    #    with, so 58.2 is used here and the Table S2 estimates are the values
    #    at FFM = 58.2 kg. Recorded in the vignette Errata.
    cl <- exp(lcl + etalcl) * (CRCL / 93)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (ffm_i / 58.2)^e_ffm_vc
    vp <- exp(lvp)          * (ffm_i / 58.2)^e_ffm_vp
    q  <- exp(lq)           * (ffm_i / 58.2)^e_ffm_q

    # 3. Micro-constants, exactly as the $PK block defines them
    #    (K10 = CL/V1, K12 = Q/V1, K21 = Q/V2).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. The source used $SUBROUTINE ADVAN5 with
    #    COMP=(CENTRAL, DEFDOSE) and COMP=(PERI); the equivalent
    #    two-compartment system is written out explicitly here.
    #
    #    Infusion duration is deliberately NOT set in the model. The
    #    estimation control stream fixes D1 = 0.083 h because every observed
    #    dose was the same 5-min bolus, but the simulation control stream --
    #    the one that produced Figures 1 and 2 -- drops D1 and takes RATE
    #    from the dataset, because the simulated regimens include 4-h
    #    extended infusions and 24-h continuous infusions. Duration therefore
    #    belongs in the event table (rxode2 `dur` or `rate`), which keeps all
    #    three published regimens reachable.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 5. Observation. S1 = V1 in the control stream, so a dose in mg over a
    #    volume in L gives mg/L. This is the TOTAL (protein-bound plus
    #    unbound) concentration, which is what the assay measured and what
    #    the model was fitted to.
    Cc <- central / vc

    # 6. Unbound concentration.
    #
    #    No unbound assay was available (Methods 2.7), so the authors
    #    predicted free cefuroxime from the total with a published saturable
    #    binding relationship rather than a single fixed free fraction. The
    #    $ERROR block of the simulation control stream implements it as:
    #      CBMAX = 23.47
    #      KBTB  = 0.02126
    #      A     = CBMAX - IPRED + (1/KBTB)
    #      B     = -IPRED/KBTB
    #      CFREE = (SQRT(A**2 - 4*B) - A)/2
    #    which is the positive root of the binding equilibrium
    #    Ctotal = Cfree + Cbmax*Cfree/(1/(kb*Tb) + Cfree), after Thonnings et
    #    al., J Med Microbiol 2020;69:387-395 (the supplement's reference 9).
    #    Cbmax is the maximal binding capacity and kb*Tb the product of the
    #    binding affinity and the mean time in the bound state; the control
    #    stream records the molecular weights it used to reach these numbers
    #    (albumin 66.4e3 g/mol, cefuroxime 424.4 g/mol).
    #
    #    The two constants are written inline rather than as ini() parameters
    #    because they are fixed constants of an EXTERNAL published formula,
    #    not parameters of this paper's model: they appear nowhere in
    #    Table S2, were not estimated here, and are applied downstream of the
    #    fit. This is the same treatment the Janmahasatian constants get in
    #    step 1 above and in Hughes_2024_vancomycin_parametric.R.
    #
    #    Cu is derived from Cc, which carries no residual error, exactly
    #    matching the source's use of IPRED. At Cc = 0 the expression returns
    #    0 rather than a division-by-zero, so pre-dose records are safe; the
    #    unbound FRACTION, where it is wanted, is Cu / Cc for Cc > 0.
    bind_a <- 23.47 - Cc + 1 / 0.02126
    bind_b <- -Cc / 0.02126
    Cu     <- (sqrt(bind_a^2 - 4 * bind_b) - bind_a) / 2

    # 7. Residual error, on the total concentration.
    Cc ~ prop(propSd)
  })
}
