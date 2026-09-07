Liu_2024_saf189s <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model with lag-time ",
    "first-order oral absorption and TIME-DEPENDENT (autoinhibited) ",
    "apparent clearance for SAF-189s, an investigational ",
    "second-generation ALK/ROS1 tyrosine kinase inhibitor, in Chinese ",
    "adults with ALK-positive or ROS1-positive advanced non-small cell ",
    "lung cancer plus healthy Chinese volunteers (Liu 2024, n = 317 ",
    "subjects and 3,173 concentrations pooled across the phase I food-",
    "effect study STL31147 and the phase I/II study SAF001, ",
    "NCT04237805; 20-210 mg orally once daily). Total apparent ",
    "clearance decays STEPWISE from day to day toward a constant: ",
    "CL/F = CL1/F + CL2/F * exp(-Kout * (DAY - 1)) with DAY the integer ",
    "post-administration day, so CL/F falls from 104.3 L/h on day 1 to ",
    "64 L/h at steady state (Liu 2024 Eq. 2, Figure 1). That is the ",
    "PK signature of SAF-189s being both a substrate and an inhibitor ",
    "of P-gp and a CYP3A substrate; the corresponding terminal ",
    "half-life lengthens from 34.6 h after the first dose to 56.4 h at ",
    "steady state. Age enters CL/F as a power term referenced to the ",
    "cohort median 53 years, and prior anti-cancer therapy in ALK+ ",
    "patients (the paper's ALKPOT classification) enters V/F as three ",
    "multiplicative indicator effects. Interindividual variability is ",
    "carried on V/F, KA, CL/F and the absorption lag time; residual ",
    "error is combined proportional plus additive. Seven companion ",
    "exposure-response models in the Liu_2024_saf189s_* family consume ",
    "the exposure metrics this model generates."
  )
  reference <- paste(
    "Liu Y, Tan Y, Hu L, Li J, Yang J, Diao L, Yang J.",
    "Population pharmacokinetics and exposure-response analyses of",
    "SAF-189s in Chinese patients with ALK+/ROS1+ non-small cell lung",
    "cancer.",
    "Front Pharmacol. 2024;15:1418549.",
    "doi:10.3389/fphar.2024.1418549.",
    sep = " "
  )
  vignette <- "Liu_2024_saf189s"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    AGE = list(
      description        = "Age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL/F as the power term (AGE / 53)^-0.314, referenced to",
        "the cohort MEDIAN age of 53 years (Liu 2024 Results, PopPK",
        "analysis: 'mean age 51.6 years, median age 53 years'; Figure 4",
        "caption: 'Compared with the median age (53 years)'). The",
        "reference value is pinned independently by the Figure 4 forest",
        "plot: (25/53)^-0.314 = 1.266, i.e. CL/F 27% higher at the 5th",
        "percentile of 25 years, and (71/53)^-0.314 = 0.912, i.e. 8.8%",
        "lower at the 95th percentile of 71 years -- both matching the",
        "printed +27% and -8.8% exactly. The exponent is NEGATIVE, so",
        "younger subjects clear SAF-189s faster. Pooled cohort age",
        "median 53 years, range 18.0-84.0 (Supplementary Table 1); 272 of",
        "317 subjects (85.8%) were under 65 years. The healthy-volunteer",
        "arm is much younger (median 27 years) than the patient arms",
        "(51-54 years), so age and subject type are partially confounded",
        "in this cohort."
      ),
      source_name        = "Age"
    ),
    PRIOR_ALKI = list(
      description        = "Prior ALK-inhibitor therapy indicator in ALK-positive patients; 1 = the patient received ALK-inhibitor treatment before enrolment, 0 = ALK-inhibitor-naive.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no previous ALK-inhibitor treatment; the paper's ALKPOT = 1 group, n = 115, 36.2% of the pooled cohort)",
      notes              = paste(
        "Collapses levels 2, 3 and 4 of the paper's six-level ALKPOT",
        "classification, which Liu 2024 itself groups into a single",
        "estimated effect (Table 2 row 'ALKPOT = 2,3,4 on V/F'). The",
        "underlying levels are 2 = intolerant to crizotinib as the only",
        "prior treatment (n = 45), 3 = one prior second- or",
        "third-generation ALK inhibitor (n = 10), 4 = at least two prior",
        "second- or third-generation ALK inhibitors (n = 9); together",
        "n = 64, 20.2% (Supplementary Table 1). Set to 0 for ROS1+ /",
        "unknown patients and for healthy subjects -- those groups carry",
        "their own indicators, and the three ALKPOT indicators are",
        "mutually exclusive by construction. Multiplies V/F by 0.734",
        "(26.6% lower apparent volume)."
      ),
      source_name        = "ALKPOT = 2, 3, 4"
    ),
    TUM_ALK_MUT = list(
      description        = "Tumour ALK rearrangement status; 1 = ALK-positive, 0 = not ALK-positive (ROS1-positive, unknown, or a healthy participant with no tumour).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (ALK-positive; the reference cell of the V/F covariate model is an ALK-positive, ALK-inhibitor-naive patient)",
      notes              = paste(
        "Used together with DIS_HEALTHY to reconstruct the paper's",
        "ALKPOT = 5 group ('others (ROS1+ patients or unknown)',",
        "n = 114, 36.0%), which in model() is the product",
        "(1 - TUM_ALK_MUT) * (1 - DIS_HEALTHY). That group is a MIXED",
        "bucket -- Liu 2024 pools genuinely ROS1-positive patients with",
        "patients of unknown status -- so it cannot be decomposed into a",
        "clean ROS1-positive indicator, and no TUM_ROS1_MUT column is",
        "introduced here. Supplementary Table 1 'ALK Mutation' row:",
        "positive 190 (60%), negative or other 127 (40%); the healthy",
        "volunteers are counted in the negative-or-other cell, which is",
        "why DIS_HEALTHY must gate this term. Multiplies V/F by 0.917",
        "(8.3% lower apparent volume) when the subject is a",
        "non-ALK-positive patient."
      ),
      source_name        = "ALK Mutation (Positive / Negative or Other); ALKPOT = 5"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant indicator; 1 = healthy Chinese volunteer from study STL31147, 0 = NSCLC patient from study SAF001.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NSCLC patient; the reference cell is an ALK-positive, ALK-inhibitor-naive patient)",
      notes              = paste(
        "The paper's ALKPOT = 6 level (n = 24, 7.6%). Multiplies V/F by",
        "0.784 (21.6% lower apparent volume) relative to",
        "ALK-inhibitor-naive ALK+ patients. Note the confounding: the",
        "healthy arm is a single-dose 160 mg crossover in young men",
        "(median age 27 years, 95.8% male), so this effect is not a clean",
        "health-status contrast. Liu 2024 Discussion nonetheless reports",
        "that 'healthy vs cancer patient' had no clinically meaningful",
        "effect on systemic exposure, because the AUCss consequence of",
        "the V/F shift is only -2.1% (Figure 4)."
      ),
      source_name        = "ALKPOT = 6 (healthy subjects); Subject Type"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Prespecified and screened in the stepwise covariate model (Liu 2024 Methods, Covariate analysis) but NOT retained: 'None of the other intrinsic factors (i.e., healthy vs cancer patient, bodyweight, sex, preexisting mild hepatic impairment, and preexisting mild or moderate renal impairment) or extrinsic factors (i.e., concomitant medications) had clinically meaningful effects on SAF-189s systemic exposure.' No point estimate is reported anywhere on disk. Pooled cohort median 63.2 kg, range 37.3-92.5 (Supplementary Table 1)."
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened, not retained (Liu 2024 Results, Impact of the PK covariates). Supplementary Table 1 reports 153 of 317 female (48.2%); the Results narrative instead says 'approximately half of them were women (56.5%)'. The two disagree and the count-backed 48.2% is used in this model's population metadata -- see the vignette Errata."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as a renal-function marker, not retained. Pooled median 96.2 mL/min, range 48.8-199 (Supplementary Table 1); 94 subjects (29.6%) had mild and 1 (0.4%) moderate renal impairment."
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened, not retained. Pooled median 41.2 g/L, range 23.0-52.6 (Supplementary Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker, not retained. Pooled median 24.1 U/L, range 2.00-127 (Supplementary Table 1); 41 subjects (13%) had mild hepatic impairment."
    ),
    ALP = list(
      description = "Alkaline phosphatase at baseline.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker, not retained. Pooled median 109 U/L, range 46.0-493 (Supplementary Table 1)."
    ),
    TBIL = list(
      description = "Total bilirubin at baseline.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker, not retained. Supplementary Table 1 prints median 68 umol/L (range 35.0-110), which is not a credible bilirubin distribution for a cohort described as 87% hepatically normal and looks like a row misalignment in that table -- see the vignette Errata. The screening outcome is unaffected because the covariate was dropped."
    ),
    SCR = list(
      description = "Serum creatinine at baseline.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened, not retained. Supplementary Table 1 prints median 163 umol/L (range 139-190), which contradicts the same table's creatinine clearance of 96.2 mL/min and its 65% normal-renal-function count; treated as a misaligned row -- see the vignette Errata."
    ),
    CONMED_ANY = list(
      description = "Any concomitant medication present at baseline.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as an extrinsic factor, not retained (Liu 2024 Methods lists 'co-administration values' among the prespecified covariates). Supplementary Table 1: present in 215 subjects (55.6%), including metformin in 46 (13.6%)."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "SAF-189s", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "SAF-189s", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 317L,
    n_studies      = 2L,
    n_observations = "3,173 measurable plasma concentrations retained of 3,538 acquired (89.68%); 329 pre-dose BQL, 6 pre-dose non-BQL, 30 post-dose BQL samples and 6 subjects with no samples were excluded",
    age_range      = "18.0-84.0 years (pooled; mean 51.6, median 53). Healthy volunteers median 27 (18.0-45.0); SAF001 phase I median 51 (28.0-68.0); SAF001 phase II median 54.1 (20.0-84.0)",
    age_median     = "53 years",
    weight_range   = "37.3-92.5 kg",
    weight_median  = "63.2 kg",
    sex_female_pct = 48.2,
    race_ethnicity = c(Asian = 100),
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer (190 ALK+ only, 103 ROS1+ only), plus 24 healthy volunteers; 136 of 317 (43%) had brain metastases at enrolment and 221 (69.8%) were disease stage IV",
    dose_range     = "SAF-189s 20, 40, 80, 120, 160 or 210 mg orally once daily in 21-day cycles (SAF001, with a single-dose 3-day PK run-in); a single 160 mg oral dose in each of two crossover periods (STL31147, fed and fasted)",
    regions        = "China",
    renal_function = "207 normal (65.2%), 94 mild (29.6%), 1 moderate (0.4%), 15 with eGFR >= 130 (4.8%)",
    hepatic_function = "276 normal (87%), 41 mild dysfunction (13%); no moderate or severe",
    notes          = paste0(
      "Baseline characteristics from Liu 2024 Supplementary document 1, ",
      "Supplementary Table 1 (a separate publisher file, Table1.docx). ",
      "Sex is reported inconsistently: the count-backed Supplementary ",
      "Table 1 gives 153 of 317 female (48.2%) while the Results ",
      "narrative says 56.5%; the count is used here. Food had no ",
      "meaningful effect (geometric mean ratio Cmax 109.1%, AUC0-t ",
      "105.1%), so the fed and fasted STL31147 periods are pooled ",
      "without a food covariate."
    )
  )

  ini({
    # ==================================================================
    # All structural values are the final-model typical values in Liu
    # 2024 Table 2, "Final model / Typical value / Estimate" column,
    # cross-checked against the Results narrative ("the typical (+/-SE)
    # values of CL1/F, CL2/F, V/F, KA, and tlag for SAF-189s were 64
    # (+/-1.51) L/h, 40.3 (+/-6.89) L/h, 5,210 (+/-339) L, 0.501
    # (+/-0.0444) 1/h, and 0.483 (+/-0.0306) h ... Kout was estimated to
    # be 1.35 (+/-0.421) 1/day"). Nothing here is fixed: every theta
    # carries an SE and a bootstrap 95% CI in Table 2.
    #
    # UNIT ERRATA IN TABLE 2 (labels only; the values are used as
    # printed). Table 2 prints "KA (L h-1)" and "Kout (L d-1)". Both are
    # first-order rate constants, so the litre is spurious: KA is 1/h
    # and Kout is 1/day. The Results narrative confirms both
    # ("0.501 ... h-1", "1.35 ... day-1").
    # ==================================================================

    lvc <- log(5210)  ; label("Apparent central volume of distribution V/F (L)")                      # Liu 2024 Table 2 theta1 = 5,210 L, SE 339, RSE 6.5%; bootstrap 5,199 (95% CI 4,415-7,373)
    lka <- log(0.501) ; label("First-order absorption rate constant KA (1/h)")                        # Liu 2024 Table 2 theta2 = 0.501, SE 0.0444, RSE 8.9%; bootstrap 0.504 (95% CI 0.445-0.585). Table 2 unit label "L h-1" is an erratum; KA is 1/h

    # ------------------------------------------------------------------
    # Time-dependent clearance, Liu 2024 Eq. 2:
    #   CL_i = TVCL1 + TVCL2 * exp(-Kout * (DAY - 1)) * exp(eta_i)
    # Registered `cl_exp_` role names (parameter-names.md, "Time-varying
    # clearance"): the constant the curve decays TO is cl_exp_inf and
    # the decaying increment is cl_exp_component.
    #
    # WHICH ARM IS WHICH. The Results prose says CL1/F is "the initial
    # apparent clearance after single dosing" and CL2/F "the subsequent
    # time-varying apparent clearance at steady state", which reads
    # backwards relative to Eq. 2 and is contradicted by three
    # independent lines of evidence, all of which agree that CL1/F is
    # the TIME-INDEPENDENT asymptote and CL2/F the DECAYING increment:
    #   (a) Eq. 2 itself -- exp(-Kout*(DAY-1)) multiplies TVCL2, so it is
    #       TVCL2 that vanishes with time and TVCL1 that survives.
    #   (b) Figure 1 labels the two elimination arms leaving the central
    #       compartment "Time independent CL1" and "Time dependent CL2".
    #   (c) The printed half-lives close the case arithmetically. On
    #       day 1, exp(0) = 1 so CL/F = 64 + 40.3 = 104.3 L/h and
    #       t1/2 = ln(2)*5210/104.3 = 34.6 h -- the paper's "34.6 h after
    #       the first dose". As DAY grows, CL/F -> 64 L/h and
    #       t1/2 = ln(2)*5210/64 = 56.4 h -- the paper's "approximately
    #       56.4 h" at steady state. The opposite assignment reproduces
    #       neither number.
    # See the vignette Errata; the prose sentence is the error, not the
    # parameter values.
    # ------------------------------------------------------------------
    lcl_exp_inf       <- log(64)   ; label("Time-independent apparent clearance component CL1/F, the steady-state asymptote of CL/F (L/h)")  # Liu 2024 Table 2 theta3 = 64 L/h, SE 1.51, RSE 2.4%; bootstrap 64 (95% CI 61.1-66.9)
    lcl_exp_component <- log(40.3) ; label("Time-dependent (autoinhibited) apparent clearance component CL2/F, the increment present on day 1 that decays away (L/h)")  # Liu 2024 Table 2 theta4 = 40.3 L/h, SE 6.89, RSE 17.1%; bootstrap 40.5 (95% CI 25.5-57.5)
    lcl_exp_kdes      <- log(1.35) ; label("Time-varying inhibition rate constant Kout governing the day-to-day decay of CL2/F (1/day)")     # Liu 2024 Table 2 theta5 = 1.35 1/day, SE 0.421, RSE 31.2%; bootstrap 1.35 (95% CI 0.927-2.24). Table 2 unit label "L d-1" is an erratum; Kout is 1/day

    ltlag <- log(0.483) ; label("Absorption lag time ALAG1 (h)")                                      # Liu 2024 Table 2 theta6 = 0.483 h, SE 0.0306, RSE 6.3%; bootstrap 0.491 (95% CI 0.374-0.562)

    # ------------------------------------------------------------------
    # Covariate effects. Age is a power term on TOTAL CL/F -- confirmed
    # by Figure 4, which quotes +27% at age 25 and -8.8% at age 71
    # against the median 53: (25/53)^-0.314 = 1.266 and
    # (71/53)^-0.314 = 0.912. Applying the exponent to CL1/F alone would
    # give only +16% at age 25, so the multiplier is shared by both
    # clearance arms (equivalently, it scales their sum).
    #
    # The three ALKPOT effects are multiplicative FACTORS on V/F, not
    # log-scale shifts: Table 2 prints 0.734 / 0.917 / 0.784 and the
    # Results narrative reads them straight off as "V/F reductions were
    # 26.6% ... 8.3% ... 21.6%", i.e. 1 - 0.734, 1 - 0.917, 1 - 0.784.
    # They are applied in power form (factor^indicator) so that an
    # indicator of 0 leaves V/F untouched.
    # ------------------------------------------------------------------
    e_age_cl      <- -0.314 ; label("Power exponent of age on total apparent clearance CL/F, referenced to 53 years (unitless)")             # Liu 2024 Table 2 theta7 = -0.314, SE 0.0799, RSE 25.4%; bootstrap -0.314 (95% CI -0.476 to -0.167)
    e_alki_vc     <- 0.734  ; label("Multiplicative effect on V/F of prior ALK-inhibitor treatment in ALK+ patients (ALKPOT 2,3,4) vs ALK-inhibitor-naive (unitless factor)")  # Liu 2024 Table 2 theta8 = 0.734, SE 0.0649, RSE 8.8%; bootstrap 0.724 (95% CI 0.541-0.89)
    e_nonalk_vc   <- 0.917  ; label("Multiplicative effect on V/F for ROS1-positive or unknown patients (ALKPOT 5) vs ALK-inhibitor-naive ALK+ patients (unitless factor)")    # Liu 2024 Table 2 theta9 = 0.917, SE 0.0841, RSE 9.2%; bootstrap 0.905 (95% CI 0.729-1.12)
    e_healthy_vc  <- 0.784  ; label("Multiplicative effect on V/F for healthy volunteers (ALKPOT 6) vs ALK-inhibitor-naive ALK+ patients (unitless factor)")                   # Liu 2024 Table 2 theta10 = 0.784, SE 0.148, RSE 18.9%; bootstrap 0.782 (95% CI 0.557-0.939)

    # ------------------------------------------------------------------
    # Interindividual variability. Table 2 reports each entry as
    # "omega^2" together with a percent CV, and the CV column proves the
    # numbers are log-scale VARIANCES rather than SDs:
    #   sqrt(exp(0.125) - 1) = 36.5%   (printed 36.5%)
    #   sqrt(exp(0.279) - 1) = 56.7%   (printed 56.7%)
    #   sqrt(exp(0.138) - 1) = 38.5%   (printed 38.5%)
    #   sqrt(exp(0.120) - 1) = 35.7%   (printed 35.7%)
    # All four reproduce to the printed decimal, so the variances are
    # entered directly. nlmixr2 expects variances here.
    #
    # A SINGLE eta is carried on clearance. Liu 2024 Results says an
    # exponential IIV model was "implemented to describe the IIVs,
    # including CL1/F, CL2/F, V/F, KA, and tlag", but Table 2 lists only
    # four eta rows and only one of them is a clearance term (omega^2
    # CL). Eq. 2 resolves the apparent conflict: exp(eta_i) multiplies
    # the whole bracket, so one eta is shared by both clearance arms.
    # That is exactly what model() does -- etalcl enters both
    # cl_exp_inf and cl_exp_component, which is algebraically identical
    # to multiplying their sum.
    #
    # Table 2 reports no off-diagonal covariances, so the etas are
    # diagonal. Shrinkage is high on KA (53.3%) and on the lag time
    # (66.2%); see the vignette Assumptions and deviations.
    # ------------------------------------------------------------------
    etalvc   ~ 0.125 ; label("IIV on V/F (log-scale variance)")                        # Liu 2024 Table 2 eta1 omega^2 V = 0.125 (CV 36.5%), SE 0.0216, RSE 17.3%, shrinkage 41.1%
    etalka   ~ 0.279 ; label("IIV on KA (log-scale variance)")                         # Liu 2024 Table 2 eta2 omega^2 KA = 0.279 (CV 56.7%), SE 0.0882, RSE 31.6%, shrinkage 53.3%
    etalcl   ~ 0.138 ; label("IIV on total apparent clearance CL/F, shared by both clearance arms (log-scale variance)")  # Liu 2024 Table 2 eta3 omega^2 CL = 0.138 (CV 38.5%), SE 0.011, RSE 8%, shrinkage 2.2%
    etaltlag ~ 0.120 ; label("IIV on the absorption lag time ALAG1 (log-scale variance)")  # Liu 2024 Table 2 eta4 omega^2 ALAG1 = 0.12 (CV 35.7%), SE 0.0338, RSE 32.3%, shrinkage 66.2%

    # ------------------------------------------------------------------
    # Residual error: combined proportional plus additive (Liu 2024
    # Results, "Residual variability was described using a combined
    # proportional and additive residual error model"). Table 2 prints
    # 0.0468 and 0.0938 as the epsilon entries.
    #
    # THESE ARE VARIANCES, not SDs, so nlmixr2's SD-scaled propSd/addSd
    # take their square roots. Two checks:
    #   (a) Read as a variance, the proportional term is
    #       sqrt(0.0468) = 21.6% -- an ordinary proportional residual
    #       error for pooled sparse phase I/II oncology data. Read as an
    #       SD it would be 4.68%, which is not a credible residual for
    #       this design and would be tighter than the assay itself.
    #   (b) Read as a variance, the additive term is
    #       sqrt(0.0938) = 0.306 ng/mL against a validated assay range of
    #       0.5-150 ng/mL -- the textbook "additive SD near the LLOQ"
    #       result. Read as an SD it would be 0.0938 ng/mL, a fifth of
    #       the LLOQ, which no assay in that range supports.
    # Both terms are therefore square-rooted. This is the one place in
    # the file where the printed number is transformed rather than used
    # verbatim; see the vignette Errata.
    # ------------------------------------------------------------------
    propSd <- 0.216333 ; label("Proportional residual error (fraction)")   # Liu 2024 Table 2 epsilon1 = 0.0468 (variance), SE 0.0008, RSE 1.6%, shrinkage 7.7%; sqrt(0.0468) = 0.216333
    addSd  <- 0.306268 ; label("Additive residual error (ng/mL)")          # Liu 2024 Table 2 epsilon2 = 0.0938 (variance), SE 0.0281, RSE 30%, shrinkage 7.7%; sqrt(0.0938) = 0.306268
  })

  model({
    # ------------------------------------------------------------------
    # 1. Derived covariate terms.
    #    The ALKPOT classification is a single four-cell factor, encoded
    #    with three mutually exclusive indicators against the reference
    #    cell (an ALK-positive, ALK-inhibitor-naive patient):
    #      PRIOR_ALKI = 1                       -> ALKPOT 2, 3, 4
    #      nonalk     = 1                       -> ALKPOT 5
    #      DIS_HEALTHY = 1                      -> ALKPOT 6
    #    Healthy volunteers are not ALK-positive, so the (1 - DIS_HEALTHY)
    #    gate is what keeps them out of the ALKPOT 5 cell.
    # ------------------------------------------------------------------
    age_cl <- (AGE / 53) ^ e_age_cl
    nonalk <- (1 - TUM_ALK_MUT) * (1 - DIS_HEALTHY)

    # ------------------------------------------------------------------
    # 2. Individual parameters.
    # ------------------------------------------------------------------
    vc   <- exp(lvc + etalvc) *
      e_alki_vc    ^ PRIOR_ALKI *
      e_nonalk_vc  ^ nonalk *
      e_healthy_vc ^ DIS_HEALTHY
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    # ------------------------------------------------------------------
    # 3. Time-dependent clearance, Liu 2024 Eq. 2.
    #
    #    DAY is the paper's INTEGER post-administration day (day 1 is the
    #    day of the first dose), and (DAY - 1) is therefore the number of
    #    completed days on treatment. rxode2 carries `time` in hours, so
    #    the exponent is floor(time / 24) -- a STEP function that drops
    #    once per calendar day, which is precisely what the abstract
    #    describes: "time-dependent elimination by allowing the clearance
    #    to decrease STEPWISE over time". A continuous time/24 would be a
    #    different (smooth) model than the one that was fitted.
    #
    #    With Kout = 1.35 /day the increment has a half-life of
    #    ln(2)/1.35 = 0.51 days, so CL/F is within 1% of its asymptote by
    #    about day 4.
    # ------------------------------------------------------------------
    cl_exp_kdes      <- exp(lcl_exp_kdes)
    cl_exp_inf       <- exp(lcl_exp_inf       + etalcl) * age_cl
    cl_exp_component <- exp(lcl_exp_component + etalcl) * age_cl
    cl               <- cl_exp_inf + cl_exp_component * exp(-cl_exp_kdes * floor(time / 24))

    # ------------------------------------------------------------------
    # 4. One-compartment ODE system with first-order oral absorption.
    # ------------------------------------------------------------------
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ------------------------------------------------------------------
    # 5. Absorption lag. Liu 2024 selected a lag-time model over transit
    #    compartments: "the pcVPC plots indicated that the transit
    #    compartment models did not offer significant improvements over
    #    the lag-time model."
    # ------------------------------------------------------------------
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 6. Observation. Dose is in mg and vc in L, so central / vc is mg/L;
    #    multiply by 1000 to report ng/mL, the units of the validated
    #    assay (range 0.5-150 ng/mL) and of every exposure metric in the
    #    companion exposure-response models. The scale checks out against
    #    the paper's own steady-state exposure: with CL/F = 64 L/h the
    #    160 mg once-daily AUCss is 160/64 = 2.5 mg*h/L = 2,500 ng*h/mL,
    #    against the reported geometric mean of 2,374 ng*h/mL.
    # ------------------------------------------------------------------
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
