Hartmann_2026_nintedanib <- function() {
  description <- paste0(
    "Pediatric population pharmacokinetic model for nintedanib in children ",
    "and adolescents 6 to less than 18 years of age with clinically ",
    "significant fibrosing interstitial lung disease (Hartmann 2026; the ",
    "phase 3 InPedILD trial and its open-label extension InPedILD-ON). ",
    "One-compartment model with an absorption lag time, first-order ",
    "absorption and first-order elimination from the central compartment. ",
    "Allometric body weight is carried on apparent clearance and apparent ",
    "volume with the standard fixed exponents 0.75 and 1 at a 75 kg ",
    "reference. Relative bioavailability carries three covariates - ",
    "ethnicity (Korean, and a composite of all non-White non-Korean ",
    "groups), systemic-sclerosis-associated ILD, and baseline lactate ",
    "dehydrogenase - plus inter-individual variability and six-occasion ",
    "inter-occasion variability whose magnitude increases with decreasing ",
    "pediatric age. Residual error is additive on the log-transformed ",
    "concentration scale. The model was estimated in NONMEM with the ",
    "NWPRI frequentist-prior functionality using the adult nintedanib ",
    "popPK meta-model as prior; every parameter except the pediatric age ",
    "effect on inter-occasion variability and the residual error was ",
    "supported by that adult prior. A fixed relative-bioavailability ",
    "multiplier for mild (Child-Pugh class A) hepatic impairment is ",
    "included so the paper's pediatric hepatic-impairment dose-adjustment ",
    "simulations can be reproduced."
  )
  reference <- paste(
    "Hartmann S, Chan Kwong A, Ribbing J, Gahlemann M, Korell J.",
    "Population Pharmacokinetics and Exposure-Response Model-Based",
    "Bayesian Extrapolation of FVC-Based Efficacy Endpoints From Adults to",
    "Pediatric Patients Receiving Nintedanib.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70135.",
    "doi:10.1002/psp4.70135. PMCID PMC12823301.",
    "Structural parameters, covariate coefficients and variance terms are",
    "the final estimates in Table S5 of Data S2 and in the final pediatric",
    "popPK NONMEM control stream reproduced in Data S1.",
    "The adult nintedanib population PK model that supplied the",
    "frequentist prior is a separate publication; a different adult",
    "nintedanib popPK analysis is packaged as",
    "modellib('Schmid_2017_nintedanib').",
    sep = " "
  )
  vignette <- "Hartmann_2026_nintedanib"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "nintedanib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nintedanib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight. Allometric power scaling on apparent clearance and apparent volume of distribution.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: InPedILD assigned the starting dose from body weight",
        "at the baseline visit and re-assigned it during treatment as body",
        "weight changed (Hartmann 2026 Table S1), and the control stream",
        "carries both a baseline WTKGBL and a time-varying WTKG column,",
        "using the time-varying WTKG in the $PK block.",
        "Reference 75 kg, the centring constant hard-coded in the final",
        "control stream (Data S1, WTCL = (WTKG/75)**0.75 and",
        "WTV = (WTKG/75)**1); it is an adult-model constant inherited with",
        "the prior, NOT the pediatric cohort median, which is 42.9 kg",
        "(Hartmann 2026 Table S4).",
        "Exponents are FIXED at 0.75 (CL/F) and 1 (V/F) in line with the",
        "adult prior (Hartmann 2026 Table S2 footnote a and Table S5",
        "footnotes a and b).",
        "Pediatric PK cohort mean 42.9 kg (SD 17.8); 27.2 kg (SD 10.3) in",
        "the 6 to less than 12 year group and 50.3 kg (SD 15.7) in the 12",
        "to less than 18 year group. Patients below 13.5 kg were excluded",
        "from the trial (Table S1 footnote)."
      ),
      source_name        = "WTKG (time-varying); WTKGBL is the baseline value used only for stratification"
    ),
    AGE = list(
      description        = "Subject age. Scales the magnitude of the inter-occasion variability on relative bioavailability, with no effect at or above 18 years.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TIME-VARYING with an upper cut-off at 18 years, exactly as tested",
        "in Hartmann 2026 Table S2 footnote b, option 1 (continuous with",
        "cut-off at 18 years of age). The final control stream (Data S1)",
        "implements it as",
        "IF(AGEY.GE.18) IOVPEDAGE = 1 ELSE IOVPEDAGE = EXP(THETA(9)*(AGEY-18)),",
        "so the IOV scale factor is 1 for an adult and rises as pediatric",
        "age falls. The model file reproduces the two branches as the",
        "single continuous expression age_iov = min(AGE, 18) so no",
        "branching is needed.",
        "This is the ONLY covariate in the pediatric popPK model that was",
        "not supported by the adult prior (Hartmann 2026 Table S5, the one",
        "row without an asterisk besides RUV).",
        "Pediatric PK cohort mean 13.2 years (SD 3.08); 9.57 (SD 1.72) in",
        "the 6 to less than 12 year group and 14.9 (SD 1.83) in the 12 to",
        "less than 18 year group (Hartmann 2026 Table S4)."
      ),
      source_name        = "AGEY (time-varying age in years); AGEYBL is the baseline value used only for stratification"
    ),
    OCC = list(
      description        = "Integer occasion index, 1 to 6, indexing the inter-occasion variability on relative bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "The final control stream (Data S1) assigns a separate eta per",
        "occasion for OCC 1 through 6",
        "(IF(OCC.EQ.1) IOVF=ETA(4)*TVETIOF through",
        "IF(OCC.EQ.6) IOVF=ETA(9)*TVETIOF), backed by",
        "$OMEGA BLOCK(1) followed by five $OMEGA BLOCK(1) SAME blocks, so",
        "all six occasions share one variance. Six is therefore stated by",
        "the source rather than inferred.",
        "An occasion is a PK sampling visit: InPedILD sampled at the week 2",
        "and week 26 visits, and InPedILD-ON at week 2 or week 12 and at",
        "week 24 (Hartmann 2026 Section 2.3.1).",
        "Records outside any sampling occasion take OCC = 0, which zeroes",
        "every indicator and leaves relative bioavailability with",
        "inter-individual variability only."
      ),
      source_name        = "OCC"
    ),
    LDH = list(
      description        = "Baseline serum lactate dehydrogenase. Exponential effect on relative bioavailability, centred at 206 U/L.",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline value. Effect on relative bioavailability:",
        "F1 = ... * exp(e_ldh_fdepot * (LDH - 206)), from the final control",
        "stream (Data S1) F1LDHBL = EXP(THETA(8)*(LDHBL - 206)).",
        "The 206 U/L centring constant is hard-coded in the control stream",
        "and is inherited from the adult model with the prior; it is NOT",
        "the pediatric cohort mean, which is 244 U/L (SD 83.9;",
        "Hartmann 2026 Table S4). At the pediatric mean the multiplier is",
        "therefore exp(0.00155711 * 38) = 1.06.",
        "The control stream comment states there are no missing LDH values",
        "in the analysis data set.",
        "Pediatric cohort mean 276 U/L (SD 105) in the 6 to less than 12",
        "year group and 229 U/L (SD 69.1) in the 12 to less than 18 year",
        "group (Hartmann 2026 Table S4).",
        "This is a mechanistic (pre-specified) covariate, included in the",
        "starting model rather than selected by the stepwise covariate",
        "search (Hartmann 2026 Table S2)."
      ),
      source_name        = "LDHBL"
    ),
    RACE_WHITE = list(
      description        = "White / Caucasian race indicator (1 = White, 0 = otherwise). The reference category of the three-level ethnicity effect on relative bioavailability.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (White) is the reference level itself; the multiplier on relative bioavailability is exactly 1 when RACE_WHITE = 1.",
      notes              = paste(
        "Time-fixed baseline. The final control stream (Data S1) encodes a",
        "three-level ethnicity effect on relative bioavailability:",
        "IF(RACEREG1.EQ.1) F1RACEREG = 1 (White);",
        "IF(RACEREG1.EQ.0.AND.RACEREG3.NE.1) F1RACEREG = (1 + THETA(5))",
        "(not White and not Korean, labelled Other);",
        "IF(RACEREG3.EQ.1) F1RACEREG = (1 + THETA(6)) (Korean).",
        "RACE_WHITE carries RACEREG1 and RACE_KOREAN carries RACEREG3; the",
        "Other group is DERIVED inside model() as",
        "(1 - RACE_WHITE) * (1 - RACE_KOREAN) so the three levels stay",
        "mutually exclusive and exhaustive.",
        "Hartmann 2026 Table S5 names the Other group as",
        "Chinese/Taiwanese/Indian/Japanese/Other Asian/Black/American",
        "Indian/Alaska Native (footnote c).",
        "This is a mechanistic (pre-specified) covariate (Table S2).",
        "The pediatric PK cohort is 77% Caucasian, 4.5% Other Asian, 9.1%",
        "Black, 6.8% American Indian/Alaska Native and 2.3% missing",
        "(Hartmann 2026 Table S4), so NO pediatric patient was Korean and",
        "the Korean coefficient is carried entirely by the adult prior.",
        "A missing-ethnicity subject was assigned to the Other group by the",
        "control stream branch, which tests only RACEREG1 and RACEREG3."
      ),
      source_name        = "RACEREG1"
    ),
    RACE_KOREAN = list(
      description        = "Korean-heritage race indicator (1 = Korean, 0 = otherwise). Third level of the ethnicity effect on relative bioavailability.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White, or any non-Korean group; the White level is the multiplicative reference).",
      notes              = paste(
        "Time-fixed baseline. Effect on relative bioavailability:",
        "F1 = ... * (1 + e_korean_fdepot) when RACE_KOREAN = 1, from the",
        "final control stream (Data S1)",
        "IF(RACEREG3.EQ.1) F1RACEREG = (1 + THETA(6)).",
        "NO pediatric patient in InPedILD or InPedILD-ON was Korean",
        "(Hartmann 2026 Table S4 ethnicity breakdown), so this coefficient",
        "is estimated entirely from the adult prior and is carried for",
        "structural fidelity and for adult-versus-pediatric comparison.",
        "The same drug and the same Korean-versus-reference contrast on",
        "relative bioavailability appear in the adult nintedanib model",
        "modellib('Schmid_2017_nintedanib'), where the coefficient is",
        "parameterised multiplicatively (0.781) rather than as the",
        "fractional change used here (-0.144, i.e. a multiplier of 0.856)."
      ),
      source_name        = "RACEREG3"
    ),
    DIS_SSC_ILD = list(
      description        = "Systemic-sclerosis-associated interstitial lung disease indicator (1 = SSc-ILD, 0 = fibrosing ILD of any other aetiology).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-SSc-ILD; the most common group, 82% of the pediatric PK cohort).",
      notes              = paste(
        "Time-fixed baseline. Effect on relative bioavailability:",
        "F1 = ... * (1 + e_sscild_fdepot) when DIS_SSC_ILD = 1, from the",
        "final control stream (Data S1)",
        "IF(SSCSUB4.EQ.0) F1DIAG6 = 1 (labelled Most common) and",
        "IF(SSCSUB4.EQ.1) F1DIAG6 = (1 + THETA(7)) (labelled Ssc ILD).",
        "This is a mechanistic (pre-specified) covariate (Hartmann 2026",
        "Table S2), inherited from the adult model where the SENSCIS",
        "SSc-ILD population contributed a large share of the data.",
        "In the pediatric PK cohort 8 of 44 patients (18%) had SSc-ILD",
        "(Hartmann 2026 Table S4); in the larger FVC cohort 9 of 53 (17%)",
        "did (Table 1).",
        "Note that SSc-ILD is a systemic-sclerosis diagnosis and is",
        "recorded separately from the ILD-aetiology categories in Table 1",
        "and Table S4, so a patient can carry both an ILD diagnosis",
        "category and DIS_SSC_ILD = 1."
      ),
      source_name        = "SSCSUB4"
    ),
    HEPIMP_MILD = list(
      description        = "Mild (Child-Pugh class A) hepatic impairment indicator (1 = Child-Pugh A, 0 = no hepatic impairment). Multiplicative effect on relative bioavailability.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hepatic impairment; every patient in the model-development data set).",
      notes              = paste(
        "NOT a fitted covariate. No patient contributing to the pediatric",
        "popPK model had hepatic impairment (Hartmann 2026 Discussion:",
        "none of the patients included in the development of the popPK",
        "model had hepatic impairment). The multiplier is an ASSUMPTION",
        "the paper states outright and then uses to generate its pediatric",
        "hepatic-impairment dose recommendation, so it is carried here as",
        "a fixed() coefficient rather than dropped.",
        "Hartmann 2026 Section 2.3.3: For children with Child-Pugh class A,",
        "a 115% higher nintedanib bioavailability was assumed compared to",
        "children without hepatic impairment, based on the estimated",
        "difference in exposure seen in adult patient with and without",
        "Child-Pugh class A. A 115% increase is a multiplier of 2.15.",
        "The underlying adult estimate comes from the paper's reference",
        "16, a separate adult hepatic-impairment trial, not from this",
        "analysis.",
        "Set HEPIMP_MILD = 0 to reproduce every result in the paper other",
        "than Figures S5 and S6. The extrapolation was deliberately",
        "limited to Child-Pugh A: nintedanib is not recommended in",
        "Child-Pugh B or C, so there is no Child-Pugh B or C coefficient",
        "to encode."
      ),
      source_name        = "Child-Pugh class A (simulation scenario flag; no source data column)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 2L,
    n_observations = "446 nintedanib plasma concentrations",
    age_range      = "6 to less than 18 years (trial eligibility); cohort mean 13.2 years (SD 3.08)",
    age_median     = "mean 13.2 years (SD 3.08); no median reported",
    weight_range   = "13.5 kg lower eligibility bound; cohort mean 42.9 kg (SD 17.8)",
    weight_median  = "mean 42.9 kg (SD 17.8); no median reported",
    sex_female_pct = 54.5,
    race_ethnicity = c(Caucasian = 77, Black = 9.1, AmericanIndianAlaskaNative = 6.8, OtherAsian = 4.5, Missing = 2.3),
    disease_state  = "Clinically significant fibrosing interstitial lung disease of mixed aetiology: pediatric autoimmune ILD 32%, surfactant protein deficiency 30%, other fibrosing ILD 23%, toxic/radiation/drug-induced pneumonitis 9.1%, chronic hypersensitivity pneumonia 4.5%, post-HSCT fibrosis 2.3%. Systemic sclerosis-associated ILD in 18%.",
    dose_range     = "Oral nintedanib twice daily, dosed by body-weight bin: 50 mg BID for 13.5 to less than 23.0 kg, 75 mg BID for 23.0 to less than 33.5 kg, 100 mg BID for 33.5 to less than 57.5 kg and 150 mg BID at 57.5 kg and above (Hartmann 2026 Table S1). Dose reduction by one step, interruption, re-escalation and discontinuation were permitted.",
    regions        = "Multinational; the InPedILD trial and the InPedILD-ON open-label extension",
    age_group_breakdown = c(`6_to_lt12_years_n` = 14, `12_to_lt18_years_n` = 30),
    formulation_breakdown = c(`25mg_soft_capsule_pct` = 32, `100mg_soft_capsule_pct` = 55, `150mg_soft_capsule_pct` = 14),
    hepatic_function = "No patient contributing to the popPK model had hepatic impairment.",
    notes          = paste0(
      "Pediatric popPK analysis set: 44 of the 53 patients enrolled in ",
      "InPedILD (phase 3, randomised 2:1 nintedanib versus placebo over ",
      "24 weeks) and its open-label extension InPedILD-ON (to the interim ",
      "database lock at 52 weeks). The 9 remaining patients contributed ",
      "only to the exposure-response analyses. Baseline characteristics ",
      "are Hartmann 2026 Table S4 for this PK set and Table 1 for the ",
      "53-patient FVC set. Nintedanib was measured by validated HPLC-MS/MS ",
      "with a lower limit of quantification of 0.415 ng/mL in InPedILD and ",
      "0.5 ng/mL (approximately 0.927 nM) in InPedILD-ON; concentrations ",
      "below the limit of quantification were omitted from the analysis. ",
      "Because the pediatric data set is small, the model was estimated ",
      "with the adult nintedanib popPK meta-model as a frequentist prior ",
      "through the NONMEM NWPRI functionality, after an external ",
      "validation step in which the pre-specified adult model was checked ",
      "against the pediatric data by visual predictive check."
    )
  )

  # Implementation notes (see the vignette section 'Assumptions and
  # deviations' for the full justification of each item):
  # * Amount units are mg and volume units are L, so central / vc is in
  #   mg/L; the observation multiplies by 1853 nM per mg/L to reach the
  #   nM scale the paper reports. Nintedanib free-base molecular weight
  #   539.62 g/mol gives 1 mg = 1.853 umol = 1853 nmol, the same
  #   conversion factor the final control stream applies as the literal
  #   1853 in its steady-state PK-metric block and the same factor used
  #   by modellib('Schmid_2017_nintedanib'). Cross-check: at the adult
  #   reference weight of 71.5 kg this model gives
  #   AUCtau,ss = 1853 * 150 / (909.629 * (71.5/75)^0.75) = 317 nM*h and
  #   Cav,ss = 26.4 nM against the adult reference geometric means of
  #   316 nM*h and 26 nM printed in the Figure 1 caption.
  # * Parameter values are taken from the final pediatric popPK NONMEM
  #   control stream in Data S1, which carries MORE significant digits
  #   than Hartmann 2026 Table S5 and rounds to Table S5 exactly on
  #   every one of the nine thetas and five variance terms. The control
  #   stream $THETA block is therefore the FINAL estimates, not initial
  #   estimates; the parallel FVC control streams in the same supplement
  #   still carry their initial values (identical to the $THETAP prior
  #   block), which is how the two cases are told apart. Table S5 values
  #   are quoted alongside each line for cross-reference.
  # * Table S5 reports every variance term in a column headed CV. Those
  #   numbers are the SQUARE ROOTS of the control-stream $OMEGA and
  #   $SIGMA variances (0.309 = sqrt(0.0952477), 1.25 = sqrt(1.56911),
  #   0.417 = sqrt(0.173479), 0.310 = sqrt(0.0961255) and
  #   0.392 = sqrt(0.153471)), i.e. standard deviations on the log
  #   scale, NOT log(1 + CV^2) coefficients of variation. The variances
  #   below are the control-stream values used directly.
  # * Inter-occasion variability on relative bioavailability is encoded
  #   with the registered OCC idiom: six per-occasion etas built from
  #   mutually exclusive indicators, with occasions 2 to 6 declared
  #   ~ fixed() at the occasion-1 variance to reproduce NONMEM
  #   $OMEGA BLOCK(1) SAME. The whole IOV term is then scaled by the
  #   pediatric-age factor exp(e_age_iov_fdepot * (min(AGE, 18) - 18)),
  #   which reproduces the control stream's IOVPEDAGE two-branch IF
  #   without branching. Because e_age_iov_fdepot is negative and
  #   min(AGE, 18) - 18 is negative or zero, the factor is at least 1
  #   and grows as pediatric age falls, matching the paper's finding
  #   that IOV in relative bioavailability is higher in younger
  #   children.
  # * Residual error is additive on the log-transformed concentration
  #   (control stream $ERROR: IPRED = LOG(F + DEL); Y = IPRED + EPS(1)),
  #   which maps in nlmixr2 to the lognormal residual lnorm(expSd) with
  #   expSd the standard deviation on the natural-log scale.
  # * The typical relative bioavailability is fixed at 1 (control stream
  #   TVF1 = 1 * covariate factors), so CL and V are apparent values
  #   CL/F and V/F. There is no absolute-bioavailability information in
  #   this analysis.
  # * There is NO inter-individual variability on apparent clearance.
  #   That is the model as fitted, and it has a useful consequence for
  #   validation: since the relative-bioavailability eta is lognormal
  #   with a geometric mean of 1, the population geometric-mean
  #   AUCtau,ss equals the typical-value AUCtau,ss exactly.
  # * The Child-Pugh class A multiplier on relative bioavailability is a
  #   fixed() ASSUMPTION carried from the paper's simulation section,
  #   not an estimate from these data; see the HEPIMP_MILD covariate
  #   notes. Leave HEPIMP_MILD at 0 to reproduce the main results.
  ini({
    # ----- Structural parameters (Hartmann 2026 Table S5 / Data S1 final control stream $THETA) -----
    lcl   <- log(909.629);  label("Log nintedanib apparent clearance CL/F at 75 kg (L/h)")          # Data S1 $THETA 1 = 909.629; Table S5 row 'Apparent clearance (CL/F)' = 910 L/h (RSE 2.31%), supported by the adult prior
    lvc   <- log(10696.3);  label("Log nintedanib apparent central volume V/F at 75 kg (L)")        # Data S1 $THETA 2 = 10696.3; Table S5 row 'Apparent volume (V/F)' = 1.07E+04 L (RSE 4.26%), supported by the adult prior
    lka   <- log(2.73232);  label("Log nintedanib first-order absorption rate constant ka (1/h)")   # Data S1 $THETA 3 = 2.73232; Table S5 row 'Absorption rate (ka)' = 2.73 /h (RSE 11.0%), supported by the adult prior
    ltlag <- log(0.71714);  label("Log nintedanib absorption lag time (h)")                         # Data S1 $THETA 4 = 0.71714; Table S5 row 'Lag time (tlag)' = 0.717 h (RSE 3.17%), supported by the adult prior
    lfdepot <- fixed(log(1)); label("Log nintedanib typical relative bioavailability Frel (unitless; reference value 1)")  # Data S1 $PK: TVF1 = 1 * F1RACEREG * F1DIAG6 * F1LDHBL, so the typical value is fixed at 1 and CL and V are apparent

    # ----- Allometric body-weight exponents (fixed in line with the adult prior) -----
    e_wt_cl <- fixed(0.75); label("Allometric exponent of WT/75 on apparent clearance (unitless)")            # Data S1 $PK: WTCL = (WTKG/75)**0.75; Table S5 footnote a and Table S2 footnote a state the exponents are fixed in line with prior information from adults
    e_wt_vc <- fixed(1);    label("Allometric exponent of WT/75 on apparent volume of distribution (unitless)")  # Data S1 $PK: WTV = (WTKG/75)**1; Table S5 footnote b and Table S2 footnote a

    # ----- Covariate effects on relative bioavailability (Hartmann 2026 Table S5) -----
    e_other_fdepot  <- 0.3289;      label("Fractional change in Frel for the non-White non-Korean ethnicity group (unitless)")  # Data S1 $THETA 5 = 0.3289; Table S5 row 'Other ethnicities on Frel' = 0.329 fraction change (RSE 13.3%), supported by the adult prior. Table S5 footnote c: Chinese/Taiwanese/Indian/Japanese/Other Asian/Black/American Indian/Alaska Native
    e_korean_fdepot <- -0.143797;   label("Fractional change in Frel for Korean ethnicity (unitless)")                          # Data S1 $THETA 6 = -0.143797; Table S5 row 'Korean on Frel' = -0.144 fraction change (RSE 36.9%), supported by the adult prior. No pediatric patient was Korean
    e_sscild_fdepot <- -0.138784;   label("Fractional change in Frel for systemic-sclerosis-associated ILD (unitless)")         # Data S1 $THETA 7 = -0.138784; Table S5 row 'SSc-ILD on Frel' = -0.139 fraction change (RSE 27.4%), supported by the adult prior
    e_ldh_fdepot    <- 0.00155711;  label("Exponential coefficient on (LDH - 206) for Frel (per U/L)")                          # Data S1 $THETA 8 = 0.00155711 with F1LDHBL = EXP(THETA(8)*(LDHBL - 206)); Table S5 row 'LDH on Frel' = 0.00156 (RSE 17.4%), supported by the adult prior

    # ----- Pediatric age effect on the magnitude of inter-occasion variability -----
    e_age_iov_fdepot <- -0.0525537; label("Exponential coefficient on (min(AGE, 18) - 18) scaling the IOV on Frel (per year)")  # Data S1 $THETA 9 = -0.0525537 with IOVPEDAGE = EXP(THETA(9)*(AGEY-18)) below 18 years and 1 at or above; Table S5 row 'Paediatric age on IOV Frel' = -0.0526 change per years of age (RSE 43.7%). This row carries NO asterisk in Table S5: it is one of only two parameters NOT supported by the adult prior

    # ----- Mild hepatic impairment (assumption carried from the paper's simulations, not an estimate) -----
    e_hepimp_fdepot <- fixed(2.15); label("Multiplicative effect on Frel for Child-Pugh class A hepatic impairment (unitless)")  # Hartmann 2026 Section 2.3.3: a 115% higher nintedanib bioavailability was assumed for Child-Pugh class A, i.e. a multiplier of 2.15, carried from the adult hepatic-impairment trial cited as reference 16. No patient in the model-development data set had hepatic impairment, so this is fixed, not estimated

    # ----- Inter-individual variability (Hartmann 2026 Table S5 / Data S1 $OMEGA) -----
    # Table S5 reports the SQUARE ROOT of each variance in a column headed
    # CV; the variances below are the control-stream values.
    etalvc      ~ 0.0952477   # Data S1 $OMEGA 1 = 0.0952477; Table S5 row 'IIV V/F' CV = 0.309 = sqrt(0.0952477) (RSE 10.1%, shrinkage 27.1%), supported by the adult prior
    etalka      ~ 1.56911     # Data S1 $OMEGA 2 = 1.56911; Table S5 row 'IIV ka' CV = 1.25 = sqrt(1.56911) (RSE 8.76%, shrinkage 32%), supported by the adult prior
    etalfdepot  ~ 0.173479    # Data S1 $OMEGA 3 = 0.173479; Table S5 row 'IIV Frel' CV = 0.417 = sqrt(0.173479) (RSE 4.01%, shrinkage 13.6%), supported by the adult prior

    # ----- Inter-occasion variability on relative bioavailability, 6 occasions -----
    # Data S1 declares $OMEGA BLOCK(1) 0.0961255 followed by five
    # $OMEGA BLOCK(1) SAME blocks, so all six occasions share one
    # variance. Occasions 2 to 6 are fixed at the occasion-1 value to
    # reproduce SAME.
    etaiov_fdepot_1 ~ 0.0961255           # Data S1 $OMEGA BLOCK(1) = 0.0961255; Table S5 row 'IOV Frel' CV = 0.310 = sqrt(0.0961255) (RSE 3.27%, shrinkage 14.9%), supported by the adult prior
    etaiov_fdepot_2 ~ fixed(0.0961255)    # Data S1 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fixed(0.0961255)    # Data S1 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fixed(0.0961255)    # Data S1 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_5 ~ fixed(0.0961255)    # Data S1 $OMEGA BLOCK(1) SAME
    etaiov_fdepot_6 ~ fixed(0.0961255)    # Data S1 $OMEGA BLOCK(1) SAME

    # ----- Residual unexplained variability -----
    expSd <- 0.391754; label("Lognormal residual SD on nintedanib plasma concentration (log[nM])")  # Data S1 $SIGMA = 0.153471 with the comment 'add err on log-scale'; sqrt(0.153471) = 0.391754. Table S5 row 'RUV' CV = 0.392 (RSE 4.26%, shrinkage 14.8%). This row carries NO asterisk: RUV was estimated for the pediatric patients without support of the adult prior
  })
  model({
    # ----- Reference covariate values (hard-coded in the final control stream) -----
    ref_wt  <- 75      # kg; Data S1 WTCL = (WTKG/75)**0.75 and WTV = (WTKG/75)**1
    ref_ldh <- 206     # U/L; Data S1 F1LDHBL = EXP(THETA(8)*(LDHBL - 206))
    ref_age <- 18      # years; Data S1 IOVPEDAGE branches at AGEY = 18

    # ----- Concentration unit conversion -----
    # central is in mg and vc in L, so central / vc is mg/L. Nintedanib
    # free-base MW = 539.62 g/mol, so 1 mg/L = 1853 nmol/L = 1853 nM,
    # the same literal the control stream uses in its PK-metric block.
    cf_nint <- 1853

    # ----- Derived ethnicity group -----
    # Three mutually exclusive levels: White is the reference with a
    # multiplier of exactly 1; Korean is its own level; everyone else
    # (including a subject with missing ethnicity, which the control
    # stream's IF branches route here) falls into the Other group.
    race_other <- (1 - RACE_WHITE) * (1 - RACE_KOREAN)

    # ----- Occasion indicators for the IOV on relative bioavailability -----
    # Records outside a sampling occasion carry OCC = 0, which zeroes
    # every indicator and leaves only the IIV on relative bioavailability.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    occ5 <- (OCC == 5)
    occ6 <- (OCC == 6)

    # ----- Pediatric-age scaling of the IOV magnitude -----
    # Control stream: IF(AGEY.GE.18) IOVPEDAGE = 1 ELSE
    # IOVPEDAGE = EXP(THETA(9)*(AGEY-18)). Capping age at 18 collapses
    # both branches into one expression, because at AGE >= 18 the
    # exponent is exactly 0 and the factor is exactly 1.
    age_iov   <- AGE * (AGE < ref_age) + ref_age * (AGE >= ref_age)
    iov_scale <- exp(e_age_iov_fdepot * (age_iov - ref_age))

    iov_fdepot <- iov_scale *
      (occ1 * etaiov_fdepot_1 + occ2 * etaiov_fdepot_2 + occ3 * etaiov_fdepot_3 +
       occ4 * etaiov_fdepot_4 + occ5 * etaiov_fdepot_5 + occ6 * etaiov_fdepot_6)

    # ----- Individual pharmacokinetic parameters -----
    # Allometric weight on CL/F and V/F with the fixed adult exponents.
    # There is no IIV on CL/F in this model.
    cl   <- exp(lcl) * (WT / ref_wt)^e_wt_cl
    vc   <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    # ----- Relative bioavailability with covariates, IIV and IOV -----
    # Control stream:
    #   TVF1 = 1 * F1RACEREG * F1DIAG6 * F1LDHBL
    #   F1A  = TVF1 * EXP(ETA(3))
    #   F1   = F1A * EXP(IOVF)
    # The Child-Pugh class A multiplier is appended last; it is 1 unless
    # the user sets HEPIMP_MILD = 1 to reproduce Figures S5 and S6.
    f_race   <- 1 + e_other_fdepot * race_other + e_korean_fdepot * RACE_KOREAN
    f_sscild <- 1 + e_sscild_fdepot * DIS_SSC_ILD
    f_ldh    <- exp(e_ldh_fdepot * (LDH - ref_ldh))

    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) *
      f_race * f_sscild * f_ldh *
      e_hepimp_fdepot^HEPIMP_MILD

    # ----- ODE system -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    # ----- Bioavailability and absorption lag time -----
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ----- Observation in nM -----
    Cc <- (central / vc) * cf_nint

    # ----- Residual error: additive on the log scale, i.e. lognormal -----
    Cc ~ lnorm(expSd)
  })
}
