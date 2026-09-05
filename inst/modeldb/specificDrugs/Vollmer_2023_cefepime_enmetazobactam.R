Vollmer_2023_cefepime_enmetazobactam <- function() {
  description <- paste(
    "Simultaneous four-compartment (two compartments per drug) population PK",
    "model for the fixed-ratio (4:1) cefepime-enmetazobactam combination",
    "(EXBLIFEP), fitted jointly to 5,070 cefepime plasma samples from 588",
    "subjects and 6,342 enmetazobactam plasma samples from 649 subjects",
    "pooled across three Phase 1 studies, one Phase 2 study and the Phase 3",
    "ALLIUM trial in adults with complicated urinary tract infection or acute",
    "pyelonephritis (Vollmer 2023). Both drugs are given as zero-order",
    "intravenous infusions with linear clearance from their own central",
    "compartment. De-indexed (absolute) eGFR is a power covariate on the",
    "clearance of both drugs with an exponent near 1, consistent with",
    "glomerular filtration being the dominant elimination route; age and body",
    "weight are power covariates on both central volumes; enmetazobactam",
    "additionally carries a cUTI-infection effect on its central volume and",
    "sex plus de-indexed eGFR effects on its peripheral volume. Both",
    "inter-compartmental clearances are FIXED to the Phase 1+2 analysis",
    "because the sparse Phase 3 sampling could not characterise distribution.",
    "Inter-individual variability is stratified by infection status (healthy",
    "volunteer versus infected patient) and is markedly larger in patients;",
    "a four-way OMEGA block across cefepime and enmetazobactam CL and Vc in",
    "infected subjects carries cross-drug correlations of 0.93 (CL) and 0.96",
    "(Vc). Cefepime is the unsuffixed parent; enmetazobactam carries the",
    "sibling-drug suffix _enm throughout."
  )
  reference <- paste(
    "Vollmer J, Belley A, Velicitat P, Machacek M. 2529. Population",
    "pharmacokinetic models for cefepime and enmetazobactam derived from",
    "pooled Phase 1 to Phase 3 clinical studies. Open Forum Infect Dis.",
    "2023;10(Suppl 2):ofad500.2147 (IDWeek 2023, Session 241, abstract 2529).",
    "doi:10.1093/ofid/ofad500.2147.",
    "The IDWeek abstract publishes the typical values for a 70 kg cUTI",
    "patient and their %CV but not the covariate coefficients, correlations",
    "or residual error, so the complete final-model parameter set encoded",
    "here is transcribed from the two regulatory reproductions of the same",
    "Allecra population PK study report AAI101-PK-21-01: US FDA NDA 216165",
    "(EXBLIFEP) Integrated Review, 22 February 2024, Table 101 'Parameter",
    "Estimates (RSE) and Median (%CV) for the Applicant's Final Model'",
    "(sponsor page 135, reproducing study-report Table 12); and EMA EPAR",
    "EMA/63929/2024 for Exblifep, Tables 5 and 6 (assessment-report pages",
    "47-48, reproducing study-report Tables 19 and 20). The EMA report is",
    "also the only source that states the reference individual: 'a body",
    "weight of 70 kg, de-indexed eGFR of 100 mL/min, and an age of 50 years'",
    "(page 46), with Sex = Female and Infection = healthy for enmetazobactam",
    "(Table 6 note).",
    sep = " "
  )
  vignette <- "Vollmer_2023_cefepime_enmetazobactam"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Inter-individual variability is stratified by infection status, so the
  # canonical etalcl / etalvc slots consumed by the individual-parameter
  # expressions are multiplexed inside model() from stratum-specific ini()
  # magnitudes. Same construction as the study-phase-stratified residual SDs
  # in Xie_2025_aztreonam_avibactam and Cammarata_2024_sulbactam_durlobactam,
  # applied to the random effects instead of the residuals. Cefepime Vc has
  # no healthy-stratum eta (FDA Table 101 prints "-" for both the estimate
  # and its RSE), so only three of the four multiplexed slots have two arms.
  paper_specific_etas <- c(
    "etalcl_inf", "etalcl_hlth", "etalvc_inf",
    "etalcl_enm_inf", "etalcl_enm_hlth", "etalvc_enm_inf", "etalvc_enm_hlth"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Doses are in mg and volumes in L, so concentrations
  # are mg/L (= ug/mL, the unit the FDA label reports).
  compartmentData <- list(
    central         = list(analyte = "cefepime",        units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "cefepime",        units = "mg", specimen = "plasma", verified = TRUE),
    central_enm     = list(analyte = "enmetazobactam",  units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_enm = list(analyte = "enmetazobactam",  units = "mg", specimen = "plasma", verified = TRUE)
  )

  dosing <- c("central", "central_enm")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate on the CENTRAL volume of distribution of both",
        "drugs, normalised to the 70 kg reference individual named in EMA",
        "EPAR EMA/63929/2024 page 46. Exponents 0.802 (cefepime) and 0.618",
        "(enmetazobactam) - FDA Table 101 rows 'Body weight on FEP V1' and",
        "'Body weight on ENM V1'. NOTE that the EMA assessment-report PROSE",
        "on page 48 attributes 0.618 to cefepime and 0.802 to",
        "enmetazobactam, which contradicts its OWN Tables 5 and 6 on the two",
        "preceding pages as well as FDA Table 101; the two independent",
        "tables agree with each other and are used here. The 70 kg reference",
        "is independently confirmed by the IDWeek abstract, whose",
        "'typical values for a 70 kg cUTI patient' table prints cefepime Vc",
        "as 11.2 L - exactly the untransformed estimate, which can only hold",
        "if the weight reference is 70 kg. Neither drug carries a weight",
        "effect on clearance or on the peripheral volume. Time-fixed per",
        "subject; analysis-set mean 76.52 kg, range 45-135 kg (FDA Table",
        "100)."
      ),
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power covariate on the CENTRAL volume of distribution of both",
        "drugs, normalised to the 50-year reference individual (EMA EPAR",
        "EMA/63929/2024 page 46: 'an age of 50 years'). Exponents 0.176",
        "(cefepime) and 0.146 (enmetazobactam) - FDA Table 101 rows 'Age on",
        "FEP V1' and 'Age on ENM V1'. Age was tested on cefepime clearance",
        "and NOT retained (EMA Table 5 prints '-' for 'Age on Cl'). The FDA",
        "reviewer concluded the age effect is not clinically meaningful and",
        "requires no dose adjustment (FDA Table 98). Time-fixed per subject;",
        "analysis-set mean 51.62 years, range 17-94 years (FDA Table 100)."
      ),
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = paste(
        "DE-INDEXED (absolute) estimated glomerular filtration rate in",
        "mL/min - i.e. the MDRD eGFR multiplied back up by BSA/1.73, NOT the",
        "BSA-normalised mL/min/1.73 m^2 form"
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW, NON-BSA-NORMALISED renal function, in mL/min. This is the",
        "un-normalised member of the CRCL family (same normalisation as",
        "Delattre_2010_amikacin and Chen_2023_nemonoxacin); supplying a",
        "BSA-normalised mL/min/1.73 m^2 value instead would silently",
        "mis-scale every clearance. The de-indexed form is what the FDA",
        "reviewer specifically requested: BSA-indexed eGFR underestimates",
        "measured GFR in subjects whose BSA exceeds 1.73 m^2, so an",
        "information request of 15 September 2023 asked the applicant to",
        "redo the renal-impairment simulations on de-indexed eGFR (FDA",
        "Integrated Review, sponsor page 138). Power covariate on the",
        "clearance of both drugs, normalised to 100 mL/min (EMA EPAR",
        "EMA/63929/2024 page 46), with exponents 0.834 (cefepime) and 0.898",
        "(enmetazobactam); the abstract reads exponents this close to 1 as",
        "evidence that 'renal filtration is the predominant clearance",
        "mechanism'. Enmetazobactam ALSO carries a de-indexed eGFR power",
        "term of 0.259 on its peripheral volume (FDA Table 101 row",
        "'De-indexed eGFR on ENM V2'); cefepime does not (EMA Table 5 prints",
        "'-' for 'De-indexed eGFR on V2'). Renal function is the only",
        "covariate that drives a labelled dose adjustment. Time-fixed per",
        "subject in the source analysis; analysis-set mean 87.98 mL/min,",
        "range 5.21-195.84 mL/min (FDA Table 100)."
      ),
      source_name        = "de-indexed eGFR"
    ),
    SEXF = list(
      description        = "Sex, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female) - the published reference individual is FEMALE, so the enmetazobactam V2 estimate of 5.28 L is the female value and males carry the multiplicative shift",
      notes              = paste(
        "Enmetazobactam PERIPHERAL volume only, as exp(0.198) = 1.219 in",
        "males (FDA Table 101 row 'Gender on ENM V2', parameter",
        "beta_ENM,V2_GENDER_Male = 0.198). Because the coefficient is",
        "carried by the MALE level while the canonical column codes FEMALE",
        "as 1, the model() term is written on (1 - SEXF) and the ini()",
        "parameter is named e_sexmale_vp_enm rather than the usual e_sexf_*",
        "- the name follows the level that carries the effect, as in",
        "e_race_chinese_vc_avi of Xie_2025_aztreonam_avibactam. Sex was NOT",
        "a significant covariate on any cefepime parameter (EMA EPAR page",
        "46: 'sex was not identified as a significant covariate on cefepime",
        "PK and thus, was not used in defining a reference individual').",
        "Analysis set 310 female (47%) / 353 male (53%) (FDA Table 100)."
      ),
      source_name        = "GENDER"
    ),
    DIS_CUTI = list(
      description        = "Infected patient indicator, 1 = complicated urinary tract infection or acute pyelonephritis, 0 = healthy volunteer",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy Phase 1 volunteer)",
      notes              = paste(
        "SINGLE POOLED infection indicator: the source analysis carries one",
        "coefficient for the whole infected cohort rather than separate cUTI",
        "and acute-pyelonephritis levels, so per the DIS_CUTI register entry",
        "this column is used ALONE and DIS_AP is deliberately not paired",
        "with it. Set it to 1 for every cUTI or AP patient and to 0 only for",
        "healthy volunteers. FDA Table 100 splits the 663-subject analysis",
        "set into 132 healthy (20%) and 531 infected - 252 acute",
        "pyelonephritis (38%), 17 cUTI (3%), 111 cUTI with a removable",
        "source (17%) and 151 cUTI without a removable source but with other",
        "risk factors (23%). This column does TWO things in the model. (1) A",
        "fixed effect on the enmetazobactam central volume only, as",
        "exp(0.145) = 1.156 (FDA Table 101 row 'cUTI infection on ENM V1'),",
        "which takes the reference 13.2 L to the 15.3 L that both the EMA",
        "report (page 47, 'For infected subjects, the estimated central",
        "volume of distribution was slightly larger with 15.3 L') and the",
        "IDWeek abstract print for an infected patient - this agreement is",
        "what pins the effect to the exponential form exp(beta) rather than",
        "the proportional form (1 + beta), which would give 15.1 L. (2) It",
        "SWITCHES the inter-individual variances: every IIV except the two",
        "peripheral-volume terms has a separate healthy-volunteer and",
        "infected-patient magnitude, and the infected magnitudes are two- to",
        "nine-fold larger, matching the abstract's 'Variability was higher",
        "in cUTI patients although differences in mean ENM or FEP exposures",
        "between healthy subjects and cUTI patients were negligible'."
      ),
      source_name        = "INFECTION_cUTI"
    )
  )

  # Screened in the source covariate analysis and NOT retained in the final
  # model, so they carry no coefficient and are documented rather than
  # declared (FDA Integrated Review, sponsor page 134, 'Covariate Analysis':
  # "BSA, body mass index, creatinine clearance, albumin, bilirubin and race
  # were not identified as statistically significant covariates in the
  # cefepime-enmetazobactam population PK model.").
  covariatesDataExcluded <- list(
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Formally evaluated, not retained. Analysis-set mean 1.89 m^2, range 1.37-2.63 (FDA Table 100)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Formally evaluated, not retained. Analysis-set mean 26.54 kg/m^2, range 16.9-45.9 (FDA Table 100)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Formally evaluated, not retained. EMA Table 5 explicitly prints '-' for 'Albumin level on V2' of cefepime. Analysis-set mean 40.38 g/L, range 21-51 (FDA Table 100)."
    ),
    BILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Formally evaluated, not retained. Analysis-set mean 11.03 umol/L, range 3-51.5 (FDA Table 100)."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race COULD NOT BE TESTED rather than merely being rejected: FDA",
        "Table 98 records that 'Race could not be tested as a predictor of",
        "IIV in the population PK analyses because 97% of the population was",
        "White', and the US label states there was insufficient information",
        "to evaluate the effect of race. Analysis set 641/663 White (97%)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = "663 in the analysis set; 588 contributed 5,070 cefepime plasma observations and 649 contributed 6,342 enmetazobactam plasma observations",
    n_studies      = 5,
    studies        = "AT-101, AT-102 and AT-103 (Phase 1, n = 83 / 30 / 19), AT-201 (Phase 2, n = 43) and AT-301 / ALLIUM (Phase 3, n = 488)",
    age_range      = "17-94 years (mean 51.62)",
    weight_range   = "45-135 kg (mean 76.52)",
    sex_female_pct = 47,
    race_ethnicity = "97% White (641/663); 4 Asian, 3 American Indian or Alaska Native, 1 Black or African American, 14 other",
    disease_state  = "complicated urinary tract infection or acute pyelonephritis (531/663, 80%) and healthy volunteers (132/663, 20%)",
    renal_function = "de-indexed eGFR mean 87.98 mL/min, range 5.21-195.84 mL/min; BSA-indexed eGFR mean 81.03 mL/min/1.73 m^2, range 4.4-166",
    dose_range     = "cefepime 1-2 g with enmetazobactam 0.25-1 g as 2-hour or 4-hour intravenous infusions, q8h to q24h depending on renal function; the approved adult regimen is 2 g cefepime / 0.5 g enmetazobactam q8h over 2 hours",
    route          = "intravenous infusion",
    notes          = paste(
      "Cefepime and enmetazobactam are co-formulated in a fixed 4:1 mass",
      "ratio and were fitted SIMULTANEOUSLY, which is what makes the",
      "cross-drug random-effect correlations estimable; the abstract",
      "highlights that 'clearances and volumes of distribution were strongly",
      "correlated within individuals between the two compounds (correlation",
      "coefficients > 0.9)'. Plasma protein binding is 20% for cefepime and",
      "negligible for enmetazobactam (US label Table 5), so the free",
      "concentrations that the PK/PD targets are defined on are 0.80 x Cc",
      "and 1.00 x Cc_enm. The model was written in Mlxtran and estimated by",
      "SAEM in Monolix."
    )
  )

  ini({
    # =====================================================================
    # REFERENCE INDIVIDUAL for every fixed effect below: body weight 70 kg,
    # de-indexed eGFR 100 mL/min, age 50 years, Sex = Female, Infection =
    # healthy (EMA EPAR EMA/63929/2024 pages 46-48, notes to Tables 5 and
    # 6). The IDWeek abstract's table instead reports "typical values for a
    # 70 kg cUTI patient", which differ from these estimates for exactly
    # one parameter - enmetazobactam Vc, 15.3 L there versus the 13.2 L
    # reference value here, the ratio being the exp(0.145) cUTI effect.
    #
    # CEFEPIME structural parameters (FDA NDA 216165 Integrated Review
    # Table 101, "Fixed effects"; EMA EPAR Table 5). Cefepime is the
    # unsuffixed parent.
    # =====================================================================
    lcl <- log(5.95); label("Cefepime clearance at de-indexed eGFR = 100 mL/min (L/h)")                  # Table 101: Cl_FEP = 5.95 L/h (RSE 1%)
    lvc <- log(11.2); label("Cefepime central volume of distribution at WT = 70 kg, AGE = 50 y (L)")     # Table 101: V1_FEP = 11.2 L (RSE 2.2%)
    lq  <- fixed(log(7.22)); label("Cefepime inter-compartmental clearance (L/h)")                       # Table 101: Q_FEP = 7.22 L/h, RSE column reads "fixed"
    lvp <- log(5.7);  label("Cefepime peripheral volume of distribution (L)")                            # Table 101: V2_FEP = 5.7 L (RSE 2.8%)

    # =====================================================================
    # ENMETAZOBACTAM structural parameters (FDA Table 101; EMA Table 6).
    # =====================================================================
    lcl_enm <- log(7.68); label("Enmetazobactam clearance at de-indexed eGFR = 100 mL/min (L/h)")                     # Table 101: Cl_ENM = 7.68 L/h (RSE 1.1%)
    lvc_enm <- log(13.2); label("Enmetazobactam central volume of distribution at WT = 70 kg, AGE = 50 y, healthy (L)") # Table 101: V1_ENM = 13.2 L (RSE 2.1%)
    lq_enm  <- fixed(log(7.16)); label("Enmetazobactam inter-compartmental clearance (L/h)")                          # Table 101: Q_ENM = 7.16 L/h, RSE column reads "fixed"
    lvp_enm <- log(5.28); label("Enmetazobactam peripheral volume of distribution at de-indexed eGFR = 100 mL/min, female (L)") # Table 101: V2_ENM = 5.28 L (RSE 3.3%)

    # Both inter-compartmental clearances are FIXED, not estimated. EMA EPAR
    # page 46: "The sparse sampling in AT-301 did not allow to characterise
    # the distribution into the peripheral compartment. The estimate of Q
    # was therefore fixed to the final estimate from the phase I+II
    # analysis." Consistent with FDA Table 101, whose RSE column prints the
    # literal word "fixed" for both Q rows and "(-)" for their %CV.

    # =====================================================================
    # COVARIATE COEFFICIENTS (FDA Table 101, "Covariate coefficients"; EMA
    # Tables 5 and 6).
    #
    # FUNCTIONAL FORM. The analysis was run in Monolix, where a covariate
    # effect enters a log-normally distributed parameter additively on the
    # log scale. Continuous covariates are entered log-transformed (the FDA
    # parameter names carry the transform marker, e.g. beta_FEP,CL_teGFR
    # and beta_ENM,V1_tAGE), which makes each of those a POWER model on the
    # covariate/reference ratio; categorical covariates enter as an
    # indicator, which makes each a multiplicative exp(beta). Three
    # independent facts confirm this reading rather than a proportional
    # (1 + beta) one:
    #   (i)  the cUTI effect on enmetazobactam Vc reproduces the separately
    #        published infected value exactly - 13.2 * exp(0.145) = 15.26,
    #        printed as 15.3 L by both the EMA report and the abstract,
    #        whereas 13.2 * (1 + 0.145) = 15.1 L does not;
    #   (ii) the abstract reads the ~0.83-0.90 eGFR coefficients as showing
    #        "renal filtration is the predominant clearance mechanism",
    #        which is a statement about a near-proportional POWER exponent;
    #        under a linear model 0.834 would mean +83% clearance per
    #        mL/min, which is absurd;
    #   (iii) inverting the power form on FDA Table 104's simulated
    #        steady-state AUCs recovers a de-indexed eGFR reference of
    #        ~100 mL/min, which is the value the EMA report states.
    # =====================================================================
    e_crcl_cl     <- 0.834; label("Cefepime de-indexed eGFR power exponent on CL, (CRCL/100) (unitless)")             # Table 101: beta_FEP,CL_teGFR = 0.834 (RSE 2.1%, p < 2.2e-16)
    e_age_vc      <- 0.176; label("Cefepime age power exponent on Vc, (AGE/50) (unitless)")                           # Table 101: beta_FEP,V1_tAGE = 0.176 (RSE 19%, p = 1.4e-07)
    e_wt_vc       <- 0.802; label("Cefepime body-weight power exponent on Vc, (WT/70) (unitless)")                    # Table 101: beta_FEP,V1_tBW = 0.802 (RSE 8.8%, p < 2.2e-16)

    e_crcl_cl_enm <- 0.898; label("Enmetazobactam de-indexed eGFR power exponent on CL, (CRCL/100) (unitless)")       # Table 101: beta_ENM,CL_teGFR = 0.898 (RSE 2%, p < 2.2e-16)
    e_age_vc_enm  <- 0.146; label("Enmetazobactam age power exponent on Vc, (AGE/50) (unitless)")                     # Table 101: beta_ENM,V1_tAGE = 0.146 (RSE 20%, p = 3.8e-07)
    e_wt_vc_enm   <- 0.618; label("Enmetazobactam body-weight power exponent on Vc, (WT/70) (unitless)")              # Table 101: beta_ENM,V1_tBW = 0.618 (RSE 9.6%, p < 2.2e-16)
    e_cuti_vc_enm <- 0.145; label("Enmetazobactam log-scale shift in Vc for an infected (cUTI or AP) patient (unitless)") # Table 101: beta_ENM,V1_INFECTION_cUTI = 0.145 (RSE 17%, p = 3.8e-09)
    e_crcl_vp_enm <- 0.259; label("Enmetazobactam de-indexed eGFR power exponent on Vp, (CRCL/100) (unitless)")       # Table 101: beta_ENM,V2_teGFR = 0.259 (RSE 21%, p = 3.2e-06)
    e_sexmale_vp_enm <- 0.198; label("Enmetazobactam log-scale shift in Vp for a male subject (unitless)")            # Table 101: beta_ENM,V2_GENDER_Male = 0.198 (RSE 19%, p = 7.9e-08)

    # =====================================================================
    # INTER-INDIVIDUAL VARIABILITY (FDA Table 101, "Standard deviations" and
    # "Correlations").
    #
    # SCALE OF THE PRINTED COLUMN. Table 101 gives BOTH the Monolix omega
    # (an SD on the log scale, under "Standard deviations") and a %CV
    # (in parentheses beside each fixed effect), so the scale is not a
    # judgement call - it is over-determined nine times over. Every one of
    # the nine estimated omegas reproduces its printed %CV under the
    # log-normal identity CV = sqrt(exp(omega^2) - 1): 0.15 -> 15%,
    # 0.356 -> 37%, 0.0549 -> 5%, 0.474 -> 50%, 0.0593 -> 6%, 0.161 -> 16%,
    # 0.297 -> 30%, 0.42 -> 44%, 0.195 -> 20%. The variances below are
    # therefore omega^2 with the printed omegas taken verbatim.
    #
    # STRATIFICATION BY INFECTION STATUS. Cefepime and enmetazobactam CL and
    # enmetazobactam Vc each have a separate healthy and infected variance;
    # cefepime Vc has an infected variance only (Table 101 prints "-" for
    # both the estimate and the RSE of omega_FEP,V1_healthy); the two
    # peripheral volumes have one variance spanning both strata (Table 101
    # labels them "FEP V2 in healthy and infected" and "ENM V2 in healthy
    # and infected"). The etas are multiplexed on DIS_CUTI in model().
    #
    # INFECTED BLOCK: a full 4x4 across cefepime CL, cefepime Vc,
    # enmetazobactam CL and enmetazobactam Vc. All six correlations are
    # printed, so the block is complete.
    #
    # POSITIVE-DEFINITENESS REPAIR. The 4x4 built from the correlations
    # EXACTLY as printed to three decimals is very slightly INDEFINITE
    # (smallest eigenvalue -1.72e-04), so chol() fails and rxSolve() cannot
    # sample from it. This is a rounding artifact, not a contradiction: each
    # printed correlation sits within 0.0005 - i.e. inside its own rounding
    # interval - of the feasible region, so the unrounded Monolix matrix was
    # positive definite and three-decimal rounding is what broke it. The
    # covariances below use the NEAREST positive-definite correlation matrix
    # (Higham projection, Matrix::nearPD with corr = TRUE, keepDiag = TRUE),
    # which moves no correlation by more than 6.7e-05 - an order of
    # magnitude INSIDE the +/-0.0005 rounding interval - so every repaired
    # correlation rounds back to the published three-decimal value exactly.
    # The repaired correlations are 0.260949, 0.931950, 0.261058, 0.339056,
    # 0.954933 and 0.437934 against printed 0.261, 0.932, 0.261, 0.339,
    # 0.955 and 0.438. The vignette re-derives this and asserts both the
    # round-trip and positive-definiteness.
    #
    # Order of the block: etalcl_inf, etalvc_inf, etalcl_enm_inf,
    # etalvc_enm_inf. Diagonals are omega^2 from the printed omegas:
    #   omega_FEP,CL_infected = 0.297 -> 0.088209   (RSE 3.4%)
    #   omega_FEP,V1_infected = 0.42  -> 0.1764     (RSE 4.5%)
    #   omega_ENM,CL_infected = 0.356 -> 0.126736   (RSE 3.3%)
    #   omega_ENM,V1_infected = 0.474 -> 0.224676   (RSE 3.9%)
    # Off-diagonals are repaired-correlation x SD x SD, from Table 101 rows:
    #   "CEF V1, FEP Cl in infected"  0.261 (RSE 19)  -> 0.03255083
    #   "ENM CL, FEP Cl in infected"  0.932 (RSE 1.1) -> 0.09853692
    #   "ENM CL, FEP V1 in infected"  0.339 (RSE 15)  -> 0.05069571
    #   "ENM V1, FEP Cl in infected"  0.261 (RSE 19)  -> 0.03675124
    #   "ENM V1, FEP V1 in infected"  0.955 (RSE 2.4) -> 0.19010800
    #   "ENM V1, ENM Cl in infected"  0.438 (RSE 9.3) -> 0.07389876
    # (The two rows printed as 0.261 with RSE 19 are the cefepime
    # within-drug Vc-CL correlation and the enmetazobactam-Vc/cefepime-CL
    # cross-drug correlation. EMA Table 5 independently prints the first as
    # "V1, Cl in cUTI = 0.261"; the second is cross-drug and so appears only
    # in the FDA table, where it is reported with the identical value.)
    # =====================================================================
    etalcl_inf + etalvc_inf + etalcl_enm_inf + etalvc_enm_inf ~
      c(0.08820900,
        0.03255083, 0.17640000,
        0.09853692, 0.05069571, 0.12673600,
        0.03675124, 0.19010800, 0.07389876, 0.22467600)

    # HEALTHY BLOCK: enmetazobactam CL and Vc, correlation 0.983 (Table 101
    # row "ENM V1, ENM Cl in healthy", RSE 42; EMA Table 6 "V1, Cl in
    # healthy"). Positive definite as printed (eigenvalues 2.54e-02 and
    # 8.99e-05), so no repair is applied here.
    #   omega_ENM,CL_healthy = 0.15    -> 0.0225      (RSE 7.6%)
    #   omega_ENM,V1_healthy = 0.0549  -> 0.00301401  (RSE 59%)
    #   covariance = 0.983 * 0.15 * 0.0549            -> 0.00809505
    etalcl_enm_hlth + etalvc_enm_hlth ~
      c(0.02250000,
        0.00809505, 0.00301401)

    # Cefepime CL in healthy volunteers carries no reported correlation with
    # anything, so it is an independent eta.
    etalcl_hlth ~ 0.02592100   # omega_FEP,CL_healthy = 0.161 -> 0.161^2 (RSE 11%)

    # Peripheral-volume IIV, one variance each spanning BOTH strata.
    etalvp     ~ 0.03802500   # omega_FEP,V2 = 0.195   -> 0.195^2   (RSE 14%)
    etalvp_enm ~ 0.00351649   # omega_ENM,V2 = 0.0593  -> 0.0593^2  (RSE 72%)

    # NO eta on either inter-compartmental clearance: both Q values are
    # fixed, and Table 101 lists no omega_Q row for either drug.
    #
    # NO healthy-stratum eta on cefepime Vc: Table 101 row "FEP V1 in
    # healthy" prints "-" in both the estimate and the RSE column. Healthy
    # volunteers therefore take the typical cefepime central volume.

    # =====================================================================
    # RESIDUAL VARIABILITY (FDA Table 101, "Observational error"; EMA
    # Tables 5 and 6). Monolix's combined error model is
    # y = f + (a + b*f)*eps, so a is an additive SD in concentration units
    # and b is a proportional SD.
    #
    # UNIT CONVERSION. Table 101 reports the enmetazobactam constant error
    # in ng/mL (a_ENM = 40.9 ng/mL) because the analysis dataset held
    # concentrations in ng/mL - the same units as the FDA's Figure 14 axis
    # label, "Observed Versus Population Predicted ... Plasma Concentrations
    # in ng/mL". This model works in mg/L (= ug/mL), so the additive SD is
    # 40.9 ng/mL / 1000 = 0.0409 mg/L. The proportional terms are unitless
    # and carry over unchanged.
    # =====================================================================
    propSd     <- 0.253;  label("Cefepime proportional residual SD (fraction)")            # Table 101: b_FEP = 0.253 (RSE 1.2%)
    propSd_enm <- 0.236;  label("Enmetazobactam proportional residual SD (fraction)")      # Table 101: b_ENM = 0.236 (RSE 1.1%)
    addSd_enm  <- 0.0409; label("Enmetazobactam additive residual SD (mg/L)")              # Table 101: a_ENM = 40.9 ng/mL (RSE 12%) = 0.0409 mg/L

    # Cefepime has NO additive residual component: Table 101 lists a
    # constant error for enmetazobactam only, and EMA Table 5 prints "-"
    # for the cefepime "Constant error" row. Its error model is purely
    # proportional.
  })

  model({
    # ------------------------------------------------------------------
    # 1. Infection-stratum multiplexing of the random effects.
    #
    # Each of these selects the healthy or the infected eta for the subject
    # via the 0/1 DIS_CUTI indicator, rather than through a branch, so the
    # expression stays smooth and finite. Cefepime Vc has no healthy-arm
    # eta, so an uninfected subject takes the typical value exactly.
    # ------------------------------------------------------------------
    infected <- DIS_CUTI
    healthy  <- 1 - DIS_CUTI

    eta_cl     <- etalcl_inf     * infected + etalcl_hlth     * healthy
    eta_vc     <- etalvc_inf     * infected
    eta_cl_enm <- etalcl_enm_inf * infected + etalcl_enm_hlth * healthy
    eta_vc_enm <- etalvc_enm_inf * infected + etalvc_enm_hlth * healthy

    # ------------------------------------------------------------------
    # 2. Cefepime individual PK parameters. Every covariate enters on the
    #    log scale (power terms on the log-transformed continuous
    #    covariates), so each exp() argument is the Monolix individual
    #    log-parameter.
    # ------------------------------------------------------------------
    cl <- exp(lcl + eta_cl) * (CRCL / 100)^e_crcl_cl
    vc <- exp(lvc + eta_vc) * (WT / 70)^e_wt_vc * (AGE / 50)^e_age_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # ------------------------------------------------------------------
    # 3. Enmetazobactam individual PK parameters. The cUTI effect on Vc and
    #    the male effect on Vp are categorical, so they are exp(beta)
    #    multipliers; the male term is written on (1 - SEXF) because the
    #    published reference individual is FEMALE.
    # ------------------------------------------------------------------
    cl_enm <- exp(lcl_enm + eta_cl_enm) * (CRCL / 100)^e_crcl_cl_enm

    vc_enm <- exp(lvc_enm + eta_vc_enm) *
      (WT / 70)^e_wt_vc_enm *
      (AGE / 50)^e_age_vc_enm *
      exp(e_cuti_vc_enm * DIS_CUTI)

    q_enm  <- exp(lq_enm)

    vp_enm <- exp(lvp_enm + etalvp_enm) *
      (CRCL / 100)^e_crcl_vp_enm *
      exp(e_sexmale_vp_enm * (1 - SEXF))

    # ------------------------------------------------------------------
    # 4. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_enm <- cl_enm / vc_enm
    k12_enm <- q_enm  / vc_enm
    k21_enm <- q_enm  / vp_enm

    # ------------------------------------------------------------------
    # 5. Four-compartment system: two independent two-compartment IV
    #    dispositions. Cefepime and enmetazobactam are co-formulated in one
    #    vial at a fixed 4:1 mass ratio and infused together, but they do
    #    not interconvert and neither inhibits the other's disposition, so
    #    each takes its own dose into its own central compartment as a
    #    zero-order infusion.
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1

    d/dt(central_enm)     <- -(kel_enm + k12_enm) * central_enm + k21_enm * peripheral1_enm
    d/dt(peripheral1_enm) <-   k12_enm * central_enm            - k21_enm * peripheral1_enm

    # ------------------------------------------------------------------
    # 6. Observations. Doses in mg and volumes in L give concentrations in
    #    mg/L (= ug/mL). These are TOTAL plasma concentrations; multiply by
    #    the unbound fractions 0.80 (cefepime, 20% protein bound) and 1.00
    #    (enmetazobactam, binding negligible) to obtain the free
    #    concentrations that the %fT > MIC and %fT > CT targets are defined
    #    on. See population$notes.
    # ------------------------------------------------------------------
    Cc     <- central     / vc
    Cc_enm <- central_enm / vc_enm

    Cc     ~ prop(propSd)
    Cc_enm ~ add(addSd_enm) + prop(propSd_enm)
  })
}
