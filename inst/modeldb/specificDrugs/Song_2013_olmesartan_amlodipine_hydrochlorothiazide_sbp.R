Song_2013_olmesartan_amlodipine_hydrochlorothiazide_sbp <- function() {
  description <- paste(
    "Steady-state exposure-response model for the change from baseline in seated trough SYSTOLIC blood pressure",
    "under the triple fixed-dose combination CS-8635 (olmesartan / amlodipine / hydrochlorothiazide) and its",
    "component mono- and dual-therapies, in adults with hypertension. The response is the sum of a",
    "study-specific placebo effect, three independent monotherapy drug effects driven by each component's",
    "steady-state AUC (saturable Emax for olmesartan and amlodipine, linear slope for hydrochlorothiazide), and",
    "a set of pairwise and three-way interaction terms proportional to the PRODUCT of the monotherapy effects,",
    "which make the combination sub-additive. Structurally the same model as the diastolic sibling but fitted",
    "independently, and the retained covariate sets differ: the systolic model keeps a baseline effect on the",
    "hydrochlorothiazide term, a female effect on the amlodipine term, a Hispanic-ethnicity effect on the",
    "CS8635-A-U301 placebo, and the olmesartan x hydrochlorothiazide pairwise interaction, while dropping the",
    "age effect on the olmesartan term and the Black-race effect on the placebo. There are no ODE states and no",
    "dosing events -- exposure enters entirely through the AUC_OLM / AUC_AML / AUC_HCTZ covariate columns, which",
    "the companion population PK models generate as AUCss = Dose / (CL/F). Sister model files from the same",
    "paper: modellib('Song_2013_olmesartan_amlodipine_hydrochlorothiazide_dbp'),",
    "modellib('Song_2013_olmesartan'), modellib('Song_2013_amlodipine'),",
    "modellib('Song_2013_hydrochlorothiazide')."
  )
  reference <- paste(
    "Song S, Carrothers TJ, Moore H, Green M, Miller R, Rohatagi S, Lee J, Wang A, Melino M, Patel M,",
    "Heyrman R, Salazar DE.",
    "Model-supported development of CS-8635: a fixed-dose combination of olmesartan, amlodipine, and",
    "hydrochlorothiazide.",
    "Clin Pharmacol Drug Dev. 2013;2(2):103-112.",
    "doi:10.1002/cpdd.17.",
    "Structural equations from main-text Equations 1-7; parameter values from Supplemental Table S7, whose row",
    "labels supersede the Systolic column of main-text Table 3 (see the model file comments and the vignette",
    "deviations section for the mislabelled Table 3 row).",
    sep = " "
  )
  vignette <- "Song_2013_olmesartan_amlodipine_hydrochlorothiazide"
  units <- list(
    time          = "h",
    dosing        = "n/a (no drug-dosing events; exposure enters as the covariates AUC_OLM, AUC_AML, AUC_HCTZ)",
    concentration = "change from baseline in seated trough systolic blood pressure, in mmHg (not a drug concentration)"
  )

  # The additive inter-subject random effect sits on the TOTAL response
  # (NONMEM Y = placebo + drug + ETA + EPS), not on any single ini()
  # parameter, so it has no 'l<param>' partner to be named after.
  paper_specific_etas <- "etadsbp"

  covariateData <- list(
    AUC_OLM = list(
      description        = "Steady-state area under the olmesartan plasma concentration-time curve over the 24 h dosing interval.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Song 2013 M&S Methods: 'area-under-the-curve (AUCss), calculated as Dose divided by Apparent",
        "Clearance', taken from the post hoc Bayesian individual clearance of the companion population PK model",
        "modellib('Song_2013_olmesartan'). Time-fixed per subject and per regimen. Enters as the saturable Emax",
        "term emax_om * AUC_OLM / (AUC_OLM + eauc50_om) (Song 2013 Equation 5). Note eauc50_om is 1590 ng*h/mL",
        "here against 1850 in the diastolic sibling. Set AUC_OLM = 0 to remove olmesartan from the regimen.",
        "Worked typical value: 40 mg / 6.32 L/h = 6329 ng*h/mL."
      ),
      source_name        = "AUCOM"
    ),
    AUC_AML = list(
      description        = "Steady-state area under the amlodipine plasma concentration-time curve over the 24 h dosing interval.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Dose / (CL/F) from modellib('Song_2013_amlodipine'). Enters as the saturable Emax term",
        "emax_aml * AUC_AML / (AUC_AML + eauc50_aml), with eauc50_aml = 309 ng*h/mL against 453 in the",
        "diastolic sibling -- the systolic amlodipine response saturates sooner.",
        "Worked typical value: 10 mg / 23.4 L/h = 427 ng*h/mL.",
        "Set AUC_AML = 0 to remove amlodipine from the regimen."
      ),
      source_name        = "AUCAML"
    ),
    AUC_HCTZ = list(
      description        = "Steady-state area under the hydrochlorothiazide plasma concentration-time curve over the 24 h dosing interval.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Dose / (CL/F) from modellib('Song_2013_hydrochlorothiazide'). Enters LINEARLY as",
        "slope_hctz * AUC_HCTZ / 1000 (Song 2013 Equation 4); the slope is reported in mmHg per 1000 ng*h/mL,",
        "so the division is mandatory. The systolic slope of -9.38 mmHg per 1000 ng*h/mL is nearly three times",
        "the diastolic slope of -3.3, making hydrochlorothiazide proportionally more systolic-selective than",
        "either of its partners.",
        "Worked typical value: 25 mg / 20.3 L/h = 1232 ng*h/mL.",
        "Set AUC_HCTZ = 0 to remove hydrochlorothiazide from the regimen."
      ),
      source_name        = "AUCHCTZ"
    ),
    SBP = list(
      description        = "Baseline seated trough systolic blood pressure.",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (pre-treatment baseline). Measured after 5 minutes seated, as the mean of three",
        "readings (two in study 866-318) taken at 1-minute intervals at approximately the same time of day at",
        "each visit. Enters in power form on the placebo term, (SBP / 165)^4.08, and on all THREE drug-effect",
        "terms: (SBP / 165)^1.96 for olmesartan, ^3.65 for amlodipine and ^2.82 for hydrochlorothiazide. The",
        "hydrochlorothiazide baseline effect is retained for systolic pressure only -- the Supplemental Results",
        "narrative confines it explicitly to 'higher systolic baseline blood pressure'. The centering value",
        "165 mmHg is the mean baseline SBP of the pooled exposure-response dataset (Song 2013 Table 2,",
        "'All data' row, n = 4873), used as a proxy for the median that Song 2013 specifies in Equation 6 but",
        "never prints."
      ),
      source_name        = "Baseline BP"
    ),
    AGE = list(
      description        = "Subject age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters in power form on the placebo term ONLY, (AGE / 54.8)^-0.746 -- note the",
        "sign is opposite to the diastolic sibling's +1.37, so older subjects show a SMALLER systolic placebo",
        "response and a LARGER diastolic one. No age effect was retained on any systolic drug-effect term.",
        "The centering value 54.8 years is the mean age of the pooled exposure-response dataset",
        "(Song 2013 Table 2, 'All data' row)."
      ),
      source_name        = "Age"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters in power form on the amlodipine drug-effect term only,",
        "(WT / 94.9)^-0.586, i.e. 'Treatment effects of AML were lower in subjects with higher body weights'",
        "(Song 2013 Results). The centering value 94.9 kg is the mean body weight of the pooled",
        "exposure-response dataset (Song 2013 Table 2, 'All data' row), used as a proxy for the unprinted median."
      ),
      source_name        = "Weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed per subject. Enters additively-fractionally (Song 2013 Equation 7) on the AMLODIPINE",
        "drug-effect term, factor (1 + 0.301), i.e. female subjects show about 30% greater systolic",
        "blood-pressure lowering from amlodipine -- the Supplemental Results narrative confirms the reading:",
        "'the female subjects might show 30 % higher systolic blood pressure lowering effects of AML'.",
        "Main-text Table 3 prints this value in a row simply labelled 'Sex' grouped under 'Placebo effects';",
        "Table S7 labels it 'Effect of sex on Drug Effect of AML' and the narrative agrees with Table S7, so",
        "the effect is placed on the amlodipine term here. No sex effect was retained in the diastolic model."
      ),
      source_name        = "Sex"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator, 1 = Black, 0 = other.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black)",
      notes              = paste(
        "Time-fixed per subject. Enters additively-fractionally (Song 2013 Equation 7) on the olmesartan",
        "drug-effect term, factor (1 - 0.393): Black subjects have less renin-angiotensin-system activation and",
        "respond less well to angiotensin receptor blockers. The systolic attenuation (39.3%) is larger than the",
        "diastolic one (26.3%). Unlike the diastolic model there is no Black-race effect on the systolic placebo",
        "term. Song 2013's Results prose quotes '39.6% and 29.3% less OM treatment effect on SeDBP and SeSBP'",
        "and the Supplemental Results quotes '26 and 50 %' for the same two endpoints; both parameter tables",
        "agree on -0.393 (systolic) and -0.263 (diastolic), so the tabulated values are used."
      ),
      source_name        = "Black race"
    ),
    RACE_HISPANIC = list(
      description        = "Hispanic / Latino ethnicity indicator, 1 = Hispanic, 0 = non-Hispanic.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Hispanic)",
      notes              = paste(
        "Time-fixed per subject. Enters additively-fractionally (Song 2013 Equation 7) on the PLACEBO term for",
        "subjects in CS8635-A-U301 only, factor (1 - 0.554): Hispanic subjects in that study showed roughly half",
        "the placebo response of the rest of the study. Song 2013 reports race and ethnicity in a single",
        "W:B:H:A:O breakdown (Table 2), so Hispanic is a category of the same variable as White and Black rather",
        "than a separate ethnicity dimension; the canonical binary encoding is nonetheless identical. Retained",
        "in the systolic model only. n = 365 Hispanic subjects in CS8635-A-U301."
      ),
      source_name        = "Hispanic race"
    ),
    STUDY_CS8635_A_U301 = list(
      description        = "Indicator for enrolment in study CS8635-A-U301 (TRINITY), the pivotal phase III CS-8635 triple-combination trial.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the other two pooled phase III studies)",
      notes              = paste(
        "Selects the study-specific placebo effect of -4.20 mmHg AND gates the Hispanic-ethnicity effect on the",
        "placebo term, which Song 2013 retained only within this study. The three study indicators partition the",
        "pooled dataset and must sum to 1 for every subject; leaving all three at 0 produces a zero placebo",
        "effect, a state the model was never fitted to. n = 2458 in the exposure-response dataset."
      ),
      source_name        = "Study CS8635-A-U301"
    ),
    STUDY_CS8663_A_U301 = list(
      description        = "Indicator for enrolment in study CS8663-A-U301 (COACH), the pivotal phase III olmesartan / amlodipine dual-combination trial.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the other two pooled phase III studies)",
      notes              = paste(
        "Selects the study-specific placebo effect of -3.45 mmHg. Unlike the diastolic sibling, the systolic",
        "model carries no study-gated race effect on this study's placebo term.",
        "n = 1920 in the exposure-response dataset."
      ),
      source_name        = "Study CS8663-A-U301"
    ),
    STUDY_866_318 = list(
      description        = "Indicator for enrolment in study 866-318 (written SE866-318 in the parameter tables), the phase III olmesartan / hydrochlorothiazide dual-combination trial.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the other two pooled phase III studies)",
      notes              = paste(
        "Selects the study-specific placebo effect of -5.26 mmHg, the largest of the three studies.",
        "n = 495 in the exposure-response dataset; this study also had the lowest mean baseline SBP",
        "(154 mmHg) and took only two rather than three readings per visit."
      ),
      source_name        = "Study SE866-318"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 4873L,
    n_studies      = 3L,
    age_range      = "mean 54.8 (SD 11) years",
    weight_range   = "mean 94.9 (SD 22) kg",
    sex_female_pct = 46.1,
    race_ethnicity = c(White = 58.9, Black = 25.0, Hispanic = 13.4, Asian = 1.9, Other = 0.8),
    disease_state  = "Adults with hypertension; 14.1% diabetic. Mean baseline seated trough BP 165/102 mmHg",
    dose_range     = paste(
      "olmesartan medoxomil 0-40 mg, amlodipine 0-10 mg and hydrochlorothiazide 0-25 mg once daily, in the",
      "mono-, dual- and triple-combination arms of three factorial phase III trials"
    ),
    regions        = "United States and Europe",
    notes          = paste(
      "Demographics from Song 2013 Table 2 'All data' row (n = 4873 = 495 + 1920 + 2458), which is the pooled",
      "exposure-response dataset: 866-318 (n = 495), CS8663-A-U301 / COACH (n = 1920) and CS8635-A-U301 /",
      "TRINITY (n = 2458). Pharmacokinetic data were available for 1471 of these subjects; the remainder",
      "contributed blood-pressure observations only, with their exposures predicted from the population PK",
      "models. Sex and race percentages here are computed from the Table 2 M:F and W:B:H:A:O counts. A fourth",
      "study, CS866-419 (n = 145), appears in Supplemental Table S2 but is NOT part of the final",
      "exposure-response dataset. The primary endpoint visit was Week 12 Period II for CS8635-A-U301 and Week 8",
      "Period II for the other two, by which time therapeutic steady state had been reached. Estimation used",
      "FOCE in NONMEM VI level 1."
    )
  )

  ini({
    # ==================================================================
    # All values from Song 2013 Supplemental Table S7 ('Parameter
    # Estimates for Final Exposure-Response Model for SBP'). Table S7 is
    # preferred over the Systolic column of main-text Table 3 because
    # Table 3 groups the sex effect under 'Placebo effects' while Table S7
    # and the Supplemental Results narrative both place it on the
    # amlodipine drug effect; see the note on e_sexf_aml below.
    # ==================================================================

    # ---- Placebo effect, one value per contributing study --------------
    pbo_8635    <- -4.20 ; label("Placebo effect on dSBP in study CS8635-A-U301 (mmHg)")  # Table S7 (SE 23% CV)
    pbo_8663    <- -3.45 ; label("Placebo effect on dSBP in study CS8663-A-U301 (mmHg)")  # Table S7 (SE 25% CV)
    pbo_866318  <- -5.26 ; label("Placebo effect on dSBP in study 866-318 (mmHg)")        # Table S7 (SE 19% CV)

    # ---- Covariate effects on the placebo term -------------------------
    e_sbp_pbo            <-  4.08  ; label("Power exponent for baseline SBP on the placebo effect (unitless)")                                  # Table S7 'Effect of baseline on Placebo Effect' (SE 13% CV)
    e_age_pbo            <- -0.746 ; label("Power exponent for age on the placebo effect (unitless)")                                           # Table S7 'Effect of age on Placebo Effect' (SE 28% CV)
    e_hispanic_pbo_8635  <- -0.554 ; label("Fractional Hispanic-ethnicity effect on the placebo effect within CS8635-A-U301 (unitless)")        # Table S7 'Effect of Hispanic Ethnicity on Placebo Effect in study CS8635-A-U301' (SE 40% CV)

    # ---- Olmesartan monotherapy effect (saturable, Eq. 5) --------------
    emax_om    <- -18.8 ; label("Maximum olmesartan effect on dSBP (mmHg)")                                     # Table S7 'Emax, OM' (SE 12% CV)
    eauc50_om  <-  1590 ; label("Olmesartan AUCss producing half the maximum dSBP effect (ng*h/mL)")            # Table S7 'EAUC50, OM' (SE 35% CV)
    e_sbp_om   <-  1.96 ; label("Power exponent for baseline SBP on the olmesartan effect (unitless)")          # Table S7 'Effect of baseline on Drug Effect of OM' (SE 19% CV)
    e_black_om <- -0.393 ; label("Fractional Black-race effect on the olmesartan effect (unitless)")            # Table S7 'Effect of Black Race on Drug Effect of OM' (SE 15% CV)

    # ---- Amlodipine monotherapy effect (saturable, Eq. 5) --------------
    emax_aml   <- -23.1 ; label("Maximum amlodipine effect on dSBP (mmHg)")                                     # Table S7 'Emax, AML' (SE 17% CV)
    eauc50_aml <-   309 ; label("Amlodipine AUCss producing half the maximum dSBP effect (ng*h/mL)")            # Table S7 'EAUC50, AML' (SE 37% CV)
    e_sbp_aml  <-  3.65 ; label("Power exponent for baseline SBP on the amlodipine effect (unitless)")          # Table S7 'Effect of baseline on Drug Effect of AML' (SE 9% CV)
    e_wt_aml   <- -0.586 ; label("Power exponent for body weight on the amlodipine effect (unitless)")          # Table S7 'Effect of weight on Drug Effect of AML' (SE 18% CV)
    # Table S7 labels this row 'Effect of sex on Drug Effect of AML'; main-text Table 3 prints the same 0.301
    # in a row labelled only 'Sex' that sits inside the 'Placebo effects' block. Table S7 wins: the Supplemental
    # Results narrative reads 'the female subjects might show 30 % higher systolic blood pressure lowering
    # effects of AML', which is exactly (1 + 0.301) on the amlodipine term.
    e_sexf_aml <-  0.301 ; label("Fractional female effect on the amlodipine effect (unitless)")                # Table S7 'Effect of sex on Drug Effect of AML' (SE 20% CV)

    # ---- Hydrochlorothiazide monotherapy effect (linear, Eq. 4) --------
    slope_hctz <- -9.38 ; label("Hydrochlorothiazide slope on dSBP (mmHg per 1000 ng*h/mL)")                    # Table S7 'Slope HCTZ' (SE 11% CV)
    e_sbp_hctz <-  2.82 ; label("Power exponent for baseline SBP on the hydrochlorothiazide effect (unitless)") # Table S7 'Effect of baseline on Drug Effect of HCTZ' (SE 30% CV)

    # ---- Interaction terms (Eq. 3), proportional to the PRODUCT of the
    #      monotherapy effects. Unlike the diastolic sibling, all three
    #      pairwise terms were retained here.
    ia_om_aml      <- 0.0182   ; label("Olmesartan x amlodipine interaction on dSBP (1/mmHg)")                  # Table S7 'Interaction OM*AML' (SE 18% CV)
    ia_om_hctz     <- 0.0195   ; label("Olmesartan x hydrochlorothiazide interaction on dSBP (1/mmHg)")         # Table S7 'Interaction OM*HCTZ' (SE 31% CV)
    ia_aml_hctz    <- 0.0263   ; label("Amlodipine x hydrochlorothiazide interaction on dSBP (1/mmHg)")         # Table S7 'Interaction HCTZ*AML' (SE 9% CV)
    ia_om_aml_hctz <- 0.000736 ; label("Olmesartan x amlodipine x hydrochlorothiazide interaction on dSBP (1/mmHg^2)")  # Table S7 'Interaction OM*HCTZ*AML' (SE 22% CV)

    # ---- Inter-subject variability -------------------------------------
    # Table S7 'Additive Inter-subject Variability [mmHg] = 14.0' with
    # footnote b 'Square root of ETA estimate', so the ETA VARIANCE is
    # 14.0^2 = 196 mmHg^2. The random effect is ADDITIVE on the response,
    # matching the NONMEM form Y = placebo + drug + ETA + EPS.
    etadsbp ~ 196  # Table S7, 14.0^2 (SE 16% CV)

    # ---- Residual error -------------------------------------------------
    # Table S7 'Residual Intra-subject Variability [mmHg] = 6.02' with
    # footnote d 'expressed as square root of EPS', i.e. already an SD.
    addSd <- 6.02 ; label("Additive residual error on dSBP (mmHg)")  # Table S7 (SE 15% CV)
  })

  model({
    # ==================================================================
    # Covariate centering values. Song 2013 Equation 6 normalizes every
    # continuous covariate by "the median observations" of the analysis
    # dataset, but the exposure-response medians are never printed. The
    # values below are the MEANS of the pooled exposure-response dataset
    # from Table 2 ('All data' row, n = 4873), used as proxies. See the
    # vignette deviations section.
    # ==================================================================
    ref_sbp <- 165
    ref_age <- 54.8
    ref_wt  <- 94.9

    # ==================================================================
    # 1. Placebo effect (Song 2013 Eq. 1 first component).
    # ==================================================================
    pbo_study <- pbo_8635 * STUDY_CS8635_A_U301 +
      pbo_8663 * STUDY_CS8663_A_U301 +
      pbo_866318 * STUDY_866_318

    placebo_effect <- pbo_study *
      (SBP / ref_sbp)^e_sbp_pbo *
      (AGE / ref_age)^e_age_pbo *
      (1 + e_hispanic_pbo_8635 * RACE_HISPANIC * STUDY_CS8635_A_U301)

    # ==================================================================
    # 2. Monotherapy drug effects. Olmesartan and amlodipine are saturable
    #    (Eq. 5); hydrochlorothiazide is linear (Eq. 4) with the slope
    #    reported per 1000 ng*h/mL. Setting a component's AUC to 0 removes
    #    that component from the regimen.
    # ==================================================================
    mono_om <- emax_om * AUC_OLM / (AUC_OLM + eauc50_om) *
      (SBP / ref_sbp)^e_sbp_om *
      (1 + e_black_om * RACE_BLACK)

    mono_aml <- emax_aml * AUC_AML / (AUC_AML + eauc50_aml) *
      (SBP / ref_sbp)^e_sbp_aml *
      (WT / ref_wt)^e_wt_aml *
      (1 + e_sexf_aml * SEXF)

    mono_hctz <- slope_hctz * AUC_HCTZ / 1000 *
      (SBP / ref_sbp)^e_sbp_hctz

    # ==================================================================
    # 3. Interaction (Song 2013 Eq. 3).
    # ==================================================================
    interaction <- ia_om_aml * mono_om * mono_aml +
      ia_om_hctz * mono_om * mono_hctz +
      ia_aml_hctz * mono_aml * mono_hctz +
      ia_om_aml_hctz * mono_om * mono_aml * mono_hctz

    # ==================================================================
    # 4. Total response (Song 2013 Eq. 1 and 2) plus the additive
    #    inter-subject random effect. dsbp is the change from baseline in
    #    seated trough systolic blood pressure, in mmHg; negative values
    #    are blood-pressure lowering. Add the subject's baseline SBP to
    #    recover an absolute pressure.
    # ==================================================================
    drug_effect <- mono_om + mono_aml + mono_hctz + interaction

    dsbp <- placebo_effect + drug_effect + etadsbp
    dsbp ~ add(addSd)
  })
}
