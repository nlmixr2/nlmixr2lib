Song_2013_olmesartan_amlodipine_hydrochlorothiazide_dbp <- function() {
  description <- paste(
    "Steady-state exposure-response model for the change from baseline in seated trough DIASTOLIC blood pressure",
    "under the triple fixed-dose combination CS-8635 (olmesartan / amlodipine / hydrochlorothiazide) and its",
    "component mono- and dual-therapies, in adults with hypertension. The response is the sum of a",
    "study-specific placebo effect, three independent monotherapy drug effects driven by each component's",
    "steady-state AUC (saturable Emax for olmesartan and amlodipine, linear slope for hydrochlorothiazide), and",
    "a set of pairwise and three-way interaction terms proportional to the PRODUCT of the monotherapy effects,",
    "which make the combination sub-additive. Baseline diastolic pressure amplifies both the placebo and every",
    "drug effect; age amplifies the placebo effect but attenuates the olmesartan effect; Black race attenuates",
    "the olmesartan effect and (within CS8663-A-U301 only) the placebo effect; body weight attenuates the",
    "amlodipine effect. There are no ODE states and no dosing events -- exposure enters entirely through the",
    "AUC_OLM / AUC_AML / AUC_HCTZ covariate columns, which the companion population PK models generate as",
    "AUCss = Dose / (CL/F). Sister model files from the same paper:",
    "modellib('Song_2013_olmesartan_amlodipine_hydrochlorothiazide_sbp'), modellib('Song_2013_olmesartan'),",
    "modellib('Song_2013_amlodipine'), modellib('Song_2013_hydrochlorothiazide')."
  )
  reference <- paste(
    "Song S, Carrothers TJ, Moore H, Green M, Miller R, Rohatagi S, Lee J, Wang A, Melino M, Patel M,",
    "Heyrman R, Salazar DE.",
    "Model-supported development of CS-8635: a fixed-dose combination of olmesartan, amlodipine, and",
    "hydrochlorothiazide.",
    "Clin Pharmacol Drug Dev. 2013;2(2):103-112.",
    "doi:10.1002/cpdd.17.",
    "Structural equations from main-text Equations 1-7; parameter values from Supplemental Table S6, whose row",
    "labels supersede the Diastolic column of main-text Table 3 (see the model file comments and the vignette",
    "deviations section for the two mislabelled Table 3 rows).",
    sep = " "
  )
  vignette <- "Song_2013_olmesartan_amlodipine_hydrochlorothiazide"
  units <- list(
    time          = "h",
    dosing        = "n/a (no drug-dosing events; exposure enters as the covariates AUC_OLM, AUC_AML, AUC_HCTZ)",
    concentration = "change from baseline in seated trough diastolic blood pressure, in mmHg (not a drug concentration)"
  )

  # The additive inter-subject random effect sits on the TOTAL response
  # (NONMEM Y = placebo + drug + ETA + EPS), not on any single ini()
  # parameter, so it has no 'l<param>' partner to be named after.
  paper_specific_etas <- "etaddbp"

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
        "term emax_om * AUC_OLM / (AUC_OLM + eauc50_om) (Song 2013 Equation 5). Set AUC_OLM = 0 to remove",
        "olmesartan from the regimen: the olmesartan mono-effect and every interaction term containing it",
        "vanish, which is exactly how the AML/HCTZ dual-combination arm of Table 4 is simulated.",
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
        "emax_aml * AUC_AML / (AUC_AML + eauc50_aml). Song 2013 Discussion notes that the saturable form",
        "REPLACED the linear amlodipine relationship used in the group's earlier modelling, the enrichment of",
        "the dataset with CS8635-A-U301 having made the plateau identifiable.",
        "Worked typical value: 10 mg / 23.4 L/h = 427 ng*h/mL, which is close to eauc50_aml = 453, so the",
        "amlodipine arm sits near the steepest part of its own Emax curve at the highest marketed dose.",
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
        "slope_hctz * AUC_HCTZ / 1000 (Song 2013 Equation 4). The division by 1000 is mandatory: Table 3 and",
        "Table S6 report the slope in mmHg per 1000 ng*h/mL, so feeding the raw AUC to an undivided slope",
        "overstates the diuretic effect a thousand-fold.",
        "Worked typical value: 25 mg / 20.3 L/h = 1232 ng*h/mL.",
        "Set AUC_HCTZ = 0 to remove hydrochlorothiazide from the regimen."
      ),
      source_name        = "AUCHCTZ"
    ),
    DBP = list(
      description        = "Baseline seated trough diastolic blood pressure.",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (pre-treatment baseline). Measured after 5 minutes seated, as the mean of three",
        "readings (two in study 866-318) taken at 1-minute intervals at approximately the same time of day at",
        "each visit. Enters in power form on the placebo term, (DBP / 102)^3.19, and on the olmesartan and",
        "amlodipine drug-effect terms, (DBP / 102)^2.46 and (DBP / 102)^4.12 -- all three positive, so",
        "higher-baseline subjects show larger placebo AND larger treatment effects, the paper's second key",
        "covariate finding. Table S6 reports NO baseline effect on the hydrochlorothiazide term for diastolic",
        "pressure (contrast the systolic sibling, which has one). The centering value 102 mmHg is the mean",
        "baseline DBP of the pooled exposure-response dataset (Song 2013 Table 2, 'All data' row, n = 4873),",
        "used as a proxy for the median that Song 2013 specifies in Equation 6 but never prints."
      ),
      source_name        = "Baseline BP"
    ),
    AGE = list(
      description        = "Subject age.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Enters in power form on the placebo term, (AGE / 54.8)^1.37, and on the",
        "olmesartan drug-effect term, (AGE / 54.8)^-0.818 -- the two have OPPOSITE signs, so older subjects",
        "show a larger placebo response but a smaller olmesartan response. The centering value 54.8 years is the",
        "mean age of the pooled exposure-response dataset (Song 2013 Table 2, 'All data' row).",
        "The age effect on the olmesartan term is the parameter whose row label differs between the two source",
        "tables: Supplemental Table S6 labels -0.818 'Effect of age on Drug Effect of OM' while main-text",
        "Table 3 prints the same value under 'Effect of baseline BP on drug effect of HCTZ'. Table S6 is used",
        "here because the Supplemental Results narrative states the diastolic age effect explicitly for",
        "olmesartan ('diastolic blood pressure lowering effects might be attenuated in subjects with advanced",
        "age') and confines the hydrochlorothiazide baseline effect to SYSTOLIC pressure, and because a negative",
        "baseline exponent on the hydrochlorothiazide term would contradict the paper's own repeated statement",
        "that higher-baseline subjects show larger treatment effects for all three components."
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
        "(WT / 94.9)^-0.830, i.e. 'Treatment effects of AML were lower in subjects with higher body weights'",
        "(Song 2013 Results). The centering value 94.9 kg is the mean body weight of the pooled",
        "exposure-response dataset (Song 2013 Table 2, 'All data' row), used as a proxy for the unprinted median."
      ),
      source_name        = "Weight"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator, 1 = Black, 0 = other.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black)",
      notes              = paste(
        "Time-fixed per subject. Enters additively-fractionally (Song 2013 Equation 7,",
        "Response x (1 + coefficient x Covariate)) in two places: on the olmesartan drug effect for all",
        "subjects, factor (1 - 0.263), and on the PLACEBO effect for subjects in CS8663-A-U301 only, factor",
        "(1 - 0.607). The olmesartan effect is the paper's first key covariate finding -- Black subjects have",
        "less renin-angiotensin-system activation and respond less well to angiotensin receptor blockers.",
        "Note that Song 2013's Results prose quotes '39.6% and 29.3% less OM treatment effect on SeDBP and",
        "SeSBP' while the Supplemental Results quotes '26 and 50 %' for the same two endpoints; both parameter",
        "tables agree on -0.263 (diastolic) and -0.393 (systolic), so the tabulated values are used and the",
        "prose percentages are treated as erroneous."
      ),
      source_name        = "Black race"
    ),
    STUDY_CS8635_A_U301 = list(
      description        = "Indicator for enrolment in study CS8635-A-U301 (TRINITY), the pivotal phase III CS-8635 triple-combination trial.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the other two pooled phase III studies)",
      notes              = paste(
        "Selects the study-specific placebo effect of -3.80 mmHg. The three study indicators partition the",
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
        "Selects the study-specific placebo effect of -3.57 mmHg AND gates the Black-race effect on the placebo",
        "term, which Song 2013 retained only within this study. n = 1920 in the exposure-response dataset."
      ),
      source_name        = "Study CS8663-A-U301"
    ),
    STUDY_866_318 = list(
      description        = "Indicator for enrolment in study 866-318 (written SE866-318 in the parameter tables), the phase III olmesartan / hydrochlorothiazide dual-combination trial.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (the other two pooled phase III studies)",
      notes              = paste(
        "Selects the study-specific placebo effect of -6.08 mmHg, which is roughly 70% larger in magnitude than",
        "in either of the other two studies. n = 495 in the exposure-response dataset; this study also had the",
        "highest mean baseline DBP (104 mmHg) and took only two rather than three readings per visit."
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
      "exposure-response dataset -- it contributed only to the pre-TRINITY predictive modelling and had no PK",
      "samples. The primary endpoint visit was Week 12 Period II for CS8635-A-U301 and Week 8 Period II for the",
      "other two, by which time therapeutic steady state had been reached. Estimation used FOCE in NONMEM VI",
      "level 1."
    )
  )

  ini({
    # ==================================================================
    # All values from Song 2013 Supplemental Table S6 ('Parameter
    # Estimates for Final Exposure-Response Model for DBP'). Table S6 is
    # preferred over the Diastolic column of main-text Table 3 because two
    # Table 3 rows carry the right values under the wrong labels; see the
    # note on e_age_om below and the vignette deviations section.
    # ==================================================================

    # ---- Placebo effect, one value per contributing study --------------
    pbo_8635    <- -3.80 ; label("Placebo effect on dDBP in study CS8635-A-U301 (mmHg)")  # Table S6 (SE 13% CV)
    pbo_8663    <- -3.57 ; label("Placebo effect on dDBP in study CS8663-A-U301 (mmHg)")  # Table S6 (SE 13% CV)
    pbo_866318  <- -6.08 ; label("Placebo effect on dDBP in study 866-318 (mmHg)")        # Table S6 (SE 10% CV)

    # ---- Covariate effects on the placebo term -------------------------
    e_dbp_pbo         <-  3.19  ; label("Power exponent for baseline DBP on the placebo effect (unitless)")                       # Table S6 'Effect of baseline on Placebo Effect' (SE 24% CV)
    e_age_pbo         <-  1.37  ; label("Power exponent for age on the placebo effect (unitless)")                                # Table S6 'Effect of age on Placebo Effect' (SE 21% CV)
    e_black_pbo_8663  <- -0.607 ; label("Fractional Black-race effect on the placebo effect within CS8663-A-U301 (unitless)")     # Table S6 'Effect of Black Race on Placebo Effect in study CS8663-A-U301' (SE 23% CV)

    # ---- Olmesartan monotherapy effect (saturable, Eq. 5) --------------
    emax_om    <- -10.5 ; label("Maximum olmesartan effect on dDBP (mmHg)")                                     # Table S6 'Emax, OM' (SE 11% CV)
    eauc50_om  <-  1850 ; label("Olmesartan AUCss producing half the maximum dDBP effect (ng*h/mL)")            # Table S6 'EAUC50, OM' (SE 31% CV)
    e_dbp_om   <-  2.46 ; label("Power exponent for baseline DBP on the olmesartan effect (unitless)")          # Table S6 'Effect of baseline on Drug Effect of OM' (SE 18% CV)
    # Table S6 labels this row 'Effect of age on Drug Effect of OM'; main-text Table 3 prints the same -0.818
    # under 'Effect of baseline BP on drug effect of HCTZ'. Table S6 wins: the Supplemental Results narrative
    # states the diastolic age effect for olmesartan explicitly and restricts the HCTZ baseline effect to
    # systolic pressure, and a NEGATIVE baseline exponent would contradict the paper's own statement that
    # higher-baseline subjects show larger treatment effects for all three components.
    e_age_om   <- -0.818 ; label("Power exponent for age on the olmesartan effect (unitless)")                  # Table S6 'Effect of age on Drug Effect of OM' (SE 17% CV)
    e_black_om <- -0.263 ; label("Fractional Black-race effect on the olmesartan effect (unitless)")            # Table S6 'Effect of Black Race on Drug Effect of OM' (SE 24% CV)

    # ---- Amlodipine monotherapy effect (saturable, Eq. 5) --------------
    emax_aml   <- -19.3 ; label("Maximum amlodipine effect on dDBP (mmHg)")                                     # Table S6 'Emax, AML' (SE 21% CV)
    eauc50_aml <-   453 ; label("Amlodipine AUCss producing half the maximum dDBP effect (ng*h/mL)")            # Table S6 'EAUC50, AML' (SE 40% CV)
    e_dbp_aml  <-  4.12 ; label("Power exponent for baseline DBP on the amlodipine effect (unitless)")          # Table S6 'Effect of baseline on Drug Effect of AML' (SE 13% CV)
    e_wt_aml   <- -0.830 ; label("Power exponent for body weight on the amlodipine effect (unitless)")          # Table S6 'Effect of weight on Drug Effect of AML' (SE 15% CV)

    # ---- Hydrochlorothiazide monotherapy effect (linear, Eq. 4) --------
    slope_hctz <- -3.3 ; label("Hydrochlorothiazide slope on dDBP (mmHg per 1000 ng*h/mL)")                     # Table S6 'Slope HCTZ' (SE 9% CV)

    # ---- Interaction terms (Eq. 3), proportional to the PRODUCT of the
    #      monotherapy effects. The monotherapy effects are negative, so a
    #      positive pairwise coefficient makes the pair SUB-additive, while
    #      the positive three-way coefficient acting on a product of three
    #      negatives adds a little lowering back.
    ia_om_aml      <- 0.043   ; label("Olmesartan x amlodipine interaction on dDBP (1/mmHg)")                   # Table S6 'Interaction OM*AML' (SE 10% CV)
    ia_aml_hctz    <- 0.0747  ; label("Amlodipine x hydrochlorothiazide interaction on dDBP (1/mmHg)")          # Table S6 'Interaction HCTZ*AML' (SE 13% CV)
    # Table S6 reports the OM x HCTZ pairwise term as 'n.s.' (main Table 3: 'NA'), i.e. it was tested and
    # dropped from the diastolic model. It is carried here as fixed(0) rather than deleted so that the
    # diastolic and systolic siblings remain structurally identical and the absence is explicit.
    ia_om_hctz     <- fixed(0); label("Olmesartan x hydrochlorothiazide interaction on dDBP (1/mmHg)")          # Table S6 'Interaction OM*HCTZ' = n.s.; not retained in the diastolic model
    ia_om_aml_hctz <- 0.00512 ; label("Olmesartan x amlodipine x hydrochlorothiazide interaction on dDBP (1/mmHg^2)")  # Table S6 'Interaction OM*HCTZ*AML' (SE 29% CV)

    # ---- Inter-subject variability -------------------------------------
    # Table S6 'Additive Inter-subject Variability [mmHg] = 8.56' with
    # footnote b 'Square root of ETA estimate', so the ETA VARIANCE is
    # 8.56^2 = 73.2736 mmHg^2. The random effect is ADDITIVE on the
    # response, matching the NONMEM form Y = placebo + drug + ETA + EPS.
    etaddbp ~ 73.2736  # Table S6, 8.56^2 (SE 16% CV)

    # ---- Residual error -------------------------------------------------
    # Table S6 'Residual Intra-subject Variability [mmHg] = 3.62' with
    # footnote d 'expressed as square root of EPS', i.e. already an SD.
    addSd <- 3.62 ; label("Additive residual error on dDBP (mmHg)")  # Table S6 (SE 18% CV)
  })

  model({
    # ==================================================================
    # Covariate centering values. Song 2013 Equation 6 normalizes every
    # continuous covariate by "the median observations" of the analysis
    # dataset, but -- unlike the population PK equations 8-10, which print
    # their centering values -- the exposure-response medians are never
    # printed. The values below are the MEANS of the pooled
    # exposure-response dataset from Table 2 ('All data' row, n = 4873),
    # used as proxies. See the vignette deviations section; the proxy is
    # calibrated against the PK equations, where substituting the
    # computable dataset mean for the printed median agrees to 1-3%.
    # ==================================================================
    ref_dbp <- 102
    ref_age <- 54.8
    ref_wt  <- 94.9

    # ==================================================================
    # 1. Placebo effect (Song 2013 Eq. 1 first component). A separate
    #    numerical value per study, then scaled by the continuous
    #    covariates in power form (Eq. 6) and by the categorical covariate
    #    in fractional form (Eq. 7).
    # ==================================================================
    pbo_study <- pbo_8635 * STUDY_CS8635_A_U301 +
      pbo_8663 * STUDY_CS8663_A_U301 +
      pbo_866318 * STUDY_866_318

    placebo_effect <- pbo_study *
      (DBP / ref_dbp)^e_dbp_pbo *
      (AGE / ref_age)^e_age_pbo *
      (1 + e_black_pbo_8663 * RACE_BLACK * STUDY_CS8663_A_U301)

    # ==================================================================
    # 2. Monotherapy drug effects. Olmesartan and amlodipine are saturable
    #    (Eq. 5, Emax * AUC / (AUC + EAUC50)); hydrochlorothiazide is
    #    linear (Eq. 4, Slope * AUC), with the slope reported per
    #    1000 ng*h/mL. Setting a component's AUC to 0 removes that
    #    component from the regimen.
    # ==================================================================
    mono_om <- emax_om * AUC_OLM / (AUC_OLM + eauc50_om) *
      (DBP / ref_dbp)^e_dbp_om *
      (AGE / ref_age)^e_age_om *
      (1 + e_black_om * RACE_BLACK)

    mono_aml <- emax_aml * AUC_AML / (AUC_AML + eauc50_aml) *
      (DBP / ref_dbp)^e_dbp_aml *
      (WT / ref_wt)^e_wt_aml

    mono_hctz <- slope_hctz * AUC_HCTZ / 1000

    # ==================================================================
    # 3. Interaction (Song 2013 Eq. 3): each term is a coefficient times
    #    the product of the monotherapy effects it couples.
    # ==================================================================
    interaction <- ia_om_aml * mono_om * mono_aml +
      ia_om_hctz * mono_om * mono_hctz +
      ia_aml_hctz * mono_aml * mono_hctz +
      ia_om_aml_hctz * mono_om * mono_aml * mono_hctz

    # ==================================================================
    # 4. Total response (Song 2013 Eq. 1 and 2) plus the additive
    #    inter-subject random effect. ddbp is the change from baseline in
    #    seated trough diastolic blood pressure, in mmHg; negative values
    #    are blood-pressure lowering. Add the subject's baseline DBP to
    #    recover an absolute pressure.
    # ==================================================================
    drug_effect <- mono_om + mono_aml + mono_hctz + interaction

    ddbp <- placebo_effect + drug_effect + etaddbp
    ddbp ~ add(addSd)
  })
}
