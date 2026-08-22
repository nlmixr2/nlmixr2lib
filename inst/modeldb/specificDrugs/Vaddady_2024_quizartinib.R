Vaddady_2024_quizartinib <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral quizartinib and",
    "its pharmacologically active metabolite AC886 in adults, pooled across",
    "13 studies in healthy volunteers, subjects with hepatic impairment,",
    "relapsed/refractory (R/R) FLT3-ITD-positive acute myeloid leukemia",
    "(AML) patients on quizartinib monotherapy, and newly diagnosed",
    "FLT3-ITD-positive AML patients receiving quizartinib alongside standard",
    "cytarabine-anthracycline induction and cytarabine consolidation",
    "chemotherapy followed by single-agent continuation (QuANTUM-First;",
    "Vaddady 2024). Quizartinib is described by a three-compartment model",
    "with sequential zero-order (duration D1) then first-order (ka)",
    "absorption from a depot, an absorption lag time, and first-order",
    "elimination from the central compartment. AC886 is described by a",
    "two-compartment model with first-order formation from the quizartinib",
    "central compartment and first-order elimination. Because no",
    "intravenous data were available, the parent-to-metabolite conversion",
    "fraction fMET is not identifiable and is fixed at 0.5; following the",
    "authors' parameterisation, the entire parent elimination flux is",
    "routed into the metabolite compartment and all four AC886 disposition",
    "parameters are divided by fMET, which is algebraically equivalent and",
    "leaves fMET acting as a pure multiplicative scale on predicted AC886",
    "concentrations. Retained covariates are allometric body-weight scaling",
    "on all clearances and volumes of both moieties (exponents fixed at",
    "0.75 and 1), strong CYP3A inhibitor coadministration on quizartinib CL",
    "and Frel, moderate CYP3A inhibitor coadministration on Frel, non-AML",
    "subject status on Frel and ka, Black race on quizartinib CL, female sex",
    "on quizartinib Vc, age on the first peripheral volume, and, for newly",
    "diagnosed AML patients only, a treatment-phase effect (induction /",
    "consolidation / continuation) on Frel. AC886 covariates are non-AML",
    "status, Black race and strong CYP3A inhibitors on CL, strong CYP3A",
    "inhibitors on Vc, and the same treatment-phase effect on fMET.",
    "Interindividual variability is reported on quizartinib CL, Q1, Vc,",
    "Tlag, ka, D1 and Frel (with separate Frel variances for AML patients,",
    "whose eta is Box-Cox transformed, and non-AML subjects) and on AC886",
    "CL (separate variances for AML patients and non-AML subjects) and Vc.",
    "Residual error is additive on the log scale, i.e. proportional on the",
    "linear concentration scale, separately for each moiety.")
  reference <- paste(
    "Vaddady P, Glatard A, Smania G, Nakayama S, Inoue H, Kurumaddali A,",
    "Abutarif M, Zheng M. Population pharmacokinetic analysis of",
    "quizartinib in patients with newly diagnosed",
    "FLT3-internal-tandem-duplication-positive acute myeloid leukemia.",
    "Clin Transl Sci. 2024;17:e70074. doi:10.1111/cts.70074")
  vignette <- "Vaddady_2024_quizartinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Amounts are carried in mg (the dosing unit) for both
  # moieties; the AC886 formation flux therefore carries the AC886/quizartinib
  # molecular-weight ratio so that mass is conserved on a molar basis (see
  # model() and Vaddady 2024 Methods, "Model development": molecular weights
  # 560.68 g/mol for quizartinib and 576.67 g/mol for AC886).
  compartmentData <- list(
    depot              = list(analyte = "quizartinib", units = "mg", specimen = "administration site", verified = TRUE),
    central            = list(analyte = "quizartinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1        = list(analyte = "quizartinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2        = list(analyte = "quizartinib", units = "mg", specimen = "plasma", verified = TRUE),
    central_ac886      = list(analyte = "AC886", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ac886  = list(analyte = "AC886", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (baseline).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Mechanistic (allometric) covariate included without statistical testing. Vaddady 2024 Equation 6: COVEff_WT = (WT/75)^theta_WT with theta_WT FIXED to 0.75 for all clearances (CL_quiz, Q1_quiz, Q2_quiz, CL_AC886, Q_AC886) and to 1 for all volumes (Vc_quiz, Vp1_quiz, Vp2_quiz, Vc_AC886, Vp_AC886). The reference weight of 75 kg is printed in Equation 6 itself and repeated in the Figure 4 reference subject; note that the Table 2 median weight is 72.0 kg, so 75 kg is a rounded reference rather than the cohort median.",
      source_name        = "WTKGBL"
    ),
    AGE = list(
      description        = "Subject age at baseline.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Exponential (log-linear) effect on the first peripheral volume of quizartinib only. Vaddady 2024 Equation 13: Vp1_quiz = 312 * (WT/75)^1 * exp(0.0152 * (AGE - 47)). The centering value of 47 years is printed in Equation 13 and in the Figure 4 reference subject; the Table 2 median age is 50.0 years, so 47 years is the model's centering constant rather than the cohort median. Selected by the stepwise covariate model (SCM) procedure.",
      source_name        = "AGEYBL"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male; the most common category, 495/932 = 53.1% per Table 2).",
      notes              = "Fractional change of -0.169 on the quizartinib central volume for females relative to the male reference (Vaddady 2024 Table 3 and Equation 11: Vc_quiz = 371 * (WT/75)^1 * (1 - 0.169 for females)). The source NONMEM column SEXN is coded 1 = male (reference) and 2 = female; the canonical SEXF column carries the female indicator directly, so no sign inversion is needed. Selected by the SCM procedure.",
      source_name        = "SEXN"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator, 1 = Black or African American, 0 = any other race.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 ('rest of the population' per Vaddady 2024 Equations 8 and 17; White 65.4%, Asian 18.1%, Other 4.8%, Black 8.5% per Table 2).",
      notes              = "Fractional change of -0.261 on quizartinib CL (Equation 8) and +0.488 on AC886 CL (Equation 17). Both were selected by the SCM procedure. The quizartinib CL effect carries the largest uncertainty of any retained effect in the model (RSE 38.0%, the only RSE above 30%), and Vaddady 2024 Results notes the Black-race effect 'was associated with large uncertainty'; the effect on total (quizartinib + AC886) exposure largely cancels between the two moieties (Figure 4c). Missing race (29/932 = 3.1%) was imputed to the most common category, i.e. RACE_BLACK = 0 (source stream: IF(RACE3N.EQ.-99) CLRACE3N = 1).",
      source_name        = "RACE3N"
    ),
    DIS_AML = list(
      description        = "Acute myeloid leukemia patient indicator, 1 = AML patient (either relapsed/refractory or newly diagnosed), 0 = non-AML subject (healthy volunteer or subject with hepatic impairment).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (AML patients) -- NOTE the reference is the AML group, which is the inverse of this canonical column's usual 0-is-reference convention. Vaddady 2024 Methods: 'The AML patient group was set as the reference category versus which the effect of non-AML subjects was estimated.' AML patients are 659/932 (70.7%: 365 R/R + 294 newly diagnosed) and non-AML subjects 273/932 (29.3%) per Table 2.",
      notes              = "The source NONMEM column AML3 is the NON-AML indicator (AML3 = 1 for healthy volunteers and hepatic-impairment subjects). The canonical DIS_AML column carries the AML indicator, so every effect below is applied via the complement (1 - DIS_AML); this is a pure relabelling with no sign change to any reported coefficient, following the (1 - SEXF) precedent in Lahu_2010_roflumilast.R. Effects: Frel = 1.73 for non-AML subjects (Equation 9, a direct multiplier rather than a fractional change), fractional change -0.188 on ka (Equation 10), and +0.843 on AC886 CL (Equation 17). This column ALSO selects which stratum-specific interindividual-variability term applies: AML patients take etalfdepot_aml (Box-Cox transformed) and etalcl_ac886_aml, while non-AML subjects take etalfdepot_nonaml and etalcl_ac886_nonaml. Vaddady 2024 additionally reports separate residual-error magnitudes for the two groups which nlmixr2 cannot express as a covariate-switched residual (see vignette Errata).",
      source_name        = "AML3"
    ),
    CONMED_CYP3A4_INH_STRONG = list(
      description        = "Concomitant strong CYP3A inhibitor coadministration indicator, 1 = a strong CYP3A inhibitor is coadministered over the record, 0 = not coadministered.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no strong CYP3A inhibitor; 487/932 = 52.3% of subjects received no CYP3A inhibitor and 184/932 = 19.7% received a strong inhibitor per Table 2).",
      notes              = "Mechanistic covariate included without statistical testing. Time-varying: Vaddady 2024 Methods states 'concomitant medications were evaluated as covariates varying over time'. Fractional changes: -0.301 on quizartinib CL and +0.273 on quizartinib Frel (Equations 8 and 9), +0.298 on AC886 CL and +2.79 on AC886 Vc (Equations 17 and 18). This is the covariate with the largest impact on quizartinib exposure, increasing steady-state AUC roughly 1.8-fold, and is the basis for the labelled dose adjustment during strong CYP3A inhibitor coadministration.",
      source_name        = "CYPINH3"
    ),
    CONMED_CYP3A4_INH_MOD = list(
      description        = "Concomitant moderate CYP3A inhibitor coadministration indicator, 1 = a moderate CYP3A inhibitor is coadministered over the record, 0 = not coadministered.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no moderate CYP3A inhibitor; 176/932 = 18.9% of subjects received a moderate inhibitor per Table 2).",
      notes              = "Time-varying, as for the strong-inhibitor column. Fractional change of +0.116 on quizartinib relative bioavailability only (Vaddady 2024 Equation 9); selected by the SCM procedure. No moderate-inhibitor effect was retained on quizartinib CL or on either AC886 parameter.",
      source_name        = "CYPINH2"
    ),
    TRTPH_INDUCTION = list(
      description        = "QuANTUM-First induction treatment-phase indicator, 1 = the record falls in the induction phase of a newly diagnosed AML patient, 0 = any other phase or population. Member of the TRTPH_<phase> canonical family.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0. All three TRTPH_ indicators are 0 for the reference group, which is R/R AML patients (and, for the fMET effect, non-AML subjects as well) -- NOT the induction phase. See notes.",
      notes              = "Protocol: QuANTUM-First (AC220-A-U302), a Phase 3 trial in newly diagnosed FLT3-ITD-positive AML. Induction = quizartinib/placebo plus intravenous cytarabine and anthracycline (daunorubicin or idarubicin). Phase-to-column mapping: induction -> TRTPH_INDUCTION, consolidation (high-dose cytarabine plus quizartinib/placebo) -> TRTPH_CONSOLIDATION, continuation (single-agent quizartinib/placebo for up to 3 years) -> TRTPH_CONTINUATION. This differs from the TRTPH_ family default in one respect: the family preamble makes the protocol's own induction phase the all-indicators-0 reference, but Vaddady 2024 sets the reference to R/R AML patients receiving quizartinib monotherapy, 'for whom no distinct treatment phases were reported'. All three phases therefore need their own indicator. Time-varying per record. Fractional change -0.419 on quizartinib Frel (Equation 9) and +0.715 on fMET (Equation 16). The Frel and fMET phase effects act in opposite directions by construction: the parent Frel effect is carried over to the metabolite in parent-metabolite modelling, and the fMET phase effect exists to counterbalance that spurious carry-over (Vaddady 2024 Discussion).",
      source_name        = "PHASE == 1"
    ),
    TRTPH_CONSOLIDATION = list(
      description        = "QuANTUM-First consolidation treatment-phase indicator, 1 = the record falls in the consolidation phase of a newly diagnosed AML patient, 0 = any other phase or population. Member of the TRTPH_<phase> canonical family.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (reference = R/R AML patients; see TRTPH_INDUCTION notes).",
      notes              = "Protocol: QuANTUM-First (AC220-A-U302). Consolidation = standard high-dose cytarabine consolidation plus quizartinib/placebo, allogeneic hematopoietic cell transplantation, or both, in patients achieving remission. Fractional change -0.192 on quizartinib Frel (Equation 9) and +0.272 on fMET (Equation 16). Vaddady 2024 Discussion notes consolidation is the phase in which newly diagnosed patients behaved most similarly to the R/R reference, because the time elapsed since initiation of chemotherapy is comparable between the two groups.",
      source_name        = "PHASE == 2"
    ),
    TRTPH_CONTINUATION = list(
      description        = "QuANTUM-First continuation treatment-phase indicator, 1 = the record falls in the continuation phase of a newly diagnosed AML patient, 0 = any other phase or population. Member of the TRTPH_<phase> canonical family.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (reference = R/R AML patients; see TRTPH_INDUCTION notes).",
      notes              = "Protocol: QuANTUM-First (AC220-A-U302). Continuation = single-agent quizartinib/placebo for up to 3 years in patients with blood count recovery, i.e. no background chemotherapy. Fractional change +0.418 on quizartinib Frel (Equation 9) and -0.249 on fMET (Equation 16). This is the phase with the highest dose-normalised quizartinib exposure (about 1.4-fold the R/R reference); Vaddady 2024 confirmed by a VPC restricted to patients who entered continuation that the effect is phase-related rather than driven by selection of subjects who survived to that phase.",
      source_name        = "PHASE == 3"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese-heritage race indicator, 1 = Japanese, 0 = otherwise.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a structural covariate on quizartinib CL, Vc, Vp and Frel (Table 2). A Frel effect was carried in the final NONMEM control stream as THETA(17) but is FIXED to exactly 0 there (`0 FIX ; 17_F1JAP1`), i.e. the multiplier is identically 1 and the covariate has no effect in the final model; it is correspondingly absent from Table 3. Documented here so the provenance of the screen is preserved without carrying an inert covariate into model()."
    ),
    ALB = list(
      description = "Serum albumin at baseline.",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Tested on AC886 Vc (structural) and on quizartinib and AC886 CL (exploratory) per Table 2. Carried in the final AC886 control stream as THETA(33) centred on 4.1 g/dL but FIXED to 0 (`0 FIX ; 33_V5ALBUBL1`), so it has no effect in the final model and is absent from Table 3."
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A inducer coadministration indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Weak (36/932) and moderate (1/932) CYP3A inducer use was tested on quizartinib CL and Frel (Table 2) but not retained in the final model. The final model's predictive performance in inducer-treated subjects was nevertheless confirmed by a stratified pcVPC (Vaddady 2024 Figure 3b)."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested as an exploratory covariate on AC886 CL only (Table 2); not retained. Quizartinib and AC886 are eliminated by hepatic metabolism, and Vaddady 2024 Conclusions state no dose adjustment is needed by renal function."
    ),
    NCI_ODWG = list(
      description = "National Cancer Institute Organ Dysfunction Working Group hepatic-function grade (normal / mild / moderate).",
      units       = "(category)",
      type        = "categorical",
      notes       = "Tested as an exploratory covariate on quizartinib and AC886 CL (Table 2); not retained. Liver function tests (alkaline phosphatase, alanine aminotransferase, aspartate aminotransferase, total bilirubin) were likewise tested on CL and not retained, despite the analysis pooling a dedicated hepatic-impairment study (AC220-016 and AC220-A-U105)."
    ),
    FORM_QUIZ_SOLUTION = list(
      description = "Oral-solution formulation indicator, 1 = oral solution, 0 = tablet.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Tested as a structural covariate on quizartinib ka and Frel (Table 2; tablet 805/932 = 86.4%, solution 127/932 = 13.6%); not retained in the final model. Antacid and proton-pump-inhibitor coadministration were also tested on ka and Frel and not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 932L,
    n_studies      = 13L,
    age_range      = "18-91 years (median 50.0); the QuANTUM-First subset enrolled patients aged 20-75 years",
    age_median     = "50.0 years",
    weight_range   = "36.8-153 kg (median 72.0)",
    weight_median  = "72.0 kg",
    sex_female_pct = 46.9,
    race_ethnicity = c(White = 65.4, Black = 8.5, Asian = 18.1, Other = 4.8),
    disease_state  = "Pooled: 294 newly diagnosed FLT3-ITD-positive AML patients (31.5%), 365 relapsed/refractory FLT3-ITD-positive AML patients (39.2%), and 273 non-AML subjects (29.3%) comprising healthy volunteers and subjects with hepatic impairment.",
    dose_range     = "20-90 mg/day oral quizartinib as tablet (805/932 subjects) or oral solution (127/932), given as single doses (30-90 mg) or multiple daily doses (20-90 mg)",
    regions        = "Multinational; three of the 13 studies enrolled Japanese patients (AC220-A-J101, AC220-A-J102, AC220-A-J201) and 75/932 subjects (8.0%) were Japanese.",
    n_observations = "14,160 quizartinib and 13,399 AC886 plasma concentrations. Quantified by two cross-validated LC-MS/MS methods with lower limits of quantification of 2 ng/mL and 0.5 ng/mL.",
    hepatic_function = "Includes a dedicated Child-Pugh hepatic-impairment study (AC220-016, 30 subjects) and an NCI-ODWG-criteria hepatic-impairment study (AC220-A-U105, 12 subjects). NCI-ODWG grade in the pooled analysis set: normal 763/932 (81.9%), mild 139/932 (14.9%), moderate 16/932 (1.7%).",
    co_medication  = "CYP3A inhibitors: none 487/932 (52.3%), weak 85/932 (9.1%), moderate 176/932 (18.9%), strong 184/932 (19.7%). CYP3A inducers: none 895/932 (96.0%), weak 36/932 (3.9%), moderate 1/932 (0.1%). Antacids 90/932 (9.7%); proton-pump inhibitors 327/932 (35.1%). Newly diagnosed patients received background cytarabine plus daunorubicin or idarubicin during induction and high-dose cytarabine during consolidation.",
    notes          = "Baseline demographics are Vaddady 2024 Table 2; the study inventory (13 studies: nine Phase 1, two Phase 2, two Phase 3, including the Phase 3 QuANTUM-First trial AC220-A-U302 and the Phase 3 R/R trial QuANTUM-R / AC220-007) is Table 1. Seven of the 13 studies were also used in the earlier healthy-volunteer / R/R AML analysis. Missing covariates were imputed with the median (continuous) or most common category (categorical). Fitted in NONMEM 7.4.4 with FOCE+I; the quizartinib model used ADVAN12 TRANS4 and the AC886 model a six-compartment ADVAN13. The AC886 model was fitted sequentially, with individual quizartinib parameters fixed to their empirical Bayes estimates from the final quizartinib model; this model file expresses the two moieties as one joint hierarchical model so that both can be simulated together."
  )

  ini({
    # ------------------------------------------------------------------
    # QUIZARTINIB (PARENT) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # All typical values are Vaddady 2024 Table 3, "Final quizartinib
    # model" column, and are reproduced in Supporting Information
    # Equations 8-15. The final AC886 control stream (Supporting
    # Information, NONMEM code) re-states each of these as a FIX row
    # carrying more decimal places, which independently confirms that the
    # Table 3 numbers are final estimates rather than initial values (the
    # non-FIX $THETA rows in both streams are initial estimates and
    # differ, e.g. CL initial 6.85443 vs final 6.6545).
    # Reference subject: 75 kg, 47 years old, male, R/R AML patient, no
    # CYP3A coadministration, not Black (Vaddady 2024 Figure 4 caption).

    lcl <- log(6.65)
    label("Quizartinib apparent clearance CL/F at the reference subject (L/h)")   # Table 3: CL_quiz = 6.65 L/h, RSE 1.61%; Eq 8; AC886 stream FIX row 1_CL = 6.6545

    lvc <- log(371)
    label("Quizartinib apparent central volume Vc/F at the reference subject (L)")  # Table 3: Vc,quiz = 371 L, RSE 3.13%; Eq 11; AC886 stream FIX row 2_V2 = 371.079

    lq <- log(40.7)
    label("Quizartinib apparent intercompartmental clearance Q1/F to peripheral1 (L/h)")  # Table 3: Q1,quiz = 40.7 L/h, RSE 4.68%; Eq 12; AC886 stream FIX row 3_Q3 = 40.7235

    lvp <- log(312)
    label("Quizartinib apparent first peripheral volume Vp1/F at the reference subject (L)")  # Table 3: Vp1,quiz = 312 L, RSE 1.95%; Eq 13; AC886 stream FIX row 4_V3 = 312.336

    lq2 <- log(0.757)
    label("Quizartinib apparent intercompartmental clearance Q2/F to peripheral2 (L/h)")  # Table 3: Q2,quiz = 0.757 L/h, RSE 3.89%; Eq 14; AC886 stream FIX row 5_Q4 = 0.757052

    lvp2 <- log(91.9)
    label("Quizartinib apparent second peripheral volume Vp2/F (L)")             # Table 3: Vp2,quiz = 91.9 L, RSE 2.24%; Eq 15; AC886 stream FIX row 6_V4 = 91.9306

    ltlag <- log(0.196)
    label("Quizartinib absorption lag time Tlag (h)")                            # Table 3: Tlag = 0.196 h, RSE 4.24%; AC886 stream FIX row 7_ALAG1 = 0.196137

    lka <- log(1.10)
    label("Quizartinib first-order absorption rate constant ka at the AML reference (1/h)")  # Table 3: ka = 1.10 1/h, RSE 6.28%; Eq 10; AC886 stream FIX row 8_KA = 1.09723

    ld1 <- log(0.710)
    label("Quizartinib duration of zero-order input into the depot D1 (h)")       # Table 3: D1 = 0.710 h, RSE 6.37%; AC886 stream FIX row 9_D1 = 0.709816

    lfdepot <- fixed(log(1))
    label("Quizartinib relative bioavailability Frel at the R/R AML reference (unitless, non-identifiable anchor)")  # Eq 9: Frel = 1 * (covariate terms), i.e. the R/R AML reference level is fixed at 1; absolute F is not identifiable without intravenous data

    # ------------------------------------------------------------------
    # ALLOMETRIC EXPONENTS (mechanistic, fixed)
    # ------------------------------------------------------------------
    # Vaddady 2024 Equation 6: COVEff_WT = (WT/75)^theta_WT, "where
    # theta_WT was fixed to 0.75 for clearance (CL) and 1 for the volume
    # of distribution (V) parameters". Body weight is a mechanistic
    # covariate included in the base model without statistical testing.

    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on (WT/75) for all clearances of both moieties (unitless)")  # Eq 6, fixed per allometric theory (refs 21, 22)

    e_wt_vc <- fixed(1)
    label("Allometric exponent on (WT/75) for all volumes of both moieties (unitless)")    # Eq 6, fixed per allometric theory (refs 21, 22)

    # ------------------------------------------------------------------
    # QUIZARTINIB COVARIATE EFFECTS
    # ------------------------------------------------------------------
    # Categorical effects follow Vaddady 2024 Equation 5: the multiplier
    # is 1 at the reference category and (1 + theta) otherwise. The
    # continuous AGE effect follows Equation 4: exp(theta * (Cov - Covref)).

    e_cyp3a4_inh_strong_cl <- -0.301
    label("Fractional change in quizartinib CL with strong CYP3A inhibitors (unitless)")  # Table 3, RSE 7.38%; Eq 8

    e_race_black_cl <- -0.261
    label("Fractional change in quizartinib CL for Black or African American race (unitless)")  # Table 3, RSE 38.0%; Eq 8

    e_sexf_vc <- -0.169
    label("Fractional change in quizartinib Vc for females vs the male reference (unitless)")  # Table 3, RSE 16.0%; Eq 11

    e_age_vp <- 0.0152
    label("Exponential coefficient for AGE centred at 47 years on quizartinib Vp1 (1/year)")  # Table 3, RSE 5.97%; Eq 13

    e_nonaml_ka <- -0.188
    label("Fractional change in quizartinib ka for non-AML subjects vs the AML reference (unitless)")  # Table 3, RSE 28.9%; Eq 10

    e_cyp3a4_inh_strong_fdepot <- 0.273
    label("Fractional change in quizartinib Frel with strong CYP3A inhibitors (unitless)")  # Table 3, RSE 10.7%; Eq 9

    e_cyp3a4_inh_mod_fdepot <- 0.116
    label("Fractional change in quizartinib Frel with moderate CYP3A inhibitors (unitless)")  # Table 3, RSE 13.8%; Eq 9

    e_nonaml_fdepot <- 1.73
    label("Quizartinib Frel multiplier for non-AML subjects vs the AML reference (unitless, a direct multiplier and not a fractional change)")  # Table 3 row "F rel for non-AML subjects" = 1.73, RSE 2.96%; Eq 9. Unlike the other categorical rows, Table 3 reports this as a LEVEL (1.73), not as a "fractional change"; the NONMEM stream confirms this (F1AML2 = THETA(12) directly, not 1 + THETA(12)), so model() applies it as (1 + (e_nonaml_fdepot - 1) * (1 - DIS_AML))

    e_trtph_induction_fdepot <- -0.419
    label("Fractional change in quizartinib Frel during QuANTUM-First induction (unitless)")  # Table 3, RSE 4.92%; Eq 9

    e_trtph_consolidation_fdepot <- -0.192
    label("Fractional change in quizartinib Frel during QuANTUM-First consolidation (unitless)")  # Table 3, RSE 15.6%; Eq 9

    e_trtph_continuation_fdepot <- 0.418
    label("Fractional change in quizartinib Frel during QuANTUM-First continuation (unitless)")  # Table 3, RSE 12.5%; Eq 9

    boxcox_fdepot <- -1.28
    label("Box-Cox shape parameter lambda for the Frel interindividual variability in AML patients (unitless)")  # Table 3 row "Shape box-cox parameter for IIV on F rel for AML patients" = -1.28, RSE 14.4%; quizartinib stream THETA(19) = -1.46038 (initial)

    # ------------------------------------------------------------------
    # AC886 (METABOLITE) STRUCTURAL PARAMETERS
    # ------------------------------------------------------------------
    # Vaddady 2024 Table 3, "Final AC886 model" column, reproduced in
    # Supporting Information Equations 16-20. These are the values BEFORE
    # the division by fMET that the source control stream applies (the
    # stream computes e.g. TVCLM = THETA(22) * (WT/75)^0.75 / FMET, and
    # Equation 17 prints only the THETA(22) part); model() applies the
    # /fMET step explicitly. The AC886 $THETA rows are initial estimates
    # (CLM initial 4.81815 vs Table 3 final 4.61), so Table 3 governs.

    lcl_ac886 <- log(4.61)
    label("AC886 apparent clearance CL/F at the AML reference before the fMET division (L/h)")  # Table 3: CL_AC886 = 4.61 L/h, RSE 3.75%; Eq 17

    lvc_ac886 <- log(8.93)
    label("AC886 apparent central volume Vc/F at the reference before the fMET division (L)")  # Table 3: Vc,AC886 = 8.93 L, RSE 6.84%; Eq 18

    lvp_ac886 <- log(68.5)
    label("AC886 apparent peripheral volume Vp/F before the fMET division (L)")   # Table 3: Vp,AC886 = 68.5 L, RSE 1.78%; Eq 20

    lq_ac886 <- log(3.76)
    label("AC886 apparent intercompartmental clearance Q/F before the fMET division (L/h)")  # Table 3: Q_AC886 = 3.76 L/h, RSE 2.06%; Eq 19

    fmet_base <- fixed(0.5)
    label("Parent-to-metabolite conversion fraction fMET at the reference (unitless)")  # Eq 16: fMET = 0.5 * (phase terms). Vaddady 2024 Methods: fMET "was fixed to 0.5, based on the quizartinib/AC886 exposure ratio in the Phase 2 study 2689-CL-2004 ... the model fit is not affected by the value of fMET"

    # ------------------------------------------------------------------
    # AC886 COVARIATE EFFECTS
    # ------------------------------------------------------------------

    e_nonaml_cl_ac886 <- 0.843
    label("Fractional change in AC886 CL for non-AML subjects vs the AML reference (unitless)")  # Table 3, RSE 11.0%; Eq 17

    e_race_black_cl_ac886 <- 0.488
    label("Fractional change in AC886 CL for Black or African American race (unitless)")  # Table 3, RSE 20.9%; Eq 17

    e_cyp3a4_inh_strong_cl_ac886 <- 0.298
    label("Fractional change in AC886 CL with strong CYP3A inhibitors (unitless)")  # Table 3, RSE 8.01%; Eq 17

    e_cyp3a4_inh_strong_vc_ac886 <- 2.79
    label("Fractional change in AC886 Vc with strong CYP3A inhibitors (unitless)")  # Table 3, RSE 7.34%; Eq 18

    e_trtph_induction_fmet <- 0.715
    label("Fractional change in fMET during QuANTUM-First induction (unitless)")   # Table 3, RSE 6.38%; Eq 16

    e_trtph_consolidation_fmet <- 0.272
    label("Fractional change in fMET during QuANTUM-First consolidation (unitless)")  # Table 3, RSE 13.3%; Eq 16

    e_trtph_continuation_fmet <- -0.249
    label("Fractional change in fMET during QuANTUM-First continuation (unitless)")  # Table 3, RSE 8.47%; Eq 16

    # ------------------------------------------------------------------
    # INTERINDIVIDUAL VARIABILITY
    # ------------------------------------------------------------------
    # Vaddady 2024 Table 3 reports each IIV in a column headed "(CV)".
    # Those numbers are the log-scale standard deviation omega, not a
    # variance and not a true lognormal CV: each squares to the
    # corresponding $OMEGA entry of the source control streams (e.g.
    # 0.695^2 = 0.4830 vs $OMEGA 2_IIV_CL_AML = 0.481313; 0.740^2 =
    # 0.5476 vs AC886 $OMEGA 2_CLM = 0.542698; 1.36^2 = 1.8496 vs
    # 3_V5 = 1.85491). nlmixr2 takes variances, so each entry below is
    # omega^2 with the reported omega shown in the expression.
    # No IIV was supported on Vp1, Vp2 or Q2 of quizartinib (the
    # quizartinib stream carries 4_IIV_V3AML as `0 FIX`), nor on the
    # AC886 peripheral volume or intercompartmental clearance.

    etalcl   ~ 0.695^2   # Table 3: IIV CL_quiz = 0.695, RSE 2.34%, shrinkage 11.8%
    etalq    ~ 0.691^2   # Table 3: IIV Q1,quiz = 0.691, RSE 4.60%
    etalvc   ~ 0.186^2   # Table 3: IIV Vc,quiz = 0.186, RSE 11.0%, shrinkage 60.2%
    etaltlag ~ 0.647^2   # Table 3: IIV Tlag = 0.647, RSE 4.10%, shrinkage 44.5%
    etalka   ~ 0.423^2   # Table 3: IIV ka = 0.423, RSE 6.77%, shrinkage 46.6%
    etald1   ~ 0.821^2   # Table 3: IIV D1 = 0.821, RSE 4.50%, shrinkage 38.1%

    # Stratum-suffixed IIV: Vaddady 2024 estimated two separate variances
    # for the same quantity, one per disease-status subpopulation, inside
    # a single joint fit ("Two different IIV terms for F rel were
    # estimated, one for AML patients and one for non-AML subjects, to
    # account for the higher IIV observed in patient data"). DIS_AML
    # selects which eta is active in model(). The AML-patient Frel eta is
    # additionally Box-Cox transformed.
    etalfdepot_aml     ~ 0.444^2   # Table 3: IIV F rel AML patients = 0.444, RSE 4.81%, shrinkage 17.5%
    etalfdepot_nonaml  ~ 0.256^2   # Table 3: IIV F rel non-AML subjects = 0.256, RSE 4.22%, shrinkage 6.58%

    etalcl_ac886_aml    ~ 0.740^2  # Table 3: IIV CL_AC886 AML patients = 0.740, RSE 2.25%, shrinkage 4.65%
    etalcl_ac886_nonaml ~ 0.516^2  # Table 3: IIV CL_AC886 non-AML subjects = 0.516, RSE 5.31%, shrinkage 1.56%

    etalvc_ac886 ~ 1.36^2          # Table 3: IIV Vc,AC886 = 1.36, RSE 3.55%, shrinkage 13.7%

    # ------------------------------------------------------------------
    # RESIDUAL UNEXPLAINED VARIABILITY
    # ------------------------------------------------------------------
    # Vaddady 2024 Equation 3: log(y) = log(yhat) + eps_add, i.e. additive
    # on the log scale, which is proportional on the linear concentration
    # scale in nlmixr2. Table 3 reports these on the same omega (SD) scale
    # as the IIV terms: 0.440^2 = 0.1936 vs the quizartinib stream's
    # $SIGMA 1_RUV_AML = 0.194347, and 0.452^2 = 0.2043 vs the AC886
    # stream's $SIGMA 1_RUV_AML = 0.225952.
    # Each moiety has TWO reported residual SDs -- one for AML patients
    # and one for non-AML subjects. nlmixr2 cannot express a
    # covariate-switched residual error, so the model file carries the
    # AML-patient value for each moiety (the population this analysis
    # targets, and 659 of 932 subjects); the non-AML values 0.102
    # (quizartinib) and 0.312 (AC886) are recorded here and in the
    # vignette Errata.

    propSd <- 0.440
    label("Quizartinib proportional residual SD on the linear concentration scale, AML patients (fraction)")  # Table 3: Additive RUV log scale AML patients = 0.440, RSE 0.363%. Non-AML subjects: 0.102 (not encodable; see vignette Errata)

    propSd_ac886 <- 0.452
    label("AC886 proportional residual SD on the linear concentration scale, AML patients (fraction)")  # Table 3: Additive RUV log scale AML patients = 0.452, RSE 0.364%. Non-AML subjects: 0.312 (not encodable; see vignette Errata)
  })

  model({
    # ------------------------------------------------------------------
    # Reference / conversion constants
    # ------------------------------------------------------------------
    # Allometric and age centring constants, Vaddady 2024 Equations 6 and 13.
    ref_wt  <- 75
    ref_age <- 47

    # Molecular weights, Vaddady 2024 Methods: quizartinib 560.68 g/mol
    # and AC886 576.67 g/mol. Metabolite formation conserves MOLES, so
    # when both amounts are carried in mg the formation flux is scaled by
    # the AC886/quizartinib molecular-weight ratio. (The source control
    # stream instead carries amounts in umol and converts at the
    # observation: CP = A(2)*560.68/V2 and CP = A(5)*576.67/V5.)
    mw_ratio_ac886 <- 576.67 / 560.68

    # ------------------------------------------------------------------
    # QUIZARTINIB individual parameters -- Equations 8 to 15
    # ------------------------------------------------------------------
    # Every categorical multiplier is Equation 5 (1 at the reference,
    # 1 + theta otherwise). The source column AML3 is the NON-AML flag, so
    # the canonical DIS_AML column enters as (1 - DIS_AML).

    cl <- exp(lcl + etalcl) * (WT / ref_wt)^e_wt_cl *
      (1 + e_cyp3a4_inh_strong_cl * CONMED_CYP3A4_INH_STRONG) *
      (1 + e_race_black_cl * RACE_BLACK)

    vc <- exp(lvc + etalvc) * (WT / ref_wt)^e_wt_vc *
      (1 + e_sexf_vc * SEXF)

    q <- exp(lq + etalq) * (WT / ref_wt)^e_wt_cl

    vp <- exp(lvp) * (WT / ref_wt)^e_wt_vc *
      exp(e_age_vp * (AGE - ref_age))

    q2 <- exp(lq2) * (WT / ref_wt)^e_wt_cl

    vp2 <- exp(lvp2) * (WT / ref_wt)^e_wt_vc

    ka <- exp(lka + etalka) * (1 + e_nonaml_ka * (1 - DIS_AML))

    tlag <- exp(ltlag + etaltlag)

    d1 <- exp(ld1 + etald1)

    # Relative bioavailability, Equation 9. The non-AML term is a direct
    # multiplier (1.73), so it enters as 1 + (multiplier - 1) * indicator
    # to keep Table 3's reported value verbatim. The three treatment-phase
    # indicators are mutually exclusive and are 0 for every non-newly-
    # diagnosed subject, so their product reduces to the single active
    # phase factor, matching the source stream's F1PHASE switch.
    frel_cov <- (1 + e_cyp3a4_inh_strong_fdepot * CONMED_CYP3A4_INH_STRONG) *
      (1 + e_cyp3a4_inh_mod_fdepot * CONMED_CYP3A4_INH_MOD) *
      (1 + (e_nonaml_fdepot - 1) * (1 - DIS_AML)) *
      (1 + e_trtph_induction_fdepot * TRTPH_INDUCTION) *
      (1 + e_trtph_consolidation_fdepot * TRTPH_CONSOLIDATION) *
      (1 + e_trtph_continuation_fdepot * TRTPH_CONTINUATION)

    # Frel interindividual variability. AML patients carry a Box-Cox
    # transformed eta (Vaddady 2024 Equation 2, Petersson 2009 form);
    # non-AML subjects carry a plain exponential eta. The source stream
    # writes the transform as (exp(eta)^lambda - 1)/lambda, which is
    # identically (exp(lambda * eta) - 1)/lambda used here (the form in
    # GonzalezSales_2015_testosterone.R).
    iiv_fdepot_aml <- (exp(boxcox_fdepot * etalfdepot_aml) - 1) / boxcox_fdepot
    iiv_fdepot     <- DIS_AML * iiv_fdepot_aml +
      (1 - DIS_AML) * etalfdepot_nonaml

    frel <- exp(lfdepot) * frel_cov * exp(iiv_fdepot)

    # ------------------------------------------------------------------
    # AC886 individual parameters -- Equations 16 to 20
    # ------------------------------------------------------------------
    # fMET, Equation 16. Baseline 0.5 with the same three mutually
    # exclusive treatment-phase indicators; the reference group here pools
    # R/R AML patients AND non-AML subjects (no non-AML effect on fMET was
    # demonstrated), which the all-zero indicator state already gives.
    fmet <- fmet_base *
      (1 + e_trtph_induction_fmet * TRTPH_INDUCTION) *
      (1 + e_trtph_consolidation_fmet * TRTPH_CONSOLIDATION) *
      (1 + e_trtph_continuation_fmet * TRTPH_CONTINUATION)

    # AC886 CL carries one of two stratum-specific etas, selected by DIS_AML.
    iiv_cl_ac886 <- DIS_AML * etalcl_ac886_aml +
      (1 - DIS_AML) * etalcl_ac886_nonaml

    # All four AC886 disposition parameters are divided by fMET, exactly as
    # the source control stream does (TVCLM = THETA(22)*(WT/75)^0.75/FMET,
    # and likewise for V5, V6 and QM). Paired with routing 100% of the
    # parent elimination flux into central_ac886 below (the stream sets
    # K25 = CL/V2 and K20 = 0, "consider that all the parent is
    # metabolized into the metabolite"), this is algebraically identical to
    # transferring only the fMET fraction with un-inflated parameters: the
    # fMET factors cancel in every rate constant, so fMET acts purely as a
    # multiplicative scale on the predicted AC886 concentration.
    cl_ac886 <- exp(lcl_ac886 + iiv_cl_ac886) * (WT / ref_wt)^e_wt_cl *
      (1 + e_nonaml_cl_ac886 * (1 - DIS_AML)) *
      (1 + e_race_black_cl_ac886 * RACE_BLACK) *
      (1 + e_cyp3a4_inh_strong_cl_ac886 * CONMED_CYP3A4_INH_STRONG) /
      fmet

    vc_ac886 <- exp(lvc_ac886 + etalvc_ac886) * (WT / ref_wt)^e_wt_vc *
      (1 + e_cyp3a4_inh_strong_vc_ac886 * CONMED_CYP3A4_INH_STRONG) /
      fmet

    vp_ac886 <- exp(lvp_ac886) * (WT / ref_wt)^e_wt_vc / fmet

    q_ac886 <- exp(lq_ac886) * (WT / ref_wt)^e_wt_cl / fmet

    # ------------------------------------------------------------------
    # Micro-constants
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    kel_ac886 <- cl_ac886 / vc_ac886
    k56       <- q_ac886  / vc_ac886
    k65       <- q_ac886  / vp_ac886

    # ------------------------------------------------------------------
    # ODE system (source: quizartinib ADVAN12 TRANS4; AC886 six-compartment
    # ADVAN13 $DES, Supporting Information NONMEM code)
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    d/dt(central_ac886)     <-  kel * central * mw_ratio_ac886 -
      kel_ac886 * central_ac886 -
      k56 * central_ac886 + k65 * peripheral1_ac886
    d/dt(peripheral1_ac886) <-  k56 * central_ac886 - k65 * peripheral1_ac886

    # ------------------------------------------------------------------
    # Absorption modifiers: sequential zero-order then first-order input
    # ------------------------------------------------------------------
    # The dose is delivered into the depot at a zero-order rate over D1
    # hours, beginning after the lag time Tlag, and then transfers to the
    # central compartment by first-order ka (Vaddady 2024 Results:
    # "sequential zero- and first-order absorption").
    # rxode2 only applies a modelled dur() to dose records carrying
    # rate = -2; a plain bolus silently ignores it and collapses Tmax onto
    # the lag time. Set rate = -2 (or supply dur= directly on the dose
    # record) in any event table built for this model.
    alag(depot) <- tlag
    dur(depot)  <- d1
    f(depot)    <- frel

    # ------------------------------------------------------------------
    # Observations
    # ------------------------------------------------------------------
    # Amounts are in mg and volumes in L, so amount/volume is mg/L =
    # ug/mL; the factor 1000 converts to the ng/mL scale Vaddady 2024
    # reports (LLOQ 2 and 0.5 ng/mL).
    Cc       <- 1000 * central       / vc
    Cc_ac886 <- 1000 * central_ac886 / vc_ac886

    Cc       ~ prop(propSd)
    Cc_ac886 ~ prop(propSd_ac886)
  })
}
