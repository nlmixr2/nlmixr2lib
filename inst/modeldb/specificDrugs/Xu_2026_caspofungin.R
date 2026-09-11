Xu_2026_caspofungin <- function() {
  description <- "Two-compartment population PK model with first-order elimination for intravenous caspofungin in critically ill Chinese children in the paediatric intensive care unit (PICU). Clearance and intercompartmental clearance scale with body surface area through a FIXED exponent of 0.66, and both volumes through a FIXED exponent of 1, standardised to a 0.79 m^2 individual; body surface area beat body weight, lean body weight and fat-free mass on OFV and AIC. Extracorporeal membrane oxygenation (ECMO) is the only clinical covariate retained and multiplies the central volume 18.2-fold, which lowers trough concentrations without changing steady-state AUC because clearance is unaffected. Inter-individual variability on clearance and central volume is strongly correlated (72%). Estimated by FOCE-I in NONMEM 7.5 with a combined proportional-plus-additive residual error. Xu 2026, n = 29 patients, 138 plasma samples, ages 0.33-16 years."
  reference <- "Xu N, Shi Y, Ju G, Liu X, Yan G, Zheng Y, Hou S, Xiang X, Lu G, Ouyang D, Zhu X, Wang Y. Population pharmacokinetics of caspofungin in critically ill Chinese children: a prospective observational study. Antimicrob Agents Chemother. 2026;70(2):e01277-25. doi:10.1128/aac.01277-25. PMC12888871. Received 22 August 2025, accepted 5 December 2025, published online 30 December 2025. ClinicalTrials.gov NCT04961593."
  vignette <- "Xu_2026_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Caspofungin was given as a once-daily 1 h intravenous
  # infusion and TOTAL plasma caspofungin was quantified by LC-MS/MS
  # (Methods, 'Caspofungin assay and fungal cultures'; calibration range
  # 0.05-50 ug/mL, LLOQ 0.05 ug/mL). The model therefore predicts total, not
  # unbound, concentration -- the authors list this explicitly as a
  # limitation, noting caspofungin's high protein binding.
  compartmentData <- list(
    central     = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area computed with the Mosteller formula from body weight and height",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The body-size descriptor of the final model, entering as a power",
        "term on all four disposition parameters: (BSA / 0.79)^0.66 on CL and",
        "Q, and (BSA / 0.79)^1 on V1 and V2. Both exponents were FIXED, not",
        "estimated (Table 2 rows 'BSA_CL 0.66 FIX' and 'BSA_V 1 FIX').",
        "Computed by Mosteller: BSA (m^2) = sqrt(weight (kg) * height (cm) /",
        "3600), stated explicitly in Methods, 'Study design'. The formula is",
        "reproducible from the paper's own Table 1 medians -- sqrt(16.0 * 104",
        "/ 3600) = 0.680 m^2 against a tabulated median BSA of 0.660 m^2, the",
        "residual gap being the usual median-of-ratios vs ratio-of-medians",
        "difference -- so this model may be driven with a correctly computed",
        "Mosteller BSA (contrast the Yao 2025 caution in the register entry).",
        "THE STANDARDISATION CONSTANT IS 0.79 m^2, NOT the cohort median of",
        "0.660 m^2. Results, 'Population pharmacokinetics analysis' states",
        "'scaled to a 0.79 m^2 individual', and three independent",
        "cross-checks in the Discussion confirm it against the Table 2",
        "thetas: BSA-normalised CL 0.196 / 0.79 = 0.248 L/h/m^2 (Discussion:",
        "0.248), BSA-normalised V1 2.22 / 0.79 = 2.81 L/m^2 (Discussion:",
        "2.81), and elimination rate constant 0.196 / 2.22 = 0.0883 /h",
        "(Discussion: 0.088). Dividing by the median 0.660 instead gives",
        "0.297 and 3.36, both of which contradict the Discussion. Note the",
        "Fig. 2 caption separately describes its reference subject as 'BSA",
        "0.66 m^2'; that is the covariate-effects forest plot's display",
        "reference, not the model's normalisation constant.",
        "Also the dosing metric: patients received 70 mg/m^2 loading and 50",
        "mg/m^2 maintenance, each capped at 70 mg. Time-fixed here --",
        "anthropometry was recorded at baseline, and the authors list the",
        "inability to capture time-varying covariates as a limitation. Must",
        "be strictly positive; it enters a power term. Studied range",
        "0.286-1.89 m^2 (Table 1)."
      ),
      source_name        = "BSA"
    ),
    ECMO_STATUS = list(
      description        = "Extracorporeal membrane oxygenation support indicator (1 = receiving ECMO)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ECMO support)",
      notes              = paste(
        "The ONLY clinical covariate retained by stepwise covariate",
        "modelling (forward P < 0.05, backward P < 0.01), and only on the",
        "central volume: dOFV = -13.262, P < 0.001 (Results, 'Population",
        "pharmacokinetics analysis'). Encoded as the power-of-indicator",
        "multiplier e_ecmo_status_vc^ECMO_STATUS so V1 is multiplied by 18.2",
        "when the patient is cannulated, following the Kang 2020 cefpirome",
        "and Watt 2015 fluconazole precedent in the register entry.",
        "The 18.2-FOLD reading (as opposed to a 1 + 18.2 = 19.2-fold",
        "increment) is the paper's own wording: the Discussion states 'ECMO",
        "was associated with a marked 18.2-fold increase in V1'. The Abstract",
        "calls 18.2 an 'effect coefficient' without a functional form, and",
        "the supplemental control stream reproduces only the earlier",
        "weight-based design model, which carries no ECMO term -- so the",
        "Discussion sentence is the only statement of the form. The",
        "distinction is 5% on V1 against a bootstrap 95% CI of 4.00-298.30,",
        "so no validation gate can separate the two readings; see the",
        "vignette Assumptions and deviations.",
        "Subject-level, not time-varying: the source does not report",
        "cannulation or decannulation times relative to caspofungin dosing,",
        "and single-day sampling after steady state gives no within-subject",
        "ECMO transition (contrast Kang 2020, where the indicator switches at",
        "decannulation). Only 4 of 29 patients (13.8%) received ECMO, which",
        "is why the effect carries 62% RSE and a bootstrap interval spanning",
        "two orders of magnitude; the authors state plainly that 'the number",
        "of patients receiving ECMO in this study was limited' and that the",
        "conclusions 'should be interpreted with caution'. An ECMO effect on",
        "clearance is absent from the final model, consistent with the",
        "paper's finding of no significant AUC difference between the ECMO",
        "and non-ECMO groups. Note 1 of the 4 ECMO patients was also on CRRT,",
        "so the two extracorporeal circuits are not fully separable in this",
        "cohort."
      ),
      source_name        = "ECMO"
    )
  )

  # Covariates the source SCREENED in stepwise covariate modelling but did
  # NOT retain in the final model. Documentation only -- none of these is
  # referenced in model(). The screened set is listed in Methods, 'Population
  # pharmacokinetic modeling of caspofungin'; the negative result for the
  # hepatic markers is stated in the Discussion ('Other factors, including
  # AST and ALT levels, did not appear to influence the PK of caspofungin').
  #
  # Two screened items are deliberately given no entry below. 'Severe
  # malnutrition' and 'renal transplant' (Table 1) have no canonical column,
  # and inventing one for a covariate that no model uses would be an
  # unwarranted register addition; they are described in population$notes
  # instead. Lean body weight and fat-free mass appear here because they were
  # tested as ALTERNATIVE body-size descriptors against BSA (Table S6), not
  # as covariates in the usual sense.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested a priori as an allometric body-size descriptor and rejected",
        "in favour of BSA. Table S6 run 2 (weight, exponents 0.75 / 1 FIXED)",
        "gives OFV 438.633 / AIC 456.633 against run 5 (BSA, 0.66 / 1 FIXED)",
        "at OFV 433.889 / AIC 451.889. Cohort median 16.0 kg (range 4.90-74.0;",
        "Table 1). Weight IS the body-size descriptor of the companion",
        "Xu_2026_caspofungin_optimalDesign.R model from the same paper's",
        "supplement, so a user who needs a weight-driven version of this model",
        "should load that one rather than substituting weight here."
      )
    ),
    LBW = list(
      description = "Lean body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested a priori as an allometric body-size descriptor and rejected in favour of BSA. Table S6 run 3 (exponents 0.75 / 1 FIXED): OFV 438.798 / AIC 456.798. Not tabulated in Table 1. The estimating formula is not stated in the source."
    ),
    FFM = list(
      description = "Fat-free mass",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested a priori as an allometric body-size descriptor and rejected in favour of BSA. Table S6 run 4 (exponents 0.75 / 1 FIXED): OFV 438.535 / AIC 456.535 -- the best of the three weight-like descriptors, still short of BSA. Not tabulated in Table 1. The estimating formula is not stated in the source."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in stepwise covariate modelling, not retained. Cohort median 5.33 years (range 0.330-16.0; Table 1); the Abstract and Results quote the range as 0.33-16 years and a median of 4.63 years, the latter matching the n = 14 intensive-sampling subset of Table S1 rather than the full cohort."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in stepwise covariate modelling, not retained. Table 1 reports 12 male / 17 female, i.e. 58.6% female. The source column records the male/female split, so the canonical SEXF requires the transformation SEXF = 1 - SEXM; immaterial here because the term is absent from the final model."
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened as a marker of immune status, not retained. Cohort median 4.92 x 10^9/L (range 0.100-54.1; Table 1)."
    ),
    RBC = list(
      description = "Red blood cell count",
      units       = "10^12/L",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 2.98 x 10^12/L (range 1.95-5.45; Table 1)."
    ),
    HGB = list(
      description = "Blood hemoglobin concentration",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 85.5 g/L (range 54.0-168; Table 1). Already in SI units, so no g/dL conversion is needed."
    ),
    AST = list(
      description = "Serum aspartate aminotransferase activity",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a hepatic-function marker, not retained in the final",
        "model -- the Discussion states 'Other factors, including AST and ALT",
        "levels, did not appear to influence the PK of caspofungin.' Cohort",
        "median 50.0 U/L (range 18.3-10,900; Table 1). Note that AST IS",
        "retained, on Q, in the companion",
        "Xu_2026_caspofungin_optimalDesign.R model fitted to the n = 14",
        "intensive-sampling subset; the effect did not survive into the final",
        "29-patient fit."
      )
    ),
    ALT = list(
      description = "Serum alanine aminotransferase activity",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker, not retained -- see the AST note for the Discussion sentence covering both. Cohort median 24.1 U/L (range 3.51-2,510; Table 1). Five children (reported as 35.7%, a percentage taken from the n = 14 subset) had abnormal liver function tests before caspofungin."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a hepatic-function marker, not retained. Cohort median 8.50 umol/L (range 2.10-329; Table 1). Direct bilirubin was also tabulated (median 4.10, range 1.00-262 umol/L) but is not listed among the screened covariates."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as the renal-function marker, not retained. Cohort median 29.7 umol/L (range 6.50-223; Table 1). Five patients (reported as 35.7%, again an n = 14 percentage) had moderate renal impairment with eGFR 30-59 mL/min/1.73 m^2."
    ),
    URIC_ACID = list(
      description = "Serum uric acid",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened, not retained. Cohort median 140 umol/L (range 6.69-736; Table 1)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal-replacement-therapy treatment-status indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened alongside ECMO_STATUS, not retained. Only 1 of 29 patients (3.4%) underwent CRRT (Table 1), so the indicator is essentially unestimable in this cohort; that single patient was also on ECMO, which is why the retained ECMO effect cannot be fully separated from CRRT."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 29L,
    n_studies      = 1L,
    n_samples      = 138L,
    age_range      = "0.330-16.0 years (Table 1). Eligibility 3 months to 18 years.",
    age_median     = "5.33 years (Table 1). The Abstract, Results and Table S1 quote 4.63 years, which is the median of the n = 14 intensive-sampling subset.",
    weight_range   = "4.90-74.0 kg (Table 1). The Abstract quotes a range of 4.9-64 kg, matching the n = 14 subset of Table S1.",
    weight_median  = "16.0 kg (Table 1). The Abstract quotes 15.9 kg, the n = 14 subset median.",
    height_median  = "104 cm (range 54.0-173; Table 1)",
    bsa_median     = "0.660 m^2 (range 0.286-1.89; Table 1). NOT the model's standardisation constant, which is 0.79 m^2.",
    sex_female_pct = 58.6,
    race_ethnicity = "Chinese. Single-centre enrolment at the Children's Hospital of Fudan University, Shanghai; the source reports no further race or ethnicity breakdown.",
    disease_state  = paste(
      "Critically ill children admitted to the paediatric intensive care unit",
      "and treated with caspofungin for suspected or proven invasive fungal",
      "infection. The cohort was deliberately enriched for the",
      "pathophysiology that perturbs caspofungin disposition: 4 of 29 on ECMO",
      "(13.8%), 1 on CRRT (3.4%), 2 with severe malnutrition (6.9%), 2 with",
      "hypoalbuminaemia (6.9%) and 1 renal-transplant recipient (3.4%)",
      "(Table 1). Five children had abnormal liver function tests and five",
      "had moderate renal impairment (eGFR 30-59 mL/min/1.73 m^2) before the",
      "first dose. Thirteen patients had a positive fungal culture during",
      "caspofungin infusion and ALL of them converted to culture-negative;",
      "the commonest isolate was Candida parapsilosis (n = 9), then C.",
      "albicans (n = 2), C. tropicalis (n = 1) and C. guilliermondii (n = 1)."
    ),
    ecmo_support   = paste(
      "4 of 29 patients (13.8%) received ECMO during caspofungin",
      "administration. The source does not report ECMO mode (VV vs VA), flow",
      "rate, or cannulation timing relative to dosing, so the covariate is",
      "encoded as a subject-level indicator. Bayesian post-hoc exposure was",
      "LOWER in the ECMO group -- median AUC(ss,24h) 117.08 h*mg/L (range",
      "17.20-176.04) vs 193.63 h*mg/L (range 53.18-380.95) in the 25",
      "non-ECMO patients -- but the difference did not reach significance",
      "(P = 0.06664), and the authors flag both the marginal P value and the",
      "severe group imbalance. The retained 18.2-fold effect on the central",
      "volume is the paper's mechanistic explanation (circuit sequestration,",
      "increased capillary permeability, altered protein binding); because",
      "clearance is unaffected, true steady-state AUC = dose / CL is",
      "identical between groups, which is exactly what the tAUC-based target",
      "attainment showed."
    ),
    renal_function = "Serum creatinine median 29.7 umol/L (range 6.50-223; Table 1). Five patients had moderate renal impairment (eGFR 30-59 mL/min/1.73 m^2). Neither creatinine nor CRRT was retained as a covariate.",
    hepatic_function = "ALT median 24.1 U/L (range 3.51-2,510), AST median 50.0 U/L (range 18.3-10,900), total bilirubin median 8.50 umol/L (range 2.10-329), albumin median 35.8 g/L (range 24.8-49.9) (Table 1). The extreme transaminase maxima reflect the PICU setting. Neither AST nor ALT was retained in the final model.",
    dose_range     = "Once-daily 1 h intravenous infusion on a BSA-based regimen: loading dose 70 mg/m^2 on day 1 and maintenance dose 50 mg/m^2 thereafter, each capped at 70 mg (per the product label, source reference 14). Monte Carlo dose simulations additionally explored maintenance doses of 10, 20, 30, 40, 50, 60 and 70 mg/m^2 for BSA <= 1.4 m^2, and flat daily doses above the 70 mg/day cap for BSA > 1.4 m^2.",
    regions        = "China (single centre: Children's Hospital of Fudan University, National Children's Medical Center, Shanghai), 1 November 2022 to 30 December 2024",
    sampling       = "Two stages. Stage 1 (intensive, n = 14) sampled pre-dose and 1, 2, 4, 8 (if feasible) and 16 h (if feasible) after the sixth dose. Stage 2 used an optimal sparse design derived from the stage-1 model -- windows of 0-1, 0.5-1.5, 6-7 and 23-24 h after the sixth dose. 138 concentrations total, median 7.475 mg/L (range 0.155-58.300), ALL above the 0.05 ug/mL LLOQ, so the planned M1 below-quantification-limit handling never had to be applied.",
    protein_binding = "Not fitted. Caspofungin is highly protein-bound; the authors assayed and modelled TOTAL plasma concentration, list unbound concentrations as future work, and note that albumin has been reported to influence caspofungin PK elsewhere. The tAUC(ss,24h)/MIC targets the paper simulates against are likewise defined on total drug.",
    notes          = paste(
      "Structural model selection: a two-compartment model was selected,",
      "consistent with three of the four prior paediatric caspofungin popPK",
      "studies. Body-size incorporation was tested a priori across weight,",
      "lean body weight, fat-free mass and BSA, with both FIXED-exponent and",
      "estimated-exponent variants (Table S6, 9 runs); BSA with FIXED",
      "exponents 0.66 on CL/Q and 1 on V1/V2 won on OFV and AIC. Estimated by",
      "FOCE-I in NONMEM 7.5. Evaluated by convergence, OFV, AIC, parameter",
      "precision, goodness-of-fit plots, a prediction-corrected VPC (Fig. 1F,",
      "Fig. S2G) and a 1,000-sample bootstrap whose medians sit close to the",
      "point estimates with the estimates inside the 95% CIs (Table 2).",
      "Marginal covariate effects on AUC(ss,24h) and C(min,ss) were explored",
      "with the coveffectsplot R package (Fig. 2). The dosing-simulation",
      "layer (1,000-subject Monte Carlo probability of target attainment",
      "against C. albicans, C. glabrata and C. parapsilosis under both",
      "tAUC(ss,24h)/MIC and C(min,ss) targets, Figs 3-5) is an application of",
      "this model rather than a separate model. The authors' dosing",
      "conclusion is that the licensed BSA-based maintenance dose is",
      "appropriate for BSA <= 1.4 m^2 while a flat daily dose is preferable",
      "above 1.4 m^2.",
      "A SECOND, independently fitted model from this paper is extracted",
      "separately as Xu_2026_caspofungin_optimalDesign.R: the supplemental",
      "Table S2 model, fitted to the n = 14 intensive-sampling subset and used",
      "to drive the $DESIGN sparse-sampling optimisation for stage 2. It uses",
      "weight rather than BSA allometry and retains an AST effect on Q that",
      "did not survive into this final fit.",
      "Limitations the authors state: no post-first-dose sampling, so the",
      "loading dose cannot be evaluated; single-day sampling cannot capture",
      "time-varying covariates or inter-occasion variability; total rather",
      "than unbound concentrations; and only 4 ECMO patients."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural fixed effects -- Xu 2026 Table 2, 'Estimates' column. The
    # paired 'Bootstrap (n = 1,000) [95% CI]' column is quoted in the
    # comments below for confidence in the point estimates, but is NOT
    # carried into any omega: it is parameter precision, not between-subject
    # spread.
    #
    # Reference subject: BSA = 0.79 m^2, not on ECMO. See
    # covariateData$BSA$notes for the three Discussion cross-checks that
    # confirm 0.79 m^2 rather than the cohort median of 0.660 m^2.
    # -----------------------------------------------------------------------

    lcl <- log(0.196)
    label("Clearance at BSA = 0.79 m^2 (L/h)")
    # Table 2: CL = 0.196 L/h (RSE 17%); bootstrap median 0.194
    # (95% CI 0.157-0.247). Cross-check, Discussion: BSA-normalised CL
    # 0.196 / 0.79 = 0.248 L/h/m^2, matching the quoted 0.248 L/h/m^2.

    lvc <- log(2.22)
    label("Central volume of distribution V1 at BSA = 0.79 m^2, without ECMO (L)")
    # Table 2: V1 = 2.22 L (RSE 41%); bootstrap median 2.10
    # (95% CI 0.58-3.38). Cross-check, Discussion: BSA-normalised V1
    # 2.22 / 0.79 = 2.81 L/m^2, matching the quoted 2.81 L/m^2. A third
    # cross-check ties CL and V1 together: kel = 0.196 / 2.22 = 0.0883 /h
    # against the Discussion's quoted 0.088 /h.

    lq <- log(1.01)
    label("Intercompartmental clearance Q at BSA = 0.79 m^2 (L/h)")
    # Table 2: Q = 1.01 L/h (RSE 57%); bootstrap median 0.95
    # (95% CI 0.25-3.47). The Table 2 footnote glosses Q as
    # 'inter-compartment exchange rate'; it is a flow in L/h, and the
    # supplemental control stream confirms it (ADVAN3 TRANS4 takes CL, V1,
    # Q, V2).

    lvp <- log(1.63)
    label("Peripheral volume of distribution V2 at BSA = 0.79 m^2 (L)")
    # Table 2: V2 = 1.63 L (RSE 30%); bootstrap median 1.89
    # (95% CI 1.20-4.29).

    # -----------------------------------------------------------------------
    # Body-surface-area allometry. Both exponents were FIXED, hence fixed().
    # Table 2 rows 'BSA_CL 0.66 FIX' and 'BSA_V 1 FIX' carry no RSE and no
    # bootstrap interval ('/' in both columns).
    #
    # ONE exponent is SHARED between CL and Q, and one between V1 and V2 --
    # hence the shared-exponent names e_bsa_cl_q and e_bsa_vc_vp rather than
    # four separate parameters. The evidence is the Table S6 column headers,
    # which are literally 'power_CL/Q' and 'power_V1/V2', and the
    # supplemental control stream for the companion design model, which
    # applies its 0.75 to both CL and Q and its 1 to both V1 and V2. Table 2
    # abbreviates the same pair as 'BSA_CL' and 'BSA_V', and the Abstract as
    # 'exponential 1 for volume of distribution and 0.66 for clearance'.
    # -----------------------------------------------------------------------

    e_bsa_cl_q <- fixed(0.66)
    label("Allometric exponent on (BSA / 0.79) for CL and Q (unitless)")
    # Table 2: BSA_CL = 0.66 FIX. Discussion: 'our current study identified a
    # fixed BSA-exponent of 0.66, consistent with allometric theory, and
    # demonstrated stability across a wider PICU cohort (0.33-16 years, BSA
    # 0.29-1.89 m2)'. Table S6 run 5 is the winning run.

    e_bsa_vc_vp <- fixed(1)
    label("Allometric exponent on (BSA / 0.79) for V1 and V2 (unitless)")
    # Table 2: BSA_V = 1 FIX. Abstract: 'exponential 1 for volume of
    # distribution'. Table S6 run 5 column 'power_V1/V2'.

    # -----------------------------------------------------------------------
    # ECMO effect on the central volume, applied in model() as the
    # power-of-indicator multiplier e_ecmo_status_vc^ECMO_STATUS so that V1
    # is multiplied by 18.2 for a cannulated patient and left unchanged
    # otherwise. This follows the Kang_2020_cefpirome.R and
    # Watt_2015_fluconazole.R encoding in the ECMO_STATUS register entry.
    # -----------------------------------------------------------------------

    e_ecmo_status_vc <- 18.2
    label("Multiplicative factor on V1 for ECMO_STATUS = 1 (vs 0) (unitless)")
    # Table 2: ECMO_V1 = 18.2 (RSE 62%); bootstrap median 17.8
    # (95% CI 4.00-298.30). Estimated, not fixed. The multiplicative
    # ('18.2-fold') rather than incremental ('1 + 18.2') reading comes from
    # the Discussion: 'ECMO was associated with a marked 18.2-fold increase
    # in V 1'. Selected by SCM with dOFV = -13.262, P < 0.001. The very wide
    # bootstrap interval reflects n = 4 ECMO patients.

    # -----------------------------------------------------------------------
    # Inter-individual variability: exponential on CL and V1, correlated.
    # Methods: 'The inter-individual variability (IIV) was described by an
    # exponential model'. Results: 'There was a high correlation between the
    # IIV in CL and the V 1 (Corr = 0.802); therefore, a correlation
    # coefficient was introduced' -- 0.802 is the diagnosed correlation that
    # motivated the block; the ESTIMATED correlation in the final model is
    # the Table 2 row 'Cor.CL.V 1' = 72%.
    #
    # OMEGA SCALE. Table 2 heads the block 'Inter-individual variability
    # (%CV)', which is ambiguous between omega-as-SD and the log-normal
    # back-transform omega^2 = log(1 + CV^2). This paper settles it
    # DEFINITIVELY, because the supplement publishes both the %CV table and
    # the raw NONMEM $OMEGA for the companion design model:
    #
    #   Table S2 IIV CL  64.4%   $OMEGA BLOCK(2) 0.415        sqrt(0.415) = 0.6442
    #   Table S2 IIV V1  50.3%                   0.253        sqrt(0.253) = 0.5030
    #   Table S2 Cor     82.4%                   0.267 (cov)  0.267 / sqrt(0.415*0.253) = 0.8240
    #
    # All three reproduce the printed percentages EXACTLY, so these authors'
    # '%CV' is the raw omega SD times 100. The log-normal back-transform
    # would give sqrt(exp(0.415) - 1) = 71.7%, not 64.4%. Applying the same
    # convention to Table 2:
    #
    #   etalcl variance = 0.507^2                = 0.257049
    #   etalvc variance = 0.91^2                 = 0.828100
    #   covariance      = 0.72 * 0.507 * 0.91    = 0.332186
    #
    # Shrinkage was low for both etas (4% on CL, 7% on V1), so the
    # stochastic layer is well supported despite n = 29.
    # -----------------------------------------------------------------------

    etalcl + etalvc ~ c(0.257049,
                        0.332186, 0.828100)
    # Table 2: IIV CL 50.7% (RSE 24%, SHR 4%), bootstrap 47.4%
    # [24.3%-72.7%]; IIV V1 91% (RSE 24%, SHR 7%), bootstrap 92.7%
    # [46.0%-174.1%]; row 'Cor.CL.V 1' 72%, bootstrap 54.7% [24.1%-83.3%].
    # See the omega-scale derivation above.

    # -----------------------------------------------------------------------
    # Residual error: combined proportional plus additive on total
    # concentration. Methods: 'the residual error for caspofungin
    # concentration was tested by proportional error, additive error, and
    # combined error model'; Table 2 reports both components, so the
    # combined form was selected.
    #
    # nlmixr2's prop() + add() forms the residual SD as
    # sqrt((propSd * f)^2 + addSd^2), which is algebraically identical to the
    # NONMEM structure Y = IPRED + IPRED*EPS(1) + EPS(2) that the
    # supplemental control stream uses for the companion model.
    #
    # Both components are read as SDs, not variances. The proportional row
    # is decisive: it is reported as a PERCENT (17.9%), which only makes
    # sense on the SD scale, and the additive row sits in the same
    # 'Estimates' column under the same heading with units of mg/L.
    # -----------------------------------------------------------------------

    propSd <- 0.179
    label("Proportional residual error (fraction)")
    # Table 2: 'Prop.error for total concentration, %' = 17.9%
    # (RSE 13%, SHR 17%); bootstrap 17.9% [11.8%-27.6%].

    addSd <- 0.838
    label("Additive residual error (mg/L)")
    # Table 2: 'Add.error for total concentration, mg/L' = 0.838
    # (RSE 29%, SHR 17%); bootstrap 0.811 [0.180-1.222]. Comfortably above
    # the 0.05 ug/mL assay LLOQ and about 11% of the median observed
    # concentration of 7.475 mg/L.
  })

  model({
    # Body-surface-area allometry, standardised to a 0.79 m^2 individual
    # (Results: 'scaled to a 0.79 m 2 individual'). One shared exponent for
    # the two clearances and one for the two volumes -- see the ini() block
    # and covariateData$BSA$notes.
    bsa_cl_factor <- (BSA / 0.79)^e_bsa_cl_q
    bsa_v_factor  <- (BSA / 0.79)^e_bsa_vc_vp

    # Individual disposition parameters. Exponential IIV on CL and V1 only;
    # Q and V2 carry no random effect in the final model (Table 2 lists IIV
    # rows for CL and V1 alone). ECMO multiplies V1 by e_ecmo_status_vc when
    # ECMO_STATUS = 1 and by 1 when it is 0.
    cl <- exp(lcl + etalcl) * bsa_cl_factor
    vc <- exp(lvc + etalvc) * bsa_v_factor * e_ecmo_status_vc^ECMO_STATUS
    q  <- exp(lq)           * bsa_cl_factor
    vp <- exp(lvp)          * bsa_v_factor

    # Linear two-compartment disposition with intravenous input into the
    # central compartment. Caspofungin was given as a 1 h infusion
    # (Methods, 'Study design'); the infusion duration is encoded on the
    # dose record by the user via rate or dur.
    #
    # The main text describes the final model as 'a two-compartment model
    # with first-order absorption and first-order elimination', but there is
    # NO absorption parameter: Table 2 lists only CL, V1, Q and V2, the drug
    # is available only as an intravenous infusion, and the supplemental
    # control stream specifies ADVAN3 TRANS4 -- two-compartment IV, no depot.
    # The 'first-order absorption' clause is a text error and no depot
    # compartment is implemented. See the vignette Assumptions and
    # deviations.
    d/dt(central)     <- -(cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-        q  / vc * central - q / vp * peripheral1

    # Total plasma caspofungin concentration. Dose mg, vc L -> Cc mg/L
    # (= ug/mL, the unit of the assay calibration range).
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
