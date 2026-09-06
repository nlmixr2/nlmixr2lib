Wu_2024_daptomycin <- function() {
  description <- paste(
    "Two-compartment population PK model with linear elimination for",
    "intravenous daptomycin in critically ill adult Han Chinese patients in a",
    "single Wuhan ICU (64 patients, 737 serum concentrations, 500 mg q24h as a",
    "30-min infusion, April 2021 to December 2022). Clearance is piecewise on",
    "continuous renal replacement therapy status: subjects on CRRT carry a",
    "single fixed total clearance of 0.386 L/h with no creatinine-clearance",
    "term, while subjects not on CRRT carry an additive non-renal plus renal",
    "decomposition CL = 0.229 + 0.148 * (CCR / 54) L/h, where CCR is the raw",
    "Cockcroft-Gault creatinine clearance in mL/min and 54 mL/min is the",
    "cohort median. Creatinine clearance was the only covariate retained:",
    "body weight, body mass index, age, sex, serum albumin, SOFA score and",
    "APACHE II score were all screened and rejected, and a sex effect on the",
    "peripheral volume (males about 1.4-fold higher) survived forward",
    "inclusion but was dropped in backward elimination. Inter-individual",
    "variability is log-normal on total clearance, central volume and",
    "peripheral volume, with no off-diagonal elements retained; residual",
    "variability is combined additive plus proportional. The paper's Monte",
    "Carlo simulations target AUC24h/MIC >= 666 at MIC 1 mg/L and conclude",
    "that 500 mg q24h suffices on CRRT and in renal impairment, while",
    "patients with CCR >= 90 mL/min need 700 mg daily to reach 90% probability",
    "of target attainment."
  )
  reference <- paste(
    "Wu J, Zheng X, Zhang L, Wang J, Lv Y, Xi Y, Wu D. Population",
    "pharmacokinetics of intravenous daptomycin in critically ill patients:",
    "implications for selection of dosage regimens. Front Pharmacol.",
    "2024;15:1378872. doi:10.3389/fphar.2024.1378872."
  )
  vignette <- "Wu_2024_daptomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE -- the assay matrix is stated in
  # Wu 2024 Methods ("Serum daptomycin concentrations were measured using a
  # high-performance liquid chromatography-tandem mass spectrometry method",
  # LLOQ 0.05 ug/mL, daptomycin-d5 internal standard).
  compartmentData <- list(
    central     = list(analyte = "daptomycin", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "daptomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CCR. Wu 2024 Methods: 'The patient's creatinine",
        "clearance (CLCR) was calculated according to the Cockcroft-Gault",
        "formula based on serum creatinine at steady state on day 3' -- so the",
        "column is a single time-fixed day-3 value per subject, not a",
        "time-varying series. Raw mL/min, NOT BSA-normalized to",
        "mL/min/1.73 m^2. Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md, which accepts raw mL/min when",
        "the source paper does not apply BSA normalization, with the assay",
        "form recorded here (precedent: Delattre 2010 amikacin,",
        "Shekar 2014 meropenem). Enters the model only for subjects with",
        "RRT_CRRT_STATUS = 0, as the dimensionless ratio (CRCL / 54) scaling",
        "the renal clearance arm; 54 mL/min is the cohort median stated in the",
        "text following the final covariate equation ('where CCR/54 is the",
        "corresponding median standardized individual CCR calculated for the",
        "current patient population') and matches Table 1 (median 54.25",
        "mL/min, range 8.3-200.2). Because 54 is the median of the WHOLE",
        "64-patient cohort including the 39 CRRT subjects, it is not the",
        "median of the 25 non-CRRT subjects the renal arm actually applies to;",
        "the paper does not report the non-CRRT subgroup median. Cockcroft-",
        "Gault CrCL is conventionally ill-defined for CRRT-dependent subjects,",
        "and the model switches the CRCL term off entirely when",
        "RRT_CRRT_STATUS = 1.",
        collapse = " "
      ),
      source_name        = "CCR"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Subject-level binary indicator for continuous renal replacement therapy during the PK sampling period",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Source concept CRRT. 1 = subject was receiving continuous renal",
        "replacement therapy; 0 = not on CRRT. Table 1: 39/64 (60.9%) on CRRT,",
        "25/64 (39.1%) not on CRRT. Filters used were MultiFiltrate Kit Ci-Ca",
        "AV1000S (n = 17, membrane area 1.8 m^2), Ultraflux AV600S (n = 2,",
        "1.4 m^2) and Prismaflex M100AN69 (n = 20, 0.9 m^2); blood flow",
        "95-350 mL/min and effluent flow 1600-2000 mL/h. The published model",
        "treats all three filters identically as a single binary covariate and",
        "does not use effluent flow, so RRT_CRRT_EFFLUENT_FLOW is NOT a column",
        "of this model. Treated as time-fixed at the subject level: the source",
        "does not resolve per-session CRRT start/stop times, hence",
        "RRT_CRRT_STATUS rather than the time-varying RRT_CRRT_ACTIVE.",
        "The indicator is a full switch, not a multiplier -- when it is 1 the",
        "entire non-renal-plus-renal clearance expression is replaced by the",
        "single estimate 0.386 L/h, so a CRRT subject carries no separate",
        "non-renal clearance term.",
        collapse = " "
      ),
      source_name        = "CRRT"
    )
  )

  # Covariates Wu 2024 screened in the stepwise covariate search but did NOT
  # retain in the final model. Documentation only -- none is referenced in
  # model(). Supplementary Table S1 ("Model building process") gives the OFV
  # for each candidate.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Tested on CL and on volume of distribution; not retained. Table 1",
        "median 64.5 kg (range 45-170; 63 of 64 patients were 45-90 kg and one",
        "was 170 kg). Wu 2024 Discussion attributes the null weight effect on",
        "CL to weight already being embedded in the Cockcroft-Gault CCR.",
        collapse = " "
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Tested as a body-size descriptor on CL; not retained. Table 1 median 23.0 kg/m^2 (range 16.5-52.5)."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on CL; not retained. Supplementary Table S1 model 5 (base + AGE on CL): dOFV = -1.631, p > 0.05. Table 1 mean 57.5 +/- 16.5 years."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested on the volume of distribution. Supplementary Table S1 model 4",
        "(base + GEND on VP) gave dOFV = -5.098 (p < 0.05) so it entered during",
        "forward inclusion, and the Results note the peripheral volume of males",
        "was approximately 1.4-fold that of females, but backward elimination",
        "removed it (Table S1 model 8, model6 - GEND: dOFV = 5.253, p > 0.001).",
        "No point estimate for the effect is reported in Table 2, so it cannot",
        "be encoded even as a documented alternative. Table 1: 43 male (67.2%),",
        "21 female (32.8%).",
        collapse = " "
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Tested on CL; not retained (Results: 'no relationships were identified between the CL, VC, or VP of daptomycin and APACHE II score, age, or serum albumin concentration'). Table 1 median 32.8 g/L (range 19.6-46.8)."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units       = "(score)",
      type        = "continuous",
      notes       = "Tested on CL; not retained. Table 1 median 24 (range 12-47)."
    ),
    ECMO_STATUS = list(
      description = "Extracorporeal membrane oxygenation treatment-status indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Wu 2024 Discussion: 'no effect of the combined use of extracorporeal",
        "membrane oxygenation on the pharmacokinetics in critically ill adult",
        "patients was observed'. Table 1: 6/64 (9.4%) on ECMO, of whom 4 were",
        "also on CRRT and 2 were on ECMO alone. The ECMO subgroup is too small",
        "to separate from CRRT and no point estimate is reported.",
        collapse = " "
      )
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units       = "(score)",
      type        = "continuous",
      notes       = paste(
        "Tested on CL and rejected: Results state the SOFA score 'was",
        "marginally associated with individual-specific CL and was excluded",
        "from the backward elimination process', and Supplementary Table S1",
        "model 3 (base + SOFA on CL) gives dOFV = -1.224, p > 0.05. Table 1",
        "median 12 (range 4-19). NOTE: as of this extraction 'SOFA' is NOT a",
        "registered canonical in inst/references/covariate-columns.md (the only",
        "mention of the score there is incidental, inside the LACT entry's",
        "notes). It is recorded here under its plain source name purely to",
        "preserve the provenance of the paper's covariate screen; because the",
        "covariate was rejected and is never referenced in model(), no new",
        "canonical was registered for it. A future paper that RETAINS a SOFA",
        "effect should ratify a canonical rather than inherit this key.",
        collapse = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 64L,
    n_studies      = 1L,
    n_observations = 737L,
    age_range      = "Not reported as a range; mean 57.5 +/- 16.5 years (Table 1). Inclusion required age >= 18 years.",
    age_median     = "Not reported (Table 1 gives mean +/- SD, not median)",
    weight_range   = "45-170 kg (Table 1 range; 63 of 64 patients weighed 45-90 kg and one patient was extremely obese at 170 kg)",
    weight_median  = "64.5 kg",
    sex_female_pct = 32.8,
    race_ethnicity = "Han Chinese (single-centre Chinese ICU cohort; the paper describes the population as 'critically ill adult Han Chinese patients')",
    disease_state  = paste(
      "Critically ill adults in a medical/surgical ICU treated with",
      "intravenous daptomycin for proven or probable Gram-positive infection",
      "unresponsive to or intolerant of other antimicrobials. Severity at",
      "baseline: APACHE II median 24 (range 12-47), SOFA median 12 (range",
      "4-19). Exclusions were pregnancy, known or suspected hypersensitivity",
      "to daptomycin or its excipients, and incomplete or missing laboratory",
      "data before and after treatment.",
      collapse = " "
    ),
    dose_range     = "Daptomycin 500 mg every 24 h, administered as a 30-min intravenous infusion in 100 mL normal saline (Jiangsu Hengrui Medicine Co., Ltd.). All subjects received the same regimen; the 400/500/600/700 mg daily doses appearing in Table 3 are simulated, not administered.",
    regions        = "China (single centre: intensive care unit of Zhongnan Hospital of Wuhan University, Wuhan; ethics approval LYL2019018)",
    renal_function = "Creatinine clearance (Cockcroft-Gault, day-3 serum creatinine) median 54.25 mL/min (range 8.3-200.2); serum creatinine median 106 umol/L (range 21.7-845.4). 39/64 (60.9%) were on CRRT.",
    extracorporeal = "CRRT 39/64 (60.9%); ECMO 6/64 (9.4%), of which 2 ECMO alone and 4 ECMO with CRRT. CRRT filters: MultiFiltrate Kit Ci-Ca AV1000S (n = 17), Ultraflux AV600S (n = 2), Prismaflex M100AN69 (n = 20); blood flow 95-350 mL/min, effluent flow 1600-2000 mL/h.",
    hepatic_function = "Albumin median 32.8 g/L (19.6-46.8); total bilirubin median 35.3 umol/L (5.2-312); ALT median 63 U/L (8-2581); AST median 71 U/L (11-3845).",
    notes          = paste(
      "Baseline demographics per Wu 2024 Table 1; the footnote states",
      "categorical variables are n (%) and the remaining variables are median",
      "(range), except age, which is printed as mean +/- SD. Sampling: day 1 at",
      "0, 0.5, 1, 2, 4, 8, 12 and 24 h after infusion; day 3 pre-dose and at 0,",
      "0.5, 1, 2, 4, 8, 12 and 24 h; days 5 and 7 pre-dose and at 0.5 and 4 h.",
      "Assay: HPLC-MS/MS, daptomycin-d5 internal standard, LLOQ 0.05 ug/mL.",
      "No daptomycin-related CPK elevation or other adverse reaction was",
      "observed. Two internal inconsistencies in the source are worth flagging:",
      "the Results text says 'Among these 45 CRRT patients' while Table 1 and",
      "the preceding sentence both give 39 CRRT patients, and Table 1 lists",
      "total protein as '54.2(37.1 - 573)' where the upper bound is almost",
      "certainly a misprint for 57.3 g/L. Neither affects the model.",
      collapse = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural fixed effects. Wu 2024 prints the final covariate model as a
    # brace-delimited piecewise equation on page 4 (immediately after "The
    # final covariate model is as follows:"):
    #
    #   CL (L/h) = { 0.386                     in subjects undergoing CRRT
    #              { 0.229 + 0.148 * (CCR/54)  in others
    #   VC (L) = 4.14 ; VP (L) = 3.52 ; Q (L/h) = 2.09
    #
    # Every value below is taken from that equation. Four of the five also
    # appear identically in Table 2 (NR = 0.229, VC = 4.14, VP = 3.52,
    # Q = 2.09, CRRT = 0.386). The ONE disagreement is the renal-clearance
    # coefficient: the equation prints 0.148, Table 2 row "R" prints 0.152,
    # and the Discussion prints 0.14. See the lcl_renal comment below.
    # ------------------------------------------------------------------

    # Non-renal clearance arm, non-CRRT subjects only.
    lcl_nonren <- log(0.229)
    label("Non-renal clearance, subjects not on CRRT (L/h)")
    # Wu 2024 final covariate model equation, "others" branch, intercept;
    # identical to Table 2 row "NR (L/h)" = 0.229 (RSE 25.1%, bootstrap
    # median 0.231, 95% CI 0.125-0.36).

    # Renal clearance arm at the reference CCR of 54 mL/min, non-CRRT only.
    lcl_renal <- log(0.148)
    label("Renal clearance at CCR = 54 mL/min, subjects not on CRRT (L/h)")
    # Wu 2024 final covariate model equation, "others" branch, slope on
    # (CCR/54). THREE-WAY SOURCE CONFLICT, resolved in favour of the printed
    # equation: the equation gives 0.148, Table 2 row "R" gives 0.152 (RSE
    # 31.3%, bootstrap median 0.150, 95% CI 0.062-0.264), and the Discussion
    # gives "(0.14 +/- 0.035) x (CCR/54)". The Discussion set is a different
    # run entirely -- it also disagrees with Table 2 on VC (4.20 vs 4.14), VP
    # (3.67 vs 3.52), Q (2.13 vs 2.09) and the CRRT clearance (0.388 vs 0.386)
    # -- so it carries no weight for the final model. That leaves equation vs
    # Table 2, which no on-disk source resolves: the Frontiers supplement is a
    # model-building OFV table with no parameter values and no control stream,
    # and the Table 3 Monte Carlo AUCs are too noisy to arbitrate (see the
    # vignette Errata). The equation is preferred per the standing
    # printed-equation-over-other-text rule. Impact is negligible: at the
    # reference CCR = 54 the renal arm is 39% of total CL, so 0.148 vs 0.152
    # moves total non-CRRT CL by 1.1%, far inside the parameter's own 31.3%
    # RSE.

    # Total clearance for subjects on CRRT. Replaces (does not add to) the
    # non-renal + renal expression above.
    lcl_crrt <- log(0.386)
    label("Total clearance, subjects on CRRT (L/h)")
    # Wu 2024 final covariate model equation, "CRRT" branch; identical to
    # Table 2 row "CRRT (L/h)" = 0.386 (RSE 5.1%, bootstrap median 0.386,
    # 95% CI 0.351-0.428).

    lvc <- log(4.14)
    label("Central volume Vc (L)")
    # Wu 2024 final covariate model equation "VC (L) 4.14"; Table 2 row
    # "V C (L)" = 4.14 (RSE 4.8%, bootstrap median 4.14, 95% CI 3.763-4.565).

    lvp <- log(3.52)
    label("Peripheral volume Vp (L)")
    # Wu 2024 final covariate model equation "VP (L) 3.52"; Table 2 row
    # "V P (L)" = 3.52 (RSE 8.2%, bootstrap median 3.527, 95% CI 2.987-4.078).

    lq <- log(2.09)
    label("Intercompartmental clearance Q (L/h)")
    # Wu 2024 final covariate model equation "Q (L/h) 2.09"; Table 2 row
    # "Q (L/h)" = 2.09 (RSE 8.2%, bootstrap median 2.086, 95% CI 1.736-2.482).

    # ------------------------------------------------------------------
    # Inter-individual variability. Wu 2024 Methods: "The interindividual
    # variability (eta) was described using lognormal distributions".
    # NONMEM $OMEGA diagonals are VARIANCES, and the Table 2 "Omega" block is
    # read as such: 0.091 -> sqrt(exp(0.091) - 1) = 30.9% CV on CL, which is a
    # plausible ICU value, whereas reading 0.091 as an SD would imply a 9.1%
    # CV that is implausibly tight for this population and inconsistent with
    # the reported shrinkage. Off-diagonals were investigated (Methods) but
    # none is reported in Table 2, so the matrix is diagonal.
    # ------------------------------------------------------------------
    etalcl ~ 0.091 # Table 2 Omega eta[CL]  (RSE 9.1%,  shrinkage 6.8%,  bootstrap median 0.085, 95% CI 0.059-0.122); 30.9% CV
    etalvc ~ 0.114 # Table 2 Omega eta[V C] (RSE 9.6%,  shrinkage 11.4%, bootstrap median 0.110, 95% CI 0.074-0.157); 34.7% CV
    etalvp ~ 0.202 # Table 2 Omega eta[V P] (RSE 14.8%, shrinkage 29.6%, bootstrap median 0.193, 95% CI 0.088-0.350); 47.3% CV

    # ------------------------------------------------------------------
    # Residual error. Wu 2024 Methods: "a mixed additive and proportional
    # model was chosen for residual variability". NONMEM $SIGMA entries are
    # VARIANCES and nlmixr2 propSd / addSd are SDs, so both are square-rooted.
    # The additive term is decisive for the variance reading: 38.095 as an SD
    # would be 38 mg/L of additive noise, larger than most observed troughs
    # in this cohort, whereas sqrt(38.095) = 6.17 mg/L is credible against
    # peak concentrations of roughly 100-120 mg/L after a 500 mg dose.
    # ------------------------------------------------------------------
    propSd <- 0.1341641
    label("Proportional residual error (fraction)")
    # sqrt(0.018); Table 2 Sigma epsilon[PROP] = 0.018 variance (RSE 9.3%,
    # shrinkage 11.9%, bootstrap median 0.018, 95% CI 0.012-0.027) -> 13.4% CV

    addSd <- 6.1721148
    label("Additive residual error (mg/L)")
    # sqrt(38.095); Table 2 Sigma epsilon[ADD] = 38.095 variance (RSE 16.1%,
    # shrinkage 9.7%, bootstrap median 38.009, 95% CI 12.368-61.724)
  })
  model({
    # 1. Clearance arms.
    #    Non-CRRT subjects: additive non-renal + renal decomposition, with the
    #    renal arm scaled by the dimensionless ratio (CCR / 54 mL/min).
    #    CRRT subjects: a single estimated total clearance.
    cl_nonren <- exp(lcl_nonren)
    cl_renal  <- exp(lcl_renal) * (CRCL / 54)
    cl_crrt   <- exp(lcl_crrt)

    # 2. Piecewise typical clearance. RRT_CRRT_STATUS is a full switch: when it
    #    is 1 the whole non-renal + renal expression is replaced, so a CRRT
    #    subject carries no separate non-renal term (Wu 2024 final covariate
    #    model equation, brace form).
    tvcl <- RRT_CRRT_STATUS * cl_crrt +
      (1 - RRT_CRRT_STATUS) * (cl_nonren + cl_renal)

    # 3. Individual parameters. The single reported eta[CL] applies to total
    #    clearance in both branches, not to the individual arms.
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # 4. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. Two-compartment disposition. Daptomycin is given as a 30-min IV
    #    infusion directly into the central compartment; the infusion duration
    #    belongs in the event table, not the model.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 6. Observation. Dose in mg and volumes in L give central / vc in mg/L,
    #    which is the ug/mL of the paper.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
