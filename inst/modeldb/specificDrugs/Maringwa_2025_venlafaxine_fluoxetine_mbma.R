Maringwa_2025_venlafaxine_fluoxetine_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis dose-response model for the mean change from baseline in the Hamilton Depression Rating (HAMD) scale in adults with major depressive disorder, fit jointly to venlafaxine and fluoxetine study-arm summary means from 16 placebo-controlled trials published 1987-2014 (43 arms, 3,432 patients) by weighted nonlinear regression (R gnls). Venlafaxine follows an Emax relationship in total daily dose; fluoxetine was not shown to have a dose-response over the studied 20-60 mg/day range and instead carries a single constant shift versus placebo applied to every active fluoxetine arm. Both drug effects are multiplied by the shared term 1 + 0.0986 * (SCORE_HAMD - 25), so the larger the arm's mean baseline HAMD score, the larger the expected drug effect. Placebo response is UNSTRUCTURED: a separate fixed effect was estimated for each of the 16 trials (mean -8.0, range -12.0 to -3.0) and this file encodes the typical value of -8 that the paper used for every published model prediction. Per-arm reported standard errors of the change were used as fixed weights and the estimation-scale sigma was held at 1, so no residual variance was estimated and no between-trial random effect exists. IMPORTANT: the two drug-effect values in the source's Table 2 are printed against TRANSPOSED row labels and the source's displayed final equation prints the baseline centring as (B - 23) rather than the fitted (B - 25); this file encodes the parameterisation that reproduces the paper's own Figures 2-4 and its own published predicted differences from placebo, which the table and equation as printed do not (see the vignette Errata for the full arithmetic). Suitable simulation scope is the study-arm mean HAMD change from baseline at a trial's primary analysis timepoint, NOT individual-patient responses; there is no PK layer and no time course."

  reference <- paste(
    "Maringwa J, Diderichsen PM, Valiathan C.",
    "Partial Residual Plots as an Integrated Model Diagnostic Tool in",
    "Model-Based Meta-Analysis.",
    "Clin Pharmacol Ther. 2025 Jan;117(1):153-159.",
    "doi:10.1002/cpt.3418.",
    "Model structure is Equations 1 and 2 (Methods, Model formulation) and the",
    "gnls fitting script in Table S1 of the supplement; parameter estimates are",
    "in Table 2; the 43-arm study-level dataset is in Table S1.",
    sep = " "
  )

  vignette <- "Maringwa_2025_venlafaxine_fluoxetine_mbma"

  # Algebraic MBMA dose-response model: no rxode2 dose events are consumed (the
  # per-arm total daily dose enters through the CONMED_<drug>_DOSE covariate
  # columns) and the output is a change from baseline in a clinical rating
  # scale rather than a drug concentration. The `units` entries follow the
  # placeholder convention already used by Mandema_2005_gemcabene_mbma and
  # Mercier_2014_tramadol_tapentadol_mbma so that checkModelConventions() sees
  # a parseable dosing / concentration pair.
  units <- list(
    time          = "week (placeholder; the model is a time-independent per-arm dose-response evaluated once per study arm at that trial's primary analysis timepoint, not a time course. Maringwa 2025 Methods: the HAMD change from baseline 'was analyzed at the time of primary analysis of each study (primary timepoint)' and no treatment-duration effect is modelled)",
    dosing        = "mg/day (per-arm TOTAL DAILY dose of venlafaxine or fluoxetine, supplied through the CONMED_VENLAFAXINE_DOSE / CONMED_FLUOXETINE_DOSE covariate columns, NOT as rxode2 dose events)",
    concentration = "score/score (arm-mean CHANGE FROM BASELINE in the 17-item Hamilton Depression Rating scale; negative values are improvements, e.g. Cc = -11.3 is an 11.3-point HAMD reduction. Output Cc is NOT a drug concentration; the slash satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    CONMED_VENLAFAXINE_DOSE = list(
      description        = "Per-arm total daily venlafaxine dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a venlafaxine arm (placebo arms and fluoxetine arms). Maringwa 2025 Methods: venlafaxine doses ranged from 25 to 375 mg/day across ten trials. The per-arm values tabulated in Table S1 are 25, 75, 150, 182, 200, 225 and 375 mg/day; where a trial randomised patients to a titrated range (e.g. 'venlafaxine 75-225 mg/day') the supplement records the TOP of the range as the arm dose, so this column carries the target daily dose of the arm rather than an achieved average. Drives the Emax dose-response term emax_venlafaxine * dose / (ed50_venlafaxine + dose) with ED50 = 29.1 mg/day, which means every studied arm except the 25 mg/day arm sits above the ED50 and the dose-response is already close to its plateau across most of the studied range. Member of the CONMED_<drug>_DOSE family; venlafaxine is spelled out in full per that family's rule for drug names with no established unambiguous 3-4 letter abbreviation.",
      source_name        = "drug1.dose where drug1 == 'venlafaxine' (Maringwa 2025 Table S1); d_ijk in Equation 2"
    ),
    CONMED_FLUOXETINE_DOSE = list(
      description        = "Per-arm total daily fluoxetine dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a fluoxetine arm (placebo arms and venlafaxine arms). Maringwa 2025 Methods: fluoxetine doses of 20, 40 and 60 mg/day were available across eight trials. NOTE the asymmetry with venlafaxine: no dose-response could be identified for fluoxetine over this narrow range, so the final model applies a single CONSTANT shift to any arm with a positive fluoxetine dose ('a constant drug effect combining all dose levels vs. placebo was identified for fluoxetine', Abstract). This column therefore acts as an arm INDICATOR in the model equation (only dose > 0 versus dose == 0 is used) even though the numeric dose is retained here for provenance and so that a downstream user can re-fit a dose-response term. Maringwa 2025 explicitly notes 'dose-normalization was only relevant for venlafaxine since the model for fluoxetine did not include a dose-response relationship.' Member of the CONMED_<drug>_DOSE family; fluoxetine is spelled out in full.",
      source_name        = "drug1.dose where drug1 == 'fluoxetine' (Maringwa 2025 Table S1); d_ijk in Equation 2"
    ),
    SCORE_HAMD = list(
      description        = "Per-arm mean baseline (pre-randomization) total score on the Hamilton Depression Rating scale.",
      units              = "(SCORE_HAMD units, 17-item score 0-52)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Study-ARM-level mean, not an individual-patient value: each row of Maringwa 2025 Table S1 carries one `baseline` value per arm, and placebo and active arms of the same trial can differ (e.g. trial 211 has 25.0 on placebo and 26.0 on the fluoxetine arm). Range across the 43 arms is 14.1 to 29.7; the unweighted mean across arms is 23.3, which the paper's own script rounds to 23 and uses as the reference baseline for every published prediction. CENTRING: the fitted model centres this covariate at 25, NOT at the 23 printed in the paper's displayed final equation -- the supplement's ModelFuncx() hard-codes `eff.bas = (1 + bbas * (baseline - 25))`. The 23/25 distinction is load-bearing: predictions AT baseline 23 carry the factor 1 + 0.0986 * (23 - 25) = 0.8028, which is exactly what reproduces the paper's published predicted differences from placebo. The effect is MULTIPLICATIVE on the whole drug effect and does not act on the placebo response, so it cancels from a placebo arm entirely. Because the multiplier 1 + 0.0986 * (SCORE_HAMD - 25) crosses zero at SCORE_HAMD = 14.9, arms with a mean baseline below about 15 have a drug effect of essentially zero or of reversed sign; the lowest arm in the dataset (14.1, Devanand 2005 placebo) sits just below that root, so extrapolation below a mean baseline of about 15 is outside the calibrated range and should be avoided.",
      source_name        = "baseline (Maringwa 2025 Table S1); B_ij in Equations 2 and the displayed final equation"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 3432L,
    n_studies      = 16L,
    n_arms         = 43L,
    age_range      = "Adults. Per-trial age distributions are not tabulated in Maringwa 2025 or its supplement; the constituent trials are adult major-depressive-disorder studies and one (Devanand 2005) is a late-life depression trial.",
    disease_state  = "Major depressive disorder. The pooled trials are placebo-controlled antidepressant efficacy studies; the clinical end point is the change from baseline in the Hamilton Depression Rating (HAMD) scale at each trial's primary analysis timepoint, chosen in line with the 2023 European Medicines Agency draft guideline on clinical investigation of medicinal products in the treatment of depression.",
    dose_range     = "venlafaxine 25-375 mg/day (ten trials); fluoxetine 20-60 mg/day (eight trials). Two of the 16 trials investigated both drugs.",
    baseline_hamd  = "Arm-mean baseline HAMD 14.1-29.7 across the 43 arms (unweighted mean 23.3, rounded to 23 in the paper). By drug (Maringwa 2025 Table 1): fluoxetine 20.8 (15, 26), venlafaxine 25.4 (23.5, 29.4), placebo arms of fluoxetine trials 20.5 (14.1, 25), placebo arms of venlafaxine trials 25 (23.7, 29.7).",
    placebo_response = "Trial-specific unstructured placebo response, one fixed effect per trial: mean -8.0, range -12.0 to -3.0 HAMD points (Maringwa 2025 Results). The individual 16 trial estimates are not tabulated in the paper.",
    n_patients_by_arm_type = c(venlafaxine = 1289L, fluoxetine = 982L, placebo = 1161L),
    data_sources   = "Literature data published between 1987 and 2014, selected from a larger proprietary antidepressant literature database that was assembled systematically following the Cochrane Handbook for Systematic Reviews of Interventions and PRISMA reporting items. The 16 selected trials, their per-arm doses, mean baseline HAMD scores, arm sizes, mean changes from baseline and change variances are tabulated in Table S1, each with a PubMed or ClinicalTrials.gov link. Venlafaxine and fluoxetine were chosen because they were among the most frequently reported drugs in placebo-controlled studies and their dose-ranging data permitted a dose-response assessment; the paper states 'the main purpose of this analysis was to illustrate methodology rather than compare the selected drugs.'",
    notes          = "Summary-level MBMA: each modelled observation is a study-ARM mean change from baseline in HAMD, not an individual patient value. The model is intended for simulating study-arm-mean HAMD change and is NOT suitable for individual-subject simulation. There is no PK layer and no time course. The paper's primary subject is the use of partial residual plots (PRPs) as an MBMA diagnostic; the antidepressant MBMA is its worked example, and it is a complete, original, fully parameterised model. The paper additionally reports an INTENTIONALLY MISSPECIFIED variant in which venlafaxine's Emax dose-response is replaced by a constant shift (Figures S1-S2), used solely to demonstrate that PRPs can reveal misspecification; per the replicate-author-structure policy that deliberate straw-man variant is NOT extracted, only the final model."
  )

  ini({
    # ========================================================================
    # Maringwa 2025 Equations 1 and 2 (Methods, Model formulation):
    #
    #   Y_ij = f(eo, d, B) + eps_ij
    #   f(eo, d, B) = eo_i + Emax_k * (1 + beta * (B_ij - Bbar)) * d_ijk
    #                        -----------------------------------------
    #                                   ED50_k + d_ijk
    #
    # and the supplement's verbatim fitting function (Table S1, chunk
    # `modeldevelopment`):
    #
    #   ModelFuncx <- function(drug1, drug1.dose, baseline, eo, emven, emflu,
    #                          edven = -Inf, edflu = -Inf, bbas = 0) {
    #     plc <- eo
    #     eff <- emven * drug1.dose / (exp(edven) + drug1.dose) *
    #              I(drug1 %in% c("venlafaxine")) +
    #            emflu * I(drug1 %in% c("fluoxetine"))
    #     eff.bas <- (1 + bbas * (baseline - 25))
    #     eff <- eff * eff.bas
    #     plc + eff
    #   }
    #
    # Y is the arm-mean CHANGE FROM BASELINE in HAMD, so every term below is
    # in HAMD points and a negative value is a clinical improvement.
    #
    # ------------------------------------------------------------------------
    # ERRATUM -- TWO PUBLISHED ERRORS, both proven from on-disk sources.
    #
    # (1) TABLE 2's TWO DRUG-EFFECT ROWS CARRY TRANSPOSED LABELS. Re-running
    #     the supplement's own gnls script verbatim against the supplement's
    #     own 43-arm Table S1 dataset returns:
    #
    #        emven = -4.4057 (SE 0.6477, 15% RSE)   <- venlafaxine Emax
    #        emflu = -1.7535 (SE 0.5884, 34% RSE)   <- fluoxetine constant shift
    #        edven =  3.3699 (SE 0.9149, 27% RSE)
    #        bbas  =  0.0986 (SE 0.0366, 37% RSE)
    #
    #     Table 2 prints exactly these four values AND exactly these four
    #     standard errors and %RSEs -- but attaches -4.41/0.648/15% to the row
    #     labelled "Constant shift fluoxetine" and -1.75/0.588/34% to the row
    #     labelled "Emax venlafaxine". The estimate/SE pairing is what proves
    #     this is a LABEL transposition and not a value error. The "Log ED50
    #     venlafaxine" and "Multiplicative effect of baseline HAMD" rows are
    #     correctly labelled.
    #
    # (2) THE DISPLAYED FINAL EQUATION PRINTS THE CENTRING AS (B - 23); THE
    #     FITTED MODEL CENTRES AT 25. The supplement hard-codes
    #     `(baseline - 25)`. 23 is a different quantity: it is
    #     round(mean(baseline)) over the 43 arms (23.33 -> 23), the reference
    #     baseline at which the paper DISPLAYS its predictions (Figure 1-4
    #     captions, Table 2 caption). The two were conflated when the equation
    #     was typeset.
    #
    #     The values encoded below are corroborated three independent ways:
    #       a. they are the re-fit estimates and SEs above;
    #       b. they reproduce the paper's own published predicted differences
    #          from placebo at a typical baseline HAMD of 23 -- fluoxetine at
    #          60 mg/day: -1.75 * (1 + 0.0986 * (23 - 25)) = -1.408 (paper:
    #          -1.41); venlafaxine at 375 mg/day: -4.41 * 375/(29.1 + 375) *
    #          0.8028 = -3.283 (paper: -3.28);
    #       c. they reproduce the model-prediction lines in Figures 2, 3 and 4
    #          (fluoxetine flat at about -9.4, venlafaxine reaching about
    #          -11.3 at 375 mg/day, both on a -8 placebo).
    #
    #     Encoding Table 2 as printed reproduces NONE of these and reverses
    #     the paper's own headline conclusion that venlafaxine outperforms
    #     fluoxetine. Full arithmetic is in the vignette's Errata section.
    # ========================================================================

    e0 <- -8
    label("Typical placebo response: arm-mean change from baseline in HAMD on placebo (HAMD points)")  # Maringwa 2025 Table 2 caption ("The overall placebo response was -8") and every figure caption. The placebo response is UNSTRUCTURED -- one fixed effect per trial, mean -8.0 and range -12.0 to -3.0 (Results) -- and this is the typical value the paper used for all published predictions. The paper's own script forms it as round(mean(trial-specific eo)); the unrounded re-fit mean is -7.7.

    # ------------------------------------------------------------------------
    # Venlafaxine: Emax in total daily dose. ED50 was estimated on the
    # log scale, so the ini() value is the printed log-scale estimate itself
    # rather than log(untransformed).
    # ------------------------------------------------------------------------
    emax_venlafaxine <- -4.41
    label("Venlafaxine maximal effect on the change from baseline in HAMD, at the centring baseline (HAMD points, signed; at infinite dose)")  # Maringwa 2025 Table 2, value -4.41 (SE 0.648, 15% RSE) -- printed on the row LABELLED "Constant shift fluoxetine"; see the erratum block above. Confirmed as emven by re-fitting the supplement's gnls script (-4.4057, SE 0.6477).

    led50_venlafaxine <- 3.37
    label("Venlafaxine log ED50 (log mg/day)")  # Maringwa 2025 Table 2, "Log ED50 venlafaxine" = 3.37 (SE 0.915, 27% RSE); the Results text gives the untransformed value as 29.1 mg/day and exp(3.37) = 29.08. Correctly labelled in Table 2.

    # ------------------------------------------------------------------------
    # Fluoxetine: no identifiable dose-response over the studied 20-60 mg/day
    # range, so a single constant shift applies to every active arm. This is
    # the ED50 -> 0 limit of Equation 2 (d / (ED50 + d) -> 1 for d > 0), which
    # is why the paper's displayed final equation writes it with an
    # indicator I(d > 0). The supplement's script sets edflu = -Inf by
    # default, giving exactly that limit.
    # ------------------------------------------------------------------------
    shift_fluoxetine <- -1.75
    label("Fluoxetine constant shift versus placebo on the change from baseline in HAMD, at the centring baseline (HAMD points, signed; applied to every active fluoxetine arm regardless of dose)")  # Maringwa 2025 Table 2, value -1.75 (SE 0.588, 34% RSE) -- printed on the row LABELLED "Emax venlafaxine"; see the erratum block above. Confirmed as emflu by re-fitting the supplement's gnls script (-1.7535, SE 0.5884).

    # ------------------------------------------------------------------------
    # Shared multiplicative baseline-severity effect on the drug effect. A
    # POSITIVE slope means a larger mean baseline HAMD gives a larger (more
    # negative) drug effect: Maringwa 2025 Results, "the larger the baseline,
    # the larger the negative changes (improvement in HAMD scores)".
    # ------------------------------------------------------------------------
    e_score_hamd_drug <- 0.0986
    label("Multiplicative effect of the arm-mean baseline HAMD score, centred at 25, on both drug effects (1/HAMD point)")  # Maringwa 2025 Table 2, "Multiplicative effect of baseline HAMD score on Emax/Constant shift" = 0.0986 (SE 0.0366, 37% RSE); likelihood-ratio test P = 0.0122 (Results). Centring constant 25 from the supplement's ModelFuncx(), NOT the 23 printed in the displayed final equation; see the erratum block above.

    # ========================================================================
    # Residual error and random effects.
    #
    # Maringwa 2025 Methods: "Reported (or imputed if missing) standard errors
    # (SEs) associated with eps_ij were used as fixed weights in the model and
    # hence no further residual variance was estimated in the model, thus
    # technically setting the residual standard deviation (sigma) to 1 in the
    # estimation." The supplement's gnls call carries
    # `weights = varFixed(~change.var)` and `control = list(sigma = 1)`.
    #
    # So the per-arm residual variance is 1 * change.var_ij, where
    # change.var_ij is the reported variance of the arm-mean change tabulated
    # in Table S1 (0.30 to 2.57 across the 43 arms, i.e. an SD of 0.55 to
    # 1.60 HAMD points, and roughly inversely proportional to arm size).
    # addSd is therefore FIXED at the estimation-scale sigma of 1, and the
    # per-arm sqrt(change.var) weighting is applied downstream in simulation
    # code -- the same pattern as Mercier_2014_tramadol_tapentadol_mbma,
    # Boucher_2018_naproxen_mbma and Vargo_2014_statins_ezetimibe_mbma.
    #
    # There is NO between-trial random effect and NO between-subject
    # variability: between-trial heterogeneity in placebo response is absorbed
    # by the 16 unstructured per-trial fixed effects rather than by a variance
    # component, and the paper reports no omega.
    # ========================================================================
    addSd <- fixed(1)
    label("Residual SD on the arm-mean change from baseline in HAMD (HAMD points); the estimation-scale sigma of 1, so the per-arm SD is sqrt(change.var) from Table S1 and that weighting is applied downstream")  # Maringwa 2025 Methods (sigma technically set to 1; reported per-arm SEs used as fixed weights) and the supplement's gnls control = list(sigma = 1) with weights = varFixed(~change.var).
  })

  model({
    # ---- Back-transform the log-scale potency parameter ---------------------
    ed50_venlafaxine <- exp(led50_venlafaxine)

    # ---- Venlafaxine Emax dose-response (Equation 2) ------------------------
    # A placebo or fluoxetine arm sets CONMED_VENLAFAXINE_DOSE to 0, which
    # collapses this term to 0 exactly.
    eff_venlafaxine <- emax_venlafaxine * CONMED_VENLAFAXINE_DOSE /
      (ed50_venlafaxine + CONMED_VENLAFAXINE_DOSE)

    # ---- Fluoxetine constant shift ------------------------------------------
    # The paper's script uses the arm indicator I(drug1 == "fluoxetine"); every
    # fluoxetine arm in Table S1 has a positive dose and every non-fluoxetine
    # arm has 0, so testing the dose column reproduces that indicator exactly
    # while keeping the numeric dose available for provenance.
    eff_fluoxetine <- shift_fluoxetine * (CONMED_FLUOXETINE_DOSE > 0)

    # ---- Shared multiplicative baseline-severity term ------------------------
    # Centred at 25 per the supplement's ModelFuncx(); see the ini() erratum.
    # Acts on the DRUG effect only, so it cancels entirely from a placebo arm.
    eff_baseline <- 1 + e_score_hamd_drug * (SCORE_HAMD - 25)

    # ---- Equation 1 ----------------------------------------------------------
    # Cc is overloaded here as the single-output observation per nlmixr2lib
    # convention; it is NOT a drug concentration but the study-arm mean change
    # from baseline in HAMD (Cc = -11.3 is an 11.3-point improvement).
    Cc <- e0 + (eff_venlafaxine + eff_fluoxetine) * eff_baseline

    Cc ~ add(addSd)
  })
}
