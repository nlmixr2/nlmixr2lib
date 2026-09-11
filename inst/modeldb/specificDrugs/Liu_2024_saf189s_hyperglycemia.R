Liu_2024_saf189s_hyperglycemia <- function() {
  description <- paste0(
    "Binomial logistic-regression exposure-safety model for ANY-GRADE ",
    "hyperglycemia in Chinese adults with ALK-positive or ",
    "ROS1-positive advanced non-small cell lung cancer treated with the ",
    "second-generation ALK/ROS1 tyrosine kinase inhibitor SAF-189s ",
    "(Liu 2024, n = 296, the safety analysis set of study SAF001, ",
    "NCT04237805, 20-210 mg orally once daily). The probability of the ",
    "event is expit(-9.430 + 1.259 * log(AUC_SAF189S)) where ",
    "AUC_SAF189S is the individual SAF-189s steady-state daily AUC in ",
    "ng*h/mL computed at the patient's FIRST dose level and the log is ",
    "a natural log, so the printed odds ratio 3.521 is per e-fold. ",
    "Hyperglycemia is the dominant toxicity of this drug (165 of 296 ",
    "patients, 55.74%, any grade) and is mechanistically expected: ALK ",
    "belongs to the insulin-receptor tyrosine-kinase superfamily, and ",
    "SAF-189s was designed from the ceritinib scaffold with an ",
    "insulin-receptor IC50 of 0.8 nM against ceritinib's 7 nM, so ",
    "insulin resistance follows from on-target off-tumour inhibition. ",
    "There is no PK layer and no ODE: the exposure metric is supplied ",
    "as a data column, derived in the source analysis from the ",
    "individual post hoc parameters of the companion population PK ",
    "model, packaged as Liu_2024_saf189s. No between-subject random ",
    "effect and no residual error are estimated (Bernoulli ",
    "likelihood). The intercept is not printed by the paper and was ",
    "recovered by digitising the fitted curve of Figure 10A; see the ",
    "in-file note on logit_ref and the vignette Errata. Companion ",
    "grade >= 2 model in Liu_2024_saf189s_hyperglycemia_grade2; six ",
    "companion models in the Liu_2024_saf189s_* family."
  )
  reference <- paste(
    "Liu Y, Tan Y, Hu L, Li J, Yang J, Diao L, Yang J.",
    "Population pharmacokinetics and exposure-response analyses of",
    "SAF-189s in Chinese patients with ALK+/ROS1+ non-small cell lung",
    "cancer.",
    "Front Pharmacol. 2024;15:1418549.",
    "doi:10.3389/fphar.2024.1418549.",
    "Individual SAF-189s steady-state AUC values derive from the companion",
    "population pharmacokinetic model reported in the same paper; see",
    "modellib('Liu_2024_saf189s').",
    sep = " "
  )
  vignette <- "Liu_2024_saf189s"
  units <- list(
    time          = "n/a (static landmark exposure-safety regression; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the AUC_SAF189S covariate column)",
    concentration = "prob_hyperglycemia (probability of any-grade hyperglycemia, 0-1; also logit_hyperglycemia)"
  )

  covariateData <- list(
    AUC_SAF189S = list(
      description        = "Individual SAF-189s area under the plasma concentration-time curve over the 24 h dosing interval at steady state (AUCss), per subject. Supplied as data: this model has no PK layer, and the source analysis used the individual post hoc parameters of the companion population PK model together with the patient's FIRST dose level (Liu 2024 Methods, E-R analysis).",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) exposure, at STEADY STATE, computed at the",
        "FIRST dose the patient received. The first-dose basis is stated",
        "outright by the Figure 9 and Figure 10 axis labels ('AUCss",
        "(ng*h/mL) ... by first dose') and matters because SAF-189s",
        "doses were reduced on toxicity: an AUCss recomputed at the most",
        "prevalent dose would be systematically lower in exactly the",
        "patients who had events. Note that the companion",
        "exposure-EFFICACY model Liu_2024_saf189s_orr uses a different",
        "metric on a different basis (Cmin,ss at the most prevalent",
        "dose), so the two exposure columns are not interchangeable.",
        "Enters on the NATURAL LOG scale: Liu 2024 Results states that",
        "'a linear effect of log AUCss described the data well', and",
        "digitising the Figure 10A fitted curve returns a slope of",
        "1.2596 per natural-log unit against the printed odds ratio",
        "3.521 (exp(1.2596) = 3.524), with the logit linear in",
        "log(AUCss) to R^2 = 0.99999. The paper's own statement that the",
        "odds ratio corresponds to 'an increase in SAF-189s exposure of",
        "1 ng*h*mL-1*d-1' is an erratum -- an odds ratio of 3.5 per",
        "1 ng*h/mL over a range spanning thousands of ng*h/mL is",
        "arithmetically impossible, and the digitisation shows the unit",
        "is one natural-log unit.",
        "Units are load-bearing: the intercept absorbs the unit choice,",
        "so supplying ug*h/mL instead of ng*h/mL would shift the logit",
        "by 1.259 * log(1000) = 8.70.",
        "For calibration, the analysis set's median AUCss is",
        "2,233 ng*h/mL, the geometric mean at the recommended 160 mg",
        "once-daily dose is 2,374 ng*h/mL, and the Figure 10 x-axis",
        "spans 174.0-8,338.1 ng*h/mL."
      ),
      source_name        = "AUCss estimated by the first dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 296L,
    n_studies      = 1L,
    n_observations = "296 binary event records (one per patient; first-occurrence landmark analysis, no repeated measures)",
    age_range      = "SAF001: phase I median 51 years (28.0-68.0), phase II median 54.1 years (20.0-84.0) (Liu 2024 Supplementary Table 1)",
    weight_range   = "median 63.2 kg, range 37.3-92.5 across the PK-evaluable cohort",
    sex_female_pct = 51.6,
    race_ethnicity = c(Asian = 100),
    disease_state  = "ALK-positive or ROS1-positive advanced non-small cell lung cancer; hyperglycemia graded by CTCAE",
    dose_range     = "SAF-189s 20, 40, 80, 120, 160 or 210 mg orally once daily in 21-day cycles",
    regions        = "China",
    notes          = paste0(
      "The exposure-safety analysis set is the 296 SAF001 patients with ",
      "safety follow-up (Liu 2024 Results, Exposure-safety analysis), a ",
      "different set from both the 317 PK-evaluable subjects and the ",
      "244 phase II efficacy patients. Observed outcome: 165 of 296 ",
      "(55.74%) had any-grade hyperglycemia and 82 (27.70%) had grade ",
      ">= 2 (Liu 2024 Table 5). That any-grade rate is high for the ",
      "class -- crizotinib and alectinib have no effect on glucose ",
      "tolerance and ceritinib's first-in-human incidence was 49%."
    )
  )

  ini({
    # ==================================================================
    # Liu 2024 reports the SLOPE of this regression but not its
    # intercept: Results, Exposure-safety analysis gives only "the ORs
    # corresponding to an increase in SAF-189s exposure ... for
    # any-grade hyperglycemia ... were 3.521 (95% CI: 2.103-5.897;
    # p < 0.001)". There is no coefficient table for any of the three
    # exposure-safety fits, in the paper or in either supplementary
    # document.
    #
    # The intercept was therefore recovered by DIGITISING the fitted
    # blue curve of Figure 10A. Method and cross-checks:
    #
    #   - Axis calibration used the two unambiguous outer tick labels,
    #     174.0126 and 8338.1020 ng*h/mL. Validation of that
    #     calibration: the figure's own median-quartile dashed line
    #     then lands at 2,231 ng*h/mL against the 2,233 ng*h/mL the
    #     Discussion prints for the same quantity (0.07% agreement).
    #   - The curve is linear in logit versus log(AUCss) to
    #     R^2 = 0.99999 over 474 digitised pixel columns.
    #   - The FREE-fit slope is 1.2596 per natural-log unit, i.e. an
    #     odds ratio of exp(1.2596) = 3.524 against the printed 3.521 --
    #     a 0.09% agreement. Recovering the published slope to four
    #     significant figures from the figure alone is what licenses
    #     reading the intercept off the same curve.
    #   - Holding the slope at the printed log(3.521) and refitting the
    #     intercept alone gives -9.430.
    #   - INDEPENDENT CORROBORATION, from a source that is not the
    #     figure: at the analysis set's median AUCss of 2,233 ng*h/mL
    #     this pair predicts P(event) = 0.569, against the observed
    #     any-grade hyperglycemia incidence of 165/296 = 55.74% counted
    #     in Table 5. A logistic MLE reproduces the observed event
    #     proportion near the centre of the exposure distribution.
    #
    # The slope below is the PRINTED value; only the intercept is
    # figure-derived. See the vignette Errata.
    # ==================================================================

    # ----- Logit intercept -----
    # digitised from Liu 2024 Figure 10A (fitted curve); NOT printed in the paper.
    # See the block comment above for the method and the three cross-checks.
    logit_ref <- -9.430 ; label("Logit of the probability of any-grade hyperglycemia at AUC_SAF189S = 1 ng*h/mL (unitless logit)")  # figure-derived: Liu 2024 Figure 10A fitted curve, intercept recovered with the slope held at the printed log(3.521); cross-checks P(event) = 0.569 at the median 2,233 ng*h/mL vs the observed 55.74% in Table 5

    # ----- Exposure slope on the logit -----
    e_auc_saf189s_logit <- log(3.521) ; label("Log-odds of any-grade hyperglycemia per e-fold increase in SAF-189s steady-state daily AUC (unitless logit)")  # Liu 2024 Results, Exposure-safety analysis: OR 3.521 (95% CI 2.103-5.897), p < 0.001; free-fit digitisation of Figure 10A independently returns OR 3.524

    # ----- No between-subject variability, no residual error -----
    # Bernoulli-likelihood logistic regression; no omega or sigma is
    # estimated. rxode2 requires an observation declaration, so the
    # deterministic probability is emitted with a tiny placeholder
    # additive residual, mirroring the
    # Chen_2021_lorlatinib_teae_grade3.R pattern. This does not perturb
    # the predicted probability; see the vignette Assumptions and
    # deviations section.
    addSd_prob_hyperglycemia <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    # ------------------------------------------------------------------
    # Linear predictor. The exposure is NOT centred: Liu 2024 fits the
    # raw log-transformed AUCss, so logit_ref is the logit at
    # AUC_SAF189S = 1 ng*h/mL rather than at any clinically meaningful
    # patient.
    # ------------------------------------------------------------------
    logit_hyperglycemia <- logit_ref + e_auc_saf189s_logit * log(AUC_SAF189S)

    prob_hyperglycemia <- expit(logit_hyperglycemia)

    # ------------------------------------------------------------------
    # Observation. Deterministic probability of any-grade hyperglycemia.
    # Downstream callers can sample binary outcomes with
    # rbinom(n, 1, prob_hyperglycemia) on the rxSolve output.
    # ------------------------------------------------------------------
    prob_hyperglycemia ~ add(addSd_prob_hyperglycemia)
  })
}
