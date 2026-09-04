Kemal_2026_nemtabrutinib_bor <- function() {
  description <- "Logistic-regression exposure-efficacy model for best overall response (BOR) to nemtabrutinib monotherapy in patients with CLL/SLL, from Kemal 2026. BOR (investigator-assessed partial or complete response per iwCLL 2018 criteria, attained at any time on treatment) is modelled on the logit scale as a linear function of the individual average on-treatment concentration (Cavg, a post-hoc exposure metric from the companion population PK model) plus a saturable Emax-type function of time on treatment. Fitted with glm in R, not in NONMEM, on the CLL/SLL subset (n = 288) of the population PK analysis set. Sister model files from the same paper: modellib('Kemal_2026_nemtabrutinib') for the population PK model that generates Cavg, and modellib('Kemal_2026_nemtabrutinib_ae') / modellib('Kemal_2026_nemtabrutinib_hypertension') for the two exposure-safety models."
  reference <- paste(
    "Kemal CC, Zweers TJ, Krekels EHJ, Chatterjee MS.",
    "Population Pharmacokinetic Modeling and Exposure-Response Analyses of",
    "Nemtabrutinib in Patients With Hematologic Malignancies.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(5). doi:10.1002/psp4.70257.",
    "Exposure-response methods are in Methods section 2.4;",
    "the coefficients this model carries are Table S4 of the supplement,",
    "and the fitted curves are drawn in Figure 3 (p. 9).",
    "Sister model files from the same paper:",
    "modellib('Kemal_2026_nemtabrutinib'),",
    "modellib('Kemal_2026_nemtabrutinib_ae'),",
    "modellib('Kemal_2026_nemtabrutinib_hypertension').",
    sep = " "
  )
  vignette <- "Kemal_2026_nemtabrutinib"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  # ------------------------------------------------------------------
  # SOURCING NOTE
  #
  # Every coefficient below is PRINTED, in Table S4 of the supplement,
  # together with the model form stated verbatim in the Table S4
  # footnote:
  #
  #   logit(p) = b0 + b1 * Cavg
  #                 + b2 * [follow up time / (ET50 + follow up time)]
  #
  # ET50 was not estimated jointly: the footnote says it "was determined
  # by optimizing the AIC of the fitted glm model", with the lowest AIC
  # (262.34) at ET50 = 200 days, "which was used in the exposure-efficacy
  # analysis". It is therefore encoded as fixed(200) rather than as a
  # free parameter.
  #
  # Nothing here is centred: Table S4 reports the raw glm intercept, so
  # logitbor is the log-odds of response at Cavg = 0 and zero time on
  # treatment. That point is far outside the observed data (the CLL/SLL
  # cohort's Cavg runs to about 2000 ng/mL and the plotted follow-up
  # scenarios are 180 and 360 days) and is not clinically interpretable
  # on its own; it is kept uncentred so the ini() values match the
  # printed table exactly.
  #
  # The transcription is self-validating against Figure 3, which draws
  # the fitted curve at two follow-up times. Substituting the Table S4
  # values gives, at Cavg = 0 / 1000 / 2000 ng/mL, response probabilities
  # of 0.037 / 0.206 / 0.635 at 180 days and 0.146 / 0.534 / 0.885 at
  # 360 days - which is what Figure 3's orange and brown curves read.
  # The vignette asserts this against curve positions digitized from the
  # published figure.
  # ------------------------------------------------------------------

  covariateData <- list(
    CSS_NEMTA = list(
      description        = "Individual average on-treatment plasma concentration of nemtabrutinib (Cavg).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A post-hoc exposure metric, not an observation. Kemal 2026 (Methods 2.3 and 2.4) simulated each participant's concentration-time profile from the companion population PK model using that participant's own dosing history - including dose interruptions and dose reductions - and their individual post hoc PK parameter estimates, then computed Cavg as the cumulative on-treatment AUC divided by the treatment duration. Generate it for simulation with modellib('Kemal_2026_nemtabrutinib'). Enters this model UNCENTRED, matching the raw glm intercept printed in Table S4. Observed distribution in the CLL/SLL exposure-efficacy cohort (Figure 3 reference boxplots): median about 600 ng/mL at 45 mg (n = 16), about 750 ng/mL at 65 mg (n = 203), and about 850 ng/mL at 80 mg (n = 54), with the pooled range extending to about 2000 ng/mL.",
      source_name        = "Cavg"
    ),
    T_TRT = list(
      description        = "Time on nemtabrutinib treatment (the paper's 'follow up time' / 'time on treatment').",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Carried continuously, NOT dichotomized: Kemal 2026 enters it through the saturable term T_TRT / (ET50 + T_TRT) with ET50 = 200 days, because 'the effect of time on response rate is only expected to persist initially' (Methods 2.4). The paper motivates the saturable form from the clinical observation that the proportion of responders to BTK inhibitors rises rapidly over the first 3-6 months and more slowly through about 12 months. Figure 3 draws the fitted relationship at two values of this covariate, 180 and 360 days. For patients who dropped out, Kemal 2026 assessed whether they had attained PR or CR and computed their exposure over the time they remained in the study, so this column is the time actually on treatment rather than the planned course.",
      source_name        = "follow up time"
    )
  )

  # Covariates Kemal 2026 names as expected predictors of efficacy but
  # that could NOT be evaluated, so no coefficient exists to carry.
  # Documentation only; checkModelConventions() does not require these to
  # be referenced in model().
  covariatesDataExcluded <- list(
    N_PRIORTHER = list(
      description        = "Number of prior lines of anticancer therapy before starting nemtabrutinib.",
      units              = "(count)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Named in Results 3.5 as a covariate 'expected to be relevant predictors of efficacy', but not estimable in this analysis because ALL patients in the cohort had received prior treatment, leaving no untreated contrast. No coefficient is reported anywhere in the paper or supplement.",
      source_name        = "prior lines of treatment"
    ),
    SNP_TP53_MUT = list(
      description        = "TP53 tumour-suppressor gene mutation status.",
      units              = "(binary)",
      type               = "binary",
      reference_category = NULL,
      notes              = "Named in Results 3.5 and again in the Discussion as an expected prognostic predictor, but not evaluable here because 'many of them had missing TP53 gene mutation status'. Kemal 2026 states that TP53 aberrations will be explored once phase 3 data accrue. No coefficient is reported.",
      source_name        = "TP53 gene mutation status"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 288L,
    n_studies      = 2L,
    age_range      = "25-89 years (pooled population PK analysis set)",
    weight_range   = "41.2-147 kg (pooled population PK analysis set)",
    disease_state  = "Relapsed / refractory chronic lymphocytic leukemia or small lymphocytic lymphoma (CLL/SLL). This is the CLL/SLL subset of the 578-patient population PK analysis set; the two exposure-safety models from the same paper instead use all 578 patients regardless of primary diagnosis.",
    dose_range     = "5-80 mg nemtabrutinib once daily orally; the doses received by most patients, and the ones drawn as reference boxplots in Figure 3, are 45, 65 and 80 mg (n = 16, 203 and 54 respectively). 65 mg daily is the recommended phase 2 dose.",
    endpoint       = "Best overall response (BOR): investigator-assessed partial response (PR) or complete response (CR) per iwCLL 2018 criteria, attained at any time throughout treatment. Methods 2.4 defines responders only as those attaining PR or CR; the assessor and the response criteria are stated nowhere in the running text and are recoverable ONLY from the Figure 3 y-axis label, which reads 'BOR (INV per IWCLL2018 Criteria)' - text that sits inside the vector figure and so does not appear in any text extraction of the PDF. In BELLWAVE-003, CLL/SLL participants treated with 65 mg showed an objective response rate of 35.0% at a median follow-up of 9.1 months; no PR or CR was observed at doses of 30 mg or below.",
    notes          = "The exposure-efficacy analysis set is the largest single-indication subset of the population PK analysis population: patients with a primary CLL/SLL diagnosis, available response data and available PK, who completed at least one treatment cycle. Patients who dropped out were retained, scored on whether they had attained PR or CR, with exposure computed over the time they remained in the study. Fitted by multivariate logistic regression in R v4.3.3 with a final AIC of 262.34."
  )

  ini({
    # Every value below is printed in Table S4 of the Kemal 2026
    # supplement. No random effects and no residual-error parameter: the
    # source fits a fixed-effects logistic regression with glm, one row
    # per patient, so all the stochasticity lives in the Bernoulli
    # endpoint itself.
    logitbor <- -7.3978; label("Logit of the probability of best overall response at Cavg = 0 and zero time on treatment (unitless)")  # Table S4, row "Intercept" (SE 0.9854; Z -7.508; p < 0.001)

    e_css_nemta_logitbor <- 0.0019; label("Effect of average on-treatment nemtabrutinib concentration on logit(BOR), per ng/mL")       # Table S4, row "Cavg" (SE 0.0005; Z 3.636; p < 0.001)
    e_ttrt_logitbor      <- 8.7637; label("Maximum effect of time on treatment on logit(BOR), approached as time on treatment >> ET50") # Table S4, row "follow up time / (ET50 + follow up time)" (SE 1.1822; Z 7.413; p < 0.001)

    # ET50 is a grid-searched constant, not a jointly estimated
    # parameter: the Table S4 footnote reports that AIC was minimised
    # (262.34) at 200 days, and that this value "was used in the
    # exposure-efficacy analysis". No SE is reported for it.
    et50 <- fixed(200); label("Time on treatment giving half the maximal time effect, ET50 (days)")                                     # Table S4 footnote (AIC-optimised, held fixed)
  })

  model({
    # ----------------------------------------------------------------
    # 1. Saturable time-on-treatment term. Zero at treatment start and
    #    approaching 1 for times far beyond ET50, so that the response
    #    rate rises quickly early and then flattens - the shape Kemal
    #    2026 Methods 2.4 argues for from the BTK-inhibitor literature.
    # ----------------------------------------------------------------
    ttrtEffect <- T_TRT / (et50 + T_TRT)

    # ----------------------------------------------------------------
    # 2. Linear predictor on the logit scale, exactly as written in the
    #    Table S4 footnote. Cavg is UNCENTRED so that logitbor is the
    #    raw glm intercept printed in the table.
    # ----------------------------------------------------------------
    lpbor <-
      logitbor +
      e_css_nemta_logitbor * CSS_NEMTA +
      e_ttrt_logitbor * ttrtEffect

    # ----------------------------------------------------------------
    # 3. Probability of response, and the Bernoulli endpoint. pbor is
    #    exposed as an output column so the fitted probabilities can be
    #    compared against Figure 3 directly, without simulating draws.
    # ----------------------------------------------------------------
    pbor <- expit(lpbor)

    bor ~ dbinom(1, pbor)
  })
}
