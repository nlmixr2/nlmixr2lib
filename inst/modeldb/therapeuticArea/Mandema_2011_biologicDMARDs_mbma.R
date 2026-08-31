Mandema_2011_biologicDMARDs_mbma <- function() {
  description <- paste0(
    "MBMA. Dose-response meta-analysis of the relative efficacy of nine ",
    "biologic disease-modifying antirheumatic drugs (DMARDs) in adults with ",
    "rheumatoid arthritis, fit jointly to the American College of ",
    "Rheumatology (ACR) 20, 50 and 70 responder counts of 50 randomized ",
    "controlled trials (21,529 patients, five mechanisms of action) by ",
    "binomial nonlinear regression in S-PLUS 6.2. Each per-arm responder ",
    "probability is the inverse logit of a trial placebo response plus an ",
    "Emax treatment effect in the arm's own normalized dose: an Emax shared ",
    "by all drugs of one mechanism of action (anti-TNF, anti-IL-1, ",
    "anti-CD28, anti-CD20, anti-IL-6) and an ED50 estimated per drug. The ",
    "three endpoints share one ED50 and differ only by a multiplicative ",
    "scaling of the treatment effect. Methotrexate given as a randomized ",
    "treatment carries a dose-independent mean effect, and initial ",
    "combination therapy of an anti-TNF with methotrexate is less than the ",
    "sum of its parts (interaction coefficient -0.32). Trials conducted in ",
    "East Asia showed a 33 percent greater treatment effect. Per-arm dose ",
    "enters through one CONMED_<drug>_DOSE column per treatment; there are no ",
    "rxode2 dose events and no PK layer. Outputs are the per-arm ACR 20, ",
    "50 and 70 responder probabilities for both reported populations ",
    "(inadequate responders to prior methotrexate, and methotrexate-naive) ",
    "together with the three log odds ratios versus the arm's own control. ",
    "Suitable simulation scope is per-arm (study-arm-mean) responder ",
    "probability, NOT individual-subject responses. The source publishes ",
    "no parameter table, so every value in ini() was recovered from the ",
    "published Tables 2-4, the Discussion dose/ED50 ratios and the Figure 2 ",
    "dose-response curves; the recovery and its accuracy checks are set out ",
    "in the validation vignette."
  )

  reference <- paste(
    "Mandema JW, Salinger DH, Baumgartner SW, Gibbs MA.",
    "A dose-response meta-analysis for quantifying relative efficacy of",
    "biologics in rheumatoid arthritis.",
    "Clin Pharmacol Ther. 2011;90(6):828-835.",
    "doi:10.1038/clpt.2011.256.",
    "The model equations are in the Methods section 'Analysis methodology'",
    "(p. 834). The Supplementary Data is a list of the 50 included trials",
    "and contains no parameter estimates.",
    sep = " "
  )

  vignette <- "Mandema_2011_biologicDMARDs_mbma"

  # Algebraic MBMA dose-response model: no rxode2 dose events are consumed (the
  # per-arm dose enters through the CONMED_<drug>_DOSE covariate columns) and the
  # outputs are dimensionless responder probabilities in [0, 1] rather than
  # drug concentrations. The `units` entries below follow the placeholder
  # convention already used by Vargo_2014_statins_ezetimibe_mbma and
  # Mandema_2011_anticoagulants_mbma so that checkModelConventions() sees a
  # parseable dosing / concentration pair.
  units <- list(
    time          = "week (placeholder; the model is a time-independent per-arm dose-response evaluated at each trial's pre-determined primary end point, not a time course)",
    dosing        = "mg/administration (per-arm dose supplied through the CONMED_<drug>_DOSE covariate columns, NOT as rxode2 dose events. Units are per-drug and are those of that drug's standard regimen: mg/kg q4w for abatacept and tocilizumab, mg q2w for adalimumab and certolizumab, mg/day for anakinra, mg twice weekly for etanercept, mg q4w for golimumab, mg/kg at weeks 0, 2, 6 then q8w for infliximab, mg at weeks 0 and 2 for rituximab, and mg/week for methotrexate -- see each covariateData entry)",
    concentration = "%/arm (per-arm ACR responder probability expressed as a percentage of patients in the study arm; the p_* outputs are on the 0-1 fraction scale. Output is NOT a drug concentration; the slash satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    # ---- Anti-TNF (tumor necrosis factor inhibitors) ------------------------
    CONMED_ADALIMUMAB_DOSE = list(
      description        = "Per-arm adalimumab dose per administration on a once-every-2-weeks regimen.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an adalimumab arm. Mandema 2011 Table 1 reports 6 trials, 2,136 patients, median 40 mg q2w with a range of 20-160 mg q2w. Drug class: anti-TNF. The source normalized each drug's dose to its standard regimen and found no statistically significant impact of dosage regimen after that normalization (once weekly vs. once every 2 weeks was tested for adalimumab), so a 40 mg every-week arm is entered as 80 mg q2w.",
      source_name        = "Dose (Mandema 2011 Table 1 adalimumab row; Table 2 suggested starting dose 40 mg q2w; Figure 2 adalimumab panel x-axis 'mg q2w')"
    ),
    CONMED_CERTOLIZUMAB_DOSE = list(
      description        = "Per-arm certolizumab pegol dose per administration on a once-every-2-weeks regimen.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a certolizumab arm. Mandema 2011 Table 1 reports 4 trials, 1,512 patients, median 200 mg q2w with a range of 200-400 mg q2w. Drug class: anti-TNF. Certolizumab is dosed the furthest up its own dose-response curve of any drug in the analysis; the source states it is given at more than 10 times its ED50 and suggests it 'may be dosed too high'.",
      source_name        = "Dose (Mandema 2011 Table 1 certolizumab row; Table 2 suggested starting dose 200 mg q2w; Figure 2 certolizumab panel x-axis 'mg q2w')"
    ),
    CONMED_ETANERCEPT_DOSE = list(
      description        = "Per-arm etanercept dose per administration on a twice-weekly regimen.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an etanercept arm. Mandema 2011 Table 1 reports 11 trials, 2,493 patients, median 25 mg twice a week with a range of 0.5-50 mg. Drug class: anti-TNF. Once-weekly and twice-weekly etanercept regimens were compared in the source and showed no statistically significant regimen effect after normalization to a weekly dose, so a 50 mg once-weekly arm is entered as 25 mg twice weekly.",
      source_name        = "Dose (Mandema 2011 Table 1 etanercept row; Table 2 suggested starting dose 25 mg biw; Figure 2 etanercept panel x-axis 'mg biw')"
    ),
    CONMED_GOLIMUMAB_DOSE = list(
      description        = "Per-arm golimumab dose per administration on a once-every-4-weeks regimen.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a golimumab arm. Mandema 2011 Table 1 reports 4 trials, 1,231 patients, median 100 mg q4w with a range of 50-200 mg q4w. Drug class: anti-TNF. Golimumab sits lowest on its dose-response curve of the five anti-TNFs (about 1.1 times its ED50 at the 50 mg q4w starting dose), which is why the source concludes it 'may be dosed too low' with respect to efficacy.",
      source_name        = "Dose (Mandema 2011 Table 1 golimumab row; Table 2 suggested starting dose 50 mg q4w; Figure 2 golimumab panel x-axis 'mg q4w')"
    ),
    CONMED_INFLIXIMAB_DOSE = list(
      description        = "Per-arm infliximab dose per administration, weight-normalized, on the weeks 0, 2, 6 then once-every-8-weeks regimen.",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an infliximab arm. Mandema 2011 Table 1 reports 7 trials, 2,179 patients, median 6 mg/kg with a range of 2-20 mg/kg on the weeks 0, 2, 6 then q8w regimen. Drug class: anti-TNF. A weight-normalized dose, NOT a mass dose; the matching ED50 (1.42 mg/kg) carries the same units. Body weight therefore enters an infliximab arm through the dose units rather than as a covariate effect.",
      source_name        = "Dose (Mandema 2011 Table 1 infliximab row; Table 2 suggested starting dose 3 mg/kg wk 0, 2, 6, q8w; Figure 2 infliximab panel x-axis 'mg/kg wk 0.2.6.q8w')"
    ),
    # ---- Anti-IL-1 ----------------------------------------------------------
    CONMED_ANAKINRA_DOSE = list(
      description        = "Per-arm anakinra dose per administration on a once-daily regimen.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an anakinra arm. Mandema 2011 Table 1 reports 3 trials, 949 patients, median 75 mg/day with a range of 3-162 mg/day. Anakinra is the only anti-IL-1 in the analysis, so that class's Emax is identified by anakinra alone; the source found anakinra to give the smallest treatment effect of all nine biologics.",
      source_name        = "Dose (Mandema 2011 Table 1 anakinra row; Table 2 suggested starting dose 100 mg qd; Figure 2 anakinra panel x-axis 'mg qd')"
    ),
    # ---- Anti-CD28 (T-cell costimulatory blocking agents) -------------------
    CONMED_ABATACEPT_DOSE = list(
      description        = "Per-arm abatacept dose per administration, weight-normalized, on a once-every-4-weeks regimen.",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an abatacept arm. Mandema 2011 Table 1 reports 6 trials, 1,242 patients, median 10 mg/kg q4w with a range of 0.5-10 mg/kg q4w. Abatacept is the only anti-CD28 agent in the analysis, so that class's Emax is identified by abatacept alone. A weight-normalized dose, NOT a mass dose; the matching ED50 (2.53 mg/kg) carries the same units. This is a per-arm dose in an MBMA, not a subject-level per-administration dose amount.",
      source_name        = "Dose (Mandema 2011 Table 1 abatacept row; Table 2 suggested starting dose 10 mg/kg q4w; Figure 2 abatacept panel x-axis 'mg/kg q4w')"
    ),
    # ---- Anti-CD20 (B-cell-depleting agents) --------------------------------
    CONMED_RITUXIMAB_DOSE = list(
      description        = "Per-arm rituximab dose per administration on the two-infusion weeks 0 and 2 regimen.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a rituximab arm. Mandema 2011 Table 1 reports 3 trials, 707 patients, median 1,000 mg with a range of 500-1,000 mg given at weeks 0 and 2. Rituximab is the only anti-CD20 agent in the analysis, so that class's Emax is identified by rituximab alone. The source found rituximab the most efficacious of the non-anti-TNF biologics, between adalimumab and etanercept.",
      source_name        = "Dose (Mandema 2011 Table 1 rituximab row; Table 2 suggested starting dose 1,000 mg wk 0, 2; Figure 2 rituximab panel x-axis 'mg wk 0.2')"
    ),
    # ---- Anti-IL-6 ----------------------------------------------------------
    CONMED_TOCILIZUMAB_DOSE = list(
      description        = "Per-arm tocilizumab dose per administration, weight-normalized, on a once-every-4-weeks regimen.",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a tocilizumab arm. Mandema 2011 Table 1 reports 7 trials, 2,377 patients, median 8 mg/kg q4w with a range of 2-8 mg/kg q4w. Tocilizumab is the only anti-IL-6 agent in the analysis. Its theoretical Emax was significantly greater than that of the anti-TNFs while its ED50 is above the studied dose range, so doubling the tocilizumab dose from 4 to 8 mg/kg buys far more response than doubling an anti-TNF dose -- the single steepest dose-response in the analysis. A weight-normalized dose, NOT a mass dose.",
      source_name        = "Dose (Mandema 2011 Table 1 tocilizumab row; Table 2 suggested starting dose 4 mg/kg q4w; Figure 2 tocilizumab panel x-axis 'mg/kg q4w')"
    ),
    # ---- Methotrexate as a RANDOMIZED treatment -----------------------------
    CONMED_MTX_DOSE = list(
      description        = "Per-arm methotrexate dose when methotrexate is given as a RANDOMIZED treatment (methotrexate-naive trials). 0 when methotrexate is background treatment or absent.",
      units              = "mg/week",
      type               = "continuous",
      reference_category = "0 (methotrexate not a randomized treatment; either absent, or present only as background treatment and therefore absorbed into the trial placebo response)",
      notes              = "Mandema 2011 Table 1 reports 8 trials, 1,726 patients, with randomized methotrexate titrated to a median 18 mg/week (range 15-20). ONLY the presence of randomized methotrexate is used by the model, not the magnitude: 'The treatment effect of MTX when given as a randomized treatment was estimated as a simple mean effect because MTX was titrated to a similar dose range in all the trials and no dose-response data were available.' Any value greater than 0 therefore produces the same methotrexate effect, and values outside the titrated 15-20 mg/week range are outside the source's calibration. Setting this column non-zero alongside an anti-TNF dose column reproduces the source's initial-combination-therapy arms and activates the interaction coefficient. This is a DISTINCT concept from the canonical CONMED_MTX (background concomitant methotrexate), which in this model carries NO effect because it is absorbed into the trial-specific placebo response, and from DOSE_MTX_MGM2 (body-surface-area-normalized oncology methotrexate dosing).",
      source_name        = "MTX dose 18 (15 to 20) mg qw (Mandema 2011 Table 1 MTX row)"
    ),
    # ---- Trial geography ----------------------------------------------------
    REGION_EASTASIA = list(
      description        = "Indicator that the trial was carried out in East Asia (1) rather than elsewhere in the world (0).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (trial carried out outside East Asia -- Australia, Europe or North America).",
      notes              = "Study-arm-level trial-location indicator, not a subject-level covariate. Country membership in this model is Japan and Taiwan: 'Three trials were carried out in Japan (two for tocilizumab and one for infliximab) and one in Taiwan (etanercept).' The source reports that in these trials 'the reported treatment effects were significantly (33% (95% CI, 12 to 54%)) greater than those reported in trials carried out elsewhere in the world', consistent across the three drugs evaluated, and that after accounting for it there was no significant remaining between-trial heterogeneity. The effect acts on the Emax of the dose-response, per the source's covariate equation Emax,i = Emax,class*(1 + theta_c*(X_ij - Xbar) + eta_Emax,i), so it scales the whole biologic treatment effect by 1.33. It is NOT applied to the randomized-methotrexate mean effect, which has no Emax. Distinct from the RACE_ASIAN* family: this records where the trial was run, not patient race, and the source tested only trial location.",
      source_name        = "primary geographic location = Asia (Mandema 2011 Methods 'Data sources' extraction list; Results paragraph on Asian trials)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Per-arm mean patient age. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Mandema 2011 Methods lists age among the extracted patient-population characteristics tested on the parameters of the dose-response relationship. Result: 'There was also no statistically significant impact of differences in background treatment (placebo, DMARD other than MTX, or MTX), age, disease duration, gender, missing data imputation method, or treatment duration.' No point estimate is reported, so no effect parameter is included in ini(). Per-drug across-trial medians and ranges are in Table 1 (overall median 52 years, range 46-59 across trials).",
      source_name        = "Age (years) (Mandema 2011 Table 1; Methods covariate list)"
    ),
    SEXF = list(
      description        = "Per-arm proportion female. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods ('gender') and reported as non-significant in the Results; no point estimate is given. Rheumatoid arthritis trial populations are predominantly female, but the source does not tabulate the per-trial proportion.",
      source_name        = "gender (Mandema 2011 Methods covariate list)"
    ),
    WT = list(
      description        = "Per-arm mean patient body weight. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods ('patient population characteristics such as age, weight, gender, and disease duration') and reported non-significant in the Discussion ('Differences in failed prior treatments, background treatment, trial duration, age, gender, weight, and disease duration did not significantly affect the treatment effect'). Note that three drugs here are dosed per kilogram (abatacept, infliximab, tocilizumab), so body weight enters those arms through the mg/kg dose units rather than as a covariate effect.",
      source_name        = "weight (Mandema 2011 Methods covariate list)"
    ),
    DIS_DURATION_RA = list(
      description        = "Per-arm median duration of rheumatoid arthritis at trial entry. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tabulated per drug in Mandema 2011 Table 1 (overall median 9 years, range 0.4-15 across trials) and screened per the Methods covariate list. Reported non-significant in the Results and Discussion; no point estimate is given. Disease duration separates the methotrexate-naive trials (median 1 year) from the methotrexate-inadequate-responder trials (median 9 years), so it is strongly confounded with the population that the two placebo-response sets in ini() already distinguish.",
      source_name        = "Disease duration (years) (Mandema 2011 Table 1; Methods covariate list)"
    ),
    PRIOR_TX_FAILURE_GROUP = list(
      description        = "Failed prior treatments, grouped by the source as nonbiologic DMARDs other than methotrexate, nonbiologic DMARDs including methotrexate, or biologics. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Tabulated per drug in Mandema 2011 Table 1 as 'Trials failed DMARD/MTX/anti-TNF' (9/36/5 overall) and screened per the Methods covariate list. Two specific findings are reported but neither reached the P < 0.01 acceptance criterion, so neither is encoded: (i) 'There was no statistically significant difference in the odds ratio (relative to placebo) for MTX-naive patients vs. patients who had shown inadequate response to prior MTX treatment'; (ii) in patients with an inadequate response to a prior anti-TNF, the golimumab treatment effect was 'estimated to be slightly smaller (15% relative decrease) but not significantly different from the treatment effect in patients with inadequate response to MTX'. The absolute responder percentage still differs between these populations because the placebo response differs, which the two placebo-response sets in ini() capture.",
      source_name        = "Trials failed DMARD/MTX/anti-TNF (Mandema 2011 Table 1; Methods covariate list)"
    ),
    BACKGROUND_TX_GROUP = list(
      description        = "Background (non-randomized) treatment in the trial: placebo alone, a nonbiologic DMARD other than methotrexate, or methotrexate. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Screened per Mandema 2011 Methods and reported non-significant in the Results. Background methotrexate is the reason the canonical CONMED_MTX indicator carries no parameter in this model: its effect is absorbed entirely into the trial-specific placebo response, which is why Table 2 is headed 'patients on background MTX treatment'. Only randomized methotrexate has an estimated effect, entered through CONMED_MTX_DOSE.",
      source_name        = "background treatment (Mandema 2011 Methods covariate list; Table 1 'Placebo' rows)"
    ),
    TRT_DURATION = list(
      description        = "Trial treatment duration to the pre-determined primary efficacy end point. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "week",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reported in Mandema 2011 Results: 34 trials had treatment periods of 22-30 weeks, 13 trials 12-16 weeks, and 4 trials had their primary efficacy end point at 1 year. Screened and reported non-significant; no point estimate is given. This is why the model has no time axis -- every arm contributes its response at its own primary end point.",
      source_name        = "treatment duration (Mandema 2011 Results 'Data'; Methods covariate list)"
    ),
    IMPUTATION_METHOD = list(
      description        = "Missing-data imputation method used by the source trial to derive the ACR responder percentage. Screened as a dose-response covariate but NOT retained in the final model.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Extracted per Mandema 2011 Methods ('All end-point information was extracted from the publications (including imputation methods)') and reported non-significant in the Results. Relevant because the certolizumab registration trials withdrew week-12/14 ACR 20 non-attainers and counted them as week-24 non-responders; the source examined and rejected this as an explanation of the large certolizumab effect.",
      source_name        = "missing data imputation method (Mandema 2011 Methods; Results)"
    ),
    PLACEBO_RESPONSE_MAGNITUDE = list(
      description        = "Magnitude of the trial's own placebo response, screened as a modifier of the treatment effect but NOT retained in the final model.",
      units              = "(fraction)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested explicitly and reported non-significant: 'There was no statistically significant remaining impact of the magnitude of the placebo response on the treatment effect, suggesting that the relative-odds model adequately normalizes the treatment effect for differences in placebo response across the trials.' This is the finding that licenses shifting the e0_* parameters in ini() to move to a different absolute-response setting WITHOUT changing the treatment effect. Observed placebo-response ranges across the 50 trials were 8.7-42% for ACR 20, 0-29% for ACR 50 and 0-16% for ACR 70.",
      source_name        = "magnitude of the placebo response (Mandema 2011 Results; Discussion)"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 21529L,
    n_studies            = 50L,
    n_treatments         = 10L,
    age_range            = "Adults with rheumatoid arthritis; across-trial median age 52 years (range of per-trial medians 46-59) per Mandema 2011 Table 1. Age was screened as a covariate and not retained.",
    weight_range         = "Not tabulated. Body weight was screened as a covariate and not retained; abatacept, infliximab and tocilizumab are dosed per kilogram so weight enters those arms through the mg/kg dose units.",
    sex_female_pct       = NA_real_,
    disease_state        = "Rheumatoid arthritis in adults. 36 trials enrolled patients with an inadequate response to methotrexate (all or a percentage of patients), 9 trials enrolled methotrexate-naive patients treated only with older nonbiologic DMARDs or with NSAIDs and oral corticosteroids, and 5 trials enrolled patients with an inadequate response to methotrexate and/or a biologic anti-TNF. Efficacy end points are the percentage of patients attaining an ACR 20, 50 or 70 response at each trial's pre-determined primary time point; all 50 trials reported ACR 20, 49 reported ACR 50 and 45 reported ACR 70. Across-trial median disease duration 9 years (range of per-trial medians 0.4-15). Relative safety and tolerability were not addressed by this analysis.",
    dose_range           = "Nine biologic DMARDs across five mechanisms of action plus randomized methotrexate; per-drug median (range) doses from Mandema 2011 Table 1 are abatacept 10 (0.5-10) mg/kg q4w, adalimumab 40 (20-160) mg q2w, anakinra 75 (3-162) mg/day, certolizumab 200 (200-400) mg q2w, etanercept 25 (0.5-50) mg twice weekly, golimumab 100 (50-200) mg q4w, infliximab 6 (2-20) mg/kg at weeks 0, 2, 6 then q8w, rituximab 1,000 (500-1,000) mg at weeks 0 and 2, tocilizumab 8 (2-8) mg/kg q4w, and methotrexate 18 (15-20) mg/week.",
    regions              = "Predominantly Australia, Europe and North America. Four trials were carried out in East Asia: three in Japan (two tocilizumab, one infliximab) and one in Taiwan (etanercept).",
    n_trials_acr20       = 50L,
    n_trials_acr50       = 49L,
    n_trials_acr70       = 45L,
    n_trials_placebo_ctl = 40L,
    mechanisms_of_action = "Anti-TNF (adalimumab, certolizumab pegol, etanercept, golimumab, infliximab); anti-IL-1 (anakinra); T-cell costimulatory blocking agent with anti-CD28 activity (abatacept); B-cell-depleting agent with anti-CD20 activity (rituximab); anti-IL-6 (tocilizumab). Ofatumumab, ocrelizumab, baminercept and ustekinumab appear in the Methods search list but no trial of theirs met the inclusion criteria, so they have no ED50 and are not columns in this model.",
    notes                = "Summary-level MBMA: the modelled observations were per-arm ACR responder COUNTS out of the arm sample size, assumed binomial, with a separate placebo-response fixed effect estimated for EACH of the 50 trials and a trial-level random effect on the ACR 20 and ACR 70 placebo responses relative to ACR 50. Because of those per-trial intercepts the identified response variable is the log odds ratio between the active and control arms, not the absolute responder rate; the source reports none of the 50 intercepts and none of the random-effect variances. Between-trial random effects on Emax and ED50 were tested and rejected ('the absence of significant between-trial heterogeneity with respect to Emax and ED50'), so the model carries NO random effects and NO between-subject variability. Within-arm correlation across the three ACR end points measured in the same arm was handled by a compound-symmetry structure that is a property of the estimation, not of the predictive model, and is therefore not encoded. Only 1 of the 50 trials directly compared two biologics (infliximab vs. abatacept), so the comparisons this model supports are indirect ones anchored on the shared placebo and methotrexate arms."
  )

  ini({
    # ========================================================================
    # PROVENANCE FOR THE WHOLE ini() BLOCK
    #
    # Mandema 2011 publishes its model equations (Methods, p. 834) but NO
    # parameter table, and its Supplementary Data is a list of the 50 included
    # trials only. Every value below was therefore recovered by inverting the
    # paper's own printed results:
    #
    #   * The placebo responses and the treatment effect of each drug at its
    #     suggested starting dose follow in closed form from Table 2 and the
    #     typical placebo responses printed beneath it (24% ACR 20, 8.6%
    #     ACR 50, 2.7% ACR 70), because g_k = logit(placebo + difference) -
    #     logit(placebo).
    #   * The anti-TNF Emax follows from those effects plus the dose/ED50
    #     ratios printed in the Discussion (golimumab 1.1, infliximab 2.1,
    #     adalimumab 2.2, etanercept 6.8 times the ED50). Each of the four
    #     ratios independently implies the SAME shared Emax to within 1.2%,
    #     which is what makes the recovery credible rather than a fit.
    #   * The four single-drug classes have only one published dose each, so
    #     their Emax was digitized from the Figure 2 dose-response curves and
    #     their ED50 was then set so that the curve passes exactly through the
    #     published Table 2 point. Digitizing the five anti-TNF panels the same
    #     way recovers the algebraically-known Emax to +2.2%, which is the
    #     measured accuracy of the digitization.
    #   * The methotrexate-naive placebo responses and the methotrexate effect
    #     follow from Table 3 plus the printed statement that methotrexate
    #     monotherapy gives 32, 24 and 12 percentage-point absolute
    #     differences from placebo; five Table 3 rows over-determine two
    #     unknowns per end point and agree to within 0.1 percentage points.
    #
    # The vignette reproduces every step and shows that the recovered
    # parameters regenerate Tables 2, 3 and 4 of the source. All parameters
    # are marked fixed(): none of them was estimated here, and none of them
    # should be re-estimated from the aggregate data this model describes.
    # ========================================================================

    # ------------------------------------------------------------------------
    # Placebo response, patients with an inadequate response to prior
    # methotrexate and receiving background methotrexate (Table 2 population).
    # The source parameterizes the placebo response as an ACR 50 fixed effect
    # per trial plus mean logit-scale shifts theta_acr20 and theta_acr70 for
    # the other two end points ("E0,ki = E0,i + h_i,k ... theta_acr20 and
    # theta_acr70 represent the mean shifts in response on the logit scale
    # from ACR 50"), which is the structure used here.
    #
    # These describe the TYPICAL trial. Observed placebo responses ranged from
    # 8.7 to 42% for ACR 20, 0 to 29% for ACR 50 and 0 to 16% for ACR 70. The
    # source found no significant effect of the placebo-response magnitude on
    # the treatment effect, so shifting e0_acr50 to move to a different
    # absolute-response setting leaves every log odds ratio unchanged.
    # ------------------------------------------------------------------------
    e0_acr50 <- fixed(-2.363483)
    label("Control-arm log-odds of an ACR 50 response, inadequate responders to prior methotrexate (logit scale; 8.6% response)")
    # Mandema 2011 p. 830: "a typical placebo response of 24% for ACR 20, 8.6%
    # for ACR 50, and 2.7% for ACR 70"; logit(0.086) = -2.363483.

    theta_acr20 <- fixed(1.210804)
    label("Mean logit-scale shift of the ACR 20 placebo response from the ACR 50 placebo response")
    # logit(0.240) - logit(0.086) = -1.152680 - (-2.363483) = 1.210804.

    theta_acr70 <- fixed(-1.221064)
    label("Mean logit-scale shift of the ACR 70 placebo response from the ACR 50 placebo response")
    # logit(0.027) - logit(0.086) = -3.584547 - (-2.363483) = -1.221064.

    # ------------------------------------------------------------------------
    # Placebo response, methotrexate-naive population (Table 3 / Table 4).
    # Recovered by solving, per end point, the two unknowns (placebo response,
    # methotrexate effect) against six printed constraints: the five Table 3
    # biologic-minus-methotrexate differences and the printed statement that
    # "The ACR 20, 50, and 70 response rates to MTX monotherapy in a typical
    # MTX-naive patient population were estimated to be 32, 24, and 12%,
    # expressed in terms of absolute difference from placebo." The solution
    # reproduces all five Table 3 rows to within 0.1 percentage points and is
    # slightly LARGER than the methotrexate-inadequate-responder placebo
    # response for every end point, independently confirming the source's
    # statement that "the absolute effect of an anti-TNF was slightly greater
    # in MTX-naive patients because of a slightly larger placebo response in
    # these patients."
    # ------------------------------------------------------------------------
    e0_acr50_naive <- fixed(-2.281447)
    label("Control-arm log-odds of an ACR 50 response, methotrexate-naive population (logit scale; 9.27% response)")
    # Recovered from Mandema 2011 Table 3 + the 32/24/12 percentage-point
    # methotrexate statement (p. 833); logit(0.092671) = -2.281447.

    theta_acr20_naive <- fixed(1.278057)
    label("Mean logit-scale shift of the ACR 20 placebo response from ACR 50, methotrexate-naive population")
    # logit(0.268275) - logit(0.092671) = 1.278057; recovered as above.

    theta_acr70_naive <- fixed(-1.225014)
    label("Mean logit-scale shift of the ACR 70 placebo response from ACR 50, methotrexate-naive population")
    # logit(0.029129) - logit(0.092671) = -1.225014; recovered as above.

    # ------------------------------------------------------------------------
    # End-point scaling of the TREATMENT EFFECT, relative to ACR 50.
    #
    # The source tested "differences in Emax and/or ED50 between ACR 20, 50,
    # and 70" and retained a difference in Emax only: "There were no
    # statistically significant differences in dose required to achieve 50% of
    # maximum effect (ED50) between ACR 20, 50, and 70." The scaling is
    # therefore multiplicative on the treatment effect at every dose, and is
    # shared by all mechanisms of action ("There was also no statistically
    # significant difference across the mechanisms of action with respect to
    # the difference between ACR 20, 50, and 70 responses").
    #
    # The values below are the least-squares scale factors over all 27 cells
    # of Table 2. They are constant to within 0.5% across the nine drugs,
    # which is itself the evidence that the end-point difference acts on Emax
    # (a difference in ED50 would make the ratio vary with dose/ED50, and the
    # nine drugs span 0.9 to 17 times their ED50). The source's prose reports
    # 11% (95% CI 6.2 to 17%) smaller for ACR 20 and 11% (2.1 to 19%) larger
    # for ACR 70; the Table 2 cells give 14.3% and 11.3%, both inside those
    # intervals. Table 2 is used here because it also regenerates Tables 3
    # and 4 -- see the vignette's Assumptions and deviations section.
    # ------------------------------------------------------------------------
    scale_acr20 <- fixed(0.857025)
    label("Multiplicative scaling of the treatment effect for ACR 20 relative to ACR 50 (unitless)")
    # Least squares over the 27 cells of Mandema 2011 Table 2.

    scale_acr70 <- fixed(1.113197)
    label("Multiplicative scaling of the treatment effect for ACR 70 relative to ACR 50 (unitless)")
    # Least squares over the 27 cells of Mandema 2011 Table 2.

    # ------------------------------------------------------------------------
    # Emax by MECHANISM OF ACTION, on the ACR 50 logit scale. "Emax,class is
    # the maximal drug effect, reflecting the theoretical maximal difference in
    # response between placebo and active treatment. A different Emax was
    # estimated for each of the five mechanisms of action. A similar Emax was
    # assumed for drugs with the same mechanism of action."
    # ------------------------------------------------------------------------
    emax_antitnf <- fixed(2.156558)
    label("Maximal ACR 50 treatment effect of the anti-TNF class (logit units)")
    # Mean of the four values implied independently by Table 2 and the
    # Discussion dose/ED50 ratios: golimumab 2.1510, infliximab 2.1600,
    # adalimumab 2.1706, etanercept 2.1446 (spread 1.2%).

    emax_antiil1 <- fixed(0.954)
    label("Maximal ACR 50 treatment effect of the anti-IL-1 class, i.e. anakinra (logit units)")
    # Digitized from the Mandema 2011 Figure 2 anakinra panel. Much smaller
    # than the anti-TNF Emax, matching "Of all the biologics, anakinra
    # provided the smallest treatment effect"; the source reports the
    # difference from anti-TNFs as not statistically significant, which with 3
    # anakinra trials over a narrow effective dose range is unsurprising.

    emax_cd28 <- fixed(1.747)
    label("Maximal ACR 50 treatment effect of the anti-CD28 class, i.e. abatacept (logit units)")
    # Digitized from the Mandema 2011 Figure 2 abatacept panel. The source
    # reports no statistically significant difference from the anti-TNF Emax.

    emax_cd20 <- fixed(1.869)
    label("Maximal ACR 50 treatment effect of the anti-CD20 class, i.e. rituximab (logit units)")
    # Digitized from the Mandema 2011 Figure 2 rituximab panel. The source
    # reports no statistically significant difference from the anti-TNF Emax.

    emax_antiil6 <- fixed(5.604)
    label("Maximal ACR 50 treatment effect of the anti-IL-6 class, i.e. tocilizumab (logit units)")
    # Digitized from the Mandema 2011 Figure 2 tocilizumab panel. This is the
    # one class Emax the source calls out as different: "The theoretical Emax
    # of tocilizumab (anti-IL-6) was significantly greater than the Emax of
    # anti-TNFs." The recovered value also reproduces the two quantitative
    # claims made about it -- that 8 mg/kg q4w gives an effect "similar to
    # those of etanercept and certolizumab", and that doubling the tocilizumab
    # dose buys much more response than doubling an anti-TNF dose.

    # ------------------------------------------------------------------------
    # ED50 by DRUG, log-transformed. "Dose is the total daily/weekly/monthly
    # dose normalized to the standard regimen for each drug, and ED50 is the
    # dose required to achieve 50% of Emax. A different ED50 was estimated for
    # each drug." Units follow the drug and match its CONMED_<drug>_DOSE column.
    #
    # Every ED50 below is the value that makes the class Emax above reproduce
    # that drug's published Table 2 effect at its suggested starting dose
    # exactly, so all nine Table 2 rows are regenerated by construction. For
    # the five anti-TNFs this is a genuine second source of information (the
    # Emax came from the Discussion ratios, not from Table 2 alone) and the
    # resulting dose/ED50 values reproduce the published ratios; for the four
    # single-drug classes it converts the digitized curvature into a potency.
    # ------------------------------------------------------------------------
    led50_adalimumab <- fixed(log(17.80647))
    label("Log adalimumab ED50 (mg q2w); back-transform ED50 = 17.81 mg q2w")
    # 40 mg q2w / 2.2 published dose/ED50 ratio; recovered 17.81, ratio 2.25.

    led50_certolizumab <- fixed(log(11.77989))
    label("Log certolizumab pegol ED50 (mg q2w); back-transform ED50 = 11.78 mg q2w")
    # 200 mg q2w at a recovered ratio of 17.0, consistent with the source's
    # only statement about it, "certolizumab at >10 times the ED50".

    led50_etanercept <- fixed(log(3.83580))
    label("Log etanercept ED50 (mg twice weekly); back-transform ED50 = 3.836 mg biw")
    # 25 mg biw / 6.8 published dose/ED50 ratio; recovered 3.836, ratio 6.52.

    led50_golimumab <- fixed(log(45.70067))
    label("Log golimumab ED50 (mg q4w); back-transform ED50 = 45.70 mg q4w")
    # 50 mg q4w / 1.1 published dose/ED50 ratio; recovered 45.70, ratio 1.09.

    led50_infliximab <- fixed(log(1.42148))
    label("Log infliximab ED50 (mg/kg on the wk 0, 2, 6 then q8w regimen); back-transform ED50 = 1.421 mg/kg")
    # 3 mg/kg / 2.1 published dose/ED50 ratio; recovered 1.421, ratio 2.11.

    led50_anakinra <- fixed(log(11.72318))
    label("Log anakinra ED50 (mg/day); back-transform ED50 = 11.72 mg/day")
    # Set so that emax_antiil1 reproduces the Table 2 anakinra effect at
    # 100 mg/day; implied dose/ED50 = 8.53 at the suggested starting dose.

    led50_abatacept <- fixed(log(2.53154))
    label("Log abatacept ED50 (mg/kg q4w); back-transform ED50 = 2.532 mg/kg q4w")
    # Set so that emax_cd28 reproduces the Table 2 abatacept effect at
    # 10 mg/kg q4w; implied dose/ED50 = 3.95.

    led50_rituximab <- fixed(log(110.97090))
    label("Log rituximab ED50 (mg at weeks 0 and 2); back-transform ED50 = 111.0 mg")
    # Set so that emax_cd20 reproduces the Table 2 rituximab effect at
    # 1,000 mg; implied dose/ED50 = 9.01.

    led50_tocilizumab <- fixed(log(12.81558))
    label("Log tocilizumab ED50 (mg/kg q4w); back-transform ED50 = 12.82 mg/kg q4w")
    # Set so that emax_antiil6 reproduces the Table 2 tocilizumab effect at
    # 4 mg/kg q4w; implied dose/ED50 = 0.31, i.e. the whole studied 2-8 mg/kg
    # range sits BELOW the ED50. This is the quantitative content of the
    # source's observation that the tocilizumab dose-response is still rising
    # steeply at the top of its clinical range.

    # ------------------------------------------------------------------------
    # Methotrexate given as a RANDOMIZED treatment: a dose-independent mean
    # effect. "The treatment effect of MTX when given as a randomized
    # treatment was estimated as a simple mean effect because MTX was titrated
    # to a similar dose range in all the trials and no dose-response data were
    # available." It uses the same end-point scaling as the biologics: the
    # recovered ACR 20 and ACR 70 methotrexate effects are 0.858 and 1.113
    # times the ACR 50 effect, matching scale_acr20 and scale_acr70 above to
    # three decimal places without being constrained to.
    # ------------------------------------------------------------------------
    e_mtx <- fixed(1.585319)
    label("ACR 50 treatment effect of randomized methotrexate monotherapy (logit units)")
    # Recovered together with e0_acr50_naive from Mandema 2011 Table 3 and the
    # printed 32/24/12 percentage-point methotrexate statement (p. 833):
    # logit(0.092671 + 0.24) - logit(0.092671) = 1.585319.

    # ------------------------------------------------------------------------
    # Initial combination therapy of a biologic with methotrexate:
    #   g = g_biologic + g_MTX + gamma * g_biologic * g_MTX
    # Only anti-TNF plus methotrexate combinations were evaluated (adalimumab,
    # etanercept, golimumab, infliximab).
    # ------------------------------------------------------------------------
    gamma_mtx <- fixed(-0.32)
    label("Interaction coefficient for initial combination therapy of a biologic with methotrexate (unitless)")
    # Mandema 2011 p. 831: "The interaction coefficient gamma for anti-TNF and
    # MTX initial combination therapy was -0.32 (-0.37 to -0.27), significantly
    # different from zero, showing that the effect of initial combination
    # therapy is less than the sum of the effects of each component."

    # ------------------------------------------------------------------------
    # East Asian trials. Acts on the Emax of the dose-response through the
    # source's covariate relationship Emax,i = Emax,class * (1 + theta_c *
    # (X_ij - Xbar) + eta_Emax,i) with a binary X, i.e. it scales the whole
    # biologic treatment effect.
    # ------------------------------------------------------------------------
    beta_eastasia <- fixed(0.33)
    label("Fractional increase in the biologic treatment effect for trials carried out in East Asia (unitless)")
    # Mandema 2011 p. 831: "In trials carried out in Asia, the reported
    # treatment effects were significantly (33% (95% CI, 12 to 54%)) greater
    # than those reported in trials carried out elsewhere in the world."

    # ========================================================================
    # No random effects and no residual error model.
    #
    # Random effects: between-trial random effects on the dose-response
    # parameters were tested and rejected -- "This is confirmed by the absence
    # of significant between-trial heterogeneity with respect to Emax and
    # ED50" -- so the final model carries none on the treatment effect. The
    # one random effect the source does retain is h_i,k on the ACR 20 and
    # ACR 70 placebo responses (variance omega_k^2), whose variances are NOT
    # reported anywhere in the paper or its Supplementary Data. Per the
    # unreported-variance convention these are omitted rather than invented;
    # the observed spread they describe is recorded in
    # covariatesDataExcluded$PLACEBO_RESPONSE_MAGNITUDE.
    #
    # Residual error: the observations were per-arm responder COUNTS with a
    # binomial likelihood (N_event,kij ~ binomial(P(event)kij, N_kij)); there
    # is no additive or proportional residual SD to encode. The outputs below
    # are therefore deterministic per-arm probabilities. A user who wants
    # simulated counts should draw rbinom(n = 1, size = N_arm, prob = p_<end
    # point>) downstream, which is exactly what the source's likelihood
    # assumes.
    # ========================================================================
  })

  model({
    # ---- Per-drug ED50, back-transformed -----------------------------------
    ed50_adalimumab   <- exp(led50_adalimumab)
    ed50_certolizumab <- exp(led50_certolizumab)
    ed50_etanercept   <- exp(led50_etanercept)
    ed50_golimumab    <- exp(led50_golimumab)
    ed50_infliximab   <- exp(led50_infliximab)
    ed50_anakinra     <- exp(led50_anakinra)
    ed50_abatacept    <- exp(led50_abatacept)
    ed50_rituximab    <- exp(led50_rituximab)
    ed50_tocilizumab  <- exp(led50_tocilizumab)

    # ---- Saturating dose fraction Dose/(Dose + ED50), per drug -------------
    # Every trial arm in the source received at most ONE biologic, so at most
    # one term of each class sum is non-zero. A placebo or methotrexate-only
    # arm sets every CONMED_<drug>_DOSE column to 0, which makes every fraction zero
    # and collapses the biologic effect to nothing. Coding two biologics
    # simultaneously is outside the source's calibration (only one of the 50
    # trials compared two biologics head to head, and never in the same arm)
    # and would make the model ADD their effects.
    f_adalimumab   <- CONMED_ADALIMUMAB_DOSE    / (CONMED_ADALIMUMAB_DOSE    + ed50_adalimumab)
    f_certolizumab <- CONMED_CERTOLIZUMAB_DOSE  / (CONMED_CERTOLIZUMAB_DOSE  + ed50_certolizumab)
    f_etanercept   <- CONMED_ETANERCEPT_DOSE    / (CONMED_ETANERCEPT_DOSE    + ed50_etanercept)
    f_golimumab    <- CONMED_GOLIMUMAB_DOSE     / (CONMED_GOLIMUMAB_DOSE     + ed50_golimumab)
    f_infliximab   <- CONMED_INFLIXIMAB_DOSE  / (CONMED_INFLIXIMAB_DOSE  + ed50_infliximab)
    f_anakinra     <- CONMED_ANAKINRA_DOSE      / (CONMED_ANAKINRA_DOSE      + ed50_anakinra)
    f_abatacept    <- CONMED_ABATACEPT_DOSE   / (CONMED_ABATACEPT_DOSE   + ed50_abatacept)
    f_rituximab    <- CONMED_RITUXIMAB_DOSE     / (CONMED_RITUXIMAB_DOSE     + ed50_rituximab)
    f_tocilizumab  <- CONMED_TOCILIZUMAB_DOSE / (CONMED_TOCILIZUMAB_DOSE + ed50_tocilizumab)

    # ---- Emax by mechanism of action ---------------------------------------
    # The Emax is a CLASS property and the ED50 is a DRUG property, which is
    # the source's central structural finding: all five anti-TNFs share one
    # dose-response shape and differ only in where their clinical dose range
    # sits on it.
    g_antitnf <- emax_antitnf *
      (f_adalimumab + f_certolizumab + f_etanercept + f_golimumab + f_infliximab)
    g_antiil1 <- emax_antiil1 * f_anakinra
    g_cd28    <- emax_cd28    * f_abatacept
    g_cd20    <- emax_cd20    * f_rituximab
    g_antiil6 <- emax_antiil6 * f_tocilizumab

    # ---- Biologic treatment effect on the ACR 50 logit scale ---------------
    # The East Asian effect multiplies the class Emax and therefore the whole
    # biologic effect. It is NOT applied to the randomized-methotrexate mean
    # effect below, which has no Emax to modify.
    g_bio <- (g_antitnf + g_antiil1 + g_cd28 + g_cd20 + g_antiil6) *
      (1 + beta_eastasia * REGION_EASTASIA)

    # ---- Randomized methotrexate: a dose-independent mean effect -----------
    # Only presence matters; the source estimated no methotrexate dose
    # response. Background (non-randomized) methotrexate is absorbed into the
    # placebo response and contributes nothing here.
    mtx_randomized <- (CONMED_MTX_DOSE > 0)
    g_mtx          <- e_mtx * mtx_randomized

    # ---- End-point-specific component effects ------------------------------
    # The end-point scaling acts on Emax, so it applies to each COMPONENT of
    # the treatment effect before they are combined.
    g_bio20 <- scale_acr20 * g_bio
    g_bio70 <- scale_acr70 * g_bio
    g_mtx20 <- scale_acr20 * g_mtx
    g_mtx70 <- scale_acr70 * g_mtx

    # ---- Initial combination therapy ---------------------------------------
    # g(x)_k = g(x)_biologic,k + g(x)_MTX,k + gamma * g(x)_biologic,k *
    # g(x)_MTX,k. The source writes the interaction on the end-point-specific
    # functions g(x)_k, so the product term is formed AFTER the end-point
    # scaling rather than before it. The distinction only matters for
    # combination arms (the product vanishes when either component is absent)
    # but it matters a lot there: forming the product on the ACR 50 scale and
    # scaling afterwards misses the published Table 4 ACR 20 and ACR 70
    # columns by up to 2.1 percentage points, while this form reproduces them
    # to within 0.4. See the validation vignette.
    g_acr20 <- g_bio20 + g_mtx20 + gamma_mtx * g_bio20 * g_mtx20
    g_acr50 <- g_bio   + g_mtx   + gamma_mtx * g_bio   * g_mtx
    g_acr70 <- g_bio70 + g_mtx70 + gamma_mtx * g_bio70 * g_mtx70

    # ---- Log odds ratios versus the arm's own trial control -----------------
    # These are the source's stated primary response variable and are the
    # quantities identified without reference to any trial placebo intercept,
    # so they are the outputs to use when comparing treatments.
    lor_acr20 <- g_acr20
    lor_acr50 <- g_acr50
    lor_acr70 <- g_acr70

    # ---- Per-arm responder probabilities, methotrexate-inadequate responders
    # (the Table 2 population: placebo means background methotrexate).
    e0_20 <- e0_acr50 + theta_acr20
    e0_70 <- e0_acr50 + theta_acr70
    p_acr20 <- 1 / (1 + exp(-(e0_20    + g_acr20)))
    p_acr50 <- 1 / (1 + exp(-(e0_acr50 + g_acr50)))
    p_acr70 <- 1 / (1 + exp(-(e0_70    + g_acr70)))

    # ---- Per-arm responder probabilities, methotrexate-naive population -----
    # (the Table 3 and Table 4 population). The treatment effect is the SAME:
    # "There was no statistically significant difference in the odds ratio
    # (relative to placebo) for MTX-naive patients vs. patients who had shown
    # inadequate response to prior MTX treatment." Only the placebo response
    # differs, and it is slightly larger.
    e0_20_naive <- e0_acr50_naive + theta_acr20_naive
    e0_70_naive <- e0_acr50_naive + theta_acr70_naive
    p_acr20_naive <- 1 / (1 + exp(-(e0_20_naive    + g_acr20)))
    p_acr50_naive <- 1 / (1 + exp(-(e0_acr50_naive + g_acr50)))
    p_acr70_naive <- 1 / (1 + exp(-(e0_70_naive    + g_acr70)))

    # ---- Derived reporting quantity ----------------------------------------
    # Multiples of the arm's own ED50: the x-axis on which the source's
    # Figure 4 collapses all five anti-TNFs onto a single curve. With one
    # active drug per arm the total fraction f is that drug's Dose/(Dose +
    # ED50), so f/(1 - f) is its Dose/ED50; a placebo arm gives 0.
    f_total <- f_adalimumab + f_certolizumab + f_etanercept +
      f_golimumab + f_infliximab + f_anakinra +
      f_abatacept + f_rituximab + f_tocilizumab
    dose_over_ed50 <- f_total / (1 - f_total + 1e-12)
  })
}
