Mandema_2005_gemcabene_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis dose-response model for the percentage change in low-density lipoprotein cholesterol (LDL-C) from pretreatment baseline for five statins (atorvastatin, rosuvastatin, simvastatin, lovastatin, pravastatin), the cholesterol absorption inhibitor ezetimibe, the investigational lipid-altering agent gemcabene, and statin-plus-nonstatin combinations. Fit to study-arm summary means from 25 randomized controlled trials (9,886 patients) by nonlinear mixed-effects regression in S-PLUS 6.1. Each drug follows a sigmoid Emax dose-response in daily dose; the five statins share a common Emax (-78.7 percent) and a common shallow Hill coefficient (0.451) and differ only in ED50 (potency), while ezetimibe and gemcabene each carry their own Emax, ED50 and Hill coefficient. Combination therapy uses the multiplicative interaction term 0.01 * g * Estatin * Enonstatin, where g = 1 is pharmacological independence: ezetimibe was found independent (g fixed to 1, giving added LDL-C lowering across the whole statin dose range) whereas gemcabene was significantly LESS than independent (g = 1.69), so gemcabene adds almost nothing on top of a high statin dose. Per-arm daily dose enters through one CONMED_<drug>_DOSE column per drug (0 outside that drug's arm; all zero = placebo) and European trials carry an extra 4 percent LDL-C reduction on the placebo intercept. The trial-specific random-effect variance was not significantly different from zero and the between-subject residual variance is not reported numerically, so the model has no eta and its residual SD is fixed to zero. Suitable simulation scope is the study-arm mean percentage change in LDL-C, NOT individual-subject responses; there is no PK layer and no time course."

  reference <- paste(
    "Mandema JW, Hermann D, Wang W, Sheiner T, Milad M, Bakker-Arkema R, Hartman D.",
    "Model-based development of gemcabene, a new lipid-altering agent.",
    "AAPS J. 2005 Oct 7;7(3):E513-E522.",
    "doi:10.1208/aapsj070352.",
    "Model structure is Equations 1 and 2 (Materials and Methods, Statistical Analysis);",
    "all parameter estimates are in Table 2.",
    sep = " "
  )

  vignette <- "Mandema_2005_gemcabene_mbma"

  # Algebraic MBMA dose-response model: no rxode2 dose events are consumed (the
  # per-arm daily dose enters through the CONMED_<drug>_DOSE covariate columns)
  # and the output is a percentage change from baseline rather than a drug
  # concentration. The `units` entries follow the placeholder convention already
  # used by Vargo_2014_statins_ezetimibe_mbma and Mandema_2011_anticoagulants_mbma
  # so that checkModelConventions() sees a parseable dosing / concentration pair.
  units <- list(
    time          = "day (placeholder; the model is a time-independent steady-state per-arm dose-response, not a time course. Mandema 2005 required at least 4 weeks of treatment for trial inclusion and found no significant effect of treatment duration on any dose-response parameter)",
    dosing        = "mg/day (per-arm TOTAL DAILY dose of each drug, supplied through the CONMED_<drug>_DOSE covariate columns, NOT as rxode2 dose events)",
    concentration = "%/baseline (percentage change in LDL-C from the pretreatment baseline; negative values are reductions, e.g. Cc = -31.9 is a 31.9 percent LDL-C reduction. Output Cc is NOT a drug concentration; the slash satisfies checkModelConventions parsing)"
  )

  covariateData <- list(
    CONMED_ATORVASTATIN_DOSE = list(
      description        = "Per-arm total daily atorvastatin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an atorvastatin arm. Atorvastatin is the most frequent statin in the trial set and is the only statin studied in combination with gemcabene (Mandema 2005 Table 1, the unpublished 2003 phase IIA hypercholesterolemia trial). Dose strengths across Table 1 span 2.5-80 mg/day. Member of the CONMED_<drug>_DOSE family; drives the shared-Emax statin dose-response term through the atorvastatin ED50 of 13.1 mg/day (Table 2).",
      source_name        = "Drug, Dose 'A<dose>' (Mandema 2005 Table 1); ED50, Atorvastatin (Table 2)"
    ),
    CONMED_RSV_DOSE = list(
      description        = "Per-arm total daily rosuvastatin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a rosuvastatin arm. Dose strengths across Mandema 2005 Table 1 span 1-80 mg/day. Rosuvastatin is the most potent statin in the analysis (ED50 4.35 mg/day, Table 2); the paper's Discussion gives its potency relative to atorvastatin as 0.33.",
      source_name        = "Drug, Dose 'R<dose>' (Mandema 2005 Table 1); ED50, Rosuvastatin (Table 2)"
    ),
    CONMED_SMV_DOSE = list(
      description        = "Per-arm total daily simvastatin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a simvastatin arm. Dose strengths across Mandema 2005 Table 1 span 10-80 mg/day. ED50 30.5 mg/day (Table 2); potency relative to atorvastatin 2.3 (Discussion).",
      source_name        = "Drug, Dose 'S<dose>' (Mandema 2005 Table 1); ED50, Simvastatin (Table 2)"
    ),
    CONMED_LOV_DOSE = list(
      description        = "Per-arm total daily lovastatin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a lovastatin arm. Dose strengths across Mandema 2005 Table 1 span 10-80 mg/day. ED50 82.8 mg/day (Table 2); potency relative to atorvastatin 6.3 (Discussion). Unlike Vargo 2014, Mandema 2005 did not estimate a separate ED50 for twice-daily or extended-release lovastatin, so there is no FORM_LOV_BID_XR analogue in this model.",
      source_name        = "Drug, Dose 'L<dose>' (Mandema 2005 Table 1); ED50, Lovastatin (Table 2)"
    ),
    CONMED_PRV_DOSE = list(
      description        = "Per-arm total daily pravastatin dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a pravastatin arm. Dose strengths across Mandema 2005 Table 1 span 10-40 mg/day. Pravastatin is the least potent statin in the analysis (ED50 97.3 mg/day, Table 2); potency relative to atorvastatin 7.4 (Discussion). Note the narrow studied dose range relative to the ED50, which is why the ED50 confidence interval is the widest of the five statins (42.4-223).",
      source_name        = "Drug, Dose 'P<dose>' (Mandema 2005 Table 1); ED50, Pravastatin (Table 2)"
    ),
    CONMED_EZT_DOSE = list(
      description        = "Per-arm total daily ezetimibe dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside an ezetimibe arm. Dose strengths across Mandema 2005 Table 1 span 0.25-40 mg/day. Ezetimibe is one of the two 'nonstatin' drugs in Equation 1 and carries its own Emax (-19.6 percent), ED50 (0.302 mg/day) and interaction coefficient. Because the ED50 is only about 0.3 mg, the marketed 10 mg dose sits essentially at Emax (paper Discussion: 19.1 percent reduction versus a 19.6 percent maximum).",
      source_name        = "Drug, Dose 'E<dose>' (Mandema 2005 Table 1); ED50,Ezetimibe (Table 2)"
    ),
    CONMED_GEMCABENE_DOSE = list(
      description        = "Per-arm total daily gemcabene dose.",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "0 outside a gemcabene arm. Gemcabene (development code CI-1027) is the investigational lipid-altering agent whose development this paper supported; it was never marketed and the paper documents the decision to stop development. Dose strengths across the four Pfizer-sponsored trials in Mandema 2005 Table 1 span 50-900 mg/day. New member of the CONMED_<drug>_DOSE family, spelled out in full per that family's naming rule (no established 3-4 letter abbreviation, and spelling out avoids any collision). Gemcabene is the second 'nonstatin' drug in Equation 1 and carries its own Emax (-34.8 percent), ED50 (314 mg/day), Hill coefficient (2.27) and, critically, its own interaction coefficient g = 1.69.",
      source_name        = "Drug, Dose 'G<dose>' (Mandema 2005 Table 1); ED50,gemcabene (Table 2)"
    ),
    REGION_EUROPE = list(
      description        = "Indicator that the trial was conducted in Europe rather than North America.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = North American (United States) trial, the reference location whose placebo intercept is the tabulated E0 of 0.802 percent.",
      notes              = "Study-arm-level (really trial-level) covariate, not a subject-level covariate. Mandema 2005 Results, LDL-C section: 'A small but statistically significant difference between North American and European studies was found for the overall placebo effect. An additional 4% reduction in LDL-C was consistently observed in Europeans trials.' The effect acts on the intercept E0 only, so it shifts every arm of a European trial (placebo and active alike) by the same -4 percent and therefore CANCELS from any within-trial active-minus-control contrast such as Mandema 2005 Table 3. The 4 percent figure is reported only in that prose sentence -- it has no row in Table 2 and no confidence interval -- so it is recorded to one significant figure exactly as printed. Location (United States or Europe) is listed in the Methods among the explanatory variables collected for each trial; the per-trial assignment itself is not tabulated in Table 1.",
      source_name        = "location (United States or Europe) (Mandema 2005 Methods, Studies Included; effect size in Results, LDL-C section)"
    )
  )

  covariatesDataExcluded <- list(
    LDLC = list(
      description        = "Per-arm mean baseline (pretreatment) LDL-C concentration. Collected and screened but NOT retained as a dose-response covariate in the final model.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tabulated per trial in Mandema 2005 Table 1 (112-217 mg/dL across the 24 trials that report it) and listed in the Methods among the collected explanatory variables. Results, LDL-C section: 'No statistically significant effect was found for the other explanatory variables on E0, Emax, ED50, and n for the statins.' No point estimate is reported, so no effect parameter appears in ini(). Baseline LDL-C DID matter at a different level: the three gemcabene healthy-subject trials had much lower baseline LDL-C (112-120 mg/dL) than the hypercholesterolemia trials, and the paper handled this by an upstream per-trial ANOVA (factors dose, trial, baseline LDL-C, triglycerides) that found 'a small impact of baseline LDL-C values < 100 mg/dL on the mean percentage of LDL-C reduction' and NO treatment-by-baseline interaction. The ANOVA least-squares means were then used as the summary data fed to this dose-response model, so that adjustment is baked into the input data rather than into the model equations. Contrast Vargo 2014, which DID retain a log(LDL.base/180) effect on statin Emax.",
      source_name        = "Baseline LDL-C (mg/dL) (Mandema 2005 Table 1; Methods explanatory-variable list)"
    ),
    TRIG = list(
      description        = "Per-arm mean baseline (pretreatment) triglyceride concentration. Used only in the upstream gemcabene ANOVA; NOT a covariate in the dose-response model.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Mandema 2005 Results, Data section: the ANOVA over the three gemcabene healthy-subject trials carried 'factors for dose, trial, baseline LDL-C, and triglycerides'. Triglycerides do not appear in the dose-response model and no coefficient is reported. One of the four gemcabene trials enrolled 'subjects with low HDL and normal or elevated TG' (Table 1), which is why triglycerides were screened at all. Per-arm triglyceride values are not tabulated.",
      source_name        = "triglycerides (Mandema 2005 Results, Data section ANOVA factor list)"
    ),
    TRT_DURATION = list(
      description        = "Trial treatment duration. Collected and screened but NOT retained as a dose-response covariate in the final model.",
      units              = "week",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tabulated per trial in Mandema 2005 Table 1 (4-16 weeks) and listed in the Methods among the collected explanatory variables; trials shorter than 4 weeks were excluded at screening. Not retained ('No statistically significant effect was found for the other explanatory variables on E0, Emax, ED50, and n'), and no point estimate is reported. The absence of a duration effect over 4-16 weeks is what licenses treating the model as a steady-state, time-independent dose-response.",
      source_name        = "Duration (weeks) (Mandema 2005 Table 1; Methods explanatory-variable list)"
    ),
    TRIAL_START_YEAR = list(
      description        = "Year the trial was published. Collected and screened but NOT retained as a dose-response covariate in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Mandema 2005 Methods: 'Variables that were collected include treatment duration, baseline LDL-C, location (United States or Europe), and the year the trial was published.' Tabulated per trial in Table 1 (1995-2003). Not retained and no point estimate reported. Named TRIAL_START_YEAR to reuse the existing canonical from Mandema_2011_anticoagulants_mbma even though this paper collected year of PUBLICATION rather than year of trial start.",
      source_name        = "Year (Mandema 2005 Table 1; Methods explanatory-variable list, 'the year the trial was published')"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 9886L,
    n_studies        = 25L,
    disease_state    = "Predominantly adults with hypercholesterolemia / elevated LDL-C (21 trials). The gemcabene programme additionally contributed three non-hypercholesterolemic cohorts: healthy volunteers, nondiabetic healthy obese subjects, and subjects with low HDL-C and normal or elevated triglycerides. Baseline LDL-C was much lower in those three healthy-subject trials (112-120 mg/dL) than in the hypercholesterolemia trials (165-217 mg/dL), which is why their arm means were first adjusted by a per-trial ANOVA before entering the dose-response analysis.",
    dose_range       = "atorvastatin 2.5-80 mg/day; rosuvastatin 1-80 mg/day; simvastatin 10-80 mg/day; lovastatin 10-80 mg/day; pravastatin 10-40 mg/day; ezetimibe 0.25-40 mg/day; gemcabene 50-900 mg/day (Mandema 2005 Table 1). All regimens once daily.",
    baseline_ldlc    = "112-217 mg/dL across trials (Mandema 2005 Table 1; one ezetimibe trial does not report a baseline).",
    regions          = "United States and Europe. The Methods collected trial location as an explanatory variable and a 4 percent larger placebo LDL-C reduction was found in European trials; the per-trial location assignment is not tabulated.",
    n_trials_statin_ezetimibe = 21L,
    n_trials_gemcabene        = 4L,
    data_sources     = "Medline search to a 31 May 2003 cut-off for randomized controlled trials of at least 4 weeks in hypercholesterolemia, plus FDA summary-basis-of-approval documents, plus Pfizer internal study reports. All 21 statin / ezetimibe trials were randomized, double-blind, multiple-dose once-daily parallel-group studies; 13 were placebo controlled and most were multicentre. Three of the four gemcabene trials were unpublished at the time of writing (marked 'n.p.' in Table 1); one has since been published.",
    notes            = "Summary-level MBMA: each modelled observation is a study-arm mean percentage change in LDL-C from pretreatment baseline, not an individual patient value. The model is intended for simulating study-arm-mean percentage change in LDL-C and is NOT suitable for individual-subject simulation. There is no PK layer and no time course. Only the LDL-C endpoint is parameterised here: the paper's Introduction and Abstract also describe models for HDL-C, persistent alanine aminotransferase elevation, myalgia, headache and coronary-artery-disease risk reduction, but states that 'the LDL-C effect ... is the main focus of this article' and reports parameter estimates for the LDL-C dose-response only (Table 2). No parameter values are published for those other endpoints anywhere in the paper, so they cannot be extracted. Gemcabene development was stopped on the strength of this analysis."
  )

  ini({
    # ========================================================================
    # Mandema 2005 Equation 1 (Materials and Methods, Statistical Analysis):
    #
    #   Y = E0 + Estatin + Enonstatin
    #         + 0.01 * g * Estatin * Enonstatin + eta + eps
    #
    # and Equation 2, the sigmoid Emax dose-response used for each drug:
    #
    #   Edrug = Dose^n * Emax / (Dose^n + ED50^n)
    #
    # Y is the percentage change in LDL-C from pretreatment baseline, so all
    # E terms below are on the PERCENT scale (not fractional) and a negative
    # value is an LDL-C reduction. The 0.01 factor in Equation 1 is exactly
    # what converts the product of two percentages back to a percentage: with
    # g = 1 the equation reproduces pharmacological independence,
    # 1 - (1 - f_statin)(1 - f_nonstatin), on the fractional scale (paper
    # Discussion). g > 1 is LESS than independent.
    #
    # All values below: Mandema 2005 Table 2. Note that the published PDF
    # renders the minus sign of every negative Table 2 entry as a leading "2"
    # (e.g. Emax,statin prints as "278.7" for -78.7); the signs used here are
    # confirmed by reproducing every cell of the paper's own Table 3 from
    # Table 2 (see the vignette's Table 3 replication).
    # ========================================================================

    e0 <- 0.802
    label("Placebo effect: percentage change in LDL-C on placebo, North American trials (%)")  # Mandema 2005 Table 2, E0 = 0.802 [0.0598 to 1.54]

    e_region_europe_e0 <- -4
    label("Additional percentage-point change in LDL-C on the placebo intercept for European trials (%)")  # Mandema 2005 Results, LDL-C section: "An additional 4% reduction in LDL-C was consistently observed in Europeans trials." Prose-only: no Table 2 row and no confidence interval; recorded to the single significant figure printed.

    # ------------------------------------------------------------------------
    # Statins. A single Emax and a single Hill coefficient are SHARED across
    # all five statins; only the ED50 (potency) differs. Paper Results: "No
    # statistically significant difference in Emax or Hill coefficient (n) was
    # found between the statins ... This suggests that all of the statins share
    # a similar shape of the dose-response relationship with a similar maximal
    # effect of about 79% reduction in LDL-C over placebo."
    # ------------------------------------------------------------------------
    emax_statin <- -78.7
    label("Maximal LDL-C reduction shared by all five statins (%, signed; at infinite dose)")  # Mandema 2005 Table 2, Emax, statin (%) = -78.7 [-90.7 to -66.7]

    led50_atorvastatin <- log(13.1)
    label("Atorvastatin ED50 (mg/day)")  # Mandema 2005 Table 2, ED50, Atorvastatin (mg) = 13.1 [6.57 to 26.2]

    led50_rsv <- log(4.35)
    label("Rosuvastatin ED50 (mg/day)")  # Mandema 2005 Table 2, ED50, Rosuvastatin (mg) = 4.35 [2.19 to 8.62]

    led50_smv <- log(30.5)
    label("Simvastatin ED50 (mg/day)")  # Mandema 2005 Table 2, ED50, Simvastatin (mg) = 30.5 [15 to 62.1]

    led50_lov <- log(82.8)
    label("Lovastatin ED50 (mg/day)")  # Mandema 2005 Table 2, ED50, Lovastatin (mg) = 82.8 [37.1 to 185]

    led50_prv <- log(97.3)
    label("Pravastatin ED50 (mg/day)")  # Mandema 2005 Table 2, ED50, Pravastatin (mg) = 97.3 [42.4 to 223]

    ln_statin <- log(0.451)
    label("Hill coefficient shared by all five statins (unitless)")  # Mandema 2005 Table 2, n statin = 0.451 [0.366 to 0.557]

    # ------------------------------------------------------------------------
    # Ezetimibe. A SIMPLE (non-sigmoid) Emax model: "A Hill coefficient for
    # ezetimibe could not be estimated and was fixed to 1." The interaction
    # coefficient was likewise fixed: "The interaction coefficient (g) was not
    # statistically significantly different from 1 and was, therefore, fixed to
    # this pharmacological value." Both carry no confidence interval in Table 2,
    # which is the tabulated signal that they are fixed rather than estimated.
    # ------------------------------------------------------------------------
    emax_ezt <- -19.6
    label("Ezetimibe maximal LDL-C reduction (%, signed)")  # Mandema 2005 Table 2, Emax,Ezetimibe (%) = -19.6 [-20.6 to -18.6]

    led50_ezt <- log(0.302)
    label("Ezetimibe ED50 (mg/day)")  # Mandema 2005 Table 2, ED50,Ezetimibe (mg) = 0.302 [0.151 to 0.604]

    ln_ezt <- fixed(log(1))
    label("Ezetimibe Hill coefficient (unitless; held at 1 per the paper's Methods)")  # Mandema 2005 Table 2, n Ezetimibe = 1 (no CI); Results: "A Hill coefficient for ezetimibe could not be estimated and was fixed to 1."

    gamma_ezt <- fixed(1)
    label("Statin x ezetimibe interaction coefficient (unitless; 1 = pharmacologically independent)")  # Mandema 2005 Table 2, g Ezetimibe = 1 (no CI); Results: "The interaction coefficient (g) was not statistically significantly different from 1 and was, therefore, fixed to this pharmacological value." When estimated it was 0.87 [0.70 to 1.04] (Discussion).

    # ------------------------------------------------------------------------
    # Gemcabene. A full sigmoid Emax, and the one parameter that decided the
    # programme: g = 1.69 is significantly greater than 1 (p < 0.001), i.e. the
    # interaction with statins is LESS than independent, so gemcabene adds very
    # little on top of a high statin dose.
    # ------------------------------------------------------------------------
    emax_gem <- -34.8
    label("Gemcabene maximal LDL-C reduction (%, signed)")  # Mandema 2005 Table 2, Emax,gemcabene (%) = -34.8 [-45 to -24.6]

    led50_gem <- log(314)
    label("Gemcabene ED50 (mg/day)")  # Mandema 2005 Table 2, ED50,gemcabene (mg) = 314 [220 to 448]

    ln_gem <- log(2.27)
    label("Gemcabene Hill coefficient (unitless)")  # Mandema 2005 Table 2, n gemcabene = 2.27 [1.19 to 4.34]

    gamma_gem <- 1.69
    label("Statin x gemcabene interaction coefficient (unitless; >1 = less than independent)")  # Mandema 2005 Table 2, g gemcabene = 1.69 [1.49 to 1.88]; Results also quote 1.69 +/- 0.10 (mean +/- SE) and p < 0.001 versus 1.

    # ========================================================================
    # Random effects and residual error.
    #
    # Trial-specific random effect (eta, variance nu^2 in Equation 1): NOT
    # encoded. Paper Results, LDL-C section: "The variance of the trial-specific
    # random effect was found to be not statistically significantly different
    # from zero suggesting homogeneity between studies in the magnitude of LDL-C
    # reduction observed for a specific treatment strategy." No numeric estimate
    # is printed. This follows the same treatment as the sibling MBMA models
    # Vargo_2014_statins_ezetimibe_mbma and Mandema_2011_anticoagulants_mbma,
    # whose between-trial variances were likewise driven to zero; encoding a
    # zero-variance eta instead would give a singular OMEGA that rxode2's
    # Cholesky sampler cannot decompose.
    #
    # Between-subject residual (eps, variance sigma^2 in Equation 1): the
    # paper declares it but reports NO numeric value for it, anywhere -- not in
    # Table 2, not in the Results text, and there is no supplement. Rather than
    # invent a variance, addSd is fixed to zero so the endpoint is declared (a
    # downstream user can unfix and estimate it) while the model predicts the
    # typical study-arm mean exactly. This is recorded in the vignette's
    # Assumptions and deviations section.
    # ========================================================================
    addSd <- fixed(0)
    label("Additive residual SD on the percentage change in LDL-C (%); 0 because Mandema 2005 reports no value for sigma")  # Mandema 2005 Equation 1 declares eps ~ N(0, sigma^2); no numeric estimate of sigma is published.
  })

  model({
    # ---- Back-transform the log-scale shape and potency parameters ---------
    n_statin <- exp(ln_statin)
    n_ezt    <- exp(ln_ezt)
    n_gem    <- exp(ln_gem)

    ed50_atorvastatin <- exp(led50_atorvastatin)
    ed50_rsv          <- exp(led50_rsv)
    ed50_smv          <- exp(led50_smv)
    ed50_lov          <- exp(led50_lov)
    ed50_prv          <- exp(led50_prv)
    ed50_ezt          <- exp(led50_ezt)
    ed50_gem          <- exp(led50_gem)

    # ---- Statin dose-response (Equation 2) ---------------------------------
    # Equation 2 is Dose^n * Emax / (Dose^n + ED50^n), which is identically
    # Emax * r^n / (r^n + 1) with r = Dose / ED50. Because the five statins
    # share Emax and n and differ only in ED50, summing the potency-normalised
    # doses r selects whichever statin the arm actually received: every trial
    # arm in Mandema 2005 Table 1 contains at most ONE statin, so at most one
    # term is non-zero and this is exact. A placebo arm sets every
    # CONMED_<drug>_DOSE to 0, making r_statin 0 and collapsing the term to 0.
    # Coding two statins simultaneously is outside the paper's calibration and
    # would add their potency-normalised doses (the pharmacologically sensible
    # generalisation, but not something the paper fitted).
    r_statin <- CONMED_ATORVASTATIN_DOSE / ed50_atorvastatin +
      CONMED_RSV_DOSE / ed50_rsv +
      CONMED_SMV_DOSE / ed50_smv +
      CONMED_LOV_DOSE / ed50_lov +
      CONMED_PRV_DOSE / ed50_prv

    u_statin <- r_statin^n_statin
    e_statin <- emax_statin * u_statin / (u_statin + 1)

    # ---- Nonstatin dose-responses (Equation 2) -----------------------------
    # Ezetimibe and gemcabene each have their own Emax, ED50 and n, so they
    # cannot be pooled the way the statins are.
    u_ezt <- (CONMED_EZT_DOSE / ed50_ezt)^n_ezt
    e_ezt <- emax_ezt * u_ezt / (u_ezt + 1)

    u_gem <- (CONMED_GEMCABENE_DOSE / ed50_gem)^n_gem
    e_gem <- emax_gem * u_gem / (u_gem + 1)

    # ---- Placebo intercept, with the European trial shift ------------------
    e0_arm <- e0 + e_region_europe_e0 * REGION_EUROPE

    # ---- Equation 1 --------------------------------------------------------
    # Y = E0 + Estatin + Enonstatin + 0.01 * g * Estatin * Enonstatin.
    # Mandema 2005 combined a statin with at most ONE nonstatin, so the paper's
    # single Enonstatin / g pair is written out here as two terms, one per
    # nonstatin, each with its own interaction coefficient. For every arm the
    # paper actually studied (at most one nonstatin present) this reduces
    # exactly to Equation 1. An arm coding both nonstatins at once would add
    # their effects with no ezetimibe-gemcabene interaction term, a combination
    # the paper never studied and did not parameterise.
    Cc <- e0_arm + e_statin + e_ezt + e_gem +
      0.01 * gamma_ezt * e_statin * e_ezt +
      0.01 * gamma_gem * e_statin * e_gem

    # Cc is overloaded here as the single-output observation per nlmixr2lib
    # convention; it is NOT a drug concentration but the percentage change in
    # LDL-C from the pretreatment baseline (Cc = -31.9 is a 31.9 percent
    # reduction).
    Cc ~ add(addSd)
  })
}
