Hansson_2013_sunitinib_OS <- function() {
  description <- "Parametric overall-survival (Weibull TTE) model in adults with imatinib-resistant gastrointestinal stromal tumours (GIST) on sunitinib. The hazard for death is a Weibull baseline (lam_haz, alfa_haz) modulated log-linearly by the model-predicted relative change in soluble VEGFR-3 (sVEGFR-3) from individual baseline and by observed baseline tumour size (sum of longest diameters, SLD). The sVEGFR-3 time course is simulated in-model as a one-compartment indirect-response turnover with simple-Imax inhibition of Kin driven by the per-cycle exposure summary auc = DOSE / CLI. A parallel Weibull censoring hazard (lam_cens, alfa_cens) is included so the model can drive prospective Kaplan-Meier simulations with censoring per the paper's published procedure. The model has no PK ODE and consumes individual posthoc upstream-PD parameters (BAS_SVEGFR3, MRT_SVEGFR3, EC50_SVEGFR3) and posthoc upstream-PK clearance (CLI) plus observed baseline tumour size (TUMSZ, mm) as data covariates. No IIV reported in the source for the OS or censoring hazard parameters."
  reference <- paste(
    "Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE,",
    "French J, Karlsson MO, Friberg LE.",
    "PKPD modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT as predictors of",
    "tumor dynamics and overall survival following sunitinib treatment in GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e84.",
    "doi:10.1038/psp.2013.61.",
    "Companion biomarker indirect-response and tumor growth inhibition models",
    "from the same paper are available as",
    "modellib('Hansson_2013a_sunitinib') (DDMODEL00000197) and",
    "modellib('Hansson_2013b_sunitinib') (DDMODEL00000198).",
    sep = " "
  )
  vignette <- "Hansson_2013_sunitinib_OS"
  paper_specific_compartments <- c("cumhaz_cens")
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "probability (the model outputs sur_os and sur_cens are survival probabilities, not drug concentrations)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles (e.g., the 4-weeks-on / 2-weeks-off GIST regimen) or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. For typical-cohort vignette simulations the value is held at 50 mg during the 4-week on-cycles and 0 mg during the 2-week off-cycles, matching the largest study (1004) in Hansson 2013 Table 1. The per-cycle daily-AUC equivalent (mg*h/L) is the same exposure summary used by the upstream Hansson 2013a biomarker model.",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit (Houk et al. 2009 sunitinib popPK meta-analysis, the framework PK model referenced by Hansson 2013 Methods). Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Typical sunitinib CL is around 33 L/h (consistent with the value used by the upstream Hansson 2013a/b extractions in the DDMORE bundle). For a refit or new-population simulation the user must supply individual CL drawn from a sunitinib popPK model; the Hansson 2013 paper text describes the upstream PK as a 'previously developed 2-compartment model' (reference 35 = Houk 2009).",
      source_name        = "CL"
    ),
    BAS_SVEGFR3 = list(
      description        = "Individual posthoc baseline sVEGFR-3 (pg/mL) from the upstream biomarker indirect-response PD fit (Hansson 2013a, DDMODEL00000197). Per-subject, time-fixed; used both as the initial condition for the in-model svegfr3 state and as the denominator in the relative-change driver bm_svegfr3 = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 paper reports a typical sVEGFR-3 baseline of 63900 pg/mL (Table 2). For new-population simulations either (a) simulate from modellib('Hansson_2013a_sunitinib') to obtain individual posthoc baselines, or (b) set every subject to the typical 63900 pg/mL.",
      source_name        = "BAS3"
    ),
    MRT_SVEGFR3 = list(
      description        = "Individual posthoc mean residence time of sVEGFR-3 (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; appears as kout3 = 1 / MRT_SVEGFR3 inside model().",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 paper reports a typical sVEGFR-3 MRT of 16.7 days (Table 2), which equals 400.8 h.",
      source_name        = "MRT3"
    ),
    EC50_SVEGFR3 = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sVEGFR-3 (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; appears in the drug-effect term eff3 = auc / (EC50_SVEGFR3 + auc).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 paper reports a typical (common across the four biomarkers) IC50 of 1.0 mg*h/L (Table 2).",
      source_name        = "EC53"
    ),
    TUMSZ = list(
      description        = "Observed baseline tumor size (sum of longest diameters of target lesions, SLD, mm) at study entry; per-subject, time-fixed. Enters the OS hazard via a log-linear multiplier exp(e_tumbase_haz * TUMSZ).",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Median baseline SLD across the four studies in Hansson 2013 Table 1: 194 mm (study 1004), 108 mm (study 1047), 166 mm (study 1045), 255 mm (study 013). Figure 4 caption: median baseline tumor size = 195 mm.",
      source_name        = "TUMSZ"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    age_range      = "adults with imatinib-resistant GIST (Hansson 2013 Table 1 lists baseline tumor size by study but does not break out age / weight / sex / race in the trimmed PDF Methods section)",
    weight_range   = "not reported in the on-disk paper trimmed text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Imatinib-resistant gastrointestinal stromal tumours (GIST). Pooled four sunitinib studies: Demetri 2006 (study 1004; placebo-controlled phase III; 202 active + 47 placebo), George 2009 (study 1047; phase II continuous-dosing 37.5 mg QD; n=13 in this analysis subset), Shirao 2010 (study 1045; Japanese phase I/II; 25-75 mg QD on a 4/2 schedule; n=36), Maki 2005 (study 013; phase I/II 25-75 mg QD on a 2/1 or 2/2 schedule; n=52).",
    dose_range     = "Sunitinib 25-75 mg PO QD on a 4/2, 2/2, 2/1 (weeks on / weeks off) or continuous treatment schedule. The largest cohort (study 1004) used 50 mg QD on a 4/2 schedule. Placebo arm: no sunitinib (study 1004 only).",
    regions        = "Phase III multinational (study 1004); Japanese phase I/II (study 1045); other studies regions not stated in the trimmed paper text.",
    biomarkers     = "Survival endpoint: time-to-death (overall survival, OS). Time-varying covariate for OS: model-predicted relative change in sVEGFR-3 from individual baseline (driven by an in-model sVEGFR-3 indirect-response sub-model). Time-fixed covariate for OS: observed baseline SLD (mm).",
    notes          = "n_subjects = 303 reported in Hansson 2013 Methods. Figure 4 caption reports median baseline tumor size = 195 mm and median steady-state decrease in sVEGFR-3_REL = -0.32 -- these are useful anchors for verifying typical-value simulations against the published Kaplan-Meier plot. Detailed baseline demographics (age, weight, sex, race) at the cohort level are not transcribed because the trimmed paper text does not provide them at the per-study or pooled level."
  )

  ini({
    # ----------------------------------------------------------------------
    # Overall-survival Weibull model parameter estimates from Hansson 2013
    # CPT Pharmacometrics Syst Pharmacol 2013;2:e84 Table 3 (Survival model
    # column). Parameter convention: lam_haz (lambda) reported in 1/week in
    # the paper; converted to 1/h via /24/7 so all rates carry units of
    # units$time = "hour", consistent with Hansson_2013a / Hansson_2013b.
    # ----------------------------------------------------------------------

    # Weibull baseline-hazard scale (paper Table 3: lambda = 0.00596/week).
    # 0.00596 / 24 / 7 = 3.5476e-5 / h.
    llam_haz  <- log(0.00596 / 24 / 7); label("Weibull OS baseline-hazard scale lambda (1/h; paper Table 3 = 0.00596/week, RSE 49%)")  # Table 3 Survival "lambda (per week)" = 0.00596

    # Weibull baseline-hazard shape (paper Table 3: alpha = 1.23, unitless).
    lalfa_haz <- log(1.23);             label("Weibull OS baseline-hazard shape alpha (unitless; paper Table 3 = 1.23, RSE 6.9%)")     # Table 3 Survival "alpha" = 1.23

    # Hazard log-linear coefficient on the relative change in sVEGFR-3 from
    # individual baseline (paper Eq. 6: hazard ~ exp(beta1 * sVEGFR-3_REL +
    # beta2 * Tumor_base)). Positive sign means that a less-negative
    # sVEGFR-3 relative change (smaller drop) corresponds to a higher
    # hazard, i.e., worse survival -- consistent with the paper's
    # observation that "smaller decreases in sVEGFR-3 ... were associated
    # with an increased risk of death".
    e_svegfr3_haz <- 3.77;              label("Hazard log-linear coefficient on relative change in sVEGFR-3 (unitless; paper Table 3 = 3.77, RSE 16%)")  # Table 3 Survival "beta1 sVEGFR-3" = 3.77

    # Hazard log-linear coefficient on baseline tumor size (paper units:
    # 1/mm). Positive sign means that larger baseline tumor size
    # corresponds to higher hazard / worse survival -- consistent with the
    # paper's observation that "... larger baseline tumor size were
    # associated with an increased risk of death".
    e_tumbase_haz <- 0.00237;           label("Hazard log-linear coefficient on baseline tumor size (1/mm; paper Table 3 = 0.00237, RSE 28%)")        # Table 3 Survival "beta2 Tumor base (/mm)" = 0.00237

    # ----------------------------------------------------------------------
    # Censoring Weibull (paper Figure 4 caption: "Censoring was described
    # by a Weibull model (lambda = 0.0017, alpha = 1.3) and applied in the
    # simulations"; refined values in Table 3 = 0.0017 / 1.27). Lambda
    # reported in 1/week and converted to 1/h via /24/7. No covariate
    # effects on the censoring hazard.
    # ----------------------------------------------------------------------
    llam_cens  <- log(0.0017 / 24 / 7); label("Weibull censoring-hazard scale lambda_cens (1/h; paper Table 3 = 0.0017/week, RSE 46%)") # Table 3 Survival "lambda_cens (per week)" = 0.0017
    lalfa_cens <- log(1.27);            label("Weibull censoring-hazard shape alpha_cens (unitless; paper Table 3 = 1.27, RSE 6.6%)")   # Table 3 Survival "alpha_cens" = 1.27

    # No IIV reported in the source for any of the OS or censoring
    # parameters (Table 3 Survival column has no IIV cV(%) row entries).
    # No residual-error block: the model output is the survival
    # probability sur_os = exp(-cumhaz_os) (and the censoring analogue),
    # intended for forward simulation. Following the same pattern as
    # modellib("Zecchin_2016_survival"), no observation-error model is
    # attached at the model-file level; for nlmixr2 fitting the user would
    # add the TTE likelihood `~ tte(<state>)` or similar at fit time.
  })

  model({
    # 1. Per-cycle drug-exposure summary (mg*h/L AUC). mg / (L/h) = mg*h/L;
    #    matches the convention used by the upstream Hansson 2013a/b
    #    biomarker and TGI models.
    auc <- DOSE / CLI

    # 2. sVEGFR-3 indirect-response turnover with simple-Imax drug effect
    #    on Kin (IMAX = 1 fixed, no Hill exponent -- Hansson 2013 Table 2).
    #    Per-subject Kin / Kout / EC50 come from the upstream Hansson 2013a
    #    posthoc Bayes parameters carried as data covariates.
    kout3 <- 1 / MRT_SVEGFR3
    kin3  <- BAS_SVEGFR3 * kout3
    eff3  <- auc / (EC50_SVEGFR3 + auc)

    svegfr3(0)    <- BAS_SVEGFR3
    d/dt(svegfr3) <- kin3 * (1 - eff3) - kout3 * svegfr3

    # 3. Relative change in sVEGFR-3 from individual baseline. Negative
    #    under drug (sVEGFR-3 depleted); zero at baseline; positive only
    #    if the state overshoots its baseline (does not happen under the
    #    published parameter values).
    bm_svegfr3 <- (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3

    # 4. Back-transformed Weibull parameters for OS and censoring.
    lam_haz   <- exp(llam_haz)
    alfa_haz  <- exp(lalfa_haz)
    lam_cens  <- exp(llam_cens)
    alfa_cens <- exp(lalfa_cens)

    # 5. OS Weibull hazard with log-linear covariate effects (paper Eq. 6):
    #      h(t) = lam_haz * alfa_haz * (lam_haz * t)^(alfa_haz - 1)
    #             * exp(beta1 * sVEGFR-3_REL(t) + beta2 * Tumor_base)
    #    A small time offset del = 1e-6 keeps the (lam_haz * t)^(alfa-1)
    #    term finite at t = 0 when alfa_haz < 1 without affecting the
    #    integrated cumulative hazard for typical sub-cycle event times
    #    (the same numerical device used in Zecchin_2016_survival).
    del <- 1e-6
    hazard_os <- lam_haz * alfa_haz * (lam_haz * (t + del))^(alfa_haz - 1) *
      exp(e_svegfr3_haz * bm_svegfr3 + e_tumbase_haz * TUMSZ)
    d/dt(cumhaz_os) <- hazard_os
    cumhaz_os(0) <- 0
    sur_os <- exp(-cumhaz_os)

    # 6. Censoring Weibull hazard (no covariate effects). Same Weibull
    #    parameterisation as the OS hazard; integrate to cumhaz_cens for
    #    prospective Kaplan-Meier simulations that apply censoring per the
    #    paper's published procedure.
    hazard_cens <- lam_cens * alfa_cens * (lam_cens * (t + del))^(alfa_cens - 1)
    d/dt(cumhaz_cens) <- hazard_cens
    cumhaz_cens(0) <- 0
    sur_cens <- exp(-cumhaz_cens)
  })
}
