Elhefnawy_2023_recurrent_ischemic_stroke <- function() {
  description <- paste(
    "Parametric time-to-event model for recurrent ischemic stroke (IS) after a",
    "first (index) IS, developed in NONMEM 7.5 from the National Neurology",
    "Registry of Malaysia (7,697 patients, 333 recurrences, up to 7.37 years of",
    "follow-up). The hazard is Gompertz, h(t) = theta_x * exp(theta_y * t), with",
    "the scale switching at 0.5 years (theta1 -> theta3) and the shape switching",
    "at 3 years (theta2 -> theta4). Four baseline comorbidity / secondary-",
    "prevention covariates act log-linearly on the hazard: hyperlipidemia,",
    "ischemic heart disease and hypertension raise it, receiving an antiplatelet",
    "at discharge lowers it. Time is in years and there is no drug exposure",
    "term, so this is a disease-progression / event-risk model rather than a",
    "PK-PD model. The cumulative hazard is evaluated in closed form (the",
    "piecewise Gompertz integral is elementary) rather than by an ODE, which",
    "avoids integrating across the two hazard discontinuities; the model exposes",
    "the instantaneous hazard `hazard_ris`, the cumulative hazard `cumhaz_ris`",
    "and the recurrence-free survivor probability `sur`.",
    "IMPORTANT -- the published parameters do not reproduce the paper's own",
    "outputs. As printed, theta1 = 0.238/year over the first 6 months implies a",
    "6-month recurrence probability of about 11.9 percent (17 percent with the",
    "Gompertz shape included), against the 108 of 7,697 = 1.40 percent actually",
    "reported in Elhefnawy 2023 Table 1, and the model overpredicts the",
    "cumulative hazard of the paper's own Kaplan-Meier VPC (Figure 3) by roughly",
    "4- to 8-fold for a covariate-free patient, or 6- to 12-fold once the cohort",
    "covariate distribution is applied. The discrepancy is not removable by any",
    "reading of the piecewise structure, because a merely constant 0.238/year is",
    "already 8.5-fold too high. The values are transcribed exactly as published",
    "and have NOT been tuned; see the validation vignette's Errata section for",
    "the full arithmetic.",
    sep = " "
  )
  reference <- paste(
    "Elhefnawy ME, Sheikh Ghadzi SM, Albitar O, Tangiisuran B, Zainal H,",
    "Looi I, Sidek NN, Aziz ZA, Harun SN. Predictive model of recurrent ischemic",
    "stroke: model development from real-world data.",
    "Front Neurol. 2023;14:1118711. doi:10.3389/fneur.2023.1118711.",
    "Structural equations are Equations 1-5 and the unnumbered covariate",
    "equation in Methods; the piecewise scale / shape switching is defined in",
    "the Table 2 footnote; all parameter estimates are Table 3. The covariate",
    "screen (univariate testing, forward inclusion, backward elimination, each",
    "as a change in objective function value) is Supplementary Table 1_S, the",
    "single file in the EuropePMC open-access supplementary deposit for",
    "PMC10176964; it contains no control stream and no parameter values.",
    sep = " "
  )
  vignette <- "Elhefnawy_2023_recurrent_ischemic_stroke"

  units <- list(
    time          = "year",
    dosing        = "n/a (no dosing events; secondary-prevention therapy enters only through the binary CONMED_ANTIPLATELET covariate)",
    concentration = "n/a (the model outputs are a hazard in 1/year, a unitless cumulative hazard and a unitless recurrence-free survivor probability, not a drug concentration)"
  )

  covariateData <- list(
    DIS_HYPERLIP = list(
      description        = "1 = the patient carried a hyperlipidemia (HPLD) diagnosis before the index ischemic stroke; 0 = no hyperlipidemia. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hyperlipidemia before the index stroke)",
      notes              = "Ascertained by physician diagnosis, the patient's electronic record, or medication history (Elhefnawy 2023 Methods, 'Collected variables'). Prevalence in the full 7,697-patient cohort is 2,028 / 7,697 = 26.34 percent (Results text; Table 1 gives 159 of 333 recurrent and 1,869 of 7,364 non-recurrent, which sum to the same 2,028). The strongest single predictor retained: HR = exp(0.799) = 2.22 (95 percent CI 1.81-2.72).",
      source_name        = "HPLD"
    ),
    DIS_IHD = list(
      description        = "1 = the patient carried an ischemic heart disease (IHD) diagnosis before the index ischemic stroke; 0 = no IHD. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ischemic heart disease before the index stroke)",
      notes              = "Ascertained the same way as the other comorbidity flags (Elhefnawy 2023 Methods, 'Collected variables'). Prevalence in the full cohort is 879 / 7,697 = 11.42 percent (Results text; Table 1 gives 77 + 802 = 879). HR = exp(0.745) = 2.10 (95 percent CI 1.64-2.69). Elhefnawy 2023 Discussion attributes the effect to shared atherosclerotic pathophysiology between IHD and ischemic stroke.",
      source_name        = "IHD"
    ),
    DIS_HYPERT = list(
      description        = "1 = the patient carried a hypertension (HTN) diagnosis before the index ischemic stroke; 0 = no hypertension. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hypertension before the index stroke)",
      notes              = "Ascertained the same way as the other comorbidity flags. Prevalence in the full cohort is 5,506 / 7,697 = 71.5 percent (Results text; Table 1 gives 288 + 5,218 = 5,506). HR = exp(0.711) = 2.03 (95 percent CI 1.52-2.71). Elhefnawy 2023 Table 1 also stratifies hypertension duration at 5 years, but duration was not retained in the final model -- only the presence / absence flag is used here.",
      source_name        = "HTN"
    ),
    CONMED_ANTIPLATELET = list(
      description        = "1 = the patient was prescribed an antiplatelet (APLT) at discharge from the index ischemic stroke admission, for secondary prevention; 0 = no antiplatelet prescribed. Time-fixed per subject.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no antiplatelet prescribed at discharge)",
      notes              = "Elhefnawy 2023 Methods, 'Collected variables': the secondary-prevention medications were 'prescribed during discharge'. The paper does not enumerate which agents were pooled into the antiplatelet class, so the class composition is unknown here (contrast CONMED_DIURETIC, where each source paper's class membership is recorded). Prevalence in the full cohort is (285 + 6,613) / 7,697 = 89.6 percent (Table 1). This is the only protective term retained: HR = exp(-0.514) = 0.59 (the paper prints the CI in descending order as '0.79-0.44'), i.e. about a 40 percent reduction in the hazard of recurrence. Because the indicator is fixed at discharge, it encodes prescription rather than adherence, and the model carries no time-varying exposure term -- Elhefnawy 2023 Limitations flags incorporating time-varying secondary-prophylaxis effects as future work.",
      source_name        = "APLT"
    )
  )

  # Covariates that Elhefnawy 2023 collected and screened but did not retain in
  # the final model. Listed here for provenance only; none is referenced in
  # model(). The paper reports no point estimate for any of them, so no value
  # could be transcribed even if one were wanted -- the screen is reported only
  # as changes in objective function value.
  #
  # The screen itself is Supplementary Table 1_S, "Univariate and multivariate
  # analysis of covariate effects on the hazard of recurrent IS after index IS",
  # the single file in the EuropePMC open-access supplementary deposit for
  # PMC10176964. Fifteen candidates were tested univariately against a base OFV
  # of 2808.68; the significant ones went through stepwise forward inclusion
  # (base OFV 2737.615) and then backward elimination. The Table 1_S footnote
  # sets the thresholds: "Significance; p value < 0.05 in univariate analysis
  # and stepwise forward inclusion. Significance < 0.01 in backward
  # elimination." That stricter backward threshold is what reduces the model to
  # the four covariates of Table 3 -- see DIS_DIAB below.
  #
  # These names are documentation labels for the paper's screen, not registered
  # canonical covariate columns: none is used in model(), so none is added to
  # inst/references/covariate-columns.md.
  covariatesDataExcluded <- list(
    DIS_DIAB = list(
      description = "Diabetes mellitus before the index stroke (3,493 / 7,697 = 45.38 percent).",
      units = "(binary)", type = "binary",
      notes = "The one covariate that reached the final elimination step and was still dropped. Strongly significant univariately (Suppl. Table 1_S, dOFV -21.75, p < 0.0001) and retained through forward inclusion (dOFV -3.96, p = 0.046, clearing the p < 0.05 forward threshold), but removing it in backward elimination cost only dOFV +3.88 (p = 0.048), short of the stricter p < 0.01 backward criterion, so it is absent from Table 3. Table 1 additionally stratifies diabetes duration into <1, 1-5, 6-10 and >10 years; no duration effect was retained either."
    ),
    DIS_HYPERURICEMIA = list(
      description = "Hyperuricemia (HU) before the index stroke (234 / 7,697 = 3.04 percent).",
      units = "(binary)", type = "binary",
      notes = "Univariately significant (Suppl. Table 1_S, dOFV -4.65, p = 0.031) and carried into forward inclusion, where it added almost nothing (dOFV -1.055, p = 0.304) and was not retained."
    ),
    DIS_AF = list(
      description = "Atrial fibrillation before the index stroke (about 3.4 percent of the cohort).",
      units = "(binary)", type = "binary",
      notes = "Not significant univariately (Suppl. Table 1_S, dOFV -0.44, p = 0.507) and never entered the stepwise procedure. Note the direction in Table 1 is opposite to the usual clinical expectation (1.2 percent of recurrent vs 3.57 percent of non-recurrent patients)."
    ),
    SEXF = list(
      description = "Female sex (4,289 / 7,697 = 55.72 percent).",
      units = "(binary)", type = "binary",
      notes = "Screened as 'Gender' in Suppl. Table 1_S and among the weakest candidates tested (dOFV -0.435, p = 0.509); not retained. Table 1 shows near-identical proportions in the recurrent (55.85 percent) and non-recurrent (55.71 percent) groups."
    ),
    FAMHX_STROKE = list(
      description = "Family history of stroke (FHOS).",
      units = "(binary)", type = "binary",
      notes = "Screened univariately in Suppl. Table 1_S (dOFV -2.735, p = 0.09) and not carried forward. Prevalence is not reported in Table 1."
    ),
    NIHSS = list(
      description = "National Institutes of Health Stroke Scale severity of the index stroke, dichotomised by the paper into minor vs moderate/severe.",
      units = "(score)", type = "continuous",
      notes = "Tabulated in Elhefnawy 2023 Table 1, defined in the Table 3 footnote, and screened univariately in Suppl. Table 1_S (dOFV -1.103, p = 0.293); no NIHSS term appears in the final model."
    ),
    CONMED_ANTIDIABETIC = list(
      description = "Antidiabetic (ADM) prescribed at discharge from the index stroke admission.",
      units = "(binary)", type = "binary",
      notes = "Univariately significant (Suppl. Table 1_S, dOFV -5.39, p = 0.0202) and carried into forward inclusion, where it was the weakest candidate tested (dOFV -0.167, p = 0.682) and was not retained."
    ),
    CONMED_DIURETIC = list(
      description = "Diuretic (DIU) prescribed at discharge (5.9 percent of the cohort).",
      units = "(binary)", type = "binary",
      notes = "Screened univariately in Suppl. Table 1_S (dOFV -2.87, p = 0.09); did not reach the p < 0.05 threshold for forward inclusion."
    ),
    CONMED_BETABLOCKER = list(
      description = "Beta-blocker (BB) prescribed at discharge (10.6 percent of the cohort).",
      units = "(binary)", type = "binary",
      notes = "Screened univariately in Suppl. Table 1_S (dOFV -2.05, p = 0.152); not carried forward."
    ),
    CONMED_CCB = list(
      description = "Calcium-channel blocker (CCB) prescribed at discharge (20.8 percent of the cohort).",
      units = "(binary)", type = "binary",
      notes = "Screened univariately in Suppl. Table 1_S (dOFV -1.52, p = 0.217); not carried forward."
    ),
    CONMED_ACEI = list(
      description = "Angiotensin-converting-enzyme inhibitor (ACEI) prescribed at discharge (31.1 percent of the cohort).",
      units = "(binary)", type = "binary",
      notes = "The weakest candidate in the whole screen (Suppl. Table 1_S, dOFV -0.03, p = 0.862); not carried forward."
    ),
    AGE = list(
      description = "Age at the index ischemic stroke (median 63.47 years).",
      units = "year", type = "continuous",
      notes = "Named as a screened demographic covariate in Elhefnawy 2023 Methods ('Based on demographic data and concomitant diseases') and dichotomised at 60 years in Table 1, but it does not appear among the fifteen candidates tabulated in Suppl. Table 1_S, so no objective-function change is available for it. Not retained."
    ),
    SMOKER = list(
      description = "Current smoker at the index stroke (about 48 percent of the cohort).",
      units = "(binary)", type = "binary",
      notes = "Reported in Elhefnawy 2023 Results and Table 1 with a sizeable unadjusted imbalance (60.66 percent of recurrent vs 48.17 percent of non-recurrent patients), but like AGE it is absent from the Suppl. Table 1_S screen, so no objective-function change is available for it. Not retained."
    )
  )

  population <- list(
    n_subjects     = 7697L,
    n_events       = 333L,
    n_studies      = 1L,
    age_range      = "adults aged over 18 years; median 63.47 years at the index stroke; 4,623 of 7,697 (60.1 percent) were 60 years or older",
    sex_female_pct = 55.72,
    race_ethnicity = "Multiethnic Malaysian registry cohort. Elhefnawy 2023 Table 1 reports Malay, Chinese, Indian and 'Others' strata separately for the recurrent (46.54 / 2.10 / 0.90 / 50.15 percent) and non-recurrent (20.08 / 2.79 / 1.08 / 76.05 percent) groups; ethnicity was screened but not retained in the final model.",
    disease_state  = "Adults with a first (index) ischemic stroke diagnosed by WHO criteria and confirmed by brain CT or MRI. The endpoint is a subsequent ischemic stroke recorded by any participating hospital.",
    dose_range     = "n/a (no drug exposure is modelled; secondary prevention enters only as the binary antiplatelet-at-discharge indicator)",
    regions        = "Malaysia -- National Neurology Registry (NNEUR), a multicentre hospital-based registry covering 13 states; index strokes registered August 2009 to December 2016.",
    notes          = paste(
      "333 of 7,697 patients (4.32 percent) had at least one recurrent ischemic",
      "stroke within the maximum 7.37 years of follow-up; 108 of those 333 (31.43",
      "percent) recurred within the first 6 months and 36 patients went on to a",
      "second recurrence. Median time to first recurrence was 1.2 years.",
      "Baseline comorbidity prevalences in the full cohort: hypertension 71.5",
      "percent, diabetes 45.38 percent, hyperlipidemia 26.34 percent, ischemic",
      "heart disease 11.42 percent, atrial fibrillation about 3.4 percent, current",
      "smoking about 48 percent. Secondary-prevention prescribing at discharge:",
      "antiplatelet 89.6 percent, antihyperlipidemic 89.6 percent, ACE inhibitor",
      "31.1 percent, calcium-channel blocker 20.8 percent, antidiabetic 31.7",
      "percent, beta-blocker 10.6 percent, diuretic 5.9 percent (computed from the",
      "Table 1 recurrent and non-recurrent counts). Non-Malaysian citizens and",
      "non-ischemic stroke diagnoses were excluded. Elhefnawy 2023 Limitations:",
      "the first stroke captured by the registry was assumed to be the patient's",
      "first ever stroke, and comorbidities were analysed independently of one",
      "another.",
      sep = " "
    )
  )

  ini({
    # ==================================================================
    # All values are Elhefnawy 2023 Table 3, "Parameters of the final
    # developed model for recurrent IS after index IS". The hazard form
    # is Table 2 row 4 (the 4-parameter Gompertz), h(t) = theta_x *
    # exp(theta_y * t), and the Table 2 footnote fixes the switching:
    # "theta_x equals theta1 if time < 0.5 year; theta3 if time >= 0.5,
    # theta_y equals theta2 if time < 3 years, theta4 if time >= 3
    # years." The covariates enter through the unnumbered Methods
    # equation h(t) = h0(t) * exp(beta1*X1 + ... + betan*Xn).
    #
    # The two baseline-hazard scales are strictly positive rates and are
    # carried on the log scale; the two Gompertz shapes are exponents
    # that could in principle take either sign and are carried linearly.
    #
    # INTERNAL CONSISTENCY OF THE COVARIATE BLOCK (all four reproduce):
    #   exp( 0.799) = 2.223 vs the published aHR 2.22
    #   exp( 0.745) = 2.106 vs the published aHR 2.10
    #   exp( 0.711) = 2.036 vs the published aHR 2.03
    #   exp(-0.514) = 0.598 vs the published aHR 0.59
    # The Table 3 "Half-life (Ln2/alpha)" column also reproduces:
    #   ln(2)/1.63 = 0.425 year (Table 3: "0.42 (5.06 months)")
    #   ln(2)/0.23 = 3.014 year (Table 3: "3.008 years")
    #
    # KNOWN DISCREPANCY IN THE BASELINE-HAZARD SCALE -- NOT TUNED.
    # theta1 = 0.238/year cannot be reconciled with the paper's own
    # event counts. Elhefnawy 2023 Table 1 reports 108 recurrences
    # within 6 months among 7,697 patients at risk, i.e. a cumulative
    # hazard of about 0.0141 at t = 0.5 year. A flat 0.238/year gives
    # 0.238 * 0.5 = 0.119 (8.5-fold too high) and the Gompertz form
    # gives (0.238/1.63) * (exp(1.63*0.5) - 1) = 0.184 (13-fold too
    # high) before the covariate multiplier -- whose cohort average,
    # 1.66, makes it worse rather than better. The same overprediction
    # holds against the Figure 3 Kaplan-Meier VPC at every read point
    # (about 8x at 0.5 year, 6x at 3 years, 4x at 7.37 years). No
    # reading of the piecewise structure removes it, and a frailty term
    # cannot either: a zero-mean random effect on the log hazard leaves
    # at least half the cohort at or above the typical hazard, which
    # bounds the population survivor function well below the published
    # curve. The published values are transcribed verbatim regardless;
    # the vignette's Errata section carries the full arithmetic.
    # ==================================================================

    lh0_early_ris <- log(0.238)
    label("Log Gompertz baseline-hazard scale over the first 6 months after the index stroke (1/year)")
    # Elhefnawy 2023 Table 3 row 1: theta1 (<6 months) = 0.238, RSE 19.92 percent. Abstract: "Within the first 6 months after the index IS, the hazard of recurrent IS was predicted to be 0.238". See the scale discrepancy note above.

    lh0_late_ris <- log(0.0016)
    label("Log Gompertz baseline-hazard scale from 6 months after the index stroke onward (1/year)")
    # Elhefnawy 2023 Table 3 row 2: theta3 (>=6 months) = 0.0016, RSE 21.62 percent. The Abstract rounds this to 0.001 ("6 months after the index attack, it reduced to 0.001"); the 4-decimal Table 3 value is used here.

    shape_early_ris <- 1.63
    label("Gompertz shape (hazard exponent) during the first 3 years after the index stroke (1/year)")
    # Elhefnawy 2023 Table 3 row 3: alpha (<3) = theta2 = 1.63, RSE 4.81 percent, "Shape parameter in the first 3 years after index IS". Positive, so the hazard rises within each interval; Results: "the exponential increase in the hazard of recurrent IS was observed in the first 3 years after the index IS".

    shape_late_ris <- 0.23
    label("Gompertz shape (hazard exponent) from 3 years after the index stroke onward (1/year)")
    # Elhefnawy 2023 Table 3 row 4: alpha (>=3) = theta4 = 0.23, RSE 20.19 percent, "Shape parameter after 3 years of index IS". Printed as positive. Note the tension with the Results sentence "and then exponentially reduced afterward": with theta4 = +0.23 the hazard still rises after 3 years, just 7-fold more slowly, and the "reduction" is the downward jump at t = 3 that the scale / shape switching produces (exp(1.63*3) = 133 falls to exp(0.23*3) = 2.0). Encoded with the printed sign; see the vignette Errata.

    e_dis_hyperlip_ris <- 0.799
    label("Log-hazard coefficient for pre-index hyperlipidemia; HR = exp(0.799) = 2.22")
    # Elhefnawy 2023 Table 3 row 5: theta5 = 0.799, RSE 12.89 percent, aHR 2.22 (95 percent CI 1.81-2.72).

    e_dis_ihd_ris <- 0.745
    label("Log-hazard coefficient for pre-index ischemic heart disease; HR = exp(0.745) = 2.10")
    # Elhefnawy 2023 Table 3 row 6: theta6 = 0.745, RSE 16.85 percent, aHR 2.10 (95 percent CI 1.64-2.69).

    e_dis_hypert_ris <- 0.711
    label("Log-hazard coefficient for pre-index hypertension; HR = exp(0.711) = 2.03")
    # Elhefnawy 2023 Table 3 row 7: theta7 = 0.711, RSE 20.62 percent, aHR 2.03 (95 percent CI 1.52-2.71).

    e_conmed_antiplatelet_ris <- -0.514
    label("Log-hazard coefficient for an antiplatelet prescribed at discharge; HR = exp(-0.514) = 0.59")
    # Elhefnawy 2023 Table 3 row 8: theta8 = -0.514, RSE 28.41 percent, aHR 0.59 (the paper prints the interval descending, "0.79-0.44"). Discussion: "receiving APLT was found to decrease the hazard of recurrent IS by ~40%".

    # Interval boundaries. These are structural design choices, not
    # estimated quantities: Elhefnawy 2023 Table 2 counts the interval
    # Gompertz model as having 4 parameters (theta1-theta4), so the two
    # breakpoints carry no degrees of freedom and are fixed here.
    tbreak_scale_ris <- fixed(0.5)
    label("Time after the index stroke at which the baseline-hazard scale switches from theta1 to theta3 (year)")
    # Elhefnawy 2023 Table 2 footnote: "theta_x equals, theta1 if time < 0.5 year; theta3 if time >= 0.5". Table 3 labels the same split "<6 months" / ">=6 months".

    tbreak_shape_ris <- fixed(3)
    label("Time after the index stroke at which the Gompertz shape switches from theta2 to theta4 (year)")
    # Elhefnawy 2023 Table 2 footnote: "theta_y equals; theta2 if time < 3 years, theta4 if time >= 3 years". Table 3 labels the same split "alpha (<3)" / "alpha (>=3)".

    # Between-subject variability. Elhefnawy 2023 Methods, "Development
    # of the base model", states "Between-subject variability around the
    # hazard was estimated, assuming an exponential distribution for the
    # random effect", but Table 3 reports no variance, no CV percent and
    # no shrinkage for it. The supplementary material was retrieved (the
    # EuropePMC open-access deposit for PMC10176964 holds a single file,
    # Suppl. Table 1_S) and is a covariate model-building table that
    # reports only objective-function changes -- it carries no variance
    # estimate either, so the magnitude is unavailable from every source
    # in hand. The structural form is therefore known and the
    # magnitude is not, so
    # the eta is declared on the log baseline hazard (the parameter the
    # paper says carried it) and fixed at zero rather than invented.
    # Simulations from this file are typical-value trajectories; see the
    # vignette's Assumptions and deviations section.
    #
    # The source describes ONE random effect on the hazard, but the
    # hazard scale is carried by two parameters (one per time interval),
    # so one eta is declared against each to keep the eta<x> / x naming
    # pairing intact. Both are fixed at zero, so the split is immaterial
    # as shipped; anyone reinstating the paper's variability should draw
    # a single value and use it for both.
    etalh0_early_ris ~ fixed(0)
    etalh0_late_ris ~ fixed(0)

    # The source fits this model with the parametric survival (event-
    # density) likelihood under LAPLACE, so there is no observation-error
    # model to translate. This placeholder additive residual is attached
    # to the survivor-probability output so the nlmixr2 likelihood
    # machinery accepts the model for forward simulation. It is NOT from
    # the source. Same device as Lindauer_2017_lacosamide_dropout.R and
    # Knebel_2012_istradefylline_dizziness.R.
    addSd <- 0.001
    label("Placeholder additive residual error on the survivor-probability output sur (unitless); not from the source")
  })

  model({
    # --- Baseline-hazard scales, back-transformed. The declared-but-
    # --- unreported between-subject random effect sits on the log scale
    # --- and is shared by both intervals (one hazard, one eta).
    h0_early_ris <- exp(lh0_early_ris + etalh0_early_ris)
    h0_late_ris  <- exp(lh0_late_ris  + etalh0_late_ris)

    # --- Covariate multiplier. Elhefnawy 2023 Methods, unnumbered
    # --- equation after Equation 5:
    # ---   h(t) = h0(t) * exp(beta1*X1 + beta2*X2 + .... + betan*Xn)
    # --- All four retained covariates are 0/1 indicators, so exp() of
    # --- each coefficient is the adjusted hazard ratio the paper
    # --- tabulates and plots in Figure 2.
    cov_ris <-
      exp(e_dis_hyperlip_ris        * DIS_HYPERLIP +
          e_dis_ihd_ris             * DIS_IHD +
          e_dis_hypert_ris          * DIS_HYPERT +
          e_conmed_antiplatelet_ris * CONMED_ANTIPLATELET)

    # --- Instantaneous baseline hazard, Elhefnawy 2023 Table 2 row 4
    # --- combined with the Table 2 footnote. The scale switches at
    # --- tbreak_scale_ris (0.5 year) and the shape switches
    # --- independently at tbreak_shape_ris (3 years), which gives three
    # --- regimes rather than two. Both switches are downward jumps in
    # --- the hazard.
    if (t < tbreak_scale_ris) {
      h0_ris <- h0_early_ris * exp(shape_early_ris * t)
    } else if (t < tbreak_shape_ris) {
      h0_ris <- h0_late_ris * exp(shape_early_ris * t)
    } else {
      h0_ris <- h0_late_ris * exp(shape_late_ris * t)
    }
    hazard_ris <- h0_ris * cov_ris

    # --- Cumulative baseline hazard in closed form. Equation 1 defines
    # --- S(t) = exp(-integral of h over 0..t); the piecewise Gompertz
    # --- integrates elementally, so the integral is evaluated exactly
    # --- rather than by an ODE. Doing it this way keeps the two hazard
    # --- discontinuities out of the solver, where they would otherwise
    # --- be stepped over and silently smeared.
    # ---   segment 1, 0 .. 0.5 yr : h0_early * exp(shape_early * t)
    # ---   segment 2, 0.5 .. 3 yr : h0_late  * exp(shape_early * t)
    # ---   segment 3, 3 yr ..     : h0_late  * exp(shape_late  * t)
    seg1_full_ris <-
      h0_early_ris / shape_early_ris *
      (exp(shape_early_ris * tbreak_scale_ris) - 1)
    seg2_full_ris <-
      h0_late_ris / shape_early_ris *
      (exp(shape_early_ris * tbreak_shape_ris) -
         exp(shape_early_ris * tbreak_scale_ris))

    if (t < tbreak_scale_ris) {
      cumhaz0_ris <-
        h0_early_ris / shape_early_ris * (exp(shape_early_ris * t) - 1)
    } else if (t < tbreak_shape_ris) {
      cumhaz0_ris <-
        seg1_full_ris +
        h0_late_ris / shape_early_ris *
        (exp(shape_early_ris * t) - exp(shape_early_ris * tbreak_scale_ris))
    } else {
      cumhaz0_ris <-
        seg1_full_ris + seg2_full_ris +
        h0_late_ris / shape_late_ris *
        (exp(shape_late_ris * t) - exp(shape_late_ris * tbreak_shape_ris))
    }

    # --- The covariate multiplier is time-constant, so it factors
    # --- straight out of the integral (proportional hazards).
    cumhaz_ris <- cumhaz0_ris * cov_ris

    # --- Elhefnawy 2023 Equation 1: the probability of not experiencing
    # --- a recurrent ischemic stroke within [0, t].
    sur <- exp(-cumhaz_ris)

    sur ~ add(addSd)
  })
}
