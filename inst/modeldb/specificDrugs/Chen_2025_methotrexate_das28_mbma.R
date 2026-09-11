Chen_2025_methotrexate_das28_mbma <- function() {
  description <- "MBMA. Model-based meta-analysis of the DAS28 (Disease Activity Score for 28 joints) time course under methotrexate (MTX) monotherapy in adults with active rheumatoid arthritis, fitted to study-arm-mean summary data digitised from 31 randomised controlled trials (of 69 trials / 7999 patients screened into the pooled database). The endpoint is the percentage change in DAS28 from baseline, described by a three-parameter Emax-in-time model E(t) = Emax * exp(eta) * t / (ET50 * exp(eta) + t), where a SINGLE study-level random effect scales Emax and ET50 together exactly as printed in the source Model Developing equation. Emax = -54.90% and ET50 = 20.6 weeks; inter-study variability on Emax is 37.4%. No covariate was retained: MTX dose, RA disease duration, baseline CRP and baseline ESR were all screened by forward selection / backward elimination and none reduced the objective function significantly, so the final model is covariate-free (see covariatesDataExcluded). Residual error is proportional and is reported UNWEIGHTED; the source weights each study-arm observation by W = 1/sqrt(N) for arm size N, which downstream simulation code must apply. Suitable simulation scope is study-arm-mean DAS28 percentage-change trajectories, NOT individual-patient responses. Companion endpoint models from the same paper are modellib('Chen_2025_methotrexate_acr20_mbma') and modellib('Chen_2025_methotrexate_acr50_mbma')."

  reference <- paste(
    "Chen S, Wu Y, Huang W, Zhou J, Wei Z, Wu X.",
    "Can Methotrexate Monotherapy Achieve Clinical Remission in Patients with",
    "Active Rheumatoid Arthritis? A Model-Based Meta-Analysis.",
    "J Clin Pharmacol. 2025;65(10):1310-1321.",
    "doi:10.1002/jcph.70039.",
    sep = " "
  )
  vignette <- "Chen_2025_methotrexate_rheumatoid_arthritis"

  units <- list(
    time          = "week (follow-up time after starting MTX treatment; ET50 is reported in weeks)",
    dosing        = "n/a (this MBMA has no exposure driver and consumes no rxode2 dose events; MTX dose was screened as a covariate and rejected)",
    concentration = "percent/percent (Cc is the study-arm mean PERCENTAGE CHANGE FROM BASELINE in DAS28, a negative quantity; it is NOT a drug concentration. The slash satisfies checkModelConventions unit parsing.)"
  )

  # The final model retained NO covariates. All four screened covariates
  # are documented here so the paper's covariate screen is preserved
  # without triggering a "declared but not referenced" convention warning.
  covariatesDataExcluded <- list(
    DOSE_MTX_MGM2 = list(
      description        = "Per-arm average administered methotrexate dose.",
      units              = "mg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on Emax as both a power form (P = theta1 * (COV/COVmedian)^theta2) and an exponential form (P = theta1 * exp(theta2 * (COV - COVmedian))); NOT retained (Chen 2025 Results, Model Developing: 'After forward selection and backward elimination, no significant reduction in the OFV was observed, as shown in Tables S2-S6.'). The Discussion attributes the failure to non-standard MTX dosing and to most trials not reporting an average dose: 'The inability to fit the dose-effect model in this study was disappointing, likely because the majority of studies did not report the average MTX dose.' Note the register canonical DOSE_MTX_MGM2 is per body-surface area; Chen 2025 does not state the dose unit used in the covariate screen, and since no coefficient was retained the unit is immaterial to this model."
    ),
    T_DIAG_RA = list(
      description        = "Per-arm mean time since rheumatoid arthritis diagnosis (RA disease duration).",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on Emax and NOT retained (Chen 2025 Results, Model Developing; Tables S2-S6). Table 1 reports a pooled median duration of rheumatoid arthritis of 4.15 years (range 0.13-12.5). Named per the register's established T_DIAG_<disease> family (T_DIAG_DIAB, T_DIAG_CANCER)."
    ),
    CRP = list(
      description        = "Per-arm mean baseline C-reactive protein.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on Emax and NOT retained (Chen 2025 Results, Model Developing; Tables S2-S6). Table 1 reports a pooled median baseline CRP of 27.82 mg/L (range 3.1-53.9). The Discussion notes the covariate screen evaluated relative rather than absolute treatment effects, which may explain the null result."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 7999L,
    n_studies      = 31L,
    age_range      = "not reported at arm level",
    weight_range   = "not reported at arm level",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported at arm level",
    disease_state  = "adults with active rheumatoid arthritis receiving methotrexate monotherapy; trials in which MTX-inadequate responders continued on MTX were excluded, so the population is MTX-naive or MTX-responsive patients",
    dose_range     = "not tabulated; most included trials did not report an average MTX dose (Chen 2025 Discussion)",
    baseline       = "pooled medians (Chen 2025 Table 1): DAS28 5.78 (range 3.79-6.84), HAQ 1.34 (0.67-2.59), CRP 27.82 mg/L (3.1-53.9), ESR 46.56 mm/h (23-63), RA duration 4.15 years (0.13-12.5)",
    timepoints     = "arm-mean DAS28 percentage change reported at multiple follow-up times per trial; pooled median treatment duration 52.73 weeks (range 12-144)",
    regions        = "not reported at arm level; the search was restricted to English-language publications",
    notes          = "MBMA at the study-arm level: each modeled data point is one trial arm's mean DAS28 percentage change at one follow-up time. n_studies = 31 is the number of trials contributing DAS28 data; n_subjects = 7999 is the pooled patient count across all 69 trials in the database, not the DAS28 subset (Chen 2025 does not report a DAS28-specific patient count). Baseline ESR was the fourth screened covariate and was likewise NOT retained; it is documented in prose rather than in covariatesDataExcluded because the covariate register has no canonical entry for erythrocyte sedimentation rate and this file introduces no new canonical for a rejected covariate. The Results text says 'Table S1 provides an overview of the 71 studies' while the Results and Abstract both state 69 included studies; the supplement is not on disk and the discrepancy is unresolved (see vignette Errata). The model simulates study-arm mean trajectories and is NOT suitable for individual-subject simulation."
  )

  ini({
    # ============================================================
    # Emax-in-time model (Chen 2025 Methods, Model Developing; the
    # equation is rendered as a display equation in the PDF and is LOST
    # from the preprocessed markdown, recovered via pdftotext -layout):
    #
    #   E_ij = Emax * exp(eta_i) * Time_j / (ET50 * exp(eta_i) + Time_j)
    #
    # Note the single eta_i multiplies BOTH Emax and ET50 as printed.
    # Table 2 labels the random effect "eta (Emax-DAS28)", i.e. names it
    # after Emax alone, but the printed equation places exp(eta_i) on
    # ET50 as well. Per the standing "trust the printed equation over
    # the prose/label" policy this file encodes the equation as printed.
    # See the vignette Assumptions and Deviations section.
    #
    # There is NO Hill/sigmoidicity exponent: the Methods call this a
    # "sigmoid Emax" model but no gamma is reported in Table 2, and the
    # gamma = 1 form reproduces every published simulation value (see
    # the vignette source-trace table).
    #
    # All values are Chen 2025 Table 2, "Value (RSE)" column.
    # ============================================================

    emax <- -54.90
    label("Maximum percentage change from baseline in DAS28 under MTX monotherapy (paper: Emax-DAS28). Negative because DAS28 decreases (improves) from baseline. Not log-transformed because the value is negative.")  # Chen 2025 Table 2, Emax-DAS28 = -54.90% (RSE 13%); bootstrap median -52.78, 95% CI -59.43 to -41.85

    let50 <- log(20.60)
    label("Log time to half-maximal effect (paper: ET50-DAS28). Back-transformed value 20.60 weeks. Fitted in log space here to keep ET50 positive; the paper reports the value on the linear scale.")  # Chen 2025 Table 2, ET50-DAS28 = 20.60 weeks (RSE 26%); bootstrap median 19.16, 95% CI 11.47-26.02

    # ============================================================
    # Inter-study variability (ISV). Chen 2025 Methods: "Inter-study
    # variability (ISV) was described using an exponential model" and
    # "eta_i represents the random effect between studies, assumed to be
    # normally distributed with a mean of 0 and a variance of omega^2".
    # Table 2 reports eta as a percentage (37.40%), i.e. omega = 0.374
    # on the exponential scale; the ini() value is the VARIANCE.
    # Encoded as an MBMA STUDY-LEVEL eta (NOT individual between-subject
    # variability) per the SKILL Phase-1 Step-3a MBMA guidance.
    # ============================================================
    eta_study_emax ~ 0.139876  # Chen 2025 Table 2, eta(Emax-DAS28) = 37.40% (RSE 15%, shrinkage 2%); omega = 0.374, variance = 0.374^2 = 0.139876

    # ============================================================
    # Study-arm-level residual correlation (paper: ERR1, tabulated as
    # "Correlation coefficient of DAS28"). Chen 2025 Methods: the
    # residual is written as W * (ERR1 + RUV), where ERR1 is realised
    # ONCE per study arm via NONMEM's L2 data item ("Level-two (L2) data
    # item in NONMEM is used to group together the data records
    # containing observations which have the same realization of the
    # level-two effects. The RUV values were consistent across different
    # time points within the same study arm."). A term constant within an
    # arm is an eta, not an epsilon, in nlmixr2, so ERR1 is encoded here
    # as a second study-level eta on the same (proportional) scale as
    # the residual. The tabulated value carries a shrinkage figure,
    # confirming it is a NONMEM VARIANCE rather than an SD.
    # ============================================================
    eta_study_err1 ~ 0.00264  # Chen 2025 Table 2, Correlation coefficient of DAS28 = 0.00264 (RSE 28%, shrinkage 10%); NONMEM variance -> SD 0.0514 (5.14% proportional)

    # ============================================================
    # Residual error. Chen 2025 Methods lists four candidate forms
    # (additive, proportional, exponential, combined) but never states
    # which was selected per endpoint. The magnitudes in Table 2 settle
    # it: the DAS28 residual (0.00135) is ~1400x smaller than the ACR20
    # residual (1.91) although both endpoints are modelled on the same
    # percentage scale (Emax -54.90% vs 70.30%). An additive SD of
    # sqrt(0.00135) = 0.037 percentage points on a -20% prediction is
    # physically impossible for WebPlotDigitizer-digitised arm means
    # (digitisation error alone is ~1%), whereas a PROPORTIONAL CV of
    # 3.67% is reasonable; conversely the ACR20 value can only be
    # additive (a proportional reading gives CV 138%). The paper's own
    # wording supports a per-endpoint selection: the four forms were
    # candidates "with adjustments made to the model structures".
    # This file therefore encodes DAS28 as PROPORTIONAL. The alternative
    # reading -- that DAS28 was fitted on the FRACTIONAL scale
    # (Emax = -0.549) with an additive SD of 0.037, i.e. 3.7 percentage
    # points after rescaling -- is documented in the vignette Errata.
    #
    # The value is the UNWEIGHTED estimate. Chen 2025 Methods applies
    # W = 1/sqrt(N) (N = per-arm sample size) to the residual; per-arm
    # weighting is left to downstream simulation code, following the
    # Boucher_2018_naproxen_mbma / Mercier_2014_tramadol_tapentadol_mbma
    # precedent.
    # ============================================================
    propSd <- 0.0367423
    label("Proportional residual SD on the study-arm mean DAS28 percentage change (paper: epsilon-DAS28). Per-observation SD is propSd / sqrt(N_arm) per the Methods weighting W = 1/sqrt(N); the N weighting is applied downstream, not here.")  # Chen 2025 Table 2, epsilon-DAS28 = 0.00135 (RSE 93%, shrinkage 10%); NONMEM variance -> SD = sqrt(0.00135) = 0.0367423
  })

  model({
    # This MBMA has no ODE and no dosing: the endpoint is an algebraic
    # function of follow-up time only. `time` is the follow-up time in
    # weeks since the start of MTX treatment.

    # The single study-level eta scales Emax and ET50 together, exactly
    # as printed in the Chen 2025 Model Developing equation.
    emaxArm <- emax * exp(eta_study_emax)
    et50Arm <- exp(let50) * exp(eta_study_emax)

    # Study-arm mean DAS28 percentage change from baseline. Cc is the
    # canonical observation name; here it is a PERCENTAGE CHANGE (a
    # negative quantity), NOT a drug concentration -- same convention as
    # Boucher_2018_naproxen_mbma. At time 0 the fraction is 0 so Cc = 0;
    # as time grows large Cc approaches emaxArm.
    #
    # The study-arm-level residual correlation ERR1 enters proportionally,
    # on the same scale as the proportional residual, per the Methods
    # form Y = E * (1 + W * (ERR1 + RUV)).
    Cc <- emaxArm * time / (et50Arm + time) * (1 + eta_study_err1)
    Cc ~ prop(propSd)
  })
}
