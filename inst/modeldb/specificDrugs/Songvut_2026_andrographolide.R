Songvut_2026_andrographolide <- function() {
  description <- paste0(
    "Exploratory sigmoidal Emax exposure-response model relating ",
    "andrographolide exposure to SARS-CoV-2 viral-load reduction in patients ",
    "with mild COVID-19 treated with standardized Andrographis paniculata ",
    "aqueous extract capsules (30 mg andrographolide every 8 h, 90 mg/day, ",
    "for 5 consecutive days). The response is the log10 RdRp viral-load ",
    "reduction between day 1 and day 5 and the driver is the day-5 ",
    "andrographolide AUC(0-4 h), supplied as the AUC_ANDRO covariate column. ",
    "This is a cross-sectional exposure-response layer only: there are no ODE ",
    "states, no dosing events and no time dependence, so the pharmacokinetics ",
    "are NOT part of this file (the source paper characterised them by ",
    "non-compartmental analysis rather than by a population PK model). The ",
    "model is a four-parameter sigmoid, baseline e0 plus a saturable ",
    "increment: reduction = e0 + emax * AUC^hill / (auc50^hill + AUC^hill). ",
    "The nonzero baseline captures spontaneous viral clearance in this ",
    "uncontrolled single-arm study and is figure-derived, not printed (see ",
    "the ini() comment on le0 and the vignette Errata). Fit by unweighted ",
    "nonlinear regression in GraphPad Prism to n = 12 patients at a single ",
    "dose level, so the parameters are hypothesis-generating rather than ",
    "confirmatory: the paper reports r2 = 0.4 and could not estimate ",
    "confidence intervals for any parameter. There is no between-subject ",
    "variability (the source fitted a single typical-value curve to the ",
    "pooled cohort)."
  )
  reference <- paste(
    "Songvut P, Payoong P, Okada PA, Rittapai N, Suntararuks S,",
    "Akanimanee J, Rangkadilok N, Panomvana D, Puranajoti P, Satayavivad J.",
    "Exploratory pharmacokinetic-pharmacodynamic characterization and safety",
    "of standardized Andrographis paniculata aqueous extract capsules in",
    "patients with mild COVID-19.",
    "Front Pharmacol. 2026;17:1781740. doi:10.3389/fphar.2026.1781740.",
    "Model equation from Section 2.4 'Statistical analysis'; parameter values",
    "from Section 3.3.2 'PK/PD relationship'. The baseline term e0 is not",
    "printed anywhere in the paper and was recovered from Figure 5; see the",
    "ini() comment on le0 for the adjudication.",
    sep = " "
  )
  vignette <- "Songvut_2026_andrographolide"

  units <- list(
    time          = "day",
    dosing        = "(none; andrographolide exposure enters as the AUC_ANDRO covariate, not as a dose record)",
    concentration = "(log10 RdRp viral-load reduction between day 1 and day 5, -log10 copies/uL; the AUC_ANDRO exposure covariate is in ug*h/L)"
  )

  covariateData <- list(
    AUC_ANDRO = list(
      description        = "Per-subject andrographolide (AP1) plasma AUC over the 0-4 h post-dose window on day 5 of q8h dosing, consumed directly as the pharmacodynamic driver of the sigmoidal Emax exposure-response curve. Time-fixed per subject: the source analysis pairs one AUC value with one day-1-to-day-5 viral-load reduction, so this is a cross-sectional covariate rather than a time-varying one.",
      units              = "ug*h/L",
      type               = "continuous",
      reference_category = "n/a -- enters via the sigmoid AUC^hill / (auc50^hill + AUC^hill). AUC_ANDRO = 0 makes the sigmoid vanish exactly and recovers the baseline reduction e0. Reference values: the cohort mean day-5 AUC(0-4 h) is 30.12 +/- 15.83 ug*h/L (paper Table 2), which the paper notes approaches the estimated auc50 of 29.80 ug*h/L; the observed range spans roughly 8-57 ug*h/L (paper Figure 5).",
      notes              = "Computed in the source analysis by non-compartmental analysis (PK Solutions 2.0, linear trapezoidal rule) on LC-MS/MS plasma andrographolide concentrations sampled pre-dose and at 0.5, 0.75, 1, 1.5, 2 and 4 h after the day-5 morning dose (paper Sections 2.3.2.2 and 2.3.6). Andrographolide, not the more highly exposed AP3, was chosen as the exposure metric because it is the Thai FDA reference compound for standardized A. paniculata preparations and the most extensively characterised diterpenoid (paper Section 4.3). No population PK model exists for this drug in the source, so downstream users must supply AUC_ANDRO from observed concentrations or from an external andrographolide PK model. Member of the AUC_<DRUG> canonical family.",
      source_name        = "AUC(0-4h, day 5)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12L,
    n_studies      = 1L,
    age_range      = "18-60 years (inclusion criterion); mean 34.08 +/- 9.29 years, two participants classified as middle-aged adults",
    age_median     = "(not reported in Songvut 2026; mean 34.08 years)",
    weight_range   = "(not reported in Songvut 2026; body mass index 21.07 +/- 2.19 kg/m^2, and obesity with BMI > 35 kg/m^2 was an exclusion criterion)",
    weight_median  = NA_character_,
    sex_female_pct = 83.33,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Mild COVID-19 by the WHO 2020 clinical-management classification: RT-PCR-confirmed SARS-CoV-2 infection within 72 h of symptom onset, symptomatic without viral pneumonia or hypoxia, resting SpO2 at least 95%, normal chest radiograph. All participants were vaccinated with their most recent dose more than one year before enrollment, and none had significant comorbidities. Baseline RdRp viral load geometric mean 8.46e4 copies/uL (range 5.41e2 to 4.40e6), baseline RdRp Ct 18.54 +/- 3.37.",
    dose_range     = "Single dose level. Standardized A. paniculata aqueous extract capsules equivalent to a labelled 30 mg of andrographolide (three capsules of 10 mg labelled content; measured content 10.28 +/- 0.03 mg/capsule per paper Table 1) every 8 h, i.e. 90 mg/day, for 5 consecutive days, alongside as-needed symptomatic standard of care.",
    regions        = "Thailand (single centre: acute respiratory infection clinic, Chulabhorn Hospital, Chulabhorn Royal Academy, Bangkok)",
    notes          = "Baseline demographics and clinical laboratory values in paper Table 5; per-subject viral loads and log10 reductions in Table 4 (Ct) and Table 3 (copies/uL). All 12 participants were of Thai nationality by inclusion criterion, hence race_ethnicity is recorded as 100% Asian; the paper does not report a race or ethnicity breakdown. All 12 completed the protocol with 100% compliance, no dropouts and no protocol deviations, so all 12 contribute to the exposure-response fit. Open-label and single-arm: there was no placebo group, so the observed viral-load decline cannot be attributed to the extract, and the paper is explicit that the exposure-response analysis is exploratory and hypothesis-generating. The paper screened no covariates on the exposure-response model, so no covariate effects are encoded here."
  )

  ini({
    # ------------------------------------------------------------------
    # Sigmoidal Emax exposure-response parameters. Source equation
    # (Section 2.4):
    #   Effect = Emax * AUC^gamma / (EAUC50^gamma + AUC^gamma)
    # Source point estimates (Section 3.3.2): EAUC50 = 29.80 ug*h/L,
    # Emax = 2.24 log10 copies/mL, Hill coefficient "approximately 8",
    # r2 = 0.4. The paper states that confidence intervals could not be
    # reliably estimated for Emax, EAUC50 or the Hill coefficient given
    # n = 12, a narrow exposure range and a single dose level, so no
    # uncertainty is available for any parameter below.
    #
    # STRUCTURAL DEVIATION FROM THE PRINTED EQUATION -- read this before
    # changing anything. A baseline term e0 is added below, making the
    # encoded model the four-parameter sigmoid
    #   reduction = e0 + emax * AUC^hill / (auc50^hill + AUC^hill)
    # rather than the printed three-parameter form. Operator-ratified
    # (sidecar oare_PMC13022710 request-001, answered 2026-08-31, option
    # B). Evidence, all internal to the paper:
    #   * Figure 5 draws the fitted curve with a clearly nonzero lower
    #     asymptote near 1.85, a plateau near 4.05 carrying the "Emax"
    #     annotation, and the dashed EAUC50 guide line at y near 3.0 --
    #     not at Emax/2 = 1.12 where the printed form would put it.
    #   * e0 + emax = 1.81 + 2.24 = 4.05 reproduces the figure plateau and
    #     e0 + emax/2 = 2.93 reproduces the dashed guide line.
    #   * Scored against the paper's own 12 y-values (Table 3) and the
    #     digitised Figure 5 x-values, the printed three-parameter form
    #     gives r2 = -0.90 while the four-parameter form gives r2 = +0.37,
    #     which is the paper's reported 0.4. The printed form cannot
    #     reproduce the paper's own reported goodness of fit.
    # Most likely explanation: GraphPad Prism's sigmoidal dose-response
    # models carry Bottom and Top parameters, so the reported
    # "Emax = 2.24" is the Top-minus-Bottom increment rather than the
    # absolute plateau, and the typeset equation dropped the Bottom term.
    # ------------------------------------------------------------------

    # FIGURE-DERIVED, not printed anywhere in the paper: least-squares
    # baseline with emax, auc50 and hill held at their printed values,
    # fitted to Table 3's 12 log10 reductions against the digitised
    # Figure 5 AUC values. The digitisation is independently validated --
    # it reproduces the paper's separately reported day-5 mean +/- SD of
    # 30.12 +/- 15.83 ug*h/L (Table 2) as 30.33 +/- 15.79 without having
    # been fitted to it. Recovered value 1.809, encoded as 1.81; it is
    # also the value read directly off the Figure 5 lower asymptote.
    le0    <- log(1.81);  label("Log baseline log10 viral-load reduction at zero andrographolide exposure (-log10 copies/uL)")

    lemax  <- log(2.24);  label("Log maximum andrographolide-attributable increment in log10 viral-load reduction (-log10 copies/uL)")
    # Section 3.3.2: "Emax = 2.24 log10 copies/mL". Note the paper writes
    # copies/mL here while Table 3 tabulates the viral loads and the
    # log10 reductions in copies/uL; a log10 ratio is unit-free, so the
    # discrepancy is a reporting slip that does not affect the value.

    lauc50 <- log(29.80); label("Log andrographolide AUC(0-4 h) producing half the maximum effect increment (ug*h/L)")
    # Section 3.3.2: "EAUC50 = 29.80 ug h/L". Corroborated in Section 4.3,
    # which contrasts it with the day-5 mean AUC(0-4 h) of 30.12 ug*h/L.

    lhill  <- log(8);     label("Log Hill coefficient of the exposure-response sigmoid (unitless)")
    # Section 3.3.2: "The Hill coefficient was estimated to be
    # approximately 8 and approached the upper boundary of the
    # constrained parameter range, indicating substantial uncertainty in
    # slope estimation." Estimated (not fixed) but at a bound, so treat
    # the steepness as poorly identified. The paper does not report the
    # bound itself. Scoring hill over 2-10 against the paper's own data
    # confirms 8 as the best-supported value under the encoded form.

    # ------------------------------------------------------------------
    # Residual error. The paper reports no residual standard deviation,
    # but the value is recoverable exactly from two printed quantities
    # because "Nonlinear regression was performed without weighting"
    # (Section 2.4) means ordinary least squares, for which
    # r2 = 1 - SSE/SST:
    #   SST = sum((y_i - ybar)^2) over Table 3's 12 printed log10
    #         reductions = 31.0025
    #   SSE = (1 - 0.4) * 31.0025 = 18.6015   (r2 = 0.4, Section 3.3.2)
    #   sd  = sqrt(SSE / 12) = 1.245
    # The n denominator is used to match nlmixr2's maximum-likelihood SD
    # parameterisation; the degrees-of-freedom-corrected alternative with
    # four estimated parameters would be sqrt(18.6015 / 8) = 1.52. Note
    # r2 is printed to one significant figure ("an r2 of 0.4"), so this
    # SD is only determined to about +/- 0.05.
    # ------------------------------------------------------------------
    addSd  <- 1.245;      label("Additive residual error on the log10 viral-load reduction (-log10 copies/uL)")

    # No between-subject variability: the source fitted a single typical-
    # value curve to the pooled 12-patient cohort by nonlinear
    # regression, with no random effects and no hierarchical structure.
  })

  model({
    # Back-transform the log-scale fixed effects. There are no random
    # effects, so these are the typical values.
    e0    <- exp(le0)
    emax  <- exp(lemax)
    auc50 <- exp(lauc50)
    hill  <- exp(lhill)

    # Four-parameter sigmoidal Emax exposure-response. Algebraic and
    # cross-sectional: no ODE states, no dosing events and no time
    # dependence, so simulation event tables carry observation records
    # only, each with its AUC_ANDRO value. AUC_ANDRO = 0 makes the
    # sigmoid vanish and returns e0.
    viralLoadReduction <- e0 + emax * AUC_ANDRO^hill / (auc50^hill + AUC_ANDRO^hill)

    viralLoadReduction ~ add(addSd)
  })
}
