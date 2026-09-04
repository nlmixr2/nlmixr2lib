Chetty_2014_efalizumab_pasi <- function() {
  description <- paste(
    "First-order asymptotic disease-progression model of the Psoriasis Area and",
    "Severity Index (PASI) during efalizumab treatment, in adults with",
    "moderate-to-severe plaque psoriasis. The PASI score is treated as",
    "continuous and decays exponentially from the cohort baseline Y(0) towards",
    "a treated steady-state asymptote Yss with a progression half-life Tp,",
    "which Chetty 2014 fitted to 397 h. Written here in the equivalent",
    "indirect-response form d(pasi)/dt = kin - kout * pasi with",
    "kout = ln(2) / Tp and kin = Yss * kout, so that pasi(0) = Y(0) and",
    "pasi(inf) = kin / kout = Yss reproduce the paper's printed closed form",
    "Y(t) = Yss + (Y(0) - Yss) * exp(-ln(2) * t / Tp) exactly.",
    "The functional form is the asymptotic progression model of Holford et al.",
    "2006 (Chetty 2014 reference 24).",
    "NO DRUG TERM IS ENCODED, because none is published. Chetty 2014 states",
    "that the baseline score Y(0) is 'modulated in proportion to concentration",
    "of bound CD11a', but reports no proportionality constant, no functional",
    "form and no linking equation, and it supplies BOTH Y(0) and Yss as inputs",
    "read off the clinical study being reproduced rather than predicting Yss",
    "from exposure. The extracted model is therefore an on-treatment PASI time",
    "course whose magnitude of benefit is set by the user-supplied Yss, which",
    "is exactly how the paper applied it to its second cohort. The companion",
    "target-engagement model from the same paper is",
    "modellib('Chetty_2014_efalizumab_cd11a_qsp'); the two are not coupled here",
    "because the coupling constant is unpublished.",
    "The ini() defaults are the Gottlieb 2002 escalating-regimen cohort",
    "(Y(0) = 24.8, Yss = 14.8). For the Gordon 2003 1 mg/kg/week cohort set",
    "Y(0) = 19 and Yss = 9.5; the vignette does this by overriding ini().",
    "KNOWN DEVIATION, quantified in the vignette: the paper's own Figure 4",
    "predicted curves fall BELOW their stated Yss late in the time course",
    "(reaching about 13.9 against Yss = 14.8, and about 8.7 against Yss = 9.5)",
    "and step discontinuously upward at the weekly dose escalations, neither of",
    "which the printed equation can produce. The printed equation tracks the",
    "paper's own predicted curves to within about 5 percent over the first",
    "1000 h and drifts to roughly 12 to 14 percent by the end of follow-up.",
    sep = " "
  )
  reference <- paste(
    "Chetty M, Li L, Rose R, Machavaram K, Jamei M, Rostami-Hodjegan A,",
    "Gardner I. (2015). Prediction of the pharmacokinetics, pharmacodynamics,",
    "and efficacy of a monoclonal antibody, using a physiologically based",
    "pharmacokinetic FcRn model. Front Immunol 5:670.",
    "doi:10.3389/fimmu.2014.00670. PMCID PMC4283607.",
    "The progression equation is printed in the Methods section 'PBPK/PD MODEL",
    "TO SIMULATE mAb EFFICACY'; the trimmed markdown companion of this paper",
    "renders it as a 'formula-not-decoded' marker, so it was recovered from the",
    "PDF. Y(0) = 24.8 (CV% 10.8) and Yss = 14.8 (CV% 22) are in the same",
    "Methods paragraph and come from Gottlieb 2002 (Chetty 2014 reference 15);",
    "Y(0) = 19 and Yss = 9.5 for the second cohort come from Gordon 2003",
    "(reference 16) in the following paragraph; Tp = 397 h is in the Results",
    "section 'PBPK LINKED PD MODEL'. The model form is attributed to Holford",
    "NH, Chan PL, Nutt JG, Kieburtz K, Shoulson I. Disease progression and",
    "pharmacodynamics in Parkinson disease. J Pharmacokinet Pharmacodyn",
    "2006;33:281-311 (Chetty 2014 reference 24).",
    "Chetty 2014 has no supplementary material (EuropePMC reports",
    "hasSuppl = N) and no erratum.",
    sep = " "
  )
  vignette <- "Chetty_2014_efalizumab_cd11a_pasi"

  units <- list(
    time          = "h",
    dosing        = "(not applicable; no drug term is published, so the model carries no dosing events)",
    concentration = "(not applicable; the PASI score is a unitless clinical index on a 0-72 scale)"
  )

  compartmentData <- list(
    pasi = list(
      analyte  = "Psoriasis Area and Severity Index score",
      units    = "(unitless score, 0-72 scale)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 500L,
    n_studies      = 2L,
    age_range      = "25-50 years (simulated cohort)",
    sex_female_pct = 50,
    disease_state  = paste(
      "Adults with moderate-to-severe plaque psoriasis. Model-building cohort",
      "(Gottlieb 2002, Chetty 2014 reference 15): mean baseline PASI 24.8",
      "(CV% 10.8), falling to a mean of approximately 14.8 (CV% 22) over the",
      "last three consecutive weeks of assessment. Application cohort",
      "(Gordon 2003, reference 16): mean baseline PASI 19, with Yss set to 9.5",
      "on the basis that 1 mg/kg/week intravenously 'usually produces a change",
      "in PASI score from baseline of about 45-50%'."
    ),
    dose_range     = paste(
      "Model-building arm: efalizumab escalating 0.3 mg/kg in week 1, 0.4 in",
      "week 2, 0.6 in week 3 and 1 mg/kg/week for the following 4 weeks, each",
      "as a 1 h intravenous infusion. Application arm: 1 mg/kg/week",
      "intravenously. Dose does not enter this model: it acts only through the",
      "user-supplied Yss."
    ),
    regions        = "north European Caucasian virtual population",
    notes          = paste(
      "Chetty 2014 Simulations section: 'Predictive studies used 5 trials with",
      "100 virtual north European Caucasian Healthy Volunteers each, aged",
      "between 25 and 50 years, with an equal proportion of males and females',",
      "hence n_subjects = 500 across 5 trials; the demographic cohort is the",
      "same virtual healthy-volunteer population used for the PK simulations,",
      "with the PASI baseline and asymptote taken from the psoriasis studies.",
      "n_studies = 2 counts the clinical PASI datasets: Gottlieb 2002 for model",
      "construction and Gordon 2003 for the external application.",
      "The paper cites the escalating-regimen study as reference 14 in the",
      "Methods and reference 15 in the Results; the Figure 3 legend keys the",
      "observed series to 'Gottlieb 2002', which is reference 15",
      "(Arch Dermatol 138:591-600), so reference 15 is the escalating",
      "multiple-dose study.",
      "Chetty 2014 Discussion records the variability caveat explicitly: 'The",
      "variability on the parameters in the disease progression model may also",
      "be a limitation, since mean data were obtained from published studies",
      "and variability could not be accurately assessed.' The model is also",
      "'not designed to simulate placebo effects of the drug'."
    )
  )

  ini({
    # ========================================================================
    # Baseline PASI score, the paper's Y(0).
    # ========================================================================
    lrbase <- fixed(log(24.8))
    label("Baseline PASI score Y(0) (unitless, 0-72 scale)")
    # Chetty 2014 Methods, "PBPK/PD MODEL TO SIMULATE mAb EFFICACY": "patients
    # with moderate to severe psoriasis (mean baseline PASI score of 24.8 -
    # CV% 10.8)". This is the Gottlieb 2002 escalating-regimen cohort. For the
    # Gordon 2003 cohort the same paragraph's successor states "The mean
    # baseline was changed to 19"; override with
    # ini(lrbase = log(19), lkin = log(9.5 * log(2) / 397)).

    # ========================================================================
    # Progression rate constant. The paper parameterises the exponential by a
    # half-life, Tp, so kout = ln(2) / Tp.
    # ========================================================================
    lkout <- fixed(log(log(2) / 397))
    label("First-order PASI progression rate constant, ln(2)/Tp (1/h)")
    # Chetty 2014 Results, "PBPK LINKED PD MODEL": "A fitted Tp parameter value
    # of 397 h was obtained using parameter estimation module of the Simulator
    # and the clinical data by Gottlieb et al." Tp is the paper's only fitted
    # parameter. ln(2)/397 = 1.7459e-3 1/h, a progression half-life of 16.5
    # days. Chetty 2014 Discussion states that efalizumab acts on Y(0) "without
    # having any effect on the half-life of the disease progression (Tp)", so
    # the same Tp is used for both cohorts.

    # ========================================================================
    # Zero-order PASI production. In the indirect-response form the treated
    # asymptote is kin / kout, so kin = Yss * ln(2) / Tp. Written with Yss and
    # Tp visible in the expression rather than pre-multiplied, so the printed
    # source values can be read straight off this line.
    # ========================================================================
    lkin <- fixed(log(14.8 * log(2) / 397))
    label("Zero-order PASI production rate, Yss * ln(2)/Tp (1/h on the PASI scale)")
    # Chetty 2014 Methods, same paragraph: "PASI scores were assessed weekly
    # and a mean score of approximately 14.8 (CV% 22) was observed in the last
    # three consecutive weeks. This value was used as the Yss in the model."
    # For the Gordon 2003 cohort: "Since a 1 mg/kg/week given by iv infusion
    # usually produces a change in PASI score from baseline of about 45-50%,
    # Yss was given a value of 9.5" - override with
    # lkin = log(9.5 * log(2) / 397).

    # ========================================================================
    # Between-subject variability. Chetty 2014 reports a CV% alongside each of
    # the two cohort inputs it took from Gottlieb 2002: 10.8% on the baseline
    # score and 22% on the steady-state asymptote. Encoded as log-normal
    # variances, omega^2 = log(1 + CV^2). Because kout is fixed, the variance
    # on kin transfers unchanged to the asymptote Yss = kin / kout, so the
    # simulated Yss carries exactly the published 22% CV. No CV% is published
    # for Tp, so kout carries no eta - consistent with the Discussion, which
    # states the drug does not affect Tp and that the disease-progression
    # variability "could not be accurately assessed".
    # ========================================================================
    etalrbase ~ log(1 + 0.108^2)
    etalkin   ~ log(1 + 0.22^2)

    # ========================================================================
    # Residual error. Chetty 2014 is a simulation study and reports no
    # residual-error model for the PASI score, so this is encoded as zero
    # rather than invented. See the vignette Errata.
    # ========================================================================
    propSd_pasi <- fixed(0)
    label("Proportional residual standard deviation for the PASI score (fraction)")
  })

  model({
    rbase <- exp(lrbase + etalrbase)
    kin   <- exp(lkin + etalkin)
    kout  <- exp(lkout)

    # =====================================================================
    # First-order asymptotic progression. Chetty 2014 Methods prints the
    # closed form
    #        Y(t) = Yss + (Y(0) - Yss) * exp(-ln(2) / Tp * t)
    # which is the solution of d(Y)/dt = -ln(2)/Tp * (Y - Yss) started at
    # Y(0). Written here as the equivalent indirect-response balance with
    # kout = ln(2)/Tp and kin = Yss * kout so that both canonical turnover
    # parameters carry their usual meaning; the two forms are algebraically
    # identical and the vignette gates the ODE solve against the closed form.
    # No drug term appears because the paper publishes none - see description.
    # =====================================================================
    d/dt(pasi) <- kin - kout * pasi

    pasi(0) <- rbase

    # The paper's Yss, recovered as the asymptote of the balance above, and
    # the maximum attainable improvement from baseline.
    pasiSs   <- kin / kout
    pasiDrop <- rbase - pasiSs

    pasi ~ prop(propSd_pasi)
  })
}
