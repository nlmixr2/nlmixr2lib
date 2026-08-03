Viberg_2012_AZD6088_rat <- function() {
  description <- paste(
    "Preclinical (rat). Direct-effect population PK/PD model for the",
    "muscarinic agonist AZD6088 (MW 406.57) in male Sprague-Dawley",
    "rats using the spinal-nerve-ligation (SNL) heat-hyperalgesia",
    "model of neuropathic pain. One-compartment oral PK (depot ->",
    "central) with the first-order absorption rate constant ka",
    "fixed at 5 x kel (equivalent to 5 x cl/vc; Viberg 2012 Results",
    "PK-model paragraph: 'Lack of plasma concentration data in the",
    "absorption phase made it impossible to estimate rate of",
    "absorption and ka was therefore fixed to 5 time clearance');",
    "the 40 umol/kg highest-dose cohort (DOSE_HIGH = 1) gates a",
    "distinct typical-value apparent oral clearance of 3.92 L/h/kg",
    "vs 10.87 L/h/kg for all other cohorts, per the paper's",
    "OFV-supported categorical CL split (delta OFV = -32). PD is a",
    "direct-effect Emax on paw withdrawal latency (baseline 5.947 s,",
    "Emax 10.09 s, EC50 0.0433 umol/L; no time delay, no sigmoidicity,",
    "no tolerance-development term - all three alternatives increased",
    "OFV or degraded goodness-of-fit). Typical-value only: the paper",
    "reports 'exponential models were used to describe inter-individual",
    "variability' on CL and V but never publishes the omega^2 / CV%",
    "magnitudes; see the vignette Assumptions and deviations section."
  )
  reference <- paste(
    "Viberg A, Martino G, Lessard E, Laird JMA. (2012).",
    "Evaluation of an Innovative Population Pharmacokinetic-Based",
    "Design for Behavioral Pharmacodynamic Endpoints.",
    "The AAPS Journal 14(4):657-663.",
    "doi:10.1208/s12248-012-9380-3.",
    sep = " "
  )
  vignette <- "Viberg_2012_AZD6088_rat"
  units <- list(
    time          = "h",
    dosing        = "umol/kg",
    concentration = "umol/L"
  )

  covariateData <- list(
    DOSE_HIGH = list(
      description        = paste(
        "Binary indicator flagging membership in the study's",
        "highest-dose cohort (40 umol/kg oral, ~16.26 mg/kg at MW",
        "406.57). 1 = subject received 40 umol/kg AZD6088 orally,",
        "0 = subject received any lower dose (1, 2.5, 5, 10, or 20",
        "umol/kg in the initial-efficacy study; 2.5 umol/kg in the",
        "pilot study; the three drug-treated main-study cohorts",
        "each received a per-cohort fixed dose below the 40",
        "umol/kg threshold)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = paste(
        "Time-fixed per rat (each animal received exactly one dose",
        "level for the full study). Gates the categorical apparent",
        "oral clearance split in model(): DOSE_HIGH = 1 selects",
        "the typical-value CL/F = 3.92 L/h/kg (paper Table I, 'Oral",
        "clearance high dose'), DOSE_HIGH = 0 selects the typical-",
        "value CL/F = 10.87 L/h/kg (paper Table I, 'Oral",
        "clearance'). The paper attributes the shift to either",
        "higher bioavailability or saturation of elimination at the",
        "40 umol/kg dose (Results, PK-model paragraph); the",
        "categorical CL encoding accommodates either mechanism",
        "without committing to one. Threshold value 40 umol/kg =",
        "40 * 406.57 / 1000 = 16.26 mg/kg. Paper alias: the source",
        "does not name a NONMEM column for this indicator (the",
        "supplement / control stream is not on disk)."
      ),
      source_name        = "DOSE_HIGH"
    )
  )

  population <- list(
    species         = "rat (Sprague-Dawley, male)",
    n_subjects      = 42L,
    n_studies       = 1L,
    age_range       = "Adult (age not reported; body weight 125-200 g at receipt)",
    weight_range    = "125-200 g at receipt (Charles River St. Constant, Canada; Harlan Inc., Indianapolis, USA)",
    sex_female_pct  = 0,
    race_ethnicity  = NA,
    disease_state   = paste(
      "Male Sprague-Dawley rats subjected to left L5 + L6 spinal-nerve",
      "ligation (SNL) under isoflurane anaesthesia (Kim and Chung",
      "method) to induce heat hyperalgesia; testing performed on",
      "days 20-27 post-surgery. SNL rats included in the study had a",
      "baseline paw withdrawal latency <= 8 s on the affected paw",
      "(pre-treatment mean 5.87 +- 0.2 s), vs the naive-rat reference",
      "of 10.50 +- 0.5 s. Selected animals were randomly assigned to",
      "a treatment group; behavioral experimenters were blinded to",
      "the drug treatment."
    ),
    dose_range      = paste(
      "AZD6088 (MW 406.57 g/mol) administered orally as vehicle or at",
      "1, 2.5, 5, 10, 20, or 40 umol/kg (0.407-16.26 mg/kg). Initial",
      "efficacy study: n = 7-12/group at 1, 2.5, 5, 10, 20, 40",
      "umol/kg. Pilot study: 2.5 umol/kg once daily for 8 days.",
      "Main study: twice-daily dosing (08:00 and 16:00) for 8 days,",
      "PD tested on days 1, 3, 5, 8 at 1, 2, 4, 6, 7, 24 h after the",
      "first daily administration; blood sampled at 2, 4, 6, 7 h",
      "after administration on days 2 and 4 in the same animals",
      "(PD-and-PK-in-one-animal design); 3 satellite rats/group",
      "sampled at 1, 2, 4, 6, 7, 24 h on day 1 only."
    ),
    regions         = "Canada (AstraZeneca R&D Montreal; single laboratory).",
    notes           = paste(
      "PK samples below the LC-MS/MS LLOQ of 0.0005 umol/L were",
      "excluded from the analysis (five PK samples in the main study,",
      "per Results). Total sample size across studies inferred from",
      "the design description: initial-efficacy SNL cohorts 1-40",
      "umol/kg (~50 rats), pilot 2.5 umol/kg n=9 + vehicle n=9 (18),",
      "main study 3 x drug PD (n=6/group) + 3 x drug satellite PK",
      "(n=3/group) + 2 x SNL vehicle PD (n=6/group) + 1 x naive",
      "PD (n=6) (~42 rats in the main study alone; the aggregate is",
      "~110 rats across the three studies but only the main-study PK",
      "and PD data feed the reported PKPD estimates in Table I).",
      "The n_subjects value 42 records the main-study PK/PD cohort",
      "used for the final model. Body weight, satellite-vs-PD-",
      "measurement status, and time-since-surgery were tested as PK",
      "covariates but only the DOSE_HIGH (40 umol/kg) categorical",
      "was retained. Housed 6/cage, 22 +- 1.5 C, 30-80% RH, 12-h",
      "light/dark, at least 3 days acclimatisation, food and water",
      "ad libitum, CCAC and AAALAC accreditation."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters. Viberg 2012 Table I (page 4) reports
    # apparent-oral point estimates with 95% non-parametric bootstrap
    # confidence intervals. The paper uses per-kg-body-weight units
    # (L/h/kg for CL/F, L/kg for V/F); ka is not estimated but fixed
    # to 5 x kel = 5 x cl/vc (Results, PK-model paragraph): 'ka was
    # therefore fixed to 5 time clearance. Changing the ka to
    # different values did not have an effect on OFV or other PK
    # parameter estimates.' The literal reading '5 x clearance' is
    # dimensionally inconsistent (CL has units L/h/kg while ka has
    # units 1/h); the intended pharmacometric idiom is 'ka fixed 5x
    # greater than the elimination rate constant kel = cl/vc' to
    # prevent flip-flop and to lock absorption on the fast branch.
    # Both apparent oral clearance values are wrapped in fixed() -
    # style single-value log() ini() expressions because the paper's
    # estimator uses NONMEM VII with likelihood-based OFV, not a
    # penalised prior; the values here reproduce the paper's point
    # estimates exactly.
    # ------------------------------------------------------------------
    lcl <- log(10.87)
    label("Apparent oral clearance CL/F at standard dose cohorts (L/h/kg)")
    # Viberg 2012 Table I: Oral clearance = 10.87 L/h/kg, 95% CI (9.06, 13.11).

    lcl_highdose <- log(3.92)
    label("Apparent oral clearance CL/F at the 40 umol/kg highest-dose cohort (L/h/kg)")
    # Viberg 2012 Table I: Oral clearance high dose = 3.92 L/h/kg, 95% CI (2.86, 6.02).
    # Gated inside model() by DOSE_HIGH = 1.

    lvc <- log(35.3)
    label("Apparent oral volume of distribution V/F (L/kg)")
    # Viberg 2012 Table I: Oral volume = 35.3 L/kg, 95% CI (28.5, 43.8).

    # ------------------------------------------------------------------
    # PD parameters (Emax direct-effect model on paw withdrawal
    # latency, seconds). Viberg 2012 Discussion first paragraph:
    # 'The pharmacodynamic response was adequately described using a
    # direct effect with an Emax model.' Alternatives tested and
    # rejected: time-delay effect compartment (no OFV improvement),
    # sigmoidicity coefficient (no OFV improvement), tolerance
    # (EC50 or Emax changing over time; no OFV improvement).
    # ------------------------------------------------------------------
    lbaseline <- log(5.947)
    label("Baseline paw withdrawal latency in SNL rats at zero AZD6088 concentration (s)")
    # Viberg 2012 Table I: Baseline latency = 5.947 s, 95% CI (5.87, 6.01).
    # Consistent with the observed SNL-rat pre-treatment baseline
    # of 5.87 +- 0.2 s (paper Methods, SNL-model assessment).

    lemax <- log(10.09)
    label("Maximum drug-induced increase in paw withdrawal latency above baseline (s)")
    # Viberg 2012 Table I: Emax = 10.09 s, 95% CI (9.44, 11.4).
    # Maximum reachable latency at infinite AZD6088 concentration =
    # baseline + Emax = 5.947 + 10.09 = 16.037 s (below the 20 s
    # assay cutoff imposed to avoid thermal injury; paper Methods,
    # SNL-model assessment).

    lec50 <- log(0.0433)
    label("AZD6088 plasma concentration producing half of Emax (umol/L)")
    # Viberg 2012 Table I: EC50 = 0.0433 umol/L, 95% CI (0.01, 0.10).
    # An improvement over the initial-efficacy-study estimate of
    # 46.6 nM with 95% CI (-14, 107) nM (paper Results, main-study
    # paragraph).

    # ------------------------------------------------------------------
    # Residual error. Viberg 2012 Table I:
    #   Proportional PK residual (linear-space fraction) = 0.325 with
    #     95% CI (0.270, 0.375). Paper Methods, Data Analysis: 'as
    #     well as proportional or additive error models were also
    #     tested' - a proportional PK error model was selected.
    #   Additive PD residual (paw withdrawal latency, seconds) = 1.027
    #     with 95% CI (0.93, 1.12). Latency is on an absolute-seconds
    #     scale (assay range 0-20 s), so an additive error is the
    #     natural choice; the paper does not report a proportional-
    #     PD alternative.
    # ------------------------------------------------------------------
    propSd <- 0.325
    label("Proportional residual error on plasma AZD6088 concentration (fraction)")
    # Viberg 2012 Table I: Proportional error = 0.325, 95% CI (0.270, 0.375).

    addSd_latency <- 1.027
    label("Additive residual error on paw withdrawal latency (s)")
    # Viberg 2012 Table I: Additive error = 1.027 s, 95% CI (0.93, 1.12).
  })

  model({
    # ------------------------------------------------------------------
    # 1. Typical-value individual PK parameters. Per the operator-
    #    approved sidecar response (frompeople-633 request-002 q1 = B),
    #    the model is encoded as typical-value only because the
    #    source publication does not report the IIV magnitudes it
    #    claims to have estimated. Downstream users wanting
    #    stochastic simulations can add multiplicative log-normal
    #    IIV on cl and vc at a documented CV%; see the vignette
    #    Assumptions and deviations section.
    #
    #    DOSE_HIGH is a per-rat binary indicator (1 at 40 umol/kg,
    #    else 0); the ifelse() selects between the two published
    #    typical-value apparent oral clearances (Table I).
    # ------------------------------------------------------------------
    cl <- exp(lcl + (lcl_highdose - lcl) * DOSE_HIGH)
    vc <- exp(lvc)

    # ka fixed at 5 x kel (= 5 x cl/vc). Paper Results:
    # 'ka was therefore fixed to 5 time clearance'; standard
    # pharmacometric idiom for 'absorption locked to the fast branch,
    # 5x faster than elimination'.
    kel <- cl / vc
    ka  <- 5 * kel

    # ------------------------------------------------------------------
    # 2. One-compartment first-order oral PK. Doses are entered per
    #    kg body weight (umol/kg) into depot; central then holds
    #    umol/kg and central / vc yields umol/L directly (no MW
    #    conversion needed because dose and concentration are on
    #    the same molar scale). This preserves the paper's native
    #    reporting units and lets the EC50 = 0.0433 umol/L compare
    #    directly against Cc without conversion.
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ------------------------------------------------------------------
    # 3. Observation variables.
    #
    #    Cc: plasma concentration of AZD6088 in umol/L (paper LC-MS/MS
    #    reporting unit; LLOQ 0.0005 umol/L).
    #
    #    latency: paw withdrawal latency in seconds, direct-effect
    #    Emax on Cc:
    #        latency = baseline + Emax * Cc / (EC50 + Cc)
    #    at Cc = 0 the latency is baseline (5.947 s, SNL pre-treatment);
    #    at Cc -> Inf the latency asymptotes to baseline + Emax
    #    (5.947 + 10.09 = 16.04 s), below the 20 s thermal-safety
    #    cutoff.
    # ------------------------------------------------------------------
    Cc <- central / vc

    baseline <- exp(lbaseline)
    emax     <- exp(lemax)
    ec50     <- exp(lec50)

    latency <- baseline + emax * Cc / (ec50 + Cc)

    Cc      ~ prop(propSd)
    latency ~ add(addSd_latency)
  })
}
