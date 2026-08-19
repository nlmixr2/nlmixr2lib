Nguyen_2025_valbenazine_tmc <- function() {
  description <- paste(
    "Longitudinal exposure-efficacy model of the change from baseline in the",
    "UHDRS Total Maximal Chorea (TMC) score under valbenazine in adults with",
    "Huntington's-disease chorea (Nguyen 2025, KINECT-HD).",
    "The response is DECOUPLED: a time-dependent placebo term applies only to",
    "placebo-arm subjects, and a nonlinear Emax term in the daily",
    "[+]-alpha-dihydrotetrabenazine ([+]-alpha-HTBZ) AUC applies to everyone.",
    "The placebo term is an asymptotic exponential onset,",
    "1 - exp(-KPBO * t), with an onset half-life of 9.18 days, multiplied",
    "after day TATT = 59 by an attenuation factor exp(-KATT * (t - TATT))",
    "with an attenuation half-life of 36.4 days -- reproducing the observed",
    "placebo response that peaks near week 8 and then returns toward",
    "baseline. The drug term is Emax * AUC / (EC50 + AUC) with",
    "Emax = -14.2 TMC points and EC50 = exp(7.50) ng*h/mL.",
    "Inter-individual variability is ADDITIVE on Pmax and log-normal on",
    "EC50 (125% CV), and the residual error is additive on the TMC scale.",
    "The model has NO ODE states and no dosing: it is an algebraic function",
    "of study time, the placebo-arm indicator, and the exposure covariate",
    "AUC_HTBZ, which must be supplied per observation record -- typically",
    "from the companion population PK model,",
    "modellib('Nguyen_2025_valbenazine').",
    sep = " "
  )
  reference <- paste(
    "Nguyen HQ, Crass RL, Chapel S, Kuan HYS, Loewen G, Brar S.",
    "Population pharmacokinetic and exposure-efficacy analyses of",
    "valbenazine in patients with Huntington's disease: supporting dose",
    "selection for chorea management.",
    "The Journal of Clinical Pharmacology. 2025;65(12):1777-1788.",
    "doi:10.1002/jcph.70092.",
    "Parameter estimates from Table 3; placebo and drug-effect equations",
    "from the Methods 'Exposure-Efficacy Modeling' subsection.",
    "Underlying phase 3 study KINECT-HD, reported in",
    "Furr Stimming E, Claassen DO, Kayson E, et al. Lancet Neurology.",
    "2023;22(6):494-504. doi:10.1016/S1474-4422(23)00127-8.",
    "Companion population PK model from the same paper:",
    "modellib('Nguyen_2025_valbenazine').",
    sep = " "
  )
  vignette <- "Nguyen_2025_valbenazine_huntington_chorea"

  depends <- c("AUC_HTBZ")

  units <- list(
    time          = "day",
    dosing        = "N/A (PD-only; the [+]-alpha-HTBZ daily AUC is a required input covariate)",
    concentration = "TMC points (observation: change from baseline in the UHDRS Total Maximal Chorea score); ng*h/mL (AUC_HTBZ input covariate)"
  )

  covariateData <- list(
    AUC_HTBZ = list(
      description        = "Daily (24 h) steady-state area under the plasma concentration-time curve of the active valbenazine metabolite [+]-alpha-dihydrotetrabenazine ([+]-alpha-HTBZ, NBI-98782), driving the Emax drug effect on the TMC score.",
      units              = "ng*h/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Supplied per observation record; 0 for placebo-arm records. Nguyen 2025 Methods: 'drug exposure as measured by the daily AUC of [+]-alpha-HTBZ'. Results: AUC24 was chosen over Cmax, Cmin and Cavg because the E-R profiles were consistent across all four metrics. In the paper these values are Bayesian post-hoc predictions from the companion population PK model, modellib('Nguyen_2025_valbenazine'), computed over a 24 h steady-state interval; Table 2 gives the HD-cohort medians (5th; 95th percentile) as 176 (107; 328), 402 (244; 749), 657 (399; 1225) and 931 (566; 1736) ng*h/mL at 20, 40, 60 and 80 mg once daily. Range explored in the E-R simulations: 0-3000 ng*h/mL. The Discussion notes that most observed data fall in the linear part of the curve (at or below 1500 ng*h/mL), so extrapolation above ~2000 ng*h/mL leans on the Emax asymptote rather than on data.",
      source_name        = "AUC[+]-alpha-HTBZ"
    ),
    PLACEBO = list(
      description        = "Placebo-arm membership indicator gating the time-dependent placebo term. 1 = placebo, 0 = valbenazine.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (valbenazine active-treatment arm)",
      notes              = "Time-fixed per subject. Nguyen 2025 Methods defines the paper's variable PLB as 'the indicator variable for placebo treatment (1 = placebo; 0 = valbenazine)' and the response equation as R = PBO * PLB + Emax * AUC / (EC50 + AUC). The polarity matches the PLACEBO canonical (1 = placebo). Note the consequence of this DECOUPLED parameterisation: valbenazine-treated subjects carry NO placebo component at all, because the two effects showed different time courses in the exploratory analysis ('Given the differing patterns in TMC scores over time, the placebo and drug effects were modeled separately using decoupled models'). KINECT-HD randomised 1:1; the E-R analysis used 61 placebo and 64 valbenazine subjects.",
      source_name        = "PLB"
    )
  )

  # Covariates that Nguyen 2025 pre-specified and SCREENED in the E-R
  # covariate model but did not retain: none reached statistical
  # significance, and no point estimates are reported for them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age. Screened on Pmax, Emax and EC50.",
      units       = "years",
      type        = "continuous",
      notes       = "Nguyen 2025 Results: 'no statistically significant covariate-parameter relationships were identified among the factors tested (i.e., sex, age, baseline body weight, baseline CGI-S, and baseline PGI-S on Pmax, Emax, and EC50)'. Backward elimination threshold was dOFV > 10.8 (P < .001). KINECT-HD full analysis set: mean (SD) 53.7 (10.8) years, median 55.0 (range 25-74) -- Supplemental Table S6."
    ),
    SEXF = list(
      description = "Sex. Screened on Pmax, Emax and EC50.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Not retained (see AGE note). KINECT-HD full analysis set: 54.4% female (Supplemental Table S5)."
    ),
    WT = list(
      description = "Baseline body weight. Screened on Pmax, Emax and EC50.",
      units       = "kg",
      type        = "continuous",
      notes       = "Not retained (see AGE note). KINECT-HD full analysis set: mean (SD) 77.0 (18.0) kg, median 75.8 (range 41.0-134) -- Supplemental Table S6."
    ),
    CGI_S_BASE = list(
      description = "Baseline Clinical Global Impression of Severity (CGI-S), a 7-point clinician-rated chorea severity scale (1 = normal, not at all ill; 7 = among the most extremely ill patients). Screened on Pmax, Emax and EC50.",
      units       = "score (1-7)",
      type        = "continuous",
      notes       = "Not retained (see AGE note). No canonical register entry is required because the covariate never enters model(). KINECT-HD full analysis set dichotomised at 4: 51.2% below 4, 48.8% at or above 4 (Supplemental Table S5)."
    ),
    PGI_S_BASE = list(
      description = "Baseline Patient Global Impression of Severity (PGI-S), a 5-point patient-rated chorea severity scale (1 = none; 5 = very severe). Screened on Pmax, Emax and EC50.",
      units       = "score (1-5)",
      type        = "continuous",
      notes       = "Not retained (see AGE note). KINECT-HD full analysis set dichotomised at 3: 55.2% below 3, 44.8% at or above 3 (Supplemental Table S5)."
    ),
    TMC_BASE = list(
      description = "Baseline UHDRS Total Maximal Chorea score (0-28; the sum of seven body-region subscores each scored 0-4). Screened specifically on the maximum placebo effect Pmax and the maximum drug effect Emax.",
      units       = "TMC points (0-28)",
      type        = "continuous",
      notes       = "Not retained. Nguyen 2025 Results: 'The effects of baseline TMC on the maximum placebo and maximum drug effects were tested in the base E-R model. However, no improvement in model fit was observed, and the baseline covariate effects were poorly estimated, indicating that the extent of placebo and drug effects on TMC scores is not related to the baseline TMC score.' Because the endpoint is the CHANGE from baseline, the baseline score is not needed to evaluate the model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 125L,
    n_studies      = 1L,
    n_observations = 937L,
    age_range      = "25-74 years (mean (SD) 53.7 (10.8); median 55.0)",
    age_median     = "55.0 years",
    weight_range   = "41.0-134 kg (mean (SD) 77.0 (18.0); median 75.8)",
    weight_median  = "75.8 kg",
    sex_female_pct = 54.4,
    race_ethnicity = c(White = 96.0, Black = 0.8, Asian = 0.8, Other = 2.4),
    hispanic_pct   = 6.4,
    disease_state  = "Adults with chorea associated with Huntington's disease (KINECT-HD, protocol HD3005; a phase 3 randomized double-blind placebo-controlled trial)",
    dose_range     = "Placebo, or valbenazine titrated from 40 mg once daily in 20 mg increments at the end of weeks 2, 4 and 6 to a target of 80 mg once daily; individual maintenance doses of 20, 40, 60 or 80 mg once daily as tolerated",
    cyp2d6_status  = c(extensive = 48.0, intermediate = 37.6, poor = 7.2, ultrarapid = 4.8, missing = 2.4),
    notes          = paste(
      "128 subjects received study treatment (64 placebo, 64 valbenazine);",
      "three placebo subjects had no post-treatment TMC assessment, so the",
      "E-R analysis dataset comprised 61 placebo and 64 valbenazine",
      "subjects contributing 937 TMC assessments including 125 baseline",
      "assessments (Nguyen 2025 Results). Twelve weeks of treatment with a",
      "follow-up visit two weeks after the final dose; TMC assessed at",
      "screening, baseline and weeks 2, 4, 6, 8, 10 and 12 (Supplemental",
      "Table S1). By the maintenance period most subjects were at 80 mg",
      "(45 of 57 at week 12), leaving very small samples at 20 and 40 mg --",
      "the paper attributes the erratic observed means at those doses to",
      "this imbalance. Baseline characteristics from Supplemental Tables S5",
      "and S6.",
      sep = " "
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Placebo model -- Nguyen 2025 Methods, 'Exposure-Efficacy Modeling':
    #
    #   PBO = (Pmax + eta1_Pmax) * { 1 - exp(-KPBO * t)                          if t <= TATT
    #                              { (1 - exp(-KPBO * t)) * exp(-KATT*(t-TATT))  if t >  TATT
    #
    #   KPBO = 0.693 / HLPLB ;  KATT = 0.693 / HLATT
    #
    # Note the model uses the paper's literal 0.693 rather than log(2) =
    # 0.6931472, because the printed equation is the definition being
    # reproduced; the difference is 0.02% and immaterial.
    # ---------------------------------------------------------------------
    pmax <- -2.35; label("Maximum placebo effect on the change from baseline in TMC score (TMC points)")  # Table 3 theta 1: -2.35 (%RSE 14.8; 95% CI -3.04, -1.67). Kept on the linear scale because the value is negative (an improvement in chorea) and cannot be log-transformed.

    lthalf_onset <- log(9.18); label("Placebo-effect onset half-life HLPLB (day)")                        # Table 3 theta 2: 9.18 (%RSE 26.5; 95% CI 4.42, 13.9)
    ltswitch     <- log(59.0); label("Start time for attenuation of the placebo effect, TATT (day)")      # Table 3 theta 3: 59.0 (%RSE 4.29; 95% CI 54.0, 64.0). Results: 'the start time of placebo attenuation precisely estimated (RSE ~4%)'.
    lthalf_att   <- log(36.4); label("Placebo-effect attenuation half-life HLATT (day)")                  # Table 3 theta 4: 36.4 (%RSE 22.6; 95% CI 20.3, 52.6)

    # ---------------------------------------------------------------------
    # Drug effect -- Nguyen 2025 Methods:
    #
    #   R = PBO * PLB + Emax * AUC / (exp(log(EC50) + eta2_EC50) + AUC)
    #
    # The paper parameterises EC50 on the log scale with the eta INSIDE the
    # exponential, which is exactly the canonical lec50 + etalec50 form.
    # ---------------------------------------------------------------------
    emax  <- -14.2; label("Maximum drug effect on the change from baseline in TMC score (TMC points)")            # Table 3 theta 5: -14.2 (%RSE 10.9; 95% CI -17.3, -11.2). Linear scale: the value is negative.
    lec50 <-   7.50; label("Log of the daily [+]-alpha-HTBZ AUC producing 50% of Emax, log(EC50) (log ng*h/mL)")  # Table 3 theta 6: 7.50 (%RSE 2.92; 95% CI 7.08, 7.93). See the exp(7.50) note in model().

    # ---------------------------------------------------------------------
    # Inter-individual variability -- Table 3.
    # ---------------------------------------------------------------------
    etapmax  ~ 5.11   # Table 3 IIV Pmax: omega^2 = 5.11 (shrinkage 36.4%; 95% CI 2.64, 7.57) -> SD 2.26 TMC points = sqrt(5.11). ADDITIVE on Pmax per the printed equation (Pmax + eta1_Pmax), so this variance is on the TMC scale, not a log-scale variance.
    etalec50 ~ 0.941  # Table 3 IIV EC50: omega^2 = 0.941 (shrinkage 37.3%; 95% CI 0.442, 1.44) -> 125% CV = sqrt(exp(0.941)-1) = 1.250, matching the published 125% (74.6, 179) exactly and confirming the log scale.

    # ---------------------------------------------------------------------
    # Residual error -- additive on the TMC change-from-baseline scale.
    # ---------------------------------------------------------------------
    addSd <- 1.41; label("Additive residual SD on the change from baseline in TMC score (TMC points)")  # Table 3 theta 7: the reported estimate 1.98 is the VARIANCE and the transformed estimate 1.41 is its SD (sqrt(1.98) = 1.407). The reported 95% CI on the variance, (1.88, 2.07), back-transforms to (1.371, 1.439) and matches the printed SD CI (1.37, 1.44) exactly.
  })

  model({
    # ---------------------------------------------------------------------
    # Placebo time course. `time` is time since the first dose, in days
    # (units$time = "day").
    #
    # The piecewise definition is written with max(time - tswitch, 0) rather
    # than a branch: below TATT the exponent is 0 so the attenuation factor
    # is exactly 1, and above TATT it decays as exp(-KATT * (time - TATT)).
    # This is algebraically identical to the printed two-case form and avoids
    # a conditional inside the model block.
    # ---------------------------------------------------------------------
    thalf_onset <- exp(lthalf_onset)
    tswitch     <- exp(ltswitch)
    thalf_att   <- exp(lthalf_att)

    kpbo <- 0.693 / thalf_onset
    katt <- 0.693 / thalf_att

    onset <- 1 - exp(-kpbo * time)
    att   <- exp(-katt * max(time - tswitch, 0))
    pbo   <- (pmax + etapmax) * onset * att

    # ---------------------------------------------------------------------
    # Emax drug effect in the daily [+]-alpha-HTBZ AUC.
    #
    # exp(7.50) = 1808 ng*h/mL. Nguyen 2025 Table 3 prints the
    # back-transformed estimate as 1820 (1190, 2780); 1820 corresponds to
    # log(EC50) = 7.5066, so the printed 7.50 is a rounded THETA. The
    # reported CI back-transforms exactly (exp(7.08) = 1188, exp(7.93) =
    # 2781). This file encodes the reported THETA 7.50, which makes EC50
    # 0.7% lower than the printed 1820.
    # ---------------------------------------------------------------------
    ec50 <- exp(lec50 + etalec50)
    drug <- emax * AUC_HTBZ / (ec50 + AUC_HTBZ)

    # ---------------------------------------------------------------------
    # Decoupled response. The placebo term is gated by PLACEBO, so
    # valbenazine-treated subjects (PLACEBO = 0) carry the drug term only.
    # ---------------------------------------------------------------------
    tmccfb <- pbo * PLACEBO + drug

    tmccfb ~ add(addSd)
  })
}
