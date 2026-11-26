Mohamed_2020_upadacitinib_clinrem <- function() {
  description <- paste(
    "Continuous-time Markov exposure-response model for clinical remission",
    "2.8/1.0 (very soft/liquid stool frequency <= 2.8 and abdominal pain score",
    "<= 1.0, neither worse than baseline) in adults with moderately-to-severely",
    "active Crohn's disease during the 16-week CELEST induction period on",
    "upadacitinib (Mohamed 2020). Three states (no remission, remission,",
    "absorbing dropout) with first-order transitions. The individual",
    "steady-state minimum upadacitinib plasma concentration over the dosing",
    "interval (CTROUGH, constant per patient) acts on the no-remission ->",
    "remission rate as 1 + Emax * C / (EC50 + C), and 1.57 is added to that",
    "factor during the first two weeks."
  )
  reference <- paste(
    "Mohamed M-EF, Klunder B, Lacerda AP, Othman AA. Exposure-Response",
    "Analyses for Upadacitinib Efficacy and Safety in the Crohn's Disease",
    "CELEST Study and Bridging to the Extended-Release Formulation. Clin",
    "Pharmacol Ther. 2020;107(3):639-649. doi:10.1002/cpt.1668"
  )
  vignette <- "Mohamed_2020_upadacitinib"

  # Markov state-occupancy probabilities (s_<state> namespace), not drug
  # compartments.
  paper_specific_compartments <- c("s_noresp", "s_resp", "s_dropout")
  units <- list(
    time = "h",
    dosing = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "ng/mL (CTROUGH covariate); prob_clinrem is a probability (0-1)"
  )

  compartmentData <- list(
    s_noresp = list(
      analyte = "probability of not being in clinical remission 2.8/1.0 and on study",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_resp = list(
      analyte = "probability of being in clinical remission 2.8/1.0 and on study",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_dropout = list(
      analyte = "probability of having dropped out",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Individual steady-state minimum (trough) upadacitinib plasma concentration over the dosing interval.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "TOTAL plasma concentration. In the source, the steady-state Cmin of each CELEST patient predicted from the empirical-Bayes PK parameters of the population PK model (Supplement Methods), held constant over the whole induction period from the first dose; 0 for placebo. For a new regimen derive it from Mohamed_2020_upadacitinib, e.g. as the minimum of a steady-state dosing interval. Table 2 gives the medians: 1.2, 2.6, 4.9 and 7.8 ng/mL for 3, 6, 12 and 24 mg IR twice daily, 1.7 ng/mL for 24 mg IR once daily.",
      source_name = "CMIN"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 216L,
    n_studies = 1L,
    age_range = "mean 40.8 years (SD 12.9)",
    weight_range = "mean 75.7 kg (SD 20.1)",
    sex_female_pct = 57,
    race_ethnicity = c(White = 86, Black = 6, Hispanic = 7, Asian = 1, Other = 0.5),
    disease_state = "Moderately-to-severely active Crohn's disease (CDAI 220-450, SES-CD >= 6 or >= 4 for ileal-only disease) with inadequate response or intolerance to an immunomodulator or anti-TNF therapy; 96% anti-TNF experienced. Baseline CDAI 300 (SD 59.9), SES-CD 14.9 (SD 8.4).",
    dose_range = "Placebo or immediate-release upadacitinib 3, 6, 12 or 24 mg twice daily or 24 mg once daily for 16 weeks (induction).",
    regions = "Multinational",
    notes = "CELEST (NCT02365649) induction period, Table 1. Remission assessed among patients with baseline stool frequency >= 4.0 or abdominal pain score >= 2.0, at baseline and weeks 2, 4, 8, 12 and 16. NONMEM, Laplace; the full control stream is printed for this model (Table S7)."
  )

  ini({
    # ---- Markov transition rates (Table S5, clinical remission 2.8/1.0, 1/h) --
    # Placebo K01 and K10 were fixed from the placebo-only fit (Supplement
    # Methods). Table S5 footnote b: every rate was estimated as exp(theta);
    # the control stream writes 10**THETA.
    llambda01 <- fixed(log(0.0000603))
    label("Placebo rate no remission -> remission (1/h)") # Table S5 'K01 [1/hour]' clinical remission = 0.0000603 fix
    llambda10 <- fixed(log(0.000324))
    label("Placebo rate remission -> no remission (1/h)") # Table S5 'K10 [1/hour]' clinical remission = 0.000324 fix
    llambda0d <- log(0.0000813)
    label("Rate no remission -> dropout (1/h)") # Table S5 'K0D [1/hour]' clinical remission = 0.0000813
    llambda1d <- log(0.0000158)
    label("Rate remission -> dropout (1/h)") # Table S5 'K1D [1/hour]' clinical remission = 0.0000158
    lscale_lambda <- log(1.80)
    label("Multiplicative factor on both K01 and K10 at all times (unitless)") # Table S5 'Multiplicative factor for K01 and K10 for all times' = 1.80 (footnote f)

    # ---- Drug effect on K01 (Supplement Methods Equation 1) -------------------
    lemax <- log(4.57)
    label("Maximum fractional increase of K01 by upadacitinib (unitless)") # Table S5 'Emax' clinical remission = 4.57
    lec50 <- log(6.11)
    label("Steady-state trough upadacitinib concentration giving half of Emax (ng/mL)") # Table S5 'EC50 (ng/mL)' clinical remission = 6.11
    lhill <- fixed(log(1))
    label("Hill coefficient of the concentration-effect relationship (unitless)") # Table S7 control stream 'HILL=THETA(9)' with THETA(9) = 1 FIX
    e_early_lambda01 <- 1.57
    label("Constant added to the K01 drug-effect factor during the first two weeks (unitless)") # Table S5 'Increase parameter for K01 in first 2 weeks' = 1.57; footnote d: (1 + Emax model + theta2) x K01

    # Placeholder residual: the Markov likelihood is categorical (observed state
    # at each visit) and has no residual error.
    addSd_prob_clinrem <- fixed(0.001)
    label("Placeholder additive residual SD on prob_clinrem (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # ---- Drug effect ------------------------------------------------------------
    # Table S7 control stream: EFF = 1 + EMAX*(CMIN+minc)**HILL /
    # ((IC50+minc)**HILL + (CMIN+minc)**HILL), minc = 0.001 ng/mL;
    # pboresp = tshift for T <= 336 h, else 0; K01 = 10**TK01 * (EFF + pboresp)
    # * pref; K10 = 10**TK10 * pref. Time 0 is the first dose.
    minc <- 0.001
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    eff <- 1 + emax * (CTROUGH + minc)^hill / ((ec50 + minc)^hill + (CTROUGH + minc)^hill)

    early <- 0
    if (t <= 336) early <- 1
    scale_lambda <- exp(lscale_lambda)
    lambda01 <- exp(llambda01) * (eff + e_early_lambda01 * early) * scale_lambda
    lambda10 <- exp(llambda10) * scale_lambda
    lambda0d <- exp(llambda0d)
    lambda1d <- exp(llambda1d)

    # ---- State-occupancy probabilities (Figure S3) -------------------------------
    # The control stream integrates the transition probabilities from the
    # previous observed state; integrating the occupancy probabilities from
    # 'everyone starts out of remission' gives the population curves the paper
    # simulates. Dropout is absorbing.
    d/dt(s_noresp) <- -(lambda01 + lambda0d) * s_noresp + lambda10 * s_resp
    d/dt(s_resp) <- lambda01 * s_noresp - (lambda10 + lambda1d) * s_resp
    d/dt(s_dropout) <- lambda0d * s_noresp + lambda1d * s_resp
    s_noresp(0) <- 1

    # Remitters who have not dropped out (nonresponder imputation).
    prob_clinrem <- s_resp
    prob_clinrem ~ add(addSd_prob_clinrem)
  })
}
