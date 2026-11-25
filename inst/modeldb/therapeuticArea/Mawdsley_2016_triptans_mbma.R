Mawdsley_2016_triptans_mbma <- function() {
  description <- paste0(
    "MBMA (model-based network meta-analysis, MBNMA). Landmark logistic Emax ",
    "dose-response of the proportion of migraine patients headache-free 2 h ",
    "after a single dose of one of seven oral triptans (eletriptan, ",
    "sumatriptan, frovatriptan, almotriptan, zolmitriptan, naratriptan, ",
    "rizatriptan) or placebo, fitted in a Bayesian network meta-analysis to ",
    "182 arms of 70 randomized controlled trials. On the logit scale an arm's ",
    "response is the placebo log-odds of a typical trial plus an agent-specific ",
    "Emax function of the dose standardized to that agent's common (US-label) ",
    "single dose: logit(p) = E0 + Emax_k * x / (ED50_k + x), x = dose / ",
    "common dose. The final model ('Emax and ED50 exchangeable') treats each ",
    "agent's Emax and log ED50 as exchangeable within the triptan class. ",
    "Variability is BETWEEN-TRIAL only - a random placebo intercept and a ",
    "random relative treatment effect - so the simulation scope is the ",
    "trial-arm response rate, NOT an individual patient's response. The ",
    "per-agent Emax and ED50 values are not printed in the source; they were ",
    "digitized from the Emax-MBNMA prediction curves of Figure 4 (see the ",
    "vignette), and are consistent with the printed class means of Table 2."
  )

  reference <- paste(
    "Mawdsley D, Bennetts M, Dias S, Boucher M, Welton NJ.",
    "Model-Based Network Meta-Analysis: A Framework for Evidence Synthesis",
    "of Clinical Trial Data.",
    "CPT Pharmacometrics Syst Pharmacol. 2016 Aug;5(8):393-401.",
    "doi:10.1002/psp4.12091.",
    "Structural model: Eqs. 6-10 and the Supplementary Material JAGS code",
    "(EmaxClassED50.jags). Placebo response of a typical trial: constants",
    "m.Pbo and sd.Pbo in the same JAGS code. Between-study SD: Table 1, row",
    "'Emax (Emax & ED50 exch.)'. Per-agent Emax and ED50: digitized from the",
    "Emax curves of Figure 4.",
    sep = " "
  )

  vignette <- "Mawdsley_2016_triptans_mbma"

  # The random relative treatment effect (Eq. 8) is a BETWEEN-STUDY effect
  # added to the logit-scale treatment contrast, not a variability on a named
  # fixed-effect parameter, so it cannot follow the eta<parameter> convention
  # (Takechi_2025_nemolizumab_mbma_iga precedent).
  paper_specific_etas <- c("eta_study_reff")

  units <- list(
    time = paste0(
      "not applicable (LANDMARK meta-analysis: every modelled observation is ",
      "the proportion of an arm headache-free 2 h after a single dose, so the ",
      "model has no time term and its prediction is constant in time)"
    ),
    dosing = paste0(
      "mg (single oral dose of the arm's triptan, supplied through the ",
      "DOSE_<TRIPTAN>_MG covariate columns and NOT as rxode2 dose events; ",
      "each dose is divided inside model() by the agent's common dose)"
    ),
    concentration = paste0(
      "fraction/arm (probability that a patient in the trial arm is ",
      "headache-free 2 h after dosing; prob_headache_free_2h is NOT a drug ",
      "concentration - the slash only satisfies checkModelConventions unit ",
      "parsing)"
    )
  )

  covariateData <- list(
    DOSE_ELETRIPTAN_MG = list(
      description = "Single oral eletriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving eletriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "eletriptan common dose of 40 mg. Doses in the analysed trials: 20, ",
        "40 and 80 mg (0.5x-2x common dose)."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 1 (TriptansReformatted.xlsx, sheet 'result')"
    ),
    DOSE_SUMATRIPTAN_MG = list(
      description = "Single oral sumatriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving sumatriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "sumatriptan common dose of 50 mg, which the source states explicitly ",
        "('following Thorlund et al., a common dose of 50 mg was used'). ",
        "Doses in the analysed trials: 25, 50, 85, 100 and 200 mg ",
        "(0.5x-4x common dose)."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 2 (TriptansReformatted.xlsx, sheet 'result')"
    ),
    DOSE_FROVATRIPTAN_MG = list(
      description = "Single oral frovatriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving frovatriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "frovatriptan common dose of 2.5 mg. Doses in the analysed trials: ",
        "2.5 and 5 mg (1x-2x common dose)."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 3 (TriptansReformatted.xlsx, sheet 'result')"
    ),
    DOSE_ALMOTRIPTAN_MG = list(
      description = "Single oral almotriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving almotriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "almotriptan common dose of 12.5 mg. Doses in the analysed trials: ",
        "6.25, 12.5 and 25 mg (0.5x-2x common dose)."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 4 (TriptansReformatted.xlsx, sheet 'result')"
    ),
    DOSE_ZOLMITRIPTAN_MG = list(
      description = "Single oral zolmitriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving zolmitriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "zolmitriptan common dose of 2.5 mg. Doses in the analysed trials: ",
        "1, 2.5, 5, 10 and 25 mg (0.4x-10x common dose) - the widest ",
        "standardized range of the seven agents, which is why the source ",
        "plots zolmitriptan a second time over its full dose range."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 5 (TriptansReformatted.xlsx, sheet 'result')"
    ),
    DOSE_NARATRIPTAN_MG = list(
      description = "Single oral naratriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving naratriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "naratriptan common dose of 2.5 mg. Doses in the analysed trials: ",
        "2.5 and 5 mg (1x-2x common dose)."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 6 (TriptansReformatted.xlsx, sheet 'result')"
    ),
    DOSE_RIZATRIPTAN_MG = list(
      description = "Single oral rizatriptan dose assigned to the trial arm, in mg; 0 for any arm not receiving rizatriptan.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "MBMA study-ARM-level covariate. Standardized inside model() by the ",
        "rizatriptan common dose of 10 mg. Doses in the analysed trials: ",
        "2.5, 5, 10 and 40 mg (0.25x-4x common dose)."
      ),
      source_name = "Dose1Txt..Dose4Txt with Drug code 7 (TriptansReformatted.xlsx, sheet 'result')"
    )
  )

  population <- list(
    species = "human",
    n_studies = 70L,
    n_arms = 182L,
    disease_state = paste0(
      "Adults with migraine (IHS criteria) treating a single migraine attack ",
      "with an oral triptan or placebo. Endpoint: headache-free at 2 hours."
    ),
    dose_range = paste0(
      "Placebo and single oral doses of eletriptan 20-80 mg, sumatriptan ",
      "25-200 mg, frovatriptan 2.5-5 mg, almotriptan 6.25-25 mg, ",
      "zolmitriptan 1-25 mg, naratriptan 2.5-5 mg and rizatriptan 2.5-40 mg"
    ),
    design = paste0(
      "Network of 70 randomized controlled trials (placebo-controlled and ",
      "head-to-head, two- and multi-arm) assembled by Thorlund et al. 2014. ",
      "Patients with at least one migraine attack, not lost to follow-up and ",
      "without protocol violation; modified-ITT data where per-protocol data ",
      "were unavailable."
    ),
    notes = paste0(
      "MBMA AT THE TRIAL-ARM LEVEL: each modelled data point is one arm's ",
      "count of headache-free patients out of the arm size (binomial ",
      "likelihood), and the model is intended for simulating arm-level ",
      "response rates - it is NOT suitable for individual-patient ",
      "simulation. Patient demographics are not modelled. Common doses (the single dose indicated on the US ",
      "label; 50 mg for sumatriptan) are those the dose axis of Figure 4 ",
      "implies for the mg doses in the deposited trial table."
    )
  )

  ini({
    # ========================================================================
    # STRUCTURE (Eqs. 6, 8, 10 and the JAGS code EmaxClassED50.jags):
    #   logit(p[i,k]) = mu[i] + delta[i,k]
    #   delta[i,k] ~ N(f(x[i,k], t[i,k]) - f(x[i,1], t[i,1]), sd^2)
    #   f(x, t)     = Emax[t] * x / (ED50[t] + x),  x = dose / common dose
    # mu[i] is a per-trial nuisance parameter in the fit. Predictions for a
    # 'typical trial' (Figure 4) replace it by the mean of a separate
    # random-effects meta-analysis of the placebo arms (m.Pbo); that is the
    # form encoded here, with the placebo between-trial SD (sd.Pbo) as a
    # between-STUDY random intercept.
    # ========================================================================

    e0 <- fixed(-2.11)
    label("Placebo log-odds of being headache-free at 2 h in a typical trial (logit units); expit(-2.11) = 10.8%")
    # JAGS code EmaxClassED50.jags: 'm.Pbo <- -2.11', the mean of the
    # separate random-effects model for the placebo arms ('We obtain the
    # placebo effect estimate using a random effects model for the placebo
    # arms ... calculated separately to the MBNMA'). Fixed: a constant input
    # to the MBNMA predictions, not estimated by it. Independently confirmed:
    # every Emax curve of Figure 4 starts at -2.07 to -2.13 at dose 0 when its
    # intercept is left free in the digitization fit.

    # ---- Per-agent Emax (logit units) --------------------------------------
    # Not printed in the source. Each value is the least-squares fit of
    # -2.11 + Emax * x / (ED50 + x) to the red 'Emax' median prediction curve
    # of that agent's Figure 4 panel, digitized from the figure raster at
    # every pixel column (100-147 points per panel, residual RMSE 0.01-0.03
    # logit). Class check, Table 2 row 'Emax and ED50 exch.': mu_Emax =
    # 2.114 (1.446-2.838); the seven digitized values average 2.05.
    emax_eletriptan <- 2.52
    label("Maximal eletriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, eletriptan panel
    emax_sumatriptan <- 1.93
    label("Maximal sumatriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, sumatriptan panel
    emax_frovatriptan <- 2.06
    label("Maximal frovatriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, frovatriptan panel
    emax_almotriptan <- 1.90
    label("Maximal almotriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, almotriptan panel
    emax_zolmitriptan <- 2.12
    label("Maximal zolmitriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, 'zolmitriptan (full dose range)' panel (0-10x); the 0-2.2x panel gives 2.10
    emax_naratriptan <- 1.41
    label("Maximal naratriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, naratriptan panel
    emax_rizatriptan <- 2.42
    label("Maximal rizatriptan effect on the log-odds of headache-free at 2 h (logit units)") # digitized from Figure 4, rizatriptan panel

    # ---- Per-agent ED50 (log of multiples of the common dose) ---------------
    # The source estimates ED50 on the log scale ('Because ED50 is a dose, it
    # is necessarily positive, so we model this parameter on the log-scale';
    # JAGS 'exp(dED50[t])'). Same digitization as Emax. Class check, Table 2:
    # mu_ED50 = 0.676 (0.315-1.552) on the natural scale; the geometric mean
    # of the seven digitized values is 0.66.
    led50_eletriptan <- log(0.50)
    label("Log eletriptan dose giving half the maximal effect (log multiples of the 40 mg common dose)") # digitized from Figure 4, eletriptan panel
    led50_sumatriptan <- log(0.70)
    label("Log sumatriptan dose giving half the maximal effect (log multiples of the 50 mg common dose)") # digitized from Figure 4, sumatriptan panel
    led50_frovatriptan <- log(0.68)
    label("Log frovatriptan dose giving half the maximal effect (log multiples of the 2.5 mg common dose)") # digitized from Figure 4, frovatriptan panel
    led50_almotriptan <- log(0.78)
    label("Log almotriptan dose giving half the maximal effect (log multiples of the 12.5 mg common dose)") # digitized from Figure 4, almotriptan panel
    led50_zolmitriptan <- log(0.72)
    label("Log zolmitriptan dose giving half the maximal effect (log multiples of the 2.5 mg common dose)") # digitized from Figure 4, 'zolmitriptan (full dose range)' panel; the 0-2.2x panel gives 0.68
    led50_naratriptan <- log(1.01)
    label("Log naratriptan dose giving half the maximal effect (log multiples of the 2.5 mg common dose)") # digitized from Figure 4, naratriptan panel
    led50_rizatriptan <- log(0.53)
    label("Log rizatriptan dose giving half the maximal effect (log multiples of the 10 mg common dose)") # digitized from Figure 4, rizatriptan panel

    # ---- Between-STUDY variability ------------------------------------------
    # Placebo log-odds varies between trials. JAGS code: 'sd.Pbo <- 0.78'
    # with 'E0Pred ~ dnorm(m.Pbo, p.Pbo)' as the predictive distribution for
    # a new trial. ini() takes the variance: 0.78^2 = 0.6084.
    eta_study_e0 ~ fixed(0.6084)

    # Random relative treatment effect (Eq. 8): delta ~ N(f_k - f_1, sd^2).
    # Table 1, row 'Emax (Emax & ED50 exch.)', between-study SD 0.242
    # (0.16-0.335), posterior median. Variance 0.242^2 = 0.058564.
    eta_study_reff ~ 0.058564

    # ---- Residual error ------------------------------------------------------
    addSd_prob_headache_free_2h <- fixed(0.001)
    label("Placeholder additive residual SD on the arm headache-free probability; the source likelihood is binomial on the arm count, so there is no source residual")
    # NOT FROM SOURCE. The source likelihood is 'resp[i,k] ~ dbin(p[i,k],
    # npt[i,k])', so arm-level sampling error follows from the arm size. Draw
    # rbinom(1, N_arm, prob_headache_free_2h) / N_arm for an observed rate.
  })

  model({
    # ---- Dose standardized to each agent's common dose -----------------------
    # Common dose = the single dose indicated on the US label (Methods,
    # 'Illustrative example'); sumatriptan uses 50 mg following Thorlund 2014.
    x_eletriptan <- DOSE_ELETRIPTAN_MG / 40
    x_sumatriptan <- DOSE_SUMATRIPTAN_MG / 50
    x_frovatriptan <- DOSE_FROVATRIPTAN_MG / 2.5
    x_almotriptan <- DOSE_ALMOTRIPTAN_MG / 12.5
    x_zolmitriptan <- DOSE_ZOLMITRIPTAN_MG / 2.5
    x_naratriptan <- DOSE_NARATRIPTAN_MG / 2.5
    x_rizatriptan <- DOSE_RIZATRIPTAN_MG / 10

    # ---- Agent-specific Emax dose-response (Eq. 10) --------------------------
    # One arm receives at most one triptan, so at most one term is nonzero;
    # placebo (all doses 0) gives f = 0 and recovers the placebo intercept.
    f_eletriptan <- emax_eletriptan * x_eletriptan / (exp(led50_eletriptan) + x_eletriptan)
    f_sumatriptan <- emax_sumatriptan * x_sumatriptan / (exp(led50_sumatriptan) + x_sumatriptan)
    f_frovatriptan <- emax_frovatriptan * x_frovatriptan / (exp(led50_frovatriptan) + x_frovatriptan)
    f_almotriptan <- emax_almotriptan * x_almotriptan / (exp(led50_almotriptan) + x_almotriptan)
    f_zolmitriptan <- emax_zolmitriptan * x_zolmitriptan / (exp(led50_zolmitriptan) + x_zolmitriptan)
    f_naratriptan <- emax_naratriptan * x_naratriptan / (exp(led50_naratriptan) + x_naratriptan)
    f_rizatriptan <- emax_rizatriptan * x_rizatriptan / (exp(led50_rizatriptan) + x_rizatriptan)
    f_triptan <- f_eletriptan + f_sumatriptan + f_frovatriptan + f_almotriptan +
      f_zolmitriptan + f_naratriptan + f_rizatriptan

    # ---- Random relative effect on active arms -------------------------------
    # In a placebo-controlled trial the active arms are exactly the arms that
    # carry a relative effect versus arm 1 (placebo).
    active <- 0
    if (x_eletriptan + x_sumatriptan + x_frovatriptan + x_almotriptan +
      x_zolmitriptan + x_naratriptan + x_rizatriptan > 0) {
      active <- 1
    }

    logit_headache_free_2h <- e0 + eta_study_e0 + f_triptan + eta_study_reff * active

    # ---- Arm probability of being headache-free at 2 h -----------------------
    prob_headache_free_2h <- expit(logit_headache_free_2h)
    prob_headache_free_2h ~ add(addSd_prob_headache_free_2h)
  })
}
