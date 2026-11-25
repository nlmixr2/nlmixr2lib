Jakobsson_2026_velsecorat <- function() {
  description <- paste(
    "Two-state mixed-effects hidden Markov model (MEHMM) of twice-daily",
    "home-measured peak expiratory flow (PEF, L/min) in adults with asthma",
    "treated for 12 weeks with once-daily inhaled velsecorat (AZD7594) or",
    "placebo in the phase 2b dose-finding trial NCT03622112 (Jakobsson 2026,",
    "final 'mixed' dose-response model of supplement Table S4). A latent",
    "discrete-time Markov chain switches each subject between a high-PEF",
    "state (mean m0) and a low-PEF state (mean d * m0); each PEF measurement",
    "is normal around the current state's mean with a subject-specific",
    "within-state SD. The high-to-low and low-to-high transition",
    "probabilities apply per MEASUREMENT (every record is one twice-daily",
    "PEF reading) and are capped at 0.5 by a logit(2p) transform. Dose (ug)",
    "enters as a covariate: a constant shift on logit(d) for any active",
    "dose, and Emax dose-response terms on log(sigma^2) and on the",
    "low-to-high transition. SIMULATION-ONLY: the latent state is drawn with",
    "rxunif() and carried record to record with lag(), so the model",
    "simulates PEF series and latent paths but cannot be used to re-estimate",
    "parameters (its nlmixr2 likelihood would not marginalise the hidden",
    "state as the paper's forward recursion does)."
  )

  reference <- paste(
    "Jakobsson L, Baaz M, Leander J, Gerlee P, Jirstrand M.",
    "Modeling Peak Expiratory Flow in Patients With Asthma and Quantifying",
    "Treatment Effects Using a Mixed-Effects Hidden Markov Model.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70281.",
    "doi:10.1002/psp4.70281. PMID 42249718; PMCID PMC13245144.",
    "Authors' code: github.com/fcc-ludvig-jakobsson/pef-mehmm (commit 0297f92).",
    sep = " "
  )

  vignette <- "Jakobsson_2026_velsecorat"

  units <- list(
    time = "day (PEF is read twice daily; the model has no explicit time dependence and advances the latent chain once per record)",
    dosing = "ug (once-daily nominal inhaled velsecorat dose supplied as the DOSE_VELSECORAT_UG covariate; the model consumes no rxode2 dose events)",
    concentration = "L/min (peak expiratory flow, a flow rate rather than a drug concentration)"
  )

  covariateData <- list(
    DOSE_VELSECORAT_UG = list(
      description = "Nominal once-daily inhaled velsecorat dose of the randomised arm",
      units = "ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The dose x_i of Jakobsson 2026 Section 2.3 ('the patients are",
        "assumed to be given either a placebo drug or an active drug at some",
        "dose x_i'); 0 for placebo. Trial arms 50, 90, 180, 360 and 720 ug",
        "(supplement Table S3 categorical-effect subscripts). Constant within",
        "a subject: the model applies the dose effect instantaneously from",
        "randomisation. Enters as an indicator (DOSE_VELSECORAT_UG > 0) on",
        "logit(d) and as an Emax term DOSE / (ED50 + DOSE) on log(sigma^2)",
        "and on logit(2 p10)."
      ),
      source_name = "x_i (dose; Jakobsson 2026 Section 2.3, dose-response equations)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "Asthma, symptomatic on low-dose inhaled corticosteroid (trial title of NCT03622112)",
    dose_range = "Placebo or velsecorat (AZD7594) DPI 50, 90, 180, 360 or 720 ug once daily for 12 weeks",
    regions = "Multi-centre (not detailed in the source)",
    notes = paste(
      "Patient-level data from the AstraZeneca phase 2b dose-finding trial",
      "NCT03622112 (Jakobsson 2026 Section 2.5, reference 23): 5 active arms",
      "and placebo, twice-daily home-measured PEF over a 12-week treatment",
      "period. The paper does not report the number of patients analysed,",
      "the number of PEF records, or baseline demographics. It notes",
      "substantial dropout in the placebo arm (Discussion)."
    )
  )

  ini({
    # Typical values (supplement Table S4, final 'mixed' dose-response model;
    # transformed scale phi = mu + beta(x) + eta, Jakobsson 2026 Section 2.3)
    lpef_high <- 5.755; label("Log high-state mean PEF, log(m0) (L/min)") # Table S4 mu1 = 5.755
    logitd_low <- 1.984; label("Logit of the low/high PEF ratio d at placebo (unitless)") # Table S4 mu2 = 1.984
    lvar_pef <- 6.837; label("Log within-state PEF variance at placebo, log(sigma^2) ((L/min)^2)") # Table S4 mu3 = 6.837
    logit2p_hl <- -2.716; label("logit(2 p01): high-to-low transition probability per measurement (unitless)") # Table S4 mu4 = -2.716
    logit2p_lh <- -3.008; label("logit(2 p10): low-to-high transition probability per measurement at placebo (unitless)") # Table S4 mu5 = -3.008

    # Dose-response terms (Table S4; forms from Section 2.3 Eqs. 'Constant' and 'Emax')
    e_active_logitd <- 0.212; label("Constant shift on logit(d) for any active velsecorat dose (unitless)") # Table S4 beta2 = 0.212
    emax_lvar <- -0.447; label("Emax of the dose effect on log(sigma^2) (unitless)") # Table S4 Emax,3 = -0.447
    led50_lvar <- 5.039; label("Log ED50 of the dose effect on log(sigma^2) (log ug)") # Table S4 log(ED50,3) = 5.039
    emax_logit2p_lh <- 0.705; label("Emax of the dose effect on logit(2 p10) (unitless)") # Table S4 Emax,5 = 0.705
    led50_logit2p_lh <- 3.811; label("Log ED50 of the dose effect on logit(2 p10) (log ug)") # Table S4 log(ED50,5) = 3.811

    # IIV: diagonal Omega, variances on the transformed scale (Table S4)
    etalpef_high ~ 0.129 # Table S4 omega1^2 = 0.129
    etalogitd_low ~ 0.390 # Table S4 omega2^2 = 0.390
    etalvar_pef ~ 0.681 # Table S4 omega3^2 = 0.681
    etalogit2p_hl ~ 0.839 # Table S4 omega4^2 = 0.839
    etalogit2p_lh ~ 1.038 # Table S4 omega5^2 = 1.038
  })

  model({
    # Dose-response terms (Section 2.3): constant effect for any active dose,
    # Emax with ED50 = exp(log ED50) in ug (authors' code simulate_two_state_model)
    active <- 0
    if (DOSE_VELSECORAT_UG > 0) active <- 1
    eff_lvar <- emax_lvar * DOSE_VELSECORAT_UG / (exp(led50_lvar) + DOSE_VELSECORAT_UG)
    eff_logit2p_lh <- emax_logit2p_lh * DOSE_VELSECORAT_UG / (exp(led50_logit2p_lh) + DOSE_VELSECORAT_UG)

    # Individual parameters, psi = inv_g(phi) with
    # inv_g = (exp, invlogit, exp, invlogit / 2, invlogit / 2) (authors' code)
    pef_high <- exp(lpef_high + etalpef_high)
    d_low <- expit(logitd_low + e_active_logitd * active + etalogitd_low)
    var_pef <- exp(lvar_pef + eff_lvar + etalvar_pef)
    sd_pef <- sqrt(var_pef)
    p_hl <- expit(logit2p_hl + etalogit2p_hl) / 2
    p_lh <- expit(logit2p_lh + eff_logit2p_lh + etalogit2p_lh) / 2
    pef_low <- d_low * pef_high

    # Latent two-state Markov chain (Section 2.3, P(Z_t = l | Z_t-1 = j) = p_jl). Every subject starts in the
    # high state (state 0) and one transition is drawn before each
    # measurement, so every record must be one PEF reading. lag() returns NA
    # on a subject's first record, which seeds the chain at state 0.
    u_trans <- rxunif(0, 1)
    state_prev <- lag(state_low)
    if (is.na(state_prev)) state_prev <- 0
    p_switch <- p_hl
    if (state_prev == 1) p_switch <- p_lh
    state_low <- state_prev
    if (u_trans < p_switch) state_low <- 1 - state_prev

    # Observation model (Section 2.3): Y_t = m_Z + eps, eps ~ N(0, sigma^2)
    PEF <- pef_high
    if (state_low == 1) PEF <- pef_low
    PEF ~ add(sd_pef)
  })
}
