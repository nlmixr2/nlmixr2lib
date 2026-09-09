Parasrampuria_2025_temsavir_resp1log_day8 <- function() {
  description <- paste0(
    "Logistic-regression exposure-response model with an Emax linear ",
    "predictor, relating the steady-state trough plasma concentration of ",
    "temsavir (Ctau) to the probability of achieving a GREATER THAN ",
    "1.0 log10 c/mL decrease in plasma HIV-1 RNA from Day 1 to Day 8 ",
    "during fostemsavir functional monotherapy, in heavily ",
    "treatment-experienced adults with multidrug-resistant HIV-1 ",
    "(Parasrampuria 2025, n = 258 from the randomised cohort of the ",
    "phase 3 BRIGHTE study AI438047 / NCT02362503: 193 on fostemsavir ",
    "600 mg BID and 65 on placebo). The probability is ",
    "expit(E0 - Emax * Ctau/(EC50 + Ctau) + 0.848 * ",
    "(log10(HIV_VLOAD) - 4.65) - 0.781 * (CD4_ABS < 20)) with ",
    "E0 = -1.99, Emax = -2.67 and EC50 = 78.8 ng/mL (Parasrampuria 2025 ",
    "Table 4 and its Note). Emax is NEGATIVE and is SUBTRACTED, so the ",
    "exposure term raises the logit: higher trough gives a higher ",
    "response probability. This is the stricter of the paper's two ",
    "responder thresholds and the one carried into the Table 5 dosing ",
    "scenarios, where the simulated proportion of responders moves only ",
    "from 0.368 to 0.478 across coadministration with a moderate CYP3A ",
    "inducer, a strong CYP3A inhibitor, fasted versus fed dosing, and ",
    "body weights from 40 to 150 kg -- the quantitative basis for the ",
    "paper's conclusion that no dose adjustment is needed. EC50 and the ",
    "CD4+ effect are imprecise (RSE 81.4% and 45.4%), and the paper ",
    "flags subjects who contradict the relationship at both extremes. ",
    "There is no PK layer, no ODE and no between-subject random effect ",
    "(Bernoulli likelihood): the exposure metric is supplied as a data ",
    "column, derived post hoc from the companion population PK model ",
    "modellib('Parasrampuria_2025_temsavir'). The sibling model ",
    "Parasrampuria_2025_temsavir_resp05log_day8 gives the parallel fit ",
    "for the 0.5 log10 threshold."
  )
  reference <- paste(
    "Parasrampuria R, Thakkar N, Moore K, Ackerman P, Magee M.",
    "Population pharmacokinetics and exposure-response relationship for",
    "temsavir following fostemsavir administration in",
    "treatment-experienced HIV patients.",
    "Pharmacol Res Perspect. 2025;13(3):e70023.",
    "doi:10.1002/prp2.70023.",
    "Individual exposure metrics derive from the companion population",
    "pharmacokinetic model reported in the same paper; see",
    "modellib('Parasrampuria_2025_temsavir').",
    sep = " "
  )
  vignette <- "Parasrampuria_2025_temsavir"
  units <- list(
    time          = "n/a (static landmark exposure-response regression at Day 8; no time dimension)",
    dosing        = "n/a (no dose events; exposure enters as the CTROUGH covariate column)",
    concentration = "prob_hivrna_decr1log (probability of a > 1.0 log10 c/mL decrease in plasma HIV-1 RNA at Day 8, 0-1; also logit_hivrna_decr1log)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual steady-state plasma temsavir concentration at the end of the dosing interval (Ctau), per subject. Supplied as data: this model has no PK layer.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma temsavir at STEADY STATE on the",
        "approved fostemsavir 600 mg BID regimen (12 h dosing interval).",
        "Parasrampuria 2025 Methods 2.3: derived from post-hoc individual",
        "parameter estimates of the final population PK model and computed",
        "from the simulated profiles with the PKNCA R package, so these",
        "are empirical-Bayes model predictions rather than observed",
        "troughs. Reproduce them with",
        "modellib('Parasrampuria_2025_temsavir') solved at 600 mg BID with",
        "rate = -2 on the dose records. Table 5 gives the simulated",
        "medians this model was exercised at: 433 ng/mL dosed alone,",
        "205 with a moderate CYP3A inducer, 775 with a strong CYP3A",
        "inhibitor, 414 with both, 257 fasted, 599 at 40 kg and 296 at",
        "150 kg. Results 3.2.2 note the relationship is contradicted at",
        "both tails -- some subjects with Ctau below 10 ng/mL achieved",
        "the decline while some above 4000 ng/mL did not."
      ),
      source_name        = "Ctau"
    ),
    HIV_VLOAD = list(
      description        = "Baseline (Day 1) plasma HIV-1 RNA concentration.",
      units              = "copies/mL (c/mL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried LINEAR in c/mL as the <VIRUS>_VLOAD family requires, and",
        "converted to log10 and centred INSIDE model(). The Table 4 Note",
        "writes the term as theta1 * (BHIVRNA - 4.65); the centring",
        "constant 4.65 is the cohort median on the log10 scale, which",
        "Table S2 gives directly as 4.65 log10 c/mL (range 1.59-6.91,",
        "N = 258) and which equals log10(44940) = 4.65263. The companion",
        "continuous-endpoint model consumes the SAME column on the LINEAR",
        "scale as (BHIVRNA/44940)^theta1. The positive coefficient means",
        "a higher baseline viral load raises the probability of clearing",
        "the 1 log10 bar, consistent with the continuous model, where",
        "higher baseline load gives a larger decline."
      ),
      source_name        = "BHIVRNA"
    ),
    CD4_ABS = list(
      description        = "Baseline (Day 1) absolute peripheral-blood CD4+ T-lymphocyte count.",
      units              = "cells/mm^3",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried as the CONTINUOUS count and DICHOTOMISED INSIDE model()",
        "at 20 cells/mm3, following the register's HCV_VLOAD / ALP",
        "pattern; a pre-binarised dataset column would hide the threshold.",
        "The Table 4 abbreviations define BSL as the 'baseline CD4+ flag",
        "(>= 20 cells/mm3 or < 20 cells/mm3)'. BSL = 1 for CD4+ BELOW",
        "20 cells/mm3, so the negative coefficient LOWERS the response",
        "probability in the severely-immunosuppressed group. The",
        "coefficient here (-0.781) is roughly half the magnitude of the",
        "0.5 log10 model's (-1.27) and is the less precisely estimated of",
        "the two (RSE 45.4% vs 29.1%). Day 8 cohort (N = 258): 69",
        "subjects (27%) below 20 cells/mm3; median 98.5 cells/mm3, range",
        "0-1160 (Table S2)."
      ),
      source_name        = "BSL (derived from baseline CD4+ count)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 258L,
    n_studies      = 1L,
    n_observations = "258 binary Day 8 responder records (one per subject; landmark analysis, no repeated measures)",
    age_range      = "18-73 years (median 48) (Table S2)",
    weight_range   = "38-146 kg (median 70) (Table S2)",
    sex_female_pct = 26.0,
    race_ethnicity = c(White = 66.0, `Black or African American` = 23.0, Asian = 1.0, Other = 10.0),
    disease_state  = "heavily treatment-experienced (HTE) adults with multidrug-resistant HIV-1 infection failing their current antiretroviral regimen, from the randomised cohort of the phase 3 BRIGHTE study; median baseline plasma HIV-1 RNA 4.65 log10 c/mL (range 1.59-6.91), median baseline CD4+ count 98.5 cells/mm3 (range 0-1160) with 69 subjects (27%) below 20 cells/mm3",
    dose_range     = "fostemsavir extended-release 600 mg twice daily (193 subjects) or matching placebo (65 subjects), both added to the failing regimen during the 8-day functional-monotherapy period",
    regions        = "North America 105 (41%), South America 97 (38%), Europe 49 (19%), other 7 (3%) (Table S2)",
    notes          = paste0(
      "Same analysis set as the companion continuous-endpoint model ",
      "Parasrampuria_2025_temsavir_hivrna_day8 and as the sibling ",
      "0.5 log10 model. Fitted with NONMEM 7.2 (FOCE-I). The published ",
      "visual predictive check (Figure 3, right panel) bins the observed ",
      "data into 7 Ctau bins and overlays the median and 95% prediction ",
      "interval from 500 trial simulations. This is the endpoint carried ",
      "into the Table 5 dosing-scenario simulations, in which 500 ",
      "replicates of each phase 3 subject were simulated under eight ",
      "scenarios."
    )
  )

  ini({
    # ==================================================================
    # Parasrampuria 2025 Table 4, RIGHT pair of columns ("Proportion of
    # subjects with > 1 log10 decrease in plasma HIV-1 RNA on Day 8"),
    # Estimate (%RSE) with the 500-replicate bootstrap 95% CI alongside.
    # The equation is the Table 4 Note, verbatim:
    #
    #   Logit = E0 - (Emax * Ctau)/(EC50 + Ctau)
    #              + theta1 * (BHIVRNA - 4.65)
    #              + (theta2) * (BSL)
    #   P = exp(Logit)/(1 + exp(Logit))
    #
    # BHIVRNA is on the LOG10 scale here (4.65 = log10 of the 44,940 c/mL
    # cohort median), unlike the linear-scale (BHIVRNA/44940)^theta
    # form used by the companion continuous-endpoint model.
    #
    # SIGN. Emax is NEGATIVE (-2.67) and the equation SUBTRACTS the Emax
    # term, so the net exposure contribution is +2.67 * Ctau/(EC50+Ctau)
    # and the logit rises with exposure. Do not "correct" the sign of
    # either the estimate or the operator.
    #
    # Which column is which: Results 3.2.2 quote the > 1.0 log10 RSEs as
    # 81.4% (EC50) and 45.4% (baseline CD4+), which are the RIGHT
    # column's values; the > 0.5 log10 RSEs quoted there (33.3%, 71.2%,
    # 29.1%) are the LEFT column's. The merged Table 4 header repeats the
    # "> 0.5 log10" caption over both column pairs; the prose
    # disambiguates it. See the vignette Errata.
    #
    # UNITS. Table 4 labels the E0 and Emax rows "(log 10 c/mL)". That is
    # a carry-over from the Table 3 layout: on a logit scale these
    # coefficients are unitless. See the vignette Errata.
    # ==================================================================

    # ----- Emax linear predictor -----
    # E0 and Emax are both negative and so are NOT log-transformed;
    # EC50 is positive and is.
    e0    <- -1.99      ; label("Logit intercept E0 at zero temsavir exposure, median baseline plasma HIV-1 RNA and a baseline CD4+ count of at least 20 cells/mm3 (unitless logit)")  # Parasrampuria 2025 Table 4, > 1 log10 column, E0 = -1.99 (RSE 21.4%), bootstrap 95% CI -2.98 to -1.37
    emax  <- -2.67      ; label("Maximal exposure contribution Emax entering the logit as -(Emax * Ctau)/(EC50 + Ctau); negative, so the net effect raises the logit (unitless logit)")  # Parasrampuria 2025 Table 4, > 1 log10 column, Emax = -2.67 (RSE 20.5%), bootstrap 95% CI -4.24 to -1.98
    lec50 <- log(78.8)  ; label("Steady-state temsavir Ctau producing half of the maximal logit shift, EC50 (ng/mL)")  # Parasrampuria 2025 Table 4, > 1 log10 column, EC50 = 78.8 ng/mL (RSE 81.4%), bootstrap 95% CI 18.6-408

    # ----- Covariate effects on the logit -----
    e_hiv_vload_logit <- 0.848  ; label("Change in the logit per 1 log10 c/mL increase in baseline plasma HIV-1 RNA above 4.65 log10 c/mL (unitless logit)")  # Parasrampuria 2025 Table 4, > 1 log10 column, "Effect of Baseline Plasma HIV-1 RNA" = 0.848 (RSE 19.5%), bootstrap 95% CI 0.608-1.17
    e_cd4_abs_logit   <- -0.781 ; label("Change in the logit for a baseline CD4+ count below 20 cells/mm3 (unitless logit)")  # Parasrampuria 2025 Table 4, > 1 log10 column, "Effect of Baseline CD4 +" = -0.781 (RSE 45.4%), bootstrap 95% CI -1.44 to -0.213

    # ----- No between-subject variability, no residual error -----
    # The source is a logistic regression fitted to binary Day 8
    # responder status; a Bernoulli likelihood has no sigma, and Table 4
    # reports no random effect. The tiny fixed additive residual below
    # exists only so rxode2 has an error model to attach to the
    # typical-value probability; it is NOT a published quantity. See the
    # vignette's Assumptions and deviations.
    addSd_prob_hivrna_decr1log <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
  })

  model({
    ec50 <- exp(lec50)

    # BSL from the Table 4 Note: 1 when the baseline CD4+ count is below
    # 20 cells/mm3. Dichotomised here rather than in the dataset so the
    # threshold stays visible.
    bslCd4Low <- (CD4_ABS < 20)

    # Baseline viral load enters on the log10 scale, centred at the
    # cohort median 4.65 log10 c/mL. log10() is rxode2's base-10 log.
    hivRnaLog10Ctr <- log10(HIV_VLOAD) - 4.65

    # Table 4 Note, verbatim.
    logit_hivrna_decr1log <- e0 -
      (emax * CTROUGH) / (ec50 + CTROUGH) +
      e_hiv_vload_logit * hivRnaLog10Ctr +
      e_cd4_abs_logit * bslCd4Low

    prob_hivrna_decr1log <- expit(logit_hivrna_decr1log)

    prob_hivrna_decr1log ~ add(addSd_prob_hivrna_decr1log)
  })
}
