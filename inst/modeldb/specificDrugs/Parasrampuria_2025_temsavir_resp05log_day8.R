Parasrampuria_2025_temsavir_resp05log_day8 <- function() {
  description <- paste0(
    "Logistic-regression exposure-response model with an Emax linear ",
    "predictor, relating the steady-state trough plasma concentration of ",
    "temsavir (Ctau) to the probability of achieving a GREATER THAN ",
    "0.5 log10 c/mL decrease in plasma HIV-1 RNA from Day 1 to Day 8 ",
    "during fostemsavir functional monotherapy, in heavily ",
    "treatment-experienced adults with multidrug-resistant HIV-1 ",
    "(Parasrampuria 2025, n = 258 from the randomised cohort of the ",
    "phase 3 BRIGHTE study AI438047 / NCT02362503: 193 on fostemsavir ",
    "600 mg BID and 65 on placebo). The probability is ",
    "expit(E0 - Emax * Ctau/(EC50 + Ctau) + 0.836 * ",
    "(log10(HIV_VLOAD) - 4.65) - 1.27 * (CD4_ABS < 20)) with ",
    "E0 = -1.07, Emax = -2.92 and EC50 = 88 ng/mL (Parasrampuria 2025 ",
    "Table 4 and its Note). Emax is NEGATIVE and is SUBTRACTED, so the ",
    "exposure term raises the logit: higher trough gives a higher ",
    "response probability. E0, EC50 and the CD4+ effect are imprecise ",
    "(RSE 33.3%, 71.2% and 29.1%), and the paper flags subjects who ",
    "contradict the relationship at both extremes -- some with observed ",
    "Ctau below 10 ng/mL achieved the decline while some above ",
    "4000 ng/mL did not. There is no PK layer, no ODE and no ",
    "between-subject random effect (Bernoulli likelihood): the exposure ",
    "metric is supplied as a data column, derived post hoc from the ",
    "companion population PK model ",
    "modellib('Parasrampuria_2025_temsavir'). The sibling model ",
    "Parasrampuria_2025_temsavir_resp1log_day8 gives the parallel fit ",
    "for the stricter 1 log10 threshold."
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
    concentration = "prob_hivrna_decr05log (probability of a > 0.5 log10 c/mL decrease in plasma HIV-1 RNA at Day 8, 0-1; also logit_hivrna_decr05log)"
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
        "rate = -2 on the dose records. Simulated median for the",
        "600 mg BID phase 3 cohort 433 ng/mL (95% CI 33.6-2400)",
        "(Table 5). Results 3.2.2 note that the relationship is",
        "contradicted at both tails -- 'some subjects with low observed",
        "Ctau (< 10 ng/mL) displayed a > 0.5 log10 reduction while some",
        "subjects with high observed Ctau (> 4000 ng/mL) did not'."
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
        "scale as (BHIVRNA/44940)^theta1 -- the two scales in one paper",
        "are exactly why the register prescribes carrying the linear",
        "concentration and transforming inside model(). The positive",
        "coefficient means a higher baseline viral load raises the",
        "probability of clearing the 0.5 log10 bar, consistent with the",
        "continuous model, where higher baseline load gives a larger",
        "decline."
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
        "probability in the severely-immunosuppressed group. That",
        "orientation is the same one the companion continuous-endpoint",
        "model requires -- where it is pinned quantitatively by two",
        "reported percentage contrasts -- and it matches Results 3.2.1",
        "('baseline CD4+ counts >= 20 cells/mm3 resulted in greater",
        "reductions in plasma HIV-1 RNA'). Day 8 cohort (N = 258): 69",
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
      "1 log10 model. Fitted with NONMEM 7.2 (FOCE-I). The published ",
      "visual predictive check (Figure 3, left panel) bins the observed ",
      "data into 7 Ctau bins and overlays the median and 95% prediction ",
      "interval from 500 trial simulations. The candidate covariate set ",
      "for this endpoint was narrower than for the continuous endpoint: ",
      "Table 1 lists only baseline plasma HIV-1 RNA and baseline CD4+, ",
      "and both were retained."
    )
  )

  ini({
    # ==================================================================
    # Parasrampuria 2025 Table 4, LEFT pair of columns ("Proportion of
    # subjects with > 0.5 log10 decrease in plasma HIV-1 RNA on Day 8"),
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
    # form used by the companion continuous-endpoint model. Both read the
    # same HIV_VLOAD column; the transform differs.
    #
    # SIGN. Emax is NEGATIVE (-2.92) and the equation SUBTRACTS the Emax
    # term, so the net exposure contribution is +2.92 * Ctau/(EC50+Ctau)
    # and the logit rises with exposure. Do not "correct" the sign of
    # either the estimate or the operator; they are consistent as
    # printed, and flipping one alone inverts the exposure-response.
    #
    # Which column is which: Results 3.2.2 quote the > 0.5 log10 RSEs as
    # 33.3% (E0), 71.2% (EC50) and 29.1% (baseline CD4+), which are the
    # LEFT column's values, and the > 1.0 log10 RSEs as 81.4% (EC50) and
    # 45.4% (CD4+), which are the RIGHT column's. The merged Table 4
    # header repeats the "> 0.5 log10" caption over both column pairs;
    # the prose disambiguates it. See the vignette Errata.
    #
    # UNITS. Table 4 labels the E0 and Emax rows "(log 10 c/mL)". That is
    # a carry-over from the Table 3 layout: on a logit scale these
    # coefficients are unitless. See the vignette Errata.
    # ==================================================================

    # ----- Emax linear predictor -----
    # E0 and Emax are both negative and so are NOT log-transformed;
    # EC50 is positive and is.
    e0    <- -1.07    ; label("Logit intercept E0 at zero temsavir exposure, median baseline plasma HIV-1 RNA and a baseline CD4+ count of at least 20 cells/mm3 (unitless logit)")  # Parasrampuria 2025 Table 4, > 0.5 log10 column, E0 = -1.07 (RSE 33.3%), bootstrap 95% CI -1.78 to -0.491
    emax  <- -2.92    ; label("Maximal exposure contribution Emax entering the logit as -(Emax * Ctau)/(EC50 + Ctau); negative, so the net effect raises the logit (unitless logit)")  # Parasrampuria 2025 Table 4, > 0.5 log10 column, Emax = -2.92 (RSE 18.0%), bootstrap 95% CI -4.11 to -2.22
    lec50 <- log(88)  ; label("Steady-state temsavir Ctau producing half of the maximal logit shift, EC50 (ng/mL)")  # Parasrampuria 2025 Table 4, > 0.5 log10 column, EC50 = 88 ng/mL (RSE 71.2%), bootstrap 95% CI 17.2-289

    # ----- Covariate effects on the logit -----
    e_hiv_vload_logit <- 0.836 ; label("Change in the logit per 1 log10 c/mL increase in baseline plasma HIV-1 RNA above 4.65 log10 c/mL (unitless logit)")  # Parasrampuria 2025 Table 4, > 0.5 log10 column, "Effect of Baseline Plasma HIV-1 RNA" = 0.836 (RSE 21.4%), bootstrap 95% CI 0.576-1.19
    e_cd4_abs_logit   <- -1.27 ; label("Change in the logit for a baseline CD4+ count below 20 cells/mm3 (unitless logit)")  # Parasrampuria 2025 Table 4, > 0.5 log10 column, "Effect of Baseline CD4 +" = -1.27 (RSE 29.1%), bootstrap 95% CI -1.96 to -0.672

    # ----- No between-subject variability, no residual error -----
    # The source is a logistic regression fitted to binary Day 8
    # responder status; a Bernoulli likelihood has no sigma, and Table 4
    # reports no random effect. The tiny fixed additive residual below
    # exists only so rxode2 has an error model to attach to the
    # typical-value probability; it is NOT a published quantity. See the
    # vignette's Assumptions and deviations.
    addSd_prob_hivrna_decr05log <- fixed(0.001) ; label("Placeholder additive residual SD on the typical-value response probability; the source likelihood is Bernoulli (no source residual)")  # not from source; see vignette Assumptions and deviations
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
    logit_hivrna_decr05log <- e0 -
      (emax * CTROUGH) / (ec50 + CTROUGH) +
      e_hiv_vload_logit * hivRnaLog10Ctr +
      e_cd4_abs_logit * bslCd4Low

    prob_hivrna_decr05log <- expit(logit_hivrna_decr05log)

    prob_hivrna_decr05log ~ add(addSd_prob_hivrna_decr05log)
  })
}
