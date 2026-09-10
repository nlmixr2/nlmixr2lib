Parasrampuria_2025_temsavir_hivrna_day8 <- function() {
  description <- paste0(
    "Inhibitory Emax exposure-response model relating the steady-state ",
    "trough plasma concentration of temsavir (Ctau) to the change in ",
    "plasma HIV-1 RNA from Day 1 to Day 8 during fostemsavir functional ",
    "monotherapy, in heavily treatment-experienced adults with ",
    "multidrug-resistant HIV-1 (Parasrampuria 2025, n = 258 from the ",
    "randomised cohort of the phase 3 BRIGHTE study AI438047 / ",
    "NCT02362503: 193 on fostemsavir 600 mg BID and 65 on placebo). This ",
    "is the phase 3 primary endpoint. The change in plasma HIV-1 RNA ",
    "(log10 c/mL, negative = a decline) is ",
    "E0 - Emax * Ctau/(EC50 + Ctau) * (HIV_VLOAD/44940)^0.150 * ",
    "0.481^(CD4_ABS < 20), with E0 = -0.129, Emax = 1.00 and ",
    "EC50 = 64.3 ng/mL (Parasrampuria 2025 Table 3 and its Note). Fitted ",
    "naive-pooled with additive residual error; there is no ",
    "between-subject random effect. E0 and EC50 are IMPRECISE ",
    "(RSE 66.1% and 98.0%) because almost every phase 3 patient sat at or ",
    "near the plateau of the exposure-response curve, so the model has ",
    "little predictive power below the observed exposure range -- Emax, ",
    "by contrast, is well estimated (RSE 17.1%). Higher baseline plasma ",
    "HIV-1 RNA and a baseline CD4+ count of 20 cells/mm3 or more each ",
    "give a larger decline; every other virologic, immunologic and ",
    "demographic covariate tested was rejected. There is no PK layer and ",
    "no ODE: the exposure metric is supplied as a data column, derived ",
    "post hoc from the companion population PK model ",
    "modellib('Parasrampuria_2025_temsavir'). Two companion logistic ",
    "exposure-response models cover the corresponding responder-rate ",
    "endpoints."
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
    concentration = "dviralLoad (change in plasma HIV-1 RNA from Day 1 to Day 8, log10 c/mL; negative = a decline)"
  )

  covariateData <- list(
    CTROUGH = list(
      description        = "Individual steady-state plasma temsavir concentration at the end of the dosing interval (Ctau), per subject. Supplied as data: this model has no PK layer.",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "TOTAL (not unbound) plasma temsavir at STEADY STATE, on the",
        "approved fostemsavir 600 mg BID regimen, so the dosing interval",
        "is 12 h. Parasrampuria 2025 Methods 2.3: post-hoc PK parameter",
        "estimates from the final population PK model were used to derive",
        "individual concentration-time profiles, and the steady-state",
        "exposure metrics Cavg, Cmax and Ctau were then computed from",
        "those simulated profiles with the PKNCA R package -- i.e. these",
        "are empirical-Bayes model predictions, not observed troughs.",
        "Reproduce them with modellib('Parasrampuria_2025_temsavir')",
        "solved at 600 mg BID with rate = -2 on the dose records.",
        "The paper's Ctau and Cavg gave very similar fits ('graphical",
        "analysis and model fits were very similar between Ctau and Cavg,",
        "data not shown'); Ctau was chosen so the results could be",
        "compared with the earlier phase 2 analysis. Normalising Ctau by",
        "the protein-binding-adjusted baseline IC50 or IC90 (Ctau/PBIC50,",
        "Ctau/PBIC90) did NOT improve the fit, so the bare concentration",
        "is what the final model uses. Simulated median for the",
        "600 mg BID phase 3 cohort: 433 ng/mL (95% CI 33.6-2400)",
        "(Table 5). Note the observed Ctau spans values at which the",
        "model does not describe the data well: some subjects with",
        "Ctau < 10 ng/mL achieved a > 0.5 log10 decline while some with",
        "Ctau > 4000 ng/mL did not."
      ),
      source_name        = "Ctau"
    ),
    HIV_VLOAD = list(
      description        = "Baseline (Day 1) plasma HIV-1 RNA concentration.",
      units              = "copies/mL (c/mL)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Carried LINEAR in c/mL, as the <VIRUS>_VLOAD family requires, and",
        "normalised inside model() to the cohort median. The Table 3 Note",
        "writes the term as (BHIVRNA/44940)^theta1, so the centring",
        "constant is 44,940 c/mL -- the median quoted in Results 3.2.1.",
        "Table S2 prints the same median as 44,943 c/mL (range",
        "39-8,207,000, N = 258); the 3 c/mL discrepancy is a rounding",
        "difference and the equation's 44,940 is what is encoded, because",
        "the fitted intercept absorbs the centring constant. Note that the",
        "SAME column is consumed on the log10 scale by the two companion",
        "logistic models, where the median is written as 4.65 log10 c/mL",
        "(log10(44940) = 4.65263); this is exactly why the register",
        "prescribes carrying the linear concentration and transforming",
        "inside model(). Higher baseline viral load gives a larger Day 8",
        "decline: Results 3.2.1 report that for subjects with baseline",
        "CD4+ below 20 cells/mm3, virologic response was 33.4% lower at a",
        "baseline of 1000 c/mL than at the 44,940 c/mL median."
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
        "at 20 cells/mm3, following the pattern the register prescribes",
        "for HCV_VLOAD and ALP: a pre-binarised dataset column would hide",
        "the threshold. The Table 3 Note's BSL is the resulting flag and",
        "the Table 3 abbreviations define it as 'baseline CD4+ flag",
        "(>= 20 cells/mm3 or < 20 cells/mm3)'.",
        "ORIENTATION. BSL = 1 for CD4+ BELOW 20 cells/mm3 and the 0.481",
        "multiplier therefore SHRINKS the drug effect in the",
        "severely-immunosuppressed group. The paper does not state the",
        "coding direction, but its own two reported percentages fix it and",
        "cross-check each other: at Ctau = 510 ng/mL this encoding",
        "reproduces BOTH the 33.4% baseline-viral-load contrast AND the",
        "45.3% CD4 contrast to three significant figures (see the",
        "vignette). The reversed coding gives the wrong sign for the CD4",
        "contrast, contradicting Results 3.2.1 ('baseline CD4+ counts",
        ">= 20 cells/mm3 resulted in greater reductions in plasma HIV-1",
        "RNA'). The two companion logistic models carry the same flag with",
        "a negative logit coefficient, which is the same orientation.",
        "Day 8 cohort (N = 258): 69 subjects (27%) below 20 cells/mm3;",
        "median 98.5 cells/mm3, range 0-1160 (Table S2)."
      ),
      source_name        = "BSL (derived from baseline CD4+ count)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened for the Day 8 efficacy exposure-response model (Table 1) and rejected: Results 3.2.1 state that 'other tested covariates (virologic, immunologic, and demographic factors) did not affect exposure-virologic response relationships.' Day 8 cohort median 48 years, range 18-73 (Table S2)."
    ),
    SEXF = list(
      description = "Female sex indicator; 1 = female, 0 = male.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Table 1, as 'gender') and rejected. Day 8 cohort 67 female (26%), 191 male (74%) (Table S2)."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (Table 1, as 'race') and rejected. Day 8 cohort 60 (23%) (Table S2)."
    ),
    WT = list(
      description = "Baseline body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (Table 1) and rejected. Day 8 cohort median 70 kg, range 38-146 (Table S2). Body weight does act on exposure through the companion population PK model; it simply adds nothing once Ctau is in the model."
    ),
    REGION = list(
      description = "Geographic region of the enrolling site.",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened (Table 1, as 'geographic region') and rejected. Day 8 cohort: North America 105 (41%), South America 97 (38%), Europe 49 (19%), other 7 (3%) (Table S2). Not a register canonical -- listed here for provenance only, since covariatesDataExcluded is documentation."
    ),
    IC50_TMR = list(
      description = "Baseline phenotypic 50% inhibitory concentration of temsavir against the subject's own virus.",
      units       = "nM",
      type        = "continuous",
      notes       = "Screened as a baseline virologic-sensitivity covariate and, in the exposure-metric screen, as the denominators PBIC50 / PBIC90 / IC50-fold-change; none improved the fit (Results 3.2.1: 'Evaluating Ctau by IC50 and IC90 adjusted for protein binding (Ctau/PBIC50, Ctau/PBIC90) did not improve model fits'). Day 8 cohort median 0.88 nM, range 0.04-6000 (0.417 ng/mL, range 0.019-2841). Table 1 gives the protein-binding adjustment: PBIC50 = 473.48 * (IC50/fu), where 473.48 g/mol is the molecular weight of free-base temsavir and fu = 0.12 is the mean estimated unbound fraction in vivo. Not a register canonical -- documentation only."
    ),
    GP160_SUBS = list(
      description = "Number of pre-defined genotypic substitutions of interest within the HIV-1 gp160 domain at baseline.",
      units       = "(count)",
      type        = "count",
      notes       = "Screened (Table 1) and rejected. Day 8 cohort: 0 substitutions 141 (55%), 1 substitution 102 (40%), 2 substitutions 15 (6%) (Table S2). Not a register canonical -- documentation only."
    ),
    CD8_ABS = list(
      description = "Baseline absolute CD8+ T-lymphocyte count.",
      units       = "cells/mm^3",
      type        = "continuous",
      notes       = "Screened (Table 1) and rejected. Day 8 cohort median 653 cells/mm3, range 61-2700 (Table S2). Not a register canonical -- documentation only."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 258L,
    n_studies      = 1L,
    n_observations = "258 Day 8 change-from-Day-1 plasma HIV-1 RNA records (one per subject; landmark analysis, no repeated measures)",
    age_range      = "18-73 years (median 48) (Table S2)",
    weight_range   = "38-146 kg (median 70) (Table S2)",
    sex_female_pct = 26.0,
    race_ethnicity = c(White = 66.0, `Black or African American` = 23.0, Asian = 1.0, Other = 10.0),
    disease_state  = "heavily treatment-experienced (HTE) adults with multidrug-resistant HIV-1 infection failing their current antiretroviral regimen, from the randomised cohort of the phase 3 BRIGHTE study; median baseline plasma HIV-1 RNA 44,943 c/mL (4.65 log10 c/mL, range 1.59-6.91), median baseline CD4+ count 98.5 cells/mm3 (range 0-1160) with 69 subjects (27%) below 20 cells/mm3",
    dose_range     = "fostemsavir extended-release 600 mg twice daily (193 subjects) or matching placebo (65 subjects), both added to the failing regimen during the 8-day functional-monotherapy period",
    regions        = "North America 105 (41%), South America 97 (38%), Europe 49 (19%), other 7 (3%) (Table S2)",
    notes          = paste0(
      "The exposure-response analysis was deliberately restricted to the ",
      "phase 3 study because that is the target population. Subjects ",
      "needed both a PK sample and Day 8 plasma HIV-1 RNA data. Fitted ",
      "with NONMEM 7.2 (FOCE-I) by a naive-pooled approach with additive ",
      "residual error; residual variability is high, about 65% of Emax. ",
      "Simulated virologic responses were consistent with the phase 3 ",
      "primary endpoint, a 0.791 log10 c/mL adjusted mean decline at ",
      "Day 8, and the estimated E0 of -0.129 log10 c/mL is consistent ",
      "with the observed placebo-plus-failing-regimen adjusted mean ",
      "decline of -0.1666 log10 c/mL. The corresponding Week 24 ",
      "endpoints (change in plasma HIV-1 RNA, change in CD4+ count, and ",
      "the proportions below 40 / 200 / 400 c/mL) showed NO relationship ",
      "with temsavir Ctau on graphical analysis and no model was fitted, ",
      "so no Week 24 model is packaged."
    )
  )

  ini({
    # ==================================================================
    # Parasrampuria 2025 Table 3 ("Parameter estimates for final
    # exposure-response model of plasma temsavir Ctau and change in
    # plasma HIV-1 RNA from Day 1 to Day 8"), Estimate (%RSE) column,
    # with the 500-replicate bootstrap 95% CI alongside. The equation is
    # the Table 3 Note, verbatim:
    #
    #   Change in plasma HIV-1 RNA from Day 1 to Day 8 =
    #     E0 - Emax * (Ctau)/(EC50 + Ctau)
    #        * (BHIVRNA/44940)^theta1
    #        * theta2^BSL
    #   Residual error: Y = IPRED + (eps1)
    #
    # Sign convention: the minus sign in front of Emax is explicit in the
    # printed equation and Emax is POSITIVE, so a larger Emax term means
    # a LARGER DECLINE (a more negative change). E0 is the placebo
    # intercept and is itself negative.
    #
    # WHAT SETTLES THE BSL ORIENTATION. Table 3 defines BSL only as the
    # "baseline CD4+ flag (>= 20 cells/mm3 or < 20 cells/mm3)" and never
    # says which level scores 1. The paper's own two reported contrasts
    # settle it jointly -- two equations, one unknown Ctau -- and they
    # are mutually consistent only under BSL = 1 for CD4+ < 20:
    #
    #   Results 3.2.1, contrast A: "At a baseline CD4+ count
    #     < 20 cells/mm3, the Day 8 virologic response was 33.4% lower
    #     for a baseline plasma HIV-1 RNA level of 1000 c/mL compared
    #     with the median baseline plasma HIV-1 RNA of 44 940 c/mL."
    #   Results 3.2.1, contrast B: "At a baseline plasma HIV-1 RNA value
    #     of 44 940 c/mL, Day 8 virologic response was 45.3% higher for
    #     subjects with baseline CD4+ count >= 20 cells/mm3 compared with
    #     subjects with baseline CD4+ count < 20 cells/mm3."
    #
    # Solving contrast A for the Ctau at which it holds exactly gives
    # 510 ng/mL, and at that same Ctau contrast B reproduces as 45.3%.
    # Both percentages are taken relative to the LARGER response. Under
    # the reversed coding contrast B comes out with the wrong sign,
    # contradicting the neighbouring sentence that CD4+ counts of
    # 20 cells/mm3 or more "resulted in greater reductions". The
    # vignette runs this as a gate.
    # ==================================================================

    # ----- Emax structure -----
    # E0 is the placebo effect and is negative, so it is NOT
    # log-transformed. Emax and EC50 are positive and are.
    e0    <- -0.129     ; label("Placebo intercept E0: change in plasma HIV-1 RNA from Day 1 to Day 8 at zero temsavir exposure (log10 c/mL)")  # Parasrampuria 2025 Table 3, E0 = -0.129 (RSE 66.1%), bootstrap 95% CI -0.268 to -0.00665. Imprecise by the authors' own account; the Discussion notes it agrees with the observed placebo adjusted mean decline of -0.1666 log10 c/mL
    lemax <- log(1.00)  ; label("Maximal drug-attributable decline Emax in plasma HIV-1 RNA at Day 8 (log10 c/mL)")                             # Parasrampuria 2025 Table 3, Emax = 1.00 (RSE 17.1%), bootstrap 95% CI 0.825-1.28. The only well-estimated parameter of the three
    lec50 <- log(64.3)  ; label("Steady-state temsavir Ctau producing half of Emax, EC50 (ng/mL)")                                              # Parasrampuria 2025 Table 3, EC50 = 64.3 ng/mL (RSE 98.0%), bootstrap 95% CI 16.6-250. Very imprecise: the Discussion attributes this to most phase 3 patients sitting at or near the plateau of the exposure-efficacy curve

    # ----- Covariate effects on the drug-effect term -----
    # Results 3.2.1 identify these as covariates "of Emax"; in the
    # printed equation both multiply the whole Emax term, which is the
    # same thing.
    e_hiv_vload_emax <- 0.150 ; label("Power exponent on (baseline plasma HIV-1 RNA / 44940 c/mL) scaling the Emax term (unitless)")  # Parasrampuria 2025 Table 3, "Effect of Baseline Plasma HIV-1 RNA" = 0.150 (RSE 29.6%), bootstrap 95% CI 0.113-0.197
    e_cd4_abs_emax   <- 0.481 ; label("Multiplicative effect on the Emax term for a baseline CD4+ count below 20 cells/mm3 (unitless)")  # Parasrampuria 2025 Table 3, "Effect of Baseline CD4 +" = 0.481 (RSE 19.9%), bootstrap 95% CI 0.293-0.673; applied as theta2^BSL

    # ----- Residual error -----
    # Additive on the response scale (log10 c/mL), per the Table 3 Note
    # "Residual error: Y = IPRED + (eps1)". Naive-pooled fit: there is
    # no between-subject random effect to estimate.
    addSd_dviralLoad <- 0.653 ; label("Additive residual error SD on the Day 8 change in plasma HIV-1 RNA (log10 c/mL)")  # Parasrampuria 2025 Table 3, epsilon = 0.653 (RSE 4.10%), bootstrap 95% CI 0.594-0.699. Results 3.2.1 describe this as "high residual variability (approximately 65% of Emax)"
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # BSL from the Table 3 Note: 1 when the baseline CD4+ count is below
    # 20 cells/mm3. Dichotomised here rather than in the dataset so the
    # threshold stays visible.
    bslCd4Low <- (CD4_ABS < 20)

    # Baseline-viral-load scaling, normalised to the cohort median in
    # the printed equation.
    hivRnaEff <- (HIV_VLOAD / 44940)^e_hiv_vload_emax

    # Table 3 Note, verbatim. Negative values are declines.
    dviralLoad <- e0 -
      emax * CTROUGH / (ec50 + CTROUGH) * hivRnaEff * e_cd4_abs_emax^bslCd4Low

    dviralLoad ~ add(addSd_dviralLoad)
  })
}
