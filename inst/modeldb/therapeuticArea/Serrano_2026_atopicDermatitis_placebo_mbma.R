Serrano_2026_atopicDermatitis_placebo_mbma <- function() {
  description <- paste0(
    "MBMA. Longitudinal model-based meta-analysis of the PLACEBO-arm EASI-75 ",
    "responder rate in moderate-to-severe atopic dermatitis (AD), fitted to ",
    "41 placebo arms from 40 randomised Phase 2 and Phase 3 trials published ",
    "2014-2024 (4827 randomised patients). There is no drug, no dose and no ",
    "PK layer: the model describes how the placebo response itself rises and ",
    "plateaus over the double-blind period, so that a trial designer can ",
    "predict the control arm of a future AD trial. The number of EASI-75 ",
    "responders in arm i at week t is binomial(N_it, P_it) and the logit of ",
    "P_it is an exponential-approach-to-plateau in time: ",
    "logit(P) = Ebase + Em_i * (1 - exp(-k_i * t)), with Ebase = -4.97 the ",
    "day-zero logit, a typical maximum placebo effect Em = 3.21 logit units ",
    "and a typical onset rate k = 0.291 /week. Two trial-level covariates ",
    "were retained on Em by forward selection / backward elimination, both ",
    "acting on the PLATEAU rather than on the onset rate: arms that permitted ",
    "concomitant topical corticosteroids reach a 0.774 logit-unit higher ",
    "plateau (a 1.84-fold higher EASI-75 rate at Week 12), and each one-point ",
    "increase in the arm's MEAN BASELINE EASI score above the across-trial ",
    "mean of 29 lowers the plateau by 0.0486 logit units (a 0.961-fold lower ",
    "EASI-75 rate). No covariate reached significance on the onset rate, ",
    "which the source explicitly flags as unexplained between-trial ",
    "variability. Between-trial random effects are additive on Em ",
    "(variance 0.0396) and exponential on k (variance 0.138). Simulation ",
    "scope is STUDY-ARM-MEAN placebo trajectories, NOT individual patients: ",
    "the source's residual is the binomial sampling variance of an arm ",
    "proportion, P*(1-P)/N, which depends on the arm size N and is therefore ",
    "left to downstream simulation code (same convention as ",
    "modellib('Chen_2025_methotrexate_acr20_mbma') and ",
    "modellib('Boucher_2018_naproxen_mbma')). The first-order autoregressive ",
    "within-trial residual correlation (phi = 0.759) is likewise not ",
    "representable in the nlmixr2 residual model; both are reported in ini() ",
    "and in the validation vignette."
  )

  reference <- paste(
    "Serrano JC, Maringwa J, Straetemans R, Willems W, Liva SG, Verhoeven J,",
    "Ford JL, Huang K-HG, Hubbard JJ, French JL, Devineni D, Vermeulen A,",
    "Valiathan C. A Model-Based Meta-Analysis Framework Quantifying Drivers",
    "of Placebo Response in Atopic Dermatitis Trials.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70150.",
    "doi:10.1002/psp4.70150. PMC12896390.",
    "The structural model is the unnumbered display equation of Methods",
    "section 2.4 'MBMA Framework'; the final covariate model, including its",
    "fully numeric form, is the unnumbered display equation of Results",
    "section 3.2; all parameter values are Table 4.",
    sep = " "
  )

  vignette <- "Serrano_2026_atopicDermatitis_placebo_mbma"

  # This model has no drug, no dosing events and no concentration. The `units`
  # entries below follow the placeholder convention already used by the
  # arm-level responder-rate MBMA siblings (Chen_2025_methotrexate_acr20_mbma,
  # Dodds_2013_psoriasis_biologics_mbma) so that checkModelConventions() sees a
  # parseable dosing / concentration pair.
  units <- list(
    time          = "week (weeks since randomisation; the onset rate k is reported in week^-1 and the source's primary read-outs are Weeks 12 and 16)",
    dosing        = "n/a (placebo arms only; this model consumes NO rxode2 dose events and has no exposure driver)",
    concentration = "probability/arm (prob_easi75 is the STUDY-ARM probability that a placebo patient achieves a 75% reduction from baseline in the Eczema Area and Severity Index, on a 0-1 scale; it is NOT a drug concentration. The slash satisfies checkModelConventions unit parsing.)"
  )

  covariateData <- list(
    SCORE_EASI = list(
      description        = "Study-arm MEAN baseline Eczema Area and Severity Index score, as reported by the source trial.",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "TRIAL-ARM-LEVEL, not subject-level: this is the arm's reported mean baseline EASI, which is the grain at which a model-based meta-analysis of published aggregate data operates. Entered CENTERED on 29, the across-trial mean baseline EASI (Serrano 2026 Results 3.2 numeric equation, 'EASI_i - 29', and Results 3.3, 'a baseline EASI score of 29 (mean EASI value across trials)'). Note that 29 is the MEAN and is distinct from the median of 29.6 reported in Table 2; the centering constant is the mean, as Methods 2.4 states the continuous-covariate form is theta_X * (X - Xbar). Effect acts on the PLATEAU (Em), not on the onset rate. Observed across-arm range 11.1-34.5 (Table 2); the source's own tertile mid-values of 19 / 29 / 34 (Figure 2B) bound the range over which the linear-in-logit effect was exercised, and extrapolation far outside 11-35 is unsupported. Two arms of one trial (Study 203121) had mean baseline EASI < 16 and satisfied the moderate-to-severe inclusion criteria via the BSA and IGA thresholds instead; a leave-one-out sensitivity analysis (Results 3.2, Figure S8) showed excluding that trial widened the confidence interval on this effect by about 50% without shifting the point estimates more than 15%, so those low-severity arms IMPROVE rather than bias the estimate.",
      source_name        = "mean baseline EASI score / EASI_i (Serrano 2026 Table 2 'Baseline EASI (score)'; Results 3.2 final equation)"
    ),
    CONMED_STEROID_TOPICAL = list(
      description        = "1 = the trial protocol permitted concomitant TOPICAL corticosteroid (TCS) therapy in the placebo arm, 0 = no concomitant therapy was permitted.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant therapy). The typical maximum placebo effect emax therefore refers to a monotherapy-design placebo arm.",
      notes              = "TRIAL-DESIGN flag at the ARM level, not a record of any individual patient's steroid use: it encodes whether the protocol allowed background TCS, which the source treats as the design decision a trial designer controls. 14 of the 41 placebo arms permitted TCS and 27 did not (Table 2). This is the single strongest driver in the analysis: dropping it costs 28 AIC points (Table 3, p = 6.80e-8). Distinct from the register's CONMED_STEROID, which is SYSTEMIC corticosteroid administration recorded per subject; the route and the grain both differ and the two columns must not be substituted for one another. The source's Discussion argues the effect is pharmacological rather than a nuisance: permitting TCS 'transform[s] what constitutes a placebo response from natural disease fluctuation to active management effects'.",
      source_name        = "I_TCS / TCS concomitant therapy (Serrano 2026 Results 3.2 final equation; Table 4 'TCS therapy effect')"
    )
  )

  # Screened in the source's forward-selection covariate analysis but NOT
  # retained in the final model. Documented here so the covariate screen is
  # preserved without triggering a 'declared but not referenced' warning.
  # Only covariates that already have a register canonical are listed; the
  # remaining rejected covariates (affected body surface area, disease
  # duration, prior therapy, trial start year, study phase, proportion of
  # males) are described in population$notes, following the
  # Chen_2025_methotrexate_acr20_mbma precedent of not minting new canonicals
  # for covariates a paper rejected.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Study-arm mean patient age. Screened as an MBMA covariate but NOT retained in the final model.",
      units              = "year",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed in Serrano 2026 Methods 2.2 among the demographics screened ('mean: age, weight, proportion of males'). Results 3.2: 'Trial study phase, mean age, mean weight and mean disease duration showed no significant covariate effects and were excluded from the final model.' No point estimate is reported. Table 2 gives a median across arms of 36.6 years (range 27.9-44.1).",
      source_name        = "Age (years) (Serrano 2026 Table 2; Methods 2.2 covariate list)"
    ),
    WT = list(
      description        = "Study-arm mean patient body weight. Screened as an MBMA covariate but NOT retained in the final model.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Listed in Serrano 2026 Methods 2.2 and rejected in Results 3.2 alongside age, disease duration and study phase. No point estimate is reported. Table 2 gives a median across arms of 75.0 kg (range 65.6-81.1). Weight would not be expected to matter here in any case: there is no drug and therefore no exposure for weight to scale.",
      source_name        = "Weight (kg) (Serrano 2026 Table 2; Methods 2.2 covariate list)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 4827L,
    n_studies      = 40L,
    n_arms         = 41L,
    age_range      = "study-arm means 27.9-44.1 years, median across arms 36.6 (Table 2). Trials primarily enrolled adults.",
    weight_range   = "study-arm means 65.6-81.1 kg, median across arms 75.0 (Table 2)",
    sex_female_pct = 44.9,
    race_ethnicity = "Not reported at arm level; the source lists geographic region among the factors it could NOT evaluate for lack of data (Discussion limitations).",
    disease_state  = "Moderate-to-severe atopic dermatitis. Trials were required to enrol patients meeting at least two of: EASI >= 16, affected body surface area >= 10%, Investigator Global Assessment >= 3. Study-arm mean baseline EASI ranged 11.1-34.5 (median 29.6) and mean affected BSA 14.5-62.0% (median 48.6%); arms with mean baseline EASI < 16 qualified via the BSA and IGA thresholds.",
    dose_range     = "n/a (placebo arms only)",
    timepoints     = "EASI-75 responder counts at the timepoints each trial reported over its double-blind period; the source's primary read-outs are Week 12 and Week 16. Observed placebo EASI-75 rates ranged 4.6-36.8% at Week 12 and 6.1-39.4% at Week 16 (Results 3.1).",
    regions        = "International; 40 randomised trials published 2014-2024 (search window 2000-2024), 18 Phase 2 and 22 Phase 3. Table 1 lists every trial with its NCT number and data source.",
    notes          = "MBMA at the STUDY-ARM level: each modelled observation is one placebo arm's EASI-75 responder count at one timepoint, so the random effects are BETWEEN-TRIAL, not between-subject, and this model must not be used to simulate individual patients. sex_female_pct is derived as 100 - 55.1, the complement of the Table 2 median 'Percent of males (%)' of 55.1 (range 35.7-82.2), and is therefore a median across arms rather than a pooled patient proportion. One trial contributed TWO placebo arms with different dosing frequencies, which is why 40 trials give 41 arms. Covariate values missing from a publication were filled by random-forest imputation (Methods 2.2), with Table S4 confirming arm-mean trends were preserved; the supplement is NOT on disk, so the per-covariate imputation fractions are unknown (see vignette Errata). Six further screened covariates were rejected and are recorded here in prose because the register has no canonical for them and this file mints none for a rejected covariate: mean affected body surface area (strongly collinear with baseline EASI, r = 0.96, and therefore not separately identifiable), mean disease duration, prior therapy, trial start year, study phase (Phase 2 vs Phase 3), and the proportion of male patients. Of these, affected BSA DID show a significant univariate inverse correlation with the Week 12 placebo rate (r = -0.39, p < 0.05) but lost to baseline EASI in the stepwise selection. Fitted in R 4.2.1 with nlme 3.1-157 by maximum likelihood, not in NONMEM."
  )

  ini({
    # ========================================================================
    # PROVENANCE FOR THE WHOLE ini() BLOCK
    #
    # Structural model, Serrano 2026 Methods 2.4 (display equation; recovered
    # with pdftotext because the preprocessed markdown drops it):
    #
    #   Y_ij ~ binomial(N_ij, P_ij)
    #   log(P_ij / (1 - P_ij)) = Ebase + Em_i * (1 - exp(-k_i * t_ij))
    #   Em_i = Em + f_Em(X, beta) + eta_i,   eta_i ~ N(0, omega^2_eta)
    #   k_i  = k * exp(f_k(X, beta) + gamma_i), gamma_i ~ N(0, omega^2_gamma)
    #   f(X, theta) = theta_X * (X - Xbar)      [continuous covariates]
    #
    # Final covariate model, Serrano 2026 Results 3.2, printed BOTH
    # symbolically and with every number substituted:
    #
    #   log(P_ij/(1-P_ij))
    #     = Ebase + (Emax + theta_TCS*I_TCS,i + theta_EASI*(EASI_i - EASIbar))
    #       * (1 - exp(-k * t_ij))
    #     = -4.97 + (3.21 + 0.774*I_TCS,i - 0.0486*(EASI_i - 29))
    #       * (1 - exp(-0.291 * t_ij))
    #
    # NOTE that the final model carries NO covariate on k, so k_i reduces to
    # k * exp(gamma_i). The source is explicit that this is a finding, not an
    # omission (Discussion): 'This limitation is particularly evident in the
    # placebo response onset rate, where none of the evaluated covariates
    # showed significant effects despite observed inter-trial variability,
    # indicating that key determinants remain unidentified.'
    #
    # All point estimates below are Serrano 2026 Table 4. Nothing is fixed()
    # except the placeholder residual: every structural value was estimated.
    #
    # SOURCE-TRACE CONFIRMATION. The parameters below were confirmed against
    # SIX published numbers that are NOT parameter estimates, every one of
    # which they reproduce to the printed precision:
    #   Week 16, EASI 29, no TCS     -> 14.3%  (Results 3.3; model 14.30%)
    #   Week 16, EASI 29, TCS        -> 26.4%  (Results 3.3; model 26.43%)
    #   Week 12 TCS / no-TCS ratio   -> 1.84   (Results 3.3; model 1.840)
    #   Week 16 TCS / no-TCS ratio   -> 1.85   (Results 3.3; model 1.848)
    #   Per +1 EASI point rate ratio -> 0.961  (Results 3.4; model 0.9604)
    #   EASI 34 vs 24 ratio, Week 12 -> 0.665  (Results 3.6; model 0.6654)
    # The vignette shows the arithmetic for all six.
    # ========================================================================

    e0 <- -4.97
    label("Baseline placebo effect at time zero, on the EASI-75 LOGIT scale (paper: E_base; unitless log-odds)")
    # Serrano 2026 Table 4, 'Baseline response (E base)' = -4.97 (%RSE 5.0;
    # 95% CI -5.47, -4.45). expit(-4.97) = 0.687%, i.e. the model places a
    # near-zero EASI-75 rate at randomisation, as it must: EASI-75 is a 75%
    # reduction from each patient's own baseline, so nobody is a responder at
    # t = 0 by construction.

    emax <- 3.21
    label("Typical maximum placebo effect, i.e. the plateau increment on the EASI-75 LOGIT scale for a reference arm (paper: E_max; unitless log-odds)")
    # Serrano 2026 Table 4, 'Maximum placebo effect (E max)' = 3.21 (%RSE 7.6;
    # 95% CI 2.71, 3.73). The reference arm is one with NO concomitant therapy
    # and a mean baseline EASI of exactly 29, because both covariates enter
    # centered / indicator-coded to vanish there.

    lkp_easi75 <- log(0.291)
    label("Log of the typical placebo-response onset rate constant (paper: k; back-transform k = 0.291 /week)")
    # Serrano 2026 Table 4, 'Rate constant (k)' = 0.291 (%RSE 11.0; 95% CI
    # 0.227, 0.353), stated in Methods 2.4 to have 'units of week -1'. Held in
    # log space here because the source's own random effect is EXPONENTIAL on
    # k (Methods 2.4: k_i = k * exp(f_k(X, beta) + gamma_i)), so log(k) is the
    # scale on which gamma is additive and normal. Corresponding half-time to
    # the plateau is log(2)/0.291 = 2.38 weeks, consistent with the source's
    # finding that the response has essentially plateaued by Week 12
    # (1 - exp(-0.291*12) = 96.96% of Em).

    # ------------------------------------------------------------------------
    # Covariate effects on the PLATEAU (emax). Both are additive on the logit
    # scale; neither acts on the onset rate.
    # ------------------------------------------------------------------------
    e_conmed_steroid_topical_emax <- 0.774
    label("Additive shift in the maximum placebo effect for a trial arm permitting concomitant topical corticosteroids (paper: theta_TCS; logit units)")
    # Serrano 2026 Table 4, 'TCS therapy effect (theta TCS)' = 0.774 (%RSE
    # 13.7; 95% CI 0.567, 0.983). The strongest effect in the analysis:
    # Table 3 shows dropping it moves the log-likelihood from 395 to 381 and
    # AIC from -775 to -747, LRT p = 6.80e-8.

    e_score_easi_emax <- -0.0486
    label("Additive shift in the maximum placebo effect per one-point increase in the arm's mean baseline EASI above 29 (paper: theta_EASI; logit units per EASI point)")
    # Serrano 2026 Table 4, 'EASI baseline effect (theta EASI)' = -0.0486
    # (%RSE 25.3; 95% CI -0.0721, -0.0260). NEGATIVE: sicker arms respond LESS
    # to placebo. Table 3 backward elimination gives LRT p = 0.000509. Its
    # %RSE of 25.3 is the largest in the table but still inside the source's
    # own stated 30% acceptance threshold (Methods 2.5).

    # ========================================================================
    # BETWEEN-TRIAL random effects (Serrano 2026 Methods 2.4). These are
    # BETWEEN-STUDY variances of an aggregate-data meta-analysis, NOT popPK
    # between-subject variability: one draw describes one whole trial arm.
    # Table 4 labels both rows 'Variance of random effect on ...', so the
    # tabulated numbers are VARIANCES and go into ini() unchanged. Neither
    # carries an RSE in Table 4.
    # ========================================================================
    eta_study_emax ~ 0.0396
    # Serrano 2026 Table 4, 'Variance of random effect on E max (omega^2 eta)'
    # = 0.0396; ADDITIVE on emax per Methods 2.4 (Em_i = Em + f_Em + eta_i).
    # SD = 0.199 logit units.

    eta_study_lkp_easi75 ~ 0.138
    # Serrano 2026 Table 4, 'Variance of random effect on k (omega^2 gamma)'
    # = 0.138; the source places gamma INSIDE an exponential
    # (k_i = k * exp(... + gamma_i)), which is exactly an additive eta on
    # log(k), hence the name. SD = 0.371 on the log scale, i.e. a 1-SD trial
    # has an onset rate 1.45-fold faster or 0.69-fold slower than typical.
    # This is the larger of the two random effects and the source flags it as
    # unexplained: no screened covariate reached significance on k.

    # ========================================================================
    # RESIDUAL ERROR. The source has NO estimated residual variance: the
    # observation is a responder COUNT and the residual is the binomial
    # sampling variance of the arm proportion, stated in Methods 2.4 as
    #
    #   var(Y_ij / N_ij | eta_i, gamma_i) = P_ij * (1 - P_ij) / N_ij
    #
    # which 'reflects the binomial sampling variability and ensures that
    # observations are weighted by their sample size'. Because that variance
    # depends on the ARM SIZE N_ij, which is data rather than a parameter, it
    # cannot live in ini(); reproducing it is left to downstream simulation
    # code, the same convention Chen_2025_methotrexate_acr20_mbma and
    # Boucher_2018_naproxen_mbma use for their 1/sqrt(N) arm weights. The
    # validation vignette shows the two-line recipe.
    #
    # Serrano 2026 additionally fits a FIRST-ORDER AUTOREGRESSIVE within-trial
    # residual correlation, Table 4 'Autoregressive coefficient (phi)' =
    # 0.759, to account for the fact that the same arm is observed repeatedly
    # over its double-blind period. nlmixr2's residual-error models are
    # independent across records and cannot express an AR(1) structure, so
    # phi is recorded here and in the vignette Errata but is NOT encoded.
    # Neither omission changes any TYPICAL-VALUE prediction: both affect only
    # the spread and the correlation of simulated arm proportions about the
    # trajectory, not the trajectory itself, which is why every published
    # number in the source-trace confirmation above is reproduced exactly.
    #
    # The placeholder below exists only so the nlmixr2 likelihood machinery
    # accepts the model for forward simulation. It is NOT from the source --
    # same device as Bhatnagar_2024_upadacitinib_asas20_as.R and
    # Chen_2021_lorlatinib_icorr.R.
    # ========================================================================
    addSd_prob_easi75 <- fixed(0.001)
    label("Placeholder additive residual SD on the arm-level probability output prob_easi75 (unitless); NOT from the source, whose residual is the binomial sampling variance P*(1-P)/N")
  })

  model({
    # ---- Trial-specific maximum placebo effect (the PLATEAU) ---------------
    # Serrano 2026 Results 3.2. Both covariates are additive on the logit
    # scale and both vanish at the reference arm (no concomitant therapy,
    # mean baseline EASI of 29), where emax_i reduces to emax + eta.
    #
    # The centering constant 29 is the across-trial MEAN baseline EASI, hard
    # coded because the source prints it inside its own final numeric equation
    # ('- 0.0486 * (EASI_i - 29)') and restates it in Results 3.3 as 'a
    # baseline EASI score of 29 (mean EASI value across trials)'. It is NOT
    # the Table 2 median of 29.6; Methods 2.4 defines the continuous-covariate
    # form as theta_X * (X - Xbar) with Xbar the mean.
    emax_i <-
      emax +
      e_conmed_steroid_topical_emax * CONMED_STEROID_TOPICAL +
      e_score_easi_emax * (SCORE_EASI - 29) +
      eta_study_emax

    # ---- Trial-specific onset rate ----------------------------------------
    # k_i = k * exp(gamma_i). No covariate was retained on the onset rate, so
    # the only trial-to-trial movement here is the random effect.
    kp_easi75 <- exp(lkp_easi75 + eta_study_lkp_easi75)

    # ---- Exponential approach to the plateau, on the logit scale ----------
    # At time 0 the bracket is 0 and the logit is e0 alone, i.e. essentially
    # nobody has yet achieved a 75% EASI reduction. As time grows the bracket
    # approaches 1 and the logit approaches e0 + emax_i.
    lp_easi75 <- e0 + emax_i * (1 - exp(-kp_easi75 * time))

    # ---- Study-arm EASI-75 placebo responder probability ------------------
    prob_easi75 <- expit(lp_easi75)

    prob_easi75 ~ add(addSd_prob_easi75)
  })
}
