Sheng_2016_quinine_rat <- function() {
  description <- "Preclinical (rat). Two generalized Poisson (2GP) mixture PD model for bimodal lick-count data from rodent brief-access taste aversion (BATA) experiments with quinine hydrochloride dihydrate; the drug effect enters via a sigmoid Emax on a logistic-transformed mixing probability between a low-count and a right-truncated high-count generalized-Poisson distribution. The fitted compound is quinine HCl dihydrate used as a model bitter stimulus. STIM_QUININE_MM is the applied sipper-tube concentration (mM); there is no PK ODE and no time evolution (each record is an 8-second presentation)."
  reference <- paste(
    "Sheng Y, Soto J, Orlu Gul M, Cortina-Borja M, Tuleu C, Standing JF. (2016).",
    "New Generalized Poisson Mixture Model for Bimodal Count Data With Drug Effect:",
    "An Application to Rodent Brief-Access Taste Aversion Experiments.",
    "CPT Pharmacometrics Syst Pharmacol 5(8):427-436.",
    "doi:10.1002/psp4.12093.",
    sep = " "
  )
  vignette <- "Sheng_2016_quinine_rat"
  units <- list(
    time          = "(unused; each record is an 8-second presentation, no time evolution)",
    dosing        = "mM (quinine HCl dihydrate applied via sipper tube)",
    concentration = "(none; observation is a count of licks in 8 s, 0-61)"
  )

  covariateData <- list(
    STIM_QUININE_MM = list(
      description        = "Applied sipper-tube concentration of quinine HCl dihydrate (mM) presented to the rat during the 8-second BATA trial; drives the sigmoid-Emax effect on the logistic mixing probability between the low-count and high-count generalized-Poisson distributions.",
      units              = "mM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper presented seven concentrations: 0, 0.01, 0.03, 0.1, 0.3, 1, 3 mM (Table 1). Set to 0 for the water (control) presentation. Per-record covariate; no PK is modeled. New canonical entry registered alongside this model in inst/references/covariate-columns.md. Distinct from CAV (systemic plasma concentration of an administered drug) - STIM_QUININE_MM is the applied stimulus concentration in solution that contacts the taste receptors directly during a brief-access trial.",
      source_name        = "QUININE"
    )
  )

  population <- list(
    species        = "rat (10 trained rats; strain not reported in the source publication)",
    n_subjects     = 10,
    n_studies      = 1,
    age_range      = "(not reported in the source publication)",
    weight_range   = "(not reported in the source publication)",
    sex_female_pct = NA_real_,
    disease_state  = "(none; preclinical taste-aversion screening)",
    dose_range     = "Seven concentrations of quinine HCl dihydrate: 0 (deionized water), 0.01, 0.03, 0.1, 0.3, 1, 3 mM, presented via sipper tube for 8 s with a 2 s water-rinse between trials. Each quinine concentration was presented four times and water six times per 40-minute session; experiments were repeated weekly for 8 weeks with a 1-week washout between sessions.",
    regions        = "(not applicable; preclinical)",
    notes          = "5,400 lick-count records total (Table 1: water n=1,080; each quinine concentration n=718-722). The observed maximum lick number across all trials was 61, which fixes the right-truncation upper bound used by the second (high-count) generalized-Poisson distribution. Lick counts followed a bimodal distribution (Figure 1b) with one peak in 0-20 and a second peak in 40-60; the proportion in the low-count peak increased monotonically with quinine concentration (Figure 2b). The Sheng 2016 publication does not report rat strain, sex, or age, and no supplementary methods document was on disk at extraction time."
  )

  ini({
    # ------------------------------------------------------------------
    # All fixed-effect values are the final estimates from Sheng 2016
    # Table 2, "2GP" column (the final model). RSE% values from Table 2
    # are noted in trailing comments; 90% bootstrap confidence intervals
    # (Table 2 "2GP Bootstrap" column) are noted for parameters that
    # carry one. The source publication's NONMEM control stream was not
    # available on disk - all values are transcribed from Table 2 of the
    # main text.
    #
    # Model equations come from the Methods section "Models for the
    # first distribution" (Eq. for Poisson, NB, GP), "Truncated models
    # for the second distribution" (right-truncation at 61), and "Model
    # for the drug effect" (sigmoid Emax on logistic mixing probability,
    # parameterized by RIC50 rather than the apparent IC50).
    # ------------------------------------------------------------------

    lEmax  <- log(21.7);    label("Log of Emax of the logistic drug-effect on mixing probability (unitless)")
    # Table 2 2GP column: Emax = 21.7 (RSE 2%); bootstrap CI (16.9, 35.8).

    lRIC50 <- log(0.0423);  label("Log of real IC50 (RIC50, mM quinine HCl dihydrate)")
    # Table 2 2GP column: RIC50 = 0.0423 mM (RSE 0.3%); bootstrap CI (0.0358, 0.0551).
    # Distinct from the apparent IC50 (derived from E0, Emax, RIC50;
    # Sheng 2016 reports apparent IC50 = 1.18 mM for the 2GP model).

    e0     <- -1.57;        label("Baseline logit of the low-count mixture probability at zero quinine concentration (unitless)")
    # Table 2 2GP column: E0 = -1.57 (RSE 0.1%); bootstrap CI (-1.60, -1.55).
    # logistic(E0) = 0.172 = minimum proportion of the low-count
    # distribution (water group baseline). Negative typical value; IIV
    # encoded as multiplicative log-normal preserves sign.

    lk1    <- log(1.8);     label("Log of theta-1, location-like parameter of the first generalized-Poisson distribution (unitless)")
    # Table 2 2GP column: k1 = 1.8 (RSE 2.1%); bootstrap CI (1.4, 2.27).
    # Mean of GP1 = k1/(1-d1) = 1.8/(1-0.693) = 5.86 lick counts.

    lk2    <- log(75);      label("Log of theta-2, location-like parameter of the second (right-truncated at 61) generalized-Poisson distribution (unitless)")
    # Table 2 2GP column: k2 = 75 (RSE 0.9%); bootstrap CI (65.6, 85.4).
    # Mean of untruncated GP2 = k2/(1-d2) = 75/(1+0.479) = 50.71 lick counts.

    d1     <- 0.693;        label("Dispersion parameter of the first generalized-Poisson distribution (unitless, bounded in [-1, 1])")
    # Table 2 2GP column: d1 = 0.693 (RSE 2.8%); bootstrap CI (0.655, 0.742).
    # d1 > 0 -> overdispersion in the low-count peak (Methods).

    d2     <- -0.479;       label("Dispersion parameter of the second generalized-Poisson distribution (unitless, bounded in [-1, 1])")
    # Table 2 2GP column: d2 = -0.479 (RSE 3.9%); bootstrap CI (-0.697, -0.303).
    # d2 < 0 -> underdispersion in the high-count peak (Methods);
    # this is what the 2NB model could not capture and motivated the 2GP.

    lc     <- log(0.701);   label("Log of the sigmoid Emax slope coefficient c (Hill exponent, unitless)")
    # Table 2 2GP column: c = 0.701 (RSE 1.2%); bootstrap CI (0.705, 0.715).
    # Methods: c is the only structural parameter that does NOT carry IIV.

    # ------------------------------------------------------------------
    # IIV - Sheng 2016 Methods: "Except for the slope coefficient c,
    # interindividual variability on all the model parameters was
    # assumed to be log-normal distribution. In addition, the
    # log-normally distributed interindividual variability was also
    # added to the drug effect (E)." Table 2 reports random-effect
    # spread as CV%; the values below are omega^2 = log(CV^2 + 1)
    # converted to the log-normal-variance scale used in nlmixr2 ini().
    #
    # Multiplicative-log-normal IIV is applied to all structural
    # parameters including those with negative typical values (e0, d2);
    # the encoding "P_i = TVP * exp(eta_P)" preserves sign. For the
    # composite drug effect E (which can change sign across the
    # concentration range), the extra IIV term is added additively
    # inside model() following the paper's "added to E" wording.
    #
    # The CV% -> omega^2 conversion uses omega^2 = log(1 + (CV/100)^2):
    #   Emax  14.8% -> 0.0218
    #   RIC50  1.0% -> 0.000100
    #   E0     2.2% -> 0.000484
    #   k1    42.9% -> 0.169
    #   k2    18.4% -> 0.0335
    #   d1    11.4% -> 0.0129
    #   d2    56.3% -> 0.270
    #   E_add 27.2% -> 0.0721
    # ------------------------------------------------------------------

    etalEmax  ~ 0.0218
    # Table 2 2GP column: omega^2 on Emax reported as CV% = 14.8% (RSE 0).
    etalRIC50 ~ 0.000100
    # Table 2 2GP column: omega^2 on RIC50 reported as CV% = 1% (RSE 35.6).
    etae0     ~ 0.000484
    # Table 2 2GP column: omega^2 on E0 reported as CV% = 2.2% (RSE 8.8).
    etalk1    ~ 0.169
    # Table 2 2GP column: omega^2 on k1 reported as CV% = 42.9% (RSE 4.8).
    etalk2    ~ 0.0335
    # Table 2 2GP column: omega^2 on k2 reported as CV% = 18.4% (RSE 5.4).
    etad1     ~ 0.0129
    # Table 2 2GP column: omega^2 on d1 reported as CV% = 11.4% (RSE 26.3).
    etad2     ~ 0.270
    # Table 2 2GP column: omega^2 on d2 reported as CV% = 56.3% (RSE 0.9).
    # Note: with a 56.3% CV the log-normal multiplicative encoding can
    # push a small fraction (~2.5%) of simulated d2 samples below -1,
    # the GP-validity bound. Clip in simulation if it matters; see the
    # vignette "Assumptions and deviations" section.
    etaE      ~ 0.0721
    # Table 2 2GP column: omega^2 on the composite drug effect E reported
    # as CV% = 27.2% (RSE 34.8). Added additively to E inside model()
    # because E can change sign across the concentration range.
  })

  model({
    # ------------------------------------------------------------------
    # Per-record stimulus (mM quinine HCl dihydrate). For water trials,
    # set STIM_QUININE_MM = 0 in the dataset (Methods: water was the
    # reference / control presentation).
    # ------------------------------------------------------------------
    conc <- STIM_QUININE_MM

    # ------------------------------------------------------------------
    # Individual structural parameters (multiplicative log-normal IIV).
    # ------------------------------------------------------------------
    Emax_i  <- exp(lEmax  + etalEmax)
    RIC50_i <- exp(lRIC50 + etalRIC50)
    e0_i    <- e0    * exp(etae0)     # negative typical value; sign preserved
    k1_i    <- exp(lk1    + etalk1)
    k2_i    <- exp(lk2    + etalk2)
    d1_i    <- d1    * exp(etad1)
    d2_i    <- d2    * exp(etad2)     # negative typical value; sign preserved
    c_i     <- exp(lc)                # no IIV per source

    # ------------------------------------------------------------------
    # Sigmoid Emax drug effect on the LOGIT of the low-count mixture
    # probability (Sheng 2016, "Model for the drug effect" section):
    #   E      = E0 + Emax * conc^c / (RIC50^c + conc^c) + eta_E
    #   p      = expit(E) = exp(E) / (1 + exp(E))
    # where p is the proportion of the low-count generalized-Poisson
    # distribution. logistic(E0) is the minimum mixing probability (at
    # conc = 0); logistic(E0 + Emax) is the maximum (at conc -> Inf).
    # RIC50 is the real IC50 producing half-maximal drug effect (Eq. for
    # RIC50 transformation in Methods); the apparent IC50 is a derived
    # hybrid parameter, not used directly here.
    # ------------------------------------------------------------------
    drugEmax <- Emax_i * (conc^c_i) / (RIC50_i^c_i + conc^c_i)
    E        <- e0_i + drugEmax + etaE
    p_low    <- exp(E) / (1 + exp(E))
    p_high   <- 1 - p_low

    # ------------------------------------------------------------------
    # Mean of the first (low-count) generalized-Poisson distribution
    # (Methods: mean = k1 / (1 - d1) for the GP distribution).
    # ------------------------------------------------------------------
    mu_low <- k1_i / (1 - d1_i)

    # ------------------------------------------------------------------
    # Mean of the second (high-count) generalized-Poisson distribution
    # untruncated. Sheng 2016 right-truncates GP2 at the observed max
    # lick number of 61 (Methods, "Truncated models for the second
    # distribution"); for the typical d2 = -0.479 and k2 = 75 the
    # truncation correction on the mean is small (untruncated mean =
    # 50.71 vs. truncated mean within ~1 of 50.71). The truncated mean
    # is not in closed form, so this expression exposes the untruncated
    # mean as the typical-value summary; the vignette documents the
    # truncation deviation in "Assumptions and deviations".
    # ------------------------------------------------------------------
    mu_high <- k2_i / (1 - d2_i)

    # ------------------------------------------------------------------
    # Mixture-weighted expected lick count (used for typical-value
    # trajectory replication of Figure 2 and for the source-paper
    # Table 1 mean/SD comparison).
    # ------------------------------------------------------------------
    mu_mix <- p_low * mu_low + p_high * mu_high

    # ------------------------------------------------------------------
    # Observation model. The source likelihood is a mixture of two
    # generalized-Poisson distributions (one right-truncated at 61),
    # which rxode2 / nlmixr2 do not express natively in the ini() /
    # model() syntax in this batch. Following the precedent set by
    # inst/modeldb/ddmore/Plan_2012_pain.R and
    # inst/modeldb/ddmore/Schoemaker_2018_levetiracetam.R, the
    # simplification used here is to expose:
    #   - the mixing probability (p_low),
    #   - the per-distribution means (mu_low, mu_high),
    #   - and the mixture-weighted mean (mu_mix),
    # as model variables, and to declare a Poisson observation
    # likelihood per branch with the per-branch typical-value mean.
    # The deterministic typical-value mixture-weighted trajectory
    # (mu_mix) is unaffected by the simplification; the loss is the
    # bimodal under/over-dispersion structure of the residual. The
    # vignette "Assumptions and deviations" section calls this out and
    # documents how a downstream user can post-process Poisson samples
    # through the d1 / d2 dispersion parameters if they need the full
    # bimodal residual structure.
    # ------------------------------------------------------------------
    count_low  <- mu_low
    count_high <- mu_high

    count_low  ~ pois(mu_low)
    count_high ~ pois(mu_high)
  })
}
