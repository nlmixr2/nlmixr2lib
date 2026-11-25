Lee_2013_fingolimod_arr <- function() {
  description <- paste(
    "Biomarker-relapse model (Model 2 of Lee 2013) for the number of",
    "confirmed relapses in adults with relapsing-remitting multiple",
    "sclerosis treated with placebo, fingolimod or interferon beta-1a. The",
    "count over each patient's time on study follows a negative binomial",
    "with mean ARR * t, where t is years on study and the log annualized",
    "relapse rate is log(ARR) = beta0 * (1 - f), with f a steep sigmoid Imax",
    "function of the patient's mean steady-state peripheral lymphocyte",
    "count. ARR rises from 0.20 at very low lymphocyte counts to",
    "0.20^0.5 = 0.45 above about 1.2 x 10^9/L. No between-subject random",
    "effect. Fitted with SAS PROC NLMIXED. rxode2 exposes no",
    "negative-binomial endpoint, so the observation is declared Poisson and",
    "the overdispersion is exposed as the model variable ovdp",
    "(variance = mean * (1 + ovdp * mean)). The lymphocyte input is predicted",
    "from exposure by the companion model Lee_2013_fingolimod_lymphocytes."
  )
  reference <- paste(
    "Lee JY, Wang Y. Use of a biomarker in exposure-response analysis to",
    "support dose selection for fingolimod. CPT Pharmacometrics Syst",
    "Pharmacol. 2013;2(8):e67. doi:10.1038/psp.2013.44"
  )
  vignette <- "Lee_2013_fingolimod"
  units <- list(
    time = "year (time on study; the observation time is the offset t)",
    dosing = "n/a (the drug acts through the LYMPH_ABS covariate)",
    concentration = "n/a",
    relapse_count = "count of confirmed relapses",
    arr = "relapses/year"
  )

  covariateData <- list(
    LYMPH_ABS = list(
      description = "Patient's mean steady-state absolute peripheral-blood lymphocyte count",
      units = "10^9 cells/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reported in 10^9 cells/L (numerically equal to K/uL), NOT the",
        "register's preferred cells/uL: supply 0.9 for 900 cells/uL. In the",
        "fit it is the mean of each patient's 3-6 counts drawn after month",
        "2 of treatment (Lee 2013 Results). For prediction at an unstudied",
        "dose the paper supplies instead the count predicted from",
        "fingolimod-phosphate exposure by Lee_2013_fingolimod_lymphocytes",
        "(Results, 'used as a bridge'); the supplement code adds that",
        "model's residual noise to the prediction and truncates it to",
        "0.11-5.15 x 10^9/L before passing it here. Time-fixed per patient."
      ),
      source_name = "dv (supplement SAS code); 'lymphocytes'"
    )
  )

  population <- list(
    species = "human",
    n_subjects = "approximately 2500 (Methods 'Data'; exact count with evaluable PK not reported)",
    n_studies = 2L,
    age_range = "18-55 years",
    weight_range = "not reported",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "Relapsing-remitting multiple sclerosis (RRMS)",
    dose_range = paste(
      "Fingolimod 0.5 or 1.25 mg orally once daily; placebo; interferon",
      "beta-1a (Avonex) intramuscularly once weekly"
    ),
    regions = "multicenter (regions not reported)",
    notes = paste(
      "Pooled from two double-blind phase III studies: study 1, 24 months,",
      "fingolimod 0.5 and 1.25 mg vs placebo; study 2, 12 months,",
      "fingolimod 0.5 and 1.25 mg vs interferon beta-1a (Methods 'Data').",
      "The dependent variable is each patient's number of confirmed",
      "relapses over the time on study; the annualized relapse rate is the",
      "total relapses divided by total years on study. Analysis performed",
      "at the US FDA."
    )
  )

  ini({
    # Lee 2013 Table 2, Model 2 ('ARR-lymphocyte relationship'), and the
    # Model 2 equations:
    #   Y_i ~ NB(mu_i, k)
    #   log(mu_i) = beta0 * (1 - f_i(lymphocytes)) + log(t_i)
    #   f_i = IEmax * L^h / (EL50^h + L^h)
    # Table 2 prints beta0 back-transformed as 'ARR at lymphocyte count = 0'
    # = exp(beta0); the code estimates beta0 unconstrained. Read this way,
    # ARR is 0.20 below about 0.7 x 10^9/L and 0.20^(1 - 0.5) = 0.45 above
    # about 1.2 x 10^9/L, the two plateaus of Figure 4 (0.203 and 0.441).
    larr0 <- log(0.20)
    label("Log annualized relapse rate at a lymphocyte count of zero (log relapses/year)") # Table 2 Model 2: 'ARR at lymphocyte count = 0' = 0.20 (SE 0.05)
    lic50 <- log(0.9)
    label("Lymphocyte count at half-maximal effect, EL50 (10^9 cells/L)") # Table 2 Model 2: 'EL50' = 0.9 (SE 0.09)
    imax <- 0.5
    label("Maximum fractional change of the log relapse rate, IEmax (unitless)") # Table 2 Model 2: 'IEmax' = 0.5 (SE 0.04)
    lhill <- log(16.6)
    label("Hill coefficient of the lymphocyte-relapse sigmoid (unitless)") # Table 2 Model 2: 'Hill coefficient (h)' = 16.6 (SE 14.9)

    # Negative-binomial dispersion k. The supplement log-likelihood
    #   lgamma(Y + 1/k) - lgamma(Y + 1) - lgamma(1/k) + Y * log(k * mean)
    #     - (Y + 1/k) * log(1 + k * mean)
    # is the NB2 form with variance = mean + k * mean^2, so k is ovdp.
    lovdp <- log(0.7)
    label("Negative-binomial overdispersion k (unitless)") # Table 2 Model 2: 'Dispersion parameter (k)' = 0.7 (SE 0.11)
  })

  model({
    arr0 <- exp(larr0)
    el50 <- exp(lic50)
    hill <- exp(lhill)
    ovdp <- exp(lovdp)

    f_arr <- imax * LYMPH_ABS^hill / (el50^hill + LYMPH_ABS^hill)
    arr <- exp(larr0 * (1 - f_arr))

    # Offset log(t_i): the expected count over t years on study. Place the
    # observation record at the patient's time on study in years.
    lambda <- arr * t
    relapse_count <- lambda
    relapse_count ~ pois(lambda)
  })
}
