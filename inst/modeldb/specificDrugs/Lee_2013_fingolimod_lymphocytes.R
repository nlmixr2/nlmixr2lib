Lee_2013_fingolimod_lymphocytes <- function() {
  description <- paste(
    "Steady-state exposure-biomarker model (Model 1 of Lee 2013) for the",
    "mean peripheral-blood lymphocyte count after month 2 of treatment in",
    "adults with relapsing-remitting multiple sclerosis given placebo,",
    "fingolimod 0.5 or 1.25 mg once daily, or intramuscular interferon",
    "beta-1a once weekly. The count is normally distributed around",
    "rbase * (1 - f), where f is a sigmoid Imax function of the steady-state",
    "pre-dose fingolimod-phosphate concentration in fingolimod-treated",
    "patients, a constant fractional reduction in interferon beta-1a-treated",
    "patients, and 0 on placebo. No time course and no between-subject",
    "random effect: one averaged count per patient. Fitted with SAS PROC",
    "NLMIXED. The predicted count is the input of the companion",
    "relapse-rate model Lee_2013_fingolimod_arr, which the paper chains to",
    "this one to bridge exposure to the annualized relapse rate."
  )
  reference <- paste(
    "Lee JY, Wang Y. Use of a biomarker in exposure-response analysis to",
    "support dose selection for fingolimod. CPT Pharmacometrics Syst",
    "Pharmacol. 2013;2(8):e67. doi:10.1038/psp.2013.44"
  )
  vignette <- "Lee_2013_fingolimod"
  units <- list(
    time = "n/a (static steady-state model; time is not used)",
    dosing = "n/a (exposure enters as the CTROUGH covariate)",
    concentration = "ng/mL (fingolimod-phosphate, CTROUGH)",
    steady_state_lymphocyte_count = "10^9 cells/L"
  )

  covariateData <- list(
    CTROUGH = list(
      description = "Steady-state pre-dose (trough) fingolimod-phosphate concentration in blood",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured (not model-predicted) concentration of the ACTIVE",
        "METABOLITE fingolimod-phosphate (FTY720-P), the exposure measure",
        "named in the Lee 2013 Results; samples were drawn pre-dose at",
        "multiple visits in the two phase III studies (Methods 'Data'), so",
        "the value is a steady-state trough. Total (not unbound)",
        "concentration. Must be 0 for placebo and interferon beta-1a",
        "patients: the sigmoid term then vanishes, which reproduces the",
        "paper's f = 0 for placebo, and the interferon arm is handled by",
        "TRT_IFNB1A. Median steady-state values read from the Lee 2013",
        "Figure 5 boxplots are about 1.25 ng/mL at 0.5 mg and 3.15 ng/mL",
        "at 1.25 mg (dose proportional)."
      ),
      source_name = "CP (supplement SAS code); 'fingolimod-p concentration at SS'"
    ),
    TRT_IFNB1A = list(
      description = "Randomised treatment arm indicator: 1 = intramuscular interferon beta-1a once weekly, 0 = placebo or fingolimod",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo or fingolimod)",
      notes = paste(
        "The active-comparator arm of the 12-month study ('study 2'). In",
        "the supplement SAS code the arm is coded PLNDOSE = 30 (the 30 ug",
        "weekly interferon beta-1a dose); PLNDOSE = 0 is placebo and",
        "PLNDOSE = 0.5 or 1.25 is fingolimod. For TRT_IFNB1A = 1 the",
        "fractional lymphocyte reduction is the constant",
        "e_trt_ifnb1a_rbase instead of the concentration-driven sigmoid."
      ),
      source_name = "PLNDOSE = 30"
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
      "Patients with evaluable PK were pooled (Results). The dependent",
      "variable is each patient's mean of the 3-6 lymphocyte counts drawn",
      "after month 2, by which time the count had reached steady state",
      "(Results). Analysis performed at the US FDA."
    )
  )

  ini({
    # Lee 2013 Table 2, Model 1 ('Lymphocyte-FTY720-p relationship'), and
    # the Model 1 equations. The table mixes scales: it prints the
    # placebo mean, the interferon effect and EC50 on the natural scale,
    # but Emax as the logit-scale estimate (the supplement's ltImax),
    # because Emax was constrained to (0, 1) by a logit transform
    # (Methods 'Model'). expit(1.38) = 0.799 reproduces the
    # 1.76 * (1 - 0.799) = 0.35 x 10^9/L high-concentration plateau of
    # Figure 3, and EC50 = 0.28 ng/mL reproduces its 0.55 x 10^9/L at
    # 0.8 ng/mL (exp(0.28) = 1.32 ng/mL would give 1.33 x 10^9/L there).
    lrbase <- log(1.76)
    label("Mean steady-state lymphocyte count on placebo (10^9 cells/L)") # Table 2 Model 1: 'Mean lymphocyte counts for placebo (beta0)' = 1.76 (SE 0.18)
    e_trt_ifnb1a_rbase <- 0.04
    label("Fractional reduction of the lymphocyte count on interferon beta-1a (unitless)") # Table 2 Model 1: 'Mean lymphocyte counts for interferon beta-1a (beta1)' = 0.04 (SE 0.01); enters as f = beta1 in the Model 1 equations
    logitimax <- 1.38
    label("Logit of the maximum fractional lymphocyte reduction by fingolimod-phosphate (logit)") # Table 2 Model 1: 'Emax' = 1.38 (SE 0.14), the logit-scale ltImax of the supplement code
    lic50 <- log(0.28)
    label("Fingolimod-phosphate concentration at half-maximal lymphocyte reduction (ng/mL)") # Table 2 Model 1: 'EC50' = 0.28 (SE 0.17)
    lhill <- log(1.64)
    label("Hill coefficient of the lymphocyte-reduction sigmoid (unitless)") # Table 2 Model 1: 'Hill coefficient (h)' = 1.64 (SE 0.72)

    # Residual error. Table 2 labels the row sigma^2, but the supplement
    # code declares 'model DV ~ normal(mlym, sigma**2)', so the estimated
    # parameter is the SD, and its simulation step draws the residual as
    # rannor * 0.348, the same value unrounded.
    addSd <- 0.35
    label("Additive residual SD of the steady-state lymphocyte count (10^9 cells/L)") # Table 2 Model 1: 'sigma^2' = 0.35 (SE 0.01), an SD per the supplement code (0.348)
  })

  model({
    rbase <- exp(lrbase)
    imax <- expit(logitimax)
    ic50 <- exp(lic50)
    hill <- exp(lhill)

    # Model 1: f = 0 on placebo (CTROUGH = 0), the sigmoid Imax on
    # fingolimod, and beta1 on interferon beta-1a (CTROUGH = 0 there too).
    f_lymph <- imax * CTROUGH^hill / (ic50^hill + CTROUGH^hill) +
      e_trt_ifnb1a_rbase * TRT_IFNB1A

    steady_state_lymphocyte_count <- rbase * (1 - f_lymph)
    steady_state_lymphocyte_count ~ add(addSd)
  })
}
