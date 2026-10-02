Zhou_2022_alisertib_febrile_neutropenia <- function() {
  description <- paste0(
    "Landmark logistic-regression exposure-SAFETY model for febrile ",
    "neutropenia as a linear function of the average steady-state plasma ",
    "concentration (Css,avg) of the Aurora A kinase inhibitor alisertib, in ",
    "children and adolescents aged 2 to 21 years with advanced malignancies ",
    "(Zhou 2022, n = 146 from ADVL0812 and ADVL0921; P = .01). This is the ",
    "paper's BASE model for the overall population: stepwise selection also ",
    "found cancer type significant, but the authors present the base model ",
    "because only 19 patients had haematological malignancies. There is no ",
    "PK layer and no ODE: Css,avg is supplied as the CSS_ALIS data column in ",
    "umol/L from the companion population PK model Zhou_2022_alisertib. The ",
    "paper prints no regression coefficients; both values are recovered from ",
    "the vector curve of Figure 6B, and they reproduce the two printed ",
    "model-predicted probabilities (0.131 at 2.218 umol/L, 0.122 at ",
    "1.66 umol/L) to within rounding (see the vignette). The companion ",
    "endpoint is Zhou_2022_alisertib_stomatitis."
  )
  reference <- paste(
    "Zhou X, Mould DR, Yuan Y, Fox E, Greengard E, Faller DV,",
    "Venkatakrishnan K. Population pharmacokinetics and exposure-safety",
    "relationships of alisertib in children and adolescents with advanced",
    "malignancies. J Clin Pharmacol. 2022;62(2):206-219.",
    "doi:10.1002/jcph.1958. Exposure-safety methods in Methods",
    "'Exposure-Safety Analysis'; fitted relationship in Figure 6B;",
    "model-predicted probabilities in Results 'Exposure-Safety Analysis'.",
    sep = " "
  )
  vignette <- "Zhou_2022_alisertib"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CSS_ALIS covariate column)",
    concentration = "prob_febrile_neutropenia (probability of febrile neutropenia, 0-1; also logit_febrile_neutropenia)"
  )

  covariateData <- list(
    CSS_ALIS = list(
      description = "Individual average steady-state plasma concentration of alisertib, Dose x 1000 / (CL/F x 518.92 x tau). Supplied as data: this model has no PK layer.",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Zhou 2022 Methods: Css,avg (uM) = Dose (mg) x 1000 / (CL/F (L/h) x",
        "518.92 x tau), with Dose the administered STARTING dose, tau 24 h",
        "(once daily) or 12 h (twice daily), and CL/F the individual post hoc",
        "value from Zhou_2022_alisertib including the tablet relative",
        "bioavailability of 0.671. Enters the logit linearly and uncentred.",
        "Reference exposures printed in Results: geometric mean 2.218 uM at",
        "80 mg/m^2 once daily (tablet) and 1.66 uM at 60 mg/m^2. Figure 6B",
        "plots the analysis range 0.7 to about 24 uM."
      ),
      source_name = "Css,avg"
    )
  )

  covariatesDataExcluded <- list(
    TUMTP_SOLID = list(
      description = "Solid-tumour (versus haematological malignancy) indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Cancer type was significant in the stepwise full model (P < .0001):",
        "the Css,avg relationship was significant in solid tumours (P = .01)",
        "but not in leukaemia (P = .36). The authors present only the base",
        "model, for the overall population, because n = 19 (13 percent) had",
        "haematological malignancies; no coefficients for the cancer-type",
        "model are reported."
      ),
      source_name = "Cancer type"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 146L,
    n_studies = 2L,
    age_range = "2 to 21 years",
    sex_female_pct = 46,
    disease_state = "Relapsed or refractory solid tumours, neuroblastoma, or acute leukaemias (19 haematological)",
    dose_range = "45 to 100 mg/m^2/day alisertib on days 1-7 of 21-day cycles",
    regions = "United States and Canada (Children's Oncology Group sites)",
    notes = paste0(
      "The exposure-safety population is the 146-patient PK population; 32 ",
      "(22 percent) had febrile neutropenia (Supplemental Table S2; the ",
      "Methods list the endpoint as 'grade >= 2 febrile neutropenia')."
    )
  )

  ini({
    # logit(p) = logit_ref + e_css_logit * CSS_ALIS
    #
    # NOT TABULATED IN THE SOURCE. Recovered from the solid curve of
    # Figure 6B, which the publisher PDF stores as vector Bezier
    # segments: the curve, calibrated on the axis tick marks, was
    # sampled and logit(p) regressed on Css,avg (residual SD 0.0098
    # logit over 0.75 to 24.5 uM). Check against the printed Results
    # predictions: 0.131 at 2.218 uM and 0.122 at 1.66 uM (model 0.1313
    # and 0.1228).
    logit_ref <- -2.194; label("Logit of the febrile-neutropenia probability at zero Css,avg (unitless logit)") # recovered from the Figure 6B vector curve
    e_css_logit <- 0.1364; label("Change in the log-odds of febrile neutropenia per 1 umol/L of alisertib Css,avg (1/(umol/L))") # recovered from the Figure 6B vector curve

    # The source likelihood is Bernoulli; this placeholder exists only
    # so rxode2 has an error model to attach. Not a published quantity.
    addSd_prob_febrile_neutropenia <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source
  })

  model({
    logit_febrile_neutropenia <- logit_ref + e_css_logit * CSS_ALIS
    prob_febrile_neutropenia <- expit(logit_febrile_neutropenia)

    prob_febrile_neutropenia ~ add(addSd_prob_febrile_neutropenia)
  })
}
