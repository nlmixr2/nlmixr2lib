Zhou_2022_alisertib_stomatitis <- function() {
  description <- paste0(
    "Landmark logistic-regression exposure-SAFETY model for CTCAE grade 2 or ",
    "higher stomatitis as a linear function of the average steady-state ",
    "plasma concentration (Css,avg) of the Aurora A kinase inhibitor ",
    "alisertib, in children and adolescents aged 2 to 21 years with advanced ",
    "malignancies (Zhou 2022, n = 146 from ADVL0812 and ADVL0921; P = .0002). ",
    "There is no PK layer and no ODE: Css,avg is supplied as the CSS_ALIS data ",
    "column in umol/L, computed in the source as starting dose / (individual ",
    "CL/F x dosing interval) from the companion population PK model ",
    "Zhou_2022_alisertib. No covariate was retained. The paper prints no ",
    "regression coefficients; the two values below are back-solved exactly ",
    "from the two model-predicted probabilities the Results print (0.159 at ",
    "2.218 umol/L and 0.143 at 1.66 umol/L), and the Figure 6A curve confirms ",
    "the logit-linear-in-Css form (see the vignette). The companion endpoint ",
    "is Zhou_2022_alisertib_febrile_neutropenia."
  )
  reference <- paste(
    "Zhou X, Mould DR, Yuan Y, Fox E, Greengard E, Faller DV,",
    "Venkatakrishnan K. Population pharmacokinetics and exposure-safety",
    "relationships of alisertib in children and adolescents with advanced",
    "malignancies. J Clin Pharmacol. 2022;62(2):206-219.",
    "doi:10.1002/jcph.1958. Exposure-safety methods in Methods",
    "'Exposure-Safety Analysis'; fitted relationship in Figure 6A;",
    "model-predicted probabilities in Results 'Exposure-Safety Analysis'.",
    sep = " "
  )
  vignette <- "Zhou_2022_alisertib"
  units <- list(
    time = "n/a (static landmark exposure-response regression; no time dimension)",
    dosing = "n/a (no dose events; exposure enters as the CSS_ALIS covariate column)",
    concentration = "prob_stomatitis_grade2 (probability of grade 2 or higher stomatitis, 0-1; also logit_stomatitis_grade2)"
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
        "80 mg/m^2 once daily (tablet) and 1.66 uM at 60 mg/m^2. Figure 6A",
        "plots the analysis range 0.7 to about 24 uM."
      ),
      source_name = "Css,avg"
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
      "The exposure-safety population is the 146-patient PK population; 47 ",
      "(32 percent) had grade 2 or higher stomatitis (Supplemental Table S2). ",
      "Age group, sex, performance score, number of prior regimens, cancer ",
      "type and number of cycles were screened in a stepwise full model and ",
      "none was retained."
    )
  )

  ini({
    # logit(p) = logit_ref + e_css_logit * CSS_ALIS
    #
    # NOT TABULATED IN THE SOURCE. Back-solved from the two printed
    # model predictions (Results): p = 0.159 at Css,avg 2.218 uM and
    # p = 0.143 at 1.66 uM. slope = (logit(0.159) - logit(0.143)) /
    # (2.218 - 1.66); intercept = logit(0.159) - slope * 2.218. The
    # printed rounding of the probabilities bounds the slope to
    # 0.21-0.24 per uM. The vector curve of Figure 6A, extracted
    # from the publisher PDF, is exactly logit-linear in Css,avg
    # (residual SD 0.011 logit) with slope 0.207 and intercept -2.07;
    # it sits about 0.008 above the printed probabilities at the
    # clinical exposures. The vignette shows both.
    logit_ref <- -2.162; label("Logit of the grade 2 or higher stomatitis probability at zero Css,avg (unitless logit)") # back-solved from Results text predictions 0.159 at 2.218 uM and 0.143 at 1.66 uM
    e_css_logit <- 0.2238; label("Change in the log-odds of grade 2 or higher stomatitis per 1 umol/L of alisertib Css,avg (1/(umol/L))") # back-solved from Results text predictions 0.159 at 2.218 uM and 0.143 at 1.66 uM

    # The source likelihood is Bernoulli; this placeholder exists only
    # so rxode2 has an error model to attach. Not a published quantity.
    addSd_prob_stomatitis_grade2 <- fixed(0.001); label("Placeholder additive residual SD on the typical-value event probability; the source likelihood is Bernoulli (no source residual)") # not from source
  })

  model({
    logit_stomatitis_grade2 <- logit_ref + e_css_logit * CSS_ALIS
    prob_stomatitis_grade2 <- expit(logit_stomatitis_grade2)

    prob_stomatitis_grade2 ~ add(addSd_prob_stomatitis_grade2)
  })
}
