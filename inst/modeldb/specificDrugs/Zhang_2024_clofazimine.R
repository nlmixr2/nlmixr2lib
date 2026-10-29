Zhang_2024_clofazimine <- function() {
  description <- "Two-compartment population PK model for oral clofazimine in HIV-infected adults with and without cryptosporidiosis-associated diarrhea: first-order absorption with a lag time, linear elimination, relative bioavailability fixed at 1 in participants without diarrhea and reduced by the participant's overall maximum diarrhea grade over the study (grade 1 and grade 2 shifts on log F), inter-individual variability on ka, V/F, Q/F and F, and inter-occasion (per study day) variability on F."
  reference <- paste(
    "Zhang CX, Conrad TM, Hermann D, Gordon MA, Houpt E, Iroh Tam PY, Jere KC,",
    "Nedi W, Operario DJ, Phulusa J, Quinnan GV Jr, Sawyer LA, Barrett LK, Thole H,",
    "Toto N, Van Voorhis WC, Arnold SLM.",
    "Clofazimine pharmacokinetics in HIV-infected adults with diarrhea: Implications",
    "of diarrheal disease on absorption of orally administered therapeutics.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(3):410-423.",
    "doi:10.1002/psp4.13092"
  )
  vignette <- "Zhang_2024_clofazimine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    DIARRHEA_MAX_GRADE_1 = list(
      description = "Overall maximum diarrhea grade over the study is 1 (mild): 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (overall maximum grade 0, no diarrhea during the study)",
      notes = "Subject-level, time-fixed. Stool was graded before each dose on study days 1-5 and once on day 6 (0 = no diarrhea, 1 = mild, 2 = severe); the overall maximum grade is the highest grade over all assessments (Zhang 2024 Methods, Data collection and analysis). Level indicator, mutually exclusive with DIARRHEA_MAX_GRADE_2. Shifts log F by -1.80, F = 0.165 (Table 3). 8 of 23 participants (all in Part A) had overall maximum grade 1 (Table 2).",
      source_name = "Overall max diarrhea grade = 1"
    ),
    DIARRHEA_MAX_GRADE_2 = list(
      description = "Overall maximum diarrhea grade over the study is 2 (severe): 1 = yes, 0 = no",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (overall maximum grade 0, no diarrhea during the study)",
      notes = "Subject-level, time-fixed; see DIARRHEA_MAX_GRADE_1 for the grading. Shifts log F by -3.08, F = 0.0460 (Table 3). Estimated from a single participant with overall maximum grade 2 (Table 2; Discussion).",
      source_name = "Overall max diarrhea grade = 2"
    ),
    OCC = list(
      description = "Occasion index (study day of the dose) for the inter-occasion variability on bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Zhang 2024 defines eight occasions, one per PK visit: visits 1-6 are study days 1-6 and visits 7-8 the two follow-up visits (Methods, Population PK structural model development). Clofazimine was dosed only on days 1-5, and the inter-occasion eta acts on bioavailability, which is evaluated at dose records, so only occasions 1-5 can change a prediction; they are encoded as OCC = 1..5 with a shared variance. Any other OCC value (for example 6-8 or 0) switches the inter-occasion eta off. OCC must be set on the dose records.",
      source_name = "visit"
    )
  )

  covariatesDataExcluded <- list(
    HIV_VL = list(
      description = "HIV viral load at screening",
      units = "copies/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on Q/F and F in the first forward-addition round (Table 1, models 2-3) but no longer significant once overall maximum diarrhea grade was on F (Results, Covariate model). No coefficient reported."
    ),
    CD4 = list(
      description = "CD4+ T-cell count at screening",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on F in the first round (Table 1, model 17); not retained (Results, Covariate model). No coefficient reported."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on V/F and ka (Table 1, models 16, 28, 29); not significant. No coefficient reported."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested on V/F and CL/F (Table 1, models 22, 32, 33); not significant over the narrow 36.6-53 kg range. No coefficient reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "clofazimine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "clofazimine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "clofazimine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 23L,
    n_studies = 1L,
    age_range = "30-64 years",
    age_median = "40.5 (Part A) and 43 (Part B) years",
    weight_range = "36.6-53 kg",
    weight_median = "45.45 (Part A) and 47 (Part B) kg",
    sex_female_pct = 34.8,
    race_ethnicity = "not reported (single site in Blantyre, Malawi)",
    disease_state = "HIV-infected adults on antiretroviral therapy: 12 with Cryptosporidium infection and >= 14 days of diarrhea at enrolment (Part A), 11 matched without cryptosporidiosis or diarrhea (Part B)",
    dose_range = "Oral clofazimine 50 mg (< 50 kg) or 100 mg (>= 50 kg) three times daily (about 5 a.m., 11 a.m. and 6 p.m.) with food for 5 days; 20 of 23 participants received 50 mg",
    regions = "Malawi",
    notes = "Phase IIa randomized, double-blind, placebo-controlled trial (NCT03341767); 428 plasma concentrations from the 23 clofazimine recipients. Overall maximum diarrhea grade 0 / 1 / 2 in 14 / 8 / 1 participants (Zhang 2024 Table 2)."
  )

  ini({
    # Structural parameters (Zhang 2024 Table 3, final model).
    lka <- log(0.625)
    label("Absorption rate constant ka (1/h)") # Table 3 final 'ka (1/h)' = 0.625 (RSE 21.0%)
    ltlag <- log(1.83)
    label("Absorption lag time Tlag (h)") # Table 3 final 'Tlag (h)' = 1.83 (RSE 2.82%)
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F, no diarrhea during the study (fraction)") # Table 3 final 'F' = 1 (fixed)
    lcl <- log(3.71)
    label("Apparent clearance CL/F (L/h)") # Table 3 final 'CL/F (L/h)' = 3.71 (RSE 10.0%)
    lvc <- log(473)
    label("Apparent central volume V/F (L)") # Table 3 final 'V/F (L)' = 473 (RSE 10.7%)
    lq <- log(18.2)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 final 'Q/F (L/h)' = 18.2 (RSE 13.2%)
    lvp <- log(3434)
    label("Apparent peripheral volume V2/F (L)") # Table 3 final 'V2/F (L)' = 3434 (RSE 8.07%)

    # Covariate effects on log F (Equation 6).
    e_diarrhea_max_grade_1_fdepot <- -1.80
    label("Shift in log F for overall maximum diarrhea grade 1 (unitless)") # Table 3 'Overall max diarrhea grade = 1' = -1.80 (RSE 11.2%); Results: Frel 0.165
    e_diarrhea_max_grade_2_fdepot <- -3.08
    label("Shift in log F for overall maximum diarrhea grade 2 (unitless)") # Table 3 'Overall max diarrhea grade = 2' = -3.08 (RSE 3.55%); Results: Frel 0.0460

    # Inter-individual variability. Table 3 reports %CV; variances are
    # omega^2 = log(CV^2 + 1).
    etalka ~ 1.169412 # Table 3 final 'IIV on ka (%)' = 149
    etalvc ~ 0.001196444 # Table 3 final 'IIV on V/F (%)' = 3.46
    etalq ~ 0.1844027 # Table 3 final 'IIV on Q/F (%)' = 45.0
    etalfdepot ~ 0.1315764 # Table 3 final 'IIV on F (%)' = 37.5

    # Inter-occasion variability on F, one eta per dosing day with a shared
    # variance (occasions 2-5 fixed to the occasion-1 value).
    etaiov_fdepot_1 ~ 0.1270061 # Table 3 final 'Interoccasion variability on F (%)' = 36.8
    etaiov_fdepot_2 ~ fixed(0.1270061) # same variance as occasion 1
    etaiov_fdepot_3 ~ fixed(0.1270061) # same variance as occasion 1
    etaiov_fdepot_4 ~ fixed(0.1270061) # same variance as occasion 1
    etaiov_fdepot_5 ~ fixed(0.1270061) # same variance as occasion 1

    # Residual error: combined proportional and additive (Equation 1). Table 3
    # prints a single sigma; in a Phoenix additive + multiplicative model that
    # is the additive SD, in the mg/L units of Figure 2. The proportional SD is
    # not reported and is fixed to 0.
    propSd <- fixed(0)
    label("Proportional residual error (fraction); not reported") # Equation 1 term Cpred x eps1; no value in Table 3
    addSd <- 0.00268
    label("Additive residual error (mg/L)") # Table 3 final 'sigma' = 0.00268 (RSE 42.0%)
  })

  model({
    # Occasion (dosing-day) indicators for the inter-occasion variability.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5

    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)

    # Equation 6: F = 1 x exp(eta_F + theta1 (grade = 1) + theta2 (grade = 2)).
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot +
      e_diarrhea_max_grade_1_fdepot * DIARRHEA_MAX_GRADE_1 +
      e_diarrhea_max_grade_2_fdepot * DIARRHEA_MAX_GRADE_2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
