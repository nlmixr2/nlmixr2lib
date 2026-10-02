Chan_2021_pregabalin_lsr28 <- function() {
  description <- paste(
    "Exposure-response (Emax) model for the natural log-transformed",
    "28-day seizure rate (LSR28) during the 12-week double-blind",
    "treatment phase in pediatric (4-16 years) and adult patients with",
    "focal onset seizures taking adjunctive pregabalin (Chan 2021).",
    "LSR28 = Intercept - (Intercept - Emax) * Cav,ss / (EC50 + Cav,ss) +",
    "Slope_baseline * baseline LSR28 (Appendix Equations II), with a",
    "common Emax (-0.924, the asymptotic intercept under maximal drug",
    "effect) and EC50 (4.69 ug/mL) and population-specific intercepts",
    "and baseline slopes for children and adults (Table 3). The drug",
    "effect therefore differs between populations only through the",
    "intercept, i.e. the placebo response. Fitted by nonlinear least",
    "squares; no between-subject or residual variance was reported, so",
    "the model returns the typical prediction only. Exposure enters as",
    "the per-patient column CAV, the individual predicted average",
    "steady-state concentration from modellib('Chan_2021_pregabalin')."
  )
  reference <- paste(
    "Chan PLS, Marshall SF, McFadyen L, Liu J.",
    "Pregabalin Population Pharmacokinetic and Exposure-Response Analyses",
    "for Focal Onset Seizures in Children (4-16 years) and Adults, to",
    "Support Dose Recommendations in Children.",
    "Clin Pharmacol Ther. 2021;110(1):132-140. doi:10.1002/cpt.2132.",
    "Model equations in Supplementary Information Appendix Equations II;",
    "parameter estimates in Table 3 (row 'Children (C) + Adult (A)')."
  )
  vignette <- "Chan_2021_pregabalin"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    CAV = list(
      description = "Average steady-state pregabalin plasma concentration (Cav,ss)",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual predicted Cav,ss from the final population PK model",
        "(modellib('Chan_2021_pregabalin')) for the patient's",
        "double-blind-phase dose, i.e. daily dose / (24 * individual",
        "CL/F). Set to 0 for placebo patients."
      ),
      source_name = "Cav,ss"
    ),
    LSR28_BL = list(
      description = "Observed baseline natural log-transformed 28-day seizure rate",
      units = "log(seizures per 28 days)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Natural log of the baseline-period 28-day rate of all focal",
        "onset (partial onset) seizures, observed per patient. Enters",
        "linearly through a population-specific slope. Medians 3.00",
        "(children) and 2.40 (adults) (Table 3 footnote a)."
      ),
      source_name = "Baseline"
    ),
    CHILD = list(
      description = "Pediatric population indicator (1 = child, 0 = adult)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "1 for pediatric patients aged 4-16 years (the PERIWINKLE study,",
        "A0081041); 0 for patients from the three adult studies. The adult",
        "population included 8 adolescents aged 13-16 years from study",
        "1008-034, who were analysed with the adult parameters (Table 3",
        "footnote b), so CHILD = 0 for them."
      ),
      source_name = "C / A population"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1138L,
    n_studies = 4L,
    age_range = "4-82 years",
    age_median = "10 years (children); 38 years (adults)",
    weight_range = "11-180 kg",
    weight_median = "35.6 kg (children); 74.9 kg (adults)",
    sex_female_pct = 48.9,
    race_ethnicity = c(White = 82.2, Black = 4.2, Asian = 8.3, Other = 5.4),
    disease_state = "Focal onset seizures, adjunctive pregabalin or placebo",
    dose_range = paste(
      "Adults 150-600 mg/day; children 2.5 or 10 mg/kg/day (>= 30 kg) or",
      "3.5 or 14 mg/kg/day (< 30 kg), b.i.d.; placebo arms included"
    ),
    regions = "USA 61.8%, European Union 21.4%, Asia-Pacific 1.5%, other 15.3%",
    notes = paste(
      "280 pediatric patients from PERIWINKLE (A0081041) and 858 adult",
      "patients (including 8 adolescents aged 13-16 years) from three",
      "adult phase III studies (Chan 2021 Table S2). Concomitant",
      "antiepileptic drugs: 1 (28.2%), 2 (46.2%), 3 (25.0%), 4 (1.3%)."
    )
  )

  ini({
    # Chan 2021 Table 3, row 'Children (C) + Adult (A)' (final model:
    # common Emax and EC50, separate baseline and placebo effects). Values
    # are estimate +/- standard error. The adult-only row is the
    # model-development step and is not encoded.
    e0 <- 0.110
    label("Intercept of LSR28 for adults (log seizures per 28 days)") # Table 3 'A: 0.110 +/- 0.066'
    e0_child <- -0.409
    label("Intercept of LSR28 for children (log seizures per 28 days)") # Table 3 'C: -0.409 +/- 0.106'
    e_lsr28_bl <- 0.945
    label("Slope of on-treatment LSR28 on baseline LSR28 for adults (unitless)") # Table 3 'Slope baseline A: 0.945 +/- 0.022'
    e_lsr28_bl_child <- 1.03
    label("Slope of on-treatment LSR28 on baseline LSR28 for children (unitless)") # Table 3 'Slope baseline C: 1.03 +/- 0.026'
    rmax_lsr28 <- -0.924
    label("Asymptotic intercept of LSR28 under maximal drug effect, both populations (log seizures per 28 days)") # Table 3 'Emax -0.924 +/- 0.214'; Appendix Equations II 'Emax ... (Intercept + maximum drug effect)'
    lec50 <- log(4.69)
    label("Cav,ss producing half the maximal drug effect, both populations (ug/mL)") # Table 3 'EC50, ug/mL 4.69 +/- 2.17'
  })

  model({
    # Population-specific intercept and baseline slope (Table 3)
    int_lsr28 <- e0 * (1 - CHILD) + e0_child * CHILD
    slope_bl <- e_lsr28_bl * (1 - CHILD) + e_lsr28_bl_child * CHILD
    ec50 <- exp(lec50)

    # Appendix Equations II, 'Emax treatment effect response model' with
    # linear baseline effect:
    # Response_i = Intercept - (Intercept - Emax) * Cav_i / (EC50 + Cav_i)
    #              + Slope_baseline * Baseline_i
    lsr28 <- int_lsr28 -
      (int_lsr28 - rmax_lsr28) * CAV / (ec50 + CAV) +
      slope_bl * LSR28_BL

    # Change from baseline in LSR28 (the quantity plotted in Figure 2)
    dlsr28 <- lsr28 - LSR28_BL
  })
}
