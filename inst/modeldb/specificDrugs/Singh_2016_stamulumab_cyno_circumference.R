Singh_2016_stamulumab_cyno_circumference <- function() {
  description <- "Preclinical (cynomolgus monkey). Landmark direct-effect Emax exposure-response model for the percentage growth of left-quadriceps circumference from predose in female cynomolgus monkeys given weekly IV stamulumab (MYO-029), a human IgG1 anti-myostatin (GDF-8) monoclonal antibody, in a 39-week toxicology study (Singh 2016). Two independent fits, at day 135 (week 20) and day 271 (week 39-40): percent change = E0 * (1 + Emax/100 * AUC / (AUC50 + AUC)), where AUC is the steady-state weekly AUC(0-168 h) supplied as the AUC_STAMULUMAB covariate. Typical-value model: no IIV or residual error is reported."
  reference <- "Singh P, Rong H, Gordi T, Bosley J, Bhattacharya I. Translational Pharmacokinetic/Pharmacodynamic Analysis of MYO-029 Antibody for Muscular Dystrophy. Clin Transl Sci. 2016;9(6):302-310. doi:10.1111/cts.12420. PMID 27700008. PMCID PMC5351001. MYO-029 is the development code of the antibody assigned the INN stamulumab."
  vignette <- "Singh_2016_stamulumab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "percent change in left-quadriceps circumference from predose (%)"
  )

  covariateData <- list(
    AUC_STAMULUMAB = list(
      description = "Steady-state serum AUC of stamulumab over the weekly (0-168 h) dosing interval",
      units = "ug*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "0 for vehicle-treated monkeys, which gives percent change = E0.",
        "In the 39-week study the AUC was measured over the day 169 dosing",
        "interval; Figure 3d places the 10, 30 and 100 mg/kg/week groups at",
        "about 0.45e5, 1.4e5 and 2.8e5 ug*h/mL. The cynomolgus PK model",
        "Singh_2016_stamulumab_cyno gives Dose/CL = 43900, 131600 and",
        "438600 ug*h/mL for the same doses."
      ),
      source_name = "AUC0-168hr"
    )
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 16,
    n_studies = 1,
    age_range = NA,
    weight_range = NA,
    sex_female_pct = 100,
    race_ethnicity = NA,
    disease_state = "Healthy female cynomolgus monkeys (39-week toxicology study)",
    dose_range = "Vehicle or 10, 30 and 100 mg/kg IV once weekly for 39 weeks (n = 4 females per group)",
    regions = NA,
    notes = paste(
      "Singh 2016 Methods 'Thirty-nine week monkey study' and Table 3: muscle",
      "circumference measured by tape predose and at weeks 20 (day 135) and",
      "39-40 (day 271) on six muscles in both sexes. Only the female",
      "left-quadriceps cohorts showed a robust exposure-response and were",
      "fitted (Results 'Muscle growth in monkeys'); the fits used group means.",
      "Table 2 'CV%' is the estimation precision of each parameter, not",
      "between-animal variability."
    )
  )

  ini({
    # Singh 2016 Table 2, 'Monkey' rows: 'Female quad-L (135d)'.
    e0_d135 <- 4.1; label("Vehicle percent change in left-quadriceps circumference at day 135 (%)") # Table 2, Female quad-L (135d) E0 = 4.1% (CV 40%)
    emax_d135 <- 274.1; label("Maximum increase over vehicle at day 135 (% of E0)") # Table 2, Female quad-L (135d) Emax = 274.1% (CV 57%)
    lauc50_d135 <- log(32092); label("Steady-state weekly AUC giving half the maximum day-135 effect (ug*h/mL)") # Table 2, Female quad-L (135d) AUC50 = 32,092 ug*h/mL (CV 98%)

    # Singh 2016 Table 2, 'Monkey' rows: 'Female quad-L (271d)'.
    e0_d271 <- 6.1; label("Vehicle percent change in left-quadriceps circumference at day 271 (%)") # Table 2, Female quad-L (271d) E0 = 6.1% (CV 1%)
    emax_d271 <- 200.5; label("Maximum increase over vehicle at day 271 (% of E0)") # Table 2, Female quad-L (271d) Emax = 200.5% (CV 2%)
    lauc50_d271 <- log(32428); label("Steady-state weekly AUC giving half the maximum day-271 effect (ug*h/mL)") # Table 2, Female quad-L (271d) AUC50 = 32428 ug*h/mL (CV 4%)
  })

  model({
    # Table 2 footnote: direct-effect Emax model of steady-state AUC(0-168h)
    # on the rate of growth in muscle circumference; Emax is the '% maximal
    # change over vehicle', so it multiplies the vehicle response E0 (the text
    # quotes Emax as '274% and 201% over vehicle').
    pct_circ_quadl_d135 <- e0_d135 * (1 + emax_d135 / 100 * AUC_STAMULUMAB / (exp(lauc50_d135) + AUC_STAMULUMAB))
    pct_circ_quadl_d271 <- e0_d271 * (1 + emax_d271 / 100 * AUC_STAMULUMAB / (exp(lauc50_d271) + AUC_STAMULUMAB))
  })
}
