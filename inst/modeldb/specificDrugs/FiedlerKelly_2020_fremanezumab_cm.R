FiedlerKelly_2020_fremanezumab_cm <- function() {
  description <- "Population PD exposure-response model relating fremanezumab average plasma concentration (Cav) to monthly moderate-to-severe headache days in adults with chronic migraine. Placebo time-course is a Hill (sigmoid) function in months and the drug effect is a power function of Cav centered on the population median Cav. Fitted to 5312 monthly observations from 1361 chronic-migraine patients pooled across the LBR-101-021 phase 2b and TV48125-CNS-30049 phase 3 studies (Fiedler-Kelly 2020)."
  reference <- "Fiedler-Kelly JB, Passarell J, Ludwig E, Levi M, Cohen-Barak O. Effect of Fremanezumab Monthly and Quarterly Doses on Efficacy Responses. Headache. 2020 Jul;60(7):1376-1391. doi:10.1111/head.13855. PMID: 32445498."
  vignette <- "FiedlerKelly_2020_fremanezumab_cm"
  units <- list(time = "month", dosing = "n/a (PD-only model; no dose events)", concentration = "ug/mL", response = "moderate-to-severe headache days/month")

  covariateData <- list(
    CAV = list(
      description        = "Average fremanezumab plasma concentration over the dosing interval used as the exposure metric in the exposure-response model",
      units              = "ug/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed within each 28-day dosing-interval aggregation. In Fiedler-Kelly 2020 the per-subject Cav values are derived from individual empirical-Bayes PK estimates from the previously-published Fiedler-Kelly 2019 population PK model (Fiedler-Kelly_2019_fremanezumab in this library). Set to 0 for placebo periods.",
      source_name        = "CAV"
    ),
    ACUTE_MED_DAYS = list(
      description        = "Baseline number of days/month of acute migraine medication use (triptans / ergot compounds)",
      units              = "days/month",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (counted during the 28-day run-in period prior to first dose). Enters as a piecewise-linear shift on baseline M/S headache days with breakpoint at 5 d/mo: contributes 0 for ACUTE_MED_DAYS <= 5 and slope_AM * (ACUTE_MED_DAYS - 5) for ACUTE_MED_DAYS > 5 (Fiedler-Kelly 2020 Results, Monthly Headache Days of at Least Moderate Severity in Patients With CM).",
      source_name        = "Baseline days/month of acute medications"
    )
  )

  population <- list(
    n_subjects        = 1361L,
    n_observations    = 5312L,
    n_studies         = 2L,
    age_range         = "18-71 years",
    age_median        = "42 years",
    weight_range      = "43.5-131.8 kg",
    weight_median     = "70.80 kg",
    sex_female_pct    = 87.2,
    race_ethnicity    = c(White = 79.8, Black = 9.0, Asian = 9.0, AmIndAlaskaNative = 0.4, NativeHawaiianPacIsl = 0.2, Other = 1.5),
    ethnicity_hispanic_pct = 9.1,
    disease_state     = "Adults with chronic migraine (headaches on at least 15 days/month with at least 8 migraine days/month per ICHD-3 criteria).",
    dose_range        = "Fremanezumab 225 mg monthly with a 675 mg starting dose, 675 mg quarterly, 900 mg monthly (phase 2b only), or placebo SC for 3 months.",
    regions           = "Multinational (LBR-101-021 phase 2b and TV48125-CNS-30049 phase 3 chronic-migraine studies).",
    notes             = "Demographics from Supplementary Table S2 of Fiedler-Kelly 2020. Concomitant analgesic-medication use 10.5% and concomitant migraine-preventive medication use 24.2% across the pooled CM cohort. Observation unit is one 28-day month."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final exposure-response parameters for monthly moderate-to-severe
    # (M/S) headache days in chronic migraine. Values from Fiedler-Kelly
    # 2020 Supplementary Table S4. Time-axis is in MONTHS (28-day periods)
    # to match the paper's data aggregation.
    # ----------------------------------------------------------------------
    bl_cm        <- 10.2;                  label("Typical baseline M/S headache days/month at acute-med use <= 5 d/mo")  # Fiedler-Kelly 2020 Table S4
    slope_AM     <- 0.460;                 label("Slope on baseline acute-med days >5 d/mo (M/S days / day acute)")       # Fiedler-Kelly 2020 Table S4
    maxPLC_cm    <- fixed(-6.24);          label("Maximum placebo response in M/S headache days (negative; days/month)")  # Fiedler-Kelly 2020 Table S4 (FIXED)
    T50_PLC      <- fixed(1.76);           label("Time to half-maximum placebo response (months, FIXED)")                 # Fiedler-Kelly 2020 Table S4 (FIXED, NE)
    lhill_PLC    <- fixed(log(0.486));     label("Log Hill coefficient for placebo time-course (unitless)")               # Fiedler-Kelly 2020 Table S4 (FIXED)
    # Logit transform of typical drug-effect intercept (typical 0.157):
    # logit(0.157) = log(0.157 / (1 - 0.157)) = -1.6816. IIV is logit-normal
    # so that the individual fractional reduction stays in (0, 1).
    logitDrugInt <- log(0.157 / (1 - 0.157)); label("Logit of drug-effect intercept (logit-transformed fraction, unitless): fractional reduction from baseline at median Cav") # Fiedler-Kelly 2020 Table S4 (typical 0.157)
    ldrugExp     <- log(0.328);            label("Log exponent for fremanezumab Cav effect (unitless)")                   # Fiedler-Kelly 2020 Table S4

    # ----------------------------------------------------------------------
    # Inter-individual variability (Fiedler-Kelly 2020 Table S4 with
    # footnotes a and b explaining the variance-to-CV% conversions):
    #   a: %CV = 100 * sqrt(omega^2)              (log-normal IIV)
    #   b: %CV = 100 * sqrt(omega^2) * (1 - typ)  (logit-normal on bounded fraction)
    # bl_cm and slope_AM share an additive eta (S4 reports identical
    # magnitude 4.69 SD on both rows), implemented as a single additive
    # eta on the composite individual baseline expression.
    # maxPLC_cm carries an additive eta (SD 6.66) per S4 row 3.
    # T50_PLC has no IIV (NE).
    # ----------------------------------------------------------------------
    etabl_cm        ~ 21.9961  # SD 4.69, variance 4.69^2 = 21.9961; Fiedler-Kelly 2020 Table S4
    etamaxPLC_cm    ~ 44.3556  # SD 6.66, variance 6.66^2 = 44.3556; Fiedler-Kelly 2020 Table S4
    etalhill_PLC    ~ 1.69     # log-normal omega^2 (130 %CV per footnote a: sqrt(1.69) = 1.30); Fiedler-Kelly 2020 Table S4
    etalogitDrugInt ~ 1.21     # logit-normal omega^2 (92.6 %CV per footnote b: 100*sqrt(1.21)*(1-0.157)); Fiedler-Kelly 2020 Table S4
    etaldrugExp     ~ 0.945    # log-normal omega^2 (97.2 %CV per footnote a: sqrt(0.945) = 0.972); Fiedler-Kelly 2020 Table S4

    # Residual error - additive on monthly M/S headache days.
    # Source reports variance 7.09 with SD column 2.66 (sqrt(7.09) = 2.663).
    # The output-prefixed name follows the multi-output residual-error
    # convention (`<output>addSd`).
    msHeadacheDaysaddSd <- 2.66; label("Additive residual error on monthly M/S headache days")  # Fiedler-Kelly 2020 Table S4 (SD of variance 7.09)
  })

  model({
    # --------------------------------------------------------------
    # Population-median Cav centering value used in the drug-effect
    # power function. Fiedler-Kelly 2020 reports the drug intercept
    # as the "fractional reduction at median Cav" but does not list
    # the median Cav numerically in Supplementary Table S4 nor in the
    # main text. The value 69 ug/mL was visually inferred from the
    # x-axis position of the per-regimen median markers in
    # Figure 2B and reproduces the narrative drug-effect ranges
    # (12-16% across Cav 28-70 ug/mL; ~18% at Cav 120 ug/mL).
    # Documented in the validation vignette under Assumptions and
    # deviations.
    # --------------------------------------------------------------
    CavMedian <- 69

    BL_i      <- bl_cm + slope_AM * max(0, ACUTE_MED_DAYS - 5) + etabl_cm
    maxPLC_i  <- maxPLC_cm + etamaxPLC_cm
    hill_i    <- exp(lhill_PLC + etalhill_PLC)
    drugInt_i <- 1 / (1 + exp(-(logitDrugInt + etalogitDrugInt)))
    drugExp_i <- exp(ldrugExp + etaldrugExp)

    # Placebo Hill function in time (Fiedler-Kelly 2020 Figure 2B form):
    #   placebo_eff = maxPLC * t^Hill / (T50^Hill + t^Hill)
    # maxPLC_i is negative, so adding to BL gives a reduction over time.
    placebo_eff <- maxPLC_i * (time^hill_i) / (T50_PLC^hill_i + time^hill_i)

    # Drug effect (multiplicative on individual baseline; power function of Cav).
    # CAV = 0 in placebo periods => 0^drugExp_i = 0 => drug_eff = 0.
    drug_eff <- BL_i * drugInt_i * (CAV / CavMedian)^drugExp_i

    # Predicted monthly M/S headache days = baseline + placebo time effect - drug effect.
    msHeadacheDays <- BL_i + placebo_eff - drug_eff

    msHeadacheDays ~ add(msHeadacheDaysaddSd)
  })
}
