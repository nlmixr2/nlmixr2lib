FiedlerKelly_2020_fremanezumab_em <- function() {
  description <- "Population PD exposure-response model relating fremanezumab average plasma concentration (Cav) to monthly migraine days in adults with episodic migraine. Placebo time-course is an exponential growth in months (predicted reduction = exp(exponent * t)) and the drug effect is an Emax/EC50 of Cav scaled by individual baseline migraine days. Fitted to 4444 monthly observations from 1142 episodic-migraine patients pooled across the LBR-101-022 phase 2b and TV48125-CNS-30050 phase 3 studies (Fiedler-Kelly 2020)."
  reference <- "Fiedler-Kelly JB, Passarell J, Ludwig E, Levi M, Cohen-Barak O. Effect of Fremanezumab Monthly and Quarterly Doses on Efficacy Responses. Headache. 2020 Jul;60(7):1376-1391. doi:10.1111/head.13855. PMID: 32445498."
  vignette <- "FiedlerKelly_2020_fremanezumab_em"
  units <- list(time = "month", dosing = "n/a (PD-only model; no dose events)", concentration = "ug/mL", response = "migraine days/month")

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
      notes              = "Time-fixed per subject (counted during the 28-day run-in period prior to first dose). Enters as a piecewise-linear shift on baseline migraine days with breakpoint at 5 d/mo: contributes 0 for ACUTE_MED_DAYS <= 5 and slope_AM * (ACUTE_MED_DAYS - 5) for ACUTE_MED_DAYS > 5 (Fiedler-Kelly 2020 Results, Monthly Migraine Days in Patients With EM).",
      source_name        = "Baseline days/month of acute medications"
    )
  )

  population <- list(
    n_subjects        = 1142L,
    n_observations    = 4444L,
    n_studies         = 2L,
    age_range         = "18-70 years",
    age_median        = "42 years",
    weight_range      = "43.1-120.0 kg",
    weight_median     = "72.26 kg",
    sex_female_pct    = 85.7,
    race_ethnicity    = c(White = 79.9, Black = 11.2, Asian = 7.0, AmIndAlaskaNative = 0.6, NativeHawaiianPacIsl = 0.1, Other = 1.2),
    ethnicity_hispanic_pct = 12.3,
    disease_state     = "Adults with episodic migraine (headaches on 6-14 days/month with at least 4 migraine days/month per ICHD-3 criteria; phase 2b cohort required 8-14 days/month).",
    dose_range        = "Fremanezumab 225 mg monthly, 675 mg monthly, 675 mg quarterly, or 225 mg monthly with a 675 mg starting dose, all SC for 3 months; placebo arms pooled.",
    regions           = "Multinational (LBR-101-022 phase 2b and TV48125-CNS-30050 phase 3 episodic-migraine studies).",
    notes             = "Demographics from Supplementary Table S1 of Fiedler-Kelly 2020. Concomitant analgesic-medication use 6.7% and concomitant migraine-preventive medication use 21.4% across the pooled EM cohort. Observation unit is one 28-day month."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final exposure-response parameters for monthly migraine days in
    # episodic migraine. Values from Fiedler-Kelly 2020 Supplementary
    # Table S3 (final parameter estimates and precisions). Time-axis is in
    # MONTHS (28-day periods) to match the paper's data aggregation; the
    # `units` field above sets the documented time unit.
    # ----------------------------------------------------------------------
    bl_em       <- 8.35;          label("Typical baseline migraine days/month at acute-med use <= 5 d/mo")     # Fiedler-Kelly 2020 Table S3
    slope_AM    <- 0.438;         label("Slope on baseline acute-med days >5 d/mo (days migraine / day acute)") # Fiedler-Kelly 2020 Table S3
    exp_PLC     <- fixed(0.360);  label("Exponent for placebo time-course (1/month, FIXED in source)")          # Fiedler-Kelly 2020 Table S3
    # Logit transform of typical maximum fractional Cav response (typical 0.252):
    # logit(0.252) = log(0.252 / (1 - 0.252)) = -1.0883. IIV is logit-normal
    # so that individual Emax stays in (0, 1) under any random effect.
    logitEmax   <- log(0.252 / (1 - 0.252));  label("Logit of maximum fractional response in migraine days due to Cav (logit-transformed fraction, unitless)") # Fiedler-Kelly 2020 Table S3 (typical 0.252)
    EC50_drug   <- 3.60;          label("Fremanezumab Cav at half-maximum response (ug/mL)")                    # Fiedler-Kelly 2020 Table S3 (NE for IIV)

    # ----------------------------------------------------------------------
    # Inter-individual variability. Magnitude and form per Fiedler-Kelly
    # 2020 Table S3 with footnote a explaining the logit-normal CV%
    # calculation for the maximum fractional Cav response:
    #   43.3 %CV = 100 * sqrt(0.335) * (1 - 0.252)  =>  omega^2 = 0.335
    # bl_em and slope_AM share the same eta (S3 reports identical magnitude
    # 1.61 SD on both rows), implemented here as a single additive eta on
    # the composite individual baseline expression.
    # ----------------------------------------------------------------------
    etabl_em      ~ 2.5921    # SD 1.61, variance 1.61^2 = 2.5921; Fiedler-Kelly 2020 Table S3
    etaexp_PLC    ~ 8.5264    # SD 2.92, variance 2.92^2 = 8.5264; Fiedler-Kelly 2020 Table S3
    etalogitEmax  ~ 0.335     # logit-normal omega^2; Fiedler-Kelly 2020 Table S3 footnote a

    # ----------------------------------------------------------------------
    # Residual error - additive on monthly migraine days. Source reports
    # variance 5.52 (days^2) with SD column 2.35 (sqrt(5.52) = 2.349).
    # nlmixr2 expects SD, so the SD form is used directly here. The
    # output-prefixed name follows the multi-output residual-error
    # convention (`addSd_<output>`).
    # ----------------------------------------------------------------------
    addSd_migraineDays <- 2.35; label("Additive residual error on monthly migraine days")  # Fiedler-Kelly 2020 Table S3 (SD of variance 5.52)
  })

  model({
    # Individual baseline migraine days with piecewise-linear acute-med-days shift
    # (breakpoint 5 d/mo) and shared additive eta on the composite baseline.
    BL_i <- bl_em + slope_AM * max(0, ACUTE_MED_DAYS - 5) + etabl_em

    # Individual placebo-time-course exponent.
    exp_i <- exp_PLC + etaexp_PLC

    # Individual maximum fractional Cav response (logit-normal IIV).
    Emax_i <- 1 / (1 + exp(-(logitEmax + etalogitEmax)))

    # Placebo time-course (Fiedler-Kelly 2020 Figure 2A form): predicted
    # reduction from baseline = exp(exponent * month). At month 0 the
    # reduction is 1 day (model anchors near, not exactly at, BL); at
    # month 3 the reduction is exp(0.360*3) = 2.94 d (matches paper's
    # narrative "approximately 3 days" placebo reduction at month 3).
    placebo_red <- exp(exp_i * time)

    # Drug-effect reduction (Emax/EC50 of Cav, multiplied by individual baseline).
    # CAV = 0 in placebo periods => drug_red = 0.
    drug_red <- BL_i * Emax_i * CAV / (EC50_drug + CAV)

    # Predicted monthly migraine days = baseline minus placebo time effect minus drug effect.
    migraineDays <- BL_i - placebo_red - drug_red

    migraineDays ~ add(addSd_migraineDays)
  })
}
