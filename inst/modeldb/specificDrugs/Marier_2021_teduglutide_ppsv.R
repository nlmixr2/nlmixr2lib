Marier_2021_teduglutide_ppsv <- function() {
  description <- "Time- and exposure-response model for the change from baseline in weekly prescribed parenteral support volume (PPSV) in adult and pediatric patients with short bowel syndrome treated with the GLP-2 analog teduglutide (Marier 2021). The change approaches a maximum reduction Emax with a hyperbolic time course (ET50 = 168 days, fixed); Emax scales as a power function of the individual steady-state teduglutide Cmax, with separate multipliers on Emax for the placebo arm and the 0.10 mg/kg dose level. Algebraic model with no dose events: the steady-state Cmax is supplied as a covariate, for example from Marier_2021_teduglutide."
  reference <- "Marier JF, Jomphe C, Peyret T, Wang Y. Population pharmacokinetics and exposure-response analyses of teduglutide in adult and pediatric patients with short bowel syndrome. Clin Transl Sci. 2021;14(6):2497-2509. doi:10.1111/cts.13117"
  vignette <- "Marier_2021_teduglutide"
  units <- list(
    time = "day (time since the start of randomized treatment)",
    dosing = "n/a (no dose events; exposure enters as the covariate CMAX)",
    concentration = "ppsvcfb (change from baseline in prescribed parenteral support volume, L/week)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual steady-state maximum plasma teduglutide concentration (Cmax,ss) at the subject's randomized dose",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives Emax as (CMAX / 30)^0.684 in actively treated subjects below the 0.10 mg/kg dose level (Table S10 'Exponent for the effect of Cmax on Emax'). In the source analysis it is the individual Cmax,ss derived from the final popPK model (Methods, 'Population PK analysis of teduglutide': 'Rich concentration-time profiles were simulated with the final PopPK model to derive exposure metrics'; 'Exposure-response analysis': 'teduglutide steady-state exposure (Cmax)'). The normalizing concentration 30 ng/mL is not printed with Table S10; it is taken from the authors' raw-PPSV sensitivity model of the same analysis (Table S14, '(Cmax/30)^0.507'). Not used for PLACEBO = 1 or DOSE_HIGH = 1 subjects (supply any finite value, e.g. 0).",
      source_name = "Cmax"
    ),
    PLACEBO = list(
      description = "Placebo-arm indicator: 1 = placebo with standard of care, 0 = teduglutide",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (teduglutide)",
      notes = "Placebo subjects take Emax x (1 - 0.276) in place of the Cmax term (Table S10 'Placebo effect on Emax' = -0.276, a proportional shift of the 1 + theta form used throughout the analysis).",
      source_name = "dose = placebo"
    ),
    DOSE_HIGH = list(
      description = "Highest dose-level indicator of the exposure-response data: 1 = teduglutide 0.10 mg/kg/day, 0 = any lower dose level (0.0125, 0.025 or 0.05 mg/kg/day) or placebo",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (dose below 0.10 mg/kg/day)",
      notes = "Subjects at 0.10 mg/kg/day take Emax x (1 - 0.225) in place of the Cmax term (Table S10 'Dose effect (0.10 mg/kg) on Emax' = -0.225; Results: 'no residual effect of dose on Emax, with the exception of the 0.1 mg/kg dose level, which was included as a covariate'). The 0.10 mg/kg arm is the highest dose of the exposure-response studies (study CL0600-004, Table S3; the 0.15 mg/kg arm of ALX-0600-92001 collected no PD samples).",
      source_name = "dose = 0.10 mg/kg"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 249L,
    n_studies = 10L,
    age_range = "4 months - 79 years (adult and pediatric patients)",
    disease_state = "Adult and pediatric patients (Japanese and non-Japanese) with short bowel syndrome dependent on parenteral support (parenteral nutrition / intravenous fluids).",
    dose_range = "Placebo with standard of care, or subcutaneous teduglutide 0.0125, 0.025, 0.05 or 0.10 mg/kg once daily for 12 weeks to more than 2 years.",
    regions = "Multinational, including Japan.",
    notes = "4918 weekly prescribed PS volume measurements from 249 patients with a baseline value across 10 SBS studies (Results, 'Exposure-response analysis'; Table S3). The long-term extension CL0600-021 (0.05 mg/kg for 24 months) informed ET50, which was estimated at 165 days and fixed to 168 days before the remaining parameters were estimated on all studies (Table S10 note). Pediatric studies recorded daily PN/IV volume in mL/kg/day (Table S3); the analysis variable is L/week. Estimation in Phoenix NLME 8.0. No between-subject variability is reported for the final model (Table S10)."
  )

  ini({
    # ---- Time-course and drug-effect parameters (Marier 2021 Table S10) ----
    # Delta PPSV(t) = Emax * t / (ET50 + t) * DrugEffect (Methods,
    # 'Exposure-response analysis'). Emax is negative (a reduction in
    # prescribed volume), so it is estimated on the linear scale.
    emax <- -5.76; label("Maximum change from baseline in prescribed PS volume at CMAX = 30 ng/mL (L/week)") # Table S10: Emax = -5.76 L/week (RSE 8.71%)
    let50 <- fixed(log(168)); label("Time to 50% of the maximum change from baseline, ET50 (day)") # Table S10: ET50 = 168, Fixed (estimated 165 days on the long-term studies)
    e_cmax_emax <- 0.684; label("Power exponent of CMAX/30 on Emax (unitless)") # Table S10: exponent for the effect of Cmax on Emax = 0.684 (RSE 3.95%)
    e_placebo_emax <- -0.276; label("Proportional shift of Emax in the placebo arm (unitless)") # Table S10: placebo effect on Emax = -0.276 (RSE 8.44%)
    e_dosehigh_emax <- -0.225; label("Proportional shift of Emax at 0.10 mg/kg (unitless)") # Table S10: dose effect (0.10 mg/kg) on Emax = -0.225 (RSE 28.9%)

    # ---- Residual error (Marier 2021 Table S10) ----
    addSd <- 1.11; label("Additive residual error (L/week)") # Table S10: additive error = 1.11 L/week (RSE 0.198%)
  })

  model({
    # ---- 1. Drug-effect multiplier on Emax ----
    # Mutually exclusive strata: placebo; 0.10 mg/kg; otherwise the Cmax power
    # term. The stratified form follows the authors' raw-PPSV model of the same
    # analysis (Table S14), which states the Cmax term applies 'if dose < 0.1
    # mg/kg and not placebo'.
    cmax_term <- (CMAX / 30)^e_cmax_emax
    drug_effect <- PLACEBO * (1 + e_placebo_emax) +
      (1 - PLACEBO) * (DOSE_HIGH * (1 + e_dosehigh_emax) + (1 - DOSE_HIGH) * cmax_term)

    # ---- 2. Time course ----
    et50 <- exp(let50)
    ppsvcfb <- emax * drug_effect * t / (et50 + t)

    # ---- 3. Observation and error ----
    ppsvcfb ~ add(addSd)
  })
}
