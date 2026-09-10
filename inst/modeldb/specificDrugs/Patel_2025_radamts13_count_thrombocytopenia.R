Patel_2025_radamts13_count_thrombocytopenia <- function() {
  description <- "Exposure-response count model for thrombocytopenia events in patients with congenital thrombotic thrombocytopenic purpura (cTTP) treated with recombinant ADAMTS13 (rADAMTS13; TAK-755) or plasma-based therapy. Poisson model with a between-subject random effect on the baseline count parameter and a sigmoid Emax inhibitory drug effect driven by the average plasma ADAMTS13 activity over a dosing interval (CAV, IU/mL). Fitted to event counts within prophylaxis Periods 1 and 2 of the pivotal phase III crossover study (NCT03393975) in 41 patients of all ages. Thrombocytopenia is defined as a platelet count decreased by at least 25% from baseline or a platelet count below 150,000/uL. The drug effect is steep and close to saturation at rADAMTS13 exposures: Emax is a 91.8% reduction in the expected count and ECave50 is 0.0149 IU/mL (about 1.5% of normal ADAMTS13 activity). Exposure comes from the companion population PK model modellib('Patel_2025_radamts13')."
  reference <- paste(
    "Patel M, Xu H, Barriere O, Diderichsen P, Patwari P, Zhu AZX,",
    "Marier JF, Peyret T, Wang LT, Mellgard B, Wang W, Bhattacharya I.",
    "Use of PopPK and E-R Analyses toward Explaining Causal Link Between",
    "ADAMTS13 in Recombinant vs. Plasma-Based Therapies and Clinical",
    "Effects in cTTP. Clin Pharmacol Ther. 2025;118(4):813-822.",
    "doi:10.1002/cpt.3720.",
    "Parameter values from Table S5 of the Supplementary Information.",
    "Exposure driver from the companion PopPK model in Table 1 of the same",
    "paper; see modellib('Patel_2025_radamts13').",
    sep = " "
  )
  vignette <- "Patel_2025_radamts13_exposure_response"

  # `etab1` carries the between-subject variability on a paper-mechanistic
  # count-model baseline rather than on a standard PK structural parameter.
  paper_specific_etas <- c("etab1")

  units <- list(
    time          = "n/a (the count is accumulated over a whole prophylaxis period, roughly 6 months; the model has no time argument and no ODE states)",
    dosing        = "n/a (no drug-dosing events; the drug input is the CAV data covariate, in IU/mL)",
    concentration = "IU/mL (the CAV covariate; 1 IU/mL = 100% of normal plasma ADAMTS13 activity)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Mean average plasma ADAMTS13 activity over a dosing interval for the patient in the prophylaxis period being scored (IU/mL).",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. One value per patient per prophylaxis period, not a time-varying column: the paper derives a Cave after each dose in each patient from the companion PopPK model's individual post hoc parameters, then averages those within Period 1 and within Period 2 to give the single mean Cave that scores that period's count (paper Methods, 'Exposure-response count modeling'). 1 IU/mL = 100% of normal ADAMTS13 activity, so the fitted ECave50 of 0.0149 IU/mL corresponds to about 1.5% of normal. The companion PopPK model modellib('Patel_2025_radamts13') returns Cc in IU/L, so a solve of that model must be divided by 1000 before it is used here; the validation vignette does this explicitly. Reference exposures from Table S2: 0.0308 / 0.0613 IU/mL for PBT 10 IU/kg Q2W / Q1W and 0.203 / 0.405 IU/mL for rADAMTS13 40 IU/kg Q2W / Q1W.",
      source_name        = "Cave"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "All ages. The count analysis was run first in patients aged >=12 years and then repeated without age restriction; the parameters extracted here are the all-ages final model (paper Methods and Results).",
    weight_range   = "Not reported separately for the 41-patient exposure-response set; the parent PK analysis set spans 18.3-130.0 kg (Table 2).",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported separately for the 41-patient exposure-response set.",
    disease_state  = "Congenital thrombotic thrombocytopenic purpura (cTTP). Endpoint: thrombocytopenia, defined as a platelet count decreased by at least 25% from baseline or a platelet count <150,000/uL (paper Methods, 'Exposure-response analyses').",
    dose_range     = "rADAMTS13 40 IU/kg IV Q1W or Q2W; PBT approximately 10 IU/kg IV.",
    regions        = "Multinational phase III crossover study NCT03393975.",
    notes          = "The exposure-response analyses use only the pivotal phase III study, Periods 1 and 2, giving N = 41 (Table S5 and Table S7 headers). Two patients included in the exploratory data analysis and the repeated time-to-event modeling were not counted as part of the prophylaxis cohort in the count exposure-response analysis (Figure 1 footnote a), which is why the count-model N (41) is smaller than the RTTE cohort and much smaller than the 65-patient PK analysis set."
  )

  ini({
    # ---- Baseline count ----
    # Table S5 reports B1 on the log scale (estimate 1.40) alongside its
    # back-transform in the "Untransformed values of estimate" column (4.05 =
    # exp(1.40)), which fixes the scale unambiguously. Because the reported
    # value is already a log, the canonical `l<name> <- log(<value>)` form
    # collapses to a bare assignment here.
    b1 <- 1.40; label("Log baseline expected thrombocytopenia count per period at zero ADAMTS13 activity (log counts)")  # Table S5, Thrombocytopenia counts: B1 = 1.40 (RSE 35.1%), untransformed 4.05

    # ---- Sigmoid Emax drug effect ----
    emax_count <- 0.918;  label("Maximum fractional reduction in the expected thrombocytopenia count (fraction)")  # Table S5: Emax = 0.918 (RSE 4.25%). Quoted in Results as "Emax: 91.8% reduction in thrombocytopenia event counts"
    ec50       <- 0.0149; label("Average ADAMTS13 activity giving half the maximum count reduction (IU/mL)")       # Table S5: EC50 = 0.0149 (RSE 68.1%). Quoted in Results as "ECave50: 0.0149 IU/mL"
    gamma      <- 2.58;   label("Hill coefficient of the ADAMTS13 activity-response relationship (unitless)")      # Table S5: Gamma = 2.58 (RSE 55.2%). Quoted in Results as "gamma: 2.58"

    # ---- Between-subject variability ----
    # Table S5 gives "BSV on B1" = 3.67 with an RSE of 35.0% and, unlike the
    # repeated time-to-event Table S6, carries no footnote declaring a
    # percent-CV back-transform. The value is therefore read as the NONMEM
    # OMEGA variance on the log-count scale (omega^2 = 3.67, omega = 1.92).
    # The alternative readings are excluded on magnitude: 3.67 as a percent CV
    # would be a 3.67% between-subject spread, which cannot produce the "count
    # zero as well as the long-tailed distribution of PBT" that the paper
    # credits this random effect with describing (Results, 'Exposure-response
    # count modeling'). See the vignette's Assumptions and deviations section.
    etab1 ~ 3.67  # Table S5, Thrombocytopenia counts: BSV on B1 = 3.67 (RSE 35.0%)
  })

  model({
    # 1. Subject-specific baseline expected count at zero ADAMTS13 activity.
    lam0 <- exp(b1 + etab1)

    # 2. Sigmoid Emax inhibition by average ADAMTS13 activity. CAV = 0 gives
    #    drug = 0 exactly (0^gamma = 0), so a placebo or untreated record
    #    collapses to the baseline count.
    drug <- emax_count * CAV^gamma / (CAV^gamma + ec50^gamma)

    # 3. Expected count over the prophylaxis period.
    lambda <- lam0 * (1 - drug)

    # 4. Poisson count observation.
    thrombocytopenia_count ~ pois(lambda)
  })
}
