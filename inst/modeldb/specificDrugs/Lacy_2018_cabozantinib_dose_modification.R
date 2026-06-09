Lacy_2018_cabozantinib_dose_modification <- function() {
  description <- "Repeated time-to-event (RTTE) hazard model for the 'dose modification of any kind' (DMAK) endpoint in adults with advanced renal cell carcinoma (RCC) treated with oral cabozantinib in the phase III METEOR study (Lacy 2018 exposure-response analysis, n=317 patients with 0-52 events per patient). The instantaneous risk of a dose modification (interruption, reduction, or escalation) depends on whether the patient is currently on an active dose or on a dose interruption. When on active dose (DOSE > 0) the hazard increases log-linearly with the time-varying average cabozantinib plasma concentration CAV. When on a dose interruption (DOSE = 0) the hazard is governed by a separate, larger baseline log-hazard with no cabozantinib effect (cabozantinib effect on hold-state hazard was tested and dropped during base-model development). The drug input Cavg is the individual predicted daily average plasma cabozantinib concentration (ng/mL) derived from the upstream Lacy 2018 popPK model `Lacy_2018_cabozantinib`. Forward simulation exposes `hazard` (instantaneous DMAK rate per day) and `sur` (survival = probability of no DMAK event since t = 0) as derived outputs."
  reference <- paste(
    "Lacy S, Nielsen J, Yang B, Miles D, Nguyen L, Hutmacher M.",
    "Population exposure-response analysis of cabozantinib efficacy",
    "and safety endpoints in patients with renal cell carcinoma.",
    "Cancer Chemother Pharmacol. 2018;81(6):1061-1070.",
    "doi:10.1007/s00280-018-3579-7.",
    "Upstream popPK (driver of CAV): Lacy S et al.,",
    "Cancer Chemother Pharmacol. 2018;81(6):1071-1082;",
    "doi:10.1007/s00280-018-3581-0; see modellib('Lacy_2018_cabozantinib').",
    sep = " "
  )
  vignette <- "Lacy_2018_cabozantinib_exposure_response"

  paper_specific_etas <- c("etalhaz_base")

  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; the time-varying drug input is supplied as the CAV data covariate, in ng/mL, derived from the upstream Lacy 2018 popPK)",
    concentration = "probability (the model output `sur` is the dose-modification-free survival probability, not a drug concentration)"
  )

  covariateData <- list(
    CAV = list(
      description        = "Time-varying individual predicted daily average plasma cabozantinib concentration (ng/mL).",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Per-record time-varying column. Set to 0 during dose interruptions (DOSE = 0). Derived in the source paper from individual empirical-Bayes post-hoc parameters of the upstream Lacy 2018 popPK model (Lacy S et al., Cancer Chemother Pharmacol 2018;81(6):1071-1082; doi:10.1007/s00280-018-3581-0). Source paper's predicted steady-state Cavg values: 375 ng/mL at 20 mg/day, 750 ng/mL at 40 mg/day, 1125 ng/mL at 60 mg/day in RCC patients. The validation vignette shows the recommended cohort construction using the upstream popPK model.",
      source_name        = "Cavg"
    ),
    DOSE = list(
      description        = "Current administered cabozantinib daily dose (mg).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Per-record time-varying column carrying the current daily cabozantinib dose. Set to 0 during dose interruptions; otherwise 20, 40, or 60 mg per the METEOR protocol's stepped-reduction rules. Only used inside model() to derive the binary dose-hold indicator `on_hold = 1 - (DOSE > 0)` that switches between the active-dose and hold-state hazards.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 317L,
    n_studies        = 1L,
    age_range        = "ages reported only at cohort level for METEOR (paper Methods cites enrolment criteria of >=18 years; demographic breakdown for the 317-patient dose-modification subset is in ER Supplemental Table 2)",
    weight_range     = "not separately reported for the dose-modification subset",
    sex_female_pct   = NA_real_,
    race_ethnicity   = "not separately reported for the dose-modification subset",
    disease_state    = "Advanced or metastatic renal cell carcinoma (RCC) with clear-cell histology and measurable disease per RECIST. Subjects had received at least one prior VEGFR-TKI therapy. METEOR randomised 658 patients 1:1 to cabozantinib 60-mg tablet QD vs everolimus 10-mg QD; the dose-modification analysis uses 317 patients on the cabozantinib arm with at least one dose record. Number of DMAK events per patient ranged 0-52 (paper Methods Dose Modification subsection).",
    dose_range       = "Cabozantinib 60-mg tablet (Cabometyx) QD starting dose; protocol allowed reductions to 40 mg and 20 mg and dose interruptions for AE management.",
    regions          = "Multinational phase III trial (NCT01865747); regional breakdown not reported in the Lacy 2018 ER paper.",
    notes            = "Baseline cohort characteristics for the 317-patient subset are in ER Supplemental Table 2 (categorical demographics: baseline ECOG score, MSKCC risk factors, sum of tumor diameters relative to the median, visceral and bone metastases, lung metastases, liver metastases, prior number of VEGF-targeted TKI therapies, organs involved, time to PD on most recent prior TKI). 79% of the 210 MTC patients in the upstream Lacy 2018 popPK cohort dose-reduced from 140 mg; the analogous metric for the RCC METEOR cohort is approximately 60% dose-reduced from 60 mg (paper Discussion)."
  )

  ini({
    # ---- Structural hazard parameters (Lacy 2018 ER Supplemental Table 3, "All Dose Modification Repeated Time to Event (DMAK model)" rows) ----
    # The paper parameterises the active-dose hazard on the log scale as
    #   log h_active(t) = theta_base + theta_drug * Cavg(t) + eta_baseline
    # and the dose-hold hazard as
    #   log h_hold = theta_base_hold
    # (no cabozantinib effect was retained on the hold-state hazard per the
    # paper Methods: "Cabozantinib exposure was tested on the hazard for dose
    # interruption, but this did not result in a reduction in the OFV and was
    # not included in the final model").
    #
    # The published theta_base = -5.4 is the log hazard at Cavg = 0 (units
    # log(1/day)); back-transformed this corresponds to a baseline rate of
    # exp(-5.4) = 0.00452 events/day. Because the paper-reported value is
    # already a log hazard, the canonical `l<base> <- log(<value>)` form
    # collapses to `l<base> <- <paper value>` here.
    lhaz_base       <- -5.4      ; label("Active-dose log baseline hazard for DMAK at CAV = 0 (log(1/day))")          # ER Supplemental Table 3: theta_base = -5.4 (90% CI -5.6, -5.2)
    e_cav_haz       <- 0.000807  ; label("Effect of CAV on active-dose log hazard (per ng/mL)")                       # ER Supplemental Table 3: theta_drug = 0.000807 (90% CI 0.000644, 0.000969)
    lhaz_base_hold  <- -2.7      ; label("Dose-hold log baseline hazard for DMAK (log(1/day))")                       # ER Supplemental Table 3: theta_base-hold = -2.7 (90% CI -2.82, -2.57)

    # ---- Inter-individual variability ----
    # Supplemental Table 3 reports a single "IIV baseline (omega)" value
    # 0.655 in variance units (omega^2), matching the convention used by the
    # companion popPK Lacy_2018_cabozantinib (footnote d, omega^2 reported
    # under the "omega" symbol). The paper does not enumerate which state's
    # baseline the IIV attaches to; this extraction applies the eta only to
    # the active-dose baseline lhaz_base (the source row is labelled
    # "baseline", which is the symbol used for the active-state log hazard
    # theta_base; the hold-state log hazard is the separately labelled
    # theta_base-hold). The vignette Assumptions and deviations section
    # records this interpretation and discusses the alternative (shared
    # subject-level eta on both states) that is mathematically defensible
    # but not explicitly stated in the paper.
    #
    # `etalhaz_base` is declared paper-specific (see `paper_specific_etas`
    # above) because the IIV attaches to a paper-mechanistic log-hazard
    # parameter rather than to a standard PK structural parameter.
    etalhaz_base ~ 0.655                                  # ER Supplemental Table 3: IIV baseline (omega) = 0.655 (90% CI 0.507, 0.803)

    # No residual error. As a parametric hazard model the likelihood is the
    # survival / event density itself ($EST ... LIKE in NONMEM), not an
    # observation-error model. Forward simulation in rxode2 exposes
    # `hazard` and `sur` as derived outputs.
  })

  model({
    # ---- Subject-specific active-dose log baseline hazard ----
    lhaz_active <- lhaz_base + etalhaz_base

    # ---- State indicator ----
    # `on_hold` is a binary indicator that is 1 during dose interruptions
    # (DOSE = 0) and 0 during active dosing (DOSE > 0). The (DOSE > 0)
    # comparison in rxode2 yields 0 / 1 in the underlying C, so the
    # arithmetic mixture below selects the active-dose hazard when DOSE > 0
    # and the hold-state hazard when DOSE <= 0.
    on_hold <- 1 - (DOSE > 0)

    # ---- Instantaneous hazard ----
    # log h_active(t) = lhaz_base + eta_baseline + theta_drug * Cavg(t)
    # log h_hold      = lhaz_base_hold
    # Cavg (ng/mL) enters linearly on the log-hazard scale in the active
    # state; the dose-hold hazard is concentration-independent.
    hazard_active <- exp(lhaz_active + e_cav_haz * CAV)
    hazard_hold   <- exp(lhaz_base_hold)
    hazard        <- hazard_active * (1 - on_hold) + hazard_hold * on_hold

    # ---- Cumulative hazard and survival ----
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)
  })
}
