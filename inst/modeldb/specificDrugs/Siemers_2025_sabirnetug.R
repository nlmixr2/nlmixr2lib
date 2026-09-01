Siemers_2025_sabirnetug <- function() {
  description <- paste(
    "Direct Emax exposure-response model of central target engagement for sabirnetug",
    "(ACU193), a humanized IgG2 monoclonal antibody selective for soluble globular",
    "amyloid-beta oligomers (AbetaOs), in participants with mild cognitive impairment",
    "or mild dementia due to Alzheimer's disease (INTERCEPT-AD phase 1, NCT04931459).",
    "Target engagement is the cerebrospinal-fluid sabirnetug-AbetaO complex",
    "concentration measured by an anti-idiotype capture / AbetaO detection immunoassay",
    "and reported in arbitrary units per mL; it is driven by the CSF sabirnetug",
    "concentration supplied as the covariate CEFFECT.",
    "PD-only model: Siemers 2025 characterised sabirnetug serum PK non-compartmentally",
    "(Table 3) and never fitted a compartmental population PK model, so no sabirnetug",
    "PK model is packaged here and users must supply their own CSF concentration",
    "trajectory. The companion population-modelling analysis is deferred to a separate",
    "publication (Trame 2023 conference abstract, reference [30] of the paper).",
    sep = " "
  )

  reference <- paste(
    "Siemers E, Feaster T, Sethuraman G, Sundell K, Skljarevski V, Cline EN, Zhang H,",
    "Jerecic J, Honig LS, Salloway S, Sperling R, Trame MN, Dodds MG, Johnson K.",
    "INTERCEPT-AD, a phase 1 study of intravenous sabirnetug in participants with mild",
    "cognitive impairment or mild dementia due to Alzheimer's disease.",
    "J Prev Alzheimers Dis. 2025;12:100005. doi:10.1016/j.tjpad.2024.100005.",
    sep = " "
  )

  vignette <- "Siemers_2025_sabirnetug"

  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external CSF sabirnetug-concentration covariate)",
    concentration = paste(
      "(observation targetEngagement is the CSF sabirnetug-AbetaO complex concentration",
      "in AU/mL; driving covariate CEFFECT is the CSF sabirnetug concentration in ng/mL)"
    )
  )

  covariateData <- list(
    CEFFECT = list(
      description = paste(
        "Instantaneous cerebrospinal-fluid sabirnetug concentration at the time of each",
        "target-engagement observation, supplied as a time-varying covariate from",
        "observed lumbar-puncture samples or an upstream PK source."
      ),
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row. Drives the direct Emax expression",
        "targetEngagement = emax * CEFFECT / (CEFFECT + ec50).",
        "CSF is the biophase for this model: the AbetaO target and the measured",
        "sabirnetug-AbetaO complex both reside in CSF, and Siemers 2025 fits the",
        "exposure-response against the CSF drug concentration rather than the serum",
        "concentration (Siemers 2025 Section 2.9 and Fig. 3). The two differ by more",
        "than an order of magnitude -- the paper reports mean CSF-to-serum percent",
        "ratios of 1.65% to 3.25% across cohorts (Section 3.7) -- so a serum",
        "concentration must NOT be substituted here. Note also that the CSF measurement",
        "is total (bound plus unbound) drug whereas the serum measurement is free drug",
        "(Section 3.7). Observed per-cohort median (range) CSF sabirnetug concentrations",
        "were 15.7 (6.7-29.2), 98.0 (36.8-130.7), 169.5 (65.0-334.7), 282.7 (26.1-455.9),",
        "148.9 (5.2-255.9), 1161.6 (48.1-1722.2) and 869.8 (419.0-1474.1) ng/mL for",
        "cohorts 1-7 respectively (Section 3.6). Set to 0 for the drug-free reference,",
        "which collapses target engagement to 0.",
        "This model reuses the general-scope CEFFECT canonical for the effect-site /",
        "biophase PD-driver concept; the sabirnetug identity and the ng/mL scale are",
        "conveyed by this per-model entry, following the same pattern as the vitreous",
        "-humour pegcetacoplan driver in Crass_2025_pegcetacoplan_ga_exposureresponse.R."
      ),
      source_name = "CSF [Sabirnetug] (Siemers 2025 Fig. 3 x-axis; C in the Section 2.9 Emax equation)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 1L,
    age_range = "55-90 years (protocol eligibility)",
    age_median = "72.3 years (mean, SD 7.9, all sabirnetug-treated participants)",
    weight_range = "41-113 kg (protocol eligibility; individual weights not reported)",
    weight_median = NA_character_,
    sex_female_pct = 55.1,
    race_ethnicity = c(
      White = 93.9,
      Black = 4.1,
      `Native American/Alaskan` = 2.0,
      Asian = 0
    ),
    disease_state = paste(
      "Early symptomatic Alzheimer's disease -- mild cognitive impairment or mild",
      "dementia due to AD by National Institute on Aging - Alzheimer's Association",
      "criteria, with a Global Clinical Dementia Rating of 0.5 or 1.0, Mini-Mental State",
      "Examination 18-30, and a positive amyloid PET scan (composite SUVr > 1.2,",
      "> 26.2 Centiloids). Mean (SD) baseline MMSE 24.1 (3.7), CDR-GS 0.6 (0.3),",
      "CDR-SB 3.6 (1.9), BMI 28.0 (5.4) kg/m2. APOE e4 carriers: 36.6% heterozygous,",
      "12.5% homozygous."
    ),
    dose_range = paste(
      "Part A (single ascending dose): one 2, 10, 25 or 60 mg/kg IV infusion.",
      "Part B (multiple ascending dose): three infusions of 10 or 60 mg/kg every four",
      "weeks, or 25 mg/kg every two weeks."
    ),
    regions = "United States (15 study centers)",
    notes = paste(
      "INTERCEPT-AD (NCT04931459), a randomized, double-blind, placebo-controlled",
      "first-in-human phase 1 study conducted 23 June 2021 to 12 June 2023. Sixty-five",
      "participants were enrolled (Part A N = 32, Part B N = 33); 48 received sabirnetug",
      "and constitute the pharmacokinetic population against which Fig. 3 is plotted.",
      "Baseline demographics are Siemers 2025 Table 1. CSF was sampled by lumbar",
      "puncture at baseline and on day 21 (Part A) and on days 70, 63 and 35 for Part B",
      "cohorts 5, 6 and 7 respectively (Section 2.8), so each participant contributes",
      "essentially one post-dose target-engagement observation and the exposure-response",
      "is fitted across participants rather than over a within-participant time course.",
      "No between-participant variability, no residual-error magnitude and no covariate",
      "effects are reported for this Emax analysis, so the model is typical-value only.",
      "Siemers 2025 Section 3.8 states that target engagement showed no correlation with",
      "APOE e4 genotype, presence of ARIA, or baseline amyloid burden."
    )
  )

  ini({
    # Direct Emax exposure-response of CSF target engagement on CSF sabirnetug
    # concentration, E = Emax * C / (C + EC50) (Siemers 2025 Section 2.9). No Hill
    # coefficient appears in the published equation, so the model is the plain
    # hyperbolic Emax form with an implicit Hill exponent of 1.
    #
    # Both point estimates are PRINTED, not digitised: they are annotated inside the
    # Fig. 3 plot panel ("Emax = 22.71 AU/mL Complex" / "EC50 = 136 ng/mL ACU193").
    # The Emax value is independently corroborated in the Section 3.8 body text
    # ("Emax = 22.71 AU/mL sabirnetug-AbetaO complex"); EC50 appears only in the
    # figure panel. Both were estimated by the authors, so neither is fixed().
    lemax <- log(22.71); label("Maximum CSF target engagement, sabirnetug-AbetaO complex (AU/mL)")  # Siemers 2025 Fig. 3 panel annotation; corroborated in Section 3.8 text

    lec50 <- log(136); label("CSF sabirnetug concentration eliciting half-maximal target engagement (ng/mL)")  # Siemers 2025 Fig. 3 panel annotation
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)

    # CSF sabirnetug-AbetaO complex (AU/mL). CEFFECT is the per-record CSF
    # sabirnetug concentration in ng/mL.
    targetEngagement <- emax * CEFFECT / (CEFFECT + ec50)
  })
}
