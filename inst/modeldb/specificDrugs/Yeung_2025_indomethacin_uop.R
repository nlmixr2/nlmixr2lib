Yeung_2025_indomethacin_uop <- function() {
  description <- paste(
    "Linear mixed-effects exposure-toxicity model for the urine output",
    "(UOP) rate nadir within 12 h after an indomethacin dose in preterm",
    "neonates at Children's Mercy Hospital: ln(UOP nadir) falls linearly",
    "with the model-predicted per-dose indomethacin Cmax (8.4% lower UOP",
    "nadir per 1 ug/mL) and rises with weight-normalised 24-h fluid",
    "intake relative to the population median, with a subject random",
    "intercept. The median UOP nadir reaches 1 and 0.5 mL/kg/h at Cmax 14",
    "and 22 ug/mL. Static per-dose regression: no ODE, no dosing events.",
    sep = " "
  )
  reference <- paste(
    "Yeung CHT, Sekulich DC, Scott A, Nolte WM, Gibson K, Su R,",
    "Alrifai MW, Lopata SM, Lewis T, Reese J, Maharaj A (2025).",
    "The Relationship of Indomethacin Exposure With Efficacy and Renal",
    "Toxicity Outcomes for Preterm Infants in the Neonatal Intensive Care",
    "Unit. Clin Transl Sci 18(5):e70251. doi:10.1111/cts.70251.",
    sep = " "
  )
  vignette <- "Yeung_2025_indomethacin"
  units <- list(
    time = "n/a (static per-dose regression; the outcome is the lowest UOP rate within 12 h after each dose)",
    dosing = "n/a (no dose events; exposure enters through the CMAX covariate)",
    concentration = "urine_output_nadir (UOP rate nadir within 12 h after a dose, mL/kg/h)"
  )

  covariateData <- list(
    CMAX = list(
      description = "Individual model-predicted peak plasma indomethacin concentration after the dose (Cmax per dose).",
      units = "ug/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "(a) Total plasma concentration. (b) Peak after each individual",
        "dose, not a steady-state value; accumulation over a course is",
        "included. (c) Simulated from empirical Bayes estimates of the IND",
        "popPK model of Krzyzanski 2023 (Yeung 2025 Methods 'Evaluation of",
        "Renal Toxicity'); Yeung_2025_indomethacin_cmh is the closest",
        "library model. Enters linearly and uncentred. Observed values in",
        "the dataset were below about 5 ug/mL (Figure 4)."
      ),
      source_name = "Cmax"
    ),
    FLUID_IN_24H = list(
      description = "Total fluid intake over the 24 h after the dose.",
      units = "mL/24h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The paper's covariate is the weight-normalised fluid intake in the",
        "24 h after the dose minus the population median ('delta fluid",
        "input'). model() forms it as FLUID_IN_24H / WT - 150, where 150",
        "mL/kg/day is the CMH median of Table 1. Set FLUID_IN_24H = 150 *",
        "WT to reproduce the paper's median-fluid simulation."
      ),
      source_name = "delta fluid intake (with WT)"
    ),
    WT = list(
      description = "Body weight at the time of the dose.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Used only to weight-normalise FLUID_IN_24H.",
      source_name = "Neonate weight"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 53,
    n_studies = 1,
    age_range = "postnatal age at first dose 0.23 h to 27.4 days",
    weight_range = "360-2630 g at first dose",
    ga_range = "22.6-34 weeks gestational age at birth",
    disease_state = "Preterm neonates receiving IV indomethacin for IVH prophylaxis or PDA treatment",
    dose_range = "0.1 mg/kg IV q24h x3 (IVH) or q12h x3 per course (PDA); CMH median 0.107 (0.085-0.40) mg/kg",
    regions = "USA (Children's Mercy Hospital, Kansas City)",
    notes = paste(
      "Renal toxicity was quantified in the CMH dataset only (Yeung 2025",
      "Methods). Table 1 lists 45 IVH and 18 PDA infants, 10 of whom are",
      "in both, giving 53 unique infants; 1820 UOP rates in total. Median",
      "UOP nadir 3.3 (0.95-9.5) mL/kg/h and median fluid intake 150",
      "(64.5-227.7) mL/kg/day in the 24 h after a dose."
    )
  )

  ini({
    lurine_output_nadir <- 1.23; label("ln UOP rate nadir at Cmax = 0 and median fluid intake (ln mL/kg/h)") # Table 2, UOP block 'Intercept' = 1.23 (95% CI 1.09, 1.37)
    e_cmax_urine_output_nadir <- -0.088; label("Change in ln UOP nadir per ug/mL of Cmax (mL/ug)") # Table 2, UOP block 'Cmax (ug/mL)' = -0.088 (95% CI -0.17, -0.00027), p = 0.044
    e_fluid_urine_output_nadir <- 0.0034; label("Change in ln UOP nadir per mL/kg/day of fluid intake above the median") # Table 2, UOP block 'delta fluid intake' printed 0.00034 (mL/kg/h) with 95% CI 0.00068, 0.0061 and p = 0.014; value and unit taken from the CI (see vignette)

    etalurine_output_nadir ~ fixed(0) # subject random intercept (Methods); variance not reported, see expSd_urine_output_nadir

    expSd_urine_output_nadir <- 0.53; label("Total SD of ln UOP nadir, between-subject plus residual (unitless)") # back-solved from Results: 5th percentile reaches 1 and 0.5 mL/kg/h at Cmax 5 and 12 ug/mL
  })

  model({
    dfluid <- FLUID_IN_24H / WT - 150
    urine_output_nadir <- exp(lurine_output_nadir + etalurine_output_nadir +
      e_cmax_urine_output_nadir * CMAX +
      e_fluid_urine_output_nadir * dfluid)
    urine_output_nadir ~ lnorm(expSd_urine_output_nadir)
  })
}
