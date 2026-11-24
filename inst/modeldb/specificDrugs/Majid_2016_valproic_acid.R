Majid_2016_valproic_acid <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for valproic acid ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 349). Css = R0 / (CL/F), R0 = daily dose / 24; no ",
    "volume or absorption term. CL/F = basal * (1 + weight + perampanel ",
    "daily-dose terms); exponential IIV and five-occasion IOV on CL/F; ",
    "proportional residual error. The covariate equation is not printed by ",
    "the paper; this form reproduces Table 2 -- see the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S10.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_VPA_MGD = list(
      description = "Total daily valproic acid dose (the drug whose concentration is modelled).",
      units = "mg/d",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_VPA_MGD / 24 (mg/h) in Css = R0 / CL; it does not act on CL/F in this model. Dose levels are not reported by the paper.",
      source_name = "valproic acid daily dose (R0 x 24)"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Majid 2016 Table S10 'Weight effect' 0.490, entered as ",
        "e_wt_cl * (WT / 70) inside the relative bracket (uncentred). Table 2 ",
        "fixes this form: 0.401 x (1 + 0.490) = 0.597 -> 0.60 L/h, whereas a ",
        "power or centred weight term returns 0.401 L/h at the reference. ",
        "Reference weight not printed; 70 kg used."
      ),
      source_name = "weight"
    ),
    DOSE_PERAMPANEL_MGD = list(
      description = "Total daily perampanel dose (once daily at bedtime); 0 for placebo.",
      units = "mg/day",
      type = "continuous",
      reference_category = "0 mg/day (placebo)",
      notes = paste0(
        "Majid 2016 Table S10 'Effect of perampanel dose' 0.00761 per mg, ",
        "additive inside the relative bracket next to the weight term. ",
        "Reproduces Table 2 (0.60 -> 0.62 -> 0.63 L/h at 0 / 8 / 12 mg); the ",
        "paper's '5.0% at 12 mg' is 0.63 / 0.60 of the rounded Table 2 ",
        "values (unrounded 6.1%)."
      ),
      source_name = "perampanel dose"
    ),
    OCC = list(
      description = "Sampling-occasion index for inter-occasion variability: 1 = visit 1 (start of baseline), 2 = visit 2 (start of titration, before the first perampanel dose), 3 = visit 6 (week 10), 4 = visit 7 (week 14), 5 = visit 8 (week 19).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Majid 2016 Methods: five distinct but identical occasions. Decomposed inside model() into oc1..oc5.",
      source_name = "visit"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 349L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 32.5 (SD 12.6) years",
    weight_median = "mean 68.6 (SD 16.9) kg",
    sex_female_pct = 43.8,
    race_ethnicity = c(White = 67.0, Asian = 15.5, Chinese = 13.2, Black = 1.1, Other = 3.2),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (246 subjects) or placebo (103); valproic acid dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (VPA row)."
  )

  ini({
    # Majid 2016 Supplemental Table S10 (valproic acid final model).
    lcl <- log(0.401); label("Basal valproic acid CL/F (L/h)") # Table S10 'Valproic acid basal clearance (L/h)' 0.401
    e_wt_cl <- 0.490; label("Additive relative change in CL/F per unit WT/70 (unitless)") # Table S10 'Weight effect' 0.490
    e_dose_perampanel_cl <- 0.00761; label("Additive relative change in CL/F per mg/day of perampanel (1/(mg/day))") # Table S10 'Effect of perampanel dose' 0.00761

    etalcl ~ 0.316 # Table S10 'IIV' 0.316 (56.2 %CV)
    etaiov_cl_1 ~ 0.0308 # Table S10 'IOV' 0.0308 (17.5 %CV)
    etaiov_cl_2 ~ fix(0.0308) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0308) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0308) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0308) # same variance as occasion 1

    propSd <- 0.1609; label("Proportional residual error (fraction)") # Table S10 'proportional' 0.0259 read as a variance: sqrt(0.0259) = 0.1609 (the printed 2.6 %CV is 100 x 0.0259; see vignette)
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    tvcl <- exp(lcl) * (1 + e_wt_cl * (WT / 70) +
      e_dose_perampanel_cl * DOSE_PERAMPANEL_MGD)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_VPA_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
