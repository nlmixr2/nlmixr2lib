Majid_2016_carbamazepine <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for carbamazepine ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 379). The observed concentration is modelled as the ",
    "average steady-state concentration Css = R0 / (CL/F), R0 = daily dose ",
    "/ 24; there is no volume or absorption term. CL/F rises linearly with ",
    "the carbamazepine daily dose (autoinduction) and is scaled by ",
    "concomitant valproic acid and by the perampanel daily dose. ",
    "Exponential IIV and five-occasion IOV on CL/F; proportional residual ",
    "error. The covariate equation is not printed by the paper; its form ",
    "was deduced from Table 2 and Figure 1 -- see the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S1.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_CBZ_MGD = list(
      description = "Total daily carbamazepine dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters twice: as the input rate R0 = DOSE_CBZ_MGD / 24 (mg/h) in ",
        "Css = R0 / CL, and as an additive linear term on CL/F (Majid 2016 ",
        "Table S1 'Carbamazepine dose effect' 0.00347 L/h per mg/day). ",
        "Dose levels are not reported by the paper."
      ),
      source_name = "carbamazepine dose"
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant valproic acid)",
      notes = "Majid 2016 Table S1 'Valproic acid effect' 0.125, a relative (1 + theta * CONMED_VPA) factor on CL/F.",
      source_name = "valproic acid co-administration"
    ),
    DOSE_PERAMPANEL_MGD = list(
      description = "Total daily perampanel dose (once daily at bedtime); 0 for placebo.",
      units = "mg/day",
      type = "continuous",
      reference_category = "0 mg/day (placebo)",
      notes = paste0(
        "Majid 2016 Table S1 'Perampanel dose effect' 0.00357 per mg, a ",
        "relative (1 + theta * DOSE_PERAMPANEL_MGD) factor on CL/F: +2.9% ",
        "at 8 mg and +4.3% at 12 mg (Table 2: 3.50 -> 3.60 -> 3.65 L/h). ",
        "Perampanel was titrated to 2-12 mg/day."
      ),
      source_name = "perampanel dose"
    ),
    OCC = list(
      description = "Sampling-occasion index for inter-occasion variability: 1 = visit 1 (start of baseline), 2 = visit 2 (start of titration, before the first perampanel dose), 3 = visit 6 (week 10), 4 = visit 7 (week 14), 5 = visit 8 (week 19).",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Majid 2016 Methods, PK model and parameters: 'Inter-occasion variability was defined for visits 1, 2, 6, 7 and 8 as five distinct but identical occasions.' Decomposed inside model() into oc1..oc5.",
      source_name = "visit"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 379L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 33.5 (SD 12.8) years",
    weight_median = "mean 68.2 (SD 17.7) kg",
    sex_female_pct = 51.7,
    race_ethnicity = c(White = 71.0, Asian = 16.1, Chinese = 9.5, Black = 0.5, Other = 2.9),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (269 subjects) or placebo (110); carbamazepine dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (CBZ row). Concentrations from visits 1, 2 (pre-perampanel) and 6, 7, 8 (maintenance weeks 10, 14, 19); two samples 1-2 h apart at each maintenance visit, sampling time relative to the AED dose unknown."
  )

  ini({
    # Majid 2016 Supplemental Table S1 (carbamazepine final model).
    lcl <- log(0.375); label("Dose-independent part of carbamazepine CL/F (L/h)") # Table S1 'Carbamazepine basal clearance (L/h)' 0.375
    e_dose_cbz_cl <- 0.00347; label("Additive increase in CL/F per mg/day of carbamazepine (L/h per mg/day)") # Table S1 'Carbamazepine dose effect' 0.00347
    e_conmed_vpa_cl <- 0.125; label("Relative change in CL/F with concomitant valproic acid (unitless)") # Table S1 'Valproic acid effect' 0.125
    e_dose_perampanel_cl <- 0.00357; label("Relative change in CL/F per mg/day of perampanel (1/(mg/day))") # Table S1 'Perampanel dose effect' 0.00357

    etalcl ~ 0.0512 # Table S1 'IIV' 0.0512 (22.6 %CV)
    # IOV: one variance shared by five identical occasions (Methods).
    etaiov_cl_1 ~ 0.0246 # Table S1 'IOV' 0.0246 (15.7 %CV)
    etaiov_cl_2 ~ fix(0.0246) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0246) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0246) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0246) # same variance as occasion 1

    propSd <- 0.1292; label("Proportional residual error (fraction)") # Table S1 'proportional' 0.0167 is a variance: sqrt(0.0167) = 0.1292, matching the printed 12.9 %CV
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # CL/F = (basal + slope * daily CBZ dose) * (1 + VPA + perampanel dose)
    # with exponential IIV and IOV. The additive-absolute dose term is the
    # only one of the candidate forms consistent with both Table 2 (3.50 L/h
    # at a typical ~900 mg/day) and the flat ~10 ug/mL median of Figure 1A.
    tvcl <- (exp(lcl) + e_dose_cbz_cl * DOSE_CBZ_MGD) *
      (1 + e_conmed_vpa_cl * CONMED_VPA + e_dose_perampanel_cl * DOSE_PERAMPANEL_MGD)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL with R0 the daily dose spread over 24 h; mg/L = ug/mL.
    Cc <- (DOSE_CBZ_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
