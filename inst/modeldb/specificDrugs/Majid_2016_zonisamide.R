Majid_2016_zonisamide <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for zonisamide ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 93). Css = R0 / (CL/F), R0 = daily dose / 24; no volume ",
    "or absorption term. CL/F = basal * (1 + phenytoin + phenobarbital + ",
    "clobazam terms); perampanel had no significant effect. Exponential IIV ",
    "and five-occasion IOV on CL/F; proportional residual error."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S11.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_ZNS_MGD = list(
      description = "Total daily zonisamide dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_ZNS_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "zonisamide daily dose (R0 x 24)"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Majid 2016 Table S11 'Phenytoin effect' 0.581, additive inside the relative bracket.",
      source_name = "phenytoin co-administration"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes = "Majid 2016 Table S11 'Phenobarbital effect' 0.491, additive inside the relative bracket.",
      source_name = "phenobarbital co-administration"
    ),
    CONMED_CLB = list(
      description = "Concomitant clobazam indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant clobazam)",
      notes = "Majid 2016 Table S11 'Clobazam effect' -0.334, additive inside the relative bracket.",
      source_name = "clobazam co-administration"
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

  covariatesDataExcluded <- list(
    DOSE_PERAMPANEL_MGD = list(
      description = "Total daily perampanel dose.",
      units = "mg/day",
      type = "continuous",
      notes = "Perampanel presence, dose and concentration were tested on zonisamide CL/F and none was significant (Majid 2016 Results and Table 2, 'Not altered by perampanel')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 93L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 33.3 (SD 13.0) years",
    weight_median = "mean 70.4 (SD 19.0) kg",
    sex_female_pct = 57.0,
    race_ethnicity = c(White = 68.8, Asian = 22.6, Chinese = 1.1, Black = 3.2, Other = 4.3),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (60 subjects) or placebo (33); zonisamide dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (ZNS row)."
  )

  ini({
    # Majid 2016 Supplemental Table S11 (zonisamide final model).
    lcl <- log(0.558); label("Basal zonisamide CL/F (L/h)") # Table S11 'Zonisamide basal clearance (L/h)' 0.558
    e_conmed_pht_cl <- 0.581; label("Additive relative change in CL/F with concomitant phenytoin (unitless)") # Table S11 'Phenytoin effect' 0.581
    e_conmed_pb_cl <- 0.491; label("Additive relative change in CL/F with concomitant phenobarbital (unitless)") # Table S11 'Phenobarbital effect' 0.491
    e_conmed_clb_cl <- -0.334; label("Additive relative change in CL/F with concomitant clobazam (unitless)") # Table S11 'Clobazam effect' -0.334

    etalcl ~ 0.235 # Table S11 'IIV' 0.235 (48.5 %CV)
    etaiov_cl_1 ~ 0.0216 # Table S11 'IOV' 0.0216 (14.7 %CV)
    etaiov_cl_2 ~ fix(0.0216) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0216) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0216) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0216) # same variance as occasion 1

    propSd <- 0.1187; label("Proportional residual error (fraction)") # Table S11 'proportional' 0.0141 is a variance: sqrt(0.0141) = 0.1187, matching the printed 11.9 %CV
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    tvcl <- exp(lcl) * (1 + e_conmed_pht_cl * CONMED_PHT +
      e_conmed_pb_cl * CONMED_PB + e_conmed_clb_cl * CONMED_CLB)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_ZNS_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
