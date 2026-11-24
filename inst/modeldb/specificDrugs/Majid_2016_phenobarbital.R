Majid_2016_phenobarbital <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for phenobarbital ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 54). Css = R0 / (CL/F), R0 = daily dose / 24; no volume ",
    "or absorption term. CL/F = basal * (1 + transaminase-elevation + ",
    "lamotrigine + oxcarbazepine terms); perampanel had no significant ",
    "effect. Exponential IIV and five-occasion IOV on CL/F; proportional ",
    "residual error."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S7.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_PB_MGD = list(
      description = "Total daily phenobarbital dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_PB_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "phenobarbital daily dose (R0 x 24)"
    ),
    TRANSAM_GT2ULN = list(
      description = "Transaminase-elevation indicator; 1 = AST or ALT greater than 2 x the upper limit of normal, 0 = otherwise.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (AST and ALT both <= 2 x ULN)",
      notes = "Majid 2016 Table S7 'Effect of liver function enzymes' 0.24, additive inside the relative bracket; footnote a: 'aspartate aminotransferase or alanine aminotransferase greater than 2-times the upper limit of normal'. The Results confirm CL was greater in these subjects.",
      source_name = "liver function enzymes"
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant lamotrigine)",
      notes = "Majid 2016 Table S7 'Lamotrigine effect' -0.24, additive inside the relative bracket.",
      source_name = "lamotrigine co-administration"
    ),
    CONMED_OXC = list(
      description = "Concomitant oxcarbazepine indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant oxcarbazepine)",
      notes = "Majid 2016 Table S7 'Oxcarbazepine effect' -0.325, additive inside the relative bracket.",
      source_name = "oxcarbazepine co-administration"
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
      notes = "Perampanel presence, dose and concentration were tested on phenobarbital CL/F and none was significant (Majid 2016 Results and Table 2, 'Not altered by perampanel')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 54L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 36.6 (SD 14.2) years",
    weight_median = "mean 72.0 (SD 22.8) kg",
    sex_female_pct = 53.7,
    race_ethnicity = c(White = 68.5, Asian = 20.4, Chinese = 7.4, Black = 1.9, Other = 1.9),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (37 subjects) or placebo (17); phenobarbital dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (PHB row). Supplemental Table S7 is a base = final model."
  )

  ini({
    # Majid 2016 Supplemental Table S7 (phenobarbital base/final model).
    lcl <- log(0.243); label("Basal phenobarbital CL/F (L/h)") # Table S7 'Phenobarbital basal clearance (L/h)' 0.243
    e_transam_gt2uln_cl <- 0.24; label("Additive relative change in CL/F with AST or ALT > 2 x ULN (unitless)") # Table S7 'Effect of liver function enzymes' 0.24
    e_conmed_lamotrigine_cl <- -0.24; label("Additive relative change in CL/F with concomitant lamotrigine (unitless)") # Table S7 'Lamotrigine effect' -0.24
    e_conmed_oxc_cl <- -0.325; label("Additive relative change in CL/F with concomitant oxcarbazepine (unitless)") # Table S7 'Oxcarbazepine effect' -0.325

    etalcl ~ 0.281 # Table S7 'IIV' 0.281 (53.0 %CV)
    etaiov_cl_1 ~ 0.0189 # Table S7 'IOV' 0.0189 (13.7 %CV)
    etaiov_cl_2 ~ fix(0.0189) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0189) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0189) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0189) # same variance as occasion 1

    propSd <- 0.1058; label("Proportional residual error (fraction)") # Table S7 'proportional' 0.0112 read as a variance: sqrt(0.0112) = 0.1058 (the printed 1.1 %CV is 100 x 0.0112; see vignette)
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    tvcl <- exp(lcl) * (1 + e_transam_gt2uln_cl * TRANSAM_GT2ULN +
      e_conmed_lamotrigine_cl * CONMED_LAMOTRIGINE + e_conmed_oxc_cl * CONMED_OXC)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_PB_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
