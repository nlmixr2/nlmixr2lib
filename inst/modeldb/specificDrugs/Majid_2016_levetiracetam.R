Majid_2016_levetiracetam <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for levetiracetam ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 330). Css = R0 / (CL/F), R0 = daily dose / 24; no ",
    "volume or absorption term. CL/F = basal * (1 + sex + weight + ",
    "phenytoin + valproic acid terms); perampanel had no significant ",
    "effect. Exponential IIV and five-occasion IOV on CL/F; proportional ",
    "residual error. The covariate equation is not printed by the paper; ",
    "the weight term follows the form Table 2 fixes for clobazam and ",
    "valproic acid -- see the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S5.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_LEV_MGD = list(
      description = "Total daily levetiracetam dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_LEV_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "levetiracetam daily dose (R0 x 24)"
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Majid 2016 Table S5 'Sex effect' -0.244, additive inside the relative bracket; the Results say CL 'was lower in female subjects', so 1 = female.",
      source_name = "sex"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Majid 2016 Table S5 'Weight effect' 0.389, entered as e_wt_cl * (WT / 70) inside the relative bracket (uncentred), the form Table 2 fixes for the clobazam and valproic-acid weight terms of the same analysis. Reference weight not printed; 70 kg used.",
      source_name = "weight"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Majid 2016 Table S5 'Phenytoin effect' -0.252, additive inside the relative bracket.",
      source_name = "phenytoin co-administration"
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant valproic acid)",
      notes = "Majid 2016 Table S5 'Valproic acid effect' -0.154, additive inside the relative bracket.",
      source_name = "valproic acid co-administration"
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
      notes = "Perampanel presence, dose and concentration were tested on levetiracetam CL/F and none was significant (Majid 2016 Results and Table 2, 'Not altered by perampanel')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 330L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 35.3 (SD 13.6) years",
    weight_median = "mean 72.3 (SD 17.9) kg",
    sex_female_pct = 52.1,
    race_ethnicity = c(White = 80.3, Asian = 8.8, Chinese = 5.8, Black = 3.6, Other = 1.5),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (235 subjects) or placebo (95); levetiracetam dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (LEV row). Supplemental Table S5 is a base = final model."
  )

  ini({
    # Majid 2016 Supplemental Table S5 (levetiracetam base/final model).
    lcl <- log(2.28); label("Basal levetiracetam CL/F (L/h)") # Table S5 'Levetiracetam basal clearance (L/h)' 2.28
    e_sexf_cl <- -0.244; label("Additive relative change in CL/F for females (unitless)") # Table S5 'Sex effect' -0.244
    e_wt_cl <- 0.389; label("Additive relative change in CL/F per unit WT/70 (unitless)") # Table S5 'Weight effect' 0.389
    e_conmed_pht_cl <- -0.252; label("Additive relative change in CL/F with concomitant phenytoin (unitless)") # Table S5 'Phenytoin effect' -0.252
    e_conmed_vpa_cl <- -0.154; label("Additive relative change in CL/F with concomitant valproic acid (unitless)") # Table S5 'Valproic acid effect' -0.154

    etalcl ~ 0.23 # Table S5 'IIV' 0.23 (48.0 %CV)
    etaiov_cl_1 ~ 0.034 # Table S5 'IOV' 0.034 (18.4 %CV)
    etaiov_cl_2 ~ fix(0.034) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.034) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.034) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.034) # same variance as occasion 1

    propSd <- 0.2678; label("Proportional residual error (fraction)") # Table S5 'proportional' 0.0717 read as a variance: sqrt(0.0717) = 0.2678
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    tvcl <- exp(lcl) * (1 + e_sexf_cl * SEXF + e_wt_cl * (WT / 70) +
      e_conmed_pht_cl * CONMED_PHT + e_conmed_vpa_cl * CONMED_VPA)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_LEV_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
