Majid_2016_clonazepam <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for clonazepam ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 81). Css = R0 / (CL/F), R0 = daily dose / 24; no volume ",
    "or absorption term. CL/F = basal * (1 + phenytoin + valproic acid + ",
    "clobazam terms); perampanel had no significant effect. Exponential IIV ",
    "and five-occasion IOV on CL/F; proportional residual error."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S3.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_CZP_MGD = list(
      description = "Total daily clonazepam dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_CZP_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "clonazepam daily dose (R0 x 24)"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Majid 2016 Table S3 'Effect of phenytoin' 1.89, additive inside the relative bracket (CL/F x 2.89).",
      source_name = "phenytoin co-administration"
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant valproic acid)",
      notes = "Majid 2016 Table S3 'Effect of valproic acid' 0.785, additive inside the relative bracket. The Results state CL 'was increased when co-administered with phenytoin, valproic acid and clobazam', which rules out reading 0.785 as a multiplicative factor (a decrease).",
      source_name = "valproic acid co-administration"
    ),
    CONMED_CLB = list(
      description = "Concomitant clobazam indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant clobazam)",
      notes = "Majid 2016 Table S3 'Effect of clobazam' 1.01, additive inside the relative bracket.",
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
      notes = "Perampanel presence, dose and concentration were tested on clonazepam CL/F and none was significant (Majid 2016 Results and Table 2, 'Not altered by perampanel'), so no perampanel term is in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 81L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 35.6 (SD 13.3) years",
    weight_median = "mean 71.2 (SD 19.0) kg",
    sex_female_pct = 59.3,
    race_ethnicity = c(White = 84.0, Asian = 4.9, Chinese = 7.4, Black = 0, Other = 3.7),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (44 subjects) or placebo (37); clonazepam dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (CLN row). Supplemental Table S3 is a base = final model."
  )

  ini({
    # Majid 2016 Supplemental Table S3 (clonazepam base/final model). The
    # SEE, 95% CI and %CV columns of S3 are copies of Table S2 (clobazam),
    # including the impossible CI 3.041-6.479 around 2.53; only the estimate
    # column is used. See the vignette Errata.
    lcl <- log(2.53); label("Basal clonazepam CL/F (L/h)") # Table S3 'Clonazepam basal clearance (L/h)' 2.53
    e_conmed_pht_cl <- 1.89; label("Additive relative change in CL/F with concomitant phenytoin (unitless)") # Table S3 'Effect of phenytoin' 1.89
    e_conmed_vpa_cl <- 0.785; label("Additive relative change in CL/F with concomitant valproic acid (unitless)") # Table S3 'Effect of valproic acid' 0.785
    e_conmed_clb_cl <- 1.01; label("Additive relative change in CL/F with concomitant clobazam (unitless)") # Table S3 'Effect of clobazam' 1.01

    etalcl ~ 0.388 # Table S3 'IIV' 0.388 (62.3 %CV; the printed 70.5 is copied from Table S2)
    etaiov_cl_1 ~ 0.0513 # Table S3 'IOV' 0.0513 (22.6 %CV; the printed 16.1 is copied from Table S2)
    etaiov_cl_2 ~ fix(0.0513) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0513) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0513) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0513) # same variance as occasion 1

    propSd <- 0.1308; label("Proportional residual error (fraction)") # Table S3 'proportional' 0.0171 read as a variance: sqrt(0.0171) = 0.1308
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
      e_conmed_vpa_cl * CONMED_VPA + e_conmed_clb_cl * CONMED_CLB)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL in mg/L, x 1000 for ng/mL.
    Cc <- 1000 * (DOSE_CZP_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
