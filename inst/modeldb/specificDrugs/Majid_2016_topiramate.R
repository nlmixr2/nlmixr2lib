Majid_2016_topiramate <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for topiramate ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 226). Css = R0 / (CL/F), R0 = daily dose / 24; no ",
    "volume or absorption term. CL/F = basal * (1 + weight + phenytoin + ",
    "zonisamide terms); perampanel had no significant effect. Exponential ",
    "IIV and five-occasion IOV on CL/F; proportional residual error. The ",
    "covariate equation is not printed by the paper; the weight term ",
    "follows the form Table 2 fixes for clobazam and valproic acid -- see ",
    "the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S9.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_TPM_MGD = list(
      description = "Total daily topiramate dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_TPM_MGD / 24 (mg/h) in Css = R0 / CL; it does not act on CL/F in this model. Dose levels are not reported by the paper.",
      source_name = "topiramate daily dose (R0 x 24)"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Majid 2016 Table S9 'Effect of weight' 0.878, entered as e_wt_cl * (WT / 70) inside the relative bracket (uncentred), the form Table 2 fixes for the clobazam and valproic-acid weight terms of the same analysis. Reference weight not printed; 70 kg used.",
      source_name = "weight"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Majid 2016 Table S9 'Phenytoin effect' 0.714, additive inside the relative bracket.",
      source_name = "phenytoin co-administration"
    ),
    CONMED_ZNS = list(
      description = "Concomitant zonisamide indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant zonisamide)",
      notes = "Majid 2016 Table S9 'Zonisamide effect' 1.12, additive inside the relative bracket.",
      source_name = "zonisamide co-administration"
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
      notes = "Perampanel presence, dose and concentration were tested on topiramate CL/F and none was significant (Majid 2016 Results and Table 2, 'Not altered by perampanel')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 226L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 33.6 (SD 13.7) years",
    weight_median = "mean 68.0 (SD 18.0) kg",
    sex_female_pct = 50.4,
    race_ethnicity = c(White = 77.0, Asian = 10.2, Chinese = 8.4, Black = 1.8, Other = 2.7),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (162 subjects) or placebo (64); topiramate dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (TPM row). Supplemental Table S9 is a base = final model."
  )

  ini({
    # Majid 2016 Supplemental Table S9 (topiramate base/final model).
    lcl <- log(0.953); label("Basal topiramate CL/F (L/h)") # Table S9 'Topiramate basal clearance (L/h)' 0.953
    e_wt_cl <- 0.878; label("Additive relative change in CL/F per unit WT/70 (unitless)") # Table S9 'Effect of weight' 0.878
    e_conmed_pht_cl <- 0.714; label("Additive relative change in CL/F with concomitant phenytoin (unitless)") # Table S9 'Phenytoin effect' 0.714
    e_conmed_zns_cl <- 1.12; label("Additive relative change in CL/F with concomitant zonisamide (unitless)") # Table S9 'Zonisamide effect' 1.12

    etalcl ~ 0.212 # Table S9 'IIV' 0.212 (46.0 %CV)
    etaiov_cl_1 ~ 0.0292 # Table S9 'IOV' 0.0292 (17.1 %CV)
    etaiov_cl_2 ~ fix(0.0292) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0292) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0292) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0292) # same variance as occasion 1

    propSd <- 0.1667; label("Proportional residual error (fraction)") # Table S9 'proportional' 0.0278 read as a variance: sqrt(0.0278) = 0.1667 (the printed 2.8 %CV is 100 x 0.0278; see vignette)
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
      e_conmed_pht_cl * CONMED_PHT + e_conmed_zns_cl * CONMED_ZNS)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_TPM_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
