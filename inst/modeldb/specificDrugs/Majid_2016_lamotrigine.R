Majid_2016_lamotrigine <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for lamotrigine ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 356). Css = R0 / (CL/F), R0 = daily dose / 24; no ",
    "volume or absorption term. CL/F = basal * (1 + carbamazepine + ",
    "phenobarbital + valproic acid + slope * ln(perampanel daily dose)); ",
    "exponential IIV and five-occasion IOV on CL/F; proportional residual ",
    "error. The covariate equation is not printed by the paper; this form ",
    "reproduces Table 2 -- see the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S4.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_LTG_MGD = list(
      description = "Total daily lamotrigine dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_LTG_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "lamotrigine daily dose (R0 x 24)"
    ),
    CONMED_CBZ = list(
      description = "Concomitant carbamazepine indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes = "Majid 2016 Table S4 'Effect of carbamazepine' 0.584, additive inside the relative bracket.",
      source_name = "carbamazepine co-administration"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes = "Majid 2016 Table S4 'Effect of phenobarbital' 0.742, additive inside the relative bracket.",
      source_name = "phenobarbital co-administration"
    ),
    CONMED_VPA = list(
      description = "Concomitant valproic acid indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant valproic acid)",
      notes = "Majid 2016 Table S4 'Effect of valproic acid' -0.51, additive inside the relative bracket (CL/F x 0.49).",
      source_name = "valproic acid co-administration"
    ),
    DOSE_PERAMPANEL_MGD = list(
      description = "Total daily perampanel dose (once daily at bedtime); 0 for placebo.",
      units = "mg/day",
      type = "continuous",
      reference_category = "0 mg/day (placebo)",
      notes = paste0(
        "Majid 2016 Table S4 'Effect of perampanel dose: slope log(dose)' ",
        "0.037; the Results say CL increased 'proportionally to the loge of ",
        "perampanel daily dose'. Enters as e * ln(DOSE_PERAMPANEL_MGD) inside ",
        "the relative bracket, and is 0 for placebo: model() evaluates ",
        "ln(max(DOSE_PERAMPANEL_MGD, 1)), which equals ln(dose) for every ",
        "studied dose (2-12 mg) and 0 at 0 mg. Reproduces Table 2 (1.18 -> ",
        "1.27 -> 1.29 L/h at 0 / 8 / 12 mg)."
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
    n_subjects = 356L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 34.7 (SD 12.5) years",
    weight_median = "mean 73.1 (SD 17.3) kg",
    sex_female_pct = 53.4,
    race_ethnicity = c(White = 79.5, Asian = 9.0, Chinese = 7.0, Black = 2.0, Other = 2.5),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (257 subjects) or placebo (99); lamotrigine dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (LTG row)."
  )

  ini({
    # Majid 2016 Supplemental Table S4 (lamotrigine final model).
    lcl <- log(1.18); label("Basal lamotrigine CL/F (L/h)") # Table S4 'Lamotrigine basal clearance (L/h)' 1.18
    e_conmed_cbz_cl <- 0.584; label("Additive relative change in CL/F with concomitant carbamazepine (unitless)") # Table S4 'Effect of carbamazepine' 0.584
    e_conmed_pb_cl <- 0.742; label("Additive relative change in CL/F with concomitant phenobarbital (unitless)") # Table S4 'Effect of phenobarbital' 0.742
    e_conmed_vpa_cl <- -0.51; label("Additive relative change in CL/F with concomitant valproic acid (unitless)") # Table S4 'Effect of valproic acid' -0.51
    e_dose_perampanel_cl <- 0.037; label("Additive relative change in CL/F per unit ln(perampanel mg/day) (unitless)") # Table S4 'Effect of perampanel dose: slope log(dose)' 0.037

    etalcl ~ 0.249 # Table S4 'IIV' 0.249 (49.9 %CV)
    etaiov_cl_1 ~ 0.0296 # Table S4 'IOV' 0.0296 (17.2 %CV)
    etaiov_cl_2 ~ fix(0.0296) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0296) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0296) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0296) # same variance as occasion 1

    propSd <- 0.1389; label("Proportional residual error (fraction)") # Table S4 'proportional' 0.0193 read as a variance: sqrt(0.0193) = 0.1389 (the printed 1.9 %CV is 100 x 0.0193; see vignette)
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    # ln(dose) for perampanel doses >= 1 mg; 0 for placebo.
    ldose_perampanel <- log(max(DOSE_PERAMPANEL_MGD, 1))
    tvcl <- exp(lcl) * (1 + e_conmed_cbz_cl * CONMED_CBZ +
      e_conmed_pb_cl * CONMED_PB + e_conmed_vpa_cl * CONMED_VPA +
      e_dose_perampanel_cl * ldose_perampanel)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_LTG_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
