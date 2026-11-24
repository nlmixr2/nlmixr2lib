Majid_2016_oxcarbazepine <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for oxcarbazepine ",
    "(parent drug; the active monohydroxy metabolite was not measured) ",
    "given with adjunctive perampanel (or placebo) in adolescents and ",
    "adults with refractory partial-onset seizures (Majid 2016, pooled ",
    "Phase III; n = 200). Css = R0 / (CL/F), R0 = daily dose / 24; no ",
    "volume or absorption term. CL/F = basal * (1 + sex + phenytoin + ",
    "perampanel yes/no terms); perampanel lowers CL/F by 26% of the basal ",
    "value. Exponential IIV and five-occasion IOV on CL/F; proportional ",
    "residual error. The covariate equation is not printed by the paper; ",
    "this form reproduces all eight Table 2 cells -- see the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S6.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    DOSE_OXC_MGD = list(
      description = "Total daily oxcarbazepine dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_OXC_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "oxcarbazepine daily dose (R0 x 24)"
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Majid 2016 Table S6 'Sex effect' -0.263, additive inside the relative bracket (Table 2: males 63.4, females 46.7 L/h).",
      source_name = "sex"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Majid 2016 Table S6 'Phenytoin effect' 0.823, additive inside the relative bracket (Table 2: males with PHT 115.6 L/h).",
      source_name = "phenytoin co-administration"
    ),
    CONMED_PERAMPANEL = list(
      description = "Concomitant perampanel indicator; 1 = on perampanel (any dose), 0 = placebo or before perampanel start.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no perampanel)",
      notes = "Majid 2016 Table S6 'Perampanel effect' -0.261, additive inside the relative bracket, 'independent of perampanel dose or concentration' (Results). Because the term subtracts 0.261 x basal, the relative reduction is 26% in males without phenytoin and 35.3%, 14.4% and 16.8% in the other three groups (Table 3).",
      source_name = "perampanel (yes/no)"
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
    n_subjects = 200L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 33.8 (SD 13.8) years",
    weight_median = "mean 73.4 (SD 17.6) kg",
    sex_female_pct = 51.0,
    race_ethnicity = c(White = 78.0, Asian = 11.0, Chinese = 5.0, Black = 3.0, Other = 3.0),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (133 subjects) or placebo; oxcarbazepine dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (OXC row); Table 2 splits the cohort into males / females without phenytoin (91 / 96) and with phenytoin (7 / 6). The goodness-of-fit plots of this model showed over-prediction of low concentrations (Supplement)."
  )

  ini({
    # Majid 2016 Supplemental Table S6 (oxcarbazepine final model).
    lcl <- log(63.4); label("Basal oxcarbazepine CL/F, male, no phenytoin, no perampanel (L/h)") # Table S6 'Oxcarbazepine basal clearance (L/h)' 63.4
    e_sexf_cl <- -0.263; label("Additive relative change in CL/F for females (unitless)") # Table S6 'Sex effect' -0.263
    e_conmed_pht_cl <- 0.823; label("Additive relative change in CL/F with concomitant phenytoin (unitless)") # Table S6 'Phenytoin effect' 0.823
    e_conmed_perampanel_cl <- -0.261; label("Additive relative change in CL/F with perampanel (unitless)") # Table S6 'Perampanel effect' -0.261

    etalcl ~ 0.305 # Table S6 'IIV' 0.305 (55.2 %CV)
    etaiov_cl_1 ~ 0.195 # Table S6 'IOV' 0.195 (44.2 %CV)
    etaiov_cl_2 ~ fix(0.195) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.195) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.195) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.195) # same variance as occasion 1

    propSd <- 0.638; label("Proportional residual error (fraction)") # Table S6 'proportional' 0.407 read as a variance: sqrt(0.407) = 0.638 (the printed 40.7 %CV is 100 x 0.407; see vignette)
  })

  model({
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5

    tvcl <- exp(lcl) * (1 + e_sexf_cl * SEXF + e_conmed_pht_cl * CONMED_PHT +
      e_conmed_perampanel_cl * CONMED_PERAMPANEL)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL; mg/L = ug/mL.
    Cc <- (DOSE_OXC_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
