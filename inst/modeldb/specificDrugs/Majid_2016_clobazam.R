Majid_2016_clobazam <- function() {
  description <- paste0(
    "Steady-state apparent-clearance population PK model for clobazam given ",
    "with adjunctive perampanel (or placebo) in adolescents and adults with ",
    "refractory partial-onset seizures (Majid 2016, pooled Phase III; ",
    "n = 114). Css = R0 / (CL/F), R0 = daily dose / 24; no volume or ",
    "absorption term. CL/F = basal * (1 + sex + weight + phenytoin + ",
    "perampanel-concentration terms), all additive inside one relative ",
    "bracket; exponential IIV and five-occasion IOV on CL/F; proportional ",
    "residual error. The covariate equation is not printed by the paper; ",
    "this form reproduces Table 2 exactly -- see the vignette."
  )
  reference <- paste(
    "Majid O, Laurenza A, Ferry J, Hussein Z. Impact of perampanel on",
    "pharmacokinetics of concomitant antiepileptics in patients with",
    "partial-onset seizures: pooled analysis of clinical trials.",
    "Br J Clin Pharmacol. 2016;82(2):422-430. doi:10.1111/bcp.12951.",
    "Parameter estimates from Supplemental Table S2.",
    sep = " "
  )
  vignette <- "Majid_2016_perampanel_concomitant_aeds"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_CLB_MGD = list(
      description = "Total daily clobazam dose (the drug whose concentration is modelled).",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Input rate R0 = DOSE_CLB_MGD / 24 (mg/h) in Css = R0 / CL. Dose levels are not reported by the paper.",
      source_name = "clobazam daily dose (R0 x 24)"
    ),
    SEXF = list(
      description = "Sex indicator; 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste0(
        "Majid 2016 Table S2 'Sex effect' -0.387, additive inside the ",
        "relative bracket. The female coding is fixed by Table 2 (males ",
        "3.27, females 1.43 L/h): coding male = 1 cannot reproduce it. ",
        "48.9% male in the clobazam cohort (Table 1)."
      ),
      source_name = "sex"
    ),
    WT = list(
      description = "Body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Majid 2016 Table S2 'Weight effect' -0.312, entered as ",
        "e_wt_cl * (WT / 70) inside the relative bracket (uncentred). The ",
        "uncentred linear form is the one that reproduces Table 2 (it needs ",
        "WT/WTref = 1.003 at the Table 2 typical value; a power or centred ",
        "form cannot reach 3.27 L/h from 4.76 L/h at any plausible weight). ",
        "The reference weight is not printed; 70 kg is used. The bracket ",
        "becomes negative for females heavier than about 137 kg, outside ",
        "the cohort (mean 67.9, SD 17.3 kg)."
      ),
      source_name = "weight"
    ),
    CONMED_PHT = list(
      description = "Concomitant phenytoin indicator; 1 = co-administered, 0 = not.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes = "Majid 2016 Table S2 'Phenytoin effect' 0.545, additive inside the relative bracket.",
      source_name = "phenytoin co-administration"
    ),
    CP_PERAMPANEL_NGML = list(
      description = "Perampanel plasma concentration (mean of the two samples taken at the visit); 0 for placebo and before perampanel start.",
      units = "ng/mL",
      type = "continuous",
      reference_category = "0 ng/mL (no perampanel)",
      notes = paste0(
        "Majid 2016 Table S2 'Perampanel concentration effect' 0.0000678, ",
        "additive inside the relative bracket (the paper: 'proportional to ",
        "perampanel concentrations and to clobazam clearance'). The data set ",
        "carried the per-visit mean of the two perampanel samples ",
        "(Supplement section A). The concentration unit is not printed; ",
        "ng/mL is assumed. Table 2's +0.09 / +0.11 L/h increments at 8 / 12 ",
        "mg (identical in both sexes, as this form requires) imply about ",
        "280 / 340 ng/mL."
      ),
      source_name = "perampanel concentration"
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
    n_subjects = 114L,
    n_studies = 3L,
    age_range = "12 years and older",
    age_median = "mean 32.8 (SD 12.9) years",
    weight_median = "mean 67.9 (SD 17.3) kg",
    sex_female_pct = 50.9,
    race_ethnicity = c(White = 66.7, Asian = 27.2, Chinese = 2.6, Black = 0, Other = 3.5),
    disease_state = "Refractory partial-onset seizures, uncontrolled despite >= 2 AEDs in the previous 2 years; up to three stable concomitant AEDs.",
    dose_range = "Adjunctive perampanel 2-12 mg once daily at bedtime (79 subjects) or placebo (35); clobazam dose not reported.",
    regions = "Multinational (three Phase III trials).",
    notes = "Majid 2016 Table 1 (CLB row); Table 2 gives 56 males (39 on perampanel) and 58 females (40 on perampanel)."
  )

  ini({
    # Majid 2016 Supplemental Table S2 (clobazam final model).
    lcl <- log(4.76); label("Basal clobazam CL/F (L/h)") # Table S2 'Clobazam basal clearance (L/h)' 4.76
    e_sexf_cl <- -0.387; label("Additive relative change in CL/F for females (unitless)") # Table S2 'Sex effect' -0.387
    e_wt_cl <- -0.312; label("Additive relative change in CL/F per unit WT/70 (unitless)") # Table S2 'Weight effect' -0.312
    e_conmed_pht_cl <- 0.545; label("Additive relative change in CL/F with concomitant phenytoin (unitless)") # Table S2 'Phenytoin effect' 0.545
    e_cp_perampanel_cl <- 0.0000678; label("Additive relative change in CL/F per ng/mL of perampanel (1/(ng/mL))") # Table S2 'Perampanel concentration effect' 0.0000678 (text: 95% CI 0.0000249-0.0001107)

    etalcl ~ 0.497 # Table S2 'IIV' 0.497 (70.5 %CV)
    etaiov_cl_1 ~ 0.0258 # Table S2 'IOV' 0.0258 (16.1 %CV)
    etaiov_cl_2 ~ fix(0.0258) # same variance as occasion 1
    etaiov_cl_3 ~ fix(0.0258) # same variance as occasion 1
    etaiov_cl_4 ~ fix(0.0258) # same variance as occasion 1
    etaiov_cl_5 ~ fix(0.0258) # same variance as occasion 1

    propSd <- 0.2429; label("Proportional residual error (fraction)") # Table S2 'proportional' 0.059 read as a variance: sqrt(0.059) = 0.2429 (the printed 5.9 %CV is 100 x 0.059; see vignette)
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
      e_conmed_pht_cl * CONMED_PHT + e_cp_perampanel_cl * CP_PERAMPANEL_NGML)
    cl <- tvcl * exp(etalcl + iov_cl)

    # Css = R0 / CL in mg/L, x 1000 for ng/mL.
    Cc <- 1000 * (DOSE_CLB_MGD / 24) / cl
    Cc ~ prop(propSd)
  })
}
