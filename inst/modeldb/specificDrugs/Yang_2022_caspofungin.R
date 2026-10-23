Yang_2022_caspofungin <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous caspofungin in hospitalised Chinese",
    "adults with confirmed or probable candidiasis (ICU, solid-organ transplant, haematology and",
    "respiratory patients; the 'general PPK model' of Yang 2022). Clearance and volume each carry",
    "two dichotomised covariates, hypoalbuminaemia (serum albumin < 35 g/L) and body weight >= 70",
    "kg, entering as (1 + theta x indicator) factors. Log-normal IIV on CL and V and a",
    "proportional residual error. A separately fitted ICU-subgroup model from the same paper is",
    "Yang_2022_caspofungin_icu.",
    sep = " "
  )
  reference <- paste(
    "Yang Q, Zhang T, Zhang Y, Sun D, Zheng X, Du Q, Wang X, Cheng X, Xing J, Dong Y. The",
    "recommended dosage regimen for caspofungin in patients with higher body weight or",
    "hypoalbuminaemia will result in low exposure: five years of data based on a population",
    "pharmacokinetic model and Monte-Carlo simulations. Front Pharmacol. 2022;13:993330.",
    "doi:10.3389/fphar.2022.993330. PMCID: PMC9669616.",
    sep = " "
  )
  vignette <- "Yang_2022_caspofungin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "caspofungin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only DICHOTOMISED: model() forms alb_low <- (ALB < 35), the paper's ALB* (Table 3",
        "footnote c: ALB* = 1 if albumin < 35 g/L, otherwise 0). The reference (ALB* = 0) is",
        "albumin >= 35 g/L. Tested as a continuous covariate it did not reduce the OFV",
        "(Discussion). Cohort 35.6 +/- 11.6 g/L (15.1-64.5); 161 of 299 patients < 35 g/L",
        "(Table 1, Results 3.1.1).",
        sep = " "
      ),
      source_name = "ALB"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Used only DICHOTOMISED: model() forms wt_high <- (WT >= 70), the paper's WT* (Table 3",
        "footnote c: WT* = 1 if weight >= 70 kg, otherwise 0). The reference (WT* = 0) is weight",
        "< 70 kg. Methods 2.2 states the same < 70 / >= 70 kg split; the dosing simulations",
        "(Figure 5, Tables S1-S2) label the groups <= 70 / > 70 kg. Cohort 62.3 +/- 11.5 kg",
        "(30-100; Table 1).",
        sep = " "
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_CSA = list(
      description = "Concomitant cyclosporine",
      units = "(binary)",
      type = "categorical",
      notes = "Significant on CL in the forward step (Table 2, dOFV -20.41) and kept in the full model (Supplementary 3.2.1) but removed in backward elimination (dOFV -9.20 < 10.83). 24 of 299 patients (Table 1)."
    ),
    CONMED_METHYLPREDNISOLONE = list(
      description = "Concomitant methylprednisolone",
      units = "(binary)",
      type = "categorical",
      notes = "Significant on CL in the forward step (Table 2, dOFV -19.35) but removed in backward elimination (dOFV -5.89). 110 of 299 patients (Table 1)."
    ),
    DIS_CRITILL = list(
      description = "ICU patient group indicator",
      units = "(binary)",
      type = "categorical",
      notes = "Significant on CL in the forward step (Table 2, dOFV -8.49) but removed in backward elimination (dOFV -5.27). 51 of 299 patients were ICU patients (Table 1)."
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy",
      units = "(binary)",
      type = "categorical",
      notes = "Collected (Table 1, 50 of 299 patients) but not reported as a significant covariate in Table 2."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Collected and tested (Methods 2.2) but not retained. Cohort 44 +/- 17 years (18-99; Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "categorical",
      notes = "Collected and tested (Methods 2.2) but not retained. 92 of 299 patients female (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 299L,
    n_studies = 1L,
    n_concentrations = 921L,
    age_range = "18-99 years",
    age_mean = "44 years (SD 17)",
    weight_range = "30-100 kg",
    weight_mean = "62.3 kg (SD 11.5)",
    sex_female_pct = 30.8,
    race_ethnicity = "Chinese (single centre in Xi'an, China)",
    disease_state = paste(
      "Hospitalised adults with confirmed or probable candidiasis (EORTC/MSG criteria) treated",
      "with caspofungin: 51 ICU, 163 solid-organ transplant, 42 haematology and 43 other",
      "(respiratory department) patients. 161 of 299 had hypoalbuminaemia (< 35 g/L).",
      sep = " "
    ),
    dose_range = paste(
      "1 h intravenous infusion; mostly 70 mg loading then 50 mg daily (70/50 mg); 70/35 mg in",
      "Child-Pugh B hepatic impairment; 70/70 mg for weight > 75 kg.",
      sep = " "
    ),
    regions = "China (First Affiliated Hospital of Xi'an Jiaotong University)",
    notes = paste(
      "Retrospective TDM data collected June 2014 - June 2019. 921 plasma concentrations (median 3",
      "per patient): 242 troughs (22-24 h post-dose), 247 at 0-12 h and 432 at 12-22 h. Assay",
      "LC-MS/MS. Mycophenolate mofetil, intravenous infusion pump vs drip, surgery and",
      "solid-organ transplant status were also significant in the forward step on CL (Table 2)",
      "but were removed in backward elimination.",
      sep = " "
    )
  )

  ini({
    # Final model (Table 3, footnote c):
    #   CL = theta1 * (1 + theta3 * ALB*) * (1 + theta5 * WT*) * exp(eta1)
    #   V  = theta2 * (1 + theta4 * ALB*) * (1 + theta6 * WT*) * exp(eta2)
    lcl <- log(0.32); label("Clearance for albumin >= 35 g/L and weight < 70 kg (L/h)") # Table 3, theta1 = 0.32 L/h (bootstrap 95% CI 0.29-0.34)
    lvc <- log(13.31); label("Volume of distribution for albumin >= 35 g/L and weight < 70 kg (L)") # Table 3, theta2 = 13.31 L (bootstrap 95% CI 10.18-16.42)
    e_alb_cl <- 0.46; label("Fractional change in CL for albumin < 35 g/L (unitless)") # Table 3, theta3 = 0.46 (bootstrap 95% CI 0.33-0.59)
    e_alb_vc <- 0.49; label("Fractional change in V for albumin < 35 g/L (unitless)") # Table 3, theta4 = 0.49 (bootstrap 95% CI 0.25-0.73)
    e_wt_cl <- 0.98; label("Fractional change in CL for weight >= 70 kg (unitless)") # Table 3, theta5 = 0.98 (bootstrap 95% CI 0.30-1.65)
    e_wt_vc <- 0.24; label("Fractional change in V for weight >= 70 kg (unitless)") # Table 3, theta6 = 0.24 (bootstrap 95% CI 0.10-0.37)

    # IIV: Table 3 'Inter-individual variability (%)' read as CV%, omega^2 = log(CV^2 + 1)
    etalcl ~ 0.081825 # Table 3, IIV CL 29.2 -> log(0.292^2 + 1)
    etalvc ~ 0.300448 # Table 3, IIV V 59.2 -> log(0.592^2 + 1)

    propSd <- 0.193; label("Proportional residual error (fraction)") # Table 3, residual variability 'Proportional (%)' = 19.3
  })
  model({
    # Dichotomised covariates (Table 3 footnote c)
    alb_low <- (ALB < 35)
    wt_high <- (WT >= 70)

    cl <- exp(lcl + etalcl) * (1 + e_alb_cl * alb_low) * (1 + e_wt_cl * wt_high)
    vc <- exp(lvc + etalvc) * (1 + e_alb_vc * alb_low) * (1 + e_wt_vc * wt_high)

    kel <- cl / vc

    # Caspofungin is given as a 1 h intravenous infusion into central.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
