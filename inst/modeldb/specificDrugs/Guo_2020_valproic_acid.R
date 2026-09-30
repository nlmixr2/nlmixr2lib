Guo_2020_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for total serum valproic acid (VPA) in Chinese adult inpatients with seizures receiving oral (immediate-release tablet / solution) or intravenous VPA (Guo 2020). Clearance falls with serum albumin through a power term centred on the 38.7 g/L cohort median and is lower in CYP2C19 non-extensive metabolizers (*1/*2, *1/*3, *2/*2, *2/*3, *3/*3) than in *1/*1 extensive metabolizers; the volume of distribution is larger in men than in women. Ka was FIXED at 2.38 1/h from the literature because the sparse steady-state sampling carried no absorption information. Residual error is proportional. Fit with FOCE-ELS in Phoenix NLME to 98 concentrations from 60 patients."
  reference <- "Guo J, Huo Y, Li F, Li Y, Guo Z, Han H, Zhou Y. Impact of gender, albumin, and CYP2C19 polymorphisms on valproic acid in Chinese patients: a population pharmacokinetic model. J Int Med Res. 2020;48(8):300060520952281. doi:10.1177/0300060520952281. PMCID PMC7469748. Final-model equations 3-6 and Table 2; cohort demographics from Table 1."
  vignette <- "Guo_2020_valproic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL as (ALB/38.7)^-1.06 (equations 3-4; Table 2 row 'f ALB-CL' = -1.06). 38.7 g/L is the cohort median albumin (Table 1, 'Albumin (g/L)' median). The negative exponent means clearance of total VPA rises as albumin falls, which the authors attribute to the higher unbound fraction of this highly albumin-bound drug (Discussion). Cohort range 25.6-53.6 g/L.",
      source_name = "ALB"
    ),
    CYP2C19_NON_EM = list(
      description = "CYP2C19 non-extensive-metabolizer indicator (1 = *1/*2, *1/*3, *2/*2, *2/*3 or *3/*3; 0 = *1/*1)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Enters CL as exp(-0.45 * CYP2C19_NON_EM) (equation 4; Table 2 row 'f CYP2C19-CL' = -0.45), i.e. 36% lower clearance than in *1/*1 extensive metabolizers. The paper's Figure 4 legend codes the covariate as CYP2C19 = 1 for *1/*1 and 2 for the pooled loss-of-function group; the canonical indicator is CYP2C19_NON_EM = (source code - 1). Only *2 (rs4244285) and *3 (rs4986893) were used for grouping; *17 was genotyped but not used in the grouping. 36 of 60 patients were *1/*1 (Table 1).",
      source_name = "CYP2C19"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "The paper codes GNDR = 1 for male and 0 for female (Figure 4 legend) and puts the effect on the MALE group: V = 22.15 * exp(0.78) for men and 22.15 for women (equations 5-6; Table 2 row 'f GNDR-V' = 0.78). Encoded here as exp(e_sexm_vc * (1 - SEXF)) so women are the typical value, exactly as printed. 16 of 60 patients were female (Table 1).",
      source_name = "GNDR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 60,
    n_studies = 1,
    n_observations = 98,
    age_range = "22-88 years (Table 1); adults 18 years or older by inclusion criterion.",
    age_median = "59 years (mean 60, SD 11.8; Table 1).",
    weight_range = "40-90 kg (Table 1).",
    weight_median = "66.5 kg (mean 66.5, SD 12.1; Table 1).",
    sex_female_pct = 26.7,
    race_ethnicity = c(Asian = 100),
    disease_state = "Inpatients with seizures on valproic acid therapeutic drug monitoring. Hepatic dysfunction, pregnancy, traditional Chinese medicine and co-medication known to change VPA concentrations (phenobarbital, carbamazepine) were exclusion criteria, although Table 1 lists one patient each on carbamazepine, lamotrigine, meropenem and imipenem.",
    dose_range = "Standard regimens were 500 mg oral twice daily (immediate-release tablets / solution) or 400 mg intravenous twice daily; individual doses 200-1200 mg, median 500 mg (Table 1).",
    regions = "China (single centre: The General Hospital of Taiyuan Iron & Steel (Group) Corporation, Taiyuan, Shanxi; January-December 2018).",
    albumin = "Serum albumin 25.6-53.6 g/L, median 38.7 g/L (mean 38.9, SD 6.4; Table 1).",
    cyp2c19 = "36 *1/*1 extensive metabolizers and 24 carriers of a *2 or *3 loss-of-function allele (Table 1).",
    notes = "Prospective sparse steady-state therapeutic-drug-monitoring sampling (98 concentrations from 60 patients), measured by homogeneous enzyme immunoassay (Roche Cobas c311; calibration range 2.8-150 ug/mL)."
  )

  ini({
    # Structural parameters - final-model estimates, Guo 2020 Table 2 and equations 3-6.
    # Ka FIXED: "it is necessary to fix Ka to 2.38 hour-1, in accordance with the
    # references" (Methods, Base model); Table 2 prints '2.38 (FIXED)'.
    lka <- fixed(log(2.38))
    label("Absorption rate constant (1/h)") # Table 2: Ka = 2.38 /hour (FIXED)
    lcl <- log(0.64)
    label("Clearance for a CYP2C19 *1/*1 patient with albumin 38.7 g/L (L/h)") # Table 2: CL = 0.64 L/hour (RSE 7.37%); equation 3
    lvc <- log(22.15)
    label("Volume of distribution for a female patient (L)") # Table 2: V = 22.15 L (RSE 10.68%); equation 5

    # Covariate effects (Table 2; minus signs confirmed from the 95% CIs,
    # -1.87 to -0.26 and -0.66 to -0.25).
    e_alb_cl <- -1.06
    label("Power exponent of albumin/38.7 on CL (unitless)") # Table 2: f ALB-CL = -1.06 (RSE 38.11%)
    e_cyp2c19nonem_cl <- -0.45
    label("Log-scale effect of CYP2C19 non-extensive-metabolizer status on CL (unitless)") # Table 2: f CYP2C19-CL = -0.45 (RSE 22.46%); equation 4
    e_sexm_vc <- 0.78
    label("Log-scale effect of male sex on V (unitless)") # Table 2: f GNDR-V = 0.78 (RSE 20.98%); equation 6

    # IIV. Table 2 reports IIV as CV% of log-normal etas (P = P_TV * exp(eta),
    # equation 1); variance = log(CV^2 + 1).
    etalcl ~ 0.18186 # Table 2: IIV CL = 44.66 CV% -> log(1 + 0.4466^2)
    etalvc ~ 0.091814 # Table 2: IIV V = 31.01 CV% -> log(1 + 0.3101^2)

    # Residual error: Cobs = Cpred * (1 + eps) (equation 2).
    propSd <- 0.1175
    label("Proportional residual error (fraction)") # Table 2: Proportional error = 11.75 (%)
  })

  model({
    ka <- exp(lka)

    # Equations 3-4: CL = 0.64 * (ALB/38.7)^-1.06 * exp(-0.45 * nonEM) * exp(eta)
    cl <- exp(lcl + etalcl) *
      (ALB / 38.7)^e_alb_cl *
      exp(e_cyp2c19nonem_cl * CYP2C19_NON_EM)

    # Equations 5-6: V = 22.15 * exp(0.78 * male) * exp(eta)
    vc <- exp(lvc + etalvc) * exp(e_sexm_vc * (1 - SEXF))

    kel <- cl / vc

    # Figure 1: oral doses enter the depot, intravenous doses enter central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, V in L -> mg/L (= ug/mL, the paper's concentration unit).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
