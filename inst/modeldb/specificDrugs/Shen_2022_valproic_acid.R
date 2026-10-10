Shen_2022_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for total serum valproic acid (VPA) in Chinese children with epilepsy on oral syrup or conventional tablets (Shen 2022). Apparent clearance rises with age as a power function centred on 5 years and differs by ABCB1 rs3789243 genotype (AG 4.7% lower, GG 8% higher than the AA wild type); apparent volume has no covariates. Ka was FIXED at 1.9 1/h from the literature because the steady-state trough-only sampling carried no absorption information. IIV on CL/F only; additive residual error. Fit with FOCE-I in NONMEM 7.5 to 376 steady-state trough concentrations from 103 patients."
  reference <- "Shen X, Chen X, Lu J, Chen Q, Li W, Zhu J, He Y, Guo H, Xu C, Fan X. Pharmacogenetics-based population pharmacokinetic analysis and dose optimization of valproic acid in Chinese southern children with epilepsy: Effect of ABCB1 gene polymorphism. Front Pharmacol. 2022;13:1037239. doi:10.3389/fphar.2022.1037239. PMCID PMC9733833. Final-model equations 1-2 and Table 4; cohort demographics from Tables 1-2."
  vignette <- "Shen_2022_valproic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as (AGE/5)^0.357 (equation 1; Table 4 row 'Age on CL/F' = 0.357). The 5-year reference is printed in equation 1 and sits near the cohort mean of 5.30 years (Table 1). Cohort range 0.5-15 years. Age and body weight were strongly correlated (r = 0.915); age was retained because it lowered the OFV more (Table 3 models 2-3), so weight does not enter the model.",
      source_name = "Age"
    ),
    SNP_ABCB1_RS3789243_AG = list(
      description = "ABCB1 rs3789243 heterozygous AG genotype indicator (1 = AG, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Enters CL/F as 0.953^SNP_ABCB1_RS3789243_AG (equation 1; Table 4 row 'ABCB1 rs3789243 AG on CL/F' = 0.953), i.e. 4.7% lower CL/F than the AA wild type. Paired with SNP_ABCB1_RS3789243_GG; both are 0 for the AA reference group, which the paper calls the wild type (text below equation 1). 50 of 103 patients were AG (Table 2).",
      source_name = "ABCB1 AG"
    ),
    SNP_ABCB1_RS3789243_GG = list(
      description = "ABCB1 rs3789243 homozygous GG genotype indicator (1 = GG, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Enters CL/F as 1.08^SNP_ABCB1_RS3789243_GG (equation 1; Table 4 row 'ABCB1 rs3789243 GG on CL/F' = 1.08), i.e. 8% higher CL/F than the AA wild type. 39 of 103 patients were GG and 14 were AA (Table 2).",
      source_name = "ABCB1 GG"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Significant on CL/F in univariate forward inclusion (Table 3 model 2, dOFV -159.5) but not carried forward because of its 0.915 correlation with age, which lowered the OFV more."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 103,
    n_studies = 1,
    n_observations = 376,
    age_range = "0.5-15 years (Table 1); younger than 16 years by inclusion criterion.",
    age_mean = "5.30 years (SD 3.39; Table 1).",
    weight_range = "6.5-52.0 kg (Table 1).",
    weight_mean = "19.9 kg (SD 10.6; Table 1).",
    sex_female_pct = 45.6,
    race_ethnicity = c(Asian = 100),
    disease_state = "Children with epilepsy or an epileptic syndrome (ILAE classification) on valproic acid monotherapy or combination antiepileptic therapy; impaired hepatic or renal function was an exclusion criterion. Co-medications (percent of samples): levetiracetam 13.8%, oxcarbazepine 11.2%, topiramate 4.3%, clonazepam 3.7%, phenobarbital 3.2%, midazolam 2.7%, ibuprofen 1.6% (Table 1); none were retained as covariates.",
    dose_range = "Oral VPA two to three times daily as syrup (Depakine) or conventional tablets; 23.8 +/- 5.7 mg/kg/day, range 9.9-45.7 mg/kg/day (Table 1).",
    regions = "China (single centre: Shenzhen Baoan Women's and Children's Hospital, September 2016 - January 2022).",
    abcb1_rs3789243 = "AA 14, AG 50, GG 39 of 103 patients (Table 2; Hardy-Weinberg P = 0.748).",
    notes = "Retrospective steady-state morning trough therapeutic-drug-monitoring samples (at least 1 week on a stable regimen), total VPA measured by homogeneous enzyme immunoassay (Siemens Viva-E; LLOQ 1 mg/L). Observed concentrations 60.54 +/- 19.32 mg/L (range 14.67-110.99; Table 1)."
  )

  ini({
    # Structural parameters - final-model estimates, Shen 2022 Table 4 and equations 1-2.
    # Ka fixed: 'Ka was fixed at 1.9 h-1, in accordance with the references
    # (Ding et al., 2015)' (Section 2.4.1).
    lka <- fixed(log(1.9))
    label("Absorption rate constant (1/h)") # Section 2.4.1: Ka fixed at 1.9 h-1
    lcl <- log(0.214)
    label("Apparent clearance CL/F for a 5-year-old ABCB1 rs3789243 AA patient (L/h)") # Table 4: theta CL = 0.214 L/h (RSE 7.4%); equation 1
    lvc <- log(3.63)
    label("Apparent volume of distribution V/F (L)") # Table 4: theta V = 3.63 L (RSE 23.8%); equation 2

    # Covariate effects on CL/F (Table 4; equation 1).
    e_age_cl <- 0.357
    label("Power exponent of AGE/5 on CL/F (unitless)") # Table 4: theta Age = 0.357 (RSE 9.6%)
    e_snp_abcb1_rs3789243_ag_cl <- 0.953
    label("Multiplicative factor on CL/F for ABCB1 rs3789243 AG vs AA (unitless)") # Table 4: theta ABCB1 AG = 0.953 (RSE 4.7%)
    e_snp_abcb1_rs3789243_gg_cl <- 1.08
    label("Multiplicative factor on CL/F for ABCB1 rs3789243 GG vs AA (unitless)") # Table 4: theta ABCB1 GG = 1.08 (RSE 5.5%)

    # IIV on CL/F only (IIV on V/F was dropped for poor precision, Section 3.2).
    # Table 4 footnote: omega = coefficient of variation of IIV, so 0.169 is a
    # CV fraction and the log-normal variance is log(1 + 0.169^2).
    etalcl ~ 0.028161 # Table 4: omega CL = 0.169 (RSE 13.0%) -> log(1 + 0.169^2)

    # Residual error: additive (Section 3.2).
    addSd <- 11.9
    label("Additive residual error (mg/L)") # Table 4: sigma additive = 11.9 mg/L (RSE 5.9%)
  })

  model({
    ka <- exp(lka)

    # Equation 1: CL/F = 0.214 * (Age/5)^0.357 * 0.953^ABCB1AG * 1.08^ABCB1GG
    cl <- exp(lcl + etalcl) *
      (AGE / 5)^e_age_cl *
      e_snp_abcb1_rs3789243_ag_cl^SNP_ABCB1_RS3789243_AG *
      e_snp_abcb1_rs3789243_gg_cl^SNP_ABCB1_RS3789243_GG

    # Equation 2: V/F = 3.63
    vc <- exp(lvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg, V/F in L -> mg/L (the paper's concentration unit).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
