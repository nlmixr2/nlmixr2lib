Liu_2022_tacrolimus_cyp3a5 <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for intravenous ",
    "(continuous 24 h infusion) and oral tacrolimus in the CYP3A5-genotyped ",
    "subpopulation of pediatric hematopoietic stem cell transplant (HSCT) ",
    "recipients (Liu 2022, model 2, n = 24): the model 1 structure ",
    "re-estimated with an added exponential CYP3A5 expresser (*1 carrier) ",
    "effect on CL; ka fixed at 4.48 1/h with an estimated oral ",
    "bioavailability."
  )
  reference <- paste0(
    "Liu XL, Guan YP, Wang Y, Huang K, Jiang FL, Wang J, Yu QH, Qiu KF, ",
    "Huang M, Wu JY, Zhou DH, Zhong GP, Yu XX. Population Pharmacokinetics ",
    "and Initial Dosage Optimization of Tacrolimus in Pediatric Hematopoietic ",
    "Stem Cell Transplant Patients. Front Pharmacol. 2022;13:891648. ",
    "doi:10.3389/fphar.2022.891648"
  )
  vignette <- "Liu_2022_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power effect on CL, (WT / 17.8)^0.48 (Liu 2022 Eq. 6). The centring ",
        "value 17.8 kg is the one printed in Eq. 6 (the same as model 1)."
      ),
      source_name = "WT"
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-varying. Liu 2022 reports hematocrit as a volume fraction; the ",
        "canonical HCT column is percent, so the centring values are rescaled ",
        "by 100 inside model(): CL uses (HCT / 29.6)^-1.06 (Eq. 6 prints ",
        "Hct/0.296) and V uses (HCT / 28.9)^-0.60 (Eq. 7 prints Hct/0.289). ",
        "Pass HCT in percent (e.g. 28.9), not as a fraction."
      ),
      source_name = "Hct"
    ),
    CONMED_AZOLE = list(
      description = "Concomitant azole antifungal therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes = paste0(
        "Time-varying. Pooled voriconazole / itraconazole / posaconazole ",
        "variable 'CZ'. Effect exp(-0.39 * CONMED_AZOLE) on CL (Eq. 6, ",
        "Table 3 subpopulation theta CZ,CL). Only 6 of the 24 genotyped ",
        "patients (all CYP3A5*1/*3, all on IV tacrolimus) received an azole ",
        "(Results and Discussion)."
      ),
      source_name = "CZ"
    ),
    CONMED_CASPOFUNGIN = list(
      description = "Concomitant caspofungin therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant caspofungin)",
      notes = "Time-varying. Effect exp(0.04 * CONMED_CASPOFUNGIN) on CL (Eq. 6, Table 3 subpopulation theta CPFG,CL).",
      source_name = "CPFG"
    ),
    POD = list(
      description = "Post-transplant day",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-varying days since the stem cell transplant. Enters only as ",
        "the indicator (POD >= 28), the paper's PTD = 2 category (footnote ",
        "to Eqs. 4-7), with effect exp(-0.33 * indicator) on CL (Eq. 6)."
      ),
      source_name = "PTD"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser status (carrier of at least one CYP3A5*1 allele)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 nonexpresser)",
      notes = paste0(
        "Time-fixed germline genotype (rs776746). Liu 2022 codes CYP = 1 for ",
        "*1/*1 or *1/*3 and CYP = 0 for *3/*3 (Methods); same orientation ",
        "as the canonical column, no value transformation. Effect ",
        "exp(0.32 * CYP3A5_EXPR) on CL, a 1.38-fold higher CL in expressers ",
        "(Eq. 6; Table 3 row 'theta CYP3A5*3, CL'). Genotyped cohort: one ",
        "*1/*1, 11-12 *1/*3 and 11-12 *3/*3 (the Results give both splits)."
      ),
      source_name = "CYP"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    disease_state = paste0(
      "CYP3A5-genotyped subpopulation of the Liu 2022 pediatric allogeneic ",
      "HSCT cohort (whole-blood EDTA sample genotyped by MassARRAY ",
      "MALDI-TOF for rs776746)."
    ),
    dose_range = paste0(
      "Tacrolimus continuous 24 h IV infusion followed by oral tacrolimus ",
      "(the full cohort received median 0.024 mg/kg/day IV and 0.057 ",
      "mg/kg/day oral)."
    ),
    regions = "China (Sun Yat-sen Memorial Hospital, Guangzhou)",
    cyp3a5_genotypes = paste0(
      "*1/*1 n = 1; *1/*3 n = 12 and *3/*3 n = 11 per the Results ",
      "'Participants' Characteristics' paragraph, or *1/*3 n = 11 and *3/*3 ",
      "n = 12 per the 'Pharmacokinetic Modeling' paragraph; Supplementary ",
      "Table S2 lists 1 TT, 11 CT and 12 CC at rs776746."
    ),
    notes = paste0(
      "Subset of the 86-patient model 1 cohort (Liu_2022_tacrolimus) with ",
      "CYP3A5 genotype available. Age, weight and sex of the subset are not ",
      "reported separately. Retrospective TDM trough data measured by EMIT; ",
      "FOCE-I in Phoenix NLME 7.0."
    )
  )

  ini({
    # Structural parameters (Liu 2022 Table 3, 'Subpopulation: Pharmacogenomic Dataset (n = 24)' column)
    lka <- fixed(log(4.48)); label("Absorption rate constant (1/h)") # Table 3 subpopulation ka = 4.48 h-1, no RSE (fixed from Jusko 1995 / Wallin 2009, Results)
    lcl <- log(2.41); label("Clearance at WT 17.8 kg, hematocrit 29.6% in CYP3A5*3/*3 (L/h)") # Table 3 subpopulation CL = 2.41 L/h (RSE 10.87%); Eq. 6
    lvc <- log(92.9); label("Volume of distribution at hematocrit 28.9% (L)") # Table 3 subpopulation V = 92.9 L (RSE 24.07%); Eq. 7
    lfdepot <- log(0.25); label("Oral bioavailability (fraction)") # Table 3 subpopulation F = 0.25 (RSE 20.12%)

    # Covariate effects (Liu 2022 Table 3 subpopulation column and Eqs. 6-7)
    e_wt_cl <- 0.48; label("Power exponent of body weight on CL (unitless)") # Table 3 subpopulation theta WT,CL = 0.48 (RSE 26.00%)
    e_hct_cl <- -1.06; label("Power exponent of hematocrit on CL (unitless)") # Table 3 subpopulation theta Hct,CL = -1.06 (RSE 18.08%)
    e_conmed_azole_cl <- -0.39; label("Exponential effect of concomitant azole antifungal on CL (unitless)") # Table 3 subpopulation theta CZ,CL = -0.39 (RSE 32.09%)
    e_conmed_caspofungin_cl <- 0.04; label("Exponential effect of concomitant caspofungin on CL (unitless)") # Table 3 subpopulation theta CPFG,CL = 0.04 (RSE 37.81%)
    e_pod_cl <- -0.33; label("Exponential effect of post-transplant day >= 28 on CL (unitless)") # Table 3 subpopulation theta PTD,CL = -0.33 (RSE 26.75%)
    e_cyp3a5_expr_cl <- 0.32; label("Exponential effect of CYP3A5 expresser (*1 carrier) on CL (unitless)") # Table 3 subpopulation theta CYP3A5*3,CL = 0.32 (RSE 27.08%); Eq. 6 applies it to 'CYP3A5*1 carriers'
    e_hct_vc <- -0.60; label("Power exponent of hematocrit on V (unitless)") # Table 3 subpopulation theta Hct,V = -0.60 (RSE 26.46%)

    # IIV: exponential (Methods Eq. 1); Table 3 reports omega^2 (variances)
    etalcl ~ 0.02 # Table 3 subpopulation omega2 CL = 0.02 (RSE 27.97%)
    etalvc ~ 0.66 # Table 3 subpopulation omega2 V = 0.66 (RSE 43.27%)
    etalfdepot ~ 0.25 # Table 3 subpopulation omega2 F = 0.25 (RSE 51.38%)

    # Residual error
    propSd <- 0.368; label("Proportional residual error (fraction)") # Table 3 subpopulation sigma proportional = 36.8% (RSE 9.90%)
  })

  model({
    # Post-transplant day indicator: PTD = 2 (>= 28 days) vs PTD = 1 (Eqs. 4-7 footnote)
    ptd_late <- POD >= 28

    # Individual parameters (Eqs. 6 and 7); HCT in percent, so the printed
    # fraction centring values 0.296 and 0.289 become 29.6 and 28.9.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (WT / 17.8)^e_wt_cl *
      (HCT / 29.6)^e_hct_cl *
      exp(e_conmed_azole_cl * CONMED_AZOLE) *
      exp(e_conmed_caspofungin_cl * CONMED_CASPOFUNGIN) *
      exp(e_pod_cl * ptd_late) *
      exp(e_cyp3a5_expr_cl * CYP3A5_EXPR)
    vc <- exp(lvc + etalvc) * (HCT / 28.9)^e_hct_vc
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc

    # Oral doses go to depot; the continuous 24 h IV infusion goes to central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- fdepot

    # Dose in mg and V in L give mg/L; x1000 gives ng/mL (ug/L).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
