Komenkul_2026_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and",
    "linear elimination for oral immediate-release tacrolimus (whole-blood",
    "concentrations) in Thai adult liver transplant recipients during the",
    "first 6 months after transplantation (Komenkul 2026, Phoenix NLME",
    "FOCE-ELS). Apparent oral clearance CL/F carries a median-normalised",
    "power effect of the total daily tacrolimus dose (reference 4 mg/day)",
    "and of haemoglobin (reference 11 g/dL), linear fractional effects of",
    "fluconazole (-33%) and prednisolone (+8%) co-administration, and a",
    "four-level combined recipient/donor CYP3A5 genotype multiplier",
    "(expressor/expressor reference; 0.85 expressor recipient + non-expressor",
    "donor, 0.84 non-expressor recipient + expressor donor, 0.76",
    "non-expressor/non-expressor). Exponential IIV on CL/F and Vc/F",
    "(correlated) and on Vp/F; no IIV on ka or Q/F; proportional residual",
    "error."
  )
  reference <- paste(
    "Komenkul V, Chantharit P, Komolmit P, Nonthasoot B, Vorasittha A,",
    "Sanpavat A, Suksawatamnuay S, Radhakrishnan C, Wattanavijitkul T",
    "(2026). Precision Dosing of Tacrolimus in Liver Transplantation:",
    "Integrating Donor-Recipient CYP3A5 Pharmacogenomics and Drug",
    "Interactions. CPT Pharmacometrics Syst Pharmacol 15(9):e70339.",
    "doi:10.1002/psp4.70339."
  )
  vignette <- "Komenkul_2026_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Komenkul 2026 Methods 'Pharmacokinetic Sampling and Bioanalytical
  # Assays': whole-blood tacrolimus by CMIA (Abbott), LLOQ 1.2 ng/mL.
  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    DOSE_TAC_MGD = list(
      description = "Total daily tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Komenkul 2026 final-model equation: (DD/4)^0.22 on CL/F, DD = 'daily",
        "dose of TAC (mg/day)'; 4 mg/day is the Table 1 median (range 1-13).",
        "The authors' Phoenix code writes (DD/4000)^dCldDD, i.e. the dataset",
        "carried the dose in micrograms; the mg/day form here is the same",
        "ratio. Time-varying: updated whenever the dose is titrated. For the",
        "twice-daily regimen DOSE_TAC_MGD is twice the per-administration amt;",
        "keep the two consistent when simulating. The Discussion attributes",
        "the positive exponent partly to dose-adaptation bias in a",
        "trough-titrated TDM cohort (see the register entry)."
      ),
      source_name = "DD"
    ),
    HGB = list(
      description = "Blood hemoglobin concentration",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Komenkul 2026 final-model equation: (HB/11)^-0.15 on CL/F; Table 1",
        "median 10.9 g/dL (Q1-Q3 9.5-12.5, range 5.7-16.7). The equation and",
        "the authors' Phoenix code both normalise by 11 g/dL, which is used",
        "here. Time-varying in the source (post-operative anaemia)."
      ),
      source_name = "HB"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole (antifungal prophylaxis) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no fluconazole)",
      notes = paste(
        "Komenkul 2026: (1 + (-0.33) * FLU) on CL/F. Time-varying: all",
        "patients started fluconazole prophylaxis after transplantation and",
        "stopped it after the prophylactic course; Table 1 reports 59% of",
        "concentrations on fluconazole (76% in months 0-3, 17% in months 3-6)."
      ),
      source_name = "FLU"
    ),
    CONMED_PREDNISOLONE = list(
      description = "Concomitant prednisolone indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prednisolone)",
      notes = paste(
        "Komenkul 2026: (1 + 0.08 * PSL) on CL/F (the authors' code names the",
        "column CPRE / CPSL). Time-varying: prednisolone was part of the",
        "initial regimen and tapered off; Table 1 reports 88% of",
        "concentrations on prednisolone. Coded as an on/off indicator; the",
        "prednisolone dose did not enter the model."
      ),
      source_name = "PSL"
    ),
    CYP3A5_EXPR = list(
      description = "Recipient CYP3A5 expressor: 1 = at least one CYP3A5*1 allele (*1/*1 or *1/*3), 0 = CYP3A5*3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 non-expressor)",
      notes = paste(
        "Recipient germline rs776746 genotype (peripheral-blood DNA).",
        "Combined with CYP3A5_EXPR_DONOR inside model() to rebuild the paper's",
        "four-level COMB covariate (recipient/donor Ex/Ex reference, Ex/Nex,",
        "Nex/Ex, Nex/Nex). Table 1: 66 of 114 recipients (58%) are expressors."
      ),
      source_name = "COMB (recipient half)"
    ),
    CYP3A5_EXPR_DONOR = list(
      description = "Liver-graft donor CYP3A5 expressor: 1 = at least one CYP3A5*1 allele, 0 = CYP3A5*3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (donor CYP3A5*3/*3 non-expressor)",
      notes = paste(
        "Donor germline rs776746 genotype (DNA from archived FFPE donor liver",
        "tissue). Combined with CYP3A5_EXPR inside model() to rebuild the",
        "four-level COMB covariate. Table 1: 58 of 114 grafts (51%) came",
        "from expressor donors."
      ),
      source_name = "COMB (donor half)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods 'Pharmacokinetic Modeling') but not retained. Table 1 median 57 (18-75)."
    ),
    WT = list(
      description = "Actual body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened but not retained. Table 1 median 61.4 (34.6-102.0)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Excluded before stepwise selection for collinearity with body weight (r > 0.8). Table 1 median 22.5."
    ),
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = "Excluded before stepwise selection for collinearity with hemoglobin (r > 0.8). Table 1 median 33.6."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened but not retained. Table 1 median 0.98 (0.21-4.23)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened but not retained. Table 1 median 3.9 (1.6-5.5)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened but not retained. Table 1 median 0.8 (0.2-32.3)."
    ),
    DBIL = list(
      description = "Direct bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Excluded before stepwise selection for collinearity with total bilirubin (r > 0.8). Table 1 median 0.4."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Table 1 median 23 (3-4042)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested last in forward inclusion (Table S1 FW6, dOFV -2.50) and not retained."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Table 1 median 87 (28-1666)."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained. Table 1: 43 of 114 (37.7%) female."
    ),
    CONMED_MPA = list(
      description = "Concomitant mycophenolic acid",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained. Table 1: 87% of concentrations."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained. Table 1: 61% of concentrations."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 114L,
    n_studies = 2L,
    n_concentrations = 1989L,
    age_range = "18-75 years (median 57)",
    weight_range = "34.6-102.0 kg (median 61.4)",
    sex_female_pct = 37.7,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Adult primary orthotopic liver transplant recipients during the first",
      "6 months after transplantation (median post-operative day 43, range",
      "0-193); multi-organ transplantation and severe renal impairment",
      "excluded."
    ),
    dose_range = paste(
      "Oral immediate-release tacrolimus (Prograf) twice daily, starting at",
      "0.02-0.03 mg/kg per dose and titrated to whole-blood troughs of 8-12",
      "ng/mL (months 0-3) and 5-8 ng/mL thereafter; daily dose median 4 mg",
      "(range 1-13)."
    ),
    regions = "Thailand (King Chulalongkorn Memorial Hospital, Bangkok)",
    sampling_design = paste(
      "Pooled retrospective TDM cohort (97 patients, troughs, 2014-2023) and",
      "prospective intensive cohort (17 patients, steady-state samples at 0,",
      "0.5, 1, 2, 4, 6 and 8 h, 2023-2024); 93.5% of observations are",
      "troughs."
    ),
    cyp3a5_genotype = c(
      `Ex/Ex (recipient/donor)` = 27.2,
      `Ex/Nex` = 30.7,
      `Nex/Ex` = 23.7,
      `Nex/Nex` = 18.4
    ),
    notes = paste(
      "Table 1: hemoglobin median 10.9 g/dL (5.7-16.7); trough median 6.3",
      "ng/mL (1.2-31.3); CYP3A5*3 allele frequency 65.55%. Co-medication by",
      "concentration record: mycophenolate 87%, prednisolone 88%, fluconazole",
      "59%, proton pump inhibitor 61%."
    )
  )

  ini({
    # ---- Structural parameters (Komenkul 2026 Table 3, final model) -----
    # Reference subject: recipient and donor CYP3A5 expressors, no
    # fluconazole, no prednisolone, daily dose 4 mg/day, hemoglobin 11 g/dL.
    lka <- log(0.95); label("First-order absorption rate constant (1/h)")  # Table 3 'Ka (h-1)' = 0.95 (RSE 10.32%)
    lcl <- log(23.61); label("Apparent oral clearance CL/F for the reference subject (L/h)")  # Table 3 'CL/F (L/h)' = 23.61 (RSE 4.68%); final-model equation intercept 23.61
    lvc <- log(81.88); label("Apparent central volume Vc/F (L)")  # Table 3 'Vc/F (L)' = 81.88 (RSE 7.45%)
    lq <- log(36.99); label("Apparent intercompartmental clearance Q/F (L/h)")  # Table 3 'Q/F (L/h)' = 36.99 (RSE 13.22%)
    lvp <- log(367.71); label("Apparent peripheral volume Vp/F (L)")  # Table 3 'Vp/F (L)' = 367.71 (RSE 12.07%)

    # ---- Covariate effects on CL/F (Komenkul 2026 final-model equation) --
    # CL/F = 23.61 x (DD/4)^0.22 x (HB/11)^-0.15 x (1 + (-0.33 x FLU))
    #        x (1 + 0.08 x PSL) x COMB
    e_dose_tac_cl <- 0.22; label("Power exponent of (DOSE_TAC_MGD / 4 mg/day) on CL/F (unitless)")  # Table 3 'Daily dose' = 0.22 (RSE 16.49%)
    e_hgb_cl <- -0.15; label("Power exponent of (HGB / 11 g/dL) on CL/F (unitless)")  # Table 3 'Hemoglobin' = -0.15 (RSE 25.74%)
    e_conmed_fluconazole_cl <- -0.33; label("Fractional change in CL/F with fluconazole co-administration (unitless)")  # Table 3 'Fluconazole' = -0.33 (RSE 4.44%)
    e_conmed_prednisolone_cl <- 0.08; label("Fractional change in CL/F with prednisolone co-administration (unitless)")  # Table 3 'Prednisolone' = 0.08 (RSE 22.96%)
    e_cyp3a5_exnex_cl <- -0.15; label("Fractional change in CL/F, recipient expressor + donor non-expressor vs Ex/Ex (unitless)")  # Table 3 'Expressors/non-expressors' = -0.15 (RSE 23.01%); equation COMB = 0.85
    e_cyp3a5_nexex_cl <- -0.16; label("Fractional change in CL/F, recipient non-expressor + donor expressor vs Ex/Ex (unitless)")  # Table 3 'Non-expressors/expressors' = -0.16 (RSE 20.56%); equation COMB = 0.84
    e_cyp3a5_nexnex_cl <- -0.24; label("Fractional change in CL/F, recipient and donor non-expressors vs Ex/Ex (unitless)")  # Table 3 'Non-expressors/non-expressors' = -0.24 (RSE 15.68%); equation COMB = 0.76

    # ---- Inter-individual variability (Komenkul 2026 Table 3) -----------
    # Phoenix NLME omega estimates (variances of exponential etas); IIV on
    # Ka and Q/F was not estimated (Table 3 'Fixed'; the authors' code has
    # no eta on Ka or Cl2). Correlation CL-Vc 94.54% gives the covariance
    # 0.9454 x sqrt(0.07 x 0.24) = 0.12254.
    etalcl + etalvc ~ c(0.07, 0.12254, 0.24) # Table 3 'IIV CL/F' = 0.07, 'IIV Vc/F' = 0.24, 'Corr CL-Vc, %' = 94.54
    etalvp ~ 1.22 # Table 3 'IIV Vp/F' = 1.22 (RSE 20.88%)

    # ---- Residual error (Komenkul 2026 Table 3) -------------------------
    propSd <- 0.2776; label("Proportional residual error (fraction)")  # Table 3 'Proportional residual error (%)' = 27.76 (RSE 2.72%)
  })

  model({
    # ---- 1. Combined recipient/donor CYP3A5 genotype (COMB) -------------
    # Reference Ex/Ex (both expressors) has multiplier 1; exactly one of the
    # three indicators below is 1 for any other combination.
    ind_exnex <- CYP3A5_EXPR * (1 - CYP3A5_EXPR_DONOR)
    ind_nexex <- (1 - CYP3A5_EXPR) * CYP3A5_EXPR_DONOR
    ind_nexnex <- (1 - CYP3A5_EXPR) * (1 - CYP3A5_EXPR_DONOR)
    comb_factor <- (1 + e_cyp3a5_exnex_cl * ind_exnex) *
      (1 + e_cyp3a5_nexex_cl * ind_nexex) *
      (1 + e_cyp3a5_nexnex_cl * ind_nexnex)

    # ---- 2. Individual parameters -----------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (DOSE_TAC_MGD / 4)^e_dose_tac_cl *
      (HGB / 11)^e_hgb_cl *
      (1 + e_conmed_fluconazole_cl * CONMED_FLUCONAZOLE) *
      (1 + e_conmed_prednisolone_cl * CONMED_PREDNISOLONE) *
      comb_factor
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 3. ODE system ----------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 4. Observation and error -----------------------------------------
    # central in mg, vc in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
