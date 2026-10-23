Kim_2026_lenalidomide <- function() {
  description <- "Two-compartment oral population PK model with first-order absorption (Ka fixed), first-order elimination and linear (concentration-independent) plasma protein binding for lenalidomide in Korean adults with multiple myeloma spanning normal renal function to severe renal impairment, fitted simultaneously to total and unbound plasma concentrations (Kim 2026). The structural parameters are apparent UNBOUND parameters: CL/F depends on CKD-EPI 2021 eGFR (exponential, centred at 46.92 mL/min/1.73 m^2) and on body surface area (linear, centred at 1.64 m^2); V1/F depends on body surface area (exponential). The total concentration is the unbound concentration divided by a constant unbound fraction of 0.614."
  reference <- paste(
    "Kim H, Park S-S, Choi S, Han S, Han S, Lee Y, Lee JH, Kim D, Kim H, Yang S, Shugg T, Chung EK.",
    "Population Pharmacokinetics and Model-Informed Precision Dosing of Lenalidomide Incorporating",
    "Total and Unbound Plasma Concentrations in Renally Impaired Patients with Multiple Myeloma.",
    "Drug Des Devel Ther. 2026;20. doi:10.2147/DDDT.S605741"
  )
  vignette <- "Kim_2026_lenalidomide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "lenalidomide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lenalidomide unbound", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "lenalidomide unbound", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Estimated glomerular filtration rate by the CKD-EPI 2021 (creatinine) equation, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject (baseline). Enters unbound CL/F as the exponential",
        "effect exp(e_crcl_cl * (CRCL - 46.92)) (Kim 2026 final-model equation, Results",
        "'Covariate Analysis'). The 46.92 reference is printed in the equation; the",
        "cohort median in Table 1 is 47 (range 11.1-94.7). Cockcroft-Gault creatinine",
        "clearance (mL/min) was tested as an alternative renal marker and gave a smaller",
        "OFV drop (-38.11 versus -48.72 for eGFR), so it was not used; do not supply a",
        "Cockcroft-Gault value to this model."
      ),
      source_name = "eGFR"
    ),
    BSA = list(
      description = "Body surface area by the Du Bois method",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject (baseline). Enters unbound CL/F as the LINEAR factor",
        "(1 + e_bsa_cl * (BSA - 1.64)) and unbound V1/F as the exponential factor",
        "exp(e_bsa_vc * (BSA - 1.64)). 1.64 m^2 is the cohort median (Table 1, range",
        "1.41-2.01). The linear factor reaches zero at BSA = 0.73 m^2, far below the",
        "studied range. See the vignette for why the CL/F term is linear rather than",
        "the additive offset the printed equation shows."
      ),
      source_name = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened in stepwise covariate modelling; not retained. Cohort 12 of 30 (40%) male (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 68 years (range 47-84) (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened as a body-size descriptor; BSA was retained instead. Cohort median 62 kg (range 47-86) (Table 1)."
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 171 (range 39-381) (Table 1)."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 5.4 (range 2.8-8.4) (Table 1)."
    ),
    NEUT = list(
      description = "Absolute neutrophil count",
      units = "10^9/L",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 3.0 (range 1.4-6.1) (Table 1)."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/dL",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 6.9 g/dL (range 5.3-9.3) (Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened; not retained, including on the unbound fraction. Cohort median 4.2 g/dL (range 3.0-5.2); no patient was hypoalbuminaemic (Table 1)."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened; not retained. Cohort median 0.2 mg/dL (range 0.01-11.7) (Table 1; the Methods state mg/L, Table 1 mg/dL)."
    ),
    SNP_ABCB1_RS1128503 = list(
      description = "ABCB1 1236C>T (rs1128503) genotype",
      units = "(genotype)",
      type = "categorical",
      notes = "Tested under dominant, recessive and additive coding; did not improve the model fit as a covariate on any unbound PK parameter (Results 'Covariate Analysis')."
    ),
    SNP_ABCB1_RS2032582 = list(
      description = "ABCB1 2677G>T/A (rs2032582) genotype",
      units = "(genotype)",
      type = "categorical",
      notes = "Tested under dominant, recessive and additive coding; not retained (Results 'Covariate Analysis')."
    ),
    SNP_ABCB1_RS1045642 = list(
      description = "ABCB1 3435C>T (rs1045642) genotype",
      units = "(genotype)",
      type = "categorical",
      notes = "Tested under dominant, recessive and additive coding; not retained (Results 'Covariate Analysis'). Compare Liang_2024_lenalidomide, where this SNP was retained on V/F."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    n_observations = "345 plasma concentrations (171 total, 174 unbound)",
    age_range = "47-84 years",
    age_median = "68 years",
    weight_range = "47-86 kg",
    weight_median = "62 kg",
    bsa_range = "1.41-2.01 m^2",
    bsa_median = "1.64 m^2",
    sex_female_pct = 60,
    race_ethnicity = c(Korean = 100),
    disease_state = "Multiple myeloma starting lenalidomide-based maintenance therapy, lenalidomide-naive.",
    renal_function = "Stratified by CKD-EPI 2021 eGFR: 11 severe renal impairment (< 30 mL/min/1.73 m^2), 11 moderate (30 to < 60), 8 control (>= 60). Overall eGFR median 47 (range 11.1-94.7). No renal replacement therapy.",
    hepatic_function = "Serum albumin median 4.2 g/dL (range 3.0-5.2); no clinically significant hypoalbuminaemia.",
    dose_range = "5, 7.5 or 10 mg orally once daily on days 1-28 of each cycle, adjusted for renal function. Serial sampling at pre-dose and 2, 4, 8, 12 and 24 h after the dose on day 3 of cycle 1 (treated as steady state), fasting 10 h before and 4 h after the dose.",
    regions = "Republic of Korea (Seoul St. Mary's Hospital, single centre).",
    notes = "Baseline demographics from Kim 2026 Table 1. NONMEM 7.5.1, Laplacian estimation with interaction. Total and unbound concentrations measured by validated HPLC-MS/MS (LLOQ 5 ng/mL for both) after centrifugal ultrafiltration."
  )

  ini({
    # Structural parameters: Kim 2026 Table 3, 'NONMEM Final Model' column.
    # CL/F, V1/F, Q/F and V2/F are apparent UNBOUND parameters (Table 3
    # abbreviations: 'apparent unbound clearance', 'apparent central unbound
    # volume of distribution', ...).
    lka <- fixed(log(1.16)); label("First-order absorption rate constant Ka (1/h)") # Table 3 theta-Ka = 1.16 (N/E); Results text: 'Ka was fixed at 1.16 h-1 ... based on our NCA analysis'
    lcl <- log(8.52); label("Apparent unbound clearance CL/F at eGFR 46.92 and BSA 1.64 m^2 (L/h)") # Table 3 theta-CL/F = 8.52 (RSE 3.9%)
    lvc <- log(48.8); label("Apparent central unbound volume V1/F at BSA 1.64 m^2 (L)") # Table 3 theta-V1/F = 48.8 (RSE 10.7%)
    lq <- log(8.05); label("Apparent unbound intercompartmental clearance Q/F (L/h)") # Table 3 theta-Q/F = 8.05 (RSE 27.9%)
    lvp <- log(21.2); label("Apparent peripheral unbound volume V2/F (L)") # Table 3 theta-V2/F = 21.2 (RSE 13.8%); the page-10 equation prints 21.1, the bootstrap median
    lfu <- log(0.614); label("Unbound fraction of lenalidomide in plasma (fraction)") # Table 3 theta-Fu = 61.4 percent (RSE 1.32%); page-10 equation Fu = 61.4 percent

    # Covariate effects: Table 3 and the page-10 final-model equations.
    e_crcl_cl <- 0.0131; label("Exponential effect of eGFR on CL/F (per mL/min/1.73 m^2)") # Table 3 theta-CLeGFR = 0.0131 (RSE 15.1%)
    e_bsa_cl <- 1.1; label("Linear effect of BSA on CL/F (per m^2)") # Table 3 theta-CLBSA = 1.1 (RSE 29.6%)
    e_bsa_vc <- 1.16; label("Exponential effect of BSA on V1/F (per m^2)") # Table 3 theta-V1BSA = 1.16 (RSE 26.2%)

    # IIV: Table 3 reports omega as a percentage under an 'omega' (not
    # 'omega^2') header; read as the log-scale SD, variance = SD^2.
    etalcl ~ 0.035721 # Table 3 omega CL/F = 18.9 percent (RSE 13.3%, shrinkage 2%) -> 0.189^2
    etalvc ~ 0.094864 # Table 3 omega V1/F = 30.8 percent (RSE 24.7%, shrinkage 8%) -> 0.308^2

    # Residual error: combined additive + proportional on each output (Table 3).
    propSd <- 0.126; label("Proportional residual error, total concentration (fraction)") # Table 3 sigma total-proportional = 12.6 percent (RSE 16.3%)
    addSd <- 2.19; label("Additive residual error, total concentration (ng/mL)") # Table 3 sigma total-additive = 2.19 ng/mL (RSE 25.2%)
    propSd_Cu <- 0.125; label("Proportional residual error, unbound concentration (fraction)") # Table 3 sigma unbound-proportional = 12.5 percent (RSE 16.8%)
    addSd_Cu <- 1.11; label("Additive residual error, unbound concentration (ng/mL)") # Table 3 sigma unbound-additive = 1.11 ng/mL (RSE 25.2%)
  })

  model({
    # Page-10 final-model equations. The BSA effect on CL/F is encoded as the
    # linear factor (1 + theta * (BSA - 1.64)); the printed equation places
    # '+ 1.1 * (BSA - 1.64)' outside the exponential as an additive offset,
    # which contradicts the paper's own Figure 6 dosing heatmaps (see vignette).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + e_crcl_cl * (CRCL - 46.92)) * (1 + e_bsa_cl * (BSA - 1.64))
    vc <- exp(lvc + etalvc + e_bsa_vc * (BSA - 1.64))
    q <- exp(lq)
    vp <- exp(lvp)
    fu <- exp(lfu)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Central and peripheral hold the unbound-equivalent amount (mg): the
    # volumes are unbound volumes, so amount / volume is the unbound concentration.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg/L -> ng/mL; linear protein binding C_total = C_unbound / Fu (page 5).
    Cu <- central / vc * 1000
    Cc <- Cu / fu

    Cc ~ add(addSd) + prop(propSd)
    Cu ~ add(addSd_Cu) + prop(propSd_Cu)
  })
}
