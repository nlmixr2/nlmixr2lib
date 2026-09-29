Chen_2021_tacrolimus <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption and",
    "elimination for oral tacrolimus in Chinese adult kidney transplant",
    "recipients during the first 90 postoperative days, pooled from two",
    "centres (Chen 2021 Pharmgenomics Pers Med, NONMEM 7.4). Apparent oral",
    "clearance CL/F carries median-normalised power effects of Cockcroft-Gault",
    "creatinine clearance, haematocrit (inverse ratio) and the daily tacrolimus",
    "dose, a 1.29-fold CYP3A5*1-carrier (expresser) multiplier, and a",
    "four-level multiplier for the 48-hour cumulative Wuzhi capsule dose",
    "(0 mg reference; below 45 mg 0.566; exactly 45 mg 0.783; above 45 mg",
    "0.598). Absorption rate constant fixed at 3.09 1/h from a previous",
    "analysis because only trough concentrations were available; no",
    "covariate on Vd/F. Exponential IIV on CL/F and Vd/F; exponential",
    "(log-normal) residual error estimated separately for each centre."
  )
  reference <- paste(
    "Chen L, Yang Y, Wang X, Wang C, Lin W, Jiao Z, Wang Z (2021).",
    "Wuzhi Capsule Dosage Affects Tacrolimus Elimination in Adult Kidney",
    "Transplant Recipients, as Determined by a Population Pharmacokinetics",
    "Analysis. Pharmgenomics Pers Med 14:1093-1106.",
    "doi:10.2147/PGPM.S321997."
  )
  vignette <- "Chen_2021_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Chen 2021 Methods 'Sample Collection and Bioanalytical Assay': whole
  # blood trough concentrations by CMIA (Changhai) and EMIT (Huashan); the
  # Huashan values were converted to CMIA equivalents before modelling.
  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  # The paper estimated one exponential residual error per centre (Table 2
  # rows 'CH exponential error' and 'HS exponential error'); the single
  # canonical expSd cannot hold both, so each carries a centre suffix.
  paper_specific_residual_sds <- c("expSd_changhai", "expSd_huashan")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation; raw mL/min, NOT BSA-normalised",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Chen 2021 Results: 'CCR is the creatinine clearance rate by",
        "Cockcroft-Gault equation'. Enters CL/F as the median-normalised power",
        "(CRCL / 45.5)^0.179; 45.5 mL/min is the pooled Table 1 median (range",
        "4.9-123.9). Raw Cockcroft-Gault scale, following the raw-CrCl",
        "precedent documented under the canonical CRCL register entry. The",
        "exponent is POSITIVE, i.e. CL/F rises with renal function; the",
        "Discussion's sentence 'patients with a lower creatinine clearance rate",
        "had a higher CL/F' contradicts the printed equation and Table 2, which",
        "are followed here. Time-varying within subject in the source data",
        "(post-transplant graft recovery)."
      ),
      source_name = "CCR"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "% (volume fraction times 100)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as the INVERSE median-normalised power (28.9 / HCT)^0.503",
        "exactly as printed in the Results equation; 28.9 % is the pooled",
        "Table 1 median (range 17.6-48.2). Equivalent to (HCT / 28.9)^-0.503:",
        "lower haematocrit gives higher CL/F, which the Discussion attributes",
        "to tacrolimus's ~95 % erythrocyte binding (a lower haematocrit leaves",
        "more unbound drug in plasma for hepatic clearance). Percent scale in",
        "the source, matching the canonical column. Time-varying within",
        "subject (early post-operative anaemia recovers over the 90 days)."
      ),
      source_name = "HCT"
    ),
    DOSE_TAC_MGD = list(
      description = "Total daily tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Chen 2021 'DOSE is the daily dose of tacrolimus'. Enters CL/F as the",
        "median-normalised power (DOSE_TAC_MGD / 5)^0.351; 5 mg/day is the",
        "pooled Table 1 median (range 1-11). Time-varying: updated whenever",
        "the dose is titrated to the trough target. For a q12h regimen",
        "DOSE_TAC_MGD is twice the per-administration amt; keep the two",
        "consistent when simulating. See the register entry for why a",
        "dose-on-clearance effect in a trough-titrated TDM dataset is partly a",
        "surrogate for unmeasured fast-clearance characteristics."
      ),
      source_name = "DOSE"
    ),
    CYP3A5_EXPR = list(
      description = "CYP3A5 expresser: 1 = at least one CYP3A5*1 allele (*1/*1 or *1/*3), 0 = CYP3A5*3/*3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CYP3A5*3/*3 non-expresser)",
      notes = paste(
        "Time-fixed germline genotype at rs776746. Chen 2021 multiplies CL/F",
        "by 1.29 'if CYP3A5*1*3 and *1*1 carriers'. Table 1 pooled genotype",
        "counts (AA/AG/GG) 12/48/82, so 60 of 142 recipients are expressers."
      ),
      source_name = "CYP3A5*3"
    ),
    DOSE_WUZHI_MG48H = list(
      description = "Cumulative Wuzhi capsule dose over the 48 h preceding the observation",
      units = "mg per 48 h",
      type = "continuous",
      reference_category = "0 (no Wuzhi capsule)",
      notes = paste(
        "Chen 2021 Methods 'Covariate Analysis': 'co-therapy with Wuzhi",
        "capsule at cumulative doses for 48 h'. The model does NOT use the",
        "dose continuously: it bins it into four categories with separate",
        "CL/F multipliers -- 0 mg (reference, 1), below 45 mg (0.566),",
        "exactly 45 mg (0.783) and above 45 mg (0.598). The Changhai regimens",
        "(11.25 mg per capsule strength) map to the bins as: 11.25 mg once",
        "daily = 22.5 mg / 48 h (below 45); 11.25 mg twice daily = 45 mg",
        "(exactly 45); 11.25 mg three times daily = 67.5 mg, 22.5 mg twice",
        "daily = 90 mg and 22.5 mg three times daily = 135 mg (above 45).",
        "Table 1 Changhai median 45 mg (range 0-135). Huashan patients",
        "received no Wuzhi capsule (value 0). The '= 45' bin is matched with a",
        "0.01 mg tolerance in model() so that a floating-point 45 lands in it."
      ),
      source_name = "WZ"
    ),
    STUDY_HUASHAN = list(
      description = "Recruiting-centre indicator: 1 = Huashan Hospital cohort, 0 = Changhai Hospital cohort",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Changhai Hospital)",
      notes = paste(
        "Selects the centre-specific exponential residual error (Table 2",
        "'CH exponential error' 0.0606 vs 'HS exponential error' 0.0887). It",
        "carries no structural effect. Huashan concentrations were measured",
        "by EMIT and converted to CMIA equivalents (CMIA = 0.93 * EMIT +",
        "0.36) before modelling, so model predictions are CMIA-scale for both",
        "centres."
      ),
      source_name = "CH / HS"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods 'Covariate Analysis') but not retained. Table 1 pooled median 60.1 kg (range 36-86.7)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened but not retained. Table 1 pooled median 169 cm (range 150-182)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened but not retained. Table 1 pooled median 15 U/L (range 6-187)."
    ),
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened (Methods 'Covariate Analysis') but not retained; not summarised in Table 1."
    ),
    POD = list(
      description = "Postoperative day",
      units = "days",
      type = "continuous",
      notes = "Screened but not retained. Table 1 pooled median 16 days (range 2-90)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 142L,
    n_studies = 2L,
    n_concentrations = 1378L,
    age_range = "20-67 years (Table 1 pooled median 40.0)",
    weight_range = "36-86.7 kg (Table 1 pooled median 60.1)",
    sex_female_pct = 32.4,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Chinese adult kidney transplant recipients followed to post-operative",
      "day 90 on triple immunosuppression (tacrolimus + mycophenolic acid +",
      "corticosteroid). Excluded were patients with missing required",
      "information, dialysis, acute rejection, second transplantation or a",
      "switch of immunosuppressant."
    ),
    dose_range = paste(
      "Oral tacrolimus twice daily, titrated to trough targets; daily dose",
      "median 5 mg (range 1-11). Changhai initial 0.05-1.0 mg/kg/day (trough",
      "10-15 ng/mL month 1, 8-10 ng/mL months 1-3); Huashan initial",
      "0.1-0.15 mg/kg/day (8-10 then 6-8 ng/mL). Wuzhi capsule co-medication",
      "at Changhai only (709 of 758 Changhai records): 11.25 mg once, twice",
      "or three times daily, or 22.5 mg twice or three times daily."
    ),
    regions = "China (Shanghai: Changhai Hospital and Huashan Hospital)",
    sampling_design = paste(
      "1378 whole-blood trough concentrations (C0), 758 from 90 Changhai",
      "patients (March-September 2016, CMIA, LLOQ 1.5 ng/mL) and 620 from 52",
      "Huashan patients (May 2009-December 2013, EMIT converted to CMIA",
      "equivalents), average 37 per subject."
    ),
    cyp3a5_genotype = c(`*1/*1 (AA)` = 8.45, `*1/*3 (AG)` = 33.80, `*3/*3 (GG)` = 57.75),
    notes = paste(
      "Table 1 (pooled): 96 male / 46 female; haematocrit median 28.9 %",
      "(17.6-48.2); Cockcroft-Gault creatinine clearance median 45.5 mL/min",
      "(4.9-123.9); observed trough median 9.9 ng/mL (2.3-42.6). NONMEM 7.4;",
      "stepwise covariate modelling with forward dOFV 10.83 and backward",
      "dOFV 10.83 (P < 0.001). Evaluated by goodness-of-fit plots, a",
      "bootstrap (95.6 % of 1000 runs successful), a pc-VPC and NPDE."
    )
  )

  ini({
    # ---- Structural parameters (Chen 2021 Table 2, final model) ---------
    # Ka fixed: Methods 'Base Model Development': 'The Ka was fixed at 3.09
    # h-1 with no BSV based on the findings of a previous study' (Zuo 2013,
    # Pharmacogenet Genomics, same population); Table 2 footnote '*Ka was
    # fixed to the published value.'
    lka <- fixed(log(3.09)); label("First-order absorption rate constant (1/h), from Zuo 2013")  # Table 2 'Ka* (h-1)' = 3.09 (fixed)

    # Reference subject: CYP3A5*3/*3, no Wuzhi capsule, CRCL 45.5 mL/min,
    # HCT 28.9 %, daily tacrolimus dose 5 mg (Discussion: 'an adult patient
    # with CYP3A5*3/*3 would require a daily tacrolimus dose of 5 mg without
    # co-treatment with Wuzhi capsule with a median creatinine clearance rate
    # of 45.5 mL/min and hematocrit of 29%; the typical CL/F was estimated to
    # be 14.4 L/h').
    lcl <- log(14.4); label("Apparent oral clearance CL/F for the reference subject (L/h)")  # Table 2 'CL/F (L/h)' final = 14.4 (RSE 7%); bootstrap median 14.30 (11.84-16.38)
    lvc <- log(275); label("Apparent volume of distribution Vd/F (L)")  # Table 2 'Vd/F (L)' final = 275 (RSE 29%); Results 'Vd/F = 275'

    # ---- Covariate effects on CL/F (Chen 2021 Results equation) --------
    # CL/F = 14.4 x (CCR/45.5)^0.179 x (28.9/HCT)^0.503 x (DOSE/5)^0.351
    #        x 1.29 (CYP3A5*1 carriers) x Wuzhi multiplier
    e_crcl_cl <- 0.179; label("Power exponent of (CRCL / 45.5 mL/min) on CL/F (unitless)")  # Table 2 'Creatinine clearance rate' = 0.179 (RSE 18%); bootstrap 0.176 (0.097-0.236)
    e_hct_cl <- 0.503; label("Power exponent of the inverse ratio (28.9 % / HCT) on CL/F (unitless)")  # Table 2 'Haematocrit' = 0.503 (RSE 14%); bootstrap 0.508 (0.385-0.657)
    e_dose_tac_cl <- 0.351; label("Power exponent of (DOSE_TAC_MGD / 5 mg/day) on CL/F (unitless)")  # Table 2 'DOSE' = 0.351 (RSE 13%); bootstrap 0.343 (0.250-0.432)
    e_cyp3a5_expr_cl <- 1.29; label("CL/F multiplier for CYP3A5 expressers (*1/*1 or *1/*3) (unitless)")  # Table 2 'CYP3A5*1/*1 and *1/*3' = 1.29 (RSE 5%); bootstrap 1.288 (1.160-1.424)

    # Wuzhi capsule, 48-h cumulative dose, four bins; 0 mg is the reference
    # with multiplier 1 (Table 2 'WZ=0mg' = 1).
    e_wuzhi_lt45_cl <- 0.566; label("CL/F multiplier for 48-h Wuzhi capsule dose below 45 mg (unitless)")  # Table 2 'WZ<45mg' = 0.566 (RSE 27%); bootstrap 0.575 (0.366-0.810)
    e_wuzhi_eq45_cl <- 0.783; label("CL/F multiplier for 48-h Wuzhi capsule dose of exactly 45 mg (unitless)")  # Table 2 'WZ=45mg' = 0.783 (RSE 8%); bootstrap 0.778 (0.671-0.921)
    e_wuzhi_gt45_cl <- 0.598; label("CL/F multiplier for 48-h Wuzhi capsule dose above 45 mg (unitless)")  # Table 2 'WZ>45mg' = 0.598 (RSE 6%); bootstrap 0.6 (0.527-0.691)

    # ---- Inter-individual variability (Chen 2021 Table 2) ---------------
    # Exponential BSV (Methods Eq. 1). Table 2 prints BSV as a percent
    # (CL/F 25.4 %, V/F 51.5 %); converted with omega^2 = log(CV^2 + 1).
    # The alternative reading omega = CV/100 would give 0.0645 and 0.265.
    etalcl ~ 0.0625 # Table 2 BSV 'CL/F(%)' final = 25.4 (RSE 15%) [shrinkage 12%]; log(0.254^2 + 1) = 0.0625
    etalvc ~ 0.2354 # Table 2 BSV 'V/F(%)' final = 51.5 (RSE 35%) [shrinkage 25%]; log(0.515^2 + 1) = 0.2354

    # ---- Residual error (Chen 2021 Table 2) ------------------------------
    # 'The residual error model was selected using an exponential method',
    # one per centre. The printed values are NONMEM SIGMA variances: the
    # estimate +/- 1.96 x RSE x estimate reproduces the bootstrap 95 % CI
    # (CH 0.0440-0.0772 vs 0.048-0.075; HS 0.0748-0.1026 vs 0.075-0.102),
    # so each log-scale SD is sqrt(variance).
    expSd_changhai <- 0.2462; label("Exponential residual SD, Changhai Hospital (log scale)")  # Table 2 'CH exponential error' = 0.0606 (RSE 14%) [shrinkage 8%]; sqrt(0.0606) = 0.2462
    expSd_huashan <- 0.2978; label("Exponential residual SD, Huashan Hospital (log scale)")  # Table 2 'HS exponential error' = 0.0887 (RSE 8%) [shrinkage 4%]; sqrt(0.0887) = 0.2978
  })

  model({
    # ---- 1. Covariate factors on CL/F (Chen 2021 Results equation) ------
    crcl_factor <- (CRCL / 45.5)^e_crcl_cl
    hct_factor <- (28.9 / HCT)^e_hct_cl
    dose_factor <- (DOSE_TAC_MGD / 5)^e_dose_tac_cl
    cyp3a5_factor <- 1 + (e_cyp3a5_expr_cl - 1) * CYP3A5_EXPR

    # Wuzhi capsule 48-h cumulative dose bins: 0 (reference), < 45, = 45,
    # > 45 mg. The '= 45' bin uses a 0.01 mg tolerance.
    wz_eq45 <- (abs(DOSE_WUZHI_MG48H - 45) < 0.01)
    wz_lt45 <- (DOSE_WUZHI_MG48H > 0) * (DOSE_WUZHI_MG48H < 45) * (1 - wz_eq45)
    wz_gt45 <- (DOSE_WUZHI_MG48H > 45) * (1 - wz_eq45)
    wuzhi_factor <- 1 + (e_wuzhi_lt45_cl - 1) * wz_lt45 +
      (e_wuzhi_eq45_cl - 1) * wz_eq45 +
      (e_wuzhi_gt45_cl - 1) * wz_gt45

    # ---- 2. Individual parameters ----------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * crcl_factor * hct_factor * dose_factor *
      cyp3a5_factor * wuzhi_factor
    vc <- exp(lvc + etalvc)

    # ---- 3. ODE system (one compartment, first-order absorption) ---------
    kel <- cl / vc
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 4. Observation and error ----------------------------------------
    # central in mg, vc in L -> mg/L; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    expSd_i <- expSd_changhai * (1 - STUDY_HUASHAN) + expSd_huashan * STUDY_HUASHAN
    Cc ~ lnorm(expSd_i)
  })
}
