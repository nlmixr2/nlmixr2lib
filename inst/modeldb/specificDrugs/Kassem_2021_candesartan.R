Kassem_2021_candesartan <- function() {
  description <- paste(
    "One-compartment population PK model for oral candesartan (dosed as the",
    "prodrug candesartan cilexetil) in white adults with chronic heart failure",
    "and reduced ejection fraction (Kassem 2021). First-order absorption with",
    "a lag time feeds a one-compartment disposition model. Apparent clearance",
    "carries median-normalised power effects of body weight and eGFR and a",
    "multiplicative effect of diabetes; interindividual variability is",
    "estimated on apparent clearance only. Combined additive + proportional",
    "residual error.",
    sep = " "
  )
  reference <- paste(
    "Kassem I, Sanche S, Li J, Bonnefois G, Dube MP, Rouleau JL, Tardif JC,",
    "White M, Turgeon J, Nekka F, de Denus S. Population pharmacokinetics of",
    "candesartan in patients with chronic heart failure. Clin Transl Sci.",
    "2021 Jan;14(1):194-203. doi:10.1111/cts.12842. PMCID: PMC7877833.",
    sep = " "
  )
  vignette <- "Kassem_2021_candesartan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "candesartan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "candesartan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as the median-normalised power function of the Results",
        "final-model equation, (Weight/82.45)^0.963. 82.45 kg is the cohort",
        "median weight (Table 2 footnote: 'for a median weight of 82.45 kg');",
        "the cohort mean (SD) is 84.0 (19.1) kg (Table 1).",
        sep = " "
      ),
      source_name = "Weight"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL/F as the median-normalised power function of the Results",
        "final-model equation, (eGFR/74)^0.56. 74 mL/min/1.73 m^2 is the",
        "cohort median (Table 2 footnote); the cohort mean (SD) is 74.2",
        "(22.3) (Table 1). The estimating equation (MDRD vs CKD-EPI) is not",
        "stated in Kassem 2021; the unit is BSA-normalised per Table 1.",
        sep = " "
      ),
      source_name = "eGFR"
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus comorbidity indicator (1 = diabetic, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Enters CL/F multiplicatively as 0.682^Diabetes (Results final-model",
        "equation; Supplementary Table S1 footnote: 'Diabetes = 0 (No), 1",
        "(Yes)'), i.e. a 31.8% lower CL/F in patients with diabetes. Type 1",
        "vs type 2 is not distinguished. Cohort prevalence 32.7% (Table 1).",
        sep = " "
      ),
      source_name = "Diabetes"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Passed the preliminary univariate screen (Results) but failed the",
        "forward step: Supplementary Table S1 univariate model",
        "CL/F = theta1 * (Age/66.26)^theta7, delta OF -0.242 (not",
        "significant). Cohort mean (SD) 65.6 (10.0) years (Table 1).",
        sep = " "
      ),
      source_name = "Age"
    ),
    NTPROBNP = list(
      description = "N-terminal pro-B-type natriuretic peptide",
      units = "ng/L",
      type = "continuous",
      notes = paste(
        "Passed the preliminary univariate screen but failed the forward",
        "step: Supplementary Table S1 univariate model",
        "CL/F = theta1 * (NT_proBNP/726)^theta8, delta OF +1.132. Reported",
        "in ng/L (= pg/mL; 726 ng/L = 0.726 ng/mL on the register's ng/mL",
        "scale). Cohort mean (SD) 1,291.9 (1,818.1) ng/L (Table 1).",
        sep = " "
      ),
      source_name = "NT_proBNP"
    ),
    CONMED_FUROSEMIDE = list(
      description = "Concomitant furosemide indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Retained in the forward step (Supplementary Table S1, delta OF",
        "-5.775) but removed in the backward step (P > 0.001; Results).",
        "71.9% of the cohort received furosemide (Table 1).",
        sep = " "
      ),
      source_name = "Furosemide"
    ),
    SEXF = list(
      description = "Female-sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Retained in the forward step (Supplementary Table S1, 'Sex = 0",
        "(Male), 1 (Female)', delta OF -12.887) but removed in the backward",
        "step (Results). The structural model gave CL/F 7.96 L/h in men vs",
        "5.9 L/h in women; the paper attributes the difference to the lower",
        "weight and eGFR of women (Discussion). 17% of the cohort was",
        "female (Table 1).",
        sep = " "
      ),
      source_name = "Sex"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 281,
    n_studies = 1,
    n_observations = 1455,
    age_mean = "65.6 years (SD 10.0)",
    weight_mean = "84.0 kg (SD 19.1); median 82.45 kg",
    sex_female_pct = 17,
    race_ethnicity = "White (100%; analysis restricted to white patients)",
    disease_state = paste(
      "Symptomatic chronic heart failure with left ventricular ejection",
      "fraction <= 40% and NYHA class II-IV (78.3% class II); 32.7% with",
      "diabetes mellitus, 56.2% hypertension, 27.0% atrial fibrillation",
      sep = " "
    ),
    renal_function = paste(
      "eGFR mean 74.2 (SD 22.3) mL/min/1.73 m^2, median 74; 20.6% normal,",
      "48.4% mild, 20.6% mild-to-moderate, 10.3% moderate-to-severe renal",
      "dysfunction (Table 1)",
      sep = " "
    ),
    dose_range = "4-32 mg candesartan cilexetil orally once daily (titrated 4, 8, 16, 32 mg)",
    regions = "Canada (16 centres)",
    notes = paste(
      "Population PK sub-study of a prospective, open-label, nonrandomized",
      "pharmacogenomic study (de Denus 2018, Pharmacogenomics 19:599-612) with",
      "visits at weeks 0, 2, 4, 6, 8 and 16. First sample 2 h after the first",
      "4 mg dose; later samples at variable times after the previous dose",
      "(majority 0-4 h post-dose; Figure S2). Only subjects with at least one",
      "concentration within 30 h of a documented dose were included. Assay",
      "LC-MS, 1.00-250 ng/mL. Baseline demographics in Table 1.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- Kassem 2021 Table 2 (final model); bootstrap
    # means/CIs in Table 3 bracket every value. The typical CL/F is for a
    # non-diabetic patient at the median weight 82.45 kg and median eGFR 74
    # mL/min/1.73 m^2 (Table 2 footnote).
    lcl <- log(8.63)
    label("Apparent clearance CL/F, non-diabetic, WT 82.45 kg, eGFR 74 (L/h)")
    # Table 2 row 'CL/F - theta1 (L/h)' = 8.63 (RSE 4%); Results equation.
    lvc <- log(12.5)
    label("Apparent volume of distribution Vd/F (L)")
    # Table 2 row 'Vd/F - theta2 (L)' = 12.5 (RSE 10%).
    lka <- log(0.131)
    label("First-order absorption rate constant (1/h)")
    # Table 2 row 'Ka - theta3 (h-1)' = 0.131 (RSE 6%).
    ltlag <- log(0.165)
    label("Absorption lag time (h)")
    # Table 2 row 'TLAG - theta4 (h)' = 0.165 (RSE 3%).

    # Covariate effects on CL/F -- Results final-model equation:
    # CL/F = 8.63 * (Weight/82.45)^0.963 * (eGFR/74)^0.56 * 0.682^Diabetes * exp(eta)
    e_wt_cl <- 0.963
    label("Power exponent on (WT / 82.45) for CL/F (unitless)")
    # Table 2 row 'Weight effect on CL/F - theta5' = 0.963 (RSE 15%).
    e_crcl_cl <- 0.56
    label("Power exponent on (CRCL / 74) for CL/F (unitless)")
    # Table 2 row 'eGFR effect on CL/F - theta6' = 0.56 (RSE 18%).
    e_diab_cl <- 0.682
    label("Multiplicative factor on CL/F for patients with diabetes (unitless)")
    # Table 2 row 'Diabetes effect on CL/F - theta7' = 0.682 (RSE 8%).

    # Interindividual variability: exponential, on CL/F only (Results:
    # diagonal omega; only omega2 CL/F is reported in Table 2).
    etalcl ~ 0.138
    # Table 2 row 'omega2 CL/F' = 0.138 (RSE 7%), a variance.

    # Residual error: mixed (additive + proportional) model (Results; Methods
    # form iii, epsilon = epsilon1 + F * epsilon2). Table 2 reports VARIANCES
    # (sigma2), so the SDs are their square roots.
    addSd <- 2.345208
    label("Additive residual error SD (ng/mL)")
    # Table 2 row 'sigma2 (additive)' = 5.5 (RSE 20%); sqrt(5.5) = 2.3452 ng/mL.
    propSd <- 0.646529
    label("Proportional residual error SD (fraction)")
    # Table 2 row 'sigma2 (proportional)' = 0.418 (RSE 3%); sqrt(0.418) = 0.6465.
  })

  model({
    # Individual PK parameters. Normalising constants 82.45 kg and 74
    # mL/min/1.73 m^2 are the cohort medians (Table 2 footnote).
    cl <- exp(lcl + etalcl) *
      (WT / 82.45)^e_wt_cl *
      (CRCL / 74)^e_crcl_cl *
      e_diab_cl^DIS_DIAB
    vc <- exp(lvc)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    kel <- cl / vc

    # One compartment with first-order absorption and lag (Results: 'best
    # described by a one-compartment model with first-order absorption and
    # first-order elimination. Adding absorption lag time ... improved model
    # fit').
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L; x1000 gives ng/mL, the assay unit
    # (Methods: 1.00-250 ng/mL) and the unit of the additive residual error.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
