Zhang_2024_aripiprazole <- function() {
  description <- "One-compartment first-order absorption population PK model for aripiprazole in Chinese adults with schizophrenia (Zhang 2024), built from routine therapeutic-drug-monitoring trough concentrations. Apparent oral clearance is allometrically scaled on body weight (exponent 0.75, 70 kg reference) and reduced 28.6% by concomitant fluoxetine, a CYP2D6 inhibitor, giving a with:without fluoxetine clearance ratio of 0.714:1; the apparent volume of distribution scales linearly with weight and the absorption rate constant is held at 1.06 1/h from Kim 2008 because the dataset contained trough samples only."
  reference <- paste(
    "Zhang C, Jiang L, Hu K, Zhang Y-J, Han J, Chen J, Bulubu, Dong B, Shi H-Z,",
    "He S-M, Yu T-T, Chen X, Wang D-D. Drug-drug interaction and initial dosage",
    "optimization of aripiprazole in patients with schizophrenia based on",
    "population pharmacokinetics. Front Psychiatry. 2024 Jun 18;15:1377268.",
    "doi:10.3389/fpsyt.2024.1377268.",
    "The fixed absorption rate constant ka = 1.06 1/h is Zhang 2024 reference 14,",
    "Kim JR, Seo HB, Cho JY, et al. Population pharmacokinetic modelling of",
    "aripiprazole and its active metabolite, dehydroaripiprazole, in psychiatric",
    "patients. Br J Clin Pharmacol. 2008;66(6):802-810.",
    "doi:10.1111/j.1365-2125.2008.03223.x; see modellib('Kim_2008_aripiprazole').",
    sep = " "
  )
  vignette <- "Zhang_2024_aripiprazole"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only continuous covariate retained in the final model, and it enters both CL/F",
        "and V/F. Zhang 2024 Methods 2.2 Eq. (3) applies allometric scaling",
        "W_i = W_std * (X_i / X_std)^R with X_std = 70 kg and R = 0.75 for CL/F and 1 for V/F,",
        "citing Anderson & Holford as the source of the exponents; both exponents are",
        "structural assumptions rather than estimates and carry no standard error in Table 3.",
        "Weight was applied to every subject before the stepwise covariate search, so it is",
        "not one of the covariates that had to clear the OFV thresholds. Cohort weight was",
        "66.77 +/- 11.68 kg, median 67.00 kg, range 41.00-115.00 kg (Table 1); the dosing",
        "simulations of Figures 4 and 5 span 40-120 kg, slightly beyond the observed range.",
        sep = " "
      ),
      source_name        = "weight"
    ),
    CONMED_FLUOXETINE = list(
      description        = "Concomitant fluoxetine, 1 = receiving fluoxetine",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant fluoxetine)",
      notes              = paste(
        "The only concomitant medication retained in the final model. Zhang 2024 Results 3.2:",
        "of the 33 concomitant drugs screened in Table 2, only fluoxetine changed the objective",
        "function value enough to survive forward inclusion (OFV drop > 6.63) and backward",
        "elimination (OFV rise > 10.8). Entered on the paper's categorical-covariate form",
        "Eq. (5), Y_i = TV(Y) * (1 + theta * Cov_i), as printed in the final-model Eq. (6):",
        "CL/F = 3.23 * (weight/70)^0.75 * (1 - 0.286 * FLU), FLU = 1 when fluoxetine is",
        "co-administered and 0 otherwise. The paper's own summary of this effect - 'the",
        "aripiprazole clearance rates were 0.714:1 in patients with or without fluoxetine'",
        "(Abstract, Results 3.4, Discussion) - reproduces exactly as 1 - 0.286 = 0.714, which",
        "fixes both the sign and the multiplicative-fraction form of the term.",
        "Mechanism per Discussion: aripiprazole is cleared mainly by CYP2D6 and fluoxetine is",
        "a CYP2D6 inhibitor, so clearance falls and concentrations rise.",
        "CAUTION ON PRECISION: Table 2 records only 3 of the 119 patients as taking fluoxetine",
        "hydrochloride capsules, so theta_FLU rests on a very small exposed subgroup despite",
        "the 6.6% standard error reported in Table 3. Figure 3B shows the two observed",
        "concentration distributions (median roughly 290 ng/mL without fluoxetine versus",
        "roughly 550 ng/mL with, p < 0.01).",
        sep = " "
      ),
      source_name        = "FLU"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "aripiprazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "aripiprazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 119L,
    n_studies      = 1L,
    age_range      = "19.00-69.38 years",
    age_median     = "46.84 years",
    weight_range   = "41.00-115.00 kg",
    weight_median  = "67.00 kg",
    sex_female_pct = 52.1,
    disease_state  = "schizophrenia",
    dose_range     = "not reported; routine clinical aripiprazole dosing, oral tablet (95 patients), orally disintegrating tablet (35) or oral solution (1), 12 patients used two dosage forms",
    regions        = "China (Xuzhou Oriental Hospital Affiliated to Xuzhou Medical University, Jiangsu)",
    notes          = paste(
      "Retrospective analysis of the hospital therapeutic-drug-monitoring database, July 2020",
      "to June 2022 (Methods 2.1). 57 men and 62 women. Demographic and laboratory data are",
      "Table 1 and the 33 screened concomitant medications are Table 2. Samples were sparse",
      "real-world trough concentrations, which the authors themselves flag as limiting the",
      "predictive ability of the model (Discussion). Dosage form and every laboratory index",
      "screened (albumin, globulin, alanine and aspartate transaminase, creatinine, urea,",
      "total protein, total cholesterol, triglyceride, direct and total bilirubin, hematocrit,",
      "hemoglobin, mean corpuscular hemoglobin and its concentration) were rejected by the",
      "stepwise covariate search; see covariatesDataExcluded.",
      sep = " "
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Collected (Methods 2.1, Table 1: 57 men / 62 women) and screened, but not retained:",
        "Results 3.2 states that apart from weight and fluoxetine, 'the aripiprazole dosage",
        "form, or the physiological and biochemical indices, or other concomitant medications",
        "were not included in the final model'. No point estimate is reported for sex.",
        sep = " "
      )
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort age 44.29 +/- 13.03 years,",
        "median 46.84, range 19.00-69.38 years (Table 1). No point estimate is reported.",
        sep = " "
      )
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort albumin 41.96 +/- 2.83 g/L,",
        "median 41.75, range 33.70-50.40 g/L (Table 1). No point estimate is reported.",
        sep = " "
      )
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort creatinine 64.41 +/- 14.15 umol/L,",
        "median 62.00, range 4.03-112.00 umol/L (Table 1). No point estimate is reported.",
        sep = " "
      )
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort ALT 26.16 +/- 19.36 IU/L,",
        "median 19.00, range 7.00-141.00 IU/L (Table 1). No point estimate is reported.",
        sep = " "
      )
    ),
    AST = list(
      description        = "Aspartate aminotransferase",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort AST 22.71 +/- 10.79 IU/L,",
        "median 20.00, range 10.00-96.00 IU/L (Table 1). No point estimate is reported.",
        sep = " "
      )
    ),
    TBIL = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort total bilirubin 10.58 +/- 4.89",
        "umol/L, median 9.50, range 3.50-34.10 umol/L; direct bilirubin was screened",
        "alongside it (3.79 +/- 1.79 umol/L, Table 1). No point estimate is reported.",
        sep = " "
      )
    ),
    HCT = list(
      description        = "Hematocrit",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened and not retained (Results 3.2). Cohort hematocrit 39.50 +/- 4.83%,",
        "median 39.05, range 29.90-53.10% (Table 1). Hemoglobin, mean corpuscular",
        "hemoglobin and mean corpuscular hemoglobin concentration were screened alongside",
        "it and were also rejected. No point estimate is reported.",
        sep = " "
      )
    ),
    CONMED_CLOZAPINE = list(
      description        = "Concomitant clozapine, 1 = receiving clozapine",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant clozapine)",
      notes              = paste(
        "The most frequent concomitant antipsychotic in the cohort (Table 2: 42/119 on",
        "clozapine tablets plus 9/119 on clozapine dispersible tablets) and the largest",
        "exposed subgroup of any screened comedication, yet it did not survive the stepwise",
        "covariate search (Results 3.2). Recorded here to preserve the provenance of that",
        "negative result: fluoxetine was retained on 3 exposed patients while clozapine was",
        "rejected on 51. No point estimate is reported.",
        sep = " "
      )
    )
  )

  ini({
    # Structural parameters: Zhang 2024 Table 3 'Estimate' column, cross-checked
    # against the printed final-model equations of Results 3.2,
    #   Eq. (6)  CL/F = 3.23 * (weight/70)^0.75 * (1 - 0.286 * FLU)
    #   Eq. (7)  V/F  = 157  * (weight/70)
    # Reported in L/h and L with dose in mg, so central/vc is mg/L and the
    # observation is scaled to the ng/mL of the TDM assay in model() below.
    lka <- fixed(log(1.06))
    label("First-order absorption rate constant ka (1/h); literature value from Kim 2008")  # Zhang 2024 Table 3 'Ka (h-1) 1.06 (fixed)'; Methods 2.2 fixes ka to the value of reference 14 (Kim 2008) because the TDM dataset held trough samples only
    lcl <- log(3.23)
    label("Apparent oral clearance CL/F at 70 kg without fluoxetine (L/h)")  # Zhang 2024 Table 3: CL/F = 3.23 L/h, SE 2.8%; bootstrap median 3.22 (95% CI 3.04-3.40); also the leading coefficient of Eq. (6)
    lvc <- log(157)
    label("Apparent volume of distribution V/F at 70 kg (L)")  # Zhang 2024 Table 3: V/F = 157 L, SE 15.3%; bootstrap median 160 (95% CI 113-235); also the leading coefficient of Eq. (7)

    # Allometric exponents. Methods 2.2 Eq. (3): "R is the allometric
    # coefficient: 0.75 for CL/F and 1 for V/F (15)" -- structural values taken
    # from the cited reference, not estimated, so neither appears in Table 3.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F for weight (unitless)")  # Zhang 2024 Methods 2.2 Eq. (3): R = 0.75 for CL/F, citing reference 15; the exponent is also visible in Eq. (6)
    e_wt_vc <- fixed(1)
    label("Allometric exponent on V/F for weight (unitless)")  # Zhang 2024 Methods 2.2 Eq. (3): R = 1 for V/F; Eq. (7) is written with the ratio to the first power

    # Covariate effect. Categorical form of Methods 2.2 Eq. (5),
    # Y_i = TV(Y) * (1 + theta * Cov_i), instantiated in Eq. (6) with the
    # minus sign shown explicitly: (1 - 0.286 * FLU). The paper's own summary
    # statistic, a with:without fluoxetine clearance ratio of 0.714:1,
    # reproduces exactly as 1 - 0.286 = 0.714.
    e_conmed_fluoxetine_cl <- -0.286
    label("Fractional change in CL/F with concomitant fluoxetine (unitless)")  # Zhang 2024 Table 3: theta_FLU = -0.286, SE 6.6%; bootstrap median -0.286 (95% CI -0.383 to -0.239); printed in Eq. (6)

    # Inter-individual variability, on CL/F only -- Table 3 reports no omega for
    # V/F and none for ka, which was fixed. Methods 2.2 Eq. (1) is the
    # exponential model J_i = TV(J) * exp(eta_i) and defines eta as having
    # "zero mean and variance omega^2", so the tabulated omega_CL/F = 0.233 is
    # the standard deviation on the log scale and the variance rxode2 needs is
    # 0.233^2 = 0.054289. See the vignette Errata for the two independent
    # checks that exclude reading 0.233 as the variance itself: the paper's own
    # Figure 5 target-attainment curves reproduce closely on the SD reading and
    # are unreachable on the variance reading, which caps attainment near 55%
    # against published values as high as 99%.
    etalcl ~ 0.233^2  # Zhang 2024 Table 3: omega_CL/F = 0.233, SE 12.2%; bootstrap median 0.233 (95% CI 0.172-0.291); squared per the Eq. (1) definition of omega as the SD

    # Residual error. Methods 2.2 Eq. (2) is the combined model
    # Q_i = P_i + P_i * eps1 + eps2, i.e. proportional plus additive, and
    # defines eps as having "zero mean and variance sigma^2" -- so the
    # tabulated sigma_1 and sigma_2 are standard deviations, read on the same
    # convention as omega above. sigma_2 is in the ng/mL of the assay.
    propSd <- 0.123
    label("Proportional residual error (fraction)")  # Zhang 2024 Table 3: sigma_1 = 0.123, SE 25.6%; bootstrap median 0.120 (95% CI 0.043-0.176); Table 3 footnote 'sigma_1, residual variability, proportional error'
    addSd <- 49.498
    label("Additive residual error (ng/mL)")  # Zhang 2024 Table 3: sigma_2 = 49.498, SE 19.4%; bootstrap median 49.498 (95% CI 14.614-64.440); Table 3 footnote 'sigma_2, residual variability, additive error'
  })

  model({
    # Individual parameters. Weight is allometric on both disposition
    # parameters (Eq. 3); fluoxetine multiplies clearance by the categorical
    # factor of Eq. (5), written out in Eq. (6).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl *
      (1 + e_conmed_fluoxetine_cl * CONMED_FLUOXETINE)
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # One-compartment model with first-order absorption and first-order
    # elimination (Methods 2.2: CL/F, V/F and a fixed Ka).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give central/vc in mg/L, whereas the TDM assay and
    # the 120-270 ng/mL therapeutic window of Methods 2.4 are in ng/mL, so
    # scale by 1000.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
