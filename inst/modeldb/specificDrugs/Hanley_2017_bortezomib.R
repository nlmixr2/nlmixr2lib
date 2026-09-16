Hanley_2017_bortezomib <- function() {
  description <- "Pediatric population pharmacokinetics of intravenous bortezomib in children and adolescents aged 2-16 years with relapsed acute lymphoblastic leukemia or de novo acute myelogenous leukemia (Hanley 2017). Three-compartment model with linear elimination after twice-weekly IV bolus repeat dosing at 1.3 mg/m2; body surface area is the only retained covariate, acting as a power function on clearance (exponent 0.97, i.e. approximately linear, supporting BSA-based dosing) and on the second intercompartmental clearance (exponent 0.75); interindividual variability on CL, V1 and Q3; study-specific log-scale residual error reflecting the different sparse-sampling schemes of the two contributing Children's Oncology Group trials."
  reference <- "Hanley MJ, Mould DR, Taylor TJ, Gupta N, Suryanarayan K, Neuwirth R, Esseltine DL, Horton TM, Aplenc R, Alonzo TA, Lu X, Milton A, Venkatakrishnan K. Population Pharmacokinetic Analysis of Bortezomib in Pediatric Leukemia Patients: Model-Based Support for Body Surface Area-Based Dosing Over the 2- to 16-Year Age Range. The Journal of Clinical Pharmacology. 2017;57(9):1183-1193. doi:10.1002/jcph.906"
  vignette <- "Hanley_2017_bortezomib"

  # The paper estimates one log-scale residual SD per contributing study to
  # absorb the different sparse-sampling schemes (Results, 'Final model
  # development'); neither is the canonical single `expSd`.
  paper_specific_residual_sds <- c("expSdAall07p1", "expSdAaml1031")

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "bortezomib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bortezomib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "bortezomib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the final model, entering as a power function on CL",
        "(exponent 0.97) and on Q3 (exponent 0.75); Results 'Final model development' and Table 2.",
        "Normalisation constant encoded here as 1.30 m^2, the total-cohort mean BSA of Table 1.",
        "The paper's Methods are self-contradictory on this constant: the displayed covariate",
        "equation defines cov_i as 'the individual value for the covariate normalized for the",
        "population mean', while the sentence that follows states 'All body size metrics were",
        "referenced to the size of an average adult (weight, 70 kg; BSA, 1.8 m2)'. The",
        "population-mean reading is the one consistent with the paper's own reported derived",
        "quantities and the 1.8 m^2 reading is refuted by them -- see the vignette's",
        "'Assumptions and deviations' section for the four-way arithmetic check.",
        "The BSA formula (DuBois / Mosteller / Haycock) is not stated in the source."
      ),
      source_name = "BSA"
    ),
    STUDY_AALL07P1 = list(
      description = "Children's Oncology Group study AALL07P1 cohort indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study AAML1031, de novo AML)",
      notes = paste(
        "1 = subject enrolled in AALL07P1 (NCT00873093; phase 2; relapsed acute lymphoblastic",
        "leukemia; N = 51), 0 = subject enrolled in AAML1031 (NCT01371981; phase 3; de novo acute",
        "myelogenous leukemia; N = 53). Subject-level and time-fixed. Selects the residual-error",
        "magnitude only: the paper estimated separate residual variability per study to absorb the",
        "different sparse-sampling schemes (AALL07P1 lacks the 18-30 h postdose sample), which",
        "dropped the OFV by more than 52 points. Disease type (ALL vs AML) was screened as a",
        "covariate on the structural parameters and was NOT significant, so this indicator carries",
        "no structural effect."
      ),
      source_name = "Study"
    )
  )

  # Screened during covariate analysis but not retained in the final model, so
  # they carry no published point estimate and are documentation only
  # (Results 'Final model development': 'No additional statistically significant
  # covariates were identified').
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened on CL, V1 and Q3; the trends in the covariate-vs-eta plots were removed once BSA was on CL, and no residual age relationship remained (Figures 4 and 6).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and decreased the OFV, but BSA on CL gave the largest OFV drop and removed the weight-vs-eta trend, so weight was not retained.",
      source_name = "WT"
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or Other)",
      notes = "Race (White 65% / Black 18% / Other 17%, Table 1) was screened by forward addition and was not statistically significant.",
      source_name = "Race"
    ),
    RACE_OTHER = list(
      description = "Race-category 'Other' indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White or Black)",
      notes = "Screened with RACE_BLACK; not statistically significant.",
      source_name = "Race"
    ),
    DIS_AML = list(
      description = "Disease-type indicator, acute myelogenous leukemia versus acute lymphoblastic leukemia",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (relapsed ALL)",
      notes = "Screened as a covariate on CL, V1 and Q3 and not significant; VPCs stratified by disease type showed no heterogeneity (Results, 'The final population PK model was evaluated by VPC'). Distinct from STUDY_AALL07P1, which is retained but only selects the residual-error magnitude.",
      source_name = "Disease type (ALL or AML)"
    ),
    RISK_AML_HIGH = list(
      description = "High-risk AML indicator within study AAML1031",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (low-risk AML)",
      notes = "Risk group within the AML population (low-risk 75% / high-risk 25% of AAML1031, Table 1); screened and not significant.",
      source_name = "Risk group (Study AAML1031 only)"
    ),
    TXSTRAT_ALL = list(
      description = "Treatment-plan stratum within study AALL07P1",
      units = "(categorical)",
      type = "categorical",
      reference_category = "Pre-B ALL relapsing within 18 months of diagnosis",
      notes = "Three strata (pre-B ALL relapse within 18 months 39%, pre-B ALL relapse 18-36 months 51%, T-cell ALL / T-cell lymphoblastic lymphoma 10%; Table 1); screened and not significant.",
      source_name = "Treatment plan stratum (Study AALL07P1 only)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 104L,
    n_studies = 2L,
    n_observations = 571L,
    age_range = "2-16 years",
    age_mean = "10 years (study AALL07P1 8.5 y, study AAML1031 11.4 y)",
    age_groups = c(`2-11 years` = 60, `12-16 years` = 40),
    weight_range = "13.9-139.7 kg",
    weight_mean = "45.2 kg (study AALL07P1 40.6 kg, study AAML1031 49.6 kg)",
    bsa_range = "0.60-2.53 m^2",
    bsa_mean = "1.30 m^2 (study AALL07P1 1.20 m^2, study AAML1031 1.40 m^2)",
    sex_female_pct = 42,
    race_ethnicity = c(White = 65, Black = 18, Other = 17),
    disease_state = "Pediatric acute leukemia: relapsed acute lymphoblastic leukemia (study AALL07P1, N = 51) and de novo acute myelogenous leukemia (study AAML1031, N = 53). Bortezomib was given with multiagent chemotherapy backbones (ALL: vincristine, prednisone, asparaginase, doxorubicin, cyclophosphamide, etoposide; AML: cytarabine, daunorubicin, etoposide, mitoxantrone).",
    dose_range = "1.3 mg/m^2 IV bolus twice weekly, on days 1, 4 and 8 of the relevant blocks or courses of both trials plus day 11 of block 1 of AALL07P1",
    regions = "North America (Children's Oncology Group)",
    notes = paste(
      "Two Children's Oncology Group trials: AALL07P1 (NCT00873093, phase 2, relapsed ALL) and",
      "AAML1031 (NCT01371981, phase 3, de novo AML). PK sampling was on day 8, over 0-72 h",
      "postdose; AALL07P1 block 1 sampled predose, 5-15 min, 30-60 min, 4-8 h and ~72 h, while",
      "AAML1031 induction 2 additionally sampled 18-30 h. Median 6 concentrations per patient",
      "(range 1-6); no concentration was below the 0.1 ng/mL limit of quantification.",
      "IMPORTANT SCOPE LIMIT stated by the authors: the model was fit to REPEAT-dose data only",
      "and bortezomib clearance is known to fall between the first dose and steady state, so",
      "these parameters must not be used to describe first-dose exposure. The analysis also",
      "excludes patients under 2 years of age. Analysis in NONMEM 7 level 2 (FOCE with",
      "log-transform both sides); parameter precision from a 5000-run study-stratified",
      "nonparametric bootstrap, of which 4665 runs converged."
    )
  )

  ini({
    # --------------------------------------------------------------------
    # Structural parameters. Table 2 'Population Median (95%CI)' -- bootstrap
    # medians over the 4665 converged runs of a 5000-run bootstrap (Table 2
    # footnote a). Time h, doses mg, plasma concentrations ng/mL.
    #
    # Reference subject: BSA = 1.30 m^2 (the total-cohort mean BSA of Table 1).
    # See covariateData$BSA$notes and the vignette for why the population-mean
    # reading, not the 'average adult 1.8 m2' sentence, is the correct one.
    # --------------------------------------------------------------------
    lcl <- log(9.59)
    label("Clearance CL at the reference BSA (L/h)") # Table 2: CL 9.59 L/h (95%CI 8.79-10.37)
    e_bsa_cl <- 0.97
    label("Power exponent on (BSA/1.30) for CL (unitless)") # Table 2: 'BSA effect on CL' 0.97 (95%CI 0.72-1.25) -- approximately linear, the paper's support for BSA-based dosing

    lvc <- log(10.0)
    label("Central volume of distribution V1 (L)") # Table 2: V1 10.0 L (95%CI 6.09-13.4)

    lq <- log(25.8)
    label("Intercompartmental clearance Q2, central to peripheral1 (L/h)") # Table 2: Q2 25.8 L/h (95%CI 18.9-31.9)

    lvp <- log(32.5)
    label("First peripheral volume of distribution V2 (L)") # Table 2: V2 32.5 L (95%CI 23.1-43.1)

    lq2 <- log(26.6)
    label("Intercompartmental clearance Q3 at the reference BSA, central to peripheral2 (L/h)") # Table 2: Q3 26.6 L/h (95%CI 21.3-30.7)
    e_bsa_q2 <- 0.75
    label("Power exponent on (BSA/1.30) for Q3 (unitless)") # Table 2: 'BSA effect on Q3' 0.75 (95%CI 0.43-0.99)

    lvp2 <- log(975)
    label("Second peripheral volume of distribution V3 (L)") # Table 2: V3 975 L (95%CI 792-1190)

    # --------------------------------------------------------------------
    # Interindividual variability. The paper used an exponential (log-normal)
    # IIV model and reports the magnitudes as percentage coefficients of
    # variation (Results: 'The IIV percentage coefficient of variation (CV)
    # for all 3 parameters in the final model (CL 29.7%, V1 34.6%, Q3
    # 29.8%)'). Converted to the internal log-scale variance with the exact
    # log-normal relation omega^2 = log(CV^2 + 1).
    # --------------------------------------------------------------------
    etalcl ~ log(1 + 0.297^2) # Table 2 CL '%CV Interindividual Variance' 29.7 (95%CI 23.1-36.8); shrinkage ~14%
    etalvc ~ log(1 + 0.346^2) # Table 2 V1 '%CV Interindividual Variance' 34.6 (95%CI 13.5-59.9); shrinkage ~30%
    etalq2 ~ log(1 + 0.298^2) # Table 2 Q3 '%CV Interindividual Variance' 29.8 (95%CI 14.3-43.3); shrinkage ~30%

    # --------------------------------------------------------------------
    # Residual error. 'Modeling was performed using ... the log-transform both
    # sides approach for residual variability' and 'The residual error model
    # was a homoscedastic model on the log scale, proportional after back
    # transformation' -- i.e. additive on the log scale, which is lnorm() in
    # nlmixr2. One magnitude per study, because AALL07P1 lacks the 18-30 h
    # postdose sample. Reported as %CV and converted with the same exact
    # log-normal relation used for the IIV above: sd = sqrt(log(CV^2 + 1)).
    # --------------------------------------------------------------------
    expSdAall07p1 <- sqrt(log(1 + 0.468^2))
    label("Log-scale residual SD, study AALL07P1 (relapsed ALL)") # Table 2: 'Residual error for study AALL07P1,%CV' 46.8% (95%CI 37.5-57.7)
    expSdAaml1031 <- sqrt(log(1 + 0.219^2))
    label("Log-scale residual SD, study AAML1031 (de novo AML)") # Table 2: 'Residual error for study AAML1031,%CV' 21.9% (95%CI 16.0-29.1)
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters.
    #
    # BSA enters as a power function normalised to the 1.30 m^2 cohort-mean
    # reference, on CL and on Q3 only (Results 'Final model development';
    # Table 2). No covariate acts on V1, Q2, V2 or V3.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (BSA / 1.30)^e_bsa_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q2 <- exp(lq2 + etalq2) * (BSA / 1.30)^e_bsa_q2
    vp2 <- exp(lvp2)

    # ------------------------------------------------------------------
    # 2. Micro-constants and the three-compartment system with linear
    #    elimination from the central compartment. Dosing is IV bolus into
    #    `central`.
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ------------------------------------------------------------------
    # 3. Observation. `central` is an amount in mg and `vc` a volume in L, so
    #    central/vc is mg/L = ug/mL; the factor 1000 converts to ng/mL, the
    #    unit of the paper's LC-MS/MS assay (dynamic range 0.1-25.0 ng/mL,
    #    Methods 'Assessments').
    # ------------------------------------------------------------------
    Cc <- 1000 * central / vc

    # Study-specific log-scale residual SD.
    expSd <- STUDY_AALL07P1 * expSdAall07p1 + (1 - STUDY_AALL07P1) * expSdAaml1031
    Cc ~ lnorm(expSd)
  })
}
