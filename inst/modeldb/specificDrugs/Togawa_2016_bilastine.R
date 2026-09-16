Togawa_2016_bilastine <- function() {
  description <- "Two-compartment population pharmacokinetic model with first-order absorption for oral bilastine in healthy adult Japanese male volunteers. Fit in NONMEM VI to pooled single-dose (10, 20, 50 mg) and 14-day once-daily multiple-dose (20, 50 mg) plasma data from a single-centre phase I study (45 bilastine-treated subjects, 1022 plasma observations). Parameters are apparent oral values (CL/F, Vc/F, Q/F, Vp/F) because no intravenous arm was studied; the absolute oral bioavailability of bilastine is reported elsewhere as 60.67%. The authors describe the structure as a 'two-compartment, semi-physiological parameter' model, meaning the clearance / volume parameterisation rather than disposition micro-constants. Inter-individual variability was carried as a full 4 x 4 random-effects block across CL, Vc, Q and Vp plus a separate random effect on ka; only the diagonal elements are published, so the off-diagonal covariances cannot be reproduced here. Residual error is proportional-only, selected over an additive model on objective function (5487.4 vs 9669.6). No covariate effect reached significance on any pharmacokinetic parameter."
  reference <- paste(
    "Togawa M, Yamaya H, Rodriguez M, Nagashima H (2016).",
    "Pharmacokinetics, pharmacodynamics and population pharmacokinetic/pharmacodynamic",
    "modelling of bilastine, a second-generation antihistamine, in healthy Japanese subjects.",
    "Clin Drug Investig 36(12):1011-1021.",
    "doi:10.1007/s40261-016-0447-2.",
    sep = " "
  )
  vignette <- "Togawa_2016_bilastine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
    # Doses are administered as mg tablets and plasma bilastine is assayed in
    # ng/mL (LC-MS/MS calibration range 0.2-400 ng/mL), so the amount unit that
    # makes `Cc <- central / vc` come out in ng/mL is ug: 20 mg = 20000 ug
    # divided by Vc/F = 51.2 L gives ug/L = ng/mL. Event tables must therefore
    # dose in ug (see the vignette), which is why compartmentData below records
    # the state amount unit as ug rather than mg.
  )

  compartmentData <- list(
    depot = list(analyte = "bilastine", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bilastine", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "bilastine", units = "ug", specimen = "plasma", verified = TRUE)
  )

  # No covariate was retained: "No covariate effects were found in any of the
  # pharmacokinetic parameters" (Togawa 2016 Results 3.4). The ten covariates
  # the paper screened are documented below so the provenance of the screen
  # survives, without triggering a "declared but not referenced" warning.
  covariateData <- list()

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years at screening. Study cohort 20-39 years by inclusion criterion; per-arm means 22.9-29.8 years (Togawa 2016 Table 1).",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Togawa 2016 Methods 2.5) against empirical-Bayes individual parameters, then for significance within the population model at p < 0.05; not retained in the final model."
    ),
    WT = list(
      description = "Body weight in kg at screening. Per-arm means 61.5-64.2 kg across the bilastine arms (Togawa 2016 Table 1); inclusion required weight >= 50 kg and BMI 18.5 to < 25 kg/m^2.",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Togawa 2016 Methods 2.5); not retained in the final model. Continuous covariates were centred on their median values so that the intercept represents the parameter at median covariate values, but no continuous covariate survived selection."
    ),
    HT = list(
      description = "Standing height in cm at screening. Per-arm means 170.2-173.9 cm (Togawa 2016 Table 1).",
      units = "cm",
      type = "continuous",
      notes = "Screened in the stepwise covariate search (Togawa 2016 Methods 2.5); not retained in the final model."
    ),
    ALB = list(
      description = "Serum albumin at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "g/L",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5 ('Albumin'); not retained in the final model. The paper reports no numeric albumin distribution, so no units can be confirmed from the source -- the canonical SI unit is recorded here for documentation only."
    ),
    TBILI = list(
      description = "Total bilirubin at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "umol/L",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5 ('Bilirubin', unqualified, read here as total bilirubin); not retained in the final model. No numeric distribution is published, so the canonical SI unit is recorded for documentation only."
    ),
    AST = list(
      description = "Aspartate transaminase at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "U/L",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5; not retained in the final model."
    ),
    ALT = list(
      description = "Alanine aminotransferase at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "U/L",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5; not retained in the final model."
    ),
    BUN = list(
      description = "Blood urea nitrogen at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "mg/dL",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5; not retained in the final model. No numeric distribution is published, so the unit is recorded for documentation only."
    ),
    ALP = list(
      description = "Alkaline phosphatase at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "U/L",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5; not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine at screening. Screened as a candidate covariate; the cohort values are not tabulated in the paper.",
      units = "umol/L",
      type = "continuous",
      notes = "Listed among the screened covariates in Togawa 2016 Methods 2.5; not retained in the final model. Relevant because bilastine is largely eliminated unchanged (44.8-51.9% of the dose recovered in urine over 72 h; Togawa 2016 Table 2), but no renal-function covariate reached significance in this healthy cohort with a narrow creatinine range."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    age_range = "20-39 years by inclusion criterion; per-arm means 22.9-29.8 years (SD 2.1-6.9)",
    weight_range = ">= 50 kg by inclusion criterion; per-arm means 61.5-64.2 kg (SD 5.19-8.22)",
    height_range = "per-arm means 170.2-173.9 cm (SD 4.70-7.55)",
    bmi_range = "18.5 to < 25 kg/m^2 by inclusion criterion; per-arm means 20.6-21.8 kg/m^2",
    sex_female_pct = 0,
    race_ethnicity = c(Japanese = 100),
    disease_state = "Healthy volunteers (no acute or chronic disease; no history of drug allergy; non-smokers for >= 90 days; normal histamine skin-prick reactivity at screening).",
    dose_range = "Oral bilastine tablets. Part I: single doses of 10, 20 or 50 mg under fasting conditions (n = 9 per dose). Part II: 20 or 50 mg once daily for 14 days (n = 9 per dose). Placebo arms are excluded from the PK analysis.",
    regions = "Japan (single centre, Tokyo)",
    n_observations = c(plasma = 1022L),
    notes = paste(
      "Togawa 2016 randomised 60 healthy Japanese male subjects (36 in Part I,",
      "24 in Part II); 45 of them received bilastine and contribute the 1022",
      "plasma observations in the population PK analysis (Table 5 footnote a).",
      "One Part II 20 mg subject withdrew on Day 4 with moderate gastroenteritis",
      "judged unrelated to study drug. Demographics are from Table 1, which",
      "reports them per dose arm rather than pooled; the ranges above are the",
      "span of the per-arm means for the six bilastine arms. The paper reports",
      "no pooled min-max for age, weight or height."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Structural parameters (Togawa 2016 Table 5, 'Japanese' column).
    # All are APPARENT oral values: the study had no intravenous arm,
    # and the table footnote defines CL as "apparent total clearance of
    # the drug from plasma". Dose / CL reproduces the published NCA
    # AUC0-inf to within 2% at every dose level (see vignette).
    # ---------------------------------------------------------------
    lka <- log(1.7); label("First-order absorption rate constant ka (1/h)") # Togawa 2016 Table 5: ka = 1.7 1/h (%SEE 7)
    lcl <- log(14.4); label("Apparent oral clearance CL/F (L/h)") # Togawa 2016 Table 5: CL = 14.4 L/h (%SEE 4)
    lvc <- log(51.2); label("Apparent central volume of distribution Vc/F (L)") # Togawa 2016 Table 5: Vc = 51.2 L (%SEE 5)
    lq <- log(1.55); label("Apparent intercompartmental clearance Q/F (L/h)") # Togawa 2016 Table 5: Q = 1.55 L/h (%SEE 8)
    lvp <- log(20.2); label("Apparent peripheral volume of distribution Vp/F (L)") # Togawa 2016 Table 5: Vp = 20.2 L (%SEE 9)

    # ---------------------------------------------------------------
    # Inter-individual variability (Togawa 2016 Table 5, rows
    # 'omega CL (%)' ... 'omega ka (%)'). The table reports omega itself
    # as a percentage, not omega^2 and not a back-transformed CV, so the
    # variance encoded here is (percent / 100)^2. Equation 1 of the
    # paper defines the IIV as exponential, under which omega = 0.28
    # corresponds to a coefficient of variation of 28.6% -- the two
    # readings differ by under 1 percentage point at these magnitudes.
    #
    # The published model used a FULL 4 x 4 random-effects block across
    # CL, Vc, Q and Vp plus a separate random effect on ka (Results 3.4).
    # Only the five diagonal elements are published; the six
    # off-diagonal covariances of the systemic block are not reported
    # anywhere in the paper or its Electronic Supplementary Material, so
    # they cannot be reproduced. The etas below are therefore encoded as
    # independent, which is a documented deviation (see vignette Errata)
    # and will under-state the correlation of CL with Vc, Q and Vp.
    # ---------------------------------------------------------------
    etalcl ~ 0.0784 # Togawa 2016 Table 5: omega CL = 28% (%SEE 22) -> variance 0.28^2
    etalvc ~ 0.1156 # Togawa 2016 Table 5: omega Vc = 34% (%SEE 19) -> variance 0.34^2
    etalq ~ 0.25 # Togawa 2016 Table 5: omega Q = 50% (%SEE 27) -> variance 0.50^2
    etalvp ~ 0.3721 # Togawa 2016 Table 5: omega Vp = 61% (%SEE 21) -> variance 0.61^2
    etalka ~ 0.0784 # Togawa 2016 Table 5: omega ka = 28% (%SEE 40) -> variance 0.28^2

    # ---------------------------------------------------------------
    # Residual error (Togawa 2016 Table 5, row 'sigma (%)'). Equation 2
    # of the paper defines a proportional-only model, and Results 3.4
    # records that it was selected over an additive model on objective
    # function (5487.4 proportional vs 9669.6 additive).
    # ---------------------------------------------------------------
    propSd <- 0.21; label("Proportional residual error (fraction)") # Togawa 2016 Table 5: sigma = 21% (%SEE 9)
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Disposition micro-constants from the clearance / volume
    # parameterisation the paper reports.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
