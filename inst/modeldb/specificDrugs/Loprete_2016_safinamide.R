Loprete_2016_safinamide <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model with first-order absorption ",
    "and first-order elimination for oral safinamide in patients with Parkinson's ",
    "disease on stable dopamine-agonist or levodopa therapy, pooled from two phase 3 ",
    "randomized placebo-controlled trials (Study 015 and Study 016; Loprete 2016 ",
    "model 110a). Apparent oral clearance and apparent volume of distribution are ",
    "allometrically scaled on body weight about a 70 kg reference with structural ",
    "exponents of 0.75 and 1. Safinamide exposure was about 30% lower in Study 016 ",
    "(tablets of safinamide free base) than in Study 015 (gelatin capsules of ",
    "safinamide methanesulfonate); the authors carried that difference as a relative ",
    "bioavailability factor of 0.724 on the depot rather than as separate clearance ",
    "and volume estimates, switched on by the binary STUDY_016 indicator, so the ",
    "default ini() values are the Study 015 reference parameters. Interindividual ",
    "variability is estimated on CL/F and Vd/F; the authors set the IIV on KA and the ",
    "additive residual error term to zero in the final model because both were poorly ",
    "estimated. Residual error is proportional. Age, sex, creatinine clearance, race ",
    "and levodopa exposure were screened and not retained."
  )
  reference <- paste(
    "Loprete L, Leuratti C, Cattaneo C, Thapar MM, Farrell C, Sardina M.",
    "Population pharmacokinetic and pharmacodynamic analyses of safinamide",
    "in subjects with Parkinson's disease.",
    "Pharmacol Res Perspect 2016;4(5):e00251.",
    "doi:10.1002/prp2.251.",
    sep = " "
  )
  vignette <- "Loprete_2016_safinamide"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-varying in the source dataset (recorded at each study visit; Loprete 2016 ",
        "Methods, PK analysis dataset, describes the imputation rules used when a visit ",
        "value was missing). Enters both CL/F and Vd/F as the allometric power term ",
        "(WT/70)^exponent with a 70 kg reference. Loprete 2016 Methods, Covariate model, ",
        "states the reference values were 'the values regarded as reference or normal in ",
        "the general population', and that 70 kg was used for body weight even though the ",
        "cohort median was 64 kg."
      ),
      source_name = "WGT"
    ),
    STUDY_016 = list(
      description = paste0(
        "Binary indicator of enrolment in phase 3 Study 016 (Borgohain 2014; patients with ",
        "motor fluctuations on stable levodopa, dosed with tablets containing safinamide ",
        "free base). 0 = Study 015 (Stocchi 2012; early Parkinson's disease on a single ",
        "dopamine agonist, dosed with gelatin capsules containing safinamide ",
        "methanesulfonate), which is the reference study."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Study 015, the reference study, whose relative bioavailability is 1)",
      notes = paste0(
        "Time-fixed per subject. Enters the relative bioavailability of the depot as ",
        "Loprete 2016 equation 2, TVP = theta1 * theta2^IND, i.e. F = 1 * 0.724^STUDY_016. ",
        "The authors first tested the study effect separately on CL/F and Vd/F (model 112) ",
        "and found similar interstudy ratios (CL/F 1.3, Vd/F 1.6); because a single ",
        "bioavailability factor was simpler and more physiologically explainable they ",
        "adopted model 110 / 110a instead. Loprete 2016 Discussion notes the ~30% ",
        "difference matches the molecular-weight difference between safinamide free base ",
        "and safinamide methanesulfonate and suggests a systematic error, so the indicator ",
        "is confounded with formulation and with analytical laboratory; it is NOT a ",
        "population difference. Setting STUDY_016 = 1 reproduces the Study 016 apparent ",
        "parameters quoted in the Loprete 2016 abstract, CL/F = 3.59 / 0.724 = 4.96 L/h ",
        "and Vd/F = 120 / 0.724 = 166 L."
      ),
      source_name = "STUD"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = paste0(
        "Screened on CL/F and Vd/F against base model 112 with a 60-year reference (the ",
        "cohort median) and not retained; Loprete 2016 Results state the maximum OFV drop ",
        "across all screened covariates was 4 points, below the 6.63-point threshold. No ",
        "point estimate is published."
      ),
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste0(
        "Screened on CL/F and Vd/F and not retained (Loprete 2016 Results). The source ",
        "column GEND used the same polarity (0 = male, 1 = female), so no value ",
        "transformation would be needed. No point estimate is published."
      ),
      source_name = "GEND"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = paste0(
        "Screened on CL/F with a 75 mL/min reference (the cohort median) and not retained ",
        "(Loprete 2016 Results). The Discussion reads this as confirmation that renal ",
        "clearance contributes little to safinamide elimination, about 90% of the drug ",
        "being metabolized, and concludes no dose adjustment is needed in mild or moderate ",
        "renal impairment. No point estimate is published."
      ),
      source_name = "CRLR"
    ),
    DOSE_LEVODOPA_MGD = list(
      description = "Levodopa dose rate",
      units = "mg/24 h",
      type = "continuous",
      notes = paste0(
        "Screened as three separate columns and none retained (Loprete 2016 Results): BLEV ",
        "(levodopa exposure at baseline), LEVO (levodopa dose rate per 24 h at each visit, ",
        "reference 500 mg/24 h, the cohort median) and LEVR (rate of change from baseline, ",
        "1 at weeks 1 and 4 by construction). All three were 0 for Study 015, where ",
        "levodopa was not administered. No point estimate is published; not registered in ",
        "inst/references/covariate-columns.md because it is screened-only and never ",
        "referenced in model()."
      ),
      source_name = "BLEV / LEVO / LEVR"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = paste0(
        "Not entered in the covariate model at all because race was recorded only in Study ",
        "016 (Loprete 2016 Methods, Covariate model; Table 1 shows 177/177 Study 015 ",
        "patients with race not available). The authors did inspect plots of the ",
        "individual random effects against race and saw no visible trend. No point estimate ",
        "is published."
      ),
      source_name = "RACE"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "safinamide",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "safinamide",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 623L,
    n_studies = 2L,
    age_range = "31.0-82.0 years",
    age_median = "60.0 years",
    weight_range = "33.5-102 kg",
    weight_median = "64.0 kg",
    sex_female_pct = 30.8,
    race_ethnicity = c(Asian = 57.5, White = 13.9, Other = 0.2, NotAvailable = 28.4),
    disease_state = paste0(
      "Parkinson's disease. Study 015 enrolled early idiopathic Parkinson's disease on a ",
      "stable dose of a single dopamine agonist; Study 016 enrolled Parkinson's disease ",
      "with motor fluctuations on a stable dose of levodopa."
    ),
    dose_range = paste0(
      "50-200 mg/day oral. Study 015: low dose titrated 50 to 100 mg/day, or high dose ",
      "titrated 100 to 150 to 200 mg/day, over 24 weeks. Study 016: 50 or 100 mg/day over ",
      "24 weeks. Placebo-arm records were excluded from the PK dataset."
    ),
    renal_function = paste0(
      "No severe renal impairment (creatinine clearance < 10 mL/min). 57 (9%) moderate ",
      "(< 50 mL/min), 307 (49%) mild (50-80 mL/min), 259 (42%) normal. Baseline creatinine ",
      "clearance median 75.5 mL/min, range 25.3-201 mL/min."
    ),
    regions = "Not reported in Loprete 2016. Study 016 was predominantly Asian (80.3%).",
    notes = paste0(
      "Demographics from Loprete 2016 Table 1 (All, n = 623). The analysis dataset held ",
      "624 patients contributing 2785 concentration records; 623 patients contributing ",
      "2719 records entered the population PK analysis (177 patients / 1099 concentrations ",
      "from Study 015 and 446 patients / 1620 concentrations from Study 016). Sex split ",
      "431 male (69%) / 192 female (31%). Race was available only for Study 016; the ",
      "Table 1 'All' race percentages therefore carry 28.4% not available. Bioanalysis was ",
      "by validated LC-MS/MS with a lower limit of quantification of 20 ng/mL; all ",
      "below-quantification-limit records were excluded rather than imputed."
    )
  )

  ini({
    # Structural parameters. Loprete 2016 Table 3 reports the final model
    # (110a) estimates for STUDY 015, i.e. with relative bioavailability
    # anchored at 1, for a 70 kg reference patient (Table 3 footnote:
    # 'The reference population for PK parameters CL/F and Vd/F is a 70 kg
    # patient').
    lka <- log(0.582)
    label("Absorption rate constant (KA, 1/h)")
    # Table 3 row 'KA (h-1)': 0.582, %RSE 21.6, 95% CI 0.335-0.829. Bootstrap median 0.572.

    lcl <- log(3.59)
    label("Apparent oral clearance, 70 kg patient in Study 015 (CL/F, L/h)")
    # Table 3 row 'CL/F (L/h)': 3.59, %RSE 2.67, 95% CI 3.40-3.78. Bootstrap median 3.58.

    lvc <- log(120)
    label("Apparent volume of distribution, 70 kg patient in Study 015 (Vd/F, L)")
    # Table 3 row 'Vd/F (L)': 120, %RSE 4.38, 95% CI 110-130. Bootstrap median 120.

    # Relative bioavailability. Loprete 2016 Methods, Structural model:
    # 'Inclusion of a relative bioavailability factor (Frel), where Frel was
    # set to 1 for the Study 015 and an estimate for Study 016, was
    # evaluated.' Study 015 is therefore the structural anchor and only the
    # Study 016 value was estimated.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability in Study 015, the reference study (unitless)")
    # Methods, Structural model: Frel set to 1 for Study 015; anchor, not an estimate.

    e_study_016_fdepot <- log(0.724)
    label("Log relative bioavailability of Study 016 versus Study 015 (unitless)")
    # Table 3 row 'F (STUD016)': 0.724, %RSE 2.49, 95% CI 0.689-0.759. Bootstrap median 0.725. Enters via Loprete 2016 equation 2, TVP = theta1 * theta2^IND.

    # Allometric exponents. Loprete 2016 Results, PK analysis and final
    # model: the base model started from 'an allometrically scaled CL/F
    # (with a scaling factor (WGT/70)^0.75) and Vd/F (with a scaling factor
    # (WGT/70))'; the Discussion restates this as 'linearly for Vd/F and
    # with a power coefficient of 0.75 for CL/F'. Neither exponent appears
    # in Table 3 with a point estimate, %RSE or confidence interval, so both
    # are structural rather than estimated.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on body weight for CL/F (unitless)")
    # Results, PK analysis and final model: scaling factor (WGT/70)^0.75 on CL/F.

    e_wt_vc <- fixed(1)
    label("Allometric exponent on body weight for Vd/F (unitless)")
    # Results, PK analysis and final model: scaling factor (WGT/70) on Vd/F, i.e. linear; restated in the Discussion.

    # Interindividual variability. Loprete 2016 equation 3 is the
    # exponential (log-normal) error model P_i = TVP * exp(eta_i), so the
    # Table 3 'x2' column holds variances on the log scale. The Table 3 CV%
    # column is the square root of that variance expressed as a percentage
    # (sqrt(0.0772) = 0.278 and sqrt(0.0892) = 0.299), which fixes the
    # column as a variance rather than a standard deviation.
    etalcl ~ 0.0772
    # Table 3 row 'x2 CL': 0.0772, %RSE 11.2, 95% CI 0.0602-0.0942, CV 27.8%. Bootstrap median 0.0759.
    etalvc ~ 0.0892
    # Table 3 row 'x2 Vd': 0.0892, %RSE 17.2, 95% CI 0.0592-0.119, CV 29.9%. Bootstrap median 0.0855.
    etalka ~ fixed(0)
    # Results, model 110a: IIV on KA was poorly estimated (%RSE 71, 95% CI included the null, shrinkage 60%) and was set to zero in the final model. Not reported in Table 3.

    # Residual error. Loprete 2016 equation 4 is combined proportional plus
    # additive; the additive component was set to zero in the final model,
    # leaving a proportional-only error. Table 3 reports the variance
    # (sqrt(0.0885) = 0.297, matching the printed CV of 29.7%), so the
    # nlmixr2 standard deviation is its square root.
    propSd <- sqrt(0.0885)
    label("Proportional residual error standard deviation (fraction)")
    # Table 3 row 'r2 prop': 0.0885, %RSE 5.11, 95% CI 0.0796-0.0974, CV 29.7%. Bootstrap median 0.0892. propSd = sqrt(0.0885) = 0.2975.

    addSd <- fixed(0)
    label("Additive residual error standard deviation (ng/mL)")
    # Results, model 110a: the additive residual error was poorly estimated (%RSE 102, 95% CI included the null) and was set to zero in the final model. Retained here so the published equation-4 combined error structure stays visible.
  })

  model({
    # Individual parameters. Body weight enters both disposition parameters
    # as the equation-1 power form (COV/COVST)^theta with COVST = 70 kg.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Relative bioavailability, Loprete 2016 equation 2 on the depot:
    # F = 1 * 0.724^STUDY_016. Study 015 (STUDY_016 = 0) gets F = 1, which
    # is the anchor the Table 3 CL/F and Vd/F estimates are conditioned on.
    f(depot) <- exp(lfdepot + e_study_016_fdepot * STUDY_016)

    # central is in mg and vc in L, so central/vc is in mg/L; multiply by
    # 1000 to report ng/mL, the unit of the source assay (lower limit of
    # quantification 20 ng/mL).
    Cc <- central / vc * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
