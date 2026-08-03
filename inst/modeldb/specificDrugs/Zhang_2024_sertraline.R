Zhang_2024_sertraline <- function() {
  description <- "One-compartment first-order absorption population PK model for sertraline in Chinese psychiatric inpatients (Zhang 2024). Apparent oral clearance decreases linearly with age around the 22-year cohort median (CL/F = 76.1 * [1 - 0.0068 * (AGE - 22)] L/h); the absorption rate constant is held at 0.098 1/h taken from Li 2013 because the therapeutic-drug-monitoring dataset contained trough samples only and the absorption phase was not identifiable."
  reference   <- paste(
    "Zhang Z, Guo Z, Tan Y, Li L, Wang Z, Wen Y, Huang S, Shang D.",
    "Population pharmacokinetic approach to guide personalized sertraline",
    "treatment in Chinese patients. Heliyon. 2024 Feb 1;10(3):e25231.",
    "doi:10.1016/j.heliyon.2024.e25231.",
    "The fixed absorption rate constant ka = 0.098 1/h is taken from",
    "Li CH, Pollock BG, Lyketsos CG, et al. Population pharmacokinetic modeling",
    "of sertraline treatment in patients with Alzheimer disease: the DIADS-2",
    "study. J Clin Pharmacol. 2013;53(2):234-239.",
    "doi:10.1177/0091270012445793 (Zhang 2024 reference 20).",
    sep = " "
  )
  vignette    <- "Zhang_2024_sertraline"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at the therapeutic-drug-monitoring sample",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model. Zhang 2024 Results 3.2: adding age",
        "to CL/F in the forward-inclusion step dropped the objective function value by 7.081",
        "(P < 0.01, df = 1) and age was the only covariate to clear the retention threshold.",
        "Entered on the median-centred linear scale of the paper's generic continuous-covariate",
        "form (Eq. 3), with the 22-year cohort median as the centring value (Table 1).",
        "Reported cohort range 11-79 years; 52/140 patients were younger than 18 years,",
        "75/140 were 18-65 years and 13/140 were older than 65 years.",
        "NOTE ON SIGN: the printed final-model equation in Results 3.2 is",
        "CL/F = 76.1 * [1 - 0.0068 * (AGE - 22)], i.e. clearance FALLS with increasing age,",
        "whereas the generic Eq. (3) is written with a plus sign and Table 2 reports",
        "theta_CL-AGE as the unsigned magnitude 0.0068. The minus form is the one encoded here:",
        "it is the equation the authors printed for this model, it matches the prose",
        "('sertraline clearance decreased progressively with aging', Abstract and Results 3.2),",
        "and it matches Fig. 3, where simulated 15-year-old adolescents reach LOWER trough",
        "concentrations than 45-year-old adults at every dose.",
        sep = " "
      ),
      source_name        = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Sex, 1 = female",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened on CL/F and V/F in the forward-inclusion step and not retained",
        "(Zhang 2024 Results 3.2: 'Gender, height, weight, BMI, and combined medications were",
        "not found to have a significant effect on the decrease in OFV values').",
        "The paper's categorical-covariate form (Eq. 4) codes 0 = male and 1 = female, which",
        "matches the SEXF orientation directly. Cohort was balanced at 70 male / 70 female",
        "(Table 1).",
        sep = " "
      ),
      source_name        = "Gender (Zhang 2024 Table 1)"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened via the median-centred linear form of Eq. (3) and not retained",
        "(Zhang 2024 Results 3.2). Cohort median (range) 60 (40-110) kg, mean (SD)",
        "61.54 (14.58) kg (Table 1). No allometric scaling was applied in the final model.",
        sep = " "
      ),
      source_name        = "Weight (Zhang 2024 Table 1)"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened via the median-centred linear form of Eq. (3) and not retained",
        "(Zhang 2024 Results 3.2). Cohort median (range) 166 (145-185) cm, mean (SD)",
        "165.46 (7.98) cm (Table 1).",
        sep = " "
      ),
      source_name        = "Height (Zhang 2024 Table 1)"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened via the median-centred linear form of Eq. (3) and not retained",
        "(Zhang 2024 Results 3.2). Table 1 reports median (range) 36 (25-62) and mean (SD)",
        "37 (7.92); these values are not internally consistent with the tabulated weight and",
        "height (60 kg at 166 cm gives 21.8 kg/m^2), so the Table 1 BMI row appears to be",
        "mislabelled or to hold a different quantity. Recorded here for provenance only -",
        "no BMI term enters the final model, so the inconsistency does not affect any",
        "encoded value.",
        sep = " "
      ),
      source_name        = "BMI (Zhang 2024 Table 1)"
    ),
    CONMED_LAMOTRIGINE = list(
      description        = "Concomitant lamotrigine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant lamotrigine)",
      notes              = paste(
        "Screened via the paper's categorical-covariate form (Eq. 4; 0 = no combination",
        "therapy, 1 = combination therapy) and not retained (Zhang 2024 Results 3.2).",
        "The most common comedication in the cohort: 99 records (33.22 %) per Table 1.",
        sep = " "
      ),
      source_name        = "Lamotrigine (Zhang 2024 Table 1)"
    ),
    CONMED_QUETIAPINE = list(
      description        = "Concomitant quetiapine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant quetiapine)",
      notes              = paste(
        "Screened via the paper's categorical-covariate form (Eq. 4) and not retained",
        "(Zhang 2024 Results 3.2). 12 records (4.03 %) per Table 1.",
        sep = " "
      ),
      source_name        = "Quetiapine (Zhang 2024 Table 1)"
    ),
    CONMED_VENLAFAXINE = list(
      description        = "Concomitant venlafaxine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant venlafaxine)",
      notes              = paste(
        "Screened via the paper's categorical-covariate form (Eq. 4) and not retained",
        "(Zhang 2024 Results 3.2). 4 records (1.34 %) per Table 1.",
        sep = " "
      ),
      source_name        = "Venlafaxine (Zhang 2024 Table 1)"
    ),
    CYP2C19_PHENO = list(
      description        = "CYP2C19 metaboliser phenotype group",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "Normal / extensive metaboliser",
      notes              = paste(
        "Zhang 2024 Discussion: CYP2C19 was included as a candidate covariate at the start of",
        "the analysis but produced no effect on clearance and the result is not shown in the",
        "paper; the authors attribute this to the small genotyped subset (n = 38 of 140).",
        "No point estimate is published, so nothing is encoded. Recorded here to preserve the",
        "screen. N-desmethyl sertraline concentrations were likewise unavailable and were not",
        "screened at all.",
        sep = " "
      ),
      source_name        = "CYP2C19 (Zhang 2024 Discussion)"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "sertraline", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sertraline", units = "mg", specimen = "serum",               verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 140L,
    n_studies        = 1L,
    n_concentrations = 298L,
    age_range        = "11-79 years",
    age_median       = "22 years",
    weight_range     = "40-110 kg",
    weight_median    = "60 kg",
    height_range     = "145-185 cm",
    height_median    = "166 cm",
    sex_female_pct   = 50,
    race_ethnicity   = "Chinese (single-centre Guangzhou cohort; not further stratified by the source)",
    disease_state    = "Hospitalised patients with psychiatric disorders (depression, panic disorder, generalized social anxiety disorder, obsessive-compulsive disorder and related indications) receiving oral sertraline; 52/140 (37%) were adolescents younger than 18 years",
    dose_range       = "Not reported for the analysis dataset; the source states sertraline is routinely dosed at 25-200 mg/d and simulates 25-250 mg/d QD and 50-200 mg/d BID",
    regions          = "China (The Affiliated Brain Hospital of Guangzhou Medical University, Guangzhou, Guangdong)",
    co_medication    = "Lamotrigine 99 records (33.22%); quetiapine 12 records (4.03%); venlafaxine 4 records (1.34%)",
    sampling         = "Retrospective therapeutic drug monitoring; all 298 samples are elimination-phase trough concentrations. Each patient contributed non-zero concentrations at at least two different dose levels.",
    notes            = paste(
      "Baseline demographics from Zhang 2024 Table 1. Retrospective TDM data collected",
      "2018-2022; IRB approval 2021027. Serum sertraline quantified by LC-MS/MS over a",
      "5-500 ng/mL calibrated range. Because only troughs were available the absorption",
      "phase could not be estimated and ka was held at a literature value (Methods 2.3).",
      sep = " "
    )
  )

  ini({
    # Structural parameters (Zhang 2024 Table 2 'Final Model - Estimate' column,
    # cross-checked against the printed final-model equations in Results 3.2:
    #   CL/F = 76.1 * [1 - 0.0068 * (AGE - 22)];  V/F = 803;  Ka = 0.098).
    # The paper reports CL/F in L/h and V/F in L: Discussion quotes CL/V = 0.095
    # for this model, and 76.1 / 803 = 0.0948 1/h reproduces it.
    lka <- fixed(log(0.098))
    label("First-order absorption rate constant ka (1/h); literature value taken from Li 2013 DIADS-2")  # Zhang 2024 Table 2 'K = 0.098 FIX'; Methods 2.3 fixes ka to the value of reference 20 (Li 2013) because only trough samples were available
    lcl <- log(76.1)
    label("Apparent oral clearance CL/F at the 22-year median age (L/h)")  # Zhang 2024 Table 2: CL/F = 76.1, RSE 7%; bootstrap median 74.68 (95% CI 67.22-87.71)
    lvc <- log(803)
    label("Apparent volume of distribution V/F (L)")  # Zhang 2024 Table 2: V/F = 803, RSE 26%; bootstrap median 809.443 (95% CI 372.19-1369.47)

    # Covariate effect. Encoded with the minus sign of the printed final-model
    # equation (Results 3.2), which agrees with the prose and with Fig. 3; the
    # generic Eq. (3) is written with a plus sign and Table 2 reports the
    # unsigned magnitude. See the AGE entry in covariateData for the full
    # sign reconciliation.
    e_age_cl <- 0.0068
    label("Linear fractional decrease in CL/F per year of age above the 22-year median (1/year)")  # Zhang 2024 Table 2: theta_CL-AGE = 0.0068, RSE 24%; bootstrap median 0.0063 (95% CI 0.004-0.0104)

    # Inter-individual variability. Table 2 reports IIV as CV% (10% on CL/F,
    # 57% on V/F) for the exponential IIV model of Eq. (1); converted to the
    # log-scale variance with omega^2 = log(CV^2 + 1), giving
    # log(0.10^2 + 1) = 0.00995 and log(0.57^2 + 1) = 0.28134.
    etalcl ~ 0.00995  # Zhang 2024 Table 2 'IIV (CV%)' = 10 on CL/F
    etalvc ~ 0.28134  # Zhang 2024 Table 2 'IIV (CV%)' = 57 on V/F

    # Residual error. The final model is proportional only: Results 3.2 states
    # "We fixed the summation error to 0 and used a proportional type error
    # model", i.e. the additive term eps2 of Eq. (2) is fixed at zero.
    # Table 2 gives the proportional term as 0.129 (RSE 12%, bootstrap median
    # 0.128, 95% CI 0.088-0.166). The row is labelled "PRO (CV%)" but the value
    # sits in the raw-estimate column alongside 76.1 / 803 / 0.098 / 0.0068, so
    # it is the NONMEM $SIGMA variance rather than a percentage; the residual SD
    # is therefore sqrt(0.129) = 0.359. Cross-checked against the paper's own
    # simulations: Fig. 4 shows a simulated adolescent CV of roughly 50-53%
    # across 25-250 mg/d, and Fig. 3 shows an adult 200 mg/d inter-quartile
    # ratio of about 2.15. Reproducing the model with the IIV above gives
    # CV 44% and an IQ ratio of 1.8 at propSd = 0.359, versus CV 28% and 1.44 at
    # propSd = 0.129 - the variance reading is the one consistent with both
    # published figures.
    propSd <- sqrt(0.129)
    label("Proportional residual error (fraction)")  # Zhang 2024 Table 2 'PRO (CV%)' = 0.129 read as the NONMEM $SIGMA variance; see the note above
  })

  model({
    # Individual parameters. Age enters CL/F on the median-centred linear scale
    # of Zhang 2024 Eq. (3), with the sign of the printed final-model equation.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (1 - e_age_cl * (AGE - 22))
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment first-order absorption with first-order elimination.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg and vc in L give central/vc in mg/L; the assay and the AGNP
    # therapeutic reference range (10-150 ng/mL) are in ng/mL, so scale by 1000.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
