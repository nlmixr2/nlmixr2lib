Franke_2015_oxycodone_food <- function() {
  description <- "One-compartment population PK model with first-order absorption for oxycodone after a single oral dose of two biphasic immediate-release/extended-release oxycodone/acetaminophen (IR/ER OC/APAP) 7.5/325 mg tablets, fitted to the 31-subject FED-vs-FASTED food-effect study (Franke 2015 Table 5). This is a separate fit from the 151-subject pooled fasted analysis in modellib('Franke_2015_oxycodone') and has a different covariate model: apparent clearance scales with body HEIGHT referenced to the 168.5 cm cohort median, while the apparent volume scales with body weight referenced to 72.3 kg and carries a meal-type factor, and the absorption rate constant carries its own meal-type factor. A low-fat or high-fat meal cuts the absorption rate constant to roughly half of its fasted value. Residual error was not reported for this analysis and is packaged as fixed(0); see the vignette Errata. Concentrations are in ng/mL, so model() scales amount/volume by 1000 to convert mg/L to ng/mL."
  reference <- "Franke RM, Morton T, Devarakonda K. Pooled post hoc analysis of population pharmacokinetics of oxycodone and acetaminophen following a single oral dose of biphasic immediate-release/extended-release oxycodone/acetaminophen tablets. Drug Des Devel Ther. 2015;9:4587-4597. doi:10.2147/DDDT.S79499"
  vignette <- "Franke_2015_oxycodone_acetaminophen"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "oxycodone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oxycodone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline; linear ratio (HT/168.5) on CL/F, the reference being the 168.5 cm median height of the 31-subject food-effect cohort (Table 2). Height rather than weight was retained on CL/F in this fit, in contrast to the pooled fasted oxycodone model where weight was retained.",
      source_name = "HT"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline; linear ratio (WT/72.3) on V/F, the reference being the 72.3 kg median weight of the 31-subject food-effect cohort (Table 2, range 59.1-101.8).",
      source_name = "WT"
    ),
    FED_LOWFAT = list(
      description = "Low-fat, low-calorie meal at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted; FED_LOWFAT and FED_HIGHFAT are mutually exclusive and both 0 under fasting)",
      notes = "Protocol meal: approximately 25%-30% of kilocalories from fat, total 800 +/- 80 kcal (Methods 'Treatments'). Per dose record, not per subject: study 4 was a three-period six-sequence crossover in which each participant was dosed fasted, after a low-fat meal, and after a high-fat meal. Corresponds to the authors' VDIET = theta4 and KADIET = theta6 levels.",
      source_name = "DIET"
    ),
    FED_HIGHFAT = list(
      description = "High-fat, high-calorie meal at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted; FED_LOWFAT and FED_HIGHFAT are mutually exclusive and both 0 under fasting)",
      notes = "Protocol meal: approximately 50% of kilocalories from fat, total 1,000 +/- 100 kcal (Methods 'Treatments'), consistent with the FDA high-fat high-calorie definition. Per dose record, as for FED_LOWFAT. Corresponds to the authors' VDIET = theta5 and KADIET = theta7 levels.",
      source_name = "DIET"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained for oxycodone in the food-effect analysis (Table 5). Sex WAS retained on CL/F and V2/F for acetaminophen in the same analysis."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained in the food-effect oxycodone model (Table 5), although race was retained in the larger 151-subject pooled fasted fit. 8 of 31 participants were Black (25.8%, Table 2)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened but not retained (Table 5). Cohort mean 30.8 years (SD 10.1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained (Table 5). Cohort median 26.3 kg/m^2 (range 19.4-29.8)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31,
    n_studies = 1,
    age_range = "18-55 years (protocol inclusion range)",
    age_median = "30.8 years (mean, SD 10.1)",
    weight_range = "59.1-101.8 kg",
    weight_median = "72.3 kg",
    height_median = "168.5 cm (range 155.1-185.5)",
    bmi_median = "26.3 kg/m^2 (range 19.4-29.8)",
    sex_female_pct = 32.3,
    race_ethnicity = c(White = 74.2, Black = 25.8, Other = 0),
    disease_state = "Healthy adults enrolled in study 4, a single-center, open-label, randomized, phase 1, three-period, six-sequence crossover food-effect study.",
    dose_range = "A single oral dose of two intact IR/ER OC/APAP 7.5/325 mg tablets (15 mg oxycodone) on each of three occasions: fasted, after a low-fat meal, and after a high-fat meal.",
    regions = "United States",
    notes = "Baseline demographics from Franke 2015 Table 2 (completers, n = 31, of 48 enrolled). Protocol number COV15000244. Unlike the pooled fasted analysis, the food-effect data were analysed with first-order conditional estimation WITH INTERACTION."
  )

  ini({
    # Structural parameters - typical values at HT = 168.5 cm, WT = 72.3 kg,
    # fasted (FED_LOWFAT = FED_HIGHFAT = 0)
    lcl <- log(77.6);  label("Apparent clearance CL/F at 168.5 cm (L/h)")  # Table 5 final model, OC row 'CL/F (L/hour)'; same value in the printed CL/F equation
    lvc <- log(640);   label("Apparent volume of distribution V/F at 72.3 kg, fasted (L)")  # Table 5 final model, OC row 'V/F (L)'; same value in the printed V/F equation
    lka <- log(0.555); label("First-order absorption rate constant, fasted (1/h)")  # Table 5 final model, OC row 'Ka (hour-1)'; same value in the printed Ka equation

    # Meal-type factors. The authors define VDIET = 1 for fasting, theta4 for
    # low fat, theta5 for high fat; and KADIET = 1 for fasting, theta6 for low
    # fat, theta7 for high fat (legend beneath the food-effect equations).
    # These are dimensionless multipliers. Table 5 labels the two Ka rows
    # '(hour-1)', which cannot be right for a multiplier: read as absolute
    # rates they would imply only a 7%-10% food effect, whereas read as
    # multipliers the paired acetaminophen values 0.802 and 0.825 reproduce
    # that model's stated 20% and 18% reductions exactly.
    e_fed_lowfat_vc  <- 0.608; label("Low-fat-meal multiplicative factor on V/F (unitless)")  # printed V/F equation, '(1 or 0.608 or 0.648)'
    e_fed_highfat_vc <- 0.648; label("High-fat-meal multiplicative factor on V/F (unitless)")  # printed V/F equation, '(1 or 0.608 or 0.648)'
    e_fed_lowfat_ka  <- 0.518; label("Low-fat-meal multiplicative factor on Ka (unitless)")  # Table 5 final model, OC row 'Ka low fat'; same value in the printed Ka equation
    e_fed_highfat_ka <- 0.497; label("High-fat-meal multiplicative factor on Ka (unitless)")  # Table 5 final model, OC row 'Ka high fat'; same value in the printed Ka equation

    # IIV. Table 5 reports no variability rows, but the authors' filled-in
    # equations substitute the log-scale variance into the exp() that carries
    # each eta, the same print convention the variance-plus-CV% pairs of
    # Tables 3 and 4 establish for this paper.
    etalcl ~ 0.0517  # printed CL/F equation for the food-effect OC model, 'exp(0.0517)'
    etalvc ~ 0.0208  # printed V/F equation for the food-effect OC model, 'exp(0.0208)'
    etalka ~ 0.104   # printed Ka equation for the food-effect OC model, 'exp(0.104)'

    # Residual error. Table 5 footnote a states that a one-compartment additive
    # and proportional error model was used for OC, but NO residual magnitude
    # is reported anywhere for the food-effect analysis - not in the base model
    # column, not in the final model column, and not in the equations. The
    # structure is therefore declared with both magnitudes held at zero rather
    # than carrying over the unrelated 151-subject pooled-fasted estimates.
    addSd  <- fixed(0); label("Additive residual error (ng/mL); magnitude not reported for the food-effect analysis")
    propSd <- fixed(0); label("Proportional residual error (fraction); magnitude not reported for the food-effect analysis")
  })

  model({
    # 1. Derived covariate terms
    # Mutually exclusive meal indicators; both 0 selects the fasted reference,
    # reproducing the authors' VDIET and KADIET three-level factors.
    vdiet  <- 1 + (e_fed_lowfat_vc - 1) * FED_LOWFAT + (e_fed_highfat_vc - 1) * FED_HIGHFAT
    kadiet <- 1 + (e_fed_lowfat_ka - 1) * FED_LOWFAT + (e_fed_highfat_ka - 1) * FED_HIGHFAT

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (HT / 168.5)
    vc <- exp(lvc + etalvc) * (WT / 72.3) * vdiet
    ka <- exp(lka + etalka) * kadiet

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observation and error
    # central is in mg and vc in L, so central/vc is in mg/L = ug/mL;
    # multiply by 1000 to express Cc in ng/mL, the paper's assay unit.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
