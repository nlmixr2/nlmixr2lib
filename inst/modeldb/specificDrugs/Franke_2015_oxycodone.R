Franke_2015_oxycodone <- function() {
  description <- "One-compartment population PK model with first-order absorption for oxycodone after a single oral dose of biphasic immediate-release/extended-release oxycodone/acetaminophen (IR/ER OC/APAP) 7.5/325 mg tablets under FASTED conditions, pooled post hoc from four phase 1 studies in 151 healthy adults and nondependent recreational users of prescription opioids (Franke 2015 Table 3). Apparent clearance and apparent volume both carry a body-weight power term referenced to the 73.45 kg cohort median (exponent 0.75 on CL/F, 1 on V/F) and a multiplicative race factor written by the authors as theta^RACE, where the analysis dataset codes RACE as White = 0, Black = 1, Asian = 2 (Figure 4 note). Concentrations are in ng/mL, so model() scales amount/volume by 1000 to convert mg/L to ng/mL. The companion acetaminophen model from the same pooled fasted dataset is modellib('Franke_2015_acetaminophen'); the fed-vs-fasted oxycodone model from the separate 31-subject food-effect cohort is modellib('Franke_2015_oxycodone_food')."
  reference <- "Franke RM, Morton T, Devarakonda K. Pooled post hoc analysis of population pharmacokinetics of oxycodone and acetaminophen following a single oral dose of biphasic immediate-release/extended-release oxycodone/acetaminophen tablets. Drug Des Devel Ther. 2015;9:4587-4597. doi:10.2147/DDDT.S79499"
  vignette <- "Franke_2015_oxycodone_acetaminophen"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot   = list(analyte = "oxycodone", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oxycodone", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline; power scaling referenced to the 73.45 kg cohort median named in the Franke 2015 Results text and in the printed CL/F and V/F equations. Table 1 gives median 73.5 kg (range 50.0-120.9). Exponent 0.75 on CL/F and 1 on V/F, both structural rather than estimated.",
      source_name        = "WT"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White reference)",
      notes              = "38 of 151 participants (25.2%, Table 1). Enters the exponent of the multiplicative race factor as RACE_BLACK + 2 * RACE_ASIAN, reproducing the authors' theta^RACE form with the White = 0 / Black = 1 / Asian = 2 dataset coding of the Figure 4 note.",
      source_name        = "RACE"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White reference)",
      notes              = "Only 1 of 151 participants (0.6%, Table 1). Because the authors parameterised race as theta^RACE over a 0/1/2 integer code (Figure 4 note: White = 0, Black = 1, Asian = 2), an Asian participant receives the SQUARE of the Black multiplier. That is an artifact of the power parameterisation over an ordinal code, not a fitted Asian effect; the paper's Discussion states explicitly that PK variability in Asian participants cannot be ruled out from these data.",
      source_name        = "RACE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator, 1 = female",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened in the stepwise forward-addition / backward-elimination covariate search but not retained in the final oxycodone model; Table 3 shows only weight and race. No point estimate is reported, so no effect can be implemented."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained (Table 3). Cohort mean 29.1 years (SD 8.6); enrolment was restricted to 18-55 years."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained in the pooled FASTED oxycodone model (Table 3). Height WAS retained on CL/F in the separate 31-subject food-effect oxycodone model; see modellib('Franke_2015_oxycodone_food')."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened but not retained (Table 3). Enrolment was restricted to 19 to <33 kg/m^2."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 151,
    n_studies      = 4,
    age_range      = "18-55 years (protocol inclusion range)",
    age_median     = "29.1 years (mean, SD 8.6)",
    weight_range   = "50.0-120.9 kg",
    weight_median  = "73.5 kg",
    height_median  = "172.7 cm (range 153.6-194.0)",
    bmi_median     = "25.1 kg/m^2 (range 18.6-32.9)",
    sex_female_pct = 30.5,
    race_ethnicity = c(White = 74.2, Black = 25.2, Asian = 0.6),
    disease_state  = "Healthy adults (studies 1, 2 and 4) and healthy nondependent recreational users of prescription opioids (study 3). Protocol inclusion required BMI 19 to <33 kg/m^2 and body weight >=59 kg.",
    dose_range     = "Single oral doses of one, two or four intact IR/ER OC/APAP 7.5/325 mg tablets, i.e. 7.5, 15 or 30 mg oxycodone, all under fasted conditions.",
    regions        = "United States",
    notes          = "Baseline demographics from Franke 2015 Table 1 (completers, n = 151, of 251 enrolled). Pooled from four phase 1 studies (protocol numbers COV15000170, COV15000172, COV15000255, COV15000244); only the fasted arms of study 4 contributed to this pooled analysis. Table 1 reports a minimum weight of 50.0 kg although the protocols specified a >=59 kg inclusion criterion."
  )

  ini({
    # Structural parameters - typical values at WT = 73.45 kg, White race (RACE = 0)
    lcl <- log(92.4); label("Apparent clearance CL/F at 73.45 kg, White (L/h)")  # Table 3 final model, 'CL/F (L/hour)'; same value in the printed CL/F equation
    lvc <- log(772);  label("Apparent volume of distribution V/F at 73.45 kg, White (L)")  # Table 3 final model, 'V/F (L)'; same value in the printed V/F equation
    lka <- log(1.15); label("First-order absorption rate constant (1/h)")  # Table 3 final model, 'Ka (hour-1)'; estimated (only the acetaminophen Ka was fixed)

    # Body-weight power exponent. Printed as '**0.75' inside the CL/F equation
    # with no standard error or confidence interval, i.e. held at the canonical
    # allometric value rather than estimated. V/F carries an exponent of 1,
    # written in the equation as a bare (WT/73.45) ratio and therefore coded
    # directly in model() rather than as a parameter.
    e_wt_cl <- fixed(0.75); label("Body-weight power exponent on CL/F (unitless)")  # printed CL/F equation, '((WT/73.45)**0.75)'

    # Race effect. The authors wrote theta^RACE over the 0/1/2 integer race code
    # of the Figure 4 note (White = 0, Black = 1, Asian = 2), so these are the
    # per-unit-RACE-code multipliers relative to the White reference.
    e_race_black_cl <- 0.831; label("Multiplicative race factor per unit RACE code on CL/F (unitless)")  # Table 3 final model, 'CL/F ~ race'; same value in the printed CL/F equation
    e_race_black_vc <- 0.827; label("Multiplicative race factor per unit RACE code on V/F (unitless)")  # Table 3 final model, 'V/F ~ race'; same value in the printed V/F equation

    # IIV. Table 3 prints the variance and, in parentheses, 100*sqrt(variance),
    # which the table header labels 'CV %'. The printed pairs self-pin the
    # scale: sqrt(0.0673) = 0.259, sqrt(0.0694) = 0.263 and sqrt(0.631) = 0.794
    # reproduce the printed 25.9, 26.3 and 79.4 exactly, so the first number of
    # each pair is a log-scale VARIANCE and is used here directly.
    etalcl ~ 0.0673  # Table 3 final model, 'interindividual variation' row 'CL/F' = 0.0673 (25.9)
    etalvc ~ 0.0694  # Table 3 final model, 'interindividual variation' row 'V/F' = 0.0694 (26.3)
    etalka ~ 0.631   # Table 3 final model, 'interindividual variation' row 'Ka' = 0.631 (79.4)

    # Residual error. Table 3 prints variance (SD or CV%) in the same style:
    # sqrt(6.94) = 2.63 and sqrt(0.0156) = 0.125 reproduce the printed
    # parenthetical 2.63 and 12.5, so the SDs below are the printed values.
    addSd  <- 2.63;  label("Additive residual error (ng/mL)")  # Table 3 final model, 'intraindividual variation' row 'additive (SD)' = 6.94 (2.63)
    propSd <- 0.125; label("Proportional residual error (fraction)")  # Table 3 final model, 'intraindividual variation' row 'proportional (CV %)' = 0.0156 (12.5)
  })

  model({
    # 1. Derived covariate terms
    # The authors' race term is theta^RACE with RACE coded White = 0,
    # Black = 1, Asian = 2 (Figure 4 note). Rebuilding that integer code from
    # the canonical binary indicators keeps the printed power form exact.
    raceCode <- RACE_BLACK + 2 * RACE_ASIAN

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * (WT / 73.45)^e_wt_cl * e_race_black_cl^raceCode
    vc <- exp(lvc + etalvc) * (WT / 73.45) * e_race_black_vc^raceCode
    ka <- exp(lka + etalka)

    # 3. Micro-constants
    kel <- cl / vc

    # 4. ODE system
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Observation and error
    # central is in mg and vc in L, so central/vc is in mg/L = ug/mL;
    # multiply by 1000 to express Cc in ng/mL, the assay unit in which the
    # paper reports its additive residual error.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
