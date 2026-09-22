Franke_2015_acetaminophen <- function() {
  description <- "Two-compartment population PK model with first-order absorption for acetaminophen after a single oral dose of biphasic immediate-release/extended-release oxycodone/acetaminophen (IR/ER OC/APAP) 7.5/325 mg tablets under FASTED conditions, pooled post hoc from four phase 1 studies in 151 healthy adults and nondependent recreational users of prescription opioids (Franke 2015 Table 4). Body weight is the only retained covariate, acting on apparent clearance (power 0.75) and on the apparent central volume (power 1) referenced to the 73.45 kg cohort median; intercompartmental clearance and the peripheral volume carry no covariate. The absorption rate constant was FIXED at 5.4 1/h to stabilise the fit. Concentrations are in ng/mL, so model() scales amount/volume by 1000 to convert mg/L to ng/mL. The companion oxycodone model from the same pooled fasted dataset is modellib('Franke_2015_oxycodone')."
  reference <- "Franke RM, Morton T, Devarakonda K. Pooled post hoc analysis of population pharmacokinetics of oxycodone and acetaminophen following a single oral dose of biphasic immediate-release/extended-release oxycodone/acetaminophen tablets. Drug Des Devel Ther. 2015;9:4587-4597. doi:10.2147/DDDT.S79499"
  vignette <- "Franke_2015_oxycodone_acetaminophen"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "acetaminophen", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "acetaminophen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "acetaminophen", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline; power scaling referenced to the 73.45 kg cohort median named in the printed CL/F and V2/F equations (the Results text for acetaminophen quotes the same median rounded to 73.5). Table 1 gives median 73.5 kg (range 50.0-120.9). Exponent 0.75 on CL/F and 1 on V2/F, both structural rather than estimated. Q and V3/F carry no weight term.",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate search but not retained in the pooled FASTED acetaminophen model (Table 4). Sex WAS retained on CL/F and V2/F in the separate 31-subject food-effect acetaminophen model, which is not packaged because Table 5 omits Q and V3/F for that fit."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened but not retained (Table 4). Cohort mean 29.1 years (SD 8.6)."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened but not retained (Table 4)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained (Table 4)."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained for acetaminophen (Table 4), in contrast to oxycodone where race was retained on both CL/F and V/F."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 151,
    n_studies = 4,
    age_range = "18-55 years (protocol inclusion range)",
    age_median = "29.1 years (mean, SD 8.6)",
    weight_range = "50.0-120.9 kg",
    weight_median = "73.5 kg",
    height_median = "172.7 cm (range 153.6-194.0)",
    bmi_median = "25.1 kg/m^2 (range 18.6-32.9)",
    sex_female_pct = 30.5,
    race_ethnicity = c(White = 74.2, Black = 25.2, Asian = 0.6),
    disease_state = "Healthy adults (studies 1, 2 and 4) and healthy nondependent recreational users of prescription opioids (study 3). Protocol inclusion required BMI 19 to <33 kg/m^2 and body weight >=59 kg.",
    dose_range = "Single oral doses of one, two or four intact IR/ER OC/APAP 7.5/325 mg tablets, i.e. 325, 650 or 1300 mg acetaminophen, all under fasted conditions.",
    regions = "United States",
    notes = "Baseline demographics from Franke 2015 Table 1 (completers, n = 151, of 251 enrolled). Same pooled fasted dataset as modellib('Franke_2015_oxycodone'); the two analytes were fitted separately and with different structural models."
  )

  ini({
    # Structural parameters - typical values at WT = 73.45 kg
    lcl <- log(20.5); label("Apparent clearance CL/F at 73.45 kg (L/h)")  # Table 4 final model, 'CL/F (L/hour)'; same value in the printed CL/F equation
    lvc <- log(58.9); label("Apparent central volume of distribution V2/F at 73.45 kg (L)")  # Table 4 final model, 'V2/F (L)'; same value in the printed V2/F equation
    lq  <- log(29.9); label("Apparent intercompartmental clearance Q (L/h)")  # Table 4 final model, 'Q (L/hour)'
    lvp <- log(99.7); label("Apparent peripheral volume of distribution V3/F (L)")  # Table 4 final model, 'V3/F (L)'

    # Ka was held constant, not estimated: Methods 'statistical analyses'
    # states 'To stabilize the resulting models, the Ka was fixed at 5.4
    # hour-1 for the single-dose population PK of APAP.'
    lka <- fixed(log(5.4)); label("First-order absorption rate constant (1/h)")  # Table 4, 'Ka (hour-1) (fixed)' = 5.4

    # Body-weight power exponent, printed as '**0.75' inside the CL/F equation
    # with no standard error or confidence interval. V2/F carries an exponent
    # of 1, written as a bare (WT/73.45) ratio and so coded directly in model().
    e_wt_cl <- fixed(0.75); label("Body-weight power exponent on CL/F (unitless)")  # printed CL/F equation, '((WT/73.45)**0.75)'

    # IIV. Table 4 prints the variance and, in parentheses, 100*sqrt(variance)
    # under a 'CV %' header. The printed pairs self-pin the scale:
    # sqrt(0.0554) = 0.235, sqrt(0.453) = 0.673, sqrt(4.35) = 2.09 and
    # sqrt(0.75) = 0.866 reproduce the printed 23.5, 67.3, 209.0 and 86.3, so
    # the first number of each pair is a log-scale VARIANCE. The CL/F and V2/F
    # variances are taken from the printed equations, which carry one more
    # significant figure than the table's rounded 0.06 and 0.45.
    etalcl ~ 0.0554  # printed CL/F equation 'exp(0.0554)'; Table 4 rounds this to 0.06 (23.5)
    etalvc ~ 0.453   # printed V2/F equation 'exp(0.453)'; Table 4 rounds this to 0.45 (67.3)
    etalq  ~ 4.35    # Table 4 final model, 'interindividual variation' row 'Q' = 4.35 (209.0)
    etalvp ~ 0.75    # Table 4 final model, 'interindividual variation' row 'V3/F' = 0.75 (86.3)

    # Residual error. Same variance (SD or CV%) print style:
    # sqrt(1.59e+3) = 39.9 reproduces the printed 39.9 exactly. The
    # proportional variance is printed to only one significant figure (0.03),
    # so the SD is taken from the more precise parenthetical CV% of 17.1.
    addSd  <- 39.9;  label("Additive residual error (ng/mL)")  # Table 4 final model, 'intraindividual variation' row 'additive (SD)' = 1.59e+3 (39.9)
    propSd <- 0.171; label("Proportional residual error (fraction)")  # Table 4 final model, 'intraindividual variation' row 'proportional (CV %)' = 0.03 (17.1)
  })

  model({
    # 1. Individual parameters
    cl <- exp(lcl + etalcl) * (WT / 73.45)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 73.45)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Observation and error
    # central is in mg and vc in L, so central/vc is in mg/L = ug/mL;
    # multiply by 1000 to express Cc in ng/mL, the assay unit in which the
    # paper reports its additive residual error.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
