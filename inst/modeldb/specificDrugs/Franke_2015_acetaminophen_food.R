Franke_2015_acetaminophen_food <- function() {
  description <- "Two-compartment population PK model with first-order absorption for acetaminophen after a single oral dose of two biphasic immediate-release/extended-release oxycodone/acetaminophen (IR/ER OC/APAP) 7.5/325 mg tablets, fitted to the 31-subject FED-vs-FASTED food-effect study (Franke 2015 Table 5). This is a separate fit from the 151-subject pooled fasted analysis in modellib('Franke_2015_acetaminophen') and has a different covariate model: apparent clearance and apparent central volume both carry a sex factor and scale with body weight referenced to the 72.3 kg cohort median (power 0.75 on CL/F, power 1 on V2/F), while the absorption rate constant carries a meal-type factor. A low-fat or high-fat meal cuts the absorption rate constant by 20% and 18% respectively. IMPORTANT: Table 5 reports no intercompartmental clearance and no peripheral volume for this two-compartment fit, in neither the base nor the final model column, so Q and V3/F are carried in as fixed() from the SAME paper's pooled fasted acetaminophen fit (Table 4) by decision of the maintainers; they are not food-effect estimates, and the two fits differ 2.4-fold in central volume. Residual error was not reported for this analysis and is packaged as fixed(0). See the vignette Errata for both gaps and for the sex-effect magnitude, which the printed equations make far larger than the paper's prose describes. Concentrations are in ng/mL, so model() scales amount/volume by 1000 to convert mg/L to ng/mL."
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
      notes = "Baseline; power scaling referenced to the 72.3 kg median weight of the 31-subject food-effect cohort (Table 2, range 59.1-101.8). Exponent 0.75 on CL/F, written '((WT/72.3)**0.75)' in the printed equation, and exponent 1 on V2/F, written as a bare (WT/72.3) ratio. Q and V3/F carry no weight term.",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Sex indicator, 1 = female",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "The authors write the sex term as a power, 'theta9**SEX' on CL/F and 'theta8**SEX' on V2/F, with the legend beneath the equations stating 'men =0; women =1'. Because SEX is binary the power collapses to a multiplier, which is how the authors themselves fill the equations in: '(1 or 0.278)' and '(1 or 0.295)'. Sex was retained for acetaminophen only in this food-effect fit; it was screened and rejected in the 151-subject pooled fasted acetaminophen model, and rejected for oxycodone in both analyses.",
      source_name = "SEX"
    ),
    FED_LOWFAT = list(
      description = "Low-fat, low-calorie meal at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted; FED_LOWFAT and FED_HIGHFAT are mutually exclusive and both 0 under fasting)",
      notes = "Protocol meal: approximately 25%-30% of kilocalories from fat, total 800 +/- 80 kcal (Methods 'Treatments'). Per dose record, not per subject: study 4 was a three-period six-sequence crossover in which each participant was dosed fasted, after a low-fat meal, and after a high-fat meal. Corresponds to the authors' KADIET = theta6 level. Acetaminophen, unlike oxycodone in the same analysis, carries no meal term on its volume.",
      source_name = "DIET"
    ),
    FED_HIGHFAT = list(
      description = "High-fat, high-calorie meal at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted; FED_LOWFAT and FED_HIGHFAT are mutually exclusive and both 0 under fasting)",
      notes = "Protocol meal: approximately 50% of kilocalories from fat, total 1,000 +/- 100 kcal (Methods 'Treatments'), consistent with the FDA high-fat high-calorie definition. Per dose record, as for FED_LOWFAT. Corresponds to the authors' KADIET = theta7 level.",
      source_name = "DIET"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Screened but not retained for acetaminophen in the food-effect analysis (Table 5). Height WAS retained on CL/F for oxycodone in the same 31-subject analysis. Cohort median 168.5 cm (range 155.1-185.5, Table 2)."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened but not retained (Table 5). 8 of 31 participants were Black (25.8%, Table 2)."
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
    dose_range = "A single oral dose of two intact IR/ER OC/APAP 7.5/325 mg tablets (650 mg acetaminophen) on each of three occasions: fasted, after a low-fat meal, and after a high-fat meal.",
    regions = "United States",
    notes = "Baseline demographics from Franke 2015 Table 2 (completers, n = 31, of 48 enrolled). Protocol number COV15000244. Same 31 participants and same dosing occasions as modellib('Franke_2015_oxycodone_food'); the two analytes were fitted separately and with different structural and covariate models. Unlike the pooled fasted analysis, the food-effect data were analysed with first-order conditional estimation WITH INTERACTION."
  )

  ini({
    # Structural parameters - typical values at WT = 72.3 kg, male
    # (SEXF = 0), fasted (FED_LOWFAT = FED_HIGHFAT = 0)
    lcl <- log(22.9); label("Apparent clearance CL/F at 72.3 kg, male (L/h)")  # Table 5 final model, APAP row 'CL/F (L/hour)'; same value in the printed CL/F equation
    lvc <- log(140); label("Apparent central volume of distribution V2/F at 72.3 kg, male (L)")  # Table 5 final model, APAP row 'V2/F (L)'; same value in the printed V2/F equation
    lka <- log(3.17); label("First-order absorption rate constant, fasted (1/h)")  # Table 5 final model, APAP row 'Ka (hour-1)'; same value in the printed Ka equation

    # Disposition parameters NOT reported for this fit. Table 5 footnote a
    # states 'a two-compartment additive and proportional error model was used
    # for APAP', but the APAP block of Table 5 has no Q row and no V3/F row in
    # EITHER the base model column or the final model column, and the printed
    # equations give only CL/F, V2/F and Ka. A two-compartment system cannot be
    # solved without them. By decision of the maintainers (2026-09-21), both
    # are carried in from the SAME paper's pooled fasted acetaminophen fit,
    # Table 4, and held fixed:
    # they are NOT food-effect estimates. The two fits differ 2.4-fold in
    # central volume (58.9 L fasted vs 140 L here), so the borrowed peripheral
    # compartment is not quantitatively coherent with the rest of this model;
    # the early profile that the food-effect analysis is about (Cmax, tmax,
    # AUC) is governed entirely by the paper-sourced CL/F, V2/F and Ka.
    lq <- fixed(log(29.9)); label("Apparent intercompartmental clearance Q (L/h); carried from Table 4, not a food-effect estimate")  # Table 4 final model, 'Q (L/hour)' = 29.9 -- borrowed within paper, see vignette Errata
    lvp <- fixed(log(99.7)); label("Apparent peripheral volume of distribution V3/F (L); carried from Table 4, not a food-effect estimate")  # Table 4 final model, 'V3/F (L)' = 99.7 -- borrowed within paper, see vignette Errata

    # Body-weight power exponent, printed as '**0.75' inside the CL/F equation
    # with no standard error or confidence interval. V2/F carries an exponent
    # of 1, written as a bare (WT/72.3) ratio and so coded directly in model().
    e_wt_cl <- fixed(0.75); label("Body-weight power exponent on CL/F (unitless)")  # printed CL/F equation, '((WT/72.3)**0.75)'

    # Sex factors. The authors write 'theta9**SEX' and 'theta8**SEX' with
    # men = 0 and women = 1, then fill the equations in as '(1 or 0.278)' and
    # '(1 or 0.295)', so for a binary covariate the power is a multiplier
    # applied to women. Table 5 labels these rows '(L/hour)' and '(L)', which
    # cannot be right for a dimensionless multiplier -- the same units-carried-
    # down artefact that affects the Ka rows below.
    e_sexf_cl <- 0.278; label("Female multiplicative factor on CL/F (unitless)")  # Table 5 final model, APAP row 'CL/F sex'; same value in the printed CL/F equation
    e_sexf_vc <- 0.295; label("Female multiplicative factor on V2/F (unitless)")  # Table 5 final model, APAP row 'V/F sex'; same value in the printed V2/F equation

    # Meal-type factors on Ka. The authors define KADIET = 1 for fasting,
    # theta6 for low fat, theta7 for high fat (legend beneath the food-effect
    # equations). Read as multipliers these reproduce the Results text exactly:
    # 1 - 0.802 = 19.8% and 1 - 0.825 = 17.5%, against the stated 'the Ka of OC
    # and APAP decreased 39% and 20% ... low-fat' and '48% and 18% ...
    # high-fat'. Table 5's '(hour-1)' label on these rows is the same units-
    # carried-down artefact noted above.
    e_fed_lowfat_ka <- 0.802; label("Low-fat-meal multiplicative factor on Ka (unitless)")  # Table 5 final model, APAP row 'Ka low fat'; same value in the printed Ka equation
    e_fed_highfat_ka <- 0.825; label("High-fat-meal multiplicative factor on Ka (unitless)")  # Table 5 final model, APAP row 'Ka high fat'; same value in the printed Ka equation

    # IIV. Table 5 reports no variability rows, but the authors' filled-in
    # equations substitute the log-scale variance into the exp() that carries
    # each eta -- the general forms print 'exp(eta1)', 'exp(eta2)', 'exp(eta5)'
    # and the filled forms print the numbers below. Table 4 pins that this
    # substituted number is a VARIANCE and not an SD: its printed pair
    # 0.0554 (23.5) satisfies 100*sqrt(0.0554) = 23.5.
    # No variability is reported for Q or V3/F in this analysis, so neither
    # carries an eta; the Table 4 variances are NOT carried over with the point
    # estimates, because inventing variability is worse than omitting it.
    etalcl ~ 0.0302  # printed CL/F equation for the food-effect APAP model, 'exp(0.0302)'
    etalvc ~ 0.0114  # printed V2/F equation for the food-effect APAP model, 'exp(0.0114)'
    etalka ~ 0.526   # printed Ka equation for the food-effect APAP model, 'exp(0.526)'

    # Residual error. Table 5 footnote a states that a two-compartment additive
    # and proportional error model was used for APAP, but NO residual magnitude
    # is reported anywhere for the food-effect analysis - not in the base model
    # column, not in the final model column, and not in the equations. The
    # structure is therefore declared with both magnitudes held at zero rather
    # than carrying over the unrelated 151-subject pooled-fasted estimates.
    addSd <- fixed(0); label("Additive residual error (ng/mL); magnitude not reported for the food-effect analysis")
    propSd <- fixed(0); label("Proportional residual error (fraction); magnitude not reported for the food-effect analysis")
  })

  model({
    # 1. Derived covariate terms
    # Mutually exclusive meal indicators; both 0 selects the fasted reference,
    # reproducing the authors' three-level KADIET factor. Acetaminophen carries
    # a meal term on Ka only, unlike oxycodone which also carries one on V/F.
    kadiet <- 1 + (e_fed_lowfat_ka - 1) * FED_LOWFAT + (e_fed_highfat_ka - 1) * FED_HIGHFAT

    # The authors' 'theta**SEX' power form, written directly. SEXF is 0 for men
    # and 1 for women, so this is 1 for men and the estimate for women.
    sexcl <- e_sexf_cl^SEXF
    sexvc <- e_sexf_vc^SEXF

    # 2. Individual parameters
    cl <- exp(lcl + etalcl) * sexcl * (WT / 72.3)^e_wt_cl
    vc <- exp(lvc + etalvc) * sexvc * (WT / 72.3)
    ka <- exp(lka + etalka) * kadiet
    q <- exp(lq)
    vp <- exp(lvp)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation and error
    # central is in mg and vc in L, so central/vc is in mg/L = ug/mL;
    # multiply by 1000 to express Cc in ng/mL, the paper's assay unit.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
