Ishibashi_2024_ensitrelvir <- function() {
  description <- "Two-compartment population PK model with first-order absorption for ensitrelvir in healthy adults and participants with SARS-CoV-2 infection (Ishibashi 2024)"
  reference <- "Ishibashi T, Shimizu R, Kubota R. Population Pharmacokinetics of Ensitrelvir in Healthy Participants and Participants with SARS-CoV-2 Infection in the SCORPIO-SR Study. Clin Pharmacokinet. 2024;63(12):1723-1734. doi:10.1007/s40262-024-01446-4"
  vignette <- "Ishibashi_2024_ensitrelvir"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Ensitrelvir is dosed orally as a suspension or tablet
  # and measured in plasma (Sect. 2.1); all disposition parameters are
  # apparent (CL/F, Vc/F, Q/F, Vp/F), so the depot carries the administered
  # amount and F is absorbed into the apparent volumes and clearances.
  compartmentData <- list(
    depot       = list(analyte = "ensitrelvir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ensitrelvir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ensitrelvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power (allometric-style) covariate model on both CL/F and Vc/F, normalized to the 62.6 kg population median body weight of Ishibashi 2024 Table 2. Exponents were estimated (0.521 on CL/F, 1.04 on Vc/F), not fixed to the 0.75/1 allometric defaults. Observed range 35.0-156.0 kg; the paper reports the resulting Cmax and AUC ratios over that range as 0.410-1.741 and 0.623-1.353 relative to the median.",
      source_name        = "BW"
    ),
    FED = list(
      description        = "Fed-state indicator at the time of dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted / any other prandial state)",
      notes              = "Ishibashi 2024 defines this covariate operationally as 'Food: 1 = administration within 2 h after a meal, 0 = other' (Table 3 footnote block), i.e. a generic fed-vs-fasted flag keyed to a 2 h post-meal dosing window rather than a standardized meal-composition challenge. The paper describes the estimate as being 'for fed to fasted' (Sect. 4). The dedicated phase I food-effect arms used a high-fat/high-calorie meal, but the covariate is applied as a general fed flag across all three pooled studies, so the general FED canonical applies rather than FED_HIGHFAT (same reasoning as Chen 2023 nemonoxacin and Comisar 2025 rimegepant). Per dose record: the food-effect arms were two-way crossovers, so a participant contributes both levels.",
      source_name        = "Food"
    ),
    FORM_TABLET = list(
      description        = "Tablet-formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral suspension)",
      notes              = "Ishibashi 2024 Table 3 footnote block: 'Formulation: 1 = tablet, 0 = suspension'. The non-tablet comparator for this model is an oral suspension. Suspension was used only in part of the phase I programme (62 of 2060 participants, 3.0%, Table 2); the phase II/III SCORPIO-SR study and all clinical use employ the tablet, so FORM_TABLET = 1 is the clinically relevant setting.",
      source_name        = "Formulation"
    )
  )

  # Covariates that Ishibashi 2024 screened but did not retain in the final
  # model (no. 511). Documentation only -- these are deliberately NOT
  # referenced in model(). Two distinct exclusion mechanisms are recorded in
  # the paper and noted per entry below:
  #   (a) collinearity pruning before forward selection -- covariate pairs
  #       with |r| >= 0.7 were reduced to the member with the larger drop in
  #       the objective function value (Sect. 3.1);
  #   (b) inferential assessment after backward deletion -- backward deletion
  #       removed nothing, but effects whose 90% CI for the parameter ratio
  #       fell entirely inside the 0.80-1.25 bioequivalence range were judged
  #       not clinically meaningful and were dropped step by step (Sect. 3.1,
  #       Figs. S1-S3).
  # Country of enrolment (Japan / Korea / Vietnam) was also tested on CL/F
  # and Vc/F and removed by mechanism (b); it has no canonical register entry
  # and is therefore described here rather than given an entry below.
  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened on CL/F and Vc/F. Excluded by mechanism (a): correlated with body weight at |r| >= 0.7, and body weight gave the larger drop in objective function value."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on CL/F and Vc/F; not selected. Observed range 12-76 years, median 35."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F and Vc/F; entered the full model (no. 302) on Vc/F but was removed by mechanism (b)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F and Vc/F; not selected. 98.5% of the analysis population was Asian (Table 2), which the paper names as its principal limitation; a race-stratified prediction-corrected VPC (Fig. 3) was used instead of a race covariate."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Reported in g/dL by Ishibashi 2024 Table 2 (mean 4.4, range 0.5-5.8); note the register canonical is g/L. Screened on CL/F and Vc/F; entered the full model (no. 302) on Vc/F but was removed by mechanism (b)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL/F. Excluded by mechanism (a): correlated with ALT at |r| >= 0.7, and ALT gave the larger drop in objective function value."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened on CL/F and carried into forward selection as the retained member of the ALT/AST correlated pair; not retained in the final model."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Reported in mg/dL by Ishibashi 2024 Table 2 (mean 0.5, range 0.1-2.1); note the register canonical is umol/L. Screened on CL/F; not selected."
    ),
    CRCL = list(
      description = "Renal function (creatinine clearance and estimated glomerular filtration rate)",
      units       = "mL/min (CrCL, eGFRabs) or mL/min/1.73 m^2 (eGFR)",
      type        = "continuous",
      notes       = "Ishibashi 2024 screened four renal-function covariates on CL/F: creatinine clearance, serum creatinine, eGFR and absolute eGFR. Creatinine clearance was carried into forward selection as the retained member of the correlated renal cluster by mechanism (a) (CrCL was correlated with both BW and eGFRabs, and eGFR with eGFRabs, at |r| >= 0.7); none was retained in the final model. Ensitrelvir is only partly renally eliminated -- 12.9-21.8% of a dose is recovered in urine as unchanged drug (Sect. 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Reported in mg/dL by Ishibashi 2024 Table 2 (mean 0.76, range 0.37-1.43). Screened on CL/F as part of the renal-function cluster; not selected."
    ),
    DIS_HEALTHY = list(
      description = "Health status: healthy participant vs participant infected with SARS-CoV-2",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F, Vc/F and relative bioavailability F1; entered the full model (no. 302) on CL/F and F1 but was removed by mechanism (b). Its exclusion is the paper's justification for pooling the 175 healthy phase I participants with the 1885 infected SCORPIO-SR participants under one set of typical values."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2060,
    n_studies      = 3,
    age_range      = "12-76 years",
    age_median     = "35 years",
    weight_range   = "35.0-156.0 kg",
    weight_median  = "62.6 kg",
    sex_female_pct = 42.6,
    race_ethnicity = c(Asian = 98.5, White = 1.2, Other = 0.2),
    disease_state  = "healthy adults (175 participants, 8.5%) pooled with participants with mild-to-moderate SARS-CoV-2 infection (1885 participants, 91.5%)",
    dose_range     = "phase I: 20, 70, 250, 500, 1000 or 2000 mg single oral dose, or 375 mg day 1 then 125 mg once daily days 2-5, or 750 mg day 1 then 250 mg once daily days 2-5; phase II/III (SCORPIO-SR): 375/125 mg or 750/250 mg",
    regions        = "Japan (72.0%), Vietnam (21.7%), Korea (6.3%)",
    notes          = "Baseline demographics from Ishibashi 2024 Table 2; study inventory from Table 1. 8034 plasma ensitrelvir concentrations were pooled from two phase I studies in healthy participants (including food-effect and drug-drug-interaction assessments) and the phase II/III SCORPIO-SR study, contributing 4341 and 3693 concentrations respectively. Adolescents aged 12 to <18 years made up 1.7% of the population. Formulation split was 3.0% suspension and 97.0% tablet."
  )

  ini({
    # Structural parameters. The typical values below are the reference
    # condition: a 62.6 kg participant (the Table 2 population median body
    # weight) dosed with the SUSPENSION in the FASTED state, i.e.
    # FED = 0 and FORM_TABLET = 0.
    lka <- log(1.50);  label("Absorption rate constant, suspension/fasted reference (1/h)")   # Table 3 'Ka (1/h)' = 1.50 (95% CI 1.25-1.75, %RSE 8.6); final-model equation in the Sect. 4 legend block: 'Ka = 1.50 x (0.594 for food) x (0.362 for formulation)'
    lcl <- log(0.211); label("Apparent clearance CL/F at 62.6 kg (L/h)")                      # Table 3 'CL/F (L/h)' = 0.211 (95% CI 0.208-0.214, %RSE 0.8); Sect. 4 legend block: 'CL/F = 0.211 x (body weight/62.6)^0.521'
    lvc <- log(14.7);  label("Apparent central volume Vc/F at 62.6 kg (L)")                   # Table 3 'Vc/F (L)' = 14.7 (95% CI 13.9-15.5, %RSE 2.7); Sect. 4 legend block: 'Vc/F = 14.7 x (body weight/62.6)^1.04'
    lq  <- log(0.539); label("Apparent intercompartmental clearance Q/F (L/h)")               # Table 3 'Q/F (L/h)' = 0.539 (95% CI 0.321-0.757, %RSE 20.6); Sect. 4 legend block: 'Q/F = 0.539' (no covariates)
    lvp <- log(2.50);  label("Apparent peripheral volume Vp/F (L)")                           # Table 3 'Vp/F (L)' = 2.50 (95% CI 1.80-3.20, %RSE 14.3); Sect. 4 legend block: 'Vp/F = 2.50' (no covariates)

    # Covariate effects. Categorical covariates enter multiplicatively
    # (Sect. 2.2.2), so the estimate is the ratio applied when the indicator
    # is 1 and the power form X^IND reproduces the printed equation exactly.
    e_fed_ka         <- 0.594; label("Fed-state multiplicative factor on ka (power-form base, unitless)")            # Table 3 'Effect of food on Ka' = 0.594 (95% CI 0.398-0.790, %RSE 16.8); a 41% slower absorption when dosed within 2 h after a meal
    e_form_tablet_ka <- 0.362; label("Tablet-formulation multiplicative factor on ka (power-form base, unitless)")   # Table 3 'Effect of formulation on Ka' = 0.362 (95% CI 0.262-0.462, %RSE 14.1); a 64% slower absorption from the tablet than from the suspension
    e_wt_cl          <- 0.521; label("Body-weight power exponent on CL/F (unitless)")                                # Table 3 'Effect of body weight on CL/F' = 0.521 (95% CI 0.456-0.586, %RSE 6.4); estimated, not fixed at the 0.75 allometric default
    e_wt_vc          <- 1.04;  label("Body-weight power exponent on Vc/F (unitless)")                                # Table 3 'Effect of body weight on Vc/F' = 1.04 (95% CI 0.960-1.12, %RSE 3.9); estimated, not fixed at the 1 allometric default

    # Interindividual variability, log-normal / exponential (Sect. 2.2.1).
    # Table 3 reports the IIV rows as percentages. Those percentages are the
    # omega STANDARD DEVIATIONS on the log scale, NOT %CV -- the covariance
    # row settles it. Table 3 gives cov(CL/F, Vc/F) = 0.0216 with R = 0.691;
    # reading the percentages as raw SDs gives 0.691 * 0.213 * 0.147 =
    # 0.021636, which rounds to the printed 0.0216, whereas reading them as
    # %CV (omega = sqrt(log(1 + CV^2))) gives 0.021282, which rounds to
    # 0.0213 and is outside the printed value's rounding interval.
    # Variances below are therefore the squared percentages:
    #   CL/F 0.213^2 = 0.045369, Vc/F 0.147^2 = 0.021609, Ka 0.729^2 = 0.531441.
    # IIV on Q/F was removed at model no. 510 (omega^2 = 1.03 with 84.3%
    # shrinkage) and IIV on Vp/F at model no. 301 (estimate almost 0, complete
    # shrinkage), so neither carries an eta here.
    etalcl + etalvc ~ c(0.045369,
                        0.0216, 0.021609)                                                    # Table 3 IIV 'CL/F (%)' = 21.3, 'Vc/F (%)' = 14.7, 'Covariance between CL/F and Vc/F' = 0.0216 (R = 0.691); covariance added at model no. 509
    etalka ~ 0.531441                                                                        # Table 3 IIV 'Ka (%)' = 72.9 (95% CI 64.3-80.6, %RSE 11.3); shrinkage 68.4%, retained because the large omega with small %RSE was judged to contribute to describing absorption (Sect. 4)

    # Residual error: combination (additive + proportional) model, selected
    # over additive-only and proportional-only on objective function value
    # (Sect. 3.1). Shrinkage in epsilon was 10.1%.
    addSd  <- 0.0317; label("Additive residual error (ug/mL)")            # Table 3 'Additive residual error (ug/mL)' = 0.0317 (95% CI 0.0224-0.0410, %RSE 14.9)
    propSd <- 0.199;  label("Proportional residual error (fraction)")     # Table 3 'Proportional residual error (%)' = 19.9 (95% CI 18.5-21.3, %RSE 3.6)
  })

  model({
    # Absorption. Categorical covariates enter multiplicatively (Sect. 2.2.2);
    # the final-model equation printed in the Sect. 4 legend block is
    #   Ka = 1.50 x (0.594 for food) x (0.362 for formulation)
    # with Food = 1 for administration within 2 h after a meal and
    # Formulation = 1 for tablet. The X^IND power form reproduces that
    # exactly: each factor is applied when its indicator is 1 and is unity
    # otherwise. The tablet fasted value is 1.50 * 0.362 = 0.543 1/h, which
    # is the setting used throughout SCORPIO-SR and in clinical use.
    ka <- exp(lka + etalka) * e_fed_ka^FED * e_form_tablet_ka^FORM_TABLET

    # Disposition. Continuous covariates enter as a power model (Sect. 2.2.2)
    # normalized to the 62.6 kg Table 2 median body weight; body weight on
    # CL/F and on Vc/F are the only continuous covariates in the final model.
    cl <- exp(lcl + etalcl) * (WT / 62.6)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 62.6)^e_wt_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # All disposition parameters are apparent (CL/F, Vc/F, Q/F, Vp/F), so
    # bioavailability is absorbed into them and there is no f(depot) hook.
    # Relative bioavailability F1 was tested as a covariate target for food,
    # formulation and health status, but every effect was removed by
    # inferential assessment (Sect. 3.1), leaving F1 with no retained effect.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
