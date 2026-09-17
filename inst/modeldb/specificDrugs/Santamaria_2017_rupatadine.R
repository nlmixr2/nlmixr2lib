Santamaria_2017_rupatadine <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for oral rupatadine 1 mg/mL solution in 6-11 year old children with allergic rhinitis, with a linear-additive body-weight effect on apparent clearance (Santamaria 2017)"
  reference <- "Santamaria E, Estevez JA, Riba J, Izquierdo I, Valle M. Population pharmacokinetic modelling of rupatadine solution in 6-11 year olds and optimisation of the experimental design in younger children. PLoS ONE. 2017;12(4):e0176091. doi:10.1371/journal.pone.0176091"
  vignette <- "Santamaria_2017_rupatadine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Santamaria 2017 Material and methods
  # ("Population", "Analytical method"): rupatadine was given as an oral
  # solution and measured in plasma by LC-MS/MS.
  compartmentData <- list(
    depot = list(analyte = "rupatadine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rupatadine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rupatadine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "The only covariate retained in the final model. Effect form is linear-additive on apparent clearance and NOT allometric: the paper states 'Although an allometric scaling model was tested, the function that best described the weight-CL/F relationship was a linear model' and prints CL = theta1 + theta2 * WEIGHT/38.5 (Santamaria 2017 Results 'Covariate inclusion' and the Table 3 footnote). Weight is normalised to 38.5 kg, the cohort median body weight (Table 2), so theta2 is the clearance increment contributed at the median weight. Cohort range 22.0-68.5 kg (Table 2); the model was subsequently used to simulate 10-25 kg children aged 2-5 years, which is below the range it was fitted over. Adding weight decreased the objective function by 3 points (P = 0.08) and explained 11.1 percent of the interindividual variability in CL/F. Weight was also evaluated on V/F, where the results argued against its inclusion.",
      source_name = "WEIGHT"
    )
  )

  # Covariates screened by the authors but NOT retained in the final model.
  # Documentation only: these are deliberately absent from model(). Santamaria
  # 2017 Material and methods ("Population analysis") lists the investigated
  # covariates as age, sex, weight, height and BMI; Results ("Covariate
  # inclusion") reports the disposition of each.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened and tested in NONMEM but excluded. On CL/F 'the estimation error increased considerably'; on V/F 'there was no change in the OF and it was deemed appropriate to exclude age from the final model following the criterion of simplicity'. The Discussion attributes the absent age effect to the narrow enrolled age range and to CYP3A4 reaching maturation at about 1.3 years of age."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 = male",
      notes = "Screened graphically but not carried into NONMEM: 'Therefore sex and height were not further evaluated in NONMEM'. The Discussion notes rupatadine clearance does depend on sex in adults, and attributes the null result here to the small sample and to weight capturing male-female differences."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened but not carried into NONMEM because of its clear correlation with weight; weight was retained as the size descriptor with more promising clinical application. Source column is in metres (Table 2 reports 1.18-1.59 m); the canonical HT column is in cm."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested for significance in NONMEM but 'no clear improvement of the model was detected respect to the base model (OFV = -114 vs -112, P > 0.05)'. Cohort range 13.4-27.1 kg/m^2 (Table 2)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    n_observations = 84L,
    age_range = "7.94-11.93 years",
    age_median = "10.41 years",
    weight_range = "22.0-68.5 kg",
    weight_median = "38.5 kg",
    height_range = "1.18-1.59 m",
    height_median = "1.44 m",
    bmi_range = "13.4-27.1 kg/m^2",
    bmi_median = "18.54 kg/m^2",
    sex_female_pct = 54.5,
    race_ethnicity = "Not reported in source",
    disease_state = "Allergic rhinitis, otherwise in good health",
    dose_range = "Single oral dose of rupatadine 1 mg/mL solution; 2.5 mg (2.5 mL) for children weighing more than 10 and less than 25 kg, 5 mg (5 mL) for children weighing 25 kg or more. Nine of the eleven children received 5 mg.",
    regions = "Australia (Royal Children's Hospital, Melbourne; Peninsula Private Hospital and Peninsula Clinical Research Centre, Rivercity)",
    co_medication = "Participants who had taken any medication that could significantly interact with CYP3A4 were excluded, as rupatadine is mainly metabolised by that enzyme.",
    notes = "Open-label single-dose study in children 6-11 years of age weighing at least 16 kg with a history of allergic rhinitis (Santamaria 2017 Material and methods, 'Population'; demographics in Table 2). A full concentration-time profile was obtained for each child, with 8 samples per child at predose and 0.5, 1, 2, 4, 8, 12 and 24 h postdose. 84 plasma concentrations entered the analysis; concentrations below the limit of quantification in the elimination phase (12 percent of observations) were discarded. Observed rupatadine concentration range 0.1-4.9 ng/mL. Bioanalysis by validated LC-MS/MS with clomipramine as internal standard; within- and between-run precision error below 12.71 percent and accuracy-related errors within plus or minus 12.83 percent. Estimation by FOCE in NONMEM. The same final model was then used to select a 2.5 mg dose and to optimise the sampling design for a planned study in 2-5 year olds; that design optimisation (WinPOPT, paper Tables 1 and 4) is a study-design exercise rather than a model and is not represented here."
  )

  ini({
    # Structural PK parameters, all apparent (secondary to bioavailability F,
    # which rupatadine's large first-pass metabolism makes small). Santamaria
    # 2017 Table 3, "Final Model" column; the parenthesised number in the table
    # is the RSE in percent.
    lka <- log(0.53)
    label("First-order absorption rate constant (1/h)") # Table 3 final model, ka = 0.53 (RSE 15)
    ltlag <- log(0.22)
    label("Absorption lag time (h)") # Table 3 final model, 'Lag time' = 0.22 (RSE 13)
    lvc <- log(108)
    label("Apparent central volume of distribution (L)") # Table 3 final model, Vc/F = 108 (RSE 52)
    lq <- log(209)
    label("Apparent inter-compartmental clearance (L/h)") # Table 3 final model, CLd/F = 209 (RSE 30)
    lvp <- log(1430)
    label("Apparent peripheral volume of distribution (L)") # Table 3 final model, Vp/F = 1430 (RSE 56)

    # Apparent clearance is a LINEAR-ADDITIVE function of body weight, not an
    # allometric power function. The paper prints
    #   CL = theta1 + theta2 * WEIGHT/38.5
    # (Santamaria 2017 Table 3 footnote and the display equation in Results,
    # 'Covariate inclusion'), having explicitly tested and rejected allometric
    # scaling. theta1 is therefore the weight-independent intercept and is
    # carried by the canonical lcl; theta2 is the slope on the weight ratio and
    # is carried by e_wt_cl, which in this model has units of L/h rather than
    # the dimensionless exponent it carries in allometric models. This
    # positive-intercept-plus-covariate-term encoding is the form prescribed by
    # inst/references/parameter-names.md ('Covariate-expression shape tokens',
    # <param>_int entry); the precedent in the registry is
    # Blair_2004_raltitrexed.R (cl <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)).
    lcl <- log(225)
    label("Weight-independent intercept of the apparent clearance regression (L/h)") # Table 3 final model, theta1 = 225 (RSE 63)
    e_wt_cl <- 333
    label("Body-weight slope on apparent clearance, per unit of weight normalised to the cohort median (L/h)") # Table 3 final model, theta2 = 333 (RSE 44)

    # Interindividual variability. The paper's Material and methods
    # ("Population analysis") gives the IIV model as CL_i = CL_pop * exp(eta_i)
    # with eta of mean 0 and variance omega^2, i.e. log-normal. Table 3 reports
    # IIV as a coefficient of variation in percent (table footnote: 'IIV,
    # interpatient variability expressed as coefficient of variation'), which
    # under the usual NONMEM reporting convention is 100 * sqrt(omega^2). The
    # variances below are therefore the squared tabulated fractions. The
    # alternative log-normal conversion omega^2 = log(1 + CV^2) is discussed in
    # the vignette Errata; it cannot be discriminated from the paper's own
    # simulation results, and the convention used here is the one consistent
    # with the paper computing its 11.1 percent explained-variability figure
    # directly from the tabulated numbers as (45 - 40)/45.
    # The data supported IIV on CL/F and Vc/F only; the paper reports that
    # interindividual variability in the absorption process was not estimable.
    # No covariance between the two random effects is reported.
    etalcl ~ 0.16 # Table 3 final model, 'IIV CL/F (%)' = 40 (RSE 25); 0.40^2 = 0.16
    etalvc ~ 0.8798 # Table 3 final model, 'IIV Vc/F (%)' = 93.8 (RSE 38); 0.938^2 = 0.8798

    # Residual error. "Residual variability was explained by means of an
    # additive error model" (Results, 'Population analysis'), and Table 3
    # reports the residual error in ng/mL, the same unit as the observations.
    addSd <- 0.18
    label("Additive residual error (ng/mL)") # Table 3 final model, 'Residual error (ng/mL)' = 0.18 (RSE 41)
  })

  model({
    # Individual parameters. The weight term is added on the linear scale to
    # reproduce the paper's printed regression exactly; the log-normal IIV
    # multiplies the whole typical clearance, matching CL_i = CL_pop * exp(eta_i)
    # where CL_pop is the covariate-adjusted typical value.
    ka <- exp(lka)
    tlag <- exp(ltlag)
    cl <- (exp(lcl) + e_wt_cl * WT / 38.5) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg and vc in L give mg/L = ug/mL; multiply by 1000 to obtain the
    # ng/mL used for the observations and for the additive residual error.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
