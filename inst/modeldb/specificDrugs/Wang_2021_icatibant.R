Wang_2021_icatibant <- function() {
  description <- "Two-compartment population pharmacokinetic model for subcutaneous icatibant, a bradykinin B2 receptor antagonist, in healthy adults and in adult and pediatric patients with hereditary angioedema (HAE) (Wang 2021). First-order absorption with a lag time; estimated allometric body-weight exponents shared by CL/F and Q/F and by Vc/F and Vp/F (70 kg reference); a linear age effect on CL/F centred at 25 years; female sex on CL/F and Vc/F; an acute HAE attack at dosing on CL/F; non-White race on Vc/F; proportional residual error."
  reference <- "Wang Y, Jomphe C, Marier JF, Martin P. Population Pharmacokinetics and Exposure-Response Analyses to Guide Dosing of Icatibant in Pediatric Patients With Hereditary Angioedema. J Clin Pharmacol. 2021;61(4):555-564. doi:10.1002/jcph.1768"
  vignette <- "Wang_2021_icatibant"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling (WT/70) on all four disposition parameters: one estimated exponent (0.516) shared by CL/F and Q/F, and a second (0.671) shared by Vc/F and Vp/F (Table 2; the shared exponents carry identical RSEs, 12.7% and 9.7%). Cohort range 12.3-102 kg (Table 1).",
      source_name = "weight"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear effect on CL/F, 1 + (-0.0107 * (AGE - 25)), applied after the body-weight effect (Table 2; Table S5 footnote a). CL/F is 16% higher at 10 years and about 20% lower at 45 years than at 25 years (Results). Cohort range 3.42-54 years (Table 1).",
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Multiplicative factors 0.882 on CL/F and 0.855 on Vc/F for females (Table 2), modelled as exp(effect) (Table S5 footnote a).",
      source_name = "sex"
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (White; the typical-value reference in this model)",
      notes = "The paper dichotomises race as White vs non-White (Table 2 'x 1.11 if nonwhite' on Vc/F). The non-White group pools Black or African American (20.9% of the cohort) and Other (2.3%) (Table 1). White is the reference, so the effect enters on (1 - RACE_WHITE).",
      source_name = "race (white / nonwhite)"
    ),
    DIS_HAE_ACUTE = list(
      description = "Acute HAE attack at the time of dosing (attack onset within 12 hours before icatibant administration)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no HAE attack at dosing: healthy volunteers, and pediatric patients with HAE dosed outside an attack)",
      notes = "Multiplicative factor 0.911 on CL/F (Table 2 'x 0.911 if HAE attack'), modelled as exp(effect) (Table S5 footnote a). Table S2 footnote a defines the attack as occurring within 12 hours before icatibant administration; by that table the flag is 1 for all 8 adult patients with HAE (study JE049-2101) and 21 of 31 pediatric patients (HGT-FIR-086), and 0 for all 133 healthy adults and the other 10 pediatric patients. The effect of HAE attack on Vc/F was tested in the full model and removed because its bootstrap 95% CI included 0 (Table S5 footnote b).",
      source_name = "HAE attack status"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "icatibant", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "icatibant", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "icatibant", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 172,
    n_studies = 6,
    age_range = "3.42-54.0 years",
    age_median = "25.0 years",
    weight_range = "12.3-102 kg",
    weight_median = "69.5 kg",
    sex_female_pct = 41.3,
    race_ethnicity = c(White = 76.7, Black = 20.9, Other = 2.3),
    disease_state = "133 healthy adults (77.3%), 8 adults with hereditary angioedema (4.7%), and 31 children and adolescents aged 2-17 years with hereditary angioedema (18.0%); 29 subjects were dosed during an acute HAE attack",
    dose_range = "Subcutaneous icatibant: 30 mg single dose, 3 x 30 mg every 6 h, 0.05-0.4 mg/kg single ascending doses, 30 or 45 mg single dose in adult patients, 0.4 mg/kg (capped at 30 mg) single dose in pediatric patients",
    regions = "Not reported",
    notes = "2172 measurable plasma icatibant concentrations (523 BLQ samples excluded; Table S3) pooled from four phase 1 studies in healthy adults (HGT-FIR-061, HGT-FIR-065, JE049-1102, JE049-1103), a phase 2 study in adult patients with HAE (JE049-2101) and a phase 3 pediatric study (HGT-FIR-086) (Table S1). Baseline demographics are in Table 1 and, by study, Table S2. The 90 mg dose arm of HGT-FIR-061 was excluded as supra-therapeutic (Table S1 footnote). NONMEM 7.3."
  )

  ini({
    # Structural parameters -- typical values for a 70 kg, 25-year-old White
    # male without an HAE attack at dosing (Table 2, final model).
    lka <- log(3.27); label("Absorption rate constant ka (1/h)") # Table 2: ka 3.27 1/h (RSE 3.50%)
    ltlag <- log(0.0426); label("Absorption lag time tlag (h)") # Table 2: tlag 0.0426 h (RSE 10.6%)
    lcl <- log(15.4); label("Apparent clearance CL/F (L/h)") # Table 2: Cl/F 15.4 L/h (RSE 2.50%)
    lvc <- log(20.4); label("Apparent central volume of distribution Vc/F (L)") # Table 2: Vc/F 20.4 L (RSE 3.10%)
    lq <- log(0.398); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2: Clp/F 0.398 L/h (RSE 9.40%)
    lvp <- log(1.75); label("Apparent peripheral volume of distribution Vp/F (L)") # Table 2: Vp/F 1.75 L (RSE 5.50%)

    # Allometric exponents -- estimated; one shared by CL/F and Q/F, one by Vc/F and Vp/F
    e_wt_cl_q <- 0.516; label("Allometric exponent on CL/F and Q/F (unitless)") # Table 2: '(weight/70)^0.516' on Cl/F and on Clp/F (RSE 12.7%)
    e_wt_vc_vp <- 0.671; label("Allometric exponent on Vc/F and Vp/F (unitless)") # Table 2: '(weight/70)^0.671' on Vc/F and on Vp/F (RSE 9.7%)

    # Covariate effects. Categorical effects are modelled as exp(effect)
    # (Table S5 footnote a); Table 2 prints the multiplicative factor exp(effect).
    e_age_cl <- -0.0107; label("Linear slope of age on CL/F, per year from 25 years (1/year)") # Table 2: '1 + (-0.0107 x [age - 25])' (RSE 10.3%)
    e_sexf_cl <- log(0.882); label("Log multiplicative effect of female sex on CL/F (unitless)") # Table 2: 'x 0.882 if female' (RSE 29.2%)
    e_dis_hae_acute_cl <- log(0.911); label("Log multiplicative effect of an HAE attack at dosing on CL/F (unitless)") # Table 2: 'x 0.911 if HAE attack' (RSE 50.2%)
    e_sexf_vc <- log(0.855); label("Log multiplicative effect of female sex on Vc/F (unitless)") # Table 2: 'x 0.855 if female' (RSE 26.8%)
    e_nonwhite_vc <- log(1.11); label("Log multiplicative effect of non-White race on Vc/F (unitless)") # Table 2: 'x 1.11 if nonwhite' (RSE 26.2%)

    # Between-subject variability -- exponential (log-normal) random effects,
    # no covariances reported. Table 2 reports BSV as CV%; converted with
    # omega^2 = log(CV^2 + 1).
    etalka ~ 0.1174 # Table 2: ka BSV 35.3%; log(0.353^2 + 1) = 0.1174
    etaltlag ~ 0.2694 # Table 2: tlag BSV 55.6%; log(0.556^2 + 1) = 0.2694
    etalcl ~ 0.05025 # Table 2: Cl/F BSV 22.7%; log(0.227^2 + 1) = 0.05025
    etalvc ~ 0.06986 # Table 2: Vc/F BSV 26.9%; log(0.269^2 + 1) = 0.06986
    etalq ~ 0.7711 # Table 2: Clp/F BSV 107.8%; log(1.078^2 + 1) = 0.7711
    etalvp ~ 0.2550 # Table 2: Vp/F BSV 53.9%; log(0.539^2 + 1) = 0.2550

    # Residual error -- proportional only in the final model
    propSd <- 0.130; label("Proportional residual error (fraction)") # Table 2: proportional error 13.0%
  })

  model({
    # 1. Covariate terms. Weight is normalised to 70 kg and age centred at
    #    25 years; White is the race reference, so the race effect enters on
    #    (1 - RACE_WHITE).
    age_cl <- 1 + e_age_cl * (AGE - 25)
    cov_cl <- (WT / 70)^e_wt_cl_q * age_cl *
      exp(e_sexf_cl * SEXF + e_dis_hae_acute_cl * DIS_HAE_ACUTE)
    cov_vc <- (WT / 70)^e_wt_vc_vp *
      exp(e_sexf_vc * SEXF + e_nonwhite_vc * (1 - RACE_WHITE))

    # 2. Individual parameters
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + etalvc) * cov_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Absorption lag
    alag(depot) <- tlag

    # 6. Observation and error. Dose in mg and volumes in L give mg/L; x 1000
    #    reports ng/mL, the unit of the source tables.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
