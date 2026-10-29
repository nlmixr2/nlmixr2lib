Liu_2024_cetirizine <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for oral cetirizine in Chinese children (about 0.5-12 years) with allergic diseases (Liu 2024). The structure is the Pitsiu 2004 paediatric cetirizine model re-fitted to sparse real-world Chinese data: CL/F is linear in age with sex-specific intercepts (CL/F = 0.28 * AGE + 0.771 [male] or 0.592 [female] L/h), V/F is linear in age (V/F = 23.1 + 1.42 * AGE L), and ka = theta_ka + CL/V guards against flip-flop. The age slope on CL/F, the V/F intercept, the IIV on CL/F and V/F and the proportional residual error were re-estimated; the sex intercepts, the V/F age slope, theta_ka and its IIV were held fixed at the Pitsiu 2004 values."
  reference <- paste(
    "Liu W, Tan Z, Yang P, Yu Z, Yao X, Zhou P, Liu L, Zhou W.",
    "Model-based exploration of the rationality of off-label use of",
    "cetirizine in Chinese pediatric patients: a prospective cohort study.",
    "Front Pharmacol. 2024;15:1322788. doi:10.3389/fphar.2024.1322788.",
    "Final parameter estimates are Supplementary Table S1. The covariate",
    "equations are not printed in Liu 2024 and are taken from the base model",
    "it re-fitted: Pitsiu M, Hussein Z, Majid O, Aarons L, de Longueville M,",
    "Stockis A. Retrospective population pharmacokinetic analysis of",
    "cetirizine in children aged 6 months to 12 years. Br J Clin Pharmacol.",
    "2004;57(4):402-411. doi:10.1046/j.1365-2125.2003.02017.x (Table 4).",
    sep = " "
  )
  vignette <- "Liu_2024_cetirizine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F and V/F as an uncentred linear term (Pitsiu 2004 Table 4 equations, re-fitted by Liu 2024): CL/F = 0.28 * AGE + sex intercept, V/F = 23.1 + 1.42 * AGE. Liu 2024 enrolled children under 14 years; the age groups of Table 1 span 0.5-12 years and the Table 2 simulations span 1-12 years.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selects the CL/F intercept: 0.771 L/h for males and 0.592 L/h for females (Liu 2024 Table S1 'CL_SEX on Male' / 'CL_SEX on Female', both fixed). Pitsiu 2004 coded SEX = 1 for female, which is the canonical SEXF directly.",
      source_name = "SEX"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cetirizine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cetirizine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 63L,
    n_studies = 1L,
    age_range = "Under 14 years at enrolment; Table 1 age groups 0.5-2, 2-6 and 6-12 years",
    age_mean = "6.34 +/- 3.7 years",
    weight_mean = "26.4 +/- 15.7 kg",
    height_mean = "120.1 +/- 27.7 cm",
    sex_female_pct = 39.7,
    race_ethnicity = "Chinese (100%)",
    disease_state = "Children with allergic diseases (allergic rhinitis, urticaria, asthma and others) prescribed oral cetirizine; normal hepatic and renal function by inclusion criteria.",
    dose_range = "Label doses 2.5 mg bid (0.5-6 years) or 5 mg bid (6-12 years); 24 of 63 children received off-label doses up to twice the label dose (3-7.5 mg bid; Table 1), five of them above 0.25 mg/kg per administration.",
    regions = "China (Peking University Third Hospital, Beijing)",
    samples = "Up to three plasma samples per child; 85 concentrations in the NPDE evaluation (Supplementary Figure S4). Plasma cetirizine by HPLC-MS/MS.",
    estimation = "NONMEM 7.4.4, FOCE-I",
    notes = "Prospective cohort study, April 2020 - April 2021: 38 male and 25 female children, 39 at label doses and 24 at off-label doses (Liu 2024 Results 3.1). Table 1 lists 10 + 22 + 33 = 65 children across age groups and the safety analysis refers to 66 patients; the modelling text states 63."
  )

  ini({
    # Re-estimated by Liu 2024 (Supplementary Table S1). The structure and
    # every value marked FIXED are the Pitsiu 2004 final model (Table 4):
    #   CL/F (L/h) = e_age_cl * AGE + CL_M * (1 - SEXF) + CL_F * SEXF
    #   V/F  (L)   = V + e_age_vc * AGE
    #   ka   (1/h) = theta_ka + CL/V
    lcl <- fixed(log(0.771)); label("CL/F intercept at AGE = 0, males (L/h)") # Table S1 'CL_SEX on Male' 0.771 FIXED (Pitsiu 2004 Table 4 theta_CL,M)
    lcl_female <- fixed(log(0.592)); label("CL/F intercept at AGE = 0, females (L/h)") # Table S1 'CL_SEX on Female' 0.592 FIXED (Pitsiu 2004 Table 4 theta_CL,F)
    e_age_cl <- 0.28; label("Linear increase of CL/F per year of age (L/h/year)") # Table S1 'CL_AGE' 0.28 (RSE 17%)
    lvc <- log(23.1); label("V/F intercept at AGE = 0 (L)") # Table S1 'V' 23.1 (RSE 40%)
    e_age_vc <- fixed(1.42); label("Linear increase of V/F per year of age (L/year)") # Table S1 'V_AGE' 1.42 FIXED (Pitsiu 2004 Table 4 theta_V,AGE)
    lka <- fixed(log(2.14)); label("Additive term theta_ka in ka = theta_ka + CL/V (1/h)") # Table S1 'Ka' 2.14 FIXED (Pitsiu 2004 Table 4 theta_ka)

    # IIV variances (exponential model). Table S1 reports OMEGA variances:
    # the fixed ka entry 2.1 equals the Pitsiu 2004 145% CV squared.
    etalcl ~ 0.378 # Table S1 'IIV_CL/F' 0.378 (RSE 34%, shrinkage 21.2%)
    etalvc ~ 0.445 # Table S1 'IIV_V/F' 0.445 (RSE 41%, shrinkage 47.1%)
    etalka ~ fixed(2.1) # Table S1 'IIV_Ka' 2.1, held at the Pitsiu 2004 value

    # Proportional error coded as a THETA with SIGMA fixed to 1, so the
    # THETA is the proportional SD.
    propSd <- 0.36; label("Proportional residual error (fraction)") # Table S1 'Proportional error' 0.36 (RSE 16%); 'Sigma' 1 FIXED
  })

  model({
    cl <- (e_age_cl * AGE + exp(lcl) * (1 - SEXF) + exp(lcl_female) * SEXF) * exp(etalcl)
    vc <- (exp(lvc) + e_age_vc * AGE) * exp(etalvc)
    ka <- exp(lka + etalka) + cl / vc
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and V/F in L give mg/L; x 1000 gives ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
