Damle_2021_dalteparin <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for subcutaneous dalteparin in pediatric patients (1 month to 19 years) with venous thromboembolism, fitted to plasma anti-factor Xa activity as the surrogate concentration (Damle 2021 full covariate model). CL/F and V/F scale allometrically with body weight (exponents fixed at 0.75 and 1, reference 43 kg); CL/F additionally has a power effect of age (reference 12 years) and multiplicative male-sex and no-cancer factors (reference: a female patient with cancer). The IIV of V/F is the IIV of CL/F multiplied by an estimated scaling factor (perfect correlation); combined proportional plus additive residual error."
  reference <- "Damle B, Jen F, Sherman N, Jani D, Sweeney K. Population Pharmacokinetic Analysis of Dalteparin in Pediatric Patients With Venous Thromboembolism. J Clin Pharmacol. 2021;61(2):172-180. doi:10.1002/jcph.1716"
  vignette <- "Damle_2021_dalteparin"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  compartmentData <- list(
    depot = list(
      analyte = "dalteparin (anti-Xa activity)",
      units = "IU",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "dalteparin (anti-Xa activity)", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power on CL/F (exponent fixed at 0.75) and V/F (exponent fixed at 1), reference weight 43 kg (Damle 2021 model equations, Results 'Population Pharmacokinetic Analysis'; cohort median 43.4 kg, Table 1).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F, (AGE / 12)^e_age_cl with e_age_cl = -0.0687, reference 12 years (cohort median, Table 1). The youngest subject in the analysis data set was about 0.04 years (15 days); age must be strictly positive.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female)",
      notes = "Damle 2021 codes SEX = 1 for male and applies theta14 = 1.03 to CL/F as theta14^I(SEX = male); encoded here as e_sexm_cl^(1 - SEXF). The bootstrap 95% CI of theta14 (0.908-1.2) includes 1 and the sex effect was dropped from the reduced model used for the paper's dose simulations.",
      source_name = "SEX"
    ),
    DIS_CANCER_PED = list(
      description = "Pediatric cancer status, 1 = patient with cancer, 0 = patient without cancer",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (with cancer)",
      notes = "Damle 2021 codes CANCERST = 1 for patients WITHOUT cancer and applies theta15 = 0.885 to CL/F as theta15^I(CANCERST = without cancer); encoded here as e_nocancer_cl^(1 - DIS_CANCER_PED), so the reference patient has cancer. 48.3% of the cohort had cancer (Table 1). The bootstrap 95% CI of theta15 (0.744-1.02) includes 1 and the cancer effect was dropped from the reduced model used for the paper's dose simulations.",
      source_name = "CANCERST"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 89L,
    n_studies = 3L,
    n_observations = 266L,
    age_range = "0.04-19.5 years (15 days to 19.5 years)",
    age_median = "12.0 years",
    weight_range = "2.3-161 kg",
    weight_median = "43.4 kg",
    sex_female_pct = 33.7,
    race_ethnicity = c(White = 45),
    disease_state = "Pediatric patients with acute venous thromboembolism requiring therapeutic anticoagulation, with (48.3%) or without cancer.",
    dose_range = "Subcutaneous dalteparin twice daily, starting doses 100-150 IU/kg by age group, titrated in 25 IU/kg steps (or by 10-20%) to a 4-6 h post-dose anti-Xa target of 0.5-1.0 IU/mL; the Mayo Clinic cohort received 100 IU/kg twice daily or 200 IU/kg once daily.",
    regions = "North America and Europe (15 sites in the Pfizer study NCT00952380), plus the multicenter Kids-DOTT pilot and a Mayo Clinic retrospective chart review.",
    age_groups = "0 to <8 weeks n = 6; 8 weeks to <2 years n = 13; 2 to <8 years n = 14; 8 to <12 years n = 11; 12 to <19 years n = 45 (Table 1).",
    notes = "Pooled from the Pfizer-sponsored open-label study (n = 37), the Kids-DOTT dose-finding pilot (n = 18) and a Mayo Clinic retrospective chart analysis (n = 34). The dependent variable is chromogenic plasma anti-Xa activity (IU/mL), measured centrally for the Pfizer study and locally for the external studies; below-quantitation samples were excluded. Most samples were collected 3-6 h post-dose."
  )

  ini({
    # Structural parameters: Damle 2021 Table 2 (full model), converted from
    # mL to L so that Cc = central / (vc * 1000) is in IU/mL.
    lka <- log(1.04); label("Absorption rate constant ka (1/h)") # Table 2 Ka theta3 = 1.04 1/h (%RSE 71.8)
    lcl <- log(0.929); label("Apparent clearance CL/F for a 43 kg, 12-year-old female with cancer (L/h)") # Table 2 CL/F theta1 = 929 mL/h
    lvc <- log(7.18); label("Apparent volume V/F for a 43 kg patient (L)") # Table 2 V/F theta2 = 7180 mL

    # Covariate effects
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)") # Table 2 WT in CL/F theta6 = 0.75 (fixed)
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on V/F (unitless)") # Table 2 WT in V/F theta7 = 1 (fixed)
    e_age_cl <- -0.0687; label("Power exponent of age (years / 12) on CL/F (unitless)") # Table 2 AGE in CL/F theta8 = -0.0687
    e_sexm_cl <- 1.03; label("Multiplicative factor on CL/F for male sex (unitless)") # Table 2 SEX = 1 in CL/F theta14 = 1.03
    e_nocancer_cl <- 0.885; label("Multiplicative factor on CL/F for patients without cancer (unitless)") # Table 2 CANCERST = 1 in CL/F theta15 = 0.885

    # eta(V/F) = theta5 * eta(CL/F): perfect correlation, SD ratio theta5.
    # Encoded as a structural scaling of etalcl (same form as Hirt 2009
    # efavirenz) because a correlation-1 OMEGA block is singular.
    vc_eta_scale <- 1.73; label("Scaling factor relating eta V/F to eta CL/F (eta V/F = vc_eta_scale * etalcl)") # Table 2 Scaling of IIV for V/F theta5 = 1.73

    etalcl ~ 0.0369 # Table 2 omega2 CL/F = 0.0369 (19% CV)

    # Residual error: Table 2 reports both terms as variances (sigma2).
    propSd <- 0.2278; label("Proportional residual error (fraction)") # Table 2 sigma2 P = 0.0519; sqrt = 0.228 (text: 23% CV)
    addSd <- 0.1265; label("Additive residual error (IU/mL)") # Table 2 sigma2 A = 0.016; sqrt = 0.1265 IU/mL
  })

  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 43)^e_wt_cl * (AGE / 12)^e_age_cl *
      e_sexm_cl^(1 - SEXF) * e_nocancer_cl^(1 - DIS_CANCER_PED)
    vc <- exp(lvc + vc_eta_scale * etalcl) * (WT / 43)^e_wt_vc

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in IU, vc in L -> IU/L; divide by 1000 for anti-Xa in IU/mL.
    Cc <- central / (vc * 1000)
    Cc ~ add(addSd) + prop(propSd)
  })
}
