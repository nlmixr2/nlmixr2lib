Yu_2026_tacrolimus <- function() {
  description <- "Two-compartment population PK model with first-order absorption and an absorption lag time for twice-daily oral immediate-release tacrolimus in adult kidney transplant recipients during the first 7 days post-transplant (Yu 2026 Bayesian PopPK comparator model, University of Michigan cohort). The structure and every disposition and absorption parameter (lag time, ka, V1/F, Q/F, V2/F) are fixed to the Francke 2025 model; only the apparent oral clearance CL/F, its inter-individual variability and a proportional residual error were re-estimated in Monolix on 4311 trough concentrations from 1624 recipients. The re-estimated CL/F carries four covariates selected by stepwise forward inclusion and backward elimination: African American race (exp(0.401) = 1.49-fold higher CL/F than non-African American recipients), body weight (power exponent 0.692), age (power exponent -0.307) and haematocrit (power exponent -0.394). The paper does not print the centring values of the three continuous covariates; this implementation centres them on the study-cohort means of Table 1 (87.4 kg, 51.9 years, 27.9 percent). The paper's mechanism-informed GRU and XGBoost machine-learning models are not population PK models and are not included."
  reference <- "Yu H, Qin Z, Smith LS, Park JM, Zhu H-J. Mechanism-informed machine learning for individualized tacrolimus dose adjustment in the early post-kidney transplant period. Br J Clin Pharmacol. 2026;92(6):1730-1743. doi:10.1002/bcp.70448"
  vignette <- "Yu_2026_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as the power term (WT/87.4)^0.692 (Supplementary Table S3 'Weight (kg)' = 0.692, the coefficient of a log-transformed covariate on the log-normal CL/F in Monolix). The centring value is not printed in the paper; 87.4 kg is the study-cohort mean in Table 1 (87.4 +/- 21.1 kg). Collected per patient from the electronic medical record (Methods, Data collection).",
      source_name = "Weight"
    ),
    AGE = list(
      description = "Recipient age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as the power term (AGE/51.9)^-0.307 (Supplementary Table S3 'Age (years)' = -0.307). The centring value is not printed in the paper; 51.9 years is the study-cohort mean in Table 1 (51.9 +/- 13.8 years). All recipients were adults (>= 18 years at transplantation, Methods, Patient population).",
      source_name = "Age"
    ),
    HCT = list(
      description = "Haematocrit, packed red blood cell volume fraction",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters CL/F as the power term (HCT/27.9)^-0.394 (Supplementary Table S3 'Hematocrit (%)' = -0.394). Supplied on the percent scale, as in the paper (Table 1 reports 27.9 +/- 4.3 percent). The centring value is not printed in the paper; 27.9 percent is the study-cohort mean in Table 1. Methods (Data collection) collected haematocrit as a time-varying laboratory value, taking the value closest in time to each dose and trough and interpolating linearly between draws; the column may therefore be supplied time-varying.",
      source_name = "Hematocrit"
    ),
    RACE_BLACK = list(
      description = "African American race indicator: 1 = African American, 0 = Caucasian or Other.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Caucasian or Other race pooled)",
      notes = "Race was recorded as Caucasian, African American or Other (Methods, Data collection; Table 1: 66.1 / 23.7 / 10.2 percent). Supplementary Table S3 prints a single race coefficient, 'Race (AA)' = 0.401, a log-additive effect on CL/F, so African American recipients have exp(0.401) = 1.49-fold higher CL/F. No coefficient is printed for the Other group, so Other is pooled with Caucasian as the reference. The paper's Results attribute the effect to the higher prevalence of CYP3A5 expressors among African American recipients, CYP3A5 genotype itself being unavailable for most of the cohort.",
      source_name = "Race"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1624L,
    n_studies = 1L,
    n_observations = 4311L,
    age_mean_sd = "51.9 +/- 13.8 years (Table 1)",
    weight_mean_sd = "87.4 +/- 21.1 kg (Table 1)",
    height_mean_sd = "171.9 +/- 10.2 cm (Table 1)",
    bmi_mean_sd = "29.5 +/- 6.2 kg/m^2 (Table 1)",
    bsa_mean_sd = "2.03 +/- 0.28 m^2 (Table 1)",
    sex_female_pct = 37.3,
    race_ethnicity = c(Caucasian = 66.1, `African American` = 23.7, Other = 10.2),
    hematocrit_mean_sd = "27.9 +/- 4.3 percent (Table 1)",
    albumin_mean_sd = "3.49 +/- 0.36 g/dL (Table 1)",
    serum_creatinine_mean_sd = "3.54 +/- 2.61 mg/dL (Table 1)",
    co_medication = "Standard regimen of tacrolimus, mycophenolate mofetil and a steroid taper. CYP3A4/5 inhibitor use within 3 days before a trough in 20.4 percent of patients and inducer use in 5.2 percent (Table 1); neither was retained as a covariate.",
    disease_state = "Adult (>= 18 years) kidney transplant recipients in the first 7 days after transplantation, excluding multi-organ transplantation and re-transplantation (Methods, Patient population).",
    dose_range = "Oral immediate-release tacrolimus (Prograf or generic) every 12 h, started within 24 h of transplantation at 0.075 mg/kg per dose and then adjusted by trough-guided therapeutic drug monitoring; administered daily dose 11.8 +/- 3.6 mg/day (Table 1).",
    sampling_window = "Trough concentrations only, within the first 7 days post-transplant; 2.7 +/- 1.0 troughs per patient, most in the 24-72 h interval (Supplementary Table S1). Observed trough 9.7 +/- 4.0 ng/mL (Table 1).",
    assay = "Whole-blood tacrolimus by LC-MS/MS (LLOQ 1.0 ng/mL), Michigan Medicine Laboratories.",
    regions = "United States (University of Michigan Hospital, Ann Arbor, MI), transplants January 2001 to January 2025.",
    external_validation = "MIMIC-IV (Beth Israel Deaconess Medical Center, 2008-2019): 164 kidney transplant recipients, 611 troughs (Table 1).",
    notes = "Retrospective single-centre cohort. CYP3A5 genotype was available for only 23.2 percent of patients and was excluded, so the genotype covariates of the Francke 2025 parent model were removed (Methods, Population pharmacokinetic modelling)."
  )

  ini({
    # Supplementary Table S3. The structural parameters lag time, ka, V1/F,
    # Q/F and V2/F are fixed to the published Francke 2025 estimates (Methods,
    # Population pharmacokinetic modelling); Table S3 prints them without a
    # standard error. Only CL/F, its covariate effects, its IIV and the
    # residual error were estimated (Monolix 2024R1, SAEM). Bioavailability is
    # not identifiable from oral data, so all clearances and volumes are
    # apparent values (CL/F, V1/F, Q/F, V2/F).
    ltlag <- fixed(log(0.375)); label("Absorption lag time tlag (h)") # Table S3 Tlag_pop = 0.375 (no S.E.: fixed to Francke 2025)
    lka <- fixed(log(6.59)); label("First-order absorption rate constant ka (1/h)") # Table S3 ka_pop = 6.59 1/h (no S.E.: fixed to Francke 2025)
    lcl <- log(24.3); label("Apparent oral clearance CL/F at the reference covariate set (L/h)") # Table S3 CL_pop = 24.3 L/h (S.E. 0.434, RSE 1.78%)
    lvc <- fixed(log(705)); label("Apparent central volume V1/F (L)") # Table S3 V1_pop = 705 L (no S.E.: fixed to Francke 2025)
    lq <- fixed(log(8.54)); label("Apparent inter-compartmental clearance Q/F (L/h)") # Table S3 Q_pop = 8.54 L/h (no S.E.: fixed to Francke 2025)
    lvp <- fixed(log(7670)); label("Apparent peripheral volume V2/F (L)") # Table S3 V2_pop = 7670 L (no S.E.: fixed to Francke 2025)

    # Covariate effects on log(CL/F). Monolix coefficients of log-transformed
    # continuous covariates are power exponents; the categorical race
    # coefficient is additive on log(CL/F).
    e_age_cl <- -0.307; label("Power exponent of age on CL/F (unitless)") # Table S3 Age (years) = -0.307 (RSE 15.4%)
    e_wt_cl <- 0.692; label("Power exponent of body weight on CL/F (unitless)") # Table S3 Weight (kg) = 0.692 (RSE 9.37%)
    e_race_black_cl <- 0.401; label("Log-additive African American effect on CL/F (unitless)") # Table S3 Race (AA) = 0.401 (RSE 8.30%)
    e_hct_cl <- -0.394; label("Power exponent of haematocrit on CL/F (unitless)") # Table S3 Hematocrit (%) = -0.394 (RSE 24.3%)

    # Table S3 omega_CL = 0.514 is the SD of the log-normal random effect
    # (Table S3 C.V. 54.9% = sqrt(exp(0.514^2) - 1)); variance = 0.514^2.
    etalcl ~ 0.264196 # Table S3 omega_CL = 0.514 (SD), CV 54.9%

    propSd <- 0.281; label("Proportional residual error (fraction)") # Table S3 b = 0.281 (RSE 1.39%)
  })

  model({
    # Continuous covariate power terms on CL/F, centred on the study-cohort
    # means of Table 1 (the paper does not print the Monolix centring values).
    # HCT is on the percent scale used by the paper.
    f_wt <- (WT / 87.4)^e_wt_cl
    f_age <- (AGE / 51.9)^e_age_cl
    f_hct <- (HCT / 27.9)^e_hct_cl

    tlag <- exp(ltlag)
    ka <- exp(lka)
    cl <- exp(lcl + e_race_black_cl * RACE_BLACK + etalcl) * f_wt * f_age * f_hct
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg and vc in L give mg/L; multiply by 1000 for ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
