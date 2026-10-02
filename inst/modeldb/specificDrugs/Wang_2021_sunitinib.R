Wang_2021_sunitinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral sunitinib (parent drug) in",
    "65 children, adolescents and young adults (3-21 years) with",
    "gastrointestinal stromal tumors or other solid / CNS tumors (Wang 2021).",
    "First-order absorption with a lag time; apparent clearance CL/F and",
    "central volume Vc/F scale with body surface area as power functions",
    "normalised to 1.44 m^2. Exponential IIV on CL/F, Vc/F and ka;",
    "proportional residual error. The active metabolite SU012662 was fitted",
    "as a separate model and ships as Wang_2021_sunitinib_su12662."
  )
  reference <- paste(
    "Wang E, DuBois SG, Wetmore C, Verschuur AC, Khosravan R (2021).",
    "Population Pharmacokinetics of Sunitinib and its Active Metabolite",
    "SU012662 in Pediatric Patients with Gastrointestinal Stromal Tumors or",
    "Other Solid Tumors. Eur J Drug Metab Pharmacokinet 46(3):343-352.",
    "doi:10.1007/s13318-021-00671-7."
  )
  vignette <- "Wang_2021_sunitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline BSA. Enters CL/F and Vc/F as power functions normalised to",
        "1.44 m^2 (Results 3.2: CL/F = 24 l/h * (BSA/1.44)^0.733 and",
        "Vc/F = 1030 l * (BSA/1.44)^1.46). 1.44 m^2 is close to the cohort",
        "median of 1.4 m^2 (Table 2). The BSA formula is not stated in the",
        "source. Selected over baseline body weight on objective function."
      ),
      source_name = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in place of BSA on CL/F and Vc/F in a separate stepwise run; the BSA model had the lower objective function and was retained (Results 3.2)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F and Vc/F in the stepwise covariate model; not significant (P > 0.001; Results 3.2)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and Vc/F; not significant (Results 3.2)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F (Asian vs non-Asian); not significant (Results 3.2)."
    ),
    TUMTP_GIST = list(
      description = "Tumor type gastrointestinal stromal tumor (vs other solid tumor)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and Vc/F; not significant (Results 3.2). Post hoc CL/F and Vc/F were similar in GIST and other tumors (Discussion)."
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status > 0 (vs 0)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F as 0 vs > 0 (Karnofsky-extrapolated where needed); not significant (Results 3.2)."
    ),
    FORM_SUNITINIB_SPRINKLE = list(
      description = "Formulation: capsule contents sprinkled on yogurt or applesauce (vs intact capsule)",
      units = "(binary)",
      type = "binary",
      notes = "The only covariate predefined for ka (Methods 2.3); not retained in the final model."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sunitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sunitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 65L,
    n_studies = 3L,
    age_range = "3-21 years (studies enrolled 18 months to 22 years)",
    age_median = "13 years",
    weight_range = "16.2-100 kg",
    weight_median = "49.1 kg",
    bsa_range = "0.7-2.1 m^2",
    bsa_median = "1.4 m^2",
    sex_female_pct = 50.8,
    race_ethnicity = c(Asian = 6.2, NonAsian = 89.2, Unknown = 4.6),
    disease_state = "Pediatric gastrointestinal stromal tumor (n = 6) or other refractory solid / CNS tumors (n = 59; primarily high-grade glioma, ependymoma, sarcoma)",
    dose_range = "Oral sunitinib 15 or 20 mg/m^2 once daily on schedule 4/2 (4 weeks on, 2 weeks off), as intact capsule or capsule contents sprinkled on yogurt / applesauce",
    regions = "North America and Europe",
    notes = paste(
      "Pooled studies ADVL0612 (NCT00387920, phase 1, n = 35), ACNS1021",
      "(NCT01462695, phase 2, n = 24) and A6181196 (NCT01396148, phase 1/2",
      "GIST, n = 6) -- Table 1. Demographics by age group in Table 2. 439",
      "post-baseline sunitinib observations; four |CWRES| > 6 outliers",
      "excluded. NONMEM 7.1.2, FOCE-I."
    )
  )

  ini({
    # Structural parameters -- Table 3 sunitinib 'Results, mean (RSE %)'
    # column; typical values at BSA = 1.44 m^2.
    lcl <- log(24.0)
    label("Apparent clearance CL/F at BSA 1.44 m^2 (L/h)") # Table 3 CL/F (theta1) 24.0 (RSE 5.8%)
    lvc <- log(1030)
    label("Apparent central volume Vc/F at BSA 1.44 m^2 (L)") # Table 3 Vc/F (theta2) 1030 (RSE 9.8%)
    lka <- log(0.37)
    label("Absorption rate constant ka (1/h)") # Table 3 ka (theta3) 0.37 (RSE 28.3%)
    ltlag <- log(0.76)
    label("Absorption lag time (h)") # Table 3 tlag (theta4) 0.76 (RSE 3.6%)
    lvp <- log(81.3)
    label("Apparent peripheral volume Vp/F (L)") # Table 3 Vp/F (theta5) 81.3 (RSE 22.9%)
    lq <- log(0.39)
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 Q/F (theta6) 0.39 (RSE 79.6%)

    # BSA power exponents -- Results 3.2 printed final equations (Table 3
    # rounds the CL/F exponent to 0.73).
    e_bsa_cl <- 0.733
    label("Power exponent of BSA/1.44 on CL/F (unitless)") # Results 3.2 CL/F = 24 l/h * (BSA/1.44)^0.733; Table 3 theta9 0.73 (RSE 25.6%)
    e_bsa_vc <- 1.46
    label("Power exponent of BSA/1.44 on Vc/F (unitless)") # Results 3.2 Vc/F = 1030 l * (BSA/1.44)^1.46; Table 3 theta8 1.46 (RSE 19.9%)

    # IIV -- Table 3 omega rows reported as %; variances via
    # omega^2 = log(1 + CV^2). No omega block (Results 3.2: eta correlations
    # judged weak).
    etalcl ~ 0.1034 # Table 3 omega CL/F 33% -> log(1 + 0.33^2)
    etalvc ~ 0.06204 # Table 3 omega Vc/F 25.3% -> log(1 + 0.253^2)
    etalka ~ 0.7271 # Table 3 omega ka 103.4% -> log(1 + 1.034^2)

    # Residual error -- Table 3 sigma (theta7) 32.2%, a THETA-scaled
    # proportional error.
    propSd <- 0.322
    label("Proportional residual error (fraction)") # Table 3 sigma (theta7) 32.2% (RSE 2.3%)
  })

  model({
    # Individual parameters (theta_i = theta * exp(eta_i), Methods 2.3)
    cl <- exp(lcl + etalcl) * (BSA / 1.44)^e_bsa_cl
    vc <- exp(lvc + etalvc) * (BSA / 1.44)^e_bsa_vc
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag)
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
