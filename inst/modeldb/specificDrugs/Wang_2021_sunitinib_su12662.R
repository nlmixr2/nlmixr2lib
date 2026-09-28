Wang_2021_sunitinib_su12662 <- function() {
  description <- paste(
    "Two-compartment population PK model for SU012662 (N-desethyl",
    "sunitinib, the primary active metabolite) after oral sunitinib in 65",
    "children, adolescents and young adults (3-21 years) with",
    "gastrointestinal stromal tumors or other solid / CNS tumors (Wang 2021).",
    "Fitted separately from the parent: the sunitinib dose enters a",
    "first-order absorption depot with a lag time, with a fixed 21% of the",
    "sunitinib dose assumed converted to SU012662, so the metabolite",
    "clearances and volumes are apparent values relative to that fraction",
    "and the unknown oral bioavailability. CL/F and Vc/F scale with body",
    "surface area as power functions normalised to 1.44 m^2. Exponential",
    "IIV on CL/F, Vc/F and ka; proportional residual error. The parent",
    "model ships as Wang_2021_sunitinib."
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
        "Baseline BSA. Enters SU012662 CL/F and Vc/F as power functions",
        "normalised to 1.44 m^2 (Results 3.3: CL/F = 11.1 l/h *",
        "(BSA/1.44)^0.87 and Vc/F = 1060 l * (BSA/1.44)^1.61). 1.44 m^2 is",
        "close to the cohort median of 1.4 m^2 (Table 2). The BSA formula is",
        "not stated in the source. Selected over baseline body weight on",
        "objective function."
      ),
      source_name = "BSA"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in place of BSA on SU012662 CL/F and Vc/F in a separate stepwise run; the BSA model had the lower objective function and was retained (Results 3.3)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F and Vc/F in the stepwise covariate model; not significant (P > 0.001; Results 3.3)."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and Vc/F; not significant (Results 3.3)."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F (Asian vs non-Asian); not significant (Results 3.3)."
    ),
    TUMTP_GIST = list(
      description = "Tumor type gastrointestinal stromal tumor (vs other solid tumor)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F and Vc/F; not significant (Results 3.3)."
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status > 0 (vs 0)",
      units = "(binary)",
      type = "binary",
      notes = "Screened on CL/F as 0 vs > 0 (Karnofsky-extrapolated where needed); not significant (Results 3.3)."
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
    central_su12662 = list(analyte = "SU12662", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_su12662 = list(analyte = "SU12662", units = "mg", specimen = "plasma", verified = TRUE)
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
      "GIST, n = 6) -- Table 1. Demographics by age group in Table 2. 417",
      "post-baseline SU012662 observations; two |CWRES| > 6 outliers",
      "excluded. NONMEM 7.1.2, FOCE-I."
    )
  )

  ini({
    # Fraction of the sunitinib dose converted to SU012662 -- assumed, not
    # estimated (Results 3.3: 'a conversion of 21% sunitinib to SU012662
    # was assumed', citing the adult analysis of Houk 2009).
    fm <- fixed(0.21)
    label("Fraction of sunitinib dose converted to SU012662 (fraction)") # Results 3.3 '21%' (assumed)

    # Structural parameters -- Table 3 SU012662 'Results, mean (RSE %)'
    # column; typical values at BSA = 1.44 m^2.
    lcl_su12662 <- log(11.1)
    label("SU012662 apparent clearance CL/F at BSA 1.44 m^2 (L/h)") # Table 3 CL/F (theta1) 11.1 (RSE 6.9%)
    lvc_su12662 <- log(1060)
    label("SU012662 apparent central volume Vc/F at BSA 1.44 m^2 (L)") # Table 3 Vc/F (theta2) 1060 (RSE 14%)
    lka_su12662 <- log(0.28)
    label("SU012662 model absorption rate constant ka (1/h)") # Table 3 ka (theta3) 0.28 (RSE 36.7%)
    ltlag_su12662 <- log(0.64)
    label("SU012662 model absorption lag time (h)") # Table 3 tlag (theta4) 0.64 (RSE 26.3%)
    lvp_su12662 <- log(63.1)
    label("SU012662 apparent peripheral volume Vp/F (L)") # Table 3 Vp/F (theta5) 63.1 (RSE 141%)
    lq_su12662 <- log(6.7)
    label("SU012662 apparent intercompartmental clearance Q/F (L/h)") # Table 3 Q/F (theta6) 6.7 (RSE 319%)

    # BSA power exponents -- Results 3.3 printed final equations.
    e_bsa_cl_su12662 <- 0.87
    label("Power exponent of BSA/1.44 on SU012662 CL/F (unitless)") # Results 3.3 CL/F = 11.1 l/h * (BSA/1.44)^0.87; Table 3 theta9 0.87 (RSE 26%)
    e_bsa_vc_su12662 <- 1.61
    label("Power exponent of BSA/1.44 on SU012662 Vc/F (unitless)") # Results 3.3 Vc/F = 1060 l * (BSA/1.44)^1.61; Table 3 theta8 1.61 (RSE 20%)

    # IIV -- Table 3 omega rows reported as %; variances via
    # omega^2 = log(1 + CV^2). No omega block (Results 3.3: eta correlations
    # judged weak).
    etalcl_su12662 ~ 0.1770 # Table 3 omega CL/F 44% -> log(1 + 0.44^2)
    etalvc_su12662 ~ 0.1625 # Table 3 omega Vc/F 42% -> log(1 + 0.42^2)
    etalka_su12662 ~ 0.6502 # Table 3 omega ka 95.7% -> log(1 + 0.957^2)

    # Residual error -- Table 3 sigma (theta7) 26%, a THETA-scaled
    # proportional error.
    propSd_su12662 <- 0.26
    label("SU012662 proportional residual error (fraction)") # Table 3 sigma (theta7) 26% (RSE 3.24%)
  })

  model({
    # Individual parameters (theta_i = theta * exp(eta_i), Methods 2.3)
    cl_su12662 <- exp(lcl_su12662 + etalcl_su12662) * (BSA / 1.44)^e_bsa_cl_su12662
    vc_su12662 <- exp(lvc_su12662 + etalvc_su12662) * (BSA / 1.44)^e_bsa_vc_su12662
    ka_su12662 <- exp(lka_su12662 + etalka_su12662)
    tlag_su12662 <- exp(ltlag_su12662)
    vp_su12662 <- exp(lvp_su12662)
    q_su12662 <- exp(lq_su12662)

    kel_su12662 <- cl_su12662 / vc_su12662
    k12_su12662 <- q_su12662 / vc_su12662
    k21_su12662 <- q_su12662 / vp_su12662

    # The depot receives the sunitinib dose (mg); fm of it enters the
    # SU012662 system (mass basis, no molecular-weight correction).
    d/dt(depot) <- -ka_su12662 * depot
    d/dt(central_su12662) <- ka_su12662 * depot -
      kel_su12662 * central_su12662 -
      k12_su12662 * central_su12662 +
      k21_su12662 * peripheral1_su12662
    d/dt(peripheral1_su12662) <- k12_su12662 * central_su12662 -
      k21_su12662 * peripheral1_su12662

    f(depot) <- fm
    alag(depot) <- tlag_su12662

    # mg / L -> ng/mL
    Cc_su12662 <- 1000 * central_su12662 / vc_su12662
    Cc_su12662 ~ prop(propSd_su12662)
  })
}
