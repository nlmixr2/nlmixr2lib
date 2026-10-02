Clements_2020_blinatumomab <- function() {
  description <- "One-compartment population PK model with linear first-order elimination for continuous intravenous blinatumomab (CD19/CD3 bispecific T cell engager) in pediatric and adult patients with hematological malignancies (Clements 2020). Fit by NONMEM 7.2 FOCE to serum concentrations from 674 patients (628 adults, 46 children aged 7 months to 16 years) pooled from eight phase I-III studies in relapsed/refractory B-precursor ALL (Philadelphia-negative and -positive), MRD-positive B-lineage ALL and relapsed NHL. Typical CL is 2.22 L/h and V 5.98 L at the reference body surface area of 1.876 m2; body surface area enters CL as a power function (exponent 0.620). Inter-individual variability is carried on CL only; the residual error is additive on the natural-log scale (transform-both-sides) and itself carries inter-individual variability on its magnitude."

  reference <- paste(
    "Clements JD, Zhu M, Kuchimanchi M, Terminello B, Doshi S. (2020).",
    "Population Pharmacokinetics of Blinatumomab in Pediatric and Adult",
    "Patients with Hematological Malignancies.",
    "Clinical Pharmacokinetics 59(4):463-474.",
    "doi:10.1007/s40262-019-00823-8.",
    sep = " "
  )

  vignette <- "Clements_2020_blinatumomab"

  # CL is reported in L/h and V in L (Table 3). Doses in ug give ug/L in the
  # central compartment; the observation multiplies by 1000 so Cc is in pg/mL,
  # the unit of the bioassay (LLOQ 50-100 pg/mL) and of every concentration
  # the paper reports (Table 2, Figure 2, Figure 3).
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")

  covariateData <- list(
    BSA = list(
      description = "Body surface area. Enters CL as the power function (BSA / 1.876)^0.620.",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 1.876 m2 is printed in Table 3 footnote a ('CL individual = CL * (BSA/1.876)^Effect of BSA on CL'); the Results text rounds it to the 'median BSA (1.88 m2)'. Cohort median 1.8 m2 (range 0.4-2.7; 0.37-2.70 per Discussion), Electronic Supplementary Material (ESM) Table 1. BSA was the only covariate retained; age, creatinine clearance, sex, AST, ALT, total bilirubin, albumin, LDH, hemoglobin, dose level and treatment cycle were screened graphically against the CL empirical Bayes estimates (no r2 above 3 percent) and never tested formally in NONMEM. The BSA formula used by the study sites is not stated.",
      source_name = "BSA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "blinatumomab", units = "ug", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 674L,
    n_studies = 8L,
    n_observations = 3629L,
    age_range = "0.6-80 years",
    age_median = "41.0 years",
    weight_range = "7.5-148.7 kg",
    weight_median = "70.7 kg",
    bsa_range = "0.37-2.70 m2",
    bsa_median = "1.8 m2 (reference 1.876 m2)",
    sex_female_pct = 39.7,
    race_ethnicity = "85.9 percent White, 3.2 percent Asian, 2.0 percent Black or African American, 0.5 percent American Indian or Alaska Native, 7.8 percent other (ESM Table 2)",
    disease_state = "Hematological malignancies: relapsed/refractory Philadelphia-negative B-precursor ALL in adults (n = 472) and children (n = 46), relapsed/refractory Philadelphia-positive ALL (n = 37), MRD-positive B-lineage ALL (n = 52) and relapsed non-Hodgkin lymphoma (n = 67).",
    age_groups = "628 adults (93.1 percent), 3 adolescents, 33 children, 10 infants (ESM Table 2)",
    renal_function = "Creatinine clearance median 121 mL/min (range 36-150; Cockcroft-Gault in adults, Schwartz in children)",
    dose_range = "Continuous IV infusion over 4-8 weeks per cycle, 0.5-90 ug/m2/day (BSA-based) or 9-28 ug/day (fixed), up to ten cycles",
    regions = "Multinational (Amgen studies MT103-104, -202, -203, -205, -206, -211, 20120216 and 00103311)",
    notes = "Demographics from Results section 3.1 and ESM Tables 1-2. Of 4841 serum samples, 548 were excluded (missing time, more than 90 days after infusion start, or more than 1 day after infusion end) and 664 were below the LLOQ, leaving 3629 in the analysis dataset. (The abstract quotes 2417 concentrations; Results 3.1 is used here.) Concentrations were measured with a CD69 T cell activation bioassay (LLOQ 50-100 pg/mL, ULOQ 1000 pg/mL)."
  )

  ini({
    # =========================================================================
    # Structural fixed effects -- Clements 2020 Table 3 (final model), at the
    # reference BSA of 1.876 m2 (Table 3 footnote a).
    # =========================================================================
    lcl <- log(2.22); label("Clearance CL at the reference BSA of 1.876 m2 (L/h)") # Table 3 'CL (L/h)' = 2.22 (RSE 2.95%, 95% CI 2.08-2.35)
    lvc <- log(5.98); label("Central volume of distribution V (L)") # Table 3 'Volume (L)' = 5.98 (RSE 8.86%, 95% CI 5.14-6.98)

    # Table 3 footnote a: CL_individual = CL * (BSA/1.876)^theta_BSA
    e_bsa_cl <- 0.620; label("Power exponent of body surface area on CL, normalised to 1.876 m2 (unitless)") # Table 3 'Effect of BSA on CL (theta)' = 0.620 (RSE 12.7%, 95% CI 0.46-0.76)

    # =========================================================================
    # Inter-individual variability -- Table 3 reports %CV only. ESM
    # 'Pharmacostatistical Modeling': IIV is exponential (log-normal) and 'the
    # magnitude of IIV and [residual variability] was expressed approximately
    # as the percent coefficient of variation', i.e. %CV/100 = sqrt(omega^2),
    # so omega^2 = (CV/100)^2. No IIV was estimated on V (Results 3.3).
    # =========================================================================
    etalcl ~ 0.2266 # Table 3 'omega CL' = 47.6 %CV (RSE 16.1%); variance = 0.476^2

    # =========================================================================
    # Residual error -- ESM: 'additive error model after natural logarithmic
    # transformation of the measured blinatumomab concentrations and model
    # predictions (the transform at both sides approach)', i.e. nlmixr2
    # lnorm() with the canonical expSd. Results 3.3: 'Residual variability was
    # modeled using an additive error model in the log-domain with IIV', i.e.
    # NONMEM W = THETA * EXP(ETA), so the per-subject log-scale residual SD is
    # expSd * exp(etaexpSd). Both follow the approximate-%CV convention above.
    # =========================================================================
    expSd <- 0.559; label("Residual error SD, additive on the natural-log scale (log-scale SD)") # Table 3 'Residual variability (%CV)' = 55.9 (RSE 3.99%)
    etaexpSd ~ 0.4134 # Table 3 'omega EPS' = 64.3 %CV (RSE 14.5%), inter-subject variability in residual variability; variance = 0.643^2
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual PK parameters (Table 3 footnote a).
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (BSA / 1.876)^e_bsa_cl
    vc <- exp(lvc)

    # -----------------------------------------------------------------------
    # 2. One compartment, linear elimination; dosing is a continuous IV
    # infusion into central.
    # -----------------------------------------------------------------------
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # -----------------------------------------------------------------------
    # 3. Observation (ug/L x 1000 = pg/mL) and residual error with its own
    # per-subject eta.
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc
    expSdInd <- expSd * exp(etaexpSd)
    Cc ~ lnorm(expSdInd)
  })
}
