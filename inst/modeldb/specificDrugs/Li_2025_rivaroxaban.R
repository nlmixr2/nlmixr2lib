Li_2025_rivaroxaban <- function() {
  description <- "One-compartment population PK model for oral rivaroxaban in elderly (65 years and older) Chinese patients with non-valvular atrial fibrillation, with a creatinine-clearance power effect and a concomitant-amiodarone effect on CL/F (Li 2025)"
  reference <- paste(
    "Li TY, Liu X, Tang L, Cai W, Zhang J, Su Y, Yu C, Liu BJ, Jiang LF, Liu C,",
    "Yu YX. Dose Optimization of Rivaroxaban in Elderly Chinese Patients with",
    "Non-Valvular Atrial Fibrillation: Analysis of Populations Pharmacokinetics",
    "and Exposure-Response Relationships. Drug Des Devel Ther. 2025;19:10703.",
    "doi:10.2147/DDDT.S563521.",
    sep = " "
  )
  vignette <- "Li_2025_rivaroxaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Rivaroxaban 10, 15 or 20 mg orally once daily (Methods 'Study Design'),
  # measured in plasma by a validated LC-MS/MS method (Methods 'Rivaroxaban
  # Quantification').
  compartmentData <- list(
    depot = list(analyte = "rivaroxaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rivaroxaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated with the Cockcroft-Gault equation",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min -- NOT BSA-normalized to",
        "mL/min/1.73 m^2 (Methods 'Study Design': 'estimated creatinine clearance",
        "[CrCL] ... calculated using the Cockcroft-Gault equation'). Enters CL/F as a",
        "power of the ratio to 47.4 mL/min, the reference value printed in the Table 2",
        "footnote equation 'CL (L/h) = 4.5 * (CrCL/47.4)^0.18 * exp(amiodarone *",
        "-0.12)*exp(etaCL)'. Table 1 gives the cohort median as 48.6 mL/min (range",
        "14.3-91.7); the printed equation's 47.4 is used as the centering value.",
        "73.1% of patients had CrCL < 50 mL/min (Results 'Population",
        "Characteristics'); patients with CrCL < 15 mL/min were excluded."
      ),
      source_name = "CrCL"
    ),
    CONMED_AMIO = list(
      description = "Concomitant amiodarone coadministration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant amiodarone)",
      notes = paste(
        "1 = the patient was taking amiodarone with rivaroxaban, 0 = not. 45 of 93",
        "patients (48.4%) were positive (Table 1); 'The standard dosing regimen for",
        "amiodarone was 200 mg qd for long-term rhythm control' (Results 'Population",
        "Characteristics'). Treated as time-fixed: every patient had been on a stable",
        "rivaroxaban regimen for at least 4 weeks before sampling. Enters CL/F as",
        "exp(-0.124 * CONMED_AMIO), i.e. CL/F 0.883-fold with amiodarone (Discussion:",
        "'reduction in CL to 0.88-fold with concomitant amiodarone use')."
      ),
      source_name = "amiodarone"
    )
  )

  # Screened in the stepwise covariate search (Methods 'Population
  # Pharmacokinetic Analysis and Model Evaluation') but NOT retained: 'none of the
  # other covariates were found to significantly affect the deposition of
  # rivaroxaban' (Results 'PopPK Model of Rivaroxaban'). The paper names these
  # covariates as examples ('eg') and reports no estimate for any of them.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened ('gender'); not retained. 43 of 93 patients female (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Mean 79 +/- 7 years (Table 1), range 66-94 years (Discussion). Age is an input to the retained Cockcroft-Gault CRCL."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained. Mean 61.4 +/- 12.1 kg (Table 1). Weight is an input to the retained Cockcroft-Gault CRCL."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened; not retained. Median 16 umol/L (range 6-40, Table 1)."
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus comorbidity indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diabetes)",
      notes = "Screened ('comorbidities (eg, diabetes, hyperlipidemia)'); not retained. 27 of 93 patients (29.0%, Table 1)."
    ),
    CONMED_ATORVASTATIN = list(
      description = "Concomitant atorvastatin indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant atorvastatin)",
      notes = "Screened ('concomitant medications (eg, amiodarone, atorvastatin)'); not retained. 37 of 93 patients (39.8%, Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 93,
    n_studies = 1,
    n_observations = 256,
    age_range = "66-94 years (Discussion); 65 years and older by eligibility",
    age_median = "not reported; mean 79 +/- 7 years",
    weight_range = "not reported; mean 61.4 +/- 12.1 kg",
    weight_median = "not reported; mean 61.4 kg",
    sex_female_pct = 46.2,
    race_ethnicity = c(Asian = 100),
    disease_state = "non-valvular atrial fibrillation, on rivaroxaban for at least 4 weeks",
    dose_range = "10, 15 or 20 mg orally once daily",
    regions = "China (single center: Suzhou Municipal Hospital, Nanjing Medical University)",
    renal_function = "Cockcroft-Gault CrCL median 48.6 mL/min, range 14.3-91.7 mL/min; 73.1% below 50 mL/min. CrCL < 15 mL/min excluded.",
    hepatic_function = "Child-Pugh B or C excluded",
    co_medication = "Amiodarone 45 (48.4%), atorvastatin 37 (39.8%), rosuvastatin 20 (21.5%), PPI 23 (24.7%), aspirin 12 (12.9%), clopidogrel 4 (4.3%), NSAIDs 3 (3.2%). Strong CYP3A4 or P-gp inhibitors excluded.",
    notes = paste(
      "Prospective single-center observational study, January 2021 to March 2022",
      "(Methods 'Study Design'); Table 1 baseline demographics. 50 male / 43",
      "female; height 162.3 +/- 8.2 cm; BMI 23.2 +/- 3.7 kg/m^2; albumin 39.64 +/-",
      "5.34 g/L. Median CHA2DS2-VASc 5, median HAS-BLED 2. AF type: paroxysmal",
      "40.9%, persistent 26.9%, permanent 32.3%. Comorbidities: hypertension 86.0%,",
      "congestive heart failure 63.4%, hyperlipoproteinemia 39.8%, diabetes 29.0%,",
      "stroke 18.3%. Sampling at steady state (at least five consecutive doses):",
      "a trough 0.5 h before the next dose and a peak 3 +/- 1 h post-dose, 256",
      "concentrations in total. Fitted in Phoenix NLME 8.3 with FOCE-ELS."
    )
  )

  ini({
    # Structural parameters. Dose in mg, CL/F in L/h, V/F in L.
    # ka was not estimated: 'Due to limited plasma concentration data points
    # during the absorption phase, the absorption rate constant (Ka) could not be
    # reliably estimated and was therefore fixed to the reported value
    # (0.821 h-1)' (Results 'PopPK Model of Rivaroxaban', citing Li 2025
    # references 28 and 29).
    lka <- fixed(log(0.821)); label("Absorption rate constant (ka, 1/h), from the literature (Li 2025 references 28-29)") # Table 2 final model 'tvKa 0.821, fixed'

    # Typical values at the covariate reference: CrCL = 47.4 mL/min, no amiodarone.
    lcl <- log(4.501); label("Apparent clearance (CL/F, L/h) at CrCL 47.4 mL/min without amiodarone") # Table 2 final model 'tvCL (L/h) 4.501 (2.84)'
    lvc <- log(30.714); label("Apparent volume of distribution (V/F, L)") # Table 2 final model 'tvV(L) 30.714 (2.19)'

    # Covariate effects on CL/F, in the form of the Table 2 footnote equation
    #   CL = tvCL * (CrCL/47.4)^dCLdCrCL * exp(dCLamiodarone * amiodarone) * exp(etaCL)
    # The footnote prints the rounded values (4.5, 0.18, -0.12); the Table 2
    # estimate column carries three decimals, used here.
    e_crcl_cl <- 0.179; label("Creatinine-clearance power exponent on CL/F (unitless)") # Table 2 final model 'dCLdCrCL 0.179 (18.61)'
    e_conmed_amio_cl <- -0.124; label("Log-scale shift in CL/F with concomitant amiodarone (unitless)") # Table 2 final model 'dCLamiodarone -0.124 (30.51)'

    # IIV, exponential on CL/F and V/F (Results: 'An exponential model ... used to
    # describe the interindividual variability'). Table 2 rows 'omega2 CL (%)' and
    # 'omega2 V (%)' print 3.6 and 3.3; read as variances 0.036 and 0.033 (the
    # percent sign marks a variance x 100). The paper's own Monte Carlo
    # simulation confirms this reading: Figure 3 boxplots within a single CrCL
    # band have an interquartile ratio of about 1.3, which matches
    # exp(2 * 0.674 * sqrt(0.036)) = 1.29 and not the 1.05 a 3.6% CV would give.
    # No eta correlation was reported.
    etalcl ~ 0.036 # Table 2 final model 'omega2 CL (%) 3.6 [shrinkage 21]'
    etalvc ~ 0.033 # Table 2 final model 'omega2 V (%) 3.3 [shrinkage 19]'

    # Residual error: additive ('An exponential model and an additive model were
    # used to describe the interindividual variability and residual variability,
    # respectively'). Phoenix 'stdev0' is the additive SD in ng/mL.
    addSd <- 18.025; label("Additive residual error (ng/mL)") # Table 2 final model 'stdev0 18.025 (7.54)'
  })
  model({
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (CRCL / 47.4)^e_crcl_cl *
      exp(e_conmed_amio_cl * CONMED_AMIO)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and vc in L give mg/L; x1000 converts to ng/mL.
    Cc <- central / vc * 1000
    Cc ~ add(addSd)
  })
}
