Jin_2026_colistinSulfate <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous colistin sulfate in",
    "critically ill adults with carbapenem-resistant Gram-negative pneumonia",
    "(Jin 2026; n = 10 Chinese ICU patients, 130 plasma concentrations).",
    "Linear elimination from the central compartment with an intravenous",
    "infusion input. Cockcroft-Gault creatinine clearance enters clearance as",
    "a power function centred on the cohort median of 116.3 mL/min (exponent",
    "0.39), so clearance rises with renal function. Inter-individual",
    "variability is exponential on all four disposition parameters, with a",
    "correlation of 0.757 between the central volume and clearance carried as",
    "a non-diagonal (block) random effect; residual error is proportional.",
    "IMPORTANT: the model was fitted to UNBOUND (free) colistin",
    "concentrations measured by LC-MS/MS after the fifth dose, while doses",
    "are TOTAL drug. Cc is therefore the free plasma concentration and CL / V",
    "are apparent values on the unbound scale -- they are not comparable with",
    "total-drug parameters without dividing by the unbound fraction (the",
    "paper assumes 0.5 for its PK/PD work but does NOT apply it in the PK",
    "model). Colistin sulfate is the active drug and must not be confused",
    "with colistimethate sodium (CMS), the inactive prodrug modelled in",
    "Plachouras 2009, Mohamed 2012, Jacobs 2016 and Karaiskos 2015. Dose",
    "units: the paper expresses every dose in million international units",
    "(MIU) and never states an MIU-to-mg potency; this model takes dose in",
    "mg. Jin's own Table 4 back-solves to about 44.6 mg/MIU, consistent with",
    "the 44 mg/MU used by modellib('Huang_2025_colistinSulfate') and",
    "modellib('Ma_2026_colistinSulfate'). See notes and the vignette Errata."
  )
  reference <- paste(
    "Jin X, Zhao D, Yang J, Yang S, Hu C (2026).",
    "Population pharmacokinetics of intravenous colistin sulfate in",
    "critically ill patients with pneumonia.",
    "Infection and Drug Resistance 19:611711.",
    "doi:10.2147/IDR.S611711.",
    sep = " "
  )
  vignette <- "Jin_2026_colistinSulfate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance computed by the Cockcroft-Gault method, RAW and",
        "NOT BSA-normalised. Methods Equation (unnumbered, page 3):",
        "CrCL (mL/min) = (140 - age in years) * body weight (kg) /",
        "(0.814 * serum creatinine in umol/L), multiplied by 0.85 for female",
        "subjects. The 0.814 divisor is the standard unit conversion",
        "that makes the Cockcroft-Gault equation accept serum creatinine in",
        "umol/L rather than mg/dL. The weight entering the numerator is",
        "BMI-dependent (Methods, immediately below the equation): actual body",
        "weight if BMI < 18.5 kg/m2, ideal body weight (IBW) if BMI is",
        "18.5-24.9 kg/m2, and adjusted body weight (ABW) if BMI >= 25 kg/m2,",
        "where IBW (kg) = Constant + 0.91 * (height in cm - 152.4) with",
        "Constant = 50 for men and 45.5 for women, and",
        "ABW (kg) = IBW + 0.4 * (actual body weight - IBW). Note that the",
        "typeset equation itself names 'ideal body weight' in its numerator",
        "while the surrounding prose specifies the three-way BMI rule; the",
        "prose is the operative description because it is the one that",
        "defines all three cases."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL centred on 116.3 mL/min, stated in the Results",
        "sentence that introduces the final model ('where 116.3 mL/min was",
        "the median of CrCL'). This is the cohort MEDIAN and is distinct from",
        "the Table 1 cohort MEAN of 103.28 +/- 41.96 mL/min -- do not",
        "substitute the mean. CrCL was the only covariate tested in the final",
        "model: the Methods prespecify that no stepwise covariate search",
        "would be run below 20 patients, so with n = 10 the model",
        "'will exclusively incorporate the highly correlated covariate,",
        "CrCL'. The 22 screened-but-never-searched laboratory and",
        "demographic candidates are recorded in covariatesDataExcluded.",
        "The simulation range used by the paper is 10-120 mL/min (Table 4).",
        "Because this column is raw mL/min and NOT normalised to",
        "1.73 m2 -- the default convention for the CRCL canonical -- supplying",
        "a BSA-normalised value to this model silently rescales clearance by",
        "(BSA / 1.73)^0.39. Same raw-mL/min convention as",
        "modellib('Sun_2025_colistinSulfate') and",
        "modellib('Delattre_2010_amikacin')."
      ),
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    # Methods 'Population Pharmacokinetic Analysis' lists 22 candidate
    # covariates, then states that a stepwise covariate model would NOT be run
    # because the sample size is below 20 patients. None of these was
    # therefore screened or retained; they are recorded here for provenance
    # only and are not referenced in model().
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened (n = 10 < 20). Cohort 48.80 +/- 18.37 years (Table 1)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened. Used inside the CrCL derivation to choose actual vs ideal vs adjusted body weight, but not as a covariate in its own right."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened. Cohort 65.60 +/- 50.59 U/L (Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened. Cohort 33.71 +/- 3.62 g/L (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened as a covariate in its own right; it enters the model only through the CrCL derivation. Cohort 93.50 +/- 102.17 umol/L (Table 1)."
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units = "mg/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened. Retained instead of CrCL by the CVVHDF sibling model modellib('Huang_2025_colistinSulfate')."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened."
    ),
    WBC = list(
      description = "White blood cell count",
      units = "10^9/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened."
    ),
    PLT = list(
      description = "Platelet count",
      units = "10^9/L",
      type = "continuous",
      notes = "Listed as a candidate covariate but never screened."
    )
    # The remaining named candidates -- hemoglobin, neutrophil percentage,
    # neutrophil count, globulin, prothrombin time, activated partial
    # thromboplastin time, fibrinogen, d-dimer, lactate, procalcitonin and
    # interleukin-6 -- have no canonical column registered and, like the ten
    # above, were never screened. They are listed in the Methods text and in
    # the vignette's Assumptions section rather than given placeholder names
    # here.
  )

  compartmentData <- list(
    # Methods 'Colistin Sulfate Administration, Sample Collection and
    # Concentration Deamination': plasma was separated by centrifugation and
    # 'The free drug concentrations were determined using a validated
    # high-performance LC-MS method'. Table 4 labels its outputs fAUC24h and
    # Cavg,free, confirming the modelled analyte is the UNBOUND fraction.
    # Colistin A and colistin B were also assayed separately (Figure 1,
    # Table 2) but the PK model was fitted to total colistin sulfate, i.e.
    # their sum.
    central = list(
      analyte = "colistin sulfate, unbound (sum of colistin A and colistin B)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "colistin sulfate, unbound (sum of colistin A and colistin B)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10,
    n_studies = 1,
    n_observations = 130,
    age_mean = "48.80 +/- 18.37 years",
    weight_mean = "64.55 +/- 11.01 kg",
    sex_female_pct = 20,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Critically ill adults (>= 18 years) in a medical ICU with pulmonary",
      "infection caused by carbapenem-resistant Gram-negative bacteria, all",
      "treated with intravenous colistin sulfate for at least 72 h. Severe",
      "illness: APACHE II 17 (12, 35), SOFA 10 (6, 18), 60% on vasoactive",
      "agents, mechanical ventilation 11.5 (6, 40) days. All 10 had pulmonary",
      "involvement and 6 (60%) had two or more infection sites (bloodstream",
      "20%, urinary tract 20%, abdomen 10%). Pathogens were",
      "carbapenem-resistant Acinetobacter baumannii 90%, Klebsiella",
      "pneumoniae 10% and Pseudomonas aeruginosa 10%. Comorbidities:",
      "hypertension, coronary heart disease and malignancy, 20% each. Length",
      "of ICU stay 24.40 +/- 11.16 days, hospital stay 40.00 +/- 20.53 days."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 103.28 +/- 41.96 mL/min (mean",
      "+/- SD, Table 1), cohort median 116.3 mL/min (the model's centring",
      "constant). Baseline serum creatinine 93.50 +/- 102.17 umol/L. Renal",
      "function was therefore largely preserved-to-augmented in this cohort;",
      "the paper's simulations extrapolate down to 10 mL/min, which is well",
      "outside the observed range and is the main reason the CrCL exponent",
      "should not be trusted far below the observed values. NO nephrotoxicity",
      "events occurred during the study."
    ),
    dose_range = paste(
      "Intravenous colistin sulfate (Shanghai SPH New Asia Pharmaceutical",
      "Co. Ltd.). Table 1 daily dose 1.75 (1.50, 2.00) MIU; the package",
      "insert recommends 1.0-1.5 MIU maintenance given two or three times",
      "daily, with a 1.0-1.5 MIU loading dose used in local practice.",
      "Treatment lasted 14.30 +/- 7.28 days. Five subjects (50%) also",
      "received adjunctive nebulised colistin sulfate, which is NOT",
      "represented in this model. The paper's simulations use a 2 h infusion.",
      "IMPORTANT: doses are stated only in MIU; this model takes mg (see",
      "notes)."
    ),
    sampling = paste(
      "Rich sampling across one full dosing interval at steady state: a",
      "pre-dose sample immediately before the FIFTH dose (0 h), then 0.5, 1,",
      "2, 3, 4, 5, 6, 7, 8, 9, 10 and 11 h after the END of the infusion --",
      "13 nominal times per subject, 130 concentrations in total. All data",
      "come from this single interval, which the Discussion flags as a",
      "limitation ('model applicability to other treatment periods remains",
      "unknown'). Observed Cmin range 0.195-1.099 mg/L."
    ),
    regions = "People's Republic of China (single centre; West China Hospital of Sichuan University, Chengdu).",
    notes = paste(
      "Baseline demographics from Jin 2026 Table 1. Single-centre prospective",
      "study, ethics approval No. [2023]1220. The final model was fitted in",
      "Phoenix NLME 8.1.0 by first-order conditional estimation with extended",
      "least squares (FOCE-ELS) and evaluated with a 1000-sample bootstrap",
      "(Table 3), goodness-of-fit plots (Figure 2) and a",
      "prediction-corrected VPC of 1000 simulations (Figure 3). The",
      "two-compartment model was selected over one-compartment on OFV",
      "(-128.93 vs -75.79), AIC (-110.93 vs -65.79) and BIC (-85.13 vs",
      "-51.45).",
      "",
      "SOURCE RECOVERY. The published PDF of this article DROPS its own",
      "Equations 1-4 during typesetting: page 6 ends the Results paragraph",
      "with 'the final population PK model was shown in Equations 1 to 4,",
      "where 116.3 mL/min was the median of CrCL' and then goes straight to",
      "Figure 1, with no equations rendered anywhere in the file (the Methods",
      "equations on page 3 render normally, so this is specific to Equations",
      "1-4, not a general math-rendering failure). The four equations were",
      "recovered from the EuropePMC JATS full text for PMC13281950, whose",
      "<disp-formula> elements carry the publisher's LaTeX source:",
      "V (L) = 18.15 * exp(eta V); V2 (L) = 10.82 * exp(eta V2);",
      "CL (L/h) = 3.20 * (CrCL/116.3)^0.39 * exp(eta CL);",
      "CL2 (L/h) = 3.23 * exp(eta CL2).",
      "Every coefficient in them matches Table 3 independently, so the",
      "equations add the covariate FORM (a power function on CrCL centred at",
      "116.3) and the IIV FORM (exponential on all four parameters) that",
      "Table 3 alone cannot supply.",
      "",
      "DOSE UNITS. The paper's dosing is entirely in MIU and no MIU-to-mg",
      "potency appears anywhere in the article, while its parameters (CL in",
      "L/h, V in L) and observations (mg/L) are on a mass basis. The",
      "conversion is therefore back-solved from the paper's own output:",
      "Table 4's simulated steady-state fAUC24h of 13.78 mg*h/L for 0.5 MIU",
      "Q12h at CrCL 120 mL/min, against the model's typical clearance of",
      "3.2393 L/h at that CrCL, gives 1 MIU = 13.78 * 3.2393 = 44.6 mg. That",
      "row is the right anchor because it is the only CrCL level at which day",
      "4 (72-96 h) is unambiguously steady state -- the terminal half-life is",
      "7.3 h at CrCL 120 but 17.3 h at CrCL 10, and the implied mg/MIU falls",
      "monotonically (44.6, 44.3, 43.8, 43.0, 40.4 at CrCL 120, 80, 50, 30,",
      "10) exactly as incomplete accumulation predicts. Two independent",
      "routes agree: modellib('Huang_2025_colistinSulfate') and",
      "modellib('Ma_2026_colistinSulfate') both use 44 mg per 10^6 IU from",
      "unrelated cohorts, and modellib('Sun_2025_colistinSulfate')",
      "back-solves 45-46 mg/MU. All of these are consistent with the potency",
      "of the colistin SULFATE salt (~22,000 IU/mg) and NOT with colistin",
      "BASE activity (30,000 IU/mg, i.e. 33.3 mg/MU); the two differ by ~35%.",
      "The conversion is deliberately NOT encoded as a model parameter -- it",
      "is a property of the marketed product, not of the pharmacokinetics --",
      "so a user supplying mg-denominated doses is unaffected by any residual",
      "uncertainty in it. The vignette derives it and validates it against all",
      "50 rows of Table 4.",
      "",
      "UNBOUND vs TOTAL. The assay measured free drug, so this model predicts",
      "UNBOUND concentrations from TOTAL doses. The unbound fraction of 0.5",
      "quoted in the PK/PD section is a literature assumption used only for",
      "the fAUC/MIC arithmetic, not a fitted parameter, and is not encoded",
      "here. Note the paper's own internal tension on this point: the Methods",
      "state 'As no specific procedure was performed to release the",
      "bounded-drug from plasma protein, we treated the observed AUC (AUCobs)",
      "as fAUC' and then separately set f = 0.5. Table 4 is labelled fAUC24h",
      "and Cavg,free throughout, and its values reconcile with the model",
      "WITHOUT any further 0.5 factor, so no unbound-fraction scaling is",
      "applied in this model.",
      "",
      "COVARIATE SIGNIFICANCE. CrCL was retained on CL despite NOT meeting",
      "the paper's own prespecified forward-inclusion threshold: Results",
      "report dOFV = -2.32 (p = 0.1281) against a stated criterion of",
      "dOFV > 3.84 (p < 0.05). The Methods justify this in advance -- with",
      "fewer than 20 patients no stepwise search is run and CrCL is included",
      "by prior knowledge rather than by test. The V-CL correlation was",
      "likewise added at dOFV = -11.30 (p = 0.0795). Both are reproduced here",
      "as the authors specified them; a user should treat the exponent 0.39",
      "(95% CI 0.29-0.50, but bootstrap 95% CI -0.94 to 1.20) as weakly",
      "identified. The bootstrap interval spanning zero is the sharpest",
      "single statement of that uncertainty in the paper.",
      "",
      "The nebulised colistin sulfate received by 50% of subjects, and the",
      "PK/PD exposure analysis (Table 5, median fAUC/MIC 39.56 against a",
      "combination MIC of 0.5 mg/L), are outside the scope of the PK model."
    )
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters -- Jin 2026 Table 3 'Final Model / Estimate',
    # each confirmed independently by the paper's Equations 1-4 as recovered
    # from the EuropePMC JATS <disp-formula> LaTeX for PMC13281950 (the
    # published PDF drops those equations; see population$notes 'SOURCE
    # RECOVERY').
    #
    # The CV (%) column of Table 3 is the RELATIVE STANDARD ERROR of the
    # estimate, not an inter-individual CV. Confirmed arithmetically on the
    # first row: 18.15 * (1 +/- 1.96 * 0.2005) = 10.02 to 26.28 against the
    # printed 95% CI of 10.94-25.36, and on the third row:
    # 3.20 * (1 +/- 1.96 * 0.1016) = 2.56 to 3.84 against the printed
    # 2.55-3.84. The IIV terms are carried separately below.
    # ----------------------------------------------------------------------
    lvc <- log(18.15); label("Central volume of distribution, V (L)")                                        # Table 3 tvV  = 18.15 L   (RSE 20.05%, 95% CI 10.94-25.36; bootstrap mean 16.86, 95% CI 9.24-31.65); Eq. 1
    lvp <- log(10.82); label("Peripheral volume of distribution, V2 (L)")                                    # Table 3 tvV2 = 10.82 L   (RSE  9.51%, 95% CI  8.78-12.86; bootstrap mean 12.20, 95% CI 8.71-17.04); Eq. 2
    lcl <- log(3.20);  label("Clearance at the reference creatinine clearance of 116.3 mL/min (L/h)")        # Table 3 tvCl = 3.20 L/h  (RSE 10.16%, 95% CI  2.55-3.84;  bootstrap mean 3.19,  95% CI 2.47-4.51);  Eq. 3
    lq  <- log(3.23);  label("Intercompartmental clearance, CL2 (L/h)")                                      # Table 3 tvCl2 = 3.23 L/h (RSE 27.03%, 95% CI  1.50-4.96;  bootstrap mean 3.30,  95% CI 0.79-7.44);  Eq. 4. Table 3 mislabels the unit of this row as '(L)'; CL2 is an intercompartmental CLEARANCE and Eq. 4 prints 'CL2 (L/h)'.

    # ----------------------------------------------------------------------
    # Covariate effect -- the ONLY covariate in the final model.
    #
    # Equation 3 prints the form unambiguously as a POWER function centred on
    # the cohort median CrCL:
    #     CL (L/h) = 3.20 * (CrCL/116.3)^0.39 * exp(eta CL)
    # The centring constant is given in the Results sentence that introduces
    # the equations ('where 116.3 mL/min was the median of CrCL'), not in
    # Table 3.
    # ----------------------------------------------------------------------
    e_crcl_cl <- 0.39; label("Power exponent on (CRCL / 116.3 mL/min) for CL (unitless)")                    # Table 3 dCldCrCL = 0.39 (RSE 13.18%, 95% CI 0.29-0.50; bootstrap mean 0.39, 95% CI -0.94 to 1.20 -- the bootstrap interval spans zero); Eq. 3

    # ----------------------------------------------------------------------
    # Inter-individual variability -- EXPONENTIAL on all four disposition
    # parameters (Equations 1-4 each carry an explicit exp(eta ...) factor;
    # Results: 'Interindividual variability was represented by an exponential
    # model').
    #
    # SCALE: Table 3 labels these rows 'omega^2 V', 'omega^2 CL', 'omega^2
    # V2', 'omega^2 CL2' -- they are VARIANCES on the log scale, which is
    # what nlmixr2's ini() expects, so no squaring is applied. Reading 0.807
    # as an SD instead would imply a 95% CV on V, against the 111% CV that
    # the variance reading gives; the variance reading is corroborated by the
    # Table 2 non-compartmental Vd,ss of 35.80 +/- 31.25 L, an 87% observed
    # CV on volume that a 95%-CV model could not generate together with the
    # residual error.
    #
    # The V-CL block: Results report 'a significant correlation was observed
    # between CL and V, which was subsequently integrated into non-diagonal
    # random effects (dOFV = -11.30, p = 0.0795). The correlation between V
    # and CL was 0.757.' The off-diagonal below is that correlation converted
    # to a covariance, 0.757 * sqrt(0.807 * 0.151) = 0.264254. The resulting
    # 2x2 block has determinant 0.807*0.151 - 0.264254^2 = 0.0520 > 0, so it
    # is comfortably positive definite and needs no nudge.
    # ----------------------------------------------------------------------
    etalvc + etalcl ~ c(0.807,
                        0.264254, 0.151)  # Table 3 'omega 2 V' = 0.807 (RSE 51.05%, eta-shrinkage 3.88%; bootstrap 0.875) and 'omega 2 CL' = 0.151 (RSE 37.09%, eta-shrinkage 0.94%; bootstrap 0.141), with off-diagonal = 0.757 * sqrt(0.807 * 0.151) from the Results correlation of 0.757
    etalvp ~ 0.025                        # Table 3 'omega 2 V2'  = 0.025 (RSE  1.68%, eta-shrinkage  5.15%; bootstrap 0.012)
    etalq  ~ 0.208                        # Table 3 'omega 2 CL2' = 0.208 (RSE 12.50%, eta-shrinkage 29.65%; bootstrap 0.878 -- the bootstrap mean is 4x the final estimate and the eta-shrinkage is the highest in the model, so this term is poorly identified)

    # ----------------------------------------------------------------------
    # Residual error -- PROPORTIONAL.
    #
    # Results: 'residual variability was characterized by a proportional
    # error model'. Table 3 reports it under 'Residual variability (sigma)'
    # as 'stdev0' = 0.1374, which is Phoenix NLME's name for the standard
    # deviation of the first residual-error term. In Phoenix a proportional
    # error is observe(CObs = C * (1 + CEps)) with error(CEps = stdev0), so
    # 0.1374 is directly the proportional SD (13.74%) and maps one-to-one
    # onto nlmixr2's prop().
    #
    # The 95% CI is on the same scale (0.0876-0.1872, i.e. 0.1374 * (1 +/-
    # 1.96 * 0.1830) = 0.0881-0.1867), which confirms the printed CV (%) of
    # 18.30 is again an RSE and that 0.1374 is an SD rather than a variance.
    # ----------------------------------------------------------------------
    propSd <- 0.1374; label("Proportional residual error (fraction)")                                        # Table 3 stdev0 = 0.1374 (RSE 18.30%, 95% CI 0.0876-0.1872; bootstrap mean 0.1318, 95% CI 0.0862-0.1841)
  })

  model({
    # 1. Individual parameters. CrCL enters CL only; V, V2 and CL2 carry no
    #    covariate in the final model (Equations 1, 2 and 4).
    cl <- exp(lcl + etalcl) * (CRCL / 116.3)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system. Colistin sulfate is given as an intravenous infusion
    #    directly into `central` (the paper's simulations use a 2 h infusion).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # 4. Observation and error. Cc is the UNBOUND plasma concentration in
    #    mg/L -- the assay measured free drug (see compartmentData).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
