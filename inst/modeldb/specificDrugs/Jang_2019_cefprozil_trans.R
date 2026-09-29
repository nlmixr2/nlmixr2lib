Jang_2019_cefprozil_trans <- function() {
  description <- "Population pharmacokinetic model for TRANS-cefprozil, the ~10% isomer of the commercial cis/trans mixture, after a single 1000 mg oral cefprozil dose in healthy adult Korean males: one compartment with first-order absorption, an absorption lag time and first-order elimination. No covariate was retained: creatinine clearance, total protein, albumin and body surface area were all screened and rejected, so this isomer's disposition is covariate-free where cis- and total cefprozil carry a creatinine-clearance effect on clearance. Unlike the cis and total models the kinetics are NOT flip-flop (Ka 0.829 1/h exceeds Kel 0.511 1/h), and the apparent volume is 2.4-fold larger. Dose the trans content of the tablet (100 mg of a 1000 mg dose); see population$notes. One of three independently fitted models in Jang 2019; see also modellib('Jang_2019_cefprozil_cis') and modellib('Jang_2019_cefprozil_total')."
  reference <- paste(
    "Jang JH, Jeong SH, Cho HY, Lee YB. (2019).",
    "Population Pharmacokinetics of Cis-, Trans-, and Total Cefprozil",
    "in Healthy Male Koreans.",
    "Pharmaceutics 11(10):531.",
    "doi:10.3390/pharmaceutics11100531.",
    sep = " "
  )
  vignette <- "Jang_2019_cefprozil"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # No covariate entered the trans-cefprozil final model. Jang 2019 Results 3.3:
  # 'no covariate affected PK parameters of trans-cefprozil', and Table 3 marks
  # the BASE model as the final model for this analyte.
  covariateData <- list()

  # Screened by Jang 2019 but NOT retained in the trans final model (Table 3
  # stepwise search plus Methods 2.3 candidate list). Documented here so the
  # covariate screen's provenance is preserved without declaring unused covariates.
  # Note that CRCL appears here, not in covariateData: it is retained for cis and
  # total cefprozil but rejected for the trans isomer.
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation",
      units = "mL/min (raw Cockcroft-Gault, NOT BSA-normalized)",
      type = "continuous",
      notes = "Table 3 'CrCl on clearance' for trans-cefprozil: dOFV +3.162 (the OFV rose, from -251.70 to -248.54), so the effect was rejected. This is the one covariate that IS retained in the cis and total models -- see modellib('Jang_2019_cefprozil_cis'). Jang 2019 Discussion attributes the difference to a possible difference in the elimination pathway of the two isomers."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; cohort 21-27 years (Table 1). Not carried into the stepwise search reported in Table 3."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate. Unlike cis and total cefprozil, 'Weight on volume' is not listed in the Table 3 trans-cefprozil block at all. Not retained."
    ),
    BSA = list(
      description = "Body surface area (Mosteller equation)",
      units = "m^2",
      type = "continuous",
      notes = "Table 3 'BSA on clearance' for trans-cefprozil: dOFV -0.676, short of the -3.84 forward-selection threshold. Not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/dL as reported by Jang 2019 Table 1 (register canonical unit is g/L)",
      type = "continuous",
      notes = "Table 3 'Total protein on clearance' for trans-cefprozil: dOFV +1.292. Not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL as reported by Jang 2019 Table 1 (register canonical unit is g/L)",
      type = "continuous",
      notes = "Table 3 'Albumin on clearance' for trans-cefprozil: dOFV -2.763, short of the -3.84 threshold. Not retained."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; Discussion reports no significant effect on any PK parameter. Not carried into Table 3."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; Discussion reports no significant effect on any PK parameter. Not carried into Table 3."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate; Discussion reports no significant effect on any PK parameter. Not carried into Table 3."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "mg/dL",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate. Not carried into Table 3."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate. Not carried into Table 3."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Methods 2.3 candidate covariate. Not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "trans-cefprozil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "trans-cefprozil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 35L,
    n_studies = 1L,
    age_range = "21-27 years",
    age_median = "24 years",
    weight_range = "53.1-91.8 kg",
    weight_median = "69.5 kg",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state = "Healthy volunteers",
    dose_range = "Single 1000 mg oral dose of cefprozil with 240 mL of water after an overnight fast (Jang 2019 Methods 2.1). For THIS model the dose amount is the TRANS CONTENT of that tablet, 100 mg at the ~9:1 cis:trans ratio; see notes.",
    regions = "Republic of Korea (Chonnam National University, Gwangju)",
    renal_function = "Normal to supranormal; CrCl 86.57-159.05 mL/min, median 124.41 mL/min (Cockcroft-Gault, Table 1). Not retained as a covariate for this analyte.",
    notes = paste(
      "35 healthy Korean males from the reference arm of a single-dose, randomized,",
      "two-way, open-label, crossover bioequivalence study (Bioequivalence Test No. 611).",
      "420 plasma samples per analyte; sampling at 0, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2, 3, 4, 8",
      "and 12 h post dose. Trans-cefprozil was assayed by UPLC-ESI-MS/MS with an LLOQ of",
      "15 ng/mL and a calibrated range of 0.015-2.5 ug/mL.",
      "DOSE BASIS: Jang 2019 never states the dose amount entered into each analyte's dataset.",
      "It is back-solved here from the paper's own reported noncompartmental results",
      "(Discussion, p. 13) and is the isomer content of the 1000 mg dose at the ~9:1",
      "cis:trans ratio quoted in the Introduction: 100 mg for trans, 900 mg for cis,",
      "1000 mg for total. Inverting Dose = AUC * CL on the three reported mean NCA AUCs",
      "gives 1022 mg (cis), 113 mg (trans) and 1137 mg (total) -- a back-solved cis:trans",
      "ratio of 9.05:1, and a cis + trans sum (1135 mg) that matches the total (1137 mg)",
      "to 0.2%. The competing reading that all three datasets were dosed at 1000 mg is",
      "falsified for THIS model in particular: it would predict a trans AUC of",
      "56.5 ug*h/mL against the 6.38 ug*h/mL reported, a 9-fold error.",
      "See the vignette Errata for the full arithmetic.",
      "PPK analysis was run in Phoenix NLME 8.1 (FOCE with extended least squares),",
      "not NONMEM."
    )
  )

  ini({
    # Structural parameters -- Jang 2019 Table 4, 'Trans-cefprozil / Final model'.
    # Volume and clearance are printed in mL and mL/h; they are divided by 1000
    # here so that a dose in mg gives a concentration in mg/L = ug/mL, the
    # concentration unit the paper reports.
    lka <- log(0.829); label("First-order absorption rate constant Ka (1/h)") # Table 4, Trans-cefprozil Final model: tvKa = 0.829 1/h (RSE 10.98%)
    lvc <- log(34617.50 / 1000); label("Apparent central volume of distribution V/F (L)") # Table 4, Trans-cefprozil Final model: tvV = 34,617.50 mL (RSE 5.77%)
    lcl <- log(17701.80 / 1000); label("Apparent clearance CL/F (L/h)") # Table 4, Trans-cefprozil Final model: tvCl = 17,701.80 mL/h (RSE 3.38%)
    ltlag <- log(0.352); label("Absorption lag time (h)") # Table 4, Trans-cefprozil Final model: tvTlag = 0.352 h (RSE 5.09%)

    # IIV. Jang 2019 reports omega^2 (variances) directly, so these are entered
    # as variances without transformation. The trans-cefprozil IIV structure is
    # the MIRROR of the cis / total structure: Table 2 selects step 02-03-01
    # 'Remove IIV V' here, so there is IIV on Ka but NONE on V, whereas the cis
    # and total models have IIV on V but none on Ka.
    etalka ~ 0.076 # Table 4, Trans-cefprozil Final model: omega^2 ka = 0.076 (shrinkage 10.95%)
    etalcl ~ 0.017 # Table 4, Trans-cefprozil Final model: omega^2 Cl = 0.017 (shrinkage 9.28%)
    etaltlag ~ 0.094 # Table 4, Trans-cefprozil Final model: omega^2 Tlag = 0.094 (shrinkage 12.06%); the Discussion quotes omega = 36.1% for this term, which does not equal sqrt(0.094) = 30.7% -- see vignette Errata

    # Residual error. Table 2 step 02-01 'Proportional' is the SELECTED residual
    # model for trans-cefprozil, unlike the log-additive model chosen for the cis
    # and total analytes. Jang 2019 Discussion: 'Proportional residual variability
    # was 23.2% for trans-cefprozil'.
    propSd <- 0.232; label("Proportional residual error (fraction)") # Table 4, Trans-cefprozil Final model: sigma = 0.232 (RSE 6.61%)
  })

  model({
    # Individual parameters. No covariate effects: the trans-cefprozil final
    # model is the base model (Jang 2019 Table 3).
    ka <- exp(lka + etalka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
