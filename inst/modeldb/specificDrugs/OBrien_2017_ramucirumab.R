OBrien_2017_ramucirumab <- function() {
  description <- "Two-compartment population PK model for ramucirumab in patients with advanced solid tumours (O'Brien 2017)"
  reference <- "O'Brien L, Westwood P, Gao L, Heathman M. Population pharmacokinetic meta-analysis of ramucirumab in cancer patients. Br J Clin Pharmacol. 2017;83(12):2741-2751. doi:10.1111/bcp.13403"
  vignette <- "OBrien_2017_ramucirumab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Ramucirumab was assayed in SERUM by a validated ELISA
  # (Methods, 'Data': 'Ramucirumab serum concentrations were evaluated using a
  # validated enzyme-linked immunosorbent assay method'), hence specimen =
  # 'serum' rather than 'plasma'.
  compartmentData <- list(
    central = list(analyte = "ramucirumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "ramucirumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline (not time-varying) body weight. The only covariate retained in the final model.",
        "Applied as a power function on both CL and V1 normalised to the population median of 68 kg",
        "(Table 4 footnotes b and c). Cohort range 31.9-143 kg, mean 70.5 kg (CV 23%; Table 2).",
        "Note that 68 kg is the MEDIAN used for centering, while Table 2 reports the MEAN of 70.5 kg;",
        "the model reference is 68 kg. Body weight also sets the dose, which was administered on a",
        "mg/kg basis (8 mg/kg Q2W or 10 mg/kg Q3W).",
        sep = " "
      ),
      source_name = "body weight"
    )
  )

  # Covariates screened by O'Brien 2017 but NOT retained in the final model.
  # Documented here (not in covariateData) so the provenance of the paper's
  # covariate screen is preserved without declaring unused model inputs.
  # Retention required all three of: OFV drop >= 6.635 (P < 0.01, 1 df) on
  # univariate screening, >= 10% relative reduction in the IIV of the affected
  # parameter, and >= 20% influence on the parameter value; survivors then had
  # to withstand backward elimination at OFV >= 10.828 (P < 0.001).
  covariatesDataExcluded <- list(
    SEXF = list(
      description = paste(
        "Female sex indicator. Statistically significant on V1 (OFV drop 43 points, IIV reduction 14%)",
        "with females having 10% lower V1 than males, and on CL (females ~10% lower CL) but with only",
        "a 7% relative IIV reduction. Results: 'Due to the small effect and the minimal impact on",
        "ramucirumab exposure, sex was not found to be clinically relevant and therefore not retained",
        "in the final model.'",
        sep = " "
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Cohort 587/1639 female (36%), 1052/1639 male (64%); Table 2. Screened, not retained."
    ),
    ALB = list(
      description = paste(
        "Baseline serum albumin. Inclusion on CL gave a statistically significant OFV decrease but",
        "reduced the relative IIV of CL by only 5%. Discussion: a patient at the 5th percentile",
        "(28 g/L) was predicted to have 17% greater CL than a patient of the same body weight at the",
        "population median (37 g/L) -- below the 20% clinical-relevance threshold.",
        sep = " "
      ),
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort mean 37.0 g/L (CV 14%), range 16.0-64.8 g/L; Table 2. Screened, not retained."
    ),
    TUM_SLD = list(
      description = paste(
        "Baseline sum of the longest diameters of target tumour lesions, a measure of tumour burden.",
        "Weak effect on CL: OFV reduced 56 points but relative IIV lowered by only 4%. Discussion: a",
        "patient at the 95th percentile (SLD 198 mm) was predicted to have 15% greater CL than a",
        "patient of the same body weight at the median (SLD 71 mm).",
        sep = " "
      ),
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort mean 84.0 mm (CV 69%), range 10.0-438 mm; Table 3. Screened, not retained."
    ),
    CRCL = list(
      description = paste(
        "Cockcroft-Gault estimated creatinine clearance. Renal status had no significant influence on",
        "CL or V1. Predicted average steady-state concentration was graphically similar across normal,",
        "mild and moderate renal impairment (Figure 4), supporting the label statement that no dose",
        "adjustment is required for mild-to-moderate renal impairment.",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cohort mean 89.8 mL/min (CV 35%), range 25.3-303 mL/min; Table 2. Categorised as normal",
        "> 90 (42%), mild 60-90 (42%), moderate 30-60 (15%), severe 15-30 (< 1%). Screened, not retained.",
        sep = " "
      )
    ),
    LDH = list(
      description = paste(
        "Baseline lactate dehydrogenase, modelled on the log-transformed scale (Table 2 footnote a).",
        "Screened as a measure of tumour burden; did not meet the predefined retention criteria.",
        sep = " "
      ),
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort mean 298 U/L (CV 99%), range 80-5034 U/L; Table 2. Screened, not retained."
    ),
    WHO_PS = list(
      description = paste(
        "Eastern Cooperative Oncology Group (ECOG) performance-status grade, screened as a measure of",
        "tumour burden (Methods; Table 3). Did not meet the predefined retention criteria.",
        sep = " "
      ),
      units = "(grade 0-2)",
      type = "continuous",
      reference_category = NULL,
      notes = "Cohort grade 0 730 (44%), grade 1 904 (55%), grade 2 4 (< 1%); Table 3. Screened, not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1639L,
    n_studies = 11L,
    n_observations = 6427L,
    age_range = "19-87 years",
    age_mean = "60.7 years (CV 18%)",
    weight_range = "31.9-143 kg",
    weight_mean = "70.5 kg (CV 23%)",
    weight_median = "68 kg (the covariate-model reference weight; Table 4 footnotes b and c)",
    sex_female_pct = 36,
    race_ethnicity = c(White = 69, Asian = 26, Other = 5),
    disease_state = paste(
      "Advanced/metastatic cancer: colorectal (27%), nonsmall cell lung (27%), gastric (24%),",
      "hepatocellular (19%), metastatic breast (< 1%) and other solid tumours (2%).",
      sep = " "
    ),
    dose_range = "8 mg/kg IV Q2W or 10 mg/kg IV Q3W, each as an approximately 1 h infusion",
    regions = "Global (11 Phase 1b/2/3 trials, including REGARD, RAINBOW, REACH, REVEL and RAISE)",
    renal_function = "Normal 42%, mild 42%, moderate 15%, severe < 1% by Cockcroft-Gault CLcr (Table 2)",
    hepatic_function = "Normal 64%, mild 32%, moderate 1% by NCI-ODWG classification (Table 2)",
    notes = paste(
      "Baseline demographics and laboratory covariates in Table 2; tumour-burden measures in Table 3;",
      "per-study breakdown in Table 1. Concentrations below the assay LLOQ (1900 or 2500 ng/mL",
      "depending on study) were excluded from the analysis, so the model was fit to quantifiable",
      "observations only. Early clinical-development PK data (doses below 8 mg/kg) were excluded",
      "because of a bioanalytical assay change, so the model is only supported over 8-10 mg/kg.",
      sep = " "
    )
  )

  ini({
    # Structural parameters -- typical values at the reference body weight of 68 kg
    # (the population median baseline body weight; Table 4 footnotes b and c).
    lcl <- log(0.0148); label("Clearance at 68 kg (L/h)") # Table 4, row 'Clearance (CL), l h-1' = 0.0148 (1.97% SEE)
    lvc <- log(3.26); label("Central volume of distribution at 68 kg (L)") # Table 4, row 'Central volume of distribution (V1), l' = 3.26 (0.880% SEE)
    lq <- log(0.0102); label("Intercompartmental clearance (L/h)") # Table 4, row 'Intercompartmental clearance (Q), l h-1' = 0.0102 (17.8% SEE)
    lvp <- log(2.04); label("Peripheral volume of distribution (L)") # Table 4, row 'Peripheral volume of distribution (V2), l' = 2.04 (5.20% SEE)

    # Body-weight effects: power functions centred on the 68 kg median. Exponents
    # were ESTIMATED (not fixed at allometric 0.75 / 1), hence no fixed().
    e_wt_cl <- 0.499; label("Power exponent on (WT/68) for CL (unitless)") # Table 4, row 'Effect of body weight on CL' = 0.499 (7.64% SEE); footnote b: CL = 0.0148 * (body weight/68)^0.499
    e_wt_vc <- 0.556; label("Power exponent on (WT/68) for Vc (unitless)") # Table 4, row 'Effect of body weight on V1' = 0.556 (6.26% SEE); footnote c: V1 = 3.26 * (body weight/68)^0.556

    # Interpatient variability. Table 4's 'Interpatient variability' block prints
    # the four diagonals as percentages but the CL-V1 off-diagonal as a bare
    # covariance (0.0478), with no footnote defining the percentage. A covariance
    # has no 'CV' form, so it can only be on the raw OMEGA scale; if the diagonals
    # were exact log-normal CVs they would be incommensurable with it and no
    # CL-V1 correlation could be recovered from the table at all. The block is
    # therefore raw OMEGA throughout: omega = percent / 100, variance = omega^2.
    # The residual block corroborates the same reporting habit -- 'Additive' is
    # printed in concentration units (4.80 ug/mL, an SD) and 'Proportional' as
    # 22.5%, i.e. 100 * sqrt(variance), not a back-transformed CV.
    # See the vignette section 'Reading the interpatient-variability column'.
    # Exponential (log-normal) IIV per Results: 'Exponential interindividual
    # variability (IIV) terms were included for CL, V1, V2, and Q, with
    # covariance between CL and V1.'
    etalcl + etalvc ~ c(
      0.104329, # 0.323^2; Table 4, row 'Clearance (CL)' = 32.3% (5.97% SEE)
      0.0478, # Table 4, row 'Covariance (CL and V1)' = 0.0478 (8.79% SEE), as printed
      0.052441 # 0.229^2; Table 4, row 'Central volume of distribution (V1)' = 22.9% (9.33% SEE)
    )
    etalq ~ 0.6724 # 0.820^2; Table 4, row 'Intercompartmental clearance (Q)' = 82.0% (26.3% SEE)
    etalvp ~ 0.2916 # 0.540^2; Table 4, row 'Peripheral volume of distribution (V2)' = 54.0% (21.1% SEE)

    # Combined additive + proportional residual error per Results: 'Residual
    # variability was accounted for by a combined additive and proportional
    # error structure.'
    propSd <- 0.225; label("Proportional residual error (fraction)") # Table 4, row 'Proportional' = 22.5% (5.25% SEE)
    addSd <- 4.80; label("Additive residual error (ug/mL)") # Table 4, row 'Additive (ug ml-1)' = 4.80 (9.39% SEE)
  })
  model({
    # Individual parameters. Body weight is centred on the 68 kg population
    # median; Q and V2 carry no covariate effect (Table 4 lists none).
    cl <- exp(lcl + etalcl) * (WT / 68)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 68)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Linear two-compartment disposition with zero-order IV infusion (the
    # infusion is supplied by the event table via rate/dur on `central`) and
    # first-order elimination, per Results 'PPK model development'.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> mg/L = ug/mL, matching the reported units.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
