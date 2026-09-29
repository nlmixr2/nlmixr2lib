Alghamdi_2019_cycloserine <- function() {
  description <- paste(
    "One-compartment population PK model for oral cycloserine in adults",
    "treated for drug-resistant tuberculosis plus healthy volunteers",
    "(Alghamdi 2019). First-order absorption with a lag time feeds a",
    "one-compartment disposition model. Apparent clearance carries an",
    "exponential shift for the TB/NTM patient stratum relative to healthy",
    "subjects and a power effect of Cockcroft-Gault creatinine clearance;",
    "apparent volume is scaled linearly by body weight (exponent fixed to",
    "1). Between-occasion variability is carried on apparent clearance.",
    sep = " "
  )
  reference <- paste(
    "Alghamdi WA, Alsultan A, Al-Shaer MH, An G, Ahmed S, Alkabab Y, Banu S,",
    "Barbakadze K, Houpt E, Kipiani M, Mikiashvili L, Schmidt S, Heysell SK,",
    "Kempker RR, Cegielski JP, Peloquin CA. Cycloserine population",
    "pharmacokinetics and pharmacodynamics in patients with tuberculosis.",
    "Antimicrob Agents Chemother. 2019 Apr 25;63(5):e00055-19.",
    "doi:10.1128/AAC.00055-19. PMCID: PMC6496076.",
    sep = " "
  )
  vignette <- "Alghamdi_2019_cycloserine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "cycloserine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cycloserine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Alghamdi 2019 Table 2 row 'beta V, wt' = 1.00, fixed: weight enters",
        "the apparent volume of distribution as the median-normalised power",
        "function of Methods equation 2 with the exponent held at 1 (Results,",
        "'Population pharmacokinetic analysis': 'The addition of weight on the",
        "apparent volume of distribution (V/F) followed allometric scaling,",
        "with the exponent fixed to 1'). The normalising constant is the",
        "median weight of the 247-subject analysis population, 59.0 kg",
        "(Results, 'Population demographics': median 59.0 kg, IQR 51.4-68.6);",
        "the by-stratum medians of Table 1 are 58.0 kg in the 235 patients",
        "and 77.3 kg in the 12 healthy subjects. Weight was the only body-size",
        "descriptor retained; BMI was screened and dropped (see",
        "covariatesDataExcluded).",
        sep = " "
      ),
      source_name = "wt"
    ),
    CRCL = list(
      description = "Creatinine clearance, Cockcroft-Gault, raw mL/min (NOT BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Alghamdi 2019 Methods, 'Population pharmacokinetic modeling and",
        "Monte Carlo simulations': 'CrCL was calculated using the",
        "Cockcroft-Gault equation', so the values are raw mL/min and are NOT",
        "BSA-normalised; supplying a mL/min/1.73 m^2 value would silently",
        "rescale the renal term. Enters CL/F as the median-normalised power",
        "function of Methods equation 2 with beta = 0.413 (Table 2). The",
        "paper motivates the effect mechanistically: cycloserine is",
        "'approximately 70% renally cleared' (Discussion).",
        "",
        "The pooled-cohort median CrCL is not printed anywhere in the paper.",
        "The model uses 89.1 mL/min, the median of the 235-patient stratum",
        "(Table 1, IQR 68.8-111.9); the 12 healthy subjects had a median of",
        "108.8 mL/min. Patients are 95% of the analysis population, so the",
        "patient median is a close proxy for the pooled median: applying the",
        "same rank shift to the weight column recovers 58.8 kg against the",
        "printed pooled median of 59.0 kg, i.e. about 1.5% high, which at",
        "the exponent 0.413 moves CL/F by under 1%. See the vignette Errata.",
        sep = " "
      ),
      source_name = "CrCL"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator: 1 = healthy volunteer, 0 = patient with tuberculosis or nontuberculous mycobacterial disease",
      units = "(binary)",
      type = "binary",
      reference_category = NULL,
      notes = paste(
        "Alghamdi 2019 Table 2 row 'beta CL, patients (vs HS)' = -0.660",
        "(RSE 18.7%, P < 0.0001), entering CL/F through the categorical",
        "exponential model of Methods equation 1,",
        "CL = CL_POP * exp(beta * indicator). The paper's indicator is the",
        "PATIENT level with healthy subjects as reference, so CL_POP = 2.00",
        "L/h is the healthy-subject clearance and 2.00 * exp(-0.660) = 1.03",
        "L/h is the patient clearance -- exactly the two values quoted in",
        "Results ('The CL/F of cycloserine was estimated to be 2.00 liter/h",
        "in healthy subjects and 1.03 liter/h in patients').",
        "",
        "The canonical column has the OPPOSITE polarity (1 = healthy), so",
        "model() rebuilds the paper's patient indicator as (1 - DIS_HEALTHY)",
        "and applies the published -0.660 verbatim. Set DIS_HEALTHY = 0 to",
        "simulate the TB/NTM patient population the paper's Monte Carlo",
        "analysis and Table 3 are based on.",
        sep = " "
      ),
      source_name = "presence or absence of disease (healthy subjects versus patients with TB)"
    ),
    OCC = list(
      description = "Integer dosing/sampling occasion index, 1-4",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Alghamdi 2019 Table 2 row 'gamma, CL/F' = 0.190 (RSE 21.1%):",
        "interoccasion variability was estimated on apparent clearance",
        "(Methods: 'Interindividual (omega) and interoccasion variabilities",
        "(gamma) were also estimated, assuming log-normal distribution').",
        "The paper never states how many occasions were defined. The richest",
        "contributing data set is the Bangladesh cohort, sampled at 2, 4 and",
        "8 weeks after treatment initiation (Methods, 'Study data sets and",
        "subjects'), i.e. three occasions; a fourth slot is provided here as",
        "headroom and occasions 2-4 repeat occasion 1's variance, the",
        "analogue of a NONMEM $OMEGA BLOCK(1) SAME. Set OCC = 1 throughout to",
        "reproduce the paper's single-interval steady-state simulations.",
        sep = " "
      ),
      source_name = "occasion"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = paste(
        "Screened as a covariate on the PK parameters (Methods: 'Age, sex,",
        "body weight, body mass index, absence or presence of disease ...,",
        "type of disease ..., CrCL, and site also were tested as covariates",
        "on the PK parameters') but not retained in the final model. Cohort",
        "median 41.0 years, IQR 28.9-52.0 (Results, 'Population",
        "demographics'); Table 1 gives 36.1 (27.9-43.9) in healthy subjects",
        "and 41.0 (29.0-52.7) in patients. Age nonetheless matters for",
        "simulation because it is a Cockcroft-Gault input and therefore",
        "enters the model indirectly through CRCL.",
        sep = " "
      ),
      source_name = "age"
    ),
    SEXF = list(
      description = "Female-sex indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened and not retained (Methods covariate list). Table 1: 6 of",
        "12 healthy subjects (50.0%) and 179 of 235 patients (76.2%) were",
        "male, i.e. about 25% of the 247-subject analysis population was",
        "female. As with AGE, sex is a Cockcroft-Gault input and so enters",
        "indirectly through CRCL.",
        sep = " "
      ),
      source_name = "sex"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = paste(
        "Screened and not retained (Methods covariate list); total body",
        "weight was the body-size descriptor kept, on V/F. Table 1 medians",
        "25.7 kg/m^2 (IQR 23.0-28.2) in healthy subjects and 20.4 (18.4-22.9)",
        "in patients.",
        sep = " "
      ),
      source_name = "BMI"
    ),
    DIS_TB_XDR = list(
      description = "Drug-resistance stratum indicator ((pre-)XDR versus MDR tuberculosis)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "'Type of disease (i.e., DS-TB, MDR-TB, pre-XDR-TB, and XDR-TB)' was",
        "screened and not retained (Methods covariate list); only the coarser",
        "healthy-versus-patient contrast survived, and it is carried by",
        "DIS_HEALTHY. Table 1 stratum counts among the 235 patients: NTM 14",
        "(6.0%), DS-TB 16 (6.8%), RR/MDR-TB 160 (68.1%), pre-XDR-TB 36",
        "(15.3%), XDR-TB 9 (3.8%).",
        sep = " "
      ),
      source_name = "type of disease"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 247,
    n_studies = 5,
    n_observations = 1069,
    age_median = "41.0 years (IQR 28.9-52.0) overall; 36.1 (27.9-43.9) in healthy subjects and 41.0 (29.0-52.7) in patients",
    weight_median = "59.0 kg (IQR 51.4-68.6) overall; 77.3 (74.9-83.2) in healthy subjects and 58.0 (50.6-67.0) in patients",
    sex_female_pct = 25.1,
    disease_state = paste(
      "235 patients (160 RR/MDR-TB, 36 pre-XDR-TB, 9 XDR-TB, 16 DS-TB, 14",
      "nontuberculous mycobacterial disease) plus 12 healthy volunteers",
      sep = " "
    ),
    renal_function = "Median CrCL 108.8 mL/min (IQR 98.9-139.9) in healthy subjects and 89.1 (68.8-111.9) in patients; median serum creatinine 0.90 mg/dL in both groups",
    dose_range = "250-1,000 mg oral cycloserine",
    regions = "Republic of Georgia (Tbilisi), Bangladesh (Dhaka), and the United States (Arizona, Colorado, Florida, Texas)",
    notes = paste(
      "Table 1 lists baseline demographics by stratum and Results,",
      "'Population demographics' gives the pooled medians. Five pooled data",
      "sets (Methods, 'Study data sets and subjects'): 12 intensively",
      "sampled healthy subjects given a single 500 mg dose (17 samples over",
      "48 h); 69 MDR-TB patients from Tbilisi, Georgia (semirich, 4-6 weeks",
      "after treatment start); 42 MDR-TB patients from Bangladesh (semirich,",
      "at 2, 4 and 8 weeks); 54 sparse-sampled MDR-TB / NTM patients from",
      "National Jewish Health, Denver; and 70 sparse-sampled patients from",
      "three U.S. TB centres. Median (range) observed peak concentration",
      "26.5 (7.5-97.9) mg/L. Protein binding was assumed to be zero",
      "(Discussion), so plasma concentrations are treated as free drug in",
      "the paper's target-attainment analysis.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters -- Alghamdi 2019 Table 2, "Final model" column.
    # The base-model column of the same table is not extracted; per the
    # replicate-the-author's-structure policy a model-development paper
    # contributes its FINAL model only.
    #
    # Internal consistency check on the transcription: the paper quotes a
    # "relatively long half-life of 16.8 h" (Discussion). With the typical
    # patient values below, ln(2) * V/F / (CL/F) = 0.693 * 24.9 / 1.0338 =
    # 16.7 h, which reproduces it.
    # ------------------------------------------------------------------------
    ltlag <- log(0.326)
    label("Absorption lag time (h)")
    # Table 2 final model, row 'T lag (h)' = 0.326 (RSE 1.47%). Results: 'Adding
    # a parameter for a lag time resulted in a better fit during the absorption
    # phase (delta -2LL = -1,035.1)'.
    lka <- log(6.61)
    label("First-order absorption rate constant (1/h)")
    # Table 2 final model, row 'k a (h-1)' = 6.61 (RSE 17.1%).
    lvc <- log(24.9)
    label("Apparent central volume of distribution V/F at 59.0 kg (L)")
    # Table 2 final model, row 'V / F (liter)' = 24.9 (RSE 2.92%). Results:
    # 'V/F was estimated to be 24.9 liters'.
    lcl <- log(2.00)
    label("Apparent clearance CL/F in healthy subjects at CrCL 89.1 mL/min (L/h)")
    # Table 2 final model, row 'CL/ F (liter/h)' = 2.00 (RSE 11.9%). Results:
    # 'The CL/F of cycloserine was estimated to be 2.00 liter/h in healthy
    # subjects and 1.03 liter/h in patients'. This is the HEALTHY-subject
    # value; the patient shift is e_patient_cl below.

    # ------------------------------------------------------------------------
    # Covariate effects. Methods gives the two functional forms explicitly;
    # both are reproduced verbatim here from the published display equations
    # (rendered as images in the publisher's XML, recovered from the EuropePMC
    # supplementaryFiles bundle as AAC.00055-19-m0001.jpg / m0002.jpg):
    #   equation 1 (categorical): CL = CL_POP * [if sex = male, e^beta_male]
    #   equation 2 (continuous):  CL = CL_POP * (age / age_median)^beta_age
    # i.e. an exponential shift on a 0/1 indicator, and a power function of the
    # covariate normalised to its MEDIAN.
    # ------------------------------------------------------------------------
    e_wt_vc <- fixed(1.00)
    label("Power exponent on (WT / 59.0) for V/F (unitless)")
    # Table 2 final model, row 'beta V, wt' = '1.00, fixed'. Results: 'followed
    # allometric scaling, with the exponent fixed to 1'.
    e_patient_cl <- -0.660
    label("Exponential shift on CL/F for TB/NTM patients vs healthy subjects (unitless)")
    # Table 2 final model, row 'beta CL, patients (vs HS)' = -0.660 (RSE 18.7%,
    # P < 0.0001). Applied to the PATIENT indicator (1 - DIS_HEALTHY); see
    # covariateData$DIS_HEALTHY for the polarity note.
    e_crcl_cl <- 0.413
    label("Power exponent on (CRCL / 89.1) for CL/F (unitless)")
    # Table 2 final model, row 'beta CL, CrCL' = 0.413 (RSE 18.1%, P < 0.0001).

    # ------------------------------------------------------------------------
    # Interindividual variability, exponential (Methods: 'Interindividual
    # (omega) and interoccasion variabilities (gamma) were also estimated,
    # assuming log-normal distribution').
    #
    # Monolix 2018R1 -- the estimation tool named in Methods -- reports omega
    # on the STANDARD-DEVIATION scale, so the Table 2 omega rows are SDs and
    # the variances below are their squares. Two independent checks against the
    # paper's own Monte Carlo output in Table 3 confirm the SD reading and rule
    # out a variance reading:
    #   * AUC0-24h for 250 mg once daily, 259.5 (97.9) mg.h/L, is a CV of
    #     37.7%. AUC scales as 1/CL, so on the SD reading the predicted CV is
    #     sqrt(exp(0.353^2) - 1) = 36.4% before the CrCL spread is added; on a
    #     variance reading it would be sqrt(exp(0.353) - 1) = 65%.
    #   * Cmax for the same regimen, 16.4 (4.3) mg/L, is a CV of 26.2%. Cmax
    #     scales as 1/V, and V carries both omega_V and the weight term, giving
    #     sqrt(0.174^2 + 0.21^2) = 27.3% on the SD reading against 48% on a
    #     variance reading.
    # ------------------------------------------------------------------------
    etaltlag ~ 0.167281
    # Table 2 final model, row 'omega , T lag' = 0.409 (RSE 22.5%); 0.409^2.
    etalka ~ 2.3104
    # Table 2 final model, row 'omega , k a' = 1.52 (RSE 13.6%); 1.52^2.
    etalvc ~ 0.030276
    # Table 2 final model, row 'omega , V / F' = 0.174 (RSE 36.6%); 0.174^2.
    # Results: interindividual variability in V/F fell 'from 0.24 to 0.17'
    # when the covariates were added.
    etalcl ~ 0.124609
    # Table 2 final model, row 'omega , CL/ F' = 0.353 (RSE 9.29%); 0.353^2.
    # Results: interindividual variability in CL/F fell 'from 0.49 to 0.35'.

    # ------------------------------------------------------------------------
    # Interoccasion variability on CL/F. Implemented by the occasion-indicator
    # expansion rather than rxode2's `~ var | OCC` syntax, which parses but
    # cannot be simulated from an rxUi. Occasions after the first repeat the
    # same variance -- the analogue of a NONMEM $OMEGA BLOCK(1) SAME -- and are
    # therefore fixed(). See covariateData$OCC for the occasion count.
    # ------------------------------------------------------------------------
    etaiov_cl_1 ~ 0.0361
    etaiov_cl_2 ~ fixed(0.0361)
    etaiov_cl_3 ~ fixed(0.0361)
    etaiov_cl_4 ~ fixed(0.0361)
    # Table 2 final model, row 'gamma , CL/ F' = 0.190 (RSE 21.1%); 0.190^2, on
    # the same SD scale as the omega rows above.

    # ------------------------------------------------------------------------
    # Residual error. Results: 'The proportional model was selected to estimate
    # the residual error'; no additive term was retained.
    # ------------------------------------------------------------------------
    propSd <- 0.190
    label("Proportional residual error (fraction)")
    # Table 2 final model, row 'Proportional' = 0.190 (RSE 3.37%).
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Derived covariate terms.
    #
    # The paper's categorical covariate is the PATIENT level (healthy subjects
    # are the reference), while the canonical column DIS_HEALTHY runs the other
    # way. Rebuilding the patient indicator here lets the published -0.660 and
    # the published healthy-subject CL/F of 2.00 L/h both be carried verbatim.
    # ----------------------------------------------------------------------
    dis_patient <- 1 - DIS_HEALTHY

    # Interoccasion-variability term on CL/F, selected by the occasion index.
    occ1 <- (OCC == 1)
    occ2 <- (OCC == 2)
    occ3 <- (OCC == 3)
    occ4 <- (OCC == 4)
    iov_cl <- occ1 * etaiov_cl_1 + occ2 * etaiov_cl_2 +
      occ3 * etaiov_cl_3 + occ4 * etaiov_cl_4

    # ----------------------------------------------------------------------
    # 2. Individual PK parameters. Normalising constants: 59.0 kg is the median
    #    weight of the 247-subject analysis population (Results, 'Population
    #    demographics'); 89.1 mL/min is the median CrCL of the 235-patient
    #    stratum (Table 1), the closest printed proxy for the unpublished
    #    pooled median -- see covariateData$CRCL.
    # ----------------------------------------------------------------------
    tlag <- exp(ltlag + etaltlag)
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / 59.0)^e_wt_vc
    cl <- exp(lcl + etalcl + iov_cl) *
      exp(e_patient_cl * dis_patient) *
      (CRCL / 89.1)^e_crcl_cl

    # 3. Micro-constant.
    kel <- cl / vc

    # ----------------------------------------------------------------------
    # 4. ODE system. One compartment with first-order absorption (Results:
    #    'The entire data set was best described by a one-compartment model,
    #    with a first-order absorption and lag phase'). Written out explicitly
    #    rather than via linCmt() so the lag term and the ODE states stay
    #    visible.
    # ----------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 5. Absorption lag applied to the depot.
    alag(depot) <- tlag

    # 6. Observation and error.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
