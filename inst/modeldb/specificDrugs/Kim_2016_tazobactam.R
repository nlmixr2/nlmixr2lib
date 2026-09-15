Kim_2016_tazobactam <- function() {
  description <- paste(
    "Two-compartment population PK model for tazobactam in 33 Korean adult",
    "inpatients with acute infections (Kim 2016), sampled at steady state after",
    "the fourth 1-h IV infusion of 0.25 g or 0.5 g tazobactam (with",
    "piperacillin) every 8 h. Zero-order IV input into the central compartment,",
    "first-order elimination, an additive non-renal-plus-renal clearance linear",
    "in Cockcroft-Gault creatinine clearance, and a linear body-weight effect on",
    "the central volume; the intercompartmental clearance and the peripheral",
    "volume were fixed.",
    sep = " "
  )
  reference <- paste(
    "Kim YK, Jung JA, Choi HK, Bae IG, Choi WS, Hur J, Jin SJ, Kim SW, Kwon KT,",
    "Lee SR, Shin JG, Kiem S.",
    "Population pharmacokinetic analysis of piperacillin/tazobactam in Korean",
    "patients with acute infections.",
    "Infect Chemother. 2016;48(3):209-215.",
    "doi:10.3947/ic.2016.48.3.209.",
    "Structural and covariate equations: Results section 3 ('the final model was",
    "as follows') and the Table 2 structural-model header rows.",
    "All parameter estimates: Table 2, tazobactam columns.",
    "The a priori linear body-weight adjustment of V1 follows reference [16] of",
    "the paper (Li C, Kuti JL, Nightingale CH, Mansfield DL, Dana A, Nicolau DP.",
    "J Antimicrob Chemother), cited in Methods section 4.",
    "The piperacillin counterpart fitted in the same paper is",
    "modellib('Kim_2016_piperacillin').",
    sep = " "
  )
  vignette <- "Kim_2016_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The paper puts a single log-normal eta on total CL (Table 2, 'omega CL,
  # interindividual variability of CL'), but the typical CL is a two-arm
  # additive structural equation (TVCL = theta1 + theta5 * CLcr/47), so there is
  # no single `lcl` ini parameter for `etalcl` to pair with. Declared here so
  # checkModelConventions() recognises the pairing.
  paper_specific_etas <- "etalcl"

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Kim 2016 Methods section 3 (plasma
  # samples assayed by LC-MS/MS) and Results section 3 (two-compartment
  # structural model with first-order elimination).
  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault formula",
        "(Kim 2016 Methods section 2). Absolute clearance in mL/min, NOT",
        "normalised to 1.73 m^2 body surface area."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters total clearance additively and linearly as the renal arm:",
        "TVCL (L/h) = 1.76 + 4.81 * CLcr/47 (Kim 2016 Results section 3 and",
        "Table 2 structural-model header). Adding CLcr on CL lowered the",
        "objective function value by 17.736 points. Cohort CLcr was",
        "61.27 +/- 36.67 mL/min, range 14.45-146.01 (Table 1).",
        "The normalising constant 47 mL/min is printed in the equation but the",
        "paper never states its provenance; Table 1 reports the cohort mean",
        "(61.27), not the median, so 47 cannot be confirmed as the cohort",
        "median from the text. It is reproduced exactly as printed."
      ),
      source_name = "CLcr"
    ),
    WT = list(
      description = "Total body weight at enrollment",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the central volume linearly (exponent fixed at 1, not",
        "estimated): TVV1 (L) = 22.6 * weight/60 (Kim 2016 Results section 3",
        "and Table 2 structural-model header). The paper states that V1 was",
        "'adjusted for a median body weight of 60 kg with a linear",
        "relationship' a priori, following its reference [16], rather than",
        "being selected by the covariate search. Cohort weight was",
        "58.17 +/- 10.08 kg, range 36.30-75.40 (Table 1)."
      ),
      source_name = "weight"
    )
  )

  # Screened as candidate covariates in the stepwise covariate search (Kim 2016
  # Methods section 4) but not retained in the final model. Documented here for
  # provenance; none is referenced in model().
  # The Glasgow Coma Scale score was screened alongside these but has no
  # canonical covariate-column name in inst/references/covariate-columns.md, and
  # this extraction does not introduce one because no model uses it; it is
  # recorded in population$notes and in the vignette instead.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at enrollment",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened for linear, exponential and power covariate models",
        "(Methods section 4) but not retained. Cohort 68.79 +/- 10.97 years,",
        "range 46-88 (Table 1)."
      ),
      source_name = "age"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 (male)",
      notes = paste(
        "Screened (Methods section 4) but not retained. The source reports sex",
        "as male/female counts (17 men, 16 women; Table 1), so SEXF = 1 - SEXM",
        "relative to the source column."
      ),
      source_name = "sex"
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened (Methods section 4) but not retained. No summary statistic",
        "for BUN is reported in Table 1."
      ),
      source_name = "blood urea nitrogen"
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened (Methods section 4) but not retained. Cohort",
        "13.48 +/- 8.48, range 3-38 (Table 1)."
      ),
      source_name = "APACHE II score"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 33,
    n_studies = 1,
    age_range = "46-88 years",
    age_mean = "68.79 years (SD 10.97)",
    weight_range = "36.30-75.40 kg",
    weight_mean = "58.17 kg (SD 10.08)",
    sex_female_pct = 48.5,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "acute infection requiring intravenous piperacillin/tazobactam;",
      "site of infection lungs 13, urinary tract 12, soft tissue 6,",
      "bloodstream 2"
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 61.27 +/- 36.67 mL/min,",
      "range 14.45-146.01; serum creatinine 1.18 +/- 0.80 mg/dL,",
      "range 0.34-4.70"
    ),
    dose_range = paste(
      "piperacillin/tazobactam 2/0.25 g (14 patients, 42.4%) for creatinine",
      "clearance <= 50 mL/min or 4/0.5 g (19 patients, 57.6%) for creatinine",
      "clearance > 50 mL/min, each infused intravenously over 1 h every 8 h for",
      "at least four consecutive doses"
    ),
    regions = "Korea (six university-affiliated hospitals)",
    notes = paste(
      "Baseline demographics: Table 1. Enrolled April 2013 to April 2015 in",
      "patients aged over 20 years with, or expected to have, TZP-susceptible",
      "pathogens. 35 patients completed the study; 2 were excluded from the",
      "analysis because sampling at the injection site was suspected (their",
      "concentrations were roughly 1,000-fold higher than the others), leaving",
      "33. Blood was drawn pre-dose and at 0 min, 30 min and 4-6 h after the",
      "end of the fourth infusion, i.e. at steady state.",
      "The Glasgow Coma Scale score (13.09 +/- 2.88, range 7-15; Table 1) was",
      "also screened as a candidate covariate and not retained; it has no",
      "canonical covariate-column name and is therefore not listed in",
      "covariatesDataExcluded."
    )
  )

  ini({
    # ===== Structural parameters (Kim 2016 Table 2, tazobactam columns) =====
    # Total clearance is additive in a non-renal intercept and a renal arm
    # linear in creatinine clearance:
    #   TVCL (L/h) = theta1 + theta5 * CLcr/47
    # mapped onto the registered multi-component clearance canonicals
    # lcl_nonren (intercept) and lcl_renal (slope), following the additive
    # renal-plus-non-renal precedent in Bulitta_2011_cefpirome.R and
    # Tong_2026_vancomycin_carreno.R.
    lcl_nonren <- log(1.76)
    label("Non-renal clearance intercept (L/h)")
    # Table 2 theta1 = 1.76 L/h (RSE 51.82%; bootstrap 1.72, 95% CI 0.01-3.67)
    lcl_renal <- log(4.81)
    label("Renal clearance arm at CLcr = 47 mL/min (L/h)")
    # Table 2 theta5 = 4.81 L/h (RSE 23.49%; bootstrap 4.79, 95% CI 2.86-7.71)

    lvc <- log(22.6)
    label("Central volume V1 at 60 kg (L)")
    # Table 2 theta2 = 22.6 L (RSE 10.84%; bootstrap 24.10, 95% CI 20.06-30.79)

    # Q and V2 are the only two structural rows in Table 2 with no %RSE and no
    # bootstrap median/CI ('-' in both columns for both drugs), while every
    # other row carries both. NONMEM RSEs come from the covariance matrix, so a
    # successful covariance step -- which the other rows prove -- would have
    # produced RSEs for Q and V2 had they been estimated. They were therefore
    # held fixed, and are encoded with fixed() accordingly. The paper does not
    # print a FIX flag or state the source of the two values; see the vignette
    # Errata.
    lq <- fixed(log(1.18))
    label("Intercompartmental clearance Q (L/h)")
    # Table 2 Q = 1.18 L/h (no RSE, no bootstrap: fixed)
    lvp <- fixed(log(4.3))
    label("Peripheral volume V2 (L)")
    # Table 2 V2 = 4.3 L (no RSE, no bootstrap: fixed)

    # Reference covariate values printed in the Table 2 structural-model
    # equations.
    crcl_ref <- 47
    label("Reference creatinine clearance in the CL equation (mL/min)")
    wt_ref <- 60
    label("Reference body weight in the V1 equation (kg; stated cohort median)")

    # ===== IIV (Kim 2016 Table 2, 'Random variability (CV,%)') =====
    # Kim 2016 Methods section 4: 'The inter-individual variability of each
    # parameter was described using a log-normal variance model.' No block
    # covariance is reported, so the two etas are independent.
    #
    # SCALE. Table 2 prints these as 'omega CL' / 'omega V1' -- omega is the
    # NONMEM symbol for the SD, not the variance -- under a block header reading
    # '(CV, %)', and omega is approximately the CV for a log-normal eta this
    # small, so the symbol and the header agree only on the SD reading. The
    # decisive check is the residual row in the same block: sigma = 0.393 with
    # RSE 12.26% on roughly 33 x 4 = 132 observations, and sqrt(2/132) = 12.31%
    # is the floor a variance-scale RSE column could not beat. 12.26% is below
    # that floor, so the block is printed on the SD scale; a variance reading
    # would also make the proportional residual error 62.7%, implausible for the
    # validated LC-MS/MS assay of Methods section 3. Variances below are
    # therefore omega^2. See the vignette Errata for the discarded readings.
    etalcl ~ 0.044521 # Table 2 'omega CL' = 0.211 (RSE 28.67%; bootstrap 0.198, 95% CI 0.099-0.500); 0.211^2
    etalvc ~ 0.022801 # Table 2 'omega V1' = 0.151 (RSE 68.87%; bootstrap 0.134, 95% CI 0.001-0.418); 0.151^2

    # ===== Residual error (Kim 2016 Table 2: proportional only) =====
    # Methods section 4: additive, proportional and combined error models were
    # tested; Results section 3: 'Residual variability was most effectively
    # explained by a proportional error model.'
    propSd <- 0.393
    label("Proportional residual error (fraction)")
    # Table 2 sigma, proportional error = 0.393 (RSE 12.26%; bootstrap 0.387, 95% CI 0.297-0.479)
  })

  model({
    # ----- Individual PK parameters (Kim 2016 Results section 3 equations) ----
    # TVCL (L/h) = 1.76 + 4.81 * CLcr/47, with a single log-normal eta on the
    # total, matching Table 2's one 'omega CL' row.
    cl_nonren <- exp(lcl_nonren)
    cl_renal <- exp(lcl_renal) * CRCL / crcl_ref
    cl <- (cl_nonren + cl_renal) * exp(etalcl)

    # TVV1 (L) = 22.6 * weight/60.
    vc <- exp(lvc + etalvc) * (WT / wt_ref)

    q <- exp(lq)
    vp <- exp(lvp)

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ----- ODE system -----
    # Tazobactam is given as a 1-h zero-order IV infusion into the central
    # compartment (infusion rate or duration supplied through the data-level
    # RATE / DUR column); there is no depot compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ----- Output -----
    # Total plasma tazobactam concentration: dose in mg, vc in L -> mg/L,
    # which is the ug/mL scale of the Methods section 3 calibration curve
    # (0.2-50.0 ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
