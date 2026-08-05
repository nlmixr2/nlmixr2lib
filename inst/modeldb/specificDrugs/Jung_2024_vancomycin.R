Jung_2024_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in critically ill children (90 days to <18 years) in the PICU who are not on extracorporeal therapy (Jung 2024). Clearance and intercompartmental clearance scale allometrically with body weight (exponent 0.75, reference 20 kg) and clearance additionally scales as a power function of bedside-Schwartz estimated GFR (exponent 0.5259, reference 141 mL/min/1.73 m^2); central and peripheral volumes scale linearly with body weight (exponent 1, reference 20 kg). Residual variability is a log error model."
  reference <- "Jung D, Kishk OA, Bhutta AT, Cummings GE, El Sahly HM, Virk MK, Moffett BS, Morris Daniel JL, Watanabe A, Fishbane N, Kotloff KL, Gu K, Ghazaryan V, Gobburu JVS, Akcan-Arikan A, Campbell JD. Evaluation of Vancomycin Dose Needed to Achieve 24-Hour Area Under the Concentration-Time Curve to Minimum Inhibitory Concentration Ratio Greater Than or Equal to 400 Using Pharmacometric Approaches in Pediatric Intensive Care Patients. Crit Care Explor. 2024;6(10):e1159. doi:10.1097/CCE.0000000000001159"
  vignette <- "Jung_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin is given as an IV infusion, so the dose
  # enters `central` directly and there is no depot state. `central` is
  # verified: Jung 2024 Methods states that SERUM vancomycin levels were
  # assayed by validated LC-MS/MS. `peripheral1` is left unverified because
  # the paper never states what matrix the peripheral distribution
  # compartment represents; "plasma" follows the repository default for a
  # mathematical peripheral compartment and is not a paper-sourced claim.
  compartmentData <- list(
    central     = list(analyte = "vancomycin", units = "mg", specimen = "serum",  verified = TRUE),
    peripheral1 = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying in the source analysis (eTable 1 lists body weight as time-varying, taken at the closest time prior to the vancomycin dose for which a sample was drawn). Jung 2024 Table 1 baseline median 20.1 kg (IQR 10.7-41.6). Reference weight 20 kg, used in all four allometric terms of Jung 2024 Table 2.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the bedside Schwartz equation, BSA-normalized: eGFR (mL/min/1.73 m^2) = 0.413 * height (cm) / serum creatinine (mg/dL)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying in the source analysis (eTable 1). Jung 2024 Table 1 footnote a gives the bedside Schwartz formula; Table 1 reports the first-chemistry median 138.3 mL/min/1.73 m^2 (IQR 101.5-179.8). Reference value 141 mL/min/1.73 m^2 is the normalizing constant in the Jung 2024 Table 2 clearance equation and is the value the paper quotes for its typical subject (CL = 2.76 L/hr at eGFR 141). Renal-impairment strata (Table 1, FDA 2020 categories): normal >=90 in 81.8%, mild 60-89 in 11.9%, moderate 30-59 in 5.0%, severe 15-29 in 1.0%. Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates; the assay form here is the bedside Schwartz creatinine-based estimate. The paper also refit the final model using the CKiD cystatin-C-based eGFR equation in the 248-subject subset with both serum creatinine and cystatin C and found equivalent predictive performance, but the reported final model (Table 2) uses bedside Schwartz.",
      source_name        = "eGFR"
    )
  )

  # Screened in the Jung 2024 covariate search (eTable 1 / eTable 2) but NOT
  # retained in the final model, so they are documentation only and are not
  # referenced in model(). eGFR was the sole statistically significant
  # predictor of clearance at the prespecified alpha = 0.01.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and Vc (eTable 1). In the eTable 2 step-2 forward selection, after bedside-Schwartz eGFR was already in the model, a power effect of age on CL gave a VOF change of -15.381 (p < 0.01) but was dropped in backward elimination at alpha = 0.001. Not retained in the Jung 2024 Table 2 final model.",
      source_name        = "Age"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL and Vc (eTable 1); highly correlated with body weight, so only one of the pair could enter a single covariate sub-model (eMethods section 5). Not retained.",
      source_name        = "BMI"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Screened on CL and Vc (eTable 1, coded 0 = male, 1 = female, the same orientation as canonical SEXF). Not retained.",
      source_name        = "Sex"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL (eTable 1). A power effect gave a VOF change of -16.057 in eTable 2 step 1, less than the -74.51 achieved by bedside-Schwartz eGFR, which uses serum creatinine as an input. Not retained as a standalone covariate.",
      source_name        = "Scr"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL (eTable 1). A power effect gave a VOF change of -36.237 in eTable 2 step 1, second only to bedside-Schwartz eGFR; not retained after eGFR entered the model.",
      source_name        = "BUN"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 301L,
    n_studies        = 1L,
    n_sites          = 2L,
    n_concentrations = 1027L,
    age_range        = "90 days to <18 years (eligibility); Table 1 median 6.0 years, IQR 1.6-13.0",
    age_median       = "6.0 years (IQR 1.6-13.0)",
    weight_range     = "IQR 10.7-41.6 kg (full range not reported)",
    weight_median    = "20.1 kg (IQR 10.7-41.6)",
    sex_female_pct   = 41.7,
    race_ethnicity   = c(White = 72.8, `Black or African American` = 16.9, Other = 10.3, `Hispanic or Latino` = 40.7),
    disease_state    = "Critically ill children admitted to a quaternary-care PICU and prescribed IV vancomycin for any indication. Excluded: any form of extracorporeal therapy (ECMO, CRRT, extracorporeal liver support), cardiopulmonary bypass within 7 days of starting vancomycin, end-stage renal disease on chronic dialysis, single-dose-only vancomycin. Day-of-admission PRISM-3 median 5.0 (IQR 2.0-10.0); PELOD median 11.0 (IQR 2.0-21.0); at least one nephrotoxic comedication in 96.4%.",
    dose_range       = "Clinician-chosen IV regimens; actual amounts 5-31 mg/kg q6h, most frequently 1-hour infusions of 15 mg/kg every 6-9 hours",
    regions          = "United States (University of Maryland Children's Hospital, Baltimore MD; Texas Children's Hospital, Houston TX)",
    renal_function   = "Bedside-Schwartz eGFR median 138.3 mL/min/1.73 m^2 (IQR 101.5-179.8); normal 81.8%, mild 11.9%, moderate 5.0%, severe 1.0%, unknown 0.3%",
    notes            = "Prospective multicenter population PK study conducted 2018-04-15 to 2020-02-06 (Jung 2024 Table 1; eMethods). 302 subjects were enrolled into the non-extracorporeal-therapy group and 1 was excluded for unavailable PK data, leaving 301 in the analysis; Table 1 baseline percentages are reported over the 302 enrolled. A median of 4 (IQR 3-6) samples per subject were drawn over up to 7 days, giving 1027 concentration records (the n for the time-varying characteristics in Table 1). Serum vancomycin was assayed by validated LC-MS/MS. A separate cohort of 33 subjects on extracorporeal therapy was recruited but excluded from this analysis. Model fit in NONMEM 7.3.0 with FOCE-I; covariate selection by forward selection at alpha = 0.01 (dVOF 6.64) then backward elimination at alpha = 0.001 (dVOF 10.83). Reported eta shrinkage 13.2% on CL and 48.2% on Vc."
  )

  ini({
    # Structural parameters (Jung 2024 Table 2). The reference subject weighs
    # 20 kg and has a bedside-Schwartz eGFR of 141 mL/min/1.73 m^2.
    lcl <- log(2.763); label("Clearance at WT=20 kg and eGFR=141 mL/min/1.73 m^2 (L/h)")  # Jung 2024 Table 2: 2.763 L/hr (RSE 2.746%)
    lvc <- log(8.63);  label("Central volume at WT=20 kg (L)")                            # Jung 2024 Table 2: 8.63 L (RSE 14.39%)
    lq  <- log(3.808); label("Intercompartmental clearance at WT=20 kg (L/h)")             # Jung 2024 Table 2: 3.808 L/hr (RSE 20.45%)
    lvp <- log(9.566); label("Peripheral volume at WT=20 kg (L)")                          # Jung 2024 Table 2: 9.566 L (RSE 14.46%)

    # Allometric exponents. Jung 2024 Table 2 prints them inside the parameter
    # equations without an RSE (0.75 on CL and Q; an implicit 1 on Vc and Vp),
    # and eTable 2 names the reference model "2-compartment model with
    # allometric scaling", so both are canonical fixed allometric exponents.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on (WT/20) for CL and Q (unitless)")   # Jung 2024 Table 2: (WT/20)^0.75 on CL and on Q
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on (WT/20) for Vc and Vp (unitless)")  # Jung 2024 Table 2: (WT/20) on Vc and on Vp

    # Renal-function effect on clearance. Estimated, so not wrapped in fixed().
    e_crcl_cl <- 0.5259; label("Power exponent on (CRCL/141) for CL (unitless)")  # Jung 2024 Table 2: (eGFR/141)^0.5259 (RSE 14.92%)

    # Between-subject variability (Jung 2024 Table 2, reported as %CV).
    # eMethods section 4.1 defines the exponential BSV model X_i = Xtilde_i *
    # exp(eta) and states %CV = sqrt(exp(omega^2) - 1) * 100, so
    # omega^2 = log(CV^2 + 1).
    etalcl ~ 0.135678  # log(0.3812^2 + 1); Jung 2024 Table 2: 38.12% CV on CL (RSE 15.94%)
    etalvc ~ 0.719144  # log(1.026^2  + 1); Jung 2024 Table 2: 102.6% CV on Vc (RSE 28.64%)
    # BSV on Q and Vp is reported as "NE" (not estimated) in Jung 2024 Table 2,
    # so no etalq / etalvp is included.

    # Residual variability. eMethods section 4.2.4 specifies the log error
    # model log(Y) = log(Yhat) + eps with eps ~ N(0, sigma^2), reported as a
    # standard deviation in log-transformed concentration units. That is
    # nlmixr2's lnorm() error structure, parameterized by expSd.
    expSd <- 0.3310; label("Log-scale residual SD (log ug/mL)")  # Jung 2024 Table 2: variance 0.1096 (RSE 12.69%), sd 0.3310 = sqrt(0.1096)
  })
  model({
    # Individual PK parameters (Jung 2024 Table 2):
    #   CL = 2.763 * (WT/20)^0.75 * (eGFR/141)^0.5259
    #   Vc = 8.63  * (WT/20)
    #   Q  = 3.808 * (WT/20)^0.75
    #   Vp = 9.566 * (WT/20)
    cl <- exp(lcl + etalcl) * (WT / 20)^e_wt_cl_q * (CRCL / 141)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 20)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 20)^e_wt_cl_q
    vp <- exp(lvp)          * (WT / 20)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L, so central/vc is mg/L == ug/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
