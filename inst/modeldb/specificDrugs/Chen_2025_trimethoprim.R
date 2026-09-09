Chen_2025_trimethoprim <- function() {
  description <- paste(
    "One-compartment population PK model for INTRAVENOUS trimethoprim in",
    "adults treated with co-trimoxazole for Pneumocystis jirovecii pneumonia",
    "(Chen 2025). Both the typical volume and the typical clearance are",
    "reported per kilogram of body weight (L/kg and L/kg/h), so body weight",
    "enters both parameters linearly; this is how the 'WT on CL' step",
    "retained in the covariate search is carried, and it is why the printed",
    "model equations show no explicit weight term. Creatinine clearance",
    "scales clearance as a power function normalized to the cohort median of",
    "75.7 mL/min. Unlike the companion sulfamethoxazole model, continuous",
    "renal replacement therapy is NOT a covariate here: trimethoprim has a",
    "large volume of distribution and is about 60% protein bound, so little",
    "is removed by ultrafiltration. NOTE the typical volume is taken as",
    "2.22 L/kg from the Discussion rather than the 8.22 L/kg printed in",
    "Table 6; four independent lines of evidence refute 8.22 and the vignette",
    "Errata sets them out in full. Doses are the TRIMETHOPRIM component of",
    "the combination product: a co-trimoxazole dose of X mg/kg/day delivers",
    "X/6 mg/kg/day of trimethoprim and 5X/6 mg/kg/day of sulfamethoxazole.",
    "The sulfamethoxazole model of the same paper is a separate file; see",
    "modellib('Chen_2025_sulfamethoxazole').",
    sep = " "
  )
  reference <- paste(
    "Chen B, Chen Y, Chen M, Mao Y, Huang Y, Zhou L, Wu W, Li X, Wu X,",
    "Cheng Y, Qiu H. Population pharmacokinetics and Monte Carlo-based",
    "dosing optimization of trimethoprim-sulfamethoxazole.",
    "Antimicrob Agents Chemother. 2025;69(11):e00519-25.",
    "doi:10.1128/aac.00519-25.",
    "Structural equations from Eq. 3 and Eq. 4 of the Results section",
    "('Population pharmacokinetic analysis'); fixed effects, interindividual",
    "variability and residual error from Table 6, EXCEPT the typical volume,",
    "which is taken from the Discussion ('the typical apparent volumes of V",
    "were 0.32 L/kg for SMX and 2.22 L/kg for TMP'). The equations are",
    "rendered as images in the publisher PDF and were recovered with",
    "'pdftotext -layout'.",
    sep = " "
  )
  vignette <- "Chen_2025_cotrimoxazole"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters BOTH clearance and volume LINEARLY (exponent exactly 1),",
        "because Table 6 reports the typical values per kilogram: tv V in",
        "L/kg and tv CL in L/kg/h. Weight is therefore not visible as a",
        "separate term in the printed equations (Eq. 3 and Eq. 4), but it is",
        "the 'WT on CL' covariate that the forward/backward search retained",
        "(Table 4: adding WT on CL dropped the objective function value by",
        "10.48, and removing it from the full model raised it by 8.74,",
        "P < 0.05). Because the paper's dosing recommendations are all in",
        "mg/kg and both parameters are per-kg, simulated concentrations are",
        "independent of the weight supplied when the dose is scaled by the",
        "same weight.",
        "Cohort median 60.0 kg (Table 2)."
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation.",
        "RAW mL/min, NOT normalized to 1.73 m^2 body surface area."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power term (CRCL / 75.7)^0.29 on clearance (Eq. 4). The text after",
        "the equations states explicitly that '75.7 represents the median",
        "CrCL value', matching the cohort median in Table 2 (75.7 mL/min,",
        "51.0-93.7). Methods: 'CrCL ... calculated using the Cockcroft-Gault",
        "equation', which returns raw mL/min; do NOT supply a BSA-normalized",
        "value here. The exponent is steeper than the 0.17 of the companion",
        "sulfamethoxazole model, consistent with the Discussion's statement",
        "that both agents undergo renal elimination but that trimethoprim",
        "clearance sits within the range of native renal clearance.",
        "The Monte Carlo simulations stratify on CrCL bands of <15, 15-29,",
        "30-49, 50-79 and 80-120 mL/min (Table 7).",
        sep = " "
      ),
      source_name        = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    RRT_CRRT_STATUS = list(
      description        = "Continuous renal replacement therapy during co-trimoxazole treatment",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "Screened and NOT carried into the final trimethoprim model. Table 4",
        "shows CRRT entering the forward chain as model 8 (dOFV -4.12,",
        "P < 0.05) and its removal costing 3.92 (P < 0.05) in the backward",
        "step, but the final model expression (Eq. 4) has no CRRT term and",
        "Table 6 reports no dCLdCRRT coefficient. The Results state it",
        "plainly: 'CRRT was identified as a significant covariate exclusively",
        "for the SMX model' and 'the model did not reveal a significant",
        "correlation between CRRT and TMP exposure; therefore, Monte Carlo",
        "simulations for TMP were not performed in this patient population'",
        "- Table 7 accordingly leaves the TMP cells of the CRRT row blank.",
        "The Discussion gives the mechanism: trimethoprim has a large volume",
        "of distribution, so plasma holds only a small fraction of the body",
        "burden, and about 60% protein binding further limits removal by",
        "ultrafiltration. No point estimate is reported, so nothing can be",
        "encoded. The companion sulfamethoxazole model DOES carry this",
        "covariate; see modellib('Chen_2025_sulfamethoxazole').",
        sep = " "
      ),
      source_name        = "CRRT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      notes              = paste(
        "Screened and NOT carried into the final trimethoprim model, despite",
        "surviving the search chain. Table 4 shows age entering as model 7",
        "(dOFV -7.19, P < 0.01) and its removal costing 5.91 (P < 0.01), yet",
        "the final model expression (Eq. 4) has no age term and Table 6",
        "reports no age coefficient. The Results narrative mentions only CrCL",
        "and CRRT as retained covariates, and the Conclusion repeats that",
        "'CrCL was identified as a significant covariate influencing the CL",
        "of both SMX and TMP'. With no point estimate anywhere in the paper",
        "the effect cannot be encoded; it is recorded here so the omission is",
        "traceable rather than silent. Cohort median 64 years (Table 2). The",
        "Discussion separately notes that renal clearance of trimethoprim is",
        "reduced in elderly subjects, which is the effect this term would",
        "have captured.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    HT = list(
      description        = "Height",
      units              = "cm",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 170.0 cm (Table 2)."
    ),
    SEXF = list(
      description        = "Female sex",
      units              = "(binary)",
      type               = "binary",
      notes              = "Screened as 'gender' but not retained (Methods covariate list). 18 of 79 patients female (Table 2)."
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 34.4 g/L (Table 2)."
    ),
    TBILI = list(
      description        = "Total bilirubin",
      units              = "umol/L",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 21.7 umol/L (Table 2)."
    ),
    DBIL = list(
      description        = "Direct bilirubin",
      units              = "umol/L",
      type               = "continuous",
      notes              = "Screened but not retained (Methods covariate list). Cohort median 10.9 umol/L (Table 2)."
    ),
    HEPIMP_SEV = list(
      description        = "Hepatic impairment severity (Child-Pugh class)",
      units              = "(categorical)",
      type               = "categorical",
      notes              = "Screened as the Child-Pugh classification (A/B/C) but not retained (Methods covariate list; Results 'Liver function ... did not exhibit statistically significant effects'). Child-Pugh A 54, B 18, C 7 (Table 2). Total protein, ALT, AST, GGT and ALP were screened alongside it and are likewise absent from the final model; NAT2 acetylator phenotype and CYP2C9 metabolizer phenotype were genotyped and screened with no significant effect. None carries a reported point estimate, so none can be encoded."
    )
  )

  compartmentData <- list(
    central = list(analyte = "trimethoprim", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 79,
    n_studies      = 1,
    age_median     = "64 years (54-73)",
    height_median  = "170.0 cm (165-175)",
    weight_median  = "60.0 kg (55-70)",
    sex_female_pct = 22.8,
    race_ethnicity = "Not reported; single-center Chinese cohort.",
    disease_state  = paste(
      "Adults (>= 18 years) with a confirmed diagnosis of Pneumocystis",
      "jirovecii pneumonia receiving intravenous co-trimoxazole, in intensive",
      "care units and general inpatient wards. Patients taking hepatic enzyme",
      "inducers (carbamazepine, rifampin) were excluded.",
      sep = " "
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCL median 75.7 mL/min (51.0-93.7); serum creatinine",
      "median 77.8 umol/L (59.0-135.0). 47 of 79 patients (59.5%) had",
      "CrCL < 80 mL/min and 32 (40.5%) had CrCL >= 80 mL/min. 19 patients",
      "(24.1%) received continuous renal replacement therapy, which was not a",
      "covariate in this model.",
      sep = " "
    ),
    hepatic_function = paste(
      "Child-Pugh A 54 (68.3%), B 18 (22.8%), C 7 (8.9%). Total bilirubin",
      "median 21.7 umol/L, albumin 34.4 g/L, ALT 34.5 U/L, AST 32.0 U/L.",
      sep = " "
    ),
    dose_range     = paste(
      "Routine care rather than protocol-assigned; intravenous infusions of",
      "approximately 1 h given every 6, 8 or 12 h. The Monte Carlo",
      "simulations explore co-trimoxazole 50, 55, 65, 70 and 90 mg/kg/day",
      "given twice, three times or four times daily.",
      sep = " "
    ),
    regions        = "China (Fujian Medical University Union Hospital, Fuzhou), March 2023 to October 2024.",
    n_observations = paste(
      "232 post-dose plasma concentrations from 79 patients: two to three",
      "samples per patient (an end-of-infusion peak, a pre-dose trough and/or",
      "an intermediate sample; Table 1). Trimethoprim was assayed by LC-MS/MS",
      "over a calibrated range of 0.20-25.0 mg/L.",
      sep = " "
    ),
    notes          = paste(
      "Prospective single-center study; demographics from Table 2. Estimation",
      "was by first-order conditional estimation with extended least squares",
      "in Phoenix NLME 8.0, and the reported parameter precision comes from",
      "1,000-sample nonparametric bootstrapping. The efficacy/toxicity window",
      "the paper simulates against is a steady-state peak of 5-10 mg/L for",
      "trimethoprim. Several arithmetic inconsistencies in the published",
      "tables are catalogued in the vignette Errata; the one that changes a",
      "parameter is the typical volume - see the lvc comment below.",
      sep = " "
    )
  )

  ini({
    # Structural parameters, Table 6 ("Population pharmacokinetic model
    # estimates and bootstrap results for TMP"), "Final model / Estimate"
    # column, EXCEPT lvc. Both are per kilogram of body weight; model()
    # multiplies by WT.
    #
    # lvc: the paper reports TWO different typical volumes for trimethoprim.
    # Table 6 prints 8.22 L/kg (bootstrap median 8.24, 95% CI 7.37-9.08); the
    # Discussion prints "the typical apparent volumes of V were 0.32 L/kg for
    # SMX and 2.22 L/kg for TMP, consistent with prior adult pharmacokinetic
    # studies (SMX: 0.17-0.34 L/kg, TMP: 1.0-2.4 L/kg)". 2.22 is used here
    # because four independent lines of evidence refute 8.22:
    #  (1) Absolute volume. At the cohort median weight of 60 kg, 2.22 L/kg is
    #      133 L and 8.22 L/kg is 493 L. The paper's own cited literature
    #      range of 1.0-2.4 L/kg brackets 2.22 and excludes 8.22.
    #  (2) Half-life. With the reported CL of 0.11 L/kg/h, 2.22 L/kg gives a
    #      terminal half-life of 12.9 h at the median CrCL, close to
    #      trimethoprim's well-established 8-12 h; 8.22 L/kg gives 47.8 h.
    #  (3) The paper's own dose-frequency simulations. The Discussion reports
    #      that at 90 mg/kg/day the fraction of patients exceeding the 10 mg/L
    #      toxicity threshold rises from 20% on a q6h schedule to 45.8% on a
    #      q12h schedule. That shift requires the steady-state peak to rise
    #      about 1.35-fold when the same daily dose is given twice rather than
    #      four times a day. 2.22 L/kg predicts 1.16-fold; 8.22 L/kg predicts
    #      1.04-fold, i.e. essentially no frequency effect at all, because a
    #      47.8 h half-life flattens the profile within a dosing interval. The
    #      same calculation applied to sulfamethoxazole, whose parameters are
    #      not in dispute, predicts 1.19-fold against a required 1.30-fold, so
    #      the method reads about 8% low and 2.22 sits within that bias while
    #      8.22 does not.
    #  (4) The Discussion sentence is internally coherent - an author who had
    #      fitted 8.22 L/kg could not describe it as consistent with a
    #      1.0-2.4 L/kg literature range.
    # See the vignette Errata for the full worked arbitration.
    lvc <- log(2.22); label("Typical weight-normalized central volume of distribution (L/kg)")               # Discussion, not Table 6: "the typical apparent volumes of V were 0.32 L/kg for SMX and 2.22 L/kg for TMP". Table 6 "tv V (L/kg)" prints 8.22 (RSE 5.28%, bootstrap median 8.24, 95% CI 7.37-9.08) and is treated as erroneous - see the block comment above and the vignette Errata.
    lcl <- log(0.11); label("Typical weight-normalized clearance at CrCL 75.7 mL/min (L/kg/h)")              # Table 6 "tv CL (L/kg/h)": 0.11, RSE 5.55%, bootstrap median 0.1 (95% CI 0.05-0.21). Discussion: "typical population CL values of ... 0.11 L/kg/h for TMP, aligning with previously reported ranges (TMP: 0.071-0.11 L/kg/h)". At 60 kg this is 6.6 L/h = 110 mL/min, matching trimethoprim's literature plasma clearance.

    # Covariate effect on clearance, Eq. 4:
    #   CL = tvCL * (CrCL/75.7)^dCLdCrCL * exp(etaCL)
    # There is no CRRT term in the trimethoprim model; see
    # covariatesDataExcluded above.
    e_crcl_cl <- 0.29; label("Power exponent on (CRCL/75.7) for CL (unitless)")                              # Table 6 "dCLdCrCL": 0.29, RSE 23.91%, bootstrap median 0.29 (95% CI 0.15-0.42)

    # Interindividual variability. Table 6 labels these rows "omega^2" and the
    # footnote confirms "omega^2, variance of interindividual variability", so
    # the printed numbers are variances on the log scale and go into ini()
    # unchanged. Methods: "Interindividual variability was characterised using
    # an exponential residual model", matching the exp(eta) form of Eq. 3-4.
    etalvc ~ 0.19                                                                                            # Table 6 "omega^2 V": 0.19, RSE 16.78%, bootstrap median 0.18 (95% CI 0.13-0.25); omega = 43.6% on the log scale
    etalcl ~ 0.25                                                                                            # Table 6 "omega^2 CL": 0.25, RSE 3.80%, bootstrap median 0.25 (95% CI 0.14-0.25); omega = 50.0% on the log scale

    # Residual error. Results: "residual variability was appropriately
    # characterised using a proportional error model" for both analytes;
    # Table 4 confirms the proportional model was selected on OFV (557.15 vs
    # 559.47 additive and 560.66 combined).
    propSd <- 0.27; label("Proportional residual error SD (fraction)")                                       # Table 6 "Proportional error", "Final model / Estimate" column: 0.27, RSE 7.84%. The Table 6 bootstrap cells for this row (median 0.06, 95% CI 0.05-0.07) are on a different scale - 0.27^2 = 0.073 - i.e. the bootstrap column appears to report the VARIANCE where the estimate column reports the SD. The estimate column is used, for consistency with the sulfamethoxazole model where both columns agree at 0.09. See the vignette Errata.
  })

  model({
    # Cohort median creatinine clearance used as the normalizing constant.
    # Stated in the sentence following Eq. 1-4: "75.7 represents the median
    # CrCL value"; matches the Table 2 median of 75.7 mL/min.
    ref_crcl <- 75.7

    # Eq. 3: V = tvV * exp(etaV), with tvV per kg, so V (L) = tvV * WT.
    vc <- exp(lvc + etalvc) * WT

    # Eq. 4: CL = tvCL * (CrCL/75.7)^dCLdCrCL * exp(etaCL), with tvCL per kg,
    # so CL (L/h) = tvCL * WT * (CrCL/75.7)^dCLdCrCL * exp(etaCL).
    cl <- exp(lcl + etalcl) * WT * (CRCL / ref_crcl)^e_crcl_cl

    kel <- cl / vc

    # One compartment with first-order elimination (Results: "The PopPK
    # characteristics of both SMX and TMP were best described by a
    # one-compartment model with first-order elimination kinetics"; Table 4
    # model 1 OFV 559.47 versus 579.49 for two compartments). Intravenous
    # only - the analysis excludes oral administration - so there is no depot
    # and no bioavailability term; dose directly into central.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
